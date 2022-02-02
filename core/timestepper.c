/**************************************************************************
 *                                                                        *
 *  TIMESTEPPER.C                                                         *
 *                                                                        *
 *  UPDATES PRIMITIVES USING THE GRIM ALGO:                               *
 *  1. COMPUTES THE RESIDUAL                                              *
 *  2. NUMERICALLY COMPUTES THE JACOBIAN                                  *
 *  3. SOLVES THE LINEAR SYSTEM FOR THE DELTA_P                           *
 *  4. FINDS THE RIGHT LINE SEARCH PARAMETER (TODO)                       *
 *  5. ITERATES TILL TOLERANCE CRITERION IS MET                           *
 *                                                                        *
 **************************************************************************/

// TODO: Put timer calls for grim calls

#include "decs.h"

#if GRIM_TIMESTEPPER

#include "mkl.h"

// Function declaration(s): residual calculation, jacobian, linear solve, L2 norm
void residual_calc(struct GridGeom *G, struct FluidState *Stmp, double U_old[NVAR], double divF[NVAR], 
  double sources[NVAR], double dt, int i, int j, int k double residual[NVAR]);

void jacobian(struct GridGeom *G, struct FluidState *Sf, double U_old[NVAR],
  double divF[NVAR], double sources[NVAR], double dt, int i, int j, int k, double J[NVAR*NVAR]);

void solve(double A[NVAR*NVAR], double b[NVAR]);

void grim_timestep(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidState *Sf, struct FluidFlux *F, double dt) {

  int nonlinear_iter = 0;

  // TODO: Source calc term
  while (nonlinear_iter < MAX_NONLINEAR_ITER){

    double max_norm = 0;

  #pragma omp parallel for reduction(max:max_norm) collapse(2)
    ZLOOP {

      // Variable declarations0
      static struct FluidState *S_guess; // Need this since prim_to_flux and mhd_calc require a FluidState object as parameter
      double prim_guess[NVAR] = {0};
      double residual[NVAR] = {0};
      double jac[NVAR*NVAR] = {0};
      double delta_prim[NVAR] = {0};
      double norm = 0;

      static int firstc = 1;
      if (firstc){
        S_guess = calloc(1, sizeof(struct FluidState));
        firstc = 0;
      }

      // Initialize U_old with Si->U
      double U_old[NVAR] = {0};
      PLOOP U_old[ip] = Si->U[ip][k][j][i];

      // Compute flux divergence
      double divF[NVAR] = {0}; 
      PLOOP divF[ip] = (F->X1[ip][k][j][i+1] - F->X1[ip][k][j][i])/dx[1]
            + (F->X2[ip][k][j+1][i] - F->X2[ip][k][j][i])/dx[2]
            + (F->X3[ip][k+1][j][i] - F->X3[ip][k][j][i])/dx[3];

      // Compute source term
      double sources[NVAR] = {0};
      get_fluid_source(G, Ss, sources, i, j, k);

      // Jacobian calculation
      jacobian(G, Sf, U_old, divF, sources, dt, i, j, k, jac);

      // Initialize prim_guess
      PLOOP {
        prim_guess[ip] = Sf->P[ip][k][j][i];
        S_guess->P[ip][k][j][i] = prim_guess[ip];
      }

      // Assign delta_prim = -residual(P)
      residual_calc(G, S_guess, U_old, divF, sources, dt, i, j, k, residual);
      PLOOP delta_prim[ip] = -residual[ip];

      // Linear solve
      solve(jac, delta_prim);

      // Update prim_guess
      PLOOP { // TODO: Add linesearch: prim_guess = prim_guess + lambda*delta_prim, lambda in (0,1]
        prim_guess[ip] += delta_prim[ip];
        S_guess->P[ip][k][j][i] = prim_guess[ip];
      }

      residual_calc(G, S_guess, U_old, divF, sources, dt, i, j, k, residual);

      // L2 norm
      PLOOP norm += pow(residual[ip], 2);
      norm = sqrt(norm);
      if (norm > max_norm) max_norm = norm;

      // Update Sf->P
      PLOOP Sf->P[ip][k][j][i] = prim_guess[ip];
      
    }// END OF ZLOOP

    max_norm = mpi_max(max_norm);

    if (mpi_io_proc()) fprintf(stdout, "\n\t\t\t\tNonlinear iter = %d. Max L2 norm: %g\n", nonlinear_iter, max_norm);

    nonlinear_iter += 1;

  } // END OF nonlinear solver while loop. prim_guess is now the solution

}

// Residual calculation
void residual_calc(struct GridGeom *G, struct FluidState *Stmp, double U_old[NVAR], double divF[NVAR], 
  double sources[NVAR], double dt, int i, int j, int k, double residual[NVAR]) {

  // Compute conserved variables from prims
  get_state(G, Stmp, i, j, k, CENT);
  prim_to_flux(G, Stmp, i, j, k, 0, CENT, Stmp->U);
  // Compute residual 
  PLOOP residual[ip] = (Stmp->U[ip][k][j][i] - U_old[ip])/dt + divF[ip] - sources[ip];
}

// Evaluate Jacobian (per zone)
void jacobian(struct GridGeom *G, struct FluidState *Sf, double U_old[NVAR],
  double divF[NVAR], double sources[NVAR], double dt, int i, int j, int k, double J[NVAR*NVAR]) {

  // Declarations
  int is_small;
  static struct FluidState *S_eps; // Need this since prim_to_flux and mhd_calc require a FluidState object as parameter
  double prims[NVAR];
  double prims_eps[NVAR];
  double residual[NVAR];
  double residual_eps[NVAR];

  // Initialize prims and residual vectors
  static int firstc = 1;
  if (firstc){
    S_eps = calloc(1, sizeof(struct FluidState));
    firstc = 0;
  }

  PLOOP {
    prims[ip] = Sf->P[ip][k][j][i];
    prims_eps[ip] = Sf->P[ip][k][j][i];
    S_eps->P[ip][k][j][i] = prims_eps[ip];
    residual[ip] = 0;
    residual_eps[ip] = 0;
  }
  
  // Calculate residual for Sf->P
  residual_calc(G, Sf, U_old, divF, sources, dt, i, j, k, residual);

  // Numerically evaluate the Jacobian
  for (int col = 0; col < NVAR; col++) {

    // Evaluate small(P)
    is_small = 0;
    if (abs(prims[col]) < 0.5 * JACOBIAN_EPS) is_small = 1;

    // Compute P_eps
    prims_eps[col] = prims[col] + JACOBIAN_EPS*prims[col]*(1 - is_small) + JACOBIAN_EPS*is_small;
    S_eps->P[col][k][j][i] = prims_eps[col];

    // Compute the residual for P_eps
    residual_calc(G, S_eps, U_old, divF, sources, dt, i, j, k, residual_eps);

    for (int row = 0; row < NVAR; row++) {
      J[NVAR*row + col] = (residual_eps[row] - residual[row]) / (prims_eps[col] - prims[col]);
    } // END of row loop

    // Reset P_eps
    prims_eps[col] = prims[col];
    S_eps->P[col][k][j][i] = prims_eps[col];

  } // END OF col loop
}    

// Linear solve (per zone)
void solve(double A[NVAR*NVAR], double b[NVAR]){

  static int ipiv[NVAR];

  LAPACKE_dgesv(LAPACK_ROW_MAJOR, NVAR, 1, A, NVAR, ipiv, b, 1);
}

#endif
