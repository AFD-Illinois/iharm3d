/**************************************************************************
 *                                                                        *
 *  TIMESTEPPER.C                                                         *
 *                                                                        *
 *  UPDATES PRIMITIVES USING THE GRIM ALGO:                               *
 *  1. COMPUTES THE RESIDUAL                                              *
 *  2. NUMERICALLY COMPUTES THE JACOBIAN                                  *
 *  3. SOLVES THE LINEAR SYSTEM FOR THE DELTA_P                           *
 *  4. FINDS THE RIGHT LINE SEARCH PARAMETER                              *
 *  5. ITERATES TILL TOLERANCE CRITERION IS MET                           *
 *                                                                        *
 **************************************************************************/

// TODO: Put timer calls for grim calls

#include "decs.h"

#if IMEX

#include "mkl.h"

// Function declaration(s): residual calculation, jacobian, and linear solve
void residual_calc(struct GridGeom *G, struct FluidState *Stmp, struct FluidState *Si, struct FluidState *Ss,
  double U_old[NFVAR], double divF[NFVAR], double sources_explicit[NFVAR], double sources_implicit_old[NFVAR], 
  double dt, int i, int j, int k, double residual[NFVAR]);

void jacobian(struct GridGeom *G, struct FluidState *S_solver, struct FluidState *Si, struct FluidState *Ss, 
  double U_old[NFVAR], double divF[NFVAR], double sources_explicit[NFVAR], double sources_implicit_old[NFVAR],
  double dt, int i, int j, int k, double J[NFVAR*NFVAR]);

void solve(double A[NFVAR*NFVAR], double b[NFVAR]);

void imex_timestep(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidState *S_solver, struct FluidFlux *F, double dt) {

  int nonlinear_iter = 0;

  static struct FluidState *S_guess; // Need this since prim_to_flux and mhd_calc require a FluidState object as parameter
  #if LINESEARCH
  static struct FluidState *S_linesearch; // Need this since prim_to_flux and mhd_calc require a FluidState object as parameter
  #endif
  static int firstc = 1;
  if (firstc){
    S_guess = calloc(1, sizeof(struct FluidState));
    #if LINESEARCH
    S_linesearch = calloc(1, sizeof(struct FluidState));
    #endif
    firstc = 0;
  }

  while (nonlinear_iter < max_nonlinear_iter){

    double max_norm = 0;

  #pragma omp parallel for reduction(max:max_norm) collapse(2)
    ZLOOP {

      if (pflag[k][j][i] == 0.) {

        // First off, compute explicit source term, old implicit source term and store old conserved variables
        double sources_explicit[NFVAR]     = {0};
        double sources_implicit_old[NFVAR] = {0};
        double U_old[NFVAR] = {0};
        explicit_sources(G, Ss, CENT, i, j, k, sources_explicit);
        implicit_sources(G, Si, Ss, CENT, i, j, k, sources_implicit_old);
        FLOOP U_old[ip] = Si->U[ip][k][j][i]; // Initialize U_old with Si->U

        // Compute flux divergence
        double divF[NFVAR] = {0}; 
        FLOOP divF[ip] = (F->X1[ip][k][j][i+1] - F->X1[ip][k][j][i])/dx[1]
              + (F->X2[ip][k][j+1][i] - F->X2[ip][k][j][i])/dx[2]
              + (F->X3[ip][k+1][j][i] - F->X3[ip][k][j][i])/dx[3];

        // Variable declarations
        double prim_guess[NVAR] = {0};

        double residual[NFVAR] = {0};
        double jac[NFVAR*NFVAR] = {0};
        double delta_prim[NFVAR] = {0};
        double norm = 0;

        // Update Ss magnetic field primitives now that explicit sources have been computed
        BLOOP Ss->P[ip][k][j][i] = S_solver->P[ip][k][j][i];

        // Initialize prim_guess
        PLOOP {
          prim_guess[ip] = S_solver->P[ip][k][j][i];
          S_guess->P[ip][k][j][i] = prim_guess[ip];
          #if LINESEARCH
          S_linesearch->P[ip][k][j][i] = prim_guess[ip];
          #endif
        }

        // Assign delta_prim to -residual(P)
        residual_calc(G, S_guess, Si, Ss, U_old, divF, sources_explicit, sources_implicit_old, dt, i, j, k, residual);
        FLOOP delta_prim[ip] = -residual[ip];

        // Jacobian calculation
        jacobian(G, S_solver, Si, Ss, U_old, divF, sources_explicit, sources_implicit_old, dt, i, j, k, jac);

        // Linear solve
        solve(jac, delta_prim);

        // Linesearch
        double step_length  = 1.;
        #if LINESEARCH
        // L2 norm
        FLOOP norm += pow(residual[ip], 2.);
        norm = sqrt(norm);
        
        double f0      = 0.5 * norm;
        double fprime0 = -2. * f0;

        int linesearch_iter = 0;
        while (linesearch_iter < max_linesearch_iter) {
          // Take step
          FLOOP S_linesearch->P[ip][k][j][i] = S_guess->P[ip][k][j][i] + (step_length * delta_prim[ip]);

          // Compute norm of residual (loss fuction)
          residual_calc(G, S_linesearch, Si, Ss, U_old, divF, sources_explicit, sources_implicit_old, dt, i, j, k, residual);
          norm        = 0.;
          FLOOP norm += pow(residual[ip], 2);
          norm        = sqrt(norm);
          double f1   = 0.5 * norm;

          // Compute new step length
          int condition = f1 > (f0 * (1. - linesearch_eps * step_length) + SMALL);
          double denom  = (f1 - f0 - (fprime0 * step_length)) * condition + (1 - condition);
          double step_length_new = -fprime0 * step_length * step_length / denom / 2.;
          step_length = step_length * (1 - condition) + (condition * step_length_new);

          // Check if new solution has converged to required tolerance
          if (condition == 0) break;

          // Update linesearch counter
          linesearch_iter += 1;
        }
        #endif

        // Update prim_guess
        FLOOP prim_guess[ip] += (step_length * delta_prim[ip]);
        PLOOP S_guess->P[ip][k][j][i] = prim_guess[ip];

        residual_calc(G, S_guess, Si, Ss, U_old, divF, sources_explicit, sources_implicit_old, dt, i, j, k, residual);

        // L2 norm
        norm = 0.;
        FLOOP norm += pow(residual[ip], 2);
        norm = sqrt(norm);
        if (norm > max_norm) max_norm = norm;
        // Update imex_errors
        imex_errors[k][j][i] = norm;

        // Update Sf->P
        FLOOP S_solver->P[ip][k][j][i] = prim_guess[ip];

        // Check for convergence
        if (norm < rootfind_tol)
          pflag[k][j][i] = 1.;
      }

      
    }// END OF ZLOOP

    max_norm = mpi_max(max_norm);

    if (mpi_io_proc()) fprintf(stdout, "\t\tNonlinear iter = %d. Max L2 norm: %g\n", nonlinear_iter, max_norm);

    nonlinear_iter += 1;

  } // END OF nonlinear solver while loop. prim_guess is now the solution

}

// Residual calculation
void residual_calc(struct GridGeom *G, struct FluidState *Stmp, struct FluidState *Si, struct FluidState *Ss,
  double U_old[NFVAR], double divF[NFVAR], double sources_explicit[NFVAR], double sources_implicit_old[NFVAR], 
  double dt, int i, int j, int k, double residual[NFVAR]) {

  // Compute Stmp->U . These are new conserved variables
  get_state(G, Stmp, i, j, k, CENT);
  prim_to_flux(G, Stmp, i, j, k, 0, CENT, Stmp->U);

  // Compute new implicit source terms and time derivative source terms
  double sources_implicit_new[NFVAR]    = {0};
  double sources_time_derivative[NFVAR] = {0};
  implicit_sources(G, Stmp, Ss, CENT, i, j, k, sources_implicit_new);
  time_derivative_sources(G, Stmp, Si, Ss, dt, CENT, i, j, k, sources_time_derivative);
  
  // Compute residual
  FLOOP residual[ip] = (Stmp->U[ip][k][j][i] - U_old[ip])/dt + divF[ip] 
    - sources_explicit[ip] - 0.5*(sources_implicit_new[ip] + sources_implicit_old[ip]) - sources_time_derivative[ip];

  #if EMHD
  #if CONDUCTION
  // Normalize the residuals
  if (higher_order_terms_conduction == 1){

    double rho      = Ss->P[RHO][k][j][i];
    double Theta    = Ss->Theta[k][j][i];
    double tau      = Ss->tau[k][j][i];
    double chi_emhd = Ss->chi_emhd[k][j][i];

    residual[Q_TILDE] *= sqrt(rho * chi_emhd * tau * pow(Theta, 2));
  }
  else {

    double tau = Ss->tau[k][j][i];

    residual[Q_TILDE] *= tau;
  }
  #endif
  #if VISCOSITY
  if (higher_order_terms_viscosity == 1){

    double rho      = Ss->P[RHO][k][j][i];
    double Theta    = Ss->Theta[k][j][i];
    double tau      = Ss->tau[k][j][i];
    double nu_emhd  = Ss->nu_emhd[k][j][i];

    residual[DELTA_P_TILDE] *= sqrt(rho * nu_emhd * tau * Theta);
  }
  else {

    double tau = Ss->tau[k][j][i];

    residual[DELTA_P_TILDE] *= tau;
  }
  #endif
  
  #endif
}

// Evaluate Jacobian (per zone)
void jacobian(struct GridGeom *G, struct FluidState *S_solver, struct FluidState *Si, struct FluidState *Ss,
  double U_old[NFVAR], double divF[NFVAR], double sources_explicit[NFVAR], double sources_implicit_old[NFVAR],
  double dt, int i, int j, int k, double J[NFVAR*NFVAR]) {

  // Declarations
  int is_small;
  static struct FluidState *S_eps; // Need this since prim_to_flux and mhd_calc require a FluidState object as parameter
  double prims[NVAR];
  double prims_eps[NVAR];

  double residual[NFVAR]     = {0};
  double residual_eps[NFVAR] = {0};

  // Initialize prims and residual vectors
  static int firstc = 1;
  if (firstc){
    S_eps = calloc(1, sizeof(struct FluidState));
    firstc = 0;
  }

  PLOOP {
    prims[ip] = S_solver->P[ip][k][j][i];
    prims_eps[ip] = S_solver->P[ip][k][j][i];
    S_eps->P[ip][k][j][i] = prims_eps[ip];
  }
  
  // Calculate residual for Sf->P
  residual_calc(G, S_solver, Si, Ss, U_old, divF, sources_explicit, sources_implicit_old, dt, i, j, k, residual);

  // Numerically evaluate the Jacobian
  for (int col = 0; col < NFVAR; col++) {

    // Evaluate small(P)
    is_small = 0;
    if (fabs(prims[col]) < 0.5 * jacobian_eps) is_small = 1;

    // Compute P_eps
    prims_eps[col] = prims[col] + jacobian_eps * prims[col] * (1 - is_small) + jacobian_eps * is_small;
    S_eps->P[col][k][j][i] = prims_eps[col];

    // Compute the residual for P_eps
    residual_calc(G, S_eps, Si, Ss, U_old, divF, sources_explicit, sources_implicit_old, dt, i, j, k, residual_eps);

    for (int row = 0; row < NFVAR; row++) {
      J[NFVAR*row + col] = (residual_eps[row] - residual[row]) / (prims_eps[col] - prims[col]);
    } // END of row loop

    // Reset P_eps
    prims_eps[col] = prims[col];
    S_eps->P[col][k][j][i] = prims_eps[col];

  } // END OF col loop
}    

// Linear solve (per zone)
void solve(double A[NFVAR*NFVAR], double b[NFVAR]){

  static int ipiv[NFVAR];

  LAPACKE_dgesv(LAPACK_ROW_MAJOR, NFVAR, 1, A, NFVAR, ipiv, b, 1);
}

#endif
