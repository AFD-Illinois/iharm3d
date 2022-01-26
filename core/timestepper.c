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

// Function declaration(s): residual calculation, L2 norm
void residual_calc(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidState *Sf, struct FluidFlux *F, double residual[NVAR], double dt, int i, int j, int k);

void residual_calc_vec(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidState *Sf, struct FluidFlux *F, GridPrim *residual, double dt);

int l2norm(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidFlux *F, double delta[NVAR], double dt, int i, int j, int k);


void grim_timestep(struct GridGeom *G, struct FluidState *Si, 
  struct FluidState *Ss, struct FluidState *Sf, struct FluidFlux *F, double dt)
{

  // Declarations and memory allocation
  static struct FluidState *S_eps;
  static GridPrim *residual;
  static GridPrimMatrix *jacobian;
  static double residual_eps[NVAR];

  static int firstc = 1;
  if (firstc) {
    S_eps = calloc(1, sizeof(struct FluidState));
    residual = calloc(1, sizeof(GridPrim));
    jacobian = calloc(1, sizeof(GridPrimMatrix));

    firstc = 0;
  }

  // Initial guess for Sf->P (right now it is set to Ss->P)
  #if INTEL_WORKAROUND
    memcpy(&(Sf->P), &(Ss->P), sizeof(GridPrim));
  #else
  #pragma omp parallel for simd collapse(3)
    PLOOP ZLOOPALL Sf->P[ip][k][j][i] = Ss->P[ip][k][j][i];
  #endif

  track_solver_iterations = 0;

  while (track_solver_iterations < MAX_ITER_SOLVER){

    // Compute conserved variables for U^(n+1/2)
    get_state_vec(G, Sf, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
    prim_to_flux_vec(G, Sf, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Sf->U);

    // Compute residual for Sf->P
    get_state_vec(G, Si, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
    residual_calc_vec(G, Si, Ss, Sf, F, residual, dt);

    // Initialize P_eps^(n+1/2)
    #if INTEL_WORKAROUND
      memcpy(&(S_eps->P), &(Sf->P), sizeof(GridPrim));
    #else
    #pragma omp parallel for simd collapse(3)
      PLOOP ZLOOPALL S_eps->P[ip][k][j][i] = Sf->P[ip][k][j][i];
    #endif

    #pragma omp parallel for simd collapse(2)
    ZLOOP {
      // Only unconverged zones
      if (!pflag[k][j][i]){
        // Numerically evaluate the Jacobian
        for (int col = 0; col < NVAR; col++) {
          // Find small(P)
          static int is_small = 0;
          if (abs(Sf->P[col][k][j][i]) < 0.5 * JACOBIAN_EPS) {
            is_small = 1;
          }

          // Compute P_eps
          S_eps->P[col][k][j][i] = Sf->P[col][k][j][i] + JACOBIAN_EPS*Sf->P[col][k][j][i]*(1 - is_small) + JACOBIAN_EPS*is_small;

          // Compute residual for P_eps
          get_state(G, S_eps, i, j, k, CENT);
          prim_to_flux(G, S_eps, i, j, k, 0, CENT, S_eps->U);
          residual_calc(G, Si, Ss, S_eps, F, residual_eps, dt, i, j, k);

          for (int row = 0; row < NVAR; row++) {
            // Jacobian computation
             (*jacobian)[NVAR*row + col][k][j][i] = (residual_eps[row] - (*residual)[row][k][j][i])
             / (S_eps->P[col][k][j][i] - Sf->P[col][k][j][i]);
          }

          // Reset P_eps
          S_eps->P[col][k][j][i] = Sf->P[col][k][j][i];
        }

        // Linear solve
        // Linear solve identifiers
        double *A_per_zone = (double*)calloc(NVAR*NVAR, sizeof(double));
        double *b_per_zone = (double*)calloc(NVAR, sizeof(double));
        int ipiv[NVAR];

        // Set zone-wise A, b
        for (int col = 0; col < NVAR; col++){
          for (int row = 0; row < NVAR; row++){
            A_per_zone[NVAR*row + col] = (*jacobian)[NVAR*row+col][k][j][i];
          }
          b_per_zone[col] = -*residual[col][k][j][i];
        }

        // Linear solve (zone-wise)
        LAPACKE_dgesv(LAPACK_ROW_MAJOR, NVAR, 1, A_per_zone, NVAR, ipiv, b_per_zone, 1);

        // update pflags and Sf->P if root finding tolerance is met
        pflag[k][j][i] = l2norm(G, Si, Ss, F, b_per_zone, dt, i, j, k);
        PLOOP Sf->P[ip][k][j][i] = Ss->P[ip][k][j][i] + b_per_zone[ip];
      } // pflag condition
    } // ZLOOP
    track_solver_iterations += 1;
  } // Solver while loop
}

void residual_calc(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidState *Sf, struct FluidFlux *F, double residual[NVAR], double dt, int i, int j, int k) {

  // Get fluid source
  static double dU[NVAR];
  
  get_fluid_source(G, Ss, dU, i, j, k);

// Compute residual
  PLOOP {
    residual[ip] = (Sf->U[ip][k][j][i] - Si->U[ip][k][j][i]) / dt
      + (F->X1[ip][k][j][i+1] - F->X1[ip][k][j][i])/dx[1]
      + (F->X2[ip][k][j+1][i] - F->X2[ip][k][j][i])/dx[2]
      + (F->X3[ip][k+1][j][i] - F->X3[ip][k][j][i])/dx[3]
      - dU[ip];
  }
}

void residual_calc_vec(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidState *Sf, struct FluidFlux *F, GridPrim *residual, double dt) {

  // Get fluid source
  static GridPrim *dU;
  static int firstc = 1;
  if (firstc) {
    dU = calloc(1, sizeof(GridPrim));
    firstc = 0;
  }
  get_fluid_source_vec(G, Ss, dU);

// Compute residual
#pragma omp parallel for collapse(3)
  PLOOP ZLOOP {
    (*residual)[ip][k][j][i] = (Sf->U[ip][k][j][i] - Si->U[ip][k][j][i]) / dt
      + (F->X1[ip][k][j][i+1] - F->X1[ip][k][j][i])/dx[1]
      + (F->X2[ip][k][j+1][i] - F->X2[ip][k][j][i])/dx[2]
      + (F->X3[ip][k+1][j][i] - F->X3[ip][k][j][i])/dx[3]
      - (*dU)[ip][k][j][i];
  }
}

// Compute L2 norm and update Sf->P if tolerance is met
int l2norm(struct GridGeom *G, struct FluidState *Si, struct FluidState *Ss,
  struct FluidFlux *F, double delta[NVAR], double dt, int i, int j, int k) {
  // Declarations
  double norm = 0;
  static struct FluidState *S_updated;
  static double residual_updated[NVAR];

  static int firstc = 1;
  if (firstc) {
    S_updated = calloc(1, sizeof(struct FluidState));

    firstc = 0;
  }

  // Updated primitives
  PLOOP {
    S_updated->P[ip][k][j][i] = Ss->P[ip][k][j][i] + delta[ip]; // will add linesearch parameter
  }

  // Compute residual
  get_state(G, S_updated, i, j, k, CENT);
  prim_to_flux(G, S_updated, i, j, k, 0, CENT, S_updated->U);
  residual_calc(G, Si, Ss, S_updated, F, residual_updated, dt, i, j, k);

  // L2 norm
  PLOOP{
    norm += pow(residual_updated[ip], 2);
  }
  norm = sqrt(norm);

  // Check against tolerance
  if (norm < ROOTFIND_TOL) return 1;
  else return 0;
}

#endif