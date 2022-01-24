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

// Function declaration(s): residual calculation
void residual_calc(struct GridGeom *G, struct FluidState *Si,
  struct FluidState *Sf, struct FluidFlux *F, GridPrim *residual, double dt);


void grim_timestep(struct GridGeom *G, struct FluidState *Si, 
  struct FluidState *Ss, struct FluidState *Sf, struct FluidFlux *F, double dt)
{

  // Declarations and memory allocation
  static struct FluidState *S_eps;
  static GridPrim *residual;
  static GridPrim *residual_eps;
  static GridPrimMatrix *jacobian;

  // Linear solve identifiers
  double *A_per_zone = calloc(NVAR*NVAR, sizeof(double));
  double *b_per_zone = calloc(NVAR, sizeof(double));
  static int ipiv[NVAR];

  static int firstc = 1;
  if (firstc) {
    S_eps = calloc(1, sizeof(struct FluidState));
    residual = calloc(1, sizeof(GridPrim));
    residual_eps = calloc(1, sizeof(GridPrim));
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

  // LOOP FOR ITERATIVE SOLVER //

  // Compute conserved variables for U^(n+1/2)
  get_state_vec(G, Sf, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, Sf, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Sf->U);

  // Compute residual for Sf->P
  get_state_vec(G, Si, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  residual_calc(G, Si, Sf, F, residual, dt);

  // Initialize P_EPS^(n+1/2)
  #if INTEL_WORKAROUND
    memcpy(&(S_eps->P), &(Sf->P), sizeof(GridPrim));
  #else
  #pragma omp parallel for simd collapse(3)
    PLOOP ZLOOPALL S_eps->P[ip][k][j][i] = Sf->P[ip][k][j][i];
  #endif

  ZLOOP {
    // Numerically evaluate the Jacobian
    for (int col = 0; col < NVAR; col++) {
      // Find small(P)
      static int is_small = 0;
      if (abs(Sf->P[col][k][j][i]) < 0.5 * JACOBIAN_EPS) {
        is_small = 1;
      }

      // Compute P_eps
      S_eps->P[col][k][j][i] = Sf->P[col][k][j][i] + JACOBIAN_EPS*Sf->P[col][k][j][i]*(1 - is_small) + JACOBIAN_EPS*is_small;

      // Compute residual for prims_eps
      prim_to_flux(G, S_eps, i, j, k, 0, CENT, S_eps->U);
      residual_calc(G, Si, S_eps, F, residual_eps, dt);

      for (int row = 0; row < NVAR; row++) {
        // Jacobian update
        *jacobian[NVAR*row + col][k][j][i] = (*residual_eps[row][k][j][i] - *residual[row][k][j][i]) 
        / (S_eps->P[col][k][j][i] - Sf->P[col][k][j][i]);
      }

      // Reset P_eps
      #if INTEL_WORKAROUND
        memcpy(&(S_eps->P), &(Sf->P), sizeof(GridPrim));
      #else
      #pragma omp parallel for simd collapse(3)
        PLOOP S_eps->P[ip][k][j][i] = Sf->P[ip][k][j][i];
      #endif
    }
  }

  // Linear solve
  #pragma omp parallel for collapse(3)
  ZLOOP{
    // Linear solve identifiers
    double *A_per_zone = calloc(NVAR*NVAR, sizeof(double));
    double *b_per_zone = calloc(NVAR, sizeof(double));
    int ipiv[NVAR];

    // Set zone-wise A, b
    for (int col = 0; col < NVAR; col++){
      for (int row = 0; row < NVAR; row++){
        A_per_zone[NVAR*row + col] = *jacobian[NVAR*row + col][k][j][i];
      }
      b_per_zone[col] = -*residual[col][k][j][i];
    }

    // Linear solve (zone-wise)
    LAPACKE_dgesv(LAPACK_ROW_MAJOR, NVAR, 1, A_per_zone, NVAR, ipiv, b_per_zone, 1);
  }

  // L2 NORM CHECK

}

void residual_calc(struct GridGeom *G, struct FluidState *Si,
  struct FluidState *Sf, struct FluidFlux *F, GridPrim *residual, double dt) {

  // Get fluid source
  static GridPrim *dU;
  static int firstc = 1;
  if (firstc) {
    dU = calloc(1, sizeof(GridPrim));
    firstc = 0;
  }
  get_fluid_source(G, Si, dU);

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

#endif