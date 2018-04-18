/******************************************************************************
 *                                                                            *
 * STEP.C                                                                     *
 *                                                                            *
 * ADVANCES SIMULATION BY ONE TIMESTEP                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Declarations
double advance_fluid(struct GridGeom *G, struct FluidState *Si,
  struct FluidState *Ss, struct FluidState *Sf, double Dt);
static struct FluidState *Stmp;
static struct FluidFlux *F;
static int first_call = 1;

void step(struct GridGeom *G, struct FluidState *S)
{
  if (first_call) {
    Stmp = (struct FluidState*)calloc(1,sizeof(struct FluidState));
    F = (struct FluidFlux*)calloc(1,sizeof(struct FluidFlux));
    first_call = 0;
  }

  // Predictor setup
  advance_fluid(G, S, S, Stmp, 0.5*dt);

  #if ELECTRONS
  heat_electrons(P, Ph, 0.5*dt);
  #endif

  fixup(G, Stmp);

  fixup_utoprim(G, Stmp);

  #if ELECTRONS
  fixup_electrons(Ph);
  #endif

  set_bounds(G, Stmp);

  // Radiation step
  #if RADIATION
  make_superphotons(Ph, t, dt);
  push_superphotons(t, dt);
  #endif

  // Corrector step
  double ndt = advance_fluid(G, S, Stmp, S, dt);

  #if ELECTRONS
  heat_electrons(Ph, P, dt);
  #endif

  fixup(G, S);
  fixup_utoprim(G, S);

  #if ELECTRONS
  fixup_electrons(P);
  #endif

  set_bounds(G, S);


  #if RADIATION
  // Apply radiation four-force to fluid
  memset((void*)&radG[0][0][0][0], 0,
    (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*NDIM*sizeof(double));
  #endif

  // Increment time
  t += dt;

  // Set next timestep
  if (ndt > SAFE * dt) {
    ndt = SAFE * dt;
  }
  dt = mpi_min(ndt);
}

double advance_fluid(struct GridGeom *G, struct FluidState *Si,
  struct FluidState *Ss, struct FluidState *Sf, double Dt)
{
  static GridPrim dU;

  // TODO this crashes ICC. Why.
//  #pragma omp parallel for collapse(3)
//  PLOOP ZLOOPALL Sf->P[ip][k][j][i] = Si->P[ip][k][j][i];

  // In the meantime...
  memcpy(&(Sf->P),&(Si->P),sizeof(GridPrim));


  double ndt = get_flux(G, Ss, F);

  #if METRIC == MKS
  fix_flux(F);
  #endif

  //Constrained transport for B
  flux_ct(F);

  // Flux diagnostic globals
  diag_flux(F);

  // Update Si to Sf
  timer_start(TIMER_UPDATE_U);
  #pragma omp parallel for collapse(3)
  ZLOOP {
    get_fluid_source(G, Ss, i, j, k, dU);
  }

  // Can remove this later after appropriate initialization
//  #pragma omp parallel for collapse(3)
//  ZLOOP {
//    get_state(G, Si, i, j, k, CENT);
//  }

  get_state_vec(G, Si, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);

  prim_to_flux_vec(G, Si, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Si->U);

  PLOOP {
    #pragma omp parallel for collapse(3)
    ZLOOP {
      Sf->U[ip][k][j][i] = Si->U[ip][k][j][i] +
        Dt*((F->X1[ip][k][j][i] - F->X1[ip][k][j][i+1])/dx[1] +
            (F->X2[ip][k][j][i] - F->X2[ip][k][j+1][i])/dx[2] +
            (F->X3[ip][k][j][i] - F->X3[ip][k+1][j][i])/dx[3] +
            dU[ip][k][j][i]);
    }
  }
  timer_stop(TIMER_UPDATE_U);

  timer_start(TIMER_U_TO_P);
  #pragma omp parallel for collapse(3)
  ZLOOP {
    pflag[k][j][i] = U_to_P(G, Sf, i, j, k, CENT);
  }
  timer_stop(TIMER_U_TO_P);

  // Not complete without setting four-vectors
//  #pragma omp parallel for collapse(3)
//  ZLOOP {
//    get_state(G, Sf, i, j, k, CENT);
//  }
  get_state_vec(G, Sf, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);

  #pragma omp parallel for collapse(3)
  ZLOOP {
    fail_save[k][j][i] = pflag[k][j][i];
  }
  //timer_stop(TIMER_UPDATE);

  return ndt;
}
