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
  //static struct FluidState Stmp;
  if (first_call) {
    Stmp = (struct FluidState*)malloc(sizeof(struct FluidState));
    F = (struct FluidFlux*)malloc(sizeof(struct FluidFlux));
    first_call = 0;
  }
  
  double ndt;

  //printf("\n");
  //PLOOP printf("P[%i][NG][64][64] = %.10e\n", ip, S->P[ip][NG][64][64]);
  //printf("P[RHO][3][64][6] = %.10e\n", S->P[RHO][3][4][6]);
  //printf("P[U1][3][64][6] = %.10e\n", S->P[U1][3][4][6]);
  /*printf("P[RHO][3][3][6] = %.10e\n", S->P[RHO][3][3][6]);
  printf("P[U1][3][3][6] = %.10e\n", S->P[U1][3][3][6]);
  printf("P[RHO][3][5][6] = %.10e\n", S->P[RHO][3][5][6]);
  printf("P[U1][3][5][6] = %.10e\n", S->P[U1][3][5][6]);
  printf("P[RHO][3][4][5] = %.10e\n", S->P[RHO][3][4][5]);
  printf("P[U1][3][4][5] = %.10e\n", S->P[U1][3][4][5]);
  printf("P[RHO][3][4][7] = %.10e\n", S->P[RHO][3][4][7]);
  printf("P[U1][3][4][7] = %.10e\n", S->P[U1][3][4][7]);*/
  //if (nstep > 10) exit(-1);

  // Predictor setup
  //ndt = advance(P, P, 0.5*dt, Ph, 0);
  ndt = advance_fluid(G, S, S, Stmp, 0.5*dt);
  /*printf("past advance\n");
  printf("P[RHO][3][4][6] = %.10e\n", Stmp->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", Stmp->P[U1][3][4][6]);*/
  //exit(-1);
  #if ELECTRONS
  heat_electrons(P, Ph, 0.5*dt);
  #endif
  fixup(G, Stmp);
  /*printf("past base fixup\n");
  printf("P[RHO][3][4][6] = %.10e\n", Stmp->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", Stmp->P[U1][3][4][6]);*/
  fixup_utoprim(G, Stmp);
  /*printf("past fixup\n");
  printf("P[RHO][3][4][6] = %.10e\n", Stmp->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", Stmp->P[U1][3][4][6]);*/
  #if ELECTRONS
  fixup_electrons(Ph);
  #endif
  set_bounds(G, Stmp);
  /*printf("past bounds\n");
  printf("P[RHO][3][4][6] = %.10e\n", Stmp->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", Stmp->P[U1][3][4][6]);*/

  // Radiation step
  #if RADIATION
  make_superphotons(Ph, t, dt);
  push_superphotons(t, dt);
  #endif

  // Corrector step
  //ndt = advance(P, Ph, dt, P, 1);
  ndt = advance_fluid(G, S, Stmp, S, dt);
  /*printf("past advance\n");
  printf("P[RHO][3][4][6] = %.10e\n", S->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", S->P[U1][3][4][6]);*/
  #if ELECTRONS
  heat_electrons(Ph, P, dt);
  #endif
  fixup(G, S);
  fixup_utoprim(G, S);
  /*printf("past fixup\n");
  printf("P[RHO][3][4][6] = %.10e\n", S->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", S->P[U1][3][4][6]);*/
  #if ELECTRONS
  fixup_electrons(P);
  #endif
  set_bounds(G, S);
  /*printf("past bounds\n");
  printf("P[RHO][3][4][6] = %.10e\n", S->P[RHO][3][4][6]);
  printf("P[U1][3][4][6] = %.10e\n", S->P[U1][3][4][6]);*/

  #if RADIATION
  // Apply radiation four-force to fluid
  memset((void*)&radG[0][0][0][0], 0,
    (N1+2*NG)*(N2+2*NG)*(N3+2*NG)*NDIM*sizeof(double));
  
  // Control Monte Carlo resolution
  #endif

  // Increment time
  t += dt;

  // Set next timestep
  if (ndt > SAFE * dt) {
    ndt = SAFE * dt;
  }
  dt = ndt;
  dt = mpi_min(dt);
}

double advance_fluid(struct GridGeom *G, struct FluidState *Si, 
  struct FluidState *Ss, struct FluidState *Sf, double Dt)
{
  static GridPrim dU;

  #pragma omp parallel for collapse(3)
  PLOOP ZLOOP Sf->P[ip][k][j][i] = Si->P[ip][k][j][i];
  //memcpy(Sf->P, Si->P, NVAR*(N3+2*NG)*(N2+2*NG)*(N1+2*NG)*sizeof(double));
  //memcpy(Sf->P, Si->P, sizeof(GridPrim);
  
  //ZLOOP PLOOP Pf[i][j][k][ip] = Pi[i][j][k][ip];
  
  double ndt = get_flux(G, Ss, F);

  #if METRIC == MKS
  fix_flux(F);
  #endif

  flux_ct(F);

  // Evaluate diagnostics based on fluxes
  //timer_start(TIMER_DIAG);
  //diag_flux(F1, F2, F3);
  //timer_stop(TIMER_DIAG);

  // Update Si to Sf
  timer_start(TIMER_UPDATE_U);
  #pragma omp parallel for collapse(3)
  ZLOOP {
    get_fluid_source(G, Ss, i, j, k, dU);
  }

  // Can remove this later after appropriate initialization
  #pragma omp parallel for collapse(3)
  ZLOOP {
    get_state(G, Si, i, j, k, CENT);
    /*if (i == 6 && j == 4 && k == 3) {
      printf("ucon[] = %e %e %e %e gdet = %e\n",
      Si->ucon[0][k][j][i], Si->ucon[1][k][j][i], Si->ucon[2][k][j][i], Si->ucon[3][k][j][i],
      G->gdet[CENT][j][i]);
    }*/
  }

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
  #pragma omp parallel for collapse(3)
  ZLOOP {
    get_state(G, Sf, i, j, k, CENT);
  }

  #pragma omp parallel for collapse(3)
  ZLOOP {
    fail_save[k][j][i] = pflag[k][j][i];
  }
  //timer_stop(TIMER_UPDATE);

  return ndt;
}

