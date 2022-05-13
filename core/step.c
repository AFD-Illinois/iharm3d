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

void step(struct GridGeom *G, struct FluidState *S)
{
  static struct FluidState *Stmp;
  static struct FluidState *Ssave;

  static int first_call = 1;
  if (first_call) {
    Stmp = calloc(1,sizeof(struct FluidState));
    Ssave = calloc(1,sizeof(struct FluidState));
    first_call = 0;
  }

  // Save radial ghost zone values to Stmp->P_BOUND if Dirichlet bc
#if X1L_BOUND == DIRICHLET
  ISLOOP(-NG, -1)
    PLOOP 
      Stmp->P_BOUND[ip][i] = S->P_BOUND[ip][i];
#endif

#if X1R_BOUND == DIRICHLET
  ISLOOP(N1, N1 - 1 + NG)
    PLOOP 
      Stmp->P_BOUND[ip][i-N1] = S->P_BOUND[ip][i-N1];
#endif

  // Need both P_n and P_n+1 to calculate current
  // Work around ICC 18.0.2 bug in assigning to pointers to structs
  // TODO use pointer tricks to avoid deep copy on both compilers
#if INTEL_WORKAROUND
  memcpy(&(Ssave->P),&(S->P),sizeof(GridPrim));
#else
#pragma omp parallel for simd collapse(3)
  PLOOP ZLOOPALL Ssave->P[ip][k][j][i] = S->P[ip][k][j][i];
#endif
  LOGN("Step %d",nstep);
  FLAG("Start step");
  // TODO add back well-named flags /after/ events

  #if DEBUG_GRIM
  fprintf(stdout, "\n----------PREDICTOR STEP----------\n");
  fprintf(stdout, "Fluid state before predictor step: \n");
  fprintf(stdout, "q_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", S->P[Q_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\ndP_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", S->P[DELTA_P_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  #endif

  // Predictor setup
  advance_fluid(G, S, S, Stmp, 0.5*dt);
  FLAG("Advance Fluid Tmp");

#if ELECTRONS
  heat_electrons(G, S, Stmp);
  FLAG("Heat Electrons Tmp");
#endif

  // Fixup routines: smooth over outlier zones
  fixup(G, Stmp);
  FLAG("Fixup Tmp");
#if ELECTRONS
  fixup_electrons(Stmp);
  FLAG("Fixup e- Tmp");
#endif
  // Need an MPI call _before_ fixup_utop to obtain correct pflags
  set_bounds(G, Stmp);
  FLAG("First bounds Tmp");
  fixup_utoprim(G, Stmp);
  FLAG("Fixup U_to_P Tmp");
  set_bounds(G, Stmp);
  FLAG("Second bounds Tmp");

  #if DEBUG_GRIM
  fprintf(stdout, "\n----------CORRECTOR STEP----------\n");
  fprintf(stdout, "Fluid state before corrector step (NOTE: This means boundary syncs have been applied): \n");
  fprintf(stdout, "q_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", S->P[Q_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\ndP_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", S->P[DELTA_P_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  #endif

  // Corrector step
  double ndt = advance_fluid(G, S, Stmp, S, dt);
  FLAG("Advance Fluid Full");

#if ELECTRONS
  heat_electrons(G, Stmp, S);
  FLAG("Heat Electrons Full");
#endif

  fixup(G, S);
  FLAG("Fixup Full");
#if ELECTRONS
  fixup_electrons(S);
  FLAG("Fixup e- Full");
#endif
  set_bounds(G, S);
  FLAG("First bounds Full");
  fixup_utoprim(G, S);
  FLAG("Fixup U_to_P Full");
  set_bounds(G, S);
  FLAG("Second bounds Full");

  #if DEBUG_GRIM
  fprintf(stdout, "\nFluid state after corrector step and boundary syncs: \n");
  fprintf(stdout, "q_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", S->P[Q_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\ndP_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", S->P[DELTA_P_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  #endif

  // Increment time
  t += dt;

  // If we're dumping this step, update the current
  if (t >= tdump) {
    current_calc(G, S, Ssave, dt);
  }

  // New dt proxy to choose fluid or light timestep
  double max_dt = 0, fake_dt = 0;
#if STATIC_TIMESTEP
  if(DEBUG) fake_dt = mpi_min(ndt);
  max_dt = cour*dt_light;
#else
  if(DEBUG) fake_dt = cour*dt_light;
  max_dt = mpi_min(ndt);
#endif

  // Set next timestep
  if (max_dt > SAFE * dt) {
    dt = SAFE * dt;
  } else {
    dt = max_dt;
  }

  LOGN("dt would have been %f",fake_dt);
  LOGN("Instead it is %f",dt);

}

inline double advance_fluid(struct GridGeom *G, struct FluidState *Si,
  struct FluidState *Ss, struct FluidState *Sf, double Dt)
{
  static GridPrim *dU;
  static struct FluidFlux *F;
  #if GRIM_TIMESTEPPER
  // Temporary fluid struct for nonlinear solver
  static struct FluidState *S_solver;
  #endif

  static int firstc = 1;
  if (firstc) {
    dU = calloc(1,sizeof(GridPrim));
    F = calloc(1,sizeof(struct FluidFlux));
    #if GRIM_TIMESTEPPER
    S_solver = calloc(1, sizeof(struct FluidState));
    #endif
    firstc = 0;
  }

  // Work around ICC 18.0.2 bug in assigning to pointers to structs
#if INTEL_WORKAROUND
  memcpy(&(Sf->P),&(Si->P),sizeof(GridPrim));
#else
#pragma omp parallel for simd collapse(3)
  PLOOP ZLOOPALL Sf->P[ip][k][j][i] = Si->P[ip][k][j][i];
#endif

	double ndt = get_flux(G, Ss, F);	

#if METRIC == MKS
  fix_flux(F);
#endif

  //Constrained transport for B
  flux_ct(F);

  // Flux diagnostic globals
  diag_flux(F);

// GRIM vs HARM time-stepper
#if GRIM_TIMESTEPPER
  
  // Set zero pflags and fail_save to zero
  zero_arrays();

  // Obtain Si->U
  get_state_vec(G, Si, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, Si, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Si->U);

  // Obtain state for Ss and Ss->U (needed for source terms)
  get_state_vec(G, Ss, CENT, -NG, N3 + NG - 1, -NG, N2 + NG - 1, -NG, N1 + NG - 1);
  prim_to_flux_vec(G, Ss, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Ss->U);

  #if DEBUG_GRIM
  fprintf(stdout, "Si (in advance fluid): \n");
  fprintf(stdout, "q:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", Si->q[k][j][i]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\ndP:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", Si->delta_p[k][j][i]);
      fprintf(stdout, "\n");
    }
  #endif

  #if DEBUG_GRIM
  fprintf(stdout, "Ss (in advance fluid): \n");
  fprintf(stdout, "q:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", Ss->q[k][j][i]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\ndP:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", Ss->delta_p[k][j][i]);
      fprintf(stdout, "\n");
    }
  #endif

  // Initial guess for S_solver->P
  #pragma omp parallel for simd collapse(3)
  PLOOP ZLOOP S_solver->P[ip][k][j][i] = Ss->P[ip][k][j][i];

  // time-step by root-finding the residual
  grim_timestep(G, Si, Ss, S_solver, F, Dt);

  // compute new conserved variables
  #pragma omp parallel for simd collapse(3)
  PLOOP ZLOOP Sf->P[ip][k][j][i] = S_solver->P[ip][k][j][i];
  get_state_vec(G, Sf, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, Sf, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Sf->U);

  #if DEBUG_GRIM
  fprintf(stdout, "Sf (in advance fluid, after timestep): \n");
  fprintf(stdout, "q_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", Sf->P[Q_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\ndP_tilde:\n");
  KLOOP_DEBUG_GRIM
    JLOOP_DEBUG_GRIM {
      ILOOP_DEBUG_GRIM
        fprintf(stdout, "%g ", Sf->P[DELTA_P_TILDE][k][j][i]);
      fprintf(stdout, "\n");
    }
  #endif
  
  // update failures
  // NOTE: These are no longer U_to_P failures but rather zones where convergence was not achieved for the nonlinear solver
  #pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    fail_save[k][j][i] = pflag[k][j][i];
  }
  
#else
  // Update Si to Sf
  timer_start(TIMER_UPDATE_U);
  get_state_vec(G, Ss, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  get_ideal_fluid_source_vec(G, Ss, dU);

  get_state_vec(G, Si, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, Si, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, Si->U);

#pragma omp parallel for collapse(3)
  PLOOP ZLOOP {
    Sf->U[ip][k][j][i] = Si->U[ip][k][j][i] +
      Dt*((F->X1[ip][k][j][i] - F->X1[ip][k][j][i+1])/dx[1] +
          (F->X2[ip][k][j][i] - F->X2[ip][k][j+1][i])/dx[2] +
          (F->X3[ip][k][j][i] - F->X3[ip][k+1][j][i])/dx[3] +
          (*dU)[ip][k][j][i]);
  }
  timer_stop(TIMER_UPDATE_U);

  timer_start(TIMER_U_TO_P);
#pragma omp parallel for collapse(3)
  ZLOOP {
    pflag[k][j][i] = U_to_P(G, Sf, i, j, k, CENT);
  }
  timer_stop(TIMER_U_TO_P);

#pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    fail_save[k][j][i] = pflag[k][j][i];
  }

#endif
  return ndt;
}
