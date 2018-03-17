/******************************************************************************
 *                                                                            *
 * TIMING.C                                                                   *
 *                                                                            *
 * PERFORMANCE TIMING AND REPORTING                                           *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static double timers[NUM_TIMERS];
static double times[NUM_TIMERS];

void time_init()
{
  for (int n = 0; n < NUM_TIMERS; n++) {
    times[n] = 0.;
  }
}

void timer_start(int timerCode)
{
  #pragma omp barrier

  #pragma omp single
  {
    timers[timerCode] = omp_get_wtime();
  }
}

void timer_stop(int timerCode)
{
  #pragma omp barrier

  #pragma omp single
  {
    times[timerCode] += (omp_get_wtime() - timers[timerCode]);
  }
}

// Report a running average of performance data
void report_performance()
{
  if (mpi_myrank() == 0) {
    fprintf(stdout, "\n********** PERFORMANCE **********\n");
    fprintf(stdout, "   RECON:    %8.4g s (%.4g %%)\n",
      times[TIMER_RECON]/nstep, 100.*times[TIMER_RECON]/times[TIMER_ALL]);
    fprintf(stdout, "   LR_TO_F:  %8.4g s (%.4g %%)\n",
      times[TIMER_LR_TO_F]/nstep, 100.*times[TIMER_LR_TO_F]/times[TIMER_ALL]);
    fprintf(stdout, "   CMAX:     %8.4g s (%.4g %%)\n",
      times[TIMER_CMAX]/nstep, 100.*times[TIMER_CMAX]/times[TIMER_ALL]);
    fprintf(stdout, "   FLUX_CT:  %8.4g s (%.4g %%)\n",
      times[TIMER_FLUX_CT]/nstep, 100.*times[TIMER_FLUX_CT]/times[TIMER_ALL]);
    fprintf(stdout, "   UPDATE_U: %8.4g s (%.4g %%)\n",
      times[TIMER_UPDATE_U]/nstep, 100.*times[TIMER_UPDATE_U]/times[TIMER_ALL]);
    fprintf(stdout, "   U_TO_P:   %8.4g s (%.4g %%)\n",
      times[TIMER_U_TO_P]/nstep, 100.*times[TIMER_U_TO_P]/times[TIMER_ALL]);
    fprintf(stdout, "   FIXUP:    %8.4g s (%.4g %%)\n",
      times[TIMER_FIXUP]/nstep, 100.*times[TIMER_FIXUP]/times[TIMER_ALL]);
    fprintf(stdout, "   BOUND:    %8.4g s (%.4g %%)\n",
      times[TIMER_BOUND]/nstep, 100.*times[TIMER_BOUND]/times[TIMER_ALL]);
    fprintf(stdout, "   BOUND_COMMS:    %8.4g s (%.4g %%)\n",
						times[TIMER_BOUND_COMMS]/nstep, 100.*times[TIMER_BOUND_COMMS]/times[TIMER_ALL]);
    fprintf(stdout, "   DIAG:     %8.4g s (%.4g %%)\n",
      times[TIMER_DIAG]/nstep, 100.*times[TIMER_DIAG]/times[TIMER_ALL]);
    fprintf(stdout, "   LR_STATE:     %8.4g s (%.4g %%)\n",
      times[TIMER_LR_STATE]/nstep, 100.*times[TIMER_LR_STATE]/times[TIMER_ALL]);
    fprintf(stdout, "   LR_PTOF:     %8.4g s (%.4g %%)\n",
      times[TIMER_LR_PTOF]/nstep, 100.*times[TIMER_LR_PTOF]/times[TIMER_ALL]);
    fprintf(stdout, "   LR_VCHAR:     %8.4g s (%.4g %%)\n",
      times[TIMER_LR_VCHAR]/nstep, 100.*times[TIMER_LR_VCHAR]/times[TIMER_ALL]);
    fprintf(stdout, "   LR_CMAX:     %8.4g s (%.4g %%)\n",
      times[TIMER_LR_CMAX]/nstep, 100.*times[TIMER_LR_CMAX]/times[TIMER_ALL]);
    fprintf(stdout, "   LR_FLUX:     %8.4g s (%.4g %%)\n",
      times[TIMER_LR_FLUX]/nstep, 100.*times[TIMER_LR_FLUX]/times[TIMER_ALL]);
    #if ELECTRONS
    fprintf(stdout, "   E_HEAT:   %8.4g s (%.4g %%)\n",
      times[TIMER_ELECTRON_HEAT]/nstep,
      100.*times[TIMER_ELECTRON_HEAT]/times[TIMER_ALL]);
    fprintf(stdout, "   E_FIXUP:  %8.4g s (%.4g %%)\n",
      times[TIMER_ELECTRON_FIXUP]/nstep,
      100.*times[TIMER_ELECTRON_FIXUP]/times[TIMER_ALL]);
    #endif
    fprintf(stdout, "   ALL:      %8.4g s\n", times[TIMER_ALL]/nstep);
    fprintf(stdout, "   ZONE CYCLES PER\n");
    fprintf(stdout, "     CORE-SECOND: %e\n",
      N1TOT*N2TOT*N3TOT/(times[TIMER_ALL]*nthreads/nstep));
    fprintf(stdout, "     NODE-SECOND: %e\n",
      N1TOT*N2TOT*N3TOT/(times[TIMER_ALL]/nstep));
  }
}
