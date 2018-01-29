/******************************************************************************
 *                                                                            *
 * MAIN.C                                                                     *
 *                                                                            *
 * ESTABLISH MPI COMMUNICATION, LOOP OVER TIME, COMPLETE                      *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "defs.h"
#include <time.h>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
  //omp_set_num_threads(1);

  // Check for minimal required MPI thread support
  int threadSafety;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &threadSafety);
  if (threadSafety < MPI_THREAD_FUNNELED) {
    fprintf(stderr, "Thread support < MPI_THREAD_FUNNELED. Unsafe.\n");
    exit(1);
  }

  init_mpi();

  if (mpi_io_proc()) {
    fprintf(stdout, "\n          ************************************************************\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          BHLIGHT                         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *          Ryan, Dolence & Gammie ApJ 807:31, 2015         *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *  B R Ryan                                                *\n");
    fprintf(stdout, "          *  J C Dolence                                             *\n");
    fprintf(stdout, "          *  C F Gammie                                              *\n");
    fprintf(stdout, "          *  S M Ressler                                             *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *                          SYNTAX                          *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          *    -o /path/to/output/directory/                         *\n");
    fprintf(stdout, "          *    -m [black hole mass in solar masses]                  *\n");
    fprintf(stdout, "          *    -M [mass unit in grams]                               *\n");
    fprintf(stdout, "          *                                                          *\n");
    fprintf(stdout, "          ************************************************************\n\n");
  }


  // Read command line arguments

  #if RADIATION
  mbh = -1.;
  M_unit = -1.;
  #endif
  for (int n = 0; n < argc; n++) {
    // Check for argv[n] of the form '-*'
    if (*argv[n] == '-' && *(argv[n]+1) != '\0' && *(argv[n]+2) == '\0' &&
        n < argc-1) {
      if (*(argv[n]+1) == 'o') { // Set output directory path
        strcpy(outputdir, argv[++n]);

        if( chdir(outputdir) != 0) {
            fprintf(stderr, "Output directory does not exist!\n");
            exit(2);
        }
      }
      #if RADIATION
      else if (*(argv[n]+1) == 'm') { // Set black hole mass in solar masses
        mbh = strtof(argv[++n], NULL)*MSUN;
      } else if (*(argv[n]+1) == 'M') { // Set mass unit in grams
        M_unit = strtof(argv[++n], NULL);
      }
      #endif
    }
  }

  strcpy(dumpdir, "dumps/");
  strcpy(restartdir, "restarts/");
  if (mpi_io_proc()) {
    mkdir(dumpdir, 0777);
    mkdir(restartdir, 0777);
  }

  // Sanity checks
  if ((RECONSTRUCTION == PPM || RECONSTRUCTION == WENO || RECONSTRUCTION == MP5)
      && NG < 3) {
    fprintf(stderr, "not enough ghost zones!\n") ;
    fprintf(stderr, "PPM/WENO/MP5 + NG < 3\n") ;
    exit(1);
  }
  #if (RADIATION && METRIC == MKS)
  if (mbh < 0. || M_unit < 0.) {
    fprintf(stderr, "Error! Bad parameters! mbh = %e and M_unit = %e\n",
      mbh, M_unit);
    exit(-1);
  }
  #endif

  #pragma omp parallel
  {
    #pragma omp master
    {
      nthreads = omp_get_num_threads();
    } // omp master
  } // omp parallel

  // Initialize global variables and arrays
  reset_log_variables();
  nstep = 0;
  struct GridGeom *G = (struct GridGeom*)malloc(sizeof(struct GridGeom));
  struct FluidState *S = (struct FluidState*)malloc(sizeof(struct FluidState));

  // Perform initializations, either directly or via checkpoint
  init_random(1); // NEED RANDOM SEED BASED ON TIME() HERE
  is_restart = restart_init(G, S);
  if (!is_restart) {
    time_init();
    init(G, S);
    #if RADIATION
    init_rad(P);
    #endif
    tdump = DTd;
    tlog  = DTl;
    if (mpi_io_proc())
      fprintf(stdout, "Initial conditions generated\n\n");
  }

  // Initial diagnostics
  diag(G, S, DIAG_INIT);

  if (mpi_io_proc())
    fprintf(stdout, "t = %e tf = %e\n", t, tf);
  //if (!is_restart)
  //  diag(DIAG_DUMP);

/*******************************************************************************
    MAIN LOOP
*******************************************************************************/
  if (mpi_io_proc())
    fprintf(stdout, "\nEntering main loop\n");
  int dumpThisStep = 0;
  while (t < tf) {
    dumpThisStep = 0;
    timer_start(TIMER_ALL);

    // Step variables forward in time
    step(G, S);
    nstep++;

    // Don't step beyond end of run
    if (t + dt > tf) {
      dt = tf - t;
    }

    if (mpi_io_proc()) {
      fprintf(stdout, "t = %10.5g dt = %10.5g n = %8d\n", t, dt, nstep);
      #if RADIATION
      fprintf(stdout, "made = %d abs = %d scatt = %d lost = %d tot = %d\n",
        step_made, step_abs, step_scatt, step_lost, step_tot);
      #endif
    }

    // File I/O with set frequencies
    if (t < tf) {
      if (t >= tdump) {
        dumpThisStep = 1;
        diag(G, S, DIAG_DUMP);
        tdump += DTd;
      }
      if (t >= tlog) {
        diag(G, S, DIAG_LOG);
        tlog += DTl;
      }
      if (nstep % DTr == 0)
        restart_write(S);
    }

    timer_stop(TIMER_ALL);

    if (nstep % DTp == 0)
      report_performance();

    reset_log_variables();
  }
/*******************************************************************************
    END MAIN LOOP
*******************************************************************************/
  if (dumpThisStep == 0) diag(G, S, DIAG_FINAL);

  return 0;
}
