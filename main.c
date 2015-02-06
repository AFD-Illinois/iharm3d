
/* M1: some updates including automatic measurement of
   zone cycles/second 
   
   cfg 12.24.14 */

#include "decs.h"
#include "defs.h"	/* defining declarations */
#include <time.h>
#include "mpi.h"

#define  TIMER_NSTEP	(512)

/*****************************************************************/
/*****************************************************************
   main():
   ------

     -- Initializes, time-steps, and concludes the simulation. 
     -- Handles timing of output routines;
     -- Main is main, what more can you say.  

-*****************************************************************/
int main(int argc, char *argv[])
{
	double tdump, tpdump, timage, tlog, zcps;
	int nfailed = 0;
	double ti_t,t0_t,t1_t,tf_t ;
	void init_mpi();

	MPI_Init(&argc, &argv);

	init_mpi();

	flag_of_convenience=0;

#if ALLOW_COOL
        if (argc > 1) {
                sscanf(argv[1], "%lf", &M_unit);
                fprintf(stderr, "M_unit: %g\n", M_unit);
        } else {
                fprintf(stderr, "usage: harm M_unit\n");
                exit(0);
        }
#endif

	ti_t = MPI_Wtime() ;	/* start timer */	
	t0_t = ti_t ;

    /* report on switches */
    fprintf(stderr,"Wind source: ") ;
    if(WIND) fprintf(stderr,"ENABLED\n") ;
    else fprintf(stderr,"DISABLED\n") ;

    fprintf(stderr,"Parabolic interpolation: ") ;
    if(RECON_PARA) fprintf(stderr,"ENABLED\n") ;
    else fprintf(stderr,"DISABLED\n") ;

    fprintf(stderr,"WENO-5 interpolation: ") ;
    if(RECON_WENO) fprintf(stderr,"ENABLED\n") ;
    else fprintf(stderr,"DISABLED\n") ;

	if((RECON_PARA || RECON_WENO) && NG == 2) {
		fprintf(stderr,"not enough ghost zones!\n") ;
		fprintf(stderr,"WENO or PARA + NG=2\n") ;
		exit(1);
	}

    fprintf(stderr,"Cooling: ") ;
    if(ALLOW_COOL) fprintf(stderr,"ENABLED\n") ;
    else fprintf(stderr,"DISABLED\n") ;

    fprintf(stderr,"Lagrangian tracers: ") ;
    if(NPTOT <= 0) fprintf(stderr,"DISABLED\n") ;
    else fprintf(stderr,"ENABLED\n") ;

	nstep = 0;

	/* Perform Initializations, either directly or via checkpoint */
	system("mkdir -p dumps images");
	if (!restart_init()) {
		init_ranc(1) ;
		init();
		if(NPTOT > 0) init_particles();		
	}

	/* do initial diagnostics */
	diag(INIT_OUT);

	tdump = t + DTd;
	if(NPTOT > 0) tpdump = t + DTp;
	timage = t + DTi;
	tlog = t + DTl;

	defcon = 1.;
	fprintf(stderr,"t, tf: %g %g\n",t,tf) ;
	dump();
	while (t < tf) {


		/* step variables forward in time */
		step_ch();

		/* uncomment to get step-by-step diagnostics */
		//diag(DUMP_OUT);
		//diag(IMAGE_OUT);

		//fprintf(stderr, "%10.5g %10.5g %8d\n", t, dt, nstep) ;

		/* Handle output frequencies: */
		if (t >= tdump) {
			//diag(DUMP_OUT);
			dump();
			tdump += DTd;
		}
		if (t >= timage) {
			//diag(IMAGE_OUT);
			timage += DTi;
		}
		if (t >= tlog) {
			diag(LOG_OUT);
			tlog += DTl;
		}
		if (NPTOT > 0 && t >= tpdump) {
			diag(PDUMP_OUT);
			tpdump += DTp;
		}

		/* restart dump */
		nstep++;
		if (nstep % DTr == 0)
			restart_write();

		/* regular timing */
		//if(nstep == 128) exit(0) ;
		if(nstep%TIMER_NSTEP == 0) {
			t1_t = MPI_Wtime() ;
			zcps = TIMER_NSTEP*N1*N2*N3/(t1_t-t0_t) ;
			fprintf(stderr,"Zone cycles/second: %g\n", zcps) ;
			t0_t = t1_t ;
		}

		/* deal with failed timestep, though we usually exit upon failure */
		if (failed) {
			//restart_init();
			failed = 0;
			nfailed = nstep;
			defcon = 0.3;
		}
		if (nstep > nfailed + DTr * 4. * (1 + 1. / defcon))
			defcon = 1.;


	}
	fprintf(stderr, "ns,ts: %d %d\n", nstep, nstep * N1 * N2 * N3);
	tf_t = MPI_Wtime() ;
	zcps = nstep*N1*N2*N3/(tf_t-ti_t) ;
	fprintf(stderr,"Zone cycles/second: %g\n", zcps) ;

	/* do final diagnostics */
	diag(FINAL_OUT);

	return (0);
}


/*****************************************************************/
/*****************************************************************
  zero_arrays():
  ----------

       -- sets to zero all arrays

 *****************************************************************/
void zero_arrays()
{
	int i, j, k;

	/* everything must be initialized to zero */
	ZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG, -NG, N3-1 + NG) {
		PLOOP {
			p[i][j][k][ip] = 0.;
			ph[i][j][k][ip] = 0.;
			dq[i][j][k][ip] = 0.;
			F1[i][j][k][ip] = 0.;
			F2[i][j][k][ip] = 0.;
			F3[i][j][k][ip] = 0.;
		}
		pflag[i][j][k] = 0;
	}

	//k = 0;
	//IMAGELOOP {
		//failimage[0][k] = failimage[1][k] = failimage[2][k] =
		    //failimage[3][k] = failimage[4][k++] = 0;
	//}
}


/*****************************************************************/
/*****************************************************************
  set_grid():
  ----------

       -- calculates all grid functions that remain constant 
          over time, such as the metric (gcov), inverse metric 
          (gcon), connection coefficients (conn), and sqrt of 
          the metric's determinant (gdet).

 *****************************************************************/
void set_grid()
{
	int i, j;
	double X[NDIM];

	/* set up boundaries, steps in coordinate grid */
	set_points();
	dV = dx[1] * dx[2] * dx[3];

	DLOOPA X[j] = 0.;

	GZSLOOP(-NG, N1-1 + NG, -NG, N2-1 + NG) {

		/* zone-centered */
		coord(i, j, CENT, X);
		gcov_func(X, ggeom[i][j][CENT].gcov);
		ggeom[i][j][CENT].g = gdet_func(ggeom[i][j][CENT].gcov);
		gcon_func(ggeom[i][j][CENT].gcov, ggeom[i][j][CENT].gcon);
		ggeom[i][j][CENT].alpha = 1.0 / sqrt(-ggeom[i][j][CENT].gcon[0][0]);

		/* only required in zone center... */
		conn_func(X, &ggeom[i][j][CENT], conn[i][j]);

		/* corner-centered */
		coord(i, j, CORN, X);
		gcov_func(X, ggeom[i][j][CORN].gcov);
		ggeom[i][j][CORN].g = gdet_func(ggeom[i][j][CORN].gcov);
		gcon_func(ggeom[i][j][CORN].gcov, ggeom[i][j][CORN].gcon);
		ggeom[i][j][CORN].alpha = 1.0 / sqrt(-ggeom[i][j][CORN].gcon[0][0]);

		/* r-face-centered */
		coord(i, j, FACE1, X);
		gcov_func(X, ggeom[i][j][FACE1].gcov);
		ggeom[i][j][FACE1].g = gdet_func(ggeom[i][j][FACE1].gcov);
		gcon_func(ggeom[i][j][FACE1].gcov, ggeom[i][j][FACE1].gcon);
		ggeom[i][j][FACE1].alpha = 1.0 / sqrt(-ggeom[i][j][FACE1].gcon[0][0]);

		/* theta-face-centered */
		coord(i, j, FACE2, X);
		gcov_func(X, ggeom[i][j][FACE2].gcov);
		ggeom[i][j][FACE2].g = gdet_func(ggeom[i][j][FACE2].gcov);
		gcon_func(ggeom[i][j][FACE2].gcov, ggeom[i][j][FACE2].gcon);
		ggeom[i][j][FACE2].alpha = 1.0 / sqrt(-ggeom[i][j][FACE2].gcon[0][0]);
	}

	/* done! */
}
