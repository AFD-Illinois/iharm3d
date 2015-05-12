
/* M1: no changes yet */

#include "decs.h"

/* all diagnostics subroutine */

void diag(int call_code)
{
	char dfnam[100] ;
	int i, j, k;
	//FILE *dump_file;
	FILE *pdump_file;
	double U[NPR], divb, e_fin, m_fin ;
	double pp = 0. ;
	double divbmax = 0.;
	double rmed = 0.;
	double e = 0.;
	struct of_geom *geom;
	struct of_state q;
	//int imax = 0;
	//int jmax = 0;
	//int kmax = 0;

	static double e_init, m_init;
	static FILE *ener_file;

	if (call_code == INIT_OUT) {
		/* set things up */
		if(mpi_io_proc()) {
			ener_file = fopen("ener.out", "a");
			if (ener_file == NULL) {
				fprintf(stderr,
					"error opening energy output file\n");
				exit(1);
			}
		}
	}

	/* calculate conserved quantities */
	if ((call_code == INIT_OUT || call_code == LOG_OUT || call_code == FINAL_OUT) && !failed) {
		pp = 0.;
		e = 0.;
		rmed = 0.;
		divbmax = 0.;
		//imax = 0;
		//jmax = 0;
		//kmax = 0;
		ZSLOOP(0, N1 - 1, 0, N2 - 1, 0, N3 - 1) {
			geom = get_geometry(i, j, CENT) ;
			get_state(p[i][j][k], geom, &q);
			primtoU(p[i][j][k], &q, geom, U);

			rmed += U[RHO] * dV;
			pp += U[U3] * dV;
			e += U[UU] * dV;

			divb = flux_ct_divb(i, j, k);

			if (divb > divbmax) {
				//imax = i;
				//jmax = j;
				//kmax = k;
				divbmax = divb;
			}
		}
	}

	rmed = mpi_io_reduce(rmed);
	pp = mpi_io_reduce(pp);
	e = mpi_io_reduce(e);
	divbmax = mpi_io_max(divbmax);

	if (call_code == INIT_OUT) {
		e_init = e;
		m_init = rmed;
	}

	if (call_code == FINAL_OUT) {
		if(mpi_io_proc()) {
			e_fin = e;
			m_fin = rmed;
			if(mpi_io_proc()) {
				fprintf(stderr, "\n\nEnergy: ini,fin,del: %g %g %g\n",
					e_init, e_fin, (e_fin - e_init) / e_init);
				fprintf(stderr, "mass: ini,fin,del: %g %g %g\n",
					m_init, m_fin, (m_fin - m_init) / m_init);
			}
		}
	}

	if (call_code == INIT_OUT ||
	    call_code == LOG_OUT || call_code == FINAL_OUT) {
		if(mpi_io_proc()) {
			fprintf(stderr, "LOG      t=%g \t divbmax: %g\n",
				t,divbmax);
			fprintf(ener_file, "%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ", 
				t,rmed,pp,e,
				p[N1/2][N2/2][N3/2][UU]*pow(p[N1/2][N2/2][N3/2][RHO], -gam),
				p[N1/2][N2/2][N3/2][UU]);
			fprintf(ener_file, "%15.8g %15.8g %15.8g ",mdot,edot,ldot);
			fprintf(ener_file, "\n");
			fflush(ener_file);
		}
	}


	/* dump at regular intervals */
	if (call_code == INIT_OUT || call_code == DUMP_OUT) {
		/* make regular dump file */
		//sprintf(dfnam, "dumps/dump%03d", dump_cnt);
		//fprintf(stderr, "DUMP     file=%s\n", dfnam);
		//dump_file = fopen(dfnam, "w");

		//if (dump_file == NULL) {
		//	fprintf(stderr, "error opening dump file\n");
		//	exit(2);
		//}

		//dump(dump_file);
		//fclose(dump_file);

		dump_cnt++;
	}

	if (NPTOT > 0 && (call_code == INIT_OUT || call_code == PDUMP_OUT)) {
		/* make lagrangian particle dump file */
		sprintf(dfnam, "dumps/pdump%03d", pdump_cnt);
		fprintf(stderr, "PDUMP     file=%s\n", dfnam);
		pdump_file = fopen(dfnam, "w");

		if (pdump_file == NULL) {
			fprintf(stderr, "error opening pdump file\n");
			exit(2);
		}

		pdump(pdump_file);
		fclose(pdump_file);

		pdump_cnt++;
	}

	/* image dump at regular intervals */
	if (call_code == IMAGE_OUT ||
	    call_code == INIT_OUT || call_code == FINAL_OUT) {

		//image_all(image_cnt);

		image_cnt++;
	}
}


/** some diagnostic routines **/
void fail(int fail_type)
{

	failed = 1;

	fprintf(stderr, "\n\nfail: %d %d %d\n", icurr, jcurr, fail_type);

	area_map(icurr, jcurr, 0, p);

	fprintf(stderr, "fail\n");

	diag(FINAL_OUT);

	/* for diagnostic purposes */
	exit(0);
}



/* map out region around failure point */
void area_map(int i, int j, int k, grid_prim_type prim)
{
	fprintf(stderr, "area map\n");

	PLOOP {
		fprintf(stderr, "variable %d \n", ip);
		fprintf(stderr, "i = \t %12d %12d %12d\n", i - 1, i,
			i + 1);
		fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j + 1,
			prim[i - 1][j + 1][k][ip], prim[i][j + 1][k][ip],
			prim[i + 1][j + 1][k][ip]);
		fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j,
			prim[i - 1][j][k][ip], prim[i][j][k][ip],
			prim[i + 1][j][k][ip]);
		fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j - 1,
			prim[i - 1][j - 1][k][ip], prim[i][j - 1][k][ip],
			prim[i + 1][j - 1][k][ip]);
	}

	/* print out other diagnostics here */

}

/* evaluate flux based diagnostics; put results in
 * global variables */
void diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3)
{
    int j,k;
	mdot = edot = ldot = 0.;
#pragma omp parallel for \
 private(j,k) \
 reduction(+:mdot) reduction(-:edot) reduction(+:ldot) \
 collapse(2)    
	for (j = 0; j < N2; j++) {
	for (k = 0; k < N3; k++) {
	
		mdot += F1[0+START1][j][k][RHO] * dx[2] * dx[3];
		edot -= (F1[0+START1][j][k][UU] - F1[0+START1][j][k][RHO]) * dx[2] * dx[3];
		ldot += F1[0+START1][j][k][U3] * dx[2] * dx[3];
	}
	}
}


double flux_ct_divb(int i, int j, int k)
{

	if(i > 0+NG && j > 0+NG && k > 0+NG && i < N1+NG && j < N2+NG && k < N3+NG) {
		return fabs(0.25*(
					p[i][j][k][B1]*ggeom[i][j][CENT].g 
					+ p[i][j-1][k][B1]*ggeom[i][j-1][CENT].g
					+ p[i][j][k-1][B1]*ggeom[i][j][CENT].g
					+ p[i][j-1][k-1][B1]*ggeom[i][j-1][CENT].g
					- p[i-1][j][k][B1]*ggeom[i-1][j][CENT].g
					- p[i-1][j-1][k][B1]*ggeom[i-1][j-1][CENT].g
					- p[i-1][j][k-1][B1]*ggeom[i-1][j][CENT].g
					- p[i-1][j-1][k-1][B1]*ggeom[i-1][j-1][CENT].g
					)/dx[1] +
					0.25*(
					p[i][j][k][B2]*ggeom[i][j][CENT].g
					+ p[i-1][j][k][B2]*ggeom[i-1][j][CENT].g
					+ p[i][j][k-1][B2]*ggeom[i][j][CENT].g
					+ p[i-1][j][k-1][B2]*ggeom[i-1][j][CENT].g
					- p[i][j-1][k][B2]*ggeom[i][j-1][CENT].g
					- p[i-1][j-1][k][B2]*ggeom[i-1][j-1][CENT].g
					- p[i][j-1][k-1][B2]*ggeom[i][j-1][CENT].g
					- p[i-1][j-1][k-1][B2]*ggeom[i-1][j-1][CENT].g
					)/dx[2] + 
					0.25*(
					p[i][j][k][B3]*ggeom[i][j][CENT].g
					+ p[i][j-1][k][B3]*ggeom[i][j-1][CENT].g
					+ p[i-1][j][k][B3]*ggeom[i-1][j][CENT].g
					+ p[i-1][j-1][k][B3]*ggeom[i-1][j-1][CENT].g
					- p[i][j][k-1][B3]*ggeom[i][j][CENT].g
					- p[i][j-1][k-1][B3]*ggeom[i][j-1][CENT].g
					- p[i-1][j][k-1][B3]*ggeom[i-1][j][CENT].g
					- p[i-1][j-1][k-1][B3]*ggeom[i-1][j-1][CENT].g
					)/dx[3]);
	} else {
		return 0.;
	}
}

