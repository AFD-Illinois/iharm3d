
/* M1: 
	added floor for erad
	cfg 12.27.13
*/

#include "decs.h"

/* apply floors to density, internal energy */

void fixup(grid_prim_type pv)
{
	int i, j, k;

#pragma omp parallel for \
 private(i,j,k) \
 collapse(3)
	ZLOOP fixup1zone(i, j, k, pv[i][j][k]);
}

void fixup1zone(int i, int j, int k, double pv[NPR])
{
	double r, th, X[NDIM], uuscal, rhoscal, rhoflr, uuflr;
	double f, gamma;
	struct of_geom *geom;

	coord(i, j, CENT, X);
	bl_coord(X, &r, &th);

	rhoscal = pow(r, -1.5);
	uuscal = rhoscal/r ;

	rhoflr = RHOMIN * rhoscal;
	uuflr = UUMIN * uuscal;

	if (rhoflr < RHOMINLIMIT) {
		//if(nstep > 0) fprintf(stdout,"applying rho floor1 in zone %d %d %d\n", i, j, k);
		rhoflr = RHOMINLIMIT;
	}
	if (uuflr < UUMINLIMIT) {
		//if(nstep > 0) fprintf(stdout,"applying u floor1 in zone %d %d %d\n", i, j, k);
		uuflr = UUMINLIMIT;
	}

	/* floor on density and internal energy density (momentum *not* conserved) */
	if (pv[RHO] < rhoflr || isnan(pv[RHO])) {
		pv[RHO] = rhoflr;
		//if(nstep > 0) fprintf(stdout,"applying rho floor2 in zone %d %d %d\n", i, j, k);
	}
	if (pv[UU] < uuflr || isnan(pv[UU])) {
		pv[UU] = uuflr;
		//if(nstep > 0) fprintf(stdout,"applying u floor2 in zone %d %d %d\n", i, j, k);
	}

	/* limit gamma wrt normal observer */
	geom = get_geometry(i, j, CENT) ;

	if (mhd_gamma_calc(pv, geom, &gamma)) {
		/* Treat gamma failure here as "fixable" for fixup_utoprim() */
		//if(nstep > 0) fprintf(stdout,"gamma failure in zone %d %d %d\n", i, j, k);
		pflag[i][j][k] = -333;
		//failimage[3][i + j * N1]++;
	} else {
		if (gamma > GAMMAMAX) {
			f = sqrt((GAMMAMAX * GAMMAMAX - 1.) /
				 (gamma * gamma - 1.)
			    );
			pv[U1] *= f;
			pv[U2] *= f;
			pv[U3] *= f;
		//	if(nstep > 0) fprintf(stdout,"fixing gamma in zone %d %d %d\n", i, j, k);
		}
	}

	return;
}



/*******************************************************************************************
  fixup_utoprim(): 

    -- figures out (w/ pflag[]) which stencil to use to interpolate bad point from neighbors;

 *******************************************************************************************/

#define FLOOP for(int ip=0;ip<B1;ip++)

void fixup_utoprim(grid_prim_type pv)
{
	int i, j, k;
	double sum[B1],wsum;

	/* Flip the logic of the pflag[] so that it now indicates which cells are good  */
	ZSLOOP(-NG, (N1-1 + NG), -NG, (N2-1 + NG), -NG, (N3-1 + NG)) {
		pflag[i][j][k] = !pflag[i][j][k];
	}

	// make sure we are not using ill defined corner regions
	for(i=0;i<NG;i++) {
	for(j=0;j<NG;j++) {
	for(k=0;k<NG;k++) {
		pflag[i][j][k] = 0;
		pflag[i+N1+NG][j][k] = 0;
		pflag[i][j+N2+NG][k] = 0;
		pflag[i][j][k+N3+NG] = 0;
		pflag[i+N1+NG][j+N2+NG][k] = 0;
		pflag[i+N1+NG][j][k+N3+NG] = 0;
		pflag[i][j+N2+NG][k+N3+NG] = 0;
		pflag[i+N1+NG][j+N2+NG][k+N3-1+NG] = 0;
	}}}

	/* Fix the interior points first */
//#pragma omp parallel for \
// private(i,j,k,sum,wsum) \
// schedule(dynamic) \
// collapse(3)
	ZSLOOP(0, (N1 - 1), 0, (N2 - 1), 0, (N3 - 1)) {
		if (pflag[i][j][k] == 0) {
			//if(nstep > 0) fprintf(stdout,"fixing utoprim %d %d %d\n", i, j, k);
			wsum = 0.;
			FLOOP sum[ip] = 0.;
			for(int l=-1;l<2;l++) {
			for(int m=-1;m<2;m++) {
			for(int n=-1;n<2;n++) {
				double w = 1./(abs(l)+abs(m)+abs(n) + 1) * pflag[i+l][j+m][k+n];
				wsum += w;
				FLOOP sum[ip] += w*pv[i+l][j+m][k+n][ip];
			}}}
			if(wsum < 1.e-10) {
				fprintf(stderr,"fixup_utoprim problem: wsum == 0!  no useable neighbors!\n");
				exit(7);
			}
			FLOOP pv[i][j][k][ip] = sum[ip]/wsum;

			pflag[i][j][k] = 1;	/* The cell has been fixed so we can use it for interpolation elsewhere */
			fixup1zone(i, j, k, pv[i][j][k]);	/* Floor and limit gamma the interpolated value */
		}
	}

	return;
}

#undef FLOOP
