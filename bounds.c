
/* M1: no change */

#include "decs.h"

void inflow_check(double *pr, int ii, int jj, int type);

void bound_prim(grid_prim_type prim)
{
	int i, j, k, iz, jrefl;
    double rescale_fac;

	sync_mpi_boundaries(prim);

	if(global_start[0] == 0) {
	/* inner r boundary condition: copy last value off grid */
#pragma omp parallel for \
 private(i,j,k,iz,rescale_fac) \
 collapse(2)
	JSLOOP(0,N2-1) {
	KSLOOP(0,N3-1) {
		iz = 0+START1 ;
		ISLOOP(-NG,-1) {
			PLOOP prim[i][j][k][ip] = prim[iz][j][k][ip];
			pflag[i][j][k] = pflag[iz][j][k];

			rescale_fac = ggeom[i][j][CENT].g/ggeom[iz][j][CENT].g ;
			prim[i][j][k][B1] *= rescale_fac ;
			prim[i][j][k][B2] *= rescale_fac ;
			prim[i][j][k][B3] *= rescale_fac ;
		}
	}
	}
	}

	if(global_stop[0] == N1TOT) {
	/* outer r BC: copy last value off grid */
#pragma omp parallel for \
 private(i,j,k,iz,rescale_fac) \
 collapse(2)
	JSLOOP(0,N2-1) {
	KSLOOP(0,N3-1) {
		iz = N1-1+START1 ;
		ISLOOP(N1,N1-1+NG) {
			PLOOP prim[i][j][k][ip] = prim[iz][j][k][ip];
			pflag[i][j][k] = pflag[iz][j][k];

			rescale_fac = ggeom[i][j][CENT].g/ggeom[iz][j][CENT].g ;
			prim[i][j][k][B1] *= rescale_fac ;
			prim[i][j][k][B2] *= rescale_fac ;
			prim[i][j][k][B3] *= rescale_fac ;
			/*
			*/
		}
	}
	}
	}

	if(global_start[1] == 0) {
	/* polar BCs */
#pragma omp parallel for \
 private(i,j,k,jrefl) \
 collapse(2)
	ISLOOP(-NG,N1-1+NG) {
	KSLOOP(-NG,N3-1+NG) {
		JSLOOP(-NG,-1) {
			jrefl = -j+2*NG-1 ;
			PLOOP prim[i][j][k][ip] = prim[i][jrefl][k][ip];
			pflag[i][j][k] = pflag[i][jrefl][k];
			prim[i][j][k][U2] *= -1.;
			prim[i][j][k][B2] *= -1.;
		}
	}
	}
	}

	if(global_stop[1] == N2TOT) {
#pragma omp parallel for \
 private(i,j,k,jrefl) \
 collapse(2)
	ISLOOP(-NG,N1-1+NG) {
	KSLOOP(-NG,N3-1+NG) {
		JSLOOP(N2,N2-1+NG) {
			jrefl = -j+2*(N2+START2)-1 ;
			PLOOP prim[i][j][k][ip] = prim[i][jrefl][k][ip];
			pflag[i][j][k] = pflag[i][jrefl][k];
			prim[i][j][k][U2] *= -1.;
			prim[i][j][k][B2] *= -1.;

		}

	}
	}
	}



#if(N3==1)
	if(global_start[2] == 0 && global_stop[2] == N3TOT) {
	/* periodic boundary in phi */
#pragma omp parallel for \
 private(i,j,k) \
 collapse(2)
	ISLOOP(-NG,N1-1+NG) {
	JSLOOP(-NG,N2-1+NG) {
		KSLOOP(-NG,-1) {
			PLOOP prim[i][j][k][ip] = prim[i][j][START3][ip];
		}
		KSLOOP(N3,N3-1+NG) {
			PLOOP prim[i][j][k][ip] = prim[i][j][START3][ip];
		}
	}
	}
	}
#else 

	if(global_start[2] == 0 && global_stop[2] == N3TOT) {
	/* periodic boundary in phi */
#pragma omp parallel for \
 private(i,j,k) \
 collapse(2)
	ISLOOP(-NG,N1-1+NG) {
	JSLOOP(-NG,N2-1+NG) {
		KSLOOP(-NG,-1) {
			int kper = N3+k;
			PLOOP prim[i][j][k][ip] = prim[i][j][kper][ip];
		}
		KSLOOP(N3,N3-1+NG) {
			int kper = k-N3;
			PLOOP prim[i][j][k][ip] = prim[i][j][kper][ip];
		}
	}
	}
	}
#endif

	if(global_start[0] == 0) {
	/* make sure there is no inflow at the inner boundary */
	ISLOOP(-NG,-1) 
		JSLOOP(-NG,N2-1+NG) 
		KSLOOP(-NG,N3-1+NG)
			inflow_check(prim[i][j][k], i, j, 0);
	}

	if(global_stop[0] == N1TOT) {
	/* make sure there is no inflow at the outer boundary */
	ISLOOP(N1,N1-1+NG)
		JSLOOP(-NG,N2-1+NG)
		KSLOOP(-NG,N3-1+NG) 
			inflow_check(prim[i][j][k], i, j, 1);
	}

}

void inflow_check(double *pr, int ii, int jj, int type)
{
	struct of_geom *geom;
	double ucon[NDIM];
	int j, k;
	double alpha, beta1, gamma, vsq;

	geom = get_geometry(ii, jj, CENT);
	ucon_calc(pr, geom, ucon);

	if (((ucon[1] > 0.) && (type == 0))
	    || ((ucon[1] < 0.) && (type == 1))) {
		/* find gamma and remove it from primitives */
		if (mhd_gamma_calc(pr, geom, &gamma)) {
			fflush(stderr);
			fprintf(stderr,
				"\ninflow_check(): gamma failure \n");
			fflush(stderr);
			fail(FAIL_GAMMA);
		}
		pr[U1] /= gamma;
		pr[U2] /= gamma;
		pr[U3] /= gamma;
		alpha = geom->alpha ;
		beta1 = geom->gcon[0][1] * alpha * alpha;

		/* reset radial velocity so radial 4-velocity
		 * is zero */
		pr[U1] = beta1 / alpha;

		/* now find new gamma and put it back in */
		vsq = 0.;
		SLOOP vsq +=
		    geom->gcov[j][k] * pr[U1 + j - 1] * pr[U1 + k - 1];
		if (fabs(vsq) < 1.e-13)
			vsq = 1.e-13;
		if (vsq >= 1.) {
			vsq = 1. - 1. / (GAMMAMAX * GAMMAMAX);
		}
		gamma = 1. / sqrt(1. - vsq);
		pr[U1] *= gamma;
		pr[U2] *= gamma;
		pr[U3] *= gamma;

		/* done */
	} else
		return;

}
