
/* M1: no changes yet */

#include "decs.h"

void dump(FILE * fp)
{
	int i, j, k;
	double divb;
	double X[NDIM];
	double r, th, vmin, vmax;
	struct of_geom *geom;
	struct of_state q;

	/***************************************************************
	  Write header information : 
	***************************************************************/

	fprintf(fp, FMT_DBL_OUT, t);
	fprintf(fp, FMT_INT_OUT, N1);
	fprintf(fp, FMT_INT_OUT, N2);
	fprintf(fp, FMT_DBL_OUT, startx[1]);
	fprintf(fp, FMT_DBL_OUT, startx[2]);
	fprintf(fp, FMT_DBL_OUT, dx[1]);
	fprintf(fp, FMT_DBL_OUT, dx[2]);
	fprintf(fp, FMT_DBL_OUT, tf);
	fprintf(fp, FMT_INT_OUT, nstep);
	fprintf(fp, FMT_DBL_OUT, a);
	fprintf(fp, FMT_DBL_OUT, gam);
	fprintf(fp, FMT_DBL_OUT, cour);
	fprintf(fp, FMT_DBL_OUT, DTd);
	fprintf(fp, FMT_DBL_OUT, DTl);
	fprintf(fp, FMT_DBL_OUT, DTi);
	fprintf(fp, FMT_INT_OUT, DTr);
	fprintf(fp, FMT_INT_OUT, dump_cnt);
	fprintf(fp, FMT_INT_OUT, image_cnt);
	fprintf(fp, FMT_INT_OUT, rdump_cnt);
	fprintf(fp, FMT_DBL_OUT, dt);
	fprintf(fp, FMT_INT_OUT, lim);
	fprintf(fp, FMT_INT_OUT, failed);
	fprintf(fp, FMT_DBL_OUT, Rin);
	fprintf(fp, FMT_DBL_OUT, Rout);
	fprintf(fp, FMT_DBL_OUT, hslope);
	fprintf(fp, FMT_DBL_OUT, R0);

	fprintf(fp, "\n");

	/* calculate charge density */
	current_calc() ;

	/***************************************************************
	  Write grid data:
	***************************************************************/

	ZSLOOP(0, N1 - 1, 0, N2 - 1) {
		coord(i, j, CENT, X);
		bl_coord(X, &r, &th);

		fprintf(fp, FMT_DBL_OUT, X[1]);
		fprintf(fp, FMT_DBL_OUT, X[2]);
		fprintf(fp, FMT_DBL_OUT, r);
		fprintf(fp, FMT_DBL_OUT, th);
		PLOOP fprintf(fp, FMT_DBL_OUT, p[i][j][k]);

		/* divb flux-ct defn; corner-centered.  Use
		   only interior corners */
		if (i > 0+NG && j > 0+NG && i < N1+NG && j < N2+NG) {
			divb = fabs(0.5*(p[i][j][B1]*ggeom[i][j][CENT].g
				+ p[i][j-1][B1]*ggeom[i][j-1][CENT].g
				- p[i-1][j][B1]*ggeom[i-1][j][CENT].g
				- p[i-1][j-1][B1]*ggeom[i-1][j-1][CENT].g
				)/dx[1] +
				0.5*(p[i][j][B2]*ggeom[i][j][CENT].g
				+ p[i-1][j][B2]*ggeom[i-1][j][CENT].g
				- p[i][j-1][B2]*ggeom[i][j-1][CENT].g
				- p[i-1][j-1][B2]*ggeom[i-1][j-1][CENT].g
				)/dx[2]);
		} else
			divb = 0.;

		fprintf(fp, FMT_DBL_OUT, divb);


		if (!failed) {
			geom = get_geometry(i, j, CENT) ;
			get_state(p[i][j], geom, &q);

			for (k = 0; k < NDIM; k++)
				fprintf(fp, FMT_DBL_OUT, q.ucon[k]);
			for (k = 0; k < NDIM; k++)
				fprintf(fp, FMT_DBL_OUT, q.ucov[k]);
			for (k = 0; k < NDIM; k++)
				fprintf(fp, FMT_DBL_OUT, q.bcon[k]);
			for (k = 0; k < NDIM; k++)
				fprintf(fp, FMT_DBL_OUT, q.bcov[k]);

			mhd_vchar(p[i][j], &q, geom, 1, &vmax, &vmin);
			fprintf(fp, FMT_DBL_OUT, vmin);
			fprintf(fp, FMT_DBL_OUT, vmax);

			mhd_vchar(p[i][j], &q, geom, 2, &vmax, &vmin);
			fprintf(fp, FMT_DBL_OUT, vmin);
			fprintf(fp, FMT_DBL_OUT, vmax);

			fprintf(fp, FMT_DBL_OUT, geom->g);

			/* current density */
			for (k = 0; k < NDIM; k++)
			        fprintf(fp, FMT_DBL_OUT, Jcon[i][j][k]) ;
                        double jcov[NDIM] ;
                        lower(Jcon[i][j], geom, jcov) ;
			for (k = 0; k < NDIM; k++)
			        fprintf(fp, FMT_DBL_OUT, jcov[k]) ;
                        
		}

		fprintf(fp, "\n");
	}
}
