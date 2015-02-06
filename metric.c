
/*

	M1: changed to eliminate use of old scn routine lu.c,
	use gsl routines instead.
	eliminates responsibility for some code
	cfg 12.24.14

*/

#include "decs.h"

/***************************************************************************/
/***************************************************************************
    coord():
    -------
       -- given the indices i,j and location in the cell, return with 
          the values of X1,X2 there;  
       -- the locations are defined by : 
           -----------------------
           |                     |
           |                     |
           |FACE1   CENT         |
           |                     |
           |CORN    FACE2        |
           ----------------------
***************************************************************************/
void coord(int i, int j, int loc, double *X)
{
	i += global_start[0];
	j += global_start[1];
	if (loc == FACE1) {
		X[1] = startx[1] + (i - START1) * dx[1];
		X[2] = startx[2] + (j + 0.5 - START2) * dx[2];
	} else if (loc == FACE2) {
		X[1] = startx[1] + (i + 0.5 - START1) * dx[1];
		X[2] = startx[2] + (j - START2) * dx[2];
	} else if (loc == CENT) {
		X[1] = startx[1] + (i + 0.5 - START1) * dx[1];
		X[2] = startx[2] + (j + 0.5 - START2) * dx[2];
	} else {
		X[1] = startx[1] + (i - START1) * dx[1];
		X[2] = startx[2] + (j - START2) * dx[2];
	}

	return;
}

/* assumes gcov has been set first; returns determinant */
/* 
   switched to use gsl
   cfg 12.24.14 
*/
double gdet_func(double gcov[][NDIM])
{
	int i,j,signum ;
	double detg ;
	gsl_matrix *gsl_gcov ;
	gsl_permutation *perm ;

	gsl_gcov = gsl_matrix_alloc(NDIM,NDIM) ;
	perm = gsl_permutation_alloc(NDIM) ;

	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++) gsl_matrix_set(gsl_gcov,i,j,gcov[i][j]) ;

	if(gsl_linalg_LU_decomp(gsl_gcov, perm, &signum) != 0) {
		fprintf(stderr,
			"gdet_func(): singular matrix encountered! \n");
		fail(FAIL_METRIC);
	}

	detg = gsl_linalg_LU_det(gsl_gcov, signum) ;

	/* clean up after yourself */
	gsl_matrix_free(gsl_gcov) ;
	gsl_permutation_free(perm) ;

	return (sqrt(fabs(detg)));
}

/* invert gcov to get gcon */
/*

	switched to use gsl
	cfg 12.24.14

*/
void gcon_func(double gcov[][NDIM], double gcon[][NDIM])
{
	int i,j,signum ;
	gsl_matrix *gsl_gcov,*gsl_gcon ;
	gsl_permutation *perm ;

	gsl_gcov = gsl_matrix_alloc(NDIM,NDIM) ;
	gsl_gcon = gsl_matrix_alloc(NDIM,NDIM) ;
	perm = gsl_permutation_alloc(NDIM) ;

	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++) gsl_matrix_set(gsl_gcov,i,j,gcov[i][j]) ;

	if(gsl_linalg_LU_decomp(gsl_gcov, perm, &signum) != 0) {
		fprintf(stderr,
			"gcon_func(): singular matrix encountered! \n");
		fail(FAIL_METRIC);
	}

	if(gsl_linalg_LU_invert(gsl_gcov, perm, gsl_gcon) != 0) {
		fprintf(stderr,
			"gcon_func(): problem inverting matrix! \n");
		fail(FAIL_METRIC);
	}

	for(i=0;i<NDIM;i++)
	for(j=0;j<NDIM;j++) gcon[i][j] = gsl_matrix_get(gsl_gcon,i,j) ;

	/* clean up after yourself */
	gsl_matrix_free(gsl_gcov) ;
	gsl_matrix_free(gsl_gcon) ;
	gsl_permutation_free(perm) ;

	/* done! */
}

/***************************************************************************/
/***************************************************************************
  conn_func():
  -----------

   -- this gives the connection coefficient
	\Gamma^{i}_{j,k} = conn[..][i][j][k]
   --  where i = {1,2,3,4} corresponds to {t,r,theta,phi}

***************************************************************************/

/* Sets the spatial discretization in numerical derivatives : */
#define DELTA 1.e-5

/* NOTE: parameter hides global variable */
void conn_func(double *X, struct of_geom *geom, double conn[][NDIM][NDIM])
{
	int i, j, k, l;
	double tmp[NDIM][NDIM][NDIM];
	double Xh[NDIM], Xl[NDIM];
	double gh[NDIM][NDIM];
	double gl[NDIM][NDIM];

	for (k = 0; k < NDIM; k++) {
		for (l = 0; l < NDIM; l++)
			Xh[l] = X[l];
		for (l = 0; l < NDIM; l++)
			Xl[l] = X[l];
		Xh[k] += DELTA;
		Xl[k] -= DELTA;
		gcov_func(Xh, gh);
		gcov_func(Xl, gl);

		for (i = 0; i < NDIM; i++)
			for (j = 0; j < NDIM; j++)
				conn[i][j][k] =
				    (gh[i][j] - gl[i][j]) / (Xh[k] -
							     Xl[k]);
	}

	/* now rearrange to find \Gamma_{ijk} */
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
			for (k = 0; k < NDIM; k++)
				tmp[i][j][k] =
				    0.5 * (conn[j][i][k] + conn[k][i][j] -
					   conn[k][j][i]);

	/* finally, raise index */
	for (i = 0; i < NDIM; i++)
		for (j = 0; j < NDIM; j++)
			for (k = 0; k < NDIM; k++) {
				conn[i][j][k] = 0.;
				for (l = 0; l < NDIM; l++)
					conn[i][j][k] +=
					    geom->gcon[i][l] *
					    tmp[l][j][k];
			}

	/* done! */
}

/* Lowers a contravariant rank-1 tensor to a covariant one */
void lower(double *ucon, struct of_geom *geom, double *ucov)
{

	ucov[0] = geom->gcov[0][0] * ucon[0]
	    + geom->gcov[0][1] * ucon[1]
	    + geom->gcov[0][2] * ucon[2]
	    + geom->gcov[0][3] * ucon[3];
	ucov[1] = geom->gcov[1][0] * ucon[0]
	    + geom->gcov[1][1] * ucon[1]
	    + geom->gcov[1][2] * ucon[2]
	    + geom->gcov[1][3] * ucon[3];
	ucov[2] = geom->gcov[2][0] * ucon[0]
	    + geom->gcov[2][1] * ucon[1]
	    + geom->gcov[2][2] * ucon[2]
	    + geom->gcov[2][3] * ucon[3];
	ucov[3] = geom->gcov[3][0] * ucon[0]
	    + geom->gcov[3][1] * ucon[1]
	    + geom->gcov[3][2] * ucon[2]
	    + geom->gcov[3][3] * ucon[3];

	return;
}

/* Raises a covariant rank-1 tensor to a contravariant one */
void raise(double *ucov, struct of_geom *geom, double *ucon)
{

	ucon[0] = geom->gcon[0][0] * ucov[0]
	    + geom->gcon[0][1] * ucov[1]
	    + geom->gcon[0][2] * ucov[2]
	    + geom->gcon[0][3] * ucov[3];
	ucon[1] = geom->gcon[1][0] * ucov[0]
	    + geom->gcon[1][1] * ucov[1]
	    + geom->gcon[1][2] * ucov[2]
	    + geom->gcon[1][3] * ucov[3];
	ucon[2] = geom->gcon[2][0] * ucov[0]
	    + geom->gcon[2][1] * ucov[1]
	    + geom->gcon[2][2] * ucov[2]
	    + geom->gcon[2][3] * ucov[3];
	ucon[3] = geom->gcon[3][0] * ucov[0]
	    + geom->gcon[3][1] * ucov[1]
	    + geom->gcon[3][2] * ucov[2]
	    + geom->gcon[3][3] * ucov[3];

	return;
}

/* load local geometry into structure geom */
struct of_geom *get_geometry(int ii, int jj, int kk)
{
	icurr = ii;
	jcurr = jj;
	pcurr = kk;

	return( &(ggeom[ii][jj][kk]) ) ;
}

#undef DELTA

/* Minkowski metric; signature +2 */
double mink(int i, int j)
{
	if (i == j) {
		if (i == 0)
			return (-1.);
		else
			return (1.);
	} else
		return (0.);
}

/* Boyer-Lindquist ("bl") metric functions */
void blgset(int i, int j, struct of_geom *geom)
{
	double r, th, X[NDIM];

	coord(i, j, CENT, X);
	bl_coord(X, &r, &th);

	if (th < 0)
		th *= -1.;
	if (th > M_PI)
		th = 2. * M_PI - th;

	geom->g = bl_gdet_func(r, th);
	bl_gcov_func(r, th, geom->gcov);
	bl_gcon_func(r, th, geom->gcon);
}

double bl_gdet_func(double r, double th)
{
	double a2, r2;

	a2 = a * a;
	r2 = r * r;
	return (r * r * fabs(sin(th)) *
		(1. + 0.5 * (a2 / r2) * (1. + cos(2. * th)))
	    );
}

void bl_gcov_func(double r, double th, double gcov[][NDIM])
{
	int j, k;
	double sth, cth, s2, a2, r2, DD, mu;

	DLOOP gcov[j][k] = 0.;

	sth = fabs(sin(th));
	s2 = sth * sth;
	cth = cos(th);
	a2 = a * a;
	r2 = r * r;
	DD = 1. - 2. / r + a2 / r2;
	mu = 1. + a2 * cth * cth / r2;

	gcov[TT][TT] = -(1. - 2. / (r * mu));
	gcov[TT][3] = -2. * a * s2 / (r * mu);
	gcov[3][TT] = gcov[TT][3];
	gcov[1][1] = mu / DD;
	gcov[2][2] = r2 * mu;
	gcov[3][3] =
	    r2 * sth * sth * (1. + a2 / r2 + 2. * a2 * s2 / (r2 * r * mu));

}

void bl_gcon_func(double r, double th, double gcon[][NDIM])
{
	int j, k;
	double sth, cth, a2, r2, r3, DD, mu;

	DLOOP gcon[j][k] = 0.;

	sth = sin(th);
	cth = cos(th);

#if(COORDSINGFIX)
	if (fabs(sth) < SINGSMALL) {
		if (sth >= 0)
			sth = SINGSMALL;
		if (sth < 0)
			sth = -SINGSMALL;
	}
#endif

	a2 = a * a;
	r2 = r * r;
	r3 = r2 * r;
	DD = 1. - 2. / r + a2 / r2;
	mu = 1. + a2 * cth * cth / r2;

	gcon[TT][TT] = -1. - 2. * (1. + a2 / r2) / (r * DD * mu);
	gcon[TT][3] = -2. * a / (r3 * DD * mu);
	gcon[3][TT] = gcon[TT][3];
	gcon[1][1] = DD / mu;
	gcon[2][2] = 1. / (r2 * mu);
	gcon[3][3] = (1. - 2. / (r * mu)) / (r2 * sth * sth * DD);

}
