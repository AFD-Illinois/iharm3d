
/* M1: no change */

#include "decs.h"

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

#define LOGR	0
#define SINHR	1
#define SINTH	0
#define POLYTH	1

//#define RCOORD	LOGR
#define RCOORD	SINHR

//#define THCOORD	SINTH
#define THCOORD	POLYTH

#if(RCOORD==SINHR)
static double sinhr_r0,sinhr_rtrans;
#endif
#if(THCOORD==POLYTH)
static double poly_alpha,poly_xt,poly_norm;
#endif

double r_of_x(double x)
{
	#if(RCOORD==LOGR)
	return exp(x) + R0;
	#elif(RCOORD==SINHR)
	return sinhr_r0 + sinhr_rtrans*sinh(x/sinhr_rtrans);
	#endif
}

double dr_dx(double x)
{
	#if(RCOORD==LOGR)
	return r_of_x(x) - R0;
	#elif(RCOORD==SINHR)
	return cosh(x/sinhr_rtrans);
	#endif
}

double th_of_x(double x)
{
	#if(THCOORD==SINTH)
	return M_PI * x + ((1. - hslope) / 2.) * sin(2. * M_PI * x);
	#elif(THCOORD==POLYTH)
	double y = 2*x - 1.;
	return poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*M_PI;
	#endif
}

double dth_dx(double x)
{
	#if(THCOORD==SINTH)
	return M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * x);
	#elif(THCOORD==POLYTH)
	double y = 2*x - 1.;
	return 2.*poly_norm*(1. + pow(y/poly_xt,poly_alpha));
	#endif
}

/* should return boyer-lindquist coordinate of point */
void bl_coord(double *X, double *r, double *th)
{
	*r = r_of_x(X[1]);
	*th = th_of_x(X[2]);

	// avoid singularity at polar axis
#if(COORDSINGFIX)
	if (fabs(*th) < SINGSMALL) {
		if ((*th) >= 0)
			*th = SINGSMALL;
		if ((*th) < 0)
			*th = -SINGSMALL;
	}
	if (fabs(M_PI - (*th)) < SINGSMALL) {
		if ((*th) >= M_PI)
			*th = M_PI + SINGSMALL;
		if ((*th) < M_PI)
			*th = M_PI - SINGSMALL;
	}
#endif

	return;
}

/* insert metric here */
void gcov_func(double *X, double gcov[][NDIM])
{
	int j, k;
	double sth, cth, s2, rho2;
	double r, th;
	double tfac, rfac, hfac, pfac;

	DLOOP gcov[j][k] = 0.;

	bl_coord(X, &r, &th);

	cth = cos(th);

	//-orig        sth = fabs(sin(th)) ;
	//-orig        if(sth < SMALL) sth = SMALL ;

	sth = sin(th);

	s2 = sth * sth;
	rho2 = r * r + a * a * cth * cth;

	tfac = 1.;
	rfac = dr_dx(X[1]);
	hfac = dth_dx(X[2]);
	pfac = 1.;

	gcov[TT][TT] = (-1. + 2. * r / rho2) * tfac * tfac;
	gcov[TT][1] = (2. * r / rho2) * tfac * rfac;
	gcov[TT][3] = (-2. * a * r * s2 / rho2) * tfac * pfac;

	gcov[1][TT] = gcov[TT][1];
	gcov[1][1] = (1. + 2. * r / rho2) * rfac * rfac;
	gcov[1][3] = (-a * s2 * (1. + 2. * r / rho2)) * rfac * pfac;

	gcov[2][2] = rho2 * hfac * hfac;

	gcov[3][TT] = gcov[TT][3];
	gcov[3][1] = gcov[1][3];
	gcov[3][3] =
	    s2 * (rho2 + a * a * s2 * (1. + 2. * r / rho2)) * pfac * pfac;
}

/* some grid location, dxs */
void set_points()
{
	#if(RCOORD==LOGR)
	startx[1] = log(Rin - R0);
	#elif(RCOORD==SINHR)
	startx[1] = 0.;
	#endif
	startx[2] = 0.;
	startx[3] = 0.;
	#if(RCOORD==LOGR)
	dx[1] = log((Rout - R0) / (Rin - R0)) / N1TOT;
	#elif(RCOORD==SINHR)
	dx[1] = 0.075;
	Rhor = 1. + sqrt(1. - a*a);
	sinhr_r0 = Rhor - 2.1*dx[1];
	sinhr_rtrans = 1.9;
	#endif
	dx[2] = 1. / N2TOT;
	dx[3] = 2.*M_PI / N3TOT;

	#if(THCOORD==POLYTH)
	poly_alpha = 6.;
	poly_xt = 0.65;
	poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha+1.) * 1./pow(poly_xt,poly_alpha));
	#endif
}

void fix_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3)
{
	int i, j, k;

	if(global_start[1] == 0) {
	ISLOOP(-1,N1)
	KSLOOP(-1,N3) {
		F1[i][-1+START2][k][B2] = -F1[i][0+START2][k][B2];
		F3[i][-1+START2][k][B2] = -F3[i][0+START2][k][B2];
		PLOOP F2[i][0+START2][k][ip] = 0.;
	}
	}
	if(global_stop[1] == N2TOT) {
	ISLOOP(-1,N1)
	KSLOOP(-1,N3) {
		F1[i][N2+START2][k][B2] = -F1[i][N2 - 1+START2][k][B2];
		F3[i][N2+START2][k][B2] = -F3[i][N2 - 1+START2][k][B2];
		PLOOP F2[i][N2+START2][k][ip] = 0.;
	}
	}

	if(global_start[0] == 0) {
	JSLOOP(-1,N2-1) 
	KSLOOP(-1,N3-1) {
		if (F1[0+START1][j][k][RHO] > 0.) F1[0+START1][j][k][RHO] = 0.;
	}
	}
	if(global_stop[0] == N1TOT) {
	JSLOOP(-1,N2-1)
	KSLOOP(-1,N3-1) {
		if (F1[N1+START1][j][k][RHO] < 0.) F1[N1+START1][j][k][RHO] = 0.;
	}
	}

	return;
}

