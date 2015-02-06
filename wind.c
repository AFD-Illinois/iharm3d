
#include "decs.h"

void wind_source(double *prim, struct of_geom *geom, struct of_state *q,
			int ii, int jj, double *dU)
{
        double X[NDIM],r,th,mu;
	double NORM,sig_pair,rISCO,dndt,ir,rst;
	double hfac,rfac,Gcon[NDIM],GconKS[NDIM],Gcov[NDIM];

        /* need coordinates to evaluate particle addtn rate */
        coord(ii, jj, CENT, X);
        bl_coord(X, &r, &th);
	mu = cos(th);

	/** based on Mosc & Gammie, eq 26,37 **/
	NORM = 1.e-0 ;
	sig_pair = 0.3 ;	/* pair production scale heigh */
	rISCO = 6.*(1. - 3.*a); /* approximate ISCO radius */
	/* invariant production rate density */
	dndt = NORM*pow(r,-6)*exp(-0.5*mu*mu/(sig_pair*sig_pair)) ;
	ir = 1./r ;

	GconKS[0] = 300.*dndt*ir ;
	rst = rISCO*(1.+0.5*mu*mu);
	GconKS[1] = 20.*(r - rst)*dndt*ir*ir ;
	GconKS[2] = mu*dndt*ir*ir;
	GconKS[3] = 150.*dndt*ir*ir;

	/* convert to MKS */
	rfac = r - R0;
	hfac = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);

	Gcon[0] = GconKS[0];
	Gcon[1] = GconKS[1]/rfac;
	Gcon[2] = GconKS[2]/hfac;
	Gcon[3] = GconKS[3];

	lower(Gcon,geom,Gcov);

	dU[0] += dndt*q->ucon[0] ;
	dU[1] += Gcov[0];
	dU[2] += Gcov[1];
	dU[3] += Gcov[2];
	dU[4] += Gcov[3];

}


