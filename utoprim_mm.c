
#include "decs.h"

#define ERRTOL	1.e-10
#define ITERMAX 8
#define DEL	1.e-5

/*

utoprim_mm.c: 
---------------

    Uses Mignone \& McKinney (2007) approach to variable inversion for a
    	general equation of state.

    vastly simplified variable inversion:
    do not use derivatives of EOS functions
    instead: Halley-Muller-Bailey step (numerical 2nd derivative)
    then set of secant steps.

    CFG 2-25-10

    added radiation sector in M1 closure

    CFG 12-23-13

*/

int Utoprim_mm(double U[NPR], struct of_geom *geom, double prim[NPR])
{
	int i, iter, eflag ;
	double Bcon[NDIM],Bcov[NDIM],Qcov[NDIM],Qcon[NDIM],
		ncov[NDIM],ncon[NDIM],Qtcon[NDIM];
	double lapse,Bsq,D,Ep,QdB,Qtsq,Qsq,Qdotn,rho0,u,gamma,Wp,Wp1,W,w,P,
		err,err1,dW ;
	double Pressure_rho0_u(double rho0, double u), 
		Pressure_rho0_w(double rho0, double w),
		err_eqn(double Bsq, double D, double Ep, double QdB, double Qtsq, double Wp, int *eflag),
		gamma_func(double Bsq, double D, double QdB, double
		Qtsq, double Wp, int *eflag),
		Wp_func(double *prim, struct of_geom *geom, int *eflag) ;
	//static long int cnt = 0 ;
	//static long int itercnt[ITERMAX] = {0} ;

	eflag = 0 ;

	/* CATCH: negative density */
	if (U[RHO] <= 0.) {
		eflag = -100 ;
		return (eflag) ;
	}

	/* Update the primitive B-fields */
	prim[B1] = U[B1]/geom->g  ;
	prim[B2] = U[B2]/geom->g  ;
	prim[B3] = U[B3]/geom->g  ;

	lapse = geom->alpha ;

	/* convert from conserved variables to four-vectors */
        D = U[RHO] * lapse / geom->g ;
        Bcon[0] = 0.;
	Bcon[1] = U[B1] * lapse / geom->g ;
	Bcon[2] = U[B2] * lapse / geom->g ;
	Bcon[3] = U[B3] * lapse / geom->g ;
	Qcov[0] = (U[UU] - U[RHO]) * lapse / geom->g ;
	Qcov[1] = U[U1] * lapse / geom->g ;
	Qcov[2] = U[U2] * lapse / geom->g ;
	Qcov[3] = U[U3] * lapse / geom->g ;

        lower(Bcon, geom, Bcov);
        Bsq = 0.; for (i = 0; i < NDIM; i++) Bsq += Bcon[i] * Bcov[i];
        raise(Qcov, geom, Qcon);
        QdB = 0.; for (i = 0; i < NDIM; i++) QdB += Qcov[i] * Bcon[i];

	ncov[0] = -lapse ; 
	ncov[1] = 0. ;
	ncov[2] = 0. ;
	ncov[3] = 0. ;
        raise(ncov, geom, ncon);
        Qdotn = 0.; for(i = 0 ; i < NDIM; i++) Qdotn += Qcon[i]*ncov[i] ;
        Qsq = 0.;  for (i = 0; i < NDIM; i++) Qsq += Qcov[i] * Qcon[i];

	for (i = 0; i < 4; i++) Qtcon[i] = Qcon[i] + ncon[i]*Qdotn ;
        Qtsq = Qsq + Qdotn * Qdotn;

	/* now set up eqtn for W'; this is the energy density */
	Ep = -Qdotn - D ;

	/** root finding section **/
	/* first find guesses from primitives.  These should be
	   OK because primitives have been through the repair
	   process */
	Wp = Wp_func(prim, geom, &eflag) ;
	if(eflag) return(eflag) ;

	double dedW,dedW2,f,errp,errm,Wpm,Wpp,h ;

	/* step around the guess & evaluate errors */
	Wpm = (1. - DEL)*Wp ; //heuristic
	h = Wp - Wpm ;
	Wpp = Wp + h ;
	errp = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wpp, &eflag) ;
	err = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wp, &eflag) ;
	errm = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wpm, &eflag) ;

	/* attempt a Halley/Muller/Bailey/Press step */
	dedW = (errp - errm)/(Wpp - Wpm) ;
	dedW2 = (errp - 2.*err + errm)/(h*h) ;
	f = 0.5*err*dedW2/(dedW*dedW) ;
		/* limits size of 2nd deriv correction */
	dW = -err/dedW/(1. - MY_MIN(MY_MAX(-0.3, f),0.3)) ;
	Wp1 = Wp ;
	err1 = err ;
		/* limit size of step */
	Wp += MY_MAX( MY_MIN(dW, 2.0*Wp), -0.5*Wp) ;
	err = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wp, &eflag) ;

	/* not good enough?  apply secant method */
	for(iter = 0 ; 
		(iter < ITERMAX) && (fabs(err/Wp) > ERRTOL); 
		iter++) {
		
		dW = (Wp1 - Wp)*err / (err - err1) ;

		Wp1 = Wp ;
		err1 = err ;

		/* normal secant increment is dW.  This expression
		   limits new guess to between 2 and 0.5 times 
		   current value */
		Wp += MY_MAX( MY_MIN(dW, 2.0*Wp), -0.5*Wp) ;
		//Wp += dW;

		if(fabs(dW/Wp) < ERRTOL) break ;

		err = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wp, &eflag) ;
	}

	/* failure to converge; do not set primitives other than B */
	if(iter == ITERMAX) return(1) ;

	/* final Halley/Muller/Bailey/Press polishing step */
	/* does not significantly improve performance */
	/*
	Wpm = (1. - DEL)*Wp ; //heuristic
	h = Wp - Wpm ;
	Wpp = Wp + h ;
	errp = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wpp, &eflag) ;
	err = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wp, &eflag) ;
	errm = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wpm, &eflag) ;
	dedW = (errp - errm)/(Wpp - Wpm) ;
	dedW2 = (errp - 2.*err + errm)/(h*h) ;
	f = 0.5*err*dedW2/(dedW*dedW) ;
	dW = -err/dedW/(1. - MY_MIN(MY_MAX(-0.2, f),0.2)) ;
	Wp += MY_MAX( MY_MIN(dW, 2.0*Wp), -0.5*Wp) ;
	*/
	//err = err_eqn(Bsq, D, Ep, QdB, Qtsq, Wp, &eflag) ;
	
#if 0
	cnt++ ;
	itercnt[iter]++ ;

	/* dump the iteration distribution periodically */
	if(cnt%(2*N1*N2*256) == 0) {
		for(i=0;i<ITERMAX;i++) 
			fprintf(stderr,"%d %lf\n",i,((double)itercnt[i])/((double)cnt)) ;


		for(i=0;i<ITERMAX;i++) itercnt[i] = 0 ;
		cnt = 0 ;
	}
#endif

	/* find utsq, gamma, rho0 from Wp */
	gamma = gamma_func(Bsq, D, QdB, Qtsq, Wp, &eflag) ;
	if(gamma < 1.) {
		fprintf(stderr,"gamma < 1 failure. STOP.\n") ;
		exit(1) ;
	}

	/* find the scalars */
	rho0 = D/gamma ;
	W = Wp + D ;
	w = W/(gamma*gamma) ;
	P = Pressure_rho0_w(rho0, w);
	u = w - (rho0 + P) ;

	/** set primitives **/

	/* return w/o updating non-B primitives */
	if(rho0 < 0 || u < 0) {
		return(5) ;
	}

	prim[RHO] = rho0 ;
	prim[UU] = u ;

	/* find u(tilde); eq. 31 of Noble et al. */
	prim[U1] = (gamma/(W + Bsq))*(Qtcon[1] + QdB*Bcon[1]/W) ;
	prim[U2] = (gamma/(W + Bsq))*(Qtcon[2] + QdB*Bcon[2]/W) ;
	prim[U3] = (gamma/(W + Bsq))*(Qtcon[3] + QdB*Bcon[3]/W) ;

	/* done! */
	return(0) ;
}

double err_eqn(
	double Bsq,
	double D,
	double Ep,
	double QdB,
	double Qtsq,
	double Wp,
	int *eflag
	) 
{
	double p,err,W,w,rho0,gamma ;
	double Pressure_rho0_w(double rho0, double w) ;
	double gamma_func(double Bsq, double D, double QdB, double Qtsq,
		double Wp, int *eflag) ;

	W = Wp + D ;
	gamma = gamma_func(Bsq,D,QdB,Qtsq,Wp,eflag) ;
	w = W/(gamma*gamma) ;
	rho0 = D/gamma ;
	p = Pressure_rho0_w(rho0,w) ;

	err = -Ep + Wp - p + 0.5*Bsq + 0.5*(Bsq*Qtsq - QdB*QdB)/((Bsq + W)*(Bsq + W)) ;
		
	return(err) ;
}

double gamma_func(
	double Bsq,
	double D,
	double QdB,
	double Qtsq,
	double Wp,
	int *eflag) 
{
	double QdBsq, W, utsq, gamma, W2, WB ;

	QdBsq = QdB*QdB ;
	W = D + Wp ;
	W2 = W*W ;
	WB = W + Bsq ;

	/* this is basically inversion of eq. A7 of MM */
	utsq = -( (W + WB)*QdBsq + W2*Qtsq )/(QdBsq*(W + WB) + W2*(Qtsq - WB*WB)) ;
	gamma = sqrt(1. + fabs(utsq)) ;

	/* CATCH: utsq < 0 */
	if(utsq < 0. || utsq > 1.e3*GAMMAMAX*GAMMAMAX) {
		*eflag = 2 ;
	}
	
	return(gamma) ;
}

double Wp_func(double *prim, struct of_geom *geom, int *eflag)
{
	double rho0, u, utcon[NDIM], utcov[NDIM], utsq, gamma, Wp ;
	double Pressure_rho0_u(double rho0, double u) ;
	int i ;

	rho0 = prim[RHO] ;
	u    = prim[UU] ;

	utcon[0] = 0. ;
	utcon[1] = prim[U1] ;
	utcon[2] = prim[U2] ;
	utcon[3] = prim[U3] ;

	lower(utcon, geom, utcov) ;
	utsq = 0. ; for(i=0;i<NDIM;i++) utsq += utcon[i]*utcov[i] ;

	/* CATCH: utsq < 0 */
	if((utsq < 0.) && (fabs(utsq) < 1.e-13)) {
		utsq = fabs(utsq) ;
	}
	if(utsq < 0. || utsq > 1.e3*GAMMAMAX*GAMMAMAX) {
		*eflag = 2 ;
		return(rho0 + u) ;	/* not sure what to do here... */
	}

	gamma = sqrt(1. + fabs(utsq)) ;

	Wp = (rho0 + u + Pressure_rho0_u(rho0,u))*gamma*gamma - rho0*gamma ;

	return(Wp) ;
}

/** EOS **/
double Pressure_rho0_u(double rho0, double u)
{
	return((gam - 1.)*u) ;
}

double Pressure_rho0_w(double rho0, double w)
{
	return((w - rho0)*(gam - 1.)/gam) ;
}

/* Hotaka EOS */
#if 0
double Pressure_rho0_u(double rho0, double u)
{

}
#endif

