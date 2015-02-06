
/* M1: no changes needed */

/***

includes all interpolation routines 
linear (MC or other limiter)
parabolic (from collela and woodward)
weno

***/

#include "decs.h"

/* performs the slope-limiting for the numerical flux calculation */

double slope_lim(double y1, double y2, double y3)
{
	double Dqm, Dqp, Dqc, s;

	/* woodward, or monotonized central, slope limiter */
	if (lim == MC) {
		Dqm = 2. * (y2 - y1);
		Dqp = 2. * (y3 - y2);
		Dqc = 0.5 * (y3 - y1);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else {
			if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
				return (Dqm);
			else if (fabs(Dqp) < fabs(Dqc))
				return (Dqp);
			else
				return (Dqc);
		}
	}
	/* van leer slope limiter */
	else if (lim == VANL) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else
			return (2. * s / (Dqm + Dqp));
	}

	/* minmod slope limiter (crude but robust) */
	else if (lim == MINM) {
		Dqm = (y2 - y1);
		Dqp = (y3 - y2);
		s = Dqm * Dqp;
		if (s <= 0.)
			return 0.;
		else if (fabs(Dqm) < fabs(Dqp))
			return Dqm;
		else
			return Dqp;
	}

	fprintf(stderr, "unknown slope limiter\n");
	exit(10);

	return (0.);
}

void linear_mc(double x1, double x2, double x3, double *lout, double *rout) 
{
	double Dqm,Dqp,Dqc,s;

	Dqm = 2. * (x2 - x1);
	Dqp = 2. * (x3 - x2);
	Dqc = 0.5 * (x3 - x1);

	s = Dqm * Dqp;

	if (s <= 0.)
		s = 0.;
	else {
		if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
			s = Dqm;
		else if (fabs(Dqp) < fabs(Dqc))
			s = Dqp;
		else
			s = Dqc;
	}

	/* reconstruct left, right */
	*lout = x2 - 0.5*s;
	*rout = x2 + 0.5*s;
}

/*
 * parabolic interpolation subroutin  
 * ref. Colella && Woodward's PPM paper
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 * 
 */
/* author: Xiaoyue Guan */

void para(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
         int i ;
         double y[5], dq[5];
         double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

         y[0]=x1;
         y[1]=x2;
         y[2]=x3;
         y[3]=x4;
         y[4]=x5;

         /*CW1.7 */
         for(i=1 ; i<4 ; i++) {
               Dqm = 2. *(y[i]-y[i-1]);
               Dqp = 2. *(y[i+1]-y[i]);
               Dqc = 0.5 *(y[i+1]-y[i-1]);
               aDqm = fabs(Dqm) ;
               aDqp = fabs(Dqp) ;
               aDqc = fabs(Dqc) ;
               s = Dqm*Dqp;

               if (s <=0.) dq[i]=0.;       //CW1.8
               else dq[i]=MY_MIN(aDqc,MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
         }

         /* CW1.6 */
         l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
         r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;

         qa=(r-y[2])*(y[2]-l);
         qd=(r-l);
         qe=6.0*(y[2]-0.5*(l+r));

         if (qa <=0. ) {
                l=y[2];
                r=y[2];
         }

         if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
         else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;

         *lout=l;   //a_L,j
         *rout=r;
}

/*
 * left and right state reconsitruction using WENO-5,
 * all numbers from WHAM paper
 *
 * using zone-centered value of n=5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 */

/* author: Monika Moscibrodzka */
/* modified and sign error corrected by CFG 12.28.13 */

void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{

    /* based on shu scholarpedia article, eqs 1,2,3 */
    /* 3rd order interpolations */
    double vr[3],vl[3];
    vr[0] =  (3./8.)*x1-(5./4.)*x2+(15./8.)*x3;
    vr[1] = (-1./8.)*x2+(3./4.)*x3+(3./8.)*x4;
    vr[2] =  (3./8.)*x3+(3./4.)*x4-(1./8.)*x5;

    vl[0] =  (3./8.)*x5-(5./4.)*x4+(15./8.)*x3;
    vl[1] = (-1./8.)*x4+(3./4.)*x3+(3./8.)*x2;
    vl[2] =  (3./8.)*x3+(3./4.)*x2-(1./8.)*x1;

    /* smoothness indicators, from Tchekh et al. A18, equiv to Shu eq. 8 */
    double beta[3];
    /*
    beta[0] = (4*x1*x1)/3.-(19*x1*x2)/3.+(25*x2*x2)/3.+(11*x1*x3)/3.-(31*x2*x3)/3.+(10*x3*x3)/3.;
    beta[1]= (4*x2*x2)/3.-(13*x2*x3)/3.+(13*x3*x3)/3.+(5*x2*x4)/3.-(13*x3*x4)/3.+(4*x4*x4)/3.;
    beta[2] = (4*x5*x5)/3.-(19*x5*x4)/3.+(25*x4*x4)/3.+(11*x5*x3)/3.-(31*x4*x3)/3.+(10*x3*x3)/3.;
    */
    beta[0]=(13./12.)*pow(x1-2.*x2+x3,2)+(1./4.)*pow(x1-4.*x2+3.*x3,2);
    beta[1]=(13./12.)*pow(x2-2.*x3+x4,2)+(1./4.)*pow(x4-x2,2);
    beta[2]=(13./12.)*pow(x3-2.*x4+x5,2)+(1./4.)*pow(x5-4.*x4+3.*x3,2);
    
    /* nonlinear weights, after shu eq. 9 */
    double den,wtr[3],Wr,wr[3],wtl[3],Wl,wl[3],eps;
    eps=1.e-26;

    den = eps+beta[0]; den *= den; wtr[0] = (1./16.)/den;
    den = eps+beta[1]; den *= den; wtr[1] = (5./8.)/den;
    den = eps+beta[2]; den *= den; wtr[2] = (5./16.)/den;
    Wr = wtr[0]+wtr[1]+wtr[2];
    wr[0] = wtr[0]/Wr ;
    wr[1] = wtr[1]/Wr ;
    wr[2] = wtr[2]/Wr ;

    den = eps+beta[2]; den *= den; wtl[0] = (1./16.)/den;
    den = eps+beta[1]; den *= den; wtl[1] = (5./8.)/den;
    den = eps+beta[0]; den *= den; wtl[2] = (5./16.)/den;
    Wl = wtl[0]+wtl[1]+wtl[2];
    wl[0] = wtl[0]/Wl ;
    wl[1] = wtl[1]/Wl ;
    wl[2] = wtl[2]/Wl ;

    *lout = vl[0]*wl[0]+vl[1]*wl[1]+vl[2]*wl[2];
    *rout = vr[0]*wr[0]+vr[1]*wr[1]+vr[2]*wr[2];
}


/***

Reconstruction routines for linear, parabolic, and weno

***/

void reconstruct_lr_lin(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
        double dqtmp[NMAX+2*NG][NPR] ;
        int i;

	ISLOOP(-1,N) PLOOP dqtmp[i][ip] = 0. ;

        /* calculate slopes */
	ISLOOP(-1,N) PLOOP dqtmp[i][ip] = slope_lim(ptmp[i-1][ip],ptmp[i][ip],ptmp[i+1][ip]) ;

        /* reconstruct left */
	ISLOOP(0,N) PLOOP p_l[i][ip] = ptmp[i][ip] - 0.5*dqtmp[i][ip];

        /* reconstruct right */
	ISLOOP(-1,N-1) PLOOP p_r[i][ip] = ptmp[i][ip] + 0.5*dqtmp[i][ip];
}


void reconstruct_lr_par(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
	int i ;

	ISLOOP(-1,N) {
		PLOOP {
			para(   ptmp[i-2][ip], 
				ptmp[i-1][ip], 
				ptmp[i][ip], 
				ptmp[i+1][ip], 
				ptmp[i+2][ip],
				&p_l[i][ip], 
				&p_r[i][ip]) ;
		}
	}
}

void reconstruct_lr_weno(double ptmp[NMAX+2*NG][NPR], int N, 
	double p_l[NMAX+2*NG][NPR], double p_r[NMAX+2*NG][NPR]) 
{
	int i ;

	ISLOOP(-1,N) {
		PLOOP {
			weno(   ptmp[i-2][ip], 
				ptmp[i-1][ip], 
				ptmp[i][ip], 
				ptmp[i+1][ip], 
				ptmp[i+2][ip],
				&p_l[i][ip], 
				&p_r[i][ip]) ;
		}
	}
}


