
#include "decs.h"

/*

special diagnostic routine to obtain current
and charge density in fluid frame

ASSUMES: p, psave, and dtsave set.

RETURNS: charge values through global array Jcon 

Does not require BCs to be set correctly; 
charge density calculated only at interior points.

Written to be readable rather than efficient.

cfg 25 june 10

- fixed bug in calculation of DF, cfg 18 june 12

*/

void current_calc()
{
    static grid_prim_type pa ;
	double Fcon_calc(double *prim, int i1, int i2, int i, int j) ;
	double gF0p[NDIM],gF0m[NDIM],gF1p[NDIM],gF1m[NDIM],gF2p[NDIM],gF2m[NDIM],gF3p[NDIM],gF3m[NDIM] ;
	struct of_geom *geom ;
	int i,j,k,d ;
	
	/* calculate a time-centered p */
	ZSLOOP(-1,N1-1+NG,-1,N2-1+NG,-1,N3-1+NG) PLOOP pa[i][j][k][ip] = 0.5*(p[i][j][k][ip] + psave[i][j][k][ip]) ;

	/* calculate J using centered differences; interior zones only */
	ZLOOP for(d=0;d<NDIM;d++) Jcon[i][j][k][d] = 0. ;
	ZLOOP {

		/* get gdet * Fmunu at neighboring points */
		/* first time direction */
		for(d = 0 ; d < NDIM ; d++) gF0p[d] = Fcon_calc(p[i][j][k],  0, d, i, j) ;
		for(d = 0 ; d < NDIM ; d++) gF0m[d] = Fcon_calc(psave[i][j][k], 0, d, i, j) ;

		/* x1 direction */
		for(d = 0 ; d < NDIM ; d++) gF1p[d] = Fcon_calc(pa[i+1][j][k], 1, d, i+1, j) ;
		for(d = 0 ; d < NDIM ; d++) gF1m[d] = Fcon_calc(pa[i-1][j][k], 1, d, i-1, j) ;

		/* x2 direction */
		for(d = 0 ; d < NDIM ; d++) gF2p[d] = Fcon_calc(pa[i][j+1][k], 2, d, i, j+1) ;
		for(d = 0 ; d < NDIM ; d++) gF2m[d] = Fcon_calc(pa[i][j-1][k], 2, d, i, j-1) ;

		/* x3 direction */
		for(d = 0 ; d < NDIM ; d++) gF3p[d] = Fcon_calc(pa[i][j][k+1], 3, d, i, j);
		for(d = 0 ; d < NDIM ; d++) gF3m[d] = Fcon_calc(pa[i][j][k-1], 3, d, i, j);

		/* get gdet at point */
		geom = get_geometry(i, j, CENT) ;

		/* difference ; use Maxwell in the form D_b F^{ab} = 4\pi J^a,
		   assuming symmetry along the 3-axis */
		for(d = 0 ; d < NDIM ; d++) {
			Jcon[i][j][k][d] = (1./(4.*M_PI*geom->g))*(
				(gF0p[d] - gF0m[d])/dtsave +
				(gF1p[d] - gF1m[d])/(2.*dx[1]) +
				(gF2p[d] - gF2m[d])/(2.*dx[2]) +
				(gF3p[d] - gF3m[d])/(2.*dx[3]) 
				) ;
		}
	}

	return ;
}

/* return single component of the contravariant maxwell tensor at position i,j
	component i1,i2, constructed from primitives prim */
double Fcon_calc(double *prim, int i1, int i2, int i, int j) 
{
	struct of_geom *geom ;
	double ucon[NDIM],ucov[NDIM],bcon[NDIM],bcov[NDIM] ;
	double Fcon,gFcon,dFcon ;
	int k,l ;
	int antisym(int a, int b, int c, int d) ;

	if(i1 == i2) return(0.) ;

	geom = get_geometry(i,j, CENT) ;
        ucon_calc(prim, geom, ucon);
        lower(ucon, geom, ucov);
        bcon_calc(prim, ucon, ucov, bcon);
        lower(bcon, geom, bcov);

	Fcon = 0. ;
	for(k=0;k<4;k++)
	for(l=0;l<4;l++) {
		dFcon = (-1./geom->g)*antisym(i1,i2,k,l)*ucov[k]*bcov[l] ;
		Fcon += dFcon ;
	}

	gFcon = Fcon * geom->g ;

	return(gFcon) ;
}


/* completely antisymmetric symbol in 4D.
   verified against mathematica */
int antisym(int a, int b, int c, int d)
{
	int pp(int n, int *P) ;

        /** check for a valid permutation **/
        /* range? */
        if(a < 0 || a > 3) return(100) ;
        if(b < 0 || b > 3) return(100) ;
        if(c < 0 || c > 3) return(100) ;
        if(d < 0 || d > 3) return(100) ;

        /* entries different? */
        if(a == b) return(0) ;
        if(a == c) return(0) ;
        if(a == d) return(0) ;
        if(b == c) return(0) ;
        if(b == d) return(0) ;
        if(c == d) return(0) ;

        /* determine parity of permutation */
        int p[4] = {a,b,c,d} ;
        return(pp(4,p)) ;
}

/* algorithm tracks back to Norm Hardy; good for general n */
int pp(int n, int P[n])
{
        int j,x ;
        int p = 0;
        int v[n] ;

        for(j=0;j<n;j++) v[j] = 0 ;

        for(j=0;j<n;j++) {
                if (v[j]) p++;
                else {
                        x = j;
                        do {
                                x = P[x];
                                v[x] = 1;
                        } while (x != j);
                }
        }

        if(p%2 == 0) return(1) ;
        else return(-1) ;
}

