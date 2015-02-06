
/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"

/* local functions */
void coord_transform(double *pr, int i, int j);
double lfish_calc(double rmax) ;

void init()
{
	int i, j;
	double r, th, sth, cth;
	double ur, uh, up, u, rho;
	double X[NDIM];
	struct of_geom *geom;

	/* for disk interior */
	double l, rin, lnh, expm2chi, up1;
	double DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
	double kappa, hm1;

	/* for magnetic field */
	static double A[N1 + 2*NG][N2 + 2*NG];
	double rho_av, rhomax, umax, beta, bsq_ij, bsq_max, norm, q,
	    beta_act;
	double rmax ;

	/* some physics parameters */
	gam = 13. / 9.;
	//gam = 4. / 3.;

	/* disk parameters (use fishbone.m to select new solutions) */
	a = 0.9375;
	rin = 6.;
	rmax = 12.;
	l = lfish_calc(rmax);

	kappa = 1.e-3;
	beta = 1.e2;

	/* some numerical parameters */
	lim = MC;
	failed = 0;		/* start slow */
	cour = 0.1;
	dtsave = dt = 1.e-8;
	R0 = 0.0;
	Rin = 0.98 * (1. + sqrt(1. - a * a));
	Rout = 40.;

	t = 0.;
	hslope = 0.3;

	zero_arrays();
	set_grid();

	coord(-2, 0, CENT, X);
	bl_coord(X, &r, &th);
	fprintf(stderr, "rmin: %g\n", r);
	fprintf(stderr, "rmin/rm: %g\n", r / (1. + sqrt(1. - a * a)));

	/* output choices */
	tf = 2000.0;

	DTd = 50.;		/* dumping frequency, in units of M */
	DTl = 0.1;		/* logfile frequency, in units of M */
	DTi = 0.1;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;
	defcon = 1.;

	rhomax = 0.;
	umax = 0.;
	ZLOOP {
		coord(i, j, CENT, X);
		bl_coord(X, &r, &th);

		sth = sin(th);
		cth = cos(th);

		/* nominal values; real value set by fixup */
		rho = 1.e-7 * RHOMIN;
		u = 1.e-7 * UUMIN;

		ur = 0.;
		uh = 0.;
		up = 0.;

		p[i][j][RHO] = rho;
		p[i][j][UU] = u;
		p[i][j][U1] = ur;
		p[i][j][U2] = uh;
		p[i][j][U3] = up;
		p[i][j][B1] = 0.;
		p[i][j][B2] = 0.;
		p[i][j][B3] = 0.;
	}

	/* enforce boundary conditions */
	fixup(p);
	bound_prim(p);

}

/* this version starts w/ BL 4-velocity and
 * converts to 3-velocities in modified
 * Kerr-Schild coordinates */

void coord_transform(double *pr, int ii, int jj)
{
	double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
	double AA, BB, CC, discr;
	double alpha, gamma, beta[NDIM];
	struct of_geom *geom, blgeom;
	int j, k;

	coord(ii, jj, CENT, X);
	bl_coord(X, &r, &th);
	blgset(ii, jj, &blgeom);

	ucon[1] = pr[U1];
	ucon[2] = pr[U2];
	ucon[3] = pr[U3];

	AA = blgeom.gcov[TT][TT];
	BB = 2. * (blgeom.gcov[TT][1] * ucon[1] +
		   blgeom.gcov[TT][2] * ucon[2] +
		   blgeom.gcov[TT][3] * ucon[3]);
	CC = 1. +
	    blgeom.gcov[1][1] * ucon[1] * ucon[1] +
	    blgeom.gcov[2][2] * ucon[2] * ucon[2] +
	    blgeom.gcov[3][3] * ucon[3] * ucon[3] +
	    2. * (blgeom.gcov[1][2] * ucon[1] * ucon[2] +
		  blgeom.gcov[1][3] * ucon[1] * ucon[3] +
		  blgeom.gcov[2][3] * ucon[2] * ucon[3]);

	discr = BB * BB - 4. * AA * CC;
	ucon[TT] = (-BB - sqrt(discr)) / (2. * AA);
	/* now we've got ucon in BL coords */

	/* transform to Kerr-Schild */
	/* make transform matrix */
	DLOOP trans[j][k] = 0.;
	DLOOPA trans[j][j] = 1.;
	trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
	trans[3][1] = a / (r * r - 2. * r + a * a);

	/* transform ucon */
	DLOOPA tmp[j] = 0.;
	DLOOP tmp[j] += trans[j][k] * ucon[k];
	DLOOPA ucon[j] = tmp[j];
	/* now we've got ucon in KS coords */

	/* transform to KS' coords */
	ucon[1] *= (1. / (r - R0));
	ucon[2] *=
	    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));

	/* now solve for v-- we can use the same u^t because
	 * it didn't change under KS -> KS' */
	geom = get_geometry(ii, jj, CENT) ;
	alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
	gamma = ucon[TT] * alpha;

	beta[1] = alpha * alpha * geom->gcon[0][1];
	beta[2] = alpha * alpha * geom->gcon[0][2];
	beta[3] = alpha * alpha * geom->gcon[0][3];

	pr[U1] = ucon[1] + beta[1] * gamma / alpha;
	pr[U2] = ucon[2] + beta[2] * gamma / alpha;
	pr[U3] = ucon[3] + beta[3] * gamma / alpha;

	/* done! */
}

double lfish_calc(double r)
{
	return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
		 ((-2. * a * r *
		   (pow(a, 2) - 2. * a * sqrt(r) +
		    pow(r,
			2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
		  ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) * 
		  (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
		/ (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
		   (pow(a, 2) + (-2. + r) * r))
	    );
}
