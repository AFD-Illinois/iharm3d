
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
	int i, j, k;
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
	cour = 0.3;
	dtsave = dt = 1.e-8;
	R0 = 0.0;
	Rhor = (1. + sqrt(1. - a*a));
	Rin = 0.98 * (1. + sqrt(1. - a * a));
	Rout = 40.;

	t = 0.;
	hslope = 0.3;

	zero_arrays();
	set_grid();

	//coord(-2, 0, CENT, X);
	//bl_coord(X, &r, &th);
	//fprintf(stderr, "rmin: %g\n", r);
	//fprintf(stderr, "rmin/rm: %g\n", r / (1. + sqrt(1. - a * a)));

	/* output choices */
	tf = 2000.0;

	DTd = 5.;		/* dumping frequency, in units of M */
	DTl = 0.5;		/* logfile frequency, in units of M */
	DTi = 2.0;		/* image file frequ., in units of M */
	DTr = 512;		/* restart file frequ., in timesteps */

	/* start diagnostic counters */
	dump_cnt = 0;
	image_cnt = 0;
	rdump_cnt = 0;
	defcon = 1.;

	rhomax = 0.;
	umax = 0.;
//	ZLOOP {
	ZSLOOP(-1, N1, -1, N2, -1, N3) {
		coord(i, j, CENT, X);
		bl_coord(X, &r, &th);

		sth = sin(th);
		cth = cos(th);

		/* calculate lnh */
		DD = r * r - 2. * r + a * a;
		AA = (r * r + a * a) * (r * r + a * a) -
		    DD * a * a * sth * sth;
		SS = r * r + a * a * cth * cth;

		thin = M_PI / 2.;

		sthin = sin(thin);
		cthin = cos(thin);
		DDin = rin * rin - 2. * rin + a * a;
		AAin = (rin * rin + a * a) * (rin * rin + a * a)
		    - DDin * a * a * sthin * sthin;
		SSin = rin * rin + a * a * cthin * cthin;

		if (r >= rin) {
			lnh =
			    0.5 *
			    log((1. +
				 sqrt(1. +
				      4. * (l * l * SS * SS) * DD / (AA *
								     sth *
								     AA *
								     sth)))
				/ (SS * DD / AA))
			    - 0.5 * sqrt(1. +
					 4. * (l * l * SS * SS) * DD /
					 (AA * AA * sth * sth))
			    - 2. * a * r * l / AA -
			    (0.5 *
			     log((1. +
				  sqrt(1. +
				       4. * (l * l * SSin * SSin) * DDin /
				       (AAin * AAin * sthin * sthin))) /
				 (SSin * DDin / AAin))
			     - 0.5 * sqrt(1. +
					  4. * (l * l * SSin * SSin) *
					  DDin / (AAin * AAin * sthin *
						  sthin))
			     - 2. * a * rin * l / AAin);
		} else
			lnh = 1.;


		/* regions outside torus */
		if (lnh < 0. || r < rin) {
			/* nominal values; real value set by fixup */
			rho = 1.e-7 * RHOMIN;
			u = 1.e-7 * UUMIN;

			/* these values are demonstrably physical
			   for all values of a and r */
			/*
			   ur = -1./(r*r) ;
			   uh = 0. ;
			   up = 0. ;
			 */

			ur = 0.;
			uh = 0.;
			up = 0.;

			/*
			   geom = get_geometry(i,j,CENT) ;
			   ur = geom->gcon[0][1]/geom->gcon[0][0] ;
			   uh = geom->gcon[0][2]/geom->gcon[0][0] ;
			   up = geom->gcon[0][3]/geom->gcon[0][0] ;
			 */

			p[i][j][k][RHO] = rho;
			p[i][j][k][UU] = u;
			p[i][j][k][U1] = ur;
			p[i][j][k][U2] = uh;
			p[i][j][k][U3] = up;
		}
		/* region inside magnetized torus; u^i is calculated in
		 * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
		 * so it needs to be transformed at the end */
		else {
			hm1 = exp(lnh) - 1.;
			rho = pow(hm1 * (gam - 1.) / (kappa * gam),
				  1. / (gam - 1.));
			u = kappa * pow(rho, gam) / (gam - 1.);
			ur = 0.;
			uh = 0.;

			/* calculate u^phi */
			expm2chi = SS * SS * DD / (AA * AA * sth * sth);
			up1 =
			    sqrt((-1. +
				  sqrt(1. + 4. * l * l * expm2chi)) / 2.);
			up = 2. * a * r * sqrt(1. +
					       up1 * up1) / sqrt(AA * SS *
								 DD) +
			    sqrt(SS / AA) * up1 / sth;


			p[i][j][k][RHO] = rho;
			if (rho > rhomax) rhomax = rho;
			p[i][j][k][UU] = u * (1. + 4.e-2 * (ranc() - 0.5));
			if (u > umax && r > rin) umax = u;
			p[i][j][k][U1] = ur;
			p[i][j][k][U2] = uh;
			p[i][j][k][U3] = up;

			/* convert from 4-vel to 3-vel */
			coord_transform(p[i][j][k], i, j);
		}

		p[i][j][k][B1] = 0.;
		p[i][j][k][B2] = 0.;
		p[i][j][k][B3] = 0.;
	}

	/* Normalize the densities so that max(rho) = 1 */
	umax = mpi_max(umax);
	rhomax = mpi_max(rhomax);

	if(mpi_io_proc()) fprintf(stderr, "rhomax: %g\n", rhomax);
	ZSLOOP(-1, N1, -1, N2, -1, N3) {
		p[i][j][k][RHO] /= rhomax;
		p[i][j][k][UU] /= rhomax;
	}
	umax /= rhomax;
	rhomax = 1.;
	fixup(p);
	bound_prim(p);

	/* first find corner-centered vector potential */
	ZSLOOP(0, N1, 0, N2, 0, 0) A[i][j] = 0.;
	ZSLOOP(0, N1, 0, N2, 0, 0) {
		/* vertical field version */
		/*
		   coord(i,j,CORN,X) ;
		   bl_coord(X,&r,&th) ;

		   A[i][j] = 0.5*r*sin(th) ;
		 */

		/* field-in-disk version */
		/* flux_ct */
		rho_av = 0.25 * (p[i][j][NG][RHO] + p[i - 1][j][NG][RHO] +
				 p[i][j - 1][NG][RHO] + p[i - 1][j - 1][NG][RHO])
				* (1. + 0.0*(ranc()-0.5)) 
				;

		q = rho_av / rhomax - 0.2;
		if (q > 0.) A[i][j] = q;
	}

	/* now differentiate to find cell-centered B,
	   and begin normalization */
	bsq_max = 0.;
	ZLOOP {
		geom = get_geometry(i, j, CENT) ;

		/* flux-ct */
		p[i][j][k][B1] = -(A[i][j] - A[i][j + 1]
				+ A[i + 1][j] - A[i + 1][j + 1]) / 
				(2. * dx[2] * geom->g);
		p[i][j][k][B2] = (A[i][j] + A[i][j + 1]
			       - A[i + 1][j] - A[i + 1][j + 1]) / 
			       (2. * dx[1] * geom->g);

		p[i][j][k][B3] = 0.;

		bsq_ij = bsq_calc(p[i][j][k], geom);
		if (bsq_ij > bsq_max) bsq_max = bsq_ij;
	}
	bsq_max = mpi_max(bsq_max);
	if(mpi_io_proc()) fprintf(stderr, "initial bsq_max: %g\n", bsq_max);

	/* finally, normalize to set field strength */
	beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
	if(mpi_io_proc()) fprintf(stderr, "initial beta: %g (should be %g)\n", beta_act,
		beta);
	norm = sqrt(beta_act / beta);
	bsq_max = 0.;
	ZLOOP {
		p[i][j][k][B1] *= norm;
		p[i][j][k][B2] *= norm;

		geom = get_geometry(i, j, CENT) ;
		bsq_ij = bsq_calc(p[i][j][k], geom);
		if (bsq_ij > bsq_max)
			bsq_max = bsq_ij;
	}
	bsq_max = mpi_max(bsq_max);
	beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
	if(mpi_io_proc()) fprintf(stderr, "final beta: %g (should be %g)\n", beta_act, beta);

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
	memset(trans, 0, 16*sizeof(double));
	DLOOPA trans[j][j] = 1.;
	trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
	trans[3][1] = a / (r * r - 2. * r + a * a);

	/* transform ucon */
	DLOOPA tmp[j] = 0.;
	DLOOP tmp[j] += trans[j][k] * ucon[k];
	DLOOPA ucon[j] = tmp[j];
	/* now we've got ucon in KS coords */

	/* transform to KS' coords */
	//ucon[1] *= (1. / (r - R0));
	ucon[1] /= dr_dx(X[1]);
	//ucon[2] *=
	//    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
	ucon[2] /= dth_dx(X[2]);

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
