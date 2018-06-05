/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR FISHBONE-MONCRIEF TORUS                             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#include <gsl/gsl_randist.h>

// Local declarations
// TODO clean this not to use old of_geom style
struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
	double alpha;
};
void coord_transform(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k);
double lfish_calc(double rmax);
void blgset(int i, int j, struct of_geom *geom);
double bl_gdet_func(double r, double th);
void bl_gcov_func(double r, double th, double gcov[][NDIM]);
void bl_gcon_func(double r, double th, double gcon[][NDIM]);


static int MAD, NORM_WITH_MAXR = 0; //TODO clean up these parameters
static double BHflux, beta;
static double rin, rmax;
static double rBstart, rBend;
void set_problem_params() {
  set_param("rin", &rin);
  set_param("rmax", &rmax);

  set_param("MAD", &MAD);
  set_param("BHflux", &BHflux);
  set_param("beta", &beta);

  set_param("rBstart", &rBstart);
  set_param("rBend", &rBend);
}

void init(struct GridGeom *G, struct FluidState *S)
{
  // Magnetic field
  double (*A)[N2 + 2*NG] = malloc(sizeof(*A) * (N1 + 2*NG));

  // Fishbone-Moncrief parameters
  double l = lfish_calc(rmax);
  double kappa = 1.e-3;

  // Grid parameters
  R0 = 0.0;
  Rhor = (1. + sqrt(1. - a*a));
  double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  Risco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));

  zero_arrays();
  set_grid(G);
  if (mpi_io_proc()) printf("grid set\n");

  // Initialize counters and such
  t = 0.; // TODO set these in main?
  failed = 0;
  dump_cnt = 0; // TODO track this in dump
  rdump_cnt = 0;

  double rhomax = 0.;
  double umax = 0.;
  ZSLOOP(-1, N3, -1, N2, -1, N1) {
    double X[NDIM];
    coord(i, j, k, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);

    double sth = sin(th);
    double cth = cos(th);

    // Calculate lnh
    double DD = r * r - 2. * r + a * a;
    double AA = (r * r + a * a) * (r * r + a * a) -
             DD * a * a * sth * sth;
    double SS = r * r + a * a * cth * cth;

    double thin = M_PI / 2.;
    double sthin = sin(thin);
    double cthin = cos(thin);

    double DDin = rin * rin - 2. * rin + a * a;
    double AAin = (rin * rin + a * a) * (rin * rin + a * a)
             - DDin * a * a * sthin * sthin;
    double SSin = rin * rin + a * a * cthin * cthin;

    double lnh;
    if (r >= rin) {
      lnh =
          0.5 *
          log((1. +
         sqrt(1. +
              4. * (l * l * SS * SS) * DD / (AA * AA * sth * sth)))
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
            4. * (l * l * SSin * SSin) * DDin / (AAin * AAin * sthin * sthin))
           - 2. * a * rin * l / AAin);
    } else {
      lnh = 1.;
    }

    // regions outside torus
    if (lnh < 0. || r < rin) {
      // Nominal values; real value set by fixup

      S->P[RHO][k][j][i] = 1.e-7 * RHOMIN;
      S->P[UU][k][j][i] = 1.e-7 * UUMIN;
      S->P[U1][k][j][i] = 0.;
      S->P[U2][k][j][i] = 0.;
      S->P[U3][k][j][i] = 0.;
    }
    /* region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
     * so it needs to be transformed at the end */
    else {
      double hm1 = exp(lnh) - 1.;
      double rho = pow(hm1 * (gam - 1.) / (kappa * gam),
               1. / (gam - 1.));
      double u = kappa * pow(rho, gam) / (gam - 1.);

      // Calculate u^phi
      double expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      double up1 =
          sqrt((-1. +
          sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      double up = 2. * a * r * sqrt(1. +
                 up1 * up1) / sqrt(AA * SS *
                 DD) +
          sqrt(SS / AA) * up1 / sth;


      S->P[RHO][k][j][i] = rho;
      if (rho > rhomax) rhomax = rho;
      //u *= (1. + 4.e-2 * (get_random() - 0.5));
      S->P[UU][k][j][i] = u;
      if (u > umax && r > rin) umax = u;
      S->P[U1][k][j][i] = 0.;
      S->P[U2][k][j][i] = 0.;
      S->P[U3][k][j][i] = up;

      // Convert from 4-velocity to 3-velocity
      coord_transform(G, S, i, j, k);
    }

    S->P[B1][k][j][i] = 0.;
    S->P[B2][k][j][i] = 0.;
    S->P[B3][k][j][i] = 0.;
  } // ZSLOOP

  // Find the zone in which rBend of Narayan condition resides
  // This just uses the farthest process in R
  // For /very/ large N1CPU it might fail
  int iend_global = 0;
  if (global_stop[0] == N1TOT && global_start[1] == 0 && global_start[2] == 0) {
    int iend = NG;
    double r_iend = 0.0;
    while (r_iend < rBend) {
      iend++;
      double Xend[NDIM];
      coord(iend,N2/2+NG,NG,CORN,Xend);
      double thend;
      bl_coord(Xend,&r_iend,&thend);
    }
    iend--;
    iend_global = global_start[0] + iend - NG; //Translate to coordinates for uu_mid below
    printf("[MPI %d] Furthest torus zone is %d (locally %d), at r = %f\n", mpi_myrank(), iend_global, iend, r_iend);
  }
  iend_global = mpi_reduce_int(iend_global); //TODO This should be a broadcast I know.

  // Normalize the densities so that max(rho) = 1
  umax = mpi_max(umax);
  rhomax = mpi_max(rhomax);

  //ZSLOOP(-1, N3, -1, N2, -1, N1) {
  ZLOOPALL {
    S->P[RHO][k][j][i] /= rhomax;
    S->P[UU][k][j][i] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;
  fixup(G, S);
  set_bounds(G, S);

  // Calculate UU along midplane, propagate to all processes
  double *uu_plane_send = calloc(N1TOT,sizeof(double));

  // This relies on an even N2TOT /and/ N2CPU
  if ((global_start[1] == N2TOT/2 || N2CPU == 1) && global_start[2] == 0) {
    int j_mid = N2TOT/2 - global_start[1] + NG;
    int k = NG; // Axisymmetric
    ILOOP {
      int i_global = global_start[0] + i - NG;
      uu_plane_send[i_global] = 0.25*(S->P[UU][k][j_mid][i] + S->P[UU][k][j_mid][i-1] +
                          S->P[UU][k][j_mid-1][i] + S->P[UU][k][j_mid-1][i-1]);
    }
  }

  double *uu_plane = calloc(N1TOT,sizeof(double));
  mpi_reduce_vector(uu_plane_send, uu_plane, N1TOT);
  free(uu_plane_send);
  //printf ("UU in plane is "); for (int i =0; i < N1TOT; i++) printf("%.10e ", uu_plane[i]);

  // first find corner-centered vector potential
  ZSLOOP(0, 0, -NG, N2+NG-1, -NG, N1+NG-1) A[i][j] = 0.;
  ZSLOOP(0, 0, -NG+1, N2+NG-1, -NG+1, N1+NG-1) {
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r, th;
    bl_coord(X,&r,&th);

    double q;

    // Field in disk
    double rho_av = 0.25*(S->P[RHO][k][j][i] + S->P[RHO][k][j][i-1] +
                    S->P[RHO][k][j-1][i] + S->P[RHO][k][j-1][i-1]);
    double uu_av = 0.25*(S->P[UU][k][j][i] + S->P[UU][k][j][i-1] +
                         S->P[UU][k][j-1][i] + S->P[UU][k][j-1][i-1]);

    int i_global = global_start[0] + i - NG;
    double uu_plane_av = uu_plane[i_global];
    double uu_end = uu_plane[iend_global];

    double b_buffer = 0.0; //Minimum rho at which there will be B field
    if (N3 > 1) {
      if (MAD == 1) { // Classic
        q = rho_av/rhomax;
      } else if (MAD == 2) { // Vertical threaded
        q = pow(sin(th),3)*pow(r/rin,3.)*exp(-r/400)*rho_av/rhomax;
      } else if (MAD == 3) { // Additional r^3 term
        q = pow(r/rin,3.)*rho_av/rhomax;
      } else if (MAD == 4) { // T,N,M (2011) uses phi = r^5 rho^2
        q = pow(r/rin,5.)*pow(rho_av/rhomax, 2);
        //OLD_MAD = 1;
      } else if (MAD == 5) { // T,N,M (2011) without renormalization step
        q = pow(r/rin,5.)*pow(rho_av/rhomax, 2);
        NORM_WITH_MAXR = 1;
      } else if (MAD == 6) { // Gaussian-strength vertical threaded field
        double wid = 2; //Radius of half-maximum. Units of rin
        q = gsl_ran_gaussian_pdf((r/rin)*sin(th), wid/sqrt(2*log(2)));
      } else if (MAD == 7) {
        // Narayan '12, Penna '12 conditions
        // Former uses rstart=25, rend=810, lam_B=25
        double uc = uu_av - uu_end;
        double ucm = uu_plane_av - uu_end;
        q = pow(sin(th),3)*(uc/(ucm+SMALL) - 0.2) / 0.8; // Note builtin buffer of 0.2
        if ( r > rBend || r < rBstart || fabs(q) > 1.e2 ) q = 0; //Exclude q outside torus and large q resulting from division by SMALL
        NORM_WITH_MAXR = 0; //?

        //if (q != 0 && th > M_PI/2-0.1 && th < M_PI/2+0.1) printf("q_mid is %.10e\n", q);
      } else {
        printf("MAD = %i not supported!\n", MAD);
        exit(-1);
      }
    } else { // TODO How about 2D?
      q = rho_av/rhomax;
    }

    // Apply floor
    q -= b_buffer;

    A[i][j] = 0.;
    //if (q == 0) printf("q is 0 at %d %d %d\n",k,j,i);
    //if (q < 0) q = -q;
    if (q > 0.) {
      if (MAD == 7) { // Narayan limit for MAD
        double lam_B = 25;
        double pre_norm = 1;
        double flux_correction = sin( 1/lam_B * (pow(r/rin,2./3) + 15./8*pow(r/rin,-2./5) - pow(rBstart/rin,2./3) - 15./8*pow(rBstart/rin,-2./5)));
        double q_mod = q*pre_norm*flux_correction;
        //printf("Correction is %.10e\n", flux_correction);
        A[i][j] = q_mod;
      } else {
        A[i][j] = q;
      }
    }
  } // ZSLOOP

  // Calculate B-field and find bsq_max
  double bsq_max = 0.;
  double beta_min = 1e100;
  ZLOOP {
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r, th;
    bl_coord(X,&r,&th);

    // Flux-ct
    S->P[B1][k][j][i] = -(A[i][j] - A[i][j + 1]
	+ A[i + 1][j] - A[i + 1][j + 1]) /
	(2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][k][j][i] = (A[i][j] + A[i][j + 1]
	     - A[i + 1][j] - A[i + 1][j + 1]) /
	     (2. * dx[1] * G->gdet[CENT][j][i]);

    S->P[B3][k][j][i] = 0.;

    get_state(G, S, i, j, k, CENT);
    double bsq_ij = bsq_calc(S, i, j, k);
    if (bsq_ij > bsq_max) bsq_max = bsq_ij;
    if (r > rBstart && r < rBend) {
      double beta_ij = (gam - 1.)*(S->P[UU][k][j][i])/(0.5*(bsq_ij+SMALL)) ;
      if(beta_ij < beta_min) beta_min = beta_ij ;
    }
  }
  bsq_max = mpi_max(bsq_max);
  beta_min = mpi_min(beta_min);

  double umax_plane = 0;
  for (int i; i < N1TOT; i++) {
      double X[NDIM];
      coord(i,NG,NG,CORN,X);
      double r, th;
      bl_coord(X,&r,&th);
      if (r > rBstart && r < rBend) {
	if (uu_plane[i] > umax_plane) umax_plane = uu_plane[i];
      }
  }

  double norm = 0;
  if (!NORM_WITH_MAXR) {
    // Ratio of max UU, beta
    //double beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
    double beta_act = (gam - 1.) * umax_plane / (0.5 * bsq_max);
    if (mpi_io_proc()) printf("Umax is %f, bsq_max is %f, beta is %f\n", umax_plane, bsq_max, beta_act);
    norm = sqrt(beta_act / beta);
  } else {
    // Beta_min = 100 normalization
    if (mpi_io_proc()) printf("Min beta in torus is %f\n", beta_min);
    norm = sqrt(beta_min / beta) ;
  }

  // Apply normalization
  if (mpi_io_proc()) printf("Normalization is %f\n", norm);
  ZLOOP {
    S->P[B1][k][j][i] *= norm ;
    S->P[B2][k][j][i] *= norm ;
  }

  // This adds a central flux based on specifying some BHflux
  // Initialize a net magnetic field inside the initial torus
  ZSLOOP(0, 0, 0, N2, 0, N1) {
    double X[NDIM];
    coord(i,j,k,CORN,X);
    double r,th;
    bl_coord(X, &r, &th);

    A[i][j] = 0.;

    double x = r*sin(th);
    double z = r*cos(th);
    double a_hyp = 20.;
    double b_hyp = 60.;
    double x_hyp = a_hyp*sqrt(1. + pow(z/b_hyp,2));

    double q = (pow(x,2) - pow(x_hyp,2))/pow(x_hyp,2);
    if (x < x_hyp) {
      A[i][j] = 10.*q;
    }
  }

  // Evaluate net flux
  double Phi_proc = 0.;
  ISLOOP(5, N1-1) {
    JSLOOP(0, N2-1) {
      int jglobal = j - NG + global_start[1];
      //int j = N2/2+NG;
      int k = NG;
      if (jglobal == N2TOT/2) {
	double X[NDIM];
	coord(i, j, k, CENT, X);
	double r,th;
	bl_coord(X, &r, &th);

	if (r < rin) {
	  double B2net =  (A[i][j] + A[i][j+1] - A[i+1][j] - A[i+1][j+1]);
	    // / (2.*dx[1]*G->gdet[CENT][j][i]);
	  Phi_proc += fabs(B2net)*M_PI/N3CPU; // * 2.*dx[1]*G->gdet[CENT][j][i]
	}
      }
    }
  }

  //If left bound in X1.  Note different convention from bhlight!
  if (global_start[0] == 0) {
    JSLOOP(0, N2/2-1) {
      int i = 5 + NG;

      double B1net = -(A[i][j] - A[i][j+1] + A[i+1][j] - A[i+1][j+1]); // /(2.*dx[2]*G->gdet[CENT][j][i]);
      Phi_proc += fabs(B1net)*M_PI/N3CPU;  // * 2.*dx[2]*G->gdet[CENT][j][i]
    }
  }
  double Phi = mpi_reduce(Phi_proc);

  norm = BHflux/(Phi + SMALL);

  ZLOOP {
    // Flux-ct
    S->P[B1][k][j][i] += -norm*(A[i][j] - A[i][j + 1]
				+ A[i + 1][j] - A[i + 1][j + 1]) /
      (2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][k][j][i] += norm*(A[i][j] + A[i][j + 1]
			       - A[i + 1][j] - A[i + 1][j + 1]) /
      (2. * dx[1] * G->gdet[CENT][j][i]);
  }


#if ELECTRONS
  init_electrons();
#endif

  // Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);

  if(mpi_io_proc()) fprintf(stderr, "Finished init()\n");

}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k)
{
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct blgeom;
  struct of_geom blgeom;

  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);
  blgset(i, j, &blgeom);

  ucon[1] = S->P[U1][k][j][i];
  ucon[2] = S->P[U2][k][j][i];
  ucon[3] = S->P[U3][k][j][i];

  AA = blgeom.gcov[0][0];
  BB = 2.*(blgeom.gcov[0][1]*ucon[1] +
           blgeom.gcov[0][2]*ucon[2] +
           blgeom.gcov[0][3]*ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1]*ucon[1]*ucon[1] +
      blgeom.gcov[2][2]*ucon[2]*ucon[2] +
      blgeom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(blgeom.gcov[1][2]*ucon[1]*ucon[2] +
          blgeom.gcov[1][3]*ucon[1]*ucon[3] +
          blgeom.gcov[2][3]*ucon[2]*ucon[3]);

  discr = BB*BB - 4.*AA*CC;
  ucon[0] = (-BB - sqrt(discr))/(2.*AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a);
  trans[3][1] = a/(r*r - 2.*r + a*a);

  // Transform ucon
  for (int mu = 0; mu < NDIM; mu++) {
    tmp[mu] = 0.;
  }
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) {
    ucon[mu] = tmp[mu];
  }

  // This is ucon in KS coords

  // Transform to KS' coords
  //ucon[1] *= (1. / (r - R0));
  ucon[1] /= dr_dx(X[1]);
  //ucon[2] *=
  //    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
  ucon[2] /= dth_dx(X[2]);

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  //geom = get_geometry(ii, jj, CENT);
  alpha = G->lapse[CENT][j][i];
  //alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
  gamma = ucon[0]*alpha;

  beta[1] = alpha*alpha*G->gcon[CENT][0][1][j][i];//geom->gcon[0][1];
  beta[2] = alpha*alpha*G->gcon[CENT][0][2][j][i];//geom->gcon[0][2];
  beta[3] = alpha*alpha*G->gcon[CENT][0][3][j][i];//geom->gcon[0][3];

  S->P[U1][k][j][i] = ucon[1] + beta[1]*gamma/alpha;
  S->P[U2][k][j][i] = ucon[2] + beta[2]*gamma/alpha;
  S->P[U3][k][j][i] = ucon[3] + beta[3]*gamma/alpha;
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

// Boyer-Lindquist metric functions
void blgset(int i, int j, struct of_geom *geom)
{
  double r, th, X[NDIM];

  coord(i, j, 0, CENT, X);
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
    (1. + 0.5 * (a2 / r2) * (1. + cos(2. * th))));
}

void bl_gcov_func(double r, double th, double gcov[NDIM][NDIM])
{
  double sth, cth, s2, a2, r2, DD, mu;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      gcov[mu][nu] = 0.;
    }
  }

  sth = fabs(sin(th));
  s2 = sth*sth;
  cth = cos(th);
  a2 = a*a;
  r2 = r*r;
  DD = 1. - 2./r + a2/r2;
  mu = 1. + a2*cth*cth/r2;

  gcov[0][0] = -(1. - 2./(r*mu));
  gcov[0][3]  = -2.*a*s2/(r*mu);
  gcov[3][0]  = gcov[0][3];
  gcov[1][1]   = mu/DD;
  gcov[2][2]   = r2*mu;
  gcov[3][3]   = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu));

}

void bl_gcon_func(double r, double th, double gcon[NDIM][NDIM])
{
  double sth, cth, a2, r2, r3, DD, mu;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      gcon[mu][nu] = 0.;
    }
  }

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

  a2 = a*a;
  r2 = r*r;
  r3 = r2*r;
  DD = 1. - 2./r + a2/r2;
  mu = 1. + a2*cth*cth/r2;

  gcon[0][0] = -1. - 2.*(1. + a2/r2)/(r*DD*mu);
  gcon[0][3] = -2.*a/(r3*DD*mu);
  gcon[3][0] = gcon[0][3];
  gcon[1][1] = DD/mu;
  gcon[2][2] = 1./(r2*mu);
  gcon[3][3] = (1. - 2./(r*mu))/(r2*sth*sth*DD);
}
