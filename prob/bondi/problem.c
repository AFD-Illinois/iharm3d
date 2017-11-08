/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR BONDI INFLOW                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Rootfinding for analytic Bondi solution
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

void set_problem_params()
{
}

double C1, C2, T_bondi[N1+2*NG][N2+2*NG], n;

// TODO TODO TODO
double Reh = 1.0;

struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
	double alpha;
};
struct params {
  double r, mdot, n, C1, C2;
};

// Local functions
void coord_transform(double *Pr, struct GridGeom *G, int i, int j);
double lfish_calc(double rmax);
// Boyer-Lindquist metric functions
void blgset(int i, int j, struct of_geom *geom);
double bl_gdet_func(double r, double th);
void bl_gcov_func(double r, double th, double gcov[][NDIM]);
void bl_gcon_func(double r, double th, double gcon[][NDIM]);
void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM]);

// Old looping mechanics
#define DLOOP1 for (int mu = 0; mu < NDIM; mu++)
#define DLOOP2 for (int mu = 0; mu < NDIM; mu++) \
               for (int nu = 0; nu < NDIM; nu++)

double rfunc(double T, void *params)
{
  struct params *p = (struct params *) params;
  //double r = p->r;
  double r = MY_MAX(2.5, p->r); // Solution breaks inside event horizon
  double mdot = p->mdot;
  double n = p->n;
  double C1 = p->C1;
  double C2 = p->C2;

  double resid = (pow(1. + (1. + n)*T,2)*(1. - 2.*mdot/r + pow(C1,2)/
         (pow(r,4)*pow(T,2*n))) - C2);

  return resid;
}

// Adapted from M. Chandra
double get_Tfunc(double T, double r)
{
  return pow(1.+(1.+n)*T,2.)*(1.-2./r+pow(C1/r/r/pow(T,n),2.))-C2;
}

double get_T(double r)
{
  double rtol = 1.e-12;
  double ftol = 1.e-14;
  double Tmin = 0.6*(sqrt(C2) - 1.)/(n + 1);
  double Tmax = pow(C1*sqrt(2./r/r/r),1./n);
  double f0, f1, fh;
  double T0, T1, Th;
  T0 = 0.6*Tmin;
  f0 = get_Tfunc(T0, r);
  T1 = Tmax;
  f1 = get_Tfunc(T1, r);

  if (f0*f1 > 0.) {
    printf("Failed solving for T at r = %e C1 = %e C2 = %e\n", r, C1, C2);
    exit(-1);
  }

  Th = (f1*T0 - f0*T1)/(f1 - f0);
  fh = get_Tfunc(Th, r);
  double epsT = rtol*(Tmin + Tmax);
  while (fabs(Th - T0) > epsT && fabs(Th - T1) > epsT && fabs(fh) > ftol) {
    if (fh*f0 < 0.) {
      T0 = Th;
      f0 = fh;
    } else {
      T1 = Th;
      f1 = fh;
    }

    Th = (f1*T0 - f0*T1)/(f1 - f0);
    fh = get_Tfunc(Th, r);
  }

  return Th;
}

void bl_to_ks(double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double trans[NDIM][NDIM];
  DLOOP2 trans[mu][nu] = 0.;
  DLOOP1 trans[mu][mu] = 1.;
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a);
  trans[3][1] = a/(r*r - 2.*r + a*a);

  DLOOP1 ucon_ks[mu] = 0.;
  DLOOP2 ucon_ks[mu] += trans[mu][nu]*ucon_bl[nu];
}

void fourvel_to_prim(double ucon[NDIM], GridPrim P,
  struct GridGeom *G, int i, int j, int k)
{
  double alpha, beta[NDIM], gamma;

  alpha = 1.0/sqrt(-G->gcon[CENT][0][0][i][j]);
  beta[1] = alpha*alpha*G->gcon[CENT][0][1][i][j];
  beta[2] = alpha*alpha*G->gcon[CENT][0][2][i][j];
  beta[3] = alpha*alpha*G->gcon[CENT][0][3][i][j];
  gamma = ucon[0]*alpha;

  P[U1][k][j][i] = ucon[1] + beta[1]*gamma/alpha;
  P[U2][k][j][i] = ucon[2] + beta[2]*gamma/alpha;
  P[U3][k][j][i] = ucon[3] + beta[3]*gamma/alpha;
}

void set_ut(double ucon[NDIM], struct of_geom *geom)
{
  double AA, BB, CC;

  AA = geom->gcov[0][0];
  BB = 2.*(geom->gcov[0][1]*ucon[1] +
           geom->gcov[0][2]*ucon[2] +
           geom->gcov[0][3]*ucon[3]);
  CC = 1. + geom->gcov[1][1]*ucon[1]*ucon[1] +
            geom->gcov[2][2]*ucon[2]*ucon[2] +
            geom->gcov[3][3]*ucon[3]*ucon[3] +
       2. *(geom->gcov[1][2]*ucon[1]*ucon[2] +
            geom->gcov[1][3]*ucon[1]*ucon[3] +
            geom->gcov[2][3]*ucon[2]*ucon[3]);

  double discr = BB*BB - 4.*AA*CC;
  ucon[0] = (-BB - sqrt(discr))/(2.*AA);
}

void get_prim_bondi(int i, int j, int k, GridPrim P, struct GridGeom *G)
{
  double r, th, X[NDIM];
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);

  while (r < Reh) {
    i++;
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);
  }

  double T = T_bondi[i][j];
  double ur = -C1/(pow(T,n)*pow(r,2));
  double rho = pow(T,n);
  double u = rho*T/(gam - 1.);
  double ucon_bl[NDIM], ucon_ks[NDIM], ucon_mks[NDIM];
  struct of_geom geom_bl;

  blgset(i, j, &geom_bl);

  DLOOP1 {
    ucon_bl[mu] = 0.;
    ucon_ks[mu] = 0.;
    ucon_mks[mu] = 0.;
  }
  ucon_bl[1] = ur;

  set_ut(ucon_bl, &geom_bl);
  bl_to_ks(X, ucon_bl, ucon_ks);

  double dxdX[NDIM][NDIM], dXdx[NDIM][NDIM];
  set_dxdX(X, dxdX);
  invert(&dxdX[0][0], &dXdx[0][0]);
  DLOOP2 {
    ucon_mks[mu] += dXdx[mu][nu]*ucon_ks[nu];
  }

  //geom_mks = get_geometry(i, j, k, CENT);
  fourvel_to_prim(ucon_mks, P, G, i, j, k);

  P[RHO][k][j][i] = rho;
  P[UU][k][j][i] = u;
  P[B1][k][j][i] = 0.;
  P[B2][k][j][i] = 0.;
  P[B3][k][j][i] = 0.;

}

void init(struct GridGeom *G, struct FluidState *S)
{
  double r, th, X[NDIM];

  tf = 100.0000;
  dt = 1.000000e-06;
  Rout = 2.000000e+01;
  gam = 1.666667e+00;
  cour = 4.000000e-01;

  DTd = 1.;
  DTl = 5.000000e-01;
  DTr = 10000;
  DTp = 10;

  a = 0.000000e+00;
  hslope = 3.000000e-01;

  double mdot = 1.;
  double rs = 8.;
  n = 1./(gam - 1.);

  // Solution constants
  double uc = sqrt(mdot/(2.*rs));
  double Vc = -sqrt(pow(uc,2)/(1. - 3.*pow(uc,2)));
  double Tc = -n*pow(Vc,2)/((n + 1.)*(n*pow(Vc,2) - 1.));
  C1 = uc*pow(rs,2)*pow(Tc,n);
  C2 = pow(1. + (1. + n)*Tc,2)*(1. - 2.*mdot/rs + pow(C1,2)/
       (pow(rs,4)*pow(Tc,2*n)));

  printf("a = %e Reh = %e\n", a, Reh);

  printf("mdot = %e\n", mdot);
  printf("rs   = %e\n", rs);
  printf("n    = %e\n", n);
  printf("uc   = %e\n", uc);
  printf("Vc   = %e\n", Vc);
  printf("Tc   = %e\n", Tc);
  printf("C1   = %e\n", C1);
  printf("C2   = %e\n", C2);

  ZSLOOP(-NG,N1+NG-1,-NG,N2+NG-1,-NG,N3+NG-1) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    T_bondi[i][j] = get_T(r);
  }

  ZLOOP {
    get_prim_bondi(i, j, k, S->P, G);
  } // ZLOOP
}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(double *Pr, struct GridGeom *G, int ii, int jj)
{
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct of_geom blgeom;

  coord(ii, jj, 0, CENT, X);
  bl_coord(X, &r, &th);
  blgset(ii, jj, &blgeom);

  ucon[1] = Pr[U1];
  ucon[2] = Pr[U2];
  ucon[3] = Pr[U3];

  AA = blgeom.gcov[0][0];
  BB = 2. * (blgeom.gcov[0][1] * ucon[1] +
       blgeom.gcov[0][2] * ucon[2] +
       blgeom.gcov[0][3] * ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1] * ucon[1] * ucon[1] +
      blgeom.gcov[2][2] * ucon[2] * ucon[2] +
      blgeom.gcov[3][3] * ucon[3] * ucon[3] +
      2. * (blgeom.gcov[1][2] * ucon[1] * ucon[2] +
      blgeom.gcov[1][3] * ucon[1] * ucon[3] +
      blgeom.gcov[2][3] * ucon[2] * ucon[3]);

  discr = BB * BB - 4. * AA * CC;
  ucon[0] = (-BB - sqrt(discr)) / (2. * AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2. * r / (r * r - 2. * r + a * a);
  trans[3][1] = a / (r * r - 2. * r + a * a);

  // Transform from BL to KS coordinates
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  // Transform from KS to MKS coordinates
  set_dxdX(X, trans);
  for (int mu = 0; mu < NDIM; mu++) tmp[mu] = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
     tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) ucon[mu] = tmp[mu];

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  //geom = get_geometry(ii, jj, 0, CENT) ;
  alpha = 1.0 / sqrt( -G->gcon[CENT][0][0][ii][jj] ) ;
  gamma = ucon[0] * alpha;

  beta[1] = alpha * alpha * G->gcon[CENT][0][1][ii][jj];
  beta[2] = alpha * alpha * G->gcon[CENT][0][2][ii][jj];
  beta[3] = alpha * alpha * G->gcon[CENT][0][3][ii][jj];

  Pr[U1] = ucon[1] + beta[1] * gamma / alpha;
  Pr[U2] = ucon[2] + beta[2] * gamma / alpha;
  Pr[U3] = ucon[3] + beta[3] * gamma / alpha;
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

void bound_gas_prob_x1r(int i, int j, int k, GridPrim  P, struct GridGeom *G)
{
  get_prim_bondi(i, j, k, P, G);
}

void set_dxdX(double X[NDIM], double dxdX[NDIM][NDIM])
{
  memset(dxdX, 0, NDIM*NDIM*sizeof(double));
  #if METRIC == MINKOWSKI
  for (int mu = 0; mu < NDIM; mu++) {
    dxdX[mu][mu] = 1.;
  }
  #elif METRIC == MKS
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  #if DEREFINE_POLES
  dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
    M_PI/2. -
    M_PI*X[2] +
    poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
    1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
    );
  dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) +
    exp(mks_smooth*(startx[1]-X[1]))*(
      -M_PI +
      2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
      (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) -
      (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
      );
  #else
  dxdX[2][2] = M_PI - (hslope - 1.)*M_PI*cos(2.*M_PI*X[2]);
  #endif
  dxdX[3][3] = 1.;
  #endif
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

