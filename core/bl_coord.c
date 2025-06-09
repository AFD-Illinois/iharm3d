/*
 * Utility functions for Boyer-Lindquist coordinates
 * Provided for problem setups but not otherwise used in core functions
 */

// TODO cleanup/minimize this file, it duplicates some of coord.c


#define met_rn       0
#define met_kerrsen  0
#define met_hayward  0
#define met_bardeen  0
#define met_jmn1     1 
#define met_jnw      0 
#define met_sv       0
#define met_deform   0
#define met_dilaton  0
#define met_lamy     0 

#define def_Rb      2.8
#define def_Q       0.
#define def_epsilon 0.
#define def_l1      0.
#define def_l2      0.
#define def_l3      0.
#define def_gamma   0.95
#define def_d       0.


#include "bl_coord.h"

void get_coefficients(double r, double th, double gvec[3]){

  gvec[0] = (1. - 2./r );
  gvec[1] = 1.;
  gvec[2] = r*r;

    if (met_jmn1){
      double def_Rbb = log(def_Rb);
      double sigg = 1./(def_Rbb - 1.);
      double rb = (1.-sigg)*def_Rbb;
      gvec[0] =  (1. - 2./r );
      gvec[1] =  1.;
      gvec[2] =  r*r;

      if(r<def_Rbb){
        gvec[0] =  (1. - 2./def_Rbb)*pow((r/rb), 2*sigg);
        gvec[1] =  1.;
        gvec[2] =  def_Rbb*def_Rb*pow(r/rb, 2. - 2.*sigg);
        fprintf(stderr ,"here: %f\n", gvec[0]);
      }
  
    }
    else if  (met_rn){
      gvec[0] = (1. - 2./r + def_Q*def_Q/r/r);
      gvec[1] =  1.;
      gvec[2] =  r*r;

    }

    else if (met_deform){

      double hr = def_epsilon/r/r/r;
      gvec[0] = (1. - 2./r)*(1. + hr);
      gvec[1] = (1. + hr)*(1. + hr);
      gvec[2] =  r*r;
    }

    else if (met_jnw){

      gvec[0] = pow((1. - 2./r/def_gamma), def_gamma);
      gvec[1] = 1.;
      gvec[2] = pow((1. - 2./r/def_gamma), 1.-def_gamma)*r*r;
    }

    else if (met_sv){

      gvec[0] = (1. - 2./sqrt(r*r + def_l1*def_l1));
      gvec[1] = 1.;
      gvec[2] = r*r + def_l1*def_l1;
    }
    else if (met_bardeen){
     
      gvec[0] = (1. - 2.*r*r/pow((r*r + def_l2*def_l2),1.5));
      gvec[1] = 1.;
      gvec[2] = r*r;
    }
    else if (met_hayward){
      double Mr = r*r*r/(r*r*r + 2.*def_l3*def_l3);
      gvec[0] =  (1. - 2.*Mr/r);
      gvec[1] =  1.;
      gvec[2] =  r*r;
    }
    else if (met_dilaton){
      double t1 = 4*r*r + def_d*def_d*def_d*def_d;
      double metG = 4.*r*r/t1;

      gvec[0] = (1. - (sqrt(t1) - def_d*def_d)/r/r);
      gvec[1] = metG;
      gvec[2] = r*r;
    }
    else if (met_lamy){

      double Mr = r*r*r/(r*r*r + 2.*def_l3*def_l3);

      if (r<0){
        Mr = fabs(r*r*r)/(fabs(r*r*r) + 2.*def_l3*def_l3);
      }
      gvec[0] =  (1. - 2.*Mr/r);
      gvec[1] =  1.;
      gvec[2] =  r*r;
      
    }
    else if (met_kerrsen){
      gvec[0] =  1.- 2./r;
      gvec[1] =  1.;
      gvec[2] =  r*(r - def_Q*def_Q);
      
    }

}

void bl_gcov_func(double r, double th, double gcov[NDIM][NDIM])
{
  DLOOP2 gcov[mu][nu] = 0.;

  double sth, cth, s2, a2, r2, DD, mu;
  sth = fabs(sin(th));
  s2 = sth*sth;
  cth = cos(th);
  // a2 = a*a;
  // r2 = r*r;
  // DD = 1. - 2./r + a2/r2;
  // mu = 1. + a2*cth*cth/r2;
  // gcov[0][0] = -(1. - 2./(r*mu));
  // gcov[0][3]  = -2.*a*s2/(r*mu);
  // gcov[3][0]  = gcov[0][3];
  // gcov[1][1]   = mu/DD;
  // gcov[2][2]   = r2*mu;
  // gcov[3][3]   = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu));

  double gtt, gg, gthth;
  double gvec[3];
  get_coefficients(r, th, gvec);
  gtt = gvec[0];
  gg = gvec[1];
  gthth = gvec[2];

  double twoF = (1./sqrt(gg) - (gtt/gg))*gthth;
  double Delta = gthth*(gtt/gg) + a*a;
  double Sigma = gthth/sqrt(gg) + a*a*cth*cth;
  double PII = (gthth/sqrt(gg) + a*a)*(gthth/sqrt(gg) + a*a) - Delta*a*a*s2;

  gcov[0][0] = -(1. - twoF/Sigma);
  gcov[0][3] = -twoF/Sigma*a*s2;
  gcov[3][0] = -twoF/Sigma*a*s2;
  gcov[1][1] = Sigma/Delta;
  gcov[2][2] = Sigma;
  gcov[3][3] = PII/Sigma*s2;

}
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
  double gcov[4][4];
  bl_gcov_func(r, M_PI/2., gcov);

  double gtt = gcov[0][0];
  double gphph = gcov[3][3];
  return sqrt(fabs(gcov[0][3]*gcov[0][3] - gtt * gphph)); 
}

void bl_gcon_func(double r, double th, double gcon[NDIM][NDIM])
{
  double sth, cth, a2, r2, r3, DD, mu;
  double gtt, gtph, gphph, grr, gthth;
  double gcov[NDIM][NDIM];
  bl_gcov_func(r, th, gcov);
  DLOOP2 gcon[mu][nu] = 0.;

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

  // a2 = a*a;
  // r2 = r*r;
  // r3 = r2*r;
  // DD = 1. - 2./r + a2/r2;
  // mu = 1. + a2*cth*cth/r2;

  // gcon[0][0] = -1. - 2.*(1. + a2/r2)/(r*DD*mu);
  // gcon[0][3] = -2.*a/(r3*DD*mu);
  // gcon[3][0] = gcon[0][3];
  // gcon[1][1] = DD/mu;
  // gcon[2][2] = 1./(r2*mu);
  // gcon[3][3] = (1. - 2./(r*mu))/(r2*sth*sth*DD);

    gtt   = gcov[0][0];
    grr   = gcov[1][1];
    gthth = gcov[2][2];
    gtph  = gcov[0][3];
    gphph = gcov[3][3];
  
    double Delta = gtt * gphph - gtph * gtph;
    gcon[0][0] =  gphph / Delta;    
    gcon[0][3] = -gtph  / Delta; 
    gcon[3][0] =  gcon[0][3];      
    gcon[3][3] =  gtt   / Delta; 
    gcon[1][1] =  1.0 / grr;    
    gcon[2][2] =  1.0 / gthth;       

}
void bl_to_ks(double X[NDIM], double ucon_bl[NDIM], double ucon_ks[NDIM])
{
  double r, th;
  bl_coord(X, &r, &th);

  double trans[NDIM][NDIM];
  DLOOP2 trans[mu][nu] = 0.;
  DLOOP1 trans[mu][mu] = 1.;
  double gtt, gg, gthth;
  double gvec[3];
  get_coefficients(r, th, gvec);
  gtt = gvec[0];
  gg = gvec[1];
  gthth = gvec[2];

  double twoF = (1./sqrt(gg) - (gtt/gg))*gthth;
  double Delta = gthth*(gtt/gg) + a*a;

  trans[0][1] = twoF/Delta;
  trans[3][1] = a/Delta;

  DLOOP1 ucon_ks[mu] = 0.;
  DLOOP2 ucon_ks[mu] += trans[mu][nu]*ucon_bl[nu];
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
  double gtt, gg, gthth;
  double gvec[3];
  get_coefficients(r, th, gvec);
  gtt = gvec[0];
  gg = gvec[1];
  gthth = gvec[2];

  double twoF = (1./sqrt(gg) - (gtt/gg))*gthth;
  double Delta = gthth*(gtt/gg) + a*a;

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
    trans[0][1] = twoF/Delta;
    trans[3][1] = a/Delta;

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

  // Transform to MKS or MMKS
  double invtrans[NDIM][NDIM];
  set_dxdX(X, invtrans);
  invert(&invtrans[0][0], &trans[0][0]);

  DLOOP1 tmp[mu] = 0.;
  DLOOP2 {
     tmp[mu] += trans[mu][nu]*ucon[nu];
  }
  DLOOP1 ucon[mu] = tmp[mu];

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  alpha = G->lapse[CENT][j][i];
  gamma = ucon[0]*alpha;

  beta[1] = alpha*alpha*G->gcon[CENT][0][1][j][i];
  beta[2] = alpha*alpha*G->gcon[CENT][0][2][j][i];
  beta[3] = alpha*alpha*G->gcon[CENT][0][3][j][i];

  S->P[U1][k][j][i] = ucon[1] + beta[1]*gamma/alpha;
  S->P[U2][k][j][i] = ucon[2] + beta[2]*gamma/alpha;
  S->P[U3][k][j][i] = ucon[3] + beta[3]*gamma/alpha;
}




