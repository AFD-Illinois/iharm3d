/******************************************************************************
 *                                                                            *
 * JNU.C                                                                      *
 *                                                                            *
 * FORMULAS FOR SPECIFIC AND ANGLE-INTEGRATED EMISSIVITIES                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
// jnu: dE/(dt dnu dV dOmega)
// Jnu: dE/(dt dnu dV)

double jnu_gray(double nu, double Ne, double Thetae);
double jnu_brem(double nu, double Ne, double Thetae);
double jnu_synch(double nu, double Ne, double Thetae, double B, double theta);
double Jnu_gray(double nu, double Ne, double Thetae);
double Jnu_brem(double nu, double Ne, double Thetae);
double Jnu_synch(double nu, double Ne, double Thetae, double B);

void init_emiss_tables();

void init_emissivity() {
  #if SYNCHROTRON
  init_emiss_tables();
  #endif
}

double jnu(double nu, double Ne, double Thetae, double B, double theta)
{
  double jnu = 0.;

  #if GRAYABSORPTION
    jnu += jnu_gray(nu, Ne, Thetae);
  #endif
  #if BREMSSTRAHLUNG
    jnu += jnu_brem(nu, Ne, Thetae);
  #endif
  #if SYNCHROTRON
    jnu += jnu_synch(nu, Ne, Thetae, B, theta);
  #endif

  return jnu;
}

double Jnu(double nu, double Ne, double Thetae, double B)
{
  double Jnu = 0.;

  #if GRAYABSORPTION
    Jnu += Jnu_gray(nu, Ne, Thetae);
  #endif
  #if BREMSSTRAHLUNG
    Jnu += Jnu_brem(nu, Ne, Thetae);
  #endif
  #if SYNCHROTRON
    Jnu += Jnu_synch(nu, Ne, Thetae, B);
  #endif

  return Jnu;
}

double jnu_brem(double nu, double Ne, double Thetae)
{
  double Te = Thetae*ME*CL*CL/KBOL;
  double x,efac;

  x = HPL*nu/(KBOL*Te);

  if(x < 1.e-3) {
    efac = (24 - 24*x + 12*x*x - 4.*x*x*x + x*x*x*x)/24.;
  } else {
    efac = exp(-x);
  }

  return 5.4e-39*Ne*Ne*1./sqrt(Te)*efac;
}

double Jnu_brem(double nu, double Ne, double Thetae)
{
  return 4.*M_PI*jnu_brem(nu, Ne, Thetae);
}

double jnu_integrand_synch(double th, void *params);
double F_eval_synch(double Thetae, double Bmag, double nu);
double F_synch[NU_BINS+1], K2[NU_BINS+1];
static double lK_min, dlK;
static double lT_min, dl_T;
#define THETAE_MIN  0.3 // Synchrotron fitting function invalid for low Thetae
#define EPSABS 0.
#define EPSREL 1.e-6
#define KMIN (0.002)
#define KMAX (1.e7)
#define TMIN (THETAE_MIN)
#define TMAX (1.e2)
void init_emiss_tables()
{
  int k;
  double result,err,K,T;
  gsl_function func;
  gsl_integration_workspace *w;

  func.function = &jnu_integrand_synch;
  func.params = &K;

  lK_min = log(KMIN);
  dlK = log(KMAX/KMIN)/(NU_BINS);

  lT_min = log(TMIN);
  dl_T = log(TMAX/TMIN)/(NU_BINS);

  //  build table for F(K) where F(K) is given by
  // \int_0^\pi ( (K/\sin\theta)^{1/2} + 2^{11/12}(K/\sin\theta)^{1/6})^2 \exp[-(K/\sin\theta)^{1/3}]
  // so that J_{\nu} = const.*F(K)

  w = gsl_integration_workspace_alloc(1000);
  for(k = 0; k <= NU_BINS; k++){
    K = exp(k*dlK + lK_min);
    gsl_integration_qag(&func, 0., M_PI/2., EPSABS, EPSREL, 1000,
      GSL_INTEG_GAUSS61, w, &result, &err);
    F_synch[k] = log(4*M_PI*result);
  }
  gsl_integration_workspace_free(w);

  // Build table for quick evaluation of the bessel function K2 for emissivity
  for(k = 0; k <= NU_BINS; k++){
    T = exp(k*dl_T + lT_min);
    K2[k] = log(gsl_sf_bessel_Kn(2,1./T));
  }
}

double linear_interp_K2_synch(double Thetae)
{
  int i;
  double di,lT;

  lT = log(Thetae);

  di = (lT - lT_min)/dl_T;
  i = (int)di;
  di = di - i;

  return exp((1.-di)*K2[i] + di*K2[i+1]);
}

double linear_interp_F_synch(double K)
{
  if (K < KMIN || K > KMAX) return 0.;

  int i;
  double di,lK;

  lK = log(K);

  di = (lK - lK_min)/dlK;
  i = (int)di;
  di = di-i;

  return exp((1.-di)*F_synch[i] + di*F_synch[i+1]);
}

double K2_eval_synch(double Thetae)
{
  if(Thetae < THETAE_MIN) return 0.;
  if(Thetae > TMAX) return 2.*Thetae*Thetae;

  return linear_interp_K2_synch(Thetae);
}

#define CST 1.88774862536 // 2^{11/12}
double jnu_synch(double nu, double Ne, double Thetae, double B, double theta)
{
  double K2,nuc,nus,x,f,j,sth,xp1, xx ;

  if(Thetae < THETAE_MIN) return 0.;

  K2 = K2_eval_synch(Thetae);

  nuc = EE*B/(2.*M_PI*ME*CL) ;
  sth = sin(theta) ;
  nus = (2./9.)*nuc*Thetae*Thetae*sth;
  if(nu > 1.e12*nus) return(0.) ;
  x = nu/nus ;
  xp1 = pow(x,1./3.);
  xx = sqrt(x) + CST*sqrt(xp1);
  f = xx*xx;
  j = (M_SQRT2*M_PI*EE*EE*Ne*nus/(3.*CL*K2)) * f * exp(-xp1) ;

  return (j);
}
#undef CST

double Jnu_synch(double nu, double Ne, double Thetae, double Bmag)
{
  double nuC = EE*Bmag/(2*M_PI*ME*CL);

  double F_interp  = F_eval_synch(Thetae, Bmag, nu);
  double K2_interp = K2_eval_synch(Thetae);

  double Jnu = sqrt(2.)*M_PI*EE*EE*Ne*(2./9.*nuC)*Thetae*Thetae;
  Jnu       /= (3.*CL*K2_interp);
  Jnu       *= F_interp;

  if (isnan(Jnu) || isinf(Jnu))
    return 0.;

  if (Thetae < THETAE_MIN)
    return 0.;

  return Jnu;
}

#define CST 1.88774862536 /* 2^{11/12} */
double jnu_integrand_synch(double th, void *params)
{
  double K = *(double *)params;
  double sth = sin(th);
  double x = K/sth;

  if(sth < 1.e-150 || x > 2.e8) return 0.;

  return sth*sth*pow(sqrt(x) + CST*pow(x,1./6.),2.)*exp(-pow(x,1./3.));
}
#undef CST

#define KFAC  (9*M_PI*ME*CL/EE)
double F_eval_synch(double Thetae, double Bmag, double nu)
{
  double K, x;
  double linear_interp_F_synch(double);

  K = KFAC*nu/(Bmag*Thetae*Thetae);

  if(K > KMAX) {
    return 0.;
  } else if(K < KMIN) {
    /* use a good approximation */
    x = pow(K, 0.333333333333333333);
    return (x*(37.67503800178 + 2.240274341836*x));
  } else if (isnan(K)) {
    return 0;
  } else {
    return linear_interp_F_synch(K);
  }
}
#undef KFAC
#undef KMIN
#undef KMAX
#undef EPSABS
#undef EPSREL
#undef THETAE_MIN
#undef TMIN
#undef TMAX

double jnu_gray(double nu, double Ne, double Thetae)
{
  double opacity = MP*Ne*kappa;

  return opacity*Bnu_inv(nu, Thetae)*nu*nu*nu;
}

double Jnu_gray(double nu, double Ne, double Thetae)
{
  return 4.*M_PI*jnu_gray(nu, Ne, Thetae);
}
#endif // RADIATION
