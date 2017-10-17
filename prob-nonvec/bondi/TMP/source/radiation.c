/******************************************************************************
 *                                                                            *
 * RADIATION.C                                                                *
 *                                                                            *
 * MODEL-INDEPENDENT RADIATION QUANTITIES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
double Bnu_inv(double nu, double Thetae)
{
  double x;
  x = HPL*nu/(ME*CL*CL*Thetae);

  if(x < 1.e-3) // Taylor-expand small arguments for numerical accuracy
    return ((2.*HPL/(CL*CL))/( x/24. * (24. + x*(12. + x*(4. + x)))));
  else
    return ((2.*HPL/(CL*CL))/(exp(x) - 1.));
}

double jnu_inv(double nu, double Thetae, double Ne, double B, double theta)
{
  double j;

  j = jnu(nu,Ne,Thetae,B,theta);

  return (j/(nu*nu));
}

// Invariant scattering opacity
double alpha_inv_scatt(double nu, double Thetae, double Ne)
{
  //double Eg = HPL*nu/(ME*CL*CL);
  double kappa = 0.;//total_scatt_cross(Eg, Thetae)/MP;

  return nu*kappa*Ne*MP;
}

// Invariant absorption opacity
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
  double theta)
{
  // Otherwise synchrotron emissivity will fail in this region
  if (isnan(jnu_inv(nu, Thetae, Ne, B, theta)))
    return 0.;
  else
    return(jnu_inv(nu, Thetae, Ne, B, theta)/(Bnu_inv(nu,Thetae) + 1.e-100));
}

// Get frequency in fluid frame, in Hz
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM])
{
  double ener, nu;

  // Energy in electron rest-mass units
  ener = -( K[0] * Ucov[0] +
            K[1] * Ucov[1] +
            K[2] * Ucov[2] +
            K[3] * Ucov[3] );

  nu = ener*ME*CL*CL/HPL;

  return nu;
}

// Return angle between magnetic field and wavevector
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
  double Bcov[NDIM], double B)
{
  double k,mu;

  if(B == 0.)
    return(M_PI/2.);

  k = fabs(K[0]*Ucov[0] + K[1]*Ucov[1] + K[2]*Ucov[2] + K[3]*Ucov[3]);

  mu = (K[0]*Bcov[0] + K[1]*Bcov[1] + K[2]*Bcov[2] + K[3]*Bcov[3])/fabs(k*B);

  #define SMALLEPS (1.e-10)

  if (fabs(mu) > 1.)
  {
    mu /= fabs(mu);

    // We need |mu| < 1 or else acos(mu) will return nan.
    if ( mu > 0. )
      mu -= SMALLEPS;
    else
      mu += SMALLEPS;
  }

  return (acos(mu));
}
#endif // RADIATION

