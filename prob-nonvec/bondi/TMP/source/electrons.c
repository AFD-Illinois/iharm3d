/******************************************************************************
 *                                                                            *
 * ELECTRONS.C                                                                *
 *                                                                            *
 * ELECTRON THERMODYNAMICS                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if ELECTRONS

void heat_electrons_zone(int i, int j, int k, double Ph[NVAR], double P[NVAR],  
  double Dt);
void fixup_electrons_1zone(double P[NVAR]);

void init_electrons()
{
  double uel;
  
  ZSLOOP(-NG, N1 + NG - 1, -NG, NG + N2 - 1, -NG, NG + N3 - 1) {
    // Set electron internal energy to constant fraction of internal energy
    uel = fel0*P[i][j][k][UU];

    // Initialize entropies
    P[i][j][k][KTOT] =  (gam-1.)*P[i][j][k][UU]*pow(P[i][j][k][RHO],-gam);
    P[i][j][k][KEL] = (game-1.)*uel*pow(P[i][j][k][RHO],-game);
  }

  bound_prim(P);
}

void heat_electrons(grid_prim_type Ph, grid_prim_type P, double Dt)
{
  timer_start(TIMER_ELECTRON_HEAT);
  #pragma omp parallel for collapse(3)
  ZLOOP {
    heat_electrons_zone(i, j, k, Ph[i][j][k], P[i][j][k], Dt);
  }
  timer_stop(TIMER_ELECTRON_HEAT);
}

void heat_electrons_zone(int i, int j, int k, double Ph[NVAR], double P[NVAR],
  double Dt)
{
	double kHarm, ugHat, ugHarm, uel, fel;

  // Actual entropy at final time
  kHarm = (gam-1.)*P[UU]/pow(P[RHO],gam);

  uel = 1./(game-1.)*P[KEL]*pow(P[RHO],game);

  #if BETA_HEAT
  fel = get_fel(i, j, k, Ph);
  #else
  fel = fel0;
  #endif

  ugHat = P[KTOT]*pow(P[RHO],gam)/(gam-1.);
  ugHarm = P[UU];

  // Update electron internal energy
  uel += fel*(ugHarm - ugHat)*pow(Ph[RHO]/P[RHO],gam-game);

  // Convert back to electron entropy
  P[KEL] = uel*(game-1.)*pow(P[RHO],-game);

  // Reset total entropy
  P[KTOT] = kHarm;
}

double get_fel(int i, int j, int k, double P[NVAR])
{
  struct of_geom geom;
  double beta, fel, c1, Trat, Tpr, mbeta, qrat;
  double pres, bsq;
  double c2, c3, c22, c32;

  c1 = 0.92;

  Tpr = (gam-1.)*P[UU]/P[RHO];
  double uel = 1./(game-1.)*P[KEL]*pow(P[RHO],game);
  double Tel = (game-1.)*uel/P[RHO];
  if (Tel <= 0.) Tel = SMALL;
  if (Tpr <= 0.) Tpr = SMALL;
  Trat = fabs(Tpr/Tel);

  pres = P[RHO]*Tpr; // Proton pressure

  //get_geometry(i, j, CENT, &geom);
  geom = *get_geometry(i, j, k);
  bsq = bsq_calc(P, &geom);

  beta = pres/bsq*2.;
  if (beta > 1.e20) beta = 1.e20;
  mbeta = 2. - 0.2*log10(Trat);

  if (Trat <= 1.) {
    c2 = 1.6/Trat;
    c3 = 18. + 5.*log10(Trat);
  } else {
    c2 = 1.2/Trat;
    c3 = 18.;
  }
  c22 = pow(c2, 2.);
  c32 = pow(c3, 2.);

  qrat = c1*(c22+pow(beta,mbeta))/(c32 + pow(beta,mbeta))*exp(-1./beta)*
         pow(MP/ME*Trat,.5);
  fel = 1./(1. + qrat);

  #if SUPPRESS_HIGHB_HEAT
  if (bsq/P[RHO] > 1.) fel = 0.;
  #endif

  return fel;
}

void fixup_electrons(grid_prim_type P)
{
  timer_start(TIMER_ELECTRON_FIXUP);
  #pragma omp parallel for collapse(3)
  ZLOOP {
    fixup_electrons_1zone(P[i][j][k]);
  }
  timer_stop(TIMER_ELECTRON_FIXUP);
}

void fixup_electrons_1zone(double P[NVAR])
{
  double kelmax = P[KTOT]*pow(P[RHO],gam-game)/(TPTEMIN + (gam-1.)/(game-1.));
  double kelmin = P[KTOT]*pow(P[RHO],gam-game)/(TPTEMAX + (gam-1.)/(game-1.));

  // Replace NANs with cold electrons
  if (isnan(P[KEL])) P[KEL] = kelmin;

  // Enforce maximum Tp/Te
  P[KEL] = MY_MAX(P[KEL], kelmin);

  // Enforce minimum Tp/Te
  P[KEL] = MY_MIN(P[KEL], kelmax);
}
#endif // ELECTRONS

