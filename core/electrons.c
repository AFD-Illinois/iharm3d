/******************************************************************************
 *                                                                            *
 * ELECTRONS.C                                                                *
 *                                                                            *
 * ELECTRON THERMODYNAMICS                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if ELECTRONS

void fixup_electrons_1zone(struct FluidState *S, int i, int j, int k);
void heat_electrons_1zone(struct GridGeom *G, struct FluidState *Sh, struct FluidState *S, int i, int j, int k);
double get_fel(struct GridGeom *G, struct FluidState *S, int i, int j, int k);

void init_electrons(struct FluidState *S)
{
  ZLOOPALL {
    // Set electron internal energy to constant fraction of internal energy
    double uel = fel0*S->P[UU][k][j][i];

    // Initialize entropies
    S->P[KTOT][k][j][i] =  (gam-1.)*S->P[UU][k][j][i]*pow(S->P[RHO][k][j][i],-gam);
    S->P[KEL][k][j][i] = (game-1.)*uel*pow(S->P[RHO][k][j][i],-game);
  }

  // Necessary?  Usually called right afterward
  //set_bounds(S);
}

void heat_electrons(struct GridGeom *G, struct FluidState *Sh, struct FluidState *S)
{
  timer_start(TIMER_ELECTRON_HEAT);

#pragma omp parallel for collapse(3)
  ZLOOP {
    heat_electrons_1zone(G, Sh, S, i, j, k);
  }

  timer_stop(TIMER_ELECTRON_HEAT);
}

inline void heat_electrons_1zone(struct GridGeom *G, struct FluidState *Sh, struct FluidState *S, int i, int j, int k)
{
  // Actual entropy at final time
  double kHarm = (gam-1.)*S->P[UU][k][j][i]/pow(S->P[RHO][k][j][i],gam);

  double uel = 1./(game-1.)*S->P[KEL][k][j][i]*pow(S->P[RHO][k][j][i],game);

#if BETA_HEAT
  double fel = get_fel(G, S, i, j, k);
#else
  double fel = fel0;
#endif

  double ugHat = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam)/(gam-1.);
  double ugHarm = S->P[UU][k][j][i];

  // Update electron internal energy
  uel += fel*(ugHarm - ugHat)*pow(Sh->P[RHO][k][j][i]/S->P[RHO][k][j][i],gam-game);

  // Convert back to electron entropy
  S->P[KEL][k][j][i] = uel*(game-1.)*pow(S->P[RHO][k][j][i],-game);

  // Reset total entropy
  S->P[KTOT][k][j][i] = kHarm;
}

inline double get_fel(struct GridGeom *G, struct FluidState *S, int i, int j, int k)
{
  double c1 = 0.92;

  double Tpr = (gam-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
  double uel = 1./(game-1.)*S->P[KEL][k][j][i]*pow(S->P[RHO][k][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][k][j][i];
  Tel = (Tel <= 0.) ? Tel : SMALL;
  Tpr = (Tpr <= 0.) ? Tpr : SMALL;

  double Trat = fabs(Tpr/Tel);

  double pres = S->P[RHO][k][j][i]*Tpr; // Proton pressure

  //TODO can I prevent this call?
  get_state(G, S, i, j, k, CENT);
  double bsq = bsq_calc(S, i, j, k);

  double beta = MY_MIN(pres/bsq*2,1.e20);

  double mbeta = 2. - 0.2*log10(Trat);

  double c2 = (Trat <= 1.) ? 1.6/Trat : 1.2/Trat;
  double c3 = (Trat <= 1.) ? 18. + 5.*log10(Trat) : 18.;
  double c22 = pow(c2, 2.);
  double c32 = pow(c3, 2.);

  double qrat = c1*(c22+pow(beta,mbeta))/(c32 + pow(beta,mbeta))*exp(-1./beta)*pow(MP/ME*Trat,.5);
  double fel = 1./(1. + qrat);

#if SUPPRESS_HIGHB_HEAT
  fel = (bsq/S->P[RHO][k][j][i] > 1.) ? 0. : fel;
#endif

  return fel;
}

void fixup_electrons(struct FluidState *S)
{
  timer_start(TIMER_ELECTRON_FIXUP);

#pragma omp parallel for collapse(3)
  ZLOOP {
    fixup_electrons_1zone(S, i, j, k);
  }

  timer_stop(TIMER_ELECTRON_FIXUP);
}

inline void fixup_electrons_1zone(struct FluidState *S, int i, int j, int k)
{
  double kelmax = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam-game)/(tptemin + (gam-1.)/(game-1.));
  double kelmin = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam-game)/(tptemax + (gam-1.)/(game-1.));

  // Replace NANs with cold electrons
  if (isnan(S->P[KEL][k][j][i])) S->P[KEL][k][j][i] = kelmin;

  // Enforce maximum Tp/Te
  S->P[KEL][k][j][i] = MY_MAX(S->P[KEL][k][j][i], kelmin);

  // Enforce minimum Tp/Te
  S->P[KEL][k][j][i] = MY_MIN(S->P[KEL][k][j][i], kelmax);
}
#endif // ELECTRONS

