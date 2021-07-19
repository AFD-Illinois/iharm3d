/******************************************************************************
 *                                                                            *
 * ELECTRONS.C                                                                *
 *                                                                            *
 * ELECTRON THERMODYNAMICS                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if ELECTRONS

// TODO put these in options with a default in decs.h
// Defined as in decs.h, CONSTANT not included in ALLMODELS version
#define KAWAZURA  9
#define WERNER    10
#define ROWAN     11
#define SHARMA    12
#define CONSTANT 5
#define FE_MODEL KAWAZURA

void fixup_electrons_1zone(struct FluidState *S, int i, int j, int k);
void heat_electrons_1zone(struct GridGeom *G, struct FluidState *Sh, struct FluidState *S, int i, int j, int k);
double get_fel(struct GridGeom *G, struct FluidState *S, int i, int j, int k);
double get_fels(struct GridGeom *G, struct FluidState *S, int i, int j, int k, int model);

void init_electrons(struct GridGeom *G, struct FluidState *S)
{
  ZLOOPALL {
    // Set electron internal energy to constant fraction of internal energy
    double uel = fel0*S->P[UU][k][j][i];

    // Initialize entropies
    S->P[KTOT][k][j][i] = (gam-1.)*S->P[UU][k][j][i]*pow(S->P[RHO][k][j][i],-gam);

    // Initialize different models' entropies
    #if ALLMODELS
    for (int idx = KEL0; idx < NVAR ; idx++) {
      S->P[idx][k][j][i] = (game-1.)*uel*pow(S->P[RHO][k][j][i],-game);
    }
    #else
    S->P[KEL][k][j][i] = (game-1.)*uel*pow(S->P[RHO][k][j][i],-game);
    #endif
  }

  // Necessary?  Usually called right afterward
  set_bounds(G, S);
}

// TODO merge these
void heat_electrons(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf)
{
  timer_start(TIMER_ELECTRON_HEAT);

#pragma omp parallel for collapse(3)
  ZLOOP {
    heat_electrons_1zone(G, Ss, Sf, i, j, k);
  }

  timer_stop(TIMER_ELECTRON_HEAT);
}

inline void heat_electrons_1zone(struct GridGeom *G, struct FluidState *Ss, struct FluidState *Sf, int i, int j, int k)
{
  // Actual entropy at final time
  double kHarm = (gam-1.)*Sf->P[UU][k][j][i]/pow(Sf->P[RHO][k][j][i],gam);

  //double uel = 1./(game-1.)*S->P[KEL][k][j][i]*pow(S->P[RHO][k][j][i],game);

  // Evolve different models' entropies
  #if ALLMODELS
  for (int idx = KEL0; idx < NVAR ; idx++) {
    double fel = get_fels(G, Ss, i, j, k, idx);
    Sf->P[idx][k][j][i] += (game-1.)/(gam-1.)*pow(Ss->P[RHO][k][j][i],gam-game)*fel*(kHarm - Sf->P[KTOT][k][j][i]);
  }
  #else
  double fel = get_fel(G, Ss, i, j, k);
  Sf->P[KEL][k][j][i] += (game-1.)/(gam-1.)*pow(Ss->P[RHO][k][j][i],gam-game)*fel*(kHarm - Sf->P[KTOT][k][j][i]);
  #endif

  // TODO bhlight calculates Qvisc here instead of this
  //double ugHat = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam)/(gam-1.);
  //double ugHarm = S->P[UU][k][j][i];

  // Update electron internal energy
  //uel += fel*(ugHarm - ugHat)*pow(Sh->P[RHO][k][j][i]/S->P[RHO][k][j][i],gam-game);

  // Convert back to electron entropy
  //S->P[KEL][k][j][i] = uel*(game-1.)*pow(S->P[RHO][k][j][i],-game);

  // Reset total entropy
  Sf->P[KTOT][k][j][i] = kHarm;
}

// New function for ALLMODELS runs.
#if ALLMODELS
inline double get_fels(struct GridGeom *G, struct FluidState *S, int i, int j, int k, int model)
{
  get_state(G, S, i, j, k, CENT);
  double bsq = bsq_calc(S, i, j, k);
  double fel = 0.0;
if (model == KAWAZURA) {
  double Tpr = (gam-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
  double uel = 1./(game-1.)*S->P[model][k][j][i]*pow(S->P[RHO][k][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][k][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;

  double Trat = fabs(Tpr/Tel);
  double pres = S->P[RHO][k][j][i]*Tpr; // Proton pressure
  double beta = pres/bsq*2;
  if(beta > 1.e20) beta = 1.e20;
  
  double QiQe = 35./(1. + pow(beta/15.,-1.4)*exp(-0.1/Trat));
  fel = 1./(1. + QiQe);
} else if (model == WERNER) {
  double sigma = bsq/S->P[RHO][k][j][i];
  fel = 0.25*(1+pow(((sigma/5.)/(2+(sigma/5.))), .5));
} else if (model == ROWAN) {
  double pres = (gam-1.)*S->P[UU][k][j][i]; // Proton pressure
  double beta = pres/bsq*2;
  double sigma = bsq/(S->P[RHO][k][j][i]+S->P[UU][k][j][i]+pres);
  double betamax = 0.25/sigma;
  fel = 0.5*exp(-pow(1-beta/betamax, 3.3)/(1+1.2*pow(sigma, 0.7)));
} else if (model == SHARMA) {
  double Tpr = (gam-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
  double uel = 1./(game-1.)*S->P[model][k][j][i]*pow(S->P[RHO][k][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][k][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;

  double Trat = fabs(Tel/Tpr);
  fel = 0.33 * pow(Trat, 0.5);
}

#if SUPPRESS_HIGHB_HEAT
  if(bsq/S->P[RHO][k][j][i] > 1.) fel = 0;
#endif

  return fel;
}
//////////////////////////////////
#else
//////////////////////////////////
inline double get_fel(struct GridGeom *G, struct FluidState *S, int i, int j, int k)
{
  get_state(G, S, i, j, k, CENT);
  double bsq = bsq_calc(S, i, j, k);

#if FE_MODEL == KAWAZURA
  double Tpr = (gam-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
  double uel = 1./(game-1.)*S->P[KEL][k][j][i]*pow(S->P[RHO][k][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][k][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;

  double Trat = fabs(Tpr/Tel);
  double pres = S->P[RHO][k][j][i]*Tpr; // Proton pressure
  double beta = pres/bsq*2;
  if(beta > 1.e20) beta = 1.e20;
  
  double QiQe = 35./(1. + pow(beta/15.,-1.4)*exp(-0.1/Trat));
  double fel = 1./(1. + QiQe);
#elif FE_MODEL == WERNER
  double sigma = bsq/S->P[RHO][k][j][i];
  double fel = 0.25*(1+pow(((sigma/5.)/(2+(sigma/5.))), .5));
#elif FE_MODEL == ROWAN
  double pres = (gam-1.)*S->P[UU][k][j][i]; // Proton pressure
  double beta = pres/bsq*2;
  double sigma = bsq/(S->P[RHO][k][j][i]+S->P[UU][k][j][i]+pres);
  double betamax = 0.25/sigma;
  double fel = 0.5*exp(-pow(1-beta/betamax, 3.3)/(1+1.2*pow(sigma, 0.7)));
#elif FE_MODEL == SHARMA
  double Tpr = (gam-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
  double uel = 1./(game-1.)*S->P[model][k][j][i]*pow(S->P[RHO][k][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][k][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;

  double Trat = fabs(Tel/Tpr);
  double fel = 0.33 * pow(Trat, 0.5);
#elif FE_MODEL == CONSTANT
  double fel = fel0;
#endif

#if SUPPRESS_HIGHB_HEAT
  if(bsq/S->P[RHO][k][j][i] > 1.) fel = 0;
#endif

  return fel;
}
#endif

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

  double kelmax = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam-game)/(tptemin*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));
  double kelmin = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam-game)/(tptemax*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));

  #if ALLMODELS
  for (int idx = KEL0; idx < NVAR ; idx++) {
    if (isnan(S->P[idx][k][j][i])) S->P[idx][k][j][i] = kelmin;
    S->P[idx][k][j][i] = MY_MAX(S->P[idx][k][j][i], kelmin);
    S->P[idx][k][j][i] = MY_MIN(S->P[idx][k][j][i], kelmax);
  }
  #else
  // Replace NANs with cold electrons
  if (isnan(S->P[KEL][k][j][i])) S->P[KEL][k][j][i] = kelmin;

  // Enforce maximum Tp/Te
  S->P[KEL][k][j][i] = MY_MAX(S->P[KEL][k][j][i], kelmin);

  // Enforce minimum Tp/Te
  S->P[KEL][k][j][i] = MY_MIN(S->P[KEL][k][j][i], kelmax);
  #endif
}
#endif // ELECTRONS

