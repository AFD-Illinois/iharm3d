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
double get_fels(struct GridGeom *G, struct FluidState *S, int i, int j, int k, int model);

void init_electrons(struct GridGeom *G, struct FluidState *S)
{
  ZLOOPALL {
    // Set electron internal energy to constant fraction of internal energy
    double uel = fel0*S->P[UU][k][j][i];

    // Initialize entropies
    S->P[KTOT][k][j][i] = (gam-1.)*S->P[UU][k][j][i]*pow(S->P[RHO][k][j][i],-gam);

    // Initialize model entropy(ies)
    for (int ip = KTOT + 1; ip < NVAR ; ip++) {
      S->P[ip][k][j][i] = (game-1.)*uel*pow(S->P[RHO][k][j][i],-game);
    }
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
  // To keep track of electron heating model
  int etracker = E_MODELS;

  // Actual entropy at final time
  double kHarm = (gam-1.)*Sf->P[UU][k][j][i]/pow(Sf->P[RHO][k][j][i],gam);

  //double uel = 1./(game-1.)*S->P[KEL][k][j][i]*pow(S->P[RHO][k][j][i],game);

  // Evolve model entropy(ies)
  for (int ip = KTOT + 1; ip < NVAR ; ip++) {
    for (int eFlagIndex = 0; eFlagIndex < sizeof(eFlagsArray) / sizeof(eFlagsArray[0]); eFlagIndex++) {
      if (etracker & eFlagsArray[eFlagIndex]) {
        double fel = get_fels(G, Ss, i, j, k, eFlagsArray[eFlagIndex]);
        Sf->P[ip][k][j][i] += fel * (game-1.)/(gam-1.) * pow(Ss->P[RHO][k][j][i],gam-game) * (kHarm - Sf->P[KTOT][k][j][i]);
        
        // update the tracker
        etracker -= eFlagsArray[eFlagIndex];
        break;
      }
    } 
  }

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
inline double get_fels(struct GridGeom *G, struct FluidState *S, int i, int j, int k, int model)
{
  get_state(G, S, i, j, k, CENT);
  double bsq = bsq_calc(S, i, j, k);
  double fel = 0.0;

if (model & CONSTANT) {
#if (E_MODELS & CONSTANT) // Special treatment for constant fel since it involves a runtime parameter
  fel = fel_constant;
#endif
}

else if (model & KAWAZURA) {
	// Equation (2) in http://www.pnas.org/lookup/doi/10.1073/pnas.1812491116
  double Tpr = (gamp-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
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
} 

else if (model & WERNER) {
	// Equation (3) in http://academic.oup.com/mnras/article/473/4/4840/4265350
  double sigma = bsq/S->P[RHO][k][j][i];
  fel = 0.25*(1+pow(((sigma/5.)/(2+(sigma/5.))), .5));
} 

else if (model & ROWAN) {
	// Equation (34) in https://iopscience.iop.org/article/10.3847/1538-4357/aa9380
  double pres = (gamp-1.)*S->P[UU][k][j][i]; // Proton pressure
  double pg = (gam-1)*S->P[UU][k][j][i];
  double beta = pres/bsq*2;
  double sigma = bsq/(S->P[RHO][k][j][i]+S->P[UU][k][j][i]+pg);
  double betamax = 0.25/sigma;
  fel = 0.5*exp(-pow(1-beta/betamax, 3.3)/(1+1.2*pow(sigma, 0.7)));
} 

else if (model & SHARMA) {
	// Equation for \delta on  pg. 719 (Section 4) in https://iopscience.iop.org/article/10.1086/520800
  double Tpr = (gamp-1.)*S->P[UU][k][j][i]/S->P[RHO][k][j][i];
  double uel = 1./(game-1.)*S->P[model][k][j][i]*pow(S->P[RHO][k][j][i],game);
  double Tel = (game-1.)*uel/S->P[RHO][k][j][i];
  if(Tel <= 0.) Tel = SMALL;
  if(Tpr <= 0.) Tpr = SMALL;
  double Trat_inv = fabs(Tel/Tpr); //Inverse of the temperature ratio in KAWAZURA
  double QeQi = 0.33 * pow(Trat_inv, 0.5);
	fel = 1./(1.+1./QeQi);
}

#if SUPPRESS_HIGHB_HEAT
  if(bsq/S->P[RHO][k][j][i] > 1.) fel = 0;
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

  double kelmax = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam-game)/(tptemin*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));
  double kelmin = S->P[KTOT][k][j][i]*pow(S->P[RHO][k][j][i],gam-game)/(tptemax*(gam-1.)/(gamp-1.) + (gam-1.)/(game-1.));

  // Replace NANs with cold electrons
  for (int ip = KTOT + 1; ip < NVAR ; ip++) {
    if (isnan(S->P[ip][k][j][i])) S->P[ip][k][j][i] = kelmin;
	// Enforce maximum Tp/Te
    S->P[ip][k][j][i] = MY_MAX(S->P[ip][k][j][i], kelmin);
	// Enforce minimum Tp/Te
    S->P[ip][k][j][i] = MY_MIN(S->P[ip][k][j][i], kelmax);
  }

}
#endif // ELECTRONS