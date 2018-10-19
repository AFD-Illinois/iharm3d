/******************************************************************************
 *                                                                            *
 * FIXUP.C                                                                    *
 *                                                                            *
 * REPAIR INTEGRATION FAILURES                                                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Allow to specify old fluid-frame floors for stable problems
// These impose no maximum magnetization
#ifndef FLUID_FRAME_FLOORS
#define FLUID_FRAME_FLOORS 0
#endif
// Can use drift-frame floors to stabilize high magnetization
#define DRIFT_FLOORS 0

// Point in m, around which to steepen floor prescription, eventually toward r^-3
#define FLOOR_R_CHAR 10

#if DEBUG
int nfixed = 0, nfixed_b = 0;
#endif

static struct FluidState *Stmp;

void fixup1zone(struct GridGeom *G, struct FluidState *S, int i, int j, int k);

// Apply floors to density, internal energy
// TODO can this be made faster?  I may be calling get_state too much
void fixup(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);

  static int firstc = 1;
  if (firstc) {Stmp = calloc(1,sizeof(struct FluidState)); firstc = 0;}

  ZLOOPALL fflag[k][j][i] = 0;

  if (!FLUID_FRAME_FLOORS) get_state_vec(G, S, CENT, 0, N3-1, 0, N2-1, 0, N1-1);

  // TODO rewrite this in light of fflag
#if DEBUG
#pragma omp parallel for collapse(3) reduction(+:nfixed) reduction(+:nfixed_b)
#else
#pragma omp parallel for collapse(3)
#endif
  ZLOOP fixup1zone(G, S, i, j, k);

#if DEBUG
  nfixed = mpi_reduce_int(nfixed);
  nfixed_b = mpi_reduce_int(nfixed_b);

  LOGN("Fixed: %d", nfixed);
  LOGN("From field: %d", nfixed_b);
  LOGN("Proportion: %f", ((double) nfixed_b) / nfixed);

  nfixed = 0;
  nfixed_b = 0;
#endif

  LOG("End fixup");

  timer_stop(TIMER_FIXUP);
}

inline void ucon_to_utcon(struct GridGeom *G, int i, int j, double ucon[NDIM], double utcon[NDIM])
{
  double beta[NDIM];

  double alpha = 1./sqrt(-G->gcon[CENT][0][0][j][i]);
  for (int mu = 1; mu < NDIM; mu++) {
    beta[mu] = G->gcon[CENT][0][mu][j][i]*alpha*alpha;
  }
  double gamma = alpha*ucon[0];

  utcon[0] = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    utcon[mu] = ucon[mu] + gamma*beta[mu]/alpha;
  }
}

inline double ut_calc_3vel(struct GridGeom *G, int i, int j, double vcon[NDIM])
{
  double AA, BB, CC, DD;

  AA = G->gcov[CENT][0][0][j][i];
  BB = 2.*(G->gcov[CENT][0][1][j][i]*vcon[1] +
      G->gcov[CENT][0][2][j][i]*vcon[2] +
      G->gcov[CENT][0][3][j][i]*vcon[3]);
  CC = G->gcov[CENT][1][1][j][i]*vcon[1]*vcon[1] +
      G->gcov[CENT][2][2][j][i]*vcon[2]*vcon[2] +
      G->gcov[CENT][3][3][j][i]*vcon[3]*vcon[3] +
       2.*(G->gcov[CENT][1][2][j][i]*vcon[1]*vcon[2] +
	   G->gcov[CENT][1][3][j][i]*vcon[1]*vcon[3] +
	   G->gcov[CENT][2][3][j][i]*vcon[2]*vcon[3]);
  DD = 1./(AA + BB + CC);

  if (DD < -G->gcon[CENT][0][0][j][i]) {
    DD = -G->gcon[CENT][0][0][j][i];
  }

  return sqrt(DD);
}

inline void fixup1zone(struct GridGeom *G, struct FluidState *S, int i, int j, int k)
{
  double rhoflr, uflr;
  double rhoscal, uscal;

  if(METRIC == MKS) {
    double r, th, X[NDIM];
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    // New, steeper floor in rho
    rhoscal = pow(r, -2.) * 1 / (1 + r/FLOOR_R_CHAR);
    // Classic harm floors
//    rhoscal = pow(r, -2.);
    uscal = pow(rhoscal, gam);
    // Very classic floors
//    rhoscal = pow(r,-1.5);
//    uscal = rhoscal/r;
    rhoflr = RHOMIN*rhoscal;
    uflr = UUMIN*uscal;
  } else if (METRIC == MINKOWSKI) {
    rhoscal = 1.e-2;
    uscal = 1.e-2;
    rhoflr = RHOMIN*rhoscal;
    uflr = UUMIN*uscal;
  }

  rhoflr = MY_MAX(rhoflr, RHOMINLIMIT);
  uflr = MY_MAX(uflr, UUMINLIMIT);

#if FLUID_FRAME_FLOORS
  // Fluid frame floors: don't conserve momentum
  if (S->P[RHO][k][j][i] < rhoflr) {
    fflag[k][j][i] = HIT_FLOOR_GEOM;
    S->P[RHO][k][j][i] = rhoflr;
  }
  if (S->P[UU][k][j][i] < uflr) {
    fflag[k][j][i] = HIT_FLOOR_GEOM;
    S->P[UU][k][j][i] = uflr;
  }
#else
  // Do this call if get_state_vec not called above
  // TODO thinner bsq calculation w/o full call -- also check whether obs frame floors should be this high
  //get_state(G, S, i, j, k, CENT);
  double bsq = bsq_calc(S, i, j, k);
  double rhoflr_b = bsq/BSQORHOMAX;
  double uflr_b = bsq/BSQOUMAX;
  // TODO apply this after updates to UU (or play equivalent tricks)
  rhoflr = MY_MAX(rhoflr, S->P[UU][k][j][i]/UORHOMAX);

  if (rhoflr > S->P[RHO][k][j][i] || uflr > S->P[UU][k][j][i] ||
      rhoflr_b > S->P[RHO][k][j][i] || uflr_b > S->P[UU][k][j][i]) { // Apply floors
#if DEBUG
    nfixed++;
#endif

    // If we would not have hit the floors w/o magnetic floors, record this
    if(S->P[RHO][k][j][i] > rhoflr && S->P[UU][k][j][i] > uflr) {
      fflag[k][j][i] = HIT_FLOOR_SIGMA;

#if DEBUG
      nfixed_b++;
#endif
    } else {
      fflag[k][j][i] = HIT_FLOOR_GEOM;
    }

    // Set single consistent floor
    rhoflr = MY_MAX(rhoflr, rhoflr_b);
    uuflr = MY_MAX(uflr, uflr_b);

#if DRIFT_FLOORS
    double trans = 10.*bsq/MY_MIN(S->P[RHO][k][j][i], S->P[UU][k][j][i]) - 1.;
    if (trans > 0.) { // Strongly magnetized region; use drift frame floors
      // Preserve pre-floor primitives
      double pv_prefloor[NVAR];
      PLOOP pv_prefloor[ip] = S->P[ip][k][j][i];

      // Update according to floor
      S->P[RHO][k][j][i] = MY_MAX(S->P[RHO][k][j][i], rhoflr);
      S->P[UU][k][j][i] = MY_MAX(S->P[UU][k][j][i], uflr);

      trans = MY_MIN(trans, 1.);

      // Set velocity to drift velocity
      double betapar = -S->bcon[0][k][j][i]/((bsq + SMALL)*S->ucon[0][k][j][i]);
      double betasq = MY_MIN((betapar*betapar*bsq),
			     (1. - 1./(GAMMAMAX*GAMMAMAX)));

      double gamma = 1./sqrt(1. - betasq);

      double ucon[NDIM], ucon_dr[NDIM]; // TODO ucon name is reused for ut*vcon below
      DLOOP1 {
        ucon[mu] = S->ucon[mu][k][j][i];
        ucon_dr[mu] = gamma
            * (S->ucon[mu][k][j][i] + betapar * S->bcon[mu][k][j][i]);
      }

      double bcon_local[NDIM], bcov_local[NDIM];
      DLOOP1 {
        bcon_local[mu] = (mu == 0) ? 0 : S->P[B1 - 1 + mu][k][j][i];
      }

      double gcov_local[NDIM][NDIM];
      get_gcov(G, i, j, CENT, gcov_local);
      lower(bcon_local, gcov_local, bcov_local);

      double bsq_local = dot(bcon_local, bcov_local);
      double b_local = sqrt(bsq_local);

      double udotB = dot(ucon, bcov_local);

      // Enthalpy before floors are applied
      double wold = pv_prefloor[RHO] + pv_prefloor[UU]*gam;
      double QdotB = udotB*wold*ucon[0];

      // Apply floors to enthalpy and recompute parallel velocity
      double wnew = S->P[RHO][k][j][i] + S->P[UU][k][j][i]*gam;
      double x = 2.*QdotB/(b_local*wnew*ucon_dr[0] + SMALL);
      double vpar = x/(ucon_dr[0]*(1. + sqrt(1. + x*x)));

      double vcon[NDIM];
      vcon[0]  = 1.;
      for (int mu = 1; mu < NDIM; mu++) {
        vcon[mu] = vpar*bcon_local[mu]/(b_local + SMALL) + ucon_dr[mu]/ucon_dr[0];
      }

      double ut = 0, utcon[NDIM];
      ut = ut_calc_3vel(G, i, j, vcon);

      DLOOP1 ucon[mu] = ut*vcon[mu];
      ucon_to_utcon(G, i, j, ucon, utcon);

      // Convert 3-velocity to relative 4-velocity and store in primitives
      for (int mu = 1; mu < NDIM; mu++) {
        S->P[mu+UU][k][j][i] = utcon[mu]*trans + pv_prefloor[mu+UU]*(1. - trans);
      }
    } else { // Weakly magnetized region; use normal observer frame floors
#endif

      PLOOP {Stmp->P[ip][k][j][i] = 0; Stmp->U[ip][k][j][i] = 0;}
      Stmp->P[RHO][k][j][i] = MY_MAX(0., rhoflr - S->P[RHO][k][j][i]);
      Stmp->P[UU][k][j][i] = MY_MAX(0., uflr - S->P[UU][k][j][i]);

      get_state(G, Stmp, i, j, k, CENT);
      prim_to_flux(G, Stmp, i, j, k, 0, CENT, Stmp->U);

      get_state(G, S, i, j, k, CENT);
      prim_to_flux(G, S, i, j, k, 0, CENT, S->U);

      PLOOP {
        S->U[ip][k][j][i] += Stmp->U[ip][k][j][i];
        S->P[ip][k][j][i] += Stmp->P[ip][k][j][i];
      }

      U_to_P(G, S, i, j, k, CENT);
      get_state(G, S, i, j, k, CENT); // TODO needed?

#if DRIFT_FLOORS
    } // if (trans > 0)
#endif
  }
#endif

#if ELECTRONS
  // Reset entropy after floors
  S->P[KTOT][k][j][i] = (gam - 1.)*S->P[UU][k][j][i]/pow(S->P[RHO][k][j][i],gam);

  // Keep to KTOTMAX by controlling u, to avoid anomalous cooling from funnel wall
  if (S->P[KTOT][k][j][i] > KTOTMAX) {
    fflag[k][j][i] = HIT_FLOOR_KTOT;
    S->P[UU][k][j][i] = KTOTMAX*pow(S->P[RHO][k][j][i],gam)/(gam-1.);
    S->P[KTOT][k][j][i] = KTOTMAX;
  }
#endif // ELECTRONS

  // Limit gamma with respect to normal observer
  // TODO check for fail here
  double gamma = mhd_gamma_calc(G, S, i, j, k, CENT);

  if (gamma > GAMMAMAX) {
    fflag[k][j][i] = HIT_FLOOR_GAMMA;

    double f = sqrt((GAMMAMAX*GAMMAMAX - 1.)/(gamma*gamma - 1.));
    S->P[U1][k][j][i] *= f;
    S->P[U2][k][j][i] *= f;
    S->P[U3][k][j][i] *= f;
    get_state(G, S, i, j, k, CENT);
  }

}

// Replace bad points with values interpolated from neighbors
#define FLOOP for(int ip=0;ip<B1;ip++)
void fixup_utoprim(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);

  // Flip the logic of the pflag[] so that it now indicates which cells are good
#pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    pflag[k][j][i] = !pflag[k][j][i];
  }

#if DEBUG
  int nbad_utop = 0;
#pragma omp parallel for simd collapse(2) reduction (+:nbad_utop)
  ZLOOP {
    // Count the 0 = bad cells
    nbad_utop += !pflag[k][j][i];
  }
  LOGN("Fixing %d bad cells", nbad_utop);
#endif

  // Make sure we are not using ill defined physical corner regions
  // TODO find a way to do this once, or put it in bounds at least?
  for (int k = 0; k < NG; k++) {
    for (int j = 0; j < NG; j++) {
      for (int i = 0; i < NG; i++) {
        if(global_start[2] == 0 && global_start[1] == 0 && global_start[0] == 0) pflag[k][j][i] = 0;
        if(global_start[2] == 0 && global_start[1] == 0 && global_stop[0] == N1TOT) pflag[k][j][i+N1+NG] = 0;
        if(global_start[2] == 0 && global_stop[1] == N2TOT && global_start[0] == 0) pflag[k][j+N2+NG][i] = 0;
        if(global_stop[2] == N3TOT && global_start[1] == 0 && global_start[0] == 0) pflag[k+N3+NG][j][i] = 0;
        if(global_start[2] == 0 && global_stop[1] == N2TOT && global_stop[0] == N1TOT) pflag[k][j+N2+NG][i+N1+NG] = 0;
        if(global_stop[2] == N3TOT && global_start[1] == 0 && global_stop[0] == N1TOT) pflag[k+N3+NG][j][i+N1+NG] = 0;
        if(global_stop[2] == N3TOT && global_stop[1] == N2TOT && global_start[0] == 0) pflag[k+N3+NG][j+N2+NG][i] = 0;
        if(global_stop[2] == N3TOT && global_stop[1] == N2TOT && global_stop[0] == N1TOT) pflag[k+N3+NG][j+N2+NG][i+N1+NG] = 0;
      }
    }
}

#if DEBUG
  // Keep track of how many points we fix
  int nfixed_utop = 0;
#endif

  // TODO is parallelizing this version okay?
//#pragma omp parallel for collapse(3) reduction(+:bad)
  ZLOOP {
    if (pflag[k][j][i] == 0) {
      double wsum = 0.;
      double sum[B1];
      FLOOP sum[ip] = 0.;
      for (int l = -1; l < 2; l++) {
        for (int m = -1; m < 2; m++) {
          for (int n = -1; n < 2; n++) {
            double w = 1./(abs(l) + abs(m) + abs(n) + 1)*pflag[k+n][j+m][i+l];
            wsum += w;
            FLOOP sum[ip] += w*S->P[ip][k+n][j+m][i+l];
          }
        }
      }
      if(wsum < 1.e-10) {
        fprintf(stderr, "fixup_utoprim problem: No usable neighbors!\n");
        // TODO set to something okay here.  This happens /very rarely/
        continue;
      }
      FLOOP S->P[ip][k][j][i] = sum[ip]/wsum;

      // Cell is fixed, could now use for other interpolations
      // However, this is harmful to MPI-determinism
      //pflag[k][j][i] = 1;

#if DEBUG
      nfixed_utop++;
#endif

      if (!FLUID_FRAME_FLOORS) get_state(G, S, i, j, k, CENT);
      fixup1zone(G, S, i, j, k);
      // TODO need final call?
      if (!FLUID_FRAME_FLOORS) get_state(G, S, i, j, k, CENT);
    }
  }

#if DEBUG
  int nleft_utop = nbad_utop - nfixed_utop;
  if(nleft_utop > 0) fprintf(stderr,"Cells STILL BAD after fixup_utoprim: %d\n", nleft_utop);
#endif

  // Re-initialize the pflag
#pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    pflag[k][j][i] = 0;
  }

  timer_stop(TIMER_FIXUP);
}
#undef FLOOP
