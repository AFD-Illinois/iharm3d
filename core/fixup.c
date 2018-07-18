/******************************************************************************
 *                                                                            *
 * FIXUP.C                                                                    *
 *                                                                            *
 * REPAIR INTEGRATION FAILURES                                                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#define OLD_FLOORS 0
#define DRIFT_FLOORS 0

// Apply floors to density, internal energy
// TODO can this be made faster?  I may be calling get_state too much
void fixup(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);

#if !OLD_FLOORS
  get_state_vec(G, S, CENT, 0, N3-1, 0, N2-1, 0, N1-1);
#endif

#pragma omp parallel for collapse(3)
  ZLOOP fixup1zone(G, S, i, j, k);

  timer_stop(TIMER_FIXUP);
}

inline void ucon_to_utcon(struct GridGeom *G, int i, int j, double ucon[NDIM], double utcon[NDIM])
{
  double alpha, beta[NDIM], gamma;

  alpha = 1./sqrt(-G->gcon[CENT][0][0][j][i]);
  for (int mu = 1; mu < NDIM; mu++) {
    beta[mu] = G->gcon[CENT][0][mu][j][i]*alpha*alpha;
  }
  gamma = alpha*ucon[0];

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

#define BSQORHOMAX (50.)
#define BSQOUMAX (2500.)
#define UORHOMAX (50.)
inline void fixup1zone(struct GridGeom *G, struct FluidState *S, int i, int j, int k)
{
  double rhoflr, uflr, f;
  double rhoscal, uscal;

  if(METRIC == MKS) {
    double r, th, X[NDIM];
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    // Classic harm floors
    rhoscal = pow(r, -2.);
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

#if OLD_FLOORS
  // Fluid frame floors: don't conserve momentum
  if (S->P[RHO][k][j][i] < rhoflr) {
    S->P[RHO][k][j][i] = rhoflr;
  }
  if (S->P[UU][k][j][i] < uflr) {
    S->P[UU][k][j][i] = uflr;
  }
#else
  // Do this call if get_state_vec not called above
  // TODO thinner bsq calculation w/o full call
  //get_state(G, S, i, j, k, CENT);
  double bsq = bsq_calc(S, i, j, k);
  rhoflr = MY_MAX(rhoflr, bsq/BSQORHOMAX);
  uflr = MY_MAX(uflr, bsq/BSQOUMAX);
  rhoflr = MY_MAX(rhoflr, S->P[UU][k][j][i]/UORHOMAX);

  if (rhoflr > S->P[RHO][k][j][i] || uflr > S->P[UU][k][j][i]) { // Apply floors

    double trans = 10.*bsq/MY_MIN(S->P[RHO][k][j][i], S->P[UU][k][j][i]) - 1.;

    if (trans > 0. && DRIFT_FLOORS == 1) { // Strongly magnetized region; use drift frame floors
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

      double Padd[NVAR];
      PLOOP Padd[ip] = 0.;
      Padd[RHO] = MY_MAX(0., rhoflr - S->P[RHO][k][j][i]);
      Padd[UU] = MY_MAX(0., uflr - S->P[UU][k][j][i]);

      // Pull a /pretty bad/ hack to get 1-zone calculations
      // TODO come back to this and re-think larger architecture
      double P_bkp[NVAR], U_bkp[NVAR];
      PLOOP {
        P_bkp[ip] = S->P[ip][k][j][i];
        U_bkp[ip] = S->U[ip][k][j][i];

        S->P[ip][k][j][i] = Padd[ip];
        S->U[ip][k][j][i] = 0;
      }

      double Uadd[NVAR];
      get_state(G, S, i, j, k, CENT);
      prim_to_flux(G, S, i, j, k, 0, CENT, S->U);

      PLOOP {
        Uadd[ip] = S->U[ip][k][j][i];

        S->P[ip][k][j][i] = P_bkp[ip];
        S->U[ip][k][j][i] = U_bkp[ip];
      }

      get_state(G, S, i, j, k, CENT);
      prim_to_flux(G, S, i, j, k, 0, CENT, S->U);

      PLOOP {
        S->U[ip][k][j][i] += Uadd[ip];
        S->P[ip][k][j][i] += Padd[ip];
      }

      U_to_P(G, S, i, j, k, CENT);
      get_state(G, S, i, j, k, CENT);
    }
  }
#endif

#if ELECTRONS
  // Reset entropy after floors
  S->P[KTOT][k][j][i] = (gam - 1.)*S->P[UU][k][j][i]/pow(S-P[RHO][k][j][i],gam);

  // Set KTOTMAX to 3 by controlling u, to avoid anomalous cooling from funnel
  // wall
  double KTOTMAX = 3.;
  if (S->P[KTOT][k][j][i] > KTOTMAX) {
    S->P[UU][k][j][i] = KTOTMAX*pow(S-P[RHO][k][j][i],gam)/(gam-1.);
    S->P[KTOT][k][j][i] = KTOTMAX;
  }
#endif // ELECTRONS

  // Limit gamma with respect to normal observer
  double gamma = mhd_gamma_calc(G, S, i, j, k, CENT);

  if (gamma > GAMMAMAX) {
    double f = sqrt((GAMMAMAX*GAMMAMAX - 1.)/(gamma*gamma - 1.));
    S->P[U1][k][j][i] *= f;
    S->P[U2][k][j][i] *= f;
    S->P[U3][k][j][i] *= f;
    // Get state since we changed primitives
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

  int nfixed_utop = 0;
#pragma omp parallel for simd collapse(2) reduction (+:nfixed_utop)
  ZLOOP {
    // Count the 0 = bad cells
    nfixed_utop += !pflag[k][j][i];
  }

  // Make sure we are not using ill defined corner regions
  // I don't think 27 iterations is worth OpenMP
  for (int k = 0; k < NG; k++) {
    for (int j = 0; j < NG; j++) {
      for (int i = 0; i < NG; i++) {
        pflag[k][j][i] = 0;
        pflag[k][j][i+N1+NG] = 0;
        pflag[k][j+N2+NG][i] = 0;
        pflag[k+N3+NG][j][i] = 0;
        pflag[k][j+N2+NG][i+N1+NG] = 0;
        pflag[k+N3+NG][j][i+N1+NG] = 0;
        pflag[k+N3+NG][j+N2+NG][i] = 0;
        pflag[k+N3+NG][j+N2+NG][i+N1+NG] = 0;
      }
    }
  }

  // Fix the interior points first
  int bad = 0;
  do {
    bad = 0;
    // TODO is this okay? Not /quite/ deterministic...
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
          bad++;
          continue;
        }
        FLOOP S->P[ip][k][j][i] = sum[ip]/wsum;

        // Cell is fixed, can now use for other interpolations
        pflag[k][j][i] = 1;

#if !OLD_FLOORS
        get_state(G, S, i, j, k, CENT);
#endif
        fixup1zone(G, S, i, j, k);
      }
    }
  } while (bad > 0);

  LOGN("Fixed %d cells", nfixed_utop);

  // Flip back pflag
#pragma omp parallel for simd collapse(2)
  ZLOOPALL {
    pflag[k][j][i] = !pflag[k][j][i];
  }

  timer_stop(TIMER_FIXUP);
}
#undef FLOOP
