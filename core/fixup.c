/******************************************************************************
 *                                                                            *
 * FIXUP.C                                                                    *
 *                                                                            *
 * REPAIR INTEGRATION FAILURES                                                *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static int nfixed = 0;
static int nfixed_b = 0;

// Apply floors to density, internal energy
void fixup(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);
  #pragma omp parallel for collapse(3)
  ZLOOP fixup1zone(G, S, i, j, k);
  timer_stop(TIMER_FIXUP);
  printf("Fixed %d zones, %d from new floor\n", nfixed, nfixed_b);
  nfixed = 0;
  nfixed_b = 0;
}

void fixup1zone(struct GridGeom *G, struct FluidState *S, int i, int j, int k)
{
  double uuscal, rhoscal, rhoflr, uuflr;
  double f, gamma;

#if METRIC == MINKOWSKI
  rhoscal = 1.e-3;
  uuscal = 1.e-3;
#elif METRIC == MKS
  double r, th, X[NDIM];
  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);

  rhoscal = pow(r, -1.5);
  uuscal = rhoscal/r;
#endif // METRIC

  // Rho_min from T,N,M '10
  // TODO check factors, insert mass in local frame, retire solid floor
  get_state(G, S, i, j, k, CENT);
  double rho_min = bsq_calc(S, i, j, k) / (8*M_PI) / GAMMAAGN;

  rhoflr = MY_MAX(RHOMIN*rhoscal, rho_min*rhoscal);
  uuflr = UUMIN*uuscal;

  rhoflr = MY_MAX(rhoflr, RHOMINLIMIT);

  uuflr = MY_MAX(uuflr, UUMINLIMIT);

  // Floor on density and internal energy density (momentum *not* conserved)
  if (S->P[RHO][k][j][i] < rhoflr || isnan(S->P[RHO][k][j][i])) {
#if METRIC == MKS
    if (S->P[RHO][k][j][i] > RHOMIN*rhoscal) {
      nfixed_b++;
      //printf("Adding from B floor at r = %f, th = %f (%d %d %d): ", r, th, i, j, k);
    }
    nfixed++;
#endif
    S->P[RHO][k][j][i] = rhoflr;
  }
  if (S->P[UU][k][j][i] < uuflr || isnan(S->P[UU][k][j][i])) {
    S->P[UU][k][j][i] = uuflr;
  }

  // Limit gamma with respect to normal observer

  if (mhd_gamma_calc(G, S, i, j, k, CENT, &gamma)) {
    // Treat gamma failure here as "fixable" for fixup_utoprim()
    pflag[k][j][i] = -333;
  } else {
    if (gamma > GAMMAMAX) {
      f = sqrt((GAMMAMAX*GAMMAMAX - 1.) /(gamma*gamma - 1.));
      S->P[U1][k][j][i] *= f;
      S->P[U2][k][j][i] *= f;
      S->P[U3][k][j][i] *= f;
    }
  }
}

// Replace bad points with values interpolated from neighbors
#define FLOOP for(int ip=0;ip<B1;ip++)
void fixup_utoprim(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_FIXUP);
  int bad;
  double sum[B1], wsum;

  // Flip the logic of the pflag[] so that it now indicates which cells are good
  ZSLOOP(-NG, N3 - 1 + NG, -NG, N2 - 1 + NG, -NG, N1 - 1 + NG) {
    pflag[k][j][i] = !pflag[k][j][i];
  }

  // Make sure we are not using ill defined corner regions
  #pragma omp parallel for collapse(3)
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
        pflag[k+N3-1+NG][j+N2+NG][i+N1+NG] = 0;
      }
    }
  }

  // Fix the interior points first
  do {
    bad = 0;
    ZSLOOP(0, N3 - 1, 0, N2 - 1, 0, N1 - 1) {
      //if i == 6 && j == 4 
      //printf("[%i %i %i] N1 N2 N3 = %i %i %i pflag = %i\n", i,j,k,N1,N2,N3,pflag[k][j][i]);
      if (pflag[k][j][i] == 0) {
        wsum = 0.;
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

        fixup1zone(G, S, i, j, k);
      }
    }
  } while (bad > 0);
  timer_stop(TIMER_FIXUP);
}
#undef FLOOP

