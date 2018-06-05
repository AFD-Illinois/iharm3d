/******************************************************************************
 *                                                                            *
 * CURRENT.C                                                                  *
 *                                                                            *
 * CALCULATE CURRENT FROM FLUID VARIABLES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

double Fcon_calc(struct GridGeom *G, struct FluidState *S, int mu, int nu, int i, int j, int k);
int antisym(int a, int b, int c, int d);
int pp(int n, int *P);


static struct FluidState *Sa;

// The current calculation takes a /lot/ of time.  This is for cutting it out easily
void current_calc_nop(struct GridGeom *G, struct FluidState *S, struct FluidState *Ssave, double dtsave)
{}

void current_calc(struct GridGeom *G, struct FluidState *S, struct FluidState *Ssave, double dtsave)
{
  double gF0p[NDIM], gF0m[NDIM], gF1p[NDIM], gF1m[NDIM], gF2p[NDIM], gF2m[NDIM];
  double gF3p[NDIM], gF3m[NDIM];

  timer_start(TIMER_CURRENT);

  if (nstep == 0) {
#pragma omp parallel for simd collapse(2)
    ZLOOP {
      for (int mu = 0; mu < NDIM; mu++) S->jcon[mu][k][j][i] = 0.;
    }
    return;
  }

  static int first_run = 1;
  if (first_run) {
      //We only need the primitives, but this is fast
      Sa = calloc(1,sizeof(struct FluidState));
      first_run = 0;
  }

  // Calculate time-centered P
#pragma omp parallel for simd collapse(3)
  PLOOP {
    ZLOOPALL {
      Sa->P[ip][k][j][i] = 0.5*(S->P[ip][k][j][i] + Ssave->P[ip][k][j][i]);
    }
  }

  // Calculate j^{\mu} using centered differences for active zones
  // TODO rewrite this vector-style
#pragma omp parallel for collapse(3)
  ZLOOP {
    for (int mu = 0; mu < NDIM; mu++) S->jcon[mu][k][j][i] = 0.;

    // Get sqrt{-g}*F^{mu nu} at neighboring points

    // X0
    for (int mu = 0; mu < NDIM; mu++) {
      gF0p[mu] = Fcon_calc(G, S,  0, mu, i, j, k);
      gF0m[mu] = Fcon_calc(G, Ssave, 0, mu, i, j, k);
    }

    // X1
    for (int mu = 0; mu < NDIM; mu++) {
      gF1p[mu] = Fcon_calc(G, Sa,  1, mu, i+1, j, k);
      gF1m[mu] = Fcon_calc(G, Sa, 1, mu, i-1, j, k);
    }

    // X2
    for (int mu = 0; mu < NDIM; mu++) {
      gF2p[mu] = Fcon_calc(G, Sa,  2, mu, i, j+1, k);
      gF2m[mu] = Fcon_calc(G, Sa, 2, mu, i, j-1, k);
    }

    // X3
    for (int mu = 0; mu < NDIM; mu++) {
      gF3p[mu] = Fcon_calc(G, Sa,  3, mu, i, j, k+1);
      gF3m[mu] = Fcon_calc(G, Sa, 3, mu, i, j, k-1);
    }

    // Difference: D_mu F^{mu nu} = 4 \pi j^nu
    for (int mu = 0; mu < NDIM; mu++) {
	S->jcon[mu][k][j][i] = (1./(4.*M_PI*G->gdet[CENT][j][i]))*(
                           (gF0p[mu] - gF0m[mu])/dtsave +
                           (gF1p[mu] - gF1m[mu])/(2.*dx[1]) +
                           (gF2p[mu] - gF2m[mu])/(2.*dx[2]) +
                           (gF3p[mu] - gF3m[mu])/(2.*dx[3]));
    }
  }

  timer_stop(TIMER_CURRENT);
}

// Return mu, nu component of contravarient Maxwell tensor at grid zone i, j, k
inline double Fcon_calc(struct GridGeom *G, struct FluidState *S, int mu, int nu, int i, int j, int k)
{
  double Fcon, gFcon, dFcon;

  if (mu == nu) return 0.;

  get_state(G,S,i,j,k, CENT);

  Fcon = 0.;
  for (int kap = 0; kap < NDIM; kap++) {
    for (int lam = 0; lam < NDIM; lam++) {
      dFcon = (-1./G->gdet[CENT][j][i])*antisym(mu,nu,kap,lam)*S->ucov[kap][k][j][i]*S->bcov[lam][k][j][i];
      Fcon += dFcon;
    }
  }

  gFcon = Fcon*G->gdet[CENT][j][i];

  return gFcon;
}

// Completely antisymmetric 4D symbol
inline int antisym(int a, int b, int c, int d)
{
  // Check for valid permutation
  if (a < 0 || a > 3) return 100;
  if (b < 0 || b > 3) return 100;
  if (c < 0 || c > 3) return 100;
  if (d < 0 || d > 3) return 100;

  // Entries different? 
  if (a == b) return 0;
  if (a == c) return 0;
  if (a == d) return 0;
  if (b == c) return 0;
  if (b == d) return 0;
  if (c == d) return 0;

  // Determine parity of permutation        
  int p[4] = {a, b, c, d};
  
  return pp(4, p);
}

// Due to Norm Hardy; good for general n
inline int pp(int n, int P[n])
{
  int x;
  int p = 0;
  int v[n];

  for (int j = 0; j < n; j++) v[j] = 0;

  for (int j = 0; j < n; j++) {
    if (v[j]) {
      p++;
    } else {
      x = j;
      do {
        x = P[x];
        v[x] = 1;
      } while (x != j);
    }
  }

  if (p % 2 == 0) {
    return 1;
  } else {
    return -1;
  }
}

