/******************************************************************************
 *                                                                            *
 * METRIC.C                                                                   *
 *                                                                            *
 * HELPER FUNCTIONS FOR METRIC TENSORS                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// MHD stress-energy tensor with first index up, second index down. A factor of
// sqrt(4 pi) is absorbed into the definition of b.
inline void mhd_calc(struct FluidState *S, int i, int j, int k, int dir, double *mhd)
{
  double r, u, pres, w, bsq, eta, ptot;

  r = S->P[RHO][k][j][i];
  u = S->P[UU][k][j][i];
  pres = (gam - 1.)*u;
  w = pres + r + u;
  bsq = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    bsq += S->bcon[mu][k][j][i]*S->bcov[mu][k][j][i];
  }
  eta = w + bsq;
  ptot = pres + 0.5*bsq;

  for (int mu = 0; mu < NDIM; mu++) {
    mhd[mu] = eta*S->ucon[dir][k][j][i]*S->ucov[mu][k][j][i] +
              ptot*delta(dir, mu) -
              S->bcon[dir][k][j][i]*S->bcov[mu][k][j][i];
  }
}

// TODO OLD only used in fixup.c and even then hacked to hell
void prim_to_flux(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int dir, int loc, GridPrim flux)
{
  double mhd[NDIM];

  // Particle number flux
  flux[RHO][k][j][i] = S->P[RHO][k][j][i]*S->ucon[dir][k][j][i];

  mhd_calc(S, i, j, k, dir, mhd);

  // MHD stress-energy tensor w/ first index up, second index down
  flux[UU][k][j][i] = mhd[0] + flux[RHO][k][j][i];
  flux[U1][k][j][i] = mhd[1];
  flux[U2][k][j][i] = mhd[2];
  flux[U3][k][j][i] = mhd[3];

  // Dual of Maxwell tensor
  flux[B1][k][j][i] = S->bcon[1][k][j][i]*S->ucon[dir][k][j][i] -
                      S->bcon[dir][k][j][i]*S->ucon[1][k][j][i];
  flux[B2][k][j][i] = S->bcon[2][k][j][i]*S->ucon[dir][k][j][i] -
                      S->bcon[dir][k][j][i]*S->ucon[2][k][j][i];
  flux[B3][k][j][i] = S->bcon[3][k][j][i]*S->ucon[dir][k][j][i] -
                      S->bcon[dir][k][j][i]*S->ucon[3][k][j][i];

  PLOOP flux[ip][k][j][i] *= G->gdet[loc][j][i];
}

// Calculate fluxes in direction dir, over given range.
// Note backward indices convention, consistent with ZSLOOP's arguments
void prim_to_flux_vec(struct GridGeom *G, struct FluidState *S, int dir, int loc,
  int kstart, int kstop, int jstart, int jstop, int istart, int istop,
  GridPrim flux)
{
  //Only collapse over first two dimensions, to preserve vectorization
#pragma omp parallel for collapse(2)
  ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
    double mhd[NDIM];

    flux[RHO][k][j][i] = S->P[RHO][k][j][i] * S->ucon[dir][k][j][i];

    mhd_calc(S, i, j, k, dir, mhd);

    // MHD stress-energy tensor w/ first index up, second index down
    flux[UU][k][j][i] = mhd[0] + flux[RHO][k][j][i];
    flux[U1][k][j][i] = mhd[1];
    flux[U2][k][j][i] = mhd[2];
    flux[U3][k][j][i] = mhd[3];

    // Dual of Maxwell tensor
    flux[B1][k][j][i] = S->bcon[1][k][j][i] * S->ucon[dir][k][j][i]
		    - S->bcon[dir][k][j][i] * S->ucon[1][k][j][i];
    flux[B2][k][j][i] = S->bcon[2][k][j][i] * S->ucon[dir][k][j][i]
		    - S->bcon[dir][k][j][i] * S->ucon[2][k][j][i];
    flux[B3][k][j][i] = S->bcon[3][k][j][i] * S->ucon[dir][k][j][i]
		    - S->bcon[dir][k][j][i] * S->ucon[3][k][j][i];

    PLOOP flux[ip][k][j][i] *= G->gdet[loc][j][i];
  }
}

// calculate magnetic field four-vector
inline void bcon_calc(struct FluidState *S, int i, int j, int k)
{
  S->bcon[0][k][j][i] = S->P[B1][k][j][i]*S->ucov[1][k][j][i] +
                        S->P[B2][k][j][i]*S->ucov[2][k][j][i] +
                        S->P[B3][k][j][i]*S->ucov[3][k][j][i];
  for (int mu = 1; mu < 4; mu++) {
    S->bcon[mu][k][j][i] = (S->P[B1-1+mu][k][j][i] +
      S->bcon[0][k][j][i]*S->ucon[mu][k][j][i])/S->ucon[0][k][j][i];
  }
}

// Find gamma-factor wrt normal observer
inline int mhd_gamma_calc(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k, int loc, double *gamma)
{
  double qsq;

  qsq = G->gcov[loc][1][1][j][i]*S->P[U1][k][j][i]*S->P[U1][k][j][i]
      + G->gcov[loc][2][2][j][i]*S->P[U2][k][j][i]*S->P[U2][k][j][i]
      + G->gcov[loc][3][3][j][i]*S->P[U3][k][j][i]*S->P[U3][k][j][i]
      + 2.*(G->gcov[loc][1][2][j][i]*S->P[U1][k][j][i]*S->P[U2][k][j][i]
          + G->gcov[loc][1][3][j][i]*S->P[U1][k][j][i]*S->P[U3][k][j][i]
          + G->gcov[loc][2][3][j][i]*S->P[U2][k][j][i]*S->P[U3][k][j][i]);


  /*if (qsq < 0.) {
    if (fabs(qsq) > 1.E-10) { // Then assume not just machine precision
      fprintf(stderr,
        "gamma_calc():  failed: [%i %i %i] qsq = %28.18e \n",
        i, j, k, qsq);
      fprintf(stderr,
        "v[1-3] = %28.18e %28.18e %28.18e  \n",
        S->P[U1][k][j][i], S->P[U2][k][j][i], S->P[U3][k][j][i]);
      *gamma = 1.;
      return 1;
    } else {
      qsq = 1.E-10; // Set floor
    }
  }*/

  *gamma = sqrt(1. + qsq);

  return 0;

}

// Find contravariant four-velocity
inline void ucon_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc)
{
  double alpha, gamma;

  mhd_gamma_calc(G, S, i, j, k, loc, &gamma);

  alpha = G->lapse[loc][j][i];
  S->ucon[0][k][j][i] = gamma/alpha;
  for (int mu = 1; mu < NDIM; mu++) {
	S->ucon[mu][k][j][i] = S->P[U1+mu-1][k][j][i] -
						   gamma*alpha*G->gcon[loc][0][mu][j][i];
  }
}

// Calculate ucon, ucov, bcon, bcov from primitive variables
inline void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc)
{
    ucon_calc(G, S, i, j, k, loc);
    lower_grid(S->ucon, S->ucov, G, i, j, k, loc);
    //lower(q->ucon, geom, q->ucov);
    bcon_calc(S, i, j, k);
    //lower(q->bcon, geom, q->bcov);
    lower_grid(S->bcon, S->bcov, G, i, j, k, loc);
}

// Calculate ucon, ucov, bcon, bcov from primitive variables, over given range
// Note same range convention as ZSLOOP and other *_vec functions
void get_state_vec(struct GridGeom *G, struct FluidState *S, int loc,
  int kstart, int kstop, int jstart, int jstop, int istart, int istop)
{
#pragma omp parallel
  {
#pragma omp for simd collapse(2)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      ucon_calc(G, S, i, j, k, loc);
    }

#pragma omp for simd collapse(2)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      lower_grid(S->ucon, S->ucov, G, i, j, k, loc);
    }

#pragma omp for simd collapse(2)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      bcon_calc(S, i, j, k);
    }

#pragma omp for simd collapse(2)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      lower_grid(S->bcon, S->bcov, G, i, j, k, loc);
    }
  }
}

// Calculate components of magnetosonic velocity from primitive variables
//void mhd_vchar(double *Pr, struct of_state *q, struct of_geom *geom, int js,
//  double *vmax, double *vmin)
inline void mhd_vchar(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc, int dir, GridDouble cmax, GridDouble cmin)
{
  double discr, vp, vm, bsq, ee, ef, va2, cs2, cms2, rho, u;
  double Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM];
  double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;

  for (int mu = 0; mu < NDIM; mu++) {
    Acov[mu] = 0.;
  }
  Acov[dir] = 1.;
  //raise(Acov, geom, Acon);
  //raise_grid(Acov, Acon, G, i, j, k, loc);

  for (int mu = 0; mu < NDIM; mu++) {
    Bcov[mu] = 0.;
  }
  Bcov[0] = 1.;
  //raise(Bcov, geom, Bcon);
  //raise_grid(Bcov, Bcon, G, i, j, k, loc);
  for (int mu = 0; mu < NDIM; mu++) {
    Acon[mu] = 0.;
    Bcon[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      Acon[mu] += G->gcon[loc][mu][nu][j][i]*Acov[nu];
      Bcon[mu] += G->gcon[loc][mu][nu][j][i]*Bcov[nu];
    }
  }

  // Find fast magnetosonic speed
  //bsq = dot(q->bcon, q->bcov);
  bsq = dot_grid(S->bcon, S->bcov, i, j, k);
  //rho = fabs(Pr[RHO]);
  rho = fabs(S->P[RHO][k][j][i]);
  //u = fabs(Pr[UU]);
  u = fabs(S->P[UU][k][j][i]);
  ef = rho + gam*u;
  ee = bsq + ef;
  va2 = bsq/ee;
  cs2 = gam*(gam - 1.)*u/ef;

  cms2 = cs2 + va2 - cs2*va2;

  // Sanity checks
  if (cms2 < 0.) {
    fprintf(stderr, "\n\ncms2: %g %g %g\n\n", gam, u, ef);
    fail(FAIL_COEFF_NEG, i, j, k);
    cms2 = SMALL;
  }
  if (cms2 > 1.) {
    fail(FAIL_COEFF_SUP, i, j, k);
    cms2 = 1.;
  }

  // Require that speed of wave measured by observer q->ucon is cms2
  Asq = dot(Acon, Acov);
  Bsq = dot(Bcon, Bcov);
  Au = Bu = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    Au += Acov[mu]*S->ucon[mu][k][j][i];
    Bu += Bcov[mu]*S->ucon[mu][k][j][i];
  }
  //Au = dot_grid(Acov, S->ucon, i, j, k);
  //Bu = dot_grid(Bcov, S->ucon, i, j, k);
  AB = dot(Acon, Bcov);
  Au2 = Au*Au;
  Bu2 = Bu*Bu;
  AuBu = Au*Bu;

  A = Bu2 - (Bsq + Bu2)*cms2;
  B = 2.*(AuBu - (AB + AuBu)*cms2);
  C = Au2 - (Asq + Au2)*cms2;

  discr = B*B - 4.*A*C;
  if ((discr < 0.0) && (discr > -1.e-10)) {
    discr = 0.0;
  } else if (discr < -1.e-10) {
    fprintf(stderr, "\n\t %g %g %g %g %g\n", A, B, C, discr, cms2);
    fprintf(stderr, "\n\t S->ucon: %g %g %g %g\n", S->ucon[0][k][j][i],
      S->ucon[1][k][j][i], S->ucon[2][k][j][i], S->ucon[3][k][j][i]);
    fprintf(stderr, "\n\t S->bcon: %g %g %g %g\n", S->bcon[0][k][j][i],
      S->bcon[1][k][j][i], S->bcon[2][k][j][i], S->bcon[3][k][j][i]);
    fprintf(stderr, "\n\t Acon: %g %g %g %g\n", Acon[0], Acon[1], Acon[2],
      Acon[3]);
    fprintf(stderr, "\n\t Bcon: %g %g %g %g\n", Bcon[0], Bcon[1], Bcon[2],
      Bcon[3]);
    fail(FAIL_VCHAR_DISCR, i, j, k);
    discr = 0.;
  }

  discr = sqrt(discr);
  vp = -(-B + discr)/(2.*A);
  vm = -(-B - discr)/(2.*A);

  if (vp > vm) {
    cmax[k][j][i] = vp;
    cmin[k][j][i] = vm;
  } else {
    cmax[k][j][i] = vm;
    cmin[k][j][i] = vp;
  }
}

// Source terms for equations of motion
inline void get_fluid_source(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k, GridPrim dU)
{
  double mhd[NDIM][NDIM];

  mhd_calc(S, i, j, k, 0, mhd[0]);
  mhd_calc(S, i, j, k, 1, mhd[1]);
  mhd_calc(S, i, j, k, 2, mhd[2]);
  mhd_calc(S, i, j, k, 3, mhd[3]);

  // Contract mhd stress tensor with connection
  PLOOP dU[ip][k][j][i] = 0.;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      dU[UU][k][j][i] += mhd[mu][nu]*G->conn[nu][0][mu][j][i];
      dU[U1][k][j][i] += mhd[mu][nu]*G->conn[nu][1][mu][j][i];
      dU[U2][k][j][i] += mhd[mu][nu]*G->conn[nu][2][mu][j][i];
      dU[U3][k][j][i] += mhd[mu][nu]*G->conn[nu][3][mu][j][i];
    }
  }

  PLOOP dU[ip][k][j][i] *= G->gdet[CENT][j][i];
}

// Returns b.b (twice magnetic pressure)
double bsq_calc(struct FluidState *S, int i, int j, int k)
{

  double bsq = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    bsq += S->bcon[mu][k][j][i]*S->bcov[mu][k][j][i];
  }

  return bsq;
}
