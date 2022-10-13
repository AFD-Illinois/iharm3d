/******************************************************************************
 *                                                                            *
 * PHYS.C                                                                     *
 *                                                                            *
 * COMPUTES PHYSICAL QUANTITIES: STRESS-ENERGY TENSOR, U(P), FLUXES,          * 
 * 4-VECTORS, LORENTZ FACTOR, MAGNETOSONIC VELOCITY, SOURCE TERMS             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if EMHD
void gradient_calc(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k,
  double grad_ucov[NDIM][NDIM], double grad_Theta[NDIM]);
#endif

// MHD stress-energy tensor with first index up, second index down. A factor of
// sqrt(4 pi) is absorbed into the definition of b.
inline void mhd_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int k, int dir, double *mhd)
{
  double u, pres, w, bsq, eta, ptot;

  u = S->P[UU][k][j][i];
  pres = (gam - 1.)*u;
  w = pres + S->P[RHO][k][j][i] + u;
  bsq = MY_MAX(bsq_calc(S, i, j, k), SMALL);
  eta = w + bsq;
  ptot = pres + 0.5*bsq;

  #if EMHD
  // Use closure relations to compute tau, chi and nu
  set_emhd_parameters(G, S, i, j, k);
  #endif

  DLOOP1 {
    mhd[mu] = eta*S->ucon[dir][k][j][i]*S->ucov[mu][k][j][i] +
              ptot*delta(dir, mu) -
              S->bcon[dir][k][j][i]*S->bcov[mu][k][j][i];
  }

  #if EMHD
  double ucon    = S->ucon[dir][k][j][i];
  double bcon    = S->bcon[dir][k][j][i];
  #if CONDUCTION
  double q = S->q[k][j][i];
  DLOOP1 {
    double bcov = S->bcov[mu][k][j][i];
    double ucov = S->ucov[mu][k][j][i];
    mhd[mu] += (q / sqrt(bsq)) * ((ucon * bcov) + (bcon * ucov));
  }
  #endif
  #if VISCOSITY
  double delta_p = S->delta_p[k][j][i];
  DLOOP1 {
    double bcov = S->bcov[mu][k][j][i];
    double ucov = S->ucov[mu][k][j][i];
    mhd[mu] += (-delta_p) * ((bcon * bcov / bsq) - (1./3.) * (delta(dir, mu) + ucon * ucov));
  }
  #endif
  #endif
}

// TODO OLD only used in fixup.c and even then hacked to hell
void prim_to_flux(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int dir, int loc, GridPrim flux)
{
  double mhd[NDIM];

  // Particle number flux
  flux[RHO][k][j][i] = S->P[RHO][k][j][i]*S->ucon[dir][k][j][i];

  mhd_calc(G, S, i, j, k, dir, mhd);

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

  #if EMHD
  #if CONDUCTION
  flux[Q_TILDE][k][j][i] = S->P[Q_TILDE][k][j][i] * S->ucon[dir][k][j][i];
  #endif
  #if VISCOSITY
  flux[DELTA_P_TILDE][k][j][i] = S->P[DELTA_P_TILDE][k][j][i] * S->ucon[dir][k][j][i];
  #endif
  #endif

#if ELECTRONS
  for (int ip = KTOT  + 1; ip < NVAR ; ip++) {
    flux[ip][k][j][i] = flux[RHO][k][j][i]*S->P[ip][k][j][i];
  }
  flux[KTOT][k][j][i] = flux[RHO][k][j][i]*S->P[KTOT][k][j][i];
#endif

  PLOOP flux[ip][k][j][i] *= G->gdet[loc][j][i];
}

// Calculate fluxes in direction dir, over given range.
// Note backward indices convention, consistent with ZSLOOP's arguments
void prim_to_flux_vec(struct GridGeom *G, struct FluidState *S, int dir, int loc,
  int kstart, int kstop, int jstart, int jstop, int istart, int istop,
  GridPrim flux)
{
  // TODO reintroduce simd pragma to see where it messes things up
#pragma omp parallel
{
#pragma omp for collapse(3) nowait
  ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
    double mhd[NDIM];

    flux[RHO][k][j][i] = S->P[RHO][k][j][i] * S->ucon[dir][k][j][i] * G->gdet[loc][j][i];

    mhd_calc(G, S, i, j, k, dir, mhd);

    // MHD stress-energy tensor w/ first index up, second index down
    flux[UU][k][j][i] = mhd[0] * G->gdet[loc][j][i] + flux[RHO][k][j][i];
    flux[U1][k][j][i] = mhd[1] * G->gdet[loc][j][i];
    flux[U2][k][j][i] = mhd[2] * G->gdet[loc][j][i];
    flux[U3][k][j][i] = mhd[3] * G->gdet[loc][j][i];
  }

#pragma omp for collapse(3) nowait
  ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
    // Dual of Maxwell tensor
    flux[B1][k][j][i] = (S->bcon[1][k][j][i] * S->ucon[dir][k][j][i]
        - S->bcon[dir][k][j][i] * S->ucon[1][k][j][i]) * G->gdet[loc][j][i];
    flux[B2][k][j][i] = (S->bcon[2][k][j][i] * S->ucon[dir][k][j][i]
        - S->bcon[dir][k][j][i] * S->ucon[2][k][j][i]) * G->gdet[loc][j][i];
    flux[B3][k][j][i] = (S->bcon[3][k][j][i] * S->ucon[dir][k][j][i]
        - S->bcon[dir][k][j][i] * S->ucon[3][k][j][i]) * G->gdet[loc][j][i];
  }

  #if EMHD
  #pragma omp for collapse(3) nowait
  ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
    #if CONDUCTION
    flux[Q_TILDE][k][j][i] = G->gdet[loc][j][i] * (S->P[Q_TILDE][k][j][i] * S->ucon[dir][k][j][i]);
    #endif
    #if VISCOSITY
    flux[DELTA_P_TILDE][k][j][i] = G->gdet[loc][j][i] * (S->P[DELTA_P_TILDE][k][j][i] * S->ucon[dir][k][j][i]);
    #endif
  }
  #endif

#if ELECTRONS
#pragma omp for collapse(3)
  ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
    // RHO already includes a factor of gdet!
    for (int ip = KTOT + 1; ip < NVAR ; ip++) {
      flux[ip][k][j][i] = flux[RHO][k][j][i]*S->P[ip][k][j][i];
    }
    flux[KTOT][k][j][i] = flux[RHO][k][j][i]*S->P[KTOT][k][j][i];
  }
#endif
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
inline double mhd_gamma_calc(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k, int loc)
{
  double qsq = G->gcov[loc][1][1][j][i]*S->P[U1][k][j][i]*S->P[U1][k][j][i]
      + G->gcov[loc][2][2][j][i]*S->P[U2][k][j][i]*S->P[U2][k][j][i]
      + G->gcov[loc][3][3][j][i]*S->P[U3][k][j][i]*S->P[U3][k][j][i]
      + 2.*(G->gcov[loc][1][2][j][i]*S->P[U1][k][j][i]*S->P[U2][k][j][i]
          + G->gcov[loc][1][3][j][i]*S->P[U1][k][j][i]*S->P[U3][k][j][i]
          + G->gcov[loc][2][3][j][i]*S->P[U2][k][j][i]*S->P[U3][k][j][i]);


#if DEBUG
  if (qsq < 0.) {
    if (fabs(qsq) > 1.E-10) { // Then assume not just machine precision
      fprintf(stderr,
        "gamma_calc():  failed: [%i %i %i] qsq = %28.18e \n",
        i, j, k, qsq);
      fprintf(stderr,
        "v[1-3] = %28.18e %28.18e %28.18e  \n",
        S->P[U1][k][j][i], S->P[U2][k][j][i], S->P[U3][k][j][i]);
      return 1.0;
    } else {
      qsq = 1.E-10; // Set floor
    }
  }
#endif

  return sqrt(1. + qsq);

}

// Find contravariant four-velocity
inline void ucon_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc)
{
  double gamma = mhd_gamma_calc(G, S, i, j, k, loc);

  double alpha = G->lapse[loc][j][i];
  S->ucon[0][k][j][i] = gamma/alpha;
  for (int mu = 1; mu < NDIM; mu++) {
    S->ucon[mu][k][j][i] = S->P[U1+mu-1][k][j][i] -
        gamma*alpha*G->gcon[loc][0][mu][j][i];
  }
}

// Calculate ucon, ucov, bcon, bcov from primitive variables
// If IMEX is enabled, compute q, delta_p, bsq, Theta
inline void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc)
{
  ucon_calc(G, S, i, j, k, loc);
  lower_grid(S->ucon, S->ucov, G, i, j, k, loc);
  bcon_calc(S, i, j, k);
  lower_grid(S->bcon, S->bcov, G, i, j, k, loc);

  // Need auxillary scalar fields to (i) set emhd parameters (ii) compute components of stress-energy tensor
  #if EMHD
  S->bsq[k][j][i]     = MY_MAX(bsq_calc(S, i, j, k), SMALL);
  S->Theta[k][j][i]   = MY_MAX((gam - 1.) * S->P[UU][k][j][i] / S->P[RHO][k][j][i], SMALL);
  
  #if CONDUCTION
  S->q[k][j][i] = S->P[Q_TILDE][k][j][i];

  if (higher_order_terms_conduction == 1) {
    
    // Initializations
    double rho   = S->P[RHO][k][j][i];
    double u     = S->P[UU][k][j][i];
    double pg     = (gam - 1.)*u;
    double cs2   = gam*pg/(rho + gam*u);
    double Theta = S->Theta[k][j][i]; // temperature has been updated above

    S->q[k][j][i] *= sqrt(rho * conduction_alpha * cs2 * pow(Theta, 2));

  }
  #endif
  #if VISCOSITY
  S->delta_p[k][j][i] = S->P[DELTA_P_TILDE][k][j][i];
  
  if (higher_order_terms_viscosity == 1) {
    
    // Initializations
    double rho   = S->P[RHO][k][j][i];
    double u     = S->P[UU][k][j][i];
    double pg     = (gam - 1.)*u;
    double cs2   = gam*pg/(rho + gam*u);
    double Theta = S->Theta[k][j][i]; // temperature has been updated above

    S->delta_p[k][j][i] *= sqrt(rho * viscosity_alpha * cs2 * Theta);

  }
  #endif
  #endif

}

// Calculate ucon, ucov, bcon, bcov from primitive variables, over given range
// Note same range convention as ZSLOOP and other *_vec functions
void get_state_vec(struct GridGeom *G, struct FluidState *S, int loc,
  int kstart, int kstop, int jstart, int jstop, int istart, int istop)
{
#pragma omp parallel
  {
#pragma omp for collapse(3)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      ucon_calc(G, S, i, j, k, loc);
      //lower_grid(S->ucon, S->ucov, G, i, j, k, loc);
    }

#pragma omp for collapse(3)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      lower_grid(S->ucon, S->ucov, G, i, j, k, loc);
    }
    //lower_grid_vec(S->ucon, S->ucov, G, kstart, kstop, jstart, jstop, istart, istop, loc);

#pragma omp for collapse(3)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      bcon_calc(S, i, j, k);
      //lower_grid(S->bcon, S->bcov, G, i, j, k, loc);
    }

#pragma omp for collapse(3)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      lower_grid(S->bcon, S->bcov, G, i, j, k, loc);
    }

    #if EMHD
    #pragma omp for collapse(3)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      S->bsq[k][j][i]     = MY_MAX(bsq_calc(S, i, j, k), SMALL);
      S->Theta[k][j][i]   = MY_MAX((gam - 1.) * S->P[UU][k][j][i] / S->P[RHO][k][j][i], SMALL);
      
      #if CONDUCTION
      S->q[k][j][i] = S->P[Q_TILDE][k][j][i];

      if (higher_order_terms_conduction == 1) {

        // Initializations
        double rho   = S->P[RHO][k][j][i];
        double u     = S->P[UU][k][j][i];
        double pg    = (gam - 1.)*u;
        double cs2   = gam*pg/(rho + gam*u);
        double Theta = S->Theta[k][j][i];

        S->q[k][j][i] *= sqrt(rho * conduction_alpha * cs2 * pow(Theta, 2));
      }
      #endif
      #if VISCOSITY
      S->delta_p[k][j][i] = S->P[DELTA_P_TILDE][k][j][i];

      if (higher_order_terms_viscosity == 1) {

        // Initializations
        double rho   = S->P[RHO][k][j][i];
        double u     = S->P[UU][k][j][i];
        double pg    = (gam - 1.)*u;
        double cs2   = gam*pg/(rho + gam*u);
        double Theta = S->Theta[k][j][i];

        S->delta_p[k][j][i] *= sqrt(rho * viscosity_alpha * cs2 * Theta);
      }
      #endif      
    }
    #endif

  }
}

// Calculate components of magnetosonic velocity from primitive variables
// TODO this is a primary candidate for splitting/vectorizing
inline void mhd_vchar(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc, int dir, GridDouble cmax, GridDouble cmin)
{
  double discr, vp, vm, bsq, ee, ef, va2, cs2, cms2, rho, u;
  double Acov[NDIM], Bcov[NDIM], Acon[NDIM], Bcon[NDIM]; 
  double Asq, Bsq, Au, Bu, AB, Au2, Bu2, AuBu, A, B, C;

  DLOOP1 {
    Acov[mu] = 0.;
  }
  Acov[dir] = 1.;

  DLOOP1 {
    Bcov[mu] = 0.;
  }
  Bcov[0] = 1.;

  DLOOP1 {
    Acon[mu] = 0.;
    Bcon[mu] = 0.;
  }
  DLOOP2 {
    Acon[mu] += G->gcon[loc][mu][nu][j][i]*Acov[nu];
    Bcon[mu] += G->gcon[loc][mu][nu][j][i]*Bcov[nu];
  }

  // Find fast magnetosonic speed
  bsq = MY_MAX(bsq_calc(S, i, j, k), SMALL);
  rho = S->P[RHO][k][j][i];
  u = S->P[UU][k][j][i];
  ef = MY_MAX(rho + gam*u, SMALL);
  ee = bsq + ef;
  va2 = bsq/ee;
  cs2 = gam*(gam - 1.)*u/ef;

  #if EMHD
  double tau    = S->tau[k][j][i];
  double ccond2 = 0.;
  double cvis2  = 0.;

  #if CONDUCTION
  double chi_emhd = S->chi_emhd[k][j][i];
  // wave speed contribution due to q
  ccond2 = (gam - 1.) * chi_emhd / tau;
  #endif
  #if VISCOSITY
  double nu_emhd  = S->nu_emhd[k][j][i];
  // wave speed contribution due to dP
  cvis2 = (4. / 3.) / (rho + (gam * u)) * rho * nu_emhd / tau;
  #endif

  ccond2 = 0.5*(cs2 + ccond2 + sqrt(cs2*cs2 + ccond2*ccond2));
  cs2 = ccond2 + cvis2;
  #endif

  cms2 = cs2 + va2 - cs2*va2;

  cms2 = (cms2 < 0) ? SMALL : cms2;
  cms2 = (cms2 > 1) ? 1 : cms2;

  // Require that speed of wave measured by observer q->ucon is cms2
  Asq = dot(Acon, Acov);
  Bsq = dot(Bcon, Bcov);
  Au = Bu = 0.;
  DLOOP1 {
    Au += Acov[mu]*S->ucon[mu][k][j][i];
    Bu += Bcov[mu]*S->ucon[mu][k][j][i];
  }
  AB = dot(Acon, Bcov);
  Au2 = Au*Au;
  Bu2 = Bu*Bu;
  AuBu = Au*Bu;

  A = Bu2 - (Bsq + Bu2)*cms2;
  B = 2.*(AuBu - (AB + AuBu)*cms2);
  C = Au2 - (Asq + Au2)*cms2;

  discr = B*B - 4.*A*C;
  discr = (discr < 0.) ? 0. : discr;
  discr = sqrt(discr);

  vp = -(-B + discr)/(2.*A);
  vm = -(-B - discr)/(2.*A);

  cmax[k][j][i] = (vp > vm) ? vp : vm;
  cmin[k][j][i] = (vp > vm) ? vm : vp;
}

// Wrapper for HARM algo call to compute ideal GRMHD source terms. Calls explicit_sources
// Having said that, this function contains a wind source term (if enabled), which is not present in the explicit_sources function
inline void get_ideal_fluid_source_vec(struct GridGeom *G, struct FluidState *S, GridPrim *dU)
{
#if WIND_TERM
  static struct FluidState *dS;
  static int firstc = 1;
  if (firstc) {dS = calloc(1, sizeof(struct FluidState)); firstc = 0;}
#endif

#pragma omp parallel for collapse(3)
  ZLOOP {
    double dU_per_zone[NVAR] = {0};
    explicit_sources(G, S, CENT, i, j, k, dU_per_zone);
    PLOOP (*dU)[ip][k][j][i] = dU_per_zone[ip];
  }

  // Add a small "wind" source term in RHO,UU
#if WIND_TERM
#pragma omp parallel for simd collapse(2)
  ZLOOP {
    // Stolen shamelessly from iharm2d_v3

    /* need coordinates to evaluate particle addtn rate */
    double X[NDIM];
    coord(i, j, k, CENT, X);
    double r, th;
    bl_coord(X, &r, &th);
    double cth = cos(th) ;

    /* here is the rate at which we're adding particles */
    /* this function is designed to concentrate effect in the
     funnel in black hole evolutions */
    double drhopdt = 2.e-4*cth*cth*cth*cth/pow(1. + r*r,2) ;

    dS->P[RHO][k][j][i] = drhopdt ;

    double Tp = 10. ;  /* temp, in units of c^2, of new plasma */
    dS->P[UU][k][j][i] = drhopdt*Tp*3. ;

    /* Leave P[U{1,2,3}]=0 to add in particles in normal observer frame */
    /* Likewise leave P[BN]=0 */
  }

  /* add in plasma to the T^t_a component of the stress-energy tensor */
  /* notice that U already contains a factor of sqrt{-g} */
  get_state_vec(G, dS, CENT, 0, N3-1, 0, N2-1, 0, N1-1);
  prim_to_flux_vec(G, dS, 0, CENT, 0, N3-1, 0, N2-1, 0, N1-1, dS->U);

#pragma omp parallel for simd collapse(3)
  PLOOP ZLOOP {
    (*dU)[ip][k][j][i] += dS->U[ip][k][j][i] ;
  }
#endif
}

// Compute explicit source terms
// NOTE: This function is outside the IMEX conditional directive since the HARM algo calls it (through get_ideal_fluid_source_vec)
double explicit_sources(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k, double dU[]) {
  
  // Ideal MHD components
  double mhd[NDIM][NDIM];

  DLOOP1 mhd_calc(G, S, i, j, k, mu, mhd[mu]);

  FLOOP dU[ip] = 0.; // explicit_sources should be the first source call
  DLOOP2 {
    for (int gam = 0; gam < NDIM; gam++)
      dU[UU+gam] += mhd[mu][nu]*G->conn[nu][gam][mu][j][i];
  }
  FLOOP dU[ip] *= G->gdet[CENT][j][i];

  // Extended MHD components
  #if EMHD
  // Initializations

  double rho      = S->P[RHO][k][j][i];
  double Theta    = S->Theta[k][j][i];
  double bsq      = S->bsq[k][j][i];
  double tau      = S->tau[k][j][i];
  double gdet     = G->gdet[loc][j][i];

  // Compute gradient of ucov and Theta
  double grad_ucov[NDIM][NDIM] = {0}, grad_Theta[NDIM] = {0};
  gradient_calc(G, S, loc, i, j, k, grad_ucov, grad_Theta);

  // Compute div of ucon (all terms but the time-derivative ones are nonzero)
  double div_ucon = 0;
  DLOOP2 {
    double gcon_mu_nu = G->gcon[loc][mu][nu][j][i];

    div_ucon += gcon_mu_nu * grad_ucov[mu][nu];
  }

  // Compute q0 and deltaP0 (everything but the time-derivative terms)
  #if CONDUCTION
  double chi_emhd = S->chi_emhd[k][j][i];
  double q_tilde  = S->P[Q_TILDE][k][j][i];
  
  double q0 = 0.;

  DLOOP1 {
    double bcon_mu = S->bcon[mu][k][j][i];

    q0 -= rho * chi_emhd * (bcon_mu / sqrt(bsq)) * grad_Theta[mu];
  }

  DLOOP2 {
    double bcon_mu = S->bcon[mu][k][j][i];
    double ucon_nu = S->ucon[nu][k][j][i];

    q0 -= rho * chi_emhd * (bcon_mu / sqrt(bsq)) * Theta * ucon_nu * grad_ucov[nu][mu];
  }

  // Higher order correction to the relaxed value
  double q0_tilde = 0.;

  if (higher_order_terms_conduction == 1)
    q0_tilde = q0 * sqrt(tau / (chi_emhd * rho * pow(Theta, 2)));
  else
    q0_tilde = q0;

  // Add explicit source term (conduction)
  dU[Q_TILDE] = gdet * (q0_tilde / tau);

  // Higher order correction to the source term
  if (higher_order_terms_conduction == 1)
    dU[Q_TILDE] += gdet * (q_tilde / 2) * div_ucon;

  #endif

  #if VISCOSITY
  double nu_emhd       = S->nu_emhd[k][j][i];
  double delta_p_tilde = S->P[DELTA_P_TILDE][k][j][i];

  double deltaP0 = 0.;

  deltaP0 = -rho * nu_emhd * div_ucon;
  DLOOP2  {
    double bcon_mu = S->bcon[mu][k][j][i];
    double bcon_nu = S->bcon[nu][k][j][i];
    
    deltaP0 += 3. * rho * nu_emhd * (bcon_mu * bcon_nu / bsq) * grad_ucov[mu][nu];
  }

  //Higher order correction to the relaxed value
  double deltaP0_tilde = 0.;

  if (higher_order_terms_viscosity == 1) 
    deltaP0_tilde = deltaP0 * sqrt(tau / (nu_emhd * rho * Theta));
  else
    deltaP0_tilde = deltaP0;

  // Add explicit source term (viscosity)
  dU[DELTA_P_TILDE] = gdet * (deltaP0_tilde / tau);

  // Higher order correction to the source term
  if (higher_order_terms_viscosity == 1)
    dU[DELTA_P_TILDE] += gdet * (delta_p_tilde / 2) * div_ucon;

  #endif

  #endif
  return 0.0;
}

// Compute source terms with time derivatives
#if IMEX
double time_derivative_sources(struct GridGeom *G, struct FluidState *S_new, struct FluidState *S_old,
  struct FluidState *S, double dt, int loc, int i, int j, int k, double dU[]) {

  #if EMHD
  // Initializations
  double rho   = S->P[RHO][k][j][i];
  double Theta = S->Theta[k][j][i];
  double bsq   = S->bsq[k][j][i];
  double tau   = S->tau[k][j][i];
  double gdet  = G->gdet[loc][j][i];

  // Compute partial derivative of ucov
  double dt_ucov[NDIM] = {0};
  DLOOP1 {
    double ucov_new = S_new->ucov[mu][k][j][i];
    double ucov_old = S_old->ucov[mu][k][j][i];

    dt_ucov[mu] = (ucov_new - ucov_old) / dt;
  }

  // Compute div of ucon (only temporal part is nonzero)
  double div_ucon = 0.;
  DLOOP1 {
    double gcon_t_mu = G->gcon[loc][0][mu][j][i];

    div_ucon += gcon_t_mu * dt_ucov[mu];
  }

  double Theta_new, Theta_old, dt_Theta;
  Theta_new = S_new->Theta[k][j][i];
  Theta_old = S_old->Theta[k][j][i];

  dt_Theta = (Theta_new - Theta_old) / dt;

  double ucon_t  = S->ucon[0][k][j][i];
  double bcon_t  = S->bcon[0][k][j][i];

  // Compute q0 and delta_P0 (temporal terms)
  #if CONDUCTION
  double chi_emhd = S->chi_emhd[k][j][i];
  double q_tilde  = S->P[Q_TILDE][k][j][i];

  double q0 = 0.;

  q0 = -rho * chi_emhd * (bcon_t / sqrt(bsq)) * dt_Theta;
  DLOOP1 {
    double bcon_mu = S->bcon[mu][k][j][i];

    q0 -= rho * chi_emhd * (bcon_mu / sqrt(bsq)) * Theta * ucon_t * dt_ucov[mu];
  }

  // Higher order correction to the relaxed value
  double q0_tilde = 0.;
  if (higher_order_terms_conduction == 1)    
    q0_tilde = q0 * sqrt(tau / (chi_emhd * rho * pow(Theta, 2)));
  else
    q0_tilde = q0;

  // Add the time derivative source term (conduction)
  dU[Q_TILDE] = gdet * (q0_tilde / tau);

  // Higher order correction to the source term
  if (higher_order_terms_conduction == 1)
    dU[Q_TILDE] += gdet * (q_tilde / 2.) * div_ucon;

  #endif
  
  #if VISCOSITY
  double nu_emhd  = S->nu_emhd[k][j][i];
  double delta_p_tilde = S->P[DELTA_P_TILDE][k][j][i];

  double deltaP0 = 0.;

  deltaP0 = -rho * nu_emhd * div_ucon;
  DLOOP1 {
    double bcon_mu = S->bcon[mu][k][j][i];

    deltaP0 += 3. * rho * nu_emhd * (bcon_t * bcon_mu / bsq) * dt_ucov[mu];
  }

  // Higher order correction to the relaxed value
  double deltaP0_tilde = 0.;
  if (higher_order_terms_viscosity == 1)    
    deltaP0_tilde = deltaP0 * sqrt(tau / (nu_emhd * rho * Theta));
  else
    deltaP0_tilde = deltaP0;

  // Add the time derivative source term (viscosity)
  dU[DELTA_P_TILDE] = gdet * (deltaP0_tilde / tau);

  // Higher order correction to the source term
  if (higher_order_terms_viscosity == 1)
    dU[DELTA_P_TILDE] += gdet * (delta_p_tilde / 2.) * div_ucon;

  #endif
  
  #endif
  return 0.0;
}

// Compute implicit source terms
double implicit_sources(struct GridGeom *G, struct FluidState *S, struct FluidState *S_tau, int loc, int i, int j, int k, double dU[]){

  #if EMHD
  // Initializations
  double tau  = S_tau->tau[k][j][i];
  double gdet = G->gdet[loc][j][i];

  #if CONDUCTION
  double q_tilde = S->P[Q_TILDE][k][j][i];
  dU[Q_TILDE]    = -gdet * (q_tilde / tau);
  #endif

  #if VISCOSITY
  double delta_p_tilde = S->P[DELTA_P_TILDE][k][j][i];
  dU[DELTA_P_TILDE]    = -gdet * (delta_p_tilde / tau);
  #endif

  #endif
  return 0.0;
}

// Compute gradient of four velocities and temperature
// Called by emhd_explicit_sources
void gradient_calc(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k,
  double grad_ucov[NDIM][NDIM], double grad_Theta[NDIM]) {

  // Compute gradient of ucov
  DLOOP1 {
    grad_ucov[0][mu] = 0;

    grad_ucov[1][mu] = slope_calc_four_vec(S->ucov, mu, 1, i, j, k);

    grad_ucov[2][mu] = slope_calc_four_vec(S->ucov, mu, 2, i, j, k);

    grad_ucov[3][mu] = slope_calc_four_vec(S->ucov, mu, 3, i, j, k);
  }

  DLOOP2 {
    for (int gam = 0; gam < NDIM; gam++)
      grad_ucov[mu][nu] -= G->conn[gam][mu][nu][j][i] * S->ucov[gam][k][j][i];
  }

  #if EMHD
  // Compute temperature gradient (need only if EMHD is enabled)
  // Time derivative component computed in time_derivative_sources
  grad_Theta[0] = 0;

  grad_Theta[1] = slope_calc_scalar(S->Theta, 1, i, j, k);

  grad_Theta[2] = slope_calc_scalar(S->Theta, 2, i, j, k);

  grad_Theta[3] = slope_calc_scalar(S->Theta, 3, i, j, k);
  #endif
}
#endif

// Returns b.b (twice magnetic pressure)
inline double bsq_calc(struct FluidState *S, int i, int j, int k)
{

  double bsq = 0.;
  DLOOP1 {
    bsq += S->bcon[mu][k][j][i]*S->bcov[mu][k][j][i];
  }

  return bsq;
}
