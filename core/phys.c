/******************************************************************************
 *                                                                            *
 * PHYS.C                                                                     *
 *                                                                            *
 * COMPUTES PHYSICAL QUANTITIES: STRESS-ENERGY TENSOR, U(P), FLUXES,          * 
 * 4-VECTORS, LORENTZ FACTOR, MAGNETOSONIC VELOCITY, SOURCE TERMS             *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if GRIM_TIMESTEPPER
void gradient_calc(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k,
  double grad_ucov[NDIM][NDIM], double grad_Theta[NDIM]);
#endif

// MHD stress-energy tensor with first index up, second index down. A factor of
// sqrt(4 pi) is absorbed into the definition of b.
inline void mhd_calc(struct FluidState *S, int i, int j, int k, int dir, double *mhd)
{
  double u, pres, w, bsq, eta, ptot;

  u = S->P[UU][k][j][i];
  pres = (gam - 1.)*u;
  w = pres + S->P[RHO][k][j][i] + u;
  bsq = bsq_calc(S, i, j, k);
  eta = w + bsq;
  ptot = pres + 0.5*bsq;

  DLOOP1 {
    mhd[mu] = eta*S->ucon[dir][k][j][i]*S->ucov[mu][k][j][i] +
              ptot*delta(dir, mu) -
              S->bcon[dir][k][j][i]*S->bcov[mu][k][j][i];
  }

  #if GRIM_TIMESTEPPER
  double q       = S->q[k][j][i];
  double delta_p = S->delta_p[k][j][i];
  double ucon    = S->ucon[dir][k][j][i];
  double bcon    = S->bcon[dir][k][j][i];
  
  DLOOP1 {
    double bcov = S->bcov[mu][k][j][i];
    double ucov = S->ucov[mu][k][j][i];

    mhd[mu] += (q / sqrt(bsq)) * ((ucon * bcov) + (bcon * ucov)) 
            + (-delta_p) * ((bcon * bcov / bsq) - (1./3.) * (delta(dir, mu) + ucon * ucov));
  }
  #endif
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

  #if GRIM_TIMESTEPPER
  flux[Q_TILDE][k][j][i]       = S->ucon[dir][k][j][i] * S->P[Q_TILDE][k][j][i];
  flux[DELTA_P_TILDE][k][j][i] = S->ucon[dir][k][j][i] * S->P[DELTA_P_TILDE][k][j][i];
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

    mhd_calc(S, i, j, k, dir, mhd);

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

  #if GRIM_TIMESTEPPER
  #pragma omp for collapse(3) nowait
  ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
    flux[Q_TILDE][k][j][i] = G->gdet[loc][j][i] * (S->ucon[dir][k][j][i] * S->P[Q_TILDE][k][j][i]);
    flux[DELTA_P_TILDE][k][j][i] = G->gdet[loc][j][i] * (S->ucon[dir][k][j][i] * S->P[DELTA_P_TILDE][k][j][i]);
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
// If GRIM_TIMESTEPPER is enabled, compute q, delta_p, bsq, Theta
inline void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc)
{
    ucon_calc(G, S, i, j, k, loc);
    lower_grid(S->ucon, S->ucov, G, i, j, k, loc);
    bcon_calc(S, i, j, k);
    lower_grid(S->bcon, S->bcov, G, i, j, k, loc);

    #if GRIM_TIMESTEPPER
    S->q[k][j][i]       = S->P[Q_TILDE][k][j][i]; // NOTE: Will have to edit this when considering higher order terms
    S->delta_p[k][j][i] = S->P[DELTA_P_TILDE][k][j][i];
    S->bsq[k][j][i]     = bsq_calc(S, i, j, k);
    S->Theta[k][j][i]   = (gam - 1.) * S->P[UU][k][j][i] / S->P[RHO][k][j][i];
    
    set_emhd_parameters(G, S, i, j, k);
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

    #if GRIM_TIMESTEPPER
    #pragma omp for collapse(3)
    ZSLOOP(kstart, kstop, jstart, jstop, istart, istop) {
      S->q[k][j][i]       = S->P[Q_TILDE][k][j][i]; // NOTE: Will have to edit this when considering higher order terms
      S->delta_p[k][j][i] = S->P[DELTA_P_TILDE][k][j][i];
      S->bsq[k][j][i]     = bsq_calc(S, i, j, k);
      S->Theta[k][j][i]   = (gam - 1.) * S->P[UU][k][j][i] / S->P[RHO][k][j][i];
      
      set_emhd_parameters(G, S, i, j, k);
    }
    #endif
  }

    //lower_grid_vec(S->bcon, S->bcov, G, kstart, kstop, jstart, jstop, istart, istop, loc);
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
  bsq = bsq_calc(S, i, j, k);
  rho = fabs(S->P[RHO][k][j][i]);
  u = fabs(S->P[UU][k][j][i]);
  ef = rho + gam*u;
  ee = bsq + ef;
  va2 = bsq/ee;
  cs2 = gam*(gam - 1.)*u/ef;

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

/*------SOURCE TERMS CALLED DURING AN IDEAL MHD RUN------*/

// Source terms for equations of motion at grid zone
// Ignoring wind term for now
inline void get_ideal_fluid_source(struct GridGeom *G, struct FluidState *S, int i, int j, int k, double dU[NVAR])
{
  double mhd[NDIM][NDIM];

  DLOOP1 mhd_calc(S, i, j, k, mu, mhd[mu]); // TODO make an mhd_calc_vec?

  // Contract mhd stress tensor with connection
  PLOOP dU[ip] = 0.;
  DLOOP2 {
    for (int gam = 0; gam < NDIM; gam++)
      dU[UU+gam] += mhd[mu][nu]*G->conn[nu][gam][mu][j][i];
  }

  PLOOP dU[ip] *= G->gdet[CENT][j][i];
}

// Source terms for equations of motion over all physical zones
inline void get_ideal_fluid_source_vec(struct GridGeom *G, struct FluidState *S, GridPrim *dU)
{
#if WIND_TERM
  static struct FluidState *dS;
  static int firstc = 1;
  if (firstc) {dS = calloc(1, sizeof(struct FluidState)); firstc = 0;}
#endif

#pragma omp parallel for collapse(3)
  ZLOOP {
    double mhd[NDIM][NDIM];

    DLOOP1 mhd_calc(S, i, j, k, mu, mhd[mu]); // TODO make an mhd_calc_vec?

    // Contract mhd stress tensor with connection
    // TODO this is scattered memory access.  Precompute mu,nu sums
    PLOOP (*dU)[ip][k][j][i] = 0.;
    DLOOP2 {
      for (int gam = 0; gam < NDIM; gam++)
        (*dU)[UU+gam][k][j][i] += mhd[mu][nu]*G->conn[nu][gam][mu][j][i];
    }

    PLOOP (*dU)[ip][k][j][i] *= G->gdet[CENT][j][i];
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

/*------END OF SOURCE TERMS CALLED DURING AN IDEAL MHD RUN------*/

/*------SOURCE TERMS CALLED DURING AN EXTENDED MHD (EMHD) RUN------*/

// Compute source terms with time derivatives
// TODO: Load all variables
#if GRIM_TIMESTEPPER
void emhd_time_derivative_sources(struct GridGeom *G, struct FluidState *S_new, struct FluidState *S_old,
  struct FluidState *S, double dt, int loc, int i, int j, int k, double dU[NVAR]) {

  // Initializations
  double rho      = S->P[RHO][k][j][i];
  double Theta    = S->Theta[k][j][i];
  double bsq      = S->bsq[k][j][i];
  double chi_emhd = S->chi_emhd[k][j][i];
  double nu_emhd  = S->nu_emhd[k][j][i];
  double tau      = S->tau[k][j][i];

  double gdet = G->gdet[loc][j][i];

  // Compute partial derivative of ucov
  double dt_ucov[NDIM];
  DLOOP1 {
    double ucov_new = S_new->ucov[mu][k][j][i];
    double ucov_old = S_old->ucov[mu][k][j][i];

    dt_ucov[mu] = (ucov_new - ucov_old) / dt;
  }

  // Compute div of ucon (only temporal part is nonzero)
  double div_ucon = 0;
  DLOOP1 {
    double gcon_t_mu = G->gcon[loc][0][mu][j][i];

    div_ucon += gcon_t_mu * dt_ucov[mu];
  }

  // Compute q0 and delta_P0 (temporal terms)
  double Theta_new, Theta_old, dt_Theta;
  Theta_new = S_new->Theta[k][j][i];
  Theta_old = S_old->Theta[k][j][i];

  dt_Theta = (Theta_new - Theta_old) / dt;

  double q0, deltaP0;
  double bcon_t  = S->bcon[0][k][j][i];

  q0 = -rho * chi_emhd * (bcon_t / sqrt(bsq)) * dt_Theta;
  DLOOP1 {
    double ucon_t  = S->ucon[0][k][j][i];
    double bcon_mu = S->bcon[mu][k][j][i];

    q0 -= rho * chi_emhd * (bcon_mu / sqrt(bsq)) * Theta * ucon_t * dt_ucov[mu];
  }

  deltaP0 = -rho * nu_emhd * div_ucon;
  DLOOP1 {
    double bcon_mu = S->bcon[mu][k][j][i];

    deltaP0 += 3. * rho * nu_emhd * (bcon_t * bcon_mu / bsq) * dt_ucov[mu];
  }
  
  // Add the time derivative source terms (conduction and viscosity)
  // NOTE: Will have to edit this when higher order terms are considered
  dU[Q_TILDE]       += gdet * (q0 / tau);
  dU[DELTA_P_TILDE] += gdet * (deltaP0 / tau);
}

// Compute implicit source terms
void emhd_implicit_sources(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k, double dU[NVAR]){

  // Initializations
  double q_tilde       = S->P[Q_TILDE][k][j][i];
  double delta_p_tilde = S->P[DELTA_P_TILDE][k][j][i];
  double tau           = S->tau[k][j][i];

  double gdet = G->gdet[loc][j][i];

  dU[Q_TILDE]       -= gdet * (q_tilde / tau);
  dU[DELTA_P_TILDE] -= gdet * (delta_p_tilde / tau);
}

// Compute explicit source terms
void emhd_explicit_sources(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k, double dU[NVAR]) {
  
  // Ideal MHD components
  double mhd[NDIM][NDIM];

  DLOOP1 mhd_calc(S, i, j, k, mu, mhd[mu]);

  PLOOP dU[ip] = 0.; // explicit_sources should be the first source call
  DLOOP2 {
    for (int gam = 0; gam < NDIM; gam++)
      dU[UU+gam] += mhd[mu][nu]*G->conn[nu][gam][mu][j][i];
  }
  FLOOP dU[ip] *= G->gdet[CENT][j][i];

  // Extended MHD components

  // Initializations

  double rho      = S->P[RHO][k][j][i];
  double Theta    = S->Theta[k][j][i];
  double bsq      = S->bsq[k][j][i];
  double chi_emhd = S->chi_emhd[k][j][i];
  double nu_emhd  = S->nu_emhd[k][j][i];
  double tau      = S->tau[k][j][i];

  double gdet = G->gdet[loc][j][i];

  double grad_ucov[NDIM][NDIM], grad_Theta[NDIM];

  // Compute gradient of ucov and Theta
  gradient_calc(G, S, loc, i, j, k, grad_ucov, grad_Theta);

  // Compute div of ucon (all terms but the time-derivative ones are nonzero)
  double div_ucon = 0;
  DLOOP2 {
    double gcon_mu_nu = G->gcon[loc][mu][nu][j][i];

    div_ucon += gcon_mu_nu * grad_ucov[mu][nu];
  }

  // Compute q0 and deltaP0 (everything but the time-derivative terms)
  double q0, deltaP0;

  DLOOP1 {
    double bcon_mu = S->bcon[mu][k][j][i];

    q0 = -rho * chi_emhd * (bcon_mu / sqrt(bsq)) * grad_Theta[mu];
  }

  DLOOP2 {
    double bcon_mu = S->bcon[mu][k][j][i];
    double ucon_nu = S->ucon[nu][k][j][i];

    q0 -= rho * chi_emhd * (bcon_mu / sqrt(bsq)) * Theta * ucon_nu * grad_ucov[nu][mu];
  }

  deltaP0 = -rho * nu_emhd * div_ucon;
  DLOOP2  {
    double bcon_mu = S->bcon[mu][k][j][i];
    double bcon_nu = S->bcon[nu][k][j][i];
    
    deltaP0 += 3. * rho * nu_emhd * (bcon_mu * bcon_nu / bsq) * grad_ucov[mu][nu];
  }

  // Add explicit source terms (conduction and viscosity)
  // NOTE: Will have to edit this when higher order terms are considered
  dU[Q_TILDE]       += gdet * (q0 / tau);
  dU[DELTA_P_TILDE] += gdet * (deltaP0) / tau;
}

// Compute gradient of four velocities and temperature
// Called by emhd_explicit_sources
void gradient_calc(struct GridGeom *G, struct FluidState *S, int loc, int i, int j, int k,
  double grad_ucov[NDIM][NDIM], double grad_Theta[NDIM]) {

  // Compute gradient of ucov
  DLOOP1 {
    grad_ucov[0][mu] = 0;

    slope_calc_four_vec(S->ucov, mu, 1, i, j, k, grad_ucov[1][mu]);

    slope_calc_four_vec(S->ucov, mu, 2, i, j, k, grad_ucov[2][mu]);

    slope_calc_four_vec(S->ucov, mu, 3, i, j, k, grad_ucov[3][mu]);
  }

  DLOOP2 {
    for (int gam = 0; gam < NDIM; gam++)
      grad_ucov[mu][nu] -= G->conn[gam][mu][nu][j][i] * S->ucov[gam][k][j][i];
  }

  // Compute temperature gradient
  // Time derivative component computed in emhd_time_derivative_sources
  grad_Theta[0] = 0;

  slope_calc_scalar(S->Theta, 1, i, j, k, grad_Theta[1]);

  slope_calc_scalar(S->Theta, 2, i, j, k, grad_Theta[2]);

  slope_calc_scalar(S->Theta, 3, i, j, k, grad_Theta[3]);
}
#endif

/*------END OF SOURCE TERMS CALLED DURING AN EXTENDED MHD (EMHD) RUN------*/

// Returns b.b (twice magnetic pressure)
inline double bsq_calc(struct FluidState *S, int i, int j, int k)
{

  double bsq = 0.;
  DLOOP1 {
    bsq += S->bcon[mu][k][j][i]*S->bcov[mu][k][j][i];
  }

  return bsq;
}