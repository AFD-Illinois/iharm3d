/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR FISHBONE-MONCRIEF TORUS                             *
 *                                                                            *
 ******************************************************************************/

#include "bl_coord.h"
#include "decs.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_rng *rng;

// Local declarations
double lfish_calc(double rmax);

static double beta;
static double rin, rmax;
void set_problem_params()
{
  set_param("rin", &rin);
  set_param("rmax", &rmax);

  set_param("beta", &beta);
}

void init(struct GridGeom *G, struct FluidState *S)
{
  double r, th, sth, cth;
  double ur, uh, up, u, rho;
  double X[NDIM];

  // Initialize RNG
  rng = gsl_rng_alloc(gsl_rng_mt19937);
  // Could set this via time, some increasing int, whatever
  // TODO include seed in dump? Likely relies on details of problem.c so not so useful
  gsl_rng_set(rng, mpi_myrank());

  // Disk interior
  double l, lnh, expm2chi, up1;
  double DD, AA, SS, thin, sthin, cthin, DDin, AAin, SSin;
  double kappa, hm1;

  // Magnetic field
  static double A[N1 + 2*NG][N2 + 2*NG];
  double rho_av, rhomax, umax, bsq_ij, bsq_max, norm, q,
      beta_act;

  // Adiabatic index
#if !ELECTRONS
  tp_over_te = 3.;
#endif

  // Fishbone-Moncrief parameters
  l = lfish_calc(rmax);
  kappa = 1.e-3;

  // Numerical parameters

  R0 = 0.0;
  Rhor = (1. + sqrt(1. - a*a));
  //double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  //double z2 = sqrt(3*a*a + z1*z1);
  //Risco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));

  zero_arrays();
  set_grid(G);

  if (mpi_io_proc()) printf("grid set\n");

  // Diagnostic counters
  t = 0.;
  dump_cnt = 0;

  rhomax = 0.;
  umax = 0.;
  ZSLOOP(-1, N3, -1, N2, -1, N1) {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);

    sth = sin(th);
    cth = cos(th);

    // Calculate lnh
    DD = r * r - 2. * r + a * a;
    AA = (r * r + a * a) * (r * r + a * a) -
        DD * a * a * sth * sth;
    SS = r * r + a * a * cth * cth;

    thin = M_PI / 2.;

    sthin = sin(thin);
    cthin = cos(thin);
    DDin = rin * rin - 2. * rin + a * a;
    AAin = (rin * rin + a * a) * (rin * rin + a * a)
        - DDin * a * a * sthin * sthin;
    SSin = rin * rin + a * a * cthin * cthin;

    if (r >= rin) {
      lnh =
          0.5 *
          log((1. +
         sqrt(1. +
              4. * (l * l * SS * SS) * DD / (AA *
                     sth *
                     AA *
                     sth)))
        / (SS * DD / AA))
          - 0.5 * sqrt(1. +
           4. * (l * l * SS * SS) * DD /
           (AA * AA * sth * sth))
          - 2. * a * r * l / AA -
          (0.5 *
           log((1. +
          sqrt(1. +
               4. * (l * l * SSin * SSin) * DDin /
               (AAin * AAin * sthin * sthin))) /
         (SSin * DDin / AAin))
           - 0.5 * sqrt(1. +
            4. * (l * l * SSin * SSin) *
            DDin / (AAin * AAin * sthin *
              sthin))
           - 2. * a * rin * l / AAin);
    } else {
      lnh = 1.;
    }

    // regions outside torus
    if (lnh < 0. || r < rin) {
      // Nominal values; real value set by fixup
      rho = 1.e-7 * RHOMIN;
      u = 1.e-7 * UUMIN;

      ur = 0.;
      uh = 0.;
      up = 0.;

      S->P[RHO][k][j][i] = rho;
      S->P[UU][k][j][i] = u;
      S->P[U1][k][j][i] = ur;
      S->P[U2][k][j][i] = uh;
      S->P[U3][k][j][i] = up;
    }
    /* region inside magnetized torus; u^i is calculated in
     * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
     * so it needs to be transformed at the end */
    else {
      hm1 = exp(lnh) - 1.;
      rho = pow(hm1 * (gam - 1.) / (kappa * gam),
          1. / (gam - 1.));
      u = kappa * pow(rho, gam) / (gam - 1.);
      ur = 0.;
      uh = 0.;

      // Calculate u^phi
      expm2chi = SS * SS * DD / (AA * AA * sth * sth);
      up1 =
          sqrt((-1. +
          sqrt(1. + 4. * l * l * expm2chi)) / 2.);
      up = 2. * a * r * sqrt(1. +
                 up1 * up1) / sqrt(AA * SS *
                 DD) +
          sqrt(SS / AA) * up1 / sth;


      S->P[RHO][k][j][i] = rho;
      if (rho > rhomax) rhomax = rho;
      S->P[UU][k][j][i] = u * (1. + 4.e-2 * (gsl_rng_uniform(rng) - 0.5));
      if (u > umax && r > rin) umax = u;
      S->P[U1][k][j][i] = ur;
      S->P[U2][k][j][i] = uh;
      S->P[U3][k][j][i] = up;

      // Convert from 4-velocity to 3-velocity
      coord_transform(G, S, i, j, k);
    }

    S->P[B1][k][j][i] = 0.;
    S->P[B2][k][j][i] = 0.;
    S->P[B3][k][j][i] = 0.;
  } // ZSLOOP

  // Normalize the densities so that max(rho) = 1
  umax = mpi_max(umax);
  rhomax = mpi_max(rhomax);

  ZSLOOP(-1, N3, -1, N2, -1, N1) {
    S->P[RHO][k][j][i] /= rhomax;
    S->P[UU][k][j][i] /= rhomax;
  }
  umax /= rhomax;
  rhomax = 1.;
  LOG("Starting first Fixup");
  fixup(G, S);
  LOG("Starting first bounds");
  set_bounds(G, S);
  LOG("Finished first Fixup/bounds");

  // first find corner-centered vector potential
  ZSLOOP(0, 0, -NG, N2+NG-1, -NG, N1+NG-1) A[i][j] = 0.;
  ZSLOOP(0, 0, 0, N2, 0, N1) {
    /* vertical field version */
    /*
       coord(i,j,CORN,X) ;
       bl_coord(X,&r,&th) ;

       A[i][j] = 0.5*r*sin(th) ;
     */

    /* field-in-disk version */
    /* flux_ct */
    rho_av = 0.25*(S->P[RHO][NG][j][i] + S->P[RHO][NG][j][i-1] +
                   S->P[RHO][NG][j-1][i] + S->P[RHO][NG][j-1][i-1])
      *(1. + 0.0*(gsl_rng_uniform(rng) - 0.5));

    q = rho_av/rhomax - 0.2;
    if (q > 0.) A[i][j] = q;
  }

  // Differentiate to find cell-centered B, and begin normalization
  bsq_max = 0.;
  ZLOOP {
    //geom = get_geometry(i, j, CENT) ;

    // Flux-ct
    S->P[B1][k][j][i] = -(A[i][j] - A[i][j + 1]
        + A[i + 1][j] - A[i + 1][j + 1]) /
        (2. * dx[2] * G->gdet[CENT][j][i]);
    S->P[B2][k][j][i] = (A[i][j] + A[i][j + 1]
             - A[i + 1][j] - A[i + 1][j + 1]) /
             (2. * dx[1] * G->gdet[CENT][j][i]);

    S->P[B3][k][j][i] = 0.;

    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i, j, k);
    if (bsq_ij > bsq_max) bsq_max = bsq_ij;
  }
  bsq_max = mpi_max(bsq_max);

  // Normalize to set field strength
  beta_act = (gam - 1.) * umax / (0.5 * bsq_max);
  norm = sqrt(beta_act / beta);
  bsq_max = 0.;
  ZLOOP {
    S->P[B1][k][j][i] *= norm;
    S->P[B2][k][j][i] *= norm;

    //geom = get_geometry(i, j, CENT) ;
    //bsq_ij = bsq_calc(P[i][j][k], geom);
    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i, j, k);
    if (bsq_ij > bsq_max)
      bsq_max = bsq_ij;
  }
  bsq_max = mpi_max(bsq_max);
  beta_act = (gam - 1.)*umax/(0.5*bsq_max);

#if ELECTRONS
  init_electrons(G, S);
#endif

  // Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);

  LOG("Finished init()");
}

double lfish_calc(double r)
{
  return (((pow(a, 2) - 2. * a * sqrt(r) + pow(r, 2)) *
     ((-2. * a * r *
       (pow(a, 2) - 2. * a * sqrt(r) +
        pow(r,
      2))) / sqrt(2. * a * sqrt(r) + (-3. + r) * r) +
      ((a + (-2. + r) * sqrt(r)) * (pow(r, 3) + pow(a, 2) *
      (2. + r))) / sqrt(1 + (2. * a) / pow (r, 1.5) - 3. / r)))
    / (pow(r, 3) * sqrt(2. * a * sqrt(r) + (-3. + r) * r) *
       (pow(a, 2) + (-2. + r) * r))
      );
}
