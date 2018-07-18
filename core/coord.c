/******************************************************************************
 *                                                                            *
 * COORD.C                                                                    *
 *                                                                            *
 * SIMULATION VOLUME COORDINATES                                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

/*
 *      -- given the indices i,j and location in the cell, return with
 *         the values of X1,X2 there;
 *      -- the locations are defined by :
 *          -----------------------
 *          |                     |
 *          |                     |
 *          |FACE1   CENT         |
 *          |                     |
 *          |CORN    FACE2        |
 *          ----------------------
 *
 */

inline void coord(int i, int j, int k, int loc, double *X)
{
  i += global_start[0];
  j += global_start[1];
  k += global_start[2];

  X[0] = 0; // Make sure all memory passed in is initialized
  if (loc == FACE1) {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG) * dx[3];
  } else if (loc == FACE2) {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG) * dx[3];
  } else if (loc == FACE3) {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
    X[3] = startx[3] + (k - NG) * dx[3];
  } else if (loc == CENT) {
    X[1] = startx[1] + (i + 0.5 - NG) * dx[1];
    X[2] = startx[2] + (j + 0.5 - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG) * dx[3];
  } else if (loc == CORN) {
    X[1] = startx[1] + (i - NG) * dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
    X[3] = startx[3] + (k - NG) * dx[3];
  }
#if DEBUG
  else {
    fprintf(stderr, "Invalid coordinate location!\n");
    exit(-1);
  }
#endif
}

// These are MKS specific
inline double r_of_x(double x)
{
  return exp(x); // TODO R0?
}

inline double dr_dx(double x)
{
  return r_of_x(x);
}

inline double th_of_x(double x)
{
#if POLYTH
  double y = 2*x - 1.;
  return poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) +
         0.5*M_PI;
#else
  return M_PI*x + ((1. - hslope)/2.)*sin(2.*M_PI*x);
#endif
}

inline double dth_dx(double x)
{
#if POLYTH
  double y = 2*x - 1.;
  return 2.*poly_norm*(1. + pow(y/poly_xt,poly_alpha));
#else
  return M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*x);
#endif
}

// Insert metric here
inline void gcov_func(double *X, double gcov[NDIM][NDIM])
{
  memset(gcov, 0, NDIM*NDIM*sizeof(double));

  #if METRIC == MINKOWSKI
  gcov[0][0] = -1.;
  for (int j = 1; j < NDIM; j++) {
    gcov[j][j] = 1.;
  }
  #elif METRIC == MKS
  double sth, cth, s2, rho2;
  double r, th;
  double tfac, rfac, hfac, pfac;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  // Exploit diagonal transformation matrix between KS and MKS
  tfac = 1.;
  rfac = dr_dx(X[1]);
  hfac = dth_dx(X[2]);
  pfac = 1.;

  gcov[0][0] = (-1. + 2.*r/rho2)*tfac*tfac;
  gcov[0][1] = (2. * r/rho2)*tfac*rfac;
  gcov[0][3] = (-2.*a*r*s2/rho2)*tfac*pfac;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = (1. + 2.*r/rho2)*rfac*rfac;
  gcov[1][3] = (-a*s2*(1. + 2.*r/rho2))*rfac*pfac;

  gcov[2][2] = rho2*hfac*hfac;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2))*pfac*pfac;
  #endif // METRIC
}

// Establish X coordinates
void set_points()
{
  #if METRIC == MINKOWSKI
  startx[1] = x1Min;
  startx[2] = x2Min;
  startx[3] = x3Min;
  dx[1] = (x1Max - x1Min)/N1TOT;
  dx[2] = (x2Max - x2Min)/N2TOT;
  dx[3] = (x3Max - x3Min)/N3TOT;
  #elif METRIC == MKS
  // Set Rin such that we have 5 zones completely inside the event horizon
  // TODO Error out if there are not enough zones
  Rin = exp((N1TOT*log(Rhor)/5.5 - log(Rout))/(1. + N1TOT/5.5));

  startx[1] = log(Rin);
  startx[2] = 0.;
  startx[3] = 0.;

  dx[1] = log(Rout/Rin)/N1TOT;
  dx[2] = 1./N2TOT;
  dx[3] = 2.*M_PI/N3TOT;

  #if POLYTH
  poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*
                           1./pow(poly_xt, poly_alpha));
  #endif
  #endif // METRIC
}

void set_grid(struct GridGeom *G)
{

  // Set up boundaries, steps in coordinate grid
  set_points();
  dV = dx[1]*dx[2]*dx[3];

#pragma omp parallel for collapse(2)
  JSLOOP(-NG, N2 - 1 + NG) {
    ISLOOP(-NG, N1 - 1 + NG) {
      set_grid_loc(G, i, j, 0, CENT);
      set_grid_loc(G, i, j, 0, CORN);
      set_grid_loc(G, i, j, 0, FACE1);
      set_grid_loc(G, i, j, 0, FACE2);
      set_grid_loc(G, i, j, 0, FACE3);

      // Connection only needed at zone center
      conn_func(G, i, j, 0);
    }
  }
}

inline void set_grid_loc(struct GridGeom *G, int i, int j, int k, int loc)
{
  double X[NDIM];
  double gcov[NDIM][NDIM], gcon[NDIM][NDIM];

  coord(i, j, k, loc, X);
  gcov_func(X, gcov);
  G->gdet[loc][j][i] = gcon_func(gcov, gcon);
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      G->gcov[loc][mu][nu][j][i] = gcov[mu][nu];
      G->gcon[loc][mu][nu][j][i] = gcon[mu][nu];
    }
  }
  G->lapse[loc][j][i] = 1./sqrt(-G->gcon[loc][0][0][j][i]);
}

void zero_arrays()
{
  ZLOOPALL {
    pflag[k][j][i] = 0;
    fail_save[k][j][i] = 0;
  }
}
