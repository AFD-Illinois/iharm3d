/******************************************************************************
 *                                                                            *
 * COORD.C                                                                    *
 *                                                                            *
 * SIMULATION VOLUME COORDINATES                                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

/*      ASSUMING X3 SYMMETRY IN COORDINATES AND SPACETIME
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

#define POLYTH (0)
#if POLYTH
double poly_norm, poly_xt, poly_alpha; 
#endif

void coord(int i, int j, int k, int loc, double *X)
{
  i += global_start[0];
  j += global_start[1];
  k += global_start[2];
  if (loc == FACE1) {
    X[1] = startx[1] + (i - NG)*dx[1];
    X[2] = startx[2] + (j + 0.5 - NG)*dx[2];
    X[3] = startx[3] + (k + 0.5 - NG)*dx[3];
  } else if (loc == FACE2) {
    X[1] = startx[1] + (i + 0.5 - NG)*dx[1];
    X[2] = startx[2] + (j - NG) * dx[2];
    X[3] = startx[3] + (k + 0.5 - NG)*dx[3];
  } else if (loc == CENT) {
    X[1] = startx[1] + (i + 0.5 - NG)*dx[1];
    X[2] = startx[2] + (j + 0.5 - NG)*dx[2];
    X[3] = startx[3] + (k + 0.5 - NG)*dx[3];
  } else {
    X[1] = startx[1] + (i - NG)*dx[1];
    X[2] = startx[2] + (j - NG)*dx[2];
    X[3] = startx[3] + (k - NG)*dx[3];
  }
}

double r_of_x(double x)
{
  return exp(x);
}

double dr_dx(double x)
{
  return r_of_x(x);
}

double th_of_x(double x)
{
#if POLYTH
  double y = 2*x - 1.;
  return poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 
         0.5*M_PI;
#else
  return M_PI*x + ((1. - hslope)/2.)*sin(2.*M_PI*x);
#endif
}

double dth_dx(double x)
{
#if POLYTH
  double y = 2*x - 1.;
  return 2.*poly_norm*(1. + pow(y/poly_xt,poly_alpha));
#else
  return M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*x);
#endif
}

// Boyer-Lindquist coordinate of point X
void bl_coord(double *X, double *r, double *th)
{
  *r = r_of_x(X[1]);
  *th = th_of_x(X[2]);

  // Avoid singularity at polar axis
#if(COORDSINGFIX)
  if (fabs(*th) < SINGSMALL) {
    if ((*th) >= 0)
      *th = SINGSMALL;
    if ((*th) < 0)
      *th = -SINGSMALL;
  }
  if (fabs(M_PI - (*th)) < SINGSMALL) {
    if ((*th) >= M_PI)
      *th = M_PI + SINGSMALL;
    if ((*th) < M_PI)
      *th = M_PI - SINGSMALL;
  }
#endif

  return;
}

// Insert metric here
void gcov_func(double *X, double gcov[NDIM][NDIM])
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
/*
void get_gcon(double X[NDIM], double gcon[NDIM][NDIM]) {
  memset(gcon, 0, NDIM*NDIM*sizeof(double));
  
  #if METRIC == MINKOWSKI
  gcov[0][0] = -1.;
  for (int j = 1; j < NDIM; j++) {
    gcov[j][j] = 1.;
  }
  #elif METRIC == MKS
  double r, th;
  double tfac, rfac, hfac, pfac;

  bl_coord(X, &r, &th);

  cth = cos(th);
  sth = sin(th);

  s2 = sth*sth;
  rho2 = r*r + a*a*cth*cth;

  tfac = 1.;
  rfac = 1./dr_dx(X[1]);
  hfac = 1./dth_dx(X[2]);
  pfac = 1.;

  gcon[0][0] = rho2/(2.*r - rho2)*tfac*tfac;
  gcon[0][1] = -2.*r*rho2/(4.*r*r - rho2*rho2);


  #endif // METRIC 
}*/

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
  #if RADIATION
  startx_rad[1] = startx[1];
  startx_rad[2] = startx[2];
  startx_rad[3] = startx[3];
  stopx_rad[1] = startx_rad[1] + N1TOT*dx[1];
  stopx_rad[2] = startx_rad[2] + N1TOT*dx[2];
  stopx_rad[3] = startx_rad[3] + N1TOT*dx[3];
  #endif
  #elif METRIC == MKS
  // Set Rin such that we have 5 zones completely inside the event horizon
  Rin = exp((N1TOT*log(Rhor)/5.5 - log(Rout))/(1. + N1TOT/5.5));

  startx[1] = log(Rin);
  startx[2] = 0.;
  startx[3] = 0.;

  dx[1] = log(Rout/Rin)/N1TOT;
  dx[2] = 1./N2TOT;
  dx[3] = 2.*M_PI/N3TOT;

  #if RADIATION
  startx_rad[1] = log(Rhor);
  startx_rad[2] = startx[2];
  startx_rad[3] = startx[3];
  stopx_rad[1] = log(Rout_rad);
  stopx_rad[2] = startx[2] + N2TOT*dx[2];
  stopx_rad[3] = startx[3] + N3TOT*dx[3];
  #endif

  #if POLYTH
  poly_alpha = 8.;
  poly_xt = 0.75;
  poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*
                           1./pow(poly_xt, poly_alpha));
  #endif
  #endif // METRIC
}

void set_grid(struct GridGeom *G)
{
  //double X[NDIM];

  // Set up boundaries, steps in coordinate grid
  set_points();
  dV = dx[1]*dx[2]*dx[3];

  //for (int mu = 0; mu < NDIM; mu++) X[mu] = 0.;

  #pragma omp parallel for collapse(2)
  JSLOOP(-NG, N2 - 1 + NG) {
    ISLOOP(-NG, N1 - 1 + NG) {
      set_grid_loc(G, i, j, 0, CENT);
      set_grid_loc(G, i, j, 0, CORN);
      set_grid_loc(G, i, j, 0, FACE1);
      set_grid_loc(G, i, j, 0, FACE2);
      
      // Connection only needed at zone center
      conn_func(G, i, j, 0);
    }
  }
}

void set_grid_loc(struct GridGeom *G, int i, int j, int k, int loc)
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
  ZSLOOP(-NG, N3-1 + NG, -NG, N2-1 + NG, -NG, N1-1 + NG) {
    pflag[k][j][i] = 0;
    fail_save[k][j][i] = 0;
  }
}

