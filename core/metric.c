/******************************************************************************
 *                                                                            *
 * METRIC.C                                                                   *
 *                                                                            *
 * HELPER FUNCTIONS FOR METRIC TENSORS                                        *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

double gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM])
{
  double gdet = invert(&gcov[0][0],&gcon[0][0]);
  return sqrt(fabs(gdet));
}

inline void get_gcov(struct GridGeom *G, int i, int j, int loc, double gcov[NDIM][NDIM]) {
  DLOOP2 gcov[mu][nu] = G->gcov[loc][mu][nu][j][i];
}

inline void get_gcon(struct GridGeom *G, int i, int j, int loc, double gcon[NDIM][NDIM])
{
  DLOOP2 gcon[mu][nu] = G->gcon[loc][mu][nu][j][i];
}

// Calculate connection coefficient \Gamma^{i}_{j,k} = conn[..][i][j][k]
void conn_func(struct GridGeom *G, int i, int j, int k)
{
  double tmp[NDIM][NDIM][NDIM];
  double X[NDIM], Xh[NDIM], Xl[NDIM];
  double gh[NDIM][NDIM];
  double gl[NDIM][NDIM];
  coord(i, j, k, CENT, X);

  for (int mu = 0; mu < NDIM; mu++) {
    for (int kap = 0; kap < NDIM; kap++) {
      Xh[kap] = X[kap];
    }
    for (int kap = 0; kap < NDIM; kap++) {
      Xl[kap] = X[kap];
    }
    Xh[mu] += DELTA;
    Xl[mu] -= DELTA;
    gcov_func(Xh, gh);
    gcov_func(Xl, gl);

    for (int lam = 0; lam < NDIM; lam++) {
      for (int nu = 0; nu < NDIM; nu++) {
        G->conn[lam][nu][mu][j][i] = (gh[lam][nu] - gl[lam][nu])/(Xh[mu] - 
                                                                  Xl[mu]);
      }
    }
  }

  // Rearrange to find \Gamma_{lam nu mu}
  for (int lam = 0; lam < NDIM; lam++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int mu = 0; mu < NDIM; mu++) {
        tmp[lam][nu][mu] = 0.5 * (G->conn[nu][lam][mu][j][i] + 
                                  G->conn[mu][lam][nu][j][i] - 
                                  G->conn[mu][nu][lam][j][i]);
      }
    }
  }

  // now mu nu kap

  // Raise index to get \Gamma^lam_{nu mu}
  for (int lam = 0; lam < NDIM; lam++) {
    for (int nu = 0; nu < NDIM; nu++) {
      for (int mu = 0; mu < NDIM; mu++) {
        G->conn[lam][nu][mu][j][i] = 0.;
        for (int kap = 0; kap < NDIM; kap++)
          G->conn[lam][nu][mu][j][i] += G->gcon[CENT][lam][kap][j][i]*
                                        tmp[kap][nu][mu];
      }
    }
  }
}

// Lower a contravariant rank-1 tensor to a covariant one
inline void lower_grid(GridVector vcon, GridVector vcov, struct GridGeom *G, int i,
  int j, int k, int loc)
{
  for (int mu = 0; mu < NDIM; mu++) {
    vcov[mu][k][j][i] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      vcov[mu][k][j][i] += G->gcov[loc][mu][nu][j][i]*vcon[nu][k][j][i];
    }
  }
}

void raise_grid(GridVector vcov, GridVector vcon, struct GridGeom *G, int i, 
  int j, int k, int loc)
{
  for (int mu = 0; mu < NDIM; mu++) {
    vcon[mu][k][j][i] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      vcon[mu][k][j][i] += G->gcon[loc][mu][nu][j][i]*vcov[nu][k][j][i];
    }
  }
}

// TODO revise this out of the following: Fcon_calc, fixup1zone, get_state, Utoprim, Wp_func
// And while you're at it revise out get_state
void lower(double ucon[NDIM], double gcov[NDIM][NDIM], double ucov[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++) {
    ucov[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      ucov[mu] += gcov[mu][nu]*ucon[nu];
    }
  }
}

// Raise a covariant rank-1 tensor to a contravariant one
void raise(double ucov[NDIM], double gcon[NDIM][NDIM], double ucon[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++) {
    ucon[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      ucon[mu] += gcon[mu][nu]*ucov[nu];
    }
  }
}

double dot_grid(GridVector vcon, GridVector vcov, int i, int j, int k)
{
  double dot = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    dot += vcon[mu][k][j][i]*vcov[mu][k][j][i];
  }
  return dot;
}

double dot(double vcon[NDIM], double vcov[NDIM])
{
  double dot = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    dot += vcon[mu]*vcov[mu];
  }
  return dot;
}

// TODO replace all the following w/GSL

double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2)
{
  return m[4*r0+c0]*(m[4*r1+c1]*m[4*r2+c2] - m[4*r2+c1]*m[4*r1+c2]) -
         m[4*r0+c1]*(m[4*r1+c0]*m[4*r2+c2] - m[4*r2+c0]*m[4*r1+c2]) +
         m[4*r0+c2]*(m[4*r1+c0]*m[4*r2+c1] - m[4*r2+c0]*m[4*r1+c1]);
}

void adjoint(double m[16], double adjOut[16])
{
  adjOut[ 0] =  MINOR(m,1,2,3,1,2,3);
  adjOut[ 1] = -MINOR(m,0,2,3,1,2,3);
  adjOut[ 2] =  MINOR(m,0,1,3,1,2,3);
  adjOut[ 3] = -MINOR(m,0,1,2,1,2,3);

  adjOut[ 4] = -MINOR(m,1,2,3,0,2,3);
  adjOut[ 5] =  MINOR(m,0,2,3,0,2,3);
  adjOut[ 6] = -MINOR(m,0,1,3,0,2,3);
  adjOut[ 7] =  MINOR(m,0,1,2,0,2,3);

  adjOut[ 8] =  MINOR(m,1,2,3,0,1,3);
  adjOut[ 9] = -MINOR(m,0,2,3,0,1,3);
  adjOut[10] =  MINOR(m,0,1,3,0,1,3);
  adjOut[11] = -MINOR(m,0,1,2,0,1,3);

  adjOut[12] = -MINOR(m,1,2,3,0,1,2);
  adjOut[13] =  MINOR(m,0,2,3,0,1,2);
  adjOut[14] = -MINOR(m,0,1,3,0,1,2);
  adjOut[15] =  MINOR(m,0,1,2,0,1,2);
}

double determinant(double m[16])
{
  return m[0]*MINOR(m,1,2,3,1,2,3) -
         m[1]*MINOR(m,1,2,3,0,2,3) +
         m[2]*MINOR(m,1,2,3,0,1,3) -
         m[3]*MINOR(m,1,2,3,0,1,2);
}

double invert(double *m, double *invOut)
{
  adjoint(m, invOut);

  double det = determinant(m);
  double inv_det = 1. / det;
  for (int i = 0; i < 16; ++i) {
    invOut[i] = invOut[i]*inv_det;
  }

  return det;
}

