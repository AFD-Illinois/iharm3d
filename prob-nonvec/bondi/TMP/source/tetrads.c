/******************************************************************************  
 *                                                                            *  
 * TETRADS.C                                                                  *  
 *                                                                            *  
 * CONSTRUCT AND APPLY TRANSFORMATION MATRICES BETWEEN FLUID AND COORDINATE   *
 * FRAMES                                                                     *  
 *                                                                            *  
 ******************************************************************************/

#include "decs.h"

#if RADIATION
void normalize(double *Vcon, double Gcov[4][4]);
void project_out(double *Vcona, double *Vconb, double Gcov[4][4]);

void coord_to_tetrad(double Ecov[NDIM][NDIM], double Kcoord[NDIM], 
  double Ktetrad[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++) {
    Ktetrad[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      Ktetrad[mu] += Ecov[mu][nu]*Kcoord[nu];
    }
  }
}

void tetrad_to_coord(double Econ[NDIM][NDIM], double Ktetrad[NDIM],
  double Kcoord[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++) {
    Kcoord[mu] = 0.;
    for (int nu = 0; nu < NDIM; nu++) {
      Kcoord[mu] += Econ[nu][mu]*Ktetrad[nu];
    }
  }
}

#define SMALL_VECTOR (1.e-30)
void make_tetrad(int i, int j, int k, double Ucon[NDIM], double trial[NDIM], 
  double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM])
{
  // Time component parallel to u^{\mu}
  for (int mu = 0; mu < NDIM; mu++) {
    Econ[0][mu] = Ucon[mu];
  }
  normalize(Econ[0], Gcov);

  // Normalize trial vector
  double norm = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      norm += trial[mu]*trial[nu]*Gcov[mu][nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) {
    trial[mu] /= norm;
  }
  
  // If trial vector bad, default to X1 direction
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      norm += trial[mu]*trial[nu]*Gcov[mu][nu];
    }
  }
  if (norm <= SMALL_VECTOR) {
    for (int mu = 0; mu < NDIM; mu++) {
      trial[mu] = delta(mu, 1);
    }
  }

  for (int mu = 0; mu < NDIM; mu++) {
    Econ[1][mu] = trial[mu];
  }
  project_out(Econ[1], Econ[0], Gcov);
  normalize(Econ[1], Gcov);

  // Use X2 direction
  for (int mu = 0; mu < NDIM; mu++) {
    Econ[2][mu] = delta(mu, 2);
  }
  project_out(Econ[2], Econ[0], Gcov);
  project_out(Econ[2], Econ[1], Gcov);
  normalize(Econ[2], Gcov);

  // Use X3 direction
  for (int mu = 0; mu < NDIM; mu++) {
    Econ[3][mu] = delta(mu, 3);
  }
  project_out(Econ[3], Econ[0], Gcov);
  project_out(Econ[3], Econ[1], Gcov);
  project_out(Econ[3], Econ[2], Gcov);
  normalize(Econ[3], Gcov);

  // Make covariant version
  for (int mu = 0; mu < NDIM; mu++) {
    lower(Econ[mu], &(ggeom[i][j][CENT]), Ecov[mu]);
  }
  for (int mu = 0; mu < NDIM; mu++) {
    Ecov[0][mu] *= -1.;
  }
}

void normalize(double *Vcon, double Gcov[NDIM][NDIM])
{
  double norm = 0.;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      norm += Vcon[mu]*Vcon[nu]*Gcov[mu][nu];
    }
  }

  norm = sqrt(fabs(norm));
  for (int mu = 0; mu < NDIM; mu++) {
    Vcon[mu] /= norm;
  }
}

void project_out(double *Vcona, double *Vconb, double Gcov[NDIM][NDIM])
{
  double Vconb_sq = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      Vconb_sq += Vconb[mu]*Vconb[nu]*Gcov[mu][nu];
    }
  }

  double adotb = 0.;
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      adotb += Vcona[mu]*Vconb[nu]*Gcov[mu][nu];
    }
  }

  for (int mu = 0; mu < NDIM; mu++) {
    Vcona[mu] -= Vconb[mu]*adotb/Vconb_sq;
  }
}

// Enforce K.K = 0
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM])
{
  double A = Gcov[0][0];
  double B = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    B += 2.*Gcov[mu][0]*K[mu];
  }
  double C = 0.;
  for (int mu = 1; mu < NDIM; mu++) {
    for (int nu = 1; nu < NDIM; nu++) {
      C += Gcov[mu][nu]*K[mu]*K[nu];
    }
  }

  K[0] = (-B - sqrt(fabs(B*B - 4.*A*C)))/(2.*A);
}
#undef SMALL_VECTOR
#endif // RADIATION
