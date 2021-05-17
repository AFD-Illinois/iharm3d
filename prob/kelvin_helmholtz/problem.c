
/**
 * Kelvin-Helmholtz instability problem
 * 
 * follows initial conditions from Lecoanet et al. 2015,
 * MNRAS 455, 4274.
 */

#include "decs.h"
#include "hdf5_utils.h"

static double tscale;

void set_problem_params()
{
  // TODO parameters
  set_param("tscale", &tscale);
}

void save_problem_data(hid_t string_type){
        hdf5_write_single_val("kelvin_helmholtz", "PROB", string_type);
        hdf5_write_single_val(&tscale, "tscale", H5T_IEEE_F64LE);
}

void init(struct GridGeom *G, struct FluidState *S)
{
  // Set up grid functions
  set_grid(G);

  LOG("Set grid");

  // Set tf based on the scaling
  tf = 6.0 / tscale;
  DTd = tf/10.;

  // follows notation of Lecoanet et al. eq. 8 et seq.
  double rho0 = 1.;
  double Drho = 0.1;
  double P0 = 10.;
  double uflow = 1.;
  double a = 0.05;
  double sigma = 0.2;
  double A = 0.01;
  double z1 = 0.5;
  double z2 = 1.5;

  ZLOOP {
    double X[NDIM];
    coord(i, j, k, CENT, X);

    /* Lecoanet's x <-> x1; z <-> x2 */
    double x = X[1];
    double z = X[2];

    S->P[RHO][k][j][i] = 1. + (Drho / rho0) * 0.5 * (tanh((z - z1) / a) - tanh((z - z2) / a));
    S->P[UU][k][j][i] = P0 / (gam - 1.);

    S->P[U1][k][j][i] = uflow * (tanh((z - z1) / a) -
                                 tanh((z - z2) / a) - 1.);
    S->P[U2][k][j][i] = A * sin(2. * M_PI * x) * (exp(-(z - z1) * (z - z1) / (sigma * sigma)) +
                        exp(-(z - z2) * (z - z2) / (sigma * sigma)));
    S->P[U3][k][j][i] = 0.;

    S->P[B1][k][j][i] = 0.;
    S->P[B2][k][j][i] = 0.;
    S->P[B3][k][j][i] = 0.;

    /* rescale */
    S->P[U1][k][j][i] *= tscale;
    S->P[U2][k][j][i] *= tscale;
    S->P[U3][k][j][i] *= tscale;
    S->P[B1][k][j][i] *= tscale;
    S->P[B2][k][j][i] *= tscale;
    S->P[B3][k][j][i] *= tscale;
    S->P[UU][k][j][i] *= tscale * tscale;
  }

  // Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);
}
