/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ORSZAG TANG                                         *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include "hdf5_utils.h"

static double tscale;
static double phase;

void set_problem_params()
{
  set_param("tscale", &tscale);
  set_param("phase", &phase);
}

void save_problem_data(hid_t string_type){
        hdf5_write_single_val("orszag_tang", "PROB", string_type);
        hdf5_write_single_val(&tscale, "tscale", H5T_IEEE_F64LE);
        hdf5_write_single_val(&phase, "phase", H5T_IEEE_F64LE);
}

// Set chi, nu, tau. Problem dependent
 void set_emhd_parameters(struct GridGeom *G, struct FluidState *S, int i, int j, int k){
     
  // Initializations
  double rho = S->P[RHO][k][j][i];
  double u   = S->P[UU][k][j][i];
  double P   = (gam - 1.) * u;
  
  //sound speed
  double cs2 = (gam * P) / (rho + (gam * u));

  // set EMHD parameters based on closure relations
  double tau   = 1.0;
  S->tau[k][j][i]      = tau;
  S->chi_emhd[k][j][i] = conduction_alpha * cs2 * tau;;
  S->nu_emhd[k][j][i]  = viscosity_alpha * cs2 * tau;;
 }

void init(struct GridGeom *G, struct FluidState *S)
{
  double X[NDIM];

  set_grid(G);
  LOG("Set grid");

  ZLOOP {

    coord(i, j, k, CENT, X);

    S->P[RHO][k][j][i]           = 25./9.;
    S->P[UU][k][j][i]            = 5./3./(gam - 1.);
    S->P[U1][k][j][i]            = -sin(X[2] + phase);
    S->P[U2][k][j][i]            = sin(X[1] + phase);
    S->P[U3][k][j][i]            = 0.;
    S->P[B1][k][j][i]            = -sin(X[2] + phase);
    S->P[B2][k][j][i]            = sin(2.*(X[1] + phase));
    S->P[B3][k][j][i]            = 0.;
    S->P[Q_TILDE][k][j][i]       = 0.0;
    S->P[DELTA_P_TILDE][k][j][i] = 0.0;

    // rescale by tscale. 0(non-relativistic) < tscale < 1(fully relativistic)
    S->P[UU][k][j][i] *= tscale * tscale;
    S->P[U1][k][j][i] *= tscale;
    S->P[U2][k][j][i] *= tscale;
    S->P[U3][k][j][i] *= tscale;
    S->P[B1][k][j][i] *= tscale;
    S->P[B2][k][j][i] *= tscale;
    S->P[B3][k][j][i] *= tscale;

  } // ZLOOP

  //Enforce boundary conditions
  set_bounds(G, S);
}
