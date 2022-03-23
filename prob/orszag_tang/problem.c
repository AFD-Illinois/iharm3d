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

void init(struct GridGeom *G, struct FluidState *S)
{
  double X[NDIM];

  set_grid(G);
  LOG("Set grid");

  ZLOOP {

    coord(i, j, k, CENT, X);

    S->P[RHO][k][j][i] = 25./9.;
    S->P[UU][k][j][i]  = 5./3./(gam - 1.);
    S->P[U1][k][j][i]  = -sin(X[2] + phase);
    S->P[U2][k][j][i]  = sin(X[1] + phase);
    S->P[U3][k][j][i]  = 0.;
    S->P[B1][k][j][i]  = -sin(X[2] + phase);
    S->P[B2][k][j][i]  = sin(2.*(X[1] + phase));
    S->P[B3][k][j][i]  = 0.;

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