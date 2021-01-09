/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR ORSZAG-TANG VORTEX                                  *
 *                                                                            *
 ******************************************************************************/

/**
  relativistic version of the Orszag-Tang vortex
  Orszag & Tang 1979, JFM 90, 129-143.
  original OT problem was incompressible
  this is based on compressible version given
  in Toth 2000, JCP 161, 605.

  in the limit tscale -> 0 the problem is identical
  to the nonrelativistic problem; as tscale increases
  the problem becomes increasingly relativistic 

  grid is fixed in coordinates.c, runs from
  -\pi < x1 <= \pi
  -\pi < x2 <= \pi
*/

#include "decs.h"
#include <complex.h>
#include "hdf5_utils.h"

static double tscale, phase;

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
  /* physical parameters */
  gam = 5./3.;

  /* set up grid functions */
  set_grid(G);

  /* end of simulation */
  t = 0.;
  tf = 1.5*M_PI/tscale;
  dt = 1.e-4;
  cour = 0.9;

  /* output frequencies */
  DTd = tf/10.;		/* dumping frequency, in units of M */
  DTl = tf/100.;		/* logfile frequency, in units of M */
  DTr = 512;		/* restart file frequ., in timesteps */

  /* start diagnostic counters */
  nstep = 0;
  dump_cnt = 0;

  /* relativistic version of the Orszag-Tang vortex 
    Orszag & Tang 1979, JFM 90, 129-143. 
    original OT problem was incompressible 
    this is based on compressible version given
    in Toth 2000, JCP 161, 605.

    in the limit tscale -> 0 the problem is identical
    to the nonrelativistic problem; as tscale increases
    the problem becomes increasingly relativistic */
	ZLOOP {
    double X[NDIM];
		coord(i, j, k, CENT, X);

    S->P[U1][k][j][i] = -sin(X[2] + phase);
    S->P[U2][k][j][i] =  sin(X[1] + phase) ;
    S->P[U3][k][j][i] = 0.;

    S->P[B1][k][j][i] = -sin(X[2] + phase);
    S->P[B2][k][j][i] =  sin(2.*(X[1] + phase));
    S->P[B3][k][j][i] = 0. ;

    S->P[RHO][k][j][i] = 25./9.;
    S->P[UU][k][j][i] = 5./(3.*(gam - 1.));

    S->P[U1][k][j][i] *= tscale ;
    S->P[U2][k][j][i] *= tscale ;
    S->P[U3][k][j][i] *= tscale ;
    S->P[B1][k][j][i] *= tscale ;
    S->P[B2][k][j][i] *= tscale ;
    S->P[B3][k][j][i] *= tscale ;
    S->P[UU][k][j][i] *= tscale*tscale ;
	}

  /* enforce boundary conditions */
  fixup(G, S, CENT);
  set_bounds(G, S);
}

