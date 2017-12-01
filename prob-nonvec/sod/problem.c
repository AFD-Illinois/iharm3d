/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR SOD SHOCKTUBE                                       *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static double tscale;
void set_problem_params()
{
  set_param("tscale", &tscale);
}

void init_prob()
{
	double X[NDIM];

  // Make problem nonrelativistic
  double tscale = 1.e-2;
  tf /= tscale;
  dt /= tscale;
  DTd /= tscale;
  DTl /= tscale;

  ZLOOP {
    coord(i, j, k, CENT, X);

    P[i][j][k][RHO] = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.125;

    double pgas = (X[1] < 0.5 || X[1] > 1.5) ? 1.0 : 0.1;
    P[i][j][k][UU] = pgas/(gam - 1.);

    P[i][j][k][U1] = 0.;
    P[i][j][k][U2] = 0.;
    P[i][j][k][U3] = 0.;
    P[i][j][k][B1] = 0.;
    P[i][j][k][B2] = 0.;
    P[i][j][k][B3] = 0.;
  } // ZLOOP

  // Rescale to make problem nonrelativistic
  ZLOOP {
    P[i][j][k][UU] *= tscale*tscale;
    P[i][j][k][U1] *= tscale;
    P[i][j][k][U2] *= tscale;
    P[i][j][k][U3] *= tscale;
    P[i][j][k][B1] *= tscale;
    P[i][j][k][B2] *= tscale;
    P[i][j][k][B3] *= tscale;
  }
}

