
/**
 * Initial conditions for Komissarov cylindrical
 * explosion problem.
 * 
 * see Komissarov 1999, MNRAS 303, 343
 * explosion problem is described on p. 360.
 */

#include "decs.h"
#include "hdf5_utils.h"

static double Bx;

double fr(double R,double min,double max);

void set_problem_params()
{
  set_param("Bx", &Bx);
}

void save_problem_data(hid_t string_type){
        hdf5_write_single_val("explosion", "PROB", string_type);
        hdf5_write_single_val(&Bx, "Bx", H5T_IEEE_F64LE);
}

void init(struct GridGeom *G, struct FluidState *S)
{
  double X[NDIM];

  /* physical parameters */
  gam = 4./3.;

  /* set up grid functions */
  set_grid(G);

  LOG("Set grid");

  ZLOOP {
    double X[NDIM];
    coord(i, j, k, CENT, X);
    double r = sqrt(X[1]*X[1] + X[2]*X[2]) ;

    S->P[RHO][k][j][i] = fr(r,1.e-4,1.e-2);
    S->P[UU][k][j][i] = fr(r,3.e-5,1.)/(gam - 1.);
    S->P[U1][k][j][i] = 0;
    S->P[U2][k][j][i] = 0;
    S->P[U3][k][j][i] = 0;
    S->P[B1][k][j][i] = Bx;
    S->P[B2][k][j][i] = 0;
    S->P[B3][k][j][i] = 0;
  }

  // Enforce boundary conditions
  fixup(G, S);
  set_bounds(G, S);
}

/* this is the initial density, internal energy profile
   described on Komissarov p. 360 */
double fr(double R,double min,double max)
{
  double Rout,Rin,dR ;

  Rin = 0.8 ;
  Rout = 1.0 ;
  dR = (Rout - Rin) ;

  if(R > Rout) return(min) ;
  if(R < Rin) return(max) ;
  else {
    return( exp(
            log(max)*(Rout - R)/dR +
            log(min)*(R - Rin)/dR)) ;
  }
}

