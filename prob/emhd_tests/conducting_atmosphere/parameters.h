/******************************************************************************
 *                                                                            *
 * PARAMETERS.H                                                               *
 *                                                                            *
 * PROBLEM-SPECIFIC CHOICES                                                   *
 *                                                                            *
 ******************************************************************************/

/* GLOBAL RESOLUTION */
#define N1TOT 128
#define N2TOT 128
#define N3TOT 1

/* MPI DECOMPOSITION */
/* DECOMPOSE IN N3 FIRST! Small leading array sizes for linear access */
#define N1CPU 1
#define N2CPU 1
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MKS

/*
 * FLOORS
 * Wind term is a small source for torii only
 * Maximum magnetization parameters should be set high for most problems
 */
#define WIND_TERM 0
#define BSQORHOMAX (100.)
#define UORHOMAX (100.)

/* ELECTRONS AND OPTIONS
 *   SUPPRESS_MAG_HEAT - (0,1) NO ELECTRON HEATING WHEN SIGMA > 1
 *   BETA_HEAT         - (0,1) BETA-DEPENDENT HEATING
 */
#define ELECTRONS           0
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1

/* TIMESTEPPER ALGORITHM
 *  GRIM_TIMESTEPPER - (0,1)
 */

#define GRIM_TIMESTEPPER 1

/* SET_RADIAL_BOUNDS - (0,1)
 * (NEED THIS IF THE DEFINITION OF RIN AND ROUT
 * MUST BE EDITED BASED ON USER INPUT FOR MKS OR FMKS.
 * DEFINE R_inner AS RUNTIME PARAMETER IF YOU'VE SET THIS PARAMETER TO 1)
 */

#define SET_RADIAL_BOUNDS 1

/* RECONSTRUCTION ALGORITHM
 *   LINEAR, PPM, WENO, MP5
 *   SLOPE LIMITERS OPTIONS (USED TO COMPUTE GRADIENTS USING PIECEWISE LINEAR RECONSTRUCTION): 
 *   MC, VAN_LEER
 */
#define RECONSTRUCTION LINEAR
#define SLOPE_LIMITER MC 

/* BOUNDARY CONDITIONS
 *   OUTFLOW PERIODIC POLAR USER
 */
#define X1L_BOUND DIRICHLET
#define X1R_BOUND DIRICHLET
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR
#define X3L_BOUND PERIODIC
#define X3R_BOUND PERIODIC

#define X1L_INFLOW 1
#define X1R_INFLOW 1
