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
/* Careful NXCPU < NXTOT!! */
#define N1CPU 1
#define N2CPU 1
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MKS
#define DEREFINE_POLES 0

/* FLOORS
 *   Wind term is a small source for torii only
 *   Maximum magnetization parameters should be set high for most problems
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
 *  IMEX - (0,1)
 *  IF IMEX IS SET TO ZERO, THE CODE DEFAULTS TO THE STANDARD HARM ALGORITHM
 */

#define IMEX 1

/* EXTENDED MHD
 *  EMHD - (0,1), CONDUCTION - (0,1), VISCOSITY - (0,1)
 *  SET EMHD TO 1 IF SOLVING AN EXTENDED MHD PROBLEM
 *  NOTE: IF EMHD IS SET TO 1, EITHER CONDUCTION OR VISCOSITY MUST BE SET TO 1
 *  SET CONDUCTION TO 1 IF HEAT CONDUCTION IS PRESENT
 *  SET VISCOSITY TO 1 IF VISCOSITY IS PRESENT
 */

#define EMHD 1
#define CONDUCTION 1
#define VISCOSITY 1

/* DEBUG_EMHD mode: PRINT VALUES IN A N3D*N2D*N1D STENCIL NEAR THE UPPER LEFT CORNER (PHYSICAL ZONES)
 *  (HOPEFULLY) USEFUL TO DEBUG ISSUES RELATED TO NEW EMHD PROBLEMS
 *  ASSUMES ALL ZONES ARE WITHIN A SINGLE MPI TASk
 */

#define DEBUG_EMHD 0
#if DEBUG_EMHD
#define N1D 8
#define N2D 8
#define N3D 1
#endif

/* RECONSTRUCTION ALGORITHM
 *   LINEAR, PPM, WENO, MP5
 */
#define RECONSTRUCTION WENO

/* BOUNDARY CONDITIONS
 *   OUTFLOW PERIODIC POLAR USER
 */
#define X1L_BOUND OUTFLOW
#define X1R_BOUND USER
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR
#define X3L_BOUND PERIODIC
#define X3R_BOUND PERIODIC

#define X1L_INFLOW 0
#define X1R_INFLOW 1
