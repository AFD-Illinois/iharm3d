/******************************************************************************
 *                                                                            *
 * PARAMETERS.H                                                               *
 *                                                                            *
 * PROBLEM-SPECIFIC CHOICES                                                   *
 *                                                                            *
 ******************************************************************************/

/* GLOBAL RESOLUTION */
#define N1TOT 256
#define N2TOT 256
#define N3TOT 1

/* MPI DECOMPOSITION */
/* COUNTERINTUITIVE: Split N3, N2, N1 order to keep k smaller than i,j*/
#define N1CPU 1
#define N2CPU 1
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MKS
#define DEREFINE_POLES 0

/*
 * FLOORS
 * Wind term is a small source for torii only
 * Maximum magnetization parameters should be set high for most problems
 * Add floors in drift frame?
 */
#define WIND_TERM 0
#define BSQORHOMAX (10.)
#define BSQOUMAX (500.)

#define DRIFT_FRAME 0

/* ELECTRONS AND OPTIONS
 *   ELECTRONS          - (0,1)
 *   MODELS             - (CONSTANT + KAWAZURA + WERNER + ROWAN + SHARMA) //bit field for electron heating models
 *   WHERE, CONSTANT    - (0,1)
 *          KAWAZURA    - (0,2)
 *          WERNER      - (0,4)
 *          ROWAN       - (0,8)
 *          SHARMA      - (0,16)
 *   SUPPRESS_MAG_HEAT  - (0,1) NO ELECTRON HEATING WHEN SIGMA > 1
 *   BETA_HEAT          - (0,1) BETA-DEPENDENT HEATING
 */
#define ELECTRONS           0
#define E_MODELS            1
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1

/* TIMESTEPPER ALGORITHM
 *  IMEX - (0,1)
 *  IF IMEX IS SET TO ZERO, THE CODE DEFAULTS TO THE STANDARD HARM ALGORITHM
 *  LINESEARCH - (0,1)
 */

#define IMEX 0
#define LINESEARCH 0

/* EXTENDED MHD
 *  EMHD - (0,1), CONDUCTION - (0,1), VISCOSITY - (0,1)
 *  SET EMHD TO 1 IF SOLVING AN EXTENDED MHD PROBLEM
 *  NOTE: IF EMHD IS SET TO 1, EITHER CONDUCTION OR VISCOSITY MUST BE SET TO 1
 */

#define EMHD 0
#define CONDUCTION 0
#define VISCOSITY 0

/* DEBUG_EMHD: PRINT VALUES IN A N3D*N2D*N1D STENCIL NEAR THE UPPER LEFT CORNER (PHYSICAL ZONES)
 *  (HOPEFULLY) USEFUL TO DEBUG ISSUES RELATED TO NEW EMHD PROBLEMS
 *  ASSUMES ALL ZONES ARE WITHIN A SINGLE MPI TASK
 */

#define DEBUG_EMHD 0
#if DEBUG_EMHD
#define N1D 223
#define N2D 131
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
#define X1R_BOUND OUTFLOW
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR
#define X3L_BOUND PERIODIC
#define X3R_BOUND PERIODIC

#define X1L_INFLOW 0
#define X1R_INFLOW 0
