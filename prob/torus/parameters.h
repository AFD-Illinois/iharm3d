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
#define N3TOT 128

/* MPI DECOMPOSITION */
/* COUNTERINTUITIVE: Split N3, N2, N1 order to keep k smaller than i,j*/
#define N1CPU 1
#define N2CPU 1
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MKS
#define DEREFINE_POLES 1

/*
 * FLOORS
 * Wind term is a small source for torii only
 * Maximum magnetization parameters should be set high for most problems
 */
#define WIND_TERM 0
#define BSQORHOMAX (400.)
#define UORHOMAX (3.)

// Apply the temp ceiling by decreasing internal energy,
// rather than adding mass 
#define TEMP_ADJUST_U 0
// Apply floors in fluid frame
#define FLUID_FRAME_FLOORS 0

/* ELECTRONS AND OPTIONS
 *   SUPPRESS_MAG_HEAT - (0,1) NO ELECTRON HEATING WHEN SIGMA > 1
 *   BETA_HEAT         - (0,1) BETA-DEPENDENT HEATING
 */
#define ELECTRONS           0
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1

/* RECONSTRUCTION ALGORITHM
 *   LINEAR, PPM, WENO, MP5
 */
#define RECONSTRUCTION WENO
// Lower reconstruction order to linear, then
// donor cell, to avoid reconstructing from polar ghosts
#define LOWER_ORDER_POLES 1

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
