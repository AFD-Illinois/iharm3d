/******************************************************************************  
 *                                                                            *  
 * PARAMETERS.H                                                               *  
 *                                                                            *  
 * PROBLEM-SPECIFIC CHOICES                                                   *  
 *                                                                            *  
 ******************************************************************************/

/* GLOBAL RESOLUTION */
#define N1TOT 64
#define N2TOT 64
#define N3TOT 64

/* MPI DECOMPOSITION */
#define N1CPU 2
#define N2CPU 4
#define N3CPU 4

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MINKOWSKI

#define NMODE 3

#define ELECTRONS           0

#define RADIATION           0

/* RECONSTRUCTION ALGORITHM
 *   LINEAR, PPM, WENO, MP5
 */
#define RECONSTRUCTION WENO

/* BOUNDARY CONDITIONS
 *   OUTFLOW PERIODIC POLAR USER
 */
#define X1L_BOUND PERIODIC
#define X1R_BOUND PERIODIC
#define X2L_BOUND PERIODIC
#define X2R_BOUND PERIODIC
#define X3L_BOUND PERIODIC
#define X3R_BOUND PERIODIC
