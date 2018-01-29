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
#define N1CPU 2
#define N2CPU 2
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MKS

#define POLYTH 0

#define ELECTRONS           0

#define RADIATION           0

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
