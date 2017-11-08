/******************************************************************************  
 *                                                                            *  
 * PARAMETERS.H                                                               *  
 *                                                                            *  
 * PROBLEM-SPECIFIC CHOICES                                                   *  
 *                                                                            *  
 ******************************************************************************/

/* GLOBAL RESOLUTION */
#define N1TOT 16
#define N2TOT 16
#define N3TOT 16

/* MPI DECOMPOSITION */
#define N1CPU 1
#define N2CPU 1
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MKS

/* ELECTRONS AND OPTIONS
 *   SUPPRESS_MAG_HEAT - (0,1) NO ELECTRON HEATING WHEN SIGMA > 1
 *   BETA_HEAT         - (0,1) BETA-DEPENDENT HEATING
 *   TPTEMIN           - MINIMUM TP/TE
 */
#define ELECTRONS           0
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1
#define TPTEMIN             0.001
#define TPTEMAX             1000.

/* RADIATION AND OPTIONS
 *   EMISSION   - (0,1)
 *   ABSORPTION - (0,1)
 *   SCATTERING - (0,1)
 *   NU_BINS    - NUMBER OF SAMPLES IN LOG(NU)
 *   THETAE_MAX - MAXIMUM ELECTRON TEMPERATURE SEEN BY RADIATION SECTOR
 *   SIGMA_MAX  - MAXIMUM B.B/RHO ABOVE WHICH FLOW DOES NOT RADIATE
 *   GRAYABSORPTION - (0,1)
 *   BREMSSTRAHLUNG - (0,1)
 *   SYNCHROTRON    - (0,1)
 */
#define RADIATION           0
#define EMISSION            1
#define ABSORPTION          1
#define SCATTERING          1
#define NU_BINS             100
#define THETAE_MAX          1000.
#define SIGMA_MAX           1.
#define GRAYABSORPTION      0
#define BREMSSTRAHLUNG      0
#define SYNCHROTRON         1

/* RECONSTRUCTION ALGORITHM
 *   LINEAR, PPM, WENO, MP5
 */
#define RECONSTRUCTION PPM

/* BOUNDARY CONDITIONS
 *   OUTFLOW PERIODIC POLAR USER
 */
#define X1L_BOUND OUTFLOW
#define X1R_BOUND OUTFLOW
#define X2L_BOUND POLAR
#define X2R_BOUND POLAR
#define X3L_BOUND PERIODIC
#define X3R_BOUND PERIODIC

