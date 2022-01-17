/******************************************************************************
 *                                                                            *
 * PARAMETERS.H                                                               *
 *                                                                            *
 * PROBLEM-SPECIFIC CHOICES                                                   *
 *                                                                            *
 ******************************************************************************/

/* GLOBAL RESOLUTION */
#define N1TOT 2000
#define N2TOT 1
#define N3TOT 1

/* MPI DECOMPOSITION */
/* DECOMPOSE IN N3 FIRST! Small leading array sizes for linear access */
#define N1CPU 1
#define N2CPU 1
#define N3CPU 1

/* METRIC
 *   MINKOWSKI, MKS
 */
#define METRIC MINKOWSKI

/*
 * FLOORS
 * Wind term is a small source for torii only
 * Maximum magnetization parameters should be set high for most problems
 */
#define WIND_TERM 0
#define BSQORHOMAX (100.)
#define UORHOMAX (100.)

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
#define ELECTRONS           1
#define E_MODELS            1
#define SUPPRESS_HIGHB_HEAT 1
#define BETA_HEAT           1

/* RECONSTRUCTION ALGORITHM
 *   LINEAR, PPM, WENO, MP5
 */
#define RECONSTRUCTION WENO

/* BOUNDARY CONDITIONS
 *   OUTFLOW PERIODIC POLAR USER
 */
#define X1L_BOUND OUTFLOW
#define X1R_BOUND OUTFLOW
#define X2L_BOUND PERIODIC
#define X2R_BOUND PERIODIC
#define X3L_BOUND PERIODIC
#define X3R_BOUND PERIODIC
