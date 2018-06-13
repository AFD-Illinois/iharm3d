/******************************************************************************
 *                                                                            *
 * DECS.H                                                                     *
 *                                                                            *
 * GLOBAL MACROS, FUNCTION DEFINITIONS, INCLUDES, AND DECLARATIONS            *
 *                                                                            *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <omp.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

// Hopefully can eliminate this dependency at some point
#include <hdf5.h>

#include "constants.h"
#include "parameters.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169164
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242
#endif

/*******************************************************************************
      COMPILE-TIME PARAMETERS :
*******************************************************************************/

#define VERSION "vharm-alpha-0.1"

// Number of active zones on each MPI process
#define N1       (N1TOT/N1CPU)
#define N2       (N2TOT/N2CPU)
#define N3       (N3TOT/N3CPU)

// Max size for 1D slice is NMAX
#define N12      (N1 > N2 ? N1 : N2)
#define NMAX     (N12 > N3 ? N12 : N3)

#define NDIM       (4)    // Number of total dimensions
#define NPG        (5)    // Number of positions on grid for grid functions
#define NG         (3)    // Number of ghost zones

// Fixup parameters
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#define RHOMIN  (1.e-5)
#define UUMIN (1.e-7)

// Numerical convenience to represent a small (<< 1) non-zero quantity
#define SMALL (1.e-20)
// Set the spatial discretization in numerical derivatives
#define DELTA 1.e-5

// Maximum value of gamma, the Lorentz factor
#define GAMMAMAX (50.)

// Maximum fractional increase in timestep per timestep
#define SAFE  (1.3)

// Whether to move polar axis slightly off of coordinate singularity
#define COORDSINGFIX 1
#define SINGSMALL (1.E-20)

// I/O format strings
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"
#define STRLEN (2048)

// Reconstruction algorithms
#define LINEAR (0)
#define PPM    (1)
#define WENO   (2)
#define MP5    (3)

// Primitive and conserved variables
#define RHO (0)
#define UU  (1)
#define U1  (2)
#define U2  (3)
#define U3  (4)
#define B1  (5)
#define B2  (6)
#define B3  (7)
#if ELECTRONS
#define KEL  (8)
#define KTOT (9)
#define NVAR (10)
#else
#define NVAR (8)
#endif

// Centering of grid functions (assumes axisymmetry in X3)
#define FACE1 (0)
#define FACE2 (1)
//#define FACE3 (2)
#define CORN  (2)
#define CENT  (3)
// For non-axisymmetric (Minkowski spacetimes etc)
#define FACE3 (4)

// Slope limiter
#define MC   (0)
#define VANL (1)
#define MINM (2)

// Boundaries
#define OUTFLOW  (0)
#define PERIODIC (1)
#define POLAR    (2)
#define USER     (3)

// Metric
#define MINKOWSKI (0)
#define MKS       (1)

// Diagnostic calls
#define DIAG_INIT  (0)
#define DIAG_DUMP  (1)
#define DIAG_LOG   (2)
#define DIAG_FINAL (3)

// Failure modes
#define FAIL_UTOPRIM     (0)
#define FAIL_VCHAR_DISCR (1)
#define FAIL_COEFF_NEG   (2)
#define FAIL_COEFF_SUP   (3)
#define FAIL_GAMMA       (4)
#define FAIL_METRIC      (5)

// Timers
#define TIMER_RECON    (1)
#define TIMER_LR_TO_F  (2)
#define TIMER_CMAX     (3)
#define TIMER_FLUX_CT  (4)
#define TIMER_UPDATE_U (5)
#define TIMER_U_TO_P   (6)
#define TIMER_FIXUP    (7)
#define TIMER_BOUND    (8)
#define TIMER_BOUND_COMMS (9)
#define TIMER_DIAG     (10)
#define TIMER_LR_STATE      (11)
#define TIMER_LR_PTOF      (12)
#define TIMER_LR_VCHAR      (13)
#define TIMER_LR_CMAX      (14)
#define TIMER_LR_FLUX      (15)
#define TIMER_IO 	(16)
#define TIMER_RESTART 	(17)
#define TIMER_CURRENT   (18)
#define TIMER_ALL      (19)
#if ELECTRONS
#define TIMER_ELECTRON_FIXUP (21)
#define TIMER_ELECTRON_HEAT  (22)
#define NUM_TIMERS           (23)
#else
#define NUM_TIMERS     (20)
#endif

/*******************************************************************************
    GLOBAL ARRAYS
*******************************************************************************/
typedef int    GridInt[N3+2*NG][N2+2*NG][N1+2*NG];
typedef double GridDouble[N3+2*NG][N2+2*NG][N1+2*NG];
typedef double GridVector[NDIM][N3+2*NG][N2+2*NG][N1+2*NG];
typedef double GridPrim[NVAR][N3+2*NG][N2+2*NG][N1+2*NG];

struct GridGeom {
  double gcov[NPG][NDIM][NDIM][N2+2*NG][N1+2*NG];
  double gcon[NPG][NDIM][NDIM][N2+2*NG][N1+2*NG];
  double gdet[NPG][N2+2*NG][N1+2*NG];
  double lapse[NPG][N2+2*NG][N1+2*NG];
  double conn[NDIM][NDIM][NDIM][N2+2*NG][N1+2*NG];
};

struct FluidState {
  GridPrim P;
  GridPrim U;
  GridVector ucon;
  GridVector ucov;
  GridVector bcon;
  GridVector bcov;
  GridVector jcon;
};

struct FluidFlux {
  GridPrim X1;
  GridPrim X2;
  GridPrim X3;
};

struct FluidEMF {
  GridDouble X1;
  GridDouble X2;
  GridDouble X3;
};

//struct FluidFail {
extern GridInt pflag;
extern GridInt fail_save;
//};

/*******************************************************************************
    GLOBAL VARIABLES SECTION
*******************************************************************************/
// Command line arguments
extern char outputdir[2048], dumpdir[2048], restartdir[2048];

// Physics parameters
extern double a;
extern double gam;
extern double M_unit;
extern double Rhor;
extern double Risco;
extern double tp_over_te;

// Numerical parameters
extern double Rin, Rout, hslope, R0;
extern double cour;
extern double dV, dx[NDIM], startx[NDIM];
extern double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
extern double dt;
extern double t, tf;
extern int nstep;
extern int is_restart;

// Output parameters
extern double DTd;
extern double DTl;
extern int DTr;
extern int DTp;
extern int dump_cnt;
extern int rdump_cnt;
extern double tdump, tlog;

// Global flags
extern int failed;
extern int lim;

// Diagnostics
extern double mdot, edot, ldot;
extern double mdot_eh, edot_eh, ldot_eh;
extern int icurr, jcurr, kcurr;
extern double mass_added;

// Parallelism
extern int nthreads;

// Electrons
#if ELECTRONS
extern double game, gamp;
extern double fel0;
#endif

#if POLYTH
extern double poly_norm, poly_xt, poly_alpha, mks_smooth;
#endif

extern gsl_rng *r; // TODO single-letter practically unused variable


// MPI-specific stuff
extern int global_start[3];
extern int global_stop[3];

/*******************************************************************************
    MACROS
*******************************************************************************/
#define ILOOP	\
  for (int i = 0 + NG; i < N1 + NG; i++)
#define ILOOPALL \
  for (int i = 0; i < N1 + 2*NG; i++)
#define JLOOP	\
  for (int j = 0 + NG; j < N2 + NG; j++)
#define JLOOPALL \
  for (int j = 0; j < N2 + 2*NG; j++)
#define KLOOP	\
  for (int k = 0 + NG; k < N3 + NG; k++)
#define KLOOPALL \
  for (int k = 0; k < N3 + 2*NG; k++)
#define ZLOOP	\
  KLOOP JLOOP ILOOP
#define ZLOOPALL \
  KLOOPALL JLOOPALL ILOOPALL
// Transpose loops for forward-index output
#define ZLOOP_OUT \
  ILOOP JLOOP KLOOP
#define ZLOOP_TRANSPOSE \
  ILOOPALL JLOOPALL KLOOPALL

#define ISLOOP(istart,istop) \
  for (int i = (istart) + NG; i <= (istop) + NG; i++)
#define JSLOOP(jstart,jstop) \
  for (int j = (jstart) + NG; j <= (jstop) + NG; j++)
#define KSLOOP(kstart,kstop) \
  for (int k = (kstart) + NG; k <= (kstop) + NG; k++)
#define ZSLOOP(kstart,kstop,jstart,jstop,istart,istop) \
  for (int k = (kstart) + NG; k <= (kstop) + NG; k++) \
  for (int j = (jstart) + NG; j <= (jstop) + NG; j++) \
  for (int i = (istart) + NG; i <= (istop) + NG; i++)
#define ZSLOOP_REVERSE(kstart,kstop,jstart,jstop,istart,istop) \
  for (int k = (kstop) + NG; k >= (kstart) + NG; k--) \
  for (int j = (jstop) + NG; j >= (jstart) + NG; j--) \
  for (int i = (istop) + NG; i >= (istart) + NG; i--)
#define ZSLOOP_OUT(kstart,kstop,jstart,jstop,istart,istop) \
  ISLOOP(istart,istop) JSLOOP(jstart,jstop) KSLOOP(kstart,kstop)

// Loop over primitive variables
#define PLOOP for(int ip = 0; ip < NVAR; ip++)

// Loop over spacetime indices
#define DLOOP1 for (int mu = 0; mu < NDIM; mu++)
#define DLOOP2 for (int mu = 0; mu < NDIM; mu++)	\
		for (int nu = 0; nu < NDIM; nu++)

#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))
#define MY_SIGN(fval) ( ((fval) <0.) ? -1. : 1. )

#define delta(i,j) ((i == j) ? 1. : 0.)

/*******************************************************************************
    FUNCTION DECLARATIONS
*******************************************************************************/
// bounds.c
void set_bounds(struct GridGeom *geom, struct FluidState *state);
void fix_flux(struct FluidFlux *F);

// coord.c
void coord(int i, int j, int k, int loc, double *X);
double r_of_x(double x);
double dr_dx(double x);
double th_of_x(double x);
double dth_dx(double x);
void bl_coord(const double *X, double *r, double *th);
void cart_coord(const double X[NDIM], double Xcart[NDIM]);
void gcov_func(double *X, double gcov[NDIM][NDIM]);
void set_points();
void set_grid(struct GridGeom *G);
void set_grid_loc(struct GridGeom *G, int i, int j, int k, int loc);
void zero_arrays();

// current.c
void current_calc(struct GridGeom *G, struct FluidState *S, struct FluidState *Ssave, double dtsave);

// diag.c
void reset_log_variables();
void diag(struct GridGeom *G, struct FluidState *S, int call_code);
void fail(struct GridGeom *G, struct FluidState *S, int fail_type, int i, int j, int k);
void area_map(int i, int j, int k, GridPrim prim);
void diag_flux(struct FluidFlux *F);
double flux_ct_divb(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k);

// electrons.c
#if ELECTRONS
void init_electrons();
void heat_electrons(grid_prim_type ph, grid_prim_type p, double Dt);
double get_fel(int i, int j, int k, double p[NVAR]);
void fixup_electrons(grid_prim_type p);
#endif

// fixup.c
void fixup(struct GridGeom *G, struct FluidState *S);
void fixup1zone(struct GridGeom *G, struct FluidState *S, int i, int j, int k);
void fixup_utoprim(struct GridGeom *G, struct FluidState *S);

// fluxes.c
double get_flux(struct GridGeom *G, struct FluidState *S, struct FluidFlux *F);
void lr_to_flux(struct GridGeom *G, struct FluidState *Sl,
  struct FluidState *Sr, int dir, int loc, GridPrim flux, GridDouble ctop);
void flux_ct(struct FluidFlux *F);

// io.c
void init_io();
void dump(struct GridGeom *G, struct FluidState *S);
void dump_grid(struct GridGeom *G);

// metric.c
double gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
void get_gcov(struct GridGeom *G, int i, int j, int loc, double gcov[NDIM][NDIM]);
void get_gcon(struct GridGeom *G, int i, int j, int loc, double gcon[NDIM][NDIM]);
void conn_func(struct GridGeom *G, int i, int j, int k);
void lower_grid(GridVector vcon, GridVector vcov, struct GridGeom *G, int i,
  int j, int k, int loc);
void raise_grid(GridVector vcov, GridVector vcon, struct GridGeom *G, int i,
  int j, int k, int loc);
void lower(double ucon[NDIM], double gcov[NDIM][NDIM], double ucov[NDIM]);
void raise(double ucov[NDIM], double gcon[NDIM][NDIM], double ucon[NDIM]);
double dot_grid(GridVector vcon, GridVector vcov, int i, int j, int k);
double dot(double vcon[NDIM], double vcov[NDIM]);
//double bl_gdet_func(double r, double th);
//void bl_gcov_func(double r, double th, double gcov[][NDIM]);
//void bl_gcon_func(double r, double th, double gcon[][NDIM]);
double MINOR(double m[16], int r0, int r1, int r2, int c0, int c1, int c2);
void adjoint(double m[16], double adjOut[16]);
double determinant(double m[16]);
double invert(double *m, double *invOut);

// mpi.c
void init_mpi();
int sync_mpi_bound_X1(struct FluidState *S);
int sync_mpi_bound_X2(struct FluidState *S);
int sync_mpi_bound_X3(struct FluidState *S);
void mpi_barrier();
int mpi_nprocs();
int mpi_myrank();
double mpi_max(double f);
double mpi_min(double f);
double mpi_reduce(double f);
int mpi_reduce_int(int f);
void mpi_reduce_vector(double *vec_send, double *vec_recv, int len);
int mpi_io_proc();
void mpi_int_broadcast(int *val);
void mpi_dbl_broadcast(double *val);

// params.c
void set_core_params();
void set_param(char *key, void *data);
void read_params(char *pfname);

// phys.c
void prim_to_flux(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int dir, int loc, GridPrim flux);
void prim_to_flux_vec(struct GridGeom *G, struct FluidState *S, int dir,
  int loc, int kstart, int kstop, int jstart, int jstop, int istart, int istop, GridPrim flux);
void bcon_calc(struct FluidState *S, int i, int j, int k);
void mhd_calc(struct FluidState *S, int i, int j, int k, int dir, double *mhd);
void get_fluid_source(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k, GridPrim dU);
double bsq_calc(struct FluidState *S, int i, int j, int k);
void get_state(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc);
void get_state_vec(struct GridGeom *G, struct FluidState *S, int loc,
  int kstart, int kstop, int jstart, int jstop, int istart, int istop);
void ucon_calc(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc);
int mhd_gamma_calc(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k, int loc, double *gamma);
void mhd_vchar(struct GridGeom *G, struct FluidState *Sr, int i, int j, int k,
  int loc, int dir, GridDouble cmax, GridDouble cmin);

// problem.c
void set_problem_params();
void init(struct GridGeom *G, struct FluidState *S);
// Boundary condition (currently used for Bondi flow)
void bound_gas_prob_x1r(int i, int j, int k, GridPrim  P, struct GridGeom *G);

// random.c
void init_random(int seed);
double get_random();

// reconstruction.c
void reconstruct(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir);

// restart.c
void restart_write(struct FluidState *S);
void restart_read(char *fname, struct FluidState *S);
int restart_init(struct GridGeom *G, struct FluidState *S);

// step.c
void step(struct GridGeom *G, struct FluidState *state);

// timing.c
void time_init();
void timer_start(int timerCode);
void timer_stop(int timerCode);
void report_performance();

// u_to_p.c
int U_to_P(struct GridGeom *G, struct FluidState *S, int i, int j, int k,
  int loc);

// xdmf_output.c
void write_xml_closing(FILE *xml);
FILE *write_xml_head(char *fname, double t);
void write_scalar(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml);
void write_scalar_dbl(double *data, hid_t file_id, char *name, hid_t filespace,
  hid_t memspace, FILE *xml);
void write_scalar_int(int *data, hid_t file_id, const char *name,
  hid_t filespace, hid_t memspace, FILE *xml);
void write_vector(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml);
void add_int_value(int val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace);
void add_dbl_value(double val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace);
void add_str_value(const char* val, const char *name, hid_t file_id,
  hid_t filespace, hid_t memspace);

void get_dbl_value(double *val, const char *name, hid_t file_id,
                   hid_t filespace, hid_t memspace);
void get_int_value(int *val, const char *name, hid_t file_id, hid_t filespace,
                   hid_t memspace);

