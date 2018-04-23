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
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

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

// Maximum value of gamma, the Lorentz factor
#define GAMMAMAX (50.)
// Average MAD value for gamma (for floors)
#define GAMMAAGN (25.)

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
#define TIMER_DIAG     (9)
#define TIMER_LR_STATE      (10)
#define TIMER_LR_PTOF      (11)
#define TIMER_LR_VCHAR      (12)
#define TIMER_LR_CMAX      (13)
#define TIMER_LR_FLUX      (14)
#define TIMER_BOUND_COMMS (15)
#define TIMER_ALL      (16)
#if ELECTRONS
#define TIMER_ELECTRON_FIXUP (21)
#define TIMER_ELECTRON_HEAT  (22)
#define NUM_TIMERS           (23)
#else
#define NUM_TIMERS     (17)
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

#if RADIATION
extern grid_fourvector_type radG; // Radiation four-force
extern struct of_photon **photon_lists;
#endif

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
#if RADIATION
extern double mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
extern double Ne_unit, Thetae_unit, kphys_to_num;
#endif

// Numerical parameters
extern double Rin, Rout, hslope, R0;
#if RADIATION
extern double Rout_rad, tune_emiss, tune_scatt;
extern double numin, numax;
extern double kappa;
extern double startx_rad[NDIM], stopx_rad[NDIM];
extern double wgtC;
extern int step_made, step_abs, step_scatt, step_lost, step_rec, step_tot;
#endif
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
extern double mdot;
extern double edot;
extern double ldot;
extern int icurr, jcurr, kcurr;

// Parallelism
extern int nthreads;

// Electrons
#if ELECTRONS
extern double game, gamp;
extern double fel0;
#endif

#if RADIATION
struct of_photon {
  double X[NDIM];
  double Kcov[NDIM];
  double Xprev[NDIM];
  double Kcovprev[NDIM];
  double w;
  double dl;
  int nscatt;
  int origin[4];
  struct of_photon *next;
};
#endif

#if POLYTH
extern double poly_norm, poly_xt, poly_alpha;
#endif

extern gsl_rng *r;


// MPI-specific stuff
extern int global_start[3];
extern int global_stop[3];

/*******************************************************************************
    MACROS
*******************************************************************************/
#define ZLOOP \
  for (int k = 0 + NG; k < N3 + NG; k++) \
  for (int j = 0 + NG; j < N2 + NG; j++) \
  for (int i = 0 + NG; i < N1 + NG; i++)
#define ZLOOPALL \
  for (int k = 0; k < N3 + 2*NG; k++) \
  for (int j = 0; j < N2 + 2*NG; j++) \
  for (int i = 0; i < N1 + 2*NG; i++)
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
// Loop over primitive variables
#define PLOOP for(int ip = 0; ip < NVAR; ip++)

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
void bl_coord(double *X, double *r, double *th);
void gcov_func(double *X, double gcov[NDIM][NDIM]);
void set_points();
void set_grid(struct GridGeom *G);
void set_grid_loc(struct GridGeom *G, int i, int j, int k, int loc);
void zero_arrays();

// diag.c
void reset_log_variables();
void diag(struct GridGeom *G, struct FluidState *S, int call_code);
void fail(int fail_type, int i, int j, int k);
//void area_map(int i, int j, int k, grid_prim_type prim);
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

// emissivity.c
#if RADIATION
double jnu(double nu, double Ne, double Thetae, double B, double theta);
double Jnu(double nu, double Ne, double Thetae, double B);
void init_emissivity();
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
void set_core_params();
void set_param(char *key, void *data);
void read_params(char *pfname);
void init_io(char *pfname);
void dump(struct GridGeom *G, struct FluidState *S);
void restart_write(struct FluidState *S);
void restart_read(char *fname, struct FluidState *S);
int restart_init(struct GridGeom *G, struct FluidState *S);
void safe_system(const char *command);
void safe_fscanf(FILE *stream, const char *format, ...);

// make_superphotons.c
#if RADIATION
void make_superphotons(grid_prim_type Prad, double t, double dt);
void set_weight(grid_prim_type Prad);
void get_dnz(grid_prim_type Prad);
#endif

// metric.c
double gcon_func(double gcov[NDIM][NDIM], double gcon[NDIM][NDIM]);
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
MPI_Status* sync_mpi_bound_X1(struct FluidState *S);
MPI_Status* sync_mpi_bound_X2(struct FluidState *S);
MPI_Status* sync_mpi_bound_X3(struct FluidState *S);
void mpi_barrier();
int mpi_nprocs();
double mpi_max(double f);
double mpi_min(double f);
int mpi_io_proc();
void mpi_int_broadcast(int *val);
void mpi_dbl_broadcast(double *val);
double mpi_io_reduce(double val);
double mpi_io_max(double val);
int mpi_myrank();

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

// push_superphotons.c
#if RADIATION
void push_superphotons(double t, double dt);
#endif

// rad_utils.c
#if RADIATION
void init_rad(grid_prim_type Prad);
void init_superphoton_resolution();
double linear_interp_log(double x, double *table, double lmin, double dl);
void set_units();
void list_remove(struct of_photon **ph, struct of_photon **ph_head,
  struct of_photon **ph_prev);
void get_fluid_zone(int i, int j, int k, grid_prim_type Prad, double *Ne,
  double *Thetae, double *B, double Ucon[NDIM], double Ucov[NDIM],
  double Bcon[NDIM], double Bcov[NDIM]);
void *my_malloc(int cnt, int size);
#endif

// radiation.c
#if RADIATION
double Bnu_inv(double nu, double Thetae);
double jnu_inv(double nu, double Thetae, double Ne, double B, double theta);
double alpha_inv_scatt(double nu, double Thetae, double Ne);
double alpha_inv_abs(double nu, double Thetae, double Ne, double B,
  double theta);
double get_fluid_nu(double X[4], double K[4], double Ucov[NDIM]);
double get_bk_angle(double X[NDIM], double K[NDIM], double Ucov[NDIM],
  double Bcov[NDIM], double B);
#endif

// random.c
void init_random(int seed);
double get_random();

// reconstruction.c
void reconstruct(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir);

// step.c
void step(struct GridGeom *G, struct FluidState *state);

// tetrads.c
#if RADIATION
void coord_to_tetrad(double Ecov[NDIM][NDIM], double Kcoord[NDIM],
  double Ktetrad[NDIM]);
void tetrad_to_coord(double Econ[NDIM][NDIM], double Ktetrad[NDIM],
  double Kcoord[NDIM]);
void make_tetrad(int i, int j, int k, double Ucon[NDIM], double trial[NDIM],
  double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM], double Ecov[NDIM][NDIM]);
void normalize_null(double Gcov[NDIM][NDIM], double K[NDIM]);
#endif

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
void dump_grid();
