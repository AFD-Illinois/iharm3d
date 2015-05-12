
/* M1: stripped conduction, modified dtsave, psave */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <string.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>
//#include <gsl/gsl_linalg.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923132169164
#endif
#ifndef M_SQRT2
#define M_SQRT2 1.4142135623730950488016887242
#endif

extern int flag_of_convenience ;

/*************************************************************************
      COMPILE-TIME PARAMETERS : 
*************************************************************************/
/** SWITCHES **/

/* include wind source terms */
#define WIND        (0)

/* use linear reconstruction */
#define RECON_LIN   (1)

/* use parabolic reconstruction */
#define RECON_PARA  (0)

/* use WENO reconstruction */
#define RECON_WENO  (0)

/* use optically thin cooling */
#define ALLOW_COOL  (0)

/** END SWITCHES **/

/** here are the few things that we change frequently **/
#define N1TOT    (96)
#define N2TOT    (64)
#define N3TOT    (64)

#define N1CPU    (1)
#define N2CPU    (1)
#define N3CPU    (1)


#define N1       (N1TOT/N1CPU)		/* number of physical zones in X1-direction */
#define N2       (N2TOT/N2CPU)		/* number of physical zones in X2-direction */
#define N3       (N3TOT/N3CPU)
#define N12      (N1 > N2 ? N1 : N2) 
#define NMAX     (N12 > N3 ? N12 : N3)	/* this sizes 1D slices */

#define NPR        (8) 		/* number of primitive variables */
#define NDIM       (4)		/* number of total dimensions.  Never changes */
#define NPG        (4)		/* number of positions on grid for grid functions */
#define NG         (3)		/* number of ghost zones */

/** FIXUP PARAMETERS, magnitudes of rho and u, respectively, in the floor : **/
#define RHOMINLIMIT (1.e-20)
#define UUMINLIMIT  (1.e-20)
#define RHOMIN	(1.e-5)
#define UUMIN	(1.e-7)

/* A numerical convenience to represent a small non-zero quantity compared to unity:*/
#define SMALL	(1.e-20)

/* Max. value of gamma, the lorentz factor */
#define GAMMAMAX (50.)

/* maximum fractional increase in timestep per timestep */
#define SAFE	(1.3)

#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
#define SINGSMALL (1.E-20)

/* I/O format strings used herein : */
#define FMT_DBL_OUT "%28.18e"
#define FMT_INT_OUT "%10d"


/*************************************************************************
    MNEMONICS SECTION 
*************************************************************************/
/* boundary condition mnemonics */
#define OUTFLOW	 (0)
#define SYMM	 (1)
#define ASYMM	 (2)
#define FIXED	 (3)
#define PERIODIC (4)

/* mnemonics for primitive vars; conserved vars */
#define RHO	(0)
#define UU	(1)
#define U1	(2)
#define U2	(3)
#define U3	(4)
#define B1	(5)
#define B2	(6)
#define B3	(7)

/* mnemonics for dimensional indices */
#define TT	(0)
#define RR	(1)
#define TH	(2)
#define PH	(3)

/* mnemonics for centering of grid functions */
#define FACE1	(0)
#define FACE2	(1)
#define CORN	(2)
#define CENT	(3)
#define FACE3	(CENT)	// appropriate for axisymmetric metrics (Kerr!)

/* mnemonics for slope limiter */
#define MC	(0)
#define VANL	(1)
#define MINM	(2)

/* mnemonics for diagnostic calls */
#define INIT_OUT	(0)
#define DUMP_OUT	(1)
#define IMAGE_OUT	(2)
#define LOG_OUT		(3)
#define FINAL_OUT	(4)
#define PDUMP_OUT	(5)

/* Directional Mnemonics */
// -------------> r
// |         3    
// |        1-0   
// |         2    
// v            
// theta      
#define X1UP    (0)
#define X1DN    (1)
#define X2UP    (2)
#define X2DN    (3)
#define X3UP	(4)
#define X3DN	(5)


/* failure modes */
#define FAIL_UTOPRIM        (1)
#define FAIL_VCHAR_DISCR    (2)
#define FAIL_COEFF_NEG	    (3)
#define FAIL_COEFF_SUP	    (4)
#define FAIL_GAMMA          (5)
#define FAIL_METRIC         (6)

/*************************************************************************
    GLOBAL ARRAY SECTION 
*************************************************************************/
typedef double	grid_prim_type[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG][NPR];

extern grid_prim_type p;	/* space for primitive vars */
extern grid_prim_type dq;	/* slopes */
extern grid_prim_type F1;	/* fluxes */
extern grid_prim_type F2;	/* fluxes */
extern grid_prim_type F3;
extern grid_prim_type ph;	/* half-step primitives */
extern int    pflag[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];	/* identifies failure points */
extern int    fail_save[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];

/* for debug & diagnostics */
extern grid_prim_type psave;   /* stores old data for time derivatives */
extern double dtsave ;
extern double Jcon[N1+2*NG][N2+2*NG][N3+2*NG][NDIM];

/* particle-related global variables */
#define NPTOT	0
extern double xp[NPTOT][NDIM];

/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
extern double a;
extern double gam;
extern double M_unit;
extern double Rhor;
extern double Risco;

/* numerical parameters */
extern double Rin, Rout, hslope, R0;
extern double cour;
extern double dV, dx[NPR], startx[NPR];
extern double dt;
extern double t, tf;
extern double x1curr, x2curr;
extern int nstep;

/* output parameters */
extern double DTd;
extern double DTl;
extern double DTi;
extern double DTp;
extern int DTr;
extern int dump_cnt;
extern int pdump_cnt;
extern int image_cnt;
extern int rdump_cnt;
extern int nstroke;
extern double t_last_dump;
extern double t_next_dump;

/* global flags */
extern int failed;
extern int lim;
extern double defcon;

/* diagnostics */
extern double mdot;
extern double edot;
extern double ldot;

/* set global variables that indicate current local metric, etc. */
extern int icurr, jcurr, pcurr;
struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
	double alpha;
};

struct of_state {
	double ucon[NDIM];
	double ucov[NDIM];
	double bcon[NDIM];
	double bcov[NDIM];
};

//extern gsl_rng *r ;

/* more grid functions */
// -- for now assume axisymmetry (JCD)
extern double conn[N1 + 2*NG][N2 + 2*NG][NDIM][NDIM][NDIM];
extern struct of_geom ggeom[N1+2*NG][N2+2*NG][NPG] ;

/* mpi specific stuff */
extern int global_start[3];
extern int global_stop[3];

/*************************************************************************
    MACROS
*************************************************************************/
#define START1	(0+NG)
#define START2	(0+NG)
#define START3	(0+NG)

/* loop over all active zones */
#define ZLOOP for(i=0+START1;i<N1+START1;i++)for(j=0+START2;j<N2+START2;j++)for(k=0+START3;k<N3+START3;k++)

/* loop over all active zones */
#define IMAGELOOP for(j=0+START2;j<N2+START2;j++)for(i=0+START1;i<N1+START1;i++)

/* specialty loop */
extern int istart, istop, jstart, jstop, kstart, kstop;
#define ISLOOP(istart,istop) \
        for(i=istart+START1;i<=istop+START1;i++)
#define JSLOOP(jstart,jstop) \
	for(j=jstart+START2;j<=jstop+START2;j++)
#define KSLOOP(kstart,kstop) \
	for(k=kstart+START3;k<=kstop+START3;k++)
#define ZSLOOP(istart,istop,jstart,jstop,kstart,kstop) \
        for(i=istart+START1;i<=istop+START1;i++)\
	for(j=jstart+START2;j<=jstop+START2;j++)\
	for(k=kstart+START3;k<=kstop+START3;k++)
#define GZSLOOP(istart,istop,jstart,jstop) \
        for(i=istart+START1;i<=istop+START1;i++)\
	for(j=jstart+START2;j<=jstop+START2;j++)

/* loop over Primitive variables */
#define PLOOP  for(int ip=0;ip<NPR;ip++)
/* loop over all Dimensions; second rank loop */
#define DLOOP  for(j=0;j<NDIM;j++) for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP  for(j=1;j<NDIM;j++) for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)


extern double fval1, fval2;
#define MY_MIN(fval1,fval2) ( ((fval1) < (fval2)) ? (fval1) : (fval2))
#define MY_MAX(fval1,fval2) ( ((fval1) > (fval2)) ? (fval1) : (fval2))
#define MY_SIGN(fval) ( ((fval) <0.) ? -1. : 1. )

#define delta(i,j) ( (i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])


/*************************************************************************
    FUNCTION DECLARATIONS 
*************************************************************************/
double bl_gdet_func(double r, double th);
double bsq_calc(double *pr, struct of_geom *geom);
int mhd_gamma_calc(double *pr, struct of_geom *geom, double *gamma);
int rad_gamma_calc(double *pr, struct of_geom *geom, double *gamma);
double gdet_func(double lgcov[][NDIM]);
double mink(int j, int k);
double ranc() ;
double slope_lim(double y1, double y2, double y3);
double synchrotron_cooling_func(double rho, double u, double bsq, double r);

int restart_init();

struct of_geom *get_geometry(int i, int j, int loc);

void advance_particles(grid_prim_type prim, double Dt);
void area_map(int i, int j, int k, grid_prim_type prim);
void bcon_calc(double *pr, double *ucon, double *ucov, double *bcon);
void blgset(int i, int j, struct of_geom *geom);
void bl_coord(double *X, double *r, double *th);
double r_of_x(double x);
double th_of_x(double x);
double dr_dx(double x);
double dth_dx(double x);
void bl_gcon_func(double r, double th, double gcov[][NDIM]);
void bl_gcov_func(double r, double th, double gcov[][NDIM]);
void bound_prim(grid_prim_type pr);
void current_calc() ;
void conn_func(double *X, struct of_geom *geom,
	       double lconn[][NDIM][NDIM]);
void cool_down(grid_prim_type, double Dt);
void coord(int i, int j, int loc, double *X);
void diag(int call_code);
void diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
double flux_ct_divb(int i, int j, int k);
void dump();

void fail(int fail_type);
void fixup(grid_prim_type pv);
void fixup1zone(int i, int j, int k, double prim[NPR]);
void fixup_utoprim(grid_prim_type pv);
void fix_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3);
double invert(double *A, double *invA);	// inverts a 4x4
double gcon_func(double lgcov[][NDIM], double lgcon[][NDIM]);
void gcov_func(double *X, double lgcov[][NDIM]);
void get_state(double *pr, struct of_geom *geom, struct of_state *q);
void image_all(int image_count);
void init(void);
void init_particles(void);
void init_ranc(int seed) ;
void linear_mc(double x1, double x2, double x3, double *lout, double *rout) ;
void lower(double *a, struct of_geom *geom, double *b);
void lr_to_flux(double p_l[NPR], double p_r[NPR], struct of_geom *geom, int dir, 
	double Flux[NPR], double *max_speed) ;
void ludcmp(double **a, int n, int *indx, double *d);
void mhd_calc(double *pr, int dir, struct of_state *q, double *mhd);
void misc_source(double *ph, int ii, int jj, struct of_geom *geom,
		 struct of_state *q, double *dU, double Dt);
void para(double x1, double x2, double x3, double x4, 
	double x5, double *lout, double *rout) ;
void pdump(FILE * fp) ;
void primtoflux(double *pa, struct of_state *q, int dir,
		struct of_geom *geom, double *fl);
void primtoU(double *p, struct of_state *q, struct of_geom *geom,
	     double *U);
void rad_calc(double *pr, int dir, struct of_state *q, double *rad);
void raise(double *v1, struct of_geom *geom, double *v2);
void reconstruct_lr_lin(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void reconstruct_lr_par(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void reconstruct_lr_weno(double ptmp[][NPR], int N, double p_l[][NPR], double p_r[][NPR]);
void restart_write(void);
void restart_read(char *fname);
void zero_arrays(void);
void set_grid(void);
void set_points(void);
void step_ch(void);
void source(double *pa, struct of_geom *geom, int ii, int jj, double *Ua,
	    double Dt);
void timestep(void);
void u_to_v(double *pr, int i, int j);
void ucon_calc(double *pr, struct of_geom *geom, double *ucon);
void uRcov_calc(double *pr, struct of_geom *geom, double *uRcov);
void usrfun(double *pr, int n, double *beta, double **alpha);
void weno(double x1, double x2, double x3, double x4, 
	double x5, double *lout, double *rout) ;
int Utoprim_mm(double *Ua, struct of_geom *geom, double *pa);
int Utoprim_sn(double *Ua, struct of_geom *geom, double *pa);

void mhd_vchar(double *pr, struct of_state *q, struct of_geom *geom,
	   int dir, double *cmax, double *cmin);
void rad_vchar(double *pr, struct of_state *q, struct of_geom *geom,
	   int dir, double *cmax, double *cmin);
void wind_source(double *ph, struct of_geom *geom, struct of_state *q,
	int ii, int jj, double *dU);

void sync_mpi_boundaries(grid_prim_type pr);
double mpi_max(double f);
double mpi_min(double f);
int mpi_io_proc();
void mpi_int_broadcast(int *val);
void mpi_dbl_broadcast(double *val);
double mpi_io_reduce(double val);
double mpi_io_max(double val);
int mpi_nprocs();
int mpi_myrank();

void write_xml_closing(FILE *xml);
FILE *write_xml_head(int dump_id, double t);
void dump_grid();
