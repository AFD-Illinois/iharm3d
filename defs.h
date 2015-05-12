
/* M1: modified op, odt after removing conduction */

int flag_of_convenience;

/*************************************************************************
    GLOBAL ARRAYS SECTION 
*************************************************************************/
grid_prim_type p;	/* space for primitive vars */
grid_prim_type dq;	/* slopes */
grid_prim_type F1;	/* fluxes */
grid_prim_type F2;	/* fluxes */
grid_prim_type F3;
grid_prim_type ph;	/* half-step primitives */
int pflag[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];	/* identifies failure points */
int fail_save[N1 + 2*NG][N2 + 2*NG][N3 + 2*NG];

/* for debug & diagnostics */
grid_prim_type psave;
double Jcon[N1+2*NG][N2+2*NG][N3+2*NG][NDIM];
//double Ccov[N1+2*NG][N2+2*NG][NDIM];
double dtsave ;

/* grid functions */
double conn[N1 + 2*NG][N2 + 2*NG][NDIM][NDIM][NDIM];
struct of_geom ggeom[N1+2*NG][N2+2*NG][NPG] ;

/* particles */
double xp[NPTOT][NDIM];

/*************************************************************************
    GLOBAL VARIABLES SECTION 
*************************************************************************/
/* physics parameters */
double a;
double gam;
double M_unit;
double Rhor;
double Risco;

/* numerical parameters */
double Rin, Rout, hslope, R0;
double cour;
double dV, dx[NPR], startx[NPR];
double dt;
double t, tf;
double rcurr, hcurr;
int istart, istop, jstart, jstop;
int icurr, jcurr, pcurr, ihere, jhere, phere;
double dminarg1, dminarg2;
int nstep;
double fval1, fval2;

/* output parameters */
double DTd;
double DTp;
double DTl;
double DTi;
int DTr;
int dump_cnt;
int pdump_cnt;
int image_cnt;
int rdump_cnt;
int nstroke;
double t_last_dump;
double t_next_dump;

/* global flags */
int failed;
int lim;
double defcon;

/* diagnostics */
double mdot = 0.;
double edot = 0.;
double ldot = 0.;

/* current local position */
int icurr, jcurr, pcurr;

//gsl_rng *r ;

/* mpi specific stuff */
int global_start[3];
int global_stop[3];

