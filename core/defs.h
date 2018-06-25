/******************************************************************************
 *                                                                            *
 * DEFS.H                                                                     *
 *                                                                            *
 * GLOBAL VARIABLE DEFINITIONS                                                *
 *                                                                            *
 ******************************************************************************/

/*******************************************************************************
    GLOBAL VARIABLES
*******************************************************************************/
// TODO encapsulate better maybe

char dumpdir[2048], restartdir[2048];

GridInt pflag;
GridInt fail_save;

double a;
double gam;
double M_unit;
double Rhor;
double Risco;
double tp_over_te;

double Rin, Rout, hslope, R0;
#if POLYTH
double poly_norm, poly_xt, poly_alpha, mks_smooth;
#endif
double cour;
double dV, dx[NDIM], startx[NDIM];
double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
double dt;
double t, tf;
double rcurr, hcurr;
int istart, istop, jstart, jstop;
int nstep;
int is_restart;

double DTd;
double DTf;
double DTl;
int DTr;
int DTp;
int dump_cnt;
int rdump_cnt;
double tdump, tlog;

int failed;
int lim;

double mdot, edot, ldot;
double mdot_eh, edot_eh, ldot_eh;
int icurr, jcurr, kcurr;

int nthreads;

#if ELECTRONS
double game, gamp;
double fel0;
#endif

int global_start[3];
int global_stop[3];

