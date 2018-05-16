/******************************************************************************
 *                                                                            *
 * DEFS.H                                                                     *
 *                                                                            *
 * GLOBAL VARIABLE DEFINITIONS                                                *
 *                                                                            *
 ******************************************************************************/

/*******************************************************************************
    GLOBAL ARRAYS
*******************************************************************************/
#if RADIATION
grid_fourvector_type radG; // Radiation four-force
struct of_photon **photon_lists;
#endif

/*******************************************************************************
    GLOBAL VARIABLES
*******************************************************************************/
char outputdir[2048], dumpdir[2048], restartdir[2048];

GridInt pflag;
GridInt fail_save;

double a;
double gam;
double M_unit;
double Rhor;
double Risco;
double tp_over_te;
#if RADIATION
double mbh, L_unit, T_unit, M_unit, RHO_unit, U_unit, B_unit;
double Ne_unit, Thetae_unit, kphys_to_num;
#endif

double Rin, Rout, hslope, R0;
#if POLYTH
double poly_norm, poly_xt, poly_alpha, mks_smooth;
#endif
#if RADIATION
double Rout_rad, tune_emiss, tune_scatt;
double numin, numax;
double kappa;
double startx_rad[NDIM], stopx_rad[NDIM];
double wgtC;
int step_made, step_abs, step_scatt, step_lost, step_rec, step_tot;
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

gsl_rng *r;

int global_start[3];
int global_stop[3];

