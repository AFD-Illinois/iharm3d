/******************************************************************************
 *                                                                            *
 * DEFS.H                                                                     *
 *                                                                            *
 * GLOBAL VARIABLE DEFINITIONS                                                *
 *                                                                            *
 ******************************************************************************/

#pragma once

// Zone flags.  TODO move these to the heap
GridInt pflag;
GridInt fail_save;
GridInt fflag;

#if DEBUG
struct FluidFlux preserve_F;
GridPrim preserve_dU;
#endif

// Parameters
// physical
double a;
double gam;
double Rhor;
double tp_over_te;

// geometry
double Rin, Rout, hslope;
double poly_norm, poly_xt, poly_alpha, mks_smooth;
double cour;
double dV, dx[NDIM], startx[NDIM];
double x1Min, x1Max, x2Min, x2Max, x3Min, x3Max;
double dt, dt_light;
double t, tf;
double rcurr, hcurr;
int istart, istop, jstart, jstop;
int nstep;
int is_restart;

// fluid dumps
double DTd;
double DTf;
double DTl;
int DTr;
int DTp;
int dump_cnt;
double tdump, tlog;

// derived logged output
double mdot, edot, ldot;
double mdot_eh, edot_eh, ldot_eh;
int icurr, jcurr, kcurr;

int nthreads;

#if ELECTRONS
double game, gamp;
double fel0;
double tptemin, tptemax;
// Array of electron heating model flags (need it to call get_fels())
int eFlagsArray[] = {CONSTANT, KAWAZURA, WERNER, ROWAN, SHARMA};
// Array of electron heating model names (need it to assign primNames)
char eNamesArray[NUM_E_MODELS][20] = {"CONSTANT", "KAWAZURA", "WERNER", "ROWAN", "SHARMA"};
#if (E_MODELS & CONSTANT)
double fel_constant;
#endif
#endif

int global_start[3];
int global_stop[3];

#if IMEX
int max_nonlinear_iter;
double jacobian_eps;
double rootfind_tol;
#endif

#if EMHD
int higher_order_terms_conduction;
int higher_order_terms_viscosity;
double conduction_alpha;
double viscosity_alpha;
#endif

#if SET_RADIAL_BOUNDS
double R_inner;
#endif
