/******************************************************************************
 *                                                                            *
 * IO.C                                                                       *
 *                                                                            *
 * HDF5 OUTPUT AND RESTART                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#include "hdf5_utils.h"

#include <sys/stat.h>
#include <ctype.h>

// TODO move bsq to this scope?

#if (DEBUG || TEST)
#define OUT_H5_TYPE H5T_IEEE_F64LE
#define OUT_TYPE double
#else
#define OUT_H5_TYPE H5T_IEEE_F32LE
#define OUT_TYPE float
#endif

#define HDF_STR_LEN 30

void dump(struct GridGeom *G, struct FluidState *S)
{
  dump_backend(G, S, IO_REGULAR);
}

void dump_backend(struct GridGeom *G, struct FluidState *S, int type)
{
  timer_start(TIMER_IO);

  static GridDouble *data;
  char fname[80];

  char varNames[NVAR][HDF_STR_LEN] =  {"RHO", "UU", "U1", "U2", "U3"};
  #if EMHD
  #if CONDUCTION
  strcpy(varNames[Q_TILDE], "Q_TILDE");
  #endif
  #if VISCOSITY
  strcpy(varNames[DELTA_P_TILDE], "DELTA_P_TILDE");
  #endif
  #endif
  strcpy(varNames[B1], "B1");
  strcpy(varNames[B2], "B2");
  strcpy(varNames[B3], "B3");
  #if ELECTRONS
  strcpy(varNames[KTOT], "KTOT");

  // Append electron heating model names (preferable as compared to KEL0, KEL1, etc..)
  // To keep track of electron heating model
  int etracker = E_MODELS;
  for (int ip = KTOT + 1; ip < NVAR ; ip++) {
    for (int eFlagIndex = 0; eFlagIndex < sizeof(eFlagsArray) / sizeof(eFlagsArray[0]); eFlagIndex++) {
      if (etracker & eFlagsArray[eFlagIndex]) {
				char kelName[HDF_STR_LEN] = "KEL_";
        strcpy(varNames[ip], strcat(kelName, eNamesArray[eFlagIndex]));
        // update the tracker
        etracker -= eFlagsArray[eFlagIndex];
        break;
      }
    } 
  }
  #endif

  static int firstc = 1;
  if(firstc) {
    data = calloc(1,sizeof(GridDouble));
    firstc = 0;
  }

  //Don't re-dump the grid after a restart
  if (dump_cnt == 0) dump_grid(G);

  if (type == IO_REGULAR) {
    sprintf(fname, "dumps/dump_%08d.h5", dump_cnt);
  } else if (type == IO_ABORT) {
    sprintf(fname, "dumps/dump_abort.h5");    
  }

  if(mpi_io_proc()) fprintf(stdout, "DUMP %s\n", fname);

  hdf5_create(fname);
  //Can use H5T_VARIABLE for any-length strings. But not compatible with parallel IO
  hid_t string_type = hdf5_make_str_type(HDF_STR_LEN);

  // Write header
  hdf5_set_directory("/");
  hdf5_make_directory("header");
  hdf5_set_directory("/header/");
  
  //Adding problem-specific data (problem type, initial data)
  hdf5_make_directory("problem");
  hdf5_set_directory("/header/problem/");
  save_problem_data(string_type);
  
  hdf5_set_directory("/header/");
  hdf5_write_single_val(VERSION, "version", string_type);
  int has_electrons = ELECTRONS;
  hdf5_write_single_val(&has_electrons, "has_electrons", H5T_STD_I32LE);

#if METRIC == MINKOWSKI
  hdf5_write_single_val("MINKOWSKI", "metric", string_type);
#elif METRIC == MKS
#if DEREFINE_POLES
  hdf5_write_single_val("MMKS", "metric", string_type);
#else
  hdf5_write_single_val("MKS", "metric", string_type);
#endif //DEREFINE
#endif //MKS
  char gridfile[HDF_STR_LEN] = "grid.h5"; // TODO match below instead of hard-coding
  hdf5_write_single_val(&gridfile, "gridfile", string_type);

#if RECONSTRUCTION == LINEAR
  hdf5_write_single_val("LINEAR", "reconstruction", string_type);
#elif RECONSTRUCTION == PPM
  hdf5_write_single_val("PPM", "reconstruction", string_type);
#elif RECONSTRUCTION == WENO
  hdf5_write_single_val("WENO", "reconstruction", string_type);
#elif RECONSTRUCTION == MP5
  hdf5_write_single_val("MP5", "reconstruction", string_type);
#endif

  int n1 = N1TOT, n2 = N2TOT, n3 = N3TOT;
  hdf5_write_single_val(&n1, "n1", H5T_STD_I32LE);
  hdf5_write_single_val(&n2, "n2", H5T_STD_I32LE);
  hdf5_write_single_val(&n3, "n3", H5T_STD_I32LE);

  int n_prims = NVAR;
  hdf5_write_single_val(&n_prims, "n_prim", H5T_STD_I32LE);
  // In case we do passive variables
  int n_prims_passive = 0;
  hdf5_write_single_val(&n_prims_passive, "n_prims_passive", H5T_STD_I32LE);
  hdf5_write_str_list(varNames, "prim_names", HDF_STR_LEN, n_prims);

  hdf5_write_single_val(&gam, "gam", H5T_IEEE_F64LE);

#if ELECTRONS
  hdf5_write_single_val(&game, "gam_e", H5T_IEEE_F64LE);
  hdf5_write_single_val(&gamp, "gam_p", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tptemin, "tptemin", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tptemax, "tptemax", H5T_IEEE_F64LE);
  hdf5_write_single_val(&fel0, "fel0", H5T_IEEE_F64LE);
  #if (E_MODELS & CONSTANT)
  hdf5_write_single_val(&fel_constant, "fel_constant", H5T_IEEE_F64LE);
  #endif
#endif

  hdf5_write_single_val(&cour, "cour", H5T_IEEE_F64LE);
  hdf5_write_single_val(&tf, "tf", H5T_IEEE_F64LE);
  hdf5_add_units("tf", "code");

  // EMHD parameters that aren't problem specific
  #if IMEX
  hdf5_make_directory("imex");
  hdf5_set_directory("/header/imex/");
  hdf5_write_single_val(&max_nonlinear_iter, "max_nonlinear_iterations", H5T_STD_I32LE);
  hdf5_write_single_val(&jacobian_eps, "jacobian_eps", H5T_IEEE_F64LE);
  hdf5_write_single_val(&rootfind_tol, "rootfinder_tolerance", H5T_IEEE_F64LE);
  #if LINESEARCH
  hdf5_write_single_val(&max_linesearch_iter, "max_linesearch_iterations", H5T_STD_I32LE);
  hdf5_write_single_val(&linesearch_eps, "linesearch_eps", H5T_IEEE_F64LE);
  #endif
  #if EMHD
  hdf5_make_directory("emhd");
  hdf5_set_directory("/header/imex/emhd/");
  int conduction_enabled = CONDUCTION;
  int viscosity_enabled  = VISCOSITY;
  hdf5_write_single_val(&conduction_enabled, "conduction", H5T_STD_I32LE);
  hdf5_write_single_val(&viscosity_enabled,  "viscosity", H5T_STD_I32LE);
  hdf5_write_single_val(&higher_order_terms_conduction, "higher_order_terms_conduction", H5T_STD_I32LE);
  hdf5_write_single_val(&higher_order_terms_viscosity, "higher_order_terms_viscosity", H5T_STD_I32LE);
  hdf5_write_single_val(&conduction_alpha,   "conduction_alpha", H5T_IEEE_F64LE);
  hdf5_write_single_val(&viscosity_alpha,    "viscosity_alpha", H5T_IEEE_F64LE);
  pack_write_scalar(S->tau, "tau", OUT_H5_TYPE);
  #if CONDUCTION
  pack_write_scalar(S->chi_emhd, "chi_emhd", OUT_H5_TYPE);
  #endif
  #if VISCOSITY
  pack_write_scalar(S->nu_emhd, "nu_emhd", OUT_H5_TYPE);
  #endif
  #endif
  #endif


  //Geometry
  hdf5_set_directory("/header/");
  hdf5_make_directory("geom");
  hdf5_set_directory("/header/geom/");

  hdf5_write_single_val(&(startx[1]), "startx1", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(startx[2]), "startx2", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(startx[3]), "startx3", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(dx[1]), "dx1", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(dx[2]), "dx2", H5T_IEEE_F64LE);
  hdf5_write_single_val(&(dx[3]), "dx3", H5T_IEEE_F64LE);
  int n_dim = NDIM;
  hdf5_write_single_val(&n_dim, "n_dim", H5T_STD_I32LE);
#if METRIC == MKS
#if DEREFINE_POLES
  hdf5_make_directory("mmks");
  hdf5_set_directory("/header/geom/mmks/");
  hdf5_write_single_val(&poly_xt, "poly_xt", H5T_IEEE_F64LE);
  hdf5_write_single_val(&poly_alpha, "poly_alpha", H5T_IEEE_F64LE);
  hdf5_write_single_val(&mks_smooth, "mks_smooth", H5T_IEEE_F64LE);
#else
  hdf5_make_directory("mks");
  hdf5_set_directory("/header/geom/mks/");
#endif
  hdf5_write_single_val(&Rin, "r_in", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Rout, "r_out", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Rhor, "r_eh", H5T_IEEE_F64LE);
  // I don't need this in code but it's in the spec to output it
  double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  double Risco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));
  hdf5_write_single_val(&Risco, "r_isco", H5T_IEEE_F64LE);
  hdf5_write_single_val(&hslope, "hslope", H5T_IEEE_F64LE);
  hdf5_write_single_val(&a, "a", H5T_IEEE_F64LE);
#endif

  hdf5_set_directory("/");

  int is_full_dump = 1; // TODO do partial/full dumps
  hdf5_write_single_val(&is_full_dump, "is_full_dump", H5T_STD_I32LE);
  hdf5_write_single_val(&t, "t", H5T_IEEE_F64LE);
  hdf5_add_units("t", "code");
  hdf5_write_single_val(&dt, "dt", H5T_IEEE_F64LE);
  hdf5_add_units("dt", "code");
  hdf5_write_single_val(&nstep, "n_step", H5T_STD_I32LE);
  hdf5_write_single_val(&dump_cnt, "n_dump", H5T_STD_I32LE);

  hdf5_write_single_val(&DTd, "dump_cadence", H5T_IEEE_F64LE);
  hdf5_write_single_val(&DTf, "full_dump_cadence", H5T_IEEE_F64LE);

  // Write primitive variables
  pack_write_vector(S->P, NVAR, "prims", OUT_H5_TYPE);
  hdf5_add_units("prims", "code");

  // Write jcon (not recoverable from prims)
  pack_write_vector(S->jcon, NDIM, "jcon", OUT_H5_TYPE);
  hdf5_add_units("jcon", "code");

  ZLOOP (*data)[k][j][i] = mhd_gamma_calc(G, S, i, j, k, CENT);
  pack_write_scalar((*data), "gamma", OUT_H5_TYPE);

  // Space for any extra items
  // Currently debug/diagnostic output, on full dumps only
  hdf5_make_directory("extras");
  hdf5_set_directory("/extras/");

  // Preserve git commit or tag of the run -- see Makefile
#ifdef GIT_VERSION
  hdf5_write_single_val(QUOTE(GIT_VERSION), "git_version", string_type);
#endif

  if (is_full_dump) {
    // TODO need sync here for consistent output?
    ZLOOP (*data)[k][j][i] = flux_ct_divb(G, S, i, j, k);
    pack_write_scalar((*data), "divB", OUT_H5_TYPE);

    pack_write_int(fail_save, "fail");
    ZLOOP fail_save[k][j][i] = 0;

    pack_write_int(fflag, "fixup");

    // Write out the L2 norm of the residuals.
    #if IMEX
    pack_write_scalar(imex_errors, "imex_errors", OUT_H5_TYPE);
    ZLOOP imex_errors[k][j][i] = 0.;
    #endif
  }

#if DEBUG
    pack_write_vector(S->U, NVAR, "U", OUT_H5_TYPE);

    pack_write_vector(preserve_F.X1, NVAR, "X1", OUT_H5_TYPE);
    pack_write_vector(preserve_F.X2, NVAR, "X2", OUT_H5_TYPE);
    pack_write_vector(preserve_F.X3, NVAR, "X3", OUT_H5_TYPE);
    pack_write_vector(preserve_dU, NVAR, "dU", OUT_H5_TYPE);
#endif

    // Extra physical variables.  These are just recalculated in post
//  ZLOOP (*data)[k][j][i] = bsq_calc(G, S, i, j, k);
//  pack_write_scalar((*data), "bsq", OUT_H5_TYPE);

  // Could write extra vectors here, but they take up too much space
//  pack_write_vector(S->bcon, NDIM, "bcon", OUT_H5_TYPE);
//  hdf5_add_units("bcon", "code");

  hdf5_close();

  timer_stop(TIMER_IO);
}

#define NGRIDVARS 9
void dump_grid(struct GridGeom *G)
{
  GridDouble *x[NGRIDVARS];
  for (int d = 0; d < NGRIDVARS; d++) x[d] = calloc(1,sizeof(GridDouble));
  const char *coordNames[] = {"X", "Y", "Z", "r", "th", "phi", "X1", "X2", "X3"};

  char *fname = "dumps/grid.h5";
  if(mpi_io_proc()) fprintf(stdout, "GRID %s\n", fname);

  hdf5_create(fname);

  hdf5_set_directory("/");

  // Batch fill grid var buffers since lots of them are the same
  ZLOOP {
	double xp[4];
    coord(i, j, k, CENT, xp);
    #if METRIC == MINKOWSKI
    (*x[0])[k][j][i] = xp[1];
    (*x[1])[k][j][i] = xp[2];
    (*x[2])[k][j][i] = xp[3];
    (*x[3])[k][j][i] = 0;
    (*x[4])[k][j][i] = 0;
    (*x[5])[k][j][i] = 0;
    #elif METRIC == MKS
    double r, th;
    bl_coord(xp, &r, &th);
    (*x[0])[k][j][i] = r*cos(xp[3])*sin(th);
    (*x[1])[k][j][i] = r*sin(xp[3])*sin(th);
    (*x[2])[k][j][i] = r*cos(th);
    (*x[3])[k][j][i] = r;
    (*x[4])[k][j][i] = th;
    (*x[5])[k][j][i] = xp[3];
    #endif

    (*x[6])[k][j][i] = xp[1];
    (*x[7])[k][j][i] = xp[2];
    (*x[8])[k][j][i] = xp[3];
  }

  for (int d = 0; d < NGRIDVARS; d++){
    pack_write_scalar((*x[d]), coordNames[d], OUT_H5_TYPE);
  }

  pack_write_axiscalar(G->gdet[CENT], "gdet", OUT_H5_TYPE);
  pack_write_axiscalar(G->lapse[CENT], "lapse", OUT_H5_TYPE);

  pack_write_Gtensor(G->gcon[CENT], "gcon", OUT_H5_TYPE);
  pack_write_Gtensor(G->gcov[CENT], "gcov", OUT_H5_TYPE);

  for (int d = 0; d < NGRIDVARS; d++) free(x[d]);

  hdf5_close();
}
