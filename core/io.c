/******************************************************************************
 *                                                                            *
 * IO.C                                                                       *
 *                                                                            *
 * HDF5 OUTPUT AND RESTART                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"
#include <hdf5.h>
#include <sys/stat.h>
#include <ctype.h>

// Variables kept for output only
GridDouble *omega;

// TODO also need to change which packing routines are called
#define OUT_TYPE float
#define OUT_H5_TYPE H5T_NATIVE_FLOAT

#define HDF_STR_LEN 20

void init_io()
{
  omega = calloc(1,sizeof(GridDouble));
}

void dump(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_IO);

  static int firstc = 1;
  static OUT_TYPE *data;
  static int *idata;
  char fname[80];

  #if ELECTRONS
  const char varNames[NVAR][HDF_STR_LEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3",
                            "KEL", "KTOT"};
  #else
  const char varNames[NVAR][HDF_STR_LEN] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"}; //Reserve some extra
  #endif

  if(firstc) {
    data = calloc(N1*N2*N3*NVAR,sizeof(OUT_TYPE));
    idata = calloc(N1*N2*N3,sizeof(int));
    if(data == NULL) {
      fprintf(stderr,"failed to allocate data in dump\n");
      exit(8);
    }
    firstc = 0;
  }

  //Don't re-dump the grid on restart
  if (dump_cnt == 0) dump_grid(G);

  sprintf(fname, "dumps/dump_%08d.h5", dump_cnt);
  if(mpi_io_proc()) {
    fprintf(stdout, "DUMP %s\n", fname);
  }

  hid_t file_id = hdf5_open(fname);
  hid_t string_type = hdf5_make_str_type(HDF_STR_LEN); //H5T_VARIABLE for any-length strings. But not compat with parallel IO

  // Write header
  hdf5_make_directory("header", file_id);
  hdf5_set_directory("/header/");

  hdf5_write_single_val(VERSION, "version", file_id, string_type);
  int has_electrons = ELECTRONS;
  hdf5_write_single_val(&has_electrons, "has_electrons", file_id, H5T_NATIVE_INT);

#if METRIC == MINKOWSKI
  hdf5_write_single_val("MINKOWSKI", "metric_name", file_id, string_type);
#elif METRIC == MKS
#if POLYTH // Morty Maxwell's Massively Modified Kerr-Schild Coordinates
  hdf5_write_single_val("MMKS", "metric_name", file_id, string_type);
#else
  hdf5_write_single_val("MKS", "metric_name", file_id, string_type);
#endif //POLYTH
#endif //MKS
  char *gridfile_name = "grid.h5"; // TODO follow grid below?
  hdf5_write_single_val(&gridfile_name, "gridfile_name", file_id, string_type);

  // TODO caps?
#if RECONSTRUCTION == LINEAR
  hdf5_write_single_val("LINEAR", "reconstruction_name", file_id, string_type);
#elif RECONSTRUCTION == PPM
  hdf5_write_single_val("PPM", "reconstruction_name", file_id, string_type);
#elif RECONSTRUCTION == WENO
  hdf5_write_single_val("WENO", "reconstruction_name", file_id, string_type);
#elif RECONSTRUCTION == MP5
  hdf5_write_single_val("MP5", "reconstruction_name", file_id, string_type);
#endif

  int n1 = N1TOT, n2 = N2TOT, n3 = N3TOT;
  hdf5_write_single_val(&n1, "n1", file_id, H5T_NATIVE_INT);
  hdf5_write_single_val(&n2, "n2", file_id, H5T_NATIVE_INT);
  hdf5_write_single_val(&n3, "n3", file_id, H5T_NATIVE_INT);

  int n_prims = NVAR;
  hdf5_write_single_val(&n_prims, "n_prims", file_id, H5T_NATIVE_INT);
  // In case we do passive variables
  int n_prims_passive = 0;
  hdf5_write_single_val(&n_prims_passive, "n_prims_passive", file_id, H5T_NATIVE_INT);
  hdf5_write_str_list(varNames, "prim_names", file_id, HDF_STR_LEN, n_prims);

  hdf5_write_single_val(&gam, "gam", file_id, H5T_NATIVE_DOUBLE);
  #if ELECTRONS
  hdf5_write_single_val(&game, "gam_e", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&gamp, "gam_p", file_id, H5T_NATIVE_DOUBLE);
  #endif
  hdf5_write_single_val(&cour, "cour", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&tf, "tf", file_id, H5T_NATIVE_DOUBLE);
  hdf5_add_units("tf", "code", file_id);

  //Geometry
  hdf5_make_directory("geom", file_id);
  hdf5_set_directory("/header/geom/");

  hdf5_write_single_val(&(startx[1]), "startx1", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&(startx[2]), "startx2", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&(startx[3]), "startx3", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&(dx[1]), "dx1", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&(dx[2]), "dx2", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&(dx[3]), "dx3", file_id, H5T_NATIVE_DOUBLE);
  int n_dim = NDIM;
  hdf5_write_single_val(&n_dim, "n_dim", file_id, H5T_NATIVE_INT);
  #if METRIC == MKS
  hdf5_make_directory("mks", file_id);
  hdf5_set_directory("/header/geom/mks/");
  hdf5_write_single_val(&Rin, "Rin", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&Rout, "Rout", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&Rhor, "Reh", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&hslope, "hslope", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&a, "a", file_id, H5T_NATIVE_DOUBLE);
  #if POLYTH
  hdf5_write_single_val(&poly_xt, "poly_xt", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&poly_alpha, "poly_alpha", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&mks_smooth, "mks_smooth", file_id, H5T_NATIVE_DOUBLE);
  #endif
  #endif

  hdf5_set_directory("/");

  int is_full_dump = 1; // TODO do full dumps
  hdf5_write_single_val(&is_full_dump, "is_full_dump", file_id, H5T_NATIVE_INT);
  hdf5_write_single_val(&t, "t", file_id, H5T_NATIVE_DOUBLE);
  hdf5_add_units("t", "code", file_id);
  hdf5_write_single_val(&dt, "dt", file_id, H5T_NATIVE_DOUBLE);
  hdf5_add_units("dt", "code", file_id);
  hdf5_write_single_val(&nstep, "n_step", file_id, H5T_NATIVE_INT);
  hdf5_write_single_val(&dump_cnt, "n_dump", file_id, H5T_NATIVE_INT);

  hdf5_write_single_val(&DTd, "dump_cadence", file_id, H5T_NATIVE_DOUBLE);
  //hdf5_write_single_val(&DTf, "full_dump_cadence", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_single_val(&failed, "failed", file_id, H5T_NATIVE_INT);

  // Write primitive variables
  pack_vector_float(S->P, data, NVAR);
  hdf5_write_vector(data, "prims", file_id, NVAR, OUT_H5_TYPE);
  hdf5_add_units("prims", "code", file_id);

  // Write jcon (not recoverable from prims)
  pack_vector_float(S->jcon, data, NDIM);
  hdf5_write_vector(data, "jcon", file_id, NDIM, OUT_H5_TYPE);

  // Write divB and fail diagnostics
  int ind = 0;
  ZLOOP_OUT {
    data[ind] = flux_ct_divb(G, S, i, j, k);
    ind++;
  }
  hdf5_write_scalar(data, "divB", file_id, OUT_H5_TYPE);

  pack_scalar_int(fail_save, idata);
  ZLOOP fail_save[k][j][i] = 0;
  hdf5_write_scalar(idata, "fail", file_id, H5T_NATIVE_INT);

  // Write some extras
  hdf5_make_directory("extras", file_id);
  hdf5_set_directory("/extras/");

  ind = 0;
  ZLOOP_OUT {
    get_state(G, S, i, j, k, CENT);
    data[ind] = bsq_calc(S, i, j, k);
    ind++;
  }
  hdf5_write_scalar(data, "bsq", file_id, OUT_H5_TYPE);

  omega_calc(G, S, omega);
  pack_scalar_float(*omega, data);
  hdf5_write_scalar(data, "omega", file_id, OUT_H5_TYPE);

  // These are far too much space to continue
  pack_vector_float(S->bcon, data, NDIM);
  hdf5_write_vector(data, "bcon", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->bcov, data, NDIM);
  hdf5_write_vector(data, "bcov", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->ucon, data, NDIM);
  hdf5_write_vector(data, "ucon", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->ucov, data, NDIM);
  hdf5_write_vector(data, "ucov", file_id, NDIM, OUT_H5_TYPE);

  hdf5_close(file_id);

  timer_stop(TIMER_IO);
}

#define NGRIDVARS 11
void dump_grid(struct GridGeom *G)
{
  OUT_TYPE *x[NGRIDVARS];
  for (int d = 0; d < NGRIDVARS; d++) x[d] = calloc(N1*N2*N3,sizeof(OUT_TYPE));
  const char *coordNames[] = {"X", "Y", "Z", "r", "th", "phi", "X1", "X2", "X3", "gdet", "lapse"};

  OUT_TYPE *g = calloc(N1*N2*N3*NDIM*NDIM,sizeof(OUT_TYPE));

  hid_t file_id = hdf5_open("dumps/grid.h5");

  // Batch fill grid var buffers since lots of them are the same
  int ind = 0;
  ZLOOP_OUT { //Change if reinstating full grid out
	double xp[4];
    coord(i, j, k, CENT, xp);
    #if METRIC == MINKOWSKI
    x[0][ind] = xp[1];
    x[1][ind] = xp[2];
    x[2][ind] = xp[3];
    x[3][ind] = 0;
    x[4][ind] = 0;
    x[5][ind] = 0;
    #elif METRIC == MKS
    double r, th;
    bl_coord(xp, &r, &th);
    x[0][ind] = r*cos(xp[3])*sin(th);
    x[1][ind] = r*sin(xp[3])*sin(th);
    x[2][ind] = r*cos(th);
    x[3][ind] = r;
    x[4][ind] = th;
    x[5][ind] = xp[3];
    #endif

    x[6][ind] = xp[1];
    x[7][ind] = xp[2];
    x[8][ind] = xp[3];

    x[9][ind] = G->gdet[CENT][j][i];
    x[10][ind] = G->lapse[CENT][j][i];

    ind++;
  }

  for (int d = 0; d < NGRIDVARS; d++){
    hdf5_write_scalar(x[d], coordNames[d], file_id, OUT_H5_TYPE);
  }

  pack_Gtensor_float(G->gcon[CENT], g);
  hdf5_write_tensor(g, "gcon", file_id, NDIM, NDIM, OUT_H5_TYPE);

  pack_Gtensor_float(G->gcov[CENT], g);
  hdf5_write_tensor(g, "gcov", file_id, NDIM, NDIM, OUT_H5_TYPE);

  for (int d = 0; d < NGRIDVARS; d++) free(x[d]);

  hdf5_close(file_id);
}
