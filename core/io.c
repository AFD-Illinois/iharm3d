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
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3",
                            "KEL", "KTOT"};
  #else
  //const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
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

  // Set our string type size.  None of our names are long so this should suffice
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, 20);

  // Write header
  hdf5_write_header_val(VERSION, "version_string", file_id, string_type);
  #if METRIC == MINKOWSKI
  hdf5_write_header_val("MINKOWSKI", "metric_name", file_id, string_type);
  #elif METRIC == MKS
  hdf5_write_header_val("MKS", "metric_name", file_id, string_type);
  #endif
  int has_electrons = ELECTRONS;
  hdf5_write_header_val(&has_electrons, "has_electrons", file_id, H5T_NATIVE_INT);
  hdf5_write_header_val(&t, "t", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&tf, "tf", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&nstep, "n_step", file_id, H5T_NATIVE_INT);
  int n1 = N1TOT, n2 = N2TOT, n3 = N3TOT;
  hdf5_write_header_val(&n1, "n1", file_id, H5T_NATIVE_INT);
  hdf5_write_header_val(&n2, "n2", file_id, H5T_NATIVE_INT);
  hdf5_write_header_val(&n3, "n3", file_id, H5T_NATIVE_INT);

  hdf5_write_header_val(&gam, "gam", file_id, H5T_NATIVE_DOUBLE);
  #if ELECTRONS
  hdf5_write_header_val(&game, "gam_e", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&gamp, "gam_p", file_id, H5T_NATIVE_DOUBLE);
  #endif
  hdf5_write_header_val(&cour, "cour", file_id, H5T_NATIVE_DOUBLE);
//  hdf5_write_header_val(&DTd, "DTd", file_id, H5T_NATIVE_DOUBLE);
//  hdf5_write_header_val(&DTl, "DTl", file_id, H5T_NATIVE_DOUBLE);
//  hdf5_write_header_val(&DTr, "DTr", file_id, H5T_NATIVE_INT);
//  hdf5_write_header_val(&DTp, "DTp", file_id, H5T_NATIVE_INT);
  hdf5_write_header_val(&dump_cnt, "dump_cnt", file_id, H5T_NATIVE_INT);
  hdf5_write_header_val(&dt, "dt", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&failed, "failed", file_id, H5T_NATIVE_INT);

  //Geometry
  hdf5_write_header_val(&(startx[1]), "startx1", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&(startx[2]), "startx2", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&(startx[3]), "startx3", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&(dx[1]), "dx1", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&(dx[2]), "dx2", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&(dx[3]), "dx3", file_id, H5T_NATIVE_DOUBLE);
  #if METRIC == MKS
  hdf5_write_header_val(&Rin, "Rin", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&Rout, "Rout", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&Rhor, "Reh", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&hslope, "hslope", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&a, "a", file_id, H5T_NATIVE_DOUBLE);
  #if POLYTH
  hdf5_write_header_val(&poly_xt, "poly_xt", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&poly_alpha, "poly_alpha", file_id, H5T_NATIVE_DOUBLE);
  hdf5_write_header_val(&mks_smooth, "mks_smooth", file_id, H5T_NATIVE_DOUBLE);
  #endif
  #endif

  // Write primitive variables
  pack_vector_float(S->P, data, NVAR);
  hdf5_write_vector(data, "prims", file_id, NVAR, OUT_H5_TYPE);

  // Write derived scalars from functions
  int ind = 0;
  ZLOOP_OUT {
    get_state(G, S, i, j, k, CENT);
    data[ind] = bsq_calc(S, i, j, k);
    ind++;
  }
  hdf5_write_scalar(data, "bsq", file_id, OUT_H5_TYPE);

  omega_calc(G, S, omega);
  pack_scalar_float(*omega, data);
  hdf5_write_scalar(data, "omega", file_id, OUT_H5_TYPE);

  ind = 0;
  ZLOOP_OUT {
    data[ind] = mhd_gamma_calc(G, S, i, j, k, CENT);
    ind++;
  }
  hdf5_write_scalar(data, "gamma", file_id, OUT_H5_TYPE);

  // TODO be able to output divb in double separately?
  // Cadence divb/fail differently from other vars?
  ind = 0;
  ZLOOP_OUT {
    data[ind] = flux_ct_divb(G, S, i, j, k);
    ind++;
  }
  hdf5_write_scalar(data, "divb", file_id, OUT_H5_TYPE);

  pack_scalar_int(fail_save, idata);
  ZLOOP_OUT fail_save[k][j][i] = 0;
  hdf5_write_scalar(idata, "fail", file_id, H5T_NATIVE_INT);

  // Write vector quantities
  pack_vector_float(S->bcon, data, NDIM);
  hdf5_write_vector(data, "bcon", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->bcov, data, NDIM);
  hdf5_write_vector(data, "bcov", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->ucon, data, NDIM);
  hdf5_write_vector(data, "ucon", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->ucov, data, NDIM);
  hdf5_write_vector(data, "ucov", file_id, NDIM, OUT_H5_TYPE);

  pack_vector_float(S->jcon, data, NDIM);
  hdf5_write_vector(data, "jcon", file_id, NDIM, OUT_H5_TYPE);

  hdf5_close(file_id);

  timer_stop(TIMER_IO);
}

// TODO delete this when possible
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
