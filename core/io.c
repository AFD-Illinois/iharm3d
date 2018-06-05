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

static int dump_id = 0;

void init_io()
{
}

void dump(struct GridGeom *G, struct FluidState *S)
{
  timer_start(TIMER_IO);

  static int firstc = 1;
  static float *data; // TODO option to output double
  static int *idata;
  char fname[80];
  FILE *xml = NULL;
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT};
  hsize_t fdims_vec[] = {N1TOT, N2TOT, N3TOT, NDIM};
  hsize_t mem_copy_dims[] = {N1, N2, N3};
  hsize_t mem_copy_dims_vec[] = {N1, N2, N3, NDIM};
  #if ELECTRONS
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3",
                            "KEL", "KTOT"};
  #else
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
  #endif

  hid_t file_id, plist_id;//, dset_id;
  hsize_t mem_start[3], file_start[3], one, zero;
  hsize_t mem_start_vec[4], file_start_vec[4];

  one = 1;
  zero = 0;

  if(firstc) {
    dump_grid(G);
    data = calloc(N1*N2*N3*NVAR,sizeof(float));
    idata = calloc(N1*N2*N3,sizeof(int));
    if(data == NULL) {
      fprintf(stderr,"failed to allocate data in dump\n");
      exit(8);
    }
    firstc = 0;
  }

  if(mpi_io_proc()) {
    sprintf(fname, "dumps/dump_%08d.xmf", dump_id); // TODO define filenames elsewhere?
    xml = write_xml_head(fname,t);
  }

  sprintf(fname, "dumps/dump_%08d.h5", dump_id);
  if(mpi_io_proc()) {
    fprintf(stdout, "DUMP %s\n", fname);
  }

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // Write header
  hid_t filespace = H5Screate_simple(1, &one, NULL);
  hsize_t single_var_count = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &zero, NULL, &single_var_count,
    NULL);
  hid_t memspace  = H5Screate_simple(1, &single_var_count, NULL);
  add_str_value(VERSION, "VERSION", file_id, filespace, memspace);
  #if METRIC == MINKOWSKI
  add_str_value("MINKOWSKI", "METRIC", file_id, filespace, memspace);
  #elif METRIC == MKS
  add_str_value("MKS", "METRIC", file_id, filespace, memspace);
  #endif
  add_int_value(ELECTRONS, "ELECTRONS", file_id, filespace, memspace);
  add_dbl_value(t, "t", file_id, filespace, memspace);
  add_dbl_value(tf, "tf", file_id, filespace, memspace);
  add_int_value(nstep, "nstep", file_id, filespace, memspace);
  add_int_value(N1TOT, "N1", file_id, filespace, memspace);
  add_int_value(N2TOT, "N2", file_id, filespace, memspace);
  add_int_value(N3TOT, "N3", file_id, filespace, memspace);
  add_dbl_value(startx[1], "startx1", file_id, filespace, memspace);
  add_dbl_value(startx[2], "startx2", file_id, filespace, memspace);
  add_dbl_value(startx[3], "startx3", file_id, filespace, memspace);
  add_dbl_value(dx[1], "dx1", file_id, filespace, memspace);
  add_dbl_value(dx[2], "dx2", file_id, filespace, memspace);
  add_dbl_value(dx[3], "dx3", file_id, filespace, memspace);
  #if METRIC == MKS
  add_dbl_value(Rin, "Rin", file_id, filespace, memspace);
  add_dbl_value(Rout, "Rout", file_id, filespace, memspace);
  add_dbl_value(Rhor, "Reh", file_id, filespace, memspace);
  add_dbl_value(hslope, "hslope", file_id, filespace, memspace);
  add_dbl_value(a, "a", file_id, filespace, memspace);
  #if POLYTH
  add_dbl_value(poly_xt, "poly_xt", file_id, filespace, memspace);
  add_dbl_value(poly_alpha, "poly_alpha", file_id, filespace, memspace);
  add_dbl_value(mks_smooth, "mks_smooth", file_id, filespace, memspace);
  #endif
  #endif
  add_dbl_value(gam, "gam", file_id, filespace, memspace);
  #if ELECTRONS
  add_dbl_value(game, "game", file_id, filespace, memspace);
  add_dbl_value(gamp, "gamp", file_id, filespace, memspace);
  #endif
  add_dbl_value(cour, "cour", file_id, filespace, memspace);
  add_dbl_value(DTd, "DTd", file_id, filespace, memspace);
  add_dbl_value(DTl, "DTl", file_id, filespace, memspace);
  add_int_value(DTr, "DTr", file_id, filespace, memspace);
  add_int_value(DTp, "DTp", file_id, filespace, memspace);
  add_int_value(dump_cnt, "dump_cnt", file_id, filespace, memspace);
  add_dbl_value(dt, "dt", file_id, filespace, memspace);
  add_int_value(failed, "failed", file_id, filespace, memspace);

  H5Sclose(filespace);
  H5Sclose(memspace);

  // Tell HDF how data is laid out in memory and in file
  filespace = H5Screate_simple(3, fdims, NULL);
  for (int d = 0; d < 3; d++) {
    file_start[d] = global_start[d];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, mem_copy_dims,
    NULL);

  memspace = H5Screate_simple(3, mem_copy_dims, NULL);
  for (int d = 0; d < 3; d++)
    mem_start[d] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, mem_copy_dims,
    NULL);


//  printf("Size of memspace: %lld x %lld x %lld\n",
//       mem_copy_dims[0], mem_copy_dims[1], mem_copy_dims[2]);
//  printf("Location of memspace: %lld %lld %lld\n",
//	file_start[0], file_start[1], file_start[2]);
//  printf("Size of filespace: %lld x %lld x %lld\n",
//       fdims[0], fdims[1], fdims[2]);

  // Write primitive variables
  PLOOP {
    int ind = 0;
    ZLOOP_OUT {
      data[ind] = (float)(S->P[ip][k][j][i]);
      ind++;
    }
    write_scalar(data, file_id, varNames[ip], filespace, memspace, xml);
  }

  // Write other scalars
  int ind = 0;
  ZLOOP_OUT {
    get_state(G, S, i, j, k, CENT);
    double bsq = bsq_calc(S, i, j, k);
    data[ind] = bsq;
    ind++;
  }
  write_scalar(data, file_id, "bsq", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    double gamma = 1.;
    mhd_gamma_calc(G, S, i, j, k, CENT, &gamma);
    data[ind] = gamma;
    ind++;
  }
  write_scalar(data, file_id, "gamma", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    double divb = flux_ct_divb(G, S, i, j, k);
    data[ind] = divb;
    ind++;
  }
  write_scalar(data, file_id, "divb", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    idata[ind] = fail_save[k][j][i];
    fail_save[k][j][i] = 0;
    ind++;
  }
  write_scalar_int(idata, file_id, "fail", filespace, memspace, xml);

  H5Sclose(filespace);
  H5Sclose(memspace);

  // Write vector quantities
  filespace = H5Screate_simple(4, fdims_vec, NULL);
  for (int d = 0; d < 3; d++)
    file_start_vec[d] = global_start[d];
  file_start_vec[3] = 0;

  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start_vec, NULL, mem_copy_dims_vec,
    NULL);

  memspace = H5Screate_simple(4, mem_copy_dims_vec, NULL);
  for (int d = 0; d < 4; d++)
    mem_start_vec[d] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start_vec, NULL, mem_copy_dims_vec,
    NULL);

  ind = 0;
  ZLOOP_OUT {
    DLOOP1 {
      data[ind] = S->bcon[mu][k][j][i];
      ind++;
    }
  }
  write_vector(data, file_id, "bcon", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    DLOOP1 {
      data[ind] = S->bcov[mu][k][j][i];
      ind++;
    }
  }
  write_vector(data, file_id, "bcov", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    DLOOP1 {
      data[ind] = S->ucon[mu][k][j][i];
      ind++;
    }
  }
  write_vector(data, file_id, "ucon", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    DLOOP1 {
      data[ind] = S->ucov[mu][k][j][i];
      ind++;
    }
  }
  write_vector(data, file_id, "ucov", filespace, memspace, xml);

  ind = 0;
  ZLOOP_OUT {
    DLOOP1 {
      data[ind] = S->jcon[mu][k][j][i];
      ind++;
    }
  }
  write_vector(data, file_id, "jcon", filespace, memspace, xml);

  // Close file
  if(mpi_io_proc()) {
    write_xml_closing(xml);
  }

  H5Sclose(memspace);
  H5Sclose(filespace);
  // TODO
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  dump_id++;

  timer_stop(TIMER_IO);
}

// TODO delete this when possible
#define NGRIDVARS 11
void dump_grid(struct GridGeom *G)
{
  float *x[NGRIDVARS];
  double xp[4];
  hid_t dset_id;
  hsize_t file_start[3], file_count[3];
  hsize_t mem_start[] = {0, 0, 0};
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT};
  hsize_t dims[] = {N1, N2, N3};

  const char *coordNames[] = {"X", "Y", "Z", "r", "th", "phi", "X1", "X2", "X3", "gdet", "lapse"};

  // This was for extra corner zones at the end.  We don't need them since we want center values
  // TODO ... right?
//  for (int d = 0; d < 3; d++) {
//    if (global_stop[d] == fdims[d]-1)
//      dims[d]++;
//  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate("dumps/grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  hid_t filespace = H5Screate_simple(3, fdims, NULL);
  for (int d = 0; d < 3; d++) {
    file_start[d] = global_start[d];
    file_count[d] = dims[d];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

  hid_t memspace = H5Screate_simple(3, dims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

  int total_size = dims[0]*dims[1]*dims[2];
  for(int d = 0; d < NGRIDVARS; d++) {
    x[d] = calloc(total_size,sizeof(float));
    if(x[d] == NULL) {
      fprintf(stderr,"Failed to allocate x[d] in dump_grid\n");
      exit(5);
    }
  }

  int ind = 0;
  ZLOOP_OUT { //Change if reinstating full grid out
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
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    dset_id = H5Dcreate(file_id, coordNames[d], H5T_NATIVE_FLOAT, filespace,
      H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x[d]);
    H5Dclose(dset_id);
    H5Pclose(plist_id);
  }

  H5Sclose(filespace);
  H5Sclose(memspace);

  hsize_t file_start_tens[5];
  hsize_t mem_start_tens[5] = {0, 0, 0, 0, 0};
  hsize_t fdims_tens[] = {N1TOT, N2TOT, N3TOT, NDIM, NDIM};
  hsize_t dims_tens[] = {N1, N2, N3, NDIM, NDIM};

  filespace = H5Screate_simple(5, fdims_tens, NULL);
  for (int d = 0; d < 5; d++) {
    file_start_tens[d] = (d < 3) ? global_start[d] : 0;
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start_tens, NULL, dims_tens, NULL);

  memspace = H5Screate_simple(5, dims_tens, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start_tens, NULL, dims_tens, NULL);

  float *g = calloc(N1*N2*N3*NDIM*NDIM,sizeof(float));
  if(g == NULL) {
    fprintf(stderr,"Failed to allocate x[d] in dump_grid\n");
    exit(5);
  }

  ind = 0;
  ZLOOP_OUT {
    DLOOP2 {
      g[ind] = G->gcon[CENT][mu][nu][j][i];
    }
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "gcon", H5T_NATIVE_FLOAT, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, g);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  ind = 0;
  ZLOOP_OUT {
    DLOOP2 {
      g[ind] = G->gcon[CENT][mu][nu][j][i];
    }
  }

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, "gcov", H5T_NATIVE_FLOAT, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, g);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(filespace);
  H5Sclose(memspace);

  for (int d = 0; d < NGRIDVARS; d++) free(x[d]);

  H5Fclose(file_id);
}
