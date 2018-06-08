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

static int restart_id = 0;

void restart_write(struct FluidState *S)
{
  timer_start(TIMER_RESTART);
  hsize_t file_start[4],mem_start[4],one,zero;
  hsize_t fdims[4] = {NVAR, N3TOT, N2TOT, N1TOT};
  hsize_t mem_copy_dims[4] = {NVAR, N3, N2, N1};
  hsize_t mem_full_dims[4] = {NVAR, N3 + 2*NG, N2 + 2*NG, N1 + 2*NG};

  // Pass by reference, baby
  one = 1;
  zero = 0;

  // Keep track of our own index
  restart_id++;

  char fname[2048];
  sprintf(fname, "restarts/restart_%08d.h5", restart_id);

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // Write header
  hid_t filespace = H5Screate_simple(1, &one, NULL);
  hsize_t single_var_count = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &zero, NULL, &single_var_count,
    NULL);
  hid_t memspace  = H5Screate_simple(1, &single_var_count, NULL);

  // TODO add version
  add_dbl_value(t, "t", file_id, filespace, memspace);
  add_int_value(nstep, "nstep", file_id, filespace, memspace);
  add_dbl_value(tf, "tf", file_id, filespace, memspace);
  add_dbl_value(a, "a", file_id, filespace, memspace);
  add_dbl_value(gam, "gam", file_id, filespace, memspace);
  #if ELECTRONS
  add_dbl_value(game, "game", file_id, filespace, memspace);
  add_dbl_value(gamp, "gamp", file_id, filespace, memspace);
  add_dbl_value(fel0, "fel0", file_id, filespace, memspace);
  #endif
  add_dbl_value(cour, "cour", file_id, filespace, memspace);
  add_dbl_value(DTd, "DTd", file_id, filespace, memspace);
  add_dbl_value(DTl, "DTl", file_id, filespace, memspace);
  add_int_value(DTr, "DTr", file_id, filespace, memspace);
  add_int_value(DTp, "DTp", file_id, filespace, memspace);
  add_int_value(restart_id, "restart_id", file_id, filespace, memspace);
  add_dbl_value(dt, "dt", file_id, filespace, memspace);
  add_int_value(lim, "lim", file_id, filespace, memspace);
  add_int_value(failed, "failed", file_id, filespace, memspace);
  add_dbl_value(Rin, "Rin", file_id, filespace, memspace);
  add_dbl_value(Rout, "Rout", file_id, filespace, memspace);
  add_dbl_value(hslope, "hslope", file_id, filespace, memspace);
  add_dbl_value(R0, "R0", file_id, filespace, memspace);
  add_dbl_value(Rhor, "Rhor", file_id, filespace, memspace);
  add_dbl_value(Risco, "Risco", file_id, filespace, memspace);
  add_dbl_value(tdump, "tdump", file_id, filespace, memspace);
  add_dbl_value(tlog, "tlog", file_id, filespace, memspace);

  H5Sclose(filespace);
  H5Sclose(memspace);

  // Write data
  filespace = H5Screate_simple(4, fdims, NULL);
  for (int d = 0; d < 4; d++) {
    file_start[d] = (d == 0) ? 0 : global_start[3-d];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, mem_copy_dims, NULL);

  memspace = H5Screate_simple(4, mem_full_dims, NULL);

  for (int d = 0; d < 4; d++)
    mem_start[d] = (d == 0) ? 0 : NG;

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, mem_copy_dims,
    NULL);

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, "p", H5T_NATIVE_DOUBLE, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id,
      &(S->P[0][0][0][0]));
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(memspace);
  H5Sclose(filespace);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);
  if(mpi_io_proc()) {
    fprintf(stdout, "RESTART %s\n", fname);

    // Symlink when we're done writing (link to last good file)
    char fname_nofolder[80];
    sprintf(fname_nofolder, "restart_%08d.h5", restart_id);

    //Chained OS functions: switch to restart directory,
    // remove current link, link last file, switch back
    int errcode;
    errcode = chdir("restarts");
    if ( access("restart.last", F_OK) != -1 ) {
      errcode = errcode || remove("restart.last");
    }
    errcode = errcode || symlink(fname_nofolder, "restart.last");
    errcode = errcode || chdir("..");
    if(errcode != 0) {
      printf("Symlink failed: errno %d\n", errno);
      exit(-1);
    }
  }
  timer_stop(TIMER_RESTART);
}

void restart_read(char *fname, struct FluidState *S)
{
  hsize_t file_start[4],mem_start[4];
  hsize_t zero,one;
  hsize_t fdims[4] = {NVAR, N3TOT, N2TOT, N1TOT};
  hsize_t mem_copy_dims[4] = {NVAR, N3, N2, N1};
  hsize_t mem_full_dims[4] = {NVAR, N3 + 2*NG, N2 + 2*NG, N1 + 2*NG};

  if(mpi_io_proc()) {
    fprintf(stderr, "Restarting from %s\n\n", fname);
  }

  zero = 0;
  one = 1;

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  // Read header
  hid_t filespace = H5Screate_simple(1, &one, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &zero, NULL, &one,
    NULL);
  hid_t memspace  = H5Screate_simple(1, &one, NULL);

  get_dbl_value(&t, "t", file_id, filespace, memspace);
  get_int_value(&nstep, "nstep", file_id, filespace, memspace);
  get_int_value(&restart_id, "restart_id", file_id, filespace, memspace);
  get_dbl_value(&tf, "tf", file_id, filespace, memspace);
  get_dbl_value(&a, "a", file_id, filespace, memspace);
  get_dbl_value(&gam, "gam", file_id, filespace, memspace);
  #if ELECTRONS
  get_dbl_value(&game, "game", file_id, filespace, memspace);
  get_dbl_value(&gamp, "gamp", file_id, filespace, memspace);
  get_dbl_value(&fel0, "fel0", file_id, filespace, memspace);
  #endif
  get_dbl_value(&cour, "cour", file_id, filespace, memspace);
  get_dbl_value(&DTd, "DTd", file_id, filespace, memspace);
  get_dbl_value(&DTl, "DTl", file_id, filespace, memspace);
  get_int_value(&DTr, "DTr", file_id, filespace, memspace);
  get_int_value(&DTp, "DTp", file_id, filespace, memspace);
  get_int_value(&restart_id, "restart_id", file_id, filespace, memspace);
  get_dbl_value(&dt, "dt", file_id, filespace, memspace);
  get_int_value(&lim, "lim", file_id, filespace, memspace);
  get_int_value(&failed, "failed", file_id, filespace, memspace);
  get_dbl_value(&Rin, "Rin", file_id, filespace, memspace);
  get_dbl_value(&Rout, "Rout", file_id, filespace, memspace);
  get_dbl_value(&hslope, "hslope", file_id, filespace, memspace);
  get_dbl_value(&R0, "R0", file_id, filespace, memspace);
  get_dbl_value(&Rhor, "Rhor", file_id, filespace, memspace);
  get_dbl_value(&Risco, "Risco", file_id, filespace, memspace);
  get_dbl_value(&tdump, "tdump", file_id, filespace, memspace);
  get_dbl_value(&tlog, "tlog", file_id, filespace, memspace);

  H5Sclose(filespace);
  H5Sclose(memspace);

  // Read data
  filespace = H5Screate_simple(4, fdims, NULL);

  for (int d = 0; d < 4; d++) {
    file_start[d] = (d == 0) ? 0 : global_start[3-d];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, mem_copy_dims,
    NULL);

  memspace = H5Screate_simple(4, mem_full_dims, NULL);

  for (int d = 0; d < 4; d++)
    mem_start[d] = (d == 0) ? 0 : NG;

  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, mem_copy_dims,
    NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dset_id = H5Dopen(file_id, "p", plist_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id,
	  &(S->P[0][0][0][0]));
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  mpi_barrier();
}

int restart_init(struct GridGeom *G, struct FluidState *S)
{
  char fname[512];
  sprintf(fname, "restarts/restart.last");

  FILE *fp = fopen(fname,"rb");
  if (fp == NULL) {
    if (mpi_io_proc())
      fprintf(stdout, "No restart file: %d\n\n", errno);
    return 0;
  }
  fclose(fp);

  if (mpi_io_proc())
    fprintf(stdout, "Loading restart file %s\n\n", fname);
  zero_arrays();

  restart_read(fname, S);

  set_grid(G);

  get_state_vec(G, S, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
  prim_to_flux_vec(G, S, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, S->U);

  set_bounds(G, S);

  return 1;
}

void add_int_value(int val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace)
{
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &val);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
}

void add_dbl_value(double val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace)
{
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &val);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
}

void add_str_value(const char* val, const char *name, hid_t file_id,
  hid_t filespace, hid_t memspace)
{
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen(val));
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, name, string_type, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Dwrite(dset_id, string_type, memspace, filespace, plist_id, val);
  H5Dclose(dset_id);
  H5Pclose(plist_id);
}

void get_int_value(int *val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace)
{
  hid_t plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dset_id = H5Dopen(file_id, name, plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, val);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  mpi_int_broadcast(val);
}

void get_dbl_value(double *val, const char *name, hid_t file_id,
  hid_t filespace, hid_t memspace)
{
  hid_t plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dset_id = H5Dopen(file_id, name, plist_id);
  H5Pclose(plist_id);
  plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, val);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  mpi_dbl_broadcast(val);
}

