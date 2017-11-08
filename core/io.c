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

void write_scalar(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml);
void write_scalar_int(int *data, hid_t file_id, const char *name,
  hid_t filespace, hid_t memspace, FILE *xml);
void add_int_value(int val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace);
void add_dbl_value(double val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace);
void add_str_value(const char* val, const char *name, hid_t file_id, 
  hid_t filespace, hid_t memspace);

static int dump_id = 0;
//static int restart_id = 0;

void dump(struct GridGeom *G, struct FluidState *S)
{
  static int firstc = 1;
  static float *data;
  static int *idata;
  char name[80];
  FILE *xml = NULL;
  hsize_t fdims[] = {N3TOT, N2TOT, N1TOT};
  hsize_t mdims[] = {N3, N2, N1};
  #if ELECTRONS
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3", 
                            "KEL", "KTOT"};
  #else
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
  #endif
  //struct of_state q;

  hid_t file_id, filespace, memspace, plist_id;//, dset_id;
  hsize_t mem_start[3], file_start[3], one, zero;
  hsize_t file_count[3], file_grid_dims[3];

  one = 1;
  zero = 0;

  if(firstc) {
    mkdir("dumps", 0777);
    dump_grid();
    data = malloc(N1*N2*N3*NVAR*sizeof(float));
    idata = malloc(N1*N2*N3*sizeof(int));
    if(data == NULL) {
      fprintf(stderr,"failed to allocate data in dump\n");
      exit(8);
    }
    firstc = 0;
  }

  if(mpi_io_proc()) {
    xml = write_xml_head(dump_id,t);
  }

  char fname[2048];
  sprintf(fname, "dump_%08d.h5", dump_id);
  strcpy(name, dumpdir);
  strcat(name, fname);
  if(mpi_io_proc()) {
    fprintf(stdout, "DUMP %s\n", name); 
  }

  plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  
  // Write header
  file_grid_dims[0] = 1;
  filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, 
    NULL);
  memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);
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
  for (int d = 0; d < 3; d++)
    file_grid_dims[d] = fdims[d];
  filespace = H5Screate_simple(3, file_grid_dims, NULL);
  for (int d = 0; d < 3; d++) {
    file_start[d] = global_start[d];
    file_count[d] = mdims[d];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, 
    NULL);

  memspace = H5Screate_simple(3, file_count, NULL);
  for (int d = 0; d < 3; d++)
    mem_start[d] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, 
    NULL);

  // Write primitive variables
  PLOOP {
    int ind = 0;
    ZLOOP {
      //data[ind] = (float)P[i][j][k][ip]; 
      data[ind] = (float)(S->P[ip][k][j][i]);
      ind++;
    }
    write_scalar(data, file_id, varNames[ip], filespace, memspace, xml);
  }

  // Write everything else
  double X[NDIM];
  int ind = 0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    data[ind] = X[1];
    ind++;
  }
  write_scalar(data, file_id, "X1", filespace, memspace, xml);
  ind = 0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    data[ind] = X[2];
    ind++;
  }
  write_scalar(data, file_id, "X2", filespace, memspace, xml);
  ind = 0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    data[ind] = X[3];
    ind++;
  }
  write_scalar(data, file_id, "X3", filespace, memspace, xml);
  
  #if METRIC == MKS
  double r, theta, phi;
  ind = 0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &theta);
    phi = X[3];
    data[ind] = r;
    ind++;
  }
  write_scalar(data, file_id, "r", filespace, memspace, xml);
  ind = 0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &theta);
    phi = X[3];
    data[ind] = theta;
    ind++;
  }
  write_scalar(data, file_id, "theta", filespace, memspace, xml);
  ind = 0;
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &theta);
    phi = X[3];
    data[ind] = phi;
    ind++;
  }
  write_scalar(data, file_id, "phi", filespace, memspace, xml);
  #endif
  
  ind = 0;
  ZLOOP {
    get_state(G, S, i, j, k, CENT);
    double bsq = dot_grid(S->bcon, S->bcov, i, j, k);
    //struct of_geom *geom = get_geometry(i, j, CENT) ;
    //get_state(P[i][j][k], geom, &q);
    //double bsq = 0.;
    //for (int d = 0; d < NDIM; d++) bsq += q.bcon[d]*q.bcov[d];
    data[ind] = bsq;
    ind++;
  }
  write_scalar(data, file_id, "bsq", filespace, memspace, xml);

  ind = 0;
  ZLOOP {
    //struct of_geom *geom = get_geometry(i, j, CENT) ;
    double gamma = 1.;
    mhd_gamma_calc(G, S, i, j, k, CENT, &gamma);
    data[ind] = gamma;
    ind++;
  }
  write_scalar(data, file_id, "gamma", filespace, memspace, xml);

  ind = 0;
  ZLOOP {
    double divb = flux_ct_divb(G, S, i, j, k);
    data[ind] = divb;
    ind++;
  }
  write_scalar(data, file_id, "divb", filespace, memspace, xml);

  ind = 0;
  ZLOOP {
    idata[ind] = fail_save[k][j][i];
    fail_save[k][j][i] = 0;
    ind++;
  }
  write_scalar_int(idata, file_id, "fail", filespace, memspace, xml);

  // Close file
  if(mpi_io_proc()) {
    write_xml_closing(xml);
  }

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  dump_id++;
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

void restart_write(struct FluidState *S)
{
  /*
  int fdims[3] = {N1TOT, N2TOT, N3TOT};
  int mdims[3] = {N1, N2, N3};
  hsize_t file_grid_dims[4],file_start[4],file_count[4];
  hsize_t mem_grid_dims[4],mem_start[4],one,zero;

  char name[80];

  one = 1;
  zero = 0;
  mkdir("restarts", 0777);

  char fname[2048];
  sprintf(fname, "restart_%08d.h5", restart_id);
  strcpy(name, restartdir);
  strcat(name, fname);
 
  char lastname[2048];
  sprintf(fname, "restart.last");
  strcpy(lastname, restartdir);
  strcat(lastname, fname);

  restart_id++;
  FILE *fp = fopen(lastname,"w");
  fprintf(fp,"%s\n", name);
  fclose(fp);
  if(mpi_io_proc()) {
    fprintf(stdout, "RESTART %s\n", name);
  }

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // Write header
  file_grid_dims[0] = 1;
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count,
    NULL);
  hid_t memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

  add_dbl_value(t, "t", file_id, filespace, memspace);
  add_int_value(nstep, "nstep", file_id, filespace, memspace);
  add_int_value(dump_id, "dump_id", file_id, filespace, memspace);
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
  for (int d = 0; d < 3; d++)
    file_grid_dims[d] = fdims[d];
  file_grid_dims[3] = NVAR;
  filespace = H5Screate_simple(4, file_grid_dims, NULL);
  for (int d = 0; d < 3; d++) {
    file_start[d] = global_start[d];
    file_count[d] = mdims[d];
  }
  file_start[3] = 0;
  file_count[3] = NVAR;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, 
    NULL);

  for (int d = 0; d < 3; d++) {
    mem_grid_dims[d] = mdims[d] + 2*NG;
  }
  mem_grid_dims[3] = NVAR;
  memspace = H5Screate_simple(4, mem_grid_dims, NULL);
  for (int d = 0; d < 3; d++)
    mem_start[d] = 0;
  mem_start[3] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, 
    NULL);

  sprintf(name,"p");
  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace,
    H5P_DEFAULT, plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id,
    &P[NG][NG][NG][0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(memspace);
  H5Sclose(filespace);

  H5Fflush(file_id, H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
  */
}

void restart_read(char *fname, struct FluidState *S)
{
  /*hsize_t file_grid_dims[4],file_start[4],file_count[4];
  hsize_t mem_grid_dims[4],mem_start[4],one;
  char *name = "p";
  int fdims[3] = {N1TOT, N2TOT, N3TOT};
  int mdims[3] = {N1, N2, N3};

  if(mpi_io_proc()) {
    fprintf(stderr, "Restarting from %s\n\n", fname);
  }

  one = 1;

  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  // Read header
  file_grid_dims[0] = 1;
  hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
  file_start[0] = 0;
  file_count[0] = 1;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, 
    NULL);
  hid_t memspace  = H5Screate_simple(1, &one, NULL);

  get_dbl_value(&t, "t", file_id, filespace, memspace);
  get_int_value(&nstep, "nstep", file_id, filespace, memspace);
  get_int_value(&dump_id, "dump_id", file_id, filespace, memspace);
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
  for (int d = 0; d < 3; d++)
    file_grid_dims[d] = fdims[d];
  file_grid_dims[3] = NVAR;  // For vectors
  filespace = H5Screate_simple(4, file_grid_dims, NULL);
  for (int d = 0; d < 3; d++) {
    file_start[d] = global_start[d];
    file_count[d] = mdims[d];
  }
  file_start[3] = 0;
  file_count[3] = NVAR;
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, 
    NULL);

  for (int d = 0; d < 3; d++) {
    mem_grid_dims[d] = mdims[d] + 2*NG;
  }
  mem_grid_dims[3] = NVAR;
  memspace = H5Screate_simple(4, mem_grid_dims, NULL);
  for (int d = 0; d < 3; d++)
    mem_start[d] = NG;
  mem_start[3] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, 
    NULL);

  plist_id = H5Pcreate(H5P_DATASET_ACCESS);
  hid_t dset_id = H5Dopen(file_id, name, plist_id);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, 
    &P[0][0][0][0]);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  H5Sclose(memspace);
  H5Sclose(filespace);
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);

  MPI_Barrier(MPI_COMM_WORLD);*/
}

int restart_init(struct GridGeom *G, struct FluidState *S)
{
  char fname[2048], lastname[2048];
  sprintf(fname, "restart.last");
  strcpy(lastname, restartdir);
  strcat(lastname, fname);

  FILE *fp = fopen(lastname,"r");
  if (fp == NULL) {
    if (mpi_io_proc())
      fprintf(stdout, "No restart file\n\n");
    return 0;
  }

  if (mpi_io_proc())
    fprintf(stdout, "Loading restart file\n\n");
  zero_arrays();

  safe_fscanf(fp, "%s\n", fname);
  restart_read(fname, S);

  //ZSLOOP(-NG, N3-1+NG, -NG, N2-1+NG, -NG, N1-1+NG) {
  //  PLOOP {
  //    Ph[i][j][k][ip] = P[i][j][k][ip];
  //  }
  //}

  set_grid(G);

  set_bounds(G, S);

  fprintf(stderr, "RESTART INIT NOT SUPPORTED WITH NEW ARRAY INDEXING\n");
  exit(-1);

  return 1;
}

// Error-handling wrappers for standard C functions
void safe_system(const char *command)
{
  int systemReturn = system(command);
  if (systemReturn == -1) {
    fprintf(stderr, "system() call %s failed! Exiting!\n", command);
    exit(-1);
  }
}

void safe_fscanf(FILE *stream, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  int vfscanfReturn = vfscanf(stream, format, args);
  va_end(args);
  if (vfscanfReturn == -1) {
    fprintf(stderr, "fscanf() call failed! Exiting!\n");
    exit(-1);
  }
}

