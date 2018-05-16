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

void write_scalar(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml);
void write_scalar_int(int *data, hid_t file_id, const char *name,
  hid_t filespace, hid_t memspace, FILE *xml);

// Defined end of this file
void add_int_value(int val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace);
void add_dbl_value(double val, const char *name, hid_t file_id, hid_t filespace,
  hid_t memspace);
void add_str_value(const char* val, const char *name, hid_t file_id,
  hid_t filespace, hid_t memspace);

void get_dbl_value(double *val, const char *name, hid_t file_id,
                   hid_t filespace, hid_t memspace);
void get_int_value(int *val, const char *name, hid_t file_id, hid_t filespace,
                   hid_t memspace);

static int dump_id = 0;
static int restart_id = 0;

// O(n) dictionary for converting strings to pointers to global variables
// TODO I think this relies on leaking memory.  Not so great
struct paramtable_entry {
  char *key;
  void *data;
};

#define MAXPARAMS (1000)
static struct paramtable_entry paramtable[MAXPARAMS];
static int nparam = 0;
static int nparamset = 0;

void set_param(char *key, void *data) {
  paramtable[nparam].key = key;
  paramtable[nparam].data = data;
  nparam++;
}

int get_param(char *key, void **data) {
  int n = 0;
  while (strncmp(key, paramtable[n].key, strlen(key)) != 0) {
    n++;
    if (n >= nparam) {
      *data = NULL;
      return 0;
    }
  }
  *data = paramtable[n].data;

  return 1;
}

void set_core_params() {
  set_param("tf", &tf);
  set_param("dt", &dt);
#if METRIC == MINKOWSKI
  set_param("x1Min", &x1Min);
  set_param("x1Max", &x1Max);
  set_param("x2Min", &x2Min);
  set_param("x2Max", &x2Max);
  set_param("x3Min", &x3Min);
  set_param("x3Max", &x3Max);
#elif METRIC == MKS
  set_param("a", &a);
  set_param("hslope", &hslope);
#if POLYTH
  set_param("poly_xt", &poly_xt);
  set_param("poly_alpha", &poly_alpha);
  set_param("mks_smooth", &mks_smooth);
#endif
  if (N2 < NG) hslope = 1.;
  set_param("Rout", &Rout);
  #if RADIATION
  set_param("Rout_rad", &Rout_rad);
  #endif
#endif

  set_param("cour", &cour);
  set_param("gam", &gam);

  #if RADIATION
  #if METRIC == MINKOWSKI
  set_param("L_unit", &L_unit);
  set_param("M_unit", &M_unit);
  #elif METRIC == MKS
  set_param("mbh", &mbh);
  set_param("M_unit", &M_unit);
  #endif
  #endif

  #if ELECTRONS
  set_param("game", &game);
  set_param("gamp", &gamp);
  set_param("fel0", &fel0);
  set_param("tptemin", &tptemin);
  set_param("tptemax", &tptemax);
  #endif

  #if RADIATION
  #if !ELECTRONS
  set_param("tp_over_te", &tp_over_te);
  #endif
  set_param("nph_per_proc", &nph_per_proc);
  set_param("numin", &numin);
  set_param("numax", &numax);
  set_param("tune_emiss", &tune_emiss);
  set_param("tune_scatt", &tune_scatt);
  set_param("t0_tune_emiss", &t_tune_emiss);
  set_param("t0_tune_scatt", &t_tune_scatt);
  set_param("thetae_max", &thetae_max);
  set_param("sigma_max", &sigma_max);
  set_param("kdotk_tol", &kdotk_tol);
  set_param("Nph_to_track", &Nph_to_track);
  #endif

  set_param("DTd", &DTd);
  set_param("DTl", &DTl);
  set_param("DTr", &DTr);
  set_param("DTp", &DTp);
  // I set output dir via command line, since it's machine-specific
}

// Set runtime parameters from file
void read_params(char* pfname)
{
  void *ptr;

  FILE *fp = fopen(pfname, "r");
  if (fp == NULL) {
    fprintf(stderr, "Cannot open parameter file: %s\n", pfname);
    exit(-1);
  }

  char line[STRLEN];
  while (fgets(line, STRLEN, fp)) {
    // Ignore comments, newlines, and leading whitespace
    if (line[0] == '#' || line[0] == '\n' || isspace(line[0]))
      continue;

    // Is key in dictionary, and is variable empty?
    char test[STRLEN], key[STRLEN];
    test[0] = '\0';
    sscanf(line, "%*s %s %*s %s", key, test);
    char *word = test;
    while(isspace(*word)) {
      word++;
    }
    if (word[0] == '\0') {
      continue;
    }

    // Read in parameter depending on datatype
    char type[6];
    strncpy(type, line, 5);
    type[5] = 0;
    if (get_param(key, &ptr)) {
      if (!strncmp(type, "[int]", 5)) {
        int buf;
        sscanf(line, "%*s %s %*s %d", key, &buf);
        *((int*)ptr) = buf;
        nparamset++;
      } else if (!strncmp(type, "[dbl]", 5)) {
        double buf;
        sscanf(line, "%*s %s %*s %lf", key, &buf);
        *((double*)ptr) = buf;
        nparamset++;
      } else if (!strncmp(type, "[str]", 5)) {
        char buf[STRLEN];
        sscanf(line, "%*s %s %*s %s", key, buf);
        strcpy((char*)ptr, buf);
        nparamset++;
      }
    }
  }

  #if (METRIC == MKS) && RADIATION
  Mbh = mbh*MSUN;
  #endif

  if (nparamset != nparam && mpi_io_proc()) {
    fprintf(stderr, "Set %i parameters, needed %i!\n", nparamset, nparam);
    exit(-1);
  }

  fclose(fp);

  if (mpi_io_proc()) fprintf(stdout, "Parameter file read\n\n");
}

void init_io(char *pfname)
{
  set_core_params();
  set_problem_params();
  read_params(pfname);
}

void dump(struct GridGeom *G, struct FluidState *S)
{
  static int firstc = 1;
  static float *data;
  static int *idata;
  char fname[80];
  FILE *xml = NULL;
  hsize_t fdims[] = {N3TOT, N2TOT, N1TOT};
  hsize_t mem_copy_dims[] = {N3, N2, N1};
  #if ELECTRONS
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3",
                            "KEL", "KTOT"};
  #else
  const char *varNames[] = {"RHO", "UU", "U1", "U2", "U3", "B1", "B2", "B3"};
  #endif

  hid_t file_id, filespace, memspace, plist_id;//, dset_id;
  hsize_t mem_start[3], file_start[3], one, zero;

  one = 1;
  zero = 0;

  if(firstc) {
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
  filespace = H5Screate_simple(1, &one, NULL);
  hsize_t single_var_count = (mpi_io_proc() ? one : zero);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &zero, NULL, &single_var_count,
    NULL);
  memspace  = H5Screate_simple(1, &single_var_count, NULL);
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
    file_start[d] = global_start[2-d];
  }
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, mem_copy_dims,
    NULL);

  memspace = H5Screate_simple(3, mem_copy_dims, NULL);
  for (int d = 0; d < 3; d++)
    mem_start[d] = 0;
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, mem_copy_dims,
    NULL);

  /* printf("Size of memspace: %lld x %lld x %lld\n", */
  /*       mem_copy_dims[0], mem_copy_dims[1], mem_copy_dims[2]); */
  /* printf("Location of memspace: %lld %lld %lld\n", */
  /*        file_start[0], file_start[1], file_start[2]); */
  /* printf("Size of filespace: %lld x %lld x %lld\n", */
  /*       fdims[0], fdims[1], fdims[2]); */

  // Write primitive variables
  PLOOP {
    int ind = 0;
    ZLOOP {
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
    data[ind] = bsq;
    ind++;
  }
  write_scalar(data, file_id, "bsq", filespace, memspace, xml);

  ind = 0;
  ZLOOP {
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

void restart_write(struct FluidState *S)
{
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

  MPI_Barrier(MPI_COMM_WORLD);
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

#if METRIC == MINKOWSKI
  // TODO these are not written to restart files, but /tend/ to be 0,1
  x1Min = 0.;
  x1Max = 1.;
  x2Min = 0.;
  x2Max = 1.;
  x3Min = 0.;
  x3Max = 1.;
#endif

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
