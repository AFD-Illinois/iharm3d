/******************************************************************************
 *                                                                            *
 * HDF5_UTILS.C                                                               *
 *                                                                            *
 * FUNCTIONS FOR DEALING WITH HDF5                                            *
 *                                                                            *
 ******************************************************************************/

#include <hdf5.h>
#include "decs.h"

static char hdf5_cur_dir[STRLEN] = "/";

hid_t hdf5_open(char *fname)
{
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL); // TODO tune HDF with an MPI info object
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  return file_id;
}

void hdf5_close(hid_t file_id)
{
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
}

void hdf5_make_directory(const char *name, hid_t file_id)
{
  // Add current directory to group name
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strnlen(path, STRLEN));

  hid_t group_id = H5Gcreate2(file_id, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}

void hdf5_set_directory(const char *path)
{
  strncpy(hdf5_cur_dir, path, STRLEN);
}

// Set our string type size.
// H5T_VARIABLE indicates any string but isn't compatible w/parallel IO
hid_t hdf5_make_str_type(hsize_t len)
{
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, len);
  return string_type;
}

void hdf5_add_att(const void *att, const char *att_name, const char *data_name, hid_t file_id, hsize_t hdf5_type)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, data_name, STRLEN - strnlen(path, STRLEN));

  hid_t attribute_id = H5Acreate_by_name(file_id, path, att_name, hdf5_type, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, hdf5_type, att);
  H5Aclose(attribute_id);
}

void hdf5_add_units(const char *name, const char *unit, hid_t file_id)
{
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strnlen(unit,STRLEN)+1);
  hdf5_add_att(unit, "units", name, file_id, string_type);
  H5Tclose(string_type);
}

// Write a 1D list of strings to a file (for labeling primitives array)
// Must be a constant-length array, i.e. strs[len][strlen]
void hdf5_write_str_list(const void *data, const char *name, hid_t file_id, hsize_t strlen, hsize_t len)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strnlen(path, STRLEN));

  // Taken (stolen) from https://support.hdfgroup.org/ftp/HDF5/examples/C/
  hsize_t dims_of_char_array[] = {len};
  hsize_t dims_of_char_dataspace[] = {1};

  hid_t vlstr_h5t = H5Tcopy(H5T_C_S1);
  H5Tset_size(vlstr_h5t, strlen);

  hid_t mem_h5t = H5Tarray_create(vlstr_h5t, 1, dims_of_char_array);
  hid_t dataspace = H5Screate_simple(1, dims_of_char_dataspace, NULL); // use same dims as int ds
  hid_t dataset = H5Dcreate(file_id, path, mem_h5t, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, mem_h5t, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Tclose(mem_h5t);
  H5Tclose(vlstr_h5t);
}

// Low level function for any array/tensor size
// Use convenience functions below instead
void hdf5_write_array(const void *data, hid_t file_id, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *mdims, hsize_t *mstart, hsize_t type)
{
  hid_t filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, mdims,
	NULL);
  hid_t memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, mdims,
	NULL);

  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strnlen(path, STRLEN));

  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, path, type, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, type, memspace, filespace, plist_id, data);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);

}

// Write a packed buffer of scalars to a file
void hdf5_write_scalar(const void *data, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2]};
  hsize_t mdims[] = {N1, N2, N3};
  hsize_t mstart[] = {0, 0, 0};

  hdf5_write_array(data, file_id, name, 3,
    fdims, fstart, mdims, mstart, hdf5_type);
}

void hdf5_write_vector(const void *data, const char *name, hid_t file_id, hsize_t len, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT, len};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2], 0};
  hsize_t mdims[] = {N1, N2, N3, len};
  hsize_t mstart[] = {0, 0, 0, 0};

  hdf5_write_array(data, file_id, name, 4,
    fdims, fstart, mdims, mstart, hdf5_type);
}

void hdf5_write_tensor(const void *data, const char *name, hid_t file_id, hsize_t n1, hsize_t n2, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT, n1, n2};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2], 0, 0};
  hsize_t mdims[] = {N1, N2, N3, n1, n2};
  hsize_t mstart[] = {0, 0, 0, 0, 0};

  hdf5_write_array(data, file_id, name, 5,
    fdims, fstart, mdims, mstart, hdf5_type);
}

void hdf5_write_single_val(const void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strnlen(path, STRLEN));

  hid_t scalarspace = H5Screate(H5S_SCALAR);
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, path, hdf5_type, scalarspace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, hdf5_type, scalarspace, scalarspace, plist_id, val);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(scalarspace);
}

void hdf5_read_single_val(const void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{

}

// Following double/float defs differ only in signature. C is a bear.
void pack_scalar_double(double in[N3+2*NG][N2+2*NG][N1+2*NG], double *out)
{
  int ind = 0;
  ZLOOP_OUT {
    out[ind] = in[k][j][i];
    ind++;
  }
}

void pack_scalar_float(double in[N3+2*NG][N2+2*NG][N1+2*NG], float *out)
{
  int ind = 0;
  ZLOOP_OUT {
    out[ind] = in[k][j][i];
    ind++;
  }
}

void pack_scalar_int(int in[N3+2*NG][N2+2*NG][N1+2*NG], int *out)
{
  int ind = 0;
  ZLOOP_OUT {
    out[ind] = in[k][j][i];
    ind++;
  }
}

void pack_vector_double(double in[][N3+2*NG][N2+2*NG][N1+2*NG], double *out, int vector_len)
{
  int ind = 0;
  for (int mu=0; mu < vector_len; mu++) {
    ZLOOP_OUT {
      out[ind] = in[mu][k][j][i];
      ind++;
    }
  }
}

void pack_vector_float(double in[][N3+2*NG][N2+2*NG][N1+2*NG], float *out, int vector_len)
{
  int ind = 0;
  ZLOOP_OUT {
    for (int mu=0; mu < vector_len; mu++) {
      out[ind] = in[mu][k][j][i];
      ind++;
    }
  }
}

// TODO more functions if I want to output F or other tensors
// Also: this wastes quite a bit of space writing all phi
void pack_Gtensor_double(double in[NDIM][NDIM][N2+2*NG][N1+2*NG], double *out)
{
  int ind = 0;
  ZLOOP_OUT {
    DLOOP2 {
	  out[ind] = in[nu][mu][j][i];
	  ind++;
    }
  }
}

void pack_Gtensor_float(double in[NDIM][NDIM][N2+2*NG][N1+2*NG], float *out)
{
  int ind = 0;
  ZLOOP_OUT {
    DLOOP2 {
	  out[ind] = in[nu][mu][j][i];
	  ind++;
    }
  }
}
