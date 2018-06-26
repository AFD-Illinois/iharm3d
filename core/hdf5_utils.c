/******************************************************************************
 *                                                                            *
 * HDF5_UTILS.C                                                               *
 *                                                                            *
 * FUNCTIONS FOR DEALING WITH HDF5                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#include <hdf5.h>

// The library remembers a "current directory" for convenience
// This is stateful, so usual caveats about resetting it
static char hdf5_cur_dir[STRLEN] = "/";

void hdf5_set_directory(const char *path)
{
  strncpy(hdf5_cur_dir, path, STRLEN);
}

// Create a new HDF file (or overwrite whatever file exists)
hid_t hdf5_create(char *fname)
{
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL); // TODO tune HDF with an MPI info object
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);
  return file_id;
}

// Open an existing file for reading
hid_t hdf5_open(char *fname)
{
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);
  return file_id;
}

// Close a file
void hdf5_close(hid_t file_id)
{
  H5Fflush(file_id,H5F_SCOPE_GLOBAL);
  H5Fclose(file_id);
}

// Make a directory (in the current directory) with given name
// This doesn't take a full path, just a name
void hdf5_make_directory(const char *name, hid_t file_id)
{
  // Add current directory to group name
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  hid_t group_id = H5Gcreate2(file_id, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}

// Return a fixed-size string type
// H5T_VARIABLE indicates any string but isn't compatible w/parallel IO
hid_t hdf5_make_str_type(size_t len)
{
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, len);
  return string_type;
}

// Add the named attribute to the named variable
void hdf5_add_att(const void *att, const char *att_name, const char *data_name, hid_t file_id, hsize_t hdf5_type)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, data_name, STRLEN - strlen(path));

  hid_t attribute_id = H5Acreate_by_name(file_id, path, att_name, hdf5_type, H5Screate(H5S_SCALAR), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attribute_id, hdf5_type, att);
  H5Aclose(attribute_id);
}

// Add an attribute named "units"
void hdf5_add_units(const char *name, const char *unit, hid_t file_id)
{
  hid_t string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, strlen(unit)+1);
  hdf5_add_att(unit, "units", name, file_id, string_type);
  H5Tclose(string_type);
}

// Write a 1D list of strings (used for labeling primitives array)
// Must be an array of constant-length strings, i.e. char strs[len][str_len] = etc.
void hdf5_write_str_list(const void *data, const char *name, hid_t file_id, size_t str_len, size_t len)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  // Taken (stolen) from https://support.hdfgroup.org/ftp/HDF5/examples/C/
  hsize_t dims_of_char_array[] = {len};
  hsize_t dims_of_char_dataspace[] = {1};

  hid_t vlstr_h5t = H5Tcopy(H5T_C_S1);
  H5Tset_size(vlstr_h5t, str_len);

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
void hdf5_write_array(const void *data, hid_t file_id, const char *name, size_t rank,
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
  strncat(path, name, STRLEN - strlen(path));

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

// Write a C-order, N{1,2,3}-order array of hdf5_type to a file
void hdf5_write_scalar(const void *data, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2]};
  hsize_t mdims[] = {N1, N2, N3};
  hsize_t mstart[] = {0, 0, 0};

  hdf5_write_array(data, file_id, name, 3,
    fdims, fstart, mdims, mstart, hdf5_type);
}

// Write a C-order, N{1,2,3},len-order array of hdf5_type to a file
void hdf5_write_vector(const void *data, const char *name, hid_t file_id, size_t len, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT, len};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2], 0};
  hsize_t mdims[] = {N1, N2, N3, len};
  hsize_t mstart[] = {0, 0, 0, 0};

  hdf5_write_array(data, file_id, name, 4,
    fdims, fstart, mdims, mstart, hdf5_type);
}

// Write a C-order, N{1,2,3},n1,n2-order array of hdf5_type to a file
void hdf5_write_tensor(const void *data, const char *name, hid_t file_id, size_t n1, size_t n2, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT, n1, n2};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2], 0, 0};
  hsize_t mdims[] = {N1, N2, N3, n1, n2};
  hsize_t mstart[] = {0, 0, 0, 0, 0};

  hdf5_write_array(data, file_id, name, 5,
    fdims, fstart, mdims, mstart, hdf5_type);
}

// Write a single value of hdf5_type to a file
void hdf5_write_single_val(const void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

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

// HARM-specific code for restarts
// These are separate because when not packed, we have to tell HDF about ghost zones
// TODO pack/unpack prims for restarts instead? HDF's a bear here

// Write a C-order, len,N{3,2,1}-order array of doubles to a file
// Note this also skips ghost zones automatically
void hdf5_write_restart_prims(const void *data, const char *name, hid_t file_id)
{
  hsize_t fdims[] = {NVAR, N3TOT, N2TOT, N1TOT};
  hsize_t fstart[] = {0, global_start[2], global_start[1], global_start[0]};
  hsize_t mdims[] = {NVAR, N3, N2, N1};
  hsize_t mfulldims[] = {NVAR, N3+2*NG, N2+2*NG, N1+2*NG};
  hsize_t mstart[] = {0, NG, NG, NG};

  hid_t filespace = H5Screate_simple(4, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, mdims,
    NULL);
  hid_t memspace = H5Screate_simple(4, mfulldims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, mdims,
    NULL);

  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, path, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
}

// These are very like above but there's not a good way to share code...
void hdf5_read_single_val(void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  hid_t scalarspace = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(file_id, path, H5P_DEFAULT);

  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, hdf5_type, scalarspace, scalarspace, plist_id, val);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(scalarspace);

}

void hdf5_read_restart_prims(void *data, const char *name, hid_t file_id)
{
  hsize_t fdims[] = {NVAR, N3TOT, N2TOT, N1TOT};
  hsize_t fstart[] = {0, global_start[2], global_start[1], global_start[0]};
  hsize_t mdims[] = {NVAR, N3, N2, N1};
  hsize_t mfulldims[] = {NVAR, N3+2*NG, N2+2*NG, N1+2*NG};
  hsize_t mstart[] = {0, NG, NG, NG};

  hid_t filespace = H5Screate_simple(4, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, mdims,
    NULL);
  hid_t memspace = H5Screate_simple(4, mfulldims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, mdims,
    NULL);

  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  hid_t dset_id = H5Dopen(file_id, path, H5P_DEFAULT);

  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
}
