/******************************************************************************
 *                                                                            *
 * HDF5_UTILS.C                                                               *
 *                                                                            *
 * FUNCTIONS FOR DEALING WITH HDF5                                            *
 *                                                                            *
 ******************************************************************************/

#include "hdf5_utils.h"

#include <string.h>
#include <hdf5.h>

// The library remembers a "current directory" for convenience
// This is stateful, so usual caveats about resetting it
#define STRLEN 2048
static char hdf5_cur_dir[STRLEN] = "/";

#define DEBUG 0

// Create a new HDF file (or overwrite whatever file exists)
hid_t hdf5_create(char *fname)
{
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL); // TODO tune HDF with an MPI info object
  hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  // Everyone expects directory to be root after open
  hdf5_set_directory("/");

  return file_id;
}

// Open an existing file for reading
hid_t hdf5_open(char *fname)
{
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
  H5Pclose(plist_id);

  // Everyone expects directory to be root after open
  hdf5_set_directory("/");

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

  if(DEBUG) printf("Adding dir %s\n", path);

  hid_t group_id = H5Gcreate2(file_id, path, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(group_id);
}

// Set the current directory
void hdf5_set_directory(const char *path)
{
  strncpy(hdf5_cur_dir, path, STRLEN);
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

  if(DEBUG) printf("Adding att %s\n", path);

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

  if(DEBUG) printf("Adding str list %s\n", path);

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

// Write the section 'mdims_copy' starting at 'mstart' of a C-order array of rank 'rank' and size 'mdims_full'
// To the section 'fdims' at 'fstart' of the file 'file_id'
void hdf5_write_array(const void *data, hid_t file_id, const char *name, size_t rank,
                      hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type)
{
  // Declare spaces of the right size
  hid_t filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount,
    NULL);
  hid_t memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount,
    NULL);

  // Add our current path to the dataset name
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  if(DEBUG) printf("Adding array %s\n", path);

  // Create the dataset in the file
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, path, hdf5_type, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  // Conduct the transfer
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, hdf5_type, memspace, filespace, plist_id, data);

  // Close spaces (TODO could keep open for speed writing arrays of same size?)
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
}

// Write a single value of hdf5_type to a file
void hdf5_write_single_val(const void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  // Add current path to the dataset name
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  if(DEBUG) printf("Adding val %s\n", path);

  // Declare scalar spaces
  hid_t scalarspace = H5Screate(H5S_SCALAR);
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, path, hdf5_type, scalarspace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  // Conduct transfer
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, hdf5_type, scalarspace, scalarspace, plist_id, val);

  // Close spaces (TODO could definitely keep these open instead of re-declaring)
  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(scalarspace);
}

// These are very like above but there's not a good way to share code...
void hdf5_read_single_val(void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  if(DEBUG) printf("Reading val %s\n", path);

  hid_t scalarspace = H5Screate(H5S_SCALAR);
  hid_t dset_id = H5Dopen(file_id, path, H5P_DEFAULT);

  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, hdf5_type, scalarspace, scalarspace, plist_id, val);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(scalarspace);

}

void hdf5_read_array(void *data, hid_t file_id, const char *name, size_t rank,
                      hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type)
{
  hid_t filespace = H5Screate_simple(4, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, fcount,
    NULL);
  hid_t memspace = H5Screate_simple(4, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, fcount,
    NULL);

  char path[STRLEN];
  strncpy(path, hdf5_cur_dir, STRLEN);
  strncat(path, name, STRLEN - strlen(path));

  if(DEBUG) printf("Reading arr %s\n", path);

  hid_t dset_id = H5Dopen(file_id, path, H5P_DEFAULT);

  hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dread(dset_id, hdf5_type, memspace, filespace, plist_id, data);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(filespace);
  H5Sclose(memspace);
}
