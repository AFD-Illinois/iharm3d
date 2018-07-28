/*
 * hdf5_utils.h
 */

#pragma once

#include <hdf5.h>

// This lib uses a global debug flag if one exists
#ifndef DEBUG
#define DEBUG 0
#endif

// Use global MPI switch
#ifndef USE_MPI
// Or set it intelligently here
#ifdef MPI_COMM_WORLD
#define USE_MPI 1
#else
#define USE_MPI 0
#endif

#endif

// Blob "copy" utility
typedef hid_t hdf5_blob;
hdf5_blob hdf5_get_blob(const char *name);
int hdf5_write_blob(hdf5_blob blob, const char *name);
int hdf5_close_blob(hdf5_blob blob);

// File
int hdf5_create(const char *fname);
int hdf5_open(const char *fname);
int hdf5_close();

// Directory
int hdf5_make_directory(const char *name);
void hdf5_set_directory(const char *path);

// Write
int hdf5_write_single_val(const void *val, const char *name, hsize_t hdf5_type);
int hdf5_write_array(const void *data, const char *name, size_t rank,
                      hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type);

// Read
int hdf5_read_single_val(void *val, const char *name, hsize_t hdf5_type);
int hdf5_read_array(void *data, const char *name, size_t rank,
                      hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type);

// Convenience and annotations
hid_t hdf5_make_str_type(size_t len);
int hdf5_write_str_list(const void *data, const char *name, size_t strlen, size_t len);
int hdf5_add_attr(const void *att, const char *att_name, const char *data_name, hsize_t hdf5_type);
int hdf5_add_units(const char *name, const char *unit);

