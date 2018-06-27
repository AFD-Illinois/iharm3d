/*
 * hdf5_utils.h
 */

#pragma once

#include <hdf5.h>

// File
hid_t hdf5_create(char *fname);
hid_t hdf5_open(char *fname);
void hdf5_close(hid_t file_id);

// Directory
void hdf5_make_directory(const char *name, hid_t file_id);
void hdf5_set_directory(const char *path);

// Write
void hdf5_write_single_val(const void *val, const char *name, hid_t file_id, hsize_t hdf5_type);
void hdf5_write_array(const void *data, hid_t file_id, const char *name, size_t rank,
                      hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type);

// Read
void hdf5_read_single_val(void *val, const char *name, hid_t file_id, hsize_t hdf5_type);
void hdf5_read_array(void *data, hid_t file_id, const char *name, size_t rank,
                      hsize_t *fdims, hsize_t *fstart, hsize_t *fcount, hsize_t *mdims, hsize_t *mstart, hsize_t hdf5_type);

// Convenience and annotations
hid_t hdf5_make_str_type(size_t len);
void hdf5_write_str_list(const void *data, const char *name, hid_t file_id, size_t strlen, size_t len);
void hdf5_add_att(const void *att, const char *att_name, const char *data_name, hid_t file_id, hsize_t hdf5_type);
void hdf5_add_units(const char *name, const char *unit, hid_t file_id);
