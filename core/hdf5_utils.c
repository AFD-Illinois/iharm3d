/******************************************************************************
 *                                                                            *
 * HDF5_UTILS.C                                                               *
 *                                                                            *
 * FUNCTIONS FOR DEALING WITH HDF5                                            *
 *                                                                            *
 ******************************************************************************/

#include <hdf5.h>
#include "decs.h"

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

// Low level function without pre-filled
void hdf5_write_array(void *data, hid_t file_id, const char *name, hsize_t rank,
  hsize_t *fdims, hsize_t *fstart, hsize_t *mdims, hsize_t *mstart, hsize_t type)
{
  hid_t filespace = H5Screate_simple(rank, fdims, NULL);
  H5Sselect_hyperslab(filespace, H5S_SELECT_SET, fstart, NULL, mdims,
	NULL);
  hid_t memspace = H5Screate_simple(rank, mdims, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mstart, NULL, mdims,
	NULL);

  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, name, type, filespace, H5P_DEFAULT,
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

// Writing is independent of type, so I define the following
void hdf5_write_scalar(void *data, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2]};
  hsize_t mdims[] = {N1, N2, N3};
  hsize_t mstart[] = {0, 0, 0};

  hdf5_write_array(data, file_id, name, 3,
    fdims, fstart, mdims, mstart, hdf5_type);
}

void hdf5_write_vector(void *data, const char *name, hid_t file_id, int len, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT, len};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2], 0};
  hsize_t mdims[] = {N1, N2, N3, len};
  hsize_t mstart[] = {0, 0, 0, 0};

  hdf5_write_array(data, file_id, name, 4,
    fdims, fstart, mdims, mstart, hdf5_type);
}

void hdf5_write_tensor(void *data, const char *name, hid_t file_id, int n1, int n2, hsize_t hdf5_type)
{
  hsize_t fdims[] = {N1TOT, N2TOT, N3TOT, n1, n2};
  hsize_t fstart[] = {global_start[0], global_start[1], global_start[2], 0, 0};
  hsize_t mdims[] = {N1, N2, N3, n1, n2};
  hsize_t mstart[] = {0, 0, 0, 0, 0};

  hdf5_write_array(data, file_id, name, 5,
    fdims, fstart, mdims, mstart, hdf5_type);
}

void hdf5_write_header_val(void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{
  // TODO add to /header
  hid_t scalarspace = H5Screate(H5S_SCALAR);
  hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
  hid_t dset_id = H5Dcreate(file_id, name, hdf5_type, scalarspace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, hdf5_type, scalarspace, scalarspace, plist_id, val);

  H5Dclose(dset_id);
  H5Pclose(plist_id);
  H5Sclose(scalarspace);
}

void hdf5_write_attributes(void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
{

}

void hdf5_read_header_val(void *val, const char *name, hid_t file_id, hsize_t hdf5_type)
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
