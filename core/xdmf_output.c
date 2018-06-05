/******************************************************************************
 *                                                                            *
 * XDMF_OUTPUT.C                                                              *
 *                                                                            *
 * XDMF OUTPUT TO ALLOW FOR SEPARATE GRID OUTPUT FILE                         *
 *                                                                            *
 ******************************************************************************/

#include <hdf5.h>
#include "decs.h"

void write_xml_closing(FILE *xml)
{
  fprintf(xml, "   </Grid>\n");
  fprintf(xml, " </Domain>\n");
  fprintf(xml, "</Xdmf>\n");

  fflush(xml);

  fclose(xml);
}

FILE *write_xml_head(char *fname, double t)
{
  FILE *fp;

  fp = fopen(fname,"w");
  if(fp == NULL) {
    fprintf(stderr, "Failed to open xmf file %s\n", fname);
    exit(123);
  }

  fprintf(fp, "<?xml version=\"1.0\" ?>\n");
  fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
  fprintf(fp, "<Xdmf Version=\"2.0\">\n");
  fprintf(fp, " <Domain>\n");
  fprintf(fp, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
  fprintf(fp, "     <Time Value=\"%16.14e\"/>\n", t);
  fprintf(fp, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n",
    N1TOT+1, N2TOT+1, N3TOT+1);
  fprintf(fp, "     <Geometry GeometryType=\"X_Y_Z\">\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
    N1TOT+1, N2TOT+1, N3TOT+1);
  fprintf(fp, "        grid.h5:/X\n");
  fprintf(fp, "       </DataItem>\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
    N1TOT+1, N2TOT+1, N3TOT+1);
  fprintf(fp, "        grid.h5:/Y\n");
  fprintf(fp, "       </DataItem>\n");
  fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1TOT+1, N2TOT+1, N3TOT+1);
  fprintf(fp, "        grid.h5:/Z\n");
  fprintf(fp, "       </DataItem>\n");

  fprintf(fp, "     </Geometry>\n");

  return fp;
}

void write_scalar(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml)
{
  char fname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if(sname != NULL) sname++;
  else sname = fname;

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(mpi_io_proc()) {
    fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
      name);
    fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT);
    fprintf(xml, "        %s:%s\n", sname, name);
    fprintf(xml, "       </DataItem>\n");
    fprintf(xml, "     </Attribute>\n");
  }
}

void write_vector(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml)
{
  char fname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if(sname != NULL) sname++;
  else sname = fname;

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(mpi_io_proc()) {
    fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
      name);
    fprintf(xml, "       <DataItem Dimensions=\"%d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT, NDIM);
    fprintf(xml, "        %s:%s\n", sname, name);
    fprintf(xml, "       </DataItem>\n");
    fprintf(xml, "     </Attribute>\n");
  }
}

void write_tensor(float *data, hid_t file_id, const char *name, hid_t filespace,
  hid_t memspace, FILE *xml)
{
  char fname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if(sname != NULL) sname++;
  else sname = fname;

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(mpi_io_proc()) {
    fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
      name);
    fprintf(xml, "       <DataItem Dimensions=\"%d %d %d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT, NDIM, NDIM);
    fprintf(xml, "        %s:%s\n", sname, name);
    fprintf(xml, "       </DataItem>\n");
    fprintf(xml, "     </Attribute>\n");
  }
}

void write_scalar_dbl(double *data, hid_t file_id, char *name, hid_t filespace,
  hid_t memspace, FILE *xml)
{
  char fname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if (sname != NULL) {
    sname++;
  }
  else sname = fname;

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(mpi_io_proc()) {
    fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
      name);
    fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT);
    fprintf(xml, "        %s:%s\n", sname, name);
    fprintf(xml, "       </DataItem>\n");
    fprintf(xml, "     </Attribute>\n");
  }
}

void write_scalar_int(int *data, hid_t file_id, const char *name,
  hid_t filespace, hid_t memspace, FILE *xml)
{
  char fname[80];
  hid_t plist_id, dset_id;

  H5Fget_name(file_id, fname, 80);
  char *sname = strrchr(fname,'/');
  if(sname != NULL) sname++;
  else sname = fname;

  plist_id = H5Pcreate(H5P_DATASET_CREATE);
  dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT,
    plist_id, H5P_DEFAULT);
  H5Pclose(plist_id);

  plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, data);
  H5Dclose(dset_id);
  H5Pclose(plist_id);

  if(mpi_io_proc()) {
    fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n",
      name);
    fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Int\" Precision=\"4\" Format=\"HDF\">\n",
      N1TOT, N2TOT, N3TOT);
    fprintf(xml, "        %s:%s\n", sname, name);
    fprintf(xml, "       </DataItem>\n");
    fprintf(xml, "     </Attribute>\n");
  }
}
