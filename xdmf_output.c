
#include <hdf5.h>
#include "decs.h"

/*#################################
  #################################

    Generic XDMF output functions

  #################################
  #################################*/

void write_xml_closing(FILE *xml) {

    fprintf(xml, "   </Grid>\n");
    fprintf(xml, " </Domain>\n");
    fprintf(xml, "</Xdmf>\n");

    fclose(xml);
}

FILE *write_xml_head(int dump_id, double t) {

	char name[80];
	FILE *fp;

	sprintf(name,"dumps/dump_%05d.xmf", dump_id);
	fp = fopen(name,"w");
	if(fp == NULL) {
		fprintf(stderr,"Failed to open xmf file...perhaps there is no \"dumps\" directory???\n");
		exit(123);
	}

    fprintf(fp, "<?xml version=\"1.0\" ?>\n");
    fprintf(fp, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(fp, "<Xdmf Version=\"2.0\">\n");
    fprintf(fp, " <Domain>\n");
    fprintf(fp, "   <Grid Name=\"mesh\" GridType=\"Uniform\">\n");
	fprintf(fp, "     <Time Value=\"%16.14e\"/>\n", t);
    fprintf(fp, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", N1TOT+1, N2TOT+1, N3TOT+1);
	fprintf(fp, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1TOT+1, N2TOT+1, N3TOT+1);
    fprintf(fp, "        grid.h5:/X\n");
    fprintf(fp, "       </DataItem>\n");
    fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1TOT+1, N2TOT+1, N3TOT+1);
    fprintf(fp, "        grid.h5:/Y\n");
    fprintf(fp, "       </DataItem>\n");
    fprintf(fp, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1TOT+1, N2TOT+1, N3TOT+1);
    fprintf(fp, "        grid.h5:/Z\n");
    fprintf(fp, "       </DataItem>\n");

    fprintf(fp, "     </Geometry>\n");

	return fp;
}

void write_scalar(float *data, hid_t file_id, const char *name, hid_t filespace, hid_t memspace, FILE *xml) {

	char fname[80];
	hid_t plist_id, dset_id;

	H5Fget_name(file_id, fname, 80);
	char *sname = strrchr(fname,'/');
	if(sname != NULL) sname++;
	else sname = fname;

	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	dset_id = H5Dcreate(file_id, name, H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	if(mpi_io_proc()) {
		fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", name);
		fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", N1TOT, N2TOT, N3TOT);
		fprintf(xml, "        %s:%s\n", sname, name);
		fprintf(xml, "       </DataItem>\n");
		fprintf(xml, "     </Attribute>\n");
	}

	return;
}

void write_scalar_dbl(double *data, hid_t file_id, char *name, hid_t filespace, hid_t memspace, FILE *xml) {

	char fname[80];
	hid_t plist_id, dset_id;

	H5Fget_name(file_id, fname, 80);
	char *sname = strrchr(fname,'/');
	if(sname != NULL) sname++;
	else sname = fname;

	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	if(mpi_io_proc()) {
		fprintf(xml, "     <Attribute Name=\"%s\" AttributeType=\"Scalar\" Center=\"Cell\">\n", name);
		fprintf(xml, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n", N1TOT, N2TOT, N3TOT);
		fprintf(xml, "        %s:%s\n", sname, name);
		fprintf(xml, "       </DataItem>\n");
		fprintf(xml, "     </Attribute>\n");
	}

	return;
}

void dump_grid() {

	int i,j,k,d;
	float *x[3];
	double xp[4],r,th;
	hid_t dset_id;
	hsize_t file_start[3], file_count[3];
	hsize_t mem_start[3];
	hsize_t fdims[] = {N1TOT+1, N2TOT+1, N3TOT+1};
	hsize_t dims[] = {N1, N2, N3};
	const char *coordNames[] = {"/X", "/Y", "/Z"};

	for(d = 0; d < 3; d++) {
		if(global_stop[d] == fdims[d]-1) dims[d]++;
	}

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	hid_t file_id = H5Fcreate("dumps/grid.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	hid_t filespace = H5Screate_simple(3, fdims, NULL);
	for(d = 0; d < 3; d++) {
		file_start[d] = global_start[d];
		file_count[d] = dims[d];
	}
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

	hid_t memspace = H5Screate_simple(3, dims, NULL);
	for(d = 0; d < 3; d++) mem_start[d] = 0;
	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

	int total_size = 1;
	for(d = 0; d < 3; d++) total_size *= dims[d];
	for(d = 0; d < 3; d++) {
		x[d] = malloc(total_size*sizeof(float));
		if(x[d] == NULL) {
			fprintf(stderr,"Failed to allocate x[d] in dump_grid\n");
			exit(5);
		}
	}

	int ind = 0;
	ZSLOOP(0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1) {
			coord(i,j,CORN,xp);
			xp[3] = (k + global_start[2] - START3)*dx[3];
			bl_coord(xp, &r, &th);
			x[0][ind] = r*cos(xp[3])*sin(th);
			x[1][ind] = r*sin(xp[3])*sin(th);
			x[2][ind] = r*cos(th);
			ind++;
	}

	for(d = 0; d < 3; d++){
		plist_id = H5Pcreate(H5P_DATASET_CREATE);
		dset_id = H5Dcreate(file_id, coordNames[d], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		H5Pclose(plist_id);

		plist_id = H5Pcreate(H5P_DATASET_XFER);
	    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
		H5Dwrite(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, x[d]);
		H5Dclose(dset_id);
		H5Pclose(plist_id);
	}

	H5Sclose(filespace);
	H5Sclose(memspace);

	for(d = 0; d < 3; d++) free(x[d]);

	H5Fclose(file_id);

	return;
}
