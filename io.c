
#include "decs.h"
#include <hdf5.h>
#include <sys/stat.h>

void write_scalar(float *data, hid_t file_id, const char *name, hid_t filespace, hid_t memspace, FILE *xml);

static int dump_id = 0;
static int restart_id = 0;


void dump() 
{
	static int firstc = 1;
	static float *data;
	int i,j,k,d;
	char name[80];
	FILE *xml = NULL;
	hsize_t fdims[] = {N1TOT, N2TOT, N3TOT};
	hsize_t mdims[] = {N1, N2, N3};
	const char *varNames[] = {"rho", "u", "U1", "U2", "U3", "B1", "B2", "B3"};
	struct of_state q;
	
	hid_t file_id,filespace,memspace,plist_id,dset_id;
	hsize_t mem_start[3], file_start[3], one, zero;
	hsize_t file_count[3], file_grid_dims[3];

	one = 1;
	zero = 0;

	if(firstc) {
		mkdir("dumps", 0777);
		dump_grid();
		data = malloc(N1*N2*N3*NPR*sizeof(float));
		if(data == NULL) {
			fprintf(stderr,"failed to allocate data in dump\n");
			exit(8);
		}
		firstc = 0;
	}

	if(mpi_io_proc()) {
		xml = write_xml_head(dump_id,t);
	}

	sprintf(name,"dumps/dump_%05d.h5", dump_id);

	plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	/****** This is crazy.  All this to write out the time! *******/

	file_grid_dims[0] = 1;
	filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = (mpi_io_proc() ? one : zero);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	dset_id = H5Dcreate(file_id, "Time", H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &t);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
	H5Sclose(filespace);
	H5Sclose(memspace);

	/****** End dumping of time...seriously.  This is totally nuts *******/


	/****** Tell the HDF library about how the data is laid out in memory and in the file ******/
	for(d = 0; d < 3; d++) file_grid_dims[d] = fdims[d];
	filespace = H5Screate_simple(3, file_grid_dims, NULL);
	for(d = 0; d < 3; d++) {
		file_start[d] = global_start[d];
		file_count[d] = mdims[d];
	}
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

	memspace = H5Screate_simple(3, file_count, NULL);
	for(d = 0; d < 3; d++) mem_start[d] = 0;
	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

	/****** Done with memory/file layout ********/



	// write out the primitive variables
	PLOOP {
		int ind = 0;
		ZLOOP {
			data[ind] = (float)p[i][j][k][ip]; ind++;
		}
		write_scalar(data, file_id, varNames[ip], filespace, memspace, xml);
	}

	// example of writing something else out
	int ind = 0;
	ZLOOP {
		struct of_geom *geom = get_geometry(i, j, CENT) ;
		get_state(p[i][j][k], geom, &q);
		double bsq = 0.;
		for(int d = 0; d < NDIM; d++) bsq += q.bcon[d]*q.bcov[d];
		data[ind] = bsq;
		ind++;
	}
	write_scalar(data, file_id, "bsq", filespace, memspace, xml);

	ind = 0;
	ZLOOP {
		struct of_geom *geom = get_geometry(i, j, CENT) ;
		double gamma = 1.;
		mhd_gamma_calc(p[i][j][k], geom, &gamma);
		data[ind] = gamma;
		ind++;
	}
	write_scalar(data, file_id, "gamma", filespace, memspace, xml);
	

	// always make sure the following is at the end
	if(mpi_io_proc()) {
		write_xml_closing(xml);
	}

	H5Sclose(memspace);
	H5Sclose(filespace);
	H5Fclose(file_id);

	dump_id++;
}

void add_int_value(int val, const char *name, hid_t file_id, hid_t filespace, hid_t memspace)
{
	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
	hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_INT, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Dwrite(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &val);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
}


void add_dbl_value(double val, const char *name, hid_t file_id, hid_t filespace, hid_t memspace)
{
	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
	hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &val);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
}

void get_int_value(int *val, const char *name, hid_t file_id, hid_t filespace, hid_t memspace)
{
	hid_t plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	hid_t dset_id = H5Dopen(file_id, name, plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dread(dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, &val);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
}


void get_dbl_value(double *val, const char *name, hid_t file_id, hid_t filespace, hid_t memspace)
{
	hid_t plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	hid_t dset_id = H5Dopen(file_id, name, plist_id);
	H5Pclose(plist_id);
	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &val);
	H5Dclose(dset_id);
	H5Pclose(plist_id);
}


void restart_write() 
{
	int d;
	int fdims[3] = {N1TOT, N2TOT, N3TOT};
	int mdims[3] = {N1, N2, N3};
	hsize_t file_grid_dims[4],file_start[4],file_count[4];
	hsize_t mem_grid_dims[4],mem_start[4],one,zero;
	
	char name[80];

	one = 1;
	zero = 0;
	mkdir("restarts", 0777);

	sprintf(name,"restarts/restart_%05d.h5", restart_id);
	restart_id++;
	FILE *fp = fopen("restarts/restart.last","w");
	fprintf(fp,"%s\n", name);
	fclose(fp);
	if(mpi_io_proc()) {
		fprintf(stderr,"writing restart file %s\n", name);
	}
	
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	hid_t file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	// first write some info about the current state of the run (e.g. time, istep, ...)
	file_grid_dims[0] = 1;
	hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = (mpi_io_proc() ? one : zero);
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	hid_t memspace  = H5Screate_simple(1, (mpi_io_proc() ? &one : &zero), NULL);

	add_dbl_value(t, "Time", file_id, filespace, memspace);
	add_int_value(nstep, "Step", file_id, filespace, memspace);
	add_int_value(dump_id, "Next_Dump_ID", file_id, filespace, memspace);
	add_dbl_value(tf, "Final_Time", file_id, filespace, memspace);
	add_dbl_value(a, "a", file_id, filespace, memspace);
	add_dbl_value(gam, "gam", file_id, filespace, memspace);
	add_dbl_value(cour, "cour", file_id, filespace, memspace);
	add_dbl_value(DTd, "DTd", file_id, filespace, memspace);
	add_dbl_value(DTl, "DTl", file_id, filespace, memspace);
	add_dbl_value(DTi, "DTi", file_id, filespace, memspace);
	add_dbl_value(DTp, "DTp", file_id, filespace, memspace);		
	add_int_value(DTr, "DTr", file_id, filespace, memspace);
	add_int_value(restart_id, "restart_id", file_id, filespace, memspace);
	add_dbl_value(dt, "dt", file_id, filespace, memspace);
	add_int_value(lim, "lim", file_id, filespace, memspace);
	add_int_value(failed, "failed", file_id, filespace, memspace);
	add_dbl_value(Rin, "Rin", file_id, filespace, memspace);
	add_dbl_value(Rout, "Rout", file_id, filespace, memspace);
	add_dbl_value(hslope, "hslope", file_id, filespace, memspace);
	add_dbl_value(R0, "R0", file_id, filespace, memspace);

	H5Sclose(filespace);
	H5Sclose(memspace);

	for(d = 0; d < 3; d++) file_grid_dims[d] = fdims[d];
	file_grid_dims[3] = NPR;
	filespace = H5Screate_simple(4, file_grid_dims, NULL);
	for(d = 0; d < 3; d++) {
		file_start[d] = global_start[d];
		file_count[d] = mdims[d];
	}
	file_start[3] = 0;
	file_count[3] = NPR;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

	for(d = 0; d < 3; d++) {
		mem_grid_dims[d] = mdims[d] + 2*NG;
	}
	mem_grid_dims[3] = NPR;
	memspace = H5Screate_simple(4, mem_grid_dims, NULL);
	for(d = 0; d < 3; d++) mem_start[d] = 0;
	mem_start[3] = 0;
	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

	sprintf(name,"p");
	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	hid_t dset_id = H5Dcreate(file_id, name, H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &p[NG][NG][NG][0]);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	H5Sclose(memspace);
	H5Sclose(filespace);

	H5Fflush(file_id, H5F_SCOPE_GLOBAL);
	H5Fclose(file_id);

	return;
}

void restart_read(char *fname) 
{

	//static int restart_id = 0;
	int d;
	hsize_t file_grid_dims[4],file_start[4],file_count[4];
	hsize_t mem_grid_dims[4],mem_start[4],one;
	char *name = "p";
	int fdims[3] = {N1TOT, N2TOT, N3TOT};
	int mdims[3] = {N1, N2, N3};

	if(mpi_io_proc()) {
		fprintf(stderr,"restarting from %s...", fname);
	}

	one = 1;

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
	hid_t file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
	H5Pclose(plist_id);

	// first write some info about the current state of the run (e.g. time, istep, ...)
	file_grid_dims[0] = 1;
	hid_t filespace = H5Screate_simple(1, file_grid_dims, NULL);
	file_start[0] = 0;
	file_count[0] = 1;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);
	hid_t memspace  = H5Screate_simple(1, &one, NULL);

	get_dbl_value(&t, "Time", file_id, filespace, memspace);
	get_int_value(&nstep, "Step", file_id, filespace, memspace);
	get_int_value(&dump_id, "Next_Dump_ID", file_id, filespace, memspace);
	get_dbl_value(&tf, "Final_Time", file_id, filespace, memspace);
	get_dbl_value(&a, "a", file_id, filespace, memspace);
	get_dbl_value(&gam, "gam", file_id, filespace, memspace);
	get_dbl_value(&cour, "cour", file_id, filespace, memspace);
	get_dbl_value(&DTd, "DTd", file_id, filespace, memspace);
	get_dbl_value(&DTl, "DTl", file_id, filespace, memspace);
	get_dbl_value(&DTi, "DTi", file_id, filespace, memspace);
	get_dbl_value(&DTp, "DTp", file_id, filespace, memspace);		
	get_int_value(&DTr, "DTr", file_id, filespace, memspace);
	get_int_value(&restart_id, "restart_id", file_id, filespace, memspace);
	get_dbl_value(&dt, "dt", file_id, filespace, memspace);
	get_int_value(&lim, "lim", file_id, filespace, memspace);
	get_int_value(&failed, "failed", file_id, filespace, memspace);
	get_dbl_value(&Rin, "Rin", file_id, filespace, memspace);
	get_dbl_value(&Rout, "Rout", file_id, filespace, memspace);
	get_dbl_value(&hslope, "hslope", file_id, filespace, memspace);
	get_dbl_value(&R0, "R0", file_id, filespace, memspace);

	H5Sclose(filespace);
	H5Sclose(memspace);

	for(d = 0; d < 3; d++) file_grid_dims[d] = fdims[d];
	file_grid_dims[3] = NPR;	// for vectors
	filespace = H5Screate_simple(4, file_grid_dims, NULL);
	for(d = 0; d< 3; d++) {
		file_start[d] = global_start[d];
		file_count[d] = mdims[d];
	}
	file_start[3] = 0;		
	file_count[3] = NPR;
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, file_start, NULL, file_count, NULL);

	for(d = 0; d < 3; d++) {
		mem_grid_dims[d] = mdims[d] + 2*NG;
	}
	mem_grid_dims[3] = NPR;
	memspace = H5Screate_simple(4, mem_grid_dims, NULL);
	for(d = 0; d< 3; d++) mem_start[d] = 0;
	mem_start[3] = 0;
	H5Sselect_hyperslab(memspace, H5S_SELECT_SET, mem_start, NULL, file_count, NULL);

	plist_id = H5Pcreate(H5P_DATASET_ACCESS);
	hid_t dset_id = H5Dopen(file_id, name, plist_id);
	H5Pclose(plist_id);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, &p[NG][NG][NG][0]);
	H5Dclose(dset_id);
	H5Pclose(plist_id);

	H5Sclose(memspace);
	H5Sclose(filespace);
	H5Fclose(file_id);

	MPI_Barrier(MPI_COMM_WORLD);
	if(mpi_io_proc()) {
		fprintf(stderr,"done\n");
	}

	return;
}

int restart_init()
{
	int i,j,k;
	char fname[256];
	FILE *fp = fopen("restarts/restart.last","r");
	if(fp == NULL) return 0;

	fscanf(fp,"%s\n", fname);
	restart_read(fname);

	set_grid();

	bound_prim(p);

	ZSLOOP(-NG, N1-1+NG, -NG, N2-1+NG, -NG, N3-1+NG) PLOOP ph[i][j][k][ip] = p[i][j][k][ip];

	return 1;
}
