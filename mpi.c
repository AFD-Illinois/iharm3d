
#include "decs.h"
#include "mpi.h"

static MPI_Comm comm;
static int neighbors[3][3][3];
static MPI_Datatype face_type[3];
static MPI_Datatype edge_type[3];
static MPI_Datatype corn_type;
static int rank;
static int numprocs;

void init_mpi() 
{
	//int i,j,k;
	int coord[3];

	numprocs = N1CPU*N2CPU*N3CPU;

	int cpudims[3] = {N1CPU, N2CPU, N3CPU};
	int periodic[3] = {0, 0, 1};								// this makes the MPI communication periodic in direction 3
	MPI_Cart_create(MPI_COMM_WORLD, 3, cpudims, periodic, 1, &comm);	// set up a communicator for a Cartesian processor topology

	MPI_Comm_rank(comm, &rank);
	MPI_Cart_coords(comm, rank, 3, coord);		// find my coordinates in the Cartesian processor topology

	// find the ranks of neighbors --- need edge/corner neighbors, too.  Can't use MPI_Cart_shift for that.
	int n[3];
	for(int i = -1; i < 2; i++) {
		n[0] = coord[0] + i;
		for(int j = -1; j < 2; j++) {
			n[1] = coord[1] + j;
			for(int k = -1; k < 2; k++) {
				n[2] = coord[2] + k;
				if((n[0] < 0 || n[0] >= N1CPU)
					|| (n[1] < 0 || n[1] >= N2CPU)) {
					neighbors[i+1][j+1][k+1] = MPI_PROC_NULL;
				} else {
					MPI_Cart_rank(comm, n, &neighbors[i+1][j+1][k+1]);
				}
			}
		}
	}

	int sdims[3] = {N1TOT, N2TOT, N3TOT};
	for(int d = 0; d < 3; d++) {
		global_start[d] = coord[d] * sdims[d]/cpudims[d];		// start and stop in a global index space
		global_stop[d] = (coord[d]+1) * sdims[d]/cpudims[d];
	}

	/* make MPI datatypes */
	int max_count = MY_MAX(N2*NG, MY_MAX(N1*NG, N1*N2));
	int *blocks = malloc(max_count*sizeof(int));
	MPI_Aint *offsets = malloc(max_count*sizeof(MPI_Aint));
	MPI_Datatype *types = malloc(max_count*sizeof(MPI_Datatype));

	// first spatial direction face
	int count = N2*NG;
	for(int i = 0; i < count ; i++) {
		blocks[i] = N3*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < NG; i++) {
		for(int j = 0; j < N2; j++) {
			offsets[i*N2 + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &face_type[0]);
	MPI_Type_commit(&face_type[0]);

	// next spatial direction face
	count = N1*NG;
	for(int i = 0; i < count ; i++) {
		blocks[i] = N3*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < N1; i++) {
		for(int j = 0; j < NG; j++) {
			offsets[NG*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &face_type[1]);
	MPI_Type_commit(&face_type[1]);

	// last spatial direction face
	count = N1*N2;
	for(int i = 0; i < count; i++) {
		blocks[i] = NG*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < N1; i++) {
		for(int j = 0; j < N2; j++) {
			offsets[N2*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &face_type[2]);
	MPI_Type_commit(&face_type[2]);

	// now edges so fluxes can be calculated that are used in CT
	// first [N1][NG][NG]
	count = N1*NG;
	for(int i = 0; i < count; i++) {
		blocks[i] = NG*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < N1; i++) {
		for(int j = 0; j < NG; j++) {
			offsets[NG*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &edge_type[0]);
	MPI_Type_commit(&edge_type[0]);

	count = N2*NG;
	for(int i = 0; i < count; i++) {
		blocks[i] = NG*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < NG; i++) {
		for(int j = 0; j < N2; j++) {
			offsets[N2*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &edge_type[1]);
	MPI_Type_commit(&edge_type[1]);


	count = NG*NG;
	for(int i = 0; i < count; i++) {
		blocks[i] = N3*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < NG; i++) {
		for(int j = 0; j < NG; j++) {
			offsets[NG*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &edge_type[2]);
	MPI_Type_commit(&edge_type[2]);

	// finally we need to share corner data so make its datatype
	count = NG*NG;
	for(int i = 0; i < count; i++) {
		blocks[i] = NG*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < NG; i++) {
		for(int j = 0; j < NG; j++) {
			offsets[NG*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &corn_type);
	MPI_Type_commit(&corn_type);

	free(blocks);
	free(offsets);
	free(types);

}

void sync_mpi_boundaries(grid_prim_type pr)
{
	MPI_Status status;

	/* Share face data */
	// first send right/receive left
	MPI_Sendrecv(&pr[N1][NG][NG][0], 1, face_type[0], neighbors[2][1][1], 0,
					&pr[0][NG][NG][0], 1, face_type[0], neighbors[0][1][1], 0, comm, &status);
	// now the other way
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, face_type[0], neighbors[0][1][1], 1,
					&pr[N1+NG][NG][NG][0], 1, face_type[0], neighbors[2][1][1], 1, comm, &status);

	// other directions
	MPI_Sendrecv(&pr[NG][N2][NG][0], 1, face_type[1], neighbors[1][2][1], 2,
					&pr[NG][0][NG][0], 1, face_type[1], neighbors[1][0][1], 2, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, face_type[1], neighbors[1][0][1], 3,
					&pr[NG][N2+NG][NG][0], 1, face_type[1], neighbors[1][2][1], 3, comm, &status);

	MPI_Sendrecv(&pr[NG][NG][N3][0], 1, face_type[2], neighbors[1][1][2], 4,
					&pr[NG][NG][0][0], 1, face_type[2], neighbors[1][1][0], 4, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, face_type[2], neighbors[1][1][0], 5,
					&pr[NG][NG][N3+NG][0], 1, face_type[2], neighbors[1][1][2], 5, comm, &status);

	/* Share edge data (needed for CT) */
	
	// first [*][NG][NG]
	MPI_Sendrecv(&pr[NG][N2][N3][0], 1, edge_type[0], neighbors[1][2][2], 6,
					&pr[NG][0][0][0], 1, edge_type[0], neighbors[1][0][0], 6, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, edge_type[0], neighbors[1][0][0], 7,
					&pr[NG][N2+NG][N3+NG][0], 1, edge_type[0], neighbors[1][2][2], 7, comm, &status);
	MPI_Sendrecv(&pr[NG][N2][NG][0], 1, edge_type[0], neighbors[1][2][0], 8,
					&pr[NG][0][N3+NG][0], 1, edge_type[0], neighbors[1][0][2], 8, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][N3][0], 1, edge_type[0], neighbors[1][0][2], 9,
					&pr[NG][N2+NG][0][0], 1, edge_type[0], neighbors[1][2][0], 9, comm, &status);

	// now [NG][*][NG]
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, edge_type[1], neighbors[0][1][0], 10, 
					&pr[N1+NG][NG][N3+NG][0], 1, edge_type[1], neighbors[2][1][2], 10, comm, &status);
	MPI_Sendrecv(&pr[N1][NG][N3][0], 1, edge_type[1], neighbors[2][1][2], 11, 
					&pr[0][NG][0][0], 1, edge_type[1], neighbors[0][1][0], 11, comm, &status);
	MPI_Sendrecv(&pr[N1][NG][NG][0], 1, edge_type[1], neighbors[2][1][0], 12,
					&pr[0][NG][N3+NG][0], 1, edge_type[1], neighbors[0][1][2], 12, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][N3][0], 1, edge_type[1], neighbors[0][1][2], 13,
					&pr[N1+NG][NG][0][0], 1, edge_type[1], neighbors[2][1][0], 13, comm, &status);

	// now [NG][NG][*]
	MPI_Sendrecv(&pr[N1][NG][NG][0], 1, edge_type[2], neighbors[2][0][1], 14,
					&pr[0][N2+NG][NG][0], 1, edge_type[2], neighbors[0][2][1], 14, comm, &status);
	MPI_Sendrecv(&pr[NG][N2][NG][0], 1, edge_type[2], neighbors[0][2][1], 15,
					&pr[N1+NG][0][NG][0], 1, edge_type[2], neighbors[2][0][1], 15, comm, &status);
	MPI_Sendrecv(&pr[N1][N2][NG][0], 1, edge_type[2], neighbors[2][2][1], 16,
					&pr[0][0][NG][0], 1, edge_type[2], neighbors[0][0][1], 16, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, edge_type[2], neighbors[0][0][1], 17,
					&pr[N1+NG][N2+NG][NG][0], 1, edge_type[2], neighbors[2][2][1], 17, comm, &status);


	/* Finally, share corner data (again for CT) */
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, corn_type, neighbors[0][0][0], 18,
					&pr[N1+NG][N2+NG][N3+NG][0], 1, corn_type, neighbors[2][2][2], 18, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][N3][0], 1, corn_type, neighbors[0][0][2], 19,
					&pr[N1+NG][N2+NG][0][0], 1, corn_type, neighbors[2][2][0], 19, comm, &status);
	MPI_Sendrecv(&pr[NG][N2][NG][0], 1, corn_type, neighbors[0][2][0], 20,
					&pr[N1+NG][0][N3+NG][0], 1, corn_type, neighbors[2][0][2], 20, comm, &status);
	MPI_Sendrecv(&pr[NG][N2][N3][0], 1, corn_type, neighbors[0][2][2], 21,
					&pr[N1+NG][0][0][0], 1, corn_type, neighbors[2][0][0], 21, comm, &status);
	MPI_Sendrecv(&pr[N1][NG][NG][0], 1, corn_type, neighbors[2][0][0], 22,
					&pr[0][N2+NG][N3+NG][0], 1, corn_type, neighbors[0][2][2], 22, comm, &status);
	MPI_Sendrecv(&pr[N1][NG][N3][0], 1, corn_type, neighbors[2][0][2], 23,
					&pr[0][N2+NG][0][0], 1, corn_type, neighbors[0][2][0], 23, comm, &status);
	MPI_Sendrecv(&pr[N1][N2][NG][0], 1, corn_type, neighbors[2][2][0], 24,
					&pr[0][0][N3+NG][0], 1, corn_type, neighbors[0][0][2], 24, comm, &status);
	MPI_Sendrecv(&pr[N1][N2][N3][0], 1, corn_type, neighbors[2][2][2], 25,
					&pr[0][0][0][0], 1, corn_type, neighbors[0][0][0], 25, comm, &status);


	/* done!  26 messages = 6 (large) face messages + 12 (medium) edge messages + 8 (small) corner messages*/

	
}

int mpi_nprocs() {
	return numprocs;
}

double mpi_max(double f) 
{
	double fmax;
	MPI_Allreduce(&f, &fmax, 1, MPI_DOUBLE, MPI_MAX, comm);
	return fmax;
}

double mpi_min(double f) 
{
	double fmin;
	MPI_Allreduce(&f, &fmin, 1, MPI_DOUBLE, MPI_MIN, comm);
	return fmin;
}

int mpi_io_proc() 
{
	return (rank == 0 ? 1 : 0);
}

void mpi_int_broadcast(int *val)
{
	MPI_Bcast(val, 1, MPI_INT, 0, comm);
}

void mpi_dbl_broadcast(double *val)
{
	MPI_Bcast(val, 1, MPI_DOUBLE, 0, comm);
}

double mpi_io_reduce(double val) {

	double local;
	MPI_Reduce(&val, &local, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
	return local;
}

double mpi_io_max(double val) {

	double local;
	MPI_Reduce(&val, &local, 1, MPI_DOUBLE, MPI_MAX, 0, comm);
	return local;
}

int mpi_myrank() {
	return rank;
}
