
#include "decs.h"
#include "mpi.h"

static MPI_Comm comm;
static int neighbors[3][2];
static MPI_Datatype bound_type[3];
static int rank;

void init_mpi() 
{
	//int i,j,k;
	int src,coord[3];

	int cpudims[3] = {N1CPU, N2CPU, N3CPU};
	int periodic[3] = {0, 0, 1};								// this makes the MPI communication periodic in direction 3
	MPI_Cart_create(MPI_COMM_WORLD, 3, cpudims, periodic, 1, &comm);	// set up a communicator for a Cartesian processor topology

	for(int d = 0; d < 3; d++) {
		MPI_Cart_shift(comm, d, -1, &src, &neighbors[d][0]);	// neighbor to the left in direction d
		MPI_Cart_shift(comm, d, 1, &src, &neighbors[d][1]);		// right neighbor
	}

	MPI_Comm_rank(comm, &rank);
	MPI_Cart_coords(comm, rank, 3, coord);		// find my coordinates in the Cartesian processor topology

	int sdims[3] = {N1TOT, N2TOT, N3TOT};
	for(int d = 0; d < 3; d++) {
		global_start[d] = coord[d] * sdims[d]/cpudims[d];		// start and stop in a global index space
		global_stop[d] = (coord[d]+1) * sdims[d]/cpudims[d];
	}

	/* make MPI datatypes */
	// first spatial direction
	int count = N2*NG;
	int *blocks = malloc(count*sizeof(int));
	MPI_Aint *offsets = malloc(count*sizeof(MPI_Aint));
	MPI_Datatype *types = malloc(count*sizeof(MPI_Datatype));
	for(int i = 0; i < count ; i++) {
		blocks[i] = N3*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < NG; i++) {
		for(int j = 0; j < N2; j++) {
			offsets[i*N2 + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &bound_type[0]);
	MPI_Type_commit(&bound_type[0]);
	free(blocks);
	free(offsets);
	free(types);

	// next spatial direction
	count = N1*NG;
	blocks = malloc(count*sizeof(int));
	offsets = malloc(count*sizeof(MPI_Aint));
	types = malloc(count*sizeof(MPI_Datatype));
	for(int i = 0; i < count ; i++) {
		blocks[i] = N3*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < N1; i++) {
		for(int j = 0; j < NG; j++) {
			offsets[NG*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &bound_type[1]);
	MPI_Type_commit(&bound_type[1]);
	free(blocks);
	free(offsets);
	free(types);

	// last spatial direction
	count = N1*N2;
	blocks = malloc(count*sizeof(int));
	offsets = malloc(count*sizeof(MPI_Aint));
	types = malloc(count*sizeof(MPI_Datatype));
	for(int i = 0; i < count; i++) {
		blocks[i] = NG*NPR;
		types[i] = MPI_DOUBLE;
	}
	for(int i = 0; i < N1; i++) {
		for(int j = 0; j < N2; j++) {
			offsets[N2*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NPR*sizeof(double) + j*(N3 + 2*NG)*NPR*sizeof(double);
		}
	}
	MPI_Type_struct(count, blocks, offsets, types, &bound_type[2]);
	MPI_Type_commit(&bound_type[2]);
	free(blocks);
	free(offsets);
	free(types);

/*	ZSLOOP(-NG, N1-1+NG, -NG, N2-1+NG, -NG, N3-1+NG) PLOOP p[i][j][k][ip] = 10*rank + ip;

	sync_mpi_boundaries(p);

	int n = neighbors[1][1];
	if(n != MPI_PROC_NULL) {
	for(i = NG; i < N1+NG; i++) {
	for(k = NG; k < N3+NG; k++) {
		for(j=N2+NG;j<N2+2*NG;j++) PLOOP {
			if(p[i][j][k][ip] != 10*n + ip) {
				fprintf(stderr,"problem with MPI on %d  %d %g %d\n", rank, n, p[i][j][k][ip], 10*n + ip);
			}
		}
	}
	}
	}
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
	exit(123);
*/
}

void sync_mpi_boundaries(grid_prim_type pr)
{
	MPI_Status status;

	// first send right/receive left
	MPI_Sendrecv(&pr[N1][NG][NG][0], 1, bound_type[0], neighbors[0][1], 0,
					&pr[0][NG][NG][0], 1, bound_type[0], neighbors[0][0], 0, comm, &status);
	// now the other way
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, bound_type[0], neighbors[0][0], 1,
					&pr[N1+NG][NG][NG][0], 1, bound_type[0], neighbors[0][1], 1, comm, &status);

	// other directions
	MPI_Sendrecv(&pr[NG][N2][NG][0], 1, bound_type[1], neighbors[1][1], 2,
					&pr[NG][0][NG][0], 1, bound_type[1], neighbors[1][0], 2, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, bound_type[1], neighbors[1][0], 3,
					&pr[NG][N2+NG][NG][0], 1, bound_type[1], neighbors[1][1], 3, comm, &status);

	MPI_Sendrecv(&pr[NG][NG][N3][0], 1, bound_type[2], neighbors[2][1], 4,
					&pr[NG][NG][0][0], 1, bound_type[2], neighbors[2][0], 4, comm, &status);
	MPI_Sendrecv(&pr[NG][NG][NG][0], 1, bound_type[2], neighbors[2][0], 5,
					&pr[NG][NG][N3+NG][0], 1, bound_type[2], neighbors[2][1], 5, comm, &status);
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
