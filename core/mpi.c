/******************************************************************************
 *                                                                            *
 * MPI.C                                                                      *
 *                                                                            *
 * HANDLES COMMUNICATION ACROSS MPI NODES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static MPI_Comm comm;
static int neighbors[3][3][3];
static MPI_Datatype face_type[3];
static int rank;
static int numprocs;

void init_mpi()
{
  int coord[3];
  numprocs = N3CPU*N2CPU*N1CPU;

  int cpudims[3] = {N3CPU, N2CPU, N1CPU};

  // Make MPI communication periodic if required
  int periodic[3] = {X3L_BOUND == PERIODIC && X3R_BOUND == PERIODIC,
      X2L_BOUND == PERIODIC && X2R_BOUND == PERIODIC,
      X1L_BOUND == PERIODIC && X1R_BOUND == PERIODIC};

  // Set up communicator for Cartesian processor topology
  // Use X3,2,1 ordering
  MPI_Cart_create(MPI_COMM_WORLD, 3, cpudims, periodic, 1, &comm);

  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 3, coord);

  // Find the ranks of neighbors, including edge/corner neighbors
  int n[3];
  for (int k = -1; k < 2; k++) {
    n[0] = coord[0] + k;
    for (int j = -1; j < 2; j++) {
      n[1] = coord[1] + j;
      for (int i = -1; i < 2; i++) {
        n[2] = coord[2] + i;
        if (((n[0] < 0 || n[0] >= N3CPU) && !periodic[0]) ||
            ((n[1] < 0 || n[1] >= N2CPU) && !periodic[1]) ||
            ((n[2] < 0 || n[2] >= N1CPU) && !periodic[2])) {
          neighbors[k+1][j+1][i+1] = MPI_PROC_NULL;
        } else {
          MPI_Cart_rank(comm, n, &neighbors[k+1][j+1][i+1]);
        }
      }
    }
  }

  // Diagnostic for processor topology
//  int me = mpi_myrank();
//  fprintf(stderr,"Process %d topology:\n", me);
//  fprintf(stderr,"%d in X:[%d]\t[%d]\t[%d]\n", me, neighbors[1][1][0], neighbors[1][1][1], neighbors[1][1][2]);
//  fprintf(stderr,"%d in Y:[%d]\t[%d]\t[%d]\n", me, neighbors[1][0][1], neighbors[1][1][1], neighbors[1][2][1]);
//  fprintf(stderr,"%d in Z:[%d]\t[%d]\t[%d]\n", me, neighbors[0][1][1], neighbors[1][1][1], neighbors[2][1][1]);

  // Start and stop in global index space
  // These are in the usual X1,2,3 order, or things would get /very confusing/
  int sdims[3] = {N3TOT, N2TOT, N1TOT};
  for (int d = 0; d < 3; d++) {
    global_start[2-d] = coord[d] * sdims[d]/cpudims[d];
    global_stop[2-d] = (coord[d]+1) * sdims[d]/cpudims[d];
  }
  printf("Process %d has X,Y,Z space [%d-%d, %d-%d, %d-%d]\n", rank,
         global_start[0], global_stop[0],
         global_start[1], global_stop[1],
         global_start[2], global_stop[2]);

  // Make MPI datatypes
  int max_count = MY_MAX(NVAR, MY_MAX(NVAR*(N3+2*NG), NVAR*(N3+2*NG)*(N2+2*NG)));
  int *blocks = malloc(max_count*sizeof(int));
  MPI_Aint *offsets = malloc(max_count*sizeof(MPI_Aint));
  MPI_Datatype *types = malloc(max_count*sizeof(MPI_Datatype));

  // First spatial direction face
  // Need slice P[0-NVAR][i:i+NG][:][:]

  int count = NVAR;
  for (int i = 0; i < count ; i++) {
    blocks[i] = NG*(N2+2*NG)*(N1+2*NG);
    types[i] = MPI_DOUBLE;
  }

  PLOOP {
    offsets[ip] = ip*(N3+2*NG)*(N2+2*NG)*(N1+2*NG)*sizeof(double);
  }
  MPI_Type_create_struct(count, blocks, offsets, types, &face_type[0]);
  MPI_Type_commit(&face_type[0]);

  // Second spatial direction face
  // Slice P[0-NVAR][:][i:i+NG][:]

  count = NVAR*(N3+2*NG);
  for (int i = 0; i < count ; i++) {
    blocks[i] = NG*(N1+2*NG);
    types[i] = MPI_DOUBLE;
  }
  PLOOP {
    for (int j = 0; j < (N3+2*NG); j++) {
      offsets[ip*(N3 + 2*NG) + j] = ip*(N3+2*NG)*(N2+2*NG)*(N1+2*NG)*sizeof(double) +
        j*(N2+2*NG)*(N1+2*NG)*sizeof(double);
    }
  }
  MPI_Type_create_struct(count, blocks, offsets, types, &face_type[1]);
  MPI_Type_commit(&face_type[1]);

  // Third spatial direction face
  // Slice P[0-NVAR][:][:][i:i+NG]

  count = NVAR*(N3+2*NG)*(N2+2*NG);
  for (int i = 0; i < count; i++) {
    blocks[i] = NG;
    types[i] = MPI_DOUBLE;
  }
  PLOOP {
    for (int j = 0; j < (N3+2*NG); j++) {
      for (int k = 0; k < (N2+2*NG); k++) {
        offsets[ip*(N3+2*NG)*(N2+2*NG) + j*(N2+2*NG) + k] = ip*(N3 + 2*NG)*(N2 + 2*NG)*(N1 + 2*NG)*sizeof(double) +
          j*(N2 + 2*NG)*(N1+2*NG)*sizeof(double) +
          k*(N1+2*NG)*sizeof(double);
      }
    }
  }
  MPI_Type_create_struct(count, blocks, offsets, types, &face_type[2]);
  MPI_Type_commit(&face_type[2]);

  free(blocks);
  free(offsets);
  free(types);
}

// Share face data
MPI_Status* sync_mpi_boundaries(struct FluidState *S)
{
  MPI_Status *status = MPI_STATUS_IGNORE;

  // These are rather slow, so I use bounds.c implementations where possible
#if (N1 > 1) && (N1CPU > 1)
  // First send right/receive left
  MPI_Sendrecv(&(S->P[0][0][0][N1]), 1, face_type[2], neighbors[1][1][2], 0,
	       &(S->P[0][0][0][0]), 1, face_type[2], neighbors[1][1][0], 0, comm, status); // TODO check the status
  // And back
  MPI_Sendrecv(&(S->P[0][0][0][NG]), 1, face_type[2], neighbors[1][1][0], 1,
	       &(S->P[0][0][0][N1+NG]), 1, face_type[2], neighbors[1][1][2], 1, comm, status);
#endif

  // Other directions
#if (N2 > 1) && (N2CPU > 1)
  MPI_Sendrecv(&(S->P[0][0][N2][0]), 1, face_type[1], neighbors[1][2][1], 2,
	       &(S->P[0][0][0][0]), 1, face_type[1], neighbors[1][0][1], 2, comm, status);
  MPI_Sendrecv(&(S->P[0][0][NG][0]), 1, face_type[1], neighbors[1][0][1], 3,
	       &(S->P[0][0][N2+NG][0]), 1, face_type[1], neighbors[1][2][1], 3, comm, status);
#endif

#if (N3 > 1) && (N3CPU > 1)
  MPI_Sendrecv(&(S->P[0][N3][0][0]), 1, face_type[0], neighbors[2][1][1], 4,
	       &(S->P[0][0][0][0]), 1, face_type[0], neighbors[0][1][1], 4, comm, status);
  MPI_Sendrecv(&(S->P[0][NG][0][0]), 1, face_type[0], neighbors[0][1][1], 5,
	       &(S->P[0][N3+NG][0][0]), 1, face_type[0], neighbors[2][1][1], 5, comm, status);
#endif

  return status;
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
