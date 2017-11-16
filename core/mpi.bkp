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

  numprocs = N1CPU*N2CPU*N3CPU;

  int cpudims[3] = {N1CPU, N2CPU, N3CPU};
/*
  // Make MPI communication periodic if required
  int isPeriodicX1, isPeriodicX2, isPeriodicX3;
  #if X1L_BOUND == PERIODIC && X1R_BOUND == PERIODIC
  isPeriodicX1 = 1;
  #else
  isPeriodicX1 = 0;
  #endif
  #if X2L_BOUND == PERIODIC && X2R_BOUND == PERIODIC
  isPeriodicX2 = 1;
  #else
  isPeriodicX2 = 0;
  #endif
  #if X3L_BOUND == PERIODIC && X3R_BOUND == PERIODIC
  isPeriodicX3 = 1;
  #else
  isPeriodicX3 = 0;
  #endif
  int periodic[3] = {isPeriodicX1, isPeriodicX2, isPeriodicX3};
*/ // THIS BREAKS MPI INITIALIZATION FOR SOME REASON

  int periodic[3] = {0,0,1};

  // Set up communicator for Cartesian processor topology
  MPI_Cart_create(MPI_COMM_WORLD, 3, cpudims, periodic, 1, &comm);

  MPI_Comm_rank(comm, &rank);
  MPI_Cart_coords(comm, rank, 3, coord);

  // Find the ranks of neighbors, including edge/corner neighbors
  int n[3];
  for (int i = -1; i < 2; i++) {
    n[0] = coord[0] + i;
    for (int j = -1; j < 2; j++) {
      n[1] = coord[1] + j;
      for (int k = -1; k < 2; k++) {
        n[2] = coord[2] + k;
        if ((n[0] < 0 || n[0] >= N1CPU) || (n[1] < 0 || n[1] >= N2CPU)) {
          neighbors[i+1][j+1][k+1] = MPI_PROC_NULL;
        } else {
          MPI_Cart_rank(comm, n, &neighbors[i+1][j+1][k+1]);
        }
      }
    }
  }

  // Start and stop in global index space
  int sdims[3] = {N1TOT, N2TOT, N3TOT};
  for (int d = 0; d < 3; d++) {
    global_start[d] = coord[d] * sdims[d]/cpudims[d];
    global_stop[d] = (coord[d]+1) * sdims[d]/cpudims[d];
  }

  // Make MPI datatypes
  int max_count = MY_MAX(N2*NG, MY_MAX((N1+2*NG)*NG, (N1+2*NG)*(N2+2*NG)));
  int *blocks = malloc(max_count*sizeof(int));
  MPI_Aint *offsets = malloc(max_count*sizeof(MPI_Aint));
  MPI_Datatype *types = malloc(max_count*sizeof(MPI_Datatype));

  // First spatial direction face
  int count = N2*NG;
  for (int i = 0; i < count ; i++) {
    blocks[i] = N3*NVAR;
    types[i] = MPI_DOUBLE;
  }
  for (int i = 0; i < NG; i++) {
    for (int j = 0; j < N2; j++) {
      offsets[i*N2 + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NVAR*sizeof(double) +
                          j*(N3 + 2*NG)*NVAR*sizeof(double);
    }
  }
  MPI_Type_create_struct(count, blocks, offsets, types, &face_type[0]);
  MPI_Type_commit(&face_type[0]);

  // Second spatial direction face
  count = (N1+2*NG)*NG;
  for (int i = 0; i < count ; i++) {
    blocks[i] = (N3+2*NG)*NVAR;
    types[i] = MPI_DOUBLE;
  }
  for (int i = 0; i < N1+2*NG; i++) {
    for (int j = 0; j < NG; j++) {
      offsets[NG*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NVAR*sizeof(double) +
                          j*(N3 + 2*NG)*NVAR*sizeof(double);
    }
  }
  MPI_Type_create_struct(count, blocks, offsets, types, &face_type[1]);
  MPI_Type_commit(&face_type[1]);

  // Third spatial direction face
  count = (N1+2*NG)*(N2+2*NG);
  for (int i = 0; i < count; i++) {
    blocks[i] = NG*NVAR;
    types[i] = MPI_DOUBLE;
  }
  for (int i = 0; i < (N1+2*NG); i++) {
    for (int j = 0; j < (N2+2*NG); j++) {
      offsets[N2*i + j] = i*(N2 + 2*NG)*(N3 + 2*NG)*NVAR*sizeof(double) +
                          j*(N3 + 2*NG)*NVAR*sizeof(double);
    }
  }
  MPI_Type_create_struct(count, blocks, offsets, types, &face_type[2]);
  MPI_Type_commit(&face_type[2]);

  free(blocks);
  free(offsets);
  free(types);
}

// Share face data
void sync_mpi_boundaries(struct FluidState *S)
{
  #if 0
  MPI_Status status;

  // First send right/receive left
  MPI_Sendrecv(&Pr[N1][NG][NG][0], 1, face_type[0], neighbors[2][1][1], 0,
    &Pr[0][NG][NG][0], 1, face_type[0], neighbors[0][1][1], 0, comm, &status);
  // Now the other way
  MPI_Sendrecv(&Pr[NG][NG][NG][0], 1, face_type[0], neighbors[0][1][1], 1,
    &Pr[N1+NG][NG][NG][0], 1, face_type[0], neighbors[2][1][1], 1, comm, 
    &status);

  #if N2 > 1
  // Other directions
  MPI_Sendrecv(&Pr[0][N2][0][0], 1, face_type[1], neighbors[1][2][1], 2,
    &Pr[0][0][0][0], 1, face_type[1], neighbors[1][0][1], 2, comm, &status);
  MPI_Sendrecv(&Pr[0][NG][0][0], 1, face_type[1], neighbors[1][0][1], 3,
    &Pr[0][N2+NG][0][0], 1, face_type[1], neighbors[1][2][1], 3, comm, &status);
  #endif

  #if N3 > 1
  MPI_Sendrecv(&Pr[0][0][N3][0], 1, face_type[2], neighbors[1][1][2], 4,
    &Pr[0][0][0][0], 1, face_type[2], neighbors[1][1][0], 4, comm, &status);
  MPI_Sendrecv(&Pr[0][0][NG][0], 1, face_type[2], neighbors[1][1][0], 5,
    &Pr[0][0][N3+NG][0], 1, face_type[2], neighbors[1][1][2], 5, comm, &status);
  #endif
#endif
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

