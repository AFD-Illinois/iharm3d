/******************************************************************************
 *                                                                            *
 * MPI.C                                                                      *
 *                                                                            *
 * HANDLES COMMUNICATION ACROSS MPI NODES                                     *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

void init_mpi()
{
  int sdims[3] = {N1TOT, N2TOT, N3TOT};
  for (int d = 0; d < 3; d++) {
    global_start[d] = 0;
    global_stop[d] = sdims[d];
  }
}

// Share face data
void sync_mpi_boundaries(struct FluidState *S)
{
}

int mpi_nprocs() {
  return 1;
}

double mpi_max(double f)
{
  return f;
}

double mpi_min(double f)
{
  return f;
}

int mpi_io_proc()
{
  return 1;
}

void mpi_int_broadcast(int *val)
{
}

void mpi_dbl_broadcast(double *val)
{
}

double mpi_io_reduce(double val)
{
  return val;
}

double mpi_io_max(double val)
{
  return val;
}

int mpi_myrank() {
  return 0;
}

