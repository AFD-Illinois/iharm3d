#!/bin/bash

#rm -rf dumps/ restarts/

. ../../../mpi-load.sh

python build.py

# Number of MPI Processes
NMPI=32

if (( $NMPI == 1 ))
then
  export OMP_NUM_THREADS=$(( $(nproc) / 2 ))
  export GOMP_CPU_AFFINITY=0-$(( $OMP_NUM_THREADS - 1 ))
  export OMP_PROC_BIND=true
  numactl --interleave=all ./bhlight
else
  export OMP_NUM_THREADS=1
  mpiexec -n $NMPI ./bhlight
  #mpiexec -n $NMPI valgrind ./bhlight
fi
