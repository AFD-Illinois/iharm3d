#!/bin/bash


OUT_DIR=${1:-.}
NMPI=${NMPI:-1}

echo "Outputting to $OUT_DIR"

if (( $NMPI == 1 ))
then
  export OMP_NUM_THREADS=$(( $(nproc --all) / 2 ))
  export GOMP_CPU_AFFINITY=0-$(( $OMP_NUM_THREADS - 1 ))
  echo "Using first $OMP_NUM_THREADS threads"

  numactl --interleave=all ./bhlight -o $OUT_DIR

else
  export OMP_NUM_THREADS=$(( $(nproc --all) / $NMPI ))
  echo "Using $NMPI local MPI processes with $OMP_NUM_THREADS threads"

  mpiexec -n $NMPI ./bhlight -o $OUT_DIR
fi
