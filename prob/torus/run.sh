#!/bin/bash

OUT_DIR=${1:-/mnt/data2/bprather/default_dumps/}
#OUT_DIR=${1:-.}
NMPI=8

#rm -rf $OUT_DIR ; mkdir $OUT_DIR

. ../../../mpi-load.sh

python build.py

if (( $NMPI == 1 ))
then
  export OMP_NUM_THREADS=$(( $(nproc --all) / 2 ))
  export GOMP_CPU_AFFINITY=0-$(( $OMP_NUM_THREADS - 1 ))
  export OMP_PROC_BIND=true
  echo "Using first $OMP_NUM_THREADS threads"

  numactl --interleave=all ./bhlight -o $OUT_DIR
else
  export OMP_NUM_THREADS=4
  echo "Using $NMPI local MPI processes"

  mpiexec -n $NMPI ./bhlight -o $OUT_DIR
fi
