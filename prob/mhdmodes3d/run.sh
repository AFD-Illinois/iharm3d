#!/bin/bash

#python build.py

NMPI=32

if (( $NMPI == 1 ))
then
  export OMP_NUM_THREADS=$(( $(nproc) / 2 ))
  export GOMP_CPU_AFFINITY=0-$(( $OMP_NUM_THREADS - 1 ))
  export OMP_PROC_BIND=true
  numactl --interleave=all ./bhlight 2>&1 >out.txt
else
  export OMP_NUM_THREADS=1
  mpiexec -n $NMPI ./bhlight 2>&1 >out.txt
fi
