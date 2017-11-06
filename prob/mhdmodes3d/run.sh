#!/bin/bash

if [ -d build ]; then
  cd build
fi

# Make absolutely sure we're setting affinities correctly!
#NP=$(( $(nproc) / 2 ))
NP=28
export OMP_NUM_THREADS=$NP
export GOMP_CPU_AFFINITY=0-$(( $NP - 1 ))
export OMP_PROC_BIND=true

numactl --interleave=all ./bhlight 2>&1 >out.txt
