#!/bin/bash

source ../test_common.sh

# Must be just a name for now
OUT_DIR=results_bondi

# Initial clean and make of work area
rm -rf build_archive param.dat harm
make_harm_here bondi

# In case we didn't clean up after test_restart_diffmpi
set_cpu_topo 2 2 1
# Try for some efficiency
#export OMP_NUM_THREADS=$(( $(nproc --all) / 4 ))

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 32 64 128 256
do

  set_problem_size $n $n 1

  sleep 1

  make_harm_here bondi
  echo "Running $n-size bondi problem w/4 MPI procs, $OMP_NUM_THREADS threads each"
  mpirun -n 4 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/output_${n}.txt
  echo "Done!"

  mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}
  rm -rf $OUT_DIR/restarts
done

# Run analysis automatically
./ana_convergence.sh bondi
