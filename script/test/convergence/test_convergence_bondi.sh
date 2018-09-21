#!/bin/bash

source ../test_common.sh

# Must be just a name for now
OUT_DIR=results_bondi

# Initial clean and make of work area
BASEDIR=../../..
rm -rf build_archive param.dat harm
make -f $BASEDIR/makefile -j4 PROB=bondi

# In case we didn't clean up after test_restart_diffmpi
set_compile_int N1CPU 2
set_compile_int N2CPU 2
set_compile_int N3CPU 1
# Try for some efficiency
export OMP_NUM_THREADS=$(( $(nproc --all) / 4 ))

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 16 32 64 128
do

  set_compile_int N1TOT $n
  set_compile_int N2TOT $n

  sleep 1

  make -f $BASEDIR/makefile -j4 PROB=bondi
  mpirun -n 4 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/output_${n}.txt

  mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}
  rm -rf $OUT_DIR/restarts
done

# Run analysis automatically
./ana_convergence.sh bondi
