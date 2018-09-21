#!/bin/bash

source ../test_common.sh

# Must be just a name for now
OUT_DIR=results_modes

# Initial clean and make of work area
BASEDIR=../../..
rm -rf build_archive param.dat harm
make -f $BASEDIR/makefile -j4 PROB=mhdmodes

# In case we didn't clean up after test_restart_diffmpi
set_compile_int N1CPU 2
set_compile_int N2CPU 2
set_compile_int N3CPU 4
# Gotta test OpenMP somewhere... sorry laptops
export OMP_NUM_THREADS=2

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 16 32 64
do

  set_compile_int N1TOT $n
  set_compile_int N2TOT $n
  set_compile_int N3TOT $n

  sleep 1

  make -f $BASEDIR/makefile -j4 PROB=mhdmodes

  for i in 1 2 3
  do

    set_run_int nmode $i

    sleep 1

    echo "Running size $n mode $i..."
    mpirun -n 16 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/output_${n}_${i}.txt
    echo "Done!"

    mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}_${i}
    rm -rf $OUT_DIR/restarts

  done
done

# Run analysis automatically
./ana_convergence.sh modes
