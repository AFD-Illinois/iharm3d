#!/bin/bash

# Must be just a name for now
OUT_DIR=test_convergence

# In case we didn't clean up after test_restart_diffmpi
sed -i -e "s/N1CPU [0-9]\+/N1CPU 2/g" parameters.h
sed -i -e "s/N2CPU [0-9]\+/N2CPU 2/g" parameters.h
export OMP_NUM_THREADS=8

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 16 32 64 128
do

  # Look, if anyone's got bright ideas for parameter management,
  # I'm all ears
  sed -i -e "s/N1TOT [0-9]\+/N1TOT $n/g" parameters.h
  sed -i -e "s/N2TOT [0-9]\+/N2TOT $n/g" parameters.h

  sleep 1

  make -f ../../makefile -j4 PROB=$(basename "$PWD")
  mpirun -n 4 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/output_${n}.txt

  mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}
  rm -rf $OUT_DIR/restarts
done

# Analysis and plots
mkdir $OUT_DIR/plots
cp plot_convergence.py $OUT_DIR/plots
cd $OUT_DIR/plots
python plot_convergence.py
