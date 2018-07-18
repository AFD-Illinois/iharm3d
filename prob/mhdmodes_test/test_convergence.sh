#!/bin/bash

# Must be just a name for now
OUT_DIR=test_convergence

# In case we didn't clean up after test_restart_diffmpi
sed -i -e "s/N1CPU [0-9]\+/N1CPU 2/g" parameters.h
sed -i -e "s/N2CPU [0-9]\+/N2CPU 2/g" parameters.h
sed -i -e "s/N3CPU [0-9]\+/N3CPU 4/g" parameters.h
export OMP_NUM_THREADS=2

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 16 32 64
do

  # Look, if anyone's got bright ideas for parameter management,
  # I'm all ears
  sed -i -e "s/N1TOT [0-9]\+/N1TOT $n/g" parameters.h
  sed -i -e "s/N2TOT [0-9]\+/N2TOT $n/g" parameters.h
  sed -i -e "s/N3TOT [0-9]\+/N3TOT $n/g" parameters.h

  sleep 1

  make -f ../../makefile -j4 PROB=$(basename "$PWD")

  for i in 1 2 3
  do

    sed -i -e "s/nmode = [0-3]/nmode = $i/g" param.dat

    sleep 1

    mpirun -n 16 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/output_${n}_${i}.txt

    mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}_${i}
    rm -rf $OUT_DIR/restarts

  done
done

# Analysis and plots
mkdir $OUT_DIR/plots
cp plot_convergence.py $OUT_DIR/plots
cd $OUT_DIR/plots
python plot_convergence.py
