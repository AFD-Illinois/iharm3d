#!/bin/bash

# Must be just a name for now
OUT_DIR=test_convergence

# In case we didn't clean up after test_restart_diffmpi
sed -i -e "s/N1CPU 1/N1CPU 2/g" parameters.h
sed -i -e "s/N2CPU 1/N2CPU 2/g" parameters.h
sed -i -e "s/N3CPU 1/N3CPU 4/g" parameters.h
export OMP_NUM_THREADS=2

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 16 32 64
do

  for i in 1 2 3
  do

    # Look, if anyone's got bright ideas for parameter management,
    # I'm all ears
    sed -i -e "s/512/$n/g" parameters.h
    sed -i -e "s/256/$n/g" parameters.h
    sed -i -e "s/128/$n/g" parameters.h
    sed -i -e "s/64/$n/g" parameters.h
    sed -i -e "s/32/$n/g" parameters.h
    sed -i -e "s/16/$n/g" parameters.h

    sed -i -e "s/nmode = 0/nmode = $i/g" param.dat
    sed -i -e "s/nmode = 1/nmode = $i/g" param.dat
    sed -i -e "s/nmode = 2/nmode = $i/g" param.dat
    sed -i -e "s/nmode = 3/nmode = $i/g" param.dat

    sleep 1

    make -f ../../makefile -j4 PROB=mpitest_modes
    mpirun -n 16 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/output_${n}_${i}.txt

    mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}_${i}
    rm -rf $OUT_DIR/restarts

  done
done

# Analysis and plots
cp -r test $OUT_DIR/
cd $OUT_DIR/test
python test_3D.py
