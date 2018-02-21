#!/bin/bash

source ../../../mpi-load.sh

# In case we didn't clean up after test_restart_diffmpi
#sed -i -e "s/N2CPU 1/N2CPU 2/g" parameters.h
#sed -i -e "s/N3CPU 1/N3CPU 4/g" parameters.h
export NMPI=1

# Must be just a name for now
OUT_DIR=test_convergence

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

for n in 16 32 64
do

  for i in 0 1 2 3
  do

    # I can't believe I'm doing this again
    sed -i -e "s/512/$n/g" parameters.h
    sed -i -e "s/256/$n/g" parameters.h
    sed -i -e "s/128/$n/g" parameters.h
    sed -i -e "s/64/$n/g" parameters.h
    sed -i -e "s/32/$n/g" parameters.h
    sed -i -e "s/16/$n/g" parameters.h

    sed -i -e "s/NMODE 0/NMODE $i/g" parameters.h
    sed -i -e "s/NMODE 1/NMODE $i/g" parameters.h
    sed -i -e "s/NMODE 2/NMODE $i/g" parameters.h
    sed -i -e "s/NMODE 3/NMODE $i/g" parameters.h

    until ./run.sh $OUT_DIR > $OUT_DIR/output_${n}_${i}.txt
    do
	echo "Retrying run"
    done

    mv $OUT_DIR/dumps $OUT_DIR/dumps_${n}_${i}
    rm -rf $OUT_DIR/restarts

  done
done
