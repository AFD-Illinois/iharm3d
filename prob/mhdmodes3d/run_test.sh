#!/bin/bash

source ../../../mpi-load.sh

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

    python build.py
    if [ $? -eq 0 ]; then
      ./run.sh

      rm -rf dumps_${n}_${i}
      mv dumps dumps_${n}_${i}
    fi

  done
done
