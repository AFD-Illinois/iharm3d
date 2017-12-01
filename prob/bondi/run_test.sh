#!/bin/bash

source ../../../mpi-load.sh

for n in 32 64 128
do

  # I can't believe I'm doing this again
  sed -i -e "s/512/$n/g" parameters.h
  sed -i -e "s/256/$n/g" parameters.h
  sed -i -e "s/128/$n/g" parameters.h
  sed -i -e "s/64/$n/g" parameters.h
  sed -i -e "s/32/$n/g" parameters.h
  sed -i -e "s/16/$n/g" parameters.h

  python build.py
  if [ $? -eq 0 ]; then
    ./run.sh

    rm -rf dumps_${n}
    mv dumps dumps_${n}
  fi

done
