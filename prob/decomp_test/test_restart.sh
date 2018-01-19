#!/bin/bash

# Must be just a name for now
OUT_DIR=restart_test

mkdir -p $OUT_DIR

./run.sh $OUT_DIR > $OUT_DIR/out_firsttime.txt

cd $OUT_DIR

cp restarts/restart.last ./last_restart_gold.h5
cp restarts/restart_00000001.h5 .

rm -rf restarts
mkdir restarts
mv restart_00000001.h5 restarts/

cd ..

./run.sh $OUT_DIR > $OUT_DIR/out_secondtime.txt

h5diff --delta=1e-10 last_restart_gold.h5 restarts/restart.last
