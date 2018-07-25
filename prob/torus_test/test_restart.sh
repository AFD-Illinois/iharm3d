#!/bin/bash

# Must be just a name for now
OUT_DIR=restart_test

#TODO set NMPI properly instead of changing here
export NMPI=1

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

./run.sh $OUT_DIR > $OUT_DIR/out_firsttime.txt

cd $OUT_DIR

cp dumps/dump_00000005.h5 ./last_dump_gold.h5
cp restarts/restart_00000001.h5 .

sleep 1
rm -rf restarts
mkdir restarts
cd restarts
mv ../restart_00000001.h5 .
ln -s restart_00000001.h5 restart.last

cd ../..

./run.sh $OUT_DIR > $OUT_DIR/out_secondtime.txt

cd $OUT_DIR

# Verification
set -x

grep restart out_firsttime.txt

grep restart out_secondtime.txt

h5diff --delta=1e-10 last_dump_gold.h5 dumps/dump_00000005.h5
