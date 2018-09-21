#!/bin/bash

source ../test_common.sh

PROB=torus

# Must be just a name for now
OUT_DIR=results_$PROB

# Initial clean and make of work area
BASEDIR=../../..
rm -rf build_archive param.dat harm
make -f $BASEDIR/makefile -j4 PROB=$PROB

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

# Give the system a reasonable size to limit runtime
set_compile_int N1TOT 96
set_compile_int N2TOT 32
if [ "$PROB" == "bondi" ]; then
  set_compile_int N3TOT 1
else
  set_compile_int N3TOT 32
fi

# Give a relatively short endpoint
# We're testing init and basic propagation
set_run_dbl tf 1.0
set_run_dbl DTd 1.0
set_run_dbl u_jitter 0.0

MADS="0 1 2 3 4"

for i in $MADS
do

rm -rf $OUT_DIR/dumps $OUT_DIR/restarts $OUT_DIR/*.h5

set_compile_int N1CPU 1
set_compile_int N2CPU 1
set_compile_int N3CPU 1

set_run_int mad_type $i

make -f $BASEDIR/makefile -j4 PROB=$PROB
./harm -p param.dat -o $OUT_DIR > $OUT_DIR/out_firsttime.txt

cd $OUT_DIR
mv dumps/dump_00000000.h5 ./first_dump_gold.h5
mv dumps/dump_00000001.h5 ./last_dump_gold.h5
mv restarts/restart_00000001.h5 ./first_restart_gold.h5
rm -rf dumps restarts
cd ..

sleep 1

set_compile_int N1CPU 2
set_compile_int N2CPU 2
set_compile_int N3CPU 4
export OMP_NUM_THREADS=2

make -f $BASEDIR/makefile -j4 PROB=$PROB
mpirun -n 16 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/out_secondtime.txt

./verify.sh $PROB

mv $OUT_DIR/verification_torus.txt $OUT_DIR/verification_torus_$i.txt

done
