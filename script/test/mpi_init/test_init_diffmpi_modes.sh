#!/bin/bash

source ../test_common.sh

PROB=mhdmodes

# Must be just a name for now
OUT_DIR=results_$PROB

# Initial clean and make of work area
BASEDIR=../../..
rm -rf build_archive param.dat harm
make -f $BASEDIR/makefile -j4 PROB=$PROB

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

# Give the system a reasonable size to limit runtime
set_compile_int N1TOT 64
set_compile_int N2TOT 64
set_compile_int N3TOT 64

for i in 0 1 2 3
do

    rm -rf $OUT_DIR/dumps $OUT_DIR/restarts $OUT_DIR/*.h5

    set_compile_int N1CPU 1
    set_compile_int N2CPU 1
    set_compile_int N3CPU 1

    set_run_int nmode $i

    make -f $BASEDIR/makefile -j4 PROB=$PROB
    echo "First run at size $n"
    ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/out_firsttime.txt
    echo "Done!"

    cd $OUT_DIR
    mv dumps/dump_00000000.h5 ./first_dump_gold.h5
    mv dumps/dump_00000005.h5 ./last_dump_gold.h5
    mv restarts/restart_00000001.h5 ./first_restart_gold.h5
    rm -rf dumps restarts
    cd ..

    sleep 1

    set_compile_int N1CPU 2
    set_compile_int N2CPU 2
    set_compile_int N3CPU 4
    export OMP_NUM_THREADS=2

    make -f $BASEDIR/makefile -j4 PROB=$PROB
    echo "Second run at size $n"
    mpirun -n 16 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/out_secondtime.txt
    echo "Done!"

    ./verify.sh $PROB

mv $OUT_DIR/verification_$PROB.txt $OUT_DIR/verification_${PROB}_$i.txt

done
