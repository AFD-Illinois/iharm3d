#!/bin/bash

source ../test_common.sh

PROB=${1:-torus}
SIZE1=96
SIZE2=64
SIZE3=32
TF=1.0

# Must be just a name for now
OUT_DIR=results_$PROB

# Initial clean and make of work area
BASEDIR=../../..
rm -rf build_archive param.dat harm
make -f $BASEDIR/makefile -j4 PROB=$PROB

# Cleanup any previous invocation or convergence test
set_compile_int N1CPU 1
set_compile_int N2CPU 1
set_compile_int N3CPU 1

# Disable electrons for fluid test
set_compile_int ELECTRONS 0

rm -rf $OUT_DIR
mkdir -p $OUT_DIR

# Set the size.  Bondi flow is 2D.
set_compile_int N1TOT $SIZE1
set_compile_int N2TOT $SIZE2
if [ "$PROB" == "bondi" ]; then
  set_compile_int N3TOT 1
else
  set_compile_int N3TOT $SIZE3
fi

set_run_dbl tf $TF
set_run_dbl DTd $TF
# TODO test setting DTr real high
# and not worrying about games below

make -f $BASEDIR/makefile -j4 PROB=$PROB
./harm -p param.dat -o $OUT_DIR > $OUT_DIR/out_firsttime.txt

cd $OUT_DIR

cp dumps/dump_00000001.h5 ./last_dump_gold.h5
cp restarts/restart_00000001.h5 .

sleep 1
# Keep the restart dir and grid/log file
rm -rf restarts/* dumps/dump_*
cd restarts
mv ../restart_00000001.h5 .
ln -s restart_00000001.h5 restart.last

cd ../..

set_compile_int N1CPU 2
set_compile_int N2CPU 2
set_compile_int N3CPU 4
export OMP_NUM_THREADS=2

make -f $BASEDIR/makefile -j4 PROB=$PROB
mpirun -n 16 ./harm -p param.dat -o $OUT_DIR > $OUT_DIR/out_secondtime.txt

./verify.sh $PROB
