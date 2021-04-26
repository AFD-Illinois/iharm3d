#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

RESULTS_DIR=$PWD/test-results

# Keep a separate folder of just results
rm -rf $RESULTS_DIR
mkdir -p $RESULTS_DIR

# In CI, errors mean something, and we do not care about copies
cpalways () { cp $1 $2 2>/dev/null ; return 0 ; }

for test_problem in *
do
    # Skip non-directories, otherwise cd in
    [ ! -d $test_problem ] && continue
    cd $test_problem

    bash run_all.sh

    mkdir -p $RESULTS_DIR/$test_problem
    cpalways -r results_*/plots $RESULTS_DIR/$test_problem
    cpalways results_*/*.txt $RESULTS_DIR/$test_problem
    cpalways results_*/*.png $RESULTS_DIR/$test_problem

    cd ..
done
