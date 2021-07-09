#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

for test_problem in *
do
    # Skip non-directories, otherwise cd in
    [ ! -d $test_problem ] && continue
    cd $test_problem

    bash run_all.sh

    cd ..
done
