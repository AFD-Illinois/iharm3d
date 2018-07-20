#!/bin/bash

OUT_DIR=test_convergence

# Analysis and plots
mkdir $OUT_DIR/plots
cp plot_convergence.py $OUT_DIR/plots
cd $OUT_DIR/plots
python plot_convergence.py

