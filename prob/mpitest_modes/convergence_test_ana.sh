#!/bin/bash

OUT_DIR=test_convergence

# Analysis and plots
cp -r test $OUT_DIR/
cd $OUT_DIR/test
python test_3D.py
