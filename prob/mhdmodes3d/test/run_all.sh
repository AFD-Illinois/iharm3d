#!/bin/bash

cd ..

./run_test.sh

cd test

python test_faux2D.py

xdg-open mhdmodes3d_ALFVEN.png
