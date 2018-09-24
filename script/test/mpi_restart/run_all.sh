#!/bin/bash

./test_restart_diffmpi.sh torus 0
./test_restart_diffmpi.sh torus 1
# More MADs?
./test_restart_diffmpi.sh bondi
./test_restart_diffmpi.sh mhdmodes

