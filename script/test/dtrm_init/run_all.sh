#!/bin/bash

# Output consistency test over multiple invocations
./test_init_dtrm.sh torus 0
./test_init_dtrm.sh torus 1
# Could put more MADs
./test_init_dtrm.sh mhdmodes
./test_init_dtrm.sh bondi

# Every-dump tests for particular bug
./test_out_everystep.sh torus 0
./test_out_everystep.sh torus 1
./test_out_everystep.sh mhdmodes
./test_out_everystep.sh bondi
