#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

./test_init_diffmpi.sh torus "0 1"
./test_init_diffmpi.sh bondi
./test_init_diffmpi.sh mhdmodes

