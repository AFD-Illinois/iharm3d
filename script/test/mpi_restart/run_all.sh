#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

#./test_restart_diffmpi.sh torus "0 1"
./test_restart_diffmpi.sh bondi
./test_restart_diffmpi.sh mhdmodes

