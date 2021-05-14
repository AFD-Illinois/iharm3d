#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Output consistency test over multiple invocations
./test_init_dtrm.sh torus
./test_init_dtrm.sh mhdmodes
./test_init_dtrm.sh bondi
