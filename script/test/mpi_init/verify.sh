# Verification

PROB=$1

cd results_$PROB

exec > verification_$PROB.txt 2>&1
set -x

grep restart out_firsttime.txt

grep restart out_secondtime.txt

h5diff --delta=1e-12 first_dump_gold.h5 dumps/dump_00000000.h5

h5diff --delta=1e-12 first_restart_gold.h5 restarts/restart_00000001.h5

if [ $PROB == "torus" ]
then
h5diff --delta=1e-12 last_dump_gold.h5 dumps/dump_00000001.h5
h5diff --delta=1e-6 last_dump_gold.h5 dumps/dump_00000001.h5
h5diff --delta=1e-3 last_dump_gold.h5 dumps/dump_00000001.h5
else
h5diff --delta=1e-12 last_dump_gold.h5 dumps/dump_00000005.h5
h5diff --delta=1e-6 last_dump_gold.h5 dumps/dump_00000005.h5
h5diff --delta=1e-3 last_dump_gold.h5 dumps/dump_00000005.h5
fi
