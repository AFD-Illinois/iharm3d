# Verification

PROB=$1

cd results_$PROB

if [ $PROB == "mhdmodes" ]
then
  LAST_DUMP=dumps/dump_00000005.h5
else
  LAST_DUMP=dumps/dump_00000001.h5
fi

cp ../../../analysis/*.py .
python plot_diff.py first_dump_gold.h5 $LAST_DUMP dumps/grid.h5 differences_$PROB.png

exec > verification_$PROB.txt 2>&1
set -x

grep restart out_firsttime.txt
grep restart out_secondtime.txt

# Diff first dumps for a sanity check
h5diff --delta=1e-12 first_dump_gold.h5 dumps/dump_00000000.h5
h5diff --delta=1e-12 first_restart_gold.h5 restarts/restart_00000001.h5

# Diff last dumps
h5diff --delta=1e-12 last_dump_gold.h5 $LAST_DUMP
h5diff --delta=1e-8 last_dump_gold.h5 $LAST_DUMP
h5diff --delta=1e-4 last_dump_gold.h5 $LAST_DUMP
