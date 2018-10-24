# Common functions for test scripts
# Mostly changing compile-time and runtime parameters

# TODO move more common restart-test functions here

# Usage: set_compile_int <param> <value>
set_compile_int () {
  sed -i -e "s/$1 [0-9]\+/$1 $2/g" build_archive/parameters.h
}

# TODO accommodate more values/spaces/etc.
# Usage: set_run_dbl <param> <value>
set_run_dbl () {
  sed -i -e "s/$1 = [0-9]\+\\.[0-9]\+/$1 = $2/g" param.dat
}

# Usage: set_run_int <param> <value>
set_run_int () {
  sed -i -e "s/$1 = [0-9]\+/$1 = $2/g" param.dat
}

set_problem_size () {
  set_compile_int N1TOT $1
  set_compile_int N2TOT $2
  set_compile_int N3TOT $3
}

set_cpu_topo () {
  set_compile_int N1CPU $1
  set_compile_int N2CPU $2
  set_compile_int N3CPU $3
  export OMP_NUM_THREADS=$(( $(nproc --all) / ($1 * $2 * $3) ))
  #echo "Exporting $OMP_NUM_THREADS threads"
}

# This allows a separate make target if we want
make_harm_here () {
  [ -z ${HARM_BASE_DIR+x} ] && HARM_BASE_DIR=../../..
  [ -z ${HARM_MAKE_JOBS+x} ] && HARM_MAKE_JOBS=$(nproc --all)
  make -f $HARM_BASE_DIR/makefile -j$HARM_MAKE_JOBS PROB=$1
  # Use default param.dat if none is present in test dir
  if [ ! -f param.dat ] cp $HARM_BASE_DIR/prob/$1/param.dat .
}
