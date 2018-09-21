# Common functions for test scripts
# Mostly changing compile-time and runtime parameters

# Usage: set_compile_int <param> <value>
set_compile_int () {
  sed -i -e "s/$1 [0-9]\+/$1 $2/g" build_archive/parameters.h
}

# TODO accommodate more values/spaces/etc.
# Usage: set_run_dbl <param> <value>
set_run_dbl () {
  sed -i -e "s/$1 = [0-9]\+\\.[0-9]\+/$1 = $2/g" param.dat
}

# Usage: set_run_dbl <param> <value>
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
}