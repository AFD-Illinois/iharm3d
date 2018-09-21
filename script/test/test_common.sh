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
