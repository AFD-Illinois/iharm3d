# TODO specify consistent h5pcc/libraries/flags combinations
CFLAGS = -std=gnu99 -O3 -march=core-avx2 -mtune=core-avx2 -flto -fopenmp -funroll-loops -mcmodel=large

#CC = /opt/phdf5-intel/bin/h5pcc
#CFLAGS = -std=gnu99 -xCORE-AVX2 -O3 -Wall -Werror -qopenmp -ipo
#GSL_DIR = /opt/gsl-intel/
#MATH_LIB =
