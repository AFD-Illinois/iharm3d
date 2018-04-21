# TODO specify consistent h5pcc/libraries/flags combinations
CFLAGS = -std=gnu99 -O3 -march=core-avx2 -mtune=core-avx2 -flto -fopenmp -funroll-loops -Wall -Werror

#CC = /opt/phdf5-intel/bin/h5pcc
#CFLAGS = -xCORE-AVX2 -Ofast -fstrict-aliasing -Wall -Werror -ipo -qopenmp
#GSL_DIR = /opt/gsl-intel/
