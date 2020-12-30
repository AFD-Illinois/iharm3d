HDF5_DIR = /home/bprather/libs/hdf5-oneapi
GSL_DIR = /home/bprather/libs/gsl-oneapi

CC=/home/bprather/libs/hdf5-oneapi/bin/h5pcc

CFLAGS = -xCORE-AVX2 -Ofast -fstrict-aliasing -Wall -Werror -ipo -qopenmp
MATH_LIB =
