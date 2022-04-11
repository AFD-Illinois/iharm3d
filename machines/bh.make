CC = /data/bh-bd1/phdf5-oneapi/bin/h5pcc
# Intel compiler flags. Use . /opt/intel/oneapi/setvars.sh to load it.
CFLAGS = -shlib -Ofast -fstrict-aliasing -Wall -Werror -ipo -qopenmp -qmkl
