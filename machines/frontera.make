#GSL_DIR = /opt/apps/intel19/gsl/2.5
#CFLAGS = -xCORE-AVX512 -std=gnu11 -O3 -funroll-loops -ipo -qopenmp -qopt-prefetch=5
GSL_DIR = /opt/apps/gcc9_1/gsl/2.6/
CFLAGS = -march=native -mtune=native -std=gnu11 -O3 -flto -fopenmp -funroll-loops -mkl
