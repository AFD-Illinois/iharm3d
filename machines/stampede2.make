GSL_DIR = /opt/apps/intel19/gsl/2.6
CFLAGS = -xMIC-AVX512 -std=gnu11 -O3 -funroll-loops -ipo -qopenmp -qopt-prefetch=5 -mkl
#GSL_DIR = /opt/apps/gcc9_1/gsl/2.6/
#CFLAGS = -march=knl -mtune=knl -std=gnu11 -O3 -flto -fopenmp -funroll-loops -mkl
