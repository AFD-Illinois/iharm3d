GSL_DIR = /opt/apps/intel18/gsl/2.3

# For SKX
#CFLAGS = -xCORE-AVX512
# For KNL
CFLAGS = -xMIC-AVX512

CFLAGS += -Ofast -funroll-loops -Wall -Werror -ipo -qopenmp -qopt-zmm-usage=high

MATH_LIB =

# Additional arguments that have been tried
#-fargument-noalias -qopt-threads-per-core=4'
#-vec-threshold0' 
# Report vectorization
#-qopt-report-phase=vec -qopt-report-file=vec.txt

#note: mcmodel is bad and I should fix the code not to use it

#GSL_DIR = /opt/apps/gcc7_1/gsl/2.3/

#CC = h5pcc -shlib

# FOR KNL
#CFLAGS = -march=knl -mtune=knl
# FOR SKX
#CFLAGS = -march=skylake-avx512 -mtune=skylake-avx512

#CFLAGS += -std=gnu99 -O3 -flto -fopenmp -funroll-loops -mcmodel=large
