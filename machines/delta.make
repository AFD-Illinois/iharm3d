# Default variables for running IMEX problem with gcc. Modules: 
# 1) cue-login-env/1.0   4) intel-oneapi-mkl/2022.0.2   7) rdma-core/32.0  10) libevent/2.1.8       13) openssh/8.0p1  16) hdf5/1.13.1
# 2) default             5) gcc/11.2.0                  8) ucx/1.11.2      11) libfabric/1.14.0     14) pmix/3.2.3     17) anaconda3_cpu/4.13.0
# 3) modtree/cpu         6) openblas/0.3.20             9) knem/1.1.4      12) lustre/2.14.0_ddn23  15) openmpi/4.1.2  18) gsl/2.7

CFLAGS = -march=native -mtune=native -std=gnu11 -O3 -flto -fopenmp -funroll-loops
CFLAGS += -DMKL_ILP64  -m64  -I"${MKLROOT}/include"
MKL_DIR = /sw/spack/deltacpu-2022-03/apps/intel-oneapi-mkl/2022.0.2-gcc-8.4.1-pes6zb6/mkl/2022.0.2/
MKL_LIB = -lmkl_scalapack_ilp64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_ilp64 -lgomp -lpthread -lm -ldl
HDF5_DIR = /sw/spack/deltacpu-2022-03/apps/hdf5/1.13.1-gcc-11.2.0-dschtbv

#HDF5_DIR = /sw/spack/deltacpu-2022-03/apps/hdf5/1.13.1-aocc-3.2.0-o3l6ma2
#GSL_DIR = /u/vdhruv2/gsl
