GSL_DIR = /opt/apps/intel18/gsl/2.3

CFLAGS = -xMIC-AVX512 -Ofast -funroll-loops -Wall -Werror -ipo -qopenmp -qopt-zmm-usage=high

# Additional arguments that have been tried
#-fargument-noalias -qopt-threads-per-core=4'
#-vec-threshold0' 
# Report vectorization
#-qopt-report-phase=vec -qopt-report-file=vec.txt

