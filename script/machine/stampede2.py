################################################################################
#                                                                              #
#  MACHINE-SPECIFIC FUNCTIONS                                                  #
#                                                                              #
#    OPTIONS:                                                                  #
#      COMPILER   : PATH TO COMPILER EXECUTABLE                                #
#      GSL_DIR    : PATH TO GSL INSTALLATION                                   #
#      MPI_DIR    : PATH TO MPI INSTALLATION                                   #
#      HDF5_DIR   : PATH TO HDF5 INSTALLATION                                  #
#      EXECUTABLE : BINARY WRAPPER USED TO LAUNCH BHLIGHT                      #
#                                                                              #
#    MPI_DIR AND HDF5_DIR ARE NOT REQUIRED IF COMPILER HANDLES HEADERS AND     #
#    LIBRARIES FOR THESE DEPENDENCIES                                          #
#                                                                              #
################################################################################

import util
import sys
import os

# module load phdf5
# module load gsl

def matches_host():
  host = os.uname()[1]
  if '.stampede2.tacc.utexas.edu' in host:
    return True
  return False

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = 'h5pcc'
  host['COMPILER_FLAGS'] = '-xMIC-AVX512 -Ofast -funroll-loops -fargument-noalias -Wall -Werror -ipo -qopenmp -qopt-threads-per-core=4' #-vec-threshold0' #-qopt-report-phase=vec -qopt-report-file=vec.txt'
  host['HDF5_DIR']       = '/opt/apps/intel18/impi18_0/phdf5/1.8.16/x86_64/'
  host['GSL_DIR']        = '/opt/apps/intel18/gsl/2.3'
  host['EXECUTABLE']     = 'mpirun'

  return host

