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
from subprocess import call
from psutil import cpu_count

def matches_host():                                                              
  host = os.uname()[1]                                                           
  if host == 'lmc':                                                            
    return True                                                                  
  return False

def get_options():
  host = {}

  host['NAME']           = os.uname()[1]
  host['COMPILER']       = '/bhhome/bryan10/software/hdf5/bin/h5pcc'
  host['COMPILER_FLAGS'] = '-Ofast -Wall -fdiagnostics-color -fopenmp'
  host['GSL_DIR']        = '/bhhome/bryan10/software/gsl'
  host['EXECUTABLE']     = '/bhhome/bryan10/software/openmpi/bin/mpirun'

  return host

