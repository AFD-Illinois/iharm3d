################################################################################
#                                                                              #
#  UTILITY FUNCTIONS                                                           #
#                                                                              #
################################################################################

import subprocess
import glob
import os

import signal
import multiprocessing
import psutil

import numpy as np

# TODO fns to process arguments

# Run a function in parallel with Python's multiprocessing
# 'function' must take only a number
def run_parallel(function, nmax, nthreads, debug=False):
  #original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
  pool = multiprocessing.Pool(nthreads)
  #signal.signal(signal.SIGINT, original_sigint_handler)
  try:
    pool.map_async(function, range(nmax)).get(720000)
  except KeyboardInterrupt:
    print 'Caught interrupt!'
    pool.terminate()
    exit(1)
  else:
    pool.close()
  pool.join()

# Calculate ideal # threads
def calc_nthreads(hdr, n_mkl=4):
  # Limit threads for 192^3+ problem due to memory
  if hdr['n1'] * hdr['n2'] * hdr['n3'] >= 192 * 192 * 192:
    # Try to add some parallelism MKL
    try:
      import ctypes
    
      mkl_rt = ctypes.CDLL('libmkl_rt.so')
      mkl_set_num_threads = mkl_rt.MKL_Set_Num_Threads
      mkl_get_max_threads = mkl_rt.MKL_Get_Max_Threads
      mkl_set_num_threads(n_mkl)
      print "Using", mkl_get_max_threads(), "MKL threads"
    except Error as e:
      print(e)
    
    # Roughly compute memory and leave some generous padding for multiple copies and Python games
    # TODO depend on success above?
    return int(0.11 * psutil.virtual_memory().total/(hdr['n1']*hdr['n2']*hdr['n3']*10*8))
  else:
    return psutil.cpu_count(logical=False)

# COLORIZED OUTPUT
class color:
  BOLD    = '\033[1m'
  WARNING = '\033[1;31m'
  BLUE    = '\033[94m'
  NORMAL  = '\033[0m'

def get_files(PATH, NAME):                                                       
  return np.sort(glob.glob(os.path.join(PATH,'') + NAME))

# PRINT ERROR MESSAGE
def warn(mesg):
  print(color.WARNING + "\n  ERROR: " + color.NORMAL + mesg + "\n")

# APPEND '/' TO PATH IF MISSING
def sanitize_path(path):
  return os.path.join(path, '')

# SEND OUTPUT TO LOG FILE AS WELL AS TERMINAL
def log_output(sys, logfile_name):
  import re
  f = open(logfile_name, 'w')
  class split(object):
    def __init__(self, *files):
      self.files = files
    def write(self, obj):
      n = 0
      ansi_escape = re.compile(r'\x1b[^m]*m')
      for f in self.files:
        if n > 0:
          f.write(ansi_escape.sub('', obj))
        else:
          f.write(obj)
        f.flush()
        n += 1
    def flush(self):
      for f in self.files:
        f.flush()
  sys.stdout = split(sys.stdout, f)
  sys.stderr = split(sys.stderr, f)

# CREATE DIRECTORY
def make_dir(path):
  if not os.path.exists(path):
    os.makedirs(path)

# CALL rm -rf ON RELATIVE PATHS ONLY
def safe_remove(path):
  import sys
  from subprocess import call
  
  # ONLY ALLOW RELATIVE PATHS
  if path[0] == '/':
    warn("DIRECTORY " + path + " IS NOT A RELATIVE PATH! DANGER OF DATA LOSS")
    sys.exit()
  elif os.path.exists(path):
    call(['rm', '-rf', path])

