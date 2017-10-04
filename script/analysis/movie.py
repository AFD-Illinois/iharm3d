################################################################################
#                                                                              # 
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      # 
#                                                                              # 
################################################################################

import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../')
import numpy as np
import hdf5_to_dict as io
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import util
import glob
import os
import plot as bplt

FIGX = 13
FIGY = 10
SIZE = 40

if len(sys.argv) != 2:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit()

path = sys.argv[1]

files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))

FRAMEDIR = 'FRAMES'
util.make_dir(FRAMEDIR)

def plot(args):
  n = args
  print '%08d / ' % (n+1) + '%08d' % len(files) 
  dump = io.load_dump(files[n])
  fig = plt.figure(figsize=(FIGX, FIGY))
  
  ax = plt.subplot(2,2,1)
  bplt.plot_xz(ax, np.log10(dump['RHO']), dump, 
    vmin=-4, vmax = 0, label='RHO')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  
  ax = plt.subplot(2,2,2)
  bplt.plot_xz(ax, np.log10(dump['Thetae']), dump,
    vmin=-2, vmax=2, label='Thetae', cmap='RdBu_r')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
 
  ax = plt.subplot(2,2,3)
  bplt.plot_xy(ax, np.log10(dump['RHO']), dump,
     vmin=-4, vmax=0, label='RHO')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  
  ax = plt.subplot(2,2,4)
  bplt.plot_xy(ax, np.log10(dump['Thetae']), dump,
     vmin=-2, vmax=2, label='Thetae', cmap='RdBu_r')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

  #ax.pcolormesh(dump['X1'][:,:,0], dump['X2'][:,:,0], dump['RHO'][:,:,0])
  plt.savefig(os.path.join(FRAMEDIR, 'frame_%08d.png' % n), 
      bbox_inches='tight', dpi=100)
  plt.close(fig)

import multiprocessing
import signal
import psutil

#nthreads = psutil.cpu_count(logical=False)
print psutil.cpu_count(logical=False)
nthreads = 10

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(nthreads)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  res = pool.map_async(plot, range(len(files)))
  res.get(720000)
except KeyboardInterrupt:
  print 'Caught interrupt!'
  pool.terminate()
else:
  pool.close()
pool.join()

