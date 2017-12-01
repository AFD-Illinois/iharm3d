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

FIGX = 20
FIGY = 10
SIZE = 40

# TODO this is a dirty hack
args_bad = False
if sys.argv[1] == '-d':
    debug = True
    path = sys.argv[2]
    if len(sys.argv) != 3:
        args_bad = True
else:
    debug = False
    path = sys.argv[1]
    if len(sys.argv) != 2:
        args_bad = True

if args_bad:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit(1)

files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))

if len(files) == 0:
    util.warn("INVALID PATH TO DUMP FOLDER")
    sys.exit(1)

FRAMEDIR = 'FRAMES'
util.make_dir(FRAMEDIR)

hdr = io.load_hdr(files[0])
geom = io.load_geom(hdr, files[0])

diag = io.load_log(os.path.join(path, "log.out"))

def plot(args):
  n = args
  imname = 'frame_%08d.png' % n
  imname = os.path.join(FRAMEDIR, imname)
  print '%08d / ' % (n+1) + '%08d' % len(files) 
 
  # Ignore if frame already exists
  if os.path.isfile(imname):
    return

  dump = io.load_dump(hdr, geom, files[n])
  fig = plt.figure(figsize=(FIGX, FIGY))
  
  ax = plt.subplot(2,3,1)
  bplt.plot_xz(ax, geom, np.log10(dump['RHO']), dump, 
    vmin=-4, vmax = 0, label='RHO')
  bplt.overlay_field(ax, geom, dump)
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

  # I change this a lot
  var2_str = 'bsq'
  var2_data = np.log10(dump[var2_str])
  
  ax = plt.subplot(2,3,2)
  bplt.plot_xz(ax, geom, var2_data, dump,
    label=var2_str, cmap='RdBu_r', vmin=-8, vmax=2)
  bplt.overlay_field(ax, geom, dump)
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
 
  ax = plt.subplot(2,3,4)
  bplt.plot_xy(ax, geom, np.log10(dump['RHO']), dump,
     vmin=-4, vmax=0, label='RHO')
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  
  ax = plt.subplot(2,3,5)
  bplt.plot_xy(ax, geom, var2_data, dump,
     label=var2_str, cmap='RdBu_r', vmin=-8, vmax=2)
  ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])
  
  ax = plt.subplot(2,3,3)
  ax.plot(diag['t'], diag['mdot'], color='k')
  ax.axvline(dump['t'], color='r') # Trace on finished plot, don't draw it as we go
  ax.set_xlim([0, dump['tf']])
  ax.set_xlabel('t/M')
  ax.set_ylabel('mdot')
  
  # This fucks with video because it generates an odd image size
  #plt.subplots_adjust(hspace=0.50) # Avoid crowding

  #ax.pcolormesh(dump['X1'][:,:,0], dump['X2'][:,:,0], dump['RHO'][:,:,0])
  plt.savefig(imname, bbox_inches='tight', dpi=100)
  plt.close(fig)

# Test plot so that backtraces work
# And for quicker turnaround
if debug:
    plot(0)
    plot(100)
    exit(0)

import multiprocessing
import signal
import psutil

nthreads = psutil.cpu_count(logical=False)
print 'Number of CPUs: %i' % psutil.cpu_count(logical=False)

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(nthreads)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  res = pool.map_async(plot, range(len(files)))
  res.get(720000)
except KeyboardInterrupt:
  print 'Caught interrupt!'
  pool.terminate()
  exit(1)
else:
  pool.close()
pool.join()

