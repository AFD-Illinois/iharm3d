################################################################################
#                                                                              #
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('Agg')

import sys; sys.dont_write_bytecode = True
import numpy as np
import hdf5_to_dict as io
import matplotlib.pyplot as plt
import util
import glob
import os
import plot as bplt
import pickle

import multiprocessing
import signal
import psutil

FIGX = 20
FIGY = 10
SIZE = 40
NLINES = 10

# For plotting debug, "array-space" plots
USEARRSPACE = False

MAD = True
if MAD:
  LOG_MDOT = False
  MAX_MDOT = 80
  LOG_PHI = False
  MAX_PHI = 100
else:
  # SANE
  LOG_MDOT = False
  MAX_MDOT = 1
  LOG_PHI = False
  MAX_PHI = 5
  


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

# Override above option
#debug = True

if args_bad:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit(1)

files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))
gridfile = os.path.join(path,"grid.h5")

if len(files) == 0:
    util.warn("INVALID PATH TO DUMP FOLDER")
    sys.exit(1)

FRAMEDIR = 'FRAMES'
util.make_dir(FRAMEDIR)

hdr = io.load_hdr(files[0])
geom = io.load_geom(gridfile)

if hdr['n1'] >= 256 or hdr['n2'] >= 256 or hdr['n3'] >= 256:
  #Roughly compute memory and leave some generous padding for multiple copies and Python games
  nthreads = int(0.1 * psutil.virtual_memory().total/(hdr['n1']*hdr['n2']*hdr['n3']*10*8))
else:
  nthreads = psutil.cpu_count(logical=False)

print 'Number of threads: %i' % nthreads


# Load diagnostics from HARM itself
#diag = io.load_log(hdr, os.path.join(path, "log.out"))

# Or from  Python post-analysis
# TODO make this an option
diag = pickle.load(open("eht_out.p", 'rb'))

def plot(args):
  n = args
  imname = 'frame_%08d.png' % n
  imname = os.path.join(FRAMEDIR, imname)
  print '%08d / ' % (n+1) + '%08d' % len(files)

  # Ignore if frame already exists
  if os.path.isfile(imname):
    return

  dump = io.load_dump(files[n], geom, hdr, calc_omega=False)
  fig = plt.figure(figsize=(FIGX, FIGY))

  ax = plt.subplot(2,4,1)
  bplt.plot_xz(ax, geom, np.log10(dump['RHO']), dump,
               vmin=-3, vmax=2, label='RHO', arrayspace=USEARRSPACE)
  if (USEARRSPACE):
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    bplt.overlay_field(ax, geom, dump, NLINES)
    ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

  ax = plt.subplot(2,4,2)
  bplt.plot_xy(ax, geom, np.log10(dump['RHO']), dump,
               vmin=-3, vmax=2, label='RHO', arrayspace=USEARRSPACE)
  if (USEARRSPACE):
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])


  # I change this a lot
  #var2_str = 'ucon_r'
  #var2_data = dump['RHO'][:,:,:]*dump['ucon'][:,:,:,1]*dump['gdet'][:,:,None]*hdr['dx2']*hdr['dx3']
  var2_str = 'beta'
  var2_data = np.log10(dump[var2_str])
  var2_min = -2
  var2_max = 2

#  var2_str = 'magnetization'
#  var2_data = dump['bsq']/dump['RHO']
#  var2_min = 0
#  var2_max = 1000
#  var2_str = 'omega'
#  var2_data = np.log10(dump['omega']*4/hdr['a'])
#  var2_min = -2
#  var2_max = 2


  ax = plt.subplot(2,4,5)
  bplt.plot_xz(ax, geom, var2_data, dump, #cmap='RdBu_r'
               label=var2_str, vmin=var2_min, vmax=var2_max, arrayspace=USEARRSPACE)
  if (USEARRSPACE):
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    bplt.overlay_field(ax, geom, dump, NLINES)
    ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

  ax = plt.subplot(2,4,6)
  bplt.plot_xy(ax, geom, var2_data, dump, #cmap='RdBu_r'
               label=var2_str, vmin=var2_min, vmax=var2_max, arrayspace=USEARRSPACE)
  if (USEARRSPACE):
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

#  ax = plt.subplot(2,4,7)
#  bplt.plot_xz(ax, geom, var3_data, dump,
#               label=var3_str, cmap='RdBu_r', vmin=var3_min, vmax=var3_max, arrayspace=USEARRSPACE)
#  if (USEARRSPACE):
#    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
#  else:
#    bplt.overlay_field(ax, geom, dump, NLINES)
#    ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

#  ax = plt.subplot(2,4,8)
#  bplt.plot_xy(ax, geom, var3_data, dump,
#               label=var3_str, cmap='RdBu_r', vmin=var3_min, vmax=var3_max, arrayspace=USEARRSPACE)
#  if (USEARRSPACE):
#    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
#  else:
#    ax.set_xlim([-SIZE, SIZE]); ax.set_ylim([-SIZE, SIZE])

  plt.subplots_adjust(hspace=0.15, wspace=0.4)

  # Don't plot time-based variables for initial conditions
  if len(diag['t'].shape) > 0:
    ax = plt.subplot(4,2,2)
    bplt.diag_plot(ax, diag, dump, 'mdot', 'mdot', logy=LOG_MDOT)

    ax = plt.subplot(4,2,4)
    bplt.diag_plot(ax, diag, dump, 'phi', 'phi_BH', logy=LOG_PHI)
    
    #ax = plt.subplot(4,2,6)
    #bplt.diag_plot(ax, diag, dump, 'sigma_max', 'Max magnetization', ylim=None, logy=False)
    
    #ax = plt.subplot(4,2,8)
    #bplt.diag_plot(ax, diag, dump, 'egas', 'Gas energy', ylim=None, logy=False)

  # TODO enlarge plots w/o messing up even pixel count
  # Maybe share axes, even?
  plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95) # Avoid crowding
  #plt.tight_layout()

  plt.savefig(imname, dpi=100) #, bbox_inches='tight')
  plt.close(fig)

# BEGIN SCRIPT

# Test-run a couple plots directly so that backtraces work
if debug:
    plot(0)
    plot(100)
    exit(0)

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(nthreads)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  pool.map_async(plot, range(len(files))).get(720000)
except KeyboardInterrupt:
  print 'Caught interrupt!'
  pool.terminate()
  exit(1)
else:
  pool.close()
pool.join()
