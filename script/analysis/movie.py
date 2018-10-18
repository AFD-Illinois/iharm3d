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
geom = io.load_geom(hdr, gridfile)

if hdr['n1'] * hdr['n2'] * hdr['n3'] >= 192*192*192:
  #Roughly compute memory and leave some generous padding for multiple copies and Python games
  nthreads = int(0.1 * psutil.virtual_memory().total/(hdr['n1']*hdr['n2']*hdr['n3']*10*8))
else:
  nthreads = psutil.cpu_count(logical=False)

print 'Number of threads: %i' % nthreads

nplotsx = 4
nplotsy = 2

# Load diagnostics from HARM itself
#diag = io.load_log(hdr, os.path.join(path, "log.out"))

# Or from  Python post-analysis
# TODO make this a runtime option
diag = pickle.load(open("eht_out.p", 'rb'))

def plot_slices(name, data, dump, min, max, subplot, avg=False, int=False, window=[-SIZE,SIZE,-SIZE,SIZE], arrspace=USEARRSPACE,
                overlay_field=True, cmap='jet', xlabel=True, ylabel=True):
  # Switch to RdBu_r default?
  if int:
    # Average multiplied values for the integral
    data_xz = data*data.shape[2] #N1,2,(3)
    data_xy = data*data.shape[1] #N1,(2),3
    avg = True
  else:
    data_xz = data
    data_xy = data
  
  ax = plt.subplot(nplotsy, nplotsx, subplot)
  bplt.plot_xz(ax, geom, data_xz, window=window, cbar=False, cmap=cmap, xlabel=xlabel, ylabel=ylabel,
               label=name, vmin=min, vmax=max, arrayspace=arrspace, average=avg)
  if overlay_field and not USEARRSPACE:
    bplt.overlay_field(ax, geom, dump, NLINES)

  ax = plt.subplot(nplotsy, nplotsx, subplot+1)
  bplt.plot_xy(ax, geom, data_xy, window=window, cmap=cmap, xlabel=xlabel, ylabel=ylabel,
               label=name, vmin=min, vmax=max, arrayspace=arrspace, average=avg)


def plot(n):
  imname = 'frame_%08d.png' % n
  imname = os.path.join(FRAMEDIR, imname)
  print '%08d / ' % (n+1) + '%08d' % len(files)

  # Don't calculate b/ucon/cov/e- stuff unless we need it below
  # Only skip this if no bsq/beta/etc
  dump = io.load_dump(files[n], geom, hdr, derived_vars = True)
  
  fig = plt.figure(figsize=(FIGX, FIGY))

  # Subplots 1 & 2
  plot_slices('RHO', np.log10(dump['RHO']), dump, -3, 2, 1)
  #plot_slices('beta', np.log10(dump['beta']), dump, -2, 2, 1)
  #plot_slices('UU/RHO', np.log10(dump['UU']/dump['RHO']), dump, -3, 3, 1, avg=True)

  # Subplots 3 & 4
  #plot_slices('sigma', np.log10(dump['bsq']/dump['RHO']), dump, -3, 3, 3, avg=True)

  # Subplots 5 & 6
  #plot_slices('inverse beta', np.log10(1/dump['beta']), dump, -3, 3, 5, avg=True)
  #plot_slices('magnetization', dump['bsq']/dump['RHO'], dump, 0, 1000, 5)
  plot_slices('beta', np.log10(dump['beta']), dump, -2, 2, 5)

  # Subplots 7 & 8
  # Zoomed in RHO
  plot_slices('RHO', np.log10(dump['RHO']), dump, -3, 2, 7, window=[-SIZE/4,SIZE/4,-SIZE/4,SIZE/4], overlay_field=False)
  # Bsq
  #plot_slices('bsq', np.log10(dump['bsq']), dump, -5, 0, 7)
  # Failures: all failed zones, one per nonzero pflag
  #plot_slices('fails', dump['fail'] != 0, dump, 0, 20, 7, cmap='Reds', int=True) #, arrspace=True)

  plt.subplots_adjust(hspace=0.15, wspace=0.4)

  # Fluxes as top right corner pair of frames
  # Don't plot time-based variables for initial conditions
  if len(diag['t'].shape) > 0:
    ax = plt.subplot(nplotsy*2,nplotsx/2,nplotsx/2)
    bplt.diag_plot(ax, diag, dump, 'mdot', 'mdot', logy=LOG_MDOT)
 
    ax = plt.subplot(nplotsy*2,nplotsx/2,nplotsx)
    bplt.diag_plot(ax, diag, dump, 'phi', 'phi_BH', logy=LOG_PHI)
    
    # Alternative to 7 & 8: more fluxes
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
# TODO if name = main

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
