################################################################################
#                                                                              #
#  GENERATE MOVIES COMPARING 2 SIMULATIONS' OUTPUT                             #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('Agg')

import plot as bplt
from analysis_fns import *

import sys; sys.dont_write_bytecode = True
import numpy as np
import hdf5_to_dict as io
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import util
import glob
import os
import pickle

import multiprocessing
import signal
import psutil

# Movie size in ~inches. Keep 16/9 for standard-size movies
FIGX = 12
FIGY = 27/4
# For plotting debug, "array-space" plots
# Default for all plots, can be overridden below
USEARRSPACE = False

MAD = False
LOG_MDOT = False
LOG_PHI = False

# Choose between several predefined layouts below
# For keeping around lots of possible movies with same infrastructure
movie_type = "simplest"

# ARGUMENTS
path1 = sys.argv[1]
path2 = sys.argv[2]

# LOAD FILES
files1 = np.sort(glob.glob(os.path.join(path1, "dump*.h5")))
files2 = np.sort(glob.glob(os.path.join(path2, "dump*.h5")))

if len(files1) == 0 or len(files2) == 0:
    util.warn("INVALID PATH TO DUMP FOLDER")
    sys.exit(1)

FRAMEDIR = 'FRAMES'
util.make_dir(FRAMEDIR)

# TODO TODO support different shapes/geoms
hdr = io.load_hdr(files1[0])
gridfile = os.path.join(path1,"grid.h5")
geom = io.load_geom(hdr, gridfile)

bplt.init_plotting(hdr, geom)
init_analysis(hdr, geom)

# DECIDE THREADS
if hdr['n1'] * hdr['n2'] * hdr['n3'] >= 192*192*192:
  #Roughly compute memory and leave some generous padding for multiple copies and Python games
  nthreads = int(0.05 * psutil.virtual_memory().total/(hdr['n1']*hdr['n2']*hdr['n3']*10*8))
else:
  nthreads = psutil.cpu_count(logical=False)

print 'Number of threads: %i' % nthreads

diag_post = True
if diag_post:
  # Load fluxes from post-analysis: more flexible
  diag = pickle.load(open("eht_out.p", 'rb'))
else:
  # Load diagnostics from HARM itself
  diag = io.load_log(hdr, os.path.join(path, "log.out"))

def plot(n):
  imname = 'frame_%08d.png' % n
  imname = os.path.join(FRAMEDIR, imname)
  print '%08d / ' % (n+1) + '%08d' % len(files)

  # Don't calculate b/ucon/cov/e- stuff unless we need it below
  dump1 = io.load_dump(files1[n], hdr, geom, derived_vars = False, extras = False)
  dump2 = io.load_dump(files2[n], hdr, geom, derived_vars = False, extras = False) 

  fig = plt.figure(figsize=(FIGX, FIGY))

  # Zoom in for SANEs
  if MAD:
    window = [-40,40,-40,40]
    nlines = 20
    rho_l, rho_h = -3, 2
  else:
    window = [-20,20,-20,20]
    nlines = 10
    rho_l, rho_h = -3, 0

  if movie_type == "simplest":
    # Simplest movie: just RHO
    gs = gridspec.GridSpec(1, 2, width_ratios=[16,17])
    ax_slc = [plt.subplot(gs[0]), plt.subplot(gs[1])]
    bplt.plot_xz(ax_slc[0], np.log10(dump1['RHO']), label=r"$\log_{10}(\rho)$, MAD", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
    bplt.plot_xz(ax_slc[1], np.log10(dump2['RHO']), label=r"$\log_{10}(\rho)$, SANE", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
  elif movie_type == "simpler":
    # Simpler movie: RHO and phi
    gs = gridspec.GridSpec(2, 3, height_ratios=[6, 1, 1], width_ratios=[16,17])
    ax_slc = [plt.subplot(gs[0,0]), plt.subplot(gs[0,1])]
    ax_flux = [plt.subplot(gs[1,:]), plt.subplot(gs[2,:])]
    bplt.plot_xz(ax_slc[0], np.log10(dump1['RHO']), label=r"$\log_{10}(\rho)$, MAD", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
    bplt.plot_xz(ax_slc[1], np.log10(dump2['RHO']), label=r"$\log_{10}(\rho)$, SANE", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
    bplt.diag_plot(ax_flux[0], diag, dump, 'phi', r"$\phi_{BH}$ MAD", logy=LOG_PHI, xlabel=False)
    bplt.diag_plot(ax_flux[1], diag, dump, 'phi', r"$\phi_{BH}$ SANE", logy=LOG_PHI, xlabel=False)

  pad = 0.05
  plt.subplots_adjust(left=2*pad, right=1-2*pad, bottom=pad, top=1-pad) # Avoid crowding
  #plt.tight_layout()

  plt.savefig(imname, dpi=1920/FIGX) #, bbox_inches='tight')
  plt.close(fig)

# Test-run a couple plots directly so that backtraces work
if debug:
    plot(0)
    plot(100)
    exit(0)

#original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(nthreads)
#signal.signal(signal.SIGINT, original_sigint_handler)
try:
  pool.map_async(plot, range(len(files))).get(720000)
except KeyboardInterrupt:
  print 'Caught interrupt!'
  pool.terminate()
  exit(1)
else:
  pool.close()
pool.join()
