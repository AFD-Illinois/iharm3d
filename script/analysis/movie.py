################################################################################
#                                                                              #
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      #
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

# Movie size in matplotlib units
FIGX = 16
FIGY = 9
# For plotting debug, "array-space" plots
# Default for all plots, can be overridden below
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

# Choose between several predefined layouts below
# For keeping around lots of possible movies with same infrastructure
movie_type = "simpler"

# ARGUMENTS
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

# LOAD FILES
files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))
gridfile = os.path.join(path,"grid.h5")

if len(files) == 0:
    util.warn("INVALID PATH TO DUMP FOLDER")
    sys.exit(1)

FRAMEDIR = 'FRAMES'
util.make_dir(FRAMEDIR)

hdr = io.load_hdr(files[0])
geom = io.load_geom(hdr, gridfile)

# Avoid passing globals in and out everywhere
bplt.init_plotting(hdr, geom)
init_analysis(hdr, geom)

# DECIDE THREADS
if hdr['n1'] * hdr['n2'] * hdr['n3'] >= 192*192*192:
  #Roughly compute memory and leave some generous padding for multiple copies and Python games
  nthreads = int(0.1 * psutil.virtual_memory().total/(hdr['n1']*hdr['n2']*hdr['n3']*10*8))
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
  # Only skip this if no bsq/beta/etc
  # Obvs add extras if plotting them
  dump = io.load_dump(files[n], hdr, geom, derived_vars = False, extras = False)
  
  fig = plt.figure(figsize=(FIGX, FIGY))
  fig.suptitle("t = %d"%dump['t']) # TODO put this at end...
  # Usual layout: 
  #ax_slc = plt.subplots(nplotsy, nplotsx)
  #ax_flux = plt.subplots(nplotsy*2, nplotsx/2)
  # Just two large slices
  gs = gridspec.GridSpec(3, 2, height_ratios=[4, 1, 1])
  #gs = gridspec.GridSpec(2, 2, height_ratios=[7, 2])
  ax_slc = [plt.subplot(gs[0,0]), plt.subplot(gs[0,1])]
  ax_flux = [plt.subplot(gs[1,:]), plt.subplot(gs[2,:])]

  if movie_type == "simplest":
    # Simplest movie: just RHO
    ax_slc = plt.subplots(1,2)
    bplt.plot_slices(ax_slc[0], ax_slc[1], r"$\log_{10}(\rho)$", np.log10(dump['RHO']), dump, -3, 2, cmap='YlOrRd')
  elif movie_type == "simpler":
    # Simpler movie: RHO and mdot
    gs = gridspec.GridSpec(2, 2, height_ratios=[5, 1])
    ax_slc = [plt.subplot(gs[0,0]), plt.subplot(gs[0,1])]
    ax_flux = [plt.subplot(gs[1,:])]
    bplt.plot_slices(ax_slc[0], ax_slc[1], r"$\log_{10}(\rho)$", np.log10(dump['RHO']), dump, -3, 2, cmap='YlOrRd')
    bplt.diag_plot(ax_flux[0], diag, dump, 'mdot', r"$\dot{M}$", logy=LOG_MDOT)
  elif movie_type == "simple":
    # Simple movie: RHO mdot phi
    gs = gridspec.GridSpec(3, 2, height_ratios=[4, 1, 1])
    ax_slc = [plt.subplot(gs[0,0]), plt.subplot(gs[0,1])]
    ax_flux = [plt.subplot(gs[1,:]), plt.subplot(gs[2,:])]
    bplt.plot_slices(ax_slc[0], ax_slc[1], r"$\log_{10}(\rho)$", np.log10(dump['RHO']), dump, -3, 2, cmap='YlOrRd')
    bplt.diag_plot(ax_flux[0], diag, dump, 'mdot', r"$\dot{M}$", logy=LOG_MDOT)
    bplt.diag_plot(ax_flux[1], diag, dump, 'phi', r"$\phi_BH$", logy=LOG_PHI)
  else: # All other movie types share a layout
    ax_slc = plt.subplots(nplotsy, nplotsx)
    ax_flux = plt.subplots(nplotsy*2, nplotsx/2)
    if movie_type == "traditional":
      # Usual movie: RHO beta fluxes
      # CUTS
      bplt.plot_slices(ax_slc[0], ax_slc[1], 'RHO', np.log10(dump['RHO']), dump, -3, 2)
      bplt.plot_slices(ax_slc[4], ax_slc[5], 'beta', np.log10(dump['beta']), dump, -2, 2, cmap="RdBu_r")
      # FLUXES
      bplt.diag_plot(ax_flux[0], diag, dump, 'mdot', 'mdot', logy=LOG_MDOT)
      bplt.diag_plot(ax_flux[1], diag, dump, 'phi', 'phi_BH', logy=LOG_PHI)
      # Mixins:
      # Zoomed in RHO
      bplt.plot_slices(ax_slc[6], ax_slc[7], 'RHO', np.log10(dump['RHO']), dump, -3, 2, 7, window=[-10,10,-10,10], overlay_field=False)
      # Bsq
      #bplt.plot_slices('bsq', np.log10(dump['bsq']), dump, -5, 0, 7)
      # Failures: all failed zones, one per nonzero pflag
      #bplt.plot_slices('fails', dump['fail'] != 0, dump, 0, 20, 7, cmap='Reds', int=True) #, arrspace=True)
      # 2D histograms
      #bplt.hist_2d(ax1, np.log10(dump['RHO']), np.log10(dump['UU']),"RHO", "UU", logcolor=True)
      #bplt.hist_2d(ax2, np.log10(dump['UU']), np.log10(dump['bsq']),"UU", "bsq", logcolor=True)
      
      # Extra fluxes:
      bplt.diag_plot(ax_flux[1], diag, dump, 'edot', 'Edot', logy=LOG_PHI)
    elif movie_type == "e_ratio":
      # Energy ratios: difficult places to integrate, with failures
      bplt.plot_slices(ax_slc[0], ax_slc[1], 'UU/RHO', np.log10(dump['UU']/dump['RHO']), dump, -3, avg=True)
      bplt.plot_slices(ax_slc[2], ax_slc[3], 'sigma', np.log10(dump['bsq']/dump['RHO']), dump, -3, 3, 3)
      bplt.plot_slices(ax_slc[4], ax_slc[5], 'inverse beta', np.log10(1/dump['beta']), dump, -3, 3, 5, avg=True)
      bplt.plot_slices(ax_slc[6], ax_slc[7], 'fails', dump['fail'] != 0, dump, 0, 20, 7, cmap='Reds', int=True) #, arrspace=True)
    elif movie_type == "conservation":
      # Continuity plots to verify local conservation of energy, angular + linear momentum
      # Integrated T01: continuity for momentum conservation
      bplt.plot_slices(ax_slc[0], ax_slc[1], 'T^1_0 integrated', Tmixed(dump, 1, 0), dump, 0, 600, 1, arrspace=True, int=True)
      # integrated T00: continuity plot for energy conservation
      bplt.plot_slices(ax_slc[4], ax_slc[5], 'T^0_0 integrated', np.abs(Tmixed(dump, 0, 0)), dump, 0, 3000, 5, arrspace=True, int=True)
      # Radial conservation plots
      E_r = sum_shell(geom,Tmixed(geom,dump,0,0))
      Ang_r = sum_shell(geom,Tmixed(geom,dump,0,3))
      mass_r = sum_shell(geom,dump['ucon'][:,:,:,0]*dump['RHO'])
      
      ax = plt.subplot(nplotsy*2,nplotsx/2,nplotsx/2*3)
      bplt.radial_plot(ax, np.abs(E_r), 'Conserved vars at R', ylim=(0,1000), col='b', rlim=(0,20), arrayspace=True)
      bplt.radial_plot(ax, np.abs(Ang_r)/10, '', ylim=(0,1000), col='r', rlim=(0,20), arrayspace=True)
      bplt.radial_plot(ax, np.abs(mass_r), '', ylim=(0,1000), rlim=(0,20), arrayspace=True)
      
      # Radial energy accretion rate
      Edot_r = sum_shell(geom,Tmixed(geom,dump,1,0))
      bplt.radial_plot(ax, np.abs(Edot_r), 'Edot at R', ylim=(0,200), rlim=(0,20), arrayspace=True)

      # Radial integrated failures
      bplt.radial_plot(ax, (dump['fail'] != 0).sum(axis=(1,2)), 'Fails at R', arrayspace=True, rlim=[0,50], ylim=[0,1000])

    elif movie_type == "ceilings":
      # TODO add measures of all ceiling efficacy.  Record ceilings in header or extras?
      bplt.plot_slices('sigma ceiling', dump['bsq']/dump['RHO'] - 100, dump, -100, 100, 5, cmap='RdBu_r')
      bplt.diag_plot(ax, diag, dump, 'sigma_max', 'sigma_max')

  # TODO enlarge plots w/o messing up even pixel count
  # Maybe share axes, even?
  plt.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95) # Avoid crowding
  #plt.tight_layout()

  plt.savefig(imname, dpi=100) #, bbox_inches='tight')
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
