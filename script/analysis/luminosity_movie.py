################################################################################
#                                                                              #
#  MOVIE OF LUMINOSITY COMPARISON                                              #
#                                                                              #
################################################################################

import os, sys
import pickle
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import util
import hdf5_to_dict as io
import plot as bplt
from analysis_fns import *

run_name = sys.argv[1]

if "SANE" in run_name:
  SIZE = 50
  AT_R = 40
else:
  SIZE = 400
  AT_R = 100

window=[0,SIZE/2,0,SIZE]
FIGX = 15
FIGY = 15

plotfile = os.path.join("/work/03002/bprather/stampede2/movies",run_name,"eht_out.p")
avg = pickle.load(open(plotfile, "rb"))

path = os.path.join("/scratch/03002/bprather/pharm_dumps/M87SimulationLibrary/GRMHD",run_name,"dumps")
dumps = io.get_dumps_list(path)
hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr, path)

os.mkdir("frames_lmovie"+run_name.replace("/","_"))

def frame(n):
  dump = io.load_dump(dumps[n], hdr, geom)
  print("Frame {} of {}".format(n,len(dumps)))
  
  for var in ['sigma', 'Be_nob', 'betagamma', 'FE_EM', 'FE']:
    dump[var] = d_fns[var](dump)
  dump['gamma'] = get_gamma(geom, dump)
  
  fig = plt.figure(figsize=(FIGX, FIGY))
  gs = gridspec.GridSpec(2, 2, width_ratios=[1,2])
  
  ax = plt.subplot(gs[0,0])
  bplt.plot_xz(ax, geom, np.log10(dump['FE_EM']), average=True, window=window)
  ax.set_title(r"$\log_{10}( -{{T_{EM}}^r}_t )$")
  
  bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')
  
  bplt.overlay_contours(ax, geom, dump['sigma'], [1.0], color='C2')
  bplt.overlay_contours(ax, geom, dump['Be_nob'], [0.02], color='C3')
  bplt.overlay_contours(ax, geom, dump['Be_nob'], [1.0], color='C4')
  bplt.overlay_contours(ax, geom, dump['betagamma'], [1.0], color='C5')
  
  ax = plt.subplot(gs[1,0])
  bplt.plot_xz(ax, geom, np.log10(dump['FE']), average=True, window=window)
  ax.set_title(r"$\log_{10}( -{T^r}_t - \rho u^r )$")
  
  bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')
  
  bplt.overlay_contours(ax, geom, dump['sigma'], [1.0], color='C2')
  bplt.overlay_contours(ax, geom, dump['Be_nob'], [0.02], color='C3')
  bplt.overlay_contours(ax, geom, dump['Be_nob'], [1.0], color='C4')
  bplt.overlay_contours(ax, geom, dump['betagamma'], [1.0], color='C5')
  
  ax = plt.subplot(gs[0,1])
  ax.plot(avg['r'], avg['LBZ_sigma1_rt'][n], label=r"$L_{BZ}$ (sigma > 1 cut)", color='C2')
  ax.plot(avg['r'], avg['LBZ_Be_nob0_rt'][n], label=r"$L_{BZ}$ ($Be > 0.02$ cut)", color='C3')
  ax.plot(avg['r'], avg['LBZ_Be_nob1_rt'][n], label=r"$L_{BZ}$ ($Be > 1.0$ cut)", color='C4')
  ax.plot(avg['r'], avg['LBZ_bg1_rt'][n], label=r"$L_{BZ}$ ($\beta\gamma > 1.0$ cut)", color='C5')
  ax.plot(avg['r'], avg['LBZ_allp_rt'][n], label=r"$L_{BZ,tot}$", color='C6')
  
  ax.set_title(r"$L_{BZ} = \int -{{T_{EM}}^r}_t \sqrt{-g} dx^{\theta} dx^{\phi}$")
  ax.set_xlim([0,SIZE])
  ax.set_xlabel("$r$ (M)")
  ax.axvline(AT_R, color='k')
  
  if "SANE" in run_name:
    ax.set_ylim([1e-4,1e-1])
    ax.set_yscale('log')
  else:
    ax.set_ylim([1,100])
    #ax.set_yscale('log')
  
  ax.legend(loc='upper right')
  
  ax = plt.subplot(gs[1,1])
  ax.plot(avg['r'], avg['Lj_sigma1_rt'][n], label=r"$L_{jet}$ (sigma > 1 cut)", color='C2')
  #ax.plot(avg['r'], avg['Lj_Be_nob0_rt'][n], label=r"$L_{jet}$ ($Be > 0.02$ cut)", color='C3')
  ax.plot(avg['r'], avg['Lj_Be_nob1_rt'][n], label=r"$L_{jet}$ ($Be > 1.0$ cut)", color='C4')
  # TODO Gamma?
  ax.plot(avg['r'], avg['Lj_bg1_rt'][n], label=r"$L_{jet}$ ($\beta\gamma > 1.0$ cut)", color='C5')
  ax.plot(avg['r'], avg['Lj_allp_rt'][n], label=r"$L_{tot}$", color='C6')
  
  ax.set_title(r"$L_{tot} = \int (-{T^r}_t - \rho u^r) \sqrt{-g} dx^{\theta} dx^{\phi}$")
  ax.set_xlim([0,SIZE])
  ax.set_xlabel("$r$ (M)")
  ax.axvline(AT_R, color='k')
  
  if "SANE" in run_name:
    ax.set_ylim([1e-4,1e-1])
    ax.set_yscale('log')
  else:
    ax.set_ylim([1,100])
    #ax.set_yscale('log')
  ax.legend(loc='lower right')
  
  plt.tight_layout()
  plt.savefig(os.path.join("frames_lmovie"+run_name.replace("/","_"),"frame_{:04}.png".format(n)), dpi=100)
  plt.close(fig)

if __name__ == "__main__":
  util.run_parallel(frame, len(dumps), util.calc_nthreads(hdr, pad=0.1))



