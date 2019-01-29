#!/usr/bin/env python3

import os, sys
import pickle
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

import util
import hdf5_to_dict as io
import plot as bplt
from analysis_fns import *

from defs import Met
from coordinates import dxdX_KS_to

FIGX=15
FIGY=15

# Decide where to measure fluxes
def i_of(rcoord):
  i = 0
  while geom['r'][i,hdr['n2']//2,0] < rcoord:
    i += 1
  i -= 1
  return i

def cut_pos(var, cut):
  var_X2 = bplt.flatten_xz(var, average=True)[iBZ,:hdr['n2']//2]
  var_cut = np.where( np.logical_or(
    np.logical_and(var_X2[:-1] > cut, var_X2[1:] < cut),
    np.logical_and(var_X2[:-1] < cut, var_X2[1:] > cut)))
  return var_cut

if __name__ == "__main__":
  run_name = sys.argv[1]
  dumpfile = os.path.join("/scratch/03002/bprather/pharm_dumps/M87SimulationLibrary/GRMHD",run_name,"dumps/dump_00001500.h5")
  hdr,geom,dump = io.load_all(dumpfile)
  
  plotfile = os.path.join("/work/03002/bprather/stampede2/movies",run_name,"eht_out.p")
  avg = pickle.load(open(plotfile, "rb"))
  
  # BZ luminosity; see eht_analysis
  if hdr['r_out'] < 100:
    iBZ = i_of(40) # most SANEs
    rBZ = 40
    rstring = "40"
  else:
    iBZ = i_of(100) # most MADs
    rBZ = 100
    rstring = "100"

  avg['X2'] = geom['X2'][iBZ,:hdr['n2']//2,0]
  
  ucon_cut = cut_pos(dump['ucon'][:,:,:,1], 0.0)
  
  sigma_cut1 = cut_pos(dump['bsq']/dump['RHO'], 1.0)
  sigma_cut10 = cut_pos(dump['bsq']/dump['RHO'], 10.0)

  be_b0_cut = cut_pos(bernoulli(dump, with_B=True), 0.02)
  be_b1_cut = cut_pos(bernoulli(dump, with_B=True), 1.0)
  be_nob0_cut = cut_pos(bernoulli(dump, with_B=False), 0.02)
  be_nob1_cut = cut_pos(bernoulli(dump, with_B=False), 1.0)

  rur_cut = cut_pos(geom['r']*dump['ucon'][:,:,:,1], 1.0)
  gamma_cut = cut_pos(get_gamma(geom, dump), 1.5)
  
  # For converting to theta
  startx1, hslope, poly_norm, poly_alpha, poly_xt, mks_smooth = geom['startx1'], geom['hslope'], geom['poly_norm'], geom['poly_alpha'], geom['poly_xt'], geom['mks_smooth']
  to_th = np.pi + (1. - hslope) * np.pi * np.cos(2. * np.pi * avg['X2']) + np.exp(
    mks_smooth * (startx1 - geom['X1'][iBZ,:hdr['n2']//2,0])) * (-np.pi +
                                             2. * poly_norm * (1. + np.power((2. * avg['X2'] - 1.) / poly_xt,
                                                                             poly_alpha) / (poly_alpha + 1.)) +
                                             (2. * poly_alpha * poly_norm * (2. * avg['X2'] - 1.) * np.power(
                                                 (2. * avg['X2'] - 1.) / poly_xt, poly_alpha - 1.)) / (
                                                         (1. + poly_alpha) * poly_xt) -
                                             (1. - hslope) * np.pi * np.cos(2. * np.pi * avg['X2']))
  # For averaging
  # What a god-awful hack.  Better to have preserved it with metadata thru the pipeline but whatever
  # NOTE: a-0.94 384 MAD should be 7000-10000, but due to a bug images were made from 5k, so we average where the images are
  qui = { "MAD" : { "384x192x192_IHARM" : { "a-0.94" : [5000,10000], "a-0.5" : [5000,9000], "a0" : [5000,10000], "a+0.5" : [5000,10000], "a+0.75" : [5000,9000], "a+0.94" : [5000,10000] },
                  "192x96x96_IHARM" : { "a-0.94" : [5000,10000], "a-0.5" : [5500,10000], "a-0.25" : [6000,9000], "a0" : [5000,8500], "a+0.25" : [6000,10000], "a+0.5" : [5000,10000], "a+0.75" : [6000,10000], "a+0.94" : [5000,10000] },
                  "288x128x128_gamma53" : { "a+0.94" : [6000,10000] },
                  "384x192x192_gamma53" : { "a+0.94" : [6000,9990] } # Last dump was not output
                  },
         "SANE" : { "288x128x128_IHARM" : { "a-0.94" : [6000,9000], "a-0.5" : [5000,8000], "a0" : [5000,10000], "a+0.5" : [3000,5500], "a+0.94" : [3000,6000] },
                  "192x192x192_IHARM" : { "a+0.5" : [6000,8000], "a+0.94" : [8000,10000] }
                  }
      }

  ND = avg['LBZ_sigma1_rt'].shape[0]
  rtype,rspin,rname = run_name.split("/")
  start, end = qui[rtype][rname][rspin]
  # I can rely on this for now
  start = int(start)//5
  end = int(end)//5

  
  # Perform average and conversion
  new_avg = {}
  for key in avg.keys():
    if "_tht" in key:
      new_avg[key[:-1]] = np.mean(avg[key][start:end,:], axis=0)*to_th*geom['gdet'][iBZ,:hdr['n2']//2]
  avg = {**avg, **new_avg}

  # L_th
  fig, axes = plt.subplots(2,2, figsize=(FIGX, FIGY))

  # Plot Luminosity contribution (incl. KE) as a fn of theta at r=100
  ax = axes[0,0]
  ax.plot(avg['th'], avg['Ltot_th'], color='xkcd:pink', label=r"$L_{tot}$ (r = "+rstring+")")
  ax.plot(avg['th'], avg['LBZ_th'], color='xkcd:cyan', label=r"$L_{BZ}$ (r = "+rstring+")")
  
  prop = np.sum(avg['Ltot_th'][np.where(avg['Ltot_th'] > 0)])/np.max(avg['Ltot_th'])
  Ltot_th_acc = [np.sum(avg['Ltot_th'][:n])/prop for n in range(hdr['n2']//2)]
  LBZ_th_acc = [np.sum(avg['LBZ_th'][:n])/prop for n in range(hdr['n2']//2)]
  ax.plot(avg['th'], Ltot_th_acc, 'C3', label=r"Accumulated $L_{tot}$ (r = "+rstring+")")
  ax.plot(avg['th'], LBZ_th_acc, 'C8', label=r"Accumulated $L_{BZ}$ (r = "+rstring+")")
  
  ax.axhline(0.0, color='k', linestyle=':')
  
  # Add axis marked with theta
#   def x2_to_theta(x2): return 180*geom['th'][iBZ,int(x2*hdr['n2']),0]/np.pi
#   formatterX = FuncFormatter(lambda x, pos: '{0:g}'.format(x2_to_theta(x)) if x2_to_theta(x) > 5 else '')
# 
#   ax2 = ax.twiny()
#   ax2.set_xlim(ax.get_xlim())
#   ax2.xaxis.set_major_formatter(formatterX)

  ymax = 1.1*max(np.max(np.abs(avg['Ltot_th'])), np.max(np.abs(avg['LBZ_th'])))
  #ymin = 1e-4*ymax
  ymin = -0.2*ymax
  #ymax = 0.1
  #ymin = -0.1

  ax.vlines(avg['th'][sigma_cut1], ymin, ymax, colors='C2', label="Sigma > 1 Cut")

  ax.vlines(avg['th'][be_nob0_cut], ymin, ymax, colors='C6', label="Be > 0.02 Cut")
  ax.vlines(avg['th'][be_nob1_cut], ymin, ymax, colors='C7', label="Be > 1.0 Cut")

  #ax.vlines(avg['th'][rur_cut], ymin, ymax, colors='C8')

  # Legend
  ax.legend(loc='upper right')
  
  ax.set_ylim([ymin, ymax])
  #ax.set_yscale('log')
  
  ax = axes[0,1]
  ax.plot(avg['th'], avg['RHO_g_th'], color='xkcd:purple', label=r"$\rho$ (r = "+rstring+")")
  ax.legend(loc='upper left')
  
  ax = axes[1,0]
  ax.plot(avg['th'], avg['FM_g_th'], color='xkcd:green', label=r"Mass flux (r = "+rstring+")")
  ax.plot(avg['th'], avg['FE_g_th'], color='xkcd:orange', label=r"Energy flux (r = "+rstring+")")
  ax.legend(loc='upper left')
  
  ax = axes[1,1]
  ax.plot(avg['th'], np.mean(avg['FL_g_tht'][1500:,:], axis=0), color='xkcd:blue', label=r"L flux (r = "+rstring+")")
  ax.legend(loc='upper left')

  plt.savefig(run_name.replace("/", "_") + '_L_th.png')
  plt.close(fig)

