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
from coordinates import dxdX_to_KS, dxdX_KS_to

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
  if var.ndim > 2:
    var_X2 = bplt.flatten_xz(var, average=True)[iBZ,:hdr['n2']//2]
  else:
    var_X2 = var
  var_cut = np.where( np.logical_or(
    np.logical_and(var_X2[:-1] > cut, var_X2[1:] < cut),
    np.logical_and(var_X2[:-1] < cut, var_X2[1:] > cut)))
  return var_cut

def overlay_thphi_contours(ax, geom, r):
  s = "_" + str(r) + "_thphi"
  r_i = i_of(r)
  max_th = geom['n2']//2
  x = bplt.loop_phi(geom['x'][r_i,:max_th,:])
  y = bplt.loop_phi(geom['y'][r_i,:max_th,:])
  prep = lambda var : bplt.loop_phi(var[:max_th,:])
  ax.contour(x,y, prep(avg['sigma'+s]), [1.0], colors='xkcd:blue')
  ax.contour(x,y, prep(avg['Be_nob'+s]), [0.02], colors='xkcd:purple')
  ax.contour(x,y, prep(avg['Be_nob'+s]), [1.0], colors='xkcd:green')

def overlay_rth_contours(ax, geom):
  s = '_rth'
  bplt.overlay_contours(ax, geom, geom['th'][:,:,0], [1.0, np.pi-1.0], color='k')
  bplt.overlay_contours(ax, geom, avg['sigma'+s], [1.0], color='xkcd:blue')
  bplt.overlay_contours(ax, geom, avg['Be_nob'+s], [0.02], color='xkcd:purple')
  bplt.overlay_contours(ax, geom, avg['Be_nob'+s], [1.0], color='xkcd:green')

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

  print(avg.keys())
  avg['sigma_100_th'] = avg['bsq_100_th']/avg['rho_100_th']
  sigma_cut1 = cut_pos(avg['sigma_100_th'][:hdr['n2']//2]/hdr['n3'], 1.0)
  sigma_cut10 = cut_pos(avg['sigma_100_th'][:hdr['n2']//2]/hdr['n3'], 10.0)
 
  be_b0_cut = cut_pos(avg['Be_b_100_th'][:hdr['n2']//2]/hdr['n3'], 0.02)
  be_b1_cut = cut_pos(avg['Be_b_100_th'][:hdr['n2']//2]/hdr['n3'], 1.0)
  be_nob0_cut = cut_pos(avg['Be_nob_100_th'][:hdr['n2']//2]/hdr['n3'], 0.02)
  be_nob1_cut = cut_pos(avg['Be_nob_100_th'][:hdr['n2']//2]/hdr['n3'], 1.0)
 
  rur_cut = cut_pos(avg['rur_100_th'][:hdr['n2']//2]/hdr['n3'], 1.0)
  gamma_cut = cut_pos(avg['gamma_100_th'][:hdr['n2']//2]/hdr['n3'], 1.5)
  
  # For converting to theta someday
#   Xgeom = np.zeros((4,1,geom['n2']))
#   Xgeom[1] = geom['X1'][iBZ,:,0]
#   Xgeom[2] = geom['X2'][iBZ,:,0]
#   to_dth = dxdX_KS_to(Xgeom, Met.FMKS, geom)[2,2,0]

  # For converting to theta now
  avg['X2'] = geom['X2'][0,:,0]
  startx1, hslope, poly_norm, poly_alpha, poly_xt, mks_smooth = geom['startx1'], geom['hslope'], geom['poly_norm'], geom['poly_alpha'], geom['poly_xt'], geom['mks_smooth']
  to_dth_bz = np.pi + (1. - hslope) * np.pi * np.cos(2. * np.pi * avg['X2']) + np.exp(
    mks_smooth * (startx1 - geom['X1'][iBZ,:,0])) * (-np.pi +
                                             2. * poly_norm * (1. + np.power((2. * avg['X2'] - 1.) / poly_xt,
                                                                             poly_alpha) / (poly_alpha + 1.)) +
                                             (2. * poly_alpha * poly_norm * (2. * avg['X2'] - 1.) * np.power(
                                                 (2. * avg['X2'] - 1.) / poly_xt, poly_alpha - 1.)) / (
                                                         (1. + poly_alpha) * poly_xt) -
                                                         (1. - hslope) * np.pi * np.cos(2. * np.pi * avg['X2']))
  to_dth_5 = np.pi + (1. - hslope) * np.pi * np.cos(2. * np.pi * avg['X2']) + np.exp(
    mks_smooth * (startx1 - geom['X1'][i_of(5),:,0])) * (-np.pi +
                                             2. * poly_norm * (1. + np.power((2. * avg['X2'] - 1.) / poly_xt,
                                                                             poly_alpha) / (poly_alpha + 1.)) +
                                             (2. * poly_alpha * poly_norm * (2. * avg['X2'] - 1.) * np.power(
                                                 (2. * avg['X2'] - 1.) / poly_xt, poly_alpha - 1.)) / (
                                                         (1. + poly_alpha) * poly_xt) -
                                                         (1. - hslope) * np.pi * np.cos(2. * np.pi * avg['X2']))

  to_dth_bz = 1/to_dth_bz
  to_dth_5 = 1/to_dth_5
  
  

  ND = avg['t'].shape[0]
  # I can rely on this for now
  start = int(avg['avg_start'])//5
  end = int(avg['avg_end'])//5

  # Perform average and conversion
  for key in avg.keys():
    if "100_th" in key:
      avg[key] *= geom['gdet'][iBZ,:]*hdr['dx3']*to_dth_bz
    if "100_thphi" in key:
      for i in range(hdr['n3']):
        avg[key][:,i] *= to_dth_bz
    if "5_th" in key:
      avg[key] *= geom['gdet'][iBZ,:]*hdr['dx3']*to_dth_5
    if "5_thphi" in key:
      for i in range(hdr['n3']):
        avg[key][:,i] *= to_dth_5

  avg['th100'] = geom['th'][iBZ,:,0]
  avg['hth100'] = geom['th'][iBZ,:hdr['n2']//2,0]
  avg['th5'] = geom['th'][i_of(5),:,0]
  
  with open("average_"+run_name.replace("/","_")+".dat", "w") as datf:
    datf.write("# x2 theta dx2/dtheta gdet rho bsq Fem_t Ffl_t F_mass\n")
    for i in range(hdr['n2']):
      datf.write("{} {} {} {} {} {} {} {} {}\n".format(avg['X2'][i], avg['th100'][i], to_dth_bz[i], geom['gdet'][iBZ,i],
                                                       avg['RHO_100_th'][i], avg['bsq_100_th'][i],
                                                       avg['FE_EM_100_th'][i], avg['FE_Fl_100_th'][i], avg['FM_100_th'][i]))
  
#  dth = hdr['dx2']/to_dth_bz
#   sigmax = sigma_cut1[0][0]
#   sigmin = sigma_cut1[0][1]
#   for flux,lum in [['FE_EM', 'LBZ'], ['FE', 'Ltot']]:
#     print("Compare {} Sigma > 1 cut {} to integrated value {}".format(lum,
#       np.mean(avg[lum+'_sigma1_rt'][start:end,iBZ]),
#       (np.sum(avg[flux+'_100_th'][:sigmax]*dth[:sigmax]) + np.sum(avg[flux+'_100_th'][sigmin:]*dth[sigmin:]))))

  # L_th
  fig, axes = plt.subplots(2,2, figsize=(FIGX, FIGY))

  # Plot Luminosity contribution (incl. KE) as a fn of theta at r=100
  ax = axes[0,0]
  ax.plot(avg['hth100'], avg['FE_100_th'][:hdr['n2']//2], color='C1', label=r"$FE_{tot}$ (r = "+rstring+")")
  ax.plot(avg['hth100'], avg['FE_EM_100_th'][:hdr['n2']//2], color='C2', label=r"$FE_{EM}$ (r = "+rstring+")")
  ax.plot(avg['hth100'], avg['FE_Fl_100_th'][:hdr['n2']//2], color='C3', label=r"$FE_{Fl}$ (r = "+rstring+")")

#   prop = np.sum(avg['Ltot_th'][np.where(avg['Ltot_th'] > 0)])/np.max(avg['Ltot_th'])
#   Ltot_th_acc = [np.sum(avg['Ltot_th'][:n])/prop for n in range(hdr['n2']//2)]
#   LBZ_th_acc = [np.sum(avg['LBZ_th'][:n])/prop for n in range(hdr['n2']//2)]
#   ax.plot(avg['th100'], Ltot_th_acc, 'C3', label=r"Accumulated $L_{tot}$ (r = "+rstring+")")
#   ax.plot(avg['th100'], LBZ_th_acc, 'C8', label=r"Accumulated $L_{BZ}$ (r = "+rstring+")")

  ax.axhline(0.0, color='k', linestyle=':')

  ymax = 1.1*max(np.max(avg['FE_100_th']), np.max(avg['FE_EM_100_th']))
  ymin = 1e-2
  ax.set_ylim([ymin, ymax])
  ax.set_yscale('log')

  ax.vlines(avg['th100'][sigma_cut1], ymin, ymax, colors='xkcd:blue', label="Sigma > 1 Cut")
  ax.vlines(avg['th100'][be_nob0_cut], ymin, ymax, colors='xkcd:purple', label="Be > 0.02 Cut")
  ax.vlines(avg['th100'][be_nob1_cut], ymin, ymax, colors='xkcd:green', label="Be > 1.0 Cut")

  ax.legend(loc='upper right')
  
  ax = axes[0,1]
  ax.plot(avg['hth100'], avg['FL_100_th'][:hdr['n2']//2], color='C1', label=r"FL (r = "+rstring+")")
  ax.plot(avg['hth100'], avg['FL_EM_100_th'][:hdr['n2']//2], color='C2', label=r"FL (r = "+rstring+")")
  ax.plot(avg['hth100'], avg['FL_Fl_100_th'][:hdr['n2']//2], color='C3', label=r"FL (r = "+rstring+")")
  
  ymax = 1.1*max(np.max(avg['FL_100_th']), np.max(avg['FL_EM_100_th']))
  ymin = 1e-2
  ax.set_ylim([ymin, ymax])
  ax.set_yscale('log')

  ax.vlines(avg['th100'][sigma_cut1], ymin, ymax, colors='xkcd:blue', label="Sigma > 1 Cut")
  ax.vlines(avg['th100'][be_nob0_cut], ymin, ymax, colors='xkcd:purple', label="Be > 0.02 Cut")
  ax.vlines(avg['th100'][be_nob1_cut], ymin, ymax, colors='xkcd:green', label="Be > 1.0 Cut")
  
  
  ax.legend(loc='upper left')
  
  ax = axes[1,0]
  ax.plot(avg['th100'], avg['RHO_100_th'], color='xkcd:purple', label=r"$\rho$ (r = "+rstring+")")
  ax.set_yscale('log')
  ax.legend(loc='upper left')

  plt.savefig(run_name.replace("/", "_") + '_L_th.png')
  plt.close(fig)
  
  fig, ax = plt.subplots(2,2,figsize=(FIGX, FIGY))
  bplt.plot_thphi(ax[0,0], geom, np.log10(avg['FE_100_thphi']), iBZ, label = "FE 2D Slice r="+rstring)
  overlay_thphi_contours(ax[0,0], geom, 100)
  bplt.plot_thphi(ax[0,1], geom, np.log10(avg['FM_100_thphi']), iBZ, label = "FM 2D Slice r="+rstring)
  overlay_thphi_contours(ax[0,1], geom, 100)
  bplt.plot_thphi(ax[1,0], geom, np.log10(avg['FL_100_thphi']), iBZ, label = "FL 2D Slice r="+rstring)
  overlay_thphi_contours(ax[1,0], geom, 100)
  bplt.plot_thphi(ax[1,1], geom, np.log10(avg['RHO_100_thphi']), iBZ, label = "RHO 2D Slice r="+rstring)
  overlay_thphi_contours(ax[1,1], geom, 100)
  
  plt.savefig(run_name.replace("/", "_") + '_L_100_thphi.png')
  plt.close(fig)
  
  fig, ax = plt.subplots(2,2,figsize=(FIGX, FIGY))
  bplt.plot_thphi(ax[0,0], geom, np.log10(avg['FE_5_thphi']), i_of(5), label = "FE 2D Slice r=5")
  overlay_thphi_contours(ax[0,0], geom, 5)
  bplt.plot_thphi(ax[0,1], geom, np.log10(avg['FM_5_thphi']), i_of(5), label = "FM 2D Slice r=5")
  overlay_thphi_contours(ax[0,1], geom, 5)
  bplt.plot_thphi(ax[1,0], geom, np.log10(avg['FL_5_thphi']), i_of(5), label = "FL 2D Slice r=5")
  overlay_thphi_contours(ax[1,0], geom, 5)
  bplt.plot_thphi(ax[1,1], geom, np.log10(avg['RHO_5_thphi']), i_of(5), label = "RHO 2D Slice r=5")
  overlay_thphi_contours(ax[1,1], geom, 5)
  
  plt.savefig(run_name.replace("/", "_") + '_L_5_thphi.png')
  plt.close(fig)
  
  fig, ax = plt.subplots(2,2,figsize=(FIGX, FIGY))
  bplt.plot_xz(ax[0,0], geom, np.log10(avg['FE_rth']), label = "FE 2D Slice")
  overlay_rth_contours(ax[0,0], geom)
  bplt.plot_xz(ax[0,1], geom, np.log10(avg['FM_rth']), label = "FM 2D Slice")
  overlay_rth_contours(ax[0,1], geom)
  bplt.plot_xz(ax[1,0], geom, np.log10(avg['FL_rth']), label = "FL 2D Slice")
  overlay_rth_contours(ax[1,0], geom)
  bplt.plot_xz(ax[1,1], geom, np.log10(avg['RHO_rth']), label = "RHO 2D Slice")
  overlay_rth_contours(ax[1,1], geom)
  
  plt.savefig(run_name.replace("/", "_") + '_L_rth.png')
  plt.close(fig)
  
  
  
  
  
  

