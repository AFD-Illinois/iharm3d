################################################################################
#                                                                              #
#  PLOTS OF VARIABLES COMPUTED IN eht_analysis.py                              #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('Agg')

import util

import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle

# For radials
FIGX = 14
FIGY = 10
# For flux plots
PLOTY = 3
SIZE = 40

RADS = True
OMEGA = True
FLUXES = True
EXTRAS = True
DIAGS = True

if len(sys.argv) < 3:
  util.warn('Format: python eht_plot.py analysis_output [analysis_output ...] [labels_list] image_name')
  sys.exit()

# All interior arguments are files to overplot
avgs = []
for filename in sys.argv[1:-2]:
  avgs.append(pickle.load(open(filename,'rb')))

labels = sys.argv[-2].split(",")

if len(labels) < len(avgs):
  util.warn("Too few labels!")
  sys.exit()

fname_out = sys.argv[-1]

# For time plots.  Also take MAD/SANE for axes?
ti = avgs[0]['t'][0]
tf = avgs[0]['t'][-1]

# Default styles
styles = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']

def plot_all(ax, iname, varname, varname_pretty, logx=False, logy=False, ylim=None, timelabels=False):
  for i,avg in enumerate(avgs):
    ax.plot(avg[iname], avg[varname], styles[i], label=labels[i])
  if logx: ax.set_xscale('log')
  if logx: ax.set_yscale('log')
  if ylim is not None: ax.set_ylim(ylim)
  ax.grid()
  if iname == 't':
    ax.set_xlim([ti,tf])
    if timelabels:
      ax.set_xlabel('t/M')
    else:
      ax.set_xticklabels([])
  elif iname == 'th':
    ax.set_xlabel(r"$\theta$")
  elif iname == 'r':
    ax.set_xlabel("r")
    

if RADS:
  fig, ax = plt.subplots(2,3, figsize=(FIGX, FIGY))
  plot_all(ax[0], 'r', 'rho_r', r"$<\rho>$", logy=True, ylim=[1.e-2, 1.e0])
  plot_all(ax, 'r', 'Pg_r', '<Pg>', logy=True, ylim=[1.e-6, 1.e-2])
  plot_all(ax, 'r', 'B_r', '<|B|>', logy=True, ylim=[1.e-4, 1.e-1])
  if len(avgs) > 1: ax.legend(loc=1)
  plot_all(ax, 'r', 'uphi_r', '<u^phi>', logy=True, ylim=[1.e-3, 1.e1])
  plot_all(ax, 'r', 'Ptot_r', '<Ptot>', logy=True, ylim=[1.e-6, 1.e-2])
  plot_all(ax, 'r', 'betainv_r', '<beta^-1>', logy=True, ylim=[1.e-2, 1.e1])

  plt.savefig(fname_out + '_ravgs.png')
  plt.close(fig)

if OMEGA:
  # Omega
  fig, ax = plt.subplots(2,2, figsize=(FIGX, FIGY))
  # Renormalize omega to omega/Omega_H for plotting
  for avg in avgs:
    avg['omega_th'] *= 4/avg['a']
    avg['omega_th_av'] *= 4/avg['a']
  plot_all(ax[0], 'th', 'omega_th', r"$\omega_f$ (EH, single shell)", ylim=[-1,2])
  if len(avgs) > 1: ax[0].legend(loc=3)
  
  plot_all(ax[1], 'th', 'omega_th_av', r"$\omega_f$ (EH, 5-zone average)")
  
  # Horizontal guidelines
  for a in ax:
    a.axhline(0.5, linestyle='--', color='k')
  
  plt.savefig(fname_out + '_omega.png')
  plt.close(fig)

if FLUXES:
  fig,ax = plt.subplots(5,1, figsize=(FIGX, 5*PLOTY))
  
  plot_all(ax[0],'t','Mdot', r"$|\dot{M}|$")
  if len(avgs) > 1: ax.legend(loc=2)
  
  plot_all(ax[1],'t','Phi', r"$\Phi$")
  
  for avg in avgs:
    avg['abs_Ldot'] = np.abs(avg['Ldot'])
  plot_all(ax[2],'t','abs_Ldot', r"$|\dot{L}|$")
  
  for avg in avgs:
    avg['EmM'] = np.fabs(np.fabs(avg['Edot']) - np.fabs(avg['Mdot']))
  plot_all(ax[3], 't', 'EmM', r"$|\dot{E} - \dot{M}|$")
  
  plot_all(ax[4], 't', 'Lum', "Lum", timelabels=True)
  
  plt.savefig(fname_out + '_fluxes.png')
  plt.close(fig)

  fig, ax = plt.subplots(5,1, figsize=(FIGX, 5*PLOTY))
  plot_all(ax[0], 't', 'Mdot', r"$|\dot{M}|$")
  if len(avgs) > 1: ax.legend(loc=2)
  
  for avg in avgs:
    avg['Phi_norm'] = avg['Phi']/np.sqrt(np.fabs(avg['Mdot']))
  plot_all(ax[1], 't', 'Phi_norm', r"$\frac{\Phi}{\sqrt{|\dot{M}|}}$")
  
  for avg in avgs:
    avg['Ldot_norm'] = np.fabs(avg['Ldot'])/np.fabs(avg['Mdot'])
  plot_all(ax[2], 't', 'Ldot_norm', r"$\frac{|Ldot|}{|\dot{M}|}$")
  
  for avg in avgs:
    avg['Edot_norm'] = np.fabs(np.fabs(avg['Edot']) - np.fabs(avg['Mdot']))/(np.fabs(avg['Mdot']))
  plot_all(ax[3], 't', 'Edot_norm', r"$\frac{|\dot{E} - \dot{M}|}{|\dot{M}|}$")
  
  for avg in avgs:
    avg['Lum_norm'] = np.fabs(avg['Lum'])/np.fabs(avg['Mdot'])
  plot_all(ax[4], 't', 'Lum_norm', r"$\frac{Lum}{|\dot{M}|}$", timelabels=True)
  
  plt.savefig(fname_out + '_normfluxes.png')
  plt.close(fig)

if EXTRAS:
  fig, ax = plt.subplots(3,1, figsize=(FIGX, 3*PLOTY))
  
  # Efficiency over mean mdot: just use the back half as <>
  for avg in avgs:
    mdot_mean = np.mean(avg['Mdot'][len(avg['Mdot'])/2:])
    avg['Eff'] = (avg['Mdot'] - avg['Edot'])/mdot_mean*100
  plot_all(ax[0], 't', 'Eff', "Efficiency (%)", ylim=[-10,200])

  ax.plot(ax[1], 't', 'Edot', r"\dot{E}")
  
  for avg in avgs:
    avg['aLBZ'] = np.abs(avg['LBZ'])
  ax.plot(ax[2], 't', 'aLBZ', "BZ Luminosity")

  plt.savefig(fout + '_extras.png')
  plt.close(fig)

if DIAGS:
  fig, ax = plt.subplots(3,1, figsize=(FIGX, FIGY/2))

  plot_all(ax[0], 't', 'Etot', "Total E")
  plot_all(ax[1], 't', 'sigma_max', r"$\sigma_{max}$")
  plot_all(ax[2], 't_d', 'divbmax_d', "max divB")

  plt.savefig(fout + '_diagnostics.png')
  plt.close(fig)
