################################################################################
#                                                                              #
#  PLOTS OF VARIABLES COMPUTED IN eht_analysis.py                              #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('Agg')

import sys; sys.dont_write_bytecode=True
import numpy as np
import matplotlib.pyplot as plt
import util
import pickle

FIGX = 14
FIGY = 14
EXTRAS = True
DIAGS = False
RADS = True

if len(sys.argv) != 3:
  util.warn('Format: python eht_plot.py [averages] [output]')
  sys.exit()

favg = sys.argv[1]
fout = sys.argv[2]

avg = pickle.load(open(favg, 'rb'))

ti = avg['t'][0]
tf = avg['t'][-1]

if RADS:
  # Radial averages
  fig = plt.figure(figsize=(FIGX, FIGY))

  ax = plt.subplot(2,3,1)
  ax.plot(avg['r'], avg['rho_r'], color='k', linewidth=2)
  ax.set_yscale('log')
  ax.set_ylabel('<rho>')
  ax.set_ylim([1.e-3, 1.e0])

  ax = plt.subplot(2,3,2)
  ax.plot(avg['r'], avg['Pg_r'], color='k', linewidth=2)
  ax.set_yscale('log')
  ax.set_ylabel('<Pg>')
  ax.set_ylim([1.e-6, 1.e-2])

  ax = plt.subplot(2,3,3)
  ax.plot(avg['r'], avg['B_r'], color='k', linewidth=2)
  ax.set_yscale('log')
  ax.set_ylabel('<|B|>')
  ax.set_ylim([1.e-4, 1.e-1])

  ax = plt.subplot(2,3,4)
  ax.plot(avg['r'], avg['uphi_r'], color='k', linewidth=2)
  ax.set_yscale('log')
  ax.set_ylabel('<u^phi>')
  ax.set_ylim([1.e-3, 1.e1])

  ax = plt.subplot(2,3,5)
  ax.plot(avg['r'], avg['Ptot_r'], color='k', linewidth=2)
  ax.set_yscale('log')
  ax.set_ylabel('<Ptot>')
  ax.set_ylim([1.e-6, 1.e-2])

  ax = plt.subplot(2,3,6)
  ax.plot(avg['r'], avg['betainv_r'], color='k', linewidth=2)
  ax.set_yscale('log')
  ax.set_ylabel('<beta^-1>')
  ax.set_ylim([1.e-2, 1.e1])

  plt.savefig(fout + '_ravgs.png')
  plt.close(fig)

  # SADW
  fig = plt.figure(figsize=(FIGX, FIGY))
  ax = plt.subplot(2,3,1)
  ax.plot(avg['r'], avg['rho_SADW'], color='k', linewidth=2)
  ax.set_yscale('log')
  plt.savefig(fout + '_sadw.png')
  plt.close(fig)

if EXTRAS:
  # RADIAL PROFILES OF FLUXES
  fig = plt.figure(figsize=(FIGX, FIGY))
  ax = plt.subplot(2,1,1)
  ax.plot(avg['r'], avg['FE_r'], color='k', linewidth=2)
  ax.set_xlim([2,20])
  ax.set_ylabel('FE_r')
  ax.set_ylim([-150, 100])
  
  ax = plt.subplot(2,1,2)
  ax.plot(avg['r'], avg['FM_r'], color='k', linewidth=2)
  ax.set_xlim([2,20])
  ax.set_ylabel('FM_r')
  ax.set_ylim([0, 100])
  
  plt.savefig(fout + '_fluxr.png')
  plt.close(fig)
  
  # Omega
  #print avg['th'].shape, avg['omega_th'].shape
  fig = plt.figure(figsize=(FIGX, FIGY))
  ax = plt.subplot(2,2,1)
  ax.plot(avg['th'], avg['omega_th']*4/avg['a'], color='k')
  ax.set_ylim([-1,2])
  ax.set_ylabel("omega (using t r component)")
  ax.set_xlabel("theta")
  ax.axhline(0.5, linestyle='--', color='k')
  
  ax = plt.subplot(2,2,2)
  ax.plot(avg['th'], avg['omega_th_av']*4/avg['a'], color='k')
  ax.set_ylim([-1,2])
  ax.set_ylabel("omega (using t r component, averaged)")
  ax.set_xlabel("theta")
  ax.axhline(0.5, linestyle='--', color='k')
  
#  ax = plt.subplot(2,2,3)
#  ax.plot(avg['th'], avg['omega_th_alt']*4/avg['a'], color='k')
#  ax.set_ylim([-1,2])
#  ax.set_ylabel("omega (using t theta component)")
#  ax.set_xlabel("theta")
#  ax.axhline(0.5, linestyle='--', color='k')
  
#  ax = plt.subplot(2,2,4)
#  ax.plot(avg['th'], avg['omega_th_alt_av']*4/avg['a'], color='k')
#  ax.set_ylim([-1,2])
#  ax.set_ylabel("omega (using t theta component, averaged)")
#  ax.set_xlabel("theta")
#  ax.axhline(0.5, linestyle='--', color='k')
  
  plt.savefig(fout + '_omega.png')
  plt.close(fig)

## FLUXES
if EXTRAS:
  fig = plt.figure(figsize=(FIGX, FIGY*1.4))
  nplots = 7
else:
  fig = plt.figure(figsize=(FIGX, FIGY))
  nplots = 5

ax = plt.subplot(nplots,1,1)
ax.plot(avg['t'], np.fabs(avg['Mdot']), color='k')
if DIAGS: ax.plot(avg['t_d'], np.fabs(avg['Mdot_d']), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('|Mdot|')

ax = plt.subplot(nplots,1,2)
ax.plot(avg['t'], avg['Phi'], color='k')
if DIAGS: ax.plot(avg['t_d'], avg['Phi_d'], linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('Phi')

ax = plt.subplot(nplots,1,3)
ax.plot(avg['t'], np.fabs(avg['Ldot']), color='k')
if DIAGS: ax.plot(avg['t_d'], np.fabs(avg['Ldot_d']), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('|Ldot|')

ax = plt.subplot(nplots,1,4)
ax.plot(avg['t'], np.fabs(avg['Edot'] - np.fabs(avg['Mdot'])), color='k')
if DIAGS: ax.plot(avg['t_d'], np.fabs(avg['Edot_d'] - np.fabs(avg['Mdot_d'])), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('|Edot - Mdot|')

ax = plt.subplot(nplots,1,5)
ax.plot(avg['t'], avg['Lum'], color='k')
if DIAGS: ax.plot(avg['t_d'], avg['Lum_d'], linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.grid(axis='y')
ax.set_ylabel('Lum')

if not EXTRAS:
  ax.set_xlabel('t/M')
else:
  ax.set_xticklabels([])

  ax = plt.subplot(nplots,1,6)
  ax.plot(avg['t'], avg['Edot'], color='k')
  if DIAGS: ax.plot(avg['t_d'], avg['Edot_d'], linestyle='--', color='r')
  ax.set_xlim([ti,tf])
  ax.set_xlabel('t/M')
  ax.grid(axis='y')
  ax.set_ylabel('BZ Luminosity')
  
  ax = plt.subplot(nplots,1,7)
  ax.plot(avg['t'], -avg['LBZ'], color='k')
  ax.set_xlim([ti,tf])
  ax.set_xlabel('t/M')
  ax.grid(axis='y')
  ax.set_ylabel('BZ Luminosity')

plt.savefig(fout + '_fluxes.png')
plt.close(fig)

## NORMALIZED FLUXES

if EXTRAS:
  fig = plt.figure(figsize=(FIGX, FIGY*1.2))
  nplots = 6
else:
  fig = plt.figure(figsize=(FIGX, FIGY))
  nplots = 5

ax = plt.subplot(nplots,1,1)
ax.plot(avg['t'], np.fabs(avg['Mdot']), color='k')
if DIAGS: ax.plot(avg['t_d'], np.fabs(avg['Mdot_d']), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('|Mdot|')

ax = plt.subplot(nplots,1,2)
ax.plot(avg['t'], avg['Phi']/np.sqrt(np.fabs(avg['Mdot'])), color='k')
if DIAGS: ax.plot(avg['t_d'], avg['Phi_d']/np.sqrt(np.fabs(avg['Mdot_d'])), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('Phi/sqrt(|Mdot|)')

ax = plt.subplot(nplots,1,3)
ax.plot(avg['t'], np.fabs(avg['Ldot'])/(np.fabs(avg['Mdot'])), color='k')
if DIAGS: ax.plot(avg['t_d'], np.fabs(avg['Ldot_d'])/(np.fabs(avg['Mdot_d'])), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('|Ldot|/|Mdot|')

ax = plt.subplot(nplots,1,4)
ax.plot(avg['t'], np.fabs(avg['Edot'] - np.fabs(avg['Mdot']))/(np.fabs(avg['Mdot'])), color='k')
if DIAGS: ax.plot(avg['t_d'], np.fabs(avg['Edot_d'] - np.fabs(avg['Mdot_d']))/(np.fabs(avg['Mdot_d'])), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.set_xticklabels([])
ax.grid(axis='y')
ax.set_ylabel('|Edot - Mdot|/|Mdot|')

ax = plt.subplot(nplots,1,5)
ax.plot(avg['t'], avg['Lum']/(np.fabs(avg['Mdot'])), color='k')
if DIAGS: ax.plot(avg['t_d'], avg['Lum_d']/(np.fabs(avg['Mdot_d'])), linestyle='--', color='r')
ax.set_xlim([ti,tf])
ax.grid(axis='y')
ax.set_ylabel('Lum/|Mdot|')

if not EXTRAS:
  ax.set_xlabel('t/M')
else:
  ax.set_xticklabels([])

  ax = plt.subplot(nplots,1,6)
  # Just use the back half as <>
  mdot_mean = np.mean(avg['Mdot'][len(avg['Mdot'])/2:])
  ax.plot(avg['t'], (avg['Mdot'] - avg['Edot'])/mdot_mean*100, color='k')
  if DIAGS:
    mdot_mean_d = np.mean(avg['Mdot'][len(avg['Mdot_d'])/2:])
    ax.plot(avg['t_d'], (avg['Mdot_d'] - np.fabs(avg['Edot_d']))/mdot_mean_d*100, linestyle='--', color='r')
  ax.set_xlim([ti,tf])
  ax.set_xlabel('t/M')
  ax.set_ylim([-10,200])
  ax.grid(axis='y')
  ax.set_ylabel('Efficiency (%)')

plt.savefig(fout + '_normfluxes.png')
plt.close(fig)

## Values which should be boring

if EXTRAS:
  fig = plt.figure(figsize=(FIGX, FIGY/2))
  nplots = 3

  ax = plt.subplot(nplots,1,1)
  ax.plot(avg['t'], avg['Etot'], color='k')
  ax.set_xlim([ti,tf])
  ax.set_xticklabels([])
  ax.grid(axis='y')
  ax.set_ylabel('Etot')

  ax = plt.subplot(nplots,1,2)
  ax.plot(avg['t'], np.fabs(avg['sigma_max']), color='k')
  ax.set_xlim([ti,tf])
  ax.set_xticklabels([])
  ax.grid(axis='y')
  ax.set_ylabel('sigma_max')

  ax = plt.subplot(nplots,1,3)
  ax.plot(avg['t_d'], avg['divbmax_d'], color='k')
  ax.set_xlim([ti,tf])
  ax.set_xlabel('t/M')
  ax.grid(axis='y')
  ax.set_ylabel('max divB')

  plt.savefig(fout + '_diagnostics.png')
  plt.close(fig)
