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
FIGY = 10
SIZE = 40

if len(sys.argv) < 4:
  util.warn('Format: python overplot.py [averages] [output name]')
  sys.exit()

# All interior arguments are files to overplot
avgs = []
for filename in sys.argv[1:-1]:
  avgs.append(pickle.load(open(filename,'rb')))

fname_out = sys.argv[-1]
fig = plt.figure(figsize=(FIGX, FIGY))

# For time plots.  Also take MAD/SANE for axes?
tf = 1e4

styles = ['k','r--','g-.','b:'][:len(sys.argv)-2] # TODO extend
def plot_all(ax, iname, varname, timelabels=False):
  for i,avg in enumerate(avgs):
    ax.plot(avg[iname], avg[varname], styles[i])
  if iname == 't':
    ax.set_xlim([0,tf])
    if not timelabels:
      ax.set_xticklabels([])
    
  # TODO legend?


ax = plt.subplot(2,3,1)
plot_all(ax, 'r', 'rho_r')
ax.set_yscale('log')
ax.set_ylabel('<rho>')
ax.set_ylim([1.e-2, 1.e0])

ax = plt.subplot(2,3,2)
plot_all(ax, 'r', 'Pg_r')
ax.set_yscale('log')
ax.set_ylabel('<Pg>')
ax.set_ylim([1.e-6, 1.e-2])

ax = plt.subplot(2,3,3)
plot_all(ax, 'r', 'B_r')
ax.set_yscale('log')
ax.set_ylabel('<|B|>')
ax.set_ylim([1.e-4, 1.e-1])

ax = plt.subplot(2,3,4)
plot_all(ax, 'r', 'uphi_r')
ax.set_yscale('log')
ax.set_ylabel('<u^phi>')
ax.set_ylim([1.e-3, 1.e1])

ax = plt.subplot(2,3,5)
plot_all(ax, 'r', 'Ptot_r')
ax.set_yscale('log')
ax.set_ylabel('<Ptot>')
ax.set_ylim([1.e-6, 1.e-2])

ax = plt.subplot(2,3,6)
plot_all(ax, 'r', 'betainv_r')
ax.set_yscale('log')
ax.set_ylabel('<beta^-1>')
ax.set_ylim([1.e-2, 1.e1])

plt.savefig(fname_out + '_ravgs.png')
plt.close(fig)

# SADW
fig = plt.figure(figsize=(FIGX, FIGY))
ax = plt.subplot(2,3,1)
plot_all(ax, 'r', 'rho_SADW')
ax.set_yscale('log')
plt.savefig(fname_out + '_sadw.png')
plt.close(fig)

# Omega
fig = plt.figure(figsize=(FIGX, FIGY))
ax = plt.subplot(2,2,1)
# Renormalize omega to omega/Omega_H for plotting
for avg in avgs:
  avg['omega_th'] *= 4/avg['a']
  avg['omega_th_av'] *= 4/avg['a']
plot_all(ax, 'th', 'omega_th')
ax.set_ylim([-1,2])
ax.set_ylabel("omega (EH, time/hemi only)")
ax.set_xlabel("theta")
ax.axhline(0.5, linestyle='--', color='k')

ax = plt.subplot(2,2,2)
plot_all(ax, 'th', 'omega_th_av')
ax.set_ylim([-1,2])
ax.set_ylabel("omega (EH, 5-zone average in R)")
ax.set_xlabel("theta")
ax.axhline(0.5, linestyle='--', color='k')

plt.savefig(fname_out + '_omega.png')
plt.close(fig)

fig = plt.figure(figsize=(FIGX, FIGY))

ax = plt.subplot(5,1,1)
plot_all(ax,'t','Mdot')
ax.set_ylabel('|Mdot|')

ax = plt.subplot(5,1,2)
for avg in avgs:
  if np.max(avg['Phi']) < 2:
    avg['Phi'] *= np.sqrt(4*np.pi)
plot_all(ax,'t','Phi')
ax.set_ylabel('Phi')

ax = plt.subplot(5,1,3)
for avg in avgs:
  avg['abs_Ldot'] = np.abs(avg['Ldot'])
plot_all(ax,'t','abs_Ldot')
ax.set_ylabel('|Ldot|')

ax = plt.subplot(5,1,4)
for avg in avgs:
  avg['EmM'] = np.fabs(np.fabs(avg['Edot']) - np.fabs(avg['Mdot']))
plot_all(ax, 't', 'EmM')
ax.set_ylabel('|Edot - Mdot|')

ax = plt.subplot(5,1,5)
plot_all(ax, 't', 'Lum', timelabels=True)
ax.set_xlabel('t/M')
ax.set_ylabel('Lum')

plt.savefig(fname_out + '_fluxes.png')
plt.close(fig)


fig = plt.figure(figsize=(FIGX, FIGY))

ax = plt.subplot(5,1,1)
plot_all(ax,'t','Mdot')
ax.set_ylabel('|Mdot|')

ax = plt.subplot(5,1,2)
for avg in avgs:
  avg['Phi_norm'] = avg['Phi']/np.sqrt(np.fabs(avg['Mdot']))
plot_all(ax, 't', 'Phi_norm')
ax.set_ylabel('Phi/sqrt(|Mdot|)')

ax = plt.subplot(5,1,3)
for avg in avgs:
  avg['Ldot_norm'] = np.fabs(avg['Ldot'])/np.fabs(avg['Mdot'])
plot_all(ax, 't', 'Ldot_norm')
ax.set_ylabel('|Ldot|/|Mdot|')

ax = plt.subplot(5,1,4)
for avg in avgs:
  avg['Edot_norm'] = np.fabs(np.fabs(avg['Edot']) - np.fabs(avg['Mdot']))/(np.fabs(avg['Mdot']))
plot_all(ax, 't', 'Edot_norm')
ax.set_ylabel('|Edot - Mdot|/|Mdot|')

ax = plt.subplot(5,1,5)
for avg in avgs:
  avg['Lum_norm'] = np.fabs(avg['Lum'])/np.fabs(avg['Mdot'])
plot_all(ax, 't', 'Lum_norm', timelabels=True)
ax.set_xlabel('t/M')
ax.set_ylabel('Lum/|Mdot|')

plt.savefig(fname_out + '_normfluxes.png')
plt.close(fig)
