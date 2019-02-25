################################################################################
#                                                                              #
#  PLOTS OF VARIABLES COMPUTED IN eht_analysis.py                              #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('Agg')

import util
import units

import sys
import numpy as np
import matplotlib.pyplot as plt
import pickle

# For radials
FIGX = 14
FIGY = 10
# For flux plots; per-plot Y dim
PLOTY = 3
SIZE = 40

RADS = True
FLUXES = True
EXTRAS = True
DIAGS = True
OMEGA = False
FLUX_PROF = False
TEMP = True
BSQ = True

def i_of(avg, rcoord):
  i = 0
  while avg['r'][i] < rcoord:
    i += 1
  i -= 1
  return i

def plot_multi(ax, iname, varname, varname_pretty, logx=False, logy=False, xlim=None, ylim=None, timelabels=False, label_list=None, linestyle='-'):
  if label_list is None: label_list = labels
  for i,avg in enumerate(avgs):
    if varname in avg.keys():
      if iname[:2] == "th":
        # We plot half-theta, from 0 to pi/2
        ax.plot(avg[iname][:avg[iname].size//2], avg[varname], styles[i]+linestyle, label=label_list[i])
      else:
        ax.plot(avg[iname], avg[varname], styles[i]+linestyle, label=label_list[i])
  # Plot additions
  if logx: ax.set_xscale('log')
  if logy: ax.set_yscale('log')
  if ylim is not None: ax.set_ylim(ylim)
  if xlim is not None: ax.set_xlim(xlim)
  ax.grid(True)
  ax.set_ylabel(varname_pretty)
  # Defaults and labels for different plot types:
  if iname == 't':
    if xlim is None:
      ax.set_xlim([ti,tf])
    if timelabels:
      ax.set_xlabel('t/M')
    else:
      ax.set_xticklabels([])
  elif iname == 'th':
    ax.set_xlabel(r"$\theta$")
  elif iname == 'r':
    ax.set_xlabel("r")
    if xlim is None:
      ax.set_xlim([0,50]) # For EHT comparison
  if logy:
    ylim = ax.get_ylim()
    if ylim[0] < 1e-5*ylim[1]:
      ax.set_ylim([1e-5*ylim[1], ylim[1]])

def plot_temp():
  fig, ax = plt.subplots(1,1, figsize=(FIGX, FIGY))
  if avgs[0]['r'][-1] > 50:
    txlim = [1e0,1e3]
  else:
    txlim = [1e0,1e2]
  
  fit_labs = []
  for i,avg in enumerate(avgs):
    cgs = units.get_cgs()
    avg['Tp_r'] = cgs['MP'] * avg['Pg_r'] / (cgs['KBOL'] * avg['rho_r']) * cgs['CL']**2
  
    # Add the fits. Aaaaaalll the fits
    x = avg['r'][i_of(avg, 3):i_of(avg, 30)]
    y = avg['Tp_r'][i_of(avg, 3):i_of(avg, 30)]
    logx=np.log(x)
    logy=np.log(y)
    # Place a guideline at 1/x starting near the top of the plot
    #ax.plot(logx, (0.1*ax.get_ylim()[1]) * 1/logx, 'k--')
    coeffs = np.polyfit(logx,logy,deg=1)
    poly = np.poly1d(coeffs)
    yfit = lambda x: np.exp(poly(np.log(x)))
    # TODO any if statements to cut down on fits
    avg['r_fit'] = x
    avg['Tp_r_fit'] = yfit(x)
    fit_lab = r"{:.2g} * r^{:.2g}".format(np.exp(coeffs[1]), coeffs[0])
    print(labels[i]," Ti fit: ",fit_lab)
    fit_labs.append(fit_lab)

  # Plot the profiles themselves
  plot_multi(ax, 'r', 'Tp_r', r"$<T_{i}>$", logx=True, xlim=txlim, logy=True)
  plot_multi(ax, 'r_fit', 'Tp_r_fit', r"$<T_{i}>$", logx=True, xlim=txlim, logy=True, label_list=fit_labs, linestyle='--')

  ax.legend(loc='lower right')
  plt.savefig(fname_out + "_Ti.png")
  plt.close(fig)

def plot_bsq_rise():
  fig, ax = plt.subplots(4,1, figsize=(FIGX, FIGY))
  for avg in avgs:
    if 'B_rt' in avg:
      avg['MagE'] = np.mean(avg['B_rt']**2, axis=-1)
      avg['MagE_close'] = np.mean(avg['B_rt'][:,:100]**2, axis=-1)
      avg['rho_close'] = np.mean(avg['rho_rt'][:,:100], axis=-1)
      n2 = avg['rho_100_tht'].shape[1]
      avg['rho_mid'] = np.mean(avg['rho_100_tht'][:,n2//2-5:n2//2+5], axis=-1)
    plot_multi(ax[0], 't', 'MagE_close', r"$<B^2>$ (low radii)", logy=True, xlim=[0,10000])
    plot_multi(ax[1], 't', 'rho_close', r"$<\rho>$ (low radii)", logy=True, xlim=[0,10000])
    plot_multi(ax[2], 't', 'MagE', r"$<B^2>$", logy=True, xlim=[0,10000])
    plot_multi(ax[3], 't', 'rho_mid', r"$<\rho>$ (midplane r=40)", logy=True, timelabels=True, xlim=[0,10000])
  
  #ax.legend(loc='lower left')
  plt.savefig(fname_out + '_bsq_rise.png')
  plt.close(fig)

def plot_rads():
  fig, ax = plt.subplots(2,4, figsize=(4/3*FIGX, FIGY))
  for avg in avgs:
    if 'beta_r' in avg:
      avg['betainv_r'] = 1/avg['beta_r']
    
    avg['Tp_r'] = avg['Pg_r'] / avg['rho_r']

  plot_multi(ax[0,0], 'r', 'rho_r', r"$<\rho>$", logy=True) #, ylim=[1.e-2, 1.e0])
  plot_multi(ax[0,1], 'r', 'Pg_r', r"$<P_g>$", logy=True) #, ylim=[1.e-6, 1.e-2])
  plot_multi(ax[0,2], 'r', 'Ptot_r', r"$<P_{tot}>$", logy=True) #, ylim=[1.e-6, 1.e-2])
  plot_multi(ax[0,3], 'r', 'B_r', r"$<|B|>$", logy=True) #, ylim=[1.e-4, 1.e-1])
  plot_multi(ax[1,0], 'r', 'u^phi_r', r"$<u^{\phi}>$", logy=True) #, ylim=[1.e-3, 1.e1])
  plot_multi(ax[1,1], 'r', 'u_phi_r', r"$<u_{\phi}>$", logy=True) #, ylim=[1.e-3, 1.e1])
  plot_multi(ax[1,2], 'r', 'Tp_r', r"$<T_{i}>$", logy=True) #, ylim=[1.e-6, 1.e-2])
  plot_multi(ax[1,3], 'r', 'betainv_r', r"$<\beta^{-1}>$", logy=True) #, ylim=[1.e-2, 1.e1])

  ax[0,3].legend(loc='upper right')

  pad = 0.05
  plt.subplots_adjust(left=pad, right=1-pad/2, bottom=pad, top=1-pad)
  plt.savefig(fname_out + '_ravgs.png')
  plt.close(fig)

def plot_fluxes():
  fig,ax = plt.subplots(4,1, figsize=(FIGX, 5*PLOTY))

  plot_multi(ax[0],'t','Mdot', r"$|\dot{M}|$")
  plot_multi(ax[1],'t','Phi_b', r"$\Phi$")

  for avg in avgs:
    if 'Ldot' in avg.keys():
      avg['abs_Ldot'] = np.abs(avg['Ldot'])
  plot_multi(ax[2],'t','abs_Ldot', r"$|\dot{L}|$")

  for avg in avgs:
    if 'Edot' in avg.keys() and 'Mdot' in avg.keys():
      avg['EmM'] = np.fabs(np.fabs(avg['Edot']) - np.fabs(avg['Mdot']))
  plot_multi(ax[3], 't', 'EmM', r"$|\dot{E} - \dot{M}|$")

#  plot_multi(ax[4], 't', 'Lum', "Lum", timelabels=True)

  ax[0].legend(loc='upper left')

  plt.savefig(fname_out + '_fluxes.png')
  plt.close(fig)

  # Don't print the norm fluxes if you can't norm
  if 'Mdot' not in avg.keys():
    return

  fig, ax = plt.subplots(4,1, figsize=(FIGX, 5*PLOTY))
  plot_multi(ax[0], 't', 'Mdot', r"$|\dot{M}|$")

  for avg in avgs:
    if 'Phi_b' in avg.keys():
      avg['phi_b'] = avg['Phi_b']/np.sqrt(np.fabs(avg['Mdot']))
  plot_multi(ax[1], 't', 'phi_b', r"$\frac{\Phi}{\sqrt{|\dot{M}|}}$")

  for avg in avgs:
    if 'Ldot' in avg.keys():
      avg['ldot'] = np.fabs(avg['Ldot'])/np.fabs(avg['Mdot'])
  plot_multi(ax[2], 't', 'ldot', r"$\frac{|Ldot|}{|\dot{M}|}$")

  for avg in avgs:
    if 'Edot' in avg.keys():
      avg['edot'] = np.fabs(-avg['Edot'] - avg['Mdot'])/(avg['Mdot'])
  plot_multi(ax[3], 't', 'edot', r"$\frac{|-\dot{E} - \dot{M}|}{|\dot{M}|}$")

#  for avg in avgs:
#    if 'Lum' in avg.keys():
#      avg['lum'] = np.fabs(avg['Lum'])/np.fabs(avg['Mdot'])
#  plot_multi(ax[4], 't', 'lum', r"$\frac{Lum}{|\dot{M}|}$", timelabels=True)
  
  ax[0].legend(loc='upper left')

  plt.savefig(fname_out + '_normfluxes.png')
  plt.close(fig)

def plot_extras():
  fig, ax = plt.subplots(6,1, figsize=(FIGX, 6*PLOTY))

  # Efficiency over mean mdot: just use the back half as <>
  for avg in avgs:
    if 'Edot' in avg.keys():
      mdot_mean = np.abs(np.mean(avg['Mdot'][len(avg['Mdot'])//2:]))
      avg['Eff'] = np.abs(avg['Mdot'] + avg['Edot'])/mdot_mean*100
  plot_multi(ax[0], 't', 'Eff', "Efficiency (%)", ylim=[-10,200])

  plot_multi(ax[1], 't', 'Edot', r"$\dot{E}$")

  for avg in avgs:
    if 'LBZ_bg1' in avg.keys():
      avg['aLBZ'] = np.abs(avg['LBZ_bg1'])
  plot_multi(ax[2], 't', 'aLBZ', "BZ Luminosity", timelabels=True)

  for avg in avgs:
    if 'aLBZ' in avg.keys() and 'Mdot' in avg.keys():
      avg['alBZ'] = avg['aLBZ']/np.fabs(avg['Mdot'])
  plot_multi(ax[3], 't', 'alBZ', r"$\frac{L_{BZ}}{\dot{M}}$", timelabels=True, ylim=[0, 4])

  for avg in avgs:
    if 'Lj_bg1' in avg.keys():
      avg['aLj'] = np.abs(avg['Lj_bg1'])
  plot_multi(ax[4], 't', 'aLj', "Jet Luminosity", timelabels=True)

  for avg in avgs:
    if 'aLBZ' in avg.keys() and 'Mdot' in avg.keys():
      avg['alj'] = avg['aLj']/np.fabs(avg['Mdot'])
  plot_multi(ax[5], 't', 'alj', r"$\frac{L_{jet}}{\dot{M}}$", timelabels=True, ylim=[0, 4])

  for i in range(len(avgs)):
    if 'alj' in avgs[i]:
      l_to_avg = avgs[i]['alj'][np.where(avgs[i]['t'] > 6000)]
      print("{} average (normalized) jet power: {}, std dev {}".format(labels[i], np.mean(l_to_avg), np.std(l_to_avg)))
#  for i in range(len(avgs)):
#    if 'alBZ' in avgs[i]:
#      print("{} average (normalized) BZ power: {}".format(labels[i], np.mean(avgs[i]['alBZ'][np.where(avgs[i]['t'] > 6000)])))

  ax[0].legend(loc='upper left')

  plt.savefig(fname_out + '_extras.png')
  plt.close(fig)

def plot_diags():
  fig, ax = plt.subplots(3,1, figsize=(FIGX, FIGY/2))
  ax = ax.flatten()

  plot_multi(ax[0], 't', 'Etot', "Total E")
  plot_multi(ax[1], 't', 'sigma_max', r"$\sigma_{max}$")
  # TODO include HARM's own diagnostics somehow? Re-insert just this one?
  plot_multi(ax[2], 't_d', 'divbmax_d', "max divB", timelabels=True)
  
  ax[0].legend(loc='lower left')

  plt.savefig(fname_out + '_diagnostics.png')
  plt.close(fig)

def plot_omega():
  # Omega
  fig, ax = plt.subplots(2,1, figsize=(FIGX, FIGY))
  # Renormalize omega to omega/Omega_H for plotting
  for avg in avgs:
    if 'omega_hth' in avg.keys(): #Then both are
      avg['omega_hth'] *= 4/avg['a']
      avg['omega_av_hth'] *= 4/avg['a']
  plot_multi(ax[0], 'th_5', 'omega_hth', r"$\omega_f$/$\Omega_H$ (EH, single shell)", ylim=[-1,2])
  plot_multi(ax[1], 'th_5', 'omega_av_hth', r"$\omega_f$/$\Omega_H$ (EH, 5-zone average)", ylim=[-1,2])

  # Legend
  ax[0].legend(loc='lower left')

  # Horizontal guidelines
  for a in ax.flatten():
    a.axhline(0.5, linestyle='--', color='k')

  plt.savefig(fname_out + '_omega.png')
  plt.close(fig)

def plot_flux_profs():
  # For converting to theta
  Xgeom = np.zeros((4,1,geom['n2']))
  Xgeom[1] = avg['r'][iBZ]
  Xgeom[2] = avg['th_100']
  to_th = 1/dxdX_to_KS(Xgeom, Met.FMKS, geom)[2,2,1]

  for avg in avgs:
    if 'FE_100_th' in avg.keys(): # Then all are
      avg['FE_100_th'] *= to_th
      avg['FE_Fl_100_th'] *= to_th
      avg['FE_EM_100_th'] *= to_th

  plot_multi(ax[0,0], 'th_100', 'FE_100_th', r"$\frac{dFE}{d\theta}$ ($r = 100$)")
  plot_multi(ax[0,1], 'th_100', 'FE_Fl_100_th', r"$\frac{dFE_{Fl}}{d\theta}$ ($r = 100$)")
  plot_multi(ax[1,0], 'th_100', 'FE_EM_100_th', r"$\frac{dFE_{EM}}{d\theta}$ ($r = 100$)")

  # Legend
  ax[0,0].legend(loc='lower left')

  plt.savefig(fname_out + '_flux_profs.png')
  plt.close(fig)

if __name__ == "__main__":
  if len(sys.argv) < 3:
    util.warn('Format: python eht_plot.py analysis_output [analysis_output ...] [labels_list] image_name')
    sys.exit()

  # All interior arguments are files to overplot, except possibly the last
  if len(sys.argv) < 4:
    last_file = -1
  else:
    last_file = -2


  avgs = []
  for filename in sys.argv[1:last_file]:
    # Encoding arg is for python2 numpy bytestrings
    avgs.append(pickle.load(open(filename,'rb'), encoding = 'latin1'))

  # Split the labels, or use the filename as a label
  # TODO something about filenames if this isn't present...
  labels = sys.argv[-2].split(",")

  if len(labels) < len(avgs):
    util.warn("Too few labels!")
    sys.exit()

  fname_out = sys.argv[-1]

  # For time plots.  Also take MAD/SANE for axes?
  ti = avgs[0]['t'][0]
  tf = avgs[0]['t'][-1]

  # Default styles
  if len(avgs) > 1:
    styles = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
  else:
    styles = ['k']
    
  if RADS: plot_rads()
  if FLUXES: plot_fluxes()
  if EXTRAS: plot_extras()
  if DIAGS: plot_diags()
  if OMEGA: plot_omega()
  if BSQ: plot_bsq_rise()
  if TEMP: plot_temp()
  if len(avgs) == 1:
    if FLUX_PROF: plot_flux_profs()


