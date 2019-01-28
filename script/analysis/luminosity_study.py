################################################################################
#                                                                              #
#  LUMINOSITY COMPARISON                                                       #
#                                                                              #
################################################################################

import os, sys
import pickle
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import hdf5_to_dict as io
import plot as bplt
from analysis_fns import *


USEARRSPACE=False

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

dumpfile = os.path.join("/scratch/03002/bprather/pharm_dumps/M87SimulationLibrary/GRMHD",run_name,"dumps/dump_00001500.h5")
hdr,geom,dump = io.load_all(dumpfile)

plotfile = os.path.join("/work/03002/bprather/stampede2/movies",run_name,"eht_out.p")
avg = pickle.load(open(plotfile, "rb"))

dump['sigma'] = dump['bsq']/dump['RHO']
dump['gamma'] = get_gamma(geom, dump)
dump['be_b'] = bernoulli(dump, with_B=True)
dump['be_nob'] = bernoulli(dump, with_B=False)
dump['TEMrt'] = TEM_mixed(dump, 1, 0)
dump['Trt'] = T_mixed(dump, 1, 0)

fig = plt.figure(figsize=(FIGX, FIGY))
gs = gridspec.GridSpec(2, 2, width_ratios=[1,2])

ax = plt.subplot(gs[0,0])
bplt.plot_xz(ax, geom, np.log10(-dump['TEMrt']), arrayspace=USEARRSPACE, average=True, window=window)
ax.set_title(r"$-{{T_{EM}}^r}_t$")

bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')

#bplt.overlay_contours(ax, geom, dump['ucon'][:,:,:,1], [0.0], color='k')

bplt.overlay_contours(ax, geom, dump['sigma'], [1.0], color='C2')
#bplt.overlay_contours(ax, geom, dump['sigma'], [10.0], color='C3')

#bplt.overlay_contours(ax, geom, dump['be_b'], [0.02], color='C4')
#bplt.overlay_contours(ax, geom, dump['be_b'], [1.0], color='C5')
bplt.overlay_contours(ax, geom, dump['be_nob'], [0.02], color='C3')
bplt.overlay_contours(ax, geom, dump['be_nob'], [1.0], color='C4')

#bplt.overlay_contours(ax, geom, geom['r']*dump['ucon'][:,:,:,1], [1.0], color='C8')
#bplt.overlay_contours(ax, geom, dump['gamma'], [1.5], color='C9')


#bplt.overlay_contours(ax, geom, geom['r']*dump['ucon'][:,:,:,1], [1.0], color='C4')
#bplt.overlay_contours(ax, geom, dump['gamma'], [1.5], color='C5')

ax = plt.subplot(gs[1,0])
bplt.plot_xz(ax, geom, np.log10(-dump['Trt']-dump['RHO']*dump['ucon'][:,:,:,1]), arrayspace=USEARRSPACE, average=True, window=window)
ax.set_title(r"$-{T^r}_t - \rho u^r$")

bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')

#bplt.overlay_contours(ax, geom, dump['ucon'][:,:,:,1], [0.0], color='k')

bplt.overlay_contours(ax, geom, dump['sigma'], [1.0], color='C2')
#bplt.overlay_contours(ax, geom, dump['sigma'], [10.0], color='C3')

#bplt.overlay_contours(ax, geom, dump['be_b'], [0.02], color='C4')
#bplt.overlay_contours(ax, geom, dump['be_b'], [1.0], color='C5')
bplt.overlay_contours(ax, geom, dump['be_nob'], [0.02], color='C3')
bplt.overlay_contours(ax, geom, dump['be_nob'], [1.0], color='C4')

#bplt.overlay_contours(ax, geom, geom['r']*dump['ucon'][:,:,:,1], [1.0], color='C8')
#bplt.overlay_contours(ax, geom, dump['gamma'], [1.5], color='C9')

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
#start, end = qui[rtype][rname][rspin]
# I can rely on this for now
#start = int(start)//5
#end = int(end)//5

# Average over quiescence
#mdav = np.mean(np.abs(avg['mdot'][start:end]))
#ab_av = lambda var : np.mean(np.abs(var[start:end,:]), axis=0)/mdav

# For pure radial versions of variables
mdav = np.mean(np.abs(avg['mdot'][1600:]))
ab_av = lambda var : np.abs(var[:])/mdav

ax = plt.subplot(gs[0,1])
ax.plot(avg['r'], ab_av(avg['LBZ_sigma1_r']), label=r"$L_{BZ}$ (sigma > 1 cut)", color='C2')
#ax.plot(avg['r'], ab_av(avg['LBZ_sigma10_r']), label=r"$L_{BZ}$ (sigma > 10 cut)", color='C3')

#ax.plot(avg['r'], ab_av(avg['LBZ_be_b0_r']), label=r"$L_{BZ}$ ($Be_B > 0.02$ cut)", color='C4')
#ax.plot(avg['r'], ab_av(avg['LBZ_be_b1_r']), label=r"$L_{BZ}$ ($Be_B > 1.0$ cut)", color='C5')
ax.plot(avg['r'], ab_av(avg['LBZ_be_nob0_r']), label=r"$L_{BZ}$ ($Be > 0.02$ cut)", color='C3')
ax.plot(avg['r'], ab_av(avg['LBZ_be_nob1_r']), label=r"$L_{BZ}$ ($Be > 1.0$ cut)", color='C4')

#ax.plot(avg['r'], ab_av(avg['LBZ_rur_r']), label=r"$L_{BZ}$ (r*u^r cut)", color='C8')
#ax.plot(avg['r'], ab_av(avg['LBZ_gamma_r']), label=r"$L_{BZ}$ (gamma cut)", color='C9')

ax.set_title(r"$L_{BZ} = \int -{{T_{EM}}^r}_t \sqrt{-g} dx^{\theta} dx^{\phi}$")
ax.set_xlim([0,SIZE])
ax.axvline(AT_R, color='k')

maxes = [np.max(ab_av(avg['LBZ_'+tag+'_r'])[hdr['n1']//4:]) for tag in ['sigma1', 'be_nob1', 'be_nob0']]
mins = [np.min(ab_av(avg['LBZ_'+tag+'_r'])[hdr['n1']//4:]) for tag in ['sigma1', 'be_nob1', 'be_nob0']]
yhi = max(maxes); ylow = min(mins)
ax.set_ylim([ylow ,yhi])
if "SANE" in run_name:
  ax.set_yscale('log')

ax.legend(loc='upper right')

ax = plt.subplot(gs[1,1])
ax.plot(avg['r'], ab_av(avg['Ltot_sigma1_r']), label=r"$L_{tot}$ (sigma > 1 cut)", color='C2')
#ax.plot(avg['r'], ab_av(avg['Ltot_sigma10_r']), label=r"$L_{tot}$ (sigma > 10 cut)", color='C3')

#ax.plot(avg['r'], ab_av(avg['Ltot_be_b0_r']), label=r"$L_{tot}$ ($Be_B > 0.02$ cut)", color='C4')
#ax.plot(avg['r'], ab_av(avg['Ltot_be_b1_r']), label=r"$L_{tot}$ ($Be_B > 1.0$ cut)", color='C5')
ax.plot(avg['r'], ab_av(avg['Ltot_be_nob0_r']), label=r"$L_{tot}$ ($Be_{gas} > 0.02$ cut)", color='C3')
ax.plot(avg['r'], ab_av(avg['Ltot_be_nob1_r']), label=r"$L_{tot}$ ($Be_{gas} > 1.0$ cut)", color='C4')

#ax.plot(avg['r'], ab_av(avg['Ltot_rur_r']), label=r"$L_{tot}$ (r*u^r cut)", color='C8')
#ax.plot(avg['r'], ab_av(avg['Ltot_gamma_r']), label=r"$L_{tot}$ (gamma cut)", color='C9')

ax.set_title(r"$L_{tot} = \int (-{T^r}_t - \rho u^r) \sqrt{-g} dx^{\theta} dx^{\phi}$")
ax.set_xlim([0,SIZE])
ax.axvline(AT_R, color='k')

maxes = [np.max(ab_av(avg['Ltot_'+tag+'_r'])[hdr['n1']//4:]) for tag in  ['sigma1', 'be_nob1', 'be_nob0']]
mins = [np.min(ab_av(avg['Ltot_'+tag+'_r'])[hdr['n1']//4:]) for tag in  ['sigma1', 'be_nob1', 'be_nob0']]
yhi = max(maxes); ylow = min(mins)
ax.set_ylim([ylow,yhi])
if "SANE" in run_name:
  ax.set_yscale('log')

ax.legend(loc='lower right')

plt.tight_layout()
plt.savefig(run_name.replace("/","_")+"_L_study.png", dpi=100)
plt.close(fig)
