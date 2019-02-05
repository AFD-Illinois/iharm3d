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
ax.set_title(r"$\log_{10}( -{{T_{EM}}^r}_t )$")

bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')

bplt.overlay_contours(ax, geom, dump['sigma'], [1.0], color='C2')
bplt.overlay_contours(ax, geom, dump['be_nob'], [0.02], color='C3')
bplt.overlay_contours(ax, geom, dump['be_nob'], [1.0], color='C4')

ax = plt.subplot(gs[1,0])
bplt.plot_xz(ax, geom, np.log10(-dump['Trt']-dump['RHO']*dump['ucon'][:,:,:,1]), arrayspace=USEARRSPACE, average=True, window=window)
ax.set_title(r"$\log_{10}( -{T^r}_t - \rho u^r )$")

bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')

bplt.overlay_contours(ax, geom, dump['sigma'], [1.0], color='C2')
bplt.overlay_contours(ax, geom, dump['be_nob'], [0.02], color='C3')
bplt.overlay_contours(ax, geom, dump['be_nob'], [1.0], color='C4')
#bplt.overlay_contours(ax, geom, dump['be_nob'], [1.0], color='C5')

ND = avg['LBZ_sigma1_rt'].shape[0]
rtype,rspin,rname = run_name.split("/")
start, end = [avg['avg_start'], avg['avg_end']]
# I can rely on this for now
start = int(start)//5
end = int(end)//5

# Average over quiescence
mdav = np.mean(np.abs(avg['mdot'][start:end]))
ab_av = lambda var : np.mean(np.abs(var[start:end,:]), axis=0)/mdav

ax = plt.subplot(gs[0,1])
ax.plot(avg['r'], ab_av(avg['LBZ_sigma1_rt']), label=r"$L_{BZ}$ (sigma > 1 cut)", color='C2')
ax.plot(avg['r'], ab_av(avg['LBZ_Be_nob0_rt']), label=r"$L_{BZ}$ ($Be > 0.02$ cut)", color='C3')
ax.plot(avg['r'], ab_av(avg['LBZ_Be_nob1_rt']), label=r"$L_{BZ}$ ($Be > 1.0$ cut)", color='C4')
ax.plot(avg['r'], ab_av(avg['LBZ_bg1_rt']), label=r"$L_{BZ}$ ($\beta\gamma > 1.0$ cut)", color='C5')

ax.set_title(r"$L_{BZ} = \int -{{T_{EM}}^r}_t \sqrt{-g} dx^{\theta} dx^{\phi}$")
ax.set_xlim([0,SIZE])
ax.set_xlabel("$r$ (M)")
ax.axvline(AT_R, color='k')

#maxes = [np.max(ab_av(avg['LBZ_'+tag+'_r'])[hdr['n1']//4:]) for tag in ['sigma1', 'be_nob1', 'be_nob0']]
#mins = [np.min(ab_av(avg['LBZ_'+tag+'_r'])[hdr['n1']//4:]) for tag in ['sigma1', 'be_nob1', 'be_nob0']]
#yhi = max(maxes); ylow = max(min(mins),1e-4*yhi)
#print(yhi, ylow)
#ax.set_ylim([ylow ,yhi])
if "SANE" in run_name:
  ax.set_yscale('log')

ax.legend(loc='upper right')

ax = plt.subplot(gs[1,1])
ax.plot(avg['r'], ab_av(avg['Lj_sigma1_rt']), label=r"$L_{tot}$ (sigma > 1 cut)", color='C2')
ax.plot(avg['r'], ab_av(avg['Lj_Be_nob0_rt']), label=r"$L_{tot}$ ($Be > 0.02$ cut)", color='C3')
ax.plot(avg['r'], ab_av(avg['Lj_Be_nob1_rt']), label=r"$L_{tot}$ ($Be > 1.0$ cut)", color='C4')
ax.plot(avg['r'], ab_av(avg['Lj_bg_rt']), label=r"$L_{tot}$ ($\beta\gamma > 1.0$ cut)", color='C5')

ax.set_title(r"$L_{tot} = \int (-{T^r}_t - \rho u^r) \sqrt{-g} dx^{\theta} dx^{\phi}$")
ax.set_xlim([0,SIZE])
ax.set_xlabel("$r$ (M)")
ax.axvline(AT_R, color='k')

#maxes = [np.max(ab_av(avg['Ltot_'+tag+'_r'])[hdr['n1']//4:]) for tag in  ['sigma1', 'be_nob1', 'be_nob0']]
#mins = [np.min(ab_av(avg['Ltot_'+tag+'_r'])[hdr['n1']//4:]) for tag in  ['sigma1', 'be_nob1', 'be_nob0']]
#yhi = max(maxes); ylow = max(min(mins),1e-4*yhi)
#print(yhi, ylow)
#ax.set_ylim([ylow,yhi])
if "SANE" in run_name:
  ax.set_yscale('log')

ax.legend(loc='lower right')

plt.tight_layout()
plt.savefig(run_name.replace("/","_")+"_L_study.png", dpi=100)
plt.close(fig)
