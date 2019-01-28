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
else:
  SIZE = 400

window=[0,SIZE/2,0,SIZE]
FIGX = 15
FIGY = 15

dumpdir = os.path.join("/scratch/03002/bprather/pharm_dumps/M87SimulationLibrary/GRMHD",run_name,"dumps")
hdr = io.load_hdr(os.path.join(dumpdir,"dump_00001500.h5"))
geom = io.load_geom(hdr, dumpdir)

plotfile = os.path.join("/work/03002/bprather/stampede2/movies",run_name,"eht_out.p")
avg = pickle.load(open(plotfile, "rb"))

fig, ax = plt.subplots(1,1,figsize=(FIGX, FIGY))

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

# Average over quiescence
mdav = np.mean(np.abs(avg['mdot'][start:end]))
ab_av = lambda var : np.mean(np.abs(var[start:end,:]), axis=0)/mdav

fig, ax = plt.subplots(1,1,figsize=(FIGX, FIGY))

ax.plot(avg['r'], ab_av(avg['LBZ_sigma1_rt']), label=r"$L_{BZ}$ (sigma > 1 cut)", color='C2')
#ax.plot(avg['r'], ab_av(avg['LBZ_be_b0_r']), label=r"$L_{BZ}$ ($Be_B > 0.02$ cut)", color='C4')
#ax.plot(avg['r'], ab_av(avg['LBZ_be_b1_r']), label=r"$L_{BZ}$ ($Be_B > 1.0$ cut)", color='C5')
#ax.plot(avg['r'], ab_av(avg['LBZ_be_nob0_r']), label=r"$L_{BZ}$ ($Be_{gas} > 0.02$ cut)", color='C6')
#ax.plot(avg['r'], ab_av(avg['LBZ_be_nob1_rt']), label=r"$L_{BZ}$ ($Be_{gas} > 1.0$ cut)", color='C7')

ax.plot(avg['r'], ab_av(avg['Ltot_sigma1_rt']), label=r"$L_{tot}$ (sigma > 1 cut)", color='C3')
#ax.plot(avg['r'], ab_av(avg['Ltot_be_b0_r']), label=r"$L_{tot}$ ($Be_B > 0.02$ cut)", color='C4')
#ax.plot(avg['r'], ab_av(avg['Ltot_be_b1_r']), label=r"$L_{tot}$ ($Be_B > 1.0$ cut)", color='C5')
#ax.plot(avg['r'], ab_av(avg['Ltot_be_nob0_r']), label=r"$L_{tot}$ ($Be_{gas} > 0.02$ cut)", color='C6')
#ax.plot(avg['r'], ab_av(avg['Ltot_be_nob1_rt']), label=r"$L_{tot}$ ($Be_{gas} > 1.0$ cut)", color='C8')

ax.axvline(geom['r_eh'], color='k')
print(geom['r_eh'])

ax.set_title(r"$L_{BZ}$ vs $L_{tot}$ comparison")
ax.set_xlim([0,100])
maxes = [np.max(ab_av(avg['LBZ_'+tag+'_rt'])[50:]) for tag in ['sigma1']]
mins = [np.min(ab_av(avg['LBZ_'+tag+'_rt'])[50:]) for tag in ['sigma1']]
yhi = max(maxes); ylow = min(mins)
ax.set_ylim([ylow ,yhi])
if "SANE" in run_name:
  ax.set_yscale('log')

ax.legend(loc='upper right')

plt.tight_layout()
plt.savefig(run_name.replace("/","_")+"_L_overplot.png", dpi=100)
plt.close(fig)
