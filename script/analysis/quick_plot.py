################################################################################
#                                                                              #
#  PLOT ONE PRIMITIVE                                                          #
#                                                                              #
################################################################################

import sys

import numpy as np
import matplotlib.pyplot as plt

import units
import util
import hdf5_to_dict as io
import plot as bplt
from analysis_fns import *

USEARRSPACE=False
UNITS=True
SIZE = 40

FIGX = 20
FIGY = 20

dumpfile = sys.argv[1]
gridfile = sys.argv[2]
var = sys.argv[3]
if UNITS:
  M_unit = float(sys.argv[4])

hdr = io.load_hdr(dumpfile)
geom = io.load_geom(hdr, gridfile)
dump = io.load_dump(dumpfile, hdr, geom)

init_analysis(hdr, geom)
bplt.init_plotting(hdr, geom)

N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3']

nplotsx = 2
nplotsy = 2

fig = plt.figure(figsize=(FIGX, FIGY))

if var == 'jcov':
  dump['jcov'] = np.zeros_like(dump['jcon'])
  for n in range(hdr['n_dim']):
    dump['jcov'][:,:,:,n] = np.sum(dump['jcon']*geom['gcov'][:,:,None,:,n], axis=3)
  dump['jsq'] = np.sum(dump['jcon']*dump['jcov'], axis=-1)
elif var == 'sigma':
  dump['sigma'] = dump['bsq']/dump['RHO']
elif var == 'bernoulli':
  dump['bernoulli'] = -Tmixed(dump,0,0) /(dump['RHO']*dump['ucon'][:,:,:,0]) - 1
elif var == 'B':
  dump['B'] = np.sqrt(dump['bsq'])

# Add units after all calculations, manually
# M87
if UNITS:
  L_unit = 9.15766e+14
  T_unit = 30546.6

  ref = units.get_cgs()
  RHO_unit = M_unit / (L_unit ** 3)
  U_unit = RHO_unit*ref['CL']**2;
  B_unit = ref['CL']*np.sqrt(4. * np.pi * RHO_unit)
  Ne_unit = RHO_unit/(ref['MP'] + ref['ME'])
  # TODO option for const
  tp_over_te = 3
  Thetae_unit = (hdr['gam']-1.)*ref['MP']/ref['ME']/(1. + tp_over_te)

  if var in ['bsq']:
    dump[var] *= B_unit**2
  elif var in ['B']:
    dump[var] *= B_unit
  elif var in ['Ne']: # Meaningless w/o units
    dump[var] = dump['RHO'] * Ne_unit
  elif var in ['Te']:
    dump[var] = ref['ME'] * ref['CL']**2 * Thetae_unit * dump['UU']/dump['RHO']
  elif var in ['Thetae']:
    dump[var] = Thetae_unit * dump['UU']/dump['RHO']
  # TODO the others


# Plot XY differently for vectors, scalars
if var in ['jcon','ucon','ucov','bcon','bcov']:
  for n in range(4):
    ax = plt.subplot(nplotsy, nplotsx, n+1)
    bplt.plot_xy(ax, geom, np.log10(np.abs(dump[var][:,:,:,n])), dump, arrayspace=USEARRSPACE)
elif var in ['sigma']:
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xy(ax, dump[var], vmin=0, vmax=10, arrayspace=USEARRSPACE)
else:
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xy(ax, dump[var], arrayspace=USEARRSPACE) 

plt.tight_layout()

plt.savefig(var+"_xy.png", dpi=100)
plt.close(fig)

fig = plt.figure(figsize=(FIGX, FIGY))

# Plot XZ
if var in ['jcon','ucon','ucov','bcon','bcov']:
  for n in range(4):
    ax = plt.subplot(nplotsx, nplotsy, n+1)
    bplt.plot_xz(ax, np.log10(np.abs(dump[var][:,:,:,n])), arrayspace=USEARRSPACE)
elif var in ['sigma']:
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xz(ax, dump[var], vmin=0, vmax=10, arrayspace=USEARRSPACE)
elif var in ['bernoulli']:
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xz(ax, dump[var], arrayspace=USEARRSPACE, average=True)
  bplt.overlay_contour(ax, dump[var], [0.05])
else:
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xz(ax, dump[var], arrayspace=USEARRSPACE)

plt.tight_layout()

plt.savefig(var+"_xz.png", dpi=100)
plt.close(fig)
