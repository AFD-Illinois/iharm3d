################################################################################
#                                                                              #
#  PLOT ONE PRIMITIVE                                                          #
#                                                                              #
################################################################################

import hdf5_to_dict as io
import plot as bplt
from analysis_fns import *
import units

import matplotlib
import matplotlib.pyplot as plt

import sys
import numpy as np

# TODO parse lots of options I set here
USEARRSPACE=False
UNITS=False

SIZE = 100
window=[-SIZE,SIZE,-SIZE,SIZE]
FIGX = 10
FIGY = 10

dumpfile = sys.argv[1]
if len(sys.argv) > 3:
  gridfile = sys.argv[2]
  var = sys.argv[3]
elif len(sys.argv) > 2:
  gridfile = None
  var = sys.argv[2]

if UNITS and var not in ['Tp']:
  M_unit = float(sys.argv[-1])

if gridfile is not None:
  hdr = io.load_hdr(dumpfile)
  geom = io.load_geom(hdr, gridfile)
  dump = io.load_dump(dumpfile, hdr, geom)
else:
  # Assumes gridfile in same directory
  hdr,geom,dump = io.load_all(dumpfile)

fig = plt.figure(figsize=(FIGX, FIGY))

# If we're plotting a derived variable, calculate + add it
if var in ['jcov', 'jsq']:
  dump['jcov'] = np.einsum("...i,...ij->...j", dump['jcon'], geom['gcov'][:,:,None,:,n])
  for n in range(hdr['n_dim']):
    dump['jcov'][:,:,:,n] = np.sum(dump['jcon']*geom['gcov'][:,:,None,:,n], axis=3)
  dump['jsq'] = np.sum(dump['jcon']*dump['jcov'], axis=-1)
elif var not in dump:
  dump[var] = d_fns[var](dump)


# Add units after all calculations, manually
if UNITS:
  if var in ['Tp']:
    cgs = units.get_cgs()
    dump[var] /= cgs['KBOL']
  else:
    unit = units.get_units_M87(M_unit, tp_over_te=3)

    if var in ['bsq']:
      dump[var] *= unit['B_unit']**2
    elif var in ['B']:
      dump[var] *= unit['B_unit']
    elif var in ['Ne']:
      dump[var] = dump['RHO'] * unit['Ne_unit']
    elif var in ['Te']:
      dump[var] = ref['ME'] * ref['CL']**2 * unit['Thetae_unit'] * dump['UU']/dump['RHO']
    elif var in ['Thetae']:
      # TODO non-const te
      dump[var] = unit['Thetae_unit'] * dump['UU']/dump['RHO']

# Plot XY differently for vectors, scalars
if var in ['jcon','ucon','ucov','bcon','bcov']:
  axes = [plt.subplot(2, 2, i) for i in range(1,5)]
  for n in range(4):
    bplt.plot_xy(axes[n], geom, np.log10(dump[var][:,:,:,n]), arrayspace=USEARRSPACE, window=window)
else:
  # TODO allow specifying vmin/max, average from command line or above
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xy(ax, geom, dump[var], arrayspace=USEARRSPACE, window=window, vmin=1e10, vmax=1e12)

plt.tight_layout()

plt.savefig(var+"_xy.png", dpi=100)
plt.close(fig)

fig = plt.figure(figsize=(FIGX, FIGY))

# Plot XZ
if var in ['jcon','ucon','ucov','bcon','bcov']:
  axes = [plt.subplot(2, 2, i) for i in range(1,5)]
  for n in range(4):
    bplt.plot_xz(axes[n], geom, np.log10(dump[var][:,:,:,n]), arrayspace=USEARRSPACE, window=window)
else:
  ax = plt.subplot(1, 1, 1)
  bplt.plot_xz(ax, geom, np.log10(dump[var]), arrayspace=USEARRSPACE, window=window, vmin=-3, vmax=2)
  #JE1 = -T_mixed(dump, 1,0)
  #JE2 = T_mixed(dump, 2,0)
  #JE1 = dump['ucon'][:,:,:,1]
  #JE2 = dump['ucon'][:,:,:,2]
  #bplt.overlay_flowlines(ax, geom, JE1, JE2, nlines=1000, arrayspace=USEARRSPACE)
  #bplt.overlay_eflow_quiver(ax, geom, dump)
  bplt.overlay_field(ax, geom, dump, nlines=30, arrayspace=USEARRSPACE)

plt.tight_layout()

plt.savefig(var+"_xz.png", dpi=100)
plt.close(fig)
