################################################################################
#                                                                              #
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      #
#                                                                              #
################################################################################

import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../')
import numpy as np
import hdf5_to_dict as io
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import util
import glob
import os
import plot as bplt

FIGX = 15
FIGY = 15
NPLOTSX = 3 # TODO switch these names
NPLOTSY = 2
NLINES = 20
SIZE = 600

# TODO this is a dirty hack
args_bad = False
if sys.argv[1] == '-d':
    debug = True
    path = sys.argv[2]
    if len(sys.argv) != 3:
        args_bad = True
else:
    debug = False
    path = sys.argv[1]
    if len(sys.argv) != 2:
        args_bad = True

if args_bad:
  util.warn('PATH TO DUMP FOLDER NEEDED AS ARGUMENT')
  sys.exit(1)

files = np.sort(glob.glob(os.path.join(path, "dump*.h5")))
gridfile = os.path.join(path,"grid.h5")

if len(files) == 0:
    util.warn("INVALID PATH TO DUMP FOLDER")
    sys.exit(1)

hdr = io.load_hdr(files[0])
geom = io.load_geom(hdr, gridfile)
dump = io.load_dump(files[0], geom, hdr)

# Plot the first dump, specifically init as in Narayan '12

imname = 'initial_conditions.png'

fig = plt.figure(figsize=(FIGX, FIGY))

N1 = hdr['N1']; N2 = hdr['N2']; N3 = hdr['N3']

ax = plt.subplot(NPLOTSX,NPLOTSY,1)
bplt.plot_r(ax, geom, dump['RHO'], N2/2, N3/2, 'RHO', logx=True, logy=True)
ax.set_xlim([8, 2*10**3]); ax.set_ylim([10**(-4), 2])

flux = np.sum(dump['B2']*geom['gdet'][:,:,None]*hdr['dx1']*hdr['dx3'],axis=-1)

ax = plt.subplot(NPLOTSX,NPLOTSY,2)
bplt.plot_r(ax, geom, flux, N2/2, None, 'flux', logx=False, logy=False)
ax.set_xlim([0, SIZE]); ax.set_ylim([-200, 10])

ax = plt.subplot(NPLOTSX,NPLOTSY,3)
bplt.plot_xz(ax, geom, np.log10(dump['RHO']), dump,
             vmin=-4, vmax = 0, label='RHO')
ax.set_xlim([0, SIZE]); ax.set_ylim([-SIZE/2, SIZE/2])

ax = plt.subplot(NPLOTSX,NPLOTSY,4)
bplt.plot_xz(ax, geom, np.log10(dump['UU']), dump,
             vmin=-4, vmax = 0, label='UU')
ax.set_xlim([0, SIZE]); ax.set_ylim([-SIZE/2, SIZE/2])

# I change this a lot
var2_str = 'beta'
var2_data = np.log10(dump[var2_str])

ax = plt.subplot(NPLOTSX,NPLOTSY,5)
bplt.plot_xz(ax, geom, var2_data, dump,
             label=var2_str, cmap='RdBu_r', vmin=1, vmax=4)
bplt.overlay_field(ax, geom, dump, NLINES)
ax.set_xlim([0, SIZE]); ax.set_ylim([-SIZE/2, SIZE/2])

# I change this a lot
var2_str = 'bsq'
var2_data = np.log10(dump[var2_str])

ax = plt.subplot(NPLOTSX,NPLOTSY,6)
bplt.plot_xz(ax, geom, var2_data, dump,
             label=var2_str, cmap='RdBu_r', vmin=-8, vmax=2)
bplt.overlay_field(ax, geom, dump, NLINES)
ax.set_xlim([0, SIZE]); ax.set_ylim([-SIZE/2, SIZE/2])

# TODO enlarge plots w/o messing up even pixel count
# Maybe share axes, even?
plt.subplots_adjust(left=0.15, right=0.85, bottom=0.05, top=0.95) # Avoid crowding
#plt.tight_layout()

plt.savefig(imname, dpi=100) #, bbox_inches='tight')
plt.close(fig)
