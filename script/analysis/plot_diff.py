################################################################################
#                                                                              #
#  PLOT DIFFERENCES IN TWO FILES                                               #
#                                                                              #
################################################################################

import matplotlib
matplotlib.use('Agg')

import sys; sys.dont_write_bytecode = True
import numpy as np
import hdf5_to_dict as io
import matplotlib.pyplot as plt
import util
import glob
import os
import plot as bplt

USEARRSPACE=True
NLINES = 20
SIZE = 600

FIGX = 16
FIGY = 16

dump1file = sys.argv[1]
dump2file = sys.argv[2]
gridfile = sys.argv[3]

hdr = io.load_hdr(dump1file)
geom = io.load_geom(gridfile)
dump1 = io.load_dump(dump1file, geom, hdr)
dump2 = io.load_dump(dump2file, geom, hdr) #Hopefully this fails for dumps that shouldn't be compared

imname = 'differences.png'

fig = plt.figure(figsize=(FIGX, FIGY))

N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3']

def plot_diff(var, rel=False, lim=None):
    if rel:
        if lim is not None:
            bplt.plot_xy(ax, geom, np.abs((dump1[var] - dump2[var])/dump1[var]), dump1, vmin=0, vmax=lim, label=var, arrayspace=USEARRSPACE)
        else:
            bplt.plot_xy(ax, geom, np.abs((dump1[var] - dump2[var])/dump1[var]), dump1, label=var, arrayspace=USEARRSPACE)
    else:
        if lim is not None:
            bplt.plot_xy(ax, geom, np.log10(np.abs(dump1[var] - dump2[var])), dump1, vmin=0, vmax=lim, label=var, arrayspace=USEARRSPACE)
        else:
            bplt.plot_xy(ax, geom, np.log10(np.abs(dump1[var] - dump2[var])), dump1, vmin=-10, vmax=0, label=var, arrayspace=USEARRSPACE)

# Plot the difference
ax = plt.subplot(2,2,1)
plot_diff('RHO')#, lim=10e-3)
ax = plt.subplot(2,2,2)
plot_diff('UU')#, lim=1.1*10e-5)
ax = plt.subplot(2,2,3)
plot_diff('beta', lim=100)
ax = plt.subplot(2,2,4)
plot_diff('bsq')#, lim=5*10e-9)

#plt.subplots_adjust(left=0.1, right=0.9, bottom=0.05, top=0.95) # Avoid crowding
plt.tight_layout()

plt.savefig(imname, dpi=100) #, bbox_inches='tight')
plt.close(fig)
