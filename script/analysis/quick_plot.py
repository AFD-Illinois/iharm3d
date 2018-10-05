################################################################################
#                                                                              #
#  PLOT ONE PRIMITIVE                                                          #
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

FIGX = 20
FIGY = 20

dumpfile = sys.argv[1]
gridfile = sys.argv[2]
var = sys.argv[3]

hdr = io.load_hdr(dumpfile)
geom = io.load_geom(gridfile)
dump = io.load_dump(dumpfile, geom, hdr)

N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3']

# TODO properly handle vectors
# Plot the variable

nplotsx = 2
nplotsy = 2

fig = plt.figure(figsize=(FIGX, FIGY))
for n in range(4):
  ax = plt.subplot(nplotsy, nplotsx, n+1)
  #bplt.plot_xy(ax, var)
  bplt.plot_xy(ax, geom, np.log10(np.abs(dump['jcon'][:,:,:,n])), dump, arrayspace=USEARRSPACE)

plt.tight_layout()

plt.savefig(var+"_xy.png", dpi=100)
plt.close(fig)

fig = plt.figure(figsize=(FIGX, FIGY))
for n in range(4):
  ax = plt.subplot(nplotsx, nplotsy, n+1)
  #bplt.plot_xz(ax, var)
  bplt.plot_xz(ax, geom, np.log10(np.abs(dump['jcon'][:,:,:,n])), dump, arrayspace=USEARRSPACE)

plt.tight_layout()

plt.savefig(var+"_xz.png", dpi=100)
plt.close(fig)
