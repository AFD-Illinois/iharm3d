################################################################################
#                                                                              #
# BONDI INFLOW CONVERGENCE PLOTS                                               #
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True
from subprocess import call
from shutil import copyfile
import glob
import numpy as np
#import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl

sys.path.insert(0, '../../../../analysis-scripts/')
import util
import hdf5_to_dict as io

AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

RES = [32, 64, 128, 256]

NVAR = 8

L1 = np.zeros(len(RES))

# RUN PROBLEM FOR EACH RESOLUTION AND ANALYZE RESULT
for m in xrange(len(RES)):
  os.chdir('../dumps_' + str(RES[m]))

  dfiles = np.sort(glob.glob('dump*.h5'))
  hdr = io.load_hdr(dfiles[0])
  geom = io.load_geom(hdr, "grid.h5")
  dump0 = io.load_dump(dfiles[0], geom, hdr)
  dump1 = io.load_dump(dfiles[-1], geom, hdr)
  
  r = geom['r'][:,hdr['n2']/2,0]
  
#   print "Reh is ", str(hdr['Reh'])

  imin = 0
  while r[imin] < hdr['Reh']:
    imin += 1
  rho0 = np.mean(dump0['RHO'][imin:,:,0], axis=1)
  rho1 = np.mean(dump1['RHO'][imin:,:,0], axis=1)

  L1[m] = np.mean(np.fabs(rho1 - rho0))
  
  # PLOT
  fig = plt.figure(figsize=(16.18,10))

  ax = fig.add_subplot(1,1,1)
  ax.plot(r[imin:], rho0, marker='s', label='Initial')
  ax.plot(r[imin:], rho1, marker='s', label='Final')
  plt.xlabel('r (M)'); plt.ylabel('RHO')
  plt.title("Radial profile")
  plt.legend(loc=1)
  plt.savefig('../plots/bondi_comp_' + str(RES[m]) + '.png', bbox_inches='tight')

# MEASURE CONVERGENCE
powerfit = np.polyfit(np.log(RES), np.log(L1), 1)[0]
print "Powerfit: ", powerfit, "L1: ", L1

os.chdir('../plots/')

if not AUTO:
  # MAKE PLOTS
  fig = plt.figure(figsize=(16.18,10))

  ax = fig.add_subplot(1,1,1)
  ax.plot(RES, L1, marker='s', label='RHO')

  amp = 1.0e-3
  ax.plot([RES[0]/2., RES[-1]*2.], 
    10.*amp*np.asarray([RES[0]/2., RES[-1]*2.])**-2.,
    color='k', linestyle='--', label='N^-2')
  plt.xscale('log', basex=2); plt.yscale('log')
  plt.xlim([RES[0]/np.sqrt(2.), RES[-1]*np.sqrt(2.)])
  plt.xlabel('N'); plt.ylabel('L1')
  plt.title("BONDI")
  plt.legend(loc=1)
  plt.savefig('bondi.png', bbox_inches='tight')

if AUTO:
  data = {}
  #data['SOL'] = -2.*np.zeros([len(MODES), NVAR])  
  data['CODE'] = powerfit
  import pickle
  pickle.dump(data, open('data.p', 'wb'))

