################################################################################
#                                                                              #
# ONE-ZONE OPTICALLY THIN BREMSSTRAHLUNG COOLING                               #
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True
from subprocess import call
import glob
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
sys.path.insert(0, '../../analysis/')
import hdf5_to_dict as io
import units_cgs
cgs = units_cgs.get_dict()

TMP_DIR = 'TMP'
#util.safe_remove(TMP_DIR)

os.chdir('../prob/brem/')

print os.getcwd()

# COMPILE CODE
call(['python', 'build_mpi.py'])
os.chdir('../../test/auto/')
call(['mv', '../../prob/brem/bhlight', TMP_DIR])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
call(['mpirun','-np','8','./bhlight'])
os.chdir('../')

# READ SIMULATION OUTPUT
dfiles = np.sort(glob.glob(os.path.join(TMP_DIR,'')+'/dumps/dump*.h5'))
Nd = len(dfiles)
hdr = io.load_hdr(dfiles[0])
t_code = np.zeros(Nd)
Te_code = np.zeros(Nd)
for n in xrange(Nd):
  dump = io.load_dump(dfiles[n])
  t_code[n] = dump['t']*hdr['T_unit']
  
  Te_avg = np.mean(dump['Thetae']*cgs['ME']*cgs['CL']**2/cgs['KBOL'])
  Te_code[n] = Te_avg

# GET ANALYTIC SOLUTION
tf = 1.e8
Te0 = 1.e8
N  = 5.4e-39 # cm^3 K^1/2 s^-1 Sr^-1 Hz^-1
t_sol = np.linspace(0, tf, 1024)
Te_sol = Te0*(1. - t_sol/tf)**2.

# MAKE FIGURE
code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(1,1,1)
ax.plot(t_code, Te_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(t_sol, Te_sol, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('t (s)'); plt.ylabel('Te (K)')
plt.xlim([0,256*hdr['T_unit']]); plt.ylim([0, 1.1e8])

plt.savefig('brem_mpi.png', bbox_inches='tight')

# CLEAN UP
call(['rm', '-rf', TMP_DIR])

