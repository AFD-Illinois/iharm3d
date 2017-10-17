################################################################################
#                                                                              #
# ONE-ZONE OPTICALLY THIN BREMSSTRAHLUNG COOLING                               #
#                                                                              #
################################################################################

import os
import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../script/')
sys.path.insert(0, '../script/analysis')
from subprocess import call
import glob
import numpy as np
import hdf5_to_dict as io
import units
cgs = units.get_cgs()
import util

TMP_DIR = 'TMP'
util.safe_remove(TMP_DIR)
PROBLEM = 'brem'
AUTO = False
for arg in sys.argv:
  if arg == '-auto':
    AUTO = True

os.chdir('../prob/' + PROBLEM)

# COMPILE CODE
call(['python', 'build.py', '-dir', TMP_DIR])
os.chdir('../../test')
call(['mv', '../prob/' + PROBLEM + '/' + TMP_DIR, './'])

# RUN EXECUTABLE
os.chdir(TMP_DIR)
call(['./bhlight', '-p', 'param_template.dat'])
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
  Te_code[n] = dump['Thetae'][0][0][0]*cgs['ME']*cgs['CL']**2/cgs['KBOL']

# GET ANALYTIC SOLUTION
tf = 1.e8
Te0 = 1.e8
N  = 5.4e-39 # cm^3 K^1/2 s^-1 Sr^-1 Hz^-1
t_sol = np.linspace(0, tf, 1024)
Te_sol = Te0*(1. - t_sol/tf)**2.

if AUTO:
  data = {}
  data['SOL'] = [t_sol, Te_sol]
  data['CODE'] = [t_code[:3./4.*len(t_code)], Te_code[:3./4.*len(Te_code)]]
  import pickle
  pickle.dump(data, open('data.p', 'wb'))
  sys.exit()

# MAKE FIGURE
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab as pl
code_col = 'r'; code_ls = ''; code_mrk = '.'
sol_col = 'k'; sol_ls = '-'; sol_mrk = ''
fig = plt.figure(figsize=(16.18,10))

ax = fig.add_subplot(1,1,1)
ax.plot(t_code, Te_code, color=code_col, linestyle=code_ls, marker=code_mrk)
ax.plot(t_sol, Te_sol, color=sol_col, linestyle=sol_ls, marker=sol_mrk)
plt.xlabel('t (s)'); plt.ylabel('Te (K)')
plt.xlim([0,256*hdr['T_unit']]); plt.ylim([0, 1.1e8])

plt.savefig('brem.png', bbox_inches='tight')

# CLEAN UP
util.safe_remove(TMP_DIR)

