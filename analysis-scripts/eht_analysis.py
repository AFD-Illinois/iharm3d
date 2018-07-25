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
import psutil
#import plot as bplt

FIGX = 14
FIGY = 10
#SIZE = 40
THMIN = np.pi/3.
THMAX = 2.*np.pi/3.

if len(sys.argv) != 3:
  util.warn('Format: python eht_analysis.py [dump path] [tavg]')
  sys.exit()

path = sys.argv[1]
tavg = float(sys.argv[2])

dumps = io.get_dumps_reduced(path)

hdr = io.load_hdr(dumps[0])
geom = io.load_geom(os.path.join(path,'grid.h5'))

# Limit threads for 256^3 problem due to memory
# TODO guess NTHREADS better (96e9/N^3*8?)
if hdr['n1'] >= 256:
  NTHREADS = 12
else:
  NTHREADS = psutil.cpu_count(logical=False)

# Calculate jmin, jmax
ths = geom['th'][-1,:,0]
for n in xrange(len(ths)):
  if ths[n] > THMIN:
    jmin = n
    break
for n in xrange(len(ths)):
  if ths[n] > THMAX:
    jmax = n
    break

diag = io.load_diag(hdr, path)

dx1 = hdr['dx1']
dx2 = hdr['dx2']
dx3 = hdr['dx3']
N1 = hdr['n1']
N2 = hdr['n2']
N3 = hdr['n3']

r = geom['r'][:,N2/2,0]
th = geom['th'][0,:N2/2,0]

# Option to calculate fluxes at (just inside) r = 5
# This reduces interference from floors
flux_out = False
iF = 0
if flux_out:
  while r[iF] < 5:
    iF += 1
  iF -= 1
else:
  iF = 5

# Some variables (Phi) should only be computed at EH
iEH = 5

ND = len(dumps)

avg_keys = ['rho_r', 'Theta_r', 'B_r', 'Pg_r', 'Ptot_r', 'betainv_r', 'uphi_r', 'rho_SADW', 'omega_th', 'omega_th_av', 'omega_th_5']


# EVALUATE DIAGNOSTICS
  
vol = (dx2*2.*np.pi*geom['gdet'][:,:]).sum(axis=-1)

def INT(var):
  return (dx2*dx3*geom['gdet'][:,:,None]*var[:,:,:]).sum(axis=-1).sum(axis=-1)

def WAVG(var, w):
  return INT(w*var)/INT(w)

def avg_dump(n):
  out = {}

  print '%08d / ' % (n+1) + '%08d' % len(dumps) 
  dump = io.load_dump(dumps[n], geom, hdr)
  out['t'] = dump['t']
  print dump['t']

  out['rho_SADW'] = WAVG(dump['RHO'], dump['RHO'])

  # SHELL AVERAGES
  #vol = (dx2*dx3*geom['gdet'][:,:]).sum(axis=-1)
  out['rho_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*dump['RHO'][:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)
  
  Theta = (hdr['gam']-1.)*dump['UU'][:,:,:]/dump['RHO'][:,:,:]
  out['Theta_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Theta[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)
  
  B = np.sqrt(dump['bsq'][:,:,:])
  out['B_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*B[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  Pg = (hdr['gam']-1.)*dump['UU'][:,:,:]
  out['Pg_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Pg[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  Ptot = Pg + dump['bsq'][:,:,:]/2
  out['Ptot_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Ptot[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  betainv = (dump['bsq'][:,:,:]/2)/Pg
  out['betainv_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*betainv[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  uphi = (dump['ucon'][:,:,:,3])
  out['uphi_r'] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*uphi[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1) 

  out['rho_r'] /= vol
  out['Theta_r'] /= vol
  out['B_r'] /= vol
  out['Pg_r'] /= vol
  out['Ptot_r'] /= vol
  out['betainv_r'] /= vol
  out['uphi_r'] /= vol

  # FLUXES
  #Phi_bcon[n] = 0.5*(np.fabs(dump['bcon'][iF,:,:,1])*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  # The HARM B_unit is sqrt(4pi)*c*sqrt(rho) which has caused issues:
  norm = np.sqrt(4*np.pi) # This is what I believe matches Narayan '12 in CGS
  #norm = 1 # This is what BR used for EHT comparison + what fits that
  out['Phi'] = 0.5*(np.fabs(norm*dump['B1'][iEH,:,:])*geom['gdet'][iEH,:,None]*dx2*dx3).sum()
  
  out['omega_th'] = dump['omega'][iEH,:N2/2,:].mean(axis=-1) + dump['omega'][iEH,:N2/2-1:-1,:].mean(axis=-1)
  out['omega_th_av'] = dump['omega'][iEH:iEH+5,:N2/2,:].mean(axis=-1).mean(axis=0)# + dump['omega'][iEH:iEH+5,:N2/2-1:-1,:].sum(axis=-1).sum(axis=0)
  #out['omega_th_av'] /= N3*2*5
  out['omega_th_5'] = dump['omega'][:10,:N2/2,:].mean(axis=-1).mean(axis=0)# + dump['omega'][:10,:N2/2-1:-1,:].sum(axis=-1).sum(axis=0)
  #out['omega_th_5'] /= N3*2*10

  out['sigma_max'] = np.max(dump['bsq']/dump['RHO'])
  
  out['Mdot'] = np.abs((dump['RHO'][iF,:,:]*dump['ucon'][iF,:,:,1]*geom['gdet'][iF,:,None]*dx2*dx3).sum())
  Trphi = (dump['RHO'] + dump['UU'] + (hdr['gam']-1.)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,1]*dump['ucov'][:,:,:,3]
  Trphi -= dump['bcon'][:,:,:,1]*dump['bcov'][:,:,:,3]
  Trt = (dump['RHO'] + dump['UU'] + (hdr['gam']-1.)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,1]*dump['ucov'][:,:,:,0]
  Trt -= dump['bcon'][:,:,:,1]*dump['bcov'][:,:,:,0]
  out['Ldot'] = (Trphi[iF,:,:]*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  out['Edot'] = -(Trt[iF,:,:]*geom['gdet'][iF,:,None]*dx2*dx3).sum()

  rho = dump['RHO']
  P = (hdr['gam']-1.)*dump['UU']
  B = np.sqrt(dump['bsq'])
  C = 0.2
  j = rho**3*P**(-2)*np.exp(-C*(rho**2/(B*P**2))**(1./3.))
  out['Lum'] = (j*geom['gdet'][:,:,None]*dx1*dx2*dx3).sum()
  
  return out

out_full = {}

# SERIAL
#for n in xrange(len(dumps)):
#  out_full.append(avg_dump(n))


# PARALLEL
import multiprocessing
import signal
import psutil
original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(NTHREADS)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  out_list = pool.map(avg_dump, range(len(dumps)))  # This is /real big/ in memory
  print out_list[0].keys()
  for key in out_list[0].keys():
    if key in avg_keys:
      out_full[key] = np.zeros((ND,out_list[0][key].size))
      for n in range(len(out_list)):
        out_full[key][n,:] = out_list[n][key]
    else:
      out_full[key] = np.zeros(ND)
      for n in range(len(out_list)):
        out_full[key][n] = out_list[n][key]
  pool.close()
except KeyboardInterrupt:
  pool.terminate()
else:
  pool.close()
pool.join()

# Toss in the common geom lists
out_full['r'] = r
out_full['th'] = th

# Time average
n = 0
for n in xrange(ND):
  if out_full['t'][n] >= tavg:
    break

print 'nmin = %i' % n

for key in avg_keys:
  out_full[key] = (out_full[key][n:,:]).mean(axis=0)

# Names for compatibility with hdf5_to_dict
out_full['mdot'] = out_full['Mdot']
out_full['phi'] = out_full['Phi']/np.sqrt(np.abs(out_full['Mdot']))

# Variables that can be read once
out_full['a'] = hdr['a']
out_full['t_d'] = diag['t']
out_full['Mdot_d'] = diag['mdot']
out_full['Phi_d'] = diag['Phi']
out_full['Ldot_d'] = diag['ldot']
out_full['Edot_d'] = -diag['edot']
out_full['Lum_d'] = diag['lum_eht']

# OUTPUT
import pickle

pickle.dump(out_full, open('eht_out.p', 'w'))
