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
NTHREADS = psutil.cpu_count(logical=False)
#NTHREADS = 1
THMIN = np.pi/3.
THMAX = 2.*np.pi/3.

if len(sys.argv) != 3:
  util.warn('Format: python eht_analysis.py [dump path] [tavg]')
  sys.exit()

path = sys.argv[1]
tavg = float(sys.argv[2])

dumps = io.get_dumps_reduced(path)

hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr, os.path.join(path,'grid.h5'))

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
N1 = hdr['N1']
N2 = hdr['N2']
N3 = hdr['N3']

r = geom['r'][:,N2/2,0]
th = geom['th'][0,:N2/2,0]

# Calculate fluxes at (just inside) r = 5
# Except Phi at EH, zone 5
iF = 0
while r[iF] < 5:
  iF += 1
iF -= 1

iEH = 5  

ND = len(dumps)

t = np.zeros(ND)
rho_r = np.zeros([ND, N1])
Theta_r = np.zeros([ND, N1])
B_r = np.zeros([ND, N1])
Pg_r = np.zeros([ND, N1])
Ptot_r = np.zeros([ND, N1])
betainv_r = np.zeros([ND, N1])
uphi_r = np.zeros([ND, N1])
rho_SADW = np.zeros([ND, N1])
omega_th = np.zeros([ND, N2/2])

Mdot = np.zeros(ND)
Phi = np.zeros(ND)
Ldot = np.zeros(ND)
Edot = np.zeros(ND)
Lum = np.zeros(ND)

# Pointers are magical things
out = {}
out['rho_SADW'] = rho_SADW
out['r'] = r
out['th'] = th
out['omega_th'] = omega_th
out['rho_r'] = rho_r
out['Theta_r'] = Theta_r
out['B_r'] = B_r
out['Pg_r'] = Pg_r
out['Ptot_r'] = Ptot_r
out['betainv_r'] = betainv_r
out['uphi_r'] = uphi_r
out['t'] = t
out['Mdot'] = Mdot
out['Phi'] = Phi
out['Ldot'] = Ldot
out['Edot'] = Edot
out['Lum'] = Lum

# EVALUATE DIAGNOSTICS
  
vol = (dx2*2.*np.pi*geom['gdet'][:,:]).sum(axis=-1)

def INT(var):
  return (dx2*dx3*geom['gdet'][:,:,None]*var[:,:,:]).sum(axis=-1).sum(axis=-1)

def WAVG(var, w):
  return INT(w*var)/INT(w)

def avg_dump(n):
  print '%08d / ' % (n+1) + '%08d' % len(dumps) 
  dump = io.load_dump(dumps[n], geom, hdr)
  t[n] = dump['t']
  print dump['t']

  rho_SADW[n,:] = WAVG(dump['RHO'], dump['RHO'])

  # SHELL AVERAGES
  #vol = (dx2*dx3*geom['gdet'][:,:]).sum(axis=-1)
  rho_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*dump['RHO'][:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)
  
  Theta = (hdr['gam']-1.)*dump['UU'][:,:,:]/dump['RHO'][:,:,:]
  Theta_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Theta[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)
  
  B = np.sqrt(dump['bsq'][:,:,:])
  B_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*B[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  Pg = (hdr['gam']-1.)*dump['UU'][:,:,:]
  Pg_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Pg[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  Ptot = Pg + dump['bsq'][:,:,:]/2
  Ptot_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*Ptot[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  betainv = (dump['bsq'][:,:,:]/2)/Pg
  betainv_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*betainv[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1)

  uphi = (dump['ucon'][:,:,:,3])
  uphi_r[n,:] = (dx2*dx3*geom['gdet'][:,jmin:jmax,None]*uphi[:,jmin:jmax,:]).sum(axis=-1).sum(axis=-1) 

  rho_r[n,:] /= vol
  Theta_r[n,:] /= vol
  B_r[n,:] /= vol
  Pg_r[n,:] /= vol
  Ptot_r[n,:] /= vol
  betainv_r[n,:] /= vol
  uphi_r[n,:] /= vol

  # FLUXES
  #Phi_bcon[n] = 0.5*(np.fabs(dump['bcon'][iF,:,:,1])*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  # The HARM B_unit is sqrt(4pi)*c*sqrt(rho) so we need to multiply that in here
  Phi[n] = 0.5*(np.fabs(np.sqrt(4*np.pi)*dump['B1'][iEH,:,:])*geom['gdet'][iEH,:,None]*dx2*dx3).sum()
  omega_th[n,:] = (dump['omega'][iEH,:,:].sum(axis=-1) + dump['omega'][iEH,::-1,:].sum(axis=-1))[:N2/2]
  
  Mdot[n] = np.abs((dump['RHO'][iF,:,:]*dump['ucon'][iF,:,:,1]*geom['gdet'][iF,:,None]*dx2*dx3).sum())
  Trphi = (dump['RHO'] + dump['UU'] + (hdr['gam']-1.)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,1]*dump['ucov'][:,:,:,3]
  Trphi -= dump['bcon'][:,:,:,1]*dump['bcov'][:,:,:,3]
  Trt = (dump['RHO'] + dump['UU'] + (hdr['gam']-1.)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,1]*dump['ucov'][:,:,:,0]
  Trt -= dump['bcon'][:,:,:,1]*dump['bcov'][:,:,:,0]
  Ldot[n] = (Trphi[iF,:,:]*geom['gdet'][iF,:,None]*dx2*dx3).sum()
  Edot[n] = -(Trt[iF,:,:]*geom['gdet'][iF,:,None]*dx2*dx3).sum()

  rho = dump['RHO']
  P = (hdr['gam']-1.)*dump['UU']
  B = np.sqrt(dump['bsq'])
  C = 0.2
  j = rho**3*P**(-2)*np.exp(-C*(rho**2/(B*P**2))**(1./3.))
  Lum[n] = (j*geom['gdet'][:,:,None]*dx1*dx2*dx3).sum()
  
  return out

# SERIAL
#for n in xrange(len(dumps)):
#  avg_dump(n)

# PARALLEL
import multiprocessing
import signal
import psutil
original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(NTHREADS)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  out_list = pool.map(avg_dump, range(len(dumps)))  # This is /real big/ in memory
  for key in out.keys():
    for i in range(len(out_list)):
      slc = np.where(out[key] == 0) # Avoid double-fill from threads that got 2+ indices
      out[key][slc] += out_list[i][key][slc]
  pool.close()
except KeyboardInterrupt:
  pool.terminate()
else:
  pool.close()
pool.join()

# Sort
order = np.argsort(out['t'])
for key in out.keys():
  if key in ['rho_r', 'Theta_r', 'B_r', 'Pg_r', 'Ptot_r', 'betainv_r', 'uphi_r', 'omega_th']:
    out[key] = out[key][order,:]
  elif key in ['r', 'th']:
    pass
  else:
    out[key] = out[key][order]

# Time average
n = 0
for n in xrange(ND):
  if t[n] >= tavg:
    break

print 'nmin = %i' % n

for key in ['rho_SADW', 'rho_r', 'Theta_r', 'B_r', 'Pg_r', 'Ptot_r', 'betainv_r', 'uphi_r', 'omega_th']:
  out[key] = (out[key][n:,:]).mean(axis=0)

# Names for compatibility with hdf5_to_dict
out['mdot'] = out['Mdot']
out['phi'] = out['Phi']/np.sqrt(np.abs(out['Mdot']))

# Variables that can be read once
out['a'] = hdr['a']
out['t_d'] = diag['t']
out['Mdot_d'] = diag['mdot']
out['Phi_d'] = diag['Phi']
out['Ldot_d'] = diag['ldot']
out['Edot_d'] = -diag['edot']
out['Lum_d'] = diag['lum_eht']

# OUTPUT
import pickle

pickle.dump(out, open('eht_out.p', 'w'))
