################################################################################
#                                                                              # 
#  INTEGRATE QUANTITIES OF INTEREST                                            #
#                                                                              # 
################################################################################

# TODO hide behind a namespace?
from analysis_fns import *
import hdf5_to_dict as io

import numpy as np
import util
import glob
import os
import sys
import psutil

# Option to calculate fluxes at (just inside) r = 5
# This reduces interference from floors
floor_workaround_flux = False
# Option to ignore accretion at high magnetization (funnel)
# This also reduces interference from floors
floor_workaround_funnel = False

if len(sys.argv) < 2:
  util.warn('Format: python eht_analysis.py DUMP_PATH [tavg]')
  sys.exit()

path = sys.argv[1]

dumps = io.get_dumps_reduced(path)
ND = len(dumps)

hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr, os.path.join(path,'grid.h5'))

if len(sys.argv) > 2:
  tavg = float(sys.argv[2])
else:
  tavg = io.load_dump(dumps[ND/2], hdr)['t'] - 1.0

init_analysis(hdr, geom)


# Limit threads for 192^3+ problem due to memory
# TODO guess NTHREADS better (96e9/N^3*8?)
if hdr['n1'] * hdr['n2'] * hdr['n3'] >= 192 * 192 * 192:
  # Roughly compute memory and leave some generous padding for multiple copies and Python games
  NTHREADS = int(0.11 * psutil.virtual_memory().total/(hdr['n1']*hdr['n2']*hdr['n3']*10*8))
  # Leave the rest of the parallelism to MKL
  try:
    import ctypes

    mkl_rt = ctypes.CDLL('libmkl_rt.so')
    mkl_set_num_threads = mkl_rt.MKL_Set_Num_Threads
    mkl_get_max_threads = mkl_rt.MKL_Get_Max_Threads
    mkl_set_num_threads(4)
    print "Using", mkl_get_max_threads(), "MKL threads"
  except Error as e:
    print(e)
else:
  NTHREADS = psutil.cpu_count(logical=False)

diag = io.load_log(hdr, os.path.join(path, "log.out"))

dx1 = hdr['dx1']
dx2 = hdr['dx2']
dx3 = hdr['dx3']
N1 = hdr['n1']
N2 = hdr['n2']
N3 = hdr['n3']
gam = hdr['gam']

r = geom['r'][:,N2/2,0]
th = geom['th'][0,:N2/2,0]

def i_of(rcoord):
  i = 0
  while r[i] < rcoord:
    i += 1
  i -= 1
  return i

iF = 5 # Zone 5 = rEH
if floor_workaround_flux:
  iF = i_of(5) # Measure fluxes at r=5M

iEmax = i_of(40)

# BZ luminosity
debug_lbz = False
iBZ = i_of(10) # TODO is there a standard measuring spot?

# Some variables (Phi) should only be computed at EH=zone 5
iEH = 5

avg_keys = ['rho_r', 'Theta_r', 'B_r', 'Pg_r', 'Ptot_r', 'betainv_r', 'uphi_r', 'FE_r', 'FM_r', 'rho_SADW', 'omega_th', 'omega_th_av' ] #, 'omega_th_alt', 'omega_th_alt_av']

def avg_dump(n):
  out = {}

  dump = io.load_dump(dumps[n], hdr, geom, extras=False)
  out['t'] = dump['t']
  print "Loaded ",(n+1),"/",len(dumps),": ",out['t']

  # SHELL AVERAGES (only for t > tavg usu. tmax/2)
  if out['t'] > tavg:

    out['rho_SADW'] = WAVG(dump['RHO'], dump['RHO'])

    out['rho_r'] = eht_profile(dump['RHO'])
    
    Theta = (hdr['gam']-1.)*dump['UU']/dump['RHO']
    out['Theta_r'] = eht_profile(Theta)
    
    B = np.sqrt(dump['bsq'])
    out['B_r'] = eht_profile(B)
  
    Pg = (hdr['gam']-1.)*dump['UU']
    out['Pg_r'] = eht_profile(Pg)
  
    Ptot = Pg + dump['bsq']/2
    out['Ptot_r'] = eht_profile(Ptot)
  
    betainv = (dump['bsq']/2)/Pg
    out['betainv_r'] = eht_profile(betainv)
  
    uphi = (dump['ucon'][:,:,:,3])
    out['uphi_r'] = eht_profile(uphi)
  
    # THETA AVERAGES
    out['omega_th'] = theta_av(Fcov(dump, 0, 1), iEH, 1) / theta_av(Fcov(dump, 1, 3), iEH, 1)
    out['omega_th_av'] = theta_av(Fcov(dump, 0, 1), iEH-2, 5) / theta_av(Fcov(dump, 1, 3), iEH-2, 5)

    # This produces much worse results
    #out['omega_th_alt'] = theta_av(Fcov(dump, 0, 2), iEH, 1) / theta_av(Fcov(dump, 2, 3), iEH, 1)
    #out['omega_th_alt_av'] = theta_av(Fcov(dump, 0, 2), iEH-2, 5) / theta_av(Fcov(dump, 2, 3), iEH-2, 5)

  else:
    for key in avg_keys:
      if '_r' in key or '_SADW' in key:
        out[key] = np.zeros_like(r)
      elif '_th' in key:
        out[key] = np.zeros_like(th)

  # The HARM B_unit is sqrt(4pi)*c*sqrt(rho) which has caused issues:
  norm = np.sqrt(4*np.pi) # This is what I believe matches T,N,M '11 and Narayan '12
  #norm = 1 # This is what the EHT comparison used for a while
  out['Phi'] = 0.5*norm*sum_shell_at(np.fabs(dump['B1']), iEH)

  # FLUXES
  # Radial profiles of Mdot and Edot, and their particular values
  out['FE_r'] = sum_shell(Tmixed(dump, 1,0))
  out['Edot'] = out['FE_r'][iF]

  out['FM_r'] = -sum_shell(dump['RHO']*dump['ucon'][:,:,:,1])
  if floor_workaround_funnel:
    mdot_full = dump['RHO'][iF,:,:]*dump['ucon'][iF,:,:,1]*geom['gdet'][iF,:,None]*dx2*dx3
    sigma_shaped = dump['bsq'][iF,:,:]/dump['RHO'][iF,:,:]
    out['Mdot'] = (mdot_full[np.where(sigma_shaped < 10)]).sum()
  else:
    out['Mdot'] = out['FM_r'][iF]

  out['Ldot'] = sum_shell_at(Tmixed(dump, 1,3), iF)

  # Maximum magnetization (and allow re-use of the variable)
  sigma = dump['bsq']/dump['RHO']
  out['sigma_max'] = np.max(sigma)

  # Blandford-Znajek Luminosity L_BZ
  # TODO define T_EM, and use sum_shell_at function
  LBZ = lambda i: (dx2*dx3*geom['gdet'][i,:,None]*(dump['bsq'][i,:,:]*dump['ucon'][i,:,:,1]*dump['ucov'][i,:,:,0] - dump['bcon'][i,:,:,1]*dump['bcov'][i,:,:,0])[np.where(sigma[i,:,:]>1)] ).sum()

  out['LBZ'] = LBZ(iBZ)

  if debug_lbz:
    out['LBZ_10'] = LBZ(10)
    out['LBZ_30'] = LBZ(30)
    out['LBZ_50'] = LBZ(50)
    if N1 > 80:
      out['LBZ_80'] = LBZ(80)

    #print "L_BZ at ",out['t']," is ",[LBZ(i) for i in range(10,100,10)]

  rho = dump['RHO']
  P = (hdr['gam']-1.)*dump['UU']
  B = np.sqrt(dump['bsq'])
  C = 0.2
  j = rho**3*P**(-2)*np.exp(-C*(rho**2/(B*P**2))**(1./3.))
  out['Lum'] = (j*geom['gdet'][:,:,None]*dx1*dx2*dx3).sum()

  # TODO define a sum_grid for this and above?
  out['Etot'] = np.sum(Tmixed(dump, 0,0)[:iEmax,:,:]*geom['gdet'][:iEmax,:,None]*dx1*dx2*dx3)
  #print "Energy on grid: ",out['Etot']

  # For an averaged energy profile
  #out['E_r'] = radial_sum(Tmixed(dump, 0,0))

  return out

out_full = {}

# SERIAL (very slow)
#out_list = [avg_dump(n) for n in range(len(dumps))]

# PARALLEL
# TODO runtime parallel/serial option
import multiprocessing
import psutil

pool = multiprocessing.Pool(NTHREADS)
try:
  # Map the above function to the dump numbers, returning a list of 'out' dicts
  out_list = pool.map_async(avg_dump, range(len(dumps))).get(99999999)
  #print out_list[0].keys()
except KeyboardInterrupt:
  pool.terminate()
  pool.join()
else:
  pool.close()
  pool.join()

# Merge the output dicts
for key in out_list[0].keys():
  if key in avg_keys:
    out_full[key] = np.zeros((ND,out_list[0][key].size))
    for n in range(len(out_list)):
      out_full[key][n,:] = out_list[n][key]
  else:
    out_full[key] = np.zeros(ND)
    for n in range(len(out_list)):
      out_full[key][n] = out_list[n][key]

# Toss in the common geom lists
out_full['r'] = r
out_full['th'] = th

# Time average the radial profiles
n = 0
for n in xrange(ND):
  if out_full['t'][n] >= tavg:
    break

print "nmin = ",n

# Todo specify radial or profile in key name?
for key in avg_keys:
  out_full[key] = (out_full[key][n:,:]).mean(axis=0)

# Names for compatibility with hdf5_to_dict
out_full['mdot'] = out_full['Mdot']
out_full['phi'] = out_full['Phi']/np.sqrt(np.abs(out_full['Mdot']))

# Pass along HARM's own diagnostics for comparison
out_full['a'] = hdr['a']
out_full['t_d'] = diag['t']
out_full['Mdot_d'] = diag['mdot']
out_full['Phi_d'] = diag['Phi']
out_full['Ldot_d'] = diag['ldot']
out_full['Edot_d'] = -diag['edot']
out_full['Lum_d'] = diag['lum_eht']
out_full['divbmax_d'] = diag['divbmax']

# OUTPUT
import pickle

pickle.dump(out_full, open('eht_out.p', 'w'))
