################################################################################
#                                                                              # 
#  CALCULATE TIME-DEPENDENT AND TIME-AVERAGED QUANTITIES                       #
#                                                                              # 
################################################################################

from __future__ import print_function, division

from analysis_fns import *
import hdf5_to_dict as io
import util

import sys
import multiprocessing
import psutil
import pickle

import numpy as np

# Option to calculate fluxes at (just inside) r = 5
# This reduces interference from floors
floor_workaround_flux = False
# Option to ignore accretion at high magnetization (funnel)
# This also reduces interference from floors
floor_workaround_funnel = False

if len(sys.argv) < 2:
  util.warn('Format: python eht_analysis.py /path/to/dumps [tavg]')
  sys.exit()

# This doesn't seem like the _right_ way to do optional args
tavg = None
if sys.argv[1] == "-d":
  debug = True
  path = sys.argv[2]
  if len(sys.argv) > 3:
    tavg = float(sys.argv[3])
else:
  debug = False
  path = sys.argv[1]
  if len(sys.argv) > 2:
    tavg = float(sys.argv[2])

dumps = io.get_dumps_list(path)
ND = len(dumps)

hdr = io.load_hdr(dumps[0])
geom = io.load_geom(hdr, path)

# If the time after which to average wasn't given, just use the back half of dumps
if tavg is None:
  tavg = io.load_dump(dumps[ND//2], hdr, geom, derived_vars=False, extras=False)['t'] - 1.0

# Decide where to measure fluxes
def i_of(rcoord):
  i = 0
  while geom['r'][i,0,0] < rcoord:
    i += 1
  i -= 1
  return i

iF = 5 # Zone 5 = rEH
if floor_workaround_flux:
  iF = i_of(5) # Measure fluxes at r=5M

# Max radius when computing "total" energy
iEmax = i_of(40)

# BZ luminosity
debug_lbz = False
iBZ = i_of(10) # TODO is there a standard measuring spot?

# Some variables (Phi) should only be computed at EH=zone 5
iEH = 5

THMIN = np.pi/3.
THMAX = 2.*np.pi/3.
# Calculate jmin, jmax for EHT radial profiles
ths = geom['th'][-1,:,0]
for n in range(len(ths)):
  if ths[n] > THMIN:
    jmin = n
    break
for n in range(len(ths)):
  if ths[n] > THMAX:
    jmax = n
    break

# Variables which should be averaged
avg_keys = ['rho_r', 'Theta_r', 'B_r', 'Pg_r', 'Ptot_r', 'betainv_r', 'uphi_r', 'FE_r', 'FM_r', 'omega_th', 'omega_th_av' ] #, 'omega_th_alt', 'omega_th_alt_av']

def avg_dump(n):
  out = {}

  # We obviously need the derived variables, but not the extras
  dump = io.load_dump(dumps[n], hdr, geom, extras=False)

  out['t'] = dump['t']
  print("Loaded {} / {}: {}".format((n+1), len(dumps), out['t']))

  # SHELL AVERAGES (only for t > tavg usu. tmax/2)
  if out['t'] > tavg:

    out['rho_r'] = eht_profile(geom, dump['RHO'], jmin, jmax)

    Theta = (hdr['gam']-1.)*dump['UU']/dump['RHO']
    out['Theta_r'] = eht_profile(geom, Theta, jmin, jmax)

    B = np.sqrt(dump['bsq'])
    out['B_r'] = eht_profile(geom, B, jmin, jmax)

    Pg = (hdr['gam']-1.)*dump['UU']
    out['Pg_r'] = eht_profile(geom, Pg, jmin, jmax)

    Ptot = Pg + dump['bsq']/2
    out['Ptot_r'] = eht_profile(geom, Ptot, jmin, jmax)

    betainv = (dump['bsq']/2)/Pg
    out['betainv_r'] = eht_profile(geom, betainv, jmin, jmax)

    uphi = (dump['ucon'][:,:,:,3])
    out['uphi_r'] = eht_profile(geom, uphi, jmin, jmax)

    # THETA AVERAGES
    Fcov01, Fcov13 = Fcov(geom, dump, 0, 1), Fcov(geom, dump, 1, 3)
    out['omega_th'] = theta_av(geom, Fcov01, iEH, 1) / theta_av(geom, Fcov13, iEH, 1)
    out['omega_th_av'] = theta_av(geom, Fcov01, iEH-2, 5) / theta_av(geom, Fcov13, iEH-2, 5)

    # This produces much worse results
    #out['omega_th_alt'] = theta_av(Fcov(dump, 0, 2), iEH, 1) / theta_av(Fcov(dump, 2, 3), iEH, 1)
    #out['omega_th_alt_av'] = theta_av(Fcov(dump, 0, 2), iEH-2, 5) / theta_av(Fcov(dump, 2, 3), iEH-2, 5)

  else:
    for key in avg_keys:
      if '_r' in key:
        out[key] = np.zeros((hdr['n1']))
      elif '_th' in key:
        out[key] = np.zeros((hdr['n2']//2))

  # The HARM B_unit is sqrt(4pi)*c*sqrt(rho) which has caused issues:
  #norm = np.sqrt(4*np.pi) # This is what I believe matches T,N,M '11 and Narayan '12
  norm = 1 # This is what the EHT comparison uses?
  out['Phi'] = 0.5*norm*sum_shell(geom, np.fabs(dump['B1']), at_zone=iEH)

  # FLUXES
  # Radial profiles of Mdot and Edot, and their particular values
  # EHT normalization has both these values positive
  out['FE_r'] = -sum_shell(geom, Tmixed(geom, dump, 1,0))
  out['Edot'] = out['FE_r'][iF]

  out['FM_r'] = -sum_shell(geom, dump['RHO']*dump['ucon'][:,:,:,1])
  if floor_workaround_funnel:
    mdot_full = dump['RHO'][iF,:,:]*dump['ucon'][iF,:,:,1]*geom['gdet'][iF,:,None]*hdr['dx2']*hdr['dx3']
    sigma_shaped = dump['bsq'][iF,:,:]/dump['RHO'][iF,:,:]
    out['Mdot'] = (mdot_full[np.where(sigma_shaped < 10)]).sum()
  else:
    out['Mdot'] = out['FM_r'][iF]

  out['Ldot'] = sum_shell(geom, Tmixed(geom, dump, 1,3), at_zone=iF)

  # Maximum magnetization (and allow re-use of the variable)
  sigma = dump['bsq']/dump['RHO']
  out['sigma_max'] = np.max(sigma)

  # Blandford-Znajek Luminosity L_BZ
  # TODO define T_EM, and use sum_shell_at function
  LBZ = lambda i: (hdr['dx2']*hdr['dx3']*geom['gdet'][i,:,None]*(dump['bsq'][i,:,:]*dump['ucon'][i,:,:,1]*dump['ucov'][i,:,:,0] - dump['bcon'][i,:,:,1]*dump['bcov'][i,:,:,0])[np.where(sigma[i,:,:]>1)] ).sum()

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
  j = rho**3 * P**(-2) * np.exp(-C*(rho**2 / (B*P**2))**(1./3.))
  out['Lum'] = eht_vol(geom, j, jmin, jmax, outside=iEH)

  T00 = Tmixed(geom, dump, 0,0)
  out['Etot'] = sum_vol(geom, T00, within=iEmax)
  #print "Energy on grid: ",out['Etot']

  # For an averaged energy profile
  #out['E_r'] = radial_sum(geom, T00)

  return out

out_full = {}

if debug:
  # SERIAL (very slow)
  out_list = [avg_dump(n) for n in range(len(dumps))]
else:
  # PARALLEL
  NTHREADS = util.calc_nthreads(hdr)
  pool = multiprocessing.Pool(NTHREADS)
  try:
    # Map the above function to the dump numbers, returning a list of 'out' dicts
    out_list = pool.map_async(avg_dump, list(range(len(dumps)))).get(99999999)
    #print out_list[0].keys()
  except KeyboardInterrupt:
    pool.terminate()
    pool.join()
  else:
    pool.close()
    pool.join()

# Merge the output dicts
for key in list(out_list[0].keys()):
  if key in avg_keys:
    out_full[key] = np.zeros((ND,out_list[0][key].size))
    for n in range(len(out_list)):
      out_full[key][n,:] = out_list[n][key]
  else:
    out_full[key] = np.zeros(ND)
    for n in range(len(out_list)):
      out_full[key][n] = out_list[n][key]

# Toss in the common geom lists
N2 = hdr['n2']
out_full['r'] = geom['r'][:,N2//2,0]
out_full['th'] = geom['th'][0,:N2//2,0]

# Time average the radial profiles
n = 0
for n in range(ND):
  if out_full['t'][n] >= tavg:
    break

print("nmin = ",n)

# Todo specify radial or profile in key name?
for key in avg_keys:
  out_full[key] = (out_full[key][n:,:]).mean(axis=0)

# Names for compatibility with hdf5_to_dict
out_full['mdot'] = out_full['Mdot']
out_full['phi'] = out_full['Phi']/np.sqrt(np.abs(out_full['Mdot']))

# Pass along HARM's own diagnostics for comparison
diag = io.load_log(path)

out_full['a'] = hdr['a']
out_full['t_d'] = diag['t']
out_full['Mdot_d'] = diag['mdot']
out_full['Phi_d'] = diag['Phi']
out_full['Ldot_d'] = diag['ldot']
out_full['Edot_d'] = diag['edot']
out_full['Lum_d'] = diag['lum_eht']
out_full['divbmax_d'] = diag['divbmax']

# OUTPUT
pickle.dump(out_full, open("eht_out.p", "wb"))
