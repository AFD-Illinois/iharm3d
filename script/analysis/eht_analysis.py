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

# Sometimes we don't know times but want averages
if tavg < 0:
  tavg = ND//2

# Decide where to measure fluxes
def i_of(rcoord):
  i = 0
  while geom['r'][i,hdr['n2']//2,0] < rcoord:
    i += 1
  i -= 1
  return i

# Leave several extra zones if using MKS3 coordinates
if geom['metric'] == "MKS3":
  iEH = i_of(hdr['r_eh'])+2
else:
  iEH = i_of(hdr['r_eh'])

if floor_workaround_flux:
  iF = i_of(5) # Measure fluxes at r=5M
else:
  iF = iEH

# Max radius when computing "total" energy
iEmax = i_of(40)

# BZ luminosity
# 100M seems like the standard measuring spot (or at least, BHAC does it that way)
# L_BZ seems constant* after that, but much higher within ~50M
iBZ = i_of(40)

jmin, jmax = get_j_vals(geom)

print("Using EH at zone {}, Fluxes at zone {}, Emax within zone {}, L_BZ at zone {}".format(iEH, iF, iEmax, iBZ))

def avg_dump(n):
  out = {}

  # We obviously need the derived variables, but not the extras
  dump = io.load_dump(dumps[n], hdr, geom, extras=False)

  out['t'] = dump['t']
  # When we don't know times, fudge
  if out['t'] == 0 and n != 0:
    out['t'] = n

  print("Loaded {} / {}: {}".format((n+1), len(dumps), out['t']))

  # SHELL AVERAGES (only for t >= tavg usu. tmax/2)
  if out['t'] >= tavg:

    out['rho_r'] = eht_profile(geom, dump['RHO'], jmin, jmax)
    out['Theta_r'] = eht_profile(geom, (hdr['gam']-1.)*dump['UU']/dump['RHO'], jmin, jmax)
    out['B_r'] = eht_profile(geom, np.sqrt(dump['bsq']), jmin, jmax)

    Pg = (hdr['gam']-1.)*dump['UU']
    Pb = dump['bsq']/2
    out['Pg_r'] = eht_profile(geom, Pg, jmin, jmax)
    out['Ptot_r'] = eht_profile(geom, Pg + Pb, jmin, jmax)
    out['betainv_r'] = eht_profile(geom, Pb/Pg, jmin, jmax)

    out['uphi_r'] = eht_profile(geom, dump['ucon'][:,:,:,3], jmin, jmax)

    # THETA AVERAGES
    Fcov01, Fcov13 = Fcov(geom, dump, 0, 1), Fcov(geom, dump, 1, 3)
    out['omega_th'] = theta_av(geom, Fcov01, iEH, 1) / theta_av(geom, Fcov13, iEH, 1)
    out['omega_av_th'] = theta_av(geom, Fcov01, iEH-2, 5) / theta_av(geom, Fcov13, iEH-2, 5)

    # This produces much worse results
    #out['omega_alt_th'] = theta_av(Fcov(dump, 0, 2), iEH, 1) / theta_av(Fcov(dump, 2, 3), iEH, 1)
    #out['omega_alt_av_th'] = theta_av(Fcov(dump, 0, 2), iEH-2, 5) / theta_av(Fcov(dump, 2, 3), iEH-2, 5)

    # For demonstrating jet stuff
    out['Ltot_th'] = theta_av(geom, -T_mixed(dump, 1, 0) - dump['RHO']*dump['ucon'][:,:,:,1], iBZ, 11)

  # The HARM B_unit is sqrt(4pi)*c*sqrt(rho) which has caused issues:
  #norm = np.sqrt(4*np.pi) # This is what I believe matches T,N,M '11 and Narayan '12
  norm = 1 # This is what the EHT comparison uses?
  
  if geom['mixed_metrics']:
    # B1 will be in the _vector_ coordinates.  Must perform the integral in those instead of zone coords
    # Some gymnastics were done to keep in-memory size small
    dxEH = np.einsum("i,...ij->...j", np.array([0, geom['dx1'], geom['dx2'], geom['dx3']]), np.linalg.inv(geom['vec_to_grid'][iEH,:,:,:]))
    out['Phi'] = 0.5*norm * np.sum( np.fabs(dump['B1'][iEH,:,:]) * geom['gdet_vec'][iEH,:,None]*dxEH[:,None,2]*dxEH[:,None,3], axis=(0,1) )
  else:
    out['Phi'] = 0.5*norm*sum_shell(geom, np.fabs(dump['B1']), at_zone=iEH)

  # FLUXES
  # Radial profiles of Mdot and Edot, and their particular values
  # EHT normalization has both these values positive
  if out['t'] >= tavg:
    out['FE_r'] = sum_shell(geom, T_mixed(dump, 1,0))
    out['Edot'] = out['FE_r'][iF]
  else:
    out['Edot'] = sum_shell(geom, T_mixed(dump, 1,0), at_zone=iF)

  # Variable useful in several contexts below
  rho_ur = dump['RHO']*dump['ucon'][:,:,:,1]
  sigma = dump['bsq']/dump['RHO']
  
  if floor_workaround_funnel:
    # TODO implement all of this with 'mask='?
    mdot_full = rho_ur[iF,:,:]*geom['gdet'][iF,:,None]*hdr['dx2']*hdr['dx3']
    sigma_shaped = dump['bsq'][iF,:,:]/dump['RHO'][iF,:,:]
    out['Mdot'] = (mdot_full[np.where(sigma_shaped < 10)]).sum()
    if out['t'] >= tavg:
      out['FM_r'] = -sum_shell(geom, rho_ur, mask=(sigma_shaped < 10))
  else:
    if out['t'] >= tavg:
      out['FM_r'] = -sum_shell(geom, rho_ur)
      out['Mdot'] = out['FM_r'][iF]
    else:
      out['Mdot'] = -sum_shell(geom, rho_ur, at_zone=iF)

  out['Ldot'] = sum_shell(geom, T_mixed(dump, 1,3), at_zone=iF)

  out['sigma_max'] = np.max(sigma)

  # Blandford-Znajek Luminosity L_BZ
  temm = TEM_mixed(dump, 1, 0)
  tfull = T_mixed(dump, 1, 0)
  bernoulli = -T_mixed(dump,0,0) /(dump['RHO']*dump['ucon'][:,:,:,0]) - 1
  if debug:
    # A bunch of radial profiles to test consistency
    out['LBZ_r'] = sum_shell(geom, temm, mask=(sigma > 1))
    out['LBZ'] = out['LBZ_r'][iBZ]
    
    mu = (-tfull + rho_ur) / rho_ur
    out['LBZ_mu2_r'] = sum_shell(geom, temm, mask=np.logical_or(sigma > 1, mu > 2))
    out['LBZ_mu3_r'] = sum_shell(geom, temm, mask=np.logical_or(sigma > 1, mu > 3))
    out['LBZ_mu4_r'] = sum_shell(geom, temm, mask=np.logical_or(sigma > 1, mu > 4))
  else:
    if out['t'] >= tavg:
      out['LBZ_sigma_r'] = sum_shell(geom, temm, mask=(sigma > 1))
      out['LBZ_sigma'] = out['LBZ_r'][iBZ]
      out['Ltot_sigma_r'] = sum_shell(geom, tfull+rho_ur, mask=(sigma > 1))
      out['Ltot_bernoulli'] = out['Ltot_r'][iBZ]

      #ucon_mean = np.mean(dump['ucon'][:,:,:,1], axis=-1)
      #out['Ltot_r'] = sum_shell(geom, tfull+rho_ur, mask=( (ucon_mean > 0)[:,:,None] ))
      out['Ltot_bernoulli_r'] = sum_shell(geom, tfull+rho_ur, mask=(bernoulli > 0.02))
      out['Ltot_bernoulli'] = out['Ltot_r'][iBZ]
      out['LBZ_bernoulli_r'] = sum_shell(geom, temm, mask=(sigma > 1))
      out['LBZ_bernoulli'] = out['LBZ_r'][iBZ]
    else:
      out['LBZ_sigma'] = sum_shell(geom, temm, at_zone=iBZ, mask=(sigma > 1))
      out['Ltot_sigma'] = sum_shell(geom, tfull+rho_ur, at_zone=iBZ, mask=(sigma > 1))
      #ucon_mean = np.mean(dump['ucon'][:,:,:,1], axis=-1)
      #out['Ltot'] = sum_shell(geom, tfull+rho_ur, at_zone=iBZ, mask=( (ucon_mean > 1)[:,:,None] ))
      out['Ltot_bernoulli'] = sum_shell(geom, tfull+rho_ur, at_zone=iBZ, mask=(bernoulli > 0.02))
      out['LBZ_bernoulli'] = sum_shell(geom, temm, at_zone=iBZ, mask=(sigma > 1))

  rho = dump['RHO']
  P = (hdr['gam']-1.)*dump['UU']
  B = np.sqrt(dump['bsq'])
  C = 0.2
  j = rho**3 * P**(-2) * np.exp(-C*(rho**2 / (B*P**2))**(1./3.))
  out['Lum'] = eht_vol(geom, j, jmin, jmax, outside=iEH)

  T00 = T_mixed(dump, 0,0)
  out['Etot'] = sum_vol(geom, T00, within=iEmax)
  #print "Energy on grid: ",out['Etot']

  # For an averaged energy profile
  #out['E_r'] = radial_sum(geom, T00)

  return out

if debug:
  # SERIAL (very slow)
  out_list = [avg_dump(n) for n in range(len(dumps))]
else:
  # PARALLEL
  NTHREADS = util.calc_nthreads(hdr, pad=0.3)
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

# Compute minimum dump for the radial averages
nmin = 0
for n in range(ND):
  if out_list[n]['t'] >= tavg:
    nmin = n
    break

print("nmin = ",nmin)

# Make a dict for merged variables
out_full = {}
# Toss in the common geom lists
N2 = hdr['n2']
out_full['r'] = geom['r'][:,N2//2,0]
out_full['th'] = geom['th'][0,:N2//2,0]

# Merge the output dicts
for key in list(out_list[-1].keys()):
  if key[-2:] == '_r':
    out_full[key] = np.zeros((ND, out_full['r'].size)) # 1D only trick
    for n in range(nmin,ND):
      out_full[key][n,:] = out_list[n][key]
  elif key[-3:] == '_th':
    out_full[key] = np.zeros((ND, out_full['th'].size))
    for n in range(nmin,ND):
      out_full[key][n,:] = out_list[n][key]
  else:
    out_full[key] = np.zeros(ND)
    for n in range(ND):
      out_full[key][n] = out_list[n][key]

# Todo specify radial or profile in key name?
for key in out_full:
  if key[-2:] == '_r' or key[-3:] == '_th':
    out_full[key] = (out_full[key][nmin:,:]).mean(axis=0)

# Compat/completeness stuff
out_full['mdot'] = out_full['Mdot']
out_full['phi'] = out_full['Phi']/np.sqrt(out_full['Mdot'])
out_full['a'] = hdr['a']

# Pass along HARM's own diagnostics for comparison
# TODO implement this separately for HARM pipeline, for new overplotting semantics
#diag = io.load_log(path)

#out_full['t_d'] = diag['t']
#out_full['Mdot_d'] = diag['mdot']
#out_full['Phi_d'] = diag['Phi']
#out_full['Ldot_d'] = diag['ldot']
#out_full['Edot_d'] = diag['edot']
#out_full['Lum_d'] = diag['lum_eht']
#out_full['divbmax_d'] = diag['divbmax']

# OUTPUT
pickle.dump(out_full, open("eht_out.p", "wb"))
