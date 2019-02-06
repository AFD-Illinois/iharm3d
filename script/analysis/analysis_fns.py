
# Convenient analysis functions for physical calculations and averages
# Meant to be imported "from analysis_fns import *" for convenience

import numpy as np

# Define a dict of names, coupled with the functions required to obtain their variables.
# That way, we only need to specify lists and final operations below,
# AND don't need to cart all these things around in memory
d_fns = {'rho': lambda dump: dump['RHO'],
               'bsq' : lambda dump: dump['bsq'],
               'sigma' : lambda dump: dump['bsq'] / dump['RHO'],
               'U' : lambda dump: dump['UU'],
               'Theta' : lambda dump: (dump['hdr']['gam'] - 1.) * dump['UU']/dump['RHO'],
               'u_t' : lambda dump: dump['ucov'][:, :, :, 0],
               'u_phi' : lambda dump: dump['ucov'][:, :, :, 3],
               'u^phi' : lambda dump: dump['ucon'][:, :, :, 3],
               'FM' : lambda dump: dump['RHO'] * dump['ucon'][:, :, :, 1],
               'FE' : lambda dump: -T_mixed(dump, 1, 0),
               'FE_EM' : lambda dump: -TEM_mixed(dump, 1, 0),
               'FE_Fl' : lambda dump: -TFl_mixed(dump, 1, 0),
               'FL' : lambda dump: T_mixed(dump, 1, 3),
               'FL_EM' : lambda dump: TEM_mixed(dump, 1, 3),
               'FL_Fl' : lambda dump: TFl_mixed(dump, 1, 3),
               'Be_b' : lambda dump: bernoulli(dump, with_B=True),
               'Be_nob' : lambda dump: bernoulli(dump, with_B=False),
               'Pg' : lambda dump: (dump['hdr']['gam'] - 1.) * dump['UU'],
               'Pb' : lambda dump: dump['bsq'] / 2,
               'Ptot' : lambda dump: d_fns['Pg'](dump) + d_fns['Pb'](dump),
               'beta' : lambda dump: dump['beta'],
               'B' : lambda dump: np.sqrt(dump['bsq']),
               'betagamma' : lambda dump: np.sqrt((d_fns['FE_EM'](dump) + d_fns['FE_Fl'](dump))/d_fns['FM'](dump) - 1)
               }
               #'rur' : lambda dump: geom['r']*dump['ucon'][:,:,:,1],
               #'gamma' : lambda dump: get_gamma(geom, dump)}

## Physics functions ##

def T_con(geom, dump, i, j):
  gam = dump['hdr']['gam']
  return ( (dump['RHO'] + gam*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucon'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcon'][:,:,None,i,j] - dump['bcon'][:,:,:,i]*dump['bcon'][:,:,:,j] )

def T_cov(geom, dump, i, j):
  gam = dump['hdr']['gam']
  return ( (dump['RHO'] + gam*dump['UU'] + dump['bsq'])*dump['ucov'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcov'][:,:,None,i,j] - dump['bcov'][:,:,:,i]*dump['bcov'][:,:,:,j] )

def T_mixed(dump, i, j):
  gam = dump['hdr']['gam']
  gmixedij = (i == j)
  return ( (dump['RHO'] + gam*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*gmixedij - dump['bcon'][:,:,:,i]*dump['bcov'][:,:,:,j] )

# TODO These only works for i != j.  Need new code path to keep speed?
def TEM_mixed(dump, i, j):
  if i == j: raise ValueError("TEM Not implemented for i==j.")
  return dump['bsq'][:,:,:]*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j] - dump['bcon'][:,:,:,i]*dump['bcov'][:,:,:,j]

def TFl_mixed(dump, i, j):
  gam = dump['hdr']['gam']
  if i == j: raise ValueError("TEM Not implemented for i==j.")
  return (dump['RHO'] + dump['hdr']['gam']*dump['UU'])*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j]

# Return the i,j component of contravarient Maxwell tensor
# TODO there's a computationally easier way to do this:
# Pre-populate an antisym ndarray and einsum
# Same below
def Fcon(geom, dump, i, j):
  NDIM = dump['hdr']['n_dim']

  Fconij = np.zeros_like(dump['RHO'])
  if i != j:
    for mu in range(NDIM):
      for nu in range(NDIM):
        Fconij[:,:,:] += _antisym(i,j,mu,nu) * dump['ucov'][:,:,:,mu] * dump['bcov'][:,:,:,nu]

  # Specify we want gdet in the vectors' coordinate system (this matters for KORAL dump files)
  # TODO is normalization correct?
  return Fconij*geom['gdet'][:,:,None]

def Fcov(geom, dump, i, j):
  NDIM = dump['hdr']['n_dim']

  Fcovij = np.zeros_like(dump['RHO'])
  for mu in range(NDIM):
    for nu in range(NDIM):
      Fcovij += Fcon(geom, dump, mu, nu)*geom['gcov'][:,:,None,mu,i]*geom['gcov'][:,:,None,nu,j]
  
  return Fcovij

def bernoulli(dump, with_B=False):
  if with_B:
    return -T_mixed(dump,0,0) /(dump['RHO']*dump['ucon'][:,:,:,0]) - 1
  else:
    return -(1 + dump['hdr']['gam']*dump['UU']/dump['RHO'])*dump['ucov'][:,:,:,0] - 1

# This is in zone metric!
def lower(geom, vec):
  return np.einsum("...i,...ij->...j", vec, geom['gcov'][:,:,None,:,:])

def to_zone_coords(geom, vec):
  return np.einsum("...i,...ij->...j", vec, geom['vec_to_grid'][:,:,None,:,:])

# Compute 4-vectors given fluid state
# Always returns vectors in the _grid_ coordinate system, to simplify analysis
def get_state(hdr, geom, dump, return_gamma=False):
  ucon = np.zeros([hdr['n1'],hdr['n2'],hdr['n3'],hdr['n_dim']])
  ucov = np.zeros_like(ucon)
  bcon = np.zeros_like(ucon)
  bcov = np.zeros_like(ucon)

  # Aliases to make the below more readable
  if geom['mixed_metrics']:
    # Make sure these are in the vector metric if mixing
    gcov = geom['gcov_vec']
    gcon = geom['gcon_vec']
    alpha = geom['lapse_vec']
  else:
    gcov = geom['gcov']
    gcon = geom['gcon']
    alpha = geom['lapse']

  B1 = dump['B1']
  B2 = dump['B2']
  B3 = dump['B3']

  gamma = get_gamma(geom, dump)

  ucon[:,:,:,0] = gamma/(alpha[:,:,None])
  ucon[:,:,:,1] = dump['U1'] - gamma*alpha[:,:,None]*gcon[:,:,None,0,1]
  ucon[:,:,:,2] = dump['U2'] - gamma*alpha[:,:,None]*gcon[:,:,None,0,2]
  ucon[:,:,:,3] = dump['U3'] - gamma*alpha[:,:,None]*gcon[:,:,None,0,3]

  ucov = np.einsum("...i,...ij->...j", ucon, gcov[:,:,None,:,:])
  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = (B1 + bcon[:,:,:,0]*ucon[:,:,:,1])/ucon[:,:,:,0]
  bcon[:,:,:,2] = (B2 + bcon[:,:,:,0]*ucon[:,:,:,2])/ucon[:,:,:,0]
  bcon[:,:,:,3] = (B3 + bcon[:,:,:,0]*ucon[:,:,:,3])/ucon[:,:,:,0]

  if geom['mixed_metrics']:
    # Convert all 4-vectors to zone coordinates
    ucon = np.einsum("...i,...ij->...j", ucon, geom['vec_to_grid'][:,:,None,:,:])
    ucov = np.einsum("...i,...ij->...j", ucon, geom['gcov'][:,:,None,:,:]) # Lower with _zone_ metric
    bcon = np.einsum("...i,...ij->...j", bcon, geom['vec_to_grid'][:,:,None,:,:])
    bcov = np.einsum("...i,...ij->...j", bcon, geom['gcov'][:,:,None,:,:])
  else:
    # Already have ucov in this case
    bcov = np.einsum("...i,...ij->...j", bcon, gcov[:,:,None,:,:])

  if return_gamma:
    return ucon, ucov, bcon, bcov, gamma
  else:
    return ucon, ucov, bcon, bcov

def get_gamma(geom, dump):
  # Aliases to make the below more readable
  if geom['mixed_metrics']:
    # Make sure this is in the vector metric if mixing
    gcov = geom['gcov_vec']
  else:
    gcov = geom['gcov']

  U1 = dump['U1']
  U2 = dump['U2']
  U3 = dump['U3']

  qsq = (gcov[:,:,None,1,1]*U1**2 + gcov[:,:,None,2,2]*U2**2 +
         gcov[:,:,None,3,3]*U3**2 + 2.*(gcov[:,:,None,1,2]*U1*U2 +
                                        gcov[:,:,None,1,3]*U1*U3 +
                                        gcov[:,:,None,2,3]*U2*U3))
  return np.sqrt(1. + qsq)

# Decide where to measure fluxes
def i_of(geom, rcoord):
  i = 0
  while geom['r'][i,geom['n2']//2,0] < rcoord:
    i += 1
  i -= 1
  return i

## Sums and Averages ##
  
# Var must be a 3D array i.e. a grid scalar
# TODO could maybe be made faster with 'where' but also harder to get right
def sum_shell(geom, var, at_zone=None, mask=None):
  if mask is not None:
    integrand = (var * geom['gdet'][:,:,None]*geom['dx2']*geom['dx3'])*(mask)
  else:
    integrand = var * geom['gdet'][:,:,None]*geom['dx2']*geom['dx3']

  if at_zone is not None:
    return np.sum(integrand[at_zone,:,:], axis=(0,1))
  else:
    return np.sum(integrand, axis=(1,2))

# TODO just pass slices here & below?
def sum_vol(geom, var, within=None):
  if within is not None:
    return np.sum(var[:within,:,:] * geom['gdet'][:within,:,None]*geom['dx1']*geom['dx2']*geom['dx3'])
  else:
    return np.sum(var * geom['gdet'][:,:,None]*geom['dx1']*geom['dx2']*geom['dx3'])

def eht_vol(geom, var, jmin, jmax, outside=None):
  if outside is not None:
    return np.sum(var[outside:,jmin:jmax,:] * geom['gdet'][outside:,jmin:jmax,None]*geom['dx1']*geom['dx2']*geom['dx3'])
  else:
    return np.sum(var[:,jmin:jmax,:] * geom['gdet'][:,jmin:jmax,None]*geom['dx1']*geom['dx2']*geom['dx3'])

# TODO can I cache the volume instead of passing these?
def get_j_vals(geom):
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

  return jmin, jmax

# TODO can I cache the volume instead of passing these?
def eht_profile(geom, var, jmin, jmax):
  return ( (var[:,jmin:jmax,:] * geom['gdet'][:,jmin:jmax,None]*geom['dx2']*geom['dx3']).sum(axis=(1,2)) /
           ((geom['gdet'][:,jmin:jmax]*geom['dx2']).sum(axis=1)*2*np.pi) )

def theta_av(geom, var, start, zones_to_av=1, use_gdet=False):
  # Sum theta from each pole to equator and take overall mean
  N2 = geom['n2']
  if use_gdet:
    return (var[start:start+zones_to_av,:N2//2,:] * geom['gdet'][start:start+zones_to_av,:N2//2,None]*geom['dx1']*geom['dx3'] +
              var[start:start+zones_to_av,:N2//2-1:-1,:] * geom['gdet'][start:start+zones_to_av,:N2//2-1:-1,None]*geom['dx1']*geom['dx3']).sum(axis=(0,2))\
           /((geom['gdet'][start:start+zones_to_av,:N2//2]*geom['dx1']).sum(axis=0)*2*np.pi)
  else:
    return (var[start:start+zones_to_av,:N2//2,:].mean(axis=(0,2)) + var[start:start+zones_to_av,:N2//2-1:-1,:].mean(axis=(0,2))) / 2

## Internal functions ##

# Completely antisymmetric 4D symbol
# TODO cache?
def _antisym(a, b, c, d):
  # Check for valid permutation
  if (a < 0 or a > 3): return 100
  if (b < 0 or b > 3): return 100
  if (c < 0 or c > 3): return 100
  if (d < 0 or d > 3): return 100

  # Entries different? 
  if (a == b): return 0
  if (a == c): return 0
  if (a == d): return 0
  if (b == c): return 0
  if (b == d): return 0
  if (c == d): return 0

  return _pp([a,b,c,d])

# Due to Norm Hardy; good for general n
def _pp(P):
  v = np.zeros_like(P)

  p = 0
  for j in range(len(P)):
    if (v[j]):
      p += 1
    else:
      x = j
      while True:
        x = P[x]
        v[x] = 1
        if x == j:
          break

  if p % 2 == 0:
    return 1
  else:
    return -1

