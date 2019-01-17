
# Convenience functions for physical calculations and averages
# Meant to be imported "from analysis_fns import *" for convenience

import numpy as np


## Physics functions ##

def T_con(geom, dump, i, j):
  gam = dump['hdr']['gam']
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucon'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcon'][:,:,None,i,j] - dump['bcon'][:,:,:,i]*dump['bcon'][:,:,:,j] )

def T_cov(geom, dump, i, j):
  gam = dump['hdr']['gam']
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucov'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcov'][:,:,None,i,j] - dump['bcov'][:,:,:,i]*dump['bcov'][:,:,:,j] )

def T_mixed(dump, i, j):
  gam = dump['hdr']['gam']
  gmixedij = (i == j)
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*gmixedij - dump['bcon'][:,:,:,i]*dump['bcov'][:,:,:,j] )

# TODO Only works for i != j
def TEM_mixed(dump, i, j):
  if i == j: raise ValueError("TEM Not implemented for i==j.")
  return dump['bsq'][:,:,:]*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j] - dump['bcon'][:,:,:,i]*dump['bcov'][:,:,:,j]

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
  return Fconij*geom['gdet_vec'][:,:,None]

def Fcov(geom, dump, i, j):
  NDIM = dump['hdr']['n_dim']

  Fcovij = np.zeros_like(dump['RHO'])
  for mu in range(NDIM):
    for nu in range(NDIM):
      Fcovij += Fcon(geom, dump, mu, nu)*geom['gcov'][:,:,None,mu,i]*geom['gcov'][:,:,None,nu,j]
  
  return Fcovij

def lower(geom, var):
  return np.einsum("...j,...ij",var,geom['gcov'][:,:,None,:,:])

# Include vectors with dumps
def get_state(hdr, geom, dump, return_gamma=False):
  N1, N2, N3, NDIM = hdr['n1'], hdr['n2'], hdr['n3'], hdr['n_dim']

  ucon = np.zeros([N1,N2,N3,NDIM])
  ucov = np.zeros_like(ucon)
  bcon = np.zeros_like(ucon)
  bcov = np.zeros_like(ucon)

  # Aliases to make the below more readable
  gcov = geom['gcov']
  gcon = geom['gcon']

  U1 = dump['U1']
  U2 = dump['U2']
  U3 = dump['U3']
  B1 = dump['B1']
  B2 = dump['B2']
  B3 = dump['B3']

  alpha = geom['lapse']
  qsq = (gcov[:,:,None,1,1]*U1**2 + gcov[:,:,None,2,2]*U2**2 +
         gcov[:,:,None,3,3]*U3**2 + 2.*(gcov[:,:,None,1,2]*U1*U2 +
                                        gcov[:,:,None,1,3]*U1*U3 +
                                        gcov[:,:,None,2,3]*U2*U3))
  gamma = np.sqrt(1. + qsq)

  ucon[:,:,:,0] = gamma/(alpha[:,:,None])
  ucon[:,:,:,1] = U1 - gamma*alpha[:,:,None]*gcon[:,:,None,0,1]
  ucon[:,:,:,2] = U2 - gamma*alpha[:,:,None]*gcon[:,:,None,0,2]
  ucon[:,:,:,3] = U3 - gamma*alpha[:,:,None]*gcon[:,:,None,0,3]

  ucov = lower(geom, ucon)

  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = (B1 + bcon[:,:,:,0]*ucon[:,:,:,1])/ucon[:,:,:,0]
  bcon[:,:,:,2] = (B2 + bcon[:,:,:,0]*ucon[:,:,:,2])/ucon[:,:,:,0]
  bcon[:,:,:,3] = (B3 + bcon[:,:,:,0]*ucon[:,:,:,3])/ucon[:,:,:,0]

  bcov = lower(geom, bcon)

  if geom['mixed_metrics']:
    to_grid = lambda vec : np.einsum("...i,...ij->...j", vec, geom['vec_to_grid'][:,:,None,:,:])
    ucon = to_grid(ucon)
    ucov = to_grid(ucov)
    bcon = to_grid(bcon)
    bcov = to_grid(bcov)

  ret = (ucon, ucov, bcon, bcov)
  if return_gamma:
    ret += gamma

  return ret

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

def theta_av(geom, var, start, av):
  # Sum theta from each pole to equator and take overall mean
  N2 = geom['n2']
  return (var[start:start+av,:N2//2,:].mean(axis=(0,2)) + var[start:start+av,:N2//2-1:-1,:].mean(axis=(0,2))) / 2

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

