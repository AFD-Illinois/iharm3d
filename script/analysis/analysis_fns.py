
# Convenience functions for physical calculations and averages
# Meant to be imported "from analysis_fns import *" for convenience

import numpy as np

THMIN = np.pi/3.
THMAX = 2.*np.pi/3.

## Physics functions ##

def Tcon(dump,i,j):
  gam = dump['hdr']['gam']
  geom = dump['geom']
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucon'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcon'][:,:,None,i,j] - dump['bcon'][:,:,:,i]*dump['bcon'][:,:,:,j] )

def Tcov(dump,i,j):
  gam = dump['hdr']['gam']
  geom = dump['geom']
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucov'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcov'][:,:,None,i,j] - dump['bcov'][:,:,:,i]*dump['bcov'][:,:,:,j] )

def Tmixed(dump,i,j):
  gam = dump['hdr']['gam']
  geom = dump['geom']
  gmixedij = (i == j)
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*gmixedij - dump['bcon'][:,:,:,i]*dump['bcov'][:,:,:,j] )

# Return mu, nu component of contravarient Maxwell tensor
# TODO there's a computationally easier way to do this:
# Pre-populate an antisym ndarray and einsum
# Same below
def Fcon(dump, i, j):
  NDIM = dump['hdr']['n_dim']
  
  Fconij = np.zeros_like(dump['RHO'])
  if i != j:
    for mu in range(NDIM):
      for nu in range(NDIM):
        Fconij[:,:,:] += _antisym(i,j,mu,nu) * dump['ucov'][:,:,:,mu] * dump['bcov'][:,:,:,nu]

  # TODO is normalization correct?
  return Fconij*dump['geom']['gdet'][:,:,None]

def Fcov(dump, i, j):
  NDIM = dump['hdr']['n_dim']
  geom = dump['geom']

  Fcovij = np.zeros_like(dump['RHO'])
  for mu in range(NDIM):
    for nu in range(NDIM):
      Fcovij += Fcon(dump, mu, nu)*geom['gcov'][:,:,None,mu,i]*geom['gcov'][:,:,None,nu,j]
  
  return Fcovij

## Sums and Averages ##

def WAVG(geom, var, w):
  return sum_shell(geom, w*var)/sum_shell(geom, w)
  
# Var must be a 3D array i.e. a grid scalar
def sum_shell(geom, var, at_zone=None):
  if at_zone is not None:

    return np.sum(var[at_zone,:,:] * geom['gdet'][at_zone,:,None]*geom['dx2']*geom['dx3'], axis=(0,1))
  else:
    return np.sum(var * geom['gdet'][:,:,None]*geom['dx2']*geom['dx3'], axis=(1,2))

def sum_vol(geom, var, within=None):
  if within is not None:
    return np.sum(var[:within,:,:] * geom['gdet'][:within,:,None]*geom['dx1']*geom['dx2']*geom['dx3'])
  else:
    return np.sum(var * geom['gdet'][:,:,None]*geom['dx1']*geom['dx2']*geom['dx3'])

# TODO can I cache the volume here without a global or object?
def eht_profile(geom, var, jmin, jmax):
  return ( np.sum(var[:,jmin:jmax,:] * geom['gdet'][:,jmin:jmax,None]*geom['dx2']*geom['dx3'], axis=(1,2)) /
           (geom['dx2']*2.*np.pi*geom['gdet'][:,:]).sum(axis=-1) )

def theta_av(var, start, av):
  # Sum theta from each pole to equator and take overall mean. N2 hack is a hack
  N2 = var.shape[1]
  return (var[start:start+av,:N2/2,:].mean(axis=-1).mean(axis=0) + var[start:start+av,:N2/2-1:-1,:].mean(axis=-1).mean(axis=0)) / 2

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

