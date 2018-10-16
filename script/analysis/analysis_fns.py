
# Convenience functions for physical calculations and averages
# I got tired of mis-typing components and sums

import numpy as np

jmin, jmax = 0,0
vol_profile = 0

THMIN = np.pi/3.
THMAX = 2.*np.pi/3.

def init_analysis(geom):
  # Setting module-wide variables
  global jmin, jmax, vol_profile
  
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

  # Constant volume for profiles
  vol_profile = (geom['dx2']*2.*np.pi*geom['gdet'][:,:]).sum(axis=-1)

def Tcon(geom, dump,i,j):
  gam = dump['hdr']['gam']
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucon'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcon'][:,:,None,i,j] - dump['bcon'][:,:,:,i]*dump['bcon'][:,:,:,j] )

def Tcov(geom, dump,i,j):
  gam = dump['hdr']['gam']
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucov'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*geom['gcov'][:,:,None,i,j] - dump['bcov'][:,:,:,i]*dump['bcov'][:,:,:,j] )

def Tmixed(geom, dump,i,j):
  gam = dump['hdr']['gam']
  gmixedij = np.sum(geom['gcon'][:,:,None,i,:]*geom['gcov'][:,:,None,:,j],axis=-1)
  return ( (dump['RHO'] + dump['UU'] + (gam-1)*dump['UU'] + dump['bsq'])*dump['ucon'][:,:,:,i]*dump['ucov'][:,:,:,j] +
           ((gam-1)*dump['UU'] + dump['bsq']/2)*gmixedij - dump['bcon'][:,:,:,i]*dump['bcov'][:,:,:,j] )

# Return mu, nu component of contravarient Maxwell tensor
def Fcon(geom, dump, i, j):
  NDIM = dump['hdr']['n_dim']
  
  Fconij = np.zeros_like(dump['RHO'])
  if i != j:
    for mu in range(NDIM):
      for nu in range(NDIM):
        Fconij[:,:,:] += _antisym(i,j,mu,nu) * dump['ucov'][:,:,:,mu] * dump['bcov'][:,:,:,nu]

  # TODO is normalization correct?
  return Fconij*geom['gdet'][:,:,None]

# TODO there's a computationally easier way to do this
def Fcov(geom, dump, i, j):
  NDIM = dump['hdr']['n_dim']
  
  Fcovij = np.zeros_like(dump['RHO'])
  for mu in range(NDIM):
    for nu in range(NDIM):
      Fcovij += Fcon(geom, dump, mu, nu)*geom['gcov'][:,:,None,mu,i]*geom['gcov'][:,:,None,nu,j]
  
  return Fcovij

# Completely antisymmetric 4D symbol
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

## SUMS AND INTEGRALS

def WAVG(geom, var, w):
  return sum_shell(geom, w*var)/sum_shell(geom, w)
  
  # Var must be a 3D array i.e. a grid scalar
def sum_shell(geom, var):
  return np.sum(var * geom['gdet'][:,:,None]*geom['dx2']*geom['dx3'], axis=(1,2))

def sum_shell_at(geom, var, i):
  return np.sum(var[i,:,:] * geom['gdet'][i,:,None]*geom['dx2']*geom['dx3'], axis=(0,1))

def eht_profile(geom, var):
  return np.sum(var[:,jmin:jmax,:] * geom['gdet'][:,jmin:jmax,None]*geom['dx2']*geom['dx3'], axis=(1,2)) / vol_profile

def theta_av(var, start, av):
  # Sum theta from each pole to equator and take overall mean. N2 hack is a hack
  N2 = var.shape[1]
  return (var[start:start+av,:N2/2,:].mean(axis=-1).mean(axis=0) + var[start:start+av,:N2/2-1:-1,:].mean(axis=-1).mean(axis=0)) / 2
