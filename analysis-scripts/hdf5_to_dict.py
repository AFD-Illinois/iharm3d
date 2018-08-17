import sys; sys.dont_write_bytecode = True
import numpy as np
import h5py
import units
import os
import glob
units = units.get_dict()

def get_dumps_reduced(path):
  return np.sort(glob.glob(os.path.join(path, "dump*.h5")))

def get_dumps_full(folder):
  alldumps = np.sort(glob.glob(folder+'dump_*.h5'))
  fulldumps = []

  for fname in alldumps:
    dfile = h5py.File(fname, 'r')
    if dfile['is_full_dump'][()] == 1:
      fulldumps.append(fname)
    dfile.close()
  return np.sort(fulldumps)

def load_hdr(fname):
  dfile = h5py.File(fname, 'r')

  hdr = {}
  try:
    # Scoop all the keys that are not folders
    for key in [key for key in dfile['header'].keys() if not key == 'geom']:
      hdr[key] = dfile['header/' + key][()]
      
    # TODO load these from grid.h5? Or is the header actually the place for them?
    for key in [key for key in dfile['header/geom'].keys() if not key in ['mks', 'mmks'] ]:
      hdr[key] = dfile['header/geom/' + key][()]
    if 'mks' in dfile['header/geom'].keys():
      for key in dfile['header/geom/mks']:
        hdr[key] = dfile['header/geom/mks/' + key][()]
    if 'mmks' in dfile['header/geom'].keys():
      for key in dfile['header/geom/mmks']:
        hdr[key] = dfile['header/geom/mmks/' + key][()]

  except KeyError, e:
    print "File is older than supported by this library. Use hdf5_to_dict_old.py"

  dfile.close()

  # Turn the version string into components
  hdr['codename'], hdr['codestatus'], hdr['vnum'] = hdr['version'].split("-")
  hdr['vnum'] = [int(x) for x in hdr['vnum'].split(".")]

  # Work around naming bug in old versions of iharm
  if hdr['vnum'][0] <= 3 and hdr['vnum'][1] <= 3:
    names = []
    for name in hdr['prim_names'][0]:
      names.append( name )
    hdr['prim_names'] = names

  print "Loading from version ", hdr['version']
  print "Size:", hdr['n1'], hdr['n2'], hdr['n3']
  print "Resolution:", hdr['startx1'], hdr['dx1'], hdr['startx2'], hdr['dx2'], hdr['startx3'], hdr['dx3']

  return hdr

def load_geom(fname):
  gfile = h5py.File(fname, 'r')

  geom = {}
  for key in gfile['/'].keys():
    geom[key] = gfile[key][()]

  # these get used interchangeably and I don't care
  geom['x'] = geom['X']
  geom['y'] = geom['Y']
  geom['z'] = geom['Z']
  
  # Compress geom in phi for normal use
  geom['gdet_full'] = geom['gdet']
  geom['gdet'] = geom['gdet'][:,:,0]
  
  geom['gcon_full'] = geom['gcon']
  geom['gcon'] = geom['gcon'][:,:,0,:,:]
  geom['gcov_full'] = geom['gcov']
  geom['gcov'] = geom['gcov'][:,:,0,:,:]

  return geom

def load_dump(fname, geom, hdr, calc_omega=True):
  dfile = h5py.File(fname)
  
  dump = {}
  dump['hdr'] = hdr # Header is a pain to carry around separately

  # TODO this necessarily grabs the /whole/ primitives array
  for key in [key for key in dfile['/'].keys() if key not in ['header', 'extras', 'prims'] ]:
    dump[key] = dfile[key][()]

  for name, num in zip(hdr['prim_names'], range(hdr['n_prim'])):
    dump[name] = dfile['prims'][:,:,:,num]

  # Load the extras. TODO option for this
  for key in dfile['extras'].keys():
    dump[key] = dfile['extras/' + key][()]
  
  dfile.close()

  if hdr['has_electrons']:
    dump['Thetae'] = units['MP']/units['ME']*dump['KEL']*dump['RHO']**(hdr['gam_e']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['gam_e'])/(hdr['gam_e']-1.)
    dump['up'] = dump['UU'] - dump['ue'] # TODO this will fail cuz I don't output ue
    dump['TpTe'] = (hdr['gam_p']-1.)*dump['up']/((hdr['gam_e']-1.)*dump['ue'])

  # This recalculates all the derived variables
  # TODO move these somewhere less disruptive where they can be fetched on demand
  dump['ucon'], dump['ucov'], dump['bcon'], dump['bcov'] = get_state(dump, geom)
  dump['bsq'] = (dump['bcon']*dump['bcov']).sum(axis=-1)
  dump['beta'] = 2.*(hdr['gam']-1.)*dump['UU']/(dump['bsq'])
  if calc_omega:
    dump['omega'] = omega_calc(dump, geom)

  return dump

# For compatibility with bhlight scripts
def load_diag(hdr, path):
  return load_log(hdr, os.path.join(path, "log.out"))

def load_log(hdr, logfile):
  dfile = np.loadtxt(logfile).transpose()
  
  # TODO include/split on log's own header rather than tying it to dumps
  diag = {}
  diag['t'] = dfile[0]
  diag['rmed'] = dfile[1]
  diag['pp'] = dfile[2]
  diag['e'] = dfile[3]
  diag['uu_rho_gam_cent'] = dfile[4]
  diag['uu_cent'] = dfile[5]
  diag['mdot'] = dfile[6]
  diag['edot'] = dfile[7]
  diag['ldot'] = dfile[8]
  diag['mass'] = dfile[9]
  diag['egas'] = dfile[10]
  diag['Phi'] = dfile[11]
  diag['phi'] = dfile[12]
  diag['jet_EM_flux'] = dfile[13]
  diag['divbmax'] = dfile[14]
  diag['lum_eht'] = dfile[15]
  diag['mdot_eh'] = dfile[16]
  diag['edot_eh'] = dfile[17]
  diag['ldot_eh'] = dfile[18]

  return diag

def log_time(diag, var, t):
  if len(diag['t'].shape) < 1:
    return diag[var]
  else:
    i = 0
    while i < len(diag['t']) and diag['t'][i] < t:
      i += 1
    return diag[var][i-1]

# Calculate field rotation rate
# Following fns are adapted from C versions
def omega_calc(dump, geom):
  hdr = dump['hdr']
  N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3'];
  NDIM = hdr['n_dim']

  Fcov01 = np.zeros((N1,N2,N3))
  Fcov13 = np.zeros((N1,N2,N3))

  F = Fcon_calc(hdr, geom, dump)
  
  for mu in range(NDIM):
    for nu in range(NDIM):
      Fcov01 += F[:,:,:,mu,nu]*geom['gcov'][:,:,None,mu,0]*geom['gcov'][:,:,None,nu,1]
      Fcov13 += F[:,:,:,mu,nu]*geom['gcov'][:,:,None,mu,1]*geom['gcov'][:,:,None,nu,3]

  return Fcov01/Fcov13

# Return mu, nu component of contravarient Maxwell tensor
def Fcon_calc(hdr, geom, dump):
  N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3'];
  NDIM = hdr['n_dim']
  
  Fcon = np.zeros((N1,N2,N3,NDIM,NDIM))
  for mu in range(NDIM):
    for nu in range(NDIM):
      if mu == nu:
        Fcon[:,:,:,mu,nu] = 0
      else:
        for kap in range(NDIM):
          for lam in range(NDIM):
            Fcon[:,:,:,mu,nu] += antisym(mu,nu,kap,lam) * dump['ucov'][:,:,:,kap] * dump['bcov'][:,:,:,lam]

  return Fcon #*geom['gdet'][:,:,None]

# Completely antisymmetric 4D symbol
def antisym(a, b, c, d):
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

  return pp([a,b,c,d])

# Due to Norm Hardy; good for general n
def pp(P):
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

def get_state(dump, geom):
  import numpy as np
  hdr = dump['hdr']
  N1 = hdr['n1']
  N2 = hdr['n2']
  N3 = hdr['n3']

  ucon = np.zeros([N1,N2,N3,4])
  ucov = np.zeros([N1,N2,N3,4])
  bcon = np.zeros([N1,N2,N3,4])
  bcov = np.zeros([N1,N2,N3,4])

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

  ucon[:,:,:,0] = gamma/alpha
  ucon[:,:,:,1] = U1 - gamma*alpha*gcon[:,:,None,0,1]
  ucon[:,:,:,2] = U2 - gamma*alpha*gcon[:,:,None,0,2]
  ucon[:,:,:,3] = U3 - gamma*alpha*gcon[:,:,None,0,3]

  for mu in range(4):
    ucov[:,:,:,mu] = (ucon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = (B1 + bcon[:,:,:,0]*ucon[:,:,:,1])/ucon[:,:,:,0]
  bcon[:,:,:,2] = (B2 + bcon[:,:,:,0]*ucon[:,:,:,2])/ucon[:,:,:,0]
  bcon[:,:,:,3] = (B3 + bcon[:,:,:,0]*ucon[:,:,:,3])/ucon[:,:,:,0]

  for mu in range(4):
    bcov[:,:,:,mu] = (bcon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  return ucon, ucov, bcon, bcov
