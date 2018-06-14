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
    if dfile['FULL_DUMP'][0]:
      fulldumps.append(fname)
  return np.sort(fulldumps)

def load_hdr(fname):
  dfile = h5py.File(fname, 'r')

  # Every version with a string outputs forward
  hdr = {}
  try:
    hdr['VERSION'] = dfile['VERSION'][0]
    hdr['reverse'] = False
    if len(hdr['VERSION'].split("-")) > 1:
      hdr['version_number'] = float(hdr['VERSION'].split("-")[-1])
    else:
      hdr['version_number'] = 0.05
  except KeyError, e:
    try:
      hdr['version'] = dfile['version_string'][()]
      hdr['reverse'] = False
      hdr['version_number'] = float(hdr['version'].split("-")[-1])
    except KeyError, e:
      hdr['VERSION'] = None
      hdr['version_number'] = 0.0
      hdr['reverse'] = True

  print "Loading from VHARM version %f" % hdr['version_number']

  if hdr['version_number'] > 3.0:
    keys = ['metric_name', 'has_electrons', #'has_radiation', 'is_full', 'n_dim'
            'n1', 'n2', 'n3',
            'startx1', 'startx2', 'startx3', 'dx1', 'dx2', 'dx3',
            'gam', 'cour', 'tf']      
    for key in keys:
      hdr[key] = dfile[key][()]
    # Fill things not output
    hdr['is_full'] = False
    hdr['has_radiation'] = False
    hdr['n_dim'] = 4
    # Compat
    hdr['N1'] = hdr['n1']; hdr['N2'] = hdr['n2']; hdr['N3'] = hdr['n3']
    hdr['RADIATION'] = hdr['has_radiation']
    hdr['METRIC'] = hdr['metric_name']
  else:
    keys = ['METRIC', 'FULL_DUMP', 'ELECTRONS', 'RADIATION',
            'N1', 'N2', 'N3',
            'startx1', 'startx2', 'startx3', 'dx1', 'dx2', 'dx3',
            'gam',
            'cour',
            'DTd', 'DTl', 'DTp', 'DTr',
            'tf']
    for key in keys:
      try:
        hdr[key] = dfile[key][0]
      except KeyError:
        hdr[key] = False
    # Fill things not output
    hdr['is_full'] = False
    hdr['has_radiation'] = False
    hdr['n_dim'] = 4
    # Compat
    hdr['has_electrons'] = hdr['ELECTRONS']
    hdr['has_radiation'] = hdr['RADIATION']
    hdr['metric_name'] = hdr['METRIC']
    hdr['n1'] = hdr['N1']; hdr['n2'] = hdr['N2']; hdr['n3'] = hdr['N3']
      
  # Set secondary keys based on properties
  keys = []
  if hdr['has_electrons']:
    keys += ['gam_e', 'gam_p']

  if hdr['has_radiation']:
    keys += ['tp_over_te']
    keys += ['L_unit', 'T_unit', 'M_unit', 'RHO_unit', 'Ne_unit', 'B_unit',
             'U_unit', 'Thetae_unit']
    keys += ['MAXNSCATT', 'NUBINS', 'numin', 'numax']

  if hdr['metric_name'] == 'MKS':
    keys += ['Rin', 'Rout', 'Reh', 'hslope', 'a'] # 'Risco',
    if False: # TODO output POLYTH
      keys += ['poly_xt', 'poly_alpha', 'mks_smooth']
    if hdr['has_radiation']:
      keys += ['Mbh', 'mbh']
  else:
    pass #TODO load Minkowski x1min/max etc here

  # Load secondary keys
  if hdr['version_number'] > 3.0:
    for key in keys:
      hdr[key] = dfile[key][()]
  else:
    for key in keys:
      hdr[key] = dfile[key][0]

  if hdr['has_radiation']:
    hdr['LEdd'] = 4.*np.pi*units['GNEWT']*hdr['Mbh']*units['MP']*units['CL']/units['THOMSON']
    hdr['nomEff'] = 0.1
    hdr['MdotEdd'] = hdr['LEdd']/(hdr['nomEff']*units['CL']**2)

  dfile.close()

  print "Size:", hdr['n1'], hdr['n2'], hdr['n3']
  print "Resolution:", hdr['startx1'], hdr['dx1'], hdr['startx2'], hdr['dx2'], hdr['startx3'], hdr['dx3']

  return hdr

def load_geom(hdr, fname):
  gfile = h5py.File(fname, 'r')

  keys = ['X1','X2','X3','X','Y','Z', 'gdet']

  if hdr['version_number'] > 0.05:
    keys += ['gcon', 'gcov']

  if hdr['metric_name'] == 'MKS':
    keys += ['r','th','phi']

  geom = {}
  for key in keys:
    geom[key] = gfile[key][()]

  # these get used interchangeably
  geom['x'] = geom['X']
  geom['y'] = geom['Y']
  geom['z'] = geom['Z']
  
  # Compress geom in phi for normal use
  geom['gdet_full'] = geom['gdet']
  geom['gdet'] = geom['gdet'][:,:,0]
  
  if hdr['version_number'] > 0.05:
    geom['gcon_full'] = geom['gcon']
    geom['gcon'] = geom['gcon'][:,:,0,:,:]
    geom['gcov_full'] = geom['gcov']
    geom['gcov'] = geom['gcov'][:,:,0,:,:]

  return geom

def load_dump(fname, geom, hdr, diag=None):
  dfile = h5py.File(fname)
  dump = {}
  
  # "Header" variables that change per dump
  # Be accomodating about these as we rarely need them
  try:
    dump['t'] = dfile['t'][0]
    dump['mass'] = dfile['mass'][0]
    dump['egas'] = dfile['egas'][0]
  except (KeyError, ValueError) as e:
    try:
      dump['t'] = dfile['t'][()]
      dump['mass'] = dfile['mass'][()]
      dump['egas'] = dfile['egas'][()]
    except KeyError, e:
      pass

  # Usual primitive variable names
  varnames = ['RHO', 'UU', 'U1', 'U2', 'U3', 'B1', 'B2', 'B3']
  if hdr['has_electrons']:
    varnames += ['KTOT', 'KEL']
  
  if hdr['version_number'] > 3.0:
    keys = ['prims']
  else:
    keys = varnames
  
  keys += ['bsq', 'divb', 'gamma', 'fail']

  if hdr['is_full']:
    if hdr['has_electrons']:
      keys += ['Qvisc']
      if hdr['has_radiation']:
        keys += ['Qcoul']

    if hdr['has_radiation']:
      keys += ['Rmunu', 'Nsph', 'nph', 'nuLnu']

  for key in keys:
    if hdr['reverse']:
      dump[key] = (dfile[key][()]).transpose()
    else:
      dump[key] = dfile[key][()]
      
  if hdr['version_number'] > 3.0:
    for name, num in zip(varnames, range(len(varnames))):
      dump[name] = dump['prims'][:,:,:,num]

  # Not all VHARM/bhlight output all variables
  # These are the optional ones so we don't mind
  ext_keys = ['bcon', 'bcov', 'ucon', 'ucov', 'jcon']

  for key in ext_keys:
    try:
      dump[key] = dfile[key][()]
    except KeyError, e:
      pass

  if hdr['has_radiation']:
    dump['ur'] = -dfile['erad'][0]

  if hdr['has_electrons']:
    dump['Thetae'] = units['MP']/units['ME']*dump['KEL']*dump['RHO']**(hdr['game']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['up'] = dump['UU'] - dump['ue']
    dump['TpTe'] = (hdr['gamp']-1.)*dump['up']/((hdr['game']-1.)*dump['ue'])
  elif hdr['has_radiation']:
    dump['Thetae'] = (hdr['gam']-1.)*units['MP']/units['ME']*(
                     1./(1. + hdr['tp_over_te'])*dump['UU']/dump['RHO'])

  N1 = hdr['n1']
  N2 = hdr['n2']
  N3 = hdr['n3']

  # This recalculates bsq
  #ucon, ucov, bcon, bcov = get_state(dump, geom)
  #dump['bsq_py'] = (bcon*bcov).sum(axis=-1)

  dump['beta'] = 2.*(hdr['gam']-1.)*dump['UU']/(dump['bsq'])

  # TODO move this calculation somewhere less disruptive (analysis.py?)
  dump['omega'] = omega_calc(hdr, geom, dump)
  
  dfile.close()

  dump['hdr'] = hdr

  return dump

# For compatibility with bhlight scripts
def load_diag(hdr, path):
  return load_log(hdr, os.path.join(path, "log.out"))

def load_log(hdr, logfile):
  diag = {}
  dfile = np.loadtxt(logfile).transpose()
  
  diag['t'] = dfile[0]
  diag['rmed'] = dfile[1]
  diag['pp'] = dfile[2]
  diag['e'] = dfile[3]
  
  # TODO find differentiator for /very/ old versions?
  if hdr['version_number'] >= 0.0:
  
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
    
    if hdr['version_number'] >= 0.1:
      diag['mass_added'] = dfile[16]
  
      diag['mdot_eh'] = dfile[17]
      diag['edot_eh'] = dfile[18]
      diag['ldot_eh'] = dfile[19]
    else:
      diag['mdot_eh'] = dfile[16]
      diag['edot_eh'] = dfile[17]
      diag['ldot_eh'] = dfile[18]
  else: # old format
    
    diag['mdot'] = dfile[4]
    diag['edot'] = dfile[5]
    diag['ldot'] = dfile[6]

    diag['phi'] = dfile[7]

    diag['divbmax'] = dfile[8]

  return diag

def log_time(diag, var, t):
  i = 0
  if len(diag['t'].shape) < 1:
    return diag[var]
  else:
    while i < len(diag['t']) and diag['t'][i] < t:
      i += 1
    return diag[var][i-1]

# Calculate field rotation rate
# Following fns are adapted from C versions
def omega_calc(hdr, geom, dump):
  N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3'];
  NDIM = hdr['n_dim']

  Fcov01 = np.zeros((N1,N2,N3))
  Fcov13 = np.zeros((N1,N2,N3))
  for mu in range(NDIM):
    for nu in range(NDIM):
      Fmunu = Fcon_calc(hdr, geom, dump, mu, nu);
      Fcov01 += Fmunu*geom['gcov'][:,:,None,mu,0]*geom['gcov'][:,:,None,nu,1]
      Fcov13 += Fmunu*geom['gcov'][:,:,None,mu,1]*geom['gcov'][:,:,None,nu,3]

    return Fcov01/Fcov13

# Return mu, nu component of contravarient Maxwell tensor at grid zone i, j, k
def Fcon_calc(hdr, geom, dump, mu, nu):
  N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3'];
  NDIM = hdr['n_dim']
  
  if (mu == nu):
    return np.zeros((N1,N2,N3))
  
  Fcon = np.zeros((N1,N2,N3))
  for kap in range(NDIM):
    for lam in range(NDIM):
      Fcon[:,:,:] += (-1./geom['gdet'][:,:,None]) * \
      antisym(mu,nu,kap,lam) * dump['ucov'][:,:,:,kap] * dump['bcov'][:,:,:,lam]

  return Fcon*geom['gdet'][:,:,None]

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

# TODO store gcon/gcov so this can be used
def get_state(hdr, geom, dump):
  N1 = hdr['N1']
  N2 = hdr['N2']
  N3 = hdr['N3']

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

# TODO vectorize this
  alpha = np.zeros([N1,N2,N3])
  for i in xrange(N1):
    for j in xrange(N2):
      alpha[i,j,:] = 1./np.sqrt(-gcon[i,j,0,0])
  qsq = np.zeros([N1,N2,N3])
  for i in xrange(N1):
    for j in xrange(N2):
      for k in xrange(N3):
        qsq[i,j,k] = (gcov[i,j,1,1]*U1[i,j,k]**2 + gcov[i,j,2,2]*U2[i,j,k]**2 +
                      gcov[i,j,3,3]*U3[i,j,k]**2 + 2.*(gcov[i,j,1,2]*U1[i,j,k]*U2[i,j,k] +
                                                gcov[i,j,1,3]*U1[i,j,k]*U3[i,j,k] + 
                                                gcov[i,j,2,3]*U2[i,j,k]*U3[i,j,k]))
  gamma = np.sqrt(1. + qsq)

  ucon[:,:,:,0] = gamma/alpha
  ucon[:,:,:,1] = U1 - gamma*alpha*gcon[:,:,None,0,1]
  ucon[:,:,:,2] = U2 - gamma*alpha*gcon[:,:,None,0,2]
  ucon[:,:,:,3] = U3 - gamma*alpha*gcon[:,:,None,0,3]

  for mu in xrange(4):
    ucov[:,:,:,mu] = (ucon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = B1 + bcon[:,:,:,0]*ucon[:,:,:,1]/ucon[:,:,:,0]
  bcon[:,:,:,2] = B2 + bcon[:,:,:,0]*ucon[:,:,:,2]/ucon[:,:,:,0]
  bcon[:,:,:,3] = B3 + bcon[:,:,:,0]*ucon[:,:,:,3]/ucon[:,:,:,0]

  for mu in xrange(4):
    bcov[:,:,:,mu] = (bcon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  bcon[:,:,:,0] = dump['B1'][:,:,:]*ucov[:,:,:,1]

  return ucon, ucov, bcon, bcov
