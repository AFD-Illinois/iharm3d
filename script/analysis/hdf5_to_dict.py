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
    hdr['VERSION'] = dfile['VERSION']
    hdr['reverse'] = False
  except KeyError, e:
    hdr['reverse'] = True

  keys = ['METRIC', 'FULL_DUMP', 'ELECTRONS', 'RADIATION',
          'N1', 'N2', 'N3',
          'startx1', 'startx2', 'startx3', 'dx1', 'dx2', 'dx3',
          'gam',
          'cour',
          'DTd', 'DTl', 'DTp', 'DTr',
          'tf']
  if dfile['ELECTRONS'][0]:
    keys += ['game', 'gamp']

  try:
    if dfile['RADIATION'][0]:
      keys += ['tp_over_te']
      keys += ['L_unit', 'T_unit', 'M_unit', 'RHO_unit', 'Ne_unit', 'B_unit',
               'U_unit', 'Thetae_unit']
      keys += ['MAXNSCATT', 'NUBINS', 'numin', 'numax']
      has_radiation = True
  except KeyError:
    has_radiation = False

  if dfile['METRIC'][0] == 'MKS':
    keys += ['Rin', 'Rout', 'Reh', 'Risco', 'hslope', 'a', 'poly_xt',
             'poly_alpha', 'mks_smooth']
    if has_radiation:
      keys += ['Mbh', 'mbh']

  for key in keys:
    try:
      hdr[key] = dfile[key][0]
    except KeyError:
      hdr[key] = False

  if hdr['METRIC'] == 'MKS' and hdr['RADIATION'] == True:
    hdr['LEdd'] = 4.*np.pi*units['GNEWT']*hdr['Mbh']*units['MP']*units['CL']/units['THOMSON']
    hdr['nomEff'] = 0.1
    hdr['MdotEdd'] = hdr['LEdd']/(hdr['nomEff']*units['CL']**2)

  dfile.close()

  print "Size:", hdr['N1'], hdr['N2'], hdr['N3']
  print "Resolution:", hdr['startx1'], hdr['dx1'], hdr['startx2'], hdr['dx2'], hdr['startx3'], hdr['dx3']
  if hdr['METRIC'] == 'MKS':
    print "MKS a, hslope, poly_xt, poly_alpha, mks_smooth:", hdr['a'], hdr['hslope'], hdr['poly_xt'], hdr['poly_alpha'], hdr['mks_smooth']

  return hdr

def load_geom(hdr, fname):
  gfile = h5py.File(fname, 'r')

  # TODO add gcon,gcov to dumps
  keys = ['X1','X2','X3','gdet','X','Y','Z']

  if hdr['METRIC'] == 'MKS':
    keys += ['r','th','phi']

  geom = {}
  for key in keys:
    geom[key] = gfile[key][()]

  # these get used interchangeably
  geom['x'] = geom['X']
  geom['y'] = geom['Y']
  geom['z'] = geom['Z']
  
  # Compress gdet for normal use
  geom['gdet_full'] = geom['gdet']
  geom['gdet'] = geom['gdet'][:,:,0]

  return geom

def load_dump(fname, geom, hdr, diag=None):
  dfile = h5py.File(fname, 'r')

  keys = ['RHO', 'UU', 'U1', 'U2', 'U3', 'B1', 'B2', 'B3']
  keys += ['bsq', 'divb', 'gamma', 'fail']
  if hdr['ELECTRONS']:
    keys += ['KTOT', 'KEL']

  if hdr['FULL_DUMP']:
    if hdr['ELECTRONS']:
      keys += ['Qvisc']
      if hdr['RADIATION']:
        keys += ['Qcoul']

    if hdr['RADIATION']:
      keys += ['Rmunu', 'Nsph', 'nph', 'nuLnu']

  dump = {}
  dump['hdr'] = hdr
  for key in keys:
    if hdr['reverse']:
      dump[key] = (dfile[key][()]).transpose()
    else:
      dump[key] = dfile[key][()]
      
  dump['t'] = dfile['t'][0]

  try:
    dump['mass'] = dfile['mass'][0]
    dump['egas'] = dfile['egas'][0]
  except KeyError, e:
    pass

  if hdr['RADIATION']:
    dump['ur'] = -dfile['erad'][0]

  if hdr['ELECTRONS']:
    dump['Thetae'] = units['MP']/units['ME']*dump['KEL']*dump['RHO']**(hdr['game']-1.)
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['up'] = dump['UU'] - dump['ue']
    dump['TpTe'] = (hdr['gamp']-1.)*dump['up']/((hdr['game']-1.)*dump['ue'])
  elif hdr['RADIATION']:
    dump['Thetae'] = (hdr['gam']-1.)*units['MP']/units['ME']*(
                     1./(1. + hdr['tp_over_te'])*dump['UU']/dump['RHO'])

  N1 = hdr['N1']
  N2 = hdr['N2']
  N3 = hdr['N3']

  # This recalculates bsq
  #ucon, ucov, bcon, bcov = get_state(dump, geom)
  #dump['bsq_py'] = (bcon*bcov).sum(axis=-1)

  dump['beta'] = 2.*(hdr['gam']-1.)*dump['UU']/(dump['bsq'])
  
  if diag is not None:
    dump['mdot'] = log_time(diag, 'mdot', dump['t'])
    dump['phi'] = log_time(diag, 'phi', dump['t'])
    dump['phi_calc'] = log_time(diag, 'phi_calc', dump['t'])
  
    dump['Phi_py'] = np.sum(np.abs(dump['B1'][5,:,:]*geom['gdet'][5,:,:]*hdr['dx2']*hdr['dx3']))
    dump['phi_py'] = dump['Phi_py']/(np.sqrt(dump['mdot'])) # *2pi/sqrt(4pi) ? Just sqrt?
  
    # TODO this is not normalized or anything
    dump['Phi_disk'] = np.sum(np.abs(dump['B2'][:,N2/2,:]*geom['gdet'][:,N2/2,:]*hdr['dx1']*hdr['dx3']))
  
    err = (dump['phi_calc'] - dump['phi_py'])/dump['phi_calc']
    if err > 1e-3:
      print "Phi calculation is wrong! Error: %f" % err
  
    # Diagnostics
    #print "From Log: t: %f mdot: %f Phi_BH: %f" % (dump['t'], dump['mdot'], dump['phi'])
    #print "Calculated: phi_BH: %f Phi_disk: %f" % (dump['phi_py'], dump['Phi_disk'])

  dump.update(geom)
  dump.update(hdr)

  dfile.close()

  return dump

# For compatibility with bhlight scripts
def load_diag(path):
  return load_log(os.path.join(path, "log.out"))

def load_log(logfile):
  diag = {}
  dfile = np.loadtxt(logfile).transpose()
  
  # TODO need version strings to parse here
  if True: # New VHARM format
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
    diag['phi_calc'] = diag['Phi'] / np.sqrt(diag['mdot'])
    diag['jet_EM_flux'] = dfile[13]
  
    diag['divbmax'] = dfile[14]
  
    diag['lum_eht'] = dfile[15]
  
    diag['mdot_eh'] = dfile[16]
    diag['edot_eh'] = dfile[17]
    diag['ldot_eh'] = dfile[18]
  else: # old format
    diag['t'] = dfile[0]
    diag['rmed'] = dfile[1]
    diag['pp'] = dfile[2]
    diag['e'] = dfile[3]

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

# TODO store gcon/gcov so this can be used
def get_state(dump, geom):
  hdr = dump['hdr']
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
