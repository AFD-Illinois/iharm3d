import sys; sys.dont_write_bytecode = True
import numpy as np
import h5py
import units
import os
import glob
units = units.get_dict()

def get_dumps_reduced(folder):
  return np.sort(glob.glob(folder+'dump_*.h5'))

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

  hdr = {}
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

  return hdr

def load_geom(hdr, dump):
  N1 = hdr['N1']
  N2 = hdr['N2']
  N3 = hdr['N3']

  r = np.zeros([N1, N2, N3])
  th = np.zeros([N1, N2, N3])
  phi = np.zeros([N1, N2, N3])
  #gcov = np.zeros([N1, N2, 4, 4])
  #gcon = np.zeros([N1, N2, 4, 4])
  #gdet = np.zeros([N1, N2, N3])

  # Old way of getting coords
  # for i in xrange(N1):
  #   for j in xrange(N2):
  #     for k in xrange(N3):
  #       if hdr['METRIC'] == 'MKS':
  #         X = harm_coord(hdr, i, j, k)
  #         r[i,j,k], th[i,j,k], phi[i,j,k] = bl_coord(hdr, X)
  #         gdet[i,j,k] = np.sqrt(-np.linalg.det(get_gcov(hdr, X)))
  #     gcov[i,j,:,:] = get_gcov(hdr, X)
  #     gcon[i,j,:,:] = get_gcon(gcov[i,j,:,:])

  # New way
  r[:,:,:]= dump['r'][:,:,:]
  th[:,:,:] = dump['theta'][:,:,:]
  phi[:,:,:] = dump['phi'][:,:,:]

  #gcon[:,:,:] = dump['gcon'][:,:,:]
  #gcov[:,:,:] = dump['gcov'][:,:,:]
  #gdet[:,:,:] = dump['gdet'][:,:,:]

  if hdr['METRIC'] == 'MKS':
    if hdr['N3'] == 1:
      phi = np.zeros([N1, N2, N3])

    x = r*np.sin(th)*np.cos(phi)
    y = r*np.sin(th)*np.sin(phi)
    z = r*np.cos(th)

    if hdr['N3'] == 1:
      x[:,0,:] = 0.
      x[:,-1,:] = 0.

  geom = {}
  geom['r'] = r
  geom['th'] = th
  geom['phi'] = phi
  geom['x'] = x
  geom['y'] = y
  geom['z'] = z
  #geom['gcov'] = gcov
  #geom['gcon'] = gcon
  #geom['gdet'] = gdet

  return geom

def load_dump(fname):
  hdr = load_hdr(fname)

  dfile = h5py.File(fname, 'r')

  geom = load_geom(hdr, dfile)

  keys = ['RHO', 'UU', 'U1', 'U2', 'U3', 'B1', 'B2', 'B3']
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
    dump[key] = dfile[key][()]
  dump['t'] = dfile['t'][0]

  try:
    dump['mass'] = dfile['mass'][0]
    dump['egas'] = dfile['egas'][0]
  except KeyError, e:
    print e

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

  #ucon, ucov, bcon, bcov = get_state(dump, geom)

  #dump['bsq'] = (bcon*bcov).sum(axis=-1)

  #dump['beta'] = 2.*(hdr['gam']-1.)*dump['UU']/(dump['bsq'])

  dump.update(geom)
  dump.update(hdr)

  dfile.close()

  return dump

# TODO TODO match this with the real logfile
def load_log(path):
  diag = {}
  dfile = np.loadtxt(os.path.join(path, 'log.out')).transpose()
  diag['t'] = dfile[0]
  diag['rmed'] = dfile[1]
  diag['pp'] = dfile[2]
  diag['e'] = dfile[3]
  diag['mdot'] = dfile[4]
  diag['edot'] = dfile[5]
  diag['ldot'] = dfile[6]
  diag['divbmax'] = dfile[7]

  return diag

def load_diag(path):
  diag = {}
  dfile = np.loadtxt(os.path.join(path, 'diag.out')).transpose()
  diag['t'] = dfile[0]
  diag['mdot'] = dfile[6]
  diag['edot'] = dfile[7]
  diag['ldot'] = dfile[8]
  diag['mass'] = dfile[9]
  diag['egas'] = dfile[10]
  diag['Phi'] = dfile[11]
  diag['phi'] = dfile[12]
  diag['jet_EM_flux'] = dfile[13]
  diag['divbmax'] = dfile[14]
  if len(dfile) > 15:
    diag['step_made']  = dfile[12]
    diag['step_abs']   = dfile[13]
    diag['step_scatt'] = dfile[14]
    diag['step_lost']  = dfile[15]
    diag['step_rec']   = dfile[16]
    diag['step_tot']   = dfile[17]
    diag['step_sent']  = dfile[18]
    diag['step_rcvd']  = dfile[19]
    diag['step_made_all']  = dfile[20]
    diag['step_abs_all']   = dfile[21]
    diag['step_scatt_all'] = dfile[22]
    diag['step_lost_all']  = dfile[23]
    diag['step_rec_all']   = dfile[24]
    diag['step_tot_all']   = dfile[25]
    diag['step_sent_all']  = dfile[26]
    diag['step_rcvd_all']  = dfile[27]
    diag['tune_emiss']  = dfile[28]
    diag['tune_scatt']  = dfile[29]
    diag['erad'] = dfile[30]
    diag['lum'] = dfile[31]
    diag['eff'] = dfile[32]

  return diag

def harm_coord(hdr, i, j, k):
  startx = [0, hdr['startx1'], hdr['startx2'], hdr['startx3']]
  dx = [0, hdr['dx1'], hdr['dx2'], hdr['dx3']]
  N1 = hdr['N1']
  N2 = hdr['N2']
  N3 = hdr['N3']
  X = np.zeros(4)

  X[1] = startx[1] + (i + 0.5)*dx[1]
  X[2] = startx[2] + (j + 0.5)*dx[2]
  X[3] = startx[3] + (k + 0.5)*dx[3]

  return X

def bl_coord(hdr, X):
  hslope = hdr['hslope']
  poly_xt = hdr['poly_xt']
  poly_alpha = hdr['poly_alpha']
  mks_smooth = hdr['mks_smooth']
  poly_norm = 0.5*np.pi*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));
  startx3 = hdr['startx3']

  # Switched 1,3 due to coord change in harm
  r = np.exp(X[3])
  y = 2*X[2] - 1
  thG = np.pi*X[2] + ((1. - hslope)/2.)*np.sin(2.*np.pi*X[2]);
  thJ = poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*np.pi;
  th = thG + np.exp(mks_smooth*(startx3 - X[3]))*(thJ - thG);

  phi = X[1]

  return r, th, phi

def get_gcov(hdr, X):
  gcov = np.zeros([4,4])
  if hdr['METRIC'] == 'MINKOWSKI':
    gcov[0,0] = -1
    gcov[1,1] = 1
    gcov[2,2] = 1
    gcov[3,3] = 1
  elif hdr['METRIC'] == 'MKS':
    a = hdr['a']
    hslope = hdr['hslope']
    poly_alpha = hdr['poly_alpha']
    poly_xt = hdr['poly_xt']
    mks_smooth = hdr['mks_smooth']
    startx3 = hdr['startx3']
    poly_norm = 0.5*np.pi*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));

    r, th, phi, = bl_coord(hdr, X)
    cth = np.cos(th)
    sth = np.sin(th)
    s2 = sth**2
    rho2 = r**2 + a**2*cth**2

    # KS coordinates
    gcovKS = np.zeros([4,4])
    gcovKS[0,0] = -1. + 2.*r/rho2
    gcovKS[0,1] = 2.*r/rho2
    gcovKS[0,3] = -2.*a*r*s2/rho2

    gcovKS[1,0] = gcovKS[0,1]
    gcovKS[1,1] = 1. + 2.*r/rho2
    gcovKS[1,3] = -a*s2*(1. + 2.*r/rho2)

    gcovKS[2,2] = rho2

    gcovKS[3,0] = gcovKS[0,3]
    gcovKS[3,1] = gcovKS[1,3]
    gcovKS[3,3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2))

    # Convert from KS to MKS
    trans = np.zeros([4,4])

    trans[0,0] = 1.
    trans[1,1] = np.exp(X[3])
    trans[2,1] = -np.exp(mks_smooth*(startx3-X[3]))*mks_smooth*(
      np.pi/2. -
      np.pi*X[2] +
      poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
      1./2.*(1. - hslope)*np.sin(2.*np.pi*X[2]))
    trans[2,2] = (np.pi + (1. - hslope)*np.pi*np.cos(2.*np.pi*X[2]) +
      np.exp(mks_smooth*(startx3-X[3]))*(
        -np.pi +
        2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
        (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)    *poly_xt) -
        (1.-hslope)*np.pi*np.cos(2.*np.pi*X[2])))
    trans[3,3] = 1.

    for mu in xrange(4):
      for nu in xrange(4):
        for lam in xrange(4):
          for kap in xrange(4):
            gcov[mu,nu] += gcovKS[lam,kap]*trans[lam,mu]*trans[kap,nu]

  return gcov

def get_gcon(gcov):
  return np.linalg.inv(gcov)

def lower(vcon, gcov):
  vcov = np.zeros(4)

  for mu in xrange(4):
    for nu in xrange(4):
      vcov[mu] += gcov[mu,nu]*vcon[nu]

  return vcov

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

  alpha = np.zeros([N1,N2,N3])
  alpha = 1./np.sqrt(-gcon[:,:,None,0,0])
  qsq = np.zeros([N1,N2,N3])
  qsq = (gcov[:,:,None,1,1]*U1**2 + gcov[:,:,None,2,2]*U2**2 +
         gcov[:,:,None,3,3]*U3**2 + 2.*(gcov[:,:,None,1,2]*U1*U2 +
                                        gcov[:,:,None,1,3]*U1*U3 +
                                        gcov[:,:,None,2,3]*U2*U3))
  gamma = np.ones([N1,N2,N3])
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
