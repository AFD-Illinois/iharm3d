import h5py
import numpy as np
import units
cgs = units.get_dict()

def load_header(filename):
  dumpfile = h5py.File(filename, "r")
  
  hdr = {}

  def add_var(name):
    hdr[name] = dumpfile[name][0]

  add_var('METRIC')
  add_var('ELECTRONS')
  add_var('t')
  add_var('tf')
  add_var('nstep')
  add_var('N1')
  add_var('N2')
  add_var('N3')
  add_var('startx1')
  add_var('startx2')
  add_var('startx3')
  add_var('dx1')
  add_var('dx2')
  add_var('dx3')
  if hdr['METRIC'] == 'MKS':
    add_var('Rin')
    add_var('Rout')
    add_var('Reh')
    add_var('hslope')
    add_var('a')
  add_var('gam')
  add_var('cour')
  add_var('DTd')
  add_var('DTl')
  add_var('DTr')
  add_var('DTp')
  add_var('dump_cnt')
  add_var('dt')
  add_var('failed')
  if hdr['ELECTRONS'] == True:
    add_var('game')
    add_var('gamp')

  dumpfile.close()
  
  return hdr

def load_dump(filename):
  dumpfile = h5py.File(filename, "r")
  hdr = load_header(filename)

  dump = {}
  
  def add_var(name):
    dump[name] = dumpfile[name]

  add_var('X1')
  add_var('X2')
  add_var('X3')
  if hdr['METRIC'] == 'MKS':
    add_var('r')
    add_var('theta')
    add_var('phi')
  add_var('RHO')
  add_var('UU')
  add_var('U1')
  add_var('U2')
  add_var('U3')
  add_var('B1')
  add_var('B2')
  add_var('B3')
  if hdr['ELECTRONS'] == True:
    add_var('KTOT')
    add_var('KEL')
  add_var('bsq')
  add_var('gamma')
  add_var('divb')
  add_var('fail')

  # Lop off unused dimensions
  for key in dump:
    if hdr['N2'] == 1:
      dump[key] = dump[key][:,0,0]
    elif hdr['N3'] == 1:
      dump[key] = dump[key][:,:,0]
    else:
      dump[key] = dump[key][:,:,:]
  
  # Do additional work
  if hdr['METRIC'] == 'MKS':
    if hdr['N3'] > 1:
      dump['x'] = dump['r']*np.sin(dump['theta'])*np.cos(dump['phi'])
      dump['x'][:,0] = 0
      dump['x'][:,-1] = 0
      dump['y'] = dump['r']*np.sin(dump['theta'])*np.sin(dump['phi'])
      dump['z'] = dump['r']*np.cos(dump['theta'])
    elif hdr['N2'] > 1:
      dump['x'] = dump['r']*np.sin(dump['theta'])
      dump['x'][:,0] = 0
      dump['x'][:,-1] = 0
      dump['z'] = dump['r']*np.cos(dump['theta'])
  if hdr['ELECTRONS'] == True:
    ME = cgs['ME']
    MP = cgs['MP']
    dump['Thetae'] = dump['KEL']*dump['RHO']**(hdr['game']-1.)*MP/ME
    dump['ue'] = dump['KEL']*dump['RHO']**(hdr['game'])/(hdr['game']-1.)
    dump['up'] = dump['UU'] - dump['ue']
    dump['Thetap'] = (hdr['gamp']-1.)*dump['up']/dump['RHO']
    dump['TpTe'] = (hdr['gamp']-1.)*dump['up']/((hdr['game']-1.)*dump['ue'])
    
  # Append header variables to flatten access
  dump.update(hdr)

  dumpfile.close()

  return dump

def load_log(path):
  diag = {}
  dfile = np.loadtxt(os.path.join(path, 'diag.out')).transpose()
  diag['t'] = dfile[0]
  diag['rmed'] = dfile[1]
  diag['pp'] = dfile[2]
  diag['e'] = dfile[3]
  diag['mdot'] = dfile[4]
  diag['edot'] = dfile[5]
  diag['ldot'] = dfile[6]
  diag['divbmax'] = dfile[7]
  
  return diag

## Coordinate calculations

def add_geom(dump):
  pass

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
    startx1 = hdr['startx1']
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
    trans[1,1] = np.exp(X[1])
    trans[2,1] = -np.exp(mks_smooth*(startx1-X[1]))*mks_smooth*(
      np.pi/2. -
      np.pi*X[2] +
      poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
      1./2.*(1. - hslope)*np.sin(2.*np.pi*X[2]))
    trans[2,2] = (np.pi + (1. - hslope)*np.pi*np.cos(2.*np.pi*X[2]) +
      np.exp(mks_smooth*(startx1-X[1]))*(
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

