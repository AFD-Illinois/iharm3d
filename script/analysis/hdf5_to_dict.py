################################################################################
#                                                                              # 
#  READ HARM OUTPUT                                                            #
#                                                                              # 
################################################################################

from __future__ import print_function, division

import os, sys
from pkg_resources import parse_version

import numpy as np
import h5py
import glob

import units

def get_dumps_list(path):
  return np.sort(glob.glob(os.path.join(path, "dump_*.h5")))

def get_full_dumps_list(folder):
  alldumps = np.sort(glob.glob(folder+'dump_*.h5'))
  fulldumps = []

  for fname in alldumps:
    dfile = h5py.File(fname, 'r')
    if dfile['is_full_dump'][()] == 1:
      fulldumps.append(fname)
    dfile.close()
  return np.sort(fulldumps)

# For single plotting scripts
def load_all(fname, **kwargs):
  hdr = load_hdr(fname)
  path = os.path.dirname(fname)
  geom = load_geom(hdr, path)
  dump = load_dump(fname, hdr, geom, **kwargs)
  return hdr, geom, dump

# Function to recursively un-bytes all the dumb HDF5 strings
def decode_all(dict):
    for key in dict:
      # Decode bytes
      if type(dict[key]) == np.bytes_:
        dict[key] = dict[key].decode('UTF-8')
      # Split ndarray of bytes into list of strings
      elif type(dict[key]) == np.ndarray:
        if dict[key].dtype.kind == 'S':
          dict[key] = [el.decode('UTF-8') for el in dict[key]]
      # Recurse for any subfolders
      elif type(dict[key]) in [list, dict]:
        decode_all(dict[key])

def load_hdr(fname):
  dfile = h5py.File(fname, 'r')

  hdr = {}
  try:
    # Scoop all the keys that are not folders
    for key in [key for key in list(dfile['header'].keys()) if not key == 'geom']:
      hdr[key] = dfile['header/' + key][()]
      
    # TODO load these from grid.h5? Or is the header actually the place for them?
    for key in [key for key in list(dfile['header/geom'].keys()) if not key in ['mks', 'mmks'] ]:
      hdr[key] = dfile['header/geom/' + key][()]
    if 'mks' in list(dfile['header/geom'].keys()):
      for key in dfile['header/geom/mks']:
        hdr[key] = dfile['header/geom/mks/' + key][()]
    if 'mmks' in list(dfile['header/geom'].keys()):
      for key in dfile['header/geom/mmks']:
        hdr[key] = dfile['header/geom/mmks/' + key][()]

  except KeyError as e:
    util.warn("File is older than supported by this library. Use hdf5_to_dict_old.py")
    exit(-1)

  decode_all(hdr)

  # Turn the version string into components
  hdr['codename'], hdr['codestatus'], hdr['vnum'] = hdr['version'].split("-")
  hdr['vnum'] = [int(x) for x in hdr['vnum'].split(".")]

  # HARM-specific workarounds:
  if hdr['codename'] == "iharm":
    # Work around naming bug before output v3.4
    if hdr['vnum'] < [3,4]:
      names = []
      for name in hdr['prim_names'][0]:
        names.append( name )
      hdr['prim_names'] = names
    
    # Work around bad radius names before output v3.6
    if hdr['vnum'] < [3,6]:
      hdr['r_in'], hdr['r_out'], hdr['r_eh'] = hdr['Rin'], hdr['Rout'], hdr['Reh']
    
    # Grab the git revision if that's something we output
    if hdr['vnum'] >= [3,6]:
      hdr['git_version'] = dfile['/extras/git_version'][()].decode('UTF-8')
  
  dfile.close()
  
  if 'git_version' in hdr.keys():
    print("Loaded header from code {}, git rev {}".format(hdr['version'], hdr['git_version']))
  else:
    print("Loaded header from code {}".format(hdr['version']))

  return hdr

def load_geom(hdr, path):
  fname = os.path.join(path, hdr['gridfile'])
  gfile = h5py.File(fname, 'r')

  geom = {}
  for key in list(gfile['/'].keys()):
    geom[key] = gfile[key][()]

  # Useful stuff for direct access in geom. TODO r_isco if available
  for key in ['n1', 'n2', 'n3', 'dx1', 'dx2', 'dx3', 'startx1', 'startx2', 'startx3', 'n_dim']:
    geom[key] = hdr[key]
  if hdr['metric'] in ["MKS", "MMKS"]:
    for key in ['r_eh', 'r_in', 'r_out']:
      geom[key] = hdr[key]

  # these get used interchangeably and I don't care
  geom['x'] = geom['X']
  geom['y'] = geom['Y']
  geom['z'] = geom['Z']

  # Compress geom in phi for normal use
  #geom['gdet_full'] = geom['gdet']
  geom['gdet'] = geom['gdet'][:,:,0]
  
  #geom['gcon_full'] = geom['gcon']
  geom['gcon'] = geom['gcon'][:,:,0,:,:]
  #geom['gcov_full'] = geom['gcov']
  geom['gcov'] = geom['gcov'][:,:,0,:,:]

  return geom

def load_dump(fname, hdr, geom, derived_vars=True, extras=True):
  dfile = h5py.File(fname, 'r')
  
  dump = {}
  
  # Carry pointers to header. Saves some pain getting shapes/parameters for plots
  # Geometry, however, _must be carried separately_ due to size in memory
  dump['hdr'] = hdr

  # TODO this necessarily grabs the /whole/ primitives array
  for key in [key for key in list(dfile['/'].keys()) if key not in ['header', 'extras', 'prims'] ]:
    dump[key] = dfile[key][()]

  for name, num in zip(hdr['prim_names'], list(range(hdr['n_prim']))):
    dump[name] = dfile['prims'][:,:,:,num]

  if extras:
    # Load the extras.
    for key in list(dfile['extras'].keys()):
      dump[key] = dfile['extras/' + key][()]
  
  dfile.close()

  # Recalculate all the derived variables, if we need to
  if derived_vars:
    dump['ucon'], dump['ucov'], dump['bcon'], dump['bcov'] = get_state(hdr, geom, dump)
    dump['bsq'] = (dump['bcon']*dump['bcov']).sum(axis=-1)
    dump['beta'] = 2.*(hdr['gam']-1.)*dump['UU']/(dump['bsq'])

    if hdr['has_electrons']:
      ref = units.get_cgs()
      dump['Thetae'] = ref['MP']/ref['ME']*dump['KEL']*dump['RHO']**(hdr['gam_e']-1.)
      dump['ue'] = dump['KEL']*dump['RHO']**(hdr['gam_e']) / (hdr['gam_e']-1.)
      dump['up'] = dump['UU'] - dump['ue']
      dump['TpTe'] = (hdr['gam_p']-1.)*dump['up']/((hdr['gam_e']-1.)*dump['ue'])

  return dump

def load_log(path):
  # TODO specify log name in dumps, like grid
  logfname = os.path.join(path,"log.out")
  dfile = np.loadtxt(logfname).transpose()
  
  # TODO log should probably have a header
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

# For adding contents of the log to dumps
def log_time(diag, var, t):
  if len(diag['t'].shape) < 1:
    return diag[var]
  else:
    i = 0
    while i < len(diag['t']) and diag['t'][i] < t:
      i += 1
    return diag[var][i-1]

# Include vectors with dumps
def get_state(hdr, geom, dump):
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

  ucon[:,:,:,0] = gamma/alpha
  ucon[:,:,:,1] = U1 - gamma*alpha*gcon[:,:,None,0,1]
  ucon[:,:,:,2] = U2 - gamma*alpha*gcon[:,:,None,0,2]
  ucon[:,:,:,3] = U3 - gamma*alpha*gcon[:,:,None,0,3]

  # TODO get the einsums working as a lower() or raise() fn
  #ucov = np.einsum("ijka,ijkba->ijkb", ucon, gcov)
  for mu in range(NDIM):
    ucov[:,:,:,mu] = (ucon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  bcon[:,:,:,0] = B1*ucov[:,:,:,1] + B2*ucov[:,:,:,2] + B3*ucov[:,:,:,3]
  bcon[:,:,:,1] = (B1 + bcon[:,:,:,0]*ucon[:,:,:,1])/ucon[:,:,:,0]
  bcon[:,:,:,2] = (B2 + bcon[:,:,:,0]*ucon[:,:,:,2])/ucon[:,:,:,0]
  bcon[:,:,:,3] = (B3 + bcon[:,:,:,0]*ucon[:,:,:,3])/ucon[:,:,:,0]

  #bcov = np.einsum("ijka,ijkba->ijkb", bcon, gcov)
  for mu in range(NDIM):
    bcov[:,:,:,mu] = (bcon[:,:,:,:]*gcov[:,:,None,mu,:]).sum(axis=-1)

  return ucon, ucov, bcon, bcov
