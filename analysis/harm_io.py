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

  dumpfile.close()

  return dump

