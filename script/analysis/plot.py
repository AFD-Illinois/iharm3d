################################################################################
#                                                                              #
#  UTILITIES FOR PLOTTING                                                      #
#                                                                              #
################################################################################

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from matplotlib.colors import LinearSegmentedColormap
cm_data = [[0.2081, 0.1663, 0.5292], [0.2116238095, 0.1897809524, 0.5776761905],
 [0.212252381, 0.2137714286, 0.6269714286], [0.2081, 0.2386, 0.6770857143],
 [0.1959047619, 0.2644571429, 0.7279], [0.1707285714, 0.2919380952,
  0.779247619], [0.1252714286, 0.3242428571, 0.8302714286],
 [0.0591333333, 0.3598333333, 0.8683333333], [0.0116952381, 0.3875095238,
  0.8819571429], [0.0059571429, 0.4086142857, 0.8828428571],
 [0.0165142857, 0.4266, 0.8786333333], [0.032852381, 0.4430428571,
  0.8719571429], [0.0498142857, 0.4585714286, 0.8640571429],
 [0.0629333333, 0.4736904762, 0.8554380952], [0.0722666667, 0.4886666667,
  0.8467], [0.0779428571, 0.5039857143, 0.8383714286],
 [0.079347619, 0.5200238095, 0.8311809524], [0.0749428571, 0.5375428571,
  0.8262714286], [0.0640571429, 0.5569857143, 0.8239571429],
 [0.0487714286, 0.5772238095, 0.8228285714], [0.0343428571, 0.5965809524,
  0.819852381], [0.0265, 0.6137, 0.8135], [0.0238904762, 0.6286619048,
  0.8037619048], [0.0230904762, 0.6417857143, 0.7912666667],
 [0.0227714286, 0.6534857143, 0.7767571429], [0.0266619048, 0.6641952381,
  0.7607190476], [0.0383714286, 0.6742714286, 0.743552381],
 [0.0589714286, 0.6837571429, 0.7253857143],
 [0.0843, 0.6928333333, 0.7061666667], [0.1132952381, 0.7015, 0.6858571429],
 [0.1452714286, 0.7097571429, 0.6646285714], [0.1801333333, 0.7176571429,
  0.6424333333], [0.2178285714, 0.7250428571, 0.6192619048],
 [0.2586428571, 0.7317142857, 0.5954285714], [0.3021714286, 0.7376047619,
  0.5711857143], [0.3481666667, 0.7424333333, 0.5472666667],
 [0.3952571429, 0.7459, 0.5244428571], [0.4420095238, 0.7480809524,
  0.5033142857], [0.4871238095, 0.7490619048, 0.4839761905],
 [0.5300285714, 0.7491142857, 0.4661142857], [0.5708571429, 0.7485190476,
  0.4493904762], [0.609852381, 0.7473142857, 0.4336857143],
 [0.6473, 0.7456, 0.4188], [0.6834190476, 0.7434761905, 0.4044333333],
 [0.7184095238, 0.7411333333, 0.3904761905],
 [0.7524857143, 0.7384, 0.3768142857], [0.7858428571, 0.7355666667,
  0.3632714286], [0.8185047619, 0.7327333333, 0.3497904762],
 [0.8506571429, 0.7299, 0.3360285714], [0.8824333333, 0.7274333333, 0.3217],
 [0.9139333333, 0.7257857143, 0.3062761905], [0.9449571429, 0.7261142857,
  0.2886428571], [0.9738952381, 0.7313952381, 0.266647619],
 [0.9937714286, 0.7454571429, 0.240347619], [0.9990428571, 0.7653142857,
  0.2164142857], [0.9955333333, 0.7860571429, 0.196652381],
 [0.988, 0.8066, 0.1793666667], [0.9788571429, 0.8271428571, 0.1633142857],
 [0.9697, 0.8481380952, 0.147452381], [0.9625857143, 0.8705142857, 0.1309],
 [0.9588714286, 0.8949, 0.1132428571], [0.9598238095, 0.9218333333,
  0.0948380952], [0.9661, 0.9514428571, 0.0755333333],
 [0.9763, 0.9831, 0.0538]]
parula = LinearSegmentedColormap.from_list('parula', cm_data)

# GET XZ SLICE OF GRID DATA
def flatten_xz(array, hdr, flip=False):
    if flip:
        sign = -1
    else:
        sign = 1
    sign = 1.
    flat = np.zeros([2*hdr['N1'],hdr['N2']])
    for j in xrange(hdr['N2']):
        for i in xrange(hdr['N1']):
            flat[i,j] = sign*array[hdr['N1'] - 1 - i,j,hdr['N3']/2]
        for i in xrange(hdr['N1']):
            flat[i+hdr['N1'],j] = array[i,j,0]
    if flip:
        flat[:,0] = 0
        flat[:,-1] = 0
    return flat

# GET XY SLICE OF GRID DATA
def flatten_xy(array, hdr):
  return np.vstack((array.transpose(),array.transpose()[0])).transpose()

def plot_xz(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
            label=None, ticks=None, arrayspace=False):
  hdr = dump['hdr']; N1 = hdr['N1']; N2 = hdr['N2']; N3 = hdr['N3']
  if N3 > 1.:
    if (arrayspace):
      x = np.reshape(np.repeat(np.linspace(0,1,N2),N1),(N1,N2))
      z = np.transpose(np.reshape(np.repeat(np.linspace(0,1,N1),N2),(N2,N1)))
      var = flatten_xz(var, hdr)[N2:,:]
    else:
      x = geom['x']
      z = geom['z']
      x = flatten_xz(x, hdr, flip=True)
      z = flatten_xz(z, hdr)
      var = flatten_xz(var, hdr)

  else:
    x = x[:,:,0]
    z = z[:,:,0]
    var = var[:,:,0]

#  print 'xshape is ', x.shape, ', zshape is ', z.shape, ', varshape is ', var.shape
  mesh = ax.pcolormesh(x, z, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading='gouraud')

  if (not arrayspace):
    circle1=plt.Circle((0,0),hdr['Reh'],color='k');
    ax.add_artist(circle1)

  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, ticks=ticks).set_label(label)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  if (arrayspace):
    ax.set_xlabel('r (arbitrary)'); ax.set_ylabel('theta (arbitrary)')
  else:
    ax.set_xlabel('x/M'); ax.set_ylabel('z/M')
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)

def overlay_field(ax, geom, dump, NLEV):
  from scipy.integrate import trapz
  hdr = dump['hdr']
  N1 = hdr['N1']; N2 = hdr['N2']
  x = flatten_xz(geom['x'], hdr).transpose()
  z = flatten_xz(geom['z'], hdr).transpose()
  A_phi = np.zeros([N2, 2*N1])
  gdet = geom['gdet'].transpose()
  B1 = dump['B1'].mean(axis=-1).transpose()
  B2 = dump['B2'].mean(axis=-1).transpose()
  for j in xrange(N2):
    for i in xrange(N1):
      A_phi[j,N1-1-i] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx1']) - 
                         trapz(gdet[:j, i]*B1[:j, i], dx=hdr['dx2']))
      A_phi[j,i+N1] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx1']) - 
                         trapz(gdet[:j, i]*B1[:j, i], dx=hdr['dx2']))
  A_phi -= (A_phi[N2/2-1,-1] + A_phi[N2/2,-1])/2.
  Apm = np.fabs(A_phi).max()
  if np.fabs(A_phi.min()) > A_phi.max():
    A_phi *= -1.
  levels = np.concatenate((np.linspace(-Apm,0,NLEV)[:-1], 
                           np.linspace(0,Apm,NLEV)[1:]))
  ax.contour(x, z, A_phi, levels=levels, colors='k')

def plot_xy(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
            label=None, ticks=None, arrayspace=False):
  hdr = dump['hdr']; N1 = hdr['N1']; N2 = hdr['N2']; N3 = hdr['N3']+1

  if (arrayspace):
    x = np.reshape(np.repeat(np.linspace(0,1,N3),N1),(N1,N3))
    y = np.transpose(np.reshape(np.repeat(np.linspace(0,1,N1),N3),(N3,N1)))
  else:
    x = geom['x'][:,N2/2,:]
    y = geom['y'][:,N2/2,:]
    x = flatten_xy(x, hdr)
    y = flatten_xy(y, hdr)

  var = flatten_xy(var[:,N2/2,:], hdr)

#  print 'xshape is ', x.shape, ', yshape is ', y.shape, ', varshape is ', var.shape
  mesh = ax.pcolormesh(x, y, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading='gouraud')

  if (not arrayspace):
    circle1=plt.Circle((0,0),hdr['Reh'],color='k');
    ax.add_artist(circle1)

  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    if label:
      plt.colorbar(mesh, cax=cax, ticks=ticks).set_label(label)
    else:
      plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  if (arrayspace):
    ax.set_xlabel('r (arbitrary)'); ax.set_ylabel('phi (arbitrary)')
  else:
    ax.set_xlabel('x/M'); ax.set_ylabel('y/M')
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)
