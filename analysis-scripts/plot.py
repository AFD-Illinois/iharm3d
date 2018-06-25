################################################################################
#                                                                              #
#  UTILITIES FOR PLOTTING                                                      #
#                                                                              #
################################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Get xz slice of 3D data
def flatten_xz(array, hdr, patch_pole=False):
  N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3']
  flat = np.zeros([2*N1,N2])
  for j in xrange(N2):
      for i in xrange(N1):
          flat[i,j] = array[N1 - 1 - i,j,N3/2]
      for i in xrange(N1):
          flat[i+N1,j] = array[i,j,0]
  # Theta is technically [small,pi/2-small]
  # This patches the X coord so the plot looks nice
  if patch_pole:
      flat[:,0] = 0
      flat[:,-1] = 0
  return flat

# Get xy slice of 3D data
def flatten_xy(array, hdr):
  return np.vstack((array.transpose(),array.transpose()[0])).transpose()

def plot_xz(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
            label=None, ticks=None, arrayspace=False):
  hdr = dump['hdr']; N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3']
  if N3 > 1.:
    if (arrayspace):
      x = np.reshape(np.repeat(np.linspace(0,1,N2),N1),(N1,N2))
      z = np.transpose(np.reshape(np.repeat(np.linspace(0,1,N1),N2),(N2,N1)))
      var = flatten_xz(var, hdr)[N2:,:]
    else:
      x = flatten_xz(geom['x'], hdr, patch_pole=True)
      z = flatten_xz(geom['z'], hdr)
      var = flatten_xz(var, hdr)

  else:
    x = x[:,:,0]
    z = z[:,:,0]
    var = var[:,:,0]

  #print 'xshape is ', x.shape, ', zshape is ', z.shape, ', varshape is ', var.shape
  mesh = ax.pcolormesh(x, z, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading='gouraud')

  if not arrayspace:
    circle1=plt.Circle((0,0),hdr['Reh'],color='k');
    ax.add_artist(circle1)

  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)
  ax.set_aspect('equal')
  if (arrayspace):
    ax.set_xlabel('r (arbitrary)'); ax.set_ylabel('theta (arbitrary)')
  else:
    ax.set_xlabel('x/M'); ax.set_ylabel('z/M')
  
  if label:
    ax.set_title(label)
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)

def overlay_field(ax, geom, dump, NLEV):
  from scipy.integrate import trapz
  hdr = dump['hdr']
  N1 = hdr['n1']; N2 = hdr['n2']
  x = flatten_xz(geom['x'], hdr).transpose()
  z = flatten_xz(geom['z'], hdr).transpose()
  A_phi = np.zeros([N2, 2*N1])
  gdet = geom['gdet'][:,:].transpose()
  B1 = dump['B1'].mean(axis=-1).transpose()
  B2 = dump['B2'].mean(axis=-1).transpose()
  for j in xrange(N2):
    for i in xrange(N1):
      A_phi[j,N1-1-i] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx1']) -
                         trapz(gdet[:j,i]*B1[:j,i], dx=hdr['dx2']))
      A_phi[j,i+N1] = (trapz(gdet[j,:i]*B2[j,:i], dx=hdr['dx1']) -
                         trapz(gdet[:j,i]*B1[:j,i], dx=hdr['dx2']))
  A_phi -= (A_phi[N2/2-1,-1] + A_phi[N2/2,-1])/2.
  Apm = np.fabs(A_phi).max()
  if np.fabs(A_phi.min()) > A_phi.max():
    A_phi *= -1.
  levels = np.concatenate((np.linspace(-Apm,0,NLEV)[:-1],
                           np.linspace(0,Apm,NLEV)[1:]))
  ax.contour(x, z, A_phi, levels=levels, colors='k')

def plot_xy(ax, geom, var, dump, cmap='jet', vmin=None, vmax=None, cbar=True,
            label=None, ticks=None, arrayspace=False):
  hdr = dump['hdr']; N1 = hdr['n1']; N2 = hdr['n2']; N3 = hdr['n3']+1

  if (arrayspace):
    x = np.reshape(np.repeat(np.linspace(0,1,N3),N1),(N1,N3))
    y = np.transpose(np.reshape(np.repeat(np.linspace(0,1,N1),N3),(N3,N1)))
  else:
    x = geom['x'][:N1,N2/2,:N3-1]
    y = geom['y'][:N1,N2/2,:N3-1]
    x = flatten_xy(x, hdr)
    y = flatten_xy(y, hdr)

  var = flatten_xy(var[:,N2/2,:], hdr)

  #print 'xshape is ', x.shape, ', yshape is ', y.shape, ', varshape is ', var.shape
  mesh = ax.pcolormesh(x, y, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading='gouraud')

  if (not arrayspace):
    circle1=plt.Circle((0,0),hdr['Reh'],color='k');
    ax.add_artist(circle1)

  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)
    if label:
      ax.set_title(label)
  ax.set_aspect('equal')
  if (arrayspace):
    ax.set_xlabel('r (arbitrary)'); ax.set_ylabel('phi (arbitrary)')
  else:
    ax.set_xlabel('x/M'); ax.set_ylabel('y/M')
  #ax.grid(True, linestyle=':', color='k', alpha=0.5, linewidth=0.5)

# TODO allow coordinate x2,3? Allow average over said?
def plot_r(ax, geom, var, n2, n3, label, logx=False, logy=False):
  r = geom['r'][:,0,0]
  if var.ndim == 1:
    data = var
  elif var.ndim == 2:
    data = var[:,n2]
  elif var.ndim == 3:
    data = var[:,n2,n3]
  
  if logx and logy:
    ax.loglog(r,data)
  elif logx:
    ax.semilogx(r,data)
  elif logy:
    ax.semilogy(r,data)
  else:
    ax.plot(r,data)
  ax.set_xlabel('r (M)')
  ax.set_ylabel(label)

def diag_plot(ax, diag, dump, varname_dump, varname_pretty, ylim=None, logy=False):
  var = diag[varname_dump]
  slc = np.nonzero(var)
  if logy:
    ax.semilogy(diag['t'][slc], var[slc], color='k')
  else:
    ax.plot(diag['t'][slc], var[slc], color='k')
  ax.axvline(dump['t'], color='r') # Trace current t on finished plot
  ax.set_xlim([0, dump['hdr']['tf']])
  if ylim is not None:
    ax.set_ylim(ylim)
  ax.set_xlabel('t/M')
  ax.set_ylabel(varname_pretty)
