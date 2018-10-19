################################################################################
#                                                                              #
#  UTILITIES FOR PLOTTING                                                      #
#                                                                              #
################################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Get xz slice of 3D data
def flatten_xz(array, patch_pole=False, average=False):
  N1 = array.shape[0]; N2 = array.shape[1]; N3 = array.shape[2]
  flat = np.zeros([2*N1,N2])
  if average:
    for i in range(N1):
      # Produce identical hemispheres to get the right size output
      flat[i,:] = np.mean(array[N1 - 1 - i,:,:], axis=-1)
      flat[i+N1,:] = np.mean(array[i,:,:], axis=-1)
  else:
    for i in range(N1):
      flat[i,:] = array[N1 - 1 - i,:,N3/2]
      flat[i+N1,:] = array[i,:,0]

  # Theta is technically [small,pi/2-small]
  # This patches the X coord so the plot looks nice
  if patch_pole:
      flat[:,0] = 0
      flat[:,-1] = 0

  return flat

# Get xy slice of 3D data
def flatten_xy(array, average=False, loop=True):
  if average:
    slice = np.mean(array, axis=1)
  else:
    slice = array[:,array.shape[1]/2,:]
  
  if loop:
    return np.vstack((slice.transpose(),slice.transpose()[0])).transpose()
  else:
    return slice

def plot_xz(ax, geom, var, cmap='jet', vmin=None, vmax=None, window=None,
            cbar=True, label=None, xlabel=True, ylabel=True,
            ticks=None, arrayspace=False, average=False):

  if (arrayspace):
    x1_norm = (geom['X1'] - geom['startx1']) / (geom['n1'] * geom['dx1'])
    x2_norm = (geom['X2'] - geom['startx2']) / (geom['n2'] * geom['dx2'])
    x = flatten_xz(x1_norm)[geom['n1']:,:]
    z = flatten_xz(x2_norm)[geom['n1']:,:]
    if geom['n3'] > 1:
      var = flatten_xz(var, average=average)[geom['n1']:,:]
    else:
      var = var[:,:,0]
  else:
    x = flatten_xz(geom['x'], patch_pole=True)
    z = flatten_xz(geom['z'])
    var = flatten_xz(var, average=average)

  #print 'xshape is ', x.shape, ', zshape is ', z.shape, ', varshape is ', var.shape
  mesh = ax.pcolormesh(x, z, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading='gouraud')

  if arrayspace:
    if xlabel: ax.set_xlabel('X1 (arbitrary)')
    if ylabel: ax.set_ylabel('X2 (arbitrary)')
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    if xlabel: ax.set_xlabel('x/M')
    if ylabel: ax.set_ylabel('z/M')
    if window:
      ax.set_xlim(window[:2]); ax.set_ylim(window[2:])
    # BH silhouette
    circle1=plt.Circle((0,0),geom['hdr']['Reh'],color='k');
    ax.add_artist(circle1)

  ax.set_aspect('equal')

  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)

  if label:
    ax.set_title(label)

def plot_xy(ax, geom, var, cmap='jet', vmin=None, vmax=None, window=None,
            cbar=True, label=None, xlabel=True, ylabel=True,
            ticks=None, arrayspace=False, average=False):

  if arrayspace:
    # Flatten_xy adds a rank. TODO is this the way to handle it?
    x1_norm = (geom['X1'] - geom['startx1']) / (geom['n1'] * geom['dx1'])
    x3_norm = (geom['X3'] - geom['startx3']) / (geom['n3'] * geom['dx3'])
    x = flatten_xy(x1_norm, loop=False)
    y = flatten_xy(x3_norm, loop=False)
    var = flatten_xy(var, average=average, loop=False)
  else:
    x = flatten_xy(geom['x'])
    y = flatten_xy(geom['y'])
    var = flatten_xy(var, average=average)

  #print 'xshape is ', x.shape, ', yshape is ', y.shape, ', varshape is ', var.shape
  mesh = ax.pcolormesh(x, y, var, cmap=cmap, vmin=vmin, vmax=vmax,
      shading='gouraud')

  if arrayspace:
    if xlabel: ax.set_xlabel('X1 (arbitrary)')
    if ylabel: ax.set_ylabel('X3 (arbitrary)')
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    if xlabel: ax.set_xlabel('x/M')
    if ylabel: ax.set_ylabel('y/M')
    if window:
      ax.set_xlim(window[:2]); ax.set_ylim(window[2:])
    # BH silhouette
    circle1=plt.Circle((0,0),geom['hdr']['Reh'],color='k');
    ax.add_artist(circle1)

  ax.set_aspect('equal')

  if cbar:
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)

  if label:
    ax.set_title(label)

def overlay_field(ax, geom, dump, NLEV):
  from scipy.integrate import trapz
  hdr = dump['hdr']
  N1 = hdr['n1']; N2 = hdr['n2']
  x = flatten_xz(geom['x']).transpose()
  z = flatten_xz(geom['z']).transpose()
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

# TODO allow coordinate x2,3? Allow average over said?
def radial_plot(ax, geom, var, label, n2=0, n3=0, logx=False, logy=False, rlim=None, ylim=None):
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
  if rlim:
    ax.set_xlim(rlim)
  if ylim:
    ax.set_ylim(ylim)
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
