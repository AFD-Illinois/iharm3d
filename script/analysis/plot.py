################################################################################
#                                                                              #
#  UTILITIES FOR PLOTTING                                                      #
#                                                                              #
################################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from mpl_toolkits.axes_grid1 import make_axes_locatable

hdr = {}
geom = {}

def init_plotting(hdr_init, geom_init):
  global hdr, geom
  hdr = hdr_init
  geom = geom_init

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

def plot_xz(ax, var, cmap='jet', vmin=None, vmax=None, window=[-40,40,-40,40],
            cbar=True, label=None, xlabel=True, ylabel=True,
            ticks=None, arrayspace=False, average=False, bh=True):

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
    if xlabel: ax.set_xlabel("X1 (arbitrary)")
    if ylabel: ax.set_ylabel("X2 (arbitrary)")
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    if xlabel: ax.set_xlabel(r"$x \frac{c^2}{G M}$")
    if ylabel: ax.set_ylabel(r"$z \frac{c^2}{G M}$")
    if window:
      ax.set_xlim(window[:2]); ax.set_ylim(window[2:])

    if bh:
      # BH silhouette
      circle1=plt.Circle((0,0), hdr['Reh'], color='k');
      ax.add_artist(circle1)

  ax.set_aspect('equal')

  if cbar:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)

  if label:
    ax.set_title(label)

def plot_xy(ax, var, cmap='jet', vmin=None, vmax=None, window=[-40,40,-40,40],
            cbar=True, label=None, xlabel=True, ylabel=True,
            ticks=None, arrayspace=False, average=False, bh=True):

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
    if xlabel: ax.set_xlabel("X1 (arbitrary)")
    if ylabel: ax.set_ylabel("X3 (arbitrary)")
    ax.set_xlim([0, 1]); ax.set_ylim([0, 1])
  else:
    if xlabel: ax.set_xlabel(r"$x \frac{c^2}{G M}$")
    if ylabel: ax.set_ylabel(r"$y \frac{c^2}{G M}$")
    if window:
      ax.set_xlim(window[:2]); ax.set_ylim(window[2:])

    if bh:
      # BH silhouette
      circle1=plt.Circle((0,0), hdr['Reh'], color='k');
      ax.add_artist(circle1)

  ax.set_aspect('equal')

  if cbar:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)

  if label:
    ax.set_title(label)

def overlay_contour(ax, var, levels):
  x = flatten_xz(geom['x'])
  z = flatten_xz(geom['z'])
  var = flatten_xz(var, average=True)
  ax.contour(x, z, var, levels=levels, colors='k')

def overlay_field(ax, dump, NLEV):
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

# Usual slice plots -- plot an _xy and _xz slice to two given axes
def plot_slices(ax1, ax2, name, data, dump, min, max, avg=False, int=False, window=[-40,40,-40,40], arrspace=False,
                field_overlay=True, nlines=10, cmap='jet', xlabel=True, ylabel=True, vert_split=False):
  if int:
    # Average multiplied values for the integral
    data_xz = data*data.shape[2] #N1,2,(3)
    data_xy = data*data.shape[1] #N1,(2),3
    avg = True
  else:
    data_xz = data
    data_xy = data
  
  plot_xz(ax1, data_xz, window=window, cbar=False, cmap=cmap, xlabel=xlabel, ylabel=ylabel,
               label=name, vmin=min, vmax=max, arrayspace=arrspace, average=avg)
  if field_overlay and not arrspace:
    overlay_field(ax1, dump, nlines)

  plot_xy(ax2, data_xy, window=window, cmap=cmap, xlabel=xlabel, ylabel=ylabel,
               label=name, vmin=min, vmax=max, arrayspace=arrspace, average=avg)

# TODO allow coordinate x2,3? Allow average over said?
def radial_plot(ax, var, label, n2=0, n3=0, average=False, logx=False, logy=False, rlim=None, ylim=None, arrayspace=False, col='k'):
  r = geom['r'][:,0,0]
  if var.ndim == 1:
    data = var
  elif var.ndim == 2:
    data = var[:,n2]
  elif var.ndim == 3:
    if average:
      data = np.mean(var[:,n2,:], axis=-1)
    else:
      data = var[:,n2,n3]

  # TODO there's probably a way to add log
  if arrayspace:
    # logx doesn't make sense here
    if logy:
      ax.semilogy(range(geom['n1']), data, color=col)
    else:
      ax.plot(range(geom['n1']), data, color=col)
  else:
    if logx and logy:
      ax.loglog(r,data, color=col)
    elif logx:
      ax.semilogx(r,data, color=col)
    elif logy:
      ax.semilogy(r,data, color=col)
    else:
      ax.plot(r,data, color=col)

  if rlim:
    ax.set_xlim(rlim)
  if ylim:
    ax.set_ylim(ylim)
  ax.set_xlabel(r"$r \frac{c^2}{G M}$")
  ax.set_ylabel(label)

def diag_plot(ax, diag, dump, varname_dump, varname_pretty, ylim=None, logy=False, xlabel=True):
  var = diag[varname_dump]
  slc = np.nonzero(var)
  if logy:
    ax.semilogy(diag['t'][slc], var[slc], color='k')
  else:
    ax.plot(diag['t'][slc], var[slc], color='k')
  ax.axvline(dump['t'], color='r') # Trace current t on finished plot
  ax.set_xlim([0, hdr['tf']])
  if ylim is not None:
    ax.set_ylim(ylim)
  if xlabel:
    ax.set_xlabel(r"$t \frac{c^3}{G M}$")
  ax.set_ylabel(varname_pretty)
  
def hist_2d(ax, var_x, var_y, xlbl, ylbl, title=None, logcolor=False, bins=40, cbar=True, ticks=None):
  # Courtesy of George Wong
  var_x_flat = var_x.flatten()
  var_y_flat = var_y.flatten()
  nidx = np.isfinite(var_x_flat) & np.isfinite(var_y_flat)
  hist = np.histogram2d(var_x_flat[nidx], var_y_flat[nidx], bins=bins)
  X,Y = np.meshgrid(hist[1], hist[2])
  
  if logcolor:
    hist[0][hist[0] == 0] = np.min(hist[0][np.nonzero(hist[0])])
    mesh = ax.pcolormesh(X,Y,np.log10(hist[0]))
  else:
    mesh = ax.pcolormesh(X,Y,hist[0])

  # Add the patented Ben Ryan colorbar
  if cbar:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mesh, cax=cax, ticks=ticks)

  if title is not None: ax.set_title(title)
  ax.set_xlabel(xlbl)
  ax.set_ylabel(ylbl)
