import numpy as np
import os, h5py, glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Function to overlay field lines
# Argument must be axes object, B1, B2 and 'nlines' -> a parameter to account for density of field lines
def plotting_bfield_lines(ax, B1, B2, nlines=20):
  B1_phi_avg = B1.mean(axis=-1)
  B2_phi_avg = B2.mean(axis=-1)
  AJ_phi = np.zeros([grid['n1'],grid['n2']])
  for j in range(grid['n2']):
    for i in range(grid['n1']):
      AJ_phi[i,j] = (np.trapz(B2_phi_avg[:i,j], dx=grid['dx1']) - np.trapz(B1_phi_avg[i,:j], dx=grid['dx2']))
  AJ_phi -= AJ_phi.min()
  levels = np.linspace(0, AJ_phi.max(), nlines*2)
  # ax.contour(xp, zp, AJ_phi, levels=levels, colors='k')
  ax.contour(grid['x'][Ellipsis,0], grid['y'][Ellipsis,0], AJ_phi, levels=levels, colors='k')

grid = {}
def load_grid(dumpsdir):
	dfile = h5py.File(sorted(glob.glob(os.path.join('dumps', 'dump_000*.h5')))[0],'r')
	grid['dx1'] = dfile['/header/geom/dx1'][()]
	grid['dx2'] = dfile['/header/geom/dx2'][()]
	grid['dx3'] = dfile['/header/geom/dx3'][()]
	grid['ndim'] = dfile['/header/geom/n_dim'][()]
	grid['n1'] = dfile['header/n1'][()]
	grid['n2'] = dfile['header/n2'][()]
	grid['n3'] = dfile['header/n3'][()]
	dfile.close()
	gfile = h5py.File(os.path.join(dumpsdir, 'grid.h5'),'r')
	grid['X1'] = gfile['X1'][()]
	grid['X2'] = gfile['X2'][()]
	grid['X3'] = gfile['X3'][()]
	grid['r'] = gfile['r'][()]
	grid['th'] = gfile['th'][()]
	grid['phi'] = gfile['phi'][()]
	grid['x'] = gfile['X'][()]
	grid['y'] = gfile['Y'][()]
	grid['z'] = gfile['Z'][()]
	grid['lapse'] = gfile['lapse'][()]
	grid['gcov'] = gfile['gcov'][()]
	grid['gcon'] = gfile['gcon'][()]
	grid['gdet'] = gfile['gdet'][()]
	gfile.close()

def plot_temperature(dumpsdir, plotsdir, dumpval):
  plt.clf()
  print("Analyzing {0:04d} dump".format(dumpval))
  dfile = h5py.File(os.path.join(dumpsdir, 'dump_0000{0:04d}.h5'.format(dumpval)),'r')
  rho = dfile['prims'][()][Ellipsis,0]
  u   = dfile['prims'][()][Ellipsis,1]
  B   = np.array(dfile['prims'][()][Ellipsis,6:])
  Theta = 1./3. * u / rho
  t = dfile['t'][()]
  dfile.close()
  t = "{:.3f}".format(t)

  fig = plt.figure(figsize=(16,9))
  heights = [1,5]
  gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=heights, figure=fig)
  
  cmap = 'jet'
  vmin = 0.33
  vmax = 0.35
  domain  = [0,1,0,1]
  shading ='gouraud'

  ax0 = fig.add_subplot(gs[0,:])
  ax0.annotate('t= '+str(t),xy=(0.5,0.5),xycoords='axes fraction',va='center',ha='center',fontsize='x-large')
  ax0.axis("off")

  ax1 = fig.add_subplot(gs[1,:])
  Theta_xyplot = ax1.pcolormesh(grid['x'][:,:,0], grid['y'][:,:,0], Theta[:,:,0], vmin=vmin, vmax=vmax, cmap=cmap, shading=shading)
  plotting_bfield_lines(ax1, B[Ellipsis,0], B[Ellipsis,1], nlines=40)
  ax1.set_xlabel('$x$' ,fontsize='x-large')
  ax1.set_ylabel('$y$', fontsize='x-large')
  ax1.set_xlim(domain[:2])
  ax1.set_ylim(domain[2:])
  ax1.set_title('$\\Theta$',fontsize='xx-large')
  ax1.set_aspect('equal')
  divider = make_axes_locatable(ax1)
  cax = divider.append_axes("right", size="5%", pad=0.05)
  plt.colorbar(Theta_xyplot, cax=cax)

  plt.tight_layout()
  plt.savefig(os.path.join(plotsdir,'snake_test_{:04d}.png'.format(dumpval)))
  plt.close()

if __name__=='__main__':
  dumpsdir = './dumps'
  plotsdir = './plots'
  Ndumps   = len(glob.glob(os.path.join('dumps', 'dump_000000*.h5')))

  load_grid(dumpsdir)

  for dumpval in range(Ndumps):
    plot_temperature(dumpsdir, plotsdir, dumpval)
