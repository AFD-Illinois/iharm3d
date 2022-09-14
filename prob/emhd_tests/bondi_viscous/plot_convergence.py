import numpy as np
import os, glob, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

RES   = [32,64,128,256]
PRIMS = ['rho','u','dP']
L1_norm = np.zeros([len(RES), len(PRIMS)])

for r, res in enumerate(RES):

  # load analytic result
  rho_analytic, u_analytic, dP_analytic  = np.loadtxt(os.path.join('{}'.format(res), 'bondi_analytic_{}.txt'.format(res)), usecols=(0,1,3), unpack=True)

  # load grid file
  dumpsdir = os.path.join('{}'.format(res), 'dumps')
  gfile    = h5py.File(os.path.join(dumpsdir, 'grid.h5'), 'r')
  rcoord   = np.squeeze(gfile['r'][()])[:,0]
  gfile.close()
  1
  # load code data
  dfile = h5py.File(sorted(glob.glob(os.path.join(dumpsdir, 'dump_*.h5')))[-1], 'r')

  rho      = np.squeeze(dfile['prims'][Ellipsis,0][()])
  u        = np.squeeze(dfile['prims'][Ellipsis,1][()])
  dP_tilde = np.squeeze(dfile['prims'][Ellipsis,5][()])

  t   = dfile['t'][()]
  gam = dfile['header/gam'][()]
  rEH = dfile['header/geom/mks/r_eh'][()]
  higher_order_terms = dfile['header/imex/emhd/higher_order_terms_viscosity'][()]

  dfile.close()

  # consider only those zones outside EH
  rEH_ind = np.argmin(np.fabs(rcoord - rEH) > 0.)
  rho = rho[rEH_ind:,:];           rho_analytic = rho_analytic[rEH_ind:]
  u   = u[rEH_ind:,:];             u_analytic = u_analytic[rEH_ind:]
  dP_tilde = dP_tilde[rEH_ind:,:]; dP_analytic = dP_analytic[rEH_ind:]

  # compute dP
  tau      = 30.
  eta      = 0.01
  P        = (gam - 1.) * u
  Theta    = P / rho
  nu_emhd  = eta / rho

  if higher_order_terms:
    dP = dP_tilde * np.sqrt(nu_emhd * rho * Theta / tau)
  else:
    dP = dP_tilde

  # compute L1 norm
  L1_norm[r,0] = np.mean(np.fabs(rho-rho_analytic[:,None]))
  L1_norm[r,1] = np.mean(np.fabs(u-u_analytic[:,None]))
  L1_norm[r,2] = np.mean(np.fabs(dP[:-1,:]-dP_analytic[:-1,None]))


# plotting parameters
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
mpl.rcParams['figure.autolayout'] = True
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['axes.xmargin'] = 0.02
mpl.rcParams['axes.ymargin'] = 0.02
mpl.rcParams['legend.fontsize'] = 'medium'
colors = ['indigo', 'goldenrod', 'darkgreen', 'crimson', 'xkcd:blue']

# plot
plt.close()
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1,1,1)

# loop over prims
tracker = 0
for n in range(len(PRIMS)):
  color = colors[tracker]
  ax.loglog(RES, L1_norm[:,n], color=color, marker='o', label=PRIMS[n])
  tracker+=1

ax.loglog([RES[0], RES[-1]], 0.1*np.asarray([float(RES[0]), float(RES[-1])])**(-2), color='k', linestyle='dashed', label='$N^{-2}$')
plt.xscale('log', base=2)
ax.set_xlabel('Resolution')
ax.set_ylabel('L1 norm')
ax.legend()
plt.savefig('bondi_viscous_convergence.png', dpi=300)
