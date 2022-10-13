import numpy as np
import os, glob, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt


RES   = [32,64,128,256]
NG    = 3
CONDUCTION = 1
if CONDUCTION:
  PRIMS = ['rho','u','q']
else:
  PRIMS = ['rho','u']
L1_norm = np.zeros([len(RES), len(PRIMS)])

for r, res in enumerate(RES):

  # load analytic result
  rho_analytic = np.loadtxt(os.path.join('{}'.format(res), 'atmosphere_soln_rho.txt'))[NG:-NG]
  u_analytic   = np.loadtxt(os.path.join('{}'.format(res), 'atmosphere_soln_u.txt'))[NG:-NG]
  if CONDUCTION:
    q_analytic = np.loadtxt(os.path.join('{}'.format(res), 'atmosphere_soln_phi.txt'))[NG:-NG]

  # load code data
  dumpsdir = os.path.join('{}'.format(res), 'dumps')
  dfile = h5py.File(sorted(glob.glob(os.path.join(dumpsdir, 'dump_*.h5')))[-1], 'r')

  rho = np.squeeze(dfile['prims'][Ellipsis,0][()])
  u   = np.squeeze(dfile['prims'][Ellipsis,1][()])
  if CONDUCTION:
    q_tilde   = np.squeeze(dfile['prims'][Ellipsis,5][()])

  t   = dfile['t'][()]
  gam = dfile['header/gam'][()]

  # compute q
  tau      = 10.
  kappa    = 0.1
  P        = (gam - 1.) * u
  Theta    = P / rho
  cs2      = (gam * P) / (rho + (gam * u))
  chi_emhd = kappa / rho
  q        = q_tilde * np.sqrt(chi_emhd * rho * Theta**2 / tau)

  # compute L1 norm
  L1_norm[r,0] = np.mean(np.fabs(rho-rho_analytic[:,None]))
  L1_norm[r,1] = np.mean(np.fabs(u-u_analytic[:,None]))
  if CONDUCTION:
    L1_norm[r,2] = np.mean(np.fabs(q-q_analytic[:,None])[1:-1])


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
ax.loglog([RES[0], RES[-1]], 1.e-4*np.asarray([float(RES[0]), float(RES[-1])])**(-2), color='k', linestyle='dashed', label='$N^{-2}$')
#ax.loglog([RES[0], RES[-1]], 0.001*np.asarray([float(RES[0]), float(RES[-1])])**(-1), color='k', linestyle='dotted', label='$N^{-1}$')
plt.xscale('log', base=2)
ax.set_xlabel('Resolution')
ax.set_ylabel('L1 norm')
ax.legend()
plt.savefig('conducting_atmosphere_convergence.png', dpi=300)