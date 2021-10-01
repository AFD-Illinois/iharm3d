import numpy as np
import matplotlib.pyplot as plt

# prettier text
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif',size=12)

# prettier colors
import seaborn as sns
sns.set_palette("tab10")

# make figure
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

## strong scaling ## 
strong_frontera_256 = np.sort(np.loadtxt("strong_scaling_frontera_256.csv").T)
strong_stampede_256 = np.sort(np.loadtxt("strong_scaling_stampede_256.csv").T)

# get analytic form and adjust first point to be halfway in log space
xs = strong_stampede_256[0]
y0a = strong_frontera_256[1, 0]
y0b = strong_stampede_256[1, 0]
ys = xs/xs[0] * np.exp((np.log(y0a) + np.log(y0b))/2.)
axes[0].plot(xs, ys, c='k', linestyle='--', label="perfect")

# plot actual data
axes[0].plot(strong_frontera_256[0], strong_frontera_256[1], 'o-', label='Frontera')
axes[0].plot(strong_stampede_256[0], strong_stampede_256[1], 'o-', label='Stampede2')

# formatting
axes[0].legend(fontsize=14)
axes[0].set_xlabel(r'\# nodes', fontsize=18)
axes[0].set_ylabel('ZCPS', fontsize=18)
axes[0].set_title(r'strong scaling ($256^3$)', fontsize=20)
axes[0].set_xscale('log')
axes[0].set_yscale('log')

#axes[1].set_yticks([])

## weak scaling ##
weak_frontera_64 = np.sort(np.loadtxt("weak_scaling_frontera_64.csv").T)
weak_stampede_64 = np.sort(np.loadtxt("weak_scaling_stampede_64.csv").T)

# get analytic form and adjust first point to be halfway in log space
xs = weak_stampede_64[0]
y0a = weak_stampede_64[1, 0]# TODO
y0b = weak_stampede_64[1, 0]
ys = xs/xs[0] * np.exp((np.log(y0a) + np.log(y0b))/2.)
axes[1].plot(xs, ys, c='k', linestyle='--', label="perfect")

# plot actual data
axes[1].plot(weak_frontera_64[0], weak_frontera_64[1], 'o-', label='Frontera')
axes[1].plot(weak_stampede_64[0], weak_stampede_64[1], 'o-', label='Stampede2')

axes[1].set_title(r'weak scaling ($64^3$ per node)', fontsize=20)
axes[1].set_xscale('log')
axes[1].set_yscale('log')

# save    
plt.tight_layout(rect=[0, 0, 1, 1])
plt.savefig('scaling.pdf')
