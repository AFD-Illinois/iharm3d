import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# prettier text
plt.rc('text', usetex=True)

plt.rc('font', family='sans-serif',size=12)

# deal with unified color cycle
sns.set_palette("tab10")
color_cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']


fig, ax = plt.subplots(1,4, figsize=(10.5,3))

vars = ("N", r"$\rho$", "$u$", r"$u^1$", r"$u^2$", r"$u^3$", r"$B^1$", r"$B^2$", r"$B^3$")
for i,(name,pretty_name) in enumerate(zip(['SLOW', 'FAST', 'ALFVEN', 'bondi'], ["Slow wave", "Fast wave", "Alfvén wave", "Bondi flow"])):
    cols = (list(range(9)), list(range(9)), (0,4,5,7,8), (0,1))[i]
    data = np.loadtxt('data_'+name+'.csv', skiprows=1, usecols=cols, delimiter=',').T

    # plot analytic behavior
    minn_maxn = np.array([data[0,0]/2, data[0,-1]*2])
    ax[i].plot(minn_maxn,
        4*np.mean(data[1:,0])*minn_maxn[0]**2/minn_maxn**2,
        color='k', linestyle='--', label=r"$N^{-2}$")

    # plot convergence data
    for j,col in enumerate(cols[1:]):
        if pretty_name in ["Slow wave", "Fast wave"] and col==5:
            # overplotting u^2 and u^3
            marker_style = dict(marker='o', fillstyle='none')
            ax[i].loglog(data[0], data[j+1], 'o--', c=color_cycle[col-1], label=vars[col], ms=4, **marker_style)
        else:
            ax[i].loglog(data[0], data[j+1], 'o-', c=color_cycle[col-1], label=vars[col], ms=4)
       
    # formatting     
    ax[i].set_xscale('log', base=2)
    ax[i].set_ylim((1e-9,2e-5))
    ax[i].set_title(pretty_name)
    if i > 0:
        ax[i].set_yticklabels([])
        ax[i].set_yticks([])
        
handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='right', fontsize=14)
    
fig.text(0.02, 0.5, 'L1 norm', fontsize=20, rotation=90, horizontalalignment='center', verticalalignment='center')
fig.text(0.48, 0.01, '$N$', fontsize=20, horizontalalignment='center', verticalalignment='bottom')
  
plt.tight_layout(rect=[0.025, 0.03, 0.9, 1])
plt.savefig("convergence.pdf")
