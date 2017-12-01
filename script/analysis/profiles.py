################################################################################
#                                                                              # 
#  GENERATE INITIAL CONDITIONS EQUATORIAL CUTS                                 # 
#                                                                              # 
################################################################################

import sys; sys.dont_write_bytecode = True
sys.path.insert(0, '../')
import numpy as np
import hdf5_to_dict as io
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import util
import glob
import os

init_file = np.sort(glob.glob(os.path.join(path, "dump*.h5")))[0]

hdr = io.load_hdr(init_file)
geom = io.load_geom(hdr, init_file)
dump = io.load_dump(hdr, geom, files[n])

r = geom['r'][:,geom['N2']/2,0]

fig = plt.figure(figsize=(10, 30))

ax = plt.subplot(2,3,1)
ax.semilogy(r, dump['RHO'][:,geom['N2']/2,0])
ax.set_xlim([0, SIZE]); ax.set_ylim([1e-8, 10])

# ax = plt.subplot(2,3,2)
# bplt.plot_xz(ax, geom, var2_data, dump,
#   label=var2_str, cmap='RdBu_r', vmin=-8, vmax=2)
# ax.set_xlim([0, SIZE])
# 
# ax = plt.subplot(2,3,3)
# bplt.plot_xy(ax, geom, np.log10(dump['RHO']), dump,
#    vmin=-4, vmax=0, label='RHO')
# ax.set_xlim([0, SIZE])
# 
# ax = plt.subplot(2,3,4)
# bplt.plot_xy(ax, geom, var2_data, dump,
#    label=var2_str, cmap='RdBu_r', vmin=-8, vmax=2)
# ax.set_xlim([0, SIZE])

#ax.pcolormesh(dump['X1'][:,:,0], dump['X2'][:,:,0], dump['RHO'][:,:,0])
plt.savefig(imname, bbox_inches='tight', dpi=100)
plt.close(fig)

# Test plot so that backtraces work
# And for quicker turnaround
if debug:
    plot(0)
    plot(100)
    exit(0)

import multiprocessing
import signal
import psutil

nthreads = psutil.cpu_count(logical=False)
print 'Number of CPUs: %i' % psutil.cpu_count(logical=False)

original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
pool = multiprocessing.Pool(nthreads)
signal.signal(signal.SIGINT, original_sigint_handler)
try:
  res = pool.map_async(plot, range(len(files)))
  res.get(720000)
except KeyboardInterrupt:
  print 'Caught interrupt!'
  pool.terminate()
  exit(1)
else:
  pool.close()
pool.join()

