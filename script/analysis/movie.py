################################################################################
#                                                                              #
#  GENERATE MOVIES FROM SIMULATION OUTPUT                                      #
#                                                                              #
################################################################################

import hdf5_to_dict as io
import plot as bplt
from analysis_fns import *
import util
from luminosity_th_study import overlay_thphi_contours

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os, sys
import pickle
import numpy as np

# Movie size in inches. Keep 16/9 for standard-size movies
FIGX = 24
FIGY = FIGX*9/16

# For plotting debug, "array-space" plots
# Certain plots can override this below
USEARRSPACE = False

LOG_MDOT = False
LOG_PHI = False

# Load diagnostic data from post-processing (eht_out.p)
diag_post = True

def plot(n):
  imname = os.path.join(frame_dir, 'frame_%08d.png' % n)
  tdump = io.get_dump_time(files[n])
  if (tstart is not None and tdump < tstart) or (tend is not None and tdump > tend):
    return
  
  print("{} / {}".format((n+1),len(files)))

  fig = plt.figure(figsize=(FIGX, FIGY))

  if movie_type not in ["simplest", "simpler", "simple"]:
    dump = io.load_dump(files[n], hdr, geom, derived_vars=True, extras=False)
    #fig.suptitle("t = %d"%dump['t']) # TODO put this at the bottom somehow?
  else:
    # Simple movies don't need derived vars
    dump = io.load_dump(files[n], hdr, geom, derived_vars=False, extras=False)

  # Zoom in for small problems
  if geom['r'][-1,0,0] > 100:
    window = [-40,40,-40,40]
    nlines = 20
    rho_l, rho_h = -3, 2
  else:
    window = [-20,20,-20,20]
    nlines = 5
    rho_l, rho_h = -3, 1

  if movie_type == "simplest":
    # Simplest movie: just RHO
    ax_slc = plt.subplots(1,2)
    bplt.plot_slices(ax_slc[0], ax_slc[1], geom, dump, np.log10(dump['RHO']),
                     label=r"$\log_{10}(\rho)$", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
  elif movie_type == "simpler":
    # Simpler movie: RHO and phi
    gs = gridspec.GridSpec(2, 2, height_ratios=[6, 1], width_ratios=[16,17])
    ax_slc = [fig.subplot(gs[0,0]), fig.subplot(gs[0,1])]
    ax_flux = [fig.subplot(gs[1,:])]
    bplt.plot_slices(ax_slc[0], ax_slc[1], geom, dump, np.log10(dump['RHO']),
                     label=r"$\log_{10}(\rho)$", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
    bplt.diag_plot(ax_flux[0], diag, 'phi_b', dump['t'], ylabel=r"$\phi_{BH}$", logy=LOG_PHI, xlabel=False)
  elif movie_type == "simple":
    # Simple movie: RHO mdot phi
    gs = gridspec.GridSpec(3, 2, height_ratios=[4, 1, 1])
    ax_slc = [fig.subplot(gs[0,0]), fig.subplot(gs[0,1])]
    ax_flux = [fig.subplot(gs[1,:]), fig.subplot(gs[2,:])]
    bplt.plot_slices(ax_slc[0], ax_slc[1], geom, dump, np.log10(dump['RHO']),
                     label=r"$\log_{10}(\rho)$", vmin=rho_l, vmax=rho_h, window=window, cmap='jet')
    bplt.diag_plot(ax_flux[0], diag, 'mdot', dump['t'], ylabel=r"$\dot{M}$", logy=LOG_MDOT)
    bplt.diag_plot(ax_flux[1], diag, 'phi_b', dump['t'], ylabel=r"$\phi_{BH}$", logy=LOG_PHI)
  elif movie_type == "radial":

    # TODO just record these in analysis output...
    rho_r = eht_profile(geom, dump['RHO'], jmin, jmax)
    B_r = eht_profile(geom, np.sqrt(dump['bsq']), jmin, jmax)
    uphi_r = eht_profile(geom, dump['ucon'][:,:,:,3], jmin, jmax)
     
    Pg = (hdr['gam']-1.)*dump['UU']
    Pb = dump['bsq']/2
     
    Pg_r = eht_profile(geom, Pg, jmin, jmax)
    Ptot_r = eht_profile(geom, Pg + Pb, jmin, jmax)
    betainv_r = eht_profile(geom, Pb/Pg, jmin, jmax)

    ax_slc = lambda i: plt.subplot(2, 3, i)
    bplt.radial_plot(ax_slc(1), geom, rho_r, ylabel=r"$<\rho>$", logy=True, ylim=[1.e-2, 1.e0])
    bplt.radial_plot(ax_slc(2), geom, Pg_r, ylabel=r"$<P_g>$", logy=True, ylim=[1.e-6, 1.e-2])
    bplt.radial_plot(ax_slc(3), geom, B_r, ylabel=r"$<|B|>$", logy=True, ylim=[1.e-4, 1.e-1])
    bplt.radial_plot(ax_slc(4), geom, uphi_r, ylabel=r"$<u^{\phi}>$", logy=True, ylim=[1.e-3, 1.e1])
    bplt.radial_plot(ax_slc(5), geom, Ptot_r, ylabel=r"$<P_{tot}>$", logy=True, ylim=[1.e-6, 1.e-2])
    bplt.radial_plot(ax_slc(6), geom, betainv_r, ylabel=r"$<\beta^{-1}>$", logy=True, ylim=[1.e-2, 1.e1])
  
  elif movie_type == "fluxes_cap":
    
    if hdr['r_out'] < 100:
      iBZ = i_of(geom,40) # most SANEs
      rBZ = 40
      rstring = "40"
    else:
      iBZ = i_of(geom,100) # most MADs
      rBZ = 100
      rstring = "100"
    
    axes = [plt.subplot(2, 2, i) for i in range(1,5)]
    bplt.plot_thphi(axes[0], geom, np.log10(d_fns['FE'](dump)[iBZ,:,:]), iBZ, vmin=-8, vmax=-4, label =r"FE $\theta-\phi$ slice")
    bplt.plot_thphi(axes[1], geom, np.log10(d_fns['FM'](dump)[iBZ,:,:]), iBZ, vmin=-8, vmax=-4, label =r"FM $\theta-\phi$ slice")
    bplt.plot_thphi(axes[2], geom, np.log10(d_fns['FL'](dump)[iBZ,:,:]), iBZ, vmin=-8, vmax=-4, label =r"FL $\theta-\phi$ slice")
    bplt.plot_thphi(axes[3], geom, np.log10(dump['RHO'][iBZ,:,:]), iBZ, vmin=-4, vmax=1, label =r"\rho $\theta-\phi$ slice")
    
    for i,axis in enumerate(axes):
      if i == 0:
        overlay_thphi_contours(axis, geom, diag, legend=True)
      else:
        overlay_thphi_contours(axis, geom, diag)
      max_th = geom['n2']//2
      x = bplt.loop_phi(geom['x'][iBZ,:max_th,:])
      y = bplt.loop_phi(geom['y'][iBZ,:max_th,:])
      prep = lambda var : bplt.loop_phi(var[:max_th,:])
      
      axis.contour(x,y, prep(geom['th'][iBZ]), [1.0], colors='k')
      axis.contour(x,y, prep(d_fns['betagamma'](dump)[iBZ]), [1.0], colors='k')
      axis.contour(x,y, prep(d_fns['sigma'](dump)[iBZ]), [1.0], colors='xkcd:green')
      axis.contour(x,y, prep(d_fns['FE'](dump)[iBZ]), [0.0], colors='xkcd:pink')
      axis.contour(x,y, prep(d_fns['Be_nob'](dump)[iBZ]), [0.02], colors='xkcd:red')
      axis.contour(x,y, prep(d_fns['mu'](dump)[iBZ]), [2.0], colors='xkcd:blue')
  
  elif movie_type == "rho_cap":
    if hdr['r_out'] < 100:
      iBZ = i_of(geom,40) # most SANEs
      rBZ = 40
      rstring = "40"
    else:
      iBZ = i_of(geom,100) # most MADs
      rBZ = 100
      rstring = "100"

    bplt.plot_slices(plt.subplot(1,3,1), plt.subplot(1,3,2), geom, dump, np.log10(dump['RHO']),
                   label=r"$\log_{10}(\rho)$", vmin=-3, vmax=2, cmap='jet')
    bplt.overlay_contours(plt.subplot(1,3,1), geom, geom['r'], [rBZ], color='k')
    
    bplt.plot_thphi(plt.subplot(1,3,3), geom, np.log10(dump['RHO'][iBZ,:,:]), iBZ, vmin=-3, vmax=2, label=r"$\log_{10}(\rho)$ $\theta-\phi$ slice")
  
  elif movie_type == "lum_cuts":
      
      gs = gridspec.GridSpec(2, 2, width_ratios=[1,2])
      
      ax = plt.subplot(gs[0,0])
      bplt.plot_xz(ax, geom, np.log10(d_fns['FE_EM'](dump)), average=True, window=window)
      ax.set_title(r"$\log_{10}( -{{T_{EM}}^r}_t )$")
      
      bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')
      
      bplt.overlay_contours(ax, geom, d_fns['sigma'](dump), [1.0], color='C2')
      #bplt.overlay_contours(ax, geom, d_fns['Be_nob'](dump), [0.02], color='C3')
      bplt.overlay_contours(ax, geom, d_fns['betagamma'](dump), [1.0], color='C5')
      
      ax = plt.subplot(gs[1,0])
      bplt.plot_xz(ax, geom, np.log10(d_fns['FE']), average=True, window=window)
      ax.set_title(r"$\log_{10}( -{T^r}_t - \rho u^r )$")
      
      bplt.overlay_contours(ax, geom, geom['r'], [AT_R], color='k')
      
      bplt.overlay_contours(ax, geom, d_fns['sigma'](dump), [1.0], color='C2')
      #bplt.overlay_contours(ax, geom, d_fns['Be_nob'](dump), [0.02], color='C3')
      bplt.overlay_contours(ax, geom, d_fns['betagamma'](dump), [1.0], color='C5')
      
      ax = plt.subplot(gs[0,1])
      ax.plot(avg['r'], avg['LBZ_sigma1_rt'][n], label=r"$L_{BZ}$ (sigma > 1 cut)", color='C2')
      ax.plot(avg['r'], avg['LBZ_Be_nob0_rt'][n], label=r"$L_{BZ}$ ($Be > 0.02$ cut)", color='C3')
      ax.plot(avg['r'], avg['LBZ_Be_nob1_rt'][n], label=r"$L_{BZ}$ ($Be > 1.0$ cut)", color='C4')
      ax.plot(avg['r'], avg['LBZ_bg1_rt'][n], label=r"$L_{BZ}$ ($\beta\gamma > 1.0$ cut)", color='C5')
      ax.plot(avg['r'], avg['LBZ_allp_rt'][n], label=r"$L_{BZ,tot}$", color='C6')
      
      ax.set_title(r"$L_{BZ} = \int -{{T_{EM}}^r}_t \sqrt{-g} dx^{\theta} dx^{\phi}$")
      ax.set_xlim([0,SIZE])
      ax.set_xlabel("$r$ (M)")
      ax.axvline(AT_R, color='k')
      
      if "SANE" in run_name:
        ax.set_ylim([1e-4,1e-1])
        ax.set_yscale('log')
      else:
        ax.set_ylim([1,100])
        #ax.set_yscale('log')
      
      ax.legend(loc='upper right')
      
      ax = plt.subplot(gs[1,1])
      ax.plot(avg['r'], avg['Lj_sigma1_rt'][n], label=r"$L_{jet}$ (sigma > 1 cut)", color='C2')
      ax.plot(avg['r'], avg['Lj_Be_nob0_rt'][n], label=r"$L_{jet}$ ($Be > 0.02$ cut)", color='C3')
      ax.plot(avg['r'], avg['Lj_Be_nob1_rt'][n], label=r"$L_{jet}$ ($Be > 1.0$ cut)", color='C4')
      ax.plot(avg['r'], avg['Lj_bg1_rt'][n], label=r"$L_{jet}$ ($\beta\gamma > 1.0$ cut)", color='C5')
      ax.plot(avg['r'], avg['Lj_allp_rt'][n], label=r"$L_{tot}$", color='C6')
      
      ax.set_title(r"$L_{tot} = \int (-{T^r}_t - \rho u^r) \sqrt{-g} dx^{\theta} dx^{\phi}$")
      ax.set_xlim([0,SIZE])
      ax.set_xlabel("$r$ (M)")
      ax.axvline(AT_R, color='k')
      
      if "SANE" in run_name:
        ax.set_ylim([1e-4,1e-1])
        ax.set_yscale('log')
      else:
        ax.set_ylim([1,100])
        #ax.set_yscale('log')
      ax.legend(loc='lower right')
  
  else: # All other movie types share a layout
    ax_slc = lambda i: plt.subplot(2, 4, i)
    ax_flux = lambda i: plt.subplot(4, 2, i)
    if movie_type == "traditional":
      # Usual movie: RHO beta fluxes
      # CUTS
      bplt.plot_slices(ax_slc(1), ax_slc(2), geom, dump, np.log10(dump['RHO']),
                       label=r"$\log_{10}(\rho)$", vmin=-3, vmax=2, cmap='jet')
      bplt.plot_slices(ax_slc(5), ax_slc(6), geom, dump, np.log10(dump['beta']),
                       label=r"$\beta$", vmin=-2, vmax=2, cmap='RdBu_r')
      # FLUXES
      bplt.diag_plot(ax_flux(2), diag, 'mdot', dump['t'], ylabel=r"$\dot{M}$", logy=LOG_MDOT)
      bplt.diag_plot(ax_flux(4), diag, 'phi_b', dump['t'], ylabel=r"$\phi_{BH}$", logy=LOG_PHI)
      # Mixins:
      # Zoomed in RHO
      bplt.plot_slices(ax_slc(7), ax_slc(8), geom, dump, np.log10(dump['RHO']),
                       label=r"$\log_{10}(\rho)$", vmin=-3, vmax=2, window=[-10,10,-10,10], field_overlay=False)
      # Bsq
#       bplt.plot_slices(ax_slc[6], ax_slc[7], geom, dump, np.log10(dump['bsq']),
#                        label=r"$b^2$", vmin=-5, vmax=0, cmap='Blues')
      # Failures: all failed zones, one per nonzero pflag
#       bplt.plot_slices(ax_slc[6], ax_slc[7], geom, dump, dump['fail'] != 0,
#                        label="Failed zones", vmin=0, vmax=20, cmap='Reds', int=True) #, arrspace=True)
      # 2D histograms
#       bplt.hist_2d(ax_slc[6], np.log10(dump['RHO']), np.log10(dump['UU']),r"$\log_{10}(\rho)$", r"$\log_{10}(U)$", logcolor=True)
#       bplt.hist_2d(ax_slc[7], np.log10(dump['UU']), np.log10(dump['bsq']),r"$\log_{10}(U)$", r"$b^2$", logcolor=True)

      # Extra fluxes:
#       bplt.diag_plot(ax_flux[1], diag, dump, 'edot', r"\dot{E}", logy=LOG_PHI)
    elif movie_type == "e_ratio":
      # Energy ratios: difficult places to integrate, with failures
      bplt.plot_slices(ax_slc(0), ax_slc(1), geom, dump, np.log10(dump['UU']/dump['RHO']),
                       label=r"$\log_{10}(U / \rho)$", vmin=-3, vmax=3, average=True)
      bplt.plot_slices(ax_slc(2), ax_slc(3), geom, dump, np.log10(dump['bsq']/dump['RHO']),
                       label=r"$\log_{10}(b^2 / \rho)$", vmin=-3, vmax=3, average=True)
      bplt.plot_slices(ax_slc(4), ax_slc(5), geom, dump, np.log10(1/dump['beta']),
                       label=r"$\beta^{-1}$", vmin=-3, vmax=3, average=True)
      bplt.plot_slices(ax_slc(6), ax_slc(7), geom, dump, dump['fail'] != 0,
                       label="Failures", vmin=0, vmax=20, cmap='Reds', int=True) #, arrspace=True)
    elif movie_type == "conservation":
      # Continuity plots to verify local conservation of energy, angular + linear momentum
      # Integrated T01: continuity for momentum conservation
      bplt.plot_slices(ax_slc[0], ax_slc[1], geom, dump, Tmixed(dump, 1, 0),
                       label=r"$T^1_0$ Integrated", vmin=0, vmax=600, arrspace=True, integrate=True)
      # integrated T00: continuity plot for energy conservation
      bplt.plot_slices(ax_slc[4], ax_slc[5], geom, dump, np.abs(Tmixed(dump, 0, 0)),
                       label=r"$T^0_0$ Integrated", vmin=0, vmax=3000, arrspace=True, integrate=True)

      # Usual fluxes for reference
      bplt.diag_plot(ax_flux[1], diag, 'mdot', dump['t'], ylabel=r"$\dot{M}$", logy=LOG_MDOT)
      #bplt.diag_plot(ax_flux[3], diag, 'phi_b', dump['t'], ylabel=r"$\phi_{BH}$", logy=LOG_PHI)

      # Radial conservation plots
      E_r = sum_shell(geom,Tmixed(geom, dump, 0,0))
      Ang_r = sum_shell(geom,Tmixed(geom, dump, 0,3))
      mass_r = sum_shell(dump['ucon'][:,:,:,0]*dump['RHO'])

      # TODO arrange legend better -- add labels when radial/diag plotting
      bplt.radial_plot(ax_flux[3], geom, np.abs(E_r), 'Conserved vars at R', ylim=(0,1000), rlim=(0,20), arrayspace=True)
      bplt.radial_plot(ax_flux[3], geom, np.abs(Ang_r)/10, '', ylim=(0,1000), rlim=(0,20), col='r', arrayspace=True)
      bplt.radial_plot(ax_flux[3], geom, np.abs(mass_r),   '', ylim=(0,1000), rlim=(0,20), col='b', arrayspace=True)
      
      # Radial energy accretion rate
      Edot_r = sum_shell(geom, Tmixed(geom, dump,1,0))
      bplt.radial_plot(ax_flux[5], geom, np.abs(Edot_r), 'Edot at R', ylim=(0,200), rlim=(0,20), arrayspace=True)

      # Radial integrated failures
      bplt.radial_plot(ax_flux[7], geom, (dump['fail'] != 0).sum(axis=(1,2)), 'Fails at R', arrayspace=True, rlim=[0,50], ylim=[0,1000])

    elif movie_type == "floors":
      # TODO add measures of all floors' efficacy.  Record ceilings in header or extras?
      bplt.plot_slices(ax_flux[0], ax_flux[1], geom, dump['bsq']/dump['RHO'] - 100,
                       vmin=-100, vmax=100, cmap='RdBu_r')
      bplt.diag_plot(ax, diag, dump, 'sigma_max', 'sigma_max')

    elif movie_type == "luminosity":
      # TODO add measures of all floors' efficacy.  Record ceilings in header or extras?
      bplt.plot_slices(ax_flux[0], ax_flux[1], geom, dump['bsq']/dump['RHO'] - 100,
                       vmin=-100, vmax=100, cmap='RdBu_r')
      bplt.diag_plot(ax, diag, dump, 'sigma_max', 'sigma_max')

  # TODO enlarge plots w/o messing up even pixel count
  pad = 0.05
  plt.subplots_adjust(left=2*pad, right=1-2*pad, bottom=pad, top=1-pad)

  plt.savefig(imname, dpi=1920/FIGX)
  plt.close(fig)
  
  dump.clear()
  del dump

if __name__ == "__main__":
  # PROCESS ARGUMENTS
  if sys.argv[1] == '-d':
    debug = True
    movie_type = sys.argv[2]
    path = sys.argv[3]
    if len(sys.argv) > 4:
      tstart = float(sys.argv[4])
    if len(sys.argv) > 5:
      tend = float(sys.argv[5])
  else:
    debug = False
    movie_type = sys.argv[1]
    path = sys.argv[2]
    if len(sys.argv) > 3:
      tstart = float(sys.argv[3])
    if len(sys.argv) > 4:
      tend = float(sys.argv[4])
  
  # LOAD FILES
  files = io.get_dumps_list(path)
  if len(files) == 0:
      util.warn("INVALID PATH TO DUMP FOLDER")
      sys.exit(1)

  frame_dir = "frames_"+movie_type
  util.make_dir(frame_dir)

  hdr = io.load_hdr(files[0])
  geom = io.load_geom(hdr, path)

  jmin, jmax = get_j_vals(geom)
  #print("jmin: {} jmax: {}".format(jmin, jmax))

  if diag_post:
    # Load fluxes from post-analysis: more flexible
    diag = pickle.load(open("eht_out.p", 'rb'))
  else:
    # Load diagnostics from HARM itself
    diag = io.load_log(path)

  nthreads = util.calc_nthreads(hdr, pad=0.3)
  if debug:
    # Run sequentially to make backtraces work
    for i in range(len(files)):
      plot(i)
  else:
    util.run_parallel(plot, len(files), nthreads)
