---
title: 'iharm3D: Vectorized General Relativistic Magnetohydrodynamics'
tags:
  - C
  - magnetohydrodynamics
  - general relativity
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Ben S. Prather^[Corresponding author]
    orcid: 0000-0002-0393-7734
    affiliation: "1, 2"
  - name: George N. Wong
    orcid: 0000-0001-6952-2147
    affiliation: "1, 2"
  - name: Vedant Dhruv
    orcid: 0000-0001-6765-877X
    affiliation: "1, 2"
  - name: Benjamin R. Ryan
    orcid: 0000-0001-8939-4461
    affiliation: 3
  - name: Joshua C. Dolence
    orcid: 0000-0003-4353-8751
    affiliation: 3
  - name: Sean M. Ressler
    orcid: 0000-0003-0220-5723
    affiliation: 5
  - name: Charles F. Gammie
    orcid: 0000-0001-7451-8935
    affiliation: "1, 2, 4"
affiliations:
 - name: Physics Department, University of Illinois at Urbana--Champaign, 1110 West Green Street, Urbana, IL 61801, USA 
   index: 1
 - name: Illinois Center for Advanced Studies of the Universe
   index: 2
 - name: CCS-2, Los Alamos National Laboratory, P.O. Box 1663, Los Alamos, NM 87545, USA
   index: 3
 - name: Astronomy Department, University of Illinois at Urbana--Champaign, 1002 West Green Street, Urbana, IL 61801, USA
   index: 4
 - name: Kavli Institute for Theoretical Physics, University of California Santa Barbara, Kohn Hall, Santa Barbara, CA 93107, USA
   index: 5

date: 12 April 2021
bibliography: paper.bib

---
# `iharm3D` Functionality and Purpose

`iharm3D`^[https://github.com/AFD-Illinois/iharm3d] is an open-source C code for simulating black hole accretion systems in arbitrary stationary spacetimes using ideal general-relativistic magnetohydrodynamics (GRMHD). It is an implementation of the HARM ("High Accuracy Relativistic Magnetohydrodynamics") algorithm outlined in @gammie_harm:_2003 with updates as outlined in @mckinney_measurement_2004 and @noble_primitive_2006.  The code is most directly derived from @ryan_bhlight:_2015 but with radiative transfer portions removed.  HARM is a conservative finite-volume scheme for solving the equations of ideal GRMHD, a hyperbolic system of partial differential equations, on a logically Cartesian mesh in arbitrary coordinates.

# Statement of Need

Numerical simulations are crucial in modeling observations of active galactic nuclei, such as the recent horizon-scale results from the Event Horizon Telescope (EHTC) and GRAVITY collaborations.  The computational simplicity of ideal GRMHD enables the generation of long, high-resolution simulations and broad parameter-exploration studies that can be compared to observations for parameter inference.

Multiple codes already exist for solving the ideal GRMHD equations on regular Eulerian meshes in 3D.  Some examples of codes currently in use:

- Athena++ (@stone_athena_2020, @white_extension_2016)
- BHAC (@porth_black_2017)
- Cosmos++ (@anninos_cosmos_2005, @fragile_numerical_2012, @fragile_numerical_2014)
- ECHO (@londrillo_high-order_2000, @londrillo_divergence-free_2004)
- H-AMR (@liska_h-amr_2019, @liska_large-scale_2020)
- HARM-Noble (@noble_primitive_2006, @noble_direct_2009, @noble_circumbinary_2012, @zilhao_dynamic_2014, @bowen_quasi-periodic_2018)
- IllinoisGRMHD (@etienne_illinoisgrmhd_2015)
- KORAL (@sadowski_semi-implicit_2013, @sadowski_numerical_2014)
- GRHydro (@mosta_grhydro_2014)
- Spritz (@cipolletta_spritz_2020, @cipolletta_spritz_2021)

As the length of this list illustrates, the field of GRMHD simulation is now well established, and many codes now exist to serve different needs.  These codes can be distinguished by the trade-offs they make in prioritizing speed, simplicity, and generality, with the latter encompassing, e.g., support for dynamical spacetimes, adaptive mesh refinement, or higher-order integration schemes.

In particular, `iharm3D` development focuses on providing the simplest and fastest possible code capable of simulating the original systems of interest when designing HARM, even at the cost of features aimed at more general applicability.  It provides a fast and scalable modern implementation of HARM, but maintains the conventions and simple structure of the original described in @gammie_harm:_2003.  The result is a code relatively easy to understand and modify, yet capable of running simulations at state-of-the-art scale.

# Implementation Notes

In MHD, uncorrected discretization errors inevitably lead to violations of the no-monopoles condition $\nabla \cdot B = 0$.  As in the original HARM implementation, `iharm3D` uses the "Flux-CT" scheme for cell-centered constrained transport outlined in @toth_b0_2000.

`iharm3D` also retains numerical evaluation of all metric-dependent quantities, allowing trivial modification of the coordinate system or background spacetime so long as the line element is available in analytic form.  This can be used as a form of static mesh refinement, since the coordinates can be adapted to place resolution in areas of interest (e.g., an accretion disk midplane).

In GRMHD, "conserved" variables (energy and momentum densities) are complicated analytic functions of "primitive" variables (density, pressure, and velocity).  Conserved variables are stepped forward in time and then inversion to primitives is done numerically. `iharm3d` uses the "$1D_W$" scheme outlined in @noble_primitive_2006.  As the equations of ideal GRMHD are rescalable, any consistent set of units may be chosen to evolve them in the code.  For numerical stability, we choose units in which $GM = c = 1$, and scale the density such that the densest zone in the initial conditions has $\rho = 1$.  The mass and length are rescaled in post-processing to reflect the particular system under study.

To model a collisionless plasma, `iharm3D` implements an optional scheme that provides a means of tracking and partitioning dissipation into ions and electrons (@ressler_electron_2015). Currently the code implements five heating models: those of @howes_prescription_2010, @kawazura_thermal_2019, @werner_non-thermal_2018, @rowan_electron_2017, and @sharma_electron_2007.

To avoid catastrophic failures caused by discretization error, especially in low density regions, fluid variables are bounded at the end of each step. Typical `iharm3D` bounds in black hole accretion problems are enforced as follows:

- Density $\rho > 10^{-6} k$, for $k \equiv \frac{1}{r^2 (1 + r/10)}$, with radius $r$ in units of gravitational radius $r_g$ of the central object,
- Internal energy density $u > 10^{-8} k^{\gamma}$ where $\gamma \equiv$ adiabatic index,
- $\rho$ and $u$ are incremented until $\sigma \equiv \frac{2 P_b}{\rho} < 400$ and $\beta \equiv \frac{P_{gas}}{P_b} > 2.5 \times 10^{-5}$ where $P_b \equiv \frac{b^2}{2}$ is the magnetic pressure,
- $\rho$ is incremented until $\frac{u}{\rho} < 100$,
- When evolving electron temperatures, $u$ is decremented until $\frac{P_{gas}}{\rho^{\gamma}} < 3$,
- Velocity components are downscaled until Lorentz factor $\Gamma \equiv \frac{1}{\sqrt{1 - v^2}} < 50$.

Global disk simulations inevitably invoke these bounds, most frequently those on $\sigma$ and $\Gamma$.

# Tests

The convergence properties of HARM are well-studied in @gammie_harm:_2003.  `iharm3D` implements most of the tests presented in that paper as integration and regression tests.  Figure \ref{fig:convergence} shows convergence results for linear modes and for un-magnetized Bondi flow.

![Results of convergence tests with `iharm3d`'s main branch, plotting L1 norm of the difference between the computed solution and the analytic or stable result with increasing domain size.  Wave solutions were performed on a 3D cubic grid N zones to one side, the Bondi accretion problem was performed on a logically Cartesian 2D square grid N zones on one side. \label{fig:convergence}](figures/convergence.pdf)

`iharm3D` implements three additional tests which check that fluid evolution is identical under different domain decompositions: one which initializes a new fluid state, one which restarts from a checkpoint file, and one comparing the initialized state to an equivalent checkpoint file.

# Scaling

Key `iharm3D` routines are highly vectorized and have efficient memory access patterns. Originally developed for Intel Knights Landing (KNL) chips on the Stampede2 supercomputer at Texas Advanced Computing Center (TACC), `iharm3D` also runs efficiently on TACC's Frontera CPU nodes.

Figure \ref{fig:scaling} presents scaling results for `iharm3D` on both Stampede2 and Frontera.

![Strong scaling performance of iharm3D. Performance is measured in zones advanced by one cycle each second (Zone-Cycles per Second), when a problem with $256^3$ zones is split among N nodes \label{fig:scaling}](figures/scaling.pdf){ width=3in }

# Research projects using iharm3D

`iharm3D` is one of several GRMHD codes used by the EHT Collaboration to produce its library of fluid simulations. Images produced from this library were used for validation tests in @PaperIV and @PaperVII and for interpretation of the M87 EHT results in total intensity (@PaperV, @PaperVI) and polarization (@PaperVIII).

Papers making use of the results of `iharm3D` simulations include @porth_event_2019, @johnson_universal_2020, @gold2020, @palumbo_discriminating_2020, @lin_feature_2020, @ricarte_decomposing_2020, @Wielgus_monitoring_2020, @tiede_variational_2020, and @gelles_role_2021.

# Acknowledgements

This work was supported by National Science Foundation grants AST 17-16327, OISE 17-43747, AST 20-07936, AST 20-34306, and PHY 17-48958, by a Donald C. and F. Shirley Jones Fellowship to G.N.W., by the Gordon and Betty Moore Foundation through Grant GBMF7392, and by the US Department of Energy through Los Alamos National Laboratory. Los Alamos National Laboratory is operated by Triad National Security, LLC, for the National Nuclear Security Administration of the US Department of Energy (Contract No. 89233218CNA000001).  This work has been assigned a document release number LA-UR-21-23714.

This work used the Extreme Science and Engineering Discovery Environment (XSEDE) resource Stampede2 at the Texas Advanced Computing Center, which is supported by National Science Foundation grant number ACI-1548562.  The authors acknowledge the Texas Advanced Computing Center (TACC) at The University of Texas at Austin for providing HPC resources that have contributed to the research results reported within this paper.

# References
