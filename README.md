# iharm3D
This code implements the HARM algorithm outlined in [Gammie et al. 2003](https://doi.org/10.1086/374594), with some modifications
outlined in [McKinney & Gammie 2004](https://doi.org/10.1086/422244).  This is a second-order, conservative, shock-capturing scheme 
for general-relativistic magnetohydrodynamics (GRMHD).

## Requirements
`iharm3D` requires an MPI/Parallel HDF5 stack.  In practice, this means that the executable `h5pcc` must be in your `PATH`
and working correctly.  It also requires the GNU Scientific Library (GSL) to be installed or loaded.
Most Linux distributions package these requirements, and most supercomputers have modules for them.

## Building
Provided the above is met, `iharm3D` can be built with
```bash
$ make
```
This builds the 'torus' problem, `prob/torus`, which sets up the initial conditions in the equilibrium torus configuration of
[Fishbone & Moncrief 1976](https://doi.org/10.1086/154565).  This is by far the most common problem used for science runs.

Alternative (mostly testing) problems can be built by specifying their folder name in `prob/`, e.g.
```bash
$ make PROB=mhdmodes
```
Refer to existing problems and/or forthcoming developer documentation for details on how to add new problem definitions.

The make process or flags can be customized by adding a host-specific makefile
```bash
$ touch machines/$(hostname).make
```
which can contain any valid `make` script and is read after setting most of the default parameters, in order to override them

## Running
`iharm3D` uses both compile-time and runtime parameters, given in the problem directories as `parameters.h` and
`param.dat` respectively.  A general workflow for customization and initiation of a run is in
`script/submit/checklist.txt`.  Most problem-specific and physical parameters are specified at runtime, while the MPI process
geometry and certain operations/stability options are specified at compile time.

## Hacking
Notes that may save you time in reading the source code:
* Grid coordinates match physical coordinates i => x^1, j => x^2, k => x^3.  However, they are indexed backward
in memory `grid[k][j][i]`.  A number of loop aliases e.g. `ILOOP`, `ZLOOP` are defined to make this counter-intuitive ordering,
as well as the presence of border "ghost" zones, easier to manage
* The fluid state `S` is often modified in-place.  Rest assured the accompanying grid `G` is not.  Both are structs of arrays,
given `typedef`s in order to allocate their backing memory contiguously
* Comments are sparse, and usually concern implementation details, not algorithmic operation. See
[iharm2d_v3](https://github.com/AFD-Illinois/iharm2d_v3) for a simpler version which may prove a gentler introduction.
