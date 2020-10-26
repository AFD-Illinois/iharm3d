# iharm3D
This code implements the HARM algorithm outlined in [Gammie et al. 2003](https://doi.org/10.1086/374594), with some modifications
outlined in [McKinney & Gammie 2004](https://doi.org/10.1086/422244).  This is a second-order, conservative, shock-capturing scheme 
for general-relativistic magnetohydrodynamics (GRMHD).  Credit also to the many people who have worked on the code over the years, 
including Scott Noble, who implemented the first 3D version of the code, Josh Dolence, Ben Ryan, George Wong, and Ben Prather.

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
which can contain any valid `make` script and is read after setting most of the default parameters, in order to override them.

## Configuration and Running
Building `iharm3d` produces a directory named `build_archive` in the directory where `make` is invoked.  This archive contains
all the source files used in the build, as well as all the object files and a copy of the executable.

If `build_archive` already exists, `make` will prefer any newer/modified files in that directory, vs their equivalents in the original source.
This allows modifying the compile-time parameters in `parameters.h`, or even modifying the C code as desired, without disrupting the
original repository and potentially committing upstream whatever compile-time or runtime configuration you happen to be using.

Note that `iharm3d` also takes runtime parameters (most of the physical parameters, whereas grid size & MPI topology are compile-time).
`iharm3d` will automatically use any file called `param.dat` in the current working directory, and will output simulation data to the
working directory as well.  You can specify an alternative parameter file with `-p` or output directory with `-o`.  Sample runtime
parameters for each problem are provided in the problem directories.

Due to this extra copy, note that between building different problems (e.g. from a torus to the MHD modes problem) one must run
```bash
$ make distclean
```
which will remove the `build_archive` directory, including any customizations that had previously been applied. A simple `make clean` will
remove just the object and executable files, preserving any customizations in `build_archive`

Full details of production runs on larger machines e.g. Stampede2 are in `script/submit/checklist.txt` in this repository, along
with job submission scripts for SLURM in the TACC environment, adaptable for a lot of SLURM machines.

## Hacking
Notes that may save you time in reading the source code:
* Grid coordinates match physical coordinates i => x^1, j => x^2, k => x^3.  However, they are indexed backward
in memory `grid[k][j][i]`.  A number of loop aliases e.g. `ILOOP`, `ZLOOP` are defined to make this counter-intuitive ordering,
as well as the presence of border "ghost" zones, easier to manage
* The fluid state `S` is often modified in-place.  Rest assured the accompanying grid `G` is not.  Both are structs of arrays,
given `typedef`s in order to allocate their backing memory contiguously
* Comments are sparse, and usually concern implementation details, not algorithmic operation. See
[iharm2d_v3](https://github.com/AFD-Illinois/iharm2d_v3) for a simpler version which may prove a gentler introduction.