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

## Running a Fishbone-Moncrief torus

The [Fishbone-Moncrief](https://doi.org/10.1086/154565)(FM) torus is the ubiquitous initial condition for modelling compact radio sources such as M87* and SgrA*. The FM problem can be simulated on `iharm3d` by using the `-PROB=torus` flag while making the program, and specifying the problem-specific parameters (compile-time and run-time) in `parameters.h` (in the build archive) and an additional parameter file which, by default `iharm3d` assumes to be named `param.dat`. The following steps outline the commands to compile and execute the problem. They presuppose that all the necessary dependencies (eg: OpenMP, MPI, phdf5, GSL) are installed and the directory variables and flags in `makefile` are pointing to them correctly.

1. Invoke the make command from the output directory,

```bash
$ make -f IHARM3D_DIRECTORY/makefile PROB=torus
```
where IHARM3D_DIRECTORY is the path to your local `iharm3d` repository that contains the makefile. The output directory is where, as explained in the section  above, the harm executable is created along with the `build_archive`. `build_archive` contains the source files necessary to run `iharm3d` along with the problem-specific compile-time parameter file, `parameters.h` and problem initialization file, `problem.c`.

2. Modify compile-time parameters in `build_archive/parameters.h`. These typically include (i) the grid size `NiTOT`; (ii) number of MPI ranks `NiCPU`; (iii) density and internal energy floors, `BSQORHOMAX`, `UORHOMAX`, `BSQOUMAX`; (iv) the reconstruction scheme `RECONSTRUCTION`. NOTE: If you're running `iharm3d` on your local system, it is recommended that the FM problem is run at a low resolution or a 2D problem is executed (set `N3TOT` to 1).

3. If the compile-time parameters have been modified or the C code in any of the source files in `build_archive` has been edited, the harm executable must be remade with the same command as in (1) from the output directory.

4. Copy any of the parameter files located at `IHARM3D_DIRECTORY/prob/torus/` labelled param_sane.dat or param_mad.dat to the output directory and rename the file as `param.dat`. This contains the runtime parameters for the FM torus (eg: duration of run, domain size, output file cadence, fluid properties, FM torus size, FMKS grid geometry). NOTE: It is again recommended to set `tf` to a reasonable value if you're running the problem on your local computer.

5. Submitting the run: Once the runtime parameters have been updated, you're good to run the FM problem. The command to launch the run depends on the capabilities of your system,
   (i) If you're executing the problem on your local system or a remote computer with a single node, you do not need the MPI dependency and following command should suffice (run from output directory),
   
   ```bash
   $ ./harm -p param.dat >LOG_FILE
   ```
   where the runtime log is redirected to `LOG_FILE`. If `STDOUT` is not redirected, the runtime log will be printed on the terminal. NOTE: You can set the number of cores over which you want `iharm3d` to execute by modifying the environment variable, `OMP_NUM_THREADS` pre-compilation. If not provided, `iharm3d` by default is run across all cores available.
   
   (ii) If you're running the problem on a HPCC, you may have multiple nodes at your disposal. In this case, you can use an MPI implementation to parallelize the jobs across several nodes. The exact command to launch `harm` depends on the MPI implementation. If you arer running `iharm3d` on a TACC system (which has the SLURM job scheduler), you may find the various job submission scripts located at `IHARM3D_DIRECTORY/scripts/submit` useful. You can submit the job on any TACC machine as,
   
   ```bash
   $sbatch -N (NODES) -p (QUEUE) IHARM3D_DIRECTORY/scripts/submit/SUBMIT_SCRIPT.sb
   ```
   where `SUBMIT_SCRIPT.sb` is the job submission script that varies in accordance with the TACC system you're logged into. 

## Hacking
Notes that may save you time in reading the source code:
* Grid coordinates match physical coordinates i => x^1, j => x^2, k => x^3.  However, they are indexed backward
in memory `grid[k][j][i]`.  A number of loop aliases e.g. `ILOOP`, `ZLOOP` are defined to make this counter-intuitive ordering,
as well as the presence of border "ghost" zones, easier to manage
* The fluid state `S` is often modified in-place.  Rest assured the accompanying grid `G` is not.  Both are structs of arrays,
given `typedef`s in order to allocate their backing memory contiguously
* Comments are sparse, and usually concern implementation details, not algorithmic operation. See
[iharm2d_v3](https://github.com/AFD-Illinois/iharm2d_v3) for a simpler version which may prove a gentler introduction.
