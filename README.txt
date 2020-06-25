UCNS3D CODE Instructions.


1. Requirements

a) Linux x86-64 (Tested on Redhat, Ubuntu, Centos, Suse) or MacOS (Catalina and newer)
b) Intel Parallel Studio Version 17 or newer, or (gfotran and gcc with an MPI distribution)
c) Intel MKL library (or OpenBLAS)
d) Tecplot, Paraview or VisIt for visualisation





2. Compiling

a) Open a terminal window in the CODE directory
b) Ensure that you have selected the desirable compiling options in the Makefile, by specifying the correct fortran compiler
(ftn, ifort etc). Always compile with full debug options when developing something new, and then proceed to the more optimised
compiler options
c) Type "make -f Makefile clean all" for a clean installation or "make -f Makefile" for recompiling the changed files and their
dependencies
d) the name of the executable is ucns3d_p


3. Grids
a) The grids 2d or 3d can be generated in any grid generation software package (such as ICEM-CFD, Pointwise, Gridgen, Gambit..)
that can export the grid and the boundary conditions in the Ansys fluent format (ASCII *.msh extension).
b) The boundary conditions codes need to be applied in the grid generation process 
(apply boundaries such as wall, outflow, inflow, symmetry, velocity inlet, pressure farfield etc). 
The actual codes and their corresponding boundary condition that are used in the code can be seen in the translate.f90 module.
c) The grid must be exported with the name grid.msh.
d) Once the code starts the grid.msh file format is translated to the native format of the ucns3d code 
where three files are created namely GRID.cel, GRID.vrt, GRID.bnd.



4. Running

Representative tests can be downloaded from 
https://doi.org/10.5281/zenodo.3375432

For running ucns3d you will need the following files in a directory of your choice:
a) a grid file generated with any software packages exported in Ansys fluent format (ASCII *.msh extension), 
given the name grid.msh or their translated to native format files GRID.cel, GRID.vrt, GRID.bnd
b)the UCNS3D.DAT parameter file responsible for all the settings of the code (more on README_UCNS3D.txt file)
c) the executable ucns3d_p
d) For interactively running the code specify the number of threads to be used by typing in the terminal window:
    "export OMP_NUM_THREADS=N", N being the number of threads to be used (use 1 for MPI only mode)
    in the same terminal window run the code by typing:
    "mpirun -np M ./ucns3d_p",  M being the number of MPI processes (at least 2 are required)
e) For running at different HPC systems sample scripts are provided

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!
The structure of the mesh files, solution files, restart files etc, are 
independent of the number of cpus selected, since parmetis is invoked 
within the code to partition the mesh, and MPI-IO to put the files together
when reading and writing.
Therefore the user can submit to any number of cpus a simulation, and 
restart/resume it to a different number of cpus. There is only one file 
written from all the cpus for the output, restart files and forces etc,
therefore minimising the number of files generated.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


5.Output Files
The ucns3d code generates numerous files while running depending on the setup chosen. The most common ones are:

a)history.txt: outlines the various stages of the initialisation and indicates the timestep size, iteration, 
time, and message when writing solution or checkpoint files
b)STATISTICS.txt: outlines the cpu time taken per iteration  for various procedures in the code while the code is running,
including cputime taken for reconstruction, fluxes, solution update, mpi boundary values exhange, mpi hallo exchange for stencils etc. Note that the collection of statistics has a noticeable computational footpring and should only be used for a limited number of iterations rather than production runs for the entire simulation.
c) FORCE.dat : specifies in four columns the (iteration, time, Lift Coefficient, Drag Coefficient). You need to ensure that
the last two are normalised with respect to the reference area of your geometry.
d) RESTART.dat: A checkpoint file written when the simulation is stopped because one of the four criteria has been met 
(the maximum wall clock time limit, maximum number of iterations, residual tolerance, time of simulation). You need to ensure
that when resubmitting a simulation for reaching a new time, number of iterations or for simply restarting from an existing 
solution, you need to ensure that this file is included in the directory that you run the executable from.
e) RESTART_AV.dat: A checkpoint file written when performing the time averaging in the code. You need to ensure that this file
is included when resubmitting an unsteady simulation and you want to keep performing the time averaging from this checkpoint.
f) GRID.plt, SURF.plt : Tecplot binary or ascii files with the geometry of the volume mesh and surface mesh. (exported only 
tecplot output is selected)
g) OUT_*.plt SURF_*.plt: Tecplot binary or ascii files for the volume and surface mesh instantaneous solution respectively.
This way only the solution files are written rather the entire geometry in order to save resources.(exported only tecplot 
binary output is selected)
h) OUT_*.vtk, SURF_*.vtk: Paraview binary or aschii files for the both the geometry and instantaneous solution for the 
volume and surface mesh respectively.
i) VOL_AVER_*.plt SURF_AV*.plt: Tecplot binary or ascii files for the volume and surface mesh AVERAGE solution respectively. 
This way only the solution files are written rather the entire geometry in order to save resources.(exported only tecplot 
binary output is selected)
j) OUT_AV*.vtk, SURF_AV*.vtk: Paraview binary or aschii files for the both the geometry and average solution for the volume 
and surface mesh respectively.
k) PROBE.* : file outlining the primitive variables variation with time (column 1 time, column 2 density, column 3 u velocity
etc) for the specified probe locations.
l) residual.dat : file containing the residuals of the conserved flow variables and turbulence flow variables at every 10 
iterations for steady state flow problems
m) ENERGY.dat: file containing the time and the normalised total kinetic energy, and kinetic energy dissipation rate for
the Taylor Green vortex test problem.

All the solution files carry timestamps and description of the flow variable names, for both tecplot and paraview outputs.










