
# UCNS3D: The Open-Source High-Order Finite-Volume CFD Solver


General Information
===================

This directory contains the [ucns3d](https://ucns3d.com/) 
Computational Fluid Dynamics (CFD) solver.

UCNS3D is an open-source computational solver for compressible flows on unstructured meshes. 
State-of-the-art high-order methods are now be available in a versatile 2D and 3D unstructured
CFD framework for a wide-range of compressible flow problems. 

Presentation
============


ucns3d is portable to Linux operating systems and MacOS (Catalina or newer)
and it is a parallel CFD code employing MPI+OpenMP for distributed memory machines.
It is based on the least-squares high-order finite-volume methods that have been
routinely developed and tested under a series of test problems as they can be tracked in the
website of [ucns3d](https://ucns3d.com/research/) and more recently a new reference article
has been published on the [ucns3d](https://www.sciencedirect.com/science/article/pii/S0010465522001722)
that details all the methods and governing equations currently available in the framework.

ucns3d can deal with triangular, quadrilateral,tetrahedral,pyramidal, prismatic, and hexahedral elements.

Meshes can be imported using the fluent *.msh, STAR-CD *.cel, or UGRID *.ugrid file formats 
and post-processing output is avilable in [tecplot](https://www.tecplot.com/) or vtk format which can be 
used by [Paraview](https://www.paraview.org/), [Visit](https://wci.llnl.gov/simulation/computer-codes/visit) or other platforms.


License
=======

ucns3d is distributed under the GNU General Public Licence v3
See the LICENSE file for details.

Installation
============

Requirements
-----------------------------------------------------

• Linux x86-64 (Tested on Redhat, Ubuntu, Centos, Suse) or MacOS (Catalina and newer)
• Intel Parallel Studio Version 17 or newer, or (gfotran and gcc with an MPI distribution)
• Intel MKL library or OpenBLAS
• Tecplot, Paraview or VisIt for visualisation


Compiling
-----------------------------------------------------

• Open a terminal window in the CODE directory
• Ensure that you have selected the desirable compiling options in the Makefile, by specifying the correct fortran compiler
(ftn, ifort etc). Always compile with full debug options when developing something new, and then proceed to the more optimised
compiler options
• For a clean installation 
```
$ make -f Makefile clean all
```
• For recompiling changed files and their dependencies
```
$make -f Makefile
```
•the name of the executable is ucns3d_p


Running
============


For running ucns3d you will need the following files in a directory of your choice
• a grid file generated with any software packages exported in Ansys fluent format (ASCII *.msh extension), 
given the name grid.msh or their translated to native format files GRID.cel, GRID.vrt, GRID.bnd
• the UCNS3D.DAT parameter file responsible for all the settings of the code (more on README_UCNS3D.txt file)
• the executable ucns3d_p
• For interactively running the code specify the number of threads to be used by typing in the terminal window
```
$ export OMP_NUM_THREADS=N
```
N being the number of threads to be used (use 1 for MPI only mode)
• in the same terminal window run the code by typing
```
$ mpirun -np M ./ucns3d_p
```
M being the number of MPI processes (at least 2 are required)
• For running at different HPC systems sample scripts are provided


Examples
-----------------------------------------------------

Representative tests can be downloaded from

[tests1](https://doi.org/10.5281/zenodo.3375432)

[tests2](https://doi.org/10.5281/zenodo.6538622)


Files Note
==============

The structure of the mesh files, solution files, restart files etc, are 
independent of the number of cpus selected, since parmetis is invoked 
within the code to partition the mesh, and MPI-IO to put the files together
when reading and writing.
Therefore the user can submit to any number of cpus a simulation, and 
restart/resume it to a different number of cpus. There is only one file 
written from all the cpus for the output, restart files and forces etc,
therefore minimising the number of files generated.


Output files
==============

The ucns3d code generates numerous files while running depending on the setup chosen. The most common ones are:

• history.txt: outlines the various stages of the initialisation and indicates the timestep size, iteration, 
time, and message when writing solution or checkpoint files
• STATISTICS.txt: outlines the cpu time taken per iteration  for various procedures in the code while the code is running,
including cputime taken for reconstruction, fluxes, solution update, mpi boundary values exhange, mpi hallo exchange for stencils etc. 
Note that the collection of statistics has a noticeable computational footpring and should only be used for a limited number of 
iterations rather than production runs for the entire simulation.
• FORCE.dat : specifies in four columns the (iteration, time, Lift Coefficient, Drag Coefficient). You need to ensure that
the last two are normalised with respect to the reference area of your geometry.
•  RESTART.dat: A checkpoint file written when the simulation is stopped because one of the four criteria has been met 
(the maximum wall clock time limit, maximum number of iterations, residual tolerance, time of simulation). You need to ensure
that when resubmitting a simulation for reaching a new time, number of iterations or for simply restarting from an existing 
solution, you need to ensure that this file is included in the directory that you run the executable from.
•  RESTART_AV.dat: A checkpoint file written when performing the time averaging in the code. You need to ensure that this file
is included when resubmitting an unsteady simulation and you want to keep performing the time averaging from this checkpoint.
•  GRID.plt, SURF.plt : Tecplot binary or ascii files with the geometry of the volume mesh and surface mesh. (exported only 
tecplot output is selected)
•  OUT_*.plt SURF_*.plt: Tecplot binary or ascii files for the volume and surface mesh instantaneous solution respectively.
This way only the solution files are written rather the entire geometry in order to save resources.(exported only tecplot 
binary output is selected)
•  OUT_*.vtk, SURF_*.vtk: Paraview binary or aschii files for the both the geometry and instantaneous solution for the 
volume and surface mesh respectively.
•  VOL_AVER_*.plt SURF_AV*.plt: Tecplot binary or ascii files for the volume and surface mesh AVERAGE solution respectively. 
This way only the solution files are written rather the entire geometry in order to save resources.(exported only tecplot 
binary output is selected)
•  OUT_AV*.vtk, SURF_AV*.vtk: Paraview binary or aschii files for the both the geometry and average solution for the volume 
and surface mesh respectively.
•  PROBE.* : file outlining the primitive variables variation with time (column 1 time, column 2 density, column 3 u velocity
etc) for the specified probe locations.
•  residual.dat : file containing the residuals of the conserved flow variables and turbulence flow variables at every 10 
iterations for steady state flow problems
•  ENERGY.dat: file containing the time and the normalised total kinetic energy, and kinetic energy dissipation rate for
the Taylor Green vortex test problem.

All the solution files carry timestamps and description of the flow variable names, for both tecplot and paraview outputs.


Support
==============
ucns3d support: ucns3d@gmail.com

