# Files Description

General
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

history.txt
-----------------------------------------------------
outlines the various stages of the initialisation and indicates the timestep size, iteration, 
time, and message when writing solution or checkpoint files.

STATISTICS.txt
-----------------------------------------------------
outlines the cpu time taken per iteration  for various procedures in the code while the code is running,
including cputime taken for reconstruction, fluxes, solution update, mpi boundary values exhange, mpi hallo exchange for stencils etc. 
Note that the collection of statistics has a noticeable computational footpring and should only be used for a limited number of 
iterations rather than production runs for the entire simulation.

FORCE.dat 
-----------------------------------------------------
specifies in four columns the (iteration, time, Lift Coefficient, Drag Coefficient). You need to ensure that
the last two are normalised with respect to the reference area of your geometry.

RESTART.dat
-----------------------------------------------------
A checkpoint file written when the simulation is stopped because one of the four criteria has been met 
(the maximum wall clock time limit, maximum number of iterations, residual tolerance, time of simulation). You need to ensure
that when resubmitting a simulation for reaching a new time, number of iterations or for simply restarting from an existing 
solution, you need to ensure that this file is included in the directory that you run the executable from.

RESTART_AV.dat
-----------------------------------------------------
A checkpoint file written when performing the time averaging in the code. You need to ensure that this file
is included when resubmitting an unsteady simulation and you want to keep performing the time averaging from this checkpoint

GRID.plt, SURF.plt
-----------------------------------------------------
Tecplot binary or ascii files with the geometry of the volume mesh and surface mesh. (exported only 
tecplot output is selected)

OUT_*.plt SURF_*.plt
-----------------------------------------------------
Tecplot binary or ascii files for the volume and surface mesh instantaneous solution respectively.
This way only the solution files are written rather the entire geometry in order to save resources.(exported only tecplot 
binary output is selected)

OUT_*.vtk, SURF_*.vtk
-----------------------------------------------------
Paraview binary or aschii files for the both the geometry and instantaneous solution for the 
volume and surface mesh respectively.

VOL_AVER_*.plt SURF_AV*.plt
-----------------------------------------------------
Tecplot binary or ascii files for the volume and surface mesh AVERAGE solution respectively. 
This way only the solution files are written rather the entire geometry in order to save resources.(exported only tecplot 
binary output is selected)

OUT_AV*.vtk, SURF_AV*.vtk
-----------------------------------------------------
Paraview binary or aschii files for the both the geometry and average solution for the volume 
and surface mesh respectively.

PROBE.* 
-----------------------------------------------------
file outlining the primitive variables variation with time (column 1 time, column 2 density, column 3 u velocity
etc) for the specified probe locations.

residual.dat 
-----------------------------------------------------
file containing the residuals of the conserved flow variables and turbulence flow variables at every 10 
iterations for steady state flow problems

ENERGY.dat
-----------------------------------------------------

file containing the time and the normalised total kinetic energy, and kinetic energy dissipation rate for
the Taylor Green vortex test problem.

All the solution files carry timestamps and description of the flow variable names, for both tecplot and paraview outputs.
