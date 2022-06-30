
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

The requirements are:

* Linux x86-64 (Tested on Redhat, Ubuntu, Centos, Suse) or MacOS (Catalina and newer)
* Intel Parallel Studio Version 17 or newer, or (gfotran and gcc with an MPI distribution)
* Intel MKL library or OpenBLAS
* Tecplot, Paraview or VisIt for visualisation.


Compiling
-----------------------------------------------------
To compile:

* Open a terminal window in the CODE directory
* Ensure that you have selected the desirable compiling options in the Makefile, by specifying the correct fortran compiler
(ftn, ifort etc). Always compile with full debug options when developing something new, and then proceed to the more optimised
compiler options
* For a clean installation 
```
$ make -f Makefile clean all
```
* For recompiling changed files and their dependencies
```
$make -f Makefile
```
*the name of the executable is ucns3d_p.


Running
============


For running ucns3d you will need the following files in a directory of your choice:
* a grid file generated with any software packages exported in Ansys fluent format (ASCII *.msh extension), given the name grid.msh or their translated to native format files GRID.cel, GRID.vrt, GRID.bnd
* the UCNS3D.DAT parameter file responsible for all the settings of the code (details for the parameters of this can be found in PARAMETERS.md file)
* the executable ucns3d_p
* For interactively running the code specify the number of threads to be used by typing in the terminal window
```
$ export OMP_NUM_THREADS=N
```
N being the number of threads to be used (use 1 for MPI only mode)
* in the same terminal window run the code by typing
```
$ mpirun -np M ./ucns3d_p
```
M being the number of MPI processes (at least 2 are required)
* For running at different HPC systems sample scripts are provided.


Examples
-----------------------------------------------------

Representative tests can be downloaded from

[tests1](https://doi.org/10.5281/zenodo.3375432)

[tests2](https://doi.org/10.5281/zenodo.6538622)

and a detailed description is provided in the file TESTS.md


Support
==============
ucns3d support: ucns3d@gmail.com

