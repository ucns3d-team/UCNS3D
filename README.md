<p align="center">
<img width="1000" height="500" src="docs/ucns3d.png">
</p>



# UCNS3D: The Open-Source High-Order Finite-Volume CFD Solver


## Overview


This repository contains the source code [ucns3d](https://ucns3d.com/) 
Computational Fluid Dynamics (CFD) solver and instructions on how to use it with representative examples.

UCNS3D is an open-source computational solver for compressible flows on unstructured meshes. State-of-the-art high-order methods are are available in a versatile 2D and 3D unstructured CFD framework for a wide-range of compressible flow problems.

The 2022 ["whitepaper"](docs/whitepaper-2022.pdf) contains a detail overview of the methods, capabilities and application of the solver.

## Overview


The `ucns3d` solver is portable to `Linux` operating systems and `MacOS` (Catalina or newer) as well as on `Windows 10` using Windows Subsystem for Linux (`WLS-2`). The parallel CFD code employing MPI+OpenMP for distributed memory machines.


ucns3d can deal with triangular, quadrilateral,tetrahedral,pyramidal, prismatic, and hexahedral elements.



## Install

There are two main methods to locally install the solver and run the solver.
* Install through Docker
* Install all dependencies and compilers manually.

### Docker

1. Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) on any operating system you are working. For Windows 10 or 11 you would need WSL ideally WSL2 installed, follow Microsoft's [instructions](https://docs.microsoft.com/en-us/windows/wsl/install).

2. In your bash terminal build the ucns3d image, you would need invoke `docker build` from the repository root directory:

```
docker build . -t ucns3d -f Dockerfile
```

3. Once the image is build run the image, you can run the image interactively like so:

```
docker run -ti ucns3d
```

The current [Dockerfile](Dockerfile) contains an example case under [tests](/tests/execute-tests.sh). Alternatively, you can mount a tmp directory and copy other uses cases when you run the image like so:

```
docker run -v $PWD/tmp/:/tmp/ -ti ucns3d
```

### Manual Build

The [source code](/src/) is written in Fortran and can be compiled in various environment. The following OS options are available.

* Linux x86-64 (Tested on Redhat, Ubuntu, Centos, Suse).
* MacOS Intel based (Catalina and newer).
* Windows 10 through WSL2, Ubuntu 20.04

### Build

The source code requires compilation and linking both static and dynamic libraries. The source code can be compiled with the following compilers.

* Intel Parallel Studio Version 17 or newer.
* G-Fortran and gcc with an MPI distribution.

### Compilation Dependencies

The solver makes use of BLAS libraries that are required for the compilation and running of the solver.

* Intel MKL library
* OpenBLAS

The mesh is partitioned using metis software.

* Open a terminal window in the [src](/src/) directory
* Ensure that you have selected the desirable compiling options in the Makefile, by specifying the appropriate fortran compiler
(ftn, ifort etc). Always compile with full debug options when developing something new, and then proceed to the more optimised
compiler options.

* For a clean installation 
```
make -f Makefile clean all
```
* For recompiling changed files and their dependencies
```
make -f Makefile
```
* the name of the executable is `ucns3d_p`.


## Running


For running ucns3d you will need the following files in a directory of your choice:
* a grid file generated with any software packages exported in Ansys fluent format (ASCII *.msh extension), given the name grid.msh or their translated to native format files GRID.cel, GRID.vrt, GRID.bnd (a detailed description of the files can be found in FILES.md)
* the UCNS3D.DAT parameter file responsible for all the settings of the code (details for the parameters of this can be found in PARAMETERS.md file)
* the executable ucns3d_p
* For interactively running the code specify the number of threads to be used by typing in the terminal window
```
export OMP_NUM_THREADS=N
```
N being the number of threads to be used (use 1 for MPI only mode)
* in the same terminal window run the code by typing
```
mpirun -np M ./ucns3d_p
```
M being the number of MPI processes (at least 2 are required), for running at different HPC systems sample [scripts](/scripts) and [libraries](/bin/lib) are provided.


## Visualisation of outputs


The solver outputs to different formats enable post-processing and visualisation through the folowing software
* [Tecplot](https://www.tecplot.com/)
* [Paraview](https://www.paraview.org/)
* [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)


## Examples


Representative tests can be downloaded from

[tests1](https://doi.org/10.5281/zenodo.3375432)

[tests2](https://doi.org/10.5281/zenodo.6538622)

and a detailed description is provided in the file TESTS.md


## License


The `ucns3d` solver is distributed under the GNU General Public Licence v3
See the LICENSE file for details.

## Support


Please get in touch and let us know how we can make this project better ucns3d@gmail.com

