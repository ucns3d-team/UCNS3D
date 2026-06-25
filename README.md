<p align="center">
  <img width="1000" height="500" src="docs/ucns3d.png" alt="UCNS3D">
</p>

# UCNS3D

**UCNS3D** is an open-source high-order Computational Fluid Dynamics solver for compressible flows on unstructured meshes. It is designed for research-grade simulations on modern parallel architectures, from MPI/OpenMP CPU clusters to selected GPU/accelerator kernels through the OpenMP target/XPU build path.

UCNS3D supports two- and three-dimensional simulations on mixed-element meshes and provides a flexible framework for high-order finite-volume methods, turbulence modelling, and large-scale aerodynamic flow problems.

The 2022 [UCNS3D whitepaper](docs/whitepaper-2022.pdf) provides a detailed description of the numerical methods, solver capabilities, and representative applications.

## At a Glance

| Capability | Support |
| --- | --- |
| Flow regime | Compressible inviscid and viscous flows |
| Meshes | 2D and 3D unstructured mixed-element meshes |
| Elements | Triangles, quadrilaterals, tetrahedra, pyramids, prisms, hexahedra |
| Numerics | MUSCL and high-order finite-volume reconstruction |
| Parallelism | MPI and OpenMP |
| Acceleration | Optional OpenMP target offload, referred to in the code as the XPU path |
| Output | ParaView, Tecplot, and VisIt-compatible files |

## Main Features

- High-order compressible-flow solver for unstructured grids
- Support for mixed-element 2D and 3D meshes
- Explicit and implicit simulation workflows
- MPI domain decomposition for distributed-memory systems
- OpenMP threading for shared-memory parallelism
- Optional GPU/accelerator offload through the `xpu` preprocessor path
- Case-specific compile-time bounds for XPU work arrays
- Standard visualization output for common CFD post-processing tools

## Repository Layout

Typical repository contents include:

```text
src/          Fortran source code
docs/         documentation and whitepaper material
scripts/      example run scripts for local and HPC systems
tests/        representative test workflows
bin/lib/      optional third-party or platform-specific libraries
```

Additional documentation is provided in:

```text
FILES.md       mesh and file-format description
PARAMETERS.md  runtime parameter description
TESTS.md       example-case description
```

## Installation

UCNS3D can be used through Docker or compiled manually on the target machine.

### Docker

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/). On Windows, use Windows Subsystem for Linux 2.

Build the image from the repository root:

```bash
docker build . -t ucns3d -f Dockerfile
```

Run the container interactively:

```bash
docker run -ti ucns3d
```

To mount a local working directory:

```bash
docker run -v $PWD/tmp/:/tmp/ -ti ucns3d
```

The current [Dockerfile](Dockerfile) includes an example workflow under [tests](/tests/execute-tests.sh).

### Manual Build

The source code is written in Fortran and is intended to build on:

- Linux x86-64 systems
- macOS
- Windows through WSL2
- HPC systems using MPI compiler wrappers

Common compiler environments include:

- Intel Fortran / Intel MPI
- GNU Fortran with MPI
- Cray compiler wrappers

The solver requires mesh-partitioning libraries:

- METIS
- ParMETIS, when using distributed mesh partitioning

## CPU Build

Open a terminal in the source directory and select the appropriate compiler, library paths, and optimization flags in the Makefile.

For a clean build:

```bash
make -f Makefile clean all
```

For an incremental rebuild:

```bash
make -f Makefile
```

The executable is:

```bash
ucns3d_p
```

For new development, start with debug flags and runtime checks. For production simulations, rebuild with optimized compiler settings appropriate for the target machine.

## XPU / Accelerated Build

In UCNS3D, **XPU** refers to the accelerator build path enabled by the `xpu` preprocessor flag. It is used for selected kernels that contain OpenMP target offload regions and accelerator-specific data handling. The MPI ranks still run as normal host processes, while supported computational kernels may be offloaded to GPUs or other accelerator devices by the compiler runtime.

This is not a separate solver and it is not a full-code GPU port. It is a compile-time option that activates the accelerated versions of supported routines while the rest of the code continues to use the standard MPI/OpenMP CPU execution path.

Selected kernels can be compiled for OpenMP target acceleration by enabling the `xpu` preprocessor path. On Cray systems, this is typically done with `-Dxpu` together with OpenMP-enabled compiler flags.

Example Cray-style compilation flags:

```bash
ftn -eZ -s real64 -fbackslash -fopenmp -e 0 -e I -O2 -Dxpu
```

The XPU build uses fixed compile-time bounds for several accelerator work arrays. These bounds must match the case being compiled. They are generated from `UCNS3D.DAT` using:

```bash
python3 gpu_max_flags.py UCNS3D.DAT
```

The script emits preprocessor definitions such as:

```bash
-DGPU_MAX_DIM=3 -DGPU_MAX_IORDER=2 -DGPU_MAX_NVAR=5
```

When using the supplied accelerated Makefile, the flags are generated automatically if `UCNS3D.DAT` is present in the build directory:

```bash
make -f view_make/makefiles_script/Makefile COMPILER=cray
```

The configuration file can also be supplied explicitly:

```bash
make -f view_make/makefiles_script/Makefile COMPILER=cray GPU_MAX_CONFIG=/path/to/UCNS3D.DAT
```

Rebuild the executable whenever the runtime case changes the dimensionality, spatial order, number of equations, turbulence setting, passive scalars, or species configuration. CPU builds do not require `UCNS3D.DAT` at compile time.

## Running a Simulation

A run directory normally contains:

- a mesh file named `grid.msh`, or native files named `GRID.cel`, `GRID.vrt`, and `GRID.bnd`
- the runtime input file `UCNS3D.DAT`
- the executable `ucns3d_p`

Set the number of OpenMP threads:

```bash
export OMP_NUM_THREADS=N
```

Use `N=1` for MPI-only execution.

Launch the solver with MPI:

```bash
mpirun -np M ./ucns3d_p
```

Here `M` is the number of MPI processes. At least two MPI processes are normally required.

## Example Workflow

```bash
cd src
make -f Makefile clean all

cd ../run_case
cp ../src/ucns3d_p .
export OMP_NUM_THREADS=1
mpirun -np 32 ./ucns3d_p
```

Machine-specific run scripts for local and HPC systems may be placed under [scripts](/scripts).

## Visualization

UCNS3D output can be post-processed with:

- [ParaView](https://www.paraview.org/)
- [Tecplot](https://www.tecplot.com/)
- [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit)

## Example Cases

Representative test cases are available from:

- [tests1](https://doi.org/10.5281/zenodo.3375432)
- [tests2](https://doi.org/10.5281/zenodo.6538622)

See `TESTS.md` for details.

## License

UCNS3D is distributed under the GNU General Public License v3. See the `LICENSE` file for details.

## Support

Questions, feedback, and suggestions can be sent to:

```text
ucns3d@gmail.com
```
