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
- Case-specific compile-time bounds for fixed work arrays
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

UCNS3D can be used through Docker or compiled manually on the target machine. The canonical build entry point is [src/Makefile](src/Makefile), which provides compiler profiles selected with `COMPILER=...`.

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

The current [Dockerfile](Dockerfile) builds the GNU CPU/OpenMP executable with:

```bash
make -f Makefile COMPILER=gnu all
```

It also installs `python3`, which is used during compilation to derive case-specific fixed array bounds when `UCNS3D.DAT` is present.

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

The repository build expects Tecplot, METIS, and ParMETIS libraries under [bin/lib](bin/lib) by default. Library paths can be overridden from the make command line.

## CPU Build

Open a terminal in the source directory:

```bash
cd src
```

Select one of the Makefile compiler profiles.

GNU Fortran / MPI, CPU/OpenMP only:

```bash
make -f Makefile COMPILER=gnu clean all
```

Intel Fortran / Intel MPI, CPU/OpenMP only:

```bash
make -f Makefile COMPILER=intel clean all
```

macOS with GNU Fortran uses the same GNU profile. If `libtecio`, `libparmetis`, and `libmetis` are not in `src`, point the build at the directory containing the macOS `.dylib` or `.a` files:

```bash
make -f Makefile COMPILER=gnu MAC_LIB_ROOT=/path/to/macos/libs clean all
```

For an incremental rebuild, omit `clean`:

```bash
make -f Makefile COMPILER=gnu all
```

The executable is:

```bash
ucns3d_p
```

The folders under [bin](bin) contain matching compiler-specific Makefile templates, but [src/Makefile](src/Makefile) is the main build file.

## Compile-Time Bounds

Several fixed work arrays use compile-time `GPU_MAX_*` bounds. Despite the historical name, these bounds are used by shared code paths and are relevant to both CPU and XPU builds.

The Makefile calls [src/gpu_max_flags.py](src/gpu_max_flags.py) automatically. With no extra options, the script checks the current build directory for `UCNS3D.DAT` and companion files such as `REALGAS.DAT`, `MULTISPECIES.DAT`, and `MULTISPECIES_DIFF.DAT`.

The normal workflow is therefore:

```bash
cd src
cp /path/to/case/UCNS3D.DAT .
cp /path/to/case/REALGAS.DAT .     # if used by the case
make -f Makefile COMPILER=gnu clean all
```

If the input deck is outside the build directory, pass it explicitly:

```bash
make -f Makefile COMPILER=gnu GPU_MAX_CONFIG=/path/to/UCNS3D.DAT clean all
```

Rebuild whenever the case changes dimensionality, spatial order, number of equations, turbulence setting, passive scalars, or species configuration. If no input deck is visible at compile time, the helper emits conservative default bounds.

## XPU / Accelerated Build

In UCNS3D, **XPU** refers to the accelerator build path enabled by the `xpu` preprocessor flag. It is used for selected kernels that contain OpenMP target offload regions and accelerator-specific data handling. The MPI ranks still run as normal host processes, while supported computational kernels may be offloaded to GPUs or other accelerator devices by the compiler runtime.

This is not a separate solver and it is not a full-code GPU port. It is a compile-time option that activates the accelerated versions of supported routines while the rest of the code continues to use the standard MPI/OpenMP CPU execution path.

On Cray systems, use the Cray compiler profile:

```bash
make -f Makefile COMPILER=cray clean all
```

The `cray` profile uses `ftn`, OpenMP, real64 defaults, and the `-Dxpu` preprocessor flag. It also uses the same automatic `gpu_max_flags.py` mechanism described above.

For an input deck outside the build directory:

```bash
make -f Makefile COMPILER=cray GPU_MAX_CONFIG=/path/to/UCNS3D.DAT clean all
```

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
make -f Makefile COMPILER=gnu clean all

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
