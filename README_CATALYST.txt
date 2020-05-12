Porting to Catalyst (Work in Progress)

1. CMake
- Edited Makefile.common to use only output files, may abolish Makefile later in favor of
  CMake build generator.
- CMake currently uses Intel Fortran as well as MKL libraries, work needs to be done to
  include gfortran.
- Steps to compile using CMake

(NOTE: Steps below only serve as author's personal reference, more work needs to be done
 to generalize it for multiple users.)
-----------------------------------------------------------------------------------------
// Load Intel compiler modules
cd /home/jason/intel/parallel_studio_xe_2019.5.075/bin
source psxevars.sh

// Store Intel libraries as system variables
export intelMPI=/home/jason/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64
export intelFort=/home/jason/intel/compilers_and_libraries_2019.5.281/linux/bin/intel64

// Navigate to directory for UCNS3D
cd $HOME/Desktop/UCNS3D

// Create .o files using Makefile
cd CODE
make -f Makefile clean all

// Compile using CMake in a separate build folder
cd ../
mkdir build
cd build
cmake -DBUILD_FORTRAN_EXAMPLES:BOOL=ON \
	-DCMAKE_BUILD_TYPE:STRING=Debug \
	-DCMAKE_Fortran_COMPILER="$intelMPI/bin/mpiifort" \
	-DMPI_Fortran_INCLUDE_PATH="$intelMPI/include" \
	--debug-output ../

make

-----------------------------------------------------------------------------------------
2. Catalyst integration (WIP)

// Compile UCNS3D

cmake -DCMAKE_PREFIX_PATH="/home/jason/Desktop/paraview_build" \
	-DUSE_CATALYST=ON -DBUILD_FORTRAN_EXAMPLES:BOOL=ON \
	-DCMAKE_BUILD_TYPE:STRING=Debug \
	-DCMAKE_Fortran_COMPILER="$intelMPI/bin/mpiifort" \
	-DMPI_Fortran_INCLUDE_PATH="$intelMPI/include" \
	--debug-output ../

make

