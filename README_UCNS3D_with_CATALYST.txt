BUILD STEPS FOR PARAVIEW AND CATALYST
----------------------------------------------------------------------------------------------------
// Get CMake versions 3.12 - 3.16 for Master Build

// 1. Paraview build compilation using GCC 7.5.0 with Ubuntu 18.04 (RECOMMENDED)

cmake -DPARAVIEW_USE_PYTHON=ON \
	-DPARAVIEW_USE_MPI=ON \
	-DPARAVIEW_USE_FORTRAN=ON \
	-DPARAVIEW_BUILD_CATALYST_ADAPTORS=ON \
	-DCMAKE_Fortran_COMPILER="/usr/bin/mpifort;/usr/bin/gfortran" \
	-DCMAKE_C_COMPILER="/usr/bin/cc" \
	-DMPI_C_COMPILER="/usr/bin/mpicc" ../paraview

make -j 4

// 2. Compile UCNS3D with Catalyst

cmake -DCMAKE_PREFIX_PATH="<Your_Paraview_Build_Dir>" \
	-DCMAKE_BUILD_TYPE:STRING=Release \
	-DCMAKE_Fortran_COMPILER="/usr/bin/mpifort" ../

cmake -DCMAKE_PREFIX_PATH="/home/jason/Desktop/pv_gcc_build" \
	-DCMAKE_BUILD_TYPE:STRING=Release \
	-DCMAKE_Fortran_COMPILER="/usr/bin/mpifort" ../

// 3. For testing of UCNS3D

cd CODE
mkdir RUN && cd RUN
cp ../../../Catalyst_Example/UCNS3D_obj/* .
cp ../ucns3d_p .

export OMP_NUM_THREADS=1
mpirun -np 4 ./ucns3d_p

----------------------------------------------------------------
// Alternative build version with Intel compilers

// 1. Paraview build compilation using Intel XE 2019 with Ubuntu 18.04

cd <Your_Intel_Compiler_Install_Dir>/parallel_studio_xe_2019.5.075/bin
source psxevars.sh

cd <Your_Paraview_Build_Dir>
export intelMPI=<Your_Intel_Compiler_Install_Dir>/compilers_and_libraries_2019.5.281/linux/mpi/intel64
export intelFort=<Your_Intel_Compiler_Install_Dir>/compilers_and_libraries_2019.5.281/linux/bin/intel64

cmake -DPARAVIEW_USE_PYTHON=ON \
	-DPARAVIEW_USE_MPI=ON \
	-DPARAVIEW_USE_FORTRAN=ON \
	-DPARAVIEW_BUILD_CATALYST_ADAPTORS=ON \
	-DCMAKE_Fortran_COMPILER="$intelFort/ifort" \
	-DMPI_Fortran_INCLUDE_PATH="$intelMPI/include" \
	-DMPI_Fortran_LIBRARIES="$intelMPI/lib/libmpifort.so" ../paraview

make -j 4

// 2. Compile UCNS3D with Catalyst

cmake -DCMAKE_PREFIX_PATH="<Your_Paraview_Build_Dir>" \
	-DCMAKE_BUILD_TYPE:STRING=Debug \
	-DCMAKE_Fortran_COMPILER="$intelMPI/bin/mpiifort" \
	-DMPI_Fortran_INCLUDE_PATH="$intelMPI/include" ../

cmake -DCMAKE_PREFIX_PATH="/home/jason/Desktop/pv" \
	-DCMAKE_BUILD_TYPE:STRING=Release \
	-DCMAKE_Fortran_COMPILER="$intelMPI/bin/mpiifort" ../
