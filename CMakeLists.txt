cmake_minimum_required(VERSION 3.12)
project(UCNS3D C CXX)

set(MKLROOT /home/jason/intel/parallel_studio_xe_2019.5.075/compilers_and_libraries_2019/linux/mkl)
set(BLASFLAGS "-I ${MKLROOT}/include")
set(BLASLIBS "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")

set(LIBS "${BLASLIBS} -lstdc")

set(CMAKE_Fortran_FLAGS "-i4 -r8 -O0  -g  -debug all -xHost -traceback -qopenmp -qopenmp-link=static ${BLASFLAGS} -check bounds -fp-model strict -debug inline-debug-info -zero" CACHE STRING "")

set(UCNS3D_ROOT /home/jason/Desktop/UCNS3D/CODE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/CMake")

option(BUILD_FORTRAN_EXAMPLES "Build Fortran Catalyst Examples" OFF)
if (BUILD_FORTRAN_EXAMPLES)
  list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake")
  if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
    set(CMAKE_Fortran_COMPILER_WORKS ON)
  endif ()
  if (CMAKE_Fortran_COMPILER_WORKS)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    add_subdirectory(CODE)
  endif ()
endif ()
