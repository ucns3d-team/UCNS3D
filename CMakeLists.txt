cmake_minimum_required(VERSION 3.12)
project(UCNS3D C CXX)

include(CMakeDependentOption)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/CMake")

#set(MKLROOT /home/jason/intel/parallel_studio_xe_2019.5.075/compilers_and_libraries_2019/linux/mkl)
#set(MKLROOT "<Your_MKL_Dir>")

set(UCNS3D_ROOT /home/jason/Desktop/UCNS3D/CODE)
#set(UCNS3D_ROOT "<Your_UCNS3D_Dir>")

set(PARAVIEW_BUILD /home/jason/Desktop/pv_gcc_build) # GCC
#set(PARAVIEW_BUILD /home/jason/Desktop/pv) # Intel
#set(PARAVIEW_BUILD "<Your_Paraview_Build_Dir>")

if (NOT EXISTS ${UCNS3D_ROOT})
  message(FATAL_ERROR "Please specify a valid source directory for UCNS3D.")
  return ()
endif ()

if (NOT EXISTS ${PARAVIEW_BUILD})
  message(FATAL_ERROR "Please specify a valid build directory for Paraview.")
  return ()
endif ()

set(CMAKE_PREFIX_PATH ${PARAVIEW_BUILD} CACHE STRING "")

# Compile with GCC if MKLROOT is empty
if (NOT DEFINED MKLROOT) 

  message(STATUS "..........................................\n"
	"-- ..........................................\n"
	"-- ............Compiling with GCC............\n"
	"-- ..........................................\n"
	"-- ..........................................")

  set(BLASFLAGS "-I /usr/local/include")
  set(BLASLIBS "-L/usr/local/lib")
  set(CMAKE_Fortran_FLAGS_DEBUG "-ffree-line-length-none -fdefault-real-8 -fopenmp -Wall -Wno-unused-parameter  -fbounds-check -finit-local-zero  -fdump-parse-tree -fdump-core -fbacktrace -fdefault-double-8 -fbackslash  -O0 -fcray-pointer -Wno-lto-type-mismatch ${BLASFLAGS}" CACHE STRING "")

  set(CMAKE_Fortran_FLAGS_RELEASE "-fdefault-real-8 -fdefault-double-8 -fbackslash -fopenmp -ffree-line-length-none -fcray-pointer -Ofast  -Wno-lto-type-mismatch ${BLASFLAGS}" CACHE STRING "")
  
# Compile with Intel compilers if MKLROOT is not empty
else ()

  message(STATUS "..........................................\n"
	"-- ..........................................\n"
	"-- ......Compiling with Intel compilers......\n"
	"-- ..........................................\n"
	"-- ..........................................")

  if (NOT EXISTS ${MKLROOT})
    message(FATAL_ERROR "Please specify a valid directory for the Intel Math Kernel Library.")
    return ()
  endif ()

  set(BLASFLAGS "-I ${MKLROOT}/include")
  set(BLASLIBS "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")

  set(CMAKE_Fortran_FLAGS_DEBUG "-i4 -r8 -O0 -g -debug all -xHost -traceback -qopenmp -qopenmp-link=static ${BLASFLAGS} -check bounds -fp-model strict -debug inline-debug-info -zero -Wl,--no-as-needed -Wl,--unresolved-symbols=ignore-in-shared-libs" CACHE STRING "")

  set(CMAKE_Fortran_FLAGS_RELEASE "-i4 -r8 -O3 -xHost -fp-model precise -zero -qopenmp  -qopenmp-link=static ${BLASFLAGS}" CACHE STRING "")

endif ()

# Locate Paraview package
find_package(ParaView REQUIRED)

# Check if Catalyst is enabled for Paraview
if (NOT TARGET ParaView::PythonCatalyst)
  message(STATUS
  "Skipping example: ${CMAKE_PROJECT_NAME} requires ParaView to be built "
  "with Catalyst and Python support enabled. Please rebuild ParaView (or "
  "point to a different build of ParaView) with PARAVIEW_USE_PYTHON set "
  "to TRUE")
  return ()
endif ()

# Check if MPI is enabled for Paraview
if (NOT PARAVIEW_USE_MPI)
  message(STATUS
  "Skipping example: ${CMAKE_PROJECT_NAME} requires ParaView to be built "
  "with MPI support enabled. Please rebuild ParaView (or point to a "
  "different build of ParaView) with PARAVIEW_USE_MPI set to TRUE")
  return ()
endif ()

# Check if Fortran compiler exists
if (CMAKE_Fortran_COMPILER)
  enable_language(Fortran OPTIONAL)
endif ()

if (CMAKE_Fortran_COMPILER_WORKS)
  find_package(MPI REQUIRED COMPONENTS Fortran)

  # Check if Fortran MPI compiler is found
  if (NOT MPI_FOUND)
    message(STATUS "Skipping example: ${CMAKE_PROJECT_NAME} requires MPI support, but none was found.")
    return ()
  endif ()

  add_subdirectory(CODE)
endif ()

