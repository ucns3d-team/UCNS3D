cmake_minimum_required(VERSION 3.12)
project(UCNS3D C CXX)

include(CMakeDependentOption)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/CMake")
cmake_dependent_option(USE_CATALYST "Link the simulator with Catalyst" ON
  "NOT WIN32" OFF)

set(MKLROOT /home/jason/intel/parallel_studio_xe_2019.5.075/compilers_and_libraries_2019/linux/mkl)
set(BLASFLAGS "-I ${MKLROOT}/include")
set(BLASLIBS "-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl")

set(LIBS "${BLASLIBS} -lstdc")

set(CMAKE_Fortran_FLAGS_DEBUG "-i4 -r8 -O0  -g  -debug all -xHost -traceback -qopenmp -qopenmp-link=static ${BLASFLAGS} -check bounds -fp-model strict -debug inline-debug-info -zero -Wl,--no-as-needed -Wl,--unresolved-symbols=ignore-in-shared-libs" CACHE STRING "")

set(UCNS3D_ROOT /home/jason/Desktop/UCNS3D/CODE)

set(PARAVIEW_BUILD /home/jason/Desktop/paraview_build)

if (USE_CATALYST)
  find_package(ParaView REQUIRED)

  if (NOT TARGET ParaView::PythonCatalyst)
    message(STATUS
      "Skipping example: ${CMAKE_PROJECT_NAME} requires ParaView to be built "
      "with Catalyst and Python support enabled. Please rebuild ParaView (or "
      "point to a different build of ParaView) with PARAVIEW_USE_PYTHON set "
      "to TRUE")
    return ()
  endif()

  if (NOT PARAVIEW_USE_MPI)
    message(STATUS
      "Skipping example: ${CMAKE_PROJECT_NAME} requires ParaView to be built "
      "with MPI support enabled. Please rebuild ParaView (or point to a "
      "different build of ParaView) with PARAVIEW_USE_MPI set to TRUE")
    return ()
  endif ()

  # FIXME: This should really be fixed to instead be done per-target.
  add_definitions(-DUSE_CATALYST)
else ()
  find_package(MPI COMPONENTS C)
  if (NOT MPI_FOUND)
    message(STATUS
      "Skipping example: ${CMAKE_PROJECT_NAME} requires MPI support, but none "
      "was found.")
    return ()
  endif ()
endif ()

option(BUILD_FORTRAN_EXAMPLES "Build Fortran Catalyst Examples" OFF)
if (BUILD_FORTRAN_EXAMPLES)
  list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake")
  if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
  endif ()
  if (CMAKE_Fortran_COMPILER_WORKS)
    find_package(MPI REQUIRED COMPONENTS Fortran)
    add_subdirectory(CODE)
  endif ()
endif ()
