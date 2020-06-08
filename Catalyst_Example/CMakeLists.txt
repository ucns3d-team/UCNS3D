cmake_minimum_required(VERSION 3.10)
project(CatalystExamples C CXX)

include(CMakeDependentOption)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../../VTK/CMake")
cmake_dependent_option(USE_CATALYST "Link the simulator with Catalyst" ON
  "NOT WIN32" OFF)
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

option(BUILD_TESTING "Build Testing" OFF)
# Setup testing.
if (BUILD_TESTING)
  include(CTest)
endif()

option(BUILD_FORTRAN_EXAMPLES "Build Fortran Catalyst Examples" OFF)
if (BUILD_FORTRAN_EXAMPLES)
  list(APPEND CMAKE_MODULE_PATH
    "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake")
  if (CMAKE_Fortran_COMPILER)
    enable_language(Fortran OPTIONAL)
  endif ()
  if (CMAKE_Fortran_COMPILER_WORKS)
    find_package(MPI REQUIRED COMPONENTS Fortran)

    add_subdirectory(Fortran90FullExample)
    
  endif ()
endif ()
