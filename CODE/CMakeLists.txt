add_executable(ucns3d_p
  declarations.o
  translate.o
  grid_t.o
  flow_operations.o
  memory.o
  io.o
  mpi_p.o
  communications.o
  riemann.o
  flux_p.o
  profile.o
  init_p.o
  parts.o
  bc_p.o
  matrix.o
  grid_p.o
  basis.o
  gradients.o
  prestore.o
  reconstruct_p.o
  source.o
  implicit_fluxes.o
  implicit_time.o
  flux_pv.o
  time_p.o
  der_r.o
  svd.o
  local_pt.o
  main.o
  libparmetis.a
  libmetis.a
  libtecio.a
  ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
  ${MKLROOT}/lib/intel64/libmkl_sequential.a
  ${MKLROOT}/lib/intel64/libmkl_core.a
)

add_library(ParMetis STATIC IMPORTED)
add_library(Metis STATIC IMPORTED)
add_library(TecIo STATIC IMPORTED)

add_library(Mkl_Intel_Lp64 SHARED IMPORTED)
add_library(Mkl_Sequential SHARED IMPORTED)
add_library(Mkl_Core SHARED IMPORTED)

target_link_libraries(ucns3d_p 
  PRIVATE
    MPI::MPI_Fortran ParMetis Metis TecIo -Wl,--start-group Mkl_Intel_Lp64 Mkl_Sequential Mkl_Core -Wl,--end-group pthread m dl stdc++
)

set_target_properties(Metis PROPERTIES
  IMPORTED_LOCATION ${UCNS3D_ROOT}/libmetis.a
)

set_target_properties(ParMetis PROPERTIES
  IMPORTED_LOCATION ${UCNS3D_ROOT}/libparmetis.a
)

set_target_properties(TecIo PROPERTIES
  IMPORTED_LOCATION ${UCNS3D_ROOT}/libtecio.a
)

set_target_properties(Mkl_Intel_Lp64 PROPERTIES
  IMPORTED_LOCATION ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
)

set_target_properties(Mkl_Sequential PROPERTIES
  IMPORTED_LOCATION ${MKLROOT}/lib/intel64/libmkl_sequential.a
)

set_target_properties(Mkl_Core PROPERTIES
  IMPORTED_LOCATION ${MKLROOT}/lib/intel64/libmkl_core.a
)

set_property(TARGET ucns3d_p
  PROPERTY
    LINKER_LANGUAGE Fortran
)
