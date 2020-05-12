add_executable(ucns3d_p
  FECxxAdaptor.cxx
  declarations.f90
  translate.f90
  grid_t.f90
  flow_operations.f90
  memory.f90
  io.f90
  mpi_p.f90
  communications.f90
  riemann.f90
  flux_p.f90
  profile.f90
  init_p.f90
  parts.f90
  bc_p.f90
  matrix.f90
  grid_p.f90
  basis.f90
  gradients.f90
  prestore.f90
  reconstruct_p.f90
  source.f90
  implicit_fluxes.f90
  implicit_time.f90
  flux_pv.f90
  time_p.f90
  der_r.f90
  svd.f90
  local_pt.f90
  FEFortranAdaptor.f90
  main.f90
  libparmetis.a
  libmetis.a
  libtecio.a
)

add_library(ParMetis STATIC IMPORTED)
add_library(Metis STATIC IMPORTED)
add_library(TecIo STATIC IMPORTED)

add_library(Mkl_Intel_Lp64 SHARED IMPORTED)
add_library(Mkl_Sequential SHARED IMPORTED)
add_library(Mkl_Core SHARED IMPORTED)

target_link_libraries(ucns3d_p
  PRIVATE
    MPI::MPI_Fortran ParMetis Metis TecIo -Wl,--start-group Mkl_Intel_Lp64 Mkl_Sequential Mkl_Core -Wl,--end-group
    ParaView::PythonCatalyst
    VTK::CommonDataModel
    VTK::IOLegacy
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
