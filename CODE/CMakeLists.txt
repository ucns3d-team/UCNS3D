set(SRCS
    basis.f90
    bc_p.f90
    communications.f90
    declarations.f90
    der_r.f90
    flow_operations.f90
    flux_p.f90
    flux_pv.f90
    gradients.f90
    grid_p.f90
    grid_t.f90
    implicit_fluxes.f90
    implicit_time.f90
    init_p.f90
    io.f90
    local_pt.f90
    main.f90
    matrix.f90
    memory.f90
    mood.f90
    mpi_p.f90
    parameters.f90
    parts.f90
    prestore.f90
    profile.f90
    reconstruct_p.f90
    riemann.f90
    source.f90
    svd.f90
    time_p.f90
    translate.f90
)

find_library(BLIS_LIBRARY NAMES libblis.a)
find_library(ParMetis_LIBRARY NAMES libparmetis.a)
find_library(Metis_LIBRARY NAMES libmetis.so)

add_executable(ucns3d_p ${SRCS} )

target_link_libraries(ucns3d_p ${BLIS_LIBRARY} ${ParMetis_LIBRARY} ${Metis_LIBRARY}
    ${CMAKE_CURRENT_SOURCE_DIR}/libtecio.so
)
