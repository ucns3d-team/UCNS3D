! Two independent meshes share one implementation without derived types.
! Volume and surface exchanges remain separate and can be interleaved safely.
module vtk_ghost_output
   include 'vtk_ghost_core.inc'
end module vtk_ghost_output

module vtk_surface_ghost_output
   include 'vtk_ghost_core.inc'
end module vtk_surface_ghost_output
