module partition
   use mesh_mpi_io
   USE MPIINFO
   use DECLARATION

   use iso_c_binding, only: c_int, c_int64_t, c_float, c_ptr, c_loc
   implicit none

#if defined(UCNS3D_METIS_INDEX_BITS) && UCNS3D_METIS_INDEX_BITS != 64
#error Only 64-bit ParMETIS indices are supported
#endif
   integer, parameter :: partition_idx_kind = c_int64_t
   ! METIS/ParMETIS real_t remains float (REALTYPEWIDTH=32), independently
   ! of the solver's default real precision. MPI handles/counts remain native.
   private :: partition_library_checked
   logical, save :: partition_library_checked = .false.

   interface
      function metis_options_raw(options) bind(C, name="METIS_SetDefaultOptions") result(status)
         import c_ptr, c_int
         type(c_ptr), value :: options
         integer(c_int) :: status
      end function
      ! Use ParMETIS's Fortran bridge: it converts MPI_Fint with MPI_Comm_f2c.
      function parmetis_mesh(elmdist, eptr, eind, vwgt, wgtflag, numflag, ncon, ncommon, &
                            nparts, tpwgts, ubvec, options, edgecut, part, comm) &
                            bind(C, name="parmetis_v3_partmeshkway_") result(status)
         import partition_idx_kind, c_float, c_int
         integer(partition_idx_kind) :: elmdist(*), eptr(*), eind(*), vwgt(*)
         integer(partition_idx_kind) :: wgtflag, numflag, ncon, ncommon, nparts
         integer(partition_idx_kind) :: options(*), edgecut, part(*)
         real(c_float) :: tpwgts(*), ubvec(*)
         integer(c_int) :: comm, status
      end function
   end interface

contains

   subroutine partition_fail(message)
      implicit none
      character(*), intent(in) :: message
      integer :: ierr
      write(*, '(A,I0,2A)') 'Partitioner error on rank ', n, ': ', trim(message)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      error stop 'partitioner failure'
   end subroutine

   subroutine partition_status(status, routine)
      implicit none
      integer(c_int), intent(in) :: status
      character(*), intent(in) :: routine
      character(24) :: code
      if (status == 1_c_int) return ! METIS_OK is a C int, even with 64-bit idx_t.
      write(code, '(I0)') status
      call partition_fail(trim(routine)//' returned '//trim(code))
   end subroutine

   subroutine partition_check_library()
      implicit none
      ! METIS 5 initializes exactly 40 idx_t entries to -1. A 320-byte,
      ! 8-byte-aligned buffer is safe for either supported integer ABI.
      ! Detect the actual loaded METIS width before passing typed arrays.
      ! This does not validate a separately mismatched ParMETIS or real_t ABI;
      ! both libraries must be built together with REALTYPEWIDTH=32.
      integer(c_int64_t), target :: probe(40)
      integer(c_int) :: status
      integer :: library_bits
      if (partition_library_checked) return
      if (storage_size(0) /= 32) call partition_fail('use 32-bit default solver/MPI integers; do not use -i8')
      probe = 0_c_int64_t
      status = metis_options_raw(c_loc(probe))
      call partition_status(status, 'METIS_SetDefaultOptions ABI probe')
      if (all(probe == -1_c_int64_t)) then
         library_bits = 64
      else if (all(probe(1:20) == -1_c_int64_t) .and. all(probe(21:40) == 0_c_int64_t)) then
         library_bits = 32
      else
         call partition_fail('unrecognized METIS options layout; require a compatible METIS 5 library')
      end if
      if (library_bits /= storage_size(0_partition_idx_kind)) then
         write(*, '(A,I0,A,I0)') 'METIS integer ABI mismatch: solver interface=', &
            storage_size(0_partition_idx_kind), ', loaded METIS=', library_bits
         call partition_fail('64-bit METIS/ParMETIS libraries are required; check library paths')
      end if
      partition_library_checked = .true.
   end subroutine

   subroutine partition_check_owners(parts, base)
      implicit none
      integer(partition_idx_kind), intent(in) :: parts(:), base
      if (any(parts < base) .or. any(parts >= int(isize, partition_idx_kind) + base)) &
         call partition_fail('partitioner returned an invalid MPI owner')
   end subroutine


   Subroutine Partitioner6(n, IMAXE, imaxn, xmpie, ieshape)
      !> @brief
!> Partition the mesh using the 64-bit ParMETIS option-6 interface.
      use ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      Implicit None
!!!!!!!!!!!!!!!
      integer, ALLOCATABLE, DIMENSION(:), intent(in) :: IESHAPE
      real :: average, average2, tsize
      integer :: maxi, maxi2
      Integer :: i, j, elementid, node1, node2, node3, node4, node5, node6, node7, node8, counternodes2, k
      Integer :: vweight
      Integer, allocatable, dimension(:, :, :) :: Compweight, Commweight
      INTEGER, INTENT(IN) :: N, IMAXE, imaxn
      INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: XMPIE
      integer, allocatable, dimension(:) :: testar, pweight, pweight2
      integer(c_int) :: partition_result
      integer(partition_idx_kind) :: counternodes
      integer(partition_idx_kind) :: imaxee, imaxnn, ncommon1, isizee, objval, NUMFLAG, ncon, NCOMMONNODES, edgecut, WGTFLAG
      integer(partition_idx_kind), allocatable, dimension(:) :: AllnodesPTR, ALLNODES, parts_el
      integer, allocatable :: xmpiee(:), mpi_displs(:), parts_local(:)
      integer(partition_idx_kind), allocatable :: vwgt(:)
      integer(partition_idx_kind) :: local_capacity
      integer(c_int64_t) :: local_count_wide
      integer(partition_idx_kind), allocatable :: elmdist(:), eptr(:), eind(:), options(:)
      real(c_float), allocatable, dimension(:) :: tpwgts
      real(c_float), allocatable, dimension(:) :: UBVEC

      integer :: mesh_nodes(8), mesh_nn, mesh_shape, mesh_faces, mesh_i
      call partition_check_library()
      if (isize < 1 .or. imaxe < isize) &
         call partition_fail('ParMETIS requires at least one input cell per MPI rank')
      WGTfLAG = 0
      imaxee = imaxe
      imaxnn = imaxn
      isizee = isize
      tsize = isizee
      allocate (xmpiee(1:isize))

      xmpiee = imaxe / isize
      xmpiee(isize) = imaxe - (imaxe / isize) * (isize - 1)

      allocate (elmdist(isize + 1))

      elmdist(1) = 1
      do i = 2, isize + 1
         elmdist(i) = elmdist(i - 1) + xmpiee(i - 1)
!
      end do
!

      if (dimensiona .eq. 3) then
         local_count_wide = 8_c_int64_t * int(xmpiee(n + 1), c_int64_t)
      else
         local_count_wide = 4_c_int64_t * int(xmpiee(n + 1), c_int64_t)
      end if

      if (local_count_wide >= int(huge(0_partition_idx_kind), c_int64_t)) &
         call partition_fail('local connectivity exceeds 64-bit partition index range; use more MPI ranks')
      local_capacity = int(local_count_wide, partition_idx_kind)
      Allocate (ALLNODES(1:local_capacity))
      Allocate (ePTR(1:xmpiee(n + 1) + 1))

      counternodes = 0
      counternodes2 = 1
      ePTR(counternodes2) = 1

      if (mesh_enabled) then
         if (size(mesh_cells, 2) /= xmpiee(n + 1)) call partition_fail('distributed mesh block/count mismatch')
         do mesh_i = 1, size(mesh_cells, 2)
            call mesh_cell_nodes(mesh_cells(:, mesh_i), mesh_nodes, mesh_nn, mesh_shape, mesh_faces)
            ALLNODES(counternodes + 1:counternodes + mesh_nn) = mesh_nodes(1:mesh_nn)
            counternodes = counternodes + mesh_nn
            counternodes2 = counternodes2 + 1
            ePTR(counternodes2) = counternodes + 1
         end do
      else
      IF (BINIO .EQ. 0) THEN

         if (dimensiona .eq. 3) then

            OPEN (1112, FILE='GRID.cel', FORM='FORMATTED', STATUS='OLD', ACTION='READ')

            Do i = 1, imaxe
               if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
                  counternodes2 = counternodes2 + 1
                  Read (1112, *) elementid, node1, node2, node3, node4, node5, node6, node7, node8
                  IF ((node3 .EQ. node4) .AND. (node5 .NE. node6) .AND. (node7 .EQ. node8)) THEN ! prism
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node6
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node7

                  ! HEXA
                  Else IF ((node5 .NE. node6) .AND. (node6 .NE. node7) .AND. (node7 .NE. node8) .AND. (node3 .NE. node4)) THEN
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node4
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node6
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node7
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node8

                  ! Tetra
                  Else IF ((node3 .EQ. node4) .AND. (node5 .EQ. node6) .AND. (node6 .EQ. node7) .AND. (node7 .EQ. node8)) THEN
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5

                  Else

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node4
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5
                  End if
                  ePTR(counternodes2) = counternodes + 1
               else
                  Read (1112, *)
               end if
            End do

            close (1112)

         else
            OPEN (1112, FILE='GRID.cel', FORM='FORMATTED', STATUS='OLD', ACTION='READ')

            Do i = 1, imaxe
               if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
                  counternodes2 = counternodes2 + 1
                  Read (1112, *) elementid, node1, node2, node3, node4
                  IF ((node3 .EQ. node4)) THEN ! triangle
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3

                  Else
                     !qudrilateral

                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node4

                  End if
                  ePTR(counternodes2) = counternodes + 1
               else
                  Read (1112, *)
               end if
            End do

            close (1112)

         end if

      ELSE
         if (dimensiona .eq. 3) then

            OPEN (1112, FILE='GRID.cel', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')

            Do i = 1, imaxe
               if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
                  counternodes2 = counternodes2 + 1
                  Read (1112) elementid, node1, node2, node3, node4, node5, node6, node7, node8
                  IF ((node3 .EQ. node4) .AND. (node5 .NE. node6) .AND. (node7 .EQ. node8)) THEN ! prism
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node6
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node7

                  ! HEXA
                  Else IF ((node5 .NE. node6) .AND. (node6 .NE. node7) .AND. (node7 .NE. node8) .AND. (node3 .NE. node4)) THEN
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node4
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node6
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node7
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node8

                  ! Tetra
                  Else IF ((node3 .EQ. node4) .AND. (node5 .EQ. node6) .AND. (node6 .EQ. node7) .AND. (node7 .EQ. node8)) THEN
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5

                  Else

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node4
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node5
                  End if
                  ePTR(counternodes2) = counternodes + 1
               else
                  Read (1112)
               end if
            End do

            close (1112)

         else
            OPEN (1112, FILE='GRID.cel', FORM='UNFORMATTED', STATUS='OLD', ACTION='READ')

            Do i = 1, imaxe
               if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
                  counternodes2 = counternodes2 + 1
                  Read (1112) elementid, node1, node2, node3, node4
                  IF ((node3 .EQ. node4)) THEN ! triangle
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3

                  Else
                     !qudrilateral

                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node1
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node2
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node3
                     Counternodes = counternodes + 1
                     ALLNODES(counternodes) = node4

                  End if
                  ePTR(counternodes2) = counternodes + 1
               else
                  Read (1112)
               end if
            End do

            close (1112)

         end if

      END IF

      end if

      allocate (eind(1:counternodes))
      eind(1:counternodes) = Allnodes(1:counternodes)
      deAllocate (ALLNODES)

      allocate (vwgt(xmpiee(n + 1)))
      vwgt = 0; wgtflag = 0; numflag = 1; ncon = 1;
      if (dimensiona .eq. 3) then

         NCOMMONNODES = 3
      else
         NCOMMONNODES = 2

      end if

      allocate (tpwgts(ncon * isize))
      tpwgts = 1./float(isize)
      allocate (ubvec(ncon))
      ubvec(1) = 1.050000

      allocate (options(3))
      options = [0, 0, 0]

      allocate (parts_el(xmpiee(n + 1))); parts_el = 0

      if (n .eq. 0) then
         write (*, *) '------------------------------------------------------------------------'
         write (*, *) '                       ParMETIS Initiated                               '
         write (*, *) '------------------------------------------------------------------------'
      end if

      CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

      if (isize == 1) then
         parts_el = 1
      else
         if (n == 0) write(*, '(A,I0,A)') 'Partitioner interface: ', &
            storage_size(0_partition_idx_kind), '-bit indices; 32-bit real_t'
         partition_result = parmetis_mesh(elmdist, eptr, eind, vwgt, WGTFLAG, numflag, NCON, NCOMMONNODES, &
            isizee, TPWGTS, UBVEC, OPTIONS, EDGECUT, PARTs_el, MPI_COMM_WORLD)
         call partition_status(partition_result, 'ParMETIS_V3_PartMeshKway')
      end if
      call partition_check_owners(parts_el, 1_partition_idx_kind)
      allocate(parts_local(size(parts_el)))
      parts_local = int(parts_el, kind(0))
      deallocate(parts_el)

      if (n .eq. 0) then
         write (*, *) '------------------------------------------------------------------------'
         write (*, *) '                       ParMETIS Operation Completed                     '
         write (*, *) '------------------------------------------------------------------------'
      end if

      deAllocate (ePTR)
      deallocate (eind)
      deallocate (vwgt)
      deallocate (tpwgts)

      deallocate (UBVEC)
      deallocate (options)

      if (mesh_enabled .and. iperiodicity /= 1) then
         if (n == 0) write(*, '(A,I0)') &
            'MPI mesh storage: distributed cells/local vertices; periodicity flag ', iperiodicity
         call mesh_store_partition(parts_local - 1)
      else
      allocate(mpi_displs(isize))
      ! Global cell IDs still fit the solver's native integer. The library's
      ! 64-bit offsets must not be passed as MPI_INTEGER counts/displacements.
      mpi_displs = int(elmdist(1:isize) - 1_partition_idx_kind, kind(0))
      call MPI_Allgatherv(parts_local, xmpiee(n + 1), MPI_INTEGER, xmpie, xmpiee, mpi_displs, &
         MPI_INTEGER, MPI_COMM_WORLD, ierror)
      if (ierror /= MPI_SUCCESS) call partition_fail('MPI_Allgatherv of partition owners failed')
      deallocate(mpi_displs)

      end if

      deAllocate (parts_local)
      deallocate (xmpiee)
      deallocate (elmdist)
      XMPIE = XMPIE - 1

   end subroutine Partitioner6

end module partition
