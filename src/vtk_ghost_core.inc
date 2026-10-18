! Shared implementation; each module instance owns independent flat buffers.
   use mpi
   use iso_fortran_env, only: int8, int32, int64, real64
   implicit none
   private
   public :: vtk_ghost_initialize, vtk_ghost_write, vtk_ghost_release
   integer, parameter :: max_vertices = 8, metadata_width = 3 + 2 * max_vertices
   integer :: output_comm, rank_id, rank_count, owned_count, ghost_count, point_count, field_capacity
   logical :: prepared = .false.
   integer, allocatable :: send_count(:), recv_count(:), send_disp(:), recv_disp(:), send_cell(:)
   integer, allocatable :: send_vertex(:, :), send_nnodes(:), point_local(:)
   integer, allocatable :: point_remote_cell(:), point_remote_vertex(:)
   integer, allocatable :: payload_sc(:), payload_rc(:), payload_sd(:), payload_rd(:)
   integer(int32), allocatable :: connectivity(:), offsets(:), cell_owner(:)
   integer(int64), allocatable :: cell_global(:), point_global(:)
   integer(int8), allocatable :: cell_type(:), cell_ghost(:), point_ghost(:)
   real(real64), allocatable :: points(:, :), send_xyz(:, :), recv_xyz(:, :)
   real(real64), allocatable :: send_values(:), recv_values(:), field(:)
contains
   subroutine fail(message)
      character(*), intent(in) :: message
      integer :: ierr
      write(*, '(A,I0,2A)') 'VTU ghost output rank ', rank_id, ': ', trim(message)
      call MPI_Abort(output_comm, 91, ierr)
      error stop
   end subroutine

   subroutine checked(ierr, message)
      integer, intent(in) :: ierr
      character(*), intent(in) :: message
      if (ierr /= MPI_SUCCESS) call fail(message)
   end subroutine

   integer function count32(value) result(count)
      integer(int64), intent(in) :: value
      if (value < 0 .or. value > int(huge(0), int64)) call fail('MPI/local array count exceeds int32')
      count = int(value)
   end function

   subroutine displacements(counts, disp, total)
      integer, intent(in) :: counts(:)
      integer, intent(out) :: disp(:), total
      integer :: i
      integer(int64) :: wide
      wide = 0
      do i = 1, size(counts)
         disp(i) = count32(wide)
         wide = wide + int(counts(i), int64)
      end do
      total = count32(wide)
   end subroutine

   ! Sort integer pairs lexicographically. Used only when building the map.
   recursive subroutine sort_pairs(a, first, last)
      integer, intent(inout) :: a(:, :)
      integer, intent(in) :: first, last
      integer :: i, j, pivot(2), tmp(2)
      if (first >= last) return
      pivot = a(:, first + (last - first) / 2)
      i = first; j = last
      do
         do while (pair_less(a(:, i), pivot))
            i = i + 1
         end do
         do while (pair_less(pivot, a(:, j)))
            j = j - 1
         end do
         if (i > j) exit
         tmp = a(:, i); a(:, i) = a(:, j); a(:, j) = tmp
         i = i + 1; j = j - 1
         if (i > j) exit
      end do
      if (first < j) call sort_pairs(a, first, j)
      if (i < last) call sort_pairs(a, i, last)
   end subroutine

   logical function pair_less(a, b)
      integer, intent(in) :: a(2), b(2)
      pair_less = a(1) < b(1) .or. (a(1) == b(1) .and. a(2) < b(2))
   end function

   integer function find_id(ids, id) result(slot)
      integer, intent(in) :: ids(:), id
      integer :: low, high, middle
      low = 1; high = size(ids); slot = 0
      do while (low <= high)
         middle = low + (high - low) / 2
         if (ids(middle) == id) then
            slot = middle
            return
         else if (ids(middle) < id) then
            low = middle + 1
         else
            high = middle - 1
         end if
      end do
   end function

   subroutine exchange_records(width, outgoing, incoming, sc, rc, sd, rd)
      integer, intent(in) :: width, outgoing(:, :), sc(:), rc(:), sd(:), rd(:)
      integer, intent(out) :: incoming(:, :)
      integer :: ierr, p
      do p = 1, rank_count
         payload_sc(p) = count32(int(width, int64) * sc(p))
         payload_rc(p) = count32(int(width, int64) * rc(p))
         payload_sd(p) = count32(int(width, int64) * sd(p))
         payload_rd(p) = count32(int(width, int64) * rd(p))
      end do
      call MPI_Alltoallv(outgoing, payload_sc, payload_sd, MPI_INTEGER, incoming, payload_rc, payload_rd, &
         MPI_INTEGER, output_comm, ierr)
      call checked(ierr, 'exchanging topology records')
   end subroutine

   subroutine vtk_ghost_initialize(global_cells, shapes, nnodes, cell_nodes, global_points, xyz, max_fields, comm)
      integer, intent(in) :: global_cells(:), shapes(:), nnodes(:), cell_nodes(:, :), global_points(:)
      real(real64), intent(in) :: xyz(:, :)
      integer, intent(in) :: max_fields, comm
      integer, allocatable :: sc(:), rc(:), sd(:), rd(:), cursor(:), degree(:), ptr(:), incident(:)
      integer, allocatable :: used_ids(:), used_local(:), pairs(:, :), registrations(:, :), received(:, :)
      integer, allocatable :: replies(:, :), shared(:, :), local_owner(:), send_pairs(:, :), vertex_map(:)
      integer, allocatable :: metadata(:, :), remote_metadata(:, :), all_ids(:), vtk_types(:)
      integer :: i, j, k, p, first, last, owner, idx, nused, ns, nr, unique_count, nsend, nconn, ierr
      integer :: nlocal_points, dest, cell, slot, total_cells, raw_count
      integer(int64) :: wide
      if (prepared) call vtk_ghost_release()
      output_comm = comm
      call MPI_Comm_rank(comm, rank_id, ierr); call checked(ierr, 'reading rank')
      call MPI_Comm_size(comm, rank_count, ierr); call checked(ierr, 'reading communicator size')
      if (storage_size(0) /= 32) call fail('32-bit default integers are required for solver mesh IDs')
      owned_count = size(global_cells); nlocal_points = size(global_points)
      if (size(shapes) /= owned_count .or. size(nnodes) /= owned_count) call fail('invalid cell metadata')
      if (size(cell_nodes, 2) /= owned_count .or. any(nnodes > size(cell_nodes, 1))) &
         call fail('invalid cell connectivity dimensions')
      if (size(xyz, 2) /= nlocal_points .or. size(xyz, 1) < 2 .or. size(xyz, 1) > 3) &
         call fail('invalid point coordinate dimensions')
      if (max_fields < 0) call fail('negative output field count')
      field_capacity = max_fields
      allocate(payload_sc(rank_count), payload_rc(rank_count), payload_sd(rank_count), payload_rd(rank_count))
      allocate(sc(rank_count), rc(rank_count), sd(rank_count), rd(rank_count), cursor(rank_count))
      allocate(degree(nlocal_points), ptr(nlocal_points + 1), local_owner(nlocal_points))
      allocate(vtk_types(owned_count))
      degree = 0; local_owner = -1
      do i = 1, owned_count
         select case(shapes(i))
         case(1)
            vtk_types(i) = 12; k = 8
         case(2)
            vtk_types(i) = 10; k = 4
         case(3)
            vtk_types(i) = 14; k = 5
         case(4)
            vtk_types(i) = 13; k = 6
         case(5)
            vtk_types(i) = 9; k = 4
         case(6)
            vtk_types(i) = 5; k = 3
         case(7)
            vtk_types(i) = 3; k = 2
         case default
            call fail('unsupported cell shape')
         end select
         if (nnodes(i) /= k) call fail('cell shape and vertex count disagree')
         do j = 1, k
            idx = cell_nodes(j, i)
            if (idx < 1 .or. idx > nlocal_points) call fail('connectivity is not local point indexing')
            degree(idx) = degree(idx) + 1
         end do
      end do
      ptr(1) = 1
      do i = 1, nlocal_points
         ptr(i + 1) = count32(int(ptr(i), int64) + degree(i))
      end do
      allocate(incident(ptr(nlocal_points + 1) - 1))
      degree = 0
      do i = 1, owned_count
         do j = 1, nnodes(i)
            idx = cell_nodes(j, i)
            incident(ptr(idx) + degree(idx)) = i
            degree(idx) = degree(idx) + 1
         end do
      end do
      nused = count(degree > 0)
      allocate(pairs(2, nused), used_ids(nused), used_local(nused))
      k = 0
      do i = 1, nlocal_points
         if (degree(i) == 0) cycle
         k = k + 1; pairs(:, k) = [global_points(i), i]
      end do
      call sort_pairs(pairs, 1, nused)
      used_ids = pairs(1, :); used_local = pairs(2, :)
      do i = 2, nused
         if (used_ids(i) == used_ids(i - 1)) call fail('duplicate local global-point ID')
      end do
      deallocate(pairs)

      ! Distributed vertex directory: each vertex is handled by one hash rank.
      ! This finds edge/corner neighbours without replicating a global mesh.
      sc = 0
      do i = 1, nused
         p = modulo(used_ids(i), rank_count) + 1
         sc(p) = sc(p) + 1
      end do
      call MPI_Alltoall(sc, 1, MPI_INTEGER, rc, 1, MPI_INTEGER, comm, ierr)
      call checked(ierr, 'counting vertex registrations')
      call displacements(sc, sd, ns); call displacements(rc, rd, nr)
      allocate(registrations(2, ns), received(2, nr)); cursor = sd
      do i = 1, nused
         p = modulo(used_ids(i), rank_count) + 1
         cursor(p) = cursor(p) + 1
         registrations(:, cursor(p)) = [used_ids(i), rank_id]
      end do
      call exchange_records(2, registrations, received, sc, rc, sd, rd)
      deallocate(registrations)
      call sort_pairs(received, 1, nr)
      sc = 0; first = 1
      do while (first <= nr)
         last = first
         do while (last < nr)
            if (received(1, last + 1) /= received(1, first)) exit
            last = last + 1
         end do
         do i = first, last
            p = received(2, i) + 1
            sc(p) = count32(int(sc(p), int64) + last - first + 1)
         end do
         first = last + 1
      end do
      call MPI_Alltoall(sc, 1, MPI_INTEGER, rc, 1, MPI_INTEGER, comm, ierr)
      call checked(ierr, 'counting vertex-sharing replies')
      call displacements(sc, sd, ns); call displacements(rc, rd, unique_count)
      allocate(replies(3, ns), shared(3, unique_count)); cursor = sd
      first = 1
      do while (first <= nr)
         last = first
         do while (last < nr)
            if (received(1, last + 1) /= received(1, first)) exit
            last = last + 1
         end do
         owner = received(2, first) ! lexicographically sorted: minimum owning rank
         do i = first, last
            p = received(2, i) + 1
            do j = first, last
               cursor(p) = cursor(p) + 1
               replies(:, cursor(p)) = [received(1, first), received(2, j), owner]
            end do
         end do
         first = last + 1
      end do
      call exchange_records(3, replies, shared, sc, rc, sd, rd)
      deallocate(received, replies)
      wide = 0
      do i = 1, size(shared, 2)
         slot = find_id(used_ids, shared(1, i))
         if (slot == 0) call fail('vertex directory returned an unknown point')
         idx = used_local(slot)
         local_owner(idx) = shared(3, i)
         if (shared(2, i) /= rank_id) wide = wide + degree(idx)
      end do
      raw_count = count32(wide)
      allocate(send_pairs(2, raw_count)); k = 0
      do i = 1, size(shared, 2)
         dest = shared(2, i)
         if (dest == rank_id) cycle
         idx = used_local(find_id(used_ids, shared(1, i)))
         do j = ptr(idx), ptr(idx + 1) - 1
            k = k + 1; send_pairs(:, k) = [dest, incident(j)]
         end do
      end do
      deallocate(shared, incident, ptr, degree)
      call sort_pairs(send_pairs, 1, raw_count)
      nsend = 0
      do i = 1, raw_count
         if (nsend > 0) then
            if (all(send_pairs(:, i) == send_pairs(:, nsend))) cycle
         end if
         nsend = nsend + 1
         send_pairs(:, nsend) = send_pairs(:, i)
      end do
      allocate(send_count(rank_count), recv_count(rank_count), send_disp(rank_count), recv_disp(rank_count))
      send_count = 0
      allocate(send_cell(nsend), send_vertex(max_vertices, nsend), send_nnodes(nsend))
      send_vertex = 0
      do i = 1, nsend
         p = send_pairs(1, i) + 1; cell = send_pairs(2, i)
         send_count(p) = send_count(p) + 1; send_cell(i) = cell
         send_nnodes(i) = nnodes(cell)
         send_vertex(1:nnodes(cell), i) = cell_nodes(1:nnodes(cell), cell)
      end do
      deallocate(send_pairs)
      call MPI_Alltoall(send_count, 1, MPI_INTEGER, recv_count, 1, MPI_INTEGER, comm, ierr)
      call checked(ierr, 'counting output ghost cells')
      call displacements(send_count, send_disp, ns)
      call displacements(recv_count, recv_disp, ghost_count)
      allocate(metadata(metadata_width, ns), remote_metadata(metadata_width, ghost_count))
      metadata = 0
      do i = 1, ns
         cell = send_cell(i); k = nnodes(cell)
         metadata(1:3, i) = [global_cells(cell), vtk_types(cell), k]
         do j = 1, k
            idx = cell_nodes(j, cell)
            metadata(3 + j, i) = global_points(idx)
            metadata(3 + max_vertices + j, i) = local_owner(idx)
         end do
      end do
      call exchange_records(metadata_width, metadata, remote_metadata, send_count, recv_count, send_disp, recv_disp)
      deallocate(metadata)

      ! Sorted global IDs provide a single point for every shared vertex within
      ! a piece. Retain direct coordinate-source maps for all subsequent writes.
      wide = int(nused, int64) + sum(int(remote_metadata(3, :), int64))
      allocate(pairs(2, count32(wide)))
      pairs(1, 1:nused) = used_ids; pairs(2, :) = 0; k = nused
      do i = 1, ghost_count
         j = remote_metadata(3, i)
         if (j < 2 .or. j > max_vertices) call fail('invalid ghost vertex count')
         pairs(1, k + 1:k + j) = remote_metadata(4:3 + j, i)
         k = k + j
      end do
      call sort_pairs(pairs, 1, k)
      point_count = 0
      do i = 1, k
         if (point_count > 0) then
            if (pairs(1, i) == pairs(1, point_count)) cycle
         end if
         point_count = point_count + 1; pairs(1, point_count) = pairs(1, i)
      end do
      allocate(all_ids(point_count)); all_ids = pairs(1, 1:point_count)
      deallocate(pairs)
      allocate(point_global(point_count), point_ghost(point_count), point_local(point_count))
      allocate(point_remote_cell(point_count), point_remote_vertex(point_count), points(3, point_count))
      allocate(vertex_map(nlocal_points)); vertex_map = 0
      point_global = int(all_ids, int64); point_ghost = 1_int8; point_local = 0
      point_remote_cell = 0; point_remote_vertex = 0; points = 0.0_real64
      do i = 1, nused
         slot = find_id(all_ids, used_ids(i)); idx = used_local(i)
         vertex_map(idx) = slot; point_local(slot) = idx
         if (local_owner(idx) == rank_id) point_ghost(slot) = 0_int8
      end do
      total_cells = count32(int(owned_count, int64) + ghost_count)
      allocate(cell_global(total_cells), cell_owner(total_cells), cell_type(total_cells), cell_ghost(total_cells))
      allocate(offsets(total_cells))
      nconn = count32(sum(int(nnodes, int64)) + sum(int(remote_metadata(3, :), int64)))
      allocate(connectivity(nconn)); k = 0
      cell_global(1:owned_count) = int(global_cells, int64)
      cell_owner(1:owned_count) = int(rank_id, int32)
      cell_type(1:owned_count) = int(vtk_types, int8)
      cell_ghost(1:owned_count) = 0_int8
      cell_ghost(owned_count + 1:) = 1_int8
      do i = 1, owned_count
         do j = 1, nnodes(i)
            k = k + 1; connectivity(k) = int(vertex_map(cell_nodes(j, i)) - 1, int32)
         end do
         offsets(i) = int(k, int32)
      end do
      do p = 1, rank_count
         do i = recv_disp(p) + 1, recv_disp(p) + recv_count(p)
            cell = owned_count + i
            cell_global(cell) = int(remote_metadata(1, i), int64)
            cell_owner(cell) = int(p - 1, int32); cell_type(cell) = int(remote_metadata(2, i), int8)
            do j = 1, remote_metadata(3, i)
               slot = find_id(all_ids, remote_metadata(3 + j, i))
               if (slot == 0) call fail('missing ghost point')
               k = k + 1; connectivity(k) = int(slot - 1, int32)
               if (point_local(slot) == 0 .and. point_remote_cell(slot) == 0) then
                  point_remote_cell(slot) = i; point_remote_vertex(slot) = j
               end if
               if (remote_metadata(3 + max_vertices + j, i) == rank_id) point_ghost(slot) = 0_int8
            end do
            offsets(cell) = int(k, int32)
         end do
      end do
      allocate(send_xyz(3 * max_vertices, nsend), recv_xyz(3 * max_vertices, ghost_count))
      allocate(send_values(count32(int(nsend, int64) * max_fields)))
      allocate(recv_values(count32(int(ghost_count, int64) * max_fields)), field(total_cells))
      prepared = .true.
      call refresh_points(xyz)
   end subroutine

   subroutine scale_payload(width)
      integer, intent(in) :: width
      integer :: p
      do p = 1, rank_count
         payload_sc(p) = count32(int(width, int64) * send_count(p))
         payload_rc(p) = count32(int(width, int64) * recv_count(p))
         payload_sd(p) = count32(int(width, int64) * send_disp(p))
         payload_rd(p) = count32(int(width, int64) * recv_disp(p))
      end do
   end subroutine

   subroutine refresh_points(xyz)
      real(real64), intent(in) :: xyz(:, :)
      integer :: i, j, idx, d, ierr
      d = size(xyz, 1); send_xyz = 0.0_real64
      do i = 1, size(send_cell)
         do j = 1, send_nnodes(i)
            idx = 3 * (j - 1)
            send_xyz(idx + 1:idx + d, i) = xyz(:, send_vertex(j, i))
         end do
      end do
      call scale_payload(3 * max_vertices)
      call MPI_Alltoallv(send_xyz, payload_sc, payload_sd, MPI_DOUBLE_PRECISION, recv_xyz, payload_rc, payload_rd, &
         MPI_DOUBLE_PRECISION, output_comm, ierr)
      call checked(ierr, 'exchanging ghost coordinates')
      do i = 1, point_count
         if (point_local(i) > 0) then
            points(:, i) = 0.0_real64
            points(1:d, i) = xyz(:, point_local(i))
         else
            j = point_remote_cell(i); idx = 3 * (point_remote_vertex(i) - 1)
            if (j == 0) call fail('ghost point has no coordinate source')
            points(:, i) = recv_xyz(idx + 1:idx + 3, j)
         end if
      end do
   end subroutine

   function number(value) result(text)
      integer(int64), intent(in) :: value
      character(:), allocatable :: text
      character(32) :: buffer
      write(buffer, '(I0)') value
      text = trim(buffer)
   end function

   function xml_name(name) result(escaped)
      character(*), intent(in) :: name
      character(:), allocatable :: escaped
      integer :: i
      escaped = ''
      do i = 1, len_trim(name)
         select case(name(i:i))
         case('&')
            escaped = escaped//'&amp;'
         case('"')
            escaped = escaped//'&quot;'
         case('<')
            escaped = escaped//'&lt;'
         case('>')
            escaped = escaped//'&gt;'
         case default
            escaped = escaped//name(i:i)
         end select
      end do
   end function

   subroutine line(unit, text)
      integer, intent(in) :: unit
      character(*), intent(in) :: text
      write(unit) text//achar(10)
   end subroutine

   subroutine data_array(unit, name, kind_name, components, items, item_bytes, offset, extra)
      integer, intent(in) :: unit, components, item_bytes
      integer(int64), intent(in) :: items
      integer(int64), intent(inout) :: offset
      character(*), intent(in) :: name, kind_name
      character(*), intent(in), optional :: extra
      character(:), allocatable :: tail
      tail = ''
      if (present(extra)) tail = ' '//extra
      call line(unit, '<DataArray type="'//kind_name//'" Name="'//xml_name(name)//'" NumberOfComponents="'// &
         number(int(components, int64))//'" format="appended" offset="'//number(offset)//'"'//tail//'/>')
      offset = offset + 8_int64 + items * int(components, int64) * int(item_bytes, int64)
   end subroutine

   subroutine parallel_array(unit, name, kind_name, components)
      integer, intent(in) :: unit, components
      character(*), intent(in) :: name, kind_name
      call line(unit, '<PDataArray type="'//kind_name//'" Name="'//xml_name(name)//'" NumberOfComponents="'// &
         number(int(components, int64))//'"/>')
   end subroutine

   subroutine vtk_ghost_write(piece_prefix, meta_prefix, iteration, time, names, values, xyz, moving)
      character(*), intent(in) :: piece_prefix, meta_prefix, names(:)
      integer, intent(in) :: iteration
      real(real64), intent(in) :: time, values(:, :), xyz(:, :)
      logical, intent(in) :: moving
      integer :: nf, i, j, ierr, unit, status, p
      integer(int64) :: nc, np, offset
      integer(int8) :: endian_test(4)
      character(:), allocatable :: filename, suffix, byte_order
      character(32) :: time_text
      if (.not. prepared) call fail('output halo has not been initialized')
      nf = size(names)
      if (size(values, 1) /= owned_count .or. size(values, 2) /= nf .or. nf > field_capacity) &
         call fail('output field dimensions do not match the prepared mesh')
      if (moving) call refresh_points(xyz)
      do i = 1, size(send_cell)
         do j = 1, nf
            send_values((i - 1) * nf + j) = values(send_cell(i), j)
         end do
      end do
      call scale_payload(nf)
      call MPI_Alltoallv(send_values, payload_sc, payload_sd, MPI_DOUBLE_PRECISION, recv_values, payload_rc, payload_rd, &
         MPI_DOUBLE_PRECISION, output_comm, ierr)
      call checked(ierr, 'exchanging ghost output fields')
      nc = size(cell_global, kind=int64); np = int(point_count, int64)
      suffix = number(int(iteration, int64))
      filename = piece_prefix//suffix//'_'//number(int(rank_id, int64))//'.vtu'
      endian_test = transfer(1_int32, endian_test)
      byte_order = 'LittleEndian'
      if (endian_test(1) /= 1_int8) byte_order = 'BigEndian'
      open(newunit=unit, file=filename, access='stream', form='unformatted', status='replace', iostat=status)
      if (status /= 0) call fail('cannot open '//filename)
      call line(unit, '<VTKFile type="UnstructuredGrid" version="2.2" byte_order="'//byte_order// &
         '" header_type="UInt64">')
      call line(unit, '<UnstructuredGrid>')
      offset = 0_int64
      call line(unit, '<FieldData>')
      call data_array(unit, 'TimeValue', 'Float64', 1, 1_int64, 8, offset, 'NumberOfTuples="1"')
      call line(unit, '</FieldData>')
      call line(unit, '<Piece NumberOfPoints="'//number(np)//'" NumberOfCells="'//number(nc)//'">')
      call line(unit, '<CellData GlobalIds="GlobalCellId">')
      do j = 1, nf
         call data_array(unit, names(j), 'Float64', 1, nc, 8, offset)
      end do
      call data_array(unit, 'vtkGhostType', 'UInt8', 1, nc, 1, offset)
      call data_array(unit, 'GlobalCellId', 'Int64', 1, nc, 8, offset)
      call data_array(unit, 'OwnerRank', 'Int32', 1, nc, 4, offset)
      call line(unit, '</CellData>')
      call line(unit, '<PointData GlobalIds="GlobalPointId">')
      call data_array(unit, 'vtkGhostType', 'UInt8', 1, np, 1, offset)
      call data_array(unit, 'GlobalPointId', 'Int64', 1, np, 8, offset)
      call line(unit, '</PointData>')
      call line(unit, '<Points>')
      call data_array(unit, 'Coordinates', 'Float64', 3, np, 8, offset)
      call line(unit, '</Points>')
      call line(unit, '<Cells>')
      call data_array(unit, 'connectivity', 'Int32', 1, size(connectivity, kind=int64), 4, offset)
      call data_array(unit, 'offsets', 'Int32', 1, nc, 4, offset)
      call data_array(unit, 'types', 'UInt8', 1, nc, 1, offset)
      call line(unit, '</Cells>')
      call line(unit, '</Piece>')
      call line(unit, '</UnstructuredGrid>')
      call line(unit, '<AppendedData encoding="raw">')
      write(unit) '_'
      write(unit) 8_int64, time
      do j = 1, nf
         field(1:owned_count) = values(:, j)
         do i = 1, ghost_count
            field(owned_count + i) = recv_values((i - 1) * nf + j)
         end do
         write(unit) 8_int64 * nc, field
      end do
      write(unit) nc, cell_ghost
      write(unit) 8_int64 * nc, cell_global
      write(unit) 4_int64 * nc, cell_owner
      write(unit) np, point_ghost
      write(unit) 8_int64 * np, point_global
      write(unit) 24_int64 * np, points
      write(unit) 4_int64 * size(connectivity, kind=int64), connectivity
      write(unit) 4_int64 * nc, offsets
      write(unit) nc, cell_type
      call line(unit, achar(10)//'</AppendedData>')
      call line(unit, '</VTKFile>')
      close(unit, iostat=status)
      if (status /= 0) call fail('cannot finish '//filename)
      call MPI_Barrier(output_comm, ierr); call checked(ierr, 'completing VTU pieces')
      if (rank_id == 0) then
         filename = meta_prefix//suffix//'.pvtu'
         open(newunit=unit, file=filename, access='stream', form='unformatted', status='replace', iostat=status)
         if (status /= 0) call fail('cannot open '//filename)
         call line(unit, '<VTKFile type="PUnstructuredGrid" version="2.2" byte_order="'//byte_order// &
            '" header_type="UInt64">')
         call line(unit, '<PUnstructuredGrid GhostLevel="1">')
         call line(unit, '<FieldData>')
         write(time_text, '(ES24.16E3)') time
         call line(unit, '<DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'// &
            trim(adjustl(time_text))//'</DataArray>')
         call line(unit, '</FieldData>')
         call line(unit, '<PCellData GlobalIds="GlobalCellId">')
         do j = 1, nf
            call parallel_array(unit, names(j), 'Float64', 1)
         end do
         call parallel_array(unit, 'vtkGhostType', 'UInt8', 1)
         call parallel_array(unit, 'GlobalCellId', 'Int64', 1)
         call parallel_array(unit, 'OwnerRank', 'Int32', 1)
         call line(unit, '</PCellData>')
         call line(unit, '<PPointData GlobalIds="GlobalPointId">')
         call parallel_array(unit, 'vtkGhostType', 'UInt8', 1)
         call parallel_array(unit, 'GlobalPointId', 'Int64', 1)
         call line(unit, '</PPointData>')
         call line(unit, '<PPoints>')
         call parallel_array(unit, 'Coordinates', 'Float64', 3)
         call line(unit, '</PPoints>')
         do p = 0, rank_count - 1
            call line(unit, '<Piece Source="'//piece_prefix//suffix//'_'//number(int(p, int64))//'.vtu"/>')
         end do
         call line(unit, '</PUnstructuredGrid>')
         call line(unit, '</VTKFile>')
         close(unit, iostat=status)
         if (status /= 0) call fail('cannot finish '//filename)
      end if
      call MPI_Barrier(output_comm, ierr); call checked(ierr, 'completing PVTU metadata')
   end subroutine

   subroutine vtk_ghost_release()
      if (.not. prepared) return
      deallocate(send_count, recv_count, send_disp, recv_disp, send_cell, send_vertex, send_nnodes)
      deallocate(point_local, point_remote_cell, point_remote_vertex, payload_sc, payload_rc, payload_sd, payload_rd)
      deallocate(connectivity, offsets, cell_owner, cell_global, point_global, cell_type, cell_ghost, point_ghost)
      deallocate(points, send_xyz, recv_xyz, send_values, recv_values, field)
      prepared = .false.
   end subroutine
