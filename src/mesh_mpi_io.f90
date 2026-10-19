! Collective mesh input for native sequential-unformatted GRID files.
! One contiguous block per rank; cached records are routed to their consumers.
module mesh_mpi_io
   use mpi
   use iso_fortran_env, only: int8, int32, int64
   implicit none
   private
   public :: mesh_enabled, mesh_cells, mesh_boundaries, mesh_boundary_count, mesh_wallxyz, mesh_wall_count
   public :: mesh_distributed, mesh_owned_count, mesh_store_partition, mesh_query_cells, mesh_owner_release
   public :: route_integer, route_real
   public :: mesh_node_owners
   public :: mesh_start, mesh_shapes, mesh_cell_nodes, mesh_prepare, mesh_vertices, mesh_incidence
   public :: mesh_load_boundaries, mesh_release_boundaries, mesh_release_vertices, mesh_finish, mesh_fail
   logical :: mesh_enabled = .false., mesh_distributed = .false.
   integer :: mesh_owned_count = 0, cell_first
   integer, allocatable :: directory_owner(:), directory_index(:)
   integer, allocatable :: mesh_cells(:, :), mesh_boundaries(:, :)
   real, allocatable :: mesh_wallxyz(:, :)
   integer :: mesh_boundary_count = 0, mesh_wall_count = 0
   integer :: rank_id, rank_count, ndim, totals(3), markers(3), strides(3)
   integer :: vertex_first, vertex_count
   integer, allocatable :: vertex_ptr(:), vertex_cells(:)
   real, allocatable :: vertex_xyz(:, :)
   integer(int64) :: disk_records(3) = 0
   real :: read_seconds(3) = 0.0
   integer, parameter :: batch_records = 65536
   character(len=8), parameter :: names(3) = ['GRID.cel', 'GRID.vrt', 'GRID.bnd']
contains
   subroutine mesh_fail(message)
      character(len=*), intent(in) :: message
      integer :: ierr
      write(*,'(A,I0,2A)') 'MPI mesh error on rank ', rank_id, ': ', trim(message)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      error stop
   end subroutine

   subroutine checked(ierr, operation)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation
      if (ierr /= MPI_SUCCESS) call mesh_fail(operation)
   end subroutine

   subroutine block_range(total, first, count)
      integer, intent(in) :: total
      integer, intent(out) :: first, count
      integer :: base
      base = total / rank_count
      first = rank_id * base + 1
      count = base
      if (rank_id == rank_count - 1) count = total - base * (rank_count - 1)
   end subroutine

   integer function vertex_owner(id) result(owner)
      integer, intent(in) :: id
      integer :: base
      if (id < 1 .or. id > totals(2)) call mesh_fail('vertex ID outside GRID.vrt')
      base = totals(2) / rank_count
      owner = rank_count - 1
      if (base > 0) owner = min((id - 1) / base, rank_count - 1)
   end function

   integer(int64) function record_marker(bytes) result(value)
      integer(int8), intent(in) :: bytes(:)
      value = -1_int64
      if (size(bytes) == 4) value = int(transfer(bytes, 0_int32), int64)
      if (size(bytes) == 8) value = transfer(bytes, 0_int64)
   end function

   subroutine probe_file(which, payload, valid)
      integer, intent(in) :: which, payload
      logical, intent(out) :: valid
      integer :: fh, ierr, status(MPI_STATUS_SIZE), m, got, trial, nread
      integer(kind=MPI_OFFSET_KIND) :: file_bytes, offset
      integer(int8) :: head(128), tail(128)
      valid = .false.
      call MPI_File_open(MPI_COMM_WORLD, names(which), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
      call checked(ierr, 'opening '//names(which))
      call MPI_File_get_size(fh, file_bytes, ierr)
      call checked(ierr, 'getting size of '//names(which))
      markers(which) = 0
      totals(which) = 0
      if (rank_id == 0) then
         if (file_bytes == 0 .and. which == 3) then
            markers(which) = 4
         else if (file_bytes >= payload + 8) then
            nread = int(min(file_bytes, int(size(head), MPI_OFFSET_KIND)))
            call MPI_File_read_at(fh, 0_MPI_OFFSET_KIND, head, nread, MPI_BYTE, status, ierr)
            call checked(ierr, 'probing '//names(which))
            call MPI_Get_count(status, MPI_BYTE, got, ierr)
            if (got /= nread) call mesh_fail('short metadata read')
            do m = 4, 8, 4
               trial = payload + 2 * m
               if (nread < trial) cycle
               if (mod(file_bytes, int(trial, MPI_OFFSET_KIND)) /= 0) cycle
               if (record_marker(head(1:m)) /= payload) cycle
               if (record_marker(head(m + payload + 1:trial)) /= payload) cycle
               if (transfer(head(m + 1:m + 4), 0_int32) /= 1) cycle
               if (file_bytes / trial > int(huge(0), MPI_OFFSET_KIND)) cycle
               offset = file_bytes - trial
               call MPI_File_read_at(fh, offset, tail, trial, MPI_BYTE, status, ierr)
               call checked(ierr, 'probing final record in '//names(which))
               call MPI_Get_count(status, MPI_BYTE, got, ierr)
               if (got /= trial) call mesh_fail('short final metadata read')
               if (record_marker(tail(1:m)) /= payload) cycle
               if (record_marker(tail(m + payload + 1:trial)) /= payload) cycle
               if (transfer(tail(m + 1:m + 4), 0_int32) /= int(file_bytes / trial)) cycle
               markers(which) = m
               totals(which) = int(file_bytes / trial)
               exit
            end do
         end if
      end if
      call MPI_File_close(fh, ierr)
      call checked(ierr, 'closing '//names(which))
      call MPI_Bcast(markers(which), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'broadcasting record layout')
      call MPI_Bcast(totals(which), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'broadcasting mesh count')
      strides(which) = payload + 2 * markers(which)
      valid = markers(which) > 0
   end subroutine

   subroutine mesh_start(dimension, binary_io, partitioner, nc, nv, nb)
      integer, intent(in) :: dimension, binary_io, partitioner
      integer, intent(out) :: nc, nv, nb
      integer :: ierr, payloads(3), k, first, count
      logical :: valid, all_valid
      integer, allocatable :: records(:, :)
      real :: started
      started = MPI_Wtime()
      call MPI_Comm_rank(MPI_COMM_WORLD, rank_id, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)
      ! Supported meshes always use collective input; no runtime selector is required.
      mesh_enabled = .false.
      nc = 0; nv = 0; nb = 0
      valid = binary_io == 1 .and. partitioner == 6 .and. (dimension == 2 .or. dimension == 3)
      valid = valid .and. storage_size(0) == 32 .and. storage_size(0.0) == 64
      if (.not. valid) then
         if (rank_id == 0) write(*,'(A)') 'Mesh input: legacy reader (configuration unsupported by MPI reader)'
         return
      end if
      ndim = dimension
      payloads = [4 * (1 + 2**ndim), 4 + 8 * ndim, 24]
      all_valid = .true.
      do k = 1, 3
         call probe_file(k, payloads(k), valid)
         all_valid = all_valid .and. valid
      end do
      if (.not. all_valid) then
         if (rank_id == 0) write(*,'(A)') 'Mesh input: legacy reader (GRID record layout not recognised)'
         return
      end if
      mesh_enabled = .true.
      nc = totals(1); nv = totals(2); nb = totals(3)
      call read_records(1, records)
      call block_range(nc, first, count)
      allocate(mesh_cells(11, count))
      mesh_cells = 0
      mesh_cells(1:size(records, 1), :) = records
      if (any(records(2:, :) < 1) .or. any(records(2:, :) > nv)) call mesh_fail('invalid cell vertex ID')
      if (rank_id == 0) then
         write(*,'(A,3(I0,1X),A,F10.3)') 'MPI mesh cells/vertices/boundaries: ', nc, nv, nb, &
            ' metadata + cell read seconds: ', MPI_Wtime() - started
      end if
   end subroutine

   subroutine read_records(which, records)
      integer, intent(in) :: which
      integer, allocatable, intent(out) :: records(:, :)
      integer :: fh, ierr, status(MPI_STATUS_SIZE), first, count, max_count, start, take, stride, m, words
      integer :: i, j, pos, got, payload
      integer(kind=MPI_OFFSET_KIND) :: offset
      integer(int8), allocatable :: raw(:)
      real :: started
      started = MPI_Wtime()
      call block_range(totals(which), first, count)
      stride = strides(which); m = markers(which); payload = stride - 2 * m
      words = payload / 4
      allocate(records(words, count), raw(max(1, min(batch_records, count) * stride)))
      call MPI_File_open(MPI_COMM_WORLD, names(which), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
      call checked(ierr, 'opening '//names(which))
      call MPI_Allreduce(count, max_count, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'agreeing collective read rounds')
      do start = 0, max_count - 1, batch_records
         take = max(0, min(batch_records, count - start))
         offset = (int(first - 1, MPI_OFFSET_KIND) + int(min(start, count), MPI_OFFSET_KIND)) * stride
         call MPI_File_read_at_all(fh, offset, raw, take * stride, MPI_BYTE, status, ierr)
         call checked(ierr, 'collective read of '//names(which))
         call MPI_Get_count(status, MPI_BYTE, got, ierr)
         if (got /= take * stride) call mesh_fail('short collective read of '//names(which))
         do i = 1, take
            pos = (i - 1) * stride
            if (record_marker(raw(pos + 1:pos + m)) /= payload) call mesh_fail('invalid leading record marker')
            if (record_marker(raw(pos + m + payload + 1:pos + stride)) /= payload) then
               call mesh_fail('invalid trailing record marker')
            end if
            do j = 1, words
               records(j, start + i) = transfer(raw(pos + m + 4 * (j - 1) + 1:pos + m + 4 * j), 0_int32)
            end do
            if (records(1, start + i) /= first + start + i - 1) call mesh_fail('GRID IDs must be consecutive from 1')
         end do
      end do
      call MPI_File_close(fh, ierr)
      call checked(ierr, 'closing '//names(which))
      disk_records(which) = disk_records(which) + int(count, int64)
      read_seconds(which) = read_seconds(which) + MPI_Wtime() - started
   end subroutine

   subroutine mesh_cell_nodes(record, nodes, nn, shape, faces)
      integer, intent(in) :: record(:)
      integer, intent(out) :: nodes(8), nn, shape, faces
      integer :: v(8)
      v = 0; v(1:2**ndim) = record(2:1 + 2**ndim)
      shape = 0; nodes = 0
      if (ndim == 2) then
         shape = 5; nn = 4; faces = 4
         if (v(3) == v(4)) then
            shape = 6; nn = 3; faces = 3
         end if
         nodes(1:nn) = v(1:nn)
      else
         if (v(5) /= v(6) .and. v(6) /= v(7) .and. v(7) /= v(8) .and. v(3) /= v(4)) then
            shape = 1; nn = 8; faces = 6; nodes = v
         end if
         if (v(3) == v(4) .and. v(5) == v(6) .and. v(6) == v(7) .and. v(7) == v(8)) then
            shape = 2; nn = 4; faces = 4; nodes(1:4) = [v(1:3), v(5)]
         end if
         if (v(3) /= v(4) .and. v(5) == v(6) .and. v(6) == v(7) .and. v(7) == v(8)) then
            shape = 3; nn = 5; faces = 5; nodes(1:5) = v(1:5)
         end if
         if (v(3) == v(4) .and. v(5) /= v(6) .and. v(7) == v(8)) then
            shape = 4; nn = 6; faces = 5; nodes(1:6) = [v(1:3), v(5:7)]
         end if
      end if
      if (shape == 0) call mesh_fail('unsupported cell connectivity')
   end subroutine

   subroutine mesh_shapes(shapes, global_counts)
      integer, intent(out) :: shapes(:), global_counts(6)
      integer :: i, nn, faces, nodes(8), counts(6), ierr
      integer(int64) :: local_nodes, offset, total_nodes
      counts = 0; local_nodes = 0
      do i = 1, size(mesh_cells, 2)
         call mesh_cell_nodes(mesh_cells(:, i), nodes, nn, shapes(i), faces)
         counts(shapes(i)) = counts(shapes(i)) + 1
         local_nodes = local_nodes + nn
      end do
      call MPI_Allreduce(counts, global_counts, 6, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'reducing shape counts')
      call MPI_Exscan(local_nodes, offset, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'computing connectivity offsets')
      if (rank_id == 0) offset = 0
      call MPI_Allreduce(local_nodes, total_nodes, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'checking connectivity offset range')
      if (total_nodes > int(huge(0), int64)) call mesh_fail('solver connectivity offsets exceed int32 capacity')
      do i = 1, size(mesh_cells, 2)
         call mesh_cell_nodes(mesh_cells(:, i), nodes, nn, shapes(i), faces)
         mesh_cells(10, i) = int(offset)
         offset = offset + nn
         mesh_cells(11, i) = int(offset)
      end do
   end subroutine

   integer function cell_directory(id) result(owner)
      integer, intent(in) :: id
      integer :: base
      if (id < 1 .or. id > totals(1)) call mesh_fail('cell ID outside GRID.cel')
      base = totals(1) / rank_count
      owner = rank_count - 1
      if (base > 0) owner = min((id - 1) / base, rank_count - 1)
   end function

   subroutine mesh_store_partition(parts)
      integer, intent(in) :: parts(:)
      integer, allocatable :: counts(:), received(:)
      integer :: i, count, ierr
      mesh_distributed = .true.
      call block_range(totals(1), cell_first, count)
      if (size(parts) /= count) call mesh_fail('partition directory block mismatch')
      directory_owner = parts
      allocate(directory_index(count)); directory_index = 0
      allocate(counts(0:rank_count - 1), received(0:rank_count - 1)); counts = 0
      do i = 1, count
         if (parts(i) < 0 .or. parts(i) >= rank_count) call mesh_fail('partition returned invalid owner')
         counts(parts(i)) = counts(parts(i)) + 1
      end do
      call MPI_Alltoall(counts, 1, MPI_INTEGER, received, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'counting partition-owned cells')
      mesh_owned_count = sum(received)
   end subroutine

   subroutine mesh_register_cells()
      integer, allocatable :: records(:, :), destinations(:), received(:, :)
      integer :: i, slot
      allocate(records(3, size(mesh_cells, 2)), destinations(size(mesh_cells, 2)))
      do i = 1, size(mesh_cells, 2)
         records(:, i) = [mesh_cells(1, i), rank_id, i]
         destinations(i) = cell_directory(mesh_cells(1, i))
      end do
      call route_integer(records, destinations, received)
      do i = 1, size(received, 2)
         slot = received(1, i) - cell_first + 1
         if (directory_owner(slot) /= received(2, i)) call mesh_fail('registered cell owner mismatch')
         directory_index(slot) = received(3, i)
      end do
      if (any(directory_index == 0)) call mesh_fail('missing cell registration')
   end subroutine

   subroutine mesh_query_cells(ids, owners, indices)
      integer, intent(in) :: ids(:)
      integer, allocatable, intent(out) :: owners(:), indices(:)
      integer, allocatable :: requests(:, :), destinations(:), received(:, :), replies(:, :), returned(:, :)
      integer :: i, slot
      allocate(requests(3, size(ids)), destinations(size(ids)))
      do i = 1, size(ids)
         requests(:, i) = [ids(i), i, rank_id]; destinations(i) = cell_directory(ids(i))
      end do
      call route_integer(requests, destinations, received)
      deallocate(requests, destinations)
      allocate(replies(3, size(received, 2)), destinations(size(received, 2)))
      do i = 1, size(received, 2)
         slot = received(1, i) - cell_first + 1
         replies(:, i) = [received(2, i), directory_owner(slot), directory_index(slot)]
         destinations(i) = received(3, i)
      end do
      call route_integer(replies, destinations, returned)
      allocate(owners(size(ids)), indices(size(ids)))
      do i = 1, size(returned, 2)
         owners(returned(1, i)) = returned(2, i); indices(returned(1, i)) = returned(3, i)
      end do
   end subroutine

   subroutine mesh_owner_release()
      if (allocated(directory_owner)) deallocate(directory_owner, directory_index)
   end subroutine

   subroutine mesh_prepare(owners)
      integer, intent(in) :: owners(:)
      integer, allocatable :: pairs(:, :), destinations(:), received(:, :), cursor(:)
      integer :: i, j, k, nv, slot
      nv = 2**ndim
      if (int(size(mesh_cells, 2), int64) * nv > huge(0)) call mesh_fail('local incidence count exceeds int32')
      allocate(pairs(2, nv * size(mesh_cells, 2)), destinations(nv * size(mesh_cells, 2)))
      k = 0
      do i = 1, size(mesh_cells, 2)
         do j = 1, nv
            k = k + 1
            pairs(:, k) = [mesh_cells(j + 1, i), mesh_cells(1, i)]
            destinations(k) = vertex_owner(pairs(1, k))
         end do
      end do
      call route_integer(pairs, destinations, received)
      deallocate(pairs, destinations)
      call block_range(totals(2), vertex_first, vertex_count)
      allocate(vertex_ptr(vertex_count + 1), cursor(vertex_count))
      vertex_ptr = 0
      do i = 1, size(received, 2)
         slot = received(1, i) - vertex_first + 1
         vertex_ptr(slot + 1) = vertex_ptr(slot + 1) + 1
      end do
      vertex_ptr(1) = 1
      do i = 1, vertex_count
         vertex_ptr(i + 1) = vertex_ptr(i + 1) + vertex_ptr(i)
      end do
      allocate(vertex_cells(size(received, 2)))
      cursor = vertex_ptr(1:vertex_count)
      do i = 1, size(received, 2)
         slot = received(1, i) - vertex_first + 1
         vertex_cells(cursor(slot)) = received(2, i)
         cursor(slot) = cursor(slot) + 1
      end do
      deallocate(received, cursor)
      allocate(destinations(size(mesh_cells, 2)))
      do i = 1, size(mesh_cells, 2)
         if (mesh_distributed) then
            destinations(i) = directory_owner(i)
         else
            destinations(i) = owners(mesh_cells(1, i))
         end if
      end do
      call route_integer(mesh_cells, destinations, received)
      call move_alloc(received, mesh_cells)
      if (mesh_distributed) call mesh_register_cells()
      call load_vertices()
   end subroutine

   subroutine load_vertices()
      integer, allocatable :: records(:, :)
      integer :: i, d
      if (allocated(vertex_xyz)) return
      call read_records(2, records)
      call block_range(totals(2), vertex_first, vertex_count)
      allocate(vertex_xyz(3, vertex_count)); vertex_xyz = 0.0
      do i = 1, vertex_count
         do d = 1, ndim
            vertex_xyz(d, i) = transfer(records(2 * d:2 * d + 1, i), 0.0)
         end do
      end do
   end subroutine

   subroutine mesh_vertices(ids, xyz)
      integer, intent(in) :: ids(:)
      real, allocatable, intent(out) :: xyz(:, :)
      integer, allocatable :: requests(:, :), destinations(:), received(:, :), returned(:, :), indices(:, :)
      real, allocatable :: values(:, :), replies(:, :)
      integer :: i
      allocate(requests(3, size(ids)), destinations(size(ids)))
      do i = 1, size(ids)
         requests(:, i) = [ids(i), i, rank_id]
         destinations(i) = vertex_owner(ids(i))
      end do
      call route_integer(requests, destinations, received)
      deallocate(requests, destinations)
      allocate(values(3, size(received, 2)), indices(1, size(received, 2)), destinations(size(received, 2)))
      do i = 1, size(received, 2)
         values(:, i) = vertex_xyz(:, received(1, i) - vertex_first + 1)
         indices(1, i) = received(2, i)
         destinations(i) = received(3, i)
      end do
      call route_integer(indices, destinations, returned)
      call route_real(values, destinations, replies)
      allocate(xyz(3, size(ids)))
      do i = 1, size(returned, 2)
         xyz(:, returned(1, i)) = replies(:, i)
      end do
   end subroutine

   subroutine mesh_incidence(ids, ptr, cells)
      integer, intent(in) :: ids(:)
      integer, allocatable, intent(out) :: ptr(:), cells(:)
      integer, allocatable :: requests(:, :), destinations(:), received(:, :), replies(:, :), returned(:, :), cursor(:)
      integer :: i, j, k, slot, total
      integer(int64) :: wide_total
      allocate(requests(3, size(ids)), destinations(size(ids)))
      do i = 1, size(ids)
         requests(:, i) = [ids(i), i, rank_id]
         destinations(i) = vertex_owner(ids(i))
      end do
      call route_integer(requests, destinations, received)
      deallocate(requests, destinations)
      wide_total = 0
      do i = 1, size(received, 2)
         slot = received(1, i) - vertex_first + 1
         wide_total = wide_total + vertex_ptr(slot + 1) - vertex_ptr(slot)
      end do
      if (wide_total >= huge(0)) call mesh_fail('local incidence replies exceed int32')
      total = int(wide_total)
      allocate(replies(2, total), destinations(total))
      k = 0
      do i = 1, size(received, 2)
         slot = received(1, i) - vertex_first + 1
         do j = vertex_ptr(slot), vertex_ptr(slot + 1) - 1
            k = k + 1
            replies(:, k) = [received(2, i), vertex_cells(j)]
            destinations(k) = received(3, i)
         end do
      end do
      call route_integer(replies, destinations, returned)
      allocate(ptr(size(ids) + 1), cursor(size(ids)))
      ptr = 0
      do i = 1, size(returned, 2)
         ptr(returned(1, i) + 1) = ptr(returned(1, i) + 1) + 1
      end do
      ptr(1) = 1
      do i = 1, size(ids)
         ptr(i + 1) = ptr(i + 1) + ptr(i)
      end do
      cursor = ptr(1:size(ids))
      allocate(cells(size(returned, 2)))
      do i = 1, size(returned, 2)
         slot = returned(1, i)
         cells(cursor(slot)) = returned(2, i)
         cursor(slot) = cursor(slot) + 1
      end do
   end subroutine

   subroutine mesh_load_boundaries(owners, scale, need_walls)
      integer, intent(in) :: owners(:)
      real, intent(in) :: scale
      logical, intent(in) :: need_walls
      integer, allocatable :: records(:, :), ptr(:), cells(:), seen(:), destinations(:), packed(:, :)
      integer, allocatable :: ids(:), counts(:), disps(:), cell_owners(:), cell_indices(:)
      real, allocatable :: xyz(:, :), wall_local(:, :)
      integer :: i, j, k, peer, total, ierr, nw, prefix, ordinal, nall, stride, wi
      integer :: wall_local_count, dtype
      call read_records(3, records)
      if (any(records(2:1 + ndim, :) < 1) .or. any(records(2:1 + ndim, :) > totals(2))) then
         call mesh_fail('invalid boundary vertex ID')
      end if
      if (ndim == 3) then
         if (any(records(5, :) < 1) .or. any(records(5, :) > totals(2))) call mesh_fail('invalid boundary vertex ID')
      end if
      wall_local_count = count(records(6, :) == 4)
      nw = wall_local_count
      if (ndim == 2) nw = nw + count(records(6, :) == 99)
      ! Count disjoint file blocks before routing can duplicate boundary records.
      ! The global count also sizes wall output when no wall-distance cache is needed.
      call MPI_Allreduce(nw, mesh_wall_count, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'global wall count')
      call MPI_Exscan(nw, prefix, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'wall numbering prefix')
      if (rank_id == 0) prefix = 0
      call mesh_incidence(records(2, :), ptr, cells)
      if (mesh_distributed) then
         call mesh_query_cells(cells, cell_owners, cell_indices)
      else
         cell_owners = owners(cells)
      end if
      allocate(seen(0:rank_count - 1)); seen = 0
      total = 0
      do i = 1, size(records, 2)
         do j = ptr(i), ptr(i + 1) - 1
            peer = cell_owners(j)
            if (seen(peer) == i) cycle
            seen(peer) = i
            if (total == huge(total)) call mesh_fail('boundary routing count overflow')
            total = total + 1
         end do
      end do
      allocate(packed(7, total), destinations(total)); seen = 0
      k = 0; ordinal = prefix
      do i = 1, size(records, 2)
         if (records(6, i) == 4 .or. (ndim == 2 .and. records(6, i) == 99)) ordinal = ordinal + 1
         do j = ptr(i), ptr(i + 1) - 1
            peer = cell_owners(j)
            if (seen(peer) == i) cycle
            seen(peer) = i
            k = k + 1
            packed(1:6, k) = records(:, i); packed(7, k) = ordinal
            destinations(k) = peer
         end do
      end do
      call route_integer(packed, destinations, mesh_boundaries)
      mesh_boundary_count = size(mesh_boundaries, 2)
      deallocate(packed, destinations, ptr, cells, seen)
      if (.not. need_walls) return
      ! Match the existing wall-centre arithmetic and global file ordering.
      stride = 2
      if (ndim == 3) stride = 4
      allocate(ids(stride * wall_local_count))
      k = 0
      do i = 1, size(records, 2)
         if (records(6, i) /= 4) cycle
         ids(k + 1:k + stride) = records(2:stride + 1, i)
         k = k + stride
      end do
      call mesh_vertices(ids, xyz)
      if (ndim == 3) xyz = xyz / scale
      allocate(wall_local(3, wall_local_count)); wall_local = 0.0
      k = 0; wi = 0
      do i = 1, size(records, 2)
         if (records(6, i) /= 4) cycle
         wi = wi + 1
         if (ndim == 2) then
            wall_local(:, wi) = (xyz(:, k + 1) + xyz(:, k + 2)) / 2.0
         else if (records(4, i) == records(5, i)) then
            wall_local(:, wi) = (xyz(:, k + 1) + xyz(:, k + 2) + xyz(:, k + 3)) / 3.0
         else
            wall_local(:, wi) = (xyz(:, k + 1) + xyz(:, k + 2) + xyz(:, k + 3) + xyz(:, k + 4)) / 4.0
         end if
         k = k + stride
      end do
      allocate(counts(0:rank_count - 1), disps(0:rank_count - 1))
      call MPI_Allgather(wall_local_count, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'wall counts')
      call offsets(counts, disps, nall)
      allocate(mesh_wallxyz(3, nall))
      call MPI_Type_contiguous(3, MPI_DOUBLE_PRECISION, dtype, ierr)
      call MPI_Type_commit(dtype, ierr)
      call MPI_Allgatherv(wall_local, wall_local_count, dtype, mesh_wallxyz, counts, disps, dtype, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'gathering wall centres')
      call MPI_Type_free(dtype, ierr)
   end subroutine

   subroutine mesh_node_owners(ids, owners)
      integer, intent(in) :: ids(:)
      integer, allocatable, intent(out) :: owners(:)
      integer, allocatable :: requests(:, :), destinations(:), received(:, :), replies(:, :), returned(:, :), maximum(:)
      integer :: i, slot
      allocate(requests(3, size(ids)), destinations(size(ids)))
      do i = 1, size(ids)
         requests(:, i) = [ids(i), i, rank_id]; destinations(i) = vertex_owner(ids(i))
      end do
      call route_integer(requests, destinations, received)
      deallocate(requests, destinations)
      allocate(maximum(vertex_count)); maximum = -1
      do i = 1, size(received, 2)
         slot = received(1, i) - vertex_first + 1
         maximum(slot) = max(maximum(slot), received(3, i))
      end do
      allocate(replies(2, size(received, 2)), destinations(size(received, 2)))
      do i = 1, size(received, 2)
         slot = received(1, i) - vertex_first + 1
         replies(:, i) = [received(2, i), maximum(slot)]; destinations(i) = received(3, i)
      end do
      call route_integer(replies, destinations, returned)
      allocate(owners(size(ids)))
      do i = 1, size(returned, 2)
         owners(returned(1, i)) = returned(2, i)
      end do
   end subroutine

   subroutine mesh_release_boundaries()
      if (allocated(mesh_boundaries)) deallocate(mesh_boundaries)
      if (allocated(vertex_ptr)) deallocate(vertex_ptr, vertex_cells)
      mesh_boundary_count = 0
   end subroutine

   subroutine mesh_release_vertices()
      if (allocated(vertex_xyz)) deallocate(vertex_xyz)
   end subroutine

   subroutine mesh_finish()
      integer :: ierr
      integer(int64) :: records(3)
      real :: seconds(3)
      if (mesh_enabled) then
         call MPI_Reduce(disk_records, records, 3, MPI_INTEGER8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         call checked(ierr, 'reducing mesh read counters')
         call MPI_Reduce(read_seconds, seconds, 3, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
         call checked(ierr, 'reducing mesh read timings')
         if (rank_id == 0) then
            write(*, '(A,3(I0,1X))') 'MPI mesh records read, cells/vertices/boundaries: ', records
            write(*, '(A,3(F12.4,1X))') 'MPI mesh read/decode maximum seconds: ', seconds
         end if
      end if
      if (allocated(mesh_cells)) deallocate(mesh_cells)
      if (allocated(vertex_xyz)) deallocate(vertex_xyz)
      if (allocated(mesh_wallxyz)) deallocate(mesh_wallxyz)
      mesh_wall_count = 0
      call mesh_release_boundaries()
   end subroutine

   subroutine offsets(counts, disps, total)
      integer, intent(in) :: counts(0:)
      integer, intent(out) :: disps(0:), total
      integer :: i
      integer(int64) :: wide
      wide = 0
      do i = 0, size(counts) - 1
         if (wide > huge(0)) call mesh_fail('MPI record displacement exceeds int32')
         disps(i) = int(wide)
         wide = wide + counts(i)
      end do
      if (wide >= huge(0)) call mesh_fail('MPI local record count exceeds int32')
      total = int(wide)
   end subroutine

   subroutine route_integer(values, destinations, received)
      integer, intent(in) :: values(:, :)
      integer, intent(in) :: destinations(:)
      integer, allocatable, intent(out) :: received(:, :)
      integer, allocatable :: packed(:, :)
      integer, allocatable :: sc(:), rc(:), sd(:), rd(:), cursor(:)
      integer :: ierr, dtype, i, peer, ns, nr, width
      width = size(values, 1)
      allocate(sc(0:rank_count - 1), rc(0:rank_count - 1), sd(0:rank_count - 1), rd(0:rank_count - 1))
      allocate(cursor(0:rank_count - 1))
      sc = 0
      do i = 1, size(destinations)
         peer = destinations(i)
         if (peer < 0 .or. peer >= rank_count) call mesh_fail('invalid routing rank')
         sc(peer) = sc(peer) + 1
      end do
      call MPI_Alltoall(sc, 1, MPI_INTEGER, rc, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'exchanging record counts')
      call offsets(sc, sd, ns)
      call offsets(rc, rd, nr)
      allocate(packed(width, ns), received(width, nr))
      cursor = sd
      do i = 1, size(destinations)
         peer = destinations(i)
         cursor(peer) = cursor(peer) + 1
         packed(:, cursor(peer)) = values(:, i)
      end do
      call MPI_Type_contiguous(width, MPI_INTEGER, dtype, ierr)
      call checked(ierr, 'creating record datatype')
      call MPI_Type_commit(dtype, ierr)
      call checked(ierr, 'committing record datatype')
      call MPI_Alltoallv(packed, sc, sd, dtype, received, rc, rd, dtype, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'routing mesh records')
      call MPI_Type_free(dtype, ierr)
   end subroutine

   subroutine route_real(values, destinations, received)
      real, intent(in) :: values(:, :)
      integer, intent(in) :: destinations(:)
      real, allocatable, intent(out) :: received(:, :)
      real, allocatable :: packed(:, :)
      integer, allocatable :: sc(:), rc(:), sd(:), rd(:), cursor(:)
      integer :: ierr, dtype, i, peer, ns, nr, width
      width = size(values, 1)
      allocate(sc(0:rank_count - 1), rc(0:rank_count - 1), sd(0:rank_count - 1), rd(0:rank_count - 1))
      allocate(cursor(0:rank_count - 1))
      sc = 0
      do i = 1, size(destinations)
         peer = destinations(i)
         if (peer < 0 .or. peer >= rank_count) call mesh_fail('invalid routing rank')
         sc(peer) = sc(peer) + 1
      end do
      call MPI_Alltoall(sc, 1, MPI_INTEGER, rc, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'exchanging record counts')
      call offsets(sc, sd, ns)
      call offsets(rc, rd, nr)
      allocate(packed(width, ns), received(width, nr))
      cursor = sd
      do i = 1, size(destinations)
         peer = destinations(i)
         cursor(peer) = cursor(peer) + 1
         packed(:, cursor(peer)) = values(:, i)
      end do
      call MPI_Type_contiguous(width, MPI_DOUBLE_PRECISION, dtype, ierr)
      call checked(ierr, 'creating record datatype')
      call MPI_Type_commit(dtype, ierr)
      call checked(ierr, 'committing record datatype')
      call MPI_Alltoallv(packed, sc, sd, dtype, received, rc, rd, dtype, MPI_COMM_WORLD, ierr)
      call checked(ierr, 'routing mesh records')
      call MPI_Type_free(dtype, ierr)
   end subroutine
end module mesh_mpi_io
