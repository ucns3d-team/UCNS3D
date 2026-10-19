! Host-only ownership and geometry cache. MPI requests occur in setup batches;
! lookup functions remain communication-free inside OpenMP loops.
module mesh_local_cells
   use declaration, only: xmpie, xmpil, ielem_ihexgl, ielem_ineighg, ielem_ifca, ielem_xxc, ielem_yyc, ielem_zzc, &
      dims, iselemt, ilocalallelg, lowmem, centerr, xand2r, cand, cand2r, glneigh
   use mpiinfo, only: n
   use mesh_mpi_io
   use mesh_local_nodes, only: node_sort
   use mpi
   use iso_fortran_env, only: int64
   implicit none
   private
   public :: cells_seed, cells_fetch, cells_geometry, cells_stencils, cells_release
   public :: cell_cache_slot
   public :: cell_owner, cell_index, cell_centre, cell_centres, cell_neighbour, cells_count
   integer :: cells_count = 0, capacity = 0
   integer, allocatable :: ids_cache(:), owner_cache(:), index_cache(:), neighbours(:, :), buckets(:)
   logical, allocatable :: geometry_ready(:)
   real, allocatable :: centres(:, :)
contains
   integer function cell_cache_slot(id) result(slot)
      integer, intent(in) :: id
      slot = slot_of(id)
      if (slot == 0) call mesh_fail('stencil cell missing from setup cache')
   end function

   integer function slot_of(id) result(slot)
      integer, intent(in) :: id
      integer :: b, start
      slot = 0
      if (.not. allocated(buckets) .or. id <= 0) return
      b = int(modulo(int(id, int64) * 104729_int64, int(size(buckets), int64))) + 1
      start = b
      do
         slot = buckets(b)
         if (slot == 0) return
         if (ids_cache(slot) == id) return
         b = mod(b, size(buckets)) + 1
         if (b == start) exit
      end do
      slot = 0
   end function

   subroutine insert_bucket(slot)
      integer, intent(in) :: slot
      integer :: b
      b = int(modulo(int(ids_cache(slot), int64) * 104729_int64, int(size(buckets), int64))) + 1
      do while (buckets(b) /= 0)
         b = mod(b, size(buckets)) + 1
      end do
      buckets(b) = slot
   end subroutine

   subroutine reserve(count)
      integer, intent(in) :: count
      integer, allocatable :: itmp(:), matrix(:, :)
      real, allocatable :: rtmp(:, :)
      logical, allocatable :: ltmp(:)
      integer :: new_capacity, i
      integer(int64) :: wide
      if (count <= capacity) return
      wide = max(int(count, int64), 2_int64 * capacity, 16_int64)
      if (2_int64 * wide + 1 >= huge(0)) call mesh_fail('local cell cache exceeds int32')
      new_capacity = int(wide)
      allocate(itmp(new_capacity)); itmp = 0
      if (cells_count > 0) itmp(1:cells_count) = ids_cache(1:cells_count)
      call move_alloc(itmp, ids_cache)
      allocate(itmp(new_capacity)); itmp = 0
      if (cells_count > 0) itmp(1:cells_count) = owner_cache(1:cells_count)
      call move_alloc(itmp, owner_cache)
      allocate(itmp(new_capacity)); itmp = 0
      if (cells_count > 0) itmp(1:cells_count) = index_cache(1:cells_count)
      call move_alloc(itmp, index_cache)
      allocate(matrix(6, new_capacity)); matrix = 0
      if (cells_count > 0) matrix(:, 1:cells_count) = neighbours(:, 1:cells_count)
      call move_alloc(matrix, neighbours)
      allocate(rtmp(3, new_capacity)); rtmp = 0.0
      if (cells_count > 0) rtmp(:, 1:cells_count) = centres(:, 1:cells_count)
      call move_alloc(rtmp, centres)
      allocate(ltmp(new_capacity)); ltmp = .false.
      if (cells_count > 0) ltmp(1:cells_count) = geometry_ready(1:cells_count)
      call move_alloc(ltmp, geometry_ready)
      if (allocated(buckets)) deallocate(buckets)
      allocate(buckets(2 * new_capacity + 1)); buckets = 0
      capacity = new_capacity
      do i = 1, cells_count
         call insert_bucket(i)
      end do
   end subroutine

   subroutine cells_seed()
      integer :: i
      call reserve(max(1, size(ielem_ihexgl)))
      do i = 1, size(ielem_ihexgl)
         cells_count = cells_count + 1
         ids_cache(cells_count) = ielem_ihexgl(i)
         owner_cache(cells_count) = n; index_cache(cells_count) = i
         call insert_bucket(cells_count)
      end do
   end subroutine

   integer function cell_owner(id) result(value)
      integer, intent(in) :: id
      integer :: slot
      if (mesh_distributed) then
         slot = slot_of(id)
         if (slot == 0) call mesh_fail('cell owner missing from setup cache')
         value = owner_cache(slot)
      else
         value = xmpie(id)
      end if
   end function

   integer function cell_index(id) result(value)
      integer, intent(in) :: id
      integer :: slot
      if (mesh_distributed) then
         slot = slot_of(id)
         if (slot == 0) call mesh_fail('cell index missing from setup cache')
         value = index_cache(slot)
      else
         value = xmpil(id)
      end if
   end function

   real function cell_centre(id, d) result(value)
      integer, intent(in) :: id, d
      integer :: slot
      if (mesh_distributed) then
         slot = slot_of(id)
         if (slot == 0) call mesh_fail('cell centre missing from setup cache')
         if (.not. geometry_ready(slot)) call mesh_fail('cell geometry not fetched')
         value = centres(d, slot)
      else if (lowmem == 0) then
         value = xand2r(cand(xmpie(id)), xmpil(id), d)
      else
         value = centerr(id, d)
      end if
   end function

   function cell_centres(id, first, last) result(values)
      integer, intent(in) :: id, first, last
      real :: values(last - first + 1)
      integer :: d
      do d = first, last
         values(d - first + 1) = cell_centre(id, d)
      end do
   end function

   integer function cell_neighbour(id, face) result(value)
      integer, intent(in) :: id, face
      integer :: slot
      if (mesh_distributed) then
         slot = slot_of(id)
         if (slot == 0) call mesh_fail('cell neighbours missing from setup cache')
         if (.not. geometry_ready(slot)) call mesh_fail('cell neighbours not fetched')
         value = neighbours(face, slot)
      else if (lowmem == 0) then
         value = cand2r(cand(xmpie(id)), xmpil(id), face)
      else
         value = glneigh(id, face)
      end if
   end function

   subroutine cells_geometry()
      integer :: i, slot
      do i = 1, size(ielem_ihexgl)
         slot = slot_of(ielem_ihexgl(i))
         neighbours(1:ielem_ifca(i), slot) = ielem_ineighg(i, 1:ielem_ifca(i))
         centres(:, slot) = [ielem_xxc(i), ielem_yyc(i), ielem_zzc(i)]
         geometry_ready(slot) = .true.
      end do
   end subroutine

   subroutine cells_fetch(ids, with_geometry)
      integer, intent(in) :: ids(:)
      logical, intent(in) :: with_geometry
      integer, allocatable :: work(:), requested(:), owners(:), indices(:), requests(:, :), destinations(:)
      integer, allocatable :: received(:, :), replies(:, :), returned(:, :)
      real, allocatable :: values(:, :), coords(:, :)
      integer :: i, k, unique, slot, local_index
      allocate(work(size(ids))); k = 0
      do i = 1, size(ids)
         if (ids(i) <= 0) cycle
         slot = slot_of(ids(i))
         if (slot > 0) then
            if (.not. with_geometry) cycle
            if (geometry_ready(slot)) cycle
         end if
         k = k + 1; work(k) = ids(i)
      end do
      call node_sort(work, 1, k)
      unique = 0
      do i = 1, k
         if (unique > 0) then
            if (work(i) == work(unique)) cycle
         end if
         unique = unique + 1; work(unique) = work(i)
      end do
      requested = work(1:unique); deallocate(work)
      call mesh_query_cells(requested, owners, indices)
      if (int(cells_count, int64) + unique >= huge(0)) call mesh_fail('local cell cache count overflow')
      call reserve(cells_count + unique)
      do i = 1, unique
         slot = slot_of(requested(i))
         if (slot /= 0) cycle
         cells_count = cells_count + 1
         ids_cache(cells_count) = requested(i); owner_cache(cells_count) = owners(i)
         index_cache(cells_count) = indices(i)
         call insert_bucket(cells_count)
      end do
      if (.not. with_geometry) return
      allocate(requests(3, unique))
      do i = 1, unique
         requests(:, i) = [indices(i), i, n]
      end do
      call route_integer(requests, owners, received)
      allocate(replies(7, size(received, 2)), values(3, size(received, 2)), destinations(size(received, 2)))
      do i = 1, size(received, 2)
         local_index = received(1, i)
         if (local_index < 1 .or. local_index > size(ielem_ihexgl)) call mesh_fail('invalid requested cell index')
         slot = slot_of(ielem_ihexgl(local_index))
         if (.not. geometry_ready(slot)) call mesh_fail('requested geometry before neighbour setup')
         replies(1, i) = received(2, i); replies(2:7, i) = neighbours(:, slot)
         values(:, i) = centres(:, slot); destinations(i) = received(3, i)
      end do
      call route_integer(replies, destinations, returned)
      call route_real(values, destinations, coords)
      do i = 1, size(returned, 2)
         slot = slot_of(requested(returned(1, i)))
         neighbours(:, slot) = returned(2:7, i); centres(:, slot) = coords(:, i)
         geometry_ready(slot) = .true.
      end do
   end subroutine

   subroutine cells_stencils()
      integer, allocatable :: head(:), used(:), requests(:), pending(:)
      integer :: i, j, k, slot, candidate, count_needed, complete, all_complete, ierr, target, nc
      nc = size(ielem_ihexgl); target = iselemt(n)
      allocate(head(nc), used(nc), requests(nc)); head = 1; used = 1
      ilocalallelg = 0
      do i = 1, nc
         ilocalallelg(n, i, 1, 1) = ielem_ihexgl(i)
      end do
      do
         count_needed = 0; complete = 1
         do i = 1, nc
            do while (used(i) < target)
               if (head(i) > used(i)) call mesh_fail('connected component is smaller than requested stencil')
               candidate = ilocalallelg(n, i, 1, head(i)); slot = slot_of(candidate)
               if (slot == 0) then
                  count_needed = count_needed + 1; requests(count_needed) = candidate
                  exit
               end if
               if (.not. geometry_ready(slot)) then
                  count_needed = count_needed + 1; requests(count_needed) = candidate
                  exit
               end if
               do j = 1, 2 * dims
                  candidate = neighbours(j, slot)
                  if (candidate <= 0) cycle
                  if (any(ilocalallelg(n, i, 1, 1:used(i)) == candidate)) cycle
                  used(i) = used(i) + 1; ilocalallelg(n, i, 1, used(i)) = candidate
                  if (used(i) == target) exit
               end do
               head(i) = head(i) + 1
            end do
            if (used(i) < target) complete = 0
         end do
         call MPI_Allreduce(complete, all_complete, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
         if (all_complete == 1) exit
         call cells_fetch(requests(1:count_needed), .true.)
      end do
      ! Fetch terminal stencil entries too: directional selection uses their centres.
      count_needed = 0
      do i = 1, nc
         do j = 1, target
            slot = slot_of(ilocalallelg(n, i, 1, j))
            if (slot > 0) then
               if (geometry_ready(slot)) cycle
            end if
            if (count_needed == huge(0)) call mesh_fail('stencil request count overflow')
            count_needed = count_needed + 1
         end do
      end do
      allocate(pending(count_needed)); k = 0
      do i = 1, nc
         do j = 1, target
            candidate = ilocalallelg(n, i, 1, j); slot = slot_of(candidate)
            if (slot > 0) then
               if (geometry_ready(slot)) cycle
            end if
            k = k + 1; pending(k) = candidate
         end do
      end do
      call cells_fetch(pending, .true.)
   end subroutine

   subroutine cells_release()
      if (allocated(ids_cache)) deallocate(ids_cache, owner_cache, index_cache, neighbours, centres, geometry_ready, buckets)
      cells_count = 0; capacity = 0
      call mesh_owner_release()
   end subroutine
end module mesh_local_cells
