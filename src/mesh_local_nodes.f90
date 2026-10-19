! Compact host-side vertex storage for nonperiodic MPI mesh startup.
! Connectivity keeps global IDs until the existing final local renumbering.
module mesh_local_nodes
   use declaration, only: dinoder, dinoder2, dimensiona, iperiodicity, imaxn
   use mpiinfo, only: n
   use mesh_mpi_io
   use iso_fortran_env, only: int64
   implicit none
   private
   public :: node_final_index, node_output_owner
   public :: local_nodes_enabled, node_ids, node_xyz, node_owners, node_slot, node_present
   public :: node_coord, node_coords, node_degree, node_cell, node_remote_count, node_remote_cell
   public :: local_nodes_init, local_nodes_adjacency, local_nodes_remote, local_nodes_release
   public :: local_nodes_release_adjacency, node_sort
   logical :: local_nodes_enabled = .false.
   integer, allocatable :: node_ids(:), node_owners(:), buckets(:)
   integer, allocatable :: adjacent_ptr(:), adjacent_cells(:), remote_ptr(:), remote_cells(:)
   real, allocatable :: node_xyz(:, :)
contains
   recursive subroutine node_sort(a, left, right)
      integer, intent(inout) :: a(:)
      integer, intent(in) :: left, right
      integer :: i, j, pivot, swap
      if (left >= right) return
      pivot = a(left + (right - left) / 2); i = left; j = right
      do
         do while (a(i) < pivot)
            i = i + 1
         end do
         do while (a(j) > pivot)
            j = j - 1
         end do
         if (i > j) exit
         swap = a(i); a(i) = a(j); a(j) = swap
         i = i + 1; j = j - 1
         if (i > j) exit
      end do
      if (left < j) call node_sort(a, left, j)
      if (i < right) call node_sort(a, i, right)
   end subroutine

   integer function node_slot(id) result(slot)
      integer, intent(in) :: id
      integer :: b, start
      slot = 0
      if (id <= 0) return
      if (.not. local_nodes_enabled) then
         if (id <= size(dinoder)) then
            if (dinoder(id)%itor > 0) slot = id
         end if
         return
      end if
      b = int(modulo(int(id, int64) * 104729_int64, int(size(buckets), int64))) + 1
      start = b
      do
         slot = buckets(b)
         if (slot == 0) return
         if (node_ids(slot) == id) return
         b = mod(b, size(buckets)) + 1
         if (b == start) exit
      end do
      slot = 0
   end function

   integer function node_final_index(id) result(value)
      integer, intent(in) :: id
      if (local_nodes_enabled) then
         value = node_slot(id)
         if (value == 0) call mesh_fail('missing vertex at final local renumbering')
      else
         value = dinoder(id)%itor
      end if
   end function

   integer function node_output_owner(id, legacy_owners) result(value)
      integer, intent(in) :: id, legacy_owners(:)
      if (local_nodes_enabled) then
         value = node_owners(node_slot(id))
      else
         value = legacy_owners(id)
      end if
   end function

   logical function node_present(id)
      integer, intent(in) :: id
      node_present = node_slot(id) > 0
   end function

   real function node_coord(id, d) result(value)
      integer, intent(in) :: id, d
      integer :: slot
      if (local_nodes_enabled) then
         slot = node_slot(id)
         if (slot == 0) call mesh_fail('coordinate requested for an absent local vertex')
         value = node_xyz(d, slot)
      else
         value = dinoder(id)%cord(d)
      end if
   end function

   function node_coords(id) result(values)
      integer, intent(in) :: id
      real :: values(dimensiona)
      integer :: slot
      if (local_nodes_enabled) then
         slot = node_slot(id)
         if (slot == 0) call mesh_fail('coordinates requested for an absent local vertex')
         values = node_xyz(1:dimensiona, slot)
      else
         values = dinoder(id)%cord(1:dimensiona)
      end if
   end function

   subroutine local_nodes_init(connectivity, counts, scale)
      integer, intent(in) :: connectivity(:, :), counts(:)
      real, intent(in) :: scale
      integer, allocatable :: ids(:)
      integer :: i, j, k, total, unique, b
      integer(int64) :: wide
      wide = sum(int(counts, int64))
      if (wide >= huge(0)) call mesh_fail('local vertex references exceed int32')
      total = int(wide); allocate(ids(total)); k = 0
      do i = 1, size(counts)
         do j = 1, counts(i)
            k = k + 1; ids(k) = connectivity(j, i)
         end do
      end do
      call node_sort(ids, 1, total)
      unique = 0
      do i = 1, total
         if (unique > 0) then
            if (ids(i) == ids(unique)) cycle
         end if
         unique = unique + 1; ids(unique) = ids(i)
      end do
      allocate(node_ids(unique)); node_ids = ids(1:unique)
      deallocate(ids)
      wide = 2_int64 * unique + 1
      if (wide >= huge(0)) call mesh_fail('local vertex lookup exceeds int32')
      allocate(buckets(int(wide))); buckets = 0
      do i = 1, unique
         b = int(modulo(int(node_ids(i), int64) * 104729_int64, wide)) + 1
         do while (buckets(b) /= 0)
            b = mod(b, size(buckets)) + 1
         end do
         buckets(b) = i
      end do
      local_nodes_enabled = .true.
      call mesh_vertices(node_ids, node_xyz)
      node_xyz = node_xyz / scale
      call mesh_node_owners(node_ids, node_owners)
      ! Empty legacy descriptors keep existing argument/deallocation contracts.
      allocate(dinoder(0), dinoder2(0))
   end subroutine

   subroutine local_nodes_adjacency(connectivity, counts)
      integer, intent(in) :: connectivity(:, :), counts(:)
      integer, allocatable :: cursor(:)
      integer :: i, j, slot
      allocate(adjacent_ptr(size(node_ids) + 1)); adjacent_ptr = 0
      do i = 1, size(counts)
         do j = 1, counts(i)
            slot = node_slot(connectivity(j, i))
            if (slot == 0) call mesh_fail('missing vertex in local adjacency')
            adjacent_ptr(slot + 1) = adjacent_ptr(slot + 1) + 1
         end do
      end do
      adjacent_ptr(1) = 1
      do i = 1, size(node_ids)
         adjacent_ptr(i + 1) = adjacent_ptr(i + 1) + adjacent_ptr(i)
      end do
      allocate(adjacent_cells(adjacent_ptr(size(node_ids) + 1) - 1))
      cursor = adjacent_ptr(1:size(node_ids))
      do i = 1, size(counts)
         do j = 1, counts(i)
            slot = node_slot(connectivity(j, i))
            adjacent_cells(cursor(slot)) = i; cursor(slot) = cursor(slot) + 1
         end do
      end do
   end subroutine

   subroutine local_nodes_remote(owners)
      integer, intent(in) :: owners(:)
      integer, allocatable :: ptr(:), cells(:), cell_owners(:), cell_indices(:)
      integer :: i, j, k
      call mesh_incidence(node_ids, ptr, cells)
      if (mesh_distributed) then
         call mesh_query_cells(cells, cell_owners, cell_indices)
      else
         cell_owners = owners(cells)
      end if
      allocate(remote_ptr(size(node_ids) + 1)); remote_ptr = 1
      k = 0
      do i = 1, size(node_ids)
         do j = ptr(i), ptr(i + 1) - 1
            if (cell_owners(j) /= n) k = k + 1
         end do
         remote_ptr(i + 1) = k + 1
      end do
      allocate(remote_cells(k)); k = 0
      do i = 1, size(node_ids)
         do j = ptr(i), ptr(i + 1) - 1
            if (cell_owners(j) == n) cycle
            k = k + 1; remote_cells(k) = cells(j)
         end do
      end do
   end subroutine

   integer function node_degree(id) result(value)
      integer, intent(in) :: id
      integer :: slot
      value = 0
      if (local_nodes_enabled) then
         slot = node_slot(id)
         if (slot > 0) value = adjacent_ptr(slot + 1) - adjacent_ptr(slot)
      else
         value = dinoder2(id)%numberofneib
      end if
   end function

   integer function node_cell(id, position) result(value)
      integer, intent(in) :: id, position
      if (local_nodes_enabled) then
         value = adjacent_cells(adjacent_ptr(node_slot(id)) + position - 1)
      else
         value = dinoder2(id)%neibids(position)
      end if
   end function

   integer function node_remote_count(id) result(value)
      integer, intent(in) :: id
      integer :: slot
      value = 0
      if (local_nodes_enabled) then
         slot = node_slot(id)
         if (slot > 0) value = remote_ptr(slot + 1) - remote_ptr(slot)
      else
         value = dinoder2(id)%xne(1)
      end if
   end function

   integer function node_remote_cell(id, position) result(value)
      integer, intent(in) :: id, position
      if (local_nodes_enabled) then
         value = remote_cells(remote_ptr(node_slot(id)) + position - 1)
      else
         value = dinoder2(id)%xneib(position)
      end if
   end function

   subroutine local_nodes_release_adjacency()
      if (allocated(adjacent_ptr)) deallocate(adjacent_ptr, adjacent_cells)
      if (allocated(remote_ptr)) deallocate(remote_ptr, remote_cells)
   end subroutine

   subroutine local_nodes_release()
      call local_nodes_release_adjacency()
      if (allocated(node_ids)) deallocate(node_ids, node_xyz, node_owners, buckets)
      local_nodes_enabled = .false.
   end subroutine
end module mesh_local_nodes
