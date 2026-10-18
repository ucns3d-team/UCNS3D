! Shared-file, native-layout restart I/O. Host only; called by the MPI master thread.
! Payload records remain in global-cell order for repartition-independent restart.
module checkpoint_mpi_io
   use mpi
   use iso_c_binding, only: c_char, c_int, c_null_char
   implicit none
   private
   public :: checkpoint_file, checkpoint_begin, checkpoint_transfer, checkpoint_finish
   public :: checkpoint_header, checkpoint_fail, checkpoint_checked, checkpoint_needs_wide
   integer, parameter :: ck_wide = selected_int_kind(18)
   integer, parameter :: max_batch_cells = 262144, max_batch_values = 8388608 ! 64 MiB of real64
   type :: checkpoint_file
      integer :: fh, rank, nvalues, capacity, rounds, local_count
      integer, allocatable :: order(:), ids(:)
      integer(kind=MPI_OFFSET_KIND) :: header_bytes, file_bytes
      logical :: writing, wide
      character(len=:), allocatable :: name
   end type
   interface
      function checkpoint_rename(old, new) bind(C, name='rename') result(rc)
         import c_char, c_int
         character(kind=c_char), intent(in) :: old(*), new(*)
         integer(c_int) :: rc
      end function
   end interface
contains
   subroutine checkpoint_fail(message)
      character(len=*), intent(in) :: message
      integer :: ierr, rank
      call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
      write(*,'(A,I0,2A)') 'Checkpoint error on rank ', rank, ': ', trim(message)
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      error stop
   end subroutine

   subroutine checkpoint_checked(ierr, operation)
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: operation
      character(len=MPI_MAX_ERROR_STRING) :: message
      integer :: n, ierr2
      if (ierr == MPI_SUCCESS) return
      call MPI_Error_string(ierr, message, n, ierr2)
      call checkpoint_fail(trim(operation)//': '//message(:n))
   end subroutine

   pure logical function checkpoint_needs_wide(total, width) result(wide)
      integer(ck_wide), intent(in) :: total
      integer, intent(in) :: width
      ! Switch before multiplying; assigning an overflowing int32 product to int64 is too late.
      wide = total >= 300000000_ck_wide
      if (width > 0) wide = wide .or. total - 1_ck_wide > int(huge(0),ck_wide) / int(width,ck_wide)
   end function

   subroutine checkpoint_begin(cp, name, ids, width, total, header_bytes, writing, batch_cells)
      type(checkpoint_file), intent(out) :: cp
      character(len=*), intent(in) :: name
      integer, intent(in) :: ids(:), width
      integer(ck_wide), intent(in) :: total
      integer(kind=MPI_OFFSET_KIND), intent(in) :: header_bytes
      logical, intent(in) :: writing
      integer, optional, intent(in) :: batch_cells ! Smaller batches for regression tests; never raises the cap.
      integer :: ierr, i, maximum, mode, stat
      integer(kind=MPI_OFFSET_KIND) :: actual
      character(len=:), allocatable :: path
      if (storage_size(0) /= 32 .or. storage_size(0.0) /= 64) &
         call checkpoint_fail('native restart layout requires int32 and real64 compiler defaults')
      if (bit_size(0_MPI_ADDRESS_KIND) < 64 .or. bit_size(0_MPI_OFFSET_KIND) < 64) &
         call checkpoint_fail('MPI address and file-offset kinds must be at least 64 bits')
      if (total < 1 .or. width < 1 .or. width > max_batch_values .or. header_bytes < 0) &
         call checkpoint_fail('invalid checkpoint dimensions')
      if (total > (huge(actual) - header_bytes) / (8_ck_wide * int(width,ck_wide))) &
         call checkpoint_fail('checkpoint file size exceeds MPI_OFFSET_KIND')
      cp%name = name
      cp%writing = writing
      cp%header_bytes = header_bytes
      cp%file_bytes = header_bytes + int(total,MPI_OFFSET_KIND) * int(width,MPI_OFFSET_KIND) * 8_MPI_OFFSET_KIND
      cp%nvalues = width
      cp%local_count = size(ids)
      cp%wide = checkpoint_needs_wide(total, width)
      cp%capacity = min(max_batch_cells, max_batch_values / width)
      if (present(batch_cells)) cp%capacity = min(cp%capacity, max(1,batch_cells))
      call MPI_Comm_rank(MPI_COMM_WORLD, cp%rank, ierr)
      call checkpoint_checked(ierr, 'getting rank')
      call MPI_Allreduce(cp%local_count, maximum, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      call checkpoint_checked(ierr, 'getting maximum local record count')
      cp%rounds = 0
      if (maximum > 0) cp%rounds = (maximum - 1) / cp%capacity + 1
      allocate(cp%ids(size(ids)), cp%order(size(ids)), stat=stat)
      if (stat /= 0) call checkpoint_fail('allocating local checkpoint index')
      cp%ids = ids
      do i = 1, size(ids)
         if (ids(i) < 1 .or. int(ids(i),ck_wide) > total) call checkpoint_fail('cell ID outside global mesh')
         cp%order(i) = i
      end do
      call sort_order(cp%ids, cp%order)
      do i = 2, size(ids)
         if (ids(cp%order(i)) == ids(cp%order(i-1))) call checkpoint_fail('duplicate local checkpoint cell ID')
      end do
      path = name
      mode = MPI_MODE_RDONLY
      if (writing) then
         path = name//'.tmp'
         mode = MPI_MODE_WRONLY + MPI_MODE_CREATE
      end if
      call MPI_File_open(MPI_COMM_WORLD, path, mode, MPI_INFO_NULL, cp%fh, ierr)
      call checkpoint_checked(ierr, 'opening '//path)
      call MPI_File_set_errhandler(cp%fh, MPI_ERRORS_RETURN, ierr)
      call checkpoint_checked(ierr, 'setting checkpoint error handler')
      if (writing) then
         call MPI_File_set_size(cp%fh, 0_MPI_OFFSET_KIND, ierr)
         call checkpoint_checked(ierr, 'truncating temporary checkpoint')
         call MPI_File_set_size(cp%fh, cp%file_bytes, ierr)
         call checkpoint_checked(ierr, 'sizing temporary checkpoint')
      else
         call MPI_File_get_size(cp%fh, actual, ierr)
         call checkpoint_checked(ierr, 'checking checkpoint size')
         if (actual /= cp%file_bytes) call checkpoint_fail('size/layout mismatch in '//name// &
            ' (check mesh, previous turbulence/passive counts and steady/unsteady settings)')
      end if
      if (cp%rank == 0) then
         if (cp%wide) then
            write(*,'(3A,I0,A,I0)') 'Checkpoint ', name, ': wide MPI-I/O; bytes=', cp%file_bytes, ', batch cells=', cp%capacity
         else
            write(*,'(3A,I0)') 'Checkpoint ', name, ': native MPI-I/O; bytes=', cp%file_bytes
         end if
      end if
   end subroutine

   subroutine checkpoint_header(cp, iteration, values)
      type(checkpoint_file), intent(in) :: cp
      integer, intent(inout) :: iteration
      real, intent(inout) :: values(:)
      integer :: ierr, got, status(MPI_STATUS_SIZE)
      if (cp%header_bytes /= 4_MPI_OFFSET_KIND + 8_MPI_OFFSET_KIND * size(values)) &
         call checkpoint_fail('header dimensions do not match payload offset')
      if (cp%rank == 0) then
         if (cp%writing) then
            call MPI_File_write_at(cp%fh, 0_MPI_OFFSET_KIND, iteration, 1, MPI_INTEGER, status, ierr)
         else
            call MPI_File_read_at(cp%fh, 0_MPI_OFFSET_KIND, iteration, 1, MPI_INTEGER, status, ierr)
         end if
         call checkpoint_checked(ierr, 'transferring checkpoint iteration')
         call MPI_Get_count(status, MPI_INTEGER, got, ierr)
         call checkpoint_checked(ierr, 'checking iteration count')
         if (got /= 1) call checkpoint_fail('incomplete checkpoint iteration')
         if (cp%writing) then
            call MPI_File_write_at(cp%fh, 4_MPI_OFFSET_KIND, values, size(values), MPI_DOUBLE_PRECISION, status, ierr)
         else
            call MPI_File_read_at(cp%fh, 4_MPI_OFFSET_KIND, values, size(values), MPI_DOUBLE_PRECISION, status, ierr)
         end if
         call checkpoint_checked(ierr, 'transferring checkpoint header')
         call MPI_Get_count(status, MPI_DOUBLE_PRECISION, got, ierr)
         call checkpoint_checked(ierr, 'checking header count')
         if (got /= size(values)) call checkpoint_fail('incomplete checkpoint header')
      end if
      if (.not. cp%writing) then
         call MPI_Bcast(iteration, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call checkpoint_checked(ierr, 'broadcasting checkpoint iteration')
         call MPI_Bcast(values, size(values), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         call checkpoint_checked(ierr, 'broadcasting checkpoint header')
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      call checkpoint_checked(ierr, 'synchronizing checkpoint header')
   end subroutine

   subroutine checkpoint_transfer(cp, first, count, values)
      type(checkpoint_file), intent(in) :: cp
      integer, intent(in) :: first, count
      real, intent(inout) :: values(:)
      integer :: ierr, j, nb, gid, previous, dtype, got, nsend, stat
      integer :: status(MPI_STATUS_SIZE)
      integer, allocatable :: lengths(:), offsets(:)
      integer(kind=MPI_ADDRESS_KIND), allocatable :: byte_offsets(:)
      if (count < 0 .or. count > cp%capacity) call checkpoint_fail('invalid checkpoint batch size')
      nsend = count * cp%nvalues ! bounded by max_batch_values
      if (size(values) < nsend) call checkpoint_fail('checkpoint buffer is too small')
      dtype = MPI_DOUBLE_PRECISION
      if (count > 0) then
         allocate(lengths(count), offsets(count), byte_offsets(count), stat=stat)
         if (stat /= 0) call checkpoint_fail('allocating checkpoint file view')
         nb = 0
         previous = -1
         do j = first, first + count - 1
            gid = cp%ids(cp%order(j))
            if (gid /= previous + 1) then
               nb = nb + 1
               lengths(nb) = cp%nvalues
               byte_offsets(nb) = (int(gid,MPI_ADDRESS_KIND) - 1_MPI_ADDRESS_KIND) * &
                  int(cp%nvalues,MPI_ADDRESS_KIND) * 8_MPI_ADDRESS_KIND
               if (.not. cp%wide) offsets(nb) = int(byte_offsets(nb) / 8_MPI_ADDRESS_KIND)
            else
               lengths(nb) = lengths(nb) + cp%nvalues
            end if
            previous = gid
         end do
         if (cp%wide) then
            call MPI_Type_create_hindexed(nb, lengths, byte_offsets, MPI_DOUBLE_PRECISION, dtype, ierr)
         else
            call MPI_Type_indexed(nb, lengths, offsets, MPI_DOUBLE_PRECISION, dtype, ierr)
         end if
         call checkpoint_checked(ierr, 'constructing checkpoint file view')
         call MPI_Type_commit(dtype, ierr)
         call checkpoint_checked(ierr, 'committing checkpoint file view')
      end if
      call MPI_File_set_view(cp%fh, cp%header_bytes, MPI_DOUBLE_PRECISION, dtype, 'native', MPI_INFO_NULL, ierr)
      call checkpoint_checked(ierr, 'setting checkpoint file view')
      if (cp%writing) then
         call MPI_File_write_all(cp%fh, values, nsend, MPI_DOUBLE_PRECISION, status, ierr)
      else
         call MPI_File_read_all(cp%fh, values, nsend, MPI_DOUBLE_PRECISION, status, ierr)
      end if
      call checkpoint_checked(ierr, 'transferring checkpoint payload')
      call MPI_Get_count(status, MPI_DOUBLE_PRECISION, got, ierr)
      call checkpoint_checked(ierr, 'checking transferred payload count')
      if (got /= nsend) call checkpoint_fail('incomplete checkpoint payload')
      if (count > 0) then
         call MPI_Type_free(dtype, ierr)
         call checkpoint_checked(ierr, 'freeing checkpoint file view')
      end if
   end subroutine

   subroutine checkpoint_finish(cp)
      type(checkpoint_file), intent(inout) :: cp
      integer :: ierr, renamed
      if (cp%writing) then
         call MPI_File_sync(cp%fh, ierr)
         call checkpoint_checked(ierr, 'syncing temporary checkpoint')
      end if
      call MPI_File_close(cp%fh, ierr)
      call checkpoint_checked(ierr, 'closing checkpoint')
      if (cp%writing) then
         renamed = 0
         if (cp%rank == 0) renamed = checkpoint_rename(cp%name//'.tmp'//c_null_char, cp%name//c_null_char)
         call MPI_Bcast(renamed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call checkpoint_checked(ierr, 'broadcasting checkpoint commit')
         if (renamed /= 0) call checkpoint_fail('cannot rename completed checkpoint '//cp%name)
      end if
      deallocate(cp%order, cp%ids)
   end subroutine

   subroutine sort_order(ids, order)
      integer, intent(in) :: ids(:)
      integer, intent(inout) :: order(:)
      integer :: i, last, temp
      logical :: sorted
      sorted = .true.
      do i = 2, size(ids)
         if (ids(i) < ids(i-1)) sorted = .false.
      end do
      if (sorted) return
      ! In-place heap sort: bounded stack and no second O(local cells) workspace.
      do i = size(order)/2, 1, -1
         call sift(i, size(order))
      end do
      do last = size(order), 2, -1
         temp = order(1); order(1) = order(last); order(last) = temp
         call sift(1, last-1)
      end do
   contains
      subroutine sift(start, finish)
         integer, intent(in) :: start, finish
         integer :: root, child, item
         root = start
         item = order(root)
         do while (root <= finish/2)
            child = 2*root
            if (child < finish) then
               if (ids(order(child)) < ids(order(child+1))) child = child+1
            end if
            if (ids(item) >= ids(order(child))) exit
            order(root) = order(child)
            root = child
         end do
         order(root) = item
      end subroutine
   end subroutine
end module checkpoint_mpi_io
