!
! UCNS3D : An Open Source, High-Order, Finite-Volume CFD Solver
!
! Main contact:  Panagiotis (Takis) Tsoutsanis (Cranfield University)
!
! Licensed under the GNU General Public License, Version 3.0 (the "License");
! you may not use this file except in compliance with the License.
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "as is" basis,
! without warranties or conditions of any kind, either express or implied.
! see the License for the specific language governing permissions and
! limitations under the License.
!

  module benchclock

  implicit none

  logical,          save, private :: firstcall = .true.
  double precision, save, private :: ticktime = 0.0

  integer, parameter :: int32kind = selected_int_kind( 9)
  integer, parameter :: int64kind = selected_int_kind(18)

!
!  Select high resolution clock
!

  integer, parameter :: intkind = int64kind
  integer(kind = intkind) :: clkcount, clkrate

contains

double precision function benchtime()

  double precision :: dummy

! Ensure clock is initialised  

  if (firstcall) dummy = benchtick()

  call system_clock(clkcount)

  benchtime  = dble(clkcount)*ticktime

end function benchtime


double precision function benchtick()

  if (firstcall) then

     firstcall = .false.
     call system_clock(clkcount, clkrate)
     ticktime = 1.0d0/dble(clkrate)

  end if

  benchtick = ticktime

end function benchtick

end module benchclock

program prototype

  use benchclock

  implicit none

  integer :: provided, ierror

  double precision :: alpha, beta, check

  ! !!! s = l, num_stencils = ll unused, num_dof = idegfree, tot_stencils = stencil_local2
  integer :: i, j, k, s, m, num_dof, tot_stencils, num_stencils

  ! !!! num_elems = kmaxe, num_vars = nof_variables, imax = num_neighbours, kmaxn unused, dof unused
  integer :: num_elems, num_vars, num_neighbours

  ! !!! r_s = r_l, r_num_elems = r_kmaxe
  double precision :: r_s, r_num_elems, r_j, r_k

  ! !!! threadprivate requires save attribute under GNU - ?
  integer, save :: stencil_local
  double precision, allocatable, save, dimension(:) :: leftv  ! threadprivate requires SAVE attribute under GNU

  ! command-line argument for convenience
  integer :: istat
  character (len=10) :: input_num_elems

  ! output filename
  character (len=32) :: output_filename
  character(len=9) :: out_prefix
  character(len=4) :: out_suffix

  ! timing parameters
  double precision :: time_start, time_end

  type local_recon3
    integer :: local, mrf, g0
    integer, allocatable, dimension(:,:) :: ihexg  ! global index of cells
    double precision, allocatable, dimension(:) :: cond        ! dummy variable used for gradient approximation estimation
    double precision, allocatable, dimension(:,:,:) :: invmat, sol, matrix_1  ! (var, dim, i_face)
  end type local_recon3
  type(local_recon3), allocatable, dimension(:) :: ilocal_recon3

  call get_command_argument(1, input_num_elems)
  read (input_num_elems, '(i7)', iostat=istat) num_elems
  num_vars = 5
  num_dof = 10
  num_neighbours = 20
  num_stencils = 4  ! unused?

  write (*, '(a35, i8)') 'NUMBER OF ELEMENTS : ', num_elems
  write (*, '(a35, i8)') 'NUMBER OF VARIABLES : ', num_vars
  write (*, '(a35, i8)') 'NUMBER OF DEGREES OF FREEDOM : ', num_dof
  write (*, '(a35, i8)') 'NUMBER OF NEIGHBOURS : ', num_neighbours
  write (*, '(a35, i8)') 'NUMBER OF STENCILS : ', num_stencils

  allocate(ilocal_recon3(1:num_elems))

  ! TIMEIT
  time_start = benchtime()

  do i=1,num_elems

    if (i .lt. 500) then
      stencil_local = 5
    else
      stencil_local = 4
    end if

    allocate(ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, 1:stencil_local));
    allocate(ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, 1:stencil_local));
    allocate(ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, 1:stencil_local));

    do j=1,num_dof
      do k=1,num_neighbours
        do s=1,stencil_local
          r_j = j
          r_k = k
          r_s = s
          ilocal_recon3(i)%invmat(j, k, s) = atan(r_j / r_k) * r_s
        end do
      end do
    end do

    do j=1,num_neighbours
      do k=1,num_vars
        do s=1,stencil_local
          r_j = j
          r_k = k
          r_s = s
          ilocal_recon3(i)%matrix_1(j, k, s) = atan(r_j / r_k) * r_s
        end do
      end do
    end do
  end do

  do i=1,num_elems
     ilocal_recon3(i)%sol(:,:,:) = 0.0
  end do

  time_start = benchtime()

  do i=1,num_elems
    if (i .lt. 500) then
      tot_stencils = 5
    else
      tot_stencils = 4
    end if

    
    ! set DGEMM parameters alpha and beta here
    alpha = 1.0
    beta = 0.0

    do s=1,tot_stencils

      ! DGEMM calculates: C := alpha * A x B + beta * C
      ! num_dof  : number of rows of A and C
      ! num_vars : number of columns of B and C
      ! num_neighbours : number of column of A and number of rows of B
      ! C is the solution matrix
      call dgemm('n', 'n', num_dof, num_vars, num_neighbours, alpha, &
                 ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, s), num_dof,&
                 ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, s), num_neighbours, &
                 beta, ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s), num_dof)

      ! TODO:NOTE this is to fill solution matrix with non-zero entries 
      !ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) = ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) + i * 2
    end do
  end do

  i = 1

  !write(*,*) 'invmat'
!  write(*,*) '------'
!    write(*,*) ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, 1:stencil_local)
!  write(*,*) '------'
!  write(*,*)    
!    write(*,*) 'matrix_1'
!  write(*,*) '------'
!    write(*,*) ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, 1:stencil_local)
!  write(*,*) '------'
!  write(*,*)    
!    write(*,*) 'sol'
!  write(*,*) '------'
!    write(*,*) ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, 1:stencil_local)
!  write(*,*) '------'


  check = 0.0
  ! verification
  do i=1,num_elems
     check = check +  sum(ilocal_recon3(i)%sol(:,:,:)**2)
  end do
  ! TIMEIT
  time_end = benchtime() - time_start
    write (*, '(a35, es15.7)') "TIME: DGEMM (s): ", time_end
    write (*, '(a35, es15.7)') "check val      : ", check

end program prototype
