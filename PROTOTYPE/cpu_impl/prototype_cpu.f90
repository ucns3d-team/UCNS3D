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

program prototype

#ifdef CIRRUS_GNU_BUILD
  implicit none
  include "mpif.h"
#else
  use mpi
  implicit none
#endif

integer :: rank       ! rank of each process  ! rank = N
integer :: tot_ranks  ! total number of ranks ! tot_ranks = isize
integer :: ierror
integer :: status(mpi_status_size) ! unused

integer :: provided, alpha, beta
integer :: i, j, k, s, m, num_dof, tot_stencils, num_stencils ! s = l, num_stencils = ll unused, num_dof = idegfree, tot_stencils = stencil_local2
integer :: num_elems, num_vars, num_neighbours ! num_elems = kmaxe, num_vars = nof_variables, imax = num_neighbours, kmaxn unused, dof unused
integer, save :: stencil_local  ! threadprivate requires save attribute under GNU 
real :: r_s, r_num_elems, r_j, r_k
real, allocatable, save, dimension(:) :: leftv ! threadprivate requires SAVE attribute under GNU

type local_recon3
  integer :: local, mrf, g0
  integer, allocatable, dimension(:,:) :: ihexg ! global index of cells
  real, allocatable, dimension(:) :: cond       ! dummy variable used for gradient approximation estimation 
  real, allocatable, dimension(:,:,:) :: invmat, sol, matrix_1 ! (var, dim, i_face)
end type local_recon3
type(local_recon3), allocatable, dimension(:) :: ilocal_recon3


!$omp threadprivate(leftv, stencil_local)

call mpi_init_thread(mpi_thread_funneled, provided, ierror)
call mpi_comm_size(mpi_comm_world, tot_ranks, ierror)
call mpi_comm_rank(mpi_comm_world, rank, ierror)

#ifdef UCNS3D_DEBUG
#ifdef CIRRUS_GNU_BUILD
  if (rank .eq. 0) then
    write(*, *) "COMPILING UCNS3D WITH GNU ON CIRRUS"
  end if
#endif
#endif

! for all mpi processes
num_elems = 100
num_vars = 5
num_dof = 10
num_neighbours = 20
num_stencils = 4  ! unused?


! now allocate a thread private variable
!$omp parallel default(shared)
allocate(leftv(1:num_vars))
!$omp end parallel 
!$omp barrier

allocate(ilocal_recon3(1:num_elems))

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

  do j=1,num_elems
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


!$omp parallel default(shared)
!$omp do
do i=1,num_elems
  if (i .lt. 500) then
    tot_stencils = 5
  else
    tot_stencils = 4
  end if

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
    ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) = ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) + i * 2
  end do
end do
!$omp end do
!$omp end parallel 


if (rank .eq. 0) then
!$omp master
  open(20, file='OUTPUT_RESULT.dat', form='formatted', status='new', action='write')
  do i=1,num_elems
    write(20, *) ilocal_recon3(i)%sol(:,:,:)
  end do
  close(20)
!$omp end master
!$omp barrier
end if


call mpi_barrier(mpi_comm_world, ierror)
call mpi_finalize(ierror)
    
end program prototype
