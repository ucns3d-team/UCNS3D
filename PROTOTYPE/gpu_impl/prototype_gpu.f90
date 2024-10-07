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

!TODO:CHECK keeping module in same file is important?

! reference implementation of LAPACK's DGEMM
module lapck

  implicit none

  contains

  subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
    !$omp declare target
    !
    ! -- reference blas level3 routine --
    ! -- reference blas is a software package provided by univ. of tennessee,    --
    ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
    !
    ! .. scalar arguments ..
      double precision alpha,beta
      integer k,lda,ldb,ldc,m,n
      logical transa,transb
    ! ..
    ! .. array arguments ..
      double precision a(lda,*),b(ldb,*),c(ldc,*)
    ! ..
    !
    ! =====================================================================
    !
    ! ..
    ! ..
    ! .. intrinsic functions ..
      intrinsic max
    ! ..
    ! .. local scalars ..
      double precision temp
      integer i,info,j,l,nrowa,nrowb
      logical nota,notb
    ! ..
    ! .. parameters ..
      double precision one,zero
      parameter(one=1.0d+0,zero=0.0d+0)
    ! ..
    !
    ! set  nota  and  notb  as  true if  a  and  b  respectively are not
    ! transposed and set  nrowa and nrowb  as the number of rows of  a
    ! and  b  respectively.
    !
      nota = .not. transa
      notb = .not. transb
      if (nota) then
        nrowa = m
      else
        nrowa = k
      end if
      if (notb) then
        nrowb = k
      else
        nrowb = n
      end if

    !
    ! quick return if possible.
    !
      if ((m .eq. 0) .or. (n .eq. 0) .or. (((alpha .eq. zero) .or. (k .eq. 0)) .and. (beta .eq. one))) return
    !
    ! and if  alpha.eq.zero.
    !
      if (alpha.eq.zero) then
        if (beta.eq.zero) then
          do 20 j = 1,n
            do 10 i = 1,m
              c(i,j) = zero
    10      continue
    20    continue
        else
          do 40 j = 1,n
            do 30 i = 1,m
              c(i,j) = beta*c(i,j)
    30      continue
    40    continue
        end if
        return
      end if
    !
    ! start the operations.
    !
      if (notb) then
        if (nota) then
    !
    !     form  c := alpha*a*b + beta*c.
    !
          do 90 j = 1,n
            if (beta.eq.zero) then
              do 50 i = 1,m
                c(i,j) = zero
    50        continue
            else if (beta.ne.one) then
              do 60 i = 1,m
                c(i,j) = beta*c(i,j)
    60        continue
            end if
            do 80 l = 1,k
              temp = alpha*b(l,j)
              do 70 i = 1,m
                c(i,j) = c(i,j) + temp*a(i,l)
    70        continue
    80      continue
    90    continue
        else
    !
    !     form  c := alpha*a**t*b + beta*c
    !
          do 120 j = 1,n
            do 110 i = 1,m
              temp = zero
              do 100 l = 1,k
                temp = temp + a(l,i)*b(l,j)
    100       continue
              if (beta.eq.zero) then
                c(i,j) = alpha*temp
              else
                c(i,j) = alpha*temp + beta*c(i,j)
              end if
    110     continue
    120   continue
        end if
      else if (nota) then
    !
    !     form  c := alpha*a*b**t + beta*c
    !
          do 170 j = 1,n
            if (beta.eq.zero) then
              do 130 i = 1,m
                 c(i,j) = zero
    130       continue
            else if (beta.ne.one) then
              do 140 i = 1,m
                c(i,j) = beta*c(i,j)
    140       continue
            end if
            do 160 l = 1,k
              temp = alpha*b(j,l)
              do 150 i = 1,m
                c(i,j) = c(i,j) + temp*a(i,l)
    150       continue
    160     continue
    170   continue
        else
   !
   !      form  c := alpha*a**t*b**t + beta*c
   !
          do 200 j = 1,n
            do 190 i = 1,m
              temp = zero
                do 180 l = 1,k
                  temp = temp + a(l,i)*b(j,l)
    180         continue
                if (beta.eq.zero) then
                  c(i,j) = alpha*temp
                else
                  c(i,j) = alpha*temp + beta*c(i,j)
                end if
    190     continue
    200   continue
        end if
      end if
    !
      return
    !
    ! end of dgemm
    !
  end subroutine dgemm

end module lapck

program prototype

  use omp_lib
  use MPI

  use lapck

  implicit none

  integer :: rank       ! rank of each process  ! rank = N
  integer :: tot_ranks  ! total number of ranks ! tot_ranks = isize
  integer :: provided, ierror
  integer :: status(MPI_STATUS_SIZE) ! unused

  real :: alpha, beta
  logical :: transa = .false.

  ! !!! s = l, num_stencils = ll unused, num_dof = idegfree, tot_stencils = stencil_local2
  integer :: i, j, k, s, m, num_dof, tot_stencils, num_stencils

  ! !!! num_elems = kmaxe, num_vars = nof_variables, imax = num_neighbours, kmaxn unused, dof unused
  integer :: num_elems, num_vars, num_neighbours

  ! !!! r_s = r_l, r_num_elems = r_kmaxe
  real :: r_s, r_num_elems, r_j, r_k

  ! !!! threadprivate requires save attribute under GNU - ?
  integer, save :: stencil_local
  real, allocatable, save, dimension(:) :: leftv  ! threadprivate requires SAVE attribute under GNU

  ! command-line argument for convenience
  integer :: istat
  character (len=10) :: input_num_elems

  ! output filename
  character (len=32) :: output_filename
  character(len=9) :: out_prefix
  character(len=4) :: out_suffix

  ! timing parameters
  real :: time_start, time_end


  !integer :: i, j, k, l, m, idegfree, stencil_local2, ll
  !integer::kmaxe,kmaxn,nof_variables,dof,imax,stencil_local
  !real::r_l,r_kmaxe,r_j,r_k
  !real,allocatable,dimension(:)::leftv

  !nbedit
  integer :: num_devices = 0
  integer :: device_num = 0
  !nbedit

  type local_recon3
    integer :: local, mrf, g0
    real, allocatable, dimension(:) :: cond       ! dummy variable used for gradient approximation estimation 
    integer, allocatable, dimension(:,:) :: ihexg ! global index of cells
    real, allocatable, dimension(:,:,:) :: invmat, sol, matrix_1  ! (var, dim, i_face)
  end type local_recon3
  type(local_recon3), allocatable, dimension(:) :: ilocal_recon3
 
!$omp threadprivate(leftv,stencil_local)

  call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierror)
  call MPI_Comm_size(MPI_COMM_WORLD, tot_ranks, ierror)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror)

  ! for all mpi processes
  call get_command_argument(1, input_num_elems)
  read (input_num_elems, '(i7)', iostat=istat) num_elems
  !num_elems =
  num_vars = 5
  num_dof = 10
  num_neighbours = 20
  num_stencils = 4  ! unused?

  write (*, '(a35, i8)') 'NUMBER OF ELEMENTS : ', num_elems
  write (*, '(a35, i8)') 'NUMBER OF VARIABLES : ', num_vars
  write (*, '(a35, i8)') 'NUMBER OF DEGREES OF FREEDOM : ', num_dof
  write (*, '(a35, i8)') 'NUMBER OF NEIGHBOURS : ', num_neighbours
  write (*, '(a35, i8)') 'NUMBER OF STENCILS : ', num_stencils

  ! for all mpi processes
  !kmaxe = 10000
  !nof_variables = 5
  !idegfree = 10
  !imax = 20
  !ll = 4

! now allocate a thread private variable
!$omp parallel default(shared)
  allocate(leftv(1:num_vars))
!$omp end parallel 
!$omp barrier

  allocate(ilocal_recon3(1:num_elems))

  ! TIMEIT
  time_start = MPI_Wtime()

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

  ! TIMEIT
  time_end = MPI_Wtime() - time_start
  if (rank .eq. 0) then
    write(*, '(a35, es15.7)') "TIME: FILL ilocal_recon3 (s): ", time_end
  end if

  ! check some values of matrix_1 in ilocal_recon3
  do j=1,5
    do k=1,1
      do s=1,1
        i = 1
        write (*, '(a, i2, a, i2, a, i2, a, i2, a, f12.10)') 'I=', i, 'J=', j, ',K=', k,', S=', s, ',ILOCAL=', &
                                                             ilocal_recon3(i)%matrix_1(j,k,s)
      end do
    end do
  end do
  
!!$omp parallel default(shared)
!!$omp do

!nbedit
!$omp target data map(from: device_num)
!nbedit

!nbedit
  device_num = omp_get_device_num()
  num_devices = omp_get_num_devices()
!nbedit

  ! TIMEIT
  time_start = MPI_Wtime()

!$omp target teams distribute parallel do
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
      call dgemm(transa, transa, num_dof, num_vars, num_neighbours, alpha, &
                 ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, s), num_dof,&
                 ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, s), num_neighbours, &
                 beta, ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s), num_dof)

      ! TODO:CHECK where did the following line come from?
      !ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) = ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) + i * 2
    end do
  end do
!$omp end target teams distribute parallel do

!!$omp end do
!!$omp end parallel 

!nbedit
!$omp end target data
  write(*, *) 'got number of devices : ', num_devices
  write(*, *) 'got device number : ', device_num
!nbedit

  ! TIMEIT
  time_end = MPI_Wtime() - time_start
  if (rank .eq. 0) then
    write (*, '(a35, es15.7)') "TIME: DGEMM (s): ", time_end
  end if

  ! TIMEIT
  time_start = MPI_Wtime()

  if (rank .eq. 0) then
!$omp master
    out_prefix = 'OUTPUT_Ne'
    out_suffix = '.DAT'
    output_filename =  out_prefix//trim(input_num_elems)//out_suffix
    open(20, file=output_filename, form='formatted', status='new', action='write')
    do i=1,num_elems
      write(20,*)ilocal_recon3(i)%sol(:,:,:)
    end do
    close(20)
!$omp end master
!$omp barrier
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierror)
  call MPI_Finalize(ierror)
    
end program prototype
