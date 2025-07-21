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
module mylapack

  implicit none

  contains

    subroutine mydgemmgpu(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
      !$omp declare target
      !
      ! -- reference blas level3 routine --
      ! -- reference blas is a software package provided by univ. of tennessee,    --
      ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
      !
      ! .. scalar arguments ..
      double precision alpha, beta
      integer k, lda, ldb, ldc, m, n
      logical transa, transb
      ! ..
      ! .. array arguments ..
      double precision a(lda,*), b(ldb,*), c(ldc,*)
      ! ..
      !
      ! =====================================================================
      ! ..
      ! .. local scalars ..
      double precision temp
      integer i, j, l
      ! ..
      ! .. parameters ..
      ! ..

      !
      !     form  c := alpha*a*b
      !
      !$omp parallel do simd private(j,l)
      do i = 1,m

         do j = 1,n
            c(i,j) = 0.0
         end do

         do l = 1,k
            do j = 1,n
               c(i,j) = c(i,j) + a(i,l)*b(l,j)
            end do
         end do

      end do

      return
      !
      ! end of dgemm
      !
    end subroutine mydgemmgpu

  subroutine mydgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c,ldc)
    !
    ! -- reference blas level3 routine --
    ! -- reference blas is a software package provided by univ. of tennessee,    --
    ! -- univ. of california berkeley, univ. of colorado denver and nag ltd..--
    !
    ! .. scalar arguments ..
      double precision alpha, beta
      integer k, lda, ldb, ldc, m, n
      logical transa, transb
    ! ..
    ! .. array arguments ..
      double precision a(lda,*), b(ldb,*), c(ldc,*)
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
      integer i, info, j, l, nrowa, nrowb
      logical nota, notb
    ! ..
    ! .. parameters ..
      double precision one, zero
      parameter(one=1.0d+0, zero=0.0d+0)
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
      if ((m .eq. 0) .or. (n .eq. 0) .or. (((alpha .eq. zero) .or. (k .eq. 0)) .and. (beta .eq. one))) then
        return
      end if
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
   !    form  c := alpha*a**t*b**t + beta*c
   !
        do 200 j = 1,n
          do 190 i = 1,m
            temp = zero
              do 180 l = 1,k
                temp = temp + a(l,i)*b(j,l)
    180       continue
              if (beta.eq.zero) then
                c(i,j) = alpha*temp
              else
                c(i,j) = alpha*temp + beta*c(i,j)
              end if
    190   continue
    200 continue
      end if
    !
      return
    !
    ! end of dgemm
    !
  end subroutine

end module mylapack

program prototype

  use benchclock
  use mylapack

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
  double precision :: time_start1, time_start2, time_end1, time_end2
  double precision :: nflop

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
  write (*, '(a35, i8)') 'NUMBER OF STENCILS : 5 then 4'

  allocate(ilocal_recon3(1:num_elems))

  nflop = dble(min(num_elems, 500-1))*5.0*dble(num_dof)*dble(num_neighbours)*dble(num_vars)

  if (num_elems >= 500) then
     nflop = nflop + &
          dble((num_elems-500+1))*4.0*dble(num_dof)*dble(num_neighbours)*dble(num_vars)
  end if

  write(*,*) "nflop = ", nflop

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


  time_start1 = benchtime()

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
        ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) &
             = matmul(ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, s),&
             ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, s))
     end do
  end do


  time_end1 = benchtime() - time_start1

  check = 0.0
  ! verification
  do i=1,num_elems
     check = check +  sum(ilocal_recon3(i)%sol(:,:,:)**2)
  end do

  write (*, '(a35, es15.7)') "TIME: MATMUL(s): ", time_end1
  write (*, '(a35, es15.7)') "GFLOP:MATMUL   : ", 1.0e-9*(nflop/time_end1)
  write (*, '(a35, es15.7)') "check val      : ", check

  do i=1,num_elems
     ilocal_recon3(i)%sol(:,:,:) = 0.0
  end do


  time_start1 = benchtime()

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
        call mydgemm(.false.,.false., num_dof, num_vars, num_neighbours, alpha, &
             ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, s), num_dof,&
             ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, s), num_neighbours, &
             beta, ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s), num_dof)
     end do
  end do


  time_end1 = benchtime() - time_start1

  check = 0.0
  ! verification
  do i=1,num_elems
     check = check +  sum(ilocal_recon3(i)%sol(:,:,:)**2)
  end do

  write (*, '(a35, es15.7)') "TIME: MYDGEMM(s): ", time_end1
  write (*, '(a35, es15.7)') "GFLOP:MYDGEMM   : ", 1.0e-9*(nflop/time_end1)
  write (*, '(a35, es15.7)') "check val       : ", check

  do i=1,num_elems
     ilocal_recon3(i)%sol(:,:,:) = 0.0
  end do

  time_start1 = benchtime()

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

     end do
  end do


  time_end1 = benchtime() - time_start1

  check = 0.0
  ! verification
  do i=1,num_elems
     check = check +  sum(ilocal_recon3(i)%sol(:,:,:)**2)
  end do
  ! TIMEIT
  write (*, '(a35, es15.7)') "TIME: DGEMM  (s): ", time_end1
  write (*, '(a35, es15.7)') "GFLOP:DGEMM     : ", 1.0e-9*(nflop/time_end1)
  write (*, '(a35, es15.7)') "check val      : ", check

  do i=1,num_elems
     ilocal_recon3(i)%sol(:,:,:) = 0.0
  end do

  time_start1 = benchtime()

  !$omp target data map(tofrom:ilocal_recon3)

  time_start2 = benchtime()

!  !$omp target teams distribute parallel do simd private(s, tot_stencils)
  !$omp target teams distribute private(s, tot_stencils)

  do i=1,num_elems
     if (i .lt. 500) then
        tot_stencils = 5
     else
        tot_stencils = 4
     end if


!     !$omp parallel do simd
     do s=1,tot_stencils

        ! DGEMM calculates: C := alpha * A x B + beta * C
        ! num_dof  : number of rows of A and C
        ! num_vars : number of columns of B and C
        ! num_neighbours : number of column of A and number of rows of B
        ! C is the solution matrix
        ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s) &
             = matmul(ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, s),&
             ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, s))

     end do
  end do

  time_end2 = benchtime() - time_start2

  !$omp end target data

  time_end1 = benchtime() - time_start1

  check = 0.0
  ! verification
  do i=1,num_elems
     check = check +  sum(ilocal_recon3(i)%sol(:,:,:)**2)
  end do
  ! TIMEIT
  write (*, '(a35, es15.7)') "TIME: tot (s): ", time_end1
  write (*, '(a35, es15.7)') "TIME: GPU (s): ", time_end2
  write (*, '(a35, es15.7)') "GFLOP:GMATMUL: ", 1.0e-9*(nflop/time_end2)
  write (*, '(a35, es15.7)') "check val    : ", check

  do i=1,num_elems
     ilocal_recon3(i)%sol(:,:,:) = 0.0
  end do

  time_start1 = benchtime()

  !$omp target data map(tofrom:ilocal_recon3)

  time_start2 = benchtime()

!  !$omp target teams distribute parallel do simd private(s, tot_stencils)
  !$omp target teams distribute private(s, tot_stencils)

  do i=1,num_elems
     if (i .lt. 500) then
        tot_stencils = 5
     else
        tot_stencils = 4
     end if


!     !$omp parallel do simd
     do s=1,tot_stencils

        ! DGEMM calculates: C := alpha * A x B + beta * C
        ! num_dof  : number of rows of A and C
        ! num_vars : number of columns of B and C
        ! num_neighbours : number of column of A and number of rows of B
        ! C is the solution matrix
        call mydgemmgpu(.false.,.false., num_dof, num_vars, num_neighbours, alpha, &
             ilocal_recon3(i)%invmat(1:num_dof, 1:num_neighbours, s), num_dof,&
             ilocal_recon3(i)%matrix_1(1:num_neighbours, 1:num_vars, s), num_neighbours, &
             beta, ilocal_recon3(i)%sol(1:num_dof, 1:num_vars, s), num_dof)
     end do
  end do

  time_end2 = benchtime() - time_start2

  !$omp end target data

  time_end1 = benchtime() - time_start1

  check = 0.0
  ! verification
  do i=1,num_elems
     check = check +  sum(ilocal_recon3(i)%sol(:,:,:)**2)
  end do
  ! TIMEIT
  write (*, '(a35, es15.7)') "TIME: tot (s)  :", time_end1
  write (*, '(a35, es15.7)') "TIME: GPU (s)  :", time_end2
  write (*, '(a35, es15.7)') "GFLOP:GMYDGM(s):", 1.0e-9*(nflop/time_end2)
  write (*, '(a35, es15.7)') "check val    : ", check

end program prototype
