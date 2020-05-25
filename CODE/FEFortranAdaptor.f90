module tcp
  use iso_c_binding
  use DECLARATION
  USE MPIINFO
  implicit none
  public :: testcoprocessor, coprocess_grid
    contains

subroutine testcoprocessor(nxstart,nxend,n_x,n_y,n_z,step,time)
  use iso_c_binding
  implicit none

  integer :: numtasks, rank, ierr
  integer, intent(in) :: n_x,n_y,n_z,step
  integer, intent(in) :: nxstart,nxend
  real*8, intent(in) :: time
  integer :: flag
  call requestdatadescription(step,time,flag)

  if (flag .ne. 0) then
    call needtocreategrid(flag)
    if (flag .ne. 0) then
      call createcpimagedata(nxstart,nxend,n_x,n_y,n_z)
    end if
    ! adding //char(0) appends the C++ terminating character
    ! to the Fortran array
    ! call addfield(psi01,"psi01"//char(0))
    call coprocess()
  end if

  return

end subroutine

subroutine coprocess_grid(step,time)
  use iso_c_binding
  implicit none

  integer, intent(in) :: step
  real*8, intent(in) :: time
  integer :: flag
  !
  call requestdatadescription(step,time,flag)
  !
  if (flag .ne. 0) then
    call needtocreategrid(flag)
    if (flag .ne. 0) then
      call testfunction(pointSet, kmaxn, el_connect, XMPIELRANK(N), N, ISIZE)
    end if

    call addfield(scalarRU, "RU"//char(0))
    call addfield(scalarRV, "RV"//char(0))
    call addfield(scalarE, "E"//char(0))
    call coprocess()
  end if

  return
end subroutine

end module tcp
