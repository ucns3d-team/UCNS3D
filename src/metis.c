module Calling_METIS

  !use constants,  only : p2 !this is for double precision
  use iso_c_bind            !inserted later

  implicit none

  !integer                                    :: ne, nn              !modified
  integer(c_int)                              :: ne, nn 
  !integer,  dimension(:), allocatable        :: eptr, eind          !modified
  integer(c_int),  dimension(:), allocatable  :: eptr, eind
  !integer,  dimension(:), allocatable        :: vwgt, vsize         !modified
  type(c_ptr)                                 :: vwgt, vsize         
  !integer                                    :: nparts              !modified
  integer(c_int)                              :: nparts
  !real(p2), dimension(:), allocatable        :: tpwgts              !modified 
  type(c_ptr)                                 :: tpwgts      
  !integer,  dimension(0:39)                  :: opts                !modified
  integer(c_int),  dimension(0:39)            :: opts        
  !integer                                    :: objval              !modified
  integer(c_int)                              :: objval
  !integer,  dimension(:), allocatable        :: epart, npart        !modified 
  integer(c_int),  dimension(:), allocatable  :: epart, npart 

  interface
    subroutine METIS_PartMeshNodal( ne, nn, eptr, eind, vwgt, vsize, nparts, tpwgt, &
                                    opts, objval, epart, npart) bind(c)
      use intrinsic        :: iso_c_binding
      !use constants,  only  : p2

      implicit none

      integer (c_int),                  intent(in)  :: ne, nn
      integer (c_int), dimension(*),    intent(in)  :: eptr, eind
      !integer (c_int), dimension(*),    intent(in) :: vwgt, vsize  !modified
      type(c_ptr),                          value   :: vwgt, vsize   
      integer (c_int),                  intent(in)  :: nparts
      !real(c_double),  dimension(*),    intent(in) :: tpwgt        !modified
      type(c_ptr),                          value   :: tpwgt
      integer (c_int), dimension(0:39), intent(in)  :: opts
      integer (c_int),                  intent(out) :: objval
      integer (c_int), dimension(*),    intent(out) :: epart
      integer (c_int), dimension(*),    intent(out) :: npart

    end subroutine METIS_PartMeshNodal  
  end interface
end module