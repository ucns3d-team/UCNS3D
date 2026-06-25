 module lapck
 use mpiinfo
 use declaration
  implicit none

  contains
 
  real function lxnorm(xqr,pdim)
  implicit none
   integer, intent(in) :: pdim
   real,allocatable,dimension(:),intent(in)::xqr
   integer i
 
    lxnorm = 0.0d0
    do i=1,pdim
     lxnorm = lxnorm + xqr(i)**2
    enddo
    lxnorm = sqrt(lxnorm)
  end function
  
  
  ! construct a householder vector v that 
  ! annihilates all but the first component of x
  subroutine house(xqr,vqr1,pdim)
	implicit none
    integer, intent(in) :: pdim
    real,allocatable,dimension(:),intent(inout)::xqr
    real,allocatable,dimension(:),intent(inout)::vqr1

    vqr1 = xqr
    vqr1(1) = xqr(1) + sign(1.0,xqr(1))*lxnorm(xqr,pdim)
    
  end subroutine


  ! construct a householder reflection matrix
  ! from a householder vector v 
  subroutine computehousematrix(pqr,vqr,ideg)
	implicit none
    integer, intent(in) :: ideg
    real,allocatable,dimension(:,:),intent(inout)::pqr
    real,allocatable,dimension(:),intent(inout)::vqr
    real ::vnorm
    integer:: i,j

    pqr = 0.0d0
    do i=1,ideg
      pqr(i,i) = 1.0d0
    enddo
   
    vnorm = lxnorm(vqr,ideg)
    vqr = vqr/vnorm
    
    do i=1,ideg
    do j=1,ideg
      pqr(i,j) = pqr(i,j) - 2.0d0*vqr(i)*vqr(j)
    enddo
    enddo
    
  end subroutine

 !%%%%%%%%%%%%%%%%%%
  subroutine transposematrix(qff,qtff,ideg)
	implicit none
    integer, intent(in) :: ideg
    real,allocatable,dimension(:,:),intent(in) :: qff
    real,allocatable,dimension(:,:),intent(inout) :: qtff
    integer i,j
    do i=1,ideg
    do j=1,ideg
      qtff (i,j) = qff(j,i)
    enddo
    enddo
  end subroutine
  
  
  
  subroutine transposematrix_dg(qff_dg,qtff_dg,ideg)
	implicit none
    integer, intent(in) :: ideg
    real,allocatable,dimension(:,:),intent(in) :: qff_dg
    real,allocatable,dimension(:,:),intent(inout) :: qtff_dg
    integer i,j
    do i=1,ideg
    do j=1,ideg
      qtff_dg (i,j) = qff_dg(j,i)
    enddo
    enddo
  end subroutine


  ! qr decomposition
  subroutine  qrdecomposition(lscqm,qff,rff,ideg)
	implicit none
integer :: mm_i
integer :: mm_j
integer :: mm_k
real :: mm_tmp(1:ideg,1:ideg)
    integer, intent(in) :: ideg
    real,allocatable,dimension(:,:),intent(in) :: lscqm
    real,allocatable,dimension(:,:),intent(inout) ::  qff,rff
    real::identity(ideg,ideg)
    real:: test(ideg)
    integer:: i,j,l,pdim
    real,allocatable,dimension(:)::xqr
    real,allocatable,dimension(:,:)::pqr
    real,allocatable,dimension(:)::vqr
    real,allocatable,dimension(:)::vqr1
    
    identity(1:ideg,1:ideg) = 0.0d0
    qff(1:ideg,1:ideg)=0.0d0
    rff(1:ideg,1:ideg)=0.0d0
    do i=1,ideg
      identity(i,i) = 1.0d0
    enddo

    
    qff(1:ideg,1:ideg) = identity(1:ideg,1:ideg)
    rff(1:ideg,1:ideg) = lscqm(1:ideg,1:ideg)
    
    allocate(pqr(1:ideg,1:ideg),vqr(1:ideg));pqr=0.0d0;vqr=0.0d0

    do l=1,ideg
     ! allocate vector and reflection matrix
     pdim = ideg-l+1
     allocate(vqr1(1:pdim),xqr(1:pdim)) 
     xqr = rff(l:ideg,l)
     ! compute the partial vector    
     call house(xqr,vqr1,pdim)
     vqr(1:ideg) = 0.0d0
     vqr(l:ideg) = vqr1
    
     
     
     ! compute the reflection matrix
     call computehousematrix(pqr,vqr,ideg)
     ! construct the q(l) matrix
     
       do mm_i=1,ideg
         do mm_j=1,ideg
           mm_tmp(mm_i,mm_j)=zero
           do mm_k=1,ideg
             mm_tmp(mm_i,mm_j)=mm_tmp(mm_i,mm_j)+pqr(mm_i,mm_k)*rff(mm_k,mm_j)
           end do
         end do
       end do
       rff(1:ideg,1:ideg)=mm_tmp(1:ideg,1:ideg)
       do mm_i=1,ideg
         do mm_j=1,ideg
           mm_tmp(mm_i,mm_j)=zero
           do mm_k=1,ideg
             mm_tmp(mm_i,mm_j)=mm_tmp(mm_i,mm_j)+qff(mm_i,mm_k)*pqr(mm_k,mm_j)
           end do
         end do
       end do
       qff(1:ideg,1:ideg)=mm_tmp(1:ideg,1:ideg)
    
     deallocate(vqr1,xqr)
    enddo

    deallocate(pqr,vqr)



 
  end subroutine
  
  
  
  subroutine  qrdecomposition_dg(lscqm_dg,qff_dg,rff_dg,ideg)
	implicit none
integer :: mm_i
integer :: mm_j
integer :: mm_k
real :: mm_tmp(1:ideg,1:ideg)
    integer, intent(in) :: ideg
    real,allocatable,dimension(:,:),intent(in) :: lscqm_dg
    real,allocatable,dimension(:,:),intent(inout) ::  qff_dg,rff_dg
    real::identity(ideg,ideg)
    real:: test(ideg)
    integer:: i,j,l,pdim
    real,allocatable,dimension(:)::xqr
    real,allocatable,dimension(:,:)::pqr
    real,allocatable,dimension(:)::vqr
    real,allocatable,dimension(:)::vqr1

    
    identity(1:ideg,1:ideg) = 0.0d0
    qff_dg(1:ideg,1:ideg)=0.0d0
    rff_dg(1:ideg,1:ideg)=0.0d0
    do i=1,ideg
      identity(i,i) = 1.0d0
    enddo

    
    qff_dg(1:ideg,1:ideg) = identity(1:ideg,1:ideg)
    rff_dg(1:ideg,1:ideg) = lscqm_dg(1:ideg,1:ideg)
    
    allocate(pqr(1:ideg,1:ideg),vqr(1:ideg));pqr=0.0d0;vqr=0.0d0

    do l=1,ideg
     ! allocate vector and reflection matrix
     pdim = ideg-l+1
     allocate(vqr1(1:pdim),xqr(1:pdim)) 
     xqr = rff_dg(l:ideg,l)
     ! compute the partial vector    
     call house(xqr,vqr1,pdim)
     vqr(1:ideg) = 0.0d0
     vqr(l:ideg) = vqr1
    
     
     
     ! compute the reflection matrix
     call computehousematrix(pqr,vqr,ideg)
     ! construct the q(l) matrix
     
       do mm_i=1,ideg
         do mm_j=1,ideg
           mm_tmp(mm_i,mm_j)=zero
           do mm_k=1,ideg
             mm_tmp(mm_i,mm_j)=mm_tmp(mm_i,mm_j)+pqr(mm_i,mm_k)*rff_dg(mm_k,mm_j)
           end do
         end do
       end do
       rff_dg(1:ideg,1:ideg)=mm_tmp(1:ideg,1:ideg)
       do mm_i=1,ideg
         do mm_j=1,ideg
           mm_tmp(mm_i,mm_j)=zero
           do mm_k=1,ideg
             mm_tmp(mm_i,mm_j)=mm_tmp(mm_i,mm_j)+qff_dg(mm_i,mm_k)*pqr(mm_k,mm_j)
           end do
         end do
       end do
       qff_dg(1:ideg,1:ideg)=mm_tmp(1:ideg,1:ideg)
    
     deallocate(vqr1,xqr)
    enddo

    deallocate(pqr,vqr)



 
  end subroutine

 end module
