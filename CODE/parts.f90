module partition
USE MPIINFO
use DECLARATION


implicit none

 contains



subroutine Partitioner5(n,imaxe,imaxn,xmpie,ieshape)
   !> @brief
!> This subroutine partitions the mesh using Metis

	integer,intent(in)::n,imaxe,imaxn
	integer,ALLOCATABLE,DIMENSION(:),intent(in)::IESHAPE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE
	INTEGER::I,J,K,IOS,IOX,IOZ,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,IOY,IHE,ITRI,idu
	CHARACTER(LEN=12)::METFILE,CELFILE
      integer::posa
      real::xc,yc
character (len=10) :: t,f,ss
		CELFILE='GRID.cel'
		METFILE='GRID'
		
		OPEN(9,FILE=METFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	
	OPEN(8,FILE=CELFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	WRITE(9,*)IMAXE

	IHE=0
	ITRI=0


	if (dimensiona.eq.3)then
	DO J=1,IMAXE
		READ(8,*)k,I1,I2,I3,I4,i5,i6,i7,i8
 	      
		
		
		WRITE(9,"(8I10)")I1,I2,I3,I4,i5,i6,i7,i8
		
		
		!WRITE(9,"(8I10)")I1,I2,I3,I4
	END DO
	else
	DO J=1,IMAXE
		READ(8,*)k,I1,I2,I3,I4
 	      
		
		IF ((I3.NE.I4))THEN
		
		WRITE(9,"(4I10)")I1,I2,I3,I4
		Else
! 
		WRITE(9,"(3I10)")I1,I2,I3

		END IF
		
		
		!WRITE(9,"(8I10)")I1,I2,I3,I4
	END DO


	end if
	CLOSE(8)
	CLOSE(9)

	 
    if (k.ge.0)then
      do idu=1,itold
	yc=xc+yc**2
	xc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	xc=yc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	yc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if


 posa=isize
 WRITE(t,FMT='(I10)') posa
 f=TRIM(ADJUSTL(t))
 ss=TRIM(ADJUSTR(f))
 call system ('pwd')

 
   if (k.ge.0)then
      do idu=1,itold
	yc=xc+yc**2
	xc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	xc=yc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	yc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if
 
 
 call system ('./mpmetis STAR'//ss) 




 
 if (k.ge.0)then
      do idu=1,itold
	yc=xc+yc**2
	xc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	xc=yc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	yc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

  CALL system ('mv GRID.epart.* GRID.epart')

 if (k.ge.0)then
      do idu=1,itold
	yc=xc+yc**2
	xc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	xc=yc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	yc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if
  
  
 OPEN(10,FILE='GRID.epart',FORM='FORMATTED',STATUS='OLD',ACTION='READ')

do i=1,imaxe
  read(10,*)K
  XMPIE(I)=K
end do
 close(10)


if (k.ge.0)then
      do idu=1,itold
	yc=xc+yc**2
	xc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	xc=yc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
	yc=xc*idu**2+sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

CALL system ('rm -rf GRID.*part* GRID')


end subroutine Partitioner5





Subroutine Partitioner1(n,IMAXE,imaxn,xmpie,ieshape)
   !> @brief
!> This subroutine partitions the mesh using Metis
use ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
EXTERNAL METIS_PartMeshDual
EXTERNAL METIS_PartMeshNodal
integer,ALLOCATABLE,DIMENSION(:),intent(in)::IESHAPE
!!!!!!!!!!!!!!!!
TYPE::elementglobal
  INTEGER::ELEMENTGLID,NodeID1,NodeID2,NodeID3,NodeID4,NodeID5,NodeID6,NodeID7,NodeID8!,vweight!,vsize
END TYPE elementglobal
!!!!!!!!!!!!!!!
real::average
integer::maxi

Integer ::i,j,elementid,node1,node2,node3,node4,node5,node6,node7,node8,counternodes,k
Integer ::vweight
INTEGER,INTENT(IN)::N,IMAXE,imaxn
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE

! Real,allocatable,dimension(:) :: tpwgts
! type(c_ptr):: vwgt, vsize 
! type(c_ptr):: tpwgts  
integer,allocatable,dimension(:) :: testar,pweight
integer(c_int)::imaxee,imaxnn,ncommon1,isizee,objval
integer(c_int),allocatable,dimension(:)::AllnodesPTR,ALLNODES
integer(c_int),allocatable,dimension(:)::xmpiee,XMPIDUMB,vwgt, vsize
real(c_float),allocatable,dimension(:) :: tpwgts
integer(c_int),dimension(0:39)::options
TYPE(elementglobal),ALLOCATABLE,DIMENSION(:)::Elements
! use intrinsic        :: iso_c_binding
!  vwgt   = c_null_ptr    !added
!  vsize  = c_null_ptr    !added
!  tpwgts  = c_null_ptr    !added  
! allocate(tpwgts(isize))

! use intrinsic        :: iso_c_binding



! print*,'I am here in the start of the partionioner'
imaxee=imaxe
imaxnn=imaxn
isizee=isize
! allocate(options(0:39))
allocate(testar(imaxe))


! print*,'I am here in the start of the OPTIONS1'
call METIS_SetDefaultOptions(options)
! print*,'I am here in the start of the OPTIONS2'

 !options(7)=1
 counternodes=0
!  print*,imaxee
! print*,options


Allocate(Elements(imaxee))
! print*,ALLNODESGLOBALL
Allocate(ALLNODES(ALLNODESGLOBALL))
Allocate(AllnodesPTR(imaxe+1))
      OPEN(1112,FILE='GRID.cel',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
 
AllnodesPTR(1)=0



Do i=1,imaxe


	Read(1112,*)elementid,node1,node2,node3,node4,node5,node6,node7,node8
    Elements(i)%ELEMENTGLID=elementid;Elements(I)%NodeID1=node1;Elements(i)%NodeID2=node2
    Elements(I)%NodeID3=node3;Elements(I)%NodeID4=node4;Elements(I)%NodeID5=node5
    Elements(I)%NodeID6=node6;Elements(I)%NodeID7=node7;Elements(I)%NodeID8=node8
	    
	    IF ((node3.EQ.node4).AND.(node5.NE.node6).AND.(node7.EQ.node8))THEN ! prism
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     
	
		     
		     
	    Else IF ((node5.NE.node6).AND.(node6.NE.node7).AND.(node7.NE.node8).AND.(node3.NE.node4))THEN ! HEXA
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node8
		

	    Else IF ((node3.EQ.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Tetra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		

	    Else 

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		

	       
	    End if

    AllnodesPTR(i+1)= counternodes


End do
   
close(1112)

 ncommon1=1

allocate(vwgt(imaxee))
allocate(vsize(imaxee))
allocate(tpwgts(isizee))

allocate(XMPIDUMB(imaxnn))
allocate(xmpiee(imaxee))


tpwgts(:)=1.0/(1.0*isizee)
vwgt(:)=1
vsize(:)=1


! 
! 

! ! print*,vwgt
! 
! print*,'I am here in the start of the partionioner 1'

!   print*,imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB
Call METIS_PartMeshDual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB)




! print*,'I am here in the start of the partionioner2'
! read(*,*)
!  METIS PartMeshDual(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *vwgt, idx t *vsize,
!  idx t *ncommon, idx t *nparts, real t *tpwgts, idx t *options, idx t *objval,
! idx t *epart, idx t *npart)



!  METIS PartMeshDual(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *vwgt, idx t *vsize,
! idx t *ncommon, idx t *nparts, real t *tpwgts, idx t *options, idx t *objval,
! idx t *epart, idx t *npart)
! open(8841,File='GRID.epart',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
! write(8841,*)objval
! XMPIE(I)
do i=1,imaxe
XMPIE(i)=xmpiee(i)
! write(8841,*)XMPIE(i)
end do
! close(8841)




! write(8841,*)XMPIDUMB

allocate(pweight(0:isize-1))
pweight(:)=0
do i=1,imaxe
  pweight(xmpiee(i))=pweight(xmpiee(i))+1
end do
! read(*,*)
  average=pweight(0)
do i=1,isize-1
  average=pweight(i)+average
  
end do 
average=average/isize
maxi=maxval(pweight)


OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',position='append')
write(63,*)'load imbalance of partioner',maxi/average
CLOSE(63)


 
deallocate(testar)
deAllocate(Elements)
deAllocate(ALLNODES)
deAllocate(AllnodesPTR)
deallocate(vwgt)
deallocate(vsize)
deallocate(tpwgts)
deallocate(pweight)
deallocate(XMPIDUMB)
deallocate(xmpiee)



end subroutine Partitioner1

!---------------------------------------------------------------------------------------------!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Partioner with variable weights
Subroutine Partitioner2(n,IMAXE,imaxn,xmpie,ieshape)
   !> @brief
!> This subroutine partitions the mesh using Metis
use ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
EXTERNAL METIS_PartMeshDual
EXTERNAL METIS_PartMeshNodal
integer,ALLOCATABLE,DIMENSION(:),intent(in)::IESHAPE
!!!!!!!!!!!!!!!!
TYPE::elementglobal
  INTEGER::ELEMENTGLID,NodeID1,NodeID2,NodeID3,NodeID4,NodeID5,NodeID6,NodeID7,NodeID8!,vweight!,vsize
END TYPE elementglobal
!!!!!!!!!!!!!!!
real::average,average2
integer::maxi,maxi2

Integer ::i,j,elementid,node1,node2,node3,node4,node5,node6,node7,node8,counternodes,k
Integer ::vweight
Integer,allocatable,dimension(:,:,:) :: Compweight,Commweight
INTEGER,INTENT(IN)::N,IMAXE,imaxn
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE
! Real,allocatable,dimension(:) :: tpwgts
! type(c_ptr):: vwgt, vsize 
! type(c_ptr):: tpwgts  
integer,allocatable,dimension(:) :: testar,pweight,pweight2
integer(c_int)::imaxee,imaxnn,ncommon1,isizee,objval
integer(c_int),allocatable,dimension(:)::AllnodesPTR,ALLNODES
integer(c_int),allocatable,dimension(:)::xmpiee,XMPIDUMB,vwgt, vsize
real(c_float),allocatable,dimension(:) :: tpwgts
integer(c_int),dimension(0:39)::options

TYPE(elementglobal),ALLOCATABLE,DIMENSION(:)::Elements
! use intrinsic        :: iso_c_binding
!  vwgt   = c_null_ptr    !added
!  vsize  = c_null_ptr    !added
!  tpwgts  = c_null_ptr    !added  
! allocate(tpwgts(isize))

! use intrinsic        :: iso_c_binding



! print*,'I am here in the start of the partionioner'
imaxee=imaxe
imaxnn=imaxn
isizee=isize
! allocate(options(0:39))
allocate(testar(imaxe))
call METIS_SetDefaultOptions(options)

 !options(7)=1 !C or FORTRAN numbering
! options(1) = 1 !cut=0 or volume=1 objective
 counternodes=0

Allocate(Elements(imaxee))
Allocate(ALLNODES(ALLNODESGLOBALL))
Allocate(AllnodesPTR(imaxe+1))
      OPEN(1112,FILE='GRID.cel',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
 
AllnodesPTR(1)=0



Do i=1,imaxe


	Read(1112,*)elementid,node1,node2,node3,node4,node5,node6,node7,node8
    Elements(i)%ELEMENTGLID=elementid;Elements(I)%NodeID1=node1;Elements(i)%NodeID2=node2
    Elements(I)%NodeID3=node3;Elements(I)%NodeID4=node4;Elements(I)%NodeID5=node5
    Elements(I)%NodeID6=node6;Elements(I)%NodeID7=node7;Elements(I)%NodeID8=node8
	    
	    IF ((node3.EQ.node4).AND.(node5.NE.node6).AND.(node7.EQ.node8))THEN ! prism
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     
	
		     
		     
	    Else IF ((node5.NE.node6).AND.(node6.NE.node7).AND.(node7.NE.node8).AND.(node3.NE.node4))THEN ! HEXA
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node8
		

	    Else IF ((node3.EQ.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Tetra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		

	    Else 

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		

	       
	    End if

    AllnodesPTR(i+1)= counternodes


End do
close(1112)
ncommon1=1

allocate(vwgt(imaxee))
allocate(vsize(imaxee))
allocate(tpwgts(isizee))
allocate(XMPIDUMB(imaxnn))
allocate(xmpiee(imaxee))


!! 
!!! Compweight(type of scheme : order of scheme : element shape)

allocate(Compweight(1:3,1:10,1:4))
allocate(Commweight(1:3,1:10,1:4))
tpwgts(:)=1.0/(1.0*isizee)

 Compweight(:,:,:)=1 ;  Commweight(:,:,:)=1
!shapeweight(1)=12 !Hexa shapeweight(2)=4 !Tetra shapeweight(3)=6 !Pyra shapeweight(4)=8 !Prism
 Compweight(:,:,1) = 6
 Compweight(:,:,2) = 4
 Compweight(:,:,3) = 5
 Compweight(:,:,4) = 5
 
!  Compweight(2,2,1) = 1
!  Compweight(2,2,2) = 1
!  Compweight(2,2,3) = 1
!  Compweight(2,2,4) = 1 
!  
!  Compweight(2,3,1) = 1
!  Compweight(2,3,2) = 1
!  Compweight(2,3,3) = 1
!  Compweight(2,3,4) = 1

  Compweight(3,:,1) = 7
 Compweight(3,:,2) = 5
 Compweight(3,:,3) = 6
 Compweight(3,:,4) = 6
  
!  Compweight(3,2,1) = 7
!  Compweight(3,2,2) = 5
!  Compweight(3,2,3) = 6
!  Compweight(3,2,4) = 6
! 
!  Compweight(3,6,1) = 7
!  Compweight(3,6,2) = 5
!  Compweight(3,6,3) = 6
!  Compweight(3,6,4) = 6
!  
!  Compweight(3,7,1) = 7
!  Compweight(3,7,2) = 5
!  Compweight(3,7,3) = 6
!  Compweight(3,7,4) = 6
! 
!  Compweight(3,3,1) = 7
!  Compweight(3,3,2) = 5
!  Compweight(3,3,3) = 6
!  Compweight(3,3,4) = 6
! 
!  Compweight(3,4,1) = 7
!  Compweight(3,4,2) = 5
!  Compweight(3,4,3) = 6
!  Compweight(3,4,4) = 6
!  
!  Compweight(3,5,1) = 7
!  Compweight(3,5,2) = 5
!  Compweight(3,5,3) = 6
!  Compweight(3,5,4) = 6
!!!!!!!!!!!!!!!!!!!!!!!
 Commweight(:,:,1) = 6
 Commweight(:,:,2) = 4
 Commweight(:,:,3) = 5
 Commweight(:,:,4) = 5

!  Commweight(:,4,1) = 24
!  Commweight(:,4,2) = 12
!  Commweight(:,4,3) = 16
!  Commweight(:,4,4) = 18
 
!  Commweight(:,4,1) = 24
!  Commweight(:,4,2) = 12
!  Commweight(:,4,3) = 16
!  Commweight(:,4,4) = 18

!  Commweight(:,5,1) = 54
!  Commweight(:,5,2) = 24
!  Commweight(:,5,3) = 33
!  Commweight(:,5,4) = 39
! 
!  Commweight(:,6,1) = 96
!  Commweight(:,6,2) = 40
!  Commweight(:,6,3) = 56
!  Commweight(:,6,4) = 68
! 
!  Commweight(:,7,1) = 150
!  Commweight(:,7,2) = 96
!  Commweight(:,7,3) = 120
!  Commweight(:,7,4) = 150
 
Do i=1,imaxee
    vwgt(i)=Compweight(spatiladiscret,spatialorder,IESHAPE(i))
end do

Do i=1,imaxee
    vsize(i)=Commweight(1,1,IESHAPE(i))
end do
  

deallocate(Commweight)
deallocate(Compweight)

 Call METIS_PartMeshDual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB)
! Call METIS_PartMeshNodal(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB)
do i=1,imaxe
  XMPIE(i)=xmpiee(i)
end do

allocate(pweight(0:isize-1))
allocate(pweight2(0:isize-1))
pweight(:)=0;pweight2(:)=0
do i=1,imaxe
  pweight(xmpiee(i))=pweight(xmpiee(i))+1
  pweight2(xmpiee(i))=pweight2(xmpiee(i))+vwgt(i)
end do
  average=pweight(0)
  average2=pweight2(0)
do i=1,isize-1
  average=pweight(i)+average
  average2=pweight2(i)+average2
end do 
average=average/isize
average2=average2/isize
maxi=maxval(pweight)
maxi2=maxval(pweight2)

OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',position='append')
write(63,*)'load imbalance of partioner',maxi/average
write(63,*)'load imbalance of partioner based on element weights',maxi2/average2
CLOSE(63)

 
deallocate(testar)
deAllocate(Elements)
deAllocate(ALLNODES)
deAllocate(AllnodesPTR)
deallocate(vwgt)
deallocate(vsize)
deallocate(tpwgts)
deallocate(pweight)
deallocate(pweight2)
deallocate(XMPIDUMB)
deallocate(xmpiee)



end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Partioner with variable weights
Subroutine Partitioner4(n,IMAXE,imaxn,xmpie,ieshape)
   !> @brief
!> This subroutine partitions the mesh using Metis
use ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
EXTERNAL METIS_PartMeshDual
EXTERNAL METIS_PartMeshNodal
integer,ALLOCATABLE,DIMENSION(:),intent(in)::IESHAPE
!!!!!!!!!!!!!!!!
TYPE::elementglobal
  INTEGER::ELEMENTGLID,NodeID1,NodeID2,NodeID3,NodeID4,NodeID5,NodeID6,NodeID7,NodeID8!,vweight!,vsize
END TYPE elementglobal
!!!!!!!!!!!!!!!
real::average,average2
integer::maxi,maxi2

Integer ::i,j,elementid,node1,node2,node3,node4,node5,node6,node7,node8,counternodes,k
Integer ::vweight
Integer,allocatable,dimension(:,:,:) :: Compweight,Commweight
INTEGER,INTENT(IN)::N,IMAXE,imaxn
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE
! Real,allocatable,dimension(:) :: tpwgts
! type(c_ptr):: vwgt, vsize 
! type(c_ptr):: tpwgts  
integer,allocatable,dimension(:) :: testar,pweight,pweight2
integer(c_int)::imaxee,imaxnn,ncommon1,isizee,objval
integer(c_int),allocatable,dimension(:)::AllnodesPTR,ALLNODES
integer(c_int),allocatable,dimension(:)::xmpiee,XMPIDUMB,vwgt, vsize
real(c_float),allocatable,dimension(:) :: tpwgts
integer(c_int),dimension(0:39)::options

TYPE(elementglobal),ALLOCATABLE,DIMENSION(:)::Elements
! use intrinsic        :: iso_c_binding
!  vwgt   = c_null_ptr    !added
!  vsize  = c_null_ptr    !added
!  tpwgts  = c_null_ptr    !added  
! allocate(tpwgts(isize))

! use intrinsic        :: iso_c_binding



! print*,'I am here in the start of the partionioner'
imaxee=imaxe
imaxnn=imaxn
isizee=isize
! allocate(options(0:39))
allocate(testar(imaxe))
call METIS_SetDefaultOptions(options)

 !options(7)=1 !C or FORTRAN numbering
options(1) = 1 !cut=0 or volume=1 objective
 counternodes=0

Allocate(Elements(imaxee))
Allocate(ALLNODES(ALLNODESGLOBALL))
Allocate(AllnodesPTR(imaxe+1))
      OPEN(1112,FILE='GRID.cel',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
 
AllnodesPTR(1)=0



Do i=1,imaxe


	Read(1112,*)elementid,node1,node2,node3,node4,node5,node6,node7,node8
    Elements(i)%ELEMENTGLID=elementid;Elements(I)%NodeID1=node1;Elements(i)%NodeID2=node2
    Elements(I)%NodeID3=node3;Elements(I)%NodeID4=node4;Elements(I)%NodeID5=node5
    Elements(I)%NodeID6=node6;Elements(I)%NodeID7=node7;Elements(I)%NodeID8=node8
	    
	    IF ((node3.EQ.node4).AND.(node5.NE.node6).AND.(node7.EQ.node8))THEN ! prism
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     
	
		     
		     
	    Else IF ((node5.NE.node6).AND.(node6.NE.node7).AND.(node7.NE.node8).AND.(node3.NE.node4))THEN ! HEXA
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node8
		

	    Else IF ((node3.EQ.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Tetra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		

	    Else 

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		

	       
	    End if

    AllnodesPTR(i+1)= counternodes


End do
close(1112)
ncommon1=1

allocate(vwgt(imaxee))
allocate(vsize(imaxee))
allocate(tpwgts(isizee))
allocate(XMPIDUMB(imaxnn))
allocate(xmpiee(imaxee))


!! 
!!! Compweight(type of scheme : order of scheme : element shape)

allocate(Compweight(1:3,1:10,1:4))
allocate(Commweight(1:3,1:10,1:4))
tpwgts(:)=1.0/(1.0*isizee)

 Compweight(:,:,:)=1 ;  Commweight(:,:,:)=1
!shapeweight(1)=12 !Hexa shapeweight(2)=4 !Tetra shapeweight(3)=6 !Pyra shapeweight(4)=8 !Prism
!  Compweight(1,1,1) = 1
!  Compweight(1,1,2) = 1
!  Compweight(1,1,3) = 1
!  Compweight(1,1,4) = 1
! !  
!  Compweight(2,2,1) = 1
!  Compweight(2,2,2) = 1
!  Compweight(2,2,3) = 1
!  Compweight(2,2,4) = 1 
!  
!  Compweight(2,3,1) = 1
!  Compweight(2,3,2) = 1
!  Compweight(2,3,3) = 1
!  Compweight(2,3,4) = 1
! 
!  Compweight(3,3,1) = 7
!  Compweight(3,3,2) = 5
!  Compweight(3,3,3) = 6
!  Compweight(3,3,4) = 6
!  
!  Compweight(3,5,1) = 7
!  Compweight(3,5,2) = 5
!  Compweight(3,5,3) = 6
!  Compweight(3,5,4) = 6
! !!!!!!!!!!!!!!!!!!!!!!!
!  Commweight(1,1,1) = 12
!  Commweight(1,1,2) = 8
!  Commweight(1,1,3) = 10
!  Commweight(1,1,4) = 10
 
Do i=1,imaxee
    vwgt(i)=Compweight(spatiladiscret,spatialorder,IESHAPE(i))
end do

Do i=1,imaxee
    vsize(i)=Commweight(1,1,IESHAPE(i))
end do
  

deallocate(Commweight)
deallocate(Compweight)

! Call METIS_PartMeshDual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB)
 Call METIS_PartMeshNodal(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB)
do i=1,imaxe
  XMPIE(i)=xmpiee(i)
end do

allocate(pweight(0:isize-1))
allocate(pweight2(0:isize-1))
pweight(:)=0;pweight2(:)=0
do i=1,imaxe
  pweight(xmpiee(i))=pweight(xmpiee(i))+1
  pweight2(xmpiee(i))=pweight2(xmpiee(i))+vwgt(i)
end do
  average=pweight(0)
  average2=pweight2(0)
do i=1,isize-1
  average=pweight(i)+average
  average2=pweight2(i)+average2
end do 
average=average/isize
average2=average2/isize
maxi=maxval(pweight)
maxi2=maxval(pweight2)

OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',position='append')
write(63,*)'load imbalance of partioner',maxi/average
write(63,*)'load imbalance of partioner based on element weights',maxi2/average2
CLOSE(63)

 
deallocate(testar)
deAllocate(Elements)
deAllocate(ALLNODES)
deAllocate(AllnodesPTR)
deallocate(vwgt)
deallocate(vsize)
deallocate(tpwgts)
deallocate(pweight)
deallocate(pweight2)
deallocate(XMPIDUMB)
deallocate(xmpiee)



end subroutine











Subroutine Partitioner3(n,IMAXE,imaxn,xmpie,ieshape)
   !> @brief
!> This subroutine partitions the mesh using Metis
use ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
EXTERNAL METIS_PartMeshDual
EXTERNAL METIS_PartMeshNodal
integer,ALLOCATABLE,DIMENSION(:),intent(in)::IESHAPE
!!!!!!!!!!!!!!!!
TYPE::elementglobal
  INTEGER::ELEMENTGLID,NodeID1,NodeID2,NodeID3,NodeID4,NodeID5,NodeID6,NodeID7,NodeID8!,vweight!,vsize
END TYPE elementglobal
!!!!!!!!!!!!!!!
real::average
integer::maxi

Integer ::i,j,elementid,node1,node2,node3,node4,node5,node6,node7,node8,counternodes,k
Integer ::vweight
INTEGER,INTENT(IN)::N,IMAXE,imaxn
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE
! Real,allocatable,dimension(:) :: tpwgts
! type(c_ptr):: vwgt, vsize 
! type(c_ptr):: tpwgts  
integer,allocatable,dimension(:) :: testar,pweight
integer(c_int)::imaxee,imaxnn,ncommon1,isizee,objval
integer(c_int),allocatable,dimension(:)::AllnodesPTR,ALLNODES
integer(c_int),allocatable,dimension(:)::xmpiee,XMPIDUMB,vwgt, vsize
real(c_float),allocatable,dimension(:) :: tpwgts
integer(c_int),dimension(0:39)::options
TYPE(elementglobal),ALLOCATABLE,DIMENSION(:)::Elements
! use intrinsic        :: iso_c_binding
!  vwgt   = c_null_ptr    !added
!  vsize  = c_null_ptr    !added
!  tpwgts  = c_null_ptr    !added  
! allocate(tpwgts(isize))

! use intrinsic        :: iso_c_binding



! print*,'I am here in the start of the partionioner'
imaxee=imaxe
imaxnn=imaxn
isizee=isize
! allocate(options(0:39))
allocate(testar(imaxe))


! print*,'I am here in the start of the OPTIONS1'
call METIS_SetDefaultOptions(options)
! print*,'I am here in the start of the OPTIONS2'

 !options(7)=1
 counternodes=0
!  print*,imaxee
! print*,options


Allocate(Elements(imaxee))
! print*,ALLNODESGLOBALL

ALLNODESGLOBALL=8*imaxe

Allocate(ALLNODES(ALLNODESGLOBALL))
Allocate(AllnodesPTR(imaxe+1))
      OPEN(1112,FILE='GRID.cel',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
 
AllnodesPTR(1)=0



Do i=1,imaxe


	Read(1112,*)elementid,node1,node2,node3,node4,node5,node6,node7,node8
!    Elements(i)%ELEMENTGLID=elementid;Elements(I)%NodeID1=node1;Elements(i)%NodeID2=node2
!    Elements(I)%NodeID3=node3;Elements(I)%NodeID4=node4;Elements(I)%NodeID5=node5
!    Elements(I)%NodeID6=node6;Elements(I)%NodeID7=node7;Elements(I)%NodeID8=node8
	    
	!    IF ((node3.EQ.node4).AND.(node5.NE.node6).AND.(node7.EQ.node8))THEN ! prism!
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node1
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node2
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node3
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node5
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node6
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node7
		     
	
		     
		     
	 !   Else IF ((node5.NE.node6).AND.(node6.NE.node7).AND.(node7.NE.node8).AND.(node3.NE.node4))THEN ! HEXA
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node8
		

!	    Else IF ((node3.EQ.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Tetra
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node1
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node2
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node3
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node5
		

!	    Else 

!! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node1
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node2
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node3
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node4
!		     Counternodes=counternodes + 1
!		     ALLNODES(counternodes) =  node5
		

	       
!	    End if

    AllnodesPTR(i+1)= counternodes


End do
   
close(1112)

 ncommon1=1

allocate(vwgt(imaxee))
allocate(vsize(imaxee))
allocate(tpwgts(isizee))

allocate(XMPIDUMB(imaxnn))
allocate(xmpiee(imaxee))


tpwgts(:)=1.0/(1.0*isizee)
vwgt(:)=1
vsize(:)=1


! 
! 

! 
! print*,'I am here in the start of the partionioner 1'

  
Call METIS_PartMeshDual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, XMPIEe, XMPIDUMB)




! print*,'I am here in the start of the partionioner2'
! read(*,*)
!  METIS PartMeshDual(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *vwgt, idx t *vsize,
!  idx t *ncommon, idx t *nparts, real t *tpwgts, idx t *options, idx t *objval,
! idx t *epart, idx t *npart)



!  METIS PartMeshDual(idx t *ne, idx t *nn, idx t *eptr, idx t *eind, idx t *vwgt, idx t *vsize,
! idx t *ncommon, idx t *nparts, real t *tpwgts, idx t *options, idx t *objval,
! idx t *epart, idx t *npart)

! XMPIE(I)
do i=1,imaxe
XMPIE(i)=xmpiee(i)

end do







allocate(pweight(0:isize-1))
pweight(:)=0
do i=1,imaxe
  pweight(xmpiee(i))=pweight(xmpiee(i))+1
end do
! read(*,*)
  average=pweight(0)
do i=1,isize-1
  average=pweight(i)+average
  
end do 
average=average/isize
maxi=maxval(pweight)


OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',position='append')
write(63,*)'load imbalance of partioner',maxi/average
CLOSE(63)


 
deallocate(testar)
deAllocate(Elements)
deAllocate(ALLNODES)
deAllocate(AllnodesPTR)
deallocate(vwgt)
deallocate(vsize)
deallocate(tpwgts)
deallocate(pweight)
deallocate(XMPIDUMB)
deallocate(xmpiee)



end subroutine



Subroutine Partitioner6(n,IMAXE,imaxn,xmpie,ieshape)
   !> @brief
!> This subroutine partitions the mesh using ParMetis (preferred option large meshes can be partitioned)
use ISO_C_BINDING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Implicit None
EXTERNAL ParMETIS_V3_PartMeshKway
!!!!!!!!!!!!!!!
integer,ALLOCATABLE,DIMENSION(:),intent(in)::IESHAPE
real::average,average2,tsize
integer::maxi,maxi2
Integer ::i,j,elementid,node1,node2,node3,node4,node5,node6,node7,node8,counternodes,counternodes2,k
Integer ::vweight
Integer,allocatable,dimension(:,:,:) :: Compweight,Commweight
INTEGER,INTENT(IN)::N,IMAXE,imaxn
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE
integer,allocatable,dimension(:) :: testar,pweight,pweight2
integer(c_int)::imaxee,imaxnn,ncommon1,isizee,objval,NUMFLAG,ncon,NCOMMONNODES,edgecut,WGTFLAG
integer(c_int),allocatable,dimension(:)::AllnodesPTR,ALLNODES,parts_el
integer(c_int),allocatable,dimension(:)::xmpiee,XMPIDUMB,vwgt, vsize
integer(c_int), allocatable :: elmdist(:),eptr(:),eind(:),options(:)
real(c_float),allocatable,dimension(:) :: tpwgts
real(c_float),allocatable,dimension(:) :: UBVEC



WGTfLAG=0
imaxee=imaxe
imaxnn=imaxn
isizee=isize
tsize=isizee
allocate(xmpiee(1:isize))


do i=1,isize-1
  xmpiee(i)=(imaxee/isize)
!  
end do

  xmpiee(isize)=imaxe-(xmpiee(isize-1)*(isize-1))
!  



allocate(elmdist(isize+1))

    elmdist(1) = 1
    do i = 2,isize+1
      elmdist(i) = elmdist(i-1) + xmpiee(i-1)
!       
    end do
!    
  
if (dimensiona.eq.3)then
ALLNODESGLOBALL=8*xmpiee(n+1)
else
ALLNODESGLOBALL=4*xmpiee(n+1)
end if
  

Allocate(ALLNODES(1:ALLNODESGLOBALL))
Allocate(ePTR(1:xmpiee(n+1)+1))
    
 counternodes=0
 counternodes2=1
  ePTR(counternodes2)=1
  
  IF (BINIO.EQ.0)THEN
  
    if (dimensiona.eq.3)then
  
    
    OPEN(1112,FILE='GRID.cel',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    
Do i=1,imaxe
	if ((i.ge.elmdist(n+1)).and.(i.lt.elmdist(n+2)))then
	counternodes2=counternodes2+1
	Read(1112,*)elementid,node1,node2,node3,node4,node5,node6,node7,node8
	IF ((node3.EQ.node4).AND.(node5.NE.node6).AND.(node7.EQ.node8))THEN ! prism
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
   
	    Else IF ((node5.NE.node6).AND.(node6.NE.node7).AND.(node7.NE.node8).AND.(node3.NE.node4))THEN ! HEXA
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node8

	    Else IF ((node3.EQ.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Tetra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5

	    Else 

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
	    End if
	ePTR(counternodes2)= counternodes+1
	else
	Read(1112,*)
	end if
End do

close(1112)


    else
    OPEN(1112,FILE='GRID.cel',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
    
Do i=1,imaxe
	if ((i.ge.elmdist(n+1)).and.(i.lt.elmdist(n+2)))then
	counternodes2=counternodes2+1
	Read(1112,*)elementid,node1,node2,node3,node4
	IF ((node3.EQ.node4))THEN ! triangle
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		
   
	    

	    Else 
		    !qudrilateral

		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     
	    End if
	ePTR(counternodes2)= counternodes+1
	else
	Read(1112,*)
	end if
End do

close(1112)
    
    
    
    
    end if

    
    ELSE
    if (dimensiona.eq.3)then
  
    
    OPEN(1112,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
    
Do i=1,imaxe
	if ((i.ge.elmdist(n+1)).and.(i.lt.elmdist(n+2)))then
	counternodes2=counternodes2+1
	Read(1112)elementid,node1,node2,node3,node4,node5,node6,node7,node8
	IF ((node3.EQ.node4).AND.(node5.NE.node6).AND.(node7.EQ.node8))THEN ! prism
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
   
	    Else IF ((node5.NE.node6).AND.(node6.NE.node7).AND.(node7.NE.node8).AND.(node3.NE.node4))THEN ! HEXA
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node6
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node7
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node8

	    Else IF ((node3.EQ.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Tetra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5

	    Else 

! IF ((node3.NE.node4).AND.(node5.EQ.node6).AND.(node6.EQ.node7).AND.(node7.EQ.node8))THEN ! Pyra
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node5
	    End if
	ePTR(counternodes2)= counternodes+1
	else
	Read(1112)
	end if
End do

close(1112)


    else
    OPEN(1112,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
    
Do i=1,imaxe
	if ((i.ge.elmdist(n+1)).and.(i.lt.elmdist(n+2)))then
	counternodes2=counternodes2+1
	Read(1112)elementid,node1,node2,node3,node4
	IF ((node3.EQ.node4))THEN ! triangle
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		
   
	    

	    Else 
		    !qudrilateral

		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node1
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node2
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node3
		     Counternodes=counternodes + 1
		     ALLNODES(counternodes) =  node4
		     
	    End if
	ePTR(counternodes2)= counternodes+1
	else
	Read(1112)
	end if
End do

close(1112)
    
    
    
    
    end if
    
    
    
    END IF
  

 allocate(eind(1:counternodes))
 eind(1:counternodes)=Allnodes(1:counternodes)
 deAllocate(ALLNODES)

 
allocate(vwgt(xmpiee(n+1)))
vwgt=0;wgtflag=0;numflag=1;ncon=1;

 if (dimensiona.eq.3)then

NCOMMONNODES=3
else
NCOMMONNODES=2

end if


allocate(tpwgts(ncon*isize))
tpwgts = 1. / float(isize)
allocate(ubvec(ncon))
ubvec(1)=1.050000

allocate(options(3))
options=[0, 0, 0]

allocate(parts_el(xmpiee(n+1)))

   if (n.eq.0) then
      write(*,*) '------------------------------------------------------------------------'
      write(*,*) '                       ParMETIS Initiated                               '
      write(*,*) '------------------------------------------------------------------------'
    end if

    









CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)

 Call ParMETIS_V3_PartMeshKway(elmdist, eptr,eind,vwgt,WGTFLAG,numflag,NCON,NCOMMONNODES,isize,TPWGTS,UBVEC,OPTIONS,EDGECUT,PARTs_el,MPI_COMM_WORLD)
 
 
 

    if (n.eq.0) then
      write(*,*) '------------------------------------------------------------------------'
      write(*,*) '                       ParMETIS Operation Completed                     '
      write(*,*) '------------------------------------------------------------------------'
    end if







 




deAllocate(ePTR)
deallocate(eind)
deallocate(vwgt)
deallocate(tpwgts)

deallocate(UBVEC)
deallocate(options)





 CALL MPI_GATHER(PARTS_EL,XMPIEE(1),MPI_INTEGER,XMPIE,XMPIEE(1),MPI_INTEGER,ISIZE-1,MPI_COMM_WORLD,IERROR)
 
  IF (N.EQ.ISIZE-1)THEN
  XMPIE(elmdist(n+1):elmdist(n+2)-1)=PARTS_EL(1:XMPIEE(N))
  END IF
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 
 call MPI_BCAST(XMPIE,IMAXE,MPI_INTEGER,ISIZE-1,MPI_COMM_WORLD,IERROR) 

deAllocate(parts_el)
deallocate(xmpiee)
 deallocate(elmdist)
  XMPIE=XMPIE-1





end subroutine Partitioner6








end module partition

