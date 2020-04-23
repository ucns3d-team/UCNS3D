MODULE MEMORY
USE MPIINFO
USE DECLARATION


CONTAINS

SUBROUTINE ALLOCATE1(N)
!> @brief
!> This subroutine allocates memory 
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::ALLS




if (dimensiona.eq.3)then
ALLS=IGQRULES*IGQRULES*IGQRULES
ALLOCATE(VVA(3,3),VVA1(3,3),DETA(1))
allocate(sb(idegfree))
ALLOCATE(VVr1(alls),VVr2(alls),VVr3(alls),VVR4(ALLS),VVwg(alls),POX(1),POY(1),POZ(1),VVnpox(igqrules),VVnpoy(igqrules),VVnpoz(igqrules),VVwpox(igqrules),VVwpoy(igqrules),VVwpoz(igqrules))
ALLOCATE(VVnxi(8),VVneta(8),VVnzeta(8),VVxi(8),VVeta(8),VVzeta(8),VVnallx(8),VVnally(8),VVnallz(8),VVB(3),VVC(3),VVD(3),VVE(3),VVF(3),VVJACOBSURF(3),VVJACOBVOLUME(4))
else
ALLS=IGQRULES*IGQRULES
ALLOCATE(VVA(2,2),VVA1(2,2),DETA(1))
allocate(sb(idegfree))
ALLOCATE(VVr1(alls),VVr2(alls),VVr3(alls),VVR4(ALLS),VVwg(alls),VVnpox(igqrules),POX(1),POY(1),VVnpoy(igqrules),VVnpoz(igqrules),VVwpox(igqrules),VVwpoy(igqrules),VVwpoz(igqrules))
ALLOCATE(VVnxi(4),VVneta(4),VVnzeta(4),VVxi(4),VVeta(4),VVzeta(4),VVnallx(4),VVnally(4),VVnallz(4),VVB(3),VVC(3),VVD(3),VVE(3),VVF(3),VVJACOBSURF(3),VVJACOBVOLUME(4))
end if




END SUBROUTINE

subroutine allocate2
!> @brief
!> This subroutine allocates memory
implicit none
ALLOCATE(LIST(2000),INEB(6),IPERB(6),NODELIST(8))
END SUBROUTINE

subroutine allocate5
!> @brief
!> This subroutine allocates memory for the stencils
implicit none
ALLOCATE(ILOCALALLELG(N:N,xmpielrank(n),1,ISELEMT(N)))
ILOCALALLELG(:,:,:,:)=0
end subroutine


subroutine allocate6_1
!> @brief
!> This subroutine allocates memory for the stencils
implicit none
ALLOCATE (ILOCALALLELG3(1,1:ISELEM),ILOCALALLELGD(1,1:ISELEM))
ALLOCATE(STCON(N:N))
ALLOCATE(STCONC(N:N))
ALLOCATE(STCONS(N:N))
ALLOCATE(STCONG(N:N))
ALLOCATE(ISOSA(N:N))
ALLOCATE(IFSAT(N:N))
ALLOCATE(IISTART(N:N))
ALLOCATE(IX(N:N))
end subroutine

subroutine allocate6_2
!> @brief
!> This subroutine deallocates memory
implicit none
deALLOCATE(STCON)
deALLOCATE(STCONC)
deALLOCATE(STCONS)
deALLOCATE(STCONG)
deALLOCATE(ISOSA)
deALLOCATE(IFSAT)
deALLOCATE(IISTART)
deALLOCATE(IX)

deALLOCATE (ILOCALALLELG3,ILOCALALLELGD)
end subroutine

subroutine allocate7_1
!> @brief
!> This subroutine allocates memory
implicit none
ALLOCATE(BC(N:N,3))
ALLOCATE(VC(N:N,8,3))
ALLOCATE(VG(N:N,3))
ALLOCATE(IWHICHSTEN(N:N))
ALLOCATE(ISATISFIED(N:N))
ALLOCATE(ISHYAPE(N:N))
ALLOCATE(XCC(3),vgg(3))

end subroutine

subroutine allocate7_2
!> @brief
!> This subroutine deallocates memory
implicit none
deALLOCATE(BC)
deALLOCATE(VC)
deALLOCATE(VG)
deALLOCATE(IWHICHSTEN)
deALLOCATE(ISATISFIED)
deALLOCATE(ISHYAPE)
DEALLOCATE(XCC,vgg)
end subroutine


SUBROUTINE ALLOCATE3
!> @brief
!> This subroutine deallocates memory
DEALLOCATE(LIST,INEB,IPERB)
END SUBROUTINE

SUBROUTINE GLOBALDEA2(XMPIL,XMPIE)
!> @brief
!> This subroutine deallocates global lists
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIL,XMPIE
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  DEALLOCATE(XMPIL)
    DEALLOCATE(XMPIE)
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 END SUBROUTINE GLOBALDEA2

SUBROUTINE QUADALLOC(QPOINTS,QPOINTS2D,WEQUA2D,WEQUA3D,NUMBEROFPOINTS,NUMBEROFPOINTS2)
!> @brief
!> This subroutine allocates memory for the quadrature points
implicit none
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::WEQUA3D
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::WEQUA2D
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::QPOINTS2D
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::QPOINTS
INTEGER,INTENT(INout)::NUMBEROFPOINTS,NUMBEROFPOINTS2


if (dimensiona.eq.3)then
NUMBEROFPOINTS=MAX(QP_HEXA,QP_TETRA,QP_PYRA,QP_PRISM)
NUMBEROFPOINTS2=MAX(QP_QUAD,QP_TRIANGLE)
ALLOCATE(WEQUA2D(NUMBEROFPOINTS2))
ALLOCATE(QPOINTS2D(3,NUMBEROFPOINTS2))
ALLOCATE(WEQUA3D(NUMBEROFPOINTS))
ALLOCATE(QPOINTS(3,NUMBEROFPOINTS))
else
NUMBEROFPOINTS=MAX(QP_QUAD,QP_TRIANGLE)
NUMBEROFPOINTS2=QP_LINE
ALLOCATE(WEQUA2D(NUMBEROFPOINTS2))
ALLOCATE(QPOINTS2D(2,NUMBEROFPOINTS2))
ALLOCATE(WEQUA3D(NUMBEROFPOINTS))
ALLOCATE(QPOINTS(2,NUMBEROFPOINTS))

end if

WEQUA2D(:)=0.0D0
QPOINTS2D(:,:)=0.0D0
WEQUA3D(:)=0.0D0
QPOINTS(:,:)=0.0D0


ALLOCATE(weight_t2(NUMBEROFPOINTS2));weight_t2=zero

END SUBROUTINE QUADALLOC




SUBROUTINE DEQUADALLOC(QPOINTS,QPOINTS2D,WEQUA2D,WEQUA3D,NUMBEROFPOINTS,NUMBEROFPOINTS2)
!> @brief
!> This subroutine deallocates memory for the quadrature points
implicit none
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::WEQUA3D
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::WEQUA2D
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::QPOINTS2D
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::QPOINTS
INTEGER,INTENT(IN)::NUMBEROFPOINTS,NUMBEROFPOINTS2
DEALLOCATE(WEQUA2D)
DEALLOCATE(QPOINTS2D)
DEALLOCATE(WEQUA3D)
DEALLOCATE(QPOINTS)
END SUBROUTINE DEQUADALLOC

!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR FLUXES!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE SUMFLUX_ALLOCATION(N)
	!> @brief
!> This subroutine allocates memory for the fluxes
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::N
	INTEGER::I,KMAXE
	KMAXE=XMPIELRANK(N)
	ALLOCATE (RHS(KMAXE))
	
	if (dimensiona.eq.3)then
	IF (ITESTCASE.EQ.4)THEN

	  IF (TURBULENCE.EQ.1)THEN
	    ALLOCATE (RHST(KMAXE))
	  
	  END IF
	  IF ((TURBULENCE.EQ.0).AND.(PASSIVESCALAR.GT.0))THEN
	  ALLOCATE (RHST(KMAXE))
      	  END IF
	END IF
	
	DO I=1,KMAXE
		IF (ITESTCASE.LT.3)THEN
			ALLOCATE (RHS(I)%VAL(1))
		END IF
		IF (ITESTCASE.EQ.3)THEN
			ALLOCATE (RHS(I)%VAL(5))
		END IF
		IF (ITESTCASE.EQ.4)THEN
			ALLOCATE (RHS(I)%VAL(5))
			IF (TURBULENCE.EQ.1) THEN
			       ALLOCATE (RHST(I)%VAL(TURBULENCEEQUATIONS+PASSIVESCALAR))
			  

			END IF
			IF ((TURBULENCE.EQ.0).AND.(PASSIVESCALAR.GT.0)) THEN
			       ALLOCATE (RHST(I)%VAL(PASSIVESCALAR))
			END IF
			
		END IF
	END DO
	else
	IF (ITESTCASE.EQ.4)THEN

	  IF (TURBULENCE.EQ.1)THEN
	    ALLOCATE (RHST(KMAXE))
	  
	  END IF
	  IF ((TURBULENCE.EQ.0).AND.(PASSIVESCALAR.GT.0))THEN
	  ALLOCATE (RHST(KMAXE))
      	  END IF
	END IF
	
	DO I=1,KMAXE
		IF (ITESTCASE.LT.3)THEN
			ALLOCATE (RHS(I)%VAL(1))
		END IF
		IF (ITESTCASE.EQ.3)THEN
			ALLOCATE (RHS(I)%VAL(4))
		END IF
		IF (ITESTCASE.EQ.4)THEN
			ALLOCATE (RHS(I)%VAL(4))
			IF (TURBULENCE.EQ.1) THEN
			       ALLOCATE (RHST(I)%VAL(TURBULENCEEQUATIONS+PASSIVESCALAR))
			  

			END IF
			IF ((TURBULENCE.EQ.0).AND.(PASSIVESCALAR.GT.0)) THEN
			       ALLOCATE (RHST(I)%VAL(PASSIVESCALAR))
			END IF
			
		END IF
	END DO
	
	
	
	
	
	
	
	
	
	
	
	end if
	END SUBROUTINE SUMFLUX_ALLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE IMPALLOCATE(N)
!> @brief
!> This subroutine allocates memory for implicit time stepping
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE
KMAXE=XMPIELRANK(N)
if (dimensiona.eq.3)then
if (lowmemory.eq.0)then

ALLOCATE (IMPDIAG(KMAXE,1:5,1:5))
ALLOCATE (IMPOFF(KMAXE,6,1:5,1:5))
ALLOCATE (IMPdu(KMAXE,1:5+TURBULENCEEQUATIONS+PASSIVESCALAR))

IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(KMAXE,6,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF


IMPDIAG(:,:,:)=zero
IMPOFF(:,:,:,:)=zero
impdu(:,:)=zero

else

ALLOCATE (IMPDIAG(1,1:5,1:5))
ALLOCATE (IMPOFF(1,6,1:5,1:5))
ALLOCATE (IMPdu(KMAXE,1:5+TURBULENCEEQUATIONS+PASSIVESCALAR))
IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(1,6,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(1,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF
IMPDIAG(1,:,:)=zero
IMPOFF(1,:,:,:)=zero
impdu(:,:)=zero
end if

else

if (lowmemory.eq.0)then

ALLOCATE (IMPDIAG(KMAXE,1:4,1:4))
ALLOCATE (IMPOFF(KMAXE,4,1:4,1:4))
ALLOCATE (IMPdu(KMAXE,1:4+TURBULENCEEQUATIONS+PASSIVESCALAR))

IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(KMAXE,4,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF


IMPDIAG(:,:,:)=zero
IMPOFF(:,:,:,:)=zero
impdu(:,:)=zero

else

ALLOCATE (IMPDIAG(1,1:4,1:4))
ALLOCATE (IMPOFF(1,4,1:4,1:4))
ALLOCATE (IMPdu(KMAXE,1:4+TURBULENCEEQUATIONS+PASSIVESCALAR))
IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(1,4,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(1,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF
IMPDIAG(1,:,:)=zero
IMPOFF(1,:,:,:)=zero
impdu(:,:)=zero
end if

end if

END  SUBROUTINE IMPALLOCATE

SUBROUTINE VERTALLOCATION(N,vext,LEFTV,RIGHTV,VISCL,LAML)
!> @brief
!> This subroutine allocates memory for vertices
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::VeXt
	REAL,ALLOCATABLE,DIMENSION(:)::LEFTV,RIGHTV
	REAL,ALLOCATABLE,DIMENSION(:)::VISCL,LAML
	
	if (dimensiona.eq.3)then

	ALLOCATE (VEXT(8,3))
	ALLOCATE (DETERJACS(1))
	ALLOCATE (JACS(3,3))
	ALLOCATE (INVERSEJACS(3,3))
	ALLOCATE(LEFTV(NOF_VARIABLES))
	ALLOCATE(RIGHTV(NOF_VARIABLES))
	ALLOCATe(CORDS(3))
	Allocate(NODES_LIST(8,3))
	Allocate(ELEM_LISTD(6,4,3))

	else
	ALLOCATE (VEXT(4,2))
	ALLOCATE (DETERJACS(1))
	ALLOCATE (JACS(2,2))
	ALLOCATE (INVERSEJACS(2,2))
	ALLOCATE(LEFTV(NOF_VARIABLES))
	ALLOCATE(RIGHTV(NOF_VARIABLES))

	ALLOCATe(CORDS(2))
	Allocate(NODES_LIST(4,2))
	Allocate(ELEM_LISTD(2,3,2))



	end if
	
	if ( turbulence .eq. 1) Then
	    ALLOCATE(VISCL(1:4))
	    ALLOCATE(LAML(1:4))
	    ALLOCATE(ETVM(1))
	    !Modified on 19/6/2013
		ALLOCATE(TURBMV(2))	
	end if
	
	if ( turbulence .ne. 1) Then
	    ALLOCATE(VISCL(1:2))
	    ALLOCATE(LAML(1:2))
	    
	end if
	
  
   
	
	 VEXT=zero
	 
	END SUBROUTINE VERTALLOCATION


SUBROUTINE TIMING(N,CPUX1,CPUX2,CPUX3,CPUX4,CPUX5,CPUX6,TIMEX1,TIMEX2,TIMEX3,TIMEX4,TIMEX5,TIMEX6)
!> @brief
!> This subroutine allocates memory for the timers
IMPLICIT NONE
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CPUX1,CPUX2,CPUX3,CPUX4,CPUX5,CPUX6,TIMEX1,TIMEX2,TIMEX3,TIMEX4,TIMEX5,TIMEX6
INTEGER,INTENT(IN)::N
ALLOCATE (CPUX1(1))
ALLOCATE (CPUX2(1))
ALLOCATE (CPUX3(1))
ALLOCATE (CPUX4(1))
ALLOCATE (CPUX5(1))
ALLOCATE (CPUX6(1))
ALLOCATE (TIMEX1(1))
ALLOCATE (TIMEX2(1))
ALLOCATE (TIMEX3(1))
ALLOCATE (TIMEX4(1))
ALLOCATE (TIMEX5(1))
ALLOCATE (TIMEX6(1))
 CPUX1(1)=0.0; CPUX2(1)=0.0;  CPUX3(1)=0.0;  CPUX4(1)=0.0;  CPUX5(1)=0.0;  CPUX6(1)=0.0
  TIMEX1(1)=0.0; TIMEX2(1)=0.0; TIMEX3(1)=0.0;  TIMEX4(1)=0.0;  TIMEX5(1)=0.0;  TIMEX6(1)=0.0
  
END  SUBROUTINE TIMING






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR ELEMENTS!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SHALLOCATION(IESHAPE,IMAXE)
	!> @brief
!> This subroutine allocates memory for the shapes
	IMPLICIT NONE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IESHAPE
	INTEGER,INTENT(INOUT)::IMAXE
	ALLOCATE (IESHAPE(IMAXE))
	IESHAPE=0
	
	END SUBROUTINE SHALLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR ELEMENTS!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE SHDEALLOCATION(IESHAPE,IMAXE)
	!> @brief
!> This subroutine deallocates memory for the shapes
	IMPLICIT NONE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IESHAPE
	INTEGER,INTENT(INOUT)::IMAXE
	DEALLOCATE (IESHAPE)
	END SUBROUTINE SHDEALLOCATION





!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR ELEMENTS!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ELALLOCATION(N,XMPIE,XMPIELRANK,IELEM,IMAXE,IESHAPE,ITESTCASE,IMAXB,IBOUND,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)
	IMPLICIT NONE
	!> @brief
!> This subroutine allocates memory for the elements
	TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::IESHAPE
	INTEGER,INTENT(IN)::N
	INTEGER,INTENT(IN)::IMAXE,ITESTCASE,IMAXB
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:)::IBOUND
	INTEGER::I,J,K,LM,IEX,KMAXE,KK
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
	ALLOCATE(XMIN(N:N))
	ALLOCATE(YMIN(N:N))
	ALLOCATE(ZMIN(N:N))
	ALLOCATE(XMAX(N:N))
	ALLOCATE(YMAX(N:N))
	ALLOCATE(ZMAX(N:N))
	KK=0; I=0; J=0; LM=0 
	KMAXE=XMPIELRANK(N)
	ALLOCATE(IELEM(N:N,XMPIELRANK(N)))
! 	ALLOCATE(IELEM2(N:N,XMPIELRANK(N)))
	IF (ITESTCASE.LT.3)THEN
		IEX=1
	END IF
	IF (ITESTCASE.GE.3)THEN
	  IF (DIMENSIONA.EQ.3)THEN
		IEX=5
	  ELSE
	      IEX=4

	  END IF
	END IF


	END SUBROUTINE ELALLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR NODES!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE NODEALLOCATION(N,INODE,IMAXN,XMPIN,XMPINRANK,INODEN)
	!> @brief
!> This subroutine allocates memory for the nodes
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIN
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPINRANK
	INTEGER::I,J,K,LM,IEX,KMAXN,KK
	TYPE(NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
	TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INODEN
	INTEGER,INTENT(IN)::IMAXN
	
!  	ALLOCATE (INODEN(N:N,IMAXN))
!  		  inoden(:,:)%itor=0
	     
	
	END SUBROUTINE NODEALLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        SUBROUTINE NODEDEALLOCATION(N,INODE,IMAXN,XMPIN,XMPINRANK,INODEN)
        !> @brief
!> This subroutine deallocates memory from the nodes
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIN
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPINRANK
	INTEGER::I,J,K,LM,IEX,KMAXN,KK
	TYPE(NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
	TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INODEN
	INTEGER,INTENT(IN)::IMAXN
 	DEALLOCATE (INODE)
! 	DEALLOCATE (INODEN)
	END SUBROUTINE NODEDEALLOCATION
!---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






SUBROUTINE ALLOCATETURB(N,EDDYFL,EDDYFR)
!> @brief
!> This subroutine allocates memory for the turbulence equations
IMPLICIT NONE
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::EDDYFL,EDDYFR
INTEGER,INTENT(IN)::N

!What is the N used for????????
!When is this subroutine called???

IF (TURBULENCEMODEL.NE.2)THEN
ALLOCATE(EDDYFL(18))
ALLOCATE(EDDYFR(18))

!Added on 20/6/2013
!Allocate with 20 if...
ELSE
ALLOCATE(EDDYFL(20))
ALLOCATE(EDDYFR(20))
end if


END SUBROUTINE



SUBROUTINE XMPIALLOCATE(XMPIE,XMPIL,XMPIN,XMPIELRANK,XMPINRANK,IMAXE,IMAXN,NPROC)
!> @brief
!> This subroutine allocates memory for the global lists
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIELRANK
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPINRANK
INTEGER,INTENT(INOUT)::NPROC,IMAXE,IMAXN

ALLOCATE (XMPIE(IMAXE))
ALLOCATE (XMPIL(IMAXE))
ALLOCATE(XMPIELRANK(n:n))



XMPIE=0
XMPIL=0
XMPIELRANK=0

END  SUBROUTINE XMPIALLOCATE


SUBROUTINE DEALLOCATEMPI1(N)
!> @brief
!> This subroutine deallocates memory for boundary exchange
IMPLICIT NONE
INTEGER,INTENT(IN)::N

DEALLOCATE (IEXCHANGES1,IEXCHANGER1)





END SUBROUTINE DEALLOCATEMPI1

SUBROUTINE DEALLOCATEMPI2(N)
!> @brief
!> This subroutine deallocates memory for the stencil exchange
IMPLICIT NONE
INTEGER,INTENT(IN)::N

DEALLOCATE (IRECEXR1,IRECEXS1)





END SUBROUTINE DEALLOCATEMPI2


SUBROUTINE LOCAL_DELALLOCATION(ILOCAL_ELEM)
!> @brief
!> This subroutine deallocates memory for reconstruction prestoring
IMPLICIT NONE
TYPE(LOCAL_ELEMENT),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_ELEM
DEALLOCATE (ILOCAL_ELEM)
END SUBROUTINE LOCAL_DELALLOCATION
!----------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!SUBROUTINE USED FOR ALLOCATING LOCALISED!!!
!!!!!!!!!!!!!STENCIL ELEMENTS!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LOCAL_DNALLOCATION(ILOCAL_NODE)
!> @brief
!> This subroutine deallocates memory for reconstruction prestoring
IMPLICIT NONE
TYPE(LOCAL_NODE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_NODE
DEALLOCATE (ILOCAL_NODE)
END SUBROUTINE LOCAL_DNALLOCATION



SUBROUTINE LOCAL_ELALLOCATION(N,ILOCAL_ELEM)
!> @brief
!> This subroutine allocates memory for reconstruction prestoring
IMPLICIT NONE
TYPE(LOCAL_ELEMENT),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_ELEM
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I,M
KMAXE=XMPIELRANK(N)
M=TYPESTEN
! IF (TYPESTEN.EQ.13)THEN
! 	m=13
! END IF
! IF (TYPESTEN.EQ.5)THEN
! 	m=5
! END IF
! IF (TYPESTEN.EQ.1)THEN
! 	m=1
! END IF
!"I HAVE STARTED ALLOCATING ELEMENT MEMORY FOR TRANSFORMATION"
ALLOCATE (ILOCAL_ELEM(1))
I=1
	ILOCAL_ELEM(I)%IADMIS=M
	ALLOCATE (ILOCAL_ELEM(I)%IHEXG(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%IHEXL(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%IHEXB(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%IHEXN(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%ISHAPE(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%XXC(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%VOLUME(M,NUMNEIGHBOURS*iextend))
	ALLOCATE (ILOCAL_ELEM(I)%YYC(M,NUMNEIGHBOURS*iextend))
	IF (DIMENSIONA.EQ.3)THEN
	ALLOCATE (ILOCAL_ELEM(I)%ZZC(M,NUMNEIGHBOURS*iextend))	
	ILOCAL_ELEM(I)%ZZC=0.d0
	END IF
	ILOCAL_ELEM(I)%IHEXG=0
	ILOCAL_ELEM(I)%IHEXL=0
	ILOCAL_ELEM(I)%IHEXB=0
	ILOCAL_ELEM(I)%ISHAPE=0
	ILOCAL_ELEM(I)%IHEXN=0
	ILOCAL_ELEM(I)%VOLUME=0.d0
	ILOCAL_ELEM(I)%XXC=0.d0
	ILOCAL_ELEM(I)%YYC=0.d0
	

!"I HAVE FINISHED ALLOCATING ELEMENT MEMORY FOR TRANSFORMATION"

END SUBROUTINE LOCAL_ELALLOCATION

SUBROUTINE LOCAL_NALLOCATION(N,ILOCAL_NODE)
!> @brief
!> This subroutine allocates memory for reconstruction prestoring
IMPLICIT NONE
TYPE(LOCAL_NODE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_NODE
INTEGER,INTENT(IN)::N
INTEGER::I,M,KMAXE,K
KMAXE=XMPIELRANK(N)
! IF (TYPESTEN.EQ.7)THEN
! 	m=7
! END IF
! IF (TYPESTEN.EQ.5)THEN
! 	m=5
! END IF
! IF (TYPESTEN.EQ.1)THEN
! 	m=1
! END IF
M=TYPESTEN
IF (DIMENSIONA.EQ.3)THEN
  K=8
ELSE
  K=4
END IF

ALLOCATE (ILOCAL_NODE(1))
! ILOCAL_NODE(:)=0
I=1
	
		
	ILOCAL_NODE(I)%IADMIS=M
	ALLOCATE (ILOCAL_NODE(I)%NODCOUNT(M,NUMNEIGHBOURS*iextend,K))
	ALLOCATE (ILOCAL_NODE(I)%X(M,NUMNEIGHBOURS*iextend,K))
	ALLOCATE (ILOCAL_NODE(I)%Y(M,NUMNEIGHBOURS*iextend,K))
	IF (DIMENSIONA.EQ.3)THEN
	ALLOCATE (ILOCAL_NODE(I)%Z(M,NUMNEIGHBOURS*iextend,K))		
	ILOCAL_NODE(I)%Z=0.D0
	END IF
	ILOCAL_NODE(I)%NODCOUNT=0
	ILOCAL_NODE(I)%X=0.D0
	ILOCAL_NODE(I)%Y=0.D0
	


END SUBROUTINE LOCAL_NALLOCATION
! -------------------------------




SUBROUTINE LOCAL_RECONALLOCATION3(N,ILOCAL_RECON3)
!> @brief
!> This subroutine allocates memory for reconstruction prestoring
IMPLICIT NONE
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,M,IKG,ITRUE,kmaxe
INTEGER::inum_points,idum,ITARGET
REAL::PERC,PERDE,PERD,PER1,PER2,PER3,PER4,PER5,PER0,PEF0,PEF1,PEF2,PEF3,PEF4,PEF5,PERV,per01,pef01
KMAXE=XMPIELRANK(N)

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

ALLOCATE (ILOCAL_RECON3(KMAXE))

perde=0.0d0

if (fastest.ne.1)then
DO I=1,KMAXE	!for all elements
	SELECT CASE(IELEM(N,I)%ISHAPE)

	CASE(1,2,3,4)
	IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	M=IELEM(N,I)%ADMIS
	IF (EES.EQ.5)THEN
	IMAX2=NUMNEIGHBOURS2-1
	INUM2=NUMNEIGHBOURS2
	IDEG2=IDEGFREE2
	M2=IELEM(N,I)%ADMIS
	END IF
	if (fastest.ne.1)then
	    ALLOCATE (ILOCAL_RECON3(I)%INVCCJAC(3,3));ILOCAL_RECON3(I)%INVCCJAC(:,:)=0.0D0
		IDUM=0
		if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
		
		if (idum.eq.1)then
	   ALLOCATE (ILOCAL_RECON3(I)%VOLUME(1,INUM));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
	   else
	   ALLOCATE (ILOCAL_RECON3(I)%VOLUME(1,1));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
	   end if
		
		
!   	    ALLOCATE (ILOCAL_RECON3(I)%INVCTJAC(3,3));ILOCAL_RECON3(I)%INVCTJAC(:,:)=0.0D0
	    
	    ! IF (EES.EQ.5)THEN
	     ! ALLOCATE (ILOCAL_RECON3(I)%VOLUMEC(M2,INUM2));ILOCAL_RECON3(I)%VOLUMEC(:,:)=0.0D0
	    ! END IF
	    
	    ALLOCATE(ILOCAL_RECON3(I)%VEXT_REF(3));ILOCAL_RECON3(I)%VEXT_REF=0.0D0
	end if
	IF (FIRSTORDER.NE.1)THEN
	   IF (GREENGO.EQ.0)then
	  
	   IDUM=0;
                if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
                if (idum.eq.1)then
	      
				   ALLOCATE (ILOCAL_RECON3(I)%STENCILS(M,IMAX,IDEG));ILOCAL_RECON3(I)%STENCILS(:,:,:)=0.0d0
				   IF (EES.EQ.5)THEN
				   ALLOCATE (ILOCAL_RECON3(I)%STENCILSC(M2,IMAX2,IDEG2));ILOCAL_RECON3(I)%STENCILSC(:,:,:)=0.0d0
				   END IF
	   
				end if
	   end if
	   
	   
	   
! 	    ALLOCATE (ILOCAL_RECON3(I)%INVMAT(M,IDEG,IDEG));ILOCAL_RECON3(I)%INVMAT(:,:,:)=0.0D0
	    allocate (ILOCAL_RECON3(I)%invmat_stencilt(ideg,imax,M));ILOCAL_RECON3(I)%invmat_stencilt(:,:,:)=0.0d0
	    IF (EES.EQ.5)THEN
	    allocate (ILOCAL_RECON3(I)%invmat_stenciltC(ideg2,imax2,M2));ILOCAL_RECON3(I)%invmat_stenciltC(:,:,:)=0.0d0
	    END IF
	    
	END IF
	
	ALLOCATE (ILOCAL_RECON3(I)%IHEXG(M,INUM))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXL(M,INUM))
	
	IF (EES.EQ.5)THEN
	ALLOCATE (ILOCAL_RECON3(I)%IHEXGC(M2,INUM2))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXLC(M2,INUM2))
	
	
	END IF
	
	ITRUE=0
	DO J=1,TYPESTEN
		IKG=0
                    IF ((EES.NE.5).OR.(J.EQ.1))THEN
                            ITARGET=INUM
                    ELSE
                            ITARGET=INUM2
                    END IF
			DO K=1,ITARGET
				IF (ILOCALSTENCIL(N,I,J,K).GT.0)THEN
				IKG=IKG+1
				IF (XMPIE(ILOCALSTENCIL(N,I,J,K)).NE.N)THEN
				  ITRUE=1
				END IF
				END IF
			END DO

	END DO
	IF (ITRUE.EQ.0)THEN
	ILOCAL_RECON3(I)%LOCAL=1
	ELSE
	ILOCAL_RECON3(I)%LOCAL=0
	END IF
	

	IF (ILOCAL_RECON3(I)%LOCAL.EQ.0)THEN
	ALLOCATE (ILOCAL_RECON3(I)%IHEXB(M,INUM))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXN(M,INUM))
	IF (EES.EQ.5)THEN
	ALLOCATE (ILOCAL_RECON3(I)%IHEXbC(M2,INUM2))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXnC(M2,INUM2))
	
	
	END IF
	
	
	
	END IF

	if (iweno.eq.1)then
	ALLOCATE (ILOCAL_RECON3(I)%INDICATOR(IDEG,IDEG));ILOCAL_RECON3(I)%INDICATOR(:,:)=0.0D0
	IF (EES.EQ.5)THEN
	ALLOCATE (ILOCAL_RECON3(I)%INDICATORc(IDEG2,IDEG2));ILOCAL_RECON3(I)%INDICATORc(:,:)=0.0D0
	end if
	END IF
	


	CASE(5,6)
	IMAX=IELEM(N,I)%inumneighbours-1
	INUM=IELEM(N,I)%inumneighbours
	IDEG=IELEM(N,I)%iDEGFREE
	M=IELEM(N,I)%ADMIS
	
	IF (EES.EQ.5)THEN
	IMAX2=NUMNEIGHBOURS2-1
	INUM2=NUMNEIGHBOURS2
	IDEG2=IDEGFREE2
	M2=IELEM(N,I)%ADMIS
	END IF
	
	
	
	if (fastest.ne.1)then
	    ALLOCATE (ILOCAL_RECON3(I)%INVCCJAC(2,2));ILOCAL_RECON3(I)%INVCCJAC(:,:)=0.0D0
!   	    ALLOCATE (ILOCAL_RECON3(I)%INVCTJAC(2,2));ILOCAL_RECON3(I)%INVCTJAC(:,:)=0.0D0
	    !ALLOCATE (ILOCAL_RECON3(I)%VOLUME(M,INUM));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
	     !IF (EES.EQ.5)THEN
	     !ALLOCATE (ILOCAL_RECON3(I)%VOLUMEC(M2,INUM2));ILOCAL_RECON3(I)%VOLUMEC(:,:)=0.0D0
	    !END IF
	    ALLOCATE(ILOCAL_RECON3(I)%VEXT_REF(2));ILOCAL_RECON3(I)%VEXT_REF=0.0D0
	end if
	
	IDUM=0
		if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
		
		if (idum.eq.1)then
	   ALLOCATE (ILOCAL_RECON3(I)%VOLUME(1,INUM));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
	   else
	   ALLOCATE (ILOCAL_RECON3(I)%VOLUME(1,1));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
	   end if
	
	
	
	IF (FIRSTORDER.NE.1)THEN
	   IF (GREENGO.EQ.0)then
	  
	   IDUM=0;
                if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
                if (idum.eq.1)then
	      
	   ALLOCATE (ILOCAL_RECON3(I)%STENCILS(M,IMAX,IDEG));ILOCAL_RECON3(I)%STENCILS(:,:,:)=0.0d0
	   IF (EES.EQ.5)THEN
	   ALLOCATE (ILOCAL_RECON3(I)%STENCILSC(M2,IMAX2,IDEG2));ILOCAL_RECON3(I)%STENCILSC(:,:,:)=0.0d0
	   END IF
	   end if
	   end if
	   
	   
	   
	   
! 	    ALLOCATE (ILOCAL_RECON3(I)%INVMAT(M,IDEG,IDEG));ILOCAL_RECON3(I)%INVMAT(:,:,:)=0.0D0
	    allocate (ILOCAL_RECON3(I)%invmat_stencilt(ideg,imax,ielem(n,i)%admis));ILOCAL_RECON3(I)%invmat_stencilt(:,:,:)=0.0d0
	    IF (EES.EQ.5)THEN
	    allocate (ILOCAL_RECON3(I)%invmat_stenciltC(ideg2,imax2,M2));ILOCAL_RECON3(I)%invmat_stenciltC(:,:,:)=0.0d0
	    END IF
	END IF
	
	ALLOCATE (ILOCAL_RECON3(I)%IHEXG(M,INUM))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXL(M,INUM))
	if (initcond.eq.0)then
	allocate (ILOCAL_RECON3(I)%cond(7))
	end if
	IF (EES.EQ.5)THEN
	ALLOCATE (ILOCAL_RECON3(I)%IHEXGC(M2,INUM2))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXLC(M2,INUM2))
	
	
	END IF
	ITRUE=0
	DO J=1,TYPESTEN
		IKG=0
			 IF ((EES.NE.5).OR.(J.EQ.1))THEN
                            ITARGET=INUM
                    ELSE
                            ITARGET=INUM2
                    END IF
			DO K=1,ITARGET
				IF (ILOCALSTENCIL(N,I,J,K).GT.0)THEN
				IKG=IKG+1
				IF (XMPIE(ILOCALSTENCIL(N,I,J,K)).NE.N)THEN
				  ITRUE=1
				END IF
				END IF
			END DO

	END DO
	IF (ITRUE.EQ.0)THEN
	ILOCAL_RECON3(I)%LOCAL=1
	ELSE
	ILOCAL_RECON3(I)%LOCAL=0
	END IF
	

	IF (ILOCAL_RECON3(I)%LOCAL.EQ.0)THEN
	ALLOCATE (ILOCAL_RECON3(I)%IHEXB(M,INUM))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXN(M,INUM))
	IF (EES.EQ.5)THEN
	ALLOCATE (ILOCAL_RECON3(I)%IHEXbC(M2,INUM2))
	ALLOCATE (ILOCAL_RECON3(I)%IHEXnC(M2,INUM2))
	
	
	END IF
	END IF

	if (iweno.eq.1)then
	ALLOCATE (ILOCAL_RECON3(I)%INDICATOR(IDEG,IDEG));ILOCAL_RECON3(I)%INDICATOR(:,:)=0.0D0
        IF (EES.EQ.5)THEN
	ALLOCATE (ILOCAL_RECON3(I)%INDICATORc(IDEG2,IDEG2));ILOCAL_RECON3(I)%INDICATORc(:,:)=0.0D0
	end if
	END IF


      end select
      perc=ilocal_recon3(i)%local
      perd=kmaxe
      perde=perde+perc

    
	
END DO
END IF


! Q5=0 ;Q4=0 ;Q3=0 ;Q2=0 ;Q1=0; Q0=0; q01=0
! DO I=1,KMAXE
! 		PERDE=Ielem(n,I)%ADMIS
! 		PERDI=(IELEM(N,I)%IFCA+1)
! 		PERC=PERDE/PERDI
! 		IF (PERC.EQ.1.0) THEN
! 		Q5=Q5+1
! 		END IF 
! 		IF ((PERC.LT.1.0).AND.(PERC.GE.0.75)) THEN
! 		Q4=Q4+1
! 		END IF 
! 		IF ((PERC.LT.0.75).AND.(PERC.GE.0.5)) THEN
! 		Q3=Q3+1
! 		END IF 
! 		IF ((PERC.LT.0.5).AND.(PERC.GE.0.25)) THEN
! 		Q2=Q2+1
! 		END IF 
! 		IF ((PERC.LT.0.25).AND.(PERC.GT.0.0)) THEN
! 		Q1=Q1+1
! 		END IF 
! 		IF ((Ielem(n,I)%ADMIS).EQ.1) THEN
! 		Q0=Q0+1
! 		END IF
! 		IF ((Ielem(n,I)%ADMIS).ge.3) THEN
! 		Q01=Q01+1
! 		END IF
! END DO
! 	PER5=Q5	;PERV=KMAXE
! 	PER4=Q4	
! 	PER3=Q3	
! 	PER2=Q2	
! 	PER1=Q1	
! 	PER0=Q0
! 	PER01=Q01
! 	
! 	PEF5=(PER5/PERV)*100; PEF4=(PER4/PERV)*100; PEF3=(PER3/PERV)*100; PEF2=(PER2/PERV)*100
!      PEF1=(PER1/PERV)*100; PEF0=(PER0/PERV)*100; PEF01=(PER01/PERV)*100
	

	
	
	
	



END SUBROUTINE LOCAL_RECONALLOCATION3

SUBROUTINE ALLWEFF(WEFF,IDEGFREE)
!> @brief
!> This subroutine allocates memory for weno weights 
IMPLICIT NONE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::WEFF
INTEGER,INTENT(IN)::IDEGFREE
ALLOCATE (WEFF(1:IDEGFREE,1:IDEGFREE))
WEFF=zero

END SUBROUTINE ALLWEFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DEALCORDINATES1(N,IEXCORDR,IEXCORDS)
!> @brief
!> This subroutine deallocates memory for exchange of info between processes
IMPLICIT NONE
TYPE(EXCHANGE_CORD),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXCORDR
TYPE(EXCHANGE_CORD),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXCORDS
INTEGER,INTENT(IN)::N
DEALLOCATE(IEXCORDR)


END SUBROUTINE DEALCORDINATES1


SUBROUTINE DEALCORDINATES2
!> @brief
!> This subroutine deallocates memory for exchange of info between processes
 if (allocated(iexcords))DEALLOCATE(IEXCORDS)

END SUBROUTINE DEALCORDINATES2





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ALLOCATE_BASIS_FUNCTION(N,INTEG_BASIS,XMPIELRANK,IDEGFREE)
!> @brief
!> This subroutine allocates memory for basis function integrals
IMPLICIT NONE
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INTEG_BASIS
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER,INTENT(IN)::IDEGFREE,N
INTEGER::KMAXE,i
KMAXE=XMPIELRANK(N)
ALLOCATE(INTEG_BASIS(KMAXE))
do i=1,kmaxe
 allocate(INTEG_BASIS(i)%value(idegfree));INTEG_BASIS(i)%value(:)=zero
 if (ees.eq.5)then
 allocate(INTEG_BASIS(i)%valuec(idegfree2));INTEG_BASIS(i)%valuec(:)=zero
 end if
		  
end do

END SUBROUTINE ALLOCATE_BASIS_FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DEALLOCATE_BASIS_FUNCTION(N,INTEG_BASIS)
!> @brief
!> This subroutine deallocates memory for basis function integrals
IMPLICIT NONE
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INTEG_BASIS
INTEGER,INTENT(IN)::N
DEALLOCATE(INTEG_BASIS)
END SUBROUTINE DEALLOCATE_BASIS_FUNCTION



SUBROUTINE MEMORY1
!> @brief
!> This subroutine allocates memory for reconstruction matrices

ALLOCATE(MATRIX_A(1:IMAXDEGFREE,1:IDEGFREE),SOURCE_T(TURBULENCEEQUATIONS),LSQM(1:IMAXDEGFREE,1:IDEGFREE-1),VELLSQMAT(1:IDEGFREE-1,1:IDEGFREE-1),MATRIX_SOLUTION(1:IDEGFREE,1:1)&
,MATRIXFACE(1:1,1:IDEGFREE),JUSTCHECK(1:1,1:1),Q(1:IDEGFREE-1,1:IDEGFREE-1),R(1:IDEGFREE-1,1:IDEGFREE-1),QT(1:IDEGFREE-1,1:IDEGFREE-1)&
,INVR(1:IDEGFREE-1,1:IDEGFREE-1),LSCQM(1:IDEGFREE,1:IDEGFREE),QFF(1:IDEGFREE,1:IDEGFREE),RFF(1:IDEGFREE,1:IDEGFREE),QTFF(1:IDEGFREE,1:IDEGFREE)&
,INVRFF(1:IDEGFREE,1:IDEGFREE),AINVJT(1:dims,1:dims),MATRIX_B(1:IMAXDEGFREE),MATRIX_X(1:IDEGFREE),BASEFACEVAL(1:IDEGFREE)&
,BASEFACGVAL(1:IDEGFREE),MOMENT(1:IDEGFREE),WEIGHTINV(1:IMAXDEGFREE),VECTOR(1:IMAXDEGFREE),INTBS(1:IDEGFREE),PERMUTATION(1:IDEGFREE),&
PERMUTATIONG(1:IDEGFREE),XDER(1:IDEGFREE),YDER(1:IDEGFREE),ZDER(1:IDEGFREE),XXDER(1:IDEGFREE,NUMBEROFPOINTS2*6),YYDER(1:IDEGFREE,NUMBEROFPOINTS2*6),ZZDER(1:IDEGFREE,NUMBEROFPOINTS2*6))


END SUBROUTINE MEMORY1

SUBROUTINE MEMORY11
!> @brief
!> This subroutine deallocates memory used for reconstruction matrices

DEALLOCATE(MATRIX_A,LSQM,VELLSQMAT,MATRIX_SOLUTION,MATRIXFACE,JUSTCHECK,Q,R,QT&
,INVR,LSCQM,QFF,RFF,QTFF,INVRFF,MATRIX_B,MATRIX_X,BASEFACEVAL&
,BASEFACGVAL,MOMENT,WEIGHTINV,VECTOR,INTBS,PERMUTATION,&
PERMUTATIONG)



END SUBROUTINE MEMORY11

   SUBROUTINE LOCALSDEALLOCATION(N,XMPIELRANK,ILOCALSTENCIL,TYPESTEN,NUMNEIGHBOURS)
   !> @brief
!> This subroutine allocates memory for stencils
	IMPLICIT NONE
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),intent(inout)::ILOCALSTENCIL
	INTEGER,INTENT(IN)::NUMNEIGHBOURS
	INTEGER,INTENT(IN)::N
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,INTENT(IN)::TYPESTEN
	
	DEALLOCATE (ILOCALSTENCIL)
	END SUBROUTINE LOCALSDEALLOCATION
	
	
	SUBROUTINE share_ALLOCATION(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::kkd,kkd1,kkd2,i
	kkd=nof_variables+turbulenceequations+passivescalar
	kkd1=nof_variables
	kkd2=turbulenceequations+passivescalar
    !> @brief
    !> This subroutine allocates memory for solution vectors
	
	
	
	ALLOCATE (CLEFT(kkd),cleft_rot(kkd),cright(kkd),cright_rot(kkd),sl(1),sr(1),sm(1),HLLCFLUX(kkd),&
	RHLLCFLUX(kkd),ROTVL(kkd),ROTVR(kkd),SUBSON1(kkd1),SUBSON2(kkd1),SUBSON3(kkd1),CTURBL(kkd2),CTURBr(kkd2),&
	turbc1(kkd2),turbc2(kkd2),tempflux2(kkd1))
	ALLOCATE (TEMPFL(KKD),TEMPFR(KKD),FLSTAR(KKD),FRSTAR(KKD),ULSTAR(KKD),UrSTAR(KKD),TEMPUL(KKD),TEMPUR(KKD),FL(KKD),FR(KKD),RML(KKD2),RMR(KKD2))

if (dimensiona.eq.3)then
ALLOCATE(TRI(5,5),INVTRI(5,5),EIGVL(5,5),EIGVR(5,5),VECTCO(kkd),VECCOS(KKD),ROTVECT(kkd))
ALLOCATE(VEIGL(kkd1),VEIGR(kkd1),RVEIGL(kkd1),RVEIGR(kkd1))
else
ALLOCATE(EIGVL(4,4),EIGVR(4,4),VECTCO(kkd),ROTVECT(kkd),VECCOS(KKD))
ALLOCATE(VEIGL(kkd1),VEIGR(kkd1),RVEIGL(kkd1),RVEIGR(kkd1))
end if
      
	
	
	
	
	
	END SUBROUTINE share_ALLOCATION	
	
	SUBROUTINE U_C_ALLOCATION(N,XMPIELRANK,U_C,U_E,ITESTCASE,U_CT)
	   !> @brief
!> This subroutine allocates memory for solution vector
	IMPLICIT NONE
	TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C,U_CT	
	TYPE(U_EXACT),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_E
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,INTENT(IN)::ITESTCASE,N
	INTEGER::I,KMAXE
	KMAXE=XMPIELRANK(N)
	ALLOCATE (U_C(KMAXE))
	
	IF (ITESTCASE.LE.3)THEN
	ALLOCATE (U_E(KMAXE))
	end if
	if (( turbulence .eq. 1).OR.(PASSIVESCALAR.GT.0))THEN
	  Allocate(U_CT(kmaxe))
	END IF

	
	IF (DIMENSIONA.EQ.3)THEN
	DO I=1,KMAXE
		IF (ITESTCASE.Lt.3)THEN
			if (rungekutta.eq.3)then
			ALLOCATE (U_C(I)%VAL(3,1));U_C(I)%VAL=ZERO
			end if
			if (rungekutta.eq.4)then
			ALLOCATE (U_C(I)%VAL(6,1));U_C(I)%VAL=ZERO
			end if
			ALLOCATE (U_E(I)%VAL(1,1));U_E(I)%VAL=ZERO
		END IF
		
		
		
		
		IF (ITESTCASE.GE.3)THEN
                    IF((RUNGEKUTTA.EQ.10))THEN
                    ALLOCATE (U_C(I)%VAL(1,5));U_C(I)%VAL=ZERO
                    if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                        Allocate(U_CT(I)%VAL(1,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
                    Endif
                    END IF
                        IF (RUNGEKUTTA.EQ.8)THEN
                        ALLOCATE (U_C(I)%VAL(5,5));U_C(I)%VAL=ZERO
                        if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                            Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))  ;U_CT(I)%VAL=ZERO    
                        Endif
                        END IF
                        IF( (RUNGEKUTTA.EQ.11))THEN
                            IF (AVERAGING.EQ.1)THEN
                        ALLOCATE (U_C(I)%VAL(5,5));U_C(I)%VAL=ZERO
                        if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                            Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))    ;U_CT(I)%VAL=ZERO  
                        Endif
                            ELSE
                            ALLOCATE (U_C(I)%VAL(3,5));U_C(I)%VAL=ZERO
                            if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                                Allocate(U_CT(I)%VAL(3,turbulenceequations+PASSIVESCALAR))  ;U_CT(I)%VAL=ZERO    
                            Endif
                            end if
                        END IF
                        IF (RUNGEKUTTA.EQ.3)THEN
                            IF (AVERAGING.EQ.1)THEN
                            ALLOCATE (U_C(I)%VAL(5,5));U_C(I)%VAL=ZERO
                            if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                            Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))  ;U_CT(I)%VAL=ZERO    
                            Endif
                            ELSE
                            ALLOCATE (U_C(I)%VAL(3,5));U_C(I)%VAL=ZERO
                            if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                            Allocate(U_CT(I)%VAL(3,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
                            Endif

			    END IF
                            IF (ITESTCASE.LE.3)THEN
                            ALLOCATE (U_E(I)%VAL(1,1));U_E(I)%VAL=ZERO
                            end if
                        END IF
                    if (rungekutta.eq.4)then
                    IF (AVERAGING.EQ.1)THEN
                    ALLOCATE (U_C(I)%VAL(7,5));U_C(I)%VAL=ZERO
                    if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                        Allocate(U_CT(I)%VAL(7,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
                    Endif
                    else
                    ALLOCATE (U_C(I)%VAL(6,5));U_C(I)%VAL=ZERO
                    if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                        Allocate(U_CT(I)%VAL(6,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
                    Endif
                        IF (ITESTCASE.LE.3)THEN
                        ALLOCATE (U_E(I)%VAL(1,1));U_E(I)%VAL=ZERO
                        end if
		      end if
		      
		      
		      
			end if
		      
			
		      
		      
		      
		      IF (RUNGEKUTTA.EQ.1)THEN
			ALLOCATE (U_C(I)%VAL(1,5));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(1,turbulenceequations+PASSIVESCALAR))  ;U_CT(I)%VAL=ZERO    
			Endif
		      END IF
! 		      IF (RUNGEKUTTA.EQ.4)THEN
! 			ALLOCATE (U_C(I)%VAL(6,5));U_C(I)%VAL=ZERO
! 			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
! 			    Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR)) ;U_CT(I)%VAL=ZERO     
! 			Endif
! 		      END IF
		      IF ((RUNGEKUTTA.EQ.2).or.(rungekutta.eq.5))THEN
			ALLOCATE (U_C(I)%VAL(2,nof_Variables));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(2,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
		      END IF


		IF (AVERAGING.EQ.1)THEN
		ALLOCATE(U_C(I)%RMS(7))
		U_C(I)%RMS(:)=ZERO
		END IF

		END IF
		
	END DO
	ELSE
	
	DO I=1,KMAXE
		IF (ITESTCASE.LT.3)THEN
			
			if (rungekutta.eq.3)then
			ALLOCATE (U_C(I)%VAL(3,1));U_C(I)%VAL=ZERO
			end if
			if (rungekutta.eq.1)then
			ALLOCATE (U_C(I)%VAL(1,1));U_C(I)%VAL=ZERO
			end if
			if (rungekutta.eq.4)then
			ALLOCATE (U_C(I)%VAL(6,1));U_C(I)%VAL=ZERO
			end if
! 			ALLOCATE (U_E(I)%VAL(1,1));U_E(I)%VAL=ZERO

			
		END IF
		
		IF (ITESTCASE.le.3)then
		ALLOCATE (U_E(I)%VAL(1,nof_Variables));U_E(I)%VAL=ZERO
		end if
		
		
		IF (ITESTCASE.GE.3)THEN
		
		
		
		
		      IF( (RUNGEKUTTA.EQ.10))THEN
			ALLOCATE (U_C(I)%VAL(1,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(1,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
		      END IF
		      IF (RUNGEKUTTA.EQ.8)THEN
			ALLOCATE (U_C(I)%VAL(5,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
		      END IF
		      IF( (RUNGEKUTTA.EQ.11))THEN
				IF (AVERAGING.EQ.1)THEN
			ALLOCATE (U_C(I)%VAL(5,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))  ;U_CT(I)%VAL=ZERO    
			Endif
			ELSE
				ALLOCATE (U_C(I)%VAL(3,4));U_C(I)%VAL=ZERO
				if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
				    Allocate(U_CT(I)%VAL(3,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
				Endif
			      end if
		      END IF
		      IF (RUNGEKUTTA.EQ.3)THEN
			IF (AVERAGING.EQ.1)THEN
			ALLOCATE (U_C(I)%VAL(5,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
			ELSE
			ALLOCATE (U_C(I)%VAL(3,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(3,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
			


			END IF
		      END IF
		      IF (RUNGEKUTTA.EQ.1)THEN
			ALLOCATE (U_C(I)%VAL(1,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(1,turbulenceequations+PASSIVESCALAR))    ;U_CT(I)%VAL=ZERO  
			Endif
		      END IF
! 		      IF (RUNGEKUTTA.EQ.4)THEN
! 			ALLOCATE (U_C(I)%VAL(5,4));U_C(I)%VAL=ZERO
! 			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
! 			    Allocate(U_CT(I)%VAL(5,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
! 			Endif
! 		      END IF
		       IF ((RUNGEKUTTA.EQ.2).or.(rungekutta.eq.5))THEN
			ALLOCATE (U_C(I)%VAL(2,nof_Variables));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(2,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
		      END IF

		      if (rungekutta.eq.4)then
			ALLOCATE (U_C(I)%VAL(6,4));U_C(I)%VAL=ZERO
			if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
			    Allocate(U_CT(I)%VAL(6,turbulenceequations+PASSIVESCALAR))   ;U_CT(I)%VAL=ZERO   
			Endif
			end if
			
		IF (AVERAGING.EQ.1)THEN
		ALLOCATE(U_C(I)%RMS(4))
		U_C(I)%RMS(:)=ZERO
		END IF

		END IF
		
	END DO
	
	
	
	END IF
	
	
	END SUBROUTINE U_C_ALLOCATION

	

	
subroutine local_reconallocation5(n)
   !> @brief
!> This subroutine allocates memory for reconstruction (one per process since these are destroyed after each element)
implicit none
INTEGER,INTENT(IN)::N

if (dimensiona.eq.3)then
allocate (ilocal_recon5(1))
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTS(TYPESTEN,IDEGFREE,NOF_VARIABLES)) !1000
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSc(TYPESTEN,IDEGFREE2,NOF_VARIABLES)) !1000
end if

if ((turbulenceequations.gt.0).or.(PASSIVESCALAR.gt.0))then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTS2(TYPESTEN,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSC2(TYPESTEN,IDEGFREE2,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
end if
end if
if (itestcase.eq.4)Then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTURB(1,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR)) !10
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTEMP(IDEGFREE))!10
ALLOCATE (ILOCAL_RECON5(1)%VELOCITYDOF(3,IDEGFREE))!30


ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTURB_wall(6*numberofpoints2,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR)) !10
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTEMP_wall(6*numberofpoints2,IDEGFREE))!10
ALLOCATE (ILOCAL_RECON5(1)%VELOCITYDOF_wall(6*numberofpoints2,3,IDEGFREE))!30
end if



else

allocate (ilocal_recon5(1))
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTS(TYPESTEN,IDEGFREE,NOF_VARIABLES)) !1000
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSc(TYPESTEN,IDEGFREE2,NOF_VARIABLES)) !1000
end if
if ((turbulenceequations.gt.0).or.(PASSIVESCALAR.gt.0))then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTS2(TYPESTEN,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSC2(TYPESTEN,IDEGFREE2,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
end if
end if





if (itestcase.eq.4)Then
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTURB(1,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR)) !10
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTEMP(IDEGFREE))!10
ALLOCATE (ILOCAL_RECON5(1)%VELOCITYDOF(2,IDEGFREE))!30


ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTURB_wall(4*numberofpoints2,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR)) !10
ALLOCATE (ILOCAL_RECON5(1)%GRADIENTSTEMP_wall(4*numberofpoints2,IDEGFREE))!10
ALLOCATE (ILOCAL_RECON5(1)%VELOCITYDOF_wall(4*numberofpoints2,2,IDEGFREE))!30
end if

end if
end subroutine 

SUBROUTINE LOCAL_RECONALLOCATION4(N)
   !> @brief
!> This subroutine allocates memory for reconstruction 
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::K,I,J,L,M,IT,KMAXE,IDUM,ICCF,decomf,SVG,points,ii
INTEGER::Q5,Q4,Q3,Q2,Q1,Q0,q01,icnn
REAL::PERC,PERDE,PERDI,PER1,PER2,PER3,PER4,PER5,PER0,PEF0,PEF1,PEF2,PEF3,PEF4,PEF5,PERV,per01,pef01
KMAXE=XMPIELRANK(N)
IF (ITESTCASE.GE.3) THEN
	IT=5
END IF
IF (ITESTCASE.LT.3) THEN
	IT=1
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)



DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I

	select case(ielem(n,i)%ishape)
	
	case(1,3,4)
	points=QP_quad_n
	
	case(2)
	points=QP_TRIANGLE_n
	end select
		
	IF (ITESTCASE.EQ.4)THEN
		
	if (fastest.ne.1)then
	ALLOCATE (ILOCAL_RECON3(I)%ULEFTV(dims,IT-1,ielem(n,i)%ifca,points))
	
	
	
	else
	ALLOCATE (ILOCAL_RECON3(I)%ULEFTV(dims,IT-1,ielem(n,i)%ifca,1))

	end if
	ILOCAL_RECON3(I)%ULEFTV=zero
	
	if ((turbulence.eq.1).or.(passivescalar.gt.0)) then
	  
	  SVG=(TURBULENCEEQUATIONS+PASSIVESCALAR)
	  if (fastest.ne.1)then
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURBV(dims,svg,ielem(n,i)%ifca,points))	! THE DERIVATIVES OF THE TURBULENCE MODEL
	 
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURB(TURBULENCEEQUATIONS+PASSIVESCALAR,ielem(n,i)%ifca,points))

	  else
	   ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURBV(dims,svg,ielem(n,i)%ifca,1))	! THE DERIVATIVES OF THE TURBULENCE MODEL
	 
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURB(TURBULENCEEQUATIONS+PASSIVESCALAR,ielem(n,i)%ifca,1))



	  end if
	  ILOCAL_RECON3(I)%ULEFTTURB=zero;ILOCAL_RECON3(I)%ULEFTTURBv=zero
	  
	END IF
	
	
	
	 
	END IF
	ALLOCATE (ILOCAL_RECON3(I)%GRADS(4+TURBULENCEEQUATIONS+passivescalar+(QSAS_MODEL*3),3))
	
	IF ((OUTSURF.EQ.1).AND.(AVERAGING.EQ.1))THEN
	ALLOCATE (ILOCAL_RECON3(I)%GRADSAV(4,3))
	END IF
	
	 if (fastest.ne.1)then
        ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,points))
	else
	 ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,1))

	  end if
	  ILOCAL_RECON3(I)%ULEFT=zero
	
! 	END IF
	
	
	
END DO



DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I
	
	select case(ielem(n,i)%ishape)
	
	case(1,3,4)
	points=QP_quad_n
	
	case(2)
	points=QP_TRIANGLE_n
	end select
		
	IF (ITESTCASE.EQ.4)THEN
		
	if (fastest.ne.1)then
	ALLOCATE (ILOCAL_RECON3(I)%ULEFTV(dims,IT-1,ielem(n,i)%ifca,points))
	
	
	
	else
	ALLOCATE (ILOCAL_RECON3(I)%ULEFTV(dims,IT-1,ielem(n,i)%ifca,1))

	end if
	ILOCAL_RECON3(I)%ULEFTV=zero
	
	if ((turbulence.eq.1).or.(passivescalar.gt.0)) then
	  
	  SVG=(TURBULENCEEQUATIONS+PASSIVESCALAR)
	  if (fastest.ne.1)then
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURBV(dims,svg,ielem(n,i)%ifca,points))	! THE DERIVATIVES OF THE TURBULENCE MODEL
	 
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURB(TURBULENCEEQUATIONS+PASSIVESCALAR,ielem(n,i)%ifca,points))

	  else
	   ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURBV(dims,svg,ielem(n,i)%ifca,1))	! THE DERIVATIVES OF THE TURBULENCE MODEL
	 
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURB(TURBULENCEEQUATIONS+PASSIVESCALAR,ielem(n,i)%ifca,1))



	  end if
	  ILOCAL_RECON3(I)%ULEFTTURB=zero;ILOCAL_RECON3(I)%ULEFTTURBv=zero
	  
	END IF
	
	
	
	 
	END IF
	ALLOCATE (ILOCAL_RECON3(I)%GRADS(4+TURBULENCEEQUATIONS+passivescalar+(QSAS_MODEL*3),3))
	IF ((OUTSURF.EQ.1).AND.(AVERAGING.EQ.1))THEN
	ALLOCATE (ILOCAL_RECON3(I)%GRADSAV(4,3))
	END IF
	
	 if (fastest.ne.1)then
        ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,points))
	else
	 ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,1))

	  end if
	  ILOCAL_RECON3(I)%ULEFT=zero
	
	
	
END DO


END SUBROUTINE LOCAL_RECONALLOCATION4


SUBROUTINE LOCAL_RECONALLOCATION42d(N)
   !> @brief
!> This subroutine allocates memory for reconstruction in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::K,I,J,L,M,IT,KMAXE,IDUM,ICCF,decomf,SVG,points
INTEGER::Q5,Q4,Q3,Q2,Q1,Q0,q01,icnn
REAL::PERC,PERDE,PERDI,PER1,PER2,PER3,PER4,PER5,PER0,PEF0,PEF1,PEF2,PEF3,PEF4,PEF5,PERV,per01,pef01
KMAXE=XMPIELRANK(N)
IF (ITESTCASE.GE.3) THEN
	IT=5
END IF
IF (ITESTCASE.LT.3) THEN
	IT=1
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




DO I=1,KMAXE
	
	
	
	
	points=qp_line_n
	
	
	
		
	IF (ITESTCASE.EQ.4)THEN
		
	if (fastest.ne.1)then
	ALLOCATE (ILOCAL_RECON3(I)%ULEFTV(dims,IT-1,ielem(n,i)%ifca,points))
	else
	ALLOCATE (ILOCAL_RECON3(I)%ULEFTV(dims,IT-1,ielem(n,i)%ifca,1))

	end if
	ILOCAL_RECON3(I)%ULEFTV=zero
	
	if ((turbulence.eq.1).or.(passivescalar.gt.0)) then
	  
	  SVG=(TURBULENCEEQUATIONS+PASSIVESCALAR)
	  if (fastest.ne.1)then
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURBV(dims,svg,ielem(n,i)%ifca,points))	! THE DERIVATIVES OF THE TURBULENCE MODEL
	 
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURB(TURBULENCEEQUATIONS+PASSIVESCALAR,ielem(n,i)%ifca,points))

	  else
	   ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURBV(dims,svg,ielem(n,i)%ifca,1))	! THE DERIVATIVES OF THE TURBULENCE MODEL
	 
	  ALLOCATE (ILOCAL_RECON3(I)%ULEFTTURB(TURBULENCEEQUATIONS+PASSIVESCALAR,ielem(n,i)%ifca,1))



	  end if
	  ILOCAL_RECON3(I)%ULEFTTURB=zero;ILOCAL_RECON3(I)%ULEFTTURBv=zero
	  
	END IF
	
	
	ALLOCATE (ILOCAL_RECON3(I)%GRADS(3+TURBULENCEEQUATIONS+passivescalar+(QSAS_MODEL*2),2))
	 
	END IF
	
	
	
	 if (fastest.ne.1)then
        ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,points))
	else
	 ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,1))

	  end if
	  ILOCAL_RECON3(I)%ULEFT=zero
	
! 	END IF
	
	
	
END DO


END SUBROUTINE LOCAL_RECONALLOCATION42d
	

SUBROUTINE MEMORY2
implicit none
   !> @brief
!> This subroutine allocates memory for all the matrices
INTEGER::KKD

KKD=nof_variables

!sorting them out first in terms of dimensions, no of variables, turbulence parameters
ALLOCATE(GRADSV(1:dims),VORTEM1(1:dims),VERTF1(1:dims),vertex1(1:dims),DERIVTURB1(1:dims),UGRADLOC(1:dims),&
GRADVEL(1:IDEGFREE),GRADTEM(1:IDEGFREE),&
GRADSS(1:IDEGFREE,NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR),VORTET1(1:dims,1:dims),consmatrix(1:NUMBEROFPOINTS2*6,1:idegfree),consmatrixc(1:NUMBEROFPOINTS2*6,1:idegfree),LAMC(TYPESTEN),GRAD3AL(IDEGFREE),GRAD5ALc(IDEGFREE,nof_variables))



ALLOCATE(VECL(kkd),VECR(kkd),CFAST(dims),LSOLP(kkd),RSOLP(kkd),SUMGGU(1:dims),SUMGGV(1:dims),SUMGGW(1:dims),SUMGGP(1:dims),SUMGGT(1:dims),&
SUMGGT2(1:dims),MATRIX_BN3(IMAXDEGFREE),MATRIX_BN4(IMAXDEGFREE),MATRIX_BN5(IMAXDEGFREE),MATRIX_B1(IMAXDEGFREE),MATRIX_B2(IMAXDEGFREE),MATRIX_B3(IMAXDEGFREE),MATRIX_B4(IMAXDEGFREE)&
,MATRIX_B5(IMAXDEGFREE),Xxx(1:IDEGFREE),X2(1:IDEGFREE-1),MATRIX_BWC(1:IDEGFREE-1),MATRIX_BN(1:IDEGFREE),MATRIX_BN1(IMAXDEGFREE),MATRIX_BN2(IMAXDEGFREE)&
,RFFAST(1:dims,1:dims),RRFFAST(1:dims,1:dims),tempsc(1:dims,1:dims))


ALLOCATE(USOL(NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2),PSI(NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2),UTEMP(IMAXDEGFREE+1,1:kkd),UTMIN(NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR),UTMAX(NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR))


ALLOCATE(GRAD1(1:kkd,0:IDEGFREE),GRAD2(1:kkd,0:IDEGFREE),INDICATEMATRIX(1:kkd,1:IDEGFREE),GRADCHAR(1:kkd,0:IDEGFREE),LIMITEDDW(1:kkd,0:IDEGFREE)&
,WENOO(1:kkd,1:TYPESTEN),WENO(1:kkd+TURBULENCEEQUATIONS+PASSIVESCALAR,1:TYPESTEN),WENO2(1:kkd,1:TYPESTEN),GRADSSL(1:IDEGFREE,1:NOF_VARIABLES),GRADCHARV(1:kkd,1:TYPESTEN,0:IDEGFREE),FINDW(kkd,0:IDEGFREE,6,2),FINDW_CHAR(kkd,0:IDEGFREE,6,2,1:TYPESTEN),LIMITEDDW_CHAR(1:kkd,0:IDEGFREE,1:TYPESTEN),&
SMOOTHINDICATOR(1:kkd,1:TYPESTEN,1:6,1:2),LAMBDA(1:kkd,1:TYPESTEN,1:6,1:2),OMEGATILDE(1:kkd,1:TYPESTEN,1:6,1:2),OMEGA(1:kkd,1:TYPESTEN,1:6,1:2),WENOOS(1:kkd,1:TYPESTEN,1:6,1:2)&
,GRAD1AL(1:IDEGFREE),INDICATEMATRIXAL(1:IDEGFREE),SMOOTHINDICATORAL(1:TYPESTEN),LAMBDAAL(1:TYPESTEN),OMEGAATILDEL(1:TYPESTEN),&
OMEGAAL(1:TYPESTEN),SUMOMEGATILDE(1:kkd),RESSOLUTION(1:NUMBEROFPOINTS2*6,1:NOF_VARIABLES),A_CHAR(1:IDEGFREE,1:NOF_VARIABLES,1:TYPESTEN),B_CHAR(1:IDEGFREE,1:IDEGFREE),X_cHAR(1:kkd,1:TYPESTEN))

 
if (rungekutta.ge.10)then
ALLOCATE(B1_imp(nof_Variables),DU1(nof_Variables),DU2(nof_Variables),DUMMY12(nof_Variables),C1_imp(nof_Variables),lscqm1(nof_Variables,nof_Variables))
ALLOCATE(DUR(nof_Variables),DUL(nof_Variables))
ALLOCATE(DURR(nof_Variables),DULR(nof_Variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
ALLOCATE(DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(B1T(TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR))
end if

end if



END SUBROUTINE MEMORY2









END MODULE MEMORY
