MODULE MEMORY
USE MPIINFO
USE DECLARATION
IMPLICIT NONE

CONTAINS



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
ALLOCATE(ILOCALALLELGPER(N:N,xmpielrank(n),1,ISELEMT(N)))
ILOCALALLELGPER(:,:,:,:)=0
end subroutine




SUBROUTINE ALLOCATE3
IMPLICIT NONE
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

SUBROUTINE QUADALLOC(NUMBEROFPOINTS,NUMBEROFPOINTS2)
!> @brief
!> This subroutine allocates memory for the quadrature points
implicit none
INTEGER,INTENT(INout)::NUMBEROFPOINTS,NUMBEROFPOINTS2
integer::i,kmaxe


IF (DIMENSIONA.EQ.3)THEN
    NUMBEROFPOINTS=MAX(QP_HEXA,QP_TETRA,QP_PYRA,QP_PRISM)
    if (dg.eq.1)NUMBEROFPOINTS=MAX(QP_HEXA,QP_TETRA*6,QP_PYRA,QP_PRISM)
    NUMBEROFPOINTS2=MAX(QP_QUAD,QP_TRIANGLE,QP_TRIANGLE)
    

ELSE
    NUMBEROFPOINTS=MAX(QP_QUAD,QP_TRIANGLE)
    if (dg.eq.1)NUMBEROFPOINTS=MAX(QP_QUAD,QP_TRIANGLE*2)
    NUMBEROFPOINTS2=QP_LINE
   
END IF



kmaxe=xmpielrank(n)
do i=1,kmaxe

select case (ielem(n,i)%ishape)
        
        
        case(1) !hexa
        ielem(n,i)%iTOTALPOINTS=QP_Tetra*6
        
        case(2) !tetra
        ielem(n,i)%iTOTALPOINTS=QP_Tetra
        
        case(3) !pyramid
        ielem(n,i)%iTOTALPOINTS=QP_Tetra*2
        
        case(4) !prism
        ielem(n,i)%iTOTALPOINTS=QP_tetra*3
        
        case(5) !quadrilateral
        ielem(n,i)%iTOTALPOINTS=QP_TRIANGLE*2
        
        
        
        case(6)!triangle
         ielem(n,i)%iTOTALPOINTS=QP_TRIANGLE
         
         
         end select


end do




END SUBROUTINE QUADALLOC





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


	
	
	IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	ALLOCATE (RHST(KMAXE))
	END IF
	
	DO I=1,KMAXE
        IF (DG == 1) THEN
            ALLOCATE(RHS(I)%VALDG(NUM_DG_DOFS, NOF_VARIABLES))

            
        end if
            ALLOCATE (RHS(I)%VAL(nof_Variables))
            
            IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)) THEN
                ALLOCATE (RHST(I)%VAL(TURBULENCEEQUATIONS+PASSIVESCALAR))
            END IF
        
    END DO
	
END SUBROUTINE SUMFLUX_ALLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE IMPALLOCATE(N)
!> @brief
!> This subroutine allocates memory for implicit time stepping
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,interf
KMAXE=XMPIELRANK(N)


if (dimensiona.eq.3)then
interf=nof_Variables
else
interf=nof_Variables
end if

IF (RUNGEKUTTA.EQ.12)THEN
ALLOCATE (IMPdu(KMAXE,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR))
impdu(:,:)=zero
ELSE



IF (RELAX.EQ.3)THEN

ALLOCATE (IMPDIAG_MF(KMAXE))
ALLOCATE (IMPOFF_MF(KMAXE,INTERF))
ALLOCATE (IMPdu(KMAXE,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR))

IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPDIAGT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPOFFt(KMAXE,INTERF,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF

IMPDIAG_MF=zero
IMPOFF_MF=zero
impdu=zero

ELSE




if (dimensiona.eq.3)then
if (lowmemory.eq.0)then

ALLOCATE (IMPDIAG(KMAXE,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPOFF(KMAXE,6,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPdu(KMAXE,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR))

IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(KMAXE,6,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF


IMPDIAG(:,:,:)=zero
IMPOFF(:,:,:,:)=zero
impdu(:,:)=zero

else

ALLOCATE (IMPdu(KMAXE,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR))
impdu(:,:)=zero
IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF
end if

else

if (lowmemory.eq.0)then

ALLOCATE (IMPDIAG(KMAXE,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPOFF(KMAXE,4,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPdu(KMAXE,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR))

IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(KMAXE,4,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF


IMPDIAG(:,:,:)=zero
IMPOFF(:,:,:,:)=zero
impdu(:,:)=zero

else

! ALLOCATE (IMPDIAG(1,1:nof_Variables,1:nof_Variables))
! ALLOCATE (IMPOFF(1,4,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPdu(KMAXE,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR))
IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
! ALLOCATE(IMPOFFt(1,4,TURBULENCEEQUATIONS+PASSIVESCALAR))
! ALLOCATE(IMPDIAGT(1,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(SHT(KMAXE,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF
! IMPDIAG(1,:,:)=zero
! IMPOFF(1,:,:,:)=zero
impdu(:,:)=zero
end if



end if
end if
END IF




END  SUBROUTINE IMPALLOCATE



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
	allocate (nodes_offset(imaxe),nodes_offset2(imaxe))
	IESHAPE=0
	nodes_offset=0
	nodes_offset2=0
	
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
	DEALLOCATE (IESHAPE,nodes_offset,nodes_offset2)
	END SUBROUTINE SHDEALLOCATION





!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR ELEMENTS!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE ELALLOCATION(N,XMPIE,XMPIELRANK,IELEM,IMAXE,IESHAPE,ITESTCASE,IMAXB,IBOUND,XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX)
	IMPLICIT NONE
	!> @brief
!> This subroutine allocates memory for the elements
	TYPE(A_ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::IESHAPE
	INTEGER,INTENT(IN)::N
	INTEGER,INTENT(IN)::IMAXE,ITESTCASE,IMAXB
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	TYPE(A_BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:)::IBOUND
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
	IF (TECPLOT.EQ.5)THEN
	allocate(nodes_offset_local(1:kmaxe));nodes_offset_local=0
	allocate(nodes_offset_local2(1:kmaxe));nodes_offset_local2=0
	END IF
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
	TYPE(A_NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
	TYPE(A_NODE_NE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INODEN
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
	TYPE(A_NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
	TYPE(A_NODE_NE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INODEN
	INTEGER,INTENT(IN)::IMAXN
 	DEALLOCATE (INODE)
! 	DEALLOCATE (INODEN)
	END SUBROUTINE NODEDEALLOCATION
!---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







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

END SUBROUTINE XMPIALLOCATE


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






SUBROUTINE LOCAL_RECONALLOCATION3(N,ILOCAL_RECON3)
!> @brief
!> This subroutine allocates memory for reconstruction prestoring
IMPLICIT NONE
TYPE(A_LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,M,IKG,ITRUE,kmaxe,idum
INTEGER::inum_points,ITARGET
INTEGER::imax,inum,ideg,imax2,inum2,ideg2,m2
REAL::PERC,PERDE,PERD,PER1,PER2,PER3,PER4,PER5,PER0,PEF0,PEF1,PEF2,PEF3,PEF4,PEF5,PERV,per01,pef01
KMAXE=XMPIELRANK(N)

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

ALLOCATE (ILOCAL_RECON3(KMAXE))

IF (DG.EQ.1)THEN
ALLOCATE (ILOCAL_RECON6(KMAXE))
DO I=1,KMAXE	!for all elements
 ALLOCATE (ILOCAL_RECON6(I)%DG2FV(1:IDEGFREE,1:NOF_vARIABLES))
 ILOCAL_RECON6(I)%DG2FV=0.0D0
END DO

END IF

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
	   ALLOCATE (ILOCAL_RECON3(I)%VOLUME(1,INUM));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
	   end if
		
		
	    
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
				   ALLOCATE (ILOCAL_RECON3(I)%WEIGHTL(M,IMAX));ILOCAL_RECON3(I)%WEIGHTL(:,:)=0.0d0!WEIGHTL
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
	ALLOCATE (ILOCAL_RECON3(I)%PERIODICFLAG(M,INUM))
	
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
	   ALLOCATE (ILOCAL_RECON3(I)%VOLUME(1,INUM));ILOCAL_RECON3(I)%VOLUME(:,:)=0.0D0
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
	   ALLOCATE (ILOCAL_RECON3(I)%WEIGHTL(M,IMAX));ILOCAL_RECON3(I)%WEIGHTL(:,:)=0.0d0
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
!       perc=ilocal_recon3(i)%local
!       perd=kmaxe
!       perde=perde+perc

    
	
END DO
END IF



	



END SUBROUTINE LOCAL_RECONALLOCATION3



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DEALCORDINATES1(N,IEXCORDR,IEXCORDS)
!> @brief
!> This subroutine deallocates memory for exchange of info between processes
IMPLICIT NONE
TYPE(A_EXCHANGE_CORD),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXCORDR
TYPE(A_EXCHANGE_CORD),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXCORDS
INTEGER,INTENT(IN)::N
DEALLOCATE(IEXCORDR)


END SUBROUTINE DEALCORDINATES1


SUBROUTINE DEALCORDINATES2
IMPLICIT NONE
!> @brief
!> This subroutine deallocates memory for exchange of info between processes
 if (allocated(iexcords))DEALLOCATE(IEXCORDS)

END SUBROUTINE DEALCORDINATES2





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ALLOCATE_BASIS_FUNCTION(N,INTEG_BASIS,XMPIELRANK,IDEGFREE)
!> @brief
!> This subroutine allocates memory for basis function integrals
IMPLICIT NONE
TYPE(A_INTEGRALBASIS),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INTEG_BASIs
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER,INTENT(IN)::IDEGFREE,N
INTEGER::KMAXE,i
KMAXE=XMPIELRANK(N)

ALLOCATE(INTEG_BASIS(KMAXE))
if (dg.eq.1)then
ALLOCATE(INTEG_BASIS_dg(KMAXE))


do i=1,kmaxe
 allocate(INTEG_BASIS_dg(i)%value(1:idegfree));INTEG_BASIS_dg(i)%value(:)=zero
end do
end if
do i=1,kmaxe
 allocate(INTEG_BASIS(i)%value(1:idegfree));INTEG_BASIS(i)%value(:)=zero
 if (ees.eq.5)then
 allocate(INTEG_BASIS(i)%valuec(1:idegfree2));INTEG_BASIS(i)%valuec(:)=zero
 end if
		  
end do

END SUBROUTINE ALLOCATE_BASIS_FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DEALLOCATE_BASIS_FUNCTION(N,INTEG_BASIS)
!> @brief
!> This subroutine deallocates memory for basis function integrals
IMPLICIT NONE
TYPE(A_INTEGRALBASIS),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INTEG_BASIS
INTEGER,INTENT(IN)::N
DEALLOCATE(INTEG_BASIS)
END SUBROUTINE DEALLOCATE_BASIS_FUNCTION




   SUBROUTINE LOCALSDEALLOCATION(N,XMPIELRANK,ILOCALSTENCIL,ILOCALSTENCILPER,TYPESTEN,NUMNEIGHBOURS)
   !> @brief
!> This subroutine allocates memory for stencils
	IMPLICIT NONE
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),intent(inout)::ILOCALSTENCIL,ILOCALSTENCILPER
	INTEGER,INTENT(IN)::NUMNEIGHBOURS
	INTEGER,INTENT(IN)::N
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,INTENT(IN)::TYPESTEN
	
	DEALLOCATE (ILOCALSTENCIL)
	DEALLOCATE (ILOCALSTENCILPER)
	END SUBROUTINE LOCALSDEALLOCATION
	
	

	
	SUBROUTINE U_C_ALLOCATION(N,XMPIELRANK,U_C,U_E,ITESTCASE,U_CT)
	   !> @brief
!> This subroutine allocates memory for solution vector
	IMPLICIT NONE
	TYPE(A_U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C,U_CT	
	TYPE(A_U_EXACT),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_E
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,INTENT(IN)::ITESTCASE,N
	INTEGER::I,KMAXE,ISTAGE,TOTALPOINTS,totalpointsvol
	KMAXE=XMPIELRANK(N)
	ALLOCATE (U_C(KMAXE))
	IF (FILTERING.EQ.1)THEN
	ALLOCATE (U_CW(KMAXE))
	ALLOCATE (U_CS(KMAXE))
	END IF
	
	
	
	IF (ITESTCASE.LE.4)THEN
	ALLOCATE (U_E(KMAXE))
	end if
	
	if (( turbulence .GT. 0).OR.(PASSIVESCALAR.GT.0))THEN
	  Allocate(U_CT(kmaxe))
	END IF

	SELECT CASE(RUNGEKUTTA)
	
	CASE(1)
	ISTAGE=1
	
	CASE(2)
	ISTAGE=2
	
	CASE(3)
	IF (AVERAGING.EQ.1)THEN
	ISTAGE=5
	ELSE
	ISTAGE=3
	
	if (mood.eq.1)then
	ISTAGE=4
	end if
	
	END IF
	
	
	
	
	CASE(4)
	IF (AVERAGING.EQ.1)THEN
	ISTAGE=7
	ELSE
	ISTAGE=6
	END IF
	
	CASE(5)
	ISTAGE=2
	
	CASE(10)
	ISTAGE=1
	
	CASE(11)
	
	IF (AVERAGING.EQ.1)THEN
	ISTAGE=5
	ELSE
	ISTAGE=3
	END IF
	
	CASE(12)
	
	IF (AVERAGING.EQ.1)THEN
	ISTAGE=5
	ELSE
	ISTAGE=3
	END IF
	
	END SELECT
	
	
	if (DG.EQ.1)THEN
	
        allocate(M_1(kmaxe))
	END IF
	
	DO I=1,KMAXE
        ALLOCATE (U_C(I)%VAL(ISTAGE,NOF_VARIABLES));U_C(I)%VAL=ZERO
        IF (FILTERING.EQ.1)tHEN
        ALLOCATE (U_CW(I)%VAL(1,NOF_VARIABLES));U_CW(I)%VAL=ZERO
        ALLOCATE (U_CS(I)%VAL(1,NOF_VARIABLES));U_CS(I)%VAL=ZERO

        END IF
         
         
         
         
         
         
        
        IF (DG.EQ.1)THEN
            ALLOCATE (U_C(I)%VALDG(ISTAGE,NOF_VARIABLES,IELEM(N,I)%IDEGFREE+1));U_C(I)%VALDG=ZERO
            IF (FILTERING.EQ.1)tHEN
            ALLOCATE (U_CW(I)%VALDG(1,NOF_VARIABLES,IELEM(N,I)%IDEGFREE+1));U_CW(I)%VALDG=ZERO
            ALLOCATE (U_CS(I)%VALDG(1,NOF_VARIABLES,IELEM(N,I)%IDEGFREE+1));U_CS(I)%VALDG=ZERO
            END IF
            
            allocate (M_1(i)%val(1:idegfree+1,1:idegfree+1));M_1(i)%val=zero
            
            IF (ITESTCASE == 4) ALLOCATE(U_C(I)%BR2_AUX_VAR(IELEM(N,I)%IDEGFREE+1,NOF_VARIABLES,DIMENSIONA)) ! NS
            
        END IF
        
        
                    if (( turbulence .eq. 1).or.(PASSIVESCALAR.GT.0))THEN
                        Allocate(U_CT(I)%VAL(ISTAGE,turbulenceequations+PASSIVESCALAR));U_CT(I)%VAL=ZERO   
                    Endif
        IF (AVERAGING.EQ.1)THEN
		ALLOCATE(U_C(I)%RMS(7))
		U_C(I)%RMS(:)=ZERO
		END IF
		IF (ITESTCASE.LE.4)THEN
		ALLOCATE (U_E(I)%VAL(1,NOF_VARIABLES));U_E(I)%VAL=ZERO
		END IF
	END DO
	
	
	IF (MOOD.EQ.1)THEN
DO I=1,KMAXE
    IELEM(N,I)%RECALC=0
END DO
ELSE
DO I=1,KMAXE
    IELEM(N,I)%RECALC=1
END DO
END IF
	
	
	
	
	
	END SUBROUTINE U_C_ALLOCATION

	

	
subroutine local_reconallocation5(n)
   !> @brief
!> This subroutine allocates memory for reconstruction (one per process since these are destroyed after each element)
implicit none
INTEGER,INTENT(IN)::N
integer::KMAXE,I
kmaxe=xmpielrank(n)

allocate (ilocal_recon5(1:kmaxe))



do i=1,kmaxe

if (dimensiona.eq.3)then


ALLOCATE (ILOCAL_RECON5(i)%GRADIENTS(TYPESTEN,IDEGFREE,NOF_VARIABLES)) !1000
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSc(TYPESTEN,IDEGFREE2,NOF_VARIABLES)) !1000
end if



		if ((turbulenceequations.gt.0).or.(PASSIVESCALAR.gt.0))then
		ALLOCATE (ILOCAL_RECON5(i)%GRADIENTS2(TYPESTEN,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
				if (ees.eq.5)then
				ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSC2(TYPESTEN,IDEGFREE2,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
				end if
		end if


if (itestcase.eq.4)Then
if ((turbulenceequations.gt.0).or.(PASSIVESCALAR.gt.0))then
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSTURB(1,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR)) !10
END IF
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSTEMP(IDEGFREE))!10
ALLOCATE (ILOCAL_RECON5(i)%VELOCITYDOF(3,IDEGFREE))!30

end if



else


ALLOCATE (ILOCAL_RECON5(i)%GRADIENTS(TYPESTEN,IDEGFREE,NOF_VARIABLES)) !1000
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSc(TYPESTEN,IDEGFREE2,NOF_VARIABLES)) !1000
end if
if ((turbulenceequations.gt.0).or.(PASSIVESCALAR.gt.0))then
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTS2(TYPESTEN,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
if (ees.eq.5)then
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSC2(TYPESTEN,IDEGFREE2,0+TURBULENCEEQUATIONS+PASSIVESCALAR))!20
end if
end if





if (itestcase.eq.4)Then
if ((turbulenceequations.gt.0).or.(PASSIVESCALAR.gt.0))then
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSTURB(1,IDEGFREE,0+TURBULENCEEQUATIONS+PASSIVESCALAR)) !10
END IF
ALLOCATE (ILOCAL_RECON5(i)%GRADIENTSTEMP(IDEGFREE))!10
ALLOCATE (ILOCAL_RECON5(i)%VELOCITYDOF(2,IDEGFREE))!30

end if

end if



end do

end subroutine 

SUBROUTINE LOCAL_RECONALLOCATION4(N)
   !> @brief
!> This subroutine allocates memory for reconstruction 
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::K,I,J,L,M,IT,KMAXE,IDUM,ICCF,decomf,SVG,points,ii
INTEGER::Q5,Q4,Q3,Q2,Q1,Q0,q01,icnn,ICONSIDERED
REAL::PERC,PERDE,PERDI,PER1,PER2,PER3,PER4,PER5,PER0,PEF0,PEF1,PEF2,PEF3,PEF4,PEF5,PERV,per01,pef01
KMAXE=XMPIELRANK(N)
IF (ITESTCASE.GE.3) THEN
	IT=nof_variables
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
	IT=nof_variables
END IF
IF (ITESTCASE.LT.3) THEN
	IT=1
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

DO I=1,KMAXE
	
    points=qp_line_n
	
		
	IF (ITESTCASE.EQ.4)THEN !Linear step?
		
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
        IF (CODE_PROFILE.EQ.2)THEN
        ALLOCATE (ILOCAL_RECON3(I)%ULEFTx(IT,ielem(n,i)%ifca,points))
        END IF
    else
        ALLOCATE (ILOCAL_RECON3(I)%ULEFT(IT,ielem(n,i)%ifca,1))

    end if
    ILOCAL_RECON3(I)%ULEFT=zero
	

END DO


END SUBROUTINE LOCAL_RECONALLOCATION42d
	











END MODULE MEMORY
