MODULE RECON
USE DECLARATION
USE DERIVATIVES
USE LIBRARY
USE TRANSFORM
USE LOCAL
USE LAPCK
USE GRADIENTS
USE BASIS
IMPLICIT NONE


 CONTAINS








SUBROUTINE AVERAGE_STRESSES(N)
implicit none
!> @brief
!> Subroutine for calling the computation of the average shear stresses
INTEGER,INTENT(IN)::N
INTEGER::II,I,ICONSIDERED
!$OMP DO
DO II=1,NOF_INTERIOR;I=EL_INT(II);ICONSIDERED=I
      CALL ALLGRADS_INNER_AV(N,I)
END DO
!$OMP END DO 

!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I
	CALL ALLGRADS_MIX_AV(N,I)
END DO	
!$OMP END DO 	
END SUBROUTINE AVERAGE_STRESSES
	
	

SUBROUTINE MEMORY_FAST(N)
!> @brief
!> Subroutine for storing the gaussian quadrature points at the cell interfaces
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX
INTEGER::ICONSIDERED,FACEX,POINTX
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ


KMAXE=XMPIELRANK(N)

if (dimensiona.eq.3)then
    DO I=1,KMAXE
        IF (IELEM(N,I)%ISHAPE.EQ.2)THEN
            ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
			IF (SRFG.EQ.1)THEN
				ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
				ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
			END IF
			IF (MRF.EQ.1)THEN
				ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
				ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
				ALLOCATE (ILOCAL_RECON3(I)%MRF_ORIGIN(1:3))
				ALLOCATE (ILOCAL_RECON3(I)%MRF_VELOCITY(1:3))
			END IF
        ELSE
            ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
			IF (SRFG.EQ.1)THEN
				ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
				ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_QUAD,3))
			END IF
			IF (MRF.EQ.1)THEN
				ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
				ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_QUAD,3))
				ALLOCATE (ILOCAL_RECON3(I)%MRF_ORIGIN(1:3))
				ALLOCATE (ILOCAL_RECON3(I)%MRF_VELOCITY(1:3))
			END IF
        END IF
        ICONSIDERED=I
        DO L=1,IELEM(N,I)%IFCA
            IDUMMY=0
            if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
                    if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
					    IDUMMY=1
					    END IF
                END IF	
                if (ielem(n,i)%types_faces(L).eq.5)then
                iqp=qp_quad
                NND=4
                    IF (IDUMMY.EQ.0)THEN
                        do K=1,nnd
                        VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
!                        VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                        END DO
                    ELSE
                        facex=l;
                        CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)

!                        do K=1,nnd
 !                       VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
 !                       END DO
                    END IF
                call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                else
                iqp=QP_TRIANGLE
                NND=3
                    IF (IDUMMY.EQ.0)THEN
                        do K=1,nnd
                        VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
!                        VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                        END DO
                    ELSE
                        facex=l;
                        CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
 !                       do K=1,nnd
  !                      VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
!                        END DO
                    END IF
                        
                call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                end if
            else
                if (ielem(n,i)%types_faces(L).eq.5)then
                    iqp=qp_quad
                    NND=4
                    do K=1,nnd
                        VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
!                        VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                    END DO 
                    call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                else
                    iqp=QP_TRIANGLE
                    NND=3 
                    do K=1,nnd
                        VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
!                        VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                    END DO  
                    call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                end if
            end if
            
            do NGP=1,iqp			!for gqp
!                ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
				IF (SRFG.EQ.1)THEN
					ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
					POX(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)-SRF_ORIGIN(1:3)
					POY(1:3)=SRF_VELOCITY(1:3)
					ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)=VECT_FUNCTION(POX,POY)
					END IF
					
					IF (MRF.EQ.1)THEN
					ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
					POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
					POY(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)
					ICONSIDERED=I
					FACEX=L
					POINTX=NGP
					CALL MRFSWITCH(N,ICONSIDERED,FACEX,POINTX,pox,poy)
					
					END IF 

		END DO	!NGP
		END DO





		DO L=1,IELEM(N,I)%IFCA
		IDUMMY=0
			    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				    IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					    if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
					    IDUMMY=1
					    END IF
				    END IF	
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    NND=4
                                        IF (IDUMMY.EQ.0)THEN
                                            do K=1,nnd
                                            VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                            VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
            END DO	!NGP
                                        ELSE
                                            facex=l;
                                            CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
                                            do K=1,nnd
                                            VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        END DO
                                        END IF
                                    call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                                    else
                                    iqp=QP_TRIANGLE
                                    NND=3
                                        IF (IDUMMY.EQ.0)THEN
                                            do K=1,nnd
                                            VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                            VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
    END DO
    
ELSE ! 2 dimensions
                                            facex=l;
                                            CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
                                            do K=1,nnd
                                            VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                            END DO
                                        END IF
                                            
                                    call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                                    end if
			    else
					  if (ielem(n,i)%types_faces(L).eq.5)then
					    iqp=qp_quad
					    NND=4
						  do K=1,nnd
						    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
						    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						  END DO 
					    call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  else
					    iqp=QP_TRIANGLE
					    NND=3 
						  do K=1,nnd
						    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
						    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						  END DO  
					    call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  end if
			    end if
		
		do NGP=1,iqp			!for gqp
		ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
		END DO	!NGP
		END DO
		END DO
	  
	  
	  else
    
    DO I=1,KMAXE
        ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,qp_line,2))
        
        ICONSIDERED=I
    
        DO L=1,IELEM(N,I)%IFCA
            IDUMMY=0
        
            if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
                        IDUMMY=1
                    END IF
                END IF
        
                IQP=QP_LINE
                NND=2
                IF (IDUMMY.EQ.0)THEN
                    DO K=1,NND
                        VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                        !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                            VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                        !END IF
                    END DO
                ELSE
                    facex=l;
                    CALL coordinates_face_PERIOD2D1(n,iconsidered,facex,VEXT,NODES_LIST)
                    DO K=1,NND
                        !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                            VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                        !END IF
                    END DO
                END IF
                CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
            ELSE
                IQP=QP_LINE
                NND=2
                DO K=1,NND
                    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                    !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                        VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                    !END IF
                END DO
                CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
            END IF
        
            DO NGP=1,iqp !for gqp
                ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:2)=QPOINTS2D(1:2,NGP) ! Storing surface quadrature points
            END DO !NGP
        END DO
    END DO
END IF
END SUBROUTINE MEMORY_FAST

! ! !---------------------------------------------------------------------------------------------!



SUBROUTINE EXTRAPOLATE_BOUND_LINEAR(USOL,varcons,FACEX,pointx,ICONSIDERED)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces for linear advection equation
IMPLICIT NONE
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED
real,dimension(1:nof_Variables)::leftv
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN)::USOL
REAL::MP_PINFl,gammal

if (WENWRT.EQ.3)THEN
    LEFTV(1:NOF_VARIABLES)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
    CALL cons2prim(N,leftv,MP_PINFl,gammal)
    LEFTV(1:NOF_VARIABLES)=LEFTV(1:NOF_VARIABLES)+USOL(1:nof_Variables,FACEX,pointx)
    CALL PRIM2CONS(N,LEFTV)
    ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)+LEFTV(1:NOF_VARIABLES)
ELSE
    ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)&
    +(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)))
END IF


IF (TURBULENCEEQUATIONS.GE.1)THEN
    ILOCAL_RECON3(ICONSIDERED)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)&
    +(U_CT(ICONSIDERED)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)))

END IF
END SUBROUTINE EXTRAPOLATE_BOUND_LINEAR




SUBROUTINE EXTRAPOLATE_BOUND_MUSCL(USOL,varcons,FACEX,pointx,ICONSIDERED,SLOPE)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces
IMPLICIT NONE
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED
real,dimension(1:nof_Variables)::leftv
REAL,allocatable,dimension(:),INTENT(IN)::SLOPE
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN)::USOL
REAL::MP_PINFl,gammal


if (WENWRT.EQ.3)THEN
    LEFTV(1:NOF_VARIABLES)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
    CALL cons2prim(N,leftv,MP_PINFl,gammal)
    LEFTV(1:NOF_VARIABLES)=LEFTV(1:NOF_VARIABLES)+USOL(1:nof_Variables,FACEX,pointx)*SLOPE(1:nof_Variables)
    CALL PRIM2CONS(N,LEFTV)
    ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)+LEFTV(1:NOF_VARIABLES)
ELSE
    ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)&
    +(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)*SLOPE(1:nof_Variables)))
END IF


IF (TURBULENCEEQUATIONS.GE.1)THEN
    ILOCAL_RECON3(ICONSIDERED)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)&
    +(U_CT(ICONSIDERED)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)*SLOPE(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)))

END IF


END SUBROUTINE EXTRAPOLATE_BOUND_MUSCL

SUBROUTINE EXTRAPOLATE_BOUND_MUSCLX(USOL,varcons,FACEX,pointx,ICONSIDERED,SLOPE)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces
IMPLICIT NONE
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED
real,dimension(1:nof_Variables)::leftv
REAL,allocatable,dimension(:),INTENT(IN)::SLOPE
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN)::USOL
REAL::MP_PINFl,gammal


if (WENWRT.EQ.3)THEN
    LEFTV(1:NOF_VARIABLES)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
    CALL cons2prim(N,leftv,MP_PINFl,gammal)
    LEFTV(1:NOF_VARIABLES)=LEFTV(1:NOF_VARIABLES)+USOL(1:nof_Variables,FACEX,pointx)*SLOPE(1:nof_Variables)
    CALL PRIM2CONS(N,LEFTV)
    ILOCAL_RECON3(ICONSIDERED)%ULEFTX(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFTX(1:nof_Variables,FACEX,pointx)+LEFTV(1:NOF_VARIABLES)
ELSE
    ILOCAL_RECON3(ICONSIDERED)%ULEFTX(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFTX(1:nof_Variables,FACEX,pointx)&
    +(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)*SLOPE(1:nof_Variables)))
END IF


! IF (TURBULENCEEQUATIONS.GE.1)THEN
!     ILOCAL_RECON3(ICONSIDERED)%ULEFTTURBX(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFTTURBX(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)&
!     +(U_CT(ICONSIDERED)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)*SLOPE(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)))
!
! END IF

END SUBROUTINE EXTRAPOLATE_BOUND_MUSCLX



SUBROUTINE WENOWEIGHTS(N, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, IBOUND_L, U_C_L, U_Ct_L, INTEG_BASIS_L,integ_basis_dg_L, INODER4_L, IEXSOLHIR_L)
!> @brief
!> Subroutine For WENO type reconstruction in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON5_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON6_L
TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C_L, U_Ct_L
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:)::INTEG_BASIS_L,integ_basis_dg_L
TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:)::INODER4_L
TYPE(EXCHANGE_SOLHI),ALLOCATABLE,DIMENSION(:)::IEXSOLHIR_L

REAL::DIVISIONBYZERO
INTEGER::I,J,K,L,M,O,LL,IEX,IEUL,FACX,IELEME,KKD,KMAXE,JF,NGP,IQP,nnd,II,icd
INTEGER::IDUMMY,POWER,ITARGET,ICONSIDERED
REAL::SUMOMEGAATILDEL
REAL::DIVBYZERO,COMPF,checkf,tau_Weno,tempxx
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T


KMAXE=XMPIELRANK(N)

#ifdef WENOWEIGHTS_GPU_KERNEL
!$OMP target teams distribute parallel do
#else
!$OMP DO
#endif
DO II=1,NOF_INTERIOR;
I=EL_INT(II)
ICONSIDERED=I

ielem_L(n,i)%LINC=LWCI1


ILOCAL_RECON3_L(ICONSIDERED)%ULEFT(:,:,:)=ZERO


            if (poly.eq.4)then
            divbyzero=ielem_L(n,iconsidered)%totvolume**2
            else
            divbyzero=10E-12
            end if
            POWER=4

            if (ADDA.EQ.1)THEN
                CALL ADDA_FILTER(N,iconsidered, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, U_C_L, INTEG_BASIS_L, integ_basis_dg_L)
            END IF

            IF (WENWRT.EQ.2)THEN
                CALL CHARACTERISTIC_RECONSTRUCTION(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, IBOUND_L, U_C_L, INTEG_BASIS_L, integ_basis_dg_L, INODER4_L, IEXSOLHIR_L)
            ELSE
                CALL CP_RECONSTRUCTION(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, U_C_L, INTEG_BASIS_L,integ_basis_dg_L)

            END IF

            IF (((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) .and. (icoupleturb.eq.1)) THEN
                CALL CP_RECONSTRUCTION_TURB(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, U_Ct_L, INTEG_BASIS_L, integ_basis_dg_L)

            END IF

END DO

#ifdef WENOWEIGHTS_GPU_KERNEL
!$OMP end target teams distribute parallel do
#else
!$OMP END DO
#endif

!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I


	ILOCAL_RECON3(ICONSIDERED)%ULEFT(:,:,:)=ZERO

            ielem(n,i)%LINC=LWCI1
            if (poly.eq.4)then
            divbyzero=ielem(n,iconsidered)%totvolume**2
            else
            divbyzero=10E-12
            end if
            POWER=4

            if (ADDA.EQ.1)THEN
                CALL ADDA_FILTER(N,iconsidered,  IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, U_C_L, INTEG_BASIS_L, integ_basis_dg_L)
            END IF

            IF (WENWRT.EQ.2)THEN
                CALL CHARACTERISTIC_RECONSTRUCTION(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, IBOUND_L, U_C_L, INTEG_BASIS_L, integ_basis_dg_L, INODER4_L, IEXSOLHIR_L)
            Else
                CALL CP_RECONSTRUCTION(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, U_C_L, INTEG_BASIS_L,integ_basis_dg_L)

            END IF

            IF (((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) .and. (icoupleturb.eq.1)) THEN
                CALL CP_RECONSTRUCTION_TURB(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, U_Ct_L, INTEG_BASIS_L, integ_basis_dg_L)

            END IF
END DO
!$OMP END DO


END SUBROUTINE WENOWEIGHTS

SUBROUTINE CHARACTERISTIC_RECONSTRUCTION(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, IBOUND_L, U_C_L, INTEG_BASIS_L, integ_basis_dg_L, INODER4_L, IEXSOLHIR_L)
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
integer,intent(in)::iconsidered,POWER
integer,intent(inOUT)::IDUMMY
REAL,INTENT(IN)::DIVBYZERO
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L, ILOCAL_RECON5_L
TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C_L
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:)::INTEG_BASIS_L,integ_basis_dg_L
TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:)::INODER4_L
TYPE(EXCHANGE_SOLHI),ALLOCATABLE,DIMENSION(:)::IEXSOLHIR_L


real::angle1,angle2,nx,ny,nz
real,dimension(1:nof_Variables)::veigl,veigr,rveigl,rveigr,leftv,rightv
REAL,DIMENSION(1:nof_Variables,1:nof_Variables)::EIGVL,EIGVR
integer::facex,KKD,l,i,ITARGET,iqp,ngp,iCOMPWRT,icd,IEX,LL,IADMISX,N_FACES,IADMIS,K
real::LWCx1,tau_Weno,ax,ay,az
real::SUMOMEGATILDE(1:nof_Variables)
real,dimension(1:IELEM_L(N,ICONSIDERED)%ADMIS)::LAMC
real,dimension(1:Nof_Variables,0:IDEGFREE)::LIMITEDDW
real,dimension(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:IDEGFREE)::CONSMATRIX
real,dimension(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:IDEGFREE)::CONSMATRIXC
real,dimension(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:Nof_Variables)::RESSOLUTION

real,dimension(1:Nof_Variables,0:IDEGFREE,1:IELEM_L(N,ICONSIDERED)%ADMIS)::LIMITEDDW_CHAR
real,dimension(1:Nof_Variables,1:IELEM_L(N,ICONSIDERED)%ADMIS,0:IDEGFREE)::GRADCHARV

real,dimension(1:nof_Variables,1:IELEM_L(N,ICONSIDERED)%ADMIS,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::LAMBDA
real,dimension(1:nof_Variables,1:IELEM_L(N,ICONSIDERED)%ADMIS,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::SMOOTHINDICATOR
real,dimension(1:nof_Variables,1:IELEM_L(N,ICONSIDERED)%ADMIS,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::omegatilde
real,dimension(1:nof_Variables,1:IELEM_L(N,ICONSIDERED)%ADMIS,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::omega
real,dimension(1:nof_Variables,1:IELEM_L(N,ICONSIDERED)%ADMIS,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::wenoos
real,dimension(1:Nof_Variables,0:IDEGFREE,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::FINDW
real,dimension(1:Nof_Variables,0:IDEGFREE,1:IELEM_L(N,ICONSIDERED)%ADMIS,1:IELEM_L(N,ICONSIDERED)%IFCA,1:2)::FINDW_CHAR

IADMIS=IELEM_L(N,ICONSIDERED)%ADMIS
N_FACES=IELEM_L(N,ICONSIDERED)%IFCA

!ALLOCATE(LAMC(1:IADMIS))
!ALLOCATE(LAMBDA(1:nof_Variables,1:IADMIS,1:N_FACES,1:2))
!ALLOCATE(SMOOTHINDICATOR(1:nof_Variables,1:IADMIS,1:N_FACES,1:2))
!ALLOCATE(OMEGATILDE(1:nof_Variables,1:IADMIS,1:N_FACES,1:2))
!ALLOCATE(OMEGA(1:nof_Variables,1:IADMIS,1:N_FACES,1:2))
!ALLOCATE(WENOOS(1:nof_Variables,1:IADMIS,1:N_FACES,1:2))
!ALLOCATE(LIMITEDDW(1:Nof_Variables,0:IDEGFREE))
!ALLOCATE(LIMITEDDW_CHAR(1:Nof_Variables,0:IDEGFREE,1:IADMIS))
!ALLOCATE(GRADCHARV(1:Nof_Variables,1:IADMIS,0:IDEGFREE))
!ALLOCATE(FINDW(1:Nof_Variables,0:IDEGFREE,1:N_FACES,1:2))
!ALLOCATE(FINDW_CHAR(1:Nof_Variables,0:IDEGFREE,1:IADMIS,1:N_FACES,1:2))
!ALLOCATE(CONSMATRIX(1:NUMBEROFPOINTS2*N_FACES,1:IDEGFREE))
!ALLOCATE(CONSMATRIXC(1:NUMBEROFPOINTS2*N_FACES,1:IDEGFREE))
!ALLOCATE(RESSOLUTION(1:NUMBEROFPOINTS2*N_FACES,1:Nof_Variables))
!allocate(eigvl(1:nof_Variables,1:nof_Variables),EIGVR(1:nof_Variables,1:nof_Variables))



I=ICONSIDERED
lwcx1=IELEM_L(n,i)%LINC

        DO L=1,IELEM_L(N,I)%IFCA  !LOOP FACES
	    !DEFINE
            ANGLE1=IELEM_L(N,I)%FACEANGLEX(L);
            ANGLE2=IELEM_L(N,I)%FACEANGLEY(L)
            FACEX=L

            IF (DIMENSIONA.EQ.3)THEN
            NX=(COS(ANGLE1)*SIN(ANGLE2));
            NY=(SIN(ANGLE1)*SIN(ANGLE2));
            NZ=(COS(ANGLE2))
            ELSE
            NX=ANGLE1
            NY=ANGLE2
            END IF

            VEIGL(1:nof_variables)=U_C_L(I)%VAL(1,1:nof_variables);

            IF (DIMENSIONA.EQ.3)THEN
            CALL ROTATEF(N,RVEIGL,VEIGL,ANGLE1,ANGLE2)
            ELSE
            CALL ROTATEF2D(N,RVEIGL,VEIGL,ANGLE1,ANGLE2)
            END IF
            !NEIGHBOUR!

            CALL WENO_NEIGHBOUR(ICONSIDERED,FACEX,VEIGL,VEIGR,NX,NY,NZ,ANGLE1,ANGLE2,IDUMMY, IELEM_L, ILOCAL_RECON3_L, IBOUND_L, U_C_L, INODER4_L, IEXSOLHIR_L)

            IF (DIMENSIONA.EQ.3)THEN
            CALL ROTATEF(N,RVEIGR,VEIGR,ANGLE1,ANGLE2)
            ELSE
            CALL ROTATEF2D(N,RVEIGR,VEIGR,ANGLE1,ANGLE2)
            END IF

            IF (DIMENSIONA.EQ.3)THEN
            CALL COMPUTE_EIGENVECTORS(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
            ELSE
            CALL COMPUTE_EIGENVECTORS2D(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
            END IF
            LAMBDA(:,:,L,1)=ZERO;
            SMOOTHINDICATOR(:,:,L,1)=ZERO;
            OMEGATILDE(:,:,L,1)=ZERO;
            OMEGA(:,:,L,1)=ZERO;FACEX=L
            CALL compute_gradcharv_smoothindicator(ICONSIDERED,FACEX,EIGVL,GRADCHARV,SMOOTHINDICATOR, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, U_C_L)
            LAMBDA(1:nof_Variables,:,L,1)=1.0D0;LAMBDA(1:nof_Variables,1,L,1)=LWCx1


            if (ees.eq.5)then
            DO KKD=1,nof_variables
                LAMC(1)=(1.0d0-(1.0d0/LWCx1))
                lamc(2:IELEM_L(n,i)%admis)=(1.0d0-lamc(1))/(IELEM_L(N,I)%ADMIS-1)
                LAMBDA(KKD,1:IELEM_L(n,i)%admis,L,1)=lamc(1:IELEM_L(n,i)%admis)
            END DO
            end if


		DO KKD=1,nof_variables
            SUMOMEGATILDE(KKD)=ZERO
                if (ees.eq.5)then
                        tau_Weno=zero
                        if (wenoz.eq.1)then
                                    DO LL=1,IELEM_L(N,I)%ADMIS
                                        tau_Weno=tau_weno+(abs(SMOOTHINDICATOR(KKD,1,L,1)-SMOOTHINDICATOR(KKD,LL,L,1)))
                                    end do
                                        tau_weno=(tau_weno/(IELEM_L(N,I)%ADMIS-1))!**power
                                    DO LL=1,IELEM_L(N,I)%ADMIS
                                        omegatilde(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATOR(KKD,LL,L,1)))**power)
                                    end do
                        else
                                        DO LL=1,IELEM_L(N,I)%ADMIS
                                            OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
                                        END DO
                        end if
                else
                    DO LL=1,IELEM_L(N,I)%ADMIS
                        OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
                    END DO
                end if

			    DO LL=1,IELEM_L(N,I)%ADMIS
				    SUMOMEGATILDE(KKD)=SUMOMEGATILDE(KKD)+OMEGATILDE(KKD,LL,L,1)
			    END DO
			    DO LL=1,IELEM_L(N,I)%ADMIS
				    OMEGA(KKD,LL,L,1)=(OMEGATILDE(KKD,LL,L,1))/SUMOMEGATILDE(KKD)
			    END DO
			    DO LL=1,IELEM_L(N,I)%ADMIS
			    WENOOS(KKD,LL,L,1)=OMEGA(KKD,LL,L,1)
			    END DO


	      END DO       !FINISHED THE LOOP FOR ALL THE VARIABLES


	      LIMITEDDW(:,:)=ZERO
		  IF (EES.EQ.5)THEN
			LIMITEDDW_CHAR(:,:,:)=ZERO
			DO LL=1,IELEM_L(N,I)%ADMIS;IF (LL.EQ.1)THEN
			ITARGET=IELEM_L(N,I)%idegfree
			ELSE
			ITARGET=IDEGFREE2
			END IF
			DO K=0,ITARGET
			LIMITEDDW_CHAR(1:nof_variables,K,1)=LIMITEDDW_CHAR(1:nof_variables,K,1)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW_CHAR(:,:,L,1,:)=ZERO

			DO K=0,IELEM_L(N,I)%idegfree
				FINDW_CHAR(1:nof_variables,K,L,1,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW_CHAR(1:nof_variables,K,1))
			END DO
                 Else
			LIMITEDDW(:,:)=ZERO
			DO K=0,IELEM_L(N,I)%idegfree;DO LL=1,IELEM_L(N,I)%ADMIS
			LIMITEDDW(1:nof_variables,K)=LIMITEDDW(1:nof_variables,K)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW(:,:,L,1)=ZERO
			DO K=0,IELEM_L(N,I)%IDEGFREE
			FINDW(1:nof_variables,K,L,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW(1:nof_variables,K))
			END DO
                  end if



            IF (DIMENSIONA.EQ.3)THEN
            if (IELEM_L(n,i)%types_faces(L).eq.5)then
            iqp=qp_quad
            else
            iqp=qp_triangle
            end if
            Else
            iqp=qp_LINE
            END IF



            icd=0
			do NGP=1,iqp
			    AX = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,1);
			    AY = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,2);
			    IF (DIMENSIONA.EQ.3)THEN
			    AZ = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,3)
			    END IF
			    icd=icd+1;iCOMPWRT=0

			    IF (DIMENSIONA.EQ.3)THEN
			    CONSMATRIX(icd,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
			    ELSE
			    CONSMATRIX(icd,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC2D(N,AX,AY,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
			    END IF

			    iCOMPWRT=0
			    if (ees.eq.5)then;
			    iCOMPWRT=1
                    IF (DIMENSIONA.EQ.3)THEN
                    CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2,icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
                    ELSE
                    CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC2D(N,AX,AY,IORDER2,I,IDEGFREE2,icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
                    END IF

			    iCOMPWRT=0;end if
			END DO


                if (ees.eq.5)then;
                ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,:)=zero
			    DO NGP=1,iqp
			    ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)&
			    +FINDW_char(1:nof_Variables,0,L,1,1)
			    end do

! 			    call DGEMM ('N','T',ICD,nof_variables,IELEM_L(n,i)%idegfree,&
! 							ALPHA,consmatrix(1:icd,1:IELEM_L(n,i)%idegfree),Icd,&
!                             FINDW_char(1:nof_variables,1:IELEM_L(N,I)%IDEGFREE,L,1,1),nof_variables,&
!                             BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),Icd)



                RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(consmatrix(1:icd,1:IELEM_L(n,i)%idegfree),transpose(FINDW_char(1:nof_variables,1:IELEM_L(N,I)%IDEGFREE,L,1,1)))



                            icd=0;
                            do NGP=1,iqp;icd=icd+1
							ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)&
							+RESSOLUTION(icd,1:NOF_vARIABLES)
							END DO

                Else

                 DO NGP=1,iqp
			      ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)=FINDW(1:nof_Variables,0,L,1)
			      end do
! 					      call DGEMM ('N','T',ICD,nof_variables,IELEM_L(n,i)%idegfree,&
! 							ALPHA,consmatrix(1:icd,1:IELEM_L(n,i)%idegfree),Icd,&
!                             FINDW(1:nof_variables,1:IELEM_L(N,I)%IDEGFREE,L,1),nof_variables,&
!                             BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),Icd)


                             RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(consmatrix(1:icd,1:IELEM_L(n,i)%idegfree),transpose(FINDW(1:nof_variables,1:IELEM_L(N,I)%IDEGFREE,L,1)))



					      icd=0;
					      do NGP=1,iqp;
					      icd=icd+1
				      ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3_L(I)%ULEFT(1:nof_Variables,L,NGP)&
				      +RESSOLUTION(icd,1:NOF_vARIABLES)
				      END DO

                end if
            END DO			!FACES




!DEALLOCATE(LAMC,LAMBDA,SMOOTHINDICATOR,OMEGATILDE,&
!OMEGA,WENOOS,LIMITEDDW,LIMITEDDW_CHAR,GRADCHARV,FINDW,&
!FINDW_CHAR,RESSOLUTION,CONSMATRIX,CONSMATRIXC,eigvl,eigvr)





END SUBROUTINE CHARACTERISTIC_RECONSTRUCTION


SUBROUTINE WENO_NEIGHBOUR(ICONSIDERED,FACEX,VEIGL,VEIGR,NX,NY,NZ,ANGLE1,ANGLE2,IDUMMY, IELEM_L, ILOCAL_RECON3_L, IBOUND_L, U_C_L, INODER4_L, IEXSOLHIR_L)
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
real,dimension(1:nof_Variables),INTENT(INOUT)::veigr
real,dimension(1:nof_Variables),INTENT(INOUT)::VEIGL
REAL,INTENT(IN)::NX,NY,NZ,ANGLE1,ANGLE2
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
INTEGER,INTENT(INOUT)::IDUMMY

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:)::U_C_L
TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:)::INODER4_L
TYPE(EXCHANGE_SOLHI),ALLOCATABLE,DIMENSION(:)::IEXSOLHIR_L

REAL::MP_PINFl,gammal
INTEGER::I,J,K,L,var2,B_CODE,N_NODE
real,dimension(1:nof_Variables)::leftv,SRF_SPEED,SRF_SPEEDROT,RIGHTV
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ,CORDS
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(TURBULENCEEQUATIONS)::CTURBL,CTURBR
REAL,DIMENSION(1:NOF_VARIABLES)::CRIGHT_ROT,CLEFT_ROT
INTEGER::IBFC



IF (DIMENSIONA.EQ.3)THEN

L=FACEX
I=ICONSIDERED


 IF (IELEM_L(n,i)%interior.EQ.0)THEN
     VEIGR(1:nof_variables)=U_C_L(IELEM_L(N,I)%INEIGH(L))%VAL(1,1:nof_variables);

 ELSE


                IF (ILOCAL_RECON3_L(I)%MRF.EQ.1)THEN
					SRF_SPEED(2:4)=ILOCAL_RECON3_L(I)%ROTVEL(L,1,1:3)
					CALL ROTATEF(N,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
                END IF
				IF (IELEM_L(N,I)%INEIGHB(l).EQ.N)THEN	!MY CPU ONLY
				      IF (IELEM_L(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if ((IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode.eq.5).or.(IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN MY CPU
					  VEIGR(1:nof_variables)=U_C_L(IELEM_L(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
					  IDUMMY=1
					  IF(PER_ROT.EQ.1)THEN
                        VEIGR(2:4)=ROTATE_PER_1(VEIGR(2:4),IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode,angle_per)
					  END IF
					  else
					  !NOT PERIODIC ONES IN MY CPU

					   CALL coordinates_face_innerx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST, IELEM_L, INODER4_L)


					    if (IELEM_L(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
                                            N_NODE=4
                                    else
                                            N_NODE=3
                                    end if

				  CORDS(1:3)=zero
 				  CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)

				  Poy(1)=cords(2)
				  Pox(1)=cords(1)
				  poz(1)=cords(3)

 				  LEFTV(1:nof_variables)=VEIGL(1:nof_variables)
				  B_CODE=IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode
 				 CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC, IELEM_L, ILOCAL_RECON3_L)

				  VEIGR(1:nof_variables)=RIGHTV(1:nof_variables)
				      	  end if
				      ELSE
				      !FLUID NEIGHBOUR
				      VEIGR(1:nof_variables)=U_C_L(IELEM_L(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
				      END IF
				else
			      !other my cpu
				    IF (IELEM_L(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if ((IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode.eq.5).or.(IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
					  VEIGR(1:nof_variables)=(IEXSOLHIR_L(ILOCAL_RECON3_L(I)%IHEXN(1,IELEM_L(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3_L(I)%IHEXL(1,IELEM_L(N,I)%INDEXI(l)),1:nof_variables))
					  IDUMMY=1
					  IF(PER_ROT.EQ.1)THEN
                        VEIGR(2:4)=ROTATE_PER_1(VEIGR(2:4),IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode,angle_per)
					  END IF
					  end if
				    else

				      VEIGR(1:nof_variables)=(IEXSOLHIR_L(ILOCAL_RECON3_L(I)%IHEXN(1,IELEM_L(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3_L(I)%IHEXL(1,IELEM_L(N,I)%INDEXI(l)),1:nof_variables))

				    end if

				end if



 END IF

Else

L=FACEX
I=ICONSIDERED


 IF (IELEM_L(n,i)%interior.EQ.0)THEN
     VEIGR(1:nof_variables)=U_C_L(IELEM_L(N,I)%INEIGH(L))%VAL(1,1:nof_variables);

 ELSE

 IF (IELEM_L(N,I)%INEIGHB(l).EQ.N)THEN	!MY CPU ONLY
				      IF (IELEM_L(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN MY CPU
					  VEIGR(1:nof_variables)=U_C_L(IELEM_L(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
					  IDUMMY=1
					  else
					  !NOT PERIODIC ONES IN MY CPU

					   CALL coordinates_face_inner2dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST, IELEM_L, INODER4_L)
					   N_NODE=2
				  CORDS(1:2)=zero
 				  CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)


				  Pox(1)=cords(1)
                  Poy(1)=cords(2)

 				  LEFTV(1:nof_variables)=VEIGL(1:nof_variables)
				  B_CODE=IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode
 				  CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC,IELEM_L)

				  VEIGR(1:nof_variables)=RIGHTV(1:nof_variables)
				      	  end if
				      ELSE
				      !FLUID NEIGHBOUR
				      VEIGR(1:nof_variables)=U_C_L(IELEM_L(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
				      END IF
				else
			      !other my cpu
				    IF (IELEM_L(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (IBOUND_L(n,IELEM_L(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					  VEIGR(1:nof_variables)=(IEXSOLHIR_L(ILOCAL_RECON3_L(I)%IHEXN(1,IELEM_L(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3_L(I)%IHEXL(1,IELEM_L(N,I)%INDEXI(l)),1:nof_variables))
					  IDUMMY=1
					  end if
				    else

				      VEIGR(1:nof_variables)=(IEXSOLHIR_L(ILOCAL_RECON3_L(I)%IHEXN(1,IELEM_L(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3_L(I)%IHEXL(1,IELEM_L(N,I)%INDEXI(l)),1:nof_variables))

				    end if

				end if

 END IF


END IF


END SUBROUTINE WENO_NEIGHBOUR


SUBROUTINE CP_RECONSTRUCTION(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, U_C_L, INTEG_BASIS_L,integ_basis_dg_L)
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
integer,intent(in)::iconsidered,POWER
integer,intent(inOUT)::IDUMMY
REAL,INTENT(IN)::DIVBYZERO

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON5_L,ILOCAL_RECON6_L
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:)::INTEG_BASIS_L,integ_basis_dg_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:)::U_C_L

integer::facex,KKD,l,i,ITARGET,iqp,ngp,iCOMPWRT,IEX,LL,IADMIS,N_FACES
real::LWCx1,ax,ay,az,tau_weno,SUMOMEGAATILDEL
INTEGER::ICD
real,dimension(1:nof_Variables)::leftv,rightv
REAL,DIMENSION(1:IDEGFREE)::GRAD1AL
REAL,DIMENSION(1:IDEGFREE)::INDICATEMATRIXAL
REAL,DIMENSION(IDEGFREE)::GRAD3AL

REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::LAMBDAAL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::OMEGAATILDEL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::SMOOTHINDICATORAL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::LAMC,OMEGAAL
REAL,DIMENSION(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:idegfree)::CONSMATRIX
REAL,DIMENSION(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:idegfree)::CONSMATRIXC
REAL,DIMENSION(1:IDEGFREE,1:NOF_VARIABLES)::GRAD5ALc
REAL,DIMENSION(1:IDEGFREE,1:NOF_VARIABLES)::GRADSSL
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:IELEM_L(N,ICONSIDERED)%ADMIS)::WENO
REAL,DIMENSION(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:NOF_vARIABLES)::RESSOLUTION


IADMIS=IELEM_L(N,ICONSIDERED)%ADMIS
N_FACES=IELEM_L(N,ICONSIDERED)%IFCA

!ALLOCATE(GRAD1AL(1:IDEGFREE),INDICATEMATRIXAL(1:IDEGFREE))
!ALLOCATE(GRAD3AL(IDEGFREE),LAMBDAAL(1:IADMIS),OMEGAATILDEL(1:IADMIS),SMOOTHINDICATORAL(1:IADMIS))
!ALLOCATE(LAMC(1:IADMIS),OMEGAAL(1:IADMIS))
!ALLOCATE(CONSMATRIX(1:NUMBEROFPOINTS2*N_FACES,1:idegfree),CONSMATRIXC(1:NUMBEROFPOINTS2*N_FACES,1:idegfree))
!ALLOCATE(GRAD5ALc(1:IDEGFREE,1:NOF_VARIABLES),GRADSSL(1:IDEGFREE,1:NOF_VARIABLES))
!ALLOCATE(WENO(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:IADMIS))
!ALLOCATE(RESSOLUTION(1:NUMBEROFPOINTS2*N_FACES,1:NOF_vARIABLES))





I=ICONSIDERED

lwcx1=IELEM_L(n,i)%LINC


DO IEX=1,nof_variables
			LAMBDAAL=ZERO;SMOOTHINDICATORAL=ZERO;OMEGAATILDEL=ZERO;OMEGAAL=ZERO
                IF (EES.EQ.5)THEN
                            LAMC(:)=ZERO; GRAD3AL(:)=ZERO; LAMC(1)=(1.0d0-(1.0d0/lwcx1));lamc(2:IELEM_L(n,i)%admis)=(1.0d0-lamc(1))/(IELEM_L(N,I)%ADMIS-1)
                            LAMBDAAL(1:IELEM_L(n,i)%admis)=lamc(1:IELEM_L(n,i)%admis)
                            !sum the low degree polynomials first
                            DO LL=2,IELEM_L(N,I)%ADMIS
                            GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTSC(LL,1:IDEGFREE2,IEX))
                            END DO
                            !this is the zero polynomial
                            GRAD1AL(1:IELEM_L(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM_L(N,I)%IDEGFREE))
                            GRAD5ALc(1:IELEM_L(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM_L(N,I)%IDEGFREE)
                        DO LL=1,IELEM_L(N,I)%ADMIS
                            IF (LL.EQ.1)THEN

!                                 CALL DGEMV('N', IELEM_L(N,I)%IDEGFREE, IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                                 ILOCAL_RECON3(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),&
!                                 IELEM_L(N,I)%IDEGFREE,GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),1,BETA,INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE),1)


                                INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE)=matmul(ILOCAL_RECON3_L(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),GRAD1AL(1:IELEM_L(N,I)%IDEGFREE))

                                SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE))

                            ELSE
                                GRAD1AL(1:IDEGFREE2)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTSC(ll,1:IDEGFREE2,IEX)
            !

!                                 CALL DGEMV('N', IDEGFREE2, IDEGFREE2,ALPHA,&
!                                 ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),&
!                                 IDEGFREE2,GRAD1AL(1:IDEGFREE2),1,BETA,INDICATEMATRIXAL(1:IDEGFREE2),1)

                                INDICATEMATRIXAL(1:IDEGFREE2)=matmul(ILOCAL_RECON3_L(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),GRAD1AL(1:IDEGFREE2))

                                SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                            END IF
                        END DO
                ELSE
                        DO LL=1,IELEM_L(N,I)%ADMIS
                        GRAD1AL(:)=ZERO
                        INDICATEMATRIXAL(:)=ZERO
                        GRAD1AL(1:IELEM_L(N,I)%IDEGFREE)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTS(LL,1:IELEM_L(N,I)%IDEGFREE,IEX)

!                         CALL DGEMV('N', IELEM_L(N,I)%IDEGFREE, IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                             ILOCAL_RECON3(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),&
!                             IELEM_L(N,I)%IDEGFREE,GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),1,BETA,INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE),1)

                            INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE)=matmul(ILOCAL_RECON3_L(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),GRAD1AL(1:IELEM_L(N,I)%IDEGFREE))


                        SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE))
                        END DO
                END IF
			      LAMBDAAL(:)=1.0D0
			      LAMBDAAL(1)=lwcx1

                if (ees.eq.5)then
				LAMC(1)=(1.0d0-(1.0d0/lwcx1))
				lamc(2:IELEM_L(n,i)%admis)=(1.0d0-lamc(1))/(IELEM_L(N,I)%ADMIS-1)
				LAMBDAAL(1:IELEM_L(n,i)%admis)=lamc(1:IELEM_L(n,i)%admis)
				end if




			     if (ees.eq.5)then

				      if (wenoz.eq.1)then
					      tau_Weno=zero
					      DO LL=1,IELEM_L(N,I)%ADMIS
					      tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
					      end do
					      tau_weno=(tau_weno/(IELEM_L(N,I)%ADMIS-1))
					      DO LL=1,IELEM_L(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL)))**power)
					      end do
				      else
					      DO LL=1,IELEM_L(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					      END DO

				      end if
                             else
					  DO LL=1,IELEM_L(N,I)%ADMIS
					  OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					  END DO
			      end if


			      SUMOMEGAATILDEL=ZERO
			      DO LL=1,IELEM_L(N,I)%ADMIS
			      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
			      END DO
			      DO LL=1,IELEM_L(N,I)%ADMIS
			      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
			      END DO

			      DO LL=1,IELEM_L(N,I)%ADMIS
			      WENO(IEX,LL)=OMEGAAL(LL)

			      if (iex.eq.1)then
				    IELEM_L(n,i)%wcx(1)=WENO(iex,1)
				    end if

			      END DO
		  END DO


		  icd=0
		DO L=1,IELEM_L(N,I)%IFCA	!FACES
                  IF (DIMENSIONA.EQ.3)THEN
                        if (IELEM_L(n,i)%types_faces(L).eq.5)then
                            iqp=qp_quad
                        else
                            iqp=qp_triangle
                        end if
				  ELSE
                            iqp=qp_LINE
				  END IF

			  do NGP=1,iqp			!for gqp
                 icd=icd+1

				AX = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,2)
				IF (DIMENSIONA.EQ.3)THEN
				AZ = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,3)
				END IF
                Icompwrt=0

                        IF (DIMENSIONA.EQ.3)THEN
	      				CONSMATRIX(icd,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L,integ_basis_dg_L)
	      				ELSE
	      				CONSMATRIX(icd,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC2D(N,AX,AY,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L,integ_basis_dg_L)
	      				END IF

                        if (ees.eq.5)then
                                Icompwrt=1

                                IF (DIMENSIONA.EQ.3)THEN
                                CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L,integ_basis_dg_L)
                                ELSE
                                CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC2D(N,AX,AY,IORDER2,I,IDEGFREE2,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L,integ_basis_dg_L)
                                END IF
                                Icompwrt=0
                        END IF
				end do

		    END DO	!FACES

                                ILOCAL_RECON3_L(I)%ULEFT(:,:,:)=ZERO


                IF (DG.EQ.1)THEN
				ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,:)=ZERO
				END IF



				DO LL=1,IELEM_L(N,I)%ADMIS	!STENCILS

				IF (EES.EQ.5)THEN
                    IF (LL.EQ.1)THEN
                        GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:nof_variables)=GRAD5ALc(1:IELEM_L(N,I)%IDEGFREE,1:nof_variables)


!                     CALL DGEMM('N','N',ICD,nof_variables,IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),ICD,&
!                     GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES),&
!                     IELEM_L(N,I)%IDEGFREE,BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),ICD)


                    RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))



                    else
				GRADSSL(1:Idegfree2,1:nof_variables)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTSc(LL,1:idegfree2,1:nof_variables)




!                                  CALL DGEMM('N','N',ICD,nof_variables,idegfree2,ALPHA,&
!                     CONSMATRIXc(1:ICD,1:IDEGFREE2),ICD,&
!                     GRADSSL(1:IDEGFREE2,1:NOF_vARIABLES),&
!                     IDEGFREE2,BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),ICD)


                RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(CONSMATRIXc(1:ICD,1:IDEGFREE2),GRADSSL(1:IDEGFREE2,1:NOF_vARIABLES))


                    end if




				ELSE
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:nof_variables)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTS(LL,1:IELEM_L(N,I)%IDEGFREE,1:nof_variables)



!                                 CALL DGEMM('N','N',ICD,nof_variables,IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),ICD,&
!                     GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES),&
!                     IELEM_L(N,I)%IDEGFREE,BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),ICD)


                    RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))


				END IF



                                 ICD=0
                                DO L=1,IELEM_L(N,I)%IFCA
                                                    IF (DIMENSIONA.EQ.3)THEN
                                                        if (IELEM_L(n,i)%types_faces(L).eq.5)then
                                                        iqp=qp_quad;
                                                        else
                                                        iqp=qp_triangle;
                                                        end if
                                                    ELSE
                                                        iqp=qp_LINE;
                                                    END IF

                                    do NGP=1,iqp
                                        ICD=ICD+1
                                            CALL  EXTRAPOLATE_BOUND(RESSOLUTION,IEX,L,NGP,I,ICD,LL,WENO, ILOCAL_RECON3_L, U_C_L)
                                    end do
                                end do

                                        DO IEX=1,nof_variables	!COMPONENTS
                                            IF (DG.EQ.1)THEN
                                                    IF (EES.EQ.5)THEN
                                                        IF (LL.EQ.1)THEN
                                                        ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)=ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)+(GRADSSL(1:Idegfree,IEX)*WENO(IEX,LL))
                                                        ELSE
                                                        ILOCAL_RECON6_L(I)%DG2FV(1:IDEGFREE2,IEX)=ILOCAL_RECON6_L(I)%DG2FV(1:IDEGFREE2,IEX)+(GRADSSL(1:Idegfree2,IEX)*WENO(IEX,LL))
                                                        END IF
                                                    ELSE
                                                        ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)=ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)+(GRADSSL(1:Idegfree,IEX)*WENO(IEX,LL))
                                                    END IF
                                            END IF
                                        end do
                                end do  !STENCILS FINISHED



                                        if (wenwrt.eq.3)then

                                        DO L=1,IELEM_L(N,I)%IFCA
                                                            IF (DIMENSIONA.EQ.3)THEN
                                                                if (IELEM_L(n,i)%types_faces(L).eq.5)then
                                                                iqp=qp_quad;
                                                                else
                                                                iqp=qp_triangle;
                                                                end if
                                                            ELSE
                                                                iqp=qp_LINE
                                                            END IF


                                                                    do NGP=1,iqp
                                                                    leftv(1:nof_variables)=ILOCAL_RECON3_L(I)%ULEFT(1:NOF_vARIABLES,l,ngp)
                                                                    call PRIM2CONS(N,leftv)
                                                                    ILOCAL_RECON3_L(I)%ULEFT(1:NOF_vARIABLES,l,ngp)=leftv(1:nof_variables)
                                                                    end do
                                        end do

                                        end if








!DEALLOCATE(GRAD1AL,INDICATEMATRIXAL,GRAD3AL,LAMBDAAL,OMEGAATILDEL,SMOOTHINDICATORAL)
!DEALLOCATE(LAMC,OMEGAAL)
!DEALLOCATE(CONSMATRIX,CONSMATRIXC,GRAD5ALc,GRADSSL)
!DEALLOCATE(WENO)
!DEALLOCATE(RESSOLUTION)












END SUBROUTINE CP_RECONSTRUCTION



SUBROUTINE CP_RECONSTRUCTION_Turb(ICONSIDERED,IDUMMY,DIVBYZERO,POWER, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, ILOCAL_RECON6_L, U_Ct_L, INTEG_BASIS_L, integ_basis_dg_L)
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
integer,intent(in)::iconsidered,POWER
integer,intent(inOUT)::IDUMMY
REAL,INTENT(IN)::DIVBYZERO

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON5_L,ILOCAL_RECON6_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:)::U_Ct_L
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:)::INTEG_BASIS_L,integ_basis_dg_L

integer::facex,KKD,l,i,ITARGET,iqp,ngp,iCOMPWRT,icd,IEX,LL,iadmis,n_faces
real::LWCx1,tau_Weno,ax,ay,az,SUMOMEGAATILDEL
real,dimension(1:nof_Variables)::leftv,rightv
REAL,DIMENSION(1:IDEGFREE)::GRAD1AL
REAL,DIMENSION(1:IDEGFREE)::INDICATEMATRIXAL
REAL,DIMENSION(IDEGFREE)::GRAD3AL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::LAMBDAAL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::OMEGAATILDEL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::SMOOTHINDICATORAL
REAL,DIMENSION(1:IELEM_L(N,ICONSIDERED)%ADMIS)::LAMC,OMEGAAL
REAL,DIMENSION(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:idegfree)::CONSMATRIX,CONSMATRIXC
REAL,DIMENSION(1:IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR) :: GRAD5ALc,GRADSSL
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:IELEM_L(N,ICONSIDERED)%ADMIS) ::WENO
REAL,DIMENSION(1:NUMBEROFPOINTS2*IELEM_L(N,ICONSIDERED)%IFCA,1:TURBULENCEEQUATIONS+PASSIVESCALAR) :: RESSOLUTION


IADMIS=IELEM_L(N,ICONSIDERED)%ADMIS
N_FACES=IELEM_L(N,ICONSIDERED)%IFCA

!ALLOCATE(GRAD1AL(1:IDEGFREE),INDICATEMATRIXAL(1:IDEGFREE))
!ALLOCATE(GRAD3AL(IDEGFREE),LAMBDAAL(1:IADMIS),OMEGAATILDEL(1:IADMIS),SMOOTHINDICATORAL(1:IADMIS))
!ALLOCATE(LAMC(1:IADMIS),OMEGAAL(1:IADMIS))
!ALLOCATE(CONSMATRIX(1:NUMBEROFPOINTS2*N_FACES,1:idegfree),CONSMATRIXC(1:NUMBEROFPOINTS2*N_FACES,1:idegfree))
!ALLOCATE(GRAD5ALc(1:IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR),GRADSSL(1:IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))
!ALLOCATE(WENO(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:IADMIS))
!LLOCATE(RESSOLUTION(1:NUMBEROFPOINTS2*N_FACES,1:TURBULENCEEQUATIONS+PASSIVESCALAR))










I=ICONSIDERED

lwcx1=IELEM_L(n,i)%LINC



DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			LAMBDAAL=ZERO;SMOOTHINDICATORAL=ZERO;OMEGAATILDEL=ZERO;OMEGAAL=ZERO
                IF (EES.EQ.5)THEN
                            LAMC(:)=ZERO; GRAD3AL(:)=ZERO; LAMC(1)=(1.0d0-(1.0d0/lwcx1));lamc(2:IELEM_L(n,i)%admis)=(1.0d0-lamc(1))/(IELEM_L(N,I)%ADMIS-1)
                            LAMBDAAL(1:IELEM_L(n,i)%admis)=lamc(1:IELEM_L(n,i)%admis)
                            !sum the low degree polynomials first
                            DO LL=2,IELEM_L(N,I)%ADMIS
                            GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTSC2(LL,1:IDEGFREE2,IEX))
                            END DO
                            !this is the zero polynomial
                            GRAD1AL(1:IELEM_L(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTS2(1,1:IELEM_L(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM_L(N,I)%IDEGFREE))
                            GRAD5ALc(1:IELEM_L(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM_L(N,I)%IDEGFREE)
                        DO LL=1,IELEM_L(N,I)%ADMIS
                            IF (LL.EQ.1)THEN

!                                 CALL DGEMV('N', IELEM_L(N,I)%IDEGFREE, IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                                 ILOCAL_RECON3(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),&
!                                 IELEM_L(N,I)%IDEGFREE,GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),1,BETA,INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE),1)

                                INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE)=matmul(ILOCAL_RECON3_L(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),GRAD1AL(1:IELEM_L(N,I)%IDEGFREE))


                                SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE))

                            ELSE
                                GRAD1AL(1:IDEGFREE2)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTSC2(ll,1:IDEGFREE2,IEX)
            !

!                                 CALL DGEMV('N', IDEGFREE2, IDEGFREE2,ALPHA,&
!                                 ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),&
!                                 IDEGFREE2,GRAD1AL(1:IDEGFREE2),1,BETA,INDICATEMATRIXAL(1:IDEGFREE2),1)

                                SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                            END IF
                        END DO
                ELSE
                        DO LL=1,IELEM_L(N,I)%ADMIS
                        GRAD1AL(:)=ZERO
                        INDICATEMATRIXAL(:)=ZERO
                        GRAD1AL(1:IELEM_L(N,I)%IDEGFREE)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTS2(LL,1:IELEM_L(N,I)%IDEGFREE,IEX)

!                         CALL DGEMV('N', IELEM_L(N,I)%IDEGFREE, IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                             ILOCAL_RECON3(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),&
!                             IELEM_L(N,I)%IDEGFREE,GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),1,BETA,INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE),1)


                           INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE)=matmul(ILOCAL_RECON3_L(I)%INDICATOR(1:IELEM_L(N,I)%IDEGFREE,1:IELEM_L(N,I)%IDEGFREE),GRAD1AL(1:IELEM_L(N,I)%IDEGFREE))


                        SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM_L(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM_L(N,I)%IDEGFREE))
                        END DO
                END IF
			      LAMBDAAL(:)=1.0D0
			      LAMBDAAL(1)=lwcx1

                if (ees.eq.5)then
				LAMC(1)=(1.0d0-(1.0d0/lwcx1))
				lamc(2:IELEM_L(n,i)%admis)=(1.0d0-lamc(1))/(IELEM_L(N,I)%ADMIS-1)
				LAMBDAAL(1:IELEM_L(n,i)%admis)=lamc(1:IELEM_L(n,i)%admis)
				end if




			     if (ees.eq.5)then

				      if (wenoz.eq.1)then
					      tau_Weno=zero
					      DO LL=1,IELEM_L(N,I)%ADMIS
					      tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
					      end do
					      tau_weno=(tau_weno/(IELEM_L(N,I)%ADMIS-1))
					      DO LL=1,IELEM_L(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL)))**power)
					      end do
				      else
					      DO LL=1,IELEM_L(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					      END DO

				      end if
                             else
					  DO LL=1,IELEM_L(N,I)%ADMIS
					  OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					  END DO
			      end if


			      SUMOMEGAATILDEL=ZERO
			      DO LL=1,IELEM_L(N,I)%ADMIS
			      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
			      END DO
			      DO LL=1,IELEM_L(N,I)%ADMIS
			      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
			      END DO

			      DO LL=1,IELEM_L(N,I)%ADMIS
			      WENO(IEX,LL)=OMEGAAL(LL)



			      END DO
		  END DO


		  icd=0
		DO L=1,IELEM_L(N,I)%IFCA	!FACES
                  IF (DIMENSIONA.EQ.3)THEN
				  if (IELEM_L(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				  else
				    iqp=qp_triangle
				  end if
				  ELSE
                    iqp=qp_LINE
				  END IF

			  do NGP=1,iqp			!for gqp
                 icd=icd+1

				AX = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,2)
				IF (DIMENSIONA.EQ.3)THEN
				AZ = ILOCAL_RECON3_L(I)%QPOINTS(L,NGP,3)
				END IF
                Icompwrt=0

                        IF (DIMENSIONA.EQ.3)THEN
	      				CONSMATRIX(icd,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L,integ_basis_dg_L)
	      				ELSE
	      				CONSMATRIX(icd,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC2D(N,AX,AY,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L,integ_basis_dg_L)
	      				END IF

                 if (ees.eq.5)then
	      				Icompwrt=1

	      				IF (DIMENSIONA.EQ.3)THEN
	      				CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
	      				ELSE
	      				CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC2D(N,AX,AY,IORDER2,I,IDEGFREE2,Icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
	      				END IF
	      				Icompwrt=0
                 END IF
				end do

		    END DO	!FACES

                                ILOCAL_RECON3_L(I)%ULEFTturb(:,:,:)=ZERO


                IF (DG.EQ.1)THEN
				ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,:)=ZERO
				END IF



				DO LL=1,IELEM_L(N,I)%ADMIS	!STENCILS

				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=GRAD5ALc(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)


!                     CALL DGEMM('N','N',ICD,TURBULENCEEQUATIONS+PASSIVESCALAR,IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),ICD,&
!                     GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR),&
!                     IELEM_L(N,I)%IDEGFREE,BETA,RESSOLUTION(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR),ICD)

                    RESSOLUTION(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=matmul(CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))


				else
				GRADSSL(1:Idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTSc2(LL,1:idegfree2,1:TURBULENCEEQUATIONS+PASSIVESCALAR)




!                                  CALL DGEMM('N','N',ICD,TURBULENCEEQUATIONS+PASSIVESCALAR,idegfree2,ALPHA,&
!                     CONSMATRIXc(1:ICD,1:IDEGFREE2),ICD,&
!                     GRADSSL(1:IDEGFREE2,1:TURBULENCEEQUATIONS+PASSIVESCALAR),&
!                     IDEGFREE2,BETA,RESSOLUTION(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR),ICD)


                    RESSOLUTION(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=matmul(CONSMATRIXc(1:ICD,1:IDEGFREE2),GRADSSL(1:IDEGFREE2,1:TURBULENCEEQUATIONS+PASSIVESCALAR))



				end if




				ELSE
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_rECON5_L(ICONSIDERED)%GRADIENTS2(LL,1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)



!                                 CALL DGEMM('N','N',ICD,TURBULENCEEQUATIONS+PASSIVESCALAR,IELEM_L(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),ICD,&
!                     GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR),&
!                     IELEM_L(N,I)%IDEGFREE,BETA,RESSOLUTION(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR),ICD)

                    RESSOLUTION(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=matmul(CONSMATRIX(1:ICD,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))



				END IF



                                 ICD=0
                                DO L=1,IELEM_L(N,I)%IFCA
                                    IF (DIMENSIONA.EQ.3)THEN
                                                    if (IELEM_L(n,i)%types_faces(L).eq.5)then
                                                    iqp=qp_quad;
                                                    else
                                                    iqp=qp_triangle;
                                                    end if
                                                    ELSE
                                                    iqp=qp_LINE;
                                                    END IF

                                    do NGP=1,iqp
                                        ICD=ICD+1
                                            CALL  EXTRAPOLATE_BOUNDt(RESSOLUTION,IEX,L,NGP,I,ICD,LL,WENO, ILOCAL_RECON3_L, U_Ct_L)
                                    end do
                                end do

                                DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR	!COMPONENTS
                                IF (DG.EQ.1)THEN
                                    IF (EES.EQ.5)THEN
                                        IF (LL.EQ.1)THEN
                                        ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)=ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)+(GRADSSL(1:Idegfree,IEX)*WENO(IEX,LL))
                                        ELSE
                                        ILOCAL_RECON6_L(I)%DG2FV(1:IDEGFREE2,IEX)=ILOCAL_RECON6_L(I)%DG2FV(1:IDEGFREE2,IEX)+(GRADSSL(1:Idegfree2,IEX)*WENO(IEX,LL))
                                        END IF
                                    ELSE
                                        ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)=ILOCAL_RECON6_L(I)%DG2FV(1:IELEM_L(N,I)%IDEGFREE,IEX)+(GRADSSL(1:Idegfree,IEX)*WENO(IEX,LL))
                                    END IF
                                    END IF
                                end do
                                end do








! 	END IF

!    DEALLOCATE(GRAD1AL,INDICATEMATRIXAL,GRAD3AL,LAMBDAAL,OMEGAATILDEL,SMOOTHINDICATORAL)
!    DEALLOCATE(LAMC,OMEGAAL)
!   DEALLOCATE(CONSMATRIX,CONSMATRIXC,GRAD5ALc,GRADSSL)
!    DEALLOCATE(WENO)
 !   DEALLOCATE(RESSOLUTION)


END SUBROUTINE CP_RECONSTRUCTION_Turb



SUBROUTINE EXTRAPOLATE_BOUND(RESSOLUTION,varcons,FACEX,pointx,ICONSIDERED,INSTEN,LLX,WENO, ILOCAL_RECON3_L, U_C_L)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED,INSTEN,LLX
REAL,dimension(:,:),INTENT(IN)::WENO
REAL,DIMENSION(:,:),intent(in)::RESSOLUTION

TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:)::U_C_L

real,dimension(1:nof_Variables)::leftv
REAL::MP_PINFl,gammal


					if (WENWRT.EQ.3)THEN	!PRIMITIVE
					LEFTV(1:NOF_VARIABLES)=U_C_L(ICONSIDERED)%VAL(1,1:nof_Variables)
                    call CONS2PRIM(N,leftv,MP_PINFl,gammal)

                    ILOCAL_RECON3_L(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)=ILOCAL_RECON3_L(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)&
				    +((leftv(1:NOF_VARIABLES)+RESSOLUTION(INSTEN,1:NOF_vARIABLES))*WENO(1:NOF_vARIABLES,llx))
                    else
											!CONSERVATIVE

                        ILOCAL_RECON3_L(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)=ILOCAL_RECON3_L(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)&
				     +(U_C_L(ICONSIDERED)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(INSTEN,1:NOF_vARIABLES))*WENO(1:NOF_vARIABLES,LLX)

				    end if



END SUBROUTINE EXTRAPOLATE_BOUND



SUBROUTINE EXTRAPOLATE_BOUNDt(RESSOLUTION,varcons,FACEX,pointx,ICONSIDERED,INSTEN,LLX,WENO, ILOCAL_RECON3_L, U_Ct_L)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED,INSTEN,LLX
REAL,dimension(:,:),INTENT(IN)::WENO
real,dimension(1:nof_Variables)::leftv
REAL,DIMENSION(:,:),intent(in)::RESSOLUTION

TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:)::U_Ct_L

REAL::MP_PINFl,gammal


				     ILOCAL_RECON3_L(ICONSIDERED)%ULEFTturb(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)=ILOCAL_RECON3_L(ICONSIDERED)%ULEFTturb(1:TURBULENCEEQUATIONS+PASSIVESCALAR,FACEX,pointx)&
				     +(U_Ct_L(ICONSIDERED)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+RESSOLUTION(INSTEN,1:TURBULENCEEQUATIONS+PASSIVESCALAR))*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,LLX)





END SUBROUTINE EXTRAPOLATE_BOUNDt


subroutine diag_At_B_A(ICONSIDERED,A_CHAR,B_CHAR,X_CHAR, IELEM_L, ILOCAL_RECON3_L)
!> @brief
!> Subroutine For general matrix A and square matrix B, computes the vector x = diag(A' * B * A), used for characteristics reconstruction
implicit none
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
   integer, intent(in):: ICONSIDERED
   REAL,dimension(:,:),INTENT(INout)::b_char,x_char
   REAL,dimension(:,:,:),INTENT(INout)::a_char
   

   TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
   TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
   real, dimension(size(A_char,1), size(A_char,2), IELEM_L(n,iconsidered)%admis) :: BA_char

   integer:: nn, mm
   integer:: i,LL,ICS

   nn = size(A_char,1) ! = size(B,1) = size(B,2)
   mm = size(A_char,2)

   !allocate(BA_char(nn, mm,IELEM_L(n,iconsidered)%admis))

   x_char=zero

   IF (EES.EQ.5)THEN
	  DO LL=1,1

! 	  CALL DGEMM('N','N',IELEM_L(N,ICONSIDERED)%idegfree,nof_variables,IELEM_L(N,ICONSIDERED)%idegfree,ALPHA,&
! 	  B_char(1:IELEM_L(N,ICONSIDERED)%idegfree,1:IELEM_L(N,ICONSIDERED)%idegfree),IELEM_L(N,ICONSIDERED)%idegfree,&
! 	  A_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_variables,LL),&
!         IELEM_L(N,ICONSIDERED)%idegfree,BETA,BA_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_Variables,LL),&
!         IELEM_L(N,ICONSIDERED)%idegfree)


        BA_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_Variables,LL)=matmul(B_char(1:IELEM_L(N,ICONSIDERED)%idegfree,1:IELEM_L(N,ICONSIDERED)%idegfree),A_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_variables,LL))


	  END DO
	B_CHAR(1:IDEGFREE2,1:IDEGFREE2)=ILOCAL_RECON3_L(ICONSIDERED)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2)
	  DO LL=2,IELEM_L(N,ICONSIDERED)%ADMIS

! 	  call gemm(B_char(1:IDEGFREE2,1:IDEGFREE2), A_char(1:IDEGFREE2,1:nof_variables,LL), BA_CHAR(1:IDEGFREE2,1:nof_Variables,LL)) ! BA = B * A

! 	  CALL DGEMM('N','N',IDEGFREE2,nof_variables,IDEGFREE2,ALPHA,B_char(1:IDEGFREE2,1:IDEGFREE2),IDEGFREE2,&
! 	  A_CHAR(1:IDEGFREE2,1:nof_variables,LL),&
!         IDEGFREE2,BETA,BA_CHAR(1:IDEGFREE2,1:nof_Variables,LL),IDEGFREE2)


        BA_CHAR(1:IDEGFREE2,1:nof_Variables,LL)=matmul(B_char(1:IDEGFREE2,1:IDEGFREE2),A_CHAR(1:IDEGFREE2,1:nof_variables,LL))


	  END DO

	  DO LL=1,1;do i = 1, mm
	    !x_char(i,LL) = dot(A_char(:,i,LL), BA_char(:,i,LL))
	    x_char(i,LL) =DOT_PRODUCT(a_char(1:IELEM_L(N,ICONSIDERED)%idegfree,i,ll),BA_char(1:IELEM_L(N,ICONSIDERED)%idegfree,i,LL))
	  end do;END DO
	  DO LL=2,IELEM_L(N,ICONSIDERED)%ADMIS;do i = 1, mm
	    !x_char(i,LL) = dot(A_char(1:IDEGFREE2,i,LL), BA_char(1:IDEGFREE2,i,LL))
	    x_char(i,LL) =DOT_PRODUCT(a_char(1:IDEGFREE2,i,ll),BA_char(1:IDEGFREE2,i,LL))
	  end do;END DO


   ELSE

	      DO LL=1,IELEM_L(N,ICONSIDERED)%ADMIS
! 	      call gemm(B_char(:,:), A_char(:,:,LL), BA_char(:,:,LL)) ! BA = B * A

! 	      CALL DGEMM('N','N',IELEM_L(N,ICONSIDERED)%idegfree,nof_variables,IELEM_L(N,ICONSIDERED)%idegfree,ALPHA,&
! 	      B_char(1:IELEM_L(N,ICONSIDERED)%idegfree,1:IELEM_L(N,ICONSIDERED)%idegfree),IELEM_L(N,ICONSIDERED)%idegfree,&
! 	      A_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_variables,LL),&
!             IELEM_L(N,ICONSIDERED)%idegfree,BETA,BA_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_Variables,LL),&
!             IELEM_L(N,ICONSIDERED)%idegfree)

            BA_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_Variables,LL)=matmul(B_char(1:IELEM_L(N,ICONSIDERED)%idegfree,1:IELEM_L(N,ICONSIDERED)%idegfree),A_CHAR(1:IELEM_L(N,ICONSIDERED)%idegfree,1:nof_variables,LL))



	      END DO
	      DO LL=1,IELEM_L(N,ICONSIDERED)%ADMIS;do i = 1, mm
		  !x_char(i,LL) = dot(A_char(:,i,LL), BA_char(:,i,LL))
		   x_char(i,LL) =DOT_PRODUCT(a_char(1:IELEM_L(N,ICONSIDERED)%idegfree,i,ll),BA_char(1:IELEM_L(N,ICONSIDERED)%idegfree,i,LL))
	      end do;END DO
   END IF

   !deallocate(BA_CHAR)
end subroutine diag_At_B_A

subroutine compute_gradcharv_smoothindicator(ICONSIDERED,FACEX,EIGVL,GRADCHARV,SMOOTHINDICATOR, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, U_C_L)
!> @brief
!> Subroutine for characteristics reconstruction of WENO schemes
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
   integer, intent(in)::ICONSIDERED,FACEX
   integer:: LL, k,I,L,ifds
   REAL::LWCX1
   REAL,DIMENSION(1:nof_Variables,1:nof_Variables),INTENT(IN)::EIGVL
   REAL,dimension(:,:,:),INTENT(INOUT)::GRADCHARV
   REAL,dimension(:,:,:,:),INTENT(INOUT)::SMOOTHINDICATOR

   TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
   TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
   TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON5_L
   TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C_L

   real,dimension(1:TYPESTEN)::lamc
   real,dimension(1:IDEGFREE,1:nof_variables)::grad5alc
   real,dimension(1:idegfree,1:idegfree)::b_char
   real,dimension(1:nof_variables,1:TYPESTEN)::x_char
   real,dimension(1:idegfree,1:nof_variables,1:TYPESTEN)::a_char
   real,dimension(0:idegfree,1:nof_variables,1:typesten)::gradients
   real,dimension(0:idegfree,1:nof_variables,1:typesten)::gradients_eigvlt

   !allocate(LAMC(1:TYPESTEN))
   !allocate(GRAD5ALC(1:IDEGFREE,1:nof_variables))
   !allocate(A_CHAR(1:idegfree,1:nof_variables,1:TYPESTEN),B_CHAR(1:idegfree,1:idegfree))
   !allocate(X_CHAR(1:nof_variables,1:TYPESTEN))
   !allocate(GRADIENTS(0:idegfree,1:nof_variables,1:typesten))
   !allocate(GRADIENTS_EIGVLT(0:idegfree,1:nof_variables,1:typesten))





   I=ICONSIDERED
   L=FACEX
   GRADCHARV=zero
   lwcx1=IELEM_L(n,i)%LINC

   gradients(:,:,:)=ZERO;gradients_eigvlt(:,:,:)=zero
   DO LL=1,IELEM_L(N,I)%ADMIS
    gradients(0,:,ll) = U_C_L(I)%VAL(1,1:nof_variables)
   end do
   IF (EES.EQ.5)THEN
   ;LAMC(:)=ZERO;GRAD5ALC=ZERO;LAMC(1)=(1.0d0-(1.0d0/LWCx1));lamc(2:IELEM_L(n,i)%admis)=(1.0d0-lamc(1))/(IELEM_L(N,I)%ADMIS-1)
		    DO LL=2,IELEM_L(N,I)%ADMIS
			GRAD5ALC(1:IDEGFREE2,1:nof_variables)=GRAD5ALC(1:IDEGFREE2,1:nof_variables)&
			+(LAMC(LL)*ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTSC(ll,1:IDEGFREE2,1:nof_variables))
			gradients(1:IDEGFREE2,1:nof_variables,ll)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTSC(LL,1:IDEGFREE2,1:nof_variables)
		    END DO
		    DO LL=1,1
		      gradients(1:IELEM_L(N,I)%idegfree,1:nof_variables,ll)=(1.0D0/LAMC(1))*&
		      (ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,1:nof_variables)-GRAD5ALC(1:IELEM_L(N,I)%IDEGFREE,1:nof_variables))
		    END DO
		    DO LL=1,1


! 			call DGEMM ('N','T',IELEM_L(N,I)%idegfree+1,nof_variables,nof_variables,&
! 			ALPHA,gradients(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll),IELEM_L(N,I)%idegfree+1,&
!             EIGVL(1:nof_variables,1:nof_variables),nof_variables,BETA,&
!             gradients_eigvlt(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll),IELEM_L(N,I)%idegfree+1)


            gradients_eigvlt(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll)=matmul(gradients(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll),transpose(EIGVL(1:nof_variables,1:nof_variables)))




		      END DO
		    DO LL=2,IELEM_L(N,I)%ADMIS

! 			call DGEMM ('N','T',IDEGFREE2+1,nof_variables,nof_variables,&
! 			ALPHA,gradients(0:IDEGFREE2,1:nof_variables,ll),IDEGFREE2+1,&
!             EIGVL(1:nof_variables,1:nof_variables),nof_variables,BETA,&
!             gradients_eigvlt(0:IDEGFREE2,1:nof_variables,ll),IDEGFREE2+1)

            gradients_eigvlt(0:IDEGFREE2,1:nof_variables,ll)=matmul(gradients(0:IDEGFREE2,1:nof_variables,ll),transpose(EIGVL(1:nof_variables,1:nof_variables)))

		    END DO

	  ELSE

	  DO LL=1,IELEM_L(N,I)%ADMIS
	      gradients(1:IELEM_L(N,I)%idegfree,:,ll) = ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(LL,1:IELEM_L(N,I)%idegfree,1:nof_variables)
	  END DO

	  DO LL=1,IELEM_L(N,I)%ADMIS


! 	      call DGEMM ('N','T',IELEM_L(N,I)%idegfree+1,nof_variables,nof_variables,&
! 			ALPHA,gradients(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll),IELEM_L(N,I)%idegfree+1,&
!         EIGVL(1:nof_variables,1:nof_variables),nof_variables,BETA,&
!         gradients_eigvlt(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll),IELEM_L(N,I)%idegfree+1)

        gradients_eigvlt(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll)=matmul(gradients(0:IELEM_L(N,I)%idegfree,1:nof_variables,ll),transpose(EIGVL(1:nof_variables,1:nof_variables)))



	  END DO
    END IF

       IF (EES.NE.5)THEN

		      DO LL=1,IELEM_L(N,I)%ADMIS;do k=0,IELEM_L(N,I)%idegfree
			GRADCHARV(1:nof_variables,LL,k)=gradients_eigvlt(k,1:nof_variables,ll)
		      end do;END DO




		      A_CHAR(1:IELEM_L(N,I)%idegfree,1:nof_variables,1:IELEM_L(N,I)%ADMIS)=gradients_eigvlt(1:IELEM_L(N,I)%idegfree,1:nof_variables,1:IELEM_L(N,I)%ADMIS)
		      B_CHAR(1:IELEM_L(N,I)%idegfree,1:IELEM_L(N,I)%idegfree)=ILOCAL_RECON3_L(I)%INDICATOR(1:IELEM_L(N,I)%idegfree,1:IELEM_L(N,I)%idegfree)
		      CALL diag_At_B_A(ICONSIDERED,A_CHAR,B_CHAR,X_CHAR, IELEM_L, ILOCAL_RECON3_L)
!
			SMOOTHINDICATOR(1:nof_variables,1:IELEM_L(N,I)%ADMIS,L,1)=X_CHAR(1:nof_variables,1:IELEM_L(N,I)%ADMIS)
        ELSE

         DO LL=1,1;do k=0,IELEM_L(N,I)%idegfree
         GRADCHARV(1:nof_variables,LL,k)=gradients_eigvlt(k,1:nof_variables,ll)
       end do;END DO
         DO LL=2,IELEM_L(N,I)%ADMIS;do k=0,IDEGFREE2
         GRADCHARV(1:nof_variables,LL,k)=gradients_eigvlt(k,1:nof_variables,ll)
      end do;END DO


      A_CHAR(1:IELEM_L(N,I)%idegfree,1:nof_variables,1:IELEM_L(N,I)%ADMIS)=gradients_eigvlt(1:IELEM_L(N,I)%idegfree,1:nof_variables,1:IELEM_L(N,I)%ADMIS)
      B_CHAR(1:IELEM_L(N,I)%idegfree,1:IELEM_L(N,I)%idegfree)=ILOCAL_RECON3_L(I)%INDICATOR(1:IELEM_L(N,I)%idegfree,1:IELEM_L(N,I)%idegfree)

      CALL diag_At_B_A(ICONSIDERED,A_CHAR,B_CHAR,X_CHAR, IELEM_L, ILOCAL_RECON3_L)!,                                     &

        SMOOTHINDICATOR(1:nof_variables,1:IELEM_L(N,I)%ADMIS,L,1)=X_CHAR(1:nof_variables,1:IELEM_L(N,I)%ADMIS)
        END IF

!deallocate(lamc,grad5alc,a_char,b_char,x_char,gradients,gradients_eigvlt)



end subroutine


SUBROUTINE FIND_BOUNDS(ICONSIDERED,MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX,UTEMP)
IMPLICIT NONE
INTEGER::I,L,J,K,KMAXE,IQP,NGP,IEX,IK,iq,jk
INTEGER,INTENT(IN)::ICONSIDERED
REAL::leftv(1:NOF_VARIABLES),MP_PINFl,gammal
REAL,allocatable,DIMENSION(:,:),INTENT(INOUT)::UTEMP
REAL,DIMENSION(1:NOF_VARIABLES+turbulenceequations+passivescalar),INTENT(INOUT)::MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX



utemp=zero



AVER_VARS=ZERO; SUMVARS=ZERO;  MAXVARS=ZERO

I=ICONSIDERED

            if (EXTENDED_BOUNDS.EQ.0)then	!STRICT BOUNDS


                UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)


                IF (TURBULENCEEQUATIONS.GE.1)THEN
                UTEMP(1,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=&
                U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                END IF



                    K=1
                    IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
                        DO L = 1, IELEM(N,I)%IFCA
                            K=K+1
                            UTEMP(K,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:NOF_VARIABLES)
                            IF (TURBULENCEEQUATIONS.GE.1)THEN
                            UTEMP(K,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=&
                            U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                            END IF



                        END DO
                    END IF

                        IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
                            DO L=1,IELEM(N,I)%IFCA
                                IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                            K=K+1
                                            UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)

                                            IF (TURBULENCEEQUATIONS.GE.1)THEN
                                            UTEMP(K,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=&
                                            U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                                            END IF


                                            ELSE
                                            !NOT PERIODIC ONES IN MY CPU
                                            END IF
                                        ELSE
                                            K=K+1
                                            UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
                                            IF (TURBULENCEEQUATIONS.GE.1)THEN
                                            UTEMP(K,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=&
                                            U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                                            END IF
                                        END IF
                                ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS

                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                        K=K+1
                                        UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                                        (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)

                                        IF (TURBULENCEEQUATIONS.GE.1)THEN
                                            UTEMP(K,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=&
                                            IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                                        (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
                                            END IF





                                        END IF
                                    ELSE

                                    K=K+1
                                    UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                                    (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)



                                    IF (TURBULENCEEQUATIONS.GE.1)THEN
                                    UTEMP(K,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                                    (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)


                                     END IF



                                    END IF

                                END IF

                        END DO
                        END IF

                        IF (WENWRT.EQ.3)THEN
                         do jk=1,k
                         leftv(1:nof_Variables)=utemp(jk,1:nof_Variables)
                         call CONS2PRIM(N,leftv,MP_PINFl,gammal)
                         utemp(jk,1:nof_Variables)=leftv(1:nof_Variables)
                         end do
                        END IF



         end if

         if (EXTENDED_BOUNDS.EQ.1)then   !extended bounds


                        k=0
                            if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
                                    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
                                        k=k+1

                                        UTEMP(k,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)

                                         IF (TURBULENCEEQUATIONS.GE.1)THEN
                                        UTEMP(k,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)

                                         END IF

                                    END DO
                            ELSE
                                    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
                                        k=k+1

                                        IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
                                        UTEMP(k,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)

                                        IF (TURBULENCEEQUATIONS.GE.1)THEN
                                        UTEMP(k,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)

                                         END IF




                                        else
                                        UTEMP(k,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)



                                        IF (TURBULENCEEQUATIONS.GE.1)THEN
                                         UTEMP(k,NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

                                         END IF







                                        end if





                                    end do
                            END IF
                                IF (WENWRT.EQ.3)THEN
                             do jk=1,k
                             leftv(1:nof_Variables)=utemp(jk,1:nof_Variables)
                            CALL cons2prim(N,leftv,MP_PINFl,gammal)
                             utemp(jk,1:nof_Variables)=leftv(1:nof_Variables)
                             end do
                                end if

                end if



                UTMIN=ZERO;UTMAX=ZERO
                    DO IEX=1,NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR

                        UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
                        UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))

                        DO IK=2,K
                            SUMVARS(IEX)=SUMVARS(IEX)+ABS(UTEMP(IK,IEX)-UTEMP(1,IEX))
                        END DO
                        DO IK=1,K
                        AVER_VARS(IEX)=AVER_VARS(IEX)+UTEMP(IK,IEX)

                        END DO
                        AVER_VARS(IEX)=AVER_VARS(IEX)/K


                        DO IK=1,K
                        MAXVARS(IEX)=MAX(MAXVARS(IEX),ABS(UTEMP(IK,IEX)))
                        END DO
                    END DO





END SUBROUTINE FIND_BOUNDS




subroutine COMPUTE_MUSCL_RECONSTRUCTION(ICONSIDERED,UTMIN,UTMAX,UTEMP)
!> @brief
!> Subroutine for computing unlimited reconstructed solution
REAL,ALLOCATABLE,DIMENSION(:,:,:)::USOL,PSI
REAL::AX,AY,AZ,MP_PINFl,gammal,limvbg
INTEGER,INTENT(IN)::ICONSIDERED
INTEGER::ICD,ICOMPWRT,I,NGP,L,IEX,iqp
REAL,dimension(1:NOF_VARIABLES)::LEFTV
REAL,allocatable,DIMENSION(:,:),intent(in)::UTEMP
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::UTMIN,UTMAX
REAL,ALLOCATABLE,DIMENSION(:,:)::CONSMATRIX,GRADSSL,RESSOLUTION,GRADSSL2,RESSOLUTION2
REAL,allocatable,dimension(:)::SLOPE




I=ICONSIDERED
ILOCAL_RECON3(ICONSIDERED)%ULEFT(:,:,:)=ZERO




ALLOCATE(SLOPE(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(USOL(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2))
ALLOCATE(PSI(NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2))
ALLOCATE(CONSMATRIX(1:6*NUMBEROFPOINTS2,1:IDEGFREE))
ALLOCATE(GRADSSL(1:IDEGFREE,1:NOF_VARIABLES))
ALLOCATE(RESSOLUTION(1:6*NUMBEROFPOINTS2,1:NOF_VARIABLES))


IF (TURBULENCEEQUATIONS.GE.1)THEN
ALLOCATE(GRADSSL2(1:IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(RESSOLUTION2(1:6*NUMBEROFPOINTS2,1:TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF


iCOMPWRT=0

USOL(:,:,:)=ZERO
PSI=ZERO
icd=0
DO L=1,IELEM(N,I)%IFCA	!faces2

                                    if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if


            do NGP=1,iqp			!for gqp
                AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				if (DIMENSIONA.eq.3) then
				az= ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				end if

            icd=icd+1

                    if (DIMENSIONA.eq.3) then
				  CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE,ICOMPWRT, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)
				  else
				  CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2D(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE,ICOMPWRT, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)
				  end if
            end do
end do

GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)
! CALL DGEMM('N','N',ICD,nof_variables,IELEM(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),ICD,&
!                     GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),&
!                     IELEM(N,I)%IDEGFREE,BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),ICD)


                 RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES))



IF (TURBULENCEEQUATIONS.GE.1)THEN
GRADSSL2(1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTS2(1,1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
! CALL DGEMM('N','N',ICD,TURBULENCEEQUATIONS+PASSIVESCALAR,IELEM(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),ICD,&
!                     GRADSSL2(1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR),&
!                     IELEM(N,I)%IDEGFREE,BETA,RESSOLUTION2(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR),ICD)

            RESSOLUTION2(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=matmul(CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),GRADSSL2(1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))


END IF



 ICD=0              !initialise counter
 DO L=1,IELEM(N,I)%IFCA     !loop all faces
                                    if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if
            do NGP=1,iqp        !all gaussian quadrature points
			     ICD=ICD+1

             IF (WENWRT.EQ.3)THEN
                LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:nof_Variables)
                        CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
				  USOL(1:NOF_VARIABLES,L,Ngp)=((LEFTV(1:NOF_VARIABLES)+RESSOLUTION(icd,1:NOF_VARIABLES)))
				  ELSE
				  USOL(1:NOF_VARIABLES,L,Ngp)=((U_C(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(icd,1:NOF_VARIABLES)))
				  END IF


				  IF (TURBULENCEEQUATIONS.GE.1)THEN

				  USOL(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,L,Ngp)=((U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+RESSOLUTION2(icd,1:TURBULENCEEQUATIONS+PASSIVESCALAR)))


				  END IF

				  END DO
			    END DO



 DO L=1,IELEM(N,I)%IFCA	!faces2
                                    if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if
			      do NGP=1,iqp
				    DO iex=1,nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR
				CALL SLOPE_LIMITERS(N,I,IEX,L,NGP,UTMIN,UTMAX,UTEMP,USOL,PSI)
				    end do


			      END DO
			 END DO


             do iex=1,NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
                                    if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if
			      do NGP=1,iqp


				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 SLOPE(IEX)=LIMVBG
			 IF (IEX.EQ.1)THEN
			 ielem(n,i)%wcx(1)=LIMVBG
			 END IF

			  IF (DG.EQ.1)THEN
				ILOCAL_RECON6(I)%DG2FV(1:IELEM(N,I)%IDEGFREE,IEX)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTS(1,IELEM(N,I)%IDEGFREE,iex)*SLOPE(IEX)
				END IF

                end do

                 ILOCAL_RECON3(I)%ULEFT(:,:,:)=ZERO
                 IF (TURBULENCEEQUATIONS.GE.1)THEN
                 ILOCAL_RECON3(I)%ULEFTTURB(:,:,:)=ZERO

                 END IF



                 DO L=1,IELEM(N,I)%IFCA	!faces2
                                    if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if

			      do NGP=1,iqp

                        IF (WENWRT.EQ.3)THEN
                        LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:nof_Variables)
                        CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
                        USOL(1:nof_Variables,l,Ngp)=USOL(1:nof_Variables,l,Ngp)-LEFTV(1:NOF_VARIABLES)
                        ELSE
                        USOL(1:nof_Variables,l,Ngp)=USOL(1:nof_Variables,l,Ngp)-U_C(I)%VAL(1,1:nof_Variables)
                        END IF
                        IF (TURBULENCEEQUATIONS.GE.1)THEN
				        USOL(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)=USOL(nof_Variables+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)&
				        -U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                        END IF



				    CALL EXTRAPOLATE_BOUND_MUSCL(USOL,IEX,L,NGP,I,SLOPE)
			      END DO
			 END DO


DEALLOCATE(SLOPE)
DEALLOCATE(USOL)
DEALLOCATE(PSI)
DEALLOCATE(CONSMATRIX)
DEALLOCATE(GRADSSL)
DEALLOCATE(RESSOLUTION)


IF (TURBULENCEEQUATIONS.GE.1)THEN
DEALLOCATE(GRADSSL2)
DEALLOCATE(RESSOLUTION2)
END IF




END subroutine COMPUTE_MUSCL_RECONSTRUCTION





SUBROUTINE MUSCL(N)
!> @brief
!> Subroutine for MUSCL type reconstruction
 IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,IEX,IEUL,JX,LX,KMAXE,iq,LL,NGP,NND,IQP,idummy,ii,icd,ICONSIDERED
REAL::UMIN,UMAX,PSITOT,ADDC,DIVG0,LIMVBG,tempxx
REAL,allocatable,DIMENSION(:,:)::UTEMP
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX

allocate(utemp(IMAXDEGFREE+1,1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR))





KMAXE=XMPIELRANK(N)


!$OMP DO
	DO II=1,NOF_INTERIOR
	I=EL_INT(II)
	ICONSIDERED=I
        IF (((ielem(n,i)%TROUBLED.EQ.1).AND.(ielem(n,i)%REDUCE.EQ.1)).OR.((ielem(n,i)%FULL.EQ.0).AND.(ielem(n,i)%TROUBLED.EQ.1)))THEN
            IF (IELEM(N,I)%RECALC.GT.0)THEN
                if (ADDA.EQ.1)THEN
                    CALL ADDA_FILTER(N,ICONSIDERED,  IELEM, ILOCAL_RECON3, ILOCAL_RECON5, U_C, INTEG_BASIS, integ_basis_dg)
                END IF

                    CALL FIND_BOUNDS(ICONSIDERED,MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX,UTEMP)

                    CALL COMPUTE_MUSCL_RECONSTRUCTION(ICONSIDERED,UTMIN,UTMAX,UTEMP)

            END IF
        END IF
    END DO
!$OMP END DO

!$OMP DO
DO II=1,NOF_bounded
	I=EL_BND(II)
	ICONSIDERED=I

     if (ADDA.EQ.1)THEN
      CALL ADDA_FILTER(N,ICONSIDERED, IELEM, ILOCAL_RECON3, ILOCAL_RECON5, U_C, INTEG_BASIS, integ_basis_dg)
      END IF

       IF (((ielem(n,i)%TROUBLED.EQ.1).AND.(ielem(n,i)%REDUCE.EQ.1)).OR.((ielem(n,i)%FULL.EQ.0).AND.(ielem(n,i)%TROUBLED.EQ.1)))THEN
            IF (IELEM(N,I)%RECALC.GT.0)THEN
                if (ADDA.EQ.1)THEN
                    CALL ADDA_FILTER(N,ICONSIDERED, IELEM, ILOCAL_RECON3, ILOCAL_RECON5, U_C, INTEG_BASIS, integ_basis_dg)
                END IF

                    CALL FIND_BOUNDS(ICONSIDERED,MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX,UTEMP)

                    CALL COMPUTE_MUSCL_RECONSTRUCTION(ICONSIDERED,UTMIN,UTMAX,UTEMP)

            END IF
        END IF


      END DO
!$OMP END DO

DEALLOCATE(UTEMP)

END SUBROUTINE MUSCL





SUBROUTINE SOLUTIONTRIAV2(N)
!> @brief
!> Subroutine for extrapolating the unlimited reconstructed values for diffusive fluxes in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,M,PPP,IEUL,IEX,IHGT,IHGJ,KMAXE,DECOMF,ICNN,IQDR,NVAR,idummy,iqp,nnd,ngp,icd,icompwrt,ICONSIDERED
REAL::RAA1,RAA2,PAA1,PAA2,ax,ay,az
REAL::SOLX
real,dimension(1:dimensiona)::ugradloc
real,dimension(1:dimensionA,1:dimensionA)::ainvjt
real,allocatable,dimension(:)::gradtem
real,allocatable,dimension(:,:)::XXDER,YYDER,ZZDER


KMAXE=XMPIELRANK(N);


allocate(xxder(1:idegfree,1:NUMBEROFPOINTS2))
allocate(yyder(1:idegfree,1:NUMBEROFPOINTS2))
allocate(zzder(1:idegfree,1:NUMBEROFPOINTS2))
allocate(gradtem(1:idegfree))










!$OMP DO
DO I=1,kmaxe
	ICONSIDERED=I





ILOCAL_RECON3(I)%ULEFTV(:,:,:,:)=zero;
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
ILOCAL_RECON3(I)%ULEFTTURBV(:,:,:,:)=zero;ILOCAL_RECON3(I)%ULEFTTURB(:,:,:)=zero;
END IF



	    DO IHGT=1,DIMENSIONA;DO IHGJ=1,DIMENSIONA
		AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
	    END DO;END DO






	DO l=1,IELEM(N,I)%IFCA;IDUMMY=0
                                                    IF (DIMENSIONA.EQ.3)THEN
                                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                                    iqp=qp_quad;
                                                    else
                                                    iqp=qp_triangle;
                                                    end if
                                                    ELSE
                                                    iqp=qp_LINE;
                                                    END IF
		 ICD=0
                do NGP=1,iqp			!for gqp

				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1);
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2);
				IF (DIMENSIONA.EQ.3)THEN
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				end if

				icd=icd+1
                        if (dimensiona.eq.3)then
                        DO K=1,IELEM(N,I)%IDEGFREE
                            IF (POLY.EQ.1) THEN
                                XXDER(K,ICD)=DFX(AX,AY,AZ,K,i);  YYDER(K,ICD)=DFY(AX,AY,AZ,K,i);  ZZDER(K,ICD)=DFZ(AX,AY,AZ,K,i)
                            END IF
                            IF (POLY.EQ.2) THEN
                                XXDER(K,ICD)=DLX(AX,AY,AZ,K,i);  YYDER(K,ICD)=DLY(AX,AY,AZ,K,i);  ZZDER(K,ICD)=DLZ(AX,AY,AZ,K,i)
                            END IF
                            IF (POLY.EQ.4) THEN
                                XXDER(K,ICD)=TL3DX(AX,AY,AZ,K,i);  YYDER(K,ICD)=TL3DY(AX,AY,AZ,K,i);  ZZDER(K,ICD)=TL3DZ(AX,AY,AZ,K,i)
                            END IF
                        END DO
                        ELSE
                         DO K=1,IELEM(N,I)%IDEGFREE
                            IF (POLY.EQ.4)THEN
							xXDER(K,icd)=TL2dX(AX,AY,K,i);  yYDER(K,icd)=TL2dY(AX,AY,K,i);
							ELSE
						    xXDER(K,icd)=DF2dX(AX,AY,K,i);  yYDER(K,icd)=DF2dY(AX,AY,K,i);
							END IF

                        END DO
                        end if


				end do
			    ICD=0
			     do NGP=1,iqp
                                    icd=icd+1

				SELECT CASE(IELEM(N,I)%GGS)

				CASE(0)

                        !TURBULENCE FIRST
                        IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN

                            if (icoupleturb.eq.0)then
                                DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
                                ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
                                END DO

                            end if


                            DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR

                                    GRADTEM(1:IELEM(N,I)%IDEGFREE)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTSTURB(1,1:IELEM(N,I)%IDEGFREE,NVAR)

                                    UGRADLOC = ZERO


                            UGRADLOC(1)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),XXDER(1:IELEM(N,I)%IDEGFREE,ICD))
                            UGRADLOC(2)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),YYDER(1:IELEM(N,I)%IDEGFREE,ICD))
                            if (dimensiona.eq.3)then
                            UGRADLOC(3)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),ZZDER(1:IELEM(N,I)%IDEGFREE,ICD))
                            end if



                                ILOCAL_RECON3(I)%ULEFTTURBV(1:dimensiona,NVAR,L,NGP) = MATMUL(AINVJT(1:dimensiona,1:dimensiona),UGRADLOC(1:dimensiona))


                            END DO
                        END IF

                    !now temperature
					GRADTEM(1:IELEM(N,I)%IDEGFREE)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTSTEMP(1:IELEM(N,I)%IDEGFREE)
!
					UGRADLOC = ZERO






                UGRADLOC(1)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),XXDER(1:IELEM(N,I)%IDEGFREE,ICD))
                UGRADLOC(2)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),YYDER(1:IELEM(N,I)%IDEGFREE,ICD))
                 if (dimensiona.eq.3)then
                UGRADLOC(3)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),ZZDER(1:IELEM(N,I)%IDEGFREE,ICD))
                end if


					  ILOCAL_RECON3(I)%ULEFTV(1:dimensiona,1,L,NGP) = MATMUL(AINVJT(1:dimensiona,1:dimensiona),UGRADLOC(1:dimensiona))


                    !now velocities
				  DO IEX=1,dimensiona
!
					GRADTEM(1:IELEM(N,I)%IDEGFREE)=ILOCAL_rECON5(ICONSIDERED)%VELOCITYDOF(IEX,1:IELEM(N,I)%IDEGFREE)
!
					 UGRADLOC = ZERO


					     UGRADLOC(1)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),XXDER(1:IELEM(N,I)%IDEGFREE,ICD))
                UGRADLOC(2)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),YYDER(1:IELEM(N,I)%IDEGFREE,ICD))
                if (dimensiona.eq.3)then
                UGRADLOC(3)=DOT_PRODUCT(GRADTEM(1:IELEM(N,I)%IDEGFREE),ZZDER(1:IELEM(N,I)%IDEGFREE,ICD))
                end if


					   ILOCAL_RECON3(I)%ULEFTV(1:dimensiona,IEX+1,L,NGP) = MATMUL(AINVJT(1:dimensiona,1:dimensiona),UGRADLOC(1:dimensiona))





!
				  END DO






				CASE(1)


				  !TURBULENCE FIRST
				IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN


                        if (icoupleturb.eq.0)then
                            DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
                        ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
                            END DO
                        end if


                        DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR

                            ILOCAL_RECON3(I)%ULEFTTURBV(1:dimensiona,NVAR,L,NGP)=ILOCAL_RECON3(I)%GRADs(dimensiona+1+NVAR,1:dimensiona)

                            END DO
				END IF

				!MEAN FLOW GRADIENTS

				ILOCAL_RECON3(I)%ULEFTV(1:dimensiona,1,L,NGP) = ILOCAL_RECON3(I)%GRADs(dimensiona+1,1:dimensiona)

				DO IEX=1,dimensiona
				    ILOCAL_RECON3(I)%ULEFTV(1:dimensiona,IEX+1,L,NGP) = ILOCAL_RECON3(I)%GRADs(IEX,1:dimensiona)
				END DO


                END  SELECT


			END DO
	  END DO




                    IF (IELEM(N,ICONSIDERED)%GGS.EQ.0)THEN

                    CALL COMPUTE_GRADIENTS_CENTER(N,ICONSIDERED)

                    END IF







	  end do
!$OMP END DO

deallocate(xxder)
deallocate(yyder)
deallocate(zzder)
deallocate(gradtem)






END SUBROUTINE SOLUTIONTRIAV2










SUBROUTINE LEAST_SQUARES(N, IELEM_L,ILOCAL_RECON3_L,ILOCAL_RECON5_L,U_C_L, U_CT_L,IEXSOLHIR_L)
!> @brief
!> Subroutine For COMPUTING LEAST SQUARES RECONSTRUCTION
IMPLICIT NONE
INTEGER,INTENT(IN)::N
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON5_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C_L, U_CT_L
TYPE(EXCHANGE_SOLHI),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXSOLHIR_L
INTEGER::ICONSIDERED,II,I


#ifdef LEASTSQUARES_GPU_KERNEL
!$OMP target teams distribute parallel do
#else
!$OMP DO
#endif
DO II=1,NOF_INTERIOR;

I=EL_INT(II)
ICONSIDERED=I
 CALL ALLGRADS_INNER(N,ICONSIDERED,IELEM_L,ILOCAL_RECON3_L,ILOCAL_RECON5_L,U_C_L,U_CT_L,IEXSOLHIR_L)
END DO
#ifdef LEASTSQUARES_GPU_KERNEL
!$OMP end target teams distribute parallel do
#else
!$OMP END DO
#endif


!$OMP DO
DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I
CALL ALLGRADS_MIX(N,I)
END DO
!$OMP END DO

END SUBROUTINE LEAST_SQUARES





! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED TO EMPLOY THE LEAST SQARES LINEAR INTERPOLATION!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!FOR DETERMINING THE SLOPES OF EACH CELL IN EACH DIRECTION!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!IN A WEIGHTED AVERAGE WAY WITH THE INVERSE DISTANCE !!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PIECEWISE_CONSTANT(N)
!> @brief
!> Subroutine For first-order scheme
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,KMAXE,IEX,l
KMAXE=XMPIELRANK(N)



!$OMP DO
	DO I=1,KMAXE
	DO IEX=1,nof_variables
 	ILOCAL_RECON3(I)%ULEFT(IEX,:,:)=U_C(I)%VAL(1,IEX)
	END DO
	
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	DO IEX=1,turbulenceequations+passivescalar
	ILOCAL_RECON3(I)%ULEFTTURB(IEX,:,:)=U_Ct(I)%VAL(1,IEX)
	END DO
	end if
	END DO
!$OMP END DO 
 
END SUBROUTINE PIECEWISE_CONSTANT






SUBROUTINE LINEAR_SCHEME(n)
!> @brief
!> Subroutine for linear type reconstruction
 IMPLICIT NONE
 INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,IEX,IEUL,JX,LX,KMAXE,iq,LL,NGP,NND,IQP,idummy,ii,icd,ICONSIDERED
REAL::UMIN,UMAX,PSITOT,ADDC,DIVG0,LIMVBG,tempxx
KMAXE=XMPIELRANK(N)


!$OMP DO
	DO II=1,NOF_INTERIOR
	I=EL_INT(II)
	ICONSIDERED=I

					CALL COMPUTE_LINEAR_RECONSTRUCTION(ICONSIDERED)
    END DO
!$OMP END DO

!$OMP DO
DO II=1,NOF_bounded
	I=EL_BND(II)
	ICONSIDERED=I

                    CALL COMPUTE_LINEAR_RECONSTRUCTION(ICONSIDERED)

      END DO
!$OMP END DO



END SUBROUTINE LINEAR_SCHEME







subroutine COMPUTE_LINEAR_RECONSTRUCTION(ICONSIDERED)
!> @brief
!> Subroutine for computing unlimited reconstructed solution
REAL,ALLOCATABLE,DIMENSION(:,:,:)::USOL,PSI
REAL::AX,AY,AZ,MP_PINFl,gammal,limvbg
INTEGER,INTENT(IN)::ICONSIDERED
INTEGER::ICD,ICOMPWRT,I,NGP,L,IEX,IQP
REAL,dimension(1:NOF_VARIABLES)::LEFTV
REAL,ALLOCATABLE,DIMENSION(:,:)::CONSMATRIX,GRADSSL,RESSOLUTION,GRADSSL2,RESSOLUTION2



I=ICONSIDERED

ILOCAL_RECON3(ICONSIDERED)%ULEFT(:,:,:)=ZERO


ALLOCATE(USOL(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2))
ALLOCATE(CONSMATRIX(1:6*NUMBEROFPOINTS2,1:IDEGFREE))
ALLOCATE(GRADSSL(1:IDEGFREE,1:NOF_VARIABLES))
ALLOCATE(RESSOLUTION(1:6*NUMBEROFPOINTS2,1:NOF_VARIABLES))
allocate(psi(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2))


IF (TURBULENCEEQUATIONS.GE.1)THEN
ALLOCATE(GRADSSL2(1:IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(RESSOLUTION2(1:6*NUMBEROFPOINTS2,1:TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF

iCOMPWRT=0

USOL(:,:,:)=ZERO
PSI=ZERO
icd=0


DO L=1,IELEM(N,I)%IFCA	!faces2

									if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if


            do NGP=1,iqp			!for gqp
                AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				if (DIMENSIONA.eq.3) then
				az= ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				end if

            icd=icd+1

                    if (DIMENSIONA.eq.3) then
				  CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE,ICOMPWRT, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)
				  else
				  CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2D(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE,ICOMPWRT, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)
				  end if
            end do
end do

GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)



! CALL DGEMM('N','N',ICD,nof_variables,IELEM(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),ICD,&
!                     GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),&
!                     IELEM(N,I)%IDEGFREE,BETA,RESSOLUTION(1:ICD,1:NOF_vARIABLES),ICD)


        RESSOLUTION(1:ICD,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES))



IF (TURBULENCEEQUATIONS.GE.1)THEN
GRADSSL2(1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_rECON5(ICONSIDERED)%GRADIENTS2(1,1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
! CALL DGEMM('N','N',ICD,TURBULENCEEQUATIONS+PASSIVESCALAR,IELEM(N,I)%IDEGFREE,ALPHA,&
!                     CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),ICD,&
!                     GRADSSL2(1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR),&
!                     IELEM(N,I)%IDEGFREE,BETA,RESSOLUTION2(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR),ICD)

     RESSOLUTION2(1:ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR)= matmul(CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),GRADSSL2(1:IELEM(N,I)%IDEGFREE,1:TURBULENCEEQUATIONS+PASSIVESCALAR))


END IF



 ICD=0              !initialise counter

 ICD=0              !initialise counter
 DO L=1,IELEM(N,I)%IFCA     !loop all faces
                                    if (DIMENSIONA.eq.3)then
                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                    iqp=qp_quad
                                    else
                                    iqp=qp_triangle
                                    end if
                                    else
                                     iqp=qp_LINE
                                    end if
            do NGP=1,iqp        !all gaussian quadrature points
			     ICD=ICD+1
					USOL(1:NOF_VARIABLES,L,NGP)=RESSOLUTION(ICD,1:NOF_VARIABLES)
					IF (turbulenceequations.GE.1)THEN
					USOL(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,L,NGP)=RESSOLUTION2(ICD,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					END IF
					CALL EXTRAPOLATE_BOUND_LINEAR(USOL,IEX,L,NGP,I)



			END DO
END DO



DEALLOCATE(USOL)
DEALLOCATE(PSI)
DEALLOCATE(CONSMATRIX)
DEALLOCATE(GRADSSL)
DEALLOCATE(RESSOLUTION)


IF (TURBULENCEEQUATIONS.GE.1)THEN
DEALLOCATE(GRADSSL2)
DEALLOCATE(RESSOLUTION2)
END IF



END SUBROUTINE COMPUTE_LINEAR_RECONSTRUCTION








SUBROUTINE ARBITRARY_ORDER(N)
!> @brief
!> Subroutine controlling the reconstruction in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I


KMAXE=XMPIELRANK(N)
ielem(n,1:kmaxe)%REDUCE=0



  CALL LEAST_SQUARES(N,IELEM, ILOCAL_RECON3,ILOCAL_RECON5,U_C,U_CT,IEXSOLHIR)

	
 SELECT CASE(IWENO)
 
 
  CASE(1)

  CALL WENOWEIGHTS(N, IELEM, ILOCAL_RECON3, ILOCAL_RECON5, ILOCAL_RECON6, IBOUND, U_C, U_Ct, INTEG_BASIS,integ_basis_dg, INODER4, IEXSOLHIR)
  CALL CHECKSOL(N)
  CALL MUSCL(N)
  CALL CHECKSOLX(N)

  CASE(-1)

  CALL MUSCL(N)		
  CALL CHECKSOLX(N)

 
 
	CASE(0)
	IF (FIRSTORDER.EQ.1)THEN
	CALL PIECEWISE_CONSTANT(N)
	ELSE

	CALL LINEAR_SCHEME(N)
	CALL CHECKSOLX(N)
	END IF



	END SELECT




	if (Itestcase.eq.4)then

	call solutiontriav2(n)
	end if


 

 
 
 

END SUBROUTINE ARBITRARY_ORDER



SUBROUTINE VISCOUS_DG_GGS(N)
!> @brief
!> Subroutine controlling the reconstruction in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I,IEX,l,ngp,iqp,ICONSIDERED


KMAXE=XMPIELRANK(N)

IF (DIMENSIONA.EQ.3)THEN
!$OMP DO
DO I=1,KMAXE
ICONSIDERED=I


			DO L=1,IELEM(N,I)%IFCA

						if (ielem(n,i)%types_faces(L).eq.5)then
							iqp=qp_quad
						else
							iqp=QP_TRIANGLE
						end if

						do NGP=1,iqp
							ILOCAL_RECON3(I)%ULEFTV(1:3,1,l,ngp) = ILOCAL_RECON3(I)%GRADs(4,1:3)
								DO IEX=1,3
									ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,l,ngp) = ILOCAL_RECON3(I)%GRADs(IEX,1:3)
								END DO
						end do
			end do




END DO
!$OMP END DO



ELSE
!$OMP DO
DO I=1,KMAXE
ICONSIDERED=I



		DO L=1,IELEM(N,I)%IFCA

                IQP=QP_LINE_N

                DO NGP=1,IQP



			ILOCAL_RECON3(I)%ULEFTV(1:2,1,l,ngp) = ILOCAL_RECON3(I)%GRADs(3,1:2)
				DO IEX=1,2
				    ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,l,ngp) = ILOCAL_RECON3(I)%GRADs(IEX,1:2)
				END DO
				end do
				end do



END DO
!$OMP END DO












END IF


END SUBROUTINE VISCOUS_DG_GGS

























SUBROUTINE CHECKSOL(N)
!> @brief
!> Subroutine for checking the reconstructed solution
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,NGP,iqp,iex
INTEGER::REDUCE1,kmaxe,indx
real::jump_cond
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
KMAXE=XMPIELRANK(N)
jump_cond=0.85






IF (ITESTCASE.GE.3)THEN

!$OMP DO
	DO I=1,KMAXE	!ALL ELEMENTS
        IELEM(N,I)%REDUCE=0;REDUCE1=0



        if (ielem(n,i)%troubled.eq.1)then

        IF (IELEM(N,I)%FULL.EQ.0)THEN
			IELEM(N,I)%REDUCE=1
        END IF


        IF (IELEM(N,I)%FULL.EQ.1)THEN
		
		  DO L=1,IELEM(N,I)%IFCA	!faces2
				  if (dimensiona.eq.3)then
				  indx=5
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      eLSE
			      indx=4
					 iqp=QP_LINE
			      end if

				  do NGP=1,iqp
					      
						LEFTV(1:NOF_VARIABLES)=ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)
						RIGHTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
						CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						
						
															IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(jump_cond*RIGHTV(1))))then
																	REDUCE1=1
																IELEM(N,I)%REDUCE=1
															end if

															IF (((ABS(LEFTV(indx)-RIGHTV(indx))).GE.(jump_cond*RIGHTV(indx))))then
																	REDUCE1=1
																IELEM(N,I)%REDUCE=1
															end if

															IF ((LEFTV(indx).LT.0.0).OR.(LEFTV(1).LT.0.0))then
																	REDUCE1=1
																IELEM(N,I)%REDUCE=1
															end if


				
					
				  END DO
		END DO	
		
		
		if (ielem(n,i)%hybrid.eq.1)then
		reduce1=1
		IELEM(N,I)%REDUCE=1
		end if
		
!
		
		
		END IF
		
		
		end if
	end do
!$OMP END DO		
		

		
END IF		
		
		

END SUBROUTINE CHECKSOL

 
SUBROUTINE CHECKSOLX(N)
!> @brief
!> Subroutine for checking the reconstructed solution
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,NGP,iqp,iex
INTEGER::REDUCE1,kmaxe,indx
real::jump_cond
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
KMAXE=XMPIELRANK(N)
jump_cond=0.95




IF (ITESTCASE.GE.3)THEN

!$OMP DO
	DO I=1,KMAXE	!ALL ELEMENTS
        

			REDUCE1=0

        if (ielem(n,i)%troubled.eq.1)then
        

						DO L=1,IELEM(N,I)%IFCA	!faces2
								if (dimensiona.eq.3)then
									indx=5
									if (ielem(n,i)%types_faces(L).eq.5)then
										iqp=qp_quad
									else
										iqp=QP_TRIANGLE
									end if
									eLSE
									indx=4
										iqp=QP_LINE
									end if
										do NGP=1,iqp

												LEFTV(1:NOF_VARIABLES)=ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)
												RIGHTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
												CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)

															IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(jump_cond*RIGHTV(1))))then
																	REDUCE1=1
																IELEM(N,I)%REDUCE=1
															end if

															IF (((ABS(LEFTV(indx)-RIGHTV(indx))).GE.(jump_cond*RIGHTV(indx))))then
																	REDUCE1=1
																IELEM(N,I)%REDUCE=1
															end if

															IF ((LEFTV(indx).LT.0.0).OR.(LEFTV(1).LT.0.0))then
																	REDUCE1=1
																IELEM(N,I)%REDUCE=1
															end if




										END DO
						END DO


					if (ielem(n,i)%hybrid.eq.1)then
					reduce1=1
					IELEM(N,I)%REDUCE=1
					end if

					IF (REDUCE1.EQ.1)THEN
						do iex=1,NOF_VARIABLES
						ILOCAL_RECON3(I)%ULEFT(iex,:,:)=u_c(i)%val(1,iex)

						end do

						if (dg.eq.1)then
						ILOCAL_RECON6(I)%DG2FV(1:IELEM(N,I)%IDEGFREE,:)=zero

						end if

					end if


					if (turbulence.eq.1)then
							if (icoupleturb.eq.1)then
							REDUCE1=0
								DO L=1,IELEM(N,I)%IFCA	!faces2

										if (dimensiona.eq.3)then
										indx=5
										if (ielem(n,i)%types_faces(L).eq.5)then
											iqp=qp_quad
										else
											iqp=QP_TRIANGLE
										end if
										eLSE
										indx=4
											iqp=QP_LINE
										end if

										do NGP=1,iqp
											leftv(1)=ILOCAL_RECON3(I)%ULEFTTURB(1,L,ngp)
											RIGHTV(1)=U_Ct(I)%VAL(1,1)
												IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(0.6*RIGHTV(1))).OR.(leftv(1).le.zero)) THEN
													REDUCE1=1
												END IF
										end do
								end do

								if (ielem(n,i)%hybrid.eq.1)then
								reduce1=1
								IELEM(N,I)%REDUCE=1
								end if

								IF (REDUCE1.EQ.1)THEN
								do iex=1,1
								ILOCAL_RECON3(I)%ULEFTTURB(1,:,:)=u_ct(i)%val(1,1)
								end do
								end if
							end if
					end if





		end if
	end do
!$OMP END DO




END IF



END SUBROUTINE CHECKSOLX
 
 












SUBROUTINE SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S,UTMIN,UTMAX,UTEMP,USOL,PSI)
!> @brief
!> Subroutine MUSCL slope limiters
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONS_E,FACEX,ICONS_S,ICONSIDERED
REAL,allocatable,DIMENSION(:,:),INTENT(IN)::UTEMP
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR),INTENT(IN)::UTMIN,UTMAX
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN)::USOL
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::PSI
INTEGER::I,L,IEX,NGP
REAL::D2,SFD,EPSI2,DMIN,DPLUS,KAPPA_VEN,psi2,sig_1,delu,y_fun,s_y,pol_MOG
IEX=ICONS_E
L=FACEX
I=ICONSIDERED
NGP=ICONS_S
KAPPA_VEN=10.0D0
psi2=zero

					  D2=USOL(IEX,L,NGP)-UTEMP(1,IEX)


					  IF (ABS(D2).LE.ZERO)THEN

					  PSI(iex,L,ngp)=1.0D0

					  ELSE IF (D2.GT.zero)THEN
					 SFD=(UTMAX(iex)-UTEMP(1,IEX))/D2

					      
					      SELECT CASE (LIMITER)

					      CASE(1)
					      PSI(iex,L,ngp) = MIN(1.0d0,SFD)		!BARTH AND JESPERSEN
					      CASE(2)
! 					      
					      pol_MOG=-((4.0d0/27.0d0)*sfd**3)+sfd
					      if (sfd.lt.1.5d0)then
					      PSI(iex,L,ngp)=pol_MOG
					      
					      else
					      PSI(iex,L,ngp)= 1.0d0
					      
					      end if
					      
					      
					      
					        	!BARTH AND JESPERSEN MICHALAK
					      
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=(2.0d0*y_fun**3)-(3.0d0*y_fun**2)+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2



                                                CASE(9)

					      pol_MOG=-((4.0d0/27.0d0)*sfd**3)+sfd
					      if (sfd.lt.1.5d0)then
					      PSI(iex,L,ngp)=pol_MOG
					      
					      else
					      PSI(iex,L,ngp)= 1.0d0
					      
					      end if
					      
					      
					      
					        	!BARTH AND JESPERSEN MICHALAK
					      
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=(2.0d0*y_fun**3)-(3.0d0*y_fun**2)+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2
! 					      
					      
					      CASE(3)
					      pol_MOG=-((4.0d0/27.0d0)*sfd**3)+sfd
					      if (sfd.lt.1.5d0)then
					      PSI(iex,L,ngp)=pol_MOG
					      
					      else
					      PSI(iex,L,ngp)= 1.0d0
					      
					      end if	!BARTH AND JESPERSEN MICHALAK EXTENDED
					      
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2
					      
					      
					      
					      
					      CASE(4)                  !VKM
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					      dmin=sign(1.0d0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMAX(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2))
								    
					delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2    
					      
					      CASE(5)
					      PSI(iex,L,ngp)= (SFD**2 + SFD) / (SFD**2 + 1.0d0)   			! VAN ALBADA
					      CASE(6)
					      PSI(iex,L,ngp)= 2.0d0*SFD / (SFD + 1.0d0) 				! VAN LEER
					      CASE(7)								!VENKATAKRISHNAN
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					      dmin=sign(1.0d0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMAX(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2))
					      
					      
					       CASE(8)								!VENKATAKRISHNAN
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					      dmin=sign(1.0d0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMAX(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2))
								    
					delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2    
								    
					      END SELECT
					  
					  ELSE
					      SFD=((UTMIN(IEX)-UTEMP(1,IEX)))/(D2)


					      SELECT CASE (LIMITER)

					    CASE(1)
					      PSI(iex,L,ngp) = MIN(1.0d0,SFD)				!MINMOD LIMITER
					      CASE(2)
					       pol_MOG=-((4.0d0/27.0d0)*sfd**3)+sfd
					       
					       
! 					       
					      if (sfd.lt.1.5d0)then
					      
					      
					      PSI(iex,L,ngp)=pol_MOG
					      
					      else
					      PSI(iex,L,ngp)= 1.0d0
					      
					      end if					!BARTH AND JESPERSEN
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2
					      
					      
! 					      
					      
					      CASE(9)
					       pol_MOG=-((4.0d0/27.0d0)*sfd**3)+sfd
					       
					       
! 					       
					      if (sfd.lt.1.5d0)then
					      
					      
					      PSI(iex,L,ngp)=pol_MOG
					      
					      else
					      PSI(iex,L,ngp)= 1.0d0
					      
					      end if					!BARTH AND JESPERSEN
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2
					      CASE(3)
					      pol_MOG=-((4.0d0/27.0d0)*sfd**3)+sfd
					      if (sfd.lt.1.5d0)then
					      PSI(iex,L,ngp)=pol_MOG
					      
					      else
					      PSI(iex,L,ngp)= 1.0d0
					      
					      end if					!BARTH AND JESPERSEN
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2
					      
					      CASE(4)
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					       dmin=sign(1.0d0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMIN(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2)) 
					      
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2 				!BARTH AND JESPERSEN (EXTENDED STENCILS)
					     
					      CASE(5)
					      PSI(iex,L,ngp)= (SFD**2 + SFD) / (SFD**2 + 1.0d0)   			! VAN ALBADA
					      CASE(6)
					      PSI(iex,L,ngp)= 2.0d0*SFD / (SFD + 1.0d0) 					! VAN LEER
					      CASE(7)								!VENKATAKRISHNAN
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					       dmin=sign(1.0d0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMIN(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2)) 
					      CASE(8)								!VENKATAKRISHNAN
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					       dmin=sign(1.0d0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMIN(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2)) 
					      
					      delu=UTMAX(IEX)-utmin(iex)
					      
					      if (delu**2.le.((kappa_ven*IELEM(N,I)%MINEDGE)**3))then
					      
					      
					      sig_1=1.0d0
					      
					      
					      else
					      
					      if ((((kappa_ven*IELEM(N,I)%MINEDGE)**3).lt.delu**2).and.(delu**2.lt.2.0d0*((kappa_ven*IELEM(N,I)%MINEDGE)**3)))then
					      
					      y_fun=((delu**2)-((kappa_ven*IELEM(N,I)%MINEDGE)**3))/((kappa_ven*IELEM(N,I)%MINEDGE)**3)
					      
					      s_y=2.0d0*y_fun**3-3.0d0*y_fun**2+1.0d0
					      
					      sig_1=s_y
					      
					      else
					      
					      
					      sig_1=0.0d0
					      
					      end if
					      end if
					      
					      
					      
					      
					      psi2=sig_1+(1.0d0-sig_1)*PSI(iex,L,ngp)
					      
					      
					      PSI(iex,L,ngp)= psi2
					      
					      END SELECT
						
						
					 
					  END IF

END SUBROUTINE SLOPE_LIMITERS


SUBROUTINE TROUBLE_INDICATOR1
IMPLICIT NONE
INTEGER::I,L,J,K,KMAXE,IQP,NGP,iex
INTEGER::TROUBLE
INTEGER::ICONSIDERED,FACEX,POINTX
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV
REAL,DIMENSION(1:NOF_VARIABLES)::MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX
REAL,allocatable,DIMENSION(:,:,:)::USOL
REAL,allocatable,DIMENSION(:,:)::UTEMP
allocate(utemp(IMAXDEGFREE+1,1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR))


allocate(usol(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR,1:6,1:NUMBEROFPOINTS2))

KMAXE=XMPIELRANK(N)


if (code_profile.ne.102)then

!$OMP DO
DO I = 1, KMAXE
ICONSIDERED=I

    
    

    
    CALL FIND_BOUNDS(ICONSIDERED,MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX,UTEMP)

    DO L = 1, IELEM(N,I)%IFCA
            
            if (dimensiona.eq.2)then

            iqp=QP_LINE_N
            else
                if (ielem(n,I)%types_faces(L).eq.5)then
					iqp=qp_quad
				  else
					iqp=QP_TRIANGLE
				  end if
            end if
            
                DO NGP = 1,iqp! 
                    facex=l
                    pointx=ngp
                    USOL(:,facex,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(:,facex,pointx)
                    LEFTV(:)=USOL(:,facex,pointx)
                    CALL PAD_DG(ICONSIDERED,LEFTV)
                    CALL NAD_DG(ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,USOL,MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX,UTEMP)
                END DO
        
    END DO
    
    
    
   
    
    
    
END DO
!$OMP END DO    
    
end if  


deallocate(usol)
DEALLOCATE(UTEMP)


END SUBROUTINE



SUBROUTINE TROUBLE_INDICATOR2
IMPLICIT NONE
INTEGER::I,L,J,K,KMAXE,IQP,NGP,iex,NDOF
INTEGER::TROUBLE, IFREE,I_DEG
INTEGER::ICONSIDERED,FACEX,POINTX

KMAXE=XMPIELRANK(N)

if (code_profile.ne.102)then

!$OMP DO
DO I = 1, KMAXE
ICONSIDERED=I

    

    if (IELEM(N,I)%TROUBLED.eq.1)then
    
    DO L = 1, IELEM(N,I)%IFCA
            
            if (dimensiona.eq.2)then

            iqp=QP_LINE_N
            else
                if (ielem(n,I)%types_faces(L).eq.5)then
					iqp=qp_quad
				  else
					iqp=QP_TRIANGLE
				  end if
            end if
            
                DO NGP = 1,iqp! 
                    facex=l
                    pointx=ngp
                    ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(:, facex,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(:, facex,pointx)
                END DO
        
    END DO
    
      DO IEX=1,NOF_VARIABLES
     
!     
      U_C(ICONSIDERED)%VALDG(1,IEX,2:IDEGFREE+1)=ILOCAL_RECON6(Iconsidered)%DG2FV(1:IDEGFREE,IEX)

      END DO
    
    
    
    
 end if
    
END DO
!$OMP END DO    
    
end if



END SUBROUTINE












    
    
SUBROUTINE PAD_DG(ICONSIDERED,LEFTV)
IMPLICIT NONE
INTEGER::I,L,J,K,KMAXE,IQP,NGP,iex
INTEGER::TROUBLE
INTEGER,INTENT(IN)::ICONSIDERED
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::LEFTV
REAL::MP_PINFL,GAMMAL
    
    I=ICONSIDERED
                                                IF (ITESTCASE.GE.3)THEN
														IF (DIMENSIONA.EQ.3)THEN
                                                
                                                    CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
                                                    
                                                    
                                                    IF(MULTISPECIES.EQ.1)THEN
                                                    
															IF ((LEFTV(1).LE.ZERO).OR.(LEFTV(1).NE.LEFTV(1)))THEN
																IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
															END IF
															IF ((LEFTV(5).LE.ZERO).OR.(LEFTV(5).NE.LEFTV(5)))THEN
																IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
															END IF
	! 														IF ((LEFTV(5).LE.-MP_PINF(1)).OR.(LEFTV(5).LE.-MP_PINF(2)))THEN
	! 															IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
	!                                                         END IF

                                                                
														!	IF((LEFTV(NOF_VARIABLES).LT.-0.05).OR.LEFTV(NOF_VARIABLES).ne.leftv(8)) THEN
														!		IELEM(N,I)%TROUBLED =1; IELEM(N,I)%CONDITION=1
														!	END IF
													ELSE
						
						!
															IF ((LEFTV(1).LE.ZERO).OR.(LEFTV(1).NE.LEFTV(1)))THEN
															IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1


															END IF
															IF ((LEFTV(5).LE.ZERO).OR.(LEFTV(5).NE.LEFTV(5)))THEN
															IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
															END IF
													end if
                                                ELSE
												CALL cons2prim(N,leftv,MP_PINFl,gammal)
                                                IF (MULTISPECIES.EQ.1)THEN

													IF ((LEFTV(1).LE.ZERO).OR.(LEFTV(1).NE.LEFTV(1)))THEN
														IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
													END IF
													IF ((LEFTV(4).LE.ZERO))THEN
														IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
													END IF
													IF ((LEFTV(4).NE.LEFTV(4)))THEN
														IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
													END IF

! IF ((LEFTV(4).LE.ZERO).OR.(LEFTV(4).NE.LEFTV(4)))THEN
! 														IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
! 													END IF

													IF((LEFTV(NOF_VARIABLES).LT.ZERO).OR.LEFTV(NOF_VARIABLES).GT.1.0D0) THEN
														IELEM(N,I)%TROUBLED =1; IELEM(N,I)%CONDITION=1
													END IF





                                                ELSE

                                                        IF ((LEFTV(1).LE.ZERO).OR.(LEFTV(1).NE.LEFTV(1)))THEN						
                                                        IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
                                                        END IF
                                                        IF ((LEFTV(4).LE.ZERO).OR.(LEFTV(4).NE.LEFTV(4)))THEN						
                                                        IELEM(N,I)%TROUBLED=1;IELEM(N,I)%CONDITION=1
                                                        END IF
                                                        
                                                END IF
                                                
                                                END IF
                                                END IF




    
END SUBROUTINE


SUBROUTINE NAD_DG(ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,USOL,MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX,UTEMP)
IMPLICIT NONE
INTEGER::I,L,J,K,KMAXE,IQP,NGP,iex
INTEGER::TROUBLE,img
INTEGER,INTENT(IN)::ICONSIDERED,FACEX,POINTX
REAL::PAR1,PAR2,d2,minb,maxb
REAL,DIMENSION(1:NoF_vARIABLES)::NAD_DG_EL
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::LEFTV,RIGHTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::MAXVARS,AVER_VARS,SUMVARS,UTMIN,UTMAX
REAL,allocatable,DIMENSION(:,:),INTENT(IN)::UTEMP
REAL,allocatable,DIMENSION(:,:,:),INTENT(IN)::USOL
REAL::MP_PINFL,GAMMAL




! PAR1=1E-4
! PAR2=4e-1
    
        if (dimensiona.eq.3)then
    
        img=5
    else
    
        img=4
    end if
        
    
  
       SELECT CASE(INDICATOR_TYPE)
    
            
             CASE(1)       !MOOD INDICATOR
    
            DO IEX=1,NOF_VARIABLES
			NAD_DG_EL(IEX)=MAX(INDICATOR_PAR1,(INDICATOR_PAR2)*(UTMAX(IEX)-UTMIN(IEX)))
			END DO
			leftv(1:nof_Variables)=USOL(1:nof_Variables,facex,pointx)


			IF (DIMENSIONA.EQ.2)THEN
			CALL cons2prim(N,leftv,MP_PINFl,gammal)
			ELSE
			CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
			END IF
    
            DO IEX=1,NOF_VARIABLES
                IF ((leftv(iex).LT.(UTMIN(IEX)-NAD_DG_EL(IEX))).OR.(leftv(iex).GT.(UTMAX(IEX)+NAD_DG_EL(IEX))))THEN
                    IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                END IF
            END DO



             CASE(11)       !MOOD INDICATOR

            DO IEX=1,NOF_VARIABLES
			NAD_DG_EL(IEX)=MAX(INDICATOR_PAR1,(INDICATOR_PAR2)*(UTMAX(IEX)-UTMIN(IEX)))
			END DO


            DO IEX=1,NOF_VARIABLES
                IF ((USOL(iex,facex,pointx).LT.(UTMIN(IEX)-NAD_DG_EL(IEX))).OR.(USOL(iex,facex,pointx).GT.(UTMAX(IEX)+NAD_DG_EL(IEX))))THEN
                    IELEM(N,ICONSIDERED)%CONDx=1
                END IF
            END DO
            
            
            CASE(2)         !SHU INDICATOR

                        DO IEX=1,NOF_VARIABLES
                        
                        
                         IF ((SUMVARS(IEX)/MAXVARS(IEX)).GT.INDICATOR_PAR1)THEN
                            IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                        END IF

                        END DO
                        
             CASE(22)         !Shock detector, INDICATOR (only density and energy)

                        DO IEX=1,NOF_VARIABLES
                        
                        if ((iex.eq.1).or.(iex.eq.img))then
                         IF ((SUMVARS(IEX)/MAXVARS(IEX)).GT.INDICATOR_PAR1)THEN
                            IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                        END IF
                        end if

                        END DO
                        
             CASE(3)         !DMP

                        DO IEX=1,NOF_VARIABLES
                        
                        
                         IF ((USOL(iex,facex,pointx).gt.(UTMax(IEX))).or.(USOL(iex,facex,pointx).lt.(UTmin(IEX))))THEN
                            IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                        END IF

                        END DO            
             
             
             
             CASE(4)    !MINMOD
                        DO IEX=1,NOF_VARIABLES
                        
                       
                        
                         if ((iex.eq.1).or.(iex.eq.img))then
                        IF (ABS(USOL(iex,facex,pointx)-UTEMP(1,IEX)).GT.(INDICATOR_PAR1*UTEMP(1,IEX)))THEN
                        IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                        END IF
                        END IF
                        END DO
             
             
             case(5)    !all troubled
             
             
                        IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                        
                        
              CASE(6)       !MOOD INDICATOR only density & energy
    
            DO IEX=1,NOF_VARIABLES
            if ((iex.eq.1).or.(iex.eq.img))then
			NAD_DG_EL(IEX)=MAX(INDICATOR_PAR1,(INDICATOR_PAR2)*(UTMAX(IEX)-UTMIN(IEX)))
			end if
			END DO
			leftv(1:nof_Variables)=USOL(1:nof_Variables,facex,pointx)

			IF (DIMENSIONA.EQ.2)THEN
			CALL cons2prim(N,leftv,MP_PINFl,gammal)
			eLSE
			CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
			END IF
    
            DO IEX=1,NOF_VARIABLES
            if ((iex.eq.1).or.(iex.eq.img))then

                IF ((leftv(iex).LT.(UTMIN(IEX)-NAD_DG_EL(IEX))).OR.(leftv(iex).GT.(UTMAX(IEX)+NAD_DG_EL(IEX))))THEN
                    IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
                END IF
                end if
            END DO          
                        

			CASE(9)       !MOOD INDICATOR only density & energy

			 IF (IELEM(N,ICONSIDERED)%FILTERED.EQ.1)THEN
				IELEM(N,ICONSIDERED)%TROUBLED=1;IELEM(N,ICONSIDERED)%CONDITION=1
			END IF





                        
                        
                                   
            END SELECT              
                        

END SUBROUTINE NAD_DG




! SUBROUTINE FIND_BOUNDS2
! IMPLICIT NONE
! INTEGER::I,L,J,K,KMAXE,IQP,NGP,IEX,IK,iq
! REAL,DIMENSION(1:NOF_VARIABLES)::AVER_VARS
! REAL,allocatable,DIMENSION(:,:)::UTEMP
! allocate(utemp(IMAXDEGFREE+1,1:NOF_VARIABLES))
!
!
!
! 	KMAXE=XMPIELRANK(N)
!
!
!
!
! 	!$OMP DO
! 	DO I=1,KMAXE
! 	AVER_VARS=ZERO
!
!
!
!
!             UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
!
!             K=1
!             IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
!                 DO L = 1, IELEM(N,I)%IFCA
!                     K=K+1
!                     UTEMP(K,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:NOF_VARIABLES)
!                 END DO
!             END IF
!
!             IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
! 			    DO L=1,IELEM(N,I)%IFCA
!                     IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
!                             IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
!                                 if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
!                                 K=K+1
!                                 UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
!                                 ELSE
!                                 !NOT PERIODIC ONES IN MY CPU
!                                 END IF
!                             ELSE
!                                 K=K+1
!                                 UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
!                             END IF
!                     ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
!
!                         IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
!                             if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
!                             K=K+1
!                             UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
!                             (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
!                             END IF
!                         ELSE
!
!                         K=K+1
!                         UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
!                         (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
!                         END IF
!
!                     END IF
!
! 			  END DO
!          END IF
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
! 			  DO IEX=1,NOF_VARIABLES
!
!                 DO IK=1,K
!                 AVER_VARS(IEX)=AVER_VARS(IEX)+UTEMP(IK,IEX)
!
!                 END DO
!                 AVER_VARS(IEX)=AVER_VARS(IEX)/K
!
!
! 			  END DO
!
!
!
! 			  IELEM(N,I)%AVARS(1:NOF_VARIABLES)=AVER_VARS(1:NOF_VARIABLES)
!
!
! 			  END DO
! 			  !$OMP END DO
!
! 			  deallocate(utemp)
!
!
!
! END SUBROUTINE FIND_BOUNDS2




SUBROUTINE APPLY_FILTER(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,j,k

KMAXE=XMPIELRANK(N)


!$OMP DO
DO I=1,KMAXE
IF (IELEM(N,I)%FILTERED.EQ.1)THEN


DO J=1,NOF_VARIABLES
do k=1,idegfree
RHS(I)%VALDG(k+1,J)=RHS(I)%VALDG(K+1,J)*MODAL_FILTER_WEAK(k)
END DO
end do

END IF
END DO
!$OMP END DO



END SUBROUTINE APPLY_FILTER


SUBROUTINE APPLY_FILTER_DG(n)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,KMAXE
REAL::FILTERED_LOW
REAL::FILTERED_HIGH
REAL::UNFILTERED,ENERGY_RATIO
REAL,DIMENSION(1:NOF_VARIABLES)::EN_F_STRONG,EN_F_WEAK,EN_UNF,EN_AVERAGE
INTEGER::fil_i
real::filx,xorder,EX1,EX2


kmaxe=xmpielrank(n)

!$OMP DO
DO I=1,KMAXE


ielem(n,I)%FILTERED=0


EN_UNF(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)

EN_F_STRONG(1:NOF_VARIABLES)=U_CS(I)%VAL(1,1:NOF_vARIABLES)

EN_F_WEAK(1:NOF_VARIABLES)=U_CW(I)%VAL(1,1:NOF_vARIABLES)

EN_AVERAGE(1:NOF_VARIABLES)=U_C(I)%VALDG(1,1:NOF_VARIABLES,1)



ex2=(((EN_UNF(2)-EN_F_WEAK(2))**2)+((EN_UNF(3)-EN_F_WEAK(3))**2)+((EN_UNF(4)-EN_F_WEAK(4))**2))

ex1=(((EN_UNF(2)-EN_F_STRONG(2))**2)+((EN_UNF(3)-EN_F_STRONG(3))**2)+((EN_UNF(4)-EN_F_STRONG(4))**2))



	ENERGY_RATIO=(EX2+10e-32)/(EX1+10e-32)

    ielem(n,i)%er1dt=((ielem(n,i)%er1-EX1))/dt
	ielem(n,i)%er2dt=((ielem(n,i)%er2-EX2))/dt



  ielem(n,i)%er=ENERGY_RATIO
  ielem(n,i)%er1=EX1
  ielem(n,i)%er2=EX2



 if (ielem(n,i)%er.GT.(1.15))THen

 ielem(n,I)%FILTERED=1
!
 end if
!


END DO
!$OMP END DO


END SUBROUTINE APPLY_FILTER_DG



SUBROUTINE APPLY_FILTER1(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,j,k

KMAXE=XMPIELRANK(N)


!$OMP DO
DO I=1,KMAXE
DO J=1,NOF_VARIABLES
do k=1,idegfree
U_C(I)%VALDG(1,J,K+1)=U_C(I)%VALDG(1,J,K+1)*MODAL_FILTER(k)
END DO
end do
END DO
!$OMP END DO



END SUBROUTINE APPLY_FILTER1


SUBROUTINE APPLY_FILTER2(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,j,k
REAL,ALLOCATABLE,DIMENSION(:)::GRAD1AL

KMAXE=XMPIELRANK(N)
ALLOCATE(GRAD1AL(1:IDEGFREE))

!$OMP DO
DO I=1,KMAXE
DO J=1,NOF_VARIABLES
GRAD1AL(1:IDEGFREE)=0.0d0
do k=1,idegfree
GRAD1AL(k)=RHS(I)%VALDG(k+1,j)*MODAL_FILTER(k)
end do
RHS(I)%VALDG(2:IDEGFREE+1,j)=GRAD1AL(1:idegfree)
END DO
END DO
!$OMP END DO
DEALLOCATE(GRAD1AL)


END SUBROUTINE APPLY_FILTER2



SUBROUTINE FILTER(N)
implicit none
INTEGER,INTENT(IN)::N
INTEGER::fil_i,i,j
real::filx,xorder
real::rfil_alpha,rfil_nc,rfil_s,rfil_i
real,dimension(1:9)::filter2
real,dimension(1:200)::dgfr
integer,dimension(0:9)::filt2






IF (FILTER_TYPE.EQ.1)THEN

do fil_i=1,iorder
      if (fil_i.le.fil_nc)then
        filter2(fil_i)=1.0d0
    end if
     if (fil_i.gt.fil_nc)then
		 rfil_alpha=fil_alpha
		 rfil_nc=fil_nc
		 rfil_s=fil_s
		 xorder=iorder
		 rfil_i=fil_i
         filx=-rfil_alpha*(((rfil_i-rfil_nc)/(xorder-rfil_nc))**rfil_s)
         filter2(fil_i)=exp(filx)
     end if


	filt2(fil_i)=(((fil_i+1)*(fil_i+2)*(fil_i+3))/6)-1
end do
filt2(0)=0

do i=1,iorder
	do j=filt2(i-1)+1,filt2(i)
	dgfr(j)=filter2(i)
	end do
end do

do fil_i=1,IDEGFREE
        MODAL_FILTER(fil_i)=dgfr(fil_i)**(1.0d0/(1.0/dt))!
        IF ((IT.le.2).and.(n.eq.0))THEN
        WRITE(200+N,*)FIL_I,MODAL_FILTER(fil_i)
        END IF

end do




END IF

IF (FILTER_TYPE.EQ.2)THEN

do fil_i=1,iorder
      if (fil_i.lt.iorder)then
        filter2(fil_i)=1.0d0
    end if
     if (fil_i.eq.iorder)then
         filter2(fil_i)=0.0d0
     end if


	filt2(fil_i)=(((fil_i+1)*(fil_i+2)*(fil_i+3))/6)-1
end do
filt2(0)=0

do i=1,iorder
	do j=filt2(i-1)+1,filt2(i)
	dgfr(j)=filter2(i)
	end do
end do

do fil_i=1,IDEGFREE
        MODAL_FILTER(fil_i)=dgfr(fil_i)**(1.0d0/(1.0/dt))!
        IF ((IT.le.2).and.(n.eq.0))THEN
        WRITE(200+N,*)FIL_I,MODAL_FILTER(fil_i)
        END IF

end do




END IF



IF (FILTER_TYPE.EQ.3)THEN

do fil_i=1,iorder
     if (fil_i.lt.iorder)then
          filter2(fil_i)=1.0d0
     eLSE
 		filter2(fil_i)=0.0d0
     END IF

!  	if (fil_i.le.fil_nc)then
!          filter2(fil_i)=1.0d0
!      end if
!       if (fil_i.gt.fil_nc)then
!  		 rfil_alpha=fil_alpha
!  		 rfil_nc=fil_nc
!  		 rfil_s=fil_s
!  		 xorder=iorder
!  		 rfil_i=fil_i
!           filx=-rfil_alpha*(((rfil_i-rfil_nc)/(xorder-rfil_nc))**rfil_s)
!           filter2(fil_i)=exp(filx)
!       end if




	filt2(fil_i)=(((fil_i+1)*(fil_i+2)*(fil_i+3))/6)-1
end do
filt2(0)=0

do i=1,iorder
	do j=filt2(i-1)+1,filt2(i)
	dgfr(j)=filter2(i)
	end do
end do



do fil_i=1,IDEGFREE
        MODAL_FILTER_weak(fil_i)=dgfr(fil_i)!**(1.0d0/(1.0/dt))!

end do

do fil_i=1,iorder
     if (fil_i.lt.iorder-1)then
          filter2(fil_i)=1.0d0
     eLSE
 		filter2(fil_i)=0.0d0
     END IF


	filt2(fil_i)=(((fil_i+1)*(fil_i+2)*(fil_i+3))/6)-1
end do
filt2(0)=0

do i=1,iorder
	do j=filt2(i-1)+1,filt2(i)
	dgfr(j)=filter2(i)
	end do
end do



do fil_i=1,IDEGFREE
        MODAL_FILTER_strong(fil_i)=dgfr(fil_i)!**(1.0d0/(1.0/dt))!
end do




do fil_i=1,IDEGFREE
        IF ((IT.le.2).and.(n.eq.0))THEN
        WRITE(200+N,*)FIL_I,MODAL_FILTER_WEAK(fil_i),MODAL_FILTER_STRONG(fil_i)
        END IF

end do





END IF









END SUBROUTINE FILTER




SUBROUTINE ADDA_FILTER(N,iconsidered, IELEM_L, ILOCAL_RECON3_L, ILOCAL_RECON5_L, U_C_L, INTEG_BASIS_L, integ_basis_dg_L)
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
INTEGER,INTENT(IN)::N,iconsidered

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3_L, ILOCAL_RECON5_L
TYPE(U_CENTRE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::U_C_L
TYPE(INTEGRALBASIS),ALLOCATABLE,DIMENSION(:)::INTEG_BASIS_L,integ_basis_dg_L

INTEGER::I,J,K
REAL::FILTERED_LOW
REAL::FILTERED_HIGH
REAL::UNFILTERED,ENERGY_RATIO
REAL,DIMENSION(1:NOF_VARIABLES)::EN_F_STRONG,EN_F_WEAK,EN_UNF
INTEGER::fil_i,countdof,icompwrt,ngp
real::filx,xorder,EX1,EX2
real::rfil_alpha,rfil_nc,rfil_s,rfil_i
real,dimension(1:9)::filter2,filter3
real,dimension(1:((IORDER+1)*(IORDER+2)*(IORDER+3))/6)::dgfr,dgfr3
integer,dimension(0:9)::filt2,filt3
REAL::AX,AY,AZ,MP_PINFl,gammal
REAL,dimension(1:NOF_VARIABLES)::LEFTV
REAL,DIMENSION(1,1:IDEGFREE)::CONSMATRIX
REAL,DIMENSION(1:IDEGFREE,1:NOF_VARIABLES)::GRADSSL
REAL,DIMENSION(1:6*NUMBEROFPOINTS2,1:NOF_VARIABLES)::RESSOLUTION



countdof=((IORDER+1)*(IORDER+2)*(IORDER+3))/6




do fil_i=1,iorder
	 if (iorder.eq.2) then
		filter2(fil_i)=0.0d0
	else

     if (fil_i.lt.2)then  !if (fil_i.le.ADDA_1)then
        filter2(fil_i)=1.0d0
    else
		filter2(fil_i)=0.0d0


    end if
    END IF



	filt2(fil_i)=(((fil_i+1)*(fil_i+2)*(fil_i+3))/6)-1
end do
filt2(0)=0

do i=1,iorder
	do j=filt2(i-1)+1,filt2(i)
	dgfr(j)=filter2(i)
	end do
end do


do fil_i=1,iorder
	if (fil_i.le.2)then!					if (fil_i.le.ADDA_2)then
        filter3(fil_i)=1.0d0
    else
		filter3(fil_i)=0.0d0
    end if



	filt3(fil_i)=(((fil_i+1)*(fil_i+2)*(fil_i+3))/6)-1
end do
filt3(0)=0

do i=1,iorder
	do j=filt3(i-1)+1,filt3(i)
	dgfr3(j)=filter3(i)
	end do
end do


do fil_i=1,IDEGFREE
        ADDA_FILTER_STRONG(fil_i)=dgfr(fil_i)
        ADDA_FILTER_WEAK(fil_i)=dgfr3(fil_i)
        if ((it.eq.0).and.(iconsidered.eq.1).AND.(N.EQ.0))then
        write(300+n,*)fil_i,ADDA_FILTER_WEAK(fil_i),ADDA_FILTER_STRONG(fil_i)

        end if
end do


I=ICONSIDERED

IF (DG.NE.1)THEN



if (adda_type.eq.1)then

AX = 0.0D0;AY = 0.0D0;AZ = 0.0D0

icompwrt=0


				CONSMATRIX(1,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)


				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,1:NOF_VARIABLES)


                RESSOLUTION(1:1,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:1,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))



				LEFTV(1:NOF_VARIABLES)=U_C_L(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(1,1:NOF_vARIABLES)


                EN_UNF(1:NOF_VARIABLES)=LEFTV(1:NOF_vARIABLES)




                do k=1,nof_Variables
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,k)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,k)*ADDA_FILTER_strong(1:IELEM_L(N,I)%IDEGFREE)
				end do



                RESSOLUTION(1:1,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:1,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))


                LEFTV(1:NOF_VARIABLES)=U_C_L(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(1,1:NOF_vARIABLES)

                EN_F_STRONG(1:NOF_VARIABLES)=LEFTV(1:NOF_vARIABLES)



				do k=1,nof_Variables
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,k)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,k)*ADDA_FILTER_weak(1:IELEM_L(N,I)%IDEGFREE)
				end do



                RESSOLUTION(1:1,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:1,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))

                LEFTV(1:NOF_VARIABLES)=U_C_L(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(1,1:NOF_vARIABLES)

                EN_F_WEAK(1:NOF_VARIABLES)=LEFTV(1:NOF_vARIABLES)


                ex2=(((EN_UNF(2)-EN_F_WEAK(2))**2)+((EN_UNF(3)-EN_F_WEAK(3))**2)+((EN_UNF(4)-EN_F_WEAK(4))**2))
                ex1=(((EN_UNF(2)-EN_F_STRONG(2))**2)+((EN_UNF(3)-EN_F_STRONG(3))**2)+((EN_UNF(4)-EN_F_STRONG(4))**2))



                ENERGY_RATIO=(EX2+10E-32)/(EX1+10E-32)

end if

if (adda_type.eq.2)then
EN_UNF(1:NOF_VARIABLES)=zero
EN_F_STRONG(1:NOF_VARIABLES)=zero
EN_F_WEAK(1:NOF_VARIABLES)=zero

do ngp=1,IELEM_L(N,ICONSIDERED)%ITOTALPOINTS

ax=QP_ARRAY(ICONSIDERED)%X(ngp)
ay=QP_ARRAY(ICONSIDERED)%Y(ngp)
az=QP_ARRAY(ICONSIDERED)%Z(ngp)



icompwrt=0

CONSMATRIX(1,1:IELEM_L(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM_L(N,I)%IORDER,I,IELEM_L(N,I)%IDEGFREE,icompwrt, IELEM_L, ILOCAL_RECON3_L, INTEG_BASIS_L, integ_basis_dg_L)
GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,1:NOF_VARIABLES)
RESSOLUTION(1:1,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:1,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))


LEFTV(1:NOF_VARIABLES)=U_C_L(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(1,1:NOF_vARIABLES)


EN_UNF(1:NOF_VARIABLES)=EN_UNF(1:NOF_VARIABLES)+LEFTV(1:NOF_vARIABLES)*QP_ARRAY(ICONSIDERED)%QP_WEIGHT(ngp)




                do k=1,nof_Variables
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,k)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,k)*ADDA_FILTER_strong(1:IELEM_L(N,I)%IDEGFREE)
				end do



                RESSOLUTION(1:1,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:1,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))


                LEFTV(1:NOF_VARIABLES)=U_C_L(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(1,1:NOF_vARIABLES)

!                 EN_F_STRONG(1:NOF_VARIABLES)=LEFTV(1:NOF_vARIABLES)

                EN_F_STRONG(1:NOF_VARIABLES)=EN_F_STRONG(1:NOF_VARIABLES)+LEFTV(1:NOF_vARIABLES)*QP_ARRAY(ICONSIDERED)%QP_WEIGHT(ngp)



				do k=1,nof_Variables
				GRADSSL(1:IELEM_L(N,I)%IDEGFREE,k)=ILOCAL_RECON5_L(ICONSIDERED)%GRADIENTS(1,1:IELEM_L(N,I)%IDEGFREE,k)*ADDA_FILTER_weak(1:IELEM_L(N,I)%IDEGFREE)
				end do



                RESSOLUTION(1:1,1:NOF_vARIABLES)=matmul(CONSMATRIX(1:1,1:IELEM_L(N,I)%IDEGFREE),GRADSSL(1:IELEM_L(N,I)%IDEGFREE,1:NOF_vARIABLES))

                LEFTV(1:NOF_VARIABLES)=U_C_L(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(1,1:NOF_vARIABLES)

!                 EN_F_WEAK(1:NOF_VARIABLES)=LEFTV(1:NOF_vARIABLES)

                EN_F_WEAK(1:NOF_VARIABLES)=EN_F_WEAK(1:NOF_VARIABLES)+LEFTV(1:NOF_vARIABLES)*QP_ARRAY(ICONSIDERED)%QP_WEIGHT(ngp)


end do


                ex2=(((EN_UNF(2)-EN_F_WEAK(2))**2)+((EN_UNF(3)-EN_F_WEAK(3))**2)+((EN_UNF(4)-EN_F_WEAK(4))**2))
                ex1=(((EN_UNF(2)-EN_F_STRONG(2))**2)+((EN_UNF(3)-EN_F_STRONG(3))**2)+((EN_UNF(4)-EN_F_STRONG(4))**2))



                ENERGY_RATIO=(EX2+10E-32)/(EX1+10E-32)


end if















				!ENERGY_RATIO=ABS(ENERGY_RATIO-1.0D0)/(IELEM_L(N,I)%TOTVOLUME)








	IELEM_L(n,i)%er1dt=(IELEM_L(n,i)%er1-EX1)/dt
	IELEM_L(n,i)%er2dt=(IELEM_L(n,i)%er2-EX2)/dt

	!IELEM_L(n,i)%er=(IELEM_L(n,i)%erX-ENERGY_RATIO)/dt

  IELEM_L(n,i)%er=ENERGY_RATIO
  IELEM_L(n,i)%er1=EX1
  IELEM_L(n,i)%er2=EX2

! 	IELEM_L(n,i)%er1er2=0.0d0

! 	if (IELEM_L(n,i)%er2dt.gt.0)then
	IELEM_L(n,i)%er1er2=abs(IELEM_L(n,i)%er1dt)/IELEM_L(n,i)%er2dt
! 	end if



   CALL APPLY_ADDA_FILTER(N,iconsidered, IELEM_L)



END IF

END SUBROUTINE ADDA_FILTER



SUBROUTINE APPLY_ADDA_FILTER(N,iconsidered, IELEM_L)
IMPLICIT NONE
#ifdef WENOWEIGHTS_GPU_KERNEL
!$omp declare target
#endif
INTEGER,INTENT(IN)::N
INTEGER::I,ICONSIDERED
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM_L
real::LWCX1

I=ICONSIDERED

LWCX1=LWCI1

IF (RUNGEKUTTA.EQ.11)THEN
IF (ISCOUN.EQ.1)THEN

LWCX1=LWCI1





		 IF (IELEM_L(n,i)%er.gt.1.2)THEN

				LWCX1=10!INCREASE DISSIPATION


		 	end if


		 IF (IELEM_L(n,i)%er.LE.0.95)THEN

						LWCX1=1000
			
			
		       
		 END IF





		IELEM_L(n,i)%lwcx2=lwcx1


ELSE

LWCX1=IELEM_L(n,i)%lwcx2


END IF



ELSE


		 IF (IELEM_L(n,i)%er.gt.1.2)THEN

				LWCX1=10!INCREASE DISSIPATION


		 	end if

! END IF
!
		 IF (IELEM_L(n,i)%er.LE.0.95)THEN

		 	

			 LWCX1=1000
			 !LWCX1=100**(12-(4*IELEM_L(n,i)%er**0.8))


		 END IF





		IELEM_L(n,i)%lwcx2=lwcx1



END IF


if (IELEM_L(n,i)%full.eq.0)then
        IELEM_L(n,i)%lwcx2=-10
end if



IELEM_L(n,i)%LINC=LWCX1


END SUBROUTINE APPLY_ADDA_FILTER






SUBROUTINE FIX_DISSIPATION(N)
IMPLICIT NONE
REAL::CHECK1
REAL::CHECK
INTEGER::I,J,K,L,KMAXE
integer,intent(in)::n

KMAXE=XMPIELRANK(N)

!$OMP DO
DO I=1,KMAXE

	if (ielem(n,i)%full.eq.1)then
	!1)reduce dissipation
	IF (IELEM(N,I)%LWCX2.GT.10)THEN
		IF (IELEM(N,I)%WCX(1).GE.0.999)THEN
		IELEM(N,I)%DISS=max(IELEM(N,I)%DISS-0.1,0.5d0)	!reduce dissipation even more
		Else
		IELEM(N,I)%DISS=1.0d0							!increase dissipation if shock
		END IF
	ELSE
	!2) INCREASE DISSIPATION
		IF (IELEM(N,I)%WCX(1).GE.0.999)THEN
		IELEM(N,I)%DISS=min(IELEM(N,I)%DISS+0.1,1.0d0)	!increase dissipation even more
		Else
		IELEM(N,I)%DISS=1.0d0							!increase dissipation if shock
		END IF
	end if
	end if

END DO
!$OMP END DO


END SUBROUTINE FIX_DISSIPATION




SUBROUTINE FIX_DISSIPATION2(N)
IMPLICIT NONE
REAL::CHECK1
REAL::CHECK
INTEGER::I,J,K,L,KMAXE,ICONSIDERED
integer,intent(in)::n

KMAXE=XMPIELRANK(N)

!$OMP DO
DO I=1,KMAXE
	iconsidered=i

	call FIND_BOUNDS_DISS(iconsidered)

END DO
!$OMP END DO


END SUBROUTINE FIX_DISSIPATION2






SUBROUTINE FIND_BOUNDS_DISS(iconsidered)
IMPLICIT NONE
INTEGER::I,L,J,K,KMAXE,IQP,NGP,IEX,IK,iq
integer,intent(in)::iconsidered
REAL,DIMENSION(1:6,1)::UTEMP

I=ICONSIDERED


UTEMP(:,:)=1.0D0



		IELEM(N,I)%FACEDISS(:)=1.0d0


            IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
                DO L = 1, IELEM(N,I)%IFCA

                    UTEMP(L,1)=IELEM(N,IELEM(N,I)%INEIGH(L))%DISS
                END DO
            END IF

            IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
			    DO L=1,IELEM(N,I)%IFCA
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                            IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU

                                UTEMP(l,1)=IELEM(N,IELEM(N,I)%INEIGH(L))%DISS
                                ELSE
                                !NOT PERIODIC ONES IN MY CPU
                                END IF
                            ELSE

                                UTEMP(l,1)=IELEM(N,IELEM(N,I)%INEIGH(L))%DISS
                            END IF
                    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS

                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU

                            UTEMP(l,1)=IEXSOLHIRd(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                            (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1)
                            END IF
                        ELSE


                        UTEMP(l,1)=IEXSOLHIRd(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                        (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1)
                        END IF

                    END IF

			  END DO
         END IF





            IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
                DO L = 1, IELEM(N,I)%IFCA
                    IELEM(N,I)%FACEDISS(L)=MAX(UTEMP(L,1),IELEM(N,I)%DISS)

                END DO
            END IF

            IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
			    DO L=1,IELEM(N,I)%IFCA
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                            IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                IELEM(N,I)%FACEDISS(L)=MAX(UTEMP(L,1),IELEM(N,I)%DISS)
!                                 UTEMP(l,1)=IELEM(N,IELEM(N,I)%INEIGH(L))%DISS
                                ELSE
                                !NOT PERIODIC ONES IN MY CPU
                                END IF
                            ELSE
                                IELEM(N,I)%FACEDISS(L)=MAX(UTEMP(L,1),IELEM(N,I)%DISS)

                            END IF
                    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS

                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU

                            IELEM(N,I)%FACEDISS(L)=MAX(UTEMP(L,1),IELEM(N,I)%DISS)
                            END IF
                        ELSE


                        IELEM(N,I)%FACEDISS(L)=MAX(UTEMP(L,1),IELEM(N,I)%DISS)
                        END IF

                    END IF

			  END DO
         END IF








END SUBROUTINE FIND_BOUNDS_DISS



SUBROUTINE VFBP_LIMITER

IMPLICIT NONE
REAL::VF_QPSOL, VF_AVSOL, LTHRESH, HTHRESH,SCALING,SCALING1,SCALING2, PD1_QPSOL, PD1_AVSOL,PD2_QPSOL, PD2_AVSOL
INTEGER::I,II,ICD,K,NUMBER_OF_NEI,IDUMMy,L,NND,NGP, IQP,I_ELEM,I_FACE, IMG1,IMG2
REAL, DIMENSION(1:NUMBEROFPOINTS2)::SCALING_F,SCALING_F1,SCALING_F2
INTEGER::FACEX,POINTX,ICONSIDERED,NUMBER_OF_DOG

 !iqp=QP_LINE_N

 LTHRESH = 1.0e-16
 HTHRESH = 1.0d0-LTHRESH



!$OMP DO
DO I = 1, XMPIELRANK(N)

    DO L = 1, IELEM(N,I)%IFCA



    if (dimensiona.eq.2)then

            iqp=QP_LINE_N
            IMG1=5
            IMG2=6
            else
            IMG1=6
            IMG2=7
                if (ielem(n,I)%types_faces(L).eq.5)then
					iqp=qp_quad
				  else
					iqp=QP_TRIANGLE

				  end if
            end if

 !COMPUTE SCALING FACTORS AT EACH FACE AND EACH QP

     DO NGP = 1,iqp! QP_LINE_N

                FACEX=L
                POINTX=NGP
                ICONSIDERED=I
                NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE

                VF_QPSOL = ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(NOF_VARIABLES, FACEX, POINTX)
                VF_AVSOL = U_C(I)%VALDG(1,NOF_VARIABLES,1)

                PD1_QPSOL = ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(IMG1, FACEX, POINTX)
                PD1_AVSOL = U_C(I)%VALDG(1,IMG1,1)

                PD2_QPSOL = ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(IMG2, FACEX, POINTX)
                PD2_AVSOL = U_C(I)%VALDG(1,IMG2,1)



        !VOLUME FRACTION
                IF (VF_QPSOL.LT.LTHRESH) THEN

                    SCALING_F(NGP) = (VF_AVSOL - LTHRESH)/(VF_AVSOL - VF_QPSOL)

                ELSE IF ((VF_QPSOL.GT.LTHRESH).and.(VF_QPSOL.lt.HTHRESH)) THEN

                    SCALING_F(NGP) = 1.0d0

                ELSE IF (VF_QPSOL.GT.HTHRESH) THEN

                    SCALING_F(NGP) = (HTHRESH - VF_AVSOL)/(VF_QPSOL - VF_AVSOL)

                END IF




        !PARTIAL DENSITIES
                IF (PD1_QPSOL.LT.LTHRESH) THEN

                    SCALING_F1(NGP) = (PD1_AVSOL - LTHRESH)/(PD1_AVSOL - PD1_QPSOL)
                ELSE

                    SCALING_F1(NGP) = 1.0D0

                END IF



                IF (PD2_QPSOL.LT.LTHRESH) THEN

                    SCALING_F2(NGP) = (PD2_AVSOL - LTHRESH)/(PD2_AVSOL - PD2_QPSOL)

                ELSE

                    SCALING_F2(NGP) = 1.0D0

                END IF



            END DO



            SCALING = MINVAL(SCALING_F(:))
            SCALING1 = MINVAL(SCALING_F1(:))
            SCALING2 = MINVAL(SCALING_F2(:))





            DO NGP = 1,iqp! QP_LINE_N

                FACEX=L
                POINTX=NGP
                ICONSIDERED=I
                NUMBER_OF_DOG=IELEM(N,I)%IDEGFREE

                VF_QPSOL = ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(NOF_VARIABLES, FACEX, POINTX)
                VF_AVSOL = U_C(I)%VALDG(1,NOF_VARIABLES,1)

                PD1_QPSOL = ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(IMG1, FACEX, POINTX)
                PD1_AVSOL = U_C(I)%VALDG(1,IMG1,1)

                PD2_QPSOL = ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(IMG2, FACEX, POINTX)
                PD2_AVSOL = U_C(I)%VALDG(1,IMG2,1)


            !VOLUME FRACTION SCALING AT QPS
                IF((SCALING.GT.0.0d0).and.(SCALING.LT.1.0d0)) THEN
                    ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(NOF_VARIABLES, FACEX, POINTX) = VF_AVSOL + SCALING*(VF_QPSOL - VF_AVSOL)

                END IF

            !PARTIAL DENSITIES SCALING AT QPS

                IF((SCALING1.GT.0.0d0).and.(SCALING1.LT.1.0d0)) THEN
                    ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(IMG1, FACEX, POINTX) = PD1_AVSOL + SCALING1*(PD1_QPSOL - PD1_AVSOL)

                END IF

                IF((SCALING2.GT.0.0d0).and.(SCALING2.LT.1.0d0)) THEN
                    ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(IMG2, FACEX, POINTX) = PD2_AVSOL + SCALING2*(PD2_QPSOL - PD2_AVSOL)

                END IF




            END DO





    END DO

END DO
!$OMP END DO


END SUBROUTINE VFBP_LIMITER


! ! ---------------------------------------------------------------------------------------------!

! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE RECON
