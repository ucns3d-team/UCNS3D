MODULE RECON
USE BLAS95
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




subroutine check_fs
!> @brief
!> For debugging purposes
implicit none
integer::i,kmaxe
kmaxe=xmpielrank(n)

end subroutine check_fs



SUBROUTINE AVERAGE_STRESSES(N)
!> @brief
!> Subroutine for calling the computation of the average shear stresses
INTEGER,INTENT(IN)::N
INTEGER::II,I
!$OMP DO SCHEDULE (STATIC)
DO II=1,NOF_INTERIOR;I=EL_INT(II);ICONSIDERED=I
      CALL ALLGRADS_INNER_AV(N,I)
END DO
!$OMP END DO 

!$OMP DO SCHEDULE (STATIC) 
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
KMAXE=XMPIELRANK(N)

       if (dimensiona.eq.3)then
	DO I=1,KMAXE
		IF (IELEM(N,I)%ISHAPE.EQ.2)THEN
		ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
		ELSE
		ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
		END IF
		ICONSIDERED=I
		DO L=1,IELEM(N,I)%IFCA
		IDUMMY=0
			    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				    IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					    if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
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
						      END DO
					    ELSE
						      facex=l;
						      CALL coordinates_face_PERIOD(n,iconsidered,facex)
						      do K=1,nnd
						      VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						      END DO
					    END IF
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    else
				    iqp=QP_TRIANGLE
				    NND=3
					    IF (IDUMMY.EQ.0)THEN
						      do K=1,nnd
							VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
							VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						      END DO
					    ELSE
						      facex=l;
						      CALL coordinates_face_PERIOD(n,iconsidered,facex)
						      do K=1,nnd
						      VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						      END DO
					    END IF
						      
				    call QUADRATURETRIANG(N,IGQRULES)
				    end if
			    else
					  if (ielem(n,i)%types_faces(L).eq.5)then
					    iqp=qp_quad
					    NND=4
						  do K=1,nnd
						    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
						    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						  END DO 
					    call  QUADRATUREQUAD3D(N,IGQRULES)
					  else
					    iqp=QP_TRIANGLE
					    NND=3 
						  do K=1,nnd
						    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
						    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
						  END DO  
					    call QUADRATURETRIANG(N,IGQRULES)
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
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2D(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  END IF
				    call  QUADRATUREline(N,IGQRULES)	  
		else
				   iqp=qp_line
				    NND=2
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO  
				    call  QUADRATUREline(N,IGQRULES)

		end if
			     do NGP=1,iqp			!for gqp
				ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:2)=QPOINTS2D(1:2,NGP)
				END DO	!NGP
	  END DO
	  END DO
	  end if
END SUBROUTINE MEMORY_FAST

! ! !---------------------------------------------------------------------------------------------!
SUBROUTINE EXTRAPOLATE_BOUND(varcons,FACEX,pointx,ICONSIDERED,INSTEN,LLX)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces
IMPLICIT NONE
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED,INSTEN,LLX
				    IF (DIMENSIONA.EQ.2)THEN
				     if (reduce_comp.eq.1)then
				     ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,1)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,1)&
				    +((U_C(ICONSIDERED)%VAL(1,varcons)+RESSOLUTION(1,1))*WENO(varcons,INSTEN))*WEQUA2D(pointx)
				    
				    
				    else
				     ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,pointx)&
				    +((U_C(ICONSIDERED)%VAL(1,varcons)+RESSOLUTION(1,1))*WENO(varcons,INSTEN))
				    
				    
				    end if
				    ELSE
				    if (reduce_comp.eq.1)then
				    
				    ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,1)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,1)&
				    +((U_C(ICONSIDERED)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(INSTEN,1:NOF_vARIABLES))*WENO(1:NOF_vARIABLES,LLX))*WEIGHT_T2(pointx)
				    
				    else
				     ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)&
				     +(U_C(ICONSIDERED)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(INSTEN,1:NOF_vARIABLES))*WENO(1:NOF_vARIABLES,LLX)
				    
				    end if
				    
				    
				    
				    END IF



END SUBROUTINE EXTRAPOLATE_BOUND


SUBROUTINE EXTRAPOLATE_BOUND_LINEAR(varcons,FACEX,pointx,ICONSIDERED,INSTEN)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces for linear advection equation
IMPLICIT NONE
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED,INSTEN

	IF (DIMENSIONA.EQ.2)THEN 
	if (reduce_comp.eq.1)then
	  ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,1)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,1)&
	+(U_C(ICONSIDERED)%VAL(1,varcons)+RESSOLUTION(1,1))*WEQUA2D(pointx)
	  else
	ILOCAL_RECON3(ICONSIDERED)%ULEFT(varcons,FACEX,pointx)=(U_C(ICONSIDERED)%VAL(1,varcons)+RESSOLUTION(1,1))
	  end if
	ELSE
	if (reduce_comp.eq.1)then
	ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,1)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,1)&
	+(U_C(ICONSIDERED)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(INSTEN,1:NOF_VARIABLES))*WEiGHT_t2(pointx) 
	else
	ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:NOF_VARIABLES,FACEX,pointx)&
	+(U_C(ICONSIDERED)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(INSTEN,1))
	end if
	END IF
END SUBROUTINE EXTRAPOLATE_BOUND_LINEAR


SUBROUTINE EXTRAPOLATE_BOUND_MUSCL(varcons,FACEX,pointx,ICONSIDERED,INSTEN)
!> @brief
!> Subroutine for extrapolating the reconstructed solution at the cell interfaces for MUSCL schemes
IMPLICIT NONE
INTEGER,INTENT(IN)::varcons,FACEX,pointx,ICONSIDERED,INSTEN
IF (DIMENSIONA.EQ.2)THEN
if (reduce_comp.eq.1)then
ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,1)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,1)&
+(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)*WENO(1:nof_Variables,1)))*WEQUA2D(pointx)
ELSE
ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)&
+(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)*WENO(1:nof_Variables,1)))
END IF
ELSE
if (reduce_comp.eq.1)then
ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,1)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,1)&
+(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)*WENO(1:nof_Variables,1)))*WEIGHT_T2(pointx)
ELSE
ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)=ILOCAL_RECON3(ICONSIDERED)%ULEFT(1:nof_Variables,FACEX,pointx)&
+(U_C(ICONSIDERED)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,FACEX,pointx)*WENO(1:nof_Variables,1)))
END IF


END IF


END SUBROUTINE 


subroutine diag_At_B_A(ICONSIDERED,A_CHAR,B_CHAR,x_CHAR)
!> @brief
!> Subroutine For general matrix A and square matrix B, computes the vector x = diag(A' * B * A), used for characteristics reconstruction
implicit none
   real, dimension(:,:,:), intent(in):: A_char
   real, dimension(:,:), intent(inOUT)::  B_char
   real, dimension(:,:), intent(inout):: x_char
   integer, intent(in):: ICONSIDERED
   real, dimension(:,:,:), allocatable:: BA_char

   integer:: nn, mm
   integer:: i,LL,ICS

   nn = size(A_char,1) ! = size(B,1) = size(B,2)
   mm = size(A_char,2)
   
   

   allocate(BA_char(nn, mm,10))
   
   x_char=zero
   
   IF (EES.EQ.5)THEN
	  DO LL=1,1
	  call gemm(B_char(:,:), A_char(:,:,LL), BA_char(:,:,LL)) ! BA = B * A
	  END DO
	B_CHAR(1:IDEGFREE2,1:IDEGFREE2)=ILOCAL_RECON3(ICONSIDERED)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2)
	  DO LL=2,IELEM(N,ICONSIDERED)%ADMIS
	  call gemm(B_char(1:IDEGFREE2,1:IDEGFREE2), A_char(1:IDEGFREE2,1:nof_variables,LL), BA_CHAR(1:IDEGFREE2,1:nof_Variables,LL)) ! BA = B * A
	  END DO
	  DO LL=1,1;do i = 1, mm
	    x_char(i,LL) = dot(A_char(:,i,LL), BA_char(:,i,LL))
	  end do;END DO
	  DO LL=2,IELEM(N,ICONSIDERED)%ADMIS;do i = 1, mm
	    x_char(i,LL) = dot(A_char(1:IDEGFREE2,i,LL), BA_char(1:IDEGFREE2,i,LL))
	  end do;END DO
   ELSE
   
	      DO LL=1,IELEM(N,ICONSIDERED)%ADMIS
	      call gemm(B_char(:,:), A_char(:,:,LL), BA_char(:,:,LL)) ! BA = B * A
	      END DO
	      DO LL=1,IELEM(N,ICONSIDERED)%ADMIS;do i = 1, mm
		  x_char(i,LL) = dot(A_char(:,i,LL), BA_char(:,i,LL))
	      end do;END DO
   END IF
   
   deallocate(BA_CHAR)
end subroutine diag_At_B_A

subroutine compute_gradcharv_smoothindicator(ICONSIDERED, FACEX)
!> @brief
!> Subroutine for characteristics reconstruction of WENO schemes
IMPLICIT NONE
   integer, intent(inOUT):: ICONSIDERED, FACEX
   integer:: LL, k,I,L
   real, dimension(0:idegfree,nof_variables,typesten):: gradients
   real, dimension(0:idegfree,nof_variables,typesten):: gradients_eigvlt
   I=ICONSIDERED
   L=FACEX
   GRADCHARV=zero
   gradients(:,:,:)=ZERO
   DO LL=1,IELEM(N,I)%ADMIS
    gradients(0,:,ll) = U_C(I)%VAL(1,1:nof_variables)
   end do
   IF (EES.EQ.5)THEN
   ;LAMC(:)=ZERO;GRAD5ALC=ZERO;LAMC(1)=(1.0d0-(1.0d0/LWCI1));lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
		    DO LL=2,IELEM(N,I)%ADMIS
			GRAD5ALC(1:IDEGFREE2,1:nof_variables)=GRAD5ALC(1:IDEGFREE2,1:nof_variables)&
			+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC(ll,1:IDEGFREE2,1:nof_variables))
			gradients(1:IDEGFREE2,1:nof_variables,ll)=ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,1:nof_variables)
		    END DO
		    DO LL=1,1
		      gradients(1:IELEM(N,I)%idegfree,1:nof_variables,ll)=(1.0D0/LAMC(1))*&
		      (ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1:nof_variables)-GRAD5ALC(1:IELEM(N,I)%IDEGFREE,1:nof_variables))
		    END DO
		    DO LL=1,1
			    CALL GEMM(                                            &
			    gradients(0:IELEM(N,I)%idegfree,1:nof_variables,ll),                    &
			    EIGVL(1:nof_variables,1:nof_variables),                           &
			    gradients_eigvlt(0:IELEM(N,I)%idegfree,1:nof_variables,ll),             &
			    'N',                                              & ! transposition flag for gradients
			    'T'                                               & ! transposition flag for EIGVL
			)
		      END DO
		    DO LL=2,IELEM(N,I)%ADMIS
			    CALL GEMM(                                            &
			    gradients(0:IDEGFREE2,1:nof_variables,ll),                    &
			    EIGVL(1:nof_variables,1:nof_variables),                           &
			    gradients_eigvlt(0:IDEGFREE2,1:nof_variables,ll),             &
			    'N',                                              & ! transposition flag for gradients
			    'T'                                               & ! transposition flag for EIGVL
			)
		    END DO

	  ELSE
	  
	  DO LL=1,IELEM(N,I)%ADMIS
	      gradients(1:IELEM(N,I)%idegfree,:,ll) = ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%idegfree,1:nof_variables)
	  END DO
	  
	  DO LL=1,IELEM(N,I)%ADMIS
	      ! gradients_eigvlt = gradients * transpose(EIGVL)
	      CALL GEMM(                                            &
		  gradients(0:IELEM(N,I)%idegfree,1:nof_variables,ll),                    &
		  EIGVL(1:nof_variables,1:nof_variables),                           &
		  gradients_eigvlt(0:IELEM(N,I)%idegfree,1:nof_variables,ll),             &
		  'N',                                              & ! transposition flag for gradients
		  'T'                                               & ! transposition flag for EIGVL
	      )
	  END DO
    END IF
    
       IF (EES.NE.5)THEN
    
		      DO LL=1,IELEM(N,I)%ADMIS;do k=0,IELEM(N,I)%idegfree
			GRADCHARV(1:nof_variables,LL,k)=gradients_eigvlt(k,1:nof_variables,ll)
		      end do;END DO

		      ! SMOOTHINDICATOR(:,LL,L,1) = diag(transpose(gradients_eigvlt) * indicator * gradients_eigvlt)
		      A_CHAR(1:IELEM(N,I)%idegfree,1:nof_variables,1:IELEM(N,I)%ADMIS)=gradients_eigvlt(1:IELEM(N,I)%idegfree,1:nof_variables,1:IELEM(N,I)%ADMIS)
		      B_CHAR(1:IELEM(N,I)%idegfree,1:IELEM(N,I)%idegfree)=ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%idegfree,1:IELEM(N,I)%idegfree)
		      CALL diag_At_B_A(ICONSIDERED,                                     &
			A_CHAR(1:IELEM(N,I)%idegfree,1:nof_variables,1:IELEM(N,I)%ADMIS),              &
			B_CHAR(1:IELEM(N,I)%idegfree,1:IELEM(N,I)%idegfree), &
			X_CHAR(1:nof_variables,1:IELEM(N,I)%ADMIS)                          &
		      )
			SMOOTHINDICATOR(1:nof_variables,1:IELEM(N,I)%ADMIS,L,1)=X_CHAR(1:nof_variables,1:IELEM(N,I)%ADMIS)
        ELSE
        
         DO LL=1,1;do k=0,IELEM(N,I)%idegfree
         GRADCHARV(1:nof_variables,LL,k)=gradients_eigvlt(k,1:nof_variables,ll)
       end do;END DO
         DO LL=2,IELEM(N,I)%ADMIS;do k=0,IDEGFREE2
         GRADCHARV(1:nof_variables,LL,k)=gradients_eigvlt(k,1:nof_variables,ll)
      end do;END DO
      
      

      ! SMOOTHINDICATOR(:,LL,L,1) = diag(transpose(gradients_eigvlt) * indicator * gradients_eigvlt)
      
      A_CHAR(1:IELEM(N,I)%idegfree,1:nof_variables,1:IELEM(N,I)%ADMIS)=gradients_eigvlt(1:IELEM(N,I)%idegfree,1:nof_variables,1:IELEM(N,I)%ADMIS)
      B_CHAR(1:IELEM(N,I)%idegfree,1:IELEM(N,I)%idegfree)=ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%idegfree,1:IELEM(N,I)%idegfree)
      
      CALL diag_At_B_A(ICONSIDERED,                                     &
         A_CHAR(1:IELEM(N,I)%idegfree,1:nof_variables,1:IELEM(N,I)%ADMIS),              &
         B_CHAR(1:IELEM(N,I)%idegfree,1:IELEM(N,I)%idegfree), &
         X_CHAR(1:nof_variables,1:IELEM(N,I)%ADMIS)                          &
      )
        SMOOTHINDICATOR(1:nof_variables,1:IELEM(N,I)%ADMIS,L,1)=X_CHAR(1:nof_variables,1:IELEM(N,I)%ADMIS)
        END IF

end subroutine




SUBROUTINE WENOWEIGHTS(N)
!> @brief
!> Subroutine For WENO type reconstruction in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL::DIVISIONBYZERO
INTEGER::I,J,K,L,M,O,LL,IEX,IEUL,FACX,IELEME,KKD,KMAXE,JF,NGP,IQP,nnd,II,icd
INTEGER::IDUMMY,POWER,ITARGET
REAL::SUMOMEGAATILDEL
REAL::DIVBYZERO,COMPF,checkf,tau_Weno
REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T


KMAXE=XMPIELRANK(N)
DIVBYZERO=1.0e-6

POWER=4

call  QUADRATUREQUAD3D(N,IGQRULES);WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
call QUADRATURETRIANG(N,IGQRULES); WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)



!$OMP DO SCHEDULE (STATIC)
DO II=1,NOF_INTERIOR;I=EL_INT(II);ICONSIDERED=I
   IF (IELEM(N,I)%FULL.EQ.1)THEN
      CALL ALLGRADS_INNER(N,I)
      
      IF (WENWRT.EQ.2)THEN
         DO L=1,IELEM(N,I)%IFCA
	    !DEFINE
            ANGLE1=IELEM(N,I)%FACEANGLEX(L);ANGLE2=IELEM(N,I)%FACEANGLEY(L)
            NX=(COS(ANGLE1)*SIN(ANGLE2));NY=(SIN(ANGLE1)*SIN(ANGLE2));NZ=(COS(ANGLE2))
            VEIGL(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables);CALL ROTATEF(N,TRI,RVEIGL,VEIGL,ANGLE1,ANGLE2)
            VEIGR(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables);CALL ROTATEF(N,TRI,RVEIGR,VEIGR,ANGLE1,ANGLE2)
            CALL COMPUTE_EIGENVECTORS(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
            LAMBDA(:,:,L,1)=ZERO;SMOOTHINDICATOR(:,:,L,1)=ZERO;OMEGATILDE(:,:,L,1)=ZERO;OMEGA(:,:,L,1)=ZERO;FACEX=L
            CALL compute_gradcharv_smoothindicator(ICONSIDERED,facex)
	    LAMBDA(1:5,:,L,1)=1.0D0;LAMBDA(1:5,1,L,1)=LWCI1
	    
	    
		if (ees.eq.5)then;DO KKD=1,nof_variables
		LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                LAMBDA(KKD,1:ielem(n,i)%admis,L,1)=lamc(1:ielem(n,i)%admis)
                END DO;end if
		DO KKD=1,nof_variables
			    SUMOMEGATILDE(KKD)=ZERO
		    if (ees.eq.5)then
			tau_Weno=zero
			 if (wenoz.eq.1)then
							DO LL=1,IELEM(N,I)%ADMIS
							tau_Weno=tau_weno+(abs(SMOOTHINDICATOR(KKD,1,L,1)-SMOOTHINDICATOR(KKD,LL,L,1)))
							end do
							tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
							DO LL=1,IELEM(N,I)%ADMIS
							omegatilde(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATOR(KKD,LL,L,1))))
							end do
                                                    else
							DO LL=1,IELEM(N,I)%ADMIS
								OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
							END DO
						   end if
		    else
			DO LL=1,IELEM(N,I)%ADMIS
			    OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
			END DO
		    end if
				      
			    DO LL=1,IELEM(N,I)%ADMIS
				    SUMOMEGATILDE(KKD)=SUMOMEGATILDE(KKD)+OMEGATILDE(KKD,LL,L,1)
			    END DO
			    DO LL=1,IELEM(N,I)%ADMIS
				    OMEGA(KKD,LL,L,1)=(OMEGATILDE(KKD,LL,L,1))/SUMOMEGATILDE(KKD)
			    END DO
			    DO LL=1,IELEM(N,I)%ADMIS
			    WENOOS(KKD,LL,L,1)=OMEGA(KKD,LL,L,1)
			    END DO
	      END DO
		LIMITEDDW(:,:)=ZERO
		  IF (EES.EQ.5)THEN
			LIMITEDDW_CHAR(:,:,:)=ZERO
			DO LL=1,IELEM(N,I)%ADMIS;IF (LL.EQ.1)THEN
			ITARGET=IELEM(N,I)%idegfree
			ELSE
			ITARGET=IDEGFREE2
			END IF
			DO K=0,ITARGET
			LIMITEDDW_CHAR(1:nof_variables,K,1)=LIMITEDDW_CHAR(1:nof_variables,K,1)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW_CHAR(:,:,L,1,:)=ZERO
			
			DO K=0,IELEM(N,I)%idegfree     
				FINDW_CHAR(1:nof_variables,K,L,1,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW_CHAR(1:nof_variables,K,1))
			END DO
                 Else
			LIMITEDDW(:,:)=ZERO
			DO K=0,IELEM(N,I)%idegfree;DO LL=1,IELEM(N,I)%ADMIS
			LIMITEDDW(1:nof_variables,K)=LIMITEDDW(1:nof_variables,K)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW(:,:,L,1)=ZERO
			DO K=0,IELEM(N,I)%IDEGFREE
			FINDW(1:nof_variables,K,L,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW(1:nof_variables,K))
			END DO
                  end if

                  
                  
                  
                  
                  
                  
                  
                  
		IF (FASTEST_Q.EQ.1)THEN;if (ielem(n,i)%types_faces(L).eq.5)then;iqp=qp_quad;else;iqp=qp_triangle;end if
		ELSE
		    if (ielem(n,i)%types_faces(L).eq.5)then;iqp=qp_quad;NND=4
			    do K=1,nnd
			    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
			    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
			    END DO
			    call  QUADRATUREQUAD3D(N,IGQRULES)
		    else;iqp=QP_TRIANGLE;NND=3
			    do K=1,nnd
			    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
			    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
			    END DO
			    call QUADRATURETRIANG(N,IGQRULES)
		    end if
		END IF
			icd=0
			do NGP=1,iqp			
			    IF (FASTEST_Q.NE.1)THEN
			    AX = QPOINTS2D(1,NGP);AY = QPOINTS2D(2,NGP);AZ = QPOINTS2D(3,NGP)
			    ELSE
			    AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1);AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2);AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
			    END IF
			    icd=icd+1;COMPWRT=0
			    CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			    COMPWRT=0
			    if (ees.eq.5)then;COMPWRT=1
			    CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2)
			    COMPWRT=0;end if
			END DO		

			if (ees.eq.5)then;ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,:)=zero
			    DO NGP=1,iqp
			    ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)&
			    +FINDW_char(1:nof_Variables,0,L,1,1)
			    end do
! 			   
					    
							call gemm(                  &
							consmatrix(1:icd,1:ielem(n,i)%idegfree) ,     &
							FINDW_char(1:nof_variables,1:IELEM(N,I)%IDEGFREE,L,1,1),  &
							RESSOLUTION(1:ICD,1:NOF_vARIABLES),   &
							'N',                                                  & ! transposition flag for consmatrix
							    'T'                                                   & ! transposition flag for findw
								)
							icd=0;do NGP=1,iqp;icd=icd+1	    
							ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)&
							+RESSOLUTION(icd,1:NOF_vARIABLES)
							END DO
			else	
			      DO NGP=1,iqp	
			      ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=FINDW(1:nof_Variables,0,L,1)
			      end do	
					call gemm(                  &
					consmatrix(1:icd,1:ielem(n,i)%idegfree) ,     &
				      FINDW(1:nof_variables,1:IELEM(N,I)%IDEGFREE,L,1),  &
				      RESSOLUTION(1:ICD,1:NOF_vARIABLES),   &
				      'N',                                                  & ! transposition flag for consmatrix
					  'T'                                                   & ! transposition flag for findw
					      )
					      icd=0;do NGP=1,iqp;icd=icd+1	    
				      ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)&
				      +RESSOLUTION(icd,1:NOF_vARIABLES)
				      END DO	
			
			end if
		 END DO			!FACES
		 
		ELSE
		
		DO IEX=1,nof_variables
			LAMBDAAL=ZERO;SMOOTHINDICATORAL=ZERO;OMEGAATILDEL=ZERO;OMEGAAL=ZERO
                         IF (EES.EQ.5)THEN
				    LAMC(:)=ZERO; GRAD3AL(:)=ZERO; LAMC(1)=(1.0d0-(1.0d0/LWCI1));lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
				    LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
				    !sum the low degree polynomials first
				    DO LL=2,IELEM(N,I)%ADMIS
				    GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,IEX))
				    END DO
				    !this is the zero polynomial
				    GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
				    GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
				  DO LL=1,IELEM(N,I)%ADMIS
				    IF (LL.EQ.1)THEN
					call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
					      GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
					      INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
					      )
					SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
				    ELSE
					GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC(ll,1:IDEGFREE2,IEX)
					call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
					      GRAD1AL(1:IDEGFREE2),   &
					      INDICATEMATRIXAL(1:IDEGFREE2)     &
					      )
					SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
				    END IF
				  END DO
                         ELSE
			      DO LL=1,IELEM(N,I)%ADMIS
				  GRAD1AL(:)=ZERO
				  INDICATEMATRIXAL(:)=ZERO
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				  call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
				  SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
			      END DO
                         END IF
			      LAMBDAAL(:)=1.0D0
			      LAMBDAAL(1)=LWCI1
                                if (ees.eq.5)then
				LAMC(1)=(1.0d0-(1.0d0/LWCI1)); lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1);LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
				end if

			     if (ees.eq.5)then

				      if (wenoz.eq.1)then
					      tau_Weno=zero
					      DO LL=1,IELEM(N,I)%ADMIS
					      tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
					      end do
					      tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
					      DO LL=1,IELEM(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
					      end do
				      else
					      DO LL=1,IELEM(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					      END DO
				    
				      end if
                             else
					  DO LL=1,IELEM(N,I)%ADMIS
					  OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					  END DO
			      end if
			      
			      
			      SUMOMEGAATILDEL=ZERO
			      DO LL=1,IELEM(N,I)%ADMIS
			      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
			      END DO
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
			      END DO
			     
			      DO LL=1,IELEM(N,I)%ADMIS
			      WENO(IEX,LL)=OMEGAAL(LL)
			      
			      END DO
			    

			      
			      
			      
		  END DO
		  icd=0
		DO L=1,IELEM(N,I)%IFCA	!FACES
				 IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				  else
				    iqp=qp_triangle
				  end if
			    ELSE
		
		
				 if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
				  END IF
				  
			  do NGP=1,iqp			!for gqp
                                icd=icd+1
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
! 				DO K=1,IELEM(N,I)%IDEGFREE
                                        compwrt=0
	      				CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
	      				if (ees.eq.5)then
	      				compwrt=1
	      				CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2)
	      				compwrt=0
	      				end if
!                                 END DO
				end do
				
				

				
		    END DO	!FACES
!                                 DO IEX=1,nof_variables	!COMPONENTS
                                ILOCAL_RECON3(I)%ULEFT(:,:,:)=ZERO
! 				ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ZERO
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				
				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1:nof_variables)=GRAD5ALc(1:IELEM(N,I)%IDEGFREE,1:nof_variables)
				

				
                                call gemm(                                                  &
                            CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
                            GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
				
				
				else
				GRADSSL(1:Idegfree2,1:nof_variables)=ILOCAL_RECON5(1)%GRADIENTSc(LL,1:idegfree2,1:nof_variables)
				

				
                                call gemm(                                                  &
                            CONSMATRIXC(1:ICD,1:idegfree2),                &
                            GRADSSL(1:idegfree2,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
				
				
				
				end if
				
				
				
				
				ELSE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1:nof_variables)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,1:nof_variables)
				

				
                                call gemm(                                                  &
                            CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
                            GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
				
				
				END IF
				
				
                                 ICD=0
                                DO L=1,IELEM(N,I)%IFCA
                                    if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad;WEIGHT_T2(1:IQP)=WEIGHTS_Q(1:IQP)
                                    else
				    iqp=qp_triangle;WEIGHT_T2(1:IQP)=WEIGHTS_T(1:IQP)
                                    end if
                                     do NGP=1,iqp
                                        ICD=ICD+1
                                
                               CALL  EXTRAPOLATE_BOUND(IEX,L,NGP,I,ICD,LL)
                                
                                
!                                 ILOCAL_RECON3(I)%ULEFT(1:NOF_vARIABLES,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:NOF_vARIABLES,L,NGP)&
! 				    +((U_C(I)%VAL(1,1:NOF_vARIABLES)+RESSOLUTION(icd,1:NOF_vARIABLES) )*WENO(1:NOF_vARIABLES,LL))
				  end do
                                end do
                                
                                end do

		    
	END IF
		
	IF (((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) .and. (icoupleturb.eq.1)) THEN
		  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		    LAMBDAAL=ZERO
		    SMOOTHINDICATORAL=ZERO
		    OMEGAATILDEL=ZERO
		    OMEGAAL=ZERO

		   IF (EES.EQ.5)THEN
		     LAMC(:)=ZERO   
                          GRAD3AL(:)=ZERO
                            LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
		    
                            DO LL=2,IELEM(N,I)%ADMIS
                             GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:IDEGFREE2,IEX))
!                         
                             END DO
                            GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS2(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
                             GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
		    
                                DO LL=1,IELEM(N,I)%ADMIS
                                        IF (LL.EQ.1)THEN
                                            call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
                                            GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
                                            INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
                                            )
                                            SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                                        else
                                            GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC2(ll,1:IDEGFREE2,IEX)
                                        call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
                                            GRAD1AL(1:IDEGFREE2),   &
                                            INDICATEMATRIXAL(1:IDEGFREE2)     &
                                            )
                                        SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                                        
                                        
                                        
                                        end if
                                end do
                    else
		    
                                DO LL=1,IELEM(N,I)%ADMIS
                                GRAD1AL=ZERO
                                INDICATEMATRIXAL=ZERO
            ! 			  DO K=1,IELEM(N,I)%idegfree
                                    GRAD1AL(1:IELEM(N,I)%idegfree)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%idegfree,IEX)
            ! 			  END DO
                                call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
                                            GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
                                            INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
                                            )
                                SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                                END DO
		      
		      
                                LAMBDAAL(:)=1.0D0
                                LAMBDAAL(1)=LWCI1
                                
!                                             if (ees.eq.5)then
!                                         LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
!                                         end if
                                                if (ees.eq.5)then
                                                
                                                LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
                                                
                                                if (wenoz.eq.1)then
                                tau_Weno=zero
                                        DO LL=1,IELEM(N,I)%ADMIS
                                        tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
                                        end do
                                tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
                                        DO LL=1,IELEM(N,I)%ADMIS
                                        OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
                                        end do
                                else
                                DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      end if
                                                else
                                                DO LL=1,IELEM(N,I)%ADMIS
                                                OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
                                                END DO
                                            end if
                                SUMOMEGAATILDEL=ZERO
                                
                                DO LL=1,IELEM(N,I)%ADMIS
                                SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
                                END DO
                                DO LL=1,IELEM(N,I)%ADMIS
                                OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
                                END DO
                                DO LL=1,IELEM(N,I)%ADMIS
                                WENO2(IEX,LL)=OMEGAAL(LL)
                                END DO
                end if
		  END DO
		  
		  
		  
		DO L=1,IELEM(N,I)%IFCA	!FACES
		
			 IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				  else
				    iqp=qp_triangle
				  end if
			    ELSE
		
		
			if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
				  END IF
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
                                    
                                     compwrt=0
				 CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
				 
				  if (ees.eq.5)then
                                        compwrt=1
				 CONSMATRIXc(1,1:idegfree2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,idegfree2)
				 compwrt=0
	      				end if
				
! 	      				CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)

				DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR!COMPONENTS
				ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)=ZERO
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				else
				GRADSSL(1:IDEGFREE2,1)=ILOCAL_RECON5(1)%GRADIENTSc2(LL,1:IDEGFREE2,IEX)
				
				end if
				else
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				
				end if
				
! 	    			END DO
                                    IF (EES.EQ.5)THEN
                                        IF (LL.EQ.1)THEN
                                        RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                        else
                                        RESSOLUTION(1,1) = DOT(CONSMATRIXc(1,1:Idegfree2),GRADSSL(1:IDEGFREE2,1))
                                    
                                        end if
                                    else
                                     RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                    
                                    end if
                                    
                                    
	    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
	    ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)&
	    +((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO2(IEX,LL))
	 END DO		!STENCILS


				END DO	!COMPONENTS




				END DO	!NGP
				END DO	!FACES
		
	END IF
		ICONSIDERED=I
		CALL SOLUTIONTRIAV2(N,ICONSIDERED)
		
		
		
		
		
		
		
		
		
	END IF
	
	END DO
	!$OMP END DO 
	
	
	!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I
		IF (IELEM(N,I)%FULL.EQ.1)THEN
		CALL ALLGRADS_MIX(N,I)
		IF (WENWRT.EQ.2)THEN	
		  DO L=1,IELEM(N,I)%IFCA
				IDUMMY=0
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
				VEIGL(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				CALL ROTATEF(N,TRI,RVEIGL,VEIGL,ANGLE1,ANGLE2)
				
				
				IF (IELEM(N,I)%INEIGHB(l).EQ.N)THEN	!MY CPU ONLY
				      IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES	
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN MY CPU
					  VEIGR(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
					  IDUMMY=1
					  else
					  !NOT PERIODIC ONES IN MY CPU
					  facex=l;iconsidered=i
					   CALL coordinates_face_inner(N,Iconsidered,facex)
				  CORDS(1:3)=zero
 				  CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
			  
				  Poy(1)=cords(2)
				  Pox(1)=cords(1)
				  poz(1)=cords(3)
				  
 				  LEFTV(1:nof_variables)=VEIGL(1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
 				  CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
				  
				  VEIGR(1:nof_variables)=RIGHTV(1:nof_variables)
				      	  end if
				      ELSE
				      !FLUID NEIGHBOUR
				      VEIGR(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
				      END IF
				else
			      !other my cpu
				    IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					  VEIGR(1:nof_variables)=(IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_variables))
					  IDUMMY=1
					  end if
				    else
				  
				      VEIGR(1:nof_variables)=(IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_variables))
				  
				    end if
				    
				end if
				CALL ROTATEF(N,TRI,RVEIGR,VEIGR,ANGLE1,ANGLE2)
				CALL COMPUTE_EIGENVECTORS(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
				LAMBDA(:,:,L,1)=ZERO
				SMOOTHINDICATOR(:,:,L,1)=ZERO
				OMEGATILDE(:,:,L,1)=ZERO
				OMEGA(:,:,L,1)=ZERO

				
				 FACEX=L
                        CALL compute_gradcharv_smoothindicator(ICONSIDERED,facex)
				

					LAMBDA(1:5,:,L,1)=1.0D0
					LAMBDA(1:5,1,L,1)=LWCI1
				
				 if (ees.eq.5)then
				 
				 LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
				 
				 
                                    DO KKD=1,nof_variables
                                    LAMBDA(KKD,1:ielem(n,i)%admis,L,1)=lamc(1:ielem(n,i)%admis)
                                    END DO
                                  end if
				
				DO KKD=1,nof_variables
			    SUMOMEGATILDE(KKD)=ZERO
		    if (ees.eq.5)then
			tau_Weno=zero
			 if (wenoz.eq.1)then
							DO LL=1,IELEM(N,I)%ADMIS
							tau_Weno=tau_weno+(abs(SMOOTHINDICATOR(KKD,1,L,1)-SMOOTHINDICATOR(KKD,LL,L,1)))
							end do
							tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
							DO LL=1,IELEM(N,I)%ADMIS
							omegatilde(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATOR(KKD,LL,L,1))))
							end do
                                                    else
							DO LL=1,IELEM(N,I)%ADMIS
								OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
							END DO
						   end if
		    else
			DO LL=1,IELEM(N,I)%ADMIS
			    OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
			END DO
		    end if
				      
			    DO LL=1,IELEM(N,I)%ADMIS
				    SUMOMEGATILDE(KKD)=SUMOMEGATILDE(KKD)+OMEGATILDE(KKD,LL,L,1)
			    END DO
			    DO LL=1,IELEM(N,I)%ADMIS
				    OMEGA(KKD,LL,L,1)=(OMEGATILDE(KKD,LL,L,1))/SUMOMEGATILDE(KKD)
			    END DO
			    DO LL=1,IELEM(N,I)%ADMIS
			    WENOOS(KKD,LL,L,1)=OMEGA(KKD,LL,L,1)
			    END DO
	      END DO
		LIMITEDDW(:,:)=ZERO
		  IF (EES.EQ.5)THEN
			LIMITEDDW_CHAR(:,:,:)=ZERO
			DO LL=1,IELEM(N,I)%ADMIS;IF (LL.EQ.1)THEN
			ITARGET=IELEM(N,I)%idegfree
			ELSE
			ITARGET=IDEGFREE2
			END IF
			DO K=0,ITARGET
			LIMITEDDW_CHAR(1:nof_variables,K,1)=LIMITEDDW_CHAR(1:nof_variables,K,1)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW_CHAR(:,:,L,1,:)=ZERO
			
			DO K=0,IELEM(N,I)%idegfree     
				FINDW_CHAR(1:nof_variables,K,L,1,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW_CHAR(1:nof_variables,K,1))
			END DO
                 Else
			LIMITEDDW(:,:)=ZERO
			DO K=0,IELEM(N,I)%idegfree;DO LL=1,IELEM(N,I)%ADMIS
			LIMITEDDW(1:nof_variables,K)=LIMITEDDW(1:nof_variables,K)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW(:,:,L,1)=ZERO
			DO K=0,IELEM(N,I)%IDEGFREE
			FINDW(1:nof_variables,K,L,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW(1:nof_variables,K))
			END DO
                  end if
			
			IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				  ELSE
				     iqp=QP_TRIANGLE
				   END IF
		ELSE
			
			
			
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
                            end if
				  icd=0
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
			        icd=icd+1
			    
			    
!                                     DO K=1,IELEM(N,I)%IDEGFREE
                                    COMPWRT=0
				    CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
				    COMPWRT=0
				    if (ees.eq.5)then
				    COMPWRT=1
 				    CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2)
 				    COMPWRT=0
                                    end if
! 				    END DO
			    
			   
			END DO		!GAUSSIAN POINTS
			
			
			
			if (ees.eq.5)then;ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,:)=zero
			    DO NGP=1,iqp
			    ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)&
			    +FINDW_char(1:nof_Variables,0,L,1,1)
			    end do
! 			   
					    
							call gemm(                  &
							consmatrix(1:icd,1:ielem(n,i)%idegfree) ,     &
							FINDW_char(1:nof_variables,1:IELEM(N,I)%IDEGFREE,L,1,1),  &
							RESSOLUTION(1:ICD,1:NOF_vARIABLES),   &
							'N',                                                  & ! transposition flag for consmatrix
							    'T'                                                   & ! transposition flag for findw
								)
							icd=0;do NGP=1,iqp;icd=icd+1	    
							ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)&
							+RESSOLUTION(icd,1:NOF_vARIABLES)
							END DO
			else	
			      DO NGP=1,iqp	
			      ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=FINDW(1:nof_Variables,0,L,1)
			      end do	
					call gemm(                  &
					consmatrix(1:icd,1:ielem(n,i)%idegfree) ,     &
				      FINDW(1:nof_variables,1:IELEM(N,I)%IDEGFREE,L,1),  &
				      RESSOLUTION(1:ICD,1:NOF_vARIABLES),   &
				      'N',                                                  & ! transposition flag for consmatrix
					  'T'                                                   & ! transposition flag for findw
					      )
					      icd=0;do NGP=1,iqp;icd=icd+1	    
				      ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)&
				      +RESSOLUTION(icd,1:NOF_vARIABLES)
				      END DO	
			
			end if
		 END DO			!FACES
		 
		 
		 ELSE
		
		DO IEX=1,nof_variables
			LAMBDAAL=ZERO;SMOOTHINDICATORAL=ZERO;OMEGAATILDEL=ZERO;OMEGAAL=ZERO
			     IF (EES.EQ.5)THEN
				    LAMC(:)=ZERO; GRAD3AL(:)=ZERO; LAMC(1)=(1.0d0-(1.0d0/LWCI1));lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
				    LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
				    !sum the low degree polynomials first
				    DO LL=2,IELEM(N,I)%ADMIS
				    GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,IEX))
				    END DO
				    !this is the zero polynomial
				    GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
				    GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
				  DO LL=1,IELEM(N,I)%ADMIS
				    IF (LL.EQ.1)THEN
					call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
					      GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
					      INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
					      )
					SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
				    ELSE
					GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC(ll,1:IDEGFREE2,IEX)
					call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
					      GRAD1AL(1:IDEGFREE2),   &
					      INDICATEMATRIXAL(1:IDEGFREE2)     &
					      )
					SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
				    END IF
				  END DO
                         ELSE
			      DO LL=1,IELEM(N,I)%ADMIS
				  GRAD1AL(:)=ZERO
				  INDICATEMATRIXAL(:)=ZERO
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				  call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
				  SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
			      END DO
                         END IF

			     LAMBDAAL(:)=1.0D0
			      LAMBDAAL(1)=LWCI1
                                if (ees.eq.5)then
				LAMC(1)=(1.0d0-(1.0d0/LWCI1)); lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1);LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
				end if
                                
                                
                                 if (ees.eq.5)then

				      if (wenoz.eq.1)then
					      tau_Weno=zero
					      DO LL=1,IELEM(N,I)%ADMIS
					      tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
					      end do
					      tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
					      DO LL=1,IELEM(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
					      end do
				      else
					      DO LL=1,IELEM(N,I)%ADMIS
					      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					      END DO
				    
				      end if
                             else
					  DO LL=1,IELEM(N,I)%ADMIS
					  OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
					  END DO
			      end if
			      
			      SUMOMEGAATILDEL=ZERO
			      DO LL=1,IELEM(N,I)%ADMIS
			      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
			      END DO
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
			      END DO

			      DO LL=1,IELEM(N,I)%ADMIS
			      WENO(IEX,LL)=OMEGAAL(LL)
 			      
			      END DO
			      
			      
			     
			      
			      
			      
			      
			      
			      
			      
			      
			      
			      
			      
			      
		  END DO
		 ICD=0
		DO L=1,IELEM(N,I)%IFCA	!FACES
				  IF (FASTEST_Q.EQ.1)THEN
					if (ielem(n,i)%types_faces(L).eq.5)then
					  iqp=qp_quad
					ELSE
					  iqp=QP_TRIANGLE
					END IF
					
				  ELSE
		
				  IDUMMY=0
		
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
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
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    end if
				  end if
				  
			 do NGP=1,iqp			!for gqp
                                icd=icd+1
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
				compwrt=0
	      			CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
	      			if (ees.eq.5)then
	      			compwrt=1
	      				CONSMATRIXC(icd,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,IDEGFREE2)
                                compwrt=0
	      				end if
				end do
				
				

				
		    END DO	!FACES
!                                 DO IEX=1,nof_variables	!COMPONENTS
                                ILOCAL_RECON3(I)%ULEFT(:,:,:)=ZERO
! 				ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ZERO
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				
				IF (EES.EQ.5)THEN
				    IF (LL.EQ.1)THEN
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1:nof_variables)=GRAD5ALc(1:IELEM(N,I)%IDEGFREE,1:nof_variables)
				    call gemm(                                                  &
				CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
				RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
				    )
				    
				    else
				    GRADSSL(1:Idegfree2,1:nof_variables)=ILOCAL_RECON5(1)%GRADIENTSc(LL,1:idegfree2,1:nof_variables)
				    call gemm(                                                  &
				CONSMATRIXC(1:ICD,1:idegfree2),                &
				GRADSSL(1:idegfree2,1:NOF_vARIABLES),                                                &
				RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
				    )
				    end if
				
				
				
				
				ELSE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1:nof_variables)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,1:nof_variables)
				

				
                                call gemm(                                                  &
                            CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
                            GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
				
				
				END IF
                                 ICD=0
                                DO L=1,IELEM(N,I)%IFCA
                                    
                                 if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad;WEIGHT_T2(1:IQP)=WEIGHTS_Q(1:IQP)
                                    else
				    iqp=qp_triangle;WEIGHT_T2(1:IQP)=WEIGHTS_T(1:IQP)
                                    end if
                                     do NGP=1,iqp
                                        ICD=ICD+1
                                
                               CALL  EXTRAPOLATE_BOUND(IEX,L,NGP,I,ICD,LL)
                                
                               
!                                 ILOCAL_RECON3(I)%ULEFT(1:NOF_vARIABLES,L,NGP)=ILOCAL_RECON3(I)%ULEFT(1:NOF_vARIABLES,L,NGP)&
! 				    +((U_C(I)%VAL(1,1:NOF_vARIABLES)+RESSOLUTION(icd,1:NOF_vARIABLES) )*WENO(1:NOF_vARIABLES,LL))
                                    end do
                                end do
                                end do
		END IF!WENWRT2
		 
		 !TURBULENCE
		IF (((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) .and. (icoupleturb.eq.1)) THEN
		  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		    LAMBDAAL=ZERO
		    SMOOTHINDICATORAL=ZERO
		    OMEGAATILDEL=ZERO
		    OMEGAAL=ZERO
! 		    DO LL=1,IELEM(N,I)%ADMIS
! 		      GRAD1AL=ZERO
! 		      INDICATEMATRIXAL=ZERO
! 		       GRAD1AL(1:IELEM(N,I)%idegfree)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%idegfree,IEX)
! 			  
! 			  
! 			  call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
! 				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
! 				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
! 				  )
! 			  
! 			  SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
! 		      
! 		      
! 		      
! ! 			  DO K=1,IELEM(N,I)%idegfree
! ! 			  GRAD1AL(K)=ILOCAL_RECON5(1)%GRADIENTS2(LL,K,IEX)
! ! 			  END DO
! ! 		      INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE) = MATMUL(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),GRAD1AL(1:IELEM(N,I)%IDEGFREE))
! ! 		      SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
! 		    END DO
		    IF (EES.EQ.5)THEN
		     LAMC(:)=ZERO   
                          GRAD3AL(:)=ZERO
                             LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
		    
                            DO LL=2,IELEM(N,I)%ADMIS
                             GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:IDEGFREE2,IEX))
!                         
                             END DO
                            GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS2(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
                             GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
		    
                                DO LL=1,IELEM(N,I)%ADMIS
                                        IF (LL.EQ.1)THEN
                                            call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
                                            GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
                                            INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
                                            )
                                            SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                                        else
                                            GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC2(ll,1:IDEGFREE2,IEX)
                                        call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
                                            GRAD1AL(1:IDEGFREE2),   &
                                            INDICATEMATRIXAL(1:IDEGFREE2)     &
                                            )
                                        SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                                        
                                        
                                        
                                        end if
                                end do
                    else
		    
		    DO LL=1,IELEM(N,I)%ADMIS
		      GRAD1AL=ZERO
		      INDICATEMATRIXAL=ZERO
! 			  DO K=1,IELEM(N,I)%idegfree
			  GRAD1AL(1:IELEM(N,I)%idegfree)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%idegfree,IEX)
! 			  END DO
		      call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
		      SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
		    END DO
		      
		      
		      
		      LAMBDAAL(:)=1.0D0
		      LAMBDAAL(1)=LWCI1
		     
! 		   if (ees.eq.5)then
!                               LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
!                               end if
		   if (ees.eq.5)then
		   
		   
		    LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
                                IF (WENOZ.EQ.1)THEN
                                                tau_Weno=zero
                                                DO LL=1,IELEM(N,I)%ADMIS
                                                tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
                                                end do
                                                tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
                                                DO LL=1,IELEM(N,I)%ADMIS
                                                OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
                                                end do
                                                ELSE
						DO LL=1,IELEM(N,I)%ADMIS
						OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
						END DO
						ENDIF
                                else
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      end if
		    SUMOMEGAATILDEL=ZERO
		    DO LL=1,IELEM(N,I)%ADMIS
		      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
		    END DO
		    DO LL=1,IELEM(N,I)%ADMIS
		      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
		    END DO
		    DO LL=1,IELEM(N,I)%ADMIS
		      WENO2(IEX,LL)=OMEGAAL(LL)
		    END DO
		  END if
		  end do
		DO L=1,IELEM(N,I)%IFCA	!FACES
			IDUMMY=0
			
			IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				  ELSE
				     iqp=QP_TRIANGLE
				   END IF
		ELSE
			
			
		
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
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
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
				  end if
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
				
	      				  compwrt=0
				 CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
				 
				  if (ees.eq.5)then
                                        compwrt=1
				 CONSMATRIXc(1,1:IDEGFREE2)=BASIS_REC(N,AX,AY,AZ,IORDER2,I,idegfree2)
				 compwrt=0
	      				end if

				DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR!COMPONENTS
				ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)=ZERO
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				else
				GRADSSL(1:IDEGFREE2,1)=ILOCAL_RECON5(1)%GRADIENTSc2(LL,1:IDEGFREE2,IEX)
				
				end if
				else
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				
				end if
				
! 	    			END DO
                                    IF (EES.EQ.5)THEN
                                        IF (LL.EQ.1)THEN
                                        RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                        else
                                        RESSOLUTION(1,1) = DOT(CONSMATRIXc(1,1:Idegfree2),GRADSSL(1:Idegfree2,1))
                                    
                                        end if
                                    else
                                     RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                    
                                    end if
                                    
                                    
	    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
	    ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)&
	    +((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO2(IEX,LL))
	 END DO		!STENCILS


				END DO	!COMPONENTS




				END DO	!NGP
				END DO	!FACES
		
	END IF
		 
		ICONSIDERED=I
		CALL SOLUTIONTRIAV2(N,ICONSIDERED)	 
		 
	END IF
	END DO
!$OMP END DO  
	
END SUBROUTINE WENOWEIGHTS




SUBROUTINE WENOWEIGHTS2d(N)
!> @brief
!> Subroutine For WENO type reconstruction in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL::DIVISIONBYZERO
INTEGER::I,J,K,L,M,O,LL,IEX,IEUL,FACX,IELEME,KKD,KMAXE,JF,NGP,IQP,nnd,ii
INTEGER::IDUMMY,POWER,itarget
REAL::SUMOMEGAATILDEL,tau_Weno
REAL::DIVBYZERO,COMPF,checkf

KMAXE=XMPIELRANK(N)
DIVBYZERO=1.0e-6
! if (ees.eq.5)then
!     POWER=1
! ELSE
POWER=4
! END IF

	call  QUADRATURELINE(N,IGQRULES)
		
	!$OMP DO SCHEDULE (STATIC) 		
	DO II=1,NOF_INTERIOR
	I=EL_INT(II)
	ICONSIDERED=I
	IF (IELEM(N,I)%FULL.EQ.1)THEN
	
		CALL ALLGRADS_INNER2d(N,I)
		  IF (WENWRT.EQ.2)THEN		
		  DO L=1,IELEM(N,I)%IFCA
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
				VEIGL(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				CALL ROTATEF2d(N,TRI,RVEIGL,VEIGL,ANGLE1,ANGLE2)
				VEIGR(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				CALL ROTATEF2d(N,TRI,RVEIGR,VEIGR,ANGLE1,ANGLE2)
				CALL COMPUTE_EIGENVECTORS2d(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
				LAMBDA(:,:,L,1)=ZERO
				SMOOTHINDICATOR(:,:,L,1)=ZERO
				OMEGATILDE(:,:,L,1)=ZERO
				
				OMEGA(:,:,L,1)=ZERO
				  FACEX=L
                        CALL compute_gradcharv_smoothindicator(ICONSIDERED,facex)
				  
				  
				  
				  
					LAMBDA(1:4,:,L,1)=1.0D0
					LAMBDA(1:4,1,L,1)=LWCI1
					
					
					
					
                                  if (ees.eq.5)then
                                  
                                  LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            
                                  
                                            DO KKD=1,nof_variables
                                            LAMBDA(KKD,1:ielem(n,i)%admis,L,1)=lamc(1:ielem(n,i)%admis)
                                            END DO
                                  end if
					
					
					
				DO KKD=1,nof_variables
				    SUMOMEGATILDE(KKD)=ZERO
				         if (ees.eq.5)then
                                                    tau_Weno=zero
                                                    if (wenoz.eq.1)then
							DO LL=1,IELEM(N,I)%ADMIS
							tau_Weno=tau_weno+(abs(SMOOTHINDICATOR(KKD,1,L,1)-SMOOTHINDICATOR(KKD,LL,L,1)))
							end do
							tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
							DO LL=1,IELEM(N,I)%ADMIS
							omegatilde(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATOR(KKD,LL,L,1))))
							end do
                                                    else
							DO LL=1,IELEM(N,I)%ADMIS
								OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
							END DO
						   end if
                                            else
                                            DO LL=1,IELEM(N,I)%ADMIS
                                                    OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
                                            END DO
                                        end if
		      
				    DO LL=1,IELEM(N,I)%ADMIS
					    SUMOMEGATILDE(KKD)=SUMOMEGATILDE(KKD)+OMEGATILDE(KKD,LL,L,1)
				    END DO
				    DO LL=1,IELEM(N,I)%ADMIS
					    OMEGA(KKD,LL,L,1)=(OMEGATILDE(KKD,LL,L,1))/SUMOMEGATILDE(KKD)
				    END DO
				    DO LL=1,IELEM(N,I)%ADMIS
				    WENOOS(KKD,LL,L,1)=OMEGA(KKD,LL,L,1)
				    
				    END DO
				    if (kkd.eq.1)then
				    ielem(n,i)%vortex(1)=WENOOS(KKD,1,L,1)
				    end if
				END DO
				
                               IF (EES.EQ.5)THEN
			LIMITEDDW_CHAR(:,:,:)=ZERO
			DO LL=1,IELEM(N,I)%ADMIS;IF (LL.EQ.1)THEN
			ITARGET=IELEM(N,I)%idegfree
			ELSE
			ITARGET=IDEGFREE2
			END IF
			DO K=0,ITARGET
			LIMITEDDW_CHAR(1:nof_variables,K,1)=LIMITEDDW_CHAR(1:nof_variables,K,1)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW_CHAR(:,:,L,1,:)=ZERO
			
			DO K=0,IELEM(N,I)%idegfree     
				FINDW_CHAR(1:nof_variables,K,L,1,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW_CHAR(1:nof_variables,K,1))
			END DO
                 Else
			LIMITEDDW(:,:)=ZERO
			DO K=0,IELEM(N,I)%idegfree;DO LL=1,IELEM(N,I)%ADMIS
			LIMITEDDW(1:nof_variables,K)=LIMITEDDW(1:nof_variables,K)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW(:,:,L,1)=ZERO
			DO K=0,IELEM(N,I)%IDEGFREE
			FINDW(1:nof_variables,K,L,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW(1:nof_variables,K))
			END DO
                  end if
			
			      IF (FASTEST_Q.EQ.1)THEN
				  
				    iqp=qp_line
				 
		ELSE
				  
				    iqp=qp_line
				    NND=2
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  
				    call  QUADRATURELINE(N,IGQRULES)
				    
				    
				  end if  
				  
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
			  IF (EES.EQ.5)THEN
				compwrt=0
			     CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                compwrt=1
                                CONSMATRIXc(1,1:IDEGFREE2)=BASIS_REC2d(N,AX,AY,IORDER2,I,Idegfree2)
                                compwrt=0
                                    DO IEX=1,nof_variables
                                    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=zero
                                         ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)+FINDW_char(IEX,0,L,1,1)
                                        
                                         DO JF=1,IELEM(N,I)%IDEGFREE
                                        ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)+(CONSMATRIX(1,JF)*&
                                        FINDW_char(IEX,JF,L,1,1))
                                        END DO	
                                    END DO
                                else
                                    compwrt=0
			     CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                    DO IEX=1,nof_variables
                                    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=FINDW(IEX,0,L,1)
                                                    
                                    
                                        DO JF=1,IELEM(N,I)%IDEGFREE
                                        ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)+(CONSMATRIX(1,JF)*&
                                        FINDW(IEX,JF,L,1))
                                        END DO	

                                    END DO
                                
                                
                                
                                end if
			END DO		!GAUSSIAN POINTS
		 END DO			!FACES
		ELSE
		
		DO IEX=1,nof_variables
			LAMBDAAL=ZERO
			SMOOTHINDICATORAL=ZERO
			OMEGAATILDEL=ZERO
			OMEGAAL=ZERO
! 			      DO LL=1,IELEM(N,I)%ADMIS
! 				  GRAD1AL(:)=ZERO
! 				  INDICATEMATRIXAL(:)=ZERO
! 				  DO K=1,IELEM(N,I)%IDEGFREE
! 					  GRAD1AL(K)=ILOCAL_RECON5(1)%GRADIENTS(LL,K,IEX)
! 				  END DO
! 				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE) = MATMUL(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),GRAD1AL(1:IELEM(N,I)%IDEGFREE))
! 				  SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
! 			      END DO

                                 IF (EES.EQ.5)THEN
                          LAMC(:)=ZERO   
                          GRAD3AL(:)=ZERO
                            LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
                             
                             DO LL=2,IELEM(N,I)%ADMIS
                             GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,IEX))
!                         
                             END DO
                             
                             GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
                             GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
                             
                             DO LL=1,IELEM(N,I)%ADMIS
                             IF (LL.EQ.1)THEN
!                              GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX))
                             call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
!                             call gemv(ILOCAL_RECON3(I)%INDICATORc(1:IDEGFREE2,1:IDEGFREE2),          &
! 				  GRAD1AL(1:IDEGFREE2),   &
! 				  INDICATEMATRIXAL(1:IDEGFREE2)     &
! 				  )
                            
                             SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                              
                             
                            
                             ELSE
                             
                             GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC(ll,1:IDEGFREE2,IEX)
                             call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
				  GRAD1AL(1:IDEGFREE2),   &
				  INDICATEMATRIXAL(1:IDEGFREE2)     &
				  )
                             SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                            
                             END IF
                             END DO
                             
                             ELSE
			      DO LL=1,IELEM(N,I)%ADMIS
				  GRAD1AL(:)=ZERO
				  INDICATEMATRIXAL(:)=ZERO
! 				  DO K=1,IELEM(N,I)%IDEGFREE
					  GRAD1AL(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
! 				  END DO
				  call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
! 				 
				  SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
				  

			      END DO
                              END IF

			      LAMBDAAL(:)=1.0D0
			      LAMBDAAL(1)=LWCI1
                               if (ees.eq.5)then
                               
                               LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
                               end if

                                if (ees.eq.5)then
                                if (wenoz.eq.1)then
                                tau_Weno=zero
                                        DO LL=1,IELEM(N,I)%ADMIS
                                        tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
                                        end do
                                tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
                                        DO LL=1,IELEM(N,I)%ADMIS
                                        OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
                                        end do
                                else
                                DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      end if
                                    
                                else
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      end if
			      SUMOMEGAATILDEL=ZERO
			      DO LL=1,IELEM(N,I)%ADMIS
			      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
			      END DO
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
			      END DO
			      DO LL=1,IELEM(N,I)%ADMIS
			      WENO(IEX,LL)=OMEGAAL(LL)
 			    
			      END DO
			       if (iex.eq.1)then
				    ielem(n,i)%vortex(1)=WENO(iex,1)
				    end if
 		  END DO
		DO L=1,IELEM(N,I)%IFCA	!FACES
				 
				 ILOCAL_RECON3(I)%ULEFT(:,L,:)=ZERO
				    IF (FASTEST_Q.EQ.1)THEN
				  
				    iqp=qp_line
				 
		ELSE
				 
				    iqp=qp_line
				    NND=2
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    end if
				  
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
                                        compwrt=0
                                        
	      				CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                        if (ees.eq.5)then
                                        compwrt=1
	      				CONSMATRIXC(1,1:IDEGFREE2)=BASIS_REC2d(N,AX,AY,IORDER2,I,IDEGFREE2)
	      				compwrt=0
	      				end if

				DO IEX=1,nof_variables	!COMPONENTS
				
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				IF (EES.EQ.5)THEN
                                            IF (LL.EQ.1)THEN
                                            GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)
                                            else
				
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IDEGFREE2,1)=ILOCAL_RECON5(1)%GRADIENTSc(LL,1:IDEGFREE2,IEX)
				end if
				else
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				
				end if
				
				
! 	    			END DO
                                    IF (EES.EQ.5)THEN
                                    IF (LL.EQ.1)THEN

                                    
				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))

				   CALL EXTRAPOLATE_BOUND(IEX,L,NGP,I,LL,ll)
				    
				   
				    
				    
				    else
				     RESSOLUTION(1,1) = DOT(CONSMATRIXC(1,1:idegfree2),GRADSSL(1:idegfree2,1))
! 				     
				    
				    
				    
				   CALL EXTRAPOLATE_BOUND(IEX,L,NGP,I,LL,ll)
				    
				    end if
				    else
				     RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
! 				    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)&
! 				    +((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO(IEX,LL))
				    
				    
				    CALL EXTRAPOLATE_BOUND(IEX,L,NGP,I,LL,ll)
				    
				    
				    end if
				    
				END DO		!STENCILS
				END DO	!COMPONENTS
				END DO	!NGP
		    END DO	!FACES
	END IF
		
	IF (((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) .and. (icoupleturb.eq.1)) THEN
		  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		    LAMBDAAL=ZERO
		    SMOOTHINDICATORAL=ZERO
		    OMEGAATILDEL=ZERO
		    OMEGAAL=ZERO
		    
		    IF (EES.EQ.5)THEN
		     LAMC(:)=ZERO   
                          GRAD3AL(:)=ZERO
                             LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
		    
                            DO LL=2,IELEM(N,I)%ADMIS
                             GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:IDEGFREE2,IEX))
!                         
                             END DO
                            GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS2(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
                             GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
		    
                                DO LL=1,IELEM(N,I)%ADMIS
                                        IF (LL.EQ.1)THEN
                                            call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
                                            GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
                                            INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
                                            )
                                            SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                                        else
                                            GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC2(ll,1:IDEGFREE2,IEX)
                                        call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
                                            GRAD1AL(1:IDEGFREE2),   &
                                            INDICATEMATRIXAL(1:IDEGFREE2)     &
                                            )
                                        SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                                        
                                        
                                        
                                        end if
                                end do
                    else
		    
		    DO LL=1,IELEM(N,I)%ADMIS
		      GRAD1AL=ZERO
		      INDICATEMATRIXAL=ZERO
! 			  DO K=1,IELEM(N,I)%idegfree
			  GRAD1AL(1:IELEM(N,I)%idegfree)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%idegfree,IEX)
! 			  END DO
		      call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
		      SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
		    END DO
		   
		      
		      
		      LAMBDAAL(:)=1.0D0
		      LAMBDAAL(1)=LWCI1
! 		      if (ees.eq.5)then
!                               LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
!                               end if
		   if (ees.eq.5)then
		   LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
				if (wenoz.eq.1)then
                                tau_Weno=zero
                                DO LL=1,IELEM(N,I)%ADMIS
                                tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
                                end do
                                tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
                                DO LL=1,IELEM(N,I)%ADMIS
                                OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
                                end do
                                else
                                
                                DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      endif
                                else
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      end if
		    SUMOMEGAATILDEL=ZERO
		    DO LL=1,IELEM(N,I)%ADMIS
		      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
		    END DO
		    DO LL=1,IELEM(N,I)%ADMIS
		      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
		    END DO
		    DO LL=1,IELEM(N,I)%ADMIS
		      WENO2(IEX,LL)=OMEGAAL(LL)
		    END DO
		  END if
		  end do
		DO L=1,IELEM(N,I)%IFCA	!FACES
			   IF (FASTEST_Q.EQ.1)THEN
				  
				    iqp=qp_line
				 
                            ELSE
				    iqp=qp_line
				    NND=2
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				   end if 
				ILOCAL_RECON3(I)%ULEFTTURB(:,L,:)=ZERO    
				  
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
				
                                            
	      				 compwrt=0
                                        
	      				CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                        if (ees.eq.5)then
                                        compwrt=1
	      				CONSMATRIXC(1,1:IDEGFREE2)=BASIS_REC2d(N,AX,AY,IORDER2,I,IDEGFREE2)
	      				compwrt=0
	      				end if

				DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR!COMPONENTS
				
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				else
				GRADSSL(1:IDEGFREE2,1)=ILOCAL_RECON5(1)%GRADIENTSc2(LL,1:IDEGFREE2,IEX)
				
				end if
				else
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				
				end if
				
! 	    			END DO
                                    IF (EES.EQ.5)THEN
                                        IF (LL.EQ.1)THEN
                                        RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                        else
                                        RESSOLUTION(1,1) = DOT(CONSMATRIXc(1,1:Idegfree2),GRADSSL(1:Idegfree2,1))
                                    
                                        end if
                                    else
                                     RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                    
                                    end if
                                    
                                    
	    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
	   
	    
	    
			  if (reduce_comp.eq.1)then
				     ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,1)=ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,1)&
	    +((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO2(IEX,LL))*WEQUA2D(ngp)
				    
				    
				    else
				      ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)&
	    +((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO2(IEX,LL))
				    
				    
				    end if
	    
	    
	    
	    
	    
	 END DO		!STENCILS


				END DO	!COMPONENTS




				END DO	!NGP
				END DO	!FACES
		
	END IF
		ICONSIDERED=I
		
		CALL SOLUTIONTRIAV22d(N,ICONSIDERED)
		
		
		
		
		
		
		
		
		
	END IF
	END DO
	!$OMP END DO 
	
	!$OMP DO SCHEDULE (STATIC) 		
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I
	
	IF (IELEM(N,I)%FULL.EQ.1)THEN
	
		CALL ALLGRADS_MIX2d(N,I)
		IF (WENWRT.EQ.2)THEN	
		  DO L=1,IELEM(N,I)%IFCA
				IDUMMY=0
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
				VEIGL(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				CALL ROTATEF2d(N,TRI,RVEIGL,VEIGL,ANGLE1,ANGLE2)
				
				
				IF (IELEM(N,I)%INEIGHB(l).EQ.N)THEN	!MY CPU ONLY
				      IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES	
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN MY CPU
					  VEIGR(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
					  IDUMMY=1
					  else
					  !NOT PERIODIC ONES IN MY CPU
					  facex=l;iconsidered=i
					   CALL coordinates_face_inner2d(N,Iconsidered,facex)
				  CORDS(1:2)=zero
 				  CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
			  
				  Poy(1)=cords(2)
				  Pox(1)=cords(1)
				 
				  
 				  LEFTV(1:nof_variables)=VEIGL(1:nof_variables)
				  B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
 				  CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
				  
				  VEIGR(1:nof_variables)=RIGHTV(1:nof_variables)
				      	  end if
				      ELSE
				      !FLUID NEIGHBOUR
				      VEIGR(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
				      END IF
				else
			      !other my cpu
				    IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					  VEIGR(1:nof_variables)=(IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_variables))
					  IDUMMY=1
					  end if
				    else
				  
				      VEIGR(1:nof_variables)=(IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_variables))
				  
				    end if
				    
				end if
				CALL ROTATEF2d(N,TRI,RVEIGR,VEIGR,ANGLE1,ANGLE2)
				CALL COMPUTE_EIGENVECTORS2d(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
				LAMBDA(:,:,L,1)=ZERO
				SMOOTHINDICATOR(:,:,L,1)=ZERO
				OMEGATILDE(:,:,L,1)=ZERO
				OMEGA(:,:,L,1)=ZERO
				 
				  FACEX=L
                        CALL compute_gradcharv_smoothindicator(ICONSIDERED,facex)
				 
				 
					LAMBDA(1:4,:,L,1)=1.0D0
					LAMBDA(1:4,1,L,1)=LWCI1
				      if (ees.eq.5)then
					  
                                  LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                                  
                                            DO KKD=1,nof_variables
                                            LAMBDA(KKD,1:ielem(n,i)%admis,L,1)=lamc(1:ielem(n,i)%admis)
                                            END DO
				      end if					
				DO KKD=1,nof_variables
				    SUMOMEGATILDE(KKD)=ZERO
				   
				     if (ees.eq.5)then
				     IF (WENOZ.EQ.1)THEN
				      tau_Weno=zero
				      DO LL=1,IELEM(N,I)%ADMIS
				      tau_Weno=tau_weno+(abs(SMOOTHINDICATOR(KKD,1,L,1)-SMOOTHINDICATOR(KKD,LL,L,1)))
				      end do
				      tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
					  DO LL=1,IELEM(N,I)%ADMIS
					  omegatilde(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATOR(KKD,LL,L,1))))
					  end do
                                 ELSE
                                    DO LL=1,IELEM(N,I)%ADMIS
					    OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
				    END DO
				  ENDIF
                                else
                                    DO LL=1,IELEM(N,I)%ADMIS
					    OMEGATILDE(KKD,LL,L,1)=(LAMBDA(KKD,LL,L,1))/((DIVBYZERO+SMOOTHINDICATOR(KKD,LL,L,1))**POWER)
				    END DO
			      end if
		      
				    DO LL=1,IELEM(N,I)%ADMIS
					    SUMOMEGATILDE(KKD)=SUMOMEGATILDE(KKD)+OMEGATILDE(KKD,LL,L,1)
				    END DO
				    DO LL=1,IELEM(N,I)%ADMIS
					    OMEGA(KKD,LL,L,1)=(OMEGATILDE(KKD,LL,L,1))/SUMOMEGATILDE(KKD)
				    END DO
				    DO LL=1,IELEM(N,I)%ADMIS
				    WENOOS(KKD,LL,L,1)=OMEGA(KKD,LL,L,1)
				   
				    END DO
				     if (kkd.eq.1)then
				    ielem(n,i)%vortex(1)=WENOOS(KKD,1,L,1)
				    end if
				END DO
                            
                            
                          
                                            
                            IF (EES.EQ.5)THEN
			LIMITEDDW_CHAR(:,:,:)=ZERO
			DO LL=1,IELEM(N,I)%ADMIS;IF (LL.EQ.1)THEN
			ITARGET=IELEM(N,I)%idegfree
			ELSE
			ITARGET=IDEGFREE2
			END IF
			DO K=0,ITARGET
			LIMITEDDW_CHAR(1:nof_variables,K,1)=LIMITEDDW_CHAR(1:nof_variables,K,1)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW_CHAR(:,:,L,1,:)=ZERO
			
			DO K=0,IELEM(N,I)%idegfree     
				FINDW_CHAR(1:nof_variables,K,L,1,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW_CHAR(1:nof_variables,K,1))
			END DO
                 Else
			LIMITEDDW(:,:)=ZERO
			DO K=0,IELEM(N,I)%idegfree;DO LL=1,IELEM(N,I)%ADMIS
			LIMITEDDW(1:nof_variables,K)=LIMITEDDW(1:nof_variables,K)+GRADCHARV(1:nof_variables,LL,K)*WENOOS(1:nof_variables,LL,L,1)
			END DO;END DO
			FINDW(:,:,L,1)=ZERO
			DO K=0,IELEM(N,I)%IDEGFREE
			FINDW(1:nof_variables,K,L,1)=MATMUL(EIGVR(1:nof_variables,1:nof_variables),LIMITEDDW(1:nof_variables,K))
			END DO
                  end if
				
				
			   IF (FASTEST_Q.EQ.1)THEN
				  
				    iqp=qp_line
				 
		ELSE
			
				 
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    end if
				    
				 
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				 IF (EES.EQ.5)THEN
				compwrt=0
			     CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                compwrt=1
                                CONSMATRIXc(1,1:IDEGFREE2)=BASIS_REC2d(N,AX,AY,IORDER2,I,Idegfree2)
                                compwrt=0
                                    DO IEX=1,nof_variables
                                    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=zero
                                         ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)+FINDW_char(IEX,0,L,1,1)
                                        
                                         DO JF=1,IELEM(N,I)%IDEGFREE
                                        ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)+(CONSMATRIX(1,JF)*&
                                        FINDW_char(IEX,JF,L,1,1))
                                        END DO	
                                    END DO
                                else
                                    compwrt=0
			     CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                    DO IEX=1,nof_variables
                                    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=FINDW(IEX,0,L,1)
                                                    
                                    
                                        DO JF=1,IELEM(N,I)%IDEGFREE
                                        ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)+(CONSMATRIX(1,JF)*&
                                        FINDW(IEX,JF,L,1))
                                        END DO	

                                    END DO
                                
                                
                                
                                end if
			END DO		!GAUSSIAN POINTS
		 END DO			!FACES
		 
		 
		 ELSE
		
		
		DO IEX=1,nof_variables
			LAMBDAAL=ZERO
			SMOOTHINDICATORAL=ZERO
			OMEGAATILDEL=ZERO
			OMEGAAL=ZERO
! 			      DO LL=1,IELEM(N,I)%ADMIS
! 				  GRAD1AL(:)=ZERO
! 				  INDICATEMATRIXAL(:)=ZERO
! 				  DO K=1,IELEM(N,I)%IDEGFREE
! 					  GRAD1AL(K)=ILOCAL_RECON5(1)%GRADIENTS(LL,K,IEX)
! 				  END DO
! 				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE) = MATMUL(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),GRAD1AL(1:IELEM(N,I)%IDEGFREE))
! 				  SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
! 			      END DO
                                  IF (EES.EQ.5)THEN
                          LAMC(:)=ZERO   
                          GRAD3AL(:)=ZERO
                          LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
                             DO LL=2,IELEM(N,I)%ADMIS
                             GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC(LL,1:IDEGFREE2,IEX))
                             END DO
                             
                             GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
                             GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
                             
                             DO LL=1,IELEM(N,I)%ADMIS
                             IF (LL.EQ.1)THEN
                             
                             call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
                                
                             SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                             
                             
!                               call gemv(ILOCAL_RECON3(I)%INDICATORc(1:IDEGFREE2,1:IDEGFREE2),          &
! 				  GRAD1AL(1:IDEGFREE2),   &
! 				  INDICATEMATRIXAL(1:IDEGFREE2)     &
! 				  )
!                             
!                              SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                             
                             
                             
                             
                             ELSE
                             
                             GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC(ll,1:IDEGFREE2,IEX)
                             call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
				  GRAD1AL(1:IDEGFREE2),   &
				  INDICATEMATRIXAL(1:IDEGFREE2)     &
				  )
                             SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                             
                             
                             END IF
                             END DO
                             
                             ELSE
			      DO LL=1,IELEM(N,I)%ADMIS
				  GRAD1AL(:)=ZERO
				  INDICATEMATRIXAL(:)=ZERO
! 				  DO K=1,IELEM(N,I)%IDEGFREE
					  GRAD1AL(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
! 				  END DO
				  call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
! 				 
				  SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
				  
				  
				  
! 				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE) = MATMUL(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),GRAD1AL(1:IELEM(N,I)%IDEGFREE))
! 				  SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
			      END DO
                              END IF
			      LAMBDAAL(:)=1.0D0
			      LAMBDAAL(1)=LWCI1
!                                 if (ees.eq.5)then
!                               LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
!                               end if
!                                  LAMBDAAL(:)=1.0D0
! 			      LAMBDAAL(1)=LWCI1
                            if (ees.eq.5)then
                            LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
                            
                            
				IF (WENOZ.EQ.1)THEN
                                tau_Weno=zero
                                DO LL=1,IELEM(N,I)%ADMIS
                                tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
                                end do
                                tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
                                DO LL=1,IELEM(N,I)%ADMIS
                                OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
                                end do
                                ELSE
                                
                               DO LL=1,IELEM(N,I)%ADMIS
                                   OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
                                   END DO
				END IF
                                
                                else
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
			      END DO
			      
			      
			      
			      end if
			     
			      SUMOMEGAATILDEL=ZERO
			      DO LL=1,IELEM(N,I)%ADMIS
			      SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
			      END DO
			      DO LL=1,IELEM(N,I)%ADMIS
			      OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
			      END DO
			      DO LL=1,IELEM(N,I)%ADMIS
			      WENO(IEX,LL)=OMEGAAL(LL)
                                    
			      
			     
			      END DO
                                    if (iex.eq.1)then
				    ielem(n,i)%vortex(1)=WENO(IEX,1)
				    end if
		  END DO
		  
		  
		  
		DO L=1,IELEM(N,I)%IFCA	!FACES
		
		     ILOCAL_RECON3(I)%ULEFT(:,L,:)=zero
		
		
				  IDUMMY=0
				  
				     IF (FASTEST_Q.EQ.1)THEN
				  
				    iqp=qp_line
				 
		ELSE
				  
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    end if
				  
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
				
                                            compwrt=0
	      				  CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
	      				  
                                            if (ees.eq.5)then
                                            compwrt=1
	      				CONSMATRIXC(1,1:IDEGFREE2)=BASIS_REC2d(N,AX,AY,IORDER2,I,IDEGFREE2)
                                            compwrt=0
! 	      				
	      				end if

! 				DO IEX=1,nof_variables	!COMPONENTS
! 				ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ZERO
! 				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
! ! 	    			DO K=1,IELEM(N,I)%IDEGFREE
! 				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
! ! 	    			END DO
! 				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
! 				    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
! 				    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)&
! 				    +((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO(IEX,LL))
! 				END DO		!STENCILS
                                DO IEX=1,nof_variables	!COMPONENTS
				ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ZERO
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)
				else
				
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IDEGFREE2,1)=ILOCAL_RECON5(1)%GRADIENTSc(LL,1:IDEGFREE2,IEX)
				end if
				else
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				
				end if
				
				
! 	    			END DO
                                    IF (EES.EQ.5)THEN
                                    IF (LL.EQ.1)THEN

                                    
				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
				    
				    
				   CALL EXTRAPOLATE_BOUND(IEX,L,NGP,I,LL,ll)
				    

				    else
				     RESSOLUTION(1,1) = DOT(CONSMATRIXC(1,1:idegfree2),GRADSSL(1:idegfree2,1))
!                                     
				   CALL EXTRAPOLATE_BOUND(IEX,L,NGP,I,LL,ll)
				    
				    end if
				    else
				     RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
				  CALL EXTRAPOLATE_BOUND(IEX,L,NGP,I,LL,ll)
				    
				    end if
				    
				END DO		!STENCILS
				END DO	!COMPONENTS
				END DO	!NGP
		    END DO	!FACES
		END IF!WENWRT2
		 
		 !TURBULENCE
		IF (((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) .and. (icoupleturb.eq.1)) THEN
		  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
		    LAMBDAAL=ZERO
		    SMOOTHINDICATORAL=ZERO
		    OMEGAATILDEL=ZERO
		    OMEGAAL=ZERO
		    IF (EES.EQ.5)THEN
		     LAMC(:)=ZERO   
                          GRAD3AL(:)=ZERO
                            LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
		    
                            DO LL=2,IELEM(N,I)%ADMIS
                             GRAD3AL(1:IDEGFREE2)=GRAD3AL(1:IDEGFREE2)+(LAMC(LL)*ILOCAL_RECON5(1)%GRADIENTSC2(LL,1:IDEGFREE2,IEX))
!                         
                             END DO
                            GRAD1AL(1:IELEM(N,I)%IDEGFREE)=(1.0D0/LAMC(1))*(ILOCAL_RECON5(1)%GRADIENTS2(1,1:IELEM(N,I)%IDEGFREE,IEX)-GRAD3AL(1:IELEM(N,I)%IDEGFREE))
                             GRAD5ALc(1:IELEM(N,I)%IDEGFREE,iex)=GRAD1AL(1:IELEM(N,I)%IDEGFREE)
		    
                                DO LL=1,IELEM(N,I)%ADMIS
                                        IF (LL.EQ.1)THEN
                                            call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
                                            GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
                                            INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
                                            )
                                            SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
                                        else
                                            GRAD1AL(1:IDEGFREE2)=ILOCAL_RECON5(1)%GRADIENTSC2(ll,1:IDEGFREE2,IEX)
                                        call gemv(ILOCAL_RECON3(I)%INDICATORC(1:IDEGFREE2,1:IDEGFREE2),          &
                                            GRAD1AL(1:IDEGFREE2),   &
                                            INDICATEMATRIXAL(1:IDEGFREE2)     &
                                            )
                                        SMOOTHINDICATORAL(LL)= DOT(GRAD1AL(1:IDEGFREE2),INDICATEMATRIXAL(1:IDEGFREE2))
                                        
                                        
                                        
                                        end if
                                end do
                    else
		    
		    DO LL=1,IELEM(N,I)%ADMIS
		      GRAD1AL=ZERO
		      INDICATEMATRIXAL=ZERO
! 			  DO K=1,IELEM(N,I)%idegfree
			  GRAD1AL(1:IELEM(N,I)%idegfree)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%idegfree,IEX)
! 			  END DO
		      call gemv(ILOCAL_RECON3(I)%INDICATOR(1:IELEM(N,I)%IDEGFREE,1:IELEM(N,I)%IDEGFREE),          &
				  GRAD1AL(1:IELEM(N,I)%IDEGFREE),   &
				  INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE)     &
				  )
		      SMOOTHINDICATORAL(LL)= DOT_PRODUCT(GRAD1AL(1:IELEM(N,I)%IDEGFREE),INDICATEMATRIXAL(1:IELEM(N,I)%IDEGFREE))
		    END DO
		   
		      
		      
		      LAMBDAAL(:)=1.0D0
		      LAMBDAAL(1)=LWCI1
		     
!                                 if (ees.eq.5)then
!                               LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
!                               end if
                                    if (ees.eq.5)then
                                    
                                    
                                    LAMC(1)=(1.0d0-(1.0d0/LWCI1))
                               lamc(2:ielem(n,i)%admis)=(1.0d0-lamc(1))/(IELEM(N,I)%ADMIS-1)
                            LAMBDAAL(1:ielem(n,i)%admis)=lamc(1:ielem(n,i)%admis)
						IF (WENOZ.EQ.1)THEN
                                                tau_Weno=zero
                                                DO LL=1,IELEM(N,I)%ADMIS
                                                tau_Weno=tau_weno+(abs(SMOOTHINDICATORAL(1)-SMOOTHINDICATORAL(LL)))
                                                end do
                                                tau_weno=(tau_weno/(IELEM(N,I)%ADMIS-1))**power
                                                DO LL=1,IELEM(N,I)%ADMIS
                                                OMEGAATILDEL(LL)=(LAMBDAAL(LL))*(1.0d0+(tau_weno/(divbyzero+SMOOTHINDICATORAL(LL))))
                                                end do
                                                ELSE
						DO LL=1,IELEM(N,I)%ADMIS
						OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
						END DO
						ENDIF
                                    else
                                    DO LL=1,IELEM(N,I)%ADMIS
                                    OMEGAATILDEL(LL)=(LAMBDAAL(LL))/((DIVBYZERO+SMOOTHINDICATORAL(LL))**POWER)
                                    END DO
                                end if
                                SUMOMEGAATILDEL=ZERO
                                DO LL=1,IELEM(N,I)%ADMIS
                                SUMOMEGAATILDEL=SUMOMEGAATILDEL+OMEGAATILDEL(LL)
                                END DO
                                DO LL=1,IELEM(N,I)%ADMIS
                                OMEGAAL(LL)=(OMEGAATILDEL(LL))/SUMOMEGAATILDEL
                                END DO
                                DO LL=1,IELEM(N,I)%ADMIS
                                WENO2(IEX,LL)=OMEGAAL(LL)
                                END DO
                    end if
		  END DO
		DO L=1,IELEM(N,I)%IFCA	!FACES
			IDUMMY=0
			ILOCAL_RECON3(I)%ULEFTTURB(:,L,:)=ZERO
			
				   IF (FASTEST_Q.EQ.1)THEN
				  
				    iqp=qp_line
				 
		ELSE
			
		
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					  ELSE
					  facex=l;iconsidered=i
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    end if
				    
				 
			  do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
				
				
	      				 compwrt=0
                                        
	      				CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                                        if (ees.eq.5)then
                                        compwrt=1
	      				CONSMATRIXC(1,1:IDEGFREE2)=BASIS_REC2d(N,AX,AY,IORDER2,I,IDEGFREE2)
	      				compwrt=0
	      				end if

				DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR!COMPONENTS
				
				DO LL=1,IELEM(N,I)%ADMIS	!STENCILS
				IF (EES.EQ.5)THEN
				IF (LL.EQ.1)THEN
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				else
				GRADSSL(1:IDEGFREE2,1)=ILOCAL_RECON5(1)%GRADIENTSc2(LL,1:IDEGFREE2,IEX)
				
				end if
				else
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS2(LL,1:IELEM(N,I)%IDEGFREE,IEX)
				
				end if
				
! 	    			END DO
                                    IF (EES.EQ.5)THEN
                                        IF (LL.EQ.1)THEN
                                        RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                        else
                                        RESSOLUTION(1,1) = DOT(CONSMATRIXc(1,1:Idegfree2),GRADSSL(1:IDEGFREE2,1))
                                    
                                        end if
                                    else
                                     RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
                                    
                                    end if
                                    
                                    
	    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
	   
	    
	    IF (REDUCE_COMP.EQ.1)THEN
				     ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,1)=ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,1)&
	    +((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO2(IEX,LL))*WEQUA2D(NGP)
				    
				    
				    ELSE
				     ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)=ILOCAL_RECON3(I)%ULEFTTURB(IEX,L,NGP)&
	    +((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1))*WENO2(IEX,LL))
				    
				    END IF
	    
	    
	    
	 END DO		!STENCILS


				END DO	!COMPONENTS




				END DO	!NGP
				END DO	!FACES
		
! 	END IF
		 
		 
		 
	END IF
	
		 ICONSIDERED=I
	 CALL SOLUTIONTRIAV22d(N,ICONSIDERED)
	
	
	end if
	
	END DO
!$OMP END DO  
	
END SUBROUTINE WENOWEIGHTS2d
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

!$OMP DO SCHEDULE (STATIC) 		
	DO I=1,KMAXE
! 	if (ielem(n,i)%interior.ne.1)then
! 		CALL ALLGRADS_INNER(N,I)
! 	ELSE
! 		CALL ALLGRADS_MIX(N,I)
! 	END IF
	DO IEX=1,nof_variables
 	ILOCAL_RECON3(I)%ULEFT(IEX,:,:)=U_C(I)%VAL(1,IEX)
	END DO
	
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	DO IEX=1,turbulenceequations+passivescalar
	ILOCAL_RECON3(I)%ULEFTTURB(IEX,:,:)=U_Ct(I)%VAL(1,IEX)
	END DO
	end if
	ICONSIDERED=I
	 CALL SOLUTIONTRIAV2(N,ICONSIDERED)
	END DO
!$OMP END DO 
 
END SUBROUTINE PIECEWISE_CONSTANT


SUBROUTINE LINEAR_SCHEME(N)
!> @brief
!> Subroutine for unlimited linear scheme
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX,ICD
REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T
KMAXE=XMPIELRANK(N)

call  QUADRATUREQUAD3D(N,IGQRULES)
	
	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
	
	call QUADRATURETRIANG(N,IGQRULES)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)



!$OMP DO SCHEDULE (STATIC) 		
	DO I=1,KMAXE
	ICONSIDERED=I
	if (ielem(n,i)%interior.ne.1)then
		CALL ALLGRADS_INNER(N,I)
	ELSE
		CALL ALLGRADS_MIX(N,I)
	END IF
            ICD=0
	
	      DO L=1,IELEM(N,I)%IFCA
		IDUMMY=0
		ILOCAL_RECON3(I)%ULEFT(:,L,:)=ZERO
		if (fastest_q.ne.1)then
				  
		
		
		
		
		    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
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
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
		else
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
		end if
		else
		if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad;WEIGHT_T2(1:IQP)=WEIGHTS_Q(1:IQP)
		else
				    iqp=qp_triangle;;WEIGHT_T2(1:IQP)=WEIGHTS_T(1:IQP)
		end if
		
		
		
		end if
		
		
		
		
		
!                              do NGP=1,iqp			!for gqp
! 				if (fastest_q.eq.1)then
! 				ax= ilocal_recon3(i)%qpoints(l,ngp,1)
! 				ay= ilocal_recon3(i)%qpoints(l,ngp,2)
! 				az=ilocal_recon3(i)%qpoints(l,ngp,3)
! 				else
! 				AX = QPOINTS2D(1,NGP)
! 				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
! 				end if
! 				CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
! 	
! 				
! 				DO IEX=1,nof_variables	!COMPONENTS
! 				ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ZERO
! 				
! ! 	    			DO K=1,IELEM(N,I)%IDEGFREE
! 				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)
! 				
! ! 	    			END DO
! 				  
! 				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
! 				    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
! 				    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)
! 				    
! 				END DO	!COMPONENTS
! 				END DO	!NGP
		
			     do NGP=1,iqp
			     ICD=ICD+1
                                !for gqp
				if (fastest_q.eq.1)then
				ax= ilocal_recon3(i)%qpoints(l,ngp,1)
				ay= ilocal_recon3(i)%qpoints(l,ngp,2)
				az=ilocal_recon3(i)%qpoints(l,ngp,3)
				else
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				end if
				CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                            END DO
				
! 				DO IEX=1,nof_variables	!COMPONENTS
! 				ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=ZERO
! 				
! ! 	    			DO K=1,IELEM(N,I)%IDEGFREE
! 				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)
! 				
! ! 	    			END DO
! 				  
! 				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
! 				    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
! 				    ILOCAL_RECON3(I)%ULEFT(IEX,L,NGP)=U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)
! 				END DO	!COMPONENTS
! 				END DO	!NGP
	END DO
	
				
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)
                                call gemm(                                                  &
                            CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
                            GRADSSL(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
                                ICD=0
                                DO L=1,IELEM(N,I)%IFCA
                                    if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
                                    else
				    iqp=qp_triangle
                                    end if
                                     do NGP=1,iqp
			     ICD=ICD+1
				    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
! 				    ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)=U_C(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(ICD,1:NOF_VARIABLES)
				    CALL EXTRAPOLATE_BOUND_LINEAR(IEX,L,NGP,I,ICD)
				   
				    
				    
				    
				    
! 				   
				    END DO
				END DO	!COMPONENTS
! 				STOP
! 				END DO	!NGP
	
	    ICONSIDERED=I
	 CALL SOLUTIONTRIAV2(N,ICONSIDERED)
	  
	  END DO
!$OMP END DO 
 
END SUBROUTINE LINEAR_SCHEME



SUBROUTINE PIECEWISE_CONSTANT2d(N)
!> @brief
!> Subroutine For first-order scheme in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,KMAXE,IEX
KMAXE=XMPIELRANK(N)
!$OMP DO SCHEDULE (STATIC) 		
	DO I=1,KMAXE
	DO IEX=1,nof_variables
	ILOCAL_RECON3(I)%ULEFT(IEX,:,:)=U_C(I)%VAL(1,IEX)
	END DO
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	DO IEX=1,turbulenceequations+passivescalar
	ILOCAL_RECON3(I)%ULEFTTURB(IEX,:,:)=U_Ct(I)%VAL(1,IEX)
	END DO
	end if
	ICONSIDERED=I
	 CALL SOLUTIONTRIAV22d(N,ICONSIDERED)
	END DO
!$OMP END DO 
 
END SUBROUTINE PIECEWISE_CONSTANT2d


SUBROUTINE LINEAR_SCHEME2d(N)
!> @brief
!> Subroutine for unlimited linear scheme in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX,IDFR
! REAL::NDMAX,NDMIN,NFMIN,NFMAX
KMAXE=XMPIELRANK(N)

call  QUADRATUREline(N,IGQRULES)

!$OMP DO SCHEDULE (STATIC) 		
	DO I=1,KMAXE
	ICONSIDERED=I
! 	NDMAX=ZERO;NDMIN=10.0E18;NFMIN=10.0E18;NFMAX=ZERO
! 	
	if (ielem(n,i)%interior.ne.1)then
		CALL ALLGRADS_INNER2d(N,I)
	ELSE
		CALL ALLGRADS_MIX2d(N,I)
	END IF
	
	
	      DO L=1,IELEM(N,I)%IFCA
	      ILOCAL_RECON3(I)%ULEFT(:,L,:)=ZERO
	      
		if (fastest_q.ne.1)then
	      
		IDUMMY=0
		    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				 
		else
				  iqp=qp_line
				    NND=2
					  
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  
					  
				    call  QUADRATUREline(N,IGQRULES)
       
       
		end if
		
		
		
		ELSE
		
		IQP=QP_LINE
		
		END IF
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
				
! 				IF (INITCOND.EQ.0)THEN
! 				DO IDFR=1,IELEM(N,I)%IDEGFREE
! 				IF (CONSMATRIX(1,IDFR).NE.ZERO)THEN
! 				NFMIN=MIN(NFMIN,CONSMATRIX(1,IDFR))
! 				NFMAX=MAX(NFMAX,CONSMATRIX(1,IDFR))
! 				
! 				END IF
! 				END DO
				
				DO IEX=1,nof_variables	!COMPONENTS
				
				
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,IEX)
! 	    			END DO
				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
				    
				    
				    
				    
				    
				    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
				    
				     CALL EXTRAPOLATE_BOUND_LINEAR(IEX,L,NGP,I,1)
				    
				    
				END DO	!COMPONENTS
				END DO	!NGP
	
	    ICONSIDERED=I
	 CALL SOLUTIONTRIAV22d(N,ICONSIDERED)
	  END DO
	  
	  
	  
	  
	  END DO
!$OMP END DO 
 
END SUBROUTINE LINEAR_SCHEME2d


SUBROUTINE condition_store(N)
!> @brief
!> Subroutine for storing the condition number for gradient approximationt tests
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX,IDFR,sols,ihgt,ihgj
REAL::solm,NDMAX,NDMIN,NFMIN,NFMAX,SOLD,solm2,SOLC
real,dimension(90)::polyfun,grad_m
real,dimension(2)::grad_1,sold2
integer::ixg
KMAXE=XMPIELRANK(N)

! call  QUADRATUREline(N,IGQRULES)

!$OMP DO SCHEDULE (STATIC) 		
	DO I=1,KMAXE
	ICONSIDERED=I
! 	NDMAX=ZERO;NDMIN=10.0E18;NFMIN=10.0E18;NFMAX=ZERO
! 	
	if (ielem(n,i)%interior.ne.1)then
		CALL ALLGRADS_INNER2d(N,I)
	ELSE
		CALL ALLGRADS_MIX2d(N,I)
	END IF
	
	VEXT(1,1)=ielem(n,i)%xxc;VEXT(1,2)=ielem(n,i)%yyc
	
	ax=VEXT(1,1)
	ay=VEXT(1,2)
	
	
	compwrt=0
  
 compwrt=0
  solm=zero;solm2=zero
   polyfun(1:ielem(n,iconsidered)%idegfree)=basis_rec2d(N,Ax,Ay,ielem(n,iconsidered)%iorder,Iconsidered,ielem(n,iconsidered)%idegfree)
 compwrt=0
   do ixg=1,ielem(n,iconsidered)%idegfree
   solm=Solm+polyfun(ixg)
   end do
   do ixg=1,ielem(n,iconsidered)%idegfree
   polyfun(ixg)=sqrt((DF2dX(AX,AY,IXG)*DF2dX(AX,AY,IXG))+(DF2dX(AX,AY,IXG)*DF2dX(AX,AY,IXG)))
   end do
  do ixg=1,ielem(n,iconsidered)%idegfree
   solm2=Solm2+polyfun(ixg)
   end do
   
	
	VEXT(1,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2),VEXT(1,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
	ax=VEXT(1,1)
	ay=VEXT(1,2)
				compwrt=0
				CONSMATRIX(1,1:IELEM(N,ICONSIDERED)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,ICONSIDERED)%IORDER,ICONSIDERED,IELEM(N,ICONSIDERED)%IDEGFREE)
				
				
							   
					DO IHGT=1,2
		DO IHGJ=1,2
		    AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
		END DO
	    END DO   
					   
					   
				
! 				    DO K=1,IELEM(N,I)%IDEGFREE
! 					    GRADM(K)=ILOCAL_RECON5(1)%GRADIENTS(1,k,1)
! 					    
! 						    XDER(K)=DF2dX(AX,AY,K);  YDER(K)=DF2dY(AX,AY,K);  
! 					    
! 					    
! 				    END DO
					    grad_1= ZERO
					    
				    DO K=1,IELEM(N,I)%IDEGFREE
					    grad_1(1) = grad_1(1) + ILOCAL_RECON5(1)%GRADIENTS(1,k,1)*DF2dX(AX,AY,K)
					    grad_1(2) = grad_1(2) + ILOCAL_RECON5(1)%GRADIENTS(1,k,1)*DF2dY(AX,AY,K);
					    

				    END DO 
				
				
				sold2(1:2)=MATMUL(AINVJT(1:2,1:2),grad_1(1:2))
				solC=sqrt((sold2(1)**2)+(sold2(2)**2))
! 				IF (INITCOND.EQ.0)THEN
! 			
				RESSOLUTION(1,1)=ZERO
				
! 	    			DO K=1,IELEM(N,I)%IDEGFREE
				GRADSSL(1:IELEM(N,I)%IDEGFREE,1)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1)
! 	    			END DO
				    RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1))
				    
! 				    sold=solm+RESSOLUTION(1,1)
 				    sold=U_C(I)%VAL(1,1)+RESSOLUTION(1,1)
				    ilocal_Recon3(i)%cond(1)=abs(sold-solm)/abs(solm)
				    ilocal_Recon3(i)%cond(2)=abs(solC-solm2)/abs(solm2)
				    !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSSL(1:IELEM(N,I)%IDEGFREE,1:1))
				    
				     
	  
	  
	  
	  
	  END DO
!$OMP END DO 
 
END SUBROUTINE condition_store





SUBROUTINE MUSCL(N)
!> @brief
!> Subroutine for MUSCL type reconstruction in 3D
 IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,IEX,IEUL,JX,LX,KMAXE,iq,LL,NGP,NND,IQP,idummy,ii,icd
REAL::UMIN,UMAX,PSITOT,ADDC,DIVG0,LIMVBG
REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T


KMAXE=XMPIELRANK(N)
call  QUADRATUREQUAD3D(N,IGQRULES);WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
call QUADRATURETRIANG(N,IGQRULES); WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)

!$OMP DO SCHEDULE (STATIC) 		
	DO II=1,NOF_INTERIOR
	I=EL_INT(II)
	ICONSIDERED=I
	   IF (IELEM(N,I)%FULL.EQ.0)THEN
	    
! 		if (ielem(n,i)%interior.ne.1)then
		CALL ALLGRADS_INNER(N,I)
		    IF (LIMITER.EQ.3)THEN !if 2
		        if (ILOCAL_RECON3(I)%LOCAL.eq.1)then  !if 1
			  DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)
			  END DO
			ELSE
			   DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)
			    else
			      UTEMP(IQ,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)
			    end if
			   end do
			END IF                   !end if 1
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,NOF_VARIABLES
			UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			END DO
			
			
		    ELSE
			   K=0
			  UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
			  K=1
			    
			  
			  
			  DO L=1,IELEM(N,I)%IFCA
			  K=K+1
			    UTEMP(K,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:NOF_VARIABLES)
			  END DO
			  
			  
			  
			  
			  
			 
			  
			  
			  			  
			  
			  
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,NOF_VARIABLES
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    END IF   !end if 2
		
		    USOL(:,:,:)=ZERO
		   
			
			      
                                        if (fastest.ne.1)then   
                                         icd=0
                                        DO L=1,IELEM(N,I)%IFCA	!faces2
                                            idummy=0
                                                IF (FASTEST_Q.EQ.1)THEN
                                                        if (ielem(n,i)%types_faces(L).eq.5)then
                                                            iqp=qp_quad
                                                        else
                                                            iqp=qp_triangle
                                                        end if
                                                ELSE
                                                    if (ielem(n,i)%types_faces(L).eq.5)then
                                                        iqp=qp_quad
                                                        NND=4
                                                            
                                                            do K=1,nnd
                                                                VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                            END DO
                                                            
                                                            
                                                        call  QUADRATUREQUAD3D(N,IGQRULES)
                                                        
                                                        
                                                        
                                                    else
                                                        iqp=QP_TRIANGLE
                                                        NND=3
                                                            
                                                            do K=1,nnd
                                                                VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                            END DO
                                                            
                                                            
                                                        call QUADRATURETRIANG(N,IGQRULES)
                                                        
                                                        
                                                    end if
                    
                                            end if
			!end fastest
			
		
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
                                    icd=icd+1
			 
				  CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
                            end do
                            end do
				      
					  GRADSS(1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)
				      
                                           
				  call gemm(                                                  &
                            CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
                            GRADSS(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
                                ICD=0
                                DO L=1,IELEM(N,I)%IFCA
                                    if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
                                    else
				    iqp=qp_triangle
                                    end if
                                     do NGP=1,iqp
			     ICD=ICD+1
	    			
				 
				  USOL(1:NOF_VARIABLES,L,Ngp)=((U_C(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(icd,1:NOF_VARIABLES)))
				  END DO	 
			    END DO
			    
			    
			    
		       else
		       DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			    
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
			
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				AZ = QPOINTS2D(3,NGP)-ielem(n,i)%zzc
		       
			 
				      DO LX=1,3
					  GRADSS(LX,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,LX,1:NOF_VARIABLES)
				      END DO
  
				  
				
				  DO IEX=1,NOF_VARIABLES		      
				  RESSOLUTION(1,1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))+(az*GRADSS(3,IEX))
				  USOL(IEX,L,ngp)=((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       
		   END DO
		   end if
			 
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				    DO iex=1,nof_Variables
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,nof_Variables
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      ILOCAL_RECON3(I)%ULEFT(:,L,:)=ZERO
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad;WEIGHT_T2(1:QP_QUAD)=WEIGHTs_Q(1:QP_QUAD)
			      else
				    iqp=QP_TRIANGLE;WEIGHT_T2(1:QP_TRIANGLE)=WEIGHTs_T(1:QP_TRIANGLE)
			      end if
			      do NGP=1,iqp
				    USOL(1:nof_Variables,l,Ngp)=USOL(1:nof_Variables,l,Ngp)-U_C(I)%VAL(1,1:nof_Variables)
				    CALL EXTRAPOLATE_BOUND_MUSCL(IEX,L,NGP,I,1)  
			      END DO
			 END DO
			 
			 !TURBULENCE NOW
			 IF (((PASSIVESCALAR.GT.0).OR.(TURBULENCE.EQ.1)).and.(icoupleturb.eq.1))then
			 IF (LIMITER.EQ.3)THEN
		        if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			  DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			  END DO
			ELSE
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			    else
			      UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
			    end if
			    end do
			END IF
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			END DO
			
			
		    ELSE
			   K=0
			  UTEMP(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			  K=1
			  DO L=1,IELEM(N,I)%IFCA
			    K=K+1
			    UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			  END DO
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    END IF
		
		    USOL(:,:,:)=ZERO
			DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			      
			   if (fastest.ne.1)then   
				      IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
		else
				    iqp=qp_triangle
		end if
			    ELSE
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
			end if
			
		
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
			 
                                  CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO

				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! 				  GRADSS(1:IELEM(N,I)%IDEGFREE,1)=GRADSS(1:IELEM(N,I)%IDEGFREE,IEX)
				  RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX))
				 
				  !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX:IEX))
				  USOL(IEX,L,Ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       else
			    
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
			
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				AZ = QPOINTS2D(3,NGP)-ielem(n,i)%zzc
		       
			 
				      DO LX=1,3
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO
  
				  
				
				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR		      
				  RESSOLUTION(1,1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))+(az*GRADSS(3,IEX))
				  USOL(IEX,L,ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       end if
		   END DO
			 
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				    DO iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR	
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR	
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				    
				    
				    USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)=USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)-U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,ngp)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))
				    
				    
			      END DO
			 END DO
			 
			 END IF
			 
			 
			 
			 
			 
			 
			 
			 
			 
			 
			 
			 
			 
			 
				 
		
		
		
		ICONSIDERED=I
	 CALL SOLUTIONTRIAV2(N,ICONSIDERED)
		
		
	END IF
	END DO
!$OMP END DO 
	
!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_bounded
	I=EL_BND(II)
	ICONSIDERED=I
	      IF (IELEM(N,I)%FULL.EQ.0)THEN
		CALL ALLGRADS_MIX(N,I)
		    IF (LIMITER.EQ.3)THEN
			  if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
				UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)
			    END DO
			  ELSE
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
				UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)
			      else
				UTEMP(IQ,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)
			      end if
			    END DO
			  END IF
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,NOF_VARIABLES
			  UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			  UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			  END DO
			
			
		    ELSE
			    K=0
			    UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
			    K=1
			    
			    
			     DO L=1,IELEM(N,I)%IFCA
				  IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
					IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
					      if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
					      K=K+1
					      UTEMP(k,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(l))%VAL(1,1:nof_variables)
					      ELSE
					      !NOT PERIODIC ONES IN MY CPU		  
					      END IF
					ELSE
					      K=K+1
						UTEMP(K,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(l))%VAL(1,1:nof_variables)		    
					END IF
				  ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			    
					      IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
							if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
							    IF (FASTEST.EQ.1)THEN
							      K=K+1
							      UTEMP(K,1:NOF_VARIABLES)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_variables)
							    ELSE
							      K=K+1
							      UTEMP(K,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_variables)
							    END IF
							END IF
					      ELSE
					      
						      IF (FASTEST.EQ.1)THEN
							K=K+1
							UTEMP(K,1:NOF_VARIABLES)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_variables)
						      ELSE
							K=K+1
							UTEMP(K,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_variables)
! 							
						      END IF
						    
					    END IF
				  END IF
			   END DO
! 				if (utemp(k,1).ne.utemp(k,1))then
! 				
! 			      
! 			      end if
			  
			UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,NOF_VARIABLES
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    
		    END IF
		    
		    USOL(:,:,:)=ZERO
			
			      
			   if (fastest.ne.1)then !1 fastest
			   icd=0
			    DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			   
                                        IF (FASTEST_Q.EQ.1)THEN   !if 1
                                                            if (ielem(n,i)%types_faces(L).eq.5)then
                                                            iqp=qp_quad
                                                            else
                                                            iqp=qp_triangle
                                                            end if
                                        ELSE
                                                
			   
                                        if ((iperiodicity.eq.1))then		!periodicity      !if 2
                                                IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                        if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
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
                                                                                END DO
                                                                                ELSE
                                                                                facex=l;
                                                                                CALL coordinates_face_PERIOD(n,iconsidered,facex)
                                                                                do K=1,nnd
                                                                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                                                END DO
                                                                                END IF
                                                                                
                                                                        call  QUADRATUREQUAD3D(N,IGQRULES)
                                                                        
                                                                        
                                                                        
                                                                        else
                                                                        iqp=QP_TRIANGLE
                                                                        NND=3
                                                                                IF (IDUMMY.EQ.0)THEN
                                                                                do K=1,nnd
                                                                                VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                                                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                                                END DO
                                                                                ELSE
                                                                                facex=l;
                                                                                CALL coordinates_face_PERIOD(n,iconsidered,facex)
                                                                                do K=1,nnd
                                                                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                                                END DO
                                                                                END IF
                                                                                
                                                                        call QUADRATURETRIANG(N,IGQRULES)
                                                                        
                                                                        
                                                                        end if
				  
                                else
			
                                            if (ielem(n,i)%types_faces(L).eq.5)then
                                                iqp=qp_quad
                                                NND=4
                                                    
                                                    do K=1,nnd
                                                        VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                                        VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                    END DO
                                                    
                                                    
                                                call  QUADRATUREQUAD3D(N,IGQRULES)
                                                
                                                
                                                
                                            else
                                                iqp=QP_TRIANGLE
                                                NND=3
                                                    
                                                    do K=1,nnd
                                                        VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                                                        VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                                                    END DO
                                                    
                                                    
                                                call QUADRATURETRIANG(N,IGQRULES)
                                                
                                                
                                            end if
       
       
                            end if      !if 2
			
		      end if      !if 1
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				icd=icd+1
				
			 
			 
				  CONSMATRIX(icd,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
			 
                            end do
                        end do
                        
                        	 
                        
                        
			 
! 				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,1:IELEM(N,I)%IDEGFREE,1:NOF_VARIABLES)
! 					  if (GRADSS(LX,1).ne.GRADSS(LX,1))then
! 					  
! 			      
! 					  end if
				      

				 call gemm(                                                  &
                            CONSMATRIX(1:ICD,1:IELEM(N,I)%IDEGFREE),                &
                            GRADSS(1:IELEM(N,I)%IDEGFREE,1:NOF_vARIABLES),                                                &
                            RESSOLUTION(1:ICD,1:NOF_vARIABLES)                                                    &
                                )
                                ICD=0
                                DO L=1,IELEM(N,I)%IFCA
                                    if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
                                    else
				    iqp=qp_triangle
                                    end if
                                     do NGP=1,iqp
			     ICD=ICD+1
	    			
				 
				  USOL(1:NOF_VARIABLES,L,Ngp)=((U_C(I)%VAL(1,1:NOF_VARIABLES)+RESSOLUTION(icd,1:NOF_VARIABLES)))
				  END DO	 
			    END DO
			    
			    
			    
		    
		       else
		       DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
		       
			    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
			else
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
			end if
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				AZ = QPOINTS2D(3,NGP)-ielem(n,i)%zzc
		       
			 
				      DO LX=1,3
					  GRADSS(LX,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,LX,1:NOF_VARIABLES)
				      END DO
  
				  
				
				  DO IEX=1,NOF_VARIABLES		      
				  RESSOLUTION(1:1,1:1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))+(az*GRADSS(3,IEX))
				  USOL(IEX,L,ngp)=((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		END DO       
		       end if
		   
		    
		    DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				    DO iex=1,nof_Variables
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,nof_Variables
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
			 
			 
			 
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      ILOCAL_RECON3(I)%ULEFT(:,L,:)=ZERO
			       if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad;WEIGHT_T2(1:QP_QUAD)=WEIGHTs_Q(1:QP_QUAD)
			      else
				    iqp=QP_TRIANGLE;WEIGHT_T2(1:QP_TRIANGLE)=WEIGHTs_T(1:QP_TRIANGLE)
			      end if
			      do NGP=1,iqp
				    USOL(1:nof_Variables,l,Ngp)=USOL(1:nof_Variables,l,Ngp)-U_C(I)%VAL(1,1:nof_Variables)
				    CALL EXTRAPOLATE_BOUND_MUSCL(IEX,L,NGP,I,1)  
			      END DO
				    
			     
			 END DO
		    
		    
		    
		    !TURBULENCE NOW
			 IF (((PASSIVESCALAR.GT.0).OR.(TURBULENCE.EQ.1)).and.(icoupleturb.eq.1))then
		    
			     IF (LIMITER.EQ.3)THEN
			  if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
				UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			    END DO
			  ELSE
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
				UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			      else
				UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
			      end if
			    end do
			  END IF
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			  UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			  UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			  END DO
			
			
		    ELSE
			    K=0
			    UTEMP(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			    K=1
			    DO L=1,IELEM(N,I)%IFCA
			    
			    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
				IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
				      if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
				      K=K+1
				      UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      ELSE
				      !NOT PERIODIC ONES IN MY CPU			  				  
				      END IF
				ELSE
					K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				END IF
			    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			      
				IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
				    if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					IF (FASTEST.EQ.1)THEN
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(L))%SOL(IELEM(N,i)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
					ELSE
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					  (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
					END IF
				    END IF
				ELSE
				
					IF (FASTEST.EQ.1)THEN
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(L))%SOL(IELEM(N,i)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
					
					ELSE
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					  (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
					END IF
				      
			      END IF
			  END IF
			  END DO
			UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    
		    END IF
		    
		    USOL(:,:,:)=ZERO
			DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			      
			   if (fastest.ne.1)then   
			   
			   
			    IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
		else
				    iqp=qp_triangle
		end if
			    ELSE
			   
			   
			   
			   if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
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
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
			else
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
			end if
			
		end if
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
			 
			 
				   CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC(N,AX,AY,AZ,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO

				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR		      
! 				   GRADSS(1:IELEM(N,I)%IDEGFREE,1)=GRADSS(1:IELEM(N,I)%IDEGFREE,IEX)
				  RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX))
				  !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX:IEX))
				  
				  USOL(IEX,L,Ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       else
			    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  
					  END IF
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  
					  END IF
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
			else
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    NND=4
					 
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREQUAD3D(N,IGQRULES)
				    
				    
				    
				  else
				    iqp=QP_TRIANGLE
				    NND=3
					  
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  
					  
				    call QUADRATURETRIANG(N,IGQRULES)
				    
				    
				  end if
       
       
			end if
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				AZ = QPOINTS2D(3,NGP)-ielem(n,i)%zzc
		       
			 
				      DO LX=1,3
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO
  
				  
				
				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR		      
				  RESSOLUTION(1:1,1:1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))+(az*GRADSS(3,IEX))
				  USOL(IEX,L,ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       end if
		   END DO
		    
		    DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				    DO iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=QP_TRIANGLE
			      end if
			      do NGP=1,iqp
				    
				    
				    USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)=USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)-U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,ngp)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))
				    
				    
			      END DO
			 END DO
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
		    
			END IF
	
	
! 	END IF
	
	  ICONSIDERED=I
	 CALL SOLUTIONTRIAV2(N,ICONSIDERED)
	  END IF
	
	
	
	END DO
!$OMP END DO



END SUBROUTINE MUSCL



SUBROUTINE ARBITRARY_ORDER3(N)
!> @brief
!> Subroutine controlling the reconstruction in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I


KMAXE=XMPIELRANK(N)
	
 SELECT CASE(IWENO)
 
 
  CASE(1)
  CALL MUSCL(N)		
  CALL WENOWEIGHTS(N)
  CALL CHECKSOL(N)


  CASE(-1)

  CALL MUSCL(N)		
  CALL CHECKSOL(N)
 
 
 
  CASE(0)
  IF (FIRSTORDER.EQ.1)THEN
  CALL PIECEWISE_CONSTANT(N)

  ELSE
  CALL LINEAR_SCHEME(N)

  END IF
 
 END SELECT
 
 
 

END SUBROUTINE ARBITRARY_ORDER3




SUBROUTINE ARBITRARY_ORDER2(N)
!> @brief
!> Subroutine controlling the reconstruction in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I


KMAXE=XMPIELRANK(N)
	
 SELECT CASE(IWENO)
 
 
  CASE(1)
  CALL MUSCL2d(N)		
  CALL WENOWEIGHTS2d(N)
  CALL CHECKSOL2d(N)


  CASE(-1)
  CALL MUSCL2d(N)		
  CALL CHECKSOL2d(N)
 
 
 
  CASE(0)
  IF (FIRSTORDER.EQ.1)THEN
  CALL PIECEWISE_CONSTANT2d(N)
  ELSE
  if (initcond.eq.0)then
  call condition_store(n)
  
  else
  CALL LINEAR_SCHEME2d(N)
  end if
  END IF
 
 END SELECT
 
 
 

END SUBROUTINE ARBITRARY_ORDER2

SUBROUTINE MUSCL2d(N)
!> @brief
!> Subroutine for MUSCL type reconstruction in 2D
 IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,IEX,IEUL,JX,LX,KMAXE,iq,LL,NGP,NND,IQP,idummy,ii
REAL::UMIN,UMAX,PSITOT,ADDC,DIVG0,LIMVBG


KMAXE=XMPIELRANK(N)
call  QUADRATURELINE(N,IGQRULES)

!$OMP DO SCHEDULE (STATIC) 		
	DO II=1,NOF_INTERIOR
	I=EL_INT(II)
	ICONSIDERED=I
	
	   IF (IELEM(N,I)%FULL.EQ.0)THEN
	    
		
		CALL ALLGRADS_INNER2d(N,I)
		    IF ((LIMITER.EQ.3).OR.(LIMITER.EQ.9))THEN
                        IF (LIMITER.EQ.3)THEN
		        if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			  DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)
			  END DO
			ELSE
			   DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)
			    else
			      UTEMP(IQ,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)
			    end if
			   end do
			END IF
! 			DO IEX=1,NOF_VARIABLES
! 			UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)=UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)-UTEMP(1,iex)
! 			END DO
			
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,NOF_VARIABLES
			UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			END DO
			ELSE
			 if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			  DO IQ=1,13
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)
			  END DO
			ELSE
			   DO IQ=1,13
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)
			    else
			      UTEMP(IQ,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)
			    end if
			   end do
			END IF
! 			DO IEX=1,NOF_VARIABLES
! 			UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)=UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)-UTEMP(1,iex)
! 			END DO
			
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,NOF_VARIABLES
			UTMIN(IEX)=MINVAL(UTEMP(1:13,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:13,IEX))
			END DO
			
			
			
			
			END IF
			
			
			
			
		    ELSE
			   K=0
			  UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
			  K=1
			  DO L=1,IELEM(N,I)%IFCA
			    K=K+1
			    UTEMP(K,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:NOF_VARIABLES)!-UTEMP(1,1:NOF_VARIABLES)
			  END DO
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,NOF_VARIABLES
			  
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    END IF
		
		    USOL(:,:,:)=ZERO
			DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			      
			   if (fastest.ne.1)then   
			   IF (FASTEST_Q.NE.1)THEN
				   
				    iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
       
       
			
			else
				    iqp=qp_line
		        end if
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
			 CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
! 				  CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(LX,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,LX,1:NOF_VARIABLES)
				      END DO

				  DO IEX=1,NOF_VARIABLES
				  
! 				   GRADSS(1:IELEM(N,I)%IDEGFREE,1)=GRADSS(1:IELEM(N,I)%IDEGFREE,IEX)
				  RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX))
				  !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX:IEX))
				  USOL(IEX,L,Ngp)=((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       else
			    
				   iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
       
       
			
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				
		       
			 
				      DO LX=1,2
					  GRADSS(LX,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,LX,1:NOF_VARIABLES)
				      END DO
  
				  
				
				  DO IEX=1,NOF_VARIABLES		      
				  RESSOLUTION(1,1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))
				  USOL(IEX,L,ngp)=((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       end if
		   END DO
			 
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      
			      do NGP=1,qp_line
				    DO iex=1,nof_Variables
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,nof_Variables
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
! 			      if (ielem(n,i)%types_faces(L).eq.5)then
! 				    iqp=qp_quad
! 			      else
! 				    iqp=QP_TRIANGLE
! 			      end if
			      do NGP=1,qp_line
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
! 			 if (itestcase.eq.4)then
			 
			 ielem(n,i)%vortex(1)=WENO(1,1)

			 
			 
! 			 end if
			 
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			     ILOCAL_RECON3(I)%ULEFT(:,L,:)=ZERO
			      do NGP=1,qp_line
				    
				    
				    USOL(1:nof_Variables,l,Ngp)=USOL(1:nof_Variables,l,Ngp)-U_C(I)%VAL(1,1:nof_Variables)
				    CALL EXTRAPOLATE_BOUND_MUSCL(IEX,L,NGP,I,1)
				    
				    
! 				    ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,ngp)=U_C(I)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,l,Ngp)*WENO(1:nof_Variables,1))
				    
				    
			      END DO
			 END DO
			 
			 !TURBULENCE NOW
			 IF (((PASSIVESCALAR.GT.0).OR.(TURBULENCE.EQ.1)).and.(icoupleturb.eq.1))then
			 IF (LIMITER.EQ.3)THEN
		        if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			  DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			  END DO
			ELSE
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			    else
			      UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
			    end if
			    end do
			END IF
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			END DO
			
			
		    ELSE
			   K=0
			  UTEMP(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			  K=1
			  DO L=1,IELEM(N,I)%IFCA
			    K=K+1
			    UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			  END DO
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    END IF
		
		    USOL(:,:,:)=ZERO
			DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			      
			   if (fastest.ne.1)then   
			   if (fastest_q.ne.1)then
				   
				    iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				else
				iqp=qp_line
				end if
				    
				  
       
       
			
			
		
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
			 
			 
				  CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO

				      
				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				  
! 				   GRADSS(1:IELEM(N,I)%IDEGFREE,1)=GRADSS(1:IELEM(N,I)%IDEGFREE,IEX)
				  RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX))
				  !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX:IEX))
				  USOL(IEX,L,Ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       else
			    
				   
				    iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
       
       
			
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				
		       
			 
				      DO LX=1,2
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO
  
				  
				
				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR		      
				  RESSOLUTION(1,1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))
				  USOL(IEX,L,ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       end if
		   END DO
			 
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      
				    iqp=qp_line
			      
			      do NGP=1,iqp
				    DO iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR	
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			 
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      iqp=qp_line
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			  ILOCAL_RECON3(I)%ULEFTTURB(:,L,:)=zero
			      iqp=qp_line
			      do NGP=1,iqp
				   
				    
				    USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)=USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)-U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				    
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,1)=ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,1)+(U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)))*wequa2d(ngp)
				    else
				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,ngp)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))
				    
				    end if
				    
			      END DO
			 END DO
			 
			 END IF
			 
			 
		 ICONSIDERED=I
	 CALL SOLUTIONTRIAV22d(N,ICONSIDERED)	 
			 
	END IF		 
		
	END DO
	!$OMP END DO
	
	
	
	!$OMP DO SCHEDULE (STATIC) 		
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I
	
	   IF (IELEM(N,I)%FULL.EQ.0)THEN
	
		CALL ALLGRADS_MIX2d(N,I)
		   IF ((LIMITER.EQ.3).OR.(LIMITER.EQ.9))THEN
                        IF (LIMITER.EQ.3)THEN
		        if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			  DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)
			  END DO
			ELSE
			   DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)
			    else
			      UTEMP(IQ,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)
			    end if
			   end do
			END IF
! 			DO IEX=1,NOF_VARIABLES
! 			UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)=UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)-UTEMP(1,iex)
! 			END DO
			
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,NOF_VARIABLES
			UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			END DO
			ELSE
			 if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			  DO IQ=1,13
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:NOF_VARIABLES)
			  END DO
			ELSE
			   DO IQ=1,13
			    IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
			      UTEMP(IQ,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:nof_variables)
			    else
			      UTEMP(IQ,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),1:nof_variables)
			    end if
			   end do
			END IF
! 			DO IEX=1,NOF_VARIABLES
! 			UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)=UTEMP(2:IELEM(N,I)%iNUMNEIGHBOURS,IEX)-UTEMP(1,iex)
! 			END DO
			
			UTMIN=ZERO;UTMAX=ZERO
			DO IEX=1,NOF_VARIABLES
			UTMIN(IEX)=MINVAL(UTEMP(1:13,IEX))
			UTMAX(IEX)=MAXVAL(UTEMP(1:13,IEX))
			END DO
			
			
			
			
			END IF
			
			
			
			
		    ELSE
			    K=0
			    UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
			    K=1
			    DO L=1,IELEM(N,I)%IFCA
			    
			    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
				IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
				      if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
				      K=K+1
				      UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				      ELSE
				      !NOT PERIODIC ONES IN MY CPU			  				  
				      END IF
				ELSE
					K=K+1
					UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				END IF
			    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			      
				IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
				    if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					IF (FASTEST.EQ.1)THEN
					
					  K=K+1
					UTEMP(K,1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(L))%SOL(IELEM(N,i)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
					ELSE
					
					  K=K+1
					UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					  (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
					END IF
				    END IF
				ELSE
				
					IF (FASTEST.EQ.1)THEN
					
					  K=K+1
					UTEMP(K,1:nof_variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(L))%SOL(IELEM(N,i)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
					
					ELSE
					
					  K=K+1
					UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					  (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
					END IF
				      
			      END IF
			  END IF
			  END DO
			  
			  
! 			  DO IEX=1,NOF_VARIABLES
! 			UTEMP(2:k,IEX)=UTEMP(2:k,IEX)-UTEMP(1,iex)
! 			END DO
			  
			  
			UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,NOF_VARIABLES
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    
		    END IF
		    
		    USOL(:,:,:)=ZERO
			DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			      
			   if (fastest.ne.1)then   
			   if (fastest_q.ne.1)then
			   if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 iqp=qp_line
				    
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
			else
				   
				    iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
       
       
			end if
			else
			 iqp=qp_line
			
			
			end if
			
		
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
			 
			 
				  CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
			 
				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(LX,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,LX,1:NOF_VARIABLES)
				      END DO

				  DO IEX=1,NOF_VARIABLES		      
! 				   GRADSS(1:IELEM(N,I)%IDEGFREE,1)=GRADSS(1:IELEM(N,I)%IDEGFREE,IEX)
				  RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX))
				  !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX:IEX))
				  
				  USOL(IEX,L,Ngp)=((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       else
			    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
! 				 if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
			else
				   
				     iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
       
       
			end if
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				
		       
			 
				      DO LX=1,2
					  GRADSS(LX,1:NOF_VARIABLES)=ILOCAL_RECON5(1)%GRADIENTS(1,LX,1:NOF_VARIABLES)
				      END DO
  
				  
				
				  DO IEX=1,NOF_VARIABLES		      
				  RESSOLUTION(1:1,1:1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))
				  USOL(IEX,L,ngp)=((U_C(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       end if
		   END DO
		    
		    DO L=1,IELEM(N,I)%IFCA	!faces2
			      iqp=qp_line
			      do NGP=1,iqp
				    DO iex=1,nof_Variables
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,nof_Variables
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      iqp=qp_line
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
! 			  if (itestcase.eq.4)then
			 
			 ielem(n,i)%vortex(1)=WENO(1,1)

			 
! 			 end if
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			  ILOCAL_RECON3(I)%ULEFT(:,L,:)=zero
			 iqp=qp_line
			      do NGP=1,iqp
				    
				    
				    USOL(1:nof_Variables,l,Ngp)=USOL(1:nof_Variables,l,Ngp)-U_C(I)%VAL(1,1:nof_Variables)
				    
				    CALL EXTRAPOLATE_BOUND_MUSCL(IEX,L,NGP,I,1)
! 				    ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,ngp)=U_C(I)%VAL(1,1:nof_Variables)+(USOL(1:nof_Variables,l,Ngp)*WENO(1:nof_Variables,1))
				    
				    
				 
			      END DO
			 END DO
		    
		    
		    
		    !TURBULENCE NOW
			 IF (((PASSIVESCALAR.GT.0).OR.(TURBULENCE.EQ.1)).and.(icoupleturb.eq.1))then
		    
			     IF (LIMITER.EQ.3)THEN
			  if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
				UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			    END DO
			  ELSE
			    DO IQ=1,IELEM(N,I)%iNUMNEIGHBOURS
			      IF (ILOCAL_RECON3(I)%IHEXB(1,IQ).EQ.N)THEN
				UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(ILOCAL_RECON3(I)%IHEXL(1,IQ))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			      else
				UTEMP(IQ,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IQ))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IQ),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
			      end if
			    end do
			  END IF
			  UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			  UTMIN(IEX)=MINVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			  UTMAX(IEX)=MAXVAL(UTEMP(1:IELEM(N,I)%iNUMNEIGHBOURS,IEX))
			  END DO
			
			
		    ELSE
			    K=0
			    UTEMP(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_Ct(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
			    K=1
			    DO L=1,IELEM(N,I)%IFCA
			    
			    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
				IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
				      if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
				      K=K+1
				      UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      ELSE
				      !NOT PERIODIC ONES IN MY CPU			  				  
				      END IF
				ELSE
					K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				END IF
			    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
			      
				IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
				    if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					IF (FASTEST.EQ.1)THEN
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(L))%SOL(IELEM(N,i)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
					ELSE
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					  (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
					END IF
				    END IF
				ELSE
				
					IF (FASTEST.EQ.1)THEN
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(L))%SOL(IELEM(N,i)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
					
					ELSE
					
					  K=K+1
					UTEMP(K,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					  (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
					END IF
				      
			      END IF
			  END IF
			  END DO
			UTMIN=ZERO;UTMAX=ZERO
			  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
			    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
			  END DO
		    
		    END IF
		    
		    USOL(:,:,:)=ZERO
			DO L=1,IELEM(N,I)%IFCA	!faces2
			      idummy=0
			      
			   if (fastest.ne.1)then 
			   IF (FASTEST_Q.NE.1)THEN
			   if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
			else
				   iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:dims)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:dims)-ILOCAL_RECON3(I)%VEXT_REF(1:dims))
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
       
       
			end if
			else
			
			iqp=qp_line
			end if
			
		
		
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
			 
			 
				  CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE)=BASIS_REC2d(N,AX,AY,IELEM(N,I)%IORDER,I,IELEM(N,I)%IDEGFREE)
				      DO LX=1,IELEM(N,I)%IDEGFREE
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO
				  
				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR		      
! 				   GRADSS(1:IELEM(N,I)%IDEGFREE,1)=GRADSS(1:IELEM(N,I)%IDEGFREE,IEX)
				  RESSOLUTION(1,1) = DOT(CONSMATRIX(1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX))
				  !RESSOLUTION(1:1,1:1)=MATMUL(CONSMATRIX(1:1,1:IELEM(N,I)%IDEGFREE),GRADSS(1:IELEM(N,I)%IDEGFREE,IEX:IEX))
				  
				  USOL(IEX,L,Ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       else
			    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				  iqp=qp_line
				  nnd=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
			else
				   
				    iqp=qp_line
				    NND=2
					 
					  do K=1,nnd
					    VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    
					  END DO
					 
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
       
       
			end if
			
		
		
			     do NGP=1,iqp
				AX = QPOINTS2D(1,NGP)-ielem(n,i)%xxc
				AY = QPOINTS2D(2,NGP)-ielem(n,i)%yyc
				
		       
			 
				      DO LX=1,2
					  GRADSS(LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=ILOCAL_RECON5(1)%GRADIENTS2(1,LX,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				      END DO
  
				  
				
				  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR		      
				  RESSOLUTION(1:1,1:1)=(ax*GRADSS(1,IEX))+(ay*GRADSS(2,IEX))
				  USOL(IEX,L,ngp)=((U_CT(I)%VAL(1,IEX)+RESSOLUTION(1,1)))
				  END DO	 
			    END DO
		       
		       		       
		       end if
		   END DO
		    
		    DO L=1,IELEM(N,I)%IFCA	!faces2
			      
				    iqp=qp_line
			      
			      do NGP=1,iqp
				    DO iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				 
				 ICONS_E=IEX
				 FACEX=L
				 ICONS_S=NGP
				 ICONSIDERED=I  
				CALL SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)    
				    end do
			      END DO
			 END DO
			 
			 do iex=1,TURBULENCEEQUATIONS+PASSIVESCALAR
			 limvbg=tolbig
			 DO L=1,IELEM(N,I)%IFCA	!faces2
			      iqp=qp_line
			      do NGP=1,iqp
				  	  
				  
				  LIMVBG=MIN(LIMVBG,PSI(iex,L,Ngp) )
			      end do
			 end do
			 WENO(IEX,1)=LIMVBG
			 end do
			 
! 			 DO L=1,IELEM(N,I)%IFCA	!faces2
! 			       iqp=qp_line
! 				    do NGP=1,iqp
! 				     
! 				    USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)=USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)-U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,ngp)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))
! 				      
! 				   
! 				    
! 				    
! 			      END DO
! 			 END DO
		    
			DO L=1,IELEM(N,I)%IFCA	!faces2
			  ILOCAL_RECON3(I)%ULEFTTURB(:,L,:)=zero
			      iqp=qp_line
			      do NGP=1,iqp
				   
				    
				    USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)=USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)-U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
				    
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,1)=ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,1)+(U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)))*wequa2d(ngp)
				    else
				    ILOCAL_RECON3(I)%ULEFTTURB(1:TURBULENCEEQUATIONS+PASSIVESCALAR,L,ngp)=U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)+(USOL(1:TURBULENCEEQUATIONS+PASSIVESCALAR,l,Ngp)*WENO(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))
				    
				    end if
				    
			      END DO
			 END DO
		    
		    
		    
		    
		    
		    
		    
		    
		    
			END IF
	
	
	
	  ICONSIDERED=I
	 CALL SOLUTIONTRIAV22d(N,ICONSIDERED)
	
	
	
	END IF
	END DO
!$OMP END DO



END SUBROUTINE MUSCL2d











SUBROUTINE CHECKSOL(N)
!> @brief
!> Subroutine for checking the reconstructed solution
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,NGP,iqp,iex
INTEGER::REDUCE1,kmaxe

KMAXE=XMPIELRANK(N)


IF (ITESTCASE.GE.3)THEN

!$OMP DO SCHEDULE (STATIC)	
	DO I=1,KMAXE	!ALL ELEMENTS
		REDUCE1=0
		  DO L=1,IELEM(N,I)%IFCA	!faces2
			      if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad_n
			      else
				    iqp=QP_TRIANGLE_n
			      end if
				  do NGP=1,iqp
					      
						LEFTV(1:NOF_VARIABLES)=ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)
						RIGHTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
						CALL CONS2PRIM2(N)
						
						
				
					    IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(0.9*RIGHTV(1))).OR.((ABS(LEFTV(5)-RIGHTV(5))).GE.(0.9*RIGHTV(5)))) THEN
						    REDUCE1=1
					    END IF
				
				
					
				  END DO
		END DO	
		
		
		if (ielem(n,i)%hybrid.eq.1)then
		reduce1=1
		end if
		
		IF (REDUCE1.EQ.1)THEN
		do iex=1,NOF_VARIABLES
		ILOCAL_RECON3(I)%ULEFT(iex,:,:)=u_c(i)%val(1,iex)
		end do
		end if			
		
		
		if (turbulence.eq.1)then
		if (icoupleturb.eq.1)then
		REDUCE1=0
		DO L=1,IELEM(N,I)%IFCA	!faces2
			      
				   if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad
				    else
				    iqp=QP_TRIANGLE
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
		end if
		
		IF (REDUCE1.EQ.1)THEN
		do iex=1,1
		ILOCAL_RECON3(I)%ULEFTTURB(1,:,:)=u_ct(i)%val(1,1)
		end do
		end if	
		end if
		end if
		
		
		
		
		
		
	end do
!$OMP END DO		
		
! ELSE


! 
! 
! !$OMP DO SCHEDULE (STATIC)	
! 	DO I=1,KMAXE	!ALL ELEMENTS
! ! 		REDUCE1=0
! 		  DO L=1,IELEM(N,I)%IFCA	!faces2
! 			      if (ielem(n,i)%types_faces(L).eq.5)then
! 				    iqp=qp_quad
! 			      else
! 				    iqp=QP_TRIANGLE
! 			      end if
! 				  do NGP=1,iqp
! 					      
! 						LEFTV(1:NOF_VARIABLES)=ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)
! 						RIGHTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
! ! 						CALL CONS2PRIM2(N)
! 						
! 						
! 				
! 					    IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(0.7*RIGHTV(1)))) THEN
! ! 						    REDUCE1=1
! 						    ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)=u_c(i)%val(1,1:NOF_VARIABLES)
! 					    END IF
! 				
! 				
! 					
! 				  END DO
! 		END DO			
! ! 		IF (REDUCE1.EQ.1)THEN
! ! 		do iex=1,NOF_VARIABLES
! ! 		ILOCAL_RECON3(I)%ULEFT(iex,:,:)=u_c(i)%val(1,iex)
! ! 		end do
! ! 		end if			
! 	end do
! !$OMP END DO		
		
END IF		
		
		

END SUBROUTINE CHECKSOL
 

 
 
SUBROUTINE CHECKSOL2D(N)
!> @brief
!> Subroutine for checking the reconstructed solution in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,NGP,iqp,iex,kmaxe
INTEGER::REDUCE1
KMAXE=XMPIELRANK(N)


IF (ITESTCASE.GE.3)THEN

!$OMP DO SCHEDULE (STATIC)	
	DO I=1,KMAXE	!ALL ELEMENTS
		REDUCE1=0
		 DO L=1,IELEM(N,I)%IFCA	!faces2
			      
				    iqp=qp_LINE_n
			      
				  do NGP=1,iqp
					      
						LEFTV(1:NOF_VARIABLES)=ILOCAL_RECON3(I)%ULEFT(1:NOF_VARIABLES,L,NGP)
						RIGHTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
						CALL CONS2PRIM2d2(N)
						
						
				
					    IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(0.8*RIGHTV(1))).OR.((ABS(LEFTV(4)-RIGHTV(4))).GE.(0.8*RIGHTV(4)))) THEN
						    REDUCE1=1
					    END IF
				
				
					
				  END DO
		END DO		
		if (ielem(n,i)%hybrid.eq.1)then
		reduce1=1
		end if
		IF (INITCOND.EQ.30)THEN
		 IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
                    DO L=1,IELEM(N,I)%IFCA
                                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.NE.5)then	!PERIODIC IN MY CPU
                                            REDUCE1=1
                                        END IF
                                        end if
                                    END IF
                    END DO
                 END IF
                 END IF
		IF (REDUCE1.EQ.1)THEN

		do iex=1,NOF_VARIABLES
		ILOCAL_RECON3(I)%ULEFT(iex,:,:)=u_c(i)%val(1,iex)
		end do
		end if	
		
		if (turbulence.eq.1)then
		if (icoupleturb.eq.1)then
		REDUCE1=0
		DO L=1,IELEM(N,I)%IFCA	!faces2
			      
				    iqp=qp_LINE_n
			      
				  do NGP=1,iqp
				    leftv(1)=ILOCAL_RECON3(I)%ULEFTTURB(1,L,ngp)
				    RIGHTV(1)=U_Ct(I)%VAL(1,1)
					      IF (((ABS(LEFTV(1)-RIGHTV(1))).GE.(0.8*RIGHTV(1))).OR.(leftv(1).le.zero)) THEN
						    REDUCE1=1
					    END IF
				  end do
		end do
		
		if (ielem(n,i)%hybrid.eq.1)then
		reduce1=1
		end if
		
		IF (REDUCE1.EQ.1)THEN
		do iex=1,1
		ILOCAL_RECON3(I)%ULEFTTURB(1,:,:)=u_ct(i)%val(1,1)
		end do
		
		do iex=1,NOF_VARIABLES
		ILOCAL_RECON3(I)%ULEFT(iex,:,:)=u_c(i)%val(1,iex)
		end do
		
		
		end if
		
		
		
		
		end if
		end if
		
	end do
!$OMP END DO		
		
END IF		
		
		

END SUBROUTINE CHECKSOL2D 

SUBROUTINE SOLUTIONTRIAV22d(N,ICONSIDERED)
!> @brief
!> Subroutine for extrapolating the unlimited reconstructed values for diffusive fluxes in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I,J,K,L,M,PPP,IEUL,IEX,IHGT,IHGJ,KMAXE,DECOMF,ICNN,IQDR,NVAR,idummy,iqp,nnd,ngp
REAL::RAA1,RAA2,PAA1,PAA2
REAL::SOLX




KMAXE=XMPIELRANK(N)
GRADSS=ZERO
			
			
 I=ICONSIDERED
			
      


call  QUADRATUREline(N,IGQRULES)

 
 if (itestcase.eq.4)then
 
 ILOCAL_RECON3(I)%ULEFTV(:,:,:,:)=zero;
 
 IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
 ILOCAL_RECON3(I)%ULEFTTURBV(:,:,:,:)=zero;
ILOCAL_RECON3(I)%ULEFTTURB(:,:,:)=zero;
 END IF
 
		
 
 
 
    if (fastest.ne.1)then
!                    AINVJT(1:2,1:2)=ILOCAL_RECON3(I)%INVCTJAC(1:2,1:2)
	    DO IHGT=1,2
		DO IHGJ=1,2
		    AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
		END DO
	    END DO
    END IF

	DO l=1,IELEM(N,I)%IFCA	!faces
	      IDUMMY=0
		    if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
				  IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
					  if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
					      IDUMMY=1
					  END IF
				  END IF
		
				
				 
				    iqp=qp_line
				    NND=2
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD2d(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  END IF
					  
				    call  QUADRATUREline(N,IGQRULES)
				    
				    
				    
				  
		else
				    iqp=qp_line
				    NND=2
					  
					  do K=1,nnd
					    VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
					  END DO
					  
					  
				    call  QUADRATUREline(N,IGQRULES)
       
       
		end if
			     do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP)
				AY = QPOINTS2D(2,NGP)
! 				AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1)
				AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2)
! 				AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				
				
				
				END IF
				
				
				SELECT CASE(IELEM(N,I)%GGS)
				
				CASE(0)
				
				!TURBULENCE FIRST
				IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
				
				if (icoupleturb.eq.0)then
				    DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,1)=u_ct(i)%val(1,nvar)
				    
				    else
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
				    
				    end if
				    
				    END DO
				
				end if
				
				
				DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    DO K=1,IELEM(N,I)%IDEGFREE
					    GRADTEM(K)=ILOCAL_RECON5(1)%GRADIENTSTURB(1,K,NVAR)
					    
						    XDER(K)=DF2dX(AX,AY,K);  YDER(K)=DF2dY(AX,AY,K);  
					    
					    
				    END DO
					    UGRADLOC = ZERO
					    
				    DO K=1,IELEM(N,I)%IDEGFREE
					    UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
					    UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
					    

				    END DO 
				    
				     if (reduce_comp.eq.1)then
				    
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,1) =ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,1)+( MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2)))*wequa2d(ngp)
				    else
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,NGP) = MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))
				    
				    end if
				    END DO
				END IF
				
				!MEAN FLOW GRADIENTS
				DO K=1,IELEM(N,I)%IDEGFREE
					GRADTEM(K)=ILOCAL_RECON5(1)%GRADIENTSTEMP(K)
					  XDER(K)=DF2dX(AX,AY,K);  YDER(K)=DF2dY(AX,AY,K);  
			         END DO
					UGRADLOC = ZERO
					  DO K=1,IDEGFREE
					  UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
					  UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
					  
					  END DO 
					  if (reduce_comp.eq.1)then
					   ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,1) = ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,1)+MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))*wequa2d(ngp)
					  
					  else
					  ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP) = MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))
					  end if
				
				  DO IEX=1,2
				  DO K=1,IELEM(N,I)%IDEGFREE
					GRADTEM(K)=ILOCAL_RECON5(1)%VELOCITYDOF(IEX,K)
					  XDER(K)=DF2dX(AX,AY,K);  YDER(K)=DF2dY(AX,AY,K);  
				  END DO
					  UGRADLOC = ZERO
					  DO K=1,IDEGFREE
					  UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
					  UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
					  
					  END DO 
					  if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,1) = ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,1)+MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))*wequa2d(ngp)
					  else
					  ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,NGP) = MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))
					  end if
				  END DO	    
				  				
				CASE(1)
				
				
				  !TURBULENCE FIRST
				IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
				DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,1)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,1)+ILOCAL_RECON3(I)%GRADs(3+NVAR,1:2)*wequa2d(ngp)
				    else
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,NGP)=ILOCAL_RECON3(I)%GRADs(3+NVAR,1:2)
				    end if
				    
				    END DO
				    
				if (icoupleturb.eq.0)then
				    DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,1)=u_ct(i)%val(1,nvar)
				    else
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
				    end if
				    END DO
				
				end if
				END IF
				
				!MEAN FLOW GRADIENTS
				if (reduce_comp.eq.1)then
				
				ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,1)= ILOCAL_RECON3(I)%GRADs(3,1:2)
				else
				ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)= ILOCAL_RECON3(I)%GRADs(3,1:2)
				
				end if
					  
				DO IEX=1,2
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,1) = ILOCAL_RECON3(I)%GRADs(IEX,1:2)
				    else
				    ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,NGP) = ILOCAL_RECON3(I)%GRADs(IEX,1:2)
				    end if
				END DO
				
				CASE(2)
				  !TURBULENCE FIRST
				IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
				
				if (icoupleturb.eq.0)then
				    DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,1)=u_ct(i)%val(1,nvar)
				    else
				     ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
				    end if
				    END DO
				
				end if
				
				
				      DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
					  if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,1)=ILOCAL_RECON3(I)%GRADs(3+NVAR,1:2)
					  else
					   ILOCAL_RECON3(I)%ULEFTTURBV(1:2,NVAR,L,NGP)=ILOCAL_RECON3(I)%GRADs(3+NVAR,1:2)
					  end if
					  END DO
				END IF
				      DO K=1,IELEM(N,I)%IDEGFREE
					      GRADTEM(K)=ILOCAL_RECON5(1)%GRADIENTSTEMP(K)
						XDER(K)=DF2dX(AX,AY,K);  YDER(K)=DF2dY(AX,AY,K);  
				      END DO
					UGRADLOC = ZERO
					  DO K=1,IDEGFREE
					  UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
					  UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
					  
					  END DO 
					  if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,1) =  ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,1)+MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))*wequa2d(ngp)
					  else
					  ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP) = MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))
					  end if
				
					DO IEX=1,2
						DO K=1,IELEM(N,I)%IDEGFREE
						      GRADTEM(K)=ILOCAL_RECON5(1)%VELOCITYDOF(IEX,K)
							XDER(K)=DF2dX(AX,AY,K);  YDER(K)=DF2dY(AX,AY,K);  
						END DO
						UGRADLOC = ZERO
						DO K=1,IDEGFREE
						UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
						UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
						
						END DO 
						if (reduce_comp.eq.1)then
						ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,1) = ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,1)+MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))*wequa2d(ngp)
						else
						ILOCAL_RECON3(I)%ULEFTV(1:2,IEX+1,L,NGP) = MATMUL(AINVJT(1:2,1:2),UGRADLOC(1:2))
						end if
					END DO	    
					
					
				 END SELECT
			END DO
	  END DO
	  
	  END IF

					    
END SUBROUTINE SOLUTIONTRIAV22d


SUBROUTINE SOLUTIONTRIAV2(N,ICONSIDERED)
!> @brief
!> Subroutine for extrapolating the unlimited reconstructed values for diffusive fluxes in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
INTEGER::I,J,K,L,M,PPP,IEUL,IEX,IHGT,IHGJ,KMAXE,DECOMF,ICNN,IQDR,NVAR,idummy,iqp,nnd,ngp,icd
REAL::RAA1,RAA2,PAA1,PAA2
REAL::SOLX
REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T


call  QUADRATUREQUAD3D(N,IGQRULES);WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
call QUADRATURETRIANG(N,IGQRULES); WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)

KMAXE=XMPIELRANK(N);GRADSS=ZERO;I=ICONSIDERED
			
if (itestcase.eq.4)then
ILOCAL_RECON3(I)%ULEFTV(:,:,:,:)=zero;
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
ILOCAL_RECON3(I)%ULEFTTURBV(:,:,:,:)=zero;ILOCAL_RECON3(I)%ULEFTTURB(:,:,:)=zero;
END IF
if (fastest.ne.1)then
	    DO IHGT=1,3;DO IHGJ=1,3
		AINVJT(IHGT,IHGJ)=ILOCAL_RECON3(I)%INVCCJAC(IHGJ,IHGT)
	    END DO;END DO
END IF

	DO l=1,IELEM(N,I)%IFCA;IDUMMY=0
	  IF (FASTEST_Q.EQ.1)THEN
				  if (ielem(n,i)%types_faces(L).eq.5)then
				    iqp=qp_quad;WEIGHT_T2(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
				    iqp=qp_triangle;WEIGHT_T2(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
	    ELSE
	           if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
			IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
			if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5)then	!PERIODIC IN OTHER CPU
			IDUMMY=1
			END IF
			END IF
		        if (ielem(n,i)%types_faces(L).eq.5)then
			      iqp=qp_quad;NND=4
			      IF (IDUMMY.EQ.0)THEN
			      do K=1,nnd
			      VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
			      VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
			      END DO
			      ELSE
			      facex=l;
			      CALL coordinates_face_PERIOD(n,iconsidered,facex)
			      do K=1,nnd
			      VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
			      END DO
			      END IF
			      call  QUADRATUREQUAD3D(N,IGQRULES)
			else
			    iqp=QP_TRIANGLE;NND=3
					  IF (IDUMMY.EQ.0)THEN
					  do K=1,nnd
					    VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
					    VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  ELSE
					  facex=l;
					  CALL coordinates_face_PERIOD(n,iconsidered,facex)
					  do K=1,nnd
					  VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
					  END DO
					  END IF
					call QUADRATURETRIANG(N,IGQRULES)
			end if
		else
			if (ielem(n,i)%types_faces(L).eq.5)then
			iqp=qp_quad;NND=4
			do K=1,nnd
			VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
			VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
			END DO
			call  QUADRATUREQUAD3D(N,IGQRULES)
			else
			iqp=QP_TRIANGLE
			NND=3
			do K=1,nnd
			VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
			VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
			END DO
			call QUADRATURETRIANG(N,IGQRULES)
			end if
			end if
		end if
		 ICD=0
			 do NGP=1,iqp			!for gqp
				IF (FASTEST_Q.NE.1)THEN
				AX = QPOINTS2D(1,NGP);AY = QPOINTS2D(2,NGP);AZ = QPOINTS2D(3,NGP)
				ELSE
				AX = ILOCAL_RECON3(I)%QPOINTS(L,NGP,1);AY = ILOCAL_RECON3(I)%QPOINTS(L,NGP,2);AZ = ILOCAL_RECON3(I)%QPOINTS(L,NGP,3)
				END IF
				icd=icd+1
			DO K=1,IELEM(N,I)%IDEGFREE
                                IF (POLY.EQ.1) THEN
					XXDER(K,ICD)=DFX(AX,AY,AZ,K);  YYDER(K,ICD)=DFY(AX,AY,AZ,K);  ZZDER(K,ICD)=DFZ(AX,AY,AZ,K)
				 END IF
				IF (POLY.EQ.2) THEN
					XXDER(K,ICD)=DLX(AX,AY,AZ,K);  YYDER(K,ICD)=DLY(AX,AY,AZ,K);  ZZDER(K,ICD)=DLZ(AX,AY,AZ,K)
				END IF
				END DO
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
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,1)=u_ct(i)%val(1,nvar)
				    ELSE
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
				    END IF
				    END DO
				
				end if
				
				
				DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
! 				    DO K=1,IELEM(N,I)%IDEGFREE
					    GRADTEM(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%GRADIENTSTURB(1:IELEM(N,I)%IDEGFREE,K,NVAR)
				   
					    UGRADLOC = ZERO
					     UGRADLOC(1)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),XXDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    UGRADLOC(2)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),YYDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    UGRADLOC(3)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),ZZDER(1:IELEM(N,I)%IDEGFREE,ICD))
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,1) = ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,1)+MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))*WEIGHT_T2(NGP)
				    ELSE
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,NGP) = MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))
				    
				    END IF
                                END DO
				END IF
				

					GRADTEM(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%GRADIENTSTEMP(1:IELEM(N,I)%IDEGFREE)
! 				
					UGRADLOC = ZERO
					     UGRADLOC(1)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),XXDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    UGRADLOC(2)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),YYDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    UGRADLOC(3)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),ZZDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,1) =ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,1)+ MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))*WEIGHT_T2(NGP)
					  ELSE
					  ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP) = MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))
					  END IF
				
				  DO IEX=1,3
! 		
					GRADTEM(1:IELEM(N,I)%IDEGFREE)=ILOCAL_RECON5(1)%VELOCITYDOF(IEX,1:IELEM(N,I)%IDEGFREE)
! 		
					 UGRADLOC = ZERO
					     UGRADLOC(1)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),XXDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    UGRADLOC(2)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),YYDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    UGRADLOC(3)=DOT(GRADTEM(1:IELEM(N,I)%IDEGFREE),ZZDER(1:IELEM(N,I)%IDEGFREE,ICD))
					    if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,1) = ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,1)+MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))*WEIGHT_T2(NGP)
					  ELSE
					   ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,NGP) = MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))
					  END IF
! 		
				  END DO	    
				  	
				CASE(1)
				
				
				  !TURBULENCE FIRST
				IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
				
				
				if (icoupleturb.eq.0)then
				    DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				   ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,1)=u_ct(i)%val(1,nvar)
				   ELSE
				   ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
				   END IF
				    END DO
				
				end if
				
				
				DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,1)=ILOCAL_RECON3(I)%GRADs(4+NVAR,1:3)
				    ELSE
				    ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,NGP)=ILOCAL_RECON3(I)%GRADs(4+NVAR,1:3)
				    END IF
				    END DO
				END IF
				
				!MEAN FLOW GRADIENTS
				if (reduce_comp.eq.1)then
				ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,1) = ILOCAL_RECON3(I)%GRADs(4,1:3)
					  
				DO IEX=1,3
				    ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,1) = ILOCAL_RECON3(I)%GRADs(IEX,1:3)
				END DO
				ELSE
				ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP) = ILOCAL_RECON3(I)%GRADs(4,1:3)
					  
				DO IEX=1,3
				    ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,NGP) = ILOCAL_RECON3(I)%GRADs(IEX,1:3)
				END DO
				
				
				END IF
				
				CASE(2)
				  !TURBULENCE FIRST
				IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
				
				if (icoupleturb.eq.0)then
				    DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
				    if (reduce_comp.eq.1)then
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,1)=u_ct(i)%val(1,nvar)
				    ELSE
				    
				    ILOCAL_RECON3(I)%ULEFTTURB(NVAR,L,NGP)=u_ct(i)%val(1,nvar)
				    END IF
				    END DO
				
				end if
				
				
				
				
				      DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
					    if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,1)=ILOCAL_RECON3(I)%GRADs(4+NVAR,1:3)
					    ELSE
					    ILOCAL_RECON3(I)%ULEFTTURBV(1:3,NVAR,L,NGP)=ILOCAL_RECON3(I)%GRADs(4+NVAR,1:3)  
					    END IF
					  END DO
				END IF
				      DO K=1,IELEM(N,I)%IDEGFREE
					      GRADTEM(K)=ILOCAL_RECON5(1)%GRADIENTSTEMP(K)
						IF (POLY.EQ.1) THEN
						      XDER(K)=DFX(AX,AY,AZ,K);  YDER(K)=DFY(AX,AY,AZ,K);  ZDER(K)=DFZ(AX,AY,AZ,K)
						END IF
						IF (POLY.EQ.2) THEN
						  XDER(K)=DLX(AX,AY,AZ,K);  YDER(K)=DLY(AX,AY,AZ,K);  ZDER(K)=DLZ(AX,AY,AZ,K)
						END IF
				      END DO
					UGRADLOC = ZERO
					  DO K=1,IDEGFREE
					  UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
					  UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
					  UGRADLOC(3) = UGRADLOC(3) + GRADTEM(K)*ZDER(K)
					  END DO 
					      if (reduce_comp.eq.1)then
					  ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,1) =ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,1)+ MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))*WEIGHT_T2(NGP)
					  ELSE
					  
					    ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP) = MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))
					  END IF
				
					DO IEX=1,3
						DO K=1,IELEM(N,I)%IDEGFREE
						      GRADTEM(K)=ILOCAL_RECON5(1)%VELOCITYDOF(IEX,K)
							IF (POLY.EQ.1) THEN
							      XDER(K)=DFX(AX,AY,AZ,K);  YDER(K)=DFY(AX,AY,AZ,K);  ZDER(K)=DFZ(AX,AY,AZ,K)
							END IF
							IF (POLY.EQ.2) THEN
							  XDER(K)=DLX(AX,AY,AZ,K);  YDER(K)=DLY(AX,AY,AZ,K);  ZDER(K)=DLZ(AX,AY,AZ,K)
							END IF
						END DO
						UGRADLOC = ZERO
						DO K=1,IDEGFREE
						UGRADLOC(1) = UGRADLOC(1) + GRADTEM(K)*XDER(K)
						UGRADLOC(2) = UGRADLOC(2) + GRADTEM(K)*YDER(K)
						UGRADLOC(3) = UGRADLOC(3) + GRADTEM(K)*ZDER(K)
						END DO 
						 if (reduce_comp.eq.1)then
						ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,1) = ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,1) +MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))*WEIGHT_T2(NGP)
						ELSE
						ILOCAL_RECON3(I)%ULEFTV(1:3,IEX+1,L,NGP) = MATMUL(AINVJT(1:3,1:3),UGRADLOC(1:3))
						
						END IF
					END DO	    
					
					
				 END SELECT
			END DO
	  END DO
	  
	  END IF
 
      
      



END SUBROUTINE SOLUTIONTRIAV2


SUBROUTINE SLOPE_LIMITERS(N,ICONSIDERED,ICONS_E,FACEX,ICONS_S)
!> @brief
!> Subroutine MUSCL slope limiters
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONS_E,FACEX,ICONS_S,ICONSIDERED
INTEGER::I,L,IEX,NGP
REAL::D2,SFD,EPSI2,DMIN,DPLUS,KAPPA_VEN,psi2,sig_1,delu,y_fun,s_y,pol_MOG
IEX=ICONS_E
L=FACEX
I=ICONSIDERED
NGP=ICONS_S
KAPPA_VEN=10.0D0
psi2=zero

					  D2=USOL(IEX,L,NGP)-UTEMP(1,IEX)


					  IF (ABS(D2).LT.tolsmall)THEN

					  PSI(iex,L,ngp)=1.0D0

					  ELSE IF (D2.GE.tolsmall)THEN
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
					      dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
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
					      dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMAX(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2))
					      
					      
					       CASE(8)								!VENKATAKRISHNAN
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					      dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
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
					       dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
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
					       dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
					      DPLUS=UTMIN(IEX)-UTEMP(1,IEX)
					      EPSI2=(KAPPA_VEN*IELEM(N,I)%MINEDGE)**3
					      PSI(iex,L,ngp)=(1.0D0/(DMIN))*(((((DPLUS**2)+EPSI2)*DMIN)+(2.0D0*(DMIN**2)*DPLUS))/((DPLUS**2)+(2.0D0*DMIN**2)+(DMIN*DPLUS)+EPSI2)) 
					      CASE(8)								!VENKATAKRISHNAN
					      DMIN=USOL(iex,L,ngp)-UTEMP(1,IEX)
					       dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
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

! ! ---------------------------------------------------------------------------------------------!

! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE RECON
