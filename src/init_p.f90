MODULE INITIALISATION
USE DECLARATION
USE LIBRARY
USE TRANSFORM
!USE FLUXES
USE PROFILE
IMPLICIT NONE

 CONTAINS
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED TO INITIALISE THE FLOW FIELD SOLUTION!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!FOR TESTING PURPOSES ONLY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIALISE(N)
 !> @brief
!> This subroutine calls the initialisation of the computational domain
IMPLICIT NONE
integer,INTENT(IN)::N
REAL,allocatable,DIMENSION(:)::RG,ARG
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
INTEGER:: prev_turbequation,INITIAL,III,i,k,jx,QQP,INC,kmaxe,jkn,ki,iterr,JX2,ind1,KX



if (rungekutta.eq.4)then
ind1=7
else
ind1=5
end if


IF (LAMPS.EQ.1)THEN
III=1
ELSE
III=0
END IF


    


prev_turbequation=0
if (prev_turbmodel.eq.1) then
prev_turbequation=1
end if 
if (prev_turbmodel.eq.2) then
prev_turbequation=2
end if 
  IF (LAMPS.EQ.1)THEN
 
 allocate(rg(5+prev_turbequation+1))
  ELSE
 allocate(rg(5+prev_turbequation))
END IF
 

VECCOS(:)=0.0
KMAXE=XMPIELRANK(N)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




IF (RESTART.EQ.0)THEN

 
!$OMP PARALLEL DEFAULT(SHARED) private(iNITIAL,i,QQP,INC,jx,jx2,k)
	IF (ISCHEME.LT.2)THEN
		IF (ITESTCASE.EQ.0)THEN
!$OMP DO SCHEDULE (STATIC) 
		DO INITIAL=1,KMAXE
			U_C(INITIAL)%VAL(1,1)=1.0D0
			U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
		END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.EQ.1)THEN
!$OMP DO SCHEDULE (STATIC) 
		DO INITIAL=1,KMAXE
			pox(1)=IELEM(N,INITIAL)%XXC
			poy(1)=IELEM(N,INITIAL)%yyC
			poz(1)=IELEM(N,INITIAL)%zzC
			U_C(INITIAL)%VAL(1,1)=LINEAR_INIT3D(N)
			U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
		END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.EQ.2)THEN
!$OMP DO SCHEDULE (STATIC) 
		DO INITIAL=1,KMAXE
			pox(1)=IELEM(N,INITIAL)%XXC
			poy(1)=IELEM(N,INITIAL)%yyC
			poz(1)=IELEM(N,INITIAL)%zzC
			U_C(INITIAL)%VAL(1,1)=LINEAR_INIT3D(N)
			U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
		END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.GE.3)THEN
!$OMP DO SCHEDULE (STATIC) 
		DO INITIAL=1,KMAXE
			VECCOS(:)=ZERO
			pox(1)=IELEM(N,INITIAL)%XXC
			poy(1)=IELEM(N,INITIAL)%yyC
			poz(1)=IELEM(N,INITIAL)%zzC
			CALL INITIALISE_EULER3D(N)
			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(INITIAL)%VAL(1,1:nof_Variables)=VECCOS(1:nof_Variables)
			U_CT(INITIAL)%VAL(1,1:0+turbulenceequations+passivescalar)=VECCOS(6:5+turbulenceequations+passivescalar)
			else
			U_C(INITIAL)%VAL(1,:)=VECCOS(:)
			if (itestcase.eq.3)U_E(INITIAL)%VAL(1,:)=U_C(INITIAL)%VAL(1,:)
			end if
		END DO
!$OMP END DO
		END IF
	ELSE
!$OMP DO SCHEDULE (STATIC) 
	DO I=1,KMAXE
		
		iconsidered=i
		 
		 VEXT=ZERO
    NODES_LIST=ZERO
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=ZERO
        
      jx=IELEM(N,I)%NONODES
!       
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,:)=inoder(JX2)%CORD(:)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3
    
      SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
      if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREHEXA(N,IGQRULES)
      VOLTEMP=HEXAVOLUME(N)/IELEM(N,I)%totvolume
      QQP=QP_HEXA
      
	  DO INC=1,QQP
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			IF (ITESTCASE.LE.2)THEN
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N)*WEQUA3D(INC)*(VOLTEMP)
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			ELSE
			CALL INITIALISE_EULER3D(N)
			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
			U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
			VECCOS(6:5+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
			else
			U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
			
			if (itestcase.eq.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
			end if
			END IF
	  END DO
      
      
      else
      COUNT_1=0
       do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	  
	  CALL QUADRATUREtetra(N,IGQRULES)
                        if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N)
                        else
                        VOLTEMP=TETRAVOLUME(N)/IELEM(N,I)%totvolume
                        end if
	  
	 
	  
	    QQP=QP_TETRA
	  DO INC=1,QQP
	  COUNT_1=COUNT_1+1
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			IF (ITESTCASE.LE.2)THEN
			
			IF (DG.EQ.1)THEN
			basis_vector(1)=1.0d0
                compwrt=-2
             BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
              compwrt=0
                             
                             tempsol(1,:)=LINEAR_INIT3D(N)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
			
			
			
			ELSE
			
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N)*WEQUA3D(INC)*(VOLTEMP)
			
			END IF
			
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			IF (DG.EQ.1)U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
			ELSE
			
			
			
			CALL INITIALISE_EULER3D(N)
                if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
                U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
                U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
                VECCOS(6:5+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
                else
                
                
                    IF (DG.EQ.1)THEN
                    basis_vector(1)=1.0d0
                    compwrt=-2
                    BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                    compwrt=0
                                    DO KX=1,nof_Variables
                                    tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                                
                                    U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                                        END DO
                    ELSE
                    U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                    if (itestcase.eq.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                    END IF
                end if
			END IF
	  END DO
	  
	  
	  
	  
	  
       END DO
       end if
	
	! 	

      CASE(2)
      
       
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	
       CALL QUADRATUREtetra(N,IGQRULES)
       
       
       
        
                        if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N)
                        COUNT_1=0
                        else
                        VOLTEMP=TETRAVOLUME(N)/IELEM(N,I)%totvolume
                        end if
! 	  
	   QQP=QP_TETRA
	  DO INC=1,QQP
	  COUNT_1=COUNT_1+1
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			IF (ITESTCASE.LE.2)THEN
			
                
			
			IF (dg.EQ.1)THEN
			basis_vector(1)=1.0d0
			compwrt=-2
			
			
			
			
             BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                             compwrt=0
                             tempsol(1,:)=LINEAR_INIT3D(N)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
!                
                
                         
                              U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
			
			
               
			
			
			
			
			ELSE
			
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N)*WEQUA3D(INC)*(VOLTEMP)
			
			END IF
			
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			IF (DG.EQ.1)U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
			ELSE
			CALL INITIALISE_EULER3D(N)
			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
			U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
			VECCOS(6:5+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
			else
                    IF (dg.EQ.1)THEN
                    basis_vector(1)=1.0d0
                    compwrt=-2
                    BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                    compwrt=0
                                    DO KX=1,nof_Variables
                                    tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                                
                                    U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                                        END DO
                    ELSE
                    U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                    if (itestcase.eq.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                    END IF
			end if
			END IF
	  END DO
      CASE(3)
      
        
      COUNT_1=0
      do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	  
	  CALL QUADRATUREtetra(N,IGQRULES)
	  if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N)
                        else
                        VOLTEMP=TETRAVOLUME(N)/IELEM(N,I)%totvolume
                        end if
	     QQP=QP_TETRA
	  DO INC=1,QQP
	   COUNT_1=COUNT_1+1
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			IF (ITESTCASE.LE.2)THEN
			IF (DG.EQ.1)THEN
			basis_vector(1)=1.0d0
			compwrt=-2
             BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
            compwrt=0
                             tempsol(1,:)=LINEAR_INIT3D(N)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
			
			
			
			ELSE
			
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N)*WEQUA3D(INC)*(VOLTEMP)
			
			END IF
			
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			IF (DG.EQ.1)U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
			ELSE
			CALL INITIALISE_EULER3D(N)
			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
			U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
			VECCOS(6:5+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
			else
			IF (DG.EQ.1)THEN
                    basis_vector(1)=1.0d0
                    compwrt=-2
                    BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                    compwrt=0
                                    DO KX=1,nof_Variables
                                    tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                                
                                    U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                                        END DO
                    ELSE
                    U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                    if (itestcase.eq.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                    END IF
			end if
			END IF
	  END DO
       END DO
	
	
	
	
       CASE(4)
       
       if (IELEM(N,I)%MODE.eq.0)then
       CALL QUADRATUREPRISM(N,IGQRULES)
      VOLTEMP=PRISMVOLUME(N)/IELEM(N,I)%totvolume
       QQP=QP_PRISM
	  DO INC=1,QQP
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			IF (ITESTCASE.LE.2)THEN
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N)*WEQUA3D(INC)*(VOLTEMP)
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			ELSE
			CALL INITIALISE_EULER3D(N)
			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
			U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
			VECCOS(6:5+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
			else
			IF (DG.EQ.1)THEN
                    basis_vector(1)=1.0d0
                    compwrt=-2
                    BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                    compwrt=0
                                    DO KX=1,nof_Variables
                                    tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                                
                                    U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                                        END DO
                    ELSE
                    U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                    if (itestcase.eq.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                    END IF
			end if
			END IF
	  END DO
      else
      COUNT_1=0
       do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	  
	  CALL QUADRATUREtetra(N,IGQRULES)
	  if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N)
                        else
                        VOLTEMP=TETRAVOLUME(N)/IELEM(N,I)%totvolume
                        end if
	    QQP=QP_TETRA
	  DO INC=1,QQP
	  COUNT_1=COUNT_1+1
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			IF (ITESTCASE.LE.2)THEN
			IF (DG.EQ.1)THEN
			basis_vector(1)=1.0d0
			compwrt=-2
             BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                          compwrt=0   
                             tempsol(1,:)=LINEAR_INIT3D(N)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
			
			
			
			ELSE
			
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N)*WEQUA3D(INC)*(VOLTEMP)
			
			END IF
			
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			IF (DG.EQ.1)U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
			ELSE
			CALL INITIALISE_EULER3D(N)
			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
			U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
			VECCOS(6:5+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
			else
			IF (DG.EQ.1)THEN
                    basis_vector(1)=1.0d0
                    compwrt=-2
                    BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE) 
                    compwrt=0
                                    DO KX=1,nof_Variables
                                    tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                                
                                    U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                                        END DO
                    ELSE
                    U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                    if (itestcase.eq.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                    END IF
			end if
			END IF
	  END DO
       END DO
       end if
       


      END SELECT
   
    		
		
		
		
		
	END DO
!$OMP END DO
	END IF
!$OMP END PARALLEL 
	RES_TIME=ZERO
END IF




! IF (RESTART.GT.0)THEN     !IF_RESTART
! 	if ( RUNGEKUTTA .ge. 5 ) Then         !LEV0
!       ! 	
! 	  RESTFILE='RESTART.dat'
! 	  OPEN(1083,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 
! 
! !Just to assure they are zero if not read	
!          INITIALRES=ZERO
! 	
!        
! 
! 
!        if (IRES_UNSTEADY .lt.1) then
! 
! 	    
! 	    READ (1083)rescounter
! 	   read (1083)INITIALRES(1)
! 	  read (1083)INITIALRES(2)
! 	  read (1083)INITIALRES(3)
! 	  read (1083)INITIALRES(4)
! 	  read (1083)INITIALRES(5)
! 	  
! 
! 	  
! 	  !IMPORTANT: The restart depends on the PREVIOUS simulation
! 	  !Important on 22/6/2013
! 	  if (prev_turbmodel.eq.1)then
! 	  read (1083)INITIALRES(6)
! 	  end if
! 	   if (prev_turbmodel.eq.2)then
! 	  read (1083)INITIALRES(6)
! 	  read (1083)INITIALRES(7)
! 	  end if
! 	
! 
!    
!    else  !It mean lamy.gt.0, coz it was Transient before
!    READ(1083)   !We read and drop: it, restime
!    
!    rescounter=0
!    end if
! 	
! 
!                         
!        
! 		
! 	    DO I=1,IMAXE
! 			jkn=i
! 
! 	
! 			IF (TURBULENCE.EQ.1)THEN		!IF TURBULENCE
! 			
! 			IF (IRES_TURB.LT.1)THEN			!NON TURBULENT
! 			 READ (1083)RG(1:5+iii)
! 
! 		  
! 			  IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(1,1:5)=RG(1:5)
! 				  IF (TURBULENCEMODEL.EQ.1)THEN
! 				  !RG(6)=VISC*TURBINIT
! 				  U_CT(KI)%VAL(1,1)=VISC*TURBINIT
!                                   END IF
! 				  IF (TURBULENCEMODEL.EQ.2)THEN
! 
! 				  U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				   U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 
! 				
! 				  END IF
! 				  IF (LAMPS.LT.1)THEN
! 				  IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=ZERO
! 				  END IF
! 				  
! 				  ELSE
! 				  
! 				   IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=RG(5+1)  !LAMX is lower than 1
! 				  END IF
! 				  END IF
! 			   
! 			  END IF
! 
! 	!-----------------------------------------------------		  
! 			  !End of the corrections on 21/6/2013
! !-------------------------------------------------------			  
! 		!Modified on the 22/6/2013	
! 			
! 			ELSE   !For LAMX.gt.1
! 			 READ (1083)RG(1:5+prev_turbequation+III)
! 	
!                            
! 			  IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(1,1:5)=RG(1:5)
! 				if (Turbulence .eq. 1) then
! 				   if (turbulencemodel .eq.prev_turbmodel)then
! 				   U_CT(KI)%VAL(1,1:Turbulenceequations)=RG(6:5+Turbulenceequations)
! 
!      			   end if
! 				   
! 				   if ((turbulencemodel.eq.1).and.(prev_turbmodel.eq.2))then
! 				   U_CT(KI)%VAL(1,1) =alpha_starinf* RG(6)/RG(7)!Nu turbulent
! 				   end if
! 				   
! 				   if ((turbulencemodel.eq.2).and.(prev_turbmodel.eq.1))then
! 				 U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				  U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 				 
! 				 !U_CT(KI)%VAL(1,1)=(1.5*(I_turb_inlet*sqrt(RG(2)**2+RG(3)**2+RG(4)**2)/RG(1))**2)*RG(1)
! 				 !U_CT(KI)%VAL(1,2)=(alpha_starinf*max(U_CT(KI)%VAL(1,1), 1.5*RG(1)*(I_turb_inlet*ufreestream)**2)/&
! 				  !                      max(RG(6),INIT_mu_RATIO*visc)*RG(1))  !Already includes Rho for rho*om
! 				   end if
! 
!                       
! 				   IF (LAMPS.LT.1)THEN
! 				  IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=0.0
! 				  END IF
! 				  
! 				  ELSE
! 				 
! 				   IF (PASSIVESCALAR.GT.0)THEN
! 			
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=RG(5+prev_turbequation+1)!*RG(1)
! 				  END IF
! 				  END IF
! 			    end if
! 				 
! 			  END IF
! 
! 
! 
! 			 END IF
! 			 
! 			ELSE ! turbulence is zero
! 
! 			READ (1083)RG(1:5+prev_turbequation+III)
! 		
! 			IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(1,1:5)=RG(1:5)
! 			if (passivescalar .gt.0) then
! 
! 
! 			  IF (LAMPS.GE.1)THEN
! 			  U_CT(KI)%VAL(1,turbulenceequations+1)=RG(5+prev_turbequation+1)
! 			  ELSE
! 			  U_CT(KI)%VAL(1,turbulenceequations+1)=ZERO
! 			  END IF
! 
! 
! 
! 			end if
! 			  
! 			END IF
! 
! 			
! 
! 			END IF	!End if for turbulence
! 
! 			END DO
! 			
! 	
! 	
! 	
! 	
! 	CLOSE(1083)
! ! 	
! 	
! !----------------------------------
! 	else   !This else means the now is RUNGE KUTTA for transient       !LEV0
! 
! ! 	
! 	!!! NO LOCAL TIME STEP
! ! 	
! 	RESTFILE='RESTART.dat'
! 	OPEN(1083,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	
! 	IF ((IRES_UNSTEADY.LT.1))THEN 
! 	READ (1083) !resiter
!          iterr=0
!         res_time=ZERO
!         ELSE
!         READ (1083)ITERR,RES_TIME
!         if (initcond.eq.95)then
!    read(1083)taylor
!    end if
!         END IF
! 	
! 	IF ((IRES_UNSTEADY.LT.1))THEN       !LEV0
! 	
! 	DO K=1,5+prev_turbequation
! 	READ (1083)
! 	  
! 
! 
! 
! 	END DO
! 	  IF (LAMPS.GT.1.0)THEN
! 	  READ (1083)
! 
! 	  END IF
! 
! 
! 	END IF
! 	
! 	
! 	DO I=1,IMAXE
! 
! 			jkn=i
! 			IF ((TURBULENCE.EQ.1))THEN     !LEV1
! 			IF (IRES_TURB.LT.1)THEN           !LEV2
! 			READ (1083)RG(1:5+iii)
! 
! 			IF (XMPIE(jkn).EQ.N)THEN        !LEV3
! 			      KI=XMPIL(jkn)
! 			      U_C(KI)%VAL(1,1:5)=RG(1:5)
! ! 				
! 				
! 				if (TURBULENCEMODEL.eq.1) then
! 				U_CT(KI)%VAL(1,1)=VISC*TURBINIT!/RG(1)
! 				end if
! 
! 				if (TURBULENCEMODEL .eq. 2) then
! 				  U_CT(KI)%VAL(1,1)=1.5D0*(I_turb_inlet*ufreestream)**2
! 				  U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25D0))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 				  
! 				end if
! 				
! 
!                             if (PASSIVESCALAR .gt. 0) then
!                               if (LAMPS .gt. 0) then
!                               U_CT(KI)%VAL(1,turbulenceequations+1)=RG(6)
! 			   else
!                               U_CT(KI)%VAL(1,turbulenceequations+1)=ZERO
!                             end if
! 			   end if
! 
! 		
! 
! 				
! 			END IF                          !LEV3
! 			ELSE  !So LAMX is gt 1          !LEV2
! 			READ (1083)RG(1:5+prev_turbequation+III )
! 			IF (XMPIE(jkn).EQ.N)THEN        !LEV3
! 			      KI=XMPIL(jkn)
! 			      U_C(KI)%VAL(1,1:5)=RG(1:5)
! 				!BUG IN THE ORIGINAL FILE (No read of RG(6), so it always pass zero)
! 				!U_CT(KI)%VAL(1,1:+Passivescalar)=RG(6:Turbulenceequations+Passivescalar)
! 		     	   if (turbulencemodel .eq.prev_turbmodel)then
! 				   U_CT(KI)%VAL(1,1:Turbulenceequations)=RG(6:5+Turbulenceequations)
!      			   end if
! 				   
! 				   if ((turbulencemodel.eq.1).and.(prev_turbmodel.eq.2))then
! 				   U_CT(KI)%VAL(1,1) =VISC*TURBINIT
! 				   end if
! 				   
! 				   if ((turbulencemodel.eq.2).and.(prev_turbmodel.eq.1))then
! 				   U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				   U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 				  ! U_CT(KI)%VAL(1,1)=(1.5*(I_turb_inlet*sqrt(RG(2)**2+RG(3)**2+RG(4)**2)/RG(1))**2)*RG(1)
! 				  ! U_CT(KI)%VAL(1,2)=(alpha_starinf*U_CT(KI)%VAL(1,1)/RG(6))*RG(1)
! 				   end if
! 
! 			    IF (PASSIVESCALAR.GT.0)THEN
! 			    IF (LAMPS.GT.1)THEN
! 			    
! 
! 
! 			    U_CT(KI)%VAL(1,turbulenceequations+1)=RG(5+prev_turbequation+1)
! 
! 
! 
! 
! 			    ELSE
! 			    U_CT(KI)%VAL(1,turbulenceequations+1)=ZERO
! 
! 
! 
! 
! 			    END IF
! 			    END IF
! 				 
! 			END IF                  !LEV3
! 
! 
! 			
! 			END IF                  !LEV2
! 			ELSE   !So turbulence is 0    !LEV1
! 			READ (1083)RG(1:5+prev_turbequation+iii)
! 			IF (XMPIE(jkn).EQ.N)THEN  !LEV2
! 			      KI=XMPIL(jkn)
! 			      U_C(KI)%VAL(1,1:5)=RG(1:5)
! 				
! 				 
! 			IF (PASSIVESCALAR.GT.0)THEN
! 			    IF (LAMPS.GT.1.0)THEN
! 			    
! 			    U_CT(KI)%VAL(1,1)=RG(5+prev_turbequation+1)
! 
! 
! 			    ELSE
! 			    U_CT(KI)%VAL(1,1)=ZERO
! 
! 
! 
! 
! 			    END IF
! 			    END IF  
! 			END IF       
! 
! 			END IF        !LEV2
! 			
! 
! 			    
! 	END DO
! 	CLOSE(1083)
! 
!     end if    !LEV1
! END IF    !LEV0


! if (Averaging .EQ. 1) then
! 
! if (average_restart .eq. 1) then
! 	    IF (LAMPS.GE.1)THEN
! 	     allocate(Arg(6+prev_turbequation+7))
! 
! 	    END IF
! 	    IF (LAMPS.LT.1)THEN
! 	     allocate(Arg(5+prev_turbequation+6))
! 
! 	    END IF
! !if ( RUNGEKUTTA .ge. 5 ) Then
!       ! 	
! 	  RESTFILE='RESTARTav.dat'
! 	  OPEN(1084,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	  
! 	    
! 		
! 	    DO I=1,IMAXE
! 			READ(1084)JKN
! 			IF (LAMPS.GE.1)THEN
! 			READ(1084)ARG(1:6+prev_turbequation+6+MIN(PASSIVESCALAR,1))
! 
! 	    END IF
! 	    IF (LAMPS.LT.1)THEN
! 	     
! 	    READ(1084)ARG(1:5+prev_turbequation+6)
! 	    END IF
! 
! 			
! 			IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(ind1,1:5)=ARG(1:5)
! 				  
! 			      
! 			  if (Turbulence .eq. 1) then
! 				   if (turbulencemodel .eq.prev_turbmodel)then
! 				   U_CT(KI)%VAL(ind1,1:Turbulenceequations)=ARG(6:6+Turbulenceequations)
! 				    end if
! 				   
! 				   if ((turbulencemodel.eq.1).and.(prev_turbmodel.eq.2))then
! 				   U_CT(KI)%VAL(ind1,1) =VISC*TURBINIT!/ARG(1)
! 				   end if
! 				   
! 				   if ((turbulencemodel.eq.2).and.(prev_turbmodel.eq.1))then
! 				   U_CT(KI)%VAL(ind1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				   U_CT(KI)%VAL(ind1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(ind1,1))&
! 			/L_TURB_INLET*RG(1)
! 				   !U_CT(KI)%VAL(5,2)=(alpha_starinf*U_CT(KI)%VAL(1,1)/ARG(6))
! 				   end if
! 			  			
! 			ELSE ! turbulence is zero
! 			
! 			       U_C(KI)%VAL(5,1:5)=ARG(1:5)
! 				 
! 			END IF
! 			
! 			
! 			IF (PASSIVESCALAR .EQ. 1) then
!     
! 			IF (LAMPS.GE.1)THEN
! 
! 			U_CT(KI)%VAL(ind1,turbulenceequations+1)=ARG(5+prev_turbequation+1)
! 
! 
! 			ELSE
! 			U_CT(KI)%VAL(ind1,turbulenceequations+1)=ZERO
! 			END IF
! 			END IF
! 			
! 			U_C(KI)%RMS(1:6+III)=ARG(5+PREV_TURBEQUATION+III+1:11+PREV_TURBEQUATION+III+III)
! 
! 			
! 		
! 
! 
! 
! 	
! 	      END IF
! 	
! 
! 	    END DO
! 	
! 	CLOSE(1084)
! !END if
! deallocate(Arg)
! else !it means that average_restart=0
!  DO I=1,kmaxe
! U_C(i)%VAL(ind1,:)=ZERO
! U_C(i)%RMS(:)=ZERO
! if ((passivescalar.gt.0).or.(turbulence.eq.1))then
! 
! U_CT(i)%VAL(ind1,:)=ZERO
! end if   
! end do
! 
! end if   !Level_Averaging
! 
! 
! 
! end if   !Level_RESTART=1/0

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 
! deallocate(POX)
! deallocate(POy)
! deallocate(POz)
deallocate(rg)
DEALLOCATE(VECCOS)


END SUBROUTINE INITIALISE


!INITIALISATION WITH BASIS

SUBROUTINE INITIALISE2d(N)
 !> @brief
!> This subroutine calls the initialisation of the computational domain in 2D
IMPLICIT NONE
integer,INTENT(IN)::N
REAL,allocatable,DIMENSION(:)::RG,ARG
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
INTEGER:: prev_turbequation,INITIAL,III,i,k,jx,QQP,INC,kmaxe,jkn,ki,iterr,JX2,IDX,KX


IF (LAMPS.EQ.1)THEN
III=1
ELSE
III=0
END IF



prev_turbequation=0
if (prev_turbmodel.eq.1) then
prev_turbequation=1
end if 
if (prev_turbmodel.eq.2) then
prev_turbequation=2
end if 
  IF (LAMPS.EQ.1)THEN
 
 allocate(rg(4+prev_turbequation+1))
  ELSE
 allocate(rg(4+prev_turbequation))
END IF
 

VECCOS(:)=0.0
KMAXE=XMPIELRANK(N)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)



IF (RESTART.EQ.0)THEN
!$OMP PARALLEL DEFAULT(SHARED) private(iNITIAL,i,QQP,INC)
	IF (ISCHEME.LT.2)THEN
		IF (ITESTCASE.EQ.0)THEN
!$OMP DO SCHEDULE (STATIC) 
            DO INITIAL=1,KMAXE
                U_C(INITIAL)%VAL(1,1)=1.0D0
                U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
            END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.EQ.1)THEN
!$OMP DO SCHEDULE (STATIC) 
            DO INITIAL=1,KMAXE
                pox(1)=IELEM(N,INITIAL)%XXC
                poy(1)=IELEM(N,INITIAL)%yyC
                
                U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
                U_C(INITIAL)%VAL(1,1)=LINEAR_INIT2D(N)
            END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.EQ.2)THEN
!$OMP DO SCHEDULE (STATIC) 
            DO INITIAL=1,KMAXE
                pox(1)=IELEM(N,INITIAL)%XXC
                poy(1)=IELEM(N,INITIAL)%yyC
                
                U_C(INITIAL)%VAL(1,1)=LINEAR_INIT2D(N)
                U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
            END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.GE.3)THEN
!$OMP DO SCHEDULE (STATIC) 
            DO INITIAL=1,KMAXE
                VECCOS(:)=ZERO
                pox(1)=IELEM(N,INITIAL)%XXC
                poy(1)=IELEM(N,INITIAL)%yyC
            
                CALL INITIALISE_EULER2D(N)
                if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
                U_C(INITIAL)%VAL(1,1:nof_Variables)=VECCOS(1:nof_Variables)
                U_CT(INITIAL)%VAL(1,1:0+turbulenceequations+passivescalar)=VECCOS(5:4+turbulenceequations+passivescalar)
                else
                U_C(INITIAL)%VAL(1,:)=VECCOS(:)
                end if
            END DO
!$OMP END DO
		END IF
	ELSE ! Original if statement: IF (ISCHEME.LT.2)THEN
!$OMP DO SCHEDULE (STATIC)
        DO I=1,KMAXE
            VEXT=ZERO
        !     NODES_LIST=ZERO
        !     ELTYPE=IELEM(N,I)%ISHAPE
        !     ELEM_LISTD=ZERO
        !         
        !       jx=IELEM(N,I)%NONODES
        ! 
        ! 	  do K=1,jx
        ! 	    JX2=IELEM(N,I)%NODES(k)
        ! 	    NODES_LIST(k,:)=inoder(JX2)%CORD(:)
        ! 	    VEXT(K,:)=NODES_LIST(k,:)
        ! 	  END DO
            
            ELTYPE=IELEM(N,I)%ISHAPE
            ELEM_DEC=IELEM(N,I)%VDEC
            DO K=1,IELEM(N,I)%NONODES
                NODES_LIST(k,1:2)=INODER(IELEM(N,I)%NODES(K))%CORD(1:2)
                VEXT(k,1:2)=NODES_LIST(k,1:2)
            END DO
            CALL DECOMPOSE2
    

            
            
            iconsidered=i
            
            SELECT CASE(ielem(n,i)%ishape)

            CASE(5)
            
            
                IF (IELEM(N,I)%MODE.EQ.0)THEN
                    CALL QUADRATUREQUAD(N,IGQRULES)
                    
                    VOLTEMP=1.0d0
                    QQP=QP_quad
                    
                    DO INC=1,QQP
                        
                        POX(1)=QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                        POY(1)=QPOINTS(2,INC) 
                    
                            IF (ITESTCASE.LE.2)THEN
                            U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT2D(N)*WEQUA3D(INC)*(VOLTEMP)
                            U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
                            ELSE
                            CALL INITIALISE_EULER2D(N)
                            
                            if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
                            U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
                            U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
                            VECCOS(5:4+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
                            else
                            U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                            if (itestcase.eq.3)then
                            U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                            end if
                            end if
                            END IF
                    END DO
                    
                    IF (DG.EQ.1)THEN
!                     
                        
                         DO INC=1,QQP
                            tempsol(1,:)=LINEAR_INIT2D(N)*WEQUA3D(INC)*(VOLTEMP)
                                U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                         END DO 
                           U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
                           ! WRITE(200+N,*) "ELEMENT", I,"QUAD NOT DECOMPOSED"
                     ENDIF
                ELSE
                ! this is where the initialisation will be performed for every decomposed element.
                ! INITIALISE COUNTER OF GAUSSIAN QUADRATURE POINT
                    
                    COUNT_1=0
                    DO K=1,ELEM_DEC
                        VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
                    
                        CALL QUADRATUREtriangle(N,IGQRULES)
                        
                        if (dg.eq.1)then
                        VOLTEMP=TRIANGLEVOLUME(N)
                        else
                        VOLTEMP=TRIANGLEVOLUME(N)/IELEM(N,I)%totvolume
                        end if
                        QQP=QP_Triangle
                        
                        IF (DG.EQ.1)THEN
                        
                        QQP=QP_Triangle
                        
                            DO INC=1,QQP
                                COUNT_1=COUNT_1+1
                               POX(1) = QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                               POY(1) = QPOINTS(2,INC)
                               
                                basis_vector(1)=1.0d0
                                compwrt=-2
                             BASIS_VECTOR(2:idegfree+1) = BASIS_REC2D(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),IORDER,I,IDEGFREE) 
                                compwrt=0
                             IF (ITESTCASE.GE.3)THEN
                             CALL INITIALISE_EULER2D(N)
                             END IF
                             
                                IF (ITESTCASE.LE.2)THEN
                              tempsol(1,:)=LINEAR_INIT2D(N)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                              END IF
                              
                              IF (ITESTCASE.GE.3)THEN
                              
                              DO KX=1,nof_Variables
                             tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                              
                             
                             
                             END DO
                              
                              END IF
                              
                            END DO 
!                             U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
                                !WRITE(200+N,*) "ELEMENT", I,"QUAD DECOMPOSED"
                        ENDIF
                         
                            DO INC=1,QQP                        
                                POX(1)=QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                                POY(1)=QPOINTS(2,INC)
                        
                                IF (ITESTCASE.LE.2)THEN ! Linear advection
                                    U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT2D(N)*WEQUA3D(INC)*(VOLTEMP)  !numerical
                                    U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)  !exact
                                    
                                ELSE
                                    CALL INITIALISE_EULER2D(N)
                                    
                                    IF ((turbulence .eq. 1).or.(passivescalar.gt.0)) THEN
                                        U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
                                        U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
                                        VECCOS(5:4+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
                                    ELSE
                                        U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                                        
                                        IF (itestcase.eq.3)THEN
                                            U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                                        END IF
                                    END IF
                                END IF
                            END DO !INC=1,QQP  
                    END DO !K=1,ELEM_DEC
                    
                    
                    
                END IF

            CASE(6)
                CALL QUADRATURETRIANGLE(N,IGQRULES)
                
                
                
                if (dg.eq.1)then
                        VOLTEMP=TRIANGLEVOLUME(N)
                        else
                        VOLTEMP=1.0d0
                        end if
                QQP=QP_Triangle
                            
                IF (DG.EQ.1)THEN
                
                        DO INC=1,QQP
                            POX(1) = QPOINTS(1,INC)  !POX,POY required for LINEAR_INIT2D
                            POY(1) = QPOINTS(2,INC) 
                            !QP_ARRAY(I)%X(INC) = QPOINTS(1,INC) - IELEM(N,I)%XXC !POX,POY required for LINEAR_INIT2D
                            !QP_ARRAY(I)%Y(INC) = QPOINTS(2,INC) - IELEM(N,I)%YYC
                           ! QP_ARRAY(I,INC)%QP_WEIGHT = WEQUA3D(INC)
                           
                           
                           
                            
                               
                           basis_vector(1)=1.0d0
                            compwrt=-2
                             BASIS_VECTOR(2:idegfree+1) = BASIS_REC2D(N,QP_ARRAY(I)%X(INC),QP_ARRAY(I)%Y(INC),IORDER,I,IDEGFREE) 
                             compwrt=0
                             
                              IF (ITESTCASE.GE.3)THEN
                             CALL INITIALISE_EULER2D(N)
                             END IF
                             
                                IF (ITESTCASE.LE.2)THEN
                              tempsol(1,:)=LINEAR_INIT2D(N)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                              END IF
                              
                              IF (ITESTCASE.GE.3)THEN
                              
                              DO KX=1,nof_Variables
                             tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
                         
                              U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                             
                             END DO
                              
                              END IF
                              
                            END DO 
                             
                             
                             
!                         END DO 
                           
                           
                                                !U_C(I)%VALDG(1,1,1)
                           
                            !WRITE(200+N,*) "ELEMENT", I,"TRIANGLE"
                ENDIF
                
                    DO INC=1,QQP
                        POX(1)=QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                        POY(1)=QPOINTS(2,INC)
                        
                        IF (ITESTCASE.LE.2)THEN
                            U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT2D(N)*WEQUA3D(INC)*(VOLTEMP)
                            U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
                        ELSE
                            CALL INITIALISE_EULER2D(N)
                            IF ((turbulence .eq. 1).or.(passivescalar.gt.0)) THEN
                                U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
                                U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
                                VECCOS(5:4+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
                            ELSE
                                U_C(I)%VAL(1,:)=U_C(I)%VAL(1,:)+VECCOS(:)*WEQUA3D(INC)*(VOLTEMP)
                                IF (itestcase.eq.3) THEN
                                    U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)
                                END IF
                            END IF
                        END IF
                    END DO
                !END IF
            END SELECT
        END DO
!$OMP END DO
	END IF
!$OMP END PARALLEL 
	RES_TIME=ZERO
END IF



! IF (RESTART.GT.0)THEN     !IF_RESTART
! 	if ( RUNGEKUTTA .ge. 5 ) Then         !LEV0
!       ! 	
! 	  RESTFILE='RESTART.dat'
! 	  OPEN(1083,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 
! 
! !Just to assure they are zero if not read	
!          INITIALRES=ZERO
! 	
!       
! 
! 
!        if (IRES_UNSTEADY .lt.1) then
! 
! 	    
! 	    READ (1083)rescounter
! 	   read (1083)INITIALRES(1)
! 	  read (1083)INITIALRES(2)
! 	  read (1083)INITIALRES(3)
! 	  read (1083)INITIALRES(4)
! 	  
! 	  
! 
! 	  
! 	  !IMPORTANT: The restart depends on the PREVIOUS simulation
! 	  !Important on 22/6/2013
! 	  if (prev_turbmodel.eq.1)then
! 	  read (1083)INITIALRES(5)
! 	  end if
! 	   if (prev_turbmodel.eq.2)then
! 	  read (1083)INITIALRES(5)
! 	  read (1083)INITIALRES(6)
! 	  end if
! 	  
! 
!    
!    else  !It mean lamy.gt.0, coz it was Transient before
!    READ(1083)   !We read and drop: it, restime
!    rescounter=0
!    end if
! 	
! 
!                         
!        
! 		
! 	    DO I=1,IMAXE
! 			jkn=i
! 
! 	
! 			IF (TURBULENCE.EQ.1)THEN		!IF TURBULENCE
! 			
! 			IF (IRES_TURB.LT.1)THEN			!NON TURBULENT
! 			 
! 			  READ (1083)RG(1:4+iii)
! 		  
! 			  IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(1,1:4)=RG(1:4)
! 				  IF (TURBULENCEMODEL.EQ.1)THEN
! 				  !RG(6)=VISC*TURBINIT
! 				  U_CT(KI)%VAL(1,1)=VISC*TURBINIT/RG(1)
!                                   END IF
! 				  IF (TURBULENCEMODEL.EQ.2)THEN
! 
! 				  U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				   U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 
! 				
! 				  END IF
! 				  IF (LAMPS.LT.1)THEN
! 				  IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=ZERO
! 				  END IF
! 				  
! 				  ELSE
! 				  
! 				   IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=RG(5+1)  !LAMX is lower than 1
! 				  END IF
! 				  END IF
! 			   
! 			  END IF
! 
! 	!-----------------------------------------------------		  
! 			  !End of the corrections on 21/6/2013
! !-------------------------------------------------------			  
! 		!Modified on the 22/6/2013	
! 			
! 			ELSE   !For LAMX.gt.1
! ! 			  DO K=1,4+prev_turbequation+III
! ! 				READ (1083)RG(K)
! !                                 
! ! 			  END DO
! 			  READ (1083)RG(1:4+prev_turbequation+III)
!                            
! 			  IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(1,1:4)=RG(1:4)
! 				if (Turbulence .eq. 1) then
! 				   if (turbulencemodel .eq.prev_turbmodel)then
! 				   U_CT(KI)%VAL(1,1:Turbulenceequations)=RG(5:4+Turbulenceequations)
! 
!      			   end if
! 				   
! 				   if ((turbulencemodel.eq.1).and.(prev_turbmodel.eq.2))then
! 				   U_CT(KI)%VAL(1,1) =alpha_starinf* RG(5)/RG(6)!Nu turbulent
! 				   end if
! 				   
! 				   if ((turbulencemodel.eq.2).and.(prev_turbmodel.eq.1))then
! 				 U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				  U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 				 
! 				 !U_CT(KI)%VAL(1,1)=(1.5*(I_turb_inlet*sqrt(RG(2)**2+RG(3)**2+RG(4)**2)/RG(1))**2)*RG(1)
! 				 !U_CT(KI)%VAL(1,2)=(alpha_starinf*max(U_CT(KI)%VAL(1,1), 1.5*RG(1)*(I_turb_inlet*ufreestream)**2)/&
! 				  !                      max(RG(6),INIT_mu_RATIO*visc)*RG(1))  !Already includes Rho for rho*om
! 				   end if
! 
!                       
! 				   IF (LAMPS.LT.1)THEN
! 				  IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=0.0
! 				  END IF
! 				  
! 				  ELSE
! 				 
! 				   IF (PASSIVESCALAR.GT.0)THEN
! 			
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=RG(4+prev_turbequation+1)!*RG(1)
! 				  END IF
! 				  END IF
! 			    end if
! 				 
! 			  END IF
! 
! 
! 
! 			 END IF
! 			 
! 			ELSE ! turbulence is zero
! 			READ (1083)RG(1:4+prev_turbequation+III)
! ! 			DO K=1,4+prev_turbequation+III
! ! 				READ (1083)RG(K)
! ! 			END DO
! 		
! 			IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(1,1:4)=RG(1:4)
! 			if (passivescalar .gt.0) then
! 
! 
! 			  IF (LAMPS.GE.1)THEN
! 			  U_CT(KI)%VAL(1,turbulenceequations+1)=RG(4+prev_turbequation+1)
! 			  ELSE
! 			  U_CT(KI)%VAL(1,turbulenceequations+1)=ZERO
! 			  END IF
! 
! 
! 
! 			end if
! 			  
! 			END IF
! 
! 			
! 
! 			END IF	!End if for turbulence
! 
! 			END DO
! 			
! 	
! 	
! 	
! 	
! 	CLOSE(1083)
! 
! 	
! !----------------------------------
! 	else   !This else means the now is RUNGE KUTTA for transient       !LEV0
! ! 	
! ! 	
! 	!!! NO LOCAL TIME STEP
! ! 	
! 	RESTFILE='RESTART.dat'
! 	OPEN(1083,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	
! 	IF ((IRES_UNSTEADY.LT.1))THEN 
! 	READ (1083) !resiter
!          iterr=0
!         res_time=ZERO
!         ELSE
!         READ (1083)ITERR,RES_TIME
!         END IF
! 	
! 	IF ((IRES_UNSTEADY.LT.1))THEN       !LEV0
! 	
! 	DO K=1,4+prev_turbequation
! 	READ (1083)
! 	  
! 
! 
! 
! 	END DO
! 	  IF (LAMPS.GT.1.0)THEN
! 	  READ (1083)
! 
! 	  END IF
! 
! 
! 	END IF
! 	
! 	
! 	DO I=1,IMAXE
! 
! 			jkn=i
! 			IF ((TURBULENCE.EQ.1))THEN     !LEV1
! 			IF (IRES_TURB.LT.1)THEN           !LEV2
! ! 			DO K=1,4+III
! ! 				READ (1083)RG(K)
! ! 
! ! 			END DO
! 			READ (1083)RG(1:4+iii)
! 
! 			IF (XMPIE(jkn).EQ.N)THEN        !LEV3
! 			      KI=XMPIL(jkn)
! 			      U_C(KI)%VAL(1,1:4)=RG(1:4)
! ! 				
! 				
! 				if (TURBULENCEMODEL.eq.1) then
! 				U_CT(KI)%VAL(1,1)=VISC*TURBINIT/RG(1)
! 				end if
! 
! 				if (TURBULENCEMODEL .eq. 2) then
! 				  U_CT(KI)%VAL(1,1)=1.5D0*(I_turb_inlet*ufreestream)**2
! 				  U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25D0))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 				  
! 				end if
! 				
! 
!                             if (PASSIVESCALAR .gt. 0) then
!                               if (LAMPS .gt. 0) then
!                               U_CT(KI)%VAL(1,turbulenceequations+1)=RG(5)
! 			   else
!                               U_CT(KI)%VAL(1,turbulenceequations+1)=ZERO
!                             end if
! 			   end if
! 
! 		
! 
! 				
! 			END IF                          !LEV3
! 			ELSE  !So LAMX is gt 1          !LEV2
! 			READ (1083)RG(1:4+prev_turbequation+III )
! 			
! 			IF (XMPIE(jkn).EQ.N)THEN        !LEV3
! 			      KI=XMPIL(jkn)
! 			      U_C(KI)%VAL(1,1:4)=RG(1:4)
! 				!BUG IN THE ORIGINAL FILE (No read of RG(6), so it always pass zero)
! 				!U_CT(KI)%VAL(1,1:+Passivescalar)=RG(6:Turbulenceequations+Passivescalar)
! 		     	   if (turbulencemodel .eq.prev_turbmodel)then
! 				   U_CT(KI)%VAL(1,1:Turbulenceequations)=RG(5:4+Turbulenceequations)
!      			   end if
! 				   
! 				   if ((turbulencemodel.eq.1).and.(prev_turbmodel.eq.2))then
! 				   U_CT(KI)%VAL(1,1) =VISC*TURBINIT
! 				   end if
! 				   
! 				   if ((turbulencemodel.eq.2).and.(prev_turbmodel.eq.1))then
! 				   U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				   U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 			/L_TURB_INLET*RG(1)
! 				  ! U_CT(KI)%VAL(1,1)=(1.5*(I_turb_inlet*sqrt(RG(2)**2+RG(3)**2+RG(4)**2)/RG(1))**2)*RG(1)
! 				  ! U_CT(KI)%VAL(1,2)=(alpha_starinf*U_CT(KI)%VAL(1,1)/RG(6))*RG(1)
! 				   end if
! 
! 			    IF (PASSIVESCALAR.GT.0)THEN
! 			    IF (LAMPS.GT.1)THEN
! 			    
! 
! 
! 			    U_CT(KI)%VAL(1,turbulenceequations+1)=RG(4+prev_turbequation+1)
! 
! 
! 
! 
! 			    ELSE
! 			    U_CT(KI)%VAL(1,turbulenceequations+1)=ZERO
! 
! 
! 
! 
! 			    END IF
! 			    END IF
! 				 
! 			END IF                  !LEV3
! 
! 
! 			
! 			END IF                  !LEV2
! 			ELSE   !So turbulence is 0    !LEV1
! 			READ (1083)RG(1:4+prev_turbequation+iii)
! 			
! 			IF (XMPIE(jkn).EQ.N)THEN  !LEV2
! 			      KI=XMPIL(jkn)
! 			      U_C(KI)%VAL(1,1:4)=RG(1:4)
! 				
! 				 
! 			IF (PASSIVESCALAR.GT.0)THEN
! 			    IF (LAMPS.GT.1.0)THEN
! 			    
! 			    U_CT(KI)%VAL(1,1)=RG(4+prev_turbequation+1)
! 
! 
! 			    ELSE
! 			    U_CT(KI)%VAL(1,1)=ZERO
! 
! 
! 
! 
! 			    END IF
! 			    END IF  
! 			END IF       
! 
! 			END IF        !LEV2
! 			
! 
! 			    
! 	END DO
! 	CLOSE(1083)
! 
!     end if    !LEV1
! END IF    !LEV0
! 
! 
! if (Averaging .EQ. 1) then
! 
! if (average_restart .eq. 1) then
! 	    IF (LAMPS.GE.1)THEN
! 	     allocate(Arg(5+prev_turbequation+4))
! 
! 	    END IF
! 	    IF (LAMPS.LT.1)THEN
! 	     allocate(Arg(4+prev_turbequation+3))
! 
! 	    END IF
! !if ( RUNGEKUTTA .ge. 5 ) Then
!       
! 	  RESTFILE='RESTARTav.dat'
! 	  OPEN(1084,FILE=RESTFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
! 	  
! 	    
! 		
! 	    DO I=1,IMAXE
! 			READ(1084,*)JKN
! 			IF (LAMPS.GE.1)THEN
! 			
! 			READ(1084)ARG(1:6+prev_turbequation+6+MIN(PASSIVESCALAR,1))
! 			
! 
! 	    END IF
! 	    IF (LAMPS.LT.1)THEN
! 		    READ(1084)ARG(1:5+prev_turbequation+6)
! 	    
! 	    END IF
! 
! 			
! 			IF (XMPIE(jkn).EQ.N)THEN
! 			      KI=XMPIL(jkn)
! 			       U_C(KI)%VAL(5,1:4)=ARG(1:4)
! 				  
! 			      
! 			  if (Turbulence .eq. 1) then
! 				   if (turbulencemodel .eq.prev_turbmodel)then
! 				   U_CT(KI)%VAL(5,1:Turbulenceequations)=ARG(5:4+Turbulenceequations)
! 				    end if
! 				   
! 				   if ((turbulencemodel.eq.1).and.(prev_turbmodel.eq.2))then
! 				   U_CT(KI)%VAL(5,1) =VISC*TURBINIT/ARG(1)
! 				   end if
! 				   
! 				   if ((turbulencemodel.eq.2).and.(prev_turbmodel.eq.1))then
! 				   U_CT(KI)%VAL(5,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				   U_CT(KI)%VAL(5,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(5,1))&
! 			/L_TURB_INLET*RG(1)
! 				   !U_CT(KI)%VAL(5,2)=(alpha_starinf*U_CT(KI)%VAL(1,1)/ARG(6))
! 				   end if
! 			  			
! 			ELSE ! turbulence is zero
! 			
! 			       U_C(KI)%VAL(5,1:4)=ARG(1:4)
! 				 
! 			END IF
! 			
! 			
! 			IF (PASSIVESCALAR .EQ. 1) then
!     
! 			IF (LAMPS.GE.1)THEN
! 
! 			U_CT(KI)%VAL(5,turbulenceequations+1)=ARG(4+prev_turbequation+1)
! 
! 
! 			ELSE
! 			U_CT(KI)%VAL(5,turbulenceequations+1)=ZERO
! 			END IF
! 			END IF
! 			
! 			U_C(KI)%RMS(1:5+III)=ARG(4+PREV_TURBEQUATION+III+1:7+PREV_TURBEQUATION+III+III)
! 
! 			
! 		
! 
! 
! 
! 	
! 	      END IF
! 	
! 
! 	    END DO
! 	
! 	CLOSE(1084)
! !END if
! deallocate(Arg)
! else !it means that average_restart=0
!  DO I=1,kmaxe
! U_C(i)%VAL(5,:)=ZERO
! U_C(i)%RMS(:)=ZERO
! if ((passivescalar.gt.0).or.(turbulence.eq.1))then
! 
! U_CT(i)%VAL(5,:)=ZERO
! end if   
! end do
! 
! end if   !Level_Averaging
! 
! 
! 
! end if   !Level_RESTART=1/0

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 
! deallocate(POX)
! deallocate(POy)

deallocate(rg)
DEALLOCATE(VECCOS)


END SUBROUTINE INITIALISE2d

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE INITIALISATION
