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

SUBROUTINE INITIALISE(N)
 !> @brief
!> This subroutine calls the initialisation of the computational domain
IMPLICIT NONE
integer,INTENT(IN)::N
REAL,allocatable,DIMENSION(:)::RG,ARG
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
INTEGER::prev_turbequation,INITIAL,III,i,k,jx,QQP,INC,kmaxe,jkn,ki,iterr,JX2,ind1,KX,icompwrt,ICONSIDERED
INTEGER::ELTYPE,ELEM_DEC,COUNT_1
REAL::VOLTEMP
real,dimension(1:nof_Variables+turbulenceequations+passivescalar)::veccos
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,dimension(1:dimensiona)::CORDS
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D



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
 allocate(rg(nof_Variables+prev_turbequation+1))
  ELSE
 allocate(rg(nof_Variables+prev_turbequation))
END IF





 

VECCOS(1:NOF_VARIABLES)=0.0
KMAXE=XMPIELRANK(N)





IF (RESTART.EQ.0)THEN

 

	IF (ISCHEME.LT.2)THEN
		IF (ITESTCASE.EQ.0)THEN
!$OMP DO
		DO INITIAL=1,KMAXE
			U_C(INITIAL)%VAL(1,1)=1.0D0
			U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
		END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.EQ.1)THEN
!$OMP DO
		DO INITIAL=1,KMAXE
			pox(1)=IELEM(N,INITIAL)%XXC
			poy(1)=IELEM(N,INITIAL)%yyC
			if (dimensiona.eq.3)then
			poz(1)=IELEM(N,INITIAL)%zzC
			U_C(INITIAL)%VAL(1,1)=LINEAR_INIT3D(N,pox,poy,poz)
			ELSE
			U_C(INITIAL)%VAL(1,1)=LINEAR_INIT2D(N,pox,poy,poz)
			end if


			U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
		END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.EQ.2)THEN
!$OMP DO
		DO INITIAL=1,KMAXE
			pox(1)=IELEM(N,INITIAL)%XXC
			poy(1)=IELEM(N,INITIAL)%yyC
			if (dimensiona.eq.3)then
			poz(1)=IELEM(N,INITIAL)%zzC
			U_C(INITIAL)%VAL(1,1)=LINEAR_INIT3D(N,pox,poy,poz)
			ELSE
			U_C(INITIAL)%VAL(1,1)=LINEAR_INIT2D(N,pox,poy,poz)
			end if

			U_E(INITIAL)%VAL(1,1)=U_C(INITIAL)%VAL(1,1)
		END DO
!$OMP END DO
		END IF
		IF (ITESTCASE.GE.3)THEN
!$OMP DO
		DO INITIAL=1,KMAXE
			VECCOS(:)=ZERO
			pox(1)=IELEM(N,INITIAL)%XXC
			poy(1)=IELEM(N,INITIAL)%yyC
			if (dimensiona.eq.3)then
			poz(1)=IELEM(N,INITIAL)%zzC
			CALL INITIALISE_EULER3D(N,veccos,pox,poy,poz)
			ELSE
			CALL INITIALISE_EULER2D(N,veccos,pox,poy,poz)
			end if

			if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
			U_C(INITIAL)%VAL(1,1:nof_Variables)=VECCOS(1:nof_Variables)
			U_CT(INITIAL)%VAL(1,1:0+turbulenceequations+passivescalar)=VECCOS(nof_Variables+1:nof_Variables+turbulenceequations+passivescalar)
			else
			U_C(INITIAL)%VAL(1,:)=VECCOS(:)
			if (itestcase.ge.3)U_E(INITIAL)%VAL(1,:)=U_C(INITIAL)%VAL(1,:)
			end if
		END DO
!$OMP END DO
		END IF
	ELSE


!$OMP DO
	DO I=1,KMAXE
		
		iconsidered=i
		 
		 VEXT=ZERO
        NODES_LIST=ZERO
        ELTYPE=IELEM(N,I)%ISHAPE
        ELEM_DEC=IELEM(N,I)%VDEC
        ELEM_LISTD=ZERO
        VOLTEMP=ZERO
        
      jx=IELEM(N,I)%NONODES
!       
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,:)=inoder(JX2)%CORD(:)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO



	  IF (DIMENSIONA.EQ.3)THEN
	  CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
	  ELSE
	  CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
	  END IF
    
      SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
       COUNT_1=0
      if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREHEXA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      VOLTEMP=HEXAVOLUME(N,VEXT,QPOINTS,WEQUA3D)/IELEM(N,I)%totvolume
      QQP=QP_HEXA
      
	  DO INC=1,QQP
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 



            CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
	  
!
	  END DO
      
      
      else
      COUNT_1=0
       do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	  
	  CALL QUADRATUREtetra(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
                        if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N,VEXT)
                        else
                        VOLTEMP=TETRAVOLUME(N,VEXT)/IELEM(N,I)%totvolume
                        end if
	  
	 
	  
                    QQP=QP_TETRA
                    DO INC=1,QQP
                    COUNT_1=COUNT_1+1
                    POX(1)=QPOINTS(1,INC)
                    POY(1)=QPOINTS(2,INC)
                    POZ(1)=QPOINTS(3,INC)

                    CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
                    END DO
       END DO
       end if
	
	! 	

      CASE(2)

       
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	
       CALL QUADRATUREtetra(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
       
       
       
        
                        if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N,VEXT)
                        COUNT_1=0
                        else
                        VOLTEMP=TETRAVOLUME(N,VEXT)/IELEM(N,I)%totvolume
                        end if
! 	  
	   QQP=QP_TETRA
	  DO INC=1,QQP
	  COUNT_1=COUNT_1+1
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
	  END DO
      CASE(3)
      
        
      COUNT_1=0
      do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
	  
	  CALL QUADRATUREtetra(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
	  if (dg.eq.1)then
                        VOLTEMP=TETRAVOLUME(N,VEXT)
                        else
                        VOLTEMP=TETRAVOLUME(N,VEXT)/IELEM(N,I)%totvolume
                        end if
	     QQP=QP_TETRA
	  DO INC=1,QQP
	   COUNT_1=COUNT_1+1
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  
			CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
	  END DO
       END DO
	
	
	
	
       CASE(4)
       
       if (IELEM(N,I)%MODE.eq.0)then
       CALL QUADRATUREPRISM(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      VOLTEMP=PRISMVOLUME(N,VEXT,QPOINTS,WEQUA3D)/IELEM(N,I)%totvolume
       QQP=QP_PRISM
	  DO INC=1,QQP
	  POX(1)=QPOINTS(1,INC)
	  POY(1)=QPOINTS(2,INC) 
	  POZ(1)=QPOINTS(3,INC) 
	  COUNT_1=INC
			CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
	  END DO
      else
      COUNT_1=0
        do K=1,ELEM_DEC
        VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)

        CALL QUADRATUREtetra(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
        if (dg.eq.1)then
                            VOLTEMP=TETRAVOLUME(N,VEXT)
                            else
                            VOLTEMP=TETRAVOLUME(N,VEXT)/IELEM(N,I)%totvolume
                            end if
            QQP=QP_TETRA
        DO INC=1,QQP
        COUNT_1=COUNT_1+1
        POX(1)=QPOINTS(1,INC)
        POY(1)=QPOINTS(2,INC)
        POZ(1)=QPOINTS(3,INC)

                CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
        END DO
        END DO
       end if

       
        CASE(5)

         IF (IELEM(N,I)%MODE.EQ.0)THEN
                CALL QUADRATUREQUAD(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)

                            VOLTEMP=1.0d0
                            QQP=QP_quad

                            DO INC=1,QQP
                            COUNT_1=QQP
                                POX(1)=QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                                POY(1)=QPOINTS(2,INC)

                                CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)

                            END DO
            ELSE
                    COUNT_1=0

                        VOLTEMP=0.0d0
                        WEQUA3D=0.0d0
                        QPOINTS=0.0d0
                            DO K=1,ELEM_DEC
                                VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)

                            CALL QUADRATUREtriangle(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
                                if (dg.eq.1)then
                                VOLTEMP=TRIANGLEVOLUME(N,VEXT)
                                else
                                VOLTEMP=TRIANGLEVOLUME(N,VEXT)/IELEM(N,I)%totvolume
                                end if


                                QQP=QP_Triangle
                                DO INC=1,QQP
                                        COUNT_1=COUNT_1+1
                                        POX(1) = QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                                    POY(1) = QPOINTS(2,INC)
                                    CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
                                END DO
                            end do
            END IF








        CASE(6)

                COUNT_1=0


                      CALL QUADRATUREtriangle(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
                        if (dg.eq.1)then
                        VOLTEMP=TRIANGLEVOLUME(N,VEXT)
                        else
                        VOLTEMP=TRIANGLEVOLUME(N,VEXT)/IELEM(N,I)%totvolume
                        end if
                        QQP=QP_Triangle
                         DO INC=1,QQP
                                COUNT_1=COUNT_1+1
                                POX(1) = QPOINTS(1,INC) !POX,POY required for LINEAR_INIT2D
                               POY(1) = QPOINTS(2,INC)
                               CALL ASSIGN_INITIAL(I,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
                        END DO








      END SELECT
   
    		
		
		
		
		
	END DO
!$OMP END DO
	END IF

	RES_TIME=ZERO
END IF
INITIALRES=0.0d0

!$OMP BARRIER
 

deallocate(rg)



END SUBROUTINE INITIALISE



SUBROUTINE ASSIGN_INITIAL(ICONSIDERED,INC,VOLTEMP,WEQUA3D,POX,POY,POZ,COUNT_1)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,inc,COUNT_1
real,dimension(1:dimensiona),INTENT(IN)::pox,poy,poz
REAL,INTENT(IN)::VOLTEMP
real,dimension(1:nof_Variables+turbulenceequations+passivescalar)::veccos
REAL,DIMENSION(1:idegfree+1)::basis_vector
REAL,DIMENSION(1,1:idegfree+1)::tempsol
REAL,DIMENSION(1:NUMBEROFPOINTS),INTENT(IN)::WEQUA3D
INTEGER::i,kx,icompwrt

I=ICONSIDERED




        IF (ITESTCASE.LE.2)THEN

			IF (DG.EQ.1)THEN
			basis_vector(1)=1.0d0
            icompwrt=-2

            IF (DIMENSIONA.EQ.2)THEN
             BASIS_VECTOR(2:idegfree+1) = BASIS_REC2D(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),IORDER,I,IDEGFREE,icompwrt, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)

            ELSE
            BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE,icompwrt, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)
            END IF


            IF (DIMENSIONA.EQ.2)THEN

            tempsol(1,:)=LINEAR_INIT2D(N,POX,POY,POZ)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)
            Else
            tempsol(1,:)=LINEAR_INIT3D(N,POX,POY,POZ)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)

            END IF



            U_C(I)%VALDG(1,1,:)=U_C(I)%VALDG(1,1,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
			ELSE


             IF (DIMENSIONA.EQ.2)THEN

			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT2D(N,POX,POY,POZ)*WEQUA3D(INC)*(VOLTEMP)
			Else
			U_C(I)%VAL(1,1)=U_C(I)%VAL(1,1)+LINEAR_INIT3D(N,POX,POY,POZ)*WEQUA3D(INC)*(VOLTEMP)

			END IF



			END IF
			U_E(I)%VAL(1,1)=U_C(I)%VAL(1,1)
			IF (DG.EQ.1)U_E(I)%VAL(1,1)=U_C(I)%VALDG(1,1,1)
        ELSE



           IF (DIMENSIONA.EQ.2)THEN
          CALL INITIALISE_EULER2D(N,veccos,pox,poy,poz)
          Else
          CALL INITIALISE_EULER3D(N,veccos,pox,poy,poz)
          END IF




                if ((turbulence .eq. 1).or.(passivescalar.gt.0)) then
                U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP)
                U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)=U_CT(I)%VAL(1,1:0+turbulenceequations+passivescalar)+&
                VECCOS(NOF_VARIABLES+1:NOF_VARIABLES+turbulenceequations+passivescalar)*WEQUA3D(INC)*(VOLTEMP)
                else
                IF (DG.EQ.1)THEN
                    basis_vector(1)=1.0d0
                    icompwrt=-2


                    IF (DIMENSIONA.EQ.2)THEN
                    BASIS_VECTOR(2:idegfree+1) = BASIS_REC2D(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),IORDER,I,IDEGFREE,icompwrt, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)

                    Else
                     BASIS_VECTOR(2:idegfree+1) = BASIS_REC(N,QP_ARRAY(I)%X(COUNT_1),QP_ARRAY(I)%Y(COUNT_1),QP_ARRAY(I)%Z(COUNT_1),IORDER,I,IDEGFREE,icompwrt, IELEM, ILOCAL_RECON3, INTEG_BASIS,integ_basis_dg)

                    END IF


                                        DO KX=1,nof_Variables
                                    tempsol(1,:)=VECCOS(KX)*WEQUA3D(INC)*(voltemp)*basis_vector(1:idegfree+1)

                                    U_C(I)%VALDG(1,KX,:)=U_C(I)%VALDG(1,KX,:)+MATMUL(m_1(i)%val(:,:),tempsol(1,:))
                                        END DO
                ELSE



                    U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+(VECCOS(1:nof_Variables)*WEQUA3D(INC)*(VOLTEMP))



                if (itestcase.ge.3)U_E(I)%VAL(1,:)=U_C(I)%VAL(1,:)


                END IF
                end if
			END IF




END SUBROUTINE ASSIGN_INITIAL






! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE INITIALISATION
