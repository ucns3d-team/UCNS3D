MODULE MOODR
USE MPIINFO
USE TRANSLATE
use DECLARATION
USE MEMORY
USE COMMUNICATIONS
USE IO
USE PARTITION
USE LIBRARY
USE TRANSFORM
USE FLUXES
USE INITIALISATION
USE BOUNDARY
USE RECON
USE LOCAL
USE PROFILE
USE FLOW_OPERATIONS
USE GRADIENTS
USE BASIS
USE PRESTORE
USE RIEMANN
USE SOURCE
USE implicit_time
USE implicit_FLUXES
USE OMP_LIB




IMPLICIT NONE

 CONTAINS


SUBROUTINE PAD_NAD(N)

    IMPLICIT NONE

    INTEGER,INTENT(IN)::N
    INTEGER::I,L,NGP,iqp,iex,kmaxe,K,ii,iconsidered
    INTEGER::REDUCE1
    REAl,DIMENSION(NOF_vARIABLES)::NAD_DELTA1
    INTEGER::PAD_TRUE,NAD_TRUE
    REAL::RELAX_MOOD1,RELAX_MOOD2
    real,dimension(1:nof_Variables)::leftv
    real::MP_PINFL,gammal
    real,dimension(1:nof_Variables)::rightv
    real::MP_PINFr,gammar
    REAL,allocatable,DIMENSION(:,:)::UTEMP
    REAL,DIMENSION(1:NOF_VARIABLES)::UTMIN,UTMAX


    allocate(utemp(IMAXDEGFREE+1,1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR))


    KMAXE=XMPIELRANK(N)


    IF (CASCADE.EQ.1)THEN
        RELAX_MOOD1=MOOD_VAR1
        RELAX_MOOD2=MOOD_VAR2
    END IF
    IF (CASCADE.EQ.2)THEN
        RELAX_MOOD1=MOOD_VAR3
        RELAX_MOOD2=MOOD_VAR4
    END IF



    IF (ITESTCASE.GE.1) THEN

        !$OMP DO
        DO II=1,NOF_INTERIOR
            I=EL_INT(II)
            ICONSIDERED=I
            IELEM(N,I)%MOOD=0
            REDUCE1=0
            PAD_TRUE=0
            NAD_TRUE=0
            
            !1 copy candidate solution at temp variable   
            LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)

            !2 TRANSFORM CONSERVATIVE  TO PRIMITIVE AND CHECK IF PRESSURE AND DENSITY ARE PHYSICALLY ADMISSIBLE IF NOT PAD_TRUE=1

            IF (ITESTCASE.GE.3) THEN ! for advection equation there is no need to transfor to primitive
                CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
            END IF
            
            IF ((LEFTV(1).LE.ZERO).OR.(LEFTV(1).NE.LEFTV(1)))THEN						
                PAD_TRUE=1
            END IF

            IF ((LEFTV(NOF_VARIABLES).LE.ZERO).OR.(LEFTV(NOF_VARIABLES).NE.LEFTV(NOF_VARIABLES))) THEN ! LEFTV(NOF_VARIABLES) = pressure or NOF_VARIABLES = 1						
                PAD_TRUE=1
            END IF
        
        
            !3 THE ONES WITH A PHYSICALLY ADMISSIBLE SOLUTION NEED TO GET CHECKED
        
            IF (PAD_TRUE.EQ.0)THEN
                
                UTEMP=ZERO
                            
                !4 NOW ESTABLISH A TEMPORARY ARRAY WITH THE CURRENT SOLUTION FROM THE DIRECT SIDE NEIGHBOURS OF CONSIDERED CELL

                K=0
                UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(3,1:NOF_VARIABLES)
                K=1
                DO L=1,IELEM(N,I)%IFCA
                    K=K+1
                    UTEMP(K,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(3,1:NOF_VARIABLES)
                END DO
                                    
                !5 NOW ESTABLISH THE MIN AND MAX BOUNDS

                DO IEX=1,NOF_VARIABLES
                    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
                    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
                END DO
                    
                ! MOOD_VAR1=0.0001;MOOD_VAR2=0.001; MOOD_MODE=RELAXED
                
                DO IEX=1,NOF_VARIABLES
                    NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                END DO

                NAD_TRUE=0

                !6 SPECIFY RELAXED OR ORIGINAL MOOD PATTERN
                IF (MOOD_MODE.GT.0)THEN
                    DO IEX=1,NOF_VARIABLES
                
                        IF ((IEX.EQ.1).OR.(IEX.EQ.NOF_vARIABLES))then
                            IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                                NAD_TRUE=1
                            END IF
                        endif
                
                    END DO
                ELSE
                    DO IEX=1,NOF_VARIABLES
                        IF ((IEX.EQ.1).OR.(IEX.EQ.NOF_vARIABLES))then
                            IF ((U_C(I)%VAL(4,IEX).lt.(UTMIN(IEX))).or.(U_C(I)%VAL(4,IEX).gt.(UTMAX(IEX))))THEN
                                NAD_TRUE=1
                            END IF
                        end if
                    END DO
                END IF
                    
            END IF

            !7 NOW SET THE MOOD FLAG FOR EACH ELEMENT
            IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE
                
        END DO
        !$OMP END DO		


        !$OMP DO
        DO II=1,NOF_BOUNDED
            I=EL_BND(II)
            ICONSIDERED=I
            REDUCE1=0
            PAD_TRUE=0
            NAD_TRUE=0
            IELEM(N,I)%MOOD=0
            LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)

            IF (ITESTCASE.GE.3) THEN ! for advection equation there is no need to transfor to primitive
                IF (DIMENSIONA.EQ.3) THEN
                    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
                ELSE
                    CALL cons2prim(N,leftv,MP_PINFl,gammal)                   
                END IF
            END IF

            IF ((LEFTV(1).LE.ZERO).OR.(LEFTV(1).NE.LEFTV(1))) THEN						
                PAD_TRUE=1
            END IF

            IF ((LEFTV(NOF_VARIABLES).LE.ZERO).OR.(LEFTV(NOF_VARIABLES).NE.LEFTV(NOF_VARIABLES))) THEN	! LEFTV(NOF_VARIABLES) = pressure or NOF_VARIABLES = 1						
                PAD_TRUE=1
            END IF
                
            
            IF (PAD_TRUE.EQ.0)THEN
                
                UTEMP=ZERO
                K=0
                UTEMP(1,1:NOF_VARIABLES)=U_C(I)%VAL(3,1:NOF_VARIABLES)
                K=1  

                DO L=1,IELEM(N,I)%IFCA	!faces2
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N) THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0) THEN	!CHECK FOR BOUNDARIES
                            IF (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5) THEN	!PERIODIC IN MY CPU
                                K=K+1
                                UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
                            ELSE
                                !NOT PERIODIC ONES IN MY CPU			  				  
                            END IF
                        ELSE
                            K=K+1
                            UTEMP(K,1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
                        END IF

                    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
            
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                K=K+1
                                UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                                (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                            END IF
                        ELSE
                            K=K+1
                            UTEMP(K,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
                            (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables) 
                        END IF
                    END IF
                END DO
                        
                DO IEX=1,NOF_VARIABLES
                    UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
                    UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
                END DO
                    
                DO IEX=1,NOF_VARIABLES
                    NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                END DO
                
                NAD_TRUE=0
                                
                                
                IF (MOOD_MODE.GT.0) THEN
                    DO IEX=1,NOF_VARIABLES
                        IF ((IEX.EQ.1).OR.(IEX.EQ.NOF_vARIABLES)) THEN
                            IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                                NAD_TRUE=1
                            END IF
                        END IF           
                    END DO
                ELSE
                    DO IEX=1,NOF_VARIABLES
                        IF ((IEX.EQ.1).OR.(IEX.EQ.NOF_vARIABLES))then
                            IF ((U_C(I)%VAL(4,IEX).lt.(UTMIN(IEX))).or.(U_C(I)%VAL(4,IEX).gt.(UTMAX(IEX))))THEN
                                NAD_TRUE=1
                            END IF
                        END IF
                    END DO                
                END IF
                                
            END IF
                
            IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE

        END DO 
        !$OMP END DO

    END IF		
		
	DEALLOCATE(UTEMP)

END SUBROUTINE PAD_NAD






SUBROUTINE MOOD_OPERATOR_2(N)
    IMPLICIT NONE
    INTEGER::I,KMAXE
    INTEGER,INTENT(IN)::N

    KMAXE=XMPIELRANK(N)

    !$OMP DO
    DO I=1,KMAXE
        IELEM(N,I)%RECALC=0
        IELEM(N,I)%MOOD=0
    END DO
    !$OMP END DO

    CASCADE=1
    CALL PAD_NAD(N)

    CALL EXHBOUNDHIGHER_MOOD(N)
 
    CALL FIX_LIST(N)
 
    CALL MUSCL(N)
 
    CALL EXHBOUNDHIGHER(N)
 
    IF (dimensiona.eq.2) THEN
        IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN ! linear advection equation
            CALL CALCULATE_FLUXESHI2D_MOOD(N)
        ELSE ! Euler or Navier-Stokes
            CALL CALCULATE_FLUXESHI_CONVECTIVE2D_MOOD(N)
        END IF
    ELSE
        CALL CALCULATE_FLUXESHI_CONVECTIVE_MOOD(N)
    END IF
 
    !$OMP DO
    DO I=1,KMAXE
        IF (IELEM(N,I)%MOOD.EQ.0) THEN
            IELEM(N,I)%MOOD_O=IORDER+1  !TARGET POLYNOMIAL/SCHEME ADMISSIBLE
        ELSE
            IELEM(N,I)%MOOD_O=2         !SECOND ORDER MUSCL ADMISSIBLE
        END IF
    END DO
    !$OMP END DO

END SUBROUTINE MOOD_OPERATOR_2





SUBROUTINE MOOD_OPERATOR_1(N)
    IMPLICIT NONE
    INTEGER::I,KMAXE
    INTEGER,INTENT(IN)::N
    KMAXE=XMPIELRANK(N)

    !$OMP DO
    DO I=1,KMAXE
        IELEM(N,I)%RECALC=0
        IELEM(N,I)%MOOD=0
    END DO
    !$OMP END DO

    CASCADE=2
    CALL PAD_NAD(N)

    CALL EXHBOUNDHIGHER_MOOD(N)
    
    CALL FIX_LIST(N)
    
    IF (dimensiona.eq.2) THEN
        IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN ! linear advection equation
            CALL CALCULATE_FLUXESHI2D_MOOD(N)
        ELSE ! Euler or Navier-Stokes
            CALL CALCULATE_FLUXESHI_CONVECTIVE2D_MOOD(N)
        END IF
    ELSE
        CALL CALCULATE_FLUXESHI_CONVECTIVE_MOOD(N)
    END IF

    !$OMP DO
    DO I=1,KMAXE
        IF (IELEM(N,I)%MOOD.EQ.1)THEN
            IELEM(N,I)%MOOD_O=1     ! ONLY FIRST ORDER SOLUTION ADMISSIBLE
        END IF
    END DO
    !$OMP END DO

END SUBROUTINE MOOD_OPERATOR_1


SUBROUTINE FIX_LIST(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL::GODFLUX2,sum_detect
	INTEGER::I,L,NGP,KMAXE,IQP
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_TEMP
	REAL,DIMENSION(1)::CRIGHT
	KMAXE=XMPIELRANK(N)
	
	
	!$OMP DO
	DO I=1,KMAXE
        IELEM(N,I)%RECALC=0
		IF (IELEM(N,I)%MOOD.EQ.1) THEN
		
            IELEM(N,I)%RECALC=1
        ELSE
            IF (IELEM(N,I)%INTERIOR.EQ.0) THEN
                CRIGHT(1)=ZERO
                
                DO L=1,IELEM(N,I)%IFCA
                    
                    cRIGHT(1)=IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD 
                    if (cright(1).gt.0.5)then
                        IELEM(N,I)%RECALC=1
                    end if	    
                END DO
            END IF
            
            
            IF (IELEM(N,I)%INTERIOR.EQ.1) THEN

                CRIGHT(1)=ZERO
                
                DO L=1,IELEM(N,I)%IFCA
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5) then	!PERIODIC IN MY CPU
                                CRIGHT(1)=ielem(n,(IELEM(N,I)%INEIGH(L)))%mood
                            ELSE
                                                                
                            END IF
                        ELSE
                            CRIGHT(1)=ielem(n,(IELEM(N,I)%INEIGH(L)))%mood
            
                        END IF

                    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5) then	!PERIODIC IN OTHER CPU
                                CRIGHT(1)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_m(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1)  
                            END IF
                        ELSE 								
                                CRIGHT(1)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_m(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1)		
                        END IF
                    END IF
                        
                        if (cright(1).gt.0.5) then
                            IELEM(N,I)%RECALC=1
                        end if	    
                
                END DO
            END IF
 		END IF		
	END DO
	!$OMP END DO 
	END SUBROUTINE FIX_LIST

END MODULE MOODR
