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


REAL FUNCTION ENTROPY(RHO, P, local_gamma)
    IMPLICIT NONE
    REAL::RHO,P, local_gamma
    ! write(*,*) "local_gamma: ", local_gamma
    ENTROPY = RHO ** local_gamma
    ENTROPY = (local_gamma - 1.0)*ENTROPY
    ENTROPY = ENTROPY / P
    ENTROPY = LOG(ENTROPY)
    ENTROPY = RHO * ENTROPY

END FUNCTION





SUBROUTINE PAD_NAD(N)

    IMPLICIT NONE

    INTEGER,INTENT(IN)::N
    INTEGER::I,L,D,NGP,iqp,iex,kmaxe,K,num_valid_neighbrs,ii,iconsidered
    INTEGER::REDUCE1
    REAl,DIMENSION(NOF_vARIABLES)::NAD_DELTA1,DELTA
    INTEGER::PAD_TRUE,NAD_TRUE
    REAL::RELAX_MOOD1,RELAX_MOOD2
    real,dimension(1:nof_Variables)::leftv
    real::MP_PINFL,gammal
    real,dimension(1:nof_Variables)::rightv
    real::MP_PINFr,gammar
    real::helper,parameter
    REAL,allocatable,DIMENSION(:,:)::UTEMP
    REAL,DIMENSION(1:NOF_VARIABLES)::UTMIN,UTMAX
    integer:: rho_index, p_index, vf_index
    REAL::CELL_AREA,CELL_SIZE
    REAL::rho_old,rho_new,p_old,p_new,entropy_old,entropy_new,helper_value

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

    select case (ITESTCASE)
      case (1, 2) 
        rho_index = 1
        p_index = 1
        vf_index = 0

      case (3)
        rho_index = 1
        p_index = dimensiona + 2
        vf_index = 0

      case (4) 
        rho_index = 1
        p_index = dimensiona + 2
        vf_index = p_index + 1

      case DEFAULT
        rho_index = 0
        p_index = 0
        vf_index = 0
    end select

    SELECT CASE (MOOD_MODE)

      CASE (0, 1)

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
                
                IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN
                    IF ((LEFTV(rho_index).LE.(-0.001)).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                ELSE
                    IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                    IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
                        PAD_TRUE=1
                    END IF
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
                    IF (MOOD_MODE.EQ.1)THEN
                        DO IEX=1,NOF_VARIABLES
                            IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
                                IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                                    NAD_TRUE=1
                                END IF
                            endif
                        END DO
                    ELSE
                        DO IEX=1,NOF_VARIABLES
                            IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
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

                IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN
                    IF ((LEFTV(rho_index).LE.(-0.001)).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                ELSE
                    IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                    IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
                        PAD_TRUE=1
                    END IF
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
                                    
                    IF (MOOD_MODE.EQ.1) THEN
                        DO IEX=1,NOF_VARIABLES
                            IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index)) THEN
                                IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                                    NAD_TRUE=1
                                END IF
                            END IF           
                        END DO
                    ELSE
                        DO IEX=1,NOF_VARIABLES
                            IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
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
    
      CASE (2)

        IF (ITESTCASE.GE.3) THEN

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
                CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
                
                IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                    PAD_TRUE=1
                END IF

                IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
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

                    ! replace U_x with |U|
                    DO L=1,K ! all considered cells
                        UTEMP(L,2) = UTEMP(L,2)**2
                        DO D=2,DIMENSIONA
                            UTEMP(L,2) = UTEMP(L,2) + (UTEMP(L,1+D)**2)
                        END DO
                        UTEMP(L,2) = SQRT(UTEMP(L,2))
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

                    IEX = rho_index ! density
                    IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                        NAD_TRUE=1
                    END IF
                    
                    IEX = 2 ! speed
                    IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                        NAD_TRUE=1
                    END IF
                    IEx = p_index ! internal energy
                    IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                        NAD_TRUE=1
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
                
                IF (DIMENSIONA.EQ.3) THEN
                    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
                ELSE
                    CALL cons2prim(N,leftv,MP_PINFl,gammal)                   
                END IF
                
                IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index))) THEN						
                    PAD_TRUE=1
                END IF

                IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN
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

                    ! replace U_x with |U|
                    DO L=1,K ! all considered cells
                        UTEMP(L,2) = UTEMP(L,2)**2
                        DO D=2,DIMENSIONA
                            UTEMP(L,2) = UTEMP(L,2) + (UTEMP(L,1+D)**2)
                        END DO
                        UTEMP(L,2) = SQRT(UTEMP(L,2))
                    END DO
                            
                    DO IEX=1,NOF_VARIABLES
                        UTMIN(IEX)=MINVAL(UTEMP(1:K,IEX))
                        UTMAX(IEX)=MAXVAL(UTEMP(1:K,IEX))
                    END DO
                        
                    DO IEX=1,NOF_VARIABLES
                        NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                    END DO
                    
                    IEX = rho_index ! density
                    IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                        NAD_TRUE=1
                    END IF
                    
                    IEX = 2 ! speed
                    IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                        NAD_TRUE=1
                    END IF
                    IEx = p_index ! internal energy
                    IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
                        NAD_TRUE=1
                    END IF
                                    
                END IF
                    
                IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE

            END DO 
            !$OMP END DO
        ELSE
            if (n.eq.0) then
                write(*,*) 'ITESTCASE incompatible with MOOD_MODE 2 in PAD_NAD function'
                call exit(777)
            endif

        END IF	

      CASE (3)

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
            
            IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                PAD_TRUE=1
            END IF

            IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
                PAD_TRUE=1
            END IF
        
            !3 THE ONES WITH A PHYSICALLY ADMISSIBLE SOLUTION NEED TO GET CHECKED
            IF (PAD_TRUE.EQ.0)THEN
                
                !4 NOW ESTABLISH A TEMPORARY ARRAY WITH THE CURRENT SOLUTION FROM THE DIRECT SIDE NEIGHBOURS OF CONSIDERED CELL
                UTEMP=ZERO
                num_valid_neighbrs = 0
                DO L=1,IELEM(N,I)%IFCA
                    helper = U_C(IELEM(N,I)%INEIGH(L))%VAL(4,rho_index)
                    if ((helper.eq.helper).and.(helper.gt.0)) then
                        helper = U_C(IELEM(N,I)%INEIGH(L))%VAL(4,p_index)
                        if ((helper.eq.helper).and.(helper.gt.0)) then
                            num_valid_neighbrs = num_valid_neighbrs + 1
                            UTEMP(num_valid_neighbrs,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(4,1:NOF_VARIABLES)
                        end if
                    end if
                END DO

                if (num_valid_neighbrs.ge.1) then
                                    
                    !5 NOW ESTABLISH THE MIN AND MAX BOUNDS
                    DO IEX=1,NOF_VARIABLES
                        UTMIN(IEX)=MINVAL(UTEMP(1:num_valid_neighbrs,IEX))
                        UTMAX(IEX)=MAXVAL(UTEMP(1:num_valid_neighbrs,IEX))
                    END DO
                        
                    ! MOOD_VAR1=0.0001;MOOD_VAR2=0.001; MOOD_MODE=RELAXED
                    
                    DO IEX=1,NOF_VARIABLES
                        NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                    END DO

                    NAD_TRUE=0

                    DO IEX=1,NOF_VARIABLES
                        IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
                            IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                                NAD_TRUE=1
                            END IF
                        endif
                    END DO
                else
                    NAD_TRUE = 1
                end if
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

            IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index))) THEN						
                PAD_TRUE=1
            END IF

            IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
                PAD_TRUE=1
            END IF
                
            IF (PAD_TRUE.EQ.0)THEN
                
                UTEMP=ZERO
                num_valid_neighbrs = 0

                DO L=1,IELEM(N,I)%IFCA	!faces2
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N) THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0) THEN	!CHECK FOR BOUNDARIES
                            IF (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5) THEN	!PERIODIC IN MY CPU
                                helper = U_C(IELEM(N,I)%INEIGH(L))%VAL(4,rho_index)
                                if ((helper.eq.helper).and.(helper.gt.0)) then
                                    helper = U_C(IELEM(N,I)%INEIGH(L))%VAL(4,p_index)
                                    if ((helper.eq.helper).and.(helper.gt.0)) then
                                        num_valid_neighbrs = num_valid_neighbrs + 1
                                        UTEMP(num_valid_neighbrs,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(4,1:NOF_VARIABLES)
                                    end if
                                end if
                            ELSE
                                !NOT PERIODIC ONES IN MY CPU			  				  
                            END IF
                        ELSE
                            helper = U_C(IELEM(N,I)%INEIGH(L))%VAL(4,rho_index)
                            if ((helper.eq.helper).and.(helper.gt.0)) then
                                helper = U_C(IELEM(N,I)%INEIGH(L))%VAL(4,p_index)
                                if ((helper.eq.helper).and.(helper.gt.0)) then
                                    num_valid_neighbrs = num_valid_neighbrs + 1
                                    UTEMP(num_valid_neighbrs,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(4,1:NOF_VARIABLES)
                                end if
                            end if
                        END IF
                    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                helper =IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),rho_index)
                                if ((helper.eq.helper).and.(helper.gt.0.0)) then
                                    helper = IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),p_index)
                                    if ((helper.eq.helper).and.(helper.gt.0.0)) then
                                        num_valid_neighbrs = num_valid_neighbrs+1
                                        UTEMP(num_valid_neighbrs,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2&
                                            (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                                    end if
                                end if 
                            END IF
                        ELSE
                            helper = IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),rho_index)
                            if ((helper.eq.helper).and.(helper.gt.0.0)) then
                                helper = IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),p_index)
                                if ((helper.eq.helper).and.(helper.gt.0.0)) then
                                    num_valid_neighbrs = num_valid_neighbrs + 1
                                    UTEMP(num_valid_neighbrs,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2&
                                        (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                                end if
                            end if
                        END IF
                    END IF
                END DO

                if (num_valid_neighbrs.ge.1) then
                                    
                    !5 NOW ESTABLISH THE MIN AND MAX BOUNDS
                    DO IEX=1,NOF_VARIABLES
                        UTMIN(IEX)=MINVAL(UTEMP(1:num_valid_neighbrs,IEX))
                        UTMAX(IEX)=MAXVAL(UTEMP(1:num_valid_neighbrs,IEX))
                    END DO
                        
                    ! MOOD_VAR1=0.0001;MOOD_VAR2=0.001; MOOD_MODE=RELAXED
                    
                    DO IEX=1,NOF_VARIABLES
                        NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                    END DO

                    NAD_TRUE=0

                    DO IEX=1,NOF_VARIABLES
                        IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
                            IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                                NAD_TRUE=1
                            END IF
                        endif
                    END DO
                else
                    NAD_TRUE = 1
                end if
                                
            END IF
                
            IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE

        END DO 
        !$OMP END DO

    CASE (4)

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
            
            IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                PAD_TRUE=1
            END IF

            IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN						
                PAD_TRUE=1
            END IF
        
            !3 THE ONES WITH A PHYSICALLY ADMISSIBLE SOLUTION NEED TO GET CHECKED
            IF (PAD_TRUE.EQ.0)THEN
                
                !4 NOW ESTABLISH A TEMPORARY ARRAY WITH THE CURRENT SOLUTION FROM THE DIRECT SIDE NEIGHBOURS OF CONSIDERED CELL
                UTEMP=ZERO
                num_valid_neighbrs = 0
                DO L=1,IELEM(N,I)%IFCA
                    num_valid_neighbrs = num_valid_neighbrs + 1
                    UTEMP(num_valid_neighbrs,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(4,1:NOF_VARIABLES)
                END DO
                           
                !5 NOW ESTABLISH THE MIN AND MAX BOUNDS
                DO IEX=1,NOF_VARIABLES
                    UTMIN(IEX)=MINVAL(UTEMP(1:num_valid_neighbrs,IEX))
                    UTMAX(IEX)=MAXVAL(UTEMP(1:num_valid_neighbrs,IEX))
                END DO
                    
                ! MOOD_VAR1=0.0001;MOOD_VAR2=0.001; MOOD_MODE=RELAXED
                
                DO IEX=1,NOF_VARIABLES
                    NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                END DO

                NAD_TRUE=0

                DO IEX=1,NOF_VARIABLES
                    if ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
                        IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                            NAD_TRUE=1
                        END IF
                    endif
                END DO
                
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

            IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index))) THEN						
                PAD_TRUE=1
            END IF

            IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN						
                PAD_TRUE=1
            END IF
                
            IF (PAD_TRUE.EQ.0)THEN
                
                UTEMP=ZERO
                num_valid_neighbrs = 0

                DO L=1,IELEM(N,I)%IFCA	!faces2
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N) THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0) THEN	!CHECK FOR BOUNDARIES
                            IF (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5) THEN	!PERIODIC IN MY CPU
                                num_valid_neighbrs = num_valid_neighbrs + 1
                                UTEMP(num_valid_neighbrs,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(4,1:NOF_VARIABLES)
                            ELSE
                                !NOT PERIODIC ONES IN MY CPU			  				  
                            END IF
                        ELSE    
                            num_valid_neighbrs = num_valid_neighbrs + 1
                            UTEMP(num_valid_neighbrs,1:NOF_VARIABLES)=U_C(IELEM(N,I)%INEIGH(L))%VAL(4,1:NOF_VARIABLES)            
                        END IF
                    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                num_valid_neighbrs = num_valid_neighbrs+1
                                UTEMP(num_valid_neighbrs,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2&
                                    (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables) 
                            end if
                        ELSE
                            num_valid_neighbrs = num_valid_neighbrs + 1
                            UTEMP(num_valid_neighbrs,1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL2&
                                (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                        END IF
                    END IF
                END DO

                if (num_valid_neighbrs.ge.1) then
                                    
                    !5 NOW ESTABLISH THE MIN AND MAX BOUNDS
                    DO IEX=1,NOF_VARIABLES
                        UTMIN(IEX)=MINVAL(UTEMP(1:num_valid_neighbrs,IEX))
                        UTMAX(IEX)=MAXVAL(UTEMP(1:num_valid_neighbrs,IEX))
                    END DO
                        
                    ! MOOD_VAR1=0.0001;MOOD_VAR2=0.001; MOOD_MODE=RELAXED
                    
                    DO IEX=1,NOF_VARIABLES
                        NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(UTMAX(IEX)-UTMIN(IEX)))
                    END DO

                    NAD_TRUE=0

                    DO IEX=1,NOF_VARIABLES
                        if ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
                            IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                                NAD_TRUE=1
                            END IF
                        endif
                    END DO
                else
                    NAD_TRUE = 1
                end if
                                
            END IF
                
            IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE

        END DO 
        !$OMP END DO

    CASE (5)

        IF (ITESTCASE.GE.1) THEN

            if (CASCADE.eq.1) then
                parameter = MOOD_VAR5
            else
                parameter = MOOD_VAR6
            end if
            
            !$OMP DO
            DO II=1,NOF_INTERIOR
                I=EL_INT(II)
                ICONSIDERED=I
                IELEM(N,I)%MOOD=0
                REDUCE1=0
                PAD_TRUE=0
                NAD_TRUE=0

                CELL_AREA=IELEM(N,I)%TOTVOLUME
                if (dimensiona.eq.2) then
                    CELL_SIZE = SQRT(CELL_AREA)
                else
                    CELL_SIZE = CELL_AREA ** (1/3)
                end if
                ! print*,CELL_SIZE
                
                !1 copy candidate solution at temp variable   
                LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)

                !2 TRANSFORM CONSERVATIVE  TO PRIMITIVE AND CHECK IF PRESSURE AND DENSITY ARE PHYSICALLY ADMISSIBLE IF NOT PAD_TRUE=1
                IF (ITESTCASE.GE.3) THEN ! for advection equation there is no need to transfor to primitive
                    CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
                END IF
                
                IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN
                    IF ((LEFTV(rho_index).LE.(-0.001)).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                ELSE
                    IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                    IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
                        PAD_TRUE=1
                    END IF
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
                        Delta(IEX) = UTMAX(IEX)-UTMIN(IEX)
                        NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*(Delta(IEX)))
                    END DO

                    NAD_TRUE=0

                    IF (DELTA(rho_index)/CELL_SIZE > parameter*U_C(I)%VAL(4,rho_index)) THEN
                        NAD_TRUE = 1
                    END IF
                    IF (DELTA(p_index)/CELL_SIZE > parameter*U_C(I)%VAL(4,p_index)) THEN
                        NAD_TRUE = 1
                    END IF
                    ! IF (DELTA(rho_index)/CELL_SIZE > parameter) THEN
                    !     NAD_TRUE = 1
                    ! END IF
                    ! IF (DELTA(p_index)/CELL_SIZE > parameter) THEN
                    !     NAD_TRUE = 1
                    ! END IF

                    !6 SPECIFY RELAXED OR ORIGINAL MOOD PATTERN
                    IF (NAD_TRUE.EQ.0)THEN
                        DO IEX=1,NOF_VARIABLES
                            IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index))then
                                IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX))))THEN
                                    NAD_TRUE=1
                                END IF
                            endif
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

                IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN
                    IF ((LEFTV(rho_index).LE.(-0.001)).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                ELSE
                    IF ((LEFTV(rho_index).LE.ZERO).OR.(LEFTV(rho_index).NE.LEFTV(rho_index)))THEN						
                        PAD_TRUE=1
                    END IF
                    IF ((LEFTV(p_index).LE.ZERO).OR.(LEFTV(p_index).NE.LEFTV(p_index))) THEN					
                        PAD_TRUE=1
                    END IF
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
                        Delta(IEX) = UTMAX(IEX)-UTMIN(IEX)
                        NAD_DELTA1(IEX)=MAX(RELAX_MOOD1,(RELAX_MOOD2)*Delta(IEX))
                    END DO
                    
                    NAD_TRUE=0

                    IF (DELTA(rho_index)/CELL_SIZE > parameter*U_C(I)%VAL(4,rho_index)) THEN
                        NAD_TRUE = 1
                    END IF
                    IF (DELTA(p_index)/CELL_SIZE > parameter*U_C(I)%VAL(4,p_index)) THEN
                        NAD_TRUE = 1
                    END IF
                    ! IF (DELTA(rho_index)/CELL_SIZE > parameter) THEN
                    !     NAD_TRUE = 1
                    ! END IF
                    ! IF (DELTA(p_index)/CELL_SIZE > parameter) THEN
                    !     NAD_TRUE = 1
                    ! END IF
                                    
                    IF (NAD_TRUE.EQ.0) THEN
                        DO IEX=1,NOF_VARIABLES
                            IF ((IEX.EQ.rho_index).OR.(IEX.EQ.p_index)) THEN
                                IF ((U_C(I)%VAL(4,IEX).LT.(UTMIN(IEX)-NAD_DELTA1(IEX))).OR.(U_C(I)%VAL(4,IEX).GT.(UTMAX(IEX)+NAD_DELTA1(IEX)))) THEN
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
    
    CASE (6)

        IF (ITESTCASE.EQ.3) THEN

            if (CASCADE.eq.1) then
                parameter = MOOD_VAR5
            else
                parameter = MOOD_VAR6
            end if
            
            !$OMP DO REDUCTION (MAX:MAX_ENTROPY)
            DO II=1,NOF_INTERIOR
                I=EL_INT(II)
                ICONSIDERED=I
                IELEM(N,I)%MOOD=0
                REDUCE1=0
                PAD_TRUE=0
                NAD_TRUE=0

                CELL_AREA=IELEM(N,I)%TOTVOLUME
                if (dimensiona.eq.2) then
                    CELL_SIZE = SQRT(CELL_AREA)
                else
                    CELL_SIZE = CELL_AREA ** (1/3)
                end if
                ! print*,CELL_SIZE
                
                !1 copy candidate solution at temp variable   
                LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)

                !2 TRANSFORM CONSERVATIVE  TO PRIMITIVE AND CHECK IF PRESSURE AND DENSITY ARE PHYSICALLY ADMISSIBLE IF NOT PAD_TRUE=1
                CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
                rho_new = leftv(rho_index)
                p_new = leftv(p_index)
                
                IF ((rho_new.LE.ZERO).OR.(rho_new.NE.rho_new))THEN						
                    PAD_TRUE=1
                END IF
                IF ((p_new.LE.ZERO).OR.(p_new.NE.p_new)) THEN					
                    PAD_TRUE=1
                END IF
            
                !3 THE ONES WITH A PHYSICALLY ADMISSIBLE SOLUTION NEED TO GET CHECKED
                IF (PAD_TRUE.EQ.0)THEN

                    !4 copy old solution at temp variable   
                    LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(3,1:NOF_VARIABLES)

                    !5 TRANSFORM CONSERVATIVE  TO PRIMITIVE AND COMPUTE ENTROPIES
                    CALL CONS2PRIM(N,leftv,MP_PINFl,GAMMA)
                    rho_old = leftv(rho_index)
                    p_old = leftv(p_index)
                    
                    entropy_new = ENTROPY(rho_new, p_new, gamma)
                    entropy_old = ENTROPY(rho_old, p_old, gamma)

                    max_entropy = MAX(entropy_old, max_entropy)
                    
                    ! helper_value = ABS(entropy_new - entropy_old)/CELL_AREA
                    helper_value = ABS(entropy_new - entropy_old)/CELL_SIZE
                    ! helper_value = ABS(entropy_new - entropy_old)
                    ! write(*,*) helper_value
                    IF (helper_value > parameter) THEN
                        NAD_TRUE = 1
                    END IF
                        
                END IF

                !7 NOW SET THE MOOD FLAG FOR EACH ELEMENT
                IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE
                    
            END DO
            !$OMP END DO		

            !$OMP DO REDUCTION (MAX:MAX_ENTROPY)
            DO II=1,NOF_BOUNDED
                I=EL_BND(II)
                ICONSIDERED=I
                REDUCE1=0
                PAD_TRUE=0
                NAD_TRUE=0
                IELEM(N,I)%MOOD=0
                LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)

                CELL_AREA=IELEM(N,I)%TOTVOLUME
                if (dimensiona.eq.2) then
                    CELL_SIZE = SQRT(CELL_AREA)
                else
                    CELL_SIZE = CELL_AREA ** (1/3)
                end if

                IF (DIMENSIONA.EQ.3) THEN
                    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
                ELSE
                    CALL cons2prim(N,leftv,MP_PINFl,gammal)                   
                END IF
                rho_new = leftv(rho_index)
                p_new = leftv(p_index)

                IF ((rho_new.LE.ZERO).OR.(rho_new.NE.rho_new))THEN						
                    PAD_TRUE=1
                END IF
                IF ((p_new.LE.ZERO).OR.(p_new.NE.p_new)) THEN					
                    PAD_TRUE=1
                END IF
                    
                IF (PAD_TRUE.EQ.0)THEN
                    
                    LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(3,1:NOF_VARIABLES)

                    !5 TRANSFORM CONSERVATIVE  TO PRIMITIVE AND COMPUTE ENTROPIES
                    CALL CONS2PRIM(N,leftv,MP_PINFl,gamma)
                    rho_old = leftv(rho_index)
                    p_old = leftv(p_index)
                    
                    entropy_new = ENTROPY(rho_new, p_new, gamma)
                    entropy_old = ENTROPY(rho_old, p_old, gamma)

                    max_entropy = MAX(entropy_old, max_entropy)
                    
                    ! helper_value = ABS(entropy_new - entropy_old)/CELL_AREA
                    helper_value = ABS(entropy_new - entropy_old)/CELL_SIZE
                    ! helper_value = ABS(entropy_new - entropy_old)
                    ! write(*,*) helper_value
                    IF (helper_value > parameter) THEN
                        NAD_TRUE = 1
                    END IF        
                END IF
                    
                IELEM(N,I)%MOOD=NAD_TRUE+PAD_TRUE

            END DO 
            !$OMP END DO
        ELSE
            write(*,*) "Unsuported combination of MOOD trigger and test case"
            call exit(777)
        END IF		
    

      CASE DEFAULT
        
        if (n.eq.0) then
            write(*,*) 'reached unreachable case in PAD_NAD function'
            call exit(777)
        end if

    END SELECT
		
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
