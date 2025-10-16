MODULE hybridCWENO_MOOD_module
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
USE MOODR

IMPLICIT NONE

CONTAINS

SUBROUTINE hybridCWENO_MOOD_OPERATOR(N)
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

    CASCADE = 1
    CALL PAD_NAD(N)

    CALL EXHBOUNDHIGHER_MOOD(N)

    CALL FIX_LIST(N)

    ! CALL MUSCL(N)
    CALL WENOWEIGHTS_hybrid(N)
    ! CALL CHECKSOL(N)
    ! CALL MUSCL(N)
    ! CALL CHECKSOLX(N)

    CALL EXHBOUNDHIGHER(N)

    IF (dimensiona.eq.2) THEN
        IF ((ITESTCASE.EQ.1).OR.(ITESTCASE.EQ.2)) THEN ! linear advection equation
            CALL CALCULATE_FLUXESHI2D_MOOD(N)
        ELSE ! Euler or Navier-Stokes
            CALL CALCULATE_FLUXESHI_CONVECTIVE2D_MOOD(N)
        END IF
    ELSE
        if (n.eq.0) then
            write (*,*) "This mode does not support 3D simulations (yet)"
        end if
        call abort
    END IF

    !$OMP DO
    DO I=1,KMAXE
        IF (IELEM(N,I)%MOOD.EQ.0) THEN
            IELEM(N,I)%MOOD_O=IORDER+1  !TARGET POLYNOMIAL/SCHEME ADMISSIBLE
        ELSE
            IELEM(N,I)%MOOD_O=1         !SECOND ORDER MUSCL ADMISSIBLE
        END IF
    END DO
    !$OMP END DO

END SUBROUTINE hybridCWENO_MOOD_OPERATOR


SUBROUTINE RUNGE_KUTTA3_2D_hybridCWENO_MOOD(N)
    !> @brief
    !> SSP RUNGE KUTTA 3RD-ORDER SCHEME
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I,KMAXE,INDS
    REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
    KMAXE=XMPIELRANK(N)
    TO4=3.0D0/4.0D0
    OO4=1.0D0/4.0D0
    TO3=2.0D0/3.0D0
    OO3=1.0D0/3.0D0	

    INDS=4

    ! IWENO = 0 ! TODO find more elegent solution

    IF (FASTEST.EQ.1)THEN
        CALL EXCHANGE_LOWER(N)
        CALL ARBITRARY_ORDER(N)
        CALL EXHBOUNDHIGHER(N)
        
        SELECT CASE(ITESTCASE)
          CASE(1,2)
            CALL CALCULATE_FLUXESHI2D(N)
          CASE(3)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
          CASE(4)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
            CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
            IF (turbulence.eq.1)THEN
                CALL SOURCES_COMPUTATION2d(N)
            END IF
        END SELECT
    ELSE
        CALL EXCHANGE_HIGHER(N)
        CALL ARBITRARY_ORDER(N)
        CALL EXHBOUNDHIGHER(N)
        SELECT CASE(ITESTCASE)
          CASE(1,2)
            CALL CALCULATE_FLUXESHI2D(N)
          CASE(3)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
          CASE(4)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
            CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
            IF (turbulence.eq.1)THEN
                CALL SOURCES_COMPUTATION2d(N)
            END IF
        END SELECT
    END IF

    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)

        U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)

        U_C(I)%VAL(INDS,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
    END DO
    !$OMP END DO

    IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        !$OMP DO
        DO I=1,KMAXE
            OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
            U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
            U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
        END DO
        !$OMP END DO
    END IF

    ! IWENO = 1 ! TODO more elegent solution
    CALL hybridCWENO_MOOD_OPERATOR(N)
    ! IWENO = 0 ! TODO more elegent solution

    !$OMP DO
    DO I=1,KMAXE
        IF (IELEM(N,I)%RECALC.EQ.1)THEN
            OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
            U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
            IELEM(N,I)%MOOD_O=1
        ELSE
            U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
        END IF

        IELEM(N,I)%RECALC = 0
        IELEM(N,I)%MOOD = 0
    END DO
    !$OMP END DO

    IF (FASTEST.EQ.1)THEN
        CALL EXCHANGE_LOWER(N)
        CALL ARBITRARY_ORDER(N)
        CALL EXHBOUNDHIGHER(N)
        SELECT CASE(ITESTCASE)
          CASE(1,2)
            CALL CALCULATE_FLUXESHI2D(N)
          CASE(3)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
          CASE(4)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
            CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
            IF (turbulence.eq.1)THEN
                CALL SOURCES_COMPUTATION2d(N)
            END IF
        END SELECT
    ELSE
        CALL EXCHANGE_HIGHER(N)
        CALL ARBITRARY_ORDER(N)
        CALL EXHBOUNDHIGHER(N)
        SELECT CASE(ITESTCASE)
          CASE(1,2)
              CALL CALCULATE_FLUXESHI2D(N)
          CASE(3)
              CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
          CASE(4)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
            CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
            IF (turbulence.eq.1)THEN
                CALL SOURCES_COMPUTATION2d(N)
            END IF
        END SELECT
    END IF

    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(inds,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
            ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    END DO
    !$OMP END DO

    IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
      !$OMP DO
      DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
        U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(TO4*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(OO4*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((OO4))*((DT)*&
          ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
      END DO
      !$OMP END DO
    END IF

    ! IWENO = 1 ! TODO more elegent solution
    CALL hybridCWENO_MOOD_OPERATOR(N)
    ! IWENO = 0 ! TODO more elegent solution

    !$OMP DO
    DO I=1,KMAXE
        IF (IELEM(N,I)%RECALC.EQ.1)THEN
            OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
            U_C(I)%VAL(1,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
                ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
            IELEM(N,I)%MOOD_O=1
        ELSE
            U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
        END IF

        IELEM(N,I)%RECALC = 0
        IELEM(N,I)%MOOD = 0
    END DO
    !$OMP END DO
 
    IF (FASTEST.EQ.1)THEN
        CALL EXCHANGE_LOWER(N)
        CALL ARBITRARY_ORDER(N)
        CALL EXHBOUNDHIGHER(N)
        SELECT CASE(ITESTCASE)
          CASE(1,2)
            CALL CALCULATE_FLUXESHI2D(N)
          CASE(3)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
          CASE(4)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
            CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
            IF (turbulence.eq.1)THEN
                CALL SOURCES_COMPUTATION2d(N)
            END IF
        END SELECT
    ELSE
        CALL EXCHANGE_HIGHER(N)
        CALL ARBITRARY_ORDER(N)
        CALL EXHBOUNDHIGHER(N)
        SELECT CASE(ITESTCASE)
          CASE(1,2)
            CALL CALCULATE_FLUXESHI2D(N)
          CASE(3)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
          CASE(4)
            CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
            CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
            IF (turbulence.eq.1)THEN
                CALL SOURCES_COMPUTATION2d(N)
            END IF
            CALL VORTEXCALC2D(N)
        END SELECT
    END IF

    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(INDS,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
            ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    END DO
    !$OMP END DO

    ! IWENO = 1 ! TODO more elegent solution
    CALL hybridCWENO_MOOD_OPERATOR(N)
    ! IWENO = 0 ! TODO more elegent solution
    
    !$OMP DO
    DO I=1,KMAXE
        IF (IELEM(N,I)%RECALC.EQ.1) THEN
            OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
            U_C(I)%VAL(1,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
                ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
            IELEM(N,I)%MOOD_O=1
        ELSE
            U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
        END IF
        
        IELEM(N,I)%RECALC = 0
        IELEM(N,I)%MOOD = 0
    END DO
    !$OMP END DO

    IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        !$OMP DO
        DO I=1,KMAXE
            OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
            U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=((OO3)*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+((TO3)*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(((TO3))*&
                ((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
        END DO
        !$OMP END DO
    END IF

    ! IF (AVERAGING.EQ.1)THEN
    !     if (n.eq.0) write(*,*) "Averaging not supported in this mode (yet)"
    !     CALL ABORT
    ! END IF

END SUBROUTINE RUNGE_KUTTA3_2D_hybridCWENO_MOOD

END MODULE hybridCWENO_MOOD_module
