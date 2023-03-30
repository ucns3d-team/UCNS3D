MODULE DG_FUNCTIONS

USE BASIS
USE DECLARATION
USE DERIVATIVES
USE LAPCK

IMPLICIT NONE

CONTAINS




FUNCTION DG_SOL(N)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
REAL,DIMENSION(1:NUMBER_OF_DOG)::BASIS_TEMP
INTEGER,INTENT(IN)::N
INTEGER::I_DOF, I_VAR
REAL,DIMENSION(1:Nof_VARIABLES)::DG_SOL
REAL,external:: ddot

 compwrt=-2
NUMBER=IELEM(N,ICONSIDERED)%IORDER




if (dimensiona.eq.2)then
BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
else
BASIS_TEMP = BASIS_REC(N, X1, Y1,z1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
end if






DO I_VAR = 1, NOF_VARIABLES
    DG_SOL(I_VAR) = U_C(ICONSIDERED)%VALDG(1,I_VAR,1)+ DDOT(NUMBER_OF_DOG,BASIS_TEMP(1:NUMBER_OF_DOG),1,U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUMBER_OF_DOG+1),1)

END DO

 compwrt=0


END FUNCTION DG_SOL


FUNCTION DG_SOLFACE(N)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given surface point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
INTEGER::I_DOF, I_VAR
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:nof_Variables)::DG_SOLFACE
REAL,external:: ddot

 compwrt=-2

x1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,1)
 y1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,2)
if (dimensiona.eq.3)then
 z1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,3)
end if







NUMBER=IELEM(N,ICONSIDERED)%IORDER
if (dimensiona.eq.2)then
BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
else
BASIS_TEMP = BASIS_REC(N, X1, Y1,z1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
end if



DO I_VAR = 1, NOF_VARIABLES
    DG_SOLFACE(I_VAR) = U_C(ICONSIDERED)%VALDG(1,I_VAR,1) + DDOT(NUMBER_OF_DOG,BASIS_TEMP(1:NUMBER_OF_DOG),1,U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUMBER_OF_DOG+1),1)
END DO





    compwrt=0

END FUNCTION DG_SOLFACE


FUNCTION DG_SOL_DER()
IMPLICIT NONE
!> @brief
!> This function returns the derivative of the DG solution at a given point (X1, Y1, [Z1])
    INTEGER::I_DOF, I_VAR, I_DIM
    REAL,DIMENSION(NUMBER_OF_DOG, DIMENSIONA)::BASIS_TEMP
    REAL,DIMENSION(NOF_VARIABLES, DIMENSIONA)::DG_SOL_DER
    REAL,external:: DDOT

    IF (BR2_YN /= 1) THEN
        DO I_DOF = 1, NUMBER_OF_DOG
            if (dimensiona.eq.2) then
                IF( POLY == 1) THEN
                    BASIS_TEMP(I_DOF,1) = DF2DX(X1,Y1,I_DOF)
                    BASIS_TEMP(I_DOF,2) = DF2DY(X1,Y1,I_DOF)
                ELSE IF( POLY == 4) THEN
                    BASIS_TEMP(I_DOF,1) = TL2DX(X1,Y1,I_DOF)
                    BASIS_TEMP(I_DOF,2) = TL2DY(X1,Y1,I_DOF)
                END IF
            ELSE
                BASIS_TEMP(I_DOF,1) = TL3DX(X1,Y1,Z1,I_DOF)
                BASIS_TEMP(I_DOF,2) = TL3DY(X1,Y1,Z1,I_DOF)
                BASIS_TEMP(I_DOF,3) = TL3DZ(X1,Y1,Z1,I_DOF)
            end if
        END DO
    END IF

    compwrt=-2
    DO I_VAR = 1, NOF_VARIABLES
        DO I_DIM = 1, DIMENSIONA
            IF (BR2_YN == 1) THEN
                IF(DIMENSIONA.EQ.2)THEN
                    DG_SOL_DER(I_VAR, I_DIM) = U_C(ICONSIDERED)%BR2_AUX_VAR(1,I_VAR,I_DIM) + DDOT(NUMBER_OF_DOG, BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG),1, U_C(ICONSIDERED)%BR2_AUX_VAR(2:NUM_DG_DOFS,I_VAR,I_DIM),1)
                ELSE
                    DG_SOL_DER(I_VAR, I_DIM) = U_C(ICONSIDERED)%BR2_AUX_VAR(1,I_VAR,I_DIM) + DDOT(NUMBER_OF_DOG, BASIS_REC(N, X1, Y1, Z1, NUMBER, ICONSIDERED, NUMBER_OF_DOG),1, U_C(ICONSIDERED)%BR2_AUX_VAR(2:NUM_DG_DOFS,I_VAR,I_DIM),1)
                 END IF
            ELSE
                DG_SOL_DER(I_VAR, I_DIM) = DDOT(NUMBER_OF_DOG, BASIS_TEMP(1:NUMBER_OF_DOG,I_DIM),1, U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUM_DG_DOFS),1)
            END IF
        END DO
    END DO
    compwrt=0
END FUNCTION DG_SOL_DER


FUNCTION BR2_LOCAL_LIFT(N, N_QP)
IMPLICIT NONE
!> REQUIRES: ICONSIDERED, FACEX, and WEIGHTS_DG or WEQUA2D
INTEGER,INTENT(IN)::N, N_QP
INTEGER::I, L, NGP, IQP, I_DIM
REAL,DIMENSION(DIMENSIONA)::NNN
REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP
REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
REAL,DIMENSION(NOF_VARIABLES,DIMENSIONA)::BR2_LOCAL_LIFT

    I = ICONSIDERED
    L = FACEX

	IF (DIMENSIONA == 3) THEN
        ANGLE1=IELEM(N,I)%FACEANGLEX(L)
        ANGLE2=IELEM(N,I)%FACEANGLEY(L)
        NNN(1)=(COS(ANGLE1)*SIN(ANGLE2))
        NNN(2)=(SIN(ANGLE1)*SIN(ANGLE2))
        NNN(3)=(COS(ANGLE2))
    ELSE
        NNN(1)=IELEM(N,I)%FACEANGLEX(L)
        NNN(2)=IELEM(N,I)%FACEANGLEY(L)
    END IF

    BR2_LOCAL_LIFT = 0.0D0
    DO NGP = 1, N_QP

        POINTX = NGP

        COMPWRT=-2
        X1 = ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,1)
        Y1 = ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,2)
        IF (DIMENSIONA == 3) Z1 = ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,3)

        IF (DIMENSIONA .EQ. 2)THEN
            BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
        ELSE
            BASIS_TEMP = BASIS_REC(N, X1, Y1, Z1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
        END IF

        CALL GET_LEFT_RIGHT_STATES(ICONSIDERED,FACEX,POINTX)

        LEFTV = CLEFT
        RIGHTV = CRIGHT
        IF (DIMENSIONA == 3) THEN
            CALL CONS2PRIM2(N)
        ELSE
            CALL CONS2PRIM2D2(N)
        END IF
        ! Change pressure to temperature
        LEFTV(NOF_VARIABLES) = LEFTV(NOF_VARIABLES) / LEFTV(1)
        RIGHTV(NOF_VARIABLES) = RIGHTV(NOF_VARIABLES) / RIGHTV(1)

        DO I_DIM = 1, DIMENSIONA
!             RHLLCFLUX = NNN(I_DIM) * (RIGHTV - LEFTV)
            BR2_LOCAL_LIFT(1:NOF_VARIABLES,I_DIM) = BR2_LOCAL_LIFT(1:NOF_VARIABLES,I_DIM) + NNN(I_DIM) * (RIGHTV - LEFTV) * WEQUA2D(POINTX) !SUM(DG_SURF_FLUX(N),DIM=1)
!             IF( N == 3 .AND. ICONSIDERED == 20) WRITE(400+N,*) BR2_LOCAL_LIFT, "BR2_LOCAL_LIFT", DG_SURF_FLUX(N), "DG_SURF_FLUX", SUM(DG_SURF_FLUX(N),DIM=1), "SUM(DG_SURF_FLUX(N),DIM=1)"
        END DO

    END DO
    BR2_LOCAL_LIFT = OO2 * BR2_LOCAL_LIFT * IELEM(N,ICONSIDERED)%SURF(FACEX) / IELEM(N,ICONSIDERED)%TOTVOLUME ! / DG_VOL_INTEGRAL3(N)

END FUNCTION BR2_LOCAL_LIFT



FUNCTION DG_SURF_FLUX(N)
!> @brief
!> Calculates the RHS flux term to be integrated in the DG formulation
IMPLICIT NONE
REAL,DIMENSION(NUMBER_OF_DOG+1,NOF_VARIABLES)::DG_SURF_FLUX
INTEGER::I
INTEGER,INTENT(IN)::N
 compwrt=-2
x1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,1)
 y1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,2)
if (dimensiona.eq.3)then
 z1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,3)
end if



   

NUMBER=IELEM(N,ICONSIDERED)%IORDER
if (dimensiona.eq.2)Then
DO I = 1, NOF_VARIABLES
    DG_SURF_FLUX(1,I) = rHLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)
    
    DG_SURF_FLUX(2:NUMBER_OF_DOG+1,I) = rHLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC2D(N,X1,Y1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
    
END DO

else
DO I = 1, NOF_VARIABLES
    DG_SURF_FLUX(1,I) = rHLLCFLUX(I) * WEIGHTS_dg(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)
    
    DG_SURF_FLUX(2:NUMBER_OF_DOG+1,I) = rHLLCFLUX(I) * WEIGHTS_dg(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC(N,X1,Y1,z1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
END DO

end if



    compwrt=0




END FUNCTION DG_SURF_FLUX


FUNCTION DG_SURF_FLUXV(N)
!> @brief
!> Calculates the RHS flux term to be integrated in the DG formulation
IMPLICIT NONE
REAL,DIMENSION(NUMBER_OF_DOG+1,1:NOF_VARIABLES)::DG_SURF_FLUXV
INTEGER::I
INTEGER,INTENT(IN)::N
 compwrt=-2
x1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,1)
 y1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,2)
if (dimensiona.eq.3)then
 z1= ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX,POINTX,3)
end if





NUMBER=IELEM(N,ICONSIDERED)%IORDER
if (dimensiona.eq.2)Then
DO I = 1, NOF_VARIABLES
    DG_SURF_FLUXV(1,I) = HLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)

    DG_SURF_FLUXV(2:NUMBER_OF_DOG+1,I) = HLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC2D(N,X1,Y1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)

END DO

else
DO I = 1, NOF_VARIABLES
    DG_SURF_FLUXV(1,I) = HLLCFLUX(I) * WEIGHTS_dg(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)

    DG_SURF_FLUXV(2:NUMBER_OF_DOG+1,I) = HLLCFLUX(I) * WEIGHTS_dg(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC(N,X1,Y1,z1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
END DO

end if





    compwrt=0





END FUNCTION DG_SURF_FLUXV




FUNCTION DG_VOL_INTEGRAL(N)
!> @brief
!> Calculates the volume integral term in the DG RHS for scalar linear advection with speed = 1
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,NQP,I_QP
REAL,DIMENSION(IDEGFREE+1,NOF_VARIABLES)::DG_VOL_INTEGRAL



    compwrt=-2



    NQP = IELEM(N,ICONSIDERED)%ITOTALPOINTS


    NUMBER=IELEM(N,ICONSIDERED)%IORDER
    NUMBER_OF_DOG = IELEM(N,ICONSIDERED)%IDEGFREE


    DG_VOL_INTEGRAL(:,:) = 0.0D0




    DO I=1,NUMBER_OF_DOG

        DO I_QP = 1, NQP

            X1=QP_ARRAY(ICONSIDERED)%X(I_QP)
            Y1=QP_ARRAY(ICONSIDERED)%Y(I_QP)
            IF (DIMENSIONA.EQ.3)THEN
            Z1=QP_ARRAY(ICONSIDERED)%Z(I_QP)
            END IF

            if (itestcase.lt.3) then ! Linear advection
                FLUX_TERM_X=dg_sol(n)*LAMX !flux in x-axis of linear advection in 2d (lamx*sol), where lamx is the wave speed for x axis, and sol is the solution
                FLUX_TERM_Y=dg_sol(n)*LAMY !flux in y-axis of linear advection in 2d (lamy*sol), where lamy is the wave speed for y axis, and sol is the solution

                if (dimensiona.eq.3) FLUX_TERM_z = dg_sol(n)*LAMz

            ELSE IF (ITESTCASE.EQ.3) THEN ! Euler

                LEFTV=DG_SOL(N)




                IF (DIMENSIONA.EQ.2)THEN
                    CALL CONS2PRIM2D(N)

                    CALL FLUX2DX
                    CALL FLUX2DY

                ELSE

                    CALL CONS2PRIM(N)
                    CALL FLUX3DX
                    CALL FLUX3DY
                    CALL FLUX3DZ

                END IF

            ELSE IF (ITESTCASE == 4) THEN ! NS
                LEFTV = DG_SOL(N)
                LEFTV_DER = DG_SOL_DER()

                IF( BR2_YN == 2) THEN
                    CALL DCONS2DPRIM
                    LEFTV_DER = LEFTV_DER + SUM(ILOCAL_RECON3(ICONSIDERED)%BR2_LOCAL_LIFT, DIM=3)
                END IF


                IF (DIMENSIONA.EQ.2)THEN
                    CALL CONS2PRIM2D(N)
                    CALL FLUX2DX
                    CALL FLUX2DY

                    CALL FLUX_VISC2D(N)

                ELSE
                    CALL CONS2PRIM(N)
                    CALL FLUX3DX
                    CALL FLUX3DY
                    CALL FLUX3DZ
                    CALL FLUX_VISC3D(N)
                END IF




            END IF


                    if (dimensiona.eq.2)then


                    if (poly.eq.1)then

                    DG_VOL_INTEGRAL(I+1,:) = DG_VOL_INTEGRAL(I+1,:)+ QP_ARRAY(ICONSIDERED)%QP_WEIGHT(I_QP) *(FLUX_TERM_X(:)*DF2DX(X1,Y1,I)+FLUX_TERM_Y(:)*DF2DY(X1,Y1,I))

                    end if

                    if (poly.eq.4)then

                    DG_VOL_INTEGRAL(I+1,:) = DG_VOL_INTEGRAL(I+1,:)+ QP_ARRAY(ICONSIDERED)%QP_WEIGHT(I_QP)*(FLUX_TERM_X(:)*TL2DX(X1,Y1,I)+FLUX_TERM_Y(:)*TL2DY(X1,Y1,I))

                    end if


                    else
                            if (poly.eq.1)then

                            DG_VOL_INTEGRAL(I+1,:) = DG_VOL_INTEGRAL(I+1,:)+ QP_ARRAY(ICONSIDERED)%QP_WEIGHT(I_QP) *((FLUX_TERM_X(:)*DFX(X1,Y1,z1,I))&
                            +(FLUX_TERM_Y(:)*DFY(X1,Y1,z1,I))+(FLUX_TERM_z(:)*DFZ(X1,Y1,z1,I)))

                            END IF

                            if (poly.eq.4)then

                            DG_VOL_INTEGRAL(I+1,:) = DG_VOL_INTEGRAL(I+1,:)+ QP_ARRAY(ICONSIDERED)%QP_WEIGHT(I_QP)*((FLUX_TERM_X(:)*TL3DX(X1,Y1,z1,I))&
                            +(FLUX_TERM_Y(:)*TL3DY(X1,Y1,z1,I))+(FLUX_TERM_z(:)*TL3DZ(X1,Y1,z1,I)))

                            END IF

                    END IF
         END DO
     END DO




        compwrt=0
END FUNCTION DG_VOL_INTEGRAL








FUNCTION DG_VOL_INTEGRAL2(N)
!> @brief
!> Calculates the volume integral term in the DG RHS for scalar linear advection with speed = 1
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:NOF_VARIABLES)::DG_VOL_INTEGRAL2
INTEGER::I,J,K,NQP,I_QP,I_VAR
REAL::PH,INTEG
REAL,DIMENSION(1:NOF_VARIABLES)::DG_SOL2
REAL,DIMENSION(1:idegfree)::BASIS_TEMP


    compwrt=-2


    NQP = ielem(n,iconsidered)%iTOTALPOINTS
    





    NUMBER=IELEM(N,ICONSIDERED)%IORDER
    NUMBER_OF_DOG = IELEM(N,ICONSIDERED)%IDEGFREE
         DG_VOL_INTEGRAL2(:) = 0.0d0  

    
    
     
       
        DO I_QP = 1, NQP 
            X1=QP_ARRAY(ICONSIDERED)%X(I_QP)
            Y1=QP_ARRAY(ICONSIDERED)%Y(I_QP)
            if (dimensiona.eq.3)then
            z1=QP_ARRAY(ICONSIDERED)%Z(I_QP)
            end if
            
            
            
            
            if (dimensiona.eq.2)then
            
            BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)
            else
            BASIS_TEMP = BASIS_REC(N, X1, Y1, z1,NUMBER, ICONSIDERED, NUMBER_OF_DOG)
            
            end if

             

            DO I_VAR = 1, NOF_VARIABLES
                DG_SOL2(I_VAR) = DOT_PRODUCT(BASIS_TEMP(1:NUMBER_OF_DOG), U_C(ICONSIDERED)%VALDG(1,I_VAR,2:NUMBER_OF_DOG+1))
            END DO 
            
            
            
              DG_VOL_INTEGRAL2 = DG_VOL_INTEGRAL2+ (QP_ARRAY(ICONSIDERED)%QP_WEIGHT(I_QP) *DG_SOL2/IELEM(N,ICONSIDERED)%TOTVOLUME)
             

         END DO


            DG_VOL_INTEGRAL2(:) =U_C(ICONSIDERED)%VALDG(1,:,1)+DG_VOL_INTEGRAL2(:)

    compwrt=0

END FUNCTION DG_VOL_INTEGRAL2


SUBROUTINE RECONSTRUCT_DG(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I_FACE, I_ELEM, I_QP,iqp
        
        
        
!$OMP DO
DO I_ELEM = 1, XMPIELRANK(N)
    compwrt=-2

    DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA
            !SOMEWHERE PRESTORED THE GUASSIAN QUADRATURE POINTS FOR YOUR SIDES IN ANOTHER SUBROUTINE AND YOU BUILD ONLY THE cleft STATES AND VOLUME INTEGRAL
            
            if (dimensiona.eq.2)then

            iqp=QP_LINE_N
            else
                if (ielem(n,I_ELEM)%types_faces(I_FACE).eq.5)then
					iqp=qp_quad
				  else
					iqp=QP_TRIANGLE
					
				  end if
            end if
            
        DO I_QP = 1,iqp! QP_LINE_N
            
            FACEX=I_FACE
            POINTX=I_QP
            ICONSIDERED=I_ELEM
            NUMBER_OF_DOG=IELEM(N,I_ELEM)%IDEGFREE

            
            ILOCAL_RECON3(ICONSIDERED)%ULEFT_DG(:, FACEX, POINTX) = DG_SOLFACE(N)
         
            
            
            
            
        END DO
    END DO
    
     compwrt=0           
    
END DO
!$OMP END DO

END SUBROUTINE RECONSTRUCT_DG


SUBROUTINE RECONSTRUCT_BR2_DG
IMPLICIT NONE
INTEGER::I_FACE, I_ELEM, I_QP, IQP,J,K
REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q, WEIGHTS_T

! Needed for DG_SURF_FLUX, used in BR2_LOCAL_LIFT
IF( DIMENSIONA == 3) THEN
    call QUADRATUREQUAD3D(N,IGQRULES)
!      WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)

    call QUADRATURETRIANG(N,IGQRULES)
!      WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)
ELSE
    CALL QUADRATURELINE(N,IGQRULES)
END IF

!$OMP DO
DO I_ELEM = 1, XMPIELRANK(N)
    ICONSIDERED = I_ELEM

    DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA

        ! Needed for BR2_LOCAL_LIFT
        if (dimensiona.eq.2)then
            iqp=QP_LINE_N
        else
            if (ielem(n,I_ELEM)%types_faces(I_FACE).eq.5)then
                iqp=qp_quad
!                  WEIGHTS_DG = WEIGHTS_Q
            else
                iqp=QP_TRIANGLE
!                  WEIGHTS_DG = WEIGHTS_T
            end if
        end if


        FACEX = I_FACE
        ILOCAL_RECON3(ICONSIDERED)%BR2_LOCAL_LIFT(:,:,I_FACE) = BR2_LOCAL_LIFT(N, IQP)

        DO I_QP = 1, IQP
            POINTX=I_QP
            IF (BR2_YN /= 1 .AND. ITESTCASE == 4) THEN
                X1 = ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX, POINTX,1)
                Y1 = ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX, POINTX,2)
                IF (DIMENSIONA == 3) Z1 = ILOCAL_RECON3(ICONSIDERED)%QPOINTS(FACEX, POINTX,3)

                LEFTV_DER = DG_SOL_DER()
                CALL DCONS2DPRIM
                ILOCAL_RECON3(ICONSIDERED)%BR2_AUX_VAR(:,:,I_FACE,I_QP) = LEFTV_DER + ILOCAL_RECON3(ICONSIDERED)%BR2_LOCAL_LIFT(:,:,I_FACE) * BR2_DAMPING

                !WRITE(150+N,*)"FIRST",ICONSIDERED,IELEM(N,ICONSIDERED)%ISHAPE,ILOCAL_RECON3(ICONSIDERED)%BR2_AUX_VAR(:,:,I_FACE,I_QP)
                !WRITE(150+N,*)"SECOND",LEFTV_DER
                !WRITE(150+N,*)"THIRD",ILOCAL_RECON3(ICONSIDERED)%BR2_LOCAL_LIFT(:,:,I_FACE)
                !WRITE(150+N,*)BR2_DAMPING

                !NOW COPY TO OTHER ARRAY
                ILOCAL_RECON3(ICONSIDERED)%ULEFTV(1:DIMS,1,I_FACE,I_QP)=ILOCAL_RECON3(ICONSIDERED)%BR2_AUX_VAR(NOF_VARIABLES,1:DIMS,I_FACE,I_QP)

    DO J=2,NOF_VAriables-1
        ILOCAL_RECON3(ICONSIDERED)%ULEFTV(1:DIMS,J,I_FACE,I_QP)=ILOCAL_RECON3(ICONSIDERED)%BR2_AUX_VAR(J,1:DIMS,I_FACE,I_QP)
    END DO
                !END COPYING
            END IF
        END DO
    END DO








END DO
!$OMP END DO

END SUBROUTINE RECONSTRUCT_BR2_DG

SUBROUTINE GET_LEFT_RIGHT_STATES(ICONSIDERED, FACEX, POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I

    I = ICONSIDERED

    IF (DIMENSIONA == 3) THEN
        IF (IELEM(N,I)%INTERIOR == 0) THEN
            CALL GET_STATES_INTERIOR(ICONSIDERED,FACEX,POINTX)
        ELSE
            CALL GET_STATES_BOUNDS(ICONSIDERED,FACEX,POINTX)
        END IF
    ELSE
        IF (IELEM(N,I)%INTERIOR == 0) THEN
            CALL GET_STATES_INTERIOR2D(ICONSIDERED,FACEX,POINTX)
        ELSE
            CALL GET_STATES_BOUNDS2D(ICONSIDERED,FACEX,POINTX)
        END IF
    END IF

END SUBROUTINE GET_LEFT_RIGHT_STATES


SUBROUTINE  GET_STATES_INTERIOR(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP
L=FACEX
NGP=POINTX
I=ICONSIDERED


                      IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                    ELSE !FV

				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state

				      END IF

					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					END IF

					IF (SRF.EQ.1)THEN

                        SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)
                        CALL ROTATEF(N,TRI,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
					END IF



END SUBROUTINE GET_STATES_INTERIOR



SUBROUTINE  GET_STATES_INTERIOR2D(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP
L=FACEX
NGP=POINTX
I=ICONSIDERED


IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables)= ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES,IELEM(N,I)%INEIGHN(L), NGP)
                        ELSE !FV

                        CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
                        CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state

                        END IF
!
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					END IF


END SUBROUTINE GET_STATES_INTERIOR2D


SUBROUTINE CALCULATE_INTERIOR_VISCOUS(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex
L=FACEX
NGP=POINTX
I=ICONSIDERED



					IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
						ELSE !FV

				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state

				      END IF



				      LCVGRAD(1,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,2,L,NGP);LCVGRAD(2,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,3,L,NGP);
				      LCVGRAD(3,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,4,L,NGP);LCVGRAD(4,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP);
				      RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,IELEM(N,I)%INEIGHN(L),NGP);

					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)


					  do nvar=1,turbulenceequations+passivescalar
					  RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
					  LCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(I)%ULEFTTURBV(1:3,nvar,L,NGP)
					  end do


					END IF



END SUBROUTINE CALCULATE_INTERIOR_VISCOUS


SUBROUTINE CALCULATE_INTERIOR_VISCOUS2D(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex
L=FACEX
NGP=POINTX
I=ICONSIDERED




					  IF (DG.EQ.1) THEN
                        CLEFT(1:nof_Variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT(1:nof_Variables) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
						LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
					  RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP)
						ELSE !FV
							CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
				      RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP)

						end if



					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF


					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)


					  do nvar=1,turbulenceequations+passivescalar
					  RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
					  LCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,nvar,L,NGP)
					  end do


					END IF





END SUBROUTINE CALCULATE_INTERIOR_VISCOUS2D


SUBROUTINE CALCULATE_BOUNDED_VISCOUS(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex,kk

L=FACEX
NGP=POINTX
I=ICONSIDERED




						IF (DG == 1) THEN
									CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
									else
								CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
								END IF



!   CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,2,L,NGP);LCVGRAD(2,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,3,L,NGP);
				      LCVGRAD(3,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,4,L,NGP);LCVGRAD(4,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP);

					  IF (SRF.EQ.1)THEN
						SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)
						CALL ROTATEF(N,TRI,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
				END IF

					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
						do nvar=1,turbulenceequations+passivescalar
						LCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(I)%ULEFTTURBV(1:3,nvar,L,NGP)
						end do
					end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if ((ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(L))%icode.eq.50))then

								  IF (DG == 1) THEN
                                    CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                    ELSE
                                    CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                    END IF



								  RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,L,NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,L,NGP);
								  RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,L,NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,L,NGP);
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF

									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do


								    END IF
									IF(PER_ROT.EQ.1)THEN
                                        CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                        DO KK=1,4
                                            RCVGRAD(KK,1:3)=ROTATE_PER_1(RCVGRAD(KK,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                        END DO
                                        RCVGRAD_T(1,1:3)=ROTATE_PER_1(RCVGRAD_T(1,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                    END IF



								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_innerx(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)

								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)


								    RCVGRAD(:,:)=LCVGRAD(:,:)
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				  				    RCVGRAD_T(:,:)=LCVGRAD_T(:,:)
				  				    end if
				  				     if (B_CODE.eq.4)then
				  				    rightv=zero
				  				    rightv(2:4)=LCVGRAD(4,1:3)
				  				    leftv=zero

				  				    CALL ROTATEF(N,TRI,leftv,rightv,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
								    leftv(2)=-leftv(2)
								    CALL ROTATEB(N,INVTRI,rightv,leftv,ANGLE1,ANGLE2)
				  				    RCVGRAD(4,1:3)=rightv(2:4)
				  				  	end if
!
!
								  END IF
							ELSE
							      IF (DG == 1) THEN
                                CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE
							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
							      END IF
							      RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,IELEM(N,I)%INEIGHN(L),NGP);
							      RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,IELEM(N,I)%INEIGHN(L),NGP);
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
								    END IF




							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if ((ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(L))%icode.eq.50))then	!periodic in other
									  IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								ELSE

									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                END IF

									  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE
                                                        IF (IEX.EQ.1)THEN
                                                RCVGRAD(4,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
                                                        ELSE
                                                RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
                                                        END IF

										      END IF


										END DO
									  END DO






								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF



									  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO

									  IF(PER_ROT.EQ.1)THEN
										CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
										  DO KK=1,4
											  RCVGRAD(KK,1:3)=ROTATE_PER_1(RCVGRAD(KK,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
										  END DO
										  RCVGRAD_T(1,1:3)=ROTATE_PER_1(RCVGRAD_T(1,1:3),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
										END IF



								END IF
							ELSE

								   IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                ELSE
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  END IF
								  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE



										      IF (IEX.EQ.1)THEN
									  RCVGRAD(4,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      END IF

										      END IF
										END DO
									  END DO



								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF

									 DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO


								    END IF

!
							END IF
					    END IF



END SUBROUTINE CALCULATE_BOUNDED_VISCOUS



SUBROUTINE CALCULATE_BOUNDED_VISCOUS2D(ICONSIDERED,FACEX,POINTX)
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,ITTT,nvar,iex

L=FACEX
NGP=POINTX
I=ICONSIDERED

IF (DG == 1) THEN
	CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
	else
CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
END IF


				      LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
						do nvar=1,turbulenceequations+passivescalar
						LCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,nvar,L,NGP)
						end do
					end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  !CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)

									IF (DG == 1) THEN
										CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
										ELSE
										CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
										END IF


								  RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,L,NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,L,NGP);
								  RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,L,NGP);
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF

									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do


								    END IF




								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_inner2dx(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)


								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)

								    RCVGRAD(:,:)=LCVGRAD(:,:)


								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				  				    RCVGRAD_T(:,:)=LCVGRAD_T(:,:)
				  				    end if
! 				  				     if (B_CODE.eq.4)then
!
!
! 				  				    rightv=zero
! 				  				    rightv(2:3)=LCVGRAD(3,1:2)
! 				  				    leftv=zero
!
! 				  				    CALL ROTATEF2d(N,TRI,leftv,rightv,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 								    leftv(2)=-leftv(2)
! 								    CALL ROTATEB2d(N,INVTRI,rightv,leftv,ANGLE1,ANGLE2)
! 				  				    RCVGRAD(3,1:2)=rightv(2:3)
! 				  				  	end if
!
!
								  END IF
							ELSE


								IF (DG == 1) THEN
									CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
									ELSE
									  CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
									  END IF




							      RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
							      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP);
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									do nvar=1,turbulenceequations+passivescalar
									RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
								    END IF




							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									 ! CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)

									IF (DG == 1) THEN
										CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									ELSE

										  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									END IF


									  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE



										      IF (IEX.EQ.1)THEN
									  RCVGRAD(3,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      END IF

										      END IF
										END DO
									  END DO






								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF




									  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO

								END IF



								END IF
							ELSE

								  !CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								IF (DG == 1) THEN
									CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									ELSE
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									  END IF


								ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (DG.EQ.1)THEN
 										      RCVGRAD(IEX,nVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
!
										      ELSE


										      IF (IEX.EQ.1)THEN
									  RCVGRAD(3,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      END IF

										      END IF
										END DO
									  END DO



								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF

									 DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										END DO
									  END DO


								    END IF

!
							END IF
					    END IF






END SUBROUTINE CALCULATE_BOUNDED_VISCOUS2D


SUBROUTINE  GET_STATES_BOUNDS(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP
L=FACEX
NGP=POINTX
I=ICONSIDERED



						IF (DG == 1) THEN
									CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
									else
								CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
								END IF

								IF (SRF.EQ.1)THEN
									SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)
									CALL ROTATEF(N,TRI,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
							END IF



								IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
										CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
									ELSE
										CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									end if
								end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  !if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
									if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then

								 IF (DG == 1) THEN
                                    CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                    ELSE
                                    CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                    END IF

									IF(PER_ROT.EQ.1)THEN
										CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
									  END IF

								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF




								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_innerx(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)

								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode



								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)


								  END IF
							ELSE


                                   IF (DG == 1) THEN
                                CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE
							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
							      END IF

								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF




							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								!if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then
								 IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								ELSE

									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                END IF
								IF(PER_ROT.EQ.1)THEN
									CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
								  END IF

								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF



								END IF
							ELSE
                                   IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                ELSE
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  END IF

!
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF

!
							END IF
					    END IF






END SUBROUTINE GET_STATES_BOUNDS



SUBROUTINE  GET_STATES_BOUNDS2D(ICONSIDERED,FACEX,POINTX)
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED, FACEX, POINTX
INTEGER::I,L,NGP,IKAS
L=FACEX
NGP=POINTX
I=ICONSIDERED

IF (DG == 1) THEN
                        CLEFT(1:nof_Variables)= ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                    ELSE

				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)

                    END IF

					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
					end if


					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                        IF (DG == 1) THEN
                                            CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                        ELSE !FV
                                            CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                        END IF
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
                                        if (icoupleturb.eq.1)then
                                        CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
                                        (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
                                        ELSE
                                        CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
                                        END IF
								    END IF





								  IKAS=1

								  ELSE
								  !NOT PERIODIC ONES IN MY CPU


								  CALL coordinates_face_inner2dx(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)

								    Poy(1)=cords(2)
								    Pox(1)=cords(1)


								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
!
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
!
				  				  	 IKAS=2

								  END IF
							ELSE
                                    IF (DG == 1) THEN
                                        CRIGHT(1:nof_Variables)= ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
        !
                                    ELSE !FV
                                        CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
                                    END IF

								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF


							       IKAS=3

							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS






							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                        IF (DG == 1) THEN
                                            CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                        else
                                            CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
 									   END IF
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF




								    END IF



								END IF
							ELSE
                                    IF (DG == 1) THEN
                                        CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                    ELSE
                                        CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                    END IF
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF


								    END IF
								   IKAS=4
!
							END IF
					    END IF







END SUBROUTINE GET_STATES_BOUNDS2D













FUNCTION CALC_DELTA_XYZ(NUM_NODES, NUM_DIMS, NODES_IN)
!> @brief
!> Calculates the "delta x/y/z" as in Luo 2012 eq 3.12
    IMPLICIT NONE
    INTEGER,INTENT(IN)::NUM_NODES, NUM_DIMS
    REAL,DIMENSION(:,:),INTENT(IN)::NODES_IN ! (NODE, DIMENSION)
    INTEGER::I_NODES, I_DIM
    REAL::XYZ_MAX,XYZ_MIN
    REAL,DIMENSION(NUM_DIMS)::CALC_DELTA_XYZ
    
    DO I_DIM = 1, NUM_DIMS
        XYZ_MAX = -1E12
        XYZ_MIN = 1E12
        
        DO I_NODES = 1, NUM_NODES
            IF (NODES_IN(I_NODES,I_DIM) > XYZ_MAX) XYZ_MAX = NODES_IN(I_NODES,I_DIM)
            IF (NODES_IN(I_NODES,I_DIM) < XYZ_MIN) XYZ_MIN = NODES_IN(I_NODES,I_DIM)
        END DO
        
        CALC_DELTA_XYZ(I_DIM) = 0.5 * ABS(XYZ_MAX - XYZ_MIN)
    END DO
    
END FUNCTION CALC_DELTA_XYZ


SUBROUTINE ALLOCATE_DG
!> @brief
!> ALLOCATES THE GAUSSIAN QUADRATURE VOLUME POINTS
    IMPLICIT NONE
    INTEGER::I, K, I_QP, N_QP, I_FACE,nnd,iqp,idummy,loopc
    real,dimension(1:idegfree+1)::tempint
    real::tempf,tempx
    
    
    ALLOCATE(QP_ARRAY(XMPIELRANK(N))); !Allocates for 2D
	

	

    DO I = 1, XMPIELRANK(N)
    
    SELECT CASE(ielem(n,i)%ishape)
        
        
        case(1) !hexa
        
        ALLOCATE(QP_ARRAY(I)%X(QP_TETRA*6))
        ALLOCATE(QP_ARRAY(I)%Y(QP_TETRA*6))
        ALLOCATE(QP_ARRAY(I)%Z(QP_TETRA*6))
        ALLOCATE(QP_ARRAY(I)%QP_WEIGHT(QP_TETRA*6))
        
        case(2) !TETRA
        ALLOCATE(QP_ARRAY(I)%X(QP_TETRA))
        ALLOCATE(QP_ARRAY(I)%Y(QP_TETRA))
        ALLOCATE(QP_ARRAY(I)%Z(QP_TETRA))
        ALLOCATE(QP_ARRAY(I)%QP_WEIGHT(QP_TETRA))
        
         case(3) !PYRAMID
        ALLOCATE(QP_ARRAY(I)%X(QP_TETRA*2))
        ALLOCATE(QP_ARRAY(I)%Y(QP_TETRA*2))
        ALLOCATE(QP_ARRAY(I)%Z(QP_TETRA*2))
        ALLOCATE(QP_ARRAY(I)%QP_WEIGHT(QP_TETRA*2))
        
        case(4) !PRISM
        ALLOCATE(QP_ARRAY(I)%X(QP_TETRA*3))
        ALLOCATE(QP_ARRAY(I)%Y(QP_TETRA*3))
        ALLOCATE(QP_ARRAY(I)%Z(QP_TETRA*3))
        ALLOCATE(QP_ARRAY(I)%QP_WEIGHT(QP_TETRA*3))
        
        
        case(5) !quad
        ALLOCATE(QP_ARRAY(I)%X(QP_TRIANGLE*2))
        ALLOCATE(QP_ARRAY(I)%Y(QP_TRIANGLE*2))
        ALLOCATE(QP_ARRAY(I)%QP_WEIGHT(QP_TRIANGLE*2))
        
        case(6) !tetra
        ALLOCATE(QP_ARRAY(I)%X(QP_TRIANGLE))
        ALLOCATE(QP_ARRAY(I)%Y(QP_TRIANGLE))
        ALLOCATE(QP_ARRAY(I)%QP_WEIGHT(QP_TRIANGLE))
        
        END SELECT
    

    ALLOCATE(ILOCAL_RECON3(I)%ULEFT_DG(NOF_VARIABLES, IELEM(N,I)%IFCA, NUMBEROFPOINTS2))
    
     IF((ITESTCASE == 4).or.((governingequations.EQ.-1).AND.(BR2_YN.EQ.1))) THEN !NS
            ALLOCATE(ILOCAL_RECON3(I)%BR2_AUX_VAR(NOF_VARIABLES, DIMENSIONA, IELEM(N,I)%IFCA, NUMBEROFPOINTS2))
            ALLOCATE(ILOCAL_RECON3(I)%BR2_LOCAL_LIFT(NOF_VARIABLES, DIMENSIONA, IELEM(N,I)%IFCA))
            ILOCAL_RECON3(I)%BR2_AUX_VAR = ZERO
            ILOCAL_RECON3(I)%BR2_LOCAL_LIFT = ZERO

        END IF
        
    
    
    END DO
    
    
END SUBROUTINE ALLOCATE_DG



SUBROUTINE PRESTORE_DG1
!> @brief
!> Prestores IELEM(N,I)%DELTA_XYZ, QP_ARRAY, SURF_QPOINTS, mass matrix
    IMPLICIT NONE
    INTEGER::I, K, I_QP, N_QP, I_FACE,nnd,iqp,idummy,loopc
    real,dimension(1:idegfree+1)::tempint
    real::tempf,tempx
        compwrt=-2
        I=iconsidered
        !Store volume quadrature points
        ELTYPE=IELEM(N,I)%ISHAPE
        ELEM_DEC=ielem(n,i)%vdec
        DO K = 1,IELEM(N,I)%NONODES

            
            NODES_LIST(k,1:dims)=INODER(IELEM(N,I)%NODES(K))%CORD(1:dims)
            
            VEXT(k,1:dims)=NODES_LIST(k,1:dims)
        END DO
        
        
        if (dimensiona.eq.2)then
        
        CALL DECOMPOSE2
        else
        CALL DECOMPOSE3
        
        end if
       
        
        SELECT CASE(ielem(n,i)%ishape)
        
        
        case(1,2,3,4) !hexa
        COUNT_1=0
            
             DO K=1,ELEM_DEC
                 VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
            
                CALL QUADRATURETETRA(N,IGQRULES)
                
                VOLTEMP=TETRAVOLUME(N)
                
                DO I_QP = 1, QP_TETRA
                    COUNT_1=COUNT_1+1
                     QP_ARRAY(I)%X(COUNT_1) = (QPOINTS(1,I_QP)- IELEM(N,I)%XXC)
                    QP_ARRAY(I)%Y(COUNT_1) = (QPOINTS(2,I_QP)- IELEM(N,I)%yyC)
                    QP_ARRAY(I)%Z(COUNT_1) = (QPOINTS(3,I_QP)- IELEM(N,I)%zzC)
                    
                    QP_ARRAY(I)%QP_WEIGHT(COUNT_1) = WEQUA3D(I_QP) * voltemp
                
                end do
            end do
            
            
            
        
       
        
        
        
        CASE(5)
            COUNT_1=0
            
             DO K=1,ELEM_DEC
                 VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            
                CALL QUADRATUREtriangle(N,IGQRULES)
                
                VOLTEMP=TRIANGLEVOLUME(N)
                
                DO I_QP = 1, QP_Triangle
                    COUNT_1=COUNT_1+1
                    
                    QP_ARRAY(I)%X(COUNT_1) = (QPOINTS(1,I_QP)- IELEM(N,I)%XXC)
                    QP_ARRAY(I)%y(COUNT_1) = (QPOINTS(2,I_QP)- IELEM(N,I)%yyC)
                   
                    
                    
                    QP_ARRAY(I)%QP_WEIGHT(COUNT_1) = WEQUA3D(I_QP) * voltemp
                    
                    
                END DO
                
            END DO
            

            
        CASE(6)
            CALL QUADRATURETRIANGLE(N,IGQRULES)
            N_QP = QP_Triangle
            VOLTEMP=TRIANGLEVOLUME(N)
            
            DO I_QP = 1, N_QP
            
               
                QP_ARRAY(I)%X(I_QP) = (QPOINTS(1,I_QP) - IELEM(N,I)%XXC)
                    QP_ARRAY(I)%Y(I_QP) = (QPOINTS(2,I_QP) - IELEM(N,I)%YYC)
                  
                    
                    
                    QP_ARRAY(I)%QP_WEIGHT(I_QP) = WEQUA3D(I_QP) * voltemp
               
            END DO
            
            
        END SELECT
        
        
        
        
        
        !Store volume quadrature points
        INTEG_BASIS_DG(ICONSIDERED)%value(1:IDEGFREE)=zero
        tempint=zero
        
            N_QP = ielem(n,iconsidered)%iTOTALPOINTS

	
        
            

                DO I_QP = 1, N_QP
                    x1=QP_ARRAY(I)%X(I_QP); y1=QP_ARRAY(I)%Y(I_QP)
                  
            
            if (dimensiona.eq.3)then
            z1=QP_ARRAY(I)%Z(I_QP)
           
            end if
!                    
                    IXX=ICONSIDERED
                    
                    if (dimensiona.eq.2)then
                    NUMBER=IELEM(N,ICONSIDERED)%IORDER
                    tempint(1:IDEGFREE)=tempint(1:IDEGFREE)+(BASIS_REC2D(N,X1,Y1,number,IXX,IDEGFREE)*QP_ARRAY(I)%QP_WEIGHT(I_QP))
                    
                    else
                    NUMBER=IELEM(N,ICONSIDERED)%IORDER
                    tempint(1:IDEGFREE)=tempint(1:IDEGFREE)+(BASIS_REC(N,X1,Y1,z1,number,IXX,IDEGFREE)*QP_ARRAY(I)%QP_WEIGHT(I_QP))
                    
                    end if
                    
                END DO
!                 
                
                INTEG_BASIS_DG(ICONSIDERED)%value(1:IDEGFREE)=tempint(1:IDEGFREE)
               
           
        
        compwrt=0
        
                

END SUBROUTINE





SUBROUTINE PRESTORE_DG2
!> @brief
!> Prestores IELEM(N,I)%DELTA_XYZ, QP_ARRAY, SURF_QPOINTS, mass matrix
    IMPLICIT NONE
    INTEGER::I, K, I_QP, N_QP, I_FACE,nnd,iqp,idummy,COUNT_1,loopc,ngp,l
    real,dimension(1:idegfree+1)::tempint
    real::tempf,tempx
    
    
    
    
    
END SUBROUTINE 






SUBROUTINE BUILD_MASS_MATRIX(N)
!> @brief
!> Assembles the mass matrix
!> REQUIRES: Globals: IELEM, QP_QUAD, QP_TRIANGLE, MASS_MATRIX
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM, I_QP, N_QP, I_DOF, J_DOF, KMAXE
    REAL::INTEG_TEST,integ_sm1,integ_sm2,phx1,phx2,KRON,maxs,mins,higher
    
    KMAXE = XMPIELRANK(N)
    
    
    
   
   
    !$OMP DO
    DO I_ELEM = 1, KMAXE
          totalmm(:,:) = ZERO
           INVmm(:,:) = ZERO
            compwrt=-2
           n_qp=ielem(n,i_elem)%iTOTALPOINTS
       
      
        
        
        iconsidered=I_ELEM
    NUMBER_OF_DOG = IELEM(N,I_ELEM)%IDEGFREE
    
    
    
        
    DO I_DOF = 1, NUM_DG_DOFS       
        DO J_DOF = 1, NUM_DG_DOFS
            INTEG_MM = ZERO
            
            
                DO I_QP = 1, N_QP
                    Ixx = I_ELEM; 
                    
                    
                    x1=QP_ARRAY(I_ELEM)%X(I_QP); 
                    y1=QP_ARRAY(I_ELEM)%Y(I_QP)
                    if (DIMENSIONA.eq.3)then
                    Z1=QP_ARRAY(I_ELEM)%Z(I_QP);
                    end if
                    
                    
            
             kron=1.0d0
             
             
             if (i_dof.ne.j_dof)then
             kron=1.0d0
             end if
                    
                    if (dimensiona.eq.2)then
                    NUMBER=IELEM(N,ICONSIDERED)%IORDER
                    BASIS_VECTOR(1:idegfree) = BASIS_REC2D(N,X1,Y1,number,IXX,NUMBER_OF_DOG)
                    else
                    NUMBER=IELEM(N,ICONSIDERED)%IORDER
                    BASIS_VECTOR(1:idegfree) = BASIS_REC(N,X1,Y1,z1,number,IXX,NUMBER_OF_DOG)
                    
                    
                    end if
                    
                        IF (I_DOF==1.and.J_DOF==1)THEN
                            PHX = 1.0d0 
                            
                            
                        ELSE IF (I_DOF==1.and.J_DOF/=1)THEN
                            PHX = BASIS_VECTOR(J_DOF-1)
                            
                            
                        ELSE IF (I_DOF/=1.and.J_DOF==1)THEN
                            PHX = BASIS_VECTOR(I_DOF-1) 
                            
                            
                            
                        ELSE IF (I_DOF/=1.and.J_DOF/=1)THEN
                            PHX = BASIS_VECTOR(I_DOF-1)*BASIS_VECTOR(J_DOF-1)
                            
                        END IF
                        

                    
                    INTEG_MM = INTEG_MM + PHX*QP_ARRAY(I_ELEM)%QP_WEIGHT(I_QP)*kron
                    
                    
                     

                END DO
                
                
                
        
        totalmm(I_DOF, J_DOF) = totalmm(I_DOF, J_DOF) + INTEG_MM
       
!         if (ielem(n,iconsidered)%ishape.eq.2)then
       
        
!         end if
        END DO
    END DO
    
    
    
        CALL COMPMASSINV
!
        m_1(i_elem)%val(:,:)=invmm(:,:)
!
!         maxs=tolsmall
!         mins=tolbig
!         DO I_DOF = 1, NUM_DG_DOFS
!         DO J_DOF = 1, NUM_DG_DOFS
!         if (abs(m_1(i_elem)%val(i_dof,j_dof)).gt.maxs)then
!         maxs=abs(m_1(i_elem)%val(i_dof,j_dof))
!
!         end if
!
!         if (abs(m_1(i_elem)%val(i_dof,j_dof)).lt.mins)then
!         mins=abs(m_1(i_elem)%val(i_dof,j_dof))
!
!         end if
!
!         end do
!         end do
!
!
!
!         if ((maxs/mins).gt.higher)then
!         higher=maxs/mins
!         end if
        
        
        
        
        
!         if ((ielem(n,iconsidered)%ishape.eq.2).AND.(ielem(n,iconsidered)%iORDER.Ge.2).AND.( IELEM(N,iconsidered)%LUMP.GT.IEVERYAV))then
!         write(510+n,*)maxs/mins
! !         if ((m_1(i_elem)%val(1,1).ne.m_1(i_elem)%val(1,1)))then
! !         if ((maxs/mins.gt.10e20).or.(m_1(i_elem)%val(1,1).ne.m_1(i_elem)%val(1,1)))then
!         
!         
!         
!         
!         
!            
!            totalmm(:,:) = ZERO
!             compwrt=-2
!            n_qp=ielem(n,i_elem)%iTOTALPOINTS
!        
!        
!         
!         
!         
!     NUMBER_OF_DOG = IELEM(N,I_ELEM)%IDEGFREE
!     
!     
!     
!         
!     
!             DO I_DOF = 1, NUM_DG_DOFS       
!         DO J_DOF = 1, NUM_DG_DOFS
!             INTEG_MM=zero
!                 DO I_QP = 1, N_QP
!                     Ixx = I_ELEM; 
!                     
!                     
!                     x1=QP_ARRAY(I_ELEM)%X(I_QP); 
!                     y1=QP_ARRAY(I_ELEM)%Y(I_QP)
!                     if (DIMENSIONA.eq.3)then
!                     Z1=QP_ARRAY(I_ELEM)%Z(I_QP);
!                     end if
!                     
!                     
!            
!             
!              kron=1.0d0
!              
!              
!              if (i_dof.ne.j_dof)then
!              kron=0.0d0
!              end if
!                     
!                     if (dimensiona.eq.2)then
!                     NUMBER=IELEM(N,ICONSIDERED)%IORDER
!                     BASIS_VECTOR(1:idegfree) = BASIS_REC2D(N,X1,Y1,number,IXX,NUMBER_OF_DOG)
!                     else
!                     NUMBER=IELEM(N,ICONSIDERED)%IORDER
!                     BASIS_VECTOR(1:idegfree) = BASIS_REC(N,X1,Y1,z1,number,IXX,NUMBER_OF_DOG)
!                     
!                     
!                     end if
!                     
!                         IF (I_DOF==1.and.J_DOF==1)THEN
!                             PHX = 1.0d0 
!                             
!                             
!                         ELSE IF (I_DOF==1.and.J_DOF/=1)THEN
!                             PHX = BASIS_VECTOR(J_DOF-1)
!                             
!                             
!                         ELSE IF (I_DOF/=1.and.J_DOF==1)THEN
!                             PHX = BASIS_VECTOR(I_DOF-1) 
!                             
!                             
!                             
!                         ELSE IF (I_DOF/=1.and.J_DOF/=1)THEN
!                             PHX = BASIS_VECTOR(I_DOF-1)*BASIS_VECTOR(J_DOF-1)
!                             
!                         END IF
!                         
! 
!                     
!                     INTEG_MM = INTEG_MM + PHX*QP_ARRAY(I_ELEM)%QP_WEIGHT(I_QP)*kron
!                     
!                     
!                      
! 
!                 END DO
!                 
!                 
!                 
!         
!         totalmm(I_DOF, J_DOF) = totalmm(I_DOF, J_DOF) + INTEG_MM
!        
!         
!         
!         END DO
!     END DO
!     
!     
!     
!         CALL COMPMASSINV
!         
!         m_1(i_elem)%val(:,:)=invmm(:,:)
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
!         
! !         end if
!         end if
        
        
         compwrt=0
         
         
         
         
         
        
    
END DO
!$OMP END DO

    
    
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
END SUBROUTINE BUILD_MASS_MATRIX

SUBROUTINE COMPMASSINV
!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
IMPLICIT NONE
 
integer :: i,j,k,l,m,irow,num_dofs
real:: big,dum
real,DIMENSION(NUM_DG_DOFS,NUM_DG_DOFS)::a,b
num_dofs=NUM_DG_DOFS

a(:,:)=totalMM(:,:)
b(:,:)=zero




do i = 1,num_dofs
    do j = 1,num_dofs
        b(i,j) = 0.0d0
    end do
    b(i,i) = 1.0d0
end do 

do i = 1,num_dofs 
   big = a(i,i)
   do j = i,num_dofs
     if (a(j,i).gt.big) then
       big = a(j,i)
       irow = j
     end if
   end do
   ! interchange lines i with irow for both a() and b() matrices
   if (big.gt.a(i,i)) then
     do k = 1,num_dofs
       dum = a(i,k)                      ! matrix a()
       a(i,k) = a(irow,k)
       a(irow,k) = dum
       dum = b(i,k)                 ! matrix b()
       b(i,k) = b(irow,k)
       b(irow,k) = dum
     end do
   end if
   ! divide all entries in line i from a(i,j) by the value a(i,i); 
   ! same operation for the identity matrix
   dum = a(i,i)
   do j = 1,num_dofs
     a(i,j) = a(i,j)/dum
     b(i,j) = b(i,j)/dum
   end do
   ! make zero all entries in the column a(j,i); same operation for indent()
   do j = i+1,num_dofs
     dum = a(j,i)
     do k = 1,num_dofs
       a(j,k) = a(j,k) - dum*a(i,k)
       b(j,k) = b(j,k) - dum*b(i,k)               
            
     end do
   end do
end do
  
 do i = 1,num_dofs-1
   do j = i+1,num_dofs
     dum = a(i,j)
     do l = 1,num_dofs
       a(i,l) = a(i,l)-dum*a(j,l)
       b(i,l) = b(i,l)-dum*b(j,l)
     end do
   end do
 end do

 

 
 
 invMM(:,:)=b(:,:)
 
 
 END SUBROUTINE COMPMASSINV
  
!  function inv(totalMM) result(Ainv)
!     implicit none
!     real,intent(in) :: totalMM(:,:)
!     real            :: Ainv(size(totalMM,1),size(totalMM,2))
!     real            :: work(size(totalMM,1))            ! work array for LAPACK
!     integer         :: vv,info,ipiv(size(totalMM,1))     ! pivot indices
! 
!     ! Store A in Ainv to prevent it from being overwritten by LAPACK
!     Ainv = totalMM
!     vv = size(totalMM,1)
!     ! SGETRF computes an LU factorization of a general M-by-N matrix A
!     ! using partial pivoting with row interchanges.
!     call SGETRF(vv,vv,Ainv,vv,ipiv,info)
!     if (info.ne.0) stop 'Matrix is numerically singular!'
!     ! SGETRI computes the inverse of a matrix using the LU factorization
!     ! computed by SGETRF.
!     call SGETRI(vv,Ainv,vv,ipiv,work,vv,info)
!     if (info.ne.0) stop 'Matrix inversion failed!'
! end function inv
  
  
  
  
  
 
  
  
  
  
  
  
  
  
  
  








END MODULE DG_FUNCTIONS
