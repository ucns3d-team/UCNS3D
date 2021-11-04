MODULE DG_FUNCTIONS

USE BASIS
USE DECLARATION
USE DERIVATIVES

IMPLICIT NONE

 CONTAINS

FUNCTION DG_SOL(N)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
    REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
    INTEGER,INTENT(IN)::N
    INTEGER::I_VAR
    REAL,DIMENSION(Nof_VARIABLES)::DG_SOL

    NUMBER = IELEM(N,ICONSIDERED)%IORDER
    BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)

    DO I_VAR = 1, NOF_VARIABLES
        DG_SOL(I_VAR) = U_C(ICONSIDERED)%VALDG(1,I_VAR,1) + DOT_PRODUCT(BASIS_TEMP, U_C(ICONSIDERED)%VALDG(1,I_VAR,2:))
    END DO

END FUNCTION DG_SOL

FUNCTION DG_SOL_I_CENTER(N, I)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
    REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
    INTEGER,INTENT(IN)::N, I
    INTEGER::I_VAR
    REAL,DIMENSION(NOF_VARIABLES)::DG_SOL_I_CENTER

    NUMBER = IELEM(N, I)%IORDER
    X1 = IELEM(N, I)%XXC
    Y1 = IELEM(N, I)%YYC
    BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, I, NUMBER_OF_DOG)

    DO I_VAR = 1, NOF_VARIABLES
        DG_SOL_I_CENTER(I_VAR) = U_C(I)%VALDG(1, I_VAR, 1) + DOT_PRODUCT(BASIS_TEMP, U_C(I)%VALDG(1, I_VAR, 2:))
    END DO

END FUNCTION DG_SOL_I_CENTER

FUNCTION DG_SOLFACE(n)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
    REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
    INTEGER::I_DOF, I_VAR
    integer,intent(in)::n
    REAL,DIMENSION(1:nof_Variables)::DG_SOLFACE

    x1=ilocal_recon3(iconsidered)%surf_qpoints(facex,pointx,1)
    y1=ilocal_recon3(iconsidered)%surf_qpoints(facex,pointx,2)
    number=ielem(n,iconsidered)%iorder
        
    BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER, ICONSIDERED, NUMBER_OF_DOG)

    DO I_VAR = 1, NOF_VARIABLES
        DG_SOLface(I_VAR) = U_C(iconsidered)%VALDG(1,I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(1:number_of_dog), U_C(iconsidered)%VALDG(1,I_VAR,2:number_of_dog+1))
    END DO

END FUNCTION DG_SOLFACE

FUNCTION DG_RHS_INTEGRAL(N, I_ELEM, QP_X, QP_Y, QP_WEIGHT, NUM_VARS, ORDER, NUM_DOFS, CELL_VOL_OR_SURF, FLUX_TERM, VOL_OR_SURF)
!> @brief
!> Calculates the volume or surface integral term in the DG RHS for scalar linear advection with speed = 1
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, I_ELEM, ORDER, NUM_VARS, NUM_DOFS, VOL_OR_SURF
    REAL,INTENT(IN)::QP_X, QP_Y, QP_WEIGHT, CELL_VOL_OR_SURF
    REAL,DIMENSION(:),INTENT(IN)::FLUX_TERM
    INTEGER::I_VAR
    REAL,DIMENSION(NUM_DOFS+1,NUM_VARS)::DG_RHS_INTEGRAL
    
    IF (SIZE(FLUX_TERM) /= NUM_VARS) THEN
        WRITE(400+N,*) 'DG_RHS_INTEGRAL: FLUX_TERM WRONG DIMENSIONS:', SHAPE(FLUX_TERM)
        STOP
    END IF
    
    IF (VOL_OR_SURF == 1) THEN ! VOLUME INTEGRAL
        DO I_VAR = 1, NUM_VARS
            DG_RHS_INTEGRAL(1,I_VAR) = 0
            DG_RHS_INTEGRAL(2:,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * (BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS,1) + BASIS_REC2D_DERIVATIVE(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS,2)) ! For linear advection with speed = 1
        END DO
    ELSE IF (VOL_OR_SURF == 2) THEN ! SURFACE INTEGRAL
        DO I_VAR = 1, NUM_VARS
            DG_RHS_INTEGRAL(1,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF
            DG_RHS_INTEGRAL(2:,I_VAR) = FLUX_TERM(I_VAR) * QP_WEIGHT * CELL_VOL_OR_SURF * BASIS_REC2D(N,QP_X,QP_Y,ORDER,I_ELEM,NUM_DOFS)
        END DO
    END IF

END FUNCTION DG_RHS_INTEGRAL

FUNCTION DG_SOLFACEX(n)
    IMPLICIT NONE
    REAL,DIMENSION(NUMBER_OF_DOG+1,NOF_VARIABLES)::DG_SOLFACEX
    INTEGER::I
    INTEGER,INTENT(IN)::N
    
    NUMBER=IELEM(N,ICONSIDERED)%IORDER
    X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
    Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
    DO I = 1, NOF_VARIABLES
        DG_SOLFACEX(1,I) = HLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)
        DG_SOLFACEX(2:,I) = HLLCFLUX(I) * WEQUA2D(POINTX)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC2D(N,X1,Y1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
    END DO

END FUNCTION DG_SOLFACEX

FUNCTION DG_RHS_INTEGRAL_vol(N)
!> @brief
!> Calculates the volume integral term in the DG RHS for scalar linear advection with speed = 1
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_VAR,iqp
    REAL,DIMENSION(idegfree+1,NOF_VARIABLES)::DG_RHS_INTEGRAL_VOL
    real,dimension(nof_variables)::FLUX_TERM
    INTEGER::I_DOF,ngp
    
    IF (IELEM(N,Iconsidered)%ISHAPE == 5) IQP = QP_TRIANGLE * 2
    IF (IELEM(N,Iconsidered)%ISHAPE == 6) IQP = QP_TRIANGLE
    number=ielem(n,iconsidered)%iorder
    number_of_dog=ielem(n,iconsidered)%idegfree
    
    DG_RHS_INTEGRAL_VOL = 0.0d0
    DO I_VAR = 1, NOF_VARIABLES
        DO I_DOF=1,NUMBER_OF_DOG
            DO NGP = 1, IQP
                X1=qp_array(iconsidered,ngp)%x
                Y1=qp_array(iconsidered,ngp)%Y
                
                FLUX_TERM=DG_SOL(n) !Linear advection
                
                if (itestcase.eq.3)then !Euler
                    leftv=DG_SOL(n)
                    call prim2cons2d(n)
                    FLUX_TERM=FLUXEVAL2D(LEFTV)
                end if
                
                DG_RHS_INTEGRAL_VOL(I_DOF+1,I_VAR) = DG_RHS_INTEGRAL_vol(I_DOF+1,I_VAR)+FLUX_TERM(I_VAR)* qp_array(iconsidered,ngp)%QP_WEIGHT * (DF2DX(X1,Y1,I_DOF)+DF2DY(X1,Y1,I_DOF))! * ielem(n,iconsidered)%totvolume Volume multiplication moved to qp_weight because decomposition
                
!                 WRITE(500+N,*) I_DOF, FLUX_TERM, DF2DX(X1,Y1,I_DOF), DF2DY(X1,Y1,I_DOF), DG_RHS_INTEGRAL_vol!, X1, Y1, qp_array(iconsidered,ngp)%QP_WEIGHT
            END DO
        END DO
    END DO
END FUNCTION DG_RHS_INTEGRAL_vol

SUBROUTINE RECONSTRUCT_DG(N)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_FACE, I_ELEM, I_QP
    
    !$OMP DO
    DO I_ELEM = 1, XMPIELRANK(N)
        ICONSIDERED = I_ELEM
        NUMBER_OF_DOG = IELEM(N,ICONSIDERED)%IDEGFREE
        
        DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA
            !SOMEWHERE PRESTORED THE GUASSIAN QUADRATURE POINTS FOR YOUR SIDES IN ANOTHER SUBROUTINE AND YOU BUILD ONLY THE cleft STATES AND VOLUME INTEGRAL
            FACEX = I_FACE
            
            DO I_QP = 1, QP_LINE_N
                POINTX = I_QP
                ILOCAL_RECON3(I_ELEM)%ULEFT_DG(:, I_FACE, I_QP) = DG_SOLFACE(n)
                
                write(700+n,*)"look here","cell",i_elem,"face",facex,"point",pointx
                write(700+n,*)ILOCAL_RECON3(I_ELEM)%ULEFT_DG(:, I_FACE, I_QP)
                
                !STORE IT HERE (ILOCAL_RECON3(I)%ULEFT_DG(1,1:FACES,1:NGP) ! you need to allocate it in memory
            END DO
        END DO
    END DO
    !$OMP END DO

END SUBROUTINE RECONSTRUCT_DG

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

SUBROUTINE PRESTORE_AND_ALLOCATE_DG
!> @brief
!> Prestores IELEM(N,I)%DELTA_XYZ, QP_ARRAY, SURF_QPOINTS, mass matrix
    IMPLICIT NONE
    INTEGER::I, K, I_QP, N_QP, I_FACE
    
	ALLOCATE(QP_ARRAY(XMPIELRANK(N),QP_Triangle*2)); !Allocates for 2D
    
    DO I = 1, XMPIELRANK(N)    
!         WRITE(600+N,*) I, IELEM(N,I)%XXC, IELEM(N,I)%YYC, IELEM(N,I)%TOTVOLUME
        !Store volume quadrature points
        DO K = 1,IELEM(N,I)%NONODES
            NODES_LIST(k,1:2)=INODER(IELEM(N,I)%NODES(K))%CORD(1:2)
            VEXT(k,1:2)=NODES_LIST(k,1:2)
            
!             WRITE(600+N,*) VEXT(k,1:2)
        END DO
        
        !Store delta xyz (normalization factor from Luo 2012)
        IELEM(N,I)%DELTA_XYZ = CALC_DELTA_XYZ(IELEM(N,I)%NONODES, DIMENSIONA, NODES_LIST)
    
        CALL DECOMPOSE2
    
        ALLOCATE(IELEM(N,I)%TAYLOR_INTEGRAL(IDEGFREE - DIMENSIONA))
        IELEM(N,I)%TAYLOR_INTEGRAL = ZERO
    
        SELECT CASE(ielem(n,i)%ishape)
        CASE(5)
            DO K=1,ELEM_DEC
                VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
            
                CALL QUADRATUREtriangle(N,IGQRULES)
                
                VOLTEMP=TRIANGLEVOLUME(N)
                
                DO I_QP = 1, QP_Triangle
                    QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
                    QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
                    QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%QP_WEIGHT = WEQUA3D(I_QP) * VOLTEMP
                    
                    IELEM(N,I)%TAYLOR_INTEGRAL(1) = IELEM(N,I)%TAYLOR_INTEGRAL(1) + (QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%X / IELEM(N,I)%DELTA_XYZ(1)) ** 2 / 2 * WEQUA3D(I_QP) ! * CELL_VOLUME No multiplication because cancelled out
                    IELEM(N,I)%TAYLOR_INTEGRAL(2) = IELEM(N,I)%TAYLOR_INTEGRAL(2) + (QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%Y / IELEM(N,I)%DELTA_XYZ(2)) ** 2 / 2 * WEQUA3D(I_QP) ! * CELL_VOLUME No multiplication because cancelled out
                    IELEM(N,I)%TAYLOR_INTEGRAL(3) = IELEM(N,I)%TAYLOR_INTEGRAL(3) + QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%X / IELEM(N,I)%DELTA_XYZ(1) *  QP_ARRAY(I,I_QP + QP_TRIANGLE * (K - 1))%Y / IELEM(N,I)%DELTA_XYZ(2) * WEQUA3D(I_QP) ! * CELL_VOLUME No multiplication because cancelled out
                END DO
                
            END DO
        CASE(6)
            CALL QUADRATURETRIANGLE(N,IGQRULES)
            N_QP = QP_Triangle
            
            DO I_QP = 1, N_QP
                QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
                QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
                QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP) * IELEM(N,I)%totvolume
            END DO
            
        END SELECT
        
        ALLOCATE(ILOCAL_RECON3(I)%SURF_QPOINTS(IELEM(N,I)%IFCA, QP_LINE_N, DIMENSIONA))
        !Store surface quadrature points
        DO I_FACE = 1, IELEM(N,I)%IFCA
            VEXT(1,1:2) = inoder(IELEM(N,I)%NODES_FACES(I_FACE,1))%CORD(1:2)  !COPY THE COORDINATE OF THE FIRST NODE OF THID EDGE
            VEXT(2,1:2) = inoder(IELEM(N,I)%NODES_FACES(I_FACE,2))%CORD(1:2)  !COPY THE COORDINATE OF THE SECOND NODE OF THID EDGE
            CALL QUADRATURELINE(N,IGQRULES)

            DO I_QP = 1, QP_LINE_N
                ILOCAL_RECON3(I)%SURF_QPOINTS(I_FACE,I_QP,:) = QPOINTS2D(:,I_QP) - (/ IELEM(N,I)%XXC, IELEM(N,I)%YYC /)
            END DO
        END DO
        
        ALLOCATE(ILOCAL_RECON3(I)%ULEFT_DG(NOF_VARIABLES, IELEM(N,I)%IFCA, QP_LINE_N))
    END DO
    
    CALL ASS_MASS_MATRIX(N)
    
END SUBROUTINE PRESTORE_AND_ALLOCATE_DG

SUBROUTINE ASS_MASS_MATRIX(N)
!> @brief
!> Assembles the mass matrix
!> REQUIRES: Globals: IELEM, QP_QUAD, QP_TRIANGLE, MASS_MATRIX
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM, I_QP, N_QP, I_DOF, J_DOF, KMAXE
    REAL,DIMENSION(IDEGFREE)::BASIS_VECTOR
    REAL::INTEG_TEST,PHX
    
    KMAXE = XMPIELRANK(N)
    
    ALLOCATE(MASS_MATRIX_CENTERS(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); MASS_MATRIX_CENTERS(N,:,:,:) = ZERO; 
    
    DO I_ELEM = 1, KMAXE
        SELECT CASE(IELEM(N,I_ELEM)%ISHAPE)
        CASE(5) ! Quadrilateral
            N_QP = QP_TRIANGLE * 2
        CASE(6) ! Triangle
            N_QP = QP_TRIANGLE
        END SELECT
        
        MASS_MATRIX_CENTERS(N,I_ELEM,1,1) = IELEM(N,I_ELEM)%TOTVOLUME

        !TAKIS START
        DO I_DOF = 1, IDEGFREE
            INTEG_TEST = ZERO
            DO I_QP = 1, N_QP
                IXX = I_ELEM; X1 = QP_ARRAY(I_ELEM,I_QP)%X; Y1 = QP_ARRAY(I_ELEM,I_QP)%Y
                BASIS_VECTOR = BASIS_REC2D(N,X1,Y1,IORDER,IXX,IDEGFREE)
            
                PHX = BASIS_VECTOR(I_DOF)
                INTEG_TEST = INTEG_TEST + PHX*QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT !* IELEM(N,I_ELEM)%TOTVOLUME
!                 WRITE(300+N,*) X1, Y1, "BASIS_VECTOR", BASIS_VECTOR, QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
            END DO
            
            MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) + INTEG_TEST
            MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) + INTEG_TEST
            
            DO J_DOF = 1, IDEGFREE
                INTEG_TEST = ZERO
                DO I_QP = 1, N_QP
                    IXX = I_ELEM; X1 = QP_ARRAY(I_ELEM,I_QP)%X; Y1 = QP_ARRAY(I_ELEM,I_QP)%Y
                    BASIS_VECTOR = BASIS_REC2D(N,X1,Y1,IORDER,IXX,IDEGFREE)
                
                    PHX = BASIS_VECTOR(I_DOF) * BASIS_VECTOR(J_DOF)
                    INTEG_TEST = INTEG_TEST + PHX*QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT !* IELEM(N,I_ELEM)%TOTVOLUME
                END DO
                MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, J_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, J_DOF+1) + INTEG_TEST
            END DO
        END DO
        !TAKIS END
        
!         DO I_QP = 1, N_QP
!             IXX=I_ELEM;X1=I_ELEM,I_QP)%X;Y1=QP_ARRAY(I_ELEM,I_QP)%Y
!             BASIS_VECTOR = BASIS_REC2D(N,X1,Y1,IORDER,IXX,IDEGFREE)
!             
!             DO I_DOF = 1, IDEGFREE
!                 MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) + BASIS_VECTOR(I_DOF) * QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!                 MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) + BASIS_VECTOR(I_DOF) * QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!                 
!                 DO J_DOF = 1, IDEGFREE
!                     MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, J_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF, J_DOF) + BASIS_VECTOR(I_DOF) * BASIS_VECTOR(I_DOF) * QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!                 END DO
!             END DO
!         END DO
!         
!         MASS_MATRIX_CENTERS(N, I_ELEM, :, :) = MASS_MATRIX_CENTERS(N, I_ELEM, :, :) * IELEM(N,I_ELEM)%TOTVOLUME
!     
    END DO
    
    ALLOCATE(INV_MASS_MATRIX(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); INV_MASS_MATRIX(N,:,:,:) = ZERO
    
    CALL COMPMASSINV(MASS_MATRIX_CENTERS(N,:,:,:), INV_MASS_MATRIX(N,:,:,:), NUM_DG_DOFS)
    
    DO I_ELEM = 1, KMAXE
        WRITE(300+N,*) I_ELEM
!         DO I_QP = 1, N_QP
!             WRITE(300+N,*) 'QP_ARRAY', QP_ARRAY(I_ELEM,I_QP)%X, QP_ARRAY(I_ELEM,I_QP)%Y, QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT
!             WRITE(300+N,*) 'BASIS,', BASIS_REC2D(N,QP_ARRAY(I_ELEM,I_QP)%X,QP_ARRAY(I_ELEM,I_QP)%Y,IORDER,I_ELEM,IDEGFREE)
!         END DO
        WRITE(300+N,*) 'XYZ', IELEM(N,I_ELEM)%NODES
!         WRITE(300+N,*) 'DELTAXYZ', IELEM(N,I_ELEM)%DELTA_XYZ
        WRITE(300+N,*) 'MMC', MASS_MATRIX_CENTERS(N,I_ELEM,:,:)
!         WRITE(300+N,*) 'Inverse,', INV_MASS_MATRIX(N,I_ELEM,:,:)
!         WRITE(300+N,*) 'Identity', MATMUL(MASS_MATRIX_CENTERS(N,I_ELEM,:,:),INV_MASS_MATRIX(N,I_ELEM,:,:))
    END DO
    
END SUBROUTINE ASS_MASS_MATRIX

SUBROUTINE COMPMASSINV(totalMM,invMM,N_DOFS)
!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
IMPLICIT NONE
 
integer :: i,j,k,l,m,irow,P,kmaxe
real:: big,dum
real,DIMENSION(N_DOFS,N_DOFS)::a,b
integer,INTENT(IN)::N_DOFS
REAL,DIMENSION(:,:,:),INTENT(IN)::totalMM
REAL,DIMENSION(:,:,:),INTENT(INOUT)::invMM
kmaxe=xmpielrank(n)
DO P=1,kmaxe

a(:,:)=totalMM(P,:,:)
b(:,:)=zero

do i = 1,N_DOFS
    do j = 1,N_DOFS
        b(i,j) = 0.0
    end do
    b(i,i) = 1.0
end do 

do i = 1,N_DOFS   
   big = a(i,i)
   do j = i,N_DOFS
     if (a(j,i).gt.big) then
       big = a(j,i)
       irow = j
     end if
   end do
   ! interchange lines i with irow for both a() and b() matrices
   if (big.gt.a(i,i)) then
     do k = 1,N_DOFS
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
   do j = 1,N_DOFS
     a(i,j) = a(i,j)/dum
     b(i,j) = b(i,j)/dum
   end do
   ! make zero all entries in the column a(j,i); same operation for indent()
   do j = i+1,N_DOFS
     dum = a(j,i)
     do k = 1,N_DOFS
       a(j,k) = a(j,k) - dum*a(i,k)
       b(j,k) = b(j,k) - dum*b(i,k)               
            
     end do
   end do
end do
  
 do i = 1,N_DOFS-1
   do j = i+1,N_DOFS
     dum = a(i,j)
     do l = 1,N_DOFS
       a(i,l) = a(i,l)-dum*a(j,l)
       b(i,l) = b(i,l)-dum*b(j,l)
     end do
   end do
 end do
 
 invMM(P,:,:)=b(:,:)
  
END DO
 
END SUBROUTINE COMPMASSINV

FUNCTION INVERSE_MATRIX(MATRIX_IN, N_DOFS)
!Calculate the inverse of the input matrix with Gauss-Jordan Elimination
IMPLICIT NONE
    INTEGER,INTENT(IN)::N_DOFS
    REAL,DIMENSION(N_DOFS,N_DOFS),INTENT(IN)::MATRIX_IN
    integer :: i,j,k,l,m,irow
    real:: big,dum
    real,DIMENSION(N_DOFS,N_DOFS)::a,b
    REAL,DIMENSION(N_DOFS,N_DOFS)::INVERSE_MATRIX

    a(:,:)=MATRIX_IN
    b(:,:)=zero

    do i = 1,N_DOFS
        do j = 1,N_DOFS
            b(i,j) = 0.0
        end do
        b(i,i) = 1.0
    end do 

    do i = 1,N_DOFS   
    big = a(i,i)
    do j = i,N_DOFS
        if (a(j,i).gt.big) then
        big = a(j,i)
        irow = j
        end if
    end do
    ! interchange lines i with irow for both a() and b() matrices
    if (big.gt.a(i,i)) then
        do k = 1,N_DOFS
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
    do j = 1,N_DOFS
        a(i,j) = a(i,j)/dum
        b(i,j) = b(i,j)/dum
    end do
    ! make zero all entries in the column a(j,i); same operation for indent()
    do j = i+1,N_DOFS
        dum = a(j,i)
        do k = 1,N_DOFS
        a(j,k) = a(j,k) - dum*a(i,k)
        b(j,k) = b(j,k) - dum*b(i,k)               
                
        end do
    end do
    end do
    
    do i = 1,N_DOFS-1
    do j = i+1,N_DOFS
        dum = a(i,j)
        do l = 1,N_DOFS
        a(i,l) = a(i,l)-dum*a(j,l)
        b(i,l) = b(i,l)-dum*b(j,l)
        end do
    end do
    end do
    
    INVERSE_MATRIX = b(:,:)
    
END FUNCTION INVERSE_MATRIX

FUNCTION NORMAL_LSQ(A, b)
    REAL,DIMENSION(:,:),INTENT(IN)::A
    REAL,DIMENSION(:),INTENT(IN)::b
    REAL,DIMENSION(SIZE(A,1))::NORMAL_LSQ
    
    NORMAL_LSQ = MATMUL(MATMUL(INVERSE_MATRIX(MATMUL(TRANSPOSE(A), A), SIZE(A,2)), TRANSPOSE(A)), b)
    
END FUNCTION NORMAL_LSQ

SUBROUTINE LUO_LSQ_RECONSTRUCT(N)
!> @brief
!> This function reconstructs an approximation of NUM_DG_RECONSTRUCT_DOFS
!> REQUIRES: IELEM, U_C as globals
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_ELEM, I_FACE, I_QP, I_VAR
    REAL,DIMENSION(NUM_DG_RECONSTRUCT_DOFS-1)::BASIS_TEMP
    REAL,DIMENSION(NOF_VARIABLES, NUM_DG_RECONSTRUCT_DOFS)::RECONSTRUCTED_SOL
    
    !$OMP DO
    DO I_ELEM = 1, XMPIELRANK(N)
        ICONSIDERED = I_ELEM
        NUMBER = IELEM(N,ICONSIDERED)%IORDER
        
        RECONSTRUCTED_SOL(:, 1:NUM_DG_DOFS) = U_C(ICONSIDERED)%VALDG(1, :, :)
        RECONSTRUCTED_SOL(:, NUM_DG_DOFS+1:NUM_DG_RECONSTRUCT_DOFS) = LUO_LSQ_RECONSTRUCT_SOL(N)
        
        DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA ! Set ULEFT_DG
            FACEX = I_FACE
            DO I_QP = 1, QP_LINE_N
                POINTX = I_QP
                
                X1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,1)
                Y1=ILOCAL_RECON3(ICONSIDERED)%SURF_QPOINTS(FACEX,POINTX,2)
                BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER + 1, ICONSIDERED, NUM_DG_RECONSTRUCT_DOFS-1)
                
                DO I_VAR = 1, NOF_VARIABLES
                    ILOCAL_RECON3(I_ELEM)%ULEFT_DG(I_VAR, I_FACE, I_QP) = RECONSTRUCTED_SOL(I_VAR,1) + DOT_PRODUCT(BASIS_TEMP(1:NUM_DG_RECONSTRUCT_DOFS-1), RECONSTRUCTED_SOL(I_VAR, 2:NUM_DG_RECONSTRUCT_DOFS))
                END DO
                
!                 write(700+n,*)"look here","cell",i_elem,"face",facex,"point",pointx
!                 write(700+n,*)ILOCAL_RECON3(I_ELEM)%ULEFT_DG(:, I_FACE, I_QP)
            END DO
        END DO
    END DO
    !$OMP END DO
    
END SUBROUTINE LUO_LSQ_RECONSTRUCT

FUNCTION LUO_LSQ_RECONSTRUCT_SOL(N)
!> @brief
!> Returns the higher dofs reconstructed solution using the Luo LSQ method
    IMPLICIT NONE
    INTEGER, INTENT(IN)::N
    INTEGER::NEIGHBOR_INDEX, I_DIM, I_FACE, I_DG_RECONSTRUCT_DOF, I_VAR, NEIGHBOR_TYPE
    REAL,DIMENSION(NUM_DG_RECONSTRUCT_DOFS - 1)::BASIS_TEMP
    REAL,DIMENSION(NOF_VARIABLES)::DG_SOL_TEMP, DG_SOL_I_J
    REAL,DIMENSION(DIMENSIONA)::NEIGHBOR_DELTA_XYZ
    REAL,DIMENSION(NOF_VARIABLES, NUM_DG_DOFS)::NEIGHBOR_U
    REAL,DIMENSION(NOF_VARIABLES, IELEM(N, ICONSIDERED)%IFCA * NUM_DG_DOFS, NUM_DG_RECONSTRUCT_DOFS - NUM_DG_DOFS)::LHS_MATRIX ! See eq. 3.19 of Luo 2012
    REAL,DIMENSION(NOF_VARIABLES, IELEM(N, ICONSIDERED)%IFCA * NUM_DG_DOFS)::RHS_DG_RECONSTRUCT ! See eq. 3.19 of Luo 2012
    REAL,DIMENSION(NOF_VARIABLES, NUM_DG_RECONSTRUCT_DOFS - NUM_DG_DOFS)::LUO_LSQ_RECONSTRUCT_SOL
    
    LHS_MATRIX = ZERO
    
    DO I_FACE = 1, IELEM(N, ICONSIDERED)%IFCA ! Construct LHS, RHS matrices
        NEIGHBOR_INDEX = IELEM(N,ICONSIDERED)%INEIGH(I_FACE)
        OFFSET = NUM_DG_DOFS * (I_FACE - 1) + 1 ! Incrementing LHS and RHS matrix indices for each face
    
        IF (IELEM(N, ICONSIDERED)%INTERIOR == 0)THEN ! Element is interior
            NEIGHBOR_TYPE = 0
        ELSE ! Element is not interior
            IF (IELEM(N,ICONSIDERED)%INEIGHB(I_FACE) == N) THEN ! Same CPU
                IF (IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE).GT.0) THEN ! Boundary
                    IF (IBOUND(N,IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE))%ICODE.EQ.5) THEN ! PERIODIC IN MY CPU
                        NEIGHBOR_TYPE = 0
                    ELSE
                        NEIGHBOR_TYPE = 1
                    END IF
                ELSE ! Not boundary
                    NEIGHBOR_TYPE = 0
                END IF
            ELSE ! Other CPU
                IF (IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE).GT.0) THEN ! Boundary
                    IF (IBOUND(N,IELEM(N,ICONSIDERED)%IBOUNDS(I_FACE))%ICODE.EQ.5) THEN ! PERIODIC IN OTHER CPU
                        NEIGHBOR_TYPE = 2
                    END IF
                ELSE ! Not boundary
                    NEIGHBOR_TYPE = 2
                END IF
            END IF
        END IF
        
        SELECT CASE(NEIGHBOR_TYPE)
        CASE(0) ! Interior or same CPU (periodic or not boundary)
            X1 = IELEM(N, NEIGHBOR_INDEX)%XXC
            Y1 = IELEM(N, NEIGHBOR_INDEX)%YYC
            BASIS_TEMP = BASIS_REC2D(N, X1, Y1, NUMBER + 1, NEIGHBOR_INDEX, NUM_DG_RECONSTRUCT_DOFS - 1)
            NEIGHBOR_U = U_C(NEIGHBOR_INDEX)%VALDG(1,:,:)
            NEIGHBOR_DELTA_XYZ = IELEM(N, NEIGHBOR_INDEX)%DELTA_XYZ(I_DIM)
            DG_SOL_TEMP = DG_SOL_I_CENTER(N, NEIGHBOR_INDEX) ! To call dg_sol without changing ICONSIDERED
        CASE(1) ! Same CPU boundary and not periodic
        CASE(2) ! Other CPU
            ! BASIS_TEMP (basis at center of neighboring element)
            BASIS_TEMP = IEXSOLHIR(ICONSIDERED)%BASIS_NEIGHBOR_CENTER(I_FACE, :)
            ! DELTA_XYZ
            NEIGHBOR_DELTA_XYZ = IEXSOLHIR(ICONSIDERED)%DELTA_XYZ(I_FACE, :)
            ! NEIGHBOR_U
            NEIGHBOR_U = IEXSOLHIR(ICONSIDERED)%SOL_DG(I_FACE,1:NOF_VARIABLES,1:NUM_DG_DOFS)
            DO I_VAR = 1, NOF_VARIABLES
                DG_SOL_TEMP(I_VAR) = NEIGHBOR_U(I_VAR,1) + DOT_PRODUCT(BASIS_TEMP, NEIGHBOR_U(I_VAR, 2:))
            END DO
        END SELECT
        
        DO I_DG_RECONSTRUCT_DOF = 1, NUM_DG_RECONSTRUCT_DOFS - NUM_DG_DOFS
            LHS_MATRIX(:, OFFSET, I_DG_RECONSTRUCT_DOF) = BASIS_TEMP(NUMBER_OF_DOG + I_DG_RECONSTRUCT_DOF) ! Higher order terms
        END DO
        
        ! 2D, first order terms
        LHS_MATRIX(:,OFFSET+1,1) = BASIS_TEMP(1) ! x, x term
        LHS_MATRIX(:,OFFSET+1,DIMENSIONA+1) = BASIS_TEMP(2) ! x, xy term
        LHS_MATRIX(:,OFFSET+2,2) = BASIS_TEMP(2) ! y, y term
        LHS_MATRIX(:,OFFSET+2,DIMENSIONA+1) = BASIS_TEMP(1) ! y, xy term
        
        DO I_VAR = 1, NOF_VARIABLES
            DG_SOL_I_J(I_VAR) = U_C(ICONSIDERED)%VALDG(1, I_VAR, 1) + DOT_PRODUCT(BASIS_TEMP, U_C(ICONSIDERED)%VALDG(1, I_VAR, 2:)) ! DG solution using values at i and basis at center of j
        END DO
        
        DG_SOL_TEMP = DG_SOL_TEMP - DG_SOL_I_J
        DO I_VAR = 1, NOF_VARIABLES
            RHS_DG_RECONSTRUCT(I_VAR,OFFSET) =  DG_SOL_TEMP(I_VAR) ! 1st term
            DO I_DIM = 1, DIMENSIONA ! 1st order terms
                RHS_DG_RECONSTRUCT(I_VAR, OFFSET+I_DIM) = IELEM(N, ICONSIDERED)%DELTA_XYZ(I_DIM) / NEIGHBOR_DELTA_XYZ * NEIGHBOR_U(I_VAR,I_DIM+1) - U_C(ICONSIDERED)%VALDG(1,I_VAR,I_DIM+1)
            END DO
        END DO
    END DO
    
    DO I_VAR = 1, NOF_VARIABLES
        LUO_LSQ_RECONSTRUCT_SOL(I_VAR, :) = NORMAL_LSQ(LHS_MATRIX(I_VAR,:,:), RHS_DG_RECONSTRUCT(I_VAR,:)) ! Solve for reconstructed solution
    END DO
    
    WRITE(700+N,*) LHS_MATRIX, 'RHS_DG_RECONSTRUCT', RHS_DG_RECONSTRUCT, 'LUO_LSQ_RECONSTRUCT_SOL', LUO_LSQ_RECONSTRUCT_SOL
    
END FUNCTION LUO_LSQ_RECONSTRUCT_SOL

END MODULE DG_FUNCTIONS
