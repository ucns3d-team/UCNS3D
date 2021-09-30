MODULE DG_FUNCTIONS

USE BASIS
USE DECLARATION
USE DERIVATIVES

IMPLICIT NONE

 CONTAINS

FUNCTION DG_SOL(n)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
    REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
    integer,intent(in)::n
    INTEGER::I_DOF, I_VAR
    REAL,DIMENSION(Nof_VARIABLES)::DG_SOL

!     IF(ALL(SHAPE(U_C_VALDG) /= (/ NUM_VARIABLES, NUM_DOFS+1 /))) THEN
!         WRITE(400+N,*) 'DG_SOL: U_C_VALDG WRONG DIMENSIONS:', SHAPE(U_C_VALDG)
!         STOP
!     END IF
    number=ielem(n,iconsidered)%iorder
    BASIS_TEMP = BASIS_REC2D(N, X1, Y1, number, ICONSIDERED, NUMBER_OF_DOG)

    DO I_VAR = 1, NOF_VARIABLES
        DG_SOL(I_VAR) = U_C(iconsidered)%VALDG(1,I_VAR,1) + DOT_PRODUCT(BASIS_TEMP, U_C(iconsidered)%VALDG(1,I_VAR,2:))
    END DO

END FUNCTION DG_SOL





FUNCTION DG_SOLFACE(n)
IMPLICIT NONE
!> @brief
!> This function returns the DG solution at a given point (X_IN, Y_IN)\n
!> REQUIRES: X_IN, Y_IN: coordinates of the point where the solution is requested, NUM_VARIABLES: number of solution variables, NUM_DOFS: number of basis terms
    REAL,DIMENSION(NUMBER_OF_DOG)::BASIS_TEMP
    INTEGER::I_DOF, I_VAR
    integer,intent(in)::n
    REAL,DIMENSION(1:nof_Variables)::DG_SOLFACE

!     IF(ALL(SHAPE(U_C_VALDG) /= (/ nof_Variables, NUMBER_OF_DOG+1 /))) THEN
!         WRITE(400+N,*) 'DG_SOL: U_C_VALDG WRONG DIMENSIONS:', SHAPE(U_C_VALDG)
!         STOP
!     END IF
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
REAL,DIMENSION(NUMBER_OF_DOG+1,NOf_VARIABLES)::DG_solfacex
integer::i
integer,intent(in)::n
NUMBER=IELEM(N,ICONSIDERED)%IORDER
X1=ilocal_recon3(iconsidered)%surf_qpoints(facex,pointx,1)
Y1=ilocal_recon3(iconsidered)%surf_qpoints(facex,pointx,2)
DO I = 1, NOf_VARIABLES
            DG_SOLFACEX(1,I) = HLLCFLUX(I) * WEQUA2D(pointx)*IELEM(N,ICONSIDERED)%SURF(FACEX)
            DG_SOLFACEX(2:,I) = HLLCFLUX(I) * WEQUA2D(pointx)*IELEM(N,ICONSIDERED)%SURF(FACEX)*BASIS_REC2D(N,X1,Y1,NUMBER,ICONSIDERED,NUMBER_OF_DOG)
        END DO

END FUNCTION DG_SOLFACEX

FUNCTION DG_RHS_INTEGRAL_vol(N)
!> @brief
!> Calculates the volume or surface integral term in the DG RHS for scalar linear advection with speed = 1
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I_VAR,iqp
    REAL,DIMENSION(idegfree+1,NOF_VARIABLES)::DG_RHS_INTEGRAL_VOL
    real,dimension(nof_variables)::flux_term
    INTEGER::I,J,K,ngp
    REAL::PH,INTEG
    
!     IF (SIZE(FLUX_TERM) /= NUM_VARS) THEN
!         WRITE(400+N,*) 'DG_RHS_INTEGRAL: FLUX_TERM WRONG DIMENSIONS:', SHAPE(FLUX_TERM)
!         STOP
!     END IF
    
    
            IF (IELEM(N,Iconsidered)%ISHAPE == 5) IQP = QP_QUAD
            IF (IELEM(N,Iconsidered)%ISHAPE == 6) IQP = QP_TRIANGLE
            number=ielem(n,iconsidered)%iorder
            number_of_dog=ielem(n,iconsidered)%idegfree
           
            DO I_VAR = 1, NOF_VARIABLES
            DG_RHS_INTEGRAL_VOL(1,I_VAR) = 0.0d0
            DG_RHS_INTEGRAL_VOL(2:,I_VAR)=0.0D0

             do I=1,NUMBER_OF_DOG
             DO NGP = 1, IQP 
             X1=qp_array(iconsidered,ngp)%x
             Y1=qp_array(iconsidered,ngp)%Y
             
             FLUX_TERM=DG_SOL(n)
            
             if (itestcase.eq.3)then
             leftv=DG_SOL(n)
             call prim2cons2d(n)
             flux_term=FLUXEVAL2D(LEFTV)
             
             end if
             
             
             DG_RHS_INTEGRAL_VOL(i+1,I_VAR) = DG_RHS_INTEGRAL_vol(i+1,I_VAR)+FLUX_TERM(I_VAR)* qp_array(iconsidered,ngp)%QP_WEIGHT * ielem(n,iconsidered)%totvolume *(DF2DX(X1,Y1,I)+DF2DY(X1,Y1,I))
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
        DO I_FACE = 1, IELEM(N,I_ELEM)%IFCA
            !SOMEWHERE PRESTORED THE GUASSIAN QUADRATURE POINTS FOR YOUR SIDES IN ANOTHER SUBROUTINE AND YOU BUILD ONLY THE cleft STATES AND VOLUME INTEGRAL
            
            DO I_QP = 1, QP_LINE_N
            
            FACEX=I_FACE
            pointx=I_QP
            ICONSIDERED=I_ELEM
            number_of_dog=ielem(n,iconsidered)%idegfree
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
    
	ALLOCATE(QP_ARRAY(XMPIELRANK(N),NUMBEROFPOINTS)); !Allocates for 2D
    
    DO I = 1, XMPIELRANK(N)    
        !Store volume quadrature points
        DO K = 1,IELEM(N,I)%NONODES
            NODES_LIST(k,1:2)=INODER(IELEM(N,I)%NODES(K))%CORD(1:2)
            VEXT(k,1:2)=NODES_LIST(k,1:2)
        END DO
        
        !Store delta xyz (normalization factor from Luo 2012)
        IELEM(N,I)%DELTA_XYZ = CALC_DELTA_XYZ(IELEM(N,I)%NONODES, DIMENSIONA, NODES_LIST)
    
        SELECT CASE(ielem(n,i)%ishape)
        CASE(5)
            CALL QUADRATUREQUAD(N,IGQRULES)
            N_QP = QP_quad
        CASE(6)
            CALL QUADRATURETRIANGLE(N,IGQRULES)
            N_QP = QP_Triangle
        END SELECT
                
        DO I_QP = 1, N_QP
            QP_ARRAY(I,I_QP)%X = QPOINTS(1,I_QP) - IELEM(N,I)%XXC
            QP_ARRAY(I,I_QP)%Y = QPOINTS(2,I_QP) - IELEM(N,I)%YYC
            QP_ARRAY(I,I_QP)%QP_WEIGHT = WEQUA3D(I_QP)
        END DO
        
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
    
    ALLOCATE(MASS_MATRIX_CENTERS(N:N,KMAXE,NUM_DG_DOFS,NUM_DG_DOFS)); MASS_MATRIX_CENTERS(N,:,:,:) = ZERO; MASS_MATRIX_CENTERS(N,:,1,1) = 1.0D0
    
    DO I_ELEM = 1, KMAXE
        SELECT CASE(IELEM(N,I_ELEM)%ISHAPE)
        CASE(5) ! Quadrilateral
            N_QP = QP_QUAD
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
                INTEG_TEST = INTEG_TEST + (PHX*QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT*IELEM(N,I_ELEM)%TOTVOLUME)
            END DO
            
            MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) = MASS_MATRIX_CENTERS(N, I_ELEM, 1, I_DOF+1) + INTEG_TEST
            MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) = MASS_MATRIX_CENTERS(N, I_ELEM, I_DOF+1, 1) + INTEG_TEST
            
            DO J_DOF = 1, IDEGFREE
                INTEG_TEST = ZERO
                DO I_QP = 1, N_QP
                    IXX = I_ELEM; X1 = QP_ARRAY(I_ELEM,I_QP)%X; Y1 = QP_ARRAY(I_ELEM,I_QP)%Y
                    BASIS_VECTOR = BASIS_REC2D(N,X1,Y1,IORDER,IXX,IDEGFREE)
                
                    PHX = BASIS_VECTOR(I_DOF) * BASIS_VECTOR(J_DOF)
                    INTEG_TEST = INTEG_TEST + (PHX*QP_ARRAY(I_ELEM,I_QP)%QP_WEIGHT*IELEM(N,I_ELEM)%TOTVOLUME)
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
        WRITE(300+N,*) 'Inverse,', INV_MASS_MATRIX(N,I_ELEM,:,:)
        WRITE(300+N,*) 'Identity', MATMUL(MASS_MATRIX_CENTERS(N,I_ELEM,:,:),INV_MASS_MATRIX(N,I_ELEM,:,:))
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

FUNCTION INVERSE_MATRIX(MATRIX_IN, N_COLS)
    INTEGER,INTENT(IN)::N_COLS
    REAL,DIMENSION(N_COLS,N_COLS),INTENT(IN)::MATRIX_IN
    REAL,DIMENSION(1,N_COLS,N_COLS)::DUMMY_MATRIX_IN, DUMMY_INVERSE_MATRIX
    REAL,DIMENSION(N_COLS,N_COLS)::INVERSE_MATRIX
    
    DUMMY_MATRIX_IN(1,:,:) = MATRIX_IN
    CALL COMPMASSINV(DUMMY_MATRIX_IN, DUMMY_INVERSE_MATRIX, N_COLS)
    INVERSE_MATRIX = DUMMY_INVERSE_MATRIX(1,:,:)
    
END FUNCTION INVERSE_MATRIX

FUNCTION NORMAL_LSQ(A, b)
    REAL,DIMENSION(:,:),INTENT(IN)::A
    REAL,DIMENSION(:),INTENT(IN)::b
    REAL,DIMENSION(SIZE(A,1))::NORMAL_LSQ
    
    NORMAL_LSQ = MATMUL(MATMUL(INVERSE_MATRIX(MATMUL(TRANSPOSE(A), A), SIZE(A,2)), TRANSPOSE(A)), b)
    
END FUNCTION NORMAL_LSQ

! FUNCTION LUO_LSQ_RECONSTRUCT(N, N_DIM, ORDER, DEGFREE)
! !> @brief
! !> This function reconstructs an approximation of NUM_DG_RECONSTRUCT_DOFS
! !> REQUIRES: IELEM, U_C as globals
!     INTEGER,INTENT(IN)::N, N_DIM, ORDER, DEGFREE
!     INTEGER::I_ELEM, I_FACE, NEIGHBOR_INDEX, I_DIM, KMAXE
!     REAL,DIMENSION(DEGFREE)::BASIS_TEMP
!     REAL,DIMENSION(IELEM(N, I_ELEM)%IFCA * 3, NUM_DG_RECONSTRUCT_DOFS)::LHS_MATRIX ! See eq. 3.19 of Luo 2012
!     REAL,DIMENSION(IELEM(N, I_ELEM)%IFCA * 3)::RHS_DG_RECONSTRUCT ! See eq. 3.19 of Luo 2012
!     REAL,DIMENSION(XMPIELRANK(N),NUM_DG_RECONSTRUCT_DOFS)::LUO_LSQ_RECONSTRUCT ! NUM_DG_RECONSTRUCT_DOFS
!     
!     KMAXE = XMPIELRANK(N)
!     
!     DO I_ELEM = 1, KMAXE
!         LHS_MATRIX = ZERO
!         IF (IELEM(N, I_ELEM)%INTERIOR == 0)THEN ! Element is interior
!             DO I_FACE = 1, IELEM(N, I_ELEM)%IFCA
!                 NEIGHBOR_INDEX = IELEM(N,I)%INEIGH(I_FACE)
!                 OFFSET = (N_DIM + 1) * (I_FACE - 1)
!             
!                 BASIS_TEMP = BASIS_REC2D(N,IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, IORDER,I_ELEM,IDEGFREE)
!                 
!                 DO I_DG_RECONSTRUCT_DOF = 1, NUM_DG_RECONSTRUCT_DOFS
!                     LHS_MATRIX(OFFSET+1,I_DG_RECONSTRUCT_DOF) = BASIS_TEMP(N_DIM+I_DG_RECONSTRUCT_DOF-1)
!                 END DO
!                 
!                 LHS_MATRIX(OFFSET+2,1) = BASIS_TEMP(1)
!                 LHS_MATRIX(OFFSET+2,N_DIM+1) = BASIS_TEMP(2)
!                 LHS_MATRIX(OFFSET+3,2) = BASIS_TEMP(2)
!                 LHS_MATRIX(OFFSET+3,N_DIM+1) = BASIS_TEMP(1)
!                 
!                 RHS_DG_RECONSTRUCT(OFFSET+1) = DG_SOL(N, I_ELEM, IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, NOF_VARIABLES, ORDER, DEGFREE, U_C(NEIGHBOR_INDEX)%VALDG(1,:,:)) - DG_SOL(N, I_ELEM, IELEM(N, NEIGHBOR_INDEX)%XXC, IELEM(N,NEIGHBOR_INDEX)%YYC, NOF_VARIABLES, ORDER, DEGFREE, U_C(I_ELEM)%VALDG(1,:,:))
!                 DO I_DIM = 1, N_DIM
!                     RHS_DG_RECONSTRUCT(OFFSET+I_DIM+1) = IELEM(N,I_ELEM)%DELTA_XYZ(I_DIM) / IELEM(N,NEIGHBOR_INDEX)%DELTA_XYZ(I_DIM) * U_C(I_ELEM)%VALDG(1,:,I_DIM+1) - U_C(NEIGHBOR_INDEX)%VALDG(1,:,I_DIM+1)
!                 END DO
!                 
!             END DO
!         END IF
!         
!         LUO_LSQ_RECONSTRUCT(I_ELEM,:) = NORMAL_LSQ(LHS_MATRIX, RHS_DG_RECONSTRUCT)
!         
!     END DO
!     
! END FUNCTION LUO_LSQ_RECONSTRUCT

END MODULE DG_FUNCTIONS
