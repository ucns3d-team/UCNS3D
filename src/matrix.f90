 MODULE LAPCK
 USE MPIINFO
 USE DECLARATION
  IMPLICIT NONE

  CONTAINS
 
  REAL FUNCTION LXNORM(XQR,PDIM)
  IMPLICIT NONE
   INTEGER, INTENT(IN) :: PDIM
   REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XQR
   INTEGER I
 
    LXNORM = 0.0d0
    DO I=1,PDIM
     LXNORM = LXNORM + XQR(I)**2
    ENDDO
    LXNORM = SQRT(LXNORM)
  END FUNCTION
  
  
  ! CONSTRUCT A HOUSEHOLDER VECTOR V THAT 
  ! ANNIHILATES ALL BUT THE FIRST COMPONENT OF X
  SUBROUTINE HOUSE(XQR,VQR1,PDIM)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: PDIM
    REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XQR
    REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VQR1

    VQR1 = XQR
    VQR1(1) = XQR(1) + SIGN(1.0,XQR(1))*LXNORM(XQR,PDIM)
    
  END SUBROUTINE


  ! CONSTRUCT A HOUSEHOLDER REFLECTION MATRIX
  ! FROM A HOUSEHOLDER VECTOR V 
  SUBROUTINE COMPUTEHOUSEMATRIX(PQR,VQR,IDEG)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDEG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::PQR
    REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VQR
    REAL ::VNORM
    INTEGER:: I,J

    PQR = 0.0d0
    DO I=1,IDEG
      PQR(I,I) = 1.0D0
    ENDDO
   
    VNORM = LXNORM(VQR,ideg)
    VQR = VQR/VNORM
    
    DO I=1,IDEG
    DO J=1,IDEG
      PQR(I,J) = PQR(I,J) - 2.0d0*VQR(I)*VQR(J)
    ENDDO
    ENDDO
    
  END SUBROUTINE

 !%%%%%%%%%%%%%%%%%%
  SUBROUTINE TRANSPOSEMATRIX(QFF,QTFF,IDEG)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDEG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN) :: QFF
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) :: QTFF
    INTEGER I,J
    DO I=1,IDEG
    DO J=1,IDEG
      QTFF (I,J) = QFF(J,I)
    ENDDO
    ENDDO
  END SUBROUTINE
  
  
  
  SUBROUTINE TRANSPOSEMATRIX_DG(QFF_DG,QTFF_DG,IDEG)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDEG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN) :: QFF_DG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) :: QTFF_DG
    INTEGER I,J
    DO I=1,IDEG
    DO J=1,IDEG
      QTFF_DG (I,J) = QFF_DG(J,I)
    ENDDO
    ENDDO
  END SUBROUTINE


  ! QR DECOMPOSITION
  SUBROUTINE  QRDECOMPOSITION(LSCQM,QFF,RFF,IDEG)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDEG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN) :: LSCQM
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::  QFF,RFF
    REAL::IDENTITY(IDEG,IDEG)
    REAL:: TEST(IDEG)
    INTEGER:: I,J,L
    REAL,ALLOCATABLE,DIMENSION(:)::XQR
    
    IDENTITY(1:IDEG,1:IDEG) = 0.0d0
    QFF(1:IDEG,1:IDEG)=0.0D0
    RFF(1:IDEG,1:IDEG)=0.0D0
    DO I=1,IDEG
      IDENTITY(I,I) = 1.0D0
    ENDDO

    
    QFF(1:IDEG,1:IDEG) = IDENTITY(1:IDEG,1:IDEG)
    RFF(1:IDEG,1:IDEG) = LSCQM(1:IDEG,1:IDEG)
    
    ALLOCATE(PQR(1:IDEG,1:IDEG),VQR(1:IDEG));PQR=0.0D0;VQR=0.0D0

    DO L=1,IDEG
     ! ALLOCATE VECTOR AND REFLECTION MATRIX
     PDIM = IDEG-L+1
     ALLOCATE(VQR1(1:PDIM),XQR(1:PDIM)) 
     XQR = RFF(L:IDEG,L)
     ! COMPUTE THE PARTIAL VECTOR    
     CALL HOUSE(XQR,VQR1,PDIM)
     VQR(1:IDEG) = 0.0d0
     VQR(L:IDEG) = VQR1
    
     
     
     ! COMPUTE THE REFLECTION MATRIX
     CALL COMPUTEHOUSEMATRIX(PQR,VQR,IDEG)
     ! CONSTRUCT THE Q(L) MATRIX
     
     RFF(1:IDEG,1:IDEG) = MATMUL(PQR(1:IDEG,1:IDEG),RFF(1:IDEG,1:IDEG))
     QFF(1:IDEG,1:IDEG) = MATMUL(QFF(1:IDEG,1:IDEG),PQR(1:IDEG,1:IDEG))
    
     DEALLOCATE(VQR1,XQR)
    ENDDO

    DEALLOCATE(PQR,VQR)



 
  END SUBROUTINE
  
  
  
  SUBROUTINE  QRDECOMPOSITION_DG(LSCQM_DG,QFF_DG,RFF_DG,IDEG)
	IMPLICIT NONE
    INTEGER, INTENT(IN) :: IDEG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN) :: LSCQM_DG
    REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT) ::  QFF_DG,RFF_DG
    REAL::IDENTITY(IDEG,IDEG)
    REAL:: TEST(IDEG)
    INTEGER:: I,J,L
    REAL,ALLOCATABLE,DIMENSION(:)::XQR
    
    IDENTITY(1:IDEG,1:IDEG) = 0.0d0
    QFF_dG(1:IDEG,1:IDEG)=0.0D0
    RFF_DG(1:IDEG,1:IDEG)=0.0D0
    DO I=1,IDEG
      IDENTITY(I,I) = 1.0D0
    ENDDO

    
    QFF_DG(1:IDEG,1:IDEG) = IDENTITY(1:IDEG,1:IDEG)
    RFF_DG(1:IDEG,1:IDEG) = LSCQM_DG(1:IDEG,1:IDEG)
    
    ALLOCATE(PQR(1:IDEG,1:IDEG),VQR(1:IDEG));PQR=0.0D0;VQR=0.0D0

    DO L=1,IDEG
     ! ALLOCATE VECTOR AND REFLECTION MATRIX
     PDIM = IDEG-L+1
     ALLOCATE(VQR1(1:PDIM),XQR(1:PDIM)) 
     XQR = RFF_DG(L:IDEG,L)
     ! COMPUTE THE PARTIAL VECTOR    
     CALL HOUSE(XQR,VQR1,PDIM)
     VQR(1:IDEG) = 0.0d0
     VQR(L:IDEG) = VQR1
    
     
     
     ! COMPUTE THE REFLECTION MATRIX
     CALL COMPUTEHOUSEMATRIX(PQR,VQR,IDEG)
     ! CONSTRUCT THE Q(L) MATRIX
     
     RFF_DG(1:IDEG,1:IDEG) = MATMUL(PQR(1:IDEG,1:IDEG),RFF_DG(1:IDEG,1:IDEG))
     QFF_DG(1:IDEG,1:IDEG) = MATMUL(QFF_DG(1:IDEG,1:IDEG),PQR(1:IDEG,1:IDEG))
    
     DEALLOCATE(VQR1,XQR)
    ENDDO

    DEALLOCATE(PQR,VQR)

    END SUBROUTINE


        DOUBLE PRECISION FUNCTION ddot_d(N,DX,INCX,DY,INCY)
!
!   -- Reference BLAS level1 routine --
!   -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!   -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
 !$omp declare target
!      .. Scalar Arguments ..
      INTEGER incx,incy,n
!      ..
!      .. Array Arguments ..
      DOUBLE PRECISION dx(*),dy(*)
!      ..
!
!   =====================================================================
!
!      .. Local Scalars ..
      DOUBLE PRECISION dtemp
      INTEGER i,ix,iy,m,mp1
!      ..
!      .. Intrinsic Functions ..
      INTRINSIC mod
!      ..
      ddot_d = 0.0d0
      dtemp = 0.0d0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) THEN
!
!         code for both increments equal to 1
!
!
!         clean-up loop
!
         m = mod(n,5)
         IF (m.NE.0) THEN
            DO i = 1,m
               dtemp = dtemp + dx(i)*dy(i)
            END DO
            IF (n.LT.5) THEN
               ddot_d=dtemp
            RETURN
            END IF
         END IF
         mp1 = m + 1
         DO i = mp1,n,5
          dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)*dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
         END DO
      ELSE
!
!         code for unequal increments or equal increments
!           not equal to 1
!
         ix = 1
         iy = 1
         IF (incx.LT.0) ix = (-n+1)*incx + 1
         IF (incy.LT.0) iy = (-n+1)*incy + 1
         DO i = 1,n
            dtemp = dtemp + dx(ix)*dy(iy)
            ix = ix + incx
            iy = iy + incy
         END DO
      END IF
      ddot_d = dtemp
      RETURN
!
!      End of DDOT
!
      END

 END MODULE
