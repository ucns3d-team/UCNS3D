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

     SUBROUTINE dgemm_d(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC)
    !$omp declare target
   !
   !  -- Reference BLAS level3 routine --
   !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
   !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
   !
   !     .. Scalar Arguments ..
         DOUBLE PRECISION ALPHA,BETA
         INTEGER K,LDA,LDB,LDC,M,N
         LOGICAL TRANSA,TRANSB
   !     ..
   !     .. Array Arguments ..
         DOUBLE PRECISION A(LDA,*),B(LDB,*),C(LDC,*)
   !     ..
   !
   !  =====================================================================
   !
   !     ..
   !     ..
   !     .. Intrinsic Functions ..
         INTRINSIC max
   !     ..
   !     .. Local Scalars ..
         DOUBLE PRECISION TEMP
         INTEGER I,INFO,J,L,NROWA,NROWB
         LOGICAL NOTA,NOTB
   !     ..
   !     .. Parameters ..
         DOUBLE PRECISION ONE,ZERO
         parameter(one=1.0d+0,zero=0.0d+0)
   !     ..
   !
   !     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
   !     transposed and set  NROWA and NROWB  as the number of rows of  A
   !     and  B  respectively.
   !
         nota = .NOT. transa
         notb = .NOT. transb
         IF (nota) THEN
             nrowa = m
         ELSE
             nrowa = k
         END IF
         IF (notb) THEN
             nrowb = k
         ELSE
             nrowb = n
         END IF
   
   !
   !     Quick return if possible.
   !
         IF ((m.EQ.0) .OR. (n.EQ.0) .OR.   (((alpha.EQ.zero).OR. (k.EQ.0)).AND. (beta.EQ.one))) RETURN
   !
   !     And if  alpha.eq.zero.
   !
         IF (alpha.EQ.zero) THEN
             IF (beta.EQ.zero) THEN
                 DO 20 j = 1,n
                     DO 10 i = 1,m
                         c(i,j) = zero
      10             CONTINUE
      20         CONTINUE
             ELSE
                 DO 40 j = 1,n
                     DO 30 i = 1,m
                         c(i,j) = beta*c(i,j)
      30             CONTINUE
      40         CONTINUE
             END IF
             RETURN
         END IF
   !
   !     Start the operations.
   !
         IF (notb) THEN
             IF (nota) THEN
   !
   !           Form  C := alpha*A*B + beta*C.
   !
                 DO 90 j = 1,n
                     IF (beta.EQ.zero) THEN
                         DO 50 i = 1,m
                             c(i,j) = zero
      50                 CONTINUE
                     ELSE IF (beta.NE.one) THEN
                         DO 60 i = 1,m
                             c(i,j) = beta*c(i,j)
      60                 CONTINUE
                     END IF
                     DO 80 l = 1,k
                         temp = alpha*b(l,j)
                         DO 70 i = 1,m
                             c(i,j) = c(i,j) + temp*a(i,l)
      70                 CONTINUE
      80             CONTINUE
      90         CONTINUE
             ELSE
   !
   !           Form  C := alpha*A**T*B + beta*C
   !
                 DO 120 j = 1,n
                     DO 110 i = 1,m
                         temp = zero
                         DO 100 l = 1,k
                             temp = temp + a(l,i)*b(l,j)
     100                 CONTINUE
                         IF (beta.EQ.zero) THEN
                             c(i,j) = alpha*temp
                         ELSE
                             c(i,j) = alpha*temp + beta*c(i,j)
                         END IF
     110             CONTINUE
     120         CONTINUE
             END IF
         ELSE
             IF (nota) THEN
   !
   !           Form  C := alpha*A*B**T + beta*C
   !
                 DO 170 j = 1,n
                     IF (beta.EQ.zero) THEN
                         DO 130 i = 1,m
                             c(i,j) = zero
     130                 CONTINUE
                     ELSE IF (beta.NE.one) THEN
                         DO 140 i = 1,m
                             c(i,j) = beta*c(i,j)
     140                 CONTINUE
                     END IF
                     DO 160 l = 1,k
                         temp = alpha*b(j,l)
                         DO 150 i = 1,m
                             c(i,j) = c(i,j) + temp*a(i,l)
     150                 CONTINUE
     160             CONTINUE
     170         CONTINUE
             ELSE
   !
   !           Form  C := alpha*A**T*B**T + beta*C
   !
                 DO 200 j = 1,n
                     DO 190 i = 1,m
                         temp = zero
                         DO 180 l = 1,k
                             temp = temp + a(l,i)*b(j,l)
     180                 CONTINUE
                         IF (beta.EQ.zero) THEN
                             c(i,j) = alpha*temp
                         ELSE
                             c(i,j) = alpha*temp + beta*c(i,j)
                         END IF
     190             CONTINUE
     200         CONTINUE
             END IF
         END IF
   !
         RETURN
   !
   !     End of DGEMM
   !
   
END

SUBROUTINE dgemv_dn(M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!$omp declare target
  !
  !  -- Reference BLAS level2 routine --
  !  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  !
  !     .. Scalar Arguments ..
        DOUBLE PRECISION ALPHA,BETA
        INTEGER INCX,INCY,LDA,M,N
  !     ..
  !     .. Array Arguments ..
        DOUBLE PRECISION A(LDA,*),X(*),Y(*)
  !     ..
  !
  !  =====================================================================
  !
  !     .. Parameters ..
        DOUBLE PRECISION ONE,ZERO
        parameter(one=1.0d+0,zero=0.0d+0)
  !     ..
  !     .. Local Scalars ..
        DOUBLE PRECISION TEMP
        INTEGER I,IX,IY,J,JX,JY,KX,KY,LENX,LENY
  !     ..
  !     ..
  !     ..
  !     .. Intrinsic Functions ..
        INTRINSIC max
  !     ..
  
  !
  !     Quick return if possible.
  !
        IF ((m.EQ.0) .OR. (n.EQ.0) .OR.    ((alpha.EQ.zero).AND. (beta.EQ.one))) RETURN
  !
  !     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
  !     up the start points in  X  and  Y.
  !
            lenx = n
            leny = m
        IF (incx.GT.0) THEN
            kx = 1
        ELSE
            kx = 1 - (lenx-1)*incx
        END IF
        IF (incy.GT.0) THEN
            ky = 1
        ELSE
            ky = 1 - (leny-1)*incy
        END IF
  !
  !     Start the operations. In this version the elements of A are
  !     accessed sequentially with one pass through A.
  !
  !     First form  y := beta*y.
  !
        IF (beta.NE.one) THEN
            IF (incy.EQ.1) THEN
                IF (beta.EQ.zero) THEN
                    DO 10 i = 1,leny
                        y(i) = zero
     10             CONTINUE
                ELSE
                    DO 20 i = 1,leny
                        y(i) = beta*y(i)
     20             CONTINUE
                END IF
            ELSE
                iy = ky
                IF (beta.EQ.zero) THEN
                    DO 30 i = 1,leny
                        y(iy) = zero
                        iy = iy + incy
     30             CONTINUE
                ELSE
                    DO 40 i = 1,leny
                        y(iy) = beta*y(iy)
                        iy = iy + incy
     40             CONTINUE
                END IF
            END IF
        END IF
        IF (alpha.EQ.zero) RETURN
  !
  !        Form  y := alpha*A*x + y.
  !
            jx = kx
            IF (incy.EQ.1) THEN
                DO 60 j = 1,n
                    temp = alpha*x(jx)
                    DO 50 i = 1,m
                        y(i) = y(i) + temp*a(i,j)
     50             CONTINUE
                    jx = jx + incx
     60         CONTINUE
            ELSE
                DO 80 j = 1,n
                    temp = alpha*x(jx)
                    iy = ky
                    DO 70 i = 1,m
                        y(iy) = y(iy) + temp*a(i,j)
                        iy = iy + incy
     70             CONTINUE
                    jx = jx + incx
     80         CONTINUE
            END IF
  !
        RETURN
  !
  !     End of DGEMV
  !
      END

 END MODULE
