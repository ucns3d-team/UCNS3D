MODULE LAPCK
IMPLICIT NONE
CONTAINS
SUBROUTINE dgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB, BETA,C,LDC)
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

END MODULE LAPCK

PROGRAM PROTOTYPE1
USE MPI
USE LAPCK
IMPLICIT NONE


INTEGER:: N !THE NUMBER OF RANK THAT I HAVE(FOR EACH PROCESSOR!
INTEGER:: ISIZE !THE TOTAL NUMBER OF RANKS(SIZE OF)
! INTEGER:: ICOMMUNICATOR !THE COMMUNICATOR OF COMM_wORLD
INTEGER::IERROR,provided
REAL::alpha,beta
INTEGER::STATUS(MPI_STATUS_SIZE)
INTEGER::I,J,K,L,M,IDEGFREE,stencil_local2,ll
INTEGER::KMAXE,KMAXN,NOF_VARIABLES,DOF,IMAX,STENCIL_LOCAL
REAL::R_l,R_KMAXE,r_j,r_k
REAL,ALLOCATABLE,DIMENSION(:)::LEFTV

TYPE LOCAL_RECON3
	INTEGER::LOCAL,MRF,G0
	REAL,ALLOCATABLE,DIMENSION(:)::cond  !dummy variable used for gradient approximation estimation 
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::IHEXG !GLOBAL INDEX OF CELLS
	REAL,ALLOCATABLE,DIMENSION(:,:,:)::INVMAT,SOL,MATRIX_1 ! (VAR, DIM, I_FACE)
END TYPE LOCAL_RECON3
TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:)::ILOCAL_RECON3
LOGICAL :: TRANSA = .FALSE.

!$OMP THREADPRIVATE(LEFTV,STENCIL_LOCAL)


CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,PROVIDED,IERROR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISIZE,IERROR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,N,IERROR)

!FOR ALL MPI PROCESSES
KMAXE=10000
NOF_VARIABLES=5
IDEGFREE=10
IMAX=20
LL=4


!NOW ALLOCATE A THREAD PRIVATE VARIABLE
!$OMP PARALLEL DEFAULT(SHARED)
ALLOCATE(LEFTV(1:NOF_VARIABLES))

!$OMP END PARALLEL 
!$OMP BARRIER

ALLOCATE (ILOCAL_RECON3(1:KMAXE))

DO I=1,KMAXE
    
    
    IF (I.LT.500)THEN
    STENCIL_LOCAL=5
    ELSE
    STENCIL_LOCAL=4
    END IF
    ALLOCATE(ILOCAL_RECON3(I)%INVMAT(1:IDEGFREE,1:IMAX,1:STENCIL_LOCAL));
    ALLOCATE(ILOCAL_RECON3(I)%MATRIX_1(1:imax,1:nof_variables,1:STENCIL_LOCAL));
    ALLOCATE(ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,1:STENCIL_LOCAL));
    
    DO J=1,IDEGFREE
            DO K=1,IMAX
                DO l=1,STENCIL_LOCAL
                    r_j=J
                    r_k=k
                    r_l=l
                    ILOCAL_RECON3(I)%INVMAT(j,k,l)=ATAN(R_j/R_k)*r_l
                END DO
            END DO
    END DO

    DO J=1,Imax
        DO K=1,nof_variables
            DO l=1,STENCIL_LOCAL
                r_j=J
                r_k=k
                r_l=l
                ILOCAL_RECON3(I)%MATRIX_1(j,k,l)=ATAN(R_J/R_k)*r_l
            END DO
        END DO
    END DO

END DO




!!$OMP PARALLEL DEFAULT(SHARED)
!!$OMP DO
!$OMP target teams distribute parallel do
do i=1,kmaxe
        IF (I.LT.500)THEN
            stencil_local2=5
        ELSE
            stencil_local2=4
        END IF
        do ll=1,stencil_local2
CALL DGEMM(TRANSA,TRANSA,IDEGFREE,nof_variables,imax,&
         ALPHA,ILOCAL_RECON3(I)%INVMAT(1:IDEGFREE,1:imax,LL),&
         IDEGFREE,ILOCAL_RECON3(I)%MATRIX_1(1:imax,1:nof_variables,ll),&
imax,BETA,ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,ll),IDEGFREE)
        ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,ll)=ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,ll)+i*2
        end do
end do
!$OMP END target teams distribute parallel do
!!$OMP END DO
!!$OMP END PARALLEL 

if (n.eq.0)then
!$OMP MASTER
        OPEN(20,FILE='OUTPUT.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
        do i=1,kmaxe
        write(20,*)ilocal_recon3(i)%sol(:,:,:)
        end do
        close(20)
!$OMP END MASTER
!$OMP BARRIER
end if


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
CALL MPI_FINALIZE(IERROR)
    
END PROGRAM PROTOTYPE1
