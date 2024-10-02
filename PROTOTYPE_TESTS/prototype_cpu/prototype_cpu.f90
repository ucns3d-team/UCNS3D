PROGRAM PROTOTYPE1
USE MPI
IMPLICIT NONE


INTEGER:: N !THE NUMBER OF RANK THAT I HAVE(FOR EACH PROCESSOR!
INTEGER:: ISIZE !THE TOTAL NUMBER OF RANKS(SIZE OF)
! INTEGER:: ICOMMUNICATOR !THE COMMUNICATOR OF COMM_wORLD
INTEGER::IERROR,provided,alpha,beta
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




!$OMP PARALLEL DEFAULT(SHARED)
!$OMP DO
do i=1,kmaxe
        IF (I.LT.500)THEN
            stencil_local2=5
        ELSE
            stencil_local2=4
        END IF
        do ll=1,stencil_local2
CALL DGEMM('N','N',IDEGFREE,nof_variables,imax,&
         ALPHA,ILOCAL_RECON3(I)%INVMAT(1:IDEGFREE,1:imax,LL),&
         IDEGFREE,ILOCAL_RECON3(I)%MATRIX_1(1:imax,1:nof_variables,ll),&
imax,BETA,ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,ll),IDEGFREE)
        ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,ll)=ILOCAL_RECON3(I)%SOL(1:IDEGFREE,1:nof_variables,ll)+i*2
        end do
end do
!$OMP END DO
!$OMP END PARALLEL 

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
