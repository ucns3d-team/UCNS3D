MODULE LIBRARY
!> @brief
!> This module includes all the subroutines related to stencils and neighbours establishing
USE MPIINFO
USE DECLARATION
USE IO
use TRANSFORM
IMPLICIT NONE
!**************************DEVELOPED BY PANAGIOTIS TSOUTSANIS & ANTONIS FOIVOS ANTONIADIS**************************!
!*****************************************************************UCNS3D*******************************************!
CONTAINS

subroutine tolerances
!> @brief
!> This subroutine specifies the tolerances and other fixed numbers used throughout the code
implicit none

tolsmall=1.0e-8
TOLBIG=1.0E+13
oo2=1.0D0/2.0D0
zero=0.0D0
! PI=(ACOS(zero))*2
PI=4.0D0*ATAN(1.0D0)

end subroutine tolerances



FUNCTION DETERMINA(EIGVL)
!> @brief
!> This function computes the determinant of a matrix
IMPLICIT NONE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::EIGVL
REAL,DIMENSION(5,5)::MATRIX
REAL::DETERMINA
INTEGER::I,J,K,L,NN
REAL::M,TEMP
LOGICAL :: DetExists = .TRUE.
matrix=EIGVL
NN=5
    l = 1
    !Convert to upper triangular form
    DO k = 1, nN-1
        IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, nN
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, nN
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO
                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO
            IF (DetExists .EQV. .FALSE.) THEN
                DETERMINA = 0
                return
            END IF
        ENDIF
        DO j = k+1, nN
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, nN
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
        END DO
    END DO
   
    !Calculate determinant by finding product of diagonal elements
    DETERMINA= l
    DO i = 1, nN
        DETERMINA = DETERMINA * matrix(i,i)
    END DO








END FUNCTION DETERMINA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !---------------------------------------------------------------------------------------------!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !---------------------------------------------------------------------------------------------!


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!

SUBROUTINE TIMERS(N,CPUX1,CPUX2,CPUX3,CPUX4,CPUX5,CPUX6,TIMEX1,TIMEX2,TIMEX3,TIMEX4,TIMEX5,TIMEX6)
!> @brief
!> This subroutine establishes the timers
REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)::CPUX1,CPUX2,CPUX3,CPUX4,CPUX5,CPUX6
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::TIMEX1,TIMEX2,TIMEX3,TIMEX4,TIMEX5,TIMEX6
INTEGER,INTENT(IN)::N

TIMEX6(1)=CPUX6(1)-CPUX1(1)
TIMEX5(1)=CPUX6(1)-CPUX5(1)
TIMEX4(1)=CPUX4(1)-CPUX3(1)
TIMEX3(1)=CPUX3(1)-CPUX2(1)
TIMEX2(1)=CPUX2(1)-CPUX1(1)


END SUBROUTINE TIMERS
 
! !---------------------------------------------------------------------------------------------!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




















SUBROUTINE XMPIFIND(XMPIE,XMPIN,XMPIELRANK,XMPINRANK,IMAXE,IMAXN,NPROC)
!> @brief
!> This subroutine finds the number of elements in each process
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIELRANK
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPINRANK
INTEGER::I,J,K,L,M,ITEMR
INTEGER,INTENT(IN)::NPROC,IMAXE,IMAXN



K=0
	DO I=1,IMAXE
		IF (XMPIE(I).EQ.n)THEN
		K=K+1
		
		
		END IF
	END DO
	
	XMPIELRANK(n)=K



END SUBROUTINE XMPIFIND




SUBROUTINE GLOBALIST(N,XMPIE,XMPIL,XMPIELRANK,IMAXE,ISIZE,CENTERR,GLNEIGH,IELEM)
!> @brief
!> This subroutine establishes the connectivity within each process
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::IELEM
INTEGER,INTENT(IN)::IMAXE,N,ISIZE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
REAL,ALLOCATABLE,DIMENSION(:,:)::CENTERX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHTS,GLNEIGHTR
INTEGER,ALLOCATABLE,DIMENSION(:)::GLNEIGHTOT
real,ALLOCATABLE,DIMENSION(:,:)::centerTS,centerTR
INTEGER::I,J,K,KMAXE,ICPUID,KJ,TEMPI,TEMPT
REAL::XV,YC,ZC
	IF (N.EQ.0) then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)"global number of elements",imaxe
	CLOSE(63)
	END IF

IF (DIMENSIONA.EQ.3)THEN

if ((typesten.gt.1).or.(icompact.ge.1))then
 ALLOCATE(CENTERR(1:IMAXE,1:3))
end if
ALLOCATE(GLNEIGH(1:IMAXE,1:6))
if ((typesten.gt.1).or.(icompact.ge.1))then
  CENTERR(:,:)=-TOLBIG
end if
 GLNEIGH(:,:)=0
 

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

 KMAXE=XMPIELRANK(N)
 

 DO K=1,KMAXE
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	CALL COMPUTE_CENTRE3d(N,K)
	
 	CENTERR(IELEM(N,K)%IHEXGL,1)=CORDS(1)
 	CENTERR(IELEM(N,K)%IHEXGL,2)=CORDS(2)
 	CENTERR(IELEM(N,K)%IHEXGL,3)=CORDS(3)
	end if
 	DO J=1,IELEM(N,K)%IFCA
		

 		GLNEIGH(IELEM(N,K)%IHEXGL,J)=IELEM(N,K)%INEIGHG(J)

		
	END DO
END DO

	KJ=0
	do K=1,kmaxe
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,k)%INTERIOR.EQ.0)THEN
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0))THEN
				KJ=KJ+1
			END IF
			ELSE
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0).AND.(IELEM(N,K)%IBOUNDS(J).EQ.0))THEN
				KJ=KJ+1
			END IF



			END IF
			
		END DO
	END DO
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)KJ,"=NUMBER OF UNKNOWN NEIGHBOURS BEFORE COMMUNICATION"
	CLOSE(63)
	end if
	
	
	ICPUID=N

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)



  ALLOCATE(GLNEIGHTOT(2))
  ALLOCATE(GLNEIGHTS(1:KMAXE,1:7))
  
  if ((typesten.gt.1).or.(icompact.ge.1))then
   ALLOCATE(centerts(1:KMAXE,1:3))
  end if
  GLNEIGHTS(:,:)=0
  

  DO K=1,KMAXE
	  GLNEIGHTS(K,1)=IELEM(N,K)%IHEXGL
	  DO J=1,IELEM(N,K)%IFCA
	  GLNEIGHTS(K,1+J)=GLNEIGH(IELEM(N,K)%IHEXGL,J)
	  END DO
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:3)=CENTERR(IELEM(N,K)%IHEXGL,1:3)
	  end if
  end do




!  !first

	DO I=0,ISIZE-1
		IF (I.NE.N)THEN
			
			TEMPI=KMAXE
			CALL MPI_SENDRECV(tempi,1,MPI_INTEGER,I,ICPUID,&
			tempt,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			ALLOCATE(GLNEIGHTR(tempt,1:7))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  ALLOCATE(centertR(tempt,1:3))
			  end if
			glneightr(:,:)=0
			 if ((typesten.gt.1).or.(icompact.ge.1))then
			centertR(:,:)=0.0
			end if
			TEMPI=KMAXE
			
			CALL MPI_SENDRECV(GLNEIGHTS(1:TEMPI,1:7),TEMPI*7,MPI_INTEGER,I,ICPUID,&
			GLNEIGHTR(1:TEMPT,1:7),TEMPT*7,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			CALL MPI_SENDRECV(CENTERts(1:TEMPI,1:3),TEMPI*3,MPI_DOUBLE_PRECISION,I,ICPUID,&
			CENTERtr(1:TEMPt,1:3),TEMPt*3,MPI_DOUBLE_PRECISION,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			end if
			
			do k=1,tempt
				      do j=1,6
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:3)=CENTERtr(k,1:3)
					    end if
					  end if
				      end do
				      
			end do
			deALLOCATE(GLNEIGHTR)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deALLOCATE(centertR)
			  end if


			
		end if
	  end do

    deALLOCATE(GLNEIGHTOT)
  deALLOCATE(GLNEIGHTS)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deALLOCATE(centerts)
      end if


END IF
IF (DIMENSIONA.EQ.2)THEN

if ((typesten.gt.1).or.(icompact.ge.1))then
 ALLOCATE(CENTERR(1:IMAXE,1:2))
end if
ALLOCATE(GLNEIGH(1:IMAXE,1:4))
if ((typesten.gt.1).or.(icompact.ge.1))then
  CENTERR(:,:)=-TOLBIG
end if
 GLNEIGH(:,:)=0
 

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

 KMAXE=XMPIELRANK(N)
 DO K=1,KMAXE
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	CALL COMPUTE_CENTRE2d(N,K)
	
 	CENTERR(IELEM(N,K)%IHEXGL,1)=CORDS(1)
 	CENTERR(IELEM(N,K)%IHEXGL,2)=CORDS(2)
 	
	end if
 	DO J=1,IELEM(N,K)%IFCA
		GLNEIGH(IELEM(N,K)%IHEXGL,J)=IELEM(N,K)%INEIGHG(J)
	end do	

END DO
	KJ=0
	do K=1,kmaxe
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,k)%INTERIOR.EQ.0)THEN
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0))THEN
				KJ=KJ+1
			END IF
			ELSE
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0).AND.(IELEM(N,K)%IBOUNDS(J).EQ.0))THEN
				KJ=KJ+1
			END IF



			END IF
			
		END DO
	END DO
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)KJ,"=NUMBER OF UNKNOWN NEIGHBOURS BEFORE COMMUNICATION"
	CLOSE(63)
	end if
	
	
	ICPUID=N

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)



  ALLOCATE(GLNEIGHTOT(2))
  ALLOCATE(GLNEIGHTS(1:KMAXE,1:5))
  if ((typesten.gt.1).or.(icompact.ge.1))then
   ALLOCATE(centerts(1:KMAXE,1:2))
  end if
  GLNEIGHTS(:,:)=0
  DO K=1,KMAXE
	  GLNEIGHTS(K,1)=IELEM(N,K)%IHEXGL
	  DO J=1,IELEM(N,K)%IFCA
	  GLNEIGHTS(K,1+J)=GLNEIGH(IELEM(N,K)%IHEXGL,J)
	  END DO
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:2)=CENTERR(IELEM(N,K)%IHEXGL,1:2)
	  end if
  end do




!  !first

	DO I=0,ISIZE-1
		IF (I.NE.N)THEN
			
			TEMPI=KMAXE
			CALL MPI_SENDRECV(tempi,1,MPI_INTEGER,I,ICPUID,&
			tempt,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			ALLOCATE(GLNEIGHTR(tempt,1:5))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  ALLOCATE(centertR(tempt,1:2))
			  end if
			glneightr(:,:)=0
			 if ((typesten.gt.1).or.(icompact.ge.1))then
			centertR(:,:)=0.0
			end if
			TEMPI=KMAXE
			
			CALL MPI_SENDRECV(GLNEIGHTS(1:TEMPI,1:5),TEMPI*5,MPI_INTEGER,I,ICPUID,&
			GLNEIGHTR(1:TEMPT,1:5),TEMPT*5,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			CALL MPI_SENDRECV(CENTERts(1:TEMPI,1:2),TEMPI*2,MPI_DOUBLE_PRECISION,I,ICPUID,&
			CENTERtr(1:TEMPt,1:2),TEMPt*2,MPI_DOUBLE_PRECISION,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			end if
			
			do k=1,tempt
				      do j=1,4
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:2)=CENTERtr(k,1:2)
					    end if
					  end if
				      end do
				      
			end do
			deALLOCATE(GLNEIGHTR)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deALLOCATE(centertR)
			  end if


			
		end if
	  end do

    deALLOCATE(GLNEIGHTOT)
  deALLOCATE(GLNEIGHTS)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deALLOCATE(centerts)
      end if


END IF


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

END SUBROUTINE GLOBALIST

SUBROUTINE GLOBALIST2(N,XMPIE,XMPIL,XMPIELRANK,IMAXE,ISIZE,CENTERR,GLNEIGH,IELEM)
!> @brief
!> This subroutine establishes the connectivity with elements from other processes
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::IELEM
INTEGER,INTENT(IN)::IMAXE,N,ISIZE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
REAL,ALLOCATABLE,DIMENSION(:,:)::CENTERX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHTS,GLNEIGHTR
INTEGER,ALLOCATABLE,DIMENSION(:)::GLNEIGHTOT
real,ALLOCATABLE,DIMENSION(:,:)::centerTS,centerTR
INTEGER::I,J,K,KMAXE,ICPUID,KJ,TEMPI,tempt
REAL::XV,YC,ZC


 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

  IF (DIMENSIONA.EQ.3)THEN

 KMAXE=XMPIELRANK(N)
 DO K=1,KMAXE
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	CALL COMPUTE_CENTRE3d(N,K)
	
 	CENTERR(IELEM(N,K)%IHEXGL,1)=CORDS(1)
 	CENTERR(IELEM(N,K)%IHEXGL,2)=CORDS(2)
 	CENTERR(IELEM(N,K)%IHEXGL,3)=CORDS(3)
	end if
 	DO J=1,IELEM(N,K)%IFCA
		

 		GLNEIGH(IELEM(N,K)%IHEXGL,J)=IELEM(N,K)%INEIGHG(J)

		
	END DO
END DO

	KJ=0
	do K=1,kmaxe
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,k)%INTERIOR.EQ.0)THEN
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0))THEN
				KJ=KJ+1
			END IF
			ELSE
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0).AND.(IELEM(N,K)%IBOUNDS(J).EQ.0))THEN
				KJ=KJ+1
			END IF



			END IF
			
		END DO
	END DO
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)KJ,"NUMBER OF UNKNOWN NEIGHBOURS AFTER COMMUNICATION"
      
	CLOSE(63)
	end if
	
	
	ICPUID=N
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

 





	    
 ALLOCATE(GLNEIGHTOT(2))
  ALLOCATE(GLNEIGHTS(1:KMAXE,1:7))
      if ((typesten.gt.1).or.(icompact.ge.1))then
   ALLOCATE(centerts(1:KMAXE,1:3))
      end if
  GLNEIGHTS(:,:)=0
  DO K=1,KMAXE
	  GLNEIGHTS(K,1)=IELEM(N,K)%IHEXGL
	  DO J=1,IELEM(N,K)%IFCA
	  GLNEIGHTS(K,1+J)=GLNEIGH(IELEM(N,K)%IHEXGL,J)
	  END DO
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:3)=CENTERR(IELEM(N,K)%IHEXGL,1:3)
	  end if
  end do




 
	DO I=0,ISIZE-1
		IF (I.NE.N)THEN


			TEMPI=KMAXE
			CALL MPI_SENDRECV(tempi,1,MPI_INTEGER,I,ICPUID,&
			tempt,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			ALLOCATE(GLNEIGHTR(tempt,1:7))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
 			  ALLOCATE(centertR(tempt,1:3))
			  centertR(:,:)=0.0
			  end if
			glneightr(:,:)=0
! 			
			
			CALL MPI_SENDRECV(GLNEIGHTS(1:TEMPI,1:7),TEMPI*7,MPI_INTEGER,I,ICPUID,&
			GLNEIGHTR(1:TEMPT,1:7),TEMPT*7,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			CALL MPI_SENDRECV(CENTERts(1:TEMPI,1:3),TEMPI*3,MPI_DOUBLE_PRECISION,I,ICPUID,&
			CENTERtr(1:TEMPt,1:3),TEMPt*3,MPI_DOUBLE_PRECISION,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			end if

			do k=1,tempt
				      do j=1,6
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:3)=CENTERtr(k,1:3)
					    end if
					  end if
				      end do
				      
			end do
			deALLOCATE(GLNEIGHTR)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deALLOCATE(centertR)
			  end if



			
		end if
	  end do
    deALLOCATE(GLNEIGHTOT)
  deALLOCATE(GLNEIGHTS)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deALLOCATE(centerts)
      end if

ELSE



 KMAXE=XMPIELRANK(N)
 DO K=1,KMAXE
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	CALL COMPUTE_CENTRE2d(N,K)
	
 	CENTERR(IELEM(N,K)%IHEXGL,1)=CORDS(1)
 	CENTERR(IELEM(N,K)%IHEXGL,2)=CORDS(2)
 	
	end if
 	DO J=1,IELEM(N,K)%IFCA
		

 		GLNEIGH(IELEM(N,K)%IHEXGL,J)=IELEM(N,K)%INEIGHG(J)

		
	END DO
END DO

	KJ=0
	do K=1,kmaxe
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,k)%INTERIOR.EQ.0)THEN
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0))THEN
				KJ=KJ+1
			END IF
			ELSE
			IF ((GLNEIGH(IELEM(N,K)%IHEXGL,J).EQ.0).AND.(IELEM(N,K)%IBOUNDS(J).EQ.0))THEN
				KJ=KJ+1
			END IF



			END IF
			
		END DO
	END DO
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)KJ,"NUMBER OF UNKNOWN NEIGHBOURS AFTER COMMUNICATION"
	CLOSE(63)
	end if
	
	
	ICPUID=N
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

 





	    
 ALLOCATE(GLNEIGHTOT(2))
  ALLOCATE(GLNEIGHTS(1:KMAXE,1:5))
      if ((typesten.gt.1).or.(icompact.ge.1))then
   ALLOCATE(centerts(1:KMAXE,1:2))
      end if
  GLNEIGHTS(:,:)=0
  DO K=1,KMAXE
	  GLNEIGHTS(K,1)=IELEM(N,K)%IHEXGL
	  DO J=1,IELEM(N,K)%IFCA
	  GLNEIGHTS(K,1+J)=GLNEIGH(IELEM(N,K)%IHEXGL,J)
	  END DO
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:2)=CENTERR(IELEM(N,K)%IHEXGL,1:2)
	  end if
  end do




 
	DO I=0,ISIZE-1
		IF (I.NE.N)THEN


			TEMPI=KMAXE
			CALL MPI_SENDRECV(tempi,1,MPI_INTEGER,I,ICPUID,&
			tempt,1,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			ALLOCATE(GLNEIGHTR(tempt,1:5))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
 			  ALLOCATE(centertR(tempt,1:2))
			  centertR(:,:)=0.0
			  end if
			glneightr(:,:)=0
! 			
			
			CALL MPI_SENDRECV(GLNEIGHTS(1:TEMPI,1:5),TEMPI*5,MPI_INTEGER,I,ICPUID,&
			GLNEIGHTR(1:TEMPT,1:5),TEMPT*5,MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			CALL MPI_SENDRECV(CENTERts(1:TEMPI,1:2),TEMPI*2,MPI_DOUBLE_PRECISION,I,ICPUID,&
			CENTERtr(1:TEMPt,1:2),TEMPt*2,MPI_DOUBLE_PRECISION,I,I,MPI_COMM_WORLD,STATUS,IERROR)
			end if

			do k=1,tempt
				      do j=1,4
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:2)=CENTERtr(k,1:2)
					    end if
					  end if
				      end do
				      
			end do
			deALLOCATE(GLNEIGHTR)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deALLOCATE(centertR)
			  end if



			
		end if
	  end do
    deALLOCATE(GLNEIGHTOT)
  deALLOCATE(GLNEIGHTS)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deALLOCATE(centerts)
      end if

END IF


END SUBROUTINE GLOBALIST2




SUBROUTINE GLOBALISTX(N,XMPIE,XMPIL,XMPIELRANK,IMAXE,ISIZE,CENTERR,GLNEIGH,IELEM)
!> @brief
!> This subroutine establishes the connectivity within each process without using a global list
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::IELEM
INTEGER,INTENT(IN)::IMAXE,N,ISIZE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
REAL,ALLOCATABLE,DIMENSION(:,:)::CENTERX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHTS,GLNEIGHTR
INTEGER,ALLOCATABLE,DIMENSION(:)::GLNEIGHTOT
real,ALLOCATABLE,DIMENSION(:,:)::centerTS,centerTR
INTEGER::I,J,K,KMAXE,ICPUID,KJ,TEMPI,TEMPT
REAL::XV,YC,ZC
KMAXE=XMPIELRANK(N)
KJ=0
	do K=1,kmaxe
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,k)%INTERIOR.EQ.0)THEN
			  IF (IELEM(N,K)%INEIGHG(J).EQ.0)THEN
				  KJ=KJ+1
			  END IF
			ELSE
			  IF ((IELEM(N,K)%INEIGHG(J).EQ.0).AND.(IELEM(N,K)%IBOUNDS(J).EQ.0))THEN
				  KJ=KJ+1
			  END IF
			END IF
		END DO
	END DO
	
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)KJ,"=NUMBER OF UNKNOWN NEIGHBOURS BEFORE COMMUNICATION"
	CLOSE(63)
	end if



END SUBROUTINE GLOBALISTX


SUBROUTINE GLOBALISTX2(N,XMPIE,XMPIL,XMPIELRANK,IMAXE,ISIZE,CENTERR,GLNEIGH,IELEM)
!> @brief
!> This subroutine establishes the connectivity with elements from other processes without a global list
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::IELEM
INTEGER,INTENT(IN)::IMAXE,N,ISIZE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
REAL,ALLOCATABLE,DIMENSION(:,:)::CENTERX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHX
INTEGER,ALLOCATABLE,DIMENSION(:,:)::GLNEIGHTS,GLNEIGHTR
INTEGER,ALLOCATABLE,DIMENSION(:)::GLNEIGHTOT
real,ALLOCATABLE,DIMENSION(:,:)::centerTS,centerTR
INTEGER::I,J,K,KMAXE,ICPUID,KJ,TEMPI,TEMPT,IFST2
integer:: n_requests
integer,allocatable, dimension(:) :: requests
integer::ifst,IMAX_CPU, IMAX_CPUT,SUMCENTRAL,SUMCENTRAL_T,APPROX
REAL::XV,YC,ZC
KMAXE=XMPIELRANK(N)
KJ=0
	DO K=1,kmaxe
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,k)%INTERIOR.EQ.0)THEN
			  IF (IELEM(N,K)%INEIGHG(J).EQ.0)THEN
				  KJ=KJ+1
			  END IF
			ELSE
			  IF ((IELEM(N,K)%INEIGHG(J).EQ.0).AND.(IELEM(N,K)%IBOUNDS(J).EQ.0))THEN
				  KJ=KJ+1
			  END IF
			END IF
		END DO
	END DO
	
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)KJ,"=NUMBER OF UNKNOWN NEIGHBOURS AFTER COMMUNICATION"
	CLOSE(63)
	end if
	
	allocate(cand(0:isize-1));cand=0		!1 ALLOCATE(CAND)
	DO K=1,kmaxe
		
		DO J=1,IELEM(N,K)%IFCA
			IF (IELEM(N,K)%INEIGHG(J).gt.0)then
! 			
			  if(xmpie(IELEM(N,K)%INEIGHG(J)).ne.n) THEN
				 cand(xmpie(IELEM(N,K)%INEIGHG(J)))=cand(xmpie(IELEM(N,K)%INEIGHG(J)))+1   
			  END IF
			end if
			 
		END DO
		
	END DO
! 	
	!first establish the cpus that are needed!
	
	
	
	IF (IPERIODICITY.EQ.1)THEN
	  
	SUMCENTRAL_T=MIN(IMAXE,ISELEM*ISIZE)
	APPROX=min(isize,(SUMCENTRAL_T/KMAXE)*iextend)
	if (n.eq.0)then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)APPROX,"=NUMBER OF CPUS REQUIRED FOR PERIODIC BOUNDARY CONDITIONS AND STENCILS"
	CLOSE(63)
	end if
	END IF
	
	
	
	
	ifst=0
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  ifst=ifst+1
! 		      
		  end if
	end do
	
	
	
	
	!first allocate memory for these additional ones
	allocate(candS(0:isize-1),CANDR(0:isize-1));candS=0;CANDR=0.0   !2 ALLOCATE(CANDS,CANDR,NEIX)
	allocate(neix1(ifst))
	ifst=0
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  ifst=ifst+1
		      neix1(ifst)%cpu=cand(i)
		      CANDS(I)=IFST
		  end if
	end do
	
	!now send to cpus of cpus!
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  CANDS(I)=IFST
		  END IF
	END DO
	CALL MPI_BARRIER(mpi_comm_world,ierror)
	
	
	
	n_requests = 0
      allocate(requests(2*ifst))					!3 ALLOCATE(REQUESTS)
	
	do i=0,isize-1
	
		  IF (cand(i).gt.0)THEN
		      n_requests = n_requests + 1
			      CALL MPI_ISEND(                                                     &
			    CANDS(I), 								& !sendbuf
			    1, MPI_INTEGER,       						& !sendcount, sendtype
			    I, 0,                                       			 & !destination, tag
			    MPI_COMM_WORLD, requests(n_requests), ierror                       & !communicator, request handle, error
			)
	
		      n_requests = n_requests + 1
			      CALL MPI_IRECV(                                                     &
				  CANDR(I),    & !recvbuf
				  1, MPI_INTEGER,          & !recvcount, recvtype
				  I, 0,                                        & !source, tag
				  MPI_COMM_WORLD, requests(n_requests), ierror                     & !communicator, request handle, error
			      )
	
	
		  END IF
	
	END DO


	CALL MPI_WAITALL(n_requests, requests, MPI_STATUSES_IGNORE, ierror)
	
	
	deallocate (requests)							!33 DALLOCATE(REQUESTS)
	
	
	
	  imax_cpu=0
	DO I=0,ISIZE-1
		IF (cand(i).gt.0)THEN
		imax_cpu=max(imax_cpu,cands(i),candr(i))
		END IF
	END DO
	imax_cpu=imax_cpu*1.3
	
	CALL MPI_ALLREDUCE(IMAX_CPU,IMAX_CPUT,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
	
	
	allocate(candxs(1:ifst,1:imax_cpuT),candxr(1:ifst,1:imax_cpuT))		!4 ALLOCATE(CANDXS,CANDXR)
	CANDXS=-1
	candxr=-1
! 	
	
	
	ifst=0
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  ifst=ifst+1
		  candxs(:,ifst)=i
		  end if
	end do
	
	
	
	
     
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
 	IFST=0
	do i=0,isize-1
	
		  IF (cand(i).gt.0)THEN
		      
		      
		      IFST=IFST+1
 		      CALL MPI_SENDRECV(CANDXS(ifst:ifst,1:IMAX_CPUT)&
,IMAX_CPUT,MPI_INTEGER,i,&
      n,CANDXr(ifst:ifst,1:IMAX_CPUT),&
IMAX_CPUT,MPI_INTEGER,I,&
      I,MPI_COMM_WORLD,STATUS,IERROR)
	
	
		  END IF
	
	END DO


	
	
 	DO I=1,IFST
	    DO J=1,IMAX_CPUT
	    IF ((CANDXR(I,J).GE.0).AND.(CANDXR(I,J).NE.N))THEN
	     CAND(CANDXR(I,J))=1
	    END IF
	    END DO
 	END DO
	
	IFST2=0
	DO I=0,ISIZE-1
	    IF (CAND(I).GT.0)THEN
	    IFST2=IFST2+1
! 	    
	    END IF
	END DO
	
	
	
	
	
	
	
	
	
	
	IFST2=0
	DO I=0,ISIZE-1
	    IF (CAND(I).GT.0)THEN
	    IFST2=IFST2+1
	    CAND(I)=ifst2
	    END IF
	END DO
	
	
	
	
!
	!now gather from all cpus the required ones
	
	IMAX_CPU=KMAXE
	IMAX_CPUT=0
	CALL MPI_ALLREDUCE(IMAX_CPU,IMAX_CPUT,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
	
	if (dimensiona.eq.3)then
	ALLOCATE(CAND2S(IMAX_CPUT,6))						!5 ALLOCATE(CAND2S,CAND2R)
	ALLOCATE(CAND2R(IFST2,IMAX_CPUT,6))
	CAND2S=-100
	CAND2R=-200
	else
	ALLOCATE(CAND2S(IMAX_CPUT,4))
	ALLOCATE(CAND2R(IFST2,IMAX_CPUT,4))
	CAND2S=-100
	CAND2R=-200
	
	
	end if
	DO I=1,KMAXE
	      DO J=1,IELEM(N,i)%IFCA
			  CAND2S(I,J)=IELEM(N,I)%INEIGHG(J)
	      END DO
	END DO
	
	
	

	
	
	if (dimensiona.eq.3)then
	IFST2=0
	
	do i=0,isize-1
	
		  IF (cand(i).gt.0)THEN
		  IFST2=ifst2+1
		   CALL MPI_SENDRECV(CAND2S(1:IMAX_CPUT,1:6)&
,IMAX_CPUT*6,MPI_INTEGER,i,&
      n,CAND2R(IFST2:IFST2,1:IMAX_CPUT,1:6),&
IMAX_CPUT*6,MPI_INTEGER,I,&
      I,MPI_COMM_WORLD,STATUS,IERROR)
		  END IF
	
	END DO
	
	else
	IFST2=0
	
	do i=0,isize-1
	
		  IF (cand(i).gt.0)THEN
		  IFST2=ifst2+1
		   CALL MPI_SENDRECV(CAND2S(1:IMAX_CPUT,1:4)&
,IMAX_CPUT*4,MPI_INTEGER,i,&
      n,CAND2R(IFST2:IFST2,1:IMAX_CPUT,1:4),&
IMAX_CPUT*4,MPI_INTEGER,I,&
      I,MPI_COMM_WORLD,STATUS,IERROR)
		  END IF
	
	END DO
	
	
	end if
	
! 	
	
	
	 if ((typesten.gt.1).or.(icompact.ge.1))then
	 
	if (dimensiona.eq.3)then
	ALLOCATE(xAND2S(IMAX_CPUT,3))						!6 ALLOCATE(XAND2S,XAND2R)
	ALLOCATE(xAND2R(IFST2,IMAX_CPUT,3))
	do k=1,kmaxe
	
	CALL COMPUTE_CENTRE3d(N,K)
	
 	xAND2S(k,1:3)=CORDS(1:3)
 	
	end do
	
	else
	ALLOCATE(xAND2S(IMAX_CPUT,2))
	ALLOCATE(xAND2R(IFST2,IMAX_CPUT,2))
	do k=1,kmaxe
	
	CALL COMPUTE_CENTRE2d(N,K)
	
 	xAND2S(k,1:2)=CORDS(1:2)
 	
	end do

	end if

	IFST2=0
	
	do i=0,isize-1
	
		  IF (cand(i).gt.0)THEN
		  IFST2=ifst2+1
		   CALL MPI_SENDRECV(xAND2S(1:IMAX_CPUT,1:dims)&
,IMAX_CPUT*dims,MPI_DOUBLE_PRECISION,i,&
      n,xAND2R(IFST2:IFST2,1:IMAX_CPUT,1:dims),&
IMAX_CPUT*dims,MPI_DOUBLE_PRECISION,I,&
      I,MPI_COMM_WORLD,STATUS,IERROR)
		  END IF
	
	END DO
	
	
	SUMCENTRAL=IMAX_CPUT*IFST2
	
	IF (N.EQ.0)THEN
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	WRITE(63,*)ISELEM,"=NUMBER OF TOTAL NEIGHBOURS NEEDED"
	WRITE(63,*)IMAX_CPUT,"=NUMBER OF NEIGHBOURS PER ADJACENT CPUS"
	WRITE(63,*)SUMCENTRAL,"=TOTAL NUMBER OF NEIGHBOURS PER ALL ADJACENT CPUS"
	CLOSE(63)
	END IF
	
	
	END IF
	CALL MPI_BARRIER(mpi_comm_world,ierror)
	
	
	
								    !NOW DEALLOCATE WHAT IS NOT NEEDED
	    DEALLOCATE(CANDS,CANDR,NEIX1,CANDXS,CANDXR)
								    
	
	
! 	
	    
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



END SUBROUTINE GLOBALISTX2



SUBROUTINE GLOBALDEA
!> @brief
!> This subroutine deallocates the memory used for establishing the connectivity
IMPLICIT NONE
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    if ((typesten.gt.1).or.(icompact.ge.1))then
	    IF (LOWMEM.EQ.0)deallocate(xand2s,xAND2R)
	    IF (LOWMEM.EQ.1)DEALLOCATE(CENTERR)
    end if
	    IF (LOWMEM.EQ.0)deallocate(cand2s,cand2r,cand)
	    IF (LOWMEM.EQ.1)DEALLOCATE(GLNEIGH)
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 END SUBROUTINE GLOBALDEA



SUBROUTINE XMPILOCAL
!> @brief
!> This subroutine establishes the local numbering of the cells
IMPLICIT NONE
INTEGER::I,KMAXE,count_block
INTEGER,ALLOCATABLE,DIMENSION(:)::XMPIC,BIN,VAL
ALLOCATE(XMPIC(0:isize-1))

XMPIC=0
KMAXE=XMPIELRANK(N)

ALLOCATE(VAL(1:KMAXE))
ALLOCATE(XGO(1:KMAXE))

DO I=1,KMAXE
XGO(I)=IELEM(N,I)%IHEXGL
END DO

IF (N.EQ.0)THEN
ALLOCATE(XMPI_RE(1:IMAXE))
xmpi_re=zero
END if

DO I=1,KMAXE
	XMPIL(IELEM(N,I)%IHEXGL)=I
END DO

	
do i=1,imaxe
      xmpic(xmpie(i))=xmpic(xmpie(i))+1
      xmpil(i)=xmpic(xmpie(i))
end do




DEALLOCATE(XMPIC)


ALLOCATE(XMPIALL(0:ISIZE-1),OFFSET(0:ISIZE-1))

XMPIALL=0

XMPIALL(N)=KMAXE



CALL MPI_ALLGATHER(kmaxe,1,MPI_INTEGER,XMPIALL,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)



OFFSET(0)=0
DO I=1,ISIZE-1
OFFSET(I)=OFFSET(I-1)+XMPIALL(I-1)
END DO




DO I=1,KMAXE
  VAL(I)=IELEM(N,I)%IHEXGL
END DO





call MPI_GATHERv(val(1:kmaxe),kmaxe,MPI_INTEGER,XMPI_RE,xmpiall,offset,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)






DEALLOCATE(VAL)


!  CHUNK_N=0
! DO I=2,KMAXE
!     IF (XGO(I).GT.XGO(I-1)+1)THEN
! 	CHUNK_N=CHUNK_N+1
!     END IF
! END DO
!   
!   WRITE(100+N,*)CHUNK_N,N,KMAXE
!   
! allocate(chunk_size(chunk_n))
! 
!  CHUNK_N=0
!  count_block=0
! DO I=2,KMAXE
!     IF (XGO(I).GT.XGO(I-1)+1)THEN
! 	CHUNK_N=CHUNK_N+1
! 	chunk_size(chunk_n)=XGO(I)
!     END IF
! END DO
! 
! do i=1,chunk_n
! write(100+n,*)i,chunk_size(i)
! end do

!  CHUNK_N=0
! do i=1,kmaxe-1
!       if (xgo(i+1).eq.xgo(i)+1)then
!       
!       else
!       CHUNK_N=CHUNK_N+1
!       
!       end if
! end do
! 
! 
! allocate(chunk_size(chunk_n))
!  count_block=0
! do i=1,kmaxe-1
!       if (xgo(i+1).eq.xgo(i)+1)then
!       count_block=count_block+1
!       else
!       CHUNK_N=CHUNK_N+1
!       
!       count_block=count_block+1
!       end if
! end do
!     





END SUBROUTINE XMPILOCAL



SUBROUTINE COUNT_WALLS
!> @brief
!> This subroutine allocates the appropriate memory for bounded walls indexing for writing files
IMPLICIT NONE
INTEGER::I,KMAXE,ILOOP,j
INTEGER,ALLOCATABLE,DIMENSION(:)::BIN,VAL
KMAXE=XMPIELRANK(N)
ILOOP=0
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
		  iloop=iloop+1
	      END IF
	  end if
	END DO
   end if
END DO



ALLOCATE(XMPIWALL(0:ISIZE-1),WOFFSET(0:ISIZE-1))

XMPIWALL=0

XMPIWALL(N)=ILOOP
CALL MPI_ALLGATHER(ILOOP,1,MPI_INTEGER,XMPIWALL,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)


WOFFSET(0)=0
DO I=1,ISIZE-1
WOFFSET(I)=WOFFSET(I-1)+XMPIWALL(I-1)
END DO

ALLOCATE(VAL(iloop))
IF (N.EQ.0)THEN
ALLOCATE(XMPI_wRE(1:totwalls))
END if


ILOOP=0

DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
		  iloop=iloop+1
		  VAL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
	      END IF
	  end if
	END DO
   end if
END DO



call MPI_GATHERv(val(1:iloop),iloop,MPI_INTEGER,XMPI_wRE,xmpiwall,woffset,MPI_INTEGER,0,MPI_COMM_WORLD,IERROR)




DEALLOCATE(VAL)





END SUBROUTINE COUNT_WALLS




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




SUBROUTINE INVERF(R,INVR,IVGT)
!> @brief
!> This subroutine inverts some special matrices
  IMPLICIT NONE
  INTEGER,INTENT(IN)::IVGT
  REAL,INTENT(IN) ::R(1:IVGT-1,1:IVGT-1)
  REAL,INTENT(OUT)::INVR(1:IVGT-1,1:IVGT-1)
  REAL::INVVVR2(1:IVGT-1,1:IVGT-1)
 INTEGER I,J,K,GT
 
  INVR = 0.d0
 GT=IVGT-1
  DO I=GT,1,-1
    INVR(i,i) = 1./R(i,i)
	do j=i+1,GT
	  INVR(i,j) = 0.
	  do k= 1,j-1
	   INVR(i,j) = INVR(i,j) - R(k,j)*INVR(i,k)
	  enddo
 	  INVR(i,j) =INVR(i,j) /R(j,j)
	enddo
  enddo

 End subroutine INVERF

SUBROUTINE ALLOCATEVECTORS(N,TRI,INVTRI,ROTVECT,VECTCO,VEIGL,VEIGR,RVEIGL,RVEIGR,EIGVL,EIGVR)
!> @brief
!> This subroutine allocates vectors and matrices frequently used
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::TRI
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INVTRI
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVECT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VECTCO
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VEIGL,VEIGR,RVEIGL,RVEIGR
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::EIGVL,EIGVR
ALLOCATE(TRI(5,5))
ALLOCATE(INVTRI(5,5))
ALLOCATE(EIGVL(5,5))
ALLOCATE(EIGVR(5,5))
ALLOCATE(VECTCO(5+turbulenceequations+passivescalar))
ALLOCATE(ROTVECT(5+turbulenceequations+passivescalar))
ALLOCATE(VEIGL(5))
ALLOCATE(VEIGR(5))
ALLOCATE(RVEIGL(5))
ALLOCATE(RVEIGR(5))
TRI=0.0
INVTRI=0.0
EIGVL=0.0
EIGVR=0.0
VECTCO=0.0
ROTVECT=0.0
VEIGL=0.0
VEIGR=0.0
RVEIGL=0.0
RVEIGR=0.0
END SUBROUTINE ALLOCATEVECTORS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR STENCIL!!!!!!!!!!!!!!!!!
! SUBROUTINE TIMECPU(TEMPCP)
! REAL,INTENT(INOUT)::TEMPCP
! CALL CPU_TIME(TEMPCP)
! END SUBROUTINE TIMECPU
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED INITIALLY TO ALLOCATE MEMORY FOR STENCIL!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	SUBROUTINE LOCALSTALLOCATION(N,XMPIELRANK,ILOCALSTENCIL,TYPESTEN,NUMNEIGHBOURS)
!> @brief
!> This subroutine memory for stencil allocation for all cells
	IMPLICIT NONE
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
	INTEGER,INTENT(IN)::NUMNEIGHBOURS
	INTEGER,INTENT(IN)::N
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,INTENT(IN)::TYPESTEN
	INTEGER::KMAXE,i
	KMAXE=XMPIELRANK(N)
	
	ALLOCATE (ILOCALSTENCIL(N:N,kmaxe,TYPESTEN,imaxdegfree+1))
	
	ILOCALSTENCIL(N:N,:,:,:)=0
	
	END SUBROUTINE LOCALSTALLOCATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED FOR DETERMINING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!THE NEIGHBOURING CELLS AT ALL DIRECTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!KIT IS INITIALLY USED FOR FIRST LEVEL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!NEIGHBOURS AND IT IS GOING TO BE MODIFIED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FIND_SHAPE(N,IMAXE,IESHAPE)
!> @brief
!> This subroutine finds the shape of each cell
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::IMAXE
	integer,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IESHAPE
	INTEGER::I1,I2,I3,I4,I5,I6,I7,I8,I,J,K,L,M,KK,INFO,FH,ll,mm
	CHARACTER(LEN=12)::CELFILE,PROC
	INTEGER,INTENT(IN)::N
	INFO=1
	FH=1
	
	if (BINIO.EQ.0)THEN
	 if (dimensiona.eq.3)then
		WRITE(PROC,FMT='(I10)') N
		CELFILE='GRID.cel'
		
		OPEN(8,FILE=CELFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	DO J=1,IMAXE	
		READ(8,*)I,I1,I2,I3,I4,I5,I6,I7,I8
		IF ((I5.NE.I6).AND.(I6.NE.I7).AND.(I7.NE.I8).AND.(I3.NE.I4))THEN
		IESHAPE(J)=1	!HEXAHEDRAL ELEMENT
		END IF
		IF ((I3.EQ.I4).AND.(I5.EQ.I6).AND.(I6.EQ.I7).AND.(I7.EQ.I8))THEN
		IESHAPE(J)=2	!TETRAHEDRAL ELEMENT
		END IF
		IF ((I3.NE.I4).AND.(I5.EQ.I6).AND.(I6.EQ.I7).AND.(I7.EQ.I8))THEN
		IESHAPE(J)=3	!PYRAMIDAL ELEMENT
		END IF
		IF ((I3.EQ.I4).AND.(I5.NE.I6).AND.(I7.EQ.I8))THEN
		IESHAPE(J)=4	!PRISM ELEMENT
		END IF
	END DO	
	else
	WRITE(PROC,FMT='(I10)') N
		CELFILE='GRID.cel'
		OPEN(8,FILE=CELFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	DO J=1,IMAXE	
		READ(8,*)I,I1,I2,I3,I4
		IF ((I3.NE.I4))THEN
		IESHAPE(J)=5	!quadrilateral ELEMENT
		else
		IESHAPE(J)=6	!Triangle eL ELEMENT
		END IF
		
	END DO	



	end if
	ELSE
	 if (dimensiona.eq.3)then
		WRITE(PROC,FMT='(I10)') N
		CELFILE='GRID.cel'
		
		OPEN(8,FILE=CELFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
	DO J=1,IMAXE	
		READ(8)I,I1,I2,I3,I4,I5,I6,I7,I8
		IF ((I5.NE.I6).AND.(I6.NE.I7).AND.(I7.NE.I8).AND.(I3.NE.I4))THEN
		IESHAPE(J)=1	!HEXAHEDRAL ELEMENT
		END IF
		IF ((I3.EQ.I4).AND.(I5.EQ.I6).AND.(I6.EQ.I7).AND.(I7.EQ.I8))THEN
		IESHAPE(J)=2	!TETRAHEDRAL ELEMENT
		END IF
		IF ((I3.NE.I4).AND.(I5.EQ.I6).AND.(I6.EQ.I7).AND.(I7.EQ.I8))THEN
		IESHAPE(J)=3	!PYRAMIDAL ELEMENT
		END IF
		IF ((I3.EQ.I4).AND.(I5.NE.I6).AND.(I7.EQ.I8))THEN
		IESHAPE(J)=4	!PRISM ELEMENT
		END IF
	END DO	
	else
	WRITE(PROC,FMT='(I10)') N
		CELFILE='GRID.cel'
		OPEN(8,FILE=CELFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
	DO J=1,IMAXE	
		READ(8)I,I1,I2,I3,I4
		IF ((I3.NE.I4))THEN
		IESHAPE(J)=5	!quadrilateral ELEMENT
		else
		IESHAPE(J)=6	!Triangle eL ELEMENT
		END IF
		
	END DO	



	end if
	
	
	
	END IF
	
	
	K=0
	L=0
	M=0
	KK=0
	ll=0
	mm=0
	DO J=1,IMAXE
	IF (IESHAPE(J).EQ.1)THEN
		K=K+1
	END IF
	IF (IESHAPE(J).EQ.2)THEN
		L=L+1
	END IF
	IF (IESHAPE(J).EQ.3)THEN
		M=M+1
	END IF
	IF (IESHAPE(J).EQ.4)THEN
		KK=KK+1
	END IF
	IF (IESHAPE(J).EQ.5)THEN
		ll=ll+1
	END IF
	IF (IESHAPE(J).EQ.6)THEN
		mm=mm+1
	END IF
	END DO
	

close(8)


IF (N.EQ.0)THEN
OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	
!-------------------FOR DEBUGGING ONLY -----------------------------------------!
	WRITE(63,*)"------------UNSTRUCTURED FLUID DYNAMICS SOLVER ------------"
	WRITE(63,*)"-----------------------GRID  TYPE--------------------------"
	IF (DIMENSIONA.EQ.3 )THEN
	WRITE(63,*)"-----------------------3D MODE--------------------------"
	ELSE
	WRITE(63,*)"-----------------------2D MODE--------------------------"
	END IF
	!IF ((L.GT.0).AND.((K.GT.0).OR.(M.GT.0).OR.(N.GT.0)))THEN
	WRITE(63,*)K,"HEXAHEDRAL ELEMENTS"
	WRITE(63,*)L,"TETRAHEDRAL ELEMENTS"
	WRITE(63,*)M,"PYRAMIDAL ELEMENTS"
	WRITE(63,*)KK,"PRISMATIC ELEMENTS"
	WRITE(63,*)ll,"QUADRILATERAL ELEMENTS"
	WRITE(63,*)mm,"TRIANGULAR ELEMENTS"
	!END IF
	
	

CLOSE(63)

If ((emetis .eq. 1) .or. (emetis .eq.2)) then
      ALLNODESGLOBALL=(K*8)+(L*4)+(M*5)+(KK*6)
      else
      ALLNODESGLOBALL=8*imaxe
    end if
    

END IF


	!-------------------FOR DEBUGGING ONLY -----------------------------------------!
END SUBROUTINE FIND_SHAPE
! ---------------------------------------------------------------------------------------------!
SUBROUTINE NEIGHBOURSS(N,IELEM,IMAXE,IMAXN,XMPIE,XMPIN,XMPIELRANK,RESTART,INODEr)
!> @brief
!> This subroutine finds the neighbours
IMPLICIT NONE
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::INODEr
INTEGER,INTENT(IN)::N,IMAXE,RESTART,IMAXN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER::I,J,JI,K,LM,IEX,KMAXN,KK,KMAXE,IHAX1,NUMNOD,L,ICN,L2,P,M,Q,IJ,c_n1,c_n2,c_n3,c_n4,d_n1,d_n2,d_n3,d_n4
INTEGER::IT1,IT2,IT3,IT4,IT5,IT6,IT7,IT8,ITX,INX,ITF,IFG,print_i,I1,I2,I3,I4,I5,I6,I7,I8,IMBG,ICFACE
REAL::DELTA,CPUER
   CHARACTER(LEN=20)::PROC,NEIBFILE,PROC3
   KMAXE=XMPIELRANK(N)
	
	
	

	
		inoder2(:)%NUMBEROFNEIB=0
		
	

	
        
	DO I=1,KMAXE
	
	
! 	if (itestcase.eq.4)then
	allocate(ielem(n,i)%vortex(1)); ielem(n,i)%vortex=zero
! 	end if
	
	SELECT CASE (IELEM(N,I)%ISHAPE)

		CASE(4) !prism
	
	allocate(IELEM(N,I)%NODES_FACES(5,4))
	
	allocate(IELEM(N,I)%INEIGHG(5))
	allocate(IELEM(N,I)%TYPES_FACES(5))
	
	IELEM(N,I)%TYPES_FACES(1:2)=6; IELEM(N,I)%TYPES_FACES(3:5)=5

	IELEM(N,I)%NODES_FACES(1:5,1:4)=0
	
	IELEM(N,I)%INEIGHG(1:5)=0
	IELEM(N,I)%ifca=5; 
	IELEM(N,I)%VDEC=3
	
	IELEM(N,I)%TOTVOLUME=0.0;IELEM(N,I)%MINEDGE=0.0;IELEM(N,I)%WALLDIST=0.0
	!FIRST FACE
	IELEM(N,I)%NODES_FACES(1,1)=IELEM(N,I)%NODES(5)
	IELEM(N,I)%NODES_FACES(1,2)=IELEM(N,I)%NODES(6)
	IELEM(N,I)%NODES_FACES(1,3)=IELEM(N,I)%NODES(4)
	
	!SEC FACE
	IELEM(N,I)%NODES_FACES(2,1)=IELEM(N,I)%NODES(2) 
	IELEM(N,I)%NODES_FACES(2,2)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(2,3)=IELEM(N,I)%NODES(3)
	
	!THIRD FACE
	IELEM(N,I)%NODES_FACES(3,1)=IELEM(N,I)%NODES(1)		!125,154
	IELEM(N,I)%NODES_FACES(3,2)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(3,3)=IELEM(N,I)%NODES(5)
	IELEM(N,I)%NODES_FACES(3,4)=IELEM(N,I)%NODES(4)
	!FOURTH FACE
	IELEM(N,I)%NODES_FACES(4,1)=IELEM(N,I)%NODES(2)		!236,265
	IELEM(N,I)%NODES_FACES(4,2)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(4,3)=IELEM(N,I)%NODES(6)
	IELEM(N,I)%NODES_FACES(4,4)=IELEM(N,I)%NODES(5)

	IELEM(N,I)%NODES_FACES(5,1)=IELEM(N,I)%NODES(3)		!314,346	
	IELEM(N,I)%NODES_FACES(5,2)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(5,3)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(5,4)=IELEM(N,I)%NODES(6)
! 	IELEM(N,I)%DEC_FACES(3,1)=IELEM(N,I)%NODES(1) ;IELEM(N,I)%DEC_FACES(3,2)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(3,3)=IELEM(N,I)%NODES(5)
! 	IELEM(N,I)%DEC_FACES(4,1)=IELEM(N,I)%NODES(1) ;IELEM(N,I)%DEC_FACES(4,2)=IELEM(N,I)%NODES(5) ;IELEM(N,I)%DEC_FACES(4,3)=IELEM(N,I)%NODES(4)
! 	IELEM(N,I)%DEC_FACES(5,1)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(5,2)=IELEM(N,I)%NODES(3) ;IELEM(N,I)%DEC_FACES(5,3)=IELEM(N,I)%NODES(6)
! 	IELEM(N,I)%DEC_FACES(6,1)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(6,2)=IELEM(N,I)%NODES(6) ;IELEM(N,I)%DEC_FACES(6,3)=IELEM(N,I)%NODES(5)
! 	IELEM(N,I)%DEC_FACES(7,1)=IELEM(N,I)%NODES(3) ;IELEM(N,I)%DEC_FACES(7,2)=IELEM(N,I)%NODES(1) ;IELEM(N,I)%DEC_FACES(7,3)=IELEM(N,I)%NODES(4)
! 	IELEM(N,I)%DEC_FACES(8,1)=IELEM(N,I)%NODES(3) ;IELEM(N,I)%DEC_FACES(8,2)=IELEM(N,I)%NODES(4) ;IELEM(N,I)%DEC_FACES(8,3)=IELEM(N,I)%NODES(6)
! 
! 	IELEM(N,I)%DEC_FACES(1,1:3)=IELEM(N,I)%NODES_FACES(1,1:3)
! 	IELEM(N,I)%DEC_FACES(2,1:3)=IELEM(N,I)%NODES_FACES(2,1:3)
	allocate(ielem(n,i)%reorient(IELEM(N,I)%ifca));ielem(n,i)%reorient(:)=0
	allocate(ielem(n,i)%faceanglex(IELEM(N,I)%ifca));allocate(ielem(n,i)%faceangley(IELEM(N,I)%ifca));
	if (fastest.eq.0)then
	allocate(ielem(n,i)%indexi(ielem(n,i)%ifca));ielem(n,i)%indexi(:)=0
	end if
	if (rungekutta.ge.2)then
	  allocate(ielem(n,i)%dih(ielem(n,i)%ifca)); ielem(n,i)%dih=zero
! 	  allocate(ielem(n,i)%dih2(ielem(n,i)%ifca,DIMS)); ielem(n,i)%dih2=zero
	  
	end if
	CASE(1) !HEXAHEDRAL
	allocate(IELEM(N,I)%NODES_FACES(6,4))
	allocate(IELEM(N,I)%INEIGHG(6))
	allocate(IELEM(N,I)%TYPEs_FACES(6))
	IELEM(N,I)%TYPEs_FACES(1:6)=5
	IELEM(N,I)%VDEC=6
	IELEM(N,I)%TOTVOLUME=0.0;IELEM(N,I)%MINEDGE=0.0;IELEM(N,I)%WALLDIST=0.0
	IELEM(N,I)%NODES_FACES(:,:)=0
	
	IELEM(N,I)%INEIGHG=0
	IELEM(N,I)%ifca=6; 
	!FIRST FACE
	IELEM(N,I)%NODES_FACES(1,1)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(1,2)=IELEM(N,I)%NODES(1)			!DEC214,243
	IELEM(N,I)%NODES_FACES(1,3)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(1,4)=IELEM(N,I)%NODES(3)
	!SEC FACE
	IELEM(N,I)%NODES_FACES(2,1)=IELEM(N,I)%NODES(5) 
	IELEM(N,I)%NODES_FACES(2,2)=IELEM(N,I)%NODES(6)				!567,578
	IELEM(N,I)%NODES_FACES(2,3)=IELEM(N,I)%NODES(7)
	IELEM(N,I)%NODES_FACES(2,4)=IELEM(N,I)%NODES(8)
	!THIRD FACE
	IELEM(N,I)%NODES_FACES(3,1)=IELEM(N,I)%NODES(6)
	IELEM(N,I)%NODES_FACES(3,2)=IELEM(N,I)%NODES(2)				!623,637
	IELEM(N,I)%NODES_FACES(3,3)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(3,4)=IELEM(N,I)%NODES(7)
	!FOURTH FACE
	IELEM(N,I)%NODES_FACES(4,1)=IELEM(N,I)%NODES(8)					
	IELEM(N,I)%NODES_FACES(4,2)=IELEM(N,I)%NODES(4)			!841,815
	IELEM(N,I)%NODES_FACES(4,3)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(4,4)=IELEM(N,I)%NODES(5)
	!FIFTH FACE
	IELEM(N,I)%NODES_FACES(5,1)=IELEM(N,I)%NODES(7)				!734,748
	IELEM(N,I)%NODES_FACES(5,2)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(5,3)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(5,4)=IELEM(N,I)%NODES(8)
	!SIXTH FACE
	IELEM(N,I)%NODES_FACES(6,1)=IELEM(N,I)%NODES(5)		!512,526
	IELEM(N,I)%NODES_FACES(6,2)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(6,3)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(6,4)=IELEM(N,I)%NODES(6)

! 	IELEM(N,I)%DEC_FACES(1,1)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(1,2)=IELEM(N,I)%NODES(1) ;IELEM(N,I)%DEC_FACES(1,3)=IELEM(N,I)%NODES(4)
! 	IELEM(N,I)%DEC_FACES(2,1)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(2,2)=IELEM(N,I)%NODES(4) ;IELEM(N,I)%DEC_FACES(2,3)=IELEM(N,I)%NODES(3)
! 	IELEM(N,I)%DEC_FACES(3,1)=IELEM(N,I)%NODES(5) ;IELEM(N,I)%DEC_FACES(3,2)=IELEM(N,I)%NODES(6) ;IELEM(N,I)%DEC_FACES(3,3)=IELEM(N,I)%NODES(7)
! 	IELEM(N,I)%DEC_FACES(4,1)=IELEM(N,I)%NODES(5) ;IELEM(N,I)%DEC_FACES(4,2)=IELEM(N,I)%NODES(7) ;IELEM(N,I)%DEC_FACES(4,3)=IELEM(N,I)%NODES(8)
! 	IELEM(N,I)%DEC_FACES(5,1)=IELEM(N,I)%NODES(6) ;IELEM(N,I)%DEC_FACES(5,2)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(5,3)=IELEM(N,I)%NODES(3)
! 	IELEM(N,I)%DEC_FACES(6,1)=IELEM(N,I)%NODES(6) ;IELEM(N,I)%DEC_FACES(6,2)=IELEM(N,I)%NODES(3) ;IELEM(N,I)%DEC_FACES(6,3)=IELEM(N,I)%NODES(7)
! 	IELEM(N,I)%DEC_FACES(7,1)=IELEM(N,I)%NODES(7) ;IELEM(N,I)%DEC_FACES(7,2)=IELEM(N,I)%NODES(3) ;IELEM(N,I)%DEC_FACES(7,3)=IELEM(N,I)%NODES(4)
! 	IELEM(N,I)%DEC_FACES(8,1)=IELEM(N,I)%NODES(7) ;IELEM(N,I)%DEC_FACES(8,2)=IELEM(N,I)%NODES(4) ;IELEM(N,I)%DEC_FACES(8,3)=IELEM(N,I)%NODES(8)
! 	IELEM(N,I)%DEC_FACES(9,1)=IELEM(N,I)%NODES(8) ;IELEM(N,I)%DEC_FACES(9,2)=IELEM(N,I)%NODES(4) ;IELEM(N,I)%DEC_FACES(9,3)=IELEM(N,I)%NODES(1)
! 	IELEM(N,I)%DEC_FACES(10,1)=IELEM(N,I)%NODES(8) ;IELEM(N,I)%DEC_FACES(10,2)=IELEM(N,I)%NODES(1) ;IELEM(N,I)%DEC_FACES(10,3)=IELEM(N,I)%NODES(5)
! 	IELEM(N,I)%DEC_FACES(11,1)=IELEM(N,I)%NODES(5) ;IELEM(N,I)%DEC_FACES(11,2)=IELEM(N,I)%NODES(1) ;IELEM(N,I)%DEC_FACES(11,3)=IELEM(N,I)%NODES(2)
! 	IELEM(N,I)%DEC_FACES(12,1)=IELEM(N,I)%NODES(5) ;IELEM(N,I)%DEC_FACES(12,2)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(12,3)=IELEM(N,I)%NODES(6)
	allocate(ielem(n,i)%reorient(IELEM(N,I)%ifca));ielem(n,i)%reorient(:)=0
	allocate(ielem(n,i)%faceanglex(IELEM(N,I)%ifca));allocate(ielem(n,i)%faceangley(IELEM(N,I)%ifca));
	if (fastest.eq.0)then
	allocate(ielem(n,i)%indexi(ielem(n,i)%ifca));ielem(n,i)%indexi(:)=0
	end if
	if (rungekutta.ge.2)then
	  allocate(ielem(n,i)%dih(ielem(n,i)%ifca)); ielem(n,i)%dih=zero
! 	  allocate(ielem(n,i)%dih2(ielem(n,i)%ifca,DIMS)); ielem(n,i)%dih2=zero
! 	  
	end if


	CASE(2) !TETRAHEDRAL
	allocate(IELEM(N,I)%NODES_FACES(4,3))
	allocate(IELEM(N,I)%INEIGHG(4))
	allocate(IELEM(N,I)%TYPES_FACES(4))
	IELEM(N,I)%TYPES_FACES(1:4)=6
	IELEM(N,I)%NODES_FACES(:,:)=0
	
	IELEM(N,I)%INEIGHG=0
	IELEM(N,I)%VDEC=1
	IELEM(N,I)%TOTVOLUME=0.0;IELEM(N,I)%MINEDGE=0.0;IELEM(N,I)%WALLDIST=0.0
	  IELEM(N,I)%ifca=4; 
	!FIRST FACE
	IELEM(N,I)%NODES_FACES(1,1)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(1,2)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(1,3)=IELEM(N,I)%NODES(1)
	
	!SEC FACE
	IELEM(N,I)%NODES_FACES(2,1)=IELEM(N,I)%NODES(2) 
	IELEM(N,I)%NODES_FACES(2,2)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(2,3)=IELEM(N,I)%NODES(3)
	
	!THIRD FACE
	IELEM(N,I)%NODES_FACES(3,1)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(3,2)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(3,3)=IELEM(N,I)%NODES(4)
	
	!FOURTH FACE
	IELEM(N,I)%NODES_FACES(4,1)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(4,2)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(4,3)=IELEM(N,I)%NODES(4)

	allocate(ielem(n,i)%faceanglex(IELEM(N,I)%ifca));allocate(ielem(n,i)%faceangley(IELEM(N,I)%ifca));
	allocate(ielem(n,i)%reorient(IELEM(N,I)%ifca));ielem(n,i)%reorient(:)=0
	if (fastest.eq.0)then
	allocate(ielem(n,i)%indexi(ielem(n,i)%ifca));ielem(n,i)%indexi(:)=0
	end if
	if (rungekutta.ge.2)then
	  allocate(ielem(n,i)%dih(ielem(n,i)%ifca)); ielem(n,i)%dih=zero
! 	  allocate(ielem(n,i)%dih2(ielem(n,i)%ifca,DIMS)); ielem(n,i)%dih2=zero
	  
	end if
	CASE(3) !pyramidal
	allocate(IELEM(N,I)%NODES_FACES(5,4))
	
	allocate(IELEM(N,I)%INEIGHG(5))
	allocate(IELEM(N,I)%TYPES_FACES(5))
	IELEM(N,I)%TYPES_FACES(1)=5; IELEM(N,I)%TYPES_FACES(2:5)=6
	IELEM(N,I)%NODES_FACES(:,:)=0
	
	IELEM(N,I)%INEIGHG=0
	IELEM(N,I)%VDEC=2
	IELEM(N,I)%TOTVOLUME=0.0;IELEM(N,I)%MINEDGE=0.0;IELEM(N,I)%WALLDIST=0.0
	 IELEM(N,I)%ifca=5; 
	!FIRST FACE
	IELEM(N,I)%NODES_FACES(1,1)=IELEM(N,I)%NODES(4)		!432,421
	IELEM(N,I)%NODES_FACES(1,2)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(1,3)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(1,4)=IELEM(N,I)%NODES(1)
	!SEC FACE
	IELEM(N,I)%NODES_FACES(2,1)=IELEM(N,I)%NODES(1) 
	IELEM(N,I)%NODES_FACES(2,2)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(2,3)=IELEM(N,I)%NODES(5)
	
	!THIRD FACE
	IELEM(N,I)%NODES_FACES(3,1)=IELEM(N,I)%NODES(2)
	IELEM(N,I)%NODES_FACES(3,2)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(3,3)=IELEM(N,I)%NODES(5)
	
	!FOURTH FACE
	IELEM(N,I)%NODES_FACES(4,1)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(4,2)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(4,3)=IELEM(N,I)%NODES(5)

	IELEM(N,I)%NODES_FACES(5,1)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(5,2)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(5,3)=IELEM(N,I)%NODES(5)

! 	IELEM(N,I)%DEC_FACES(1,1)=IELEM(N,I)%NODES(4) ;IELEM(N,I)%DEC_FACES(1,2)=IELEM(N,I)%NODES(3) ;IELEM(N,I)%DEC_FACES(1,3)=IELEM(N,I)%NODES(2)
! 	IELEM(N,I)%DEC_FACES(2,1)=IELEM(N,I)%NODES(4) ;IELEM(N,I)%DEC_FACES(2,2)=IELEM(N,I)%NODES(2) ;IELEM(N,I)%DEC_FACES(2,3)=IELEM(N,I)%NODES(1)
! 
! 	IELEM(N,I)%DEC_FACES(3,1:3)=IELEM(N,I)%NODES_FACES(2,1:3)
! 	IELEM(N,I)%DEC_FACES(4,1:3)=IELEM(N,I)%NODES_FACES(3,1:3)
! 	IELEM(N,I)%DEC_FACES(5,1:3)=IELEM(N,I)%NODES_FACES(4,1:3)
! 	IELEM(N,I)%DEC_FACES(6,1:3)=IELEM(N,I)%NODES_FACES(5,1:3)

	  if (rungekutta.ge.2)then
	  allocate(ielem(n,i)%dih(ielem(n,i)%ifca)); ielem(n,i)%dih=zero
! 	  allocate(ielem(n,i)%dih2(ielem(n,i)%ifca,DIMS)); ielem(n,i)%dih2=zero
	  
	end if
	
	allocate(ielem(n,i)%faceanglex(IELEM(N,I)%ifca));allocate(ielem(n,i)%faceangley(IELEM(N,I)%ifca));
	
	allocate(ielem(n,i)%reorient(IELEM(N,I)%ifca));ielem(n,i)%reorient(:)=0
	if (fastest.eq.0)then
	allocate(ielem(n,i)%indexi(ielem(n,i)%ifca));ielem(n,i)%indexi(:)=0
	end if


	


	CASE(5) !quadrilateral
	IELEM(N,I)%ifca=4; 
	allocate(IELEM(N,I)%NODES_FACES(4,2))
	allocate(IELEM(N,I)%INEIGHG(4))
	IELEM(N,I)%NODES_FACES(:,:)=0
	IELEM(N,I)%VDEC=2
	IELEM(N,I)%TOTVOLUME=0.0;IELEM(N,I)%MINEDGE=0.0;IELEM(N,I)%WALLDIST=0.0
	
	IELEM(N,I)%INEIGHG=0
	!FIRST FACE
	IELEM(N,I)%NODES_FACES(1,1)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(1,2)=IELEM(N,I)%NODES(2)
	
	
	!SEC FACE
	IELEM(N,I)%NODES_FACES(2,1)=IELEM(N,I)%NODES(2) 
	IELEM(N,I)%NODES_FACES(2,2)=IELEM(N,I)%NODES(3)
	
	
	!THIRD FACE
	IELEM(N,I)%NODES_FACES(3,1)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(3,2)=IELEM(N,I)%NODES(4)
	
	
	!FOURTH FACE
	IELEM(N,I)%NODES_FACES(4,1)=IELEM(N,I)%NODES(4)
	IELEM(N,I)%NODES_FACES(4,2)=IELEM(N,I)%NODES(1)
	allocate(ielem(n,i)%faceanglex(IELEM(N,I)%ifca));allocate(ielem(n,i)%faceangley(IELEM(N,I)%ifca));
	allocate(ielem(n,i)%reorient(IELEM(N,I)%ifca));ielem(n,i)%reorient(:)=0
	if (fastest.eq.0)then
	allocate(ielem(n,i)%indexi(ielem(n,i)%ifca));ielem(n,i)%indexi(:)=0
	end if
	if (rungekutta.ge.2)then
	  allocate(ielem(n,i)%dih(ielem(n,i)%ifca)); ielem(n,i)%dih=zero
! 	  allocate(ielem(n,i)%dih2(ielem(n,i)%ifca,DIMS)); ielem(n,i)%dih2=zero
	  
	end if
	CASE(6) !Triangular
	IELEM(N,I)%ifca=3; 
	allocate(IELEM(N,I)%NODES_FACES(3,2))
	allocate(IELEM(N,I)%INEIGHG(3))
	IELEM(N,I)%NODES_FACES(:,:)=0
	  IELEM(N,I)%VDEC=1
	IELEM(N,I)%TOTVOLUME=0.0;IELEM(N,I)%MINEDGE=0.0;IELEM(N,I)%WALLDIST=0.0

	IELEM(N,I)%INEIGHG=0
	!FIRST FACE
	IELEM(N,I)%NODES_FACES(1,1)=IELEM(N,I)%NODES(1)
	IELEM(N,I)%NODES_FACES(1,2)=IELEM(N,I)%NODES(2)
	
	
	!SEC FACE
	IELEM(N,I)%NODES_FACES(2,1)=IELEM(N,I)%NODES(2) 
	IELEM(N,I)%NODES_FACES(2,2)=IELEM(N,I)%NODES(3)
	
	
	!THIRD FACE
	IELEM(N,I)%NODES_FACES(3,1)=IELEM(N,I)%NODES(3)
	IELEM(N,I)%NODES_FACES(3,2)=IELEM(N,I)%NODES(1)
	
	
	
	





	 if (rungekutta.ge.2)then
	  allocate(ielem(n,i)%dih(ielem(n,i)%ifca)); ielem(n,i)%dih=zero
! 	  allocate(ielem(n,i)%dih2(ielem(n,i)%ifca,DIMS)); ielem(n,i)%dih2=zero
	  
	end if



	allocate(ielem(n,i)%faceanglex(IELEM(N,I)%ifca));allocate(ielem(n,i)%faceangley(IELEM(N,I)%ifca));

	allocate(ielem(n,i)%reorient(IELEM(N,I)%ifca));ielem(n,i)%reorient(:)=0
	if (fastest.eq.0)then
	allocate(ielem(n,i)%indexi(ielem(n,i)%ifca));ielem(n,i)%indexi(:)=0
	end if


	END SELECT
	end do

	
! 	
       

	DO I=1,KMAXE
		DO J=1,IELEM(N,I)%nonodes
		

		K=IELEM(N,I)%nodeS(j)
		inoder2(K)%NUMBEROFNEIB=inoder2(K)%NUMBEROFNEIB+1
		
		end do
		
		
	END DO
	
	

	
	

	
	DO I=1,IMAXN
		if (inoder2(I)%NUMBEROFNEIB.gt.0)then
		allocate(inoder2(i)%neibids(inoder2(I)%NUMBEROFNEIB))
		inoder2(i)%neibids(:)=0
		inoder2(I)%NUMBEROFNEIB=0
		end if
	END DO
	
	
	
! 	
	
	DO I=1,KMAXE
		DO Ji=1,IELEM(N,I)%nonodes
		K=IELEM(N,I)%nodeS(ji)
		inoder2(K)%NUMBEROFNEIB=inoder2(K)%NUMBEROFNEIB+1
		J=inoder2(K)%NUMBEROFNEIB
		inoder2(K)%NEIBIDS(J)=I
		end do
		
		
	END DO
	
! 	


	PRINT_I=KMAXE/5
OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',position='append')	
	if (n.eq. 0) then
	
	write(63,*)'Find Neigbours'
	end if
CLOSE(63)

 

	DO IHAX1=1,KMAXE


	
	
	
			NUMNOD=IELEM(N,IHAX1)%nonodes
	 		NODELIST(1:NUMNOD) = IELEM(N,IHAX1)%nodes(1:NUMNOD)
	 		
			
	

		
				
				LIST(:)=0
				M=0
				DO K=1,NUMNOD
					P=NODELIST(K)
						if (inoder2(P)%NUMBEROFNEIB.gt.0)then
						DO Q=1,inoder2(P)%NUMBEROFNEIB
							J=inoder2(P)%NEIBIDS(Q)
							ICN=0
							DO IJ=1,M
								IF (J.EQ.LIST(IJ))THEN
									ICN=1
								END IF
							END DO
							IF ((ICN.EQ.0).AND.(IHAX1.NE.J))THEN
								M=M+1
								LIST(M)=J
							END IF
						END DO
						end if
				END DO
				
			      
    DO L=1,IELEM(N,IHAX1)%IFCA
  if (dimensiona.eq.3)then

  DO P=1,M
  
  J=LIST(P)
	  IF (J.EQ.IHAX1) CYCLE
	DO L2=1,IELEM(N,J)%IFCA
		  
		  IF (IELEM(N,IHAX1)%TYPES_FACES(L).EQ.IELEM(N,J)%TYPES_FACES(L2))THEN
		  if (IELEM(N,IHAX1)%TYPES_FACES(L).eq.6)then
		  c_n4=0
		  c_n1=ielem(n,ihax1)%NODES_FACES(l,1)
		  c_n2=ielem(n,ihax1)%NODES_FACES(l,2)
		  c_n3=ielem(n,ihax1)%NODES_FACES(l,3)
		  d_n4=0
		  d_n1=ielem(n,j)%NODES_FACES(l2,1)
		  d_n2=ielem(n,j)%NODES_FACES(l2,2)
		  d_n3=ielem(n,j)%NODES_FACES(l2,3)
		  if (((c_n1.eq.d_n1).or.(c_n1.eq.d_n2).or.(c_n1.eq.d_n3)).and.&
		      ((c_n2.eq.d_n1).or.(c_n2.eq.d_n2).or.(c_n2.eq.d_n3)).and.&
		      ((c_n3.eq.d_n1).or.(c_n3.eq.d_n2).or.(c_n3.eq.d_n3)))then
			  IELEM(N,IHAX1)%INEIGHG(L)=IELEM(N,J)%IHEXGL
			 
			  
		  GOTO 101
		  END IF
		  else
		  c_n4=ielem(n,ihax1)%NODES_FACES(l,4)
		  c_n1=ielem(n,ihax1)%NODES_FACES(l,1)
		  c_n2=ielem(n,ihax1)%NODES_FACES(l,2)
		  c_n3=ielem(n,ihax1)%NODES_FACES(l,3)
		  d_n4=ielem(n,j)%NODES_FACES(l2,4)
		  d_n1=ielem(n,j)%NODES_FACES(l2,1)
		  d_n2=ielem(n,j)%NODES_FACES(l2,2)
		  d_n3=ielem(n,j)%NODES_FACES(l2,3)
		  if (((c_n1.eq.d_n1).or.(c_n1.eq.d_n2).or.(c_n1.eq.d_n3).or.(c_n1.eq.d_n4)).and.&
		      ((c_n2.eq.d_n1).or.(c_n2.eq.d_n2).or.(c_n2.eq.d_n3).or.(c_n2.eq.d_n4)).and.&
		      ((c_n3.eq.d_n1).or.(c_n3.eq.d_n2).or.(c_n3.eq.d_n3).or.(c_n3.eq.d_n4)).and.&
		      ((c_n4.eq.d_n1).or.(c_n4.eq.d_n2).or.(c_n4.eq.d_n3).or.(c_n4.eq.d_n4)))then
			  IELEM(N,IHAX1)%INEIGHG(L)=IELEM(N,J)%IHEXGL
			  
			  
		  GOTO 101
		  END IF
		  end if

		  
		 
		 
		  
		  
		  END IF
	END DO
  END DO
  else
  DO P=1,M
  J=LIST(P)
	  IF (J.EQ.IHAX1) CYCLE
	DO L2=1,IELEM(N,J)%IFCA
		  c_n1=ielem(n,ihax1)%NODES_FACES(l,1)
		  c_n2=ielem(n,ihax1)%NODES_FACES(l,2)
		  

		  d_n1=ielem(n,j)%NODES_FACES(l2,1)
		  d_n2=ielem(n,j)%NODES_FACES(l2,2)
		 
		 
		  if (((c_n1.eq.d_n1).or.(c_n1.eq.d_n2)).and.&
		      ((c_n2.eq.d_n1).or.(c_n2.eq.d_n2)))then
			  IELEM(N,IHAX1)%INEIGHG(L)=IELEM(N,J)%IHEXGL
			  
		  GOTO 101
		  END IF
	END DO
  END DO
  end if





  
  101 CONTINUE
			END DO
		END DO	
     
	
		
		
		

	JI=0
      do i=1,kmaxe
	    ielem(n,i)%interior=0
	    DO L=1,IELEM(N,I)%IFCA
	      if (IELEM(N,I)%INEIGHG(l).eq.0)then
	      ielem(n,i)%interior=1
	      JI=JI+1
	       end if
	    end do
      end do
     
      
		
END SUBROUTINE NEIGHBOURSS
! !---------------------------------------------------------------------------------------------!







SUBROUTINE CONS(N,ICONR,ICONS,IPERIODICITY,XMPIELRANK,ISIZE,ICONRPA,ICONRPM,ICONSPO,XPER,&
YPER,ZPER,ICONRPF,NUMNEIGHBOURS,TYPESTEN)
!> @brief
!> This subroutine establishes the connectivity across different cpus
IMPLICIT NONE
INTEGER,INTENT(IN)::IPERIODICITY,N,ISIZE,NUMNEIGHBOURS,TYPESTEN
REAL::SMALL,DIST
REAL,INTENT(IN)::XPER,YPER,ZPER
INTEGER::I,K,JJJ,KJ,J,L,IM,IO,IR,KMAXE,ICPUID,KKJ,countxsize
INTEGER,DIMENSION(6)::JX
REAL,DIMENSION(3)::DUMFACE1,DUMFACE2
INTEGER::FACEPER1,FACEPER2,p,kkj2,kkj3,kkj4,kkj5,JJ1,P1,P2,P3,P4,Q1,Q2,Q3,Q4,J1,J2,J3,J4

INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
TYPE(CONNX),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ICONR,ICONRPA,ICONRPM,ICONRPF
TYPE(CONNX),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ICONS,ICONSPO

KMAXE=XMPIELRANK(N)

SMALL=TOLSMALL					! A SMALL NUMBER!
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

           
KJ=0; KKJ=0
	DO I=1,KMAXE
	IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		DO K=1,IELEM(N,I)%IFCA
		
		  IF (((IELEM(N,I)%INEIGHG(K).EQ.0)))THEN
		      if(IELEM(N,I)%IBOUNDS(K).EQ.0)then
						KJ=KJ+1
		      else
		      if (((IBOUND(N,IELEM(N,I)%IBOUNDS(k))%ICODE)).eq.5)then
						  KJ=KJ+1	
						
		      end if
		      end if
		END IF
		
		

		END DO
	END IF
	END DO
	
	




! ALLOCATE(ICONS(1:COUNTXSIZE))   1ST


ALLOCATE(ICONR(N:N))
ALLOCATE(ICONR(N)%HOWMANYI(1))


  ICONR(N)%PROCID=N
  ICONR(N)%HOWMANYI(1)=KJ			!THIS IS THE TOTAL NUMBER OF NEIGHBOURS I NEED

	
	
	
		ALLOCATE(ICONR(N)%WHICHI(ICONR(N)%HOWMANYI(1),6))
		ICONR(N)%WHICHI(:,:)=0
! 		
		KJ=0
		DO I=1,KMAXE
		IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		IO=ielem(n,i)%ifca
		DO K=1,IO
		 IF (((IELEM(N,I)%INEIGHG(K).EQ.0)))then
		    if ((IELEM(N,I)%IBOUNDS(K).EQ.0))then
						KJ=KJ+1
				ICONR(N)%WHICHI(KJ,1)=I
				ICONR(N)%WHICHI(KJ,2)=K
		    
						
		      eLSE
		      
		      if (((IBOUND(N,IELEM(N,I)%IBOUNDS(k))%ICODE)).eq.5)then
						  KJ=KJ+1
				ICONR(N)%WHICHI(KJ,1)=I
				ICONR(N)%WHICHI(KJ,2)=K
		      end if
		      end if
		END IF
		END DO
		end if
		END DO

		 




		DO KJ=1,ICONR(N)%HOWMANYI(1)
			      
			      SELECT case(IELEM(N,ICONR(N)%WHICHI(KJ,1))%TYPES_FACES(ICONR(N)%WHICHI(KJ,2)))
			      CASE(5)
			      ICONR(N)%WHICHI(KJ,3:6)=IELEM(N,ICONR(N)%WHICHI(KJ,1))%NODES_FACES(ICONR(N)%WHICHI(KJ,2),1:4)


			      CASE(6)
			      ICONR(N)%WHICHI(KJ,3:5)=IELEM(N,ICONR(N)%WHICHI(KJ,1))%NODES_FACES(ICONR(N)%WHICHI(KJ,2),1:3)
			       ICONR(N)%WHICHI(KJ,6)=0
				end select

			
		END DO





DO IHAX1=1,ICONR(N)%HOWMANYI(1)

		      SELECT case(IELEM(N,ICONR(N)%WHICHI(ihax1,1))%TYPES_FACES(ICONR(N)%WHICHI(ihax1,2)))
		      
			      CASE(5)
			      ICONR(N)%WHICHI(ihax1,3:6)=IELEM(N,ICONR(N)%WHICHI(ihax1,1))%NODES_FACES(ICONR(N)%WHICHI(ihax1,2),1:4)
	 		       NODELIST(1:4) = ICONR(N)%WHICHI(ihax1,3:6)
	 		       
			
			P1=NODELIST(1)
			    DO Q1=1,inoder2(P1)%xne(1)
				  J1=inoder2(P1)%xneib(Q1)
				      P2=NODELIST(2)
				      DO Q2=1,inoder2(P2)%xne(1)
					  J2=inoder2(P2)%xneib(Q2)
					    IF ((J2.EQ.J1).AND.(J2.NE.IELEM(N,ICONR(N)%WHICHI(ihax1,1))%IHEXGL))THEN
					      P3=NODELIST(3)
						  DO Q3=1,inoder2(P3)%xne(1)
						      J3=inoder2(P3)%xneib(Q3)
							  IF ((J3.EQ.J2).AND.(J3.NE.IELEM(N,ICONR(N)%WHICHI(ihax1,1))%IHEXGL))THEN
							      P4=NODELIST(4)
								  DO Q4=1,inoder2(P4)%xne(1)
								    J4=inoder2(P4)%xneib(Q4)
									IF ((J4.EQ.J3).AND.(J4.NE.IELEM(N,ICONR(N)%WHICHI(ihax1,1))%IHEXGL))THEN
									IELEM(N,ICONR(N)%WHICHI(IHAX1,1))%INEIGHG(ICONR(N)%WHICHI(IHAX1,2))=J4
									GO TO 331
									END IF
								  END DO
							  END IF
						  END DO
					    END IF
				     END DO
			    END DO
		CASE(6)
			      ICONR(N)%WHICHI(ihax1,3:5)=IELEM(N,ICONR(N)%WHICHI(ihax1,1))%NODES_FACES(ICONR(N)%WHICHI(ihax1,2),1:3)
	 		       NODELIST(1:3) = ICONR(N)%WHICHI(ihax1,3:5)
	 		       
			
			P1=NODELIST(1)
			    DO Q1=1,inoder2(P1)%xne(1)
				  J1=inoder2(P1)%xneib(Q1)
				      P2=NODELIST(2)
				      DO Q2=1,inoder2(P2)%xne(1)
					  J2=inoder2(P2)%xneib(Q2)
					    IF ((J2.EQ.J1).AND.(J2.NE.IELEM(N,ICONR(N)%WHICHI(ihax1,1))%IHEXGL))THEN
					      P3=NODELIST(3)
						  DO Q3=1,inoder2(P3)%xne(1)
						      J3=inoder2(P3)%xneib(Q3)
							  IF ((J3.EQ.J2).AND.(J3.NE.IELEM(N,ICONR(N)%WHICHI(ihax1,1))%IHEXGL))THEN
									IELEM(N,ICONR(N)%WHICHI(IHAX1,1))%INEIGHG(ICONR(N)%WHICHI(IHAX1,2))=J3
									GO TO 331
							  END IF
						  END DO
					    END IF
				     END DO
			    END DO	
	      END SELECT
	      331 CONTINUE
				
END DO


            KJ=0

	      DO I=1,ICONR(N)%HOWMANYI(1)
! 		    

		       IF ((IELEM(N,ICONR(N)%WHICHI(I,1))%INEIGHG(ICONR(N)%WHICHI(I,2)).EQ.0))THEN
		      IF (IPERIODICITY.EQ.1)THEN
		      if (((IBOUND(N,IELEM(N,ICONR(N)%WHICHI(I,1))%IBOUNDS(ICONR(N)%WHICHI(I,2)))%ICODE)).eq.5)then
						  KJ=KJ+1
			
		      end if
		      else

					    kj=kj+1
		      end if
		      end if
		END DO

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		IF (IPERIODICITY.EQ.1)THEN
		    ! 		
				    
				    ALLOCATE(ICONRPM(N:N))
				    ALLOCATE(ICONRPM(N)%HOWMANYI(1))
				    ICONRPM(N)%HOWMANYI(1)=KJ
				    ALLOCATE(ICONRPM(N)%WHICHI(KJ,3))
				    ICONRPM(N)%WHICHI(:,:)=0
				    KJ=0
		    ! 		

				    


				    KJ=0
				    DO I=1,ICONR(N)%HOWMANYI(1)


					  IF ((IELEM(N,ICONR(N)%WHICHI(I,1))%INEIGHG(ICONR(N)%WHICHI(I,2)).EQ.0))THEN
					  if (((IBOUND(N,IELEM(N,ICONR(N)%WHICHI(I,1))%IBOUNDS(ICONR(N)%WHICHI(I,2)))%ICODE)).eq.5)then
								      KJ=KJ+1
					    ICONRPM(N)%WHICHI(KJ,1)=ICONR(N)%WHICHI(I,1)
					    ICONRPM(N)%WHICHI(KJ,3)=IELEM(N,ICONR(N)%WHICHI(I,1))%IHEXGL
					    ICONRPM(N)%WHICHI(KJ,2)=(ICONR(N)%WHICHI(I,2))
					  end if
					  end if
				    END DO
		    ! 		
				    
				    JJJ=0
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
				    ALLOCATE(ICONRPA(N:N))
				    ALLOCATE(ICONRPF(1:ISIZE-1))
		    
				    DO I=1,ISIZE-1
				    ALLOCATE(ICONRPF(I)%HOWMANYI(1))
				    ICONRPF(I)%HOWMANYI(1)=0
				    ICONRPF(I)%HOWMANYI(1)=ICONRPM(N)%HOWMANYI(1)-JJJ
				    ALLOCATE(ICONRPF(I)%WHICHI(ICONRPF(I)%HOWMANYI(1),1))
				    ICONRPF(I)%WHICHI(:,:)=0
				    END DO
				    ALLOCATE(ICONRPA(N)%HOWMANYI(1))
				    ICONRPA(N)%HOWMANYI(1)=0
				    ICONRPA(N)%HOWMANYI(1)=ICONRPM(N)%HOWMANYI(1)-JJJ
				    ALLOCATE(ICONRPA(N)%WHICHI(ICONRPM(N)%HOWMANYI(1)-JJJ,3))
				    ALLOCATE(ICONRPA(N)%FACX(ICONRPM(N)%HOWMANYI(1)-JJJ,3))
				    ICONRPA(N)%WHICHI(:,:)=0
				    ICONRPA(N)%FACX(:,:)=0.d0
				    KJ=0
				    DO I=1,ICONRPM(N)%HOWMANYI(1)
		   
					    IF ((IELEM(N,ICONRpm(N)%WHICHI(I,1))%INEIGHG(ICONRpm(N)%WHICHI(I,2)).EQ.0))THEN
					  if (((IBOUND(N,IELEM(N,ICONRpm(N)%WHICHI(I,1))%IBOUNDS(ICONRpm(N)%WHICHI(I,2)))%ICODE)).eq.5)then
		    ! 						  KJ=KJ+1
		    ! 			
					    KJ=KJ+1
					    ICONRPA(N)%WHICHI(KJ,1)=ICONRPM(N)%WHICHI(I,1)
					    ICONRPA(N)%WHICHI(KJ,2)=ICONRPM(N)%WHICHI(I,2)
					    ICONRPA(N)%WHICHI(KJ,3)=ICONRPM(N)%WHICHI(I,3)
					      IF (IELEM(N,ICONRPA(N)%WHICHI(KJ,1))%TYPEs_FACES(ICONRPa(N)%WHICHI(kj,2)).EQ.5)THEN

					      IXXFF=4
					      ELSE
					      IXXFF=3
					      END IF
						Iconsi=ICONRPA(N)%WHICHI(KJ,1); facex=ICONRPa(N)%WHICHI(kj,2)

					      CALL COMPUTE_CENTRE3DF(N,Iconsi,facex,IXXFF)
					    

					    ICONRPA(N)%FACX(KJ,1)=CORDS(1)
					    ICONRPA(N)%FACX(KJ,2)=CORDS(2)
					    ICONRPA(N)%FACX(KJ,3)=CORDS(3)
					    end if
					    end if
		    ! 			END IF
				    END DO
				    
				    DUMFACE1(:)=0.0
				    DUMFACE2(:)=0.0
				    FACEPER1=0
				    FACEPER2=0
				    ICPUID=N
				    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		   
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
				    ALLOCATE(ICONSPO(1:ISIZE-1))
		    
				    DO I=1,ISIZE-1
				    ALLOCATE(ICONSPO(I)%HOWMANYTHEY(1))
				    ICONSPO(I)%HOWMANYTHEY(1)=0
				    END DO
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
		    
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
				    K=0
				    DO I=0,ISIZE-1
				    IF (I.NE.N)THEN
				    K=K+1
				    ICONSPO(K)%PROCID=I
				    END IF
				    END DO
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
		    
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
				    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
		    
		    !-------------------FOR DEBUGGING ONLY -----------------------------------------!
					    ICPUID=N
				    K=0
				    DO I=0,ISIZE-1
					    IF (I.NE.N)THEN
							    K=K+1
			CALL MPI_SENDRECV(ICONRPA(N)%HOWMANYI(1),1,MPI_INTEGER,I,ICPUID,&
			ICONSPO(K)%HOWMANYTHEY(1),1,MPI_INTEGER,ICONSPO(K)%PROCID,ICONSPO(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
					    END IF
				    END DO


				    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
				    DO I=1,ISIZE-1
				    ALLOCATE(ICONSPO(I)%WHICHTHEY(ICONSPO(I)%HOWMANYTHEY(1),1))
				    ALLOCATE(ICONSPO(I)%FACX(ICONSPO(I)%HOWMANYTHEY(1),3))
				    ICONSPO(I)%WHICHTHEY(:,:)=0
				    ICONSPO(I)%FACX(:,:)=0.d0
				    END DO
				    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
					    ICPUID=N
				    K=0
		    DO I=0,ISIZE-1
		    IF (I.NE.N)THEN
		      K=K+1
		      CALL MPI_SENDRECV(ICONRPA(N)%FACX(1:ICONRPA(N)%HOWMANYI(1),1:3),&
		    ICONRPA(N)%HOWMANYI(1)*3,MPI_DOUBLE_PRECISION,I,ICPUID,&
		      ICONSPO(K)%FACX(1:ICONSPO(K)%HOWMANYTHEY(1),1:3),ICONSPO(K)%HOWMANYTHEY(1)*3,&
		    MPI_DOUBLE_PRECISION,ICONSPO(K)%PROCID,ICONSPO(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
		    END IF
		    END DO

			      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
			      K=0
			      JJ1=0
			      DO J=0,ISIZE-1
				      IF (J.NE.N)THEN
				      K=K+1
				      DO I=1,ICONSPO(K)%HOWMANYTHEY(1)
					VEXT(1,1)=ICONSPO(K)%FACX(I,1); VEXT(1,2)=ICONSPO(K)%FACX(I,2); VEXT(1,3)=ICONSPO(K)%FACX(I,3)


				      DO KJ=1,ICONRPA(N)%HOWMANYI(1)
					VEXT(2,1)=ICONRPA(N)%FACX(KJ,1); VEXT(2,2)=ICONRPA(N)%FACX(KJ,2);VEXT(2,3)=ICONRPA(N)%FACX(KJ,3)

		dist=distance3(n)


						    if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
						    ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
						    if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
						    JJ1=JJ1+1
	      ! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL).or.((abs(dist-zper)).lt.TOLSMALL))then

						    ICONSPO(K)%WHICHTHEY(I,1)=ICONRPA(N)%WHICHI(KJ,3);go to 401
						    end if
						    end if
						    if (((abs(vext(2,2)-Yper).lt.tolsmall).or.(abs((abs(vext(2,2)-Yper))-Yper).lt.tolsmall)).and.&
						    ((abs(vext(1,2)-Yper).lt.tolsmall).or.(abs((abs(vext(1,2)-Yper))-Yper).lt.tolsmall)))then
						    if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
						    JJ1=JJ1+1
	      ! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL).or.((abs(dist-zper)).lt.TOLSMALL))then

						    ICONSPO(K)%WHICHTHEY(I,1)=ICONRPA(N)%WHICHI(KJ,3);go to 401
						    end if
						    end if
						    
						    
						    if (((abs(vext(2,3)-Zper).lt.tolsmall).or.(abs((abs(vext(2,3)-Zper))-Zper).lt.tolsmall)).and.&
						    ((abs(vext(1,3)-Zper).lt.tolsmall).or.(abs((abs(vext(1,3)-Zper))-Zper).lt.tolsmall)))then
						    if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,1)-vext(1,1)).lt.tolsmall))then
						    JJ1=JJ1+1
	      ! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL).or.((abs(dist-zper)).lt.TOLSMALL))then

						  ICONSPO(K)%WHICHTHEY(I,1)=ICONRPA(N)%WHICHI(KJ,3);go to 401
						  
						    end if
						    end if










					  END DO
					  401 CONTINUE
					  END DO
				      END IF
			      END DO
! 			
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		
		      !-------------------FOR DEBUGGING ONLY -----------------------------------------!
			K=0
		!-------------------FOR DEBUGGING ONLY -----------------------------------------!
! 		
		!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
			ICPUID=N
		K=0
DO I=0,ISIZE-1
IF (I.NE.N)THEN
  K=K+1
  CALL MPI_SENDRECV(ICONSPO(K)%WHICHTHEY(1:ICONSPO(K)%HOWMANYTHEY(1),1:1),&
ICONSPO(K)%HOWMANYTHEY(1),MPI_INTEGER,I,ICPUID,&
  ICONRPF(K)%WHICHI(1:ICONRPF(K)%HOWMANYI(1),1:1),ICONRPF(K)%HOWMANYI(1),&
MPI_INTEGER,ICONSPO(K)%PROCID,ICONSPO(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
END IF
END DO

K=0
! 		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		DO I=0,ISIZE-1
			IF (I.NE.N)THEN
			K=K+1
			DO KJ=1,ICONRPF(K)%HOWMANYI(1)
			IF (ICONRPF(K)%WHICHI(KJ,1).GT.0)THEN
! 			
 			IELEM(N,ICONRPA(N)%WHICHI(KJ,1))%Ineighg(ICONRPA(N)%WHICHI(KJ,2))=ICONRPF(K)%WHICHI(KJ,1)
			END IF
			END DO
			END IF
		END DO
		KJ=0
		DO I=1,ICONRPA(N)%HOWMANYI(1)
			
			IF (IELEM(N,ICONRPA(N)%WHICHI(I,1))%Ineighg(ICONRPA(N)%WHICHI(I,2)).EQ.0)THEN
			  if (((IBOUND(N,IELEM(N,ICONRpa(N)%WHICHI(I,1))%IBOUNDS(ICONRpa(N)%WHICHI(I,2)))%ICODE)).eq.5)then
			KJ=KJ+1
				  
			      
			END IF
			end if
		END DO

		if (kj.gt.0)then
		  PRINT*,"error detected in boundary conditions matching rules",kj,ICONRPA(N)%HOWMANYI(1)
		end if
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

			DEALLOCATE(ICONRPA)
			DEALLOCATE(ICONRPF)
			DEALLOCATE(ICONRPM)
			DEALLOCATE(ICONSPO)		
		END IF
		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		
! 		DEALLOCATE(ICONS)
		DEALLOCATE(ICONR)
		
! 		

	
			

END SUBROUTINE CONS



SUBROUTINE CONS2d(N,ICONR,ICONS,IPERIODICITY,XMPIELRANK,ISIZE,ICONRPA,ICONRPM,ICONSPO,XPER,&
YPER,ZPER,ICONRPF,NUMNEIGHBOURS,TYPESTEN)
!> @brief
!> This subroutine establishes the connectivity across different cpus in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::IPERIODICITY,N,ISIZE,NUMNEIGHBOURS,TYPESTEN
REAL::SMALL,DIST
REAL,INTENT(IN)::XPER,YPER,ZPER
INTEGER::I,K,JJJ,KJ,J,L,IM,IO,IR,KMAXE,ICPUID,KKJ
INTEGER,DIMENSION(6)::JX
REAL,DIMENSION(3)::DUMFACE1,DUMFACE2
INTEGER::FACEPER1,FACEPER2,p,kkj2,kkj3,kkj4,kkj5

INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
TYPE(CONNX),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ICONR,ICONRPA,ICONRPM,ICONRPF
TYPE(CONNX),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ICONS,ICONSPO

KMAXE=XMPIELRANK(N)

SMALL=TOLSMALL


           
KJ=0; KKJ=0
	DO I=1,KMAXE
	IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		DO K=1,IELEM(N,I)%IFCA
		
		  IF (((IELEM(N,I)%INEIGHG(K).EQ.0)))THEN
		      if(IELEM(N,I)%IBOUNDS(K).EQ.0)then
						KJ=KJ+1
		      else
		      if (((IBOUND(N,IELEM(N,I)%IBOUNDS(k))%ICODE)).eq.5)then
						  KJ=KJ+1	
						
		      end if
		      end if
		END IF
		
		

		END DO
	END IF
	END DO


ALLOCATE(ICONR(N:N))
ALLOCATE(ICONS(1:ISIZE-1))
ALLOCATE(ICONR(N)%HOWMANYI(1))


DO I=1,ISIZE-1
ALLOCATE(ICONS(I)%HOWMANYTHEY(1))
ICONS(I)%HOWMANYTHEY(1)=0
END DO
ICONR(N)%PROCID=N
ICONR(N)%HOWMANYI(1)=KJ
K=0
DO I=0,ISIZE-1
	 IF (I.NE.N)THEN
	 K=K+1
	 ICONS(K)%PROCID=I
	 END IF
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
ICPUID=N
K=0

DO I=0,ISIZE-1
	IF (I.NE.N)THEN
			K=K+1
			CALL MPI_SENDRECV(ICONR(N)%HOWMANYI(1),1,MPI_INTEGER,I,ICPUID,&
			ICONS(K)%HOWMANYTHEY(1),1,MPI_INTEGER,ICONS(K)%PROCID,ICONS(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
	END IF
END DO

		ALLOCATE(ICONR(N)%WHICHI(ICONR(N)%HOWMANYI(1),4))
		ICONR(N)%WHICHI(:,:)=0
		DO I=1,ISIZE-1
		ALLOCATE(ICONS(I)%WHICHTHEY(ICONS(I)%HOWMANYTHEY(1),4))
		ALLOCATE(ICONS(I)%RET(ICONS(I)%HOWMANYTHEY(1)))
		ALLOCATE(ICONS(I)%RETM(ICONR(N)%HOWMANYI(1)))
		ICONS(I)%WHICHTHEY(:,:)=0
		ICONS(I)%RET(:)=0
		ICONS(I)%RETM(:)=0
		END DO		
! 	
		KJ=0
		DO I=1,KMAXE
		IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		IO=ielem(n,i)%ifca
		DO K=1,IO
		 IF (((IELEM(N,I)%INEIGHG(K).EQ.0)))then
		    if ((IELEM(N,I)%IBOUNDS(K).EQ.0))then
						KJ=KJ+1
				ICONR(N)%WHICHI(KJ,1)=I
				ICONR(N)%WHICHI(KJ,2)=K
		    
						
		      eLSE
		      
		      if (((IBOUND(N,IELEM(N,I)%IBOUNDS(k))%ICODE)).eq.5)then
						  KJ=KJ+1
				ICONR(N)%WHICHI(KJ,1)=I
				ICONR(N)%WHICHI(KJ,2)=K
		      end if
		      end if
		END IF
		END DO
		end if
		END DO

		 



		DO KJ=1,ICONR(N)%HOWMANYI(1)
			      
! 			     
			      ICONR(N)%WHICHI(KJ,3:4)=IELEM(N,ICONR(N)%WHICHI(KJ,1))%NODES_FACES(ICONR(N)%WHICHI(KJ,2),1:2)


			
		END DO
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		ICPUID=N
		K=0
		
		
DO I=0,ISIZE-1
IF (I.NE.N)THEN
    K=K+1
    CALL MPI_SENDRECV(ICONR(N)%WHICHI(1:ICONR(N)%HOWMANYI(1),1:4),ICONR(N)%HOWMANYI(1)*4,MPI_INTEGER,I,ICPUID,&
    ICONS(K)%WHICHTHEY(1:ICONS(K)%HOWMANYTHEY(1),1:4),ICONS(K)%HOWMANYTHEY(1)*4,&
MPI_INTEGER,ICONS(K)%PROCID,ICONS(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
END IF
END DO





DO KJ=1,ICONR(N)%HOWMANYI(1)
DO I=1,ISIZE-1
      DO K=1,ICONS(I)%HOWMANYTHEY(1)
 	      if (inoder(icons(i)%whichthey(k,3))%itor.gt.0)then
	      
	      IF ((((ICONR(N)%WHICHI(KJ,3)).EQ.(ICONS(I)%WHICHTHEY(K,3))).OR.((ICONR(N)%WHICHI(KJ,3))&
.EQ.(ICONS(I)%WHICHTHEY(K,4)))).AND.(((ICONR(N)%WHICHI(KJ,4)).EQ.(ICONS(I)%WHICHTHEY(K,3))).OR.((ICONR(N)%WHICHI(KJ,4))&
.EQ.(ICONS(I)%WHICHTHEY(K,4)))))THEN
	      ICONS(I)%RET(K)=IELEM(N,ICONR(N)%WHICHI(KJ,1))%IHEXGL
	      
	      GO TO 801
	      END IF
	      end if
      END DO
END DO
801 CONTINUE
END DO

		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

		K=0
		ICPUID=N
		K=0
DO I=0,ISIZE-1
IF (I.NE.N)THEN
  K=K+1
  CALL MPI_SENDRECV(ICONS(K)%RET(1:ICONS(k)%HOWMANYTHEY(1)),ICONS(K)%HOWMANYTHEY(1),MPI_INTEGER,ICONS(K)%PROCID,ICPUID,&
  ICONS(K)%RETM(1:ICONR(N)%HOWMANYI(1)),ICONR(N)%HOWMANYI(1),MPI_INTEGER,I,I,MPI_COMM_WORLD,STATUS,IERROR)
END IF
END DO

			
			DO I=1,ICONR(N)%HOWMANYI(1)
				K=0
				DO J=0,ISIZE-1
					IF (J.NE.N)THEN
					K=K+1
					IF (ICONS(K)%RETM(I).GT.0)THEN
					IELEM(N,ICONR(N)%WHICHI(I,1))%INEIGHG(ICONR(N)%WHICHI(I,2))=ICONS(K)%RETM(I)
					END IF
					END IF
				END DO
			END DO
		  
		  KJ=0
		DO I=1,ICONR(N)%HOWMANYI(1)


		       IF ((IELEM(N,ICONR(N)%WHICHI(I,1))%INEIGHG(ICONR(N)%WHICHI(I,2)).EQ.0))THEN
		      IF (IPERIODICITY.EQ.1)THEN
		      
! 		     
		      if (IELEM(N,ICONR(N)%WHICHI(I,1))%IBOUNDS(ICONR(N)%WHICHI(I,2)).ne.0)then
		      if (((IBOUND(N,IELEM(N,ICONR(N)%WHICHI(I,1))%IBOUNDS(ICONR(N)%WHICHI(I,2)))%ICODE)).eq.5)then
						  KJ=KJ+1
! 			
		      end if
		      end if
		      else

					    kj=kj+1
		      end if
		      end if
		END DO

	     
		
!-------------------FOR DEBUGGING ONLY -----------------------------------------!



!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		IF (IPERIODICITY.EQ.1)THEN
! 		
		
		ALLOCATE(ICONRPM(N:N))
		ALLOCATE(ICONRPM(N)%HOWMANYI(1))
		ICONRPM(N)%HOWMANYI(1)=KJ
		ALLOCATE(ICONRPM(N)%WHICHI(KJ,3))
		ICONRPM(N)%WHICHI(:,:)=0
		KJ=0
! 		
		


		 KJ=0
		DO I=1,ICONR(N)%HOWMANYI(1)


		       IF ((IELEM(N,ICONR(N)%WHICHI(I,1))%INEIGHG(ICONR(N)%WHICHI(I,2)).EQ.0))THEN
		      if (((IBOUND(N,IELEM(N,ICONR(N)%WHICHI(I,1))%IBOUNDS(ICONR(N)%WHICHI(I,2)))%ICODE)).eq.5)then
						  KJ=KJ+1
			ICONRPM(N)%WHICHI(KJ,1)=ICONR(N)%WHICHI(I,1)
			ICONRPM(N)%WHICHI(KJ,3)=IELEM(N,ICONR(N)%WHICHI(I,1))%IHEXGL
			ICONRPM(N)%WHICHI(KJ,2)=(ICONR(N)%WHICHI(I,2))
		      end if
		      end if
		END DO
! 		
		
		JJJ=0
!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		ALLOCATE(ICONRPA(N:N))
		ALLOCATE(ICONRPF(1:ISIZE-1))
! 		 
! 		ICONRPA(N:N)=0
! 		ICONRPF(1:ISIZE-1)=0
		DO I=1,ISIZE-1
		ALLOCATE(ICONRPF(I)%HOWMANYI(1))
		ICONRPF(I)%HOWMANYI(1)=0
		ICONRPF(I)%HOWMANYI(1)=ICONRPM(N)%HOWMANYI(1)-JJJ
		ALLOCATE(ICONRPF(I)%WHICHI(ICONRPF(I)%HOWMANYI(1),1))
		ICONRPF(I)%WHICHI(:,:)=0
		END DO
		ALLOCATE(ICONRPA(N)%HOWMANYI(1))
		ICONRPA(N)%HOWMANYI(1)=0
		ICONRPA(N)%HOWMANYI(1)=ICONRPM(N)%HOWMANYI(1)-JJJ
		ALLOCATE(ICONRPA(N)%WHICHI(ICONRPM(N)%HOWMANYI(1)-JJJ,3))
		ALLOCATE(ICONRPA(N)%FACX(ICONRPM(N)%HOWMANYI(1)-JJJ,2))
		ICONRPA(N)%WHICHI(:,:)=0
		ICONRPA(N)%FACX(:,:)=0.d0
		KJ=0
		DO I=1,ICONRPM(N)%HOWMANYI(1)
! 			
			 IF ((IELEM(N,ICONRpm(N)%WHICHI(I,1))%INEIGHG(ICONRpm(N)%WHICHI(I,2)).EQ.0))THEN
		      if (((IBOUND(N,IELEM(N,ICONRpm(N)%WHICHI(I,1))%IBOUNDS(ICONRpm(N)%WHICHI(I,2)))%ICODE)).eq.5)then
! 						  KJ=KJ+1
! 			
			KJ=KJ+1
			ICONRPA(N)%WHICHI(KJ,1)=ICONRPM(N)%WHICHI(I,1)
			ICONRPA(N)%WHICHI(KJ,2)=ICONRPM(N)%WHICHI(I,2)
			ICONRPA(N)%WHICHI(KJ,3)=ICONRPM(N)%WHICHI(I,3)
! 			  IF (IELEM(N,ICONRPA(N)%WHICHI(KJ,1))%TYPEs_FACES(ICONRPa(N)%WHICHI(kj,2)).EQ.5)THEN

			  IXXFF=2
! 			  ELSE
! 			  IXXFF=3
! 			  END IF
			    Iconsi=ICONRPA(N)%WHICHI(KJ,1); facex=ICONRPa(N)%WHICHI(kj,2)

			  CALL COMPUTE_CENTRE2dF(N,Iconsi,facex,IXXFF)
			

			ICONRPA(N)%FACX(KJ,1)=CORDS(1)
			ICONRPA(N)%FACX(KJ,2)=CORDS(2)
			
			end if
			end if
! 			END IF
		END DO
		
		DUMFACE1(:)=0.0
		DUMFACE2(:)=0.0
		FACEPER1=0
		FACEPER2=0
		ICPUID=N
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 		
!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		ALLOCATE(ICONSPO(1:ISIZE-1))
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		DO I=1,ISIZE-1
		ALLOCATE(ICONSPO(I)%HOWMANYTHEY(1))
		ICONSPO(I)%HOWMANYTHEY(1)=0
		END DO
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		K=0
		DO I=0,ISIZE-1
		IF (I.NE.N)THEN
		K=K+1
		ICONSPO(K)%PROCID=I
	 	END IF
		END DO
!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!-----------------------------!
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
			ICPUID=N
		K=0
		DO I=0,ISIZE-1
			IF (I.NE.N)THEN
					K=K+1
    CALL MPI_SENDRECV(ICONRPA(N)%HOWMANYI(1),1,MPI_INTEGER,I,ICPUID,&
    ICONSPO(K)%HOWMANYTHEY(1),1,MPI_INTEGER,ICONSPO(K)%PROCID,ICONSPO(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
			END IF
		END DO


		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		DO I=1,ISIZE-1
		ALLOCATE(ICONSPO(I)%WHICHTHEY(ICONSPO(I)%HOWMANYTHEY(1),1))
		ALLOCATE(ICONSPO(I)%FACX(ICONSPO(I)%HOWMANYTHEY(1),3))
		ICONSPO(I)%WHICHTHEY(:,:)=0
		ICONSPO(I)%FACX(:,:)=0.d0
		END DO
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
			ICPUID=N
		K=0
DO I=0,ISIZE-1
IF (I.NE.N)THEN
  K=K+1
  CALL MPI_SENDRECV(ICONRPA(N)%FACX(1:ICONRPA(N)%HOWMANYI(1),1:2),&
ICONRPA(N)%HOWMANYI(1)*2,MPI_DOUBLE_PRECISION,I,ICPUID,&
  ICONSPO(K)%FACX(1:ICONSPO(K)%HOWMANYTHEY(1),1:2),ICONSPO(K)%HOWMANYTHEY(1)*2,&
MPI_DOUBLE_PRECISION,ICONSPO(K)%PROCID,ICONSPO(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
END IF
END DO
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		K=0
		DO J=0,ISIZE-1
			IF (J.NE.N)THEN
			K=K+1
			DO I=1,ICONSPO(K)%HOWMANYTHEY(1)
			  VEXT(1,1)=ICONSPO(K)%FACX(I,1); VEXT(1,2)=ICONSPO(K)%FACX(I,2)


			DO KJ=1,ICONRPA(N)%HOWMANYI(1)
			  VEXT(2,1)=ICONRPA(N)%FACX(KJ,1); VEXT(2,2)=ICONRPA(N)%FACX(KJ,2)

  dist=distance2(n)

 if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
				      ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall))then


! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL))then
						
				     ICONSPO(K)%WHICHTHEY(I,1)=ICONRPA(N)%WHICHI(KJ,3);go to 901
				      end if
				      end if
				      if (((abs(vext(2,2)-yper).lt.tolsmall).or.(abs((abs(vext(2,2)-yper))-yper).lt.tolsmall)).and.&
				      ((abs(vext(1,2)-yper).lt.tolsmall).or.(abs((abs(vext(1,2)-yper))-yper).lt.tolsmall)))then
				      if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall))then



						
				      ICONSPO(K)%WHICHTHEY(I,1)=ICONRPA(N)%WHICHI(KJ,3);go to 901
				      end if
				      end if




			    END DO
			    901 CONTINUE
			    END DO
			END IF
		END DO
				
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		
		      !-------------------FOR DEBUGGING ONLY -----------------------------------------!
		K=0
		!-------------------FOR DEBUGGING ONLY -----------------------------------------!
! 		
		!-------------------FOR DEBUGGING ONLY -----------------------------------------!
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
			ICPUID=N
		K=0
DO I=0,ISIZE-1
IF (I.NE.N)THEN
  K=K+1
  CALL MPI_SENDRECV(ICONSPO(K)%WHICHTHEY(1:ICONSPO(K)%HOWMANYTHEY(1),1:1),&
ICONSPO(K)%HOWMANYTHEY(1),MPI_INTEGER,I,ICPUID,&
  ICONRPF(K)%WHICHI(1:ICONRPF(K)%HOWMANYI(1),1:1),ICONRPF(K)%HOWMANYI(1),&
MPI_INTEGER,ICONSPO(K)%PROCID,ICONSPO(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
END IF
END DO
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
K=0
! 		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		DO I=0,ISIZE-1
			IF (I.NE.N)THEN
			K=K+1
			DO KJ=1,ICONRPF(K)%HOWMANYI(1)
			IF (ICONRPF(K)%WHICHI(KJ,1).GT.0)THEN
! 			
 			IELEM(N,ICONRPA(N)%WHICHI(KJ,1))%Ineighg(ICONRPA(N)%WHICHI(KJ,2))=ICONRPF(K)%WHICHI(KJ,1)
			END IF
			END DO
			END IF
		END DO
		KJ=0
		DO I=1,ICONRPA(N)%HOWMANYI(1)
			
			IF (IELEM(N,ICONRPA(N)%WHICHI(I,1))%Ineighg(ICONRPA(N)%WHICHI(I,2)).EQ.0)THEN
			  if (((IBOUND(N,IELEM(N,ICONRpa(N)%WHICHI(I,1))%IBOUNDS(ICONRpa(N)%WHICHI(I,2)))%ICODE)).eq.5)then
			KJ=KJ+1
				  
			      
			END IF
			end if
		END DO
! 		
		if (kj.gt.0)then
		  PRINT*,"error detected in boundary conditions matching rules"
		end if
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 		
		DEALLOCATE(ICONRPA)
		DEALLOCATE(ICONRPF)
		DEALLOCATE(ICONRPM)
		DEALLOCATE(ICONSPO)		
		END IF
		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		
		DEALLOCATE(ICONS)
		DEALLOCATE(ICONR)

! 		

	
			

END SUBROUTINE CONS2d

!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DETSTENX(N,ISIZE,IELEM,ISELEMT,XMPIELRANK,ILOCALALLS,TYPESTEN,&
ILOCALALLELG,STCON,STCONC,STCONS,STCONG,ISOSA,IX,IISTART,IFSAT,PARE,DOSE,PAREEL,DOSEEL,&
PARES,SOSEEL,XMPIE,IFIN,TFIN,XMPIL,GLNEIGH)
!> @brief
!> This subroutine builds the large stencils
IMPLICIT NONE
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INout)::IELEM
INTEGER,INTENT(IN)::N,TYPESTEN,ISIZE
INTEGER::I,K,JJJ,KJ,J,L,IM,IO,IR,KMAXE,ICPUID,ITRUE,IC,PRINT_I,kx,kxk,kk,CANDID,CANDID2
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IFIN,TFIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::STCON,STCONC,STCONS,STCONG,ISOSA,IISTART,IFSAT,IX
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLS
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::PARE,DOSE,PAREEL,PARES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::DOSEEL,SOSEEL
integer::pik,jk,ibn,igd1,icomp_set,i2comp_set,extf2,test6
real::distf,x_ste,y_ste,z_ste,max_sten,min_sten,MAX_STEN2,min_sten2,rcomp_set,TESTDIST,testdiv
! CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

KMAXE=XMPIELRANK(N)






ICPUID=N

STCON(N:N)=0
STCONC(N:N)=0
STCONS(N:N)=0
STCONG(N:N)=0
ISOSA(N:N)=0
IFSAT(N:N)=0
IISTART(N:N)=0
IX(N:N)=0
PRINT_I=KMAXE/20

!$OMP DO SCHEDULE (STATIC)
DO K=1,KMAXE



! 	ILOCALALLEL(N,K,1,1)=IELEM(N,K)%IHEX
	ILOCALALLELG(N,K,1,1)=IELEM(N,K)%IHEXGL
! 	ILOCALALLS(N,K,1,1)=IELEM(N,K)%IFCA
	STCON(N)=k
	STCONC(N)=k
	STCONS(N)=0
	STCONG(N)=IELEM(N,K)%IHEXGL
	ISOSA(N)=1
	IISTART(N)=1
	IX(N)=0
	IC=N
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
	CALL ALLSx(IC,ISIZE,IELEM,STCON,STCONC,STCONS,STCONG,ISELEMT,ILOCALALLS,&
ILOCALALLELG,XMPIE,ISOSA,IFSAT,IISTART,IX,XMPIELRANK,PARE,DOSE,PAREEL,DOSEEL,PARES,SOSEEL,IFIN,TFIN,XMPIL,GLNEIGH)


!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
        !TYPE NOW!
        
       
        
! 	if (iweno.eq.1)then
	
      IF (ICOMPACT.ge.1)THEN      !IF ONE
	
                    if (dimensiona.eq.3)then        !IF DIMENSIONA
                                        ibn=ielem(n,k)%ifca
                                        ILOCALALLELG3(1,1:ISELEM)=ZERO
                                        
                                        ILOCALALLELG3(1,1:ISELEMT(N)-ielem(n,k)%ifca)=ILOCALALLELG(N,K,1,1:ISELEMT(N)-ibn)
                                        ILOCALALLELGD(:,:)=zero
                                                        do j=2,iselemt(n)-ielem(n,k)%ifca
							    CANDID2=ILOCALALLELG(N,K,1,J)
							    IF (CAND(XMPIE(CANDID2)).GT.0)THEN
                                                            
                                                            x_ste=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),1)
                                                            y_ste=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),2)
                                                            z_ste=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),3)
                                                            if (iperiodicity.eq.1)then      !IF THREE
                                                                    
                                                                        IF(ABS(x_ste-IELEM(N,K)%XXC).GT.XPER*oo2)THEN
                                                                        x_ste=x_ste+(XPER*SIGN(1.0,IELEM(N,K)%XXC-XPER*oo2))
                                                                        end if
                                                                        IF(ABS(y_ste-IELEM(N,K)%yyC).GT.yPER*oo2)THEN
                                                                        y_ste=y_ste+(yPER*SIGN(1.0,IELEM(N,K)%yyC-yPER*oo2))
                                                                        end if
                                                                        IF(ABS(z_ste-IELEM(N,K)%zzC).GT.zPER*oo2)THEN
                                                                        z_ste=z_ste+(zPER*SIGN(1.0,IELEM(N,K)%zzC-zPER*oo2))
                                                                        end if
                                                            
                                                            end if          !IF THREE
                                                        ILOCALALLELGD(1,J)=SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
                                                ((y_ste-IELEM(N,K)%YYC)**2)+((z_ste-IELEM(N,K)%zzC)**2)+(tolsmall*j))
							    ELSE
							    WRITE(99+N,*)"KILLED HERE DUE TO STENCILS REQUIRING MORE CPUS"
							    STOP
							    
							    
							    END IF
                                                        END DO
                                        
                                        max_sten2=tolsmall
                                        min_sten2=tolbig
                                        
                                        DO KK=1,ilx
                                        
                                                MIN_STEN2=MIN(MIN_STEN2,ILOCALALLELGD(1,KK))
                                        
                                                MAX_STEN2=MAX(MAX_STEN2,ILOCALALLELGD(1,KK))
                                        END DO
	
	
      
                                    igd1=1
                                    testdist=zero
                                    max_sten=tolsmall
                                    min_sten=tolbig
                                                                            do j=1,ielem(n,k)%ifca
                                                    if (ielem(n,K)%ineighg(j).gt.0)then
                                                    igd1=igd1+1
                                                    testdist=testdist+ILOCALALLELGD(1,igd1)
                                                        if (ILOCALALLELGD(1,igd1).ge.max_sten)then
                                                        max_sten=ILOCALALLELGD(1,igd1)
                                                        
                                                        end if
                                                        if (ILOCALALLELGD(1,igd1).le.min_sten)then
                                                        min_sten=ILOCALALLELGD(1,igd1)
                                                        end if
                                                    end if
                                                end do
                                      testdiv=igd1-1
	TESTDIST=(testdist)/(testdiv*2)
! 	    min_sten=2.0*IELEM(N,K)%MINEDGE
! 	        ielem(n,K)%vortex(1)=max_sten/min_sten
                                            if (icompact.eq.2)then   !IF ICOMPACT 
                                ! 	       IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS)
                                        imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
                                        
                                        kxk=0
                                        
                                                            do J=ielem(n,k)%ifca+2,iselemt(n)-ibn
                                                            do JK=iselemt(n)-ibn,J+1,-1
                                                                        if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                                        pIK= ILOCALALLELG3(1,JK-1)
                                                                        distf=ILOCALALLELGD(1,JK-1)
                                                                            ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                                        ILOCALALLELG3(1,JK) = pik
                                                                            ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                                        ILOCALALLELGD(1,JK)=distf
                                                                        endif
                                                            enddo ! j
                                                            enddo  ! i
                                                                            
                                                                            
                                                                            
                                            ILOCALALLELG(N,K,1,ielem(n,k)%ifca+2:ISELEMT(N))=ILOCALALLELG3(1,ielem(n,k)%ifca+2:ISELEMT(N))
                                            
                                            
                                            
                                            
                                            eND IF           !IF THREE
	      
	      
	      
                                                if (icompact.eq.3)then   !IF ICOMPACT 3
                                        kxk=0
                                        icomp_set=igd1+1
                                !  	  
                                                                    IF ((TESTDIST.LE.(1.4*ielem(n,k)%minedge)).AND.(MAX_STEN/MIN_STEN.LE.(1.4)))THEN
                                                                ielem(n,K)%vortex(1)=100.0
                                                                
                                                                                                
                                                                                                
                                                                                                            do J=icomp_set,iselemt(n)-ibn
                                                                                                            do JK=iselemt(n)-icomp_set,J+1,-1
                                                                                                        if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                                                                        pIK= ILOCALALLELG3(1,JK-1)
                                                                                                        distf=ILOCALALLELGD(1,JK-1)
                                                                                                            ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                                                                        ILOCALALLELG3(1,JK) = pik
                                                                                                            ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                                                                        ILOCALALLELGD(1,JK)=distf
                                                                                                        endif
                                                                                                    enddo ! j
                                                                                                    enddo  ! i
                                                                                                    
                                                                                                    
                                                                        max_sten2=((TESTDIST)*(iorder+1))  
                                                                        
                                                                            do kk=icomp_set,ISELEMT(N)-ibn
                                                                            if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
                                                                                
                                                                                kxk=kxk+1
                                                                            
                                                                            end if
                                                                            end do
                                                                            
                                                                    
                                                                    END IF     !IF THREE
                                            test6=IELEM(N,K)%iNUMNEIGHBOURS*1.2
                                            
                                            IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS,MIN(kxk+icomp_set,test6))
                                        imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
                                        ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))
                                        
                                        
                                        
                                        
                                        end if      !IF Three
                                        if ((icompact.eq.1).AND.((MAX_STEN/MIN_STEN).LE.1.2))then      !IF ICOMPACT 1
                        
                                            extf2=extf      
                                            icomp_set=(Ilx*extf2/2)+1
                                            
                                            kxk=0
                                        
                                                                do J=icomp_set,iselemt(n)-ibn
                                                                do JK=iselemt(n)-icomp_set-ibn,J+1,-1
                                                                        if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                                        pIK= ILOCALALLELG3(1,JK-1)
                                                                        distf=ILOCALALLELGD(1,JK-1)
                                                                            ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                                        ILOCALALLELG3(1,JK) = pik
                                                                            ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                                        ILOCALALLELGD(1,JK)=distf
                                                                        endif
                                                        enddo ! j
                                                        enddo  ! i
                                            
                                            
                                                            do kk=icomp_set,ISELEMT(N)-ibn
                                                                if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
                                            !                          
                                                                    kxk=kxk+1
                                                                
                                                                end if
                                                                end do
                                !  	

                                        IELEM(N,K)%iNUMNEIGHBOURS=max(IELEM(N,K)%iNUMNEIGHBOURS,min(kxk+icomp_set,ilx*extf2))
                                        imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
                                        ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))

                                        
                                !  	 END IF
                                        
                                            
                                        end if
	
	   
	   
	   
	  else
     
      
         
	 
	   ibn=ielem(n,k)%ifca
	ILOCALALLELG3(1,1:ISELEM)=ZERO
	
	ILOCALALLELG3(1,1:ISELEMT(N)-ielem(n,k)%ifca)=ILOCALALLELG(N,K,1,1:ISELEMT(N)-ibn)
! 	
	ILOCALALLELGD(:,:)=zero
	do j=2,iselemt(n)-ielem(n,k)%ifca
	
	    CANDID2=ILOCALALLELG(N,K,1,J)
		IF (CAND(XMPIE(CANDID2)).GT.0)THEN
                 
                 x_ste=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),1)
                 y_ste=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),2)
                                                            
	
	
	
	   
	    if (iperiodicity.eq.1)then
		     
			IF(ABS(x_ste-IELEM(N,K)%XXC).GT.XPER*oo2)THEN
			x_ste=x_ste+(XPER*SIGN(1.0,IELEM(N,K)%XXC-XPER*oo2))
			end if
			IF(ABS(y_ste-IELEM(N,K)%yyC).GT.yPER*oo2)THEN
			y_ste=y_ste+(yPER*SIGN(1.0,IELEM(N,K)%yyC-yPER*oo2))
			end if
	    
	    
	    end if
! 	    IF (ICOMPACT.EQ.1)THEN
	  ILOCALALLELGD(1,J)=(SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
                   ((y_ste-IELEM(N,K)%YYC)**2)))+(tolsmall)
!             ELSE
!             ILOCALALLELGD(1,J)=MAX(((SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
!                    ((y_ste-IELEM(N,K)%YYC)**2)+(tolsmall)))/(SQRT((x_ste-IELEM(N,K)%XXC)**2)+(tolsmall))),((SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
!                    ((y_ste-IELEM(N,K)%YYC)**2)+(tolsmall)))/(SQRT((Y_ste-IELEM(N,K)%YYC)**2)+(tolsmall))))
!             
!             END IF
	  ELSE
	  WRITE(99+N,*)"KILLED HERE DUE TO STENCILS REQUIRING MORE CPUS"
							    STOP
	  
	  
	  END IF
	END DO
	
	
	max_sten2=tolsmall
	min_sten2=tolbig
	
	
        
        
          DO KK=1,ILX
                MAX_STEN2=MAX(MAX_STEN2,ILOCALALLELGD(1,KK))
                MIN_STEN2=MIN(MIN_STEN2,ILOCALALLELGD(1,KK))
          END DO
!           
            
      
	igd1=1
	testdist=zero
	max_sten=tolsmall
	min_sten=tolbig
	do j=1,ielem(n,k)%ifca
	      if (ielem(n,K)%ineighg(j).gt.0)then
	      igd1=igd1+1
	      testdist=testdist+ILOCALALLELGD(1,igd1)
		if (ILOCALALLELGD(1,igd1).ge.max_sten)then
		  max_sten=ILOCALALLELGD(1,igd1)
		  
		end if
		if (ILOCALALLELGD(1,igd1).le.min_sten)then
		  min_sten=ILOCALALLELGD(1,igd1)
		end if
	      end if
	end do
	testdiv=igd1-1
	TESTDIST=(testdist)/(testdiv*2)
! 	    min_sten=2.0*IELEM(N,K)%MINEDGE
! 	        ielem(n,K)%vortex(1)=max_sten/min_sten
	      if (icompact.eq.2)then
! 	       IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS)
 	  imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
 	  
 	   kxk=0
	
                    do J=ielem(n,k)%ifca+2,iselemt(n)-ibn
                    do JK=iselemt(n)-ibn,J+1,-1
                if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                pIK= ILOCALALLELG3(1,JK-1)
                distf=ILOCALALLELGD(1,JK-1)
                    ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                ILOCALALLELG3(1,JK) = pik
                    ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                ILOCALALLELGD(1,JK)=distf
                endif
            enddo ! j
            enddo  ! i
 	  
 	  
	   ILOCALALLELG(N,K,1,ielem(n,k)%ifca+2:ISELEMT(N))=ILOCALALLELG3(1,ielem(n,k)%ifca+2:ISELEMT(N))
	      
	      
	      
	      
	      eND IF
	      
	        if (icompact.eq.3)then
	   kxk=0
	   icomp_set=igd1+1

            IF ((TESTDIST.LE.(1.4*ielem(n,k)%minedge)).AND.(MAX_STEN/MIN_STEN.LE.(1.4)))THEN
	   ielem(n,K)%vortex(1)=100.0
	   
                                        
                                        
                                                    do J=icomp_set,iselemt(n)-ibn
                                                    do JK=iselemt(n)-icomp_set,J+1,-1
                                                if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                pIK= ILOCALALLELG3(1,JK-1)
                                                distf=ILOCALALLELGD(1,JK-1)
                                                    ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                ILOCALALLELG3(1,JK) = pik
                                                    ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                ILOCALALLELGD(1,JK)=distf
                                                endif
                                            enddo ! j
                                            enddo  ! i
                                            
                                            
                max_sten2=((TESTDIST)*(iorder+1))  
                
	     do kk=icomp_set,ISELEMT(N)-ibn
                    if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
                         
                        kxk=kxk+1
                       
                    end if
                    end do
                    
 	    
	    END IF
	    test6=IELEM(N,K)%iNUMNEIGHBOURS*2.6
	    
	    IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS,MIN(kxk+icomp_set,test6))
 	  imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
	   ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))
	   
	   
	   
	   
	   end if
	      
	      
	     if ((icompact.eq.1).AND.((MAX_STEN/MIN_STEN).LE.1.2))then   

            extf2=extf          !EXTENSION
	    icomp_set=(Ilx*extf2/2)+3      
	   
	    kxk=0
	
                    do J=icomp_set,iselemt(n)-ibn
                    do JK=iselemt(n)-icomp_set-ibn,J+1,-1
                if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                pIK= ILOCALALLELG3(1,JK-1)
                distf=ILOCALALLELGD(1,JK-1)
                    ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                ILOCALALLELG3(1,JK) = pik
                    ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                ILOCALALLELGD(1,JK)=distf
                endif
            enddo ! j
            enddo  ! i
            
            
                   do kk=icomp_set,ISELEMT(N)-ibn
                    if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
!                          
                        kxk=kxk+1
                       
                    end if
                    end do


         IELEM(N,K)%iNUMNEIGHBOURS=max(IELEM(N,K)%iNUMNEIGHBOURS,min(kxk+icomp_set,ilx*extf2))
 	  imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
	   ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))

 	 
!  	 END IF
 	 
	    
	end if
	end if
        
       
    END IF


	
END DO
!$OMP END DO
    
 
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

END SUBROUTINE DETSTENX
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECURSIVE SUBROUTINE ALLSX(N,ISIZE,IELEM,STCON,STCONC,STCONS,STCONG,ISELEMT,ILOCALALLS,&
ILOCALALLELG,XMPIE,ISOSA,IFSAT,IISTART,IX,XMPIELRANK,PARE,DOSE,PAREEL,DOSEEL,PARES,SOSEEL,IFIN,TFIN,XMPIL,GLNEIGH)
!> @brief
!> This recursive subroutine builds the large stencils until the prescribed number of elements is reached
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ISIZE
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::IELEM
INTEGER::I,K,JJJ,KJ,J,L,IM,IO,IR,KMAXE,ICPUID,ITRUE,jloop,CANDID
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IFIN,TFIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::STCON,STCONC,STCONS,STCONG,ISOSA,IFSAT,IISTART,IX
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLS
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::PARE,DOSE,PAREEL,PARES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::DOSEEL,SOSEEL
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
KMAXE=XMPIELRANK(N)
      if (dimensiona.eq.3)then
			  jloop=6
		else
			  jloop=4
		end if
		ICPUID=N
		IISTART(N)=IISTART(N)+1
		IF (XMPIE(STCONG(N)).EQ.N)THEN		!global array checking if element belongs to my cpu 
			L=XMPIL(STCONG(N))		!if yes then set l=to the number of this cpu
			STCONC(N)=L			!use this stconc(n)=l for referencing after
		
		DO J=1,IELEM(N,STCONC(N))%IFCA		!loop all the sides of this element
			IF (IELEM(N,STCONC(N))%INEIGHG(J).GT.0)THEN  !if the neighbour is gt.0 then
			IX(N)=IELEM(N,STCONC(N))%INEIGHG(J)		!set the ix(n) as the global index of this element
			CALL CHECK(N,STCON,IX,ILOCALALLELG,IFSAT,ISELEMT)	!check if this element is already in the list
			IF (IFSAT(N).EQ.1)THEN					!if not then include
			if (isosa(n).le.ISELEMT(N)-1)then
			ISOSA(N)=ISOSA(N)+1
			ILOCALALLELG(N,STCON(N),1,ISOSA(N))=IELEM(N,STCONC(N))%INEIGHG(J)
			END IF
			end if
			END IF
! 			
		END DO
		END IF
		
		IF (XMPIE(STCONG(N)).NE.N)THEN	!if this element belongs to another cpu then
		      IF (CAND(XMPIE(STCONG(N))).GT.0)THEN
		      DO J=1,jloop			!loop all the sides of this particular element (i need to find the cpu that it belongs) therefore indexing to parent cpu
			      CANDID=CAND2R(CAND(XMPIE(STCONG(N))),XMPIL(STCONG(N)),J)
			
			      IF(CANDID.GT.0) THEN
			      
			      IX(N)=CANDID
			      CALL CHECK(N,STCON,IX,ILOCALALLELG,IFSAT,ISELEMT)
			      IF (IFSAT(N).EQ.1)THEN
			      if (isosa(n).le.ISELEMT(N)-1)then
			      ISOSA(N)=ISOSA(N)+1
			      ILOCALALLELG(N,STCON(N),1,ISOSA(N))=CANDID
			      END IF
			      END IF
			      END IF
		      END DO
		      ELSE
		
		      
		      
		      WRITE(99+N,*)"KILLED DUE TO STENCILS REQUIRING MORE CPUS"
		      STOP
		      END IF
		
! 		
		END IF
		IF (ISOSA(N).LT.ISELEMT(N))THEN
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
! 			
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
			STCONG(N)=ILOCALALLELG(N,STCON(N),1,IISTART(N))
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
! 			
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
			CALL ALLSX (N,ISIZE,IELEM,STCON,STCONC,STCONS,STCONG,ISELEMT,ILOCALALLS,&
ILOCALALLELG,XMPIE,ISOSA,IFSAT,IISTART,IX,XMPIELRANK,PARE,DOSE,PAREEL,DOSEEL,PARES,SOSEEL,IFIN,TFIN,XMPIL,GLNEIGH)
		END IF
		IF (ISOSA(N).eq.ISELEMT(N))THEN
			RETURN
		END IF

			
END SUBROUTINE ALLSX



SUBROUTINE DETSTEN(N,ISIZE,IELEM,ISELEMT,XMPIELRANK,ILOCALALLS,TYPESTEN,&
ILOCALALLELG,STCON,STCONC,STCONS,STCONG,ISOSA,IX,IISTART,IFSAT,PARE,DOSE,PAREEL,DOSEEL,&
PARES,SOSEEL,XMPIE,IFIN,TFIN,XMPIL,GLNEIGH)
!> @brief
!> This subroutine builds the large stencils suitable for periodic boundaries
IMPLICIT NONE
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INout)::IELEM
INTEGER,INTENT(IN)::N,TYPESTEN,ISIZE
INTEGER::I,K,JJJ,KJ,J,L,IM,IO,IR,KMAXE,ICPUID,ITRUE,IC,PRINT_I,kx,kxk,kk
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IFIN,TFIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::STCON,STCONC,STCONS,STCONG,ISOSA,IISTART,IFSAT,IX
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLS
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::PARE,DOSE,PAREEL,PARES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::DOSEEL,SOSEEL
integer::pik,jk,ibn,igd1,icomp_set,i2comp_set,extf2,test6
real::distf,x_ste,y_ste,z_ste,max_sten,min_sten,MAX_STEN2,min_sten2,rcomp_set,TESTDIST,testdiv

KMAXE=XMPIELRANK(N)






ICPUID=N

STCON(N:N)=0
STCONC(N:N)=0
STCONS(N:N)=0
STCONG(N:N)=0
ISOSA(N:N)=0
IFSAT(N:N)=0
IISTART(N:N)=0
IX(N:N)=0
PRINT_I=KMAXE/20

!$OMP DO SCHEDULE (STATIC)
DO K=1,KMAXE



! 	ILOCALALLEL(N,K,1,1)=IELEM(N,K)%IHEX
	ILOCALALLELG(N,K,1,1)=IELEM(N,K)%IHEXGL
! 	ILOCALALLS(N,K,1,1)=IELEM(N,K)%IFCA
	STCON(N)=k
	STCONC(N)=k
	STCONS(N)=0
	STCONG(N)=IELEM(N,K)%IHEXGL
	ISOSA(N)=1
	IISTART(N)=1
	IX(N)=0
	IC=N
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
	CALL ALLS(IC,ISIZE,IELEM,STCON,STCONC,STCONS,STCONG,ISELEMT,ILOCALALLS,&
ILOCALALLELG,XMPIE,ISOSA,IFSAT,IISTART,IX,XMPIELRANK,PARE,DOSE,PAREEL,DOSEEL,PARES,SOSEEL,IFIN,TFIN,XMPIL,GLNEIGH)


        
        
        
! 	if (iweno.eq.1)then
	
      IF (ICOMPACT.ge.1)THEN      !IF ONE
	
                    if (dimensiona.eq.3)then        !IF DIMENSIONA
                                        ibn=ielem(n,k)%ifca
                                        ILOCALALLELG3(1,1:ISELEM)=ZERO
                                        
                                        ILOCALALLELG3(1,1:ISELEMT(N)-ielem(n,k)%ifca)=ILOCALALLELG(N,K,1,1:ISELEMT(N)-ibn)
                                        ILOCALALLELGD(:,:)=zero
                                                        do j=2,iselemt(n)-ielem(n,k)%ifca
                                                            
                                                            x_ste=CENTERR(ILOCALALLELG(N,K,1,J),1);y_ste=CENTERR(ILOCALALLELG(N,K,1,J),2);z_ste=CENTERR(ILOCALALLELG(N,K,1,J),3)
                                                            if (iperiodicity.eq.1)then      !IF THREE
                                                                    
                                                                        IF(ABS(x_ste-IELEM(N,K)%XXC).GT.XPER*oo2)THEN
                                                                        x_ste=x_ste+(XPER*SIGN(1.0,IELEM(N,K)%XXC-XPER*oo2))
                                                                        end if
                                                                        IF(ABS(y_ste-IELEM(N,K)%yyC).GT.yPER*oo2)THEN
                                                                        y_ste=y_ste+(yPER*SIGN(1.0,IELEM(N,K)%yyC-yPER*oo2))
                                                                        end if
                                                                        IF(ABS(z_ste-IELEM(N,K)%zzC).GT.zPER*oo2)THEN
                                                                        z_ste=z_ste+(zPER*SIGN(1.0,IELEM(N,K)%zzC-zPER*oo2))
                                                                        end if
                                                            
                                                            end if          !IF THREE
                                                        ILOCALALLELGD(1,J)=SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
                                                ((y_ste-IELEM(N,K)%YYC)**2)+((z_ste-IELEM(N,K)%zzC)**2)+(tolsmall*j))
                                                        END DO
                                        
                                        max_sten2=tolsmall
                                        min_sten2=tolbig
                                        
                                        DO KK=1,ilx
                                        
                                                MIN_STEN2=MIN(MIN_STEN2,ILOCALALLELGD(1,KK))
                                        
                                                MAX_STEN2=MAX(MAX_STEN2,ILOCALALLELGD(1,KK))
                                        END DO
	
	
      
                                    igd1=1
                                    testdist=zero
                                    max_sten=tolsmall
                                    min_sten=tolbig
                                                                            do j=1,ielem(n,k)%ifca
                                                    if (ielem(n,K)%ineighg(j).gt.0)then
                                                    igd1=igd1+1
                                                    testdist=testdist+ILOCALALLELGD(1,igd1)
                                                        if (ILOCALALLELGD(1,igd1).ge.max_sten)then
                                                        max_sten=ILOCALALLELGD(1,igd1)
                                                        
                                                        end if
                                                        if (ILOCALALLELGD(1,igd1).le.min_sten)then
                                                        min_sten=ILOCALALLELGD(1,igd1)
                                                        end if
                                                    end if
                                                end do
                                      testdiv=igd1-1
	TESTDIST=(testdist)/(testdiv*2)
! 	    min_sten=2.0*IELEM(N,K)%MINEDGE
! 	        ielem(n,K)%vortex(1)=max_sten/min_sten
                                            if (icompact.eq.2)then   !IF ICOMPACT 
                                ! 	       IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS)
                                        imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
                                        
                                        kxk=0
                                        
                                                            do J=ielem(n,k)%ifca+2,iselemt(n)-ibn
                                                            do JK=iselemt(n)-ibn,J+1,-1
                                                                        if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                                        pIK= ILOCALALLELG3(1,JK-1)
                                                                        distf=ILOCALALLELGD(1,JK-1)
                                                                            ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                                        ILOCALALLELG3(1,JK) = pik
                                                                            ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                                        ILOCALALLELGD(1,JK)=distf
                                                                        endif
                                                            enddo ! j
                                                            enddo  ! i
                                                                            
                                                                            
                                                                            
                                            ILOCALALLELG(N,K,1,ielem(n,k)%ifca+2:ISELEMT(N))=ILOCALALLELG3(1,ielem(n,k)%ifca+2:ISELEMT(N))
                                            
                                            
                                            
                                            
                                            eND IF           !IF THREE
	      
	      
	      
                                                if (icompact.eq.3)then   !IF ICOMPACT 3
                                        kxk=0
                                        icomp_set=igd1+1
                                !  	   
                                                                    IF ((TESTDIST.LE.(1.4*ielem(n,k)%minedge)).AND.(MAX_STEN/MIN_STEN.LE.(1.4)))THEN
                                                                ielem(n,K)%vortex(1)=100.0
                                                                
                                                                                                
                                                                                                
                                                                                                            do J=icomp_set,iselemt(n)-ibn
                                                                                                            do JK=iselemt(n)-icomp_set,J+1,-1
                                                                                                        if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                                                                        pIK= ILOCALALLELG3(1,JK-1)
                                                                                                        distf=ILOCALALLELGD(1,JK-1)
                                                                                                            ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                                                                        ILOCALALLELG3(1,JK) = pik
                                                                                                            ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                                                                        ILOCALALLELGD(1,JK)=distf
                                                                                                        endif
                                                                                                    enddo ! j
                                                                                                    enddo  ! i
                                                                                                    
                                                                                                    
                                                                        max_sten2=((TESTDIST)*(iorder+1))  
                                                                        
                                                                            do kk=icomp_set,ISELEMT(N)-ibn
                                                                            if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
                                                                                
                                                                                kxk=kxk+1
                                                                            
                                                                            end if
                                                                            end do
                                                                            
                                                                    
                                                                    END IF     !IF THREE
                                            test6=IELEM(N,K)%iNUMNEIGHBOURS*1.2
                                            
                                            IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS,MIN(kxk+icomp_set,test6))
                                        imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
                                        ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))
                                        
                                        
                                        
                                        
                                        end if      !IF Three
                                        if ((icompact.eq.1).AND.((MAX_STEN/MIN_STEN).LE.1.2))then      !IF ICOMPACT 1
                        
                                            extf2=extf      
                                            icomp_set=(Ilx*extf2/2)+1
                                            
                                            kxk=0
                                        
                                                                do J=icomp_set,iselemt(n)-ibn
                                                                do JK=iselemt(n)-icomp_set-ibn,J+1,-1
                                                                        if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                                        pIK= ILOCALALLELG3(1,JK-1)
                                                                        distf=ILOCALALLELGD(1,JK-1)
                                                                            ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                                        ILOCALALLELG3(1,JK) = pik
                                                                            ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                                        ILOCALALLELGD(1,JK)=distf
                                                                        endif
                                                        enddo ! j
                                                        enddo  ! i
                                            
                                            
                                                            do kk=icomp_set,ISELEMT(N)-ibn
                                                                if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
                                            !                          
                                                                    kxk=kxk+1
                                                                
                                                                end if
                                                                end do
                                !  	 

                                        IELEM(N,K)%iNUMNEIGHBOURS=max(IELEM(N,K)%iNUMNEIGHBOURS,min(kxk+icomp_set,ilx*extf2))
                                        imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
                                        ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))

                                        
                                !  	 END IF
                                        
                                            
                                        end if
	
	   
	   
	   
	  else
     
      
         
	 
	   ibn=ielem(n,k)%ifca
	ILOCALALLELG3(1,1:ISELEM)=ZERO
	
	ILOCALALLELG3(1,1:ISELEMT(N)-ielem(n,k)%ifca)=ILOCALALLELG(N,K,1,1:ISELEMT(N)-ibn)
! 	
	ILOCALALLELGD(:,:)=zero
	do j=2,iselemt(n)-ielem(n,k)%ifca
	    x_ste=CENTERR(ILOCALALLELG(N,K,1,J),1);y_ste=CENTERR(ILOCALALLELG(N,K,1,J),2)
	    if (iperiodicity.eq.1)then
		     
			IF(ABS(x_ste-IELEM(N,K)%XXC).GT.XPER*oo2)THEN
			x_ste=x_ste+(XPER*SIGN(1.0,IELEM(N,K)%XXC-XPER*oo2))
			end if
			IF(ABS(y_ste-IELEM(N,K)%yyC).GT.yPER*oo2)THEN
			y_ste=y_ste+(yPER*SIGN(1.0,IELEM(N,K)%yyC-yPER*oo2))
			end if
	    
	    
	    end if
! 	    IF (ICOMPACT.EQ.1)THEN
	  ILOCALALLELGD(1,J)=(SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
                   ((y_ste-IELEM(N,K)%YYC)**2)))+(tolsmall)
!             ELSE
!             ILOCALALLELGD(1,J)=MAX(((SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
!                    ((y_ste-IELEM(N,K)%YYC)**2)+(tolsmall)))/(SQRT((x_ste-IELEM(N,K)%XXC)**2)+(tolsmall))),((SQRT(((x_ste-IELEM(N,K)%XXC)**2)+&
!                    ((y_ste-IELEM(N,K)%YYC)**2)+(tolsmall)))/(SQRT((Y_ste-IELEM(N,K)%YYC)**2)+(tolsmall))))
!             
!             END IF
	END DO
	
	
	max_sten2=tolsmall
	min_sten2=tolbig
	
	
        
        
          DO KK=1,ILX
                MAX_STEN2=MAX(MAX_STEN2,ILOCALALLELGD(1,KK))
                MIN_STEN2=MIN(MIN_STEN2,ILOCALALLELGD(1,KK))
          END DO
!            
            
            
      
	igd1=1
	testdist=zero
	max_sten=tolsmall
	min_sten=tolbig
	do j=1,ielem(n,k)%ifca
	      if (ielem(n,K)%ineighg(j).gt.0)then
	      igd1=igd1+1
	      testdist=testdist+ILOCALALLELGD(1,igd1)
		if (ILOCALALLELGD(1,igd1).ge.max_sten)then
		  max_sten=ILOCALALLELGD(1,igd1)
		  
		end if
		if (ILOCALALLELGD(1,igd1).le.min_sten)then
		  min_sten=ILOCALALLELGD(1,igd1)
		end if
	      end if
	end do
	testdiv=igd1-1
	TESTDIST=(testdist)/(testdiv*2)
! 	    min_sten=2.0*IELEM(N,K)%MINEDGE
! 	        ielem(n,K)%vortex(1)=max_sten/min_sten
	      if (icompact.eq.2)then
! 	       IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS)
 	  imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
 	  
 	   kxk=0
	
                    do J=ielem(n,k)%ifca+2,iselemt(n)-ibn
                    do JK=iselemt(n)-ibn,J+1,-1
                if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                pIK= ILOCALALLELG3(1,JK-1)
                distf=ILOCALALLELGD(1,JK-1)
                    ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                ILOCALALLELG3(1,JK) = pik
                    ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                ILOCALALLELGD(1,JK)=distf
                endif
            enddo ! j
            enddo  ! i
 	  
 	  
	   ILOCALALLELG(N,K,1,ielem(n,k)%ifca+2:ISELEMT(N))=ILOCALALLELG3(1,ielem(n,k)%ifca+2:ISELEMT(N))
	      
	      
	      
	      
	      eND IF
	      
	        if (icompact.eq.3)then
	   kxk=0
	   icomp_set=igd1+1
!  	   
            IF ((TESTDIST.LE.(1.4*ielem(n,k)%minedge)).AND.(MAX_STEN/MIN_STEN.LE.(1.4)))THEN
	   ielem(n,K)%vortex(1)=100.0
	   
                                        
                                        
                                                    do J=icomp_set,iselemt(n)-ibn
                                                    do JK=iselemt(n)-icomp_set,J+1,-1
                                                if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                                                pIK= ILOCALALLELG3(1,JK-1)
                                                distf=ILOCALALLELGD(1,JK-1)
                                                    ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                                                ILOCALALLELG3(1,JK) = pik
                                                    ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                                                ILOCALALLELGD(1,JK)=distf
                                                endif
                                            enddo ! j
                                            enddo  ! i
                                            
                                            
                max_sten2=((TESTDIST)*(iorder+1))  
                
	     do kk=icomp_set,ISELEMT(N)-ibn
                    if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
                         
                        kxk=kxk+1
                       
                    end if
                    end do
                    
 	    
	    END IF
	    test6=IELEM(N,K)%iNUMNEIGHBOURS*2.6
	    
	    IELEM(N,K)%iNUMNEIGHBOURS=MAX(IELEM(N,K)%iNUMNEIGHBOURS,MIN(kxk+icomp_set,test6))
 	  imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
	   ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))
	   
	   
	   
	   
	   end if
	      
	      
	     if ((icompact.eq.1).AND.((MAX_STEN/MIN_STEN).LE.1.2))then   

            extf2=extf          !EXTENSION
	    icomp_set=(Ilx*extf2/2)+3      
	   
	    kxk=0
	
                    do J=icomp_set,iselemt(n)-ibn
                    do JK=iselemt(n)-icomp_set-ibn,J+1,-1
                if (ILOCALALLELGD(1,JK-1) .gt. ILOCALALLELGD(1,jK)) then
                pIK= ILOCALALLELG3(1,JK-1)
                distf=ILOCALALLELGD(1,JK-1)
                    ILOCALALLELG3(1,JK-1) = ILOCALALLELG3(1,JK)
                ILOCALALLELG3(1,JK) = pik
                    ILOCALALLELGD(1,JK-1)=ILOCALALLELGD(1,JK)
                ILOCALALLELGD(1,JK)=distf
                endif
            enddo ! j
            enddo  ! i
            
            
                   do kk=icomp_set,ISELEMT(N)-ibn
                    if  ((ILOCALALLELGD(1,kk)).le.MAX_STEN2)then
!                         
                        kxk=kxk+1
                       
                    end if
                    end do
!  	

         IELEM(N,K)%iNUMNEIGHBOURS=max(IELEM(N,K)%iNUMNEIGHBOURS,min(kxk+icomp_set,ilx*extf2))
 	  imaxdegfree=max(imaxdegfree,IELEM(N,K)%iNUMNEIGHBOURS-1)
	   ILOCALALLELG(N,K,1,icomp_set:ISELEMT(N))=ILOCALALLELG3(1,icomp_set:ISELEMT(N))

 	 
!  	 END IF
 	 
	    
	end if
	end if
        
       
    END IF


	
END DO
!$OMP END DO
    
 
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

END SUBROUTINE DETSTEN
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
RECURSIVE SUBROUTINE ALLS(N,ISIZE,IELEM,STCON,STCONC,STCONS,STCONG,ISELEMT,ILOCALALLS,&
ILOCALALLELG,XMPIE,ISOSA,IFSAT,IISTART,IX,XMPIELRANK,PARE,DOSE,PAREEL,DOSEEL,PARES,SOSEEL,IFIN,TFIN,XMPIL,GLNEIGH)
!> @brief
!> This recursive subroutine builds the large stencils suitable for periodic boundaries
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ISIZE
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::IELEM
INTEGER::I,K,JJJ,KJ,J,L,IM,IO,IR,KMAXE,ICPUID,ITRUE,jloop
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IFIN,TFIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::STCON,STCONC,STCONS,STCONG,ISOSA,IFSAT,IISTART,IX
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLS
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::PARE,DOSE,PAREEL,PARES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::DOSEEL,SOSEEL
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::GLNEIGH
KMAXE=XMPIELRANK(N)
      if (dimensiona.eq.3)then
			  jloop=6
		else
			  jloop=4
		end if
		ICPUID=N
		IISTART(N)=IISTART(N)+1
		IF (XMPIE(STCONG(N)).EQ.N)THEN
			L=XMPIL(STCONG(N))
			STCONC(N)=L
		
		DO J=1,IELEM(N,STCONC(N))%IFCA
			IF (IELEM(N,STCONC(N))%INEIGHG(J).GT.0)THEN
			IX(N)=IELEM(N,STCONC(N))%INEIGHG(J)
			CALL CHECK(N,STCON,IX,ILOCALALLELG,IFSAT,ISELEMT)
			IF (IFSAT(N).EQ.1)THEN
			if (isosa(n).le.ISELEMT(N)-1)then
			ISOSA(N)=ISOSA(N)+1
			ILOCALALLELG(N,STCON(N),1,ISOSA(N))=IELEM(N,STCONC(N))%INEIGHG(J)
			END IF
			end if
			END IF
! 			
		END DO
		END IF
		IF (XMPIE(STCONG(N)).NE.N)THEN
		
		DO J=1,jloop
			IF(GLNEIGH(STCONG(N),J).GT.0) THEN
			
			IX(N)=GLNEIGH(STCONG(N),J)
			CALL CHECK(N,STCON,IX,ILOCALALLELG,IFSAT,ISELEMT)
			IF (IFSAT(N).EQ.1)THEN
			if (isosa(n).le.ISELEMT(N)-1)then
			ISOSA(N)=ISOSA(N)+1
			ILOCALALLELG(N,STCON(N),1,ISOSA(N))=GLNEIGH(STCONG(N),J)
			END IF
			END IF
			end if
		END DO
		END IF
		IF (ISOSA(N).LT.ISELEMT(N))THEN
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
! 			
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
			STCONG(N)=ILOCALALLELG(N,STCON(N),1,IISTART(N))
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
! 			
			!-------------------FOR DEBUGGING ONLY -----------------------------------------!
			CALL ALLS (N,ISIZE,IELEM,STCON,STCONC,STCONS,STCONG,ISELEMT,ILOCALALLS,&
ILOCALALLELG,XMPIE,ISOSA,IFSAT,IISTART,IX,XMPIELRANK,PARE,DOSE,PAREEL,DOSEEL,PARES,SOSEEL,IFIN,TFIN,XMPIL,GLNEIGH)
		END IF
		IF (ISOSA(N).eq.ISELEMT(N))THEN
			RETURN
		END IF

			
END SUBROUTINE ALLS



SUBROUTINE CHECK(N,STCON,IX,ILOCALALLELG,IFSAT,ISELEMT)
!> @brief
!> This subroutine checks if some candidate elements already belong to the an existing list of neighbours
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IFSAT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::STCON,ISELEMT,IX
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(IN)::ILOCALALLELG	
INTEGER::I,J,K
IFSAT(N)=1
DO I=1,ISELEMT(N)
IF (IX(N).EQ.ILOCALALLELG(N,STCON(N),1,I))THEN
	IFSAT(N)=0
END IF
END DO
END SUBROUTINE CHECK













SUBROUTINE CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
!> @brief
!> This subroutine checks which candidate cells satisfy the directionality condition for directional stencils
IMPLICIT NONE
INTEGER,INTENT(IN)::N,IPERIODICITY
REAL,INTENT(IN)::XPER,YPER,ZPER
INTEGER,INTENT(IN)::ISSF
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::IWHICHSTEN,ISHYAPE
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ISATISFIED
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::BC
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::VG
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(IN)::VC
REAL::SMALL
INTEGER::isat1,isat2,isat4,isat3








SMALL=tolsmall

ISATISFIED(N)=0

SELECT CASE(n_node)

CASE (4)
        
	 vext(1,:)=BC(N,:)
	 vext(2,1:3)=vext(6,1:3); vext(3,1:3)=vext(5,1:3);vext(4,1:3)=vext(7,1:3);  
	vgg(:)=vg(n,:)
	CALL COMPUTEJACOBIANS
	 IF (IPERIODICITY.EQ.1) THEN
	 if (abs(VG(N,1) - BC(N,1)) .ge. XPer/2.0d0)    VGg(1) = VG(N,1) + XPer*sign(1.D0,BC(N,1) - XPer/2.0d0)
          if (abs(VG(N,2) - BC(N,2)) .ge. YPer/2.0d0)    VGg(2) = VG(N,2) + YPer*sign(1.D0,BC(N,2) - yPer/2.0d0)
	   if (abs(VG(N,3) - BC(N,3)) .ge. zPer/2.0d0)    VGg(3) = VG(N,3) + zPer*sign(1.D0,BC(N,3) - zPer/2.0d0)
         
        END IF
          XCC(1:3) =  matmul(vva1(1:3,1:3),VGg(1:3) -VEXT(1,1:3))
          
           IF ((XcC(1).GE.ZERO).AND.(XcC(2).GE.ZERO).AND.(XcC(3).GE.ZERO))THEN
	    ISAT1=1
	  ELSE
	      ISAT1=0
	  END IF
	  
! 	  IF ((ISAT1.EQ.1))THEN
! 	    ISATISFIED(N)=1
! 	  ELSE
! 	      ISATISFIED(N)=0
! 	  END IF
	  
	  vext(1,:)=BC(N,:)
	 vext(2,1:3)=vext(6,1:3); vext(3,1:3)=vext(8,1:3);vext(4,1:3)=vext(7,1:3);  
	 CALL COMPUTEJACOBIANS
	 IF (IPERIODICITY.EQ.1) THEN
	 if (abs(VG(N,1) - BC(N,1)) .ge. XPer/2.0d0)    VGg(1) = VG(N,1) + XPer*sign(1.D0,BC(N,1) - XPer/2.0d0)
          if (abs(VG(N,2) - BC(N,2)) .ge. YPer/2.0d0)    VGg(2) = VG(N,2) + YPer*sign(1.D0,BC(N,2) - yPer/2.0d0)
	   if (abs(VG(N,3) - BC(N,3)) .ge. zPer/2.0d0)    VGg(3) = VG(N,3) + zPer*sign(1.D0,BC(N,3) - zPer/2.0d0)
         
        END IF
          XCC(1:3) =  matmul(vva1(1:3,1:3),VGg(1:3) -VEXT(1,1:3))
          IF ((XcC(1).GE.ZERO).AND.(XcC(2).GE.ZERO).AND.(XcC(3).GE.ZERO))THEN
	    ISAT2=1
	  ELSE
	      ISAT2=0
	  END IF
	  
	  
	  
	   vext(1,:)=BC(N,:)
	 vext(2,1:3)=vext(5,1:3); vext(3,1:3)=vext(8,1:3);vext(4,1:3)=vext(7,1:3);  
	 CALL COMPUTEJACOBIANS
	 IF (IPERIODICITY.EQ.1) THEN
	 if (abs(VG(N,1) - BC(N,1)) .ge. XPer/2.0d0)    VGg(1) = VG(N,1) + XPer*sign(1.D0,BC(N,1) - XPer/2.0d0)
          if (abs(VG(N,2) - BC(N,2)) .ge. YPer/2.0d0)    VGg(2) = VG(N,2) + YPer*sign(1.D0,BC(N,2) - yPer/2.0d0)
	   if (abs(VG(N,3) - BC(N,3)) .ge. zPer/2.0d0)    VGg(3) = VG(N,3) + zPer*sign(1.D0,BC(N,3) - zPer/2.0d0)
         
        END IF
          XCC(1:3) =  matmul(vva1(1:3,1:3),VGg(1:3) -VEXT(1,1:3))
          IF ((XcC(1).GE.ZERO).AND.(XcC(2).GE.ZERO).AND.(XcC(3).GE.ZERO))THEN
	    ISAT3=1
	  ELSE
	      ISAT3=0
	  END IF
	  
	  
          
          vext(1,:)=BC(N,:)
	 vext(2,1:3)=vext(5,1:3); vext(3,1:3)=vext(8,1:3);vext(4,1:3)=vext(6,1:3);  
	 CALL COMPUTEJACOBIANS
	 IF (IPERIODICITY.EQ.1) THEN
	 if (abs(VG(N,1) - BC(N,1)) .ge. XPer/2.0d0)    VGg(1) = VG(N,1) + XPer*sign(1.D0,BC(N,1) - XPer/2.0d0)
          if (abs(VG(N,2) - BC(N,2)) .ge. YPer/2.0d0)    VGg(2) = VG(N,2) + YPer*sign(1.D0,BC(N,2) - yPer/2.0d0)
	   if (abs(VG(N,3) - BC(N,3)) .ge. zPer/2.0d0)    VGg(3) = VG(N,3) + zPer*sign(1.D0,BC(N,3) - zPer/2.0d0)
         
        END IF
          XCC(1:3) =  matmul(vva1(1:3,1:3),VGg(1:3) -VEXT(1,1:3))
          IF ((XcC(1).GE.ZERO).AND.(XcC(2).GE.ZERO).AND.(XcC(3).GE.ZERO))THEN
	    ISAT4=1
	  ELSE
	      ISAT4=0
	  END IF
          
          
	 IF ((ISAT1.EQ.1).or.(ISAT2.EQ.1).OR.(ISAT3.EQ.1).OR.(ISAT4.EQ.1))THEN
	    ISATISFIED(N)=1
	  ELSE
	      ISATISFIED(N)=0
	  END IF
           


CASE (3)
vext(1,:)=BC(N,:)
	vgg(:)=vg(n,:)
	CALL COMPUTEJACOBIANS
	 IF (IPERIODICITY.EQ.1) THEN
	 if (abs(VG(N,1) - BC(N,1)) .ge. XPer/2.d0)    VGg(1) = VG(N,1) + XPer*sign(1.D0,BC(N,1) - XPer/2.d0)
          if (abs(VG(N,2) - BC(N,2)) .ge. YPer/2.d0)    VGg(2) = VG(N,2) + YPer*sign(1.D0,BC(N,2) - yPer/2.d0)
	   if (abs(VG(N,3) - BC(N,3)) .ge. zPer/2.d0)    VGg(3) = VG(N,3) + zPer*sign(1.D0,BC(N,3) - zPer/2.d0)
         
        END IF
          XCC =  matmul(vva1(:,:),VGg(:) -VEXT(1,:))
          IF ((XcC(1).GE.ZERO).AND.(XcC(2).GE.ZERO).AND.(XcC(3).GE.ZERO))THEN
	    ISAT1=1
	  ELSE
	      ISAT1=0
	  END IF

 IF ((ISAT1.EQ.1))THEN
	    ISATISFIED(N)=1
	  ELSE
	      ISATISFIED(N)=0
	  END IF


case(2)


       vext(1,1:2)=BC(N,1:2)
	vgg(1:2)=vg(n,1:2)
	CALL COMPUTeJACOBIANS2
	 IF (IPERIODICITY.EQ.1) THEN
	 if (abs(VG(N,1) - BC(N,1)) .ge. XPer/2.0d0)    VGg(1) = VG(N,1) + XPer*sign(1.0D0,BC(N,1) - XPer/2.0d0)
          if (abs(VG(N,2) - BC(N,2)) .ge. YPer/2.0d0)    VGg(2) = VG(N,2) + YPer*sign(1.0D0,BC(N,2) - yPer/2.0d0)
	 
         
        END IF
          XCC(1:2) =  matmul(vva1(1:2,1:2),VGg(1:2) -VEXT(1,1:2))
          IF ((XcC(1).GE.zero).AND.(XcC(2).GE.zero))THEN
	    ISAT1=1
	  ELSE
	      ISAT1=0
	  END IF

 IF ((ISAT1.EQ.1))THEN
	    ISATISFIED(N)=1
	  ELSE
	      ISATISFIED(N)=0
	  END IF


end select





! 
 END SUBROUTINE CHECK_CONDITION





subroutine sortstencils(N)
!> @brief
!> This subroutine sorts the stencils with respect to their distance from the cell-centre
implicit none
INTEGER,INTENT(IN)::N
real,dimension(iselemt(n))::rdistl
integer,dimension(iselemt(n))::idistl
real::minds,xc,yc,zc,xl,yl,zl
integer::i,j,k,l,m,kmaxe
KMAXE=XMPIELRANK(N)
DO I=1,KMAXE
	minds=0.0
	  xc=ielem(n,i)%xxc
	  yc=ielem(n,i)%yyc
	  zc=ielem(n,i)%zzc
	DO J=2,iselemt(n)

	    IF (XMPIE(ILOCALALLELG(N,I,1,J)).ne.N)THEN
	      xl=CENTERR((ILOCALALLELG(N,I,1,J)),1)
	      yl=CENTERR((ILOCALALLELG(N,I,1,J)),2)
	      zl=CENTERR((ILOCALALLELG(N,I,1,J)),3)
	    end if
	    IF (XMPIE(ILOCALALLELG(N,I,1,J)).eq.N)THEN
	      xl=ielem(n,XMPIL((ILOCALALLELG(N,I,1,J))))%xxc
	      yl=ielem(n,XMPIL((ILOCALALLELG(N,I,1,J))))%yyc
	      zl=ielem(n,XMPIL((ILOCALALLELG(N,I,1,J))))%zzc
	    end if
	      rdistl(j)=SQRT(((XL-XC)**2.0)+((YL-YC)**2.0)+((ZL-ZC)**2.0))
		
	END DO
	do k=2,iselemt(n)
	do j=2,iselemt(n)-1
	if (rdistl(j)>rdistl(j+1))then
	  m=ILOCALALLELG(N,I,1,J)
	  ILOCALALLELG(N,I,1,J)=ILOCALALLELG(N,I,1,J+1)
	  ILOCALALLELG(N,I,1,J+1)=m
	end if
   enddo 
  enddo  

end do








end subroutine sortstencils




SUBROUTINE STENCIILS_EESX(N,IELEM,ILOCALALLELG,TYPESTEN,ILOCALSTENCIL,NUMNEIGHBOURS,ISELEMT,XMPIE,&
XMPIELRANK,ISIZE,BC,VC,VG,ISATISFIED,IWHICHSTEN,IPERIODICITY,XPER,YPER,ZPER,ISSF,ISHYAPE,XMPIL,CENTERR)
!> @brief
!> This subroutine builds all the directional stencils from the large stencil
IMPLICIT NONE
INTEGER,INTENT(IN)::N,TYPESTEN,NUMNEIGHBOURS,ISIZE,IPERIODICITY,ISSF
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
REAL,INTENT(INOUT)::XPER,YPER,ZPER
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::BC,VG
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::VC
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IWHICHSTEN,ISHYAPE,ISATISFIED
INTEGER::I,J,K,L,M,O,P,KMAXE,ICOUNT,ICPUID,IATRUE,IG,IL,IFG,STNSHA,ITGH,KK,KXK,IX,IFVS,IXCZ,iadd,iadd2,iadd3,IADDX,iaddx1
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTSIZE,OPS
INTEGER::COUNTERST1,ISB1,ISB2,ISB3,TEMPSTEN3,ifno,igvd,CANDID2
INTEGER,ALLOCATABLE,DIMENSION(:)::TOTALN
REAL,ALLOCATABLE,DIMENSION(:)::IFIN,TFIN
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
TYPE::TEMPCOUNT
INTEGER::TOTAL
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTTOT
END TYPE
INTEGER,ALLOCATABLE,DIMENSION(:)::SENDTW,RECVFW
REAL,ALLOCATABLE,DIMENSION(:,:)::SENDVER,RECVER
TYPE::WENOSTENC
INTEGER,ALLOCATABLE,DIMENSION(:,:)::IDCN
REAL,ALLOCATABLE,DIMENSION(:,:,:)::VERTICES
END TYPE
! TYPE(WENOSTENC),ALLOCATABLE,DIMENSION(:)::IWENOST
TYPE(TEMPCOUNT),ALLOCATABLE,DIMENSION(:)::ICOUNT1
KMAXE=XMPIELRANK(N)

if (dimensiona.eq.3)then
! ALLOCATE(BC(N:N,3))
! ALLOCATE(VC(N:N,8,3))
! ALLOCATE(VG(N:N,3))
! ALLOCATE(IWHICHSTEN(N:N))
! ALLOCATE(ISATISFIED(N:N))
! ALLOCATE(ISHYAPE(N:N))
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP END DO
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
		
		
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		     BC(N,1)=IELEM(N,I)%XXC;	BC(N,2)=IELEM(N,I)%YYC;	BC(N,3)=IELEM(N,I)%ZZC;	
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
		
		 IL=0
			DO IADDX=1,STNSHA	!for all stencils
			
			
			if (ielem(n,i)%types_faces(iaddx).eq.5)then
			
			igvd=2
			
			
			else
			
			igvd=2
			
			
			
			end if
			
			DO IADDX1=1,igvd
			IL=IL+1
			ifno=3
			
			if (ielem(n,i)%types_faces(iaddx).eq.5)then
			
			
			IF (iADDX1.EQ.1)THEN
			
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
				
						
			ELSE
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,4))%cord(1:3)
			END IF
		
		
			else
			
			IF (iADDX1.EQ.1)THEN
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3)
			vext(4,1:3)=(inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)+inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3))/2.0D0
			
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(4,1:3))/2.0d0
			
			
			
			end if
			IF (iADDX1.EQ.2)THEN
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=(inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)+inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3))/2.0D0
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
			
			
			
			end if
			
! 			IF (iADDX1.EQ.3)THEN
! 			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
! 			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
! 			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3)
! 			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
! 			
! 			
! 			
! 			end if
			
			
			
			
			
			
			end if
					




			IF ((IELEM(N,I)%INEIGHG(iaddx).GT.0))THEN
			IWHICHSTEN(N)=IL
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN
  !DO IFG=1,KMAXE
	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	   VG(N,1)=IELEM(N,ICONSI)%XXC ;VG(N,2)=IELEM(N,I)%YYC;VG(N,3)=IELEM(N,I)%ZZC
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	  
	  CANDID2=ILOCALALLELG(N,I,1,J)
	  IF (CAND(XMPIE(CANDID2)).GT.0)THEN
                 
                 VG(N,1:3)=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),1:3)
	  
	  
! 	    VG(N,1)=CENTERR((ILOCALALLELG(N,I,1,J)),1)
! 	  VG(N,2)=CENTERR((ILOCALALLELG(N,I,1,J)),2)
! 	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

	    CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
	    IF (ISATISFIED(N).EQ.1)THEN
	    if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
	    ITGH=ITGH+1
	    ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
	    END IF
	    end if
	    IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
	    EXIT
	    END IF
	    
	    
	    ELSE
	     WRITE(99+N,*)"KILLED HERE 3 DUE TO STENCILS REQUIRING MORE CPUS"
							    STOP
	    
	    
	    END IF
  
  END IF
END IF
					END DO		!elements in stencil
			      end if
			END DO		!directional stencils
			end do
	
	END DO
	!$OMP END DO
END IF



end if

if (dimensiona.eq.2)then

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 	
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP end do
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
		
		
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		    iconsi=i
		 BC(N,1)=IELEM(N,I)%XXC;	BC(N,2)=IELEM(N,I)%YYC;
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
		  IL=0
			DO IADDX=1,STNSHA	!for all stencils
			
			DO IADDX1=1,2
			IL=IL+1
			ifno=2
			
			IF (iADDX1.EQ.1)THEN
			vext(2,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:2)
			vext(3,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:2)
			VEXT(3,1:2)=(VEXT(2,1:2)+VEXT(3,1:2))*OO2
			
			ELSE
			vext(2,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:2)
			vext(3,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:2)
			VEXT(2,1:2)=(VEXT(2,1:2)+VEXT(3,1:2))*OO2
			END IF
			
			
			
			IF ((IELEM(N,I)%INEIGHG(IADDX).GT.0))THEN
			IWHICHSTEN(N)=IL
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN

	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	    VG(N,1)=IELEM(N,ICONSI)%XXC ;VG(N,2)=IELEM(N,I)%YYC
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	 CANDID2=ILOCALALLELG(N,I,1,J)
	  IF (CAND(XMPIE(CANDID2)).GT.0)THEN
                 
                 VG(N,1:2)=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),1:2)
! 	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

	      CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
	      IF (ISATISFIED(N).EQ.1)THEN
	      if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
	      ITGH=ITGH+1
	      ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
	      END IF
	      end if
	      IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
	      EXIT
	      END IF
	      ELSE
	      
	      WRITE(99+N,*)"KILLED HERE 2 DUE TO STENCILS REQUIRING MORE CPUS"
							    STOP
	      
	      END IF
  
  END IF
END IF
					END DO		!elements in stencil
			      end if
			      
			END DO		!directional stencils
	
	
	
	END DO
	END DO
	!$OMP END DO
END IF


















end if






END SUBROUTINE STENCIILS_EESX


SUBROUTINE STENCIILS_EES(N,IELEM,ILOCALALLELG,TYPESTEN,ILOCALSTENCIL,NUMNEIGHBOURS,ISELEMT,XMPIE,&
XMPIELRANK,ISIZE,BC,VC,VG,ISATISFIED,IWHICHSTEN,IPERIODICITY,XPER,YPER,ZPER,ISSF,ISHYAPE,XMPIL,CENTERR)
!> @brief
!> This subroutine builds all the directional stencils from the large stencil (suitable for periodic boundaries)
IMPLICIT NONE
INTEGER,INTENT(IN)::N,TYPESTEN,NUMNEIGHBOURS,ISIZE,IPERIODICITY,ISSF
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
REAL,INTENT(INOUT)::XPER,YPER,ZPER
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::BC,VG
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::VC
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IWHICHSTEN,ISHYAPE,ISATISFIED
INTEGER::I,J,K,L,M,O,P,KMAXE,ICOUNT,ICPUID,IATRUE,IG,IL,IFG,STNSHA,ITGH,KK,KXK,IX,IFVS,IXCZ,iadd,iadd2,iadd3,IADDX,iaddx1
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTSIZE,OPS
INTEGER::COUNTERST1,ISB1,ISB2,ISB3,TEMPSTEN3,ifno,igvd
INTEGER,ALLOCATABLE,DIMENSION(:)::TOTALN
REAL,ALLOCATABLE,DIMENSION(:)::IFIN,TFIN
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
TYPE::TEMPCOUNT
INTEGER::TOTAL
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTTOT
END TYPE
INTEGER,ALLOCATABLE,DIMENSION(:)::SENDTW,RECVFW
REAL,ALLOCATABLE,DIMENSION(:,:)::SENDVER,RECVER
TYPE::WENOSTENC
INTEGER,ALLOCATABLE,DIMENSION(:,:)::IDCN
REAL,ALLOCATABLE,DIMENSION(:,:,:)::VERTICES
END TYPE
! TYPE(WENOSTENC),ALLOCATABLE,DIMENSION(:)::IWENOST
TYPE(TEMPCOUNT),ALLOCATABLE,DIMENSION(:)::ICOUNT1
KMAXE=XMPIELRANK(N)

if (dimensiona.eq.3)then
! ALLOCATE(BC(N:N,3))
! ALLOCATE(VC(N:N,8,3))
! ALLOCATE(VG(N:N,3))
! ALLOCATE(IWHICHSTEN(N:N))
! ALLOCATE(ISATISFIED(N:N))
! ALLOCATE(ISHYAPE(N:N))
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 	
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP END DO
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
		
		
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		      call COMPUTE_CENTRE3d(N,Iconsi)
			
		BC(N,1)=cords(1)   ;   BC(N,2)=cords(2);    BC(N,3)=cords(3)
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
		
		 IL=0
			DO IADDX=1,STNSHA	!for all stencils
			
			
			if (ielem(n,i)%types_faces(iaddx).eq.5)then
			
			igvd=2
			
			
			else
			
			igvd=2
			
			
			
			end if
			
			DO IADDX1=1,igvd
			IL=IL+1
			ifno=3
			
			if (ielem(n,i)%types_faces(iaddx).eq.5)then
			
			
			IF (iADDX1.EQ.1)THEN
			
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
				
						
			ELSE
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,4))%cord(1:3)
			END IF
		
		
			else
			
			IF (iADDX1.EQ.1)THEN
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3)
			vext(4,1:3)=(inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)+inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3))/2.0D0
			
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(4,1:3))/2.0d0
			
			
			
			end if
			IF (iADDX1.EQ.2)THEN
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
			vext(3,1:3)=(inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)+inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3))/2.0D0
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
			
			
			
			end if
			
! 			IF (iADDX1.EQ.3)THEN
! 			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,3))%cord(1:3)
! 			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:3)
! 			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:3)
! 			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
! 			
! 			
! 			
! 			end if
			
			
			
			
			
			
			end if
					

! 			if (ifno.eq.3)then
! 			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
! 			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
! 			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
! 			else
! 			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
! 			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
! 			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
! 			vext(5,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
! 			vext(6,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
! 			vext(7,1:3)=inoder(ielem(n,i)%nodes_faces(il,4))%cord(1:3)
! 
! 			end if


			IF ((IELEM(N,I)%INEIGHG(iaddx).GT.0))THEN
			IWHICHSTEN(N)=IL
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN
  !DO IFG=1,KMAXE
	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	  call COMPUTE_CENTRE3d(N,Iconsi)
	  VG(N,1)=cords(1) ;VG(N,2)=cords(2) ; VG(N,3)=cords(3)
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	    VG(N,1)=CENTERR((ILOCALALLELG(N,I,1,J)),1)
	  VG(N,2)=CENTERR((ILOCALALLELG(N,I,1,J)),2)
	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
  ITGH=ITGH+1
  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
  EXIT
  END IF
  
  END IF
END IF
					END DO		!elements in stencil
			      end if
			END DO		!directional stencils
			end do
	
	END DO
	!$OMP END DO
END IF



end if

if (dimensiona.eq.2)then

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE

	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP end do
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
		
		
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		      call COMPUTE_CENTRE2d(N,Iconsi)
			
		BC(N,1)=cords(1)   ;   BC(N,2)=cords(2);    !BC(N,3)=cords(3)
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
		  IL=0
			DO IADDX=1,STNSHA	!for all stencils
			
			DO IADDX1=1,2
			IL=IL+1
			ifno=2
			
			IF (iADDX1.EQ.1)THEN
			vext(2,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:2)
			vext(3,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:2)
			VEXT(3,1:2)=(VEXT(2,1:2)+VEXT(3,1:2))*OO2
			
			ELSE
			vext(2,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,1))%cord(1:2)
			vext(3,1:2)=inoder(ielem(n,i)%nodes_faces(IADDX,2))%cord(1:2)
			VEXT(2,1:2)=(VEXT(2,1:2)+VEXT(3,1:2))*OO2
			END IF
			
			
			
			IF ((IELEM(N,I)%INEIGHG(IADDX).GT.0))THEN
			IWHICHSTEN(N)=IL
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN

	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	  call COMPUTE_CENTRE2d(N,Iconsi)
	  VG(N,1)=cords(1) ;VG(N,2)=cords(2) 
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	    VG(N,1)=CENTERR((ILOCALALLELG(N,I,1,J)),1)
	  VG(N,2)=CENTERR((ILOCALALLELG(N,I,1,J)),2)
! 	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ielem(n,i)%iNUMNEIGHBOURS)then
  ITGH=ITGH+1
  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ielem(n,i)%iNUMNEIGHBOURS)THEN
  EXIT
  END IF
  
  END IF
END IF
					END DO		!elements in stencil
			      end if
			      
			END DO		!directional stencils
	
	
	
	END DO
	END DO
	!$OMP END DO
END IF


















end if






END SUBROUTINE STENCIILS_EES





SUBROUTINE STENCIILSX(N,IELEM,ILOCALALLELG,TYPESTEN,ILOCALSTENCIL,NUMNEIGHBOURS,ISELEMT,XMPIE,&
XMPIELRANK,ISIZE,BC,VC,VG,ISATISFIED,IWHICHSTEN,IPERIODICITY,XPER,YPER,ZPER,ISSF,ISHYAPE,XMPIL,CENTERR)
!> @brief
!> This subroutine builds all the directional stencils from the large stencil based on various aglorithms
IMPLICIT NONE
INTEGER,INTENT(IN)::N,TYPESTEN,NUMNEIGHBOURS,ISIZE,IPERIODICITY,ISSF
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
REAL,INTENT(INOUT)::XPER,YPER,ZPER
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::BC,VG
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::VC
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IWHICHSTEN,ISHYAPE,ISATISFIED
INTEGER::I,J,K,L,M,O,P,KMAXE,ICOUNT,ICPUID,IATRUE,IG,IL,IFG,STNSHA,ITGH,KK,KXK,IX,IFVS,IXCZ,iadd,iadd2,iadd3,ITARGET
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTSIZE,OPS
INTEGER::COUNTERST1,ISB1,ISB2,ISB3,TEMPSTEN3,ifno,INV,CANDID2
INTEGER,ALLOCATABLE,DIMENSION(:)::TOTALN
REAL,ALLOCATABLE,DIMENSION(:)::IFIN,TFIN
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
TYPE::TEMPCOUNT
INTEGER::TOTAL
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTTOT
END TYPE
INTEGER,ALLOCATABLE,DIMENSION(:)::SENDTW,RECVFW
REAL,ALLOCATABLE,DIMENSION(:,:)::SENDVER,RECVER
TYPE::WENOSTENC
INTEGER,ALLOCATABLE,DIMENSION(:,:)::IDCN
REAL,ALLOCATABLE,DIMENSION(:,:,:)::VERTICES
END TYPE
! TYPE(WENOSTENC),ALLOCATABLE,DIMENSION(:)::IWENOST
TYPE(TEMPCOUNT),ALLOCATABLE,DIMENSION(:)::ICOUNT1
KMAXE=XMPIELRANK(N)

if (dimensiona.eq.3)then

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 	
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP END DO
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
                 IF (EES.EQ.5)THEN
                   ITARGET=NUMNEIGHBOURS2
                ELSE 
                   ITARGET=ielem(n,i)%iNUMNEIGHBOURS
                END IF
                 
                    
		
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
! 		      call COMPUTE_CENTRE3d(N,Iconsi)
		BC(N,1)=IELEM(N,I)%XXC;	BC(N,2)=IELEM(N,I)%YYC;	BC(N,3)=IELEM(N,I)%ZZC;	
! 		BC(N,1)=cords(1)   ;   BC(N,2)=cords(2);    BC(N,3)=cords(3)
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
			DO IL=1,STNSHA	!for all stencils
			   if (ielem(n,i)%types_faces(il).eq.5)then
			ifno=4
			else
			ifno=3
			end if
			

			if (ifno.eq.3)then
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
			else
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
			vext(5,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
			vext(6,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
			vext(7,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
                        vext(8,1:3)=inoder(ielem(n,i)%nodes_faces(il,4))%cord(1:3)
			end if


			IF ((IELEM(N,I)%INEIGHG(il).GT.0))THEN
			IWHICHSTEN(N)=IL
			
! 			
			
			
			
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						
						
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN
  !DO IFG=1,KMAXE
	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	  VG(N,1)=IELEM(N,Iconsi)%XXC;	VG(N,2)=IELEM(N,Iconsi)%YYC;	VG(N,3)=IELEM(N,Iconsi)%ZZC;	
	  ! 	  call COMPUTE_CENTRE3d(N,Iconsi)
! 	  VG(N,1)=cords(1) ;VG(N,2)=cords(2) ; VG(N,3)=cords(3)
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ITARGET)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ITARGET)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	  CANDID2=ILOCALALLELG(N,I,1,J)
		IF (CAND(XMPIE(CANDID2)).GT.0)THEN
                 
                 VG(N,1:3)=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),1:3)
                

		CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
		      IF (ISATISFIED(N).EQ.1)THEN
		      if (itgh+1.le.ITARGET)then
		      ITGH=ITGH+1
		      ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
		      END IF
		      end if
		      IF (ITGH.EQ.ITARGET)THEN
		      EXIT
		      END IF
		ELSE
		
		WRITE(99+N,*)"KILLED HERE2 DUE TO STENCILS REQUIRING MORE CPUS"
							    STOP
		
		
		END IF
		END IF
END IF
                        
					END DO		!elements in stencil
			      end if
			END DO		!directional stencils
	
	END DO
	!$OMP END DO
END IF



end if

if (dimensiona.eq.2)then

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 	
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP END DO
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
		
                         IF (EES.EQ.5)THEN
                   ITARGET=NUMNEIGHBOURS2
                ELSE 
                   ITARGET=ielem(n,i)%iNUMNEIGHBOURS
                END IF
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		 BC(N,1)=IELEM(N,I)%XXC;	BC(N,2)=IELEM(N,I)%YYC;
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
			DO IL=1,STNSHA	!for all stencils
			ifno=2

			vext(2,1:2)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:2)
			vext(3,1:2)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:2)

! 			


			IF ((IELEM(N,I)%INEIGHG(il).GT.0))THEN
			IWHICHSTEN(N)=IL
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
			
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN

	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	   VG(N,1)=IELEM(N,ICONSI)%XXC ;VG(N,2)=IELEM(N,Iconsi)%YYC
! 	  call COMPUTE_CENTRE2d(N,Iconsi)
! 	  VG(N,1)=cords(1) ;VG(N,2)=cords(2) 
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ITARGET)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ITARGET)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	  CANDID2=ILOCALALLELG(N,I,1,J)
	  IF (CAND(XMPIE(CANDID2)).GT.0)THEN
                 
                 VG(N,1:2)=XAND2R(CAND(XMPIE(CANDID2)),XMPIL(CANDID2),1:2)
	  
! 	    VG(N,1)=CENTERR((ILOCALALLELG(N,I,1,J)),1)
! 	  VG(N,2)=CENTERR((ILOCALALLELG(N,I,1,J)),2)
! 	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

	  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
	  IF (ISATISFIED(N).EQ.1)THEN
	  if (itgh+1.le.ITARGET)then
	  ITGH=ITGH+1
	  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
	  END IF
	  end if
	  IF (ITGH.EQ.ITARGET)THEN
	  EXIT
	  END IF
  
	  ELSE
	  
	  WRITE(99+N,*)"KILLED HERE 2 DUE TO STENCILS REQUIRING MORE CPUS"
							    STOP
	  
	  END IF
  END IF
END IF
					END DO		!elements in stencil
			      end if
			END DO		!directional stencils
	
	END DO
	!$OMP END DO
END IF

















end if






END SUBROUTINE STENCIILSX


SUBROUTINE STENCIILS(N,IELEM,ILOCALALLELG,TYPESTEN,ILOCALSTENCIL,NUMNEIGHBOURS,ISELEMT,XMPIE,&
XMPIELRANK,ISIZE,BC,VC,VG,ISATISFIED,IWHICHSTEN,IPERIODICITY,XPER,YPER,ZPER,ISSF,ISHYAPE,XMPIL,CENTERR)
!> @brief
!> This subroutine builds all the directional stencils from the large stencil based on various algorithms (suitable for period boundaries)
IMPLICIT NONE
INTEGER,INTENT(IN)::N,TYPESTEN,NUMNEIGHBOURS,ISIZE,IPERIODICITY,ISSF
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALALLELG	
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
REAL,INTENT(INOUT)::XPER,YPER,ZPER
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::BC,VG
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::VC
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::ISELEMT
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIL
INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IWHICHSTEN,ISHYAPE,ISATISFIED
INTEGER::I,J,K,L,M,O,P,KMAXE,ICOUNT,ICPUID,IATRUE,IG,IL,IFG,STNSHA,ITGH,KK,KXK,IX,IFVS,IXCZ,iadd,iadd2,iadd3,ITARGET
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTSIZE,OPS
INTEGER::COUNTERST1,ISB1,ISB2,ISB3,TEMPSTEN3,ifno,INV
INTEGER,ALLOCATABLE,DIMENSION(:)::TOTALN
REAL,ALLOCATABLE,DIMENSION(:)::IFIN,TFIN
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::CENTERR
TYPE::TEMPCOUNT
INTEGER::TOTAL
INTEGER,ALLOCATABLE,DIMENSION(:)::COUNTTOT
END TYPE
INTEGER,ALLOCATABLE,DIMENSION(:)::SENDTW,RECVFW
REAL,ALLOCATABLE,DIMENSION(:,:)::SENDVER,RECVER
TYPE::WENOSTENC
INTEGER,ALLOCATABLE,DIMENSION(:,:)::IDCN
REAL,ALLOCATABLE,DIMENSION(:,:,:)::VERTICES
END TYPE
! TYPE(WENOSTENC),ALLOCATABLE,DIMENSION(:)::IWENOST
TYPE(TEMPCOUNT),ALLOCATABLE,DIMENSION(:)::ICOUNT1
KMAXE=XMPIELRANK(N)

if (dimensiona.eq.3)then

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 	
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP END DO
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
                 IF (EES.EQ.5)THEN
                   ITARGET=NUMNEIGHBOURS2
                ELSE 
                   ITARGET=ielem(n,i)%iNUMNEIGHBOURS
                END IF
                 
                    
		
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		      call COMPUTE_CENTRE3d(N,Iconsi)
			
		BC(N,1)=cords(1)   ;   BC(N,2)=cords(2);    BC(N,3)=cords(3)
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
			DO IL=1,STNSHA	!for all stencils
			   if (ielem(n,i)%types_faces(il).eq.5)then
			ifno=4
			else
			ifno=3
			end if
			

			if (ifno.eq.3)then
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
			else
			vext(2,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
			vext(3,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
			vext(4,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
			vext(5,1:3)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:3)
			vext(6,1:3)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:3)
			vext(7,1:3)=inoder(ielem(n,i)%nodes_faces(il,3))%cord(1:3)
                        vext(8,1:3)=inoder(ielem(n,i)%nodes_faces(il,4))%cord(1:3)
			end if


			IF ((IELEM(N,I)%INEIGHG(il).GT.0))THEN
			IWHICHSTEN(N)=IL
			
! 	
			
			
			
			
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						
						
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN
  !DO IFG=1,KMAXE
	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	  call COMPUTE_CENTRE3d(N,Iconsi)
	  VG(N,1)=cords(1) ;VG(N,2)=cords(2) ; VG(N,3)=cords(3)
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ITARGET)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ITARGET)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	    VG(N,1)=CENTERR((ILOCALALLELG(N,I,1,J)),1)
	  VG(N,2)=CENTERR((ILOCALALLELG(N,I,1,J)),2)
	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ITARGET)then
  ITGH=ITGH+1
  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ITARGET)THEN
  EXIT
  END IF
  
  END IF
END IF
                        
					END DO		!elements in stencil
			      end if
			END DO		!directional stencils
	
	END DO
	!$OMP END DO
END IF



end if

if (dimensiona.eq.2)then

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
	
!-------------------FOR DEBUGGING ONLY -----------------------------------------!
!$OMP DO SCHEDULE (STATIC)
DO I=1,KMAXE
! 	
	DO J=1,ielem(n,i)%iNUMNEIGHBOURS
	      
		ILOCALSTENCIL(N,I,1,J)=ILOCALALLELG(N,I,1,J)
! 		
	END DO
	
END DO
!$OMP END DO
IF (TYPESTEN.GT.1)THEN
	ICPUID=N
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE	!for all elements
		
                         IF (EES.EQ.5)THEN
                   ITARGET=NUMNEIGHBOURS2
                ELSE 
                   ITARGET=ielem(n,i)%iNUMNEIGHBOURS
                END IF
			STNSHA=ielem(n,i)%ifca
			    iconsi=i
		      call COMPUTE_CENTRE2d(N,Iconsi)
			
		BC(N,1)=cords(1)   ;   BC(N,2)=cords(2);    !BC(N,3)=cords(3)
	    
		ISHYAPE(N)=IELEM(N,I)%ISHAPE
			DO IL=1,STNSHA	!for all stencils
			ifno=2

			vext(2,1:2)=inoder(ielem(n,i)%nodes_faces(il,1))%cord(1:2)
			vext(3,1:2)=inoder(ielem(n,i)%nodes_faces(il,2))%cord(1:2)




			IF ((IELEM(N,I)%INEIGHG(il).GT.0))THEN
			IWHICHSTEN(N)=IL
			
			ILOCALSTENCIL(N,I,IL+1,1)=ILOCALALLELG(N,I,1,1)
			
					ITGH=1
					DO J=2,ISELEMT(N)!for all stencil elements
						IF ((ILOCALALLELG(N,I,1,J)).GT.0) THEN
						ISATISFIED(N)=0
IF (XMPIE(ILOCALALLELG(N,I,1,J)).EQ.N)THEN

	  IFG=XMPIL((ILOCALALLELG(N,I,1,J)))
	  iconsi=ifg
	  call COMPUTE_CENTRE2d(N,Iconsi)
	  VG(N,1)=cords(1) ;VG(N,2)=cords(2) 
		     n_node=ifno

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ITARGET)then
  ITGH=ITGH+1

  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ITARGET)THEN
  EXIT
  END IF
  
  END IF
  IF (XMPIE(ILOCALALLELG(N,I,1,J)).NE.N)THEN
	  n_node=ifno
	    VG(N,1)=CENTERR((ILOCALALLELG(N,I,1,J)),1)
	  VG(N,2)=CENTERR((ILOCALALLELG(N,I,1,J)),2)
! 	  VG(N,3)=CENTERR((ILOCALALLELG(N,I,1,J)),3)

  CALL CHECK_CONDITION(N,IWHICHSTEN,ISATISFIED,IPERIODICITY,XPER,YPER,ZPER,ISHYAPE,ISSF,BC,VC,VG)
  IF (ISATISFIED(N).EQ.1)THEN
  if (itgh+1.le.ITARGET)then
  ITGH=ITGH+1
  ILOCALSTENCIL(N,I,IL+1,ITGH)=ILOCALALLELG(N,I,1,J)
  END IF
  end if
  IF (ITGH.EQ.ITARGET)THEN
  EXIT
  END IF
  
  END IF
END IF
					END DO		!elements in stencil
			      end if
			END DO		!directional stencils
	
	END DO
	!$OMP END DO
END IF

















end if






END SUBROUTINE STENCIILS



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SUBROUTINE FOR DETERMINING THE NUMBER OF ELEMENTS!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!REQUIRED FOR EACH STENCIL!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!FOR VARIOUS ORDER OF ACCURACY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE DETERMINE_SIZE(N,IORDER,ISELEM,ISELEMT,IOVERST,IOVERTO,ILX,NUMNEIGHBOURS,IDEGFREE,IMAXDEGFREE,IEXTEND)
!> @brief
!> This subroutine determines the degress of freedom, neighbours and polynomial order for each stencil of each cell
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::IORDER
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ISELEMT
	INTEGER,INTENT(OUT)::IOVERST,ISELEM
	INTEGER,INTENT(OUT)::IOVERTO
	INTEGER,INTENT(OUT)::IDEGFREE
	INTEGER,INTENT(OUT)::IMAXDEGFREE
	INTEGER,INTENT(OUT)::ILX     !NUMBER OF DEGREES OF FREEDOM
	INTEGER,INTENT(INOUT)::NUMNEIGHBOURS
	INTEGER,INTENT(IN)::IEXTEND,N
	integer::i,itemd
	ALLOCATE(ISELEMT(N:N))
	if (dimensiona.eq.3)then
	
	IF (IORDER.GE.2)THEN
		ILX=((IORDER+1)*(IORDER+2)*(IORDER+3))/6
	IDEGFREE=ILX-1
	NUMNEIGHBOURS=ILX*extf
	IMAXDEGFREE=NUMNEIGHBOURS-1
	ISELEM=(ILX*extf)*IEXTEND
	ISELEMT(N:N)=ISELEM
	IOVERST=ISELEM
	IOVERTO=ISELEM
	
	SELECT CASE(IORDER)
	
	CASE (1,2,3)
	idegfree2=3
	IORDER2=1
	NUMNEIGHBOURS2=9
	
	
	CASE(4,5,6,7)
	idegfree2=3
	IORDER2=1
	NUMNEIGHBOURS2=9
	
	
	END SELECT
	
	
      
      
	END IF
	IF (IORDER.EQ.1)THEN
		ILX=((IORDER+1)*(IORDER+2)*(IORDER+3))/6
	IDEGFREE=ILX-1
	NUMNEIGHBOURS=ILX*extf
	IMAXDEGFREE=NUMNEIGHBOURS-1
	ISELEM=(ILX*extf)*IEXTEND
	ISELEMT(N:N)=ISELEM
	IOVERST=ISELEM
	IOVERTO=ISELEM
	
	idegfree2=3
	IORDER2=1
	NUMNEIGHBOURS2=9
	
	
	
	END IF
	else
	
	
	
	IF (IORDER.GE.2)THEN
		ILX=((IORDER+1)*(IORDER+2))/2
	IDEGFREE=ILX-1
	NUMNEIGHBOURS=ILX*extf
	IMAXDEGFREE=NUMNEIGHBOURS-1
	itemd=(ILX*extf)*IEXTEND
	ISELEM=min(itemd,imaxe-1)
	ISELEMT(N:N)=ISELEM
	IOVERST=ISELEM
	IOVERTO=ISELEM
	SELECT CASE(IORDER)
	
	CASE (1,2,3)
	idegfree2=2
	IORDER2=1
	NUMNEIGHBOURS2=5
	
	
	CASE(4,5,6,7)
	idegfree2=2
	IORDER2=1
	NUMNEIGHBOURS2=5
	
	
	END SELECT
	
	END IF
	IF (IORDER.EQ.1)THEN
		ILX=((IORDER+1)*(IORDER+2))/2
	IDEGFREE=ILX-1
	NUMNEIGHBOURS=ILX*extf
	IMAXDEGFREE=NUMNEIGHBOURS-1
	ISELEM=(ILX*extf)*IEXTEND
	ISELEMT(N:N)=ISELEM
	IOVERST=ISELEM
	IOVERTO=ISELEM
	IORDER2=1
	idegfree2=2
	NUMNEIGHBOURS2=5
	
	END IF





	end if
	do i=1,xmpielrank(n)
	  ielem(n,i)%inumneighbours=NUMNEIGHBOURS
	 ielem(n,i)%idegfree=idegfree
	 ielem(n,i)%iorder=iorder
        end do



	END SUBROUTINE DETERMINE_SIZE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!SUBROUTINE CALLED FOR READING DAT FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!FOR THE SIMULATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!WITH INFORMATION SUCH AS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!ORDER OF ACCURACY METHODS ETC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!	
SUBROUTINE GAUSSIANPOINTS(IGQRULES,NUMBEROFPOINTS,NUMBEROFPOINTS2)
!> @brief
!> This subroutine determines quadrature points required for each spatial order of accuracy and for each cell type
IMPLICIT NONE
INTEGER,INTENT(INOUT)::IGQRULES,NUMBEROFPOINTS,NUMBEROFPOINTS2
if (dimensiona.eq.3)then
IF (IGQRULES.EQ.1)THEN
	QP_HEXA=1;QP_TETRA=1;QP_PYRA=1;QP_PRISM=1;
	QP_QUAD=1;QP_TRIANGLE=1
END IF
IF (IGQRULES.EQ.2)THEN

	QP_HEXA=8;QP_TETRA=4;QP_PYRA=5;QP_PRISM=6;
	QP_QUAD=4;QP_TRIANGLE=3
	
	
END IF
IF (IGQRULES.EQ.3)THEN
	QP_HEXA=27;QP_TETRA=10;QP_PYRA=15;QP_PRISM=18;
	QP_QUAD=9;QP_TRIANGLE=6
	
END IF
IF (IGQRULES.EQ.4)THEN
	QP_HEXA=64;QP_TETRA=20;QP_PYRA=15;QP_PRISM=24;
	QP_QUAD=16;QP_TRIANGLE=10
	
END IF
IF (IGQRULES.EQ.5)THEN
	QP_HEXA=125;QP_TETRA=35;QP_PYRA=15;QP_PRISM=50;
	QP_QUAD=25;QP_TRIANGLE=15
	
END IF
IF (IGQRULES.EQ.6)THEN
	QP_HEXA=216;QP_TETRA=35;QP_PYRA=15;QP_PRISM=60;
	QP_QUAD=36;QP_TRIANGLE=21
	
END IF
IF (IGQRULES.GE.7)THEN
	QP_HEXA=216;QP_TETRA=35;QP_PYRA=15;QP_PRISM=60;
	QP_QUAD=36;QP_TRIANGLE=36
	
END IF
if (reduce_comp.eq.1)then
qp_quad_n=1;QP_TRIANGLE_n=1
else
qp_quad_n=qp_quad;QP_TRIANGLE_n=QP_TRIANGLE
end if

else


IF (IGQRULES.EQ.1)THEN
	QP_QUAD=1;QP_TRIANGLE=1;qp_line=1
	
END IF
IF (IGQRULES.EQ.2)THEN
	QP_QUAD=4;QP_TRIANGLE=3;qp_line=2
END IF
IF (IGQRULES.EQ.3)THEN
	QP_QUAD=9;QP_TRIANGLE=6;qp_line=3
END IF
IF (IGQRULES.EQ.4)THEN
	QP_QUAD=16;QP_TRIANGLE=10;qp_line=4
END IF
IF (IGQRULES.EQ.5)THEN
	QP_QUAD=25;QP_TRIANGLE=15;qp_line=5
END IF

IF (IGQRULES.EQ.6)THEN
	QP_QUAD=36;QP_TRIANGLE=21;qp_line=6
END IF

IF (IGQRULES.GE.7)THEN
	QP_QUAD=36;QP_TRIANGLE=36;qp_line=9
END IF


if (reduce_comp.eq.1)then
qp_line_n=1
else
qp_line_n=qp_line
end if


end if
END SUBROUTINE GAUSSIANPOINTS








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writeout
!> @brief
!> This subroutine is solely for debugging purposes
integer::ijd,ic2,kmaxe,i
real,allocatable,dimension(:)::xx,yy,zz
kmaxe=xmpielrank(n)

if (n.eq.0)then
CALL OPEN_INPUT(N,ITT)
	WRITE(97,*)'TITLE="FINAL SOLUTION"'
 	WRITE(97,*)'VARIABLES="X","Y","Z"'
 	
 	
WRITE(97,*) 'Zone N=',4,',E=',1,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
allocate(xx(imaxn),yy(imaxn),zz(imaxn))
DO I=1,kmaxe
    IF (ielem(n,i)%ishape.eq.6)then
	      DO IJD=1,IELEM(N,I)%NONODES
		WRITE(97,*)inoder(IELEM(N,I)%NODES(IJD))%CORD(1)
	    END DO
	    DO IJD=1,IELEM(N,I)%NONODES
		WRITE(97,*)inoder(IELEM(N,I)%NODES(IJD))%CORD(2)
	    END DO
! 	      DO IJD=1,IELEM(N,I)%NONODES
! 		WRITE(97,*)inoder(IELEM(N,I)%NODES(IJD))%CORD(3)
! 	    END DO
	    
    EXIT

  END IF

END DO

do i=1,kmaxe
if (ielem(n,i)%ishape.eq.6)then

write(97,*)1,2,3,3
EXIT
end if

end do

CALL CLOSE_INPUT(N,ITT)
end if



call mpi_barrier(mpi_comm_world,ierror)
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE STENCILS(N,IELEM,IMAXE,XMPIE,XMPIELRANK,ILOCALSTENCIL,TYPESTEN,NUMNEIGHBOURS,RESTART)
!> @brief
!> This subroutine is establishing which of the stencils are admissible and which cells can use the WENO algorithms
	IMPLICIT NONE
	TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
	INTEGER,INTENT(IN)::N,IMAXE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
	INTEGER,INTENT(IN)::NUMNEIGHBOURS
	INTEGER,INTENT(IN)::TYPESTEN,RESTART
	INTEGER::I,J,JI,K,LM,KMAXN,KK,KMAXE,IAA,KX,L,ITRR,ITRX,ITRY,ITARGET
	INTEGER::IT1,IT2,IT3,IT4,IT5,IT6,IT7,IT8,ITX,INX,ITF,IFG,KMH,JNG,ichg
	integer,dimension(12,1000)::items

	
	KMAXE=XMPIELRANK(N)
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	DO I=1,KMAXE
                 IELEM(N,I)%FULL=0
			KX=0
		DO LM=1,TYPESTEN
                                IF ((LM.EQ.1).OR.(EES.NE.5))THEN
                                ITARGET=ielem(n,i)%iNUMNEIGHBOURS
                                ELSE
                                ITARGET=NUMNEIGHBOURS2
                                END IF
                                
				KK=0
				DO IFG=1,ITARGET
					IF (ILOCALSTENCIL(N,I,LM,IFG).GT.0)THEN
					KK=KK+1
					END IF
				END DO
				IF (KK.EQ.ITARGET)THEN
				KX=KX+1
				END IF

			
		END DO
			IELEM(N,I)%ADMIS=KX
			IF ((EES.EQ.0))THEN
                            
                            if (initcond.eq.101)then
                            
                            IF (IELEM(N,I)%ADMIS.eq.ielem(n,i)%ifca+1)THEN
                            IELEM(N,I)%FULL=1
!                             
                            end if
                            
                            else
                            IF (IELEM(N,I)%ADMIS.Gt.3)THEN
                            IELEM(N,I)%FULL=1
                            END IF
                            
                            
                            
                            end if
			END IF
			
			IF ((EES.GE.4))THEN
                            
                            if (initcond.eq.101)then
                            
                            IF (IELEM(N,I)%ADMIS.eq.ielem(n,i)%ifca+1)THEN
                            IELEM(N,I)%FULL=1
!                             
                            end if
                            
                            else
                            IF (IELEM(N,I)%ADMIS.GE.3)THEN
                            IELEM(N,I)%FULL=1
                            END IF
                            
                            
                            
                            end if
			END IF
			IF (EES.EQ.1)THEN
                            IF (IELEM(N,I)%ADMIS.GT.(ielem(n,i)%ifca+1))THEN
                            IELEM(N,I)%FULL=1
                            END IF
			END IF
			IF (EES.EQ.2)THEN
                            IF (IELEM(N,I)%ADMIS.GT.(ielem(n,i)%ifca+1))THEN
                            IELEM(N,I)%FULL=1
                            END IF
			END IF
			
			
			
                        
			
! 			
	END DO	
	
	IF (EES.EQ.2)THEN
	  DO I=1,KMAXE
            ITRR=0

				
	    
	    IF (IELEM(N,I)%ADMIS.EQ.((ielem(n,i)%ifca*2)+1))THEN
	      ITRR=1    
	    
	    
	    END IF
	    IF (ITRR.EQ.1)THEN
	      IELEM(N,I)%ADMIS=ielem(n,i)%ifca+1
	      ITRX=1;ITRY=1
	      DO IFG=2,ielem(n,i)%iNUMNEIGHBOURS
			IF (MOD(IFG,2).eq.0)THEN
			ITRX=ITRX+1
				IF (IELEM(N,I)%ISHAPE.EQ.5)THEN !QUAD
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,4,ITRX)
				items(4,IFG)=ILOCALSTENCIL(N,I,5,ITRX)
                                END IF
				IF (IELEM(N,I)%ISHAPE.EQ.6)THEN !TRI
				
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,4,ITRX)
! 				
				
				
				END IF
                                
                                IF (IELEM(N,I)%ISHAPE.EQ.1)THEN !HEXA
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,6,ITRX)
				items(4,IFG)=ILOCALSTENCIL(N,I,7,ITRX)
				items(5,IFG)=ILOCALSTENCIL(N,I,10,ITRX)
				items(6,IFG)=ILOCALSTENCIL(N,I,11,ITRX)
                                END IF
                                 IF (IELEM(N,I)%ISHAPE.EQ.2)THEN !TETRA
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,6,ITRX)
				items(4,IFG)=ILOCALSTENCIL(N,I,8,ITRX)
                                END IF
                                 IF (IELEM(N,I)%ISHAPE.EQ.3)THEN !PYRAMIDAL
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,4,ITRX)
				items(4,IFG)=ILOCALSTENCIL(N,I,5,ITRX)
				items(5,IFG)=ILOCALSTENCIL(N,I,7,ITRX)
                                END IF
                                IF (IELEM(N,I)%ISHAPE.EQ.4)THEN !PRISM
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,6,ITRX)
				items(4,IFG)=ILOCALSTENCIL(N,I,7,ITRX)
				items(5,IFG)=ILOCALSTENCIL(N,I,8,ITRX)
                                END IF
			ELSE
			ITRY=ITRY+1
			    IF (IELEM(N,I)%ISHAPE.EQ.5)THEN !QUAD
			    items(1,IFG)=ILOCALSTENCIL(N,I,6,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,7,ITRY)
			    items(3,IFG)=ILOCALSTENCIL(N,I,8,ITRY)
			    items(4,IFG)=ILOCALSTENCIL(N,I,9,ITRY)
 			    END IF
			    IF (IELEM(N,I)%ISHAPE.EQ.6)THEN !TRI
			    
			    items(1,IFG)=ILOCALSTENCIL(N,I,5,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,6,ITRY)
			    items(3,IFG)=ILOCALSTENCIL(N,I,7,ITRY)
			    END IF
			    IF (IELEM(N,I)%ISHAPE.EQ.1)THEN !hexa
			    items(1,IFG)=ILOCALSTENCIL(N,I,4,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,5,ITRY)
			    items(3,IFG)=ILOCALSTENCIL(N,I,8,ITRY)
			    items(4,IFG)=ILOCALSTENCIL(N,I,9,ITRY)
			    items(5,IFG)=ILOCALSTENCIL(N,I,12,ITRY)
			    items(6,IFG)=ILOCALSTENCIL(N,I,13,ITRY)
 			    END IF
 			    IF (IELEM(N,I)%ISHAPE.EQ.2)THEN !tetra
			    items(1,IFG)=ILOCALSTENCIL(N,I,5,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,6,ITRY)
			    items(3,IFG)=ILOCALSTENCIL(N,I,9,ITRY)
			    items(4,IFG)=ILOCALSTENCIL(N,I,7,ITRY)
 			    END IF
 			     IF (IELEM(N,I)%ISHAPE.EQ.3)THEN !pyra
			    items(1,IFG)=ILOCALSTENCIL(N,I,10,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,6,ITRY)
			    items(3,IFG)=ILOCALSTENCIL(N,I,8,ITRY)
			    items(4,IFG)=ILOCALSTENCIL(N,I,9,ITRY)
			    items(5,IFG)=ILOCALSTENCIL(N,I,11,ITRY)
			    
 			    END IF
 			    IF (IELEM(N,I)%ISHAPE.EQ.4)THEN !PRISM
				items(1,IFG)=ILOCALSTENCIL(N,I,4,ITRy)
				items(2,IFG)=ILOCALSTENCIL(N,I,5,ITRy)
				items(3,IFG)=ILOCALSTENCIL(N,I,9,ITRy)
				items(4,IFG)=ILOCALSTENCIL(N,I,11,ITRy)
				items(5,IFG)=ILOCALSTENCIL(N,I,10,ITRy)
                                END IF
			end if
	    end do
	    IF (IELEM(N,I)%ISHAPE.EQ.5)THEN
	    ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4,2:ielem(n,i)%iNUMNEIGHBOURS)=items(3,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,5,2:ielem(n,i)%iNUMNEIGHBOURS)=items(4,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,6:TYPESTEN,:)=0
	    
	    
	    END IF
	    IF (IELEM(N,I)%ISHAPE.EQ.6)THEN
	      ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4,2:ielem(n,i)%iNUMNEIGHBOURS)=items(3,2:ielem(n,i)%iNUMNEIGHBOURS)
	     ILOCALSTENCIL(N,I,5:TYPESTEN,:)=0
	    
	    end if
	    IF (IELEM(N,I)%ISHAPE.EQ.1)THEN
	    ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4,2:ielem(n,i)%iNUMNEIGHBOURS)=items(3,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,5,2:ielem(n,i)%iNUMNEIGHBOURS)=items(4,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,6,2:ielem(n,i)%iNUMNEIGHBOURS)=items(5,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,7,2:ielem(n,i)%iNUMNEIGHBOURS)=items(6,2:ielem(n,i)%iNUMNEIGHBOURS)
	     ILOCALSTENCIL(N,I,8:TYPESTEN,:)=0
	    
	    END IF
	    IF (IELEM(N,I)%ISHAPE.EQ.2)THEN
	    ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4,2:ielem(n,i)%iNUMNEIGHBOURS)=items(3,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,5,2:ielem(n,i)%iNUMNEIGHBOURS)=items(4,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,6:TYPESTEN,:)=0
	    
	    
	    END IF
	      IF ((IELEM(N,I)%ISHAPE.EQ.3).OR.(IELEM(N,I)%ISHAPE.EQ.4))THEN
	    ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4,2:ielem(n,i)%iNUMNEIGHBOURS)=items(3,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,5,2:ielem(n,i)%iNUMNEIGHBOURS)=items(4,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,6,2:ielem(n,i)%iNUMNEIGHBOURS)=items(5,2:ielem(n,i)%iNUMNEIGHBOURS)
	     ILOCALSTENCIL(N,I,7:TYPESTEN,:)=0
	    
	    END IF
	    
	    END IF
	  
	  
	  
	  END DO
	
	
	
	END IF
	
	
	
	
	IF (EES.EQ.4)THEN
	  DO I=1,KMAXE
            ITRR=0

				
	    
	    IF (IELEM(N,I)%ADMIS.EQ.((ielem(n,i)%ifca)+1))THEN
	      ITRR=1    
	    
	    
	    END IF
	    IF ((ITRR.EQ.1).AND.((IELEM(N,I)%ISHAPE.EQ.5).OR.(IELEM(N,I)%ISHAPE.EQ.1)))THEN
              IF (IELEM(N,I)%ISHAPE.EQ.5)THEN
	      IELEM(N,I)%ADMIS=3
	      END IF
	      IF (IELEM(N,I)%ISHAPE.EQ.1)THEN
	      IELEM(N,I)%ADMIS=4
	      END IF
	      
	      ITRX=1;ITRY=1
	      DO IFG=2,ielem(n,i)%iNUMNEIGHBOURS
			IF (MOD(IFG,2).eq.0)THEN
			ITRX=ITRX+1
				IF (IELEM(N,I)%ISHAPE.EQ.5)THEN !QUAD
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,3,ITRX)
                                END IF
				
                                
                                IF (IELEM(N,I)%ISHAPE.EQ.1)THEN !HEXA
				items(1,IFG)=ILOCALSTENCIL(N,I,2,ITRX)
				items(2,IFG)=ILOCALSTENCIL(N,I,4,ITRX)
				items(3,IFG)=ILOCALSTENCIL(N,I,6,ITRX)
                                END IF
                                 
			ELSE
			ITRY=ITRY+1
			    IF (IELEM(N,I)%ISHAPE.EQ.5)THEN !QUAD
			    items(1,IFG)=ILOCALSTENCIL(N,I,4,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,5,ITRY)
 			    END IF
			   
			    IF (IELEM(N,I)%ISHAPE.EQ.1)THEN !hexa
			    items(1,IFG)=ILOCALSTENCIL(N,I,3,ITRY)
			    items(2,IFG)=ILOCALSTENCIL(N,I,5,ITRY)
			    items(3,IFG)=ILOCALSTENCIL(N,I,7,ITRY)
			    
 			    END IF
 			   
			end if
	    end do
	    IF (IELEM(N,I)%ISHAPE.EQ.5)THEN
	    ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	     ILOCALSTENCIL(N,I,4:TYPESTEN,:)=0
	    
	    END IF
	   
	    IF (IELEM(N,I)%ISHAPE.EQ.1)THEN
	    ILOCALSTENCIL(N,I,2,2:ielem(n,i)%iNUMNEIGHBOURS)=items(1,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,3,2:ielem(n,i)%iNUMNEIGHBOURS)=items(2,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4,2:ielem(n,i)%iNUMNEIGHBOURS)=items(3,2:ielem(n,i)%iNUMNEIGHBOURS)
	    ILOCALSTENCIL(N,I,4:TYPESTEN,:)=0
	    END IF
	   
	    
	    END IF
	  
	  
	  
	  END DO
	
	
	
	END IF
	
	
	
	
	
	
	
	
	
	
	
	
	
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	END SUBROUTINE STENCILS

	
	
	SUBROUTINE STENCILS3(N,IELEM,IMAXE,XMPIE,XMPIELRANK,ILOCALSTENCIL,TYPESTEN,NUMNEIGHBOURS,RESTART)
	!> @brief
!> This subroutine is establishing which of the stencils are admissible under different set of rules and which cells can use the WENO algorithms
	IMPLICIT NONE
	TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
	INTEGER,INTENT(IN)::N,IMAXE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::ILOCALSTENCIL
	INTEGER,INTENT(IN)::NUMNEIGHBOURS
	INTEGER,INTENT(IN)::TYPESTEN,RESTART
	INTEGER::I,J,JI,K,LM,KMAXN,KK,KMAXE,IAA,KX,L,ITRR,ITRX,ITRY
	INTEGER::IT1,IT2,IT3,IT4,IT5,IT6,IT7,IT8,ITX,INX,ITF,IFG,KMH,JNG,ichg
	integer,dimension(6,1000)::items

	
	KMAXE=XMPIELRANK(N)
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	DO I=1,KMAXE
		
			KX=1
		DO LM=1,TYPESTEN
				DO IFG=2,ielem(n,i)%iNUMNEIGHBOURS
					KX=KX+1
					ILOCALSTENCIL(N,I,LM,IFG)=ILOCALALLELG(N,I,1,KX)
				END DO
		END DO
		
			ILOCALSTENCIL(N,I,:,1)=ILOCALALLELG(N,I,1,1)
! 			
	END DO	
	
	
      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	END SUBROUTINE STENCILS3


SUBROUTINE ADAPT_CRITERION
!> @brief
!> This subroutine is establishing a region for which to use a very high-order discretisation and a lower one outside this region
	IMPLICIT NONE
	INTEGER::KMAXE,I,FC
	  KMAXE=XMPIELRANK(N)
	  DO I=1,KMAXE
	      FC=0
	      IF (IELEM(N,I)%XXC.LT.-0.025)THEN
	      FC=1

	      END IF
	       IF (IELEM(N,I)%XXC.GT.1.34)THEN

		FC=1
	      END IF
	       IF (IELEM(N,I)%YYC.GT.0.17)THEN
		FC=1

	      END IF

	       IF (IELEM(N,I)%YYC.LT.-0.05)THEN

		FC=1
	      END IF

	IF (FC.EQ.1)THEN
		IELEM(N,I)%ADMIS=1

	END IF


	  END DO




	END SUBROUTINE ADAPT_CRITERION



!---------------------------------------------------------------------------------------------!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE WRITE_BLOCKS(N,XMPIELRANK,XMPINRANK,XMPIE,XMPIN,IELEM,INODE,IMAXN,IMAXE,IBOUND,IMAXB,XMPINNUMBER)
!> @brief
!> This subroutine is solely for debugging purposes
IMPLICIT NONE
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
TYPE(NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::XMPINNUMBER
INTEGER,INTENT(IN)::N,IMAXN,IMAXE
TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIN
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK,XMPINRANK
INTEGER,INTENT(INOUT)::IMAXB
INTEGER::I,J,JI,K,LM,IEX,KMAXN,KK,KMAXE
INTEGER::IT1,IT2,IT3,IT4,IT5,IT6,IT7,IT8,ITX,INX
REAL::X,Y,Z
KMAXE=XMPIELRANK(N)

CALL OPEN_INPUT(N,ITT)
DO I=1,IMAXN
	 READ(9,*)INX,X,Y,Z
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
END DO
CALL CLOSE_INPUT(N,ITT)
CALL OPEN_INPUT(N,ITT)
DO I=1,IMAXN
	READ(9,*)INX,X,Y,Z
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
END DO
CALL CLOSE_INPUT(N,ITT)
CALL OPEN_INPUT(N,ITT)
DO I=1,IMAXN
	READ(9,*)INX,X,Y,Z
!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
END DO
CALL CLOSE_INPUT(N,ITT)
!-------------------FOR DEBUGGING ONLY -----------------------------------------!


END SUBROUTINE WRITE_BLOCKS	




















!---------------------------------------------------------------------------------------------!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED FOR OPENING GRID FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!AND EVERY OTHER RELATED FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED FOR OPENING GRID FILES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!AND EVERY OTHER RELATED FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





SUBROUTINE NEW_ARRAYS(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,KMAXE

J=0
K=0

KMAXE=XMPIELRANK(N)

DO I=1,KMAXE
    IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
    
    J=J+1
    ELSE
    
    K=K+1
    
    END IF
END DO
NOF_INTERIOR=K
NOF_BOUNDED=J

 ALLOCATE(EL_INT(K));EL_INT(:)=0
 ALLOCATE(EL_BND(J));EL_INT(:)=0

J=0
K=0

DO I=1,KMAXE
    IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
    
    J=J+1
    EL_BND(J)=I
    ELSE
    
    K=K+1
    EL_INT(K)=I
    END IF
END DO


END SUBROUTINE NEW_ARRAYS





!---------------------------------------------------------------------------------------------!
END MODULE LIBRARY
