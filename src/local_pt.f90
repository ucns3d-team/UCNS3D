MODULE LOCAL
USE LIBRARY
USE DECLARATION
USE TRANSFORM
IMPLICIT NONE

CONTAINS

SUBROUTINE EXCH_CORDS(N)
!> @brief
!> This subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,INEEDT,TNEEDT,INDL,TNDL,ICPUID,IXFLAG,ITEE,ITEEDUM,IAVC,IAVT,I_CNT,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
REAL,DIMENSION(1)::DUMTS,RUMTS
integer,dimension(4)::icfv1,icfv2
real::rcfv1,rcfv2

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
kmaxe=xmpielrank(n)


INDL=IEXCHANGER(1)%TOT
TNDL=IEXCHANGES(1)%TOT



IF (FASTEST.NE.1)THEN
INEEDT=IRECEXR(1)%TOT
TNEEDT=IRECEXS(1)%TOT
ALLOCATE (IEXCORDR(INEEDT))
ALLOCATE (IEXCORDS(TNEEDT))
ALLOCATE (IEXSOLHIR(INEEDT))
ALLOCATE (IEXSOLHIS(TNEEDT))
END IF
ALLOCATE (IEXBOUNDHIR(INDL))
ALLOCATE (IEXBOUNDHIS(TNDL))
ALLOCATE (IEXBOUNDHIRR(INDL))
ALLOCATE (IEXBOUNDHISS(TNDL))

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)	

if (dimensiona.eq.3)then
i_cnt2=4;i_cnt3=3;i_cnt4=8
else
i_cnt2=2;i_cnt3=2;i_cnt4=4
end if
i_cnt5=i_cnt3*i_cnt4
IF (FASTEST.NE.1)THEN
DO I=1,INEEDT

	IEXSOLHIR(I)%PROCID=IRECEXR(I)%PROCID
 	IEXCORDR(I)%PROCID=IRECEXR(I)%PROCID
 	ALLOCATE (IEXCORDR(I)%NODECORD(IRECEXR(I)%MUCHINEED(1),i_cnt4,i_cnt3))
 	IEXCORDR(I)%NODECORD(1:IRECEXR(I)%MUCHINEED(1),i_cnt4,i_cnt3)=-tolbig
    
		 ALLOCATE(IEXSOLHIR(I)%SOL(IRECEXR(I)%MUCHINEED(1),nof_variables+turbulenceequations+passivescalar))

	IEXSOLHIR(I)%SOL(:,:)=0.0d0
	
END DO
END IF
DO I=1,INDL
	IEXBOUNDHIR(I)%PROCID=IEXCHANGER(I)%PROCID
	IEXBOUNDHIRR(I)%PROCID=IEXCHANGER(I)%PROCID
	IF (ITESTCASE.Le.3)THEN
! 	ALLOCATE(IEXBOUNDHIR(I)%FACESOL(IEXCHANGER(I)%MUCHINEED(1),nof_variables))
	ALLOCATE(IEXBOUNDHIRR(I)%vertpp(IEXCHANGER(I)%MUCHINEED(1),i_cnt2))

	Else
	
	  if (dimensiona.eq.3)then
	 I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((4+TURBULENCEEQUATIONS+PASSIVESCALAR)*3)
	  else
	  I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((3+TURBULENCEEQUATIONS+PASSIVESCALAR)*2)
	  end if

! 	    ALLOCATE(IEXBOUNDHIR(I)%FACESOL(IEXCHANGER(I)%MUCHINEED(1),I_CNT))

	   ALLOCATE(IEXBOUNDHIRR(I)%vertpp(IEXCHANGER(I)%MUCHINEED(1),i_cnt2))
	END IF

! 	IEXBOUNDHIR(I)%FACESOL(:,:)=0.0d0
	IEXBOUNDHIRR(I)%vertpp(:,:)=0

END DO
IF (FASTEST.NE.1)THEN
DO I=1,TNEEDT
	
	IEXSOLHIS(I)%PROCID=IRECEXS(I)%PROCID

	IEXCORDS(I)%PROCID=IRECEXS(I)%PROCID
	ALLOCATE (IEXCORDS(I)%NODECORD(IRECEXS(I)%MUCHTHEYNEED(1),i_cnt4,i_cnt3))
	   


	IEXCORDs(I)%NODECORD(1:IRECEXs(I)%MUCHTHEYNEED(1),1:i_cnt4,1:i_cnt3)=-tolbig
	ALLOCATE (IEXSOLHIS(I)%SOL(IRECEXS(I)%MUCHTHEYNEED(1),nof_variables+turbulenceequations+passivescalar))
    
	IEXSOLHIs(I)%SOL(:,:)=0.0d0
	
END DO
END IF
DO I=1,TNDL


      IEXBOUNDHIs(I)%PROCID=IEXCHANGEs(I)%PROCID
	IEXBOUNDHIss(I)%PROCID=IEXCHANGEs(I)%PROCID
	IF (ITESTCASE.Le.3)THEN
! 	ALLOCATE(IEXBOUNDHIs(I)%FACESOL(IEXCHANGEs(I)%MUCHTHEYNEED(1),nof_variables))
	
	ALLOCATE(IEXBOUNDHIss(I)%vertpp(IEXCHANGEs(I)%MUCHTHEYNEED(1),i_cnt2))

	Else
	
	  if (dimensiona.eq.3)then
	 I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((4+TURBULENCEEQUATIONS+PASSIVESCALAR)*3)
	  else
	  I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((3+TURBULENCEEQUATIONS+PASSIVESCALAR)*2)

	  end if
! 	    ALLOCATE(IEXBOUNDHIs(I)%FACESOL(IEXCHANGEs(I)%MUCHTHEYNEED(1),I_CNT))

	   ALLOCATE(IEXBOUNDHIss(I)%vertpp(IEXCHANGEs(I)%MUCHTHEYNEED(1),i_cnt2))
	END IF

! 	IEXBOUNDHIs(I)%FACESOL(:,:)=0.0d0
	IEXBOUNDHIss(I)%vertpp(:,:)=0
END DO

IF (FASTEST.NE.1)THEN

DO I=1,TNEEDT
! 	
	DO K=1,IRECEXS(I)%MUCHTHEYNEED(1)
! 	
 	    do j=1,ielem(n,IRECEXS(I)%LOCALREF(K))%nonodes
! 		
		IEXCORDS(I)%NODECORD(K,j,1:DIMS)=inoder(ielem(n,IRECEXS(I)%LOCALREF(K))%nodes(j))%CORD(1:DIMS)
! 		
	    end do
	END DO
END DO	




CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)	
ICPUID=N
DUMTS(1:1)=tolsmall
		DO I=0,ISIZE-1
			IF (I.NE.N) THEN
				DO J=1,TNEEDT
					IAVT=10000
					IF (IRECEXS(J)%PROCID.EQ.I)THEN
					IAVT=J 
					GO TO 7001
					END IF
				END DO
				7001 CONTINUE
				DO K=1,INEEDT
					IAVC=10000
					IF (IRECEXR(K)%PROCID.EQ.I) THEN
					IAVC=K
					GO TO 8001
					END IF
				END DO
				8001 CONTINUE
				IF ((IAVT.EQ.10000).AND.(IAVC.NE.10000)) THEN
				CALL MPI_SENDRECV(DUMTS(1:1),1,MPI_DOUBLE_PRECISION,I,ICPUID,&
				IEXCORDR(IAVC)%NODECORD(1:IRECEXR(IAVC)%MUCHINEED(1),1:i_cnt4,1:i_cnt3),&
IRECEXR(IAVC)%MUCHINEED(1)*i_cnt4*i_cnt3,MPI_DOUBLE_PRECISION,IEXCORDR(IAVC)%PROCID,IEXCORDR(IAVC)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)	
				
				END IF
				IF ((IAVT.NE.10000).AND.(IAVC.EQ.10000)) THEN
				!MESSAGE 1
				CALL MPI_SENDRECV(IEXCORDS(IAVT)%NODECORD(1:IRECEXS(IAVT)%MUCHTHEYNEED(1),1:i_cnt4,1:i_cnt3),&
IRECEXS(IAVT)%MUCHTHEYNEED(1)*i_cnt4*i_cnt3,MPI_DOUBLE_PRECISION,IEXCORDS(IAVT)%PROCID,ICPUID,&
				DUMTS(1:1),1,MPI_DOUBLE_PRECISION,I,I,MPI_COMM_WORLD,STATUS,IERROR)
				END IF
				IF ((IAVT.NE.10000).AND.(IAVC.NE.10000)) THEN
				!MESSAGE 1
				CALL MPI_SENDRECV(IEXCORDS(IAVT)%NODECORD(1:IRECEXS(IAVT)%MUCHTHEYNEED(1),1:i_cnt4,1:i_cnt3),&
IRECEXS(IAVT)%MUCHTHEYNEED(1)*i_cnt4*i_cnt3,MPI_DOUBLE_PRECISION,IEXCORDS(IAVT)%PROCID,ICPUID,&
				IEXCORDR(IAVC)%NODECORD(1:IRECEXR(IAVC)%MUCHINEED(1),1:i_cnt4,1:i_cnt3),IRECEXR(IAVC)%MUCHINEED(1)*i_cnt4*i_cnt3,&
MPI_DOUBLE_PRECISION,IEXCORDR(IAVC)%PROCID,IEXCORDR(IAVC)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
! 				
				
				
				END IF
				

			!	END DO !K
			!END DO! J
			END IF	! I.NE.N
		END DO
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)











		

END SUBROUTINE EXCH_CORDS


sUBROUTINE EXCH_CORDS_opt(N)
!> @brief
!> This subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,INEEDT,TNEEDT,INDL,TNDL,ICPUID,IXFLAG,ITEE,ITEEDUM,IAVC,IAVT,I_CNT,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
REAL,DIMENSION(1)::DUMTS,RUMTS
integer,dimension(4)::icfv1,icfv2
real::rcfv1,rcfv2

!-------------------FOR DEBUGGING ONLY -----------------------------------------!

!-------------------FOR DEBUGGING ONLY -----------------------------------------!
kmaxe=xmpielrank(n)


INDL=IEXCHANGER(1)%TOT
TNDL=IEXCHANGES(1)%TOT



if (dimensiona.eq.3)then
i_cnt2=4;i_cnt3=3;i_cnt4=8
else
i_cnt2=2;i_cnt3=2;i_cnt4=4
end if
i_cnt5=i_cnt3*i_cnt4

DO I=1,INDL

	IF (ITESTCASE.Le.3)THEN
        ALLOCATE(IEXBOUNDHIR(I)%FACESOL(IEXCHANGER(I)%MUCHINEED(1),nof_variables))
    ! 	ALLOCATE(IEXBOUNDHIRR(I)%vertpp(IEXCHANGER(I)%MUCHINEED(1),i_cnt2))


        IF (DG == 1)THEN
            ALLOCATE(IEXBOUNDHIR(I)%FACESOL_DG(IEXCHANGER(I)%MUCHINEED(1),NOF_VARIABLES))
        END IF
        
        if (mood.eq.1)then
            ALLOCATE(IEXBOUNDHIR(I)%FACESOL_m(IEXCHANGER(I)%MUCHINEED(1),1))
        end if

	Else
	
	  if (dimensiona.eq.3)then
	 I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((4+TURBULENCEEQUATIONS+PASSIVESCALAR)*3)
	  else
	  I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((3+TURBULENCEEQUATIONS+PASSIVESCALAR)*2)
	  end if

	    ALLOCATE(IEXBOUNDHIR(I)%FACESOL(IEXCHANGER(I)%MUCHINEED(1),I_CNT))

            
	END IF

	IEXBOUNDHIR(I)%FACESOL(:,:)=0.0d0
! 	IEXBOUNDHIRR(I)%vertpp(:,:)=0

END DO

DO I=1,TNDL
	IF (ITESTCASE.Le.3)THEN
        ALLOCATE(IEXBOUNDHIs(I)%FACESOL(IEXCHANGEs(I)%MUCHTHEYNEED(1),nof_variables))
        
        IF (DG == 1)THEN
            ALLOCATE(IEXBOUNDHIS(I)%FACESOL_DG(IEXCHANGER(I)%MUCHINEED(1),NOF_VARIABLES))
        END IF
        
        if (mood.eq.1)then
            ALLOCATE(IEXBOUNDHIs(I)%FACESOL_m(IEXCHANGEs(I)%MUCHTHEYNEED(1),1))
        end if

	Else
	
	  if (dimensiona.eq.3)then
	 I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((4+TURBULENCEEQUATIONS+PASSIVESCALAR)*3)
	  else
	  I_CNT=(nof_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)+((3+TURBULENCEEQUATIONS+PASSIVESCALAR)*2)

	  end if
	    ALLOCATE(IEXBOUNDHIs(I)%FACESOL(IEXCHANGEs(I)%MUCHTHEYNEED(1),I_CNT))

	   
	    
	  
	END IF

	IEXBOUNDHIs(I)%FACESOL(:,:)=0.0d0
	
END DO


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)











		

END SUBROUTINE EXCH_CORDS_opt





SUBROUTINE EXCH_CORD3(N)
!> @brief
!> This subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values and mapping of gaussian quadrature points
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,INEEDT,TNEEDT,INDL,TNDL,ICPUID,IXFLAG,ITEE,ITEEDUM,IAVC,IAVT,I_CNT,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
REAL,DIMENSION(1)::DUMTS,RUMTS
integer,dimension(4)::icfv1,icfv2
real::rcfv1,rcfv2

ICPUID=N

kmaxe=xmpielrank(n)

INDL=IEXCHANGER(1)%TOT
TNDL=IEXCHANGES(1)%TOT

if (dimensiona.eq.3)then
i_cnt2=4;i_cnt3=3;i_cnt4=8
else
i_cnt2=2;i_cnt3=2;i_cnt4=4
end if
i_cnt5=i_cnt3*i_cnt4

if (dimensiona.eq.3)then



!now try to remap the local and global neighbours and mapping of the gaussian quadrature points


!first the global
DO J=1,INDL
DO I=1,TNDL
IF (IEXCHANGER(J)%PROCID.EQ.IEXCHANGES(I)%PROCID)THEN
DO K=1,IEXCHANGES(I)%MUCHTHEYNEED(1)

! 	    
	    if (ielem(n,(IEXCHANGES(I)%LOCALREF(K)))%TYPES_FACES(IEXCHANGES(I)%SIDETHEYNEED(K)).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
! 	    
	    IEXBOUNDHISS(I)%VERTPP(K,1:ixf4)=ielem(n,(IEXCHANGES(I)%LOCALREF(K)))%NODES_FACES(IEXCHANGES(I)%SIDETHEYNEED(K),1:ixf4)


	    

end do
END IF
end do
END DO

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

DO K=1,INDL
DO J=1,TNDL
      IF (IEXCHANGER(k)%PROCID.EQ.IEXCHANGES(j)%PROCID)THEN

      CALL MPI_SENDRECV(IEXBOUNDHISs(J)%vertpp(1:IEXCHANGES(J)%MUCHTHEYNEED(1),1:i_cnt2)&
,IEXCHANGES(J)%MUCHTHEYNEED(1)*I_cnt2,MPI_INTEGER,IEXCHANGES(j)%PROCID,&
      IEXCHANGES(j)%PROCID,IEXBOUNDHIRr(K)%vertpp(1:IEXCHANGER(K)%MUCHINEED(1),1:i_cnt2),&
IEXCHANGER(K)%MUCHINEED(1)*I_cnt2,MPI_INTEGER,IEXCHANGES(j)%PROCID,&
      ICPUID,MPI_COMM_WORLD,STATUS,IERROR)
      END IF
END DO
END DO







do i=1,kmaxe
  if (ielem(n,i)%interior.eq.1)then
  do k=1,ielem(n,i)%ifca
		  if (ielem(n,i)%TYPES_FACES(k).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
	    
      if (ielem(n,i)%ineighg(K).gt.0)then
	if (ielem(n,i)%ineighb(K).ne.n)then
	    
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
! 		ielem(n,i)%NODES_FACES(k,1:ixf4)=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),1:ixf4)
! 		ielem(n,i)%REORIENT(K)=1
	    if (ielem(n,i)%ibounds(k).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(k))%icode.eq.5)then
		do ixfv=1,ixf4
		inoder(IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),ixfv))%itor=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),ixfv)
		end do
	      end if
	    end if
	    end if

	END IF
	eND IF
	
      END DO
    END IF
  
end do










else







!now try to remap the local and global neighbours and mapping of the gaussian quadrature points


!first the global
DO J=1,INDL
DO I=1,TNDL
IF (IEXCHANGER(J)%PROCID.EQ.IEXCHANGES(I)%PROCID)THEN
DO K=1,IEXCHANGES(I)%MUCHTHEYNEED(1)

ixf4=2

IEXBOUNDHISs(I)%vertpp(K,1:ixf4)=ielem(n,(IEXCHANGES(I)%LOCALREF(K)))%NODES_FACES(IEXCHANGES(I)%SIDETHEYNEED(K),1:ixf4)
end do
END IF
end do
END DO


DO K=1,INDL
DO J=1,TNDL
      IF (IEXBOUNDHIRR(K)%PROCID.EQ.IEXBOUNDHISs(J)%PROCID)THEN
      CALL MPI_SENDRECV(IEXBOUNDHISs(J)%vertpp(1:IEXCHANGES(J)%MUCHTHEYNEED(1),1:i_cnt2)&
,IEXCHANGES(J)%MUCHTHEYNEED(1)*I_cnt2,MPI_INTEGER,IEXBOUNDHISs(J)%PROCID,&
      ICPUID,IEXBOUNDHIRr(K)%vertpp(1:IEXCHANGER(K)%MUCHINEED(1),1:i_cnt2),&
IEXCHANGER(K)%MUCHINEED(1)*I_cnt2,MPI_INTEGER,IEXBOUNDHIRr(K)%PROCID,&
      IEXBOUNDHIRR(K)%PROCID,MPI_COMM_WORLD,STATUS,IERROR)
      END IF
END DO
END DO

















do i=1,kmaxe
  if (ielem(n,i)%interior.eq.1)then
  do k=1,ielem(n,i)%ifca
		 
	    ixf4=2
	    
      if (ielem(n,i)%ineighg(K).gt.0)then
	if (ielem(n,i)%ineighb(K).ne.n)then
	    
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then

! 		ielem(n,i)%NODES_FACES(k,1:ixf4)=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),1:ixf4)
! 		ielem(n,i)%REORIENT(K)=1
	    if (ielem(n,i)%ibounds(k).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(k))%icode.eq.5)then
		do ixfv=1,ixf4
 		inoder(IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),ixfv))%itor=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),ixfv)
		end do
	      end if
	    end if
	    end if


	else
! 	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
! 		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
! 		ielem(n,i)%REORIENT(K)=1
! 	    end if
	end if
      end if
  end do
  else
!       do k=1,ielem(n,i)%ifca
! 
! 
! 	  ixf4=2
!       if (ielem(n,i)%ineighg(K).gt.0)then
! 	
! 	
! 	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
! 		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
! 		ielem(n,i)%REORIENT(K)=1
! 	    end if
! 	
!       end if
!       end do
  end if
end do


end if

END SUBROUTINE EXCH_CORD3


SUBROUTINE EXCH_CORDS2(N,ISIZE,IEXBOUNDHIRi,IEXBOUNDHISi,&
ITESTCASE,NUMBEROFPOINTS2,IEXCHANGER,IEXCHANGES)
!> @brief
!> This subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values in 2D
IMPLICIT NONE
TYPE(EXCHANGE_BOUNDHI),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXBOUNDHIRi
TYPE(EXCHANGE_BOUNDHI),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXBOUNDHISi
TYPE(EXCHANGE),ALLOCATABLE,DIMENSION(:),INTENT(IN)::IEXCHANGER
TYPE(EXCHANGE),ALLOCATABLE,DIMENSION(:),INTENT(IN)::IEXCHANGES
INTEGER,INTENT(IN)::ISIZE,N,ITESTCASE,NUMBEROFPOINTS2
INTEGER::I,J,K,L,INEEDT,TNEEDT,INDL,TNDL,ICPUID,IXFLAG,ITEE,ITEEDUM,IAVC,IAVT,I_CNT
REAL,DIMENSION(1:1)::DUMTS,RUMTS

INDL=IEXCHANGER(1)%TOT
TNDL=IEXCHANGES(1)%TOT


ALLOCATE (IEXBOUNDHIRi(INDL))
ALLOCATE (IEXBOUNDHISi(TNDL))


I_CNT=(nof_variables+turbulenceequations+passivescalar)
DO I=1,INDL
	IEXBOUNDHIRi(I)%PROCID=IEXCHANGER(I)%PROCID
	
	
	ALLOCATE(IEXBOUNDHIRi(I)%FACESOL(IEXCHANGER(I)%MUCHINEED(1),I_CNT))
	IEXBOUNDHIRi(I)%FACESOL(:,:)=0.0d0
	
END DO

DO I=1,TNDL
	IEXBOUNDHISi(I)%PROCID=IEXCHANGES(I)%PROCID
	
	ALLOCATE(IEXBOUNDHISi(I)%FACESOL(IEXCHANGES(I)%MUCHTHEYNEED(1),I_CNT))
	
	IEXBOUNDHISi(I)%FACESOL(:,:)=0.0d0
	
END DO


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		

END SUBROUTINE EXCH_CORDS2





SUBROUTINE FIND_ROT_ANGLES(N,ICONSI)
!> @brief
!> This subroutine determines the normal vectors for each face
IMPLICIT NONE
real::Xc1,Yc1,Zc1,Xc2,Yc2,Zc2,Xc3,Yc3,Zc3,DELXYA,DELyzA,DELzxA,DELXYb,DELyzb,DELzxb,DELXYc,DELyzc,DELzxc,nx,ny,nz
REAL::X5,X6,X7,X8,Y5,Y6,Y7,Y8,Z5,Z6,Z7,Z8,XX,YY,ZZ
REAL::DELXA,DELXB,DELYA,DELYB,DELZA,DELZB,L_ANGLE1,L_ANGLE2,lNX,LNY,LNZ
INTEGER::K,KMAXE,I,J,kk,kk2,ixf4,IXFV
INTEGER,INTENT(IN)::N,ICONSI
KMAXE=XMPIELRANK(N)

i=iconsi

			DO K=1,IELEM(N,I)%IFCA
			       if (ielem(n,i)%types_faces(k).eq.5)then
				    kk2=4;n_node=kk2
			      else
				    kk2=3;n_node=kk2
			      end if
			    IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
			    IF ((IELEM(N,I)%INEIGHG(K).GT.0).AND.(IELEM(N,I)%IBOUNDS(K).GT.0))THEN 	!PERIODIC NEIGHBOUR
		             
			    XX=IELEM(N,I)%XXC  ;YY=IELEM(N,I)%YYC; ZZ=IELEM(N,I)%ZZC
			     
				    DO Kk=1,n_node
				       IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN	
				       vext(kk,1:3)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:3)
				       ELSE
					
				       
				       
					vext(kk,1:3)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:3)
!  					
				       END IF
				      IF(ABS(vext(kk,1)-xx).GT.XPER*oo2)THEN
				      vext(kk,1)=vext(kk,1)+(XPER*SIGN(1.0,xx-XPER*oo2))
				      end if
				      IF(ABS(vext(kk,2)-yy).GT.yPER*oo2)THEN
				      vext(kk,2)=vext(kk,2)+(yPER*SIGN(1.0,yy-yPER*oo2))
				      end if
				      IF(ABS(vext(kk,3)-zz).GT.zPER*oo2)THEN
				      vext(kk,3)=vext(kk,3)+(zPER*SIGN(1.0,zz-zPER*oo2))
				      end if
				      
! 				    
				      
				      
				      
			      end do

			    Else
				  DO Kk=1,n_node
					 IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN	
				       vext(kk,1:3)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:3)
				       ELSE
					vext(KK,1:3)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:3)
				       END IF
				  END DO

			    END IF
			    ELSE
				   DO Kk=1,n_node
					 IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN	
				       vext(kk,1:3)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:3)
				       ELSE
					vext(KK,1:3)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:3)
				       END IF
				  END DO




			    END IF
					IF (KK2.EQ.3)THEN  
					Xc1=veXt(1,1); Xc2=veXt(2,1); Xc3=veXt(3,1);
					Yc1=veXt(1,2);Yc2=veXt(2,2); Yc3=veXt(3,2);
					Zc1=veXt(1,3); Zc2=veXt(2,3); Zc3=veXt(3,3);
					DELXYA=(xc1-xc2)*(yc1+yc2);DELyzA=(yc1-yc2)*(zc1+zc2);DELzxA=(zc1-zc2)*(xc1+xc2)
					DELXYb=(xc2-xc3)*(yc2+yc3);DELyzb=(yc2-yc3)*(zc2+zc3);DELzxb=(zc2-zc3)*(xc2+xc3)
					DELXYc=(xc3-xc1)*(yc3+yc1);DELyzc=(yc3-yc1)*(zc3+zc1);DELzxc=(zc3-zc1)*(xc3+xc1)
					nx=(delyza+delyzb+delyzc)
					ny=(delzxa+delzxb+delzxc)
					nz=(delxya+delxyb+delxyc)
					ROOT_ROT=SQRT((nx**2)+(ny**2)+(nz**2))
					nx=nx/root_ROT; ny=ny/root_ROT; nz=nz/root_ROT
					root_ROT=1.0D0
					a_ROT=nx
					b_ROT=ny
					c_ROT=nz
					CALL ANGLEX(ANGLEFACEX)
					CALL ANGLEY(ANGLEFACEY)
					IELEM(N,I)%FACEANGLEX(K)=anglefacex
					IELEM(N,I)%FACEANGLEY(K)=anglefacey

					ELSE

					 
					
					
					
					DELXA=VEXT(4,1)-VEXT(2,1)
					DELXB=VEXT(3,1)-VEXT(1,1)
					DELYA=VEXT(4,2)-VEXT(2,2)
					DELYB=VEXT(3,2)-VEXT(1,2)
					DELZA=VEXT(4,3)-VEXT(2,3)
					DELZB=VEXT(3,3)-VEXT(1,3)
					
					
					NX=-0.50D0*((DELYA*DELZB)-(DELZA*DELYB))
					NY=-0.50D0*((DELZA*DELXB)-(DELXA*DELZB))
					NZ=-0.50D0*((DELXA*DELYB)-(DELYA*DELXB))
					
					
					
					!newells method

					NX=ZERO;NY=ZERO;NZ=ZERO;ROOT_ROT=ZERO
 					do kk=1,n_node
 					if (kk.ne.n_node)then
 					nx=nx+(vext(kk,2)-vext(kk+1,2))*(vext(kk,3)+vext(kk+1,3))
 					ny=ny+(vext(kk,3)-vext(kk+1,3))*(vext(kk,1)+vext(kk+1,1))
 					nz=nz+(vext(kk,1)-vext(kk+1,1))*(vext(kk,2)+vext(kk+1,2))
					else
 					nx=nx+(vext(kk,2)-vext(1,2))*(vext(kk,3)+vext(1,3))
 					ny=ny+(vext(kk,3)-vext(1,3))*(vext(kk,1)+vext(1,1))
 					nz=nz+(vext(kk,1)-vext(1,1))*(vext(kk,2)+vext(1,2))
 					
 					end if
 					end do
                                            
					
					
					
! 					
					root_ROT=SQRT((nx**2)+(ny**2)+(nz**2))
					nx=nx/root_ROT; ny=ny/root_ROT; nz=nz/root_ROT
					root_ROT=1.0D0
					a_ROT=nx
					b_ROT=ny
					c_ROT=nz
					CALL ANGLEX(ANGLEFACEX)
					CALL ANGLEY(ANGLEFACEY)
					
					
					IELEM(N,I)%FACEANGLEX(K)=anglefacex
					IELEM(N,I)%FACEANGLEY(K)=anglefacey
					END IF
                    IELEM(N,I)%LUMP=0
                    if (ielem(n,i)%ishape.eq.2)then
                    L_ANGLE1=IELEM(N,I)%FACEANGLEX(K);L_ANGLE2=IELEM(N,I)%FACEANGLEY(K)
                    LNX=COS(L_ANGLE1)*SIN(L_ANGLE2)
                    LNY=SIN(L_ANGLE1)*SIN(L_ANGLE2)
                    LNZ=COS(L_ANGLE2)
                    if ((abs(LNX-1.0d0).le.10e-16).OR.(abs(LNY-1.0d0).le.10e-16).OR.(abs(LNZ-1.0d0).le.10e-16))THEN
                    IELEM(N,I)%LUMP=100
                    end if
                    END IF


			    END DO







END SUBROUTINE FIND_ROT_ANGLES



SUBROUTINE FIND_ROT_ANGLES2d(N,ICONSI)
!> @brief
!> This subroutine determines the normal vectors for each edge
IMPLICIT NONE
real::Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,DELXYA,DELyzA,DELzxA,DELXYb,DELyzb,DELzxb,DELXYc,DELyzc,DELzxc,nx,ny,nz
REAL::X5,X6,X7,X8,Y5,Y6,Y7,Y8,Z5,Z6,Z7,Z8,XX,YY,ZZ
REAL::DELXA,DELXB,DELYA,DELYB,DELZA,DELZB
INTEGER::K,KMAXE,I,J,kk,kk2,ixf4,IXFV
INTEGER,INTENT(IN)::N,ICONSI
KMAXE=XMPIELRANK(N)

i=iconsi
	
	
	i=iconsi

			DO K=1,IELEM(N,I)%IFCA
! 			       if (ielem(n,i)%types_faces(k).eq.5)then
				    kk2=2;n_node=kk2
! 			      else
! 				    kk2=3;n_node=kk2
! 			      end if
			    IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
			    IF ((IELEM(N,I)%INEIGHG(K).GT.0).AND.(IELEM(N,I)%IBOUNDS(K).GT.0))THEN 	!PERIODIC NEIGHBOUR
		             
 			    XX=IELEM(N,I)%XXC  ;YY=IELEM(N,I)%YYC; !ZZ=IELEM(N,I)%ZZC
			     
				    DO Kk=1,n_node
				       IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN	
				       vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:2)
				       ELSE
					
! 				     
					vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:2)
					
! 					
				       END IF
				       
				       
				       
				       
				       
				      IF(ABS(vext(kk,1)-xx).GT.XPER*oo2)THEN
				      vext(kk,1)=vext(kk,1)+(XPER*SIGN(1.0d0,xx-XPER/2.0D0))
				      end if
				      IF(ABS(vext(kk,2)-yy).GT.yPER*oo2)THEN
				      vext(kk,2)=vext(kk,2)+(yPER*SIGN(1.0d0,yy-yPER/2.0D0))
				      end if
				      
				      
! 				     
				      
				      
				      
			      end do

			    Else
				  DO Kk=1,n_node
					 IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN	
				       vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:2)
				       ELSE
					vext(KK,1:2)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:2)
				       END IF
				  END DO

			    END IF
			    ELSE
				   DO Kk=1,n_node
					 IF (IELEM(N,I)%REORIENT(K).EQ.0)THEN	
				       vext(kk,1:2)=inoder(ielem(n,i)%NODES_FACES(k,kk))%CORD(1:2)
				       ELSE
					vext(KK,1:2)=inoder(ielem(n,i)%NODES_FACES(k,n_node-KK+1))%CORD(1:2)
				       END IF
				  END DO




			    END IF
					
			
			
					CALL ANGLE2D(ANGLEFACEX,ANGLEFACEY)
					
					
! 					
					IELEM(N,I)%FACEANGLEX(K)=anglefacex
					IELEM(N,I)%FACEANGLEY(K)=anglefacey
			
! 					
			



			    END DO
	
	
	
	
	
	
	
	
	
	
	
	
	
! 			







END SUBROUTINE FIND_ROT_ANGLES2d








SUBROUTINE LOCALISE_STENCIL(N,Iconsi)
!> @brief
!> This subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell
IMPLICIT NONE
INTEGER,INTENT(IN)::n,iconsi
INTEGER::I,J,l,IXFF,IXSST,ikg,ismp,inv,ineedt,ikg2,IN_STEN,K,itarget
REAL::RIN_STEN
INEEDT=IRECEXR(1)%TOT
i=iconsi
	IKG2=0
	DO ISMP=1,TYPESTEN
		IKG=0
		
                        if ((ees.ne.5).or.(ismp.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        end if
                        
		
			DO L=1,itarget
				IF (ILOCALSTENCIL(N,I,ISMP,L).GT.0)THEN
				IKG=IKG+1
				END IF
			END DO
			IF (IKG.EQ.itarget)THEN
			
			
			
				if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
				IKG2=IKG2+1
				
				DO L=1,itarget
				ILOCAL_ELEM(1)%IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				j=(XMPIL(ILOCAL_ELEM(1)%IHEXG(IKG2,L)))
				ILOCAL_ELEM(1)%IHEXL(IKG2,L)=J
				ILOCAL_ELEM(1)%ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
				CALL COMPUTE_CENTRE3d(N,j)
				ilocal_elem(1)%XXC(IKG2,L)=CORDS(1)
				ilocal_elem(1)%YYC(IKG2,L)=CORDS(2)
				ilocal_elem(1)%ZZC(IKG2,L)=CORDS(3)
				
				ILOCAL_NODE(1)%NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes) 
				
				DO K=1,ielem(n,j)%nonodes    
				ILOCAL_NODE(1)%x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1) 
				ILOCAL_NODE(1)%y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2) 
				ILOCAL_NODE(1)%z(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(3) 
				END DO
				end do
				
				
				
				ELSE
				IKG2=IKG2+1
				DO L=1,itarget
				ILOCAL_ELEM(1)%IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)	
				ILOCAL_ELEM(1)%IHEXB(IKG2,L)=XMPIE(ILOCALSTENCIL(N,I,ISMP,L))
				IF (ILOCAL_ELEM(1)%IHEXB(IKG2,L).EQ.N)THEN
				      j=(XMPIL(ILOCAL_ELEM(1)%IHEXG(IKG2,L)))
						IF (ILOCAL_ELEM(1)%IHEXG(IKG2,L).EQ.IELEM(N,J)%IHEXGL)THEN
						ILOCAL_ELEM(1)%IHEXL(IKG2,L)=J
						ILOCAL_ELEM(1)%ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
						END IF
						CALL COMPUTE_CENTRE3d(N,j)
						ilocal_elem(1)%XXC(IKG2,L)=CORDS(1)
						ilocal_elem(1)%YYC(IKG2,L)=CORDS(2)
						ilocal_elem(1)%ZZC(IKG2,L)=CORDS(3)
						ILOCAL_NODE(1)%NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)  
						DO K=1,ielem(n,j)%nonodes    
				ILOCAL_NODE(1)%x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1) 
				ILOCAL_NODE(1)%y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2) 
				ILOCAL_NODE(1)%z(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(3) 
				
				
				
				END DO
				
				
				
				
				ELSE
				
				
				
					
					 DO IXFF=1,INEEDT
						IF (IEXCORDR(IXFF)%PROCID.EQ.ilocal_elem(1)%IHEXB(IKG2,L))THEN
							DO IXSST=1,IRECEXR(IXFF)%MUCHINEED(1)
							      IF ((ilocal_elem(1)%IHEXG(IKG2,L)).EQ.IRECEXR1(IXFF)%WHATINEED(IXSST))THEN
								ilocal_elem(1)%IHEXL(IKG2,L)=IXSST
								ilocal_elem(1)%IHEXN(IKG2,L)=IXFF
								ilocal_elem(1)%ISHAPE(IKG2,L)=IRECEXR1(IXFF)%ISHAPE(IXSST)
								EXIT
							       END IF
							END DO
						END IF
					 END DO	
					SELECT CASE(ilocal_elem(1)%ISHAPE(IKG2,L))
					CASE(1)
					IN_STEN=8
					RIN_STEN=8.0D0
					CASE(2)
					IN_STEN=4
					RIN_STEN=4.0D0

					CASE(3)
					IN_STEN=5
					RIN_STEN=5.0D0

					CASE(4)
					IN_STEN=6
					RIN_STEN=6.0D0
					END SELECT
					ILOCAL_NODE(1)%X(IKG2,L,1:IN_STEN)=IEXCORDR(ilocal_elem(1)%IHEXN(IKG2,L))%NODECORD(ilocal_elem(1)%IHEXL(IKG2,L),1:IN_STEN,1)
					ILOCAL_NODE(1)%Y(IKG2,L,1:IN_STEN)=IEXCORDR(ilocal_elem(1)%IHEXN(IKG2,L))%NODECORD(ilocal_elem(1)%IHEXL(IKG2,L),1:IN_STEN,2)
					ILOCAL_NODE(1)%Z(IKG2,L,1:IN_STEN)=IEXCORDR(ilocal_elem(1)%IHEXN(IKG2,L))%NODECORD(ilocal_elem(1)%IHEXL(IKG2,L),1:IN_STEN,3)
					ilocal_elem(1)%XXC(IKG2,L)=sum(ILOCAL_NODE(1)%X(IKG2,L,1:IN_STEN))/RIN_STEN
					ilocal_elem(1)%YYC(IKG2,L)=sum(ILOCAL_NODE(1)%Y(IKG2,L,1:IN_STEN))/RIN_STEN
					ilocal_elem(1)%ZZC(IKG2,L)=sum(ILOCAL_NODE(1)%Z(IKG2,L,1:IN_STEN))/RIN_STEN
					
					
					
					
					
					
! 			
					
				END IF
				END DO
				END IF
				END IF
				END DO

END SUBROUTINE LOCALISE_STENCIL


SUBROUTINE LOCALISE_STENCIL2d(N,Iconsi)
!> @brief
!> This subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::n,iconsi
INTEGER::I,J,l,IXFF,IXSST,ikg,ismp,inv,ineedt,ikg2,IN_STEN,K,itarget
REAL::RIN_STEN
INEEDT=IRECEXR(1)%TOT
i=iconsi
	IKG2=0
	DO ISMP=1,TYPESTEN
		IKG=0
                             if ((ees.ne.5).or.(ismp.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        end if
		
		
			DO L=1,itarget
				IF (ILOCALSTENCIL(N,I,ISMP,L).GT.0)THEN
				IKG=IKG+1
				END IF
			END DO
			IF (IKG.EQ.itarget)THEN
				if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
				IKG2=IKG2+1
				DO L=1,itarget
				ILOCAL_ELEM(1)%IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				j=(XMPIL(ILOCAL_ELEM(1)%IHEXG(IKG2,L)))
				ILOCAL_ELEM(1)%IHEXL(IKG2,L)=J
				ILOCAL_ELEM(1)%ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
				CALL COMPUTE_CENTRE2d(N,j)
				ilocal_elem(1)%XXC(IKG2,L)=CORDS(1)
				ilocal_elem(1)%YYC(IKG2,L)=CORDS(2)
				
				
				ILOCAL_NODE(1)%NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes) 
				DO K=1,ielem(n,j)%nonodes    
				ILOCAL_NODE(1)%x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1) 
				ILOCAL_NODE(1)%y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2) 
				
				END DO
				end do
				ELSE
				IKG2=IKG2+1
				DO L=1,itarget
				ILOCAL_ELEM(1)%IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				ILOCAL_ELEM(1)%IHEXB(IKG2,L)=XMPIE(ILOCALSTENCIL(N,I,ISMP,L))
				IF (ILOCAL_ELEM(1)%IHEXB(IKG2,L).EQ.N)THEN
				      j=(XMPIL(ILOCAL_ELEM(1)%IHEXG(IKG2,L)))
						IF (ILOCAL_ELEM(1)%IHEXG(IKG2,L).EQ.IELEM(N,J)%IHEXGL)THEN
						ILOCAL_ELEM(1)%IHEXL(IKG2,L)=J
						ILOCAL_ELEM(1)%ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
						END IF
						CALL COMPUTE_CENTRE2d(N,j)
						ilocal_elem(1)%XXC(IKG2,L)=CORDS(1)
						ilocal_elem(1)%YYC(IKG2,L)=CORDS(2)
						
						ILOCAL_NODE(1)%NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)  
						DO K=1,ielem(n,j)%nonodes    
				ILOCAL_NODE(1)%x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1) 
				ILOCAL_NODE(1)%y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2) 
				
				END DO
				ELSE
					 DO IXFF=1,INEEDT
						IF (IEXCORDR(IXFF)%PROCID.EQ.ilocal_elem(1)%IHEXB(IKG2,L))THEN
							DO IXSST=1,IRECEXR(IXFF)%MUCHINEED(1)
							      IF ((ilocal_elem(1)%IHEXG(IKG2,L)).EQ.IRECEXR1(IXFF)%WHATINEED(IXSST))THEN
								ilocal_elem(1)%IHEXL(IKG2,L)=IXSST
								ilocal_elem(1)%IHEXN(IKG2,L)=IXFF
								ilocal_elem(1)%ISHAPE(IKG2,L)=IRECEXR1(IXFF)%ISHAPE(IXSST)
								EXIT
							       END IF
							END DO
						END IF
					 END DO	
					SELECT CASE(ilocal_elem(1)%ISHAPE(IKG2,L))
					CASE(5)
					IN_STEN=4
					RIN_STEN=4.D0
					CASE(6)
					IN_STEN=3
					RIN_STEN=3.D0

					
					END SELECT
					ILOCAL_NODE(1)%X(IKG2,L,1:IN_STEN)=IEXCORDR(ilocal_elem(1)%IHEXN(IKG2,L))%NODECORD(ilocal_elem(1)%IHEXL(IKG2,L),1:IN_STEN,1)
					ILOCAL_NODE(1)%Y(IKG2,L,1:IN_STEN)=IEXCORDR(ilocal_elem(1)%IHEXN(IKG2,L))%NODECORD(ilocal_elem(1)%IHEXL(IKG2,L),1:IN_STEN,2)
					
					ilocal_elem(1)%XXC(IKG2,L)=sum(ILOCAL_NODE(1)%X(IKG2,L,1:IN_STEN))/RIN_STEN
					ilocal_elem(1)%YYC(IKG2,L)=sum(ILOCAL_NODE(1)%Y(IKG2,L,1:IN_STEN))/RIN_STEN
					
				END IF
				END DO
				END IF
				END IF
				END DO

END SUBROUTINE LOCALISE_STENCIL2d




SUBROUTINE  LOCALISE_STEN2(N,ICONSI)
!> @brief
!> This subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell
IMPLICIT NONE
INTEGER::I,J,K,L,KK,PRK,JJ,kmaxe,ineedt,jx2,jx,iivd,iivd3,facexx,IXXFFf,in1,itarget,idum
real,dimension(3)::tempcentres,TEMP_cG
real::dumv1,dumv2,detjc,dist1,MA,MB,MC,MD,ME,MF,MG,MH,MI,MDD
INTEGER,INTENT(IN)::ICONSI,N
INEEDT=IRECEXR(1)%TOT
i=iconsi




	IF (IPERIODICITY.EQ.1)THEN

	  DO JJ=1,IELEM(N,I)%ADMIS
	  
             if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        end if
	  
	  
	    DO J=2,itarget
	    IF(ABS(ilocal_elem(1)%XXC(JJ,J)-ilocal_elem(1)%XXC(JJ,1)).GT.XPER*oo2)THEN
		    ilocal_elem(1)%XXC(JJ,J)=ilocal_elem(1)%XXC(JJ,J)+(XPER*SIGN(1.0,ilocal_elem(1)%XXC(JJ,1)-XPER*oo2))
		    DO KK=1,8
		    ILOCAL_NODE(1)%X(JJ,J,KK)=ILOCAL_NODE(1)%X(JJ,J,KK)+(XPER*SIGN(1.0,ilocal_elem(1)%XXC(JJ,1)-XPER*oo2))
		    END DO
	    END IF
	    IF(ABS(ilocal_elem(1)%YYC(JJ,J)-ilocal_elem(1)%YYC(JJ,1)).GT.YPER*oo2)THEN
		    ilocal_elem(1)%YYC(JJ,J)=ilocal_elem(1)%YYC(JJ,J)+(YPER*SIGN(1.0,ilocal_elem(1)%YYC(JJ,1)-YPER*oo2))
		  
		    DO KK=1,8
		    ILOCAL_NODE(1)%Y(JJ,J,KK)=ILOCAL_NODE(1)%Y(JJ,J,KK)+(YPER*SIGN(1.0,ilocal_elem(1)%YYC(JJ,1)-YPER*oo2))
		    END DO
		    
	    END IF
	    IF(ABS(ilocal_elem(1)%ZZC(JJ,J)-ilocal_elem(1)%ZZC(JJ,1)).GT.ZPER*oo2)THEN
		    ilocal_elem(1)%ZZC(JJ,J)=ilocal_elem(1)%ZZC(JJ,J)+(ZPER*SIGN(1.0,ilocal_elem(1)%ZZC(JJ,1)-ZPER*oo2))
		      DO KK=1,8
		    ILOCAL_NODE(1)%Z(JJ,J,KK)=ILOCAL_NODE(1)%Z(JJ,J,KK)+(ZPER*SIGN(1.0,ilocal_elem(1)%ZZC(JJ,1)-ZPER*oo2))
		    END DO
		   
	    END IF
	    END DO
		
	END DO
	END IF
	
	  
      VEXT=0.0d0
      NODES_LIST=0.0d0
      ELTYPE=IELEM(N,I)%ISHAPE
      ELEM_DEC=IELEM(N,I)%VDEC
      ELEM_LISTD=0.0d0
      jx=IELEM(N,I)%NONODES
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILOCAL_NODE(1)%x(1,1,K)
	    NODES_LIST(k,2)=ILOCAL_NODE(1)%y(1,1,K)
	    NODES_LIST(k,3)=ILOCAL_NODE(1)%z(1,1,K)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3
    
      SELECT CASE(ielem(n,i)%ishape)
      CASE(1)
      if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREHEXA(N,IGQRULES)
      DUMV1=HEXAVOLUME(N)
      CALL COMPUTE_CENTRE3d(N,i)
      vext(1,1:dims)=cords(1:dims)
      else
	VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	  call COMPUTEJACOBIANS
      end if	    
      CASE(2)
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	call COMPUTEJACOBIANS
      CASE(3)
	VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	  call COMPUTEJACOBIANS
      CASE(4)
       if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREPRISM(N,IGQRULES)
      DUMV1=PRISMVOLUME(N)
      VEXT(1,1:3)=VEXT(6,1:3)
      else
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
      call COMPUTEJACOBIANS
      end if
      END SELECT
		  TEMP_CG(1:3)=VEXT(1,1:3)
		  ILOCAL_RECON3(I)%INVCCJAC(1:3,1:3)=VVa1(1:3,1:3)
		  

		  
		  
		  
		  
		  detjc=deta(1)
	    ILOCAL_RECON3(I)%VEXT_REF(1:3)=TEMP_cG(1:3)
	    
	    
	    if ((poly.eq.4).OR.(DG.EQ.1))then
		    ILOCAL_RECON3(I)%INVCCJAC(1:3,1:3)=zero
		     ILOCAL_RECON3(I)%INVCCJAC(1,1)=1.0d0
		      ILOCAL_RECON3(I)%INVCCJAC(2,2)=1.0d0
		       ILOCAL_RECON3(I)%INVCCJAC(3,3)=1.0d0
		      detjc=1.0
		      
		      TEMP_cG(1)=ielem(n,i)%xxc
		      TEMP_cG(2)=ielem(n,i)%yyc
		      TEMP_cG(3)=ielem(n,i)%zzc
		      ILOCAL_RECON3(I)%VEXT_REF(1:3)=TEMP_cG(1:3)
		    end if
	    
	    
	    
	    
	    
	    

		DO JJ=1,IELEM(N,I)%ADMIS
		
                        if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        end if
		
			DO L=1,itarget
				
				
				
				if (fastest.ne.1)then
				
				
				
				DO KK=1,8
					TEMPCENTRES=ZERO
					TEMPCENTRES(1)=ILOCAL_NODE(1)%X(JJ,L,KK)
					TEMPCENTRES(2)=ILOCAL_NODE(1)%Y(JJ,L,KK)
					TEMPCENTRES(3)=ILOCAL_NODE(1)%Z(JJ,L,KK)
					TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-TEMP_CG(:))
					ILOCAL_NODE(1)%X(JJ,L,KK)=TEMPCENTRES(1)
					ILOCAL_NODE(1)%Y(JJ,L,KK)=TEMPCENTRES(2)
					ILOCAL_NODE(1)%Z(JJ,L,KK)=TEMPCENTRES(3)
				END DO
				end if
			END DO	
		END DO	
		IDUM=0
		if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
		
		


		

		
		DO JJ=1,IELEM(N,I)%ADMIS
			if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        end if
		
			DO L=1,itarget
				if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
				ilocal_elem(1)%VOLUME(JJ,L)=(IELEM(N,ilocal_elem(1)%IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				
				
				
				
				if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ilocal_elem(1)%VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				end if
				
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				else
				
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				
				end if
				
! 				
				Else
				IF (ilocal_elem(1)%IHEXB(JJ,L).eq.N)THEN
				ilocal_elem(1)%VOLUME(JJ,L)=(IELEM(N,ilocal_elem(1)%IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ilocal_elem(1)%VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				end if
				
				
	
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
 				else
                                ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				

                                end if

				else 
				if ((EES.ne.5).or.(jj.eq.1))then
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXN(JJ,L)=ilocal_elem(1)%IHEXN(JJ,L)
				else
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXNc(JJ,L)=ilocal_elem(1)%IHEXN(JJ,L)
				
				end if
				
				
				
				ELTYPE=ilocal_elem(1)%ISHAPE(jj,L)
				ELEM_LISTD=0.0d0; VEXT=0.0d0; NODES_LIST=0.0d0
				      
				      select case(ELTYPE)
				      
				      case(1)
				      ELEM_DEC=6; jx=8
				      case(2)
				    ELEM_DEC=1; jx=4
				      case(3)
					ELEM_DEC=2; jx=5
				      case(4)
				      ELEM_DEC=3; jx=6
				      end select

					    do K=1,jx
					      NODES_LIST(k,1)=ILOCAL_NODE(1)%x(jj,l,K)
					      NODES_LIST(k,2)=ILOCAL_NODE(1)%y(jj,l,K)
					      NODES_LIST(k,3)=ILOCAL_NODE(1)%z(jj,l,K)
					      VEXT(K,:)=NODES_LIST(k,:)
					      
					    END DO
					    CALL DECOMPOSE3
					    
					    dumv2=0.0d0
					    do k=1,ELEM_DEC
					    VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
					    dumv2=dumv2+TETRAVOLUME(N)
					    
					    end do
					    
					    	
					
					    
								    
					    
					    
					    ilocal_elem(1)%VOLUME(JJ,L)=dumv2
				

				if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ilocal_elem(1)%VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				end if
				
				!ILOCAL_RECON3(I)%VOLUME(JJ,L)=ilocal_elem(1)%VOLUME(JJ,L)
				
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				!ILOCAL_RECON3(I)%VOLUMEc(JJ,L)=ilocal_elem(1)%VOLUME(JJ,L)
				end if
				
				

				end if
				END IF
! 							
			END DO	
		END DO

		if (ielem(n,i)%interior.eq.0)then
	  CALL COMPUTE_CENTRE3d(N,i)
	  vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
! 		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(1,1:DIMS)-VEXT(2,1:DIMS)
		    end if
	  end do
      else
		    CALL COMPUTE_CENTRE3d(N,i)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		  select case(ielem(n,i)%types_faces(k))
		  case(5)
		  IXXFFf=4
		  case(6)
		  IXXFFf=3
		  end select
		  call COMPUTE_CENTRE3dF(N,I,facexx,IXXFFf)
		  VEXT(2,1:dims)=cords(1:dims)
	
		  
		      dist1=distance3(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
! 		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(1,1:DIMS)-VEXT(2,1:DIMS)
		    end if
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
! 		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(1,1:DIMS)-VEXT(2,1:DIMS)
		    end if
		else						!from another cpu 
		    DO In1=1,ielem(n,i)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ilocal_elem(1)%XXC(1,In1);vext(2,2)=ilocal_elem(1)%yyC(1,In1); vext(2,3)=ilocal_elem(1)%zzC(1,In1)	      
		     dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
! 		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(1,1:DIMS)-VEXT(2,1:DIMS)
				      end if
			  end if
		    end do
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(N,j)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    IF(ABS(vext(2,3)-vext(1,3)).GT.zPER*oo2)THEN
		    vext(2,3)=vext(2,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
		    end if
		    dist1=distance3(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
! 		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(1,1:DIMS)-VEXT(2,1:DIMS)
		    end if
		else	!periodic boundaries from another cpu
		     DO In1=1,ielem(n,i)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ilocal_elem(1)%XXC(1,In1);vext(2,2)=ilocal_elem(1)%yyC(1,In1); vext(2,3)=ilocal_elem(1)%zzC(1,In1)	      
		     
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    IF(ABS(vext(2,3)-vext(1,3)).GT.zPER*oo2)THEN
		    vext(2,3)=vext(2,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
		    end if
		    dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
! 		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(1,1:DIMS)-VEXT(2,1:DIMS)
				end if
			  end if
		    end do

		end if
		end if
	  end do
	end if


        DO JJ=1,IELEM(N,I)%ADMIS
                        if ((EES.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        
                        end if
			DO L=1,itarget
				if (fastest.ne.1)then
				TEMPCENTRES=ZERO
				TEMPCENTRES(1)=ilocal_elem(1)%XXC(JJ,L)
				TEMPCENTRES(2)=ilocal_elem(1)%YYC(JJ,L)
				TEMPCENTRES(3)=ilocal_elem(1)%ZZC(JJ,L)
				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-TEMP_CG(:))
				ilocal_elem(1)%XXC(JJ,L)=TEMPCENTRES(1)
				ilocal_elem(1)%YYC(JJ,L)=TEMPCENTRES(2)
				ilocal_elem(1)%ZZC(JJ,L)=TEMPCENTRES(3)

				
				
				
				
				
				
				
				
				
                        
				end if
			END DO	
		END DO					

	

		
	
! 
END SUBROUTINE LOCALISE_STEN2



SUBROUTINE  LOCALISE_STEN2d(N,ICONSI)
!> @brief
!> This subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
IMPLICIT NONE
INTEGER::I,J,K,L,KK,PRK,JJ,kmaxe,ineedt,jx2,jx,in1,facexx,ixxfff,IHGT,IHGJ,ITARGET,IDUM
real,dimension(2)::tempcentres,rel2
real::dumv1,dumv2,detjc,dist1
INTEGER,INTENT(IN)::ICONSI,N

KMAXE=XMPIELRANK(N)
INEEDT=IRECEXR(1)%TOT
i=iconsi
	IF (IPERIODICITY.EQ.1)THEN

	  DO JJ=1,IELEM(N,I)%ADMIS
	  if ((EES.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        
                        end if
	    DO J=2,Itarget
	    IF(ABS(ilocal_elem(1)%XXC(JJ,J)-ilocal_elem(1)%XXC(JJ,1)).GT.XPER*oo2)THEN
		    ilocal_elem(1)%XXC(JJ,J)=ilocal_elem(1)%XXC(JJ,J)+(XPER*SIGN(1.0,ilocal_elem(1)%XXC(JJ,1)-XPER*oo2))
		    DO KK=1,4
		    ILOCAL_NODE(1)%X(JJ,J,KK)=ILOCAL_NODE(1)%X(JJ,J,KK)+(XPER*SIGN(1.0,ilocal_elem(1)%XXC(JJ,1)-XPER*oo2))
		    END DO
	    END IF
	    IF(ABS(ilocal_elem(1)%YYC(JJ,J)-ilocal_elem(1)%YYC(JJ,1)).GT.YPER*oo2)THEN
		    ilocal_elem(1)%YYC(JJ,J)=ilocal_elem(1)%YYC(JJ,J)+(YPER*SIGN(1.0,ilocal_elem(1)%YYC(JJ,1)-YPER*oo2))
		  
		    DO KK=1,4
		    ILOCAL_NODE(1)%Y(JJ,J,KK)=ILOCAL_NODE(1)%Y(JJ,J,KK)+(YPER*SIGN(1.0,ilocal_elem(1)%YYC(JJ,1)-YPER*oo2))
		    END DO
		    
	    END IF
	    END DO
		
	END DO
	END IF
	
	  
     VEXT=0.0d0
    NODES_LIST=0.0d0
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=0.0d0
       
      jx=IELEM(N,I)%NONODES
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILOCAL_NODE(1)%x(1,1,K)
	    NODES_LIST(k,2)=ILOCAL_NODE(1)%y(1,1,K)
	   
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE2
    
      SELECT CASE(ielem(n,i)%ishape)

      CASE(5)
      if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREQUAD(N,IGQRULES)
      DUMV1=QUADVOLUME(N)
      CALL COMPUTE_CENTRE2d(N,i)
      vext(1,1:dims)=cords(1:dims)
      else
	VEXT(1:3,1:2)=ELEM_LISTD(1,1:3,1:2)
! 	  DUMV1=TRIANGLEVOLUME(N)
	  call COMPUTEJACOBIANS2
	end if
      

      CASE(6)
      VEXT(1:3,1:2)=ELEM_LISTD(1,1:3,1:2)
! 	DUMV1=TRIANGLEVOLUME(N)
	call COMPUTEJACOBIANS2

      END SELECT
      
		
		  ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2)=VVa1(1:2,1:2)
		  
		  
		  
		  
		  
		  
		  
		 
		    
 		IDUM=0
		if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
                                IDUM=1
                            end if
                        END IF
                        END DO
                end if
		   
		   
		   
		   
		  
        
		  
		  
		  
            detjc=deta(1)
		    ILOCAL_RECON3(I)%VEXT_REF(1:2)=VEXT(1,1:2)
		    
		    if ((poly.eq.4).OR.(DG.EQ.1))then
		    ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2)=zero
		     ILOCAL_RECON3(I)%INVCCJAC(1,1)=1.0d0
		      ILOCAL_RECON3(I)%INVCCJAC(2,2)=1.0d0
		      detjc=1.0
		      VEXT(1,1)=ielem(n,i)%xxc
		      VEXT(1,2)=ielem(n,i)%yyc
		      ILOCAL_RECON3(I)%VEXT_REF(1:2)=VEXT(1,1:2)
		    end if
		    
		    
		  DO JJ=1,IELEM(N,I)%ADMIS
                            if ((EES.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        
                        end if
		  
		  
			DO L=1,Itarget
! 				TEMPCENTRES=ZERO
! 				TEMPCENTRES(1)=ilocal_elem(1)%XXC(JJ,L)
! 				TEMPCENTRES(2)=ilocal_elem(1)%YYC(JJ,L)
! 				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
! 				
! 				ilocal_elem(1)%XXC(JJ,L)=TEMPCENTRES(1)
! 				ilocal_elem(1)%YYC(JJ,L)=TEMPCENTRES(2)
				
				
				
				
				DO KK=1,4
					TEMPCENTRES=ZERO
					TEMPCENTRES(1)=ILOCAL_NODE(1)%X(JJ,L,KK)
					TEMPCENTRES(2)=ILOCAL_NODE(1)%Y(JJ,L,KK)
					
					
					TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
					
					ILOCAL_NODE(1)%X(JJ,L,KK)=TEMPCENTRES(1)
					ILOCAL_NODE(1)%Y(JJ,L,KK)=TEMPCENTRES(2)
					
				END DO
				
			END DO	
		END DO			

					
		DO JJ=1,IELEM(N,I)%ADMIS
                        if ((EES.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        
                        end if
			DO L=1,Itarget
				if (ILOCAL_RECON3(I)%LOCAL.eq.1)then
				ilocal_elem(1)%VOLUME(JJ,L)=(IELEM(N,ilocal_elem(1)%IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				 if ((EES.ne.5).or.(jj.eq.1))then
				 if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ilocal_elem(1)%VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				end if
				
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				
				
				end if
				Else
				IF (ilocal_elem(1)%IHEXB(JJ,L).eq.N)THEN
				ilocal_elem(1)%VOLUME(JJ,L)=(IELEM(N,ilocal_elem(1)%IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				 if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ilocal_elem(1)%VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				end if
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				end if
				else 
				 if ((EES.ne.5).or.(jj.eq.1))then
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXN(JJ,L)=ilocal_elem(1)%IHEXN(JJ,L)
				else
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ilocal_elem(1)%IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ilocal_elem(1)%IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ilocal_elem(1)%IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXNc(JJ,L)=ilocal_elem(1)%IHEXN(JJ,L)
				
				
				end if
				
				      ELTYPE=ilocal_elem(1)%ISHAPE(jj,L)
				      ELEM_LISTD=0.0d0; VEXT=0.0d0; NODES_LIST=0.0d0
				      select case(ELTYPE)
				      
				      case(5)
				      ELEM_DEC=2; jx=4
				      case(6)
				    ELEM_DEC=1; jx=3
				      
				      end select

					    do K=1,jx
					      NODES_LIST(k,1)=ILOCAL_NODE(1)%x(jj,l,K)
					      NODES_LIST(k,2)=ILOCAL_NODE(1)%y(jj,l,K)
					     
					      VEXT(K,:)=NODES_LIST(k,:)
					    END DO
					    CALL DECOMPOSE2
					    dumv2=0.0d0
					    do k=1,ELEM_DEC
					    VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
					    dumv2=dumv2+TRIANGLEVOLUME(N)
					    end do
					    ilocal_elem(1)%VOLUME(JJ,L)=dumv2
					     if ((EES.ne.5).or.(jj.eq.1))then
						 if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ilocal_elem(1)%VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)
				end if
				
                                            else
                                  ILOCAL_RECON3(I)%VOLUME(1,1)=ilocal_elem(1)%VOLUME(1,1)         
                                            
                                            end if
! 				
				end if
				END IF
			END DO	
		END DO

	      if (ielem(n,i)%interior.eq.0)then
	  CALL COMPUTE_CENTRE2d(N,i)
	  vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
		    end if
	  end do
      else
		    CALL COMPUTE_CENTRE2d(N,i)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		  
		  IXXFFf=2
		  
		  
		  call COMPUTE_CENTRE2dF(N,iconsi,facexx,IXXFFf)
		  VEXT(2,1:dims)=cords(1:dims)
		  
		      dist1=distance2(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
		    end if
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
		    end if
		else						!from another cpu 
		    DO In1=1,IELEM(N,I)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ilocal_elem(1)%XXC(1,In1);vext(2,2)=ilocal_elem(1)%yyC(1,In1)    
		     dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
				      end if
			  end if
		    end do
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(N,j)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    
		    dist1=distance2(n)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
		    end if
		else	!periodic boundaries from another cpu
		     DO In1=1,IELEM(N,I)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ilocal_elem(1)%XXC(1,In1);vext(2,2)=ilocal_elem(1)%yyC(1,In1);     
		     
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    
		    dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
				end if
			  end if
		    end do

		end if
		end if
	  end do
	end if

        DO JJ=1,IELEM(N,I)%ADMIS
			if ((EES.ne.5).or.(jj.eq.1))then
                        itarget=ielem(n,i)%iNUMNEIGHBOURS
                        
                        else
                        itarget=NUMNEIGHBOURS2
                        
                        
                        end if
			DO L=1,Itarget
				TEMPCENTRES=ZERO
				TEMPCENTRES(1)=ilocal_elem(1)%XXC(JJ,L)
				TEMPCENTRES(2)=ilocal_elem(1)%YYC(JJ,L)
				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
				
				ilocal_elem(1)%XXC(JJ,L)=TEMPCENTRES(1)
				ilocal_elem(1)%YYC(JJ,L)=TEMPCENTRES(2)
				
				
				
				
				
				
			END DO	
		END DO


		
	
! 
END SUBROUTINE LOCALISE_STEN2d










subroutine direct_side(n)
!> @brief
!> This subroutine establishes the distance betwen cell centres for each face
implicit none
integer,intent(in)::n
integer::i,j,k,kmaxe,facexx,ixxfff
real::dist1

kmaxe=xmpielrank(n)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,facexx,ixxfff) 
!$OMP DO SCHEDULE (STATIC)
do i=1,kmaxe
	if (ielem(n,i)%interior.eq.0)then
		    CALL COMPUTE_CENTRE3d(N,i)
		    vext(1,1:dims)=cords(1:dims)
      
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
	  end do
	
	else
		    CALL COMPUTE_CENTRE3d(N,i)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
		  
		  facexx=k
		  select case(ielem(n,i)%types_faces(k))
		  case(5)
		  IXXFFf=4
		  case(6)
		  IXXFFf=3
		  end select
		  call COMPUTE_CENTRE3dF(N,I,facexx,IXXFFf)
		  VEXT(2,1:dims)=cords(1:dims)
! 		  j=IELEM(N,i)%INEIGH(K)
		  
		      dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
		else						!from another cpu 
		
		    vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims)
		     dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(N,j)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    IF(ABS(vext(2,3)-vext(1,3)).GT.zPER*oo2)THEN
		    vext(2,3)=vext(2,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
		    end if
		    dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
		
		else	!periodic boundaries from another cpu

		     vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims) 
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    IF(ABS(vext(2,3)-vext(1,3)).GT.zPER*oo2)THEN
		    vext(2,3)=vext(2,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
		    end if
		    dist1=distance3(n)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
	  end do
	end if

end do
!$OMP END DO
!$OMP END PARALLEL


end subroutine direct_side



subroutine direct_side2d(n)
!> @brief
!> This subroutine establishes the distance betwen cell centres for each edge
implicit none
integer,intent(in)::n
integer::i,j,k,kmaxe,facexx,ixxfff
real::dist1

kmaxe=xmpielrank(n)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,facexx,ixxfff) 
!$OMP DO
do i=1,kmaxe
	if (ielem(n,i)%interior.eq.0)then
		    CALL COMPUTE_CENTRE2d(N,i)
		    vext(1,1:dims)=cords(1:dims)
      
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
	  end do
	
	else
		    CALL COMPUTE_CENTRE2d(N,i)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		 
		  IXXFFf=2
		 
		  call COMPUTE_CENTRE2dF(N,I,facexx,IXXFFf)
		  VEXT(2,1:dims)=cords(1:dims)
		  
		      dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(N,j)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
		else						!from another cpu 
		    vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims)
		     dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(N,j)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER/2.d0)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0D0,vext(1,1)-XPER/2.D0))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER/2.d0)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0D0,vext(1,2)-yPER/2.D0))
		    end if
		    
		    dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
		
		else	!periodic boundaries from another cpu

		     vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims) 
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0D0,vext(1,1)-XPER/2.0D0))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0D0,vext(1,2)-yPER/2.0D0))
		    end if
		    
		    dist1=distance2(n)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
	  end do
	end if

end do
!$OMP END DO
!$OMP END PARALLEL


end subroutine direct_side2d

SUBROUTINE CHECKGRADS(N,ICONSI)
!> @brief
!> This subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSI
REAL::DXX1,dxx2,TEMPG1,dist1,dist2,oo2,surfmin,surfmax
INTEGER::I,J,K,L,jj,icount3,nnd,ixf4,IDC,IDC2

i=iconsi
    IELEM(N,I)%GGS=greengo
    dxx1=-TOLBIG; dxx2=TOLBIG
    CALL COMPUTE_CENTRE3d(N,i)
    vext(1,1:dims)=CORDS(1:dims)
      
I=ICONSI
  if (ielem(n,i)%interior.eq.1)then
  do k=1,ielem(n,i)%ifca
		  if (ielem(n,i)%TYPES_FACES(k).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
	    
      if (ielem(n,i)%ineighg(K).gt.0)then
	if (ielem(n,i)%ineighb(K).ne.n)then
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then

		ielem(n,i)%NODES_FACES(k,1:ixf4)=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),1:ixf4)
		ielem(n,i)%REORIENT(K)=1

	    end if


	else
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then

		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
		  ielem(n,i)%REORIENT(K)=1
! 		
	    end if
	end if
      end if
      
  end do
  else
      do k=1,ielem(n,i)%ifca
		  if (ielem(n,i)%TYPES_FACES(k).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if

      if (ielem(n,i)%ineighg(K).gt.0)then
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
! 		
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
! 		
	    end if
	
      end if

      
      end do
  end if
	
























    
surfmin=1.0e16
surfmax=1.0e-16



       
      do l=1,ielem(n,i)%ifca
            
	      select case (ielem(n,i)%types_faces(l))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select

		
		  if (ielem(n,i)%interior.eq.1)then
		    if ((ielem(n,i)%ineighg(l).gt.0).and.(ielem(n,i)%ibounds(l).gt.0)) then
! 		     

                    end if
                    else
! 		
		       call COMPUTE_CENTRE3dF(N,I,l,nnd)
		  VEXT(2,1:dims)=cords(1:dims)
		
		end if
		
		
		
	      surfmin=min(surfmin,IELEM(N,I)%SURF(L))
	      surfmax=max(surfmax,IELEM(N,I)%SURF(L))
	      
	      
	      
	      
		    
	


		  dist1=distance3(n)
		if (dist1.lt.dxx2)then
		  dxx2=dist1
		end if
		if (dist1.gt.dxx1)then
		  dxx1=dist1
		end if
end do
        ielem(n,i)%condition=1.0

        if (turbulence.gt.0)then
        

        if (ielem(n,i)%ishape.eq.2)then
        ielem(n,i)%condition=surfmax/surfmin
        
        if (ielem(n,i)%condition.gt.30)then
         ielem(n,i)%hybrid=1
        end if
        end if
        
        if (ielem(n,i)%ishape.eq.3)then
        ielem(n,i)%condition=surfmax/surfmin
        
        if (ielem(n,i)%condition.gt.10)then
         ielem(n,i)%hybrid=1
        end if
        
       
        
        end if


        end if
       



	    TEMPG1=MAX((DXX1/DXX2),(DXX2/DXX1))
! 	   
			
	    
	      IF (TEMPG1.GT.GRIDAR1)THEN
	      IELEM(N,I)%GGS=1
	      end if 
   if (fastest.eq.0)then
	       dxx1=-tolbig; dxx2=tolbig
	       JJ=1
	       DO L=1,ielem(n,i)%iNUMNEIGHBOURS
		       if (ilocal_elem(1)%VOLUME(JJ,L).lt.dxx2)then
		      
		       dxx2=ilocal_elem(1)%VOLUME(JJ,L)
		       end if
		       if (ilocal_elem(1)%VOLUME(JJ,L).gt.dxx1)then
		      
		       dxx1=ilocal_elem(1)%VOLUME(JJ,L)
		       end if
	       end do
	 TEMPG1=MAX((DXX1/DXX2),(DXX2/DXX1))
	     IF (TEMPG1.GT.GRIDAR2)THEN
	       IELEM(N,I)%GGS=1
	       end if  
	      
! ! 	      
 end if

IDC=0
idc2=0
if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
	        IDC=IDC+1
	      Else
idc2=idc2+1

		end if
	  END IF
        END DO
END IF

IF ((IDC.Gt.1))THEN     !until GE is fully adaptive
IELEM(N,I)%GGS=1

END IF







END SUBROUTINE CHECKGRADS






SUBROUTINE CHECK3(N,ICONSI)
!> @brief
!> This subroutine assigns the ordering for the faces
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,INTENT(INout)::ICONSI
REAL::DXX1,dxx2,TEMPG1,dist1,dist2,oo2
INTEGER::I,J,K,L,jj,icount3,nnd,ixf4


do iconsi=1,xmpielrank(n)

i=iconsi
    IELEM(N,I)%GGS=greengo
    dxx1=-TOLBIG; dxx2=TOLBIG
    CALL COMPUTE_CENTRE3d(N,i)
    vext(1,1:dims)=CORDS(1:dims)
      
I=ICONSI
  if (ielem(n,i)%interior.eq.1)then
  do k=1,ielem(n,i)%ifca
		  if (ielem(n,i)%TYPES_FACES(k).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
	    
      if (ielem(n,i)%ineighg(K).gt.0)then
	if (ielem(n,i)%ineighb(K).ne.n)then
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
	    
	    end if


	else
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
	    end if
	end if
      end if
!       
  end do
  else
      do k=1,ielem(n,i)%ifca
		  if (ielem(n,i)%TYPES_FACES(k).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if

      if (ielem(n,i)%ineighg(K).gt.0)then
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
	    end if
	
      end if

!      

      end do
  end if
end do

end subroutine check3



SUBROUTINE CHECKGRADS2d(N,ICONSI)
!> @brief
!> This subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSI
REAL::DXX1,dxx2,TEMPG1,dist1,dist2,oo2
INTEGER::I,J,K,L,jj,icount3,nnd,ixf4
i=iconsi

tempg1=0.0; 
    IELEM(N,I)%GGS=greengo
    dxx1=-TOLBIG; dxx2=TOLBIG
    CALL COMPUTE_CENTRE2d(N,i)
    vext(1,1:dims)=CORDS(1:dims)
      
      
      
     

  if (ielem(n,i)%interior.eq.1)then
  do k=1,ielem(n,i)%ifca
		  
	    ixf4=2
	   
	    
      if (ielem(n,i)%ineighg(K).gt.0)then
	if (ielem(n,i)%ineighb(K).ne.n)then
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IEXBOUNDHIRR(IELEM(N,I)%INEIGHN(K))%VERTPP(IELEM(N,I)%Q_FACE(K)%Q_MAPL(1),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
	    
	    end if


	else
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
	    end if
	end if
      end if
  end do
  else
      do k=1,ielem(n,i)%ifca
		ixf4=2

      if (ielem(n,i)%ineighg(K).gt.0)then
	    if (ielem(n,i)%ineighg(K).gt.ielem(n,i)%ihexgl)then
		ielem(n,i)%NODES_FACES(k,1:ixf4)=IELEM(N,ielem(n,i)%ineigh(K))%NODES_FACES(IELEM(N,I)%INEIGHN(K),1:ixf4)
		ielem(n,i)%REORIENT(K)=1
	    end if
	
      end if
      end do
  end if


     
      do l=1,ielem(n,i)%ifca
      
	      
	      nnd=2
	     

		
		  if (ielem(n,i)%interior.eq.1)then
		    if ((ielem(n,i)%ineighg(l).gt.0).and.(ielem(n,i)%ibounds(l).gt.0))then
		      if (ielem(n,i)%ineighb(l).ne.n)then
           
			do K=1,nnd
			  NODES_LIST(k,1:dims)=inoder(IELEM(N,Iconsi)%NODES_FACES(l,K))%CORD(1:dims)
			END DO
			
			
			do K=1,nnd
			IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
			NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0d0,vext(1,1)-XPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
			NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0d0,vext(1,2)-yPER*oo2))
			end if
			end do
			
			CORDS=CORDINATES2(N,NODES_LIST,nnd)
			vext(2,1:dims)=CORDS(1:dims)
		    else
		      call COMPUTE_CENTRE2dF(N,I,l,nnd)
		      VEXT(2,1:dims)=cords(1:dims)
		    end if
		  else
		
		    
		      call COMPUTE_CENTRE2dF(N,I,l,nnd)
		  VEXT(2,1:dims)=cords(1:dims)
		  end if
		else
		       call COMPUTE_CENTRE2dF(N,I,l,nnd)
		  VEXT(2,1:dims)=cords(1:dims)
		
		end if
	      
		    
	


		  dist1=distance2(n)
		if (dist1.lt.dxx2)then
		  dxx2=dist1
		end if
		if (dist1.gt.dxx1)then
		  dxx1=dist1
		end if
		
		
		
end do
	    IF (TEMPG1.GT.GRIDAR1)THEN
	      IELEM(N,I)%GGS=1
	      end if 
	      
	      
	      
   if (fastest.eq.0)then
	       dxx1=-tolbig; dxx2=tolbig
	       JJ=1
	       DO L=1,ielem(n,i)%iNUMNEIGHBOURS
		       if (ilocal_elem(1)%VOLUME(JJ,L).lt.dxx2)then
		      
		       dxx2=ilocal_elem(1)%VOLUME(JJ,L)
		       end if
		       if (ilocal_elem(1)%VOLUME(JJ,L).gt.dxx1)then
		      
		       dxx1=ilocal_elem(1)%VOLUME(JJ,L)
		       end if
	       end do
	 TEMPG1=MAX((DXX1/DXX2),(DXX2/DXX1))
	     IF (TEMPG1.GT.GRIDAR2)THEN
	       IELEM(N,I)%GGS=1
	       end if 
end if


END SUBROUTINE CHECKGRADS2d







END MODULE LOCAL
