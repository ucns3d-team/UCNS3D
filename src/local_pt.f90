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


if (adda.eq.1)then
ALLOCATE (IEXSOLHIRd(INEEDT))
ALLOCATE (IEXSOLHISd(TNEEDT))
end if

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

		 if (adda.eq.1)then
		ALLOCATE(IEXSOLHIRd(I)%SOL(IRECEXR(I)%MUCHINEED(1),1))
		end if
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

	if (adda.eq.1)then
	ALLOCATE(IEXSOLHIsd(I)%SOL(IRECEXS(I)%MUCHTHEYNEED(1),1))
	end if

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




CALL EXCHANGE_CORDX(N,INEEDT,TNEEDT,i_cnt4,i_cnt3)









END SUBROUTINE EXCH_CORDS



SUBROUTINE EXCHANGE_CORDX(N,INEEDT,TNEEDT,i_cnt4,i_cnt3)
INTEGER,INTENT(IN)::N,INEEDT,TNEEDT,i_cnt4,i_cnt3
INTEGER::I,J,K,L,ICPUID,IXFLAG,ITEE,ITEEDUM,IAVC,IAVT
REAL,DIMENSION(1)::DUMTS,RUMTS

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



END SUBROUTINE EXCHANGE_CORDX


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
        IF (DG == 1)THEN
            ALLOCATE(IEXBOUNDHIR(I)%FACESOL_DG(IEXCHANGER(I)%MUCHINEED(1),I_CNT))
        END IF


            
	END IF

	IEXBOUNDHIR(I)%FACESOL(:,:)=0.0d0
	if (DG.EQ.1)THEN
	IEXBOUNDHIR(I)%FACESOL_DG(:,:)=0.0d0
	END IF
! 	IEXBOUNDHIRR(I)%vertpp(:,:)=0

END DO

DO I=1,TNDL
	IF (ITESTCASE.Le.3)THEN
        ALLOCATE(IEXBOUNDHIs(I)%FACESOL(IEXCHANGEs(I)%MUCHTHEYNEED(1),nof_variables))
        
        IF (DG == 1)THEN
            ALLOCATE(IEXBOUNDHIS(I)%FACESOL_DG(IEXCHANGEs(I)%MUCHTHEYNEED(1),NOF_VARIABLES))
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

	   IF (DG == 1)THEN
            ALLOCATE(IEXBOUNDHIS(I)%FACESOL_DG(IEXCHANGEs(I)%MUCHTHEYNEED(1),I_CNT))
        END IF
	    
	  
	END IF

	IEXBOUNDHIS(I)%FACESOL(:,:)=0.0d0
	if (DG.EQ.1)THEN
	IEXBOUNDHIS(I)%FACESOL_DG(:,:)=0.0d0
	END IF

	
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
	      if ((ibound(n,ielem(n,i)%ibounds(k))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(k))%icode.eq.50))then
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
	      if ((ibound(n,ielem(n,i)%ibounds(k))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(k))%icode.eq.50))then
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
TYPE(A_EXCHANGE_BOUNDHI),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXBOUNDHIRi
TYPE(A_EXCHANGE_BOUNDHI),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::IEXBOUNDHISi
TYPE(A_EXCHANGE),ALLOCATABLE,DIMENSION(:),INTENT(IN)::IEXCHANGER
TYPE(A_EXCHANGE),ALLOCATABLE,DIMENSION(:),INTENT(IN)::IEXCHANGES
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














SUBROUTINE LOCALISE_STENCIL(N,Iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell
IMPLICIT NONE
INTEGER,INTENT(IN)::n,iconsi
INTEGER::I,J,l,IXFF,IXSST,ikg,ismp,inv,ineedt,ikg2,IN_STEN,K,itarget
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X
REAL::RIN_STEN
REAL,dimension(1:dimensiona)::CORDS
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
				ILOX_IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				ILOX_PERIODICFLAG(IKG2,L)=ILOCALSTENCILPER(N,I,ISMP,L)
				j=(XMPIL(ILOX_IHEXG(IKG2,L)))
				ILOX_IHEXL(IKG2,L)=J
				ILOX_ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
				CALL COMPUTE_CENTRE3d(j,cords)
				ILOX_XXC(IKG2,L)=CORDS(1)
				ILOX_YYC(IKG2,L)=CORDS(2)
				ILOX_ZZC(IKG2,L)=CORDS(3)
				
				ILON_NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)
				
				DO K=1,ielem(n,j)%nonodes    
				ILON_x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1)
				ILON_y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2)
				ILON_z(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(3)
				END DO
				end do
				
				
				
				ELSE
				IKG2=IKG2+1
				DO L=1,itarget
				ILOX_IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				ILOX_PERIODICFLAG(IKG2,L)=ILOCALSTENCILPER(N,I,ISMP,L)
				ILOX_IHEXB(IKG2,L)=XMPIE(ILOCALSTENCIL(N,I,ISMP,L))
				IF (ILOX_IHEXB(IKG2,L).EQ.N)THEN
				      j=(XMPIL(ILOX_IHEXG(IKG2,L)))
						IF (ILOX_IHEXG(IKG2,L).EQ.IELEM(N,J)%IHEXGL)THEN
						ILOX_IHEXL(IKG2,L)=J
						ILOX_ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
						END IF
						CALL COMPUTE_CENTRE3d(j,cords)
						ILOX_XXC(IKG2,L)=CORDS(1)
						ILOX_YYC(IKG2,L)=CORDS(2)
						ILOX_ZZC(IKG2,L)=CORDS(3)
						ILON_NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)
						DO K=1,ielem(n,j)%nonodes    
				ILON_x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1)
				ILON_y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2)
				ILON_z(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(3)
				
				
				
				END DO
				
				
				
				
				ELSE
				
				
				
					
					 DO IXFF=1,INEEDT
						IF (IEXCORDR(IXFF)%PROCID.EQ.ILOX_IHEXB(IKG2,L))THEN
							DO IXSST=1,IRECEXR(IXFF)%MUCHINEED(1)
							      IF ((ILOX_IHEXG(IKG2,L)).EQ.IRECEXR1(IXFF)%WHATINEED(IXSST))THEN
								ILOX_IHEXL(IKG2,L)=IXSST
								ILOX_IHEXN(IKG2,L)=IXFF
								ILOX_ISHAPE(IKG2,L)=IRECEXR1(IXFF)%ISHAPE(IXSST)
								EXIT
							       END IF
							END DO
						END IF
					 END DO	
					SELECT CASE(ILOX_ISHAPE(IKG2,L))
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
					ILON_X(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,1)
					ILON_Y(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,2)
					ILON_Z(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,3)
					ILOX_XXC(IKG2,L)=sum(ILON_X(IKG2,L,1:IN_STEN))/RIN_STEN
					ILOX_YYC(IKG2,L)=sum(ILON_Y(IKG2,L,1:IN_STEN))/RIN_STEN
					ILOX_ZZC(IKG2,L)=sum(ILON_Z(IKG2,L,1:IN_STEN))/RIN_STEN
					
					
					
					
					
					
! 			
					
				END IF
				END DO
				END IF
				END IF
				END DO

END SUBROUTINE LOCALISE_STENCIL


SUBROUTINE LOCALISE_STENCIL2d(N,Iconsi,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::n,iconsi
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X
REAL,DIMENSION(1:DIMENSIONA)::CORDS
INTEGER::I,J,l,IXFF,IXSST,ikg,ismp,inv,ineedt,ikg2,IN_STEN,K,itarget,NOJCOUNT
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
				ILOX_IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				j=(XMPIL(ILOX_IHEXG(IKG2,L)))
				ILOX_IHEXL(IKG2,L)=J
				ILOX_ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
				CALL COMPUTE_CENTRE2d(j,cords)
				ILOX_XXC(IKG2,L)=CORDS(1)
				ILOX_YYC(IKG2,L)=CORDS(2)
				
				
				ILON_NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)
				DO K=1,ielem(n,j)%nonodes    
				ILON_x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1)
				ILON_y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2)
				
				END DO
				end do
				ELSE
				IKG2=IKG2+1
				DO L=1,itarget
				ILOX_IHEXG(IKG2,L)=ILOCALSTENCIL(N,I,ISMP,L)
				ILOX_IHEXB(IKG2,L)=XMPIE(ILOCALSTENCIL(N,I,ISMP,L))
				IF (ILOX_IHEXB(IKG2,L).EQ.N)THEN
				      j=(XMPIL(ILOX_IHEXG(IKG2,L)))
						IF (ILOX_IHEXG(IKG2,L).EQ.IELEM(N,J)%IHEXGL)THEN
						ILOX_IHEXL(IKG2,L)=J
						ILOX_ISHAPE(IKG2,L)=IELEM(N,J)%ISHAPE
						END IF
						CALL COMPUTE_CENTRE2d(j,cords)
						ILOX_XXC(IKG2,L)=CORDS(1)
						ILOX_YYC(IKG2,L)=CORDS(2)
						
						ILON_NODCOUNT(IKG2,L,1:ielem(n,j)%nonodes)=ielem(n,j)%nodes(1:ielem(n,j)%nonodes)
						DO K=1,ielem(n,j)%nonodes    
				ILON_x(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(1)
				ILON_y(IKG2,L,K)=inoder(ielem(n,j)%nodes(K))%cord(2)
				
				END DO
				ELSE
					 DO IXFF=1,INEEDT
						IF (IEXCORDR(IXFF)%PROCID.EQ.ILOX_IHEXB(IKG2,L))THEN
							DO IXSST=1,IRECEXR(IXFF)%MUCHINEED(1)
							      IF ((ILOX_IHEXG(IKG2,L)).EQ.IRECEXR1(IXFF)%WHATINEED(IXSST))THEN
								ILOX_IHEXL(IKG2,L)=IXSST
								ILOX_IHEXN(IKG2,L)=IXFF
								ILOX_ISHAPE(IKG2,L)=IRECEXR1(IXFF)%ISHAPE(IXSST)
								EXIT
							       END IF
							END DO
						END IF
					 END DO	
					SELECT CASE(ILOX_ISHAPE(IKG2,L))
					CASE(5)
					IN_STEN=4
					RIN_STEN=4.D0
					CASE(6)
					IN_STEN=3
					RIN_STEN=3.D0

					
					END SELECT
					ILON_X(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,1)
					ILON_Y(IKG2,L,1:IN_STEN)=IEXCORDR(ILOX_IHEXN(IKG2,L))%NODECORD(ILOX_IHEXL(IKG2,L),1:IN_STEN,2)
					
					ILOX_XXC(IKG2,L)=sum(ILON_X(IKG2,L,1:IN_STEN))/RIN_STEN
					ILOX_YYC(IKG2,L)=sum(ILON_Y(IKG2,L,1:IN_STEN))/RIN_STEN
					
				END IF
				END DO
				END IF
				END IF
				END DO

				

				

END SUBROUTINE LOCALISE_STENCIL2d




SUBROUTINE  LOCALISE_STEN2(N,ICONSI,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSI,N
INTEGER::I,J,K,L,KK,PRK,JJ,kmaxe,ineedt,jx2,jx,iivd,iivd3,facexx,IXXFFf,in1,itarget,idum,eltype,n_node,ELEM_DEC
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X
real,dimension(3)::tempcentres,TEMP_cG
real::dumv1,dumv2,detjc,dist1,MA,MB,MC,MD,ME,MF,MG,MH,MI,MDD,tempxx
REAL,dimension(1:dimensiona)::CORDS
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL,DIMENSION(1:dimensiona,1:dimensiona)::VVA1
REAL,DIMENSION(1)::DETA

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
	    IF(PER_ROT.EQ.0)THEN
	    IF(ABS(ILOX_XXC(JJ,J)-ILOX_XXC(JJ,1)).GT.XPER*oo2)THEN
		    ILOX_XXC(JJ,J)=ILOX_XXC(JJ,J)+(XPER*SIGN(1.0,ILOX_XXC(JJ,1)-XPER*oo2))
		    DO KK=1,8
		    ILON_X(JJ,J,KK)=ILON_X(JJ,J,KK)+(XPER*SIGN(1.0,ILOX_XXC(JJ,1)-XPER*oo2))
		    END DO
	    END IF
	    IF(ABS(ILOX_YYC(JJ,J)-ILOX_YYC(JJ,1)).GT.YPER*oo2)THEN
		    ILOX_YYC(JJ,J)=ILOX_YYC(JJ,J)+(YPER*SIGN(1.0,ILOX_YYC(JJ,1)-YPER*oo2))
		  
		    DO KK=1,8
		    ILON_Y(JJ,J,KK)=ILON_Y(JJ,J,KK)+(YPER*SIGN(1.0,ILOX_YYC(JJ,1)-YPER*oo2))
		    END DO
		    
	    END IF
	    IF(ABS(ILOX_ZZC(JJ,J)-ILOX_ZZC(JJ,1)).GT.ZPER*oo2)THEN
		    ILOX_ZZC(JJ,J)=ILOX_ZZC(JJ,J)+(ZPER*SIGN(1.0,ILOX_ZZC(JJ,1)-ZPER*oo2))
		      DO KK=1,8
		    ILON_Z(JJ,J,KK)=ILON_Z(JJ,J,KK)+(ZPER*SIGN(1.0,ILOX_ZZC(JJ,1)-ZPER*oo2))
		    END DO
		   
	    END IF
	    ELSE
            if (ILOX_PERIODICFLAG(jj,J).EQ.2) THEN
                    tempxx=ILOX_XXC(JJ,J)
                    ILOX_XXC(JJ,J)=tempxx*cos(angle_per)-ILOX_YYC(JJ,J)*sin(angle_per)
                    ILOX_YYC(JJ,J)=tempxx*sin(angle_per)+ILOX_YYC(JJ,J)*cos(angle_per)
                    !write(3300+n,'(6es14.6,I5)'),ILOX_XXC(JJ,1),ILOX_YYC(JJ,1),ILOX_ZZC(JJ,1),ILOX_XXC(JJ,J),ILOX_YYC(JJ,J),ILOX_ZZC(JJ,j),ILOX_PERIODICFLAG(jj,J)
                    DO KK=1,8
                    tempxx=ILON_X(JJ,J,KK)
                            ILON_X(JJ,J,KK)=tempxx*cos(angle_per)-ILON_Y(JJ,J,KK)*sin(angle_per)
                            ILON_Y(JJ,J,KK)=tempxx*sin(angle_per)+ILON_Y(JJ,J,KK)*cos(angle_per)
                    END DO
            end if
            if (ILOX_PERIODICFLAG(jj,J).EQ.1) THEN
                        tempxx=ILOX_XXC(JJ,J)
                        ILOX_XXC(JJ,J)=tempxx*cos(-angle_per)-ILOX_YYC(JJ,J)*sin(-angle_per)
                        ILOX_YYC(JJ,J)=tempxx*sin(-angle_per)+ILOX_YYC(JJ,J)*cos(-angle_per)
                        !write(3300+n,'(6es14.6,I5)'),ILOX_XXC(JJ,1),ILOX_YYC(JJ,1),ILOX_ZZC(JJ,1),ILOX_XXC(JJ,J),ILOX_YYC(JJ,J),ILOX_ZZC(JJ,j),ILOX_PERIODICFLAG(jj,J)
                DO KK=1,8
                        tempxx=ILON_X(JJ,J,KK)
                        ILON_X(JJ,J,KK)=tempxx*cos(-angle_per)-ILON_Y(JJ,J,KK)*sin(-angle_per)
                        ILON_Y(JJ,J,KK)=tempxx*sin(-angle_per)+ILON_Y(JJ,J,KK)*cos(-angle_per)
                END DO
            end if
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
	    NODES_LIST(k,1)=ILON_x(1,1,K)
	    NODES_LIST(k,2)=ILON_y(1,1,K)
	    NODES_LIST(k,3)=ILON_z(1,1,K)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  call DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
    
      SELECT CASE(ielem(n,i)%ishape)
      CASE(1)
      if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREHEXA(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      DUMV1=HEXAVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      CALL COMPUTE_CENTRE3d(i,cords)
      vext(1,1:dims)=cords(1:dims)
      else
	VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	  call COMPUTEJACOBIANS(N,VEXT,VVA1,DETA)
      end if	    
      CASE(2)
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	call COMPUTEJACOBIANS(N,VEXT,VVA1,DETA)
      CASE(3)
	VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	  call COMPUTEJACOBIANS(N,VEXT,VVA1,DETA)
      CASE(4)
       if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREPRISM(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      DUMV1=PRISMVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      VEXT(1,1:3)=VEXT(6,1:3)
      else
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
      call COMPUTEJACOBIANS(N,VEXT,VVA1,DETA)
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
					TEMPCENTRES(1)=ILON_X(JJ,L,KK)
					TEMPCENTRES(2)=ILON_Y(JJ,L,KK)
					TEMPCENTRES(3)=ILON_Z(JJ,L,KK)
					TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-TEMP_CG(:))
					ILON_X(JJ,L,KK)=TEMPCENTRES(1)
					ILON_Y(JJ,L,KK)=TEMPCENTRES(2)
					ILON_Z(JJ,L,KK)=TEMPCENTRES(3)
				END DO
				end if
			END DO	
		END DO	
		IDUM=0
		if (ielem(n,i)%interior.eq.1)then
                        DO j=1,IELEM(N,I)%IFCA
                        if (ielem(n,i)%ibounds(J).gt.0)then
                            if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
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
				ILOX_VOLUME(JJ,L)=(IELEM(N,ILOX_IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				
				
				
				
				if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				end if
				
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%PERIODICFLAG(JJ,L)=ILOX_PERIODICFLAG(JJ,L)
				else
				
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				
				end if
				
! 				
				Else
				IF (ILOX_IHEXB(JJ,L).eq.N)THEN
				ILOX_VOLUME(JJ,L)=(IELEM(N,ILOX_IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				end if
				
				
	
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%PERIODICFLAG(JJ,L)=ILOX_PERIODICFLAG(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ILOX_IHEXB(JJ,L)
 				else
                                ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ILOX_IHEXB(JJ,L)
				

                                end if

				else 
				if ((EES.ne.5).or.(jj.eq.1))then
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%PERIODICFLAG(JJ,L)=ILOX_PERIODICFLAG(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ILOX_IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXN(JJ,L)=ILOX_IHEXN(JJ,L)
				else
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ILOX_IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXNc(JJ,L)=ILOX_IHEXN(JJ,L)
				
				end if
				
				
				
				ELTYPE=ILOX_ISHAPE(jj,L)
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
					      NODES_LIST(k,1)=ILON_x(jj,l,K)
					      NODES_LIST(k,2)=ILON_y(jj,l,K)
					      NODES_LIST(k,3)=ILON_z(jj,l,K)
					      VEXT(K,:)=NODES_LIST(k,:)
					      
					    END DO
					    CALL DECOMPOSE3(n,eltype,NODES_LIST,ELEM_LISTD)
					    
					    dumv2=0.0d0
					    do k=1,ELEM_DEC
					    VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
					    dumv2=dumv2+TETRAVOLUME(N,vext)
					    
					    end do
					    
					    	
					
					    
								    
					    
					    
					    ILOX_VOLUME(JJ,L)=dumv2
				

				if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				end if
				
				!ILOCAL_RECON3(I)%VOLUME(JJ,L)=ILOX_VOLUME(JJ,L)
				
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				!ILOCAL_RECON3(I)%VOLUMEc(JJ,L)=ILOX_VOLUME(JJ,L)
				end if
				
				

				end if
				END IF
! 							
			END DO	
		END DO

		if (ielem(n,i)%interior.eq.0)then
	  CALL COMPUTE_CENTRE3d(i,cords)
	  vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
!  		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(2,1:DIMS)-VEXT(1,1:DIMS)
		    end if
	  end do
      else
		    CALL COMPUTE_CENTRE3d(i,cords)
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
		  call COMPUTE_CENTRE3dF(N,I,k,IXXFFf,cords)
		  VEXT(2,1:dims)=cords(1:dims)
	
		  
		      dist1=distance3(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
! 			IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(2,1:DIMS)-VEXT(1,1:DIMS)
		    end if
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
!  		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(2,1:DIMS)-VEXT(1,1:DIMS)
		    end if
		else						!from another cpu 
		    DO In1=1,ielem(n,i)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ILOX_XXC(1,In1);vext(2,2)=ILOX_yyC(1,In1); vext(2,3)=ILOX_zzC(1,In1)
		     dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
! 			IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(2,1:DIMS)-VEXT(1,1:DIMS)
				      end if
			  end if
		    end do
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(PER_ROT.EQ.0)THEN
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    IF(ABS(vext(2,3)-vext(1,3)).GT.zPER*oo2)THEN
		    vext(2,3)=vext(2,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
		    end if
		    ELSE
                if (ibound(n,ielem(n,i)%ibounds(k))%icode.eq.5) then
                    DO KK=1,n_node
                    tempxx=vext(kk,1)
                    vext(kk,1)=tempxx*cos(-angle_per)-sin(-angle_per)*vext(kk,2)
                    vext(kk,2)=tempxx*sin(-angle_per)+cos(-angle_per)*vext(kk,2)
                    END DO
                else
                    DO KK=1,n_node
                    tempxx=vext(kk,1)
                    vext(kk,1)=tempxx*cos(angle_per)-sin(angle_per)*vext(kk,2)
                    vext(kk,2)=tempxx*sin(angle_per)+cos(angle_per)*vext(kk,2)
                    END DO
                end if
		    END IF
		    dist1=distance3(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
!  		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(2,1:DIMS)-VEXT(1,1:DIMS)
		    end if
		else	!periodic boundaries from another cpu
		     DO In1=1,ielem(n,i)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ILOX_XXC(1,In1);vext(2,2)=ILOX_yyC(1,In1); vext(2,3)=ILOX_zzC(1,In1)
		    IF(PER_ROT.EQ.0)THEN 
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    IF(ABS(vext(2,3)-vext(1,3)).GT.zPER*oo2)THEN
		    vext(2,3)=vext(2,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
		    end if
		    ELSE
                if (ibound(n,ielem(n,i)%ibounds(k))%icode.eq.5) then
                    DO KK=1,n_node
                      tempxx=vext(kk,1)
				      vext(kk,1)=tempxx*cos(-angle_per)-sin(-angle_per)*vext(kk,2)
				      vext(kk,2)=tempxx*sin(-angle_per)+cos(-angle_per)*vext(kk,2)
				     END DO 
                else
                    DO KK=1,n_node
				      tempxx=vext(kk,1)
				      vext(kk,1)=tempxx*cos(angle_per)-sin(angle_per)*vext(kk,2)
				      vext(kk,2)=tempxx*sin(angle_per)+cos(angle_per)*vext(kk,2)
				    END DO  
                end if
		    END IF
		    dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
!  		    IELEM(N,i)%DIH2(K,1:DIMS)=VEXT(2,1:DIMS)-VEXT(1,1:DIMS)
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
				TEMPCENTRES(1)=ILOX_XXC(JJ,L)
				TEMPCENTRES(2)=ILOX_YYC(JJ,L)
				TEMPCENTRES(3)=ILOX_ZZC(JJ,L)
				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-TEMP_CG(:))
				ILOX_XXC(JJ,L)=TEMPCENTRES(1)
				ILOX_YYC(JJ,L)=TEMPCENTRES(2)
				ILOX_ZZC(JJ,L)=TEMPCENTRES(3)

				
				
				
				
				
				
				
				
				
                        
				end if
			END DO	
		END DO					

	

		
	
! 
END SUBROUTINE LOCALISE_STEN2



SUBROUTINE  LOCALISE_STEN2d(N,ICONSI,ILOX_IHEXG,ILOX_IHEXL,ILOX_IHEXB,ILOX_IHEXN,ILOX_ISHAPE,ILOX_XXC,ILOX_YYC,ILOX_ZZC,ILOX_VOLUME,ILOX_PERIODICFLAG,ILON_NODCOUNT,ILON_X,ILON_Y,ILON_Z)
!> @brief
!> This subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
IMPLICIT NONE
INTEGER::I,J,K,L,KK,PRK,JJ,kmaxe,ineedt,jx2,jx,in1,facexx,ixxfff,IHGT,IHGJ,ITARGET,IDUM,IN_STEN,NJ,ELEM_DEC,ELtype
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXG  !GLOBAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXL  !LOCAL INDEX OF CELLS
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXB  !CPU THAT THAT EACH CELL BELONGS TO
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_IHEXN  !INTERNAL INDEX FROM WHERE TO TAKE THE VALUES FROM COMMUNICATED MESSAGES
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ISHAPE !SHAPE OF EACH ELEMENT
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_XXC       !CELL CENTRE COORDINATES IN X
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_YYC       !CELL CENTRE COORDINATES IN Y
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_ZZC      !CELL CENTRE COORDINATES IN Z
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_VOLUME    !CELL VOLUME
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::ILOX_PERIODICFLAG
INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_NODCOUNT  !NUMBER OF NODES
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_X           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Y           !COORDINATES OF EACH NODE IN X
REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::ILON_Z           !COORDINATES OF EACH NODE IN X
real,dimension(2)::tempcentres,rel2
real::dumv1,dumv2,detjc,dist1,DISTFD,X1X,Y1Y,X2X,Y2Y
REAL,DIMENSION(1:DIMENSIONA)::CORDS
INTEGER,INTENT(IN)::ICONSI,N
INTEGER,DIMENSION(4)::NOJCOUNT
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:6,1:4,1:DIMENSIONA)::ELEM_LISTD
REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS)::QPOINTS
REAL,DIMENSION(1:NUMBEROFPOINTS)::WEQUA3D
REAL,DIMENSION(1:dimensiona,1:dimensiona)::VVA1
REAL,DIMENSION(1)::DETA

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
	    IF(ABS(ILOX_XXC(JJ,J)-ILOX_XXC(JJ,1)).GT.XPER*oo2)THEN
		    ILOX_XXC(JJ,J)=ILOX_XXC(JJ,J)+(XPER*SIGN(1.0,ILOX_XXC(JJ,1)-XPER*oo2))
		    DO KK=1,4
		    ILON_X(JJ,J,KK)=ILON_X(JJ,J,KK)+(XPER*SIGN(1.0,ILOX_XXC(JJ,1)-XPER*oo2))
		    END DO
	    END IF
	    IF(ABS(ILOX_YYC(JJ,J)-ILOX_YYC(JJ,1)).GT.YPER*oo2)THEN
		    ILOX_YYC(JJ,J)=ILOX_YYC(JJ,J)+(YPER*SIGN(1.0,ILOX_YYC(JJ,1)-YPER*oo2))
		  
		    DO KK=1,4
		    ILON_Y(JJ,J,KK)=ILON_Y(JJ,J,KK)+(YPER*SIGN(1.0,ILOX_YYC(JJ,1)-YPER*oo2))
		    END DO
		    
	    END IF
	    END DO
		
	END DO
	END IF
	

	!$OMP MASTER

	IF (CODE_PROFILE.EQ.30)THEN
		DO NJ=1,IELEM(N,I)%nonodes
				X1X=ILON_x(1,1,NJ)
				Y1Y=ILON_Y(1,1,NJ)
				NOJCOUNT(NJ)=0

				IF (ILOCAL_RECON3(I)%LOCAL.eq.1)then
					DO L=2,itarget
						j=(XMPIL(ILOX_IHEXG(1,L)))
						DO K=1,ielem(n,j)%nonodes    
							X2X=ILON_x(1,L,K)
							Y2Y=ILON_y(1,L,K)

							DISTFD=SQRT(((X1X-X2X)**2)+((Y1Y-Y2Y)**2))

							IF (DISTFD.LT.TOLSMALL)THEN
								NOJCOUNT(NJ)=NOJCOUNT(NJ)+1
								IELEM(N,I)%NODES_NEIGHBOURS(NJ,NOJCOUNT(NJ))=L
							END IF
						END DO
					END DO
				END IF
				!---------------------- MIXED----------!
				IF (ILOCAL_RECON3(I)%LOCAL.NE.1)then
					DO L=2,itarget
						IF (ILOX_IHEXB(1,L).EQ.N)THEN
							j=(XMPIL(ILOX_IHEXG(1,L)))
							DO K=1,ielem(n,j)%nonodes    
								X2X=ILON_x(1,L,K)
								Y2Y=ILON_y(1,L,K)

								DISTFD=SQRT(((X1X-X2X)**2)+((Y1Y-Y2Y)**2))

								IF (DISTFD.LT.TOLSMALL)THEN
									NOJCOUNT(NJ)=NOJCOUNT(NJ)+1
									IELEM(N,I)%NODES_NEIGHBOURS(NJ,NOJCOUNT(NJ))=L
								END IF
							END DO
						ELSE
							SELECT CASE(ILOX_ISHAPE(1,L))
							CASE(5)
							IN_STEN=4
							CASE(6)
							IN_STEN=3
							END SELECT
							DO K=1,IN_STEN
								X2X=ILON_X(1,L,K)
								Y2Y=ILON_Y(1,L,K)

								DISTFD=SQRT(((X1X-X2X)**2)+((Y1Y-Y2Y)**2))

								IF (DISTFD.LT.TOLSMALL)THEN
									NOJCOUNT(NJ)=NOJCOUNT(NJ)+1
									IELEM(N,I)%NODES_NEIGHBOURS(NJ,NOJCOUNT(NJ))=L
								END IF

							END DO
						END IF
					END DO
				END IF

			END DO
			WRITE(630+N,*)"ELEMENT NUMBER GLOBAL",IELEM(N,I)%IHEXGL
			ALLOCATE(IELEM(N,I)%NOJECOUNT(IELEM(N,I)%nonodes))
			
			DO NJ=1,IELEM(N,I)%nonodes
				WRITE(630+N,*)"NODE NUMBER",NJ
				IELEM(N,I)%NOJECOUNT(NJ)=NOJCOUNT(NJ)

				DO J=1,NOJCOUNT(NJ)
				!IF (IELEM(N,I)%NODES_NEIGHBOURS(NJ,J).GT.0)THEN
					WRITE(630+N,*)J,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)
				!END IF
				END DO
			END DO


	END IF
	!$OMP END MASTER




	


	  
     VEXT=0.0d0
    NODES_LIST=0.0d0
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=0.0d0
       
      jx=IELEM(N,I)%NONODES
	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,1)=ILON_x(1,1,K)
	    NODES_LIST(k,2)=ILON_y(1,1,K)
	   
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
    
      SELECT CASE(ielem(n,i)%ishape)

      CASE(5)
      if (IELEM(N,I)%MODE.eq.0)then
      CALL QUADRATUREQUAD(N,IGQRULES,VEXT,QPOINTS,WEQUA3D)
      DUMV1=QUADVOLUME(N,VEXT,QPOINTS,WEQUA3D)
      CALL COMPUTE_CENTRE2d(i,cords)
      vext(1,1:dims)=cords(1:dims)
      else
	VEXT(1:3,1:2)=ELEM_LISTD(1,1:3,1:2)
! 	  DUMV1=TRIANGLEVOLUME(N)
	  call COMPUTEJACOBIANS2(N,VEXT,VVA1,DETA)
	end if
      

      CASE(6)
      VEXT(1:3,1:2)=ELEM_LISTD(1,1:3,1:2)
! 	DUMV1=TRIANGLEVOLUME(N)
	call COMPUTEJACOBIANS2(N,VEXT,VVA1,DETA)

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
		    ILOCAL_RECON3(I)%INVCCJAC(1:2,1:2)=0.0d0
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
! 				TEMPCENTRES(1)=ILOX_XXC(JJ,L)
! 				TEMPCENTRES(2)=ILOX_YYC(JJ,L)
! 				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
! 				
! 				ILOX_XXC(JJ,L)=TEMPCENTRES(1)
! 				ILOX_YYC(JJ,L)=TEMPCENTRES(2)
				
				
				
				
				DO KK=1,4
					TEMPCENTRES=ZERO
					TEMPCENTRES(1)=ILON_X(JJ,L,KK)
					TEMPCENTRES(2)=ILON_Y(JJ,L,KK)
					
					
					TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
					
					ILON_X(JJ,L,KK)=TEMPCENTRES(1)
					ILON_Y(JJ,L,KK)=TEMPCENTRES(2)
					
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
				ILOX_VOLUME(JJ,L)=(IELEM(N,ILOX_IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				 if ((EES.ne.5).or.(jj.eq.1))then
				 if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				end if
				
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				
				
				end if
				Else
				IF (ILOX_IHEXB(JJ,L).eq.N)THEN
				ILOX_VOLUME(JJ,L)=(IELEM(N,ILOX_IHEXL(JJ,L))%TOTVOLUME)/ABS(detjc)
				 if ((EES.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				end if
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ILOX_IHEXB(JJ,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ILOX_IHEXB(JJ,L)
				end if
				else 
				 if ((EES.ne.5).or.(jj.eq.1))then
				ILOCAL_RECON3(I)%IHEXG(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXL(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXB(JJ,L)=ILOX_IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXN(JJ,L)=ILOX_IHEXN(JJ,L)
				else
				ILOCAL_RECON3(I)%IHEXGc(JJ,L)=ILOX_IHEXG(JJ,L)
				ILOCAL_RECON3(I)%IHEXLc(JJ,L)=ILOX_IHEXL(JJ,L)
				ILOCAL_RECON3(I)%IHEXBc(JJ,L)=ILOX_IHEXB(JJ,L)
				ILOCAL_RECON3(I)%IHEXNc(JJ,L)=ILOX_IHEXN(JJ,L)
				
				
				end if
				
				      ELTYPE=ILOX_ISHAPE(jj,L)
				      ELEM_LISTD=0.0d0; VEXT=0.0d0; NODES_LIST=0.0d0
				      select case(ELTYPE)
				      
				      case(5)
				      ELEM_DEC=2; jx=4
				      case(6)
				    ELEM_DEC=1; jx=3
				      
				      end select

					    do K=1,jx
					      NODES_LIST(k,1)=ILON_x(jj,l,K)
					      NODES_LIST(k,2)=ILON_y(jj,l,K)
					     
					      VEXT(K,:)=NODES_LIST(k,:)
					    END DO
					    CALL DECOMPOSE2(n,eltype,NODES_LIST,ELEM_LISTD)
					    dumv2=0.0d0
					    do k=1,ELEM_DEC
					    VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
					    dumv2=dumv2+TRIANGLEVOLUME(N,vext)
					    end do
					    ILOX_VOLUME(JJ,L)=dumv2
					     if ((EES.ne.5).or.(jj.eq.1))then
						 if (idum.eq.1)then
				ILOCAL_RECON3(I)%VOLUME(1,L)=ILOX_VOLUME(1,L)
				else
				ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
				end if
				
                                            else
                                  ILOCAL_RECON3(I)%VOLUME(1,1)=ILOX_VOLUME(1,1)
                                            
                                            end if
! 				
				end if
				END IF
			END DO	
		END DO

	      if (ielem(n,i)%interior.eq.0)then
	  CALL COMPUTE_CENTRE2d(i,cords)
	  vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
		    end if
	  end do
      else
		    CALL COMPUTE_CENTRE2d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		  
		  IXXFFf=2
		  
		  
		  call COMPUTE_CENTRE2dF(N,iconsi,facexx,IXXFFf,CORDS)
		  VEXT(2,1:dims)=cords(1:dims)
		  
		      dist1=distance2(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
		    end if
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
		    end if
		else						!from another cpu 
		    DO In1=1,IELEM(N,I)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ILOX_XXC(1,In1);vext(2,2)=ILOX_yyC(1,In1)
		     dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
				      end if
			  end if
		    end do
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    
		    dist1=distance2(n,vext)
		    IF (RUNGEKUTTA.ge.2)THEN
		    IELEM(N,i)%DIH(K)=dist1
		    end if
		else	!periodic boundaries from another cpu
		     DO In1=1,IELEM(N,I)%iNUMNEIGHBOURS
			  IF (IELEM(N,i)%INEIGHG(K).EQ.ILOCAL_RECON3(i)%IHEXG(1,In1))THEN
				  IELEM(N,i)%INDEXI(K)=In1
				      IF (RUNGEKUTTA.ge.2)THEN
		    vext(2,1)=ILOX_XXC(1,In1);vext(2,2)=ILOX_yyC(1,In1);
		     
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		    end if
		    
		    dist1=distance2(n,vext)
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
				TEMPCENTRES(1)=ILOX_XXC(JJ,L)
				TEMPCENTRES(2)=ILOX_YYC(JJ,L)
				TEMPCENTRES(:)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),TEMPCENTRES(:)-VEXT(1,:))
				
				ILOX_XXC(JJ,L)=TEMPCENTRES(1)
				ILOX_YYC(JJ,L)=TEMPCENTRES(2)
				
				
				
				
				
				
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
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA)::CORDS
real::dist1

kmaxe=xmpielrank(n)

!$OMP DO
do i=1,kmaxe
	if (ielem(n,i)%interior.eq.0)then
		    CALL COMPUTE_CENTRE3d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
      
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
	  end do
	
	else
		    CALL COMPUTE_CENTRE3d(i,cords)
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
		  call COMPUTE_CENTRE3dF(N,I,k,IXXFFf,cords)
		  VEXT(2,1:dims)=cords(1:dims)

		  
		      dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		else						!from another cpu 
		
		    vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims)
		     dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE3d(j,cords)
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
		    dist1=distance3(n,vext)
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
		    dist1=distance3(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
	  end do
	end if

end do
!$OMP END DO



end subroutine direct_side



subroutine direct_side2d(n)
!> @brief
!> This subroutine establishes the distance betwen cell centres for each edge
implicit none
integer,intent(in)::n
integer::i,j,k,kmaxe,facexx,ixxfff
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA)::CORDS
real::dist1

kmaxe=xmpielrank(n)


!$OMP DO
do i=1,kmaxe
	if (ielem(n,i)%interior.eq.0)then
		    CALL COMPUTE_CENTRE2d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
      
	  do k=1,ielem(n,i)%ifca
		  j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
	  end do
	
	else
		    CALL COMPUTE_CENTRE2d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem(n,i)%ifca
		if (ielem(n,i)%ineighg(k).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		 
		  IXXFFf=2
		 
		  call COMPUTE_CENTRE2dF(N,I,facexx,IXXFFf,cords)
		  VEXT(2,1:dims)=cords(1:dims)
		  
		      dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1*2.0d0
		 end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).eq.0))then	!non periodic boundaries 
		if (ielem(n,i)%ineighb(k).eq.n)then		!within my cpu
		 j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		else						!from another cpu 
		    vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims)
		     dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
		if ((ielem(n,i)%ineighg(k).gt.0).and.(ielem(n,i)%ibounds(k).gt.0))then	!periodic boundaries within my cpu
		if (ielem(n,i)%ineighb(k).eq.n)then	
		     j=IELEM(N,i)%INEIGH(K)
		  CALL COMPUTE_CENTRE2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER/2.d0)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0D0,vext(1,1)-XPER/2.D0))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER/2.d0)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0D0,vext(1,2)-yPER/2.D0))
		    end if
		    
		    dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		
		else	!periodic boundaries from another cpu

		     vext(2,1:dims)=SOLCHANGER(IELEM(N,I)%INEIGHN(k))%CENTRES(IELEM(N,i)%Q_FACE(k)%Q_MAPL(1),1:dims) 
		    IF(ABS(vext(2,1)-vext(1,1)).GT.XPER*oo2)THEN
		    vext(2,1)=vext(2,1)+(XPER*SIGN(1.0D0,vext(1,1)-XPER/2.0D0))
		    end if
		    IF(ABS(vext(2,2)-vext(1,2)).GT.yPER*oo2)THEN
		    vext(2,2)=vext(2,2)+(yPER*SIGN(1.0D0,vext(1,2)-yPER/2.0D0))
		    end if
		    
		    dist1=distance2(n,vext)
		    IELEM(N,i)%DIH(K)=dist1
		end if
		end if
	  end do
	end if

end do
!$OMP END DO



end subroutine direct_side2d


SUBROUTINE GRADS_ASSIGN(N)
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,jj,KMAXE

KMAXE=XMPIELRANK(N)

IF (DIMENSIONA.EQ.3)THEN

!$OMP DO
do i=1,kmaxe
	CALL CHECKGRADS(N,I)
end do
!$OMP END DO

Else

!$OMP DO
do i=1,kmaxe
	CALL CHECKGRADS2D(N,I)
end do
!$OMP END DO

END IF




END SUBROUTINE GRADS_ASSIGN


SUBROUTINE CHECKGRADS(N,ICONSI)
!> @brief
!> This subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSI
REAL::DXX1,dxx2,TEMPG1,dist1,dist2,oo2,surfmin,surfmax
INTEGER::I,J,K,L,jj,icount3,nnd,ixf4,IDC,IDC2
REAL,dimension(1:dimensiona)::CORDS
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST

i=iconsi
    IELEM(N,I)%GGS=greengo
    dxx1=-TOLBIG; dxx2=TOLBIG
    CALL COMPUTE_CENTRE3d(i,cords)
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
					call COMPUTE_CENTRE3dF(N,I,l,nnd,cords)
				VEXT(2,1:dims)=cords(1:dims)

				end if
		
		
		
	      surfmin=min(surfmin,IELEM(N,I)%SURF(L))
	      surfmax=max(surfmax,IELEM(N,I)%SURF(L))
	      
	      
	      
	      
		    
	


		  dist1=distance3(n,vext)
		if (dist1.lt.dxx2)then
		  dxx2=dist1
		end if
		if (dist1.gt.dxx1)then
		  dxx1=dist1
		end if
end do
        ielem(n,i)%condition=1.0


        ielem(n,i)%condition=surfmax/surfmin

		if (CODE_PROFILE.eq.9)then
        if ((ielem(n,i)%condition.gt.30).or.(ielem(n,i)%ishape.eq.4))then
        ielem(n,i)%full=0
		end if
		end if

        !if (turbulence.gt.0)then
        

        !if (ielem(n,i)%ishape.eq.2)then
        !ielem(n,i)%condition=surfmax/surfmin
        
        !if (ielem(n,i)%condition.gt.30)then
        ! ielem(n,i)%hybrid=1
        !end if
        !end if
        
        !if (ielem(n,i)%ishape.eq.3)then
        !ielem(n,i)%condition=surfmax/surfmin
        
        !if (ielem(n,i)%condition.gt.10)then
         !ielem(n,i)%hybrid=1
        !end if
        
       
        
        !end if


        !end if
       



	    TEMPG1=ielem(n,i)%condition!MAX((DXX1/DXX2),(DXX2/DXX1))
! 	   
			ielem(n,i)%erx=tempg1
	    

              IF (CODE_PROFILE.EQ.88)THEN
                      IF ((IELEM(N,I)%ISHAPE.EQ.3)) THEN
                                        IELEM(N,I)%FULL=0
                         END IF

              END IF


	      IF (TEMPG1.GT.GRIDAR1)THEN
	      IELEM(N,I)%GGS=1
	      IF ((IADAPT.EQ.1).or.(code_profile.eq.88).or.(code_profile.eq.98))THEN
                IELEM(N,I)%FULL=0
			END IF
	      end if 
   if (fastest.eq.0)then
	       dxx1=-tolbig; dxx2=tolbig
	       JJ=1
	       DO L=1,ielem(n,i)%iNUMNEIGHBOURS
		       if (ILOCAL_RECON3(i)%VOLUME(JJ,L).lt.dxx2)then
		      
		       dxx2=ILOCAL_RECON3(i)%VOLUME(JJ,L)
		       end if
		       if (ILOCAL_RECON3(i)%VOLUME(JJ,L).gt.dxx1)then
		      
		       dxx1=ILOCAL_RECON3(i)%VOLUME(JJ,L)
		       end if
	       end do
	 TEMPG1=MAX((DXX1/DXX2),(DXX2/DXX1))
	 !ielem(n,i)%walldist=tempg1
! 	     IF (TEMPG1.GT.GRIDAR2)THEN
! 	       IELEM(N,I)%GGS=1
! 	  	IF ((IADAPT.EQ.1).or.(code_profile.eq.88).or.(code_profile.eq.98))THEN
!                 IELEM(N,I)%FULL=0
! 			END IF
! 	       end if
	      
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
REAL,dimension(1:dimensiona)::CORDS
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT


do iconsi=1,xmpielrank(n)

i=iconsi
    IELEM(N,I)%GGS=greengo
    dxx1=-TOLBIG; dxx2=TOLBIG
    CALL COMPUTE_CENTRE3d(i,cords)
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
REAL::DXX1,dxx2,TEMPG1,dist1,dist2,oo2,surfmin,surfmax
INTEGER::I,J,K,L,jj,icount3,nnd,ixf4
REAL,dimension(1:dimensiona)::CORDS
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST

i=iconsi

tempg1=0.0; 
    IELEM(N,I)%GGS=greengo
    dxx1=-TOLBIG; dxx2=TOLBIG
    CALL COMPUTE_CENTRE2d(i,CORDS)
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




  surfmin=1.0e16
surfmax=1.0e-16

     
      do l=1,ielem(n,i)%ifca
      
	      
	      nnd=2
	     
			 surfmin=min(surfmin,IELEM(N,I)%SURF(L))
	      surfmax=max(surfmax,IELEM(N,I)%SURF(L))
		
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
		      call COMPUTE_CENTRE2dF(N,I,l,nnd,CORDS)
		      VEXT(2,1:dims)=cords(1:dims)
		    end if
		  else
		
		    
		      call COMPUTE_CENTRE2dF(N,I,l,nnd,CORDS)
		  VEXT(2,1:dims)=cords(1:dims)
		  end if
		else
		       call COMPUTE_CENTRE2dF(N,I,l,nnd,CORDS)
		  VEXT(2,1:dims)=cords(1:dims)
		
		end if
	      
		    
	


		  dist1=distance2(n,VEXT)
		if (dist1.lt.dxx2)then
		  dxx2=dist1
		end if
		if (dist1.gt.dxx1)then
		  dxx1=dist1
		end if
		
		
		
end do


			ielem(n,i)%condition=1.0


        ielem(n,i)%condition=surfmax/surfmin

        TEMPG1=ielem(n,i)%condition

	    IF (TEMPG1.GT.GRIDAR1)THEN
	      IELEM(N,I)%GGS=1

	      IF ((IADAPT.EQ.1).or.(code_profile.eq.88).or.(code_profile.eq.98))THEN
                IELEM(N,I)%FULL=0
			END IF

	      end if 
	      
	      
	      
   if (fastest.eq.0)then
	       dxx1=-tolbig; dxx2=tolbig
	       JJ=1
	       DO L=1,ielem(n,i)%iNUMNEIGHBOURS
		        if (ILOCAL_RECON3(i)%VOLUME(JJ,L).lt.dxx2)then

		       dxx2=ILOCAL_RECON3(i)%VOLUME(JJ,L)
		       end if
		       if (ILOCAL_RECON3(i)%VOLUME(JJ,L).gt.dxx1)then

		       dxx1=ILOCAL_RECON3(i)%VOLUME(JJ,L)
		       end if
	       end do
! 	 TEMPG1=MAX((DXX1/DXX2),(DXX2/DXX1))
! 	     IF (TEMPG1.GT.GRIDAR2)THEN
! 	       IELEM(N,I)%GGS=1
! 	       IF ((IADAPT.EQ.1).or.(code_profile.eq.88).or.(code_profile.eq.98))THEN
!                 IELEM(N,I)%FULL=0
! 			END IF
! 	       end if
end if


END SUBROUTINE CHECKGRADS2d







END MODULE LOCAL
