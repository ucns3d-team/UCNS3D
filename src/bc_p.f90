MODULE BOUNDARY
!> @brief
!> This module includes the subroutines for reading and establishing the boundary conditions
!> for all the cells that need to be bounded including the periodic ones
USE DECLARATION
USE LIBRARY
USE TRANSFORM
IMPLICIT NONE
CONTAINS



SUBROUTINE READ_BOUND(N,IMAXB,IBOUND,XMPIELRANK)
!> @brief
!> This subroutine reads the boundary conditions for all the bounded elements

IMPLICIT NONE
INTEGER,INTENT(IN)::N,IMAXB
TYPE(A_BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND
CHARACTER(LEN=12)::BNDFILE,ibx_code
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
INTEGER::I,J,JI,K,LM,IEX,KMAXN,KK,KMAXE,kkk,JJJ,jfx,JB,KXK,ITR1,jj,JJ1
INTEGER::IOY,IBID,IB1,IB2,IB3,IB4,IBX1,ITL,ibgw,ibgw2,IBLEED,ibdum
integer,dimension(4)::ib_n
INTEGER::NB1,NB2,NB3,NB4
	KMAXE=XMPIELRANK(N)

kkk=0
totiw=0
ibgw=0
ibgw2=0






IF (DIMENSIONA.EQ.3)THEN

BNDFILE='GRID.bnd'
		if (binio.eq.0)OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
		if (binio.eq.1)OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)

itl=0

IF (BINIO.EQ.0)THEN
DO ji=1,IMAXB
	READ(10,*)IBID,IB1,IB2,IB3,IB4
	
	if ((inoder(IB1)%itor.gt.0).and.(inoder(IB2)%itor.GT.0).and.(inoder(IB3)%itor.GT.0)&
	.and.(inoder(IB4)%itor.GT.0))then
	itl=itl+1

	end if
end do
ELSE
DO ji=1,IMAXB
	READ(10)IBID,IB1,IB2,IB3,IB4
	
	if ((inoder(IB1)%itor.gt.0).and.(inoder(IB2)%itor.GT.0).and.(inoder(IB3)%itor.GT.0)&
	.and.(inoder(IB4)%itor.GT.0))then
	itl=itl+1

	end if
end do

END IF

if (itl.gt.0)then

allocate(ibound(n:n,itl))
ibound(n:n,:)%inum=0
end if


 close(10)

itl=0
if (binio.eq.0)OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
if (binio.eq.1)OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)

IF (BINIO.EQ.0)THEN
DO ji=1,IMAXB
	READ(10,*)IBID,IB_n(1),IB_n(2),IB_n(3),IB_n(4),ibx1
	
	if (ibx1.eq.4)then
	ibgw=ibgw+1
	end if
	
	if ((inoder(IB_n(1))%itor.gt.0).and.(inoder(IB_n(2))%itor.GT.0).and.(inoder(IB_n(3))%itor.GT.0)&
	.and.(inoder(IB_n(4))%itor.GT.0))then
	itl=itl+1
	    if (ibx1.eq.4)then
	    ibound(n,itl)%inum=ibgw
	    totiw=totiw+1
	    end if
		


	    ibound(n,itl)%icode=ibx1
	    ibound(n,itl)%IBID=IBID
	    if (ib_n(3).eq.ib_n(4))then
	    
	      ibound(n,itl)%ishape=6
	      ALLOCATE( ibound(n,itl)%ibl(1:3))
			ibound(n,itl)%ibl(1:3)=ib_n(1:3)
	    else

	      ibound(n,itl)%ishape=5
	       ALLOCATE( ibound(n,itl)%ibl(1:4))
			ibound(n,itl)%ibl(1:4)=ib_n(1:4)
	    end if

	end if
end do


ELSE
	
DO ji=1,IMAXB
	READ(10)IBID,IB_n(1),IB_n(2),IB_n(3),IB_n(4),ibx1
		
	

	if (ibx1.eq.4)then
	ibgw=ibgw+1
	end if
	
	if ((inoder(IB_n(1))%itor.gt.0).and.(inoder(IB_n(2))%itor.GT.0).and.(inoder(IB_n(3))%itor.GT.0)&
	.and.(inoder(IB_n(4))%itor.GT.0))then
	itl=itl+1
	    if (ibx1.eq.4)then
	    ibound(n,itl)%inum=ibgw
	    totiw=totiw+1
	    end if

	    ibound(n,itl)%icode=ibx1
	    ibound(n,itl)%ibid=IBID
	    if (ib_n(3).eq.ib_n(4))then
	    
	      ibound(n,itl)%ishape=6
	      ALLOCATE( ibound(n,itl)%ibl(1:3))
			ibound(n,itl)%ibl(1:3)=ib_n(1:3)
	    else

	      ibound(n,itl)%ishape=5
	       ALLOCATE( ibound(n,itl)%ibl(1:4))
			ibound(n,itl)%ibl(1:4)=ib_n(1:4)
	    end if

	end if
end do

END IF
 close(10)
n_boundaries=itl

itl=0
DO JI=1,n_boundaries
  IF ((IBOUND(N,JI)%ICODE.EQ.5).or.(IBOUND(N,JI)%ICODE.EQ.50))THEN
    ALLOCATE(IBOUND(N,JI)%LOCALN(2))
    allocate(ibound(n,ji)%CPUN(2))
    itl=itl+1
    ibound(n,ji)%localn=0
    ibound(n,ji)%cpun=0
  end if
end do


DO I=1,KMAXE
  IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
  ALLOCATE(IELEM(N,I)%iBOUNDs(IELEM(N,I)%IFCA))
  IELEM(N,I)%iBOUNDs=0
  ielem(n,I)%nofbc=0
  END IF
END DO

itl=0
JJ1=0
   do ji=1,n_boundaries
	do jfx=1,inoder2(ibound(n,ji)%ibl(1))%NUMBEROFNEIB
	     j=inoder2(ibound(n,ji)%ibl(1))%NEIBIDS(jfx)
	      if (ielem(n,j)%interior.eq.1)then
		  do jj=1,ielem(n,j)%ifca
			if (ielem(n,j)%ineighg(jj).eq.0)then
		      IF (IELEM(N,j)%TYPES_fACES(Jj).EQ.IBOUND(N,ji)%ISHAPE)THEN
				IF (IBOUND(N,ji)%ISHAPE.EQ.6)THEN
				  NB1=IELEM(N,j)%NODES_FACES(jJ,1)
				  NB2=IELEM(N,j)%NODES_FACES(jJ,2)
				  NB3=IELEM(N,j)%NODES_FACES(jJ,3)
				  IB1=IBOUND(N,ji)%IBL(1)
				  IB2=IBOUND(N,ji)%IBL(2)
				  IB3=IBOUND(N,ji)%IBL(3)
				  IF (((NB1.EQ.IB1).OR.(NB1.EQ.IB2).OR.(NB1.EQ.IB3)).AND.&
				  ((NB2.EQ.IB1).OR.(NB2.EQ.IB2).OR.(NB2.EQ.IB3)).AND.&
				  ((NB3.EQ.IB1).OR.(NB3.EQ.IB2).OR.(NB3.EQ.IB3)))THEN
				 IELEM(N,J)%IBOUNDS(JJ)=ji
				 IBOUND(N,Ji)%which=j
				 IBOUND(N,Ji)%face=jj
				  ielem(n,J)%nofbc=ielem(n,J)%nofbc+1
				    if ((IBOUND(N,JI)%ICODE.EQ.5).or.(IBOUND(N,JI)%ICODE.EQ.50))then
				    IBOUND(N,Ji)%localn(1)=J;IBOUND(N,Ji)%cpun(1)=n
				    JJ1=JJ1+1
				      
				    go to 51
				    end if
				  END IF
			      END IF
				  IF (IBOUND(N,ji)%ISHAPE.EQ.5)THEN
				    NB1=IELEM(N,j)%NODES_FACES(jJ,1)
				    NB2=IELEM(N,j)%NODES_FACES(jJ,2)
				    NB3=IELEM(N,j)%NODES_FACES(jJ,3)
				    NB4=IELEM(N,j)%NODES_FACES(jJ,4)
				    IB1=IBOUND(N,ji)%IBL(1)
				    IB2=IBOUND(N,ji)%IBL(2)
				    IB3=IBOUND(N,ji)%IBL(3)
				    IB4=IBOUND(N,ji)%IBL(4)
					IF (((NB1.EQ.IB1).OR.(NB1.EQ.IB2).OR.(NB1.EQ.IB3).OR.(NB1.EQ.IB4)).AND.&
					((NB2.EQ.IB1).OR.(NB2.EQ.IB2).OR.(NB2.EQ.IB3).OR.(NB2.EQ.IB4)).AND.&
					((NB3.EQ.IB1).OR.(NB3.EQ.IB2).OR.(NB3.EQ.IB3).OR.(NB3.EQ.IB4)).AND.&
					((NB4.EQ.IB1).OR.(NB4.EQ.IB2).OR.(NB4.EQ.IB3).OR.(NB4.EQ.IB4)))THEN
					IELEM(N,J)%IBOUNDS(JJ)=ji
					IBOUND(N,Ji)%which=j
				 IBOUND(N,Ji)%face=jj
				    ielem(n,J)%nofbc=ielem(n,J)%nofbc+1
				    if ((IBOUND(N,JI)%ICODE.EQ.5).or.(IBOUND(N,JI)%ICODE.EQ.50))then
				    IBOUND(N,Ji)%localn(1)=J;IBOUND(N,Ji)%cpun(1)=n
				    JJ1=JJ1+1
				      
				    go to 51
				    end if
					END IF
				    END IF
			  END IF
		      end if
		END DO
		
	      END IF
	  END DO
      51 continue
    END DO



     


end if

IF (DIMENSIONA.EQ.2)THEN

BNDFILE='GRID.bnd'
		IF (BINIO.EQ.0)OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
		IF (BINIO.EQ.1)OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)

itl=0

IF (BINIO.EQ.0)THEN
DO ji=1,IMAXB
	READ(10,*)IBID,IB1,IB2
	
	if ((inoder(IB1)%itor.gt.0).and.(inoder(IB2)%itor.GT.0))then
	itl=itl+1

	end if
end do
ELSE
DO ji=1,IMAXB
	READ(10)IBID,IB1,IB2
	
	if ((inoder(IB1)%itor.gt.0).and.(inoder(IB2)%itor.GT.0))then
	itl=itl+1

	end if
end do

END IF

if (itl.gt.0)then

allocate(ibound(n:n,itl))
ibound(n:n,:)%inum=0

end if


close(10)

itl=0
IF (BINIO.EQ.0)OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
IF (BINIO.EQ.1)OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)

IF (BINIO.EQ.0)THEN
DO ji=1,IMAXB
	READ(10,*)IBID,IB_n(1),IB_n(2),IB_n(3),IB_n(4),ibx1
	IF (BLEED.EQ.1)THEN
		ibdum=ibx1
		if (ibx1.eq.4)then

		if ((inoder(IB_n(1))%itor.gt.0).and.(inoder(IB_n(2))%itor.GT.0))then

		!NOW CHECK IF THIS IS IN A BLEED ZONE


		DO IBLEED=1,BLEED_NUMBER
			IF (((inoder(IB_n(1))%CORD(1).GE.bleed_start(IBLEED,1)).AND.(inoder(IB_n(1))%CORD(1).LE.bleed_END(IBLEED,1))).AND.((inoder(IB_n(2))%CORD(1).GE.bleed_start(IBLEED,1)).AND.(inoder(IB_n(2))%CORD(1).LE.bleed_END(IBLEED,1))).AND.((inoder(IB_n(1))%CORD(2).GE.bleed_start(IBLEED,2)).AND.(inoder(IB_n(1))%CORD(2).LE.bleed_END(IBLEED,2))).AND.((inoder(IB_n(2))%CORD(2).GE.bleed_start(IBLEED,2)).AND.(inoder(IB_n(2))%CORD(2).LE.bleed_END(IBLEED,2))))THEN
			IBdum=99


			END IF

		END DO

		END IF
		end if
		ibx1=IBdum


	END IF


	if ((ibx1.eq.4).or.(ibx1.eq.99))then
	ibgw=ibgw+1
	end if
	if ((inoder(IB_n(1))%itor.gt.0).and.(inoder(IB_n(2))%itor.GT.0))then
	itl=itl+1
	     if ((ibx1.eq.4).or.(ibx1.eq.99))then
	    ibound(n,itl)%inum=ibgw;totiw=totiw+1
	    end if
	    ibound(n,itl)%icode=ibx1
	   ibound(n,itl)%ibid=IBID

	      ibound(n,itl)%ishape=7
	       ALLOCATE( ibound(n,itl)%ibl(1:2))
			ibound(n,itl)%ibl(1:2)=ib_n(1:2)
	   

	end if
end do
ELSE
DO ji=1,IMAXB
	READ(10)IBID,IB_n(1),IB_n(2),IB_n(3),IB_n(4),ibx1
	IF (BLEED.EQ.1)THEN
		ibdum=ibx1
		if (ibx1.eq.4)then
		if ((inoder(IB_n(1))%itor.gt.0).and.(inoder(IB_n(2))%itor.GT.0))then

		!NOW CHECK IF THIS IS IN A BLEED ZONE


		DO IBLEED=1,BLEED_NUMBER
			IF (((inoder(IB_n(1))%CORD(1).GE.bleed_start(IBLEED,1)).AND.(inoder(IB_n(1))%CORD(1).LE.bleed_END(IBLEED,1))).AND.((inoder(IB_n(2))%CORD(1).GE.bleed_start(IBLEED,1)).AND.(inoder(IB_n(2))%CORD(1).LE.bleed_END(IBLEED,1))).AND.((inoder(IB_n(1))%CORD(2).GE.bleed_start(IBLEED,2)).AND.(inoder(IB_n(1))%CORD(2).LE.bleed_END(IBLEED,2))).AND.((inoder(IB_n(2))%CORD(2).GE.bleed_start(IBLEED,2)).AND.(inoder(IB_n(2))%CORD(2).LE.bleed_END(IBLEED,2))))THEN
			IBdum=99


			END IF

		END DO

		END IF
		end if
		ibx1=IBdum


	END IF


	if ((ibx1.eq.4).or.(ibx1.eq.99))then
	ibgw=ibgw+1
	end if
	if ((inoder(IB_n(1))%itor.gt.0).and.(inoder(IB_n(2))%itor.GT.0))then
	itl=itl+1
	     if ((ibx1.eq.4).or.(ibx1.eq.99))then
	    ibound(n,itl)%inum=ibgw;totiw=totiw+1
	    end if
	    ibound(n,itl)%icode=ibx1
	 ibound(n,itl)%ibid=IBID

	      ibound(n,itl)%ishape=7
	       ALLOCATE( ibound(n,itl)%ibl(1:2))
			ibound(n,itl)%ibl(1:2)=ib_n(1:2)
	   

	end if
end do

END IF

close(10)

n_boundaries=itl




DO JI=1,n_boundaries
  IF ((IBOUND(N,JI)%ICODE.EQ.5).or.(IBOUND(N,JI)%ICODE.EQ.50))THEN
    ALLOCATE(IBOUND(N,JI)%LOCALN(2))
    allocate(ibound(n,ji)%CPUN(2))
    ibound(n,ji)%localn=0
    ibound(n,ji)%cpun=0
  end if
end do

DO I=1,KMAXE
  IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
  ALLOCATE(IELEM(N,I)%iBOUNDs(IELEM(N,I)%IFCA))

IF (BLEED.EQ.1)THEN
  ALLOCATE(IELEM(N,I)%BLEEDN(IELEM(N,I)%IFCA))
  IELEM(N,I)%BLEEDN(:)=0
  END IF
  IELEM(N,I)%iBOUNDs=0
  ielem(n,i)%nofbc=0
  END IF
END DO

itl=0

   do ji=1,n_boundaries
	do jfx=1,inoder2(ibound(n,ji)%ibl(1))%NUMBEROFNEIB
	     j=inoder2(ibound(n,ji)%ibl(1))%NEIBIDS(jfx)
	      if (ielem(n,j)%interior.eq.1)then
		  do jj=1,ielem(n,j)%ifca
		     if (ielem(n,j)%ineighg(jj).eq.0)then
				
				  NB1=IELEM(N,j)%NODES_FACES(jJ,1)
				  NB2=IELEM(N,j)%NODES_FACES(jJ,2)
				  
				  IB1=IBOUND(N,ji)%IBL(1)
				  IB2=IBOUND(N,ji)%IBL(2)
				  
				   IF (((NB1.EQ.IB1).OR.(NB1.EQ.IB2)).AND.&
				((NB2.EQ.IB1).OR.(NB2.EQ.IB2)))THEN
				 IELEM(N,J)%IBOUNDS(JJ)=ji
				 IBOUND(N,Ji)%which=j
				 IBOUND(N,Ji)%face=jj
				    ielem(n,J)%nofbc=ielem(n,J)%nofbc+1
				    if ((IBOUND(N,JI)%ICODE.EQ.5).or.(IBOUND(N,JI)%ICODE.EQ.50))then
				    IBOUND(N,Ji)%localn(1)=J;IBOUND(N,Ji)%cpun(1)=n
				   
				    itl=itl+1
				    end if
				  END IF
			      
			  end if
		END DO
	      END IF
	  END DO
    END DO




 


  end if



totwalls=ibgw





      IF (TOTIW.GT.0)THEN
	ALLOCATE(IBOUND_T(TOTIW))
	ALLOCATE(IBOUND_T2(TOTIW))
	TOTIW=0
	
	
	
	DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	     if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
		 TOTIW=TOTIW+1
				IBOUND_T(TOTIW)=I
				IBOUND_T2(TOTIW)=j
	      END IF
	  end if
	END DO
   end if
END DO
	
	end if



	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 	DEALLOCATE(IBOUND)
END SUBROUTINE READ_BOUND
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE APPLY_BOUNDARY(N,XPER,YPER,ZPER,IPERIODICITY,XMPIELRANK)
!> @brief
!> This subroutine assigns the correct boundary condition code for all the bounded elements
IMPLICIT NONE
REAL,INTENT(IN)::XPER,YPER,ZPER
INTEGER,INTENT(IN)::IPERIODICITY,N
REAL::SMALL,tolerance,dist,temp_x
INTEGER::I,K,j,kk,ii,kmaxe,jj1,jj2,ji,l,IBLEED
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK
integer::dum1,dum2,N_NODE
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST

	KMAXE=XMPIELRANK(N)
 

TOLERANCE=TOLSMALL !> The tolerance can have a significant impact on the periodic boundary conditions matching rules





jj1=0
if (dimensiona.eq.3)then

jj2=0
do i=1,n_boundaries
if ((IBOUND(N,I)%ICODE.EQ.5).or.(IBOUND(N,I)%ICODE.EQ.50))then
jj2=jj2+1
end if


end do


!$OMP DO
DO I=1,KMAXE			! For ALL ELEMENTS
    if (ielem(n,i)%interior.eq.1)then		! THAT HAVE AT LEAST ONE UNKNWON NEIGHBOUR
	    IF (ielem(n,i)%nofbc.GT.0)THEN		! THAT HAVE AT LEAST ESTABLISHED A BOUNDARY CONDITION CODE
			
		  DO J=1,ielem(n,i)%IFCA			! LOOP ALL THEIR FACES
		      if (IELEM(N,I)%IBOUNDS(J).gt.0)then
		     IF ((IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ICODE.EQ.5).or.(IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ICODE.EQ.50))THEN	!IF ANY OF THEM HAS A PERIODIC BOUNDARY CONDITION THEN
				    if (IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ISHAPE.EQ.5)then
				    N_NODE=4
				    else
				    N_NODE=3
				    end if
				    DO Kk=1,n_node
				      NODES_LIST(kk,1:3)=inoder(IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ibl(kk))%CORD(1:3)
				    END DO
				    vext(1,1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
				    
			   do ii=1,n_boundaries				! loop all the boundaries
				if ((ii.ne.IELEM(N,I)%IBOUNDS(J)).and.((IBOUND(n,ii)%icode.eq.5).or.(IBOUND(n,ii)%icode.eq.50)).and.&
(IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ishape.eq.IBOUND(n,ii)%ishape))then
				      if ((IBOUND(n,ii)%localn(1).gt.0)) then	! excluding itself, and of same shape type
 				   if (ielem(n,IBOUND(N,ii)%localn(1))%ihexgl.ne.ielem(n,i)%ihexgl)then
! 				    
				    DO Kk=1,n_node
				      NODES_LIST(kk,1:3)=inoder(ibound(N,ii)%ibl(kk))%CORD(1:3)
				    END DO
				    vext(2,1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
				    
				    dist=distance3(n,VEXT)
				      
                    IF(PER_ROT.EQ.0)THEN
				      if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
				      ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
				
				      
! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL).or.((abs(dist-zper)).lt.TOLSMALL))then

				      IBOUND(N,ii)%localn(2)=i;IBOUND(N,ii)%cpun(2)=n
				      IELEM(N,I)%INEIGHG(J)=ielem(n,IBOUND(N,ii)%localn(1))%ihexgl

			


! 				      IELEM(N,I)%INEIGHG(J)=IELEM(N,IBOUND(N,IELEM(N,I)%IBOUNDS(J))%localn(1))%IHEXGL
! 				      
				      jj1=jj1+1
 				      go to 101
 				      end if
				      end if
				       if (((abs(vext(2,2)-Yper).lt.tolsmall).or.(abs((abs(vext(2,2)-Yper))-Yper).lt.tolsmall)).and.&
				      ((abs(vext(1,2)-Yper).lt.tolsmall).or.(abs((abs(vext(1,2)-Yper))-Yper).lt.tolsmall)))then
				      if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
				      
! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL).or.((abs(dist-zper)).lt.TOLSMALL))then

				      IBOUND(N,ii)%localn(2)=i;IBOUND(N,ii)%cpun(2)=n
				      IELEM(N,I)%INEIGHG(J)=ielem(n,IBOUND(N,ii)%localn(1))%ihexgl
				      jj1=jj1+1
 				      go to 101
 				      end if
				      end if
				       if (((abs(vext(2,3)-Zper).lt.tolsmall).or.(abs((abs(vext(2,3)-Zper))-Zper).lt.tolsmall)).and.&
				      ((abs(vext(1,3)-Zper).lt.tolsmall).or.(abs((abs(vext(1,3)-Zper))-Zper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,1)-vext(1,1)).lt.tolsmall))then
				      
! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL).or.((abs(dist-zper)).lt.TOLSMALL))then

				      IBOUND(N,ii)%localn(2)=i;IBOUND(N,ii)%cpun(2)=n
				      IELEM(N,I)%INEIGHG(J)=ielem(n,IBOUND(N,ii)%localn(1))%ihexgl
				      jj1=jj1+1
 				      go to 101
 				      end if
				      end if
				    ELSE
                        vext(2,:)=ROTATE_PER(vext(2,:),IBOUND(n,ii)%icode,angle_per)
                        if ((abs(vext(1,1)-vext(2,1)).lt.tol_per).and.&
                            (abs(vext(1,2)-vext(2,2)).lt.tol_per).and.&
                            (abs(vext(1,3)-vext(2,3)).lt.tol_per)) then              
                                IBOUND(N,ii)%localn(2)=i
                                IBOUND(N,ii)%cpun(2)=n
                                IELEM(N,I)%INEIGHG(J)=ielem(n,IBOUND(N,ii)%localn(1))%ihexgl
                                jj1=jj1+1
                            go to 101
 				      end if
				    END IF
				    end if

				      
				end if
 			    end if
			  end do
			  101 continue
		      end if
		    end if
		    end do
	      end if
    end if
end do
!$OMP END DO


	 	
		    
	else
jj2=0
do i=1,n_boundaries
if (IBOUND(n,i)%icode.eq.5)then
jj2=jj2+1
end if
end do
jj1=0
!$OMP DO
DO I=1,KMAXE			!> ALL ELEMENTS
    if (ielem(n,i)%interior.eq.1)then		! THAT HAVE AT LEAST ONE UNKNWON NEIGHBOUR
	    IF (ielem(n,i)%nofbc.GT.0)THEN		! THAT HAVE AT LEAST ESTABLISHED A BOUNDARY CONDITION CODE
		  DO J=1,ielem(n,i)%IFCA			! LOOP ALL THEIR BOUNDARY FACES
		      if (IELEM(N,I)%IBOUNDS(J).gt.0)then


			  IF (BLEED.EQ.1)THEN
		      IF (IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ICODE.EQ.99)THEN


					!ASSIGN THE BLEED ZONE
					N_NODE=2
				    DO Kk=1,n_node
				      NODES_LIST(kk,1:2)=inoder(IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ibl(kk))%CORD(1:2)
				    END DO
				    vext(1,:)=CORDINATES2(N,NODES_LIST,N_NODE)

					DO IBLEED=1,BLEED_NUMBER
						IF (((vext(1,1).GE.bleed_start(IBLEED,1)).AND.(vext(1,1).LE.bleed_END(IBLEED,1))).AND.((vext(1,2).GE.bleed_start(IBLEED,2)).AND.(vext(1,2).LE.bleed_END(IBLEED,2))))THEN

						IELEM(N,I)%BLEEDN(J)=IBLEED	!ASSIGN THE BLEED NUMBER
						END IF
					END DO

				END IF
				END IF


		     IF (IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ICODE.EQ.5)THEN	! IF ANY OF THEM HAS A PERIODIC BOUNDARY CONDITION THEN
		     
				    N_NODE=2
				    DO Kk=1,n_node
				      NODES_LIST(kk,1:2)=inoder(IBOUND(N,IELEM(N,I)%IBOUNDS(J))%ibl(kk))%CORD(1:2)
				    END DO
				    vext(1,:)=CORDINATES2(N,NODES_LIST,N_NODE)
			   do ii=1,n_boundaries				! loop all the boundaries
				if (((ii.ne.IELEM(N,I)%IBOUNDS(J)).and.(IBOUND(n,ii)%icode.eq.5)))then
				      if(IBOUND(n,ii)%localn(1).gt.0)then	! excluding itself, and of same shape type
				    if (ielem(n,IBOUND(N,ii)%localn(1))%ihexgl.ne.ielem(n,i)%ihexgl)then
				    DO Kk=1,n_node
				      NODES_LIST(kk,1:2)=inoder(ibound(N,ii)%ibl(kk))%CORD(1:2)
				    END DO
				    vext(2,:)=CORDINATES2(N,NODES_LIST,N_NODE)
				      dist=distance2(n,VEXT)
				      


				      if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
				      ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall))then


! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL))then
						
				      IBOUND(N,ii)%localn(2)=i;IBOUND(N,ii)%cpun(2)=n
				      IELEM(N,I)%INEIGHG(J)=ielem(n,IBOUND(N,ii)%localn(1))%ihexgl
				      jj1=jj1+1
					go to 201
				      end if
				      end if
				      if (((abs(vext(2,2)-yper).lt.tolsmall).or.(abs((abs(vext(2,2)-yper))-yper).lt.tolsmall)).and.&
				      ((abs(vext(1,2)-yper).lt.tolsmall).or.(abs((abs(vext(1,2)-yper))-yper).lt.tolsmall)))then
				      if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall))then
				      !

! 				      if (((abs(dist-xper)).lt.TOLSMALL).or.((abs(dist-yper)).lt.TOLSMALL))then
						
				      IBOUND(N,ii)%localn(2)=i;IBOUND(N,ii)%cpun(2)=n
				      IELEM(N,I)%INEIGHG(J)=ielem(n,IBOUND(N,ii)%localn(1))%ihexgl
				      jj1=jj1+1
					go to 201
				      end if
				      end if
				end if
				end if
				end if
				
			  end do
			  201 continue
			  end if
		      end if
		    end do
	      end if
    end if
end do			      
!$OMP END DO 

  
	
end if    
	      


      



	
    

END SUBROUTINE APPLY_BOUNDARY





















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
END MODULE BOUNDARY
