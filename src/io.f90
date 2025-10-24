MODULE IO
USE MPIINFO
USE DECLARATION
USE FLOW_OPERATIONS
use ISO_C_BINDING
USE TRANSFORM
IMPLICIT NONE
contains


SUBROUTINE OUTWRITEGRIDB
 !> @brief
!> This subroutine writes the grid file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
INTEGER,DIMENSION(70)::IVALID
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,kmmg
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
     
if (n.eq.0)then

		INQUIRE (FILE='GRID.plt',EXIST=HEREV)

	if (herev)then



	else

NullPtr = 0
      Debug   = 0
      FileType = 1
      VIsDouble = 1
      IMax    = imaxN
      JMax    = IMAXE
      KMax    = 0
      ZoneType = 5
      SolTime = 360.0
      StrandID = 0
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
NULCHAR = CHAR(0)


ierr =  TecIni112('SIMPLE DATASET'//NULCHAR, &
                    'X Y Z'//NULCHAR, &
                    'GRID.plt'//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)


 ierr= TecZne112('GRID1'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    ShrConn)




	
    allocate(xbin(imaxn))
    allocate(ybin(imaxn))
    allocate(zbin(imaxn))
	if (binio.eq.0)then
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	xbin(i)=x/scaler
	ybin(i)=y/scaler
 	Zbin(i)=Z/scaler
	END DO
	CLOSE(96)
	else
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96)j,x,y,z
	xbin(i)=x/scaler
	ybin(i)=y/scaler
 	Zbin(i)=Z/scaler
	END DO
	CLOSE(96)
	end if


    ierr = TECDAT112(imaxn,xbin,1) 
   
    ierr = TECDAT112(imaxn,ybin,1)

     ierr = TECDAT112(imaxn,zbin,1)
   
    deallocate(xbin,YBIN,zbin)
    
    
    IF (BINIO.EQ.0)THEN
    OPEN(98,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(8,1))
    icon=0
     cv=0
		DO K=1,iMAXE
               
 		read(98,*)i,Icon(1,1),icon(2,1),icon(3,1),icon(4,1),Icon(5,1),icon(6,1),icon(7,1),icon(8,1)
    
		ierr = TECNODE112(8,icon)
    !cv=cv+4
        	
		END DO
 		close(98)
		!ierr = TECNOD112(icon)
		deallocate(icon)	
    ELSE
     OPEN(98,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(8,1))
    icon=0
     cv=0
		DO K=1,iMAXE
               
 		read(98)i,Icon(1:8,1)
    
		ierr = TECNODE112(8,icon)
    !cv=cv+4
        	
		END DO
 		close(98)
		!ierr = TECNOD112(icon)
		deallocate(icon)	
    
    
    
    END IF
    
    
    
  ierr = TECEND112()


END IF


end if


	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	
	

	
	

END SUBROUTINE OUTWRITEGRIDB


SUBROUTINE OUTWRITEGRIDB2D
 !> @brief
!> This subroutine writes the grid file in tecplot binary format in 2D
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
INTEGER,DIMENSION(70)::IVALID
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
     
if (n.eq.0)then
INQUIRE (FILE='GRID.plt',EXIST=HEREV)

if (herev)then



else


NullPtr = 0
      Debug   = 0
      FileType = 1
      VIsDouble = 1
      IMax    = imaxN
      JMax    = IMAXE
      KMax    = 0
      ZoneType = 3
      SolTime = 360.0
      StrandID = 0
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
NULCHAR = CHAR(0)


ierr =  TecIni112('SIMPLE DATASET'//NULCHAR, &
                    'X Y'//NULCHAR, &
                    'GRID.plt'//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)


 ierr= TecZne112('GRID1'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    ShrConn)




	
    allocate(xbin(imaxn))
    allocate(ybin(imaxn))
	IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96,*)j,x,y
	xbin(i)=x/scaler
	ybin(i)=y/scaler
	END DO
	CLOSE(96)
	ELSE
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96)j,x,y
	xbin(i)=x/scaler
	ybin(i)=y/scaler
	END DO
	CLOSE(96)
	END IF


    ierr = TECDAT112(imaxn,xbin,1)  !!! why not xbin instead of xbin(1) ??
   
    ierr = TECDAT112(imaxn,ybin,1)

    
   
    deallocate(xbin,YBIN)
    
    IF (BINIO.EQ.0)THEN
    
    OPEN(98,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(4,1))
    icon=0
     cv=0
		DO K=1,iMAXE
               
 		read(98,*)i,Icon(1,1),icon(2,1),icon(3,1),icon(4,1)
    
		ierr = TECNODE112(4,icon)
    !cv=cv+4
        	
		END DO
 		close(98)
		!ierr = TECNOD112(icon)
		deallocate(icon)
		
		
    ELSE
	   OPEN(98,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(4,1))
    icon=0
     cv=0
		DO K=1,iMAXE
               
 		read(98)i,Icon(1:4,1)
    
		ierr = TECNODE112(4,icon)
    !cv=cv+4
        	
		END DO
 		close(98)
		!ierr = TECNOD112(icon)
		deallocate(icon)
    
    
    
    END IF
         
        
  ierr = TECEND112()


END IF

end if	



	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	

	

	
	

END SUBROUTINE OUTWRITEGRIDB2D


SUBROUTINE OUTWRITE3N
 !> @brief
!> This subroutine is solely for debugging
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)


KMAXE=XMPIELRANK(N)



! DUMG=KMAXE
! call mpi_barrier(mpi_comm_world,IERROR)
! 
! CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
! IMAXP=DUML

! ALLOCATE(ICELL(IMAXE))
! ICELL=0

! DO I=1,KMAXE
!   ICELL(I)=IELEM(N,I)%IHEXGL
! END DO


! IF (N.EQ.0)THEN
! 	ALLOCATE(ICELLA(IMAXP*ISIZE))
! 	 ICELLA=0
! 
! END IF
! 
! call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)
! 
! ! if (n.eq.0)then

! ! 
! ! 
! ! end if
! 
! call mpi_barrier(mpi_comm_world,IERROR)
! deallocate (icell)




if (n.eq.0)then


allocate(xbin(imaxe))
allocate(variables(8))

 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	

end if



IF (ITESTCASE.EQ.3)THEN
	nvar1=2
if (n.eq.0)ierr =  TecIni112('sols1'//NULCHAR, &
                    'solution,sols2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
END IF


	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 5
!       if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = 0.0
!       else
!       SolTime = IT
! 
!       end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


!  ALLOCATE(VALUESA(IMAXP*ISIZE))
!   allocate(xbin(imaxe))
! 	VALUESA=0.0
! 
!  END IF
! 
!   call mpi_barrier(MPI_COMM_WORLD,IERROR)
!   ALLOCATE(VALUESS(imaxp))
!   VALUESS=0.0
    
   


!   IF (ITESTCASE.LT.3)THEN
!     DO I=1,KMAXE
!       VALUESS(I)=U_C(I)%VAL(1,1)
!     END DO
! 
!     call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
! 
!     IF (N.EQ.0)THEN
!     do i=1,imaxp*isize
! 	if (icella(i).gt.0)then
! 	xbin(icella(i))=valuesa(i)
! 	end if
!     end do
    xbin(1:imaxe)=xmpie(1:imaxe)
	  

    ierr = TECDAT112(imaxe,XBIN,1)


    ierr = TECDAT112(imaxe,XBIN,1)

  ierr = TECEND112()
    deallocate(xbin,variables, valuelocation)
    END IF

     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 
!   END IF
! 
! 
! 
! 
!     
! 
! 
! 
! 
! 
!   
! 
! 
! 
!   IF (ITESTCASE.ge.3)THEN
! 			  do j=1,nvar1
! 				  
! 				    if (j.eq.1)then
! 				    
! 				    DO I=1,KMAXE
! 				      
! 				      VALUESS(I)=U_C(I)%VAL(1,1)
! 				    END DO
! 				    else
! 				    if (j.lt.6)then
! 				    DO I=1,KMAXE
! 				      VALUESS(I)=U_C(I)%VAL(1,j)/U_C(I)%VAL(1,1)
! 				    END DO
! 				    end if
! 				    end if
! 				    if (j.eq.6)then
! 
! 				      DO I=1,KMAXE
! 											VARIABLES(1)=U_C(I)%VAL(1,1)
! 											VARIABLES(2)=(U_C(I)%VAL(1,2)/VARIABLES(1))
! 											VARIABLES(3)=(U_C(I)%VAL(1,3)/VARIABLES(1))
! 											VARIABLES(4)=(U_C(I)%VAL(1,4)/VARIABLES(1))
! 											VARIABLES(5)=U_C(I)%VAL(1,5)
! 									    VARIABLES(6)=((GAMMA-1.0))*((VARIABLES(5))-0.5*VARIABLES(1)*(((VARIABLES(2)**2.0)+(VARIABLES(3)**2.0)+(VARIABLES(4)**2.0))))
! 									    VALUESS(I)=VARIABLES(6)
! 
! 				    END DO
! 
! 				    
! 				    end if 
! 
! 				    IF (TURBULENCE.NE.1)THEN
! 					      IF (PASSIVESCALAR.GT.0)THEN
! 						      IF (J.EQ.7)THEN
! 						      DO I=1,KMAXE
! 							VALUESS(I)=U_CT(I)%VAL(1,1)
! 						      END DO
! 						      end if
! 					      END IF
! 				    END IF
!     
! 
! 					  if (turbulence.eq.1)then
! 								if (j.eq.7)then
! 										do i=1,kmaxe
! 														leftv(1:5)=U_C(I)%VAL(1,1:5)
! 														rightv(1:5)=leftv(1:5)
! 										  Call SUTHERLANDIi(N,LEFTV,RIGHTV,VISCL,LAML,PRES,RRES,GAMMA,VISC,BETAAS,SUTHER,PRANDTL)
! 												  Variables(1)=Viscl(1)/VISC
! 												  VALUESS(I)=VARIABLES(1)
! 									      end do
! 								end if
! 								if (j.eq.8)then
! 								  do i=1,kmaxe
! 				      leftv(1:5)=U_C(I)%VAL(1,1:5)
! 				      rightv(1:5)=leftv(1:5)
! 							Call SUTHERLANDII(N,LEFTV,RIGHTV,VISCL,LAML,PRES,RRES,GAMMA,VISC,BETAAS,SUTHER,PRANDTL)
! 							 EDDYFL(1)=IELEM(N,i)%WALLDIST
! 							  IF (TURBULENCEMODEL.EQ.1)THEN
! 							 EDDYFL(2)=U_CT(i)%VAL(1,1)*U_C(i)%VAL(1,1)
! 							 EDDYFL(3)=0
! 							  
! 							  END IF
! 							  IF (TURBULENCEMODEL.EQ.2)THEN
! 							 EDDYFL(2)=U_CT(i)%VAL(1,1)
! 							 EDDYFL(3)=U_CT(i)%VAL(1,2)
! 							  
! 							  END IF
! 							EDDYFL(4:6)=ILOCAL_RECON3(i)%GRADS(1,1:3)
! 							EDDYFL(7:9)=ILOCAL_RECON3(i)%GRADS(2,1:3)
! 							EDDYFL(10:12)=ILOCAL_RECON3(i)%GRADS(3,1:3)
! 							IF (TURBULENCEMODEL.EQ.2)THEN
! 							EDDYFL(13:15)=ILOCAL_RECON3(i)%GRADS(4,1:3)
! 							EDDYFL(16:18)=ILOCAL_RECON3(i)%GRADS(5,1:3)
! 							END IF
! 							EDDYFR=EDDYFL
! 							
! 							call EDDYVISCOo(N,VISCL,LAML,TURBMV,ETVM,LEFTV,RIGHTV,EDDYFL,EDDYFR)
! 							Variables(8) =  (VISCL(3))/VISC	
! 			VALUESS(I)=VARIABLES(8)
! 
! 								  end do
! 								end if
!    
! 							if (turbulencemodel.eq.2)then
! 							    if ((j.eq.9))then
! 							      DO I=1,KMAXE
! 								VALUESS(I)=U_CT(I)%VAL(1,1)
! 							      END DO
! 
! 							    end if
! 							    if ((j.eq.10))then
! 							      DO I=1,KMAXE
! 								VALUESS(I)=U_CT(I)%VAL(1,2)
! 							      END DO
! 
! 							    end if
! 							end if
! 
! 							IF (passivescalar.gt.0)THEN
! 								      IF (J.EQ.NVAR1-1)THEN
! 								    DO I=1,KMAXE
! 								      VALUESS(I)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
! 								    END DO
! 								      END IF
! 							 END IF
! 
! 
! 					  end if !turbulence
! 
! 
! 				    IF (IVORTEX.EQ.1)THEN
! 					      IF (J.EQ.NVAR1)THEN
! 						    DO I=1,KMAXE
! 						      VALUESS(I)=IELEM(N,i)%VORTEX
! 						    END DO
! 					    END IF
! 				    END IF
! 
! 
! 
! 
!       
!       CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 
!     
!     call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
! 
!      
! !       print*,"writing output var",j,nvar1,n
! 
!     
!     IF (N.EQ.0)THEN
!     do i=1,imaxp*isize
! 	if (icella(i).gt.0)then
! 	xbin(icella(i))=valuesa(i)
! 	end if
!     end do
!     
!     ierr = TECDAT112(imaxe,xbin,1)
!     END IF
!     
!       
! !       print*,"writing 1r",n
!       
!   end do
! 
! 
! 
!   END IF
!   IF (N.EQ.0)THEN
!   ierr = TECEND112()
!   DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
!   deallocate(out1)
!   END IF
!   DEALLOCATE (VALUESS)
!   
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
! 
!   CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 
! 
! 
! 
! 
! 
! ! 
! !     allocate(xbin(imaxe))
! !     ALLOCATE(FBIN(1,IMAXE))
! ! 
! !  
! !     do i=1,1
! !     fbin(i,:)=10.0
! !     end do
! !     
! !     do i=1,1
! !     xbin(:)=fbin(i,:)
! ! 
! !     ierr = TECDAT112(imaxe,xbin,1)  !!! why not xbin instead of xbin(1) ??
! !    end do
! 
! 		
!          
!         
!  deallocate(variables)









	
	
	
	
	

	
	

END SUBROUTINE OUTWRITE3N



SUBROUTINE MOVIE
 !> @brief
!> This subroutine writes only the 3D solution without the grid in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::FBIN
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
 character(LEN=:),allocatable::out1
 character*1 NULCHAR

      Integer::   Debug,III,NPts,NElm


      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)









allocate(variables(2))

KMAXE=XMPIELRANK(N)






if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="MOV_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)


end if
call mpi_barrier(mpi_comm_world,IERROR)


 IF (ITESTCASE.EQ.4)THEN
 NVAR1=2



              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'v_MAG,Q'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)




END IF




if (n.eq.0)then
allocate (Valuelocation(nvar1))


      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 5

      SolTime = T

      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)



  allocate(xbin(imaxe),xbin2(imaxe))


 eLSE
 allocate(xbin2(1))
 end if



  allocate(valuess(kmaxe))

!   call mpi_barrier(MPI_COMM_WORLD,IERROR)









		    DO I=1,KMAXE


            leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)

		    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)

		    valuess(i)=sqrt(leftv(2)**2+leftv(3)**2+leftv(4)**2)

			end do

		    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
					IF (N.EQ.0)THEN
					do i=1,imaxe
				xbin(XMPI_RE(I))=xbin2(I)
				end do
					ierr = TECDAT112(imaxe,xbin,1)
					END IF





		DO I=1,KMAXE

		  VALUESS(i)=ielem(n,i)%vortex(1)
		END DO


		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF












  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1)
  END IF

  DEALLOCATE (VALUESS,VARIABLES,xbin2)











END SUBROUTINE MOVIE





SUBROUTINE OUTWRITE3vb
 !> @brief
!> This subroutine writes only the 3D solution without the grid in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
 character(LEN=:),allocatable::out1
 character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)


 






allocate(variables(14))

KMAXE=XMPIELRANK(N)






if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	

end if
call mpi_barrier(mpi_comm_world,IERROR)
 IF (ITESTCASE.LE.2)THEN
  NVAR1=4
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'solution1,SOLUTION2,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=8+PASSIVESCALAR
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     if (multispecies.eq.1)then
     NVAR1=10

     if (dg.eq.1)then
     NVAR1=10
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,species1,species2,vfraction,TROUBLED'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)



     else

      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,species1,species2,vfraction,AUX'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)


         end if

     
     else
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     end if
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=9+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF

	

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))


      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 5
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0

 



 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(imaxe),xbin2(imaxe))
	

 eLSE
 allocate(xbin2(1))
 end if



  allocate(valuess(kmaxe))

!   call mpi_barrier(MPI_COMM_WORLD,IERROR)
 
    
   


    IF (ITESTCASE.LE.2)THEN
		DO I=1,KMAXE
		 
     VALUESS(i)=U_C(I)%VAL(1,1)!0.0
    
		END DO
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF

    
		DO I=1,KMAXE
		
      VALUESS(i)=ielem(n,i)%inumneighbours
     
		END DO
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    
		DO I=1,KMAXE
		  VALUESS(i)=ielem(n,i)%TROUBLED
		END DO
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF

    
		DO I=1,KMAXE
		  VALUESS(i)=IELEM(N,I)%ADMIS
		END DO
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 IF (N.EQ.0)THEN
		 do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
		    DO I=1,KMAXE
		    
            
            leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
            
		    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
		      VALUESS(i)=leftv(kkd)
			if (kkd.eq.5)then
			
            VALUESS(i)=U_C(I)%VAL(1,kkd)
           
			end if
		    END DO
		    
		    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF
		end do
		
		
		
    
		DO I=1,KMAXE
		      
            leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
            
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		END DO
		
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF

		
		
		
		
		
                if (multispecies.eq.1)then
                
                DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,6)
                END DO
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
                    IF (N.EQ.0)THEN
                    do i=1,imaxe
                xbin(XMPI_RE(I))=xbin2(I)
                end do
                    ierr = TECDAT112(imaxe,xbin,1)
                    END IF
			
			
                    DO I=1,KMAXE
                        VALUESS(i)=U_C(I)%VAL(1,7)
                        END DO
                    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
                    IF (N.EQ.0)THEN
                    do i=1,imaxe
                xbin(XMPI_RE(I))=xbin2(I)
                end do
                    ierr = TECDAT112(imaxe,xbin,1)
                    END IF
                    
                    
                       DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,8)
                END DO
                
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF
                    

                    IF (DG.EQ.1)THEN
                     DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%TROUBLED
                END DO

                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF
                Else
                        DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%REDUCE
                END DO

                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
                        IF (N.EQ.0)THEN
                        do i=1,imaxe
                xbin(XMPI_RE(I))=xbin2(I)
                end do
                        ierr = TECDAT112(imaxe,xbin,1)
                        END IF







                end if





                    
                    
			else
			
			
                
                IF (MOOD.EQ.1)THEN
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%MOOD_O
                END DO
                ELSE
                DO I=1,KMAXE
                IF (ADDA.EQ.1)THEN

                VALUESS(i)=IELEM(N,I)%lwcx2!DISS!IELEM(N,I)%STENCIL_DIST

                ELSE
                VALUESS(i)=IELEM(N,I)%ggs!WCX(1)!TROUBLED!FILTERED

                END IF
                END DO
                END IF
                
                
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF
            
                
               
                
                
                DO I=1,KMAXE
                
                if (adda.eq.1)then
                
                valuess(i)=ielem(n,i)%diss
                else
                
                
                VALUESS(i)=ielem(n,i)%FULL!WCX(1)!IELEM(N,I)%ADMIS
                end if
                END DO
                
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF

              end if  
		
    
		  if (passivescalar.gt.0)then
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
		  END DO
		  
		  
		  call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  DO I=1,KMAXE
		      VALUESS(i)=ielem(n,i)%VOrtex(1)!%inumneighbours
		  END DO
		  
		  
		 call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
			DO I=1,KMAXE
			    VALUESS(i)=U_CT(I)%VAL(1,kkd)
			END DO
		      
			
		      call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			      IF (N.EQ.0)THEN
			      do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			      ierr = TECDAT112(imaxe,xbin,1)
			      END IF


		 end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1)
  END IF
  
  DEALLOCATE (VALUESS,VARIABLES,xbin2)
  

  





	
	

END SUBROUTINE OUTWRITE3vb


SUBROUTINE OUTWRITEtec3dbp
 !> @brief
!> This subroutine writes only the 3D solution without the grid in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2,xbin3
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
 character(LEN=:),allocatable::out1
 character*1 NULCHAR
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)


 






allocate(variables(15),icon(8,1))

KMAXE=XMPIELRANK(N)





 NullPtr = 0
      Debug   = 0
      FileType = 0
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') N
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	
	
	
	
	
 IF (ITESTCASE.LE.2)THEN
  NVAR1=7
  ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,solution1,SOLUTION2,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=11+PASSIVESCALAR
  if (passivescalar.gt.0)then
 ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     if (multispecies.eq.1)then
     NVAR1=12
      ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,species1,species2,vfraction'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     else
    ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     end if
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=12+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF

	

	

allocate (Valuelocation(nvar1))


      IMax    = kmaxn
      JMax    = kMAXe
      KMax    = 0
      ZoneType = 5
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0
Valuelocation(1:3)=1
 



 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(Kmaxn),xbin2(kmaxn),xbin3(kmaxn))
	

 
  allocate(valuess(kmaxe))


    do i=1,kmaxn
        xbin(i)=INODER4(i)%CORD(1);
        xbin2(i)=iNODER4(i)%CORD(2);
        xbin3(i)=INODER4(i)%CORD(3)
        
    end do
    
    ierr = TECDAT112(kmaxn,xbin,1) 
   
    ierr = TECDAT112(kmaxn,xbin2,1)

     ierr = TECDAT112(kmaxn,xbin3,1)

     
    
     

    IF (ITESTCASE.LE.2)THEN
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(1,1)!0.0
		END DO
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		

    
		DO I=1,KMAXE
		  VALUESS(i)=n
		END DO
		
		
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		
    
		DO I=1,KMAXE
		  VALUESS(i)=IELEM(N,I)%iNUMNEIGHBOURS!%STENCIL_DIST
		END DO
		
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		

    
		DO I=1,KMAXE
		  VALUESS(i)=IELEM(N,I)%ADMIS
		END DO
		
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
		    DO I=1,KMAXE
		    leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
		      VALUESS(i)=leftv(kkd)
			if (kkd.eq.5)then
			VALUESS(i)=U_C(I)%VAL(1,kkd)!/U_C(I)%VAL(1,1)
			end if
		    END DO
		    
		    
		ierr = TECDAT112(kmaxe,VALUESS,1)
			
		end do
		
		
		
    
		DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		END DO
		
		
		
			ierr = TECDAT112(kmaxe,VALUESS,1)
			

		
		
		
		
		
                if (multispecies.eq.1)then
                
                DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,6)
                END DO
                
                 ierr = TECDAT112(kmaxe,VALUESS,1)
                   
			
			
                    DO I=1,KMAXE
                        VALUESS(i)=U_C(I)%VAL(1,7)
                        END DO
                    
                ierr = TECDAT112(kmaxe,VALUESS,1)
                    
                    
                    
                       DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,8)
                END DO
                
               
			ierr = TECDAT112(kmaxe,VALUESS,1)
			
                    
                    
                    
			else
			
			
                
                IF (MOOD.EQ.1)THEN
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%MOOD_O
                END DO
                ELSE
                DO I=1,KMAXE
                VALUESS(i)=ielem(n,i)%condition!IELEM(N,I)%STENCIL_DIST
                END DO
                END IF
                
                
                
			ierr = TECDAT112(kmaxe,VALUESS,1)
			
            
                
               
                
                
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%ADMIS
                END DO
                
                
			ierr = TECDAT112(kmaxe,VALUESS,1)
			

              end if  
		
    
		  if (passivescalar.gt.0)then
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
		  END DO
		  
		  
		  
			ierr = TECDAT112(kmaxe,VALUESS,1)
			
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  DO I=1,KMAXE
		      VALUESS(i)=ielem(n,i)%vortex(1)!%inumneighbours
		  END DO
		  
		  
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
			
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
			DO I=1,KMAXE
			    VALUESS(i)=U_CT(I)%VAL(1,kkd)
			END DO
		      
			
		      
		ierr = TECDAT112(kmaxe,VALUESS,1)
			      


		 end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    

    do i=1,kmaxe
    icon(1:8,1)=el_connect(I,1:8)
    ierr = TECNODE112(8,icon)
    end do
    
  
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2,xbin3,icon)
  
  
  DEALLOCATE (VALUESS,VARIABLES)
  

  





	
	

END SUBROUTINE OUTWRITEtec3dbp


SUBROUTINE OUTWRITEtec3dbpav
 !> @brief
!> This subroutine writes only the 3D solution without the grid in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2,xbin3
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
 character(LEN=:),allocatable::out1
 character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)


 






allocate(variables(15),icon(1,8))

KMAXE=XMPIELRANK(N)





 NullPtr = 0
      Debug   = 0
      FileType = 0
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') N
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	
 IF (ITESTCASE.LE.2)THEN
  NVAR1=7
  ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,solution1,SOLUTION2,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=11+PASSIVESCALAR
  if (passivescalar.gt.0)then
 ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     if (multispecies.eq.1)then
     NVAR1=12
      ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,species1,species2,vfraction'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     else
    ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     end if
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=12+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,passivescalar,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  TecIni112('sols'//NULCHAR, &
                    'X,Y,Z,Density,U,V,W,energy,Pressure,STEN1,STEN2,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF

	

	

allocate (Valuelocation(nvar1))


      IMax    = kmaxn
      JMax    = kMAXe
      KMax    = 0
      ZoneType = 5
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0
Valuelocation(1:3)=1
 



 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(Kmaxn),xbin2(kmaxn),xbin3(kmaxn))
	

 
  allocate(valuess(kmaxe))


    do i=1,kmaxn
        xbin=INODER4(i)%CORD(1);
        xbin2=iNODER4(i)%CORD(2);
        xbin3=INODER4(i)%CORD(3)
    end do
    
    ierr = TECDAT112(kmaxn,xbin,1) 
   
    ierr = TECDAT112(kmaxn,xbin2,1)

     ierr = TECDAT112(kmaxn,xbin3,1)


    IF (ITESTCASE.LE.2)THEN
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(1,1)!0.0
		END DO
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		

    
		DO I=1,KMAXE
		  VALUESS(i)=n
		END DO
		
		
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		
    
		DO I=1,KMAXE
		  VALUESS(i)=IELEM(N,I)%iNUMNEIGHBOURS!%STENCIL_DIST
		END DO
		
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		

    
		DO I=1,KMAXE
		  VALUESS(i)=IELEM(N,I)%ADMIS
		END DO
		
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
		
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
		    DO I=1,KMAXE
		    leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
		      VALUESS(i)=leftv(kkd)
			if (kkd.eq.5)then
			VALUESS(i)=U_C(I)%VAL(1,kkd)!/U_C(I)%VAL(1,1)
			end if
		    END DO
		    
		    
		ierr = TECDAT112(kmaxe,VALUESS,1)
			
		end do
		
		
		
    
		DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		END DO
		
		
		
			ierr = TECDAT112(kmaxe,VALUESS,1)
			

		
		
		
		
		
                if (multispecies.eq.1)then
                
                DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,6)
                END DO
                
                 ierr = TECDAT112(kmaxe,VALUESS,1)
                   
			
			
                    DO I=1,KMAXE
                        VALUESS(i)=U_C(I)%VAL(1,7)
                        END DO
                    
                ierr = TECDAT112(kmaxe,VALUESS,1)
                    
                    
                    
                       DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,8)
                END DO
                
               
			ierr = TECDAT112(kmaxe,VALUESS,1)
			
                    
                    
                    
			else
			
			
                
                IF (MOOD.EQ.1)THEN
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%MOOD_O
                END DO
                ELSE
                DO I=1,KMAXE
                VALUESS(i)=ielem(n,i)%condition!IELEM(N,I)%STENCIL_DIST
                END DO
                END IF
                
                
                
			ierr = TECDAT112(kmaxe,VALUESS,1)
			
            
                
               
                
                
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%ADMIS
                END DO
                
                
			ierr = TECDAT112(kmaxe,VALUESS,1)
			

              end if  
		
    
		  if (passivescalar.gt.0)then
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
		  END DO
		  
		  
		  
			ierr = TECDAT112(kmaxe,VALUESS,1)
			
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  DO I=1,KMAXE
		      VALUESS(i)=ielem(n,i)%vortex(1)!%inumneighbours
		  END DO
		  
		  
		
		ierr = TECDAT112(kmaxe,VALUESS,1)
			
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
			DO I=1,KMAXE
			    VALUESS(i)=U_CT(I)%VAL(1,kkd)
			END DO
		      
			
		      
		ierr = TECDAT112(kmaxe,VALUESS,1)
			      


		 end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    do i=1,kmaxe
    icon(1,1:8)=el_connect(I,1:8)
    ierr = TECNODE112(8,icon)
    end do
    
  
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2,xbin3,icon)
  
  
  DEALLOCATE (VALUESS,VARIABLES)
  

  





	
	

END SUBROUTINE OUTWRITEtec3dbpav


SUBROUTINE OUTWRITE3v
!> @brief
!> This subroutine writes only the 3D solution without the grid in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(10))

KMAXE=XMPIELRANK(N)

DUMG=KMAXE
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then
! write(1000+n,*)icella(:)
! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
		OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')

end if
call mpi_barrier(mpi_comm_world,IERROR)
 IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="solution"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
  
  if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if
  
      ELSE
     
     if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
   end if
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=7+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","K","OMEGA"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED, [9] = CELLCENTERED,[10] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","MU"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED, [9] = CELLCENTERED)'
	  end if
              
              end if
              if (turbulenceequations.eq.0)then
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	                if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","K","OMEGA"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                              if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","MU"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
             
              end if
              if (turbulenceequations.eq.0)then
                              if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
	  end if
              
              end if
	    
	    
	    
	    
	    end if
	   
 END IF

	IF (N.EQ.0)THEN
	    WRITE(97,*) ', SOLUTIONTIME=',T
	    END IF

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







Valuelocation(:)=0

 





 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(imaxe))
	VALUESA=ZERO

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=ZERO
    
   


    IF (ITESTCASE.LE.2)THEN
    DO I=1,KMAXE
      VALUESS(i)=U_C(I)%VAL(1,1)!0.0
    END DO
    
    call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    
			WRITE(97,*)XBIN(1:IMAXE)
			
    
    END IF

     
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(1,kkd)
		  if ((kkd.ge.2).and.(kkd.le.4))then
		  VALUESS(i)=U_C(I)%VAL(1,kkd)/U_C(I)%VAL(1,1)
		  end if
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
    
		DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    
		  if (passivescalar.gt.0)then
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:IMAXE)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  DO I=1,KMAXE
		      VALUESS(i)=ielem(n,i)%vortex(1)
		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:IMAXE)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,kkd)
		  END DO
		
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:IMAXE)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  CLOSE(97)
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







         
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITE3v


SUBROUTINE OUTWRITE3v2d
!> @brief
!> This subroutine writes only the 2D solution without the grid in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD, I_DOF
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(10))

KMAXE=XMPIELRANK(N)

DUMG=KMAXE
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
		OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')

end if
call mpi_barrier(mpi_comm_world,IERROR)
 IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="solution"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=5+PASSIVESCALAR
  if (passivescalar.gt.0)then
  
  if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
  
      ELSE
     
     if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED)'
   end if
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=7+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","K","OMEGA"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED, [9] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","MU"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
              
              end if
              if (turbulenceequations.eq.0)then
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
	  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	                if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","K","OMEGA"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                              if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","MU"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
	  end if
             
              end if
              if (turbulenceequations.eq.0)then
                              if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
	  end if
              
              end if
	    
	    
	    
	    
	    end if
	   
 END IF

	IF (N.EQ.0)THEN
	    WRITE(97,*) ', SOLUTIONTIME=',T
	    END IF

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







Valuelocation(:)=0

 





 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(imaxe))
	VALUESA=ZERO

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=ZERO
    
   


    IF (ITESTCASE.LE.2)THEN
        IF (DG == 1) THEN
!         DO I_DOF = 1, IELEM(N,ICONSIDERED)%IDEGFREE + 1
            do i=1,kmaxe
            VALUESS(i)=U_C(I)%VAL(1,1)!U_C(i)%VALDG(1,1,1)
            end do
!         END DO
            
        ELSE
         do i=1,kmaxe
            VALUESS(i)=U_C(i)%VAL(1,1)!0.0
            end do
        END IF
        
        call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

        IF (N.EQ.0)THEN
            do i=1,imaxp*isize
                if (icella(i).gt.0)then
                    xbin(icella(i))=valuesa(i)
                end if
            end do
        
            WRITE(97,*)XBIN(1:IMAXE)
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    ELSE IF (ITESTCASE.ge.3)THEN
		do kkd=1,4
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(1,kkd)
		  if ((kkd.ge.2).and.(kkd.le.3))then
		  VALUESS(i)=U_C(I)%VAL(1,kkd)/U_C(I)%VAL(1,1)
		  end if
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
    
		DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL cons2prim(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(4)
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    
		  if (passivescalar.gt.0)then
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:IMAXE)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  DO I=1,KMAXE
		      VALUESS(i)=ielem(n,i)%vortex(1)
		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:IMAXE)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,kkd)
		  END DO
		
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:IMAXE)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  CLOSE(97)
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







         
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITE3v2d

SUBROUTINE OUTWRITE3vb2d
!> @brief
!> This subroutine writes only the 2D solution without the grid in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(10))

KMAXE=XMPIELRANK(N)






if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	

end if

 IF (ITESTCASE.LE.2)THEN
  NVAR1=4
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'sol1,sol2,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=8+PASSIVESCALAR
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,STEN1,STEN2,SLOPE,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     if (multispecies.eq.1)then
     NVAR1=9
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,species1,species2,vf,aux'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     else
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,SLOPE'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     end if
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=9+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,passivescalar,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,passivescalar,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,passivescalar,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,vortex,k,omega'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,vortex,mu'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,vortex'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF

	
	
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 3
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)

allocate(xbin(1:imaxe),xbin2(1:imaxe))
 eLSE
 allocate(xbin2(1))
 end if
  


 
 
  allocate(valuess(1:kmaxe))
  
    
   


    IF (ITESTCASE.LE.2)THEN
    DO I=1,KMAXE
     VALUESS(i)=U_C(I)%VAL(1,1)!0.0
    END DO
    
    call MPI_GATHERv(valuess(1:kmaxe),kmaxe,MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    
    
    
	      
    
    
		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    

    
    
    DO I=1,KMAXE
        
      VALUESS(i)=ielem(n,i)%TROUBLED
      
    END DO
    
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    
    
    DO I=1,KMAXE
      VALUESS(i)=ielem(n,i)%wcx(1)!IELEM(N,I)%ADMIS
    END DO
    
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    
    if (initcond.eq.0)then
     DO I=1,KMAXE
      VALUESS(i)=ilocal_recon3(i)%cond(1)
    END DO
    else
    DO I=1,KMAXE
      VALUESS(i)=IELEM(N,I)%STENCIL_DIST
    END DO
    end if
    
   call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    
    
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,4
		DO I=1,KMAXE
        
        leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
        
		  CALL cons2prim(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(kkd)
          if (KKD.EQ.4)THEN
          
        
        VALUESS(i)=U_C(I)%VAL(1,KKD)
        
          
          END IF
		END DO
		
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		end do
    
		DO I=1,KMAXE
		 
        leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
        
		  CALL cons2prim(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(4)
		END DO
		
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    
                if (multispecies.eq.1)then
                DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,5)
                END DO
                else
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%full!ADMIS
                END DO
                end if
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
                if (multispecies.eq.1)then
                DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,6)
                END DO
                else
                IF (MOOD.EQ.1)THEN
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%MOOD_O
                END DO
                ELSE
                
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%TROUBLED
                END DO
                END IF
                end if
                
                
               call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF

                
                	
		
		
    
		  if (passivescalar.gt.0)then
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,turbulenceequations+passivescalar)
		  END DO
		  
		  
		 call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
		  
		  
		   if (itestcase.eq.3)then
                            if (multispecies.eq.1)then
                DO I=1,KMAXE
                VALUESS(i)=U_C(I)%VAL(1,7)
                END DO
                else
                DO I=1,KMAXE
                VALUESS(i)=ielem(n,i)%wcx(1)!vortex(1)
                END DO
                end if
                            
                            
                            
                            call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		
		
		if (multispecies.eq.1)then
		IF (MOOD.EQ.1)THEN
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%MOOD_O
                END DO
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
                IF (N.EQ.0)THEN
                do i=1,imaxe
                xbin(XMPI_RE(I))=xbin2(I)
                end do
                ierr = TECDAT112(imaxe,xbin,1)
                END IF
        else
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%REDUCE
                END DO
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
                IF (N.EQ.0)THEN
                do i=1,imaxe
                xbin(XMPI_RE(I))=xbin2(I)
                end do
                ierr = TECDAT112(imaxe,xbin,1)
                END IF
        
        
        
        end if
        end if
		  end if
    
		  if (itestcase.eq.4)then
                            DO I=1,KMAXE
                                VALUESS(i)=ielem(n,i)%ggs!vortex(1)
                            END DO
                            
                            
                            call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  DO I=1,KMAXE
		      VALUESS(i)=U_CT(I)%VAL(1,kkd)
		  END DO
		
		  
		  call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		  end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
!     read_ucns

     
    
    
    
  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1)
  END IF
  
  DEALLOCATE (VALUESS,variables,xbin2)


	
	

END SUBROUTINE OUTWRITE3vb2d


SUBROUTINE CHECKRES
!> @brief
!> This subroutine checks the presence of restart file
IMPLICIT NONE
LOGICAL::HERE
INTEGER::I,J,K,L,ITER,DIP
 CHARACTER(LEN=20)::PROC,RESTFILE

 	RESTFILE='RESTART.dat'
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	INQUIRE (FILE=RESTFILE,EXIST=HERE)
	IF (HERE) THEN
	OPEN(1083+N,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',access='stream')
	DIP=1
	if (IRES_UNSTEADY.eq.1)then
	READ(1083+N,pos=dip)ITER


	dip=dip+4
	READ(1083+N,pos=dip)RES_TIME
	dip=dip+8
	    if (initcond.eq.95)then
	    READ(1083+N,pos=dip)taylor
	    end if
	else
	  READ(1083+N,pos=dip)ITER
	end if
	RESTART=ITER
	CLOSE(1083+N)
	ELSE
	RESTART=0
        average_restart=0
	RES_TIME=0.0d0
	END IF


	IF (N.EQ.0)THEN
	PRINT*,"RESTARTING",ITER,RES_TIME,RESTART

	END IF
	

	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

END SUBROUTINE CHECKRES


SUBROUTINE OPEN_ARBITRARY(N,IMAXE,IMAXN,IMAXB)
!> @brief
!> This subroutine opens the grid files and establishes the number of cells, nodes, boundary conditions
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::IMAXE,IMAXN,IMAXB
	INTEGER::I,J,K,IOS,IOX,IOZ,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,IOY
	REAL::IX1,IX2,IX3
	CHARACTER(LEN=12)::PROC,VRTFILE,CELFILE,BNDFILE
	integer,ALLOCATABLE,dimension(:)::isent
	INTEGER,INTENT(IN)::N
		WRITE(PROC,FMT='(I10)') N
		CELFILE='GRID.cel'
		VRTFILE='GRID.vrt'
		BNDFILE='GRID.bnd'

	ALLOCATE(ISENT(3))
	
	IF (N.EQ.0)THEN
	
	if(binio.eq.0)then
	
		OPEN(8,FILE=CELFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
		OPEN(9,FILE=VRTFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOX)
		OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
	I=0
	J=0
	K=0

	
	DO 
		READ(8,*,IOSTAT=IOS)I
		IF (IOS.NE.0) THEN
			EXIT
		END IF
		
	END DO
		IMAXE=i
	DO 
		READ(9,*,IOSTAT=IOX)J
		IF (IOX.NE.0) THEN
			EXIT
		END IF
		
	END DO
		IMAXN=J
	DO 
		READ(10,*,IOSTAT=IOY)K
		IF (IOY.NE.0) THEN
			EXIT
		END IF
 		
	END DO
	 		IMAXB=K
	else
	
	
		OPEN(8,FILE=CELFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
		OPEN(9,FILE=VRTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOX)
		OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
	I=0
	J=0
	K=0

	IF (DIMENSIONA.EQ.3)THEN
	DO 
		READ(8,IOSTAT=IOS)I,I1,I2,I3,I4,I5,I6,I7,I8
		IF (IOS.NE.0) THEN
			EXIT
		END IF
		
	END DO
	
		IMAXE=i
		
	DO 
		READ(9,IOSTAT=IOX)J,IX1,IX2,IX3
		IF (IOX.NE.0) THEN
			EXIT
		END IF
		
	END DO
		IMAXN=J
	DO 
		READ(10,IOSTAT=IOY)K,I1,I2,I3,I4,I5
		IF (IOY.NE.0) THEN
			EXIT
		END IF
 		
	END DO
	 		IMAXB=K
	END IF
	
	IF (DIMENSIONA.EQ.2)THEN
	DO 
		READ(8,IOSTAT=IOS)I,I1,I2,I3,I4
		IF (IOS.NE.0) THEN
			EXIT
		END IF
		
	END DO
	
		IMAXE=i
		
	DO 
		READ(9,IOSTAT=IOX)J,IX1,IX2
		IF (IOX.NE.0) THEN
			EXIT
		END IF
		
	END DO
		IMAXN=J
	DO 
		READ(10,IOSTAT=IOY)K,I1,I2,I3,I4,I5
		IF (IOY.NE.0) THEN
			EXIT
		END IF
 		
	END DO
	 		IMAXB=K
	END IF
	
	
	end if

	CLOSE(8)
	CLOSE(9)
 	CLOSE(10)

	ISENT(1)=IMAXE
	ISENT(2)=IMAXN
	ISENT(3)=IMAXB

	

	END IF
	call MPI_Bcast( isent, 3, MPI_integer, 0 ,MPI_COMM_WORLD,IERROR )

	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	IMAXE=ISENT(1)
	IMAXN=ISENT(2)
	IMAXB=ISENT(3)
	DEALLOCATE(ISENT)
	
	
END SUBROUTINE OPEN_ARBITRARY


SUBROUTINE OPEN_INPUT1(N,ITT)
!> @brief
!> This subroutine opens the parameter file
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::ITT
	INTEGER,INTENT(IN)::N
	CHARACTER(LEN=12)::VRTFILE,CELFILE
! 		CELFILE='GRID.cel'
! 		VRTFILE='GRID.vrt'
! 		if (binio.eq.0)then
! 		OPEN(8,FILE=CELFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
! 		OPEN(9,FILE=VRTFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
! 		else
		OPEN(15,FILE='UCNS3D.DAT',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
END SUBROUTINE OPEN_INPUT1

SUBROUTINE CLOSE_INPUT1(N,ITT)
!> @brief
!> This subroutine closes the parameter file
 IMPLICIT NONE
 	INTEGER,INTENT(INOUT)::ITT
 	INTEGER,INTENT(IN)::N	
!  	CLOSE(8)
!  	CLOSE(9)
 	CLOSE(15)
END SUBROUTINE CLOSE_INPUT1

SUBROUTINE OPEN_INPUT(N,ITT)
!> @brief
!> This subroutine opens the grid files
	IMPLICIT NONE
	INTEGER,INTENT(INOUT)::ITT
	INTEGER,INTENT(IN)::N
	CHARACTER(LEN=12)::VRTFILE,CELFILE
		CELFILE='GRID.cel'
		VRTFILE='GRID.vrt'
		if (binio.eq.0)then
		OPEN(8,FILE=CELFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
		OPEN(9,FILE=VRTFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
		else
		OPEN(8,FILE=CELFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
		OPEN(9,FILE=VRTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
		END IF
	END SUBROUTINE OPEN_INPUT

SUBROUTINE CLOSE_INPUT(N,ITT)
!> @brief
!> This subroutine closes the grid files
 IMPLICIT NONE
 	INTEGER,INTENT(INOUT)::ITT
 	INTEGER,INTENT(IN)::N	
  	CLOSE(8)
  	CLOSE(9)

 END SUBROUTINE CLOSE_INPUT




SUBROUTINE READ_INPUT(N,XMPIELRANK,XMPINRANK,XMPIE,XMPIN,IELEM,INODE,IMAXN,IMAXE,IBOUND,IMAXB,XMPINNUMBER,SCALER,inoder)
!> @brief
!> This subroutine reads all the data from the grid files
	IMPLICIT NONE
	TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
	TYPE(NODE_NE),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::inoder
	TYPE(NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::XMPINNUMBER
	INTEGER,INTENT(IN)::N,IMAXN,IMAXE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIN
	TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK,XMPINRANK
	INTEGER,INTENT(INOUT)::IMAXB
	REAL,INTENT(INOUT)::SCALER
	integer,allocatable,dimension(:)::nodep,nodec
	INTEGER::I,J,JI,K,LM,IEX,KMAXN,KK,KMAXE,print_out,kk2,shap,nodal,iftrue,kxk2
	INTEGER,DIMENSION(8)::IDV
	INTEGER::IT1,IT2,IT3,IT4,IT5,IT6,IT7,IT8,ITX,INX,IT55,in,out
	REAL::X,Y,Z

	

	KK=0; I=0; J=0; LM=0 ;kk2=1
	XMIN(N)=tolbig; YMIN(N)=tolbig; ZMIN(N)=tolbig
	XMAX(N)=-tolbig; YMAX(N)=-tolbig; ZMAX(N)=-tolbig
	CALL OPEN_INPUT(N,ITT)
        !OPEN(82,FILE="GRID.epart",FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	KMAXE=XMPIELRANK(N)
	
	
 	if (dimensiona.eq.3)then
	
! 	ALLOCATE(XSIZE(0:ISIZE-1))
! 	XSIZE=-1
	ALLOCATE(iNODEr(imaxn))
	ALLOCATE(iNODEr2(imaxn))
	INODER2(1:imaxn)%NUMBEROFNEIB=0
	INODER(1:imaxn)%ITOR=0
	
	
	IF (BINIO.EQ.0)then
	
	DO J=1,IMAXE
	
	if (xmpie(j).eq.n)then

	KK=KK+1
	
	

	rEAD(8,*) ITX,IDV(1),IDV(2),IDV(3),IDV(4),IDV(5),IDV(6),IDV(7),IDV(8)

	    
		
				
	
	  do kxk2=1,8
	  inoder(idv(kxk2))%itor=idv(kxk2)
	 
	  end do
	  
	  
	  
	  
	  ielem(n,kk)%IHEXGL=itx
		IF ((IDV(5).NE.IDV(6)).AND.(IDV(6).NE.IDV(7)).AND.(IDV(7).NE.IDV(8)).AND.(IDV(3).NE.IDV(4)))THEN
		shap=1	;nodal=8
		IELEM(N,KK)%ISHAPE=shap
		ielem(n,kk)%ifca=6
		IELEM(N,KK)%nonodes=nodal
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%nodes(1:nodal)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		
		
		END IF
		IF ((IDV(3).EQ.IDV(4)).AND.(IDV(5).EQ.IDV(6)).AND.(IDV(6).EQ.IDV(7)).AND.(IDV(7).EQ.IDV(8)))THEN
		shap=2	;nodal=4!TETRAHEDRAL ELEMENT
		IELEM(N,KK)%ISHAPE=shap
		allocate(ielem(n,kk)%nodes(nodal))
		IELEM(N,KK)%nonodes=nodal
		ielem(n,kk)%ifca=4
		ielem(n,kk)%nodes(1:3)=idv(1:3)
		ielem(n,kk)%nodes(4)=idv(5)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		END IF
		IF ((IDV(3).NE.IDV(4)).AND.(Idv(5).EQ.Idv(6)).AND.(IDV(6).EQ.IDV(7)).AND.(IDV(7).EQ.IDV(8)))THEN
		shap=3	;nodal=5!PYRAMIDAL ELEMENT
		IELEM(N,KK)%ISHAPE=shap
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%ifca=5
		IELEM(N,KK)%nonodes=nodal
		ielem(n,kk)%nodes(1:nodal)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		END IF
		IF ((IDV(3).EQ.IDV(4)).AND.(IDV(5).NE.IDV(6)).AND.(IDV(7).EQ.IDV(8)))THEN
		shap=4	;nodal=6!PRISM ELEMENT
		IELEM(N,KK)%ISHAPE=shap
		ielem(n,kk)%ifca=5
		IELEM(N,KK)%nonodes=nodal
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%nodes(1:3)=idv(1:3)
		ielem(n,kk)%nodes(4:6)=idv(5:7)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		END IF
		  
		
			
	
	
	else
	rEAD(8,*)


	end if 
	end do
	ELSE
	DO J=1,IMAXE
	
	if (xmpie(j).eq.n)then

	KK=KK+1
	
	

	rEAD(8) ITX,IDV(1),IDV(2),IDV(3),IDV(4),IDV(5),IDV(6),IDV(7),IDV(8)

	    
		
				
	
	  do kxk2=1,8
	  inoder(idv(kxk2))%itor=idv(kxk2)
	 
	  end do
	  ielem(n,kk)%IHEXGL=itx
		IF ((IDV(5).NE.IDV(6)).AND.(IDV(6).NE.IDV(7)).AND.(IDV(7).NE.IDV(8)).AND.(IDV(3).NE.IDV(4)))THEN
		shap=1	;nodal=8
		IELEM(N,KK)%ISHAPE=shap
		ielem(n,kk)%ifca=6
		IELEM(N,KK)%nonodes=nodal
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%nodes(1:nodal)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		
		
		END IF
		IF ((IDV(3).EQ.IDV(4)).AND.(IDV(5).EQ.IDV(6)).AND.(IDV(6).EQ.IDV(7)).AND.(IDV(7).EQ.IDV(8)))THEN
		shap=2	;nodal=4!TETRAHEDRAL ELEMENT
		IELEM(N,KK)%ISHAPE=shap
		allocate(ielem(n,kk)%nodes(nodal))
		IELEM(N,KK)%nonodes=nodal
		ielem(n,kk)%ifca=4
		ielem(n,kk)%nodes(1:3)=idv(1:3)
		ielem(n,kk)%nodes(4)=idv(5)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		END IF
		IF ((IDV(3).NE.IDV(4)).AND.(Idv(5).EQ.Idv(6)).AND.(IDV(6).EQ.IDV(7)).AND.(IDV(7).EQ.IDV(8)))THEN
		shap=3	;nodal=5!PYRAMIDAL ELEMENT
		IELEM(N,KK)%ISHAPE=shap
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%ifca=5
		IELEM(N,KK)%nonodes=nodal
		ielem(n,kk)%nodes(1:nodal)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		END IF
		IF ((IDV(3).EQ.IDV(4)).AND.(IDV(5).NE.IDV(6)).AND.(IDV(7).EQ.IDV(8)))THEN
		shap=4	;nodal=6!PRISM ELEMENT
		IELEM(N,KK)%ISHAPE=shap
		ielem(n,kk)%ifca=5
		IELEM(N,KK)%nonodes=nodal
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%nodes(1:3)=idv(1:3)
		ielem(n,kk)%nodes(4:6)=idv(5:7)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		END IF
		  
		
			
	
	
	else
	rEAD(8)ITX,IDV(1),IDV(2),IDV(3),IDV(4),IDV(5),IDV(6),IDV(7),IDV(8)


	end if 
	end do
	
	
	
	
	
	
	
	
	
	
	END IF







	
	

	
	
	
	
 	else
 	
 	
 	
 	
	ALLOCATE(iNODEr(imaxn))
	ALLOCATE(iNODEr2(imaxn))
	INODER2(:)%NUMBEROFNEIB=0
	INODER(:)%ITOR=0
	
	if (binio.eq.0)then
	
	
	
	DO J=1,IMAXE
	
	if (xmpie(j).eq.n)then

	KK=KK+1
	
	

	rEAD(8,*) ITX,IDV(1),IDV(2),IDV(3),IDV(4)
	do kxk2=1,4
	  inoder(idv(kxk2))%itor=idv(kxk2)
	  
	  end do
	ielem(n,kk)%IHEXGL=itx
		IF ((IDV(3).NE.IDV(4)))THEN
		shap=5	;nodal=4	!QUADRILATERAL
		IELEM(N,KK)%ISHAPE=shap
		IELEM(N,KK)%nonodes=nodal
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%ifca=4
		ielem(n,kk)%nodes(1:nodal)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		Else
		
		shap=6	;nodal=3!TRIANGULAR
		IELEM(N,KK)%ISHAPE=shap
		allocate(ielem(n,kk)%nodes(nodal))
		IELEM(N,KK)%nonodes=nodal
		ielem(n,kk)%ifca=3
		ielem(n,kk)%nodes(1:NODAL)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		
		END IF
		


	 

	  else
	rEAD(8,*)



	end if 
	end do
	else
	
	DO J=1,IMAXE
	
	if (xmpie(j).eq.n)then

	KK=KK+1
	
	

	rEAD(8) ITX,IDV(1),IDV(2),IDV(3),IDV(4)
	do kxk2=1,4
	  inoder(idv(kxk2))%itor=idv(kxk2)
	  
	  end do
	ielem(n,kk)%IHEXGL=itx
		IF ((IDV(3).NE.IDV(4)))THEN
		shap=5	;nodal=4	!QUADRILATERAL
		IELEM(N,KK)%ISHAPE=shap
		IELEM(N,KK)%nonodes=nodal
		allocate(ielem(n,kk)%nodes(nodal))
		ielem(n,kk)%ifca=4
		ielem(n,kk)%nodes(1:nodal)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		Else
		
		shap=6	;nodal=3!TRIANGULAR
		IELEM(N,KK)%ISHAPE=shap
		allocate(ielem(n,kk)%nodes(nodal))
		IELEM(N,KK)%nonodes=nodal
		ielem(n,kk)%ifca=3
		ielem(n,kk)%nodes(1:NODAL)=idv(1:nodal)
		allocate(ielem(n,kk)%SURF(ielem(n,kk)%ifca))
		
		END IF
		


	 

	  else
	rEAD(8)ITX,IDV(1),IDV(2),IDV(3),IDV(4)



	end if 
	end do
	
	
	
	
	end if
	end if

	
	
	

	
	  
      


	      if (dimensiona.eq.3)then
	      
	      
	      if (binio.eq.0)then
	      DO j=1,IMAXN

	      IF (inoder(J)%ITOR.gt.0)THEN


		READ(9,*)INX,X,Y,Z
		x=x/scaler; y=y/scaler;  z=z/scaler
		XMIN(N)=MIN(XMIN(N),X)
		YMIN(N)=MIN(YMIN(N),Y)
		ZMIN(N)=MIN(ZMIN(N),Z)
		XMAX(N)=MAX(XMAX(N),X)
		YMAX(N)=MAX(YMAX(N),Y)
		ZMAX(N)=MAX(ZMAX(N),Z)
		
		
		  ALLOCATE(inoder(J)%CORD(1:3))
		  inoder(J)%CORD(1)=X
		  inoder(J)%CORD(2)=Y
		  inoder(J)%CORD(3)=Z
		  
		
		
		
		else
		read(9,*)
		end if
		
		
		!print_out = KMaxE/5


		END DO
		else
		 DO j=1,IMAXN

	      IF (inoder(J)%ITOR.gt.0)THEN


		READ(9)INX,X,Y,Z
		x=x/scaler; y=y/scaler;  z=z/scaler
		XMIN(N)=MIN(XMIN(N),X)
		YMIN(N)=MIN(YMIN(N),Y)
		ZMIN(N)=MIN(ZMIN(N),Z)
		XMAX(N)=MAX(XMAX(N),X)
		YMAX(N)=MAX(YMAX(N),Y)
		ZMAX(N)=MAX(ZMAX(N),Z)
		
		
		  ALLOCATE(inoder(J)%CORD(1:3))
		  inoder(J)%CORD(1)=X
		  inoder(J)%CORD(2)=Y
		  inoder(J)%CORD(3)=Z
		  
		
		
		
		else
		read(9)INX,X,Y,Z
		end if
		
		
		!print_out = KMaxE/5


		END DO
		
		
		
		
		
		
		
		end if
		else
		
		if (binio.eq.0)then
		 DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		  

		    
! -------------------FOR DEBUGGING ONLY -----------------------------------------!
! 		READ(9,'(I14,1X,3ES16.9)')INX,X,Y,Z
! -------------------FOR DEBUGGING ONLY -----------------------------------------!
		READ(9,*)INX,X,Y
		x=x/scaler; y=y/scaler
		XMIN(N)=MIN(XMIN(N),X)
		YMIN(N)=MIN(YMIN(N),Y)
		
		XMAX(N)=MAX(XMAX(N),X)
		YMAX(N)=MAX(YMAX(N),Y)
		
		
		
		  ALLOCATE(inoder(J)%CORD(1:2))
		  inoder(J)%CORD(1)=X
		  inoder(J)%CORD(2)=Y
		  Else

		  read(9,*)
		  end if
		
		
		
		
		
		
		!print_out = KMaxE/5


		END DO
		else
		 DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		  

		    
! -------------------FOR DEBUGGING ONLY -----------------------------------------!
! 		READ(9,'(I14,1X,3ES16.9)')INX,X,Y,Z
! -------------------FOR DEBUGGING ONLY -----------------------------------------!
		READ(9)INX,X,Y
		x=x/scaler; y=y/scaler
		XMIN(N)=MIN(XMIN(N),X)
		YMIN(N)=MIN(YMIN(N),Y)
		
		XMAX(N)=MAX(XMAX(N),X)
		YMAX(N)=MAX(YMAX(N),Y)
		
		
		
		  ALLOCATE(inoder(J)%CORD(1:2))
		  inoder(J)%CORD(1)=X
		  inoder(J)%CORD(2)=Y
		  Else

		  read(9)INX,X,Y
		  end if
		
		
		
		
		
		
		!print_out = KMaxE/5


		END DO
		
		
		
		
		


		end if
		end if
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	 CALL CLOSE_INPUT(N,ITT)






	x=xmax(n);CALL MPI_ALLREDUCE(x,XMax,1,MPI_DOUBLE_PRECISION,MPI_Max,MPI_COMM_WORLD,IERROR)
	x=Ymax(n);CALL MPI_ALLREDUCE(x,ymax,1,MPI_DOUBLE_PRECISION,MPI_Max,MPI_COMM_WORLD,IERROR)
	x=Zmax(n);CALL MPI_ALLREDUCE(x,zmax,1,MPI_DOUBLE_PRECISION,MPI_max,MPI_COMM_WORLD,IERROR)
	x=xmin(n);CALL MPI_ALLREDUCE(x,XMIN,1,MPI_DOUBLE_PRECISION,MPI_Min,MPI_COMM_WORLD,IERROR)
	x=ymin(n);CALL MPI_ALLREDUCE(x,yMIN,1,MPI_DOUBLE_PRECISION,MPI_Min,MPI_COMM_WORLD,IERROR)
	x=zmin(n);CALL MPI_ALLREDUCE(x,zMIN,1,MPI_DOUBLE_PRECISION,MPI_Min,MPI_COMM_WORLD,IERROR)
		
		
	  IF (DIMENSIONA.EQ.3)THEN
		
	  
	  
	  
	  CALL OPEN_INPUT(N,ITT)
	  DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		  ALLOCATE(INODER2(J)%XNE(1));INODER2(J)%XNE(1)=0 
		  end if
		  
	   END DO
	  
	  
	  DO J=1,IMAXE
	
	
	  
	  rEAD(8) ITX,IDV(1),IDV(2),IDV(3),IDV(4),IDV(5),IDV(6),IDV(7),IDV(8)

	  
	  
	    
		  if ((xmpie(j).NE.n))then
		      do kxk2=1,8
			  if (inoder(idv(kxk2))%itor.gt.0)then
! 			  XSIZE(XMPIE(J))=1
			  inoder2(idv(kxk2))%XNE(1)=INODER2(idv(kxk2))%XNE(1)+1
			  end if
		      end do
		  end if
	 END DO
	  
	  CALL CLOSE_INPUT(N,ITT)
	  
	  
	  
	  
	  
	  
	  
	  
	  CALL OPEN_INPUT(N,ITT)
	  DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		  IF (inoder2(J)%XNE(1).GT.0)THEN
		  ALLOCATE(INODER2(J)%XNEIB(inoder2(J)%XNE(1)));INODER2(J)%XNEIB=0 
		  inoder2(J)%XNE(1)=0
		  end if
		  END IF
	   END DO
	  
	  
	  DO J=1,IMAXE
	
	
	  
	  rEAD(8) ITX,IDV(1),IDV(2),IDV(3),IDV(4),IDV(5),IDV(6),IDV(7),IDV(8)

	  
	  
	    
		  if ((xmpie(j).NE.n))then
		      do kxk2=1,8
			  if (inoder(idv(kxk2))%itor.gt.0)then
! 			  XSIZE(XMPIE(J))=1
			  inoder2(idv(kxk2))%XNE(1)=INODER2(idv(kxk2))%XNE(1)+1
			  INODER2(IDV(KXK2))%XNEIB(inoder2(idv(kxk2))%XNE(1))=J
			  end if
		      end do
		  end if
	 END DO
	  
	  CALL CLOSE_INPUT(N,ITT)
	  
	  
	  
	  END IF
	  
	  
	  
	  
	  
	  

! DUMMYOUT=TIMEC5
! DUMMYIN=0.0
! CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
! !CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! TIMEC5=DUMMYIN





XPER = (XMAX(N))- XMIN(N)
YPER = (YMAX(N)) - YMIN(N)
ZPER = (ZMAX(N)) - ZMIN(N)
! 
! 
! XPER = abs(XMAX(N))! - XMIN(N)
! YPER = abs(YMAX(N))! - YMIN(N)
! ZPER = abs(ZMAX(N))! - ZMIN(N)



	IF (N.EQ.0) then
	OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	 WRITE(63,*)XPER,YPER,ZPER,"PERIODICS"
     	write(63,*)"MIN",XMIN(N),YMIN(N),ZMIN(N)
    	write(63,*)"MAX",XMAX(N),YMAX(N),ZMAX(N)
	CLOSE(63)
	END IF
	

	  
	END SUBROUTINE READ_INPUT





SUBROUTINE READ_INPUT_PERIOD(N,XMPIELRANK,XMPINRANK,XMPIE,XMPIN,IELEM,INODE,IMAXN,IMAXE,IBOUND,IMAXB,XMPINNUMBER,SCALER)
IMPLICIT NONE
!> @brief
!> This subroutine reads all the periodic data from the grid files
TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
	TYPE(NODE_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INODE
	INTEGER,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::XMPINNUMBER
	INTEGER,INTENT(IN)::N,IMAXN,IMAXE
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIE,XMPIN
	TYPE(BOUND_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IBOUND
	INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK,XMPINRANK
	INTEGER,INTENT(INOUT)::IMAXB
	REAL,INTENT(INOUT)::SCALER
	integer,allocatable,dimension(:)::nodep,nodec
	INTEGER::I,J,JI,K,LM,IEX,KMAXN,KK,KMAXE,print_out,kk2,shap,nodal,iftrue,kxk2
	INTEGER,DIMENSION(8)::IDV
	INTEGER::IT1,IT2,IT3,IT4,IT5,IT6,IT7,IT8,ITX,INX,IT55,in,out
	REAL::X,Y,Z


CALL OPEN_INPUT(N,ITT)

if (binio.eq.0)then
IF (DIMENSIONA.EQ.3)THEN
DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		READ(9,*)INX,X,Y,Z
		x=x/scaler; y=y/scaler ; Z=Z/scaler
		      IF (inoder2(J)%NUMBEROFNEIB.EQ.0)THEN		
		      
		      
			ALLOCATE(inoder(J)%CORD(1:3))
! 			  if ((j.eq.709).or.(j.eq.710).or.(j.eq.693).or.(j.eq.692))then

! 			  end if
			inoder(J)%CORD(1)=X
			inoder(J)%CORD(2)=Y
			inoder(J)%CORD(3)=Z
		      END IF
		  
		  else
		  
		  READ(9,*)
		  end if
END DO

ELSE
DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		READ(9,*)INX,X,Y
		x=x/scaler; y=y/scaler
		      IF (inoder2(J)%NUMBEROFNEIB.EQ.0)THEN		
		      
		      
			ALLOCATE(inoder(J)%CORD(1:2))
			inoder(J)%CORD(1)=X
			inoder(J)%CORD(2)=Y
		      END IF
		  Else
		  read(9,*)
		  end if
END DO





END IF
else
IF (DIMENSIONA.EQ.3)THEN
DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		READ(9)INX,X,Y,Z
		x=x/scaler; y=y/scaler ; Z=Z/scaler
		      IF (inoder2(J)%NUMBEROFNEIB.EQ.0)THEN		
		      
		      
			ALLOCATE(inoder(J)%CORD(1:3))
! 			  if ((j.eq.709).or.(j.eq.710).or.(j.eq.693).or.(j.eq.692))then

! 			  end if
			inoder(J)%CORD(1)=X
			inoder(J)%CORD(2)=Y
			inoder(J)%CORD(3)=Z
		      END IF
		  
		  else
		  
		  READ(9)INX,X,Y,Z
		  end if
END DO

ELSE
DO j=1,IMAXN


		  IF (inoder(J)%ITOR.gt.0)THEN
		READ(9)INX,X,Y
		x=x/scaler; y=y/scaler
		      IF (inoder2(J)%NUMBEROFNEIB.EQ.0)THEN		
		      
		      
			ALLOCATE(inoder(J)%CORD(1:2))
			inoder(J)%CORD(1)=X
			inoder(J)%CORD(2)=Y
		      END IF
		  Else
		  read(9)INX,X,Y
		  end if
END DO





END IF









end if



CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	 CALL CLOSE_INPUT(N,ITT)

! if (dimensiona.eq.3)then
!  deallocate(inoder2)

! end if




END SUBROUTINE READ_INPUT_PERIOD



SUBROUTINE stenprint(n)
IMPLICIT NONE
!> @brief
!> This subroutine prints the stencils at a selected position
integer,intent(in)::n
integer::i,j,k,inv,ismp,l,i1,i2,i3,i4,i5,i6,i7,i8,ixx
INTEGER::KMAXE,KK,KFK,ICPUID,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,DECOMF,KD,itarget
INTEGER::INX,M,O,P,Q,JK,ELEMENTSS
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
 





ICPUID=N
	
	



















if (nprobes.gt.0)then
DO INV=1,NPROBES
				  
				IF (PROBEI(N,inv).NE.0) THEN
			WRITE(PROC3,FMT='(I10)') INV
	!proc4=".plt"
	OUTFILE="STENCILS_"//TRIM(ADJUSTL(PROC3))//".dat"!//TRIM(ADJUSTL(PROC4))	
				
				
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	IF (BINIO.EQ.0)OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	IF (BINIO.EQ.1)OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	
	WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	
	IF (DIMENSIONA.EQ.3)THEN
	WRITE(97,*)'VARIABLES="X","Y","Z"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	ELSE
	WRITE(97,*)'VARIABLES="X","Y"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	
	END IF
	
	IF (BINIO.EQ.0)THEN
	DO I=1,IMAXN
		READ(96,*)j,X
		WRITE(97,*)X/SCALER
	END DO

	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96,*)j,X,Y
		WRITE(97,*)Y/SCALER
	END DO
	CLOSE(96)
	
	IF (DIMENSIONA.EQ.3)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96,*)j,X,Y,Z
		WRITE(97,*)Z/SCALER
	END DO
	CLOSE (97)
	CLOSE(96)
        END IF
        
        
	ELSE
	DO I=1,IMAXN
		READ(96)j,X
		WRITE(97,*)X/SCALER
	END DO

	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96)j,X,Y
		WRITE(97,*)Y/SCALER
	END DO
	CLOSE(96)
	
	IF (DIMENSIONA.EQ.3)THEN
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96)j,X,Y,Z
		WRITE(97,*)Z/SCALER
	END DO
	CLOSE (97)
	CLOSE(96)
        END IF
	
	
	END IF
	
	
	IF (BINIO.EQ.0)THEN
	
	OPEN(96,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
	IF (DIMENSIONA.EQ.3)THEN
	do i=1,imaxe
	      read(96,*)ix,I5,I6,I8,I7,I1,I2,I4,I3
	      write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
	end do
	else
	do i=1,imaxe
	      read(96,*)ix,I1,I2,I3,I4
	      write(97,*)I1,I2,I3,I4
	end do
	END IF
	close(96)
	close(97)
	ELSE
	OPEN(96,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
	IF (DIMENSIONA.EQ.3)THEN
	do i=1,imaxe
	      read(96)ix,I5,I6,I8,I7,I1,I2,I4,I3
	      write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
	end do
	else
	do i=1,imaxe
	      read(96)ix,I1,I2,I3,I4
	      write(97,*)I1,I2,I3,I4
	end do
	END IF
	close(96)
	close(97)
	
	END IF
				
				
				
				
				
				
				
				
				
				
                                    OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                                    
                                    if (binio.eq.0)OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                    if (binio.eq.1)OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                    WRITE(97,*)'TITLE="GRID"'
                                    WRITE(97,*)'FILETYPE=GRID'
                                    IF (DIMENSIONA.EQ.3)THEN
                                    WRITE(97,*)'VARIABLES="X","Y","Z"'
                                    WRITE(97,*)'Zone N=',IMAXN,',E=',1,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
                                    ELSE
                                    WRITE(97,*)'VARIABLES="X","Y"'
                                    WRITE(97,*)'Zone N=',IMAXN,',E=',1,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
                                    
                                    END IF
                                    
					    IF (BINIO.EQ.0)THEN
                                            DO I=1,IMAXN
                                                    READ(96,*)j,X
                                                    WRITE(97,*)X/SCALER
                                            END DO

                                            CLOSE(96)
                                            OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96,*)j,X,Y
                                                    WRITE(97,*)Y/SCALER
                                            END DO
                                            CLOSE(96)
                                            
                                            IF (DIMENSIONA.EQ.3)THEN
                                            OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96,*)j,X,Y,Z
                                                    WRITE(97,*)Z/SCALER
                                            END DO
                                            CLOSE (97)
                                            CLOSE(96)
                                            END IF
                                            ELSE
                                            DO I=1,IMAXN
                                                    READ(96)j,X
                                                    WRITE(97,*)X/SCALER
                                            END DO

                                            CLOSE(96)
                                            OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96)j,X,Y
                                                    WRITE(97,*)Y/SCALER
                                            END DO
                                            CLOSE(96)
                                            
                                            IF (DIMENSIONA.EQ.3)THEN
                                            OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96)j,X,Y,Z
                                                    WRITE(97,*)Z/SCALER
                                            END DO
                                            CLOSE (97)
                                            CLOSE(96)
                                            END IF
                                            
                                            END IF
                                            
                    
											



                                            IF (BINIO.EQ.0)THEN
                                            OPEN(96,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
                                            OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                                            IF (DIMENSIONA.EQ.3)THEN
                                            do i=1,imaxe
                                                read(96,*)ix,I5,I6,I8,I7,I1,I2,I4,I3
                                                if (ielem(n,PROBEI(N,inv))%ihexgl.eq.i)then
                                                write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
                                                end if
                                            end do
                                            else
                                            do i=1,imaxe
                                                read(96,*)ix,I1,I2,I3,I4
                                                 if (ielem(n,PROBEI(N,inv))%ihexgl.eq.i)then
                                                write(97,*)I1,I2,I3,I4
                                                end if
                                            end do
                                            END IF
                                            close(96)
                                            close(97)
					     ELSE
					      OPEN(96,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                            OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                                            IF (DIMENSIONA.EQ.3)THEN
                                            do i=1,imaxe
                                                read(96)ix,I5,I6,I8,I7,I1,I2,I4,I3
                                                if (ielem(n,PROBEI(N,inv))%ihexgl.eq.i)then
                                                write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
                                                end if
                                            end do
                                            else
                                            do i=1,imaxe
                                                read(96)ix,I1,I2,I3,I4
                                                 if (ielem(n,PROBEI(N,inv))%ihexgl.eq.i)then
                                                write(97,*)I1,I2,I3,I4
                                                end if
                                            end do
                                            END IF
                                            close(96)
                                            close(97)
					     
					     
					     END IF
                                    
                                    
                                    
                                    
                                    
				
				      DO ISMP=1,typesten
					ELEMENTSS=0
					
					if ((ismp.eq.1).or.(ees.ne.5))then
					itarget=ielem(n,PROBEI(N,inv))%inumneighbours
					else
					itarget=numneighbours2
					end if
        
					
					
					
					DO  L=2,itarget
                                            IF (ILOCALSTENCIL(N,PROBEI(N,inv),ISMP,L).GT.0)THEN
                                                     ELEMENTSS= ELEMENTSS+1 
                                            END IF
                                        END DO
                                        
                                            IF (ELEMENTSS+1.EQ.itarget)THEN
                                                                                     
                                            OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                                    IF (BINIO.EQ.0)OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                    IF (BINIO.EQ.1)OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                    WRITE(97,*)'TITLE="GRID"'
                                    WRITE(97,*)'FILETYPE=GRID'
                                    IF (DIMENSIONA.EQ.3)THEN
                                    WRITE(97,*)'VARIABLES="X","Y","Z"'
                                    WRITE(97,*)'Zone N=',IMAXN,',E=',ELEMENTSS,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
                                    ELSE
                                    WRITE(97,*)'VARIABLES="X","Y"'
                                    WRITE(97,*)'Zone N=',IMAXN,',E=',ELEMENTSS,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
                                    
                                    END IF
                                    
					    IF (BINIO.EQ.0)THEN
                                            DO I=1,IMAXN
                                                    READ(96,*)j,X
                                                    WRITE(97,*)X/SCALER
                                            END DO

                                            CLOSE(96)
                                            OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96,*)j,X,Y
                                                    WRITE(97,*)Y/SCALER
                                            END DO
                                            CLOSE(96)
                                            
                                            IF (DIMENSIONA.EQ.3)THEN
                                            OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96,*)j,X,Y,Z
                                                    WRITE(97,*)Z/SCALER
                                            END DO
                                            CLOSE (97)
                                            CLOSE(96)
                                            END IF
                                            ELSE
                                            DO I=1,IMAXN
                                                    READ(96)j,X
                                                    WRITE(97,*)X/SCALER
                                            END DO

                                            CLOSE(96)
                                            OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96)j,X,Y
                                                    WRITE(97,*)Y/SCALER
                                            END DO
                                            CLOSE(96)
                                            
                                            IF (DIMENSIONA.EQ.3)THEN
                                            OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                            DO I=1,IMAXN
                                                    READ(96)j,X,Y,Z
                                                    WRITE(97,*)Z/SCALER
                                            END DO
                                            CLOSE (97)
                                            CLOSE(96)
                                            END IF
                                            
                                            
                                            END IF
                                            
                                                                                        
                                            DO  L=2,itarget
                                            	IF (ILOCALSTENCIL(N,PROBEI(N,inv),ISMP,L).GT.0)THEN
                                            	IF (BINIO.EQ.0)THEN
                                            	 OPEN(96,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
                                            OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                                            
                                                    IF (DIMENSIONA.EQ.3)THEN
                                                        do i=1,imaxe
                                                            read(96,*)ix,I5,I6,I8,I7,I1,I2,I4,I3
                                                            IF (I.EQ.ILOCALSTENCIL(N,PROBEI(N,inv),ISMP,L))THEN
                                                            write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
                                                            end if
                                                        end do
                                                    else
                                                        do i=1,imaxe
                                                            read(96,*)ix,I1,I2,I3,I4
                                                            IF (I.EQ.ILOCALSTENCIL(N,PROBEI(N,inv),ISMP,L))THEN
                                                            write(97,*)I1,I2,I3,I4
                                                            end if
                                                        end do
                                                    END IF
                                            close(96)
                                            close(97)
						ELSE
						 OPEN(96,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                            OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                                            
                                                    IF (DIMENSIONA.EQ.3)THEN
                                                        do i=1,imaxe
                                                            read(96)ix,I5,I6,I8,I7,I1,I2,I4,I3
                                                            IF (I.EQ.ILOCALSTENCIL(N,PROBEI(N,inv),ISMP,L))THEN
                                                            write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
                                                            end if
                                                        end do
                                                    else
                                                        do i=1,imaxe
                                                            read(96)ix,I1,I2,I3,I4
                                                            IF (I.EQ.ILOCALSTENCIL(N,PROBEI(N,inv),ISMP,L))THEN
                                                            write(97,*)I1,I2,I3,I4
                                                            end if
                                                        end do
                                                    END IF
                                            close(96)
                                            close(97)
						
						
						END IF
                                                END IF
                                            END DO
					   
					  end if
					END DO
				      
				   END IF
				  
				  
				  
			    END DO
 			  

end if

 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
end subroutine stenprint









Subroutine WallDistance(N,ielem,imaxe,XMPIELRANK)
!> @brief
!> This subroutine establishes the wall distance for every cell in 3D

Implicit None

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK !,XMPINRANK
INTEGER,INTENT(IN)::N,IMAXE
Integer :: countwall,KmaxE,Countwallglobal,i,l,icpu,doyouhavewall,howmanyhavewall,counterall,j
Integer :: counterall2,wall1,wall2,wall3,wall4,wall5,wall6,IOY,wl1,wl2,wl3,wl4,k
Real,allocatable,dimension(:,:) :: WallElemArrayCord,WallElemArrayCordGlobal
Real :: Distance
CHARACTER(LEN=12)::BNDFILE,VRTFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Type(Wallboundary),Allocatable,Dimension(:):: Wallbnd

Type(NodesWall),Allocatable,Dimension(:):: Wallvrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! How many element for this block have wall
KMAXE=XMPIELRANK(N)
 countwall = 0
 Do i=1,KmaxE
    if (ielem(n,i)%interior.eq.1)then
    Do l=1,IELEM(N,I)%IFCA
	if (ielem(n,i)%ibounds(l).gt.0)then
	      if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.99))then
	  countwall = countwall + 1
	      end if
	      
	End If
    End Do
    end if
End Do
! How many for all cpu blocks
CALL MPI_ALLREDUCE(countwall,Countwallglobal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!  Print*,Countwallglobal,'total wall elements',N
Allocate(Wallbnd(Countwallglobal))

 countwall = 0
BNDFILE='GRID.bnd'
VRTFILE='GRID.vrt'
IF (BINIO.EQ.0)OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
IF (BINIO.EQ.1)OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
! Find the nodes id of the wall globally from the bnd file
 countwall = 0
IF (BINIO.EQ.0)THEN
DO I=1,IMAXB
	READ(10,*)wall1,wall2,wall3,wall4,wall5,wall6
	    If (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		Wallbnd(countwall)%wbid = wall1 ; Wallbnd(countwall)%wb1 = wall2 ; Wallbnd(countwall)%wb2 = wall3 
		Wallbnd(countwall)%wb3 = wall4  ;Wallbnd(countwall)%wb4 = wall5  
	    End If
END DO
ELSE
DO I=1,IMAXB
	READ(10)wall1,wall2,wall3,wall4,wall5,wall6
	    If (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		Wallbnd(countwall)%wbid = wall1 ; Wallbnd(countwall)%wb1 = wall2 ; Wallbnd(countwall)%wb2 = wall3 
		Wallbnd(countwall)%wb3 = wall4  ;Wallbnd(countwall)%wb4 = wall5  
	    End If
END DO

END IF
 CLOSE(10)

!  Print*,countwall,'total wall elements2',N
Allocate(Wallvrt(IMAXN))
IF (BINIO.EQ.0)OPEN(11,FILE=VRTFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
IF (BINIO.EQ.1)OPEN(11,FILE=VRTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! Read Node coordinates and store
IF (BINIO.EQ.0)THEN
Do i = 1,imaxn
  Read(11,*)Wallvrt(i)%ID,Wallvrt(i)%wnx,Wallvrt(i)%wny,Wallvrt(i)%wnz
  Wallvrt(i)%wnx=Wallvrt(i)%wnx/SCALER
  Wallvrt(i)%wnY=Wallvrt(i)%wnY/SCALER
  Wallvrt(i)%wnZ=Wallvrt(i)%wnZ/SCALER
End Do
ELSE
Do i = 1,imaxn
  Read(11)Wallvrt(i)%ID,Wallvrt(i)%wnx,Wallvrt(i)%wny,Wallvrt(i)%wnz
  Wallvrt(i)%wnx=Wallvrt(i)%wnx/SCALER
  Wallvrt(i)%wnY=Wallvrt(i)%wnY/SCALER
  Wallvrt(i)%wnZ=Wallvrt(i)%wnZ/SCALER
End Do


END IF
  Close(11)
! Compute and store wall face centers globally
Do i = 1,Countwallglobal

	wl1 = Wallbnd(i)%wb1 ;wl2 = Wallbnd(i)%wb2 ;wl3 = Wallbnd(i)%wb3 ;wl4 = Wallbnd(i)%wb4 ;
	IF (Wl4.EQ.Wl3)THEN
	WALLBND(I)%WALLX=(wallvrt(wl1)%wnx+Wallvrt(wl2)%wnx+wallvrt(wl3)%wnx)/3.0
	wallbnd(i)%wally=(wallvrt(wl1)%wny+wallvrt(wl2)%wny+wallvrt(wl3)%wny)/3.0
	wallbnd(i)%wallz=(wallvrt(wl1)%wnz+wallvrt(wl2)%wnz+wallvrt(wl3)%wnz)/3.0


	ELSE
	Wallbnd(i)%Wallx = (Wallvrt(wl1)%wnx + Wallvrt(wl2)%wnx + Wallvrt(wl3)%wnx + Wallvrt(wl4)%wnx)/4.0
	Wallbnd(i)%Wally = (Wallvrt(wl1)%wny + Wallvrt(wl2)%wny + Wallvrt(wl3)%wny + Wallvrt(wl4)%wny)/4.0
	Wallbnd(i)%Wallz = (Wallvrt(wl1)%wnz + Wallvrt(wl2)%wnz + Wallvrt(wl3)%wnz + Wallvrt(wl4)%wnz)/4.0
	END IF
End Do
! Find distance from element barycenter to the nearest wall for this block.
KMAXE=XMPIELRANK(N)
Do i=1,KmaxE
    Distance=TOLBIG 
	Do k = 1,Countwallglobal
	      If ( Distance .gt. (sqrt(((Wallbnd(k)%Wallx-Ielem(N,i)%xxc)**2) &
				     + ((Wallbnd(k)%Wally-Ielem(N,i)%yyc)**2)&
				     + ((Wallbnd(k)%Wallz-Ielem(N,i)%zzc)**2)))) Then
		Distance = sqrt(((Wallbnd(k)%Wallx-Ielem(N,i)%xxc)**2)&
			       +((Wallbnd(k)%Wally-Ielem(N,i)%yyc)**2)&
			       +((Wallbnd(k)%Wallz-Ielem(N,i)%zzc)**2))
		Ielem(N,i)%WallDist = Distance
		 IF (IELEM(N,I)%WALLDIST.LT.HYBRIDIST)THEN
		    IELEM(N,I)%HYBRID=1
		 END IF
		
	      End If
	End Do
End Do

DeAllocate(Wallbnd)
DeAllocate(Wallvrt)

End Subroutine


Subroutine WallDistance2d(N,ielem,imaxe,XMPIELRANK)
!> @brief
!> This subroutine establishes the wall distance for every cell in 2D

Implicit None

TYPE(ELEMENT_NUMBER),ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IELEM
INTEGER,ALLOCATABLE,DIMENSION(:),INTENT(IN)::XMPIELRANK !,XMPINRANK
INTEGER,INTENT(IN)::N,IMAXE
Integer :: countwall,KmaxE,Countwallglobal,i,l,icpu,doyouhavewall,howmanyhavewall,counterall,j
Integer :: counterall2,wall1,wall2,wall3,wall4,wall5,wall6,IOY,wl1,wl2,wl3,wl4,k
Real,allocatable,dimension(:,:) :: WallElemArrayCord,WallElemArrayCordGlobal
Real :: Distance
CHARACTER(LEN=12)::BNDFILE,VRTFILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Type(Wallboundary),Allocatable,Dimension(:):: Wallbnd

Type(NodesWall),Allocatable,Dimension(:):: Wallvrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! How many element for this block have wall
KMAXE=XMPIELRANK(N)
 countwall = 0
 Do i=1,KmaxE
    if (ielem(n,i)%interior.eq.1)then
    Do l=1,IELEM(N,I)%IFCA
	if (ielem(n,i)%ibounds(l).gt.0)then
	      if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.99))then
	  countwall = countwall + 1
	      end if
	       
	End If
    End Do
    end if
End Do
! How many for all cpu blocks
CALL MPI_ALLREDUCE(countwall,Countwallglobal,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

!  Print*,Countwallglobal,'total wall elements',N
Allocate(Wallbnd(Countwallglobal))

 countwall = 0
BNDFILE='GRID.bnd'
VRTFILE='GRID.vrt'
IF (BINIO.EQ.0)OPEN(10,FILE=BNDFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
IF (BINIO.EQ.1)OPEN(10,FILE=BNDFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ',IOSTAT=IOY)
! Find the nodes id of the wall globally from the bnd file
 countwall = 0
IF (BINIO.EQ.0)THEN
DO I=1,IMAXB
	READ(10,*)wall1,wall2,wall3,wall4,wall5,wall6
	    If (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		Wallbnd(countwall)%wbid = wall1 ; Wallbnd(countwall)%wb1 = wall2 ; Wallbnd(countwall)%wb2 = wall3 
		Wallbnd(countwall)%wb3 = wall4  ;Wallbnd(countwall)%wb4 = wall5  
	    End If
END DO
ELSE
DO I=1,IMAXB
	READ(10)wall1,wall2,wall3,wall4,wall5,wall6
	    If (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		Wallbnd(countwall)%wbid = wall1 ; Wallbnd(countwall)%wb1 = wall2 ; Wallbnd(countwall)%wb2 = wall3 
		Wallbnd(countwall)%wb3 = wall4  ;Wallbnd(countwall)%wb4 = wall5  
	    End If
END DO

END IF
 CLOSE(10)

 
Allocate(Wallvrt(IMAXN))
IF (BINIO.EQ.0)OPEN(11,FILE=VRTFILE,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
IF (BINIO.EQ.1)OPEN(11,FILE=VRTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! Read Node coordinates and store
IF (BINIO.EQ.0)THEN
Do i = 1,imaxn
  Read(11,*)Wallvrt(i)%ID,Wallvrt(i)%wnx,Wallvrt(i)%wny
End Do
ELSE
Do i = 1,imaxn
  Read(11)Wallvrt(i)%ID,Wallvrt(i)%wnx,Wallvrt(i)%wny
End Do

END IF
  Close(11)
! Compute and store wall face centers globally
Do i = 1,Countwallglobal

	wl1 = Wallbnd(i)%wb1 ;wl2 = Wallbnd(i)%wb2 ;wl3 = Wallbnd(i)%wb3 ;wl4 = Wallbnd(i)%wb4 ;
	
	WALLBND(I)%WALLX=(wallvrt(wl1)%wnx+Wallvrt(wl2)%wnx)/2.0
	wallbnd(i)%wally=(wallvrt(wl1)%wny+wallvrt(wl2)%wny)/2.0
	


	
End Do


! Find distance from element barycenter to the nearest wall for this block.
KMAXE=XMPIELRANK(N)
Do i=1,KmaxE
    Distance=TOLBIG 
	Do k = 1,Countwallglobal
	      If ( Distance .gt. (sqrt(((Wallbnd(k)%Wallx-Ielem(N,i)%xxc)**2) &
				     + ((Wallbnd(k)%Wally-Ielem(N,i)%yyc)**2)))) Then
		Distance = sqrt(((Wallbnd(k)%Wallx-Ielem(N,i)%xxc)**2)&
			       +((Wallbnd(k)%Wally-Ielem(N,i)%yyc)**2))
		Ielem(N,i)%WallDist = Distance
		
		 IF (IELEM(N,I)%WALLDIST.LT.HYBRIDIST)THEN
		    IELEM(N,I)%HYBRID=1
		 END IF
		
	      End If
	End Do

End Do

DeAllocate(Wallbnd)
DeAllocate(Wallvrt)

End Subroutine


SUBROUTINE OUTWRITEGRIDBs
!> @brief
!> This subroutine writes the wall surface grid in tecplotm binary format for 3D meshes
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
INTEGER,DIMENSION(70)::IVALID
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog
character*1 NULCHAR
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
     

INQUIRE (FILE='SURF.plt',EXIST=HEREv)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do

END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
deallocate(inog)
end if








else


if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2

NullPtr = 0
      Debug   = 0
      FileType = 1
      VIsDouble = 1
      IMax    = igf2
      JMax    = totwalls
      KMax    = 0
      ZoneType = 3
      SolTime = 360.0
      StrandID = 0
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
NULCHAR = CHAR(0)


ierr =  TecIni112('SIMPLE DATASET'//NULCHAR, &
                    'X Y Z'//NULCHAR, &
                    'SURF.plt'//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)


 ierr= TecZne112('WALLS'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    ShrConn)




	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
    allocate(zbin(igf2))
	IF (BINIO.EQ.0)OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	IF (BINIO.EQ.1)OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	igf2=0
	IF (BINIO.EQ.0)THEN
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	Zbin(igf2)=Z/scaler
	end if
	END DO
	ELSE
	 DO I=1,IMAXN
	READ(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	Zbin(igf2)=Z/scaler
	end if
	END DO
	
	END IF

    CLOSE(96)


    ierr = TECDAT112(igf2,xbin,1)  !!! why not xbin instead of xbin(1) ??
   
    ierr = TECDAT112(igf2,ybin,1)

     ierr = TECDAT112(igf2,zbin,1)
   
    deallocate(xbin,YBIN,zbin)
    
    IF (BINIO.EQ.0)OPEN(98,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
    IF (BINIO.EQ.1)OPEN(98,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		IF (BINIO.EQ.0)THEN
		DO K=1,iMAXB
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		ELSE
		DO K=1,iMAXB
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		
		END IF
	    
 		close(98)
		ierr = TECNOD112(icon)
		!ierr = TECNOD112(icon)
		deallocate(icon)	
		deallocate(inog)
         
        
  ierr = TECEND112()


END IF





	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	end if
	

	
	

END SUBROUTINE OUTWRITEGRIDBs



SUBROUTINE OUTWRITEGRIDs
!> @brief
!> This subroutine writes the wall surface grid in tecplotm ascii format for 3D meshes
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
INTEGER,DIMENSION(70)::IVALID
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
     

INQUIRE (FILE='SURF.plt',EXIST=HEREv)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then



IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')

allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do

END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
deallocate(inog)
end if








else


if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do

END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
	OPEN(97,FILE='SURF.plt',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	WRITE(97,*)'VARIABLES="X","Y","Z"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'







	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
    allocate(zbin(igf2))
	IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	Zbin(igf2)=Z/scaler
	end if
	END DO

    CLOSE(96)
      ELSE
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	Zbin(igf2)=Z/scaler
	end if
	END DO

    CLOSE(96)
      
      
      END IF

    
	WRITE(97,*)XBIN(1:IGF2)
	WRITE(97,*)YBIN(1:IGF2)
	WRITE(97,*)ZBIN(1:IGF2)
    
   
   
    deallocate(xbin,YBIN,zbin)
    IF (BINIO.EQ.0)OPEN(98,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
    IF (BINIO.EQ.1)OPEN(98,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		IF (BINIO.EQ.0)THEN
		DO K=1,iMAXB
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		ELSE
		DO K=1,iMAXB
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		
		END IF
	    
 		close(98)
 		
 		DO I=1,IGF2
 		WRITE(97,*)ICON(1:4,I)
 		END DO
 		
	
		deallocate(icon)	
		deallocate(inog)
         
        CLOSE(97)
  
END IF





	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	end if
	

	
	

END SUBROUTINE OUTWRITEGRIDs



SUBROUTINE OUTWRITEGRIDs2D
!> @brief
!> This subroutine writes the wall surface grid in tecplot ascii format for 2D meshes
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
INTEGER,DIMENSION(70)::IVALID
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
     

INQUIRE (FILE='SURF.plt',EXIST=HEREv)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
deallocate(inog)
end if








else


if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do

END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
OPEN(97,FILE='SURF.plt',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	WRITE(97,*)'VARIABLES="X","Y"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FELINE,','DATAPACKING = BLOCK'







	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
	IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96,*)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	END DO

    CLOSE(96)
      ELSE
      OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96,*)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	END DO

    CLOSE(96)
      
      END IF

    
	WRITE(97,*)XBIN(1:IGF2)
	WRITE(97,*)YBIN(1:IGF2)
	
    
   
   
    deallocate(xbin,YBIN)
    IF (BINIO.EQ.0)OPEN(98,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
    IF (BINIO.EQ.1)OPEN(98,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(2,totwalls))
    icon=0
     cv=0
		igf2=0
		IF (BINIO.EQ.0)THEN
		DO K=1,iMAXB
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  
    
		end if
		END DO
		ELSE
		DO K=1,iMAXB
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  
    
		end if
		END DO
		
		END IF
	    
 		close(98)
 		
 		DO I=1,IGF2
 		WRITE(97,*)ICON(1:2,I)
 		END DO
 		
	
		deallocate(icon)	
		deallocate(inog)
         
        CLOSE(97)
  
END IF





	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	end if
	

	
	

END SUBROUTINE OUTWRITEGRIDs2D



SUBROUTINE OUTWRITEGRIDBs2d
!> @brief
!> This subroutine writes the wall surface grid in tecplot binary format for 2D meshes
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
INTEGER,DIMENSION(70)::IVALID
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
     

INQUIRE (FILE='SURF.plt',EXIST=HEREv)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do

END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
deallocate(inog)
end if








else


if ((n.eq.0).and.(totwalls.gt.0))then
IF (BINIO.EQ.0)OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
IF (BINIO.EQ.1)OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
IF (BINIO.EQ.0)THEN
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
ELSE
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do

END IF
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2

NullPtr = 0
      Debug   = 0
      FileType = 1
      VIsDouble = 1
      IMax    = igf2
      JMax    = totwalls
      KMax    = 0
      ZoneType = 1
      SolTime = 360.0
      StrandID = 0
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0
NULCHAR = CHAR(0)


ierr =  TecIni112('SIMPLE DATASET'//NULCHAR, &
                    'X Y'//NULCHAR, &
                    'SURF.plt'//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)


 ierr= TecZne112('WALLS'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Null, &
                    Null, &
                    ShrConn)




	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
   
	IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96,*)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	END DO

    CLOSE(96)
      ELSE
      OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	END DO

    CLOSE(96)
      
      END IF

    ierr = TECDAT112(igf2,xbin,1)  !!! why not xbin instead of xbin(1) ??
   
    ierr = TECDAT112(igf2,ybin,1)

    
   
    deallocate(xbin,YBIN)
    IF (BINIO.EQ.0)OPEN(98,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
    IF (BINIO.EQ.1)OPEN(98,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	  allocate(icon(2,totwalls))
    icon=0
     cv=0
		igf2=0
		IF (BINIO.EQ.0)THEN
		DO K=1,iMAXB
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		 
    
		end if
    !cv=cv+4
        	
		END DO
		ELSE
		DO K=1,iMAXB
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		 
    
		end if
    !cv=cv+4
        	
		END DO
		
		END IF
	    
 		close(98)
		ierr = TECNOD112(icon)
		!ierr = TECNOD112(icon)
		deallocate(icon)	
		deallocate(inog)
         
        
  ierr = TECEND112()


END IF





	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	end if
	

	
	

END SUBROUTINE OUTWRITEGRIDBs2d


SUBROUTINE OUTWRITEGRID(N)
!> @brief
!> This subroutine writes the 3D grid file in tecplot ASCII format
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,i6,i7,i8,DECOMF,KD
INTEGER::INX,I,K,J,M,O,P,Q,JK,IXXFF
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR




if (n.eq.0)then 

ICPUID=N
	KK=0
	WRITE(PROC3,FMT='(I10)') IXXFF
	OUTFILE='GRID.dat'
! 	OUTFILE='OUT.'//TRIM(ADJUSTL(PROC3))
	
	
	IF (BINIO.EQ.0)THEN
	
	
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	
	
	WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	WRITE(97,*)'VARIABLES="X","Y","Z"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	
	
	
	DO I=1,IMAXN
		READ(96,*)j,X
		WRITE(97,*)X/SCALER
	END DO

	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96,*)j,X,Y
		WRITE(97,*)Y/SCALER
	END DO
	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96,*)j,X,Y,Z
		WRITE(97,*)Z/SCALER
	END DO
	CLOSE (97)
	CLOSE(96)

	ELSE
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	
	
	WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	WRITE(97,*)'VARIABLES="X","Y","Z"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	
	
	
	DO I=1,IMAXN
		READ(96)j,X
		WRITE(97,*)X/SCALER
	END DO

	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96)j,X,Y
		WRITE(97,*)Y/SCALER
	END DO
	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96)j,X,Y,Z
		WRITE(97,*)Z/SCALER
	END DO
	CLOSE (97)
	CLOSE(96)
	
	
	END IF
	IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
	do i=1,imaxe
	      read(96,*)ix,I5,I6,I8,I7,I1,I2,I4,I3
	      write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
	end do
	close(96)
	close(97)
	ELSE
	OPEN(96,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
	do i=1,imaxe
	      read(96)ix,I5,I6,I8,I7,I1,I2,I4,I3
	      write(97,*)I6,I2,I1,I5,I8,I4,I3,I7
	end do
	close(96)
	close(97)
	
	END IF
END IF
END SUBROUTINE OUTWRITEGRID


SUBROUTINE OUTWRITEGRID2d(N)
!> @brief
!> This subroutine writes the 2D grid file in tecplot ASCII format
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,i6,i7,i8,DECOMF,KD
INTEGER::INX,I,K,J,M,O,P,Q,JK,IXXFF
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR




if (n.eq.0)then 

ICPUID=N
	KK=0
	WRITE(PROC3,FMT='(I10)') IXXFF
	OUTFILE='GRID.dat'
! 	OUTFILE='OUT.'//TRIM(ADJUSTL(PROC3))
	
	
	IF (BINIO.EQ.0)THEN
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	
	
	WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	WRITE(97,*)'VARIABLES="X","Y"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	
	
	
	DO I=1,IMAXN
		READ(96,*)j,X
		WRITE(97,*)X/SCALER
	END DO

	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96,*)j,X,Y
		WRITE(97,*)Y/SCALER
	END DO
	CLOSE(96)
	
	CLOSE (97)
	ELSE
	
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	
	
	WRITE(97,*)'TITLE="GRID"'
	WRITE(97,*)'FILETYPE=GRID'
	WRITE(97,*)'VARIABLES="X","Y"'
	WRITE(97,*)'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	
	
	
	DO I=1,IMAXN
		READ(96)j,X
		WRITE(97,*)X/SCALER
	END DO

	CLOSE(96)
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	DO I=1,IMAXN
		READ(96)j,X,Y
		WRITE(97,*)Y/SCALER
	END DO
	CLOSE(96)
	
	CLOSE (97)

	END IF	
	IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
	do i=1,imaxe
	      read(96,*)ix,I1,I2,I3,I4
	      write(97,*)I1,I2,I3,I4
	end do
	close(96)
	close(97)
	ELSE
	OPEN(96,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
	do i=1,imaxe
	      read(96)ix,I1,I2,I3,I4
	      write(97,*)I1,I2,I3,I4
	end do
	close(96)
	close(97)
	
	END IF
	
END IF
END SUBROUTINE OUTWRITEGRID2d


SUBROUTINE OUTWRITE3vSb
!> @brief
!> This subroutine writes the 3D surface solution file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
      integer::iconsidered,facex
      REAL::SHEAR_TEMP


 KMAXE=XMPIELRANK(N)
! ! 
! DUMG=TOTIW
! call mpi_barrier(mpi_comm_world,IERROR)
! 
! CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
! IMAXP=DUML
! 
! ALLOCATE(ICELL(IMAXP))
! ICELL=0
! 
! ILOOP=0
! 
! ! IF (TOTIW.GT.0)THEN
! !     
! !    DO I=1,TOTIW
! ! 	k=IBOUND_T(I)
! ! 	do i1=k,k
! ! 	DO J=1,IELEM(N,k)%IFCA
! ! 	if (ielem(n,k)%percorg(j).eq.-4)then
! ! 	iloop=iloop+1
! ! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! ! 	go to 1043
! ! 	end if
! ! 	end do
! ! 	END DO
! ! 	1043 continue
! !   ENDDO
! ! END IF
! IF (TOTIW.GT.0)THEN
! DO I=1,KMAXE
!   if (ielem(n,i)%interior.eq.1)then
! 	DO j=1,IELEM(N,I)%IFCA
! 	  if (ielem(n,i)%ibounds(J).gt.0)then
! 	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
! 		  iloop=iloop+1
! 		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
! 	      END IF
! 	  end if
! 	END DO
!    end if
! END DO
! end if
! 
! IF (N.EQ.0)THEN
! 	ALLOCATE(ICELLA(IMAXP*ISIZE))
! 	 ICELLA=0
! 
! END IF
! 
! call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)
! 
! ! if (n.eq.0)then
! ! write(1000+n,*)icella(:)
! ! 
! ! 
! ! end if
! 
! call mpi_barrier(mpi_comm_world,IERROR)
! deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)


end if



IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'solution'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar,vortex,k,omega,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar,vortex,mu,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar,vortex,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,vortex,k,omega,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,vortex,mu,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,vortex,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = itotalb
      JMax    = totwalls
      KMax    = 0
      ZoneType = 3
      !if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .GE. 11)) Then
      SolTime = T
      !else
      !SolTime = IT

      !end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(totwalls),xbin2(totwalls))
	

 ELSE
 allocate(xbin2(1))

 END IF
 
  totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  ALLOCATE(VALUESS(xmpiwall(n)))
!   end if
    call mpi_barrier(MPI_COMM_WORLD,IERROR)
   
IF (ITESTCASE.LE.2)THEN
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(1,1)
					  ENDDO
					  END IF


		      
    
    call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
!                                                 if (kkd.eq.5)then
!                                                 if (IELEM(N,IBOUND_T(I))%ishape.eq.2)then
!                                               VALUESS(I)=-1000
!                                                 else
!                                                 VALUESS(i)=IELEM(N,IBOUND_T(I))%dih(kkd)
!                                                 end if
!                                                 else
!                                                 VALUESS(i)=IELEM(N,IBOUND_T(I))%dih(kkd)
!                                                 end if
                                                
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(1,KKD)
						if ((kkd.ge.2).and.(kkd.le.4))then
						valuess(i)=U_C(IBOUND_T(I))%VAL(1,kkd)/U_C(IBOUND_T(I))%VAL(1,1)
						end if	
					  ENDDO
					  END IF
		
		     
		
		
		call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		end do
    
					IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
                                                
                                                
                                                
                                                
                                               
					  
                                            
						leftv(1:nof_Variables)=U_C(IBOUND_T(I))%VAL(1,1:nof_Variables)
						CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
						VALUESS(i)=leftv(5)
						
					  ENDDO
					  END IF
		     
    
    
    
		call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    
		  if (passivescalar.gt.0)then
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(i)=U_CT(IBOUND_T(I))%VAL(1,turbulenceequations+passivescalar)
												
					  ENDDO
					  END IF
		  
		  
		  
		  
		  	  
		  
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  
		  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(i)=ielem(n,IBOUND_T(I))%vortex(1)
												
					  ENDDO
					  END IF
		  
		  
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  
		   IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(i)=u_ct(IBOUND_T(I))%val(1,kkd)
												
					  ENDDO
					  END IF
		  
		  
		  
				  
		   call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		
		  end do
		  end if
		  
		  
		  
		  do kkd=1,3
		  
				      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						
					 ICONSIdered=IBOUND_T(I)
					 facex=IBOUND_T2(I)
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X(ICONSIdered,facex,shear_temp)
					 CASE (2)
					 call SHEAR_y(ICONSIdered,facex,shear_temp)
					 CASE(3)
					 call SHEAR_z(ICONSIdered,facex,shear_temp)
					 END SELECT
				       		VALUESS(i)=SHEAR_TEMP					
					  ENDDO
					  END IF
		  
		  
		  
		  
		  
				  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		
		  end do
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(valuelocation,out1,xbin)
  END IF
!   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS,xbin2)
!   end if











  







	
	
	

	
	

END SUBROUTINE OUTWRITE3vSb


SUBROUTINE OUTWRITE3vSb2d
!> @brief
!> This subroutine writes the 2D surface solution file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
      integer::iconsidered,facex
      REAL::SHEAR_TEMP


 KMAXE=XMPIELRANK(N)
! 
totiw=xmpiwall(n)
! call mpi_barrier(mpi_comm_world,IERROR)
! 
! CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
! IMAXP=DUML
! 
! ALLOCATE(ICELL(IMAXP))
! ICELL=0
! 
! ILOOP=0

! IF (TOTIW.GT.0)THEN
!     
!    DO I=1,TOTIW
! 	k=IBOUND_T(I)
! 	do i1=k,k
! 	DO J=1,IELEM(N,k)%IFCA
! 	if (ielem(n,k)%percorg(j).eq.-4)then
! 	iloop=iloop+1
! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! 	go to 1043
! 	end if
! 	end do
! 	END DO
! 	1043 continue
!   ENDDO
! END IF
! IF (TOTIW.GT.0)THEN
! 
! DO I=1,KMAXE
!   if (ielem(n,i)%interior.eq.1)then
! 	DO j=1,IELEM(N,I)%IFCA
! 	  if (ielem(n,i)%ibounds(J).gt.0)then
! 	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
! 		  iloop=iloop+1
! 		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
! 	      END IF
! 	  end if
! 	END DO
!    end if
! END DO
! end if
! 
! IF (N.EQ.0)THEN
! 	ALLOCATE(ICELLA(IMAXP*ISIZE))
! 	 ICELLA=0
! 
! END IF
! 
! call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

! call mpi_barrier(mpi_comm_world,IERROR)
! deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)


end if



IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'solution'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar,vortex,k,omega,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar,vortex,mu,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar,vortex,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,vortex,k,omega,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,vortex,mu,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,vortex,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = itotalb
      JMax    = totwalls
      KMax    = 0
      ZoneType = 1
      !if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      !else
      !SolTime = IT

      !end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(totwalls),xbin2(totwalls))
	

 Else
 allocate(xbin2(1))

 end if

 totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  ALLOCATE(VALUESS(xmpiwall(n)))
!   end if
    
   
IF (ITESTCASE.LE.2)THEN
		     
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(1,1)
					  ENDDO
					  END IF
    
    
     
    
    
    
     call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,4
		     IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(1,kkd)
							    if ((kkd.ge.2).and.(kkd.le.3))then
			      valuess(i)=U_C(IBOUND_T(I))%VAL(1,kkd)/U_C(IBOUND_T(I))%VAL(1,1)
			      end if	
						
						
					  ENDDO
		  END IF
		
		
		      
		
		
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		
		end do
    
		      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					  leftv(1:nof_Variables)=U_C(IBOUND_T(I))%VAL(1,1:nof_Variables)
				    CALL cons2prim(N,leftv,MP_PINFl,gammal)
				    VALUESS(i)=leftv(4)
					end do	
			      end if	
			
    
    
    
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		
    
		  if (passivescalar.gt.0)then
		  
		  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					VALUESS(i)=U_Ct(IBOUND_T(I))%VAL(1,turbulenceequations+passivescalar)
					end do	
			      end if	
		  
		  
		   
		  
		  
		  	  
		  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					VALUESS(i)=ielem(n,IBOUND_T(I))%vortex(1)
					end do	
			      end if	
		 
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		   call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		   IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					VALUESS(i)=u_ct(IBOUND_T(I))%val(1,kkd)
					end do	
			      end if
		  
		  
				  
		   call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		end do
		  end if
		  
		  
		  
		  do kkd=1,2
			    IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					  ICONSIdered=IBOUND_T(I)
					  facex=IBOUND_t2(i)
					  select case(kkd)
					 case(1)
					 
					 call SHEAR_x2d(ICONSIdered,facex,shear_temp)
					 CASE (2)
					 call SHEAR_y2d(ICONSIdered,facex,shear_temp)
					 
					 END SELECT
					VALUESS(i)=SHEAR_TEMP
					end do	
			      end if
		 
		  
				  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  
		  
		  
		  
		  
		  
		  end do
		  
		  end if
    
    
    end if





!   END IF
!    end if

    IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1)
  END IF
  
!   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS,xbin2)
!   end if
  











  







	
	
	

	
	

END SUBROUTINE OUTWRITE3vSb2d



SUBROUTINE OUTWRITE3vS
!> @brief
!> This subroutine writes the 3D surface solution file in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
integer::iconsidered,facex
      REAL::SHEAR_TEMP

 KMAXE=XMPIELRANK(N)
! 
DUMG=TOTIW
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

ILOOP=0

! IF (TOTIW.GT.0)THEN
!     
!    DO I=1,TOTIW
! 	k=IBOUND_T(I)
! 	do i1=k,k
! 	DO J=1,IELEM(N,k)%IFCA
! 	if (ielem(n,k)%percorg(j).eq.-4)then
! 	iloop=iloop+1
! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! 	go to 1043
! 	end if
! 	end do
! 	END DO
! 	1043 continue
!   ENDDO
! END IF
IF (TOTIW.GT.0)THEN
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
		  iloop=iloop+1
		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
	      END IF
	  end if
	END DO
   end if
END DO
end if

IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)

OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
end if








		
   if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if





IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  
  
   if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="solution"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 
                    
                    
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if
     ELSE
     
     
     if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      
                    
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","k","omega","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED,[13] = CELLCENTERED)'
  end if
                    
                    
                    
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","MU","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
              
                if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","k","omega","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","M","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
         
               if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED)'
  end if
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







     

Valuelocation(:)=0





 


 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(totwalls))
	VALUESA=0.0

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=0.0
    
   
IF (ITESTCASE.LE.2)THEN
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(1,1)
				END IF
		      end do
		      end if
    
    call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    WRITE(97,*)XBIN(1:TOTWALLS)
    
    END IF

     
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(1,kkd)
				if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(icount_wall)=U_C(IBOUND(N,i)%which)%VAL(1,kkd)/U_C(IBOUND(N,i)%which)%VAL(1,1)
		  end if	
				END IF
		      end do
		      end if
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
    
    
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_Variables)=U_C(IBOUND(N,i)%which)%VAL(1,1:nof_Variables)
				    CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
				    VALUESS(icount_wall)=leftv(5)
				  
				END IF
		      end do
		      end if
    
    
    
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    
		  if (passivescalar.gt.0)then
		   IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=U_Ct(IBOUND(N,i)%which)%VAL(1,turbulenceequations+passivescalar)
				  
				END IF
		      end do
		      end if
		  
		  
		  	  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=ielem(n,IBOUND(N,i)%which)%vortex(1)
				  
				END IF
		      end do
		      end if
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=u_ct(IBOUND(N,i)%which)%val(1,kkd)
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,3
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
					ICONSIdered=IBOUND(N,i)%which
					 facex=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X(ICONSIdered,facex,shear_temp)
					 CASE (2)
					 call SHEAR_y(ICONSIdered,facex,shear_temp)
					 CASE(3)
					 call SHEAR_z(ICONSIdered,facex,shear_temp)
					 END SELECT
				       		VALUESS(icount_wall)=SHEAR_TEMP		
				       
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

  IF (N.EQ.0)THEN
  
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  











  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







	
	
	

	
	

END SUBROUTINE OUTWRITE3vS


SUBROUTINE OUTWRITE3vS2d
!> @brief
!> This subroutine writes the 2D surface solution file in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
integer::iconsidered,facex
      REAL::SHEAR_TEMP

 KMAXE=XMPIELRANK(N)
! 
DUMG=TOTIW
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

ILOOP=0

! IF (TOTIW.GT.0)THEN
!     
!    DO I=1,TOTIW
! 	k=IBOUND_T(I)
! 	do i1=k,k
! 	DO J=1,IELEM(N,k)%IFCA
! 	if (ielem(n,k)%percorg(j).eq.-4)then
! 	iloop=iloop+1
! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! 	go to 1043
! 	end if
! 	end do
! 	END DO
! 	1043 continue
!   ENDDO
! END IF
IF (TOTIW.GT.0)THEN
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
		  iloop=iloop+1
		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
	      END IF
	  end if
	END DO
   end if
END DO
end if

IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)

	OPEN(97,FILE=outfile,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	




end if



IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="SOLUTION"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
 
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 
                    
                     if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     ELSE
      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED)'
  end if
!      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
!                     'Density,U,V,energy,Pressure'//NULCHAR, &
!                     out1//NULCHAR, &
!                     '.'//NULCHAR, &
!                     FileType, &
!                     Debug, &
!                     VIsDouble)
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","K","OMEGA","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED, [11] = CELLCENTERED)'
  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                   if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","M","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
               if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED)'
  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	             if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","K","OMEGA","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
! 	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
!                     'Density,U,V,energy,Pressure,vortex,k,omega,ssx,ssy'//NULCHAR, &
!                     out1//NULCHAR, &
!                     '.'//NULCHAR, &
!                     FileType, &
!                     Debug, &
!                     VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","MU","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED)'
	
  end if
              
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = itotalb
      JMax    = totwalls
      KMax    = 0
      ZoneType = 1
     ! if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      !else
      !SolTime = IT

      !end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





!  ierr= TecZne112('GRID2'//NULCHAR, &
!                     ZoneType, &
!                     IMax, &
!                     JMax, &
!                     kmax, &
!                     ICellMax, &
!                     JCellMax, &
!                     KCellMax, &
!                     SolTime, &
!                     StrandID, &
!                     ParentZn, &
!                     IsBlock, &
!                     NFConns, &
!                     FNMode, &
!                     0, &
!                     0, &
!                     0, &
!                     Null, &
!                     Valuelocation, &
!                     Null, &
!                     ShrConn)


 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(totwalls))
	VALUESA=0.0

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=0.0
    
   
IF (ITESTCASE.LE.2)THEN
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(1,1)
				END IF
		      end do
		      end if
    
    call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    WRITE(97,*)XBIN(1:TOTWALLS)
!     ierr = TECDAT112(totwalls,xbin,1)
    END IF

     
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,4
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(1,kkd)
				if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(icount_wall)=U_C(IBOUND(N,i)%which)%VAL(1,kkd)/U_C(IBOUND(N,i)%which)%VAL(1,1)
		  end if	
				END IF
		      end do
		      end if
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
    
    
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_Variables)=U_C(IBOUND(N,i)%which)%VAL(1,1:nof_Variables)
				    CALL cons2prim(N,leftv,MP_PINFl,gammal)
				    VALUESS(icount_wall)=leftv(4)
				  
				END IF
		      end do
		      end if
    
    
    
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    
		  if (passivescalar.gt.0)then
		   IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=U_Ct(IBOUND(N,i)%which)%VAL(1,turbulenceequations+passivescalar)
				  
				END IF
		      end do
		      end if
		  
		  
		  	  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=ielem(n,IBOUND(N,i)%which)%vortex(1)
				  
				END IF
		      end do
		      end if
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=u_ct(IBOUND(N,i)%which)%val(1,kkd)
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,2
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
					ICONSIdered=IBOUND(N,i)%which
					 facex=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call SHEAR_x2d(ICONSIdered,facex,shear_temp)
					 CASE (2)
					 call SHEAR_y2d(ICONSIdered,facex,shear_temp)
					 
					 END SELECT
				       		VALUESS(icount_wall)=SHEAR_TEMP		
				       
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

  IF (N.EQ.0)THEN
 
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  











  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







	
	
	

	
	

END SUBROUTINE OUTWRITE3vS2d




SUBROUTINE OUTWRITE3vbav
!> @brief
!> This subroutine writes the 3D averaged solution file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,ind1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(13))

KMAXE=XMPIELRANK(N)





if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="VOL_AVER_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	

end if
call mpi_barrier(mpi_comm_world,IERROR)

  NVAR1=11
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'R_mean,U_mean,V_mean,W_mean,P_mean,U_rms,V_rms,W_rms,UV,UW,WV'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 

	

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 5
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0

 



 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(imaxe),xbin2(imaxe))
	

 END IF

  
  ALLOCATE(VALUESS(KMAXE))
  VALUESS=ZERO
    
   
	      if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if
	      
	      

    
    
    		do kkd=1,5
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(ind1,kkd)
		  if ((kkd.ge.2).and.(kkd.le.4))then
		  VALUESS(i)=U_C(I)%VAL(ind1,kkd)/U_C(I)%VAL(ind1,1)
		  end if
		  if (kkd.eq.5)then
		   leftv(1:nof_Variables)=U_C(I)%VAL(ind1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		  
		  end if
		  
		END DO
		
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		
		end do
		do kkd=1,6

    		DO I=1,KMAXE

		  VALUESS(i)=U_C(I)%RMS(kkd)
		END DO
		
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
		END DO
    
		  
    
    
      
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1,xbin2)
  END IF
  
  DEALLOCATE (VALUESS,VARIABLES)

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)









	
	

END SUBROUTINE OUTWRITE3vbav

SUBROUTINE OUTWRITE3vb2Dav
!> @brief
!> This subroutine writes the 2D averaged solution file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,ind1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(10))

KMAXE=XMPIELRANK(N)

DUMG=KMAXE
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="VOL_AVER_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	

end if
call mpi_barrier(mpi_comm_world,IERROR)

  NVAR1=7
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'R_mean,U_mean,V_mean,P_mean,U_rms,V_rms,UV'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 

	

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 3
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0

 



 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(imaxe))
	VALUESA=ZERO

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=ZERO
    
	      if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if
	      


    
    
    		do kkd=1,4
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(ind1,kkd)
		  if ((kkd.ge.2).and.(kkd.le.3))then
		  VALUESS(i)=U_C(I)%VAL(ind1,kkd)/U_C(I)%VAL(ind1,1)
		  end if
		  if (kkd.eq.4)then
		   leftv(1:nof_Variables)=U_C(I)%VAL(ind1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(4)
		  
		  end if
		  
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
		do kkd=1,3
    		DO I=1,KMAXE
		  
		  VALUESS(i)=U_C(I)%RMS(kkd)
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		END DO
    
		  
    
    
      
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







         
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITE3vb2Dav



SUBROUTINE OUTWRITE3vav
!> @brief
!> This subroutine writes the 3D averaged solution file in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(13))

KMAXE=XMPIELRANK(N)

DUMG=KMAXE
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="VOL_AVER_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
 WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="R_mean","U_mean","V_mean","W_mean","P_mean","U_rms","V_rms","W_rms","UV","UW","WV"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED)'
end if
call mpi_barrier(mpi_comm_world,IERROR)

  
  
      
  
  
 

	

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))









 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(imaxe))
	VALUESA=ZERO

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=ZERO
    
   


    
    
    		do kkd=1,5
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(5,kkd)
		  if ((kkd.ge.2).and.(kkd.le.4))then
		  VALUESS(i)=U_C(I)%VAL(5,kkd)/U_C(I)%VAL(5,1)
		  end if
		  if (kkd.eq.5)then
		   leftv(1:nof_Variables)=U_C(I)%VAL(5,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		  
		  end if
		  
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
		do kkd=1,6
    		DO I=1,KMAXE
		  
		  VALUESS(i)=U_C(I)%RMS(kkd)
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		END DO
    
		  
    
    
      
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







         
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITE3vav

SUBROUTINE OUTWRITE3v2Dav
!> @brief
!> This subroutine writes the 2D averaged solution file in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(10))

KMAXE=XMPIELRANK(N)

DUMG=KMAXE
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then





 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="VOL_AVER_"//TRIM(ADJUSTL(PROC3))//".plt"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)
	OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')

end if
call mpi_barrier(mpi_comm_world,IERROR)


  
  IF (N.EQ.0)THEN
        WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="R_mean","U_mean","V_mean","P_mean","U_rms","V_rms","UV"'
	WRITE(97,*) 'Zone N=',IMAXN,',E=',IMAXE,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  END IF


  NVAR1=7
  
  
  
  
 

	

	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = imaxn
      JMax    = IMAXe
      KMax    = 0
      ZoneType = 3
      
      SolTime = T
      
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0

 






 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(imaxe))
	VALUESA=ZERO

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=ZERO
    
   


    
    
    		do kkd=1,4
		DO I=1,KMAXE
		  VALUESS(i)=U_C(I)%VAL(5,kkd)
		  if ((kkd.ge.2).and.(kkd.le.3))then
		  VALUESS(i)=U_C(I)%VAL(5,kkd)/U_C(I)%VAL(5,1)
		  end if
		  if (kkd.eq.4)then
		   leftv(1:nof_Variables)=U_C(I)%VAL(5,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(4)
		  
		  end if
		  
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*) XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
		do kkd=1,3
    		DO I=1,KMAXE
		  
		  VALUESS(i)=U_C(I)%RMS(kkd)
		END DO
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*) XBIN(1:IMAXE)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		END DO
    
		  
    
    
      
    
    
    

     
    
    
    
  IF (N.EQ.0)THEN
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







         
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITE3v2Dav



SUBROUTINE OUTWRITE3vSbav
!> @brief
!> This subroutine writes the 3D surface averaged solution file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall,ind1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
integer::iconsidered,facex
      REAL::SHEAR_TEMP

 KMAXE=XMPIELRANK(N)
! 

! IF (TOTIW.GT.0)THEN
!     
!    DO I=1,TOTIW
! 	k=IBOUND_T(I)
! 	do i1=k,k
! 	DO J=1,IELEM(N,k)%IFCA
! 	if (ielem(n,k)%percorg(j).eq.-4)then
! 	iloop=iloop+1
! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! 	go to 1043
! 	end if
! 	end do
! 	END DO
! 	1043 continue
!   ENDDO
! END IF




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_AV"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)


end if



IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'solution'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar,vortex,k,omega,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar,vortex,mu,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,passivescalar,vortex,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,vortex,k,omega,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,vortex,mu,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,vortex,ssx,ssy,ssz'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = itotalb
      JMax    = totwalls
      KMax    = 0
      ZoneType = 3
!       if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
!       else
!       SolTime = IT

!       end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 
  allocate(xbin(totwalls),xbin2(totwalls))
	

 END IF

  totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  ALLOCATE(VALUESS(xmpiwall(n)))
!   end if
    if (rungekutta.eq.4)then
    ind1=7
    
    
    else
    ind1=5
    
    end if
   
IF (ITESTCASE.LE.2)THEN

		      
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(ind1,1)
					  ENDDO
					  END IF
    
   call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(ind1,KKD)
						if ((kkd.ge.2).and.(kkd.le.4))then
					      valuess(i)=U_C(IBOUND_T(I))%VAL(ind1,kkd)/U_C(IBOUND_T(I))%VAL(ind1,1)
					      end if	
						
						
					  ENDDO
					  END IF
		
		
		
		
		call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		end do
    
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					        leftv(1:nof_Variables)=U_C(IBOUND_T(I))%VAL(ind1,1:nof_Variables)
						
						CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
				    VALUESS(i)=leftv(5)
						
						
						
					  ENDDO
					  END IF
		     
    
    
    
		call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    
		  if (passivescalar.gt.0)then
		  
					IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					        
				    VALUESS(i)=U_Ct(IBOUND_T(I))%VAL(ind1,turbulenceequations+passivescalar)
						
						
						
					  ENDDO
					  END IF
		  
		  
		   
		  
		  
		  	  
		  
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					        
				    VALUESS(i)=ielem(n,IBOUND_T(I))%vortex(1)
						
						
						
					  ENDDO
					  END IF
		  
		  
		  
		 
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					        
				    VALUESS(i)=u_ct(IBOUND_T(I))%val(ind1,kkd)
						
						
						
					  ENDDO
					  END IF
		  
		  
		
		  
				  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  end do
		  end if
		  
		  
		  if (ITESTCASE.eq.4)then
		  do kkd=1,3
		  
				      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
					        
				    
						
						ICONSIdered=IBOUND_T(I)
					 facex=IBOUND_T2(I)
					select case(kkd)
					 case(1)
					 
					 call shear_x_av(iconsidered,facex,shear_temp)
					 CASE (2)
					 call shear_y_av(iconsidered,facex,shear_temp)
					 CASE(3)
					 call shear_z_av(iconsidered,facex,shear_temp)
					 END SELECT
				       		VALUESS(i)=SHEAR_TEMP		
						
					  ENDDO
					  END IF
		  
		  
		  
		  
		  
				  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

  IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2)
  END IF
  
!  IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS)
!   end if
  











 







	
	
	

	
	

END SUBROUTINE OUTWRITE3vSbav


SUBROUTINE OUTWRITE3vSb2dav
!> @brief
!> This subroutine writes the 2D surface averaged solution file in tecplot binary format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
integer::iconsidered,facex
      REAL::SHEAR_TEMP

 KMAXE=XMPIELRANK(N)
 
! 
! DUMG=TOTIW
! call mpi_barrier(mpi_comm_world,IERROR)
! 
! CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
! IMAXP=DUML
! 
! ALLOCATE(ICELL(IMAXP))
! ICELL=0
! 
! ILOOP=0
! 
! ! IF (TOTIW.GT.0)THEN
! !     
! !    DO I=1,TOTIW
! ! 	k=IBOUND_T(I)
! ! 	do i1=k,k
! ! 	DO J=1,IELEM(N,k)%IFCA
! ! 	if (ielem(n,k)%percorg(j).eq.-4)then
! ! 	iloop=iloop+1
! ! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! ! 	go to 1043
! ! 	end if
! ! 	end do
! ! 	END DO
! ! 	1043 continue
! !   ENDDO
! ! END IF
! IF (TOTIW.GT.0)THEN
! DO I=1,KMAXE
!   if (ielem(n,i)%interior.eq.1)then
! 	DO j=1,IELEM(N,I)%IFCA
! 	  if (ielem(n,i)%ibounds(J).gt.0)then
! 	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
! 		  iloop=iloop+1
! 		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
! 	      END IF
! 	  end if
! 	END DO
!    end if
! END DO
! end if
! 
! IF (N.EQ.0)THEN
! 	ALLOCATE(ICELLA(IMAXP*ISIZE))
! 	 ICELLA=0
! 
! END IF
! 
! call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)
! 
! ! if (n.eq.0)then
! 
! ! 
! ! 
! ! end if
! 
! call mpi_barrier(mpi_comm_world,IERROR)
! deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_av"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)


end if



IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'solution'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
  
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     ELSE
     
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar,vortex,k,omega,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar,vortex,mu,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,passivescalar,vortex,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,vortex,k,omega,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,vortex,mu,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,energy,Pressure,vortex,ssx,ssy'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = itotalb
      JMax    = totwalls
      KMax    = 0
      ZoneType = 1
!       if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
!       else
!       SolTime = IT

!       end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





 ierr= TecZne112('GRID2'//NULCHAR, &
                    ZoneType, &
                    IMax, &
                    JMax, &
                    kmax, &
                    ICellMax, &
                    JCellMax, &
                    KCellMax, &
                    SolTime, &
                    StrandID, &
                    ParentZn, &
                    IsBlock, &
                    NFConns, &
                    FNMode, &
                    0, &
                    0, &
                    0, &
                    Null, &
                    Valuelocation, &
                    Null, &
                    ShrConn)


 allocate(xbin(totwalls),xbin2(totwalls))
	

 END IF

 totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  ALLOCATE(VALUESS(xmpiwall(n)))
!   end if
    
   
IF (ITESTCASE.LE.2)THEN

		      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(5,1)
					  ENDDO
					  END IF
		     
    
   call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,4
		      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						VALUESS(I)=U_C(IBOUND_T(I))%VAL(5,KKD)
						if ((kkd.ge.2).and.(kkd.le.3))then
						valuess(i)=U_C(IBOUND_T(I))%VAL(5,kkd)/U_C(IBOUND_T(I))%VAL(5,1)
						end if	
					  ENDDO
					  END IF
		     
		
		
		call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		end do
    
		      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						leftv(1:nof_Variables)=U_C(IBOUND_T(I))%VAL(5,1:nof_Variables)
						CALL cons2prim(N,leftv,MP_PINFl,gammal)
					  VALUESS(i)=leftv(4)
					  ENDDO
					  END IF
		      
    
    
    
		call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
    
		  if (passivescalar.gt.0)then
		    IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						
						
					  VALUESS(i)=U_CT(IBOUND_T(I))%VAL(5,turbulenceequations+passivescalar)
					  ENDDO
					  END IF
		  
		  
		  
		  
		  	  
		  
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		   IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						
						
					  VALUESS(i)=ielem(n,IBOUND_T(I))%vortex(1)
					  ENDDO
					  END IF
		
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
					  IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						
						
					  VALUESS(i)=u_ct(IBOUND_T(I))%val(5,kkd)
					  ENDDO
					  END IF
		  
		  
				  
		 call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,2
				      IF (TOTIW.GT.0)THEN
					  DO I=1,TOTIW
						ICONSIdered=IBOUND_T(I)
					 facex=IBOUND_T2(I)
					select case(kkd)
					 case(1)
					 
					 call shear_x2d_av(iconsidered,facex,shear_temp)
					 CASE (2)
					 call shear_y2d_av(iconsidered,facex,shear_temp)
					 
					 END SELECT
				       		VALUESS(i)=SHEAR_TEMP		
						
					  
					  ENDDO
					  END IF
		 
		  
				  
		  call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		ierr = TECDAT112(totwalls,xbin,1)
		END IF
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

   IF (N.EQ.0)THEN
  ierr = TECEND112()
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2)
  END IF
  
!   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS)
!   end if






	
	
	

	
	

END SUBROUTINE OUTWRITE3vSb2dav


SUBROUTINE OUTWRITE3vSav
!> @brief
!> This subroutine writes the 3D surface averaged solution file in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
integer::iconsidered,facex
      REAL::SHEAR_TEMP

 KMAXE=XMPIELRANK(N)
! 
DUMG=TOTIW
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

ILOOP=0

! IF (TOTIW.GT.0)THEN
!     
!    DO I=1,TOTIW
! 	k=IBOUND_T(I)
! 	do i1=k,k
! 	DO J=1,IELEM(N,k)%IFCA
! 	if (ielem(n,k)%percorg(j).eq.-4)then
! 	iloop=iloop+1
! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! 	go to 1043
! 	end if
! 	end do
! 	END DO
! 	1043 continue
!   ENDDO
! END IF
IF (TOTIW.GT.0)THEN
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
		  iloop=iloop+1
		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
	      END IF
	  end if
	END DO
   end if
END DO
end if

IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_av"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)

OPEN(97,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
end if








		
   if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if





IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  
  
   if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="solution"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 
                    
                    
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if
     ELSE
     
     
     if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      
                    
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","k","omega","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED,[13] = CELLCENTERED)'
  end if
                    
                    
                    
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","MU","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
              
                if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","k","omega","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","M","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
         
               if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","ssx","ssy","ssz"'
	WRITE(97,*) 'Zone N=',Itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	WRITE(97,*) '[10] = CELLCENTERED)'
  end if
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







     

Valuelocation(:)=0





 


 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(totwalls))
	VALUESA=0.0

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=0.0
    
   
IF (ITESTCASE.LE.2)THEN
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(5,1)
				END IF
		      end do
		      end if
    
    call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    WRITE(97,*)XBIN(1:TOTWALLS)
    
    END IF

     
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,5
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(5,kkd)
				if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(icount_wall)=U_C(IBOUND(N,i)%which)%VAL(5,kkd)/U_C(IBOUND(N,i)%which)%VAL(5,1)
		  end if	
				END IF
		      end do
		      end if
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
    
    
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_Variables)=U_C(IBOUND(N,i)%which)%VAL(5,1:nof_Variables)
				    CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
				    VALUESS(icount_wall)=leftv(5)
				  
				END IF
		      end do
		      end if
    
    
    
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    
		  if (passivescalar.gt.0)then
		   IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=U_Ct(IBOUND(N,i)%which)%VAL(5,turbulenceequations+passivescalar)
				  
				END IF
		      end do
		      end if
		  
		  
		  	  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=ielem(n,IBOUND(N,i)%which)%vortex(1)
				  
				END IF
		      end do
		      end if
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=u_ct(IBOUND(N,i)%which)%val(5,kkd)
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,3
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
					ICONSIdered=IBOUND(N,i)%which
					 facex=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call shear_x_av(iconsidered,facex,shear_temp)
					 CASE (2)
					 call shear_y_av(iconsidered,facex,shear_temp)
					 CASE(3)
					 call shear_z_av(iconsidered,facex,shear_temp)
					 END SELECT
				       		VALUESS(icount_wall)=SHEAR_TEMP		
				       
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

  IF (N.EQ.0)THEN
  
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  











  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







	
	
	

	
	

END SUBROUTINE OUTWRITE3vSav


SUBROUTINE OUTWRITE3vS2dav
!> @brief
!> This subroutine writes the 2D surface averaged solution file in tecplot ascii format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD,im
REAL,DIMENSION(8)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE,proc4
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNOD112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType,ILOOP
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
      REAL::SHEAR_TEMP
integer::iconsidered,facex

 KMAXE=XMPIELRANK(N)
! 
DUMG=TOTIW
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

ILOOP=0

! IF (TOTIW.GT.0)THEN
!     
!    DO I=1,TOTIW
! 	k=IBOUND_T(I)
! 	do i1=k,k
! 	DO J=1,IELEM(N,k)%IFCA
! 	if (ielem(n,k)%percorg(j).eq.-4)then
! 	iloop=iloop+1
! 	ICELL(Iloop)=IELEM(N,IBOUND_T(k))%indexi(j)
! 	go to 1043
! 	end if
! 	end do
! 	END DO
! 	1043 continue
!   ENDDO
! END IF
IF (TOTIW.GT.0)THEN
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
		  iloop=iloop+1
		    ICELL(Iloop)=ibound(n,ielem(n,i)%ibounds(j))%inum
	      END IF
	  end if
	END DO
   end if
END DO
end if

IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)




if (n.eq.0)then



 NullPtr = 0
      Debug   = 0
      FileType = 2
VIsDouble = 1

NULCHAR = CHAR(0)

WRITE(PROC3,FMT='(I10)') IT
	!proc4=".plt"
	OUTFILE="SURF_av"//TRIM(ADJUSTL(PROC3))//'.plt'
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
! 	out1=out1//CHAR(0)

	OPEN(97,FILE=outfile,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
	




end if



IF (ITESTCASE.LE.2)THEN
  NVAR1=1
  
                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="SOLUTION"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
 
  
  
 END IF
 IF (ITESTCASE.EQ.3)THEN
 NVAR1=6+PASSIVESCALAR
  if (passivescalar.gt.0)then
 
                    
                     if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     ELSE
      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED)'
  end if
!      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
!                     'Density,U,V,energy,Pressure'//NULCHAR, &
!                     out1//NULCHAR, &
!                     '.'//NULCHAR, &
!                     FileType, &
!                     Debug, &
!                     VIsDouble)
     
     
     
     END IF
 END IF
 IF (ITESTCASE.EQ.4)THEN
 NVAR1=10+PASSIVESCALAR+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	                    if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","K","OMEGA","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED, [11] = CELLCENTERED)'
  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                   if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","M","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
               if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED)'
  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	             if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","K","OMEGA","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
! 	      if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
!                     'Density,U,V,energy,Pressure,vortex,k,omega,ssx,ssy'//NULCHAR, &
!                     out1//NULCHAR, &
!                     '.'//NULCHAR, &
!                     FileType, &
!                     Debug, &
!                     VIsDouble)
              end if
              if (turbulenceequations.eq.1)then
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","MU","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	WRITE(97,*) '[9] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
                      if (n.eq.0)then
	WRITE(97,*) 'FILETYPE=SOLUTION'
	WRITE(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","SSX","SSY"'
	WRITE(97,*) 'Zone N=',ITOTALB,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	WRITE(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	WRITE(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED)'
	
  end if
              
              
              end if
	    
	    
	    
	    
	    end if
 
 
 END IF
	
	

if (n.eq.0)then
allocate (Valuelocation(nvar1))







      IMax    = itotalb
      JMax    = totwalls
      KMax    = 0
      ZoneType = 1
!       if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
!       else
!       SolTime = IT

!       end if
      StrandID = 1
      ParentZn = 0
      IsBlock = 1
      ICellMax = 0
      JCellMax = 0
      KCellMax = 0
      NFConns = 0
      FNMode = 0
      ShrConn = 0

Valuelocation(:)=0





!  ierr= TecZne112('GRID2'//NULCHAR, &
!                     ZoneType, &
!                     IMax, &
!                     JMax, &
!                     kmax, &
!                     ICellMax, &
!                     JCellMax, &
!                     KCellMax, &
!                     SolTime, &
!                     StrandID, &
!                     ParentZn, &
!                     IsBlock, &
!                     NFConns, &
!                     FNMode, &
!                     0, &
!                     0, &
!                     0, &
!                     Null, &
!                     Valuelocation, &
!                     Null, &
!                     ShrConn)


 ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(totwalls))
	VALUESA=0.0

 END IF

  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  ALLOCATE(VALUESS(imaxp))
  VALUESS=0.0
    
   
IF (ITESTCASE.LE.2)THEN
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(5,1)
				END IF
		      end do
		      end if
    
    call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    WRITE(97,*)XBIN(1:TOTWALLS)
!     ierr = TECDAT112(totwalls,xbin,1)
    END IF

     
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    END IF
    
    IF (ITESTCASE.ge.3)THEN
		do kkd=1,4
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c(IBOUND(N,i)%which)%val(5,kkd)
				if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(icount_wall)=U_C(IBOUND(N,i)%which)%VAL(5,kkd)/U_C(IBOUND(N,i)%which)%VAL(5,1)
		  end if	
				END IF
		      end do
		      end if
		
		
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		 
		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		end do
    
    
		      IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_Variables)=U_C(IBOUND(N,i)%which)%VAL(5,1:nof_Variables)
				    CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
				    VALUESS(icount_wall)=leftv(4)
				  
				END IF
		      end do
		      end if
    
    
    
		call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		IF (N.EQ.0)THEN
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		WRITE(97,*)XBIN(1:TOTWALLS)
		END IF

		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    
		  if (passivescalar.gt.0)then
		   IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=U_Ct(IBOUND(N,i)%which)%VAL(5,turbulenceequations+passivescalar)
				  
				END IF
		      end do
		      end if
		  
		  
		  	  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=ielem(n,IBOUND(N,i)%which)%vortex(1)
				  
				END IF
		      end do
		      end if
		  
		  
		  
! 		  DO I=1,KMAXE
! 		      VALUESS(i)=ielem(n,i)%vortex(1)
! 		  END DO
		  
		  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
				       VALUESS(icount_wall)=u_ct(IBOUND(N,i)%which)%val(5,kkd)
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,2
		  
		  IF (TOTIW.GT.0)THEN
		      icount_wall=0	
		      DO I=1,n_boundaries
			      if ((ibound(n,i)%icode.eq.4).and.(IBOUND(N,i)%which.gt.0)) then
				    icount_wall=icount_wall+1
					ICONSIdered=IBOUND(N,i)%which
					 facex=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call shear_x2d_av(iconsidered,facex,shear_temp)
					 CASE (2)
					 call shear_y2d_av(iconsidered,facex,shear_temp)
					 
					 END SELECT
				       		VALUESS(icount_wall)=SHEAR_TEMP		
				       
				  
				END IF
		      end do
		      end if
		  
				  
		  call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

		  IF (N.EQ.0)THEN
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  WRITE(97,*)XBIN(1:TOTWALLS)
		  END IF

		  
		  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   END IF
!    end if

  IF (N.EQ.0)THEN
 
  DEALLOCATE(XBIN,VALUESA,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  











  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







	
	
	

	
	

END SUBROUTINE OUTWRITE3vS2dav






SUBROUTINE GRID_WRITE
!> @brief
!> This subroutine calls the appropriate grid writing subroutine based on the settings
IMPLICIT NONE



IF (TECPLOT.EQ.1)THEN		!BINARY TECPLOT
  if (dimensiona.eq.3)then
  call outwritegridb
  ELSE
  call outwritegridb2D
  END IF

END IF
IF (TECPLOT.EQ.0)THEN		!ASCII TECPLOT
    if (dimensiona.eq.3)then
      call outwritegrid(n)
  
    eLSE
      call outwritegrid2D(n)

    END IF
  
END IF

IF (TECPLOT.EQ.2)THEN		!BINARY PARAVIEW 3D ONLY

  if (dimensiona.eq.3)then

        call OUTWRITEPARA3Db
    else
    
        call OUTWRITEPARA2Db
    end if

END IF


IF (TECPLOT.EQ.3)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEPARA3DbP

END IF


IF (TECPLOT.EQ.4)THEN		!BINARY tecplot partitioned 3D ONLY

  call OUTWRITEtec3DbP
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
END IF

IF (TECPLOT.EQ.5)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS USING MPI-IO

CALL PARALLEL_VTK_COMBINE(N)

END IF

IF (TECPLOT.EQ.6)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS



CALL PARALLEL_VTK_COMBINE_PARTITIONED(N)

END IF






END SUBROUTINE GRID_WRITE


SUBROUTINE SURF_WRITE
!> @brief
!> This subroutine calls the appropriate surface writing subroutine based on the settings
IMPLICIT NONE
IF (TECPLOT.EQ.1)THEN
   if (dimensiona.eq.3)then
call OUTWRITEGRIDBs
ELSE
call OUTWRITEGRIDBs2d
END IF
END IF
IF (TECPLOT.EQ.0)THEN
 if (dimensiona.eq.3)then
call OUTWRITEGRIDs
ELSE
call OUTWRITEGRIDs2d
END IF

END IF


IF (TECPLOT.eq.3)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEPARA3Dsb

END IF


IF (TECPLOT.EQ.5)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS



CALL PARALLEL_VTK_COMBINE_WALL(N)

END IF


IF (TECPLOT.EQ.6)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS



CALL PARALLEL_VTK_COMBINE_PARTITIONED_WALL(N)

END IF




END SUBROUTINE SURF_WRITE


SUBROUTINE VOLUME_SOLUTION_WRITE
!> @brief
!> This subroutine calls the appropriate volume writing subroutine based on the settings
IMPLICIT NONE
				  

				  IF (N.EQ.0)THEN
				  OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
				  WRITE(63,*)"output1",T
				  CLOSE(63)
				  END IF
	  IF (TECPLOT.EQ.1)THEN
				if (dimensiona.eq.3)then
		
					IF (N.EQ.0)THEN
					  OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
					  WRITE(63,*)"output2",T
					  CLOSE(63)
					  END IF
					if (fastmovie.eq.1)then

					call movie

					else
					call OUTWRITE3vb

					end if
					
					
						IF (N.EQ.0)THEN
					  OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
					  WRITE(63,*)"output3",T
					  CLOSE(63)
					  END IF
				eLSE
				      IF (N.EQ.0)THEN
					  OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
					  WRITE(63,*)"output1",T
					  CLOSE(63)
					  END IF


				call OUTWRITE3vb2D
				
					IF (N.EQ.0)THEN
					  OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
					  WRITE(63,*)"output3",T
					  CLOSE(63)
					  END IF
				END IF
	END IF
	IF (TECPLOT.EQ.0)THEN
  			if (dimensiona.eq.3)then

				call OUTWRITE3v
		      eLSE

 				call OUTWRITE3v2D
		      END IF

	END IF
	
	IF (TECPLOT.EQ.2)THEN		!BINARY PARAVIEW 3D ONLY
	if (dimensiona.eq.3)then




        if (fastmovie.eq.1)then

					call movie_PARA

					else

					call OUTWRITEPARA3Db

					END IF
    else
    
        call OUTWRITEPARA2Db
    end if

END IF

IF (TECPLOT.EQ.3)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEPARA3DbP

END IF

IF (TECPLOT.EQ.4)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEtec3DbP
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
end if

IF (TECPLOT.EQ.5)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS USING MPI-IO

CALL PARALLEL_VTK_COMBINE(N)

END IF

IF (TECPLOT.EQ.6)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS

CALL PARALLEL_VTK_COMBINE_PARTITIONED(N)

END IF




					IF (N.EQ.0)THEN
				  OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
				  WRITE(63,*)"finished writing output",T
				  CLOSE(63)
				  END IF





END SUBROUTINE VOLUME_SOLUTION_WRITE



SUBROUTINE surface_SOLUTION_WRITE
!> @brief
!> This subroutine calls the appropriate surface solution writing subroutine based on the settings
IMPLICIT NONE

IF (TECPLOT.EQ.1)THEN
    if (dimensiona.eq.3)then

  call OUTWRITE3vsb
  eLSE

  call OUTWRITE3vsb2D
  END IF
END IF
IF (TECPLOT.EQ.0)THEN
      if (dimensiona.eq.3)then

    call OUTWRITE3vs
    eLSE

    call OUTWRITE3vs2D
    END IF

END IF
IF (TECPLOT.eq.2)THEN		!BINARY PARAVIEW 3D ONLY
    if (dimensiona.eq.3)then

  call OUTWRITEPARA3Dsb
    else
  !call OUTWRITEPARA2Dsb    !not implemented yet
    end if

END IF


IF (TECPLOT.eq.3)THEN		!BINARY PARAVIEW 3D ONLY
    if (dimensiona.eq.3)then

  call OUTWRITEPARA3Dsb
    else
  !call OUTWRITEPARA2Dsb    !not implemented yet
    end if

END IF


IF (TECPLOT.EQ.4)THEN
    if (dimensiona.eq.3)then

  call OUTWRITE3vsb
  eLSE

  call OUTWRITE3vsb2D
  END IF
END IF

IF (TECPLOT.eq.5)THEN		!BINARY PARAVIEW 3D ONLY

  call PARALLEL_VTK_COMBINE_WALL(N)

END IF


IF (TECPLOT.EQ.6)THEN		!FAST OUTPUT WRITTEN BY ALL PROCESSORS



CALL PARALLEL_VTK_COMBINE_PARTITIONED_WALL(N)

END IF


END SUBROUTINE surface_SOLUTION_WRITE

SUBROUTINE VOLUME_SOLUTION_WRITE_av
!> @brief
!> This subroutine calls the appropriate average volume writing subroutine based on the settings
IMPLICIT NONE

IF (TECPLOT.EQ.1)THEN


call OUTWRITE3vbav

END IF

IF (TECPLOT.EQ.0)THEN
    if (dimensiona.eq.3)then

  call OUTWRITE3vav
  eLSE

  call OUTWRITE3v2Dav
  END IF

END IF


IF (TECPLOT.EQ.2)THEN

CALL OUTWRITEPARA3Dbav

END IF


IF (TECPLOT.EQ.3)THEN

CALL OUTWRITEPARA3DbPav

END IF


IF (TECPLOT.EQ.4)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEtec3DbPav
  
end if

IF (TECPLOT.EQ.5)THEN
	IF (DIMENSIONA.EQ.3)THEN
	CALL PARALLEL_VTK_COMBINE_AV(N)
	END IF
END IF

IF (TECPLOT.EQ.6)THEN
	IF (DIMENSIONA.EQ.3)THEN
	CALL PARALLEL_VTK_COMBINE_PARTITIONED_AV(N)
	END IF
END IF




END SUBROUTINE VOLUME_SOLUTION_WRITE_av



SUBROUTINE surface_SOLUTION_WRITE_av
!> @brief
!> This subroutine calls the appropriate surface writing subroutine based on the settings
IMPLICIT NONE

IF (TECPLOT.EQ.1)THEN
   if (dimensiona.eq.3)then

call OUTWRITE3vsbav
eLSE

call OUTWRITE3vsb2Dav
END IF
END IF

IF (TECPLOT.EQ.0)THEN
  if (dimensiona.eq.3)then

call OUTWRITE3vsav
eLSE

call OUTWRITE3vs2Dav
END IF

END IF
IF (TECPLOT.eq.2)THEN

CALL OUTWRITEPARA3Dsbav

END IF

IF (TECPLOT.eq.3)THEN

CALL OUTWRITEPARA3Dsbav

END IF



IF (TECPLOT.EQ.4)THEN
   if (dimensiona.eq.3)then

call OUTWRITE3vsbav
eLSE

call OUTWRITE3vsb2Dav
END IF
END IF

IF (TECPLOT.EQ.5)THEN
	IF (DIMENSIONA.EQ.3)THEN
	CALL PARALLEL_VTK_COMBINE_WALL_AV(N)
	END IF
END IF



IF (TECPLOT.EQ.6)THEN
	IF (DIMENSIONA.EQ.3)THEN
	CALL PARALLEL_VTK_COMBINE_partitioned_wall_av(N)
	END IF
END IF





END SUBROUTINE surface_SOLUTION_WRITE_av

SUBROUTINE forces
!> @brief
!> This subroutine calls the appropriate force computation subroutine based on the dimensionality of the problem
IMPLICIT NONE


   if (dimensiona.eq.3)then

call computeforce(n)
eLSE

call computeforce2d(n)
END IF



END SUBROUTINE forces



SUBROUTINE RESIDUAL_COMPUTE
!> @brief
!> This subroutine calls the appropriate residual computation subroutine based on the dimensionality of the problem
IMPLICIT NONE


   if (dimensiona.eq.3)then

call CALCULATE_RESIDUAL(n)
eLSE

call CALCULATE_RESIDUAL2D(n)
END IF

END SUBROUTINE RESIDUAL_COMPUTE




SUBROUTINE CHECKPOINTING
!> @brief
!> This subroutine calls the appropriate checkpointing subroutine based on the dimensionality of the problem
IMPLICIT NONE


   if (dimensiona.eq.3)then

call CHECKPOINT(n)
eLSE

call CHECKPOINT2D(n)
END IF

END SUBROUTINE CHECKPOINTING


SUBROUTINE CHECKPOINTING_av
!> @brief
!> This subroutine calls the appropriate averaged checkpointing subroutine based on the dimensionality of the problem
IMPLICIT NONE


   if (dimensiona.eq.3)then

call CHECKPOINTav(n)
eLSE

call CHECKPOINTav2D(n)
END IF

END SUBROUTINE CHECKPOINTING_av


SUBROUTINE CHECKPOINT(N)
!> @brief
!> This subroutine uses MPI-IO for writing the checkpointing files
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::dispt
REAL,ALLOCATABLE,DIMENSION(:)::array2
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,size_of_real,size_of_int,dip,N_END,datatype,ifg
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
LOGICAL::HERE1
REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init
disp_in_file=0
tmp=0
disp_init=0

KMAXE=XMPIELRANK(N)

size_of_int=4
size_of_real=8
ICPUID=N
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
!  ioCPt1=MPI_WTIME()


IF (DG.EQ.1)THEN
ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1)))
ELSE
ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)))	
END IF

     if (dg.eq.1)then
     DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*((NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1))
      END DO
      

      n_end=(NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1)
     
     else

      DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*((NOF_VARIABLES+turbulenceequations+passivescalar))
      END DO
      

      n_end=NOF_VARIABLES+turbulenceequations+passivescalar
      end if

      
      if (dg.eq.1)then
      IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	    K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	      ARRAY2(K:K+turbulenceequations+passivescalar-1)=U_CT(I)%VAL(1,1:turbulenceequations+passivescalar)
	      K=K+turbulenceequations+passivescalar
	  END DO
      ELSE
	  K=1
	  DO I=1,KMAXE
        DO J=1,NOF_VARIABLES
	      ARRAY2(K:K+IDEGFREE)=U_C(I)%VALDG(1,J,1:IDEGFREE+1)
	      K=K+(IDEGFREE+1)
        END DO
	  END DO
      END IF
      
      else
      
      IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	    K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	      ARRAY2(K:K+turbulenceequations+passivescalar-1)=U_CT(I)%VAL(1,1:turbulenceequations+passivescalar)
	      K=K+turbulenceequations+passivescalar
	  END DO
      ELSE
	  K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	  END DO
      END IF
      end if
      
	RESTFILE='RESTART.dat'

       IF (N.EQ.0)THEN
	 INQUIRE (FILE=RESTFILE,EXIST=HERE1)
	  CALL MPI_FILE_DELETE(RESTFILE,MPI_INFO_NULL,IERROR)
	END IF
     
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    !CREATE TYPE FIRST OF INDEXED BLOCK
    
    CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXE,n_end,DISPT,MPI_DOUBLE_PRECISION,DATATYPE,IERROR)
    CALL MPI_TYPE_COMMIT(DATATYPE,IERROR)
   

    ALLOCATE(ARRAY(1:nof_Variables+turbulenceequations+passivescalar))
    
    
	
   
	
	
	
	!CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

	call MPI_file_open(MPI_COMM_WORLD, RESTFILE,MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, fh, ierror)
	
	
	 
	  
	
	if (n.eq.0)then
	
	    if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	          call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_write(fh, it, 1, MPI_INTEGER, MPI_STATUS_IGNORE,ierror)
		    disp_in_file = disp_in_file + size_of_int	!1
		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_write(fh, INITIALRES(1:nof_Variables+turbulenceequations),nof_Variables+turbulenceequations, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_Variables+turbulenceequations)	!3
	    ELSE
		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  
		  call MPI_file_write(fh, it, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
		    disp_in_file = disp_in_file + size_of_int 	!4
		    
		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  
		  call MPI_file_write(fh, T, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
		    disp_in_file = disp_in_file + size_of_real    !5
		      if (initcond.eq.95)then
		      call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
		      call MPI_file_write(fh, TAYLOR, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
			    disp_in_file = disp_in_file + size_of_real !6
		      end if  
	    end if
	ELSE
	      if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
		disp_in_file = disp_in_file + size_of_int 
		disp_in_file = disp_in_file + size_of_real*(nof_Variables+turbulenceequations)
	      ELSE
		  disp_in_file = disp_in_file + size_of_int 
		  disp_in_file = disp_in_file + size_of_real
		  if (initcond.eq.95)then
		  disp_in_file = disp_in_file + size_of_real 
		  END IF
	      END IF
	      
	      
	END IF
	
	
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,datatype, 'native',MPI_INFO_NULL, ierror)
	call MPI_FILE_WRITE_ALL(fh, ARRAY2, KMAXE*n_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)        
        
        call MPI_FILE_CLOSE(fh, ierror)
	CALL MPI_TYPE_FREE(DATATYPE,IERROR)
          

	
	
	
	
	DEALLOCATE(ARRAY,DISPT,ARRAY2)
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	
	

	
	


END SUBROUTINE CHECKPOINT



SUBROUTINE PREPARE_SURFACES_V(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,J,K,KMAXE,temp_cord,ILOOP,KLOOP,zloop,kloopf,temp_loop,KN
	INTEGER,ALLOCATABLE,DIMENSION(:)::LIST_NIN2,LIST_NOUT2
	KMAXE=XMPIELRANK(N)
ILOOP=0
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	      if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
		  iloop=iloop+1
	      END IF
	  end if
	END DO
   end if
END DO


iloopx=iloop

if (iloopx.gt.0)then
ALLOCATE(WALL_L(1:ILOOPX,1:4))

eLSE
ALLOCATE(WALL_L(0,0))
end if

ILOOP=0
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	     if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
		  iloop=iloop+1
				if (DIMENSIONA.EQ.2)then
					KLOOP=2
				ELSE
					if (ielem(n,i)%types_faces(J).eq.5)then
					KLOOP=4
					ELSE
					KLOOP=3
					END IF
				END IF

		  WALL_L(ILOOP,1)=I		!this contains the iloop of all the elements that are needed with local numbering
		  WALL_L(ILOOP,2)=j		!face
		  WALL_L(ILOOP,3)=KLOOP !NUMBER OF NODES
	      END IF
	  end if
	END DO
   end if
END DO


	IWMAXE=0


	!NOW FIND ALL THE ELEMENT TYPES AND COPY THEM IN A LOCAL LIST

	CALL MPI_ALLREDUCE(ILOOPX,IWMAXE,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)




	ALLOCATE(WALLSHAPE_G(1:IWMAXE));WALLSHAPE_G=0

	IF (ILOOPX.GT.0)then
	ALLOCATE(WALLSHAPE(1:ILOOPX));WALLSHAPE=0

	DO I=1,ILOOPX
		IF (WALL_L(I,3).EQ.2)THEN	!LINE
			WALLSHAPE(I)=3
		END IF

		IF (WALL_L(I,3).EQ.3)THEN	!TRIANGULAR
			WALLSHAPE(I)=5
		END IF

		IF (WALL_L(I,3).EQ.4)THEN	!QUADRILATERAL
			WALLSHAPE(I)=9
		END IF
	END DO


	eLSE
	ALLOCATE(WALLSHAPE(0:0));WALLSHAPE=0
	END IF



	!WE NEED THE TOTAL WITH AN OFFSET
	ALLOCATE(WALLCX_G(0:ISIZE-1),OFFSETWC_G(0:ISIZE-1))
WALLCX_G(:)=0
OFFSETWC_G(0:ISIZE-1)=0

WALLCX_G(N)=ILOOPX

CALL MPI_ALLGATHER(ILOOPX,1,MPI_INTEGER,WALLCX_G,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)

OFFSETWC_G(0)=0
DO I=1,ISIZE-1
		OFFSETWC_G(I)=OFFSETWC_G(I-1)+WALLCX_g(I-1)
END DO



IF (ILOOPX.GT.0)THEN
WALLSHAPE_G(OFFSETWC_G(N)+1:OFFSETWC_G(N)+WALLCX_g(N))=WALLSHAPE(1:ILOOPX)
END IF










ALLOCATE(WALLSHAPE_G2(1:IWMAXE));WALLSHAPE_G2=0




CALL MPI_ALLREDUCE(WALLSHAPE_G,WALLSHAPE_G2,IWMAXE,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)








	IF (ILOOPX.GT.0)then
	ALLOCATE(TYP_NODESN_W(1:ILOOPX));TYP_NODESN_w(:)=0
	eLSE
	ALLOCATE(TYP_NODESN_W(0:0));TYP_NODESN_W(:)=0
	END IF
	TYP_COUNTN_W=0
	IF (ILOOPX.GT.0)THEN
	temp_loop=OFFSETWC_G(N)
	DO I=1,ILOOPX
		TYP_COUNTN_W=TYP_COUNTN_w+WALL_L(I,3)
		temp_loop=temp_loop+1
		WALL_L(i,4)=temp_loop

		TYP_NODESN_W(I)=WALL_L(I,3)
	END DO
	END IF


allocate(nodes_offsetW(1:IWMAXE),nodes_offsetW2(1:IWMAXE))






CALL MPI_ALLREDUCE(TYP_COUNTN_W,TYP_COUNTN_GLOBAL_W,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)

KN=0
DO J=1,IWMAXE
	nodes_offsetW(j)=kn
		IF (WALLSHAPE_G2(J).EQ.3)THEN
			KN=KN+2
		END IF
		IF (WALLSHAPE_G2(J).EQ.5)THEN
			KN=KN+3
		END IF
		IF (WALLSHAPE_G2(J).EQ.9)THEN
			KN=KN+4
		END IF
	nodes_offsetW2(j)=kn


END DO


	IF (ILOOPX.GT.0)THEN
	allocate(nodes_offset_localW(1:ILOOPX));nodes_offset_localW=0
	allocate(nodes_offset_local2W(1:ILOOPX));nodes_offset_local2W=0




	ELSE
	Allocate(nodes_offset_localW(0:0));nodes_offset_localW=0
	allocate(nodes_offset_local2W(0:0));nodes_offset_local2W=0
	END IF

	IF (ILOOPX.GT.0)THEN
	DO I=1,ILOOPX
	NODES_OFFSET_LOCALW(I)=NODES_OFFSETW(WALL_L(i,4))
	NODES_OFFSET_LOCAL2W(I)=NODES_OFFSETW2(WALL_L(i,4))

	END DO
	END IF




END SUBROUTINE PREPARE_SURFACES_V




















SUBROUTINE PARTITION_PREPARATION_WALLv(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,J,K,temp_cord,KLOOP,varg_max,KMAXN_P
	INTEGER,DIMENSION(4)::TEMPWNODE
	KMAXN_P=XMPIALL_v(N)
	temp_cord=3

	IF (DIMENSIONA.EQ.2)THEN
	WNODES_PART=2
	ELSE
	WNODES_PART=3
	END IF
	varg_max=max(WRITE_VARIABLES_W,WRITE_VARIABLES_Av_w)

	IF (ILOOPX.GT.0)THEN
 	ALLOCATE(WDISPART1(1:ILOOPX),WrARRAY_PART1(1:ILOOPX,1:varg_max))	!

 		DO I=1,ILOOPX
 			WDISPART1(I)=(wall_l(i,4)-1)*1
 		END DO
			WPART1_end=1
	eLSE
		ALLOCATE(WDISPART1(1),WrARRAY_PART1(1,1:WRITE_VARIABLES_W))	!
		WrARRAY_PART1(1,:)=0
		WDISPART1(1)=0!
			WPART1_end=0

	END IF





	IF (ILOOPX.GT.0)THEN
 	ALLOCATE(WDISPART2(1:ILOOPX),WiARRAY_PART2(1:TYP_COUNTN_W))		!
 		DO I=1,iloopx
 			WDISPART2(I)=nodes_offset_localW(i)
 		END DO

 		WPART2_end=WNODES_PART


	eLSE
		ALLOCATE(WDISPART2(1),WiARRAY_PART2(1))	!
		WiARRAY_PART2(1)=0
		WDISPART2(1)=0!(totwallsc-1)*WNODES_PART	!MAYBE TOTAL NUMBER OF VALUES WRITTEN
		WPART2_end=0

	END IF



	!the nodes of the boundary cells

	IF (ILOOPX.GT.0)THEN




 		K=1
 	  DO I=1,iloopx

		WiARRAY_PART2(K:K+TYP_NODESN_W(I)-1)=IELEM(N,WALL_L(I,1))%NODES_FACES_V(WALL_L(I,2),1:TYP_NODESN_W(I))
		K=K+TYP_NODESN_W(I)
 	  END DO


    END IF




		IF (ILOOPX.GT.0)THEN
		ALLOCATE(WDISPART5(1:ILOOPX),WiARRAY_PART5(1:ILOOPX))


			DO I=1,ILOOPX
				WDISPART5(I)=(wall_l(i,4)-1)*1
			END DO

			DO I=1,iloopx
				WiARRAY_PART5(i)=nodes_offset_local2W(i)
			END DO


		ELSE
			ALLOCATE(WDISPART5(1),WiARRAY_PART5(1))
				WDISPART5(1)=0!totwallsc-1
				WiARRAY_PART5(1)=0


		END IF



	IF (ILOOPX.GT.0)THEN

!	!tHE TYPES of elements OF THE nodes OF  tHE boundary cells
	ALLOCATE(wDISPART3(1:ILOOPX),wiARRAY_PART3(1:ILOOPX))

			DO I=1,ILOOPX
				WDISPART3(I)=(wall_l(i,4)-1)*1
				wiARRAY_PART3(I)=WALLSHAPE(I)
			END DO
			wPART3_end=1
	ELSE
	ALLOCATE(wDISPART3(1),wiARRAY_PART3(1))
			WDISPART3(1)=0!totwallsc-1
			WiARRAY_PART3(1)=0
			wPART3_end=0
	END IF





  	ALLOCATE(WDISPART4(1:KMAXN_P),WrARRAY_PART4(1:KMAXN_P*temp_cord))		!
  	wPART4_end=temp_cord

  	DO I=1,KMAXN_P
			WDISPART4(I)=(my_nodesg(I)-1)*(temp_cord)
		END DO

		K=1
	  DO I=1,KMAXN_P
		WrARRAY_PART4(K:K+dims-1)=INODER4(my_nodesl(i))%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		WrARRAY_PART4(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO



!
!
IF (N.EQ.0)THEN
WKDUM1=1
WKDUM2=1
WKDUM3(1)=0
ELSE
WKDUM1=0
WKDUM2=0
WKDUM3(1)=1
END IF
!
			!now commit datatypes

 				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(ILOOPX,WPART1_end,WDISPART1,MPI_DOUBLE_PRECISION,WDATATYPEX,IERROR)
 				CALL MPI_TYPE_COMMIT(WDATATYPEX,IERROR)
!
 				!DUMMY TYPE FOR WRITING ONE COMPONENT ONLY FROM ONE CPU
 				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(WKDUM1,WKDUM2,WKDUM3,MPI_INTEGER,WDATATYPEINT,IERROR)
 				CALL MPI_TYPE_COMMIT(WDATATYPEINT,IERROR)
!
!
!
! 				!point coordinates
 				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXN_P,WPART4_end,WDISPART4,MPI_DOUBLE_PRECISION,WDATATYPEz,IERROR)
 				CALL MPI_TYPE_COMMIT(WDATATYPEz,IERROR)
!
! 				!connectivity
!  				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(ILOOPX,TYP_NODESN_W(:),WDISPART2,MPI_INTEGER,WDATATYPEy,IERROR)
 				CALL MPI_TYPE_INDEXED(ILOOPX,TYP_NODESN_W(:),WDISPART2,MPI_INTEGER,WDATATYPEy,IERROR)
 				CALL MPI_TYPE_COMMIT(WDATATYPEy,IERROR)
!
! 				!type of element
 				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(ILOOPX,WPART1_end,WDISPART5,MPI_INTEGER,WDATATYPEXx,IERROR)
 				CALL MPI_TYPE_COMMIT(WDATATYPEXx,IERROR)
!
! 				!nodes
 				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(ILOOPX,wPART3_end,WDISPART3,MPI_INTEGER,WDATATYPEyy,IERROR)
 				CALL MPI_TYPE_COMMIT(WDATATYPEyy,IERROR)
!
!
!


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


DEALLOCATE(WALLSHAPE_G,WALLSHAPE,WALLCX_G,OFFSETWC_G,WALLSHAPE_G2,TYP_NODESN_W,nodes_offsetW,nodes_offsetW2,NODES_OFFSET_LOCALW,NODES_OFFSET_LOCAL2W)





END SUBROUTINE PARTITION_PREPARATION_WALLv





SUBROUTINE PARTITION_PREPARATION(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,J,K,KMAXE,KMAXN_P,temp_cord,varg_max
	KMAXE=XMPIELRANK(N)
	KMAXN_P=XMPIALL_v(N)

	temp_cord=3

	IF (DIMENSIONA.EQ.2)THEN
	NODES_PART=4
	ELSE
	NODES_PART=8
	END IF
	varg_max=max(WRITE_VARIABLES,WRITE_VARIABLES_Av)
	ALLOCATE(DISPART1(1:KMAXE),rARRAY_PART1(1:KMAXE,1:varg_max))	!

		DO I=1,KMAXE
			DISPART1(I)=(XGO(I)-1)*(1)
		END DO
		PART1_end=1


	ALLOCATE(TYP_NODESN(1:kmaxe));TYP_NODESN(:)=0
	TYP_COUNTN=0
	DO I=1,KMAXE
			IF (IELEM(N,I)%ishape.EQ.1)THEN	!HEXA
				TYP_COUNTN=TYP_COUNTN+8
				TYP_NODESN(I)=8

			END IF
			IF (IELEM(N,I)%ishape.EQ.2)THEN	!TETRA
				TYP_COUNTN=TYP_COUNTN+4
				TYP_NODESN(I)=4
! 				TYP_COUNTN=TYP_COUNTN+8
! 				TYP_NODESN(I)=8

			END IF
			IF (IELEM(N,I)%ishape.EQ.3)THEN	!PYRAMID
				TYP_COUNTN=TYP_COUNTN+5
				TYP_NODESN(I)=5
! 				TYP_COUNTN=TYP_COUNTN+8
! 				TYP_NODESN(I)=8
			END IF
			IF (IELEM(N,I)%ishape.EQ.4)THEN	!PRISM
				TYP_COUNTN=TYP_COUNTN+6
				TYP_NODESN(I)=6
! 				TYP_COUNTN=TYP_COUNTN+8
! 				TYP_NODESN(I)=8
			END IF
			IF (IELEM(N,I)%ishape.EQ.5)THEN	!QUAD
				TYP_COUNTN=TYP_COUNTN+4
				TYP_NODESN(I)=4
			END IF
			IF (IELEM(N,I)%ishape.EQ.6)THEN	!TRIANGULAR
				TYP_COUNTN=TYP_COUNTN+3
				TYP_NODESN(I)=3
! 				TYP_COUNTN=TYP_COUNTN+4
! 				TYP_NODESN(I)=4
			END IF
	END DO




	CALL MPI_ALLREDUCE(TYP_COUNTN,TYP_COUNTN_GLOBAL,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)


! 		print*,TYP_COUNTN,TYP_COUNTN_GLOBAL,"check here",IMAXE*6




	ALLOCATE(DISPART2(1:KMAXE),iARRAY_PART2(1:TYP_COUNTN))		!
		DO I=1,KMAXE
! 			DISPART2(I)=(XGO(I)-1)*(NODES_PART)
!  			DISPART2(I)=(XGO(I)-1)*(TYP_NODESN(I))
 			DISPART2(I)=nodes_offset_local(i)

		END DO

		PART2_end=NODES_PART

			K=1
	  DO I=1,KMAXE
! 		iARRAY_PART2(K:K+NODES_PART-1)=IELEM(N,I)%NODES_v(1:NODES_PART)
! 	      K=K+NODES_PART
	      iARRAY_PART2(K:K+TYP_NODESN(I)-1)=IELEM(N,I)%NODES_v(1:TYP_NODESN(I))
	      K=K+TYP_NODESN(I)
	  END DO



	  ALLOCATE(DISPART5(1:KMAXE),iARRAY_PART5(1:KMAXE))
		DO I=1,KMAXE
			DISPART5(I)=(XGO(I)-1)*(1)
		END DO


	   DO I=1,KMAXE
		iARRAY_PART5(i)=nodes_offset_local2(i)!(XGO(I))*TYP_NODESN(I)
	  END DO




	ALLOCATE(DISPART3(1:KMAXe),iARRAY_PART3(1:KMAXE))
 		DO I=1,KMAXE
 			DISPART3(I)=(XGO(I)-1)*(1)


 			IF (IELEM(N,I)%ishape.EQ.1)THEN	!HEXA
				IARRAY_PART3(I)=12
			END IF
			IF (IELEM(N,I)%ishape.EQ.2)THEN	!TETRA
				IARRAY_PART3(I)=10
			END IF
			IF (IELEM(N,I)%ishape.EQ.3)THEN	!PYRAMID
				IARRAY_PART3(I)=14
			END IF
			IF (IELEM(N,I)%ishape.EQ.4)THEN	!PRISM
				IARRAY_PART3(I)=13
			END IF
			IF (IELEM(N,I)%ishape.EQ.5)THEN	!QUAD
				IARRAY_PART3(I)=9
			END IF
			IF (IELEM(N,I)%ishape.EQ.6)THEN	!TRIANGULAR
				IARRAY_PART3(I)=5
			END IF





 		END DO

		PART3_end=1





	ALLOCATE(DISPART4(1:KMAXN_P),rARRAY_PART4(1:KMAXN_P*temp_cord))

		PART4_end=temp_cord
		DO I=1,KMAXN_P
			DISPART4(I)=(my_nodesg(I)-1)*(temp_cord)
		END DO


		K=1
	  DO I=1,KMAXN_P
		rARRAY_PART4(K:K+dims-1)=INODER4(my_nodesl(i))%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		rARRAY_PART4(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO


IF (N.EQ.0)THEN
KDUM1=1
KDUM2=1
KDUM3(1)=0
ELSE
KDUM1=0
KDUM2=0
KDUM3(1)=1
END IF

 !now commit datatypes

				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXE,PART1_end,DISPART1,MPI_DOUBLE_PRECISION,DATATYPEX,IERROR)
				CALL MPI_TYPE_COMMIT(DATATYPEX,IERROR)

				!DUMMY TYPE FOR WRITING ONE COMPONENT ONLY FROM ONE CPU
				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KDUM1,KDUM2,KDUM3,MPI_INTEGER,DATATYPEINT,IERROR)
				CALL MPI_TYPE_COMMIT(DATATYPEINT,IERROR)



				!point coordinates
				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXN_P,PART4_end,DISPART4,MPI_DOUBLE_PRECISION,DATATYPEz,IERROR)
				CALL MPI_TYPE_COMMIT(DATATYPEz,IERROR)

				!connectivity
! 				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(kmaxe,nodes_part,DISPART2,MPI_INTEGER,DATATYPEy,IERROR)
				CALL MPI_TYPE_INDEXED(kmaxe,TYP_NODESN(1:kmaxe),DISPART2,MPI_INTEGER,DATATYPEy,IERROR)
				CALL MPI_TYPE_COMMIT(DATATYPEy,IERROR)

				!type of element
				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(kmaxe,PART1_end,DISPART5,MPI_INTEGER,DATATYPEXx,IERROR)
				CALL MPI_TYPE_COMMIT(DATATYPEXx,IERROR)

				!nodes
				CALL MPI_TYPE_CREATE_INDEXED_BLOCK(kmaxe,PART1_end,DISPART3,MPI_INTEGER,DATATYPEyy,IERROR)
				CALL MPI_TYPE_COMMIT(DATATYPEyy,IERROR)














END SUBROUTINE PARTITION_PREPARATION

SUBROUTINE PARTITION_PREPARATION_p(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,J,K,KMAXE,KMAXN_P,temp_cord,varg_max
	INTEGER,DIMENSION(8)::TEMP_PART_N
	KMAXE=XMPIELRANK(N)


	temp_cord=3


	varg_max=max(WRITE_VARIABLES,WRITE_VARIABLES_Av)
	ALLOCATE(sol_vtu(1:KMAXE,1:varg_max));sol_vtu=0	!

	ALLOCATE(TYP_NODESN(1:kmaxe));TYP_NODESN(:)=0
	TYP_COUNTN=0
	DO I=1,KMAXE
			IF (IELEM(N,I)%ishape.EQ.1)THEN	!HEXA
				TYP_COUNTN=TYP_COUNTN+8
				TYP_NODESN(I)=8

			END IF
			IF (IELEM(N,I)%ishape.EQ.2)THEN	!TETRA
				TYP_COUNTN=TYP_COUNTN+4
				TYP_NODESN(I)=4


			END IF
			IF (IELEM(N,I)%ishape.EQ.3)THEN	!PYRAMID
				TYP_COUNTN=TYP_COUNTN+5
				TYP_NODESN(I)=5

			END IF
			IF (IELEM(N,I)%ishape.EQ.4)THEN	!PRISM
				TYP_COUNTN=TYP_COUNTN+6
				TYP_NODESN(I)=6

			END IF
			IF (IELEM(N,I)%ishape.EQ.5)THEN	!QUAD
				TYP_COUNTN=TYP_COUNTN+4
				TYP_NODESN(I)=4
			END IF
			IF (IELEM(N,I)%ishape.EQ.6)THEN	!TRIANGULAR
				TYP_COUNTN=TYP_COUNTN+3
				TYP_NODESN(I)=3

			END IF
	END DO






	ALLOCATE(offset_vtu(1:KMAXE),connect_vtu(1:TYP_COUNTN))
	!
	K=0
		DO I=1,KMAXE
				K=k+TYP_NODESN(I)
			offset_vtu(I)=k
		END DO



			K=1
	  DO I=1,KMAXE
			TEMP_PART_N(1:TYP_NODESN(I))=IELEM(N,I)%NODES(1:TYP_NODESN(I))

			DO J=1,TYP_NODESN(I)
				TEMP_PART_N(J)=TEMP_PART_N(J)-1
			END DO

	      connect_vtu(K:K+TYP_NODESN(I)-1)=TEMP_PART_N(1:TYP_NODESN(I))
	      K=K+TYP_NODESN(I)
	  END DO





	ALLOCATE(type_vtu(1:KMAXE))
 		DO I=1,KMAXE
 			IF (IELEM(N,I)%ishape.EQ.1)THEN	!HEXA
				type_vtu(I)=12
			END IF
			IF (IELEM(N,I)%ishape.EQ.2)THEN	!TETRA
				type_vtu(I)=10
			END IF
			IF (IELEM(N,I)%ishape.EQ.3)THEN	!PYRAMID
				type_vtu(I)=14
			END IF
			IF (IELEM(N,I)%ishape.EQ.4)THEN	!PRISM
				type_vtu(I)=13
			END IF
			IF (IELEM(N,I)%ishape.EQ.5)THEN	!QUAD
				type_vtu(I)=9
			END IF
			IF (IELEM(N,I)%ishape.EQ.6)THEN	!TRIANGULAR
				type_vtu(I)=5
			END IF


 		END DO






	ALLOCATE(nodes_vtu(1:KMAXN*temp_cord))





		K=1
	  DO I=1,KMAXN
		nodes_vtu(K:K+dims-1)=INODER4(I)%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		nodes_vtu(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO




END SUBROUTINE PARTITION_PREPARATION_p


SUBROUTINE PARTITION_PREPARATION_p_WALL(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,J,K,KMAXE,KMAXN_P,temp_cord,varg_max,ILOOP,KLOOP,zloop,kloopf,temp_loop,KN,KKD
	INTEGER,DIMENSION(8)::TEMP_PART_N
	KMAXE=XMPIELRANK(N)


	temp_cord=3
				ILOOP=0
			DO I=1,KMAXE
			if (ielem(n,i)%interior.eq.1)then
				DO j=1,IELEM(N,I)%IFCA
				if (ielem(n,i)%ibounds(J).gt.0)then
					if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
					iloop=iloop+1
					END IF
				end if
				END DO
			end if
			END DO



			iloopx=iloop


	allocate(wallcount_cpu_l(0:isize-1),wallcount_cpu_g(0:isize-1));wallcount_cpu_l=0;wallcount_cpu_g=0


	wallcount_cpu_l(n)=iloopx

	call MPI_ALLREDUCE(wallcount_cpu_l(0:isize-1),wallcount_cpu_g(0:isize-1),isize,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)


	call mpi_barrier(mpi_comm_world,ierror)


if (iloopx.gt.0)then
ALLOCATE(WALL_L(1:ILOOPX,1:4))

eLSE
ALLOCATE(WALL_L(0,0))
end if

ILOOP=0
DO I=1,KMAXE
  if (ielem(n,i)%interior.eq.1)then
	DO j=1,IELEM(N,I)%IFCA
	  if (ielem(n,i)%ibounds(J).gt.0)then
	    if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).or.(ibound(n,ielem(n,i)%ibounds(j))%icode.eq.99))then
		  iloop=iloop+1
				if (DIMENSIONA.EQ.2)then
					KLOOP=2
				ELSE
					if (ielem(n,i)%types_faces(J).eq.5)then
					KLOOP=4
					ELSE
					KLOOP=3
					END IF
				END IF

		  WALL_L(ILOOP,1)=I		!this contains the iloop of all the elements that are needed with local numbering
		  WALL_L(ILOOP,2)=j		!face
		  WALL_L(ILOOP,3)=KLOOP !NUMBER OF NODES
	      END IF
	  end if
	END DO
   end if
END DO



	varg_max=max(WRITE_VARIABLES_w,WRITE_VARIABLES_Av_w)
	IF (ILOOPX.GT.0)THEN
	ALLOCATE(sol_vtu_W(1:ILOOPX,1:varg_max));sol_vtu=0	!
	ELSE
	ALLOCATE(sol_vtu_W(0:0,1:varg_max));sol_vtu=0	!
	END IF

	IF (ILOOPX.GT.0)THEN
	ALLOCATE(TYP_NODESN_W(1:ILOOPX),type_vtu_W(1:ILOOPX));TYP_NODESN_W(:)=0;type_vtu_w=0
	eLSE
	ALLOCATE(TYP_NODESN_W(0:0),type_vtu_w(0:0));TYP_NODESN_W(:)=0;type_vtu_w=0
	END IF



	TYP_COUNTN_W=0
	IF (ILOOPX.GT.0)THEN



	DO I=1,ILOOPX
			IF (WALL_L(I,3).EQ.4)THEN	!QUAD
				TYP_COUNTN_W=TYP_COUNTN+4
				TYP_NODESN_W(I)=4
				type_vtu_W(I)=9
			END IF

			IF (WALL_L(I,3).EQ.3)THEN	!TRI
				TYP_COUNTN_W=TYP_COUNTN+3
				TYP_NODESN_W(I)=3
				type_vtu_W(I)=5
			END IF

			IF (WALL_L(I,3).EQ.2)THEN	!LINE
				TYP_COUNTN_W=TYP_COUNTN+2
				TYP_NODESN_W(I)=2
				type_vtu_W(I)=3
			END IF
	END DO

	END IF



	IF (ILOOPX.GT.0)THEN
	ALLOCATE(offset_vtu_w(1:ILOOPX),connect_vtu_w(1:TYP_COUNTN_W))
	ELSE
	ALLOCATE(offset_vtu_w(0:0),connect_vtu_w(0:0))
	END IF
	!

	IF (ILOOPX.GT.0)THEN
	K=0
		DO I=1,ILOOPX
				K=k+TYP_NODESN_W(I)
			offset_vtu_W(I)=k
		END DO



			K=1
	  DO I=1,ILOOPX
			TEMP_PART_N(1:TYP_NODESN_W(I))=IELEM(N,WALL_L(I,1))%NODES_FACES(WALL_L(I,2),1:TYP_NODESN_W(I))

			DO J=1,TYP_NODESN_W(I)
				TEMP_PART_N(J)=TEMP_PART_N(J)-1
			END DO

	      connect_vtu_W(K:K+TYP_NODESN_w(I)-1)=TEMP_PART_N(1:TYP_NODESN_W(I))
	      K=K+TYP_NODESN_W(I)
	  END DO


	END IF








	IF (ILOOPX.GT.0)THEN
	ALLOCATE(nodes_vtu_w(1:KMAXN*temp_cord))

	eLSE
	ALLOCATE(nodes_vtu_w(0:0))
	END IF


	IF (ILOOPX.GT.0)THEN

		K=1
	  DO I=1,KMAXN
		nodes_vtu_W(K:K+dims-1)=INODER4(I)%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		nodes_vtu_W(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO
	END IF



END SUBROUTINE PARTITION_PREPARATION_p_WALL



SUBROUTINE SPECIFY_WRITE_VARIABLES(N)
!> @brief
!> This subroutine uses MPI-IO for writing the VTK FILES
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,J,k




!volume instantaneous variables

if (dimensiona.eq.3)then



IF (multispecies.EQ.1)THEN
WRITE_VARIABLES=NOF_VARIABLES+1
!!specify the name of the variable names!!

Variable_names(1)='density'
Variable_names(2)='U'
Variable_names(3)='V'
Variable_names(4)='W'
Variable_names(5)='Pressure'
Variable_names(6)='rho vf1'
Variable_names(7)='rho vf2'
Variable_names(8)='volume_fraction'
Variable_names(9)='Q'


Else

WRITE_VARIABLES=NOF_VARIABLES+1+TURBULENCEEQUATIONS+adda


!!specify the name of the variable names!!

Variable_names(1)='density'
Variable_names(2)='U'
Variable_names(3)='V'
Variable_names(4)='W'
Variable_names(5)='Pressure'
Variable_names(6)='Q'

	if (adda.eq.1)then
	Variable_names(NOF_VARIABLES+1+adda)='ADDA'
	end if

	if (turbulence.eq.1)then
	Variable_names(NOF_VARIABLES+1+TURBULENCEEQUATIONS+adda)='turb'
	end if

	IF (ITESTCASE.EQ.1)THEN
	Variable_names(1)='solution'
	Variable_names(2)='aux'
	end if



END IF


else




IF (multispecies.EQ.1)THEN
WRITE_VARIABLES=NOF_VARIABLES+1
!!specify the name of the variable names!!

Variable_names(1)='density'
Variable_names(2)='U'
Variable_names(3)='V'
Variable_names(4)='Pressure'
Variable_names(5)='rho vf1'
Variable_names(6)='rho vf2'
Variable_names(7)='volume_fraction'
Variable_names(8)='Q'


Else

WRITE_VARIABLES=NOF_VARIABLES+1+TURBULENCEEQUATIONS
!!specify the name of the variable names!!

Variable_names(1)='density'
Variable_names(2)='U'
Variable_names(3)='V'
Variable_names(4)='Pressure'
Variable_names(5)='Q'
if (turbulence.eq.1)then
Variable_names(6)='turb'
end if



IF (ITESTCASE.EQ.1)THEN
Variable_names(1)='solution'
Variable_names(2)='aux'
end if



end if



end if





if(Averaging.eq.1)then

if (dimensiona.eq.3)then




WRITE_VARIABLES_Av=11


!!specify the name of the variable names!!

Variable_names_AV(1)='R_mean'
Variable_names_AV(2)='U_mean'
Variable_names_AV(3)='V_mean'
Variable_names_AV(4)='W_mean'
Variable_names_AV(5)='P_mean'
Variable_names_AV(6)='U_rms'
Variable_names_AV(7)='V_rms'
Variable_names_AV(8)='W_rms'
Variable_names_AV(9)='UV'
Variable_names_AV(10)='UW'
Variable_names_AV(11)='WV'

WRITE_VARIABLES_Av_w=8

Variable_names_AV_w(1)='R_mean'
Variable_names_AV_w(2)='U_mean'
Variable_names_AV_w(3)='V_mean'
Variable_names_AV_w(4)='W_mean'
Variable_names_AV_w(5)='P_mean'
Variable_names_AV_w(6)='ssx_mean'
Variable_names_AV_w(7)='ssy_mean'
Variable_names_AV_w(8)='ssz_mean'
end if
end if




































!wall instantaneous variables



 if (dimensiona.eq.3)then


 if (itestcase.eq.4)then
 WRITE_VARIABLES_W=NOF_VARIABLES+5
 Variable_names_W(1)='density'
Variable_names_W(2)='U'
Variable_names_W(3)='V'
Variable_names_W(4)='W'
Variable_names_W(5)='Pressure'
Variable_names_W(6)='Q'
Variable_names_W(7)='ssx'
Variable_names_W(8)='ssy'
Variable_names_W(9)='ssz'
Variable_names_W(10)='q_heat'


 end if

 if (itestcase.eq.3)then
 WRITE_VARIABLES_W=NOF_VARIABLES+1
Variable_names_W(1)='density'
Variable_names_W(2)='U'
Variable_names_W(3)='V'
Variable_names_W(4)='W'
Variable_names_W(5)='Pressure'
Variable_names_W(6)='Q'
 end if



 if (itestcase.eq.-1)then
 WRITE_VARIABLES_W=NOF_VARIABLES+1
 Variable_names_W(1)='density'
Variable_names_W(2)='U'
Variable_names_W(3)='V'
Variable_names_W(4)='W'
Variable_names_W(5)='Pressure'
Variable_names_W(6)='rho vf1'
Variable_names_W(7)='rho vf2'
Variable_names_W(8)='volume_fraction'
Variable_names_W(9)='Q'
 end if




ELSE

if (itestcase.eq.4)then
 WRITE_VARIABLES_W=NOF_VARIABLES+4
 Variable_names_W(1)='density'
Variable_names_W(2)='U'
Variable_names_W(3)='V'
Variable_names_W(4)='Pressure'
Variable_names_W(5)='Q'
Variable_names_W(6)='ssx'
Variable_names_W(7)='ssy'
Variable_names_W(8)='q_heat'


 end if

 if (itestcase.eq.3)then
 WRITE_VARIABLES_W=NOF_VARIABLES+1
Variable_names_W(1)='density'
Variable_names_W(2)='U'
Variable_names_W(3)='V'
Variable_names_W(4)='Pressure'
Variable_names_W(5)='Q'

 end if

 if (itestcase.eq.-1)then
 WRITE_VARIABLES_W=NOF_VARIABLES+1
 Variable_names_W(1)='density'
Variable_names_W(2)='U'
Variable_names_W(3)='V'
Variable_names_W(4)='Pressure'
Variable_names_W(5)='rho vf1'
Variable_names_W(6)='rho vf2'
Variable_names_W(7)='volume_fraction'
Variable_names_W(8)='Q'
 end if



END IF






END SUBROUTINE SPECIFY_WRITE_VARIABLES


SUBROUTINE PARALLEL_VTK_COMBINE(N)
!> @brief
!> This subroutine uses MPI-IO for writing the VTK FILES
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord
CHARACTER(LEN=20)::PROC,FILEX,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
LOGICAL::HERE1
REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
INTEGER                     :: nbytes,eight
CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
character(LEN=:),allocatable::VTU
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
nbytes=1


 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0
KMAXE=XMPIELRANK(N)
KMAXN_P=XMPIALL_v(N)


temp_cord=3


! 						if (n.eq.0)then
                               WRITE(PROC3,FMT='(I10)') IT
                               FILEX="OUT_"//TRIM(ADJUSTL(PROC3))//".vtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)
!                             end if

		IF (MOVEMENT.EQ.1)THEN

					K=1
	  DO I=1,KMAXN_P
		rARRAY_PART4(K:K+dims-1)=INODER4(my_nodesl(i))%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		rARRAY_PART4(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO


	  END IF












				if (dimensiona.eq.3)then
				DO I=1,KMAXE
				leftv(1:nof_Variables)=U_C(I)%VAL(1,1:NOF_VARIABLES)
				call CONS2PRIM(N,leftv,MP_PINFl,gammal)
					rARRAY_PART1(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)
										do j=nof_Variables+1,write_variables-TURBULENCEEQUATIONS
                                        if (multispecies.eq.1)then
											rARRAY_PART1(i,j)=ielem(n,i)%REDUCE!ielem(n,i)%vortex(1)
                                        else

											rARRAY_PART1(i,j)=ielem(n,i)%vortex(1)
											 if (j.eq.write_variables-TURBULENCEEQUATIONS)then
											 if (adda.eq.1)then
											 rARRAY_PART1(i,j)=ielem(n,i)%diss
											 end if
											 end if
                                        end if        
                                        end do
                                        IF (TURBULENCEEQUATIONS.GT.0)THEN
                                        rARRAY_PART1(i,write_variables)=U_CT(I)%VAL(1,1)
                                        END IF

				END DO

				temp_node=8;temp_dims=3

				end if

				if (dimensiona.eq.2)then
				DO I=1,KMAXE
				leftv(1:nof_Variables)=U_C(I)%VAL(1,1:NOF_VARIABLES)
				call cons2prim(N,leftv,MP_PINFl,gammal)
					rARRAY_PART1(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)
										do j=nof_Variables+1,write_variables-TURBULENCEEQUATIONS
										if (multispecies.eq.1)then

										rARRAY_PART1(i,j)=IELEM(N,I)%REDUCE!ielem(n,i)%vortex(1)
                                        else
                                        if (mood.eq.1)then
                                         rARRAY_PART1(i,j)=ielem(n,i)%mood_o
                                        else
										if (Dg.eq.1)then
											rARRAY_PART1(i,j)=ielem(n,i)%troubled
										else
                                        rARRAY_PART1(i,j)=ielem(n,i)%vortex(1)
                                        end if
                                        end if
                                        end if
										end do
										IF (TURBULENCEEQUATIONS.GT.0)THEN
                                        rARRAY_PART1(i,write_variables)=U_CT(I)%VAL(1,1)
                                        END IF
				END DO
				temp_node=4;temp_dims=3
				end if



temp_imaxe=imaxe
temp_imaxn=imaxn


if (n.eq.0)then

	!first write the header xml file from one MPI process

	lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
	 offset_temp=0
     WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify field pieces
    WRITE(tempstamp1,'(I16)')imaxn
    WRITE(tempstamp2,'(I16)')imaxe
    Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
           &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    DO i=1,WRITE_VARIABLES
      Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names(i))//'" '// &
                       'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      WRITE(Offset_stamp,'(I16)')offset_temp
    END DO
    Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
    ! Specify necessary cell data
    Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
    ! Connectivity
    Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+TYP_COUNTN_GLOBAL*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Offsets
    Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Elem types
    Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
    Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Prepare append section
    Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
    ! Write leading data underscore
    Buffer=TRIM('_');WRITE(300) Buffer
	Bytes = size_of_real
	close(300)


end if


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


				call MPI_file_open(MPI_COMM_WORLD,VTU,MPI_MODE_WRONLY + MPI_MODE_APPEND,MPI_INFO_NULL, fh, ierror)
				call MPI_FILE_GET_POSITION(fh, disp_in_file, ierror)
				disp_init=disp_in_FILE

				!----write time stamp----!
				IF (N.EQ.0)THEN
				call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
				BYTES=size_of_real
				call MPI_file_write(fh, bytes, nbytes, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
				disp_in_file = disp_in_file + size_of_int
				call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
				call MPI_file_write(fh, T, nbytes, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
				disp_in_file=disp_in_file+size_of_Real
				else
				disp_in_file=disp_in_file+size_of_int+size_of_real
				end if
				!end time stamp




				do i=1,WRITE_VARIABLES
				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxe*size_of_real
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int
				!write variables---within loop
				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,DATATYPEX,'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh,rARRAY_PART1(1:kmaxe,i),KMAXE*PART1_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierror)
				!end write variables---within loop
				disp_in_file=disp_in_file+temp_imaxe*size_of_real
				!end loop
				end do

! 				IF (N.EQ.0)print*,"LOCATION2",disp_in_file

				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxn*size_of_real*temp_dims
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int

				call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_DOUBLE_PRECISION,DATATYPEz,'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh,rARRAY_PART4,KMAXN_P*PART4_end,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)


				disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)

! 				IF (N.EQ.0)print*,"LOCATION3",disp_in_file


				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=size_of_int*TYP_COUNTN_GLOBAL!temp_imaxe*size_of_int*temp_node
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int


				call MPI_FILE_SET_VIEW(fh,disp_in_file,MPI_INTEGER,DATATYPEy,'native',MPI_INFO_NULL,ierror)

				call MPI_FILE_WRITE_ALL(fh,iARRAY_PART2,TYP_COUNTN, MPI_INTEGER,STATUS,ierror)

				disp_in_file=disp_in_file+(size_of_int*TYP_COUNTN_GLOBAL)!(temp_imaxe*size_of_int*temp_node)

! 				IF (N.EQ.0)print*,"LOCATION4",disp_in_file

				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxe*size_of_int
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int



				call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,DATATYPEXx, 'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh, iARRAY_PART5,kmaxe*PART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)


! 				IF (N.EQ.0)print*,"LOCATION5",disp_in_file


				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxe*size_of_int
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int


				call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,DATATYPEyy, 'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh, iARRAY_PART3,kmaxe*PART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)


				disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)

				call MPI_FILE_CLOSE(fh, ierror)
				CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)





if (n.eq.0)then
OPEN(300,FILE=FILEX,ACCESS='STREAM',position='APPEND')
  lf = char(10)
  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)
end if


 DEallocate(vTU)


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







END SUBROUTINE PARALLEL_VTK_COMBINE


SUBROUTINE PARALLEL_VTK_COMBINE_AV(N)
!> @brief
!> This subroutine uses MPI-IO for writing the VTK FILES
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord,IND1
CHARACTER(LEN=40)::PROC,FILEX,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
LOGICAL::HERE1
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
INTEGER                     :: nbytes,eight
CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
character(LEN=:),allocatable::VTU
nbytes=1


 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0
KMAXE=XMPIELRANK(N)
KMAXN_P=XMPIALL_v(N)


temp_cord=3


if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if



! 						if (n.eq.0)then
                               WRITE(PROC3,FMT='(I10)') IT
                               FILEX="VOL_AVER"//TRIM(ADJUSTL(PROC3))//".vtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)
!                             end if

		IF (MOVEMENT.EQ.1)THEN

					K=1
	  DO I=1,KMAXN_P
		rARRAY_PART4(K:K+dims-1)=INODER4(my_nodesl(i))%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		rARRAY_PART4(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO


	  END IF












				if (dimensiona.eq.3)then
				DO I=1,KMAXE
				leftv(1:nof_Variables)=U_C(I)%VAL(ind1,1:NOF_VARIABLES)
				call CONS2PRIM(N,leftv,MP_PINFl,gammal)
					rARRAY_PART1(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)
					do j=nof_Variables+1,write_variables_AV
					rARRAY_PART1(i,j)=U_C(I)%RMS(J-nof_Variables)
                     end do
				END DO

				temp_node=8;temp_dims=3

				end if





temp_imaxe=imaxe
temp_imaxn=imaxn


if (n.eq.0)then

	!first write the header xml file from one MPI process

	lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
	 offset_temp=0
     WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify field pieces
    WRITE(tempstamp1,'(I16)')imaxn
    WRITE(tempstamp2,'(I16)')imaxe
    Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
           &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    DO i=1,write_variables_AV
      Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names_AV(i))//'" '// &
                       'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      WRITE(Offset_stamp,'(I16)')offset_temp
    END DO
    Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <Points>'//lf);WRITE(300)Buffer
    Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
    ! Specify necessary cell data
    Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
    ! Connectivity
    Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+TYP_COUNTN_GLOBAL*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Offsets
    Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Elem types
    Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
    Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Prepare append section
    Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
    ! Write leading data underscore
    Buffer=TRIM('_');WRITE(300) Buffer
	Bytes = size_of_real
	close(300)


end if


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


				call MPI_file_open(MPI_COMM_WORLD,VTU,MPI_MODE_WRONLY + MPI_MODE_APPEND,MPI_INFO_NULL, fh, ierror)
				call MPI_FILE_GET_POSITION(fh, disp_in_file, ierror)
				disp_init=disp_in_FILE

				!----write time stamp----!
				IF (N.EQ.0)THEN
				call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
				BYTES=size_of_real
				call MPI_file_write(fh, bytes, nbytes, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
				disp_in_file = disp_in_file + size_of_int
				call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
				call MPI_file_write(fh, T, nbytes, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
				disp_in_file=disp_in_file+size_of_Real
				else
				disp_in_file=disp_in_file+size_of_int+size_of_real
				end if
				!end time stamp




				do i=1,write_variables_AV
				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxe*size_of_real
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int
				!write variables---within loop
				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,DATATYPEX,'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh,rARRAY_PART1(1:kmaxe,i),KMAXE*PART1_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierror)
				!end write variables---within loop
				disp_in_file=disp_in_file+temp_imaxe*size_of_real
				!end loop
				end do

! 				IF (N.EQ.0)print*,"LOCATION2",disp_in_file

				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxn*size_of_real*temp_dims
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int

				call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_DOUBLE_PRECISION,DATATYPEz,'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh,rARRAY_PART4,KMAXN_P*PART4_end,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)


				disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)

! 				IF (N.EQ.0)print*,"LOCATION3",disp_in_file


				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=size_of_int*TYP_COUNTN_GLOBAL!temp_imaxe*size_of_int*temp_node
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int


				call MPI_FILE_SET_VIEW(fh,disp_in_file,MPI_INTEGER,DATATYPEy,'native',MPI_INFO_NULL,ierror)

				call MPI_FILE_WRITE_ALL(fh,iARRAY_PART2,TYP_COUNTN, MPI_INTEGER,STATUS,ierror)

				disp_in_file=disp_in_file+(size_of_int*TYP_COUNTN_GLOBAL)!(temp_imaxe*size_of_int*temp_node)

! 				IF (N.EQ.0)print*,"LOCATION4",disp_in_file

				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxe*size_of_int
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int



				call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,DATATYPEXx, 'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh, iARRAY_PART5,kmaxe*PART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)


! 				IF (N.EQ.0)print*,"LOCATION5",disp_in_file


				call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,DATATYPEINT,'native',MPI_INFO_NULL, ierror)

				IF (N.EQ.0)THEN
				BYTES=temp_imaxe*size_of_int
				nbytes=1
				Else
				BYTES=0
				nbytes=0
				end if

				call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

				disp_in_file = disp_in_file + size_of_int


				call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,DATATYPEyy, 'native',MPI_INFO_NULL, ierror)
				call MPI_FILE_WRITE_ALL(fh, iARRAY_PART3,kmaxe*PART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)


				disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)

				call MPI_FILE_CLOSE(fh, ierror)
				CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)





if (n.eq.0)then
OPEN(300,FILE=FILEX,ACCESS='STREAM',position='APPEND')
  lf = char(10)
  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)
end if


 DEallocate(vTU)


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







END SUBROUTINE PARALLEL_VTK_COMBINE_AV





SUBROUTINE PARALLEL_VTK_COMBINE_partitioned(N)
!> @brief
!> This subroutine uses MPI-IO for writing the VTK FILES
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord
CHARACTER(LEN=20)::PROC,PROC3,proc5,proc6,proc7
CHARACTER(LEN=90)::FILEX,filev
REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
LOGICAL::HERE1
REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
INTEGER                     :: nbytes,eight
CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
character(LEN=:),allocatable::VTU
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
nbytes=1

temp_cord=3
 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0
KMAXE=XMPIELRANK(N)



temp_cord=3





! 						if (n.eq.0)then
                               WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') N
                               FILEX="OUT_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".vtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)
!                             end if



	if (movement.eq.1)then
				K=1
	  DO I=1,KMAXN
		nodes_vtu(K:K+dims-1)=INODER4(I)%CORD(1:DIMS)
		if (dimensiona.eq.2)then
		nodes_vtu(K+temp_cord-1:K+temp_cord-1)=0.0D0
		end if
		K=K+TEMP_CORD
	  END DO

   end if






				if (dimensiona.eq.3)then
				DO I=1,KMAXE
				leftv(1:nof_Variables)=U_C(I)%VAL(1,1:NOF_VARIABLES)
				call CONS2PRIM(N,leftv,MP_PINFl,gammal)
					sol_vtu(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)
										do j=nof_Variables+1,write_variables-TURBULENCEEQUATIONS
                                        if (multispecies.eq.1)then
										sol_vtu(i,j)=IELEM(N,I)%REDUCE!ielem(n,i)%vortex(1)
                                        else
                                        sol_vtu(i,j)=ielem(n,i)%vortex(1)
                                        end if
                                        end do
                                        IF (TURBULENCEEQUATIONS.GT.0)THEN
                                        sol_vtu(i,write_variables)=U_CT(I)%VAL(1,1)
                                        END IF
				END DO

				temp_node=8;temp_dims=3

				end if

				if (dimensiona.eq.2)then
				DO I=1,KMAXE
				leftv(1:nof_Variables)=U_C(I)%VAL(1,1:NOF_VARIABLES)
				call cons2prim(N,leftv,MP_PINFl,gammal)
					sol_vtu(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)
										do j=nof_Variables+1,write_variables-TURBULENCEEQUATIONS
                                        if (multispecies.eq.1)then

					sol_vtu(i,j)=IELEM(N,I)%REDUCE!ielem(n,i)%vortex(1)
                                        else
                                        sol_vtu(i,j)=ielem(n,i)%vortex(1)
                                        end if
                                        end do
                                        IF (TURBULENCEEQUATIONS.GT.0)THEN
                                        sol_vtu(i,write_variables)=U_CT(I)%VAL(1,1)
                                        END IF

				END DO
				temp_node=4;temp_dims=3
				end if





temp_imaxe=kmaxe
temp_imaxn=kmaxn


	!first write the header xml file from one MPI process

	lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
	 offset_temp=0
     WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify field pieces
    WRITE(tempstamp1,'(I16)')kmaxn
    WRITE(tempstamp2,'(I16)')kmaxe
    Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
           &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    DO i=1,WRITE_VARIABLES
      Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names(i))//'" '// &
                       'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      WRITE(Offset_stamp,'(I16)')offset_temp
    END DO
    Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
    ! Specify necessary cell data
    Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
    ! Connectivity
    Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+(TYP_COUNTN)*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Offsets
    Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Elem types
    Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
    Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Prepare append section
    Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
    ! Write leading data underscore
    Buffer=TRIM('_');WRITE(300) Buffer
	Bytes = size_of_real





				!----write time stamp----!
				BYTES=size_of_real
				WRITE(300)BYTES,T

				!end time stamp



				!----write variables----!
				do j=1,WRITE_VARIABLES
				BYTES=temp_imaxe*size_of_real
					 WRITE(300)bytes,sol_vtu(1:kmaxe,j)

				end do
				!end----write variables----!


				!write nodes now!
				BYTES=kmaxn*3*size_of_real
					 WRITE(300)bytes,nodes_vtu(1:kmaxn*3)

				!end write nodes now!


				!write connectivity now!
				BYTES=TYP_COUNTN*size_of_int
					 WRITE(300)bytes, connect_vtu(1:TYP_COUNTN)

				!end write connectivity now!


				!write offsets now!
				BYTES=kmaxe*size_of_int
					 WRITE(300)bytes,offset_vtu(1:kmaxe)

				!end write nodes now!


				!write types now!
				BYTES=kmaxe*size_of_int
					 WRITE(300)bytes,type_vtu(1:kmaxe)

				!end write nodes now!

  lf = char(10)
  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)



 DEallocate(vTU)

	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


 IF (N.EQ.0)THEN


							WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') N
                               FILEX="PAR_"//TRIM(ADJUSTL(PROC3))//".pvtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)




                               lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <PUnstructuredGrid GhostLevel="0">'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf);WRITE(300) Buffer
!     Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;WRITE(300) Buffer
    WRITE(tempstamp1,'(F17.9)')T
    Buffer=TRIM('        '//TRIM(ADJUSTL(tempstamp1))//lf);WRITE(300) Buffer
	Buffer=TRIM('      </DataArray>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('    <PCellData>'//lf);WRITE(300) Buffer
    DO i=1,WRITE_VARIABLES
      Buffer=TRIM('        <PDataArray type="Float64" Name="'//TRIM(Variable_names(i))//'" '// &
                       'format="appended"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
    END DO
	Buffer=TRIM('    </PCellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    <PPoints>'//lf);WRITE(300) Buffer

     Buffer=TRIM('        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </PPoints>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    <PCells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <PDataArray type="Int32" Name="connectivity"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('        <PDataArray type="Int32" Name="offsets"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('        <PDataArray type="Int32" Name="types"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('    </PCells>'//lf);WRITE(300) Buffer
   DO procx=0,isize-1

				WRITE(PROC6,FMT='(I10)') IT
				WRITE(PROC7,FMT='(I10)') PROCX
		Buffer=TRIM('    <Piece Source="OUT_'//TRIM(ADJUSTL(PROC6))//"_"//TRIM(ADJUSTL(PROC7))//'.vtu"/>'//lf);WRITE(300) Buffer
   end do
  Buffer=TRIM('  </PUnstructuredGrid>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)



	DEallocate(vTU)



 END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




END SUBROUTINE PARALLEL_VTK_COMBINE_partitioned


SUBROUTINE PARALLEL_VTK_COMBINE_partitioned_wall(N)
!> @brief
!> This subroutine uses MPI-IO for writing the VTK FILES
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord,KKD_I,KKD
CHARACTER(LEN=20)::PROC,PROC3,proc5,proc6,proc7
CHARACTER(LEN=90)::FILEX,filev
REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
LOGICAL::HERE1
REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
INTEGER                     :: nbytes,eight
CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
character(LEN=:),allocatable::VTU
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::SHEAR_TEMP
integer::iconsidered,facex
nbytes=1

temp_cord=3
 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0


if (iloopx.gt.0)then





temp_cord=3






                               WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') N
                               FILEX="SURF_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".vtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)




			if (movement.eq.1)then
						K=1
			DO I=1,KMAXN
				nodes_vtu_W(K:K+dims-1)=INODER4(I)%CORD(1:DIMS)
				if (dimensiona.eq.2)then
				nodes_vtu_w(K+temp_cord-1:K+temp_cord-1)=0.0D0
				end if
				K=K+TEMP_CORD
			END DO

		end if








				DO I=1,iloopx
				facex=WALL_L(I,2)
				ICONSIdered=WALL_L(I,1)
				LEFTV(1:NOF_VARIABLES)=U_C(ICONSIdered)%VAL(1,1:NOF_VARIABLES)
							IF (DIMENSIONA.EQ.3)THEN
								CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
								temp_node=3;temp_dims=3

								ELSE
								CALL cons2prim(N,leftv,MP_PINFl,gammal)
								temp_node=2;temp_dims=2
								END IF



					sol_vtu_w(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)

					sol_vtu_w(I,NOF_VARIABLES+1:NOF_VARIABLES+1)=IELEM(N,ICONSIdered)%VORTEX(1)


					KKD_I=NOF_VARIABLES+1
								    IF (ITESTCASE.EQ.4)THEN

										IF (DIMENSIONA.EQ.3)THEN
											DO KKD=1,3




										select case(kkd)
											case(1)

											call SHEAR_X(ICONSIdered,facex,shear_temp)
											CASE (2)
											call SHEAR_y(ICONSIdered,facex,shear_temp)
											CASE(3)
											call SHEAR_z(ICONSIdered,facex,shear_temp)
											END SELECT

											sol_vtu_w(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF

										IF (DIMENSIONA.EQ.2)THEN
											DO KKD=1,2




										select case(kkd)
											case(1)

											call SHEAR_x2d(ICONSIdered,facex,shear_temp)
											CASE (2)
											call SHEAR_y2d(ICONSIdered,facex,shear_temp)

											END SELECT

											sol_vtu_w(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF



									END IF

								END DO






temp_imaxe=iloopx
temp_imaxn=kmaxn



	!first write the header xml file from one MPI process

	lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
	 offset_temp=0
     WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify field pieces
    WRITE(tempstamp1,'(I16)')kmaxn
    WRITE(tempstamp2,'(I16)')iloopx
    Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
           &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    DO i=1,WRITE_VARIABLES_W
      Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names_W(i))//'" '// &
                       'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      WRITE(Offset_stamp,'(I16)')offset_temp
    END DO
    Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
    ! Specify necessary cell data
    Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
    ! Connectivity
    Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+(TYP_COUNTN_w)*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Offsets
    Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Elem types
    Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
    Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Prepare append section
    Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
    ! Write leading data underscore
    Buffer=TRIM('_');WRITE(300) Buffer
	Bytes = size_of_real





				!----write time stamp----!
				BYTES=size_of_real
				WRITE(300)BYTES,T

				!end time stamp



				!----write variables----!

				do j=1,WRITE_VARIABLES_W
				BYTES=temp_imaxe*size_of_real
					 WRITE(300)bytes,sol_vtu_W(1:iloopx,j)

				end do
				!end----write variables----!


				!write nodes now!
				BYTES=kmaxn*3*size_of_real
					 WRITE(300)bytes,nodes_vtu_W(1:kmaxn*3)

				!end write nodes now!


				!write connectivity now!
				BYTES=TYP_COUNTN_w*size_of_int
					 WRITE(300)bytes, connect_vtu_W(1:TYP_COUNTN_w)

				!end write connectivity now!


				!write offsets now!
				BYTES=iloopx*size_of_int
					 WRITE(300)bytes,offset_vtu_w(1:iloopx)

				!end write nodes now!


				!write types now!
				BYTES=iloopx*size_of_int
					 WRITE(300)bytes,type_vtu_W(1:iloopx)

				!end write nodes now!

  lf = char(10)
  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)



 DEallocate(vTU)


end if


	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


 IF (N.EQ.0)THEN


							WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') N
                               FILEX="PAR_SURF_"//TRIM(ADJUSTL(PROC3))//".pvtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)




                               lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <PUnstructuredGrid GhostLevel="0">'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf);WRITE(300) Buffer
!     Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;WRITE(300) Buffer
    WRITE(tempstamp1,'(F17.9)')T
     Buffer=TRIM('        '//TRIM(ADJUSTL(tempstamp1))//lf);WRITE(300) Buffer
	Buffer=TRIM('      </DataArray>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('    <PCellData>'//lf);WRITE(300) Buffer
    DO i=1,WRITE_VARIABLES_W
      Buffer=TRIM('        <PDataArray type="Float64" Name="'//TRIM(Variable_names_W(i))//'" '// &
                       'format="appended"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
    END DO
	Buffer=TRIM('    </PCellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    <PPoints>'//lf);WRITE(300) Buffer

     Buffer=TRIM('        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </PPoints>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    <PCells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <PDataArray type="Int32" Name="connectivity"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('        <PDataArray type="Int32" Name="offsets"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('        <PDataArray type="Int32" Name="types"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('    </PCells>'//lf);WRITE(300) Buffer
   DO procx=0,isize-1
				if (wallcount_cpu_G(procx).gt.0)then
				WRITE(PROC6,FMT='(I10)') IT
				WRITE(PROC7,FMT='(I10)') PROCX
		Buffer=TRIM('    <Piece Source="SURF_'//TRIM(ADJUSTL(PROC6))//"_"//TRIM(ADJUSTL(PROC7))//'.vtu"/>'//lf);WRITE(300) Buffer
				end if
   end do
  Buffer=TRIM('  </PUnstructuredGrid>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)



	DEallocate(vTU)



 END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




END SUBROUTINE PARALLEL_VTK_COMBINE_partitioned_wall




SUBROUTINE PARALLEL_VTK_COMBINE_partitioned_wall_av(N)
!> @brief
!> This subroutine uses MPI-IO for writing the VTK FILES
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord,KKD_I,KKD,ind1
CHARACTER(LEN=20)::PROC,PROC3,proc5,proc6,proc7
CHARACTER(LEN=90)::FILEX,filev
REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
LOGICAL::HERE1
REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
INTEGER                     :: nbytes,eight
CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
CHARACTER(LEN=200)          :: Buffer
CHARACTER(LEN=1)            :: lf
character(LEN=:),allocatable::VTU
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::SHEAR_TEMP
integer::iconsidered,facex
nbytes=1

temp_cord=3
 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0

if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if



if (iloopx.gt.0)then





temp_cord=3






                               WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') N
                               FILEX="SURF_AV_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".vtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)




			if (movement.eq.1)then
						K=1
			DO I=1,KMAXN
				nodes_vtu_W(K:K+dims-1)=INODER4(I)%CORD(1:DIMS)
				if (dimensiona.eq.2)then
				nodes_vtu_w(K+temp_cord-1:K+temp_cord-1)=0.0D0
				end if
				K=K+TEMP_CORD
			END DO

		end if








				DO I=1,iloopx
				facex=WALL_L(I,2)
				ICONSIdered=WALL_L(I,1)
				LEFTV(1:NOF_VARIABLES)=U_C(ICONSIdered)%VAL(ind1,1:NOF_VARIABLES)
							IF (DIMENSIONA.EQ.3)THEN
								CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
								temp_node=3;temp_dims=3

								ELSE
								CALL cons2prim(N,leftv,MP_PINFl,gammal)
								temp_node=2;temp_dims=2
								END IF



					sol_vtu_w(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)




					KKD_I=NOF_VARIABLES
								    IF (ITESTCASE.EQ.4)THEN

										IF (DIMENSIONA.EQ.3)THEN
											DO KKD=1,3




										select case(kkd)
											case(1)

											call shear_x_av(iconsidered,facex,shear_temp)
											CASE (2)
											call shear_y_av(iconsidered,facex,shear_temp)
											CASE(3)
											call shear_z_av(iconsidered,facex,shear_temp)
											END SELECT

											sol_vtu_w(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF

										IF (DIMENSIONA.EQ.2)THEN
											DO KKD=1,2




										select case(kkd)
											case(1)

											call shear_x2d_av(iconsidered,facex,shear_temp)
											CASE (2)
											call shear_y2d_av(iconsidered,facex,shear_temp)

											END SELECT

											sol_vtu_w(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF



									END IF

								END DO






temp_imaxe=iloopx
temp_imaxn=kmaxn



	!first write the header xml file from one MPI process

	lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
	 offset_temp=0
     WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify field pieces
    WRITE(tempstamp1,'(I16)')kmaxn
    WRITE(tempstamp2,'(I16)')iloopx
    Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
           &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    DO i=1,WRITE_VARIABLES_Av_w
      Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names_AV_w(i))//'" '// &
                       'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      WRITE(Offset_stamp,'(I16)')offset_temp
    END DO
    Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    WRITE(Offset_stamp,'(I16)')offset_temp
    Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
    ! Specify necessary cell data
    Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
    ! Connectivity
    Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+(TYP_COUNTN_w)*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Offsets
    Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    WRITE(Offset_stamp,'(I16)')offset_temp
    ! Elem types
    Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
    ! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
    Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
    Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
    ! Prepare append section
    Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
    ! Write leading data underscore
    Buffer=TRIM('_');WRITE(300) Buffer
	Bytes = size_of_real





				!----write time stamp----!
				BYTES=size_of_real
				WRITE(300)BYTES,T

				!end time stamp



				!----write variables----!

				do j=1,WRITE_VARIABLES_Av_w
				BYTES=temp_imaxe*size_of_real
					 WRITE(300)bytes,sol_vtu_W(1:iloopx,j)

				end do
				!end----write variables----!


				!write nodes now!
				BYTES=kmaxn*3*size_of_real
					 WRITE(300)bytes,nodes_vtu_W(1:kmaxn*3)

				!end write nodes now!


				!write connectivity now!
				BYTES=TYP_COUNTN_w*size_of_int
					 WRITE(300)bytes, connect_vtu_W(1:TYP_COUNTN_w)

				!end write connectivity now!


				!write offsets now!
				BYTES=iloopx*size_of_int
					 WRITE(300)bytes,offset_vtu_w(1:iloopx)

				!end write nodes now!


				!write types now!
				BYTES=iloopx*size_of_int
					 WRITE(300)bytes,type_vtu_W(1:iloopx)

				!end write nodes now!

  lf = char(10)
  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)



 DEallocate(vTU)


end if


	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


 IF (N.EQ.0)THEN


							WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') N
                               FILEX="PAR_SURF_AV_"//TRIM(ADJUSTL(PROC3))//".pvtu"
                               ITRIMM=len_trim(FILEX)
                               allocate(character(LEN=ITRIMM)::VTU)
                               VTU=FILEX(1:ITRIMM)




                               lf = char(10)
   ! Write file name
    OPEN(300,FILE=VTU,ACCESS='STREAM')
    ! Write header
    Buffer=TRIM('<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf);WRITE(300) Buffer
    ! Write unstructured grid type
    Buffer=TRIM('  <PUnstructuredGrid GhostLevel="0">'//lf);WRITE(300) Buffer
    ! Write solution time type
    Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf);WRITE(300) Buffer
!     Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;WRITE(300) Buffer
    WRITE(tempstamp1,'(F17.9)')T
     Buffer=TRIM('        '//TRIM(ADJUSTL(tempstamp1))//lf);WRITE(300) Buffer
	Buffer=TRIM('      </DataArray>'//lf);WRITE(300) Buffer
	Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
    ! Specify point data
    Buffer=TRIM('    <PCellData>'//lf);WRITE(300) Buffer
    DO i=1,WRITE_VARIABLES_Av_w
      Buffer=TRIM('        <PDataArray type="Float64" Name="'//TRIM(Variable_names_AV_w(i))//'" '// &
                       'format="appended"/>'//lf);WRITE(300) Buffer
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
    END DO
	Buffer=TRIM('    </PCellData>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    <PPoints>'//lf);WRITE(300) Buffer

     Buffer=TRIM('        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    </PPoints>'//lf);WRITE(300) Buffer
    Buffer=TRIM('    <PCells>'//lf);WRITE(300) Buffer
    Buffer=TRIM('        <PDataArray type="Int32" Name="connectivity"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('        <PDataArray type="Int32" Name="offsets"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('        <PDataArray type="Int32" Name="types"/>'//lf);WRITE(300) Buffer
     Buffer=TRIM('    </PCells>'//lf);WRITE(300) Buffer
   DO procx=0,isize-1
				if (wallcount_cpu_G(procx).gt.0)then
				WRITE(PROC6,FMT='(I10)') IT
				WRITE(PROC7,FMT='(I10)') PROCX
		Buffer=TRIM('    <Piece Source="SURF_AV_'//TRIM(ADJUSTL(PROC6))//"_"//TRIM(ADJUSTL(PROC7))//'.vtu"/>'//lf);WRITE(300) Buffer
				end if
   end do
  Buffer=TRIM('  </PUnstructuredGrid>'//lf);WRITE(300) Buffer
  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
  CLOSE(300)



	DEallocate(vTU)



 END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




END SUBROUTINE PARALLEL_VTK_COMBINE_partitioned_wall_av





SUBROUTINE PARALLEL_VTK_COMBINE_partitioned_AV(N)
	!> @brief
	!> This subroutine uses MPI-IO for writing the VTK FILES
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
	INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord,ind1
	CHARACTER(LEN=20)::PROC,PROC3,proc5,proc6,proc7
	CHARACTER(LEN=90)::FILEX,filev
	REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
	LOGICAL::HERE1
	REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
	integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
	INTEGER                     :: nbytes,eight
	CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
	CHARACTER(LEN=200)          :: Buffer
	CHARACTER(LEN=1)            :: lf
	character(LEN=:),allocatable::VTU
	real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
	nbytes=1
	
	temp_cord=3
		 size_of_int=4
		 size_of_real=8
	
	offset_temp=0
	disp_in_file=0
	tmp=0
	disp_init=0
	KMAXE=XMPIELRANK(N)
	
	
	
	temp_cord=3
	

	if (rungekutta.eq.4)then
		ind1=7
		else
		ind1=5
		end if

	
	
	
	
	! 						if (n.eq.0)then
								   WRITE(PROC3,FMT='(I10)') IT
								   WRITE(PROC5,FMT='(I10)') N
								   FILEX="VOL_AVER_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".vtu"
								   ITRIMM=len_trim(FILEX)
								   allocate(character(LEN=ITRIMM)::VTU)
								   VTU=FILEX(1:ITRIMM)
	!                             end if
	
	
	
		if (movement.eq.1)then
					K=1
		  DO I=1,KMAXN
			nodes_vtu(K:K+dims-1)=INODER4(I)%CORD(1:DIMS)
			if (dimensiona.eq.2)then
			nodes_vtu(K+temp_cord-1:K+temp_cord-1)=0.0D0
			end if
			K=K+TEMP_CORD
		  END DO
	
	   end if
	
	
	
	
	
	
					if (dimensiona.eq.3)then
					DO I=1,KMAXE
					leftv(1:nof_Variables)=U_C(I)%VAL(IND1,1:NOF_VARIABLES)
					call CONS2PRIM(N,leftv,MP_PINFl,gammal)
						sol_vtu(i,1:NOF_VARIABLES)=leftv(1:nof_Variables)
						do j=nof_Variables+1,write_variables_AV
							sol_vtu(i,J)=U_C(I)%RMS(J-nof_Variables)
						end do
					END DO
	
					temp_node=8;temp_dims=3
	
					end if
	
					
	
	
	
	
	
	temp_imaxe=kmaxe
	temp_imaxn=kmaxn
	
	
		!first write the header xml file from one MPI process
	
		lf = char(10)
	   ! Write file name
		OPEN(300,FILE=VTU,ACCESS='STREAM')
		! Write header
		Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
		! Write unstructured grid type
		Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
		! Write solution time type
		Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
		 offset_temp=0
		 WRITE(Offset_stamp,'(I16)')offset_temp
		Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
						 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
		Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
		! Specify field pieces
		WRITE(tempstamp1,'(I16)')kmaxn
		WRITE(tempstamp2,'(I16)')kmaxe
		Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
			   &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
		! Specify point data
		Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
		Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
		Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
		offset_temp=offset_temp+size_of_int+size_of_real
		WRITE(Offset_stamp,'(I16)')offset_temp
		DO i=1,WRITE_VARIABLES_av
		  Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names_AV(i))//'" '// &
						   'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
		  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
		  WRITE(Offset_stamp,'(I16)')offset_temp
		END DO
		Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
		Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
		Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
						 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
		! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
		offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
		WRITE(Offset_stamp,'(I16)')offset_temp
		Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
		! Specify necessary cell data
		Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
		! Connectivity
		Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
						 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
		offset_temp=offset_temp+size_of_int+(TYP_COUNTN)*size_of_int
		WRITE(Offset_stamp,'(I16)')offset_temp
		! Offsets
		Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
						 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
		offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
		WRITE(Offset_stamp,'(I16)')offset_temp
		! Elem types
		Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
						 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
		! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
		Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
		Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
		Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
		! Prepare append section
		Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
		! Write leading data underscore
		Buffer=TRIM('_');WRITE(300) Buffer
		Bytes = size_of_real
	
	
	
	
	
					!----write time stamp----!
					BYTES=size_of_real
					WRITE(300)BYTES,T
	
					!end time stamp
	
					
	
					!----write variables----!
					do j=1,WRITE_VARIABLES_AV
					BYTES=temp_imaxe*size_of_real
						 WRITE(300)bytes,sol_vtu(1:kmaxe,j)
	
					end do
					!end----write variables----!
	
	
					!write nodes now!
					BYTES=kmaxn*3*size_of_real
						 WRITE(300)bytes,nodes_vtu(1:kmaxn*3)
	
					!end write nodes now!
	
	
					!write connectivity now!
					BYTES=TYP_COUNTN*size_of_int
						 WRITE(300)bytes, connect_vtu(1:TYP_COUNTN)
	
					!end write connectivity now!
	
	
					!write offsets now!
					BYTES=kmaxe*size_of_int
						 WRITE(300)bytes,offset_vtu(1:kmaxe)
	
					!end write nodes now!
	
	
					!write types now!
					BYTES=kmaxe*size_of_int
						 WRITE(300)bytes,type_vtu(1:kmaxe)
	
					!end write nodes now!
	
	  lf = char(10)
	  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
	  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
	  CLOSE(300)
	
	
	
	 DEallocate(vTU)
	
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	 IF (N.EQ.0)THEN
	
	
								WRITE(PROC3,FMT='(I10)') IT
								   WRITE(PROC5,FMT='(I10)') N
								   FILEX="PAR_VOL_AVER_"//TRIM(ADJUSTL(PROC3))//".pvtu"
								   ITRIMM=len_trim(FILEX)
								   allocate(character(LEN=ITRIMM)::VTU)
								   VTU=FILEX(1:ITRIMM)
	
	
	
	
								   lf = char(10)
	   ! Write file name
		OPEN(300,FILE=VTU,ACCESS='STREAM')
		! Write header
		Buffer=TRIM('<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf);WRITE(300) Buffer
		! Write unstructured grid type
		Buffer=TRIM('  <PUnstructuredGrid GhostLevel="0">'//lf);WRITE(300) Buffer
		! Write solution time type
		Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
		Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf);WRITE(300) Buffer
	!     Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;WRITE(300) Buffer
		WRITE(tempstamp1,'(F17.9)')T
		 Buffer=TRIM('        '//TRIM(ADJUSTL(tempstamp1))//lf);WRITE(300) Buffer
		Buffer=TRIM('      </DataArray>'//lf);WRITE(300) Buffer
		Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
		! Specify point data
		Buffer=TRIM('    <PCellData>'//lf);WRITE(300) Buffer
		DO i=1,WRITE_VARIABLES_AV
		  Buffer=TRIM('        <PDataArray type="Float64" Name="'//TRIM(Variable_names_AV(i))//'" '// &
						   'format="appended"/>'//lf);WRITE(300) Buffer
		  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
		END DO
		Buffer=TRIM('    </PCellData>'//lf);WRITE(300) Buffer
		Buffer=TRIM('    <PPoints>'//lf);WRITE(300) Buffer
	
		 Buffer=TRIM('        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf);WRITE(300) Buffer
		Buffer=TRIM('    </PPoints>'//lf);WRITE(300) Buffer
		Buffer=TRIM('    <PCells>'//lf);WRITE(300) Buffer
		Buffer=TRIM('        <PDataArray type="Int32" Name="connectivity"/>'//lf);WRITE(300) Buffer
		 Buffer=TRIM('        <PDataArray type="Int32" Name="offsets"/>'//lf);WRITE(300) Buffer
		 Buffer=TRIM('        <PDataArray type="Int32" Name="types"/>'//lf);WRITE(300) Buffer
		 Buffer=TRIM('    </PCells>'//lf);WRITE(300) Buffer
	   DO procx=0,isize-1
	
					WRITE(PROC6,FMT='(I10)') IT
					WRITE(PROC7,FMT='(I10)') PROCX
			Buffer=TRIM('    <Piece Source="VOL_AVER_'//TRIM(ADJUSTL(PROC6))//"_"//TRIM(ADJUSTL(PROC7))//'.vtu"/>'//lf);WRITE(300) Buffer
	   end do
	  Buffer=TRIM('  </PUnstructuredGrid>'//lf);WRITE(300) Buffer
	  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
	  CLOSE(300)
	
	
	
		DEallocate(vTU)
	
	
	
	 END IF
	
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	
	
	END SUBROUTINE PARALLEL_VTK_COMBINE_partitioned_AV





	SUBROUTINE PARALLEL_VTK_COMBINE_WALL(N)
		!> @brief
		!> This subroutine uses MPI-IO for writing the VTK FILES
		IMPLICIT NONE
		INTEGER,INTENT(IN)::N
		REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
		INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord,KKD_I,KKD,TYP_COUNTN_GLOBAL
		CHARACTER(LEN=20)::PROC,FILEX,PROC3
		REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
		LOGICAL::HERE1
		REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
		integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
		INTEGER                     :: nbytes,eight
		CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
		CHARACTER(LEN=200)          :: Buffer
		CHARACTER(LEN=1)            :: lf
		character(LEN=:),allocatable::VTU
		real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::SHEAR_TEMP
integer::iconsidered,facex
		nbytes=1
		
		
			 size_of_int=4
			 size_of_real=8
		
		offset_temp=0
		disp_in_file=0
		tmp=0
		disp_init=0
		KMAXE=XMPIELRANK(N)
		KMAXN_P=XMPIALL_v(N)
		if (dimensiona.eq.3)then
 		temp_node=3;temp_dims=3
 		else
 		temp_node=2;temp_dims=3
 		end if
		
		temp_cord=3



		
		

									   WRITE(PROC3,FMT='(I10)') IT
									   FILEX="SURF_"//TRIM(ADJUSTL(PROC3))//".vtu"
									   ITRIMM=len_trim(FILEX)
									   allocate(character(LEN=ITRIMM)::VTU)
									   VTU=FILEX(1:ITRIMM)

		
				IF (MOVEMENT.EQ.1)THEN

						k=1
						do i=1,kmaxn_P
									WrARRAY_PART4(K:k+dims-1)=INODER4(my_nodesl(i))%CORD(1:DIMS)
									if (dimensiona.eq.2)then
									wrARRAY_PART4(K+temp_cord-1:K+temp_cord-1)=0.0D0
									end if
									K=K+TEMP_CORD
						end do

				END IF
		
		
		
		
		
		
			  			!LOOP THE CORRECT NUMBER OF ELEMENTS THAT ARE BOUNDED 
			  			IF (ILOOPX.GT.0)THEN
							do i=1,ILOOPX
								facex=WALL_L(I,2)
								ICONSIdered=WALL_L(I,1)
								LEFTV(1:NOF_VARIABLES)=U_C(ICONSIdered)%VAL(1,1:NOF_VARIABLES)


								IF (DIMENSIONA.EQ.3)THEN
								CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
								temp_node=3;temp_dims=3

								ELSE
								CALL cons2prim(N,leftv,MP_PINFl,gammal)
								temp_node=2;temp_dims=3
								END IF

									WrARRAY_PART1(I,1:NOF_VARIABLES)=LEFTV(1:NOF_VARIABLES)

									!the next variable is always going to be an auxiliary
									WrARRAY_PART1(I,NOF_VARIABLES+1:NOF_VARIABLES+1)=IELEM(N,ICONSIdered)%VORTEX(1)
						
									KKD_I=NOF_VARIABLES+1
								    IF (ITESTCASE.EQ.4)THEN
									
										IF (DIMENSIONA.EQ.3)THEN
											DO KKD=1,4




										select case(kkd)
											case(1)
											
											call SHEAR_X(ICONSIdered,facex,shear_temp)
											CASE (2)
											call SHEAR_y(ICONSIdered,facex,shear_temp)
											CASE(3)
											call SHEAR_z(ICONSIdered,facex,shear_temp)

											case(4)

											call heat_X(ICONSIdered,facex,shear_temp)



											END SELECT

											WrARRAY_PART1(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF

										IF (DIMENSIONA.EQ.2)THEN
											DO KKD=1,3




										select case(kkd)
											case(1)
											
											call SHEAR_x2d(ICONSIdered,facex,shear_temp)
											CASE (2)
											call SHEAR_y2d(ICONSIdered,facex,shear_temp)
											case(3)

											call heat_x2d(ICONSIdered,facex,shear_temp)

											
											END SELECT

											WrARRAY_PART1(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF



									END IF
		
								END DO
		
							END IF




		
		
		
		
		temp_imaxe=IWMAXE !TOTAL NUMBER OF WALL ELEMENTS IN THE DOMAIN
		temp_imaxn=imaxn	 !imaxn	!WE NEED THE TOTAL NUMBER OF NODES IN THE DOMAIN
		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




		if (n.eq.0)then
		
			!first write the header xml file from one MPI process
		
			lf = char(10)
		   ! Write file name
			OPEN(300,FILE=VTU,ACCESS='STREAM')
			! Write header
			Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
			! Write unstructured grid type
			Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
			! Write solution time type
			Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
			 offset_temp=0
			 WRITE(Offset_stamp,'(I16)')offset_temp
			Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
			! Specify field pieces
			WRITE(tempstamp1,'(I16)')temp_imaxn
			WRITE(tempstamp2,'(I16)')temp_imaxe
			Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
				   &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
			! Specify point data
			Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
			Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
			Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+size_of_real
			WRITE(Offset_stamp,'(I16)')offset_temp
			DO i=1,WRITE_VARIABLES_W
			  Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names_W(i))//'" '// &
							   'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
			  WRITE(Offset_stamp,'(I16)')offset_temp
			END DO
			Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
			Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
			Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
			WRITE(Offset_stamp,'(I16)')offset_temp
			Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
			! Specify necessary cell data
			Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
			! Connectivity
			Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+TYP_COUNTN_GLOBAL_W*size_of_int
			WRITE(Offset_stamp,'(I16)')offset_temp
			! Offsets
			Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
			WRITE(Offset_stamp,'(I16)')offset_temp
			! Elem types
			Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
			Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
			Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
			Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
			! Prepare append section
			Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
			! Write leading data underscore
			Buffer=TRIM('_');WRITE(300) Buffer
			Bytes = size_of_real
			close(300)
		
		
		end if
		


		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

		
						call MPI_file_open(MPI_COMM_WORLD,VTU,MPI_MODE_WRONLY + MPI_MODE_APPEND,MPI_INFO_NULL, fh, ierror)
						call MPI_FILE_GET_POSITION(fh, disp_in_file, ierror)
						disp_init=disp_in_FILE
		
						!----write time stamp----!
						IF (N.EQ.0)THEN
						call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
						BYTES=size_of_real
						call MPI_file_write(fh, bytes, nbytes, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
						disp_in_file = disp_in_file + size_of_int
						call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
						call MPI_file_write(fh, T, nbytes, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
						disp_in_file=disp_in_file+size_of_Real
						else
						disp_in_file=disp_in_file+size_of_int+size_of_real
						end if
						!end time stamp
		


		
		
! 						do i=1,WRITE_VARIABLES
						DO I=1,WRITE_VARIABLES_W




						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,WDATATYPEINT,'native',MPI_INFO_NULL, ierror)
		
						IF (N.EQ.0)THEN
						BYTES=temp_imaxe*size_of_real
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if
		
						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		
						disp_in_file = disp_in_file + size_of_int

						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,WDATATYPEX,'native',MPI_INFO_NULL, ierror)
						call MPI_FILE_WRITE_ALL(fh,WrARRAY_PART1(:,i),iloopx*WPART1_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierror)
						!end write variables---within loop
						disp_in_file=disp_in_file+temp_imaxe*size_of_real
						!end loop
						end do
		
		! 				IF (N.EQ.0)print*,"LOCATION2",disp_in_file
		
						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)
		
						IF (N.EQ.0)THEN
						BYTES=temp_imaxn*size_of_real*temp_dims
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if
		
						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		
						disp_in_file = disp_in_file + size_of_int

						call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_DOUBLE_PRECISION,wDATATYPEz,'native',MPI_INFO_NULL, ierror)
						call MPI_FILE_WRITE_ALL(fh,wrARRAY_PART4,KMAXN_P*wPART4_end,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)

						disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)




						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)
		
						IF (N.EQ.0)THEN
						BYTES=size_of_int*TYP_COUNTN_GLOBAL_W!temp_imaxe*size_of_int*temp_node
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if
		




						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		
						disp_in_file = disp_in_file + size_of_int


						call MPI_FILE_SET_VIEW(fh,disp_in_file,MPI_INTEGER,wDATATYPEy,'native',MPI_INFO_NULL,ierror)


						call MPI_FILE_WRITE_ALL(fh,WiARRAY_PART2,TYP_COUNTN_W, MPI_INTEGER,STATUS,ierror)


		
						disp_in_file=disp_in_file+(size_of_int*TYP_COUNTN_GLOBAL_W)!(temp_imaxe*size_of_int*temp_node)
		



						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)
		
						IF (N.EQ.0)THEN
						BYTES=temp_imaxe*size_of_int
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if
		
						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		
						disp_in_file = disp_in_file + size_of_int
		
		
						CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

						call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,wDATATYPEXx, 'native',MPI_INFO_NULL, ierror)

						call MPI_FILE_WRITE_ALL(fh,WiARRAY_PART5,ILOOPX*WPART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		
						disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)
		

						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)
		
						IF (N.EQ.0)THEN
						BYTES=temp_imaxe*size_of_int
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if
		
						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		
						disp_in_file = disp_in_file + size_of_int
		
		


						call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,wDATATYPEyy, 'native',MPI_INFO_NULL, ierror)

						call MPI_FILE_WRITE_ALL(fh, wiARRAY_PART3,iloopx*WPART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)
		


		
						disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)
		
						call MPI_FILE_CLOSE(fh, ierror)
						CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		
		
		
		
		if (n.eq.0)then
		OPEN(300,FILE=FILEX,ACCESS='STREAM',position='APPEND')
		  lf = char(10)
		  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
		  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
		  CLOSE(300)
		end if
		
		
		 DEallocate(vTU)
		
		
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
		
		

		
		
		
		
		END SUBROUTINE PARALLEL_VTK_COMBINE_WALL


		SUBROUTINE PARALLEL_VTK_COMBINE_WALL_AV(N)
		!> @brief
		!> This subroutine uses MPI-IO for writing the VTK FILES
		IMPLICIT NONE
		INTEGER,INTENT(IN)::N
		REAL,ALLOCATABLE,DIMENSION(:)::array2,ARRAY3,ARRAY4
		INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,dip,N_END,ifg,kmaxn_p,ITRIMM,temp_cord,KKD_I,KKD,TYP_COUNTN_GLOBAL,ind1
		CHARACTER(LEN=20)::PROC,FILEX,PROC3
		REAL,ALLOCATABLE,DIMENSION(:)::ARRAY
		LOGICAL::HERE1
		REAL::IN1,iocpt1,iocpt2,iocpt3,iocpt4
		integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init,offset_temp,Bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
		INTEGER                     :: nbytes,eight
		CHARACTER(LEN=35)           :: Offset_stamp,tempstamp1,tempstamp2
		CHARACTER(LEN=200)          :: Buffer
		CHARACTER(LEN=1)            :: lf
		character(LEN=:),allocatable::VTU
		real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::SHEAR_TEMP
integer::iconsidered,facex
		nbytes=1


			 size_of_int=4
			 size_of_real=8

		offset_temp=0
		disp_in_file=0
		tmp=0
		disp_init=0
		KMAXE=XMPIELRANK(N)
		KMAXN_P=XMPIALL_v(N)


		temp_cord=3
	if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if





									   WRITE(PROC3,FMT='(I10)') IT
									   FILEX="SURF_AV"//TRIM(ADJUSTL(PROC3))//".vtu"
									   ITRIMM=len_trim(FILEX)
									   allocate(character(LEN=ITRIMM)::VTU)
									   VTU=FILEX(1:ITRIMM)


				IF (MOVEMENT.EQ.1)THEN

						k=1
						do i=1,kmaxn_P
									WrARRAY_PART4(K:k+dims-1)=INODER4(my_nodesl(i))%CORD(1:DIMS)
									if (dimensiona.eq.2)then
									wrARRAY_PART4(K+temp_cord-1:K+temp_cord-1)=0.0D0
									end if
									K=K+TEMP_CORD
						end do

				END IF






			  			!LOOP THE CORRECT NUMBER OF ELEMENTS THAT ARE BOUNDED
			  			IF (ILOOPX.GT.0)THEN
							do i=1,ILOOPX
								facex=WALL_L(I,2)
								ICONSIdered=WALL_L(I,1)
								LEFTV(1:NOF_VARIABLES)=U_C(ICONSIdered)%VAL(ind1,1:NOF_VARIABLES)


								IF (DIMENSIONA.EQ.3)THEN
								CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
								temp_node=3;temp_dims=3

								ELSE
								CALL cons2prim(N,leftv,MP_PINFl,gammal)
								temp_node=2;temp_dims=3
								END IF

									WrARRAY_PART1(I,1:NOF_VARIABLES)=LEFTV(1:NOF_VARIABLES)


									KKD_I=NOF_VARIABLES
								    IF (ITESTCASE.EQ.4)THEN

										IF (DIMENSIONA.EQ.3)THEN
											DO KKD=1,3




										select case(kkd)
											case(1)

											call shear_x_av(iconsidered,facex,shear_temp)
											CASE (2)
											call shear_y_av(iconsidered,facex,shear_temp)
											CASE(3)
											call shear_z_av(iconsidered,facex,shear_temp)
											END SELECT

											WrARRAY_PART1(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF

										IF (DIMENSIONA.EQ.2)THEN
											DO KKD=1,2




										select case(kkd)
											case(1)

											call shear_x2d_av(iconsidered,facex,shear_temp)
											CASE (2)
											call shear_y2d_av(iconsidered,facex,shear_temp)

											END SELECT

											WrARRAY_PART1(I,KKD_I+kkd)=SHEAR_TEMP

											END DO
										END IF



									END IF

								END DO

							END IF






		temp_imaxe=IWMAXE !TOTAL NUMBER OF WALL ELEMENTS IN THE DOMAIN
		temp_imaxn=imaxn	 !imaxn	!WE NEED THE TOTAL NUMBER OF NODES IN THE DOMAIN

		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




		if (n.eq.0)then

			!first write the header xml file from one MPI process

			lf = char(10)
		   ! Write file name
			OPEN(300,FILE=VTU,ACCESS='STREAM')
			! Write header
			Buffer=TRIM('<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf);WRITE(300) Buffer
			! Write unstructured grid type
			Buffer=TRIM('  <UnstructuredGrid>'//lf);WRITE(300) Buffer
			! Write solution time type
			Buffer=TRIM('    <FieldData>'//lf);WRITE(300) Buffer
			 offset_temp=0
			 WRITE(Offset_stamp,'(I16)')offset_temp
			Buffer=TRIM('      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			Buffer=TRIM('    </FieldData>'//lf);WRITE(300) Buffer
			! Specify field pieces
			WRITE(tempstamp1,'(I16)')temp_imaxn
			WRITE(tempstamp2,'(I16)')temp_imaxe
			Buffer=TRIM('    <Piece NumberOfPoints="'//TRIM(ADJUSTL(tempstamp1))//'" &
				   &NumberOfCells="'//TRIM(ADJUSTL(tempstamp2))//'">'//lf);WRITE(300) Buffer
			! Specify point data
			Buffer=TRIM('     <PointData>'//lf);WRITE(300) Buffer
			Buffer=TRIM('     </PointData>'//lf);WRITE(300) Buffer
			Buffer=TRIM('     <CellData>'//lf);WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+size_of_real
			WRITE(Offset_stamp,'(I16)')offset_temp
			DO i=1,WRITE_VARIABLES_Av_w
			  Buffer=TRIM('        <DataArray type="Float64" Name="'//TRIM(Variable_names_AV_w(i))//'" '// &
							   'format="appended" offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
			  WRITE(Offset_stamp,'(I16)')offset_temp
			END DO
			Buffer=TRIM('     </CellData>'//lf);WRITE(300) Buffer
			Buffer=TRIM('     <Points>'//lf);WRITE(300) Buffer
			Buffer=TRIM('        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
			WRITE(Offset_stamp,'(I16)')offset_temp
			Buffer=TRIM('     </Points>'//lf);WRITE(300) Buffer
			! Specify necessary cell data
			Buffer=TRIM('      <Cells>'//lf);WRITE(300) Buffer
			! Connectivity
			Buffer=TRIM('        <DataArray type="Int32" Name="connectivity" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+TYP_COUNTN_GLOBAL_W*size_of_int
			WRITE(Offset_stamp,'(I16)')offset_temp
			! Offsets
			Buffer=TRIM('        <DataArray type="Int32" Name="offsets" format="appended" ' // &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
			WRITE(Offset_stamp,'(I16)')offset_temp
			! Elem types
			Buffer=TRIM('        <DataArray type="Int32" Name="types" format="appended" '// &
							 'offset="'//TRIM(ADJUSTL(Offset_stamp))//'"/>'//lf);WRITE(300) Buffer
			! Buffer=TRIM('        </DataArray>'//lf;WRITE(300) Buffer
			Buffer=TRIM('      </Cells>'//lf);WRITE(300) Buffer
			Buffer=TRIM('    </Piece>'//lf);WRITE(300) Buffer
			Buffer=TRIM('  </UnstructuredGrid>'//lf);WRITE(300) Buffer
			! Prepare append section
			Buffer=TRIM('  <AppendedData encoding="raw">'//lf);WRITE(300) Buffer
			! Write leading data underscore
			Buffer=TRIM('_');WRITE(300) Buffer
			Bytes = size_of_real
			close(300)


		end if




		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


						call MPI_file_open(MPI_COMM_WORLD,VTU,MPI_MODE_WRONLY + MPI_MODE_APPEND,MPI_INFO_NULL, fh, ierror)
						call MPI_FILE_GET_POSITION(fh, disp_in_file, ierror)
						disp_init=disp_in_FILE

						!----write time stamp----!
						IF (N.EQ.0)THEN
						call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
						BYTES=size_of_real
						call MPI_file_write(fh, bytes, nbytes, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
						disp_in_file = disp_in_file + size_of_int
						call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
						call MPI_file_write(fh, T, nbytes, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
						disp_in_file=disp_in_file+size_of_Real
						else
						disp_in_file=disp_in_file+size_of_int+size_of_real
						end if
						!end time stamp





!
						DO I=1,WRITE_VARIABLES_Av_w




						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,WDATATYPEINT,'native',MPI_INFO_NULL, ierror)

						IF (N.EQ.0)THEN
						BYTES=temp_imaxe*size_of_real
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if

						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

						disp_in_file = disp_in_file + size_of_int

						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,WDATATYPEX,'native',MPI_INFO_NULL, ierror)
						call MPI_FILE_WRITE_ALL(fh,WrARRAY_PART1(:,i),iloopx*WPART1_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierror)
						!end write variables---within loop
						disp_in_file=disp_in_file+temp_imaxe*size_of_real
						!end loop
						end do

		! 				IF (N.EQ.0)print*,"LOCATION2",disp_in_file

						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)

						IF (N.EQ.0)THEN
						BYTES=temp_imaxn*size_of_real*temp_dims
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if

						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

						disp_in_file = disp_in_file + size_of_int

						call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_DOUBLE_PRECISION,wDATATYPEz,'native',MPI_INFO_NULL, ierror)
						call MPI_FILE_WRITE_ALL(fh,wrARRAY_PART4,KMAXN_P*wPART4_end,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)

						disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)




						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)

						IF (N.EQ.0)THEN
						BYTES=size_of_int*TYP_COUNTN_GLOBAL_W!temp_imaxe*size_of_int*temp_node
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if





						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

						disp_in_file = disp_in_file + size_of_int


						call MPI_FILE_SET_VIEW(fh,disp_in_file,MPI_INTEGER,wDATATYPEy,'native',MPI_INFO_NULL,ierror)


						call MPI_FILE_WRITE_ALL(fh,WiARRAY_PART2,TYP_COUNTN_W, MPI_INTEGER,STATUS,ierror)



						disp_in_file=disp_in_file+(size_of_int*TYP_COUNTN_GLOBAL_W)!(temp_imaxe*size_of_int*temp_node)




						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)

						IF (N.EQ.0)THEN
						BYTES=temp_imaxe*size_of_int
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if

						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

						disp_in_file = disp_in_file + size_of_int


						CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

						call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,wDATATYPEXx, 'native',MPI_INFO_NULL, ierror)

						call MPI_FILE_WRITE_ALL(fh,WiARRAY_PART5,ILOOPX*WPART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

						disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)


						call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_INTEGER,wDATATYPEINT,'native',MPI_INFO_NULL, ierror)

						IF (N.EQ.0)THEN
						BYTES=temp_imaxe*size_of_int
						nbytes=1
						Else
						BYTES=0
						nbytes=0
						end if

						call MPI_FILE_WRITE_ALL(fh,bytes,nbytes,MPI_INTEGER,MPI_STATUS_IGNORE, ierror)

						disp_in_file = disp_in_file + size_of_int




						call MPI_FILE_SET_VIEW(fh, disp_in_file,MPI_INTEGER,wDATATYPEyy, 'native',MPI_INFO_NULL, ierror)

						call MPI_FILE_WRITE_ALL(fh, wiARRAY_PART3,iloopx*WPART1_end, MPI_INTEGER,MPI_STATUS_IGNORE, ierror)




						disp_in_file=disp_in_file+(temp_imaxe*size_of_INT)

						call MPI_FILE_CLOSE(fh, ierror)
						CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




		if (n.eq.0)then
		OPEN(300,FILE=FILEX,ACCESS='STREAM',position='APPEND')
		  lf = char(10)
		  Buffer=TRIM(lf//'  </AppendedData>'//lf);WRITE(300) Buffer
		  Buffer=TRIM('</VTKFile>'//lf);WRITE(300) Buffer
		  CLOSE(300)
		end if


		 DEallocate(vTU)


		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)







		END SUBROUTINE PARALLEL_VTK_COMBINE_WALL_AV



SUBROUTINE CHECKPOINTv2(N)
!> @brief
!> This subroutine is writing the checkpointing files
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS
REAL,ALLOCATABLE,DIMENSION(:,:)::xbin
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT
 KMAXE=XMPIELRANK(N)


ICPUID=N
	RESTFILE='RESTART2.dat'
	IF (N.EQ.0)THEN
	OPEN(1086,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
	if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	write(1086)it,ZERO
	      write(1086)INITIALRES(1)
	      write(1086)INITIALRES(2)
	      write(1086)INITIALRES(3)
	      write(1086)INITIALRES(4)
	      write(1086)INITIALRES(5)
	      if ( turbulence .eq. 1) then
	      if (turbulencemodel.eq.1)then
	      write(1086)INITIALRES(6)
	      end if
	      if (turbulencemodel.eq.2)then
	      write(1086)INITIALRES(6)
	      write(1086)INITIALRES(7)
	      end if
	      end if
	      
	      
	ELSE
	WRITE (1086)IT,T
	
	if (initcond.eq.95)then
	write(1086)taylor
	end if
	      
	end if
	end if
	
	KMAXE=XMPIELRANK(N)
    
DUMG=KMAXE


call mpi_barrier(mpi_comm_world,IERROR) !NOT NEEDED

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


    IF (N.EQ.0)THEN
	    ALLOCATE(ICELLA(IMAXP*ISIZE))
	    ICELLA=0

    END IF

    call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

!     call mpi_barrier(mpi_comm_world,IERROR)
    deallocate (icell)

    IF (N.EQ.0) then
    ALLOCATE(VALUESA(IMAXP*ISIZE))
    allocate(xbin(imaxe,5+turbulenceequations+passivescalar))
    VALUESA=ZERO
    END IF
    
    
    ALLOCATE(VALUESS(imaxp));VALUESS=ZERO

  
  
IF (TURBULENCE.EQ.1)THEN
do jj=1,5+turbulenceequations+passivescalar
 DO I=1,KMAXE
      IF (jj.gt.5) THEN

		      VALUESS(I)=U_CT(I)%VAL(1,jj-5)
      Else
	VALUESS(I)=U_C(I)%VAL(1,jj)
      end if
 END DO
    call MPI_GATHER(VALUESS,imaxp,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i),jj)=valuesa(i)
	end if
    end do
    end if
    
end do
ELSE
  do jj=1,5
	DO I=1,KMAXE
		VALUESS(I)=U_C(I)%VAL(1,jj)
	END DO
	call MPI_GATHER(VALUESS,imaxp,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
	
	    IF (N.EQ.0)THEN
	    do i=1,imaxp*isize
		if (icella(i).gt.0)then
		xbin(icella(i),jj)=valuesa(i)
		end if
	    end do
	    end if
    
  end do

END IF

!   CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN

    DO I=1,IMAXE
    WRITE(1086)XBIN(i,1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR)
    END DO

    DEALLOCATE(XBIN,ICELLA,VALUESA)
    close(1086)
    END IF
    
DEALLOCATE(VALUESS)








END SUBROUTINE CHECKPOINTv2




SUBROUTINE CHECKPOINT2D(N)
!> @brief
!> This subroutine uses MPI-IO for writing the checkpointing files for 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA,dispt
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS,array2
REAL,ALLOCATABLE,DIMENSION(:,:)::xbin
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,fh,size_of_real,size_of_int,dip,N_END,datatype
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT,ARRAY
REAL::IN1
logical::here1
integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init
disp_in_file=0
tmp=0
disp_init=0
 KMAXE=XMPIELRANK(N)


size_of_int=4
size_of_real=8
ICPUID=N
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

IF (DG.EQ.1)THEN
ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1)))
ELSE
ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)))	!I ALLOCATE IN MEMORY THE PATTERN OF ACCESS OF DATA IN TERMS OF DISPLACEMENT, AND IN TERMS OF BLOCKLENGTH, AND FINALY AN ARRAY WITH THIS PROCESSOR DATA
END IF

       if (dg.eq.1)then
     DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*((NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1))
      END DO
      

      n_end=(NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1)
     
     else


      DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+turbulenceequations+passivescalar)
      END DO

      n_end=NOF_VARIABLES+turbulenceequations+passivescalar

      end if
      
      
      
      if (dg.eq.1)then
      IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	    K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	      ARRAY2(K:K+turbulenceequations+passivescalar-1)=U_CT(I)%VAL(1,1:turbulenceequations+passivescalar)
	      K=K+turbulenceequations+passivescalar
	  END DO
      ELSE
	  K=1
	  DO I=1,KMAXE
        DO J=1,NOF_VARIABLES
	      ARRAY2(K:K+IDEGFREE)=U_C(I)%VALDG(1,J,1:IDEGFREE+1)
	     
	      
	      K=K+(IDEGFREE+1)
        END DO
	  END DO
      END IF
      
      
      else
      
      IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	    K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	      ARRAY2(K:K+turbulenceequations+passivescalar-1)=U_CT(I)%VAL(1,1:turbulenceequations+passivescalar)
	      K=K+turbulenceequations+passivescalar
	  END DO
      ELSE
	  K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	  END DO
      END IF
      
      end if
      RESTFILE='RESTART.dat'
       IF (N.EQ.0)THEN
	  !INQUIRE (FILE=RESTFILE,EXIST=HERE1)
	  !IF (HEREss) THEN
	 INQUIRE (FILE=RESTFILE,EXIST=HERE1)
	 
	  
	  CALL MPI_FILE_DELETE(RESTFILE,MPI_INFO_NULL,IERROR)
	  
	 

	  !END IF

   END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
      
    !CREATE TYPE FIRST OF INDEXED BLOCK
    CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXE,n_end,DISPT,MPI_DOUBLE_PRECISION,DATATYPE,IERROR)
    CALL MPI_TYPE_COMMIT(DATATYPE,IERROR)
    

    ALLOCATE(ARRAY(1:nof_Variables+turbulenceequations+passivescalar))
    
	
	
	
	call MPI_file_open(MPI_COMM_WORLD, RESTFILE,MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, fh, ierror)
	
	
	 
	  
	
	if (n.eq.0)then
	
	    if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	          call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_write(fh, it, 1, MPI_INTEGER, MPI_STATUS_IGNORE,ierror)
		    disp_in_file = disp_in_file + size_of_int	!1
		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_write(fh, INITIALRES(1:nof_Variables+turbulenceequations),nof_Variables+turbulenceequations, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_Variables+turbulenceequations)	!3
	    ELSE
		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_write(fh, it, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
		    disp_in_file = disp_in_file + size_of_int 	!4
		    
		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_write(fh, T, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
		    disp_in_file = disp_in_file + size_of_real    !5
		      if (initcond.eq.95)then
		      call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
		      call MPI_file_write(fh, TAYLOR, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
			    disp_in_file = disp_in_file + size_of_real !6
		      end if  
	    end if
	ELSE
	      if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
		disp_in_file = disp_in_file + size_of_int 
		disp_in_file = disp_in_file + size_of_real*(nof_Variables+turbulenceequations)
	      ELSE
		  disp_in_file = disp_in_file + size_of_int 
		  disp_in_file = disp_in_file + size_of_real
		  if (initcond.eq.95)then
		  disp_in_file = disp_in_file + size_of_real 
		  END IF
	      END IF
	      
	      
	END IF
	
	
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,datatype, 'native',MPI_INFO_NULL, ierror)
	call MPI_FILE_WRITE_ALL(fh, ARRAY2, KMAXE*n_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)        
        call MPI_FILE_CLOSE(fh, ierror)
	CALL MPI_TYPE_FREE(DATATYPE,IERROR)
          DEALLOCATE(ARRAY,DISPT,ARRAY2)
          
          
	
	
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	
	
	



END SUBROUTINE CHECKPOINT2D




SUBROUTINE CHECKPOINTAV(N) 
!> @brief
!> This subroutine uses MPI-IO for writing the averaged checkpointing files
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA,dispt
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS,array2
REAL,ALLOCATABLE,DIMENSION(:,:)::xbin
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,IND1,fh,size_of_real,size_of_int,dip,ISTA,IEND,N_END,datatype
 CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
 REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT,array
 integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp,disp_init
 logical::here1
 KMAXE=XMPIELRANK(N)
disp_in_file=0
tmp=0
disp_init=0
 
 
 
 
 
 
 size_of_int=4
size_of_real=8
ICPUID=N
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 
IF (RUNGEKUTTA.EQ.4)THEN
IND1=7
ELSE
IND1=5
END IF


	
	
 ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar+6+passivescalar)))
    DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+turbulenceequations+passivescalar+6+passivescalar)
      END DO

 RESTFILE='RESTARTav.dat'
 
  n_end=NOF_VARIABLES+turbulenceequations+passivescalar+6+passivescalar
  
  
 
 
  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	    K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(IND1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	      ARRAY2(K:K+turbulenceequations+passivescalar-1)=U_CT(I)%VAL(IND1,1:turbulenceequations+passivescalar)
	      K=K+turbulenceequations+passivescalar
	      ARRAY2(K:K+6+PASSIVESCALAR-1)=U_C(i)%RMS(1:6+passivescalar)
	      K=K+6+PASSIVESCALAR
	  END DO
      ELSE
	  K=1
	  DO I=1,KMAXE
	      ARRAY2(K:K+NOF_VARIABLES-1)=U_C(I)%VAL(ind1,1:NOF_VARIABLES)
	      K=K+NOF_VARIABLES
	      ARRAY2(K:K+6+PASSIVESCALAR-1)=U_C(i)%RMS(1:6+passivescalar)
	      K=K+6+PASSIVESCALAR
	  END DO
      END IF
  IF (N.EQ.0)THEN
	 INQUIRE (FILE=RESTFILE,EXIST=HERE1)
	  CALL MPI_FILE_DELETE(RESTFILE,MPI_INFO_NULL,IERROR)
	END IF
 
 CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 
 
 CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXE,n_end,DISPT,MPI_DOUBLE_PRECISION,DATATYPE,IERROR)
    CALL MPI_TYPE_COMMIT(DATATYPE,IERROR)
 
 
 
	
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	call MPI_file_open(MPI_COMM_WORLD, RESTFILE,MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL, fh, ierror)
	call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,datatype, 'native',MPI_INFO_NULL, ierror)
	call MPI_FILE_WRITE_ALL(fh, ARRAY2, KMAXE*n_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)        
        call MPI_FILE_CLOSE(fh, ierror)
	CALL MPI_TYPE_FREE(DATATYPE,IERROR)
          
          
          
	DEALLOCATE(DISPT,ARRAY2)
	
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	
	
	
	
 
	

END SUBROUTINE CHECKPOINTAV


SUBROUTINE REST_READ(N)
!> @brief
!> This subroutine uses MPI-IO for reading the checkpointing files
IMPLICIT NONE
integer,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA,dispt
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS,array2
REAL,allocatable,DIMENSION(:)::RG,ARG
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
INTEGER:: prev_turbequation,INITIAL,III,i,k,j,jx,QQP,INC,kmaxe,jkn,ki,iterr,JX2,ind1,fh,size_of_real,size_of_int,dip,N_END,datatype
REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT,ARRAY
integer(kind=MPI_OFFSET_KIND) :: disp_in_file, tmp
logical::here
disp_in_file=0
tmp=0

KMAXE=XMPIELRANK(N)

if (rungekutta.eq.4)then
ind1=7
else
ind1=5
end if
prev_turbequation=0
if (prev_turbmodel.eq.1) then
prev_turbequation=1
end if 
if (prev_turbmodel.eq.2) then
prev_turbequation=2
end if 

 size_of_int=4
size_of_real=8

!$OMP MASTER
IF (DG.EQ.1)THEN
ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+prev_turbequation+passivescalar)*(IDEGFREE+1)))
ELSE
ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+prev_turbequation+passivescalar)))	
END IF

    if (dg.eq.1)then
     DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*((NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1))
      END DO
    n_end=(NOF_VARIABLES+turbulenceequations+passivescalar)*(IDEGFREE+1)
     
     else


    DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+prev_turbequation+passivescalar)
      END DO
      
      
      n_end=NOF_VARIABLES+prev_turbequation+LAMPS
      
      end if
      
      CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXE,n_end,DISPT,MPI_DOUBLE_PRECISION,DATATYPE,IERROR)
    CALL MPI_TYPE_COMMIT(DATATYPE,IERROR)
    
    RESTFILE='RESTART.dat'
    
    call MPI_file_open(MPI_COMM_WORLD, RESTFILE,MPI_MODE_RDONLY,MPI_INFO_NULL, fh, ierror)
    
    
    
	    if (IRES_UNSTEADY.eq.0)then
! 	   if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	          call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_READ(fh, it, 1, MPI_INTEGER, MPI_STATUS_IGNORE,ierror)
		    disp_in_file = disp_in_file + size_of_int	!1
		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_READ(fh, INITIALRES(1:nof_Variables+prev_turbequation),nof_Variables+prev_turbequation, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_Variables+prev_turbequation)	!3
	    ELSE
		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_READ(fh, it, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierror)
		  
		    disp_in_file = disp_in_file + size_of_int 	!4
		    
		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
		  call MPI_file_READ(fh, T, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
		 
		    disp_in_file = disp_in_file + size_of_real    !5
		      if (initcond.eq.95)then
		      call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET,ierror)
		      call MPI_file_READ(fh, TAYLOR, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
			    disp_in_file = disp_in_file + size_of_real !6
		      end if  
	    end if
      
      
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,datatype, 'native',MPI_INFO_NULL, ierror)
	call MPI_FILE_READ_ALL(fh, ARRAY2, KMAXE*n_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)        
        call MPI_FILE_CLOSE(fh, ierror)
	CALL MPI_TYPE_FREE(DATATYPE,IERROR)
	
	
	
	if (dg.eq.1)then
	
	
	IF ((prev_turbmodel.GT.0).OR.(LAMPS.GT.0))THEN
	
	    K=1
	    DO I=1,KMAXE
		U_C(I)%VAL(1,1:NOF_VARIABLES)=ARRAY2(K:K+NOF_VARIABLES-1)
		K=K+NOF_VARIABLES
		    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		    U_CT(I)%VAL(1,1:turbulenceequations+passivescalar)=ARRAY2(K:K+prev_turbequation+LAMPS-1)
		    END IF
		K=K+prev_turbmodel+LAMPS
		
	    END DO
	ELSE
	    K=1
	    DO I=1,KMAXE
	    DO J=1,NOF_VARIABLES
	    
		U_C(I)%VALdg(1,j,1:IDEGFREE+1)=ARRAY2(K:K+IDEGFREE)
		
		
		K=K+(IDEGFREE+1)
		END DO
		
		
		      IF (TURBULENCE.EQ.1)THEN
			    IF (TURBULENCEMODEL.EQ.1)THEN
				U_CT(I)%VAL(1,1)=VISC*TURBINIT
			    ELSE
				U_CT(I)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
				U_CT(I)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
					/L_TURB_INLET*RG(1)
			    END IF
			ENDIF
		K=K+prev_turbmodel+LAMPS
		
	    END DO
	END IF
	
	
	
	
	
	
	
	
	
	
	else
	IF ((prev_turbmodel.GT.0).OR.(LAMPS.GT.0))THEN
	
	    K=1
	    DO I=1,KMAXE
		U_C(I)%VAL(1,1:NOF_VARIABLES)=ARRAY2(K:K+NOF_VARIABLES-1)
		K=K+NOF_VARIABLES
		    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		    U_CT(I)%VAL(1,1:turbulenceequations+passivescalar)=ARRAY2(K:K+prev_turbequation+LAMPS-1)
		    END IF
		K=K+prev_turbmodel+LAMPS
		
	    END DO
	ELSE
	    K=1
	    DO I=1,KMAXE
		U_C(I)%VAL(1,1:NOF_VARIABLES)=ARRAY2(K:K+NOF_VARIABLES-1)
		K=K+NOF_VARIABLES
		      IF (TURBULENCE.EQ.1)THEN
			    IF (TURBULENCEMODEL.EQ.1)THEN
				U_CT(I)%VAL(1,1)=VISC*TURBINIT
			    ELSE
				U_CT(I)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
				U_CT(I)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
					/L_TURB_INLET*RG(1)
			    END IF
			ENDIF
		K=K+prev_turbmodel+LAMPS
		
	    END DO
	END IF
	
	
	end if
	
	
	
	
	
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	DEALLOCATE(DISPT,ARRAY2)
	
	
      
if (Averaging .EQ. 1) then
 
IF (Average_restart.EQ.1)THEN

disp_in_file=0
ALLOCATE(DISPT(KMAXE),ARRAY2(kmaxe*(NOF_VARIABLES+prev_turbequation+passivescalar+6+passivescalar)))
DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+prev_turbequation+passivescalar+6+passivescalar)
      END DO

      RESTFILE='RESTARTav.dat'
      n_end=NOF_VARIABLES+prev_turbequation+passivescalar+6+passivescalar
      
       CALL MPI_TYPE_CREATE_INDEXED_BLOCK(KMAXE,n_end,DISPT,MPI_DOUBLE_PRECISION,DATATYPE,IERROR)
    CALL MPI_TYPE_COMMIT(DATATYPE,IERROR)
    call MPI_file_open(MPI_COMM_WORLD, RESTFILE,MPI_MODE_RDONLY,MPI_INFO_NULL, fh, ierror)
	call MPI_FILE_SET_VIEW(fh, disp_in_file, MPI_DOUBLE_PRECISION,datatype, 'native',MPI_INFO_NULL, ierror)
	call MPI_FILE_READ_ALL(fh, ARRAY2, KMAXE*n_end, MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE, ierror)        
        call MPI_FILE_CLOSE(fh, ierror)
	CALL MPI_TYPE_FREE(DATATYPE,IERROR)
	
	IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
	    K=1
	  DO I=1,KMAXE
	      U_C(I)%VAL(IND1,1:NOF_VARIABLES)=ARRAY2(K:K+NOF_VARIABLES-1)
	      K=K+NOF_VARIABLES
	      U_CT(I)%VAL(IND1,1:turbulenceequations+passivescalar)=ARRAY2(K:K+turbulenceequations+passivescalar-1)
	      K=K+turbulenceequations+passivescalar
	      U_C(i)%RMS(1:6+passivescalar)=ARRAY2(K:K+6+PASSIVESCALAR-1)
	      K=K+6+PASSIVESCALAR
	  END DO
      ELSE
	  K=1
	  DO I=1,KMAXE
	      U_C(I)%VAL(ind1,1:NOF_VARIABLES)=ARRAY2(K:K+NOF_VARIABLES-1)
	      
	      K=K+NOF_VARIABLES
	      U_C(i)%RMS(1:6+passivescalar)=ARRAY2(K:K+6+PASSIVESCALAR-1)
	      
	      K=K+6+PASSIVESCALAR
	  END DO
      END IF
	
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	DEALLOCATE(DISPT,ARRAY2)
	
    
ELSE
DO I=1,kmaxe
      U_C(i)%VAL(ind1,:)=ZERO
      U_C(i)%RMS(:)=ZERO
	if ((passivescalar.gt.0).or.(turbulence.eq.1))then
	U_CT(i)%VAL(ind1,:)=ZERO
	end if   
END DO


END IF

END IF
!$OMP END MASTER













END SUBROUTINE REST_READ







SUBROUTINE CHECKPOINTAV2d(N)  
!> @brief
!> This subroutine writes the average checkpointing files in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS
REAL,ALLOCATABLE,DIMENSION(:,:)::xbin
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj
 CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
 REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT
 KMAXE=XMPIELRANK(N)



ICPUID=N
	RESTFILE='RESTARTav.dat'
	IF (N.EQ.0)THEN
	OPEN(1086,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
	
	
	end if
	
	KMAXE=XMPIELRANK(N)
    
DUMG=KMAXE
call mpi_barrier(mpi_comm_world,IERROR)

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


IF (N.EQ.0)THEN
	ALLOCATE(ICELLA(IMAXP*ISIZE))
	 ICELLA=0

END IF

call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)


call mpi_barrier(mpi_comm_world,IERROR)
deallocate (icell)

IF (N.EQ.0) then
ALLOCATE(VALUESA(IMAXP*ISIZE))
  allocate(xbin(imaxe,(4+turbulenceequations+passivescalar+3+passivescalar)))
	VALUESA=0.0

 END IF
ALLOCATE(VALUESS(imaxp))
  VALUESS=0.0

do jj=1,4+turbulenceequations+passivescalar+3+passivescalar
 DO I=1,KMAXE
      
      if (jj.le.4+turbulenceequations+passivescalar)then
      IF (jj.le.4) THEN
	  VALUESS(I)=U_C(I)%VAL(5,jj)
		      
      Else

	VALUESS(I)=U_CT(I)%VAL(5,jj-4)

      end if
      end if
      if (jj.gt.4+turbulenceequations+passivescalar)then

       VALUESS(I)=U_C(i)%RMS(jj-(4+turbulenceequations+passivescalar))

      end if
      
 END DO

    call MPI_GATHER(VALUESS,imaxp,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i),jj)=valuesa(i)
	end if
    end do
    end if
    
    
 end do
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN

    DO I=1,IMAXE
     WRITE(1086)I
!     DO NVAR=1,5+TURBULENCEEQUATIONS+PASSIVESCALAR
    WRITE(1086)XBIN(xmpi_re(i),1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR+3+PASSIVESCALAR)
!     END DO
    END DO
      
    DEALLOCATE(XBIN,ICELLA,VALUESA)

      close(1086)
    END IF
    

DEALLOCATE(VALUESS)
	

END SUBROUTINE CHECKPOINTAV2d


SUBROUTINE PROBING
!> @brief
!> This subroutine writes the primitve variables at the probe positions
IMPLICIT NONE
INTEGER::INV
CHARACTER(LEN=120)::PROB,PROBFILE,PROC3
LOGICAL::HERES
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
                    IF (nof_variables.GT.1)THEN

			IF (NPROBES.GT.0)THEN
			    
			    DO INV=1,NPROBES
			    IF (PROBEI(N,INV).NE.0) THEN
			      WRITE(PROB,FMT='(I10)') INV
			      PROBFILE='PROBE.'//TRIM(ADJUSTL(PROB))

			      INQUIRE (FILE=PROBFILE,EXIST=HEREs)
			    IF (HEREs.EQV..TRUE.) THEN
				OPEN(3000+N,FILE=PROBFILE,FORM='FORMATTED',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
				
				
				ELSE
				OPEN(3000+N,FILE=PROBFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
				
				END IF
				IF (PASSIVESCALAR.EQ.0)THEN
				LEFTV(1:NOF_vARIABLES)=U_C(PROBEI(N,INV))%VAL(1,1:NOF_vARIABLES)
				CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
	WRITE(3000+N,'(1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5)
				ELSE
				LEFTV(1:NOF_vARIABLES)=U_C(PROBEI(N,INV))%VAL(1,1:NOF_vARIABLES)
				CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
	WRITE(3000+N,'(1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5)&
	,U_CT(PROBEI(N,INV))%VAL(1,1)/U_C(PROBEI(N,INV))%VAL(1,1)


				END IF
				CLOSE(3000+N)
				
		    
			      END IF     
			    END DO
			  END IF
			  END IF

END SUBROUTINE PROBING


SUBROUTINE PROBING2D
!> @brief
!> This subroutine writes the primitve variables at the probe positions in 2D
IMPLICIT NONE
INTEGER::INV
CHARACTER(LEN=120)::PROB,PROBFILE,PROC3
LOGICAL::HERES
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
                        IF(nof_variables.GT.1)THEN

			IF (NPROBES.GT.0)THEN
			    
			    DO INV=1,NPROBES
			    IF (PROBEI(N,INV).NE.0) THEN
			      WRITE(PROB,FMT='(I10)') INV
			      PROBFILE='PROBE.'//TRIM(ADJUSTL(PROB))

			      INQUIRE (FILE=PROBFILE,EXIST=HEREs)
			    IF (HEREs.eqv..TRUE.) THEN
				OPEN(3000+N,FILE=PROBFILE,FORM='FORMATTED',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
				
				
				ELSE
				OPEN(3000+N,FILE=PROBFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
				
				END IF
				IF (PASSIVESCALAR.EQ.0)THEN
	WRITE(3000+N,'(1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,U_C(PROBEI(N,INV))%VAL(1,1),&
	U_C(PROBEI(N,INV))%VAL(1,2)/U_C(PROBEI(N,INV))%VAL(1,1),&
	U_C(PROBEI(N,INV))%VAL(1,3)/U_C(PROBEI(N,INV))%VAL(1,1)
				ELSE
	WRITE(3000+N,'(1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,U_C(PROBEI(N,INV))%VAL(1,1),&
	  U_C(PROBEI(N,INV))%VAL(1,2)/U_C(PROBEI(N,INV))%VAL(1,1),&
	U_C(PROBEI(N,INV))%VAL(1,3)/U_C(PROBEI(N,INV))%VAL(1,1)&
	,U_CT(PROBEI(N,INV))%VAL(1,1)/U_C(PROBEI(N,INV))%VAL(1,1)


				END IF
				CLOSE(3000+N)
				
		    
			      END IF     
			    END DO
			  END IF
			  END IF

END SUBROUTINE PROBING2D



SUBROUTINE COMPUTEFORCE(N)
!> @brief
!> This subroutine computes the forces on the wall
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,J,KMAXE,gqi_points,nnd
INTEGER:: MYSURFACE
CHARACTER(LEN=12)::PROC,RESTFILE,PROC3
REAL::DRAG,LIFT,CD,CL,RX,PX,EX,surface_temp,RTEMP,FX,FY,FZ,MX,MY,MZ
REAL::FORCEXFR,SSX,CDF,LIFTF,DRAGF,FRICTIONF,TAUYX,TAUZX,TAUZY,SSY,SSZ,CF,TAUXX,TAUYY,TAUZZ
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,angle1,angle2,nx,ny,nz
 REAL,DIMENSION(3)::CI,CO
 logical::heref
 INTEGER::IM
 REAL::TSOLR,TSOLU,TSOLE,TSOLV,TSOLW,TSOLP,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA,1:NUMBEROFPOINTS2)::QPOINTS2D
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
real,dimension(1:4)::viscl,laml
FORCEX=zero; FORCEY=zero; FORCEZ=zero;  FORCEXFR=zero
 CD=zero
 CL=zero
 CI(:)=zero
 CO(:)=zero
 FX=zero
 FY=zero
 FZ=zero
 MX=ZERO
 MY=ZERO
 MZ=ZERO
 momenty=ZERO 
 momentz=ZERO
 MOMENTX=ZERO
 KMAXE=XMPIELRANK(N)
 
!$OMP BARRIER 
!$OMP DO  REDUCTION(+:FORCEX,FORCEY,FORCEZ,momentx,momenty,momentz)
DO I=1,kmaxe
		if (ielem(n,i)%interior.eq.1)then
			IF(MRF.EQ.1)THEN
				MYSURFACE=ILOCAL_RECON3(I)%MRF
			else
				MYSURFACE=1	
			END IF	
		IF(MYSURFACE.EQ.1)THEN
		    do j=1,ielem(n,i)%ifca
		      if (ielem(n,i)%ibounds(j).gt.0)then
			  if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4))then
			      ANGLE1=IELEM(N,I)%FACEANGLEX(j)
			      ANGLE2=IELEM(N,I)%FACEANGLEY(j)
			      NX=(COS(ANGLE1)*SIN(ANGLE2))
			      NY=(SIN(ANGLE1)*SIN(ANGLE2))
			      NZ=(COS(ANGLE2))
			      
			  SSX=zero; SSP=zero; SSY=zero; SSZ=zero
			  
				select case(ielem(n,i)%types_faces(j))
				case (5)
					  gqi_points=qp_quad_n
					
					  
					  
					  if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					    NND=4
				      do K=1,nnd
					VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO
					  call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  end if
					  surface_temp=IELEM(N,I)%SURF(J)
					  
				  
				case(6)
					gqi_points=qp_triangle_n
					
					    
					if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					  NND=3
					do K=1,nnd
					  VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
					END DO
					call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					end if
 					    surface_temp=IELEM(N,I)%SURF(J)
 					    
 					    
					  
				end select
				  
				  
				  do im=1,gqi_points
				  
				  if (itestcase.eq.4)then
				  if (ielem(n,i)%ggs.eq.1)then
				  
				  VORTET1(1:3,1:3) = ILOCAL_RECON3(I)%GRADS(1:3,1:3)
				  ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
				  vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
				  wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)
				  
				  else
				  
				  vortet1(1,1:3)=ILOCAL_RECON3(i)%ULEFTV(1:3,2,J,IM)
				  vortet1(2,1:3)=ILOCAL_RECON3(i)%ULEFTV(1:3,3,J,IM)
				  vortet1(3,1:3)=ILOCAL_RECON3(i)%ULEFTV(1:3,4,J,IM)
				  ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
				  vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
				  wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)
				  
				  
				  end if
				 end if
				  
				  IF (DG.EQ.1)THEN
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT_DG(1:nof_Variables, J,IM)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT_DG(1:nof_Variables, J,IM)


				  ELSE
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  END IF
				  
				  
				  
				    call CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
				    px=leftv(5)
				    ssp=ssp+(px*WEQUA2D(im))
				    if (itestcase.eq.4)then
				    CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
				  
				  
				  TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
				  TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
				  TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
				  TAUYX=(UY + VX)
				  TAUZX=(WX + UZ)
				  TAUZY=(VZ + WY)
				  SSX=SSX-((VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))*WEQUA2D(im))
				  SSY=SSY-((VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))*WEQUA2D(im))
				  SSZ=SSZ-((VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))*WEQUA2D(im))
				 end if
				   end do
				   
				   
				  SSP=ssp-PRES	
				  
				  FORCEX=FORCEX+(((SSP)*(surface_temp)*NX))+((SSX)*surface_temp)
				  FORCEY=FORCEY+(((SSP)*(surface_temp)*NY))+((SSY)*surface_temp)
				  FORCEZ=FORCEZ+(((SSP)*(surface_temp)*NZ))+((SSZ)*surface_temp)

				  MOMENTX=MOMENTX+(((SSP)*(surface_temp)*NZ))*IELEM(N,I)%YYC-(((SSP)*(surface_temp)*NY))*IELEM(N,I)%ZZC
				  MOMENTY=MOMENTY+ (((SSP)*(surface_temp)*NX))*IELEM(N,I)%ZZC-(((SSP)*(surface_temp)*NZ))*IELEM(N,I)%XXC	
				  MOMENTZ=MOMENTZ+(((SSP)*(surface_temp)*NY))*IELEM(N,I)%XXC-(((SSP)*(surface_temp)*NX))*IELEM(N,I)%YYC
					    
			END IF
		      end if
		    end do
		END IF !MYSURFACE
		end if

			
		
END DO					 
!$OMP END DO
	

		
	
	
!$OMP BARRIER 
!$OMP MASTER 
	FORCEX=FORCEX*VECTORX
	FORCEY=FORCEY*VECTORY
	FORCEZ=FORCEZ*VECTORZ
	IF(RFRAME.EQ.0)THEN
        RTEMP=((AOA/180.0d0)*PI)
	LIFTF=(FORCEZ*COS(RTEMP))+(FORCEy*COS(RTEMP))-(FORCEX*SIN(RTEMP))
	DRAGF=(FORCEX*COS(RTEMP))+(FORCEy*SIN(RTEMP))+(FORCEZ*SIN(RTEMP))
	CL=(2.0D0*LIFTF)/((RRES)*(ufreestream**2))
	CD=(2.0D0*DRAGF)/((RRES)*(ufreestream**2))


	CO(1)=CL
	CO(2)=CD
	CALL MPI_ALLREDUCE(CO(1:2),CI(1:2),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
	CL=CI(1)
	CD=CI(2)
	ELSE
        CO(1)=FORCEX
        CO(2)=FORCEY
        CO(3)=FORCEZ
        CALL MPI_ALLREDUCE(CO(1:3),CI(1:3),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        FX=CI(1)
        FY=CI(2)
        FZ=CI(3)
		CO(1)=MOMENTX
        CO(2)=MOMENTY
        CO(3)=MOMENTZ
        CALL MPI_ALLREDUCE(CO(1:3),CI(1:3),3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        MX=CI(1)
        MY=CI(2)
        MZ=CI(3)
	END IF
	IF (N.EQ.0) THEN
	INQUIRE (FILE='FORCE.dat',EXIST=HEREf)
		IF (HEREf) THEN
		OPEN(50+N,FILE='FORCE.dat',FORM='FORMATTED',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
		ELSE
	OPEN(50+N,FILE='FORCE.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
		END IF
	INQUIRE (FILE='MOMENT.dat',EXIST=HEREf)
		IF (HEREf) THEN
		OPEN(500+N,FILE='MOMENT.dat',FORM='FORMATTED',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
		ELSE
	OPEN(500+N,FILE='MOMENT.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
		END IF	
	IF(RFRAME.EQ.0)THEN	
	WRITE(50+N,'(I14,1X,E14.7,1X,E14.7,1X,E14.7)')it,T,CL,CD
	WRITE(500+N,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')it,T,FX,FY,FZ
	ELSE
	WRITE(50+N,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')it,T,FX,FY,FZ
	WRITE(500+N,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')it,T,MX,MY,MZ
	END IF
	CLOSE(50+N)
	CLOSE(500+N)
	END IF		
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
!$OMP END MASTER 
!$OMP BARRIER 
	
	
	
	
	
	

END SUBROUTINE COMPUTEFORCE


SUBROUTINE COMPUTEFORCE2d(N)
!> @brief
!> This subroutine computes the forces on the wall in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,J,KMAXE,gqi_points,nnd
CHARACTER(LEN=12)::PROC,RESTFILE,PROC3
REAL::DRAG,LIFT,CD,CL,RX,PX,EX,surface_temp,RTEMP
REAL::FORCEXFR,SSX,CDF,LIFTF,DRAGF,FRICTIONF,TAUYX,TAUZX,TAUZY,SSY,SSZ,CF,TAUXX,TAUYY,TAUZZ
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,nx,ny,angle1,angle2
 REAL,DIMENSION(2)::CI,CO
 logical::heref
 INTEGER::IM
 REAL::TSOLR,TSOLU,TSOLE,TSOLV,TSOLW,TSOLP,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA,1:NUMBEROFPOINTS2)::QPOINTS2D
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
real,dimension(1:4)::viscl,laml
FORCEX=zero; FORCEY=zero; FORCEZ=zero;  FORCEXFR=zero
 CD=zero
 CL=zero
 CI(:)=zero
 CO(:)=zero
 KMAXE=XMPIELRANK(N)
 
!$OMP BARRIER 
!$OMP DO  REDUCTION(+:FORCEX,FORCEY,FORCEZ)
DO I=1,kmaxe
		if (ielem(n,i)%interior.eq.1)then	
		    do j=1,ielem(n,i)%ifca
		      if (ielem(n,i)%ibounds(j).gt.0)then
			  if (ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4)then
			      nx=IELEM(N,I)%FACEANGLEX(j)
			      ny=IELEM(N,I)%FACEANGLEY(j)
			      
			      
			  SSX=zero; SSP=zero; SSY=zero; 
			  
				
					  gqi_points=qp_line_n
					   if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					  NND=2
				      do K=1,nnd
					VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO
					  
					  call  QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  end if
					  surface_temp=IELEM(N,I)%SURF(J)
					  
				  
							  
				  
				  do im=1,gqi_points
				  
				  if (itestcase.eq.4)then
				  if (ielem(n,i)%ggs.eq.1)then
				  
				  VORTET1(1:2,1:2) = ILOCAL_RECON3(I)%GRADS(1:2,1:2)
				  ux = Vortet1(1,1);uy = Vortet1(1,2)
				  vx = Vortet1(2,1);vy = Vortet1(2,2)
				
				  
				  else
				  
				  vortet1(1,1:2)=ILOCAL_RECON3(i)%ULEFTV(1:2,2,J,IM)
				  vortet1(2,1:2)=ILOCAL_RECON3(i)%ULEFTV(1:2,3,J,IM)
				  
				 ux = Vortet1(1,1);uy = Vortet1(1,2)
				  vx = Vortet1(2,1);vy = Vortet1(2,2)
				
				  
				  
				  end if
				  end if
				  
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				    call cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
				    px=leftv(4)
				    
				    
				    if (itestcase.eq.4)then
				    CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
				  
				  
				  TAUXX=2.0d0*ux
				  TAUYY=2.0d0*vy
				  TAUYX=(UY + VX)
				  
				  SSX=SSX-((VISCL(1)*((NY*TAUYX)))*WEQUA2D(im))
				  SSY=SSY-((VISCL(1)*((NX*TAUYX)))*WEQUA2D(im))
				  end if
				  ssp=ssp+(px*WEQUA2D(im))
				   end do
				   
				   
				  SSP=ssp-PRES	
				  
				  FORCEX=FORCEX+(((SSP)*(surface_temp)*NX))+((SSX)*surface_temp)
				  FORCEY=FORCEY+(((SSP)*(surface_temp)*NY))+((SSY)*surface_temp)
				  

					    
					    
			END IF
		      end if
		    end do
		end if

			
		
END DO					 
!$OMP END DO
	

		
	
	
!$OMP BARRIER 
!$OMP MASTER 
	FORCEX=FORCEX*VECTORX
	FORCEY=FORCEY*VECTORY
	
        RTEMP=((AOA/180.0d0)*PI)
	LIFTF=(FORCEy*COS(RTEMP))-(FORCEX*SIN(RTEMP))
	DRAGF=(FORCEX*COS(RTEMP))+(FORCEy*SIN(RTEMP))
	CL=(2.0D0*LIFTF)/((RRES)*(ufreestream**2))
	CD=(2.0D0*DRAGF)/((RRES)*(ufreestream**2))


	CO(1)=CL
	CO(2)=CD
	CALL MPI_ALLREDUCE(CO(1:2),CI(1:2),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
	CL=CI(1)
	CD=CI(2)
	
	IF (N.EQ.0) THEN
	INQUIRE (FILE='FORCE.dat',EXIST=HEREf)
		IF (HEREf) THEN
		OPEN(50+N,FILE='FORCE.dat',FORM='FORMATTED',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
		ELSE
	OPEN(50+N,FILE='FORCE.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
		END IF
	WRITE(50+N,'(I14,1X,E14.7,1X,E14.7,1X,E14.7)')it,T,CL,CD
	
	CLOSE(50+N)
	END IF		
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
!$OMP END MASTER 
!$OMP BARRIER 
	
	
	
	
	
	

END SUBROUTINE COMPUTEFORCE2d


SUBROUTINE CALCULATE_RESIDUAL(N)
!> @brief
!> This subroutine computes and writes the residual for 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
real::suml3,dum_resi

KMAXE=XMPIELRANK(N)

ALLRES(:)=ZERO


IF ((ITESTCASE.LE.4).AND.(TURBULENCE.NE.1))THEN
!$OMP BARRIER 
!$OMP DO  REDUCTION(+:ALLRES)
DO I=1,KMAXE
	if (dg.eq.1)then
    ALLRES(1:nof_Variables)=ALLRES(1:nof_Variables)+((rhs(i)%VALdg(1,1:nof_Variables)*ielem(n,i)%totvolume)**2)
    else
    ALLRES(1:nof_Variables)=ALLRES(1:nof_Variables)+((rhs(i)%VAL(1:nof_Variables)*ielem(n,i)%totvolume)**2)

    end if
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,5
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/totalvolume

END DO


DO I=1,5
IF (INITIALRES(I).LE.ALLRES(I))THEN
INITIALRES(I)=ALLRES(I)
END IF
ALLRES(I)=ALLRES(I)/INITIALRES(I)

END DO
!$OMP END MASTER




END IF

IF (TURBULENCE.EQ.1)THEN
!$OMP BARRIER 
!$OMP DO  REDUCTION(+:ALLRES)
DO I=1,KMAXE
    ALLRES(1:nof_Variables)=ALLRES(1:nof_Variables)+((rhs(i)%VAL(1:nof_Variables)*ielem(n,i)%totvolume)**2)
    ALLRES(6:5+TURBULENCEEQUATIONS)=ALLRES(6:5+TURBULENCEEQUATIONS)+((RHST(I)%VAL(1:TURBULENCEEQUATIONS)*ielem(n,i)%totvolume)**2)
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,7
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/totalvolume

END DO



DO I=1,7
IF (INITIALRES(I).LE.ALLRES(I))THEN
INITIALRES(I)=ALLRES(I)
END IF
ALLRES(I)=ALLRES(I)/INITIALRES(I)

END DO

if (turbulenceequations.eq.1) allres(7)=1.0d0
!$OMP END MASTER

END IF



!$OMP MASTER
IF (N.EQ.0)THEN
IF ((ITESTCASE.LE.4).AND.(TURBULENCE.NE.1))THEN

	    OPEN(67,FILE='residual.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	    WRITE(67,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IT,ALLRES(1),ALLRES(2),ALLRES(3),ALLRES(4),ALLRES(5)
	    CLOSE(67)

ELSE

	    OPEN(67,FILE='residual.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	    WRITE(67,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IT,ALLRES(1),ALLRES(2),ALLRES(3),ALLRES(4),ALLRES(5),ALLRES(6),ALLRES(7)
	    CLOSE(67)



END IF

END IF 


if ((ALLRES(1).lt.reslimit).AND.(ALLRES(2).lt.reslimit).AND.(ALLRES(3).lt.reslimit).AND.(ALLRES(4).lt.reslimit).AND.(ALLRES(5).lt.reslimit))then
 kill=1
 end if



!$OMP END MASTER

 









End Subroutine


SUBROUTINE CALCULATE_RESIDUAL2D(N)
!> @brief
!> This subroutine computes and writes the residual for 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
real::suml3,dum_resi
KMAXE=XMPIELRANK(N)

ALLRES(:)=ZERO


IF ((ITESTCASE.LE.4).AND.(TURBULENCE.NE.1))THEN
!$OMP BARRIER 
!$OMP DO  REDUCTION(+:ALLRES)
DO I=1,KMAXE

    if (dg.eq.1)then
    ALLRES(1:nof_Variables)=ALLRES(1:nof_Variables)+((rhs(i)%VALdg(1,1:nof_Variables)*ielem(n,i)%totvolume)**2)
    else
    ALLRES(1:nof_Variables)=ALLRES(1:nof_Variables)+((rhs(i)%VAL(1:nof_Variables)*ielem(n,i)%totvolume)**2)

    end if
END DO
!$OMP END DO



!$OMP BARRIER
!$OMP MASTER

DO I=1,4
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/totalvolume

END DO


DO I=1,4
IF (INITIALRES(I).LE.ALLRES(I))THEN
INITIALRES(I)=ALLRES(I)
END IF
ALLRES(I)=ALLRES(I)/INITIALRES(I)

END DO

!$OMP END MASTER




END IF

IF (TURBULENCE.EQ.1)THEN
!$OMP BARRIER 
!$OMP DO  REDUCTION(+:ALLRES)
DO I=1,KMAXE
    ALLRES(1:nof_Variables)=ALLRES(1:nof_Variables)+((rhs(i)%VAL(1:nof_Variables)*ielem(n,i)%totvolume)**2)
    ALLRES(5:4+TURBULENCEEQUATIONS)=ALLRES(5:4+TURBULENCEEQUATIONS)+((RHST(I)%VAL(1:TURBULENCEEQUATIONS)*ielem(n,i)%totvolume)**2)
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,nof_variables+turbulenceequations
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/totalvolume

END DO



DO I=1,nof_variables+turbulenceequations
IF (INITIALRES(I).LE.ALLRES(I))THEN
INITIALRES(I)=ALLRES(I)
END IF
ALLRES(I)=ALLRES(I)/INITIALRES(I)

END DO
!$OMP END MASTER

END IF



!$OMP MASTER
IF (N.EQ.0)THEN
IF ((ITESTCASE.LE.4).AND.(TURBULENCE.NE.1))THEN

	    OPEN(67,FILE='residual.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	    WRITE(67,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IT,ALLRES(1),ALLRES(2),ALLRES(3),ALLRES(4)
	    CLOSE(67)

ELSE

	    OPEN(67,FILE='residual.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	    WRITE(67,'(I14,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IT,ALLRES(1),ALLRES(2),ALLRES(3),ALLRES(4),ALLRES(5)
	    CLOSE(67)



END IF

END IF 


 if ((ALLRES(1).lt.reslimit).AND.(ALLRES(2).lt.reslimit).AND.(ALLRES(3).lt.reslimit).AND.(ALLRES(4).lt.reslimit))then
 kill=1
 end if

!$OMP END MASTER







End Subroutine

SUBROUTINE CALCULATE_ERROR(N)
!> @brief
!> This subroutine computes and writes the l2,linfinity or l1 norm
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,K,KMAXE,ind_er
	REAL::EXACT,DUMMYOUT,DUMMYIN
	REAL::APROXIMATE
	real,dimension(15)::condm
	KMAXE=XMPIELRANK(N)
	L0NORM=ZERO;STENNORM=ZERO;L1NORM=ZERO
	ind_er=1
	if (multispecies.eq.1)then
	ind_er=nof_variables
	end if
	  !$OMP MASTER
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	  !$OMP END MASTER
	  !$OMP BARRIER
	  
			!$OMP DO REDUCTION (+:L1NORM)
			DO I=1,KMAXE
				IF (ITESTCASE.Le.3)THEN
				
				EXACT=U_E(I)%VAL(1,ind_er)
				
				APROXIMATE=U_C(I)%VAL(1,ind_er)
				
! 					IF ((ABS(APROXIMATE-EXACT)).GT.L0NORM(N,1))THEN
! 					L0NORM(N,1)=ABS(APROXIMATE-EXACT)
! 					END IF
					L1NORM=L1NORM+((APROXIMATE-EXACT)**2)*ielem(n,i)%TOTVOLUME
				END IF
 			END DO
 			!$OMP END DO 
 			
 			if (initcond.eq.0)then
 			!$OMP DO REDUCTION (+:L0NORM)
			DO I=1,KMAXE
				IF (ITESTCASE.Le.3)THEN
! 				condm(2)=ilocal_recon3(i)%cond(2)
! 				L0NORM=ilocal_recon3(i)%cond(1)
! 					IF (maxval(condm).GT.L0NORM)THEN
					L0NORM=L0NORM+abs(ilocal_recon3(i)%cond(1))
! 					END IF
! 					L1NORM(N,1)=L1NORM(N,1)+((ABS(APROXIMATE-EXACT)))
				END IF
 			END DO
 			!$OMP END DO 
 			ELSE
 			!$OMP DO REDUCTION (MAX:L0NORM)
			DO I=1,KMAXE
				IF (ITESTCASE.Le.3)THEN
				EXACT=U_E(I)%VAL(1,ind_er)
				
				
				APROXIMATE=U_C(I)%VAL(1,ind_er)
				
					IF ((ABS(APROXIMATE-EXACT)).GT.L0NORM)THEN
					L0NORM=ABS(APROXIMATE-EXACT)
					END IF
! 					L1NORM(N,1)=L1NORM(N,1)+((ABS(APROXIMATE-EXACT)))
				END IF
 			END DO
 			!$OMP END DO 
 			
 			
 			
 			END IF
 			IF (INITCOND.EQ.3)THEN
 			L0NORM=ZERO;L1NORM=TOLBIG
 			
 			!$OMP DO REDUCTION (MAX:L0NORM)
			DO I=1,KMAXE
					
					IF (U_C(I)%VAL(1,ind_er).GT.L0NORM)THEN
					L0NORM=U_C(I)%VAL(1,ind_er)
					END IF
! 					L1NORM(N,1)=L1NORM(N,1)+((ABS(APROXIMATE-EXACT)))
				
 			END DO
 			!$OMP END DO 
 			!$OMP DO REDUCTION (MIN:L1NORM)
			DO I=1,KMAXE
					
					IF (U_C(I)%VAL(1,ind_er).LT.L1NORM)THEN
					L1NORM=U_C(I)%VAL(1,ind_er)
					END IF
! 					L1NORM(N,1)=L1NORM(N,1)+((ABS(APROXIMATE-EXACT)))
				
 			END DO
 			!$OMP END DO 
 			
 			
 			END IF
 			
 			
 			if (initcond.eq.0)then
 			!$OMP DO REDUCTION (+:STENNORM)
 			DO I=1,KMAXE
				
				STENNORM=STENNORM+abs(ilocal_recon3(i)%cond(2))
 			END DO
 			!$OMP END DO 
 			
 			else
 			!$OMP DO REDUCTION (+:STENNORM)
 			DO I=1,KMAXE
				STENNORM=STENNORM+ielem(n,i)%STENCIL_DIST
 			END DO
 			!$OMP END DO 
 			end if
			
 			
 			
 			
 			!$OMP MASTER
 			IF (INITCOND.EQ.3)THEN
 			DUMMYOUT=L0NORM
 			CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
 			L0NORM=DUMMYIN
 			DUMMYOUT=L1NORM
 			CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
 			L1NORM=DUMMYIN
 			Else
 			
 			DUMMYOUT=L1NORM
 			CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
 			L1NORM=DUMMYIN
 			DUMMYOUT=L0NORM
 			IF (INITCOND.EQ.0)THEN
 			CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
 			DUMMYIN=DUMMYIN/IMAXE
 			ELSE
 			
 			CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
 			END IF
 			L0NORM=DUMMYIN
 			DUMMYOUT=STENNORM
 			CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
 			STENNORM=DUMMYIN
 			
 			
 			
 			
 			END IF
 			
 			
 			CPUX3(1) = MPI_Wtime()
 			if (n.eq.0)then
			OPEN(30,FILE='Errors.dat',FORM='FORMATTED',ACTION='write',position='append')
			if (initcond.eq.1)then
			WRITE(30,'(I9,1X,E14.7,1X,I4,1X,E14.7,1X,E14.7)')IMAXE,T,spatiladiscret,L0NORM,STENNORM/IMAXE
			
			else
			IF (INITCOND.NE.3)THEN
			
			WRITE(30,'(I9,1X,I4,1X,I4,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IMAXE,iorder,spatiladiscret,L0NORM,SQRT(L1NORM/TOTALVOLUME),STENNORM/IMAXE,(CPUX3(1)-CPUX2(1))*isize
			ELSE
			WRITE(30,'(I9,1X,I4,1X,I4,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IMAXE,iorder,spatiladiscret,L0NORM,L1NORM,STENNORM/IMAXE,(CPUX3(1)-CPUX2(1))*isize
			
			END IF
			end if
! 			WRITE(30,'(I9,1X,I4,1X,I4,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')IMAXE,iorder,spatiladiscret,L0NORM,SQRT(L1NORM/TOTALVOLUME),STENNORM/IMAXE,(CPUX3(1)-CPUX2(1))*isize
			close(30)
			end if
			!$OMP END MASTER
			!$OMP BARRIER

END SUBROUTINE CALCULATE_ERROR



SUBROUTINE OUTWRITEPARA3D
!> @brief
!> This subroutine writes the solution and the grid file in Ascii vtk format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 

INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml


 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(12))

KMAXE=XMPIELRANK(N)
nvar1=2





if (n.eq.0)then

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') 
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
end if



if (n.eq.0)then

OPEN(400+N,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW',ACTION='WRITE')
WRITE(400+N,'(A)')"# vtk DataFile Version 3.0"
WRITE(400+N,'(A)')"vtk output"
WRITE(400+N,'(A)')"ASCII"
WRITE(400+N,'(A)')"DATASET UNSTRUCTURED_GRID"
WRITE(400+N,'(A)')"FIELD FieldData 1"
WRITE(400+N,'(A)')"TIME 1 1 double"
WRITE(400+N,*) T
WRITE(400+N,'(A6,2X,I10,2X,A6)')"POINTS",IMAXN,"double"



allocate (Valuelocation(nvar1))


Valuelocation(:)=0
Valuelocation(1:2)=1

    
	if (binio.eq.0)then
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N,'(2x,g14.6,2x,g14.6,2x,g14.6)')x,y,z
	END DO
	CLOSE(96)
	else
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N,'(2x,g14.6,2x,g14.6,2x,g14.6)')x,y,z
	END DO
	CLOSE(96)
	end if

   
                    
                    
 WRITE(400+N,*)                   
write(400+N,'(A5,2X,I10,2X,I10)')"CELLS",IMAXE,(IMAXE*8)+IMAXE

if (binio.eq.0)then
OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97,*)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N,'(9i12)')8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)

else
OPEN(97,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N,'(9i12)')8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)
end if
                    
                    
WRITE(400+N,*)                    
WRITE(400+N,'(A10,2X,I10)')"CELL_TYPES",IMAXE
DO I=1,IMAXE
WRITE(400+N,*)12
END DO
WRITE(400+N,*)
WRITE(400+N,'(A9,2X,I10)')"CELL_DATA",IMAXE
               
                    
                    

 
  allocate(xbin(imaxe),XBIN2(IMAXE))
	

 END IF
 allocate(valuess(kmaxe))
 
  call mpi_barrier(MPI_COMM_WORLD,IERROR)
  
    
   
do j=1,NOF_VARIABLES
     
     if ((J.ge.2).and.(J.le.4))then
     
        DO I=1,KMAXE
        VALUESS(i)=U_C(I)%VAL(1,j)/U_C(I)%VAL(1,1)!0.0
        END DO
    END IF
    if (j.eq.1)then
        DO I=1,KMAXE
        VALUESS(i)=U_C(I)%VAL(1,j)
        END DO
    end if
    if (j.eq.5)then
                DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		END DO
    
    end if
    
		
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)


    IF (N.EQ.0)THEN

    IF (J.EQ.1)THEN
    WRITE(400+N,'(A)')"SCALARS  R double 1"
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N,*)
    WRITE(400+N,'(A)')"SCALARS  U double 1"
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N,*)
    WRITE(400+N,'(A)')"SCALARS  V double 1"
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N,*)
    WRITE(400+N,'(A)')"SCALARS  W double 1"
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N,*)
    WRITE(400+N,'(A)')"SCALARS  P  double 1"
    END IF
    WRITE(400+N,'(A)')"LOOKUP_TABLE default"
    
    
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
    
     WRITE(400+N,*)xbin(1:imaxe)
    END IF

    
end do

    
    
    
    

     
    
    
    
 IF (N.EQ.0)THEN
  CLOSE(400+N) 
   
   DEALLOCATE(XBIN,VALUESA,XBIN2,VALUELOCATION,ICELLA)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITEPARA3D

SUBROUTINE OUTWRITEPARA3Db
!> @brief
!> This subroutine writes the solution and the grid file in binary vtk format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(12))
nvar1=2
KMAXE=XMPIELRANK(N)




if (n.eq.0)then

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') 
	!proc4=".plt"
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
end if

lf=char(10)

if (n.eq.0)then

OPEN(400+N,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') IMAXN
WRITE(400+N)"POINTS "//str1//" double"//lf

allocate (Valuelocation(nvar1))


Valuelocation(:)=0
Valuelocation(1:2)=1

    
     
! 	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
!         DO I=1,IMAXN
! 	READ(96,*)j,x,y,z
! 	xbin(i)=x/scaler
! 	ybin(i)=y/scaler
!  	zbin(i)=z/scaler
! 	END DO
! 
!     CLOSE(96)
    
    if (binio.eq.0)then
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N)x,y,z
	END DO
	CLOSE(96)
	else
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N)x,y,z
	END DO
	CLOSE(96)
	end if

    
    
   
    
                    
WRITE(str1(1:15),'(i15)') IMAXE
WRITE(str2(1:15),'(i15)') (IMAXE*8)+IMAXE
write(400+N)"CELLS",str1//str2//lf
if (binio.eq.0)then
OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97,*)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)
ELSE
OPEN(97,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)

END IF
                    
                    
WRITE(str1(1:15),'(i15)') IMAXE                    
WRITE(400+N)"CELL_TYPES"//str1//lf	
DO I=1,IMAXE
WRITE(400+N)12
END DO

WRITE(400+N)"CELL_DATA"//str1//lf
!                
!                     
!                     
! 
 
  allocate(xbin(imaxe),XBIN2(IMAXE))
	
! 
 END IF
! 
!  
  
  ALLOCATE(VALUESS(KMAXE))
 
 call mpi_barrier(MPI_COMM_WORLD,IERROR)
!     
   
do j=1,NOF_VARIABLES
     
         DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(j)
		END DO
    
     
		
    
    
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    
    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  R double 1"//lf
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U double 1"//lf
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V double 1"//lf
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W double 1"//lf
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  P  double 1"//lf
    END IF
    IF (J.EQ.6)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  species1  double 1"//lf
    END IF
    IF (J.EQ.7)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  species2  double 1"//lf
    END IF
    IF (J.EQ.8)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  volumef1  double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
    
    
     WRITE(400+N)xbin(1:imaxe)
    END IF

     
    
    
end do
    
    if (ITESTCASE.EQ.4)THEN

	 DO I=1,KMAXE
		  
		  VALUESS(i)=IELEM(N,I)%VORTEX(1)
		END DO
    
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

	IF (N.EQ.0)THEN
    
    	WRITE(400+N)
    	WRITE(400+N)"SCALARS  Q double 1"//lf
	    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
     WRITE(400+N)xbin(1:imaxe)
    END IF
   END IF
   
   
   if (turbulence.EQ.1)THEN

	 DO I=1,KMAXE
		  
		  VALUESS(i)=U_CT(I)%VAL(1,1)
		END DO
    
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

	IF (N.EQ.0)THEN
    
    	WRITE(400+N)
    	WRITE(400+N)"SCALARS  NUT double 1"//lf
	    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
     WRITE(400+N)xbin(1:imaxe)
    END IF
   END IF
    
    

     
    
  
    
  IF (N.EQ.0)THEN
  CLOSE(400+N) 
   
   DEALLOCATE(XBIN,VALUELOCATION,XBIN2)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITEPARA3Db




SUBROUTINE movie_PARA
!> @brief
!> This subroutine writes the solution and the grid file in binary vtk format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

!
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2

      Integer::   Debug,III,NPts,NElm


      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(12))
nvar1=2
KMAXE=XMPIELRANK(N)




if (n.eq.0)then

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)')
	!proc4=".plt"
	OUTFILE="MOV_"//TRIM(ADJUSTL(PROC3))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
end if

lf=char(10)

if (n.eq.0)then

OPEN(400+N,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') IMAXN
WRITE(400+N)"POINTS "//str1//" double"//lf

allocate (Valuelocation(nvar1))


Valuelocation(:)=0
Valuelocation(1:2)=1



! 	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
!         DO I=1,IMAXN
! 	READ(96,*)j,x,y,z
! 	xbin(i)=x/scaler
! 	ybin(i)=y/scaler
!  	zbin(i)=z/scaler
! 	END DO
!
!     CLOSE(96)

    if (binio.eq.0)then
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N)x,y,z
	END DO
	CLOSE(96)
	else
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N)x,y,z
	END DO
	CLOSE(96)
	end if






WRITE(str1(1:15),'(i15)') IMAXE
WRITE(str2(1:15),'(i15)') (IMAXE*8)+IMAXE
write(400+N)"CELLS",str1//str2//lf
if (binio.eq.0)then
OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97,*)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)
ELSE
OPEN(97,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)

END IF


WRITE(str1(1:15),'(i15)') IMAXE
WRITE(400+N)"CELL_TYPES"//str1//lf
DO I=1,IMAXE
WRITE(400+N)12
END DO

WRITE(400+N)"CELL_DATA"//str1//lf
!
!
!
!

  allocate(xbin(imaxe),XBIN2(IMAXE))

!
 END IF
!
!

  ALLOCATE(VALUESS(KMAXE))

 call mpi_barrier(MPI_COMM_WORLD,IERROR)
!

do j=1,1

         DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)

		  valuess(i)=sqrt(leftv(2)**2+leftv(3)**2+leftv(4)**2)
		END DO





    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN

    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  vel double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf


		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do


     WRITE(400+N)xbin(1:imaxe)
    END IF




end do

    if (ITESTCASE.EQ.4)THEN

	 DO I=1,KMAXE

		  VALUESS(i)=IELEM(N,I)%VORTEX(1)
		END DO

    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

	IF (N.EQ.0)THEN

    	WRITE(400+N)
    	WRITE(400+N)"SCALARS  Q double 1"//lf
	    WRITE(400+N)"LOOKUP_TABLE default"//lf

		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
     WRITE(400+N)xbin(1:imaxe)
    END IF
   END IF


   if (turbulence.EQ.1)THEN

	 DO I=1,KMAXE

		  VALUESS(i)=U_CT(I)%VAL(1,1)
		END DO

    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

	IF (N.EQ.0)THEN

    	WRITE(400+N)
    	WRITE(400+N)"SCALARS  NUT double 1"//lf
	    WRITE(400+N)"LOOKUP_TABLE default"//lf

		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
     WRITE(400+N)xbin(1:imaxe)
    END IF
   END IF







  IF (N.EQ.0)THEN
  CLOSE(400+N)

   DEALLOCATE(XBIN,VALUELOCATION,XBIN2)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)


  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


 deallocate(variables)





END SUBROUTINE movie_PARA


 !added by holger foysi
                          !!-----------------------------------------------------------------------------------
                          SUBROUTINE OUTWRITEPARA2Db
                          !!-----------------------------------------------------------------------------------
                            !!-----------------------------------------------------------------------------------
                            !> @brief
                            !> This subroutine writes the 2D solution including the grid in binary vtk format
                            !!-----------------------------------------------------------------------------------
                            use ISO_C_BINDING
                            IMPLICIT NONE
                            real,dimension(1:nof_Variables)::leftv
							real::MP_PINFL,gammal
							real,dimension(1:nof_Variables)::RIGHTv
							real::MP_PINFR,gammaR
							real::angle1,angle2,nx,ny,nz
							real,dimension(1:4)::viscl,laml
                            INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
                            REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
                            REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
                            INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
                            REAL,allocatable,DIMENSION(:)::VARIABLES
                            REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
                            INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
                            LOGICAL::HEREV
                            REAL,DIMENSION(5)::TOTAL
                            CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
                            integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
                            real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
                            real,allocatable,dimension(:,:)::FBIN
                            integer,allocatable,dimension(:,:)::icon
                            INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
                            real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
                            character(LEN=:),allocatable::out1
                            character*1 NULCHAR
                            CHARACTER(LEN=1)   :: flui,lf
                            CHARACTER(LEN=15)  :: str1,str2,str1a,str2a
                            Integer::   Debug,III,NPts,NElm
                            Real::    SolTime
                            Integer:: VIsDouble, FileType, count
                            Integer:: ZoneType,StrandID,ParentZn,IsBlock,strl
                            Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
                            POINTER   (NullPtr,Null)
                            Integer:: Null(*)
                            character(len=:),allocatable  :: str_imaxn,str_imaxe,str_imaxe2,str_time
                            character(1) :: c
                            allocate(variables(12))
                            nvar1=2
                            KMAXE=XMPIELRANK(N)
                            
                            !the previous use of field didn't work, so time is added to the title

                            !create filename
                            if (n.eq.0)then
                               WRITE(PROC3,FMT='(I10)') IT
                               WRITE(PROC5,FMT='(I10)') 
                               OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//".vtk"
                               ITGFD=len_trim(OUTFILE)
                               allocate(character(LEN=itgfd) ::out1)
                               out1=OUTFILE(1:itgfd)
                            end if

                            ! char(10)  = linefeed ! Alternativeley use ",new_line(c)"
                            lf=char(10)

                            if (n.eq.0)then
                               !file
                               OPEN(400+N,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
                               !header
                               WRITE(400+N)"# vtk DataFile Version 3.0",new_line(c)
                               !create dynamic string with size equal to the number of digits of Time
                               WRITE(str1,'(f12.6)') T
                               str_time=trim(adjustl(str1))
                               WRITE(400+N)"UCNS3D vtk output 2D. Time: "//str_time,new_line(c)
                               WRITE(400+N)"BINARY",new_line(c)
                               WRITE(400+N)"DATASET UNSTRUCTURED_GRID",new_line(c)

                               !create dynamic string with size equal to the number of digits of the integer
                               WRITE(str1,'(i0)') IMAXN
                               str_imaxn=trim(adjustl(str1))

                               !Number of node points
                               WRITE(400+N)"POINTS "//str_imaxn//" double",new_line(c)

                               if (binio.eq.0)then !ascii
                                  OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
                                  DO I=1,IMAXN
                                     READ(96,*)j,x,y
                                     x=x/scaler;y=y/scaler
                                     write(400+N) x,y,0.d0  !for 2d vtk needs x,y,z, too, just with z=0
                                  END DO
                                  CLOSE(96)
                               else !binary
                                  OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                  DO I=1,IMAXN
                                     READ(96)j,x,y
                                     x=x/scaler;y=y/scaler
                                     write(400+N) x,y,0.d0  !for 2d vtk needs x,y,z, too, just with z=0
                                  END DO
                                  CLOSE(96)
                               end if

                               !create dynamic strings with size equal to the number of digits of the integer
                               !# of Cells
                               WRITE(str1,'(i0)') IMAXE
                               str_imaxe  = trim(adjustl(str1))

                               !cell list size:the total number of integer values required to represent the list
                               !changes needed if grids with different cell types are used, here every cell has the same type
                               WRITE(str2,'(i0)') IMAXE*5 ! 
                               str_imaxe2 = trim(adjustl(str2))        
                               
                               write(400+N) "CELLS "//str_imaxe//" "//str_imaxe2,new_line(c)
                               !the grid is such, that triangles are defined using quads with two points being equal (J3=J4)
                               if (binio.eq.0)then
                                  OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
                                  DO I=1,IMAXE
                                     READ(97,*) j,J1,J2,J3,J4
                                     write(400+N) 4,J1-1,J2-1,J3-1,J4-1
                                  END DO
                                  CLOSE(97)
                               ELSE
                                  OPEN(97,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
                                  DO I=1,IMAXE
                                     READ(97) j,J1,J2,J3,J4
                                     write(400+N) 4,J1-1,J2-1,J3-1,J4-1
                                  END DO
                                  CLOSE(97)
                               END IF

                               !Specify type of cells for each cell
                               !changes needed if grids with different cell types are used, here every cell has the same type
                               !hexahedron in 3D (12), quad in 2D (9), triangles (5) etc., see
                               !https://lorensen.github.io/VTKExamples/site/VTKFileFormats/
                               WRITE(400+N)"CELL_TYPES "//str_imaxe,new_line(c)
                               DO I=1,IMAXE
                                  WRITE(400+N) 9  
                               END DO

                               ! the data is cell centered. Use CellDataToPointData in Paraview for Warp by Scalar
                               WRITE(400+N)"CELL_DATA "//str_imaxe,new_line(c)

                               allocate(xbin(imaxe),XBIN2(IMAXE))
                               ! 
                            END IF

                            ALLOCATE(VALUESS(KMAXE))
                            call mpi_barrier(MPI_COMM_WORLD,IERROR)

                            do j=1,NOF_VARIABLES
                                
                              
                                  DO I=1,KMAXE
                                        leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
                                        call cons2prim(N,leftv,MP_PINFl,gammal)
                                        VALUESS(i)=leftv(j)
                                  END DO

                               call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset, &
                                                mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

                               IF (N.EQ.0)THEN
                                  IF (J.EQ.1)THEN
                                     WRITE(400+N)"SCALARS R double 1",new_line(c)
                                  END IF
                                  IF (J.EQ.2)THEN
                                     WRITE(400+N)"SCALARS U double 1",new_line(c)
                                  END IF
                                  IF (J.EQ.3)THEN
                                     WRITE(400+N)"SCALARS V double 1",new_line(c)
                                  END IF
                                  IF (J.EQ.4)THEN
                                     WRITE(400+N)"SCALARS P double 1",new_line(c)
                                  END IF
                                  IF (J.EQ.5)THEN
                                     WRITE(400+N)"SCALARS species1 double 1",new_line(c)
                                  END IF
                                  IF (J.EQ.6)THEN
                                     WRITE(400+N)"SCALARS species2 double 1",new_line(c)
                                  END IF
                                  IF (J.EQ.7)THEN
                                     WRITE(400+N)"SCALARS volumef double 1",new_line(c)
                                  END IF
                                  WRITE(400+N)"LOOKUP_TABLE default",new_line(c)
                                  do i=1,imaxe
                                     xbin(XMPI_RE(i))=xbin2(i)
                                  end do
                                  do i=1,imaxe
                                     WRITE(400+N) xbin(i)  !test, not clear of whether cell or node based
                                  end do
                               END IF
                            end do

                            
                            if (turbulence.EQ.1)THEN

                            DO I=1,KMAXE
                                
                                VALUESS(i)=U_CT(I)%VAL(1,1)
                                END DO
                            
                            call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

                            IF (N.EQ.0)THEN
                            
                                WRITE(400+N)
                                WRITE(400+N)"SCALARS  NUT double 1"//lf
                                WRITE(400+N)"LOOKUP_TABLE default"//lf
                            
                                do i=1,imaxe
                                xbin(XMPI_RE(I))=xbin2(I)
                                end do
                            WRITE(400+N)xbin(1:imaxe)
                            END IF
                        END IF
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            IF (N.EQ.0)THEN
                               ! close file
                               CLOSE(400+N) 

                               DEALLOCATE(XBIN,XBIN2)
                               deallocate(out1)
                            END IF
                            DEALLOCATE (VALUESS)

                            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
                            deallocate(variables)

    END SUBROUTINE OUTWRITEPARA2Db




SUBROUTINE OUTWRITEPARA3DbP
!> @brief
!> This subroutine writes the solution and the grid file in binary vtk format
use ISO_C_BINDING
IMPLICIT NONE
INTEGER::KMAXE
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,OUTFILE2,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
KMAXE=XMPIELRANK(N)

allocate(xbin(kmaxe))
WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') N 
	!proc4=".plt
	OUTFILE="OUT_"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)


lf=char(10)


OPEN(400+n,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') kmaxn
WRITE(400+N)"POINTS "//str1//" double"//lf


    
    DO I=1,KMAXN
	write(400+N)INODER4(i)%CORD(1),INODER4(i)%CORD(2),INODER4(i)%CORD(3)
	END DO
   
                    
WRITE(str1(1:15),'(i15)') KMAXE
WRITE(str2(1:15),'(i15)') (KMAXE*8)+KMAXE
write(400+N)"CELLS",str1//str2//lf
DO I=1,kMAXE
write(400+N)8,el_connect(I,1)-1,el_connect(I,2)-1,el_connect(I,3)-1,el_connect(I,4)-1,el_connect(I,5)-1,el_connect(I,6)-1,el_connect(I,7)-1,el_connect(I,8)-1
END DO
                   
                    
WRITE(str1(1:15),'(i15)') KMAXE                    
WRITE(400+N)"CELL_TYPES"//str1//lf	
DO I=1,KMAXE
WRITE(400+N)12

END DO

WRITE(400+N)"CELL_DATA"//str1//lf
!     
   
do j=1,NOF_VARIABLES


     DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  XBIN(i)=leftv(j)
		END DO



    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  R double 1"//lf
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U double 1"//lf
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V double 1"//lf
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W double 1"//lf
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  P  double 1"//lf
    END IF
    IF (J.EQ.6)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  species1  double 1"//lf
    END IF
    IF (J.EQ.7)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  species2 double 1"//lf
    END IF
    IF (J.EQ.8)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  volumef  double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    
		
    
    
     WRITE(400+N)xbin(1:Kmaxe)
     
    

     
    
    
end do

    if (ITESTCASE.EQ.4)THEN

	 DO I=1,KMAXE
		  
		  XBIN(i)=IELEM(N,I)%VORTEX(1)
		END DO

        WRITE(400+N)
    	WRITE(400+N)"SCALARS  Q double 1"//lf
	    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    WRITE(400+N)xbin(1:imaxe)
    END IF
    
    if (TURBULENCE.EQ.1)THEN

	 DO I=1,KMAXE
		  
		  XBIN(i)=U_CT(I)%VAL(1,1)
		END DO

        WRITE(400+N)
    	WRITE(400+N)"SCALARS  NUT double 1"//lf
	    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    WRITE(400+N)xbin(1:imaxe)
    END IF
    
    
    
    close(400+n)
    DEALLOCATE(XBIN,OUT1)
    
    
  
    
  
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        


	
	

END SUBROUTINE OUTWRITEPARA3DbP




SUBROUTINE OUTWRITEPARA3DbPav
!> @brief
!> This subroutine writes the solution and the grid file in binary vtk format
use ISO_C_BINDING
IMPLICIT NONE
INTEGER::KMAXE,ind1
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
KMAXE=XMPIELRANK(N)


if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if


allocate(xbin(kmaxe))
WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)')N 
	!proc4=".plt
	OUTFILE="OUT_AV"//TRIM(ADJUSTL(PROC3))//"_"//TRIM(ADJUSTL(PROC5))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)


lf=char(10)


OPEN(400+n,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') kmaxn
WRITE(400+N)"POINTS "//str1//" double"//lf


    
    DO I=1,KMAXN
	write(400+N)INODER4(i)%CORD(1),INODER4(i)%CORD(2),INODER4(i)%CORD(3)
	END DO
   
                    
WRITE(str1(1:15),'(i15)') KMAXE
WRITE(str2(1:15),'(i15)') (KMAXE*8)+KMAXE
write(400+N)"CELLS",str1//str2//lf
DO I=1,kMAXE
write(400+N)8,el_connect(I,1)-1,el_connect(I,2)-1,el_connect(I,3)-1,el_connect(I,4)-1,el_connect(I,5)-1,el_connect(I,6)-1,el_connect(I,7)-1,el_connect(I,8)-1
END DO
                   
                    
WRITE(str1(1:15),'(i15)') KMAXE                    
WRITE(400+N)"CELL_TYPES"//str1//lf	
DO I=1,KMAXE
WRITE(400+N)12
END DO

WRITE(400+N)"CELL_DATA"//str1//lf
!     
   
do j=1,11

    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  R_mean double 1"//lf
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U_mean double 1"//lf
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V_mean double 1"//lf
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W_mean double 1"//lf
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  P_mean  double 1"//lf
    END IF
    IF (J.EQ.6)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U_rms  double 1"//lf
    END IF
    IF (J.EQ.7)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V_rms double 1"//lf
    END IF
    IF (J.EQ.8)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W_rms  double 1"//lf
    END IF
     IF (J.EQ.9)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  UV  double 1"//lf
    END IF
     IF (J.EQ.10)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  UW  double 1"//lf
    END IF
     IF (J.EQ.11)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  WV  double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
        
        if (j.eq.1)then
		do i=1,Kmaxe
		xbin(I)=U_C(I)%VAL(ind1,J)
		end do
		end if
		if ((j.gt.1).and.(j.lt.5))then
		do i=1,Kmaxe
		xbin(I)=U_C(I)%VAL(ind1,J)/U_C(I)%VAL(ind1,1)
		end do
		end if
		if (j.eq.5)then
		do i=1,Kmaxe
		leftv(1:nof_variables)=U_C(I)%VAL(ind1,1:nof_variables)
		call CONS2PRIM(N,leftv,MP_PINFl,gammal)
		xbin(I)=leftv(5)
		end do
		end if
		if (j.ge.6)then
		do i=1,Kmaxe
		xbin(I)=U_C(I)%RMS(j-6)
		end do
		end if
    
    
     WRITE(400+N)xbin(1:Kmaxe)
   

     
    
    
end do
    
    
    
    close(400+n)
        DEALLOCATE(XBIN,OUT1)
     
    
  
    
  
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        


	
	

END SUBROUTINE OUTWRITEPARA3DbPav





SUBROUTINE OUTWRITEPARA3Dsb
!> @brief
!> This subroutine writes the solution and the surface file in binary vtk format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(12))
nvar1=2





if (n.eq.0)then

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') 
	!proc4=".plt"
	OUTFILE="SURF_"//TRIM(ADJUSTL(PROC3))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
end if

lf=char(10)



if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)then
OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
 close(96)
 else
 OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
 close(96)
 
 
 end if
 
 
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
! OPEN(400+N,FILE=OUTFILE,FORM='UNFORMATTED',STATUS='NEW',ACTION='WRITE')
 OPEN(400+N,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') igf2
WRITE(400+N)"POINTS "//str1//" double"//lf

allocate (Valuelocation(nvar1))

    IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+N)x/scaler,y/scaler,Z/scaler
	end if
	END DO

    CLOSE(96)
    ELSE
    OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+N)x/scaler,y/scaler,Z/scaler
	end if
	END DO

    CLOSE(96)
    
    
    END IF

    
   
    
    
     IF (BINIO.EQ.0)THEN
     OPEN(98,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
     END IF
      IF (BINIO.EQ.1)THEN
     OPEN(98,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
     END IF
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		IF (BINIO.EQ.0)THEN
		DO K=1,iMAXB
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		ELSE
		DO K=1,iMAXB
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		
		END IF
	    
 		close(98)
 		
 		
 		
 		
 		WRITE(str1(1:15),'(i15)') igf2
 WRITE(str2(1:15),'(i15)') (igf2*4)+igf2
 write(400+N)"CELLS",str1//str2//lf
 		DO I=1,IGF2
 		WRITE(400+N)4,ICON(1,I)-1,ICON(2,I)-1,ICON(3,I)-1,ICON(4,I)-1
 		END DO
 		
	
		deallocate(icon)	
		deallocate(inog)


                    
! WRITE(str1(1:15),'(i15)') IMAXE
! WRITE(str2(1:15),'(i15)') (IMAXE*8)+IMAXE
! write(400+N)"CELLS",str1//str2//lf
! 
! OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
! DO I=1,IMAXE
! READ(97,*)j,J1,J2,J3,J4,j5,j6,j7,j8
! write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
! END DO
! CLOSE(97)
                    
                    
                 
 WRITE(400+N)"CELL_TYPES"//str1//lf	
 DO I=1,igf2 
 WRITE(400+N)9
 END DO
! 
WRITE(400+N)"CELL_DATA"//str1//lf
!                
!                     
!                     
! 
allocate(xbin(totwalls),XBIN2(TOTWALLS))

 END IF

  totiw=xmpiwall(n)
  if (xmpiwall(n).gt.0)then
  ALLOCATE(VALUESS(xmpiwall(n)))
  end if
!     
   
do j=1,NOF_VARIABLES
     
     if ((J.ge.2).and.(J.le.4))then
     
         IF (TOTIW.GT.0)THEN
		      	
		      do I=1,TOTIW
                valuess(i)=U_C(IBOUND_T(I))%VAL(1,j)/U_C(IBOUND_T(I))%VAL(1,1)
                end do
				
		      end if
    END IF
    if (j.eq.1)then
         IF (TOTIW.GT.0)THEN
		      	
		      do I=1,TOTIW
                valuess(i)=U_C(IBOUND_T(I))%VAL(1,j)
                end do
				
		      end if
    end if
    if (j.eq.5)then
                IF (TOTIW.GT.0)THEN
		     DO I=1,TOTIW
						leftv(1:nof_Variables)=U_C(IBOUND_T(I))%VAL(1,1:nof_Variables)
						CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
						VALUESS(i)=leftv(5)
						
					  END DO
		      end if
    
    end if
    
		
    
    call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
!     call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    
    
    
    
    
    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  R double 1"//lf
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U double 1"//lf
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V double 1"//lf
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W double 1"//lf
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  P  double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    
    
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
		
    
    
     WRITE(400+N)XBIN(1:TOTWALLS)
    END IF

     
    
    
end do
    
    
    
    

     
    
 
    
  IF (N.EQ.0)THEN
  
   
   DEALLOCATE(XBIN,XBIN2,VALUELOCATION)
  deallocate(out1)
  CLOSE(400+N)  
  END IF
   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS)
  end if
  

          
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITEPARA3Dsb






SUBROUTINE OUTWRITEPARA3Dbav
!> @brief
!> This subroutine writes the averaged solution and the grid file in binary vtk format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,j1,j2,j3,j4,j5,j6,j7,j8,IND1
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(12))
nvar1=2
KMAXE=XMPIELRANK(N)

if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if


if (n.eq.0)then

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') 
	!proc4=".plt"
	OUTFILE="OUT_AV"//TRIM(ADJUSTL(PROC3))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
end if

lf=char(10)

if (n.eq.0)then

OPEN(400+N,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') IMAXN
WRITE(400+N)"POINTS "//str1//" double"//lf

allocate (Valuelocation(nvar1))


Valuelocation(:)=0
Valuelocation(1:2)=1

    
     
! 	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
!         DO I=1,IMAXN
! 	READ(96,*)j,x,y,z
! 	xbin(i)=x/scaler
! 	ybin(i)=y/scaler
!  	zbin(i)=z/scaler
! 	END DO
! 
!     CLOSE(96)
    
    if (binio.eq.0)then
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N)x,y,z
	END DO
	CLOSE(96)
	else
	OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
        DO I=1,IMAXN
	READ(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+N)x,y,z
	END DO
	CLOSE(96)
	end if

    
    
   
    
                    
WRITE(str1(1:15),'(i15)') IMAXE
WRITE(str2(1:15),'(i15)') (IMAXE*8)+IMAXE
write(400+N)"CELLS",str1//str2//lf
if (binio.eq.0)then
OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97,*)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)
ELSE
OPEN(97,FILE='GRID.cel',FORM='UNFORMATTED',STATUS='old',ACTION='read')
DO I=1,IMAXE
READ(97)j,J1,J2,J3,J4,j5,j6,j7,j8
write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
END DO
CLOSE(97)

END IF
                    
                    
WRITE(str1(1:15),'(i15)') IMAXE                    
WRITE(400+N)"CELL_TYPES"//str1//lf	
DO I=1,IMAXE
WRITE(400+N)12
END DO

WRITE(400+N)"CELL_DATA"//str1//lf
!                
!                     
!                     
! 
 
  allocate(xbin(imaxe),XBIN2(IMAXE))
	
! 
 END IF
! 
!  
  
  ALLOCATE(VALUESS(KMAXE))
 
 call mpi_barrier(MPI_COMM_WORLD,IERROR)
!     
   
do j=1,NOF_VARIABLES+6
     
     if ((J.ge.2).and.(J.le.4))then
     
        DO I=1,KMAXE
        VALUESS(i)=U_C(I)%VAL(IND1,j)/U_C(I)%VAL(IND1,1)!0.0
        END DO
    END IF
    if (j.eq.1)then
        DO I=1,KMAXE
        VALUESS(i)=U_C(I)%VAL(IND1,j)
        END DO
    end if
    if (j.eq.5)then
                DO I=1,KMAXE
		  leftv(1:nof_Variables)=U_C(I)%VAL(IND1,1:nof_Variables)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  VALUESS(i)=leftv(5)
		END DO
    
    end if
    IF (J.GT.5)THEN
    DO I=1,KMAXE
	VALUESS(i)=U_C(I)%RMS(J-5)
    END DO
    END if
    
		
    
    
    call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    
    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  R_MEAN double 1"//lf
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U_MEAN double 1"//lf
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V_MEAN double 1"//lf
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W_MEAN double 1"//lf
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  P_MEAN  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U_RMS  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V_RMS  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W_RMS  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  UV  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  UW  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  WV  double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
    
     WRITE(400+N)xbin(1:imaxe)
    END IF
    
     
    
    
end do
    
    
    
    

     
    
  
    
  IF (N.EQ.0)THEN
  CLOSE(400+N) 
   
   DEALLOCATE(XBIN,XBIN2,VALUELOCATION)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITEPARA3Dbav

SUBROUTINE OUTWRITEPARA3Dsbav
!> @brief
!> This subroutine writes the averaged solution and the surface file in Ascii vtk format
use ISO_C_BINDING
IMPLICIT NONE

! EXTERNAL TecIni112
! EXTERNAL TecZne112
! EXTERNAL TECDAT112
! EXTERNAL TECNODE112
! EXTERNAL  TECEND112

! 
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
INTEGER::KMAXE,KK,KFK,ICPUID,L,IHGT,IHGJ,kkd
REAL::X,Y,Z,DENOMINATOR,TUY,TVX,TWX,TUZ,TVZ,TWY,SNORM,ONORM
REAL,ALLOCATABLE,DIMENSION(:)::IFINT,TFINT,NDR,NDS
INTEGER::INEEDT,JJ,IX,IX1,I1,I2,I3,I4,I5,DECOMF,KD
REAL,allocatable,DIMENSION(:)::VARIABLES
REAL,DIMENSION(3,3)::AVORT,TVORT,SVORT,OVORT
INTEGER::INX,I,K,J,M,O,P,Q,JK,imax,jmax,kmax,igf,igf2,DUMG,DUML,IMAXP,nvar1,IND1,j1,j2,j3,j4,j5,j6,j7,j8,icount_wall
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=30)::PROC,OUTFILE,PROC3,SURFILE,proc4,proc5
integer::ierr,cv,TecIni112,TecZne112,TECDAT112,TECNODE112,TECEND112,ITGFD
real,allocatable,dimension(:)::xbin,ybin,zbin,XBIN2
real,allocatable,dimension(:,:)::FBIN
integer,allocatable,dimension(:,:)::icon
INTEGER,ALLOCATABLE,DIMENSION(:)::Valuelocation,inog,ICELL,ICELLA
real,ALLOCATABLE,DIMENSION(:)::valuess,VALUESA
character(LEN=:),allocatable::out1
character*1 NULCHAR
CHARACTER(LEN=1)   :: flui,lf
CHARACTER(LEN=15)  :: str1,str2
 
      Integer::   Debug,III,NPts,NElm

  
      Real::    SolTime
      Integer:: VIsDouble, FileType
      Integer:: ZoneType,StrandID,ParentZn,IsBlock
      Integer:: ICellMax,JCellMax,KCellMax,NFConns,FNMode,ShrConn
      POINTER   (NullPtr,Null)
      Integer:: Null(*)
allocate(variables(12))
NVAR1=2

if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if



if (n.eq.0)then

WRITE(PROC3,FMT='(I10)') IT
WRITE(PROC5,FMT='(I10)') 
	!proc4=".plt"
	OUTFILE="SURF_"//TRIM(ADJUSTL(PROC3))//".vtk"!//TRIM(ADJUSTL(PROC4))
	ITGFD=len_trim(OUTFILE)
	allocate(character(LEN=itgfd) ::out1)
	out1=OUTFILE(1:itgfd)
end if

lf=char(10)



if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)then
OPEN(96,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
 close(96)
 else
 OPEN(96,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
allocate (inog(imaxn))
inog(:)=0
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
 close(96)
 
 
 end if
 
 
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
! OPEN(400+N,FILE=OUTFILE,FORM='UNFORMATTED',STATUS='NEW',ACTION='WRITE')
 OPEN(400+N,FILE=OUTFILE,STATUS='REPLACE',ACCESS='STREAM',CONVERT='BIG_ENDIAN')
WRITE(400+N)"# vtk DataFile Version 3.0"//lf
WRITE(400+N)"vtk output"//lf
WRITE(400+N)"BINARY"//lf
WRITE(400+N)"DATASET UNSTRUCTURED_GRID"//lf
WRITE(400+N)"FIELD FieldData 1"//lf
WRITE(400+N)"TIME 1 1 double"//lf
WRITE(400+N) T
WRITE(str1(1:15),'(i15)') igf2
WRITE(400+N)"POINTS "//str1//" double"//lf

allocate (Valuelocation(nvar1))

    IF (BINIO.EQ.0)THEN
	OPEN(96,FILE='GRID.vrt',FORM='FORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+N)x/scaler,y/scaler,Z/scaler
	end if
	END DO

    CLOSE(96)
    ELSE
    OPEN(96,FILE='GRID.vrt',FORM='UNFORMATTED',STATUS='old',ACTION='read')
	igf2=0
        DO I=1,IMAXN
	READ(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+N)x/scaler,y/scaler,Z/scaler
	end if
	END DO

    CLOSE(96)
    
    
    END IF

    
   
    
    
     IF (BINIO.EQ.0)THEN
     OPEN(98,FILE='GRID.bnd',FORM='FORMATTED',STATUS='old',ACTION='read')
     END IF
      IF (BINIO.EQ.1)THEN
     OPEN(98,FILE='GRID.bnd',FORM='UNFORMATTED',STATUS='old',ACTION='read')
     END IF
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		IF (BINIO.EQ.0)THEN
		DO K=1,iMAXB
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		ELSE
		DO K=1,iMAXB
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		END DO
		
		END IF
	    
 		close(98)
 		
 		
 		
 		
 		WRITE(str1(1:15),'(i15)') igf2
 WRITE(str2(1:15),'(i15)') (igf2*4)+igf2
 write(400+N)"CELLS",str1//str2//lf
 		DO I=1,IGF2
 		WRITE(400+N)4,ICON(1,I)-1,ICON(2,I)-1,ICON(3,I)-1,ICON(4,I)-1
 		END DO
 		
	
		deallocate(icon)	
		deallocate(inog)


                    
! WRITE(str1(1:15),'(i15)') IMAXE
! WRITE(str2(1:15),'(i15)') (IMAXE*8)+IMAXE
! write(400+N)"CELLS",str1//str2//lf
! 
! OPEN(97,FILE='GRID.cel',FORM='FORMATTED',STATUS='old',ACTION='read')
! DO I=1,IMAXE
! READ(97,*)j,J1,J2,J3,J4,j5,j6,j7,j8
! write(400+N)8,J1-1,J2-1,J3-1,J4-1,J5-1,j6-1,j7-1,j8-1
! END DO
! CLOSE(97)
                    
                    
                 
 WRITE(400+N)"CELL_TYPES"//str1//lf	
 DO I=1,igf2 
 WRITE(400+N)9
 END DO
! 
WRITE(400+N)"CELL_DATA"//str1//lf
!                
!                     
!                     
! 
allocate(xbin(totwalls),xbin2(totwalls))

 END IF

  totiw=xmpiwall(n)
  if (xmpiwall(n).gt.0)then
  ALLOCATE(VALUESS(xmpiwall(n)))
  end if
!     
   
do j=1,NOF_VARIABLES+6
     
     if ((J.ge.2).and.(J.le.4))then
     
         IF (TOTIW.GT.0)THEN
		      	
		      do I=1,TOTIW
                valuess(i)=U_C(IBOUND_T(I))%VAL(1,j)/U_C(IBOUND_T(I))%VAL(IND1,1)
                end do
				
		      end if
    END IF
    if (j.eq.1)then
         IF (TOTIW.GT.0)THEN
		      	
		      do I=1,TOTIW
                valuess(i)=U_C(IBOUND_T(I))%VAL(IND1,j)
                end do
				
		      end if
    end if
    if (j.eq.5)then
                IF (TOTIW.GT.0)THEN
		     DO I=1,TOTIW
						leftv(1:nof_Variables)=U_C(IBOUND_T(I))%VAL(IND1,1:nof_Variables)
						CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
						VALUESS(i)=leftv(5)
						
					  END DO
		      end if
    
    end if
     IF (J.GT.5)THEN
      IF (TOTIW.GT.0)THEN
	DO I=1,TOTIW
	  VALUESS(i)=U_C(IBOUND_T(I))%RMS(J-5)
	end do
      END IF
    END if
		
    
    call MPI_GATHERv(valuess,xmpiwall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiwall,woffset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
!     call MPI_GATHER(VALUESS,IMAXP,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN
    
    IF (J.EQ.1)THEN
    WRITE(400+N)"SCALARS  R_MEAN double 1"//lf
    END IF
    IF (J.EQ.2)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U_MEAN double 1"//lf
    END IF
    IF (J.EQ.3)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V_MEAN double 1"//lf
    END IF
    IF (J.EQ.4)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W_MEAN double 1"//lf
    END IF
    IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  P_MEAN  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  U_RMS  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  V_RMS  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  W_RMS  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  UV  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  UW  double 1"//lf
    END IF
     IF (J.EQ.5)THEN
    WRITE(400+N)
    WRITE(400+N)"SCALARS  WV  double 1"//lf
    END IF
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    
		do i=1,totwalls
		xbin(XMPI_WRE(I))=xbin2(I)
		end do
    
     WRITE(400+N)xbin(1:totwalls)
    END IF

     
    
    
end do
    
    
    
    

     
    
 
    
  IF (N.EQ.0)THEN
  
   
   DEALLOCATE(XBIN,VALUELOCATION,XBIN2)
  deallocate(out1)
  CLOSE(400+N)  
  END IF
   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS)
  end if
  

          
        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITEPARA3Dsbav



SUBROUTINE FIX_NODES_LOCAL
IMPLICIT NONE
INTEGER::NDLC_COUNT1,NDLC_COUNT2,indexgt
INTEGER::KMAXE,ND_LC_NODES,I,J,K,L,M,COUNTFNODES,indfc
INTEGER,ALLOCATABLE,DIMENSION(:)::NDLC_ARRAY1
INTEGER,ALLOCATABLE,DIMENSION(:)::LIST_NIN,LIST_NOUT

KMAXE=XMPIELRANK(N)
ALLOCATE(NDLC_ARRAY1(KMAXE*8));NDLC_ARRAY1(:)=0

NDLC_COUNT1=0

ALLOCATE(LIST_NIN(1:IMAXN),LIST_NOUT(1:IMAXN))
LIST_NIN(:)=-10
LIST_NOuT(:)=-10

!FIND ONLY THE UNIQUES FIRST
COUNTFNODES=0
DO I=1,IMAXN
    IF (INODER(I)%ITOR.GT.0)THEN
    COUNTFNODES=COUNTFNODES+1
    LIST_NIN(I)=n
    END IF
END do



CALL MPI_ALLREDUCE(LIST_NIN,LIST_NOUT,IMAXN,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)

!FIND ONLY THE UNIQUES FIRST
COUNTFNODES=0
DO I=1,IMAXN
    IF (LIST_NOUT(i).EQ.N)THEN
    COUNTFNODES=COUNTFNODES+1
    END IF
END do




deallocate(LIST_NIN)



!FIND ONLY THE UNIQUES FIRST





DO I=1,IMAXN
    IF (INODER(I)%ITOR.GT.0)THEN

    NDLC_COUNT1=NDLC_COUNT1+1
    NDLC_ARRAY1(NDLC_COUNT1)=I
    END IF

END DO

KMAXN=NDLC_COUNT1



ALLOCATE(INODER4(1:KMAXN));
do i=1,kmaxn
    allocate(inoder4(i)%cord(1:dims))
    ALLOCATE(INODER4(I)%BCT(1:3));INODER4(I)%BCT(:)=0
end do

allocate(my_nodesl(1:COUNTFNODES),my_nodesg(COUNTFNODES))
COUNTFNODES=0
DO I=1,KMAXN

    INODER4(I)%ITOR=NDLC_ARRAY1(I)
    INODER4(I)%CORD(1:DIMS)=INODER(INODER4(I)%iTOR)%CORD(1:DIMS)
    INODER(INODER4(I)%ITOR)%ITOR=I
    if (LIST_NOUT(INODER4(I)%ITOR).eq.n)then
    COUNTFNODES=COUNTFNODES+1
    INODER4(I)%ITORm=COUNTFNODES
    my_nodesl(COUNTFNODES)=i
    my_nodesg(COUNTFNODES)=INODER4(I)%ITOR
    end if
END DO






ALLOCATE(XMPIALL_v(0:ISIZE-1),OFFSET_v(0:ISIZE-1))

XMPIALL_v=0

XMPIALL_v(N)=COUNTFNODES



call mpi_barrier(MPI_COMM_WORLD,IERROR)


CALL MPI_ALLGATHER(COUNTFNODES,1,MPI_INTEGER,XMPIALL_v,1,MPI_INTEGER,MPI_COMM_WORLD,IERROR)




OFFSET_v(0)=0
DO I=1,ISIZE-1
OFFSET_v(I)=OFFSET_v(I-1)+XMPIALL_v(I-1)
END DO










!NOW WE FIXED THE CORRECT NUMBERING FOR WRITING VTK OUTPUT



DO I=1,KMAXE
	if (tecplot.eq.5)then
	if (DIMENSIONA.EQ.2)then
	IF (ielem(n,i)%ishape.eq.5)then !QUAD
	indfc=4
	eLSE
	indfc=3
	end if
	allocate(IELEM(N,I)%NODES_v(1:indfc))
	eLSE
	IF (ielem(n,i)%ishape.eq.1)then !HEXA
	indfc=8
	end if
	IF (ielem(n,i)%ishape.eq.2)then !tetra
	indfc=4
	end if
	IF (ielem(n,i)%ishape.eq.3)then !pyramid
	indfc=5
	end if
	IF (ielem(n,i)%ishape.eq.4)then !prism
	indfc=6
	end if


	allocate(IELEM(N,I)%NODES_v(1:indfc))
	end if

	end if



    DO J=1,ielem(n,i)%NONODES
    IELEM(N,I)%NODES(J)=INODER(IELEM(N,I)%NODES(J))%ITOR
    END DO





    if (tecplot.eq.5)then
					if (DIMENSIONA.EQ.2)then
						IF (ielem(n,i)%ishape.eq.5)then !QUAD
						IELEM(N,I)%NODES_v(1)=INODER4(IELEM(N,I)%NODES(1))%ITOR-1
						IELEM(N,I)%NODES_v(2)=INODER4(IELEM(N,I)%NODES(2))%ITOR-1
						IELEM(N,I)%NODES_v(3)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1
						IELEM(N,I)%NODES_v(4)=INODER4(IELEM(N,I)%NODES(4))%ITOR-1
						eLSE							!TRIANGLE
						IELEM(N,I)%NODES_v(1)=INODER4(IELEM(N,I)%NODES(1))%ITOR-1
						IELEM(N,I)%NODES_v(2)=INODER4(IELEM(N,I)%NODES(2))%ITOR-1
						IELEM(N,I)%NODES_v(3)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1
						!IELEM(N,I)%NODES_v(4)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1

						END IF
					END IF
					if (DIMENSIONA.EQ.3)then
						IF (ielem(n,i)%ishape.eq.1)then !HEXA
						IELEM(N,I)%NODES_v(1)=INODER4(IELEM(N,I)%NODES(1))%ITOR-1
						IELEM(N,I)%NODES_v(2)=INODER4(IELEM(N,I)%NODES(2))%ITOR-1
						IELEM(N,I)%NODES_v(3)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1
						IELEM(N,I)%NODES_v(4)=INODER4(IELEM(N,I)%NODES(4))%ITOR-1
						IELEM(N,I)%NODES_v(5)=INODER4(IELEM(N,I)%NODES(5))%ITOR-1
						IELEM(N,I)%NODES_v(6)=INODER4(IELEM(N,I)%NODES(6))%ITOR-1
						IELEM(N,I)%NODES_v(7)=INODER4(IELEM(N,I)%NODES(7))%ITOR-1
						IELEM(N,I)%NODES_v(8)=INODER4(IELEM(N,I)%NODES(8))%ITOR-1
						END IF
						IF (ielem(n,i)%ishape.eq.2)then !TETRA
						IELEM(N,I)%NODES_v(1)=INODER4(IELEM(N,I)%NODES(1))%ITOR-1
						IELEM(N,I)%NODES_v(2)=INODER4(IELEM(N,I)%NODES(2))%ITOR-1
						IELEM(N,I)%NODES_v(3)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1
						IELEM(N,I)%NODES_v(4)=INODER4(IELEM(N,I)%NODES(4))%ITOR-1
						END IF
						IF (ielem(n,i)%ishape.eq.3)then !PYRAMID
						IELEM(N,I)%NODES_v(1)=INODER4(IELEM(N,I)%NODES(1))%ITOR-1
						IELEM(N,I)%NODES_v(2)=INODER4(IELEM(N,I)%NODES(2))%ITOR-1
						IELEM(N,I)%NODES_v(3)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1
						IELEM(N,I)%NODES_v(4)=INODER4(IELEM(N,I)%NODES(4))%ITOR-1
						IELEM(N,I)%NODES_v(5)=INODER4(IELEM(N,I)%NODES(5))%ITOR-1
						END IF
						IF (ielem(n,i)%ishape.eq.4)then !PRISM

						IELEM(N,I)%NODES_v(1)=INODER4(IELEM(N,I)%NODES(1))%ITOR-1
						IELEM(N,I)%NODES_v(2)=INODER4(IELEM(N,I)%NODES(2))%ITOR-1
						IELEM(N,I)%NODES_v(3)=INODER4(IELEM(N,I)%NODES(3))%ITOR-1
						IELEM(N,I)%NODES_v(4)=INODER4(IELEM(N,I)%NODES(4))%ITOR-1
						IELEM(N,I)%NODES_v(5)=INODER4(IELEM(N,I)%NODES(5))%ITOR-1
						IELEM(N,I)%NODES_v(6)=INODER4(IELEM(N,I)%NODES(6))%ITOR-1

						END IF

					END IF




			DO L=1,IELEM(N,I)%IFCA
					if (DIMENSIONA.EQ.2)then
					DO J=1,2
							indexgt=INODER(IELEM(N,I)%NODES_FACES(L,J))%ITOR
							IELEM(N,I)%NODES_FACES_v(L,J)=INODER4(indexgt)%ITOR-1
					END DO

					else
						if (ielem(n,i)%types_faces(L).eq.5)then
							DO J=1,4
								indexgt=INODER(IELEM(N,I)%NODES_FACES(L,J))%ITOR
								IELEM(N,I)%NODES_FACES_v(L,J)=INODER4(indexgt)%ITOR-1
							END DO
						END IF
						if (ielem(n,i)%types_faces(L).eq.6)then
							DO J=1,3
								indexgt=INODER(IELEM(N,I)%NODES_FACES(L,J))%ITOR
								IELEM(N,I)%NODES_FACES_v(L,J)=INODER4(indexgt)%ITOR-1

							END DO
						END IF
					END IF
			END DO

	end if



END DO



	DO I=1,KMAXE


    DO L=1,IELEM(N,I)%IFCA
         if (DIMENSIONA.EQ.2)then
          DO J=1,2
                IELEM(N,I)%NODES_FACES(L,J)=INODER(IELEM(N,I)%NODES_FACES(L,J))%ITOR

            END DO

            else
        if (ielem(n,i)%types_faces(L).eq.5)then
            DO J=1,4
                IELEM(N,I)%NODES_FACES(L,J)=INODER(IELEM(N,I)%NODES_FACES(L,J))%ITOR
            END DO
        END IF
        if (ielem(n,i)%types_faces(L).eq.6)then
            DO J=1,3
                IELEM(N,I)%NODES_FACES(L,J)=INODER(IELEM(N,I)%NODES_FACES(L,J))%ITOR
            END DO
        END IF



        END IF
    END DO
END DO














DEALLOCATE (INODER,NDLC_ARRAY1,LIST_NOuT)


IF (TECPLOT.EQ.3)then

if (dimensiona.eq.3)then

allocate(el_connect(1:kmaxe,1:8))
do i=1,kmaxe
    if (ielem(n,i)%ishape.eq.1)then!hexa
    el_connect(i,1:8)=ielem(n,i)%nodes(1:8)
    end if
     if (ielem(n,i)%ishape.eq.2)then!tetra
    el_connect(i,1:2)=ielem(n,i)%nodes(1:2)
    el_connect(i,3)=ielem(n,i)%nodes(3)
    el_connect(i,4)=ielem(n,i)%nodes(3)
    el_connect(i,5)=ielem(n,i)%nodes(4)
    el_connect(i,6)=ielem(n,i)%nodes(4)
    el_connect(i,7)=ielem(n,i)%nodes(4)
    el_connect(i,8)=ielem(n,i)%nodes(4)
    end if
    if (ielem(n,i)%ishape.eq.3)then!pyramid
    el_connect(i,1:5)=ielem(n,i)%nodes(1:5)
    el_connect(i,6)=ielem(n,i)%nodes(5)
    el_connect(i,7)=ielem(n,i)%nodes(5)
    el_connect(i,8)=ielem(n,i)%nodes(5)
    end if
    if (ielem(n,i)%ishape.eq.4)then!prism
    el_connect(i,1:3)=ielem(n,i)%nodes(1:3)
    el_connect(i,4)=ielem(n,i)%nodes(3)
    el_connect(i,5:7)=ielem(n,i)%nodes(4:6)
    el_connect(i,8)=ielem(n,i)%nodes(6)
    end if

end do

END IF

END IF
!array el_connect(index1,indices 1:8)
!first index is the element index
!second index varies from 1:8 for all the vertices of this element (connectivity list)
!for paraview the numbering might require switching (from 1 to number of nodes----> 0  to number of nodes-1)

!INODER4(1:number of nodes (kmaxn))%cord(1:3) holds the coordinates for each point 
!kmaxn is global across all modules
!kmaxe must be set as equal to xmpielrank(n)
!this is shown in OUTWRITEPARA3DbP (P STANDS for partitioned mesh writing)















END SUBROUTINE FIX_NODES_LOCAL




SUBROUTINE CHECKPOINTv4(N)
!> @brief
!> This subroutine is writing the checkpointing files
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS
REAL,ALLOCATABLE,DIMENSION(:,:)::xbin
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj,igfs
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV
REAL::MP_PINFl,gammal
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT
 KMAXE=XMPIELRANK(N)
igfs=t

WRITE(PROC3,FMT='(I10)') igfs
	RESTFILE="REST_"//TRIM(ADJUSTL(PROC3))//".dat"!//TRIM(ADJUSTL(PROC4))
	
ICPUID=N

	IF (N.EQ.0)THEN
	OPEN(1086,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
	end if
	
	KMAXE=XMPIELRANK(N)
    
DUMG=KMAXE


call mpi_barrier(mpi_comm_world,IERROR) !NOT NEEDED

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


    IF (N.EQ.0)THEN
	    ALLOCATE(ICELLA(IMAXP*ISIZE))
	    ICELLA=0

    END IF

    call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

    deallocate (icell)

    IF (N.EQ.0) then
    ALLOCATE(VALUESA(IMAXP*ISIZE))
    allocate(xbin(imaxe,5+turbulenceequations+passivescalar))
    VALUESA=ZERO
    END IF
    
    
    ALLOCATE(VALUESS(imaxp));VALUESS=ZERO

  
  do jj=1,5
	DO I=1,KMAXE
        LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		VALUESS(I)=leftv(jj)
	END DO
	call MPI_GATHER(VALUESS,imaxp,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
	
	    IF (N.EQ.0)THEN
	    do i=1,imaxp*isize
		if (icella(i).gt.0)then
		xbin(icella(i),jj)=valuesa(i)
		end if
	    end do
	    end if
    
  end do



!   CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN

    DO I=1,IMAXE
    WRITE(1086)XBIN(i,1:nof_Variables)
    END DO

    DEALLOCATE(XBIN,ICELLA,VALUESA)
    close(1086)
    END IF
    
DEALLOCATE(VALUESS)







END SUBROUTINE CHECKPOINTv4


SUBROUTINE CHECKPOINTv3(N)
!> @brief
!> This subroutine is writing the checkpointing files
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER,ALLOCATABLE,DIMENSION(:)::ICELL,ICELLA
REAL,ALLOCATABLE,DIMENSION(:)::VALUESA,VALUESS
REAL,ALLOCATABLE,DIMENSION(:,:)::xbin
INTEGER::I,K,KMAXE,J,JK,ICPUID,nvar,IMAXP,DUMG,DUML,jj
CHARACTER(LEN=20)::PROC,RESTFILE,PROC3
REAL,ALLOCATABLE,DIMENSION(:)::IGINT,TGINT
 KMAXE=XMPIELRANK(N)



	RESTFILE="CORD.dat"
	
ICPUID=N

	IF (N.EQ.0)THEN
	OPEN(1086,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='REPLACE',ACTION='WRITE')
	end if
	
	KMAXE=XMPIELRANK(N)
    
DUMG=KMAXE


call mpi_barrier(mpi_comm_world,IERROR) !NOT NEEDED

CALL MPI_ALLREDUCE(DUMG,DUML,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,IERROR)
IMAXP=DUML

ALLOCATE(ICELL(IMAXP))
ICELL=0

DO I=1,KMAXE
  ICELL(I)=IELEM(N,I)%IHEXGL
END DO


    IF (N.EQ.0)THEN
	    ALLOCATE(ICELLA(IMAXP*ISIZE))
	    ICELLA=0

     ELSE
		ALLOCATE(ICELLA(1))

     END IF

    call MPI_GATHER(ICELL,IMAXP,MPI_INTEGER,icella,imaxp,mpi_integer,0,MPI_COMM_WORLD,IERROR)

    deallocate (icell)

	IF (N.EQ.0) then
    ALLOCATE(VALUESA(IMAXP*ISIZE))

    allocate(xbin(imaxe,3))
    VALUESA=ZERO
    ELSE
    ALLOCATE(VALUESA(1),xbin(1,3))

    END IF
    
    
    ALLOCATE(VALUESS(imaxp));VALUESS=ZERO

  
  do jj=1,3
	DO I=1,KMAXE
        IF (JJ.EQ.1)THEN
		VALUESS(I)=IELEM(N,I)%XXC
		END IF
		IF (JJ.EQ.2)THEN
		VALUESS(I)=IELEM(N,I)%YYC
		END IF
		IF (JJ.EQ.3)THEN
		VALUESS(I)=IELEM(N,I)%ZZC
		END IF
	END DO
	call MPI_GATHER(VALUESS,imaxp,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
	
	    IF (N.EQ.0)THEN
	    do i=1,imaxp*isize
		if (icella(i).gt.0)then
		xbin(icella(i),jj)=valuesa(i)
		end if
	    end do
	    end if
    
  end do



!   CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    IF (N.EQ.0)THEN

    DO I=1,IMAXE
    WRITE(1086)XBIN(i,1:3)
    END DO


    close(1086)
    END IF

    DEALLOCATE(XBIN,ICELLA,VALUESA)


    
DEALLOCATE(VALUESS)







END SUBROUTINE CHECKPOINTv3


SUBROUTINE TROUBLED_HISTORY
INTEGER::I,J,K,TRAJ1,TRAJ2,TRAJ3,TRAJ4,kmaxe,writeid,writeconf
REAL::WIN1,WIN2,WIN3,WIN4,POST,POST1,POST2,POST3,POST4
real,dimension(3)::pos_l,pos_g
integer,dimension(3)::ipos_l,ipos_g
KMAXE=XMPIELRANK(N)
POST1=0
traj1=0
pos_l(1)=zero
pos_G(1)=zero
ipos_l(:)=0
ipos_G(:)=0

if (mood.gt.0)then

DO I=1,KMAXE
	if (IELEM(N,I)%mood_o.lt.(iorder+1))then
    ipos_l(1)=ipos_l(1)+1		!number of cells
    end if
    if (IELEM(N,I)%mood_o.EQ.1)then
    ipos_l(2)=ipos_l(2)+1
    end if
END DO


CALL MPI_ALLREDUCE(Ipos_l(1:2),Ipos_g(1:2),2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,IERROR)

POS_G(1)=IPOS_G(1)
POS_G(2)=IPOS_G(2)
POS_G(1)=(POS_G(1)/IMAXE)*100
POS_G(2)=(POS_G(2)/IMAXE)*100

IF (n.eq.0)THEN

OPEN(70,FILE='TROUBLED.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),POS_G(2)
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)






eLSE

DO I=1,KMAXE
    pos_l(1)=pos_l(1)+IELEM(N,I)%CONDITION
END DO

  
CALL MPI_ALLREDUCE(pos_l(1),pos_g(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

IF (n.eq.0)THEN

OPEN(70,FILE='TROUBLED.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7)')T,(POS_G(1)/IMAXE)*100.0
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


end if

  






END SUBROUTINE TROUBLED_HISTORY


SUBROUTINE REDUCED_HISTORY
INTEGER::I,J,K,TRAJ1,TRAJ2,TRAJ3,TRAJ4,kmaxe,writeid,writeconf
REAL::WIN1,WIN2,WIN3,WIN4,POST,POST1,POST2,POST3,POST4
real,dimension(1)::pos_l,pos_g
KMAXE=XMPIELRANK(N)
POST1=0
traj1=0
pos_l(1)=zero
pos_G(1)=zero
DO I=1,KMAXE
    pos_l(1)=pos_l(1)+IELEM(N,I)%REDUCE
END DO


CALL MPI_ALLREDUCE(pos_l(1),pos_g(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

IF (n.eq.0)THEN

OPEN(70,FILE='REDUCED.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7)')T,(POS_G(1)/IMAXE)*100.0
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)



END SUBROUTINE REDUCED_HISTORY


SUBROUTINE FILTERED_HISTORY
INTEGER::I,J,K,TRAJ1,TRAJ2,TRAJ3,TRAJ4,kmaxe,writeid,writeconf,countfd
REAL::WIN1,WIN2,WIN3,WIN4,POST,POST1,POST2,POST3,POST4
real,dimension(5)::pos_l,pos_g
KMAXE=XMPIELRANK(N)
POST1=0
traj1=0
pos_l(:)=zero
pos_G(:)=zero
! pos_L(2)=10e20
! pos_g(2)=0
! pos_L(3)=0.0d0
! pos_g(3)=0
! pos_L(4)=-10e20
! pos_g(4)=0
! countfd=0
DO I=1,KMAXE

    pos_l(1)=pos_l(1)+IELEM(N,I)%FILTERED

!     if (ielem(n,i)%er2dt.gt.0)THen
!
!
! 			if (ielem(n,i)%er1er2.gt.1)then
!
!
! 			countfd=countfd+1
!
! 			pos_l(2)=min(pos_l(2),ielem(n,i)%er1er2)
!
! 			pos_l(3)=ielem(n,i)%er1er2+pos_l(3)
!
!
!
!
!
! 			end if








! 	end if


! 	if (ielem(n,i)%er2dt.lt.0.0)THen
! 			pos_l(4)=max(pos_L(4),ielem(n,i)%er2dt)
!
! 			end if




!
!
!
!
!  			countfd=countfd+1
!

!
 			pos_l(3)=ielem(n,i)%er+pos_l(3)




!  			if (ielem(n,i)%er.gt.0)then

 			pos_l(2)=ielem(n,i)%er1+pos_l(2)
!
 			pos_l(4)=ielem(n,i)%er2+pos_l(4)

!  			end if
















END DO

! pos_l(3)=pos_l(3)/countfd


CALL MPI_ALLREDUCE(pos_l(1),pos_g(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)


post1=(POS_G(1)/IMAXE)*100.0

CALL MPI_ALLREDUCE(pos_l(2),pos_g(2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)



post2=(POS_G(2)/IMAXE)

CALL MPI_ALLREDUCE(pos_l(3),pos_g(3),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)


post3=(POS_G(3)/IMAXE)

CALL MPI_ALLREDUCE(pos_l(4),pos_g(4),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

post4=(POS_G(4)/IMAXE)



IF (n.eq.0)THEN

OPEN(70,FILE='FILTERED.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POST1,POST2,POST3,POST4
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)



END SUBROUTINE FILTERED_HISTORY



END MODULE IO
