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
     
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                    'Density,U,V,W,energy,Pressure,STEN1,STEN2'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     
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
	

 END IF
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
		  VALUESS(i)=n
		END DO
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
		IF (N.EQ.0)THEN
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
		ierr = TECDAT112(imaxe,xbin,1)
		END IF
    
		DO I=1,KMAXE
		  VALUESS(i)=IELEM(N,I)%iNUMNEIGHBOURS!%STENCIL_DIST
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
		do kkd=1,nof_variables
		    DO I=1,KMAXE
		      VALUESS(i)=U_C(I)%VAL(1,kkd)
			if ((kkd.ge.2).and.(kkd.le.nof_variables-1))then
			VALUESS(i)=U_C(I)%VAL(1,kkd)/U_C(I)%VAL(1,1)
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
		  leftv(1:5)=U_C(I)%VAL(1,1:5)
		  CALL CONS2PRIM(N)
		  VALUESS(i)=leftv(5)
		END DO
		
		
		call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF

		
		
		
		
		
		
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%GGS!%ISHAPE!IELEM(N,I)%STENCIL_DIST
                END DO
                
                call MPI_GATHERv(valuess,xmpiall(n),MPI_DOUBLE_PRECISION,xbin2,xmpiall,offset,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
			IF (N.EQ.0)THEN
			do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
			ierr = TECDAT112(imaxe,xbin,1)
			END IF

                
               
                
                
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%STENCIL_DIST!ADMIS
                END DO
                
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
    
		  if (itestcase.eq.4)then
		  DO I=1,KMAXE
		      VALUESS(i)=ielem(n,i)%vortex(1)!%inumneighbours
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
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2)
  END IF
  
  DEALLOCATE (VALUESS,VARIABLES)
  

  





	
	

END SUBROUTINE OUTWRITE3vb


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
		  leftv(1:5)=U_C(I)%VAL(1,1:5)
		  CALL CONS2PRIM(N)
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
		  leftv(1:4)=U_C(I)%VAL(1,1:4)
		  CALL CONS2PRIM2d(N)
		  VALUESS(i)=ielem(n,i)%condition!leftv(4)
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
     
     if (n.eq.0)ierr =  TecIni112('sols'//NULCHAR, &
                     'Density,U,V,energy,Pressure,STEN1,STEN2,SLOPE'//NULCHAR, &
                    out1//NULCHAR, &
                    '.'//NULCHAR, &
                    FileType, &
                    Debug, &
                    VIsDouble)
     
     
     
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


 END IF
 
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
      VALUESS(i)=IELEM(N,I)%ADMIS
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
		  VALUESS(i)=U_C(I)%VAL(1,kkd)
		  if (kkd.ge.2.and.kkd.le.3)then
		  VALUESS(i)=U_C(I)%VAL(1,kkd)/U_C(I)%VAL(1,1)
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
		  leftv(1:4)=U_C(I)%VAL(1,1:4)
		  CALL CONS2PRIM2D(N)
		  VALUESS(i)=leftv(4)
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
                
                DO I=1,KMAXE
                VALUESS(i)=IELEM(N,I)%STENCIL_DIST
                END DO
                
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
		  
		  end if
    
		  if (itestcase.eq.4)then
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
  DEALLOCATE(XBIN,xbin2,valuelocation,out1)
  END IF
  
  DEALLOCATE (VALUESS,variables)


	
	

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
	RES_TIME=0.0
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

SUBROUTINE READ_UCNS3D
!> @brief
!> This subroutine reads the parameter file

	IMPLICIT NONE

 	Integer :: INV
 	Real :: angledum
	CHARACTER(48)::STAMP1
	LOGICAL::HERE1


 	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	INQUIRE (FILE='RESTARTav.dat',EXIST=HERE1)
	IF (HERE1) THEN
	Average_restart=1
	Else
	Average_restart=0
	end if
	
	


	
	
	
	
	
			
	
	
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	
	
	
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)DIMENSIONA,STATISTICS,CODE_PROFILE
	READ(15,*)
	READ(15,*)
	READ(15,*)governingequations,INITCOND
	READ(15,*)
	READ(15,*)
	READ(15,*)turbulence,icoupleturb,PASSIVESCALAR
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)RRES,ufreestream,VVEL,WVEL,PRES
	READ(15,*)
	READ(15,*)
	READ(15,*)AOA,vectorx,vectory,vectorz
	READ(15,*)
	READ(15,*)
	READ(15,*)GAMMA,PRANDTL,Reynolds,CharLength
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)spatiladiscret,iRiemann,spatialorder,LIMITER,POLY
	READ(15,*)
	READ(15,*)
	READ(15,*)wenocnschar,EES,wenoz,wenocentralweight
	READ(15,*)
	READ(15,*)
	READ(15,*)temporder,CFL,timestep,upperlimit,reSLIMIT
	READ(15,*)
	READ(15,*)
	READ(15,*)iboundary,boundtype,SCALER
	READ(15,*)
	READ(15,*)
	READ(15,*)GREENGO,LMACH
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)OUT_TIME,NTMAX,WALLC
	READ(15,*)
	READ(15,*)
	READ(15,*)TECPLOT,IEVERY,IEVERY2,IEVERYAV,STENCIL_IO
	READ(15,*)
	READ(15,*)
	READ(15,*)Averaging
	READ(15,*)
	READ(15,*)
	READ(15,*)OUTSURF,IFORCE,surfshear
	READ(15,*)
	READ(15,*)
	READ(15,*)IRES_TURB,IRES_UNSTEADY,LAMPS,Prev_turbmodel  
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)NPROBES
	READ(15,*)
	    
	    
	    
	    
	    SELECT CASE(CODE_PROFILE)
	
	
	CASE (0)
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	icoupleturb=0	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=2.2		!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	if (dimensiona.eq.2)then
	if (tecplot.eq.2)then
	tecplot=1
	end if
	end if
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 
	 
	 DES_model=0
	 
	 
	 
	CASE (1)
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	icoupleturb=1	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=3		!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	if (dimensiona.eq.2)then
	if (tecplot.eq.2)then
	tecplot=1
	end if
	end if
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 DES_model=2
	 
	 
	 CASE (3)
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	icoupleturb=0	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=2.5	!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=10.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=30.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	if (dimensiona.eq.2)then
	if (tecplot.eq.2)then
	tecplot=1
	end if
	end if
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 
	 
	 DES_model=2
	 
	 
	  CASE (4)
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	icoupleturb=0	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=1	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=2	!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=10.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=30.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	if (dimensiona.eq.2)then
	if (tecplot.eq.2)then
	tecplot=1
	end if
	end if
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 
	 
	 DES_model=0
	 
	 
	 
	
	CASE DEFAULT
	
	LOWMEMORY=0 	 !MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	 !I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	 !GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	 !QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	icoupleturb=0	 !COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	 !HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0  !HYBRID DISTANCE
	swirl=0		 !swirling flow:0 deactivated, 1 activated
	IADAPT=0	 !ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	 !COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=2		 !STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	 !WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0 !GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	 !STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		 !RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	 !CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	 !CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	 !Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	 !TOLERANCE=n_iterations
	GRIDAR1=5.0	 ! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	 ! LIMIT VOLUME CELLS
	fastest=0	 ! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	 !0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	 !LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	if (dimensiona.eq.2)then
	if (tecplot.eq.2)then
	tecplot=1
	end if
	end if
	 if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	if (dimensiona.eq.2)then
	if (tecplot.eq.2)then
	tecplot=1
	end if
	end if
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	
	
	DES_model=0
	
	
	
	
	END SELECT
	    
	    
	    
	    
	   
	    
	    
	    
	    
	    
	    
	    
	    
	
	 IF (N.EQ.0)THEN
      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
      CALL DATE_and_time(STAMP1)
      WRITE(63,*)"UCNS3D Parallel Computation Start",STAMP1
      END IF

      
      
	    NPROC=isize-1;vorder=iorder
      
	    !--------------------------1------------------------------!
	    !PROBE POSITION ALLOCATION AND READING OF COORDINATES
	    !$OMP MASTER
	    IF (NPROBES.GT.0)THEN
	    IF (DIMENSIONA.EQ.3)THEN
		  ALLOCATE(PROBEC(1:NPROBES,1:3))
		    DO INV=1,NPROBES
		      READ(15,*)PROBEC(INV,1),PROBEC(INV,2),PROBEC(INV,3)
		    END DO
	    ELSE
		  ALLOCATE(PROBEC(1:NPROBES,1:2))
	    
		    DO INV=1,NPROBES
		      READ(15,*)PROBEC(INV,1),PROBEC(INV,2)
		    END DO
	    END IF
	    END IF
	    !$OMP END MASTER
	     !--------------------------END 1-------------------------!
	     
	     
	     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	     
	     
	     !-------------------------2---------------------------------!
	     !NUMBER OF VARIABLES SETUP
	     if (governingequations.le.2)then
		  if (dimensiona.eq.3)then
		      nof_variables=5;dims=3
		  else
		      nof_variables=4;dims=2
		  end if
	      else
		  nof_variables=1
		    if (dimensiona.eq.3)then
		    dims=3
		    else
		    dims=2
		    end if
	      end if
	    !--------------------------END 2-------------------------!


	    !-------------------------3---------------------------------!
	     !TURBULENCE
	      if (turbulence .eq.1) then
		      if (turbulencemodel .eq. 1) then
			      turbulenceequations=1
		      else if (turbulencemodel .eq. 2) then
			      turbulenceequations=2
		      end if 
	      eLSE
		      turbulenceequations=0
		      turbulencemodel=0
	      END IF
	    !------------------------- END 3---------------------------------!
	   
	   
! 	   !-------------------------4---------------------------------!
	   !EQUATIONS TYPE
	   
	   SELECT CASE(governingequations)
	   
	   
	   CASE(1)
	   !NAVIER STOKES EQUATIONS
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Navier-Stokes Solver Engaged'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 4;IVORTEX = 1
	    
	  CASE(2)
	   !EULER EQUATIONS
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Euler Solver Engaged'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 3;IVORTEX = 0
	    
	    
	    
	  CASE(3)
	   !LINEAR ADVECTION EQUATION
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Linear Solver Engaged Sinewave'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 1;IVORTEX = 0
	    
	 CASE(4)
	   !LINEAR ADVECTION EQUATION
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Linear Solver Engaged STEP FUNCTION'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 2;IVORTEX = 0
	    
	  END SELECT
	  
	  
	  
	  !-------------------------END EQUATIONS 4---------------------------------!
	  
	  !-------------------------5---------------------------------!
	  !FLOW PARAMETERS
	  UNWOU = 3
	  BETAAS=1.5d0
	  SUTHER=0.412158681d0
	  uvel=ufreestream
	   ! Set pressure
	  if ( PRES .lt. 0 ) PRES = RRES/GAMMA	
	  ! Set dynamic free-stream viscosity
	  VISC = (RRES*ufreestream*CharLength)/Reynolds
	  if (swirl.eq.1)then
	  uvel=ZERO
	  end if
	  If (AOA .NE. 0.0D0) Then
	  angledum=(AOA*PI)/180.0d0
	  UVEL = COS(angledum)*ufreestream*VECTORX
	  WVEL = SIN(angledum)*ufreestream*VECTORZ
	  VVEL = SIN(angledum)*ufreestream*VECTORY
	  end if
	  
	  IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      if (ITESTCASE .eq. 4) Then
		write(63,*)'----Reynolds Number:',(RRES*ufreestream*CharLength)/VISC
		end if
	      CLOSE(63)
	   END IF
	  
	  
	  
	  !-------------------------END FLOW PARAMETERS 5---------------------------------!
	  
	  !-------------------------6---------------------------------!
	   !SPATIAL & TEMPORAL DISCRETISATION
	   
	   IORDER = max(1,spatialorder-1)
	   firstorder = 0
! 	   if (ees.eq.5)then
! 	   icompact=1
! 	   end if
	   IEXTEND =2; IF (ICOMPACT.EQ.1)IEXTEND = 8
	   IWENO =0
	   WENWRT = 1
	   ISCHEME=3
	   TYPESTEN=1
	   LWCI1=wenocentralweight
	  if (fastest.eq.1)ischeme=1
	    
	   
	   SELECT CASE(spatiladiscret)
	   
	   CASE(1)	!NO LIMITER
	      
	   IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'----1st Order in Space|LINEAR Scheme'
	      CLOSE(63)
	   END IF
	   
	      if (spatialorder.eq.1) then 
	      firstorder = 1;iorder=1
	      else
	      firstorder = 0
	      end if

	     CASE(2)	!MUSCL TYPE
	      
	   IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'----MUSCL Scheme Engaged----'
	      CLOSE(63)
	   END IF
	  
	  IWENO = -1
	  
	  
	   CASE(3)	!WENO TYPE
	      
	   IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'----WENO Scheme Engaged----'
	      CLOSE(63)
	   END IF
	  
	  IWENO = 1
	  IWENO = 1
	  IEXTEND = 12	;
	  
	  
	  if (EES.EQ.5)THEN
	   IF (ICOMPACT.EQ.1)THEN
	  IEXTEND = 8
	  
	  
	  
	  ELSE
	  IEXTEND = 5
	  
	  END IF
	  END IF
	  WENWRT=wenocnschar
	  
	  IF (DIMENSIONA.EQ.2)THEN
	    if ((EES.EQ.0).OR.(EES.GE.4))THEN;TYPESTEN = 5;End if
	    if (EES.EQ.1)THEN;TYPESTEN = 9;end if
	    if (EES.EQ.2)THEN;TYPESTEN = 9;End if
	    if (EES.EQ.3)THEN;TYPESTEN = 5;end if
	  ELSE
	    if ((EES.EQ.0).OR.(EES.GE.4))THEN;TYPESTEN = 7;End if
	    if (EES.EQ.1)THEN ;TYPESTEN = 15;End if
	    if (EES.EQ.2)THEN;TYPESTEN = 15;End if
	    if (EES.EQ.3)THEN;TYPESTEN = 7;End if
	  END IF
	  
	  
	  
	  END SELECT
	  
	  IF (N.EQ.0)THEN
	    OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	    write(63,*)'Order of Accuracy in space:',spatialorder
	    CLOSE(63)
	END IF
	   
	  ! Temporal order
	  RUNGEKUTTA = temporder 
	  
	if ( iboundary .eq. 0 ) then 
	    IPERIODICITY = -3 
	    end if
	    if ( iboundary .eq. 1 ) then 
	    IPERIODICITY = 1 
	    end if

	    if (guassianquadra.eq.0)then
	  IGQRULES=min(iorder,6)
	else

	if (guassianquadra.gt.1)then
	IGQRULES = guassianquadra-1
	  else
	IGQRULES =guassianquadra
	end if
	End if
	    
	    
	  !-------------------------END DISCRETISATION 6---------------------------------!




	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Total Number of Processes:',isize
	      CLOSE(63)
	  END IF

	  
	   
	  
	  !TURBULENCE DEFAULT VALUES
	
	
			! 	TURBULENCE MODEL PARAMETERS:
			! ---OPTIONS Spalart Almaras---
			ISPAL=1! 1				||SPALART ALLMARAS VERSION:| 1:ORIGINAL |2: NEGATIVE MODIFICATION
			!DES_model=0! 0              			|| 1-Detached Eddy Simulation 2-Delayed DES
			! ---CONSTANTS-------
			CB1=0.1355	! 0.1355				|| Cb1 
			CB2=0.622! 0.622   			|| Cb2
			SIGMA=0.666666667! 0.666666667      		|| sigma
			KAPPA=0.41  ! 0.41     			|| kappa
			CW1=3.23886781677! 3.23886781677      		|| CW1 (0 for (CB1/kappa**2)+1/sigma(1+CB2) or predifined)
			CW2=0.3! 0.3      			|| CW2
			CW3=2.0 ! 2.0      			|| CW3
			CV1=7.1! 7.1      			|| CV1
			CT1=1.0 ! 1.0      			|| CT1
			CT2=2.0! 2.0      			|| CT2
			CT3=1.1 ! 1.1      			|| CT3
			CT4=2.0 ! 2.0      			|| CT4
			PRTU=0.9! 0.9				|| Turbulent Prandtl Number
			TWALL=0! 0				|| Wall temperature (Kelvin) leave 0 for adiabatic (q_wall =0 <=> dT/dn=0)
			TURBINIT=3.0! 3.0	  			|| Initial value for turbulence parameter (multiplyied by the freestream viscosity from given Re)
			Upturblimit=1000000! 1000000				|| Upper limit for turbulence
			residualfreq=10! 10				|| Residual compute every
			IRS=0! 0				||IMPLICIT RESIDUAL SMOOTHING (DOUBLES CFL)
			C_DES_SA=0.61	! 0.61				||C_DES_SA
			! =============================================================================
			! K-OMEGA SST: 
			! ---OPTIONS---
			VORT_MODEL=0! 0					||0:Default strain-production for k, 1:Vorticity-production for k
			QSAS_MODEL=0! 0					||Hybrid models: 1-SAS-SST  //  2-DES-SST 
			ZERO_TURB_INIT=0! 0					||Zero turbulence option (if 1, initialization for k/w is done with zero turbulence)
			! --------------K-OMEGA CONSTANTS--------------------------
			SIGMA_K1=1.176470588! 1.176470588				||sigma_k1
			SIGMA_K2=1.0! 1.0					||sigma_k2
			SIGMA_OM1=2.0	! 2.0					||sigma_om1
			SIGMA_OM2=1.168! 1.168					||sigma_om2
			AA_1=0.31! 0.31					||aa_1
			BETA_I1=0.075	! 0.075					||beta_i1
			BETA_I2=0.0828! 0.0828					||beta_i2
			ALPHA_STARINF=1.0! 1.0					||alpha_starinf
			ALPHA_0=0.111111111! 0.111111111				||alpha_0
			BETA_STARINF=0.09! 0.09					||beta_starinf
			R_BETA=8.0! 8.0					||R_beta
			R_K_SST=6.0! 6.0					||R_k_SST
			BETA_T=0.072	! 0.072					||beta_t
			KAPPA_sst=0.41! 0.41					||kappa   (Same in Spalart-Allmaras)
			R_OM_SST=2.95	! 2.95					||R_om_SST
			ZETA_STAR= 1.5	! 1.5					||zeta_star  (Only for Mach corrections)
			M_T0=0.25! 0.25					||M_t0		(Only for Mach corrections)
			C_MU_INLET=0.09! 0.09					||C_mu_inlet
			C_SMG=0.11! 0.11					||C_smg  (SAS)
			ETA2_SAS=3.51! 3.51					||eta2_SAS (SAS)
			SIGMA_PHI=0.666666667	! 0.666666667				||sigma_phi (SAS)
			C_SAS=2.0! 2.0					||C_SAS   (SAS)
			L_turb_inlet=0.0026! 0.0026					||L_turb_inlet (Default turbulence lengthscale)
			I_turb_inlet=0.1! 0.1					||I_turb_inlet (Default turbulene intensity)
			Init_mu_ratio=0.01! 0.01					||Init_mu_ratio	(Default initial ratio for turbulent viscosity)
			C_DES_SST=0.61! 0.61					||C_DES_SST
			! =============================================================================
			! PASSIVE SCALAR TRANSPORT
			SCHMIDT_LAM=10.0! 10.0					||Laminar Schmidt number
			SCHMIDT_TURB=0.7! 0.7					||Turbulent Schmidt number

		
	

	END SUBROUTINE READ_UCNS3D


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
	!WRITE(400+N,*)KMAXE,IMAXE,IMAXN
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
Type :: Wallboundary
	    Integer :: GlobalID,wbid,wb1,wb2,wb3,wb4,NumNodes!,wbdescr
	    Real :: Wallx,Wally,Wallz
End Type
Type(Wallboundary),Allocatable,Dimension(:):: Wallbnd
Type :: NodesWall
	    Integer :: ID
	    Real :: wnx,wny,wnz
End Type
Type(NodesWall),Allocatable,Dimension(:):: Wallvrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! How many element for this block have wall
KMAXE=XMPIELRANK(N)
 countwall = 0
 Do i=1,KmaxE
    if (ielem(n,i)%interior.eq.1)then
    Do l=1,IELEM(N,I)%IFCA
	if (ielem(n,i)%ibounds(l).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.4)then
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
Type :: Wallboundary
	    Integer :: GlobalID,wbid,wb1,wb2,wb3,wb4,NumNodes!,wbdescr
	    Real :: Wallx,Wally,Wallz
End Type
Type(Wallboundary),Allocatable,Dimension(:):: Wallbnd
Type :: NodesWall
	    Integer :: ID
	    Real :: wnx,wny,wnz
End Type
Type(NodesWall),Allocatable,Dimension(:):: Wallvrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! How many element for this block have wall
KMAXE=XMPIELRANK(N)
 countwall = 0
 Do i=1,KmaxE
    if (ielem(n,i)%interior.eq.1)then
    Do l=1,IELEM(N,I)%IFCA
	if (ielem(n,i)%ibounds(l).gt.0)then
	      if (ibound(n,ielem(n,i)%ibounds(l))%icode.eq.4)then
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
INTEGER::INX,I,K,J,M,O,P,Q,JK
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
 




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
INTEGER::INX,I,K,J,M,O,P,Q,JK
LOGICAL::HEREV
REAL,DIMENSION(5)::TOTAL
 CHARACTER(LEN=20)::PROC,OUTFILE,PROC3,SURFILE
 




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
      if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      else
      SolTime = IT

      end if
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
                                                
                                                
                                                
                                                
                                               
					  
                                            
						leftv(1:5)=U_C(IBOUND_T(I))%VAL(1,1:5)
						CALL CONS2PRIM(N)
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
						
					 ICONSI=IBOUND_T(I)
					 ICONSIDERED=IBOUND_T2(I)
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y(ICONSI,ICONSIDERED)
					 CASE(3)
					 call SHEAR_Z(ICONSI,ICONSIDERED)
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
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2)
  END IF
!   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS)
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
      if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      else
      SolTime = IT

      end if
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
					  leftv(1:4)=U_C(IBOUND_T(I))%VAL(1,1:4)
				    CALL CONS2PRIM2d(N)
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
					  ICONSI=IBOUND_T(I)
					  ICONSIDERED=IBOUND_t2(i)
					  select case(kkd)
					 case(1)
					 
					 call SHEAR_X2d(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y2d(ICONSI,ICONSIDERED)
					 
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
  DEALLOCATE(XBIN,valuelocation,out1,XBIN2)
  END IF
  
!   IF (TOTIW.GT.0)THEN
  DEALLOCATE (VALUESS)
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
				    leftv(1:5)=U_C(IBOUND(N,i)%which)%VAL(1,1:5)
				    CALL CONS2PRIM(N)
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
					ICONSI=IBOUND(N,i)%which
					 ICONSIDERED=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y(ICONSI,ICONSIDERED)
					 CASE(3)
					 call SHEAR_Z(ICONSI,ICONSIDERED)
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
      if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      else
      SolTime = IT

      end if
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
				    leftv(1:4)=U_C(IBOUND(N,i)%which)%VAL(1,1:4)
				    CALL CONS2PRIM2d(N)
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
					ICONSI=IBOUND(N,i)%which
					 ICONSIDERED=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X2d(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y2d(ICONSI,ICONSIDERED)
					 
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
		   leftv(1:5)=U_C(I)%VAL(ind1,1:5)
		  CALL CONS2PRIM(N)
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
		   leftv(1:4)=U_C(I)%VAL(ind1,1:4)
		  CALL CONS2PRIM(N)
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
		   leftv(1:5)=U_C(I)%VAL(5,1:5)
		  CALL CONS2PRIM(N)
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
		   leftv(1:4)=U_C(I)%VAL(5,1:4)
		  CALL CONS2PRIM(N)
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
      if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      else
      SolTime = IT

      end if
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
					        leftv(1:5)=U_C(IBOUND_T(I))%VAL(ind1,1:5)
						
						CALL CONS2PRIM(N)
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
					        
				    
						
						ICONSI=IBOUND_T(I)
					 ICONSIDERED=IBOUND_T2(I)
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X_av(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y_av(ICONSI,ICONSIDERED)
					 CASE(3)
					 call SHEAR_Z_av(ICONSI,ICONSIDERED)
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
      if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      else
      SolTime = IT

      end if
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
						leftv(1:4)=U_C(IBOUND_T(I))%VAL(5,1:4)
						CALL CONS2PRIM2d(N)
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
						ICONSI=IBOUND_T(I)
					 ICONSIDERED=IBOUND_T2(I)
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X2d_av(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y2d_av(ICONSI,ICONSIDERED)
					 
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
				    leftv(1:5)=U_C(IBOUND(N,i)%which)%VAL(5,1:5)
				    CALL CONS2PRIM(N)
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
					ICONSI=IBOUND(N,i)%which
					 ICONSIDERED=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X_av(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y_av(ICONSI,ICONSIDERED)
					 CASE(3)
					 call SHEAR_Z_av(ICONSI,ICONSIDERED)
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
      if (( RUNGEKUTTA .LT. 5).or.( RUNGEKUTTA .eq. 11)) Then
      SolTime = T
      else
      SolTime = IT

      end if
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
				    leftv(1:4)=U_C(IBOUND(N,i)%which)%VAL(5,1:4)
				    CALL CONS2PRIM(N)
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
					ICONSI=IBOUND(N,i)%which
					 ICONSIDERED=IBOUND(N,i)%FACE
					select case(kkd)
					 case(1)
					 
					 call SHEAR_X2d_av(ICONSI,ICONSIDERED)
					 CASE (2)
					 call SHEAR_Y2d_av(ICONSI,ICONSIDERED)
					 
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

  call OUTWRITEPARA3Db

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


IF (TECPLOT.EQ.2)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEPARA3Dsb

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
					call OUTWRITE3vb
					
					
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

  call OUTWRITEPARA3Db

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
IF (TECPLOT.EQ.2)THEN		!BINARY PARAVIEW 3D ONLY

  call OUTWRITEPARA3Dsb

END IF



END SUBROUTINE surface_SOLUTION_WRITE

SUBROUTINE VOLUME_SOLUTION_WRITE_av
!> @brief
!> This subroutine calls the appropriate average volume writing subroutine based on the settings
IMPLICIT NONE

IF (TECPLOT.EQ.1)THEN
   if (dimensiona.eq.3)then

call OUTWRITE3vbav
eLSE

call OUTWRITE3vb2Dav
END IF
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
IF (TECPLOT.EQ.2)THEN

CALL OUTWRITEPARA3Dsbav

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

ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)))	!I ALLOCATE IN MEMORY THE PATTERN OF ACCESS OF DATA IN TERMS OF DISPLACEMENT, AND IN TERMS OF BLOCKLENGTH, AND FINALY AN ARRAY WITH THIS PROCESSOR DATA

      DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+turbulenceequations+passivescalar)
      END DO

      n_end=NOF_VARIABLES+turbulenceequations+passivescalar

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
    WRITE(1086)XBIN(i,1:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
    END DO

    DEALLOCATE(XBIN,ICELLA,VALUESA)
    close(1086)
    END IF
    
DEALLOCATE(VALUESS)


! !NEW APPROACH WITH GATHERV
!     IF (N.EQ.0) then
!     ALLOCATE(VALUESA(ImaxE));VALUESA=ZERO
!     allocate(xbin(imaxe,5+turbulenceequations+passivescalar))
!     END IF
!     
!     
!     ALLOCATE(VALUESS(KmaxE));VALUESS=ZERO
! 
! 
! 
! 
! IF (TURBULENCE.EQ.1)THEN
! do jj=1,5+turbulenceequations+passivescalar
!  
!       IF (jj.gt.5) THEN
!       DO I=1,KMAXE
! 		      VALUESS(I)=U_CT(I)%VAL(1,jj-5)
!       END DO
!       Else
!       DO I=1,KMAXE
! 	VALUESS(I)=U_C(I)%VAL(1,jj)
!       END DO
!       end if
! 
! !       CALL MPI_GATHERV(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNTS,DISPLS, RECVTYPE, ROOT, COMM, IERROR)
!       
!       
!     call MPI_GATHERV(VALUESS,KMAXE,MPI_DOUBLE_PRECISION,VALUESA,imaxE,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
!     IF (N.EQ.0)THEN
!     do i=1,imaxp*isize
! 	if (icella(i).gt.0)then
! 	xbin(icella(i),jj)=valuesa(i)
! 	end if
!     end do
!     end if
!     
! end do
! ELSE
!   do jj=1,5
! 	DO I=1,KMAXE
! 		VALUESS(I)=U_C(I)%VAL(1,jj)
! 	END DO
! 	call MPI_GATHER(VALUESS,imaxp,MPI_DOUBLE_PRECISION,VALUESA,imaxp,mpi_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IERROR)
! 	
! 	    IF (N.EQ.0)THEN
! 	    do i=1,imaxp*isize
! 		if (icella(i).gt.0)then
! 		xbin(icella(i),jj)=valuesa(i)
! 		end if
! 	    end do
! 	    end if
!     
!   end do
! 
! END IF
! 
! !   CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 
!     IF (N.EQ.0)THEN
! 
!     DO I=1,IMAXE
!     WRITE(1086)XBIN(i,1:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
!     END DO
! 
!     DEALLOCATE(XBIN,ICELLA,VALUESA)
!     close(1086)
!     END IF
!     
! DEALLOCATE(VALUESS)





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

ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)))	!I ALLOCATE IN MEMORY THE PATTERN OF ACCESS OF DATA IN TERMS OF DISPLACEMENT, AND IN TERMS OF BLOCKLENGTH, AND FINALY AN ARRAY WITH THIS PROCESSOR DATA

      DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+turbulenceequations+passivescalar)
      END DO

      n_end=NOF_VARIABLES+turbulenceequations+passivescalar

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
INTEGER:: prev_turbequation,INITIAL,III,i,k,jx,QQP,INC,kmaxe,jkn,ki,iterr,JX2,ind1,fh,size_of_real,size_of_int,dip,N_END,datatype
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

ALLOCATE(DISPT(KMAXE),ARRAY2(KMAXE*(NOF_VARIABLES+turbulenceequations+passivescalar)))
    DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+turbulenceequations+passivescalar)
      END DO
      
      
      n_end=NOF_VARIABLES+prev_turbequation+LAMPS
      
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
		  call MPI_file_READ(fh, INITIALRES(1:nof_Variables+turbulenceequations),nof_Variables+turbulenceequations, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_Variables+turbulenceequations)	!3
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
				U_CT(KI)%VAL(1,1)=VISC*TURBINIT
			    ELSE
				U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
				U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
					/L_TURB_INLET*RG(1)
			    END IF
			ENDIF
		K=K+prev_turbmodel+LAMPS
		
	    END DO
	END IF
	call MPI_Barrier(MPI_COMM_WORLD, ierror)
	DEALLOCATE(DISPT,ARRAY2)
	
	
      
if (Averaging .EQ. 1) then
 
IF (Average_restart.EQ.1)THEN

disp_in_file=0
ALLOCATE(DISPT(KMAXE),ARRAY2(kmaxe*(NOF_VARIABLES+turbulenceequations+passivescalar+6+passivescalar)))
DO I=1,KMAXE
	DISPT(I)=(XGO(I)-1)*(NOF_VARIABLES+turbulenceequations+passivescalar+6+passivescalar)
      END DO

      RESTFILE='RESTARTav.dat'
      n_end=NOF_VARIABLES+turbulenceequations+passivescalar+6+passivescalar
      
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


!INFO_NULL, fh, ierror)
! 
! 
! 	  
! 	 if (IRES_UNSTEADY .lt.1) then
! 	    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 	    call MPI_file_write(fh, RESCOUNTER, 1, MPI_INTEGER, MPI_STATUS_IGNORE,ierror)
! 			disp_in_file = disp_in_file + size_of_int	!1
! 	    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 	    call MPI_file_write(fh, INITIALRES(nof_variables+prev_turbmodel+PASSIVESCALAR), nof_variables+prev_turbmodel+PASSIVESCALAR, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
! 			disp_in_file = disp_in_file + size_of_real*(Nof_variables+prev_turbmodel+PASSIVESCALAR)	!1	    
! 	 else
! 	    if (initcond.eq.95)then
! 		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 		    call MPI_file_write(fh, iterr, 1, MPI_INTEGER, MPI_STATUS_IGNORE,ierror)
! 		      disp_in_file = disp_in_file + size_of_int	!1
! 		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 		    call MPI_file_write(fh, RES_TIME, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
! 		      disp_in_file = disp_in_file + size_of_real   
! 		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 		    call MPI_file_write(fh, taylor, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
! 		      disp_in_file = disp_in_file + size_of_real   
! 		
! 	      eLSE
! 		  call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 		    call MPI_file_write(fh, iterr, 1, MPI_INTEGER, MPI_STATUS_IGNORE,ierror)
! 		      disp_in_file = disp_in_file + size_of_int	!1
! 		    call MPI_file_seek(fh, disp_in_file, MPI_SEEK_SET, ierror)
! 		    call MPI_file_write(fh, RES_TIME, 1, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE,ierror)
! 		      disp_in_file = disp_in_file + size_of_real   
! 	      END IF
! 	 end if
!       
!       
!       call MPI_Barrier(MPI_COMM_WORLD, ierror)
! 	  
! 	  
! 	  if (n.ne.0) then
! 	    do i=0, n-1
! 	      tmp = tmp + (xmpiall(i)*size_of_int) + (xmpiall(i)*(size_of_real*(nof_Variables+prev_turbequation+lamps)))
! 	    end do
! 	    disp_in_file = disp_in_file + tmp
! 	  end if
! 	  write(100+n,*)disp_in_file,xmpiall(n)
! 
! 
! 
!   DO I=1,IMAXE
!   
! 	READ(1083)KI,RG(1:NOF_VARIABLES+prev_turbmodel+LAMPS)
! 	
! 	IF (XMPIE(I).EQ.N)THEN
! 	    KI=XMPIL(I)
! 	    U_C(KI)%VAL(1,1:NOF_VARIABLES)=RG(1:NOF_VARIABLES)
! 	 
! 		  IF (IRES_TURB.LT.1)THEN	!FROM NON TURBULENT RESTART	
! 		      IF (TURBULENCE.EQ.1)THEN
! 			  IF (TURBULENCEMODEL.EQ.1)THEN
! 			      U_CT(KI)%VAL(1,1)=VISC*TURBINIT
! 			  ELSE
! 			      U_CT(KI)%VAL(1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 			      U_CT(KI)%VAL(1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(1,1))&
! 				      /L_TURB_INLET*RG(1)
! 			  END IF
! 		      ENDIF
! 		  
! 		  ELSE			!FROM TURBULENT RESTART
! 		  
! 		      if (turbulencemodel .eq.prev_turbmodel)then
! 			      U_CT(KI)%VAL(1,1:Turbulenceequations)=RG(nof_Variables+1:nof_variables+Turbulenceequations)
! 		      end if
! 		  END IF
! 		  
! 		  
! 		  IF (LAMPS.LT.1)THEN
! 			  IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=ZERO
! 			  END IF
! 		  ELSE
! 			  IF (PASSIVESCALAR.GT.0)THEN
! 				  U_CT(KI)%VAL(1,TURBULENCEEQUATIONS+1)=RG(NOF_VARIABLES+prev_turbmodel+LAMPS)
! 			  END IF
! 		  END IF
! 		  
! 		  
! 	  END IF
!   END DO
! 
!   CLOSE(1083)
! 
! if (Averaging .EQ. 1) then
! 
! IF (Average_restart.EQ.1)THEN
! 
! allocate(Arg(nof_variables+prev_turbmodel+lamps+6))
! 
! RESTFILE='RESTARTav.dat'
! 
! OPEN(1084,FILE=RESTFILE,FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')
! 	  
! 	    DO I=1,IMAXE
!   
! 		  READ(1084)KI,aRG(1:NOF_VARIABLES+prev_turbmodel+LAMPS)
! 	
! 		  IF (XMPIE(I).EQ.N)THEN
! 		      KI=XMPIL(I)
! 			U_C(KI)%VAL(ind1,1:NOF_VARIABLES)=aRG(1:NOF_VARIABLES)
! 	    
! 			IF (IRES_TURB.LT.1)THEN	!FROM NON TURBULENT RESTART	
! 			  IF (TURBULENCE.EQ.1)THEN
! 			      IF (TURBULENCEMODEL.EQ.1)THEN
! 				  U_CT(KI)%VAL(ind1,1)=VISC*TURBINIT
! 			      ELSE
! 				  U_CT(KI)%VAL(ind1,1)=1.5*(I_turb_inlet*ufreestream)**2
! 				  U_CT(KI)%VAL(ind1,2)=(C_MU_INLET**(-0.25))*SQRT(U_CT(KI)%VAL(ind1,1))&
! 					  /L_TURB_INLET*aRG(1)
! 			      END IF
! 			  ENDIF
! 		      
! 			ELSE			!FROM TURBULENT RESTART
! 		      
! 			  if (turbulencemodel .eq.prev_turbmodel)then
! 				  U_CT(KI)%VAL(ind1,1:Turbulenceequations)=aRG(nof_Variables+1:nof_variables+Turbulenceequations)
! 			  end if
! 		      END IF
! 		  
! 		  
! 			IF (LAMPS.LT.1)THEN
! 				IF (PASSIVESCALAR.GT.0)THEN
! 					U_CT(KI)%VAL(ind1,TURBULENCEEQUATIONS+1)=ZERO
! 				END IF
! 			ELSE
! 				IF (PASSIVESCALAR.GT.0)THEN
! 					U_CT(KI)%VAL(ind1,TURBULENCEEQUATIONS+1)=aRG(NOF_VARIABLES+prev_turbmodel+LAMPS)
! 				END IF
! 			END IF
! 		  
! 		  U_C(KI)%RMS(1:6+lamps)=ARG(5+prev_turbmodel+lamps+1:11+prev_turbmodel+lamps+lamps)
! 		  
! 		END IF
! 		
! 		
! 		
! 	   END DO
! 
! 	   CLOSE(1084)
! 	   
!   DEALLOCATE(ARG)
!   
!   
!   
!   ELSE
! 	    
!   DO I=1,kmaxe
!       U_C(i)%VAL(ind1,:)=ZERO
!       U_C(i)%RMS(:)=ZERO
! 	if ((passivescalar.gt.0).or.(turbulence.eq.1))then
! 	U_CT(i)%VAL(ind1,:)=ZERO
! 	end if   
!   END DO
!   
! 
! 
!   END IF
! 
! 
! 
! 
!  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
! 
! deallocate(rg)
! 
! 










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
    WRITE(1086)XBIN(xmpi_re(i),1:4+TURBULENCEEQUATIONS+PASSIVESCALAR+3+PASSIVESCALAR)
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
				CALL CONS2PRIM(N)
	WRITE(3000+N,'(1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5)
				ELSE
				LEFTV(1:NOF_vARIABLES)=U_C(PROBEI(N,INV))%VAL(1,1:NOF_vARIABLES)
				CALL CONS2PRIM(N)
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
CHARACTER(LEN=12)::PROC,RESTFILE,PROC3
REAL::DRAG,LIFT,CD,CL,RX,PX,EX,surface_temp,RTEMP
REAL::FORCEXFR,SSX,CDF,LIFTF,DRAGF,FRICTIONF,TAUYX,TAUZX,TAUZY,SSY,SSZ,CF,TAUXX,TAUYY,TAUZZ
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ
 REAL,DIMENSION(2)::CI,CO
 logical::heref
 INTEGER::IM
 REAL::TSOLR,TSOLU,TSOLE,TSOLV,TSOLW,TSOLP,SSP
FORCEX=zero; FORCEY=zero; FORCEZ=zero;  FORCEXFR=zero
 CD=zero
 CL=zero
 CI(:)=zero
 CO(:)=zero
 KMAXE=XMPIELRANK(N)
 
!$OMP BARRIER 
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:FORCEX,FORCEY,FORCEZ)
DO I=1,kmaxe
		if (ielem(n,i)%interior.eq.1)then	
		    do j=1,ielem(n,i)%ifca
		      if (ielem(n,i)%ibounds(j).gt.0)then
			  if ((ibound(n,ielem(n,i)%ibounds(j))%icode.eq.4).AND.(inoder(IELEM(N,I)%NODES_FACES(J,1))%CORD(3).GT.0.001))then
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
					VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO
					  call  QUADRATUREQUAD3D(N,IGQRULES)
					  end if
					  surface_temp=IELEM(N,I)%SURF(J)
					  
				  
				case(6)
					gqi_points=qp_triangle_n
					
					    
					if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					  NND=3
					do K=1,nnd
					  VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
					END DO
					call QUADRATURETRIANG(N,IGQRULES)
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
				  
				  
				  LEFTV(1:5)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  RIGHTV(1:5)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  
				  
				  
				    call cons2prim2(n)
				    px=leftv(5)
				    ssp=ssp+(px*WEQUA2D(im))
				    if (itestcase.eq.4)then
				    CALL SUTHERLAND(N,LEFTV,RIGHTV)
				  
				  
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
	FORCEZ=FORCEZ*VECTORZ
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
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ
 REAL,DIMENSION(2)::CI,CO
 logical::heref
 INTEGER::IM
 REAL::TSOLR,TSOLU,TSOLE,TSOLV,TSOLW,TSOLP,SSP
FORCEX=zero; FORCEY=zero; FORCEZ=zero;  FORCEXFR=zero
 CD=zero
 CL=zero
 CI(:)=zero
 CO(:)=zero
 KMAXE=XMPIELRANK(N)
 
!$OMP BARRIER 
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:FORCEX,FORCEY,FORCEZ)
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
					VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO
					  
					  call  QUADRATURELINE(N,IGQRULES)
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
				  
				  LEFTV(1:4)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  RIGHTV(1:4)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				    call cons2prim2d2(n)
				    px=leftv(4)
				    
				    
				    if (itestcase.eq.4)then
				    CALL SUTHERLAND2d(N,LEFTV,RIGHTV)
				  
				  
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
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:ALLRES)
DO I=1,KMAXE
    ALLRES(1:5)=ALLRES(1:5)+(rhs(i)%VAL(1:5)**2)
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,5
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/IMAXE

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
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:ALLRES)
DO I=1,KMAXE
    ALLRES(1:5)=ALLRES(1:5)+(rhs(i)%VAL(1:5)**2)
    ALLRES(6:5+TURBULENCEEQUATIONS)=ALLRES(6:5+TURBULENCEEQUATIONS)+(RHST(I)%VAL(1:TURBULENCEEQUATIONS)**2)
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,7
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/IMAXE

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
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:ALLRES)
DO I=1,KMAXE
    ALLRES(1:4)=ALLRES(1:4)+(rhs(i)%VAL(1:4)**2)
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,4
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/IMAXE

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
!$OMP DO SCHEDULE(GUIDED) REDUCTION(+:ALLRES)
DO I=1,KMAXE
    ALLRES(1:4)=ALLRES(1:4)+(rhs(i)%VAL(1:4)**2)
    ALLRES(5:4+TURBULENCEEQUATIONS)=ALLRES(5:4+TURBULENCEEQUATIONS)+(RHST(I)%VAL(1:TURBULENCEEQUATIONS)**2)
END DO
!$OMP END DO

!$OMP MASTER
DO I=1,nof_variables+turbulenceequations
SUML3=ALLRES(I)
DUM_RESI=ZERO
CALL MPI_ALLREDUCE(SUML3,DUM_RESI,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
ALLRES(I)=DUM_RESI/IMAXE

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
	INTEGER::I,K,KMAXE
	REAL::EXACT,DUMMYOUT,DUMMYIN
	REAL::APROXIMATE
	real,dimension(15)::condm
	KMAXE=XMPIELRANK(N)
	L0NORM=ZERO;STENNORM=ZERO;L1NORM=ZERO
	  !$OMP MASTER
	  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	  !$OMP END MASTER
	  !$OMP BARRIER
	  
			!$OMP DO REDUCTION (+:L1NORM)
			DO I=1,KMAXE
				IF (ITESTCASE.Le.3)THEN
				EXACT=U_E(I)%VAL(1,1)
				APROXIMATE=U_C(I)%VAL(1,1)
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
				EXACT=U_E(I)%VAL(1,1)
				APROXIMATE=U_C(I)%VAL(1,1)
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
					
					IF (U_C(I)%VAL(1,1).GT.L0NORM)THEN
					L0NORM=U_C(I)%VAL(1,1)
					END IF
! 					L1NORM(N,1)=L1NORM(N,1)+((ABS(APROXIMATE-EXACT)))
				
 			END DO
 			!$OMP END DO 
 			!$OMP DO REDUCTION (MIN:L1NORM)
			DO I=1,KMAXE
					
					IF (U_C(I)%VAL(1,1).LT.L1NORM)THEN
					L1NORM=U_C(I)%VAL(1,1)
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
			WRITE(30,'(I9,1X,I4,1X,I4,1X,E14.7,1X,E14.7)')IMAXE,iorder,spatiladiscret,L0NORM,STENNORM/IMAXE
			
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
		  leftv(1:5)=U_C(I)%VAL(1,1:5)
		  CALL CONS2PRIM(N)
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
		  leftv(1:5)=U_C(I)%VAL(1,1:5)
		  CALL CONS2PRIM(N)
		  VALUESS(i)=leftv(5)
		END DO
    
    end if
    
		
    
    
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
    WRITE(400+N)"LOOKUP_TABLE default"//lf
    
    
		do i=1,imaxe
		xbin(XMPI_RE(I))=xbin2(I)
		end do
    
    
     WRITE(400+N)xbin(1:imaxe)
    END IF

     
    
    
end do
    
    
    
    

     
    
  
    
  IF (N.EQ.0)THEN
  CLOSE(400+N) 
   
   DEALLOCATE(XBIN,VALUELOCATION,XBIN2)
  deallocate(out1)
  END IF
  DEALLOCATE (VALUESS)
  

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

        
 deallocate(variables)


	
	

END SUBROUTINE OUTWRITEPARA3Db

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
						leftv(1:5)=U_C(IBOUND_T(I))%VAL(1,1:5)
						CALL CONS2PRIM(N)
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
		  leftv(1:5)=U_C(I)%VAL(IND1,1:5)
		  CALL CONS2PRIM(N)
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
						leftv(1:5)=U_C(IBOUND_T(I))%VAL(IND1,1:5)
						CALL CONS2PRIM(N)
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

















END MODULE IO
