MODULE TRANSLATe
USE DECLARATION
implicit none
integer,allocatable,dimension(:):: interray


contains

SUBROUTINE TRANSLATE_mesh
!> @brief
!> subroutine for transforming fluent style msh file to native format
implicit none
LOGICAL::HEREs,HEREF,HEREU
 CHARACTER(LEN=20)::PROC,UCNS3DFILE,FLUENTFILE,UGRIDFILE
 	UCNS3DFILE='GRID.bnd'
 	FLUENTFILE='grid.msh'
 	UGRIDFILE='grid.ugrid'
	
	INQUIRE (FILE=UCNS3DFILE,EXIST=HEREs)
	IF (HEREs) THEN
	
	!PROCEED WITHOUT TRANSLATION


	
	else


			INQUIRE (FILE=FLUENTFILE,EXIST=HEREF)


			IF (HEREF)THEN

			call Drive(interray)

			END IF



			INQUIRE (FILE=UGRIDFILE,EXIST=HEREU)


			IF (HEREU)THEN

			call TRANSUGRID

			END IF


	END IF



END SUBROUTINE TRANSLATE_mesh


!!! http://people.sc.fsu.edu/~jburkardt/f_src/chrpak/chrpak.html
Subroutine Drive(interray)
implicit none


TYPE(A_ELEMENT_NUMBER1),ALLOCATABLE,DIMENSION(:)::IELE
TYPE(A_NODE_NUMBER1),ALLOCATABLE,DIMENSION(:)::INOD
TYPE(A_BOUNDARY_NUMBER1),ALLOCATABLE,DIMENSION(:)::IBOU
TYPE(A_FACE_Number),ALLOCATABLE,DIMENSION(:)::IFAC

Integer::ing2,jj,dimen,imaxe,imaxn,imaxb,dum,index10,zoneid,in1,dum1,bctypdum,iosx,ii,vrt1,vrt2,vrt3,vrt4,vct1,vct2,vct3,vct4,countb,eltype,ing,ibtr,nfin,icte
Integer::lexist,checkbrac,dum2,dum3,dum4,ierr,il,il1,il2,ispace1,ispace2
Integer :: countline,counto,countword,lengt,intsize,imaxnglobal,imaxeglobal,imin,imax,imaxfglobal,index1,ina,in2,ichen,iix,iiy,icountfc,corn,icorn
integer,dimension(4)::nod,nodx,cans,xcand,cane,cang,canh,iiz,canf
integer,allocatable,dimension(:),intent(inout):: interray
integer,allocatable,dimension(:)::spaces,integerarray,ishape
integer,allocatable,dimension(:,:)::ifaci,iface
character(len=3) :: checkcomm,checkdim,checknod,checkelety,checkfac,checkbcne
integer::IOS,dumhex,ioss,endline,str2int
character(len=1)::braco,bracstr
character(len=3)::dumc,bracc
character(len=3)::dumc1
character(len=2)::dumc2
real(8),allocatable,dimension(:)::x,y,z
logical:: back
character(256) :: gchar,char1
character(len=256)::gchar1,gchar2
character:: ch
integer:: i
ios=0;ioss=0
 endline=0


  

OPEN(82,FILE="grid.msh",FORM='formatted',STATUS='OLD',ACTION='READ',IOSTAT=ios)
 countline=0
do 
	read(82,"(A3)",advance='NO',IOSTAT=ios)dumc
	dumc2=dumc(2:3)
	read(dumc2,*,iostat=ioss)str2int
	if ((ios .eq. -1))goto 11
	countline=countline+1
! 	print*,'line',countlieq. "(0")ne !,dumc
	     ! if (dumc .eq. "(0")  then ! COMMENT CONDITIONS
	      if (dumc .eq. "(0 ")  then ! COMMENT CONDITIONS
		  read(82,*,IOSTAT=ios)
		  countline=countline+1
		  go to 10
	      Endif
! 	      if (dumc .eq. "(2") then
! 	      read(82,*,IOSTAT=IOS)dimen
	       if (dumc .eq. "(2 ") then
       read(82,'(I1)',IOSTAT=IOS) dimen !changed by holger foysi
		  countline=countline+1
! 	      print*,'Dimension:',dimen
	      end if
	      if (dumc .eq. "(10") then
		      read(82,'(A)',iostat=ios) gchar
			countline=countline+1
			
		      call removebrac(gchar,gchar2)
! 		      print*,gchar2!,il
		      call string2int (gchar2,interray,intsize)
		      if (interray(1) .eq. 0)then
! 			      do i=1,intsize
! 				  print*,interray(i)
! 			      enddo
			  imaxnglobal=interray(3) ! set global maximum number of nodes
			  allocate(x(imaxnglobal));allocate(y(imaxnglobal))
! 			  print*,'Max Number Nodes:',imaxnglobal
			  if (dimen.eq.3) allocate(z(imaxnglobal))
		      endif
		      if (interray(1) .ne. 0) then 
			read(82,"(A1)",advance='NO',IOSTAT=ios)dumc
			countline=countline+1
			  if (dumc .eq. "(") then
			      Do i=interray(2),interray(3)
				if (dimen.eq.2) then
				read(82,*,IOSTAT=ios)x(i),y(i)
				countline=countline+1
				endif
				if (dimen.eq.3) then
				read(82,*,IOSTAT=ios)x(i),y(i),z(i)
				countline=countline+1
				endif
			      enddo
			  else
			    backspace(82,iostat=ios)
			  Do i=interray(2),interray(3)
				if (dimen.eq.2) then
				read(82,*,IOSTAT=ios)x(i),y(i)
				countline=countline+1
				endif
				if (dimen.eq.3) then
				read (82,*,IOSTAT=ios)x(i),y(i),z(i)
				countline=countline+1
				endif
			      enddo
			  endif
			  if (binio.eq.0)OPEN(10,FILE="GRID.vrt",FORM='formatted',ACTION='WRITE',IOSTAT=iosx,position='append')
			  if (binio.eq.1)OPEN(10,FILE="GRID.vrt",FORM='unformatted',ACTION='WRITE',IOSTAT=iosx,position='append')
			  selectcase (dimen)
			    case(2)
			    Do i=interray(2),interray(3)
			      if (binio.eq.0)write(10,"(5X, I8, 2X,ES21.14,2X,ES21.14)")i,x(i),y(i)
			      if (binio.eq.1)write(10)i,x(i),y(i)
			    end do
			    close(10)
			    case(3)
			    Do i=interray(2),interray(3)
			      if (binio.eq.0)write(10,"(5X,I8,2X,ES21.14,2X,ES21.14,2X,ES21.14)")i,x(i),y(i),z(i)
			      if (binio.eq.1)write(10)i,x(i),y(i),z(i)
			    end do 
			    close(10)
			    deallocate(x,y)
                            if (dimen.eq.3) deallocate(z)
!                             
			  end select
		      endif
		deallocate(interray)
		
! 			  print*,'Max Number Nodes:',imaxnglobal
! 			  if (dimen.eq.3) allocate(z(imaxnglobal))
	      end if ! (10
	      if (dumc .eq. "(12") then
			  read(82,'(A)',iostat=ios) gchar
			  countline=countline+1
			  call removebrac(gchar,gchar2)
	  ! 		print*,gchar2,"edw1"!,il
			  call string2int (gchar2,interray,intsize)
			  if (interray(1) .eq. 0)then
			      imaxeglobal=interray(3) ! set global maximum number of elements
			      allocate(IELE(imaxeglobal))
			      allocate(ishape(imaxeglobal))
! 			      print*,'Max Number Cells:',imaxeglobal
			      iele(1:imaxeglobal)%IECOUNTER=0
			  endif

!                         
			   if ((interray(1) .ne. 0))then
			    if (interray(5).ne.0) then 
			      IELE(interray(2):interray(3))%ishape=interray(5)
			       Do i=interray(2),interray(3)
! 				      
				      iele(i)%ieindex=i
				       SELECT CASE (IELE(I)%ISHAPE)
					  Case(1)
					    iele(i)%iface=3
					    allocate(iele(i)%faces(1:3,1:2))
					    !allocate(iele(i)%nd(1:44))
					  Case(3)
					    iele(i)%iface=4
					    allocate(iele(i)%faces(1:4,1:2))
					    !allocate(iele(i)%nd(4))
					  Case(2)
					    iele(i)%iface=4
					    allocate(iele(i)%faces(1:4,1:3))
! 					    allocate(iele(i)%nd(8))
					  Case(4)
					    iele(i)%iface=6
					    allocate(iele(i)%faces(1:6,1:4))
! 					    allocate(iele(i)%nd(8))
					  Case(5)
					    iele(i)%iface=5
					    allocate(iele(i)%faces(1:5,1:4))
! 					    allocate(iele(i)%nd(8))
					  Case(6)
					    iele(i)%iface=5
					    allocate(iele(i)%faces(1:5,1:4))
! 					    allocate(iele(i)%nd(8))

				    End select

				      iele(i)%faces(:,:)=0
				  end do
			  end if
			  end if
			  if ((interray(1) .ne. 0))then
			    if (interray(5).eq.0) then 

			    read(82,"(A1)",advance='NO',IOSTAT=ios)dumc
			    countline=countline+1
			      if (dumc .eq. "(") then
				  read(82,*,IOSTAT=ios)Ishape(interray(2):interray(3))!IELE(interray(2):interray(3))%ishape
				    IELE(interray(2):interray(3))%ishape=Ishape(interray(2):interray(3))
				  countline=countline+1
			      else
				backspace (82,iostat=ios)
			      read(82,*,IOSTAT=ios)Ishape(interray(2):interray(3))!IELE(interray(2):interray(3))%ishape
				    IELE(interray(2):interray(3))%ishape=Ishape(interray(2):interray(3))
			      countline=countline+1
			      endif
				  Do i=interray(2),interray(3)
! 				      
				      iele(i)%ieindex=i
				       SELECT CASE (IELE(I)%ISHAPE)
					  Case(1)
					    iele(i)%iface=3
					    allocate(iele(i)%faces(3,2))
					    
					  Case(3)
					    iele(i)%iface=4
					    allocate(iele(i)%faces(4,2))
					   
					  Case(2)
					    iele(i)%iface=4
					    allocate(iele(i)%faces(4,3))
					   
					  Case(4)
					    iele(i)%iface=6
					    allocate(iele(i)%faces(6,4))
					   
					  Case(5)
					    iele(i)%iface=5
					    allocate(iele(i)%faces(5,4))
					    
					  Case(6)
					    iele(i)%iface=5
					    allocate(iele(i)%faces(5,4))
					    

				    End select

				      iele(i)%faces(:,:)=0
				  end do
			  endif
			  end if
! 			  
			  deallocate(interray)
	      end if
	      if (dumc .eq. "(13") then
		read(82,'(A)',advance='no',iostat=ios) gchar
		countline=countline+1
		call removebrac(gchar,gchar2)
		 call string2int (gchar2,interray,intsize)
		 index1=interray(1);imin=interray(2);imax=interray(3);bctypdum=interray(4)
		  if (index1.ne.0)eltype=interray(5)

		 deallocate(interray)
		if (index1 .eq. 0)then
		    imaxfglobal=imax ! set global maximum number of faces
		    allocate(ifac(imaxfglobal))
		endif
		if (index1 .ne. 0) then 
		    
		  read(82,"(A1)",advance='NO',IOSTAT=ios)dumc
		  countline=countline+1
! 		  print*,'faces',imin,imax,index1,dumc
  		    if (dumc .eq. "(") then
!   		    
			   Do i=imin,imax
! 				if (dimen.eq.2) then
				  read(82,'(A)',advance='no',iostat=ios) gchar
				  countline=countline+1
! 				  print*,gchar,'ssss'
				  call string2int (gchar,interray,intsize)
! 				  
				  select case (dimen)
				    case (2)! 2d
					     if (eltype.eq.0)then
					ifac(i)%ishb=2
					ifac(i)%IFACBTYPE=bctypdum
					allocate(ifac(i)%ifa(2,2))
					ifac(i)%ifa(1,1)=interray(2)
					ifac(i)%ifa(1,2)=interray(3)
					ifac(i)%ifa(2,1)=interray(4)
					ifac(i)%ifa(2,2)=interray(5)
! 					if ( ifac(i)%ifa(2,1) .ne. 0) then ! 1st element
! 					Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					end if
! 					if ( ifac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! 					Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,2)
! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,1)
! 					endif
					else
					 ifac(i)%ishb=2
					    ifac(i)%IFACBTYPE=bctypdum
					allocate(ifac(i)%ifa(2,2))
					ifac(i)%ifa(1,1)=interray(1)
					ifac(i)%ifa(1,2)=interray(2)
					ifac(i)%ifa(2,1)=interray(3)
					ifac(i)%ifa(2,2)=interray(4)
! 					if ( ifac(i)%ifa(2,1) .ne. 0) then ! 1st element
! 					Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					end if
! 					if ( ifac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! 					Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,2)
! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,1)
! 					endif


					    end if

					
				    case (3)!3d
					ifac(i)%IFACBTYPE=bctypdum

					  if (eltype.eq.0)then
					ifac(i)%ishb=Interray(1)
					select case (Interray(1))
					  case (3)! triangles
					      allocate(ifac(i)%ifa(2,3))
					      ifac(i)%ifa(1,1)=interray(2)
					      ifac(i)%ifa(1,2)=interray(3)
					      ifac(i)%ifa(1,3)=interray(4)
					      ifac(i)%ifa(2,1)=interray(5)
					      ifac(i)%ifa(2,2)=interray(6)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then ! no wall 
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      allocate(ifac(i)%ifa(2,4))
					      ifac(i)%ifa(1,1)=interray(2)
					      ifac(i)%ifa(1,2)=interray(3)
					      ifac(i)%ifa(1,3)=interray(4)
					      ifac(i)%ifa(1,4)=interray(5)
					      ifac(i)%ifa(2,1)=interray(6)
					      ifac(i)%ifa(2,2)=interray(7)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,4)=ifac(i)%ifa(1,4)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,4)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,4)=ifac(i)%ifa(1,1)
! 					      endif

					  end select
					    else
					      ifac(i)%ishb=eltype
					select case (eltype)
					  case (3)! triangles
					      allocate(ifac(i)%ifa(2,3))
					      ifac(i)%ifa(1,1)=interray(1)
					      ifac(i)%ifa(1,2)=interray(2)
					      ifac(i)%ifa(1,3)=interray(3)
					      ifac(i)%ifa(2,1)=interray(4)
					      ifac(i)%ifa(2,2)=interray(5)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then ! no wall 
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      allocate(ifac(i)%ifa(2,4))
					      ifac(i)%ifa(1,1)=interray(1)
					      ifac(i)%ifa(1,2)=interray(2)
					      ifac(i)%ifa(1,3)=interray(3)
					      ifac(i)%ifa(1,4)=interray(4)
					      ifac(i)%ifa(2,1)=interray(5)
					      ifac(i)%ifa(2,2)=interray(6)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,4)=ifac(i)%ifa(1,4)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,4)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,4)=ifac(i)%ifa(1,1)
! 					      endif



					  endselect ! quad or triangles faces (3d)
		    
					end if
					endselect ! dimensions
					
! 				endif
! 				if (dimen.eq.3) then
! 				read (82,*)x(i),y(i),z(i)
! 				endif
				  deallocate(interray)
			      enddo
		    else
!                       
		      backspace (82,iostat=ios)
     			   Do i=imin,imax
! 				if (dimen.eq.2) then
				  read(82,'(A)',advance='no',iostat=ios) gchar
				  countline=countline+1
! 				  print*,gchar,'ssss'
				  call string2int (gchar,interray,intsize)
! 				  
				  select case (dimen)
				    case (2)! 2d
					ifac(i)%IFACBTYPE=bctypdum
					     if (eltype.eq.0)then
 					ifac(i)%ishb=2
					ifac(i)%IFACBTYPE=bctypdum
					allocate(ifac(i)%ifa(2,2))
					ifac(i)%ifa(1,1)=interray(2)
					ifac(i)%ifa(1,2)=interray(3)
					ifac(i)%ifa(2,1)=interray(4)
					ifac(i)%ifa(2,2)=interray(5)
! 					if ( ifac(i)%ifa(2,1) .ne. 0) then ! 1st element
! ! 					Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! ! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! ! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! ! 					end if
! ! 					if ( ifac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! ! 					Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! ! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,2)
! ! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,1)
! ! 					endif
					else
					 ifac(i)%ishb=2
					    ifac(i)%IFACBTYPE=bctypdum
					allocate(ifac(i)%ifa(2,2))
					ifac(i)%ifa(1,1)=interray(1)
					ifac(i)%ifa(1,2)=interray(2)
					ifac(i)%ifa(2,1)=interray(3)
					ifac(i)%ifa(2,2)=interray(4)
! 					if ( ifac(i)%ifa(2,1) .ne. 0) then ! 1st element
! 					Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					end if
! 					if ( ifac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! 					Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,2)
! 					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,1)
! 					endif


					    end if
				    case (3)!3d
					ifac(i)%IFACBTYPE=bctypdum

					  if (eltype.eq.0)then
					ifac(i)%ishb=Interray(1)
					select case (Interray(1))
					  case (3)! triangles
					      allocate(ifac(i)%ifa(2,3))
					      ifac(i)%ifa(1,1)=interray(2)
					      ifac(i)%ifa(1,2)=interray(3)
					      ifac(i)%ifa(1,3)=interray(4)
					      ifac(i)%ifa(2,1)=interray(5)
					      ifac(i)%ifa(2,2)=interray(6)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then ! no wall 
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      allocate(ifac(i)%ifa(2,4))
					      ifac(i)%ifa(1,1)=interray(2)
					      ifac(i)%ifa(1,2)=interray(3)
					      ifac(i)%ifa(1,3)=interray(4)
					      ifac(i)%ifa(1,4)=interray(5)
					      ifac(i)%ifa(2,1)=interray(6)
					      ifac(i)%ifa(2,2)=interray(7)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,4)=ifac(i)%ifa(1,4)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,4)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,4)=ifac(i)%ifa(1,1)
! 					      endif
					    end select
					    else
					      ifac(i)%ishb=eltype
					select case (eltype)
					  case (3)! triangles
						ifac(i)%ishb=eltype
					      allocate(ifac(i)%ifa(2,3))
					      ifac(i)%ifa(1,1)=interray(1)
					      ifac(i)%ifa(1,2)=interray(2)
					      ifac(i)%ifa(1,3)=interray(3)
					      ifac(i)%ifa(2,1)=interray(4)
					      ifac(i)%ifa(2,2)=interray(5)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then ! no wall 
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      ifac(i)%ishb=eltype
					      allocate(ifac(i)%ifa(2,4))
					      ifac(i)%ifa(1,1)=interray(1)
					      ifac(i)%ifa(1,2)=interray(2)
					      ifac(i)%ifa(1,3)=interray(3)
					      ifac(i)%ifa(1,4)=interray(4)
					      ifac(i)%ifa(2,1)=interray(5)
					      ifac(i)%ifa(2,2)=interray(6)
! 					      if ((ifac(i)%ifa(2,1)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,4)=ifac(i)%ifa(1,4)
! 					      endif
! 					      if ((ifac(i)%ifa(2,2)).ne.0)then
! 					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,4)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,3)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,2)
! 					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,4)=ifac(i)%ifa(1,1)
! 					      endif



					  endselect ! quad or triangles faces (3d)
					end if
					endselect ! dimensions
! 				endif
! 				if (dimen.eq.3) then
! 				read (82,*)x(i),y(i),z(i)
! 				endif
				  deallocate(interray)
			      enddo
		      endif ! brac
		  endif! index1
	      end if !(13 all faces
	      
! 	      if (dumc .eq. "(39")  return
      
 10 continue
!  print*,ios,'ios'
! read(*,*)
! if (ios .lt. 0) return
 
!  if (dunc .gt. 800) stop

end do

11 continue

close(82)
do i=1,imaxfglobal
! 					
					if  (ifac(i)%ishb.eq.2)then

					if ( ifac(i)%ifa(2,1) .ne. 0) then ! 1st element
					Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
					Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
					end if
					if ( ifac(i)%ifa(2,2) .ne. 0) then ! 2nd element
					Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,2)
					Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,1)
					endif
					end if
					if  (ifac(i)%ishb.eq.4)then
					 if ((ifac(i)%ifa(2,1)).ne.0)then
					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,4)=ifac(i)%ifa(1,4)
					      endif
					      if ((ifac(i)%ifa(2,2)).ne.0)then
					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,4)
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,3)
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,2)
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,4)=ifac(i)%ifa(1,1)
					      endif

					end if
					if  (ifac(i)%ishb.eq.3)then
					 if ((ifac(i)%ifa(2,1)).ne.0)then ! no wall 
					      Iele(ifac(i)%ifa(2,1))%IECOUNTER=Iele(ifac(i)%ifa(2,1))%IECOUNTER+1
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,1)=ifac(i)%ifa(1,1)
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,2)=ifac(i)%ifa(1,2)
					      Iele(ifac(i)%ifa(2,1))%FACES(Iele(ifac(i)%ifa(2,1))%IECOUNTER,3)=ifac(i)%ifa(1,3)
					      endif
					      if ((ifac(i)%ifa(2,2)).ne.0)then
					      Iele(ifac(i)%ifa(2,2))%IECOUNTER=Iele(ifac(i)%ifa(2,2))%IECOUNTER+1
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,1)=ifac(i)%ifa(1,3)
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,2)=ifac(i)%ifa(1,2)
					      Iele(ifac(i)%ifa(2,2))%FACES(Iele(ifac(i)%ifa(2,2))%IECOUNTER,3)=ifac(i)%ifa(1,1)
					      endif

				      end if
end do



IF (BINIO.EQ.0)THEN
OPEN(10,FILE="GRID.cel",FORM='formatted',ACTION='WRITE',IOSTAT=iosx)
ELSE
OPEN(10,FILE="GRID.cel",FORM='UNFORMATTED',ACTION='WRITE',IOSTAT=iosx)

END IF
! 
!  INTEGER,ALLOCATABLE,DIMENSION(:,:)::FACES ! id of face and id of node
!     Integer,allocatable,dimension(:):: ND
icte=0

Do i=1,imaxeglobal


      SELECT CASE(IELE(I)%ISHAPE)


      CASE(1)

	VRT1=IELE(I)%FACEs(1,1);VRT2=IELE(I)%FACEs(1,2)
	DO II=2,IELE(I)%IECOUNTER
	    VCT1=IELE(I)%FACEs(II,1)
	  IF (((VRT1.EQ.VCT1)).OR.((VRT2.EQ.VCT1)))THEN


	  ELSE
! 	  WRITE(10,"(5I10)")I,VRT2,VRT1,VCT1,VCT1
	  IF (BINIO.EQ.0)WRITE(10,"(5I10)")I,Vct1,VRT1,Vrt2,Vrt2
	  IF (BINIO.EQ.1)WRITE(10)I,Vct1,VRT1,Vrt2,Vrt2
	  CYCLE
	  END IF
	END DO



      CASE(3)
	VRT1=IELE(I)%FACEs(1,1);VRT2=IELE(I)%FACEs(1,2)
	DO II=2,IELE(I)%IECOUNTER
	    VCT1=IELE(I)%FACEs(II,1);VCT2=IELE(I)%FACEs(II,2)
	  IF (((VRT1.EQ.VCT1).OR.(VRT1.EQ.VCT2)).OR.((VRT2.EQ.VCT1).OR.(VRT2.EQ.VCT2)))THEN


	  ELSE
	  IF (BINIO.EQ.0)WRITE(10,"(5I10)")I,VRT1,VRT2,VCT1,VCT2
	  IF (BINIO.EQ.1)WRITE(10)I,VRT1,VRT2,VCT1,VCT2
	 
	  CYCLE
	  END IF
	END DO


      CASE(2)
	VRT1=IELE(I)%FACEs(1,1);VRT2=IELE(I)%FACEs(1,2);VRT3=IELE(I)%FACEs(1,3)

	nod(1)=vrt1
	nod(2)=vrt2
	nod(3)=vrt3

	DO II=2,IELE(I)%IECOUNTER
	       VCT1=IELE(I)%FACEs(II,1);VCT2=IELE(I)%FACEs(II,2);VCT3=IELE(I)%FACEs(II,3)
	      nodx(1)=vct1; nodx(2)=vct2; nodx(3)=vct3
	      do ina=1,3
		    ichen=0
		    do in2=1,3
		      if (nod(in2).ne.nodx(ina))then
			      ichen=ichen+1
		      end if
		    end do

		    if (ichen.eq.3)then
		    nod(4)=nodx(ina)
		
		    end if
  
	      end do
	 end do

	    
		    
! 		icte=icte+1
	   
   	       IF (BINIO.EQ.0)WRITE(10,"(9I10)")I,VRT1,VRT2,VRT3,vrt3,nod(4),nod(4),nod(4),nod(4)		!or
		IF (BINIO.EQ.1)WRITE(10)I,VRT1,VRT2,VRT3,vrt3,nod(4),nod(4),nod(4),nod(4)
		


	CASE(4)
		ing=0
		  VRT1=IELE(I)%FACEs(1,1);VRT2=IELE(I)%FACEs(1,2);VRT3=IELE(I)%FACEs(1,3);vrt4=IELE(I)%FACEs(1,4)
	DO II=2,IELE(I)%IECOUNTER
! 	    print*,"edw1",vrt1,vrt2,vrt3,ii,IELE(I)%IECOUNTER
	    VCT1=IELE(I)%FACEs(II,1);VCT2=IELE(I)%FACEs(II,2);VCT3=IELE(I)%FACEs(II,3);vct4=IELE(I)%FACEs(II,4)
	  IF( ((VRT1.EQ.VCT1).OR.(VRT2.EQ.VCT1).OR.(VRT3.EQ.VCT1).or.(vrt4.eq.vct1)).or.&
	      ((VRT1.EQ.VCT2).OR.(VRT2.EQ.VCT2).OR.(VRT3.EQ.VCT2).or.(vrt4.eq.vct2)).or.&
	      ((VRT1.EQ.VCT3).OR.(VRT2.EQ.VCT3).OR.(VRT3.EQ.VCT3).or.(vrt4.eq.vct3)).or.&
	      ((VRT1.EQ.VCT4).OR.(VRT2.EQ.VCT4).OR.(VRT3.EQ.VCT4).or.(vrt4.eq.vct4)))then
	    
		

	  ELSE





	 
!   		  WRITE(10,"(9I10)")I,VRT1,VRT2,VRT3,vrt4,VCT4,vct3,vct2,vct1

		  ing=ii


		    xcand(1)=vrt1
		    xcand(2)=vrt2
		    xcand(3)=vrt3
		    xcand(4)=vrt4


		    cane(1)=vct1
		    cane(2)=vct2
		    cane(3)=vct3
		    cane(4)=vct4

		    icountfc=0

		    do iix=2,IELE(I)%IECOUNTER
			if (iix.ne.ing)then
			  canf(1)=IELE(I)%FACEs(IIx,1)
			  canf(2)=IELE(I)%FACEs(IIx,2)
			  canf(3)=IELE(I)%FACEs(IIx,3)
			  canf(4)=IELE(I)%FACEs(IIx,4)
			do iiy=1,4
			      if (canf(iiy).eq.xcand(1))then

				icountfc=icountfc+1
				    iiz(icountfc)=iix
			       end if
			end do
			
			end if
		    end do

		    icountfc=0

		    do iix=2,IELE(I)%IECOUNTER
			 if ((iix.eq.iiz(1)).or.(iix.eq.iiz(2)))then
			      icountfc=icountfc+1
			      if (icountfc.eq.1)then
			       cang(1)=IELE(I)%FACEs(IIx,1)
			  cang(2)=IELE(I)%FACEs(IIx,2)
			  cang(3)=IELE(I)%FACEs(IIx,3)
			  cang(4)=IELE(I)%FACEs(IIx,4)
			      end if
			      if (icountfc.eq.2)then
			       canh(1)=IELE(I)%FACEs(IIx,1)
			  canh(2)=IELE(I)%FACEs(IIx,2)
			  canh(3)=IELE(I)%FACEs(IIx,3)
			  canh(4)=IELE(I)%FACEs(IIx,4)
			      end if
			 end if
		    end do
		      

		    do iix=1,4
			if (cang(iix).ne.xcand(1))then
		      do iiy=1,4
			  if (canh(iiy).ne.xcand(1))then
			      if (cang(iix).eq.canh(iiy))then
			      icorn=canh(iiy)

			      end if
			  end if
		      end do
			end if
		    end do

		     do iix=1,4

			  if (icorn.eq.cane(iix))then
			      iiy=iix




			  end if
		      end do
		      if (iiy.eq.1)then
			    cans(:)=cane(:)

		      end if
		      if (iiy.eq.2)then
			   cans(1)=cane(2)
			    cans(2)=cane(3)
			    cans(3)=cane(4)
			      cans(4)=cane(1)
			end if
			 if (iiy.eq.3)then
			    cans(1)=cane(3)
			    cans(2)=cane(4)
			    cans(3)=cane(1)
			      cans(4)=cane(2)

			end if
			if (iiy.eq.4)then
			    cans(1)=cane(4)
			    cans(2)=cane(1)
			    cans(3)=cane(2)
			      cans(4)=cane(3)

		      end if

		   
 		  IF (BINIO.EQ.0)WRITE(10,"(9I10)")I,xcand(1),xcand(2),xcand(3),xcand(4),cans(1),cans(4),cans(3),cans(2)
		  IF (BINIO.EQ.1)WRITE(10)I,xcand(1),xcand(2),xcand(3),xcand(4),cans(1),cans(4),cans(3),cans(2)
		cycle

	    end if
	  end do

	CASE(5)

	    	 
	ing=0
	DO II=1,IELE(I)%IECOUNTER

	   

	    if (IELE(I)%FACEs(II,4).ne.0)then
	    
	      ing=ii
	      VCT1=IELE(I)%FACEs(II,1);VCT2=IELE(I)%FACEs(II,2);VCT3=IELE(I)%FACEs(II,3);vct4=IELE(I)%FACEs(II,4)
	      nod(1)=vct1; nod(2)=vct2; nod(3)=vct3; nod(4)=vct4

	    end if

	end do
	DO II=1,IELE(I)%IECOUNTER

	     


		VRT1=IELE(I)%FACEs(ii,1);VRT2=IELE(I)%FACEs(ii,2);VRT3=IELE(I)%FACEs(ii,3)
	      nodx(1)=vrt1; nodx(2)=vrt2; nodx(3)=vrt3
	      do ina=1,3
		    ichen=0
		    do in2=1,4
		      if (nod(in2).ne.nodx(ina))then
			      ichen=ichen+1
		      end if
		    end do

		    if (ichen.eq.4)then
		    nfin=nodx(ina)
		
		    end if
  
	      end do

	  
	end do
! 		  icte=icte+1
 		  IF (BINIO.EQ.0)WRITE(10,"(9I10)")i,Vct1,Vct2,Vct3,vct4,nfin,nfin,nfin,nfin
 		  IF (BINIO.EQ.1)WRITE(10)i,Vct1,Vct2,Vct3,vct4,nfin,nfin,nfin,nfin
! 		  

	

	CASE(6)
	ing=0
	ing2=0
		      
	DO II=1,IELE(I)%IECOUNTER

! 	    if VRT1=IELE(I)%FACEs(1,1);VRT2=IELE(I)%FACEs(1,2);VRT3=IELE(I)%FACEs(1,3)

	     if (IELE(I)%FACEs(II,4).eq.0)then
	      
	      ing=ii
	      VRT1=IELE(I)%FACEs(ii,1);VRT2=IELE(I)%FACEs(ii,2);VRT3=IELE(I)%FACEs(ii,3)
 	      cycle

	    end if
	end do
	do ii=1,iele(i)%iecounter
	    if (ii.ne.ing)then

	     if (IELE(i)%FACEs(II,4).eq.0)then
	     VCT1=IELE(I)%FACEs(II,1);VCT2=IELE(I)%FACEs(II,2);VCT3=IELE(I)%FACEs(II,3)
	      ing2=ii
	      end if

	    end if


	end do

	     


		    xcand(1)=vrt1
		    xcand(2)=vrt2
		    xcand(3)=vrt3
 		    xcand(4)=0


		    cane(1)=vct1
		    cane(2)=vct2
		    cane(3)=vct3
 		    cane(4)=0

		    icountfc=0

		    do iix=1,IELE(I)%IECOUNTER
			if ((iix.ne.ing).and.(iix.ne.ing2))then
			  canf(1)=IELE(I)%FACEs(IIx,1)
			  canf(2)=IELE(I)%FACEs(IIx,2)
			  canf(3)=IELE(I)%FACEs(IIx,3)
			  canf(4)=IELE(I)%FACEs(IIx,4)
			do iiy=1,4
			      if (canf(iiy).eq.xcand(1))then

				icountfc=icountfc+1
				    iiz(icountfc)=iix
			       end if
			end do
			
			end if
		    end do

		    icountfc=0

		    do iix=1,IELE(I)%IECOUNTER
			 if ((iix.eq.iiz(1)).or.(iix.eq.iiz(2)))then
			      icountfc=icountfc+1
			      if (icountfc.eq.1)then
			       cang(1)=IELE(I)%FACEs(IIx,1)
			  cang(2)=IELE(I)%FACEs(IIx,2)
			  cang(3)=IELE(I)%FACEs(IIx,3)
			  cang(4)=IELE(I)%FACEs(IIx,4)
			      end if
			      if (icountfc.eq.2)then
			       canh(1)=IELE(I)%FACEs(IIx,1)
			  canh(2)=IELE(I)%FACEs(IIx,2)
			  canh(3)=IELE(I)%FACEs(IIx,3)
			  canh(4)=IELE(I)%FACEs(IIx,4)
			      end if
			 end if
		    end do
		      

		    do iix=1,4
			if (cang(iix).ne.xcand(1))then
		      do iiy=1,4
			  if (canh(iiy).ne.xcand(1))then
			      if ((cang(iix).ne.0).and.(cang(iix).eq.canh(iiy)))then
			      icorn=canh(iiy)

			      end if
			  end if
		      end do
			end if
		    end do

		     do iix=1,3

			  if (icorn.eq.cane(iix))then
			      iiy=iix




			  end if
		      end do
		      if (iiy.eq.1)then
			    cans(:)=cane(:)

		      end if
		      if (iiy.eq.2)then
			   cans(1)=cane(2)
			    cans(2)=cane(3)
			    cans(3)=cane(1)
! 			      cans(4)=cane(1)
			end if
			 if (iiy.eq.3)then
			    cans(1)=cane(3)
			    cans(2)=cane(1)
			    cans(3)=cane(2)
! 			      cans(4)=cane(2)

			end if
! 			if (iiy.eq.4)then
! 			    cans(1)=cane(4)
! 			    cans(2)=cane(1)
! 			    cans(3)=cane(2)
! 			      cans(4)=cane(3)
! 
! 		      end if



























  icte=icte+1









   IF (BINIO.EQ.0)WRITE(10,"(9I10)")I,xcand(1),xcand(2),xcand(3),xcand(3),cans(1),cans(3),cans(2),cans(2)
   IF (BINIO.EQ.1)WRITE(10)I,xcand(1),xcand(2),xcand(3),xcand(3),cans(1),cans(3),cans(2),cans(2)

		!WRITE(10,"(9I10)")I,VRT1,VRT2,VRT3,vrt3,VCT1,vct2,vct3,vct3
!  		WRITE(10,"(9I10)")Icte,VRT1,VRT2,VRT3,vrt3,VCT2,vct1,vct3,vct3
		
	      
! 	      



      ENDSELECT



	
end do

close(10) 


IF (BINIO.EQ.0)THEN
OPEN(10,FILE="GRID.bnd",FORM='formatted',ACTION='WRITE',IOSTAT=iosx)
ELSE
OPEN(10,FILE="GRID.bnd",FORM='unformatted',ACTION='WRITE',IOSTAT=iosx)
END IF
 countb=0
do i=1,imaxfglobal
	if (ifac(i)%IFACBTYPE.ne.2)then



	  select case(ifac(i)%ifacbtype)


	  case(7) !symmetry
	  ibtr=3
	
	  case(8,24,37)!periodicity	note that we use several codes
	  ibtr=5

	  case(3)  !wall
	  ibtr=4

	  case(36) !outflow

	  ibtr=2

	  case(10)!velocity inlet

	  ibtr=1

	  case(9)!pressure far field
	  ibtr=6

	   case(4)!pressure inlet

	  ibtr=1

	    case(5)!pressure outlet

	  ibtr=2

	  case(20)!mass flow inlet

	  ibtr=1

	  case(12)!mass flow inlet

	  ibtr=5



	end select



	  countb=countb+1


	  select case(dimen)

	  case(2)
	  IF (BINIO.EQ.0)write(10,"(6i12)")countb,ifac(i)%ifa(1,1),ifac(i)%ifa(1,2),0,0,ibtr
	  IF (BINIO.EQ.1)write(10)countb,ifac(i)%ifa(1,1),ifac(i)%ifa(1,2),0,0,ibtr				
	  case(3)

	    if (ifac(i)%ishb.eq.3)then

	    IF (BINIO.EQ.0)write(10,"(6i12)")countb,ifac(i)%ifa(1,1),ifac(i)%ifa(1,2),ifac(i)%ifa(1,3),ifac(i)%ifa(1,3),ibtr
	    IF (BINIO.EQ.1)write(10)countb,ifac(i)%ifa(1,1),ifac(i)%ifa(1,2),ifac(i)%ifa(1,3),ifac(i)%ifa(1,3),ibtr

	    end if

	    if (ifac(i)%ishb.eq.4)then

	    IF (BINIO.EQ.0)write(10,"(6i12)")countb,ifac(i)%ifa(1,1),ifac(i)%ifa(1,2),ifac(i)%ifa(1,3),ifac(i)%ifa(1,4),ibtr
	    IF (BINIO.EQ.1)write(10)countb,ifac(i)%ifa(1,1),ifac(i)%ifa(1,2),ifac(i)%ifa(1,3),ifac(i)%ifa(1,4),ibtr

	    end if






	end select

	end if

end do





close(10)

deallocate(ifac)
deallocate(ishape)
deallocate(iele)


! 
! 
! 
! 
! 
! 
! 
! ! read dimensions
! ! read(82,"(A1,I1,1x,i1)",advance='no')braco,dum,dum1
! ! if (dum .eq. 2) then
! ! dimen=dum1
! ! read(82,"(1x)",advance='YES')
! ! end if
! ! 
! ! (0 " Created by : Fluent_V6 Interface Vers. 14.0.3")
! ! (2 2)
! ! (0 "Node Section")
! ! (10 (0 1 81 0 2))
! ! (10 (2 1 81 1 2)
! ! (
! 
! 
! 
! read(82,"(A1,I2)",advance='NO')braco,dum
! ! if comment preceed index = 0
! if (dum.eq.0) then
! read(82,"(1x)",advance='YES')
! endif
! 
! read(82,"(A1,I2,1x,A1,I1,1x,I1,z4)",advance='no')braco,dum,braco,zoneid,in1,imaxn
! 
! 
! ! read(82,"(A1,I1,1x,i1)")braco,dum,dimen!dum,dimen
!  print*,braco,dum,dimen
! 
! read(82,*)
! 
!   read(82,"(A1,I2,1x,A1,I1,1x,I1,z4)")braco,index10,braco,zoneid,in1,imaxn
! !  read(82,*)braco,index10,zoneid,in1
!  read(82,*)
!  read(82,*)
! 
! allocate(x(imaxn))
! allocate(y(imaxn))
! if (dimen.eq.3) allocate(z(imaxn))
! 
! do i=1,imaxn
! read(82,*)x(i),y(i)
! end do
! 
! read(82,*)
! read(82,"(A1,I2,1x,A1,I1,1x,I1,z4)")braco,index10,braco,zoneid,in1,imaxe
! 
! 
! ! (12 (0 1 a6 0 0))
! ! (12 (3 1 a6 1 0)
! 
! 
! close(82)
! 
! 
! 
! 
! 
! 
! ! print*,braco,dum,dimen
! print*,index10,zoneid,in1,imaxn
! 
! 
! ! (0 " Created by : Fluent_V6 Interface Vers. 14.0.3")
! ! (2 2)
! ! (0 "Node Section")
! ! (10 (0 1 81 0 2))
! ! (10 (2 1 81 1 2)
! ! (
! ! 
! ! 
! ! (0 " Created by : Fluent_V6 Interface Vers. 14.0.3")
! ! (2 3)
! ! (0 "Node Section")
! ! (10 (0 1 13a 0 3))
! ! (10 (a 1 13a 1 3)
! ! (
! 
! 
! ! end subroutine init
! 
! 
! 
! ! subroutine test034 ( )
! ! 
! ! !*****************************************************************************80
! ! !
! ! !! TEST034 tests CHR4_TO_8 and CHR8_TO_4.
! ! !
! ! !  Licensing:
! ! !
! ! !    This code is distributed under the GNU LGPL license.
! ! !
! ! !  Modified:
! ! !
! ! !    19 January 2007
! ! !
! ! !  Author:
! ! !
! ! !    John Burkardt
! ! !
! !   implicit none
! ! 
! !   character chrtmp
! !   character chrtmp2
! !   integer ( kind = 4 ) i
! !   integer ( kind = 4 ) ichr
! !   integer ( kind = 4 ) j
! !   character ( len = 256 ) s1
! !   character ( len = 512 ) s2
! !   character ( len = 256 ) s3
! ! 
! !   write ( *, '(a)' ) ' '
! !   write ( *, '(a)' ) 'TEST034'
! !   write ( *, '(a)' ) '  CHR8_TO_4 convert characters to pairs of hexadecimals.'
! !   write ( *, '(a)' ) '  CHR4_TO_8 converts pairs of hexadecimals to characters.'
! !   write ( *, '(a)' ) ' '
! ! 
! !   do i = 1, 256
! !     s1(i:i) = char(i-1)
! !   end do
! ! 
! !   call chr8_to_4 ( s1, s2 )
! ! 
! !   call chr4_to_8 ( s2, s3 )
! ! 
! !   write ( *, '(a)' ) ' '
! !   write ( *, '(a)' ) '  Coded characters that can''t be printed are shown as blanks.'
! !   write ( *, '(a)' ) ' '
! !   write ( *, '(a)' ) '   ASCII  Coded  Decoded'
! !   write ( *, '(a)' ) ' '
! ! 
! !   do i = 1, 256
! ! 
! !     ichr = i - 1
! !     j = 2 * i - 1
! ! 
! !     if ( 33 <= ichr .and. ichr <= 127 ) then
! !       chrtmp = s1(i:i)
! !       chrtmp2 = s3(i:i)
! !     else
! !       chrtmp = ' '
! !       chrtmp2 = ' '
! !     end if
! ! 
! !     write ( *, '(2x,i3,1x,a1,6x,a2,7x,a1)' ) ichr, chrtmp, s2(j:j+1), chrtmp2
! ! 
! !   end do
! ! 
! !   return
! ! end


End subroutine drive




subroutine hex_to_i4 ( s, i4 )

!*****************************************************************************80
!
!! HEX_TO_I4 converts a hexadecimal string to an I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 December 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) S, the string of hexadecimal digits.
!
!    Output, integer ( kind = 4 ) I4, the corresponding I4.
!
  implicit none

  integer ( kind = 4 ) first
  integer ( kind = 4 ) idig
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) isgn
  integer ( kind = 4 ) j
  character ( len = * ) s
  integer ( kind = 4 ) s_length

  s_length = len_trim ( s )
!
!  Determine if there is a plus or minus sign.
!
  isgn = 1

  first = s_length + 1

  do j = 1, s_length

    if ( s(j:j) == '-' ) then
      isgn = -1
    else if ( s(j:j) == '+' ) then
      isgn = + 1
    else if ( s(j:j) /= ' ' ) then
      first = j
      exit
    end if

  end do
!
!  Read the numeric portion of the string.
!
  i4 = 0

  do j = first, s_length
    call hex_digit_to_i4 ( s(j:j), idig )
    i4 = i4 * 16 + idig
  end do

  i4 = isgn * i4

  return
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine hex_digit_to_i4 ( ch, i )

!*****************************************************************************80
!
!! HEX_DIGIT_TO_I4 converts a hexadecimal digit to an I4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 August 2009
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character CH, the hexadecimal digit, '0'
!    through '9', or 'A' through 'F', or also 'a' through 'f'
!    are allowed.
!
!    Output, integer ( kind = 4 ) I, the corresponding integer, or -1 if
!    CH was illegal.
!
  implicit none

  character ch
  integer ( kind = 4 ) i

  i = iachar ( ch )

  if ( lle ( '0', ch ) .and. lle ( ch, '9' ) ) then

    i = i - 48

  else if ( 65 <= i .and. i <= 70 ) then

    i = i - 55

  else if ( 97 <= i .and. i <= 102 ) then

    i = i - 87

  else if ( ch == ' ' ) then

    i = 0

  else

    i = -1

  end if

  return
end subroutine



 subroutine string2int (gcharr,interray,intsizee)
implicit none

integer :: counto,i,il,intsizee,outint
integer,allocatable,dimension(:):: spaces!,interray
integer,allocatable,dimension(:),intent(inout):: interray
character (len = *) :: gcharr
 character (len =256)::char1
	il=len_trim(gcharr)
		counto=0
		Do i=1,il
		  if ( gcharr(i:i) .eq. ' ') then 
 		  counto=counto+1
		  end if
		End do
		allocate(spaces(0:counto+1))
		spaces(:)=0
		counto=0
		spaces(1)=0
		Do i=1,il
		  if ( gcharr(i:i) .eq. ' ') then 
 		  counto=counto+1
		  spaces(counto)=i
		  end if
		End do
		    spaces(counto+1)=il+1
		intsizee=counto+1
		allocate (interray(intsizee))
		interray=0
		do i=0,counto
		    if (i .eq. 1) then
		    char1 = char1(1:spaces(i))
		    end if
		  char1=gcharr(spaces(i)+1:spaces(i+1)-1)
		    call hex_to_i4 ( char1, outint )
		  interray(i+1)=outint
		end do

deallocate(spaces)

end subroutine

subroutine removebrac(ch1,ch2)
implicit none
character(len=256) ::ch1,ch2,chdum
integer :: il,dum,dum1
		
! 		il=len_trim(gchar)
! 		dum= SCAN (gchar,'(' ,BACK = .true.)		! removes first braquet
! 		gchar1=gchar(dum+1:il)
! 		dum1= SCAN (gchar1,')' ,BACK = .false.)	! removes last braquet
! 		gchar2=gchar1(1:dum1-1)


		il=len_trim(ch1)
		dum= SCAN (ch1,'(' ,BACK = .false.)		! removes first braquet
		chdum=ch1(dum+1:il)
		dum1= SCAN (chdum,')' ,BACK = .false.)	! removes last braquet
		ch2=chdum(1:dum1-1)
! print*,ch1,ch2,'hello'
end subroutine removebrac









subroutine TRANSUGRID
IMPLICIT NONE
	INTEGER::I,J,K,L,N,I1,I2,I3,I4,I5,I6,I7,I8,IOS,IOX,IOY,IMAXEu,IMAXBu,IMAXNu,ICG,KX,NBOUND,DIP
	INTEGER::afnnodesg ! = number of nodes
        INTEGER::afntface  ! = number of boundary triangles
        INTEGER::afnqface  ! = number of boundary quads
        INTEGER::afntet    ! = number of volume TETRA_4 elements
        INTEGER::afnpyr    ! = number of volume PYRA_5 elements
        INTEGER::afnprz    ! = number of volume PENTA_6 elements
        INTEGER::afnhex    ! = number of volume HEXA_8 elements
        INTEGER,ALLOCATABLE,DIMENSION(:)::IBID,IBX,IBXX,ifacetag
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::if2nt,if2nq,ic2nt,ic2np,ic2nz,ic2nh
	real,allocatable,dimension(:)::x,y,z



OPEN(180,FILE="grid.ugrid",FORM='UNFORMATTED',STATUS='OLD',ACCESS='STREAM',CONVERT="BIG_ENDIAN")
read(180)afnnodesg,afntface,afnqface, afntet, afnpyr, afnprz, afnhex
print*,afnnodesg,afntface, afnqface, afntet, afnpyr, afnprz, afnhex

allocate(x(afnnodesg),y(afnnodesg),z(afnnodesg))
allocate(if2nt(3,afntface),if2nq(4,afnqface),ifacetag(afntface+afnqface),ic2nt(4,afntet),ic2np(5,afnpyr),ic2nz(6,afnprz),ic2nh(8,afnhex))

imaxeu=afntet+afnpyr+afnprz+afnhex
imaxbu=afntface+afnqface
imaxnu=afnnodesg

			 do i=1,afnnodesg;read(180)x(i),y(i),z(i);end do


			 do i=1,afntface;do j=1,3;read(180)if2nt(j,i); end do;end do
print*,"2"
			 do i=1,afnqface;do j=1,4;read(180)if2nq(j,i); end do;end do
print*,"3"
			 do i=1,afntface+afnqface;read(180)ifacetag(i);end do
print*,"4"
			 do i=1,afntet; do j=1,4; read(180)ic2nt(j,i);end do; end do
print*,"5"
			 do i=1,afnpyr; do j=1,5; read(180)ic2np(j,i);end do; end do
print*,"6"
			 do i=1,afnprz; do j=1,6; read(180)ic2nz(j,i);end do; end do
print*,"7"
			 do i=1,afnhex; do j=1,8; read(180)ic2nh(j,i);end do; end do
print*,"8"

close(180)
OPEN(120,FILE="grid.mapbc",STATUS='OLD',FORM='FORMATTED')
READ (120,*) NBOUND

ALLOCATE (IBID(NBOUND),IBX(NBOUND),IBXX(NBOUND))
DO I=1,NBOUND
  READ(120,*)IBID(I),IBX(I)
END DO




DO I=1,NBOUND

	select case(ibx(i))

	case (5000,5050)	!farfield
	IBXX(i)=6
	case(6662,6661,6663)	!symmetry
	IBXX(i)=3
	case(4000)	!wall
	IBXX(i)=4
	case(7031)	!outflow
	IBXX(i)=2

	case(7036,7100)	!inflow
	IBXX(i)=1

	case(6100)	!periodicity
	IBXX(i)=5

	end select

end do








			  !write nodes first
			  OPEN(12,FILE="GRID.vrt",FORM='unformatted',ACTION='WRITE')
			 do i=1,afnnodesg
			write(12)i,x(i),y(i),z(i)
! 			 write(1200,"(5X,I8,2X,ES21.14,2X,ES21.14,2X,ES21.14)")i,x(i),y(i),z(i)
			 end do
			 close(12)

			 !end nodes writing


			 !write elements now
			 kx=0
			 OPEN(11,FILE="GRID.cel",FORM='unformatted',ACTION='WRITE')

			 !tetra: 1 2 3 3 4 4 4 4
			 do i=1,afntet
			    kx=kx+1
! 			    write(150,"(9I10)")kx,ic2nt(1,i),ic2nt(2,i),ic2nt(3,i),ic2nt(3,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i)
				write(11)kx,ic2nt(1,i),ic2nt(2,i),ic2nt(3,i),ic2nt(3,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i)

			 end do
			 !pyramid: 1 2 3 4 5 5 5 5


			 do i=1,afnpyr
			    kx=kx+1
 			    !write(150,*)kx,ic2np(1,i),ic2np(2,i),ic2np(5,i),ic2np(4,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
! 			    	write(150,"(9I10)")kx,ic2np(1,i),ic2np(4,i),ic2np(5,i),ic2np(2,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
				write(11)kx,ic2np(1,i),ic2np(4,i),ic2np(5,i),ic2np(2,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
			 end do
			 !prism: 1 2 3 3 4 5 6 6





			 do i=1,afnprz
			    kx=kx+1
 			    !write(11,"(9I10)")kx,ic2nz(1,i),ic2nz(2,i),ic2nz(3,i),ic2nz(3,i),ic2nz(4,i),ic2nz(5,i),ic2nz(6,i),ic2nz(6,i)
				write(11)kx,ic2nz(1,i),ic2nz(2,i),ic2nz(3,i),ic2nz(3,i),ic2nz(4,i),ic2nz(5,i),ic2nz(6,i),ic2nz(6,i)
			 end do
			 !hexa: 1 2 3 4 5 6 7 8
			  do i=1,afnhex
			    kx=kx+1
! 			    write(150,"(9I10)")kx,ic2nh(1,i),ic2nh(2,i),ic2nh(3,i),ic2nh(4,i),ic2nh(5,i),ic2nh(6,i),ic2nh(7,i),ic2nh(8,i)
				write(11)kx,ic2nh(1,i),ic2nh(2,i),ic2nh(3,i),ic2nh(4,i),ic2nh(5,i),ic2nh(6,i),ic2nh(7,i),ic2nh(8,i)
			 end do
			 close(11)
			 !end writing elements

			 !now write the boundary file
			 OPEN(10,FILE="GRID.bnd",FORM='unformatted',ACTION='WRITE')
			 kx=0
			 !triangle: 1 2 3 3
			  do i=1,afntface
			    kx=kx+1
! 			    write(1000,"(6I12)")kx,if2nt(1,i),if2nt(2,i),if2nt(3,i),if2nt(3,i),ibxx(ifacetag(kx))
			   write(10)kx,if2nt(1,i),if2nt(2,i),if2nt(3,i),if2nt(3,i),ibxx(ifacetag(kx))
			 end do
			 !quad: 1 2 3 4
			  do i=1,afnqface
			    kx=kx+1
! 			    write(1000,"(6I12)")kx,if2nq(1,i),if2nq(2,i),if2nq(3,i),if2nq(4,i),ibxx(ifacetag(kx))
				write(10)kx,if2nq(1,i),if2nq(2,i),if2nq(3,i),if2nq(4,i),ibxx(ifacetag(kx))
			 end do

			close(10)

			 DEALLOCATE(IBID,IBX,IBXX, x,y,z, if2nt,if2nq,ifacetag,ic2nt,ic2np,ic2nz,ic2nh)


END SUBROUTINE













 END MODULE translate
