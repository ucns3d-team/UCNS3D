module translate
use declaration
implicit none
integer,allocatable,dimension(:):: interray


contains

subroutine translate_mesh
!> @brief
!> subroutine for transforming fluent style msh file to native format
implicit none
logical::heres,heref,hereu
 character(len=20)::proc,ucns3dfile,fluentfile,ugridfile
 	ucns3dfile='GRID.bnd'
 	fluentfile='grid.msh'
 	ugridfile='grid.ugrid'
	
	inquire (file=ucns3dfile,exist=heres)
	if (heres) then
	
	!proceed without translation


	
	else


			inquire (file=fluentfile,exist=heref)


			if (heref)then

			call drive(interray)

			end if



			inquire (file=ugridfile,exist=hereu)


			if (hereu)then

			call transugrid

			end if


	end if



end subroutine translate_mesh


!!! http://people.sc.fsu.edu/~jburkardt/f_src/chrpak/chrpak.html
subroutine drive(interray)
implicit none

type::anode_number1	!name of type for the set of nodes
	integer::noden	!identification number that can be used as a pointer inside an array
	real::x		!coordinates in x axis
	real::y		!coordinates in y axis
	real::z
end type anode_number1
type::aelement_number1	!name of type for the set of elements
    integer::ieindex
    integer::iecounter
    integer:: iface
    integer,allocatable,dimension(:,:)::faces ! id of face and id of node
    integer,allocatable,dimension(:):: nd
    integer ::ishape ! id of shape ! 1 triangle, 2 tetra, 3 quad, 4 hexa, 5 pyramid, 6 prism
end type aelement_number1
type::aboundary_number1	!name of type for the set of elements
    integer::ibindex
    integer::ibcounter
    integer:: ibtype
    integer,allocatable,dimension(:):: ndb
    integer ::ibshape ! 1 line ! 3 triangle ! 4 quad
end type aboundary_number1
type::aface_number	!name of type for the set of elements
    integer::ifindex
    integer::ifcounter
    integer::ifacbtype,ishb
!     integer::ifmaxnode
    integer,allocatable,dimension(:,:):: ifa ! elements id (always two) , nodes id only 1s row
    integer ::ifshape ! 1 line ! 3 triangle ! 4 quad
end type aface_number
type(aelement_number1),allocatable,dimension(:)::diele
type(anode_number1),allocatable,dimension(:)::dinod
type(aboundary_number1),allocatable,dimension(:)::dibou
type(aface_number),allocatable,dimension(:)::difac

integer::ing2,jj,dimen,imaxe,imaxn,imaxb,dum,index10,zoneid,in1,dum1,bctypdum,iosx,ii,vrt1,vrt2,vrt3,vrt4,vct1,vct2,vct3,vct4,countb,eltype,ing,ibtr,nfin,icte
integer::lexist,checkbrac,dum2,dum3,dum4,ierr,il,il1,il2,ispace1,ispace2
integer :: countline,counto,countword,lengt,intsize,imaxnglobal,imaxeglobal,imin,imax,imaxfglobal,index1,ina,in2,ichen,iix,iiy,icountfc,corn,icorn
integer,dimension(4)::nod,nodx,cans,xcand,cane,cang,canh,iiz,canf
integer,allocatable,dimension(:),intent(inout):: interray
integer,allocatable,dimension(:)::spaces,integerarray,ishape
integer,allocatable,dimension(:,:)::ifaci,iface
character(len=3) :: checkcomm,checkdim,checknod,checkelety,checkfac,checkbcne
integer::ios,dumhex,ioss,endline,str2int
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


  

open(82,file="grid.msh",form='formatted',status='old',action='read',iostat=ios)
 countline=0
do 
	read(82,"(a3)",advance='no',iostat=ios)dumc
	dumc2=dumc(2:3)
	read(dumc2,*,iostat=ioss)str2int
	if ((ios .eq. -1))goto 11
	countline=countline+1
! 	print*,'line',countlieq. "(0")ne !,dumc
	     ! if (dumc .eq. "(0")  then ! comment conditions
	      if (dumc .eq. "(0 ")  then ! comment conditions
		  read(82,*,iostat=ios)
		  countline=countline+1
		  go to 10
	      endif
! 	      if (dumc .eq. "(2") then
! 	      read(82,*,iostat=ios)dimen
	       if (dumc .eq. "(2 ") then
       read(82,'(i1)',iostat=ios) dimen !changed by holger foysi
		  countline=countline+1
! 	      print*,'dimension:',dimen
	      end if
	      if (dumc .eq. "(10") then
		      read(82,'(a)',iostat=ios) gchar
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
! 			  print*,'max number nodes:',imaxnglobal
			  if (dimen.eq.3) allocate(z(imaxnglobal))
		      endif
		      if (interray(1) .ne. 0) then 
			read(82,"(a1)",advance='no',iostat=ios)dumc
			countline=countline+1
			  if (dumc .eq. "(") then
			      do i=interray(2),interray(3)
				if (dimen.eq.2) then
				read(82,*,iostat=ios)x(i),y(i)
				countline=countline+1
				endif
				if (dimen.eq.3) then
				read(82,*,iostat=ios)x(i),y(i),z(i)
				countline=countline+1
				endif
			      enddo
			  else
			    backspace(82,iostat=ios)
			  do i=interray(2),interray(3)
				if (dimen.eq.2) then
				read(82,*,iostat=ios)x(i),y(i)
				countline=countline+1
				endif
				if (dimen.eq.3) then
				read (82,*,iostat=ios)x(i),y(i),z(i)
				countline=countline+1
				endif
			      enddo
			  endif
			  if (binio.eq.0)open(10,file="GRID.vrt",form='formatted',action='write',iostat=iosx,position='append')
			  if (binio.eq.1)open(10,file="GRID.vrt",form='unformatted',action='write',iostat=iosx,position='append')
			  selectcase (dimen)
			    case(2)
			    do i=interray(2),interray(3)
			      if (binio.eq.0)write(10,"(5x, i8, 2x,es21.14,2x,es21.14)")i,x(i),y(i)
			      if (binio.eq.1)write(10)i,x(i),y(i)
			    end do
			    close(10)
			    case(3)
			    do i=interray(2),interray(3)
			      if (binio.eq.0)write(10,"(5x,i8,2x,es21.14,2x,es21.14,2x,es21.14)")i,x(i),y(i),z(i)
			      if (binio.eq.1)write(10)i,x(i),y(i),z(i)
			    end do 
			    close(10)
			    deallocate(x,y)
                            if (dimen.eq.3) deallocate(z)
!                             
			  end select
		      endif
		deallocate(interray)
		
! 			  print*,'max number nodes:',imaxnglobal
! 			  if (dimen.eq.3) allocate(z(imaxnglobal))
	      end if ! (10
	      if (dumc .eq. "(12") then
			  read(82,'(a)',iostat=ios) gchar
			  countline=countline+1
			  call removebrac(gchar,gchar2)
	  ! 		print*,gchar2,"edw1"!,il
			  call string2int (gchar2,interray,intsize)
			  if (interray(1) .eq. 0)then
			      imaxeglobal=interray(3) ! set global maximum number of elements
			      allocate(diele(imaxeglobal))
			      allocate(ishape(imaxeglobal))
! 			      print*,'max number cells:',imaxeglobal
			      diele(1:imaxeglobal)%iecounter=0
			  endif

!                         
			   if ((interray(1) .ne. 0))then
			    if (interray(5).ne.0) then 
			      diele(interray(2):interray(3))%ishape=interray(5)
			       do i=interray(2),interray(3)
! 				      
				      diele(i)%ieindex=i
				       select case (diele(i)%ishape)
					  case(1)
					    diele(i)%iface=3
					    allocate(diele(i)%faces(1:3,1:2))
					    !allocate(diele(i)%nd(1:44))
					  case(3)
					    diele(i)%iface=4
					    allocate(diele(i)%faces(1:4,1:2))
					    !allocate(diele(i)%nd(4))
					  case(2)
					    diele(i)%iface=4
					    allocate(diele(i)%faces(1:4,1:3))
! 					    allocate(diele(i)%nd(8))
					  case(4)
					    diele(i)%iface=6
					    allocate(diele(i)%faces(1:6,1:4))
! 					    allocate(diele(i)%nd(8))
					  case(5)
					    diele(i)%iface=5
					    allocate(diele(i)%faces(1:5,1:4))
! 					    allocate(diele(i)%nd(8))
					  case(6)
					    diele(i)%iface=5
					    allocate(diele(i)%faces(1:5,1:4))
! 					    allocate(diele(i)%nd(8))

				    end select

				      diele(i)%faces(:,:)=0
				  end do
			  end if
			  end if
			  if ((interray(1) .ne. 0))then
			    if (interray(5).eq.0) then 

			    read(82,"(a1)",advance='no',iostat=ios)dumc
			    countline=countline+1
			      if (dumc .eq. "(") then
				  read(82,*,iostat=ios)ishape(interray(2):interray(3))!diele(interray(2):interray(3))%ishape
				    diele(interray(2):interray(3))%ishape=ishape(interray(2):interray(3))
				  countline=countline+1
			      else
				backspace (82,iostat=ios)
			      read(82,*,iostat=ios)ishape(interray(2):interray(3))!diele(interray(2):interray(3))%ishape
				    diele(interray(2):interray(3))%ishape=ishape(interray(2):interray(3))
			      countline=countline+1
			      endif
				  do i=interray(2),interray(3)
! 				      
				      diele(i)%ieindex=i
				       select case (diele(i)%ishape)
					  case(1)
					    diele(i)%iface=3
					    allocate(diele(i)%faces(3,2))
					    
					  case(3)
					    diele(i)%iface=4
					    allocate(diele(i)%faces(4,2))
					   
					  case(2)
					    diele(i)%iface=4
					    allocate(diele(i)%faces(4,3))
					   
					  case(4)
					    diele(i)%iface=6
					    allocate(diele(i)%faces(6,4))
					   
					  case(5)
					    diele(i)%iface=5
					    allocate(diele(i)%faces(5,4))
					    
					  case(6)
					    diele(i)%iface=5
					    allocate(diele(i)%faces(5,4))
					    

				    end select

				      diele(i)%faces(:,:)=0
				  end do
			  endif
			  end if
! 			  
			  deallocate(interray)
	      end if
	      if (dumc .eq. "(13") then
		read(82,'(a)',advance='no',iostat=ios) gchar
		countline=countline+1
		call removebrac(gchar,gchar2)
		 call string2int (gchar2,interray,intsize)
		 index1=interray(1);imin=interray(2);imax=interray(3);bctypdum=interray(4)
		  if (index1.ne.0)eltype=interray(5)

		 deallocate(interray)
		if (index1 .eq. 0)then
		    imaxfglobal=imax ! set global maximum number of faces
		    allocate(difac(imaxfglobal))
		endif
		if (index1 .ne. 0) then 
		    
		  read(82,"(a1)",advance='no',iostat=ios)dumc
		  countline=countline+1
! 		  print*,'faces',imin,imax,index1,dumc
  		    if (dumc .eq. "(") then
!   		    
			   do i=imin,imax
! 				if (dimen.eq.2) then
				  read(82,'(a)',advance='no',iostat=ios) gchar
				  countline=countline+1
! 				  print*,gchar,'ssss'
				  call string2int (gchar,interray,intsize)
! 				  
				  select case (dimen)
				    case (2)! 2d
					     if (eltype.eq.0)then
					difac(i)%ishb=2
					difac(i)%ifacbtype=bctypdum
					allocate(difac(i)%ifa(2,2))
					difac(i)%ifa(1,1)=interray(2)
					difac(i)%ifa(1,2)=interray(3)
					difac(i)%ifa(2,1)=interray(4)
					difac(i)%ifa(2,2)=interray(5)
! 					if ( difac(i)%ifa(2,1) .ne. 0) then ! 1st element
! 					diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					end if
! 					if ( difac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! 					diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,2)
! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,1)
! 					endif
					else
					 difac(i)%ishb=2
					    difac(i)%ifacbtype=bctypdum
					allocate(difac(i)%ifa(2,2))
					difac(i)%ifa(1,1)=interray(1)
					difac(i)%ifa(1,2)=interray(2)
					difac(i)%ifa(2,1)=interray(3)
					difac(i)%ifa(2,2)=interray(4)
! 					if ( difac(i)%ifa(2,1) .ne. 0) then ! 1st element
! 					diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					end if
! 					if ( difac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! 					diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,2)
! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,1)
! 					endif


					    end if

					
				    case (3)!3d
					difac(i)%ifacbtype=bctypdum

					  if (eltype.eq.0)then
					difac(i)%ishb=interray(1)
					select case (interray(1))
					  case (3)! triangles
					      allocate(difac(i)%ifa(2,3))
					      difac(i)%ifa(1,1)=interray(2)
					      difac(i)%ifa(1,2)=interray(3)
					      difac(i)%ifa(1,3)=interray(4)
					      difac(i)%ifa(2,1)=interray(5)
					      difac(i)%ifa(2,2)=interray(6)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then ! no wall
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      allocate(difac(i)%ifa(2,4))
					      difac(i)%ifa(1,1)=interray(2)
					      difac(i)%ifa(1,2)=interray(3)
					      difac(i)%ifa(1,3)=interray(4)
					      difac(i)%ifa(1,4)=interray(5)
					      difac(i)%ifa(2,1)=interray(6)
					      difac(i)%ifa(2,2)=interray(7)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,4)=difac(i)%ifa(1,4)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,4)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,4)=difac(i)%ifa(1,1)
! 					      endif

					  end select
					    else
					      difac(i)%ishb=eltype
					select case (eltype)
					  case (3)! triangles
					      allocate(difac(i)%ifa(2,3))
					      difac(i)%ifa(1,1)=interray(1)
					      difac(i)%ifa(1,2)=interray(2)
					      difac(i)%ifa(1,3)=interray(3)
					      difac(i)%ifa(2,1)=interray(4)
					      difac(i)%ifa(2,2)=interray(5)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then ! no wall
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      allocate(difac(i)%ifa(2,4))
					      difac(i)%ifa(1,1)=interray(1)
					      difac(i)%ifa(1,2)=interray(2)
					      difac(i)%ifa(1,3)=interray(3)
					      difac(i)%ifa(1,4)=interray(4)
					      difac(i)%ifa(2,1)=interray(5)
					      difac(i)%ifa(2,2)=interray(6)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,4)=difac(i)%ifa(1,4)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,4)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,4)=difac(i)%ifa(1,1)
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
     			   do i=imin,imax
! 				if (dimen.eq.2) then
				  read(82,'(a)',advance='no',iostat=ios) gchar
				  countline=countline+1
! 				  print*,gchar,'ssss'
				  call string2int (gchar,interray,intsize)
! 				  
				  select case (dimen)
				    case (2)! 2d
					difac(i)%ifacbtype=bctypdum
					     if (eltype.eq.0)then
 					difac(i)%ishb=2
					difac(i)%ifacbtype=bctypdum
					allocate(difac(i)%ifa(2,2))
					difac(i)%ifa(1,1)=interray(2)
					difac(i)%ifa(1,2)=interray(3)
					difac(i)%ifa(2,1)=interray(4)
					difac(i)%ifa(2,2)=interray(5)
! 					if ( difac(i)%ifa(2,1) .ne. 0) then ! 1st element
! ! 					diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! ! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! ! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! ! 					end if
! ! 					if ( difac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! ! 					diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! ! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,2)
! ! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,1)
! ! 					endif
					else
					 difac(i)%ishb=2
					    difac(i)%ifacbtype=bctypdum
					allocate(difac(i)%ifa(2,2))
					difac(i)%ifa(1,1)=interray(1)
					difac(i)%ifa(1,2)=interray(2)
					difac(i)%ifa(2,1)=interray(3)
					difac(i)%ifa(2,2)=interray(4)
! 					if ( difac(i)%ifa(2,1) .ne. 0) then ! 1st element
! 					diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					end if
! 					if ( difac(i)%ifa(2,2) .ne. 0) then ! 2nd element
! 					diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,2)
! 					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,1)
! 					endif


					    end if
				    case (3)!3d
					difac(i)%ifacbtype=bctypdum

					  if (eltype.eq.0)then
					difac(i)%ishb=interray(1)
					select case (interray(1))
					  case (3)! triangles
					      allocate(difac(i)%ifa(2,3))
					      difac(i)%ifa(1,1)=interray(2)
					      difac(i)%ifa(1,2)=interray(3)
					      difac(i)%ifa(1,3)=interray(4)
					      difac(i)%ifa(2,1)=interray(5)
					      difac(i)%ifa(2,2)=interray(6)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then ! no wall
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      allocate(difac(i)%ifa(2,4))
					      difac(i)%ifa(1,1)=interray(2)
					      difac(i)%ifa(1,2)=interray(3)
					      difac(i)%ifa(1,3)=interray(4)
					      difac(i)%ifa(1,4)=interray(5)
					      difac(i)%ifa(2,1)=interray(6)
					      difac(i)%ifa(2,2)=interray(7)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,4)=difac(i)%ifa(1,4)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,4)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,4)=difac(i)%ifa(1,1)
! 					      endif
					    end select
					    else
					      difac(i)%ishb=eltype
					select case (eltype)
					  case (3)! triangles
						difac(i)%ishb=eltype
					      allocate(difac(i)%ifa(2,3))
					      difac(i)%ifa(1,1)=interray(1)
					      difac(i)%ifa(1,2)=interray(2)
					      difac(i)%ifa(1,3)=interray(3)
					      difac(i)%ifa(2,1)=interray(4)
					      difac(i)%ifa(2,2)=interray(5)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then ! no wall
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,1)
! 					      endif
					  case(4)! quads
					      difac(i)%ishb=eltype
					      allocate(difac(i)%ifa(2,4))
					      difac(i)%ifa(1,1)=interray(1)
					      difac(i)%ifa(1,2)=interray(2)
					      difac(i)%ifa(1,3)=interray(3)
					      difac(i)%ifa(1,4)=interray(4)
					      difac(i)%ifa(2,1)=interray(5)
					      difac(i)%ifa(2,2)=interray(6)
! 					      if ((difac(i)%ifa(2,1)).ne.0)then
! 					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,4)=difac(i)%ifa(1,4)
! 					      endif
! 					      if ((difac(i)%ifa(2,2)).ne.0)then
! 					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,4)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,3)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,2)
! 					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,4)=difac(i)%ifa(1,1)
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
					if  (difac(i)%ishb.eq.2)then

					if ( difac(i)%ifa(2,1) .ne. 0) then ! 1st element
					diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
					diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
					end if
					if ( difac(i)%ifa(2,2) .ne. 0) then ! 2nd element
					diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,2)
					diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,1)
					endif
					end if
					if  (difac(i)%ishb.eq.4)then
					 if ((difac(i)%ifa(2,1)).ne.0)then
					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,4)=difac(i)%ifa(1,4)
					      endif
					      if ((difac(i)%ifa(2,2)).ne.0)then
					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,4)
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,3)
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,2)
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,4)=difac(i)%ifa(1,1)
					      endif

					end if
					if  (difac(i)%ishb.eq.3)then
					 if ((difac(i)%ifa(2,1)).ne.0)then ! no wall
					      diele(difac(i)%ifa(2,1))%iecounter=diele(difac(i)%ifa(2,1))%iecounter+1
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,1)=difac(i)%ifa(1,1)
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,2)=difac(i)%ifa(1,2)
					      diele(difac(i)%ifa(2,1))%faces(diele(difac(i)%ifa(2,1))%iecounter,3)=difac(i)%ifa(1,3)
					      endif
					      if ((difac(i)%ifa(2,2)).ne.0)then
					      diele(difac(i)%ifa(2,2))%iecounter=diele(difac(i)%ifa(2,2))%iecounter+1
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,1)=difac(i)%ifa(1,3)
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,2)=difac(i)%ifa(1,2)
					      diele(difac(i)%ifa(2,2))%faces(diele(difac(i)%ifa(2,2))%iecounter,3)=difac(i)%ifa(1,1)
					      endif

				      end if
end do



if (binio.eq.0)then
open(10,file="GRID.cel",form='formatted',action='write',iostat=iosx)
else
open(10,file="GRID.cel",form='unformatted',action='write',iostat=iosx)

end if
! 
!  integer,allocatable,dimension(:,:)::faces ! id of face and id of node
!     integer,allocatable,dimension(:):: nd
icte=0

do i=1,imaxeglobal


      select case(diele(i)%ishape)


      case(1)

	vrt1=diele(i)%faces(1,1);vrt2=diele(i)%faces(1,2)
	do ii=2,diele(i)%iecounter
	    vct1=diele(i)%faces(ii,1)
	  if (((vrt1.eq.vct1)).or.((vrt2.eq.vct1)))then


	  else
! 	  write(10,"(5i10)")i,vrt2,vrt1,vct1,vct1
	  if (binio.eq.0)write(10,"(5i10)")i,vct1,vrt1,vrt2,vrt2
	  if (binio.eq.1)write(10)i,vct1,vrt1,vrt2,vrt2
	  cycle
	  end if
	end do



      case(3)
	vrt1=diele(i)%faces(1,1);vrt2=diele(i)%faces(1,2)
	do ii=2,diele(i)%iecounter
	    vct1=diele(i)%faces(ii,1);vct2=diele(i)%faces(ii,2)
	  if (((vrt1.eq.vct1).or.(vrt1.eq.vct2)).or.((vrt2.eq.vct1).or.(vrt2.eq.vct2)))then


	  else
	  if (binio.eq.0)write(10,"(5i10)")i,vrt1,vrt2,vct1,vct2
	  if (binio.eq.1)write(10)i,vrt1,vrt2,vct1,vct2
	 
	  cycle
	  end if
	end do


      case(2)
	vrt1=diele(i)%faces(1,1);vrt2=diele(i)%faces(1,2);vrt3=diele(i)%faces(1,3)

	nod(1)=vrt1
	nod(2)=vrt2
	nod(3)=vrt3

	do ii=2,diele(i)%iecounter
	       vct1=diele(i)%faces(ii,1);vct2=diele(i)%faces(ii,2);vct3=diele(i)%faces(ii,3)
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
	   
   	       if (binio.eq.0)write(10,"(9i10)")i,vrt1,vrt2,vrt3,vrt3,nod(4),nod(4),nod(4),nod(4)		!or
		if (binio.eq.1)write(10)i,vrt1,vrt2,vrt3,vrt3,nod(4),nod(4),nod(4),nod(4)
		


	case(4)
		ing=0
		  vrt1=diele(i)%faces(1,1);vrt2=diele(i)%faces(1,2);vrt3=diele(i)%faces(1,3);vrt4=diele(i)%faces(1,4)
	do ii=2,diele(i)%iecounter
! 	    print*,"edw1",vrt1,vrt2,vrt3,ii,diele(i)%iecounter
	    vct1=diele(i)%faces(ii,1);vct2=diele(i)%faces(ii,2);vct3=diele(i)%faces(ii,3);vct4=diele(i)%faces(ii,4)
	  if( ((vrt1.eq.vct1).or.(vrt2.eq.vct1).or.(vrt3.eq.vct1).or.(vrt4.eq.vct1)).or.&
	      ((vrt1.eq.vct2).or.(vrt2.eq.vct2).or.(vrt3.eq.vct2).or.(vrt4.eq.vct2)).or.&
	      ((vrt1.eq.vct3).or.(vrt2.eq.vct3).or.(vrt3.eq.vct3).or.(vrt4.eq.vct3)).or.&
	      ((vrt1.eq.vct4).or.(vrt2.eq.vct4).or.(vrt3.eq.vct4).or.(vrt4.eq.vct4)))then
	    
		

	  else





	 
!   		  write(10,"(9i10)")i,vrt1,vrt2,vrt3,vrt4,vct4,vct3,vct2,vct1

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

		    do iix=2,diele(i)%iecounter
			if (iix.ne.ing)then
			  canf(1)=diele(i)%faces(iix,1)
			  canf(2)=diele(i)%faces(iix,2)
			  canf(3)=diele(i)%faces(iix,3)
			  canf(4)=diele(i)%faces(iix,4)
			do iiy=1,4
			      if (canf(iiy).eq.xcand(1))then

				icountfc=icountfc+1
				    iiz(icountfc)=iix
			       end if
			end do
			
			end if
		    end do

		    icountfc=0

		    do iix=2,diele(i)%iecounter
			 if ((iix.eq.iiz(1)).or.(iix.eq.iiz(2)))then
			      icountfc=icountfc+1
			      if (icountfc.eq.1)then
			       cang(1)=diele(i)%faces(iix,1)
			  cang(2)=diele(i)%faces(iix,2)
			  cang(3)=diele(i)%faces(iix,3)
			  cang(4)=diele(i)%faces(iix,4)
			      end if
			      if (icountfc.eq.2)then
			       canh(1)=diele(i)%faces(iix,1)
			  canh(2)=diele(i)%faces(iix,2)
			  canh(3)=diele(i)%faces(iix,3)
			  canh(4)=diele(i)%faces(iix,4)
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

		   
 		  if (binio.eq.0)write(10,"(9i10)")i,xcand(1),xcand(2),xcand(3),xcand(4),cans(1),cans(4),cans(3),cans(2)
		  if (binio.eq.1)write(10)i,xcand(1),xcand(2),xcand(3),xcand(4),cans(1),cans(4),cans(3),cans(2)
		cycle

	    end if
	  end do

	case(5)

	    	 
	ing=0
	do ii=1,diele(i)%iecounter

	   

	    if (diele(i)%faces(ii,4).ne.0)then
	    
	      ing=ii
	      vct1=diele(i)%faces(ii,1);vct2=diele(i)%faces(ii,2);vct3=diele(i)%faces(ii,3);vct4=diele(i)%faces(ii,4)
	      nod(1)=vct1; nod(2)=vct2; nod(3)=vct3; nod(4)=vct4

	    end if

	end do
	do ii=1,diele(i)%iecounter

	     


		vrt1=diele(i)%faces(ii,1);vrt2=diele(i)%faces(ii,2);vrt3=diele(i)%faces(ii,3)
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
 		  if (binio.eq.0)write(10,"(9i10)")i,vct1,vct2,vct3,vct4,nfin,nfin,nfin,nfin
 		  if (binio.eq.1)write(10)i,vct1,vct2,vct3,vct4,nfin,nfin,nfin,nfin
! 		  

	

	case(6)
	ing=0
	ing2=0
		      
	do ii=1,diele(i)%iecounter

! 	    if vrt1=diele(i)%faces(1,1);vrt2=diele(i)%faces(1,2);vrt3=diele(i)%faces(1,3)

	     if (diele(i)%faces(ii,4).eq.0)then
	      
	      ing=ii
	      vrt1=diele(i)%faces(ii,1);vrt2=diele(i)%faces(ii,2);vrt3=diele(i)%faces(ii,3)
 	      cycle

	    end if
	end do
	do ii=1,diele(i)%iecounter
	    if (ii.ne.ing)then

	     if (diele(i)%faces(ii,4).eq.0)then
	     vct1=diele(i)%faces(ii,1);vct2=diele(i)%faces(ii,2);vct3=diele(i)%faces(ii,3)
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

		    do iix=1,diele(i)%iecounter
			if ((iix.ne.ing).and.(iix.ne.ing2))then
			  canf(1)=diele(i)%faces(iix,1)
			  canf(2)=diele(i)%faces(iix,2)
			  canf(3)=diele(i)%faces(iix,3)
			  canf(4)=diele(i)%faces(iix,4)
			do iiy=1,4
			      if (canf(iiy).eq.xcand(1))then

				icountfc=icountfc+1
				    iiz(icountfc)=iix
			       end if
			end do
			
			end if
		    end do

		    icountfc=0

		    do iix=1,diele(i)%iecounter
			 if ((iix.eq.iiz(1)).or.(iix.eq.iiz(2)))then
			      icountfc=icountfc+1
			      if (icountfc.eq.1)then
			       cang(1)=diele(i)%faces(iix,1)
			  cang(2)=diele(i)%faces(iix,2)
			  cang(3)=diele(i)%faces(iix,3)
			  cang(4)=diele(i)%faces(iix,4)
			      end if
			      if (icountfc.eq.2)then
			       canh(1)=diele(i)%faces(iix,1)
			  canh(2)=diele(i)%faces(iix,2)
			  canh(3)=diele(i)%faces(iix,3)
			  canh(4)=diele(i)%faces(iix,4)
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









   if (binio.eq.0)write(10,"(9i10)")i,xcand(1),xcand(2),xcand(3),xcand(3),cans(1),cans(3),cans(2),cans(2)
   if (binio.eq.1)write(10)i,xcand(1),xcand(2),xcand(3),xcand(3),cans(1),cans(3),cans(2),cans(2)

		!write(10,"(9i10)")i,vrt1,vrt2,vrt3,vrt3,vct1,vct2,vct3,vct3
!  		write(10,"(9i10)")icte,vrt1,vrt2,vrt3,vrt3,vct2,vct1,vct3,vct3
		
	      
! 	      



      endselect



	
end do

close(10) 


if (binio.eq.0)then
open(10,file="GRID.bnd",form='formatted',action='write',iostat=iosx)
else
open(10,file="GRID.bnd",form='unformatted',action='write',iostat=iosx)
end if
 countb=0
do i=1,imaxfglobal
	if (difac(i)%ifacbtype.ne.2)then



	  select case(difac(i)%ifacbtype)


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
	  if (binio.eq.0)write(10,"(6i12)")countb,difac(i)%ifa(1,1),difac(i)%ifa(1,2),0,0,ibtr
	  if (binio.eq.1)write(10)countb,difac(i)%ifa(1,1),difac(i)%ifa(1,2),0,0,ibtr
	  case(3)

	    if (difac(i)%ishb.eq.3)then

	    if (binio.eq.0)write(10,"(6i12)")countb,difac(i)%ifa(1,1),difac(i)%ifa(1,2),difac(i)%ifa(1,3),difac(i)%ifa(1,3),ibtr
	    if (binio.eq.1)write(10)countb,difac(i)%ifa(1,1),difac(i)%ifa(1,2),difac(i)%ifa(1,3),difac(i)%ifa(1,3),ibtr

	    end if

	    if (difac(i)%ishb.eq.4)then

	    if (binio.eq.0)write(10,"(6i12)")countb,difac(i)%ifa(1,1),difac(i)%ifa(1,2),difac(i)%ifa(1,3),difac(i)%ifa(1,4),ibtr
	    if (binio.eq.1)write(10)countb,difac(i)%ifa(1,1),difac(i)%ifa(1,2),difac(i)%ifa(1,3),difac(i)%ifa(1,4),ibtr

	    end if






	end select

	end if

end do





close(10)

deallocate(difac)
deallocate(ishape)
deallocate(diele)


! 
! 
! 
! 
! 
! 
! 
! ! read dimensions
! ! read(82,"(a1,i1,1x,i1)",advance='no')braco,dum,dum1
! ! if (dum .eq. 2) then
! ! dimen=dum1
! ! read(82,"(1x)",advance='yes')
! ! end if
! ! 
! ! (0 " created by : fluent_v6 interface vers. 14.0.3")
! ! (2 2)
! ! (0 "node section")
! ! (10 (0 1 81 0 2))
! ! (10 (2 1 81 1 2)
! ! (
! 
! 
! 
! read(82,"(a1,i2)",advance='no')braco,dum
! ! if comment preceed index = 0
! if (dum.eq.0) then
! read(82,"(1x)",advance='yes')
! endif
! 
! read(82,"(a1,i2,1x,a1,i1,1x,i1,z4)",advance='no')braco,dum,braco,zoneid,in1,imaxn
! 
! 
! ! read(82,"(a1,i1,1x,i1)")braco,dum,dimen!dum,dimen
!  print*,braco,dum,dimen
! 
! read(82,*)
! 
!   read(82,"(a1,i2,1x,a1,i1,1x,i1,z4)")braco,index10,braco,zoneid,in1,imaxn
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
! read(82,"(a1,i2,1x,a1,i1,1x,i1,z4)")braco,index10,braco,zoneid,in1,imaxe
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
! ! (0 " created by : fluent_v6 interface vers. 14.0.3")
! ! (2 2)
! ! (0 "node section")
! ! (10 (0 1 81 0 2))
! ! (10 (2 1 81 1 2)
! ! (
! ! 
! ! 
! ! (0 " created by : fluent_v6 interface vers. 14.0.3")
! ! (2 3)
! ! (0 "node section")
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
! ! !! test034 tests chr4_to_8 and chr8_to_4.
! ! !
! ! !  licensing:
! ! !
! ! !    this code is distributed under the gnu lgpl license.
! ! !
! ! !  modified:
! ! !
! ! !    19 january 2007
! ! !
! ! !  author:
! ! !
! ! !    john burkardt
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
! !   write ( *, '(a)' ) 'test034'
! !   write ( *, '(a)' ) '  chr8_to_4 convert characters to pairs of hexadecimals.'
! !   write ( *, '(a)' ) '  chr4_to_8 converts pairs of hexadecimals to characters.'
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
! !   write ( *, '(a)' ) '  coded characters that can''t be printed are shown as blanks.'
! !   write ( *, '(a)' ) ' '
! !   write ( *, '(a)' ) '   ascii  coded  decoded'
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


end subroutine drive




subroutine hex_to_i4 ( s, i4 )

!*****************************************************************************80
!
!! hex_to_i4 converts a hexadecimal string to an i4.
!
!  licensing:
!
!    this code is distributed under the gnu lgpl license.
!
!  modified:
!
!    07 december 2000
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, character ( len = * ) s, the string of hexadecimal digits.
!
!    output, integer ( kind = 4 ) i4, the corresponding i4.
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
!  determine if there is a plus or minus sign.
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
!  read the numeric portion of the string.
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
!! hex_digit_to_i4 converts a hexadecimal digit to an i4.
!
!  licensing:
!
!    this code is distributed under the gnu lgpl license.
!
!  modified:
!
!    31 august 2009
!
!  author:
!
!    john burkardt
!
!  parameters:
!
!    input, character ch, the hexadecimal digit, '0'
!    through '9', or 'a' through 'f', or also 'a' through 'f'
!    are allowed.
!
!    output, integer ( kind = 4 ) i, the corresponding integer, or -1 if
!    ch was illegal.
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
		do i=1,il
		  if ( gcharr(i:i) .eq. ' ') then 
 		  counto=counto+1
		  end if
		end do
		allocate(spaces(0:counto+1))
		spaces(:)=0
		counto=0
		spaces(1)=0
		do i=1,il
		  if ( gcharr(i:i) .eq. ' ') then 
 		  counto=counto+1
		  spaces(counto)=i
		  end if
		end do
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
! 		dum= scan (gchar,'(' ,back = .true.)		! removes first braquet
! 		gchar1=gchar(dum+1:il)
! 		dum1= scan (gchar1,')' ,back = .false.)	! removes last braquet
! 		gchar2=gchar1(1:dum1-1)


		il=len_trim(ch1)
		dum= scan (ch1,'(' ,back = .false.)		! removes first braquet
		chdum=ch1(dum+1:il)
		dum1= scan (chdum,')' ,back = .false.)	! removes last braquet
		ch2=chdum(1:dum1-1)
! print*,ch1,ch2,'hello'
end subroutine removebrac









subroutine transugrid
implicit none
	integer::i,j,k,l,ntg,i1,i2,i3,i4,i5,i6,i7,i8,ios,iox,ioy,imaxeu,imaxbu,imaxnu,icg,kx,nbound,dip
	integer::afnnodesg ! = number of nodes
        integer::afntface  ! = number of boundary triangles
        integer::afnqface  ! = number of boundary quads
        integer::afntet    ! = number of volume tetra_4 elements
        integer::afnpyr    ! = number of volume pyra_5 elements
        integer::afnprz    ! = number of volume penta_6 elements
        integer::afnhex    ! = number of volume hexa_8 elements
        integer,allocatable,dimension(:)::ibid,ibx,ibxx,ifacetag
	integer,allocatable,dimension(:,:)::if2nt,if2nq,ic2nt,ic2np,ic2nz,ic2nh
	real,allocatable,dimension(:)::x,y,z



open(180,file="grid.ugrid",form='unformatted',status='old',access='stream',convert="big_endian")
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
open(120,file="grid.mapbc",status='old',form='formatted')
read (120,*) nbound

allocate (ibid(nbound),ibx(nbound),ibxx(nbound))
do i=1,nbound
  read(120,*)ibid(i),ibx(i)
end do




do i=1,nbound

	select case(ibx(i))

	case (5000,5050)	!farfield
	ibxx(i)=6
	case(6662,6661,6663)	!symmetry
	ibxx(i)=3
	case(4000)	!wall
	ibxx(i)=4
	case(7031)	!outflow
	ibxx(i)=2

	case(7036,7100)	!inflow
	ibxx(i)=1

	case(6100)	!periodicity
	ibxx(i)=5

	end select

end do








			  !write nodes first
			  open(12,file="GRID.vrt",form='unformatted',action='write')
			 do i=1,afnnodesg
			write(12)i,x(i),y(i),z(i)
! 			 write(1200,"(5x,i8,2x,es21.14,2x,es21.14,2x,es21.14)")i,x(i),y(i),z(i)
			 end do
			 close(12)

			 !end nodes writing


			 !write elements now
			 kx=0
			 open(11,file="GRID.cel",form='unformatted',action='write')

			 !tetra: 1 2 3 3 4 4 4 4
			 do i=1,afntet
			    kx=kx+1
! 			    write(150,"(9i10)")kx,ic2nt(1,i),ic2nt(2,i),ic2nt(3,i),ic2nt(3,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i)
				write(11)kx,ic2nt(1,i),ic2nt(2,i),ic2nt(3,i),ic2nt(3,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i)

			 end do
			 !pyramid: 1 2 3 4 5 5 5 5


			 do i=1,afnpyr
			    kx=kx+1
 			    !write(150,*)kx,ic2np(1,i),ic2np(2,i),ic2np(5,i),ic2np(4,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
! 			    	write(150,"(9i10)")kx,ic2np(1,i),ic2np(4,i),ic2np(5,i),ic2np(2,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
				write(11)kx,ic2np(1,i),ic2np(4,i),ic2np(5,i),ic2np(2,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
			 end do
			 !prism: 1 2 3 3 4 5 6 6





			 do i=1,afnprz
			    kx=kx+1
 			    !write(11,"(9i10)")kx,ic2nz(1,i),ic2nz(2,i),ic2nz(3,i),ic2nz(3,i),ic2nz(4,i),ic2nz(5,i),ic2nz(6,i),ic2nz(6,i)
				write(11)kx,ic2nz(1,i),ic2nz(2,i),ic2nz(3,i),ic2nz(3,i),ic2nz(4,i),ic2nz(5,i),ic2nz(6,i),ic2nz(6,i)
			 end do
			 !hexa: 1 2 3 4 5 6 7 8
			  do i=1,afnhex
			    kx=kx+1
! 			    write(150,"(9i10)")kx,ic2nh(1,i),ic2nh(2,i),ic2nh(3,i),ic2nh(4,i),ic2nh(5,i),ic2nh(6,i),ic2nh(7,i),ic2nh(8,i)
				write(11)kx,ic2nh(1,i),ic2nh(2,i),ic2nh(3,i),ic2nh(4,i),ic2nh(5,i),ic2nh(6,i),ic2nh(7,i),ic2nh(8,i)
			 end do
			 close(11)
			 !end writing elements

			 !now write the boundary file
			 open(10,file="GRID.bnd",form='unformatted',action='write')
			 kx=0
			 !triangle: 1 2 3 3
			  do i=1,afntface
			    kx=kx+1
! 			    write(1000,"(6i12)")kx,if2nt(1,i),if2nt(2,i),if2nt(3,i),if2nt(3,i),ibxx(ifacetag(kx))
			   write(10)kx,if2nt(1,i),if2nt(2,i),if2nt(3,i),if2nt(3,i),ibxx(ifacetag(kx))
			 end do
			 !quad: 1 2 3 4
			  do i=1,afnqface
			    kx=kx+1
! 			    write(1000,"(6i12)")kx,if2nq(1,i),if2nq(2,i),if2nq(3,i),if2nq(4,i),ibxx(ifacetag(kx))
				write(10)kx,if2nq(1,i),if2nq(2,i),if2nq(3,i),if2nq(4,i),ibxx(ifacetag(kx))
			 end do

			close(10)

			 deallocate(ibid,ibx,ibxx, x,y,z, if2nt,if2nq,ifacetag,ic2nt,ic2np,ic2nz,ic2nh)


end subroutine













 end module translate
