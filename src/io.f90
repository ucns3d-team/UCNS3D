module io
use mpiinfo
use declaration
use flow_operations
use iso_c_binding
use transform
implicit none
contains


subroutine outwritegridb
 !> @brief
!> this subroutine writes the grid file in tecplot binary format
use iso_c_binding
implicit none

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
integer,dimension(70)::ivalid
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,kmmg
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
     
if (n.eq.0)then

		inquire (file='GRID.plt',exist=herev)

	if (herev)then



	else

nullptr = 0
      debug   = 0
      filetype = 1
      visdouble = 1
      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 5
      soltime = 360.0
      strandid = 0
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
nulchar = char(0)


ierr =  tecini112('simple dataset'//nulchar, &
                    'x y z'//nulchar, &
                    'GRID.plt'//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)


 ierr= teczne112('grid1'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    null, &
                    null, &
                    shrconn)




	
    allocate(xbin(imaxn))
    allocate(ybin(imaxn))
    allocate(zbin(imaxn))
	if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
        do i=1,imaxn
	read(96,*)j,x,y,z
	xbin(i)=x/scaler
	ybin(i)=y/scaler
 	zbin(i)=z/scaler
	end do
	close(96)
	else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
        do i=1,imaxn
	read(96)j,x,y,z
	xbin(i)=x/scaler
	ybin(i)=y/scaler
 	zbin(i)=z/scaler
	end do
	close(96)
	end if


    ierr = tecdat112(imaxn,xbin,1) 
   
    ierr = tecdat112(imaxn,ybin,1)

     ierr = tecdat112(imaxn,zbin,1)
   
    deallocate(xbin,ybin,zbin)
    
    
    if (binio.eq.0)then
    open(98,file='GRID.cel',form='formatted',status='old',action='read')
	  allocate(icon(8,1))
    icon=0
     cv=0
		do k=1,imaxe
               
 		read(98,*)i,icon(1,1),icon(2,1),icon(3,1),icon(4,1),icon(5,1),icon(6,1),icon(7,1),icon(8,1)
    
		ierr = tecnode112(8,icon)
    !cv=cv+4
        	
		end do
 		close(98)
		!ierr = tecnod112(icon)
		deallocate(icon)	
    else
     open(98,file='GRID.cel',form='unformatted',status='old',action='read')
	  allocate(icon(8,1))
    icon=0
     cv=0
		do k=1,imaxe
               
 		read(98)i,icon(1:8,1)
    
		ierr = tecnode112(8,icon)
    !cv=cv+4
        	
		end do
 		close(98)
		!ierr = tecnod112(icon)
		deallocate(icon)	
    
    
    
    end if
    
    
    
  ierr = tecend112()


end if


end if


	call mpi_barrier(mpi_comm_world,ierror)
	
	
	
	

	
	

end subroutine outwritegridb


subroutine outwritegridb2d
 !> @brief
!> this subroutine writes the grid file in tecplot binary format in 2d
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
integer,dimension(70)::ivalid
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
     
if (n.eq.0)then
inquire (file='GRID.plt',exist=herev)

if (herev)then



else


nullptr = 0
      debug   = 0
      filetype = 1
      visdouble = 1
      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 3
      soltime = 360.0
      strandid = 0
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
nulchar = char(0)


ierr =  tecini112('simple dataset'//nulchar, &
                    'x y'//nulchar, &
                    'GRID.plt'//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)


 ierr= teczne112('grid1'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    null, &
                    null, &
                    shrconn)




	
    allocate(xbin(imaxn))
    allocate(ybin(imaxn))
	if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
        do i=1,imaxn
	read(96,*)j,x,y
	xbin(i)=x/scaler
	ybin(i)=y/scaler
	end do
	close(96)
	else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
        do i=1,imaxn
	read(96)j,x,y
	xbin(i)=x/scaler
	ybin(i)=y/scaler
	end do
	close(96)
	end if


    ierr = tecdat112(imaxn,xbin,1)  !!! why not xbin instead of xbin(1) ??
   
    ierr = tecdat112(imaxn,ybin,1)

    
   
    deallocate(xbin,ybin)
    
    if (binio.eq.0)then
    
    open(98,file='GRID.cel',form='formatted',status='old',action='read')
	  allocate(icon(4,1))
    icon=0
     cv=0
		do k=1,imaxe
               
 		read(98,*)i,icon(1,1),icon(2,1),icon(3,1),icon(4,1)
    
		ierr = tecnode112(4,icon)
    !cv=cv+4
        	
		end do
 		close(98)
		!ierr = tecnod112(icon)
		deallocate(icon)
		
		
    else
	   open(98,file='GRID.cel',form='unformatted',status='old',action='read')
	  allocate(icon(4,1))
    icon=0
     cv=0
		do k=1,imaxe
               
 		read(98)i,icon(1:4,1)
    
		ierr = tecnode112(4,icon)
    !cv=cv+4
        	
		end do
 		close(98)
		!ierr = tecnod112(icon)
		deallocate(icon)
    
    
    
    end if
         
        
  ierr = tecend112()


end if

end if	



	call mpi_barrier(mpi_comm_world,ierror)
	
	

	

	
	

end subroutine outwritegridb2d


subroutine outwrite3n
 !> @brief
!> this subroutine is solely for debugging
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)


kmaxe=xmpielrank(n)



! dumg=kmaxe
! call mpi_barrier(mpi_comm_world,ierror)
! 
! call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
! imaxp=duml

! allocate(icell(imaxe))
! icell=0

! do i=1,kmaxe
!   icell(i)=ielem_ihexgl(i)
! end do


! if (n.eq.0)then
! 	allocate(icella(imaxp*isize))
! 	 icella=0
! 
! end if
! 
! call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)
! 
! ! if (n.eq.0)then

! ! 
! ! 
! ! end if
! 
! call mpi_barrier(mpi_comm_world,ierror)
! deallocate (icell)




if (n.eq.0)then


allocate(xbin(imaxe))
allocate(variables(8))

 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	

end if



if (itestcase.eq.3)then
	nvar1=2
if (n.eq.0)ierr =  tecini112('sols1'//nulchar, &
                    'solution,sols2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
end if


	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 5
!       if (( rungekutta .lt. 5).or.( rungekutta .eq. 11)) then
      soltime = 0.0
!       else
!       soltime = it
! 
!       end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


!  allocate(valuesa(imaxp*isize))
!   allocate(xbin(imaxe))
! 	valuesa=0.0
! 
!  end if
! 
!   call mpi_barrier(mpi_comm_world,ierror)
!   allocate(valuess(imaxp))
!   valuess=0.0
    
   


!   if (itestcase.lt.3)then
!     do i=1,kmaxe
!       valuess(i)=u_c_val(1,1,i)
!     end do
! 
!     call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
! 
!     if (n.eq.0)then
!     do i=1,imaxp*isize
! 	if (icella(i).gt.0)then
! 	xbin(icella(i))=valuesa(i)
! 	end if
!     end do
    xbin(1:imaxe)=xmpie(1:imaxe)
	  

    ierr = tecdat112(imaxe,xbin,1)


    ierr = tecdat112(imaxe,xbin,1)

  ierr = tecend112()
    deallocate(xbin,variables, valuelocation)
    end if

     call mpi_barrier(mpi_comm_world,ierror)
! 
!   end if
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
!   if (itestcase.ge.3)then
! 			  do j=1,nvar1
! 				  
! 				    if (j.eq.1)then
! 				    
! 				    do i=1,kmaxe
! 				      
! 				      valuess(i)=u_c_val(1,1,i)
! 				    end do
! 				    else
! 				    if (j.lt.6)then
! 				    do i=1,kmaxe
! 				      valuess(i)=u_c_val(1,j,i)/u_c_val(1,1,i)
! 				    end do
! 				    end if
! 				    end if
! 				    if (j.eq.6)then
! 
! 				      do i=1,kmaxe
! 											variables(1)=u_c_val(1,1,i)
! 											variables(2)=(u_c_val(1,2,i)/variables(1))
! 											variables(3)=(u_c_val(1,3,i)/variables(1))
! 											variables(4)=(u_c_val(1,4,i)/variables(1))
! 											variables(5)=u_c_val(1,5,i)
! 									    variables(6)=((gamma-1.0))*((variables(5))-0.5*variables(1)*(((variables(2)**2.0)+(variables(3)**2.0)+(variables(4)**2.0))))
! 									    valuess(i)=variables(6)
! 
! 				    end do
! 
! 				    
! 				    end if 
! 
! 				    if (turbulence.ne.1)then
! 					      if (passivescalar.gt.0)then
! 						      if (j.eq.7)then
! 						      do i=1,kmaxe
! 							valuess(i)=u_ct_val(1,1,i)
! 						      end do
! 						      end if
! 					      end if
! 				    end if
!     
! 
! 					  if (turbulence.eq.1)then
! 								if (j.eq.7)then
! 										do i=1,kmaxe
! 														leftv(1:5)=u_c_val(1,1:5,i)
! 														rightv(1:5)=leftv(1:5)
! 										  call sutherlandii(n,leftv,rightv,viscl,laml,pres,rres,gamma,visc,betaas,suther,prandtl)
! 												  variables(1)=viscl(1)/visc
! 												  valuess(i)=variables(1)
! 									      end do
! 								end if
! 								if (j.eq.8)then
! 								  do i=1,kmaxe
! 				      leftv(1:5)=u_c_val(1,1:5,i)
! 				      rightv(1:5)=leftv(1:5)
! 							call sutherlandii(n,leftv,rightv,viscl,laml,pres,rres,gamma,visc,betaas,suther,prandtl)
! 							 eddyfl(1)=ielem_walldist(i)
! 							  if (turbulencemodel.eq.1)then
! 							 eddyfl(2)=u_ct_val(1,1,i)*u_c_val(1,1,i)
! 							 eddyfl(3)=0
! 							  
! 							  end if
! 							  if (turbulencemodel.eq.2)then
! 							 eddyfl(2)=u_ct_val(1,1,i)
! 							 eddyfl(3)=u_ct_val(1,2,i)
! 							  
! 							  end if
! 							eddyfl(4:6)=rec_grads(1,1:3,i)
! 							eddyfl(7:9)=rec_grads(2,1:3,i)
! 							eddyfl(10:12)=rec_grads(3,1:3,i)
! 							if (turbulencemodel.eq.2)then
! 							eddyfl(13:15)=rec_grads(4,1:3,i)
! 							eddyfl(16:18)=rec_grads(5,1:3,i)
! 							end if
! 							eddyfr=eddyfl
! 							
! 							call eddyviscoo(n,viscl,laml,turbmv,etvm,leftv,rightv,eddyfl,eddyfr)
! 							variables(8) =  (viscl(3))/visc	
! 			valuess(i)=variables(8)
! 
! 								  end do
! 								end if
!    
! 							if (turbulencemodel.eq.2)then
! 							    if ((j.eq.9))then
! 							      do i=1,kmaxe
! 								valuess(i)=u_ct_val(1,1,i)
! 							      end do
! 
! 							    end if
! 							    if ((j.eq.10))then
! 							      do i=1,kmaxe
! 								valuess(i)=u_ct_val(1,2,i)
! 							      end do
! 
! 							    end if
! 							end if
! 
! 							if (passivescalar.gt.0)then
! 								      if (j.eq.nvar1-1)then
! 								    do i=1,kmaxe
! 								      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
! 								    end do
! 								      end if
! 							 end if
! 
! 
! 					  end if !turbulence
! 
! 
! 				    if (ivortex.eq.1)then
! 					      if (j.eq.nvar1)then
! 						    do i=1,kmaxe
! 						      valuess(i)=ielem_vortex(i)
! 						    end do
! 					    end if
! 				    end if
! 
! 
! 
! 
!       
!       call mpi_barrier(mpi_comm_world,ierror)
! 
!     
!     call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
! 
!      
! !       print*,"writing output var",j,nvar1,n
! 
!     
!     if (n.eq.0)then
!     do i=1,imaxp*isize
! 	if (icella(i).gt.0)then
! 	xbin(icella(i))=valuesa(i)
! 	end if
!     end do
!     
!     ierr = tecdat112(imaxe,xbin,1)
!     end if
!     
!       
! !       print*,"writing 1r",n
!       
!   end do
! 
! 
! 
!   end if
!   if (n.eq.0)then
!   ierr = tecend112()
!   deallocate(xbin,valuesa,valuelocation,icella)
!   deallocate(out1)
!   end if
!   deallocate (valuess)
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
!   call mpi_barrier(mpi_comm_world,ierror)
! 
! 
! 
! 
! 
! 
! ! 
! !     allocate(xbin(imaxe))
! !     allocate(fbin(1,imaxe))
! ! 
! !  
! !     do i=1,1
! !     fbin(i,:)=10.0
! !     end do
! !     
! !     do i=1,1
! !     xbin(:)=fbin(i,:)
! ! 
! !     ierr = tecdat112(imaxe,xbin,1)  !!! why not xbin instead of xbin(1) ??
! !    end do
! 
! 		
!          
!         
!  deallocate(variables)









	
	
	
	
	

	
	

end subroutine outwrite3n



subroutine movie
 !> @brief
!> this subroutine writes only the 3d solution without the grid in tecplot binary format
use iso_c_binding
implicit none
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
 character(len=:),allocatable::out1
 character*1 nulchar

      integer::   debug,iii,npts,nelm


      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)









allocate(variables(2))

kmaxe=xmpielrank(n)






if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="mov_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)


end if
call mpi_barrier(mpi_comm_world,ierror)


 if (itestcase.eq.4)then
 nvar1=2



              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'v_mag,q'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)




end if




if (n.eq.0)then
allocate (valuelocation(nvar1))


      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 5

      soltime = t

      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)



  allocate(xbin(imaxe),xbin2(imaxe))


 else
 allocate(xbin2(1))
 end if



  allocate(valuess(kmaxe))

!   call mpi_barrier(mpi_comm_world,ierror)









		    do i=1,kmaxe


            leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)

		    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

		    valuess(i)=sqrt(leftv(2)**2+leftv(3)**2+leftv(4)**2)

			end do

		    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
					if (n.eq.0)then
					do i=1,imaxe
				xbin(xmpi_re(i))=xbin2(i)
				end do
					ierr = tecdat112(imaxe,xbin,1)
					end if





		do i=1,kmaxe

		  valuess(i)=ielem_vortex(1,i)
		end do


		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if












  if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1)
  end if

  deallocate (valuess,variables,xbin2)











end subroutine movie





subroutine outwrite3vb
 !> @brief
!> this subroutine writes only the 3d solution without the grid in tecplot binary format
use iso_c_binding
implicit none
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
 character(len=:),allocatable::out1
 character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)


 






allocate(variables(14))

kmaxe=xmpielrank(n)






if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	

end if
call mpi_barrier(mpi_comm_world,ierror)
 if (itestcase.le.2)then
  nvar1=4
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'solution1,solution2,sten1,sten2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=8+passivescalar
	  if (passivescalar.gt.0)then
			if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
							'density,u,v,w,energy,pressure,sten1,sten2,passivescalar'//nulchar, &
							out1//nulchar, &
							'.'//nulchar, &
							filetype, &
							debug, &
							visdouble)
     else
     if (multispecies.eq.1)then
     nvar1=10

				if (dg.eq.1)then
						nvar1=10
						if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
										'density,u,v,w,energy,pressure,species1,species2,vfraction,troubled'//nulchar, &
										out1//nulchar, &
										'.'//nulchar, &
										filetype, &
										debug, &
										visdouble)



				else

						if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
										'density,u,v,w,energy,pressure,species1,species2,vfraction,aux'//nulchar, &
										out1//nulchar, &
										'.'//nulchar, &
										filetype, &
										debug, &
										visdouble)


				end if

     
     else	!no multispecies
			if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
							'density,u,v,w,energy,pressure,sten1,sten2'//nulchar, &
							out1//nulchar, &
							'.'//nulchar, &
							filetype, &
							debug, &
							visdouble)






	end if




     end if
 end if
 if (itestcase.eq.4)then
 nvar1=9+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,sten1,sten2,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,sten1,sten2,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,sten1,sten2,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if

	

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))


      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 5
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0

 



 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(imaxe),xbin2(imaxe))
	

 else
 allocate(xbin2(1))
 end if



  allocate(valuess(kmaxe))

!   call mpi_barrier(mpi_comm_world,ierror)
 
    
   


    if (itestcase.le.2)then
		do i=1,kmaxe
		 
     valuess(i)=u_c_val(1,1,i)!0.0
    
		end do
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if

    
		do i=1,kmaxe
		
      valuess(i)=ielem_inumneighbours(i)
     
		end do
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
    
		do i=1,kmaxe
		  valuess(i)=ielem_troubled(i)
		end do
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if

    
		do i=1,kmaxe
		  valuess(i)=ielem_admis(i)
		end do
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
		 if (n.eq.0)then
		 do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
		    do i=1,kmaxe
		    
            
            leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
            
		    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		      valuess(i)=leftv(kkd)
			if (kkd.eq.5)then
			
            valuess(i)=u_c_val(1,kkd,i)
           
			end if
		    end do
		    
		    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if
		end do
		
		
		
    
		do i=1,kmaxe
		      
            leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
            
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		end do
		
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if

		
		
		
		
		
                if (multispecies.eq.1)then
                
                do i=1,kmaxe
                valuess(i)=u_c_val(1,6,i)
                end do
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
                    if (n.eq.0)then
                    do i=1,imaxe
                xbin(xmpi_re(i))=xbin2(i)
                end do
                    ierr = tecdat112(imaxe,xbin,1)
                    end if
			
			
                    do i=1,kmaxe
                        valuess(i)=u_c_val(1,7,i)
                        end do
                    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
                    if (n.eq.0)then
                    do i=1,imaxe
                xbin(xmpi_re(i))=xbin2(i)
                end do
                    ierr = tecdat112(imaxe,xbin,1)
                    end if
                    
                    
                       do i=1,kmaxe
                valuess(i)=u_c_val(1,8,i)
                end do
                
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if
                    

                    if (dg.eq.1)then
                     do i=1,kmaxe
                valuess(i)=ielem_troubled(i)
                end do

                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if
                else
                        do i=1,kmaxe
                valuess(i)=ielem_reduce(i)
                end do

                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
                        if (n.eq.0)then
                        do i=1,imaxe
                xbin(xmpi_re(i))=xbin2(i)
                end do
                        ierr = tecdat112(imaxe,xbin,1)
                        end if







                end if





                    
                    
			else
			
			
                
                if (mood.eq.1)then
                do i=1,kmaxe
                valuess(i)=ielem_mood_o(i)
                end do
                else
                do i=1,kmaxe
                if (adda.eq.1)then

                valuess(i)=ielem_lwcx2(i)!diss!ielem_stencil_dist(i)

                else
                valuess(i)=ielem_ggs(i)!wcx(1)!troubled!filtered

                end if
                end do
                end if
                
                
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if
            
                
               
                
                
                do i=1,kmaxe
                
                if (adda.eq.1)then
                
                valuess(i)=ielem_diss(i)
                else
                
                
                valuess(i)=ielem_full(i)!wcx(1)!ielem_admis(i)
                end if
                end do
                
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if

              end if  
		
    
		  if (passivescalar.gt.0)then
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
		  end do
		  
		  
		  call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  do i=1,kmaxe
		      valuess(i)=ielem_vortex(1,i)!%inumneighbours
		  end do
		  
		  
		 call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			if (n.eq.0)then
			do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			ierr = tecdat112(imaxe,xbin,1)
			end if
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
			do i=1,kmaxe
			    valuess(i)=u_ct_val(1,kkd,i)
			end do
		      
			
		      call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
			      if (n.eq.0)then
			      do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
			      ierr = tecdat112(imaxe,xbin,1)
			      end if


		 end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
    

     
    
    
    
  if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1)
  end if
  
  deallocate (valuess,variables,xbin2)
  

  





	
	

end subroutine outwrite3vb


subroutine outwritetec3dbp
 !> @brief
!> this subroutine writes only the 3d solution without the grid in tecplot binary format
use iso_c_binding
implicit none
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2,xbin3
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
 character(len=:),allocatable::out1
 character*1 nulchar
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)


 






allocate(variables(15),icon(8,1))

kmaxe=xmpielrank(n)





 nullptr = 0
      debug   = 0
      filetype = 0
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') n
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	
	
	
	
	
 if (itestcase.le.2)then
  nvar1=7
  ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,solution1,solution2,sten1,sten2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=11+passivescalar
  if (passivescalar.gt.0)then
 ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     if (multispecies.eq.1)then
     nvar1=12
      ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,species1,species2,vfraction'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     else
    ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     end if
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=12+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if

	

	

allocate (valuelocation(nvar1))


      imax    = kmaxn
      jmax    = kmaxe
      kmax    = 0
      ZONETYPE = 5
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0
valuelocation(1:3)=1
 



 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(kmaxn),xbin2(kmaxn),xbin3(kmaxn))
	

 
  allocate(valuess(kmaxe))


    do i=1,kmaxn
        xbin(i)=inoder4_cord(1,i);
        xbin2(i)=inoder4_cord(2,i);
        xbin3(i)=inoder4_cord(3,i)
        
    end do
    
    ierr = tecdat112(kmaxn,xbin,1) 
   
    ierr = tecdat112(kmaxn,xbin2,1)

     ierr = tecdat112(kmaxn,xbin3,1)

     
    
     

    if (itestcase.le.2)then
		do i=1,kmaxe
		  valuess(i)=u_c_val(1,1,i)!0.0
		end do
		
		ierr = tecdat112(kmaxe,valuess,1)
		

    
		do i=1,kmaxe
		  valuess(i)=n
		end do
		
		
		
		ierr = tecdat112(kmaxe,valuess,1)
		
    
		do i=1,kmaxe
		  valuess(i)=ielem_inumneighbours(i)!%stencil_dist
		end do
		
		
		ierr = tecdat112(kmaxe,valuess,1)
		

    
		do i=1,kmaxe
		  valuess(i)=ielem_admis(i)
		end do
		
		
		ierr = tecdat112(kmaxe,valuess,1)
		
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
		    do i=1,kmaxe
		    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		      valuess(i)=leftv(kkd)
			if (kkd.eq.5)then
			valuess(i)=u_c_val(1,kkd,i)!/u_c_val(1,1,i)
			end if
		    end do
		    
		    
		ierr = tecdat112(kmaxe,valuess,1)
			
		end do
		
		
		
    
		do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		end do
		
		
		
			ierr = tecdat112(kmaxe,valuess,1)
			

		
		
		
		
		
                if (multispecies.eq.1)then
                
                do i=1,kmaxe
                valuess(i)=u_c_val(1,6,i)
                end do
                
                 ierr = tecdat112(kmaxe,valuess,1)
                   
			
			
                    do i=1,kmaxe
                        valuess(i)=u_c_val(1,7,i)
                        end do
                    
                ierr = tecdat112(kmaxe,valuess,1)
                    
                    
                    
                       do i=1,kmaxe
                valuess(i)=u_c_val(1,8,i)
                end do
                
               
			ierr = tecdat112(kmaxe,valuess,1)
			
                    
                    
                    
			else
			
			
                
                if (mood.eq.1)then
                do i=1,kmaxe
                valuess(i)=ielem_mood_o(i)
                end do
                else
!                 do i=1,kmaxe
!                 valuess(i)=ielem_condition(i)!ielem_stencil_dist(i)
!                 end do
                end if
                
                
                
			ierr = tecdat112(kmaxe,valuess,1)
			
            
                
               
                
                
                do i=1,kmaxe
                valuess(i)=ielem_admis(i)
                end do
                
                
			ierr = tecdat112(kmaxe,valuess,1)
			

              end if  
		
    
		  if (passivescalar.gt.0)then
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
		  end do
		  
		  
		  
			ierr = tecdat112(kmaxe,valuess,1)
			
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  do i=1,kmaxe
		      valuess(i)=ielem_vortex(1,i)!%inumneighbours
		  end do
		  
		  
		
		ierr = tecdat112(kmaxe,valuess,1)
			
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
			do i=1,kmaxe
			    valuess(i)=u_ct_val(1,kkd,i)
			end do
		      
			
		      
		ierr = tecdat112(kmaxe,valuess,1)
			      


		 end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    

    do i=1,kmaxe
    icon(1:8,1)=el_connect(i,1:8)
    ierr = tecnode112(8,icon)
    end do
    
  
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1,xbin2,xbin3,icon)
  
  
  deallocate (valuess,variables)
  

  





	
	

end subroutine outwritetec3dbp


subroutine outwritetec3dbpav
 !> @brief
!> this subroutine writes only the 3d solution without the grid in tecplot binary format
use iso_c_binding
implicit none
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2,xbin3
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
 character(len=:),allocatable::out1
 character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)


 






allocate(variables(15),icon(1,8))

kmaxe=xmpielrank(n)





 nullptr = 0
      debug   = 0
      filetype = 0
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') n
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	
 if (itestcase.le.2)then
  nvar1=7
  ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,solution1,solution2,sten1,sten2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=11+passivescalar
  if (passivescalar.gt.0)then
 ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     if (multispecies.eq.1)then
     nvar1=12
      ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,species1,species2,vfraction'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     else
    ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     end if
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=12+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,passivescalar,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              ierr =  tecini112('sols'//nulchar, &
                    'x,y,z,density,u,v,w,energy,pressure,sten1,sten2,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if

	

	

allocate (valuelocation(nvar1))


      imax    = kmaxn
      jmax    = kmaxe
      kmax    = 0
      ZONETYPE = 5
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0
valuelocation(1:3)=1
 



 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(kmaxn),xbin2(kmaxn),xbin3(kmaxn))
	

 
  allocate(valuess(kmaxe))


    do i=1,kmaxn
        xbin=inoder4_cord(1,i);
        xbin2=inoder4_cord(2,i);
        xbin3=inoder4_cord(3,i)
    end do
    
    ierr = tecdat112(kmaxn,xbin,1) 
   
    ierr = tecdat112(kmaxn,xbin2,1)

     ierr = tecdat112(kmaxn,xbin3,1)


    if (itestcase.le.2)then
		do i=1,kmaxe
		  valuess(i)=u_c_val(1,1,i)!0.0
		end do
		
		ierr = tecdat112(kmaxe,valuess,1)
		

    
		do i=1,kmaxe
		  valuess(i)=n
		end do
		
		
		
		ierr = tecdat112(kmaxe,valuess,1)
		
    
		do i=1,kmaxe
		  valuess(i)=ielem_inumneighbours(i)!%stencil_dist
		end do
		
		
		ierr = tecdat112(kmaxe,valuess,1)
		

    
		do i=1,kmaxe
		  valuess(i)=ielem_admis(i)
		end do
		
		
		ierr = tecdat112(kmaxe,valuess,1)
		
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
		    do i=1,kmaxe
		    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		      valuess(i)=leftv(kkd)
			if (kkd.eq.5)then
			valuess(i)=u_c_val(1,kkd,i)!/u_c_val(1,1,i)
			end if
		    end do
		    
		    
		ierr = tecdat112(kmaxe,valuess,1)
			
		end do
		
		
		
    
		do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		end do
		
		
		
			ierr = tecdat112(kmaxe,valuess,1)
			

		
		
		
		
		
                if (multispecies.eq.1)then
                
                do i=1,kmaxe
                valuess(i)=u_c_val(1,6,i)
                end do
                
                 ierr = tecdat112(kmaxe,valuess,1)
                   
			
			
                    do i=1,kmaxe
                        valuess(i)=u_c_val(1,7,i)
                        end do
                    
                ierr = tecdat112(kmaxe,valuess,1)
                    
                    
                    
                       do i=1,kmaxe
                valuess(i)=u_c_val(1,8,i)
                end do
                
               
			ierr = tecdat112(kmaxe,valuess,1)
			
                    
                    
                    
			else
			
			
                
                if (mood.eq.1)then
                do i=1,kmaxe
                valuess(i)=ielem_mood_o(i)
                end do
                else
!                 do i=1,kmaxe
!                 valuess(i)=ielem_condition(i)!ielem_stencil_dist(i)
!                 end do
                end if
                
                
                
			ierr = tecdat112(kmaxe,valuess,1)
			
            
                
               
                
                
                do i=1,kmaxe
                valuess(i)=ielem_admis(i)
                end do
                
                
			ierr = tecdat112(kmaxe,valuess,1)
			

              end if  
		
    
		  if (passivescalar.gt.0)then
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
		  end do
		  
		  
		  
			ierr = tecdat112(kmaxe,valuess,1)
			
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  do i=1,kmaxe
		      valuess(i)=ielem_vortex(1,i)!%inumneighbours
		  end do
		  
		  
		
		ierr = tecdat112(kmaxe,valuess,1)
			
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
			do i=1,kmaxe
			    valuess(i)=u_ct_val(1,kkd,i)
			end do
		      
			
		      
		ierr = tecdat112(kmaxe,valuess,1)
			      


		 end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    do i=1,kmaxe
    icon(1,1:8)=el_connect(i,1:8)
    ierr = tecnode112(8,icon)
    end do
    
  
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1,xbin2,xbin3,icon)
  
  
  deallocate (valuess,variables)
  

  





	
	

end subroutine outwritetec3dbpav


subroutine outwrite3v
!> @brief
!> this subroutine writes only the 3d solution without the grid in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(10))

kmaxe=xmpielrank(n)

dumg=kmaxe
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then
! write(1000+n,*)icella(:)
! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
		open(97,file=outfile,form='formatted',status='new',action='write')

end if
call mpi_barrier(mpi_comm_world,ierror)
 if (itestcase.le.2)then
  nvar1=1
  if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="solution"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
  
  if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if
  
      else
     
     if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
   end if
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=7+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","K","OMEGA"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED, [9] = CELLCENTERED,[10] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","MU"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED, [9] = CELLCENTERED)'
	  end if
              
              end if
              if (turbulenceequations.eq.0)then
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	                if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","K","OMEGA"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                              if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","MU"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
             
              end if
              if (turbulenceequations.eq.0)then
                              if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
	  end if
              
              end if
	    
	    
	    
	    
	    end if
	   
 end if

	if (n.eq.0)then
	    write(97,*) ', solutiontime=',t
	    end if

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







valuelocation(:)=0

 





 allocate(valuesa(imaxp*isize))
  allocate(xbin(imaxe))
	valuesa=zero

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=zero
    
   


    if (itestcase.le.2)then
    do i=1,kmaxe
      valuess(i)=u_c_val(1,1,i)!0.0
    end do
    
    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    
			write(97,*)xbin(1:imaxe)
			
    
    end if

     
    call mpi_barrier(mpi_comm_world,ierror)
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
		do i=1,kmaxe
		  valuess(i)=u_c_val(1,kkd,i)
		  if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(i)=u_c_val(1,kkd,i)/u_c_val(1,1,i)
		  end if
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
		do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
    
		  if (passivescalar.gt.0)then
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:imaxe)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  do i=1,kmaxe
		      valuess(i)=ielem_vortex(1,i)
		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:imaxe)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,kkd,i)
		  end do
		
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:imaxe)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
    

     
    
    
    
  if (n.eq.0)then
  close(97)
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)







         
        
 deallocate(variables)


	
	

end subroutine outwrite3v


subroutine outwrite3v2d
!> @brief
!> this subroutine writes only the 2d solution without the grid in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd, i_dof
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(10))

kmaxe=xmpielrank(n)

dumg=kmaxe
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
		open(97,file=outfile,form='formatted',status='new',action='write')

end if
call mpi_barrier(mpi_comm_world,ierror)
 if (itestcase.le.2)then
  nvar1=1
  if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="solution"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=5+passivescalar
  if (passivescalar.gt.0)then
  
  if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
  
      else
     
     if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED)'
   end if
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=7+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","K","OMEGA"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED, [9] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","MU"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
              
              end if
              if (turbulenceequations.eq.0)then
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
	  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	                if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","K","OMEGA"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED)'
	  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                              if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","MU"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
	  end if
             
              end if
              if (turbulenceequations.eq.0)then
                              if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
	  end if
              
              end if
	    
	    
	    
	    
	    end if
	   
 end if

	if (n.eq.0)then
	    write(97,*) ', solutiontime=',t
	    end if

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







valuelocation(:)=0

 





 allocate(valuesa(imaxp*isize))
  allocate(xbin(imaxe))
	valuesa=zero

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=zero
    
   


    if (itestcase.le.2)then
        if (dg == 1) then
!         do i_dof = 1, ielem_idegfree(iconsidered) + 1
            do i=1,kmaxe
            valuess(i)=u_c_val(1,1,i)!u_c_val(1,1,i)%valdg(1,1,1)
            end do
!         end do
            
        else
         do i=1,kmaxe
            valuess(i)=u_c_val(1,1,i)!0.0
            end do
        end if
        
        call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

        if (n.eq.0)then
            do i=1,imaxp*isize
                if (icella(i).gt.0)then
                    xbin(icella(i))=valuesa(i)
                end if
            end do
        
            write(97,*)xbin(1:imaxe)
        end if

        call mpi_barrier(mpi_comm_world,ierror)
    else if (itestcase.ge.3)then
		do kkd=1,4
		do i=1,kmaxe
		  valuess(i)=u_c_val(1,kkd,i)
		  if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(i)=u_c_val(1,kkd,i)/u_c_val(1,1,i)
		  end if
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
		do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(4)
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
    
		  if (passivescalar.gt.0)then
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:imaxe)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  do i=1,kmaxe
		      valuess(i)=ielem_vortex(1,i)
		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:imaxe)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,kkd,i)
		  end do
		
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:imaxe)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
    

     
    
    
    
  if (n.eq.0)then
  close(97)
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)







         
        
 deallocate(variables)


	
	

end subroutine outwrite3v2d

subroutine outwrite3vb2d
!> @brief
!> this subroutine writes only the 2d solution without the grid in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(10))

kmaxe=xmpielrank(n)






if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	

end if

 if (itestcase.le.2)then
  nvar1=4
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'sol1,sol2,sten1,sten2'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=8+passivescalar
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,sten1,sten2,slope,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     if (multispecies.eq.1)then
     nvar1=9
     if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,species1,species2,vf,aux'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     else
     if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,slope'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     end if
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=9+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,passivescalar,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,passivescalar,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,passivescalar,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,vortex,k,omega'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,vortex,mu'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                     'density,u,v,energy,pressure,sten1,sten2,vortex'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if

	
	
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 3
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)

allocate(xbin(1:imaxe),xbin2(1:imaxe))
 else
 allocate(xbin2(1))
 end if
  


 
 
  allocate(valuess(1:kmaxe))
  
    
   


    if (itestcase.le.2)then
    do i=1,kmaxe
     valuess(i)=u_c_val(1,1,i)!0.0
    end do
    
    call mpi_gatherv(valuess(1:kmaxe),kmaxe,mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

    
    
    
	      
    
    
		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		
		ierr = tecdat112(imaxe,xbin,1)
		end if
    

    
    
    do i=1,kmaxe
        
      valuess(i)=ielem_troubled(i)
      
    end do
    
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
    
    
    do i=1,kmaxe
      valuess(i)=ielem_wcx(i)!ielem_admis(i)
    end do
    
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
    
    if (initcond.eq.0)then
!      do i=1,kmaxe
!       valuess(i)=rec_cond(1,i)
!     end do
    else
    do i=1,kmaxe
      valuess(i)=ielem_wcx(i)!ielem_stencil_dist(i)
    end do
    end if
    
   call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
    
    
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,4
		do i=1,kmaxe
        
        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
        
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(kkd)
          if (kkd.eq.4)then
          
        
        valuess(i)=u_c_val(1,kkd,i)
        
          
          end if
		end do
		
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		end do
    
		do i=1,kmaxe
		 
        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
        
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(4)
		end do
		
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
    
                if (multispecies.eq.1)then
                do i=1,kmaxe
                valuess(i)=u_c_val(1,5,i)
                end do
                else
                do i=1,kmaxe
                valuess(i)=ielem_full(i)!admis
                end do
                end if
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
                if (multispecies.eq.1)then
                do i=1,kmaxe
                valuess(i)=u_c_val(1,6,i)
                end do
                else
                if (mood.eq.1)then
                do i=1,kmaxe
                valuess(i)=ielem_mood_o(i)
                end do
                else
                
                do i=1,kmaxe
                valuess(i)=ielem_troubled(i)
                end do
                end if
                end if
                
                
               call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if

                
                	
		
		
    
		  if (passivescalar.gt.0)then
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,i)
		  end do
		  
		  
		 call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
		  
		  
		   if (itestcase.eq.3)then
                            if (multispecies.eq.1)then
                do i=1,kmaxe
                valuess(i)=u_c_val(1,7,i)
                end do
                else
                do i=1,kmaxe
                valuess(i)=ielem_wcx(i)!vortex(1)
                end do
                end if
                            
                            
                            
                            call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		
		
		if (multispecies.eq.1)then
		if (mood.eq.1)then
                do i=1,kmaxe
                valuess(i)=ielem_mood_o(i)
                end do
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
                if (n.eq.0)then
                do i=1,imaxe
                xbin(xmpi_re(i))=xbin2(i)
                end do
                ierr = tecdat112(imaxe,xbin,1)
                end if
        else
                do i=1,kmaxe
                valuess(i)=ielem_reduce(i)
                end do
                call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)
                if (n.eq.0)then
                do i=1,imaxe
                xbin(xmpi_re(i))=xbin2(i)
                end do
                ierr = tecdat112(imaxe,xbin,1)
                end if
        
        
        
        end if
        end if
		  end if
    
		  if (itestcase.eq.4)then
                            do i=1,kmaxe
                                valuess(i)=ielem_ggs(i)!vortex(1)
                            end do
                            
                            
                            call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  do i=1,kmaxe
		      valuess(i)=u_ct_val(1,kkd,i)
		  end do
		
		  
		  call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		  end do
		  end if
		  
		  
		  
		  
		  end if
    
    
    end if
    
    
    
!     read_ucns

     
    
    
    
  if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1)
  end if
  
  deallocate (valuess,variables,xbin2)


	
	

end subroutine outwrite3vb2d


subroutine checkres
!> @brief
!> this subroutine checks the presence of restart file
implicit none
logical::here
integer::i,j,k,l,iter,dip
 character(len=20)::proc,restfile

 	restfile='RESTART.dat'
	call mpi_barrier(mpi_comm_world,ierror)
	inquire (file=restfile,exist=here)
	if (here) then
	open(1083+n,file=restfile,form='unformatted',status='old',action='read',access='stream')
	dip=1
	if (ires_unsteady.eq.1)then
	read(1083+n,pos=dip)iter


	dip=dip+4
	read(1083+n,pos=dip)res_time
	dip=dip+8
	    if (initcond.eq.95)then
	    read(1083+n,pos=dip)taylor
	    end if
	else
	  read(1083+n,pos=dip)iter
	end if
	restart=iter
	close(1083+n)
	else
	restart=0
        average_restart=0
	res_time=0.0d0
	end if


	if (n.eq.0)then
	print*,"restarting",iter,res_time,restart

	end if
	

	call mpi_barrier(mpi_comm_world,ierror)

end subroutine checkres


subroutine open_arbitrary(n,imaxe,imaxn,imaxb)
!> @brief
!> this subroutine opens the grid files and establishes the number of cells, nodes, boundary conditions
	implicit none
	integer,intent(inout)::imaxe,imaxn,imaxb
	integer::i,j,k,ios,iox,ioz,i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,ioy
	real::ix1,ix2,ix3
	character(len=12)::proc,vrtfile,celfile,bndfile
	integer,allocatable,dimension(:)::isent
	integer,intent(in)::n
		write(proc,fmt='(i10)') n
		celfile='GRID.cel'
		vrtfile='GRID.vrt'
		bndfile='GRID.bnd'

	allocate(isent(3))
	
	if (n.eq.0)then
	
	if(binio.eq.0)then
	
		open(8,file=celfile,form='formatted',status='old',action='read',iostat=ios)
		open(9,file=vrtfile,form='formatted',status='old',action='read',iostat=iox)
		open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
	i=0
	j=0
	k=0

	
	do 
		read(8,*,iostat=ios)i
		if (ios.ne.0) then
			exit
		end if
		
	end do
		imaxe=i
	do 
		read(9,*,iostat=iox)j
		if (iox.ne.0) then
			exit
		end if
		
	end do
		imaxn=j
	do 
		read(10,*,iostat=ioy)k
		if (ioy.ne.0) then
			exit
		end if
 		
	end do
	 		imaxb=k
	else
	
	
		open(8,file=celfile,form='unformatted',status='old',action='read',iostat=ios)
		open(9,file=vrtfile,form='unformatted',status='old',action='read',iostat=iox)
		open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)
	i=0
	j=0
	k=0

	if (dimensiona.eq.3)then
	do 
		read(8,iostat=ios)i,i1,i2,i3,i4,i5,i6,i7,i8
		if (ios.ne.0) then
			exit
		end if
		
	end do
	
		imaxe=i
		
	do 
		read(9,iostat=iox)j,ix1,ix2,ix3
		if (iox.ne.0) then
			exit
		end if
		
	end do
		imaxn=j
	do 
		read(10,iostat=ioy)k,i1,i2,i3,i4,i5
		if (ioy.ne.0) then
			exit
		end if
 		
	end do
	 		imaxb=k
	end if
	
	if (dimensiona.eq.2)then
	do 
		read(8,iostat=ios)i,i1,i2,i3,i4
		if (ios.ne.0) then
			exit
		end if
		
	end do
	
		imaxe=i
		
	do 
		read(9,iostat=iox)j,ix1,ix2
		if (iox.ne.0) then
			exit
		end if
		
	end do
		imaxn=j
	do 
		read(10,iostat=ioy)k,i1,i2,i3,i4,i5
		if (ioy.ne.0) then
			exit
		end if
 		
	end do
	 		imaxb=k
	end if
	
	
	end if

	close(8)
	close(9)
 	close(10)

	isent(1)=imaxe
	isent(2)=imaxn
	isent(3)=imaxb

	

	end if
	call mpi_bcast( isent, 3, mpi_integer, 0 ,mpi_comm_world,ierror )

	call mpi_barrier(mpi_comm_world,ierror)
	imaxe=isent(1)
	imaxn=isent(2)
	imaxb=isent(3)
	deallocate(isent)
	
	
end subroutine open_arbitrary


subroutine open_input1(n,itt)
!> @brief
!> this subroutine opens the parameter file
	implicit none
	integer,intent(inout)::itt
	integer,intent(in)::n
	character(len=12)::vrtfile,celfile
! 		celfile='GRID.cel'
! 		vrtfile='GRID.vrt'
! 		if (binio.eq.0)then
! 		open(8,file=celfile,form='formatted',status='old',action='read')
! 		open(9,file=vrtfile,form='formatted',status='old',action='read')
! 		else
		open(15,file='UCNS3D.DAT',form='formatted',status='old',action='read')
end subroutine open_input1

subroutine close_input1(n,itt)
!> @brief
!> this subroutine closes the parameter file
 implicit none
 	integer,intent(inout)::itt
 	integer,intent(in)::n	
!  	close(8)
!  	close(9)
 	close(15)
end subroutine close_input1

subroutine open_input(n,itt)
!> @brief
!> this subroutine opens the grid files
	implicit none
	integer,intent(inout)::itt
	integer,intent(in)::n
	character(len=12)::vrtfile,celfile
		celfile='GRID.cel'
		vrtfile='GRID.vrt'
		if (binio.eq.0)then
		open(8,file=celfile,form='formatted',status='old',action='read')
		open(9,file=vrtfile,form='formatted',status='old',action='read')
		else
		open(8,file=celfile,form='unformatted',status='old',action='read')
		open(9,file=vrtfile,form='unformatted',status='old',action='read')
		end if
	end subroutine open_input

subroutine close_input(n,itt)
!> @brief
!> this subroutine closes the grid files
 implicit none
 	integer,intent(inout)::itt
 	integer,intent(in)::n	
  	close(8)
  	close(9)

 end subroutine close_input




subroutine read_input(n,xmpielrank,xmpinrank,xmpie,xmpin,dinode,imaxn,imaxe,imaxb,xmpinnumber,scaler,dinoder)
!> @brief
!> this subroutine reads all the data from the grid files
	implicit none
	type(anode_ne),allocatable,dimension(:),intent(inout)::dinoder
	type(anode_number),allocatable,dimension(:,:),intent(inout)::dinode
	integer,allocatable,dimension(:,:,:),intent(inout)::xmpinnumber
	integer,intent(in)::n,imaxn,imaxe
	integer,allocatable,dimension(:),intent(in)::xmpie,xmpin
	integer,allocatable,dimension(:),intent(in)::xmpielrank,xmpinrank
	integer,intent(inout)::imaxb
	real,intent(inout)::scaler
	integer,allocatable,dimension(:)::nodep,nodec
	integer::i,j,ji,k,lm,iex,kmaxn,kk,kmaxe,print_out,kk2,shap,nodal,iftrue,kxk2
	integer,dimension(8)::idv
	integer::it1,it2,it3,it4,it5,it6,it7,it8,itx,inx,it55,in,out
	real::x,y,z

	

	kk=0; i=0; j=0; lm=0 ;kk2=1
	xmin(n)=tolbig; ymin(n)=tolbig; zmin(n)=tolbig
	xmax(n)=-tolbig; ymax(n)=-tolbig; zmax(n)=-tolbig
	call open_input(n,itt)
        !open(82,file="grid.epart",form='formatted',status='old',action='read')
	kmaxe=xmpielrank(n)
	
	
 	if (dimensiona.eq.3)then
	

	allocate(dinoder(imaxn))
	allocate(dinoder2(imaxn))
	dinoder2(1:imaxn)%numberofneib=0
	dinoder(1:imaxn)%itor=0
	
	
	if (binio.eq.0)then
	
	do j=1,imaxe
	
	if (xmpie(j).eq.n)then

	kk=kk+1
	
	

	read(8,*) itx,idv(1),idv(2),idv(3),idv(4),idv(5),idv(6),idv(7),idv(8)

	    
		
				
	
	  do kxk2=1,8
	  dinoder(idv(kxk2))%itor=idv(kxk2)
	 
	  end do
	  
	  
	  
	  
	  ielem_ihexgl(kk)=itx
		if ((idv(5).ne.idv(6)).and.(idv(6).ne.idv(7)).and.(idv(7).ne.idv(8)).and.(idv(3).ne.idv(4)))then
		shap=1	;nodal=8
		ielem_ishape(kk)=shap
		ielem_ifca(kk)=6
		ielem_nonodes(kk)=nodal
		ielem_nodes(1:nodal,kk)=idv(1:nodal)
		end if
		if ((idv(3).eq.idv(4)).and.(idv(5).eq.idv(6)).and.(idv(6).eq.idv(7)).and.(idv(7).eq.idv(8)))then
		shap=2	;nodal=4!tetrahedral element
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=4
		ielem_nodes(1:3,kk)=idv(1:3)
		ielem_nodes(4,kk)=idv(5)
		end if
		if ((idv(3).ne.idv(4)).and.(idv(5).eq.idv(6)).and.(idv(6).eq.idv(7)).and.(idv(7).eq.idv(8)))then
		shap=3	;nodal=5!pyramidal element
		ielem_ishape(kk)=shap
		ielem_ifca(kk)=5
		ielem_nonodes(kk)=nodal
		ielem_nodes(1:nodal,kk)=idv(1:nodal)
		end if
		if ((idv(3).eq.idv(4)).and.(idv(5).ne.idv(6)).and.(idv(7).eq.idv(8)))then
		shap=4	;nodal=6!prism element
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=5
		ielem_nodes(1:3,kk)=idv(1:3)
		ielem_nodes(4:6,kk)=idv(5:7)
		end if
		  
		
			
	
	
	else
	read(8,*)


	end if 
	end do
	else
	do j=1,imaxe
	
	if (xmpie(j).eq.n)then

	kk=kk+1
	
	

	read(8) itx,idv(1),idv(2),idv(3),idv(4),idv(5),idv(6),idv(7),idv(8)

	    
		
				
	
	  do kxk2=1,8
	  dinoder(idv(kxk2))%itor=idv(kxk2)
	 
	  end do
	  ielem_ihexgl(kk)=itx
		if ((idv(5).ne.idv(6)).and.(idv(6).ne.idv(7)).and.(idv(7).ne.idv(8)).and.(idv(3).ne.idv(4)))then
		shap=1	;nodal=8
		ielem_ishape(kk)=shap
		ielem_ifca(kk)=6
		ielem_nonodes(kk)=nodal
		ielem_nodes(1:nodal,kk)=idv(1:nodal)
		end if
		if ((idv(3).eq.idv(4)).and.(idv(5).eq.idv(6)).and.(idv(6).eq.idv(7)).and.(idv(7).eq.idv(8)))then
		shap=2	;nodal=4!tetrahedral element
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=4
		ielem_nodes(1:3,kk)=idv(1:3)
		ielem_nodes(4,kk)=idv(5)
		end if
		if ((idv(3).ne.idv(4)).and.(idv(5).eq.idv(6)).and.(idv(6).eq.idv(7)).and.(idv(7).eq.idv(8)))then
		shap=3	;nodal=5!pyramidal element
		ielem_ishape(kk)=shap
		ielem_ifca(kk)=5
		ielem_nonodes(kk)=nodal
		ielem_nodes(1:nodal,kk)=idv(1:nodal)
		end if
		if ((idv(3).eq.idv(4)).and.(idv(5).ne.idv(6)).and.(idv(7).eq.idv(8)))then
		shap=4	;nodal=6!prism element
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=5
		ielem_nodes(1:3,kk)=idv(1:3)
		ielem_nodes(4:6,kk)=idv(5:7)
		end if
		  
		
			
	
	
	else
	read(8)itx,idv(1),idv(2),idv(3),idv(4),idv(5),idv(6),idv(7),idv(8)


	end if 
	end do
	
	
	
	
	
	
	
	
	
	
	end if







	
	

	
	
	
	
 	else
 	
 	
 	
 	
	allocate(dinoder(imaxn))
	allocate(dinoder2(imaxn))
	dinoder2(:)%numberofneib=0
	dinoder(:)%itor=0
	
	if (binio.eq.0)then
	
	
	
	do j=1,imaxe
	
	if (xmpie(j).eq.n)then

	kk=kk+1
	
	

	read(8,*) itx,idv(1),idv(2),idv(3),idv(4)
	do kxk2=1,4
	  dinoder(idv(kxk2))%itor=idv(kxk2)
	  
	  end do
	ielem_ihexgl(kk)=itx
		if ((idv(3).ne.idv(4)))then
		shap=5	;nodal=4	!quadrilateral
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=4
		ielem_nodes(1:nodal,kk)=idv(1:nodal)

		else
		
		shap=6	;nodal=3!triangular
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=3
		ielem_nodes(1:nodal,kk)=idv(1:nodal)
		
		end if
		


	 

	  else
	read(8,*)



	end if 
	end do
	else
	
	do j=1,imaxe
	
	if (xmpie(j).eq.n)then

	kk=kk+1
	
	

	read(8) itx,idv(1),idv(2),idv(3),idv(4)
	do kxk2=1,4
	  dinoder(idv(kxk2))%itor=idv(kxk2)
	  
	  end do
	ielem_ihexgl(kk)=itx
		if ((idv(3).ne.idv(4)))then
		shap=5	;nodal=4	!quadrilateral
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=4
		ielem_nodes(1:nodal,kk)=idv(1:nodal)

		else

		shap=6	;nodal=3!triangular
		ielem_ishape(kk)=shap
		ielem_nonodes(kk)=nodal
		ielem_ifca(kk)=3
		ielem_nodes(1:nodal,kk)=idv(1:nodal)

		end if
		


	 

	  else
	read(8)itx,idv(1),idv(2),idv(3),idv(4)



	end if 
	end do
	
	
	
	
	end if
	end if

	
	
	

	
	  
      


	      if (dimensiona.eq.3)then
	      
	      
	      if (binio.eq.0)then
	      do j=1,imaxn

	      if (dinoder(j)%itor.gt.0)then


		read(9,*)inx,x,y,z
		x=x/scaler; y=y/scaler;  z=z/scaler
		xmin(n)=min(xmin(n),x)
		ymin(n)=min(ymin(n),y)
		zmin(n)=min(zmin(n),z)
		xmax(n)=max(xmax(n),x)
		ymax(n)=max(ymax(n),y)
		zmax(n)=max(zmax(n),z)
		
		
		  allocate(dinoder(j)%cord(1:3))
		  dinoder(j)%cord(1)=x
		  dinoder(j)%cord(2)=y
		  dinoder(j)%cord(3)=z
		  
		
		
		
		else
		read(9,*)
		end if
		
		
		!print_out = kmaxe/5


		end do
		else
		 do j=1,imaxn

	      if (dinoder(j)%itor.gt.0)then


		read(9)inx,x,y,z
		x=x/scaler; y=y/scaler;  z=z/scaler
		xmin(n)=min(xmin(n),x)
		ymin(n)=min(ymin(n),y)
		zmin(n)=min(zmin(n),z)
		xmax(n)=max(xmax(n),x)
		ymax(n)=max(ymax(n),y)
		zmax(n)=max(zmax(n),z)
		
		
		  allocate(dinoder(j)%cord(1:3))
		  dinoder(j)%cord(1)=x
		  dinoder(j)%cord(2)=y
		  dinoder(j)%cord(3)=z
		  
		
		
		
		else
		read(9)inx,x,y,z
		end if
		
		
		!print_out = kmaxe/5


		end do
		
		
		
		
		
		
		
		end if
		else
		
		if (binio.eq.0)then
		 do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		  

		    
! -------------------for debugging only -----------------------------------------!
! 		read(9,'(i14,1x,3es16.9)')inx,x,y,z
! -------------------for debugging only -----------------------------------------!
		read(9,*)inx,x,y
		x=x/scaler; y=y/scaler
		xmin(n)=min(xmin(n),x)
		ymin(n)=min(ymin(n),y)
		
		xmax(n)=max(xmax(n),x)
		ymax(n)=max(ymax(n),y)
		
		
		
		  allocate(dinoder(j)%cord(1:2))
		  dinoder(j)%cord(1)=x
		  dinoder(j)%cord(2)=y
		  else

		  read(9,*)
		  end if
		
		
		
		
		
		
		!print_out = kmaxe/5


		end do
		else
		 do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		  

		    
! -------------------for debugging only -----------------------------------------!
! 		read(9,'(i14,1x,3es16.9)')inx,x,y,z
! -------------------for debugging only -----------------------------------------!
		read(9)inx,x,y
		x=x/scaler; y=y/scaler
		xmin(n)=min(xmin(n),x)
		ymin(n)=min(ymin(n),y)
		
		xmax(n)=max(xmax(n),x)
		ymax(n)=max(ymax(n),y)
		
		
		
		  allocate(dinoder(j)%cord(1:2))
		  dinoder(j)%cord(1)=x
		  dinoder(j)%cord(2)=y
		  else

		  read(9)inx,x,y
		  end if
		
		
		
		
		
		
		!print_out = kmaxe/5


		end do
		
		
		
		
		


		end if
		end if
	
	call mpi_barrier(mpi_comm_world,ierror)
	 call close_input(n,itt)






	x=xmax(n);call mpi_allreduce(x,xmax,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
	x=ymax(n);call mpi_allreduce(x,ymax,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
	x=zmax(n);call mpi_allreduce(x,zmax,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
	x=xmin(n);call mpi_allreduce(x,xmin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
	x=ymin(n);call mpi_allreduce(x,ymin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
	x=zmin(n);call mpi_allreduce(x,zmin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
		
		
	  if (dimensiona.eq.3)then
		
	  
	  
	  
	  call open_input(n,itt)
	  do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		  allocate(dinoder2(j)%xne(1));dinoder2(j)%xne(1)=0
		  end if
		  
	   end do
	  
	  
	  do j=1,imaxe
	
	
	  
	  read(8) itx,idv(1),idv(2),idv(3),idv(4),idv(5),idv(6),idv(7),idv(8)

	  
	  
	    
		  if ((xmpie(j).ne.n))then
		      do kxk2=1,8
			  if (dinoder(idv(kxk2))%itor.gt.0)then
! 			  xsize(xmpie(j))=1
			  dinoder2(idv(kxk2))%xne(1)=dinoder2(idv(kxk2))%xne(1)+1
			  end if
		      end do
		  end if
	 end do
	  
	  call close_input(n,itt)
	  
	  
	  
	  
	  
	  
	  
	  
	  call open_input(n,itt)
	  do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		  if (dinoder2(j)%xne(1).gt.0)then
		  allocate(dinoder2(j)%xneib(dinoder2(j)%xne(1)));dinoder2(j)%xneib=0
		  dinoder2(j)%xne(1)=0
		  end if
		  end if
	   end do
	  
	  
	  do j=1,imaxe
	
	
	  
	  read(8) itx,idv(1),idv(2),idv(3),idv(4),idv(5),idv(6),idv(7),idv(8)

	  
	  
	    
		  if ((xmpie(j).ne.n))then
		      do kxk2=1,8
			  if (dinoder(idv(kxk2))%itor.gt.0)then
! 			  xsize(xmpie(j))=1
			  dinoder2(idv(kxk2))%xne(1)=dinoder2(idv(kxk2))%xne(1)+1
			  dinoder2(idv(kxk2))%xneib(dinoder2(idv(kxk2))%xne(1))=j
			  end if
		      end do
		  end if
	 end do
	  
	  call close_input(n,itt)
	  
	  
	  
	  end if
	  
	  
	  
	  
	  
	  

! dummyout=timec5
! dummyin=0.0
! call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
! !call mpi_barrier(mpi_comm_world,ierror)
! timec5=dummyin





xper = (xmax(n))- xmin(n)
yper = (ymax(n)) - ymin(n)
zper = (zmax(n)) - zmin(n)
! 
! 
! xper = abs(xmax(n))! - xmin(n)
! yper = abs(ymax(n))! - ymin(n)
! zper = abs(zmax(n))! - zmin(n)



	if (n.eq.0) then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	 write(63,*)xper,yper,zper,"periodics"
     	write(63,*)"min",xmin(n),ymin(n),zmin(n)
    	write(63,*)"max",xmax(n),ymax(n),zmax(n)
	close(63)
	end if
	

	  
	end subroutine read_input





subroutine read_input_period(n,xmpielrank,xmpinrank,xmpie,xmpin,dinode,imaxn,imaxe,imaxb,xmpinnumber,scaler)
implicit none
!> @brief
!> this subroutine reads all the periodic data from the grid files
	type(anode_number),allocatable,dimension(:,:),intent(inout)::dinode
	integer,allocatable,dimension(:,:,:),intent(inout)::xmpinnumber
	integer,intent(in)::n,imaxn,imaxe
	integer,allocatable,dimension(:),intent(in)::xmpie,xmpin
	integer,allocatable,dimension(:),intent(in)::xmpielrank,xmpinrank
	integer,intent(inout)::imaxb
	real,intent(inout)::scaler
	integer,allocatable,dimension(:)::nodep,nodec
	integer::i,j,ji,k,lm,iex,kmaxn,kk,kmaxe,print_out,kk2,shap,nodal,iftrue,kxk2
	integer,dimension(8)::idv
	integer::it1,it2,it3,it4,it5,it6,it7,it8,itx,inx,it55,in,out
	real::x,y,z


call open_input(n,itt)

if (binio.eq.0)then
if (dimensiona.eq.3)then
do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		read(9,*)inx,x,y,z
		x=x/scaler; y=y/scaler ; z=z/scaler
		      if (dinoder2(j)%numberofneib.eq.0)then
		      
		      
			allocate(dinoder(j)%cord(1:3))
! 			  if ((j.eq.709).or.(j.eq.710).or.(j.eq.693).or.(j.eq.692))then

! 			  end if
			dinoder(j)%cord(1)=x
			dinoder(j)%cord(2)=y
			dinoder(j)%cord(3)=z
		      end if
		  
		  else
		  
		  read(9,*)
		  end if
end do

else
do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		read(9,*)inx,x,y
		x=x/scaler; y=y/scaler
		      if (dinoder2(j)%numberofneib.eq.0)then
		      
		      
			allocate(dinoder(j)%cord(1:2))
			dinoder(j)%cord(1)=x
			dinoder(j)%cord(2)=y
		      end if
		  else
		  read(9,*)
		  end if
end do





end if
else
if (dimensiona.eq.3)then
do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		read(9)inx,x,y,z
		x=x/scaler; y=y/scaler ; z=z/scaler
		      if (dinoder2(j)%numberofneib.eq.0)then
		      
		      
			allocate(dinoder(j)%cord(1:3))
! 			  if ((j.eq.709).or.(j.eq.710).or.(j.eq.693).or.(j.eq.692))then

! 			  end if
			dinoder(j)%cord(1)=x
			dinoder(j)%cord(2)=y
			dinoder(j)%cord(3)=z
		      end if
		  
		  else
		  
		  read(9)inx,x,y,z
		  end if
end do

else
do j=1,imaxn


		  if (dinoder(j)%itor.gt.0)then
		read(9)inx,x,y
		x=x/scaler; y=y/scaler
		      if (dinoder2(j)%numberofneib.eq.0)then
		      
		      
			allocate(dinoder(j)%cord(1:2))
			dinoder(j)%cord(1)=x
			dinoder(j)%cord(2)=y
		      end if
		  else
		  read(9)inx,x,y
		  end if
end do





end if









end if



call mpi_barrier(mpi_comm_world,ierror)
	 call close_input(n,itt)

! if (dimensiona.eq.3)then
!  deallocate(inoder2)

! end if




end subroutine read_input_period



subroutine stenprint(n)
implicit none
!> @brief
!> this subroutine prints the stencils at a selected position
integer,intent(in)::n
integer::i,j,k,inv,ismp,l,i1,i2,i3,i4,i5,i6,i7,i8,ixx
integer::kmaxe,kk,kfk,icpuid,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,decomf,kd,itarget
integer::inx,m,o,p,q,jk,elementss
logical::herev
real,dimension(5)::total
character(len=20)::proc,outfile,proc3,surfile
 





icpuid=n
	
	



















if (nprobes.gt.0)then
do inv=1,nprobes
				  
				if (probei(n,inv).ne.0) then
			write(proc3,fmt='(i10)') inv
	!proc4=".plt"
	outfile="stencils_"//trim(adjustl(proc3))//".dat"!//trim(adjustl(proc4))	
				
				
	open(97,file=outfile,form='formatted',status='new',action='write')
	if (binio.eq.0)open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	if (binio.eq.1)open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	
	write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	
	if (dimensiona.eq.3)then
	write(97,*)'VARIABLES="x","y","z"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	else
	write(97,*)'VARIABLES="x","y"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	
	end if
	
	if (binio.eq.0)then
	do i=1,imaxn
		read(96,*)j,x
		write(97,*)x/scaler
	end do

	close(96)
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	do i=1,imaxn
		read(96,*)j,x,y
		write(97,*)y/scaler
	end do
	close(96)
	
	if (dimensiona.eq.3)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	do i=1,imaxn
		read(96,*)j,x,y,z
		write(97,*)z/scaler
	end do
	close (97)
	close(96)
        end if
        
        
	else
	do i=1,imaxn
		read(96)j,x
		write(97,*)x/scaler
	end do

	close(96)
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	do i=1,imaxn
		read(96)j,x,y
		write(97,*)y/scaler
	end do
	close(96)
	
	if (dimensiona.eq.3)then
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	do i=1,imaxn
		read(96)j,x,y,z
		write(97,*)z/scaler
	end do
	close (97)
	close(96)
        end if
	
	
	end if
	
	
	if (binio.eq.0)then
	
	open(96,file='GRID.cel',form='formatted',status='old',action='read')
	open(97,file=outfile,form='formatted',status='old',action='write',position='append')
	if (dimensiona.eq.3)then
	do i=1,imaxe
	      read(96,*)ix,i5,i6,i8,i7,i1,i2,i4,i3
	      write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
	end do
	else
	do i=1,imaxe
	      read(96,*)ix,i1,i2,i3,i4
	      write(97,*)i1,i2,i3,i4
	end do
	end if
	close(96)
	close(97)
	else
	open(96,file='GRID.cel',form='unformatted',status='old',action='read')
	open(97,file=outfile,form='formatted',status='old',action='write',position='append')
	if (dimensiona.eq.3)then
	do i=1,imaxe
	      read(96)ix,i5,i6,i8,i7,i1,i2,i4,i3
	      write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
	end do
	else
	do i=1,imaxe
	      read(96)ix,i1,i2,i3,i4
	      write(97,*)i1,i2,i3,i4
	end do
	end if
	close(96)
	close(97)
	
	end if
				
				
				
				
				
				
				
				
				
				
                                    open(97,file=outfile,form='formatted',status='old',action='write',position='append')
                                    
                                    if (binio.eq.0)open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                    if (binio.eq.1)open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                    write(97,*)'TITLE="GRID"'
                                    write(97,*)'filetype=grid'
                                    if (dimensiona.eq.3)then
                                    write(97,*)'VARIABLES="x","y","z"'
                                    write(97,*)'Zone N=',imaxn,',E=',1,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
                                    else
                                    write(97,*)'VARIABLES="x","y"'
                                    write(97,*)'Zone N=',imaxn,',E=',1,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
                                    
                                    end if
                                    
					    if (binio.eq.0)then
                                            do i=1,imaxn
                                                    read(96,*)j,x
                                                    write(97,*)x/scaler
                                            end do

                                            close(96)
                                            open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96,*)j,x,y
                                                    write(97,*)y/scaler
                                            end do
                                            close(96)
                                            
                                            if (dimensiona.eq.3)then
                                            open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96,*)j,x,y,z
                                                    write(97,*)z/scaler
                                            end do
                                            close (97)
                                            close(96)
                                            end if
                                            else
                                            do i=1,imaxn
                                                    read(96)j,x
                                                    write(97,*)x/scaler
                                            end do

                                            close(96)
                                            open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96)j,x,y
                                                    write(97,*)y/scaler
                                            end do
                                            close(96)
                                            
                                            if (dimensiona.eq.3)then
                                            open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96)j,x,y,z
                                                    write(97,*)z/scaler
                                            end do
                                            close (97)
                                            close(96)
                                            end if
                                            
                                            end if
                                            
                    
											



                                            if (binio.eq.0)then
                                            open(96,file='GRID.cel',form='formatted',status='old',action='read')
                                            open(97,file=outfile,form='formatted',status='old',action='write',position='append')
                                            if (dimensiona.eq.3)then
                                            do i=1,imaxe
                                                read(96,*)ix,i5,i6,i8,i7,i1,i2,i4,i3
                                                if (ielem_ihexgl(probei(n,inv)).eq.i)then
                                                write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
                                                end if
                                            end do
                                            else
                                            do i=1,imaxe
                                                read(96,*)ix,i1,i2,i3,i4
                                                 if (ielem_ihexgl(probei(n,inv)).eq.i)then
                                                write(97,*)i1,i2,i3,i4
                                                end if
                                            end do
                                            end if
                                            close(96)
                                            close(97)
					     else
					      open(96,file='GRID.cel',form='unformatted',status='old',action='read')
                                            open(97,file=outfile,form='formatted',status='old',action='write',position='append')
                                            if (dimensiona.eq.3)then
                                            do i=1,imaxe
                                                read(96)ix,i5,i6,i8,i7,i1,i2,i4,i3
                                                if (ielem_ihexgl(probei(n,inv)).eq.i)then
                                                write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
                                                end if
                                            end do
                                            else
                                            do i=1,imaxe
                                                read(96)ix,i1,i2,i3,i4
                                                 if (ielem_ihexgl(probei(n,inv)).eq.i)then
                                                write(97,*)i1,i2,i3,i4
                                                end if
                                            end do
                                            end if
                                            close(96)
                                            close(97)
					     
					     
					     end if
                                    
                                    
                                    
                                    
                                    
				
				      do ismp=1,typesten
					elementss=0
					
					if ((ismp.eq.1).or.(ees.ne.5))then
					itarget=ielem_inumneighbours(probei(n,inv))
					else
					itarget=numneighbours2
					end if
        
					
					
					
					do  l=2,itarget
                                            if (ilocalstencil(n,probei(n,inv),ismp,l).gt.0)then
                                                     elementss= elementss+1 
                                            end if
                                        end do
                                        
                                            if (elementss+1.eq.itarget)then
                                                                                     
                                            open(97,file=outfile,form='formatted',status='old',action='write',position='append')
                                    if (binio.eq.0)open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                    if (binio.eq.1)open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                    write(97,*)'TITLE="GRID"'
                                    write(97,*)'filetype=grid'
                                    if (dimensiona.eq.3)then
                                    write(97,*)'VARIABLES="x","y","z"'
                                    write(97,*)'Zone N=',imaxn,',E=',elementss,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
                                    else
                                    write(97,*)'VARIABLES="x","y"'
                                    write(97,*)'Zone N=',imaxn,',E=',elementss,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
                                    
                                    end if
                                    
					    if (binio.eq.0)then
                                            do i=1,imaxn
                                                    read(96,*)j,x
                                                    write(97,*)x/scaler
                                            end do

                                            close(96)
                                            open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96,*)j,x,y
                                                    write(97,*)y/scaler
                                            end do
                                            close(96)
                                            
                                            if (dimensiona.eq.3)then
                                            open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96,*)j,x,y,z
                                                    write(97,*)z/scaler
                                            end do
                                            close (97)
                                            close(96)
                                            end if
                                            else
                                            do i=1,imaxn
                                                    read(96)j,x
                                                    write(97,*)x/scaler
                                            end do

                                            close(96)
                                            open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96)j,x,y
                                                    write(97,*)y/scaler
                                            end do
                                            close(96)
                                            
                                            if (dimensiona.eq.3)then
                                            open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                            do i=1,imaxn
                                                    read(96)j,x,y,z
                                                    write(97,*)z/scaler
                                            end do
                                            close (97)
                                            close(96)
                                            end if
                                            
                                            
                                            end if
                                            
                                                                                        
                                            do  l=2,itarget
                                            	if (ilocalstencil(n,probei(n,inv),ismp,l).gt.0)then
                                            	if (binio.eq.0)then
                                            	 open(96,file='GRID.cel',form='formatted',status='old',action='read')
                                            open(97,file=outfile,form='formatted',status='old',action='write',position='append')
                                            
                                                    if (dimensiona.eq.3)then
                                                        do i=1,imaxe
                                                            read(96,*)ix,i5,i6,i8,i7,i1,i2,i4,i3
                                                            if (i.eq.ilocalstencil(n,probei(n,inv),ismp,l))then
                                                            write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
                                                            end if
                                                        end do
                                                    else
                                                        do i=1,imaxe
                                                            read(96,*)ix,i1,i2,i3,i4
                                                            if (i.eq.ilocalstencil(n,probei(n,inv),ismp,l))then
                                                            write(97,*)i1,i2,i3,i4
                                                            end if
                                                        end do
                                                    end if
                                            close(96)
                                            close(97)
						else
						 open(96,file='GRID.cel',form='unformatted',status='old',action='read')
                                            open(97,file=outfile,form='formatted',status='old',action='write',position='append')
                                            
                                                    if (dimensiona.eq.3)then
                                                        do i=1,imaxe
                                                            read(96)ix,i5,i6,i8,i7,i1,i2,i4,i3
                                                            if (i.eq.ilocalstencil(n,probei(n,inv),ismp,l))then
                                                            write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
                                                            end if
                                                        end do
                                                    else
                                                        do i=1,imaxe
                                                            read(96)ix,i1,i2,i3,i4
                                                            if (i.eq.ilocalstencil(n,probei(n,inv),ismp,l))then
                                                            write(97,*)i1,i2,i3,i4
                                                            end if
                                                        end do
                                                    end if
                                            close(96)
                                            close(97)
						
						
						end if
                                                end if
                                            end do
					   
					  end if
					end do
				      
				   end if
				  
				  
				  
			    end do
 			  

end if

 call mpi_barrier(mpi_comm_world,ierror)
end subroutine stenprint









subroutine walldistance(n,imaxe,xmpielrank)
!> @brief
!> this subroutine establishes the wall distance for every cell in 3d

implicit none

integer,allocatable,dimension(:),intent(in)::xmpielrank !,xmpinrank
integer,intent(in)::n,imaxe
integer :: countwall,kmaxe,countwallglobal,i,l,icpu,doyouhavewall,howmanyhavewall,counterall,j
integer :: counterall2,wall1,wall2,wall3,wall4,wall5,wall6,ioy,wl1,wl2,wl3,wl4,k
real,allocatable,dimension(:,:) :: wallelemarraycord,wallelemarraycordglobal
real :: distance
character(len=12)::bndfile,vrtfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type :: awallboundary
	    integer :: globalid,wbid,wb1,wb2,wb3,wb4,numnodes!,wbdescr
	    real :: wallx,wally,wallz,walltrans
end type
type(awallboundary),allocatable,dimension(:):: dwallbnd
type :: anodeswall
	    integer :: id
	    real :: wnx,wny,wnz
end type
type(anodeswall),allocatable,dimension(:):: dwallvrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! how many element for this block have wall
kmaxe=xmpielrank(n)
 countwall = 0
 do i=1,kmaxe
    if (ielem_interior(i).eq.1)then
    do l=1,ielem_ifca(i)
	if (ielem_ibounds(l,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(l,i)).eq.4).or.(ibound_icode(ielem_ibounds(l,i)).eq.99))then
	  countwall = countwall + 1
	      end if
	      
	end if
    end do
    end if
end do
! how many for all cpu blocks
call mpi_allreduce(countwall,countwallglobal,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)
!  print*,countwallglobal,'total wall elements',n
allocate(dwallbnd(countwallglobal))

 countwall = 0
bndfile='GRID.bnd'
vrtfile='GRID.vrt'
if (binio.eq.0)open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
if (binio.eq.1)open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)
! find the nodes id of the wall globally from the bnd file
 countwall = 0
if (binio.eq.0)then
do i=1,imaxb
	read(10,*)wall1,wall2,wall3,wall4,wall5,wall6
	    if (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		dwallbnd(countwall)%wbid = wall1 ; dwallbnd(countwall)%wb1 = wall2 ; dwallbnd(countwall)%wb2 = wall3
		dwallbnd(countwall)%wb3 = wall4  ;dwallbnd(countwall)%wb4 = wall5
	    end if
end do
else
do i=1,imaxb
	read(10)wall1,wall2,wall3,wall4,wall5,wall6
	    if (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		dwallbnd(countwall)%wbid = wall1 ; dwallbnd(countwall)%wb1 = wall2 ; dwallbnd(countwall)%wb2 = wall3
		dwallbnd(countwall)%wb3 = wall4  ;dwallbnd(countwall)%wb4 = wall5
	    end if
end do

end if
 close(10)

!  print*,countwall,'total wall elements2',n
allocate(dwallvrt(imaxn))
if (binio.eq.0)open(11,file=vrtfile,form='formatted',status='old',action='read')
if (binio.eq.1)open(11,file=vrtfile,form='unformatted',status='old',action='read')
! read node coordinates and store
if (binio.eq.0)then
do i = 1,imaxn
  read(11,*)dwallvrt(i)%id,dwallvrt(i)%wnx,dwallvrt(i)%wny,dwallvrt(i)%wnz
  dwallvrt(i)%wnx=dwallvrt(i)%wnx/scaler
  dwallvrt(i)%wny=dwallvrt(i)%wny/scaler
  dwallvrt(i)%wnz=dwallvrt(i)%wnz/scaler
end do
else
do i = 1,imaxn
  read(11)dwallvrt(i)%id,dwallvrt(i)%wnx,dwallvrt(i)%wny,dwallvrt(i)%wnz
  dwallvrt(i)%wnx=dwallvrt(i)%wnx/scaler
  dwallvrt(i)%wny=dwallvrt(i)%wny/scaler
  dwallvrt(i)%wnz=dwallvrt(i)%wnz/scaler
end do


end if
  close(11)
! compute and store wall face centers globally
	do i = 1,countwallglobal

		wl1 = dwallbnd(i)%wb1 ;wl2 = dwallbnd(i)%wb2 ;wl3 = dwallbnd(i)%wb3 ;wl4 = dwallbnd(i)%wb4 ;
	if (wl4.eq.wl3)then
	dwallbnd(i)%wallx=(dwallvrt(wl1)%wnx+dwallvrt(wl2)%wnx+dwallvrt(wl3)%wnx)/3.0
	dwallbnd(i)%wally=(dwallvrt(wl1)%wny+dwallvrt(wl2)%wny+dwallvrt(wl3)%wny)/3.0
	dwallbnd(i)%wallz=(dwallvrt(wl1)%wnz+dwallvrt(wl2)%wnz+dwallvrt(wl3)%wnz)/3.0


	else
	dwallbnd(i)%wallx = (dwallvrt(wl1)%wnx + dwallvrt(wl2)%wnx + dwallvrt(wl3)%wnx + dwallvrt(wl4)%wnx)/4.0
	dwallbnd(i)%wally = (dwallvrt(wl1)%wny + dwallvrt(wl2)%wny + dwallvrt(wl3)%wny + dwallvrt(wl4)%wny)/4.0
		dwallbnd(i)%wallz = (dwallvrt(wl1)%wnz + dwallvrt(wl2)%wnz + dwallvrt(wl3)%wnz + dwallvrt(wl4)%wnz)/4.0
		end if
		if (transition_model.eq.1)then
		    if (transition_axis.eq.2)then
			dwallbnd(i)%walltrans=dwallbnd(i)%wally-transition_location
		    else if (transition_axis.eq.3)then
			dwallbnd(i)%walltrans=dwallbnd(i)%wallz-transition_location
		    else
			dwallbnd(i)%walltrans=dwallbnd(i)%wallx-transition_location
		    end if
		else
		    dwallbnd(i)%walltrans=0.0
		end if
	end do
! find distance from element barycenter to the nearest wall for this block.
kmaxe=xmpielrank(n)
	do i=1,kmaxe
	    distance=tolbig 
	    ielem_walltrans(i)=0.0
		do k = 1,countwallglobal
		      if ( distance .gt. (sqrt(((dwallbnd(k)%wallx-ielem_xxc(i))**2) &
				     + ((dwallbnd(k)%wally-ielem_yyc(i))**2)&
				     + ((dwallbnd(k)%wallz-ielem_zzc(i))**2)))) then
		distance = sqrt(((dwallbnd(k)%wallx-ielem_xxc(i))**2)&
			       +((dwallbnd(k)%wally-ielem_yyc(i))**2)&
			       +((dwallbnd(k)%wallz-ielem_zzc(i))**2))
			ielem_walldist(i) = distance
			ielem_walltrans(i) = dwallbnd(k)%walltrans
			 if (ielem_walldist(i).lt.hybridist)then
		    ielem_hybrid(i)=1
		 end if
		
	      end if
	end do
end do

deallocate(dwallbnd)
deallocate(dwallvrt)

end subroutine


subroutine walldistancex(n, imaxe, xmpielrank)
!> establish wall distance for every cell in 3D
!> rewritten with flat arrays + spatial bins + OpenMP

implicit none

integer, allocatable, dimension(:), intent(in) :: xmpielrank
integer, intent(in) :: n, imaxe

!-----------------------------
! basic integers
!-----------------------------
integer :: i, j, k, l, p
integer :: kmaxe, ioy, ibin
integer :: countwall, countwallglobal
integer :: wall1, wall2, wall3, wall4, wall5, wall6
integer :: wl1, wl2, wl3, wl4

character(len=12) :: bndfile, vrtfile

!-----------------------------
! flat arrays for wall faces
!-----------------------------
integer, allocatable :: wbid(:), wb1(:), wb2(:), wb3(:), wb4(:)

!-----------------------------
! flat arrays for wall centers
!-----------------------------
real, allocatable :: wallx(:), wally(:), wallz(:), walltr(:)

!-----------------------------
! flat arrays for nodes
!-----------------------------
integer, allocatable :: vid(:)
real, allocatable :: vx(:), vy(:), vz(:)

!-----------------------------
! bin data
!-----------------------------
integer, allocatable :: bincount(:), binstart(:), binfill(:), binlist(:)
integer :: nx, ny, nz, nbins
integer :: ix, iy, iz
integer :: target_walls_per_bin
integer :: shellmax
integer :: ic, jc, kc, shell
integer :: ilo, ihi, jlo, jhi, klo, khi
integer :: ii, jj, kk, b
integer :: p1, p2, bestwall

!-----------------------------
! reals
!-----------------------------
real :: xc, yc, zc
real :: dx, dy, dz, d2, mind2
real :: xminw, xmaxw, yminw, ymaxw, zminw, zmaxw
real :: xrange, yrange, zrange, volbox
real :: hbin, dxbin, dybin, dzbin
real :: x0, x1, y0, y1, z0, z1
real :: dbox2, lower2
real :: leftd, rightd, backd, frontd, botd, topd
real :: hybridist2
real, parameter :: tinybin = 1.0e-20

!------------------------------------------------------------
! count local wall faces from local mesh info
!------------------------------------------------------------
kmaxe = xmpielrank(n)

countwall = 0
do i = 1, kmaxe
    if (ielem_interior(i) .eq. 1) then
        do l = 1, ielem_ifca(i)
            if (ielem_ibounds(l,i) .gt. 0) then
                if (ibound_icode(ielem_ibounds(l,i)) .eq. 4) then
                    countwall = countwall + 1
                end if
            end if
        end do
    end if
end do

call mpi_allreduce(countwall, countwallglobal, 1, mpi_integer, mpi_sum, mpi_comm_world, ierror)
call mpi_barrier(mpi_comm_world, ierror)

!------------------------------------------------------------
! allocate flat arrays
!------------------------------------------------------------
allocate(wbid(countwallglobal))
allocate(wb1(countwallglobal))
allocate(wb2(countwallglobal))
allocate(wb3(countwallglobal))
allocate(wb4(countwallglobal))

allocate(wallx(countwallglobal))
allocate(wally(countwallglobal))
allocate(wallz(countwallglobal))
allocate(walltr(countwallglobal))

allocate(vid(imaxn))
allocate(vx(imaxn))
allocate(vy(imaxn))
allocate(vz(imaxn))

bndfile = 'GRID.bnd'
vrtfile = 'GRID.vrt'

!------------------------------------------------------------
! read boundary file
!------------------------------------------------------------
if (binio .eq. 0) open(10, file=bndfile, form='formatted',   status='old', action='read', iostat=ioy)
if (binio .eq. 1) open(10, file=bndfile, form='unformatted', status='old', action='read', iostat=ioy)

if (ioy /= 0) then
    print *, 'Error opening ', trim(bndfile), ' on rank ', n
    call mpi_abort(mpi_comm_world, 1, ierror)
end if

countwall = 0

if (binio .eq. 0) then
    do i = 1, imaxb
        read(10,*) wall1, wall2, wall3, wall4, wall5, wall6
        if (wall6 .eq. 4) then
            countwall       = countwall + 1
            wbid(countwall) = wall1
            wb1(countwall)  = wall2
            wb2(countwall)  = wall3
            wb3(countwall)  = wall4
            wb4(countwall)  = wall5
        end if
    end do
else
    do i = 1, imaxb
        read(10) wall1, wall2, wall3, wall4, wall5, wall6
        if (wall6 .eq. 4) then
            countwall       = countwall + 1
            wbid(countwall) = wall1
            wb1(countwall)  = wall2
            wb2(countwall)  = wall3
            wb3(countwall)  = wall4
            wb4(countwall)  = wall5
        end if
    end do
end if

close(10)

if (countwall /= countwallglobal) then
    print *, 'Mismatch in wall counts on rank ', n, ': file=', countwall, ' global=', countwallglobal
    call mpi_abort(mpi_comm_world, 2, ierror)
end if

!------------------------------------------------------------
! read vertex file
!------------------------------------------------------------
if (binio .eq. 0) open(11, file=vrtfile, form='formatted',   status='old', action='read', iostat=ioy)
if (binio .eq. 1) open(11, file=vrtfile, form='unformatted', status='old', action='read', iostat=ioy)

if (ioy /= 0) then
    print *, 'Error opening ', trim(vrtfile), ' on rank ', n
    call mpi_abort(mpi_comm_world, 3, ierror)
end if

if (binio .eq. 0) then
    do i = 1, imaxn
        read(11,*) vid(i), vx(i), vy(i), vz(i)
        vx(i) = vx(i) / scaler
        vy(i) = vy(i) / scaler
        vz(i) = vz(i) / scaler
    end do
else
    do i = 1, imaxn
        read(11) vid(i), vx(i), vy(i), vz(i)
        vx(i) = vx(i) / scaler
        vy(i) = vy(i) / scaler
        vz(i) = vz(i) / scaler
    end do
end if

close(11)

!------------------------------------------------------------
! compute wall-face centers
!------------------------------------------------------------
do i = 1, countwallglobal

    wl1 = wb1(i)
    wl2 = wb2(i)
    wl3 = wb3(i)
    wl4 = wb4(i)

    if (wl1 < 1 .or. wl1 > imaxn .or. &
        wl2 < 1 .or. wl2 > imaxn .or. &
        wl3 < 1 .or. wl3 > imaxn .or. &
        wl4 < 1 .or. wl4 > imaxn) then
        print *, 'Wall node index out of range on rank ', n, ' at face ', i
        print *, 'wl1,wl2,wl3,wl4 = ', wl1, wl2, wl3, wl4
        call mpi_abort(mpi_comm_world, 4, ierror)
    end if

    if (wl4 .eq. wl3) then
        wallx(i) = (vx(wl1) + vx(wl2) + vx(wl3)) / 3.0
        wally(i) = (vy(wl1) + vy(wl2) + vy(wl3)) / 3.0
        wallz(i) = (vz(wl1) + vz(wl2) + vz(wl3)) / 3.0
    else
        wallx(i) = (vx(wl1) + vx(wl2) + vx(wl3) + vx(wl4)) / 4.0
        wally(i) = (vy(wl1) + vy(wl2) + vy(wl3) + vy(wl4)) / 4.0
        wallz(i) = (vz(wl1) + vz(wl2) + vz(wl3) + vz(wl4)) / 4.0
    end if

	end do

	if (transition_model .eq. 1) then
	    if (transition_axis .eq. 2) then
	        walltr(1:countwallglobal) = wally(1:countwallglobal) - transition_location
	    else if (transition_axis .eq. 3) then
	        walltr(1:countwallglobal) = wallz(1:countwallglobal) - transition_location
	    else
	        walltr(1:countwallglobal) = wallx(1:countwallglobal) - transition_location
	    end if
	else
	    walltr(1:countwallglobal) = 0.0
	end if

	!------------------------------------------------------------
! build 3D spatial bins for wall centers
!------------------------------------------------------------
xminw = minval(wallx)
xmaxw = maxval(wallx)
yminw = minval(wally)
ymaxw = maxval(wally)
zminw = minval(wallz)
zmaxw = maxval(wallz)

xrange = max(xmaxw - xminw, tinybin)
yrange = max(ymaxw - yminw, tinybin)
zrange = max(zmaxw - zminw, tinybin)

! tune this if needed: 16, 32, 64, 128
target_walls_per_bin = 64

volbox = xrange * yrange * zrange
hbin   = (volbox * real(target_walls_per_bin) / real(countwallglobal)) ** (1.0 / 3.0)
hbin   = max(hbin, tinybin)

dxbin = hbin
dybin = hbin
dzbin = hbin

nx = int(xrange / dxbin) + 1
ny = int(yrange / dybin) + 1
nz = int(zrange / dzbin) + 1
nbins = nx * ny * nz

allocate(bincount(nbins))
allocate(binstart(nbins+1))
allocate(binfill(nbins))
allocate(binlist(countwallglobal))

bincount = 0

!------------------------------------------------------------
! first pass: count walls in each bin
!------------------------------------------------------------
do k = 1, countwallglobal

    ix = int((wallx(k) - xminw) / dxbin) + 1
    iy = int((wally(k) - yminw) / dybin) + 1
    iz = int((wallz(k) - zminw) / dzbin) + 1

    if (ix < 1) ix = 1
    if (ix > nx) ix = nx
    if (iy < 1) iy = 1
    if (iy > ny) iy = ny
    if (iz < 1) iz = 1
    if (iz > nz) iz = nz

    ibin = ix + (iy-1)*nx + (iz-1)*nx*ny
    bincount(ibin) = bincount(ibin) + 1
end do

!------------------------------------------------------------
! prefix sum
!------------------------------------------------------------
binstart(1) = 1
do b = 1, nbins
    binstart(b+1) = binstart(b) + bincount(b)
end do

binfill = binstart(1:nbins)

!------------------------------------------------------------
! second pass: fill packed wall-index list
!------------------------------------------------------------
do k = 1, countwallglobal

    ix = int((wallx(k) - xminw) / dxbin) + 1
    iy = int((wally(k) - yminw) / dybin) + 1
    iz = int((wallz(k) - zminw) / dzbin) + 1

    if (ix < 1) ix = 1
    if (ix > nx) ix = nx
    if (iy < 1) iy = 1
    if (iy > ny) iy = ny
    if (iz < 1) iz = 1
    if (iz > nz) iz = nz

    ibin = ix + (iy-1)*nx + (iz-1)*nx*ny

    p = binfill(ibin)
    binlist(p) = k
    binfill(ibin) = p + 1
end do

shellmax   = max(nx, max(ny, nz))
hybridist2 = hybridist * hybridist
kmaxe      = xmpielrank(n)

! optional diagnostics
! print *, 'rank=', n, ' countwallglobal=', countwallglobal
! print *, 'rank=', n, ' nx,ny,nz=', nx, ny, nz, ' nbins=', nbins
! print *, 'rank=', n, ' avg walls/bin=', real(countwallglobal)/real(nbins)

!------------------------------------------------------------
! nearest-wall search using expanding bin shells + OpenMP
!------------------------------------------------------------
!$omp parallel do default(none) &
	!$omp shared(kmaxe,ielem_xxc,ielem_yyc,ielem_zzc,ielem_walldist,ielem_walltrans,ielem_hybrid, &
	!$omp        wallx,wally,wallz,xminw,yminw,zminw,dxbin,dybin,dzbin, &
	!$omp        nx,ny,nz,binstart,binlist,walltr,hybridist2,tolbig,shellmax) &
	!$omp private(i,xc,yc,zc,ic,jc,kc,mind2,shell,ilo,ihi,jlo,jhi,klo,khi, &
	!$omp         ii,jj,kk,b,p1,p2,p,k,x0,x1,y0,y1,z0,z1,dbox2,dx,dy,dz,d2, &
	!$omp         leftd,rightd,backd,frontd,botd,topd,lower2,bestwall) &
!$omp schedule(static)
do i = 1, kmaxe

    xc = ielem_xxc(i)
    yc = ielem_yyc(i)
    zc = ielem_zzc(i)

    ic = int((xc - xminw) / dxbin) + 1
    jc = int((yc - yminw) / dybin) + 1
    kc = int((zc - zminw) / dzbin) + 1

    if (ic < 1) ic = 1
    if (ic > nx) ic = nx
    if (jc < 1) jc = 1
    if (jc > ny) jc = ny
    if (kc < 1) kc = 1
    if (kc > nz) kc = nz

	    mind2 = tolbig * tolbig
	    bestwall = 0
	    ielem_hybrid(i) = 0
	    ielem_walltrans(i) = 0.0

    do shell = 0, shellmax

        ilo = max(1,  ic - shell)
        ihi = min(nx, ic + shell)
        jlo = max(1,  jc - shell)
        jhi = min(ny, jc + shell)
        klo = max(1,  kc - shell)
        khi = min(nz, kc + shell)

        ! search only the outer surface of the shell
        do kk = klo, khi
            do jj = jlo, jhi
                do ii = ilo, ihi

                    if (shell > 0) then
                        if (ii > ilo .and. ii < ihi .and. &
                            jj > jlo .and. jj < jhi .and. &
                            kk > klo .and. kk < khi) cycle
                    end if

                    b = ii + (jj-1)*nx + (kk-1)*nx*ny

                    ! empty bin
                    if (binstart(b) == binstart(b+1)) cycle

                    ! point-to-bin-box lower bound
                    x0 = xminw + real(ii-1) * dxbin
                    x1 = xminw + real(ii  ) * dxbin
                    y0 = yminw + real(jj-1) * dybin
                    y1 = yminw + real(jj  ) * dybin
                    z0 = zminw + real(kk-1) * dzbin
                    z1 = zminw + real(kk  ) * dzbin

                    if (xc < x0) then
                        dx = x0 - xc
                    else if (xc > x1) then
                        dx = xc - x1
                    else
                        dx = 0.0
                    end if

                    if (yc < y0) then
                        dy = y0 - yc
                    else if (yc > y1) then
                        dy = yc - y1
                    else
                        dy = 0.0
                    end if

                    if (zc < z0) then
                        dz = z0 - zc
                    else if (zc > z1) then
                        dz = zc - z1
                    else
                        dz = 0.0
                    end if

                    dbox2 = dx*dx + dy*dy + dz*dz
                    if (dbox2 >= mind2) cycle

                    p1 = binstart(b)
                    p2 = binstart(b+1) - 1

                    do p = p1, p2
                        k = binlist(p)

                        dx = wallx(k) - xc
                        dy = wally(k) - yc
                        dz = wallz(k) - zc

                        d2 = dx*dx + dy*dy + dz*dz

	                        if (d2 < mind2) then
	                            mind2 = d2
	                            bestwall = k
	                        end if
                    end do

                end do
            end do
        end do

        ! exact stopping rule:
        ! if the nearest possible point outside the searched box
        ! is already farther than current best, stop
        lower2 = tolbig * tolbig

        if (ilo > 1) then
            leftd  = xc - (xminw + real(ilo-1) * dxbin)
            lower2 = min(lower2, leftd*leftd)
        end if

        if (ihi < nx) then
            rightd = (xminw + real(ihi) * dxbin) - xc
            lower2 = min(lower2, rightd*rightd)
        end if

        if (jlo > 1) then
            backd  = yc - (yminw + real(jlo-1) * dybin)
            lower2 = min(lower2, backd*backd)
        end if

        if (jhi < ny) then
            frontd = (yminw + real(jhi) * dybin) - yc
            lower2 = min(lower2, frontd*frontd)
        end if

        if (klo > 1) then
            botd   = zc - (zminw + real(klo-1) * dzbin)
            lower2 = min(lower2, botd*botd)
        end if

        if (khi < nz) then
            topd   = (zminw + real(khi) * dzbin) - zc
            lower2 = min(lower2, topd*topd)
        end if

        if (mind2 <= lower2) exit

    end do

	    ielem_walldist(i) = sqrt(mind2)
	    if (bestwall .gt. 0) ielem_walltrans(i) = walltr(bestwall)
	    if (mind2 < hybridist2) ielem_hybrid(i) = 1

end do
!$omp end parallel do

!------------------------------------------------------------
! cleanup
!------------------------------------------------------------
deallocate(wbid, wb1, wb2, wb3, wb4)
deallocate(wallx, wally, wallz, walltr)
deallocate(vid, vx, vy, vz)
deallocate(bincount, binstart, binfill, binlist)

end subroutine


subroutine walldistance2d(n,imaxe,xmpielrank)
!> @brief
!> this subroutine establishes the wall distance for every cell in 2d

implicit none
integer,allocatable,dimension(:),intent(in)::xmpielrank !,xmpinrank
integer,intent(in)::n,imaxe
integer :: countwall,kmaxe,countwallglobal,i,l,icpu,doyouhavewall,howmanyhavewall,counterall,j
integer :: counterall2,wall1,wall2,wall3,wall4,wall5,wall6,ioy,wl1,wl2,wl3,wl4,k
real,allocatable,dimension(:,:) :: wallelemarraycord,wallelemarraycordglobal
real :: distance
character(len=12)::bndfile,vrtfile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type :: awallboundary
	    integer :: globalid,wbid,wb1,wb2,wb3,wb4,numnodes!,wbdescr
	    real :: wallx,wally,wallz,walltrans
end type
type(awallboundary),allocatable,dimension(:):: dwallbnd
type :: anodeswall
	    integer :: id
	    real :: wnx,wny,wnz
end type
type(anodeswall),allocatable,dimension(:):: dwallvrt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! how many element for this block have wall
kmaxe=xmpielrank(n)
 countwall = 0
 do i=1,kmaxe
    if (ielem_interior(i).eq.1)then
    do l=1,ielem_ifca(i)
	if (ielem_ibounds(l,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(l,i)).eq.4).or.(ibound_icode(ielem_ibounds(l,i)).eq.99))then
	  countwall = countwall + 1
	      end if
	       
	end if
    end do
    end if
end do
! how many for all cpu blocks
call mpi_allreduce(countwall,countwallglobal,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)
call mpi_barrier(mpi_comm_world,ierror)

!  print*,countwallglobal,'total wall elements',n
allocate(dwallbnd(countwallglobal))

 countwall = 0
bndfile='GRID.bnd'
vrtfile='GRID.vrt'
if (binio.eq.0)open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
if (binio.eq.1)open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)
! find the nodes id of the wall globally from the bnd file
 countwall = 0
if (binio.eq.0)then
do i=1,imaxb
	read(10,*)wall1,wall2,wall3,wall4,wall5,wall6
	    if (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		dwallbnd(countwall)%wbid = wall1 ; dwallbnd(countwall)%wb1 = wall2 ; dwallbnd(countwall)%wb2 = wall3
		dwallbnd(countwall)%wb3 = wall4  ;dwallbnd(countwall)%wb4 = wall5
	    end if
end do
else
do i=1,imaxb
	read(10)wall1,wall2,wall3,wall4,wall5,wall6
	    if (wall6 .eq. 4) then	! wall face
		countwall = countwall + 1
		dwallbnd(countwall)%wbid = wall1 ; dwallbnd(countwall)%wb1 = wall2 ; dwallbnd(countwall)%wb2 = wall3
		dwallbnd(countwall)%wb3 = wall4  ;dwallbnd(countwall)%wb4 = wall5
	    end if
end do

end if
 close(10)

 
allocate(dwallvrt(imaxn))
if (binio.eq.0)open(11,file=vrtfile,form='formatted',status='old',action='read')
if (binio.eq.1)open(11,file=vrtfile,form='unformatted',status='old',action='read')
! read node coordinates and store
if (binio.eq.0)then
do i = 1,imaxn
  read(11,*)dwallvrt(i)%id,dwallvrt(i)%wnx,dwallvrt(i)%wny
end do
else
do i = 1,imaxn
  read(11)dwallvrt(i)%id,dwallvrt(i)%wnx,dwallvrt(i)%wny
end do

end if
  close(11)
! compute and store wall face centers globally
do i = 1,countwallglobal

	wl1 = dwallbnd(i)%wb1 ;wl2 = dwallbnd(i)%wb2 ;wl3 = dwallbnd(i)%wb3 ;wl4 = dwallbnd(i)%wb4 ;
	
	dwallbnd(i)%wallx=(dwallvrt(wl1)%wnx+dwallvrt(wl2)%wnx)/2.0
	dwallbnd(i)%wally=(dwallvrt(wl1)%wny+dwallvrt(wl2)%wny)/2.0
	if (transition_model.eq.1)then
	    if (transition_axis.eq.2)then
		dwallbnd(i)%walltrans=dwallbnd(i)%wally-transition_location
	    else
		dwallbnd(i)%walltrans=dwallbnd(i)%wallx-transition_location
	    end if
	else
	    dwallbnd(i)%walltrans=0.0
	end if
	


	
end do


! find distance from element barycenter to the nearest wall for this block.
kmaxe=xmpielrank(n)
do i=1,kmaxe
    distance=tolbig 
    ielem_walltrans(i)=0.0
	do k = 1,countwallglobal
	      if ( distance .gt. (sqrt(((dwallbnd(k)%wallx-ielem_xxc(i))**2) &
				     + ((dwallbnd(k)%wally-ielem_yyc(i))**2)))) then
		distance = sqrt(((dwallbnd(k)%wallx-ielem_xxc(i))**2)&
		       +((dwallbnd(k)%wally-ielem_yyc(i))**2))
		ielem_walldist(i) = distance
		ielem_walltrans(i) = dwallbnd(k)%walltrans
		
		 if (ielem_walldist(i).lt.hybridist)then
		    ielem_hybrid(i)=1
		 end if
		
	      end if
	end do

end do

deallocate(dwallbnd)
deallocate(dwallvrt)

end subroutine


subroutine outwritegridbs
!> @brief
!> this subroutine writes the wall surface grid in tecplotm binary format for 3d meshes
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
integer,dimension(70)::ivalid
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog
character*1 nulchar
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
     

inquire (file='surf.plt',exist=herev)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do

end if
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
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
end if
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2

nullptr = 0
      debug   = 0
      filetype = 1
      visdouble = 1
      imax    = igf2
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 3
      soltime = 360.0
      strandid = 0
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
nulchar = char(0)


ierr =  tecini112('simple dataset'//nulchar, &
                    'x y z'//nulchar, &
                    'surf.plt'//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)


 ierr= teczne112('walls'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    null, &
                    null, &
                    shrconn)




	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
    allocate(zbin(igf2))
	if (binio.eq.0)open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	if (binio.eq.1)open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	igf2=0
	if (binio.eq.0)then
        do i=1,imaxn
	read(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	zbin(igf2)=z/scaler
	end if
	end do
	else
	 do i=1,imaxn
	read(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	zbin(igf2)=z/scaler
	end if
	end do
	
	end if

    close(96)


    ierr = tecdat112(igf2,xbin,1)  !!! why not xbin instead of xbin(1) ??
   
    ierr = tecdat112(igf2,ybin,1)

     ierr = tecdat112(igf2,zbin,1)
   
    deallocate(xbin,ybin,zbin)
    
    if (binio.eq.0)open(98,file='GRID.bnd',form='formatted',status='old',action='read')
    if (binio.eq.1)open(98,file='GRID.bnd',form='unformatted',status='old',action='read')
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		if (binio.eq.0)then
		do k=1,imaxb
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		else
		do k=1,imaxb
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		
		end if
	    
 		close(98)
		ierr = tecnod112(icon)
		!ierr = tecnod112(icon)
		deallocate(icon)	
		deallocate(inog)
         
        
  ierr = tecend112()


end if





	call mpi_barrier(mpi_comm_world,ierror)
	
	
	end if
	

	
	

end subroutine outwritegridbs



subroutine outwritegrids
!> @brief
!> this subroutine writes the wall surface grid in tecplotm ascii format for 3d meshes
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
integer,dimension(70)::ivalid
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
     

inquire (file='surf.plt',exist=herev)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then



if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')

allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do

end if
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
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	inog(l)=1
	inog(m)=1
	end if
end do

end if
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
	open(97,file='surf.plt',form='formatted',status='new',action='write')
	write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	write(97,*)'VARIABLES="x","y","z"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'







	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
    allocate(zbin(igf2))
	if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	zbin(igf2)=z/scaler
	end if
	end do

    close(96)
      else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	zbin(igf2)=z/scaler
	end if
	end do

    close(96)
      
      
      end if

    
	write(97,*)xbin(1:igf2)
	write(97,*)ybin(1:igf2)
	write(97,*)zbin(1:igf2)
    
   
   
    deallocate(xbin,ybin,zbin)
    if (binio.eq.0)open(98,file='GRID.bnd',form='formatted',status='old',action='read')
    if (binio.eq.1)open(98,file='GRID.bnd',form='unformatted',status='old',action='read')
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		if (binio.eq.0)then
		do k=1,imaxb
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		else
		do k=1,imaxb
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		
		end if
	    
 		close(98)
 		
 		do i=1,igf2
 		write(97,*)icon(1:4,i)
 		end do
 		
	
		deallocate(icon)	
		deallocate(inog)
         
        close(97)
  
end if





	call mpi_barrier(mpi_comm_world,ierror)
	
	
	end if
	

	
	

end subroutine outwritegrids



subroutine outwritegrids2d
!> @brief
!> this subroutine writes the wall surface grid in tecplot ascii format for 2d meshes
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
integer,dimension(70)::ivalid
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
     

inquire (file='surf.plt',exist=herev)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
end if
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
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do

end if
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2
open(97,file='surf.plt',form='formatted',status='new',action='write')
write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	write(97,*)'VARIABLES="x","y"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FELINE,','DATAPACKING = BLOCK'







	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
	if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96,*)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	end do

    close(96)
      else
      open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96,*)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	end do

    close(96)
      
      end if

    
	write(97,*)xbin(1:igf2)
	write(97,*)ybin(1:igf2)
	
    
   
   
    deallocate(xbin,ybin)
    if (binio.eq.0)open(98,file='GRID.bnd',form='formatted',status='old',action='read')
    if (binio.eq.1)open(98,file='GRID.bnd',form='unformatted',status='old',action='read')
	  allocate(icon(2,totwalls))
    icon=0
     cv=0
		igf2=0
		if (binio.eq.0)then
		do k=1,imaxb
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  
    
		end if
		end do
		else
		do k=1,imaxb
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  
    
		end if
		end do
		
		end if
	    
 		close(98)
 		
 		do i=1,igf2
 		write(97,*)icon(1:2,i)
 		end do
 		
	
		deallocate(icon)	
		deallocate(inog)
         
        close(97)
  
end if





	call mpi_barrier(mpi_comm_world,ierror)
	
	
	end if
	

	
	

end subroutine outwritegrids2d



subroutine outwritegridbs2d
!> @brief
!> this subroutine writes the wall surface grid in tecplot binary format for 2d meshes
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
integer,dimension(70)::ivalid
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
     

inquire (file='surf.plt',exist=herev)



if (herev)then
if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do

end if
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
if (binio.eq.0)open(96,file='GRID.bnd',form='formatted',status='old',action='read')
if (binio.eq.1)open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
allocate (inog(imaxn))
inog(:)=0
if (binio.eq.0)then
do i=1,imaxb
	read(96,*)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do
else
do i=1,imaxb
	read(96)igf,k,j,l,m,o
	if (o.eq.4)then
	inog(k)=1
	inog(j)=1
	
	end if
end do

end if
 close(96)
igf2=0
do i=1,imaxn
      if (inog(i).eq.1)then
      igf2=igf2+1
      end if

end do
	
itotalb=igf2

nullptr = 0
      debug   = 0
      filetype = 1
      visdouble = 1
      imax    = igf2
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 1
      soltime = 360.0
      strandid = 0
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0
nulchar = char(0)


ierr =  tecini112('simple dataset'//nulchar, &
                    'x y'//nulchar, &
                    'surf.plt'//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)


 ierr= teczne112('walls'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    null, &
                    null, &
                    shrconn)




	
    allocate(xbin(igf2))
    allocate(ybin(igf2))
   
	if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96,*)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	end do

    close(96)
      else
      open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96)j,x,y
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      
	xbin(igf2)=x/scaler
	ybin(igf2)=y/scaler
 	
	end if
	end do

    close(96)
      
      end if

    ierr = tecdat112(igf2,xbin,1)  !!! why not xbin instead of xbin(1) ??
   
    ierr = tecdat112(igf2,ybin,1)

    
   
    deallocate(xbin,ybin)
    if (binio.eq.0)open(98,file='GRID.bnd',form='formatted',status='old',action='read')
    if (binio.eq.1)open(98,file='GRID.bnd',form='unformatted',status='old',action='read')
	  allocate(icon(2,totwalls))
    icon=0
     cv=0
		igf2=0
		if (binio.eq.0)then
		do k=1,imaxb
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		 
    
		end if
    !cv=cv+4
        	
		end do
		else
		do k=1,imaxb
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		 
    
		end if
    !cv=cv+4
        	
		end do
		
		end if
	    
 		close(98)
		ierr = tecnod112(icon)
		!ierr = tecnod112(icon)
		deallocate(icon)	
		deallocate(inog)
         
        
  ierr = tecend112()


end if





	call mpi_barrier(mpi_comm_world,ierror)
	
	
	end if
	

	
	

end subroutine outwritegridbs2d


subroutine outwritegrid(n)
!> @brief
!> this subroutine writes the 3d grid file in tecplot ascii format
implicit none
integer,intent(in)::n
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,i6,i7,i8,decomf,kd
integer::inx,i,k,j,m,o,p,q,jk,ixxff
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar




if (n.eq.0)then 

icpuid=n
	kk=0
	write(proc3,fmt='(i10)') ixxff
	outfile='GRID.dat'
! 	outfile='out.'//trim(adjustl(proc3))
	
	
	if (binio.eq.0)then
	
	
	open(97,file=outfile,form='formatted',status='new',action='write')
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	
	
	write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	write(97,*)'VARIABLES="x","y","z"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	
	
	
	do i=1,imaxn
		read(96,*)j,x
		write(97,*)x/scaler
	end do

	close(96)
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	do i=1,imaxn
		read(96,*)j,x,y
		write(97,*)y/scaler
	end do
	close(96)
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	do i=1,imaxn
		read(96,*)j,x,y,z
		write(97,*)z/scaler
	end do
	close (97)
	close(96)

	else
	open(97,file=outfile,form='formatted',status='new',action='write')
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	
	
	write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	write(97,*)'VARIABLES="x","y","z"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEBRICK,','DATAPACKING = BLOCK'
	
	
	
	do i=1,imaxn
		read(96)j,x
		write(97,*)x/scaler
	end do

	close(96)
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	do i=1,imaxn
		read(96)j,x,y
		write(97,*)y/scaler
	end do
	close(96)
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	do i=1,imaxn
		read(96)j,x,y,z
		write(97,*)z/scaler
	end do
	close (97)
	close(96)
	
	
	end if
	if (binio.eq.0)then
	open(96,file='GRID.cel',form='formatted',status='old',action='read')
	open(97,file=outfile,form='formatted',status='old',action='write',position='append')
	do i=1,imaxe
	      read(96,*)ix,i5,i6,i8,i7,i1,i2,i4,i3
	      write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
	end do
	close(96)
	close(97)
	else
	open(96,file='GRID.cel',form='unformatted',status='old',action='read')
	open(97,file=outfile,form='formatted',status='old',action='write',position='append')
	do i=1,imaxe
	      read(96)ix,i5,i6,i8,i7,i1,i2,i4,i3
	      write(97,*)i6,i2,i1,i5,i8,i4,i3,i7
	end do
	close(96)
	close(97)
	
	end if
end if
end subroutine outwritegrid


subroutine outwritegrid2d(n)
!> @brief
!> this subroutine writes the 2d grid file in tecplot ascii format
implicit none
integer,intent(in)::n
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,i6,i7,i8,decomf,kd
integer::inx,i,k,j,m,o,p,q,jk,ixxff
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar




if (n.eq.0)then 

icpuid=n
	kk=0
	write(proc3,fmt='(i10)') ixxff
	outfile='GRID.dat'
! 	outfile='out.'//trim(adjustl(proc3))
	
	
	if (binio.eq.0)then
	open(97,file=outfile,form='formatted',status='new',action='write')
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	
	
	write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	write(97,*)'VARIABLES="x","y"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	
	
	
	do i=1,imaxn
		read(96,*)j,x
		write(97,*)x/scaler
	end do

	close(96)
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	do i=1,imaxn
		read(96,*)j,x,y
		write(97,*)y/scaler
	end do
	close(96)
	
	close (97)
	else
	
	open(97,file=outfile,form='formatted',status='new',action='write')
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	
	
	write(97,*)'TITLE="GRID"'
	write(97,*)'filetype=grid'
	write(97,*)'VARIABLES="x","y"'
	write(97,*)'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	
	
	
	do i=1,imaxn
		read(96)j,x
		write(97,*)x/scaler
	end do

	close(96)
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	do i=1,imaxn
		read(96)j,x,y
		write(97,*)y/scaler
	end do
	close(96)
	
	close (97)

	end if	
	if (binio.eq.0)then
	open(96,file='GRID.cel',form='formatted',status='old',action='read')
	open(97,file=outfile,form='formatted',status='old',action='write',position='append')
	do i=1,imaxe
	      read(96,*)ix,i1,i2,i3,i4
	      write(97,*)i1,i2,i3,i4
	end do
	close(96)
	close(97)
	else
	open(96,file='GRID.cel',form='unformatted',status='old',action='read')
	open(97,file=outfile,form='formatted',status='old',action='write',position='append')
	do i=1,imaxe
	      read(96)ix,i1,i2,i3,i4
	      write(97,*)i1,i2,i3,i4
	end do
	close(96)
	close(97)
	
	end if
	
end if
end subroutine outwritegrid2d


subroutine outwrite3vsb
!> @brief
!> this subroutine writes the 3d surface solution file in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
      integer::iconsidered,facex
      real::shear_temp


 kmaxe=xmpielrank(n)
! ! 
! dumg=totiw
! call mpi_barrier(mpi_comm_world,ierror)
! 
! call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
! imaxp=duml
! 
! allocate(icell(imaxp))
! icell=0
! 
! iloop=0
! 
! ! if (totiw.gt.0)then
! !     
! !    do i=1,totiw
! ! 	k=ibound_t(i)
! ! 	do i1=k,k
! ! 	do j=1,ielem_ifca(k)
! ! 	if (ielem_percorg(j,k).eq.-4)then
! ! 	iloop=iloop+1
! ! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! ! 	go to 1043
! ! 	end if
! ! 	end do
! ! 	end do
! ! 	1043 continue
! !   enddo
! ! end if
! if (totiw.gt.0)then
! do i=1,kmaxe
!   if (ielem_interior(i).eq.1)then
! 	do j=1,ielem_ifca(i)
! 	  if (ielem_ibounds(j,i).gt.0)then
! 	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
! 		  iloop=iloop+1
! 		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
! 	      end if
! 	  end if
! 	end do
!    end if
! end do
! end if
! 
! if (n.eq.0)then
! 	allocate(icella(imaxp*isize))
! 	 icella=0
! 
! end if
! 
! call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)
! 
! ! if (n.eq.0)then
! ! write(1000+n,*)icella(:)
! ! 
! ! 
! ! end if
! 
! call mpi_barrier(mpi_comm_world,ierror)
! deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)


end if



if (itestcase.le.2)then
  nvar1=1
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'solution'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     
     if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar,vortex,k,omega,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar,vortex,mu,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar,vortex,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,vortex,k,omega,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,vortex,mu,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,vortex,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = itotalb
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 3
      !if (( rungekutta .lt. 5).or.( rungekutta .ge. 11)) then
      soltime = t
      !else
      !soltime = it

      !end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(totwalls),xbin2(totwalls))
	

 else
 allocate(xbin2(1))

 end if
 
  totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  allocate(valuess(xmpiwall(n)))
!   end if
    call mpi_barrier(mpi_comm_world,ierror)
   
if (itestcase.le.2)then
					  if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_val(1,1,ibound_t(i))
					  enddo
					  end if


		      
    
    call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
					  if (totiw.gt.0)then
					  do i=1,totiw
!                                                 if (kkd.eq.5)then
!                                                 if (ielem_ishape(ibound_t(i)).eq.2)then
!                                               valuess(i)=-1000
!                                                 else
!                                                 valuess(i)=ielem_dih(kkd,ibound_t(i))
!                                                 end if
!                                                 else
!                                                 valuess(i)=ielem_dih(kkd,ibound_t(i))
!                                                 end if
                                                
						valuess(i)=u_c_val(1,kkd,ibound_t(i))
						if ((kkd.ge.2).and.(kkd.le.4))then
						valuess(i)=u_c_val(1,kkd,ibound_t(i))/u_c_val(1,1,ibound_t(i))
						end if	

                                               

					  enddo
					  end if
		
		     
		
		
		call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		end do
    
					if (totiw.gt.0)then
					  do i=1,totiw
                                                
                                                
                                                
                                                
                                               
					  
                                            
						leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ibound_t(i))
						call cons2prim(n,leftv,mp_pinfl,gammal)
						valuess(i)=leftv(5)
                                                
                                               ! if (ielem_ishape(ibound_t(i)).eq.1)then
                                                        
                                                !        valuess(i)=ielem_dih(6,ibound_t(i))
                                                
                                                !else

                                                 !       valuess(i)=0.0d0
                                               ! end if

						
					  enddo
					  end if
		     
    
    
    
		call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
    
		  if (passivescalar.gt.0)then
					  if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,ibound_t(i))
												
					  enddo
					  end if
		  
		  
		  
		  
		  	  
		  
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  
		  if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=ielem_vortex(1,ibound_t(i))
												
					  enddo
					  end if
		  
		  
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  
		   if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_ct_val(1,kkd,ibound_t(i))
												
					  enddo
					  end if
		  
		  
		  
				  
		   call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		
		  end do
		  end if
		  
		  
		  
		  do kkd=1,3
		  
				      if (totiw.gt.0)then
					  do i=1,totiw
						
					 iconsidered=ibound_t(i)
					 facex=ibound_t2(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x(iconsidered,facex,shear_temp)
					 case (2)

					 call shear_y(iconsidered,facex,shear_temp)
					
                                        case(3)
					 call shear_z(iconsidered,facex,shear_temp)
					 end select
				       		valuess(i)=shear_temp					
					  enddo
					  end if
		  
		  
		  
		  
		  
				  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		
		  end do
		  
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

  if (n.eq.0)then
  ierr = tecend112()
  deallocate(valuelocation,out1,xbin)
  end if
!   if (totiw.gt.0)then
  deallocate (valuess,xbin2)
!   end if











  







	
	
	

	
	

end subroutine outwrite3vsb


subroutine outwrite3vsb2d
!> @brief
!> this subroutine writes the 2d surface solution file in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
      integer::iconsidered,facex
      real::shear_temp


 kmaxe=xmpielrank(n)
! 
totiw=xmpiwall(n)
! call mpi_barrier(mpi_comm_world,ierror)
! 
! call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
! imaxp=duml
! 
! allocate(icell(imaxp))
! icell=0
! 
! iloop=0

! if (totiw.gt.0)then
!     
!    do i=1,totiw
! 	k=ibound_t(i)
! 	do i1=k,k
! 	do j=1,ielem_ifca(k)
! 	if (ielem_percorg(j,k).eq.-4)then
! 	iloop=iloop+1
! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! 	go to 1043
! 	end if
! 	end do
! 	end do
! 	1043 continue
!   enddo
! end if
! if (totiw.gt.0)then
! 
! do i=1,kmaxe
!   if (ielem_interior(i).eq.1)then
! 	do j=1,ielem_ifca(i)
! 	  if (ielem_ibounds(j,i).gt.0)then
! 	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
! 		  iloop=iloop+1
! 		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
! 	      end if
! 	  end if
! 	end do
!    end if
! end do
! end if
! 
! if (n.eq.0)then
! 	allocate(icella(imaxp*isize))
! 	 icella=0
! 
! end if
! 
! call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

! call mpi_barrier(mpi_comm_world,ierror)
! deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)


end if



if (itestcase.le.2)then
  nvar1=1
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'solution'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     
     if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar,vortex,k,omega,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar,vortex,mu,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar,vortex,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,vortex,k,omega,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,vortex,mu,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,vortex,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = itotalb
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 1
      !if (( rungekutta .lt. 5).or.( rungekutta .eq. 11)) then
      soltime = t
      !else
      !soltime = it

      !end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(totwalls),xbin2(totwalls))
	

 else
 allocate(xbin2(1))

 end if

 totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  allocate(valuess(xmpiwall(n)))
!   end if
    
   
if (itestcase.le.2)then
		     
					  if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_val(1,1,ibound_t(i))
					  enddo
					  end if
    
    
     
    
    
    
     call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,4
		     if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_val(1,kkd,ibound_t(i))
							    if ((kkd.ge.2).and.(kkd.le.3))then
			      valuess(i)=u_c_val(1,kkd,ibound_t(i))/u_c_val(1,1,ibound_t(i))
			      end if	
						
						
					  enddo
		  end if
		
		
		      
		
		
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		
		end do
    
		      if (totiw.gt.0)then
					  do i=1,totiw
					  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ibound_t(i))
				    call cons2prim(n,leftv,mp_pinfl,gammal)
				    valuess(i)=leftv(4)
					end do	
			      end if	
			
    
    
    
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		
    
		  if (passivescalar.gt.0)then
		  
		  if (totiw.gt.0)then
					  do i=1,totiw
					valuess(i)=u_ct_val(1,turbulenceequations+passivescalar,ibound_t(i))
					end do	
			      end if	
		  
		  
		   
		  
		  
		  	  
		  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  if (totiw.gt.0)then
					  do i=1,totiw
					valuess(i)=ielem_vortex(1,ibound_t(i))
					end do	
			      end if	
		 
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		   call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		   if (totiw.gt.0)then
					  do i=1,totiw
					valuess(i)=u_ct_val(1,kkd,ibound_t(i))
					end do	
			      end if
		  
		  
				  
		   call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		end do
		  end if
		  
		  
		  
		  do kkd=1,2
			    if (totiw.gt.0)then
					  do i=1,totiw
					  iconsidered=ibound_t(i)
					  facex=ibound_t2(i)
					  select case(kkd)
					 case(1)
					 
					 call shear_x2d(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y2d(iconsidered,facex,shear_temp)
					 
					 end select
					valuess(i)=shear_temp
					end do	
			      end if
		 
		  
				  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  
		  
		  
		  
		  
		  
		  end do
		  
		  end if
    
    
    end if





!   end if
!    end if

    if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1)
  end if
  
!   if (totiw.gt.0)then
  deallocate (valuess,xbin2)
!   end if
  











  







	
	
	

	
	

end subroutine outwrite3vsb2d



subroutine outwrite3vs
!> @brief
!> this subroutine writes the 3d surface solution file in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
integer::iconsidered,facex
      real::shear_temp

 kmaxe=xmpielrank(n)
! 
dumg=totiw
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

iloop=0

! if (totiw.gt.0)then
!     
!    do i=1,totiw
! 	k=ibound_t(i)
! 	do i1=k,k
! 	do j=1,ielem_ifca(k)
! 	if (ielem_percorg(j,k).eq.-4)then
! 	iloop=iloop+1
! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! 	go to 1043
! 	end if
! 	end do
! 	end do
! 	1043 continue
!   enddo
! end if
if (totiw.gt.0)then
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		  iloop=iloop+1
		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
	      end if
	  end if
	end do
   end if
end do
end if

if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)

open(97,file=outfile,form='formatted',status='new',action='write')
end if








		
   if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if





if (itestcase.le.2)then
  nvar1=1
  
  
   if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="solution"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 
                    
                    
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if
     else
     
     
     if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      
                    
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","K","OMEGA","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED,[13] = CELLCENTERED)'
  end if
                    
                    
                    
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","MU","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
              
                if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","K","OMEGA","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","m","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
         
               if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED)'
  end if
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







     

valuelocation(:)=0





 


 allocate(valuesa(imaxp*isize))
  allocate(xbin(totwalls))
	valuesa=0.0

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=0.0
    
   
if (itestcase.le.2)then
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(1,1,ibound_which(i))
				end if
		      end do
		      end if
    
    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    write(97,*)xbin(1:totwalls)
    
    end if

     
    call mpi_barrier(mpi_comm_world,ierror)
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(1,kkd,ibound_which(i))
				if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(icount_wall)=u_c_val(1,kkd,ibound_which(i))/u_c_val(1,1,ibound_which(i))
		  end if	
				end if
		      end do
		      end if
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
    
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ibound_which(i))
				    call cons2prim(n,leftv,mp_pinfl,gammal)
				    valuess(icount_wall)=leftv(5)
				  
				end if
		      end do
		      end if
    
    
    
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
    
		  if (passivescalar.gt.0)then
		   if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(1,turbulenceequations+passivescalar,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  	  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=ielem_vortex(1,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(1,kkd,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,3
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
					iconsidered=ibound_which(i)
					 facex=ibound_face(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y(iconsidered,facex,shear_temp)
					 case(3)
					 call shear_z(iconsidered,facex,shear_temp)
					 end select
				       		valuess(icount_wall)=shear_temp		
				       
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

  if (n.eq.0)then
  
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  











  call mpi_barrier(mpi_comm_world,ierror)







	
	
	

	
	

end subroutine outwrite3vs


subroutine outwrite3vs2d
!> @brief
!> this subroutine writes the 2d surface solution file in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
integer::iconsidered,facex
      real::shear_temp

 kmaxe=xmpielrank(n)
! 
dumg=totiw
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

iloop=0

! if (totiw.gt.0)then
!     
!    do i=1,totiw
! 	k=ibound_t(i)
! 	do i1=k,k
! 	do j=1,ielem_ifca(k)
! 	if (ielem_percorg(j,k).eq.-4)then
! 	iloop=iloop+1
! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! 	go to 1043
! 	end if
! 	end do
! 	end do
! 	1043 continue
!   enddo
! end if
if (totiw.gt.0)then
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
		  iloop=iloop+1
		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
	      end if
	  end if
	end do
   end if
end do
end if

if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)

	open(97,file=outfile,form='formatted',status='new',action='write')
	




end if



if (itestcase.le.2)then
  nvar1=1
  
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="solution"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
 
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 
                    
                     if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     else
      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED)'
  end if
!      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
!                     'density,u,v,energy,pressure'//nulchar, &
!                     out1//nulchar, &
!                     '.'//nulchar, &
!                     filetype, &
!                     debug, &
!                     visdouble)
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","K","OMEGA","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED, [11] = CELLCENTERED)'
  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                   if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","m","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
               if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED)'
  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	             if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","K","OMEGA","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
! 	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
!                     'density,u,v,energy,pressure,vortex,k,omega,ssx,ssy'//nulchar, &
!                     out1//nulchar, &
!                     '.'//nulchar, &
!                     filetype, &
!                     debug, &
!                     visdouble)
              end if
              if (turbulenceequations.eq.1)then
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","MU","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED)'
	
  end if
              
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = itotalb
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 1
     ! if (( rungekutta .lt. 5).or.( rungekutta .eq. 11)) then
      soltime = t
      !else
      !soltime = it

      !end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





!  ierr= teczne112('grid2'//nulchar, &
!                     zonetype, &
!                     imax, &
!                     jmax, &
!                     kmax, &
!                     icellmax, &
!                     jcellmax, &
!                     kcellmax, &
!                     soltime, &
!                     strandid, &
!                     parentzn, &
!                     isblock, &
!                     nfconns, &
!                     fnmode, &
!                     0, &
!                     0, &
!                     0, &
!                     null, &
!                     valuelocation, &
!                     null, &
!                     shrconn)


 allocate(valuesa(imaxp*isize))
  allocate(xbin(totwalls))
	valuesa=0.0

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=0.0
    
   
if (itestcase.le.2)then
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(1,1,ibound_which(i))
				end if
		      end do
		      end if
    
    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    write(97,*)xbin(1:totwalls)
!     ierr = tecdat112(totwalls,xbin,1)
    end if

     
    call mpi_barrier(mpi_comm_world,ierror)
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,4
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(1,kkd,ibound_which(i))
				if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(icount_wall)=u_c_val(1,kkd,ibound_which(i))/u_c_val(1,1,ibound_which(i))
		  end if	
				end if
		      end do
		      end if
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
    
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ibound_which(i))
				    call cons2prim(n,leftv,mp_pinfl,gammal)
				    valuess(icount_wall)=leftv(4)
				  
				end if
		      end do
		      end if
    
    
    
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
    
		  if (passivescalar.gt.0)then
		   if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(1,turbulenceequations+passivescalar,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  	  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=ielem_vortex(1,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(1,kkd,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,2
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
					iconsidered=ibound_which(i)
					 facex=ibound_face(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x2d(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y2d(iconsidered,facex,shear_temp)
					 
					 end select
				       		valuess(icount_wall)=shear_temp		
				       
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

  if (n.eq.0)then
 
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  











  call mpi_barrier(mpi_comm_world,ierror)







	
	
	

	
	

end subroutine outwrite3vs2d




subroutine outwrite3vbav
!> @brief
!> this subroutine writes the 3d averaged solution file in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(13))

kmaxe=xmpielrank(n)





if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="VOL_AVER_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	

end if
call mpi_barrier(mpi_comm_world,ierror)

  nvar1=11
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'r_mean,u_mean,v_mean,w_mean,p_mean,u_rms,v_rms,w_rms,uv,uw,wv'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 

	

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 5
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0

 



 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(imaxe),xbin2(imaxe))
	

 end if

  
  allocate(valuess(kmaxe))
  valuess=zero
    
   

	      
	      

    
    
    		do kkd=1,5
		do i=1,kmaxe
		  valuess(i)=u_c_val(ind1,kkd,i)
		  if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(i)=u_c_val(ind1,kkd,i)/u_c_val(ind1,1,i)
		  end if
		  if (kkd.eq.5)then
		   leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		  
		  end if
		  
		end do
		
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		
		end do
		do kkd=1,6

    		do i=1,kmaxe

		  valuess(i)=u_c_rms(kkd,i)
		end do
		
		
		call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if
		end do
    
		  
    
    
      
    
    
    

     
    
    
    
  if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1,xbin2)
  end if
  
  deallocate (valuess,variables)

  call mpi_barrier(mpi_comm_world,ierror)









	
	

end subroutine outwrite3vbav

subroutine outwrite3vb2dav
!> @brief
!> this subroutine writes the 2d averaged solution file in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(10))

kmaxe=xmpielrank(n)

dumg=kmaxe
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="VOL_AVER_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	

end if
call mpi_barrier(mpi_comm_world,ierror)

  nvar1=7
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'r_mean,u_mean,v_mean,p_mean,u_rms,v_rms,uv'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 

	

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 3
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0

 



 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 allocate(valuesa(imaxp*isize))
  allocate(xbin(imaxe))
	valuesa=zero

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=zero
    

	      


    
    
    		do kkd=1,4
		do i=1,kmaxe
		  valuess(i)=u_c_rms(kkd,i)
		  if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(i)=u_c_val(ind1,kkd,i)/u_c_val(ind1,1,i)
		  end if
		  if (kkd.eq.4)then
		   leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(4)
		  
		  end if
		  
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
		do kkd=1,3
    		do i=1,kmaxe
		  
		  valuess(i)=u_c_val(5,kkd,i)
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		ierr = tecdat112(imaxe,xbin,1)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
		  
    
    
      
    
    
    

     
    
    
    
  if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)







         
        
 deallocate(variables)


	
	

end subroutine outwrite3vb2dav



subroutine outwrite3vav
!> @brief
!> this subroutine writes the 3d averaged solution file in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(13))

kmaxe=xmpielrank(n)

dumg=kmaxe
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="VOL_AVER_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	open(97,file=outfile,form='formatted',status='new',action='write')
 write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="r_mean","u_mean","v_mean","w_mean","p_mean","u_rms","v_rms","w_rms","uv","uw","wv"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED)'
end if
call mpi_barrier(mpi_comm_world,ierror)

  
  
      
  
  
 

	

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))









 allocate(valuesa(imaxp*isize))
  allocate(xbin(imaxe))
	valuesa=zero

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=zero
    
   


    
    
    		do kkd=1,5
		do i=1,kmaxe
		  valuess(i)=u_c_rms(kkd,i)
		  if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(i)=u_c_val(5,kkd,i)/u_c_val(5,1,i)
		  end if
		  if (kkd.eq.5)then
		   leftv(1:nof_variables)=u_c_val(5,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		  
		  end if
		  
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
		do kkd=1,6
    		do i=1,kmaxe
		  
		  valuess(i)=u_c_val(5,kkd,i)
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
		  
    
    
      
    
    
    

     
    
    
    
  if (n.eq.0)then
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)







         
        
 deallocate(variables)


	
	

end subroutine outwrite3vav

subroutine outwrite3v2dav
!> @brief
!> this subroutine writes the 2d averaged solution file in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(10))

kmaxe=xmpielrank(n)

dumg=kmaxe
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then





 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="VOL_AVER_"//trim(adjustl(proc3))//".plt"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)
	open(97,file=outfile,form='formatted',status='new',action='write')

end if
call mpi_barrier(mpi_comm_world,ierror)


  
  if (n.eq.0)then
        write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="r_mean","u_mean","v_mean","p_mean","u_rms","v_rms","uv"'
	write(97,*) 'Zone N=',imaxn,',E=',imaxe,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if


  nvar1=7
  
  
  
  
 

	

	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = imaxn
      jmax    = imaxe
      kmax    = 0
      ZONETYPE = 3
      
      soltime = t
      
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0

 






 allocate(valuesa(imaxp*isize))
  allocate(xbin(imaxe))
	valuesa=zero

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=zero
    
   


    
    
    		do kkd=1,4
		do i=1,kmaxe
		  valuess(i)=u_c_val(5,kkd,i)
		  if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(i)=u_c_val(5,kkd,i)/u_c_val(5,1,i)
		  end if
		  if (kkd.eq.4)then
		   leftv(1:nof_variables)=u_c_val(5,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(4)
		  
		  end if
		  
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*) xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
		do kkd=1,3
    		do i=1,kmaxe
		  
		  valuess(i)=u_c_val(5,1,i)
		end do
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*) xbin(1:imaxe)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
		  
    
    
      
    
    
    

     
    
    
    
  if (n.eq.0)then
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)







         
        
 deallocate(variables)


	
	

end subroutine outwrite3v2dav



subroutine outwrite3vsbav
!> @brief
!> this subroutine writes the 3d surface averaged solution file in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
integer::iconsidered,facex
      real::shear_temp

 kmaxe=xmpielrank(n)
! 

! if (totiw.gt.0)then
!     
!    do i=1,totiw
! 	k=ibound_t(i)
! 	do i1=k,k
! 	do j=1,ielem_ifca(k)
! 	if (ielem_percorg(j,k).eq.-4)then
! 	iloop=iloop+1
! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! 	go to 1043
! 	end if
! 	end do
! 	end do
! 	1043 continue
!   enddo
! end if




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_AV"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)


end if



if (itestcase.le.2)then
  nvar1=1
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'solution'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     
     if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar,vortex,k,omega,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar,vortex,mu,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,passivescalar,vortex,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,vortex,k,omega,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,vortex,mu,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,w,energy,pressure,vortex,ssx,ssy,ssz'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = itotalb
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 3
!       if (( rungekutta .lt. 5).or.( rungekutta .eq. 11)) then
      soltime = t
!       else
!       soltime = it

!       end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 
  allocate(xbin(totwalls),xbin2(totwalls))
	

 end if

  totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  allocate(valuess(xmpiwall(n)))
!   end if

   
if (itestcase.le.2)then

		      
					  if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_rms(j,ibound_t(i))
					  enddo
					  end if
    
   call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
					  if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_val(ind1,kkd,ibound_t(i))
						if ((kkd.ge.2).and.(kkd.le.4))then
					      valuess(i)=u_c_val(ind1,kkd,ibound_t(i))/u_c_val(ind1,1,ibound_t(i))
					      end if	
						
						
					  enddo
					  end if
		
		
		
		
		call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		end do
    
					  if (totiw.gt.0)then
					  do i=1,totiw
					        leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,ibound_t(i))
						
						call cons2prim(n,leftv,mp_pinfl,gammal)
				    valuess(i)=leftv(5)
						
						
						
					  enddo
					  end if
		     
    
    
    
		call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
    
		  if (passivescalar.gt.0)then
		  
					if (totiw.gt.0)then
					  do i=1,totiw
					        
				    valuess(i)=u_ct_val(ind1,turbulenceequations+passivescalar,ibound_t(i))
						
						
						
					  enddo
					  end if
		  
		  
		   
		  
		  
		  	  
		  
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  if (totiw.gt.0)then
					  do i=1,totiw
					        
				    valuess(i)=ielem_vortex(1,ibound_t(i))
						
						
						
					  enddo
					  end if
		  
		  
		  
		 
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
					  if (totiw.gt.0)then
					  do i=1,totiw
					        
				    valuess(i)=u_ct_val(ind1,kkd,ibound_t(i))
						
						
						
					  enddo
					  end if
		  
		  
		
		  
				  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  end do
		  end if
		  
		  
		  if (itestcase.eq.4)then
		  do kkd=1,3
		  
				      if (totiw.gt.0)then
					  do i=1,totiw
					        
				    
						
						iconsidered=ibound_t(i)
					 facex=ibound_t2(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x_av(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y_av(iconsidered,facex,shear_temp)
					 case(3)
					 call shear_z_av(iconsidered,facex,shear_temp)
					 end select
				       		valuess(i)=shear_temp		
						
					  enddo
					  end if
		  
		  
		  
		  
		  
				  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

  if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1,xbin2)
  end if
  
!  if (totiw.gt.0)then
  deallocate (valuess)
!   end if
  











 







	
	
	

	
	

end subroutine outwrite3vsbav


subroutine outwrite3vsb2dav
!> @brief
!> this subroutine writes the 2d surface averaged solution file in tecplot binary format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
integer::iconsidered,facex
      real::shear_temp

 kmaxe=xmpielrank(n)
 
! 
! dumg=totiw
! call mpi_barrier(mpi_comm_world,ierror)
! 
! call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
! imaxp=duml
! 
! allocate(icell(imaxp))
! icell=0
! 
! iloop=0
! 
! ! if (totiw.gt.0)then
! !     
! !    do i=1,totiw
! ! 	k=ibound_t(i)
! ! 	do i1=k,k
! ! 	do j=1,ielem_ifca(k)
! ! 	if (ielem_percorg(j,k).eq.-4)then
! ! 	iloop=iloop+1
! ! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! ! 	go to 1043
! ! 	end if
! ! 	end do
! ! 	end do
! ! 	1043 continue
! !   enddo
! ! end if
! if (totiw.gt.0)then
! do i=1,kmaxe
!   if (ielem_interior(i).eq.1)then
! 	do j=1,ielem_ifca(i)
! 	  if (ielem_ibounds(j,i).gt.0)then
! 	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
! 		  iloop=iloop+1
! 		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
! 	      end if
! 	  end if
! 	end do
!    end if
! end do
! end if
! 
! if (n.eq.0)then
! 	allocate(icella(imaxp*isize))
! 	 icella=0
! 
! end if
! 
! call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)
! 
! ! if (n.eq.0)then
! 
! ! 
! ! 
! ! end if
! 
! call mpi_barrier(mpi_comm_world,ierror)
! deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_AV"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)


end if



if (itestcase.le.2)then
  nvar1=1
  if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'solution'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
  
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     else
     
     if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar,vortex,k,omega,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar,vortex,mu,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,passivescalar,vortex,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,vortex,k,omega,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.1)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,vortex,mu,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              end if
              if (turbulenceequations.eq.0)then
              if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
                    'density,u,v,energy,pressure,vortex,ssx,ssy'//nulchar, &
                    out1//nulchar, &
                    '.'//nulchar, &
                    filetype, &
                    debug, &
                    visdouble)
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = itotalb
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 1
!       if (( rungekutta .lt. 5).or.( rungekutta .eq. 11)) then
      soltime = t
!       else
!       soltime = it

!       end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





 ierr= teczne112('grid2'//nulchar, &
                    zonetype, &
                    imax, &
                    jmax, &
                    kmax, &
                    icellmax, &
                    jcellmax, &
                    kcellmax, &
                    soltime, &
                    strandid, &
                    parentzn, &
                    isblock, &
                    nfconns, &
                    fnmode, &
                    0, &
                    0, &
                    0, &
                    null, &
                    valuelocation, &
                    null, &
                    shrconn)


 allocate(xbin(totwalls),xbin2(totwalls))
	

 end if

 totiw=xmpiwall(n)
!   if (xmpiwall(n).gt.0)then
  allocate(valuess(xmpiwall(n)))
!   end if
    
   
if (itestcase.le.2)then

		      if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_val(5,1,ibound_t(i))
					  enddo
					  end if
		     
    
   call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,4
		      if (totiw.gt.0)then
					  do i=1,totiw
						valuess(i)=u_c_val(5,kkd,ibound_t(i))
						if ((kkd.ge.2).and.(kkd.le.3))then
						valuess(i)=u_c_val(5,kkd,ibound_t(i))/u_c_val(5,1,ibound_t(i))
						end if	
					  enddo
					  end if
		     
		
		
		call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
		end do
    
		      if (totiw.gt.0)then
					  do i=1,totiw
						leftv(1:nof_variables)=u_c_val(5,1:nof_variables,ibound_t(i))
						call cons2prim(n,leftv,mp_pinfl,gammal)
					  valuess(i)=leftv(4)
					  enddo
					  end if
		      
    
    
    
		call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
    
		  if (passivescalar.gt.0)then
		    if (totiw.gt.0)then
					  do i=1,totiw
						
						
					  valuess(i)=u_ct_val(5,turbulenceequations+passivescalar,ibound_t(i))
					  enddo
					  end if
		  
		  
		  
		  
		  	  
		  
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		   if (totiw.gt.0)then
					  do i=1,totiw
						
						
					  valuess(i)=ielem_vortex(1,ibound_t(i))
					  enddo
					  end if
		
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
					  if (totiw.gt.0)then
					  do i=1,totiw
						
						
					  valuess(i)=u_ct_val(5,kkd,ibound_t(i))
					  enddo
					  end if
		  
		  
				  
		 call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,2
				      if (totiw.gt.0)then
					  do i=1,totiw
						iconsidered=ibound_t(i)
					 facex=ibound_t2(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x2d_av(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y2d_av(iconsidered,facex,shear_temp)
					 
					 end select
				       		valuess(i)=shear_temp		
						
					  
					  enddo
					  end if
		 
		  
				  
		  call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		ierr = tecdat112(totwalls,xbin,1)
		end if
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

   if (n.eq.0)then
  ierr = tecend112()
  deallocate(xbin,valuelocation,out1,xbin2)
  end if
  
!   if (totiw.gt.0)then
  deallocate (valuess)
!   end if






	
	
	

	
	

end subroutine outwrite3vsb2dav


subroutine outwrite3vsav
!> @brief
!> this subroutine writes the 3d surface averaged solution file in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
integer::iconsidered,facex
      real::shear_temp

 kmaxe=xmpielrank(n)
! 
dumg=totiw
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

iloop=0

! if (totiw.gt.0)then
!     
!    do i=1,totiw
! 	k=ibound_t(i)
! 	do i1=k,k
! 	do j=1,ielem_ifca(k)
! 	if (ielem_percorg(j,k).eq.-4)then
! 	iloop=iloop+1
! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! 	go to 1043
! 	end if
! 	end do
! 	end do
! 	1043 continue
!   enddo
! end if
if (totiw.gt.0)then
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
		  iloop=iloop+1
		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
	      end if
	  end if
	end do
   end if
end do
end if

if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_AV"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)

open(97,file=outfile,form='formatted',status='new',action='write')
end if








		
   if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if





if (itestcase.le.2)then
  nvar1=1
  
  
   if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="solution"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 
                    
                    
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED)'
  end if
     else
     
     
     if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	      
                    
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","K","OMEGA","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED,[13] = CELLCENTERED)'
  end if
                    
                    
                    
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","MU","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
              
                if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","passivescalar","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
	    else
	    if (turbulenceequations.eq.2)then
	      
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","K","OMEGA","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED, [11] = CELLCENTERED,[12] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.1)then
              
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","m","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED,[11] = CELLCENTERED)'
  end if
              end if
              if (turbulenceequations.eq.0)then
         
               if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","W","energy","Pressure","VORTEX","ssx","ssy","ssz"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED,[7] = CELLCENTERED,[8] = CELLCENTERED,[9] = CELLCENTERED,'
	write(97,*) '[10] = CELLCENTERED)'
  end if
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







     

valuelocation(:)=0





 


 allocate(valuesa(imaxp*isize))
  allocate(xbin(totwalls))
	valuesa=0.0

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=0.0
    
   
if (itestcase.le.2)then
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(5,1,ibound_which(i))
				end if
		      end do
		      end if
    
    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    write(97,*)xbin(1:totwalls)
    
    end if

     
    call mpi_barrier(mpi_comm_world,ierror)
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,5
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(5,kkd,ibound_which(i))
				if ((kkd.ge.2).and.(kkd.le.4))then
		  valuess(icount_wall)=u_c_val(5,kkd,ibound_which(i))/u_c_val(5,1,ibound_which(i))
		  end if	
				end if
		      end do
		      end if
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
    
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_variables)=u_c_val(5,1:nof_variables,ibound_which(i))
				    call cons2prim(n,leftv,mp_pinfl,gammal)
				    valuess(icount_wall)=leftv(5)
				  
				end if
		      end do
		      end if
    
    
    
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		 write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
    
		  if (passivescalar.gt.0)then
		   if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(5,turbulenceequations+passivescalar,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  	  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=ielem_vortex(1,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(5,kkd,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,3
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
					iconsidered=ibound_which(i)
					 facex=ibound_face(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x_av(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y_av(iconsidered,facex,shear_temp)
					 case(3)
					 call shear_z_av(iconsidered,facex,shear_temp)
					 end select
				       		valuess(icount_wall)=shear_temp		
				       
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		   write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

  if (n.eq.0)then
  
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  











  call mpi_barrier(mpi_comm_world,ierror)







	
	
	

	
	

end subroutine outwrite3vsav


subroutine outwrite3vs2dav
!> @brief
!> this subroutine writes the 2d surface averaged solution file in tecplot ascii format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm,nx,ny,nz,ssx,ssy,ssz,ssp,tauyx,tauzx,tauzy
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd,im
real,dimension(8)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,icount_wall
logical::herev
real,dimension(5)::total
 character(len=20)::proc,outfile,proc3,surfile,proc4
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnod112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype,iloop
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
      real::shear_temp
integer::iconsidered,facex

 kmaxe=xmpielrank(n)
! 
dumg=totiw
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

iloop=0

! if (totiw.gt.0)then
!     
!    do i=1,totiw
! 	k=ibound_t(i)
! 	do i1=k,k
! 	do j=1,ielem_ifca(k)
! 	if (ielem_percorg(j,k).eq.-4)then
! 	iloop=iloop+1
! 	icell(iloop)=ielem_indexi(j,ibound_t(k))
! 	go to 1043
! 	end if
! 	end do
! 	end do
! 	1043 continue
!   enddo
! end if
if (totiw.gt.0)then
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
		  iloop=iloop+1
		    icell(iloop)=ibound_inum(ielem_ibounds(j,i))
	      end if
	  end if
	end do
   end if
end do
end if

if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

! if (n.eq.0)then

! 
! 
! end if

call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)




if (n.eq.0)then



 nullptr = 0
      debug   = 0
      filetype = 2
visdouble = 1

nulchar = char(0)

write(proc3,fmt='(i10)') it
	!proc4=".plt"
	outfile="SURF_AV"//trim(adjustl(proc3))//'.plt'
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
! 	out1=out1//char(0)

	open(97,file=outfile,form='formatted',status='new',action='write')
	




end if



if (itestcase.le.2)then
  nvar1=1
  
                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="solution"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED)'
  end if
 
  
  
 end if
 if (itestcase.eq.3)then
 nvar1=6+passivescalar
  if (passivescalar.gt.0)then
 
                    
                     if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED)'
  end if
     else
      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED)'
  end if
!      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
!                     'density,u,v,energy,pressure'//nulchar, &
!                     out1//nulchar, &
!                     '.'//nulchar, &
!                     filetype, &
!                     debug, &
!                     visdouble)
     
     
     
     end if
 end if
 if (itestcase.eq.4)then
 nvar1=10+passivescalar+turbulenceequations
	    if (passivescalar.gt.0)then
	      if (turbulenceequations.eq.2)then
	                    if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","K","OMEGA","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED, [11] = CELLCENTERED)'
  end if
	      
              end if
              if (turbulenceequations.eq.1)then
                   if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","m","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
               if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","passivescalar","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED)'
  end if
              
              
              end if
	    else
	    if (turbulenceequations.eq.2)then
	             if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","K","OMEGA","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED, [10] = CELLCENTERED)'
  end if
! 	      if (n.eq.0)ierr =  tecini112('sols'//nulchar, &
!                     'density,u,v,energy,pressure,vortex,k,omega,ssx,ssy'//nulchar, &
!                     out1//nulchar, &
!                     '.'//nulchar, &
!                     filetype, &
!                     debug, &
!                     visdouble)
              end if
              if (turbulenceequations.eq.1)then
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","MU","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED,'
	write(97,*) '[9] = CELLCENTERED)'
  end if
              
              end if
              if (turbulenceequations.eq.0)then
                      if (n.eq.0)then
	write(97,*) 'FILETYPE=SOLUTION'
	write(97,*)'VARIABLES="Density","U","V","energy","Pressure","VORTEX","ssx","ssy"'
	write(97,*) 'Zone N=',itotalb,',E=',totwalls,',ZONETYPE = FEQUADRILATERAL,','DATAPACKING = BLOCK'
	write(97,*) ',VARLOCATION = ([1] = CELLCENTERED,[2] = CELLCENTERED,[3] = CELLCENTERED,[4] = CELLCENTERED,'
	write(97,*) '[5] = CELLCENTERED, [6] = CELLCENTERED, [7] = CELLCENTERED, [8] = CELLCENTERED)'
	
  end if
              
              
              end if
	    
	    
	    
	    
	    end if
 
 
 end if
	
	

if (n.eq.0)then
allocate (valuelocation(nvar1))







      imax    = itotalb
      jmax    = totwalls
      kmax    = 0
      ZONETYPE = 1
!       if (( rungekutta .lt. 5).or.( rungekutta .eq. 11)) then
      soltime = t
!       else
!       soltime = it

!       end if
      strandid = 1
      parentzn = 0
      isblock = 1
      icellmax = 0
      jcellmax = 0
      kcellmax = 0
      nfconns = 0
      fnmode = 0
      shrconn = 0

valuelocation(:)=0





!  ierr= teczne112('grid2'//nulchar, &
!                     zonetype, &
!                     imax, &
!                     jmax, &
!                     kmax, &
!                     icellmax, &
!                     jcellmax, &
!                     kcellmax, &
!                     soltime, &
!                     strandid, &
!                     parentzn, &
!                     isblock, &
!                     nfconns, &
!                     fnmode, &
!                     0, &
!                     0, &
!                     0, &
!                     null, &
!                     valuelocation, &
!                     null, &
!                     shrconn)


 allocate(valuesa(imaxp*isize))
  allocate(xbin(totwalls))
	valuesa=0.0

 end if

  call mpi_barrier(mpi_comm_world,ierror)
  allocate(valuess(imaxp))
  valuess=0.0
    
   
if (itestcase.le.2)then
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(5,1,ibound_which(i))
				end if
		      end do
		      end if
    
    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i))=valuesa(i)
	end if
    end do
    write(97,*)xbin(1:totwalls)
!     ierr = tecdat112(totwalls,xbin,1)
    end if

     
    call mpi_barrier(mpi_comm_world,ierror)
    end if
    
    if (itestcase.ge.3)then
		do kkd=1,4
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    valuess(icount_wall)=u_c_val(5,kkd,ibound_which(i))
				if ((kkd.ge.2).and.(kkd.le.3))then
		  valuess(icount_wall)=u_c_val(5,kkd,ibound_which(i))/u_c_val(5,1,ibound_which(i))
		  end if	
				end if
		      end do
		      end if
		
		
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
		 
		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
		end do
    
    
		      if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				    leftv(1:nof_variables)=u_c_val(5,1:nof_variables,ibound_which(i))
				    call cons2prim(n,leftv,mp_pinfl,gammal)
				    valuess(icount_wall)=leftv(4)
				  
				end if
		      end do
		      end if
    
    
    
		call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		if (n.eq.0)then
		do i=1,imaxp*isize
		    if (icella(i).gt.0)then
		    xbin(icella(i))=valuesa(i)
		    end if
		end do
		write(97,*)xbin(1:totwalls)
		end if

		
		call mpi_barrier(mpi_comm_world,ierror)
    
		  if (passivescalar.gt.0)then
		   if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(5,turbulenceequations+passivescalar,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  	  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  
		  end if
    
		  if (itestcase.eq.4)then
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=ielem_vortex(1,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
		  
		  
! 		  do i=1,kmaxe
! 		      valuess(i)=ielem_vortex(1,i)
! 		  end do
		  
		  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  
		  if (turbulence.eq.1)then
		  do kkd=1,turbulenceequations
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
				       valuess(icount_wall)=u_ct_val(5,kkd,ibound_which(i))
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  if (turbulence.eq.1)then
		  do kkd=1,2
		  
		  if (totiw.gt.0)then
		      icount_wall=0	
		      do i=1,n_boundaries
			      if ((ibound_icode(i).eq.4).and.(ibound_which(i).gt.0)) then
				    icount_wall=icount_wall+1
					iconsidered=ibound_which(i)
					 facex=ibound_face(i)
					select case(kkd)
					 case(1)
					 
					 call shear_x2d_av(iconsidered,facex,shear_temp)
					 case (2)
					 call shear_y2d_av(iconsidered,facex,shear_temp)
					 
					 end select
				       		valuess(icount_wall)=shear_temp		
				       
				  
				end if
		      end do
		      end if
		  
				  
		  call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

		  if (n.eq.0)then
		  do i=1,imaxp*isize
		      if (icella(i).gt.0)then
		      xbin(icella(i))=valuesa(i)
		      end if
		  end do
		  write(97,*)xbin(1:totwalls)
		  end if

		  
		  call mpi_barrier(mpi_comm_world,ierror)
		  end do
		  end if
		  
		  
		  
		  
		  
		  
		  
		  
		  end if
    
    
    end if





!   end if
!    end if

  if (n.eq.0)then
 
  deallocate(xbin,valuesa,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  











  call mpi_barrier(mpi_comm_world,ierror)







	
	
	

	
	

end subroutine outwrite3vs2dav






subroutine grid_write
!> @brief
!> this subroutine calls the appropriate grid writing subroutine based on the settings
implicit none



if (tecplot.eq.1)then		!binary tecplot
  if (dimensiona.eq.3)then
  call outwritegridb
  else
  call outwritegridb2d
  end if

end if
if (tecplot.eq.0)then		!ascii tecplot
    if (dimensiona.eq.3)then
      call outwritegrid(n)
  
    else
      call outwritegrid2d(n)

    end if
  
end if

if (tecplot.eq.2)then		!binary paraview 3d only

  if (dimensiona.eq.3)then

        call outwritepara3db
    else
    
        call outwritepara2db
    end if

end if


if (tecplot.eq.3)then		!binary paraview 3d only

  call outwritepara3dbp

end if


if (tecplot.eq.4)then		!binary tecplot partitioned 3d only

  call outwritetec3dbp
call mpi_barrier(mpi_comm_world,ierror)
end if

if (tecplot.eq.5)then		!fast output written by all processors using mpi-io

call parallel_vtk_combine(n)

end if

if (tecplot.eq.6)then		!fast output written by all processors



call parallel_vtk_combine_partitioned(n)

end if






end subroutine grid_write


subroutine surf_write
!> @brief
!> this subroutine calls the appropriate surface writing subroutine based on the settings
implicit none
if (tecplot.eq.1)then
   if (dimensiona.eq.3)then
call outwritegridbs
else
call outwritegridbs2d
end if
end if
if (tecplot.eq.0)then
 if (dimensiona.eq.3)then
call outwritegrids
else
call outwritegrids2d
end if

end if


if (tecplot.eq.3)then		!binary paraview 3d only

  call outwritepara3dsb

end if


if (tecplot.eq.5)then		!fast output written by all processors



call parallel_vtk_combine_wall(n)

end if


if (tecplot.eq.6)then		!fast output written by all processors



call parallel_vtk_combine_partitioned_wall(n)

end if




end subroutine surf_write


subroutine volume_solution_write
!> @brief
!> this subroutine calls the appropriate volume writing subroutine based on the settings
implicit none


#ifdef gpu
				call gpu_to_host_vol
#endif




				  

				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"output1",t
				  close(63)
				  end if
	  if (tecplot.eq.1)then
				if (dimensiona.eq.3)then
		
					if (n.eq.0)then
					  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
					  write(63,*)"output2",t
					  close(63)
					  end if
					if (fastmovie.eq.1)then

					call movie

					else
					call outwrite3vb

					end if
					
					
						if (n.eq.0)then
					  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
					  write(63,*)"output3",t
					  close(63)
					  end if
				else
				      if (n.eq.0)then
					  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
					  write(63,*)"output1",t
					  close(63)
					  end if


				call outwrite3vb2d
				
					if (n.eq.0)then
					  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
					  write(63,*)"output3",t
					  close(63)
					  end if
				end if
	end if
	if (tecplot.eq.0)then
  			if (dimensiona.eq.3)then

				call outwrite3v
		      else

 				call outwrite3v2d
		      end if

	end if
	
	if (tecplot.eq.2)then		!binary paraview 3d only
	if (dimensiona.eq.3)then




        if (fastmovie.eq.1)then

					call movie_para

					else

					call outwritepara3db

					end if
    else
    
        call outwritepara2db
    end if

end if

if (tecplot.eq.3)then		!binary paraview 3d only

  call outwritepara3dbp

end if

if (tecplot.eq.4)then		!binary paraview 3d only

  call outwritetec3dbp
  call mpi_barrier(mpi_comm_world,ierror)
end if

if (tecplot.eq.5)then		!fast output written by all processors using mpi-io

call parallel_vtk_combine(n)

end if

if (tecplot.eq.6)then		!fast output written by all processors

call parallel_vtk_combine_partitioned(n)

end if




					if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished writing output",t
				  close(63)
				  end if





end subroutine volume_solution_write



subroutine surface_solution_write
!> @brief
!> this subroutine calls the appropriate surface solution writing subroutine based on the settings
implicit none



#ifdef gpu
call gpu_to_host_surf
#endif


if (tecplot.eq.1)then
    if (dimensiona.eq.3)then

  call outwrite3vsb
  else

  call outwrite3vsb2d
  end if
end if
if (tecplot.eq.0)then
      if (dimensiona.eq.3)then

    call outwrite3vs
    else

    call outwrite3vs2d
    end if

end if
if (tecplot.eq.2)then		!binary paraview 3d only
    if (dimensiona.eq.3)then

  call outwritepara3dsb
    else
  !call outwritepara2dsb    !not implemented yet
    end if

end if


if (tecplot.eq.3)then		!binary paraview 3d only
    if (dimensiona.eq.3)then

  call outwritepara3dsb
    else
  !call outwritepara2dsb    !not implemented yet
    end if

end if


if (tecplot.eq.4)then
    if (dimensiona.eq.3)then

  call outwrite3vsb
  else

  call outwrite3vsb2d
  end if
end if

if (tecplot.eq.5)then		!binary paraview 3d only

  call parallel_vtk_combine_wall(n)

end if


if (tecplot.eq.6)then		!fast output written by all processors



call parallel_vtk_combine_partitioned_wall(n)

end if


end subroutine surface_solution_write

subroutine volume_solution_write_av
!> @brief
!> this subroutine calls the appropriate average volume writing subroutine based on the settings
implicit none

#ifdef gpu
call gpu_to_host_vol_av
#endif



if (tecplot.eq.1)then


call outwrite3vbav

end if

if (tecplot.eq.0)then
    if (dimensiona.eq.3)then

  call outwrite3vav
  else

  call outwrite3v2dav
  end if

end if


if (tecplot.eq.2)then

call outwritepara3dbav

end if


if (tecplot.eq.3)then

call outwritepara3dbpav

end if


if (tecplot.eq.4)then		!binary paraview 3d only

  call outwritetec3dbpav
  
end if

if (tecplot.eq.5)then
	if (dimensiona.eq.3)then
	call parallel_vtk_combine_av(n)
	end if
end if

if (tecplot.eq.6)then
	if (dimensiona.eq.3)then
	call parallel_vtk_combine_partitioned_av(n)
	end if
end if




end subroutine volume_solution_write_av



subroutine surface_solution_write_av
!> @brief
!> this subroutine calls the appropriate surface writing subroutine based on the settings
implicit none


#ifdef gpu
call gpu_to_host_surf_av
#endif

if (tecplot.eq.1)then
   if (dimensiona.eq.3)then

call outwrite3vsbav
else

call outwrite3vsb2dav
end if
end if

if (tecplot.eq.0)then
  if (dimensiona.eq.3)then

call outwrite3vsav
else

call outwrite3vs2dav
end if

end if
if (tecplot.eq.2)then

call outwritepara3dsbav

end if

if (tecplot.eq.3)then

call outwritepara3dsbav

end if



if (tecplot.eq.4)then
   if (dimensiona.eq.3)then

call outwrite3vsbav
else

call outwrite3vsb2dav
end if
end if

if (tecplot.eq.5)then
	if (dimensiona.eq.3)then
	call parallel_vtk_combine_wall_av(n)
	end if
end if



if (tecplot.eq.6)then
	if (dimensiona.eq.3)then
	call parallel_vtk_combine_partitioned_wall_av(n)
	end if
end if





end subroutine surface_solution_write_av

subroutine forces
!> @brief
!> this subroutine calls the appropriate force computation subroutine based on the dimensionality of the problem
implicit none


   if (dimensiona.eq.3)then

call computeforce(n)
else

call computeforce2d(n)
end if



end subroutine forces



subroutine residual_compute
!> @brief
!> this subroutine calls the appropriate residual computation subroutine based on the dimensionality of the problem
implicit none


   if (dimensiona.eq.3)then

call calculate_residual(n)
else

call calculate_residual2d(n)
end if

end subroutine residual_compute




subroutine checkpointing
!> @brief
!> this subroutine calls the appropriate checkpointing subroutine based on the dimensionality of the problem
implicit none


#ifdef gpu
				if (allocated(u_c_val)) then
                  !$omp target update from(u_c_val(1,1:nof_variables,1:xmpielrank(n)))
                  end if
				 if (allocated(u_ct_val)) then
                  !$omp target update from(u_ct_val(ind1,1,1:xmpielrank(n)))
                  end if
#endif


   if (dimensiona.eq.3)then

call checkpoint(n)
else

call checkpoint2d(n)
end if

end subroutine checkpointing


subroutine checkpointing_av
!> @brief
!> this subroutine calls the appropriate averaged checkpointing subroutine based on the dimensionality of the problem
implicit none


#ifdef gpu
call gpu_to_host_vol_av
#endif


   if (dimensiona.eq.3)then

call checkpointav(n)
else

call checkpointav2d(n)
end if

end subroutine checkpointing_av



subroutine writethem(n)
integer,intent(in)::n
integer::i,l,ngp,k,kmaxe
character(len=20)::proc,outfile,proc3,surfile,proc4


write(proc3,fmt='(i10)') n



outfile="OUTX_"//trim(adjustl(proc3))//".dat"


open(140,file=outfile,form='formatted',status='new',action='write')

kmaxe=xmpielrank(n)
do i=1,nof_interior

		!		write(140,*) u_c_val(1,1:5,i)

! 		do l=1,ielem_ifca(i)
!
! 				!write(140,*)ielem_ihexgl(i),l
 				write(140,*)rec_uleft(1,1,1, i)
!
!
! 		end do
end do


close(140)


end subroutine







subroutine checkpoint(n)
!> @brief
!> this subroutine uses mpi-io for writing the checkpointing files
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::dispt
real,allocatable,dimension(:)::array2
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,size_of_real,size_of_int,dip,n_end,datatype,ifg
character(len=20)::proc,restfile,proc3
real,allocatable,dimension(:)::array
logical::here1
real::in1,iocpt1,iocpt2,iocpt3,iocpt4
integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init
disp_in_file=0
tmp=0
disp_init=0

kmaxe=xmpielrank(n)

size_of_int=4
size_of_real=8
icpuid=n
call mpi_barrier(mpi_comm_world,ierror)
!  iocpt1=mpi_wtime()


if (dg.eq.1)then
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+turbulenceequations+passivescalar)*(idegfree+1)))
else
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+turbulenceequations+passivescalar)))	
end if

     if (dg.eq.1)then
     do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*((nof_variables+turbulenceequations+passivescalar)*(idegfree+1))
      end do
      

      n_end=(nof_variables+turbulenceequations+passivescalar)*(idegfree+1)
     
     else

      do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*((nof_variables+turbulenceequations+passivescalar))
      end do
      

      n_end=nof_variables+turbulenceequations+passivescalar
      end if

      
      if (dg.eq.1)then
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	    k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(1,1:nof_variables,i)
	      k=k+nof_variables
	      array2(k:k+turbulenceequations+passivescalar-1)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
	      k=k+turbulenceequations+passivescalar
	  end do
      else
	  k=1
	  do i=1,kmaxe
        do j=1,nof_variables
	      array2(k:k+idegfree)=u_c_valdg(1,j,1:idegfree+1,i)
	      k=k+(idegfree+1)
        end do
	  end do
      end if
      
      else
      
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	    k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(1,1:nof_variables,i)
	      k=k+nof_variables
	      array2(k:k+turbulenceequations+passivescalar-1)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
	      k=k+turbulenceequations+passivescalar
	  end do
      else
	  k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(1,1:nof_variables,i)
	      k=k+nof_variables
	  end do
      end if
      end if
      
	restfile='RESTART.dat'

       if (n.eq.0)then
	 inquire (file=restfile,exist=here1)
	  call mpi_file_delete(restfile,mpi_info_null,ierror)
	end if
     
    call mpi_barrier(mpi_comm_world,ierror)

    !create type first of indexed block
    
    call mpi_type_create_indexed_block(kmaxe,n_end,dispt,mpi_double_precision,datatype,ierror)
    call mpi_type_commit(datatype,ierror)
   

    allocate(array(1:nof_variables+turbulenceequations+passivescalar))
    
    
	
   
	
	
	
	!call mpi_barrier(mpi_comm_world,ierror)

	call mpi_file_open(mpi_comm_world, restfile,mpi_mode_wronly + mpi_mode_create,mpi_info_null, fh, ierror)
	
	
	 
	  
	
	if (n.eq.0)then
	
	    if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	          call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_write(fh, it, 1, mpi_integer, mpi_status_ignore,ierror)
		    disp_in_file = disp_in_file + size_of_int	!1
		    call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_write(fh, initialres(1:nof_variables+turbulenceequations),nof_variables+turbulenceequations, mpi_double_precision, mpi_status_ignore, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_variables+turbulenceequations)	!3
	    else
		  call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  
		  call mpi_file_write(fh, it, 1, mpi_integer, mpi_status_ignore, ierror)
		    disp_in_file = disp_in_file + size_of_int 	!4
		    
		  call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  
		  call mpi_file_write(fh, t, 1, mpi_double_precision, mpi_status_ignore,ierror)
		    disp_in_file = disp_in_file + size_of_real    !5
		      if (initcond.eq.95)then
		      call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
		      call mpi_file_write(fh, taylor, 1, mpi_double_precision, mpi_status_ignore, ierror)
			    disp_in_file = disp_in_file + size_of_real !6
		      end if  
	    end if
	else
	      if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
		disp_in_file = disp_in_file + size_of_int 
		disp_in_file = disp_in_file + size_of_real*(nof_variables+turbulenceequations)
	      else
		  disp_in_file = disp_in_file + size_of_int 
		  disp_in_file = disp_in_file + size_of_real
		  if (initcond.eq.95)then
		  disp_in_file = disp_in_file + size_of_real 
		  end if
	      end if
	      
	      
	end if
	
	
	call mpi_barrier(mpi_comm_world, ierror)
	call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatype, 'native',mpi_info_null, ierror)
	call mpi_file_write_all(fh, array2, kmaxe*n_end, mpi_double_precision,mpi_status_ignore, ierror)        
        
        call mpi_file_close(fh, ierror)
	call mpi_type_free(datatype,ierror)
          

	
	
	
	
	deallocate(array,dispt,array2)
	call mpi_barrier(mpi_comm_world, ierror)
	
	

	
	


end subroutine checkpoint



subroutine prepare_surfaces_v(n)
	implicit none
	integer,intent(in)::n
	integer::i,j,k,kmaxe,temp_cord,iloop,kloop,zloop,kloopf,temp_loop,kn
	integer,allocatable,dimension(:)::list_nin2,list_nout2
	kmaxe=xmpielrank(n)
iloop=0
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		  iloop=iloop+1
	      end if
	  end if
	end do
   end if
end do


iloopx=iloop

if (iloopx.gt.0)then
allocate(wall_l(1:iloopx,1:4))

else
allocate(wall_l(0,0))
end if

iloop=0
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	     if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		  iloop=iloop+1
				if (dimensiona.eq.2)then
					kloop=2
				else
					if (ielem_types_faces(j,i).eq.5)then
					kloop=4
					else
					kloop=3
					end if
				end if

		  wall_l(iloop,1)=i		!this contains the iloop of all the elements that are needed with local numbering
		  wall_l(iloop,2)=j		!face
		  wall_l(iloop,3)=kloop !number of nodes
	      end if
	  end if
	end do
   end if
end do


	iwmaxe=0


	!now find all the element types and copy them in a local list

	call mpi_allreduce(iloopx,iwmaxe,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)




	allocate(wallshape_g(1:iwmaxe));wallshape_g=0

	if (iloopx.gt.0)then
	allocate(wallshape(1:iloopx));wallshape=0

	do i=1,iloopx
		if (wall_l(i,3).eq.2)then	!line
			wallshape(i)=3
		end if

		if (wall_l(i,3).eq.3)then	!triangular
			wallshape(i)=5
		end if

		if (wall_l(i,3).eq.4)then	!quadrilateral
			wallshape(i)=9
		end if
	end do


	else
	allocate(wallshape(0:0));wallshape=0
	end if



	!we need the total with an offset
	allocate(wallcx_g(0:isize-1),offsetwc_g(0:isize-1))
wallcx_g(:)=0
offsetwc_g(0:isize-1)=0

wallcx_g(n)=iloopx

call mpi_allgather(iloopx,1,mpi_integer,wallcx_g,1,mpi_integer,mpi_comm_world,ierror)

offsetwc_g(0)=0
do i=1,isize-1
		offsetwc_g(i)=offsetwc_g(i-1)+wallcx_g(i-1)
end do



if (iloopx.gt.0)then
wallshape_g(offsetwc_g(n)+1:offsetwc_g(n)+wallcx_g(n))=wallshape(1:iloopx)
end if










allocate(wallshape_g2(1:iwmaxe));wallshape_g2=0




call mpi_allreduce(wallshape_g,wallshape_g2,iwmaxe,mpi_integer,mpi_max,mpi_comm_world,ierror)








	if (iloopx.gt.0)then
	allocate(typ_nodesn_w(1:iloopx));typ_nodesn_w(:)=0
	else
	allocate(typ_nodesn_w(0:0));typ_nodesn_w(:)=0
	end if
	typ_countn_w=0
	if (iloopx.gt.0)then
	temp_loop=offsetwc_g(n)
	do i=1,iloopx
		typ_countn_w=typ_countn_w+wall_l(i,3)
		temp_loop=temp_loop+1
		wall_l(i,4)=temp_loop

		typ_nodesn_w(i)=wall_l(i,3)
	end do
	end if


allocate(nodes_offsetw(1:iwmaxe),nodes_offsetw2(1:iwmaxe))






call mpi_allreduce(typ_countn_w,typ_countn_global_w,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)

kn=0
do j=1,iwmaxe
	nodes_offsetw(j)=kn
		if (wallshape_g2(j).eq.3)then
			kn=kn+2
		end if
		if (wallshape_g2(j).eq.5)then
			kn=kn+3
		end if
		if (wallshape_g2(j).eq.9)then
			kn=kn+4
		end if
	nodes_offsetw2(j)=kn


end do


	if (iloopx.gt.0)then
	allocate(nodes_offset_localw(1:iloopx));nodes_offset_localw=0
	allocate(nodes_offset_local2w(1:iloopx));nodes_offset_local2w=0




	else
	allocate(nodes_offset_localw(0:0));nodes_offset_localw=0
	allocate(nodes_offset_local2w(0:0));nodes_offset_local2w=0
	end if

	if (iloopx.gt.0)then
	do i=1,iloopx
	nodes_offset_localw(i)=nodes_offsetw(wall_l(i,4))
	nodes_offset_local2w(i)=nodes_offsetw2(wall_l(i,4))

	end do
	end if




end subroutine prepare_surfaces_v




















subroutine partition_preparation_wallv(n)
	implicit none
	integer,intent(in)::n
	integer::i,j,k,temp_cord,kloop,varg_max,kmaxn_p
	integer,dimension(4)::tempwnode
	kmaxn_p=xmpiall_v(n)
	temp_cord=3

	if (dimensiona.eq.2)then
	wnodes_part=2
	else
	wnodes_part=3
	end if
	varg_max=max(write_variables_w,write_variables_av_w)

	if (iloopx.gt.0)then
 	allocate(wdispart1(1:iloopx),wrarray_part1(1:iloopx,1:varg_max))	!

 		do i=1,iloopx
 			wdispart1(i)=(wall_l(i,4)-1)*1
 		end do
			wpart1_end=1
	else
		allocate(wdispart1(1),wrarray_part1(1,1:varg_max))	!
		wrarray_part1(1,:)=0
		wdispart1(1)=0!
			wpart1_end=0

	end if





	if (iloopx.gt.0)then
 	allocate(wdispart2(1:iloopx),wiarray_part2(1:typ_countn_w))		!
 		do i=1,iloopx
 			wdispart2(i)=nodes_offset_localw(i)
 		end do

 		wpart2_end=wnodes_part


	else
		allocate(wdispart2(1),wiarray_part2(1))	!
		wiarray_part2(1)=0
		wdispart2(1)=0!(totwallsc-1)*wnodes_part	!maybe total number of values written
		wpart2_end=0

	end if



	!the nodes of the boundary cells

	if (iloopx.gt.0)then




 		k=1
 	  do i=1,iloopx

		wiarray_part2(k:k+typ_nodesn_w(i)-1)=ielem_nodes_faces_v(wall_l(i,2),1:typ_nodesn_w(i),wall_l(i,1))
		k=k+typ_nodesn_w(i)
 	  end do


    end if




		if (iloopx.gt.0)then
		allocate(wdispart5(1:iloopx),wiarray_part5(1:iloopx))


			do i=1,iloopx
				wdispart5(i)=(wall_l(i,4)-1)*1
			end do

			do i=1,iloopx
				wiarray_part5(i)=nodes_offset_local2w(i)
			end do


		else
			allocate(wdispart5(1),wiarray_part5(1))
				wdispart5(1)=0!totwallsc-1
				wiarray_part5(1)=0


		end if



	if (iloopx.gt.0)then

!	!the types of elements of the nodes of  the boundary cells
	allocate(wdispart3(1:iloopx),wiarray_part3(1:iloopx))

			do i=1,iloopx
				wdispart3(i)=(wall_l(i,4)-1)*1
				wiarray_part3(i)=wallshape(i)
			end do
			wpart3_end=1
	else
	allocate(wdispart3(1),wiarray_part3(1))
			wdispart3(1)=0!totwallsc-1
			wiarray_part3(1)=0
			wpart3_end=0
	end if





  	allocate(wdispart4(1:kmaxn_p),wrarray_part4(1:kmaxn_p*temp_cord))		!
  	wpart4_end=temp_cord

  	do i=1,kmaxn_p
			wdispart4(i)=(my_nodesg(i)-1)*(temp_cord)
		end do

		k=1
	  do i=1,kmaxn_p
		wrarray_part4(k:k+dims-1)=inoder4_cord(1:dims,my_nodesl(i))
		if (dimensiona.eq.2)then
		wrarray_part4(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do



!
!
if (n.eq.0)then
wkdum1=1
wkdum2=1
wkdum3(1)=0
else
wkdum1=0
wkdum2=0
wkdum3(1)=1
end if
!
			!now commit datatypes

 				call mpi_type_create_indexed_block(iloopx,wpart1_end,wdispart1,mpi_double_precision,wdatatypex,ierror)
 				call mpi_type_commit(wdatatypex,ierror)
!
 				!dummy type for writing one component only from one cpu
 				call mpi_type_create_indexed_block(wkdum1,wkdum2,wkdum3,mpi_integer,wdatatypeint,ierror)
 				call mpi_type_commit(wdatatypeint,ierror)
!
!
!
! 				!point coordinates
 				call mpi_type_create_indexed_block(kmaxn_p,wpart4_end,wdispart4,mpi_double_precision,wdatatypez,ierror)
 				call mpi_type_commit(wdatatypez,ierror)
!
! 				!connectivity
!  				call mpi_type_create_indexed_block(iloopx,typ_nodesn_w(:),wdispart2,mpi_integer,wdatatypey,ierror)
 				call mpi_type_indexed(iloopx,typ_nodesn_w(:),wdispart2,mpi_integer,wdatatypey,ierror)
 				call mpi_type_commit(wdatatypey,ierror)
!
! 				!type of element
 				call mpi_type_create_indexed_block(iloopx,wpart1_end,wdispart5,mpi_integer,wdatatypexx,ierror)
 				call mpi_type_commit(wdatatypexx,ierror)
!
! 				!nodes
 				call mpi_type_create_indexed_block(iloopx,wpart3_end,wdispart3,mpi_integer,wdatatypeyy,ierror)
 				call mpi_type_commit(wdatatypeyy,ierror)
!
!
!


call mpi_barrier(mpi_comm_world,ierror)


deallocate(wallshape_g,wallshape,wallcx_g,offsetwc_g,wallshape_g2,typ_nodesn_w,nodes_offsetw,nodes_offsetw2,nodes_offset_localw,nodes_offset_local2w)





end subroutine partition_preparation_wallv





subroutine partition_preparation(n)
	implicit none
	integer,intent(in)::n
	integer::i,j,k,kmaxe,kmaxn_p,temp_cord,varg_max
	kmaxe=xmpielrank(n)
	kmaxn_p=xmpiall_v(n)

	temp_cord=3

	if (dimensiona.eq.2)then
	nodes_part=4
	else
	nodes_part=8
	end if
	varg_max=max(write_variables,write_variables_av)
	allocate(dispart1(1:kmaxe),rarray_part1(1:kmaxe,1:varg_max))	!

		do i=1,kmaxe
			dispart1(i)=(xgo(i)-1)*(1)
		end do
		part1_end=1


	allocate(typ_nodesn(1:kmaxe));typ_nodesn(:)=0
	typ_countn=0
	do i=1,kmaxe
			if (ielem_ishape(i).eq.1)then	!hexa
				typ_countn=typ_countn+8
				typ_nodesn(i)=8

			end if
			if (ielem_ishape(i).eq.2)then	!tetra
				typ_countn=typ_countn+4
				typ_nodesn(i)=4
! 				typ_countn=typ_countn+8
! 				typ_nodesn(i)=8

			end if
			if (ielem_ishape(i).eq.3)then	!pyramid
				typ_countn=typ_countn+5
				typ_nodesn(i)=5
! 				typ_countn=typ_countn+8
! 				typ_nodesn(i)=8
			end if
			if (ielem_ishape(i).eq.4)then	!prism
				typ_countn=typ_countn+6
				typ_nodesn(i)=6
! 				typ_countn=typ_countn+8
! 				typ_nodesn(i)=8
			end if
			if (ielem_ishape(i).eq.5)then	!quad
				typ_countn=typ_countn+4
				typ_nodesn(i)=4
			end if
			if (ielem_ishape(i).eq.6)then	!triangular
				typ_countn=typ_countn+3
				typ_nodesn(i)=3
! 				typ_countn=typ_countn+4
! 				typ_nodesn(i)=4
			end if
	end do




	call mpi_allreduce(typ_countn,typ_countn_global,1,mpi_integer,mpi_sum,mpi_comm_world,ierror)


! 		print*,typ_countn,typ_countn_global,"check here",imaxe*6




	allocate(dispart2(1:kmaxe),iarray_part2(1:typ_countn))		!
		do i=1,kmaxe
! 			dispart2(i)=(xgo(i)-1)*(nodes_part)
!  			dispart2(i)=(xgo(i)-1)*(typ_nodesn(i))
 			dispart2(i)=nodes_offset_local(i)

		end do

		part2_end=nodes_part

			k=1
	  do i=1,kmaxe
! 		iarray_part2(k:k+nodes_part-1)=ielem_nodes_v(1:nodes_part,i)
! 	      k=k+nodes_part
	      iarray_part2(k:k+typ_nodesn(i)-1)=ielem_nodes_v(1:typ_nodesn(i),i)
	      k=k+typ_nodesn(i)
	  end do



	  allocate(dispart5(1:kmaxe),iarray_part5(1:kmaxe))
		do i=1,kmaxe
			dispart5(i)=(xgo(i)-1)*(1)
		end do


	   do i=1,kmaxe
		iarray_part5(i)=nodes_offset_local2(i)!(xgo(i))*typ_nodesn(i)
	  end do




	allocate(dispart3(1:kmaxe),iarray_part3(1:kmaxe))
 		do i=1,kmaxe
 			dispart3(i)=(xgo(i)-1)*(1)


 			if (ielem_ishape(i).eq.1)then	!hexa
				iarray_part3(i)=12
			end if
			if (ielem_ishape(i).eq.2)then	!tetra
				iarray_part3(i)=10
			end if
			if (ielem_ishape(i).eq.3)then	!pyramid
				iarray_part3(i)=14
			end if
			if (ielem_ishape(i).eq.4)then	!prism
				iarray_part3(i)=13
			end if
			if (ielem_ishape(i).eq.5)then	!quad
				iarray_part3(i)=9
			end if
			if (ielem_ishape(i).eq.6)then	!triangular
				iarray_part3(i)=5
			end if





 		end do

		part3_end=1





	allocate(dispart4(1:kmaxn_p),rarray_part4(1:kmaxn_p*temp_cord))

		part4_end=temp_cord
		do i=1,kmaxn_p
			dispart4(i)=(my_nodesg(i)-1)*(temp_cord)
		end do


		k=1
	  do i=1,kmaxn_p
		rarray_part4(k:k+dims-1)=inoder4_cord(1:dims,my_nodesl(i))
		if (dimensiona.eq.2)then
		rarray_part4(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do


if (n.eq.0)then
kdum1=1
kdum2=1
kdum3(1)=0
else
kdum1=0
kdum2=0
kdum3(1)=1
end if

 !now commit datatypes

				call mpi_type_create_indexed_block(kmaxe,part1_end,dispart1,mpi_double_precision,datatypex,ierror)
				call mpi_type_commit(datatypex,ierror)

				!dummy type for writing one component only from one cpu
				call mpi_type_create_indexed_block(kdum1,kdum2,kdum3,mpi_integer,datatypeint,ierror)
				call mpi_type_commit(datatypeint,ierror)



				!point coordinates
				call mpi_type_create_indexed_block(kmaxn_p,part4_end,dispart4,mpi_double_precision,datatypez,ierror)
				call mpi_type_commit(datatypez,ierror)

				!connectivity
! 				call mpi_type_create_indexed_block(kmaxe,nodes_part,dispart2,mpi_integer,datatypey,ierror)
				call mpi_type_indexed(kmaxe,typ_nodesn(1:kmaxe),dispart2,mpi_integer,datatypey,ierror)
				call mpi_type_commit(datatypey,ierror)

				!type of element
				call mpi_type_create_indexed_block(kmaxe,part1_end,dispart5,mpi_integer,datatypexx,ierror)
				call mpi_type_commit(datatypexx,ierror)

				!nodes
				call mpi_type_create_indexed_block(kmaxe,part1_end,dispart3,mpi_integer,datatypeyy,ierror)
				call mpi_type_commit(datatypeyy,ierror)














end subroutine partition_preparation

subroutine partition_preparation_p(n)
	implicit none
	integer,intent(in)::n
	integer::i,j,k,kmaxe,kmaxn_p,temp_cord,varg_max
	integer,dimension(8)::temp_part_n
	kmaxe=xmpielrank(n)


	temp_cord=3


	varg_max=max(write_variables,write_variables_av)
	allocate(sol_vtu(1:kmaxe,1:varg_max));sol_vtu=0	!

	allocate(typ_nodesn(1:kmaxe));typ_nodesn(:)=0
	typ_countn=0
	do i=1,kmaxe
			if (ielem_ishape(i).eq.1)then	!hexa
				typ_countn=typ_countn+8
				typ_nodesn(i)=8

			end if
			if (ielem_ishape(i).eq.2)then	!tetra
				typ_countn=typ_countn+4
				typ_nodesn(i)=4


			end if
			if (ielem_ishape(i).eq.3)then	!pyramid
				typ_countn=typ_countn+5
				typ_nodesn(i)=5

			end if
			if (ielem_ishape(i).eq.4)then	!prism
				typ_countn=typ_countn+6
				typ_nodesn(i)=6

			end if
			if (ielem_ishape(i).eq.5)then	!quad
				typ_countn=typ_countn+4
				typ_nodesn(i)=4
			end if
			if (ielem_ishape(i).eq.6)then	!triangular
				typ_countn=typ_countn+3
				typ_nodesn(i)=3

			end if
	end do






	allocate(offset_vtu(1:kmaxe),connect_vtu(1:typ_countn))
	!
	k=0
		do i=1,kmaxe
				k=k+typ_nodesn(i)
			offset_vtu(i)=k
		end do



			k=1
	  do i=1,kmaxe
			temp_part_n(1:typ_nodesn(i))=ielem_nodes(1:typ_nodesn(i),i)

			do j=1,typ_nodesn(i)
				temp_part_n(j)=temp_part_n(j)-1
			end do

	      connect_vtu(k:k+typ_nodesn(i)-1)=temp_part_n(1:typ_nodesn(i))
	      k=k+typ_nodesn(i)
	  end do





	allocate(type_vtu(1:kmaxe))
 		do i=1,kmaxe
 			if (ielem_ishape(i).eq.1)then	!hexa
				type_vtu(i)=12
			end if
			if (ielem_ishape(i).eq.2)then	!tetra
				type_vtu(i)=10
			end if
			if (ielem_ishape(i).eq.3)then	!pyramid
				type_vtu(i)=14
			end if
			if (ielem_ishape(i).eq.4)then	!prism
				type_vtu(i)=13
			end if
			if (ielem_ishape(i).eq.5)then	!quad
				type_vtu(i)=9
			end if
			if (ielem_ishape(i).eq.6)then	!triangular
				type_vtu(i)=5
			end if


 		end do






	allocate(nodes_vtu(1:kmaxn*temp_cord))





		k=1
	  do i=1,kmaxn
		nodes_vtu(k:k+dims-1)=inoder4_cord(1:dims,i)
		if (dimensiona.eq.2)then
		nodes_vtu(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do




end subroutine partition_preparation_p


subroutine partition_preparation_p_wall(n)
	implicit none
	integer,intent(in)::n
	integer::i,j,k,kmaxe,kmaxn_p,temp_cord,varg_max,iloop,kloop,zloop,kloopf,temp_loop,kn,kkd
	integer,dimension(8)::temp_part_n
	kmaxe=xmpielrank(n)


	temp_cord=3
				iloop=0
			do i=1,kmaxe
			if (ielem_interior(i).eq.1)then
				do j=1,ielem_ifca(i)
				if (ielem_ibounds(j,i).gt.0)then
					if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
					iloop=iloop+1
					end if
				end if
				end do
			end if
			end do



			iloopx=iloop


	allocate(wallcount_cpu_l(0:isize-1),wallcount_cpu_g(0:isize-1));wallcount_cpu_l=0;wallcount_cpu_g=0


	wallcount_cpu_l(n)=iloopx

	call mpi_allreduce(wallcount_cpu_l(0:isize-1),wallcount_cpu_g(0:isize-1),isize,mpi_integer,mpi_max,mpi_comm_world,ierror)


	call mpi_barrier(mpi_comm_world,ierror)


if (iloopx.gt.0)then
allocate(wall_l(1:iloopx,1:4))

else
allocate(wall_l(0,0))
end if

iloop=0
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	    if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		  iloop=iloop+1
				if (dimensiona.eq.2)then
					kloop=2
				else
					if (ielem_types_faces(j,i).eq.5)then
					kloop=4
					else
					kloop=3
					end if
				end if

		  wall_l(iloop,1)=i		!this contains the iloop of all the elements that are needed with local numbering
		  wall_l(iloop,2)=j		!face
		  wall_l(iloop,3)=kloop !number of nodes
	      end if
	  end if
	end do
   end if
end do



	varg_max=max(write_variables_w,write_variables_av_w)
	if (iloopx.gt.0)then
	allocate(sol_vtu_w(1:iloopx,1:varg_max));sol_vtu=0	!
	else
	allocate(sol_vtu_w(0:0,1:varg_max));sol_vtu=0	!
	end if

	if (iloopx.gt.0)then
	allocate(typ_nodesn_w(1:iloopx),type_vtu_w(1:iloopx));typ_nodesn_w(:)=0;type_vtu_w=0
	else
	allocate(typ_nodesn_w(0:0),type_vtu_w(0:0));typ_nodesn_w(:)=0;type_vtu_w=0
	end if



	typ_countn_w=0
	if (iloopx.gt.0)then



	do i=1,iloopx
			if (wall_l(i,3).eq.4)then	!quad
				typ_countn_w=typ_countn+4
				typ_nodesn_w(i)=4
				type_vtu_w(i)=9
			end if

			if (wall_l(i,3).eq.3)then	!tri
				typ_countn_w=typ_countn+3
				typ_nodesn_w(i)=3
				type_vtu_w(i)=5
			end if

			if (wall_l(i,3).eq.2)then	!line
				typ_countn_w=typ_countn+2
				typ_nodesn_w(i)=2
				type_vtu_w(i)=3
			end if
	end do

	end if



	if (iloopx.gt.0)then
	allocate(offset_vtu_w(1:iloopx),connect_vtu_w(1:typ_countn_w))
	else
	allocate(offset_vtu_w(0:0),connect_vtu_w(0:0))
	end if
	!

	if (iloopx.gt.0)then
	k=0
		do i=1,iloopx
				k=k+typ_nodesn_w(i)
			offset_vtu_w(i)=k
		end do



			k=1
	  do i=1,iloopx
			temp_part_n(1:typ_nodesn_w(i))=ielem_nodes_faces(wall_l(i,2),1:typ_nodesn_w(i),wall_l(i,1))

			do j=1,typ_nodesn_w(i)
				temp_part_n(j)=temp_part_n(j)-1
			end do

	      connect_vtu_w(k:k+typ_nodesn_w(i)-1)=temp_part_n(1:typ_nodesn_w(i))
	      k=k+typ_nodesn_w(i)
	  end do


	end if








	if (iloopx.gt.0)then
	allocate(nodes_vtu_w(1:kmaxn*temp_cord))

	else
	allocate(nodes_vtu_w(0:0))
	end if


	if (iloopx.gt.0)then

		k=1
	  do i=1,kmaxn
		nodes_vtu_w(k:k+dims-1)=inoder4_cord(1:dims,i)
		if (dimensiona.eq.2)then
		nodes_vtu_w(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do
	end if



end subroutine partition_preparation_p_wall



subroutine specify_write_variables(n)
!> @brief
!> this subroutine uses mpi-io for writing the vtk files
implicit none
integer,intent(in)::n
integer::i,j,k




!volume instantaneous variables

if (dimensiona.eq.3)then



			if (multispecies.eq.1)then
			write_variables=nof_variables+1
			!!specify the name of the variable names!!

			variable_names(1)='density'
			variable_names(2)='u'
			variable_names(3)='v'
			variable_names(4)='w'
			variable_names(5)='pressure'
			variable_names(6)='rho vf1'
			variable_names(7)='rho vf2'
			variable_names(8)='volume_fraction'
			variable_names(9)='q'


			else

			write_variables=nof_variables+1+turbulenceequations+adda


			!!specify the name of the variable names!!

			variable_names(1)='density'
			variable_names(2)='u'
			variable_names(3)='v'
			variable_names(4)='w'
			variable_names(5)='pressure'
			variable_names(6)='q'

				if (adda.eq.1)then
				variable_names(nof_variables+1+adda)='adda'
				end if

				if (turbulence.eq.1)then
				variable_names(nof_variables+1+turbulenceequations+adda)='turb'
				end if

				if (itestcase.eq.1)then
				variable_names(1)='solution'
				variable_names(2)='aux'
				end if



			end if

			if (realgas.eq.1)then
			write_variables=nof_variables+2+turbulenceequations
			variable_names(1)='density'
			variable_names(2)='u'
			variable_names(3)='v'
			variable_names(4)='w'
			variable_names(5)='ttr'
			variable_names(6)='tvb'
			variable_names(7)='n2'
			variable_names(8)='o2'
			variable_names(9)='no'
			variable_names(10)='n'
			variable_names(11)='o'
			variable_names(12)='pressure'
			variable_names(13)='q'
			if (turbulence.eq.1)then
				variable_names(nof_variables+2+turbulenceequations)='turb'
			end if
			end if














else




		if (multispecies.eq.1)then
		write_variables=nof_variables+1
		!!specify the name of the variable names!!

		variable_names(1)='density'
		variable_names(2)='u'
		variable_names(3)='v'
		variable_names(4)='pressure'
		variable_names(5)='rho vf1'
		variable_names(6)='rho vf2'
		variable_names(7)='volume_fraction'
		variable_names(8)='q'


		else

		write_variables=nof_variables+1+turbulenceequations
		!!specify the name of the variable names!!

		variable_names(1)='density'
		variable_names(2)='u'
		variable_names(3)='v'
		variable_names(4)='pressure'
		variable_names(5)='q'
		if (turbulence.eq.1)then
		variable_names(6)='turb'
		end if



		if (itestcase.eq.1)then
		variable_names(1)='solution'
		variable_names(2)='aux'
		end if



		end if



		if (realgas.eq.1)then
		write_variables=nof_variables+2+turbulenceequations
		variable_names(1)='density'
			variable_names(2)='u'
			variable_names(3)='v'
			variable_names(4)='ttr'
			variable_names(5)='tvb'
			variable_names(6)='n2'
			variable_names(7)='o2'
			variable_names(8)='no'
			variable_names(9)='n'
			variable_names(10)='o'
			variable_names(11)='pressure'
			variable_names(12)='q'
		if (turbulence.eq.1)then
		variable_names(nof_variables+2+turbulenceequations)='turb'
		end if
		end if


end if





if(averaging.eq.1)then

if (dimensiona.eq.3)then




write_variables_av=11


!!specify the name of the variable names!!

variable_names_av(1)='r_mean'
variable_names_av(2)='u_mean'
variable_names_av(3)='v_mean'
variable_names_av(4)='w_mean'
variable_names_av(5)='p_mean'
variable_names_av(6)='u_rms'
variable_names_av(7)='v_rms'
variable_names_av(8)='w_rms'
variable_names_av(9)='uv'
variable_names_av(10)='uw'
variable_names_av(11)='wv'

if (realgas.eq.1)then
write_variables_av=18



variable_names_av(1)='r_mean'
variable_names_av(2)='u_mean'
variable_names_av(3)='v_mean'
variable_names_av(4)='w_mean'
variable_names_av(5)='ttr_mean'
variable_names_av(6)='tvb_mean'
variable_names_av(7)='n2_mean'
variable_names_av(8)='o2_mean'
variable_names_av(9)='no_mean'
variable_names_av(10)='n_mean'
variable_names_av(11)='o_mean'
variable_names_av(12)='p_mean'
variable_names_av(13)='u_rms'
variable_names_av(14)='v_rms'
variable_names_av(15)='w_rms'
variable_names_av(16)='uv'
variable_names_av(17)='uw'
variable_names_av(18)='wv'






end if



write_variables_av_w=8

variable_names_av_w(1)='r_mean'
variable_names_av_w(2)='u_mean'
variable_names_av_w(3)='v_mean'
variable_names_av_w(4)='w_mean'
variable_names_av_w(5)='p_mean'
variable_names_av_w(6)='ssx_mean'
variable_names_av_w(7)='ssy_mean'
variable_names_av_w(8)='ssz_mean'


if (realgas.eq.1)then
write_variables_av_w=10

variable_names_av_w(1)='r_mean'
variable_names_av_w(2)='u_mean'
variable_names_av_w(3)='v_mean'
variable_names_av_w(4)='w_mean'
variable_names_av_w(5)='ttr_mean'
variable_names_av_w(6)='tvb_mean'
variable_names_av_w(7)='p_mean'
variable_names_av_w(8)='ssx_mean'
variable_names_av_w(9)='ssy_mean'
variable_names_av_w(10)='ssz_mean'
end if










end if
end if




































!wall instantaneous variables



 if (dimensiona.eq.3)then


 if (itestcase.eq.4)then
 write_variables_w=nof_variables+5
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='w'
variable_names_w(5)='pressure'
variable_names_w(6)='q'
variable_names_w(7)='ssx'
variable_names_w(8)='ssy'
variable_names_w(9)='ssz'
variable_names_w(10)='q_heat'

if (realgas.eq.1)then

 write_variables_w=11
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='w'
variable_names_w(5)='ttr'
variable_names_w(6)='tvb'
variable_names_w(7)='pressure'
variable_names_w(8)='ssx'
variable_names_w(9)='ssy'
variable_names_w(10)='ssz'
variable_names_w(11)='q_heat'

 end if

 end if
 if (itestcase.eq.3)then
 write_variables_w=nof_variables+1
variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='w'
variable_names_w(5)='pressure'
variable_names_w(6)='q'

if (realgas.eq.1)then

 write_variables_w=7
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='w'
variable_names_w(5)='ttr'
variable_names_w(6)='tvb'
variable_names_w(7)='pressure'

 end if







 end if



 if (itestcase.eq.-1)then
 write_variables_w=nof_variables+1
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='w'
variable_names_w(5)='pressure'
variable_names_w(6)='rho vf1'
variable_names_w(7)='rho vf2'
variable_names_w(8)='volume_fraction'
variable_names_w(9)='q'
 end if




else

if (itestcase.eq.4)then
 write_variables_w=nof_variables+4
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='pressure'
variable_names_w(5)='q'
variable_names_w(6)='ssx'
variable_names_w(7)='ssy'
variable_names_w(8)='q_heat'

if (realgas.eq.1)then

 write_variables_w=9
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='ttr'
variable_names_w(5)='tvb'
variable_names_w(6)='pressure'
variable_names_w(7)='ssx'
variable_names_w(8)='ssy'
variable_names_w(9)='q_heat'

 end if




 end if

 if (itestcase.eq.3)then
 write_variables_w=nof_variables+1
variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='pressure'
variable_names_w(5)='q'
if (realgas.eq.1)then

 write_variables_w=6
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='ttr'
variable_names_w(5)='tvb'
variable_names_w(6)='pressure'

 end if



 end if

 if (itestcase.eq.-1)then
 write_variables_w=nof_variables+1
 variable_names_w(1)='density'
variable_names_w(2)='u'
variable_names_w(3)='v'
variable_names_w(4)='pressure'
variable_names_w(5)='rho vf1'
variable_names_w(6)='rho vf2'
variable_names_w(7)='volume_fraction'
variable_names_w(8)='q'
 end if



end if






end subroutine specify_write_variables


subroutine parallel_vtk_combine(n)
!> @brief
!> this subroutine uses mpi-io for writing the vtk files
implicit none
integer,intent(in)::n
real,allocatable,dimension(:)::array2,array3,array4
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord
character(len=20)::proc,filex,proc3
real,allocatable,dimension(:)::array
logical::here1
real::in1,iocpt1,iocpt2,iocpt3,iocpt4,ptemp
integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
integer                     :: nbytes,eight
character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
character(len=200)          :: buffer
character(len=1)            :: lf
character(len=:),allocatable::vtu
real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
nbytes=1


 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0
kmaxe=xmpielrank(n)
kmaxn_p=xmpiall_v(n)


temp_cord=3


! 						if (n.eq.0)then
                               write(proc3,fmt='(i10)') it
                               filex="OUT_"//trim(adjustl(proc3))//".vtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)
!                             end if

		if (movement.eq.1)then

					k=1
	  do i=1,kmaxn_p
		rarray_part4(k:k+dims-1)=inoder4_cord(1:dims,my_nodesl(i))
		if (dimensiona.eq.2)then
		rarray_part4(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do


	  end if












				if (dimensiona.eq.3)then
				do i=1,kmaxe
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
				call cons2prim(n,leftv,mp_pinfl,gammal)	!r,u,v,w,p,y_n2,y_o2,y_no,y_n,y_o

				if (realgas.eq.1)then
				ptemp=leftv(dimensiona+2)
				vectco(1:nof_variables)=u_c_val(1,1:nof_variables,i)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
				call cons2div(n,vectco,mp_pinfl,gammal)
				leftv(1:nof_variables)=vectco(1:nof_variables)
				end if


					rarray_part1(i,1:nof_variables)=leftv(1:nof_variables)
										do j=nof_variables+1,write_variables-turbulenceequations
                                        if (multispecies.eq.1)then
											rarray_part1(i,j)=ielem_reduce(i)!ielem_vortex(1,i)
                                        else
											if (realgas.eq.1)then
											if (j.eq.nof_variables+1)then
											rarray_part1(i,j)=ptemp
											else
											rarray_part1(i,j)=ielem_reduce(i)!ielem_vortex(1,i)
											end if
											else
											rarray_part1(i,j)=ielem_vortex(1,i)
											end if
											 if (j.eq.write_variables-turbulenceequations)then
											 if (adda.eq.1)then
											 rarray_part1(i,j)=ielem_diss(i)
											 end if
											 end if
                                        end if        
                                        end do
                                        if (turbulenceequations.gt.0)then
                                        rarray_part1(i,write_variables)=u_ct_val(1,1,i)
                                        end if

				end do

				temp_node=8;temp_dims=3

				end if

				if (dimensiona.eq.2)then
				do i=1,kmaxe

				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)


				call cons2prim(n,leftv,mp_pinfl,gammal)	!r,u,v,p,y_n2,y_o2,y_no,y_n,y_o

				if (realgas.eq.1)then
				ptemp=leftv(dimensiona+2)
				vectco(1:nof_variables)=u_c_val(1,1:nof_variables,i)	!r,u,v,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
				call cons2div(n,vectco,mp_pinfl,gammal)
				leftv(1:nof_variables)=vectco(1:nof_variables)
				end if
					rarray_part1(i,1:nof_variables)=leftv(1:nof_variables)
										do j=nof_variables+1,write_variables-turbulenceequations
										if (multispecies.eq.1)then

										rarray_part1(i,j)=ielem_reduce(i)!ielem_vortex(1,i)
                                        else
											if (mood.eq.1)then
											rarray_part1(i,j)=ielem_mood_o(i)
											else
												if (dg.eq.1)then
													rarray_part1(i,j)=ielem_troubled(i)
												else
													if (realgas.eq.1)then
														if (j.eq.nof_variables+1)then
														rarray_part1(i,j)=ptemp
														else
														rarray_part1(i,j)=ielem_reduce(i)!vortex(1)
														end if
													else
														if (initcond.gt.100000)then
															rarray_part1(i,j)=u_e_val(1,1,i)
														else
															rarray_part1(i,j)=ielem_reduce(i)
														end if

													end if
												end if
											end if
                                        end if
										end do
										if (turbulenceequations.gt.0)then



                                        rarray_part1(i,write_variables)=u_ct_val(1,1,i)
                                        end if
				end do
				temp_node=4;temp_dims=3
				end if



temp_imaxe=imaxe
temp_imaxn=imaxn


if (n.eq.0)then

	!first write the header xml file from one mpi process

	lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
	 offset_temp=0
     write(offset_stamp,'(i16)')offset_temp
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify field pieces
    write(tempstamp1,'(i16)')imaxn
    write(tempstamp2,'(i16)')imaxe
    buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
           &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='     <PointData>'//lf;write(300) trim(buffer)
    buffer='     </PointData>'//lf;write(300) trim(buffer)
    buffer='     <CellData>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+size_of_real
    write(offset_stamp,'(i16)')offset_temp
    do i=1,write_variables
      buffer='        <DataArray type="Float64" Name="'//trim(variable_names(i))//'" '// &
                       'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      write(offset_stamp,'(i16)')offset_temp
    end do
    buffer='     </CellData>'//lf;write(300) trim(buffer)
    buffer='     <Points>'//lf;write(300) trim(buffer)
    buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    write(offset_stamp,'(i16)')offset_temp
    buffer='     </Points>'//lf;write(300) trim(buffer)
    ! specify necessary cell data
    buffer='      <Cells>'//lf;write(300) trim(buffer)
    ! connectivity
    buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+typ_countn_global*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! offsets
    buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! elem types
    buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    buffer='      </Cells>'//lf;write(300) trim(buffer)
    buffer='    </Piece>'//lf;write(300) trim(buffer)
    buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! prepare append section
    buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
    ! write leading data underscore
    buffer='_';write(300) trim(buffer)
	bytes = size_of_real
	close(300)


end if


call mpi_barrier(mpi_comm_world,ierror)


				call mpi_file_open(mpi_comm_world,vtu,mpi_mode_wronly + mpi_mode_append,mpi_info_null, fh, ierror)
				call mpi_file_get_position(fh, disp_in_file, ierror)
				disp_init=disp_in_file

				!----write time stamp----!
				if (n.eq.0)then
				call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
				bytes=size_of_real
				call mpi_file_write(fh, bytes, nbytes, mpi_integer, mpi_status_ignore, ierror)
				disp_in_file = disp_in_file + size_of_int
				call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
				call mpi_file_write(fh, t, nbytes, mpi_double_precision, mpi_status_ignore, ierror)
				disp_in_file=disp_in_file+size_of_real
				else
				disp_in_file=disp_in_file+size_of_int+size_of_real
				end if
				!end time stamp




				do i=1,write_variables
				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxe*size_of_real
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int
				!write variables---within loop
				call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatypex,'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh,rarray_part1(1:kmaxe,i),kmaxe*part1_end, mpi_double_precision,mpi_status_ignore,ierror)
				!end write variables---within loop
				disp_in_file=disp_in_file+temp_imaxe*size_of_real
				!end loop
				end do

! 				if (n.eq.0)print*,"location2",disp_in_file

				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxn*size_of_real*temp_dims
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int

				call mpi_file_set_view(fh, disp_in_file,mpi_double_precision,datatypez,'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh,rarray_part4,kmaxn_p*part4_end,mpi_double_precision,mpi_status_ignore, ierror)


				disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)

! 				if (n.eq.0)print*,"location3",disp_in_file


				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=size_of_int*typ_countn_global!temp_imaxe*size_of_int*temp_node
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int


				call mpi_file_set_view(fh,disp_in_file,mpi_integer,datatypey,'native',mpi_info_null,ierror)

				call mpi_file_write_all(fh,iarray_part2,typ_countn, mpi_integer,status,ierror)

				disp_in_file=disp_in_file+(size_of_int*typ_countn_global)!(temp_imaxe*size_of_int*temp_node)

! 				if (n.eq.0)print*,"location4",disp_in_file

				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxe*size_of_int
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int



				call mpi_file_set_view(fh, disp_in_file,mpi_integer,datatypexx, 'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh, iarray_part5,kmaxe*part1_end, mpi_integer,mpi_status_ignore, ierror)

				disp_in_file=disp_in_file+(temp_imaxe*size_of_int)


! 				if (n.eq.0)print*,"location5",disp_in_file


				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxe*size_of_int
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int


				call mpi_file_set_view(fh, disp_in_file,mpi_integer,datatypeyy, 'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh, iarray_part3,kmaxe*part1_end, mpi_integer,mpi_status_ignore, ierror)


				disp_in_file=disp_in_file+(temp_imaxe*size_of_int)

				call mpi_file_close(fh, ierror)
				call mpi_barrier(mpi_comm_world,ierror)





if (n.eq.0)then
open(300,file=filex,access='stream',position='append')
  lf = char(10)
  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)
end if


 deallocate(vtu)


call mpi_barrier(mpi_comm_world,ierror)







end subroutine parallel_vtk_combine


subroutine parallel_vtk_combine_av(n)
!> @brief
!> this subroutine uses mpi-io for writing the vtk files
implicit none
integer,intent(in)::n
real,allocatable,dimension(:)::array2,array3,array4
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord
character(len=40)::proc,filex,proc3
real,allocatable,dimension(:)::array
logical::here1
real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar,ptemp
real::in1,iocpt1,iocpt2,iocpt3,iocpt4
integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
integer                     :: nbytes,eight
character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
character(len=200)          :: buffer
character(len=1)            :: lf
character(len=:),allocatable::vtu
nbytes=1


 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0
kmaxe=xmpielrank(n)
kmaxn_p=xmpiall_v(n)


temp_cord=3






! 						if (n.eq.0)then
                               write(proc3,fmt='(i10)') it
                               filex="vol_aver"//trim(adjustl(proc3))//".vtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)
!                             end if

		if (movement.eq.1)then

					k=1
	  do i=1,kmaxn_p
		rarray_part4(k:k+dims-1)=inoder4_cord(1:dims,my_nodesl(i))
		if (dimensiona.eq.2)then
		rarray_part4(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do


	  end if












				if (dimensiona.eq.3)then
				do i=1,kmaxe
				leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
				call cons2prim(n,leftv,mp_pinfl,gammal)

				if (realgas.eq.1)then
				ptemp=leftv(dimensiona+2)
				vectco(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
				call cons2div(n,vectco,mp_pinfl,gammal)
				leftv(1:nof_variables)=vectco(1:nof_variables)
				end if



					rarray_part1(i,1:nof_variables)=leftv(1:nof_variables)
					if (realgas.eq.1)then
					rarray_part1(i,nof_variables+1)=ptemp
					do j=1,6
					rarray_part1(i,nof_variables+1+j)=u_c_val(1,nof_variables+1+j,i)
                     end do
					else
					do j=nof_variables+1,write_variables_av
					rarray_part1(i,j)=u_c_rms(j-nof_variables,i)
                     end do
                     end if
				end do

				temp_node=8;temp_dims=3

				end if





temp_imaxe=imaxe
temp_imaxn=imaxn


if (n.eq.0)then

	!first write the header xml file from one mpi process

	lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
	 offset_temp=0
     write(offset_stamp,'(i16)')offset_temp
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify field pieces
    write(tempstamp1,'(i16)')imaxn
    write(tempstamp2,'(i16)')imaxe
    buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
           &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='     <PointData>'//lf;write(300) trim(buffer)
    buffer='     </PointData>'//lf;write(300) trim(buffer)
    buffer='     <CellData>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+size_of_real
    write(offset_stamp,'(i16)')offset_temp
    do i=1,write_variables_av
      buffer='        <DataArray type="Float64" Name="'//trim(variable_names_av(i))//'" '// &
                       'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      write(offset_stamp,'(i16)')offset_temp
    end do
    buffer='     </CellData>'//lf;write(300) trim(buffer)
    buffer='     <Points>'//lf;write(300) trim(buffer)
    buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    write(offset_stamp,'(i16)')offset_temp
    buffer='     </Points>'//lf;write(300) trim(buffer)
    ! specify necessary cell data
    buffer='      <Cells>'//lf;write(300) trim(buffer)
    ! connectivity
    buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+typ_countn_global*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! offsets
    buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! elem types
    buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    buffer='      </Cells>'//lf;write(300) trim(buffer)
    buffer='    </Piece>'//lf;write(300) trim(buffer)
    buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! prepare append section
    buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
    ! write leading data underscore
    buffer='_';write(300) trim(buffer)
	bytes = size_of_real
	close(300)


end if


call mpi_barrier(mpi_comm_world,ierror)


				call mpi_file_open(mpi_comm_world,vtu,mpi_mode_wronly + mpi_mode_append,mpi_info_null, fh, ierror)
				call mpi_file_get_position(fh, disp_in_file, ierror)
				disp_init=disp_in_file

				!----write time stamp----!
				if (n.eq.0)then
				call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
				bytes=size_of_real
				call mpi_file_write(fh, bytes, nbytes, mpi_integer, mpi_status_ignore, ierror)
				disp_in_file = disp_in_file + size_of_int
				call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
				call mpi_file_write(fh, t, nbytes, mpi_double_precision, mpi_status_ignore, ierror)
				disp_in_file=disp_in_file+size_of_real
				else
				disp_in_file=disp_in_file+size_of_int+size_of_real
				end if
				!end time stamp




				do i=1,write_variables_av
				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxe*size_of_real
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int
				!write variables---within loop
				call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatypex,'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh,rarray_part1(1:kmaxe,i),kmaxe*part1_end, mpi_double_precision,mpi_status_ignore,ierror)
				!end write variables---within loop
				disp_in_file=disp_in_file+temp_imaxe*size_of_real
				!end loop
				end do

! 				if (n.eq.0)print*,"location2",disp_in_file

				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxn*size_of_real*temp_dims
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int

				call mpi_file_set_view(fh, disp_in_file,mpi_double_precision,datatypez,'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh,rarray_part4,kmaxn_p*part4_end,mpi_double_precision,mpi_status_ignore, ierror)


				disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)

! 				if (n.eq.0)print*,"location3",disp_in_file


				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=size_of_int*typ_countn_global!temp_imaxe*size_of_int*temp_node
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int


				call mpi_file_set_view(fh,disp_in_file,mpi_integer,datatypey,'native',mpi_info_null,ierror)

				call mpi_file_write_all(fh,iarray_part2,typ_countn, mpi_integer,status,ierror)

				disp_in_file=disp_in_file+(size_of_int*typ_countn_global)!(temp_imaxe*size_of_int*temp_node)

! 				if (n.eq.0)print*,"location4",disp_in_file

				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxe*size_of_int
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int



				call mpi_file_set_view(fh, disp_in_file,mpi_integer,datatypexx, 'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh, iarray_part5,kmaxe*part1_end, mpi_integer,mpi_status_ignore, ierror)

				disp_in_file=disp_in_file+(temp_imaxe*size_of_int)


! 				if (n.eq.0)print*,"location5",disp_in_file


				call mpi_file_set_view(fh, disp_in_file, mpi_integer,datatypeint,'native',mpi_info_null, ierror)

				if (n.eq.0)then
				bytes=temp_imaxe*size_of_int
				nbytes=1
				else
				bytes=0
				nbytes=0
				end if

				call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

				disp_in_file = disp_in_file + size_of_int


				call mpi_file_set_view(fh, disp_in_file,mpi_integer,datatypeyy, 'native',mpi_info_null, ierror)
				call mpi_file_write_all(fh, iarray_part3,kmaxe*part1_end, mpi_integer,mpi_status_ignore, ierror)


				disp_in_file=disp_in_file+(temp_imaxe*size_of_int)

				call mpi_file_close(fh, ierror)
				call mpi_barrier(mpi_comm_world,ierror)





if (n.eq.0)then
open(300,file=filex,access='stream',position='append')
  lf = char(10)
  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)
end if


 deallocate(vtu)


call mpi_barrier(mpi_comm_world,ierror)







end subroutine parallel_vtk_combine_av





subroutine parallel_vtk_combine_partitioned(n)
!> @brief
!> this subroutine uses mpi-io for writing the vtk files
implicit none
integer,intent(in)::n
real,allocatable,dimension(:)::array2,array3,array4
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord
character(len=20)::proc,proc3,proc5,proc6,proc7
character(len=90)::filex,filev
real,allocatable,dimension(:)::array
logical::here1
real::in1,iocpt1,iocpt2,iocpt3,iocpt4,ptemp
integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
integer                     :: nbytes,eight
character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
character(len=200)          :: buffer
character(len=1)            :: lf
character(len=:),allocatable::vtu
real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
nbytes=1

temp_cord=3
 	size_of_int=4
 	size_of_real=8

offset_temp=0
disp_in_file=0
tmp=0
disp_init=0
kmaxe=xmpielrank(n)



temp_cord=3





! 						if (n.eq.0)then
                               write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') n
                               filex="OUT_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".vtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)
!                             end if



	if (movement.eq.1)then
				k=1
	  do i=1,kmaxn
		nodes_vtu(k:k+dims-1)=inoder4_cord(1:dims,i)
		if (dimensiona.eq.2)then
		nodes_vtu(k+temp_cord-1:k+temp_cord-1)=0.0d0
		end if
		k=k+temp_cord
	  end do

   end if






				if (dimensiona.eq.3)then
				do i=1,kmaxe
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
				call cons2prim(n,leftv,mp_pinfl,gammal)	!r,u,v,w,p,y_n2,y_o2,y_no,y_n,y_o

				if (realgas.eq.1)then
				ptemp=leftv(dimensiona+2)
				vectco(1:nof_variables)=u_c_val(1,1:nof_variables,i)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
				call cons2div(n,vectco,mp_pinfl,gammal)
				leftv(1:nof_variables)=vectco(1:nof_variables)
				end if
					sol_vtu(i,1:nof_variables)=leftv(1:nof_variables)
										do j=nof_variables+1,write_variables-turbulenceequations
                                        if (multispecies.eq.1)then
										sol_vtu(i,j)=ielem_reduce(i)!ielem_vortex(1,i)
                                        else
                                        if (realgas.eq.1)then
											if (j.eq.nof_variables+1)then
											sol_vtu(i,j)=ptemp
											else
											sol_vtu(i,j)=ielem_vortex(1,i)
											end if
											else
											sol_vtu(i,j)=ielem_vortex(1,i)
											end if
											 if (j.eq.write_variables-turbulenceequations)then
											 if (adda.eq.1)then
											 sol_vtu(i,j)=ielem_diss(i)
											 end if
											 end if
                                        end if
                                        end do
                                        if (turbulenceequations.gt.0)then
                                        sol_vtu(i,write_variables)=u_ct_val(1,1,i)
                                        end if
				end do

				temp_node=8;temp_dims=3

				end if

				if (dimensiona.eq.2)then
				do i=1,kmaxe
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
				call cons2prim(n,leftv,mp_pinfl,gammal)	!r,u,v,p,y_n2,y_o2,y_no,y_n,y_o

				if (realgas.eq.1)then
				ptemp=leftv(dimensiona+2)
				vectco(1:nof_variables)=u_c_val(1,1:nof_variables,i)	!r,u,v,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
				call cons2div(n,vectco,mp_pinfl,gammal)
				leftv(1:nof_variables)=vectco(1:nof_variables)
				end if
					sol_vtu(i,1:nof_variables)=leftv(1:nof_variables)
										do j=nof_variables+1,write_variables-turbulenceequations
                                        if (multispecies.eq.1)then

										sol_vtu(i,j)=ielem_reduce(i)!ielem_vortex(1,i)
                                        else
                                        if (mood.eq.1)then
                                         sol_vtu(i,j)=ielem_mood_o(i)
                                        else
										if (dg.eq.1)then
											sol_vtu(i,j)=ielem_troubled(i)
										else
										if (realgas.eq.1)then
											if (j.eq.nof_variables+1)then
											sol_vtu(i,j)=ptemp
											else
											sol_vtu(i,j)=ielem_vortex(1,i)
											end if
										else
											sol_vtu(i,j)=ielem_vortex(1,i)
										end if
                                        end if
                                        end if
                                        end if
										end do
										if (turbulenceequations.gt.0)then
                                        sol_vtu(i,write_variables)=u_ct_val(1,1,i)
                                        end if

				end do
				temp_node=4;temp_dims=3
				end if





temp_imaxe=kmaxe
temp_imaxn=kmaxn


	!first write the header xml file from one mpi process

	lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
	 offset_temp=0
     write(offset_stamp,'(i16)')offset_temp
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify field pieces
    write(tempstamp1,'(i16)')kmaxn
    write(tempstamp2,'(i16)')kmaxe
    buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
           &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='     <PointData>'//lf;write(300) trim(buffer)
    buffer='     </PointData>'//lf;write(300) trim(buffer)
    buffer='     <CellData>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+size_of_real
    write(offset_stamp,'(i16)')offset_temp
    do i=1,write_variables
      buffer='        <DataArray type="Float64" Name="'//trim(variable_names(i))//'" '// &
                       'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      write(offset_stamp,'(i16)')offset_temp
    end do
    buffer='     </CellData>'//lf;write(300) trim(buffer)
    buffer='     <Points>'//lf;write(300) trim(buffer)
    buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    write(offset_stamp,'(i16)')offset_temp
    buffer='     </Points>'//lf;write(300) trim(buffer)
    ! specify necessary cell data
    buffer='      <Cells>'//lf;write(300) trim(buffer)
    ! connectivity
    buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+(typ_countn)*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! offsets
    buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! elem types
    buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    buffer='      </Cells>'//lf;write(300) trim(buffer)
    buffer='    </Piece>'//lf;write(300) trim(buffer)
    buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! prepare append section
    buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
    ! write leading data underscore
    buffer='_';write(300) trim(buffer)
	bytes = size_of_real





				!----write time stamp----!
				bytes=size_of_real
				write(300)bytes,t

				!end time stamp



				!----write variables----!
				do j=1,write_variables
				bytes=temp_imaxe*size_of_real
					 write(300)bytes,sol_vtu(1:kmaxe,j)

				end do
				!end----write variables----!


				!write nodes now!
				bytes=kmaxn*3*size_of_real
					 write(300)bytes,nodes_vtu(1:kmaxn*3)

				!end write nodes now!


				!write connectivity now!
				bytes=typ_countn*size_of_int
					 write(300)bytes, connect_vtu(1:typ_countn)

				!end write connectivity now!


				!write offsets now!
				bytes=kmaxe*size_of_int
					 write(300)bytes,offset_vtu(1:kmaxe)

				!end write nodes now!


				!write types now!
				bytes=kmaxe*size_of_int
					 write(300)bytes,type_vtu(1:kmaxe)

				!end write nodes now!

  lf = char(10)
  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)



 deallocate(vtu)

	call mpi_barrier(mpi_comm_world,ierror)


 if (n.eq.0)then


							write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') n
                               filex="PAR_"//trim(adjustl(proc3))//".pvtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)




                               lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf;write(300) trim(buffer)
!     buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;write(300) trim(buffer)
    write(tempstamp1,'(f17.9)')t
     buffer='        '//trim(adjustl(tempstamp1))//lf;write(300) trim(buffer)
	buffer='      </DataArray>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='    <PCellData>'//lf;write(300) trim(buffer)
    do i=1,write_variables
      buffer='        <PDataArray type="Float64" Name="'//trim(variable_names(i))//'" '// &
                       'format="appended"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
    end do
	buffer='    </PCellData>'//lf;write(300) trim(buffer)
    buffer='    <PPoints>'//lf;write(300) trim(buffer)

     buffer='        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf;write(300) trim(buffer)
    buffer='    </PPoints>'//lf;write(300) trim(buffer)
    buffer='    <PCells>'//lf;write(300) trim(buffer)
    buffer='        <PDataArray type="Int32" Name="connectivity"/>'//lf;write(300) trim(buffer)
     buffer='        <PDataArray type="Int32" Name="offsets"/>'//lf;write(300) trim(buffer)
     buffer='        <PDataArray type="Int32" Name="types"/>'//lf;write(300) trim(buffer)
     buffer='    </PCells>'//lf;write(300) trim(buffer)
   do procx=0,isize-1

				write(proc6,fmt='(i10)') it
				write(proc7,fmt='(i10)') procx
		buffer='    <Piece source="out_'//trim(adjustl(proc6))//"_"//trim(adjustl(proc7))//'.vtu"/>'//lf;write(300) trim(buffer)
   end do
  buffer='  </PUnstructuredGrid>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)



	deallocate(vtu)



 end if

	call mpi_barrier(mpi_comm_world,ierror)




end subroutine parallel_vtk_combine_partitioned


subroutine parallel_vtk_combine_partitioned_wall(n)
!> @brief
!> this subroutine uses mpi-io for writing the vtk files
implicit none
integer,intent(in)::n
real,allocatable,dimension(:)::array2,array3,array4
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord,kkd_i,kkd
character(len=20)::proc,proc3,proc5,proc6,proc7
character(len=90)::filex,filev
real,allocatable,dimension(:)::array
logical::here1
real::in1,iocpt1,iocpt2,iocpt3,iocpt4,ptemp
integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
integer                     :: nbytes,eight
character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
character(len=200)          :: buffer
character(len=1)            :: lf
character(len=:),allocatable::vtu
real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::shear_temp
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






                               write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') n
                               filex="SURF_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".vtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)




			if (movement.eq.1)then
						k=1
			do i=1,kmaxn
				nodes_vtu_w(k:k+dims-1)=inoder4_cord(1:dims,i)
				if (dimensiona.eq.2)then
				nodes_vtu_w(k+temp_cord-1:k+temp_cord-1)=0.0d0
				end if
				k=k+temp_cord
			end do

		end if








				do i=1,iloopx
				facex=wall_l(i,2)
				iconsidered=wall_l(i,1)
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
							if (dimensiona.eq.3)then
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=3;temp_dims=3

								else
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=2;temp_dims=2
								end if



								if (realgas.eq.1)then
								ptemp=leftv(dimensiona+2)
								vectco(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
								call cons2div(n,vectco,mp_pinfl,gammal)
								leftv(1:nof_variables)=vectco(1:nof_variables)
								end if

									if (realgas.eq.1)then
									sol_vtu_w(i,1:dimensiona+3)=leftv(1:dimensiona+3)
									sol_vtu_w(i,dimensiona+4)=ptemp

									else


									!the next variable is always going to be an auxiliary
									sol_vtu_w(i,1:nof_variables)=leftv(1:nof_variables)
									sol_vtu_w(i,nof_variables+1:nof_variables+1)=ielem_vortex(1,iconsidered)
									end if





					kkd_i=nof_variables+1

									if (realgas.eq.1)then
									kkd_i=dimensiona+4
									end if

								    if (itestcase.eq.4)then

										if (dimensiona.eq.3)then
											do kkd=1,4




										select case(kkd)
											case(1)

											call shear_x(iconsidered,facex,shear_temp)
											case (2)
											call shear_y(iconsidered,facex,shear_temp)
											case(3)
											call shear_z(iconsidered,facex,shear_temp)

											case(4)

											call heat_x(iconsidered,facex,shear_temp)

											end select

											sol_vtu_w(i,kkd_i+kkd)=shear_temp

											end do
										end if

										if (dimensiona.eq.2)then
											do kkd=1,3




										select case(kkd)
											case(1)

											call shear_x2d(iconsidered,facex,shear_temp)
											case (2)
											call shear_y2d(iconsidered,facex,shear_temp)
											case(3)

											call heat_x2d(iconsidered,facex,shear_temp)
											end select

											sol_vtu_w(i,kkd_i+kkd)=shear_temp

											end do
										end if



									end if

								end do






temp_imaxe=iloopx
temp_imaxn=kmaxn



	!first write the header xml file from one mpi process

	lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
	 offset_temp=0
     write(offset_stamp,'(i16)')offset_temp
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify field pieces
    write(tempstamp1,'(i16)')kmaxn
    write(tempstamp2,'(i16)')iloopx
    buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
           &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='     <PointData>'//lf;write(300) trim(buffer)
    buffer='     </PointData>'//lf;write(300) trim(buffer)
    buffer='     <CellData>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+size_of_real
    write(offset_stamp,'(i16)')offset_temp
    do i=1,write_variables_w
      buffer='        <DataArray type="Float64" Name="'//trim(variable_names_w(i))//'" '// &
                       'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      write(offset_stamp,'(i16)')offset_temp
    end do
    buffer='     </CellData>'//lf;write(300) trim(buffer)
    buffer='     <Points>'//lf;write(300) trim(buffer)
    buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    write(offset_stamp,'(i16)')offset_temp
    buffer='     </Points>'//lf;write(300) trim(buffer)
    ! specify necessary cell data
    buffer='      <Cells>'//lf;write(300) trim(buffer)
    ! connectivity
    buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+(typ_countn_w)*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! offsets
    buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! elem types
    buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    buffer='      </Cells>'//lf;write(300) trim(buffer)
    buffer='    </Piece>'//lf;write(300) trim(buffer)
    buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! prepare append section
    buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
    ! write leading data underscore
    buffer='_';write(300) trim(buffer)
	bytes = size_of_real





				!----write time stamp----!
				bytes=size_of_real
				write(300)bytes,t

				!end time stamp



				!----write variables----!

				do j=1,write_variables_w
				bytes=temp_imaxe*size_of_real
					 write(300)bytes,sol_vtu_w(1:iloopx,j)

				end do
				!end----write variables----!


				!write nodes now!
				bytes=kmaxn*3*size_of_real
					 write(300)bytes,nodes_vtu_w(1:kmaxn*3)

				!end write nodes now!


				!write connectivity now!
				bytes=typ_countn_w*size_of_int
					 write(300)bytes, connect_vtu_w(1:typ_countn_w)

				!end write connectivity now!


				!write offsets now!
				bytes=iloopx*size_of_int
					 write(300)bytes,offset_vtu_w(1:iloopx)

				!end write nodes now!


				!write types now!
				bytes=iloopx*size_of_int
					 write(300)bytes,type_vtu_w(1:iloopx)

				!end write nodes now!

  lf = char(10)
  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)



 deallocate(vtu)


end if


	call mpi_barrier(mpi_comm_world,ierror)


 if (n.eq.0)then


							write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') n
                               filex="PAR_SURF_"//trim(adjustl(proc3))//".pvtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)




                               lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf;write(300) trim(buffer)
!     buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;write(300) trim(buffer)
    write(tempstamp1,'(f17.9)')t
     buffer='        '//trim(adjustl(tempstamp1))//lf;write(300) trim(buffer)
	buffer='      </DataArray>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='    <PCellData>'//lf;write(300) trim(buffer)
    do i=1,write_variables_w
      buffer='        <PDataArray type="Float64" Name="'//trim(variable_names_w(i))//'" '// &
                       'format="appended"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
    end do
	buffer='    </PCellData>'//lf;write(300) trim(buffer)
    buffer='    <PPoints>'//lf;write(300) trim(buffer)

     buffer='        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf;write(300) trim(buffer)
    buffer='    </PPoints>'//lf;write(300) trim(buffer)
    buffer='    <PCells>'//lf;write(300) trim(buffer)
    buffer='        <PDataArray type="Int32" Name="connectivity"/>'//lf;write(300) trim(buffer)
     buffer='        <PDataArray type="Int32" Name="offsets"/>'//lf;write(300) trim(buffer)
     buffer='        <PDataArray type="Int32" Name="types"/>'//lf;write(300) trim(buffer)
     buffer='    </PCells>'//lf;write(300) trim(buffer)
   do procx=0,isize-1
				if (wallcount_cpu_g(procx).gt.0)then
				write(proc6,fmt='(i10)') it
				write(proc7,fmt='(i10)') procx
		buffer='    <Piece source="SURF_'//trim(adjustl(proc6))//"_"//trim(adjustl(proc7))//'.vtu"/>'//lf;write(300) trim(buffer)
				end if
   end do
  buffer='  </PUnstructuredGrid>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)



	deallocate(vtu)



 end if

	call mpi_barrier(mpi_comm_world,ierror)




end subroutine parallel_vtk_combine_partitioned_wall




subroutine parallel_vtk_combine_partitioned_wall_av(n)
!> @brief
!> this subroutine uses mpi-io for writing the vtk files
implicit none
integer,intent(in)::n
real,allocatable,dimension(:)::array2,array3,array4
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord,kkd_i,kkd
character(len=20)::proc,proc3,proc5,proc6,proc7
character(len=90)::filex,filev
real,allocatable,dimension(:)::array
logical::here1
real::in1,iocpt1,iocpt2,iocpt3,iocpt4
integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
integer                     :: nbytes,eight
character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
character(len=200)          :: buffer
character(len=1)            :: lf
character(len=:),allocatable::vtu
real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal,ptemp
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::shear_temp
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






                               write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') n
                               filex="SURF_AV_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".vtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)




			if (movement.eq.1)then
						k=1
			do i=1,kmaxn
				nodes_vtu_w(k:k+dims-1)=inoder4_cord(1:dims,i)
				if (dimensiona.eq.2)then
				nodes_vtu_w(k+temp_cord-1:k+temp_cord-1)=0.0d0
				end if
				k=k+temp_cord
			end do

		end if








				do i=1,iloopx
				facex=wall_l(i,2)
				iconsidered=wall_l(i,1)
				leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)
							if (dimensiona.eq.3)then
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=3;temp_dims=3

								else
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=2;temp_dims=2
								end if


								if (realgas.eq.1)then
								ptemp=leftv(dimensiona+2)
								vectco(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
								call cons2div(n,vectco,mp_pinfl,gammal)
								leftv(1:nof_variables)=vectco(1:nof_variables)
								end if

									if (realgas.eq.1)then
									sol_vtu_w(i,1:dimensiona+3)=leftv(1:dimensiona+3)
									sol_vtu_w(i,dimensiona+4)=ptemp
									else
									sol_vtu_w(i,1:nof_variables)=leftv(1:nof_variables)
									end if



! 					sol_vtu_w(i,1:nof_variables)=leftv(1:nof_variables)




					kkd_i=nof_variables

									if (realgas.eq.1)then
									kkd_i=dimensiona+4
									end if


								    if (itestcase.eq.4)then

										if (dimensiona.eq.3)then
											do kkd=1,3




										select case(kkd)
											case(1)

											call shear_x_av(iconsidered,facex,shear_temp)
											case (2)
											call shear_y_av(iconsidered,facex,shear_temp)
											case(3)
											call shear_z_av(iconsidered,facex,shear_temp)
											end select

											sol_vtu_w(i,kkd_i+kkd)=shear_temp

											end do
										end if

										if (dimensiona.eq.2)then
											do kkd=1,2




										select case(kkd)
											case(1)

											call shear_x2d_av(iconsidered,facex,shear_temp)
											case (2)
											call shear_y2d_av(iconsidered,facex,shear_temp)

											end select

											sol_vtu_w(i,kkd_i+kkd)=shear_temp

											end do
										end if



									end if

								end do






temp_imaxe=iloopx
temp_imaxn=kmaxn



	!first write the header xml file from one mpi process

	lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
	 offset_temp=0
     write(offset_stamp,'(i16)')offset_temp
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify field pieces
    write(tempstamp1,'(i16)')kmaxn
    write(tempstamp2,'(i16)')iloopx
    buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
           &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='     <PointData>'//lf;write(300) trim(buffer)
    buffer='     </PointData>'//lf;write(300) trim(buffer)
    buffer='     <CellData>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+size_of_real
    write(offset_stamp,'(i16)')offset_temp
    do i=1,write_variables_av_w
      buffer='        <DataArray type="Float64" Name="'//trim(variable_names_av_w(i))//'" '// &
                       'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
      write(offset_stamp,'(i16)')offset_temp
    end do
    buffer='     </CellData>'//lf;write(300) trim(buffer)
    buffer='     <Points>'//lf;write(300) trim(buffer)
    buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
    write(offset_stamp,'(i16)')offset_temp
    buffer='     </Points>'//lf;write(300) trim(buffer)
    ! specify necessary cell data
    buffer='      <Cells>'//lf;write(300) trim(buffer)
    ! connectivity
    buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+(typ_countn_w)*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! offsets
    buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
    write(offset_stamp,'(i16)')offset_temp
    ! elem types
    buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
                     'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
    ! buffer='        </DataArray>'//lf;write(300) trim(buffer)
    buffer='      </Cells>'//lf;write(300) trim(buffer)
    buffer='    </Piece>'//lf;write(300) trim(buffer)
    buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
    ! prepare append section
    buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
    ! write leading data underscore
    buffer='_';write(300) trim(buffer)
	bytes = size_of_real





				!----write time stamp----!
				bytes=size_of_real
				write(300)bytes,t

				!end time stamp



				!----write variables----!

				do j=1,write_variables_av_w
				bytes=temp_imaxe*size_of_real
					 write(300)bytes,sol_vtu_w(1:iloopx,j)

				end do
				!end----write variables----!


				!write nodes now!
				bytes=kmaxn*3*size_of_real
					 write(300)bytes,nodes_vtu_w(1:kmaxn*3)

				!end write nodes now!


				!write connectivity now!
				bytes=typ_countn_w*size_of_int
					 write(300)bytes, connect_vtu_w(1:typ_countn_w)

				!end write connectivity now!


				!write offsets now!
				bytes=iloopx*size_of_int
					 write(300)bytes,offset_vtu_w(1:iloopx)

				!end write nodes now!


				!write types now!
				bytes=iloopx*size_of_int
					 write(300)bytes,type_vtu_w(1:iloopx)

				!end write nodes now!

  lf = char(10)
  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)



 deallocate(vtu)


end if


	call mpi_barrier(mpi_comm_world,ierror)


 if (n.eq.0)then


							write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') n
                               filex="PAR_SURF_AV_"//trim(adjustl(proc3))//".pvtu"
                               itrimm=len_trim(filex)
                               allocate(character(len=itrimm)::vtu)
                               vtu=filex(1:itrimm)




                               lf = char(10)
   ! write file name
    open(300,file=vtu,access='stream')
    ! write header
    buffer='<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf;write(300) trim(buffer)
    ! write unstructured grid type
    buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;write(300) trim(buffer)
    ! write solution time type
    buffer='    <FieldData>'//lf;write(300) trim(buffer)
    buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf;write(300) trim(buffer)
!     buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;write(300) trim(buffer)
    write(tempstamp1,'(f17.9)')t
     buffer='        '//trim(adjustl(tempstamp1))//lf;write(300) trim(buffer)
	buffer='      </DataArray>'//lf;write(300) trim(buffer)
	buffer='    </FieldData>'//lf;write(300) trim(buffer)
    ! specify point data
    buffer='    <PCellData>'//lf;write(300) trim(buffer)
    do i=1,write_variables_av_w
      buffer='        <PDataArray type="Float64" Name="'//trim(variable_names_av_w(i))//'" '// &
                       'format="appended"/>'//lf;write(300) trim(buffer)
      offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
    end do
	buffer='    </PCellData>'//lf;write(300) trim(buffer)
    buffer='    <PPoints>'//lf;write(300) trim(buffer)

     buffer='        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf;write(300) trim(buffer)
    buffer='    </PPoints>'//lf;write(300) trim(buffer)
    buffer='    <PCells>'//lf;write(300) trim(buffer)
    buffer='        <PDataArray type="Int32" Name="connectivity"/>'//lf;write(300) trim(buffer)
     buffer='        <PDataArray type="Int32" Name="offsets"/>'//lf;write(300) trim(buffer)
     buffer='        <PDataArray type="Int32" Name="types"/>'//lf;write(300) trim(buffer)
     buffer='    </PCells>'//lf;write(300) trim(buffer)
   do procx=0,isize-1
				if (wallcount_cpu_g(procx).gt.0)then
				write(proc6,fmt='(i10)') it
				write(proc7,fmt='(i10)') procx
		buffer='    <Piece source="surf_av_'//trim(adjustl(proc6))//"_"//trim(adjustl(proc7))//'.vtu"/>'//lf;write(300) trim(buffer)
				end if
   end do
  buffer='  </PUnstructuredGrid>'//lf;write(300) trim(buffer)
  buffer='</VTKFile>'//lf;write(300) trim(buffer)
  close(300)



	deallocate(vtu)



 end if

	call mpi_barrier(mpi_comm_world,ierror)




end subroutine parallel_vtk_combine_partitioned_wall_av





subroutine parallel_vtk_combine_partitioned_av(n)
	!> @brief
	!> this subroutine uses mpi-io for writing the vtk files
	implicit none
	integer,intent(in)::n
	real,allocatable,dimension(:)::array2,array3,array4
	integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord
	character(len=20)::proc,proc3,proc5,proc6,proc7
	character(len=90)::filex,filev
	real,allocatable,dimension(:)::array
	logical::here1
	real::in1,iocpt1,iocpt2,iocpt3,iocpt4,ptemp
	integer:: disp_in_file,procx ,tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
	integer                     :: nbytes,eight
	character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
	character(len=200)          :: buffer
	character(len=1)            :: lf
	character(len=:),allocatable::vtu
	real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
	nbytes=1
	
	temp_cord=3
		 size_of_int=4
		 size_of_real=8
	
	offset_temp=0
	disp_in_file=0
	tmp=0
	disp_init=0
	kmaxe=xmpielrank(n)
	
	
	
	temp_cord=3
	



	
	
	
	
	! 						if (n.eq.0)then
								   write(proc3,fmt='(i10)') it
								   write(proc5,fmt='(i10)') n
								   filex="VOL_AVER_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".vtu"
								   itrimm=len_trim(filex)
								   allocate(character(len=itrimm)::vtu)
								   vtu=filex(1:itrimm)
	!                             end if
	
	
	
		if (movement.eq.1)then
					k=1
		  do i=1,kmaxn
			nodes_vtu(k:k+dims-1)=inoder4_cord(1:dims,i)
			if (dimensiona.eq.2)then
			nodes_vtu(k+temp_cord-1:k+temp_cord-1)=0.0d0
			end if
			k=k+temp_cord
		  end do
	
	   end if
	
	
	
	
	
	
					if (dimensiona.eq.3)then
					do i=1,kmaxe
					leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
					call cons2prim(n,leftv,mp_pinfl,gammal)

					if (realgas.eq.1)then
					ptemp=leftv(dimensiona+2)
					vectco(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
					call cons2div(n,vectco,mp_pinfl,gammal)
					leftv(1:nof_variables)=vectco(1:nof_variables)
					end if


						sol_vtu(i,1:nof_variables)=leftv(1:nof_variables)
						if (realgas.eq.1)then
						sol_vtu(i,nof_variables+1)=ptemp

						do j=1,6
							sol_vtu(i,nof_variables+1+j)=u_c_val(1,nof_variables+1+j,i)
						end do

						else

						do j=nof_variables+1,write_variables_av
							sol_vtu(i,j)=u_c_rms(j-nof_variables,i)
						end do
						end if
					end do
	
					temp_node=8;temp_dims=3
	
					end if
	
					
	
	
	
	
	
	temp_imaxe=kmaxe
	temp_imaxn=kmaxn
	
	
		!first write the header xml file from one mpi process
	
		lf = char(10)
	   ! write file name
		open(300,file=vtu,access='stream')
		! write header
		buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
		! write unstructured grid type
		buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
		! write solution time type
		buffer='    <FieldData>'//lf;write(300) trim(buffer)
		 offset_temp=0
		 write(offset_stamp,'(i16)')offset_temp
		buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
						 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
		buffer='    </FieldData>'//lf;write(300) trim(buffer)
		! specify field pieces
		write(tempstamp1,'(i16)')kmaxn
		write(tempstamp2,'(i16)')kmaxe
		buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
			   &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
		! specify point data
		buffer='     <PointData>'//lf;write(300) trim(buffer)
		buffer='     </PointData>'//lf;write(300) trim(buffer)
		buffer='     <CellData>'//lf;write(300) trim(buffer)
		offset_temp=offset_temp+size_of_int+size_of_real
		write(offset_stamp,'(i16)')offset_temp
		do i=1,write_variables_av
		  buffer='        <DataArray type="Float64" Name="'//trim(variable_names_av(i))//'" '// &
						   'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
		  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
		  write(offset_stamp,'(i16)')offset_temp
		end do
		buffer='     </CellData>'//lf;write(300) trim(buffer)
		buffer='     <Points>'//lf;write(300) trim(buffer)
		buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
						 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
		! buffer='        </DataArray>'//lf;write(300) trim(buffer)
		offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
		write(offset_stamp,'(i16)')offset_temp
		buffer='     </Points>'//lf;write(300) trim(buffer)
		! specify necessary cell data
		buffer='      <Cells>'//lf;write(300) trim(buffer)
		! connectivity
		buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
						 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
		offset_temp=offset_temp+size_of_int+(typ_countn)*size_of_int
		write(offset_stamp,'(i16)')offset_temp
		! offsets
		buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
						 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
		offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
		write(offset_stamp,'(i16)')offset_temp
		! elem types
		buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
						 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
		! buffer='        </DataArray>'//lf;write(300) trim(buffer)
		buffer='      </Cells>'//lf;write(300) trim(buffer)
		buffer='    </Piece>'//lf;write(300) trim(buffer)
		buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
		! prepare append section
		buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
		! write leading data underscore
		buffer='_';write(300) trim(buffer)
		bytes = size_of_real
	
	
	
	
	
					!----write time stamp----!
					bytes=size_of_real
					write(300)bytes,t
	
					!end time stamp
	
					
	
					!----write variables----!
					do j=1,write_variables_av
					bytes=temp_imaxe*size_of_real
						 write(300)bytes,sol_vtu(1:kmaxe,j)
	
					end do
					!end----write variables----!
	
	
					!write nodes now!
					bytes=kmaxn*3*size_of_real
						 write(300)bytes,nodes_vtu(1:kmaxn*3)
	
					!end write nodes now!
	
	
					!write connectivity now!
					bytes=typ_countn*size_of_int
						 write(300)bytes, connect_vtu(1:typ_countn)
	
					!end write connectivity now!
	
	
					!write offsets now!
					bytes=kmaxe*size_of_int
						 write(300)bytes,offset_vtu(1:kmaxe)
	
					!end write nodes now!
	
	
					!write types now!
					bytes=kmaxe*size_of_int
						 write(300)bytes,type_vtu(1:kmaxe)
	
					!end write nodes now!
	
	  lf = char(10)
	  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
	  buffer='</VTKFile>'//lf;write(300) trim(buffer)
	  close(300)
	
	
	
	 deallocate(vtu)
	
		call mpi_barrier(mpi_comm_world,ierror)
	
	
	 if (n.eq.0)then
	
	
								write(proc3,fmt='(i10)') it
								   write(proc5,fmt='(i10)') n
								   filex="par_vol_aver_"//trim(adjustl(proc3))//".pvtu"
								   itrimm=len_trim(filex)
								   allocate(character(len=itrimm)::vtu)
								   vtu=filex(1:itrimm)
	
	
	
	
								   lf = char(10)
	   ! write file name
		open(300,file=vtu,access='stream')
		! write header
		buffer='<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//lf;write(300) trim(buffer)
		! write unstructured grid type
		buffer='  <PUnstructuredGrid GhostLevel="0">'//lf;write(300) trim(buffer)
		! write solution time type
		buffer='    <FieldData>'//lf;write(300) trim(buffer)
		buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="ascii">'//lf;write(300) trim(buffer)
	!     buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1">'//lf;write(300) trim(buffer)
		write(tempstamp1,'(f17.9)')t
		 buffer='        '//trim(adjustl(tempstamp1))//lf;write(300) trim(buffer)
		buffer='      </DataArray>'//lf;write(300) trim(buffer)
		buffer='    </FieldData>'//lf;write(300) trim(buffer)
		! specify point data
		buffer='    <PCellData>'//lf;write(300) trim(buffer)
		do i=1,write_variables_av
		  buffer='        <PDataArray type="Float64" Name="'//trim(variable_names_av(i))//'" '// &
						   'format="appended"/>'//lf;write(300) trim(buffer)
		  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
		end do
		buffer='    </PCellData>'//lf;write(300) trim(buffer)
		buffer='    <PPoints>'//lf;write(300) trim(buffer)
	
		 buffer='        <PDataArray type="Float64" Name="Coordinates" NumberOfComponents="3"/>'//lf;write(300) trim(buffer)
		buffer='    </PPoints>'//lf;write(300) trim(buffer)
		buffer='    <PCells>'//lf;write(300) trim(buffer)
		buffer='        <PDataArray type="Int32" Name="connectivity"/>'//lf;write(300) trim(buffer)
		 buffer='        <PDataArray type="Int32" Name="offsets"/>'//lf;write(300) trim(buffer)
		 buffer='        <PDataArray type="Int32" Name="types"/>'//lf;write(300) trim(buffer)
		 buffer='    </PCells>'//lf;write(300) trim(buffer)
	   do procx=0,isize-1
	
					write(proc6,fmt='(i10)') it
					write(proc7,fmt='(i10)') procx
			buffer='    <Piece source="vol_aver_'//trim(adjustl(proc6))//"_"//trim(adjustl(proc7))//'.vtu"/>'//lf;write(300) trim(buffer)
	   end do
	  buffer='  </PUnstructuredGrid>'//lf;write(300) trim(buffer)
	  buffer='</VTKFile>'//lf;write(300) trim(buffer)
	  close(300)
	
	
	
		deallocate(vtu)
	
	
	
	 end if
	
		call mpi_barrier(mpi_comm_world,ierror)
	
	
	
	
	end subroutine parallel_vtk_combine_partitioned_av





	subroutine parallel_vtk_combine_wall(n)
		!> @brief
		!> this subroutine uses mpi-io for writing the vtk files
		implicit none
		integer,intent(in)::n
		real,allocatable,dimension(:)::array2,array3,array4
		integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord,kkd_i,kkd,typ_countn_global
		character(len=20)::proc,filex,proc3
		real,allocatable,dimension(:)::array
		logical::here1
		real::in1,iocpt1,iocpt2,iocpt3,iocpt4,ptemp
		integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
		integer                     :: nbytes,eight
		character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
		character(len=200)          :: buffer
		character(len=1)            :: lf
		character(len=:),allocatable::vtu
		real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::shear_temp
integer::iconsidered,facex
		nbytes=1
		
		
			 size_of_int=4
			 size_of_real=8
		
		offset_temp=0
		disp_in_file=0
		tmp=0
		disp_init=0
		kmaxe=xmpielrank(n)
		kmaxn_p=xmpiall_v(n)
		if (dimensiona.eq.3)then
 		temp_node=3;temp_dims=3
 		else
 		temp_node=2;temp_dims=3
 		end if
		
		temp_cord=3



		
		

									   write(proc3,fmt='(i10)') it
									   filex="SURF_"//trim(adjustl(proc3))//".vtu"
									   itrimm=len_trim(filex)
									   allocate(character(len=itrimm)::vtu)
									   vtu=filex(1:itrimm)

		
				if (movement.eq.1)then

						k=1
						do i=1,kmaxn_p
									wrarray_part4(k:k+dims-1)=inoder4_cord(1:dims,my_nodesl(i))
									if (dimensiona.eq.2)then
									wrarray_part4(k+temp_cord-1:k+temp_cord-1)=0.0d0
									end if
									k=k+temp_cord
						end do

				end if
		
		
		
		
		
		
			  			!loop the correct number of elements that are bounded 
			  			if (iloopx.gt.0)then
							do i=1,iloopx
								facex=wall_l(i,2)
								iconsidered=wall_l(i,1)
								leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)


								if (dimensiona.eq.3)then
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=3;temp_dims=3

								else
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=2;temp_dims=3
								end if



								if (realgas.eq.1)then
								ptemp=leftv(dimensiona+2)
								vectco(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
								call cons2div(n,vectco,mp_pinfl,gammal)
								leftv(1:nof_variables)=vectco(1:nof_variables)
								end if










									if (realgas.eq.1)then
									wrarray_part1(i,1:dimensiona+3)=leftv(1:dimensiona+3)
									wrarray_part1(i,dimensiona+4)=ptemp

									else


									!the next variable is always going to be an auxiliary
									wrarray_part1(i,1:nof_variables)=leftv(1:nof_variables)
									wrarray_part1(i,nof_variables+1:nof_variables+1)=ielem_vortex(1,iconsidered)
									end if
						
									kkd_i=nof_variables+1

									if (realgas.eq.1)then
									kkd_i=dimensiona+4
									end if


								    if (itestcase.eq.4)then
									
										if (dimensiona.eq.3)then
											do kkd=1,4




										select case(kkd)
											case(1)
											
											call shear_x(iconsidered,facex,shear_temp)
											case (2)
											call shear_y(iconsidered,facex,shear_temp)
											case(3)
											call shear_z(iconsidered,facex,shear_temp)

											case(4)

											call heat_x(iconsidered,facex,shear_temp)



											end select

											wrarray_part1(i,kkd_i+kkd)=shear_temp

											end do
										end if

										if (dimensiona.eq.2)then
											do kkd=1,3




										select case(kkd)
											case(1)
											
											call shear_x2d(iconsidered,facex,shear_temp)
											case (2)
											call shear_y2d(iconsidered,facex,shear_temp)
											case(3)

											call heat_x2d(iconsidered,facex,shear_temp)

											
											end select

											wrarray_part1(i,kkd_i+kkd)=shear_temp

											end do
										end if



									end if
		
								end do
		
							end if




		
		
		
		
		temp_imaxe=iwmaxe !total number of wall elements in the domain
		temp_imaxn=imaxn	 !imaxn	!we need the total number of nodes in the domain
		
		call mpi_barrier(mpi_comm_world,ierror)




		if (n.eq.0)then
		
			!first write the header xml file from one mpi process
		
			lf = char(10)
		   ! write file name
			open(300,file=vtu,access='stream')
			! write header
			buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
			! write unstructured grid type
			buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
			! write solution time type
			buffer='    <FieldData>'//lf;write(300) trim(buffer)
			 offset_temp=0
			 write(offset_stamp,'(i16)')offset_temp
			buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			buffer='    </FieldData>'//lf;write(300) trim(buffer)
			! specify field pieces
			write(tempstamp1,'(i16)')temp_imaxn
			write(tempstamp2,'(i16)')temp_imaxe
			buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
				   &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
			! specify point data
			buffer='     <PointData>'//lf;write(300) trim(buffer)
			buffer='     </PointData>'//lf;write(300) trim(buffer)
			buffer='     <CellData>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+size_of_real
			write(offset_stamp,'(i16)')offset_temp
			do i=1,write_variables_w
			  buffer='        <DataArray type="Float64" Name="'//trim(variable_names_w(i))//'" '// &
							   'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
			  write(offset_stamp,'(i16)')offset_temp
			end do
			buffer='     </CellData>'//lf;write(300) trim(buffer)
			buffer='     <Points>'//lf;write(300) trim(buffer)
			buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			! buffer='        </DataArray>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
			write(offset_stamp,'(i16)')offset_temp
			buffer='     </Points>'//lf;write(300) trim(buffer)
			! specify necessary cell data
			buffer='      <Cells>'//lf;write(300) trim(buffer)
			! connectivity
			buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+typ_countn_global_w*size_of_int
			write(offset_stamp,'(i16)')offset_temp
			! offsets
			buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
			write(offset_stamp,'(i16)')offset_temp
			! elem types
			buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			! buffer='        </DataArray>'//lf;write(300) trim(buffer)
			buffer='      </Cells>'//lf;write(300) trim(buffer)
			buffer='    </Piece>'//lf;write(300) trim(buffer)
			buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
			! prepare append section
			buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
			! write leading data underscore
			buffer='_';write(300) trim(buffer)
			bytes = size_of_real
			close(300)
		
		
		end if
		


		
		call mpi_barrier(mpi_comm_world,ierror)

		
						call mpi_file_open(mpi_comm_world,vtu,mpi_mode_wronly + mpi_mode_append,mpi_info_null, fh, ierror)
						call mpi_file_get_position(fh, disp_in_file, ierror)
						disp_init=disp_in_file
		
						!----write time stamp----!
						if (n.eq.0)then
						call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
						bytes=size_of_real
						call mpi_file_write(fh, bytes, nbytes, mpi_integer, mpi_status_ignore, ierror)
						disp_in_file = disp_in_file + size_of_int
						call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
						call mpi_file_write(fh, t, nbytes, mpi_double_precision, mpi_status_ignore, ierror)
						disp_in_file=disp_in_file+size_of_real
						else
						disp_in_file=disp_in_file+size_of_int+size_of_real
						end if
						!end time stamp
		


		
		
! 						do i=1,write_variables
						do i=1,write_variables_w




						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)
		
						if (n.eq.0)then
						bytes=temp_imaxe*size_of_real
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if
		
						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)
		
						disp_in_file = disp_in_file + size_of_int

						call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,wdatatypex,'native',mpi_info_null, ierror)
						call mpi_file_write_all(fh,wrarray_part1(:,i),iloopx*wpart1_end, mpi_double_precision,mpi_status_ignore,ierror)
						!end write variables---within loop
						disp_in_file=disp_in_file+temp_imaxe*size_of_real
						!end loop
						end do
		
		! 				if (n.eq.0)print*,"location2",disp_in_file
		
						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)
		
						if (n.eq.0)then
						bytes=temp_imaxn*size_of_real*temp_dims
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if
		
						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)
		
						disp_in_file = disp_in_file + size_of_int

						call mpi_file_set_view(fh, disp_in_file,mpi_double_precision,wdatatypez,'native',mpi_info_null, ierror)
						call mpi_file_write_all(fh,wrarray_part4,kmaxn_p*wpart4_end,mpi_double_precision,mpi_status_ignore, ierror)

						disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)




						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)
		
						if (n.eq.0)then
						bytes=size_of_int*typ_countn_global_w!temp_imaxe*size_of_int*temp_node
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if
		




						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)
		
						disp_in_file = disp_in_file + size_of_int


						call mpi_file_set_view(fh,disp_in_file,mpi_integer,wdatatypey,'native',mpi_info_null,ierror)


						call mpi_file_write_all(fh,wiarray_part2,typ_countn_w, mpi_integer,status,ierror)


		
						disp_in_file=disp_in_file+(size_of_int*typ_countn_global_w)!(temp_imaxe*size_of_int*temp_node)
		



						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)
		
						if (n.eq.0)then
						bytes=temp_imaxe*size_of_int
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if
		
						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)
		
						disp_in_file = disp_in_file + size_of_int
		
		
						call mpi_barrier(mpi_comm_world,ierror)

						call mpi_file_set_view(fh, disp_in_file,mpi_integer,wdatatypexx, 'native',mpi_info_null, ierror)

						call mpi_file_write_all(fh,wiarray_part5,iloopx*wpart1_end, mpi_integer,mpi_status_ignore, ierror)
		
						disp_in_file=disp_in_file+(temp_imaxe*size_of_int)
		

						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)
		
						if (n.eq.0)then
						bytes=temp_imaxe*size_of_int
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if
		
						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)
		
						disp_in_file = disp_in_file + size_of_int
		
		


						call mpi_file_set_view(fh, disp_in_file,mpi_integer,wdatatypeyy, 'native',mpi_info_null, ierror)

						call mpi_file_write_all(fh, wiarray_part3,iloopx*wpart1_end, mpi_integer,mpi_status_ignore, ierror)
		


		
						disp_in_file=disp_in_file+(temp_imaxe*size_of_int)
		
						call mpi_file_close(fh, ierror)
						call mpi_barrier(mpi_comm_world,ierror)
		
		
		
		
		if (n.eq.0)then
		open(300,file=filex,access='stream',position='append')
		  lf = char(10)
		  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
		  buffer='</VTKFile>'//lf;write(300) trim(buffer)
		  close(300)
		end if
		
		
		 deallocate(vtu)
		
		
		call mpi_barrier(mpi_comm_world,ierror)
		
		

		
		
		
		
		end subroutine parallel_vtk_combine_wall


		subroutine parallel_vtk_combine_wall_av(n)
		!> @brief
		!> this subroutine uses mpi-io for writing the vtk files
		implicit none
		integer,intent(in)::n
		real,allocatable,dimension(:)::array2,array3,array4
		integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,dip,n_end,ifg,kmaxn_p,itrimm,temp_cord,kkd_i,kkd,typ_countn_global
		character(len=20)::proc,filex,proc3
		real,allocatable,dimension(:)::array
		logical::here1
		real::in1,iocpt1,iocpt2,iocpt3,iocpt4,ptemp
		integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init,offset_temp,bytes,temp_imaxe,temp_imaxn,temp_node,temp_dims,size_of_real,size_of_int
		integer                     :: nbytes,eight
		character(len=35)           :: offset_stamp,tempstamp1,tempstamp2
		character(len=200)          :: buffer
		character(len=1)            :: lf
		character(len=:),allocatable::vtu
		real,dimension(1:nof_variables)::leftv,vectco
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::shear_temp
integer::iconsidered,facex
		nbytes=1


			 size_of_int=4
			 size_of_real=8

		offset_temp=0
		disp_in_file=0
		tmp=0
		disp_init=0
		kmaxe=xmpielrank(n)
		kmaxn_p=xmpiall_v(n)


		temp_cord=3
		if (dimensiona.eq.3)then
 		temp_node=3;temp_dims=3
 		else
 		temp_node=2;temp_dims=3
 		end if








									   write(proc3,fmt='(i10)') it
									   filex="SURF_AV"//trim(adjustl(proc3))//".vtu"
									   itrimm=len_trim(filex)
									   allocate(character(len=itrimm)::vtu)
									   vtu=filex(1:itrimm)


				if (movement.eq.1)then

						k=1
						do i=1,kmaxn_p
									wrarray_part4(k:k+dims-1)=inoder4_cord(1:dims,my_nodesl(i))
									if (dimensiona.eq.2)then
									wrarray_part4(k+temp_cord-1:k+temp_cord-1)=0.0d0
									end if
									k=k+temp_cord
						end do

				end if


			  			!loop the correct number of elements that are bounded
			  			if (iloopx.gt.0)then
							do i=1,iloopx
								facex=wall_l(i,2)
								iconsidered=wall_l(i,1)
								leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)


								if (dimensiona.eq.3)then
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=3;temp_dims=3

								else
								call cons2prim(n,leftv,mp_pinfl,gammal)
								temp_node=2;temp_dims=3
								end if

								if (realgas.eq.1)then
								ptemp=leftv(dimensiona+2)
								vectco(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)	!r,u,v,w,ttr,tv,y_n2,y_o2,y_no,y_n,y_o
								call cons2div(n,vectco,mp_pinfl,gammal)
								leftv(1:nof_variables)=vectco(1:nof_variables)
								end if

									if (realgas.eq.1)then
									wrarray_part1(i,1:dimensiona+3)=leftv(1:dimensiona+3)
									wrarray_part1(i,dimensiona+4)=ptemp
									else
									wrarray_part1(i,1:nof_variables)=leftv(1:nof_variables)
									end if


									kkd_i=nof_variables


									if (realgas.eq.1)then
									kkd_i=dimensiona+4
									end if

								    if (itestcase.eq.4)then

										if (dimensiona.eq.3)then
											do kkd=1,3




										select case(kkd)
											case(1)

											call shear_x_av(iconsidered,facex,shear_temp)
											case (2)
											call shear_y_av(iconsidered,facex,shear_temp)
											case(3)
											call shear_z_av(iconsidered,facex,shear_temp)
											end select

											wrarray_part1(i,kkd_i+kkd)=shear_temp

											end do
										end if

										if (dimensiona.eq.2)then
											do kkd=1,2




										select case(kkd)
											case(1)

											call shear_x2d_av(iconsidered,facex,shear_temp)
											case (2)
											call shear_y2d_av(iconsidered,facex,shear_temp)

											end select

											wrarray_part1(i,kkd_i+kkd)=shear_temp

											end do
										end if



									end if

								end do

							end if






		temp_imaxe=iwmaxe !total number of wall elements in the domain
		temp_imaxn=imaxn	 !imaxn	!we need the total number of nodes in the domain

		call mpi_barrier(mpi_comm_world,ierror)




		if (n.eq.0)then

			!first write the header xml file from one mpi process

			lf = char(10)
		   ! write file name
			open(300,file=vtu,access='stream')
			! write header
			buffer='<VTKFile type="UnstructuredGrid" version="2.2" byte_order="LittleEndian" header_type="UInt32">'//lf;write(300) trim(buffer)
			! write unstructured grid type
			buffer='  <UnstructuredGrid>'//lf;write(300) trim(buffer)
			! write solution time type
			buffer='    <FieldData>'//lf;write(300) trim(buffer)
			 offset_temp=0
			 write(offset_stamp,'(i16)')offset_temp
			buffer='      <DataArray type="Float64" Name="TimeValue" NumberOfTuples="1" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			buffer='    </FieldData>'//lf;write(300) trim(buffer)
			! specify field pieces
			write(tempstamp1,'(i16)')temp_imaxn
			write(tempstamp2,'(i16)')temp_imaxe
			buffer='    <Piece NumberOfPoints="'//trim(adjustl(tempstamp1))//'" &
				   &NumberOfCells="'//trim(adjustl(tempstamp2))//'">'//lf;write(300) trim(buffer)
			! specify point data
			buffer='     <PointData>'//lf;write(300) trim(buffer)
			buffer='     </PointData>'//lf;write(300) trim(buffer)
			buffer='     <CellData>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+size_of_real
			write(offset_stamp,'(i16)')offset_temp
			do i=1,write_variables_av_w
			  buffer='        <DataArray type="Float64" Name="'//trim(variable_names_av_w(i))//'" '// &
							   'format="appended" offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			  offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_real
			  write(offset_stamp,'(i16)')offset_temp
			end do
			buffer='     </CellData>'//lf;write(300) trim(buffer)
			buffer='     <Points>'//lf;write(300) trim(buffer)
			buffer='        <DataArray type="Float64" Name="Coordinates" NumberOfComponents="3" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			! buffer='        </DataArray>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+temp_dims*temp_imaxn*size_of_real
			write(offset_stamp,'(i16)')offset_temp
			buffer='     </Points>'//lf;write(300) trim(buffer)
			! specify necessary cell data
			buffer='      <Cells>'//lf;write(300) trim(buffer)
			! connectivity
			buffer='        <DataArray type="Int32" Name="connectivity" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+typ_countn_global_w*size_of_int
			write(offset_stamp,'(i16)')offset_temp
			! offsets
			buffer='        <DataArray type="Int32" Name="offsets" format="appended" ' // &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			offset_temp=offset_temp+size_of_int+temp_imaxe*size_of_int
			write(offset_stamp,'(i16)')offset_temp
			! elem types
			buffer='        <DataArray type="Int32" Name="types" format="appended" '// &
							 'offset="'//trim(adjustl(offset_stamp))//'"/>'//lf;write(300) trim(buffer)
			! buffer='        </DataArray>'//lf;write(300) trim(buffer)
			buffer='      </Cells>'//lf;write(300) trim(buffer)
			buffer='    </Piece>'//lf;write(300) trim(buffer)
			buffer='  </UnstructuredGrid>'//lf;write(300) trim(buffer)
			! prepare append section
			buffer='  <AppendedData encoding="raw">'//lf;write(300) trim(buffer)
			! write leading data underscore
			buffer='_';write(300) trim(buffer)
			bytes = size_of_real
			close(300)


		end if




		call mpi_barrier(mpi_comm_world,ierror)


						call mpi_file_open(mpi_comm_world,vtu,mpi_mode_wronly + mpi_mode_append,mpi_info_null, fh, ierror)
						call mpi_file_get_position(fh, disp_in_file, ierror)
						disp_init=disp_in_file

						!----write time stamp----!
						if (n.eq.0)then
						call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
						bytes=size_of_real
						call mpi_file_write(fh, bytes, nbytes, mpi_integer, mpi_status_ignore, ierror)
						disp_in_file = disp_in_file + size_of_int
						call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
						call mpi_file_write(fh, t, nbytes, mpi_double_precision, mpi_status_ignore, ierror)
						disp_in_file=disp_in_file+size_of_real
						else
						disp_in_file=disp_in_file+size_of_int+size_of_real
						end if
						!end time stamp





!
						do i=1,write_variables_av_w




						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)

						if (n.eq.0)then
						bytes=temp_imaxe*size_of_real
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if

						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

						disp_in_file = disp_in_file + size_of_int

						call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,wdatatypex,'native',mpi_info_null, ierror)
						call mpi_file_write_all(fh,wrarray_part1(:,i),iloopx*wpart1_end, mpi_double_precision,mpi_status_ignore,ierror)
						!end write variables---within loop
						disp_in_file=disp_in_file+temp_imaxe*size_of_real
						!end loop
						end do

		! 				if (n.eq.0)print*,"location2",disp_in_file

						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)

						if (n.eq.0)then
						bytes=temp_imaxn*size_of_real*temp_dims
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if

						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

						disp_in_file = disp_in_file + size_of_int

						call mpi_file_set_view(fh, disp_in_file,mpi_double_precision,wdatatypez,'native',mpi_info_null, ierror)
						call mpi_file_write_all(fh,wrarray_part4,kmaxn_p*wpart4_end,mpi_double_precision,mpi_status_ignore, ierror)

						disp_in_file=disp_in_file+(temp_imaxn*size_of_real*temp_dims)




						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)

						if (n.eq.0)then
						bytes=size_of_int*typ_countn_global_w!temp_imaxe*size_of_int*temp_node
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if





						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

						disp_in_file = disp_in_file + size_of_int


						call mpi_file_set_view(fh,disp_in_file,mpi_integer,wdatatypey,'native',mpi_info_null,ierror)


						call mpi_file_write_all(fh,wiarray_part2,typ_countn_w, mpi_integer,status,ierror)



						disp_in_file=disp_in_file+(size_of_int*typ_countn_global_w)!(temp_imaxe*size_of_int*temp_node)




						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)

						if (n.eq.0)then
						bytes=temp_imaxe*size_of_int
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if

						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

						disp_in_file = disp_in_file + size_of_int


						call mpi_barrier(mpi_comm_world,ierror)

						call mpi_file_set_view(fh, disp_in_file,mpi_integer,wdatatypexx, 'native',mpi_info_null, ierror)

						call mpi_file_write_all(fh,wiarray_part5,iloopx*wpart1_end, mpi_integer,mpi_status_ignore, ierror)

						disp_in_file=disp_in_file+(temp_imaxe*size_of_int)


						call mpi_file_set_view(fh, disp_in_file, mpi_integer,wdatatypeint,'native',mpi_info_null, ierror)

						if (n.eq.0)then
						bytes=temp_imaxe*size_of_int
						nbytes=1
						else
						bytes=0
						nbytes=0
						end if

						call mpi_file_write_all(fh,bytes,nbytes,mpi_integer,mpi_status_ignore, ierror)

						disp_in_file = disp_in_file + size_of_int




						call mpi_file_set_view(fh, disp_in_file,mpi_integer,wdatatypeyy, 'native',mpi_info_null, ierror)

						call mpi_file_write_all(fh, wiarray_part3,iloopx*wpart1_end, mpi_integer,mpi_status_ignore, ierror)




						disp_in_file=disp_in_file+(temp_imaxe*size_of_int)

						call mpi_file_close(fh, ierror)
						call mpi_barrier(mpi_comm_world,ierror)




		if (n.eq.0)then
		open(300,file=filex,access='stream',position='append')
		  lf = char(10)
		  buffer=lf//'  </AppendedData>'//lf;write(300) trim(buffer)
		  buffer='</VTKFile>'//lf;write(300) trim(buffer)
		  close(300)
		end if


		 deallocate(vtu)


		call mpi_barrier(mpi_comm_world,ierror)







		end subroutine parallel_vtk_combine_wall_av



subroutine checkpointv2(n)
!> @brief
!> this subroutine is writing the checkpointing files
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella
real,allocatable,dimension(:)::valuesa,valuess
real,allocatable,dimension(:,:)::xbin
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj
character(len=20)::proc,restfile,proc3
real,allocatable,dimension(:)::igint,tgint
 kmaxe=xmpielrank(n)


icpuid=n
	restfile='restart2.dat'
	if (n.eq.0)then
	open(1086,file=restfile,form='unformatted',status='replace',action='write')
	if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	write(1086)it,zero
	      write(1086)initialres(1)
	      write(1086)initialres(2)
	      write(1086)initialres(3)
	      write(1086)initialres(4)
	      write(1086)initialres(5)
	      if ( turbulence .eq. 1) then
	      if (turbulencemodel.eq.1)then
	      write(1086)initialres(6)
	      end if
	      if (turbulencemodel.eq.2)then
	      write(1086)initialres(6)
	      write(1086)initialres(7)
	      end if
	      end if
	      
	      
	else
	write (1086)it,t
	
	if (initcond.eq.95)then
	write(1086)taylor
	end if
	      
	end if
	end if
	
	kmaxe=xmpielrank(n)
    
dumg=kmaxe


call mpi_barrier(mpi_comm_world,ierror) !not needed

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


    if (n.eq.0)then
	    allocate(icella(imaxp*isize))
	    icella=0

    end if

    call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

!     call mpi_barrier(mpi_comm_world,ierror)
    deallocate (icell)

    if (n.eq.0) then
    allocate(valuesa(imaxp*isize))
    allocate(xbin(imaxe,5+turbulenceequations+passivescalar))
    valuesa=zero
    end if
    
    
    allocate(valuess(imaxp));valuess=zero

  
  
if (turbulence.eq.1)then
do jj=1,5+turbulenceequations+passivescalar
 do i=1,kmaxe
      if (jj.gt.5) then

		      valuess(i)=u_ct_val(1,jj-5,i)
      else
	valuess(i)=u_c_val(1,jj,i)
      end if
 end do
    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i),jj)=valuesa(i)
	end if
    end do
    end if
    
end do
else
  do jj=1,5
	do i=1,kmaxe
		valuess(i)=u_c_val(1,jj,i)
	end do
	call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
	
	    if (n.eq.0)then
	    do i=1,imaxp*isize
		if (icella(i).gt.0)then
		xbin(icella(i),jj)=valuesa(i)
		end if
	    end do
	    end if
    
  end do

end if

!   call mpi_barrier(mpi_comm_world,ierror)

    if (n.eq.0)then

    do i=1,imaxe
    write(1086)xbin(i,1:nof_variables+turbulenceequations+passivescalar)
    end do

    deallocate(xbin,icella,valuesa)
    close(1086)
    end if
    
deallocate(valuess)








end subroutine checkpointv2




subroutine checkpoint2d(n)
!> @brief
!> this subroutine uses mpi-io for writing the checkpointing files for 2d
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella,dispt
real,allocatable,dimension(:)::valuesa,valuess,array2
real,allocatable,dimension(:,:)::xbin
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,size_of_real,size_of_int,dip,n_end,datatype
character(len=20)::proc,restfile,proc3
real,allocatable,dimension(:)::igint,tgint,array
real::in1
logical::here1
integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init
disp_in_file=0
tmp=0
disp_init=0
 kmaxe=xmpielrank(n)


size_of_int=4
size_of_real=8
icpuid=n
call mpi_barrier(mpi_comm_world,ierror)

if (dg.eq.1)then
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+turbulenceequations+passivescalar)*(idegfree+1)))
else
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+turbulenceequations+passivescalar)))	!i allocate in memory the pattern of access of data in terms of displacement, and in terms of blocklength, and finaly an array with this processor data
end if

       if (dg.eq.1)then
     do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*((nof_variables+turbulenceequations+passivescalar)*(idegfree+1))
      end do
      

      n_end=(nof_variables+turbulenceequations+passivescalar)*(idegfree+1)
     
     else


      do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*(nof_variables+turbulenceequations+passivescalar)
      end do

      n_end=nof_variables+turbulenceequations+passivescalar

      end if
      
      
      
      if (dg.eq.1)then
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	    k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(1,1:nof_variables,i)
	      k=k+nof_variables
	      array2(k:k+turbulenceequations+passivescalar-1)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
	      k=k+turbulenceequations+passivescalar
	  end do
      else
	  k=1
	  do i=1,kmaxe
        do j=1,nof_variables
	      array2(k:k+idegfree)=u_c_valdg(1,j,1:idegfree+1,i)
	     
	      
	      k=k+(idegfree+1)
        end do
	  end do
      end if
      
      
      else
      
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	    k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(1,1:nof_variables,i)
	      k=k+nof_variables
	      array2(k:k+turbulenceequations+passivescalar-1)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
	      k=k+turbulenceequations+passivescalar
	  end do
      else
	  k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(1,1:nof_variables,i)
	      k=k+nof_variables
	  end do
      end if
      
      end if
      restfile='RESTART.dat'
       if (n.eq.0)then
	  !inquire (file=restfile,exist=here1)
	  !if (heress) then
	 inquire (file=restfile,exist=here1)
	 
	  
	  call mpi_file_delete(restfile,mpi_info_null,ierror)
	  
	 

	  !end if

   end if

      call mpi_barrier(mpi_comm_world,ierror)
      
    !create type first of indexed block
    call mpi_type_create_indexed_block(kmaxe,n_end,dispt,mpi_double_precision,datatype,ierror)
    call mpi_type_commit(datatype,ierror)
    

    allocate(array(1:nof_variables+turbulenceequations+passivescalar))
    
	
	
	
	call mpi_file_open(mpi_comm_world, restfile,mpi_mode_wronly + mpi_mode_create,mpi_info_null, fh, ierror)
	
	
	 
	  
	
	if (n.eq.0)then
	
	    if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	          call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_write(fh, it, 1, mpi_integer, mpi_status_ignore,ierror)
		    disp_in_file = disp_in_file + size_of_int	!1
		    call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_write(fh, initialres(1:nof_variables+turbulenceequations),nof_variables+turbulenceequations, mpi_double_precision, mpi_status_ignore, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_variables+turbulenceequations)	!3
	    else
		  call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_write(fh, it, 1, mpi_integer, mpi_status_ignore, ierror)
		    disp_in_file = disp_in_file + size_of_int 	!4
		    
		  call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_write(fh, t, 1, mpi_double_precision, mpi_status_ignore,ierror)
		    disp_in_file = disp_in_file + size_of_real    !5
		      if (initcond.eq.95)then
		      call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
		      call mpi_file_write(fh, taylor, 1, mpi_double_precision, mpi_status_ignore, ierror)
			    disp_in_file = disp_in_file + size_of_real !6
		      end if  
	    end if
	else
	      if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
		disp_in_file = disp_in_file + size_of_int 
		disp_in_file = disp_in_file + size_of_real*(nof_variables+turbulenceequations)
	      else
		  disp_in_file = disp_in_file + size_of_int 
		  disp_in_file = disp_in_file + size_of_real
		  if (initcond.eq.95)then
		  disp_in_file = disp_in_file + size_of_real 
		  end if
	      end if
	      
	      
	end if
	
	
	call mpi_barrier(mpi_comm_world, ierror)
	call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatype, 'native',mpi_info_null, ierror)
	call mpi_file_write_all(fh, array2, kmaxe*n_end, mpi_double_precision,mpi_status_ignore, ierror)        
        call mpi_file_close(fh, ierror)
	call mpi_type_free(datatype,ierror)
          deallocate(array,dispt,array2)
          
          
	
	
	call mpi_barrier(mpi_comm_world, ierror)
	
	
	



end subroutine checkpoint2d




subroutine checkpointav(n) 
!> @brief
!> this subroutine uses mpi-io for writing the averaged checkpointing files
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella,dispt
real,allocatable,dimension(:)::valuesa,valuess,array2
real,allocatable,dimension(:,:)::xbin
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,fh,size_of_real,size_of_int,dip,ista,iend,n_end,datatype
 character(len=20)::proc,restfile,proc3
 real,allocatable,dimension(:)::igint,tgint,array
 integer(kind=mpi_offset_kind) :: disp_in_file, tmp,disp_init
 logical::here1
 kmaxe=xmpielrank(n)
disp_in_file=0
tmp=0
disp_init=0
 
 
 
 
 
 
 size_of_int=4
size_of_real=8
icpuid=n
call mpi_barrier(mpi_comm_world,ierror)
 



	
	
 allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+turbulenceequations+passivescalar+6+passivescalar)))
    do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*(nof_variables+turbulenceequations+passivescalar+6+passivescalar)
      end do

 restfile='RESTARTAV.dat'
 
  n_end=nof_variables+turbulenceequations+passivescalar+6+passivescalar
  
  
 
 
  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	    k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(ind1,1:nof_variables,i)
	      k=k+nof_variables
	      array2(k:k+turbulenceequations+passivescalar-1)=u_ct_val(ind1,1:turbulenceequations+passivescalar,i)
	      k=k+turbulenceequations+passivescalar
	      array2(k:k+6+passivescalar-1)=u_c_rms(1:6+passivescalar,i)
	      k=k+6+passivescalar
	  end do
      else
	  k=1
	  do i=1,kmaxe
	      array2(k:k+nof_variables-1)=u_c_val(ind1,1:nof_variables,i)
	      k=k+nof_variables
	      array2(k:k+6+passivescalar-1)=u_c_rms(1:6+passivescalar,i)
	      k=k+6+passivescalar
	  end do
      end if
  if (n.eq.0)then
	 inquire (file=restfile,exist=here1)
	  call mpi_file_delete(restfile,mpi_info_null,ierror)
	end if
 
 call mpi_barrier(mpi_comm_world,ierror)
 
 
 call mpi_type_create_indexed_block(kmaxe,n_end,dispt,mpi_double_precision,datatype,ierror)
    call mpi_type_commit(datatype,ierror)
 
 
 
	
	
	call mpi_barrier(mpi_comm_world,ierror)
	call mpi_file_open(mpi_comm_world, restfile,mpi_mode_wronly + mpi_mode_create,mpi_info_null, fh, ierror)
	call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatype, 'native',mpi_info_null, ierror)
	call mpi_file_write_all(fh, array2, kmaxe*n_end, mpi_double_precision,mpi_status_ignore, ierror)        
        call mpi_file_close(fh, ierror)
	call mpi_type_free(datatype,ierror)
          
          
          
	deallocate(dispt,array2)
	
	call mpi_barrier(mpi_comm_world, ierror)
	
	
	
	
 
	

end subroutine checkpointav


subroutine rest_read(n)
!> @brief
!> this subroutine uses mpi-io for reading the checkpointing files
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella,dispt
real,allocatable,dimension(:)::valuesa,valuess,array2
real,allocatable,dimension(:)::rg,arg
character(len=20)::proc,restfile,proc3
integer:: prev_turbequation,initial,iii,i,k,j,jx,qqp,inc,kmaxe,jkn,ki,iterr,jx2,fh,size_of_real,size_of_int,dip,n_end,datatype
real,allocatable,dimension(:)::igint,tgint,array
integer(kind=mpi_offset_kind) :: disp_in_file, tmp
logical::here
disp_in_file=0
tmp=0

kmaxe=xmpielrank(n)


prev_turbequation=0
if (prev_turbmodel.eq.1) then
prev_turbequation=1
end if 
if (prev_turbmodel.eq.2) then
prev_turbequation=2
end if 

 size_of_int=4
size_of_real=8

!$omp master
if (dg.eq.1)then
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+prev_turbequation+passivescalar)*(idegfree+1)))
else
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+prev_turbequation+passivescalar)))	
end if

    if (dg.eq.1)then
     do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*((nof_variables+turbulenceequations+passivescalar)*(idegfree+1))
      end do
    n_end=(nof_variables+turbulenceequations+passivescalar)*(idegfree+1)
     
     else


    do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*(nof_variables+prev_turbequation+passivescalar)
      end do
      
      
      n_end=nof_variables+prev_turbequation+lamps
      
      end if
      
      call mpi_type_create_indexed_block(kmaxe,n_end,dispt,mpi_double_precision,datatype,ierror)
    call mpi_type_commit(datatype,ierror)
    
    restfile='RESTART.dat'
    
    call mpi_file_open(mpi_comm_world, restfile,mpi_mode_rdonly,mpi_info_null, fh, ierror)
    
    
    
	    if (ires_unsteady.eq.0)then
! 	   if ((rungekutta .ge. 5).and.(rungekutta .lt. 11)) then
	          call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_read(fh, it, 1, mpi_integer, mpi_status_ignore,ierror)
		    disp_in_file = disp_in_file + size_of_int	!1
		    call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_read(fh, initialres(1:nof_variables+prev_turbequation),nof_variables+prev_turbequation, mpi_double_precision, mpi_status_ignore, ierror)
		    disp_in_file = disp_in_file + size_of_real*(nof_variables+prev_turbequation)	!3
	    else
		  call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_read(fh, it, 1, mpi_integer, mpi_status_ignore, ierror)
		  
		    disp_in_file = disp_in_file + size_of_int 	!4
		    
		  call mpi_file_seek(fh, disp_in_file, mpi_seek_set, ierror)
		  call mpi_file_read(fh, t, 1, mpi_double_precision, mpi_status_ignore,ierror)
		 
		    disp_in_file = disp_in_file + size_of_real    !5
		      if (initcond.eq.95)then
		      call mpi_file_seek(fh, disp_in_file, mpi_seek_set,ierror)
		      call mpi_file_read(fh, taylor, 1, mpi_double_precision, mpi_status_ignore, ierror)
			    disp_in_file = disp_in_file + size_of_real !6
		      end if  
	    end if
      
      
	call mpi_barrier(mpi_comm_world, ierror)
	call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatype, 'native',mpi_info_null, ierror)
	call mpi_file_read_all(fh, array2, kmaxe*n_end, mpi_double_precision,mpi_status_ignore, ierror)        
        call mpi_file_close(fh, ierror)
	call mpi_type_free(datatype,ierror)
	
	
	
	if (dg.eq.1)then
	
	
	if ((prev_turbmodel.gt.0).or.(lamps.gt.0))then
	
	    k=1
	    do i=1,kmaxe
		u_c_val(1,1:nof_variables,i)=array2(k:k+nof_variables-1)
		k=k+nof_variables
		    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		    u_ct_val(1,1:turbulenceequations+passivescalar,i)=array2(k:k+prev_turbequation+lamps-1)
		    end if
		k=k+prev_turbmodel+lamps
		
	    end do
	else
	    k=1
	    do i=1,kmaxe
	    do j=1,nof_variables
	    
		u_c_valdg(1,j,1:idegfree+1,i)=array2(k:k+idegfree)
		
		
		k=k+(idegfree+1)
		end do
		
		
		      if (turbulence.eq.1)then
			    if (turbulencemodel.eq.1)then
				u_ct_val(1,1,i)=visc*turbinit
			    else
				u_ct_Val(1,1,i)=1.5*(i_turb_inlet*ufreestream)**2
				u_ct_val(1,2,i)=(c_mu_inlet**(-0.25))*sqrt(u_ct_val(1,1,ki))&
					/l_turb_inlet*rg(1)
			    end if
			endif
		k=k+prev_turbmodel+lamps
		
	    end do
	end if
	
	
	
	
	
	
	
	
	
	
	else
	if ((prev_turbmodel.gt.0).or.(lamps.gt.0))then
	
	    k=1
	    do i=1,kmaxe
		u_c_val(1,1:nof_variables,i)=array2(k:k+nof_variables-1)
		k=k+nof_variables
		    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		    u_ct_val(1,1:turbulenceequations+passivescalar,i)=array2(k:k+prev_turbequation+lamps-1)
		    end if
		k=k+prev_turbmodel+lamps
		
	    end do
	else
	    k=1
	    do i=1,kmaxe
		u_c_val(1,1:nof_variables,i)=array2(k:k+nof_variables-1)
		k=k+nof_variables
		      if (turbulence.eq.1)then
			    if (turbulencemodel.eq.1)then
				u_ct_val(1,1,i)=visc*turbinit
			    else
				u_ct_val(1,1,i)=1.5*(i_turb_inlet*ufreestream)**2
				u_ct_val(1,2,i)=(c_mu_inlet**(-0.25))*sqrt(u_ct_val(1,1,ki))&
					/l_turb_inlet*rg(1)
			    end if
			endif
		k=k+prev_turbmodel+lamps
		
	    end do
	end if
	
	
	end if
	
	
	
	
	
	call mpi_barrier(mpi_comm_world, ierror)
	deallocate(dispt,array2)
	
	
      
if (averaging .eq. 1) then
 
if (average_restart.eq.1)then

disp_in_file=0
allocate(dispt(kmaxe),array2(kmaxe*(nof_variables+prev_turbequation+passivescalar+6+passivescalar)))
do i=1,kmaxe
	dispt(i)=(xgo(i)-1)*(nof_variables+prev_turbequation+passivescalar+6+passivescalar)
      end do

      restfile='RESTARTAV.dat'
      n_end=nof_variables+prev_turbequation+passivescalar+6+passivescalar
      
       call mpi_type_create_indexed_block(kmaxe,n_end,dispt,mpi_double_precision,datatype,ierror)
    call mpi_type_commit(datatype,ierror)
    call mpi_file_open(mpi_comm_world, restfile,mpi_mode_rdonly,mpi_info_null, fh, ierror)
	call mpi_file_set_view(fh, disp_in_file, mpi_double_precision,datatype, 'native',mpi_info_null, ierror)
	call mpi_file_read_all(fh, array2, kmaxe*n_end, mpi_double_precision,mpi_status_ignore, ierror)        
        call mpi_file_close(fh, ierror)
	call mpi_type_free(datatype,ierror)
	
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	    k=1
	  do i=1,kmaxe
	      u_c_val(ind1,1:nof_variables,i)=array2(k:k+nof_variables-1)
	      k=k+nof_variables
	      u_ct_val(ind1,1:turbulenceequations+passivescalar,i)=array2(k:k+turbulenceequations+passivescalar-1)
	      k=k+turbulenceequations+passivescalar
	      u_c_rms(1:6+passivescalar,i)=array2(k:k+6+passivescalar-1)
	      k=k+6+passivescalar
	  end do
      else
	  k=1
	  do i=1,kmaxe
	      u_c_val(ind1,1:nof_variables,i)=array2(k:k+nof_variables-1)
	      
	      k=k+nof_variables
	      u_c_rms(1:6+passivescalar,i)=array2(k:k+6+passivescalar-1)
	      
	      k=k+6+passivescalar
	  end do
      end if
	
	call mpi_barrier(mpi_comm_world, ierror)
	deallocate(dispt,array2)
	
    
else
do i=1,kmaxe
      u_c_rms(:,i)=zero
      u_c_val(ind1,:,i)=zero
	if ((passivescalar.gt.0).or.(turbulence.eq.1))then
	u_ct_val(ind1,:,i)=zero
	end if   
end do


end if

end if
!$omp end master













end subroutine rest_read







subroutine checkpointav2d(n)  
!> @brief
!> this subroutine writes the average checkpointing files in 2d
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella
real,allocatable,dimension(:)::valuesa,valuess
real,allocatable,dimension(:,:)::xbin
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj
 character(len=20)::proc,restfile,proc3
 real,allocatable,dimension(:)::igint,tgint
 kmaxe=xmpielrank(n)



icpuid=n
	restfile='RESTARTAV.dat'
	if (n.eq.0)then
	open(1086,file=restfile,form='unformatted',status='replace',action='write')
	
	
	end if
	
	kmaxe=xmpielrank(n)
    
dumg=kmaxe
call mpi_barrier(mpi_comm_world,ierror)

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


if (n.eq.0)then
	allocate(icella(imaxp*isize))
	 icella=0

end if

call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)


call mpi_barrier(mpi_comm_world,ierror)
deallocate (icell)

if (n.eq.0) then
allocate(valuesa(imaxp*isize))
  allocate(xbin(imaxe,(4+turbulenceequations+passivescalar+3+passivescalar)))
	valuesa=0.0

 end if
allocate(valuess(imaxp))
  valuess=0.0

do jj=1,4+turbulenceequations+passivescalar+3+passivescalar
 do i=1,kmaxe
      
      if (jj.le.4+turbulenceequations+passivescalar)then
      if (jj.le.4) then
	  valuess(i)=u_c_val(5,jj,i)
		      
      else

		valuess(i)=u_ct_val(5,jj-turbulenceequations+passivescalar,-i)

      end if
      end if
      if (jj.gt.4+turbulenceequations+passivescalar)then

       valuess(i)=u_c_rms(jj-(4+turbulenceequations+passivescalar),i)

      end if
      
 end do

    call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    do i=1,imaxp*isize
	if (icella(i).gt.0)then
	xbin(icella(i),jj)=valuesa(i)
	end if
    end do
    end if
    
    
 end do
call mpi_barrier(mpi_comm_world,ierror)

    if (n.eq.0)then

    do i=1,imaxe
     write(1086)i
!     do nvar=1,5+turbulenceequations+passivescalar
    write(1086)xbin(xmpi_re(i),1:nof_variables+turbulenceequations+passivescalar+3+passivescalar)
!     end do
    end do
      
    deallocate(xbin,icella,valuesa)

      close(1086)
    end if
    

deallocate(valuess)
	

end subroutine checkpointav2d


subroutine probing
!> @brief
!> this subroutine writes the primitve variables at the probe positions
implicit none
integer::inv
character(len=120)::prob,probfile,proc3
logical::heres
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
                    if (nof_variables.gt.1)then

			if (nprobes.gt.0)then
			    
			    do inv=1,nprobes
			    if (probei(n,inv).ne.0) then
			      write(prob,fmt='(i10)') inv
			      probfile='probe.'//trim(adjustl(prob))

			      inquire (file=probfile,exist=heres)
			    if (heres.eqv..true.) then
				open(3000+n,file=probfile,form='formatted',status='old',action='write',position='append')
				
				
				else
				open(3000+n,file=probfile,form='formatted',status='new',action='write')
				
				end if
				if (passivescalar.eq.0)then
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,probei(n,inv))
				call cons2prim(n,leftv,mp_pinfl,gammal)
	write(3000+n,'(1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,leftv(1),leftv(2),leftv(3),leftv(4),leftv(5)
				else
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,probei(n,inv))
				call cons2prim(n,leftv,mp_pinfl,gammal)
	write(3000+n,'(1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,leftv(1),leftv(2),leftv(3),leftv(4),leftv(5)&
	,u_ct_val(1,1,probei(n,inv))/u_c_val(1,1,probei(n,inv))


				end if
				close(3000+n)
				
		    
			      end if     
			    end do
			  end if
			  end if

end subroutine probing


subroutine probing2d
!> @brief
!> this subroutine writes the primitve variables at the probe positions in 2d
implicit none
integer::inv
character(len=120)::prob,probfile,proc3
logical::heres
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
                        if(nof_variables.gt.1)then

			if (nprobes.gt.0)then
			    
			    do inv=1,nprobes
			    if (probei(n,inv).ne.0) then
			      write(prob,fmt='(i10)') inv
			      probfile='probe.'//trim(adjustl(prob))

			      inquire (file=probfile,exist=heres)
			    if (heres.eqv..true.) then
				open(3000+n,file=probfile,form='formatted',status='old',action='write',position='append')
				
				
				else
				open(3000+n,file=probfile,form='formatted',status='new',action='write')
				
				end if
				if (passivescalar.eq.0)then
	write(3000+n,'(1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,u_c_val(1,1,probei(n,inv)),&
	u_c_val(1,2,probei(n,inv))/u_c_val(1,1,probei(n,inv)),&
	u_c_val(1,3,probei(n,inv))/u_c_val(1,1,probei(n,inv))
				else
	write(3000+n,'(1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,u_c_val(1,1,probei(n,inv)),&
	  u_c_val(1,2,probei(n,inv))/u_c_val(1,1,probei(n,inv)),&
	u_c_val(1,3,probei(n,inv))/u_c_val(1,1,probei(n,inv))&
	,u_ct_val(1,1,probei(n,inv))/u_c_val(1,1,probei(n,inv))


				end if
				close(3000+n)
				
		    
			      end if     
			    end do
			  end if
			  end if

end subroutine probing2d



subroutine computeforce(n)
!> @brief
!> this subroutine computes the forces on the wall
implicit none
integer,intent(in)::n
integer::i,k,j,kmaxe,gqi_points,nnd
integer:: mysurface
character(len=12)::proc,restfile,proc3
real::drag,lift,cd,cl,rx,px,ex,surface_temp,rtemp,fx,fy,fz,mx,my,mz
real::forcexfr,ssx,cdf,liftf,dragf,frictionf,tauyx,tauzx,tauzy,ssy,ssz,cf,tauxx,tauyy,tauzz
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,angle1,angle2,nx,ny,nz
 real,dimension(3)::ci,co
 logical::heref
 integer::im
 real::tsolr,tsolu,tsole,tsolv,tsolw,tsolp,ssp
 real,dimension(1:dims,1:dims)::vortet1
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real,dimension(1:numberofpoints2)::weights_temp
real,dimension(1:4)::viscl,laml
real::fxr,fyr,fzr,mome_xcc,mome_ycc,mome_zcc

!$omp master

forcex=zero; forcey=zero; forcez=zero;  forcexfr=zero
 cd=zero
 cl=zero
 ci(:)=zero
 co(:)=zero
 fx=zero
 fy=zero
 fz=zero
 mx=zero
 my=zero
 mz=zero
 momenty=zero 
 momentz=zero
 momentx=zero
 !$omp end master
!$omp barrier


 kmaxe=xmpielrank(n)
 
#ifdef gpu
!$omp target teams distribute parallel do                              &
!$omp& map(tofrom:forcex, forcey, forcez, momentx, momenty, momentz)    &
!$omp& private(mysurface, j, k, im, nnd, gqi_points,                   &
!$omp& angle1, angle2, nx, ny, nz, ssx, ssy, ssz, ssp, surface_temp,   &
!$omp& ux, uy, uz, vx, vy, vz, wx, wy, wz, px, tauxx, tauyy, tauzz,    &
!$omp& tauyx, tauzx, tauzy, fxr, fyr, fzr, mome_xcc, mome_ycc,         &
!$omp& mome_zcc, vortet1, leftv, rightv, viscl, laml, weights_temp)    &
!$omp& reduction(+:forcex, forcey, forcez, momentx, momenty, momentz)
#else
!$omp do reduction(+:forcex,forcey,forcez,momentx,momenty,momentz)
#endif
do i=1,kmaxe
		if (ielem_interior(i).eq.1)then
			if(mrf.eq.1)then
				mysurface=rec_mrf(i)
			else
				mysurface=1	
			end if	
		if(mysurface.eq.1)then
		    do j=1,ielem_ifca(i)
		      if (ielem_ibounds(j,i).gt.0)then
			  if ((ibound_icode(ielem_ibounds(j,i)).eq.4))then
			      angle1=ielem_faceanglex(j,i)
			      angle2=ielem_faceangley(j,i)
			      nx=(cos(angle1)*sin(angle2))
			      ny=(sin(angle1)*sin(angle2))
			      nz=(cos(angle2))
			      
			  ssx=zero; ssp=zero; ssy=zero; ssz=zero
			  
				select case(ielem_types_faces(j,i))
				case (5)
					  gqi_points=qp_quad
					  weights_temp(1:gqi_points)=weights_q(1:gqi_points)
					  surface_temp=ielem_surf(j,i)
					  
				  
				case(6)
					gqi_points=qp_triangle
					weights_temp(1:gqi_points)=weights_t(1:gqi_points)
 					    surface_temp=ielem_surf(j,i)
 					    
 					    
					  
				end select
				  
				  
				  do im=1,gqi_points
				  
				  if (itestcase.eq.4)then
				  if (ielem_ggs(i).eq.1)then
				  
				  vortet1(1:3,1:3) = rec_grads(1:3,1:3,i)
				  ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
				  vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
				  wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)
				  
				  else
				  
				  vortet1(1,1:3)=rec_uleftv(1:3,2,j,im,i)
				  vortet1(2,1:3)=rec_uleftv(1:3,3,j,im,i)
				  vortet1(3,1:3)=rec_uleftv(1:3,4,j,im,i)
				  ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
				  vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
				  wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)
				  
				  
				  end if
				 end if
				  
				  if (dg.eq.1)then
				  leftv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)
				  rightv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)


				  else
				  leftv(1:nof_variables)=rec_uleft(:,j,im,i)
				  rightv(1:nof_variables)=rec_uleft(:,j,im,i)
				  end if
				  
				  
				  
				    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
				    px=leftv(5)
				    ssp=ssp+(px*weights_temp(im))
				    if (itestcase.eq.4)then
				    if (dg.eq.1)then
				  leftv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)
				  rightv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)


				  else
				  leftv(1:nof_variables)=rec_uleft(:,j,im,i)
				  rightv(1:nof_variables)=rec_uleft(:,j,im,i)
				  end if
				    call get_visc_conduct(n,leftv,rightv,viscl,laml)
				  
				  
				  tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
				  tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
				  tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
				  tauyx=(uy + vx)
				  tauzx=(wx + uz)
				  tauzy=(vz + wy)
				  ssx=ssx-((viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))*weights_temp(im))
				  ssy=ssy-((viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))*weights_temp(im))
				  ssz=ssz-((viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))*weights_temp(im))
				 end if
				   end do
				   
				   
				  ssp=ssp-pres	
				  
				  forcex=forcex+(((ssp)*(surface_temp)*nx))+((ssx)*surface_temp)
				  forcey=forcey+(((ssp)*(surface_temp)*ny))+((ssy)*surface_temp)
				  forcez=forcez+(((ssp)*(surface_temp)*nz))+((ssz)*surface_temp)

				  fxr=(((ssp)*(surface_temp)*nx))+((ssx)*surface_temp)
				  fyr=(((ssp)*(surface_temp)*ny))+((ssy)*surface_temp)
				  fzr=(((ssp)*(surface_temp)*nz))+((ssz)*surface_temp)

				  mome_xcc=ielem_xxc(i)-origin(1)
				  mome_ycc=ielem_yyc(i)-origin(2)
				  mome_zcc=ielem_zzc(i)-origin(3)


				  momentx=momentx+(fzr*mome_ycc)-(fyr*mome_zcc)
				  momenty=momenty+(fxr*mome_zcc)-(fzr*mome_xcc)
				  momentz=momentz+(fyr*mome_xcc)-(fxr*mome_ycc)

					    
			end if
		      end if
		    end do
		end if !mysurface
		end if


		
end do					 
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	

		
	
	

!$omp master
	forcex=forcex*vectorx
	forcey=forcey*vectory
	forcez=forcez*vectorz
	if(rframe.eq.0)then
        rtemp=((aoa/180.0d0)*pi)
	liftf=(forcez*cos(rtemp))+(forcey*cos(rtemp))-(forcex*sin(rtemp))
	dragf=(forcex*cos(rtemp))+(forcey*sin(rtemp))+(forcez*sin(rtemp))
	cl=(2.0d0*liftf)/((rres)*(ufreestream**2))
	cd=(2.0d0*dragf)/((rres)*(ufreestream**2))


	co(1)=cl
	co(2)=cd
	call mpi_allreduce(co(1:2),ci(1:2),2,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
	cl=ci(1)
	cd=ci(2)
		co(1)=momentx
        co(2)=momenty
        co(3)=momentz
        call mpi_allreduce(co(1:3),ci(1:3),3,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        mx=ci(1)
        my=ci(2)
        mz=ci(3)



	else
        co(1)=forcex
        co(2)=forcey
        co(3)=forcez
        call mpi_allreduce(co(1:3),ci(1:3),3,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        fx=ci(1)
        fy=ci(2)
        fz=ci(3)
		co(1)=momentx
        co(2)=momenty
        co(3)=momentz
        call mpi_allreduce(co(1:3),ci(1:3),3,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        mx=ci(1)
        my=ci(2)
        mz=ci(3)
	end if
	if (n.eq.0) then
	inquire (file='FORCE.dat',exist=heref)
		if (heref) then
		open(50+n,file='FORCE.dat',form='formatted',status='old',action='write',position='append')
		else
	open(50+n,file='FORCE.dat',form='formatted',status='new',action='write')
		end if
	inquire (file='MOMENT.dat',exist=heref)
		if (heref) then
		open(500+n,file='MOMENT.dat',form='formatted',status='old',action='write',position='append')
		else
	open(500+n,file='MOMENT.dat',form='formatted',status='new',action='write')
		end if	
	if(rframe.eq.0)then	
	write(50+n,'(i14,1x,e14.7,1x,e14.7,1x,e14.7)')it,t,cl,cd
	write(500+n,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,t,mx,my,mz
	else
	write(50+n,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,t,fx,fy,fz
	write(500+n,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,t,mx,my,mz
	end if
	close(50+n)
	close(500+n)
	end if		

	
	
!$omp end master
!$omp barrier
	
	
	
	
	
	

end subroutine computeforce


subroutine computeforce2d(n)
!> @brief
!> this subroutine computes the forces on the wall in 2d
implicit none
integer,intent(in)::n
integer::i,k,j,kmaxe,gqi_points,nnd
character(len=12)::proc,restfile,proc3
real::drag,lift,cd,cl,rx,px,ex,surface_temp,rtemp
real::forcexfr,ssx,cdf,liftf,dragf,frictionf,tauyx,tauzx,tauzy,ssy,ssz,cf,tauxx,tauyy,tauzz
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,nx,ny,nz,angle1,angle2
 real,dimension(2)::ci,co
 logical::heref
 integer::im
 real::tsolr,tsolu,tsole,tsolv,tsolw,tsolp,ssp
 real,dimension(1:dims,1:dims)::vortet1
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real,dimension(1:numberofpoints2)::weights_temp
real,dimension(1:4)::viscl,laml


!$omp master
forcex=zero; forcey=zero; forcez=zero;  forcexfr=zero
 cd=zero
 cl=zero
 ci(:)=zero
 co(:)=zero

 !$omp end master
!$omp barrier


 kmaxe=xmpielrank(n)
 
#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:forcex,forcey) &
!$omp& private(j, k, im, nnd, gqi_points, angle1, angle2, nx, ny, nz, &
!$omp&         ssx, ssy, ssz, ssp, surface_temp, ux, uy, uz, vx, vy, vz, &
!$omp&         wx, wy, wz, px, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy, &
!$omp&         vortet1, leftv, rightv, viscl, laml, weights_temp) &
!$omp& reduction(+:forcex,forcey)
#else
!$omp do reduction(+:forcex,forcey)
#endif
do i=1,kmaxe
		if (ielem_interior(i).eq.1)then	
		    do j=1,ielem_ifca(i)
		      if (ielem_ibounds(j,i).gt.0)then
			  if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
			      nx=ielem_faceanglex(j,i)
			      ny=ielem_faceangley(j,i)
			      
			      
			  ssx=zero; ssp=zero; ssy=zero; 
			  
				
					  gqi_points=qp_line_n
					  weights_temp(1:qp_line) = weights_l(1:qp_line)
					  surface_temp=ielem_surf(j,i)
					  
				  
							  
				  
				  do im=1,gqi_points
				  
				  if (itestcase.eq.4)then
				  if (ielem_ggs(i).eq.1)then
				  
				  vortet1(1:2,1:2) = rec_grads(1:2,1:2,i)
				  ux = vortet1(1,1);uy = vortet1(1,2)
				  vx = vortet1(2,1);vy = vortet1(2,2)
				
				  
				  else
				  
				  vortet1(1,1:2)=rec_uleftv(1:2,2,j,im,i)
				  vortet1(2,1:2)=rec_uleftv(1:2,3,j,im,i)
				  
				 ux = vortet1(1,1);uy = vortet1(1,2)
				  vx = vortet1(2,1);vy = vortet1(2,2)
				
				  
				  
				  end if
				  end if
				  
				  leftv(1:nof_variables)=rec_uleft(:,j,im,i)
				  rightv(1:nof_variables)=rec_uleft(:,j,im,i)
				    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
				    px=leftv(4)
				    
				    
				    if (itestcase.eq.4)then
				    leftv(1:nof_variables)=rec_uleft(:,j,im,i)
				  rightv(1:nof_variables)=rec_uleft(:,j,im,i)
				    call get_visc_conduct(n,leftv,rightv,viscl,laml)
				  
				  
				  tauxx=2.0d0*ux
				  tauyy=2.0d0*vy
				  tauyx=(uy + vx)
				  
				  ssx=ssx-((viscl(1)*((ny*tauyx)))*weights_temp(im))
				  ssy=ssy-((viscl(1)*((nx*tauyx)))*weights_temp(im))
				  end if
				  ssp=ssp+(px*weights_temp(im))
				   end do
				   
				   
				  ssp=ssp-pres	
				  
				  forcex=forcex+(((ssp)*(surface_temp)*nx))+((ssx)*surface_temp)
				  forcey=forcey+(((ssp)*(surface_temp)*ny))+((ssy)*surface_temp)
				  

					    
					    
			end if
		      end if
		    end do
		end if

			
		
end do					 
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	

		
	
	

!$omp master
	forcex=forcex*vectorx
	forcey=forcey*vectory
	
        rtemp=((aoa/180.0d0)*pi)
	liftf=(forcey*cos(rtemp))-(forcex*sin(rtemp))
	dragf=(forcex*cos(rtemp))+(forcey*sin(rtemp))
	cl=(2.0d0*liftf)/((rres)*(ufreestream**2))
	cd=(2.0d0*dragf)/((rres)*(ufreestream**2))


	co(1)=cl
	co(2)=cd
	call mpi_allreduce(co(1:2),ci(1:2),2,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
	cl=ci(1)
	cd=ci(2)
	
	if (n.eq.0) then
	inquire (file='FORCE.dat',exist=heref)
		if (heref) then
		open(50+n,file='FORCE.dat',form='formatted',status='old',action='write',position='append')
		else
	open(50+n,file='FORCE.dat',form='formatted',status='new',action='write')
		end if
	write(50+n,'(i14,1x,e14.7,1x,e14.7,1x,e14.7)')it,t,cl,cd
	
	close(50+n)
	end if		
	call mpi_barrier(mpi_comm_world,ierror)
	
	
!$omp end master
!$omp barrier
	
	
	
	
	

end subroutine computeforce2d


subroutine calculate_residual(n)
!> @brief
!> this subroutine computes and writes the residual for 3d
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::suml3,dum_resi

kmaxe=xmpielrank(n)

!$omp master
allres(:)=zero

!$omp end master
!$omp barrier


if ((itestcase.le.4).and.(turbulence.ne.1))then
#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:allres(1:nof_variables)) &
!$omp& reduction(+:allres(1:nof_variables))
#else
!$omp do reduction(+:allres(1:nof_variables))
#endif
do i=1,kmaxe
	if (dg.eq.1)then
    allres(1:nof_variables)=allres(1:nof_variables)+((rhs_valdg(1,1:nof_variables,i)*ielem_totvolume(i))**2)
    else
    allres(1:nof_variables)=allres(1:nof_variables)+((rhs_val(1:nof_variables,i)*ielem_totvolume(i))**2)

    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

!$omp master
do i=1,5
suml3=allres(i)
dum_resi=zero
call mpi_allreduce(suml3,dum_resi,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allres(i)=sqrt(dum_resi/totalvolume)

end do


do i=1,5
if (initialres(i).le.allres(i))then
initialres(i)=allres(i)
end if
allres(i)=allres(i)/initialres(i)

end do
!$omp end master




end if

if (turbulence.eq.1)then
#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:allres(1:nof_variables)) &
!$omp& reduction(+:allres(1:nof_variables))
#else
!$omp do reduction(+:allres(1:nof_variables))
#endif
do i=1,kmaxe
    allres(1:nof_variables)=allres(1:nof_variables)+((rhs_val(1:nof_variables,i)*ielem_totvolume(i))**2)
    allres(nof_variables+1:nof_variables+turbulenceequations)=allres(nof_variables+1:nof_variables+turbulenceequations)+((rhst_val(1:turbulenceequations,i)*ielem_totvolume(i))**2)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

!$omp master
do i=1,7
suml3=allres(i)
dum_resi=zero
call mpi_allreduce(suml3,dum_resi,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allres(i)=sqrt(dum_resi/totalvolume)

end do



do i=1,7
if (initialres(i).le.allres(i))then
initialres(i)=allres(i)
end if
allres(i)=allres(i)/initialres(i)

end do

if (turbulenceequations.eq.1) allres(7)=1.0d0
!$omp end master

end if



!$omp master
if (n.eq.0)then
if ((itestcase.le.4).and.(turbulence.ne.1))then

	    open(67,file='residual.dat',form='formatted',action='write',position='append')
	    write(67,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,allres(1),allres(2),allres(3),allres(4),allres(5)
	    close(67)

else

	    open(67,file='residual.dat',form='formatted',action='write',position='append')
	    write(67,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,allres(1),allres(2),allres(3),allres(4),allres(5),allres(6),allres(7)
	    close(67)



end if

end if 


if ((allres(1).lt.reslimit).and.(allres(2).lt.reslimit).and.(allres(3).lt.reslimit).and.(allres(4).lt.reslimit).and.(allres(5).lt.reslimit))then
 kill=1
 end if



!$omp end master

 









end subroutine


subroutine calculate_residual2d(n)
!> @brief
!> this subroutine computes and writes the residual for 2d
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::suml3,dum_resi
kmaxe=xmpielrank(n)


!$omp master
allres(:)=zero
!$omp end master
!$omp barrier


if ((itestcase.le.4).and.(turbulence.ne.1))then
#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:allres(1:nof_variables)) &
!$omp& reduction(+:allres(1:nof_variables))
#else
!$omp do reduction(+:allres(1:nof_variables))
#endif
do i=1,kmaxe

    if (dg.eq.1)then
    allres(1:nof_variables)=allres(1:nof_variables)+((rhs_val(1:nof_variables,i)*ielem_totvolume(i))**2)
    else
    allres(1:nof_variables)=allres(1:nof_variables)+((rhs_val(1:nof_variables,i)*ielem_totvolume(i))**2)

    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




!$omp master

do i=1,4
suml3=allres(i)
dum_resi=zero
call mpi_allreduce(suml3,dum_resi,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allres(i)=sqrt(dum_resi/totalvolume)

end do


do i=1,4
if (initialres(i).le.allres(i))then
initialres(i)=allres(i)
end if
allres(i)=allres(i)/initialres(i)

end do

!$omp end master




end if

if (turbulence.eq.1)then
#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:allres(1:nof_variables)) &
!$omp& reduction(+:allres(1:nof_variables))
#else
!$omp do reduction(+:allres(1:nof_variables))
#endif
do i=1,kmaxe
    allres(1:nof_variables)=allres(1:nof_variables)+((rhs_val(1:nof_variables,i)*ielem_totvolume(i))**2)
    allres(nof_variables+1:nof_variables+turbulenceequations)=allres(nof_variables+1:nof_variables+turbulenceequations)+((rhst_val(1:turbulenceequations,i)*ielem_totvolume(i))**2)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

!$omp master
do i=1,nof_variables+turbulenceequations
suml3=allres(i)
dum_resi=zero
call mpi_allreduce(suml3,dum_resi,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allres(i)=sqrt(dum_resi/totalvolume)

end do



do i=1,nof_variables+turbulenceequations
 if (initialres(i).le.allres(i))then
initialres(i)=allres(i)
end if
allres(i)=allres(i)/initialres(i)

end do
!$omp end master

end if



!$omp master
if (n.eq.0)then
if ((itestcase.le.4).and.(turbulence.ne.1))then

	    open(67,file='residual.dat',form='formatted',action='write',position='append')
	    write(67,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,allres(1),allres(2),allres(3),allres(4)
	    close(67)

else

	    open(67,file='residual.dat',form='formatted',action='write',position='append')
	    write(67,'(i14,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')it,allres(1),allres(2),allres(3),allres(4),allres(5)
	    close(67)



end if

end if 


 if ((allres(1).lt.reslimit).and.(allres(2).lt.reslimit).and.(allres(3).lt.reslimit).and.(allres(4).lt.reslimit))then
 kill=1
 end if

!$omp end master







end subroutine

subroutine calculate_error(n)
!> @brief
!> this subroutine computes and writes the l2,linfinity or l1 norm
	implicit none
	integer,intent(in)::n
	integer::i,k,kmaxe,ind_er
	real::exact,dummyout,dummyin
	real::aproximate
	real,dimension(15)::condm
	kmaxe=xmpielrank(n)
	l0norm=zero;stennorm=zero;l1norm=zero
	ind_er=1
	if (multispecies.eq.1)then
	ind_er=nof_variables
	end if
	  !$omp master
	  call mpi_barrier(mpi_comm_world,ierror)
	  !$omp end master
	  !$omp barrier
	  
			!$omp do reduction (+:l1norm)
			do i=1,kmaxe
				if (itestcase.le.3)then
				
				exact=u_e_val(1,ind_er,i)
				
				aproximate=u_c_val(1,ind_er,i)
				
! 					if ((abs(aproximate-exact)).gt.l0norm(n,1))then
! 					l0norm(n,1)=abs(aproximate-exact)
! 					end if
					l1norm=l1norm+((aproximate-exact)**2)*ielem_totvolume(i)
				end if
 			end do
 			!$omp end do 
 			
 			if (initcond.eq.0)then
 			!$omp do reduction (+:l0norm)
			do i=1,kmaxe
				if (itestcase.le.3)then
! 				condm(2)=rec_cond(2,i)
! 				l0norm=rec_cond(1,i)
! 					if (maxval(condm).gt.l0norm)then
! 					l0norm=l0norm+abs(rec_cond(1,i))
! 					end if
! 					l1norm(n,1)=l1norm(n,1)+((abs(aproximate-exact)))
				end if
 			end do
 			!$omp end do 
 			else
 			!$omp do reduction (max:l0norm)
			do i=1,kmaxe
				if (itestcase.le.3)then
				exact=u_e_val(1,ind_er,i)
				
				
				aproximate=u_c_val(1,ind_er,i)
				
					if ((abs(aproximate-exact)).gt.l0norm)then
					l0norm=abs(aproximate-exact)
					end if
! 					l1norm(n,1)=l1norm(n,1)+((abs(aproximate-exact)))
				end if
 			end do
 			!$omp end do 
 			
 			
 			
 			end if
 			if (initcond.eq.3)then
 			l0norm=zero;l1norm=tolbig
 			
 			!$omp do reduction (max:l0norm)
			do i=1,kmaxe
					
					if (u_c_val(1,ind_er,i).gt.l0norm)then
					l0norm=u_c_val(1,ind_er,i)
					end if
! 					l1norm(n,1)=l1norm(n,1)+((abs(aproximate-exact)))
				
 			end do
 			!$omp end do 
 			!$omp do reduction (min:l1norm)
			do i=1,kmaxe
					
					if (u_c_val(1,ind_er,i).lt.l1norm)then
					l1norm=u_c_val(1,ind_er,i)
					end if
! 					l1norm(n,1)=l1norm(n,1)+((abs(aproximate-exact)))
				
 			end do
 			!$omp end do 
 			
 			
 			end if
 			
 			
 			if (initcond.eq.0)then
 			!$omp do reduction (+:stennorm)
 			do i=1,kmaxe
				
! 				stennorm=stennorm+abs(rec_cond(2,i))
 			end do
 			!$omp end do 
 			
 			else
 			!$omp do reduction (+:stennorm)
 			do i=1,kmaxe
				stennorm=stennorm+u_c_val(1,ind_er,i)
 			end do
 			!$omp end do 
 			end if
			
 			!$omp master
			call mpi_barrier(mpi_comm_world,ierror)
			!$omp end master
			!$omp barrier
 			
 			
 			!$omp master
 			if (initcond.eq.3)then
 			dummyout=l0norm
 			call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
 			l0norm=dummyin
 			dummyout=l1norm
 			call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
 			l1norm=dummyin
 			else
 			
 			dummyout=l1norm
 			call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
 			l1norm=dummyin
 			dummyout=l0norm
 			if (initcond.eq.0)then
 			call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
 			dummyin=dummyin/imaxe
 			else
 			
 			call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
 			end if
 			l0norm=dummyin
 			dummyout=stennorm
 			call mpi_allreduce(dummyout,dummyin,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
 			stennorm=dummyin
 			
 			
 			
 			
 			end if

 			
 			cpux3(1) = mpi_wtime()
 			if (n.eq.0)then
			open(30,file='errors.dat',form='formatted',action='write',position='append')
				if (initcond.eq.1)then
				write(30,'(i9,1x,e14.7,1x,i4,1x,e14.7,1x,e14.7)')imaxe,t,spatiladiscret,l0norm,stennorm/imaxe

				else
					if (initcond.ne.3)then

					write(30,'(i9,1x,i4,1x,i4,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')imaxe,iorder,spatiladiscret,l0norm,sqrt(l1norm/totalvolume),stennorm/imaxe,(cpux3(1)-cpux2(1))*isize
					else
					write(30,'(i9,1x,i4,1x,i4,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')imaxe,iorder,spatiladiscret,l0norm,l1norm,stennorm/imaxe,(cpux3(1)-cpux2(1))*isize

					end if
				end if
! 			write(30,'(i9,1x,i4,1x,i4,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')imaxe,iorder,spatiladiscret,l0norm,sqrt(l1norm/totalvolume),stennorm/imaxe,(cpux3(1)-cpux2(1))*isize
			close(30)
			end if
			!$omp end master
			!$omp barrier


			!$omp master
			call mpi_barrier(mpi_comm_world,ierror)
			!$omp end master
			!$omp barrier



end subroutine calculate_error



subroutine outwritepara3d
!> @brief
!> this subroutine writes the solution and the grid file in ascii vtk format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 

integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml


 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(12))

kmaxe=xmpielrank(n)
nvar1=2





if (n.eq.0)then

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') 
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
end if



if (n.eq.0)then

open(400+n,file=outfile,form='formatted',status='new',action='write')
write(400+n,'(a)')"# vtk datafile version 3.0"
write(400+n,'(a)')"vtk output"
write(400+n,'(a)')"ascii"
write(400+n,'(a)')"dataset unstructured_grid"
write(400+n,'(a)')"field fielddata 1"
write(400+n,'(a)')"time 1 1 double"
write(400+n,*) t
write(400+n,'(a6,2x,i10,2x,a6)')"points",imaxn,"double"



allocate (valuelocation(nvar1))


valuelocation(:)=0
valuelocation(1:2)=1

    
	if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
        do i=1,imaxn
	read(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n,'(2x,g14.6,2x,g14.6,2x,g14.6)')x,y,z
	end do
	close(96)
	else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
        do i=1,imaxn
	read(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n,'(2x,g14.6,2x,g14.6,2x,g14.6)')x,y,z
	end do
	close(96)
	end if

   
                    
                    
 write(400+n,*)                   
write(400+n,'(a5,2x,i10,2x,i10)')"cells",imaxe,(imaxe*8)+imaxe

if (binio.eq.0)then
open(97,file='GRID.cel',form='formatted',status='old',action='read')
do i=1,imaxe
read(97,*)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n,'(9i12)')8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)

else
open(97,file='GRID.cel',form='unformatted',status='old',action='read')
do i=1,imaxe
read(97)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n,'(9i12)')8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)
end if
                    
                    
write(400+n,*)                    
write(400+n,'(a10,2x,i10)')"cell_types",imaxe
do i=1,imaxe
write(400+n,*)12
end do
write(400+n,*)
write(400+n,'(a9,2x,i10)')"cell_data",imaxe
               
                    
                    

 
  allocate(xbin(imaxe),xbin2(imaxe))
	

 end if
 allocate(valuess(kmaxe))
 
  call mpi_barrier(mpi_comm_world,ierror)
  
    
   
do j=1,nof_variables
     
     if ((j.ge.2).and.(j.le.4))then
     
        do i=1,kmaxe
        valuess(i)=u_c_val(1,j,i)/u_c_val(1,1,i)!0.0
        end do
    end if
    if (j.eq.1)then
        do i=1,kmaxe
        valuess(i)=u_c_val(1,j,i)
        end do
    end if
    if (j.eq.5)then
                do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		end do
    
    end if
    
		
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)


    if (n.eq.0)then

    if (j.eq.1)then
    write(400+n,'(a)')"scalars  r double 1"
    end if
    if (j.eq.2)then
    write(400+n,*)
    write(400+n,'(a)')"scalars  u double 1"
    end if
    if (j.eq.3)then
    write(400+n,*)
    write(400+n,'(a)')"scalars  v double 1"
    end if
    if (j.eq.4)then
    write(400+n,*)
    write(400+n,'(a)')"scalars  w double 1"
    end if
    if (j.eq.5)then
    write(400+n,*)
    write(400+n,'(a)')"scalars  p  double 1"
    end if
    write(400+n,'(a)')"lookup_table default"
    
    
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
    
     write(400+n,*)xbin(1:imaxe)
    end if

    
end do

    
    
    
    

     
    
    
    
 if (n.eq.0)then
  close(400+n) 
   
   deallocate(xbin,valuesa,xbin2,valuelocation,icella)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)

        
 deallocate(variables)


	
	

end subroutine outwritepara3d

subroutine outwritepara3db
!> @brief
!> this subroutine writes the solution and the grid file in binary vtk format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(12))
nvar1=2
kmaxe=xmpielrank(n)




if (n.eq.0)then

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') 
	!proc4=".plt"
	outfile="OUT_"//trim(adjustl(proc3))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
end if

lf=char(10)

if (n.eq.0)then

open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') imaxn
write(400+n)"points "//str1//" double"//lf

allocate (valuelocation(nvar1))


valuelocation(:)=0
valuelocation(1:2)=1

    
     
! 	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
!         do i=1,imaxn
! 	read(96,*)j,x,y,z
! 	xbin(i)=x/scaler
! 	ybin(i)=y/scaler
!  	zbin(i)=z/scaler
! 	end do
! 
!     close(96)
    
    if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
        do i=1,imaxn
	read(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n)x,y,z
	end do
	close(96)
	else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
        do i=1,imaxn
	read(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n)x,y,z
	end do
	close(96)
	end if

    
    
   
    
                    
write(str1(1:15),'(i15)') imaxe
write(str2(1:15),'(i15)') (imaxe*8)+imaxe
write(400+n)"cells",str1//str2//lf
if (binio.eq.0)then
open(97,file='GRID.cel',form='formatted',status='old',action='read')
do i=1,imaxe
read(97,*)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)
else
open(97,file='GRID.cel',form='unformatted',status='old',action='read')
do i=1,imaxe
read(97)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)

end if
                    
                    
write(str1(1:15),'(i15)') imaxe                    
write(400+n)"cell_types"//str1//lf	
do i=1,imaxe
write(400+n)12
end do

write(400+n)"cell_data"//str1//lf
!                
!                     
!                     
! 
 
  allocate(xbin(imaxe),xbin2(imaxe))
	
! 
 end if
! 
!  
  
  allocate(valuess(kmaxe))
 
 call mpi_barrier(mpi_comm_world,ierror)
!     
   
do j=1,nof_variables
     
         do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(j)
		end do
    
     
		
    
    
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    
    if (j.eq.1)then
    write(400+n)"scalars  r double 1"//lf
    end if
    if (j.eq.2)then
    write(400+n)
    write(400+n)"scalars  u double 1"//lf
    end if
    if (j.eq.3)then
    write(400+n)
    write(400+n)"scalars  v double 1"//lf
    end if
    if (j.eq.4)then
    write(400+n)
    write(400+n)"scalars  w double 1"//lf
    end if
    if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  p  double 1"//lf
    end if
    if (j.eq.6)then
    write(400+n)
    write(400+n)"scalars  species1  double 1"//lf
    end if
    if (j.eq.7)then
    write(400+n)
    write(400+n)"scalars  species2  double 1"//lf
    end if
    if (j.eq.8)then
    write(400+n)
    write(400+n)"scalars  volumef1  double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf
    
    
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
    
    
     write(400+n)xbin(1:imaxe)
    end if

     
    
    
end do
    
    if (itestcase.eq.4)then

	 do i=1,kmaxe
		  
		  valuess(i)=ielem_vortex(1,i)
		end do
    
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

	if (n.eq.0)then
    
    	write(400+n)
    	write(400+n)"scalars  q double 1"//lf
	    write(400+n)"lookup_table default"//lf
    
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
     write(400+n)xbin(1:imaxe)
    end if
   end if
   
   
   if (turbulence.eq.1)then

	 do i=1,kmaxe
		  
		  valuess(i)=u_ct_val(1,1,i)
		end do
    
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

	if (n.eq.0)then
    
    	write(400+n)
    	write(400+n)"scalars  nut double 1"//lf
	    write(400+n)"lookup_table default"//lf
    
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
     write(400+n)xbin(1:imaxe)
    end if
   end if
    
    

     
    
  
    
  if (n.eq.0)then
  close(400+n) 
   
   deallocate(xbin,valuelocation,xbin2)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)

        
 deallocate(variables)


	
	

end subroutine outwritepara3db




subroutine movie_para
!> @brief
!> this subroutine writes the solution and the grid file in binary vtk format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

!
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2

      integer::   debug,iii,npts,nelm


      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(12))
nvar1=2
kmaxe=xmpielrank(n)




if (n.eq.0)then

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)')
	!proc4=".plt"
	outfile="mov_"//trim(adjustl(proc3))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
end if

lf=char(10)

if (n.eq.0)then

open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') imaxn
write(400+n)"points "//str1//" double"//lf

allocate (valuelocation(nvar1))


valuelocation(:)=0
valuelocation(1:2)=1



! 	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
!         do i=1,imaxn
! 	read(96,*)j,x,y,z
! 	xbin(i)=x/scaler
! 	ybin(i)=y/scaler
!  	zbin(i)=z/scaler
! 	end do
!
!     close(96)

    if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
        do i=1,imaxn
	read(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n)x,y,z
	end do
	close(96)
	else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
        do i=1,imaxn
	read(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n)x,y,z
	end do
	close(96)
	end if






write(str1(1:15),'(i15)') imaxe
write(str2(1:15),'(i15)') (imaxe*8)+imaxe
write(400+n)"cells",str1//str2//lf
if (binio.eq.0)then
open(97,file='GRID.cel',form='formatted',status='old',action='read')
do i=1,imaxe
read(97,*)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)
else
open(97,file='GRID.cel',form='unformatted',status='old',action='read')
do i=1,imaxe
read(97)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)

end if


write(str1(1:15),'(i15)') imaxe
write(400+n)"cell_types"//str1//lf
do i=1,imaxe
write(400+n)12
end do

write(400+n)"cell_data"//str1//lf
!
!
!
!

  allocate(xbin(imaxe),xbin2(imaxe))

!
 end if
!
!

  allocate(valuess(kmaxe))

 call mpi_barrier(mpi_comm_world,ierror)
!

do j=1,1

         do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)

		  valuess(i)=sqrt(leftv(2)**2+leftv(3)**2+leftv(4)**2)
		end do





    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then

    if (j.eq.1)then
    write(400+n)"scalars  vel double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf


		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do


     write(400+n)xbin(1:imaxe)
    end if




end do

    if (itestcase.eq.4)then

	 do i=1,kmaxe

		  valuess(i)=ielem_vortex(1,i)
		end do

    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

	if (n.eq.0)then

    	write(400+n)
    	write(400+n)"scalars  q double 1"//lf
	    write(400+n)"lookup_table default"//lf

		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
     write(400+n)xbin(1:imaxe)
    end if
   end if


   if (turbulence.eq.1)then

	 do i=1,kmaxe

		  valuess(i)=u_ct_val(1,1,i)
		end do

    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

	if (n.eq.0)then

    	write(400+n)
    	write(400+n)"scalars  nut double 1"//lf
	    write(400+n)"lookup_table default"//lf

		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
     write(400+n)xbin(1:imaxe)
    end if
   end if







  if (n.eq.0)then
  close(400+n)

   deallocate(xbin,valuelocation,xbin2)
  deallocate(out1)
  end if
  deallocate (valuess)


  call mpi_barrier(mpi_comm_world,ierror)


 deallocate(variables)





end subroutine movie_para


 !added by holger foysi
                          !!-----------------------------------------------------------------------------------
                          subroutine outwritepara2db
                          !!-----------------------------------------------------------------------------------
                            !!-----------------------------------------------------------------------------------
                            !> @brief
                            !> this subroutine writes the 2d solution including the grid in binary vtk format
                            !!-----------------------------------------------------------------------------------
                            use iso_c_binding
                            implicit none
                            real,dimension(1:nof_variables)::leftv
							real::mp_pinfl,gammal
							real,dimension(1:nof_variables)::rightv
							real::mp_pinfr,gammar
							real::angle1,angle2,nx,ny,nz
							real,dimension(1:4)::viscl,laml
                            integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
                            real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
                            real,allocatable,dimension(:)::ifint,tfint,ndr,nds
                            integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
                            real,allocatable,dimension(:)::variables
                            real,dimension(3,3)::avort,tvort,svort,ovort
                            integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
                            logical::herev
                            real,dimension(5)::total
                            character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
                            integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
                            real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
                            real,allocatable,dimension(:,:)::fbin
                            integer,allocatable,dimension(:,:)::icon
                            integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
                            real,allocatable,dimension(:)::valuess,valuesa
                            character(len=:),allocatable::out1
                            character*1 nulchar
                            character(len=1)   :: flui,lf
                            character(len=15)  :: str1,str2,str1a,str2a
                            integer::   debug,iii,npts,nelm
                            real::    soltime
                            integer:: visdouble, filetype, count
                            integer:: zonetype,strandid,parentzn,isblock,strl
                            integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
                            pointer   (nullptr,null)
                            integer:: null(*)
                            character(len=:),allocatable  :: str_imaxn,str_imaxe,str_imaxe2,str_time
                            character(1) :: c
                            allocate(variables(12))
                            nvar1=2
                            kmaxe=xmpielrank(n)
                            
                            !the previous use of field didn't work, so time is added to the title

                            !create filename
                            if (n.eq.0)then
                               write(proc3,fmt='(i10)') it
                               write(proc5,fmt='(i10)') 
                               outfile="OUT_"//trim(adjustl(proc3))//".vtk"
                               itgfd=len_trim(outfile)
                               allocate(character(len=itgfd) ::out1)
                               out1=outfile(1:itgfd)
                            end if

                            ! char(10)  = linefeed ! alternativeley use ",new_line(c)"
                            lf=char(10)

                            if (n.eq.0)then
                               !file
                               open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
                               !header
                               write(400+n)"# vtk datafile version 3.0",new_line(c)
                               !create dynamic string with size equal to the number of digits of time
                               write(str1,'(f12.6)') t
                               str_time=trim(adjustl(str1))
                               write(400+n)"ucns3d vtk output 2d. time: "//str_time,new_line(c)
                               write(400+n)"binary",new_line(c)
                               write(400+n)"dataset unstructured_grid",new_line(c)

                               !create dynamic string with size equal to the number of digits of the integer
                               write(str1,'(i0)') imaxn
                               str_imaxn=trim(adjustl(str1))

                               !number of node points
                               write(400+n)"points "//str_imaxn//" double",new_line(c)

                               if (binio.eq.0)then !ascii
                                  open(96,file='GRID.vrt',form='formatted',status='old',action='read')
                                  do i=1,imaxn
                                     read(96,*)j,x,y
                                     x=x/scaler;y=y/scaler
                                     write(400+n) x,y,0.d0  !for 2d vtk needs x,y,z, too, just with z=0
                                  end do
                                  close(96)
                               else !binary
                                  open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
                                  do i=1,imaxn
                                     read(96)j,x,y
                                     x=x/scaler;y=y/scaler
                                     write(400+n) x,y,0.d0  !for 2d vtk needs x,y,z, too, just with z=0
                                  end do
                                  close(96)
                               end if

                               !create dynamic strings with size equal to the number of digits of the integer
                               !# of cells
                               write(str1,'(i0)') imaxe
                               str_imaxe  = trim(adjustl(str1))

                               !cell list size:the total number of integer values required to represent the list
                               !changes needed if grids with different cell types are used, here every cell has the same type
                               write(str2,'(i0)') imaxe*5 ! 
                               str_imaxe2 = trim(adjustl(str2))        
                               
                               write(400+n) "cells "//str_imaxe//" "//str_imaxe2,new_line(c)
                               !the grid is such, that triangles are defined using quads with two points being equal (j3=j4)
                               if (binio.eq.0)then
                                  open(97,file='GRID.cel',form='formatted',status='old',action='read')
                                  do i=1,imaxe
                                     read(97,*) j,j1,j2,j3,j4
                                     write(400+n) 4,j1-1,j2-1,j3-1,j4-1
                                  end do
                                  close(97)
                               else
                                  open(97,file='GRID.cel',form='unformatted',status='old',action='read')
                                  do i=1,imaxe
                                     read(97) j,j1,j2,j3,j4
                                     write(400+n) 4,j1-1,j2-1,j3-1,j4-1
                                  end do
                                  close(97)
                               end if

                               !specify type of cells for each cell
                               !changes needed if grids with different cell types are used, here every cell has the same type
                               !hexahedron in 3d (12), quad in 2d (9), triangles (5) etc., see
                               !https://lorensen.github.io/vtkexamples/site/vtkfileformats/
                               write(400+n)"cell_types "//str_imaxe,new_line(c)
                               do i=1,imaxe
                                  write(400+n) 9  
                               end do

                               ! the data is cell centered. use celldatatopointdata in paraview for warp by scalar
                               write(400+n)"cell_data "//str_imaxe,new_line(c)

                               allocate(xbin(imaxe),xbin2(imaxe))
                               ! 
                            end if

                            allocate(valuess(kmaxe))
                            call mpi_barrier(mpi_comm_world,ierror)

                            do j=1,nof_variables
                                
                              
                                  do i=1,kmaxe
                                        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
                                        call cons2prim(n,leftv,mp_pinfl,gammal)
                                        valuess(i)=leftv(j)
                                  end do

                               call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset, &
                                                mpi_double_precision,0,mpi_comm_world,ierror)

                               if (n.eq.0)then
                                  if (j.eq.1)then
                                     write(400+n)"scalars r double 1",new_line(c)
                                  end if
                                  if (j.eq.2)then
                                     write(400+n)"scalars u double 1",new_line(c)
                                  end if
                                  if (j.eq.3)then
                                     write(400+n)"scalars v double 1",new_line(c)
                                  end if
                                  if (j.eq.4)then
                                     write(400+n)"scalars p double 1",new_line(c)
                                  end if
                                  if (j.eq.5)then
                                     write(400+n)"scalars species1 double 1",new_line(c)
                                  end if
                                  if (j.eq.6)then
                                     write(400+n)"scalars species2 double 1",new_line(c)
                                  end if
                                  if (j.eq.7)then
                                     write(400+n)"scalars volumef double 1",new_line(c)
                                  end if
                                  write(400+n)"lookup_table default",new_line(c)
                                  do i=1,imaxe
                                     xbin(xmpi_re(i))=xbin2(i)
                                  end do
                                  do i=1,imaxe
                                     write(400+n) xbin(i)  !test, not clear of whether cell or node based
                                  end do
                               end if
                            end do

                            
                            if (turbulence.eq.1)then

                            do i=1,kmaxe
                                
                                valuess(i)=u_ct_val(1,1,i)
                                end do
                            
                            call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

                            if (n.eq.0)then
                            
                                write(400+n)
                                write(400+n)"scalars  nut double 1"//lf
                                write(400+n)"lookup_table default"//lf
                            
                                do i=1,imaxe
                                xbin(xmpi_re(i))=xbin2(i)
                                end do
                            write(400+n)xbin(1:imaxe)
                            end if
                        end if
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            if (n.eq.0)then
                               ! close file
                               close(400+n) 

                               deallocate(xbin,xbin2)
                               deallocate(out1)
                            end if
                            deallocate (valuess)

                            call mpi_barrier(mpi_comm_world,ierror)
                            deallocate(variables)

    end subroutine outwritepara2db




subroutine outwritepara3dbp
!> @brief
!> this subroutine writes the solution and the grid file in binary vtk format
use iso_c_binding
implicit none
integer::kmaxe
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,outfile2,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
kmaxe=xmpielrank(n)

allocate(xbin(kmaxe))
write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') n 
	!proc4=".plt
	outfile="OUT_"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)


lf=char(10)


open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') kmaxn
write(400+n)"points "//str1//" double"//lf


    
    do i=1,kmaxn
	write(400+n)inoder4_cord(1,i),inoder4_cord(2,i),inoder4_cord(3,i)
	end do
   
                    
write(str1(1:15),'(i15)') kmaxe
write(str2(1:15),'(i15)') (kmaxe*8)+kmaxe
write(400+n)"cells",str1//str2//lf
do i=1,kmaxe
write(400+n)8,el_connect(i,1)-1,el_connect(i,2)-1,el_connect(i,3)-1,el_connect(i,4)-1,el_connect(i,5)-1,el_connect(i,6)-1,el_connect(i,7)-1,el_connect(i,8)-1
end do
                   
                    
write(str1(1:15),'(i15)') kmaxe                    
write(400+n)"cell_types"//str1//lf	
do i=1,kmaxe
write(400+n)12

end do

write(400+n)"cell_data"//str1//lf
!     
   
do j=1,nof_variables


     do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  xbin(i)=leftv(j)
		end do



    if (j.eq.1)then
    write(400+n)"scalars  r double 1"//lf
    end if
    if (j.eq.2)then
    write(400+n)
    write(400+n)"scalars  u double 1"//lf
    end if
    if (j.eq.3)then
    write(400+n)
    write(400+n)"scalars  v double 1"//lf
    end if
    if (j.eq.4)then
    write(400+n)
    write(400+n)"scalars  w double 1"//lf
    end if
    if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  p  double 1"//lf
    end if
    if (j.eq.6)then
    write(400+n)
    write(400+n)"scalars  species1  double 1"//lf
    end if
    if (j.eq.7)then
    write(400+n)
    write(400+n)"scalars  species2 double 1"//lf
    end if
    if (j.eq.8)then
    write(400+n)
    write(400+n)"scalars  volumef  double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf
    
    
		
    
    
     write(400+n)xbin(1:kmaxe)
     
    

     
    
    
end do

    if (itestcase.eq.4)then

	 do i=1,kmaxe
		  
		  xbin(i)=ielem_vortex(1,i)
		end do

        write(400+n)
    	write(400+n)"scalars  q double 1"//lf
	    write(400+n)"lookup_table default"//lf
    
    write(400+n)xbin(1:imaxe)
    end if
    
    if (turbulence.eq.1)then

	 do i=1,kmaxe
		  
		  xbin(i)=u_ct_val(1,1,i)
		end do

        write(400+n)
    	write(400+n)"scalars  nut double 1"//lf
	    write(400+n)"lookup_table default"//lf
    
    write(400+n)xbin(1:imaxe)
    end if
    
    
    
    close(400+n)
    deallocate(xbin,out1)
    
    
  
    
  
  

  call mpi_barrier(mpi_comm_world,ierror)

        


	
	

end subroutine outwritepara3dbp




subroutine outwritepara3dbpav
!> @brief
!> this subroutine writes the solution and the grid file in binary vtk format
use iso_c_binding
implicit none
integer::kmaxe
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2
 real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
kmaxe=xmpielrank(n)




allocate(xbin(kmaxe))
write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)')n 
	!proc4=".plt
	outfile="out_av"//trim(adjustl(proc3))//"_"//trim(adjustl(proc5))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)


lf=char(10)


open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') kmaxn
write(400+n)"points "//str1//" double"//lf


    
    do i=1,kmaxn
	write(400+n)inoder4_cord(1,i),inoder4_cord(2,i),inoder4_cord(3,i)
	end do
   
                    
write(str1(1:15),'(i15)') kmaxe
write(str2(1:15),'(i15)') (kmaxe*8)+kmaxe
write(400+n)"cells",str1//str2//lf
do i=1,kmaxe
write(400+n)8,el_connect(i,1)-1,el_connect(i,2)-1,el_connect(i,3)-1,el_connect(i,4)-1,el_connect(i,5)-1,el_connect(i,6)-1,el_connect(i,7)-1,el_connect(i,8)-1
end do
                   
                    
write(str1(1:15),'(i15)') kmaxe                    
write(400+n)"cell_types"//str1//lf	
do i=1,kmaxe
write(400+n)12
end do

write(400+n)"cell_data"//str1//lf
!     
   
do j=1,11

    if (j.eq.1)then
    write(400+n)"scalars  r_mean double 1"//lf
    end if
    if (j.eq.2)then
    write(400+n)
    write(400+n)"scalars  u_mean double 1"//lf
    end if
    if (j.eq.3)then
    write(400+n)
    write(400+n)"scalars  v_mean double 1"//lf
    end if
    if (j.eq.4)then
    write(400+n)
    write(400+n)"scalars  w_mean double 1"//lf
    end if
    if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  p_mean  double 1"//lf
    end if
    if (j.eq.6)then
    write(400+n)
    write(400+n)"scalars  u_rms  double 1"//lf
    end if
    if (j.eq.7)then
    write(400+n)
    write(400+n)"scalars  v_rms double 1"//lf
    end if
    if (j.eq.8)then
    write(400+n)
    write(400+n)"scalars  w_rms  double 1"//lf
    end if
     if (j.eq.9)then
    write(400+n)
    write(400+n)"scalars  uv  double 1"//lf
    end if
     if (j.eq.10)then
    write(400+n)
    write(400+n)"scalars  uw  double 1"//lf
    end if
     if (j.eq.11)then
    write(400+n)
    write(400+n)"scalars  wv  double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf
    
        
        if (j.eq.1)then
		do i=1,kmaxe
		xbin(i)=u_c_val(ind1,j,i)
		end do
		end if
		if ((j.gt.1).and.(j.lt.5))then
		do i=1,kmaxe
		xbin(i)=u_c_val(ind1,j,i)/u_c_val(ind1,1,i)
		end do
		end if
		if (j.eq.5)then
		do i=1,kmaxe
		leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
		call cons2prim(n,leftv,mp_pinfl,gammal)
		xbin(i)=leftv(5)
		end do
		end if
		if (j.ge.6)then
		do i=1,kmaxe
		xbin(i)=u_c_val(1,j,i)
		end do
		end if
    
    
     write(400+n)xbin(1:kmaxe)
   

     
    
    
end do
    
    
    
    close(400+n)
        deallocate(xbin,out1)
     
    
  
    
  
  

  call mpi_barrier(mpi_comm_world,ierror)

        


	
	

end subroutine outwritepara3dbpav





subroutine outwritepara3dsb
!> @brief
!> this subroutine writes the solution and the surface file in binary vtk format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8,icount_wall
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(12))
nvar1=2





if (n.eq.0)then

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') 
	!proc4=".plt"
	outfile="SURF_"//trim(adjustl(proc3))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
end if

lf=char(10)



if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)then
open(96,file='GRID.bnd',form='formatted',status='old',action='read')
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
 open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
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
! open(400+n,file=outfile,form='unformatted',status='new',action='write')
 open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') igf2
write(400+n)"points "//str1//" double"//lf

allocate (valuelocation(nvar1))

    if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+n)x/scaler,y/scaler,z/scaler
	end if
	end do

    close(96)
    else
    open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+n)x/scaler,y/scaler,z/scaler
	end if
	end do

    close(96)
    
    
    end if

    
   
    
    
     if (binio.eq.0)then
     open(98,file='GRID.bnd',form='formatted',status='old',action='read')
     end if
      if (binio.eq.1)then
     open(98,file='GRID.bnd',form='unformatted',status='old',action='read')
     end if
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		if (binio.eq.0)then
		do k=1,imaxb
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		else
		do k=1,imaxb
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		
		end if
	    
 		close(98)
 		
 		
 		
 		
 		write(str1(1:15),'(i15)') igf2
 write(str2(1:15),'(i15)') (igf2*4)+igf2
 write(400+n)"cells",str1//str2//lf
 		do i=1,igf2
 		write(400+n)4,icon(1,i)-1,icon(2,i)-1,icon(3,i)-1,icon(4,i)-1
 		end do
 		
	
		deallocate(icon)	
		deallocate(inog)


                    
! write(str1(1:15),'(i15)') imaxe
! write(str2(1:15),'(i15)') (imaxe*8)+imaxe
! write(400+n)"cells",str1//str2//lf
! 
! open(97,file='GRID.cel',form='formatted',status='old',action='read')
! do i=1,imaxe
! read(97,*)j,j1,j2,j3,j4,j5,j6,j7,j8
! write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
! end do
! close(97)
                    
                    
                 
 write(400+n)"cell_types"//str1//lf	
 do i=1,igf2 
 write(400+n)9
 end do
! 
write(400+n)"cell_data"//str1//lf
!                
!                     
!                     
! 
allocate(xbin(totwalls),xbin2(totwalls))

 end if

  totiw=xmpiwall(n)
  if (xmpiwall(n).gt.0)then
  allocate(valuess(xmpiwall(n)))
  end if
!     
   
do j=1,nof_variables
     
     if ((j.ge.2).and.(j.le.4))then
     
         if (totiw.gt.0)then
		      	
		      do i=1,totiw
                valuess(i)=u_c_rms(j-5,ibound_t(i))/u_c_val(1,1,ibound_t(i))
                end do
				
		      end if
    end if
    if (j.eq.1)then
         if (totiw.gt.0)then
		      	
		      do i=1,totiw
                valuess(i)=u_c_val(1,j,ibound_t(i))
                end do
				
		      end if
    end if
    if (j.eq.5)then
                if (totiw.gt.0)then
		     do i=1,totiw
						leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ibound_t(i))
						call cons2prim(n,leftv,mp_pinfl,gammal)
						valuess(i)=leftv(5)
						
					  end do
		      end if
    
    end if
    
		
    
    call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)
!     call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    
    
    
    
    
    if (j.eq.1)then
    write(400+n)"scalars  r double 1"//lf
    end if
    if (j.eq.2)then
    write(400+n)
    write(400+n)"scalars  u double 1"//lf
    end if
    if (j.eq.3)then
    write(400+n)
    write(400+n)"scalars  v double 1"//lf
    end if
    if (j.eq.4)then
    write(400+n)
    write(400+n)"scalars  w double 1"//lf
    end if
    if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  p  double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf
    
    
    
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
		
    
    
     write(400+n)xbin(1:totwalls)
    end if

     
    
    
end do
    
    
    
    

     
    
 
    
  if (n.eq.0)then
  
   
   deallocate(xbin,xbin2,valuelocation)
  deallocate(out1)
  close(400+n)  
  end if
   if (totiw.gt.0)then
  deallocate (valuess)
  end if
  

          
        
 deallocate(variables)


	
	

end subroutine outwritepara3dsb






subroutine outwritepara3dbav
!> @brief
!> this subroutine writes the averaged solution and the grid file in binary vtk format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(12))
nvar1=2
kmaxe=xmpielrank(n)




if (n.eq.0)then

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') 
	!proc4=".plt"
	outfile="out_av"//trim(adjustl(proc3))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
end if

lf=char(10)

if (n.eq.0)then

open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') imaxn
write(400+n)"points "//str1//" double"//lf

allocate (valuelocation(nvar1))


valuelocation(:)=0
valuelocation(1:2)=1

    
     
! 	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
!         do i=1,imaxn
! 	read(96,*)j,x,y,z
! 	xbin(i)=x/scaler
! 	ybin(i)=y/scaler
!  	zbin(i)=z/scaler
! 	end do
! 
!     close(96)
    
    if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
        do i=1,imaxn
	read(96,*)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n)x,y,z
	end do
	close(96)
	else
	open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
        do i=1,imaxn
	read(96)j,x,y,z
	x=x/scaler;y=y/scaler;z=z/scaler
	write(400+n)x,y,z
	end do
	close(96)
	end if

    
    
   
    
                    
write(str1(1:15),'(i15)') imaxe
write(str2(1:15),'(i15)') (imaxe*8)+imaxe
write(400+n)"cells",str1//str2//lf
if (binio.eq.0)then
open(97,file='GRID.cel',form='formatted',status='old',action='read')
do i=1,imaxe
read(97,*)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)
else
open(97,file='GRID.cel',form='unformatted',status='old',action='read')
do i=1,imaxe
read(97)j,j1,j2,j3,j4,j5,j6,j7,j8
write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
end do
close(97)

end if
                    
                    
write(str1(1:15),'(i15)') imaxe                    
write(400+n)"cell_types"//str1//lf	
do i=1,imaxe
write(400+n)12
end do

write(400+n)"cell_data"//str1//lf
!                
!                     
!                     
! 
 
  allocate(xbin(imaxe),xbin2(imaxe))
	
! 
 end if
! 
!  
  
  allocate(valuess(kmaxe))
 
 call mpi_barrier(mpi_comm_world,ierror)
!     
   
do j=1,nof_variables+6
     
     if ((j.ge.2).and.(j.le.4))then
     
        do i=1,kmaxe
        valuess(i)=u_c_val(ind1,j,i)/u_c_val(ind1,1,i)!0.0
        end do
    end if
    if (j.eq.1)then
        do i=1,kmaxe
        valuess(i)=u_c_val(ind1,j,i)
        end do
    end if
    if (j.eq.5)then
                do i=1,kmaxe
		  leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
		  call cons2prim(n,leftv,mp_pinfl,gammal)
		  valuess(i)=leftv(5)
		end do
    
    end if
    if (j.gt.5)then
    do i=1,kmaxe
	valuess(i)=u_c_val(1,j,i)
    end do
    end if
    
		
    
    
    call mpi_gatherv(valuess,xmpiall(n),mpi_double_precision,xbin2,xmpiall,offset,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    
    if (j.eq.1)then
    write(400+n)"scalars  r_mean double 1"//lf
    end if
    if (j.eq.2)then
    write(400+n)
    write(400+n)"scalars  u_mean double 1"//lf
    end if
    if (j.eq.3)then
    write(400+n)
    write(400+n)"scalars  v_mean double 1"//lf
    end if
    if (j.eq.4)then
    write(400+n)
    write(400+n)"scalars  w_mean double 1"//lf
    end if
    if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  p_mean  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  u_rms  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  v_rms  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  w_rms  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  uv  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  uw  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  wv  double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf
    
    
		do i=1,imaxe
		xbin(xmpi_re(i))=xbin2(i)
		end do
    
     write(400+n)xbin(1:imaxe)
    end if
    
     
    
    
end do
    
    
    
    

     
    
  
    
  if (n.eq.0)then
  close(400+n) 
   
   deallocate(xbin,xbin2,valuelocation)
  deallocate(out1)
  end if
  deallocate (valuess)
  

  call mpi_barrier(mpi_comm_world,ierror)

        
 deallocate(variables)


	
	

end subroutine outwritepara3dbav

subroutine outwritepara3dsbav
!> @brief
!> this subroutine writes the averaged solution and the surface file in ascii vtk format
use iso_c_binding
implicit none

! external tecini112
! external teczne112
! external tecdat112
! external tecnode112
! external  tecend112

! 
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
integer::kmaxe,kk,kfk,icpuid,l,ihgt,ihgj,kkd
real::x,y,z,denominator,tuy,tvx,twx,tuz,tvz,twy,snorm,onorm
real,allocatable,dimension(:)::ifint,tfint,ndr,nds
integer::ineedt,jj,ix,ix1,i1,i2,i3,i4,i5,decomf,kd
real,allocatable,dimension(:)::variables
real,dimension(3,3)::avort,tvort,svort,ovort
integer::inx,i,k,j,m,o,p,q,jk,imax,jmax,kmax,igf,igf2,dumg,duml,imaxp,nvar1,j1,j2,j3,j4,j5,j6,j7,j8,icount_wall
logical::herev
real,dimension(5)::total
 character(len=30)::proc,outfile,proc3,surfile,proc4,proc5
integer::ierr,cv,tecini112,teczne112,tecdat112,tecnode112,tecend112,itgfd
real,allocatable,dimension(:)::xbin,ybin,zbin,xbin2
real,allocatable,dimension(:,:)::fbin
integer,allocatable,dimension(:,:)::icon
integer,allocatable,dimension(:)::valuelocation,inog,icell,icella
real,allocatable,dimension(:)::valuess,valuesa
character(len=:),allocatable::out1
character*1 nulchar
character(len=1)   :: flui,lf
character(len=15)  :: str1,str2
 
      integer::   debug,iii,npts,nelm

  
      real::    soltime
      integer:: visdouble, filetype
      integer:: zonetype,strandid,parentzn,isblock
      integer:: icellmax,jcellmax,kcellmax,nfconns,fnmode,shrconn
      pointer   (nullptr,null)
      integer:: null(*)
allocate(variables(12))
nvar1=2





if (n.eq.0)then

write(proc3,fmt='(i10)') it
write(proc5,fmt='(i10)') 
	!proc4=".plt"
	outfile="SURF_"//trim(adjustl(proc3))//".vtk"!//trim(adjustl(proc4))
	itgfd=len_trim(outfile)
	allocate(character(len=itgfd) ::out1)
	out1=outfile(1:itgfd)
end if

lf=char(10)



if ((n.eq.0).and.(totwalls.gt.0))then
if (binio.eq.0)then
open(96,file='GRID.bnd',form='formatted',status='old',action='read')
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
 open(96,file='GRID.bnd',form='unformatted',status='old',action='read')
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
! open(400+n,file=outfile,form='unformatted',status='new',action='write')
 open(400+n,file=outfile,status='replace',access='stream',convert='big_endian')
write(400+n)"# vtk datafile version 3.0"//lf
write(400+n)"vtk output"//lf
write(400+n)"binary"//lf
write(400+n)"dataset unstructured_grid"//lf
write(400+n)"field fielddata 1"//lf
write(400+n)"time 1 1 double"//lf
write(400+n) t
write(str1(1:15),'(i15)') igf2
write(400+n)"points "//str1//" double"//lf

allocate (valuelocation(nvar1))

    if (binio.eq.0)then
	open(96,file='GRID.vrt',form='formatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96,*)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+n)x/scaler,y/scaler,z/scaler
	end if
	end do

    close(96)
    else
    open(96,file='GRID.vrt',form='unformatted',status='old',action='read')
	igf2=0
        do i=1,imaxn
	read(96)j,x,y,z
	    
      if (inog(i).eq.1)then
      igf2=igf2+1
      inog(i)=igf2
      write(400+n)x/scaler,y/scaler,z/scaler
	end if
	end do

    close(96)
    
    
    end if

    
   
    
    
     if (binio.eq.0)then
     open(98,file='GRID.bnd',form='formatted',status='old',action='read')
     end if
      if (binio.eq.1)then
     open(98,file='GRID.bnd',form='unformatted',status='old',action='read')
     end if
	  allocate(icon(4,totwalls))
    icon=0
     cv=0
		igf2=0
		if (binio.eq.0)then
		do k=1,imaxb
               
 		read(98,*)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		else
		do k=1,imaxb
               
 		read(98)igf,i,j,l,m,o
		  if (o.eq.4)then
		  igf2=igf2+1
		  icon(1,igf2)=inog(i)
		  icon(2,igf2)=inog(j)
		  icon(3,igf2)=inog(l)
		  icon(4,igf2)=inog(m)
    
		end if
    !cv=cv+4
        	
		end do
		
		end if
	    
 		close(98)
 		
 		
 		
 		
 		write(str1(1:15),'(i15)') igf2
 write(str2(1:15),'(i15)') (igf2*4)+igf2
 write(400+n)"cells",str1//str2//lf
 		do i=1,igf2
 		write(400+n)4,icon(1,i)-1,icon(2,i)-1,icon(3,i)-1,icon(4,i)-1
 		end do
 		
	
		deallocate(icon)	
		deallocate(inog)


                    
! write(str1(1:15),'(i15)') imaxe
! write(str2(1:15),'(i15)') (imaxe*8)+imaxe
! write(400+n)"cells",str1//str2//lf
! 
! open(97,file='GRID.cel',form='formatted',status='old',action='read')
! do i=1,imaxe
! read(97,*)j,j1,j2,j3,j4,j5,j6,j7,j8
! write(400+n)8,j1-1,j2-1,j3-1,j4-1,j5-1,j6-1,j7-1,j8-1
! end do
! close(97)
                    
                    
                 
 write(400+n)"cell_types"//str1//lf	
 do i=1,igf2 
 write(400+n)9
 end do
! 
write(400+n)"cell_data"//str1//lf
!                
!                     
!                     
! 
allocate(xbin(totwalls),xbin2(totwalls))

 end if

  totiw=xmpiwall(n)
  if (xmpiwall(n).gt.0)then
  allocate(valuess(xmpiwall(n)))
  end if
!     
   
do j=1,nof_variables+6
     
     if ((j.ge.2).and.(j.le.4))then
     
         if (totiw.gt.0)then
		      	
		      do i=1,totiw
                valuess(i)=u_c_rms(j-5,ibound_t(i))/u_c_val(ind1,1,ibound_t(i))
                end do
				
		      end if
    end if
    if (j.eq.1)then
         if (totiw.gt.0)then
		      	
		      do i=1,totiw
                valuess(i)=u_c_val(ind1,j,ibound_t(i))
                end do
				
		      end if
    end if
    if (j.eq.5)then
                if (totiw.gt.0)then
		     do i=1,totiw
						leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,ibound_t(i))
						call cons2prim(n,leftv,mp_pinfl,gammal)
						valuess(i)=leftv(5)
						
					  end do
		      end if
    
    end if
     if (j.gt.5)then
      if (totiw.gt.0)then
	do i=1,totiw
	  valuess(i)=u_c_val(1,j,ibound_t(i))
	end do
      end if
    end if
		
    
    call mpi_gatherv(valuess,xmpiwall(n),mpi_double_precision,xbin2,xmpiwall,woffset,mpi_double_precision,0,mpi_comm_world,ierror)
!     call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)

    if (n.eq.0)then
    
    if (j.eq.1)then
    write(400+n)"scalars  r_mean double 1"//lf
    end if
    if (j.eq.2)then
    write(400+n)
    write(400+n)"scalars  u_mean double 1"//lf
    end if
    if (j.eq.3)then
    write(400+n)
    write(400+n)"scalars  v_mean double 1"//lf
    end if
    if (j.eq.4)then
    write(400+n)
    write(400+n)"scalars  w_mean double 1"//lf
    end if
    if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  p_mean  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  u_rms  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  v_rms  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  w_rms  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  uv  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  uw  double 1"//lf
    end if
     if (j.eq.5)then
    write(400+n)
    write(400+n)"scalars  wv  double 1"//lf
    end if
    write(400+n)"lookup_table default"//lf
    
    
		do i=1,totwalls
		xbin(xmpi_wre(i))=xbin2(i)
		end do
    
     write(400+n)xbin(1:totwalls)
    end if

     
    
    
end do
    
    
    
    

     
    
 
    
  if (n.eq.0)then
  
   
   deallocate(xbin,valuelocation,xbin2)
  deallocate(out1)
  close(400+n)  
  end if
   if (totiw.gt.0)then
  deallocate (valuess)
  end if
  

          
        
 deallocate(variables)


	
	

end subroutine outwritepara3dsbav



subroutine fix_nodes_local
implicit none
integer::ndlc_count1,ndlc_count2,indexgt
integer::kmaxe,nd_lc_nodes,i,j,k,l,m,countfnodes,indfc
integer,allocatable,dimension(:)::ndlc_array1
integer,allocatable,dimension(:)::list_nin,list_nout

kmaxe=xmpielrank(n)
allocate(ndlc_array1(kmaxe*8));ndlc_array1(:)=0

ndlc_count1=0

allocate(list_nin(1:imaxn),list_nout(1:imaxn))
list_nin(:)=-10
list_nout(:)=-10

!find only the uniques first
countfnodes=0
do i=1,imaxn
    if (dinoder(i)%itor.gt.0)then
    countfnodes=countfnodes+1
    list_nin(i)=n
    end if
end do



call mpi_allreduce(list_nin,list_nout,imaxn,mpi_integer,mpi_max,mpi_comm_world,ierror)

!find only the uniques first
countfnodes=0
do i=1,imaxn
    if (list_nout(i).eq.n)then
    countfnodes=countfnodes+1
    end if
end do




deallocate(list_nin)



!find only the uniques first





do i=1,imaxn
    if (dinoder(i)%itor.gt.0)then

    ndlc_count1=ndlc_count1+1
    ndlc_array1(ndlc_count1)=i
    end if

end do

kmaxn=ndlc_count1




    allocate(inoder4_cord(1:dims,1:kmaxn))
    allocate(inoder4_bct(1:3,1:kmaxn));inoder4_bct(:,:)=0
    allocate(inoder4_itor(1:kmaxn))


allocate(my_nodesl(1:countfnodes),my_nodesg(countfnodes))
countfnodes=0
do i=1,kmaxn

    inoder4_itor(i)=ndlc_array1(i)
    inoder4_cord(1:dims,i)=dinoder(inoder4_itor(i))%cord(1:dims)
    dinoder(inoder4_itor(i))%itor=i
    if (list_nout(inoder4_itor(i)).eq.n)then
    countfnodes=countfnodes+1
!     inoder4_itorm(i)=countfnodes
    my_nodesl(countfnodes)=i
    my_nodesg(countfnodes)=inoder4_itor(i)
    end if
end do






allocate(xmpiall_v(0:isize-1),offset_v(0:isize-1))

xmpiall_v=0

xmpiall_v(n)=countfnodes



call mpi_barrier(mpi_comm_world,ierror)


call mpi_allgather(countfnodes,1,mpi_integer,xmpiall_v,1,mpi_integer,mpi_comm_world,ierror)




offset_v(0)=0
do i=1,isize-1
offset_v(i)=offset_v(i-1)+xmpiall_v(i-1)
end do










!now we fixed the correct numbering for writing vtk output



do i=1,kmaxe
	if (tecplot.eq.5)then
	if (dimensiona.eq.2)then
	if (ielem_ishape(i).eq.5)then !quad
	indfc=4
	else
	indfc=3
	end if

	else
	if (ielem_ishape(i).eq.1)then !hexa
	indfc=8
	end if
	if (ielem_ishape(i).eq.2)then !tetra
	indfc=4
	end if
	if (ielem_ishape(i).eq.3)then !pyramid
	indfc=5
	end if
	if (ielem_ishape(i).eq.4)then !prism
	indfc=6
	end if



	end if

	end if



    do j=1,ielem_nonodes(i)
    ielem_nodes(j,i)=dinoder(ielem_nodes(j,i))%itor
    end do





    if (tecplot.eq.5)then
					if (dimensiona.eq.2)then
						if (ielem_ishape(i).eq.5)then !quad
						ielem_nodes_v(1,i)=inoder4_itor(ielem_nodes(1,i))-1
						ielem_nodes_v(2,i)=inoder4_itor(ielem_nodes(2,i))-1
						ielem_nodes_v(3,i)=inoder4_itor(ielem_nodes(3,i))-1
						ielem_nodes_v(4,i)=inoder4_itor(ielem_nodes(4,i))-1
						else							!triangle
						ielem_nodes_v(1,i)=inoder4_itor(ielem_nodes(1,i))-1
						ielem_nodes_v(2,i)=inoder4_itor(ielem_nodes(2,i))-1
						ielem_nodes_v(3,i)=inoder4_itor(ielem_nodes(3,i))-1
						!ielem_nodes_v(4,i)=inoder4_itor(ielem_nodes(3,i))-1

						end if
					end if
					if (dimensiona.eq.3)then
						if (ielem_ishape(i).eq.1)then !hexa
						ielem_nodes_v(1,i)=inoder4_itor(ielem_nodes(1,i))-1
						ielem_nodes_v(2,i)=inoder4_itor(ielem_nodes(2,i))-1
						ielem_nodes_v(3,i)=inoder4_itor(ielem_nodes(3,i))-1
						ielem_nodes_v(4,i)=inoder4_itor(ielem_nodes(4,i))-1
						ielem_nodes_v(5,i)=inoder4_itor(ielem_nodes(5,i))-1
						ielem_nodes_v(6,i)=inoder4_itor(ielem_nodes(6,i))-1
						ielem_nodes_v(7,i)=inoder4_itor(ielem_nodes(7,i))-1
						ielem_nodes_v(8,i)=inoder4_itor(ielem_nodes(8,i))-1
						end if
						if (ielem_ishape(i).eq.2)then !tetra
						ielem_nodes_v(1,i)=inoder4_itor(ielem_nodes(1,i))-1
						ielem_nodes_v(2,i)=inoder4_itor(ielem_nodes(2,i))-1
						ielem_nodes_v(3,i)=inoder4_itor(ielem_nodes(3,i))-1
						ielem_nodes_v(4,i)=inoder4_itor(ielem_nodes(4,i))-1
						end if
						if (ielem_ishape(i).eq.3)then !pyramid
						ielem_nodes_v(1,i)=inoder4_itor(ielem_nodes(1,i))-1
						ielem_nodes_v(2,i)=inoder4_itor(ielem_nodes(2,i))-1
						ielem_nodes_v(3,i)=inoder4_itor(ielem_nodes(3,i))-1
						ielem_nodes_v(4,i)=inoder4_itor(ielem_nodes(4,i))-1
						ielem_nodes_v(5,i)=inoder4_itor(ielem_nodes(5,i))-1
						end if
						if (ielem_ishape(i).eq.4)then !prism

						ielem_nodes_v(1,i)=inoder4_itor(ielem_nodes(1,i))-1
						ielem_nodes_v(2,i)=inoder4_itor(ielem_nodes(2,i))-1
						ielem_nodes_v(3,i)=inoder4_itor(ielem_nodes(3,i))-1
						ielem_nodes_v(4,i)=inoder4_itor(ielem_nodes(4,i))-1
						ielem_nodes_v(5,i)=inoder4_itor(ielem_nodes(5,i))-1
						ielem_nodes_v(6,i)=inoder4_itor(ielem_nodes(6,i))-1

						end if

					end if




			do l=1,ielem_ifca(i)
					if (dimensiona.eq.2)then
					do j=1,2
							indexgt=dinoder(ielem_nodes_faces(l,j,i))%itor
							ielem_nodes_faces_v(l,j,i)=inoder4_itor(indexgt)-1
					end do

					else
						if (ielem_types_faces(l,i).eq.5)then
							do j=1,4
								indexgt=dinoder(ielem_nodes_faces(l,j,i))%itor
								ielem_nodes_faces_v(l,j,i)=inoder4_itor(indexgt)-1
							end do
						end if
						if (ielem_types_faces(l,i).eq.6)then
							do j=1,3
								indexgt=dinoder(ielem_nodes_faces(l,j,i))%itor
								ielem_nodes_faces_v(l,j,i)=inoder4_itor(indexgt)-1

							end do
						end if
					end if
			end do

	end if



end do



	do i=1,kmaxe


    do l=1,ielem_ifca(i)
         if (dimensiona.eq.2)then
          do j=1,2
                ielem_nodes_faces(l,j,i)=dinoder(ielem_nodes_faces(l,j,i))%itor

            end do

            else
        if (ielem_types_faces(l,i).eq.5)then
            do j=1,4
                ielem_nodes_faces(l,j,i)=dinoder(ielem_nodes_faces(l,j,i))%itor
            end do
        end if
        if (ielem_types_faces(l,i).eq.6)then
            do j=1,3
                ielem_nodes_faces(l,j,i)=dinoder(ielem_nodes_faces(l,j,i))%itor
            end do
        end if



        end if
    end do
end do














deallocate (dinoder,ndlc_array1,list_nout)


if (tecplot.eq.3)then

if (dimensiona.eq.3)then

allocate(el_connect(1:kmaxe,1:8))
do i=1,kmaxe
    if (ielem_ishape(i).eq.1)then!hexa
    el_connect(i,1:8)=ielem_nodes(1:8,i)
    end if
     if (ielem_ishape(i).eq.2)then!tetra
    el_connect(i,1:2)=ielem_nodes(1:2,i)
    el_connect(i,3)=ielem_nodes(3,i)
    el_connect(i,4)=ielem_nodes(3,i)
    el_connect(i,5)=ielem_nodes(4,i)
    el_connect(i,6)=ielem_nodes(4,i)
    el_connect(i,7)=ielem_nodes(4,i)
    el_connect(i,8)=ielem_nodes(4,i)
    end if
    if (ielem_ishape(i).eq.3)then!pyramid
    el_connect(i,1:5)=ielem_nodes(1:5,i)
    el_connect(i,6)=ielem_nodes(5,i)
    el_connect(i,7)=ielem_nodes(5,i)
    el_connect(i,8)=ielem_nodes(5,i)
    end if
    if (ielem_ishape(i).eq.4)then!prism
    el_connect(i,1:3)=ielem_nodes(1:3,i)
    el_connect(i,4)=ielem_nodes(3,i)
    el_connect(i,5:7)=ielem_nodes(4:6,i)
    el_connect(i,8)=ielem_nodes(6,i)
    end if

end do

end if

end if
!array el_connect(index1,indices 1:8)
!first index is the element index
!second index varies from 1:8 for all the vertices of this element (connectivity list)
!for paraview the numbering might require switching (from 1 to number of nodes----> 0  to number of nodes-1)

!inoder4_cord(1:3,1:number of nodes (kmaxn)) holds the coordinates for each point 
!kmaxn is global across all modules
!kmaxe must be set as equal to xmpielrank(n)
!this is shown in outwritepara3dbp (p stands for partitioned mesh writing)















end subroutine fix_nodes_local




subroutine checkpointv4(n)
!> @brief
!> this subroutine is writing the checkpointing files
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella
real,allocatable,dimension(:)::valuesa,valuess
real,allocatable,dimension(:,:)::xbin
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj,igfs
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
character(len=20)::proc,restfile,proc3
real,allocatable,dimension(:)::igint,tgint
 kmaxe=xmpielrank(n)
igfs=t



#ifdef gpu
                  !$omp target update from(u_c_val(1,1:nof_variables,1:xmpielrank(n)))
#endif

write(proc3,fmt='(i10)') igfs
	restfile="REST_"//trim(adjustl(proc3))//".dat"!//trim(adjustl(proc4))
	
icpuid=n

	if (n.eq.0)then
	open(1086,file=restfile,form='unformatted',status='replace',action='write')
	end if
	
	kmaxe=xmpielrank(n)
    
dumg=kmaxe


call mpi_barrier(mpi_comm_world,ierror) !not needed

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


    if (n.eq.0)then
	    allocate(icella(imaxp*isize))
	    icella=0

    end if

    call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

    deallocate (icell)

    if (n.eq.0) then
    allocate(valuesa(imaxp*isize))
    allocate(xbin(imaxe,5+turbulenceequations+passivescalar))
    valuesa=zero
    end if
    
    
    allocate(valuess(imaxp));valuess=zero

  
  do jj=1,5
	do i=1,kmaxe
        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
        call cons2prim(n,leftv,mp_pinfl,gammal)
		valuess(i)=leftv(jj)
	end do
	call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
	
	    if (n.eq.0)then
	    do i=1,imaxp*isize
		if (icella(i).gt.0)then
		xbin(icella(i),jj)=valuesa(i)
		end if
	    end do
	    end if
    
  end do



!   call mpi_barrier(mpi_comm_world,ierror)

    if (n.eq.0)then

    do i=1,imaxe
    write(1086)xbin(i,1:nof_variables)
    end do

    deallocate(xbin,icella,valuesa)
    close(1086)
    end if
    
deallocate(valuess)







end subroutine checkpointv4


subroutine checkpointv3(n)
!> @brief
!> this subroutine is writing the checkpointing files
implicit none
integer,intent(in)::n
integer,allocatable,dimension(:)::icell,icella
real,allocatable,dimension(:)::valuesa,valuess
real,allocatable,dimension(:,:)::xbin
integer::i,k,kmaxe,j,jk,icpuid,nvar,imaxp,dumg,duml,jj
character(len=20)::proc,restfile,proc3
real,allocatable,dimension(:)::igint,tgint
 kmaxe=xmpielrank(n)



	restfile="CORD.dat"
	
icpuid=n

	if (n.eq.0)then
	open(1086,file=restfile,form='unformatted',status='replace',action='write')
	end if
	
	kmaxe=xmpielrank(n)
    
dumg=kmaxe


call mpi_barrier(mpi_comm_world,ierror) !not needed

call mpi_allreduce(dumg,duml,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
imaxp=duml

allocate(icell(imaxp))
icell=0

do i=1,kmaxe
  icell(i)=ielem_ihexgl(i)
end do


    if (n.eq.0)then
	    allocate(icella(imaxp*isize))
	    icella=0

     else
		allocate(icella(1))

     end if

    call mpi_gather(icell,imaxp,mpi_integer,icella,imaxp,mpi_integer,0,mpi_comm_world,ierror)

    deallocate (icell)

	if (n.eq.0) then
    allocate(valuesa(imaxp*isize))

    allocate(xbin(imaxe,3))
    valuesa=zero
    else
    allocate(valuesa(1),xbin(1,3))

    end if
    
    
    allocate(valuess(imaxp));valuess=zero

  
  do jj=1,3
	do i=1,kmaxe
        if (jj.eq.1)then
		valuess(i)=ielem_xxc(i)
		end if
		if (jj.eq.2)then
		valuess(i)=ielem_yyc(i)
		end if
		if (jj.eq.3)then
		valuess(i)=ielem_zzc(i)
		end if
	end do
	call mpi_gather(valuess,imaxp,mpi_double_precision,valuesa,imaxp,mpi_double_precision,0,mpi_comm_world,ierror)
	
	    if (n.eq.0)then
	    do i=1,imaxp*isize
		if (icella(i).gt.0)then
		xbin(icella(i),jj)=valuesa(i)
		end if
	    end do
	    end if
    
  end do



   call mpi_barrier(mpi_comm_world,ierror)

    if (n.eq.0)then

    do i=1,imaxe
    write(1086)xbin(i,1:3)
    end do


    close(1086)
    end if

    deallocate(xbin,icella,valuesa)


    
deallocate(valuess)

call mpi_barrier(mpi_comm_world,ierror) !not needed





end subroutine checkpointv3


subroutine troubled_history
implicit none
integer::i,j,k,traj1,traj2,traj3,traj4,kmaxe,writeid,writeconf
real::win1,win2,win3,win4,post,post1,post2,post3,post4



kmaxe=xmpielrank(n)

!$omp barrier
!$omp master
post1=0
traj1=0
pos_l(:)=zero
pos_g(:)=zero
ipos_l(:)=0
ipos_g(:)=0
!$omp end master
!$omp barrier

if (mood.gt.0)then

#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:ipos_l(1:2)) &
!$omp& private(i) reduction(+:ipos_l(1), ipos_l(2))
#else
!$omp do reduction(+:ipos_l(1), ipos_l(2))
#endif
do i=1,kmaxe
    if (ielem_mood_o(i).lt.(iorder+1))then
        ipos_l(1)=ipos_l(1)+1
    end if
    if (ielem_mood_o(i).eq.1)then
        ipos_l(2)=ipos_l(2)+1
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



!$omp barrier
!$omp master
call mpi_allreduce(ipos_l(1:2),ipos_g(1:2),2,mpi_integer,mpi_sum,mpi_comm_world,ierror)

pos_g(1)=ipos_g(1)
pos_g(2)=ipos_g(2)
pos_g(1)=(pos_g(1)/imaxe)*100
pos_g(2)=(pos_g(2)/imaxe)*100

if (n.eq.0)then
    open(70,file='troubled.dat',form='formatted',action='write',position='append')
    write(70,'(e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),pos_g(2)
    close(70)
end if

call mpi_barrier(mpi_comm_world,ierror)
!$omp end master
!$omp barrier

else

#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:pos_l(1)) &
!$omp& private(i) reduction(+:pos_l(1))
#else
!$omp do reduction(+:pos_l(1))
#endif
do i=1,kmaxe
    pos_l(1)=pos_l(1)+ielem_condition(i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


!$omp barrier
!$omp master

call mpi_allreduce(pos_l(1),pos_g(1),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

if (n.eq.0)then
    open(70,file='troubled.dat',form='formatted',action='write',position='append')
    write(70,'(e14.7,1x,e14.7,1x,e14.7)')t,(pos_g(1)/imaxe)*100.0
    close(70)
end if

call mpi_barrier(mpi_comm_world,ierror)
!$omp end master
!$omp barrier

end if

end subroutine troubled_history


subroutine reduced_history
implicit none
integer::i,j,k,traj1,traj2,traj3,traj4,kmaxe,writeid,writeconf,tempint
real::win1,win2,win3,win4,post,post1,post2,post3,post4


kmaxe=xmpielrank(n)
!$omp barrier
!$omp master
post1=0
traj1=0
pos_l(1)=zero
pos_g(1)=zero
!$omp end master
!$omp barrier

#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:pos_l(1)) &
!$omp& private(i,tempint) reduction(+:pos_l(1))
#else
!$omp do private(tempint) reduction(+:pos_l(1))
#endif
do i=1,kmaxe
    tempint=0
    if (ielem_reduce(i).ge.1)then
        tempint=1
    end if
    pos_l(1)=pos_l(1)+tempint
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


!$omp barrier
!$omp master

call mpi_allreduce(pos_l(1),pos_g(1),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

if (n.eq.0)then
    open(70,file='reduced.dat',form='formatted',action='write',position='append')
    write(70,'(e14.7,1x,e14.7,1x,e14.7)')t,(pos_g(1)/imaxe)*100.0
    close(70)
end if

call mpi_barrier(mpi_comm_world,ierror)


!$omp end master
!$omp barrier

end subroutine reduced_history


subroutine filtered_history
implicit none
integer::i,j,k,traj1,traj2,traj3,traj4,kmaxe,writeid,writeconf,countfd
real::win1,win2,win3,win4,post,post1,post2,post3,post4


kmaxe=xmpielrank(n)

!$omp barrier
!$omp master
post1=0
traj1=0
pos_l(:)=zero
pos_g(:)=zero
!$omp end master
!$omp barrier


#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:pos_l) &
!$omp& private(i) reduction(+:pos_l)
#else
!$omp do reduction(+:pos_l)
#endif
do i=1,kmaxe
    pos_l(1)=pos_l(1)+ielem_filtered(i)
    pos_l(3)=pos_l(3)+ielem_er(i)
    pos_l(2)=pos_l(2)+ielem_er1(i)
    pos_l(4)=pos_l(4)+ielem_er2(i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


!$omp barrier
!$omp master

call mpi_allreduce(pos_l(1),pos_g(1),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
post1=(pos_g(1)/imaxe)*100.0

call mpi_allreduce(pos_l(2),pos_g(2),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
post2=(pos_g(2)/imaxe)

call mpi_allreduce(pos_l(3),pos_g(3),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
post3=(pos_g(3)/imaxe)

call mpi_allreduce(pos_l(4),pos_g(4),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
post4=(pos_g(4)/imaxe)

if (n.eq.0)then
    open(70,file='filtered.dat',form='formatted',action='write',position='append')
    write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,post1,post2,post3,post4
    close(70)
end if

call mpi_barrier(mpi_comm_world,ierror)


!$omp end master
!$omp barrier

end subroutine filtered_history



subroutine writesols(n)
implicit none
integer,intent(in)::n
integer::i,j,k,facex,ngp
real::pix
real:: leftv(1)
real::coord_gqp(1:2),rd

  pix=4.0d0*atan(1.0d0)

  do i=1,xmpielrank(n)
	write(500+n,*)"element id",ielem_ihexgl(i)
	do facex=1,ielem_ifca(i)
	write(500+n,*)"face",facex
		do ngp=1,qp_line
		write(500+n,*)"gqp",ngp
			!get gqp coords
			coord_gqp(1:2)= rec_qpoints_p(facex,ngp,1:2,i)



			!compute function

			    leftv(1)=0.0d0
if (sqrt(((coord_gqp(1)-0.25d0)**2)+((coord_gqp(2)-0.5d0)**2)).le.0.15)then
rd=(1.0d0/0.15d0)*sqrt(((coord_gqp(1)-0.25d0)**2)+((coord_gqp(2)-0.5d0)**2))

leftv(1)=0.25d0*(1.0d0+cos(pix*min(rd,1.0d0)))
end if

if (sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.25d0)**2)).le.0.15)then

rd=(1.0d0/0.15d0)*sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.25d0)**2))
leftv(1)=1.0d0-rd
end if

    if (sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.75d0)**2)).le.0.15)then

    rd=(1.0d0/0.15d0)*sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.75d0)**2))
	  if ((abs(coord_gqp(1)-0.5).ge.0.025d0).or.(coord_gqp(2).gt.0.85))then

	 leftv(1)=1.0d0
	  else

	  leftv(1)=0.0d0

	  end if
    end if

			write(500+n,*)leftv(1),rec_uleft(1,facex,ngp,i)
			u_e_val(1,1,i)=leftv(1)
			u_c_val(1,1,i)=rec_uleft(1,facex,ngp,i)

			end do
		end do
	end do











end subroutine writesols








subroutine gpu_to_host_vol
implicit none
#ifdef gpu
                  if (allocated(u_ct_val)) then
                  !$omp target update from(u_ct_val(1,:,:))
                  end if
                  if (allocated(u_c_val)) then
                  !$omp target update from(u_c_val(1,:,:))
                  end if
                  if (allocated(ielem_mood_o)) then
                  !$omp target update from(ielem_mood_o)
                  end if
                  if (allocated(ielem_reduce)) then
                  !$omp target update from(ielem_reduce)
                  end if
                  if (allocated(ielem_troubled)) then
                  !$omp target update from(ielem_troubled)
                  end if
                  if (allocated(ielem_vortex)) then
                  !$omp target update from(ielem_vortex)
                  end if
                  if (allocated(ielem_wcx)) then
                  !$omp target update from(ielem_wcx)
                  end if


#endif


end subroutine gpu_to_host_vol

subroutine gpu_to_host_vol_av
implicit none
#ifdef gpu
                  if (allocated(u_ct_val)) then
                  !$omp target update from(u_ct_val(ind1,:,:))
                  end if
                  if (allocated(u_c_val)) then
                  !$omp target update from(u_c_val(ind1,:,:))
                  end if
                  if (allocated(u_c_rms)) then
                  !$omp target update from(u_c_rms)
                  end if

#endif
end subroutine gpu_to_host_vol_av

subroutine gpu_to_host_surf
implicit none

#ifdef gpu
                  if (allocated(rec_grads)) then
                  !$omp target update from(rec_grads)
                  end if
                  if (itestcase.eq.4)then
                  if (allocated(rec_uleft)) then
                  !$omp target update from(rec_uleft)
                  end if
                  if (allocated(rec_uleft_dg)) then
                  !$omp target update from(rec_uleft_dg)
                  end if
                  if (allocated(rec_uleftturb)) then
                  !$omp target update from(rec_uleftturb)
                  end if
                  end if

#endif

end subroutine gpu_to_host_surf

subroutine gpu_to_host_surf_av
implicit none

#ifdef gpu
                  if (allocated(rec_gradsav)) then
                  !$omp target update from(rec_gradsav)
                  end if

#endif

end subroutine gpu_to_host_surf_av





end module io
