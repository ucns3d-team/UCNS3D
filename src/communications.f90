module communications
!> @brief
!> this module includes all the subroutines related to the mpi communications for various procedures
!> including exchange of halo cells from stencils, direct-side neighbours, and the boundary extrapolated
!> values for the conserved variables and their gradients at each surface/edge gaussian quadrature point

use mpiinfo
use declaration
use transform
implicit none

contains

subroutine renumber_neighbours(n,xmpie,xmpielrank,diexchanger,diexchanges)
!> @brief
!> this subroutine renumbers the neighbours indexing for cross referencing between different cpus
!> it is a process that is performed once the beginning of each run

implicit none
integer,intent(in)::n
integer,allocatable,dimension(:),intent(in)::xmpie
integer,allocatable,dimension(:),intent(in)::xmpielrank
type(aexchange),allocatable,dimension(:),intent(inout)::diexchanger,diexchanges
integer::i,j,k,l,m,e,kmaxe,ineedt,tneedt,ix1,ix2,cnbt,ivt,itax,inum_points,inn,jjj,kk,itogg,ifdn,ifdn2,ifdn3,iii,iouf,jj1
integer::c_n1,c_n2,c_n3,c_n4,d_n1,d_n2,d_n3,d_n4,kvf,itor4
ineedt=diexchanger(1)%tot
tneedt=diexchanges(1)%tot
kmaxe=xmpielrank(n)
do i=1,kmaxe
	if (ielem_interior(i).eq.1)then

	ielem_ineigh(:,i)=0;ielem_ineighn(:,i)=0;ielem_ineighb(:,i)=n
	else

	ielem_ineighn(:,i)=0;ielem_ineigh(:,i)=0;
	end if
end do


if (dimensiona.eq.3)then
 ifdn=4
else
  ifdn=2
end if


do i=1,kmaxe
      if (ielem_interior(i).eq.0)then
	do j=1,ielem_ifca(i)
		if (ielem_ineighg(j,i).gt.0)then
			if (xmpie(ielem_ineighg(j,i)).eq.n)then
			
			k=xmpil(ielem_ineighg(j,i))
			do ivt=1,ielem_ifca(k)
							if (ielem_ineighg(ivt,k).eq.ielem_ihexgl(i)) then
							ielem_ineigh(j,i)=xmpil(ielem_ineighg(j,i))
							ielem_ineighn(j,i)=ivt
							go to 101
							end if
			end do
			end if
		end if
	101 continue
		
	end do
      end if
end do






do i=1,kmaxe
if (ielem_interior(i).eq.1)then
do j=1,ielem_ifca(i)
itax=0
	if ((ielem_ineighg(j,i).gt.0))then
	    if(xmpie(ielem_ineighg(j,i)).ne.n)then
	      ielem_ineighb(j,i)=xmpie(ielem_ineighg(j,i))
	     
	    else
			      
				    k=xmpil(ielem_ineighg(j,i))
    ! 				  
				    do ivt=1,ielem_ifca(k)
							  if (ielem_ineighg(ivt,k).eq.ielem_ihexgl(i)) then
    ! 							
							    ielem_ineigh(j,i)=k
							    ielem_ineighn(j,i)=ivt
							    ielem_ineighb(j,i)=n
							    end if
			    
				    end do  

	    end if
	end if
end do
end if
end do












jj1=0
do k=1,ineedt
do ix1=1,tneedt
do e=1,diexchanger(k)%muchineed(1)
do ix2=1,diexchanges(ix1)%muchtheyneed(1)


if ((diexchanges1(ix1)%whattheyneed(ix2).eq.diexchanger1(k)%sideineedn(e)).and.&
(diexchanges1(ix1)%sidetheyneedn(ix2).eq.diexchanger1(k)%whatineed(e)))then		!if 1

i=xmpil(diexchanger1(k)%sideineedn(e))
! do i=1,kmaxe
  if (ielem_interior(i).eq.1)then												!if 2
    do j=1,ielem_ifca(i)
	  if (ielem_ineighb(j,i).ne.n)then										!if 3
	  
	    if (ielem_ibounds(j,i).gt.0)then										!if 4
		if ((ibound_icode(ielem_ibounds(j,i)).eq.5).or.(ibound_icode(ielem_ibounds(j,i)).eq.50))then		!if 5
		    if (dimensiona.eq.3)then											!if 6
				if ( ielem_types_faces(j,i).eq.5)then							!if 7
					inum_points=qp_quad_n
				else
					inum_points=qp_triangle_n
				end if
			else
					inum_points=qp_line_n
			end if																!end if6
		    
		   
		    



if ((diexchanger(k)%procid.eq.ielem_ineighb(j,i)).and. (diexchanges(ix1)%procid.eq.ielem_ineighb(j,i)))then		!if
if ((diexchanges1(ix1)%sidetheyneedn(ix2).eq.ielem_ineighg(j,i)).and.&
(ielem_ineighg(j,i).eq.diexchanger1(k)%whatineed(e)))then
 do inn=1,inum_points
if ((diexchanges1(ix1)%whattheyneed(ix2).eq.ielem_ihexgl(i)).and.(diexchanges1(ix1)%qtheyneed(ix2).eq.diexchanger1(k)%qineed(e)).and.(inn.eq.diexchanges1(ix1)%qtheyneed(ix2)))then

ielem_qface(j,inn,ielem_inter_id(ielem_indexf(i)))=e

! ielem_q_face_q_mapl(inn,j,i)=e
ielem_ineighn(j,i)=k
diexchanges(ix1)%sidetheyneed(ix2)=j

jj1=jj1+1
end if
end do
end if
end if



      

      
end if      
else

if (dimensiona.eq.3)then
			if ( ielem_types_faces(j,i).eq.5)then
			    inum_points=qp_quad_n
			  else
			    inum_points=qp_triangle_n
			 end if
		    else
			    inum_points=qp_line_n
		    end if
		    
		    



if ((diexchanger(k)%procid.eq.ielem_ineighb(j,i)).and. (diexchanges(ix1)%procid.eq.ielem_ineighb(j,i)))then
if ((diexchanges1(ix1)%sidetheyneedn(ix2).eq.ielem_ineighg(j,i)).and.&
(ielem_ineighg(j,i).eq.diexchanger1(k)%whatineed(e)))then
do inn=1,inum_points
if ((diexchanges1(ix1)%whattheyneed(ix2).eq.ielem_ihexgl(i)).and.(diexchanges1(ix1)%qtheyneed(ix2).eq.diexchanger1(k)%qineed(e)).and.(inn.eq.diexchanges1(ix1)%qtheyneed(ix2)))then

ielem_qface(j,inn,ielem_inter_id(ielem_indexf(i)))=e

! ielem_q_face_q_mapl(inn,j,i)=e
ielem_ineighn(j,i)=k
diexchanges(ix1)%sidetheyneed(ix2)=j

end if
end do
end if
end if


end if
end if				!if 3


end do
end if					!if 2


! end do
end if			!if 1
end do
end do
end do
end do





end subroutine renumber_neighbours



subroutine solex_alloc(n)
!> @brief
!> this subroutine allocates the memory for the halo cells of the direct side neighbours only
!> and is therefore used for lower-order schemes 
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,e,kmaxe,ineedt,tneedt,icpuid
real,dimension(1:dimensiona)::cords
ineedt=diexchanger(1)%tot
tneedt=diexchanges(1)%tot
allocate(dsolchanger(ineedt))
allocate(dsolchanges(tneedt))
! dsolchanger(:)=0
! dsolchanges(:)=0
do i=1,ineedt
	dsolchanger(i)%procid=diexchanger(i)%procid
	allocate (dsolchanger(i)%centres(diexchanger(i)%muchineed(1),dims))
	
	dsolchanger(i)%centres(:,:)=0.d0
	
	allocate (dsolchanger(i)%sol(diexchanger(i)%muchineed(1),nof_variables+turbulenceequations+passivescalar))
	dsolchanger(i)%sol(:,:)=0.0d0
end do
do i=1,tneedt
	dsolchanges(i)%procid=diexchanges(i)%procid
	allocate (dsolchanges(i)%centres(diexchanges(i)%muchtheyneed(1),dims))
	
	dsolchanges(i)%centres(:,:)=0.0d0
	
	allocate (dsolchanges(i)%sol(diexchanges(i)%muchtheyneed(1),nof_variables+turbulenceequations+passivescalar))
	dsolchanges(i)%sol(:,:)=0.0d0
end do

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
do i=1,tneedt
	do k=1,diexchanges(i)%muchtheyneed(1)
		j=diexchanges(i)%localref(k)
		    if (dimensiona.eq.3)then
		    call compute_centre3d(j,cords)
		    else
		     call compute_centre2d(j,cords)
		    end if
	dsolchanges(i)%centres(k,1:dims)=cords(1:dims)
	end do
end do
call mpi_barrier(mpi_comm_world,ierror)
icpuid=n
do i=1,ineedt
		do k=1,tneedt
		if (dsolchanger(i)%procid.eq.dsolchanges(k)%procid)then
		call mpi_sendrecv(dsolchanges(k)%centres(1:diexchanges(k)%muchtheyneed(1),1:dims),&
diexchanges(k)%muchtheyneed(1)*dims,mpi_double_precision,dsolchanger(i)%procid,&
		dsolchanges(k)%procid,dsolchanger(i)%centres(1:diexchanger(i)%muchineed(1),1:dims),&
diexchanger(i)%muchineed(1)*dims,mpi_double_precision,&
		dsolchanger(i)%procid,icpuid,mpi_comm_world,status,ierror)
		end if
		end do
end do

end subroutine solex_alloc




subroutine estabexhange(n,imaxe,xmpie,xmpin,xmpielrank,ilocalstencil,diexchanger,&
diexchanges,direcexr,direcexs,numneighbours,ischeme,isize,iperiodicity,typesten,xmpil)
!> @brief
!> this subroutine is establishing the communication patterns for all the mpi processes and in particular for
!> the halo cells of the reconstruction stencils and it must be noted that each cell can have a different number of stencils
!> of different size

	implicit none
	integer,intent(in)::n,imaxe,isize,iperiodicity,typesten
	integer,allocatable,dimension(:),intent(in)::xmpie,xmpin,xmpil
	integer,allocatable,dimension(:),intent(in)::xmpielrank
	type(aexchange),allocatable,dimension(:),intent(inout)::diexchanges
	type(aexchange),allocatable,dimension(:),intent(inout)::diexchanger
	integer,allocatable,dimension(:,:,:,:),intent(in)::ilocalstencil
	integer,intent(in)::numneighbours,ischeme
	type(arecex),allocatable,dimension(:),intent(inout)::direcexr		!receive elements due to stencils
	type(arecex),allocatable,dimension(:),intent(inout)::direcexs
	integer,allocatable,dimension(:)::totstenc,totstenc2
	integer,allocatable,dimension(:)::stenred
	integer,allocatable,dimension(:)::stenredproc
	type aliststencils
	integer::procid
	integer,allocatable,dimension(:)::listsarray
	end type aliststencils
	integer,dimension(0:isize-1)::sendwh
	integer,dimension(0:isize-1)::listofpr
	integer,dimension(0:isize-1)::sendstwh
	type:: alistgog
	integer::procid,imuch,iavc,iavt
	integer,allocatable,dimension(:)::globarray
	integer,allocatable,dimension(:,:)::nodex
	end type alistgog
	type(aliststencils),allocatable,dimension(:)::diliststen
	type(alistgog),allocatable,dimension(:)::dilistgog,dilistside,dilistglo,dilistq
	integer::howmany,itee,source,tag,cpuidrec,icpuid,ixflag,iteedum,iavc,iavt
	integer,dimension(1:1)::rumts,sumts
	real,dimension(1:1)::dumts
	integer::stencounter,stencounter2,stensame,iiflag,tempsten2,iiflag2,kxk,ix,icx,icf,ifdn,ifdn2,INDEX_INT
	integer,dimension(0:isize-1)::ihmste
	integer::i,j,ji,k,lm,kmaxn,kk,kmaxe,iaa,l,ingt,ik,tempint,temp2,itarget,isty5,i86,i87,i65,inum_points,iin
	kmaxe=xmpielrank(n)
	listofpr(:)=0
	i65=0
	if (dimensiona.eq.3)then
	ifdn=4
	
	else

	ifdn=2
	end if
	index_int=0
	do k=1,kmaxe
	      if (ielem_interior(k).eq.1)then
	      index_int=index_int+1
		  ielem_indexf(k)=index_int
	      end if
	end do

	allocate(ielem_inter_id(index_int));ielem_inter_id=0






	if (dimensiona.eq.3)then
	allocate(ielem_qface(max_faces,numberofpoints2,index_int));ielem_qface=0
	else
	allocate(ielem_qface(max_faces,numberofpoints2,index_int));ielem_qface=0
	end if

	index_int=0
	do k=1,kmaxe
	      if (ielem_interior(k).eq.1)then
	      index_int=index_int+1
	      ielem_inter_id(ielem_indexf(k))=index_int

	      end if
	end do







	do k=1,kmaxe
		if (ielem_interior(k).eq.1)then
			do l=1,ielem_ifca(k)
				    if (ielem_ineighg(l,k).gt.0)then
					    if (xmpie(ielem_ineighg(l,k)).ne.n)then





					    if (dimensiona.eq.3)then
									if ( ielem_types_faces(l,k).eq.5)then
									inum_points=qp_quad_n;
									else
									inum_points=qp_triangle_n ;
									end if
								else
								inum_points=qp_line_n;
								end if




					    
					   
					      
					    listofpr(xmpie(ielem_ineighg(l,k)))=listofpr(xmpie(ielem_ineighg(l,k)))+inum_points
					    end if
				    end if
			  end do
		
		end if
	end do

	howmany=0
	do k=0,isize-1
		if (listofpr(k).gt.0)then
		howmany=howmany+1
		end if
	end do
	allocate (dilistgog(howmany))
	allocate (dilistside(howmany))
	allocate (dilistq(howmany))
	allocate (dilistglo(howmany))
	
! 	dilistgog(:)=0
! 	dilistside(:)=0
	kk=0
	do k=0,isize-1
		if (listofpr(k).gt.0)then
		kk=kk+1
		dilistgog(kk)%procid=k
		dilistgog(kk)%imuch=listofpr(k)
		dilistside(kk)%procid=k
		dilistside(kk)%imuch=listofpr(k)
		dilistq(kk)%procid=k
		dilistq(kk)%imuch=listofpr(k)
		allocate(dilistgog(kk)%globarray(listofpr(k)))
		allocate(dilistside(kk)%globarray(listofpr(k)))
		allocate(dilistside(kk)%nodex(listofpr(k),ifdn))
		allocate(dilistglo(kk)%globarray(listofpr(k)))
		allocate(dilistq(kk)%globarray(listofpr(k)))
		dilistgog(kk)%globarray(:)=0
		dilistside(kk)%globarray(:)=0;dilistside(kk)%nodex(:,:)=0
		dilistglo(kk)%globarray(:)=0
		dilistq(kk)%globarray(:)=0
		end if
	end do
	kk=0
	do ingt=0,isize-1
		if (listofpr(ingt).gt.0)then
! 		
		kk=kk+1
		iaa=0
		do k=1,kmaxe
		      if (ielem_interior(k).eq.1)then
			do l=1,ielem_ifca(k)
				    if (ielem_ineighg(l,k).gt.0)then
					    if (xmpie(ielem_ineighg(l,k)).ne.n)then
						if (xmpie(ielem_ineighg(l,k)).eq.dilistgog(kk)%procid)then
						     if (dimensiona.eq.3)then
									if ( ielem_types_faces(l,k).eq.5)then
									inum_points=qp_quad_n; ifdn2=4
									else
									inum_points=qp_triangle_n ; ifdn2=3
									end if
								else
								inum_points=qp_line_n; ifdn2=2
								end if
						 
							    do iin=iaa+1,iaa+inum_points
							    dilistgog(kk)%globarray(iin)=(ielem_ineighg(l,k))
							    dilistglo(kk)%globarray(iin)=(ielem_ihexgl(k))
							    dilistside(kk)%globarray(iin)=l
							    dilistq(kk)%globarray(iin)=iin-iaa
							    dilistside(kk)%nodex(iin,1:ifdn2)=ielem_nodes_faces(l,1:ifdn2,k)
							    
							    end do
							    iaa=iaa+inum_points
						      
						  end if
					    end if
				      end if
			end do
		      end if
		end do

	      end if
	  end do
	      


	sendwh(:)=0	
	allocate (diexchanger(howmany))
	allocate (diexchanger1(howmany))
! 	diexchanger(:)=0
	diexchanger(:)%tot=howmany
	diexchanger1(:)%tot=howmany

	kk=0
	do k=0,isize-1
		if (listofpr(k).gt.0)then
		kk=kk+1
		diexchanger(kk)%procid=k
		end if
	end do
	
	do i=1,howmany
		if (diexchanger(i)%procid.ne.n)then
		allocate(diexchanger(i)%muchineed(1))
		diexchanger(i)%muchineed(:)=0
		diexchanger(i)%muchineed(1)=listofpr(diexchanger(i)%procid)
		allocate(diexchanger1(i)%whatineed(diexchanger(i)%muchineed(1)))	!element number global
		allocate(diexchanger1(i)%localref(diexchanger(i)%muchineed(1)))
		allocate(diexchanger1(i)%sideineed(diexchanger(i)%muchineed(1)))
		allocate(diexchanger1(i)%nodex(diexchanger(i)%muchineed(1),ifdn))
		allocate(diexchanger1(i)%qineed(diexchanger(i)%muchineed(1)))
		allocate(diexchanger1(i)%sideineedn(diexchanger(i)%muchineed(1)))
! 		allocate(diexchanger1(i)%qineed(diexchanger(i)%muchineed(1)))
		diexchanger1(i)%whatineed(:)=0	!element number global
		diexchanger1(i)%localref(:)=0
		diexchanger1(i)%sideineed(:)=0;diexchanger1(i)%nodex(:,:)=0
		diexchanger1(i)%sideineedn(:)=0
 		diexchanger1(i)%qineed(:)=0
		
		end if
	end do
	tempint=0
! 	
	do i=1,howmany
		if (diexchanger(i)%procid.ne.n)then
		tempint=tempint+1
		if (diexchanger(i)%procid.eq.dilistgog(tempint)%procid)then
		diexchanger1(i)%whatineed(:)=dilistgog(tempint)%globarray(:)
		diexchanger1(i)%sideineedn(:)=dilistglo(tempint)%globarray(:)
		diexchanger1(i)%nodex(1:diexchanger(i)%muchineed(1),1:ifdn)=dilistside(tempint)%nodex(1:diexchanger(i)%muchineed(1),1:ifdn)
		diexchanger1(i)%qineed(1:diexchanger(i)%muchineed(1))=dilistq(tempint)%globarray(1:diexchanger(i)%muchineed(1))


		
		
		end if
		end if
	end do
	deallocate(dilistgog)
	call mpi_barrier(mpi_comm_world,ierror)
	icpuid=n
	do i=1,howmany
		if (diexchanger(i)%procid.ne.n)then
		call mpi_sendrecv(diexchanger(i)%muchineed(1),1,mpi_integer,diexchanger(i)%procid,diexchanger(i)%procid,&
		itee,1,mpi_integer,diexchanger(i)%procid,icpuid,mpi_comm_world,status,ierror)
		sendwh(diexchanger(i)%procid)=itee
		end if
	end do
	call mpi_barrier(mpi_comm_world,ierror)
	temp2=0
	do i=0,isize-1
		if (sendwh(i).gt.0)then
		temp2=temp2+1
		end if
	end do
	
! 
	allocate (diexchanges(temp2))
	allocate (diexchanges1(temp2))
! 	diexchanges(:)=0
	diexchanges(:)%tot=temp2
	kk=0
	do i=0,isize-1
		if (sendwh(i).gt.0)then
		kk=kk+1
		diexchanges(kk)%procid=i
		allocate(diexchanges(kk)%muchtheyneed(1))
		diexchanges(kk)%muchtheyneed(:)=0
		diexchanges(kk)%muchtheyneed(1)=sendwh(i)
		allocate(diexchanges1(kk)%whattheyneed(diexchanges(kk)%muchtheyneed(1)))
		allocate(diexchanges(kk)%localref(diexchanges(kk)%muchtheyneed(1)))
		allocate(diexchanges(kk)%sidetheyneed(diexchanges(kk)%muchtheyneed(1)))

		
		allocate(diexchanges1(kk)%sidetheyneedn(diexchanges(kk)%muchtheyneed(1)))
		allocate(diexchanges1(kk)%nodex(diexchanges(kk)%muchtheyneed(1),ifdn))
		allocate(diexchanges1(kk)%qtheyneed(diexchanges(kk)%muchtheyneed(1)))
		allocate(diexchanges(kk)%qtheyneed(diexchanges(kk)%muchtheyneed(1)))
		diexchanges1(kk)%whattheyneed(:)=0
		diexchanges(kk)%localref(:)=0
		diexchanges(kk)%sidetheyneed(:)=0
		diexchanges1(kk)%sidetheyneedn(:)=0;diexchanges1(kk)%nodex(:,:)=0
		diexchanges1(kk)%qtheyneed(:)=0
		diexchanges(kk)%qtheyneed(:)=0
		end if
	end do
!
	do i=1,howmany
		do k=1,temp2
		if (diexchanges(k)%procid.eq.dilistside(i)%procid)then
		diexchanges(k)%sidetheyneed(:)=dilistside(i)%globarray(:)
!  		diexchanges(k)%qtheyneed(:)=dilistq(i)%globarray(:)
		end if
		end do
	end do
 	call mpi_barrier(mpi_comm_world,ierror)

	do i=1,howmany
		do k=1,temp2
		if (diexchanger(i)%procid.eq.diexchanges(k)%procid)then
		call mpi_sendrecv(diexchanger1(i)%whatineed(1:diexchanger(i)%muchineed(1)),&
diexchanger(i)%muchineed(1),mpi_integer,diexchanger(i)%procid,&
		diexchanger(i)%procid,diexchanges1(k)%whattheyneed(1:diexchanges(k)%muchtheyneed(1)),&
diexchanges(k)%muchtheyneed(1),mpi_integer,diexchanger(i)%procid,&
		icpuid,mpi_comm_world,status,ierror)
		end if
		end do
	end do
	do i=1,howmany
		do k=1,temp2
		if (diexchanger(i)%procid.eq.diexchanges(k)%procid)then
		call mpi_sendrecv(diexchanger1(i)%nodex(1:diexchanger(i)%muchineed(1),1:ifdn),&
diexchanger(i)%muchineed(1)*ifdn,mpi_integer,diexchanger(i)%procid,&
		diexchanger(i)%procid,diexchanges1(k)%nodex(1:diexchanges(k)%muchtheyneed(1),1:ifdn),&
diexchanges(k)%muchtheyneed(1)*ifdn,mpi_integer,diexchanger(i)%procid,&
		icpuid,mpi_comm_world,status,ierror)
		end if
		end do
	end do
	do i=1,howmany
		do k=1,temp2
		if (diexchanger(i)%procid.eq.diexchanges(k)%procid)then
		call mpi_sendrecv(diexchanger1(i)%qineed(1:diexchanger(i)%muchineed(1)),&
diexchanger(i)%muchineed(1),mpi_integer,diexchanger(i)%procid,&
		diexchanger(i)%procid,diexchanges1(k)%qtheyneed(1:diexchanges(k)%muchtheyneed(1)),&
diexchanges(k)%muchtheyneed(1),mpi_integer,diexchanger(i)%procid,&
		icpuid,mpi_comm_world,status,ierror)
		diexchanges(k)%qtheyneed(1:diexchanges(k)%muchtheyneed(1))=diexchanges1(k)%qtheyneed(1:diexchanges(k)%muchtheyneed(1))
		
		
		end if
		end do
	end do

	call mpi_barrier(mpi_comm_world,ierror)
	do i=1,howmany
		do k=1,temp2
		if (diexchanger(i)%procid.eq.diexchanges(k)%procid)then
		call mpi_sendrecv(diexchanger1(i)%sideineedn(1:diexchanger(i)%muchineed(1)),&
diexchanger(i)%muchineed(1),mpi_integer,diexchanger(i)%procid,&
		diexchanger(i)%procid,diexchanges1(k)%sidetheyneedn(1:diexchanges(k)%muchtheyneed(1)),&
diexchanges(k)%muchtheyneed(1),mpi_integer,diexchanger(i)%procid,&
		icpuid,mpi_comm_world,status,ierror)
		end if
		end do
	end do
!make the ones below q_points!
 	call mpi_barrier(mpi_comm_world,ierror)
	do i=1,howmany
		do k=1,temp2
		if (diexchanger(i)%procid.eq.diexchanges(k)%procid)then
		call mpi_sendrecv(diexchanges(k)%sidetheyneed(1:diexchanges(k)%muchtheyneed(1)),&
diexchanges(k)%muchtheyneed(1),mpi_integer,diexchanges(k)%procid,&
		diexchanges(k)%procid,diexchanger1(i)%sideineed(1:diexchanger(i)%muchineed(1)),&
diexchanger(i)%muchineed(1),mpi_integer,diexchanges(k)%procid,&
		icpuid,mpi_comm_world,status,ierror)

		
		end if
		end do
	end do
      
	call mpi_barrier(mpi_comm_world,ierror)
	
	
      
	

	deallocate(dilistside,dilistq)
	do i=1,temp2
		do k=1,diexchanges(i)%muchtheyneed(1)

 			if (xmpie(diexchanges1(i)%whattheyneed(k)).eq.n)then

 			diexchanges(i)%localref(k)=xmpil(diexchanges1(i)%whattheyneed(k))
   			end if

			
		end do
	end do
	do i=1,howmany
		do k=1,temp2
		if (diexchanger(i)%procid.eq.diexchanges(k)%procid)then
		call mpi_sendrecv(diexchanges(k)%localref(1:diexchanges(k)%muchtheyneed(1)),&
diexchanges(k)%muchtheyneed(1),mpi_integer,diexchanges(k)%procid,&
		diexchanges(k)%procid,diexchanger1(i)%localref(1:diexchanger(i)%muchineed(1)),&
diexchanger(i)%muchineed(1),mpi_integer,diexchanger(i)%procid,&
		icpuid,mpi_comm_world,status,ierror)
		end if
		end do
	end do
	call mpi_barrier(mpi_comm_world,ierror)

	
	if (ischeme.gt.1)then
	stencounter=0
!
  	
	allocate(totstenc2(imaxe))
	totstenc2=0
	

	do i=1,kmaxe
		do k=1,typesten
			isty5=0
			if ((k.eq.1).or.(ees.ne.5))then
			itarget=ielem_inumneighbours(i)
			else
			itarget=numneighbours2
			end if
			
			
			do kk=1,itarget
				if(ilocalstencil(n,i,k,kk).gt.0)then
				isty5=isty5+1
				end if
			end do
			if (isty5.eq.itarget)then
			do kk=1,itarget
					
				iiflag=0

					if (xmpie(ilocalstencil(n,i,k,kk)).eq.n) then

						iiflag=1
					end if

					
					if (iiflag.eq.0)then
					totstenc2((ilocalstencil(n,i,k,kk)))=ilocalstencil(n,i,k,kk)
					end if
			end do
			end if
		end do
	end do

	
	
	
	kxk=0
	do i=1,imaxe
		if (totstenc2(i).gt.0)then
		kxk=kxk+1
		end if
	end do

	allocate (stenred(kxk))
	allocate (stenredproc(kxk))

	stencounter=kxk

	allocate(totstenc(stencounter))
	totstenc=0
	kxk=0
	do i=1,imaxe
		if (totstenc2(i).gt.0)then
		kxk=kxk+1
		totstenc(kxk)=totstenc2(i)
		end if
	end do



	deallocate(totstenc2)




	stenred=0
	stenredproc=0
	ihmste=0
	kxk=0
	do i=1,stencounter
		if (totstenc(i).gt.0)then
		kxk=kxk+1
		stenred(kxk)=totstenc(i)
		end if
	end do
	deallocate (totstenc)




	do i=1,kxk
		stenredproc(i)=xmpie(stenred(i))
	end do
	do i=0,isize-1
		ihmste(i)=0
		if (i.ne.n)then
		do ix=1,kxk
		if(stenredproc(ix).eq.i)then
		ihmste(i)=ihmste(i)+1
		end if
		end do
		end if
! 		
	end do
	ix=0
	do i=0,isize-1
		if (i.ne.n)then
		if (ihmste(i).gt.0)then
		ix=ix+1
		end if
		end if
	end do
	call mpi_barrier(mpi_comm_world,ierror)
	allocate (diliststen(ix))
! 	diliststen(:)=0
	ix=0
	do i=0,isize-1
		if (i.ne.n)then
		if (ihmste(i).gt.0)then
		ix=ix+1
		diliststen(ix)%procid=i
		allocate (diliststen(ix)%listsarray(ihmste(i)))
		diliststen(ix)%listsarray(:)=0
		icf=0
		do icx=1,kxk
		if (stenredproc(icx).eq.i)then
		icf=icf+1
		diliststen(ix)%listsarray(icf)=stenred(icx)
		end if
		end do

		end if
		end if
	end do
	allocate (direcexr(ix))
	allocate (direcexr1(ix))
! 	direcexr(:)=0
	direcexr(:)%tot=0
	direcexr(:)%tot=ix
	icx=0
	do i=0,isize-1
		if (ihmste(i).gt.0)then
		icx=icx+1
		direcexr(icx)%procid=i
		allocate(direcexr(icx)%muchineed(1))
		direcexr(icx)%muchineed(:)=0
		direcexr(icx)%muchineed(1)=ihmste(i)
		allocate(direcexr1(icx)%whatineed(direcexr(icx)%muchineed(1)))
		allocate(direcexr1(icx)%localref(direcexr(icx)%muchineed(1)))
		allocate(direcexr1(icx)%ishape(direcexr(icx)%muchineed(1)))
		direcexr1(icx)%whatineed(:)=0
		direcexr1(icx)%localref(:)=0
		direcexr1(icx)%ishape(:)=0
! 		
		do icf=1,ix
		if (diliststen(icf)%procid.eq.direcexr(icx)%procid)then
		direcexr1(icx)%whatineed(:)=diliststen(icf)%listsarray(:)
! 		
		end if
		end do
		end if
	end do
	

	
	sendstwh(:)=0
	icpuid=n
	do k=0,isize-1
		ixflag=0
		iteedum=0
		itee=0
		if (k.ne.n)then
		do i=1,icx
			if (direcexr(i)%procid.eq.k)then
			call mpi_sendrecv(direcexr(i)%muchineed(1),1,mpi_integer,direcexr(i)%procid,icpuid,&
			itee,1,mpi_integer,direcexr(i)%procid,direcexr(i)%procid,mpi_comm_world,status,ierror)
			sendstwh(direcexr(i)%procid)=itee
			ixflag=1
			end if
		end do
		if (ixflag.eq.0)then
			call mpi_sendrecv(iteedum,1,mpi_integer,k,icpuid,&
			itee,1,mpi_integer,k,k,mpi_comm_world,status,ierror)
			sendstwh(k)=itee
		end if
		end if
	end do	
	call mpi_barrier(mpi_comm_world,ierror)
	temp2=0
	do i=0,isize-1
		if (sendstwh(i).gt.0)then
		temp2=temp2+1
		end if
	end do

	call mpi_barrier(mpi_comm_world,ierror)
	allocate (direcexs(temp2))
	allocate (direcexs1(temp2))
! 	direcexs(:)=0
	direcexs(:)%tot=temp2
	kk=0
	

	
	do i=0,isize-1
		if (sendstwh(i).gt.0)then
		kk=kk+1
		direcexs(kk)%procid=i
		allocate(direcexs(kk)%muchtheyneed(1))
		direcexs(kk)%muchtheyneed(:)=0
		direcexs(kk)%muchtheyneed(1)=sendstwh(i)
		allocate(direcexs1(kk)%whattheyneed(direcexs(kk)%muchtheyneed(1)))
		allocate(direcexs(kk)%localref(direcexs(kk)%muchtheyneed(1)))
		allocate(direcexs1(kk)%ishape(direcexs(kk)%muchtheyneed(1)))
		direcexs1(kk)%whattheyneed(:)=0
		direcexs(kk)%localref(:)=0
		direcexs1(kk)%ishape(:)=0
		end if
	end do
	
			rumts(1:1)=0
	
	call mpi_barrier(mpi_comm_world,ierror)
	
	

	
	do i=0,isize-1
			if (i.ne.n) then
				do j=1,temp2
					iavt=10000
					if (direcexs(j)%procid.eq.i)then
					iavt=j 
					go to 7001
					end if
				end do
				7001 continue
				do k=1,icx
					iavc=10000
					if (direcexr(k)%procid.eq.i) then
					iavc=k
					go to 8001
					end if
				end do
				8001 continue
if ((iavc.eq.10000).and.(iavt.ne.10000)) then
call mpi_sendrecv(sumts(1:1),1,mpi_integer,i,icpuid,&
direcexs1(iavt)%whattheyneed(1:direcexs(iavt)%muchtheyneed(1)),direcexs(iavt)%muchtheyneed(1),&
mpi_integer,direcexs(iavt)%procid,direcexs(iavt)%procid,mpi_comm_world,status,ierror)
end if
if ((iavc.eq.10000).and.(iavt.eq.10000)) then
call mpi_sendrecv(sumts(1:1),1,mpi_integer,i,icpuid,&
rumts(1:1),1,mpi_integer,i,i,mpi_comm_world,status,ierror)
end if
 if ((iavc.ne.10000).and.(iavt.eq.10000)) then
 call mpi_sendrecv(direcexr1(iavc)%whatineed(1:direcexr(iavc)%muchineed(1)),direcexr(iavc)%muchineed(1),&
 mpi_integer,direcexr(iavc)%procid,icpuid,rumts(1:1),1,mpi_integer,i,i,mpi_comm_world,status,ierror)
 end if
if ((iavc.ne.10000).and.(iavt.ne.10000)) then
call mpi_sendrecv(direcexr1(iavc)%whatineed(1:direcexr(iavc)%muchineed(1)),direcexr(iavc)%muchineed(1),&
mpi_integer,direcexr(iavc)%procid,icpuid,&
direcexs1(iavt)%whattheyneed(1:direcexs(iavt)%muchtheyneed(1)),direcexs(iavt)%muchtheyneed(1),mpi_integer,&
direcexs(iavt)%procid,direcexs(iavt)%procid,mpi_comm_world,status,ierror)
end if
			!	end do !k
			!end do! j
			end if	! i.ne.n
	end do
	call mpi_barrier(mpi_comm_world,ierror)
! 	
	   	
	do i=1,temp2
		do k=1,direcexs(i)%muchtheyneed(1)
			!do j=1,kmaxe
			if (xmpie(direcexs1(i)%whattheyneed(k)).eq.n) then
			j=xmpil(direcexs1(i)%whattheyneed(k))
			!if (direcexs(i)%whattheyneed(k).eq.ielem_ihexgl(j))then
			direcexs(i)%localref(k)=j
			direcexs1(i)%ishape(k)=ielem_ishape(j)
			
			end if
			!end do
		end do
	end do
	call mpi_barrier(mpi_comm_world,ierror)
		do i=0,isize-1
			rumts=0
			sumts=0
			if (i.ne.n) then
				do j=1,temp2
					iavt=10000
					if (direcexs(j)%procid.eq.i)then
					iavt=j 
					go to 9001
					end if
				end do
				9001 continue
				do k=1,icx
					iavc=10000
					if (direcexr(k)%procid.eq.i) then
					iavc=k
					go to 10001
					end if
				end do
				10001 continue
			
  if ((iavt.eq.10000).and.(iavc.ne.10000)) then
  call mpi_sendrecv(sumts(1:1),1,mpi_integer,i,icpuid,&
  direcexr1(iavc)%localref(1:direcexr(iavc)%muchineed(1)),direcexr(iavc)%muchineed(1),&
mpi_integer,direcexr(iavc)%procid,direcexr(iavc)%procid,mpi_comm_world,status,ierror)
  call mpi_sendrecv(sumts(1:1),1,mpi_integer,i,icpuid,&
  direcexr1(iavc)%ishape(1:direcexr(iavc)%muchineed(1)),direcexr(iavc)%muchineed(1),&
mpi_integer,direcexr(iavc)%procid,direcexr(iavc)%procid,mpi_comm_world,status,ierror)
  end if
   if ((iavt.eq.10000).and.(iavc.eq.10000)) then
   call mpi_sendrecv(sumts(1:1),1,mpi_integer,i,icpuid,&
   rumts(1:1),1,mpi_integer,i,i,mpi_comm_world,status,ierror)
   call mpi_sendrecv(sumts(1:1),1,mpi_integer,i,icpuid,&
   rumts(1:1),1,mpi_integer,i,i,mpi_comm_world,status,ierror)
   end if
  if ((iavt.ne.10000).and.(iavc.eq.10000)) then
  call mpi_sendrecv(direcexs(iavt)%localref(1:direcexs(iavt)%muchtheyneed(1)),&
direcexs(iavt)%muchtheyneed(1),mpi_integer,direcexs(iavt)%procid,icpuid,&
  rumts(1:1),1,mpi_integer,i,i,mpi_comm_world,status,ierror)
  call mpi_sendrecv(direcexs1(iavt)%ishape(1:direcexs(iavt)%muchtheyneed(1)),&
  direcexs(iavt)%muchtheyneed(1),mpi_integer,direcexs(iavt)%procid,icpuid,&
  rumts(1:1),1,mpi_integer,i,i,mpi_comm_world,status,ierror)
  end if
  if ((iavt.ne.10000).and.(iavc.ne.10000)) then
  call mpi_sendrecv(direcexs(iavt)%localref(1:direcexs(iavt)%muchtheyneed(1)),&
direcexs(iavt)%muchtheyneed(1),mpi_integer,direcexs(iavt)%procid,icpuid,direcexr1(iavc)%localref(1:direcexr(iavc)%muchineed(1)),&
direcexr(iavc)%muchineed(1),mpi_integer,direcexr(iavc)%procid,direcexr(iavc)%procid,mpi_comm_world,status,ierror)
  call mpi_sendrecv(direcexs1(iavt)%ishape(1:direcexs(iavt)%muchtheyneed(1)),direcexs(iavt)%muchtheyneed(1),&
  mpi_integer,direcexs(iavt)%procid,icpuid,direcexr1(iavc)%ishape(1:direcexr(iavc)%muchineed(1)),direcexr(iavc)%muchineed(1)&
  ,mpi_integer,direcexr(iavc)%procid,direcexr(iavc)%procid,mpi_comm_world,status,ierror)
  end if
			
			end if	! i.ne.n
		end do
	call mpi_barrier(mpi_comm_world,ierror)
	deallocate (stenred)
	deallocate (stenredproc)
	deallocate (diliststen)
	end if	
      
	
	deallocate (dilistglo)
	




	!temporary allocate

	
end subroutine estabexhange

subroutine exchange_lower(n)
!> @brief
!> this subroutine is exchanging the variables of all the halo cells for direct side neighbours only
!> and is primarily used by lower order schemes

implicit none
integer,intent(in)::n
integer::i,k,ineedt,tneedt,icpuid,itest
ineedt=diexchanger(1)%tot
tneedt=diexchanges(1)%tot

  itest=nof_variables+turbulenceequations+passivescalar

  
if (( turbulence .gt. 0 ).or.(passivescalar.gt.0)) then
!$omp do
do i=1,tneedt
	do k=1,diexchanges(i)%muchtheyneed(1)
	      dsolchanges(i)%sol(k,1:nof_variables)=u_c_val(1,1:nof_variables,diexchanges(i)%localref(k))
	      dsolchanges(i)%sol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,diexchanges(i)%localref(k))
	end do
end do
!$omp end do

else
!$omp do
do i=1,tneedt
	do k=1,diexchanges(i)%muchtheyneed(1)
	      dsolchanges(i)%sol(k,1:nof_variables)=u_c_val(1,1:nof_variables,diexchanges(i)%localref(k))
	end do
end do
!$omp end do
end if




 
!$omp barrier
!$omp master
icpuid=n
do i=1,ineedt
		do k=1,tneedt
		if (dsolchanger(i)%procid.eq.dsolchanges(k)%procid)then
		call mpi_sendrecv(dsolchanges(k)%sol(1:diexchanges(k)%muchtheyneed(1),1:itest),&
diexchanges(k)%muchtheyneed(1)*itest,mpi_double_precision,dsolchanger(i)%procid,&
		dsolchanges(k)%procid,dsolchanger(i)%sol(1:diexchanger(i)%muchineed(1),1:itest),&
diexchanger(i)%muchineed(1)*itest,mpi_double_precision,&
		dsolchanger(i)%procid,icpuid,mpi_comm_world,status,ierror)
		end if
		end do
end do
!$omp end master 
!$omp barrier

end subroutine exchange_lower
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exchange_higher(n)
!> @brief
!> this subroutine is exchanging the variables of all the halo cells for all the reconstruction stencils
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,itest,itee,iteedum,itemp1,itemp2,iavc,iavt,r0,r1,s0,s1,rb,sb
real,dimension(1:1)::dumts,rumts
integer:: n_requests,cell,baseS,baseR,v
integer, dimension(:), allocatable:: requests



	
itest=nof_variables+turbulenceequations+passivescalar


! ineedt=ineedhalo!direcexr(1)%tot
! tneedt=ineedhalos!direcexs(1)%tot

if (( turbulence .gt. 0).or.(passivescalar.gt.0)) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,halos_total
	solhis(i,1:nof_variables)=u_c_val(1,1:nof_variables,solhi_loc(i))
	solhis(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,solhi_loc(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

else


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,halos_total
	solhis(i,1:nof_variables)=u_c_val(1,1:nof_variables,solhi_loc(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do cell = 1, halos_total

  do v = 1, ilength1
    solhis_flat((cell-1)*ilength1 + v) = solhis(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif



n_requests = 0

allocate(requests(jtotal*2))
requests(:)=0
icpuid=n

do k=1,jtotal
        

        if ((jtot(k,1).eq.-1).and.(jtot(k,2).ne.-1))then
        
        
        
          n_requests = n_requests + 1
        iavc=jtot(k,2)
!        
        call mpi_isend(                                                     &
      dumts(1:1), & !sendbuf
      1, mpi_double_precision,       & !sendcount, sendtype
      jtot(k,3), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )


         n_requests = n_requests + 1
         iavc=jtot(k,2)
         r0=halo_offset(iavc)
         baseR = (r0-1)*itest + 1
         rb=halo_len(iavc)

#ifdef gpu
		!$omp target data use_device_ptr(solhir_flat)
#endif
   call mpi_irecv(                                                     &
      solhir_flat(baseR), rb*itest,    & !recvbuf
       mpi_double_precision,          & !recvcount, recvtype
     halo_proc(iavc), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif


        end if
        
       if ((jtot(k,1).ne.-1).and.(jtot(k,2).eq.-1))then
        
        n_requests = n_requests + 1
        iavt=jtot(k,1)
!         
		s0=halos_offset(iavt)
		baseS = (s0-1)*itest + 1
		sb=halos_len(iavt)

#ifdef gpu
		!$omp target data use_device_ptr(solhis_flat)
#endif
        call mpi_isend(                                                     &
      solhis_flat(baseS), sb*itest, & !sendbuf
      mpi_double_precision,       & !sendcount, sendtype
      halos_proc(iavt), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif





         n_requests = n_requests + 1
         iavt=jtot(k,1)
   call mpi_irecv(                                                     &
      dumts(1:1),    & !recvbuf
      1, mpi_double_precision,          & !recvcount, recvtype
     jtot(k,3), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
        
        
        
        end if
        
        
        
        if ((jtot(k,1).ne.-1).and.(jtot(k,2).ne.-1))then
        
       
       
        
        
         n_requests = n_requests + 1
		 iavt=jtot(k,1)
		 s0=halos_offset(iavt)
		baseS = (s0-1)*itest + 1
		sb=halos_len(iavt)

!          
#ifdef gpu
		!$omp target data use_device_ptr(solhis_flat)
#endif
        call mpi_isend(                                                     &
      solhis_flat(baseS), sb*itest, & !sendbuf
      mpi_double_precision,       & !sendcount, sendtype
      halos_proc(iavt), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif


         n_requests = n_requests + 1
         iavc=jtot(k,2)
         r0=halo_offset(iavc)
         baseR = (r0-1)*itest + 1
         rb=halo_len(iavc)




#ifdef gpu
		!$omp target data use_device_ptr(solhir_flat)
#endif
   call mpi_irecv(                                                     &
      solhir_flat(baseR), rb*itest,   & !recvbuf
       mpi_double_precision,          & !recvcount, recvtype
     halo_proc(iavc), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif
        
        
        
       end if
       
end do

       

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

deallocate(requests)





#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, halo_total
  do v = 1, ilength1
    solhir(cell, v) = solhir_flat((cell-1)*ilength1 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif






end subroutine exchange_higher




subroutine exchange_adda_diss(n)
!> @brief
!> this subroutine is exchanging the variables of all the halo cells for all the reconstruction stencils
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,itest,itee,iteedum,itemp1,itemp2,iavc,iavt,r0,r1,s0,s1,rb,sb
real,dimension(1:1)::dumts,rumts
integer:: n_requests,baseR,baseS,cell,v
integer, dimension(:), allocatable:: requests




itest=1












#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,halos_total
	solhisd(i)=ielem_diss(solhi_loc(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif





n_requests = 0

allocate(requests(jtotal*2))
requests(:)=0
icpuid=n

do k=1,jtotal


        if ((jtot(k,1).eq.-1).and.(jtot(k,2).ne.-1))then



          n_requests = n_requests + 1
        iavc=jtot(k,2)
!
        call mpi_isend(                                                     &
      dumts(1:1), & !sendbuf
      1, mpi_double_precision,       & !sendcount, sendtype
      jtot(k,3), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )


         n_requests = n_requests + 1
         iavc=jtot(k,2)

	r0=halo_offset(iavc)
         r1=halo_offset(iavc)+halo_len(iavc)-1
         rb=halo_len(iavc)

#ifdef gpu
		!$omp target data use_device_ptr(solhird)
#endif
   call mpi_irecv(                                                     &
      solhird(r0:r1),    & !recvbuf
      rb, mpi_double_precision,          & !recvcount, recvtype
     halo_proc(iavc), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif
        end if

       if ((jtot(k,1).ne.-1).and.(jtot(k,2).eq.-1))then

        n_requests = n_requests + 1
        iavt=jtot(k,1)
!
       s0=halos_offset(iavt)
         s1=halos_offset(iavt)+halos_len(iavt)-1
         sb=halos_len(iavt)

!
#ifdef gpu
		!$omp target data use_device_ptr(solhisd)
#endif
        call mpi_isend(                                                     &
      solhisd(s0:s1), sb, & !sendbuf
      mpi_double_precision,       & !sendcount, sendtype
      halos_proc(iavt), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif


         n_requests = n_requests + 1
         iavt=jtot(k,1)
   call mpi_irecv(                                                     &
      dumts(1:1),    & !recvbuf
      1, mpi_double_precision,          & !recvcount, recvtype
     jtot(k,3), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )



        end if



        if ((jtot(k,1).ne.-1).and.(jtot(k,2).ne.-1))then





         n_requests = n_requests + 1
		 iavt=jtot(k,1)
         s0=halos_offset(iavt)
         s1=halos_offset(iavt)+halos_len(iavt)-1
         sb=halos_len(iavt)

!
#ifdef gpu
		!$omp target data use_device_ptr(solhisd)
#endif
        call mpi_isend(                                                     &
      solhisd(s0:s1), sb, & !sendbuf
      mpi_double_precision,       & !sendcount, sendtype
      halos_proc(iavt), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

         n_requests = n_requests + 1
         iavc=jtot(k,2)
         r0=halo_offset(iavc)
         r1=halo_offset(iavc)+halo_len(iavc)-1
         rb=halo_len(iavc)

#ifdef gpu
		!$omp target data use_device_ptr(solhird)
#endif
   call mpi_irecv(                                                     &
      solhird(r0:r1),    & !recvbuf
      rb, mpi_double_precision,          & !recvcount, recvtype
     halo_proc(iavc), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif


       end if

end do



call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

deallocate(requests)

#ifdef gpu

#else
!$omp end master
!$omp barrier
#endif




end subroutine exchange_adda_diss


subroutine exchange_higher_av(n)
!> @brief
!> this subroutine is exchanging the averaged variables of all the halo cells for all the reconstruction stencils
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,itest,itee,iteedum,itemp1,itemp2,iavc,iavt,r0,r1,rb,s0,s1,sb
real,dimension(1:1)::dumts,rumts
integer:: n_requests,baseR,baseS,cell,v
integer, dimension(:), allocatable:: requests



itest=nof_variables+turbulenceequations+passivescalar



if (( turbulence .gt. 0).or.(passivescalar.gt.0)) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,halos_total
	solhis(i,1:nof_variables)=u_c_val(ind1,1:nof_variables,solhi_loc(i))
	solhis(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct_val(ind1,1:turbulenceequations+passivescalar,solhi_loc(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



else
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,halos_total
	solhis(i,1:nof_variables)=u_c_val(ind1,1:nof_variables,solhi_loc(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end if


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do cell = 1, halos_total
  do v = 1, ilength1
    solhis_flat((cell-1)*ilength1 + v) = solhis(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif




n_requests = 0

allocate(requests(jtotal*2))
requests(:)=0
icpuid=n

do k=1,jtotal


        if ((jtot(k,1).eq.-1).and.(jtot(k,2).ne.-1))then



          n_requests = n_requests + 1
        iavc=jtot(k,2)
!
        call mpi_isend(                                                     &
      dumts(1:1), & !sendbuf
      1, mpi_double_precision,       & !sendcount, sendtype
      jtot(k,3), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )


         n_requests = n_requests + 1
         iavc=jtot(k,2)
         r0=halo_offset(iavc)
         baseR = (r0-1)*itest + 1
         rb=halo_len(iavc)

#ifdef gpu
		!$omp target data use_device_ptr(solhir_flat)
#endif
   call mpi_irecv(                                                     &
      solhir_flat(baseR), rb*itest,    & !recvbuf
       mpi_double_precision,          & !recvcount, recvtype
     halo_proc(iavc), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif


        end if

       if ((jtot(k,1).ne.-1).and.(jtot(k,2).eq.-1))then

        n_requests = n_requests + 1
        iavt=jtot(k,1)
!
		s0=halos_offset(iavt)
		baseS = (s0-1)*itest + 1
		sb=halos_len(iavt)

!
#ifdef gpu
		!$omp target data use_device_ptr(solhis_flat)
#endif
        call mpi_isend(                                                     &
      solhis_flat(baseS), sb*itest, & !sendbuf
      mpi_double_precision,       & !sendcount, sendtype
      halos_proc(iavt), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif


         n_requests = n_requests + 1
         iavt=jtot(k,1)
   call mpi_irecv(                                                     &
      dumts(1:1),    & !recvbuf
      1, mpi_double_precision,          & !recvcount, recvtype
     jtot(k,3), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )



        end if



        if ((jtot(k,1).ne.-1).and.(jtot(k,2).ne.-1))then





         n_requests = n_requests + 1
		 iavt=jtot(k,1)
		 s0=halos_offset(iavt)
		baseS = (s0-1)*itest + 1
		sb=halos_len(iavt)

!
#ifdef gpu
		!$omp target data use_device_ptr(solhis_flat)
#endif
        call mpi_isend(                                                     &
      solhis_flat(baseS), sb*itest, & !sendbuf
      mpi_double_precision,       & !sendcount, sendtype
      halos_proc(iavt), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

         n_requests = n_requests + 1
         iavc=jtot(k,2)
         r0=halo_offset(iavc)
         baseR = (r0-1)*itest + 1
         rb=halo_len(iavc)



#ifdef gpu
		!$omp target data use_device_ptr(solhir_flat)
#endif

   call mpi_irecv(                                                     &
      solhir_flat(baseR), rb*itest,   & !recvbuf
       mpi_double_precision,          & !recvcount, recvtype
     halo_proc(iavc), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif



       end if

end do



call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

deallocate(requests)



#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, halo_total
  do v = 1, ilength1
    solhir(cell, v) = solhir_flat((cell-1)*ilength1 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	


end subroutine exchange_higher_av



subroutine exchange_higher_pre(n)
!> @brief
!> this subroutine is establishing the communication pattern and allocating the appropriate memory for 
!> exchanging the variables of all the halo cells for all the reconstruction stencils
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,itest,itee,iteedum,itemp1,itemp2,iavc,iavt
real,dimension(1:1)::dumts,rumts

	
itest=nof_variables+turbulenceequations+passivescalar

ineedt=direcexr(1)%tot
tneedt=direcexs(1)%tot

if (( turbulence .gt. 0).or.(passivescalar.gt.0)) then
!$omp do
do i=1,tneedt
	do k=1,direcexs(i)%muchtheyneed(1)
		  diexsolhis(i)%sol(k,1:nof_variables)=u_c_val(1,1:nof_variables,direcexs(i)%localref(k))
		  diexsolhis(i)%sol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,direcexs(i)%localref(k))
		  
	end do
end do
!$omp end do

else
!  dsolchanges(i)%sol(k,1:nof_variables)=u_c_val(1,1:nof_variables,diexchanges(i)%localref(k))
!$omp do
do i=1,tneedt
	do k=1,direcexs(i)%muchtheyneed(1)
		  diexsolhis(i)%sol(k,1:nof_variables+turbulenceequations+passivescalar)=u_c_val(1,1:nof_variables+turbulenceequations+passivescalar,direcexs(i)%localref(k))
	end do
end do
!$omp end do
end if

!$omp barrier
!$omp master

jtotal=0;jtotal1=0;jtotal2=0;jtotal3=0
icpuid=n
 dumts=zero
		do i=0,isize-1
			if (i.ne.n) then
				do j=1,tneedt
					iavt=10000
					if (direcexs(j)%procid.eq.i)then
					itemp1=direcexs(j)%muchtheyneed(1)*itest
					iavt=j 
					go to 7001
                                   
					end if
				end do
				7001 continue
				do k=1,ineedt
					iavc=10000
					if (direcexr(k)%procid.eq.i) then
					iavc=k
					itemp2=direcexr(k)%muchineed(1)*itest
					go to 8001
					
					end if
				end do
				8001 continue
				if ((iavt.eq.10000).and.(iavc.ne.10000)) then
				!message 1
! 				
				jtotal1=jtotal1+1
				call mpi_sendrecv(dumts(1:1),1,mpi_double_precision,i,icpuid,&
				diexsolhir(iavc)%sol(1:direcexr(iavc)%muchineed(1),1:itest),&
			      itemp2,mpi_double_precision,diexsolhir(iavc)%procid,diexsolhir(iavc)%procid,mpi_comm_world,status,ierror)
				end if
				if ((iavt.ne.10000).and.(iavc.eq.10000)) then
				!message 1
! 			
				jtotal2=jtotal2+1
				call mpi_sendrecv(diexsolhis(iavt)%sol(1:direcexs(iavt)%muchtheyneed(1),1:itest),itemp1,&
			      mpi_double_precision,diexsolhis(iavt)%procid,icpuid,&
				dumts(1:1),1,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
				end if
				if ((iavt.ne.10000).and.(iavc.ne.10000)) then
				!message 1
				jtotal3=jtotal3+1
! 				
				
				call mpi_sendrecv(diexsolhis(iavt)%sol(1:direcexs(iavt)%muchtheyneed(1),1:itest),itemp1,&
			      mpi_double_precision,diexsolhis(iavt)%procid,icpuid,&
				diexsolhir(iavc)%sol(1:direcexr(iavc)%muchineed(1),1:itest),itemp2,mpi_double_precision,&
			      diexsolhir(iavc)%procid,diexsolhir(iavc)%procid,mpi_comm_world,status,ierror)
				end if		
			end if	
		end do
		
		
		allocate(jtot1(jtotal1,3),jtot2(jtotal2,3),jtot3(jtotal3,3),jtot(jtotal1+jtotal2+jtotal3,3))
		jtot1(:,:)=-1;jtot2(:,:)=-1;jtot3(:,:)=-1;jtot(:,:)=-1
!  		
		jtotal=0;jtotal1=0;jtotal2=0;jtotal3=0
		
icpuid=n
 dumts=zero
		do i=0,isize-1
			if (i.ne.n) then
				do j=1,tneedt
					iavt=10000
					if (direcexs(j)%procid.eq.i)then
					itemp1=direcexs(j)%muchtheyneed(1)*itest
					iavt=j 
					go to 10001
                                   
					end if
				end do
				10001 continue
				do k=1,ineedt
					iavc=10000
					if (direcexr(k)%procid.eq.i) then
					iavc=k
					itemp2=direcexr(k)%muchineed(1)*itest
					go to 11001
					
					end if
				end do
				11001 continue
				if ((iavt.eq.10000).and.(iavc.ne.10000)) then
				!message 1
! 				jtotal1=jtotal1+1
! 				jtot1(jtotal1,1)=-1
! 				jtot1(jtotal1,2)=j
				
				jtotal=jtotal+1
				jtot(jtotal,1)=-1
				jtot(jtotal,2)=k
				jtot(jtotal,3)=i
				call mpi_sendrecv(dumts(1:1),1,mpi_double_precision,i,icpuid,&
				diexsolhir(iavc)%sol(1:direcexr(iavc)%muchineed(1),1:itest),&
			      itemp2,mpi_double_precision,diexsolhir(iavc)%procid,diexsolhir(iavc)%procid,mpi_comm_world,status,ierror)
				end if
				if ((iavt.ne.10000).and.(iavc.eq.10000)) then
				!message 1
! 				jtotal2=jtotal2+1
! 				jtot1(jtotal2,1)=k
! 				jtot1(jtotal2,2)=-1
				jtotal=jtotal+1
				jtot(jtotal,1)=j
				jtot(jtotal,2)=-1
				jtot(jtotal,3)=i
				call mpi_sendrecv(diexsolhis(iavt)%sol(1:direcexs(iavt)%muchtheyneed(1),1:itest),itemp1,&
			      mpi_double_precision,diexsolhis(iavt)%procid,icpuid,&
				dumts(1:1),1,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
				end if
				if ((iavt.ne.10000).and.(iavc.ne.10000)) then
				!message 1
! 				jtotal3=jtotal3+1
! 				jtot3(jtotal3,1)=j
! 				jtot3(jtotal3,2)=k
				jtotal=jtotal+1
				jtot(jtotal,1)=j
				jtot(jtotal,2)=k
				call mpi_sendrecv(diexsolhis(iavt)%sol(1:direcexs(iavt)%muchtheyneed(1),1:itest),itemp1,&
			      mpi_double_precision,diexsolhis(iavt)%procid,icpuid,&
				diexsolhir(iavc)%sol(1:direcexr(iavc)%muchineed(1),1:itest),itemp2,mpi_double_precision,&
			      diexsolhir(iavc)%procid,diexsolhir(iavc)%procid,mpi_comm_world,status,ierror)
				end if		
			end if	
		end do
		
		
!$omp end master 
!$omp barrier
	
!                                

end subroutine exchange_higher_pre









subroutine exhboundhigher(n)
!> @brief
!> this subroutine is communicating the boundary extrapolated values for the variables and their gradients
!> for the gaussian quadrature points of direct-side neighbours between mpi processes 
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,ittt,iex,imulti,k_cnt,nvar,r0,r1,rb,s0,s1,sb
integer::itee,iteedum,jk,jjk,jjk4,jjk12,imulti2,icpe,jmnb,j76,j78,j79,j80,imulti3,i_cnt,cinout2
integer:: n_requests,baseR,baseS,cell,v
integer, dimension(:), allocatable:: requests
real::pr_t31,pr_t32,pr_t33,pr_t34,pr_t35,temp_prin,temp_prout
 cinout2=0

pr_t31=zero
pr_t32=zero
pr_t33=zero
pr_t34=zero
pr_t35=zero
temp_prin=zero
temp_prout=zero

indl=ineedbound
tndl=ineedbounds

if (dimensiona.eq.3)then


if( itestcase.eq.4)then
i_cnt=(nof_variables+turbulenceequations+passivescalar)+(((nof_variables-1)+turbulenceequations+passivescalar)*3)
else

i_cnt=nof_variables
end if
else
if( itestcase.eq.4)then
i_cnt=(nof_variables+turbulenceequations+passivescalar)+(((nof_variables-1)+turbulenceequations+passivescalar)*2)
else

i_cnt=nof_variables
end if
end if


if (itestcase.le.3) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
    do i=1,bounds_total
			boundhis(i,1:nof_variables)=rec_uleft(1:nof_variables,need_side(i),need_q(i),need_loc(i))
    end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end  if


if (itestcase.eq.4) then

if (turbulence.ne.1)then
#ifdef gpu
!$omp target teams distribute parallel do private(ittt)
#else
!$omp do
#endif
do i=1,bounds_total
			boundhis(i,1:nof_variables)=rec_uleft(1:nof_variables,need_side(i),need_q(i),need_loc(i))
			ittt=0
			do iex=1,nof_variables-1
			do nvar=1,dims
			ittt=ittt+1
			boundhis(i,nof_variables+ittt)=rec_uleftv(nvar, iex, &
             need_side(i), need_q(i), need_loc(i))
			end do
			end do

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

else
#ifdef gpu
!$omp target teams distribute parallel do private(ittt)
#else
!$omp do
#endif
do i=1,bounds_total
			boundhis(i,1:nof_variables)=rec_uleft(1:nof_variables,need_side(i),need_q(i),need_loc(i))
			boundhis(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,need_side(i),need_q(i),need_loc(i))
			ittt=0
			do iex=1,nof_variables-1
			do nvar=1,dims
			ittt=ittt+1
			boundhis(i,nof_variables+turbulenceequations+passivescalar+ittt)=rec_uleftv(nvar, iex, &
             need_side(i), need_q(i), need_loc(i))
			end do
			end do
			do iex=1,turbulenceequations+passivescalar
			do nvar=1,dims
			ittt=ittt+1
			boundhis(i,nof_variables+turbulenceequations+passivescalar+ittt)=rec_uleftturbv(nvar, iex, &
             need_side(i), need_q(i), need_loc(i))
			end do
			end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if
end  if

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do cell = 1, bounds_total
  do v = 1, ilength2
    boundhis_flat((cell-1)*ilength2 + v) = boundhis(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif








n_requests = 0
allocate(requests(2*indl))
requests(:)=0
icpuid=n



do k=1,indl

   ! search unique j such
   j = 1
   do while(bound_proc(k) .ne. bounds_proc(j))
      j = j + 1
   end do

   ! non-blocking send
   n_requests = n_requests + 1

   s0 = bounds_offset(j)
	baseS=(s0-1)*i_cnt + 1
	sb= bounds_len(j)


#ifdef gpu
		!$omp target data use_device_ptr(boundhis_flat)
#endif

   call mpi_isend(                                                     &
      boundhis_flat(baseS), & !sendbuf
      sb*i_cnt, mpi_double_precision,       & !sendcount, sendtype
      bounds_proc(j), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

   ! non-blocking receive
   n_requests = n_requests + 1

    r0 = bound_offset(k)
    baseR =(r0-1)*i_cnt + 1
	rb= bound_len(k)
#ifdef gpu
		!$omp target data use_device_ptr(boundhir_flat)
#endif
   call mpi_irecv(                                                     &
      boundhir_flat(baseR),    & !recvbuf
      rb*i_cnt, mpi_double_precision,          & !recvcount, recvtype
      bound_proc(k), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

end do

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

deallocate(requests)


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, bound_total
  do v = 1, ilength2
    boundhir(cell, v) = boundhir_flat((cell-1)*ilength2 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end subroutine exhboundhigher




subroutine exhboundhigher2(n)
!> @brief
!> this subroutine is communicating the boundary extrapolated values for the variables and their gradients
!> for the gaussian quadrature points of direct-side neighbours between mpi processes for the implicit time stepping

implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,ittt,iex,imulti,k_cnt,nvar,r0,r1,rb,s0,s1,sb
integer::itee,iteedum,jk,jjk,jjk4,jjk12,imulti2,icpe,jmnb,j76,j78,j79,j80,imulti3,n_requests,cell,baseS,baseR,v
integer, dimension(:), allocatable:: requests
indl=ineedbound
tndl=ineedbounds

if (itestcase.lt.3)then
	iex=1
      imulti=iex
      imulti2=iex
end if
if (itestcase.eq.3)then
	iex=nof_variables
    imulti2=iex
end if
if (itestcase.eq.4)then
    k_cnt=(nof_variables+turbulenceequations+passivescalar)
      iex = (nof_variables+turbulenceequations+passivescalar)
       imulti2=k_cnt
       imulti3=k_cnt
end if







imulti=iex


if (itestcase.le.3) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,bounds_total
	if (relax.eq.3)then
		if (iscoun.eq.1)then
		boundhisi(i,1:iex)=-rhs_val(1:iex,need_loc(i))/impdiag_mf(need_loc(i))
		else
		boundhisi(i,1:iex)=-(rhs_val(1:nof_variables,need_loc(i))+((((1.5*u_c_val(1,1:nof_variables,need_loc(i)))-(2.0d0*u_c_val(2,1:nof_variables,need_loc(i)))+(0.5d0*u_c_val(3,1:nof_variables,need_loc(i))))/(dt))*ielem_totvolume(need_loc(i))))/impdiag_mf(need_loc(i))
		end if
	else
		boundhisi(i,1:iex)=impdu(need_loc(i),1:iex)
	end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end  if



if (itestcase.eq.4) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,bounds_total
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		if (relax.eq.3)then
			if (iscoun.eq.1)then
			boundhisi(i,1:nof_variables)=-rhs_val(1:nof_variables,need_loc(i))/impdiag_mf(need_loc(i))
			boundhisi(i,nof_variables:nof_variables+turbulenceequations+passivescalar)=-rhst_val(:,need_loc(i))/impdiagt(:,need_loc(i))
			else
			boundhisi(i,1:nof_variables)=-(rhs_val(1:nof_variables,need_loc(i))+((((1.5*u_c_val(1,1:nof_variables,need_loc(i)))-(2.0d0*u_c_val(2,1:nof_variables,need_loc(i)))+(0.5d0*u_c_val(3,1:nof_variables,need_loc(i))))/(dt))*ielem_totvolume(need_loc(i))))/impdiag_mf(need_loc(i))
			boundhisi(i,nof_variables:nof_variables+turbulenceequations+passivescalar)=-(rhst_val(:,need_loc(i))+((((1.5*u_ct_val(1,:,need_loc(i)))-(2.0d0*u_ct_val(2,:,need_loc(i)))+(0.5d0*u_ct_val(3,:,need_loc(i))))/(dt))*ielem_totvolume(need_loc(i))))/impdiagt(:,need_loc(i))

			end if
		else
			boundhisi(i,:)=impdu(need_loc(i),:)
		end if
	end if
	if ((turbulence.eq.0).or.(passivescalar.eq.0))then
			if (relax.eq.3)then
			if (iscoun.eq.1)then
			boundhisi(i,1:iex)=-rhs_val(1:iex,need_loc(i))/impdiag_mf(need_loc(i))

			else
			boundhisi(i,1:iex)=-(rhs_val(1:nof_variables,need_loc(i))+((((1.5*u_c_val(1,1:nof_variables,need_loc(i)))-(2.0d0*u_c_val(2,1:nof_variables,need_loc(i)))+(0.5d0*u_c_val(3,1:nof_variables,need_loc(i))))/(dt))*ielem_totvolume(need_loc(i))))/impdiag_mf(need_loc(i))
			end if
		else
			boundhisi(i,1:iex)=impdu(need_loc(i),1:iex)
		end if
	end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if



#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do cell = 1, bounds_total
  do v = 1, ilength1
    boundhisi_flat((cell-1)*ilength1 + v) = boundhisi(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif



n_requests = 0
 allocate(requests(2*indl))
requests(:)=0
icpuid=n






do k=1,indl

   j = 1
   do while(bound_proc(k) .ne. bounds_proc(j))
      j = j + 1
   end do

   ! non-blocking send
   n_requests = n_requests + 1

    s0 = bounds_offset(j)
	baseS=(s0-1)*imulti2 + 1
	sb= bounds_len(j)
#ifdef gpu
		!$omp target data use_device_ptr(boundhisi_flat)
#endif
   call mpi_isend(                                                     &
      boundhisi_flat(baseS), & !sendbuf
      sb*imulti2, mpi_double_precision,       & !sendcount, sendtype
      bounds_proc(j), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

   ! non-blocking receive
   n_requests = n_requests + 1
    r0 = bound_offset(k)
    baseR =(r0-1)*imulti2 + 1
	rb= bound_len(k)
#ifdef gpu
		!$omp target data use_device_ptr(boundhiri_flat)
#endif
   call mpi_irecv(                                                     &
      boundhiri_flat(baseR),    & !recvbuf
      rb*imulti2, mpi_double_precision,          & !recvcount, recvtype
      bound_proc(k), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

end do

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

 deallocate(requests)


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, bound_total
  do v = 1, ilength1
    boundhiri(cell, v) = boundhiri_flat((cell-1)*ilength1 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine exhboundhigher2



subroutine exhboundhigherlu(n)
!> @brief
!> this subroutine is communicating the boundary extrapolated values for the variables and their gradients
!> for the gaussian quadrature points of direct-side neighbours between mpi processes for the implicit time stepping

implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,ittt,iex,imulti,k_cnt,nvar,n_requests,r0,r1,rb,s0,s1,sb
integer::itee,iteedum,jk,jjk,jjk4,jjk12,imulti2,icpe,jmnb,j76,j78,j79,j80,imulti3
integer::cell,baseR,baseS,v
integer, dimension(:), allocatable:: requests
indl=ineedbound
tndl=ineedbounds

if (itestcase.lt.3)then
	iex=1
      imulti=iex
      imulti2=iex
end if
if (itestcase.eq.3)then
	iex=nof_variables
    imulti2=iex
end if
if (itestcase.eq.4)then
    k_cnt=(nof_variables+turbulenceequations+passivescalar)
      iex = (nof_variables+turbulenceequations+passivescalar)
       imulti2=k_cnt
       imulti3=k_cnt
end if







imulti=iex


if (itestcase.le.3) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,bounds_total
boundhisi(i,1:iex)=impdu(need_loc(i),1:ilength1)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end  if

if (itestcase.eq.4) then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,bounds_total
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	  boundhisi(i,1:nof_variables)=impdu(need_loc(i),1:nof_variables)
	
	do nvar=1,0+turbulenceequations+passivescalar
	boundhisi(i,nof_variables+nvar)=impdu(need_loc(i),nof_variables+nvar)
	end do
      end if
!     
    if ((turbulence .eq. 0).and.(passivescalar.eq.0)) then
	   boundhisi(i,1:nof_variables)=impdu(need_loc(i),1:nof_variables)
    end if ! turbulence

end do 
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do cell = 1, bounds_total
  do v = 1, ilength1
    boundhisi_flat((cell-1)*ilength1 + v) = boundhisi(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif





n_requests = 0
allocate(requests(2*indl))
requests(:)=0
icpuid=n


do k=1,indl

   j = 1
   do while(bound_proc(k) .ne. bounds_proc(j))
      j = j + 1
   end do

   ! non-blocking send
   n_requests = n_requests + 1

    s0 = bounds_offset(j)
	baseS=(s0-1)*imulti2 + 1
	sb= bounds_len(j)
#ifdef gpu
		!$omp target data use_device_ptr(boundhisi_flat)
#endif
   call mpi_isend(                                                     &
      boundhisi_flat(baseS), & !sendbuf
      sb*imulti2, mpi_double_precision,       & !sendcount, sendtype
      bounds_proc(j), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

   ! non-blocking receive
   n_requests = n_requests + 1
    r0 = bound_offset(k)
    baseR =(r0-1)*imulti2 + 1
	rb= bound_len(k)
#ifdef gpu
		!$omp target data use_device_ptr(boundhiri_flat)
#endif
   call mpi_irecv(                                                     &
      boundhiri_flat(baseR),    & !recvbuf
      rb*imulti2, mpi_double_precision,          & !recvcount, recvtype
      bound_proc(k), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

end do

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
!!$omp end target data
 deallocate(requests)

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, bound_total
  do v = 1, ilength1
    boundhiri(cell, v) = boundhiri_flat((cell-1)*ilength1 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif







end subroutine exhboundhigherlu

subroutine exhboundhigher_mood(n)
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,ittt,iex,imulti,k_cnt,nvar,r0,r1,rb,s0,s1,sb
integer::itee,iteedum,jk,jjk,jjk4,jjk12,imulti2,icpe,jmnb,j76,j78,j79,j80,imulti3,i_cnt,cinout2
integer:: n_requests,cell,baseR,baseS,v
integer, dimension(:), allocatable:: requests
real::pr_t31,pr_t32,pr_t33,pr_t34,pr_t35,temp_prin,temp_prout
 cinout2=0
indl=ineedbound
tndl=ineedbounds

pr_t31=zero
pr_t32=zero
pr_t33=zero
pr_t34=zero
pr_t35=zero
temp_prin=zero
temp_prout=zero




if (dimensiona.eq.3)then


            if( itestcase.eq.4)then
            i_cnt=1
            else

            i_cnt=1
            end if
else
        if( itestcase.eq.4)then
        i_cnt=1
        else

        i_cnt=1
        end if
end if



#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,bounds_total
			boundhism(i)=ielem_mood(need_loc(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif





n_requests = 0
allocate(requests(2*indl))

icpuid=n



do k=1,indl

   j = 1
   do while(bound_proc(k) .ne. bounds_proc(j))
      j = j + 1
   end do

   ! non-blocking send
   n_requests = n_requests + 1
   s0 = bounds_offset(j)
	s1 = bounds_offset(j) + bounds_len(j) - 1
	sb= bounds_len(j)
#ifdef gpu
		!$omp target data use_device_ptr(boundhism)
#endif
   call mpi_isend(                                                     &
      boundhism(s0:s1), & !sendbuf
      sb, mpi_double_precision,       & !sendcount, sendtype
      bounds_proc(j), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

   ! non-blocking receive
   n_requests = n_requests + 1
   r0 = bound_offset(k)
   r1 = bound_offset(k) + bound_len(k) - 1
   rb= bound_len(k)
#ifdef gpu
		!$omp target data use_device_ptr(boundhism)
#endif
   call mpi_irecv(                                                     &
      boundhirm(r0:r1),    & !recvbuf
      rb, mpi_double_precision,          & !recvcount, recvtype
      bound_proc(k), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

end do

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)



deallocate(requests)

#ifdef gpu

#else
!$omp end master
!$omp barrier
#endif

end subroutine exhboundhigher_mood





subroutine exhboundhigher_dg(n)
!> @brief
!> this subroutine is communicating the boundary extrapolated values for the variables and their gradients
!> for the gaussian quadrature points of direct-side neighbours between mpi processes
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,ittt,iex,imulti,k_cnt,nvar,r0,r1,rb,s0,s1,sb
integer::itee,iteedum,jk,jjk,jjk4,jjk12,imulti2,icpe,jmnb,j76,j78,j79,j80,imulti3,i_cnt,cinout2
integer:: n_requests,cell,baseR,baseS,v
integer, dimension(:), allocatable:: requests
real::pr_t31,pr_t32,pr_t33,pr_t34,pr_t35,temp_prin,temp_prout
 cinout2=0
indl=ineedbound
tndl=ineedbounds

pr_t31=zero
pr_t32=zero
pr_t33=zero
pr_t34=zero
pr_t35=zero
temp_prin=zero
temp_prout=zero




if( itestcase.eq.4)then
i_cnt=(nof_variables+turbulenceequations+passivescalar)+(((nof_variables-1)+turbulenceequations+passivescalar)*dimensiona)
else

i_cnt=nof_variables
end if



if (itestcase.le.3)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
    do i=1,bounds_total
			boundhis_dg(i,1:nof_variables)=rec_uleft_dg(1:nof_variables,need_side(i),need_q(i),need_loc(i))
    end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
if (itestcase.eq.4) then
#ifdef gpu
!$omp target teams distribute parallel do private(ittt)
#else
!$omp do
#endif
		do i=1,bounds_total
				boundhis_dg(i,1:nof_variables)=rec_uleft_dg(1:nof_variables,need_side(i),need_q(i),need_loc(i))
				ittt=0
				do iex=1,nof_variables-1
					do nvar=1,dims
					ittt=ittt+1
					boundhis_dg(i,nof_variables+ittt)=rec_uleftv(nvar, iex, &
             need_side(i), need_q(i), need_loc(i))
					end do
				end do
		end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end  if



#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do cell = 1, bounds_total
  do v = 1, ilength1
    boundhis_dgflat((cell-1)*ilength1 + v) = boundhis_dg(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif




!call mpi_barrier(mpi_comm_world,ierror)

n_requests = 0
allocate(requests(2*indl))
requests(:)=0
icpuid=n




do k=1,indl


   j = 1
   do while(bound_proc(k) .ne. bounds_proc(j))
      j = j + 1
   end do

   ! non-blocking send
   n_requests = n_requests + 1
  s0 = bounds_offset(j)
	baseS=(s0-1)*i_cnt + 1
	sb= bounds_len(j)


#ifdef gpu
		!$omp target data use_device_ptr(boundhis_dgflat)
#endif
   call mpi_isend(                                                     &
      boundhis_dgflat(baseS), & !sendbuf
      sb*i_cnt, mpi_double_precision,       & !sendcount, sendtype
      bounds_proc(j), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif
   ! non-blocking receive
   n_requests = n_requests + 1
   r0 = bound_offset(k)
    baseR =(r0-1)*i_cnt + 1
	rb= bound_len(k)
#ifdef gpu
		!$omp target data use_device_ptr(boundhir_dgflat)
#endif
   call mpi_irecv(                                                     &
      boundhir_dgflat(baseR),    & !recvbuf
      rb*i_cnt, mpi_double_precision,          & !recvcount, recvtype
      bound_proc(k), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

end do

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

deallocate(requests)

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, bound_total
  do v = 1, ilength2
    boundhir_dg(cell, v) = boundhir_dgflat((cell-1)*ilength2 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine exhboundhigher_dg





subroutine exhboundhigher_dg2(n)
!> @brief
!> this subroutine is communicating the boundary extrapolated values for the variables and their gradients
!> for the gaussian quadrature points of direct-side neighbours between mpi processes
implicit none
integer,intent(in)::n
integer::i,j,k,l,m,o,p,q,ineedt,tneedt,indl,tndl,icpuid,ittt,iex,imulti,k_cnt,nvar
integer::itee,iteedum,jk,jjk,jjk4,jjk12,imulti2,icpe,jmnb,j76,j78,j79,j80,imulti3,i_cnt,cinout2,r0,r1,rb,s0,s1,sb
integer:: n_requests,cell,baseR,baseS,v
integer, dimension(:), allocatable:: requests
real::pr_t31,pr_t32,pr_t33,pr_t34,pr_t35,temp_prin,temp_prout
 cinout2=0
indl=ineedbound
tndl=ineedbounds

pr_t31=zero
pr_t32=zero
pr_t33=zero
pr_t34=zero
pr_t35=zero
temp_prin=zero
temp_prout=zero





if (dimensiona.eq.3)then


if( itestcase.eq.4)then
i_cnt=(nof_variables+turbulenceequations+passivescalar)+(((nof_variables-1)+turbulenceequations+passivescalar)*dimensiona)
else

i_cnt=nof_variables
end if

else
if( itestcase.eq.4)then
i_cnt=(nof_variables+turbulenceequations+passivescalar)+(((nof_variables-1)+turbulenceequations+passivescalar)*dimensiona)
else
	i_cnt=nof_variables
end if
end if


if (itestcase.le.3) then
	if((multispecies.eq.1).and.(br2_yn.eq.1)) then
#ifdef gpu
!$omp target teams distribute parallel do private (ittt)
#else
!$omp do
#endif
		do i=1,bounds_total
				boundhis_dg(i,1:nof_variables)=rec_uleft_dg(1:nof_variables,need_side(i),need_q(i),need_loc(i))
				ittt=0
				do iex=1,nof_variables-1
					do nvar=1,dims
					ittt=ittt+1
					boundhis_dg(i,nof_variables+ittt)=rec_br2_aux_var(iex,nvar,need_side(i),need_q(i),need_loc(i))
					end do
				end do
		end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

	else

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
		do i=1,bounds_total
				boundhis_dg(i,1:nof_variables)=rec_uleft_dg(1:nof_variables,need_side(i),need_q(i),need_loc(i))
		end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
end  if


if (itestcase.eq.4) then
#ifdef gpu
!$omp target teams distribute parallel do private (ittt)
#else
!$omp do
#endif
		do i=1,bounds_total
				boundhis_dg(i,1:nof_variables)=rec_uleft_dg(1:nof_variables,need_side(i),need_q(i),need_loc(i))
				ittt=0
				do iex=1,nof_variables-1
					do nvar=1,dims
					ittt=ittt+1
					boundhis_dg(i,nof_variables+ittt)=rec_uleftv(nvar, iex, &
             need_side(i), need_q(i), need_loc(i))
					end do
				end do
		end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end  if







#ifdef gpu
!$omp target teams distribute parallel do private
#else
!$omp do
#endif
do cell = 1, bounds_total
  do v = 1, ilength2
    boundhis_dgflat((cell-1)*ilength2 + v) = boundhis_dg(cell, v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
!$omp barrier
!$omp master
#endif






n_requests = 0
allocate(requests(2*indl))
requests(:)=0
icpuid=n







do k=1,indl


   j = 1
   do while(bound_proc(k) .ne. bounds_proc(j))
      j = j + 1
   end do

   ! non-blocking send
   n_requests = n_requests + 1
  s0 = bounds_offset(j)
	baseS=(s0-1)*i_cnt + 1
	sb= bounds_len(j)



#ifdef gpu
		!$omp target data use_device_ptr(boundhis_dgflat)
#endif
   call mpi_isend(                                                     &
      boundhis_dgflat(baseS), & !sendbuf
      sb*i_cnt, mpi_double_precision,       & !sendcount, sendtype
      bounds_proc(j), 0,                                        & !destination, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif
   ! non-blocking receive
   n_requests = n_requests + 1
   r0 = bound_offset(k)
    baseR =(r0-1)*i_cnt + 1
	rb= bound_len(k)
#ifdef gpu
		!$omp target data use_device_ptr(boundhir_flat)
#endif
   call mpi_irecv(                                                     &
      boundhir_flat(baseR),    & !recvbuf
      rb*i_cnt, mpi_double_precision,          & !recvcount, recvtype
      bound_proc(k), 0,                                        & !source, tag
      mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
   )
#ifdef gpu
		!$omp end target data
#endif

end do

call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)

deallocate(requests)

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp end master
!$omp barrier
!$omp do
#endif
do cell = 1, bound_total
  do v = 1, ilength2
    boundhir_dg(cell, v) = boundhir_dgflat((cell-1)*ilength2 + v)
  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end subroutine exhboundhigher_dg2


subroutine flat_comm(n)
implicit none
integer,intent(in)::n
integer::i,j,k,row



ineedhalo=direcexr(1)%tot
ineedbound=diexchanger(1)%tot

ineedhalos=direcexs(1)%tot
ineedbounds=diexchanges(1)%tot



if (itestcase.eq.4)then
ilength1=(nof_variables+turbulenceequations+passivescalar)
ilength2=(nof_variables+turbulenceequations+passivescalar)+(((nof_variables-1)+turbulenceequations+passivescalar)*dimensiona)
else
ilength1=nof_variables
ilength2=nof_variables
end if







!first allocate memory for high-order solution exchanges
allocate(halo_len(1:ineedhalo),halo_offset(1:ineedhalo),halo_proc(1:ineedhalo))
allocate(halos_len(1:ineedhalos),halos_offset(1:ineedhalos),halos_proc(1:ineedhalos))
do i=1,ineedhalo
      halo_len(i)=direcexr(i)%muchineed(1)  !store the number of elements needed per cpu
      halo_proc(i)=direcexr(i)%procid
end do
do i=1,ineedhalos
      halos_len(i)=direcexs(i)%muchtheyneed(1)  !store the number of elements needed per cpu
	  halos_proc(i)=direcexs(i)%procid
end do




halo_offset(1)=1                             !start computing the entries for each one of them

do j=2,ineedhalo
   halo_offset(j)=halo_offset(j-1)+halo_len(j-1)
!    write(140+n,*)"receive offset",halo_offset(j)
end do


halos_offset(1)=1                             !start computing the entries for each one of them

do j=2,ineedhalos
   halos_offset(j)=halos_offset(j-1)+halos_len(j-1)
!    write(140+n,*)"send offset",halos_offset(j)
end do




halo_total = halo_offset(ineedhalo) + halo_len(ineedhalo) - 1
halos_total = halos_offset(ineedhalos) + halos_len(ineedhalos) - 1

allocate(solhir(halo_total,1:ilength1));solhir=zero
allocate(solhis(halos_total,1:ilength1));solhis=zero

allocate(solhir_flat(1:halo_total*ilength1));solhir_flat=zero
allocate(solhis_flat(1:halos_total*ilength1));solhis_flat=zero





allocate(solhi_loc(halos_total))

do i=1,ineedhalos
	do k=1,halos_len(i)
		row=halos_offset(i)+k-1
		solhi_loc(row)=direcexs(i)%localref(k)
	end do
end do







if (adda.eq.1)then
allocate(solhird(halo_total));solhird=zero
allocate(solhisd(halos_total));solhisd=zero
end if




allocate(bound_len(ineedbound), bound_offset(ineedbound),bound_proc(ineedbound))


allocate(bounds_len(ineedbounds), bounds_offset(ineedbounds),bounds_proc(ineedbounds))

do i = 1, ineedbound
  bound_len(i) = diexchanger(i)%muchineed(1)
  bound_proc(i) = diexchanger(i)%procid
end do


do i = 1, ineedbounds
  bounds_len(i) = diexchanges(i)%muchtheyneed(1)
  bounds_proc(i) = diexchanges(i)%procid
end do


  bound_offset(1) = 1
  do j = 2, ineedbound
    bound_offset(j) = bound_offset(j-1) + bound_len(j-1)
  end do

  bound_total = bound_offset(ineedbound) + bound_len(ineedbound) - 1


   bounds_offset(1) = 1
  do j = 2, ineedbounds
    bounds_offset(j) = bounds_offset(j-1) + bounds_len(j-1)
  end do

  bound_total = bound_offset(ineedbound) + bound_len(ineedbound) - 1
  bounds_total = bounds_offset(ineedbounds) + bounds_len(ineedbounds) - 1


  allocate(need_side(bounds_total),need_q(bounds_total),need_loc(bound_total))

  ! Flat per-neighbor k-lists into a single row-based list
  do i = 1, ineedbounds
  do k = 1, bounds_len(i)
    row = bounds_offset(i) + k - 1
    need_side(row) = diexchanges(i)%sidetheyneed(k)
    need_q(row)    = diexchanges(i)%qtheyneed(k)
    need_loc(row)  = diexchanges(i)%localref(k)
  end do
end do



  allocate(boundhir(bound_total,1:ilength2));boundhir=zero
  allocate(boundhis(bounds_total,1:ilength2));boundhis=zero
  allocate(boundhir_flat(1:bound_total*ilength2));boundhir_flat=zero
  allocate(boundhis_flat(1:bounds_total*ilength2));boundhis_flat=zero
  if (dg.eq.1)then
  allocate(boundhir_dg(bound_total,1:ilength2)); boundhir_dg=zero
  allocate(boundhis_dg(bounds_total,1:ilength2)); boundhis_dg=zero
  allocate(boundhir_dgflat(1:bound_total*ilength2));boundhir_dgflat=zero
  allocate(boundhis_dgflat(1:bounds_total*ilength2));boundhis_dgflat=zero
  end if

  if (mood.eq.1)then
  allocate(boundhirm(bound_total)); boundhirm=zero
   allocate(boundhism(bounds_total));boundhism=zero
  end if

  if ((rungekutta.ge.10))then
  allocate(boundhiri(bound_total,1:ilength1)); boundhiri=zero
  allocate(boundhisi(bounds_total,1:ilength1)); boundhisi=zero
  allocate(boundhiri_flat(1:bound_total*ilength1)); boundhiri_flat=zero
  allocate(boundhisi_flat(1:bounds_total*ilength1)); boundhisi_flat=zero
  end if




!these are for the mpi net




if (allocated(diexboundhir))  deallocate(diexboundhir)
if (allocated(diexboundhirr)) deallocate(diexboundhirr)
if (allocated(diexboundhiri)) deallocate(diexboundhiri)
if (allocated(diexboundhis))  deallocate(diexboundhis)
if (allocated(diexboundhiss)) deallocate(diexboundhiss)
if (allocated(diexboundhisi)) deallocate(diexboundhisi)
if (allocated(diexsolhir))    deallocate(diexsolhir)
if (allocated(diexsolhis))    deallocate(diexsolhis)
if (allocated(diexsolhird))   deallocate(diexsolhird)
if (allocated(diexsolhisd))   deallocate(diexsolhisd)
if (allocated(diexcordr))     deallocate(diexcordr)
if (allocated(diexcords))     deallocate(diexcords)
if (allocated(diexchanger))   deallocate(diexchanger)
if (allocated(diexchanges))   deallocate(diexchanges)







end subroutine


end module communications
