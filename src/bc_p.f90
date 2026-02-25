module boundary
!> @brief
!> this module includes the subroutines for reading and establishing the boundary conditions
!> for all the cells that need to be bounded including the periodic ones
use declaration
use library
use transform
implicit none
contains



subroutine read_bound(n,imaxb,xmpielrank)
!> @brief
!> this subroutine reads the boundary conditions for all the bounded elements

implicit none
integer,intent(in)::n,imaxb
character(len=12)::bndfile,ibx_code
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer::i,j,ji,k,lm,iex,kmaxn,kk,kmaxe,kkk,jjj,jfx,jb,kxk,itr1,jj,jj1
integer::ioy,ibid,ib1,ib2,ib3,ib4,ibx1,itl,ibgw,ibgw2,ibleed,ibdum
integer,dimension(4)::ib_n
integer::nb1,nb2,nb3,nb4
	kmaxe=xmpielrank(n)

kkk=0
totiw=0
ibgw=0
ibgw2=0






if (dimensiona.eq.3)then

bndfile='GRID.bnd'
		if (binio.eq.0)open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
		if (binio.eq.1)open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)

itl=0

if (binio.eq.0)then
do ji=1,imaxb
	read(10,*)ibid,ib1,ib2,ib3,ib4
	
	if ((dinoder(ib1)%itor.gt.0).and.(dinoder(ib2)%itor.gt.0).and.(dinoder(ib3)%itor.gt.0)&
	.and.(dinoder(ib4)%itor.gt.0))then
	itl=itl+1

	end if
end do
else
do ji=1,imaxb
	read(10)ibid,ib1,ib2,ib3,ib4
	
	if ((dinoder(ib1)%itor.gt.0).and.(dinoder(ib2)%itor.gt.0).and.(dinoder(ib3)%itor.gt.0)&
	.and.(dinoder(ib4)%itor.gt.0))then
	itl=itl+1

	end if
end do

end if

if (itl.gt.0)then

allocate(ibound_inum(itl));ibound_inum=0
allocate(ibound_icode(itl));ibound_icode=0
allocate(ibound_ibid(itl));ibound_ibid=0
allocate( ibound_ibl(1:4,itl));ibound_ibl=0
allocate(ibound_ishape(itl));ibound_ishape=0
allocate(ibound_which(itl));ibound_which=0
allocate(ibound_face(itl));ibound_face=0
if (iperiodicity.eq.1)then
allocate(ibound_localn(2,itl));ibound_ishape=0
allocate(ibound_cpun(2,itl));ibound_ishape=0


end if


end if


 close(10)

itl=0
if (binio.eq.0)open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
if (binio.eq.1)open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)

if (binio.eq.0)then
do ji=1,imaxb
	read(10,*)ibid,ib_n(1),ib_n(2),ib_n(3),ib_n(4),ibx1
	
	if (ibx1.eq.4)then
	ibgw=ibgw+1
	end if
	
	if ((dinoder(ib_n(1))%itor.gt.0).and.(dinoder(ib_n(2))%itor.gt.0).and.(dinoder(ib_n(3))%itor.gt.0)&
	.and.(dinoder(ib_n(4))%itor.gt.0))then
	itl=itl+1
	    if (ibx1.eq.4)then
	    ibound_inum(itl)=ibgw
	    totiw=totiw+1
	    end if
		


	    ibound_icode(itl)=ibx1
	    ibound_ibid(itl)=ibid
	    if (ib_n(3).eq.ib_n(4))then
	    
	      ibound_ishape(itl)=6

			ibound_ibl(1:3,itl)=ib_n(1:3)
	    else

	      ibound_ishape(itl)=5

			ibound_ibl(1:4,itl)=ib_n(1:4)
	    end if

	end if
end do


else
	
do ji=1,imaxb
	read(10)ibid,ib_n(1),ib_n(2),ib_n(3),ib_n(4),ibx1
		
	

	if (ibx1.eq.4)then
	ibgw=ibgw+1
	end if
	
	if ((dinoder(ib_n(1))%itor.gt.0).and.(dinoder(ib_n(2))%itor.gt.0).and.(dinoder(ib_n(3))%itor.gt.0)&
	.and.(dinoder(ib_n(4))%itor.gt.0))then
	itl=itl+1
	    if (ibx1.eq.4)then
	    ibound_inum(itl)=ibgw
	    totiw=totiw+1
	    end if

	    ibound_icode(itl)=ibx1
	    ibound_ibid(itl)=ibid
	    if (ib_n(3).eq.ib_n(4))then
	    
	      ibound_ishape(itl)=6

			ibound_ibl(1:3,itl)=ib_n(1:3)
	    else

	      ibound_ishape(itl)=5

			ibound_ibl(1:4,itl)=ib_n(1:4)
	    end if

	end if
end do

end if
 close(10)
n_boundaries=itl

itl=0
do ji=1,n_boundaries
  if ((ibound_icode(ji).eq.5).or.(ibound_icode(ji).eq.50))then

    itl=itl+1
    ibound_localn(:,ji)=0
    ibound_cpun(:,ji)=0
  end if
end do


do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
  ielem_ibounds(:,i)=0
  ielem_nofbc(i)=0
  end if
end do

itl=0
jj1=0
   do ji=1,n_boundaries
	do jfx=1,dinoder2(ibound_ibl(1,ji))%numberofneib
	     j=dinoder2(ibound_ibl(1,ji))%neibids(jfx)
	      if (ielem_interior(j).eq.1)then
		  do jj=1,ielem_ifca(j)
			if (ielem_ineighg(jj,j).eq.0)then
		      if (ielem_types_faces(jj,j).eq.ibound_ishape(ji))then
				if (ibound_ishape(ji).eq.6)then
				  nb1=ielem_nodes_faces(jj,1,j)
				  nb2=ielem_nodes_faces(jj,2,j)
				  nb3=ielem_nodes_faces(jj,3,j)
				  ib1=ibound_ibl(1,ji)
				  ib2=ibound_ibl(2,ji)
				  ib3=ibound_ibl(3,ji)
				  if (((nb1.eq.ib1).or.(nb1.eq.ib2).or.(nb1.eq.ib3)).and.&
				  ((nb2.eq.ib1).or.(nb2.eq.ib2).or.(nb2.eq.ib3)).and.&
				  ((nb3.eq.ib1).or.(nb3.eq.ib2).or.(nb3.eq.ib3)))then
				 ielem_ibounds(jj,j)=ji
				 ibound_which(ji)=j
				 ibound_face(ji)=jj
				  ielem_nofbc(j)=ielem_nofbc(j)+1
				    if ((ibound_icode(ji).eq.5).or.(ibound_icode(ji).eq.50))then
				    ibound_localn(1,ji)=j;ibound_cpun(1,ji)=n
				    jj1=jj1+1
				      
				    go to 51
				    end if
				  end if
			      end if
				  if (ibound_ishape(ji).eq.5)then
				    nb1=ielem_nodes_faces(jj,1,j)
				    nb2=ielem_nodes_faces(jj,2,j)
				    nb3=ielem_nodes_faces(jj,3,j)
				    nb4=ielem_nodes_faces(jj,4,j)
				    ib1=ibound_ibl(1,ji)
				    ib2=ibound_ibl(2,ji)
				    ib3=ibound_ibl(3,ji)
				    ib4=ibound_ibl(4,ji)
					if (((nb1.eq.ib1).or.(nb1.eq.ib2).or.(nb1.eq.ib3).or.(nb1.eq.ib4)).and.&
					((nb2.eq.ib1).or.(nb2.eq.ib2).or.(nb2.eq.ib3).or.(nb2.eq.ib4)).and.&
					((nb3.eq.ib1).or.(nb3.eq.ib2).or.(nb3.eq.ib3).or.(nb3.eq.ib4)).and.&
					((nb4.eq.ib1).or.(nb4.eq.ib2).or.(nb4.eq.ib3).or.(nb4.eq.ib4)))then
					ielem_ibounds(jj,j)=ji
					ibound_which(ji)=j
				 ibound_face(ji)=jj
				    ielem_nofbc(j)=ielem_nofbc(j)+1
				    if ((ibound_icode(ji).eq.5).or.(ibound_icode(ji).eq.50))then
				    ibound_localn(1,ji)=j;ibound_cpun(1,ji)=n
				    jj1=jj1+1
				      
				    go to 51
				    end if
					end if
				    end if
			  end if
		      end if
		end do
		
	      end if
	  end do
      51 continue
    end do



     


end if

if (dimensiona.eq.2)then

bndfile='GRID.bnd'
		if (binio.eq.0)open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
		if (binio.eq.1)open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)

itl=0

if (binio.eq.0)then
do ji=1,imaxb
	read(10,*)ibid,ib1,ib2
	
	if ((dinoder(ib1)%itor.gt.0).and.(dinoder(ib2)%itor.gt.0))then
	itl=itl+1

	end if
end do
else
do ji=1,imaxb
	read(10)ibid,ib1,ib2
	
	if ((dinoder(ib1)%itor.gt.0).and.(dinoder(ib2)%itor.gt.0))then
	itl=itl+1

	end if
end do

end if

if (itl.gt.0)then

print*,"itl",itl

allocate(ibound_inum(itl));ibound_inum=0
allocate(ibound_icode(itl));ibound_icode=0
allocate(ibound_ibid(itl));ibound_ibid=0
allocate(ibound_ibl(1:2,itl));ibound_ibl=0
allocate(ibound_ishape(itl));ibound_ishape=0
allocate(ibound_which(itl));ibound_which=0
allocate(ibound_face(itl));	ibound_face=0

if (iperiodicity.eq.1)then
allocate(ibound_localn(2,itl));ibound_ishape=0
allocate(ibound_cpun(2,itl));ibound_ishape=0


end if


end if


close(10)

itl=0
if (binio.eq.0)open(10,file=bndfile,form='formatted',status='old',action='read',iostat=ioy)
if (binio.eq.1)open(10,file=bndfile,form='unformatted',status='old',action='read',iostat=ioy)

if (binio.eq.0)then
do ji=1,imaxb
	read(10,*)ibid,ib_n(1),ib_n(2),ib_n(3),ib_n(4),ibx1
	if (bleed.eq.1)then
		ibdum=ibx1
		if (ibx1.eq.4)then

		if ((dinoder(ib_n(1))%itor.gt.0).and.(dinoder(ib_n(2))%itor.gt.0))then

		!now check if this is in a bleed zone


		do ibleed=1,bleed_number
			if (((dinoder(ib_n(1))%cord(1).ge.bleed_start(ibleed,1)).and.(dinoder(ib_n(1))%cord(1).le.bleed_end(ibleed,1))).and.((dinoder(ib_n(2))%cord(1).ge.bleed_start(ibleed,1)).and.(dinoder(ib_n(2))%cord(1).le.bleed_end(ibleed,1))).and.((dinoder(ib_n(1))%cord(2).ge.bleed_start(ibleed,2)).and.(dinoder(ib_n(1))%cord(2).le.bleed_end(ibleed,2))).and.((dinoder(ib_n(2))%cord(2).ge.bleed_start(ibleed,2)).and.(dinoder(ib_n(2))%cord(2).le.bleed_end(ibleed,2))))then
			ibdum=99


			end if

		end do

		end if
		end if
		ibx1=ibdum


	end if


	if ((ibx1.eq.4).or.(ibx1.eq.99))then
	ibgw=ibgw+1
	end if
	if ((dinoder(ib_n(1))%itor.gt.0).and.(dinoder(ib_n(2))%itor.gt.0))then
	itl=itl+1
	     if ((ibx1.eq.4).or.(ibx1.eq.99))then
	    ibound_inum(itl)=ibgw;totiw=totiw+1
	    end if
	    ibound_icode(itl)=ibx1
	   ibound_ibid(itl)=ibid

	      ibound_ishape(itl)=7

			ibound_ibl(1:2,itl)=ib_n(1:2)
	   

	end if
end do
else
do ji=1,imaxb
	read(10)ibid,ib_n(1),ib_n(2),ib_n(3),ib_n(4),ibx1
	if (bleed.eq.1)then
		ibdum=ibx1
		if (ibx1.eq.4)then
		if ((dinoder(ib_n(1))%itor.gt.0).and.(dinoder(ib_n(2))%itor.gt.0))then

		!now check if this is in a bleed zone


		do ibleed=1,bleed_number
			if (((dinoder(ib_n(1))%cord(1).ge.bleed_start(ibleed,1)).and.(dinoder(ib_n(1))%cord(1).le.bleed_end(ibleed,1))).and.((dinoder(ib_n(2))%cord(1).ge.bleed_start(ibleed,1)).and.(dinoder(ib_n(2))%cord(1).le.bleed_end(ibleed,1))).and.((dinoder(ib_n(1))%cord(2).ge.bleed_start(ibleed,2)).and.(dinoder(ib_n(1))%cord(2).le.bleed_end(ibleed,2))).and.((dinoder(ib_n(2))%cord(2).ge.bleed_start(ibleed,2)).and.(dinoder(ib_n(2))%cord(2).le.bleed_end(ibleed,2))))then
			ibdum=99


			end if

		end do

		end if
		end if
		ibx1=ibdum


	end if


	if ((ibx1.eq.4).or.(ibx1.eq.99))then
	ibgw=ibgw+1
	end if
	if ((dinoder(ib_n(1))%itor.gt.0).and.(dinoder(ib_n(2))%itor.gt.0))then
	itl=itl+1
	     if ((ibx1.eq.4).or.(ibx1.eq.99))then
	    ibound_inum(itl)=ibgw;totiw=totiw+1
	    end if
	    ibound_icode(itl)=ibx1
	 ibound_ibid(itl)=ibid

	      ibound_ishape(itl)=7

			ibound_ibl(1:2,itl)=ib_n(1:2)
	   

	end if
end do

end if

close(10)

n_boundaries=itl




do ji=1,n_boundaries
  if ((ibound_icode(ji).eq.5).or.(ibound_icode(ji).eq.50))then
    ibound_localn(:,ji)=0
    ibound_cpun(:,ji)=0
  end if
end do

do i=1,kmaxe
  if (ielem_interior(i).eq.1)then


if (bleed.eq.1)then
  ielem_bleedn(:,i)=0
  end if
  ielem_ibounds(:,i)=0
  ielem_nofbc(i)=0
  end if
end do

itl=0

   do ji=1,n_boundaries
	do jfx=1,dinoder2(ibound_ibl(1,ji))%numberofneib
	     j=dinoder2(ibound_ibl(1,ji))%neibids(jfx)
	      if (ielem_interior(j).eq.1)then
		  do jj=1,ielem_ifca(j)
		     if (ielem_ineighg(jj,j).eq.0)then
				
				  nb1=ielem_nodes_faces(jj,1,j)
				  nb2=ielem_nodes_faces(jj,2,j)
				  
				  ib1=ibound_ibl(1,ji)
				  ib2=ibound_ibl(2,ji)
				  
				   if (((nb1.eq.ib1).or.(nb1.eq.ib2)).and.&
				((nb2.eq.ib1).or.(nb2.eq.ib2)))then
				 ielem_ibounds(jj,j)=ji
				 ibound_which(ji)=j
				 ibound_face(ji)=jj
				    ielem_nofbc(j)=ielem_nofbc(j)+1
				    if ((ibound_icode(ji).eq.5).or.(ibound_icode(ji).eq.50))then
				    ibound_localn(1,ji)=j;ibound_cpun(1,ji)=n
				   
				    itl=itl+1
				    end if
				  end if
			      
			  end if
		end do
	      end if
	  end do
    end do




 


  end if



totwalls=ibgw





      if (totiw.gt.0)then
	allocate(ibound_t(totiw))
	allocate(ibound_t2(totiw))
	totiw=0
	
	
	
	do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	     if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		 totiw=totiw+1
				ibound_t(totiw)=i
				ibound_t2(totiw)=j
	      end if
	  end if
	end do
   end if
end do
	
	end if



	call mpi_barrier(mpi_comm_world,ierror)
end subroutine read_bound
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine apply_boundary(n,xper,yper,zper,iperiodicity,xmpielrank)
!> @brief
!> this subroutine assigns the correct boundary condition code for all the bounded elements
implicit none
real,intent(in)::xper,yper,zper
integer,intent(in)::iperiodicity,n
real::small,tolerance,dist,temp_x
integer::i,k,j,kk,ii,kmaxe,jj1,jj2,ji,l,ibleed
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer::dum1,dum2,n_node
real,dimension(1:8,1:dimensiona)::vext,nodes_list

	kmaxe=xmpielrank(n)
 

tolerance=tolsmall !> the tolerance can have a significant impact on the periodic boundary conditions matching rules





jj1=0
if (dimensiona.eq.3)then

jj2=0
do i=1,n_boundaries
if ((ibound_icode(i).eq.5).or.(ibound_icode(i).eq.50))then
jj2=jj2+1
end if


end do


!$omp do
do i=1,kmaxe			! for all elements
    if (ielem_interior(i).eq.1)then		! that have at least one unknwon neighbour
	    if (ielem_nofbc(i).gt.0)then		! that have at least established a boundary condition code
			
		  do j=1,ielem_ifca(i)			! loop all their faces
		      if (ielem_ibounds(j,i).gt.0)then
		     if ((ibound_icode(ielem_ibounds(j,i)).eq.5).or.(ibound_icode(ielem_ibounds(j,i)).eq.50))then	!if any of them has a periodic boundary condition then
				    if (ibound_ishape(ielem_ibounds(j,i)).eq.5)then
				    n_node=4
				    else
				    n_node=3
				    end if
				    do kk=1,n_node
				      nodes_list(kk,1:3)=dinoder(ibound_ibl(kk,ielem_ibounds(j,i)))%cord(1:3)
				    end do
				    vext(1,1:3)=cordinates3(n,nodes_list,n_node)
				    
			   do ii=1,n_boundaries				! loop all the boundaries
				if ((ii.ne.ielem_ibounds(j,i)).and.((ibound_icode(ii).eq.5).or.(ibound_icode(ii).eq.50)).and.&
(ibound_ishape(ielem_ibounds(j,i)).eq.ibound_ishape(ii)))then
				      if ((ibound_localn(1,ii).gt.0)) then	! excluding itself, and of same shape type
 				   if (ielem_ihexgl(ibound_localn(1,ii)).ne.ielem_ihexgl(i))then
! 				    
				    do kk=1,n_node
				      nodes_list(kk,1:3)=dinoder(ibound_ibl(kk,ii))%cord(1:3)
				    end do
				    vext(2,1:3)=cordinates3(n,nodes_list,n_node)
				    
				    dist=distance3(n,vext)
				      
                    if(per_rot.eq.0)then
				      if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
				      ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
				
				      
! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall).or.((abs(dist-zper)).lt.tolsmall))then

				      ibound_localn(2,ii)=i;ibound_cpun(2,ii)=n
				      ielem_ineighg(j,i)=ielem_ihexgl(ibound_localn(1,ii))

			



! 				      
				      jj1=jj1+1
 				      go to 101
 				      end if
				      end if
				       if (((abs(vext(2,2)-yper).lt.tolsmall).or.(abs((abs(vext(2,2)-yper))-yper).lt.tolsmall)).and.&
				      ((abs(vext(1,2)-yper).lt.tolsmall).or.(abs((abs(vext(1,2)-yper))-yper).lt.tolsmall)))then
				      if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
				      
! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall).or.((abs(dist-zper)).lt.tolsmall))then

				      ibound_localn(2,ii)=i;ibound_cpun(2,ii)=n
				      ielem_ineighg(j,i)=ielem_ihexgl(ibound_localn(1,ii))
				      jj1=jj1+1
 				      go to 101
 				      end if
				      end if
				       if (((abs(vext(2,3)-zper).lt.tolsmall).or.(abs((abs(vext(2,3)-zper))-zper).lt.tolsmall)).and.&
				      ((abs(vext(1,3)-zper).lt.tolsmall).or.(abs((abs(vext(1,3)-zper))-zper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,1)-vext(1,1)).lt.tolsmall))then
				      
! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall).or.((abs(dist-zper)).lt.tolsmall))then

				      ibound_localn(2,ii)=i;ibound_cpun(2,ii)=n
				      ielem_ineighg(j,i)=ielem_ihexgl(ibound_localn(1,ii))
				      jj1=jj1+1
 				      go to 101
 				      end if
				      end if
				    else
                        vext(2,:)=rotate_per(vext(2,:),ibound_icode(ii),angle_per)
                        if ((abs(vext(1,1)-vext(2,1)).lt.tol_per).and.&
                            (abs(vext(1,2)-vext(2,2)).lt.tol_per).and.&
                            (abs(vext(1,3)-vext(2,3)).lt.tol_per)) then              
                                ibound_localn(2,ii)=i
                                ibound_cpun(2,ii)=n
                                ielem_ineighg(j,i)=ielem_ihexgl(ibound_localn(1,ii))
                                jj1=jj1+1
                            go to 101
 				      end if
				    end if
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
!$omp end do


	 	
		    
	else
jj2=0
do i=1,n_boundaries
if (ibound_icode(i).eq.5)then
jj2=jj2+1
end if
end do
jj1=0
!$omp do
do i=1,kmaxe			!> all elements
    if (ielem_interior(i).eq.1)then		! that have at least one unknwon neighbour
	    if (ielem_nofbc(i).gt.0)then		! that have at least established a boundary condition code
		  do j=1,ielem_ifca(i)			! loop all their boundary faces
		      if (ielem_ibounds(j,i).gt.0)then


			  if (bleed.eq.1)then
		      if (ibound_icode(ielem_ibounds(j,i)).eq.99)then


					!assign the bleed zone
					n_node=2
				    do kk=1,n_node
				      nodes_list(kk,1:2)=dinoder(ibound_ibl(kk,ielem_ibounds(j,i)))%cord(1:2)
				    end do
				    vext(1,:)=cordinates2(n,nodes_list,n_node)

					do ibleed=1,bleed_number
						if (((vext(1,1).ge.bleed_start(ibleed,1)).and.(vext(1,1).le.bleed_end(ibleed,1))).and.((vext(1,2).ge.bleed_start(ibleed,2)).and.(vext(1,2).le.bleed_end(ibleed,2))))then

						ielem_bleedn(j,i)=ibleed	!assign the bleed number
						end if
					end do

				end if
				end if


		     if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	! if any of them has a periodic boundary condition then
		     
				    n_node=2
				    do kk=1,n_node
				      nodes_list(kk,1:2)=dinoder(ibound_ibl(kk,ielem_ibounds(j,i)))%cord(1:2)
				    end do
				    vext(1,:)=cordinates2(n,nodes_list,n_node)
			   do ii=1,n_boundaries				! loop all the boundaries
				if (((ii.ne.ielem_ibounds(j,i)).and.(ibound_icode(ii).eq.5)))then
				      if(ibound_localn(1,ii).gt.0)then	! excluding itself, and of same shape type
				    if (ielem_ihexgl(ibound_localn(1,ii)).ne.ielem_ihexgl(i))then
				    do kk=1,n_node
				      nodes_list(kk,1:2)=dinoder(ibound_ibl(kk,ii))%cord(1:2)
				    end do
				    vext(2,:)=cordinates2(n,nodes_list,n_node)
				      dist=distance2(n,vext)
				      


				      if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
				      ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall))then


! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall))then
						
				      ibound_localn(2,ii)=i;ibound_cpun(2,ii)=n
				      ielem_ineighg(j,i)=ielem_ihexgl(ibound_localn(1,ii))
				      jj1=jj1+1
					go to 201
				      end if
				      end if
				      if (((abs(vext(2,2)-yper).lt.tolsmall).or.(abs((abs(vext(2,2)-yper))-yper).lt.tolsmall)).and.&
				      ((abs(vext(1,2)-yper).lt.tolsmall).or.(abs((abs(vext(1,2)-yper))-yper).lt.tolsmall)))then
				      if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall))then
				      !

! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall))then
						
				      ibound_localn(2,ii)=i;ibound_cpun(2,ii)=n
				      ielem_ineighg(j,i)=ielem_ihexgl(ibound_localn(1,ii))
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
!$omp end do 

  
	
end if    
	      


      



	
    

end subroutine apply_boundary





















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
end module boundary
