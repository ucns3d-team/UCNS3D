module local
use library
use declaration
use transform
implicit none

contains

subroutine exch_cords(n)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values
implicit none
integer,intent(in)::n
integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
real,dimension(1)::dumts,rumts
integer,dimension(4)::icfv1,icfv2
real::rcfv1,rcfv2

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
kmaxe=xmpielrank(n)


indl=diexchanger(1)%tot
tndl=diexchanges(1)%tot



if (fastest.ne.1)then
ineedt=direcexr(1)%tot
tneedt=direcexs(1)%tot
allocate (diexcordr(ineedt))
allocate (diexcords(tneedt))
allocate (diexsolhir(ineedt))
allocate (diexsolhis(tneedt))




if (adda.eq.1)then
allocate (diexsolhird(ineedt))
allocate (diexsolhisd(tneedt))
end if

end if
allocate (diexboundhir(indl))
allocate (diexboundhis(tndl))
allocate (diexboundhirr(indl))
allocate (diexboundhiss(tndl))

call mpi_barrier(mpi_comm_world,ierror)

if (dimensiona.eq.3)then
i_cnt2=4;i_cnt3=3;i_cnt4=8
else
i_cnt2=2;i_cnt3=2;i_cnt4=4
end if
i_cnt5=i_cnt3*i_cnt4
if (fastest.ne.1)then
do i=1,ineedt

	diexsolhir(i)%procid=direcexr(i)%procid
 	diexcordr(i)%procid=direcexr(i)%procid
 	allocate (diexcordr(i)%nodecord(direcexr(i)%muchineed(1),i_cnt4,i_cnt3))
 	diexcordr(i)%nodecord(1:direcexr(i)%muchineed(1),i_cnt4,i_cnt3)=-tolbig

		 allocate(diexsolhir(i)%sol(direcexr(i)%muchineed(1),nof_variables+turbulenceequations+passivescalar))

		 if (adda.eq.1)then
		allocate(diexsolhird(i)%sol(direcexr(i)%muchineed(1),1))
		end if
	diexsolhir(i)%sol(:,:)=0.0d0

end do
end if
do i=1,indl
	diexboundhir(i)%procid=diexchanger(i)%procid
	diexboundhirr(i)%procid=diexchanger(i)%procid
	if (itestcase.le.3)then
! 	allocate(diexboundhir(i)%facesol(diexchanger(i)%muchineed(1),nof_variables))
	allocate(diexboundhirr(i)%vertpp(diexchanger(i)%muchineed(1),i_cnt2))

	else

	  if (dimensiona.eq.3)then
	 i_cnt=(nof_variables+turbulenceequations+passivescalar)+((nof_variables-1+turbulenceequations+passivescalar)*3)
	  else
	  i_cnt=(nof_variables+turbulenceequations+passivescalar)+((nof_variables-1+turbulenceequations+passivescalar)*2)
	  end if

! 	    allocate(diexboundhir(i)%facesol(diexchanger(i)%muchineed(1),i_cnt))

	   allocate(diexboundhirr(i)%vertpp(diexchanger(i)%muchineed(1),i_cnt2))
	end if

! 	diexboundhir(i)%facesol(:,:)=0.0d0
	diexboundhirr(i)%vertpp(:,:)=0

end do
if (fastest.ne.1)then
do i=1,tneedt

	diexsolhis(i)%procid=direcexs(i)%procid

	diexcords(i)%procid=direcexs(i)%procid
	allocate (diexcords(i)%nodecord(direcexs(i)%muchtheyneed(1),i_cnt4,i_cnt3))



	diexcords(i)%nodecord(1:direcexs(i)%muchtheyneed(1),1:i_cnt4,1:i_cnt3)=-tolbig
	allocate (diexsolhis(i)%sol(direcexs(i)%muchtheyneed(1),nof_variables+turbulenceequations+passivescalar))

	if (adda.eq.1)then
	allocate(diexsolhisd(i)%sol(direcexs(i)%muchtheyneed(1),1))
	end if

	diexsolhis(i)%sol(:,:)=0.0d0

end do
end if
do i=1,tndl


      diexboundhis(i)%procid=diexchanges(i)%procid
	diexboundhiss(i)%procid=diexchanges(i)%procid
	if (itestcase.le.3)then
! 	allocate(diexboundhis(i)%facesol(diexchanges(i)%muchtheyneed(1),nof_variables))

	allocate(diexboundhiss(i)%vertpp(diexchanges(i)%muchtheyneed(1),i_cnt2))

	else

	  if (dimensiona.eq.3)then
	 i_cnt=(nof_variables+turbulenceequations+passivescalar)+((4+turbulenceequations+passivescalar)*3)
	  else
	  i_cnt=(nof_variables+turbulenceequations+passivescalar)+((3+turbulenceequations+passivescalar)*2)

	  end if
! 	    allocate(diexboundhis(i)%facesol(diexchanges(i)%muchtheyneed(1),i_cnt))

	   allocate(diexboundhiss(i)%vertpp(diexchanges(i)%muchtheyneed(1),i_cnt2))
	end if

! 	diexboundhis(i)%facesol(:,:)=0.0d0
	diexboundhiss(i)%vertpp(:,:)=0
end do




call exchange_cordx(n,ineedt,tneedt,i_cnt4,i_cnt3)









end subroutine exch_cords



subroutine exchange_cordx(n,ineedt,tneedt,i_cnt4,i_cnt3)
integer,intent(in)::n,ineedt,tneedt,i_cnt4,i_cnt3
integer::i,j,k,l,icpuid,ixflag,itee,iteedum,iavc,iavt
real,dimension(1)::dumts,rumts

if (fastest.ne.1)then

do i=1,tneedt
!
	do k=1,direcexs(i)%muchtheyneed(1)
!
 	    do j=1,ielem_nonodes(direcexs(i)%localref(k))
!
		diexcords(i)%nodecord(k,j,1:dims)=dinoder(ielem_nodes(j,direcexs(i)%localref(k)))%cord(1:dims)
!
	    end do
	end do
end do




call mpi_barrier(mpi_comm_world,ierror)
icpuid=n
dumts(1:1)=tolsmall
		do i=0,isize-1
			if (i.ne.n) then
				do j=1,tneedt
					iavt=10000
					if (direcexs(j)%procid.eq.i)then
					iavt=j
					go to 7001
					end if
				end do
				7001 continue
				do k=1,ineedt
					iavc=10000
					if (direcexr(k)%procid.eq.i) then
					iavc=k
					go to 8001
					end if
				end do
				8001 continue
				if ((iavt.eq.10000).and.(iavc.ne.10000)) then
				call mpi_sendrecv(dumts(1:1),1,mpi_double_precision,i,icpuid,&
				diexcordr(iavc)%nodecord(1:direcexr(iavc)%muchineed(1),1:i_cnt4,1:i_cnt3),&
direcexr(iavc)%muchineed(1)*i_cnt4*i_cnt3,mpi_double_precision,diexcordr(iavc)%procid,diexcordr(iavc)%procid,mpi_comm_world,status,ierror)

				end if
				if ((iavt.ne.10000).and.(iavc.eq.10000)) then
				!message 1
				call mpi_sendrecv(diexcords(iavt)%nodecord(1:direcexs(iavt)%muchtheyneed(1),1:i_cnt4,1:i_cnt3),&
direcexs(iavt)%muchtheyneed(1)*i_cnt4*i_cnt3,mpi_double_precision,diexcords(iavt)%procid,icpuid,&
				dumts(1:1),1,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
				end if
				if ((iavt.ne.10000).and.(iavc.ne.10000)) then
				!message 1
				call mpi_sendrecv(diexcords(iavt)%nodecord(1:direcexs(iavt)%muchtheyneed(1),1:i_cnt4,1:i_cnt3),&
direcexs(iavt)%muchtheyneed(1)*i_cnt4*i_cnt3,mpi_double_precision,diexcords(iavt)%procid,icpuid,&
				diexcordr(iavc)%nodecord(1:direcexr(iavc)%muchineed(1),1:i_cnt4,1:i_cnt3),direcexr(iavc)%muchineed(1)*i_cnt4*i_cnt3,&
mpi_double_precision,diexcordr(iavc)%procid,diexcordr(iavc)%procid,mpi_comm_world,status,ierror)
!


				end if


			!	end do !k
			!end do! j
			end if	! i.ne.n
		end do
end if
call mpi_barrier(mpi_comm_world,ierror)



end subroutine exchange_cordx


subroutine exch_cords_opt(n)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values
implicit none
integer,intent(in)::n
integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
real,dimension(1)::dumts,rumts
integer,dimension(4)::icfv1,icfv2
real::rcfv1,rcfv2

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
kmaxe=xmpielrank(n)


indl=diexchanger(1)%tot
tndl=diexchanges(1)%tot



if (dimensiona.eq.3)then
i_cnt2=4;i_cnt3=3;i_cnt4=8
else
i_cnt2=2;i_cnt3=2;i_cnt4=4
end if
i_cnt5=i_cnt3*i_cnt4

do i=1,indl

	if (itestcase.le.3)then
        allocate(diexboundhir(i)%facesol(diexchanger(i)%muchineed(1),nof_variables))
    ! 	allocate(diexboundhirr(i)%vertpp(diexchanger(i)%muchineed(1),i_cnt2))


        if (dg == 1)then
            allocate(diexboundhir(i)%facesol_dg(diexchanger(i)%muchineed(1),nof_variables))
        end if
        
        if (mood.eq.1)then
            allocate(diexboundhir(i)%facesol_m(diexchanger(i)%muchineed(1),1))
        end if

	else
	
	  if (dimensiona.eq.3)then
	 i_cnt=(nof_variables+turbulenceequations+passivescalar)+((nof_variables-1+turbulenceequations+passivescalar)*3)
	  else
	  i_cnt=(nof_variables+turbulenceequations+passivescalar)+((nof_variables-1+turbulenceequations+passivescalar)*2)
	  end if

	    allocate(diexboundhir(i)%facesol(diexchanger(i)%muchineed(1),i_cnt))
        if (dg == 1)then
            allocate(diexboundhir(i)%facesol_dg(diexchanger(i)%muchineed(1),i_cnt))
        end if


            
	end if

	diexboundhir(i)%facesol(:,:)=0.0d0
	if (dg.eq.1)then
	diexboundhir(i)%facesol_dg(:,:)=0.0d0
	end if
! 	diexboundhirr(i)%vertpp(:,:)=0

end do

do i=1,tndl
	if (itestcase.le.3)then
        allocate(diexboundhis(i)%facesol(diexchanges(i)%muchtheyneed(1),nof_variables))
        
        if (dg == 1)then
            allocate(diexboundhis(i)%facesol_dg(diexchanges(i)%muchtheyneed(1),nof_variables))
        end if
        
        if (mood.eq.1)then
            allocate(diexboundhis(i)%facesol_m(diexchanges(i)%muchtheyneed(1),1))
        end if

	else
	
	  if (dimensiona.eq.3)then
	 i_cnt=(nof_variables+turbulenceequations+passivescalar)+((nof_variables-1+turbulenceequations+passivescalar)*3)
	  else
	  i_cnt=(nof_variables+turbulenceequations+passivescalar)+((nof_variables-1+turbulenceequations+passivescalar)*2)

	  end if
	    allocate(diexboundhis(i)%facesol(diexchanges(i)%muchtheyneed(1),i_cnt))

	   if (dg == 1)then
            allocate(diexboundhis(i)%facesol_dg(diexchanges(i)%muchtheyneed(1),i_cnt))
        end if
	    
	  
	end if

	diexboundhis(i)%facesol(:,:)=0.0d0
	if (dg.eq.1)then
	diexboundhis(i)%facesol_dg(:,:)=0.0d0
	end if

	
end do


call mpi_barrier(mpi_comm_world,ierror)











		

end subroutine exch_cords_opt





subroutine exch_cord3(n)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values and mapping of gaussian quadrature points
implicit none
integer,intent(in)::n
integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
real,dimension(1)::dumts,rumts
integer,dimension(4)::icfv1,icfv2
real::rcfv1,rcfv2

icpuid=n

kmaxe=xmpielrank(n)

indl=diexchanger(1)%tot
tndl=diexchanges(1)%tot

if (dimensiona.eq.3)then
i_cnt2=4;i_cnt3=3;i_cnt4=8
else
i_cnt2=2;i_cnt3=2;i_cnt4=4
end if
i_cnt5=i_cnt3*i_cnt4

if (dimensiona.eq.3)then



!now try to remap the local and global neighbours and mapping of the gaussian quadrature points


!first the global
do j=1,indl
do i=1,tndl
if (diexchanger(j)%procid.eq.diexchanges(i)%procid)then
do k=1,diexchanges(i)%muchtheyneed(1)

! 	    
	    if (ielem_types_faces(diexchanges(i)%sidetheyneed(k),(diexchanges(i)%localref(k))).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
! 	    
	    diexboundhiss(i)%vertpp(k,1:ixf4)=ielem_nodes_faces(diexchanges(i)%sidetheyneed(k),1:ixf4,(diexchanges(i)%localref(k)))


	    

end do
end if
end do
end do

call mpi_barrier(mpi_comm_world,ierror)

do k=1,indl
do j=1,tndl
      if (diexchanger(k)%procid.eq.diexchanges(j)%procid)then

      call mpi_sendrecv(diexboundhiss(j)%vertpp(1:diexchanges(j)%muchtheyneed(1),1:i_cnt2)&
,diexchanges(j)%muchtheyneed(1)*i_cnt2,mpi_integer,diexchanges(j)%procid,&
      diexchanges(j)%procid,diexboundhirr(k)%vertpp(1:diexchanger(k)%muchineed(1),1:i_cnt2),&
diexchanger(k)%muchineed(1)*i_cnt2,mpi_integer,diexchanges(j)%procid,&
      icpuid,mpi_comm_world,status,ierror)
      end if
end do
end do







do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
  do k=1,ielem_ifca(i)
		  if (ielem_types_faces(k,i).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
	    
      if (ielem_ineighg(k,i).gt.0)then
	if (ielem_ineighb(k,i).ne.n)then
	    
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
! 		ielem_nodes_faces(k,1:ixf4,i)=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_q_face_q_mapl(1,k,i),1:ixf4)
! 		ielem_reorient(k,i)=1
	    if (ielem_ibounds(k,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(k,i)).eq.5).or.(ibound_icode(ielem_ibounds(k,i)).eq.50))then
		do ixfv=1,ixf4
! 		dinoder(diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_q_face_q_mapl(1,k,i),ixfv))%itor=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_q_face_q_mapl(1,k,i),ixfv)
		dinoder(diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),ixfv))%itor=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),ixfv)

		end do


	      end if
	    end if
	    end if

	end if
	end if
	
      end do
    end if
  
end do










else







!now try to remap the local and global neighbours and mapping of the gaussian quadrature points


!first the global
do j=1,indl
do i=1,tndl
if (diexchanger(j)%procid.eq.diexchanges(i)%procid)then
do k=1,diexchanges(i)%muchtheyneed(1)

ixf4=2

diexboundhiss(i)%vertpp(k,1:ixf4)=ielem_nodes_faces(diexchanges(i)%sidetheyneed(k),1:ixf4,(diexchanges(i)%localref(k)))
end do
end if
end do
end do


do k=1,indl
do j=1,tndl
      if (diexboundhirr(k)%procid.eq.diexboundhiss(j)%procid)then
      call mpi_sendrecv(diexboundhiss(j)%vertpp(1:diexchanges(j)%muchtheyneed(1),1:i_cnt2)&
,diexchanges(j)%muchtheyneed(1)*i_cnt2,mpi_integer,diexboundhiss(j)%procid,&
      icpuid,diexboundhirr(k)%vertpp(1:diexchanger(k)%muchineed(1),1:i_cnt2),&
diexchanger(k)%muchineed(1)*i_cnt2,mpi_integer,diexboundhirr(k)%procid,&
      diexboundhirr(k)%procid,mpi_comm_world,status,ierror)
      end if
end do
end do

















do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
  do k=1,ielem_ifca(i)
		 
	    ixf4=2
	    
      if (ielem_ineighg(k,i).gt.0)then
	if (ielem_ineighb(k,i).ne.n)then
	    
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then

! 		ielem_nodes_faces(k,1:ixf4,i)=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_q_face_q_mapl(1,k,i),1:ixf4)
! 		ielem_reorient(k,i)=1
	    if (ielem_ibounds(k,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(k,i)).eq.5).or.(ibound_icode(ielem_ibounds(k,i)).eq.50))then
		do ixfv=1,ixf4
 		dinoder(diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),ixfv))%itor=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),ixfv)



		end do
	      end if
	    end if
	    end if


	else
! 	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
! 		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
! 		ielem_reorient(k,i)=1
! 	    end if
	end if
      end if
  end do
  else
!       do k=1,ielem_ifca(i)
! 
! 
! 	  ixf4=2
!       if (ielem_ineighg(k,i).gt.0)then
! 	
! 	
! 	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
! 		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
! 		ielem_reorient(k,i)=1
! 	    end if
! 	
!       end if
!       end do
  end if
end do


end if

end subroutine exch_cord3


subroutine exch_cords2(n,isize,diexboundhiri,diexboundhisi,&
itestcase,numberofpoints2,diexchanger,diexchanges)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values in 2d
implicit none
type(aexchange_boundhi),allocatable,dimension(:),intent(inout)::diexboundhiri
type(aexchange_boundhi),allocatable,dimension(:),intent(inout)::diexboundhisi
type(aexchange),allocatable,dimension(:),intent(in)::diexchanger
type(aexchange),allocatable,dimension(:),intent(in)::diexchanges
integer,intent(in)::isize,n,itestcase,numberofpoints2
integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt
real,dimension(1:1)::dumts,rumts

indl=diexchanger(1)%tot
tndl=diexchanges(1)%tot


allocate (diexboundhiri(indl))
allocate (diexboundhisi(tndl))


i_cnt=(nof_variables+turbulenceequations+passivescalar)
do i=1,indl
	diexboundhiri(i)%procid=diexchanger(i)%procid
	
	
	allocate(diexboundhiri(i)%facesol(diexchanger(i)%muchineed(1),i_cnt))
	diexboundhiri(i)%facesol(:,:)=0.0d0
	
end do

do i=1,tndl
	diexboundhisi(i)%procid=diexchanges(i)%procid
	
	allocate(diexboundhisi(i)%facesol(diexchanges(i)%muchtheyneed(1),i_cnt))
	
	diexboundhisi(i)%facesol(:,:)=0.0d0
	
end do


call mpi_barrier(mpi_comm_world,ierror)
		

end subroutine exch_cords2














subroutine localise_stencil(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell
implicit none
integer,intent(in)::n,iconsi
integer::i,j,l,ixff,ixsst,ikg,ismp,inv,ineedt,ikg2,in_sten,k,itarget
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexg  !global index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexl  !local index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
integer,allocatable,dimension(:,:),intent(inout)::ilox_ishape !shape of each element
real,allocatable,dimension(:,:),intent(inout)::ilox_xxc       !cell centre coordinates in x
real,allocatable,dimension(:,:),intent(inout)::ilox_yyc       !cell centre coordinates in y
real,allocatable,dimension(:,:),intent(inout)::ilox_zzc      !cell centre coordinates in z
real,allocatable,dimension(:,:),intent(inout)::ilox_volume    !cell volume
integer,allocatable,dimension(:,:),intent(inout)::ilox_periodicflag
integer,allocatable,dimension(:,:,:),intent(inout)::ilon_nodcount  !number of nodes
real,allocatable,dimension(:,:,:),intent(inout)::ilon_x           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_y           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_z           !coordinates of each node in x
real::rin_sten
real,dimension(1:dimensiona)::cords
ineedt=direcexr(1)%tot
i=iconsi
	ikg2=0
	do ismp=1,typesten
		ikg=0
		
                        if ((ees.ne.5).or.(ismp.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        end if
                        
		
			do l=1,itarget
				if (ilocalstencil(n,i,ismp,l).gt.0)then
				ikg=ikg+1
				end if
			end do
			if (ikg.eq.itarget)then
			
			
			
				if (rec_local(i).eq.0)then
				ikg2=ikg2+1
				
					do l=1,itarget
					ilox_ihexg(ikg2,l)=ilocalstencil(n,i,ismp,l)
					ilox_periodicflag(ikg2,l)=ilocalstencilper(n,i,ismp,l)
					j=(xmpil(ilox_ihexg(ikg2,l)))
					ilox_ihexl(ikg2,l)=j
					ilox_ishape(ikg2,l)=ielem_ishape(j)
					call compute_centre3d(j,cords)
					ilox_xxc(ikg2,l)=cords(1)
					ilox_yyc(ikg2,l)=cords(2)
					ilox_zzc(ikg2,l)=cords(3)

					ilon_nodcount(ikg2,l,1:ielem_nonodes(j))=ielem_nodes(1:ielem_nonodes(j),j)

					do k=1,ielem_nonodes(j)
					ilon_x(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(1)
					ilon_y(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(2)
					ilon_z(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(3)
					end do
					end do
				
				
				
				else
				ikg2=ikg2+1
				do l=1,itarget
				ilox_ihexg(ikg2,l)=ilocalstencil(n,i,ismp,l)
				ilox_periodicflag(ikg2,l)=ilocalstencilper(n,i,ismp,l)
				ilox_ihexb(ikg2,l)=xmpie(ilocalstencil(n,i,ismp,l))
				if (ilox_ihexb(ikg2,l).eq.n)then
				      j=(xmpil(ilox_ihexg(ikg2,l)))
						if (ilox_ihexg(ikg2,l).eq.ielem_ihexgl(j))then
						ilox_ihexl(ikg2,l)=j
						ilox_ishape(ikg2,l)=ielem_ishape(j)
						end if
						call compute_centre3d(j,cords)
						ilox_xxc(ikg2,l)=cords(1)
						ilox_yyc(ikg2,l)=cords(2)
						ilox_zzc(ikg2,l)=cords(3)
						ilon_nodcount(ikg2,l,1:ielem_nonodes(j))=ielem_nodes(1:ielem_nonodes(j),j)
						do k=1,ielem_nonodes(j)    
				ilon_x(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(1)
				ilon_y(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(2)
				ilon_z(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(3)
				
				
				
				end do
				
				
				
				
				else
				
				
				
					
					 do ixff=1,ineedt
						if (diexcordr(ixff)%procid.eq.ilox_ihexb(ikg2,l))then
							do ixsst=1,direcexr(ixff)%muchineed(1)
							      if ((ilox_ihexg(ikg2,l)).eq.direcexr1(ixff)%whatineed(ixsst))then
								ilox_ihexl(ikg2,l)=ixsst
								ilox_ihexn(ikg2,l)=ixff
								ilox_ishape(ikg2,l)=direcexr1(ixff)%ishape(ixsst)
								exit
							       end if
							end do
						end if
					 end do	
					select case(ilox_ishape(ikg2,l))
					case(1)
					in_sten=8
					rin_sten=8.0d0
					case(2)
					in_sten=4
					rin_sten=4.0d0

					case(3)
					in_sten=5
					rin_sten=5.0d0

					case(4)
					in_sten=6
					rin_sten=6.0d0
					end select
					ilon_x(ikg2,l,1:in_sten)=diexcordr(ilox_ihexn(ikg2,l))%nodecord(ilox_ihexl(ikg2,l),1:in_sten,1)
					ilon_y(ikg2,l,1:in_sten)=diexcordr(ilox_ihexn(ikg2,l))%nodecord(ilox_ihexl(ikg2,l),1:in_sten,2)
					ilon_z(ikg2,l,1:in_sten)=diexcordr(ilox_ihexn(ikg2,l))%nodecord(ilox_ihexl(ikg2,l),1:in_sten,3)
					ilox_xxc(ikg2,l)=sum(ilon_x(ikg2,l,1:in_sten))/rin_sten
					ilox_yyc(ikg2,l)=sum(ilon_y(ikg2,l,1:in_sten))/rin_sten
					ilox_zzc(ikg2,l)=sum(ilon_z(ikg2,l,1:in_sten))/rin_sten
					
					
					
					
					
					
! 			
					
				end if
				end do
				end if
				end if
				end do














end subroutine localise_stencil


subroutine localise_stencil2d(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
implicit none
integer,intent(in)::n,iconsi
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexg  !global index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexl  !local index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
integer,allocatable,dimension(:,:),intent(inout)::ilox_ishape !shape of each element
real,allocatable,dimension(:,:),intent(inout)::ilox_xxc       !cell centre coordinates in x
real,allocatable,dimension(:,:),intent(inout)::ilox_yyc       !cell centre coordinates in y
real,allocatable,dimension(:,:),intent(inout)::ilox_zzc      !cell centre coordinates in z
real,allocatable,dimension(:,:),intent(inout)::ilox_volume    !cell volume
integer,allocatable,dimension(:,:),intent(inout)::ilox_periodicflag
integer,allocatable,dimension(:,:,:),intent(inout)::ilon_nodcount  !number of nodes
real,allocatable,dimension(:,:,:),intent(inout)::ilon_x           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_y           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_z           !coordinates of each node in x
real,dimension(1:dimensiona)::cords
integer::i,j,l,ixff,ixsst,ikg,ismp,inv,ineedt,ikg2,in_sten,k,itarget,nojcount
real::rin_sten
ineedt=direcexr(1)%tot
i=iconsi
	ikg2=0
	do ismp=1,typesten
		ikg=0
                             if ((ees.ne.5).or.(ismp.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        end if
		
		
			do l=1,itarget
				if (ilocalstencil(n,i,ismp,l).gt.0)then
				ikg=ikg+1
				end if
			end do
			if (ikg.eq.itarget)then
				if (rec_local(i).eq.0)then
				ikg2=ikg2+1
				do l=1,itarget
				ilox_ihexg(ikg2,l)=ilocalstencil(n,i,ismp,l)
				j=(xmpil(ilox_ihexg(ikg2,l)))
				ilox_ihexl(ikg2,l)=j
				ilox_ishape(ikg2,l)=ielem_ishape(j)
				call compute_centre2d(j,cords)
				ilox_xxc(ikg2,l)=cords(1)
				ilox_yyc(ikg2,l)=cords(2)
				
				
				ilon_nodcount(ikg2,l,1:ielem_nonodes(j))=ielem_nodes(1:ielem_nonodes(j),j)
				do k=1,ielem_nonodes(j)    
				ilon_x(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(1)
				ilon_y(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(2)
				
				end do
				end do
				else
				ikg2=ikg2+1
				do l=1,itarget
				ilox_ihexg(ikg2,l)=ilocalstencil(n,i,ismp,l)
				ilox_ihexb(ikg2,l)=xmpie(ilocalstencil(n,i,ismp,l))
				if (ilox_ihexb(ikg2,l).eq.n)then
				      j=(xmpil(ilox_ihexg(ikg2,l)))
						if (ilox_ihexg(ikg2,l).eq.ielem_ihexgl(j))then
						ilox_ihexl(ikg2,l)=j
						ilox_ishape(ikg2,l)=ielem_ishape(j)
						end if
						call compute_centre2d(j,cords)
						ilox_xxc(ikg2,l)=cords(1)
						ilox_yyc(ikg2,l)=cords(2)
						
						ilon_nodcount(ikg2,l,1:ielem_nonodes(j))=ielem_nodes(1:ielem_nonodes(j),j)
						do k=1,ielem_nonodes(j)    
				ilon_x(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(1)
				ilon_y(ikg2,l,k)=dinoder(ielem_nodes(k,j))%cord(2)
				
				end do
				else
					 do ixff=1,ineedt
						if (diexcordr(ixff)%procid.eq.ilox_ihexb(ikg2,l))then
							do ixsst=1,direcexr(ixff)%muchineed(1)
							      if ((ilox_ihexg(ikg2,l)).eq.direcexr1(ixff)%whatineed(ixsst))then
								ilox_ihexl(ikg2,l)=ixsst
								ilox_ihexn(ikg2,l)=ixff
								ilox_ishape(ikg2,l)=direcexr1(ixff)%ishape(ixsst)
								exit
							       end if
							end do
						end if
					 end do	
					select case(ilox_ishape(ikg2,l))
					case(5)
					in_sten=4
					rin_sten=4.d0
					case(6)
					in_sten=3
					rin_sten=3.d0

					
					end select
					ilon_x(ikg2,l,1:in_sten)=diexcordr(ilox_ihexn(ikg2,l))%nodecord(ilox_ihexl(ikg2,l),1:in_sten,1)
					ilon_y(ikg2,l,1:in_sten)=diexcordr(ilox_ihexn(ikg2,l))%nodecord(ilox_ihexl(ikg2,l),1:in_sten,2)
					
					ilox_xxc(ikg2,l)=sum(ilon_x(ikg2,l,1:in_sten))/rin_sten
					ilox_yyc(ikg2,l)=sum(ilon_y(ikg2,l,1:in_sten))/rin_sten
					
				end if
				end do
				end if
				end if
				end do

! 				if (initcond.gt.100000)then
! 				write(200+n,*)"#element number",ielem_ihexgl,i,itarget
! 				do l=1,itarget
! 				write(200+n,*)ilox_xxc(1,l),ilox_yyc(1,l)
! 				end do
! 				end if

				

end subroutine localise_stencil2d




subroutine  localise_sten2(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell
implicit none
integer :: mm_i
integer :: mm_j
real :: mm_old(1:dimensiona)
integer,intent(in)::iconsi,n
integer::i,j,k,l,kk,prk,jj,kmaxe,ineedt,jx2,jx,iivd,iivd3,facexx,ixxfff,in1,itarget,idum,eltype,n_node,elem_dec
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexg  !global index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexl  !local index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
integer,allocatable,dimension(:,:),intent(inout)::ilox_ishape !shape of each element
real,allocatable,dimension(:,:),intent(inout)::ilox_xxc       !cell centre coordinates in x
real,allocatable,dimension(:,:),intent(inout)::ilox_yyc       !cell centre coordinates in y
real,allocatable,dimension(:,:),intent(inout)::ilox_zzc      !cell centre coordinates in z
real,allocatable,dimension(:,:),intent(inout)::ilox_volume    !cell volume
integer,allocatable,dimension(:,:),intent(inout)::ilox_periodicflag
integer,allocatable,dimension(:,:,:),intent(inout)::ilon_nodcount  !number of nodes
real,allocatable,dimension(:,:,:),intent(inout)::ilon_x           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_y           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_z           !coordinates of each node in x
real,dimension(3)::tempcentres,temp_cg
real::dumv1,dumv2,detjc,dist1,ma,mb,mc,md,me,mf,mg,mh,mi,mdd,tempxx
real,dimension(1:dimensiona)::cords
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real,dimension(1:dimensiona,1:dimensiona)::vva1
real,dimension(1)::deta

ineedt=direcexr(1)%tot
i=iconsi




	if (iperiodicity.eq.1)then

	  do jj=1,ielem_admis(i)
	  
             if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        end if
	  
	  
	    do j=2,itarget
	    if(per_rot.eq.0)then
	    if(abs(ilox_xxc(jj,j)-ilox_xxc(jj,1)).gt.xper*oo2)then
		    ilox_xxc(jj,j)=ilox_xxc(jj,j)+(xper*sign(1.0,ilox_xxc(jj,1)-xper*oo2))
		    do kk=1,8
		    ilon_x(jj,j,kk)=ilon_x(jj,j,kk)+(xper*sign(1.0,ilox_xxc(jj,1)-xper*oo2))
		    end do
	    end if
	    if(abs(ilox_yyc(jj,j)-ilox_yyc(jj,1)).gt.yper*oo2)then
		    ilox_yyc(jj,j)=ilox_yyc(jj,j)+(yper*sign(1.0,ilox_yyc(jj,1)-yper*oo2))
		  
		    do kk=1,8
		    ilon_y(jj,j,kk)=ilon_y(jj,j,kk)+(yper*sign(1.0,ilox_yyc(jj,1)-yper*oo2))
		    end do
		    
	    end if
	    if(abs(ilox_zzc(jj,j)-ilox_zzc(jj,1)).gt.zper*oo2)then
		    ilox_zzc(jj,j)=ilox_zzc(jj,j)+(zper*sign(1.0,ilox_zzc(jj,1)-zper*oo2))
		      do kk=1,8
		    ilon_z(jj,j,kk)=ilon_z(jj,j,kk)+(zper*sign(1.0,ilox_zzc(jj,1)-zper*oo2))
		    end do
		   
	    end if
	    else
            if (ilox_periodicflag(jj,j).eq.2) then
                    tempxx=ilox_xxc(jj,j)
                    ilox_xxc(jj,j)=tempxx*cos(angle_per)-ilox_yyc(jj,j)*sin(angle_per)
                    ilox_yyc(jj,j)=tempxx*sin(angle_per)+ilox_yyc(jj,j)*cos(angle_per)
                    !write(3300+n,'(6es14.6,i5)'),ilox_xxc(jj,1),ilox_yyc(jj,1),ilox_zzc(jj,1),ilox_xxc(jj,j),ilox_yyc(jj,j),ilox_zzc(jj,j),ilox_periodicflag(jj,j)
                    do kk=1,8
                    tempxx=ilon_x(jj,j,kk)
                            ilon_x(jj,j,kk)=tempxx*cos(angle_per)-ilon_y(jj,j,kk)*sin(angle_per)
                            ilon_y(jj,j,kk)=tempxx*sin(angle_per)+ilon_y(jj,j,kk)*cos(angle_per)
                    end do
            else if (ilox_periodicflag(jj,j).eq.1) then
                if(abs(ilox_xxc(jj,j)-ilox_xxc(jj,1)).gt.xper*oo2)then
                    ilox_xxc(jj,j)=ilox_xxc(jj,j)+(xper*sign(1.0,ilox_xxc(jj,1)-xper*oo2))
                    do kk=1,8
                        ilon_x(jj,j,kk)=ilon_x(jj,j,kk)+(xper*sign(1.0,ilox_xxc(jj,1)-xper*oo2))
                    end do
                end if
                if(abs(ilox_yyc(jj,j)-ilox_yyc(jj,1)).gt.yper*oo2)then
                    ilox_yyc(jj,j)=ilox_yyc(jj,j)+(yper*sign(1.0,ilox_yyc(jj,1)-yper*oo2))
                    do kk=1,8
                        ilon_y(jj,j,kk)=ilon_y(jj,j,kk)+(yper*sign(1.0,ilox_yyc(jj,1)-yper*oo2))
                    end do
                end if
                if(abs(ilox_zzc(jj,j)-ilox_zzc(jj,1)).gt.zper*oo2)then
                    ilox_zzc(jj,j)=ilox_zzc(jj,j)+(zper*sign(1.0,ilox_zzc(jj,1)-zper*oo2))
                    do kk=1,8
                        ilon_z(jj,j,kk)=ilon_z(jj,j,kk)+(zper*sign(1.0,ilox_zzc(jj,1)-zper*oo2))
                    end do
                end if
            end if
	    end if
	    end do
		
	end do
	end if
	
	  
      vext=0.0d0
      nodes_list=0.0d0
      eltype=ielem_ishape(i)
      elem_dec=ielem_vdec(i)
      elem_listd=0.0d0
      jx=ielem_nonodes(i)
	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,1)=ilon_x(1,1,k)
	    nodes_list(k,2)=ilon_y(1,1,k)
	    nodes_list(k,3)=ilon_z(1,1,k)
	    vext(k,:)=nodes_list(k,:)
	  end do
	  call decompose3(n,eltype,nodes_list,elem_listd)
    
      select case(ielem_ishape(i))
      case(1)
      if (ielem_mode(i).eq.0)then
      call quadraturehexa(n,igqrules,vext,qpoints,wequa3d)
      dumv1=hexavolume(n,vext,qpoints,wequa3d)
      call compute_centre3d(i,cords)
      vext(1,1:dims)=cords(1:dims)
      else
	vext(1:4,1:3)=elem_listd(1,1:4,1:3)
	  call computejacobians(n,vext,vva1,deta)
      end if	    
      case(2)
      vext(1:4,1:3)=elem_listd(1,1:4,1:3)
	call computejacobians(n,vext,vva1,deta)
      case(3)
	vext(1:4,1:3)=elem_listd(1,1:4,1:3)
	  call computejacobians(n,vext,vva1,deta)
      case(4)
       if (ielem_mode(i).eq.0)then
      call quadratureprism(n,igqrules,vext,qpoints,wequa3d)
      dumv1=prismvolume(n,vext,qpoints,wequa3d)
      vext(1,1:3)=vext(6,1:3)
      else
      vext(1:4,1:3)=elem_listd(1,1:4,1:3)
      call computejacobians(n,vext,vva1,deta)
      end if
      end select
		  temp_cg(1:3)=vext(1,1:3)
		  rec_invccjac(1:3,1:3,i)=vva1(1:3,1:3)
		  

		  
		  
		  
		  
		  detjc=deta(1)
	    rec_vext_ref(1:3,i)=temp_cg(1:3)
	    
	    
	    if ((poly.eq.4).or.(dg.eq.1))then
		    rec_invccjac(1:3,1:3,i)=zero
		     rec_invccjac(1,1,i)=1.0d0
		      rec_invccjac(2,2,i)=1.0d0
		       rec_invccjac(3,3,i)=1.0d0
		      detjc=1.0
		      
		      temp_cg(1)=ielem_xxc(i)
		      temp_cg(2)=ielem_yyc(i)
		      temp_cg(3)=ielem_zzc(i)
		      rec_vext_ref(1:3,i)=temp_cg(1:3)
		    end if
	    
	    
	    
	    
	    
	    

		do jj=1,ielem_admis(i)
		
                        if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        end if
		
			do l=1,itarget
				
				
				
				if (fastest.ne.1)then
				
				
				
				do kk=1,8
					tempcentres=zero
					tempcentres(1)=ilon_x(jj,l,kk)
					tempcentres(2)=ilon_y(jj,l,kk)
					tempcentres(3)=ilon_z(jj,l,kk)
					  do mm_j=1,dimensiona
					    mm_old(mm_j)=tempcentres(mm_j)-temp_cg(mm_j)
					  end do
					  do mm_i=1,dimensiona
					    tempcentres(mm_i)=zero
					    do mm_j=1,dimensiona
					      tempcentres(mm_i)=tempcentres(mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
					    end do
					  end do
					ilon_x(jj,l,kk)=tempcentres(1)
					ilon_y(jj,l,kk)=tempcentres(2)
					ilon_z(jj,l,kk)=tempcentres(3)
				end do
				end if
			end do	
		end do	
		idum=0
		if (ielem_interior(i).eq.1)then
                        do j=1,ielem_ifca(i)
                        if (ielem_ibounds(j,i).gt.0)then
                            if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
                                idum=1
                            end if
                        end if
                        end do
                end if
		
		


		

		
		do jj=1,ielem_admis(i)
			if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        end if
		
			do l=1,itarget
				if (rec_local(i).eq.0)then
				ilox_volume(jj,l)=(ielem_totvolume(ilox_ihexl(jj,l)))/abs(detjc)
				
				
				
				
					if ((ees.ne.5).or.(jj.eq.1))then
						if (idum.eq.1)then
						if (itestcase.eq.4)then
						if (greengo.eq.0)rec_volume_w(1,l,rec_wall(i))=ilox_volume(1,l)
						end if
						rec_volume(1,1,i)=ilox_volume(1,1)
						else
						rec_volume(1,1,i)=ilox_volume(1,1)
						end if

					rec_ihexg(jj,l,i)=ilox_ihexg(jj,l)
					rec_ihexl(jj,l,i)=ilox_ihexl(jj,l)
					 if (iperiodicity == 1) then
					rec_periodicflag(jj,l,i)=ilox_periodicflag(jj,l)
					end if
					else

					rec_ihexgc(jj,l,i)=ilox_ihexg(jj,l)
					rec_ihexlc(jj,l,i)=ilox_ihexl(jj,l)

					end if
				
! 				
				else
				if (ilox_ihexb(jj,l).eq.n)then
				ilox_volume(jj,l)=(ielem_totvolume(ilox_ihexl(jj,l)))/abs(detjc)
				if ((ees.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				if (itestcase.eq.4)then
				if (greengo.eq.0)rec_volume_w(1,l,rec_wall(i))=ilox_volume(1,l)
				end if
				rec_volume(1,1,i)=ilox_volume(1,1)
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				end if
				
				
	
				rec_ihexg(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexl(jj,l,i)=ilox_ihexl(jj,l)
				 if (iperiodicity == 1) then
				rec_periodicflag(jj,l,i)=ilox_periodicflag(jj,l)
				end if
				rec_ihexb(jj,l,rec_local(i))=ilox_ihexb(jj,l)
 				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				rec_ihexgc(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexlc(jj,l,i)=ilox_ihexl(jj,l)
				rec_ihexbc(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				

                                end if

				else 
				if ((ees.ne.5).or.(jj.eq.1))then
				rec_ihexg(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexl(jj,l,i)=ilox_ihexl(jj,l)
				 if (iperiodicity == 1) then
				rec_periodicflag(jj,l,i)=ilox_periodicflag(jj,l)
				end if
				rec_ihexb(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				rec_ihexn(jj,l,rec_local(i))=ilox_ihexn(jj,l)
				else
				rec_ihexgc(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexlc(jj,l,i)=ilox_ihexl(jj,l)
				rec_ihexbc(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				rec_ihexnc(jj,l,rec_local(i))=ilox_ihexn(jj,l)
				
				end if
				
				
				
				eltype=ilox_ishape(jj,l)
				elem_listd=0.0d0; vext=0.0d0; nodes_list=0.0d0
				      
				      select case(eltype)
				      
				      case(1)
				      elem_dec=6; jx=8
				      case(2)
				    elem_dec=1; jx=4
				      case(3)
					elem_dec=2; jx=5
				      case(4)
				      elem_dec=3; jx=6
				      end select

					    do k=1,jx
					      nodes_list(k,1)=ilon_x(jj,l,k)
					      nodes_list(k,2)=ilon_y(jj,l,k)
					      nodes_list(k,3)=ilon_z(jj,l,k)
					      vext(k,:)=nodes_list(k,:)
					      
					    end do
					    call decompose3(n,eltype,nodes_list,elem_listd)
					    
					    dumv2=0.0d0
					    do k=1,elem_dec
					    vext(1:4,1:3)=elem_listd(k,1:4,1:3)
					    dumv2=dumv2+tetravolume(n,vext)
					    
					    end do
					    
					    	
					
					    
								    
					    
					    
					    ilox_volume(jj,l)=dumv2
				

				if ((ees.ne.5).or.(jj.eq.1))then
				if (idum.eq.1)then
				if (itestcase.eq.4)then
				if (greengo.eq.0)rec_volume_w(1,l,rec_wall(i))=ilox_volume(1,l)
				end if
				rec_volume(1,1,i)=ilox_volume(1,1)
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				end if
				
				!rec_volume(jj,l,i)=ilox_volume(jj,l)
				
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				!rec_volumec(jj,l,i)=ilox_volume(jj,l)
				end if
				
				

				end if
				end if
! 							
			end do	
		end do

		if (ielem_interior(i).eq.0)then
	  call compute_centre3d(i,cords)
	  vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem_ifca(i)
		  j=ielem_ineigh(k,i)
		  call compute_centre3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1
  		    ielem_dih2(k,1:dims,i)=vext(2,1:dims)-vext(1,1:dims)
		    end if
	  end do
      else
		    call compute_centre3d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem_ifca(i)
		if (ielem_ineighg(k,i).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		  select case(ielem_types_faces(k,i))
		  case(5)
		  ixxfff=4
		  case(6)
		  ixxfff=3
		  end select
		  call compute_centre3df(n,i,k,ixxfff,cords)
		  vext(2,1:dims)=cords(1:dims)
	
		  
		      dist1=distance3(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1*2.0d0
 			ielem_dih2(k,1:dims,i)=vext(2,1:dims)-vext(1,1:dims)
		    end if
		 end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).eq.0))then	!non periodic boundaries 
		if (ielem_ineighb(k,i).eq.n)then		!within my cpu
		 j=ielem_ineigh(k,i)
		  call compute_centre3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1
  		    ielem_dih2(k,1:dims,i)=vext(2,1:dims)-vext(1,1:dims)
		    end if
		else						!from another cpu 
		    do in1=1,ielem_inumneighbours(i)
			  if (ielem_ineighg(k,i).eq.rec_ihexg(1,in1,i))then
				  ielem_indexi(k,i)=in1
				      if (rungekutta.ge.2)then
		    vext(2,1)=ilox_xxc(1,in1);vext(2,2)=ilox_yyc(1,in1); vext(2,3)=ilox_zzc(1,in1)
		     dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
 			ielem_dih2(k,1:dims,i)=vext(2,1:dims)-vext(1,1:dims)
				      end if
			  end if
		    end do
		end if
		end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).gt.0))then	!periodic boundaries within my cpu
		if (ielem_ineighb(k,i).eq.n)then	
		     j=ielem_ineigh(k,i)
		  call compute_centre3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    if(per_rot.eq.0)then
		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
		    end if
		    if(abs(vext(2,3)-vext(1,3)).gt.zper*oo2)then
		    vext(2,3)=vext(2,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
		    end if
		    else
                if (ibound_icode(ielem_ibounds(k,i)).eq.50) then
                    do kk=1,n_node
                    tempxx=vext(kk,1)
                    vext(kk,1)=tempxx*cos(angle_per)-sin(angle_per)*vext(kk,2)
                    vext(kk,2)=tempxx*sin(angle_per)+cos(angle_per)*vext(kk,2)
                    end do
                else
                    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
                    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
                    end if
                    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
                    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
                    end if
                    if(abs(vext(2,3)-vext(1,3)).gt.zper*oo2)then
                    vext(2,3)=vext(2,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
                    end if
                end if
		    end if
		    dist1=distance3(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1
  		    ielem_dih2(k,1:dims,i)=vext(2,1:dims)-vext(1,1:dims)
		    end if
		else	!periodic boundaries from another cpu
		     do in1=1,ielem_inumneighbours(i)
			  if (ielem_ineighg(k,i).eq.rec_ihexg(1,in1,i))then
				  ielem_indexi(k,i)=in1
				      if (rungekutta.ge.2)then
		    vext(2,1)=ilox_xxc(1,in1);vext(2,2)=ilox_yyc(1,in1); vext(2,3)=ilox_zzc(1,in1)
		    if(per_rot.eq.0)then 
		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
		    end if
		    if(abs(vext(2,3)-vext(1,3)).gt.zper*oo2)then
		    vext(2,3)=vext(2,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
		    end if
		    else
                if (ibound_icode(ielem_ibounds(k,i)).eq.50) then
                    do kk=1,n_node
				      tempxx=vext(kk,1)
				      vext(kk,1)=tempxx*cos(angle_per)-sin(angle_per)*vext(kk,2)
				      vext(kk,2)=tempxx*sin(angle_per)+cos(angle_per)*vext(kk,2)
				    end do  
                else
                    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
                    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
                    end if
                    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
                    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
                    end if
                    if(abs(vext(2,3)-vext(1,3)).gt.zper*oo2)then
                    vext(2,3)=vext(2,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
                    end if
                end if
		    end if
		    dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
  		    ielem_dih2(k,1:dims,i)=vext(2,1:dims)-vext(1,1:dims)
				end if
			  end if
		    end do

		end if
		end if
	  end do
	end if


        do jj=1,ielem_admis(i)
                        if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        
                        end if
			do l=1,itarget
				if (fastest.ne.1)then
				tempcentres=zero
				tempcentres(1)=ilox_xxc(jj,l)
				tempcentres(2)=ilox_yyc(jj,l)
				tempcentres(3)=ilox_zzc(jj,l)
				  do mm_j=1,dimensiona
				    mm_old(mm_j)=tempcentres(mm_j)-temp_cg(mm_j)
				  end do
				  do mm_i=1,dimensiona
				    tempcentres(mm_i)=zero
				    do mm_j=1,dimensiona
				      tempcentres(mm_i)=tempcentres(mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
				    end do
				  end do
				ilox_xxc(jj,l)=tempcentres(1)
				ilox_yyc(jj,l)=tempcentres(2)
				ilox_zzc(jj,l)=tempcentres(3)

				
				
				
				
				
				
				
				
				
                        
				end if
			end do	
		end do					

	

		
	
! 
end subroutine localise_sten2



subroutine  localise_sten2d(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
implicit none
integer :: mm_i
integer :: mm_j
real :: mm_old(1:dimensiona)
integer::i,j,k,l,kk,prk,jj,kmaxe,ineedt,jx2,jx,in1,facexx,ixxfff,ihgt,ihgj,itarget,idum,in_sten,nj,elem_dec,eltype
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexg  !global index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexl  !local index of cells
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
integer,allocatable,dimension(:,:),intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
integer,allocatable,dimension(:,:),intent(inout)::ilox_ishape !shape of each element
real,allocatable,dimension(:,:),intent(inout)::ilox_xxc       !cell centre coordinates in x
real,allocatable,dimension(:,:),intent(inout)::ilox_yyc       !cell centre coordinates in y
real,allocatable,dimension(:,:),intent(inout)::ilox_zzc      !cell centre coordinates in z
real,allocatable,dimension(:,:),intent(inout)::ilox_volume    !cell volume
integer,allocatable,dimension(:,:),intent(inout)::ilox_periodicflag
integer,allocatable,dimension(:,:,:),intent(inout)::ilon_nodcount  !number of nodes
real,allocatable,dimension(:,:,:),intent(inout)::ilon_x           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_y           !coordinates of each node in x
real,allocatable,dimension(:,:,:),intent(inout)::ilon_z           !coordinates of each node in x
real,dimension(2)::tempcentres,rel2
real::dumv1,dumv2,detjc,dist1,distfd,x1x,y1y,x2x,y2y
real,dimension(1:dimensiona)::cords
integer,intent(in)::iconsi,n
integer,dimension(4)::nojcount
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real,dimension(1:dimensiona,1:dimensiona)::vva1
real,dimension(1)::deta

kmaxe=xmpielrank(n)
ineedt=direcexr(1)%tot
i=iconsi
	if (iperiodicity.eq.1)then

	  do jj=1,ielem_admis(i)
	  if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        
                        end if
	    do j=2,itarget
	    if(abs(ilox_xxc(jj,j)-ilox_xxc(jj,1)).gt.xper*oo2)then
		    ilox_xxc(jj,j)=ilox_xxc(jj,j)+(xper*sign(1.0,ilox_xxc(jj,1)-xper*oo2))
		    do kk=1,4
		    ilon_x(jj,j,kk)=ilon_x(jj,j,kk)+(xper*sign(1.0,ilox_xxc(jj,1)-xper*oo2))
		    end do
	    end if
	    if(abs(ilox_yyc(jj,j)-ilox_yyc(jj,1)).gt.yper*oo2)then
		    ilox_yyc(jj,j)=ilox_yyc(jj,j)+(yper*sign(1.0,ilox_yyc(jj,1)-yper*oo2))
		  
		    do kk=1,4
		    ilon_y(jj,j,kk)=ilon_y(jj,j,kk)+(yper*sign(1.0,ilox_yyc(jj,1)-yper*oo2))
		    end do
		    
	    end if
	    end do
		
	end do
	end if
	

	!$omp master

	if (code_profile.eq.30)then
		do nj=1,ielem_nonodes(i)
				x1x=ilon_x(1,1,nj)
				y1y=ilon_y(1,1,nj)
				nojcount(nj)=0

				if (rec_local(i).eq.0)then
					do l=2,itarget
						j=(xmpil(ilox_ihexg(1,l)))
						do k=1,ielem_nonodes(j)    
							x2x=ilon_x(1,l,k)
							y2y=ilon_y(1,l,k)

							distfd=sqrt(((x1x-x2x)**2)+((y1y-y2y)**2))

							if (distfd.lt.tolsmall)then
								nojcount(nj)=nojcount(nj)+1
								ielem_nodes_neighbours(nj,nojcount(nj),i)=l
							end if
						end do
					end do
				end if
				!---------------------- mixed----------!
				if (rec_local(i).ne.0)then
					do l=2,itarget
						if (ilox_ihexb(1,l).eq.n)then
							j=(xmpil(ilox_ihexg(1,l)))
							do k=1,ielem_nonodes(j)    
								x2x=ilon_x(1,l,k)
								y2y=ilon_y(1,l,k)

								distfd=sqrt(((x1x-x2x)**2)+((y1y-y2y)**2))

								if (distfd.lt.tolsmall)then
									nojcount(nj)=nojcount(nj)+1
									ielem_nodes_neighbours(nj,nojcount(nj),i)=l
								end if
							end do
						else
							select case(ilox_ishape(1,l))
							case(5)
							in_sten=4
							case(6)
							in_sten=3
							end select
							do k=1,in_sten
								x2x=ilon_x(1,l,k)
								y2y=ilon_y(1,l,k)

								distfd=sqrt(((x1x-x2x)**2)+((y1y-y2y)**2))

								if (distfd.lt.tolsmall)then
									nojcount(nj)=nojcount(nj)+1
									ielem_nodes_neighbours(nj,nojcount(nj),i)=l
								end if

							end do
						end if
					end do
				end if

			end do
			write(630+n,*)"element number global",ielem_ihexgl(i)

			
			do nj=1,ielem_nonodes(i)
				write(630+n,*)"node number",nj
				ielem_nojecount(nj,i)=nojcount(nj)

				do j=1,nojcount(nj)
				!if (ielem_nodes_neighbours(nj,j,i).gt.0)then
					write(630+n,*)j,ielem_nodes_neighbours(nj,j,i)
				!end if
				end do
			end do


	end if
	!$omp end master




	


	  
     vext=0.0d0
    nodes_list=0.0d0
    eltype=ielem_ishape(i)
    elem_dec=ielem_vdec(i)
    elem_listd=0.0d0
       
      jx=ielem_nonodes(i)
	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,1)=ilon_x(1,1,k)
	    nodes_list(k,2)=ilon_y(1,1,k)
	   
	    vext(k,:)=nodes_list(k,:)
	  end do
	  call decompose2(n,eltype,nodes_list,elem_listd)
    
      select case(ielem_ishape(i))

      case(5)
      if (ielem_mode(i).eq.0)then
      call quadraturequad(n,igqrules,vext,qpoints,wequa3d)
      dumv1=quadvolume(n,vext,qpoints,wequa3d)
      call compute_centre2d(i,cords)
      vext(1,1:dims)=cords(1:dims)
      else
	vext(1:3,1:2)=elem_listd(1,1:3,1:2)
! 	  dumv1=trianglevolume(n)
	  call computejacobians2(n,vext,vva1,deta)
	end if
      

      case(6)
      vext(1:3,1:2)=elem_listd(1,1:3,1:2)
! 	dumv1=trianglevolume(n)
	call computejacobians2(n,vext,vva1,deta)

      end select
      
		
		  rec_invccjac(1:2,1:2,i)=vva1(1:2,1:2)
		  
		  
		  
		  
		  
		  
		  
		 
		    
 		idum=0
		if (ielem_interior(i).eq.1)then
                        do j=1,ielem_ifca(i)
                        if (ielem_ibounds(j,i).gt.0)then
                            if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
                                idum=1
                            end if
                        end if
                        end do
                end if
		   
		   
		   
		   
		  
        
		  
		  
		  
            detjc=deta(1)
		    rec_vext_ref(1:2,i)=vext(1,1:2)
		    
		    if ((poly.eq.4).or.(dg.eq.1))then
		    rec_invccjac(1:2,1:2,i)=0.0d0
 		     rec_invccjac(1,1,i)=1.0d0
 		      rec_invccjac(2,2,i)=1.0d0
		      detjc=1.0
		      vext(1,1)=ielem_xxc(i)
		      vext(1,2)=ielem_yyc(i)
		      rec_vext_ref(1:2,i)=vext(1,1:2)
		    end if
		    
		    
		  do jj=1,ielem_admis(i)
                            if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        
                        end if
		  
		  
			do l=1,itarget
! 				tempcentres=zero
! 				tempcentres(1)=ilox_xxc(jj,l)
! 				tempcentres(2)=ilox_yyc(jj,l)
! 				tempcentres(:)=explicit_matrix_product(rec_invccjac(:,:,i),tempcentres(:)-vext(1,:))
! 				
! 				ilox_xxc(jj,l)=tempcentres(1)
! 				ilox_yyc(jj,l)=tempcentres(2)
				
				
				
				
				do kk=1,4
					tempcentres=zero
					tempcentres(1)=ilon_x(jj,l,kk)
					tempcentres(2)=ilon_y(jj,l,kk)
					
					
					  do mm_j=1,dimensiona
					    mm_old(mm_j)=tempcentres(mm_j)-vext(1,mm_j)
					  end do
					  do mm_i=1,dimensiona
					    tempcentres(mm_i)=zero
					    do mm_j=1,dimensiona
					      tempcentres(mm_i)=tempcentres(mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
					    end do
					  end do
					
					ilon_x(jj,l,kk)=tempcentres(1)
					ilon_y(jj,l,kk)=tempcentres(2)
					
				end do
				
			end do	
		end do			

					
		do jj=1,ielem_admis(i)
                        if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        
                        end if
			do l=1,itarget
				if (rec_local(i).eq.0)then
				ilox_volume(jj,l)=(ielem_totvolume(ilox_ihexl(jj,l)))/abs(detjc)
				 if ((ees.ne.5).or.(jj.eq.1))then
				 if (idum.eq.1)then
				 if (itestcase.eq.4)then
				if (greengo.eq.0)rec_volume_w(1,l,rec_wall(i))=ilox_volume(1,l)
				end if
				rec_volume(1,1,i)=ilox_volume(1,1)
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				end if
				
				rec_ihexg(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexl(jj,l,i)=ilox_ihexl(jj,l)
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				rec_ihexgc(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexlc(jj,l,i)=ilox_ihexl(jj,l)
				
				
				end if
				else
				if (ilox_ihexb(jj,l).eq.n)then
				ilox_volume(jj,l)=(ielem_totvolume(ilox_ihexl(jj,l)))/abs(detjc)
				 if ((ees.ne.5).or.(jj.eq.1))then
						if (idum.eq.1)then
								if (itestcase.eq.4)then
								if (greengo.eq.0)rec_volume_w(1,l,rec_wall(i))=ilox_volume(1,l)
								end if
						rec_volume(1,1,i)=ilox_volume(1,1)
						else
						rec_volume(1,1,i)=ilox_volume(1,1)
						end if
				rec_ihexg(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexl(jj,l,i)=ilox_ihexl(jj,l)
				rec_ihexb(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				rec_ihexgc(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexlc(jj,l,i)=ilox_ihexl(jj,l)
				rec_ihexbc(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				end if
				else 
				 if ((ees.ne.5).or.(jj.eq.1))then
				rec_ihexg(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexl(jj,l,i)=ilox_ihexl(jj,l)
				rec_ihexb(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				rec_ihexn(jj,l,rec_local(i))=ilox_ihexn(jj,l)
				else
				rec_ihexgc(jj,l,i)=ilox_ihexg(jj,l)
				rec_ihexlc(jj,l,i)=ilox_ihexl(jj,l)
				rec_ihexbc(jj,l,rec_local(i))=ilox_ihexb(jj,l)
				rec_ihexnc(jj,l,rec_local(i))=ilox_ihexn(jj,l)
				
				
				end if
				
				      eltype=ilox_ishape(jj,l)
				      elem_listd=0.0d0; vext=0.0d0; nodes_list=0.0d0
				      select case(eltype)
				      
				      case(5)
				      elem_dec=2; jx=4
				      case(6)
				    elem_dec=1; jx=3
				      
				      end select

					    do k=1,jx
					      nodes_list(k,1)=ilon_x(jj,l,k)
					      nodes_list(k,2)=ilon_y(jj,l,k)
					     
					      vext(k,:)=nodes_list(k,:)
					    end do
					    call decompose2(n,eltype,nodes_list,elem_listd)
					    dumv2=0.0d0
					    do k=1,elem_dec
					    vext(1:3,1:2)=elem_listd(k,1:3,1:2)
					    dumv2=dumv2+trianglevolume(n,vext)
					    end do
					    ilox_volume(jj,l)=dumv2
					     if ((ees.ne.5).or.(jj.eq.1))then
						 if (idum.eq.1)then
						 if (itestcase.eq.4)then
				if (greengo.eq.0)rec_volume_w(1,l,rec_wall(i))=ilox_volume(1,l)
				end if
				rec_volume(1,1,i)=ilox_volume(1,1)
				else
				rec_volume(1,1,i)=ilox_volume(1,1)
				end if
				
                                            else
                                  rec_volume(1,1,i)=ilox_volume(1,1)
                                            
                                            end if
! 				
				end if
				end if
			end do	
		end do

	      if (ielem_interior(i).eq.0)then
	  call compute_centre2d(i,cords)
	  vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem_ifca(i)
		  j=ielem_ineigh(k,i)
		  call compute_centre2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1
		     ielem_dih2(k,1:dimensiona,i)=vext(2,1:dimensiona)-vext(1,1:dimensiona)
		    end if
	  end do
      else
		    call compute_centre2d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem_ifca(i)
		if (ielem_ineighg(k,i).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		  
		  ixxfff=2
		  
		  
		  call compute_centre2df(n,iconsi,facexx,ixxfff,cords)
		  vext(2,1:dims)=cords(1:dims)
		  
		      dist1=distance2(n,vext)

		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1*2.0d0
		    ielem_dih2(k,1:dimensiona,i)=vext(2,1:dimensiona)-vext(1,1:dimensiona)
		    end if
		 end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).eq.0))then	!non periodic boundaries 
		if (ielem_ineighb(k,i).eq.n)then		!within my cpu
		 j=ielem_ineigh(k,i)
		  call compute_centre2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1
		     ielem_dih2(k,1:dimensiona,i)=vext(2,1:dimensiona)-vext(1,1:dimensiona)
		    end if
		else						!from another cpu 
		    do in1=1,ielem_inumneighbours(i)
			  if (ielem_ineighg(k,i).eq.rec_ihexg(1,in1,i))then
				  ielem_indexi(k,i)=in1
				      if (rungekutta.ge.2)then
		    vext(2,1)=ilox_xxc(1,in1);vext(2,2)=ilox_yyc(1,in1)
		     dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
		    ielem_dih2(k,1:dimensiona,i)=vext(2,1:dimensiona)-vext(1,1:dimensiona)
				      end if
			  end if
		    end do
		end if
		end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).gt.0))then	!periodic boundaries within my cpu
		if (ielem_ineighb(k,i).eq.n)then	
		     j=ielem_ineigh(k,i)
		  call compute_centre2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
		    end if
		    
		    dist1=distance2(n,vext)
		    if (rungekutta.ge.2)then
		    ielem_dih(k,i)=dist1
		     ielem_dih2(k,1:dimensiona,i)=vext(2,1:dimensiona)-vext(1,1:dimensiona)
		    end if
		else	!periodic boundaries from another cpu
		     do in1=1,ielem_inumneighbours(i)
			  if (ielem_ineighg(k,i).eq.rec_ihexg(1,in1,i))then
				  ielem_indexi(k,i)=in1
				      if (rungekutta.ge.2)then
		    vext(2,1)=ilox_xxc(1,in1);vext(2,2)=ilox_yyc(1,in1);
		     
		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
		    end if
		    
		    dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
		    ielem_dih2(k,1:dimensiona,i)=vext(2,1:dimensiona)-vext(1,1:dimensiona)
				end if
			  end if
		    end do

		end if
		end if
	  end do
	end if

        do jj=1,ielem_admis(i)
			if ((ees.ne.5).or.(jj.eq.1))then
                        itarget=ielem_inumneighbours(i)
                        
                        else
                        itarget=numneighbours2
                        
                        
                        end if
			do l=1,itarget
				tempcentres=zero
				tempcentres(1)=ilox_xxc(jj,l)
				tempcentres(2)=ilox_yyc(jj,l)
				  do mm_j=1,dimensiona
				    mm_old(mm_j)=tempcentres(mm_j)-vext(1,mm_j)
				  end do
				  do mm_i=1,dimensiona
				    tempcentres(mm_i)=zero
				    do mm_j=1,dimensiona
				      tempcentres(mm_i)=tempcentres(mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
				    end do
				  end do
				
				ilox_xxc(jj,l)=tempcentres(1)
				ilox_yyc(jj,l)=tempcentres(2)
				
				
				
				
				
				
			end do	
		end do


		
	
! 
end subroutine localise_sten2d










subroutine direct_side(n)
!> @brief
!> this subroutine establishes the distance betwen cell centres for each face
implicit none
integer,intent(in)::n
integer::i,j,k,kmaxe,facexx,ixxfff
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dimensiona)::cords
real::dist1

kmaxe=xmpielrank(n)

!$omp do
do i=1,kmaxe
	if (ielem_interior(i).eq.0)then
		    call compute_centre3d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
      
	  do k=1,ielem_ifca(i)
		  j=ielem_ineigh(k,i)
		  call compute_centre3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
	  end do
	
	else
		    call compute_centre3d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem_ifca(i)
		if (ielem_ineighg(k,i).eq.0)then	!boundaries except other cpus and periodics
		  
		  facexx=k
		  select case(ielem_types_faces(k,i))
		  case(5)
		  ixxfff=4
		  case(6)
		  ixxfff=3
		  end select
		  call compute_centre3df(n,i,k,ixxfff,cords)
		  vext(2,1:dims)=cords(1:dims)

		  
		      dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1*2.0d0
		 end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).eq.0))then	!non periodic boundaries 
		if (ielem_ineighb(k,i).eq.n)then		!within my cpu
		 j=ielem_ineigh(k,i)
		  call compute_centre3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
		else						!from another cpu 
		
! 		    vext(2,1:dims)=dsolchanger(ielem_ineighn(k,i))%centres(ielem_q_face_q_mapl(1,k,i),1:dims)
		    vext(2,1:dims)=dsolchanger(ielem_ineighn(k,i))%centres(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:dims)


		     dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
		end if
		end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).gt.0))then	!periodic boundaries within my cpu
		if (ielem_ineighb(k,i).eq.n)then	
		     j=ielem_ineigh(k,i)
		  call compute_centre3d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
		    end if
		    if(abs(vext(2,3)-vext(1,3)).gt.zper*oo2)then
		    vext(2,3)=vext(2,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
		    end if
		    dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
		
		else	!periodic boundaries from another cpu

! 		     vext(2,1:dims)=dsolchanger(ielem_ineighn(k,i))%centres(ielem_q_face_q_mapl(1,k,i),1:dims)
		     vext(2,1:dims)=dsolchanger(ielem_ineighn(k,i))%centres(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:dims)



		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
		    end if
		    if(abs(vext(2,3)-vext(1,3)).gt.zper*oo2)then
		    vext(2,3)=vext(2,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
		    end if
		    dist1=distance3(n,vext)
		    ielem_dih(k,i)=dist1
		end if
		end if
	  end do
	end if

end do
!$omp end do



end subroutine direct_side



subroutine direct_side2d(n)
!> @brief
!> this subroutine establishes the distance betwen cell centres for each edge
implicit none
integer,intent(in)::n
integer::i,j,k,kmaxe,facexx,ixxfff
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dimensiona)::cords
real::dist1

kmaxe=xmpielrank(n)


!$omp do
do i=1,kmaxe
	if (ielem_interior(i).eq.0)then
		    call compute_centre2d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
      
	  do k=1,ielem_ifca(i)
		  j=ielem_ineigh(k,i)
		  call compute_centre2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
	  end do
	
	else
		    call compute_centre2d(i,cords)
		    vext(1,1:dims)=cords(1:dims)
	  do k=1,ielem_ifca(i)
		if (ielem_ineighg(k,i).eq.0)then	!boundaries except other cpus and periodics
		  facexx=k
		 
		  ixxfff=2
		 
		  call compute_centre2df(n,i,facexx,ixxfff,cords)
		  vext(2,1:dims)=cords(1:dims)
		  
		      dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1*2.0d0
		 end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).eq.0))then	!non periodic boundaries 
		if (ielem_ineighb(k,i).eq.n)then		!within my cpu
		 j=ielem_ineigh(k,i)
		  call compute_centre2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)
		      dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
		else						!from another cpu 
		    vext(2,1:dims)=dsolchanger(ielem_ineighn(k,i))%centres(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:dims)

		     dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
		end if
		end if
		if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).gt.0))then	!periodic boundaries within my cpu
		if (ielem_ineighb(k,i).eq.n)then	
		     j=ielem_ineigh(k,i)
		  call compute_centre2d(j,cords)
		    vext(2,1:dims)=cords(1:dims)  
		    if(abs(vext(2,1)-vext(1,1)).gt.xper/2.d0)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0d0,vext(1,1)-xper/2.d0))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper/2.d0)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0d0,vext(1,2)-yper/2.d0))
		    end if
		    
		    dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
		
		else	!periodic boundaries from another cpu

		     vext(2,1:dims)=dsolchanger(ielem_ineighn(k,i))%centres(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:dims)
		    if(abs(vext(2,1)-vext(1,1)).gt.xper*oo2)then
		    vext(2,1)=vext(2,1)+(xper*sign(1.0d0,vext(1,1)-xper/2.0d0))
		    end if
		    if(abs(vext(2,2)-vext(1,2)).gt.yper*oo2)then
		    vext(2,2)=vext(2,2)+(yper*sign(1.0d0,vext(1,2)-yper/2.0d0))
		    end if
		    
		    dist1=distance2(n,vext)
		    ielem_dih(k,i)=dist1
		end if
		end if
	  end do
	end if

end do
!$omp end do



end subroutine direct_side2d


subroutine grads_assign(n)
integer,intent(in)::n
integer::i,j,k,l,jj,kmaxe

kmaxe=xmpielrank(n)

if (dimensiona.eq.3)then

!$omp do
do i=1,kmaxe
	call checkgrads(n,i)
end do
!$omp end do

else

!$omp do
do i=1,kmaxe
	call checkgrads2d(n,i)
end do
!$omp end do

end if




end subroutine grads_assign


subroutine checkgrads(n,iconsi)
!> @brief
!> this subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
implicit none
integer,intent(in)::n,iconsi
real::dxx1,dxx2,tempg1,dist1,dist2,oo2,surfmin,surfmax
integer::i,j,k,l,jj,icount3,nnd,ixf4,idc,idc2
real,dimension(1:dimensiona)::cords
real,dimension(1:8,1:dimensiona)::vext,nodes_list

i=iconsi
    ielem_ggs(i)=greengo
    dxx1=-tolbig; dxx2=tolbig
    call compute_centre3d(i,cords)
    vext(1,1:dims)=cords(1:dims)
      
i=iconsi
  if (ielem_interior(i).eq.1)then
			do k=1,ielem_ifca(i)
								if (ielem_types_faces(k,i).eq.5)then
								ixf4=4
								else
								ixf4=3
								end if

							if (ielem_ineighg(k,i).gt.0)then
										if (ielem_ineighb(k,i).ne.n)then
											if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then

											ielem_nodes_faces(k,1:ixf4,i)=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:ixf4)
											ielem_reorient(k,i)=1

											end if


										else
											if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then

											ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
											ielem_reorient(k,i)=1
									!
											end if
										end if
							end if

			end do
  else
      do k=1,ielem_ifca(i)
		  if (ielem_types_faces(k,i).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if

      if (ielem_ineighg(k,i).gt.0)then
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
! 		
		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
		ielem_reorient(k,i)=1
! 		
	    end if
	
      end if

      
      end do
  end if
	
























    
surfmin=1.0e16
surfmax=1.0e-16



       
      do l=1,ielem_ifca(i)
            
	      select case (ielem_types_faces(l,i))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select

		
				if (ielem_interior(i).eq.1)then
							if ((ielem_ineighg(l,i).gt.0).and.(ielem_ibounds(l,i).gt.0)) then
				!

									end if
							else
		!
					call compute_centre3df(n,i,l,nnd,cords)
				vext(2,1:dims)=cords(1:dims)

				end if
		
		
		
	      surfmin=min(surfmin,ielem_surf(l,i))
	      surfmax=max(surfmax,ielem_surf(l,i))
	      
	      
	      
	      
		    
	


		  dist1=distance3(n,vext)
		if (dist1.lt.dxx2)then
		  dxx2=dist1
		end if
		if (dist1.gt.dxx1)then
		  dxx1=dist1
		end if
end do
        ielem_condition(i)=1.0


        ielem_condition(i)=surfmax/surfmin

		if (code_profile.eq.9)then
        if ((ielem_condition(i).gt.30).or.(ielem_ishape(i).eq.4))then
        ielem_full(i)=0
		end if
		end if

        !if (turbulence.gt.0)then
        

        !if (ielem_ishape(i).eq.2)then
        !ielem_condition(i)=surfmax/surfmin
        
        !if (ielem_condition(i).gt.30)then
        ! ielem_hybrid(i)=1
        !end if
        !end if
        
        !if (ielem_ishape(i).eq.3)then
        !ielem_condition(i)=surfmax/surfmin
        
        !if (ielem_condition(i).gt.10)then
         !ielem_hybrid(i)=1
        !end if
        
       
        
        !end if


        !end if
       



	    tempg1=ielem_condition(i)!max((dxx1/dxx2),(dxx2/dxx1))
! 	   

	    

              if (code_profile.eq.88)then
                      if ((ielem_ishape(i).eq.3)) then
                                        ielem_full(i)=0
                         end if

              end if


	      if (tempg1.gt.gridar1)then
! 	      ielem_ggs(i)=1
	      if ((iadapt.eq.1).or.(code_profile.eq.88).or.(code_profile.eq.98))then
                ielem_full(i)=0
			end if
	      end if
   if (fastest.eq.0)then

! 	      dxx1=-tolbig; dxx2=tolbig
! 	       jj=1
! 	       do l=1,ielem_inumneighbours(i)
! 		       if (rec_volume(jj,l,i).lt.dxx2)then
!
! 		       dxx2=rec_volume(jj,l,i)
! 		       end if
! 		       if (rec_volume(jj,l,i).gt.dxx1)then
!
! 		       dxx1=rec_volume(jj,l,i)
! 		       end if
! 	       end do
! 	 tempg1=max((dxx1/dxx2),(dxx2/dxx1))
	 !ielem_walldist(i)=tempg1
! 	     if (tempg1.gt.gridar2)then
! 	       ielem_ggs(i)=1
! 	  	if ((iadapt.eq.1).or.(code_profile.eq.88).or.(code_profile.eq.98))then
!                 ielem_full(i)=0
! 			end if
! 	       end if

! !
 end if

idc=0
idc2=0
if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if (ibound_icode(ielem_ibounds(j,i)).eq.4)then
	        idc=idc+1
	      else
idc2=idc2+1

		end if
	  end if
        end do
end if

! if ((idc.gt.1))then     !until ge is fully adaptive
! ielem_ggs(i)=1
!
! end if

 if (code_profile.eq.888)then
 if (ielem_interior(i).eq.1)then
 	do j=1,ielem_ifca(i)
 	  if (ielem_ibounds(j,i).gt.0)then
! !                   if (ibound_icode(ielem_ibounds(j,i)).ne.4)then
 	                ielem_hybrid(i)=1
! !                   end if
 	  end if
         end do
 end if
 end if




end subroutine checkgrads






subroutine check3(n,iconsi)
!> @brief
!> this subroutine assigns the ordering for the faces
implicit none
integer,intent(in)::n
integer,intent(inout)::iconsi
real::dxx1,dxx2,tempg1,dist1,dist2,oo2
integer::i,j,k,l,jj,icount3,nnd,ixf4
real,dimension(1:dimensiona)::cords
real,dimension(1:8,1:dimensiona)::vext


do iconsi=1,xmpielrank(n)

i=iconsi
    ielem_ggs(i)=greengo
    dxx1=-tolbig; dxx2=tolbig
    call compute_centre3d(i,cords)
    vext(1,1:dims)=cords(1:dims)
      
i=iconsi
  if (ielem_interior(i).eq.1)then
  do k=1,ielem_ifca(i)
		  if (ielem_types_faces(k,i).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if
	    
      if (ielem_ineighg(k,i).gt.0)then
	if (ielem_ineighb(k,i).ne.n)then
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
		ielem_nodes_faces(k,1:ixf4,i)=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:ixf4)
		ielem_reorient(k,i)=1
	    
	    end if


	else
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
		ielem_reorient(k,i)=1
	    end if
	end if
      end if
!       
  end do
  else
      do k=1,ielem_ifca(i)
		  if (ielem_types_faces(k,i).eq.5)then
	    ixf4=4
	    else
	    ixf4=3
	    end if

      if (ielem_ineighg(k,i).gt.0)then
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
		ielem_reorient(k,i)=1
	    end if
	
      end if

!      

      end do
  end if
end do

end subroutine check3



subroutine checkgrads2d(n,iconsi)
!> @brief
!> this subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
implicit none
integer,intent(in)::n,iconsi
real::dxx1,dxx2,tempg1,dist1,dist2,oo2,surfmin,surfmax
integer::i,j,k,l,jj,icount3,nnd,ixf4
real,dimension(1:dimensiona)::cords
real,dimension(1:8,1:dimensiona)::vext,nodes_list

i=iconsi

tempg1=0.0; 
    ielem_ggs(i)=greengo
    dxx1=-tolbig; dxx2=tolbig
    call compute_centre2d(i,cords)
    vext(1,1:dims)=cords(1:dims)
      
! ielem_inx(i)=0
!   if (ielem_interior(i).eq.1)then
!         do j=1,ielem_ifca(i)
!           if (ielem_ibounds(j,i).gt.0)then
!               if (ibound_icode(ielem_ibounds(j,i)).eq.1)then
!
!               ielem_inx(i)=1
!
!                 end if
!           end if
!         end do
! end if
     

  if (ielem_interior(i).eq.1)then
  do k=1,ielem_ifca(i)
		  
	    ixf4=2
	   
	    
      if (ielem_ineighg(k,i).gt.0)then
	if (ielem_ineighb(k,i).ne.n)then
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
		ielem_nodes_faces(k,1:ixf4,i)=diexboundhirr(ielem_ineighn(k,i))%vertpp(ielem_qface(k,1,ielem_inter_id(ielem_indexf(i))),1:ixf4)
		ielem_reorient(k,i)=1
	    
	    end if


	else
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
		ielem_reorient(k,i)=1
	    end if
	end if
      end if
  end do
  else
      do k=1,ielem_ifca(i)
		ixf4=2

      if (ielem_ineighg(k,i).gt.0)then
	    if (ielem_ineighg(k,i).gt.ielem_ihexgl(i))then
		ielem_nodes_faces(k,1:ixf4,i)=ielem_nodes_faces(ielem_ineighn(k,i),1:ixf4,ielem_ineigh(k,i))
		ielem_reorient(k,i)=1
	    end if
	
      end if
      end do
  end if




  surfmin=1.0e16
surfmax=1.0e-16

     
      do l=1,ielem_ifca(i)
      
	      
	      nnd=2
	     
			 surfmin=min(surfmin,ielem_surf(l,i))
	      surfmax=max(surfmax,ielem_surf(l,i))
		
		  if (ielem_interior(i).eq.1)then
		    if ((ielem_ineighg(l,i).gt.0).and.(ielem_ibounds(l,i).gt.0))then
		      if (ielem_ineighb(l,i).ne.n)then
           
			do k=1,nnd
			  nodes_list(k,1:dims)=dinoder(ielem_nodes_faces(l,k,iconsi))%cord(1:dims)
			end do
			
			
			do k=1,nnd
			if(abs(nodes_list(k,1)-vext(1,1)).gt.xper*oo2)then
			nodes_list(k,1)=nodes_list(k,1)+(xper*sign(1.0d0,vext(1,1)-xper*oo2))
			end if
			if(abs(nodes_list(k,2)-vext(1,2)).gt.yper*oo2)then
			nodes_list(k,2)=nodes_list(k,2)+(yper*sign(1.0d0,vext(1,2)-yper*oo2))
			end if
			end do
			
			call cordinates2(n,nodes_list,nnd,cords(1:2))
			vext(2,1:dims)=cords(1:dims)
		    else
		      call compute_centre2df(n,i,l,nnd,cords)
		      vext(2,1:dims)=cords(1:dims)
		    end if
		  else
		
		    
		      call compute_centre2df(n,i,l,nnd,cords)
		  vext(2,1:dims)=cords(1:dims)
		  end if
		else
		       call compute_centre2df(n,i,l,nnd,cords)
		  vext(2,1:dims)=cords(1:dims)
		
		end if
	      
		    
	


		  dist1=distance2(n,vext)
		if (dist1.lt.dxx2)then
		  dxx2=dist1
		end if
		if (dist1.gt.dxx1)then
		  dxx1=dist1
		end if
		
		
		
end do


			ielem_condition(i)=1.0


        ielem_condition(i)=surfmax/surfmin

        tempg1=ielem_condition(i)


   if ((realgas.eq.1).or.(code_profile.eq.888))then
  if (ielem_interior(i).eq.1)then
  	do j=1,ielem_ifca(i)
  	  if (ielem_ibounds(j,i).gt.0)then
  	        ielem_hybrid(i)=1
  	  end if
          end do
  end if
   end if







! 	    if (tempg1.gt.gridar1)then
! 	      ielem_ggs(i)=1
!
! 	      if ((iadapt.eq.1).or.(code_profile.eq.88).or.(code_profile.eq.98))then
!                 ielem_full(i)=0
! 			end if
!
! 	      end if
!
!
!
!    if (fastest.eq.0)then
! 	       dxx1=-tolbig; dxx2=tolbig
! 	       jj=1
! 	       do l=1,ielem_inumneighbours(i)
! 		        if (rec_volume(jj,l,i).lt.dxx2)then
!
! 		       dxx2=rec_volume(jj,l,i)
! 		       end if
! 		       if (rec_volume(jj,l,i).gt.dxx1)then
!
! 		       dxx1=rec_volume(jj,l,i)
! 		       end if
! 	       end do
! ! 	 tempg1=max((dxx1/dxx2),(dxx2/dxx1))
! ! 	     if (tempg1.gt.gridar2)then
! ! 	       ielem_ggs(i)=1
! ! 	       if ((iadapt.eq.1).or.(code_profile.eq.88).or.(code_profile.eq.98))then
! !                 ielem_full(i)=0
! ! 			end if
! ! 	       end if
! end if


end subroutine checkgrads2d







end module local
