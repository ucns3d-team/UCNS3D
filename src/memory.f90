module memory
use mpiinfo
use declaration
implicit none

contains



subroutine allocate2
!> @brief
!> this subroutine allocates memory
implicit none
allocate(list(2000),ineb(6),iperb(6),nodelist(8))
end subroutine

subroutine allocate5
!> @brief
!> this subroutine allocates memory for the stencils
implicit none
allocate(ilocalallelg(n:n,xmpielrank(n),1,iselemt(n)))
ilocalallelg(:,:,:,:)=0
allocate(ilocalallelgper(n:n,xmpielrank(n),1,iselemt(n)))
ilocalallelgper(:,:,:,:)=0
end subroutine




subroutine allocate3
implicit none
!> @brief
!> this subroutine deallocates memory
deallocate(list,ineb,iperb)
end subroutine

subroutine globaldea2(xmpil,xmpie)
!> @brief
!> this subroutine deallocates global lists
implicit none
integer,allocatable,dimension(:),intent(inout)::xmpil,xmpie
 call mpi_barrier(mpi_comm_world,ierror)
  deallocate(xmpil)
    deallocate(xmpie)
 call mpi_barrier(mpi_comm_world,ierror)
 end subroutine globaldea2

subroutine quadalloc(numberofpoints,numberofpoints2)
!> @brief
!> this subroutine allocates memory for the quadrature points
implicit none
integer,intent(inout)::numberofpoints,numberofpoints2
integer::i,kmaxe


if (dimensiona.eq.3)then

    numberofpoints=max(qp_hexa,qp_tetra,qp_pyra,qp_prism)


    if (dg.eq.1)numberofpoints=max(qp_hexa,qp_tetra*6,qp_pyra,qp_prism)


    numberofpoints2=max(qp_quad,qp_triangle)

    

else
    numberofpoints=max(qp_quad,qp_triangle)
    if (dg.eq.1)numberofpoints=max(qp_quad,qp_triangle*2)
    numberofpoints2=qp_line
   
end if



kmaxe=xmpielrank(n)
do i=1,kmaxe

select case (ielem_ishape(i))
        
        
        case(1) !hexa
        ielem_itotalpoints(i)=qp_tetra*6
        
        case(2) !tetra
        ielem_itotalpoints(i)=qp_tetra
        
        case(3) !pyramid
        ielem_itotalpoints(i)=qp_tetra*2
        
        case(4) !prism
        ielem_itotalpoints(i)=qp_tetra*3
        
        case(5) !quadrilateral
        ielem_itotalpoints(i)=qp_triangle*2
        
        
        
        case(6)!triangle
         ielem_itotalpoints(i)=qp_triangle
         
         
         end select


end do




end subroutine quadalloc








!!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for fluxes!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sumflux_allocation(n)
!> @brief
!> this subroutine allocates memory for the fluxes
	implicit none
	integer,intent(inout)::n
	integer::i,kmaxe
	kmaxe=xmpielrank(n)
	
	allocate (rhs_val(nof_variables,kmaxe))


	
	
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	allocate (rhst_val(turbulenceequations+passivescalar,kmaxe))
	end if
	

	if (dg.eq.1)then
	allocate(rhs_valdg(num_dg_dofs, nof_variables,kmaxe))
    allocate(rhs_sol_mm_dg(1:num_dg_dofs,1:nof_variables,kmaxe))
	end if

	
end subroutine sumflux_allocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine impallocate(n)
!> @brief
!> this subroutine allocates memory for implicit time stepping
implicit none
integer,intent(in)::n
integer::kmaxe,interf
kmaxe=xmpielrank(n)


if (dimensiona.eq.3)then
interf=nof_variables
else
interf=nof_variables
end if

if (rungekutta.eq.12)then
allocate (impdu(kmaxe,1:nof_variables+turbulenceequations+passivescalar))
impdu(:,:)=zero
else



if (relax.eq.3)then

allocate (impdiag_mf(kmaxe))
allocate (impoff_mf(kmaxe,interf))
allocate (impdu(kmaxe,1:nof_variables+turbulenceequations+passivescalar))

if ((itestcase.eq.4).and.((turbulence.gt.0).or.(passivescalar.gt.0)))then
allocate(impdiagt(kmaxe,turbulenceequations+passivescalar))
allocate(impofft(kmaxe,interf,turbulenceequations+passivescalar))
allocate(sht(kmaxe,turbulenceequations+passivescalar))
end if

impdiag_mf=zero
impoff_mf=zero
impdu=zero

else




if (dimensiona.eq.3)then
if (lowmemory.eq.0)then

allocate (impdiag(kmaxe,1:nof_variables,1:nof_variables))
allocate (impoff(kmaxe,6,1:nof_variables,1:nof_variables))
allocate (impdu(kmaxe,1:nof_variables+turbulenceequations+passivescalar))

if (realgas.eq.1)then
allocate(sht_rg(kmaxe,1:nof_variables));sht_rg(:,:)=0.0d0
end if

if ((itestcase.eq.4).and.((turbulence.gt.0).or.(passivescalar.gt.0)))then
allocate(impofft(kmaxe,6,turbulenceequations+passivescalar))
allocate(impdiagt(kmaxe,turbulenceequations+passivescalar))
allocate(sht(kmaxe,turbulenceequations+passivescalar))
end if


impdiag(:,:,:)=zero
impoff(:,:,:,:)=zero
impdu(:,:)=zero

else

allocate (impdu(kmaxe,1:nof_variables+turbulenceequations+passivescalar))
impdu(:,:)=zero

if (realgas.eq.1)then
allocate(sht_rg(kmaxe,1:nof_variables));sht_rg(:,:)=0.0d0
end if

if ((itestcase.eq.4).and.((turbulence.gt.0).or.(passivescalar.gt.0)))then
allocate(sht(kmaxe,turbulenceequations+passivescalar))
end if
end if

else

if (lowmemory.eq.0)then

allocate (impdiag(kmaxe,1:nof_variables,1:nof_variables))
allocate (impoff(kmaxe,4,1:nof_variables,1:nof_variables))
allocate (impdu(kmaxe,1:nof_variables+turbulenceequations+passivescalar))

if (realgas.eq.1)then
allocate(sht_rg(kmaxe,1:nof_variables));sht_rg(:,:)=0.0d0
end if


if ((itestcase.eq.4).and.((turbulence.gt.0).or.(passivescalar.gt.0)))then
allocate(impofft(kmaxe,4,turbulenceequations+passivescalar))
allocate(impdiagt(kmaxe,turbulenceequations+passivescalar))
allocate(sht(kmaxe,turbulenceequations+passivescalar))
end if


impdiag(:,:,:)=zero
impoff(:,:,:,:)=zero
impdu(:,:)=zero

else

! allocate (impdiag(1,1:nof_variables,1:nof_variables))
! allocate (impoff(1,4,1:nof_variables,1:nof_variables))
allocate (impdu(kmaxe,1:nof_variables+turbulenceequations+passivescalar))

if (realgas.eq.1)then
allocate(sht_rg(kmaxe,1:nof_variables));sht_rg(:,:)=0.0d0
end if

if ((itestcase.eq.4).and.((turbulence.gt.0).or.(passivescalar.gt.0)))then
! allocate(impofft(1,4,turbulenceequations+passivescalar))
! allocate(impdiagt(1,turbulenceequations+passivescalar))
allocate(sht(kmaxe,turbulenceequations+passivescalar))
end if
! impdiag(1,:,:)=zero
! impoff(1,:,:,:)=zero
impdu(:,:)=zero
end if



end if
end if
end if




end  subroutine impallocate



subroutine timing(n,cpux1,cpux2,cpux3,cpux4,cpux5,cpux6,timex1,timex2,timex3,timex4,timex5,timex6)
!> @brief
!> this subroutine allocates memory for the timers
implicit none
real,allocatable,dimension(:),intent(inout)::cpux1,cpux2,cpux3,cpux4,cpux5,cpux6,timex1,timex2,timex3,timex4,timex5,timex6
integer,intent(in)::n
allocate (cpux1(1))
allocate (cpux2(1))
allocate (cpux3(1))
allocate (cpux4(1))
allocate (cpux5(1))
allocate (cpux6(1))
allocate (timex1(1))
allocate (timex2(1))
allocate (timex3(1))
allocate (timex4(1))
allocate (timex5(1))
allocate (timex6(1))
 cpux1(1)=0.0; cpux2(1)=0.0;  cpux3(1)=0.0;  cpux4(1)=0.0;  cpux5(1)=0.0;  cpux6(1)=0.0
  timex1(1)=0.0; timex2(1)=0.0; timex3(1)=0.0;  timex4(1)=0.0;  timex5(1)=0.0;  timex6(1)=0.0
  
end  subroutine timing






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for elements!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine shallocation(ieshape,imaxe)
	!> @brief
!> this subroutine allocates memory for the shapes
	implicit none
	integer,allocatable,dimension(:),intent(inout)::ieshape
	integer,intent(inout)::imaxe
	allocate (ieshape(imaxe))
	allocate (nodes_offset(imaxe),nodes_offset2(imaxe))
	ieshape=0
	nodes_offset=0
	nodes_offset2=0
	
	end subroutine shallocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for elements!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine shdeallocation(ieshape,imaxe)
	!> @brief
!> this subroutine deallocates memory for the shapes
	implicit none
	integer,allocatable,dimension(:),intent(inout)::ieshape
	integer,intent(inout)::imaxe
	deallocate (ieshape,nodes_offset,nodes_offset2)
	end subroutine shdeallocation





!!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for elements!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine elallocation(n,xmpie,xmpielrank,imaxe,ieshape,itestcase,imaxb,xmin,xmax,ymin,ymax,zmin,zmax)
	implicit none
	!> @brief
!> this subroutine allocates memory for the elements
	integer,allocatable,dimension(:),intent(in)::ieshape
	integer,intent(in)::n
	integer,intent(in)::imaxe,itestcase,imaxb
	integer,allocatable,dimension(:),intent(in)::xmpie
	integer,allocatable,dimension(:),intent(in)::xmpielrank
	integer::i,j,k,lm,iex,kmaxe,kk
	real,allocatable,dimension(:),intent(inout)::xmin,xmax,ymin,ymax,zmin,zmax
	allocate(xmin(n:n))
	allocate(ymin(n:n))
	allocate(zmin(n:n))
	allocate(xmax(n:n))
	allocate(ymax(n:n))
	allocate(zmax(n:n))
	kk=0; i=0; j=0; lm=0 
	kmaxe=xmpielrank(n)



  if (dimensiona.eq.3)then

      max_fnodes=4; max_nodes=8




      max_faces=4


      if (num_pyramids.gt.0)then

      max_faces=5;

      end if


      if (num_prisms.gt.0)then

      max_faces=5;

      end if


      if (num_hexas.gt.0)then

      max_faces=6;

      end if


else

max_faces=4; max_fnodes=2; max_nodes=4
end if


!-------------------------
! scalar integer members
!-------------------------
allocate(ielem_ihex(kmaxe) )                ; ielem_ihex = 0
allocate(ielem_ihexgl(kmaxe) )              ; ielem_ihexgl = 0
allocate(ielem_full(kmaxe) )                ; ielem_full = 0
allocate(ielem_interior(kmaxe) )            ; ielem_interior = 0
allocate(ielem_itotalpoints(kmaxe), ielem_troubled(kmaxe) )
allocate(ielem_nofbc(kmaxe) )                ;ielem_nofbc = 0
allocate(ielem_inumneighbours(kmaxe))         ;ielem_inumneighbours=0
allocate(ielem_idegfree(kmaxe))             ;ielem_idegfree=0
allocate(ielem_iorder(kmaxe))               ;ielem_iorder=0
allocate(ielem_mode(kmaxe))                 ;ielem_mode=0
allocate(ielem_indexf(kmaxe))               ;ielem_indexf=0


ielem_itotalpoints = 0                      ; ielem_troubled = 0
allocate(ielem_vdec(kmaxe) )                ; ielem_vdec = 0
allocate(ielem_ggs(kmaxe) )                 ; ielem_ggs = 0
allocate(ielem_ishape(kmaxe) )              ; ielem_ishape = 0
allocate(ielem_ifca(kmaxe) )                ; ielem_ifca = 0
allocate(ielem_admis(kmaxe) )               ; ielem_admis = 0
allocate(ielem_hybrid(kmaxe) )              ; ielem_hybrid = 0
allocate(ielem_nonodes(kmaxe) )             ; ielem_nonodes = 0

if (filtering.eq.1)then
  allocate(ielem_filtered(kmaxe))           ; ielem_filtered = 0
end if
allocate(ielem_reduce(kmaxe) )              ; ielem_reduce = 0


!-------------------------
! scalar real members
!-------------------------
allocate( ielem_totvolume(kmaxe) )          ; ielem_totvolume = 0.0
allocate( ielem_dtl(kmaxe) )                ; ielem_dtl = 0.0
allocate( ielem_minedge(kmaxe) )            ; ielem_minedge = 0.0
allocate( ielem_walldist(kmaxe) )           ; ielem_walldist = 0.0
allocate( ielem_xxc(kmaxe), ielem_yyc(kmaxe), ielem_zzc(kmaxe) )
ielem_xxc = 0.0                             ; ielem_yyc = 0.0                             ; ielem_zzc = 0.0
allocate( ielem_condition(kmaxe) )          ; ielem_condition = 0.0

if (adda.eq.1)then
  allocate(ielem_er(kmaxe), ielem_er2(kmaxe), ielem_er1(kmaxe) )
  ielem_er = 0.0                            ; ielem_er2 = 0.0                              ; ielem_er1 = 0.0
  allocate(ielem_er2dt(kmaxe), ielem_er1dt(kmaxe), ielem_er1er2(kmaxe) )
  ielem_er2dt = 0.0                         ; ielem_er1dt = 0.0                             ; ielem_er1er2 = 0.0
  allocate(ielem_lwcx2(kmaxe), ielem_diss(kmaxe), ielem_erx(kmaxe) )
  ielem_lwcx2 = 0.0                         ; ielem_diss = 0.0                              ; ielem_erx = 0.0
end if

if (mood.gt.0)then
  allocate(ielem_mood(kmaxe), ielem_mood_o(kmaxe) )
  ielem_mood = 0.0                          ; ielem_mood_o = 0.0

end if

allocate(ielem_recalc(kmaxe) )

if (code_profile.eq.30)then
  allocate( ielem_nojecount(max_nodes, kmaxe) )
  allocate(ielem_nodes_neighbours(max_nodes,30,kmaxe))
  ielem_nojecount = 0
  ielem_nodes_neighbours=0
end if

allocate(ielem_linc(kmaxe) )                ; ielem_linc = 0.0
allocate(ielem_wcx(kmaxe) )                 ; ielem_wcx = 0.0

!===========================================================
! fixed-max flattened allocatable members in element_number
!===========================================================

! 1d integer per element (fixed max sizes)

allocate( ielem_nodes(max_nodes, kmaxe) )              ; ielem_nodes = 0
allocate( ielem_types_faces(max_faces, kmaxe) )        ; ielem_types_faces = 0
allocate( ielem_reorient(max_faces, kmaxe) )           ; ielem_reorient = 0

allocate(ielem_indexi(max_faces, kmaxe) )              ; ielem_indexi = 0
allocate(ielem_ibounds(max_faces, kmaxe) )             ; ielem_ibounds = 0
allocate(ielem_ineigh(max_faces, kmaxe) )              ; ielem_ineigh = 0

allocate(ielem_ineighg(max_faces, kmaxe) )             ; ielem_ineighg = 0
allocate(ielem_ineighb(max_faces, kmaxe) )             ; ielem_ineighb = n
allocate(ielem_ineighn(max_faces, kmaxe) )             ; ielem_ineighn = 0


if (tecplot.eq.5)then
allocate(ielem_nodes_v(max_nodes, kmaxe) )             ; ielem_nodes_v = 0
end if

! 2d integer per element
allocate( ielem_nodes_faces(max_faces, max_fnodes, kmaxe) )   ; ielem_nodes_faces = 0
if (tecplot.eq.5)then
allocate( ielem_nodes_faces_v(max_faces, max_fnodes, kmaxe) ) ; ielem_nodes_faces_v = 0
end if

! 1d real per element
allocate( ielem_faceanglex(max_faces, kmaxe) )          ; ielem_faceanglex = 0.0
allocate( ielem_facediss(max_faces, kmaxe) )            ; ielem_facediss = 0.0
allocate( ielem_faceangley(max_faces, kmaxe) )          ; ielem_faceangley = 0.0
allocate( ielem_dih(max_faces,kmaxe) )                 ; ielem_dih = 0.0
allocate( ielem_dih2(max_faces,1:dimensiona, kmaxe) )                 ; ielem_dih2 = 0.0
allocate( ielem_vortex(3, kmaxe) )                      ; ielem_vortex = 0.0
allocate( ielem_surf(max_faces, kmaxe) )                ; ielem_surf = 0.0





	if (tecplot.eq.5)then
	allocate(nodes_offset_local(1:kmaxe));nodes_offset_local=0
	allocate(nodes_offset_local2(1:kmaxe));nodes_offset_local2=0
	end if
! 	allocate(ielem2(n:n,xmpielrank(n)))
	if (itestcase.lt.3)then
		iex=1
	end if
	if (itestcase.ge.3)then
	  if (dimensiona.eq.3)then
		iex=5
	  else
	      iex=4

	  end if
	end if


	end subroutine elallocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for nodes!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine nodeallocation(n,dinode,imaxn,xmpin,xmpinrank,dinoden)
	!> @brief
!> this subroutine allocates memory for the nodes
	implicit none
	integer,intent(in)::n
	integer,allocatable,dimension(:),intent(in)::xmpin
	integer,allocatable,dimension(:),intent(in)::xmpinrank
	integer::i,j,k,lm,iex,kmaxn,kk
	type(anode_number),allocatable,dimension(:,:),intent(inout)::dinode
	type(anode_ne),allocatable,dimension(:),intent(inout)::dinoden
	integer,intent(in)::imaxn
	
!  	allocate (inoden(n:n,imaxn))
!  		  inoden(:,:)%itor=0
	     
	
	end subroutine nodeallocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine nodedeallocation(n,dinode,imaxn,xmpin,xmpinrank,dinoden)
        !> @brief
!> this subroutine deallocates memory from the nodes
	implicit none
	integer,intent(in)::n
	integer,allocatable,dimension(:),intent(in)::xmpin
	integer,allocatable,dimension(:),intent(in)::xmpinrank
	integer::i,j,k,lm,iex,kmaxn,kk
	type(anode_number),allocatable,dimension(:,:),intent(inout)::dinode
	type(anode_ne),allocatable,dimension(:),intent(inout)::dinoden
	integer,intent(in)::imaxn
 	deallocate (dinode)
! 	deallocate (inoden)
	end subroutine nodedeallocation
!---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







subroutine xmpiallocate(xmpie,xmpil,xmpin,xmpielrank,xmpinrank,imaxe,imaxn,nproc)
!> @brief
!> this subroutine allocates memory for the global lists
implicit none
integer,allocatable,dimension(:),intent(inout)::xmpie,xmpil
integer,allocatable,dimension(:),intent(inout)::xmpin
integer,allocatable,dimension(:),intent(inout)::xmpielrank
integer,allocatable,dimension(:),intent(inout)::xmpinrank
integer,intent(inout)::nproc,imaxe,imaxn

allocate (xmpie(imaxe))
allocate (xmpil(imaxe))
allocate (xmpielrank(n:n))

xmpie=0
xmpil=0
xmpielrank=0

end subroutine xmpiallocate


subroutine deallocatempi1(n)
!> @brief
!> this subroutine deallocates memory for boundary exchange
implicit none
integer,intent(in)::n

deallocate (diexchanges1,diexchanger1)

end subroutine deallocatempi1

subroutine deallocatempi2(n)
!> @brief
!> this subroutine deallocates memory for the stencil exchange
implicit none
integer,intent(in)::n

deallocate (direcexr1,direcexs1)

end subroutine deallocatempi2


subroutine local_reconallocation3(n)
  !> @brief
  !> this subroutine allocates memory for reconstruction prestoring
  implicit none
  integer, intent(in) :: n
  integer :: i, j, k, m, ikg, itrue, kmaxe, idum
  integer :: inum_points, itarget, int_wall, idx, int_local
  integer :: imax, inum, ideg, imax2, inum2, ideg2, m2
  real    :: perc, perde, perd, per1, per2, per3, per4, per5, per0
  real    :: pef0, pef1, pef2, pef3, pef4, pef5, perv, per01, pef01

  kmaxe = xmpielrank(n)

  call mpi_barrier(mpi_comm_world, ierror)

  allocate(rec_local(1:kmaxe));                            rec_local = 0
  allocate(rec_invccjac(1:dimensiona,1:dimensiona,1:kmaxe)); rec_invccjac = zero
  allocate(rec_vext_ref(1:dimensiona,1:kmaxe));            rec_vext_ref = zero
  allocate(rec_wall(1:kmaxe));                             rec_wall = 0
  allocate(rec_volume(1,1,1:kmaxe))                 ;       rec_volume=0

  if (ees == 5) then
    allocate(rec_ihexg(1,numneighbours,1:kmaxe));          rec_ihexg = 0
    allocate(rec_ihexl(1,numneighbours,1:kmaxe));          rec_ihexl = 0
    if (iperiodicity == 1) then
      allocate(rec_periodicflag(1,numneighbours,1:kmaxe)); rec_periodicflag = 0
    end if
    allocate(rec_ihexgc(2:typesten,numneighbours2,1:kmaxe)); rec_ihexgc = 0
    allocate(rec_ihexlc(2:typesten,numneighbours2,1:kmaxe)); rec_ihexlc = 0
  else
    allocate(rec_ihexg(1:typesten,numneighbours,1:kmaxe));  rec_ihexg = 0
    allocate(rec_ihexl(1:typesten,numneighbours,1:kmaxe));  rec_ihexl = 0
    if (iperiodicity == 1) then
      allocate(rec_periodicflag(1:typesten,numneighbours,1:kmaxe)); rec_periodicflag = 0
    end if
  end if

  if (ees == 5) then
    allocate(rec_invmat_stencilt(idegfree,numneighbours-1,1,1:kmaxe));             rec_invmat_stencilt = zero
    allocate(rec_invmat_stenciltc(idegfree2,numneighbours2-1,2:typesten,1:kmaxe)); rec_invmat_stenciltc = zero
  else
    allocate(rec_invmat_stencilt(idegfree,numneighbours-1,1:typesten,1:kmaxe));    rec_invmat_stencilt = zero
  end if

  int_wall = 0
  do i = 1, kmaxe
    idum = 0
    if (ielem_interior(i) == 1) then
      do j = 1, ielem_ifca(i)
        if (ielem_ibounds(j,i) > 0) then
          if (ibound_icode(ielem_ibounds(j,i)) == 4) then
            idum = 1
          end if
        end if
      end do
    end if
    int_wall = int_wall + idum
  end do

  int_wall = 0
  do i = 1, kmaxe
    idum = 0
    if (ielem_interior(i) == 1) then
      do j = 1, ielem_ifca(i)
        if (ielem_ibounds(j,i) > 0) then
          if (ibound_icode(ielem_ibounds(j,i)) == 4) then
            idum = 1
          end if
        end if
      end do
    end if
    int_wall = int_wall + idum

    if (idum.eq.1)then

    rec_wall(i) = int_wall

    else

    rec_wall(i) =0
    end if




  end do

  if (dimensiona == 3) then
    idx = max_faces
  else
    idx = max_faces
  end if

  allocate(rec_grads(nof_variables-1+turbulenceequations+passivescalar,1:dimensiona,1:kmaxe)); rec_grads = zero
  allocate(rec_uleft(1:nof_variables,idx,1:numberofpoints2,1:kmaxe));           rec_uleft = zero

  if (dg.eq.1)then
  allocate(rec_uleft_dg(1:nof_variables,idx,1:numberofpoints2,1:kmaxe)); rec_uleft_dg=zero
  end if


   if ((turbulenceequations > 0) .or. (passivescalar > 0)) then
  allocate(rec_uleftturb(1:turbulenceequations+passivescalar,idx,1:numberofpoints2,1:kmaxe)); rec_uleftturb = zero
  end if

  if (ees == 5) then
    allocate(rec_gradients(1,1:idegfree,1:nof_variables,1:kmaxe));    rec_gradients = zero
    allocate(rec_gradientsc(2:typesten,1:idegfree2,1:nof_variables,1:kmaxe)); rec_gradientsc = zero
    if ((turbulenceequations > 0) .or. (passivescalar > 0)) then
      allocate(rec_gradients2(1,1:idegfree,1:turbulenceequations+passivescalar,1:kmaxe));      rec_gradients2 = zero
      allocate(rec_gradientsc2(2:typesten,1:idegfree2,1:turbulenceequations+passivescalar,1:kmaxe)); rec_gradientsc2 = zero
    end if
  else
    allocate(rec_gradients(1:typesten,1:idegfree,1:nof_variables,1:kmaxe)); rec_gradients = zero
    if ((turbulenceequations > 0) .or. (passivescalar > 0)) then
      allocate(rec_gradients2(1:typesten,1:idegfree,1:turbulenceequations+passivescalar,1:kmaxe)); rec_gradients2 = zero
    end if
  end if

  if (itestcase >= 4) then

    allocate(rec_gradf(1:nof_variables-1,1:idegfree,1:kmaxe)); rec_gradf = zero
    allocate(rec_uleftv(1:dimensiona,1:nof_variables-1,idx,1:numberofpoints2,1:kmaxe)); rec_uleftv = zero

    if (averaging == 1) then
      allocate(rec_gradsav(1:nof_variables-1,1:dimensiona,1:kmaxe)); rec_gradsav = zero
    end if

    if ((turbulence == 1) .or. (passivescalar > 0)) then
      allocate(rec_uleftturbv(1:dimensiona,1:turbulenceequations+passivescalar,idx,1:numberofpoints2,1:kmaxe)); rec_uleftturbv = zero
      allocate(rec_gradientsturb(1,1:idegfree,1:turbulenceequations+passivescalar,1:kmaxe));                   rec_gradientsturb = zero
    end if

    if (greengo == 0) then

      if (int_wall > 0) then
        if (ees == 5) then
          allocate(rec_stencils(1,numneighbours-1,1:idegfree,1:int_wall));             rec_stencils = 0
!           allocate(rec_stencilsc(2:typesten,numneighbours2-1,1:idegfree2,1:int_wall)); rec_stencilsc = 0
          allocate(rec_weightl(1,numneighbours-1,1:int_wall));                         rec_weightl = zero
        else
          allocate(rec_stencils(1:typesten,numneighbours-1,1:idegfree,1:int_wall));    rec_stencils = 0
          allocate(rec_weightl(1,numneighbours-1,1:int_wall));                         rec_weightl = zero
        end if
        allocate(rec_volume_w(1,numneighbours,1:int_wall));                 rec_volume_w = zero
        allocate(rec_k0(1:int_wall));                                      rec_k0 = 0
        allocate(rec_g0(1:int_wall));                                      rec_g0 = 0
        allocate(rec_velinvlsqmat(idegfree-1,idegfree-1,1:int_wall));      rec_velinvlsqmat = zero
        allocate(rec_wallcoeff(idegfree,1:int_wall));                      rec_wallcoeff = zero
        allocate(rec_vellsq(numneighbours-1,idegfree-1,1:int_wall));       rec_vellsq = zero
        allocate(rec_tempsqmat(idegfree-1,idegfree-1,1:int_wall));         rec_tempsqmat = zero
        allocate(rec_wallcoefg(idegfree,1:int_wall));                      rec_wallcoefg = zero
        allocate(rec_tempsq(numneighbours-1,idegfree-1,1:int_wall));       rec_tempsq = zero
      end if

    end if
  end if

  int_local = 0
  do i = 1, kmaxe
    itrue = 0
    do j = 1, typesten
      ikg = 0
      if ((ees /= 5) .or. (j == 1)) then
        itarget = numneighbours
      else
        itarget = numneighbours2
      end if
      do k = 1, itarget
        if (ilocalstencil(n,i,j,k) > 0) then
          ikg = ikg + 1
          if (xmpie(ilocalstencil(n,i,j,k)) /= n) then
            itrue = 1
          end if
        end if
      end do
    end do
    if (itrue .ne. 0) then
      int_local = int_local + 1
      rec_local(i) = int_local
    end if
  end do

  if (int_local > 0) then
    if (ees == 5) then
      allocate(rec_ihexb(1,1:numneighbours,1:int_local));                rec_ihexb = -100
      allocate(rec_ihexn(1,1:numneighbours,1:int_local));                rec_ihexn = -100
      allocate(rec_ihexbc(2:typesten,numneighbours2,1:int_local));       rec_ihexbc = -100
      allocate(rec_ihexnc(2:typesten,numneighbours2,1:int_local));       rec_ihexnc = -100
    else
      allocate(rec_ihexb(1:typesten,1:numneighbours,1:int_local));       rec_ihexb = -100
      allocate(rec_ihexn(1:typesten,1:numneighbours,1:int_local));       rec_ihexn = -100
    end if
  end if

  if (iweno == 1) then
    allocate(rec_indicator(1:idegfree,1:idegfree,1:kmaxe));              rec_indicator = zero
    if (ees == 5) then
      allocate(rec_indicatorc(1:idegfree2,1:idegfree2,1:kmaxe));         rec_indicatorc = zero
    end if
  end if

  if (dg == 1) then
    allocate(dg2fv(1:idegfree,1:nof_variables,1:kmaxe));             dg2fv = zero
  end if












end subroutine local_reconallocation3








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dealcordinates1(n,diexcordr,diexcords)
!> @brief
!> this subroutine deallocates memory for exchange of info between processes
implicit none
type(aexchange_cord),allocatable,dimension(:),intent(inout)::diexcordr
type(aexchange_cord),allocatable,dimension(:),intent(inout)::diexcords
integer,intent(in)::n
deallocate(diexcordr)


end subroutine dealcordinates1


subroutine dealcordinates2
implicit none
!> @brief
!> this subroutine deallocates memory for exchange of info between processes
 if (allocated(diexcords))deallocate(diexcords)

end subroutine dealcordinates2





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine allocate_basis_function(n,xmpielrank,idegfree)
!> @brief
!> this subroutine allocates memory for basis function integrals
implicit none
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer,intent(in)::idegfree,n
integer::kmaxe,i
kmaxe=xmpielrank(n)

allocate(integ_basis_value(1:idegfree,1:kmaxe));integ_basis_value=zero
if (ees.eq.5)then
allocate(integ_basis_valuec(1:idegfree2,1:kmaxe));integ_basis_valuec=zero
end if

if (dg.eq.1)then
allocate(integ_basis_dg_value(1:idegfree,1:kmaxe));integ_basis_dg_value=zero
end if


end subroutine allocate_basis_function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





   subroutine localsdeallocation(n,xmpielrank,ilocalstencil,ilocalstencilper,typesten,numneighbours)
   !> @brief
!> this subroutine allocates memory for stencils
	implicit none
	integer,allocatable,dimension(:,:,:,:),intent(inout)::ilocalstencil,ilocalstencilper
	integer,intent(in)::numneighbours
	integer,intent(in)::n
	integer,allocatable,dimension(:),intent(in)::xmpielrank
	integer,intent(in)::typesten
	
	deallocate (ilocalstencil)
	deallocate (ilocalstencilper)
	end subroutine localsdeallocation
	
	
subroutine u_c_allocation(n, xmpielrank, itestcase)
  implicit none
  integer, intent(in) :: n, itestcase
  integer, allocatable,dimension(:),intent(in) :: xmpielrank

  integer :: kmaxe, istage
  integer :: max_ideg

  kmaxe = xmpielrank(n)


  call mpi_barrier(mpi_comm_world, ierror)

  ! -----------------------------
  ! decide number of rk stages
  ! -----------------------------
  select case (rungekutta)
  case (1)
    istage = 1
  case (2)
    istage = 2
  case (3)
    if (averaging == 1) then
      istage = 5
    else
      istage = 3
      if (mood == 1) istage = 4
    end if
  case (4)
    if (averaging == 1) then
      istage = 7
    else
      istage = 6
    end if
  case (5)
    istage = 2
  case (10)
    istage = 1
  case (11,12)
    if (averaging == 1) then
      istage = 5
    else
      istage = 3
    end if
  case default
    istage = 1
  end select

  ! -----------------------------
  ! base solution arrays (flat)
  ! u_c_val(j,k,i) -> u_c_val(j,k,i)
  ! -----------------------------
  if (.not. allocated(u_c_val)) then
    allocate(u_c_val(istage, nof_variables, kmaxe)); u_c_val = zero
  end if

  if ((turbulence > 0) .or. (passivescalar > 0)) then
    if (.not. allocated(u_ct_val)) then
      allocate(u_ct_val(istage, turbulenceequations + passivescalar, kmaxe)); u_ct_val = zero
    end if
  end if

  if (itestcase <= 4) then
    if (.not. allocated(u_e_val)) then
      allocate(u_e_val(1, nof_variables, kmaxe)); u_e_val = zero
    end if
  end if

  if (averaging == 1) then
    if (.not. allocated(u_c_rms)) then
      allocate(u_c_rms(nof_variables, kmaxe)); u_c_rms = zero
    end if
  end if

  ! -----------------------------
  ! dg arrays (flat)
  ! u_c_valdg(j,k,l,i) -> u_c_valdg(j,k,l,i)
  !
  ! important: ideg varies per element; allocate with max degree
  ! and use slices 1:(ielem_idegfree(i)+1) at use sites.
  ! -----------------------------
  if (dg == 1) then
    max_ideg = idegfree

    if (.not. allocated(u_c_valdg)) then
      allocate(u_c_valdg(istage,nof_variables, max_ideg+1, kmaxe)); u_c_valdg = zero
    end if

    if (.not. allocated(m_1_val)) then
      ! if your m_1 is per-element matrix in modal/nodal space, this is the usual shape
      allocate(m_1_val(max_ideg+1, max_ideg+1, kmaxe)); m_1_val = zero
    end if

    if (itestcase == 4) then
      if (.not. allocated(u_c_br2_aux_var)) then
        ! your declared br2 aux var is 4d: (j,k,l,i) in the mapping you used earlier
        allocate(u_c_br2_aux_var(max_ideg+1, max_ideg+1, nof_variables, kmaxe)); u_c_br2_aux_var = zero
      end if
    end if
  end if

  ! -----------------------------
  ! filtering arrays (flat)
  ! u_cs(i)%val(...) -> u_cs_val(...)
  ! u_cw(i)%val(...) -> u_cw_val(...)
  ! and dg variants: u_cs_valdg / u_cw_valdg
  ! -----------------------------
  if (filtering == 1) then
    if (.not. allocated(u_cs_val)) then
      allocate(u_cs_val(istage, nof_variables, kmaxe)); u_cs_val = zero
    end if
    if (.not. allocated(u_cw_val)) then
      allocate(u_cw_val(istage, nof_variables, kmaxe)); u_cw_val = zero
    end if

    if (dg == 1) then
      if (.not. allocated(u_cs_valdg)) then
        allocate(u_cs_valdg(1, nof_variables, max_ideg+1, kmaxe)); u_cs_valdg = zero
      end if
      if (.not. allocated(u_cw_valdg)) then
        allocate(u_cw_valdg(1, nof_variables, max_ideg+1, kmaxe)); u_cw_valdg = zero
      end if
    end if
  end if

  ! -----------------------------
  ! set recalc flags
  ! ielem_recalc(i) already flat
  ! -----------------------------
  if (mood == 1) then
    ielem_recalc(1:kmaxe) = 0
  else
    ielem_recalc(1:kmaxe) = 1
  end if

end subroutine u_c_allocation



subroutine omp_map_first(n)
  implicit none
  integer,intent(in)::n


  
#ifdef gpu
  !$omp target enter data map(alloc: n, aa_1, adda, adda_1, adda_1_s, adda_2, adda_2_s, adda_alpha_1, adda_alpha_2)
  !$omp target update to(n, aa_1, adda, adda_1, adda_1_s, adda_2, adda_2_s, adda_alpha_1, adda_alpha_2)
  !$omp target enter data map(alloc: adda_type, allnodesgloball, allres, allresdt, alls, alpha, alpha_0, alpha_inf1)
  !$omp target update to(adda_type, allnodesgloball, allres, allresdt, alls, alpha, alpha_0, alpha_inf1)
  !$omp target enter data map(alloc: alpha_inf2, alpha_star0, alpha_starinf, angle_per, aoa, average_restart, averaging)
  !$omp target update to(alpha_inf2, alpha_star0, alpha_starinf, angle_per, aoa, average_restart, averaging)
  !$omp target enter data map(alloc: beta, beta_i1, beta_i2, beta_starinf, beta_t, betaas, binio, bleed, bleed_number)
  !$omp target update to(beta, beta_i1, beta_i2, beta_starinf, beta_t, betaas, binio, bleed, bleed_number)
  !$omp target enter data map(alloc: bleed_type, boundtype, br2_damping, br2_yn, bubble_centre, bubble_radius, c_des_sa)
  !$omp target update to(bleed_type, boundtype, br2_damping, br2_yn, bubble_centre, bubble_radius, c_des_sa)
  !$omp target enter data map(alloc: c_des_sst, c_mu_inlet, c_sas, c_smg, cascade, catalytic_wall, cavitation, cb1, cb2)
  !$omp target update to(c_des_sst, c_mu_inlet, c_sas, c_smg, cascade, catalytic_wall, cavitation, cb1, cb2)
  !$omp target enter data map(alloc: cfl, cflmax, cflramp, cfw, charlength, chunk_n, code_profile, ct1, ct2, ct3, ct4, cv1)
  !$omp target update to(cfl, cflmax, cflramp, cfw, charlength, chunk_n, code_profile, ct1, ct2, ct3, ct4, cv1)
  !$omp target enter data map(alloc: cw1, cw2, cw3, d_corr, datatypeint, datatypex, datatypexx, datatypey, datatypeyy)
  !$omp target update to(cw1, cw2, cw3, d_corr, datatypeint, datatypex, datatypexx, datatypey, datatypeyy)
  !$omp target enter data map(alloc: datatypez, des_model, dg, dimensiona, dims, dt, ees, ek_time, emetis, eta2_sas)
  !$omp target update to(datatypez, des_model, dg, dimensiona, dims, dt, ees, ek_time, emetis, eta2_sas)
  !$omp target enter data map(alloc: every_time, extended_bounds, extf, fastest, fastest_q,ind1,origin)
  !$omp target update to(every_time, extended_bounds, extf, fastest, fastest_q,ind1,origin)

  !$omp target enter data map(alloc: fastmovie, fil_alpha, fil_nc, fil_s, filter_type, filtering, firstorder, firstrese,totk,totens,totensx,kill_nan)
  !$omp target update to(fastmovie, fil_alpha, fil_nc, fil_s, filter_type, filtering, firstorder, firstrese,totk,totens,totensx,kill_nan)
  !$omp target enter data map(alloc: firstresk, firstresomega, firstrespass, firstresr, firstrest, firstresu, firstresv)
  !$omp target update to(firstresk, firstresomega, firstrespass, firstresr, firstrest, firstresu, firstresv)
  !$omp target enter data map(alloc: firstresw, forcex, forcey, forcez, gamma, governingequations, greengo, gridar1)
  !$omp target update to(firstresw, forcex, forcey, forcez, gamma, governingequations, greengo, gridar1)
  !$omp target enter data map(alloc: gridar2, guassianquadra, hybridist, i_turb_inlet, iadapt, ibcode, iboundary, ibside)
  !$omp target update to(gridar2, guassianquadra, hybridist, i_turb_inlet, iadapt, ibcode, iboundary, ibside)
  !$omp target enter data map(alloc: icarlos1, icarlos2, icompact, icong, iconimp, iconsgvq, iconsr, icoupleturb, idegfree)
  !$omp target update to(icarlos1, icarlos2, icompact, icong, iconimp, iconsgvq, iconsr, icoupleturb, idegfree)
  !$omp target enter data map(alloc: idegfree2, idegfree3, ievery, ievery2, ieveryav, iforce, igianagraps, igqrules, ihax1)
  !$omp target update to(idegfree2, idegfree3, ievery, ievery2, ieveryav, iforce, igianagraps, igqrules, ihax1)
  !$omp target enter data map(alloc: ihybrid, iloop, iloopx, ilx, imaxb, imaxdegfree, imaxdegfree2, imaxe, imaxn, in)
  !$omp target update to(ihybrid, iloop, iloopx, ilx, imaxb, imaxdegfree, imaxdegfree2, imaxe, imaxn, in)
  !$omp target enter data map(alloc: indicator_type, init_mu_ratio, initcond, initialres, inum2, inwhichel, iorder)
  !$omp target update to(indicator_type, init_mu_ratio, initcond, initialres, inum2, inwhichel, iorder)
  !$omp target enter data map(alloc: iorder2, ioverst, ioverto, iperiodicity, ires_turb, ires_unsteady, iriemann, irs)
  !$omp target update to(iorder2, ioverst, ioverto, iperiodicity, ires_turb, ires_unsteady, iriemann, irs)
  !$omp target enter data map(alloc: ischeme, iscoun, iselem, ispal, isplit, issf)
  !$omp target update to(ischeme, iscoun, iselem, ispal, isplit, issf)

  !$omp target enter data map(alloc: istn, it, itestcase, itold, itotalb, itt, ivortex, iweightlsqr, iweno, iwmaxe, jk)
  !$omp target update to(istn, it, itestcase, itold, itotalb, itt, ivortex, iweightlsqr, iweno, iwmaxe, jk)
  !$omp target enter data map(alloc: jtotal, jtotal1, jtotal2, jtotal3, jump_cond1, jump_cond2, jump_cond3, kappa,a405,nof_perturbations405)
  !$omp target update to(jtotal, jtotal1, jtotal2, jtotal3, jump_cond1, jump_cond2, jump_cond3, kappa,a405,nof_perturbations405)
  !$omp target enter data map(alloc: kappa_sst, kdum1, kdum2, kdum3, kill, kinit_srf, kloopx, kmaxn, l0norm, l1norm)
  !$omp target update to(kappa_sst, kdum1, kdum2, kdum3, kill, kinit_srf, kloopx, kmaxn, l0norm, l1norm)
  !$omp target enter data map(alloc: l2norm, l_turb_inlet, lam, lamps, lamx, lamy, lamz, limiter, lmach, lmach_style)
  !$omp target update to(l2norm, l_turb_inlet, lam, lamps, lamx, lamy, lamz, limiter, lmach, lmach_style)
  !$omp target enter data map(alloc: lowmem, lowmemory, lwci1, m_t0, modeio, momentx, momenty, momentz, mood, mood_mode)
  !$omp target update to(lowmem, lowmemory, lwci1, m_t0, modeio, momentx, momenty, momentz, mood, mood_mode)
  !$omp target enter data map(alloc: mood_var1, mood_var2, mood_var3, mood_var4, movement, mp_modelc, mrf, multispecies)
  !$omp target update to(mood_var1, mood_var2, mood_var3, mood_var4, movement, mp_modelc, mrf, multispecies)
  !$omp target enter data map(alloc: n_boundaries, nderivative, nodes_i, nodes_part, nof_bounded, nof_bubbles)
  !$omp target update to(n_boundaries, nderivative, nodes_i, nodes_part, nof_bounded, nof_bubbles)
  !$omp target enter data map(alloc: nof_interior, nof_species, nof_variables, nprobes, nproc, nrotors, ntmax, num_dg_dofs)
  !$omp target update to(nof_interior, nof_species, nof_variables, nprobes, nproc, nrotors, ntmax, num_dg_dofs)
  !$omp target enter data map(alloc: num_dg_reconstruct_dofs, numberofpoints, numberofpoints2, numneighbours)
  !$omp target update to(num_dg_reconstruct_dofs, numberofpoints, numberofpoints2, numneighbours)
  !$omp target enter data map(alloc: numneighbours2, oo2, out_time, output_freq, outsurf, part1_end, part2_end)
  !$omp target update to(numneighbours2, oo2, out_time, output_freq, outsurf, part1_end, part2_end)

  !$omp target enter data map(alloc: part3_end, part4_end, part5_end, passivescalar, per_rot, pi, poly, pr_t1, pr_t2)
  !$omp target update to(part3_end, part4_end, part5_end, passivescalar, per_rot, pi, poly, pr_t1, pr_t2)
  !$omp target enter data map(alloc: pr_t3, pr_t4, pr_t5, pr_t6, pr_t7, pr_t8, prace_t1, prace_t2, prace_t3, prace_t4)
  !$omp target update to(pr_t3, pr_t4, pr_t5, pr_t6, pr_t7, pr_t8, prace_t1, prace_t2, prace_t3, prace_t4)
  !$omp target enter data map(alloc: prace_t5, prace_t6, prace_t7, prace_t8, prace_t9, prace_tx1, prace_tx2, prace_tx3)
  !$omp target update to(prace_t5, prace_t6, prace_t7, prace_t8, prace_t9, prace_tx1, prace_tx2, prace_tx3)
  !$omp target enter data map(alloc: ccfl,prandtl, pres, press_outlet, prev_turbmodel, prevres, prtu, qp_hexa, qp_line)
  !$omp target update to(ccfl,prandtl, pres, press_outlet, prev_turbmodel, prevres, prtu, qp_hexa, qp_line)
  !$omp target enter data map(alloc: qp_line_n, qp_prism, qp_pyra, qp_quad, qp_quad_n, qp_tetra, qp_triangle)
  !$omp target update to(qp_line_n, qp_prism, qp_pyra, qp_quad, qp_quad_n, qp_tetra, qp_triangle)
  !$omp target enter data map(alloc: qp_triangle_n, qrde, qsas_model, r_beta, r_gas, r_k_sst, r_om_sst, realgas)
  !$omp target update to(qp_triangle_n, qrde, qsas_model, r_beta, r_gas, r_k_sst, r_om_sst, realgas)
  !$omp target enter data map(alloc: reduce_comp, relax, required, res_time, rescounter, rescountert, residualfreq)
  !$omp target update to(reduce_comp, relax, required, res_time, rescounter, rescountert, residualfreq)
  !$omp target enter data map(alloc: reslimit, resmax, resmaxt, restart, reynolds, rframe, rg_kf_type, rg_nof_reactions)
  !$omp target update to(reslimit, resmax, resmaxt, restart, reynolds, rframe, rg_kf_type, rg_nof_reactions)
  !$omp target enter data map(alloc: rg_nof_tv_coef, rg_relax, rg_t_inf, rg_t_ref, rg_t_wall_init, rg_ttr, rg_tve, rhc1)
  !$omp target update to(rg_nof_tv_coef, rg_relax, rg_t_inf, rg_t_ref, rg_t_wall_init, rg_ttr, rg_tve, rhc1)
  !$omp target enter data map(alloc: rhc2, rhc3, rhc4, rot_corr, rres, rungekutta, scaler, schmidt_lam)
  !$omp target update to(rhc2, rhc3, rhc4, rot_corr, rres, rungekutta, scaler, schmidt_lam)

  !$omp target enter data map(alloc: schmidt_turb, sigma, sigma_k1, sigma_k2, sigma_om1, sigma_om2, sigma_phi)
  !$omp target update to(schmidt_turb, sigma, sigma_k1, sigma_k2, sigma_om1, sigma_om2, sigma_phi)
  !$omp target enter data map(alloc: source_active, spatialorder, spatiladiscret, spkin, spos, srf_origin, srf_velocity)
  !$omp target update to(source_active, spatialorder, spatiladiscret, spkin, spos, srf_origin, srf_velocity)
  !$omp target enter data map(alloc: srfg, st_n_cpu, st_n_threads, statfile, statistics, stencil_io, stennorm, subdiv)
  !$omp target update to(srfg, st_n_cpu, st_n_threads, statfile, statistics, stencil_io, stennorm, subdiv)
  !$omp target enter data map(alloc: surfshear, suther, swirl, t, taylor, taylor_ens, taylor_ensx, tecplot, temp_model)
  !$omp target update to(surfshear, suther, swirl, t, taylor, taylor_ens, taylor_ensx, tecplot, temp_model)
  !$omp target enter data map(alloc: temporder, thermal, thread_n, timestep, tol_per, tolbig, tolsmall, totalvolume, totiw)
  !$omp target update to(temporder, thermal, thread_n, timestep, tol_per, tolbig, tolsmall, totalvolume, totiw)
  !$omp target enter data map(alloc: totwalls, totwallsc, turbinit, turbulence, turbulenceequations, turbulencemodel)
  !$omp target update to(totwalls, totwallsc, turbinit, turbulence, turbulenceequations, turbulencemodel)
  !$omp target enter data map(alloc: twall, typ_countn, typ_countn_global, typ_countn_global_w, typ_countn_w, typesten)
  !$omp target update to(twall, typ_countn, typ_countn_global, typ_countn_global_w, typ_countn_w, typesten)
  !$omp target enter data map(alloc: tz1, ufreestream, unwou, upperlimit, upturblimit, uvel, v_ref, variable_names)
  !$omp target update to(tz1, ufreestream, unwou, upperlimit, upturblimit, uvel, v_ref, variable_names)
  !$omp target enter data map(alloc: variable_names_av, variable_names_av_w, variable_names_w, vectorx, vectory, vectorz)
  !$omp target update to(variable_names_av, variable_names_av_w, variable_names_w, vectorx, vectory, vectorz)
  !$omp target enter data map(alloc: visc, viscous_s, voll, vorder, vort_model, vvel, wall_temp)
  !$omp target update to(visc, viscous_s, voll, vorder, vort_model, vvel, wall_temp)

  !$omp target enter data map(alloc: wallc, wdatatypeint, wdatatypex, wdatatypexx, wdatatypey, wdatatypeyy, wdatatypez)
  !$omp target update to(wallc, wdatatypeint, wdatatypex, wdatatypexx, wdatatypey, wdatatypeyy, wdatatypez)
  !$omp target enter data map(alloc: weight_lsqr, wenocentralweight, wenocnschar, wenoz, wenwrt, wkdum1, wkdum2, wkdum3)
  !$omp target update to(weight_lsqr, wenocentralweight, wenocnschar, wenoz, wenwrt, wkdum1, wkdum2, wkdum3)
  !$omp target enter data map(alloc: wnodes_part, wpart1_end, wpart2_end, wpart3_end, wpart4_end, wpart5_end)
  !$omp target update to(wnodes_part, wpart1_end, wpart2_end, wpart3_end, wpart4_end, wpart5_end)
  !$omp target enter data map(alloc: write_variables, write_variables_av, write_variables_av_w, write_variables_w, wvel)
  !$omp target update to(write_variables, write_variables_av, write_variables_av_w, write_variables_w, wvel)
  !$omp target enter data map(alloc: xper, yper, zero, zero_turb_init, zeta_star, zper)
  !$omp target update to(xper, yper, zero, zero_turb_init, zeta_star, zper)
  !$omp target enter data map(alloc: indicator_par1, indicator_par2, indicator_par3)
  !$omp target update to(indicator_par1, indicator_par2, indicator_par3)


  if (allocated(jtot)) then
    !$omp target enter data map(alloc: jtot)
    !$omp target update to(jtot)
  end if

  if (allocated(adda_filter_strong)) then
    !$omp target enter data map(alloc: adda_filter_strong)
    !$omp target update to(adda_filter_strong)
  end if
  if (allocated(adda_filter_weak)) then
    !$omp target enter data map(alloc: adda_filter_weak)
    !$omp target update to(adda_filter_weak)
  end if
  if (allocated(bleed_end)) then
    !$omp target enter data map(alloc: bleed_end)
    !$omp target update to(bleed_end)
  end if
  if (allocated(bleed_plenum)) then
    !$omp target enter data map(alloc: bleed_plenum)
    !$omp target update to(bleed_plenum)
  end if
  if (allocated(bleed_porosity)) then
    !$omp target enter data map(alloc: bleed_porosity)
    !$omp target update to(bleed_porosity)
  end if
  if (allocated(bleed_start)) then
    !$omp target enter data map(alloc: bleed_start)
    !$omp target update to(bleed_start)
  end if
  if (allocated(bound_len)) then
    !$omp target enter data map(alloc: bound_len)
    !$omp target update to(bound_len)
  end if
  if (allocated(bound_offset)) then
    !$omp target enter data map(alloc: bound_offset)
    !$omp target update to(bound_offset)
  end if
  if (allocated(boundhir_dg)) then
    !$omp target enter data map(alloc: boundhir_dg)
    !$omp target update to(boundhir_dg)
  end if
  if (allocated(catalytic_con)) then
    !$omp target enter data map(alloc: catalytic_con)
    !$omp target update to(catalytic_con)
  end if
  if (allocated(el_bnd)) then
    !$omp target enter data map(alloc: el_bnd)
    !$omp target update to(el_bnd)
  end if
  if (allocated(el_int)) then
    !$omp target enter data map(alloc: el_int)
    !$omp target update to(el_int)
  end if
  if (allocated(gamma_in)) then
    !$omp target enter data map(alloc: gamma_in)
    !$omp target update to(gamma_in)
  end if
  if (allocated(halo_offset)) then
    !$omp target enter data map(alloc: halo_offset)
    !$omp target update to(halo_offset)
  end if
  if (allocated(ibound_cpun)) then
    !$omp target enter data map(alloc: ibound_cpun)
    !$omp target update to(ibound_cpun)
  end if
  if (allocated(ibound_face)) then
    !$omp target enter data map(alloc: ibound_face)
    !$omp target update to(ibound_face)
  end if
  if (allocated(ibound_ibid)) then
    !$omp target enter data map(alloc: ibound_ibid)
    !$omp target update to(ibound_ibid)
  end if
  if (allocated(ibound_ibl)) then
    !$omp target enter data map(alloc: ibound_ibl)
    !$omp target update to(ibound_ibl)
  end if
  if (allocated(ibound_icode)) then
    !$omp target enter data map(alloc: ibound_icode)
    !$omp target update to(ibound_icode)
  end if
  if (allocated(ibound_inum)) then
    !$omp target enter data map(alloc: ibound_inum)
    !$omp target update to(ibound_inum)
  end if
  if (allocated(ibound_ishape)) then
    !$omp target enter data map(alloc: ibound_ishape)
    !$omp target update to(ibound_ishape)
  end if
  if (allocated(ibound_localn)) then
    !$omp target enter data map(alloc: ibound_localn)
    !$omp target update to(ibound_localn)
  end if
  if (allocated(ibound_nibl)) then
    !$omp target enter data map(alloc: ibound_nibl)
    !$omp target update to(ibound_nibl)
  end if
  if (allocated(ibound_nlocal)) then
    !$omp target enter data map(alloc: ibound_nlocal)
    !$omp target update to(ibound_nlocal)
  end if
  if (allocated(ibound_t)) then
    !$omp target enter data map(alloc: ibound_t)
    !$omp target update to(ibound_t)
  end if
  if (allocated(ibound_t2)) then
    !$omp target enter data map(alloc: ibound_t2)
    !$omp target update to(ibound_t2)
  end if
  if (allocated(ibound_which)) then
    !$omp target enter data map(alloc: ibound_which)
    !$omp target update to(ibound_which)
  end if
  if (allocated(ielem_admis)) then
    !$omp target enter data map(alloc: ielem_admis)
    !$omp target update to(ielem_admis)
  end if
  if (allocated(ielem_avars)) then
    !$omp target enter data map(alloc: ielem_avars)
    !$omp target update to(ielem_avars)
  end if
  if (allocated(ielem_bleedn)) then
    !$omp target enter data map(alloc: ielem_bleedn)
    !$omp target update to(ielem_bleedn)
  end if
  if (allocated(ielem_condition)) then
    !$omp target enter data map(alloc: ielem_condition)
    !$omp target update to(ielem_condition)
  end if
  if (allocated(ielem_condx)) then
    !$omp target enter data map(alloc: ielem_condx)
    !$omp target update to(ielem_condx)
  end if
  if (allocated(ielem_dih)) then
    !$omp target enter data map(alloc: ielem_dih)
    !$omp target update to(ielem_dih)
  end if
  if (allocated(ielem_dih2)) then
    !$omp target enter data map(alloc: ielem_dih2)
    !$omp target update to(ielem_dih2)
  end if
  if (allocated(ielem_diss)) then
    !$omp target enter data map(alloc: ielem_diss)
    !$omp target update to(ielem_diss)
  end if
  if (allocated(ielem_dtl)) then
    !$omp target enter data map(alloc: ielem_dtl)
    !$omp target update to(ielem_dtl)
  end if
  if (allocated(ielem_er)) then
    !$omp target enter data map(alloc: ielem_er)
    !$omp target update to(ielem_er)
  end if
  if (allocated(ielem_er1)) then
    !$omp target enter data map(alloc: ielem_er1)
    !$omp target update to(ielem_er1)
  end if
  if (allocated(ielem_er1dt)) then
    !$omp target enter data map(alloc: ielem_er1dt)
    !$omp target update to(ielem_er1dt)
  end if
  if (allocated(ielem_er1er2)) then
    !$omp target enter data map(alloc: ielem_er1er2)
    !$omp target update to(ielem_er1er2)
  end if
  if (allocated(ielem_er2)) then
    !$omp target enter data map(alloc: ielem_er2)
    !$omp target update to(ielem_er2)
  end if
  if (allocated(ielem_er2dt)) then
    !$omp target enter data map(alloc: ielem_er2dt)
    !$omp target update to(ielem_er2dt)
  end if
  if (allocated(ielem_erx)) then
    !$omp target enter data map(alloc: ielem_erx)
    !$omp target update to(ielem_erx)
  end if
  if (allocated(ielem_faceanglex)) then
    !$omp target enter data map(alloc: ielem_faceanglex)
    !$omp target update to(ielem_faceanglex)
  end if
  if (allocated(ielem_faceangley)) then
    !$omp target enter data map(alloc: ielem_faceangley)
    !$omp target update to(ielem_faceangley)
  end if
  if (allocated(ielem_facediss)) then
    !$omp target enter data map(alloc: ielem_facediss)
    !$omp target update to(ielem_facediss)
  end if
  if (allocated(ielem_filtered)) then
    !$omp target enter data map(alloc: ielem_filtered)
    !$omp target update to(ielem_filtered)
  end if
  if (allocated(ielem_full)) then
    !$omp target enter data map(alloc: ielem_full)
    !$omp target update to(ielem_full)
  end if
  if (allocated(ielem_ggs)) then
    !$omp target enter data map(alloc: ielem_ggs)
    !$omp target update to(ielem_ggs)
  end if
  if (allocated(ielem_hybrid)) then
    !$omp target enter data map(alloc: ielem_hybrid)
    !$omp target update to(ielem_hybrid)
  end if
  if (allocated(ielem_ibounds)) then
    !$omp target enter data map(alloc: ielem_ibounds)
    !$omp target update to(ielem_ibounds)
  end if
  if (allocated(ielem_idegfree)) then
    !$omp target enter data map(alloc: ielem_idegfree)
    !$omp target update to(ielem_idegfree)
  end if
  if (allocated(ielem_ifca)) then
    !$omp target enter data map(alloc: ielem_ifca)
    !$omp target update to(ielem_ifca)
  end if
  if (allocated(ielem_ihex)) then
    !$omp target enter data map(alloc: ielem_ihex)
    !$omp target update to(ielem_ihex)
  end if
  if (allocated(ielem_ihexgl)) then
    !$omp target enter data map(alloc: ielem_ihexgl)
    !$omp target update to(ielem_ihexgl)
  end if
  if (allocated(ielem_indexf)) then
    !$omp target enter data map(alloc: ielem_indexf)
    !$omp target update to(ielem_indexf)
  end if
  if (allocated(ielem_indexi)) then
    !$omp target enter data map(alloc: ielem_indexi)
    !$omp target update to(ielem_indexi)
  end if
  if (allocated(ielem_ineigh)) then
    !$omp target enter data map(alloc: ielem_ineigh)
    !$omp target update to(ielem_ineigh)
  end if
  if (allocated(ielem_ineighb)) then
    !$omp target enter data map(alloc: ielem_ineighb)
    !$omp target update to(ielem_ineighb)
  end if
  if (allocated(ielem_ineighg)) then
    !$omp target enter data map(alloc: ielem_ineighg)
    !$omp target update to(ielem_ineighg)
  end if
  if (allocated(ielem_ineighn)) then
    !$omp target enter data map(alloc: ielem_ineighn)
    !$omp target update to(ielem_ineighn)
  end if
  if (allocated(ielem_inter_id)) then
    !$omp target enter data map(alloc: ielem_inter_id)
    !$omp target update to(ielem_inter_id)
  end if
  if (allocated(ielem_interior)) then
    !$omp target enter data map(alloc: ielem_interior)
    !$omp target update to(ielem_interior)
  end if
  if (allocated(ielem_inumneighbours)) then
    !$omp target enter data map(alloc: ielem_inumneighbours)
    !$omp target update to(ielem_inumneighbours)
  end if
  if (allocated(ielem_iorder)) then
    !$omp target enter data map(alloc: ielem_iorder)
    !$omp target update to(ielem_iorder)
  end if
  if (allocated(ielem_ishape)) then
    !$omp target enter data map(alloc: ielem_ishape)
    !$omp target update to(ielem_ishape)
  end if
  if (allocated(ielem_itotalpoints)) then
    !$omp target enter data map(alloc: ielem_itotalpoints)
    !$omp target update to(ielem_itotalpoints)
  end if
  if (allocated(ielem_linc)) then
    !$omp target enter data map(alloc: ielem_linc)
    !$omp target update to(ielem_linc)
  end if
  if (allocated(ielem_lwcx2)) then
    !$omp target enter data map(alloc: ielem_lwcx2)
    !$omp target update to(ielem_lwcx2)
  end if
  if (allocated(ielem_minedge)) then
    !$omp target enter data map(alloc: ielem_minedge)
    !$omp target update to(ielem_minedge)
  end if
  if (allocated(ielem_mode)) then
    !$omp target enter data map(alloc: ielem_mode)
    !$omp target update to(ielem_mode)
  end if
  if (allocated(ielem_mood)) then
    !$omp target enter data map(alloc: ielem_mood)
    !$omp target update to(ielem_mood)
  end if
  if (allocated(ielem_mood_o)) then
    !$omp target enter data map(alloc: ielem_mood_o)
    !$omp target update to(ielem_mood_o)
  end if
  if (allocated(ielem_nodes)) then
    !$omp target enter data map(alloc: ielem_nodes)
    !$omp target update to(ielem_nodes)
  end if
  if (allocated(ielem_nodes_faces)) then
    !$omp target enter data map(alloc: ielem_nodes_faces)
    !$omp target update to(ielem_nodes_faces)
  end if
  if (allocated(ielem_nodes_faces_v)) then
    !$omp target enter data map(alloc: ielem_nodes_faces_v)
    !$omp target update to(ielem_nodes_faces_v)
  end if
  if (allocated(ielem_nodes_neighbours)) then
    !$omp target enter data map(alloc: ielem_nodes_neighbours)
    !$omp target update to(ielem_nodes_neighbours)
  end if
  if (allocated(ielem_nodes_v)) then
    !$omp target enter data map(alloc: ielem_nodes_v)
    !$omp target update to(ielem_nodes_v)
  end if
  if (allocated(ielem_nofbc)) then
    !$omp target enter data map(alloc: ielem_nofbc)
    !$omp target update to(ielem_nofbc)
  end if
  if (allocated(ielem_nojecount)) then
    !$omp target enter data map(alloc: ielem_nojecount)
    !$omp target update to(ielem_nojecount)
  end if
  if (allocated(ielem_nonodes)) then
    !$omp target enter data map(alloc: ielem_nonodes)
    !$omp target update to(ielem_nonodes)
  end if
  if (allocated(ielem_q_face_q_mapl)) then
    !$omp target enter data map(alloc: ielem_q_face_q_mapl)
    !$omp target update to(ielem_q_face_q_mapl)
  end if
  if (allocated(ielem_qface)) then
    !$omp target enter data map(alloc: ielem_qface)
    !$omp target update to(ielem_qface)
  end if
  if (allocated(ielem_recalc)) then
    !$omp target enter data map(alloc: ielem_recalc)
    !$omp target update to(ielem_recalc)
  end if
  if (allocated(ielem_reduce)) then
    !$omp target enter data map(alloc: ielem_reduce)
    !$omp target update to(ielem_reduce)
  end if
  if (allocated(ielem_reorient)) then
    !$omp target enter data map(alloc: ielem_reorient)
    !$omp target update to(ielem_reorient)
  end if
  if (allocated(ielem_stencil_dist)) then
    !$omp target enter data map(alloc: ielem_stencil_dist)
    !$omp target update to(ielem_stencil_dist)
  end if
  if (allocated(ielem_surf)) then
    !$omp target enter data map(alloc: ielem_surf)
    !$omp target update to(ielem_surf)
  end if
  if (allocated(ielem_totvolume)) then
    !$omp target enter data map(alloc: ielem_totvolume)
    !$omp target update to(ielem_totvolume)
  end if
  if (allocated(ielem_troubled)) then
    !$omp target enter data map(alloc: ielem_troubled)
    !$omp target update to(ielem_troubled)
  end if
  if (allocated(ielem_types_faces)) then
    !$omp target enter data map(alloc: ielem_types_faces)
    !$omp target update to(ielem_types_faces)
  end if
  if (allocated(ielem_vdec)) then
    !$omp target enter data map(alloc: ielem_vdec)
    !$omp target update to(ielem_vdec)
  end if
  if (allocated(ielem_viscx)) then
    !$omp target enter data map(alloc: ielem_viscx)
    !$omp target update to(ielem_viscx)
  end if
  if (allocated(ielem_vortex)) then
    !$omp target enter data map(alloc: ielem_vortex)
    !$omp target update to(ielem_vortex)
  end if
  if (allocated(ielem_walldist)) then
    !$omp target enter data map(alloc: ielem_walldist)
    !$omp target update to(ielem_walldist)
  end if
  if (allocated(ielem_walls)) then
    !$omp target enter data map(alloc: ielem_walls)
    !$omp target update to(ielem_walls)
  end if
  if (allocated(ielem_wcx)) then
    !$omp target enter data map(alloc: ielem_wcx)
    !$omp target update to(ielem_wcx)
  end if
  if (allocated(ielem_xxc)) then
    !$omp target enter data map(alloc: ielem_xxc)
    !$omp target update to(ielem_xxc)
  end if
  if (allocated(ielem_yyc)) then
    !$omp target enter data map(alloc: ielem_yyc)
    !$omp target update to(ielem_yyc)
  end if
  if (allocated(ielem_zzc)) then
    !$omp target enter data map(alloc: ielem_zzc)
    !$omp target update to(ielem_zzc)
  end if
  if (allocated(impdiag_mf)) then
    !$omp target enter data map(alloc: impdiag_mf)
    !$omp target update to(impdiag_mf)
  end if
  if (allocated(impoff_mf)) then
    !$omp target enter data map(alloc: impoff_mf)
    !$omp target update to(impoff_mf)
  end if
  if (allocated(inoder4_bct)) then
    !$omp target enter data map(alloc: inoder4_bct)
    !$omp target update to(inoder4_bct)
  end if
  if (allocated(inoder4_cord)) then
    !$omp target enter data map(alloc: inoder4_cord)
    !$omp target update to(inoder4_cord)
  end if
  if (allocated(inoder4_itor)) then
    !$omp target enter data map(alloc: inoder4_itor)
    !$omp target update to(inoder4_itor)
  end if
  if (allocated(integ_basis_dg_value)) then
    !$omp target enter data map(alloc: integ_basis_dg_value)
    !$omp target update to(integ_basis_dg_value)
  end if
  if (allocated(integ_basis_value)) then
    !$omp target enter data map(alloc: integ_basis_value)
    !$omp target update to(integ_basis_value)
  end if
  if (allocated(integ_basis_valuec)) then
    !$omp target enter data map(alloc: integ_basis_valuec)
    !$omp target update to(integ_basis_valuec)
  end if
  if (allocated(m_1_val)) then
    !$omp target enter data map(alloc: m_1_val)
    !$omp target update to(m_1_val)
  end if
  if (allocated(modal_filter)) then
    !$omp target enter data map(alloc: modal_filter)
    !$omp target update to(modal_filter)
  end if
  if (allocated(modal_filter_strong)) then
    !$omp target enter data map(alloc: modal_filter_strong)
    !$omp target update to(modal_filter_strong)
  end if
  if (allocated(modal_filter_weak)) then
    !$omp target enter data map(alloc: modal_filter_weak)
    !$omp target update to(modal_filter_weak)
  end if
  if (allocated(mp_a_in)) then
    !$omp target enter data map(alloc: mp_a_in)
    !$omp target update to(mp_a_in)
  end if
  if (allocated(mp_brok_a)) then
    !$omp target enter data map(alloc: mp_brok_a)
    !$omp target update to(mp_brok_a)
  end if
  if (allocated(mp_brok_b)) then
    !$omp target enter data map(alloc: mp_brok_b)
    !$omp target update to(mp_brok_b)
  end if
  if (allocated(mp_brok_c)) then
    !$omp target enter data map(alloc: mp_brok_c)
    !$omp target update to(mp_brok_c)
  end if
  if (allocated(mp_janaf)) then
    !$omp target enter data map(alloc: mp_janaf)
    !$omp target update to(mp_janaf)
  end if
  if (allocated(mp_m)) then
    !$omp target enter data map(alloc: mp_m)
    !$omp target update to(mp_m)
  end if
  if (allocated(mp_pinf)) then
    !$omp target enter data map(alloc: mp_pinf)
    !$omp target update to(mp_pinf)
  end if
  if (allocated(mp_r_in)) then
    !$omp target enter data map(alloc: mp_r_in)
    !$omp target update to(mp_r_in)
  end if
  if (allocated(mp_thigh_in)) then
    !$omp target enter data map(alloc: mp_thigh_in)
    !$omp target update to(mp_thigh_in)
  end if
  if (allocated(mp_tlow_in)) then
    !$omp target enter data map(alloc: mp_tlow_in)
    !$omp target update to(mp_tlow_in)
  end if
  if (allocated(mp_tmid_in)) then
    !$omp target enter data map(alloc: mp_tmid_in)
    !$omp target update to(mp_tmid_in)
  end if
  if (allocated(mrf_rot_gl)) then
    !$omp target enter data map(alloc: mrf_rot_gl)
    !$omp target update to(mrf_rot_gl)
  end if
  if (allocated(point1_gl)) then
    !$omp target enter data map(alloc: point1_gl)
    !$omp target update to(point1_gl)
  end if
  if (allocated(point2_gl)) then
    !$omp target enter data map(alloc: point2_gl)
    !$omp target update to(point2_gl)
  end if
  if (allocated(qp_array_qp_weight)) then
    !$omp target enter data map(alloc: qp_array_qp_weight)
    !$omp target update to(qp_array_qp_weight)
  end if
  if (allocated(qp_array_x)) then
    !$omp target enter data map(alloc: qp_array_x)
    !$omp target update to(qp_array_x)
  end if
  if (allocated(qp_array_y)) then
    !$omp target enter data map(alloc: qp_array_y)
    !$omp target update to(qp_array_y)
  end if
  if (allocated(qp_array_z)) then
    !$omp target enter data map(alloc: qp_array_z)
    !$omp target update to(qp_array_z)
  end if
  if (allocated(radius_gl)) then
    !$omp target enter data map(alloc: radius_gl)
    !$omp target update to(radius_gl)
  end if
  if (allocated(rec_br2_aux_var)) then
    !$omp target enter data map(alloc: rec_br2_aux_var)
    !$omp target update to(rec_br2_aux_var)
  end if
  if (allocated(rec_br2_local_lift)) then
    !$omp target enter data map(alloc: rec_br2_local_lift)
    !$omp target update to(rec_br2_local_lift)
  end if
  if (allocated(rec_cgradientstemp)) then
    !$omp target enter data map(alloc: rec_cgradientstemp)
    !$omp target update to(rec_cgradientstemp)
  end if
  if (allocated(rec_cond)) then
    !$omp target enter data map(alloc: rec_cond)
    !$omp target update to(rec_cond)
  end if
  if (allocated(rec_findw)) then
    !$omp target enter data map(alloc: rec_findw)
    !$omp target update to(rec_findw)
  end if
  if (allocated(rec_g0)) then
    !$omp target enter data map(alloc: rec_g0)
    !$omp target update to(rec_g0)
  end if
  if (allocated(rec_gradf)) then
    !$omp target enter data map(alloc: rec_gradf)
    !$omp target update to(rec_gradf)
  end if
  if (allocated(rec_gradients)) then
    !$omp target enter data map(alloc: rec_gradients)
    !$omp target update to(rec_gradients)
  end if
  if (allocated(rec_gradients2)) then
    !$omp target enter data map(alloc: rec_gradients2)
    !$omp target update to(rec_gradients2)
  end if
  if (allocated(rec_gradientsc)) then
    !$omp target enter data map(alloc: rec_gradientsc)
    !$omp target update to(rec_gradientsc)
  end if
  if (allocated(rec_gradientsc2)) then
    !$omp target enter data map(alloc: rec_gradientsc2)
    !$omp target update to(rec_gradientsc2)
  end if
  if (allocated(rec_gradientstemp)) then
    !$omp target enter data map(alloc: rec_gradientstemp)
    !$omp target update to(rec_gradientstemp)
  end if
  if (allocated(rec_gradientstemp_wall)) then
    !$omp target enter data map(alloc: rec_gradientstemp_wall)
    !$omp target update to(rec_gradientstemp_wall)
  end if
  if (allocated(rec_gradientsturb)) then
    !$omp target enter data map(alloc: rec_gradientsturb)
    !$omp target update to(rec_gradientsturb)
  end if
  if (allocated(rec_gradientsturb_wall)) then
    !$omp target enter data map(alloc: rec_gradientsturb_wall)
    !$omp target update to(rec_gradientsturb_wall)
  end if
  if (allocated(rec_grads)) then
    !$omp target enter data map(alloc: rec_grads)
    !$omp target update to(rec_grads)
  end if
  if (allocated(rec_gradsav)) then
    !$omp target enter data map(alloc: rec_gradsav)
    !$omp target update to(rec_gradsav)
  end if
  if (allocated(rec_ihexb)) then
    !$omp target enter data map(alloc: rec_ihexb)
    !$omp target update to(rec_ihexb)
  end if
  if (allocated(rec_ihexbc)) then
    !$omp target enter data map(alloc: rec_ihexbc)
    !$omp target update to(rec_ihexbc)
  end if
  if (allocated(rec_ihexg)) then
    !$omp target enter data map(alloc: rec_ihexg)
    !$omp target update to(rec_ihexg)
  end if
  if (allocated(rec_ihexgc)) then
    !$omp target enter data map(alloc: rec_ihexgc)
    !$omp target update to(rec_ihexgc)
  end if
  if (allocated(rec_ihexl)) then
    !$omp target enter data map(alloc: rec_ihexl)
    !$omp target update to(rec_ihexl)
  end if
  if (allocated(rec_ihexlc)) then
    !$omp target enter data map(alloc: rec_ihexlc)
    !$omp target update to(rec_ihexlc)
  end if
  if (allocated(rec_ihexn)) then
    !$omp target enter data map(alloc: rec_ihexn)
    !$omp target update to(rec_ihexn)
  end if
  if (allocated(rec_ihexnc)) then
    !$omp target enter data map(alloc: rec_ihexnc)
    !$omp target update to(rec_ihexnc)
  end if
  if (allocated(rec_indicator)) then
    !$omp target enter data map(alloc: rec_indicator)
    !$omp target update to(rec_indicator)
  end if
  if (allocated(rec_indicatorc)) then
    !$omp target enter data map(alloc: rec_indicatorc)
    !$omp target update to(rec_indicatorc)
  end if
  if (allocated(rec_invccjac)) then
    !$omp target enter data map(alloc: rec_invccjac)
    !$omp target update to(rec_invccjac)
  end if
  if (allocated(rec_invctjac)) then
    !$omp target enter data map(alloc: rec_invctjac)
    !$omp target update to(rec_invctjac)
  end if
  if (allocated(rec_invmat_stencilt)) then
    !$omp target enter data map(alloc: rec_invmat_stencilt)
    !$omp target update to(rec_invmat_stencilt)
  end if
  if (allocated(rec_invmat_stenciltc)) then
    !$omp target enter data map(alloc: rec_invmat_stenciltc)
    !$omp target update to(rec_invmat_stenciltc)
  end if
  if (allocated(rec_k0)) then
    !$omp target enter data map(alloc: rec_k0)
    !$omp target update to(rec_k0)
  end if
  if (allocated(rec_local)) then
    !$omp target enter data map(alloc: rec_local)
    !$omp target update to(rec_local)
  end if
  if (allocated(rec_mrf)) then
    !$omp target enter data map(alloc: rec_mrf)
    !$omp target update to(rec_mrf)
  end if
  if (allocated(rec_mrf_origin)) then
    !$omp target enter data map(alloc: rec_mrf_origin)
    !$omp target update to(rec_mrf_origin)
  end if
  if (allocated(rec_mrf_velocity)) then
    !$omp target enter data map(alloc: rec_mrf_velocity)
    !$omp target update to(rec_mrf_velocity)
  end if
  if (allocated(rec_periodicflag)) then
    !$omp target enter data map(alloc: rec_periodicflag)
    !$omp target update to(rec_periodicflag)
  end if
  if (allocated(rec_qpoints)) then
    !$omp target enter data map(alloc: rec_qpoints)
    !$omp target update to(rec_qpoints)
  end if
  if (allocated(rec_rotvel)) then
    !$omp target enter data map(alloc: rec_rotvel)
    !$omp target update to(rec_rotvel)
  end if
  if (allocated(rec_rpoints)) then
    !$omp target enter data map(alloc: rec_rpoints)
    !$omp target update to(rec_rpoints)
  end if
  if (allocated(rec_stencils)) then
    !$omp target enter data map(alloc: rec_stencils)
    !$omp target update to(rec_stencils)
  end if
  if (allocated(rec_stencilsc)) then
    !$omp target enter data map(alloc: rec_stencilsc)
    !$omp target update to(rec_stencilsc)
  end if
  if (allocated(rec_surf_qpoints)) then
    !$omp target enter data map(alloc: rec_surf_qpoints)
    !$omp target update to(rec_surf_qpoints)
  end if
  if (allocated(rec_tempsq)) then
    !$omp target enter data map(alloc: rec_tempsq)
    !$omp target update to(rec_tempsq)
  end if
  if (allocated(rec_tempsqmat)) then
    !$omp target enter data map(alloc: rec_tempsqmat)
    !$omp target update to(rec_tempsqmat)
  end if
  if (allocated(rec_uleft)) then
    !$omp target enter data map(alloc: rec_uleft)
    !$omp target update to(rec_uleft)
  end if
  if (allocated(rec_uleft_dg)) then
    !$omp target enter data map(alloc: rec_uleft_dg)
    !$omp target update to(rec_uleft_dg)
  end if
  if (allocated(rec_uleftturb)) then
    !$omp target enter data map(alloc: rec_uleftturb)
    !$omp target update to(rec_uleftturb)
  end if
  if (allocated(rec_uleftturbv)) then
    !$omp target enter data map(alloc: rec_uleftturbv)
    !$omp target update to(rec_uleftturbv)
  end if
  if (allocated(rec_uleftv)) then
    !$omp target enter data map(alloc: rec_uleftv)
    !$omp target update to(rec_uleftv)
  end if
  if (allocated(rec_uleftx)) then
    !$omp target enter data map(alloc: rec_uleftx)
    !$omp target update to(rec_uleftx)
  end if
  if (allocated(rec_velinvlsqmat)) then
    !$omp target enter data map(alloc: rec_velinvlsqmat)
    !$omp target update to(rec_velinvlsqmat)
  end if
  if (allocated(rec_vellsq)) then
    !$omp target enter data map(alloc: rec_vellsq)
    !$omp target update to(rec_vellsq)
  end if
  if (allocated(rec_velocitydof_wall)) then
    !$omp target enter data map(alloc: rec_velocitydof_wall)
    !$omp target update to(rec_velocitydof_wall)
  end if
  if (allocated(rec_vext_ref)) then
    !$omp target enter data map(alloc: rec_vext_ref)
    !$omp target update to(rec_vext_ref)
  end if
  if (allocated(rec_volume)) then
    !$omp target enter data map(alloc: rec_volume)
    !$omp target update to(rec_volume)
  end if
  if (allocated(rec_volume_w)) then
    !$omp target enter data map(alloc: rec_volume_w)
    !$omp target update to(rec_volume_w)
  end if
  if (allocated(rec_volumec)) then
    !$omp target enter data map(alloc: rec_volumec)
    !$omp target update to(rec_volumec)
  end if
  if (allocated(rec_wall)) then
    !$omp target enter data map(alloc: rec_wall)
    !$omp target update to(rec_wall)
  end if
  if (allocated(rec_wallcoeff)) then
    !$omp target enter data map(alloc: rec_wallcoeff)
    !$omp target update to(rec_wallcoeff)
  end if
  if (allocated(rec_wallcoefg)) then
    !$omp target enter data map(alloc: rec_wallcoefg)
    !$omp target update to(rec_wallcoefg)
  end if
  if (allocated(rec_weightl)) then
    !$omp target enter data map(alloc: rec_weightl)
    !$omp target update to(rec_weightl)
  end if
  if (allocated(rec_weno)) then
    !$omp target enter data map(alloc: rec_weno)
    !$omp target update to(rec_weno)
  end if
  if (allocated(rec_weno2)) then
    !$omp target enter data map(alloc: rec_weno2)
    !$omp target update to(rec_weno2)
  end if
  if (allocated(rec_wenos)) then
    !$omp target enter data map(alloc: rec_wenos)
    !$omp target update to(rec_wenos)
  end if
  if (allocated(rg_hzero)) then
    !$omp target enter data map(alloc: rg_hzero)
    !$omp target update to(rg_hzero)
  end if
  if (allocated(rg_molm)) then
    !$omp target enter data map(alloc: rg_molm)
    !$omp target update to(rg_molm)
  end if
  if (allocated(rg_thetag)) then
    !$omp target enter data map(alloc: rg_thetag)
    !$omp target update to(rg_thetag)
  end if
  if (allocated(rg_tv_coef)) then
    !$omp target enter data map(alloc: rg_tv_coef)
    !$omp target update to(rg_tv_coef)
  end if
  if (allocated(rg_vf)) then
    !$omp target enter data map(alloc: rg_vf)
    !$omp target update to(rg_vf)
  end if
  if (allocated(rgs_ab)) then
    !$omp target enter data map(alloc: rgs_ab)
    !$omp target update to(rgs_ab)
  end if
  if (allocated(rgs_bb)) then
    !$omp target enter data map(alloc: rgs_bb)
    !$omp target update to(rgs_bb)
  end if
  if (allocated(rgs_cb)) then
    !$omp target enter data map(alloc: rgs_cb)
    !$omp target update to(rgs_cb)
  end if
  if (allocated(rgs_eps_over_k)) then
    !$omp target enter data map(alloc: rgs_eps_over_k)
    !$omp target update to(rgs_eps_over_k)
  end if
  if (allocated(rgs_mg)) then
    !$omp target enter data map(alloc: rgs_mg)
    !$omp target update to(rgs_mg)
  end if
  if (allocated(rgs_sigmaa)) then
    !$omp target enter data map(alloc: rgs_sigmaa)
    !$omp target update to(rgs_sigmaa)
  end if
  if (allocated(rhs_sol_mm_dg)) then
    !$omp target enter data map(alloc: rhs_sol_mm_dg)
    !$omp target update to(rhs_sol_mm_dg)
  end if
  if (allocated(rhs_val)) then
    !$omp target enter data map(alloc: rhs_val)
    !$omp target update to(rhs_val)
  end if
  if (allocated(rhs_valdg)) then
    !$omp target enter data map(alloc: rhs_valdg)
    !$omp target update to(rhs_valdg)
  end if
  if (allocated(rhst_val)) then
    !$omp target enter data map(alloc: rhst_val)
    !$omp target update to(rhst_val)
  end if
  if (allocated(sht_rg)) then
    !$omp target enter data map(alloc: sht_rg)
    !$omp target update to(sht_rg)
  end if
  if (allocated(u_c_br2_aux_var)) then
    !$omp target enter data map(alloc: u_c_br2_aux_var)
    !$omp target update to(u_c_br2_aux_var)
  end if
  if (allocated(u_c_rms)) then
    !$omp target enter data map(alloc: u_c_rms)
    !$omp target update to(u_c_rms)
  end if
  if (allocated(u_c_val)) then
    !$omp target enter data map(alloc: u_c_val)
    !$omp target update to(u_c_val)
  end if
  if (allocated(u_c_valdg)) then
    !$omp target enter data map(alloc: u_c_valdg)
    !$omp target update to(u_c_valdg)
  end if
  if (allocated(u_cs_val)) then
    !$omp target enter data map(alloc: u_cs_val)
    !$omp target update to(u_cs_val)
  end if
  if (allocated(u_cs_valdg)) then
    !$omp target enter data map(alloc: u_cs_valdg)
    !$omp target update to(u_cs_valdg)
  end if
  if (allocated(u_ct_val)) then
    !$omp target enter data map(alloc: u_ct_val)
    !$omp target update to(u_ct_val)
  end if
  if (allocated(u_cw_val)) then
    !$omp target enter data map(alloc: u_cw_val)
    !$omp target update to(u_cw_val)
  end if
  if (allocated(u_cw_valdg)) then
    !$omp target enter data map(alloc: u_cw_valdg)
    !$omp target update to(u_cw_valdg)
  end if
  if (allocated(u_e_val)) then
    !$omp target enter data map(alloc: u_e_val)
    !$omp target update to(u_e_val)
  end if

  if (allocated(dg2fv)) then
    !$omp target enter data map(alloc: dg2fv)
    !$omp target update to(dg2fv)
  end if
  if (allocated(impdiag)) then
    !$omp target enter data map(alloc: impdiag)
    !$omp target update to(impdiag)
  end if
  if (allocated(impdiagt)) then
    !$omp target enter data map(alloc: impdiagt)
    !$omp target update to(impdiagt)
  end if
  if (allocated(impdu)) then
    !$omp target enter data map(alloc: impdu)
    !$omp target update to(impdu)
  end if
  if (allocated(impoff)) then
    !$omp target enter data map(alloc: impoff)
    !$omp target update to(impoff)
  end if
  if (allocated(impofft)) then
    !$omp target enter data map(alloc: impofft)
    !$omp target update to(impofft)
  end if
  if (allocated(nodelist)) then
    !$omp target enter data map(alloc: nodelist)
    !$omp target update to(nodelist)
  end if
  if (allocated(sht)) then
    !$omp target enter data map(alloc: sht)
    !$omp target update to(sht)
  end if

  if (allocated(xmpielrank)) then
    !$omp target enter data map(alloc: xmpielrank)
    !$omp target update to(xmpielrank)
  end if



    !$omp target enter data map(alloc: ineedhalo, ineedhalos,ineedbound,ineedbounds,bound_total,bounds_total, halo_total,halos_total,ilength1,ilength2)
    !$omp target update to(ineedhalo, ineedhalos,ineedbound,ineedbounds,bound_total,bounds_total, halo_total,halos_total,ilength1,ilength2)

!----------------------------
! Integers
!----------------------------
if (allocated(halo_len)       .and. size(halo_len)       > 0) then
  !$omp target enter data map(alloc: halo_len)
  !$omp target update to(halo_len)
end if
if (allocated(halos_len)      .and. size(halos_len)      > 0) then
  !$omp target enter data map(alloc: halos_len)
  !$omp target update to(halos_len)
end if

if (allocated(halo_offset)    .and. size(halo_offset)    > 0) then
  !$omp target enter data map(alloc: halo_offset)
  !$omp target update to(halo_offset)
end if
if (allocated(halos_offset)   .and. size(halos_offset)   > 0) then
  !$omp target enter data map(alloc: halos_offset)
  !$omp target update to(halos_offset)
end if

if (allocated(bound_len)      .and. size(bound_len)      > 0) then
  !$omp target enter data map(alloc: bound_len)
  !$omp target update to(bound_len)
end if
if (allocated(bounds_len)     .and. size(bounds_len)     > 0) then
  !$omp target enter data map(alloc: bounds_len)
  !$omp target update to(bounds_len)
end if

if (allocated(need_side)      .and. size(need_side)      > 0) then
  !$omp target enter data map(alloc: need_side)
  !$omp target update to(need_side)
end if
if (allocated(need_q)         .and. size(need_q)         > 0) then
  !$omp target enter data map(alloc: need_q)
  !$omp target update to(need_q)
end if
if (allocated(need_loc)       .and. size(need_loc)       > 0) then
  !$omp target enter data map(alloc: need_loc)
  !$omp target update to(need_loc)
end if

if (allocated(halo_proc)      .and. size(halo_proc)      > 0) then
  !$omp target enter data map(alloc: halo_proc)
  !$omp target update to(halo_proc)
end if
if (allocated(halos_proc)     .and. size(halos_proc)     > 0) then
  !$omp target enter data map(alloc: halos_proc)
  !$omp target update to(halos_proc)
end if

if (allocated(bound_proc)     .and. size(bound_proc)     > 0) then
  !$omp target enter data map(alloc: bound_proc)
  !$omp target update to(bound_proc)
end if
if (allocated(bounds_proc)    .and. size(bounds_proc)    > 0) then
  !$omp target enter data map(alloc: bounds_proc)
  !$omp target update to(bounds_proc)
end if

if (allocated(bound_offset)   .and. size(bound_offset)   > 0) then
  !$omp target enter data map(alloc: bound_offset)
  !$omp target update to(bound_offset)
end if
if (allocated(bounds_offset)  .and. size(bounds_offset)  > 0) then
  !$omp target enter data map(alloc: bounds_offset)
  !$omp target update to(bounds_offset)
end if


if (allocated(bounds_offset)  .and. size(bounds_offset)  > 0) then
  !$omp target enter data map(alloc: bounds_offset)
  !$omp target update to(bounds_offset)
end if

!----------------------------
! Reals (2D halo arrays)
!----------------------------
if (allocated(solhir)  .and. size(solhir)  > 0) then
  !$omp target enter data map(alloc: solhir)
   !$omp target update to(solhir)
end if
if (allocated(solhis)  .and. size(solhis)  > 0) then
  !$omp target enter data map(alloc: solhis)
  !$omp target update to(solhis)
end if

if (allocated(solhird) .and. size(solhird) > 0) then
  !$omp target enter data map(alloc: solhird)
  !$omp target update to(solhird)
end if
if (allocated(solhisd) .and. size(solhisd) > 0) then
  !$omp target enter data map(alloc: solhisd)
  !$omp target update to(solhisd)
end if


!----------------------------
! Reals (boundary arrays)
!----------------------------
if (allocated(boundhiri) .and. size(boundhiri) > 0) then
  !$omp target enter data map(alloc: boundhiri)
  !$omp target update to(boundhiri)
end if

if (allocated(boundhisi) .and. size(boundhisi) > 0) then
  !$omp target enter data map(alloc: boundhisi)
  !$omp target update to(boundhisi)
end if

if (allocated(boundhir) .and. size(boundhir) > 0) then
  !$omp target enter data map(alloc: boundhir)
  !$omp target update to(boundhir)
end if

if (allocated(boundhis) .and. size(boundhis) > 0) then
  !$omp target enter data map(alloc: boundhis)
  !$omp target update to(boundhis)
end if

if (allocated(boundhir_dg) .and. size(boundhir_dg) > 0) then
  !$omp target enter data map(alloc: boundhir_dg)
  !$omp target update to(boundhir_dg)
end if

if (allocated(boundhis_dg) .and. size(boundhis_dg) > 0) then
  !$omp target enter data map(alloc: boundhis_dg)
  !$omp target update to(boundhis_dg)
end if

if (allocated(boundhirm) .and. size(boundhirm) > 0) then
  !$omp target enter data map(alloc: boundhirm)
  !$omp target update to(boundhirm)
end if

if (allocated(boundhism) .and. size(boundhism) > 0) then
  !$omp target enter data map(alloc: boundhism)
  !$omp target update to(boundhism)
end if

if (allocated(solhi_loc) .and. size(solhi_loc) > 0) then
  !$omp target enter data map(alloc: solhi_loc)
  !$omp target update to(solhi_loc)
end if

if (allocated(solhir_flat) .and. size(solhir_flat) > 0) then
  !$omp target enter data map(alloc: solhir_flat)
  !$omp target update to(solhir_flat)
end if

if (allocated(solhis_flat) .and. size(solhis_flat) > 0) then
  !$omp target enter data map(alloc: solhis_flat)
  !$omp target update to(solhis_flat)
end if

if (allocated(boundhir_flat) .and. size(boundhir_flat) > 0) then
  !$omp target enter data map(alloc: boundhir_flat)
  !$omp target update to(boundhir_flat)
end if

if (allocated(boundhis_flat) .and. size(boundhis_flat) > 0) then
  !$omp target enter data map(alloc: boundhis_flat)
  !$omp target update to(boundhis_flat)
end if

if (allocated(boundhir_dgflat) .and. size(boundhir_dgflat) > 0) then
  !$omp target enter data map(alloc: boundhir_dgflat)
  !$omp target update to(boundhir_dgflat)
end if

if (allocated(boundhis_dgflat) .and. size(boundhis_dgflat) > 0) then
  !$omp target enter data map(alloc: boundhis_dgflat)
  !$omp target update to(boundhis_dgflat)
end if

if (allocated(boundhiri_flat) .and. size(boundhiri_flat) > 0) then
  !$omp target enter data map(alloc: boundhiri_flat)
  !$omp target update to(boundhiri_flat)
end if

if (allocated(boundhisi_flat) .and. size(boundhisi_flat) > 0) then
  !$omp target enter data map(alloc: boundhisi_flat)
  !$omp target update to(boundhisi_flat)
end if

if (allocated(weights_l) .and. size(weights_l) > 0) then
  !$omp target enter data map(to: weights_l)
end if
if (allocated(weights_t) .and. size(weights_t) > 0) then
  !$omp target enter data map(to: weights_t)
end if

if (allocated(weights_q) .and. size(weights_q) > 0) then
  !$omp target enter data map(to: weights_q)
end if





#endif
end subroutine omp_map_first



end module memory
