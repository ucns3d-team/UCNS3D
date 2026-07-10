module prestore
use declaration
use derivatives
use library
use basis
use transform
use local
use lapck
use dg_functions
implicit none


 contains
 
 
 subroutine prestore_1(n)
!> @brief
!> this subroutine calls other subroutines to prestore the pseudoinverse reconstruction least square matrices
implicit none
integer,intent(in)::n
integer::kmaxe,i,m,k
integer,allocatable,dimension(:,:)::ilox_ihexg  !global index of cells
integer,allocatable,dimension(:,:)::ilox_ihexl  !local index of cells
integer,allocatable,dimension(:,:)::ilox_ihexb  !cpu that that each cell belongs to
integer,allocatable,dimension(:,:)::ilox_ihexn  !internal index from where to take the values from communicated messages
integer,allocatable,dimension(:,:)::ilox_ishape !shape of each element
real,allocatable,dimension(:,:)::ilox_xxc       !cell centre coordinates in x
real,allocatable,dimension(:,:)::ilox_yyc       !cell centre coordinates in y
real,allocatable,dimension(:,:)::ilox_zzc      !cell centre coordinates in z
real,allocatable,dimension(:,:)::ilox_volume    !cell volume
integer,allocatable,dimension(:,:)::ilox_periodicflag
integer,allocatable,dimension(:,:,:)::ilon_nodcount  !number of nodes
real,allocatable,dimension(:,:,:)::ilon_x           !coordinates of each node in x
real,allocatable,dimension(:,:,:)::ilon_y           !coordinates of each node in x
real,allocatable,dimension(:,:,:)::ilon_z           !coordinates of each node in x



kmaxe=xmpielrank(n)
m=typesten




	allocate (ilox_ihexg(m,numneighbours*iextend))
	allocate (ilox_ihexl(m,numneighbours*iextend))
	allocate (ilox_ihexb(m,numneighbours*iextend))
	allocate (ilox_ihexn(m,numneighbours*iextend))
	allocate (ilox_periodicflag(m,numneighbours*iextend))
	allocate (ilox_ishape(m,numneighbours*iextend))
	allocate (ilox_xxc(m,numneighbours*iextend))
	allocate (ilox_volume(m,numneighbours*iextend))
	allocate (ilox_yyc(m,numneighbours*iextend))
	if (dimensiona.eq.3)then
	allocate (ilox_zzc(m,numneighbours*iextend))
	ilox_zzc=0.d0
	end if
	ilox_ihexg=0
	ilox_ihexl=0
	ilox_ihexb=0
	ilox_ishape=0
	ilox_ihexn=0
	ilox_volume=0.d0
	ilox_xxc=0.d0
	ilox_yyc=0.d0


if (dimensiona.eq.3)then
  k=8
else
  k=4
end if





	allocate (ilon_nodcount(m,numneighbours*iextend,k))
	allocate (ilon_x(m,numneighbours*iextend,k))
	allocate (ilon_y(m,numneighbours*iextend,k))
	if (dimensiona.eq.3)then
	allocate (ilon_z(m,numneighbours*iextend,k))
	ilon_z=0.d0
	end if
	ilon_nodcount=0
	ilon_x=0.d0
	ilon_y=0.d0






kmaxe=xmpielrank(n)

if (dimensiona.eq.3)then
	!!$omp do
do i=1,kmaxe

call localise_stencil(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)


call localise_sten2(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
call checkgrads(n,i)
call find_rot_angles(n,i)

call prestore_reconstruction3(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)

 if ((dg.eq.1).or.(adda_type.eq.2))then
 
 call prestore_dg1(i)
    end if
 

end do
	!!$omp end do




else





	!!$omp do
do i=1,kmaxe
call localise_stencil2d(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
call localise_sten2d(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
 call checkgrads2d(n,i)
call find_rot_angles2d(n,i)
 call prestore_reconstruction2(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
 if ((dg.eq.1).or.(adda_type.eq.2))then
 
 call prestore_dg1(i)
 
 
    end if
end do
	!!$omp end do






end if


deallocate(ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y)

	if (dimensiona.eq.3)then
	deallocate(ilox_zzc,ilon_z)
	end if




end subroutine prestore_1






subroutine prestore_reconstruction3(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine prestores the pseudoinverse reconstruction least square matrices
implicit none
integer :: kbasis
integer :: mm_i
integer :: mm_j
integer :: mm_k
real :: mm_old(1:3)
integer,intent(in)::n,iconsi
integer::i,j,k,llco,ll,ii,igf,igf2,ifd2,idum,idum2,iq,jq,lq,ihgt,ihgj,iqp,iqp2,nnd,k0,g0,lcou,lcc,iqqq,icond1,icond2,n_node
integer::ideg,ideg2,imax,imax2,ivgt,jxx,ixx,lxx1,kxx,icompwrt,number_of_dog,eltype,inumo,inumo2,inum,facex,iconsidered
real::ssss,gggg,uptemp,lotemp,x_stencil,y_stencil,z_stencil,dist_sten,dist_sten2
real::ax,ay,az,angle1,angle2,nx,ny,nz,nnx,nnz,nny,x1,y1,z1
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dims,1:dims)::ainvjt
real,dimension(1:dimensiona)::cords
real,allocatable,dimension(:)::intbs,basefaceval,basefacgval,permutation,permutationg,xder,yder,zder
real,allocatable,dimension(:,:)::stencil,invmat,wlsqr,lsqm,lscqm,qff,rff,qtff,invrff,vellsqmat
real,allocatable,dimension(:,:,:)::stencils
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



allocate(lscqm(1:idegfree,1:idegfree),qff(1:idegfree,1:idegfree),rff(1:idegfree,1:idegfree),qtff(1:idegfree,1:idegfree),invrff(1:idegfree,1:idegfree))
allocate(vellsqmat(1:idegfree-1,1:idegfree-1))
allocate(lsqm(1:imaxdegfree,1:idegfree-1))	
allocate(intbs(1:idegfree),basefaceval(1:idegfree),basefacgval(1:idegfree),permutation(1:idegfree),permutationg(1:idegfree),xder(1:idegfree),yder(1:idegfree),zder(1:idegfree))
allocate(invmat(1:idegfree,1:idegfree))
allocate(wlsqr(1:20,1:numneighbours-1))
allocate(stencil(1:numneighbours-1,1:idegfree))
allocate(stencils(1:20,1:numneighbours-1,1:idegfree))	



i=iconsi


idum=0;
                if (ielem_interior(i).eq.1)then
                        do j=1,ielem_ifca(i)
                        if (ielem_ibounds(j,i).gt.0)then
                            if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
                                idum=1
                            end if
                        end if
                        end do
                end if
                

iconsidered=i


		
        !firstly compute the integrals of the basis functions
		
		intbs=zero;jxx=1;ixx=i;lxx1=1;number_of_dog=ielem_idegfree(i);kxx=ielem_iorder(i);eltype=ielem_ishape(i)
		icompwrt=0

		intbs=calintbasis(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,eltype,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)

        ! store integration-basis vector intbs into all entries 1:number_of_dog for element i

            integ_basis_value(1:number_of_dog,i) = intbs(1:number_of_dog)

                
        !the indicator matrix for the smaller polynomial
		if (iweno.eq.1)then
		call indicatormatrix(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
		end if
		
					!the indicator matrix for the lower order polynomial
					if (ees.eq.5)then
					intbs=zero;jxx=1;ixx=i;lxx1=1;number_of_dog=idegfree2;kxx=iorder2;eltype=ielem_ishape(i)
					icompwrt=1
					intbs=calintbasis(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,eltype,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
                    ! store lower-order integration-basis vector

                        integ_basis_valuec(1:number_of_dog,i) = intbs(1:number_of_dog)

						if (iweno.eq.1)then
						call indicatormatrix2(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
						icompwrt=0
						end if
					end if

		
		llco=ielem_admis(i);imax=ielem_inumneighbours(i)-1;inum=ielem_inumneighbours(i);ideg=ielem_idegfree(i)
	   inumo=ielem_iorder(i);imax2=numneighbours2-1;inum2=numneighbours2;ideg2=idegfree2;inumo2=iorder2


	   !this section computes the ratio of volumes inside each stencil
                dist_sten=zero
				
				dist_sten=-tolbig; dist_sten2=tolbig
	       	do ll=1,ielem_inumneighbours(i)
				if (ilox_volume(1,ll).lt.dist_sten2)then
				dist_sten2=ilox_volume(1,ll)
				end if
				if (ilox_volume(1,ll).gt.dist_sten)then
				dist_sten=ilox_volume(1,ll)
				end if
			end do
			
! 			ielem_stencil_dist(i)=max((dist_sten/dist_sten2),(dist_sten2/dist_sten))
				
		


		!now we start with looping all the admissible stencils
		do ll=1,llco	!for all stencils
									if((ees.ne.5).or.(ll.eq.1))then
										imax=ielem_inumneighbours(i)-1;inum=ielem_inumneighbours(i);ideg=ielem_idegfree(i);inumo=ielem_iorder(i)
										icompwrt=0;number_of_dog=ideg
									else
										imax=numneighbours2-1;inum=numneighbours2;ideg=idegfree2;inumo=iorder2;number_of_dog=ideg
										icompwrt=1
									end if

				do k=1,imax	!for all neighbours
								ixx=i; kxx=inumo


								x_stencil=(ilox_xxc(ll,k+1)-ilox_xxc(ll,1))**2
								y_stencil=(ilox_yyc(ll,k+1)-ilox_yyc(ll,1))**2

								z_stencil=(ilox_zzc(ll,k+1)-ilox_zzc(ll,1))**2

								dist_sten2=sqrt(x_stencil+y_stencil+z_stencil)

									if (weight_lsqr.eq.1)then
									wlsqr(ll,k)=1.0d0/((dist_sten2))
									else
									wlsqr(ll,k)=1.0d0
									end if






					if (fastest.eq.1)then	!this is when transformation is not active (rarely used)
					x1 = ilox_xxc(ll,k+1)-ilox_xxc(ll,1)
					y1 = ilox_yyc(ll,k+1)-ilox_yyc(ll,1)
					z1 = ilox_zzc(ll,k+1)-ilox_zzc(ll,1)


					if (greengo.eq.0)then	!for the least squares green gauss gradients
							if (idum.eq.1)then
								if((ees.ne.5).or.(ll.eq.1))then
									icompwrt=0
									if (itestcase.eq.4)then
									rec_stencils(ll,k,1:ielem_idegfree(i),rec_wall(i))=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,ielem_iorder(i),i,ielem_idegfree(i),icompwrt)
									rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
									stencils(ll,k,1:ideg)=rec_stencils(ll,k,1:ielem_idegfree(i),rec_wall(i))
									else
									stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,ielem_iorder(i),i,ielem_idegfree(i),icompwrt)
									end if

								else
									icompwrt=1
								if (itestcase.eq.4)then
! 								rec_stencilsc(ll,k,1:ideg,rec_wall(i))=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,i,ideg,icompwrt)
! 								rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
 								stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,i,ideg,icompwrt)
								else
								stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,i,ideg,icompwrt)
								end if
								end if
							else
								if((ees.ne.5).or.(ll.eq.1))then

								icompwrt=0
								stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,ixx,ideg,icompwrt)
								else
								icompwrt=1
								stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,ixx,ideg,icompwrt)

								end if
							end if

					else
							if((ees.ne.5).or.(ll.eq.1))then
							icompwrt=0
							stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,ixx,ideg,icompwrt)
							else
							icompwrt=1
							stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec(n,x1,y1,z1,inumo,ixx,ideg,icompwrt)


							end if
					end if



					else		!with coordinate transformation


					ixx=i;jxx=k+1;lxx1=ll
					eltype=ilox_ishape(ll,k+1)
						if (greengo.eq.0)then
								if (idum.eq.1)then	!smaller memory footprint for non boundary elements
										if((ees.ne.5).or.(ll.eq.1))then
										icompwrt=0



										if (itestcase.eq.4)then
										rec_stencils(ll,k,1:ideg,rec_wall(i))=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
										rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
										stencils(ll,k,1:ideg)=rec_stencils(ll,k,1:ideg,rec_wall(i))
										else
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)

										end if



										else
										icompwrt=1
										if (itestcase.eq.4)then
! 										rec_stencilsc(ll,k,1:ideg,rec_wall(i))=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
! 										rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
										else
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
										end if

										end if
								else
										if((ees.ne.5).or.(ll.eq.1))then
										icompwrt=0
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
										else
										icompwrt=1
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)

										end if
								end if
						else
								if((ees.ne.5).or.(ll.eq.1))then
							icompwrt=0
								stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
								else
								icompwrt=1
								stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)

								end if

						end if


					end if




		end do	!close loop for all neighbours








		if ((greengo.eq.0))then
				if (idum.eq.1)then
					invmat=zero;lscqm=zero
					do iq=1,ideg;do jq=1,ideg;do lq=1,imax
					if((ees.ne.5).or.(ll.eq.1))then
					lscqm(jq,iq)=lscqm(jq,iq)&
					+((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
	!				lscqm(jq,iq)=lscqm(jq,iq)&
	!				+((rec_stencils(ll,lq,jq,rec_wall(i))*rec_stencils(ll,lq,iq,rec_wall(i))))
					else
					lscqm(jq,iq)=lscqm(jq,iq)&
					+((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
				!	lscqm(jq,iq)=lscqm(jq,iq)&
				!	+((rec_stencilsc(ll,lq,jq,rec_wall(i))*rec_stencilsc(ll,lq,iq,rec_wall(i))))
					end if
					end do;end do;end do
				else
					invmat=zero;lscqm=zero
					do iq=1,ideg;do jq=1,ideg;do lq=1,imax
					lscqm(jq,iq)=lscqm(jq,iq)&
					+((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
					end do;end do;end do
				end if
		else

				invmat=zero;lscqm=zero
				do iq=1,ideg;do jq=1,ideg;do lq=1,imax
				lscqm(jq,iq)=lscqm(jq,iq)&
				+((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
				end do;end do;end do

		end if




		qff(:,:)=zero; rff(:,:)=zero; qtff(:,:)=zero; rff(:,:)=zero;  invrff(:,:)=zero
		call qrdecomposition(lscqm,qff,rff,ideg)
		call transposematrix(qff,qtff,ideg)
		ivgt=ideg+1
		call invert(rff,invrff,ivgt)
		  do mm_i=1,ideg
		    do mm_j=1,ideg
		      invmat(mm_i,mm_j)=zero
		      do mm_k=1,ideg
		        invmat(mm_i,mm_j)=invmat(mm_i,mm_j)+invrff(mm_i,mm_k)*qtff(mm_k,mm_j)
		      end do
		    end do
		  end do


		if (greengo.eq.0)then
			if (idum.eq.1)then
				if((ees.ne.5).or.(ll.eq.1))then
				stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)!rec_stencils(ll,1:imax,1:ideg,rec_wall(i))
				else
				stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)!rec_stencilsc(ll,1:imax,1:ideg,rec_wall(i))
				end if
			else
			stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)
			end if
		else
			stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)
		end if

		if((ees.ne.5).or.(ll.eq.1))then

! 		call dgemm ('n','t',ideg,imax,ideg,alpha,invmat(1:ideg,1:ideg),ideg,&
! 		stencil(1:imax,1:ideg),imax,beta,rec_invmat_stencilt(1:ideg,1:imax,ll,i),ideg)


		  do mm_i=1,ideg
		    do mm_j=1,imax
		      rec_invmat_stencilt(mm_i,mm_j,ll,i)=zero
		      do mm_k=1,ideg
		        rec_invmat_stencilt(mm_i,mm_j,ll,i)=rec_invmat_stencilt(mm_i,mm_j,ll,i)+invmat(mm_i,mm_k)*stencil(mm_j,mm_k)
		      end do
		    end do
		  end do



					do iq=1,imax
					rec_invmat_stencilt(:,iq,ll,i)=rec_invmat_stencilt(:,iq,ll,i)&
					*ilox_volume(ll,iq+1)*wlsqr(ll,iq)
					end do


					if (rec_invmat_stencilt(1,1,ll,i).ne.rec_invmat_stencilt(1,1,ll,i))then	!prevents non-invertible nan stencils from being deployed

					ielem_full(i)=0
					end if



		else
! 					call dgemm ('n','t',ideg,imax,ideg,alpha,invmat(1:ideg,1:ideg),ideg,&
! 				stencil(1:imax,1:ideg),imax,beta,rec_invmat_stenciltc(1:ideg,1:imax,ll,i),ideg)


				  do mm_i=1,ideg
				    do mm_j=1,imax
				      rec_invmat_stenciltc(mm_i,mm_j,ll,i)=zero
				      do mm_k=1,ideg
				        rec_invmat_stenciltc(mm_i,mm_j,ll,i)=rec_invmat_stenciltc(mm_i,mm_j,ll,i)+invmat(mm_i,mm_k)*stencil(mm_j,mm_k)
				      end do
				    end do
				  end do


					if (rec_invmat_stenciltc(1,1,ll,i).ne.rec_invmat_stenciltc(1,1,ll,i))then	!!prevents non-invertible nan stencils from being deployed

					ielem_full(i)=0
					end if


					do iq=1,imax
					rec_invmat_stenciltc(:,iq,ll,i)=rec_invmat_stenciltc(:,iq,ll,i)&
					*ilox_volume(ll,iq+1)*wlsqr(ll,iq)
					end do
		end if


	
		if (ielem_ggs(i).ne.1)then	!ggs

		if (itestcase.eq.4)then	!test

	 

	 

if (ll.eq.1)then		!stencils
! 

  
  if (idum.eq.1)then		!for wall only
    
	 
	
!     	
    do ihgt=1,3; do ihgj=1,3
    ainvjt(ihgt,ihgj)=rec_invccjac(ihgj,ihgt,i)
   end do; end do

	 idum2=0
	 
	 
				basefaceval=zero
				basefacgval=zero
				permutation=zero
				permutationg=zero
			
	 
	 
	 
	 do j=1,ielem_ifca(i)		!for all faces
	    if (ielem_ibounds(j,i).gt.0)then		!for bounded only
		      if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then		!!for bounded only 2                         
					facex=j;iconsidered=i
								  call coordinates_face_inner(n,iconsidered,facex,vext,nodes_list)


								   if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if



								    cords(1:3)=zero
								    call cordinates3(n,nodes_list,n_node,cords(1:3))
							    
								    ay=cords(2)
								    ax=cords(1)
								    az=cords(3)
				
					
					    vext(1,1)=ax;vext(1,2)=ay;vext(1,3)=az
					      do mm_j=1,3
					        mm_old(mm_j)=vext(1,mm_j)-rec_vext_ref(mm_j,i)
					      end do
					      do mm_i=1,3
					        vext(1,mm_i)=zero
					        do mm_j=1,3
					          vext(1,mm_i)=vext(1,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
					        end do
					      end do
					  
					    ax=vext(1,1);ay=vext(1,2);az=vext(1,3)
				
				
				
				
				
				
				angle1=ielem_faceanglex(j,i)
				angle2=ielem_faceangley(j,i)
				
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
				nnx=(nx*ainvjt(1,1))+(ny*ainvjt(2,1))+(nz*ainvjt(3,1))
				nny=(nx*ainvjt(1,2))+(ny*ainvjt(2,2))+(nz*ainvjt(3,2))
				nnz=(nx*ainvjt(1,3))+(ny*ainvjt(2,3))+(nz*ainvjt(3,3))
				
				do iq=1, ideg
				if (poly.eq.1) then
				xder(iq)=dfx(ax,ay,az,iq,i);  yder(iq)=dfy(ax,ay,az,iq,i);  zder(iq)=dfz(ax,ay,az,iq,i)
				end if
				if (poly.eq.2) then
				xder(iq)=dlx(ax,ay,az,iq,i);  yder(iq)=dly(ax,ay,az,iq,i);  zder(iq)=dlz(ax,ay,az,iq,i)
				end if
				if (poly.eq.4) then
				xder(iq)=tl3dx(ax,ay,az,iq,i);  yder(iq)=tl3dy(ax,ay,az,iq,i);  zder(iq)=tl3dz(ax,ay,az,iq,i)
				end if
				
				end do

				icompwrt=0

				basefaceval(1:ielem_idegfree(i))=basis_rec(n,ax,ay,az,ielem_iorder(i),i,ielem_idegfree(i),icompwrt)

				!if (thermal.eq.0)then
				basefacgval(1:ielem_idegfree(i))=((nnx*xder(1:ielem_idegfree(i)))+(nny*yder(1:ielem_idegfree(i)))+(nnz*zder(1:ielem_idegfree(i))))
				!else
				!icompwrt=0
				!basefacgval(1:ielem_idegfree(i))=basis_rec(n,ax,ay,az,ielem_iorder(i),i,ielem_idegfree(i),icompwrt)

! 				end if


				do iq=1,ideg
				rec_wallcoeff(iq,rec_wall(i))=basefaceval(iq)
				rec_wallcoefg(iq,rec_wall(i))=basefacgval(iq)
				permutation(iq)=iq
				permutationg(iq)=iq
				end do
		
				
		
		
		
		
		
		end if!for bounded only 2        
		end if!for bounded only
		end do!for all faces
		
		gggg=-tolbig
g0=0
do iq=1,ideg
if (abs(basefacgval(iq)) >gggg)then
gggg=abs(basefacgval(iq))
g0=iq
end if
end do




ssss=-tolbig
k0=0
do iq=1,ideg
if (abs(basefaceval(iq)) >ssss)then
ssss=abs(basefaceval(iq))
k0=iq
end if
end do

!  

! now swap basis functions and thus coefficients
permutation(1)=k0; permutation(k0)=1; ssss=basefaceval(1)
basefaceval(1)=basefaceval(k0); basefaceval(k0)=ssss
permutationg(1)=g0; permutation(g0)=1; gggg=basefacgval(1)
basefacgval(1)=basefacgval(g0); basefacgval(g0)=gggg
rec_k0(rec_wall(i))=k0
rec_g0(rec_wall(i))=g0
		
		
		
		
lsqm = zero
do lq=1,imax
lcou=0
    do iq=1,ideg
	if (iq.eq.g0) cycle
	    lcou=lcou+1
lsqm(lq,lcou)=rec_stencils(ll,lq,iq,rec_wall(i))&
-rec_stencils(ll,lq,g0,rec_wall(i))*rec_wallcoefg(iq,rec_wall(i))/rec_wallcoefg(g0,rec_wall(i))
    end do
end do
rec_tempsq(1:imax,1:ideg-1,rec_wall(i))=lsqm(1:imax,1:ideg-1)





vellsqmat=zero
do iq=1,ideg-1; do jq=1,ideg-1;	do lcc=1,imax
!now store the least square matrix
vellsqmat(jq,iq)= vellsqmat(jq,iq)+(lsqm(lcc,jq)*lsqm(lcc,iq))
end do;	end do;	end do

lscqm=zero
lscqm(1:ideg-1,1:ideg-1)=vellsqmat(1:ideg-1,1:ideg-1)
qff=zero; rff=zero; qtff=zero; invrff=zero
	  ivgt=ideg
	  
    call qrdecomposition(lscqm,qff,rff,ivgt-1)
    call transposematrix(qff,qtff,ivgt-1)
    
    call invert(rff,invrff,ivgt)
!   final inverted r^(-1)*q^(-1)
  do mm_i=1,ideg-1
    do mm_j=1,ideg-1
      rec_tempsqmat(mm_i,mm_j,rec_wall(i))=zero
      do mm_k=1,ideg-1
        rec_tempsqmat(mm_i,mm_j,rec_wall(i))=rec_tempsqmat(mm_i,mm_j,rec_wall(i))+invrff(mm_i,mm_k)*qtff(mm_k,mm_j)
      end do
    end do
  end do

! definy the matrix defining the least-square reconstruction in this case for velocity
lsqm = zero
do lq=1,imax
lcou=0
    do iq=1,ideg
	if (iq.eq.k0) cycle
	    lcou=lcou+1
	     lsqm(lq,lcou)=rec_stencils(ll,lq,iq,rec_wall(i))&
 -rec_stencils(ll,lq,k0,rec_wall(i))*rec_wallcoeff(iq,rec_wall(i))/rec_wallcoeff(k0,rec_wall(i))
    end do
end do
rec_vellsq(1:imax,1:ideg-1,rec_wall(i))=lsqm(1:imax,1:ideg-1)





		vellsqmat=zero
		do iq=1,ideg-1; do jq=1,ideg-1;	do lcc=1,imax
		!now store the least square matrix
		vellsqmat(jq,iq)= vellsqmat(jq,iq)+(lsqm(lcc,jq)*lsqm(lcc,iq))
		end do;	end do;	end do

		lscqm=zero
		lscqm(1:ideg-1,1:ideg-1)=vellsqmat(1:ideg-1,1:ideg-1)

		qff=zero; rff=zero; qtff=zero; invrff=zero
			ivgt=ideg
			call qrdecomposition(lscqm,qff,rff,ivgt-1)
			call transposematrix(qff,qtff,ivgt-1)

			call invert(rff,invrff,ivgt)
		!   final inverted r^(-1)*q^(-1)
		  do mm_i=1,ideg-1
		    do mm_j=1,ideg-1
		      rec_velinvlsqmat(mm_i,mm_j,rec_wall(i))=zero
		      do mm_k=1,ideg-1
		        rec_velinvlsqmat(mm_i,mm_j,rec_wall(i))=rec_velinvlsqmat(mm_i,mm_j,rec_wall(i))+invrff(mm_i,mm_k)*qtff(mm_k,mm_j)
		      end do
		    end do
		  end do
			
		
		
		
		
		
		
		end if!for wall only
		end if!stencils
		end if!for test
		end if!ggs
		end do !for all stencils
		

		deallocate(lsqm,qff,rff,qtff,invrff)	
		deallocate(vellsqmat)
		deallocate(lscqm)	
		deallocate(intbs,basefaceval,basefacgval,permutation,permutationg,xder,yder,zder)	
		deallocate(wlsqr)	
		deallocate(stencil)	
		deallocate(invmat)	
		deallocate(stencils)	




end subroutine prestore_reconstruction3





subroutine prestore_reconstruction2(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine prestores the pseudoinverse reconstruction least square matrices in 2d
implicit none
integer :: kbasis
integer :: mm_i
integer :: mm_j
integer :: mm_k
real :: mm_old(1:2)
integer,intent(in)::n,iconsi
integer::i,j,k,llco,ll,ii,igf,igf2,ifd2,idum,idum2,iq,jq,lq,ihgt,ihgj,iqp,iqp2,nnd,k0,g0,lcou,lcc,iqqq,icond1,icond2,n_node
integer::ideg,ideg2,imax,imax2,ivgt,jxx,ixx,lxx1,kxx,icompwrt,number_of_dog,eltype,inumo,inumo2,iconsidered,facex,inum,ai,aj
real::ssss,gggg,uptemp,lotemp,x_stencil,y_stencil,z_stencil,dist_sten,dist_sten2
real::ax,ay,az,angle1,angle2,nx,ny,nz,nnx,nnz,nny,x1,y1,z1,maxai,minai
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dims,1:dims)::ainvjt
real,dimension(1:dimensiona)::cords
real,allocatable,dimension(:)::intbs,basefaceval,basefacgval,permutation,permutationg,xder,yder,zder
real,allocatable,dimension(:,:)::stencil,invmat,wlsqr,lsqm,lscqm,qff,rff,qtff,invrff,vellsqmat
real,allocatable,dimension(:,:,:)::stencils
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



allocate(lscqm(1:idegfree,1:idegfree),qff(1:idegfree,1:idegfree),rff(1:idegfree,1:idegfree),qtff(1:idegfree,1:idegfree),invrff(1:idegfree,1:idegfree))
allocate(vellsqmat(1:idegfree-1,1:idegfree-1))
allocate(lsqm(1:imaxdegfree,1:idegfree-1))
allocate(intbs(1:idegfree),basefaceval(1:idegfree),basefacgval(1:idegfree),permutation(1:idegfree),permutationg(1:idegfree),xder(1:idegfree),yder(1:idegfree),zder(1:idegfree))
allocate(invmat(1:idegfree,1:idegfree))
allocate(wlsqr(1:20,1:numneighbours-1))
allocate(stencil(1:numneighbours-1,1:idegfree))
allocate(stencils(1:20,1:numneighbours-1,1:idegfree))























i=iconsi

        
		idum=0;
                if (ielem_interior(i).eq.1)then
                        do j=1,ielem_ifca(i)
                        if (ielem_ibounds(j,i).gt.0)then
                            if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
                                idum=1
                            end if
                        end if
                        end do
                end if
		
		
		

		
		
		
		
		
		intbs=zero;jxx=1;ixx=i;lxx1=1;number_of_dog=ielem_idegfree(i);kxx=ielem_iorder(i);eltype=ielem_ishape(i)
		icompwrt=0

		intbs=calintbasis(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,eltype,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
        ! store integration-basis vector intbs into all entries 1:number_of_dog for element i

            integ_basis_value(1:number_of_dog,i) = intbs(1:number_of_dog)

                
                
        !the indicator matrix for the large polynomial
		if (iweno.eq.1)then
		call indicatormatrix(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
		end if

		if (ees.eq.5)then
		intbs=zero;jxx=1;ixx=i;lxx1=1
		number_of_dog=idegfree2;kxx=iorder2;eltype=ielem_ishape(i);

		icompwrt=1

		intbs=calintbasis(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,eltype,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
		
		integ_basis_valuec(1:number_of_dog,i)=intbs(1:number_of_dog)

		 !the indicator matrix for the smaller polynomial
		if (iweno.eq.1)then
		call indicatormatrix2(n,i,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
		icompwrt=0
		end if
		
		end if
		llco=ielem_admis(i)

                
		imax=ielem_inumneighbours(i)-1;inum=ielem_inumneighbours(i);ideg=ielem_idegfree(i);inumo=ielem_iorder(i)
         imax2=numneighbours2-1;inum2=numneighbours2;ideg2=idegfree2;inumo2=iorder2
			dist_sten=zero
			!now we start with looping all the admissible stencils
		do ll=1,llco	!admis

			

				if((ees.ne.5).or.(ll.eq.1))then
				imax=ielem_inumneighbours(i)-1;inum=ielem_inumneighbours(i);ideg=ielem_idegfree(i);inumo=ielem_iorder(i)
				icompwrt=0;number_of_dog=ideg

				else
				imax=numneighbours2-1;inum=numneighbours2;ideg=idegfree2;inumo=iorder2;number_of_dog=ideg
				icompwrt=1
				end if




					do k=1,imax !for all neighbours
							ixx=i;kxx=inumo

							if (weight_lsqr.eq.1)then
							wlsqr(ll,k)=1.0d0/((sqrt(((ilox_xxc(ll,k+1)-ilox_xxc(ll,1))**2)+((ilox_yyc(ll,k+1)-ilox_yyc(ll,1))**2))))
							else
							wlsqr(ll,k)=1.0d0
							end if
							x_stencil=(ilox_xxc(ll,k+1)-ilox_xxc(ll,1))**2
							y_stencil=(ilox_yyc(ll,k+1)-ilox_yyc(ll,1))**2


							dist_sten2=sqrt(x_stencil+y_stencil)

							dist_sten=max(dist_sten,dist_sten2)


							if (weight_lsqr.eq.1)then
							wlsqr(ll,k)=1.0d0/sqrt(x_stencil+y_stencil)
							else
							wlsqr(ll,k)=1.0d0
							end if

! 							ielem_stencil_dist(i)=dist_sten/(ilox_volume(1,1)**(1/2))

								if (fastest.eq.1)then
									x1 = ilox_xxc(ll,k+1)-ilox_xxc(ll,1)
									y1 = ilox_yyc(ll,k+1)-ilox_yyc(ll,1)

										if((ees.ne.5).or.(ll.eq.1))then
										icompwrt=0
											if (itestcase.eq.4)then
											rec_stencils(ll,k,1:ielem_idegfree(i),rec_wall(i))=wlsqr(ll,k)*basis_rec2d(n,x1,y1,ielem_iorder(i),ixx,ielem_idegfree(i),icompwrt)
											rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
											stencils(ll,k,1:ideg)=rec_stencils(ll,k,1:ielem_idegfree(i),rec_wall(i))
											else
											stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec2d(n,x1,y1,ielem_iorder(i),ixx,ielem_idegfree(i),icompwrt)
											end if
										else
										icompwrt=1
										if (itestcase.eq.4)then
! 										rec_stencilsc(ll,k,1:ideg,rec_wall(i))=wlsqr(ll,k)*basis_rec2d(n,x1,y1,inumo,ixx,ideg,icompwrt)
! 										rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec2d(n,x1,y1,inumo,ixx,ideg,icompwrt)
										else
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*basis_rec2d(n,x1,y1,inumo,ixx,ideg,icompwrt)
										end if

										end if

								else
									ixx=i;jxx=k+1;lxx1=ll
									eltype=ilox_ishape(ll,k+1)

									if (greengo.eq.0)then

										if (idum.eq.1)then
											if((ees.ne.5).or.(ll.eq.1))then
											icompwrt=0
											if (itestcase.eq.4)then
											rec_stencils(ll,k,1:ielem_idegfree(i),rec_wall(i))=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
											rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
											stencils(ll,k,1:ideg)=rec_stencils(ll,k,1:ielem_idegfree(i),rec_wall(i))
											else
											stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)



											end if

											else
											icompwrt=1
											if (itestcase.eq.4)then
! 											rec_stencilsc(ll,k,1:ideg,rec_wall(i))=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
! 											rec_weightl(ll,k,rec_wall(i))=wlsqr(ll,k)
											stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
											else
											stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
											end if


											icompwrt=0
											end if
										else
											if((ees.ne.5).or.(ll.eq.1))then
											icompwrt=0
											else
											icompwrt=1
											end if
										stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
										icompwrt=0
										end if

									else
									if((ees.ne.5).or.(ll.eq.1))then
									icompwrt=0
									else
									icompwrt=1
									end if
									stencils(ll,k,1:ideg)=wlsqr(ll,k)*compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
									icompwrt=0
									end if
								end if		!fastest
							end do	!imaxedegfree





     
      if ((greengo.eq.0))then
		if (idum.eq.1)then

		      invmat=zero
		      lscqm(:,:)=zero
		      do iq=1,ideg;do jq=1,ideg;do lq=1,imax
		      if((ees.ne.5).or.(ll.eq.1))then
		      lscqm(jq,iq)=lscqm(jq,iq)&
		      +((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
		      !lscqm(jq,iq)=lscqm(jq,iq)&
		      !+((rec_stencils(ll,lq,jq,rec_wall(i))*rec_stencils(ll,lq,iq,rec_wall(i))))
		      else
		      lscqm(jq,iq)=lscqm(jq,iq)&
		      +((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
		      !lscqm(jq,iq)=lscqm(jq,iq)&
		      !+((rec_stencilsc(ll,lq,jq,rec_wall(i))*rec_stencilsc(ll,lq,iq,rec_wall(i))))
		      end if
		      end do;end do;end do  
		else

		      
		      invmat=zero
		      lscqm(:,:)=zero
		      do iq=1,ideg;do jq=1,ideg;do lq=1,imax
		      lscqm(jq,iq)=lscqm(jq,iq)&
		      +((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
		      end do;end do;end do 
		end if


      else

		
		invmat=zero
		lscqm(:,:)=zero
		do iq=1,ideg;do jq=1,ideg;do lq=1,imax
		lscqm(jq,iq)=lscqm(jq,iq)&
		+((stencils(ll,lq,jq)*stencils(ll,lq,iq)))
		end do;end do;end do 

	end if 








			qff(:,:)=zero; rff(:,:)=zero; qtff(:,:)=zero; rff(:,:)=zero;  invrff(:,:)=zero
			call qrdecomposition(lscqm,qff,rff,ideg)
			call transposematrix(qff,qtff,ideg)
			ivgt=ideg+1
			call invert(rff,invrff,ivgt)
			  do mm_i=1,ideg
			    do mm_j=1,ideg
			      invmat(mm_i,mm_j)=zero
			      do mm_k=1,ideg
			        invmat(mm_i,mm_j)=invmat(mm_i,mm_j)+invrff(mm_i,mm_k)*qtff(mm_k,mm_j)
			      end do
			    end do
			  end do




! 
! 
! 

if (greengo.eq.0)then
    if (idum.eq.1)then
            if((ees.ne.5).or.(ll.eq.1))then
            
            stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)!rec_stencils(ll,1:imax,1:ideg,rec_wall(i))
            else
           
            stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)!rec_stencilsc(ll,1:imax,1:ideg,rec_wall(i))

            end if
    else
stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)
    end if
else

stencil(1:imax,1:ideg)=stencils(ll,1:imax,1:ideg)

end if

if((ees.ne.5).or.(ll.eq.1))then
! call gemm(                                               &
!                invmat,                                               &
!                stencil,                                              &
!                rec_invmat_stencilt(:,:,ll,i),             &
!                'n',                                                  & ! transposition flag for invmat
!                't'                                                   & ! transposition flag for stencil
!             )
            
            
! call dgemm ('n','t',ideg,imax,ideg,alpha,invmat(1:ideg,1:ideg),ideg,&
! stencil(1:imax,1:ideg),imax,beta,rec_invmat_stencilt(1:ideg,1:imax,ll,i),ideg)

            
              do mm_i=1,ideg
                do mm_j=1,imax
                  rec_invmat_stencilt(mm_i,mm_j,ll,i)=zero
                  do mm_k=1,ideg
                    rec_invmat_stencilt(mm_i,mm_j,ll,i)=rec_invmat_stencilt(mm_i,mm_j,ll,i)+invmat(mm_i,mm_k)*stencil(mm_j,mm_k)
                  end do
                end do
              end do





           do iq=1,imax
			rec_invmat_stencilt(:,iq,ll,i)=rec_invmat_stencilt(:,iq,ll,i)&
			*ilox_volume(ll,iq+1)*wlsqr(ll,iq)
			end do


			if (rec_invmat_stencilt(1,1,ll,i).ne.rec_invmat_stencilt(1,1,ll,i))then	!prevents non-invertible nan stencils from being deployed

					ielem_full(i)=0
					end if



else
! call gemm(                                               &
!                invmat(1:ideg,1:ideg),                                               &
!                stencil(1:imax,1:ideg),                                              &
!                rec_invmat_stenciltc(1:ideg,1:imax,ll,i),             &
!                'n',                                                  & ! transposition flag for invmat
!                't'                                                   & ! transposition flag for stencil
!             )
            
! call dgemm ('n','t',ideg,imax,ideg,alpha,invmat(1:ideg,1:ideg),ideg,&
! stencil(1:imax,1:ideg),imax,beta,rec_invmat_stenciltc(1:ideg,1:imax,ll,i),ideg)

			  do mm_i=1,ideg
			    do mm_j=1,imax
			      rec_invmat_stenciltc(mm_i,mm_j,ll,i)=zero
			      do mm_k=1,ideg
			        rec_invmat_stenciltc(mm_i,mm_j,ll,i)=rec_invmat_stenciltc(mm_i,mm_j,ll,i)+invmat(mm_i,mm_k)*stencil(mm_j,mm_k)
			      end do
			    end do
			  end do



		do iq=1,imax
			rec_invmat_stenciltc(:,iq,ll,i)=rec_invmat_stenciltc(:,iq,ll,i)&
			*ilox_volume(ll,iq+1)*wlsqr(ll,iq)
			end do


			if (rec_invmat_stenciltc(1,1,ll,i).ne.rec_invmat_stenciltc(1,1,ll,i))then	!prevents non-invertible nan stencils from being deployed

					ielem_full(i)=0
					end if

end if





if (ielem_ggs(i).ne.1)then	!ggs
    
if (itestcase.eq.4)then	!test

	

	 

if (ll.eq.1)then		!stencils
! 

  if (idum.eq.1)then		!for wall only
    
	 
	 

    	
    do ihgt=1,2; do ihgj=1,2
    ainvjt(ihgt,ihgj)=rec_invccjac(ihgj,ihgt,i)
   end do; end do

	 idum2=0
	 
	 
				basefaceval=zero
				basefacgval=zero
				permutation=zero
				permutationg=zero
			
	 
	 
	 
	 do j=1,ielem_ifca(i)		!for all faces
	    if (ielem_ibounds(j,i).gt.0)then		!for bounded only
		      if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then		!!for bounded only 2                         
					facex=j;iconsidered=i
								  call coordinates_face_inner2d(n,iconsidered,facex,vext,nodes_list)
								  n_node=2
								    cords(1:2)=zero
								    call cordinates2(n,nodes_list,n_node,cords(1:2))
							    
								    ay=cords(2)
								    ax=cords(1)

					    vext(1,1)=ax;vext(1,2)=ay;
					      do mm_j=1,2
					        mm_old(mm_j)=vext(1,mm_j)-rec_vext_ref(mm_j,i)
					      end do
					      do mm_i=1,2
					        vext(1,mm_i)=zero
					        do mm_j=1,2
					          vext(1,mm_i)=vext(1,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
					        end do
					      end do
					  
					    ax=vext(1,1);ay=vext(1,2)
				



				
				
				
				
				angle1=ielem_faceanglex(j,i)
				angle2=ielem_faceangley(j,i)
				
				nx=angle1
				ny=angle2
				
				nnx=(nx*ainvjt(1,1))+(ny*ainvjt(2,1))
				nny=(nx*ainvjt(1,2))+(ny*ainvjt(2,2))
				
				if (poly.eq.4)then


				do iq=1, ideg

				xder(iq)=tl2dx(ax,ay,iq,i);  yder(iq)=tl2dy(ax,ay,iq,i);


				end do


				else



				do iq=1, ideg
				
				xder(iq)=df2dx(ax,ay,iq,i);  yder(iq)=df2dy(ax,ay,iq,i);
				
				
				end do

				end if



				iconsidered=i
				icompwrt=0
				basefaceval(1:ielem_idegfree(i))=basis_rec2d(n,ax,ay,ielem_iorder(i),iconsidered,ielem_idegfree(i),icompwrt)

				!if (thermal.eq.0)then
				icompwrt=0

				basefacgval(1:ielem_idegfree(i))=((nnx*xder(1:ielem_idegfree(i)))+(nny*yder(1:ielem_idegfree(i))))
				!else
				!icompwrt=0
				!basefacgval(1:ielem_idegfree(i))=basis_rec2d(n,ax,ay,ielem_iorder(i),iconsidered,ielem_idegfree(i),icompwrt)

				!end if


				do iq=1,ideg
				rec_wallcoeff(iq,rec_wall(i))=basefaceval(iq)
				rec_wallcoefg(iq,rec_wall(i))=basefacgval(iq)
				permutation(iq)=iq
				permutationg(iq)=iq
				end do
		
				
		
		
		
		
		
		end if!for bounded only 2        
		end if!for bounded only
		end do!for all faces
		
		gggg=-tolbig
g0=0
do iq=1,ideg
if (abs(basefacgval(iq)) >gggg)then
gggg=abs(basefacgval(iq))
g0=iq
end if
end do




ssss=-tolbig
k0=0
do iq=1,ideg
if (abs(basefaceval(iq)) >ssss)then
ssss=abs(basefaceval(iq))
k0=iq
end if
end do



! now swap basis functions and thus coefficients
permutation(1)=k0; permutation(k0)=1; ssss=basefaceval(1)
basefaceval(1)=basefaceval(k0); basefaceval(k0)=ssss

permutationg(1)=g0; permutation(g0)=1; gggg=basefacgval(1)
basefacgval(1)=basefacgval(g0); basefacgval(g0)=gggg
! rec_k0(ii)=k0
! rec_g0(ii)=g0
rec_k0(rec_wall(i))=k0
rec_g0(rec_wall(i))=g0
		
		
		
		
lsqm = zero
do lq=1,imax
lcou=0
    do iq=1,ideg
	if (iq.eq.g0) cycle
	    lcou=lcou+1
lsqm(lq,lcou)=rec_stencils(ll,lq,iq,rec_wall(i))&
-rec_stencils(ll,lq,g0,rec_wall(i))*rec_wallcoefg(iq,rec_wall(i))/rec_wallcoefg(g0,rec_wall(i))
    end do
end do
rec_tempsq(1:imax,1:ideg-1,rec_wall(i))=lsqm(1:imax,1:ideg-1)






vellsqmat=zero
do iq=1,ideg-1; do jq=1,ideg-1;	do lcc=1,imax
!now store the least square matrix
vellsqmat(jq,iq)= vellsqmat(jq,iq)+(lsqm(lcc,jq)*lsqm(lcc,iq))
end do;	end do;	end do
lscqm=zero
lscqm(1:ideg-1,1:ideg-1)=vellsqmat(1:ideg-1,1:ideg-1)

qff=zero; rff=zero; qtff=zero; invrff=zero
	  ivgt=ideg
    call qrdecomposition(lscqm,qff,rff,ivgt-1)
    call transposematrix(qff,qtff,ivgt-1)
    
    call invert(rff,invrff,ivgt)
!   final inverted r^(-1)*q^(-1)
  do mm_i=1,ideg-1
    do mm_j=1,ideg-1
      rec_tempsqmat(mm_i,mm_j,rec_wall(i))=zero
      do mm_k=1,ideg-1
        rec_tempsqmat(mm_i,mm_j,rec_wall(i))=rec_tempsqmat(mm_i,mm_j,rec_wall(i))+invrff(mm_i,mm_k)*qtff(mm_k,mm_j)
      end do
    end do
  end do

! definy the matrix defining the least-square reconstruction in this case for velocity
lsqm = zero
do lq=1,imax
lcou=0
    do iq=1,ideg
	if (iq.eq.k0) cycle
	    lcou=lcou+1
	     lsqm(lq,lcou)=rec_stencils(ll,lq,iq,rec_wall(i))&
 -rec_stencils(ll,lq,k0,rec_wall(i))*rec_wallcoeff(iq,rec_wall(i))/rec_wallcoeff(k0,rec_wall(i))
    end do
end do
rec_vellsq(1:imax,1:ideg-1,rec_wall(i))=lsqm(1:imax,1:ideg-1)
vellsqmat=zero




do iq=1,ideg-1; do jq=1,ideg-1;	do lcc=1,imax
!now store the least square matrix
vellsqmat(jq,iq)= vellsqmat(jq,iq)+(lsqm(lcc,jq)*lsqm(lcc,iq))
end do;	end do;	end do

lscqm=zero
lscqm(1:ideg-1,1:ideg-1)=vellsqmat(1:ideg-1,1:ideg-1)


qff=zero; rff=zero; qtff=zero; invrff=zero
	  ivgt=ideg
    call qrdecomposition(lscqm,qff,rff,ivgt-1)
    call transposematrix(qff,qtff,ivgt-1)
    
    call invert(rff,invrff,ivgt)
!   final inverted r^(-1)*q^(-1)
  do mm_i=1,ideg-1
    do mm_j=1,ideg-1
      rec_velinvlsqmat(mm_i,mm_j,rec_wall(i))=zero
      do mm_k=1,ideg-1
        rec_velinvlsqmat(mm_i,mm_j,rec_wall(i))=rec_velinvlsqmat(mm_i,mm_j,rec_wall(i))+invrff(mm_i,mm_k)*qtff(mm_k,mm_j)
      end do
    end do
  end do
			
		
		
		
		
		
		
		end if!for wall only
		end if!stencils
		end if!for test
		end if!ggs
		end do
		
		
deallocate(lsqm,qff,rff,qtff,invrff)
		deallocate(vellsqmat)
		deallocate(lscqm)
		deallocate(intbs,basefaceval,basefacgval,permutation,permutationg,xder,yder,zder)
		deallocate(wlsqr)
		deallocate(stencil)
		deallocate(invmat)
		deallocate(stencils)


end subroutine prestore_reconstruction2


subroutine indicatormatrix(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine computes the indicator matrices for weno reconstructions
implicit none
integer,intent(in)::n,iconsi
integer::i,j,k,l,m,jx,jx2,imax,inum,ideg,inumo,eltype,elem_dec,inump,iconsidered
real::voltemp
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
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
real,allocatable,dimension(:,:)::weff
allocate(weff(1:idegfree,1:idegfree))
i=iconsi


weff=zero

	imax=ielem_inumneighbours(i)-1
	inum=ielem_inumneighbours(i)
	ideg=ielem_idegfree(i)
	 inumo=ielem_iorder(i)
	 
	iconsidered=i
	 

	vext=zero
    nodes_list=zero
    eltype=ielem_ishape(i)
    elem_dec=ielem_vdec(i)
    elem_listd=zero
          rec_indicator(1:ideg,1:ideg,i)=zero
      jx=ielem_nonodes(i)
     
	  if (dimensiona.eq.3)then
	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,1)=ilon_x(1,1,k)
	    nodes_list(k,2)=ilon_y(1,1,k)
	    nodes_list(k,3)=ilon_z(1,1,k)
	    vext(k,:)=nodes_list(k,:)
	  end do
	  call decompose3(n,eltype,nodes_list,elem_listd)
	  
	  else

	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,1)=ilon_x(1,1,k)
	    nodes_list(k,2)=ilon_y(1,1,k)
	    vext(k,1:2)=nodes_list(k,1:2)
	  end do
	  call decompose2(n,eltype,nodes_list,elem_listd)
	  	  end if

	 select case(ielem_ishape(i))

      case(1)
      call quadraturehexa(n,igqrules,vext,qpoints,wequa3d)
      voltemp=hexavolume(n,vext,qpoints,wequa3d)
      inump=qp_hexa
      case(2)
      call quadraturetetra(n,igqrules,vext,qpoints,wequa3d)
      voltemp=tetravolume(n,vext)
      inump=qp_tetra
      case(3)
      call quadraturepyra(n,igqrules,vext,qpoints,wequa3d)
      voltemp=pyravolume(n,vext,qpoints,wequa3d)
      inump=qp_pyra
      case(4)
      call quadratureprism(n,igqrules,vext,qpoints,wequa3d)
      voltemp=prismvolume(n,vext,qpoints,wequa3d)
      inump=qp_prism

      case(5)
      call quadraturequad(n,igqrules,vext,qpoints,wequa3d)
      voltemp=quadvolume(n,vext,qpoints,wequa3d)
      inump=qp_quad


      case(6)
      call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)
      voltemp=trianglevolume(n,vext)
      inump=qp_triangle

      end select      
      
	if (dimensiona.eq.3)then
	if (ielem_mode(i).eq.1)then
	  do k=1,elem_dec
	      vext(1:4,1:3)=elem_listd(k,1:4,1:3)
	      call quadraturetetra(n,igqrules,vext,qpoints,wequa3d)
	      voltemp=tetravolume(n,vext)
	      inump=qp_tetra
	      call wenotet(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
	      rec_indicator(1:ideg,1:ideg,i)=rec_indicator(1:ideg,1:ideg,i)+weff(1:ideg,1:ideg)
	
	  end do
	else
	  

	    call wenotet(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
	      rec_indicator(1:ideg,1:ideg,i)=rec_indicator(1:ideg,1:ideg,i)+weff(1:ideg,1:ideg)






	end if
	else

	    do k=1,elem_dec
            vext(1:3,1:2)=elem_listd(k,1:3,1:2)
            call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)
            voltemp=trianglevolume(n,vext)
            inump=qp_triangle
           
	    call wenotet2d(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
             
	      rec_indicator(1:ideg,1:ideg,i)=rec_indicator(1:ideg,1:ideg,i)+weff(1:ideg,1:ideg)
            end do


	end if

deallocate(weff)


end subroutine indicatormatrix


subroutine indicatormatrix2(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine computes the indicator matrices for cweno reconstructions of the lower-order polynomials
implicit none
integer,intent(in)::n,iconsi
integer::i,j,k,l,m,jx,jx2,imax,inum,ideg,inumo,eltype,elem_dec,inump,iconsidered
real::voltemp
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
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
real,allocatable,dimension(:,:)::weff
allocate(weff(1:idegfree2,1:idegfree2))
i=iconsi

	imax=numneighbours2-1
	inum=numneighbours2
	ideg=idegfree2
	 inumo=iorder2
	 weff=zero
	
	 iconsidered=i

	vext=zero
    nodes_list=zero
    eltype=ielem_ishape(i)
    elem_dec=ielem_vdec(i)
    elem_listd=zero
          rec_indicatorc(1:ideg,1:ideg,i)=zero
      jx=ielem_nonodes(i)
     
	  if (dimensiona.eq.3)then
	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,1)=ilon_x(1,1,k)
	    nodes_list(k,2)=ilon_y(1,1,k)
	    nodes_list(k,3)=ilon_z(1,1,k)
	    vext(k,:)=nodes_list(k,:)
	  end do
	  call decompose3(n,eltype,nodes_list,elem_listd)
	  
	  else

	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,1)=ilon_x(1,1,k)
	    nodes_list(k,2)=ilon_y(1,1,k)
	    vext(k,:)=nodes_list(k,:)
	  end do
	  call decompose2(n,eltype,nodes_list,elem_listd)
	  	  end if

	 select case(ielem_ishape(i))

      case(1)
      call quadraturehexa(n,igqrules,vext,qpoints,wequa3d)
      voltemp=hexavolume(n,vext,qpoints,wequa3d)
      inump=qp_hexa
      case(2)
      call quadraturetetra(n,igqrules,vext,qpoints,wequa3d)
      voltemp=tetravolume(n,vext)
      inump=qp_tetra
      case(3)
      call quadraturepyra(n,igqrules,vext,qpoints,wequa3d)
      voltemp=pyravolume(n,vext,qpoints,wequa3d)
      inump=qp_pyra
      case(4)
      call quadratureprism(n,igqrules,vext,qpoints,wequa3d)
      voltemp=prismvolume(n,vext,qpoints,wequa3d)
      inump=qp_prism

      case(5)
      call quadraturequad(n,igqrules,vext,qpoints,wequa3d)
      voltemp=quadvolume(n,vext,qpoints,wequa3d)
      inump=qp_quad


      case(6)
      call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)
      voltemp=trianglevolume(n,vext)
      inump=qp_triangle

      end select      
      
	if (dimensiona.eq.3)then
	if (ielem_mode(i).eq.1)then
	  do k=1,elem_dec
	      vext(1:4,1:3)=elem_listd(k,1:4,1:3)
	      call quadraturetetra(n,igqrules,vext,qpoints,wequa3d)
	      voltemp=tetravolume(n,vext)
	      inump=qp_tetra
	      call wenotet(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
	      rec_indicatorc(1:ideg,1:ideg,i)=rec_indicatorc(1:ideg,1:ideg,i)+weff(1:ideg,1:ideg)
	
	  end do
	else
	  

	    call wenotet(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
	      rec_indicatorc(1:ideg,1:ideg,i)=rec_indicatorc(1:ideg,1:ideg,i)+weff(1:ideg,1:ideg)






	end if
	else
            
            do k=1,elem_dec
            vext(1:3,1:2)=elem_listd(k,1:3,1:2)
            call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)
            voltemp=trianglevolume(n,vext)
            inump=qp_triangle

	    call wenotet2d(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

	      rec_indicatorc(1:ideg,1:ideg,i)=rec_indicatorc(1:ideg,1:ideg,i)+weff(1:ideg,1:ideg)
            end do
            
	end if

	
	
deallocate(weff)


end subroutine indicatormatrix2

	



subroutine wenotet(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n,inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real,allocatable,dimension(:,:),intent(inout)::weff
weff=zero


		if (poly.eq.1)then
           select case(inumo)
		 case(1)
		call ph1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(2)
		call ph2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(3)
		call ph3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(4)
		call ph4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(5)
		call ph5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(6)
		call ph6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
		
	end select
	end if
	if (poly.eq.2)then
           select case(inumo)
		 case(1)
		call pl1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(2)
		call pl2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(3)
		call pl3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(4)
		call pl4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(5)
		call pl5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(6)
		call pl6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
		
	end select
	end if
	
	if (poly.eq.4)then
           select case(inumo)
		 case(1)
		call tl3d1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(2)
		call tl3d2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(3)
		call tl3d3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(4)
		call tl3d4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(5)
		call tl3d5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(6)
		call tl3d6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
		
	end select
	end if
	
	
	

end subroutine wenotet

subroutine ph1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
			integ =zero
                scalerx=1.0d0
                
                
       			 do k=1,inump
	
	    ph=dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

	    integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine ph1
subroutine ph2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
  do k=1,inump
	ph=dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)
	
	integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
  enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine ph2
subroutine ph3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
  do k=1,inump
	  ph=dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
	dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)
	

	integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
  enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine ph3
subroutine ph4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
do k=1,inump
    ph=dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
    dfz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

    integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine ph4
subroutine ph5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
do i=1,ideg
do j=1,ideg
scalerx=1.0d0
                
	integ =zero
	  do k=1,inump
ph=dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
  dfz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

		integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
	  enddo
	weff(i,j)=weff(i,j)+integ
enddo
enddo

end subroutine ph5
subroutine ph6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real::ph,integ,scalerx
integer::i,j,k
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
			ph=dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		    dfx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		    dfx6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx5y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx5y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx4y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx4yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2y2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2y2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2yz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2yz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfx2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfx2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfy2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfy2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfxz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfxz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfyz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfyz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dfz6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dfz6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine ph6

subroutine pl1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
	
				ph=dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				   dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				   dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine pl1
subroutine pl2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)
				
				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine pl2
subroutine pl3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				 ph=dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)
				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine pl3
subroutine pl4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dly2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				dlz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine pl4
subroutine pl5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
		ph=dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		  dlx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine pl5
subroutine pl6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
			     ph=dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		    dlx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		    dlx6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx5y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx5y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2y2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2y2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2yz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2yz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlx2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlx2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dly2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dly2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlxz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlxz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlyz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlyz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  dlz6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*dlz6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine pl6


subroutine tl3d1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
	
				ph=tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				   tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				   tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl3d1
subroutine tl3d2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)
				
				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl3d2
subroutine tl3d3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				 ph=tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)
				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl3d3
subroutine tl3d4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
				tl3dz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl3d4
subroutine tl3d5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
		ph=tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		  tl3dx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl3d5
subroutine tl3d6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
			     ph=tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		    tl3dx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
		    tl3dx6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx5y(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx5y(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4y2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4yz(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3y2z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y3z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy4z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy5z(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3yz2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2y2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2y2z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy3z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy4z2(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2yz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2yz3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxy2z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy3z3(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dx2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dx2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxyz4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dy2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dy2z4(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dxz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dxz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dyz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dyz5(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered) + &
                  tl3dz6(qpoints(1,k),qpoints(2,k),qpoints(3,k),i,iconsidered)*tl3dz6(qpoints(1,k),qpoints(2,k),qpoints(3,k),j,iconsidered)

				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl3d6





subroutine wenotet2d(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n,inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real,allocatable,dimension(:,:),intent(inout)::weff
weff=zero

 
		if (poly.eq.1)then
           select case(inumo)
		 case(1)
		call p2dh1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(2)
		call p2dh2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(3)
		call p2dh3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(4)
		call p2dh4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(5)
		call p2dh5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(6)
		call p2dh6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
		
	end select
	
	end if
	
	
	
	if (poly.eq.4)then
           select case(inumo)
		 case(1)
		call tl2dh1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(2)
		call tl2dh2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(3)
		call tl2dh3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(4)
		call tl2dh4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(5)
		call tl2dh5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)

		case(6)
		call tl2dh6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
		
	end select
	
	end if
	
	
	

end subroutine wenotet2d
! 
subroutine p2dh1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
	
				ph=df2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				   df2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy(qpoints(1,k),qpoints(2,k),j,iconsidered)
				   

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
!           		
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine p2dh1
subroutine p2dh2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=df2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered)
				
				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine p2dh2
subroutine p2dh3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
        	
        	
        	
        	
			integ =zero
       			 do k=1,inump
				 ph=df2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered)
				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine p2dh3
subroutine p2dh4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=df2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx3y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dx2y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dxy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				df2dy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy4(qpoints(1,k),qpoints(2,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine p2dh4
subroutine p2dh5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
		ph=df2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
		  df2dx5(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx4y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx4y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy5(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy5(qpoints(1,k),qpoints(2,k),j,iconsidered)
                 

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine p2dh5
subroutine p2dh6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
			ph=df2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
		    df2dx5(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx4y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx4y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy5(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
		    df2dx6(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx6(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx5y(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx5y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx4y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx4y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx3y3(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx3y3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dx2y4(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dx2y4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dxy5(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dxy5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  df2dy6(qpoints(1,k),qpoints(2,k),i,iconsidered)*df2dy6(qpoints(1,k),qpoints(2,k),j,iconsidered)

				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine p2dh6

subroutine tl2dh1(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
	
				ph=tl2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				   tl2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy(qpoints(1,k),qpoints(2,k),j,iconsidered)
				   

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
!           		
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl2dh1
subroutine tl2dh2(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=tl2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered)
				
				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl2dh2
subroutine tl2dh3(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
        	
        	
        	
        	
			integ =zero
       			 do k=1,inump
				 ph=tl2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered)
				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl2dh3
subroutine tl2dh4(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
				ph=tl2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx3y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dx2y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dxy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
				tl2dy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy4(qpoints(1,k),qpoints(2,k),j,iconsidered)

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl2dh4
subroutine tl2dh5(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
		ph=tl2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
		  tl2dx5(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx4y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx4y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy5(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy5(qpoints(1,k),qpoints(2,k),j,iconsidered)
                 

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl2dh5
subroutine tl2dh6(n,weff,inumo,inump,voltemp,qpoints,wequa3d,ideg,iconsidered)
!> @brief
!> this subroutine computes derivative operator for the smoothness indicator
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::weff
integer,intent(in)::inumo,inump,ideg,iconsidered
real,intent(in)::voltemp
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
integer::i,j,k
real::ph,integ,scalerx
weff=zero
		
 	do i=1,ideg
        	do j=1,ideg
        	scalerx=1.0d0
                
			integ =zero
       			 do k=1,inump
			ph=tl2dx(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
		    tl2dx5(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx4y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx4y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy5(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
		    tl2dx6(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx6(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx5y(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx5y(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx4y2(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx4y2(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx3y3(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx3y3(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dx2y4(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dx2y4(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dxy5(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dxy5(qpoints(1,k),qpoints(2,k),j,iconsidered) + &
                  tl2dy6(qpoints(1,k),qpoints(2,k),i,iconsidered)*tl2dy6(qpoints(1,k),qpoints(2,k),j,iconsidered)

				

				integ=integ+(ph*wequa3d(k)*voltemp)*(scalerx**(i+j-1))
          		 enddo
         		weff(i,j)=weff(i,j)+integ
        	enddo
        enddo

end subroutine tl2dh6



function calintbasis(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,eltype,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine computes basis functions for each element
implicit none
integer,intent(in)::n,number_of_dog,icompwrt,eltype
integer,intent(in):: ixx,jxx,kxx,lxx1
integer::k
real,dimension(1:number_of_dog)::calintbasis
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


! 	kxx=ielem_iorder(ixx)
	calintbasis = zero
     	
	calintbasis(1:number_of_dog)=compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
	
   end function


function compbasel(n,eltype,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine computes basis functions for each element
implicit none
integer,intent(in)::n,lxx1,jxx,ixx,icompwrt,kxx
integer,intent(in)::eltype,number_of_dog
integer::jx,k,elem_dec
real,dimension(1:number_of_dog)::compbasel,s1
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
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
s1=zero

 select case(eltype)
      case(1)
      jx=8;   elem_dec=6

		 if ((jxx.eq.1))then
		  
		      do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			  nodes_list(k,3)=ilon_z(lxx1,jxx,k)
			vext(k,:)=nodes_list(k,:)
		      end do
		  if (ielem_mode(ixx).eq.0)then
		  
		  compbasel=compbashex(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)

		    else
		  call decompose3(n,eltype,nodes_list,elem_listd)
		  do k=1,elem_dec
		    vext(1:4,1)=elem_listd(k,1:4,1)
		    vext(1:4,2)=elem_listd(k,1:4,2)
		    vext(1:4,3)=elem_listd(k,1:4,3)
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		  end do
		    compbasel=s1
		    end if
		  
		  else

		   do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			  nodes_list(k,3)=ilon_z(lxx1,jxx,k)
			vext(k,:)=nodes_list(k,:)
			
		      end do
		      
		      
		  call decompose3(n,eltype,nodes_list,elem_listd)
		  do k=1,elem_dec
		    vext(1:4,1)=elem_listd(k,1:4,1)
		    vext(1:4,2)=elem_listd(k,1:4,2)
		    vext(1:4,3)=elem_listd(k,1:4,3)
		    
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		    
		  end do
		    compbasel=s1
		 end if

      case(2)	

		jx=4;   elem_dec=1
		do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			  nodes_list(k,3)=ilon_z(lxx1,jxx,k)
			vext(k,:)=nodes_list(k,:)
		      end do
		  call decompose3(n,eltype,nodes_list,elem_listd)
		  do k=1,elem_dec
		    vext(1:4,1)=elem_listd(k,1:4,1)
		    vext(1:4,2)=elem_listd(k,1:4,2)
		    vext(1:4,3)=elem_listd(k,1:4,3)
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		  end do
		    compbasel=s1
		
      case(3)

		jx=5;   elem_dec=2
		do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			  nodes_list(k,3)=ilon_z(lxx1,jxx,k)
			vext(k,:)=nodes_list(k,:)
		      end do
		  call decompose3(n,eltype,nodes_list,elem_listd)
		  do k=1,elem_dec
		    vext(1:4,1)=elem_listd(k,1:4,1)
		    vext(1:4,2)=elem_listd(k,1:4,2)
		    vext(1:4,3)=elem_listd(k,1:4,3)
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		  end do
		    compbasel=s1
	case(4)
		jx=6;   elem_dec=3

	      if ((jxx.eq.1))then
		  
		      do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			  nodes_list(k,3)=ilon_z(lxx1,jxx,k)
			vext(k,:)=nodes_list(k,:)
		      end do
		      call decompose3(n,eltype,nodes_list,elem_listd)
		  if (ielem_mode(ixx).eq.0)then
		  
		  compbasel=compbaspr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)

		    else
		  
		  do k=1,elem_dec
		   vext(1:4,1)=elem_listd(k,1:4,1)
		    vext(1:4,2)=elem_listd(k,1:4,2)
		    vext(1:4,3)=elem_listd(k,1:4,3)
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
! 		   
		  end do
		    compbasel=s1
		    end if
		  
		  else
			    
		   do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			  nodes_list(k,3)=ilon_z(lxx1,jxx,k)
			vext(k,:)=nodes_list(k,:)
			
		      end do
		  call decompose3(n,eltype,nodes_list,elem_listd)
		  
		  do k=1,elem_dec
		    vext(1:4,1)=elem_listd(k,1:4,1)
		    vext(1:4,2)=elem_listd(k,1:4,2)
		    vext(1:4,3)=elem_listd(k,1:4,3)
		   
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		    
		  end do
		    compbasel=s1
		 end if

	      case(5)
		jx=4;   elem_dec=2
		     if ((jxx.eq.1))then
		     
		     
		   do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
		    vext(k,1:2)=nodes_list(k,1:2)
		    
		    end do
		    call decompose2(n,eltype,nodes_list,elem_listd)
		    if (ielem_mode(ixx).eq.0)then
		  compbasel=compbasquad(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		    else
		    
		  do k=1,elem_dec
		    vext(1:3,1)=elem_listd(k,1:3,1)
		    vext(1:3,2)=elem_listd(k,1:3,2)
		    
		    
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastri(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		    
		  end do
		    compbasel=s1
		    
		    
		    end if
		    
		    
		    
		else
		  do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
			 
			vext(k,:)=nodes_list(k,:)
			
		      end do
		      
		      
		  call decompose2(n,eltype,nodes_list,elem_listd)
		  do k=1,elem_dec
		    vext(1:3,1)=elem_listd(k,1:3,1)
		    vext(1:3,2)=elem_listd(k,1:3,2)
		    
		    
		    s1(1:number_of_dog)=s1(1:number_of_dog)+compbastri(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
		    
		  end do
		    compbasel=s1	
		
		end if

	      case(6)

		jx=3;
		
		
		   do k=1,jx
			  nodes_list(k,1)=ilon_x(lxx1,jxx,k)
			  nodes_list(k,2)=ilon_y(lxx1,jxx,k)
		    vext(k,1:2)=nodes_list(k,1:2)
		    
		    end do
		    
		  compbasel=compbastri(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)


     end select

end function

function compbastr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
!> @brief
!> this subroutine computes basis functions for each triangle
implicit none
integer :: kbasis
integer,intent(in)::n,number_of_dog,ixx,jxx,kxx,lxx1,icompwrt
real::vol
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real::x1,y1,z1
integer::lc
real,dimension(1:number_of_dog)::integ,compbastr


qpoints(1:3,1:qp_tetra)=zero
wequa3d(1:qp_tetra)=zero

vol=tetravolume(n,vext)

call quadraturetetra(n,igqrules,vext,qpoints,wequa3d)

integ=zero
do lc=1,qp_tetra
	  x1=qpoints(1,lc);y1=qpoints(2,lc);z1=qpoints(3,lc)
	integ(1:number_of_dog)=integ(1:number_of_dog)+(basis_rec(n,x1,y1,z1,kxx,ixx,number_of_dog,icompwrt)*wequa3d(lc)*vol)
end do
	  compbastr = integ

end function


function compbashex(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
!> @brief
!> this subroutine computes basis functions for hexahedral
implicit none
integer :: kbasis
integer,intent(in)::n,number_of_dog,ixx,jxx,kxx,lxx1,icompwrt
real::vol
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real::x1,y1,z1
integer::lc
real,dimension(number_of_dog)::integ,compbashex

qpoints(1:3,1:qp_hexa)=zero
wequa3d(1:qp_hexa)=zero

call quadraturehexa(n,igqrules,vext,qpoints,wequa3d)


 vol=hexavolume(n,vext,qpoints,wequa3d)

integ=zero
do lc=1,qp_hexa
	x1=qpoints(1,lc);y1=qpoints(2,lc);z1=qpoints(3,lc)
	integ(1:number_of_dog)=integ(1:number_of_dog)+(basis_rec(n,x1,y1,z1,kxx,ixx,number_of_dog,icompwrt)*wequa3d(lc)*vol)
end do
	  compbashex = integ

end function

function compbaspr(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
!> @brief
!> this subroutine computes basis functions for prisms
implicit none
integer :: kbasis
integer,intent(in)::n,number_of_dog,ixx,jxx,kxx,lxx1,icompwrt
real::vol
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real::x1,y1,z1
integer::lc
real,dimension(number_of_dog)::integ,compbaspr
qpoints(1:3,1:qp_prism)=zero
wequa3d(1:qp_prism)=zero

call quadratureprism(n,igqrules,vext,qpoints,wequa3d)


 vol=prismvolume(n,vext,qpoints,wequa3d)

integ=zero
do lc=1,qp_prism
	x1=qpoints(1,lc);y1=qpoints(2,lc);z1=qpoints(3,lc)
	integ(1:number_of_dog)=integ(1:number_of_dog)+(basis_rec(n,x1,y1,z1,kxx,ixx,number_of_dog,icompwrt)*wequa3d(lc)*vol)
end do
	  compbaspr= integ

end function

function compbasquad(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
!> @brief
!> this subroutine computes basis functions for quadrilateral
implicit none
integer :: kbasis
integer,intent(in)::n,number_of_dog,ixx,jxx,kxx,lxx1,icompwrt
real::vol
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real::x1,y1,z1
integer::lc
real,dimension(number_of_dog)::integ,compbasquad

qpoints(1:2,1:qp_quad)=zero
wequa3d(1:qp_quad)=zero

call quadraturequad(n,igqrules,vext,qpoints,wequa3d)


 vol=quadvolume(n,vext,qpoints,wequa3d)

integ=zero
do lc=1,qp_quad
	x1=qpoints(1,lc);y1=qpoints(2,lc)
	integ(1:number_of_dog)=integ(1:number_of_dog)+(basis_rec2d(n,x1,y1,kxx,ixx,number_of_dog,icompwrt)*wequa3d(lc)*vol)
end do
	  compbasquad= integ

end function

function compbastri(n,ixx,jxx,kxx,lxx1,number_of_dog,icompwrt,vext)
!> @brief
!> this subroutine computes basis functions for each triangle
implicit none
integer :: kbasis
integer,intent(in)::n,number_of_dog,ixx,jxx,kxx,lxx1,icompwrt
real::vol
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d
real::x1,y1,z1
integer::lc
real,dimension(number_of_dog)::integ,compbastri

qpoints(1:2,1:qp_quad)=zero
wequa3d(1:qp_quad)=zero

call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)


 vol=trianglevolume(n,vext)

integ=zero
do lc=1,qp_triangle
	x1=qpoints(1,lc);y1=qpoints(2,lc)
	integ(1:number_of_dog)=integ(1:number_of_dog)+(basis_rec2d(n,x1,y1,kxx,ixx,number_of_dog,icompwrt)*wequa3d(lc)*vol)
end do
	  compbastri= integ

end function


subroutine invert(rff,invrff,ivgt)
!> @brief
!> this subroutine inverts a special matrix
  implicit none
  integer,intent(in)::ivgt
  real,allocatable,dimension(:,:),intent(in) ::rff
  real,allocatable,dimension(:,:),intent(inout)::invrff
  integer i,j,k,gt
 
  invrff(1:ivgt-1,1:ivgt-1) = zero
 gt=ivgt-1
  do i=gt,1,-1
    invrff(i,i) = 1./rff(i,i)
	do j=i+1,gt
	  invrff(i,j) = zero
	  do k= 1,j-1
	   invrff(i,j) = invrff(i,j) - rff(k,j)*invrff(i,k)
	  enddo
 	  invrff(i,j) =invrff(i,j) /rff(j,j)
	enddo
  enddo

 end subroutine invert
 


 
 
end module prestore
