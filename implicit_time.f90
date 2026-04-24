module implicit_time
use library
use transform
use flow_operations
use implicit_fluxes
use communications
use declaration
implicit none

 contains



subroutine relaxation(n)
 !> @brief
!> this subroutine solves the linear system for implicit time stepping either through jacobian or lu-sgs in 3d
implicit none
integer,intent(in)::n
integer::i,l,k,ii,sweeps,kmaxe,nvar,igoflux, icaseb,n_node,ibfc,srf,j
real::impres1,impres2,impres3,tempxx
real:: w1,w2,w3,denx
real,dimension(1:nof_variables,1:nof_variables)::lscqm1
real::b1_imp(1:nof_variables),du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::b1t(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::iconsidered, facex, pointx
integer::ngp
integer::b_code,nfx,lfx,rowfx
real::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv,srf_speedrot,srf_speed
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords

sweeps=10
kmaxe=xmpielrank(n)

impdu(:,:)=zero


du1=zero
b1_imp=zero
lscqm1=zero
dur=zero; dul=zero
durr=zero; dulr=zero



call calculate_jacobian(n)
if (rframe.eq.0) then
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(i,1:nof_variables,1:nof_variables)
  do j=1,nof_variables
	impdiag(i,j,j)=1.0d0/lscqm1(j,j)
  end do
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
if (srfg.eq.1)then
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  lscqm1(1:5,1:5)=impdiag(i,1:5,1:5)
    w1=lscqm1(4,3)
    w2=lscqm1(2,4)
    w3=lscqm1(3,2)
    !inverse of the matrix in case of source term
    denx=(lscqm1(2,2)*lscqm1(3,3)*lscqm1(4,4)+lscqm1(2,2)*w1**2+lscqm1(3,3)*w2**2+lscqm1(4,4)*w3**2)
    impdiag(i,1,1)=1.0d0/lscqm1(1,1)
    impdiag(i,2,2)=(lscqm1(3,3)*lscqm1(4,4)+w1**2)/denx            
    impdiag(i,2,3)=(lscqm1(4,4)*w3+w1*w2)/denx
    impdiag(i,2,4)=(-lscqm1(3,3)*w2+w1*w3)/denx            
    impdiag(i,3,2)=(-lscqm1(4,4)*w3+w1*w2)/denx
    impdiag(i,3,3)=(lscqm1(2,2)*lscqm1(4,4)+w2**2)/denx
    impdiag(i,3,4)=(lscqm1(2,2)*w1+w2*w3)/denx
    impdiag(i,4,2)=(lscqm1(3,3)*w2+w1*w3)/denx            
    impdiag(i,4,3)=(-lscqm1(2,2)*w1+w2*w3)/denx
    impdiag(i,4,4)=(lscqm1(2,2)*lscqm1(3,3)+w3**2)/denx
    impdiag(i,5,5)=1.0d0/lscqm1(5,5)
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
if(mrf.eq.1)then
 srf=0
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
	srf=rec_mrf(i)
	if (rec_mrf(i).eq.0)then
        lscqm1(1:5,1:5)=impdiag(i,1:5,1:5)
        impdiag(i,1,1)=1.0d0/lscqm1(1,1)
        impdiag(i,2,2)=1.0d0/lscqm1(2,2)
        impdiag(i,3,3)=1.0d0/lscqm1(3,3)
        impdiag(i,4,4)=1.0d0/lscqm1(4,4)
        impdiag(i,5,5)=1.0d0/lscqm1(5,5)
    else
        lscqm1(1:5,1:5)=impdiag(i,1:5,1:5)
        w1=lscqm1(4,3)
        w2=lscqm1(2,4)
        w3=lscqm1(3,2)
        !inverse of the matrix in case of source term
        denx=(lscqm1(2,2)*lscqm1(3,3)*lscqm1(4,4)+lscqm1(2,2)*w1**2+lscqm1(3,3)*w2**2+lscqm1(4,4)*w3**2)
        impdiag(i,1,1)=1.0d0/lscqm1(1,1)
        impdiag(i,2,2)=(lscqm1(3,3)*lscqm1(4,4)+w1**2)/denx            
        impdiag(i,2,3)=(lscqm1(4,4)*w3+w1*w2)/denx
        impdiag(i,2,4)=(-lscqm1(3,3)*w2+w1*w3)/denx            
        impdiag(i,3,2)=(-lscqm1(4,4)*w3+w1*w2)/denx
        impdiag(i,3,3)=(lscqm1(2,2)*lscqm1(4,4)+w2**2)/denx
        impdiag(i,3,4)=(lscqm1(2,2)*w1+w2*w3)/denx
        impdiag(i,4,2)=(lscqm1(3,3)*w2+w1*w3)/denx            
        impdiag(i,4,3)=(-lscqm1(2,2)*w1+w2*w3)/denx
        impdiag(i,4,4)=(lscqm1(2,2)*lscqm1(3,3)+w3**2)/denx
        impdiag(i,5,5)=1.0d0/lscqm1(5,5)
    end if
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
impdiagt(i,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(i,1:turbulenceequations+passivescalar),1.0d-30)
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if


if (relax.eq.1)then
do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
    if(mrf.eq.1)then
        srf=rec_mrf(i)
    end if
if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))
end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)
end if
end if
if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if



		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)

		if ((turbulence.eq.1).or.(passivescalar.gt.0))then
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if

		b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
		if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		b1t(:)=b1t(:)-(impofft(i,l,:)*dut1(:))
		end if
end do	!loop f


else

do l=1,ielem_ifca(i)	!loop3
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
               if (rec_mrf(i).eq.1)then
                    !retrieve rotational velocity in case of rotating reference frame to calculate the correct value of the boundary condition
                    srf_speed(2:4)=rec_rotvel(l,1,1:3,i)
                    call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
                end if

					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
								    if(per_rot.eq.1)then
                                        du1(2:4)=rotate_per_1(du1(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
								    end if
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then

								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

									end if



								  else
								  !not periodic ones in my cpu

								  facex=l;iconsidered=i
								  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

								   if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if

								    cords(1:3)=zero
								    cords(1:3)=cordinates3(n,nodes_list,n_node)

								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    poz(1)=cords(3)

								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))




								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if


								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if

								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if


								    case(6,9,99)

								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if


								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if



								    end if


								    end select




								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then

								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

									end if




							end if
					    else	!in other cpus they can only be periodic or mpi neighbours



							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then	!periodic in other cpu

! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
                                if(per_rot.eq.1)then
                                    du1(2:4)=rotate_per_1(du1(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
                                end if
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then

! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)

								 end if

								end if
							else
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)

								if ((turbulence.gt.0).or.(passivescalar.gt.0))then

! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)

								 end if


!
							end if
					    end if



b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(:)=b1t(:)-(impofft(i,l,:)*dut1(:))
end if
end do	!loop f
end if
impdu(i,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(i,1:turbulenceequations+passivescalar)*b1t(1:turbulenceequations+passivescalar)
end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)


end do!sweeps


else


do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(:)=zero
end if

if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))
! dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))

end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)
! dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)
! dummy12t(:)=zero
end if
end if



if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_reorient(l,i).eq.0)then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3	
                                 igoflux=0
		     if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                        if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                                icaseb=1        !periodic mine
                                                        else
                                                                icaseb=3        !physical
                                                        end if
                                                   
                                                    else
                                                    
                                                                icaseb=2!no boundaries interior
                                                    
                                                    end if
                                else
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                                                icaseb=4
                                                                end if
                                                    else
                                                                icaseb=5
                                                    
                                                    end if
                                
                                
                                end if
				  if (icaseb.le.2)then
                                        if (ielem_reorient(l,i).eq.0)then
                                            igoflux=1
                                        else
                                            igoflux=0
                                        end if
                                else
                                        igoflux=2
                                end if
                                 if (ielem_reorient(l,i).eq.0)then
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

								   if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if

								    cords(1:3)=zero
								    cords(1:3)=cordinates3(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    

								    
								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f
end if
impdu(i,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(i,1:turbulenceequations+passivescalar)*&
(b1t(1:turbulenceequations+passivescalar)-dummy12t(1:turbulenceequations+passivescalar))


end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)
 


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe,-1	!loop2
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(:)=zero
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_reorient(l,i).eq.1)then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
+(impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
end if
		end if
end do	!loop f
			
else

do l=1,ielem_ifca(i)	!loop3	

                          igoflux=0
		     if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                        if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                                icaseb=1        !periodic mine
                                                        else
                                                                icaseb=3        !physical
                                                        end if
                                                   
                                                    else
                                                    
                                                                icaseb=2!no boundaries interior
                                                    
                                                    end if
                                else
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                                                icaseb=4
                                                                end if
                                                    else
                                                                icaseb=5
                                                    
                                                    end if
                                
                                
                                end if
				  if (icaseb.le.2)then
                                        if (ielem_reorient(l,i).eq.1)then
                                            igoflux=1
                                        else
                                            igoflux=0
                                        end if
                                else
                                        igoflux=2
                                end if
                                 if (ielem_reorient(l,i).eq.1)then

				
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

 dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
+(impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
end if
		end if
end do	!loop f
end if

impdu(i,1:nof_variables)=impdu(i,1:nof_variables)-matmul(impdiag(i,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then

impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)-&
(impdiagt(i,1:turbulenceequations+passivescalar)*dummy12t(1:turbulenceequations+passivescalar))

end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 
 call exhboundhigher2(n)
 

end do!sweeps

end if







end subroutine relaxation


subroutine relaxation_lm(n)
 !> @brief
!> this subroutine solves the linear system for implicit time stepping either through jacobian or lu-sgs in 3d with low memory footprint
implicit none
integer,intent(in)::n
integer::i,l,k,ii,sweeps,kmaxe,nvar,ibfc,igoflux
real::impres1,impres2,impres3
real,dimension(1:nof_variables,1:nof_variables)::lscqm1
real::b1_imp(1:nof_variables),du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::b1t(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::iconsidered, facex, pointx
integer::ngp,n_node
integer::b_code,nfx,lfx,rowfx
real::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv,srf_speedrot,srf_speed
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords
real::impdiagt(1,turbulenceequations+passivescalar)
real::impdiag(1,1:nof_variables,1:nof_variables),impofft(1,6,turbulenceequations+passivescalar)
real::impoff(1,6,1:nof_variables,1:nof_variables)









sweeps=10
kmaxe=xmpielrank(n)

impdu(iconsidered,:)=zero


du1=zero
b1_imp=zero
lscqm1=zero
dur=zero; dul=zero
durr=zero; dulr=zero







if (relax.eq.1)then
do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
iconsidered=i
call calculate_jacobianlm(n,iconsidered,impdiag,impdiagt,impoff,impofft)
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(1,1:nof_variables,1:nof_variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)
impdiag(1,5,5)=1.0d0/lscqm1(5,5)



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdiagt(1,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(1,1:turbulenceequations+passivescalar),1.0d-30)
end if

if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))

end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)

end if
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
		if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		b1t(:)=b1t(:)-(impofft(1,l,:)*dut1(:))
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3			
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

								   if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if

								    cords(1:3)=zero
								    cords(1:3)=cordinates3(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    

								    
								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1,6)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(:)=b1t(:)-(impofft(1,l,:)*dut1(:))
end if
end do	!loop f
end if

impdu(i,1:nof_variables)=matmul(impdiag(1,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(1,1:turbulenceequations+passivescalar)*b1t(1:turbulenceequations+passivescalar)
end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)


end do!sweeps


else


do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
iconsidered=i
call calculate_jacobianlm(n,iconsidered,impdiag,impdiagt,impoff,impofft)
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(1,1:nof_variables,1:nof_variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)
impdiag(1,5,5)=1.0d0/lscqm1(5,5)



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdiagt(1,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(1,1:turbulenceequations+passivescalar),1.0d-30)
end if




if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))
dummy12t(:)=zero
end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)
dummy12t(:)=zero
end if
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_reorient(l,i).eq.0)then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3	
				if (ielem_reorient(l,i).eq.0)then
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

								   if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if


								    cords(1:3)=zero
								    cords(1:3)=cordinates3(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    

								    
								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								 select case(b_code)
								    case(1)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f
end if

impdu(i,1:nof_variables)=matmul(impdiag(1,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(1,1:turbulenceequations+passivescalar)*&
(b1t(1:turbulenceequations+passivescalar)-dummy12t(1:turbulenceequations+passivescalar))


end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)
 


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe,-1	!loop2

iconsidered=i
call calculate_jacobianlm(n,iconsidered,impdiag,impdiagt,impoff,impofft)
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(1,1:nof_variables,1:nof_variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)
impdiag(1,5,5)=1.0d0/lscqm1(5,5)



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdiagt(1,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(1,1:turbulenceequations+passivescalar),1.0d-30)
end if







dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(:)=zero
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_reorient(l,i).eq.1)then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
+(impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
end if
		end if
end do	!loop f

				
else

do l=1,ielem_ifca(i)	!loop3	
				if (ielem_reorient(l,i).eq.1)then
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

 dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
+(impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
end if
		end if
end do	!loop f
end if
impdu(i,1:nof_variables)=impdu(i,1:nof_variables)-matmul(impdiag(1,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then

impdu(1,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)&
-((impdiagt(1,1:turbulenceequations+passivescalar)*dummy12t(1:turbulenceequations+passivescalar)))

end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 
 call exhboundhigher2(n)
 

end do!sweeps

end if







end subroutine relaxation_lm


subroutine relaxation2d(n)
 !> @brief
!> this subroutine solves the linear system for implicit time stepping either through jacobian or lu-sgs in 2d
implicit none
integer,intent(in)::n
integer::i,l,k,ii,sweeps,kmaxe,nvar,igoflux, icaseb,n_node,ibfc,j
real::impres1,impres2,impres3
real,dimension(1:nof_variables,1:nof_variables)::lscqm1
real::b1_imp(1:nof_variables),du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::b1t(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::iconsidered, facex, pointx
integer::ngp
integer::b_code,nfx,lfx,rowfx
real::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv,srf_speedrot,srf_speed
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords

sweeps=10
kmaxe=xmpielrank(n)

impdu(:,:)=zero


du1=zero
b1_imp=zero
lscqm1=zero
dur=zero; dul=zero
durr=zero; dulr=zero



call calculate_jacobian_2d(n)

#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(i,1:nof_variables,1:nof_variables)
  do j=1,nof_variables

	impdiag(i,j,j)=1.0d0/lscqm1(j,j)
  end do
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
impdiagt(i,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(i,1:turbulenceequations+passivescalar),1.0d-30)
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if


if (relax.eq.1)then
do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2

if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))
end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)
end if
end if
if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
		if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		b1t(:)=b1t(:)-(impofft(i,l,:)*dut1(:))
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3			
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=angle1
				ny=angle2
				
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
								  n_node=2
								    cords(1:2)=zero
								    cords(1:2)=cordinates2(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    

								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    



								    
								     select case(b_code)
								    case(1)
								    du1(:)=zero
								    
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								   
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if


								

								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	

b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(:)=b1t(:)-(impofft(i,l,:)*dut1(:))
end if
end do	!loop f
end if
impdu(i,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(i,1:turbulenceequations+passivescalar)*b1t(1:turbulenceequations+passivescalar)
end if
	!loop elements

end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)

 

end do!sweeps


else


do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(:)=zero
end if

if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))
! dummy12t(:)=zero
end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)
! dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)
! dummy12t(:)=zero
end if
end if



if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_ineigh(l,i).lt.i) then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3	
				  igoflux=0
		     if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                        if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                                icaseb=1        !periodic mine
                                                        else
                                                                icaseb=3        !physical
                                                        end if
                                                   
                                                    else
                                                    
                                                                icaseb=2!no boundaries interior
                                                    
                                                    end if
                                else
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                                                icaseb=4
                                                                end if
                                                    else
                                                                icaseb=5
                                                    
                                                    end if
                                
                                
                                end if
				  if (icaseb.le.2)then
                                        if (ielem_ineigh(l,i).lt.i)then
                                            igoflux=1
                                        else
                                            igoflux=0
                                        end if
                                else
                                        igoflux=2
                                end if
                                 if (ielem_ineigh(l,i).lt.i)then
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=angle1
				ny=angle2
				
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
								  n_node=2
								    cords(1:2)=zero
								    cords(1:2)=cordinates2(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    

								    
								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f
end if
impdu(i,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(i,1:turbulenceequations+passivescalar)*&
(b1t(1:turbulenceequations+passivescalar)-dummy12t(1:turbulenceequations+passivescalar))


end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)
 


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(:)=zero
end if

if (ielem_interior(i).eq.0)then
                        do l=1,ielem_ifca(i)	!loop3
                                    if (ielem_ineigh(l,i).gt.i)then
                                                    du1(1:nof_variables)=zero
                                                    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                    dut1(:)=zero
                                                    end if
                                                        
                                                        
                                
                                        du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
                                                                                            
                                                    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
                                                    dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                    end if	
                                        
                                                                    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1))
                                                if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
                                                +(impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
                                                end if
                                        end if
                        end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3	
					  igoflux=0
                                    
                                    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                                        if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                                                icaseb=1        !periodic mine
                                                                        else
                                                                                icaseb=3        !physical
                                                                        end if
                                                                
                                                                    else
                                                                    
                                                                                icaseb=2!no boundaries interior
                                                                    
                                                                    end if
                                                else
                                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                                                if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                                                                icaseb=4
                                                                                end if
                                                                    else
                                                                                icaseb=5
                                                                    
                                                                    end if
                                                
                                                
                                                end if
                                                if (icaseb.le.2)then
                                                        if (ielem_reorient(l,i).eq.1)then
                                                            igoflux=1
                                                        else
                                                            igoflux=0
                                                        end if
                                                else
                                                        igoflux=2
                                                end if
                                                                    if (ielem_ineigh(l,i).gt.i)then
                                                                    angle1=ielem_faceanglex(l,i)
                                                                    angle2=ielem_faceangley(l,i)
                                                                    nx=angle1
                                                                    ny=angle2
                                                                    
                                            
                                                                                        if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                                                                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                                                                                        if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                                                                                            du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
                                                                                                                                if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                                                                                            
                                                                                                                            dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                                                                                                
                                                                                                                                end if
                                                                                                                                                                
                                                                                                                        
                                                                                                                        
                                                                                                                        else
                                                                                                                        
                                                                                                                            
                                                                                                                            
                                                                                                                                                                                        
                                                                                                                            
                                                                                                                        end if
                                                                                                                else
                                                                                                                    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
                                                                                                                                if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                                                                                            
                                                                                                                            dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                                                                                                
                                                                                                                                end if
                                                                                                                    
                                                                                                                    
                                                                                                                    
                                                                                                                    
                                                                                                                end if
                                                                                            else	!in other cpus they can only be periodic or mpi neighbours
                                                                                            
                                                                                            
                                                                                                
                                                                                                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                                                                                if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                                                                                                
!                                                                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
                                                                                                                nfx  = ielem_ineighn(l,i)
                                                                                                                lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
                                                                                                                rowfx = bound_offset(nfx) + lfx - 1
                                                                                                                du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
                                                                    
                                                                                                                if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                                                                                
!                                                                                                                 dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                                                                                nfx  = ielem_ineighn(l,i)
                                                                                                                lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
                                                                                                                rowfx = bound_offset(nfx) + lfx - 1
                                                                                                                dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                                                                                
                                                                                                                end if
                                                                                        
                                                                                                                end if
                                                                                                        else 			
!                                                                                                             du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
                                                                                                            nfx  = ielem_ineighn(l,i)
                                                                                                            lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
                                                                                                            rowfx = bound_offset(nfx) + lfx - 1
                                                                                                            du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
                                                                    
                                                                                                                if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                                                                                
!                                                                                                                 dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                                                                                nfx  = ielem_ineighn(l,i)
                                                                                                                lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
                                                                                                                rowfx = bound_offset(nfx) + lfx - 1
                                                                                                                dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                                                                                                
                                                                                                                end if
                                                                                                        
                                                                                                                
                                                ! 								   
                                                                                                                end if
                                                                                                    end if
                                                    
                                                    

                                            dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(i,l,1:nof_variables,1:nof_variables),du1))
                                                                                            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                                                                                            dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
                                                                                            +(impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
                                                                                            end if
                                                            end if
end do	!loop f
end if

! impdu(i,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),impdu(i,1:nof_variables)-dummy12(1:nof_variables))


impdu(i,1:nof_variables)=impdu(i,1:nof_variables)-matmul(impdiag(i,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))


if ((turbulence.gt.0).or.(passivescalar.gt.0))then

impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)-&
(impdiagt(i,1:turbulenceequations+passivescalar)*dummy12t(1:turbulenceequations+passivescalar))

end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 
 call exhboundhigher2(n)
 

end do!sweeps

end if







end subroutine relaxation2d


subroutine relaxation_lm2d(n)
 !> @brief
!> this subroutine solves the linear system for implicit time stepping either through jacobian or lu-sgs in 2d with low memory footprint
implicit none
integer,intent(in)::n
integer::i,l,k,ii,sweeps,kmaxe,nvar,ibfc,igoflux
real::impres1,impres2,impres3
real,dimension(1:nof_variables,1:nof_variables)::lscqm1
real::b1_imp(1:nof_variables),du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::b1t(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::iconsidered, facex, pointx,n_node
integer::ngp
integer::b_code,nfx,lfx,rowfx
real::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv,srf_speedrot,srf_speed
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords
real::impdiagt(1,turbulenceequations+passivescalar)
real::impdiag(1,1:nof_variables,1:nof_variables),impofft(1,4,turbulenceequations+passivescalar)
real::impoff(1,4,1:nof_variables,1:nof_variables)




sweeps=10
kmaxe=xmpielrank(n)

impdu(iconsidered,:)=zero


du1=zero
b1_imp=zero
lscqm1=zero
dur=zero; dul=zero
durr=zero; dulr=zero







if (relax.eq.1)then
do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
iconsidered=i
call calculate_jacobian_2dlm(n,iconsidered,impdiag,impdiagt,impoff,impofft)
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(1,1:nof_variables,1:nof_variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdiagt(1,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(1,1:turbulenceequations+passivescalar),1.0d-30)
end if

if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))

end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)

end if
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
		if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		b1t(:)=b1t(:)-(impofft(1,l,:)*dut1(:))
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3			
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=angle1
				ny=angle2
				
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
								  n_node=2
								    cords(1:2)=zero
								    cords(1:2)=cordinates2(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								   
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    
								    
								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								     select case(b_code)
								    case(1)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    case(6)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

b1_imp(1:nof_variables)=b1_imp(1:nof_variables)-matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(:)=b1t(:)-(impofft(1,l,:)*dut1(:))
end if
end do	!loop f
end if
impdu(i,1:nof_variables)=matmul(impdiag(1,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(1,1:turbulenceequations+passivescalar)*b1t(1:turbulenceequations+passivescalar)
end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)


end do!sweeps


else


do ii=1,sweeps	!loop1
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe	!loop2
iconsidered=i
call calculate_jacobian_2dlm(n,iconsidered,impdiag,impdiagt,impoff,impofft)
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(1,1:nof_variables,1:nof_variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)




if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdiagt(1,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(1,1:turbulenceequations+passivescalar),1.0d-30)
end if




if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
do nvar=1,turbulenceequations+passivescalar
b1t(nvar)=-(rhst_val(nvar,i)+((((1.5d0*u_ct_val(1,nvar,i))-(2.0d0*u_ct_val(2,nvar,i))+(0.5d0*u_ct_val(3,nvar,i)))/(dt))*ielem_totvolume(i)))
dummy12t(:)=zero
end do
end if
else
b1_imp(1:nof_variables)=-rhs_val(1:nof_variables,i)
dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)
dummy12t(:)=zero
end if
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_reorient(l,i).eq.0)then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3	
				if (ielem_reorient(l,i).eq.0)then
				angle1=ielem_faceanglex(l,i)
				angle2=ielem_faceangley(l,i)
				nx=angle1
				ny=angle2
				
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
								  n_node=2
								    cords(1:2)=zero
								    cords(1:2)=cordinates2(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								   
								    
								    leftv(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    cturbl(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    b_code=ibound_icode(ielem_ibounds(l,i))
								   

								    
								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    
								    du1(1:nof_variables)=rightv(1:nof_variables)
				  				    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    case(6)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
		  dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)+&
		  (impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
		  end if
		end if
end do	!loop f
end if

impdu(i,1:nof_variables)=matmul(impdiag(1,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdiagt(1,1:turbulenceequations+passivescalar)*&
(b1t(1:turbulenceequations+passivescalar)-dummy12t(1:turbulenceequations+passivescalar))


end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigher2(n)
 


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe,-1	!loop2

iconsidered=i
call calculate_jacobian_2dlm(n,iconsidered,impdiag,impdiagt,impoff,impofft)
  lscqm1(1:nof_variables,1:nof_variables)=impdiag(1,1:nof_variables,1:nof_variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)




if ((turbulence.gt.0).or.(passivescalar.gt.0))then
impdiagt(1,1:turbulenceequations+passivescalar)=1.0d0/max(impdiagt(1,1:turbulenceequations+passivescalar),1.0d-30)
end if







dummy12(:)=zero
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(:)=zero
end if

if (ielem_interior(i).eq.0)then
do l=1,ielem_ifca(i)	!loop3
	      if (ielem_reorient(l,i).eq.1)then
			    du1(1:nof_variables)=zero
			    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
			    dut1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
				     				      
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
		dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		end if	
		
		    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
+(impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
end if
		end if
end do	!loop f


				
else

do l=1,ielem_ifca(i)	!loop3	
				if (ielem_reorient(l,i).eq.1)then
				
	
					if (ielem_ineighb(l,i).eq.n)then	!my cpu only
						if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								    du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
								  					  
								  
								  
								  else
								  
								    
								    
				  				  				  				  
								    
								  end if
							else
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									
									end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
								
! 								du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
								nfx  = ielem_ineighn(l,i)
								lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								rowfx = bound_offset(nfx) + lfx - 1
								du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 			
! 							      du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
							      nfx  = ielem_ineighn(l,i)
							      lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
							      rowfx = bound_offset(nfx) + lfx - 1
							      du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 								  dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
								  
								 end if
							
								  
! 								   
							end if
					    end if
	
	

 dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(impoff(1,l,1:nof_variables,1:nof_variables),du1))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
dummy12t(1:turbulenceequations+passivescalar)=dummy12t(1:turbulenceequations+passivescalar)&
+(impofft(1,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar))
end if
		end if
end do	!loop f
end if
impdu(i,1:nof_variables)=impdu(i,1:nof_variables)-matmul(impdiag(1,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))
if ((turbulence.gt.0).or.(passivescalar.gt.0))then

impdu(1,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)&
-((impdiagt(1,1:turbulenceequations+passivescalar)*dummy12t(1:turbulenceequations+passivescalar)))

end if
end do	!loop elements
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 
 call exhboundhigher2(n)
 

end do!sweeps

end if







end subroutine relaxation_lm2d






subroutine relaxation_lumfree(n)
implicit none
!> @brief
!> this subroutine solves the linear system for implicit time stepping either through matrix free lu-sgs low memory footprint
integer,intent(in)::n
integer::i,l,k,ii,sweeps,kmaxe,nvar,igoflux,icaseb,inds,modeu
real,dimension(1:nof_variables,1:nof_variables)::lscqm1
real::b1_imp(1:nof_variables),du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::b1t(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::iconsidered, facex, pointx
integer::ngp
integer::b_code,nfx,lfx,rowfx
real::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot,srf_speedrot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords

sweeps=1
kmaxe=xmpielrank(n)
impdu(:,:)=zero

if (turbulence.gt.0)then
dut1=zero
end if


if (dimensiona.eq.3)then
inds=5
call calculate_jacobian_3d_mf(n)

else
inds=4
call calculate_jacobian_2d_mf(n)
end if




do ii=1,sweeps	!loop1
du1=zero
call exhboundhigher2(n)

#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe   !forward sweep !loop 1
iconsidered=i
                    b1_imp=zero
                    if (turbulence.gt.0)b1t=zero
                    modeu=1
                     if (ielem_interior(i).eq.0)then  
                    call lusgs_interior(n,iconsidered,modeu,b1_imp,b1t)
                     else
                    
                    call lusgs_out(n,iconsidered,modeu,b1_imp,b1t)
                     end if
                     
                    if (iscoun.ne.1)then
b1_imp(1:nof_variables)=-(rhs_val(1:nof_variables,i)+((((1.5*u_c_val(1,1:nof_variables,i))-(2.0d0*u_c_val(2,1:nof_variables,i))+(0.5d0*u_c_val(3,1:nof_variables,i)))/(dt))*ielem_totvolume(i)))-b1_imp(1:nof_variables)
                    else
			
                    b1_imp(1:nof_variables)=(-rhs_val(1:nof_variables,i))-b1_imp(1:nof_variables)
                    end if
                    
                    
                    impdu(i,1:nof_variables)=(b1_imp(1:nof_variables)/impdiag_mf(i))
                    
                    
                    
                    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                    
                    if (iscoun.ne.1)then
                    b1t(1:turbulenceequations+passivescalar)=-(rhst_val(1:turbulenceequations+passivescalar,i)+((((1.5d0*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(2.0d0*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.5d0*u_ct_val(3,1:turbulenceequations+passivescalar,i)))/(dt))*ielem_totvolume(i)))-b1t(1:turbulenceequations+passivescalar)
                    impdu(i,inds+1:inds+turbulenceequations+passivescalar)=b1t(:)/impdiagt(i,1:turbulenceequations+passivescalar)
                    
                    
                    else
                    b1t(1:turbulenceequations+passivescalar)=-rhst_val(1:turbulenceequations+passivescalar,i)-b1t(1:turbulenceequations+passivescalar)
                    impdu(i,inds+1:inds+turbulenceequations+passivescalar)=b1t(:)/impdiagt(i,1:turbulenceequations+passivescalar)
                    end if
                    
                    
                        end if
                    
                    
                  
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 call exhboundhigherlu(n)

#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=kmaxe,1,-1     !forward sweep !loop 1
iconsidered=i
                     b1_imp=zero
                    modeu=0
                    if (turbulence.gt.0)b1t=zero
                     if (ielem_interior(i).eq.0)then  !if 1
                        call lusgs_interior(n,iconsidered,modeu,b1_imp,b1t)
                     else
                    
                    call lusgs_out(n,iconsidered,modeu,b1_imp,b1t)
                     end if

                   
                    
!                     impdu(i,1:nof_variables)=impdu(i,1:nof_variables)-(b1_imp(1:nof_variables)/impdiag_mf(i))
                    
                    
                     impdu(i,1:nof_variables)=impdu(i,1:nof_variables)-(b1_imp(1:nof_variables)/impdiag_mf(i))
                    
                    
                     if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                    impdu(i,inds+1:inds+turbulenceequations+passivescalar)=impdu(i,inds+1:inds+turbulenceequations+passivescalar)-b1t(1:turbulenceequations+passivescalar)/impdiagt(i,1:turbulenceequations+passivescalar)
                        end if
                    
                    
                 
end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end do!sweeps
                                       
end subroutine relaxation_lumfree                  
                    
                    
subroutine lusgs_interior(n,iconsidered,modeu,b1_imp,b1t)
implicit none
!> @brief
!> this subroutine solves the linear system for implicit time stepping either through matrix free lu-sgs low memory footprint
integer,intent(in)::n,iconsidered,modeu
integer::i,l,k,ii,inds
real::impres1,impres2,impres3
real,intent(inout)::b1_imp(1:nof_variables),b1t(turbulenceequations+passivescalar)
real::du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::facex, pointx,igoflux
integer::ngp
integer::b_code,nfx,lfx,rowfx
real::angle1,angle2,nx,ny,nz,velnormal
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot,srf_speedrot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords
real::deltaf(1:nof_variables+turbulenceequations+passivescalar)




i=iconsidered

!     modeu=1!lower
!     modeu=0!upper
if (dimensiona.eq.3)then
inds=5
else
inds=4
end if

do l=1,ielem_ifca(i)	!loop2     
                                    if (modeu.eq.0)then
                                    
							       if (ielem_ineigh(l,i).lt.iconsidered) cycle
							       
							       else
							       
							        if (ielem_ineigh(l,i).gt.iconsidered) cycle
							        
							       end if
        
        
        
         du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
         cright(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
              if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
                  dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),inds+1:inds+turbulenceequations+passivescalar)   
                  cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
                  cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
               end if	
        angle1=ielem_faceanglex(l,i)
        angle2=ielem_faceangley(l,i)
        nx=angle1
        ny=angle2
        facex=l;
        call deltafx(iconsidered,facex,deltaf,velnormal)
        
        
        
        
         b1_imp(1:nof_variables)=b1_imp(1:nof_variables)+(0.5d0*(deltaf(1:nof_variables)-velnormal*du1(1:nof_variables))*ielem_surf(l,i))
        
!         
       
        
        
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        
!         b1t(1:turbulenceequations+passivescalar)=b1t(1:turbulenceequations+passivescalar)+impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar)
        
        b1t(1:turbulenceequations+passivescalar)=b1t(1:turbulenceequations+passivescalar)+0.5d0*(deltaf(inds+1:inds+turbulenceequations+passivescalar)-velnormal*dut1(1:turbulenceequations+passivescalar))*ielem_surf(l,i)
        end if
        
        

        
end do
        
        

        
        
        
end subroutine lusgs_interior

                    
subroutine lusgs_out(n,iconsidered,modeu,b1_imp,b1t)
implicit none
!> @brief
!> this subroutine solves the linear system for implicit time stepping either through matrix free lu-sgs low memory footprint
integer,intent(in)::n,iconsidered,modeu
integer::i,l,k,ii,inds,interfaces,n_node
real,intent(inout)::b1_imp(1:nof_variables),b1t(turbulenceequations+passivescalar)
real::du1(1:nof_variables),du2(1:nof_variables),dummy12(1:nof_variables),c1_imp(1:nof_variables)
real::dur(nof_variables),dul(nof_variables)
real::durr(nof_variables),dulr(nof_variables)
real::dut1(turbulenceequations+passivescalar)
real::dummy12t(turbulenceequations+passivescalar)
integer::facex, pointx,igoflux
integer::ngp
integer::b_code
real::angle1,angle2,nx,ny,nz,velnormal
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:nof_variables)::leftv,srf_speedrot,srf_speed
real,dimension(1:nof_variables)::rightv
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:8,1:dimensiona)::vext,nodes_list
real,dimension(1:dimensiona)::cords
integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf
real::deltaf(1:nof_variables+turbulenceequations+passivescalar)
i=iconsidered

if (dimensiona.eq.3)then
inds=nof_variables
else
inds=nof_variables
end if

do l=1,ielem_ifca(i)	!cond1   
        b_code=0
        
        
        
        if (ielem_ineighb(l,i).eq.n)then	!my cpu only !cond2
        
		if (ielem_ibounds(l,i).gt.0)then	!check for boundaries !cond3
		if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu !cond4
		
		
                                    if (modeu.eq.0)then
							       if (ielem_ineigh(l,i).lt.iconsidered)cycle
							       else
							        if (ielem_ineigh(l,i).gt.iconsidered)cycle
							       end if

		
		
		
		
		 du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
		 
		 cright(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
		 
		 
		 
		 
		if ((turbulence.gt.0).or.(passivescalar.gt.0))then  
        dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),inds+1:inds+turbulenceequations+passivescalar)
        cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
        end if      
        else                                !cond4
		!not periodic ones in my cpu
        facex=l;
        
        
        
        if (dimensiona.eq.2)then
		call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)

		n_node=2
		cords(1:2)=zero
		cords(1:2)=cordinates2(n,nodes_list,n_node);pox(1)=cords(1);poy(1)=cords(2)
		else
		call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

		 if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if

		cords(1:3)=zero
		cords(1:3)=cordinates3(n,nodes_list,n_node);pox(1)=cords(1);poy(1)=cords(2);poz(1)=cords(3);
		
		end if
		
		
        leftv(1:nof_variables)=impdu(i,1:nof_variables)!u_c_val(1,1:nof_variables,iconsidered)						           
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then  !cond5
		 cturbl(1:turbulenceequations+passivescalar)=impdu(i,inds+1:inds+turbulenceequations+passivescalar)
		end if                                                !end cond5
		b_code=ibound_icode(ielem_ibounds(l,i))
		
		if (dimensiona.eq.2)then
			


		call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
		else
		
		call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
		end if
		
		du1(1:nof_variables)=rightv(1:nof_variables)
		
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then !cond6
		 dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
		end if                                              !end cond6  
                                    select case(b_code)     !cond 7
                                    case(1)
                                    du1(:)=zero
								    
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								   
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,inds+1:inds+turbulenceequations+passivescalar)
								    end if
								    
								    case(4)
        
								    
								    
								    case(6)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=impdu(i,1:nof_variables)
								    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								    dut1(1:turbulenceequations+passivescalar)=impdu(i,inds+1:inds+turbulenceequations+passivescalar)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select   !end cond7  
								    
								    
				  				  				  				  
								    
								  end if            !end cond4
							else     !cond3
                                    
                                    
                                     if (modeu.eq.0)then
							       if (ielem_ineigh(l,i).lt.iconsidered)cycle
							       else
							        if (ielem_ineigh(l,i).gt.iconsidered)cycle
							       end if
                                    
                                    
							       du1(1:nof_variables)=impdu(ielem_ineigh(l,i),1:nof_variables)
							       cright(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
							       
							       if (modeu.eq.0)then
							       if (ielem_ineigh(l,i).lt.iconsidered)cycle
							       else
							        if (ielem_ineigh(l,i).gt.iconsidered)cycle
							       end if
							       
							       
							       
									if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								      
								      dut1(1:turbulenceequations+passivescalar)=impdu(ielem_ineigh(l,i),inds+1:inds+turbulenceequations+passivescalar)
									cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
							      
							      
							      
							      
							end if   !end cond3
					    else	!in other cpus they can only be periodic or mpi neighbours !cond2
					    
					    
						
	if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
	
			if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                
                                    
!
! 			cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 			(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)
			nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
			lf=rec_ihexl(1,ielem_indexi(L,i),i)
			rowf=halo_offset(nf) + lf - 1
			cright(1:nof_variables)=solhir(rowf,1:nof_variables)

! 			du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
			nfx  = ielem_ineighn(l,i)
			lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
			rowfx = bound_offset(nfx) + lfx - 1
			du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)
	      	      
			if ((turbulence.gt.0).or.(passivescalar.gt.0))then
! 			cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 							      (rec_ihexl(1,ielem_indexi(l,i)),inds+1:inds+turbulenceequations+passivescalar)

			  nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
			lf=rec_ihexl(1,ielem_indexi(L,i),i)
				rowf=halo_offset(nf) + lf - 1
				cturbr(1:turbulenceequations+passivescalar)=solhir(rowf,nof_variables+1:nof_variables+turbulenceequations+passivescalar)


! 			dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),inds+1:inds+turbulenceequations+passivescalar)
			nfx  = ielem_ineighn(l,i)
			lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
			rowfx = bound_offset(nfx) + lfx - 1
			dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,inds+1:inds+turbulenceequations+passivescalar)
								  
								 end if
					
								end if
							else 	
							
                                    
                                    
! 				du1(1:nof_variables)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),1:nof_variables)
				nfx  = ielem_ineighn(l,i)
				lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
				rowfx = bound_offset(nfx) + lfx - 1
				du1(1:nof_variables) = boundhiri(rowfx,1:nof_variables)

				nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
			lf=rec_ihexl(1,ielem_indexi(L,i),i)
			rowf=halo_offset(nf) + lf - 1
			cright(1:nof_variables)=solhir(rowf,1:nof_variables)

!
! 	      	      cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 							      (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)
								if ((turbulence.gt.0).or.(passivescalar.gt.0))then
								  
! 				dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(l,i))%facesol(ielem_qface(l,1,i),inds+1:inds+turbulenceequations+passivescalar)
				nfx  = ielem_ineighn(l,i)
				lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
				rowfx = bound_offset(nfx) + lfx - 1
				dut1(1:turbulenceequations+passivescalar) = boundhiri(rowfx,inds+1:inds+turbulenceequations+passivescalar)
! 					cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 							      (rec_ihexl(1,ielem_indexi(l,i)),inds+1:inds+turbulenceequations+passivescalar)

					 nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								cturbr(1:turbulenceequations+passivescalar)=solhir(rowf,nof_variables+1:nof_variables+turbulenceequations+passivescalar)

								 end if
							
								  
! 								   
							end if
					    end if       
               
               
               
               
               
               
               
               
               
        angle1=ielem_faceanglex(l,i)
        angle2=ielem_faceangley(l,i)
        nx=angle1
        ny=angle2
        facex=l;
        call deltafx(iconsidered,facex,deltaf,velnormal)
        
       
        if (b_code.gt.0)then
        if (b_code.ne.5)then
        deltaf=0.0;velnormal=0.0
        end if
        end if

        

        
        
         b1_imp(1:nof_variables)=b1_imp(1:nof_variables)+(0.5d0*(deltaf(1:nof_variables)-velnormal*du1(1:nof_variables))*ielem_surf(l,i))
        

        
        
        
       
        
        
                            
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
!          b1t(1:turbulenceequations+passivescalar)=b1t(1:turbulenceequations+passivescalar)+impofft(i,l,1:turbulenceequations+passivescalar)*dut1(1:turbulenceequations+passivescalar)
        b1t(1:turbulenceequations+passivescalar)=b1t(1:turbulenceequations+passivescalar)+0.5d0*(deltaf(inds+1:inds+turbulenceequations+passivescalar)-velnormal*dut1(1:turbulenceequations+passivescalar))*ielem_surf(l,i)
        end if
        
end do
        

        
end subroutine lusgs_out  



subroutine deltafx(iconsidered,facex,deltaf,velnormal)
implicit none
integer,intent(in)::iconsidered,facex
real::uul,uur,angle1,angle2
real,dimension(1:nof_variables)::cdx1,cdx2
real,intent(inout)::deltaf(1:nof_variables+turbulenceequations+passivescalar),velnormal
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot,srf_speedrot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr,dut1
real,dimension(1:nof_variables)::leftv,du1
real,dimension(1:nof_variables)::rightv
real::mp_pinfl,mp_pinfr,gammal,gammar
integer::inds,nfx,lfx,rowfx


                            if (dimensiona.eq.2)then
                            inds=nof_variables
                            !2 rotate the solution of the neighbour along the face normal and store in cright_rot
                            call rotatef2d(n,cright_rot,cright,angle1,angle2)
                            !3 copy to cleft the solution of the neighbour plus the change du
                            cleft(1:nof_variables)=cright(1:nof_variables)+du1(1:nof_variables)
                            !4 rotate the cleft along the face normal and store in cleft_rot
                            call rotatef2d(n,cleft_rot,cleft,angle1,angle2)
                            !5 compute the spectral radius
                            
                            velnormal=impoff_mf(iconsidered,facex)
                            !copy to temporary variables
                            rightv(1:nof_variables)=cright_rot(1:nof_variables)
                            leftv(1:nof_variables)=cleft_rot(1:nof_variables)
                            
                            !transform leftv and rightv from cons to prim variables
                            call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
                            
                            
                            uul=leftv(2);uur=rightv(2)
                            !evaluate the fluxes
                           
                            cdx2(1:nof_variables)=fluxeval2d(leftv)
                            
                            
                            leftv=rightv
                             cdx1(1:nof_variables)=fluxeval2d(leftv)
                            
                           
                            
                            deltaf(1:nof_variables)=cdx2(1:nof_variables)-cdx1(1:nof_variables)!-velnormal*du1(1:nof_variables)
                           
                           
                            
                            
                            cleft_rot(1:nof_variables)=deltaf(1:nof_variables)
                            call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
                            deltaf(1:nof_variables)=cleft(1:nof_variables)-velnormal*du1(1:nof_variables)
                            
                            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                            
                            deltaf(inds+1:inds+turbulenceequations+passivescalar)=(((cturbr(1:turbulenceequations+passivescalar)+dut1(1:turbulenceequations+passivescalar))*uul)-(cturbr(1:turbulenceequations+passivescalar)*uur))-&
                             velnormal*dut1(1:turbulenceequations+passivescalar)
                            end if
                            
                            else
                            
                            
                            
                            inds=nof_variables
                            !2 rotate the solution of the neighbour along the face normal and store in cright_rot
                            call rotatef(n,cright_rot,cright,angle1,angle2)
                            !3 copy to cleft the solution of the neighbour plus the change du
                            cleft(1:nof_variables)=cright(1:nof_variables)+du1(1:nof_variables)
                            !4 rotate the cleft along the face normal and store in cleft_rot
                            call rotatef(n,cleft_rot,cleft,angle1,angle2)
                            !5 compute the spectral radius
                            
                            velnormal=impoff_mf(iconsidered,facex)
                            !copy to temporary variables
                            rightv(1:nof_variables)=cright_rot(1:nof_variables)
                            leftv(1:nof_variables)=cleft_rot(1:nof_variables)
                            
                            !transform leftv and rightv from cons to prim variables
                            call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
                            
                            
                            uul=leftv(2);uur=rightv(2)
                            !evaluate the fluxes
                           
                            cdx2(1:nof_variables)=fluxeval3d(leftv)
                            
                            
                            leftv=rightv
                             cdx1(1:nof_variables)=fluxeval3d(leftv)
                            
                           
                            
                            deltaf(1:nof_variables)=cdx2(1:nof_variables)-cdx1(1:nof_variables)!-velnormal*du1(1:nof_variables)
                           
                           
                            
                            
                            cleft_rot(1:nof_variables)=deltaf(1:nof_variables)
                            call rotateb(n,cleft,cleft_rot,angle1,angle2)
                            deltaf(1:nof_variables)=cleft(1:nof_variables)-velnormal*du1(1:nof_variables)
                            
                            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
                            
                            deltaf(inds+1:inds+turbulenceequations+passivescalar)=(((cturbr(1:turbulenceequations+passivescalar)+dut1(1:turbulenceequations+passivescalar))*uul)-(cturbr(1:turbulenceequations+passivescalar)*uur))-&
                             velnormal*dut1(1:turbulenceequations+passivescalar)
                            end if
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            end if
                            
                            
                            

end subroutine deltafx



















end module implicit_time
