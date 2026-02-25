module transform
use mpiinfo
use declaration


implicit none

 contains
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine alloc_weights(n)
integer,intent(in)::n
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dimensiona,1:numberofpoints2)::qpoints2d
real,dimension(1:numberofpoints2)::wequa2d



if( dimensiona == 3) then
	allocate(weights_q(1:numberofpoints2),weights_t(1:numberofpoints2))
    call quadraturequad3d(n,igqrules,vext,qpoints2d,wequa2d)
      weights_q(1:qp_quad)=wequa2d(1:qp_quad)

    call quadraturetriang(n,igqrules,vext,qpoints2d,wequa2d)
      weights_t(1:qp_triangle)=wequa2d(1:qp_triangle)
else
	allocate(weights_l(1:numberofpoints2))
    call quadratureline(n,igqrules,vext,qpoints2d,wequa2d)
      weights_l(1:qp_line_n)=wequa2d(1:qp_line_n)
end if

end subroutine alloc_weights

! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 real function trianglearea(n,vext)
!> @brief
!> this function computes the area of a triangle in 3d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:3,1:3)::vva
real,dimension(1:3)::vvb,vvc,vve,vvd,vvjacobsurf
	vvb(1:3)=vext(1,1:3)
	vvc(1:3)=vext(2,1:3)
	vvd(1:3)=vext(3,1:3)
	
	vva(1,1)=vvb(2);vva(2,1)=vvc(2);vva(3,1)=vvd(2)
	vva(1,2)=vvb(3);vva(2,2)=vvc(3);vva(3,2)=vvd(3)
	vva(1,3)=1.0d0;vva(2,3)=1.0d0;vva(3,3)=1.0d0
		vvjacobsurf(1)=vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3)))
	vva(1,1)=vvb(3);vva(2,1)=vvc(3);vva(3,1)=vvd(3)
	vva(1,2)=vvb(1);vva(2,2)=vvc(1);vva(3,2)=vvd(1)
	vva(1,3)=1.0d0;vva(2,3)=1.0d0;vva(3,3)=1.0d0
		vvjacobsurf(2)=vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3)))
	vva(1,1)=vvb(1);vva(2,1)=vvc(1);vva(3,1)=vvd(1)
	vva(1,2)=vvb(2);vva(2,2)=vvc(2);vva(3,2)=vvd(2)
	vva(1,3)=1.0d0; vva(2,3)=1.0d0;	vva(3,3)=1.0d0
		vvjacobsurf(3)=vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3)))
		trianglearea=((oo2)*(sqrt((vvjacobsurf(1)**2)+(vvjacobsurf(2)**2)+(vvjacobsurf(3)**2))))

end function trianglearea



 real function quadarea(n,vext)
 !> @brief
!> this function computes the area of a quadrilateral in 3d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:3)::vve,vvd
	
	vve(1:3)=vext(4,1:3)-vext(2,1:3)
	vvd(1:3)=vext(3,1:3)-vext(1,1:3)
	
	
	quadarea=oo2*sqrt((((vve(2)*vvd(3))-(vve(3)*vvd(2)))**2)+(((vve(3)*vvd(1))-(vve(1)*vvd(3)))**2)+(((vve(1)*vvd(2))-(vve(2)*vvd(1)))**2))
	
	

end function quadarea


 real function linearea(n,vext)
 !> @brief
!> this function computes the length of an edge in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:3)::vve
	
	
	vve(1:2)=vext(2,1:2)-vext(1,1:2)
	
	
	
	linearea=sqrt((vve(1)**2)+(vve(2)**2))
	
	

end function linearea

real function quadvolume(n,vext,qpoints,wequa3d)
!> @brief
!> this function computes the area of quad in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real::s,tx,r,vol
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real,dimension(1:4)::vvxi,vveta,vvnallx,vvnally,vvnxi,vvneta
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1:3)::vvb,vvc,vve,vvd
real,dimension(1)::deta
integer::kk,ii
	
vvxi(1)=-1.0d0; vveta(1)=-1.0d0;
vvxi(2)=1.0d0; vveta(2)=-1.0d0;
vvxi(3)=1.0d0; vveta(3)=1.0d0; 
vvxi(4)=-1.0d0; vveta(4)=1.0d0; 



vvnallx(:)=0.0d0;vvnally(:)=0.0


do kk=1,qp_quad
r=qpoints(1,kk)
s=qpoints(2,kk)


do ii=1,4
vvnxi(1)=-(0.25d0)*(1.0d0-s); vvneta(1)=-(0.25d0)*(1.d0-r);
vvnxi(2)=(0.25d0)*(1.0d0-s); vvneta(2)=-(0.25d0)*(1.d0+r);
vvnxi(3)=(0.25d0)*(1.0d0+s); vvneta(3)=(0.25d0)*(1.d0+r); 
vvnxi(4)=-(0.25d0)*(1.0d0+s); vvneta(4)=(0.25d0)*(1.d0-r); 


vvnallx(ii)=vvnallx(ii)+(vvnxi(ii)*wequa3d(kk))
vvnally(ii)=vvnally(ii)+(vvneta(ii)*wequa3d(kk))

end do
end do



vva=0.0d0
vva1=0.0d0

do ii=1,4

vvnxi(ii)=vvnallx(ii)
vvneta(ii)=vvnally(ii)

    vva(1,1)=vva(1,1)+vvnxi(ii)*vext(ii,1); vva(1,2)=vva(1,2)+vvnxi(ii)*vext(ii,2)
    vva(2,1)=vva(2,1)+vvneta(ii)*vext(ii,1); vva(2,2)=vva(2,2)+vvneta(ii)*vext(ii,2)
    
end do

deta(1)=(vva(1,1)*vva(2,2))-(vva(1,2)*vva(2,1))

vva1(1,1)=(vva(2,2))
vva1(1,2)=-(vva(1,2))
vva1(2,1)=-(vva(2,1))
vva1(2,2)=(vva(1,1))

vol=deta(1)*4.0d0
vva1=vva1/deta(1)


quadvolume=vol



end function quadvolume



real function trianglevolume(n,vext)
!> @brief
!> this function computes the area of triangle in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real::s,tx,r,vol
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1:3)::vvb,vvc,vve,vvd
real,dimension(1)::deta
vol=0.0d0

vva(1,1)=vext(1,1)-vext(3,1)
vva(1,2)=vext(1,2)-vext(3,2)
vva(2,1)=vext(2,1)-vext(3,1)
vva(2,2)=vext(2,2)-vext(3,2)	
vol=(vva(1,1)*vva(2,2))-(vva(2,1)*vva(1,2))


vva1(1,1)=(vva(2,2))
vva1(1,2)=-(vva(1,2))
vva1(2,1)=-(vva(2,1))
vva1(2,2)=(vva(1,1))
vol=vol*0.50d0
vva1=vva1/vol
deta(1)=vol
trianglevolume=vol



end function trianglevolume


! ! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function tetravolume(n,vext)
!> @brief
!> this function computes the volume of a tetrahedrals 

implicit none
#ifdef gpu
!$omp declare target
#endif
!!$omp threadprivate(tetravolume)
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:3,1:3)::vva
real,dimension(1:3)::vvb,vvc,vve,vvd
real,dimension(1:4)::vvjacobvolume



	vvb(1:3)=vext(1,1:3)
	vvc(1:3)=vext(2,1:3)
	vvd(1:3)=vext(3,1:3)
	vve(1:3)=vext(4,1:3)
	
	

	vva(1,1)=vvc(2);vva(2,1)=vvd(2);vva(3,1)=vve(2)
	vva(1,2)=vvc(3);vva(2,2)=vvd(3);vva(3,2)=vve(3)
	vva(1,3)=1.0d0;vva(2,3)=1.0d0;vva(3,3)=1.0d0
	vvjacobvolume(1)=vvb(1)*(vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3))))
	vva(1,1)=vvc(1);vva(2,1)=vvd(1);vva(3,1)=vve(1)
	vva(1,2)=vvc(3);vva(2,2)=vvd(3);vva(3,2)=vve(3)
	vva(1,3)=1.0d0;vva(2,3)=1.0d0;vva(3,3)=1.0d0
	vvjacobvolume(2)=(-vvb(2))*(vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3))))
	vva(1,1)=vvc(1);vva(2,1)=vvd(1);vva(3,1)=vve(1)
	vva(1,2)=vvc(2);vva(2,2)=vvd(2);vva(3,2)=vve(2)
	vva(1,3)=1.0d0;vva(2,3)=1.0d0;vva(3,3)=1.0
	vvjacobvolume(3)=vvb(3)*(vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3))))
	vva(1,1)=vvc(1);vva(2,1)=vvd(1);vva(3,1)=vve(1)
	vva(1,2)=vvc(2);vva(2,2)=vvd(2);vva(3,2)=vve(2)
	vva(1,3)=vvc(3);vva(2,3)=vvd(3);vva(3,3)=vve(3)
	vvjacobvolume(4)=(-1.0)*(vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3)))-vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3)))+vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3))))

	tetravolume=abs((0.166666666666666)*(vvjacobvolume(1)+vvjacobvolume(2)+vvjacobvolume(3)+vvjacobvolume(4)))

	
	
end function tetravolume



subroutine find_angles(n)
integer,intent(in)::n
integer::i,j,k,l,jj,kmaxe

kmaxe=xmpielrank(n)

if (dimensiona.eq.3)then

!$omp do
do i=1,kmaxe
	call find_rot_angles(n,i)
end do
!$omp end do

else

!$omp do
do i=1,kmaxe
	call find_rot_angles2d(n,i)
end do
!$omp end do

end if




end subroutine find_angles


subroutine find_rot_angles(n,iconsi)
!> @brief
!> this subroutine determines the normal vectors for each face
implicit none
real::xc1,yc1,zc1,xc2,yc2,zc2,xc3,yc3,zc3,delxya,delyza,delzxa,delxyb,delyzb,delzxb,delxyc,delyzc,delzxc,nx,ny,nz
real::x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,xx,yy,zz
real::delxa,delxb,delya,delyb,delza,delzb,l_angle1,l_angle2,lnx,lny,lnz,a_rot,b_rot,c_rot,root_rot,anglefacex,anglefacey
integer::k,kmaxe,i,j,kk,kk2,ixf4,ixfv,n_node
integer,intent(in)::n,iconsi
real,dimension(1:8,1:dimensiona)::vext
real:: tempxx
kmaxe=xmpielrank(n)

i=iconsi

			do k=1,ielem_ifca(i)
			       if (ielem_types_faces(k,i).eq.5)then
				    kk2=4;n_node=kk2
			      else
				    kk2=3;n_node=kk2
			      end if
			    if (ielem_interior(i).eq.1)then
			    if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).gt.0))then 	!periodic neighbour

			    xx=ielem_xxc(i)  ;yy=ielem_yyc(i); zz=ielem_zzc(i)

				    do kk=1,n_node
				       if (ielem_reorient(k,i).eq.0)then
				       vext(kk,1:3)=dinoder(ielem_nodes_faces(k,kk,i))%cord(1:3)
				       else



					vext(kk,1:3)=dinoder(ielem_nodes_faces(k,n_node-kk+1,i))%cord(1:3)
!
				       end if
                    if(per_rot.eq.0)then
				      if(abs(vext(kk,1)-xx).gt.xper*oo2)then
				      vext(kk,1)=vext(kk,1)+(xper*sign(1.0,xx-xper*oo2))
				      end if
				      if(abs(vext(kk,2)-yy).gt.yper*oo2)then
				      vext(kk,2)=vext(kk,2)+(yper*sign(1.0,yy-yper*oo2))
				      end if
				      if(abs(vext(kk,3)-zz).gt.zper*oo2)then
				      vext(kk,3)=vext(kk,3)+(zper*sign(1.0,zz-zper*oo2))
				      end if

                    else
                        if (ielem_reorient(k,i).eq.1) then
                            if (ibound_icode(ielem_ibounds(k,i)).eq.5) then
                                tempxx=vext(kk,1)
                                vext(kk,1)=tempxx*cos(-angle_per)-sin(-angle_per)*vext(kk,2)
                                vext(kk,2)=tempxx*sin(-angle_per)+cos(-angle_per)*vext(kk,2)

                            else
                                tempxx=vext(kk,1)
                                vext(kk,1)=tempxx*cos(angle_per)-sin(angle_per)*vext(kk,2)
                                vext(kk,2)=tempxx*sin(angle_per)+cos(angle_per)*vext(kk,2)
                            end if
				      end if
				    end if

			      end do

			    else
				  do kk=1,n_node
					 if (ielem_reorient(k,i).eq.0)then
				       vext(kk,1:3)=dinoder(ielem_nodes_faces(k,kk,i))%cord(1:3)
				       else
					vext(kk,1:3)=dinoder(ielem_nodes_faces(k,n_node-kk+1,i))%cord(1:3)
				       end if
				  end do

			    end if
			    else
				   do kk=1,n_node
					 if (ielem_reorient(k,i).eq.0)then
				       vext(kk,1:3)=dinoder(ielem_nodes_faces(k,kk,i))%cord(1:3)
				       else
					vext(kk,1:3)=dinoder(ielem_nodes_faces(k,n_node-kk+1,i))%cord(1:3)
				       end if
				  end do




			    end if
					if (kk2.eq.3)then
					xc1=vext(1,1); xc2=vext(2,1); xc3=vext(3,1);
					yc1=vext(1,2);yc2=vext(2,2); yc3=vext(3,2);
					zc1=vext(1,3); zc2=vext(2,3); zc3=vext(3,3);
					delxya=(xc1-xc2)*(yc1+yc2);delyza=(yc1-yc2)*(zc1+zc2);delzxa=(zc1-zc2)*(xc1+xc2)
					delxyb=(xc2-xc3)*(yc2+yc3);delyzb=(yc2-yc3)*(zc2+zc3);delzxb=(zc2-zc3)*(xc2+xc3)
					delxyc=(xc3-xc1)*(yc3+yc1);delyzc=(yc3-yc1)*(zc3+zc1);delzxc=(zc3-zc1)*(xc3+xc1)
					nx=(delyza+delyzb+delyzc)
					ny=(delzxa+delzxb+delzxc)
					nz=(delxya+delxyb+delxyc)
					root_rot=sqrt((nx**2)+(ny**2)+(nz**2))
					nx=nx/root_rot; ny=ny/root_rot; nz=nz/root_rot
					root_rot=1.0d0
					a_rot=nx
					b_rot=ny
					c_rot=nz

					call anglex(a_rot,b_rot,anglefacex)
					call angley(c_rot,root_rot,anglefacey)
					ielem_faceanglex(k,i)=anglefacex
					ielem_faceangley(k,i)=anglefacey

					else





					delxa=vext(4,1)-vext(2,1)
					delxb=vext(3,1)-vext(1,1)
					delya=vext(4,2)-vext(2,2)
					delyb=vext(3,2)-vext(1,2)
					delza=vext(4,3)-vext(2,3)
					delzb=vext(3,3)-vext(1,3)


					nx=-0.50d0*((delya*delzb)-(delza*delyb))
					ny=-0.50d0*((delza*delxb)-(delxa*delzb))
					nz=-0.50d0*((delxa*delyb)-(delya*delxb))



					!newells method

					nx=zero;ny=zero;nz=zero;root_rot=zero
 					do kk=1,n_node
 					if (kk.ne.n_node)then
 					nx=nx+(vext(kk,2)-vext(kk+1,2))*(vext(kk,3)+vext(kk+1,3))
 					ny=ny+(vext(kk,3)-vext(kk+1,3))*(vext(kk,1)+vext(kk+1,1))
 					nz=nz+(vext(kk,1)-vext(kk+1,1))*(vext(kk,2)+vext(kk+1,2))
					else
 					nx=nx+(vext(kk,2)-vext(1,2))*(vext(kk,3)+vext(1,3))
 					ny=ny+(vext(kk,3)-vext(1,3))*(vext(kk,1)+vext(1,1))
 					nz=nz+(vext(kk,1)-vext(1,1))*(vext(kk,2)+vext(1,2))

 					end if
 					end do




!
					root_rot=sqrt((nx**2)+(ny**2)+(nz**2))
					nx=nx/root_rot; ny=ny/root_rot; nz=nz/root_rot
					root_rot=1.0d0
					a_rot=nx
					b_rot=ny
					c_rot=nz
					call anglex(a_rot,b_rot,anglefacex)
					call angley(c_rot,root_rot,anglefacey)


					ielem_faceanglex(k,i)=anglefacex
					ielem_faceangley(k,i)=anglefacey
					end if

                    if (ielem_ishape(i).eq.2)then
                    l_angle1=ielem_faceanglex(k,i);l_angle2=ielem_faceangley(k,i)
                    lnx=cos(l_angle1)*sin(l_angle2)
                    lny=sin(l_angle1)*sin(l_angle2)
                    lnz=cos(l_angle2)
                    if ((abs(lnx-1.0d0).le.10e-16).or.(abs(lny-1.0d0).le.10e-16).or.(abs(lnz-1.0d0).le.10e-16))then

                    end if
                    end if


			    end do







end subroutine find_rot_angles



subroutine find_rot_angles2d(n,iconsi)
!> @brief
!> this subroutine determines the normal vectors for each edge
implicit none
real::y1,z1,x2,y2,z2,x3,y3,z3,delxya,delyza,delzxa,delxyb,delyzb,delzxb,delxyc,delyzc,delzxc,nx,ny,nz,anglefacex,anglefacey
real::x5,x6,x7,x8,y5,y6,y7,y8,z5,z6,z7,z8,xx,yy,zz
real::delxa,delxb,delya,delyb,delza,delzb
integer::k,kmaxe,i,j,kk,kk2,ixf4,ixfv,n_node
integer,intent(in)::n,iconsi
real,dimension(1:8,1:dimensiona)::vext

i=iconsi


	i=iconsi

			do k=1,ielem_ifca(i)
! 			       if (ielem_types_faces(k,i).eq.5)then
				    kk2=2;n_node=kk2
! 			      else
! 				    kk2=3;n_node=kk2
! 			      end if
			    if (ielem_interior(i).eq.1)then
			    if ((ielem_ineighg(k,i).gt.0).and.(ielem_ibounds(k,i).gt.0))then 	!periodic neighbour

 			    xx=ielem_xxc(i)  ;yy=ielem_yyc(i); !zz=ielem_zzc(i)

				    do kk=1,n_node
				       if (ielem_reorient(k,i).eq.0)then
				       vext(kk,1:2)=dinoder(ielem_nodes_faces(k,kk,i))%cord(1:2)
				       else

!
					vext(kk,1:2)=dinoder(ielem_nodes_faces(k,n_node-kk+1,i))%cord(1:2)

!
				       end if





				      if(abs(vext(kk,1)-xx).gt.xper*oo2)then
				      vext(kk,1)=vext(kk,1)+(xper*sign(1.0d0,xx-xper/2.0d0))
				      end if
				      if(abs(vext(kk,2)-yy).gt.yper*oo2)then
				      vext(kk,2)=vext(kk,2)+(yper*sign(1.0d0,yy-yper/2.0d0))
				      end if


!



			      end do

			    else
				  do kk=1,n_node
					 if (ielem_reorient(k,i).eq.0)then
				       vext(kk,1:2)=dinoder(ielem_nodes_faces(k,kk,i))%cord(1:2)
				       else
					vext(kk,1:2)=dinoder(ielem_nodes_faces(k,n_node-kk+1,i))%cord(1:2)
				       end if
				  end do

			    end if
			    else
				   do kk=1,n_node
					 if (ielem_reorient(k,i).eq.0)then
				       vext(kk,1:2)=dinoder(ielem_nodes_faces(k,kk,i))%cord(1:2)
				       else
					vext(kk,1:2)=dinoder(ielem_nodes_faces(k,n_node-kk+1,i))%cord(1:2)
				       end if
				  end do




			    end if



					call angle2d(vext,anglefacex,anglefacey)


!
					ielem_faceanglex(k,i)=anglefacex
					ielem_faceangley(k,i)=anglefacey

!




			    end do













!







end subroutine find_rot_angles2d




subroutine computejacobians(n,vext,vva1,deta)
!> @brief
!> this function computes the jacobian of a tetrahedral
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:3,1:3)::vva
real,dimension(1:3,1:3),intent(inout)::vva1
real,dimension(1),intent(inout)::deta


     vva(:,1) = vext(2,:)-vext(1,:)
     vva(:,2) = vext(3,:)-vext(1,:)
     vva(:,3) = vext(4,:)-vext(1,:)
	 deta(1) =  vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3))) - vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3))) + vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3)))


	 vva1(1,1) =  vva(3,3)*vva(2,2) - vva(3,2)*vva(2,3)
	 vva1(1,2) = -(vva(3,3)*vva(1,2) - vva(3,2)*vva(1,3))
	 vva1(1,3) =  vva(2,3)*vva(1,2) - vva(2,2)*vva(1,3)

	 vva1(2,1) = -(vva(3,3)*vva(2,1)-vva(3,1)*vva(2,3) )   
	 vva1(2,2) = vva(3,3)*vva(1,1) -vva(3,1)*vva(1,3)
	 vva1(2,3) = -(vva(2,3)*vva(1,1)-vva(2,1)*vva(1,3)) 

	 vva1(3,1) = vva(3,2)*vva(2,1)-vva(3,1)*vva(2,2)
	 vva1(3,2) =  -(vva(3,2)*vva(1,1)-vva(3,1)*vva(1,2))
	 vva1(3,3) =   vva(2,2)*vva(1,1)-vva(2,1)*vva(1,2)

	 vva1 = vva1/deta(1)

	
	
end subroutine computejacobians

subroutine computejacobians2(n,vext,vva1,deta)
!> @brief
!> this function computes the volume of a triangle
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:2,1:2)::vva
real,dimension(1:2,1:2),intent(inout)::vva1
real,dimension(1),intent(inout)::deta

vva(1,1) = vext(2,1) - vext(1,1); 	 vva(1,2) = vext(3,1) - vext(1,1)
	  vva(2,1) = vext(2,2) - vext(1,2); 	 vva(2,2) = vext(3,2) - vext(1,2)
      deta(1) = vva(1,1)*vva(2,2) - vva(1,2)*vva(2,1)
	  vva1(1,1) = vext(3,2) - vext(1,2);  vva1(1,2) = -(vext(3,1) - vext(1,1))
	  vva1(2,1) = -(vext(2,2) - vext(1,2));   vva1(2,2) = vext(2,1) - vext(1,1)
      vva1(:,:) = vva1(:,:)/deta(1)




end subroutine computejacobians2



real function hexavolume(n,vext,qpoints,wequa3d)
!> @brief
!> this function computes the volume of a hexahedral
implicit none
#ifdef gpu
!$omp declare target
#endif
!!$omp threadprivate(hexavolume)
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real::s,tx,r,vol
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real,dimension(1:8)::vvxi,vveta,vvzeta,vvnallx,vvnally,vvnallz,vvnxi,vvneta,vvnzeta
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1:3)::vvb,vvc,vve,vvd
real,dimension(1)::deta


integer::kk,ii


	
vvxi(1)=-1.0d0; vveta(1)=-1.0d0; vvzeta(1)=-1.0d0
vvxi(2)=1.0d0; vveta(2)=-1.0d0; vvzeta(2)=-1.0d0
vvxi(3)=1.0d0; vveta(3)=1.0d0; vvzeta(3)=-1.0d0
vvxi(4)=-1.0d0; vveta(4)=1.0d0; vvzeta(4)=-1.0d0
vvxi(5)=-1.0d0; vveta(5)=-1.0d0; vvzeta(5)=1.0d0
vvxi(6)=1.0d0; vveta(6)=-1.0d0; vvzeta(6)=1.0d0
vvxi(7)=1.0d0; vveta(7)=1.0d0; vvzeta(7)=1.0d0
vvxi(8)=-1.0d0; vveta(8)=1.0d0; vvzeta(8)=1.0d0


vvnallx(:)=0.0d0;vvnally(:)=0.0d0;vvnallz(:)=0.0d0


do kk=1,qp_hexa
r=qpoints(1,kk)
s=qpoints(2,kk)
tx=qpoints(3,kk)

do ii=1,8
vvnxi(1)=-(1.0d0/8.0d0)*(1.0-s)*(1.0d0-tx); vvneta(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-tx); vvnzeta(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-s);
vvnxi(2)=(1.0d0/8.0d0)*(1.0-s)*(1.0d0-tx); vvneta(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-tx); vvnzeta(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-s);
vvnxi(3)=(1.0d0/8.0d0)*(1.0+s)*(1.0d0-tx); vvneta(3)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-tx); vvnzeta(3)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+s);
vvnxi(4)=-(1.0d0/8.0d0)*(1.0+s)*(1.0d0-tx); vvneta(4)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-tx); vvnzeta(4)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+s);
vvnxi(5)=-(1.0d0/8.0d0)*(1.0-s)*(1.0d0+tx); vvneta(5)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+tx); vvnzeta(5)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-s);
vvnxi(6)=(1.0d0/8.0d0)*(1.0-s)*(1.0d0+tx); vvneta(6)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+tx); vvnzeta(6)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-s);
vvnxi(7)=(1.0d0/8.0d0)*(1.0+s)*(1.0d0+tx); vvneta(7)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+tx); vvnzeta(7)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+s);
vvnxi(8)=-(1.0d0/8.0d0)*(1.0+s)*(1.0d0+tx); vvneta(8)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+tx); vvnzeta(8)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+s);

vvnallx(ii)=vvnallx(ii)+(vvnxi(ii)*wequa3d(kk))
vvnally(ii)=vvnally(ii)+(vvneta(ii)*wequa3d(kk))
vvnallz(ii)=vvnallz(ii)+(vvnzeta(ii)*wequa3d(kk))
end do
end do



vva=0.0d0
vva1=0.0d0

do ii=1,8

!  r=xi(ii)
!  s=eta(ii)
!  t=zeta(ii)
vvnxi(ii)=vvnallx(ii)
vvneta(ii)=vvnally(ii)
vvnzeta(ii)=vvnallz(ii)





    vva(1,1)=vva(1,1)+vvnxi(ii)*vext(ii,1); vva(1,2)=vva(1,2)+vvnxi(ii)*vext(ii,2); vva(1,3)=vva(1,3)+vvnxi(ii)*vext(ii,3)
    vva(2,1)=vva(2,1)+vvneta(ii)*vext(ii,1); vva(2,2)=vva(2,2)+vvneta(ii)*vext(ii,2); vva(2,3)=vva(2,3)+vvneta(ii)*vext(ii,3)
    vva(3,1)=vva(3,1)+vvnzeta(ii)*vext(ii,1); vva(3,2)=vva(3,2)+vvnzeta(ii)*vext(ii,2); vva(3,3)=vva(3,3)+vvnzeta(ii)*vext(ii,3)
end do

vol=(vva(1,1)*vva(2,2)*vva(3,3))-(vva(1,1)*vva(2,3)*vva(3,2))-(vva(1,2)*vva(2,1)*vva(3,3))+(vva(1,2)*vva(2,3)*vva(3,1))+(vva(1,3)*vva(2,1)*vva(3,2))-(vva(1,3)*vva(2,2)*vva(3,1))

! vol=vol*8.0d0




 vva1(1,1)=(vva(2,2)*vva(3,3))-(vva(2,3)*vva(3,2));vva1(1,2)=((vva(1,3)*vva(3,2))-(vva(3,3)*vva(1,2)));vva1(1,3)=(vva(1,2)*vva(2,3))-(vva(2,2)*vva(1,3));
vva1(2,1)=((vva(2,3)*vva(3,1))-(vva(3,3)*vva(2,1)));vva1(2,2)=((vva(1,1)*vva(3,3))-(vva(3,1)*vva(1,3)));vva1(2,3)=((vva(1,3)*vva(2,1))-(vva(2,3)*vva(1,1)));
vva1(3,1)=(vva(2,1)*vva(3,2))-(vva(3,1)*vva(2,2));vva1(3,2)=((vva(1,2)*vva(3,1))-(vva(3,2)*vva(1,1)));vva1(3,3)=(vva(1,1)*vva(2,2))-(vva(2,1)*vva(1,2));





deta(1)=(vva(1,1)*vva1(1,1))+(vva(1,2)*vva1(2,1))+(vva(1,3)*vva1(3,1))


vol=deta(1)*8.0d0
deta(1)=deta(1)

vva1=vva1/deta(1)




hexavolume=vol




end function hexavolume



real function pyravolume(n,vext,qpoints,wequa3d)
!> @brief
!> this function computes the volume of a pyramid 
implicit none
#ifdef gpu
!$omp declare target
#endif
!!$omp threadprivate(pyravolume)
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real::s,tx,r,vol
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real,dimension(1:8)::vvxi,vveta,vvzeta,vvnallx,vvnally,vvnallz,vvnxi,vvneta,vvnzeta
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1:3)::vvb,vvc,vve,vvd
real,dimension(1)::deta
integer::kk,ii

	
vvxi(1)=-1.0d0; vveta(1)=-1.0d0; vvzeta(1)=-1.0d0
vvxi(2)=1.0d0; vveta(2)=-1.0d0; vvzeta(2)=-1.0d0
vvxi(3)=1.0d0; vveta(3)=1.0d0; vvzeta(3)=-1.0d0
vvxi(4)=-1.0d0; vveta(4)=1.0d0; vvzeta(4)=-1.0d0
vvxi(5)=0.0d0; vveta(5)=0.0d0; vvzeta(5)=1.0d0



vvnallx(:)=0.0d0;vvnally(:)=0.0d0;vvnallz(:)=0.0d0


do kk=1,qp_pyra
r=qpoints(1,kk)
s=qpoints(2,kk)
tx=qpoints(3,kk)

do ii=1,5
vvnxi(1)=-(1.0d0/8.0d0)*(1.0-s)*(1.0d0-tx); vvneta(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-tx); vvnzeta(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-s);
vvnxi(2)=(1.0d0/8.0d0)*(1.0-s)*(1.0d0-tx); vvneta(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-tx); vvnzeta(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-s);
vvnxi(3)=(1.0d0/8.0d0)*(1.0+s)*(1.0d0-tx); vvneta(3)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-tx); vvnzeta(3)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+s);
vvnxi(4)=-(1.0d0/8.0d0)*(1.0+s)*(1.0d0-tx); vvneta(4)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-tx); vvnzeta(4)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+s);
vvnxi(5)=0.0d0; vvneta(5)=0.0d0; vvnzeta(5)=0.5d0;


vvnallx(ii)=vvnallx(ii)+(vvnxi(ii)*wequa3d(kk))
vvnally(ii)=vvnally(ii)+(vvneta(ii)*wequa3d(kk))
vvnallz(ii)=vvnallz(ii)+(vvnzeta(ii)*wequa3d(kk))
end do
end do



vva=0.0d0
vva1=0.0d0

do ii=1,5


vvnxi(ii)=vvnallx(ii)
vvneta(ii)=vvnally(ii)
vvnzeta(ii)=vvnallz(ii)



    vva(1,1)=vva(1,1)+vvnxi(ii)*vext(ii,1); vva(1,2)=vva(1,2)+vvnxi(ii)*vext(ii,2); vva(1,3)=vva(1,3)+vvnxi(ii)*vext(ii,3)
    vva(2,1)=vva(2,1)+vvneta(ii)*vext(ii,1); vva(2,2)=vva(2,2)+vvneta(ii)*vext(ii,2); vva(2,3)=vva(2,3)+vvneta(ii)*vext(ii,3)
    vva(3,1)=vva(3,1)+vvnzeta(ii)*vext(ii,1);vva(3,2)=vva(3,2)+vvnzeta(ii)*vext(ii,2); vva(3,3)=vva(3,3)+vvnzeta(ii)*vext(ii,3)
end do


 vva1(1,1)=(vva(2,2)*vva(3,3))-(vva(2,3)*vva(3,2));vva1(1,2)=((vva(1,3)*vva(3,2))-(vva(3,3)*vva(1,2)));vva1(1,3)=(vva(1,2)*vva(2,3))-(vva(2,2)*vva(1,3));
vva1(2,1)=((vva(2,3)*vva(3,1))-(vva(3,3)*vva(2,1)));vva1(2,2)=((vva(1,1)*vva(3,3))-(vva(3,1)*vva(1,3)));vva1(2,3)=((vva(1,3)*vva(2,1))-(vva(2,3)*vva(1,1)));
vva1(3,1)=(vva(2,1)*vva(3,2))-(vva(3,1)*vva(2,2));vva1(3,2)=((vva(1,2)*vva(3,1))-(vva(3,2)*vva(1,1)));vva1(3,3)=(vva(1,1)*vva(2,2))-(vva(2,1)*vva(1,2));





deta(1)=(vva(1,1)*vva1(1,1))+(vva(1,2)*vva1(2,1))+(vva(1,3)*vva1(3,1))

deta(1)=deta(1)

vva1=vva1/deta(1)

vol=deta(1)

pyravolume=vol










end function pyravolume


real function prismvolume(n,vext,qpoints,wequa3d)
!> @brief
!> this function computes the volume of a prism 
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real::s,tx,r,vol
real,dimension(1:dimensiona,1:numberofpoints),intent(in)::qpoints
real,dimension(1:numberofpoints),intent(in)::wequa3d
real,dimension(1:8)::vvxi,vveta,vvzeta,vvnallx,vvnally,vvnallz,vvnxi,vvneta,vvnzeta
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1:3)::vvb,vvc,vve,vvd
real,dimension(1)::deta
integer::kk,ii




	
vvxi(1)=1.0d0; vveta(1)=0.0d0; vvzeta(1)=-1.0d0
vvxi(2)=0.0d0; vveta(2)=1.0d0; vvzeta(2)=-1.0d0
vvxi(3)=0.0d0; vveta(3)=0.0d0; vvzeta(3)=-1.0d0
vvxi(4)=1.0d0; vveta(4)=0.0d0; vvzeta(4)=0.0d0
vvxi(5)=0.0d0; vveta(5)=1.0d0; vvzeta(5)=0.0d0
vvxi(6)=0.0d0; vveta(6)=0.0d0; vvzeta(6)=0.0d0

 
vvnallx(:)=0.0d0;vvnally(:)=0.0d0;vvnallz(:)=0.0d0
do ii=1,6

do kk=1,qp_prism
r=qpoints(1,kk)
s=qpoints(2,kk)
tx=qpoints(3,kk)



vvnxi(1)=0.5d0*(1.0d0-tx); vvneta(1)=0.0d0; vvnzeta(1)=-0.5d0*r;
vvnxi(2)=0.0d0; vvneta(2)=0.5d0*(1.0d0-tx); vvnzeta(2)=-0.5d0*s;
vvnxi(3)=-0.5d0*(1.0d0-tx); vvneta(3)=-0.5d0*(1.0d0-tx); vvnzeta(3)=-0.5d0*(1.0-r-s);
vvnxi(4)=0.5d0*(1.0d0+tx); vvneta(4)=0.0d0; vvnzeta(4)=0.5d0*r;
vvnxi(5)=0.0d0; vvneta(5)=0.5d0*(1.0d0+tx); vvnzeta(5)=0.5d0*s;
vvnxi(6)=-0.5d0*(1.0d0+tx); vvneta(6)=-0.5d0*(1.0d0+tx); vvnzeta(6)=0.5d0*(1.0-r-s);


vvnallx(ii)=vvnallx(ii)+(vvnxi(ii)*wequa3d(kk))
vvnally(ii)=vvnally(ii)+(vvneta(ii)*wequa3d(kk))
vvnallz(ii)=vvnallz(ii)+(vvnzeta(ii)*wequa3d(kk))
end do
end do



vva=0.0d0
vva1=0.0d0

do ii=1,6

!  r=xi(ii)
!  s=eta(ii)
!  t=zeta(ii)
vvnxi(ii)=vvnallx(ii)
vvneta(ii)=vvnally(ii)
vvnzeta(ii)=vvnallz(ii)





    vva(1,1)=vva(1,1)+vvnxi(ii)*vext(ii,1); vva(1,2)=vva(1,2)+vvnxi(ii)*vext(ii,2); vva(1,3)=vva(1,3)+vvnxi(ii)*vext(ii,3)
    vva(2,1)=vva(2,1)+vvneta(ii)*vext(ii,1); vva(2,2)=vva(2,2)+vvneta(ii)*vext(ii,2); vva(2,3)=vva(2,3)+vvneta(ii)*vext(ii,3)
    vva(3,1)=vva(3,1)+vvnzeta(ii)*vext(ii,1); vva(3,2)=vva(3,2)+vvnzeta(ii)*vext(ii,2); vva(3,3)=vva(3,3)+vvnzeta(ii)*vext(ii,3)
end do


 vva1(1,1)=(vva(2,2)*vva(3,3))-(vva(2,3)*vva(3,2));vva1(1,2)=((vva(1,3)*vva(3,2))-(vva(3,3)*vva(1,2)));vva1(1,3)=(vva(1,2)*vva(2,3))-(vva(2,2)*vva(1,3));
vva1(2,1)=((vva(2,3)*vva(3,1))-(vva(3,3)*vva(2,1)));vva1(2,2)=((vva(1,1)*vva(3,3))-(vva(3,1)*vva(1,3)));vva1(2,3)=((vva(1,3)*vva(2,1))-(vva(2,3)*vva(1,1)));
vva1(3,1)=(vva(2,1)*vva(3,2))-(vva(3,1)*vva(2,2));vva1(3,2)=((vva(1,2)*vva(3,1))-(vva(3,2)*vva(1,1)));vva1(3,3)=(vva(1,1)*vva(2,2))-(vva(2,1)*vva(1,2));





deta(1)=(vva(1,1)*vva1(1,1))+(vva(1,2)*vva1(2,1))+(vva(1,3)*vva1(3,1))

deta(1)=deta(1)

vva1=vva1/deta(1)

vol=deta(1)




prismvolume=vol



end function prismvolume


 function cordinates3(n,nodes_list,n_node)
 !> @brief
!> this function computes the centre of 3d element 
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,n_node
 real,dimension(3)::cordinates3
real::rnode
real,dimension(1:8,1:dimensiona),intent(in)::nodes_list
  rnode=n_node
 cordinates3(1)=sum(nodes_list(1:n_node,1))/rnode
 cordinates3(2)=sum(nodes_list(1:n_node,2))/rnode
 cordinates3(3)=sum(nodes_list(1:n_node,3))/rnode


end function cordinates3


 function distance3(n,vext)
 !> @brief
!> this function computes the distance between two points in 3d
implicit none
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
 real::distance3
real::rnode

  distance3=sqrt(((vext(1,1)-vext(2,1))**2)+((vext(1,2)-vext(2,2))**2)+((vext(1,3)-vext(2,3))**2))
  

end function distance3

 function distance2(n,vext)
 !> @brief
!> this function computes the distance between two points in 2d
implicit none
integer,intent(in)::n
real,dimension(1:8,1:dimensiona),intent(in)::vext
 real::distance2
real::rnode

  distance2=sqrt(((vext(1,1)-vext(2,1))**2)+((vext(1,2)-vext(2,2))**2))
  

end function distance2

 function cordinates2(n,nodes_list,n_node)
  !> @brief
!> this function computes the centre of 2d element 
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,n_node
 real,dimension(2)::cordinates2
real::rnode
real,dimension(1:8,1:dimensiona),intent(in)::nodes_list
  rnode=n_node
 cordinates2(1)=sum(nodes_list(1:n_node,1))/rnode
 cordinates2(2)=sum(nodes_list(1:n_node,2))/rnode
 


end function cordinates2














subroutine decompose3(n,eltype,nodes_list,elem_listd)
 !> @brief
!> this function decomposes element into tetrahedrals (counterclockwise numbering)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,eltype
real,dimension(1:8,1:dimensiona),intent(in)::nodes_list
real,dimension(1:6,1:4,1:dimensiona),intent(inout)::elem_listd


elem_listd(:,:,:)=zero


select case(eltype)
    
    case(1)
      elem_listd(1,1,:)=nodes_list(1,:);elem_listd(1,2,:)=nodes_list(6,:);elem_listd(1,3,:)=nodes_list(8,:);elem_listd(1,4,:)=nodes_list(5,:)
      elem_listd(2,1,:)=nodes_list(1,:);elem_listd(2,2,:)=nodes_list(2,:);elem_listd(2,3,:)=nodes_list(8,:);elem_listd(2,4,:)=nodes_list(6,:)
      elem_listd(3,1,:)=nodes_list(2,:);elem_listd(3,2,:)=nodes_list(7,:);elem_listd(3,3,:)=nodes_list(8,:);elem_listd(3,4,:)=nodes_list(6,:)
      elem_listd(4,1,:)=nodes_list(1,:);elem_listd(4,2,:)=nodes_list(8,:);elem_listd(4,3,:)=nodes_list(3,:);elem_listd(4,4,:)=nodes_list(4,:)
      elem_listd(5,1,:)=nodes_list(1,:);elem_listd(5,2,:)=nodes_list(8,:);elem_listd(5,3,:)=nodes_list(2,:);elem_listd(5,4,:)=nodes_list(3,:)
      elem_listd(6,1,:)=nodes_list(2,:);elem_listd(6,2,:)=nodes_list(8,:);elem_listd(6,3,:)=nodes_list(7,:);elem_listd(6,4,:)=nodes_list(3,:)

    case(2)
      elem_listd(1,1,1:3)=nodes_list(1,1:3);elem_listd(1,2,1:3)=nodes_list(2,1:3);elem_listd(1,3,1:3)=nodes_list(3,1:3);elem_listd(1,4,1:3)=nodes_list(4,1:3)
    case(3)
      elem_listd(1,1,:)=nodes_list(1,:);elem_listd(1,2,:)=nodes_list(2,:);elem_listd(1,3,:)=nodes_list(3,:);elem_listd(1,4,:)=nodes_list(5,:)
      elem_listd(2,1,:)=nodes_list(1,:);elem_listd(2,2,:)=nodes_list(3,:);elem_listd(2,3,:)=nodes_list(4,:);elem_listd(2,4,:)=nodes_list(5,:)
     

    case(4)
    
    
    

      elem_listd(1,1,:)=nodes_list(1,:);elem_listd(1,2,:)=nodes_list(2,:);elem_listd(1,3,:)=nodes_list(3,:);elem_listd(1,4,:)=nodes_list(6,:)
      elem_listd(2,1,:)=nodes_list(1,:);elem_listd(2,2,:)=nodes_list(2,:);elem_listd(2,3,:)=nodes_list(6,:);elem_listd(2,4,:)=nodes_list(5,:)
      elem_listd(3,1,:)=nodes_list(1,:);elem_listd(3,2,:)=nodes_list(5,:);elem_listd(3,3,:)=nodes_list(6,:);elem_listd(3,4,:)=nodes_list(4,:)
     
 


  end select



 

end subroutine decompose3



subroutine decompose2(n,eltype,nodes_list,elem_listd)
 !> @brief
!> this function writes decomposed triangle element nodes into elem_listd from nodes_list (counterclockwise numbering)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,eltype
real,dimension(1:8,1:dimensiona),intent(in)::nodes_list
real,dimension(1:6,1:4,1:dimensiona),intent(inout)::elem_listd




 select case(eltype)
    
    case(5)
      elem_listd(1,1,1:2)=nodes_list(1,1:2);elem_listd(1,2,1:2)=nodes_list(2,1:2);elem_listd(1,3,1:2)=nodes_list(3,1:2)
      elem_listd(2,1,1:2)=nodes_list(1,1:2);elem_listd(2,2,1:2)=nodes_list(3,1:2);elem_listd(2,3,1:2)=nodes_list(4,1:2)
     

    case(6)
      elem_listd(1,1,1:2)=nodes_list(1,1:2);elem_listd(1,2,1:2)=nodes_list(2,1:2);elem_listd(1,3,1:2)=nodes_list(3,1:2)
   
     
 


  end select


end subroutine decompose2


subroutine edge_calculator3d(iconsidered)
 !> @brief
!> this subroutine computes the radius of inscribed sphere or circle
implicit none
integer,intent(in)::iconsidered
integer::i,l,facex,n_node
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real::edgel,dist

	i=iconsidered

    
 	ielem_minedge(i)=(3.0d0*ielem_totvolume(i))/(sum(ielem_surf(1:ielem_ifca(i),i)))
	
	do l=1,ielem_ifca(i)
	facex=l

	 select case (ielem_types_faces(facex,i))
	      case(5)
	      n_node=4
	      case(6)
	      n_node=3
	      end select


				  call coordinates_face_inner(n,i,facex,vext,nodes_list)

 				  vext(2,1:3)=cordinates3(n,nodes_list,n_node)
				  vext(1,1)=ielem_xxc(i);vext(1,2)=ielem_yyc(i); vext(1,3)=ielem_zzc(i)
				  dist=distance3(n,vext)

  				   ielem_minedge(i)=min(dist,ielem_minedge(i))
 	end do
	


end subroutine edge_calculator3d






subroutine edge_calculator2d(iconsidered)
 !> @brief
!> this subroutine computes the radius of inscribed sphere or circle
implicit none
integer,intent(in)::iconsidered
integer::i,l,facex,n_node
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real::edgel,dist



    i=iconsidered
    
	ielem_minedge(i)=(2.0d0*ielem_totvolume(i))/(sum(ielem_surf(1:ielem_ifca(i),i)))
	
	do l=1,ielem_ifca(i)
	facex=l
	n_node=2
				  call coordinates_face_inner2d(n,i,facex,vext,nodes_list)

 				  vext(2,1:2)=cordinates2(n,nodes_list,n_node)
				  vext(1,1)=ielem_xxc(i);vext(1,2)=ielem_yyc(i);
				  dist=distance2(n,vext)


				  ielem_minedge(i)=min(dist,ielem_minedge(i))


	end do
	





end subroutine edge_calculator2d


subroutine geometry_calc
 !> @brief
!> this subroutine computes the volume, surface, centre and min edge for each element
implicit none
integer::kmaxe,i
real::dumv5

kmaxe=xmpielrank(n)


if (dimensiona.eq.3)then

!$omp do 
do i=1,kmaxe

	call volume_calculator3(i)
	call surface_calculator3(i)
	call centre3d(i)
	call edge_calculator3d(i)
end do
!$omp end do

else
!$omp do
do i=1,kmaxe
	call volume_calculator2(i)
	call surface_calculator2(i)
	call centre2d(i)
	call edge_calculator2d(i)
end do
!$omp end do 

end if





!$omp barrier 
!$omp master
dumv5=zero
do i=1,kmaxe
    dumv5=dumv5+ielem_totvolume(i)
end do
 call mpi_allreduce(dumv5,totalvolume,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
!$omp end master
!$omp barrier 



end subroutine geometry_calc




subroutine volume_calculator3(iconsidered)
 !> @brief
!> this subroutine computes the volume of elements
implicit none
integer,intent(in)::iconsidered
integer::i,k,kmaxe,jx,jx2,eltype,elem_dec
real::dumv1,dumv2,dumv3,dumv5
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d



 	i=iconsidered

    
    vext=0.0d0
    nodes_list=0.0d0
    eltype=ielem_ishape(i)
    elem_dec=ielem_vdec(i)
    elem_listd=0.0d0
     ielem_totvolume(i)=0.0d0
      jx=ielem_nonodes(i)
      
	  do k=1,jx
	    jx2=ielem_nodes(k,i)
	    nodes_list(k,:)=dinoder(jx2)%cord(:)
	    vext(k,:)=nodes_list(k,:)
	  end do
	  call decompose3(n,eltype,nodes_list,elem_listd)




    
      select case(ielem_ishape(i))

      case(1)	!hexa
      call quadraturehexa(n,igqrules,vext,qpoints,wequa3d)
      
      dumv1=hexavolume(n,vext,qpoints,wequa3d)
      dumv2=0.0d0
       do k=1,elem_dec
		vext(1:4,1:3)=elem_listd(k,1:4,1:3)

	  dumv2=dumv2+tetravolume(n,vext)
	
		end do
	
! 		if (abs(dumv2-dumv1).le.(0.001d0*abs(dumv2)))then
! 		ielem_totvolume(i)=dumv1
! 		ielem_mode(i)=0
! 		else
! 		ielem_totvolume(i)=dumv2
! 		ielem_mode(i)=1
! 		end if
!
! 		if (dumv1.le.zero)then
! 		ielem_mode(i)=1
! 		ielem_totvolume(i)=dumv2
! 		end if
	
		ielem_mode(i)=1
		ielem_totvolume(i)=dumv2

      case(2)	!tetra
      vext(1:4,1:3)=elem_listd(1,1:4,1:3)
	
		ielem_totvolume(i)=tetravolume(n,vext)
		ielem_mode(i)=1
	
    
      case(3)	!pyramid
      
      call quadraturepyra(n,igqrules,vext,qpoints,wequa3d)
      dumv1=pyravolume(n,vext,qpoints,wequa3d)
    
      dumv2=0.0d0
       do k=1,elem_dec
		vext(1:4,1:3)=elem_listd(k,1:4,1:3)
        dumv3=tetravolume(n,vext)     
	  dumv2=dumv2+tetravolume(n,vext)
		end do
		ielem_totvolume(i)=dumv2
		ielem_mode(i)=1

	
       case(4)	!prism
      call quadratureprism(n,igqrules,vext,qpoints,wequa3d)
      
      dumv1=prismvolume(n,vext,qpoints,wequa3d)

      dumv2=0.0d0
       do k=1,elem_dec
		vext(1:4,1:3)=elem_listd(k,1:4,1:3)
				dumv3=tetravolume(n,vext)

		dumv2=dumv2+tetravolume(n,vext)
		
		end do

! 	if (abs(dumv2-dumv1).le.(0.001d0*abs(dumv2)))then
! 	ielem_totvolume(i)=dumv1
! 	ielem_mode(i)=0
! 	else
! 	ielem_totvolume(i)=dumv2
! 	ielem_mode(i)=1
! 	end if
! 	if (dumv1.le.zero)then
! 	ielem_mode(i)=1
! 	end if
	
 	
 	ielem_mode(i)=1
        ielem_totvolume(i)=dumv2

      end select
   
   
    


end subroutine volume_calculator3



subroutine volume_calculator2(iconsidered)
 !> @brief
!> this subroutine computes the volume of elements in 2d
implicit none
integer,intent(in)::iconsidered
!$ integer::omp_in_parallel,omp_get_thread_num
integer::i,k,kmaxe,jx,jx2,eltype,elem_dec
real::dumv1,dumv2,dumv3,dumv5
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
real,dimension(1:dimensiona,1:numberofpoints)::qpoints
real,dimension(1:numberofpoints)::wequa3d

i=iconsidered


    eltype=ielem_ishape(i)
    elem_dec=ielem_vdec(i)
     ielem_totvolume(i)=0.0d0
	  do k=1,ielem_nonodes(i)
	    nodes_list(k,1:2)=dinoder(ielem_nodes(k,i))%cord(1:2)
	    vext(k,1:2)=nodes_list(k,1:2)
	  end do
	 call decompose2(n,eltype,nodes_list,elem_listd)

	
	    select case(ielem_ishape(i))

      case(5)

      call quadraturequad(n,igqrules,vext,qpoints,wequa3d)
      dumv1=quadvolume(n,vext,qpoints,wequa3d)

      

      
      dumv2=0.0d0
      

      
      
       do k=1,elem_dec
	vext(1:3,1:2)=elem_listd(k,1:3,1:2)
	  dumv2=dumv2+trianglevolume(n,vext)
    
	end do


! 	if (abs(dumv2-dumv1).le.(0.001*dumv2))then
! 	ielem_totvolume(i)=dumv1
! 	ielem_mode(i)=0
! 	else
! 	ielem_totvolume(i)=dumv2
! 	ielem_mode(i)=1
! 	end if



 	ielem_mode(i)=1
 	ielem_totvolume(i)=dumv2
     
      case(6)

      dumv1=trianglevolume(n,vext)

      

      dumv2=0.0d0
       do k=1,elem_dec
	vext(1:3,1:2)=elem_listd(k,1:3,1:2)
	  dumv2=dumv2+trianglevolume(n,vext)
    
	end do
! 	if (abs(dumv2-dumv1).le.(0.001d0*dumv2))then
! 	ielem_totvolume(i)=dumv1
! 	ielem_mode(i)=0
! 	else
! 	ielem_totvolume(i)=dumv2
! 	ielem_mode(i)=1
! 	end if
        ielem_mode(i)=1
 	ielem_totvolume(i)=dumv2
     
    end select
   
    


 
    


end subroutine volume_calculator2




subroutine surface_calculator3(iconsidered)
 !> @brief
!> this subroutine computes the surface area of elements in 3d
implicit none
integer,intent(in)::iconsidered
integer::i,k,kmaxe,jx,jx2,eltype,elem_dec,nnd,j
real::dumv1,dumv2,dumv3,dumv5
real,dimension(1:8,1:dimensiona)::vext


	i=iconsidered


    
    do j=1,ielem_ifca(i)
				select case(ielem_types_faces(j,i))
				case (5)
					 
					  nnd=4
				      do k=1,nnd
					vext(k,1:dims)=dinoder(ielem_nodes_faces(j,k,i))%cord(1:dims)
				      end do
					  
					  
					  ielem_surf(j,i)=quadarea(n,vext)
					  
				  
				case(6)
					
					nnd=3
					do k=1,nnd
					  vext(k,1:dims)=dinoder(ielem_nodes_faces(j,k,i))%cord(1:dims)
					end do
					    
					
 					    ielem_surf(j,i)=trianglearea(n,vext)
 					    
 					    
					  
				end select
    
    
    
    end do
    
    
    
    do j=1,ielem_ifca(i)
				select case(ielem_types_faces(j,i))
                        case (5)
                                dumv2=zero
				
					 
					
				     
					vext(1,1:dims)=dinoder(ielem_nodes_faces(j,1,i))%cord(1:dims)
					vext(2,1:dims)=dinoder(ielem_nodes_faces(j,2,i))%cord(1:dims)
					vext(3,1:dims)=dinoder(ielem_nodes_faces(j,3,i))%cord(1:dims)
				     
					  
					  
					  dumv2=dumv2+trianglearea(n,vext)
					  
					  vext(1,1:dims)=dinoder(ielem_nodes_faces(j,1,i))%cord(1:dims)
					vext(2,1:dims)=dinoder(ielem_nodes_faces(j,3,i))%cord(1:dims)
					vext(3,1:dims)=dinoder(ielem_nodes_faces(j,4,i))%cord(1:dims)
					  dumv2=dumv2+trianglearea(n,vext)
					  
					  if (abs((dumv2-ielem_surf(j,i))/ielem_surf(j,i))*100.0d0.gt.10.0d0)then
! 					  
					  
					  
					  ielem_surf(j,i)=dumv2
						end if
                                          
                                          

    
    end select
    end do
    
    

    


end subroutine surface_calculator3



subroutine surface_calculator2(iconsidered)
 !> @brief
!> this subroutine computes the length of edges of elements in 2d
implicit none
integer,intent(in)::iconsidered
integer::i,k,jx,jx2,eltype,elem_dec,nnd,j
real::dumv1,dumv2,dumv3,dumv5,dumr
real,dimension(1:8,1:dimensiona)::vext



i=iconsidered
    do j=1,ielem_ifca(i)
        nnd=2
        
        do k=1,nnd
            vext(k,1:dims)=dinoder(ielem_nodes_faces(j,k,i))%cord(1:dims)
        end do
        
        ielem_surf(j,i)=linearea(n,vext)
    end do


end subroutine surface_calculator2











subroutine compute_centre3df(n,iconsidered,facex,n_node,cords)
 !> @brief
!> this subroutine retrieve the nodes of faces of elements in 3d
implicit none
integer,intent(in)::n,iconsidered,facex
integer,intent(inout)::n_node
real,dimension(1:dimensiona),intent(inout)::cords
real,dimension(1:8,1:dimensiona)::nodes_list
integer::k,i

i=iconsidered
  
    do k=1,n_node

      nodes_list(k,1:3)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:3)
    end do
   

    cords=cordinates3(n,nodes_list,n_node)
   


end subroutine


subroutine coordinates_face_inner(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieve the nodes of interior faces of elements in 3d
implicit none
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k
i=iconsidered


	      select case (ielem_types_faces(facex,i))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      
	      do k=1,nnd
		  nodes_list(k,1:3)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:3)
		  vext(k,1:3)=nodes_list(k,1:3)
	      end do
	      
	      

	      
	     
end subroutine coordinates_face_inner


subroutine coordinates_face_inner2d(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieves the nodes of edges of elements in 2d
implicit none
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k
i=iconsidered


	      nnd=2
	      
	      
	      do k=1,nnd
		  nodes_list(k,1:2)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:2)
		  vext(k,1:2)=nodes_list(k,1:2)
	      end do
	      
	      
	      
	     
end subroutine coordinates_face_inner2d


subroutine coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieve the nodes of interior faces of elements in 3d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k


i=iconsidered
	      select case (ielem_types_faces(facex,i))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      
	      do k=1,nnd
		  nodes_list(k,1:3)=inoder4_cord(1:3,ielem_nodes_faces(facex,k,i))
		  vext(k,1:3)=nodes_list(k,1:3)
	      end do
	      
	      
	     
end subroutine coordinates_face_innerx


subroutine coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieves the nodes of edges of elements in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k

i=iconsidered


	      nnd=2
	      
	      
	      do k=1,nnd
		  nodes_list(k,1:2)=inoder4_cord(1:2,ielem_nodes_faces(facex,k,i))
		  vext(k,1:2)=nodes_list(k,1:2)
	      end do
	      
	      

	      
	     
end subroutine coordinates_face_inner2dx


subroutine coordinates_face_period1(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieve the nodes of periodic faces of elements in 3d
implicit none
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k
real::tempxx
i=iconsidered


	      select case (ielem_types_faces(facex,i))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      vext(1,1)=ielem_xxc(i)
	      vext(1,2)=ielem_yyc(i)
	      vext(1,3)=ielem_zzc(i)
	      
			do k=1,nnd
			  nodes_list(k,1:3)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:dims)
			end do
			do k=1,nnd
			if(per_rot.eq.0)then
			if(abs(nodes_list(k,1)-vext(1,1)).gt.xper*oo2)then
			nodes_list(k,1)=nodes_list(k,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
			end if
			if(abs(nodes_list(k,2)-vext(1,2)).gt.yper*oo2)then
			nodes_list(k,2)=nodes_list(k,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
			end if
			if(abs(nodes_list(k,3)-vext(1,3)).gt.zper*oo2)then
			nodes_list(k,3)=nodes_list(k,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
			end if
			else
                if (ielem_reorient(facex,i).eq.1) then
                    if (ibound_icode(ielem_ibounds(facex,i)).eq.5) then
                        tempxx=nodes_list(k,1)
                        nodes_list(k,1)=tempxx*cos(-angle_per)-sin(-angle_per)*nodes_list(k,2)
                        nodes_list(k,2)=tempxx*sin(-angle_per)+cos(-angle_per)*nodes_list(k,2)
                    else
                        tempxx=nodes_list(k,1)
                        nodes_list(k,1)=tempxx*cos(angle_per)-sin(angle_per)*nodes_list(k,2)
                        nodes_list(k,2)=tempxx*sin(angle_per)+cos(angle_per)*nodes_list(k,2)
                    end if
                end if
			end if
			end do
	      
	      
	      do k=1,nnd
		  vext(k,1:3)=nodes_list(k,1:3)
	      end do
	      

end subroutine coordinates_face_period1


subroutine coordinates_face_period2d1(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieve the nodes of periodic edges of elements in 2d
implicit none
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k
i=iconsidered


	     nnd=2
	      
	      vext(1,1)=ielem_xxc(i)
	      vext(1,2)=ielem_yyc(i)
	      
	      
			do k=1,nnd
			  nodes_list(k,1:2)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:dims)
			end do
			do k=1,nnd
			if(abs(nodes_list(k,1)-vext(1,1)).gt.xper*oo2)then
			nodes_list(k,1)=nodes_list(k,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
			end if
			if(abs(nodes_list(k,2)-vext(1,2)).gt.yper*oo2)then
			nodes_list(k,2)=nodes_list(k,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
			end if
			end do
	      
	      
	      do k=1,nnd
		  vext(k,1:2)=nodes_list(k,1:2)
	      end do
	      

end subroutine coordinates_face_period2d1



subroutine coordinates_face_period(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieve the nodes of periodic faces of elements in 3d
implicit none
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k
real::tempxx
i=iconsidered


	      select case (ielem_types_faces(facex,i))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      vext(1,1)=ielem_xxc(i)
	      vext(1,2)=ielem_yyc(i)
	      vext(1,3)=ielem_zzc(i)
	      
			do k=1,nnd
			  nodes_list(k,1:3)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:dims)
			end do
            if(per_rot.eq.0)then
			do k=1,nnd
			if(abs(nodes_list(k,1)-vext(1,1)).gt.xper*oo2)then
			nodes_list(k,1)=nodes_list(k,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
			end if
			if(abs(nodes_list(k,2)-vext(1,2)).gt.yper*oo2)then
			nodes_list(k,2)=nodes_list(k,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
			end if
			if(abs(nodes_list(k,3)-vext(1,3)).gt.zper*oo2)then
			nodes_list(k,3)=nodes_list(k,3)+(zper*sign(1.0,vext(1,3)-zper*oo2))
			end if
			end do
			else
                if (ielem_reorient(facex,i).eq.1) then
                do k=1,nnd
                    if (ibound_icode(ielem_ibounds(facex,i)).eq.5) then
                        tempxx=nodes_list(k,1)
                        nodes_list(k,1)=tempxx*cos(-angle_per)-sin(-angle_per)*nodes_list(k,2)
                        nodes_list(k,2)=tempxx*sin(-angle_per)+cos(-angle_per)*nodes_list(k,2)
                    else
                        tempxx=nodes_list(k,1)
                        nodes_list(k,1)=tempxx*cos(angle_per)-sin(angle_per)*nodes_list(k,2)
                        nodes_list(k,2)=tempxx*sin(angle_per)+cos(angle_per)*nodes_list(k,2)
                    end if
                end do
                end if
			end if
	      
	      
	      do k=1,nnd
		  vext(k,1:3)=nodes_list(k,1:3)
	      end do
	      
end subroutine coordinates_face_period


subroutine coordinates_face_period2d(n,iconsidered,facex,vext,nodes_list)
 !> @brief
!> this subroutine retrieve the nodes of periodic edges of elements in 2d
implicit none
integer,intent(in)::n,iconsidered,facex
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1:8,1:dimensiona),intent(inout)::nodes_list
integer::nnd
integer::i,k
i=iconsidered


	     nnd=2
	      
	      vext(1,1)=ielem_xxc(i)
	      vext(1,2)=ielem_yyc(i)
	      
	      
			do k=1,nnd
			  nodes_list(k,1:2)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:2)
			end do
			do k=1,nnd
			if(abs(nodes_list(k,1)-vext(1,1)).gt.xper*oo2)then
			nodes_list(k,1)=nodes_list(k,1)+(xper*sign(1.0,vext(1,1)-xper*oo2))
			end if
			if(abs(nodes_list(k,2)-vext(1,2)).gt.yper*oo2)then
			nodes_list(k,2)=nodes_list(k,2)+(yper*sign(1.0,vext(1,2)-yper*oo2))
			end if
			end do
	      
	      
	      do k=1,nnd
		  vext(k,1:2)=nodes_list(k,1:2)
	      end do
	      

end subroutine coordinates_face_period2d



subroutine compute_centre2df(n,iconsidered,facex,n_node,cords)
 !> @brief
!> this subroutine retrieves the nodes of the vertices of edges of 2d elements
implicit none
integer,intent(in)::n,iconsidered,facex
integer,intent(inout)::n_node
real,dimension(1:dimensiona),intent(inout)::cords
real,dimension(1:8,1:dimensiona)::nodes_list
integer::k,i
i=iconsidered

    
    do k=1,n_node
      nodes_list(k,1:2)=dinoder(ielem_nodes_faces(facex,k,i))%cord(1:2)
    end do
    cords=cordinates2(n,nodes_list,n_node)
   


end subroutine



subroutine compute_centre3d(iconsidered,cords)
 !> @brief
!> this subroutine computes the cell centre of elements in 3d
implicit none
integer,intent(in)::iconsidered
integer::k,i,j,n_node
real,dimension(1:dimensiona),intent(inout)::cords
real,dimension(1:8,1:dimensiona)::nodes_list


i=iconsidered

    n_node=ielem_nonodes(i)
    do k=1,ielem_nonodes(i)
      nodes_list(k,1:3)=dinoder(ielem_nodes(k,i))%cord(1:3)
    end do
    cords=cordinates3(n,nodes_list,n_node)




end subroutine



subroutine compute_centre2d(iconsidered,cords)
 !> @brief
!> this subroutine retrieves the nodes of the vertices of 2d elements
implicit none
integer,intent(in)::iconsidered
integer::k,i,j,n_node
real,dimension(1:dimensiona),intent(inout)::cords
real,dimension(1:8,1:dimensiona)::nodes_list


i=iconsidered

    n_node=ielem_nonodes(i)
    do k=1,ielem_nonodes(i)
      nodes_list(k,1:2)=dinoder(ielem_nodes(k,i))%cord(1:2)
    end do
    cords=cordinates2(n,nodes_list,n_node)
   


end subroutine


subroutine centre3d(iconsidered)
 !> @brief
!> this subroutine computes the cell centres
implicit none
integer,intent(in)::iconsidered
real,dimension(1:dimensiona)::cords
integer::i
i=iconsidered



    call compute_centre3d(i,cords)
    ielem_xxc(i)=cords(1);ielem_yyc(i)=cords(2);ielem_zzc(i)=cords(3);





end subroutine centre3d


subroutine centre2d(iconsidered)
 !> @brief
!> this subroutine computes the cell centres
implicit none
integer,intent(in)::iconsidered
real,dimension(1:dimensiona)::cords
integer::i
i=iconsidered


    call compute_centre2d(i,cords)
    ielem_xxc(i)=cords(1);ielem_yyc(i)=cords(2)




end subroutine centre2d

subroutine quadraturetriang(n,igqrules,vext,qpoints2d,wequa2d)
 !> @brief
!> this subroutine computes the quadrature points and weights for triangle in 3d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints2),intent(inout)::qpoints2d
real,dimension(1:numberofpoints2),intent(inout)::wequa2d
real,dimension(1:dimensiona,1:dimensiona)::vva,vva1
real,dimension(1)::deta
real,dimension(1:4)::vvnxi
real,dimension(1:alls)::vvwg
real,dimension(1:alls)::vvr1,vvr2,vvr3
integer::kk

wequa2d=0.0d0
qpoints2d=0.0d0


select case(igqrules)

case(1)
		
	    vvwg(1) = 1.0d0
	    vvr1(1)=1.0d0/3.0d0;	vvr2(1)=1.0d0/3.0d0;	vvr3(1)=1.0d0/3.0d0
case(2)
		vvwg(1)=0.33333333333333333333
  		vvwg(2)=0.33333333333333333333
  		vvwg(3)=0.33333333333333333333

		vvr1(1)=0.666666666666667 ;vvr2(1)=0.166666666666667 ;vvr3(1)=0.166666666666667 
		vvr1(2)=0.166666666666667 ;vvr2(2)=0.666666666666667 ;vvr3(2)=0.166666666666667 
		vvr1(3)=0.166666666666667 ;vvr2(3)=0.166666666666667 ;vvr3(3)=0.666666666666667 


case(3)
		vvr1(1)=0.816847572980440 ;vvr2(1)=0.091576213509780 ;vvr3(1)=0.091576213509780 ;vvwg(1)=0.109951743655333
		vvr1(2)=0.091576213509780 ;vvr2(2)=0.816847572980440 ;vvr3(2)=0.091576213509780 ;vvwg(2)=0.109951743655333
		vvr1(3)=0.091576213509780 ;vvr2(3)=0.091576213509780 ;vvr3(3)=0.816847572980440 ;vvwg(3)=0.109951743655333
		vvr1(4)=0.445948490915964 ;vvr2(4)=0.445948490915964 ;vvr3(4)=0.108103018168071 ;vvwg(4)=0.223381589678000
		vvr1(5)=0.445948490915964 ;vvr2(5)=0.108103018168071 ;vvr3(5)=0.445948490915964 ;vvwg(5)=0.223381589678000
		vvr1(6)=0.108103018168071 ;vvr2(6)=0.445948490915964 ;vvr3(6)=0.445948490915964 ;vvwg(6)=0.223381589678000

case(4)
		
vvr1(1)=0.888871894660413 ;vvr2(1)=0.055564052669793 ;vvr3(1)=0.055564052669793 ;vvwg(1)=0.041955512996649
vvr1(2)=0.055564052669793 ;vvr2(2)=0.888871894660413 ;vvr3(2)=0.055564052669793 ;vvwg(2)=0.041955512996649
vvr1(3)=0.055564052669793 ;vvr2(3)=0.055564052669793 ;vvr3(3)=0.888871894660413 ;vvwg(3)=0.041955512996649
vvr1(4)=0.295533711735893 ;vvr2(4)=0.634210747745723 ;vvr3(4)=0.070255540518384 ;vvwg(4)=0.112098412070887
vvr1(5)=0.295533711735893 ;vvr2(5)=0.070255540518384 ;vvr3(5)=0.634210747745723 ;vvwg(5)=0.112098412070887
vvr1(6)=0.070255540518384 ;vvr2(6)=0.295533711735893 ;vvr3(6)=0.634210747745723 ;vvwg(6)=0.112098412070887
vvr1(7)=0.634210747745723 ;vvr2(7)=0.295533711735893 ;vvr3(7)=0.070255540518384 ;vvwg(7)=0.112098412070887
vvr1(8)=0.634210747745723 ;vvr2(8)=0.070255540518384 ;vvr3(8)=0.295533711735893 ;vvwg(8)=0.112098412070887
vvr1(9)=0.070255540518384 ;vvr2(9)=0.634210747745723 ;vvr3(9)=0.295533711735893 ;vvwg(9)=0.112098412070887
vvr1(10)=0.333333333333333 ;vvr2(10)=0.333333333333333 ;vvr3(10)=0.333333333333333 ;vvwg(10)=0.201542988584730


case(5)

vvr1(1)=0.928258244608533;vvr2(1)= 0.035870877695734 ;vvr3(1)=0.035870877695734 ;vvwg(1)=0.017915455012303
vvr1(2)=0.035870877695734;vvr2(2)= 0.928258244608533 ;vvr3(2)=0.035870877695734 ;vvwg(2)=0.017915455012303
vvr1(3)=0.035870877695734;vvr2(3)= 0.035870877695734 ;vvr3(3)=0.928258244608533 ;vvwg(3)=0.017915455012303
vvr1(4)=0.516541208464066;vvr2(4)= 0.241729395767967 ;vvr3(4)=0.241729395767967 ;vvwg(4)=0.127712195881265
vvr1(5)=0.241729395767967;vvr2(5)= 0.516541208464066 ;vvr3(5)=0.241729395767967 ;vvwg(5)=0.127712195881265
vvr1(6)=0.241729395767967;vvr2(6)= 0.241729395767967 ;vvr3(6)=0.516541208464066 ;vvwg(6)=0.127712195881265
vvr1(7)=0.474308787777079;vvr2(7)= 0.474308787777079 ;vvr3(7)=0.051382424445843 ;vvwg(7)=0.076206062385535
vvr1(8)=0.474308787777079;vvr2(8)= 0.051382424445843 ;vvr3(8)=0.474308787777079 ;vvwg(8)=0.076206062385535
vvr1(9)=0.051382424445843;vvr2(9)= 0.474308787777079 ;vvr3(9)=0.474308787777079 ;vvwg(9)=0.076206062385535
vvr1(10)=0.201503881881800;vvr2(10)= 0.751183631106484 ;vvr3(10)=0.047312487011716 ;vvwg(10)=0.055749810027115
vvr1(11)=0.201503881881800;vvr2(11)= 0.047312487011716 ;vvr3(11)=0.751183631106484 ;vvwg(11)=0.055749810027115
vvr1(12)=0.047312487011716;vvr2(12)= 0.201503881881800 ;vvr3(12)=0.751183631106484 ;vvwg(12)=0.055749810027115
vvr1(13)=0.751183631106484;vvr2(13)= 0.201503881881800 ;vvr3(13)=0.047312487011716 ;vvwg(13)=0.055749810027115
vvr1(14)=0.751183631106484;vvr2(14)= 0.047312487011716 ;vvr3(14)=0.201503881881800 ;vvwg(14)=0.055749810027115
vvr1(15)=0.047312487011716;vvr2(15)= 0.751183631106484 ;vvr3(15)=0.201503881881800 ;vvwg(15)=0.055749810027115








case(6)

vvr1(1)=0.943774095634672   ;vvr2(1)=0.028112952182664  ;vvr3(1)=0.028112952182664 ;vvwg(1)=0.010359374696538
vvr1(2)=0.028112952182664   ;vvr2(2)=0.943774095634672  ;vvr3(2)=0.028112952182664 ;vvwg(2)=0.010359374696538
vvr1(3)=0.028112952182664   ;vvr2(3)=0.028112952182664  ;vvr3(3)=0.943774095634672 ;vvwg(3)=0.010359374696538
vvr1(4)=0.645721803061365   ;vvr2(4)=0.177139098469317  ;vvr3(4)=0.177139098469317 ;vvwg(4)=0.075394884326738
vvr1(5)=0.177139098469317   ;vvr2(5)=0.645721803061365  ;vvr3(5)=0.177139098469317 ;vvwg(5)=0.075394884326738
vvr1(6)=0.177139098469317   ;vvr2(6)=0.177139098469317  ;vvr3(6)=0.645721803061365 ;vvwg(6)=0.075394884326738
vvr1(7)=0.405508595867433   ;vvr2(7)=0.405508595867433  ;vvr3(7)=0.188982808265134 ;vvwg(7)=0.097547802373242
vvr1(8)=0.405508595867433   ;vvr2(8)=0.188982808265134  ;vvr3(8)=0.405508595867433 ;vvwg(8)=0.097547802373242
vvr1(9)=0.188982808265134   ;vvr2(9)=0.405508595867433  ;vvr3(9)=0.405508595867433 ;vvwg(9)=0.097547802373242
vvr1(10)=0.148565812270887 ;vvr2(10)=0.817900980028499 ;vvr3(10)=0.033533207700614 ;vvwg(10)=0.028969269372473
vvr1(11)=0.148565812270887 ;vvr2(11)=0.033533207700614 ;vvr3(11)=0.817900980028499 ;vvwg(11)=0.028969269372473
vvr1(12)=0.033533207700614 ;vvr2(12)=0.148565812270887 ;vvr3(12)=0.817900980028499 ;vvwg(12)=0.028969269372473
vvr1(13)=0.817900980028499 ;vvr2(13)=0.148565812270887 ;vvr3(13)=0.033533207700614 ;vvwg(13)=0.028969269372473
vvr1(14)=0.817900980028499 ;vvr2(14)=0.033533207700614 ;vvr3(14)=0.148565812270887 ;vvwg(14)=0.028969269372473
vvr1(15)=0.033533207700614 ;vvr2(15)=0.817900980028499 ;vvr3(15)=0.148565812270887 ;vvwg(15)=0.028969269372473
vvr1(16)=0.357196298615681 ;vvr2(16)=0.604978911775132 ;vvr3(16)=0.037824789609186 ;vvwg(16)=0.046046366595935
vvr1(17)=0.357196298615681 ;vvr2(17)=0.037824789609186 ;vvr3(17)=0.604978911775132 ;vvwg(17)=0.046046366595935
vvr1(18)=0.037824789609186 ;vvr2(18)=0.357196298615681 ;vvr3(18)=0.604978911775132 ;vvwg(18)=0.046046366595935
vvr1(19)=0.604978911775132 ;vvr2(19)=0.357196298615681 ;vvr3(19)=0.037824789609186 ;vvwg(19)=0.046046366595935
vvr1(20)=0.604978911775132 ;vvr2(20)=0.037824789609186 ;vvr3(20)=0.357196298615681 ;vvwg(20)=0.046046366595935
vvr1(21)=0.037824789609186 ;vvr2(21)=0.604978911775132 ;vvr3(21)=0.357196298615681 ;vvwg(21)=0.046046366595935




case(7,8,9)

vvr1(1)=0.957657154441070
vvr1(2)=0.021171422779465
vvr1(3)=0.021171422779465
vvr1(4)=0.798831205208225
vvr1(5)=0.100584397395888
vvr1(6)=0.100584397395888
vvr1(7)=0.457923384576135
vvr1(8)=0.271038307711932
vvr1(9)=0.271038307711932
vvr1(10)=0.440191258403832
vvr1(11)=0.440191258403832
vvr1(12)=0.119617483192335
vvr1(13)=0.101763679498021
vvr1(14)=0.101763679498021
vvr1(15)=0.018256679074748
vvr1(16)=0.879979641427232
vvr1(17)=0.879979641427232
vvr1(18)=0.018256679074748
vvr1(19)=0.394033271669987
vvr1(20)=0.394033271669987
vvr1(21)=0.023404705466341
vvr1(22)=0.582562022863673
vvr1(23)=0.582562022863673
vvr1(24)=0.023404705466341
vvr1(25)=0.226245530909229
vvr1(26)=0.226245530909229
vvr1(27)=0.022223854547989
vvr1(28)=0.751530614542782
vvr1(29)=0.751530614542782
vvr1(30)=0.022223854547989
vvr1(31)=0.635737183263105
vvr1(32)=0.635737183263105
vvr1(33)=0.115183589115563
vvr1(34)=0.249079227621332
vvr1(35)=0.249079227621332
vvr1(36)=0.115183589115563



vvr2(1)=0.021171422779465
vvr2(2)=0.957657154441070
vvr2(3)=0.021171422779465
vvr2(4)=0.100584397395888
vvr2(5)=0.798831205208225
vvr2(6)=0.100584397395888
vvr2(7)=0.271038307711932
vvr2(8)=0.457923384576135
vvr2(9)=0.271038307711932
vvr2(10)=0.440191258403832
vvr2(11)=0.119617483192335
vvr2(12)=0.440191258403832
vvr2(13)=0.879979641427232
vvr2(14)=0.018256679074748
vvr2(15)=0.101763679498021
vvr2(16)=0.101763679498021
vvr2(17)=0.018256679074748
vvr2(18)=0.879979641427232
vvr2(19)=0.582562022863673
vvr2(20)=0.023404705466341
vvr2(21)=0.394033271669987
vvr2(22)=0.394033271669987
vvr2(23)=0.023404705466341
vvr2(24)=0.582562022863673
vvr2(25)=0.751530614542782
vvr2(26)=0.022223854547989
vvr2(27)=0.226245530909229
vvr2(28)=0.226245530909229
vvr2(29)=0.022223854547989
vvr2(30)=0.751530614542782
vvr2(31)=0.249079227621332
vvr2(32)=0.115183589115563
vvr2(33)=0.635737183263105
vvr2(34)=0.635737183263105
vvr2(35)=0.115183589115563
vvr2(36)=0.249079227621332




vvr3(1)=0.021171422779465
vvr3(2)=0.021171422779465
vvr3(3)=0.957657154441070
vvr3(4)=0.100584397395888
vvr3(5)=0.100584397395888
vvr3(6)=0.798831205208225
vvr3(7)=0.271038307711932
vvr3(8)=0.271038307711932
vvr3(9)=0.457923384576135
vvr3(10)=0.119617483192335
vvr3(11)=0.440191258403832
vvr3(12)=0.440191258403832
vvr3(13)=0.018256679074748
vvr3(14)=0.879979641427232
vvr3(15)=0.879979641427232
vvr3(16)=0.018256679074748
vvr3(17)=0.101763679498021
vvr3(18)=0.101763679498021
vvr3(19)=0.023404705466341
vvr3(20)=0.582562022863673
vvr3(21)=0.582562022863673
vvr3(22)=0.023404705466341
vvr3(23)=0.394033271669987
vvr3(24)=0.394033271669987
vvr3(25)=0.022223854547989
vvr3(26)=0.751530614542782
vvr3(27)=0.751530614542782
vvr3(28)=0.022223854547989
vvr3(29)=0.226245530909229
vvr3(30)=0.226245530909229
vvr3(31)=0.115183589115563
vvr3(32)=0.249079227621332
vvr3(33)=0.249079227621332
vvr3(34)=0.115183589115563
vvr3(35)=0.635737183263105
vvr3(36)=0.635737183263105



vvwg(1)=0.005639123786910
vvwg(2)=0.005639123786910
vvwg(3)=0.005639123786910
vvwg(4)=0.027148968192278
vvwg(5)=0.027148968192278
vvwg(6)=0.027148968192278
vvwg(7)=0.063100912533359
vvwg(8)=0.063100912533359
vvwg(9)=0.063100912533359
vvwg(10)=0.051752795679899
vvwg(11)=0.051752795679899
vvwg(12)=0.051752795679899
vvwg(13)=0.009866753574646
vvwg(14)=0.009866753574646
vvwg(15)=0.009866753574646
vvwg(16)=0.009866753574646
vvwg(17)=0.009866753574646
vvwg(18)=0.009866753574646
vvwg(19)=0.022008204800147
vvwg(20)=0.022008204800147
vvwg(21)=0.022008204800147
vvwg(22)=0.022008204800147
vvwg(23)=0.022008204800147
vvwg(24)=0.022008204800147
vvwg(25)=0.016644570076736
vvwg(26)=0.016644570076736
vvwg(27)=0.016644570076736
vvwg(28)=0.016644570076736
vvwg(29)=0.016644570076736
vvwg(30)=0.016644570076736
vvwg(31)=0.044326238118914
vvwg(32)=0.044326238118914
vvwg(33)=0.044326238118914
vvwg(34)=0.044326238118914
vvwg(35)=0.044326238118914
vvwg(36)=0.044326238118914















		
end select



		do kk=1,qp_triangle
			wequa2d(kk)=vvwg(kk)
			qpoints2d(:,kk)=(vvr1(kk)*vext(1,:))+(vvr2(kk)*vext(2,:))+(vvr3(kk)*vext(3,:))
			
		end do




end subroutine quadraturetriang

subroutine quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)
!> @brief
!> this subroutine computes the quadrature points and weights for triangle in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints),intent(inout)::qpoints
real,dimension(1:numberofpoints),intent(inout)::wequa3d
real,dimension(1:dimensiona,1:dimensiona)::vva,vva1
real,dimension(1)::deta
real,dimension(1:4)::vvnxi
real,dimension(1:alls)::vvwg,vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
integer::kk

wequa3d=0.0d0
qpoints=0.0d0

select case(igqrules)


case(1)
		
	    vvwg(1) = 1.0d0
	    vvr1(1)=1.0d0/3.0d0;	vvr2(1)=1.0d0/3.0d0;	vvr3(1)=1.0d0/3.0d0

case(2)
		vvwg(1)=0.33333333333333333333
  		vvwg(2)=0.33333333333333333333
  		vvwg(3)=0.33333333333333333333

		vvr1(1)=0.666666666666667 ;vvr2(1)=0.166666666666667 ;vvr3(1)=0.166666666666667 
		vvr1(2)=0.166666666666667 ;vvr2(2)=0.666666666666667 ;vvr3(2)=0.166666666666667 
		vvr1(3)=0.166666666666667 ;vvr2(3)=0.166666666666667 ;vvr3(3)=0.666666666666667 

case(3)
		vvr1(1)=0.816847572980440 ;vvr2(1)=0.091576213509780 ;vvr3(1)=0.091576213509780 ;vvwg(1)=0.109951743655333
		vvr1(2)=0.091576213509780 ;vvr2(2)=0.816847572980440 ;vvr3(2)=0.091576213509780 ;vvwg(2)=0.109951743655333
		vvr1(3)=0.091576213509780 ;vvr2(3)=0.091576213509780 ;vvr3(3)=0.816847572980440 ;vvwg(3)=0.109951743655333
		vvr1(4)=0.445948490915964 ;vvr2(4)=0.445948490915964 ;vvr3(4)=0.108103018168071 ;vvwg(4)=0.223381589678000
		vvr1(5)=0.445948490915964 ;vvr2(5)=0.108103018168071 ;vvr3(5)=0.445948490915964 ;vvwg(5)=0.223381589678000
		vvr1(6)=0.108103018168071 ;vvr2(6)=0.445948490915964 ;vvr3(6)=0.445948490915964 ;vvwg(6)=0.223381589678000

case(4)
		
vvr1(1)=0.888871894660413 ;vvr2(1)=0.055564052669793 ;vvr3(1)=0.055564052669793 ;vvwg(1)=0.041955512996649
vvr1(2)=0.055564052669793 ;vvr2(2)=0.888871894660413 ;vvr3(2)=0.055564052669793 ;vvwg(2)=0.041955512996649
vvr1(3)=0.055564052669793 ;vvr2(3)=0.055564052669793 ;vvr3(3)=0.888871894660413 ;vvwg(3)=0.041955512996649
vvr1(4)=0.295533711735893 ;vvr2(4)=0.634210747745723 ;vvr3(4)=0.070255540518384 ;vvwg(4)=0.112098412070887
vvr1(5)=0.295533711735893 ;vvr2(5)=0.070255540518384 ;vvr3(5)=0.634210747745723 ;vvwg(5)=0.112098412070887
vvr1(6)=0.070255540518384 ;vvr2(6)=0.295533711735893 ;vvr3(6)=0.634210747745723 ;vvwg(6)=0.112098412070887
vvr1(7)=0.634210747745723 ;vvr2(7)=0.295533711735893 ;vvr3(7)=0.070255540518384 ;vvwg(7)=0.112098412070887
vvr1(8)=0.634210747745723 ;vvr2(8)=0.070255540518384 ;vvr3(8)=0.295533711735893 ;vvwg(8)=0.112098412070887
vvr1(9)=0.070255540518384 ;vvr2(9)=0.634210747745723 ;vvr3(9)=0.295533711735893 ;vvwg(9)=0.112098412070887
vvr1(10)=0.333333333333333 ;vvr2(10)=0.333333333333333 ;vvr3(10)=0.333333333333333 ;vvwg(10)=0.201542988584730


case(5)

vvr1(1)=0.928258244608533;vvr2(1)= 0.035870877695734 ;vvr3(1)=0.035870877695734 ;vvwg(1)=0.017915455012303
vvr1(2)=0.035870877695734;vvr2(2)= 0.928258244608533 ;vvr3(2)=0.035870877695734 ;vvwg(2)=0.017915455012303
vvr1(3)=0.035870877695734;vvr2(3)= 0.035870877695734 ;vvr3(3)=0.928258244608533 ;vvwg(3)=0.017915455012303
vvr1(4)=0.516541208464066;vvr2(4)= 0.241729395767967 ;vvr3(4)=0.241729395767967 ;vvwg(4)=0.127712195881265
vvr1(5)=0.241729395767967;vvr2(5)= 0.516541208464066 ;vvr3(5)=0.241729395767967 ;vvwg(5)=0.127712195881265
vvr1(6)=0.241729395767967;vvr2(6)= 0.241729395767967 ;vvr3(6)=0.516541208464066 ;vvwg(6)=0.127712195881265
vvr1(7)=0.474308787777079;vvr2(7)= 0.474308787777079 ;vvr3(7)=0.051382424445843 ;vvwg(7)=0.076206062385535
vvr1(8)=0.474308787777079;vvr2(8)= 0.051382424445843 ;vvr3(8)=0.474308787777079 ;vvwg(8)=0.076206062385535
vvr1(9)=0.051382424445843;vvr2(9)= 0.474308787777079 ;vvr3(9)=0.474308787777079 ;vvwg(9)=0.076206062385535
vvr1(10)=0.201503881881800;vvr2(10)= 0.751183631106484 ;vvr3(10)=0.047312487011716 ;vvwg(10)=0.055749810027115
vvr1(11)=0.201503881881800;vvr2(11)= 0.047312487011716 ;vvr3(11)=0.751183631106484 ;vvwg(11)=0.055749810027115
vvr1(12)=0.047312487011716;vvr2(12)= 0.201503881881800 ;vvr3(12)=0.751183631106484 ;vvwg(12)=0.055749810027115
vvr1(13)=0.751183631106484;vvr2(13)= 0.201503881881800 ;vvr3(13)=0.047312487011716 ;vvwg(13)=0.055749810027115
vvr1(14)=0.751183631106484;vvr2(14)= 0.047312487011716 ;vvr3(14)=0.201503881881800 ;vvwg(14)=0.055749810027115
vvr1(15)=0.047312487011716;vvr2(15)= 0.751183631106484 ;vvr3(15)=0.201503881881800 ;vvwg(15)=0.055749810027115








case(6)

vvr1(1)=0.943774095634672 ;vvr2(1)=0.028112952182664 ;vvr3(1)=0.028112952182664 ;vvwg(1)=0.010359374696538
vvr1(2)=0.028112952182664 ;vvr2(2)=0.943774095634672 ;vvr3(2)=0.028112952182664 ;vvwg(2)=0.010359374696538
vvr1(3)=0.028112952182664 ;vvr2(3)=0.028112952182664 ;vvr3(3)=0.943774095634672 ;vvwg(3)=0.010359374696538
vvr1(4)=0.645721803061365 ;vvr2(4)=0.177139098469317 ;vvr3(4)=0.177139098469317 ;vvwg(4)=0.075394884326738
vvr1(5)=0.177139098469317 ;vvr2(5)=0.645721803061365 ;vvr3(5)=0.177139098469317 ;vvwg(5)=0.075394884326738
vvr1(6)=0.177139098469317 ;vvr2(6)=0.177139098469317 ;vvr3(6)=0.645721803061365 ;vvwg(6)=0.075394884326738
vvr1(7)=0.405508595867433 ;vvr2(7)=0.405508595867433 ;vvr3(7)=0.188982808265134 ;vvwg(7)=0.097547802373242
vvr1(8)=0.405508595867433 ;vvr2(8)=0.188982808265134 ;vvr3(8)=0.405508595867433 ;vvwg(8)=0.097547802373242
vvr1(9)=0.188982808265134 ;vvr2(9)=0.405508595867433 ;vvr3(9)=0.405508595867433 ;vvwg(9)=0.097547802373242
vvr1(10)=0.148565812270887 ;vvr2(10)=0.817900980028499 ;vvr3(10)=0.033533207700614 ;vvwg(10)=0.028969269372473
vvr1(11)=0.148565812270887 ;vvr2(11)=0.033533207700614 ;vvr3(11)=0.817900980028499 ;vvwg(11)=0.028969269372473
vvr1(12)=0.033533207700614 ;vvr2(12)=0.148565812270887 ;vvr3(12)=0.817900980028499 ;vvwg(12)=0.028969269372473
vvr1(13)=0.817900980028499 ;vvr2(13)=0.148565812270887 ;vvr3(13)=0.033533207700614 ;vvwg(13)=0.028969269372473
vvr1(14)=0.817900980028499 ;vvr2(14)=0.033533207700614 ;vvr3(14)=0.148565812270887 ;vvwg(14)=0.028969269372473
vvr1(15)=0.033533207700614 ;vvr2(15)=0.817900980028499 ;vvr3(15)=0.148565812270887 ;vvwg(15)=0.028969269372473
vvr1(16)=0.357196298615681 ;vvr2(16)=0.604978911775132 ;vvr3(16)=0.037824789609186 ;vvwg(16)=0.046046366595935
vvr1(17)=0.357196298615681 ;vvr2(17)=0.037824789609186 ;vvr3(17)=0.604978911775132 ;vvwg(17)=0.046046366595935
vvr1(18)=0.037824789609186 ;vvr2(18)=0.357196298615681 ;vvr3(18)=0.604978911775132 ;vvwg(18)=0.046046366595935
vvr1(19)=0.604978911775132 ;vvr2(19)=0.357196298615681 ;vvr3(19)=0.037824789609186 ;vvwg(19)=0.046046366595935
vvr1(20)=0.604978911775132 ;vvr2(20)=0.037824789609186 ;vvr3(20)=0.357196298615681 ;vvwg(20)=0.046046366595935
vvr1(21)=0.037824789609186 ;vvr2(21)=0.604978911775132 ;vvr3(21)=0.357196298615681 ;vvwg(21)=0.046046366595935

case(7,8,9)

vvr1(1)=0.957657154441070
vvr1(2)=0.021171422779465
vvr1(3)=0.021171422779465
vvr1(4)=0.798831205208225
vvr1(5)=0.100584397395888
vvr1(6)=0.100584397395888
vvr1(7)=0.457923384576135
vvr1(8)=0.271038307711932
vvr1(9)=0.271038307711932
vvr1(10)=0.440191258403832
vvr1(11)=0.440191258403832
vvr1(12)=0.119617483192335
vvr1(13)=0.101763679498021
vvr1(14)=0.101763679498021
vvr1(15)=0.018256679074748
vvr1(16)=0.879979641427232
vvr1(17)=0.879979641427232
vvr1(18)=0.018256679074748
vvr1(19)=0.394033271669987
vvr1(20)=0.394033271669987
vvr1(21)=0.023404705466341
vvr1(22)=0.582562022863673
vvr1(23)=0.582562022863673
vvr1(24)=0.023404705466341
vvr1(25)=0.226245530909229
vvr1(26)=0.226245530909229
vvr1(27)=0.022223854547989
vvr1(28)=0.751530614542782
vvr1(29)=0.751530614542782
vvr1(30)=0.022223854547989
vvr1(31)=0.635737183263105
vvr1(32)=0.635737183263105
vvr1(33)=0.115183589115563
vvr1(34)=0.249079227621332
vvr1(35)=0.249079227621332
vvr1(36)=0.115183589115563



vvr2(1)=0.021171422779465
vvr2(2)=0.957657154441070
vvr2(3)=0.021171422779465
vvr2(4)=0.100584397395888
vvr2(5)=0.798831205208225
vvr2(6)=0.100584397395888
vvr2(7)=0.271038307711932
vvr2(8)=0.457923384576135
vvr2(9)=0.271038307711932
vvr2(10)=0.440191258403832
vvr2(11)=0.119617483192335
vvr2(12)=0.440191258403832
vvr2(13)=0.879979641427232
vvr2(14)=0.018256679074748
vvr2(15)=0.101763679498021
vvr2(16)=0.101763679498021
vvr2(17)=0.018256679074748
vvr2(18)=0.879979641427232
vvr2(19)=0.582562022863673
vvr2(20)=0.023404705466341
vvr2(21)=0.394033271669987
vvr2(22)=0.394033271669987
vvr2(23)=0.023404705466341
vvr2(24)=0.582562022863673
vvr2(25)=0.751530614542782
vvr2(26)=0.022223854547989
vvr2(27)=0.226245530909229
vvr2(28)=0.226245530909229
vvr2(29)=0.022223854547989
vvr2(30)=0.751530614542782
vvr2(31)=0.249079227621332
vvr2(32)=0.115183589115563
vvr2(33)=0.635737183263105
vvr2(34)=0.635737183263105
vvr2(35)=0.115183589115563
vvr2(36)=0.249079227621332




vvr3(1)=0.021171422779465
vvr3(2)=0.021171422779465
vvr3(3)=0.957657154441070
vvr3(4)=0.100584397395888
vvr3(5)=0.100584397395888
vvr3(6)=0.798831205208225
vvr3(7)=0.271038307711932
vvr3(8)=0.271038307711932
vvr3(9)=0.457923384576135
vvr3(10)=0.119617483192335
vvr3(11)=0.440191258403832
vvr3(12)=0.440191258403832
vvr3(13)=0.018256679074748
vvr3(14)=0.879979641427232
vvr3(15)=0.879979641427232
vvr3(16)=0.018256679074748
vvr3(17)=0.101763679498021
vvr3(18)=0.101763679498021
vvr3(19)=0.023404705466341
vvr3(20)=0.582562022863673
vvr3(21)=0.582562022863673
vvr3(22)=0.023404705466341
vvr3(23)=0.394033271669987
vvr3(24)=0.394033271669987
vvr3(25)=0.022223854547989
vvr3(26)=0.751530614542782
vvr3(27)=0.751530614542782
vvr3(28)=0.022223854547989
vvr3(29)=0.226245530909229
vvr3(30)=0.226245530909229
vvr3(31)=0.115183589115563
vvr3(32)=0.249079227621332
vvr3(33)=0.249079227621332
vvr3(34)=0.115183589115563
vvr3(35)=0.635737183263105
vvr3(36)=0.635737183263105



vvwg(1)=0.005639123786910
vvwg(2)=0.005639123786910
vvwg(3)=0.005639123786910
vvwg(4)=0.027148968192278
vvwg(5)=0.027148968192278
vvwg(6)=0.027148968192278
vvwg(7)=0.063100912533359
vvwg(8)=0.063100912533359
vvwg(9)=0.063100912533359
vvwg(10)=0.051752795679899
vvwg(11)=0.051752795679899
vvwg(12)=0.051752795679899
vvwg(13)=0.009866753574646
vvwg(14)=0.009866753574646
vvwg(15)=0.009866753574646
vvwg(16)=0.009866753574646
vvwg(17)=0.009866753574646
vvwg(18)=0.009866753574646
vvwg(19)=0.022008204800147
vvwg(20)=0.022008204800147
vvwg(21)=0.022008204800147
vvwg(22)=0.022008204800147
vvwg(23)=0.022008204800147
vvwg(24)=0.022008204800147
vvwg(25)=0.016644570076736
vvwg(26)=0.016644570076736
vvwg(27)=0.016644570076736
vvwg(28)=0.016644570076736
vvwg(29)=0.016644570076736
vvwg(30)=0.016644570076736
vvwg(31)=0.044326238118914
vvwg(32)=0.044326238118914
vvwg(33)=0.044326238118914
vvwg(34)=0.044326238118914
vvwg(35)=0.044326238118914
vvwg(36)=0.044326238118914
		
		
end select



		do kk=1,qp_triangle
			wequa3d(kk)=vvwg(kk)
			qpoints(:,kk)=(vvr1(kk)*vext(1,:))+(vvr2(kk)*vext(2,:))+(vvr3(kk)*vext(3,:))
			
		end do




end subroutine quadraturetriangle


subroutine quadraturequad(n,igqrules,vext,qpoints,wequa3d)
 !> @brief
!> this subroutine computes the quadrature points and weights for quadrilateral in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints),intent(inout)::qpoints
real,dimension(1:numberofpoints),intent(inout)::wequa3d
real,dimension(1:2,1:2)::vva,vva1
real,dimension(1)::deta
real,dimension(1:4)::vvnxi
real,dimension(1:alls)::vvwg,vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
real::r,s,tx,a,b,c,d,e,f
real::a1,b1,c1,d1,e1,f1
integer::kk,j,ii,ij,ik,count1


 wequa3d=0.0d0
  qpoints=0.0d0

select case(igqrules)
 

 case(1)

		vvwg(1) = 4.0d0
	    vvr1(1)=0.0d0	;vvr2(1)=0.0d0	


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
  vvnpoy(1)= a	;vvnpoy(2)=b	
  vvnpoz(1)= a	;vvnpoz(2)=b

  vvwpox(1)=a1	;vvwpox(2)=b1	
  vvwpoy(1)= a1	;vvwpoy(2)=b1	
  vvwpoz(1)= a1	;vvwpoz(2)=b1
  
  count1=0
 do ii=1,2
    do ij=1,2
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
  	


  	
 case(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1
  count1=0
 do ii=1,3
    do ij=1,3
     
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
 

			
 case(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1
  count1=0
 do ii=1,4
    do ij=1,4
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
     
    end do
end do


			
case(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  e1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1
  count1=0
 do ii=1,5
    do ij=1,5
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=(vvwpox(ii)*vvwpox(ij))
      
    end do
end do
	
			
case(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  f=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  f1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e;vvnpoy(6)=f
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e;vvnpoz(6)=f
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1 ;vvwpoy(6)=f1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1 ;vvwpoz(6)=f1
  count1=0
 do ii=1,6
    do ij=1,6
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
	
			
end select
		qpoints(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.25d0
! 		  wequa3d(:)=vvwg(:)
		do kk=1,qp_quad
			 wequa3d(kk)=vvwg(kk)
			r=vvr1(kk); s=vvr2(kk);
			vvnxi(1)=(0.25d0)*(1.0d0-r)*(1.0d0-s)
			vvnxi(2)=(0.25d0)*(1.0d0+r)*(1.0d0-s)
			vvnxi(3)=(0.25d0)*(1.0d0+r)*(1.0d0+s)
			vvnxi(4)=(0.25d0)*(1.0d0-r)*(1.0d0+s)
			
			do j=1,4
			qpoints(1:2,kk)=qpoints(1:2,kk)+(vvnxi(j)*vext(j,1:2))
			end do
! 			

		end do
		




end subroutine quadraturequad




subroutine quadraturequad3d(n,igqrules,vext,qpoints2d,wequa2d)
 !> @brief
!> this subroutine computes the quadrature points and weights for quadrilateral in 3d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,igqrules
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints2),intent(inout)::qpoints2d
real,dimension(1:numberofpoints2),intent(inout)::wequa2d
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1)::deta
real,dimension(1:4)::vvnxi
real,dimension(1:alls)::vvwg
real,dimension(1:alls)::vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
real::r,s,a,b,c,d,e,f
real::a1,b1,c1,d1,e1,f1
integer::kk,j,ii,ij,ik,count1



 wequa2d=0.0d0
  qpoints2d=0.0d0

select case(igqrules)
 

 case(1)

		vvwg(1) = 4.0d0
	    vvr1(1)=0.0d0	;vvr2(1)=0.0d0	


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
  vvnpoy(1)= a	;vvnpoy(2)=b	
  vvnpoz(1)= a	;vvnpoz(2)=b

  vvwpox(1)=a1	;vvwpox(2)=b1	
  vvwpoy(1)= a1	;vvwpoy(2)=b1	
  vvwpoz(1)= a1	;vvwpoz(2)=b1
  
  count1=0
 do ii=1,2
    do ij=1,2
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
  	


  	
 case(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1
  count1=0
 do ii=1,3
    do ij=1,3
     
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
 

			
 case(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1
  count1=0
 do ii=1,4
    do ij=1,4
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
     
    end do
end do


			
case(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  e1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1
  count1=0
 do ii=1,5
    do ij=1,5
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=(vvwpox(ii)*vvwpox(ij))
      
    end do
end do
	
			
case(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  f=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  f1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e;vvnpoy(6)=f
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e;vvnpoz(6)=f
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1 ;vvwpoy(6)=f1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1 ;vvwpoz(6)=f1
  count1=0
 do ii=1,6
    do ij=1,6
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
	
			
end select
		qpoints2d(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.25d0
! 		  wequa2d(:)=vvwg(:)
		do kk=1,qp_quad
			wequa2d(kk)=vvwg(kk)
			r=vvr1(kk); s=vvr2(kk);
			vvnxi(1)=(0.25d0)*(1.0d0-r)*(1.0d0-s)
			vvnxi(2)=(0.25d0)*(1.0d0+r)*(1.0d0-s)
			vvnxi(3)=(0.25d0)*(1.0d0+r)*(1.0d0+s)
			vvnxi(4)=(0.25d0)*(1.0d0-r)*(1.0d0+s)
			
			do j=1,4
			qpoints2d(1:3,kk)=qpoints2d(1:3,kk)+(vvnxi(j)*vext(j,1:3))
			end do
! 			

		end do
		




end subroutine quadraturequad3d


subroutine quadratureline(n,igqrules,vext,qpoints2d,wequa2d)
 !> @brief
!> this subroutine computes the quadrature points for a line and returns it in qpoints2d(dim,qp)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,igqrules
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints2),intent(inout)::qpoints2d
real,dimension(1:numberofpoints2),intent(inout)::wequa2d
real,dimension(1:2,1:2)::vva,vva1
real,dimension(1)::deta
real,dimension(1:4)::vvnxi
real,dimension(1:alls)::vvwg
real,dimension(1:alls)::vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
real::r,s,tx,a,b,c,d,e,f,g,h,k
real::a1,b1,c1,d1,e1,f1,g1,h1,k1
integer::kk,j,ii,ij,ik,count1



 wequa2d=0.0d0
  qpoints2d=0.0d0

select case(igqrules)
 

 case(1)

		vvwg(1) = 2.0d0
	    vvr1(1)=0.0d0	;vvr2(1)=0.0d0	


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
 

  vvwpox(1)=a1	;vvwpox(2)=b1	
  
  
  count1=0
 do ii=1,2
         
	count1=count1+1
	vvr1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
    
end do
  	


  	
 case(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  
  count1=0
 do ii=1,3
    
     
	count1=count1+1
	vvr1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
    
end do
 

			
 case(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  
  count1=0
 do ii=1,4
   
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
     
   
end do


			
case(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  e1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
 
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  
  count1=0
 do ii=1,5
   
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii)
	vvwg(count1)=(vvwpox(ii))
      
    
end do
	
			
case(6)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  f=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  f1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  
  count1=0
 do ii=1,6
    
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
   
end do
	
	
	
case(7,8,9)

  
  
  a1=0.3302393550012598	         ;a=0.0000000000000000
	b1=0.1806481606948574	;b=-0.8360311073266358
	c1=0.1806481606948574	;c=0.8360311073266358
	d1=0.0812743883615744	;d=-0.9681602395076261
	e1=0.0812743883615744	;e=0.9681602395076261
	f1=0.3123470770400029	;f=-0.3242534234038089
	g1=0.3123470770400029	;g=0.3242534234038089
	h1=0.2606106964029354	;h=-0.6133714327005904
	k1=0.2606106964029354	;k=0.6133714327005904
  
  
  
  
  
  
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e; vvnpox(6)=f ; vvnpox(7)=g ; vvnpox(8)=h ; vvnpox(9)=k
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1 ;vvwpox(7)=g1 ;vvwpox(8)=h1 ;vvwpox(9)=k1
  
  count1=0
 do ii=1,9
    
      
	count1=count1+1
	vvr1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
   
end do	
	
	
			
end select
		qpoints2d(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.5d0

		do kk=1,qp_line
			wequa2d(kk)=vvwg(kk)
			r=vvr1(kk); 
				
			
			
			qpoints2d(:,kk)=((vext(1,1:2)+vext(2,1:2))/2.0d0)+(r*(vext(2,1:2)-vext(1,1:2))/2.0d0)
			
			
! 		
		end do
		




end subroutine quadratureline



subroutine quadraturetetra(n,igqrules,vext,qpoints,wequa3d)
 !> @brief
!> this subroutine computes the quadrature points and weights for a tetrahedral
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints),intent(inout)::qpoints
real,dimension(1:numberofpoints),intent(inout)::wequa3d
integer::kk
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1)::deta
real,dimension(1:8)::vvnxi
real,dimension(1:alls)::vvwg,vvr1,vvr2,vvr3,vvr4
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz



wequa3d=0.0d0
qpoints=0.0d0


select case(igqrules)
case(1)

		vvwg(1) = 1.0d0
	    vvr1(1)=0.25d0	;vvr2(1)=0.25d0	;vvr3(1)=0.25d0; vvr4(1)=0.25d0
case(2)
	
  	

vvr1(1)=0.5854101966249680 ;vvr2(1)=0.1381966011250110 ;vvr3(1)=0.1381966011250110 ;vvr4(1)=0.1381966011250110 ;vvwg(1)=0.2500000000000000
vvr1(2)=0.1381966011250110 ;vvr2(2)=0.5854101966249680 ;vvr3(2)=0.1381966011250110 ;vvr4(2)=0.1381966011250110 ;vvwg(2)=0.2500000000000000
vvr1(3)=0.1381966011250110 ;vvr2(3)=0.1381966011250110 ;vvr3(3)=0.5854101966249680 ;vvr4(3)=0.1381966011250110 ;vvwg(3)=0.2500000000000000
vvr1(4)=0.1381966011250110 ;vvr2(4)=0.1381966011250110 ;vvr3(4)=0.1381966011250110 ;vvr4(4)=0.5854101966249680 ;vvwg(4)=0.2500000000000000

  	
case(3)

  
vvr1( 1)= 0.3108859192633006e+00;vvr2( 1)= 0.3108859192633006e+00;vvr3( 1)= 0.6734224221009816e-01;vvr4( 1)= 0.3108859192633006e+00;vvwg( 1)= 0.1126879257180159e+00
vvr1( 2)= 0.3108859192633006e+00;vvr2( 2)= 0.6734224221009816e-01;vvr3( 2)= 0.3108859192633006e+00;vvr4( 2)= 0.3108859192633005e+00;vvwg( 2)= 0.1126879257180159e+00
vvr1( 3)= 0.6734224221009816e-01;vvr2( 3)= 0.3108859192633006e+00;vvr3( 3)= 0.3108859192633006e+00;vvr4( 3)= 0.3108859192633005e+00;vvwg( 3)= 0.1126879257180159e+00
vvr1( 4)= 0.3108859192633006e+00;vvr2( 4)= 0.3108859192633006e+00;vvr3( 4)= 0.3108859192633006e+00;vvr4( 4)= 0.6734224221009810e-01;vvwg( 4)= 0.1126879257180159e+00
vvr1( 5)= 0.9273525031089125e-01;vvr2( 5)= 0.9273525031089125e-01;vvr3( 5)= 0.7217942490673264e+00;vvr4( 5)= 0.9273525031089114e-01;vvwg( 5)= 0.7349304311636196e-01
vvr1( 6)= 0.9273525031089125e-01;vvr2( 6)= 0.7217942490673264e+00;vvr3( 6)= 0.9273525031089125e-01;vvr4( 6)= 0.9273525031089114e-01;vvwg( 6)= 0.7349304311636196e-01
vvr1( 7)= 0.7217942490673264e+00;vvr2( 7)= 0.9273525031089125e-01;vvr3( 7)= 0.9273525031089125e-01;vvr4( 7)= 0.9273525031089114e-01;vvwg( 7)= 0.7349304311636196e-01
vvr1( 8)= 0.9273525031089125e-01;vvr2( 8)= 0.9273525031089125e-01;vvr3( 8)= 0.9273525031089125e-01;vvr4( 8)= 0.7217942490673263e+00;vvwg( 8)= 0.7349304311636196e-01
vvr1( 9)= 0.4550370412564964e-01;vvr2( 9)= 0.4544962958743504e+00;vvr3( 9)= 0.4544962958743504e+00;vvr4( 9)= 0.4550370412564964e-01;vvwg( 9)= 0.4254602077708147e-01
vvr1(10)= 0.4544962958743504e+00;vvr2(10)= 0.4550370412564964e-01;vvr3(10)= 0.4544962958743504e+00;vvr4(10)= 0.4550370412564964e-01;vvwg(10)= 0.4254602077708147e-01
vvr1(11)= 0.4550370412564964e-01;vvr2(11)= 0.4550370412564964e-01;vvr3(11)= 0.4544962958743504e+00;vvr4(11)= 0.4544962958743504e+00;vvwg(11)= 0.4254602077708147e-01
vvr1(12)= 0.4550370412564964e-01;vvr2(12)= 0.4544962958743504e+00;vvr3(12)= 0.4550370412564964e-01;vvr4(12)= 0.4544962958743504e+00;vvwg(12)= 0.4254602077708147e-01
vvr1(13)= 0.4544962958743504e+00;vvr2(13)= 0.4550370412564964e-01;vvr3(13)= 0.4550370412564964e-01;vvr4(13)= 0.4544962958743504e+00;vvwg(13)= 0.4254602077708147e-01
vvr1(14)= 0.4544962958743504e+00;vvr2(14)= 0.4544962958743504e+00;vvr3(14)= 0.4550370412564964e-01;vvr4(14)= 0.4550370412564964e-01;vvwg(14)= 0.4254602077708147e-01






  












			
case(4)
vvr1( 1)= 0.4067395853461137e-01;vvr2( 1)= 0.4067395853461137e-01;vvr3( 1)= 0.8779781243961660e+00;vvr4( 1)= 0.4067395853461120e-01;vvwg( 1)= 0.1007721105532064e-01
vvr1( 2)= 0.4067395853461137e-01;vvr2( 2)= 0.8779781243961660e+00;vvr3( 2)= 0.4067395853461137e-01;vvr4( 2)= 0.4067395853461125e-01;vvwg( 2)= 0.1007721105532064e-01
vvr1( 3)= 0.8779781243961660e+00;vvr2( 3)= 0.4067395853461137e-01;vvr3( 3)= 0.4067395853461137e-01;vvr4( 3)= 0.4067395853461131e-01;vvwg( 3)= 0.1007721105532064e-01
vvr1( 4)= 0.4067395853461137e-01;vvr2( 4)= 0.4067395853461137e-01;vvr3( 4)= 0.4067395853461137e-01;vvr4( 4)= 0.8779781243961657e+00;vvwg( 4)= 0.1007721105532064e-01
vvr1( 5)= 0.3223378901422755e+00;vvr2( 5)= 0.3223378901422755e+00;vvr3( 5)= 0.3298632957317349e-01;vvr4( 5)= 0.3223378901422754e+00;vvwg( 5)= 0.5535718154365472e-01
vvr1( 6)= 0.3223378901422755e+00;vvr2( 6)= 0.3298632957317349e-01;vvr3( 6)= 0.3223378901422755e+00;vvr4( 6)= 0.3223378901422754e+00;vvwg( 6)= 0.5535718154365472e-01
vvr1( 7)= 0.3298632957317349e-01;vvr2( 7)= 0.3223378901422755e+00;vvr3( 7)= 0.3223378901422755e+00;vvr4( 7)= 0.3223378901422754e+00;vvwg( 7)= 0.5535718154365472e-01
vvr1( 8)= 0.3223378901422755e+00;vvr2( 8)= 0.3223378901422755e+00;vvr3( 8)= 0.3223378901422755e+00;vvr4( 8)= 0.3298632957317338e-01;vvwg( 8)= 0.5535718154365472e-01
vvr1( 9)= 0.2146028712591520e+00;vvr2( 9)= 0.2146028712591520e+00;vvr3( 9)= 0.3561913862225439e+00;vvr4( 9)= 0.2146028712591521e+00;vvwg( 9)= 0.3992275025816749e-01
vvr1(10)= 0.2146028712591520e+00;vvr2(10)= 0.3561913862225439e+00;vvr3(10)= 0.2146028712591520e+00;vvr4(10)= 0.2146028712591521e+00;vvwg(10)= 0.3992275025816749e-01
vvr1(11)= 0.3561913862225439e+00;vvr2(11)= 0.2146028712591520e+00;vvr3(11)= 0.2146028712591520e+00;vvr4(11)= 0.2146028712591521e+00;vvwg(11)= 0.3992275025816749e-01
vvr1(12)= 0.2146028712591520e+00;vvr2(12)= 0.2146028712591520e+00;vvr3(12)= 0.2146028712591520e+00;vvr4(12)= 0.3561913862225440e+00;vvwg(12)= 0.3992275025816749e-01
vvr1(13)= 0.6030056647916491e+00;vvr2(13)= 0.6366100187501750e-01;vvr3(13)= 0.2696723314583158e+00;vvr4(13)= 0.6366100187501761e-01;vvwg(13)= 0.4821428571428571e-01
vvr1(14)= 0.6030056647916491e+00;vvr2(14)= 0.6366100187501750e-01;vvr3(14)= 0.6366100187501750e-01;vvr4(14)= 0.2696723314583159e+00;vvwg(14)= 0.4821428571428571e-01
vvr1(15)= 0.6366100187501750e-01;vvr2(15)= 0.6366100187501750e-01;vvr3(15)= 0.6030056647916491e+00;vvr4(15)= 0.2696723314583160e+00;vvwg(15)= 0.4821428571428571e-01
vvr1(16)= 0.2696723314583158e+00;vvr2(16)= 0.6030056647916491e+00;vvr3(16)= 0.6366100187501750e-01;vvr4(16)= 0.6366100187501761e-01;vvwg(16)= 0.4821428571428571e-01
vvr1(17)= 0.6366100187501750e-01;vvr2(17)= 0.2696723314583158e+00;vvr3(17)= 0.6030056647916491e+00;vvr4(17)= 0.6366100187501766e-01;vvwg(17)= 0.4821428571428571e-01
vvr1(18)= 0.6366100187501750e-01;vvr2(18)= 0.6030056647916491e+00;vvr3(18)= 0.6366100187501750e-01;vvr4(18)= 0.2696723314583160e+00;vvwg(18)= 0.4821428571428571e-01
vvr1(19)= 0.2696723314583158e+00;vvr2(19)= 0.6366100187501750e-01;vvr3(19)= 0.6030056647916491e+00;vvr4(19)= 0.6366100187501766e-01;vvwg(19)= 0.4821428571428571e-01
vvr1(20)= 0.6366100187501750e-01;vvr2(20)= 0.2696723314583158e+00;vvr3(20)= 0.6366100187501750e-01;vvr4(20)= 0.6030056647916493e+00;vvwg(20)= 0.4821428571428571e-01
vvr1(21)= 0.6366100187501750e-01;vvr2(21)= 0.6366100187501750e-01;vvr3(21)= 0.2696723314583158e+00;vvr4(21)= 0.6030056647916493e+00;vvwg(21)= 0.4821428571428571e-01
vvr1(22)= 0.6366100187501750e-01;vvr2(22)= 0.6030056647916491e+00;vvr3(22)= 0.2696723314583158e+00;vvr4(22)= 0.6366100187501766e-01;vvwg(22)= 0.4821428571428571e-01
vvr1(23)= 0.2696723314583158e+00;vvr2(23)= 0.6366100187501750e-01;vvr3(23)= 0.6366100187501750e-01;vvr4(23)= 0.6030056647916493e+00;vvwg(23)= 0.4821428571428571e-01
vvr1(24)= 0.6030056647916491e+00;vvr2(24)= 0.2696723314583158e+00;vvr3(24)= 0.6366100187501750e-01;vvr4(24)= 0.6366100187501761e-01;vvwg(24)= 0.4821428571428571e-01




			
case(5)


vvr1( 1)= 0.2500000000000000e+00;vvr2( 1)= 0.2500000000000000e+00;vvr3( 1)= 0.2500000000000000e+00;vvr4( 1)= 0.2500000000000000e+00;vvwg( 1)= 0.9548528946413085e-01
vvr1( 2)= 0.3157011497782028e+00;vvr2( 2)= 0.3157011497782028e+00;vvr3( 2)= 0.5289655066539162e-01;vvr4( 2)= 0.3157011497782028e+00;vvwg( 2)= 0.4232958120996703e-01
vvr1( 3)= 0.3157011497782028e+00;vvr2( 3)= 0.5289655066539162e-01;vvr3( 3)= 0.3157011497782028e+00;vvr4( 3)= 0.3157011497782029e+00;vvwg( 3)= 0.4232958120996703e-01
vvr1( 4)= 0.5289655066539162e-01;vvr2( 4)= 0.3157011497782028e+00;vvr3( 4)= 0.3157011497782028e+00;vvr4( 4)= 0.3157011497782029e+00;vvwg( 4)= 0.4232958120996703e-01
vvr1( 5)= 0.3157011497782028e+00;vvr2( 5)= 0.3157011497782028e+00;vvr3( 5)= 0.3157011497782028e+00;vvr4( 5)= 0.5289655066539167e-01;vvwg( 5)= 0.4232958120996703e-01
vvr1( 6)= 0.5048982259839635e-01;vvr2( 6)= 0.4495101774016036e+00;vvr3( 6)= 0.4495101774016036e+00;vvr4( 6)= 0.5048982259839652e-01;vvwg( 6)= 0.3189692783285758e-01
vvr1( 7)= 0.4495101774016036e+00;vvr2( 7)= 0.5048982259839635e-01;vvr3( 7)= 0.4495101774016036e+00;vvr4( 7)= 0.5048982259839641e-01;vvwg( 7)= 0.3189692783285758e-01
vvr1( 8)= 0.5048982259839635e-01;vvr2( 8)= 0.5048982259839635e-01;vvr3( 8)= 0.4495101774016036e+00;vvr4( 8)= 0.4495101774016038e+00;vvwg( 8)= 0.3189692783285758e-01
vvr1( 9)= 0.5048982259839635e-01;vvr2( 9)= 0.4495101774016036e+00;vvr3( 9)= 0.5048982259839635e-01;vvr4( 9)= 0.4495101774016038e+00;vvwg( 9)= 0.3189692783285758e-01
vvr1(10)= 0.4495101774016036e+00;vvr2(10)= 0.5048982259839635e-01;vvr3(10)= 0.5048982259839635e-01;vvr4(10)= 0.4495101774016036e+00;vvwg(10)= 0.3189692783285758e-01
vvr1(11)= 0.4495101774016036e+00;vvr2(11)= 0.4495101774016036e+00;vvr3(11)= 0.5048982259839635e-01;vvr4(11)= 0.5048982259839646e-01;vvwg(11)= 0.3189692783285758e-01
vvr1(12)= 0.5751716375870000e+00;vvr2(12)= 0.1888338310260010e+00;vvr3(12)= 0.4716070036099790e-01;vvr4(12)= 0.1888338310260011e+00;vvwg(12)= 0.3720713072833462e-01
vvr1(13)= 0.5751716375870000e+00;vvr2(13)= 0.1888338310260010e+00;vvr3(13)= 0.1888338310260010e+00;vvr4(13)= 0.4716070036099795e-01;vvwg(13)= 0.3720713072833462e-01
vvr1(14)= 0.1888338310260010e+00;vvr2(14)= 0.1888338310260010e+00;vvr3(14)= 0.5751716375870000e+00;vvr4(14)= 0.4716070036099795e-01;vvwg(14)= 0.3720713072833462e-01
vvr1(15)= 0.4716070036099790e-01;vvr2(15)= 0.5751716375870000e+00;vvr3(15)= 0.1888338310260010e+00;vvr4(15)= 0.1888338310260010e+00;vvwg(15)= 0.3720713072833462e-01
vvr1(16)= 0.1888338310260010e+00;vvr2(16)= 0.4716070036099790e-01;vvr3(16)= 0.5751716375870000e+00;vvr4(16)= 0.1888338310260012e+00;vvwg(16)= 0.3720713072833462e-01
vvr1(17)= 0.1888338310260010e+00;vvr2(17)= 0.5751716375870000e+00;vvr3(17)= 0.1888338310260010e+00;vvr4(17)= 0.4716070036099795e-01;vvwg(17)= 0.3720713072833462e-01
vvr1(18)= 0.4716070036099790e-01;vvr2(18)= 0.1888338310260010e+00;vvr3(18)= 0.5751716375870000e+00;vvr4(18)= 0.1888338310260010e+00;vvwg(18)= 0.3720713072833462e-01
vvr1(19)= 0.1888338310260010e+00;vvr2(19)= 0.4716070036099790e-01;vvr3(19)= 0.1888338310260010e+00;vvr4(19)= 0.5751716375870001e+00;vvwg(19)= 0.3720713072833462e-01
vvr1(20)= 0.1888338310260010e+00;vvr2(20)= 0.1888338310260010e+00;vvr3(20)= 0.4716070036099790e-01;vvr4(20)= 0.5751716375870000e+00;vvwg(20)= 0.3720713072833462e-01
vvr1(21)= 0.1888338310260010e+00;vvr2(21)= 0.5751716375870000e+00;vvr3(21)= 0.4716070036099790e-01;vvr4(21)= 0.1888338310260011e+00;vvwg(21)= 0.3720713072833462e-01
vvr1(22)= 0.4716070036099790e-01;vvr2(22)= 0.1888338310260010e+00;vvr3(22)= 0.1888338310260010e+00;vvr4(22)= 0.5751716375870000e+00;vvwg(22)= 0.3720713072833462e-01
vvr1(23)= 0.5751716375870000e+00;vvr2(23)= 0.4716070036099790e-01;vvr3(23)= 0.1888338310260010e+00;vvr4(23)= 0.1888338310260011e+00;vvwg(23)= 0.3720713072833462e-01
vvr1(24)= 0.8108302410985486e+00;vvr2(24)= 0.2126547254148325e-01;vvr3(24)= 0.1466388138184849e+00;vvr4(24)= 0.2126547254148320e-01;vvwg(24)= 0.8110770829903342e-02
vvr1(25)= 0.8108302410985486e+00;vvr2(25)= 0.2126547254148325e-01;vvr3(25)= 0.2126547254148325e-01;vvr4(25)= 0.1466388138184849e+00;vvwg(25)= 0.8110770829903342e-02
vvr1(26)= 0.2126547254148325e-01;vvr2(26)= 0.2126547254148325e-01;vvr3(26)= 0.8108302410985486e+00;vvr4(26)= 0.1466388138184849e+00;vvwg(26)= 0.8110770829903342e-02
vvr1(27)= 0.1466388138184849e+00;vvr2(27)= 0.8108302410985486e+00;vvr3(27)= 0.2126547254148325e-01;vvr4(27)= 0.2126547254148325e-01;vvwg(27)= 0.8110770829903342e-02
vvr1(28)= 0.2126547254148325e-01;vvr2(28)= 0.1466388138184849e+00;vvr3(28)= 0.8108302410985486e+00;vvr4(28)= 0.2126547254148314e-01;vvwg(28)= 0.8110770829903342e-02
vvr1(29)= 0.2126547254148325e-01;vvr2(29)= 0.8108302410985486e+00;vvr3(29)= 0.2126547254148325e-01;vvr4(29)= 0.1466388138184849e+00;vvwg(29)= 0.8110770829903342e-02
vvr1(30)= 0.1466388138184849e+00;vvr2(30)= 0.2126547254148325e-01;vvr3(30)= 0.8108302410985486e+00;vvr4(30)= 0.2126547254148325e-01;vvwg(30)= 0.8110770829903342e-02
vvr1(31)= 0.2126547254148325e-01;vvr2(31)= 0.1466388138184849e+00;vvr3(31)= 0.2126547254148325e-01;vvr4(31)= 0.8108302410985485e+00;vvwg(31)= 0.8110770829903342e-02
vvr1(32)= 0.2126547254148325e-01;vvr2(32)= 0.2126547254148325e-01;vvr3(32)= 0.1466388138184849e+00;vvr4(32)= 0.8108302410985486e+00;vvwg(32)= 0.8110770829903342e-02
vvr1(33)= 0.2126547254148325e-01;vvr2(33)= 0.8108302410985486e+00;vvr3(33)= 0.1466388138184849e+00;vvr4(33)= 0.2126547254148320e-01;vvwg(33)= 0.8110770829903342e-02
vvr1(34)= 0.1466388138184849e+00;vvr2(34)= 0.2126547254148325e-01;vvr3(34)= 0.2126547254148325e-01;vvr4(34)= 0.8108302410985486e+00;vvwg(34)= 0.8110770829903342e-02
vvr1(35)= 0.8108302410985486e+00;vvr2(35)= 0.1466388138184849e+00;vvr3(35)= 0.2126547254148325e-01;vvr4(35)= 0.2126547254148320e-01;vvwg(35)= 0.8110770829903342e-02


	
	
case(6,7,8,9)
vvr1( 1)= 0.9551438045408220e+00;vvr2( 1)= 0.1495206515305920e-01;vvr3( 1)= 0.1495206515305920e-01;vvr4( 1)= 0.1495206515305920e-01;vvwg( 1)= 0.1037311233614000e-02
vvr1( 2)= 0.1495206515305920e-01;vvr2( 2)= 0.9551438045408220e+00;vvr3( 2)= 0.1495206515305920e-01;vvr4( 2)= 0.1495206515305920e-01;vvwg( 2)= 0.1037311233614000e-02
vvr1( 3)= 0.1495206515305920e-01;vvr2( 3)= 0.1495206515305920e-01;vvr3( 3)= 0.9551438045408220e+00;vvr4( 3)= 0.1495206515305920e-01;vvwg( 3)= 0.1037311233614000e-02
vvr1( 4)= 0.1495206515305920e-01;vvr2( 4)= 0.1495206515305920e-01;vvr3( 4)= 0.1495206515305920e-01;vvr4( 4)= 0.9551438045408220e+00;vvwg( 4)= 0.1037311233614000e-02
vvr1( 5)= 0.7799760084415400e+00;vvr2( 5)= 0.1518319491659370e+00;vvr3( 5)= 0.3409602119626150e-01;vvr4( 5)= 0.3409602119626150e-01;vvwg( 5)= 0.9601664539948001e-02
vvr1( 6)= 0.1518319491659370e+00;vvr2( 6)= 0.7799760084415400e+00;vvr3( 6)= 0.3409602119626150e-01;vvr4( 6)= 0.3409602119626150e-01;vvwg( 6)= 0.9601664539948001e-02
vvr1( 7)= 0.7799760084415400e+00;vvr2( 7)= 0.3409602119626150e-01;vvr3( 7)= 0.1518319491659370e+00;vvr4( 7)= 0.3409602119626150e-01;vvwg( 7)= 0.9601664539948001e-02
vvr1( 8)= 0.1518319491659370e+00;vvr2( 8)= 0.3409602119626150e-01;vvr3( 8)= 0.7799760084415400e+00;vvr4( 8)= 0.3409602119626150e-01;vvwg( 8)= 0.9601664539948001e-02
vvr1( 9)= 0.7799760084415400e+00;vvr2( 9)= 0.3409602119626150e-01;vvr3( 9)= 0.3409602119626150e-01;vvr4( 9)= 0.1518319491659370e+00;vvwg( 9)= 0.9601664539948001e-02
vvr1(10)= 0.1518319491659370e+00;vvr2(10)= 0.3409602119626150e-01;vvr3(10)= 0.3409602119626150e-01;vvr4(10)= 0.7799760084415400e+00;vvwg(10)= 0.9601664539948001e-02
vvr1(11)= 0.3409602119626150e-01;vvr2(11)= 0.7799760084415400e+00;vvr3(11)= 0.1518319491659370e+00;vvr4(11)= 0.3409602119626150e-01;vvwg(11)= 0.9601664539948001e-02
vvr1(12)= 0.3409602119626150e-01;vvr2(12)= 0.1518319491659370e+00;vvr3(12)= 0.7799760084415400e+00;vvr4(12)= 0.3409602119626150e-01;vvwg(12)= 0.9601664539948001e-02
vvr1(13)= 0.3409602119626150e-01;vvr2(13)= 0.7799760084415400e+00;vvr3(13)= 0.3409602119626150e-01;vvr4(13)= 0.1518319491659370e+00;vvwg(13)= 0.9601664539948001e-02
vvr1(14)= 0.3409602119626150e-01;vvr2(14)= 0.1518319491659370e+00;vvr3(14)= 0.3409602119626150e-01;vvr4(14)= 0.7799760084415400e+00;vvwg(14)= 0.9601664539948001e-02
vvr1(15)= 0.3409602119626150e-01;vvr2(15)= 0.3409602119626150e-01;vvr3(15)= 0.7799760084415400e+00;vvr4(15)= 0.1518319491659370e+00;vvwg(15)= 0.9601664539948001e-02
vvr1(16)= 0.3409602119626150e-01;vvr2(16)= 0.3409602119626150e-01;vvr3(16)= 0.1518319491659370e+00;vvr4(16)= 0.7799760084415400e+00;vvwg(16)= 0.9601664539948001e-02
vvr1(17)= 0.3549340560639790e+00;vvr2(17)= 0.5526556431060170e+00;vvr3(17)= 0.4620515041500170e-01;vvr4(17)= 0.4620515041500170e-01;vvwg(17)= 0.1644939767982320e-01
vvr1(18)= 0.5526556431060170e+00;vvr2(18)= 0.3549340560639790e+00;vvr3(18)= 0.4620515041500170e-01;vvr4(18)= 0.4620515041500170e-01;vvwg(18)= 0.1644939767982320e-01
vvr1(19)= 0.3549340560639790e+00;vvr2(19)= 0.4620515041500170e-01;vvr3(19)= 0.5526556431060170e+00;vvr4(19)= 0.4620515041500170e-01;vvwg(19)= 0.1644939767982320e-01
vvr1(20)= 0.5526556431060170e+00;vvr2(20)= 0.4620515041500170e-01;vvr3(20)= 0.3549340560639790e+00;vvr4(20)= 0.4620515041500170e-01;vvwg(20)= 0.1644939767982320e-01
vvr1(21)= 0.3549340560639790e+00;vvr2(21)= 0.4620515041500170e-01;vvr3(21)= 0.4620515041500170e-01;vvr4(21)= 0.5526556431060170e+00;vvwg(21)= 0.1644939767982320e-01
vvr1(22)= 0.5526556431060170e+00;vvr2(22)= 0.4620515041500170e-01;vvr3(22)= 0.4620515041500170e-01;vvr4(22)= 0.3549340560639790e+00;vvwg(22)= 0.1644939767982320e-01
vvr1(23)= 0.4620515041500170e-01;vvr2(23)= 0.3549340560639790e+00;vvr3(23)= 0.5526556431060170e+00;vvr4(23)= 0.4620515041500170e-01;vvwg(23)= 0.1644939767982320e-01
vvr1(24)= 0.4620515041500170e-01;vvr2(24)= 0.5526556431060170e+00;vvr3(24)= 0.3549340560639790e+00;vvr4(24)= 0.4620515041500170e-01;vvwg(24)= 0.1644939767982320e-01
vvr1(25)= 0.4620515041500170e-01;vvr2(25)= 0.3549340560639790e+00;vvr3(25)= 0.4620515041500170e-01;vvr4(25)= 0.5526556431060170e+00;vvwg(25)= 0.1644939767982320e-01
vvr1(26)= 0.4620515041500170e-01;vvr2(26)= 0.5526556431060170e+00;vvr3(26)= 0.4620515041500170e-01;vvr4(26)= 0.3549340560639790e+00;vvwg(26)= 0.1644939767982320e-01
vvr1(27)= 0.4620515041500170e-01;vvr2(27)= 0.4620515041500170e-01;vvr3(27)= 0.3549340560639790e+00;vvr4(27)= 0.5526556431060170e+00;vvwg(27)= 0.1644939767982320e-01
vvr1(28)= 0.4620515041500170e-01;vvr2(28)= 0.4620515041500170e-01;vvr3(28)= 0.5526556431060170e+00;vvr4(28)= 0.3549340560639790e+00;vvwg(28)= 0.1644939767982320e-01
vvr1(29)= 0.5381043228880020e+00;vvr2(29)= 0.2281904610687610e+00;vvr3(29)= 0.2281904610687610e+00;vvr4(29)= 0.5514754974477500e-02;vvwg(29)= 0.1537477665133100e-01
vvr1(30)= 0.2281904610687610e+00;vvr2(30)= 0.5381043228880020e+00;vvr3(30)= 0.2281904610687610e+00;vvr4(30)= 0.5514754974477500e-02;vvwg(30)= 0.1537477665133100e-01
vvr1(31)= 0.2281904610687610e+00;vvr2(31)= 0.2281904610687610e+00;vvr3(31)= 0.5381043228880020e+00;vvr4(31)= 0.5514754974477500e-02;vvwg(31)= 0.1537477665133100e-01
vvr1(32)= 0.5381043228880020e+00;vvr2(32)= 0.2281904610687610e+00;vvr3(32)= 0.5514754974477500e-02;vvr4(32)= 0.2281904610687610e+00;vvwg(32)= 0.1537477665133100e-01
vvr1(33)= 0.2281904610687610e+00;vvr2(33)= 0.5381043228880020e+00;vvr3(33)= 0.5514754974477500e-02;vvr4(33)= 0.2281904610687610e+00;vvwg(33)= 0.1537477665133100e-01
vvr1(34)= 0.2281904610687610e+00;vvr2(34)= 0.2281904610687610e+00;vvr3(34)= 0.5514754974477500e-02;vvr4(34)= 0.5381043228880020e+00;vvwg(34)= 0.1537477665133100e-01
vvr1(35)= 0.5381043228880020e+00;vvr2(35)= 0.5514754974477500e-02;vvr3(35)= 0.2281904610687610e+00;vvr4(35)= 0.2281904610687610e+00;vvwg(35)= 0.1537477665133100e-01
vvr1(36)= 0.2281904610687610e+00;vvr2(36)= 0.5514754974477500e-02;vvr3(36)= 0.5381043228880020e+00;vvr4(36)= 0.2281904610687610e+00;vvwg(36)= 0.1537477665133100e-01
vvr1(37)= 0.2281904610687610e+00;vvr2(37)= 0.5514754974477500e-02;vvr3(37)= 0.2281904610687610e+00;vvr4(37)= 0.5381043228880020e+00;vvwg(37)= 0.1537477665133100e-01
vvr1(38)= 0.5514754974477500e-02;vvr2(38)= 0.5381043228880020e+00;vvr3(38)= 0.2281904610687610e+00;vvr4(38)= 0.2281904610687610e+00;vvwg(38)= 0.1537477665133100e-01
vvr1(39)= 0.5514754974477500e-02;vvr2(39)= 0.2281904610687610e+00;vvr3(39)= 0.5381043228880020e+00;vvr4(39)= 0.2281904610687610e+00;vvwg(39)= 0.1537477665133100e-01
vvr1(40)= 0.5514754974477500e-02;vvr2(40)= 0.2281904610687610e+00;vvr3(40)= 0.2281904610687610e+00;vvr4(40)= 0.5381043228880020e+00;vvwg(40)= 0.1537477665133100e-01
vvr1(41)= 0.1961837595745600e+00;vvr2(41)= 0.3523052600879940e+00;vvr3(41)= 0.3523052600879940e+00;vvr4(41)= 0.9920572024945300e-01;vvwg(41)= 0.2935201183752300e-01
vvr1(42)= 0.3523052600879940e+00;vvr2(42)= 0.1961837595745600e+00;vvr3(42)= 0.3523052600879940e+00;vvr4(42)= 0.9920572024945300e-01;vvwg(42)= 0.2935201183752300e-01
vvr1(43)= 0.3523052600879940e+00;vvr2(43)= 0.3523052600879940e+00;vvr3(43)= 0.1961837595745600e+00;vvr4(43)= 0.9920572024945300e-01;vvwg(43)= 0.2935201183752300e-01
vvr1(44)= 0.1961837595745600e+00;vvr2(44)= 0.3523052600879940e+00;vvr3(44)= 0.9920572024945300e-01;vvr4(44)= 0.3523052600879940e+00;vvwg(44)= 0.2935201183752300e-01
vvr1(45)= 0.3523052600879940e+00;vvr2(45)= 0.1961837595745600e+00;vvr3(45)= 0.9920572024945300e-01;vvr4(45)= 0.3523052600879940e+00;vvwg(45)= 0.2935201183752300e-01
vvr1(46)= 0.3523052600879940e+00;vvr2(46)= 0.3523052600879940e+00;vvr3(46)= 0.9920572024945300e-01;vvr4(46)= 0.1961837595745600e+00;vvwg(46)= 0.2935201183752300e-01
vvr1(47)= 0.1961837595745600e+00;vvr2(47)= 0.9920572024945300e-01;vvr3(47)= 0.3523052600879940e+00;vvr4(47)= 0.3523052600879940e+00;vvwg(47)= 0.2935201183752300e-01
vvr1(48)= 0.3523052600879940e+00;vvr2(48)= 0.9920572024945300e-01;vvr3(48)= 0.1961837595745600e+00;vvr4(48)= 0.3523052600879940e+00;vvwg(48)= 0.2935201183752300e-01
vvr1(49)= 0.3523052600879940e+00;vvr2(49)= 0.9920572024945300e-01;vvr3(49)= 0.3523052600879940e+00;vvr4(49)= 0.1961837595745600e+00;vvwg(49)= 0.2935201183752300e-01
vvr1(50)= 0.9920572024945300e-01;vvr2(50)= 0.1961837595745600e+00;vvr3(50)= 0.3523052600879940e+00;vvr4(50)= 0.3523052600879940e+00;vvwg(50)= 0.2935201183752300e-01
vvr1(51)= 0.9920572024945300e-01;vvr2(51)= 0.3523052600879940e+00;vvr3(51)= 0.1961837595745600e+00;vvr4(51)= 0.3523052600879940e+00;vvwg(51)= 0.2935201183752300e-01
vvr1(52)= 0.9920572024945300e-01;vvr2(52)= 0.3523052600879940e+00;vvr3(52)= 0.3523052600879940e+00;vvr4(52)= 0.1961837595745600e+00;vvwg(52)= 0.2935201183752300e-01
vvr1(53)= 0.5965649956210169e+00;vvr2(53)= 0.1344783347929940e+00;vvr3(53)= 0.1344783347929940e+00;vvr4(53)= 0.1344783347929940e+00;vvwg(53)= 0.3662913664051080e-01
vvr1(54)= 0.1344783347929940e+00;vvr2(54)= 0.5965649956210169e+00;vvr3(54)= 0.1344783347929940e+00;vvr4(54)= 0.1344783347929940e+00;vvwg(54)= 0.3662913664051080e-01
vvr1(55)= 0.1344783347929940e+00;vvr2(55)= 0.1344783347929940e+00;vvr3(55)= 0.5965649956210169e+00;vvr4(55)= 0.1344783347929940e+00;vvwg(55)= 0.3662913664051080e-01
vvr1(56)= 0.1344783347929940e+00;vvr2(56)= 0.1344783347929940e+00;vvr3(56)= 0.1344783347929940e+00;vvr4(56)= 0.5965649956210169e+00;vvwg(56)= 0.3662913664051080e-01
 	
	
	
	
	
	
			
end select



		do kk=1,qp_tetra
			wequa3d(kk)=vvwg(kk)
			qpoints(:,kk)=(vvr1(kk)*vext(1,:))+(vvr2(kk)*vext(2,:))+(vvr3(kk)*vext(3,:))+(vvr4(kk)*vext(4,:))
			
		
		end do


end subroutine quadraturetetra

subroutine quadratureprism(n,igqrules,vext,qpoints,wequa3d)
 !> @brief
!> this subroutine computes the quadrature points and weights for a prism
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints),intent(inout)::qpoints
real,dimension(1:numberofpoints),intent(inout)::wequa3d
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1)::deta
real,dimension(1:8)::vvnxi
real,dimension(1:alls)::vvwg,vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
real::r,s,tx,a,b,c,d,e,f
real::a1,b1,c1,d1,e1,f1,sumwe
integer::kk,j,ii,ij,ik,count1




 
 wequa3d=0.0d0
  qpoints=0.0d0
sumwe=0.0d0

select case(igqrules)
 
  case(1)
		vvr1(1)=0.666666666666667/2.0d0 ;vvr2(1)=0.666666666666667/2.0d0 ;vvr3(1)=0.0d0
		
		  vvwg(1)=2.0d0

 
  	


 

 
 case(2)
		vvr1(1)=0.666666666666667 ;vvr2(1)=0.166666666666667 ;vvr3(1)=0.5773502691896257;vvwg(1)=0.33333333333333333333*1.0d0
		vvr1(2)=0.166666666666667 ;vvr2(2)=0.666666666666667 ;vvr3(2)=0.5773502691896257;vvwg(2)=0.33333333333333333333*1.0d0
		vvr1(3)=0.166666666666667 ;vvr2(3)=0.166666666666667 ;vvr3(3)=0.5773502691896257;vvwg(3)=0.33333333333333333333*1.0d0
		vvr1(4)=0.666666666666667 ;vvr2(4)=0.166666666666667 ;vvr3(4)=-0.5773502691896257;vvwg(4)=0.33333333333333333333*1.0d0
		vvr1(5)=0.166666666666667 ;vvr2(5)=0.666666666666667 ;vvr3(5)=-0.5773502691896257;vvwg(5)=0.33333333333333333333*1.0d0
		vvr1(6)=0.166666666666667 ;vvr2(6)=0.166666666666667 ;vvr3(6)=-0.5773502691896257;vvwg(6)=0.33333333333333333333*1.0d0
 

 
  	


  	
 case(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
 

		vvr1(1)=0.816847572980440 ;vvr2(1)=0.091576213509780 ;vvr3(1)=0.0 ;vvwg(1)=0.109951743655333*0.8888888888888888
		vvr1(2)=0.091576213509780 ;vvr2(2)=0.816847572980440 ;vvr3(2)=0.0 ;vvwg(2)=0.109951743655333*0.8888888888888888
		vvr1(3)=0.091576213509780 ;vvr2(3)=0.091576213509780 ;vvr3(3)=0.0 ;vvwg(3)=0.109951743655333*0.8888888888888888
		vvr1(4)=0.445948490915964 ;vvr2(4)=0.445948490915964 ;vvr3(4)=0.0 ;vvwg(4)=0.223381589678000*0.8888888888888888
		vvr1(5)=0.445948490915964 ;vvr2(5)=0.108103018168071 ;vvr3(5)=0.0 ;vvwg(5)=0.223381589678000*0.8888888888888888
		vvr1(6)=0.108103018168071 ;vvr2(6)=0.445948490915964 ;vvr3(6)=0.0 ;vvwg(6)=0.223381589678000*0.8888888888888888
		vvr1(7)=0.816847572980440 ;vvr2(7)=0.091576213509780 ;vvr3(7)=-0.7745966692414834 ;vvwg(7)=0.109951743655333*0.5555555555555556
		vvr1(8)=0.091576213509780 ;vvr2(8)=0.816847572980440 ;vvr3(8)=-0.7745966692414834 ;vvwg(8)=0.109951743655333*0.5555555555555556
		vvr1(9)=0.091576213509780 ;vvr2(9)=0.091576213509780 ;vvr3(9)=-0.7745966692414834 ;vvwg(9)=0.109951743655333*0.5555555555555556
		vvr1(10)=0.445948490915964 ;vvr2(10)=0.445948490915964 ;vvr3(10)=-0.7745966692414834 ;vvwg(10)=0.223381589678000*0.5555555555555556
		vvr1(11)=0.445948490915964 ;vvr2(11)=0.108103018168071 ;vvr3(11)=-0.7745966692414834 ;vvwg(11)=0.223381589678000*0.5555555555555556
		vvr1(12)=0.108103018168071 ;vvr2(12)=0.445948490915964 ;vvr3(12)=-0.7745966692414834 ;vvwg(12)=0.223381589678000*0.5555555555555556
		vvr1(13)=0.816847572980440 ;vvr2(13)=0.091576213509780 ;vvr3(13)=0.7745966692414834 ;vvwg(13)=0.109951743655333*0.5555555555555556
		vvr1(14)=0.091576213509780 ;vvr2(14)=0.816847572980440 ;vvr3(14)=0.7745966692414834 ;vvwg(14)=0.109951743655333*0.5555555555555556
		vvr1(15)=0.091576213509780 ;vvr2(15)=0.091576213509780 ;vvr3(15)=0.7745966692414834 ;vvwg(15)=0.109951743655333*0.5555555555555556
		vvr1(16)=0.445948490915964 ;vvr2(16)=0.445948490915964 ;vvr3(16)=0.7745966692414834 ;vvwg(16)=0.223381589678000*0.5555555555555556
		vvr1(17)=0.445948490915964 ;vvr2(17)=0.108103018168071 ;vvr3(17)=0.7745966692414834 ;vvwg(17)=0.223381589678000*0.5555555555555556
		vvr1(18)=0.108103018168071 ;vvr2(18)=0.445948490915964 ;vvr3(18)=0.7745966692414834 ;vvwg(18)=0.223381589678000*0.5555555555555556











			
 case(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
 
		vvr1(1)=0.816847572980440 ;vvr2(1)=0.091576213509780 ;vvr3(1)=-0.3399810435848563 ;vvwg(1)=0.109951743655333*0.6521451548625461
		vvr1(2)=0.091576213509780 ;vvr2(2)=0.816847572980440 ;vvr3(2)=-0.3399810435848563 ;vvwg(2)=0.109951743655333*0.6521451548625461
		vvr1(3)=0.091576213509780 ;vvr2(3)=0.091576213509780 ;vvr3(3)=-0.3399810435848563 ;vvwg(3)=0.109951743655333*0.6521451548625461
		vvr1(4)=0.445948490915964 ;vvr2(4)=0.445948490915964 ;vvr3(4)=-0.3399810435848563 ;vvwg(4)=0.223381589678000*0.6521451548625461
		vvr1(5)=0.445948490915964 ;vvr2(5)=0.108103018168071 ;vvr3(5)=-0.3399810435848563 ;vvwg(5)=0.223381589678000*0.6521451548625461
		vvr1(6)=0.108103018168071 ;vvr2(6)=0.445948490915964 ;vvr3(6)=-0.3399810435848563 ;vvwg(6)=0.223381589678000*0.6521451548625461
		vvr1(7)=0.816847572980440 ;vvr2(7)=0.091576213509780 ;vvr3(7)=0.3399810435848563 ;vvwg(7)=0.109951743655333*0.6521451548625461
		vvr1(8)=0.091576213509780 ;vvr2(8)=0.816847572980440 ;vvr3(8)=0.3399810435848563 ;vvwg(8)=0.109951743655333*0.6521451548625461
		vvr1(9)=0.091576213509780 ;vvr2(9)=0.091576213509780 ;vvr3(9)=0.3399810435848563 ;vvwg(9)=0.109951743655333*0.6521451548625461
		vvr1(10)=0.445948490915964 ;vvr2(10)=0.445948490915964 ;vvr3(10)=0.3399810435848563 ;vvwg(10)=0.223381589678000*0.6521451548625461
		vvr1(11)=0.445948490915964 ;vvr2(11)=0.108103018168071 ;vvr3(11)=0.3399810435848563 ;vvwg(11)=0.223381589678000*0.6521451548625461
		vvr1(12)=0.108103018168071 ;vvr2(12)=0.445948490915964 ;vvr3(12)=0.3399810435848563 ;vvwg(12)=0.223381589678000*0.6521451548625461
		vvr1(13)=0.816847572980440 ;vvr2(13)=0.091576213509780 ;vvr3(13)=-0.8611363115940526 ;vvwg(13)=0.109951743655333*0.3478548451374538
		vvr1(14)=0.091576213509780 ;vvr2(14)=0.816847572980440 ;vvr3(14)=-0.8611363115940526 ;vvwg(14)=0.109951743655333*0.3478548451374538
		vvr1(15)=0.091576213509780 ;vvr2(15)=0.091576213509780 ;vvr3(15)=-0.8611363115940526 ;vvwg(15)=0.109951743655333*0.3478548451374538
		vvr1(16)=0.445948490915964 ;vvr2(16)=0.445948490915964 ;vvr3(16)=-0.8611363115940526 ;vvwg(16)=0.223381589678000*0.3478548451374538
		vvr1(17)=0.445948490915964 ;vvr2(17)=0.108103018168071 ;vvr3(17)=-0.8611363115940526 ;vvwg(17)=0.223381589678000*0.3478548451374538
		vvr1(18)=0.108103018168071 ;vvr2(18)=0.445948490915964 ;vvr3(18)=-0.8611363115940526 ;vvwg(18)=0.223381589678000*0.3478548451374538
		vvr1(19)=0.816847572980440 ;vvr2(19)=0.091576213509780 ;vvr3(19)=0.8611363115940526 ;vvwg(19)=0.109951743655333*0.3478548451374538
		vvr1(20)=0.091576213509780 ;vvr2(20)=0.816847572980440 ;vvr3(20)=0.8611363115940526 ;vvwg(20)=0.109951743655333*0.3478548451374538
		vvr1(21)=0.091576213509780 ;vvr2(21)=0.091576213509780 ;vvr3(21)=0.8611363115940526 ;vvwg(21)=0.109951743655333*0.3478548451374538
		vvr1(22)=0.445948490915964 ;vvr2(22)=0.445948490915964 ;vvr3(22)=0.8611363115940526 ;vvwg(22)=0.223381589678000*0.3478548451374538
		vvr1(23)=0.445948490915964 ;vvr2(23)=0.108103018168071 ;vvr3(23)=0.8611363115940526 ;vvwg(23)=0.223381589678000*0.3478548451374538
		vvr1(24)=0.108103018168071 ;vvr2(24)=0.445948490915964 ;vvr3(24)=0.8611363115940526 ;vvwg(24)=0.223381589678000*0.3478548451374538


			
case(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  e1=0.2369268850561891
				vvr1(1)=0.888871894660413 ;vvr2(1)=0.055564052669793 ;vvr3(1)=0.0 ;vvwg(1)=0.041955512996649*0.5688888888888889
				vvr1(2)=0.055564052669793 ;vvr2(2)=0.888871894660413 ;vvr3(2)=0.0 ;vvwg(2)=0.041955512996649*0.5688888888888889
				vvr1(3)=0.055564052669793 ;vvr2(3)=0.055564052669793 ;vvr3(3)=0.0 ;vvwg(3)=0.041955512996649*0.5688888888888889
				vvr1(4)=0.295533711735893 ;vvr2(4)=0.634210747745723 ;vvr3(4)=0.0 ;vvwg(4)=0.112098412070887*0.5688888888888889
				vvr1(5)=0.295533711735893 ;vvr2(5)=0.070255540518384 ;vvr3(5)=0.0 ;vvwg(5)=0.112098412070887*0.5688888888888889
				vvr1(6)=0.070255540518384 ;vvr2(6)=0.295533711735893 ;vvr3(6)=0.0 ;vvwg(6)=0.112098412070887*0.5688888888888889
				vvr1(7)=0.634210747745723 ;vvr2(7)=0.295533711735893 ;vvr3(7)=0.0 ;vvwg(7)=0.112098412070887*0.5688888888888889
				vvr1(8)=0.634210747745723 ;vvr2(8)=0.070255540518384 ;vvr3(8)=0.0 ;vvwg(8)=0.112098412070887*0.5688888888888889
				vvr1(9)=0.070255540518384 ;vvr2(9)=0.634210747745723 ;vvr3(9)=0.0 ;vvwg(9)=0.112098412070887*0.5688888888888889
				vvr1(10)=0.333333333333333 ;vvr2(10)=0.333333333333333 ;vvr3(10)=0.0 ;vvwg(10)=0.201542988584730*0.5688888888888889
				vvr1(11)=0.888871894660413 ;vvr2(11)=0.055564052669793 ;vvr3(11)=-0.5384693101056831 ;vvwg(11)=0.041955512996649*0.4786286704993665
				vvr1(12)=0.055564052669793 ;vvr2(12)=0.888871894660413 ;vvr3(12)=-0.5384693101056831 ;vvwg(12)=0.041955512996649*0.4786286704993665
				vvr1(13)=0.055564052669793 ;vvr2(13)=0.055564052669793 ;vvr3(13)=-0.5384693101056831 ;vvwg(13)=0.041955512996649*0.4786286704993665
				vvr1(14)=0.295533711735893 ;vvr2(14)=0.634210747745723 ;vvr3(14)=-0.5384693101056831 ;vvwg(14)=0.112098412070887*0.4786286704993665
				vvr1(15)=0.295533711735893 ;vvr2(15)=0.070255540518384 ;vvr3(15)=-0.5384693101056831 ;vvwg(15)=0.112098412070887*0.4786286704993665
				vvr1(16)=0.070255540518384 ;vvr2(16)=0.295533711735893 ;vvr3(16)=-0.5384693101056831 ;vvwg(16)=0.112098412070887*0.4786286704993665
				vvr1(17)=0.634210747745723 ;vvr2(17)=0.295533711735893 ;vvr3(17)=-0.53846931010568314 ;vvwg(17)=0.112098412070887*0.4786286704993665
				vvr1(18)=0.634210747745723 ;vvr2(18)=0.070255540518384 ;vvr3(18)=-0.5384693101056831 ;vvwg(18)=0.112098412070887*0.4786286704993665
				vvr1(19)=0.070255540518384 ;vvr2(19)=0.634210747745723 ;vvr3(19)=-0.5384693101056831 ;vvwg(19)=0.112098412070887*0.4786286704993665
				vvr1(20)=0.333333333333333 ;vvr2(20)=0.333333333333333 ;vvr3(20)=-0.5384693101056831 ;vvwg(20)=0.201542988584730*0.4786286704993665
				vvr1(21)=0.888871894660413 ;vvr2(21)=0.055564052669793 ;vvr3(21)=0.5384693101056831 ;vvwg(21)=0.041955512996649*0.4786286704993665
				vvr1(22)=0.055564052669793 ;vvr2(22)=0.888871894660413 ;vvr3(22)=0.5384693101056831 ;vvwg(22)=0.041955512996649*0.4786286704993665
				vvr1(23)=0.055564052669793 ;vvr2(23)=0.055564052669793 ;vvr3(23)=0.5384693101056831 ;vvwg(23)=0.041955512996649*0.4786286704993665
				vvr1(24)=0.295533711735893 ;vvr2(24)=0.634210747745723 ;vvr3(24)=0.5384693101056831;vvwg(24)=0.112098412070887*0.4786286704993665
				vvr1(25)=0.295533711735893 ;vvr2(25)=0.070255540518384 ;vvr3(25)=0.5384693101056831 ;vvwg(25)=0.112098412070887*0.4786286704993665
				vvr1(26)=0.070255540518384 ;vvr2(26)=0.295533711735893 ;vvr3(26)=0.5384693101056831 ;vvwg(26)=0.112098412070887*0.4786286704993665
				vvr1(27)=0.634210747745723 ;vvr2(27)=0.295533711735893 ;vvr3(27)=0.5384693101056831 ;vvwg(27)=0.112098412070887*0.4786286704993665
				vvr1(28)=0.634210747745723 ;vvr2(28)=0.070255540518384 ;vvr3(28)=0.5384693101056831 ;vvwg(28)=0.112098412070887*0.4786286704993665
				vvr1(29)=0.070255540518384 ;vvr2(29)=0.634210747745723 ;vvr3(29)=0.5384693101056831 ;vvwg(29)=0.112098412070887*0.4786286704993665
				vvr1(30)=0.333333333333333 ;vvr2(30)=0.333333333333333 ;vvr3(30)=0.5384693101056831 ;vvwg(30)=0.201542988584730*0.4786286704993665
				vvr1(31)=0.888871894660413 ;vvr2(31)=0.055564052669793 ;vvr3(31)=-0.9061798459386640 ;vvwg(31)=0.041955512996649*0.2369268850561891
				vvr1(32)=0.055564052669793 ;vvr2(32)=0.888871894660413 ;vvr3(32)=-0.9061798459386640 ;vvwg(32)=0.041955512996649*0.2369268850561891
				vvr1(33)=0.055564052669793 ;vvr2(33)=0.055564052669793 ;vvr3(33)=-0.9061798459386640 ;vvwg(33)=0.041955512996649*0.2369268850561891
				vvr1(34)=0.295533711735893 ;vvr2(34)=0.634210747745723 ;vvr3(34)=-0.9061798459386640 ;vvwg(34)=0.112098412070887*0.2369268850561891
				vvr1(35)=0.295533711735893 ;vvr2(35)=0.070255540518384 ;vvr3(35)=-0.9061798459386640 ;vvwg(35)=0.112098412070887*0.2369268850561891
				vvr1(36)=0.070255540518384 ;vvr2(36)=0.295533711735893 ;vvr3(36)=-0.9061798459386640 ;vvwg(36)=0.112098412070887*0.2369268850561891
				vvr1(37)=0.634210747745723 ;vvr2(37)=0.295533711735893 ;vvr3(37)=-0.9061798459386640 ;vvwg(37)=0.112098412070887*0.2369268850561891
				vvr1(38)=0.634210747745723 ;vvr2(38)=0.070255540518384 ;vvr3(38)=-0.9061798459386640 ;vvwg(38)=0.112098412070887*0.2369268850561891
				vvr1(39)=0.070255540518384 ;vvr2(39)=0.634210747745723 ;vvr3(39)=-0.9061798459386640 ;vvwg(39)=0.112098412070887*0.2369268850561891
				vvr1(40)=0.333333333333333 ;vvr2(40)=0.333333333333333 ;vvr3(40)=-0.9061798459386640 ;vvwg(40)=0.201542988584730*0.2369268850561891
				vvr1(41)=0.888871894660413 ;vvr2(41)=0.055564052669793 ;vvr3(41)=0.9061798459386640 ;vvwg(41)=0.041955512996649*0.2369268850561891
				vvr1(42)=0.055564052669793 ;vvr2(42)=0.888871894660413 ;vvr3(42)=0.9061798459386640;vvwg(42)=0.041955512996649*0.2369268850561891
				vvr1(43)=0.055564052669793 ;vvr2(43)=0.055564052669793 ;vvr3(43)=0.9061798459386640 ;vvwg(43)=0.041955512996649*0.2369268850561891
				vvr1(44)=0.295533711735893 ;vvr2(44)=0.634210747745723 ;vvr3(44)=0.9061798459386640;vvwg(44)=0.112098412070887*0.2369268850561891
				vvr1(45)=0.295533711735893 ;vvr2(45)=0.070255540518384 ;vvr3(45)=0.9061798459386640 ;vvwg(45)=0.112098412070887*0.2369268850561891
				vvr1(46)=0.070255540518384 ;vvr2(46)=0.295533711735893 ;vvr3(46)=0.9061798459386640 ;vvwg(46)=0.112098412070887*0.2369268850561891
				vvr1(47)=0.634210747745723 ;vvr2(47)=0.295533711735893 ;vvr3(47)=0.9061798459386640 ;vvwg(47)=0.112098412070887*0.2369268850561891
				vvr1(48)=0.634210747745723 ;vvr2(48)=0.070255540518384 ;vvr3(48)=0.9061798459386640 ;vvwg(48)=0.112098412070887*0.2369268850561891
				vvr1(49)=0.070255540518384 ;vvr2(49)=0.634210747745723 ;vvr3(49)=0.9061798459386640 ;vvwg(49)=0.112098412070887*0.2369268850561891
				vvr1(50)=0.333333333333333 ;vvr2(50)=0.333333333333333 ;vvr3(50)=0.9061798459386640 ;vvwg(50)=0.201542988584730*0.2369268850561891
	
			
case(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  f=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  f1=0.1713244923791704
 
			      vvr1(1)=0.888871894660413 ;vvr2(1)=0.055564052669793 ;vvr3(1)=0.6612093864662645 ;vvwg(1)=0.041955512996649*0.3607615730481386
				vvr1(2)=0.055564052669793 ;vvr2(2)=0.888871894660413 ;vvr3(2)=0.6612093864662645 ;vvwg(2)=0.041955512996649*0.3607615730481386
				vvr1(3)=0.055564052669793 ;vvr2(3)=0.055564052669793 ;vvr3(3)=0.6612093864662645 ;vvwg(3)=0.041955512996649*0.3607615730481386
				vvr1(4)=0.295533711735893 ;vvr2(4)=0.634210747745723 ;vvr3(4)=0.6612093864662645 ;vvwg(4)=0.112098412070887*0.3607615730481386
				vvr1(5)=0.295533711735893 ;vvr2(5)=0.070255540518384 ;vvr3(5)=0.6612093864662645 ;vvwg(5)=0.112098412070887*0.3607615730481386
				vvr1(6)=0.070255540518384 ;vvr2(6)=0.295533711735893 ;vvr3(6)=0.6612093864662645 ;vvwg(6)=0.112098412070887*0.3607615730481386
				vvr1(7)=0.634210747745723 ;vvr2(7)=0.295533711735893 ;vvr3(7)=0.6612093864662645 ;vvwg(7)=0.112098412070887*0.3607615730481386
				vvr1(8)=0.634210747745723 ;vvr2(8)=0.070255540518384 ;vvr3(8)=0.6612093864662645 ;vvwg(8)=0.112098412070887*0.3607615730481386
				vvr1(9)=0.070255540518384 ;vvr2(9)=0.634210747745723 ;vvr3(9)=0.6612093864662645 ;vvwg(9)=0.112098412070887*0.3607615730481386
				vvr1(10)=0.333333333333333 ;vvr2(10)=0.333333333333333 ;vvr3(10)=0.6612093864662645 ;vvwg(10)=0.201542988584730*0.3607615730481386
				vvr1(11)=0.888871894660413 ;vvr2(11)=0.055564052669793 ;vvr3(11)=-0.6612093864662645 ;vvwg(11)=0.041955512996649*0.3607615730481386
				vvr1(12)=0.055564052669793 ;vvr2(12)=0.888871894660413 ;vvr3(12)=-0.6612093864662645 ;vvwg(12)=0.041955512996649*0.3607615730481386
				vvr1(13)=0.055564052669793 ;vvr2(13)=0.055564052669793 ;vvr3(13)=-0.6612093864662645 ;vvwg(13)=0.041955512996649*0.3607615730481386
				vvr1(14)=0.295533711735893 ;vvr2(14)=0.634210747745723 ;vvr3(14)=-0.6612093864662645 ;vvwg(14)=0.112098412070887*0.3607615730481386
				vvr1(15)=0.295533711735893 ;vvr2(15)=0.070255540518384 ;vvr3(15)=-0.6612093864662645 ;vvwg(15)=0.112098412070887*0.3607615730481386
				vvr1(16)=0.070255540518384 ;vvr2(16)=0.295533711735893 ;vvr3(16)=-0.6612093864662645;vvwg(16)=0.112098412070887*0.3607615730481386
				vvr1(17)=0.634210747745723 ;vvr2(17)=0.295533711735893 ;vvr3(17)=-0.6612093864662645 ;vvwg(17)=0.112098412070887*0.3607615730481386
				vvr1(18)=0.634210747745723 ;vvr2(18)=0.070255540518384 ;vvr3(18)=-0.6612093864662645;vvwg(18)=0.112098412070887*0.3607615730481386
				vvr1(19)=0.070255540518384 ;vvr2(19)=0.634210747745723 ;vvr3(19)=-0.6612093864662645 ;vvwg(19)=0.112098412070887*0.3607615730481386
				vvr1(20)=0.333333333333333 ;vvr2(20)=0.333333333333333 ;vvr3(20)=-0.6612093864662645 ;vvwg(20)=0.201542988584730*0.3607615730481386
				vvr1(21)=0.888871894660413 ;vvr2(21)=0.055564052669793 ;vvr3(21)=-0.2386191860831969 ;vvwg(21)=0.041955512996649*0.4679139345726910
				vvr1(22)=0.055564052669793 ;vvr2(22)=0.888871894660413 ;vvr3(22)=-0.2386191860831969 ;vvwg(22)=0.041955512996649*0.4679139345726910
				vvr1(23)=0.055564052669793 ;vvr2(23)=0.055564052669793 ;vvr3(23)=-0.2386191860831969;vvwg(23)=0.041955512996649*0.4679139345726910
				vvr1(24)=0.295533711735893 ;vvr2(24)=0.634210747745723 ;vvr3(24)=-0.2386191860831969;vvwg(24)=0.112098412070887*0.4679139345726910
				vvr1(25)=0.295533711735893 ;vvr2(25)=0.070255540518384 ;vvr3(25)=-0.2386191860831969 ;vvwg(25)=0.112098412070887*0.4679139345726910
				vvr1(26)=0.070255540518384 ;vvr2(26)=0.295533711735893 ;vvr3(26)=-0.2386191860831969 ;vvwg(26)=0.112098412070887*0.4679139345726910
				vvr1(27)=0.634210747745723 ;vvr2(27)=0.295533711735893 ;vvr3(27)=-0.2386191860831969 ;vvwg(27)=0.112098412070887*0.4679139345726910
				vvr1(28)=0.634210747745723 ;vvr2(28)=0.070255540518384 ;vvr3(28)=-0.2386191860831969 ;vvwg(28)=0.112098412070887*0.4679139345726910
				vvr1(29)=0.070255540518384 ;vvr2(29)=0.634210747745723 ;vvr3(29)=-0.2386191860831969 ;vvwg(29)=0.112098412070887*0.4679139345726910
				vvr1(30)=0.333333333333333 ;vvr2(30)=0.333333333333333 ;vvr3(30)=-0.2386191860831969 ;vvwg(30)=0.201542988584730*0.4679139345726910
				vvr1(31)=0.888871894660413 ;vvr2(31)=0.055564052669793 ;vvr3(31)=0.2386191860831969 ;vvwg(31)=0.041955512996649*0.4679139345726910
				vvr1(32)=0.055564052669793 ;vvr2(32)=0.888871894660413 ;vvr3(32)=0.2386191860831969 ;vvwg(32)=0.041955512996649*0.4679139345726910
				vvr1(33)=0.055564052669793 ;vvr2(33)=0.055564052669793 ;vvr3(33)=0.2386191860831969 ;vvwg(33)=0.041955512996649*0.4679139345726910
				vvr1(34)=0.295533711735893 ;vvr2(34)=0.634210747745723 ;vvr3(34)=0.2386191860831969 ;vvwg(34)=0.112098412070887*0.4679139345726910
				vvr1(35)=0.295533711735893 ;vvr2(35)=0.070255540518384 ;vvr3(35)=0.2386191860831969 ;vvwg(35)=0.112098412070887*0.4679139345726910
				vvr1(36)=0.070255540518384 ;vvr2(36)=0.295533711735893 ;vvr3(36)=0.2386191860831969 ;vvwg(36)=0.112098412070887*0.4679139345726910
				vvr1(37)=0.634210747745723 ;vvr2(37)=0.295533711735893 ;vvr3(37)=0.2386191860831969 ;vvwg(37)=0.112098412070887*0.4679139345726910
				vvr1(38)=0.634210747745723 ;vvr2(38)=0.070255540518384 ;vvr3(38)=0.2386191860831969 ;vvwg(38)=0.112098412070887*0.4679139345726910
				vvr1(39)=0.070255540518384 ;vvr2(39)=0.634210747745723 ;vvr3(39)=0.2386191860831969 ;vvwg(39)=0.112098412070887*0.4679139345726910
				vvr1(40)=0.333333333333333 ;vvr2(40)=0.333333333333333 ;vvr3(40)=0.2386191860831969 ;vvwg(40)=0.201542988584730*0.4679139345726910
				vvr1(41)=0.888871894660413 ;vvr2(41)=0.055564052669793 ;vvr3(41)=-0.9324695142031521 ;vvwg(41)=0.041955512996649*0.1713244923791704
				vvr1(42)=0.055564052669793 ;vvr2(42)=0.888871894660413 ;vvr3(42)=-0.9324695142031521;vvwg(42)=0.041955512996649*0.1713244923791704
				vvr1(43)=0.055564052669793 ;vvr2(43)=0.055564052669793 ;vvr3(43)=-0.9324695142031521 ;vvwg(43)=0.041955512996649*0.1713244923791704
				vvr1(44)=0.295533711735893 ;vvr2(44)=0.634210747745723 ;vvr3(44)=-0.9324695142031521;vvwg(44)=0.112098412070887*0.1713244923791704
				vvr1(45)=0.295533711735893 ;vvr2(45)=0.070255540518384 ;vvr3(45)=-0.9324695142031521 ;vvwg(45)=0.112098412070887*0.1713244923791704
				vvr1(46)=0.070255540518384 ;vvr2(46)=0.295533711735893 ;vvr3(46)=-0.9324695142031521 ;vvwg(46)=0.112098412070887*0.1713244923791704
				vvr1(47)=0.634210747745723 ;vvr2(47)=0.295533711735893 ;vvr3(47)=-0.9324695142031521 ;vvwg(47)=0.112098412070887*0.1713244923791704
				vvr1(48)=0.634210747745723 ;vvr2(48)=0.070255540518384 ;vvr3(48)=-0.9324695142031521 ;vvwg(48)=0.112098412070887*0.1713244923791704
				vvr1(49)=0.070255540518384 ;vvr2(49)=0.634210747745723 ;vvr3(49)=-0.9324695142031521 ;vvwg(49)=0.112098412070887*0.1713244923791704
				vvr1(50)=0.333333333333333 ;vvr2(50)=0.333333333333333 ;vvr3(50)=-0.9324695142031521 ;vvwg(50)=0.201542988584730*0.1713244923791704
				vvr1(51)=0.888871894660413 ;vvr2(51)=0.055564052669793 ;vvr3(51)=0.9324695142031521 ;vvwg(51)=0.041955512996649*0.1713244923791704
				vvr1(52)=0.055564052669793 ;vvr2(52)=0.888871894660413 ;vvr3(52)=0.9324695142031521;vvwg(52)=0.041955512996649*0.1713244923791704
				vvr1(53)=0.055564052669793 ;vvr2(53)=0.055564052669793 ;vvr3(53)=0.9324695142031521 ;vvwg(53)=0.041955512996649*0.1713244923791704
				vvr1(54)=0.295533711735893 ;vvr2(54)=0.634210747745723 ;vvr3(54)=0.9324695142031521;vvwg(54)=0.112098412070887*0.1713244923791704
				vvr1(55)=0.295533711735893 ;vvr2(55)=0.070255540518384 ;vvr3(55)=0.9324695142031521 ;vvwg(55)=0.112098412070887*0.1713244923791704
				vvr1(56)=0.070255540518384 ;vvr2(56)=0.295533711735893 ;vvr3(56)=0.9324695142031521 ;vvwg(56)=0.112098412070887*0.1713244923791704
				vvr1(57)=0.634210747745723 ;vvr2(57)=0.295533711735893 ;vvr3(57)=0.9324695142031521 ;vvwg(57)=0.112098412070887*0.1713244923791704
				vvr1(58)=0.634210747745723 ;vvr2(58)=0.070255540518384 ;vvr3(58)=0.9324695142031521 ;vvwg(58)=0.112098412070887*0.1713244923791704
				vvr1(59)=0.070255540518384 ;vvr2(59)=0.634210747745723 ;vvr3(59)=0.9324695142031521 ;vvwg(59)=0.112098412070887*0.1713244923791704
				vvr1(60)=0.333333333333333 ;vvr2(60)=0.333333333333333 ;vvr3(60)=0.9324695142031521 ;vvwg(60)=0.201542988584730*0.1713244923791704




	
			
end select
		qpoints(:,:)=0.0d0
! 		
! 		  wequa3d(:)=vvwg(:)*0.5d0
		do kk=1,qp_prism
			wequa3d(kk)=vvwg(kk)*0.5d0
			
			r=vvr1(kk); s=vvr2(kk); tx=vvr3(kk)
			vvnxi(1)=(0.5d0)*r*(1.0d0-tx)
			vvnxi(2)=(0.5d0)*(s)*(1.0d0-tx)
			vvnxi(3)=(0.5d0)*(1.0-r-s)*(1.0d0-tx)
			vvnxi(4)=(0.5d0)*r*(1.0d0+tx)
			vvnxi(5)=(0.5d0)*(s)*(1.0d0+tx)
			vvnxi(6)=(0.5d0)*(1.0-r-s)*(1.0d0+tx)
			
			do j=1,6
			qpoints(:,kk)=qpoints(:,kk)+(vvnxi(j)*vext(j,:))
			end do
		    
! 			

		end do
! 		



end subroutine quadratureprism


subroutine quadraturepyra(n,igqrules,vext,qpoints,wequa3d)
 !> @brief
!> this subroutine computes the quadrature points and weights for a pyramid
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints),intent(inout)::qpoints
real,dimension(1:numberofpoints),intent(inout)::wequa3d
real::r,s,tx,a,b,c,d,e,f,g
real::a1,b1,c1,d1,e1,f1,sumwe
integer::kk,j,ii,ij,ik,count1
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1)::deta
real,dimension(1:8)::vvnxi
real,dimension(1:alls)::vvwg,vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz




 
 wequa3d=0.0d0
  qpoints=0.0d0
sumwe=0.0d0

select case(igqrules)
 
  case(1)
		vvr1(1)=0.0d0 ;vvr2(1)=0.0d0 ;vvr3(1)=-0.5d0
		
		  vvwg(1)=8.0d0

 
  	


 

 
 case(2)
		vvr1(1)=-0.584237394672177188;vvr2(1)=-0.58423739467217718 ;vvr3(1)=-0.6666666666666666;vvwg(1)=0.81
		vvr1(2)=0.58423739467217718 ;vvr2(2)=-0.58423739467217718 ;vvr3(2)=-0.6666666666666666;vvwg(2)=0.81
		vvr1(3)=0.58423739467217718 ;vvr2(3)=0.58423739467217718 ;vvr3(3)=-0.6666666666666666;vvwg(3)=0.81
		vvr1(4)=-0.58423739467217718 ;vvr2(4)=0.58423739467217718 ;vvr3(4)=-0.6666666666666666;vvwg(4)=0.81
		vvr1(5)=0.0d0 ;vvr2(5)=0.0d0 ;vvr3(5)=0.4d0;vvwg(5)=3.76d0
		
 

 
  	


  	
 case(3)
  a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006d0*0.515003019323671498
  b1=1.104848006d0*0.2571837452420646589
  c1=1.104848006d0*2.474004977113405936
  d1=1.104848006d0*0.419515737191525950
 
 


		vvr1(1)=-a ;vvr2(1)=-a ;vvr3(1)=d ;vvwg(1)=a1
		vvr1(2)=a ;vvr2(2)=-a ;vvr3(2)=d ;vvwg(2)=a1
		vvr1(3)=a ;vvr2(3)=a ;vvr3(3)=d ;vvwg(3)=a1
		vvr1(4)=-a ;vvr2(4)=a ;vvr3(4)=d ;vvwg(4)=a1
		vvr1(5)=-b ;vvr2(5)=0.0 ;vvr3(5)=e ;vvwg(5)=b1
		vvr1(6)=b ;vvr2(6)=0.0 ;vvr3(6)=e ;vvwg(6)=b1
		vvr1(7)=0.0 ;vvr2(7)=-b ;vvr3(7)=e ;vvwg(7)=b1
		vvr1(8)=0.0 ;vvr2(8)=b ;vvr3(8)=e ;vvwg(8)=b1
		vvr1(9)=0.0 ;vvr2(9)=0.0 ;vvr3(9)=f ;vvwg(9)=c1
		vvr1(10)=-c ;vvr2(10)=-c ;vvr3(10)=g ;vvwg(10)=d1
		vvr1(11)=c ;vvr2(11)=-c ;vvr3(11)=g ;vvwg(11)=d1
		vvr1(12)=c ;vvr2(12)=c ;vvr3(12)=g ;vvwg(12)=d1
		vvr1(13)=-c ;vvr2(13)=c ;vvr3(13)=g ;vvwg(13)=d1










			
 case(4)
  a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006*0.515003019323671498
  b1=1.104848006*0.2571837452420646589
  c1=1.104848006*2.474004977113405936
  d1=1.104848006*0.419515737191525950
 
 


		vvr1(1)=-a ;vvr2(1)=-a ;vvr3(1)=d ;vvwg(1)=a1
		vvr1(2)=a ;vvr2(2)=-a ;vvr3(2)=d ;vvwg(2)=a1
		vvr1(3)=a ;vvr2(3)=a ;vvr3(3)=d ;vvwg(3)=a1
		vvr1(4)=-a ;vvr2(4)=a ;vvr3(4)=d ;vvwg(4)=a1
		vvr1(5)=-b ;vvr2(5)=0.0 ;vvr3(5)=e ;vvwg(5)=b1
		vvr1(6)=b ;vvr2(6)=0.0 ;vvr3(6)=e ;vvwg(6)=b1
		vvr1(7)=0.0 ;vvr2(7)=-b ;vvr3(7)=e ;vvwg(7)=b1
		vvr1(8)=0.0 ;vvr2(8)=b ;vvr3(8)=e ;vvwg(8)=b1
		vvr1(9)=0.0 ;vvr2(9)=0.0 ;vvr3(9)=f ;vvwg(9)=c1
		vvr1(10)=-c ;vvr2(10)=-c ;vvr3(10)=g ;vvwg(10)=d1
		vvr1(11)=c ;vvr2(11)=-c ;vvr3(11)=g ;vvwg(11)=d1
		vvr1(12)=c ;vvr2(12)=c ;vvr3(12)=g ;vvwg(12)=d1
		vvr1(13)=-c ;vvr2(13)=c ;vvr3(13)=g ;vvwg(13)=d1


			
case(5)

 a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006*0.515003019323671498
  b1=1.104848006*0.2571837452420646589
  c1=1.104848006*2.474004977113405936
  d1=1.104848006*0.419515737191525950
 
 


		vvr1(1)=-a ;vvr2(1)=-a ;vvr3(1)=d ;vvwg(1)=a1
		vvr1(2)=a ;vvr2(2)=-a ;vvr3(2)=d ;vvwg(2)=a1
		vvr1(3)=a ;vvr2(3)=a ;vvr3(3)=d ;vvwg(3)=a1
		vvr1(4)=-a ;vvr2(4)=a ;vvr3(4)=d ;vvwg(4)=a1
		vvr1(5)=-b ;vvr2(5)=0.0 ;vvr3(5)=e ;vvwg(5)=b1
		vvr1(6)=b ;vvr2(6)=0.0 ;vvr3(6)=e ;vvwg(6)=b1
		vvr1(7)=0.0 ;vvr2(7)=-b ;vvr3(7)=e ;vvwg(7)=b1
		vvr1(8)=0.0 ;vvr2(8)=b ;vvr3(8)=e ;vvwg(8)=b1
		vvr1(9)=0.0 ;vvr2(9)=0.0 ;vvr3(9)=f ;vvwg(9)=c1
		vvr1(10)=-c ;vvr2(10)=-c ;vvr3(10)=g ;vvwg(10)=d1
		vvr1(11)=c ;vvr2(11)=-c ;vvr3(11)=g ;vvwg(11)=d1
		vvr1(12)=c ;vvr2(12)=c ;vvr3(12)=g ;vvwg(12)=d1
		vvr1(13)=-c ;vvr2(13)=c ;vvr3(13)=g ;vvwg(13)=d1
	
			
case(6,7,8,9)

 a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006*0.515003019323671498
  b1=1.104848006*0.2571837452420646589
  c1=1.104848006*2.474004977113405936
  d1=1.104848006*0.419515737191525950
 
 


		vvr1(1)=-a ;vvr2(1)=-a ;vvr3(1)=d ;vvwg(1)=a1
		vvr1(2)=a ;vvr2(2)=-a ;vvr3(2)=d ;vvwg(2)=a1
		vvr1(3)=a ;vvr2(3)=a ;vvr3(3)=d ;vvwg(3)=a1
		vvr1(4)=-a ;vvr2(4)=a ;vvr3(4)=d ;vvwg(4)=a1
		vvr1(5)=-b ;vvr2(5)=0.0 ;vvr3(5)=e ;vvwg(5)=b1
		vvr1(6)=b ;vvr2(6)=0.0 ;vvr3(6)=e ;vvwg(6)=b1
		vvr1(7)=0.0 ;vvr2(7)=-b ;vvr3(7)=e ;vvwg(7)=b1
		vvr1(8)=0.0 ;vvr2(8)=b ;vvr3(8)=e ;vvwg(8)=b1
		vvr1(9)=0.0 ;vvr2(9)=0.0 ;vvr3(9)=f ;vvwg(9)=c1
		vvr1(10)=-c ;vvr2(10)=-c ;vvr3(10)=g ;vvwg(10)=d1
		vvr1(11)=c ;vvr2(11)=-c ;vvr3(11)=g ;vvwg(11)=d1
		vvr1(12)=c ;vvr2(12)=c ;vvr3(12)=g ;vvwg(12)=d1
		vvr1(13)=-c ;vvr2(13)=c ;vvr3(13)=g ;vvwg(13)=d1



	
			
end select
		qpoints(:,:)=0.0d0


		do kk=1,qp_pyra
			wequa3d(kk)=vvwg(kk)*0.1250000000000
			
			r=vvr1(kk); s=vvr2(kk); tx=vvr3(kk)
			vvnxi(1)=(0.1250000000000)*(1.0-r)*(1.0d0-s)*(1.0d0-tx)
			vvnxi(2)=(0.1250000000000)*(1.0+r)*(1.0d0-s)*(1.0d0-tx)
			vvnxi(3)=(0.1250000000000)*(1.0+r)*(1.0d0+s)*(1.0d0-tx)
			vvnxi(4)=(0.1250000000000)*(1.0-r)*(1.0d0+s)*(1.0d0-tx)
			vvnxi(5)=0.5d0*(1.0d0+tx)
			

			do j=1,5
			qpoints(:,kk)=qpoints(:,kk)+(vvnxi(j)*vext(j,:))
			end do
		    
			

		end do




end subroutine quadraturepyra


subroutine quadraturehexa(n,igqrules,vext,qpoints,wequa3d)
 !> @brief
!> this subroutine computes the quadrature points and weights for a hexahedral
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::igqrules,n
real,dimension(1:8,1:dimensiona),intent(in)::vext
real,dimension(1:dimensiona,1:numberofpoints),intent(inout)::qpoints
real,dimension(1:numberofpoints),intent(inout)::wequa3d
real::r,s,tx,a,b,c,d,e,f
real::a1,b1,c1,d1,e1,f1
integer::kk,j,ii,ij,ik,count1
real,dimension(1:3,1:3)::vva,vva1
real,dimension(1)::deta
real,dimension(1:8)::vvnxi
real,dimension(1:alls)::vvwg,vvr1,vvr2,vvr3
real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz



 wequa3d=0.0d0
  qpoints=0.0d0

select case(igqrules)
 

 case(1)

		vvwg(1) = 8.0d0
	    vvr1(1)=0.0d0	;vvr2(1)=0.0d0	;vvr3(1)=0.0d0


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
  vvnpoy(1)= a	;vvnpoy(2)=b	
  vvnpoz(1)= a	;vvnpoz(2)=b

  vvwpox(1)=a1	;vvwpox(2)=b1	
  vvwpoy(1)= a1	;vvwpoy(2)=b1	
  vvwpoz(1)= a1	;vvwpoz(2)=b1
  
  count1=0
 do ii=1,2
    do ij=1,2
      do ik=1,2
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) ;vvr3(count1)=vvnpoz(ik)   
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do
  	

  	
 case(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1
  count1=0
 do ii=1,3
    do ij=1,3
      do ik=1,3
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) ;vvr3(count1)=vvnpoz(ik)   
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do
 

			
 case(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1
  count1=0
 do ii=1,4
    do ij=1,4
      do ik=1,4
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) ;vvr3(count1)=vvnpoz(ik)  
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do


			
case(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  e1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1
  count1=0
 do ii=1,5
    do ij=1,5
      do ik=1,5
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) ;vvr3(count1)=vvnpoz(ik) 
	vvwg(count1)=(vvwpox(ii)*vvwpox(ij)*vvwpox(ik))
      end do
    end do
end do
	
			
case(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  f=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  f1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e;vvnpoy(6)=f
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e;vvnpoz(6)=f
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1 ;vvwpoy(6)=f1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1 ;vvwpoz(6)=f1
  count1=0
 do ii=1,6
    do ij=1,6
      do ik=1,6
	count1=count1+1
	vvr1(count1)=vvnpox(ii);vvr2(count1)=vvnpoy(ij) ;vvr3(count1)=vvnpoz(ik) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do
	
			
end select
		qpoints(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.125d0
		do kk=1,qp_hexa
			wequa3d(kk)=vvwg(kk)
			r=vvr1(kk); s=vvr2(kk); tx=vvr3(kk)
			vvnxi(1)=(0.125d0)*(1.0d0-r)*(1.0d0-s)*(1.0d0-tx)
			vvnxi(2)=(0.125d0)*(1.0d0+r)*(1.0d0-s)*(1.0d0-tx)
			vvnxi(3)=(0.125d0)*(1.0d0+r)*(1.0d0+s)*(1.0d0-tx)
			vvnxi(4)=(0.125d0)*(1.0d0-r)*(1.0d0+s)*(1.0d0-tx)
			vvnxi(5)=(0.125d0)*(1.0d0-r)*(1.0d0-s)*(1.0d0+tx)
			vvnxi(6)=(0.125d0)*(1.0d0+r)*(1.0d0-s)*(1.0d0+tx)
			vvnxi(7)=(0.125d0)*(1.0d0+r)*(1.0d0+s)*(1.0d0+tx)
			vvnxi(8)=(0.125d0)*(1.0d0-r)*(1.0d0+s)*(1.0d0+tx)
			do j=1,8
			qpoints(:,kk)=qpoints(:,kk)+(vvnxi(j)*vext(j,:))
			end do
			

		end do




 

end subroutine quadraturehexa

subroutine rotatef(n,rotvect,vectco,angle1,angle2)
 !> @brief
!> this subroutine rotates the vector of fluxes in the directions normal to the face in 3d 
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:5,1:5)::tri
real,dimension(1:nof_variables),intent(inout)::rotvect
real,dimension(1:nof_variables),intent(inout)::vectco
real,intent(in)::angle1,angle2
real::sia1,coa1,coa2,sia2


!build matrix of rotation!
 coa1=cos(angle1)
 sia1=sin(angle1)
 coa2=cos(angle2)
 sia2=sin(angle2)
 
tri=zero

tri(1,1)=1.0d0

tri(2,2)=coa1*sia2	!cos(angle1)*sin(angle2)
tri(2,3)=sia1*sia2	!sin(angle1)*sin(angle2)
tri(2,4)=coa2		!cos(angle2)

tri(3,2)=coa1*coa2	!cos(angle1)*cos(angle2)
tri(3,3)=sia1*coa2	!sin(angle1)*cos(angle2)
tri(3,4)=-sia2		!-sin(angle2)

tri(4,2)=-sia1		!-sin(angle1)
tri(4,3)=coa1		!cos(angle1)
tri(5,5)=1.0d0

rotvect(1:nof_variables)=vectco(1:nof_variables)

rotvect(1:5)=matmul(tri(1:5,1:5),vectco(1:5))





end subroutine rotatef


subroutine rotateb(n,rotvect,vectco,angle1,angle2)
 !> @brief
!> this subroutine rotates back the vector of fluxes from the directions normal to the face to cartesian coordinates
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:5,1:5)::invtri
real,dimension(1:nof_variables),intent(inout)::rotvect
real,dimension(1:nof_variables),intent(inout)::vectco
real,intent(in)::angle1,angle2
real::sia1,coa1,coa2,sia2

!build matrix of rotation!
 coa1=cos(angle1)
 sia1=sin(angle1)
 coa2=cos(angle2)
 sia2=sin(angle2)




invtri=zero


!build matrix of rotation!
invtri(1,1)=1.0d0

invtri(2,2)=coa1*sia2!cos(angle1)*sin(angle2)
invtri(2,3)=coa1*coa2!cos(angle1)*cos(angle2)
invtri(2,4)=-sia1!-sin(angle1)

invtri(3,2)=sia1*sia2!sin(angle1)*sin(angle2)
invtri(3,3)=sia1*coa2!sin(angle1)*cos(angle2)
invtri(3,4)=coa1!cos(angle1)

invtri(4,2)=coa2!cos(angle2)
invtri(4,3)=-sia2!-sin(angle2)
invtri(5,5)=1.0d0

rotvect(1:nof_variables)=vectco(1:nof_variables)

rotvect(1:5)=matmul(invtri(1:5,1:5),vectco(1:5))



end subroutine rotateb


subroutine rotatef2d(n,rotvect,vectco,angle1,angle2)
 !> @brief
!> this subroutine rotates the vector of fluxes in the directions normal to the edge in 2d
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::rotvect
real,dimension(1:nof_variables),intent(inout)::vectco
real,intent(in)::angle1,angle2

rotvect(1:nof_variables)=vectco(1:nof_variables)


rotvect(1)=vectco(1)
rotvect(2)=(angle1*vectco(2))+(angle2*vectco(3))
rotvect(3)=-(angle2*vectco(2))+(angle1*vectco(3))
rotvect(4)=vectco(4)






end subroutine rotatef2d


subroutine rotateb2d(n,rotvect,vectco,angle1,angle2)
 !> @brief
!> this subroutine rotates back the vector of fluxes from the directions normal to the edge to cartesian coordinates
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::rotvect
real,dimension(1:nof_variables),intent(inout)::vectco
real,intent(in)::angle1,angle2


rotvect(1:nof_variables)=vectco(1:nof_variables)

!build matrix of rotation!

rotvect(1)=vectco(1)
rotvect(2)=(angle1*vectco(2))-(angle2*vectco(3))
rotvect(3)=(angle2*vectco(2))+(angle1*vectco(3))
rotvect(4)=vectco(4)



end subroutine rotateb2d




subroutine probepos(n,probei)
 !> @brief
!> this subroutine establishes the cells where the probe positions belong to
implicit none
integer,allocatable,dimension(:,:),intent(inout)::probei
integer,intent(in)::n
integer::i,j,k,l,kmaxe,inv
real,dimension(1:8,1:dimensiona)::vext
real::dist
real::dumin,dumout,delta
real,dimension(1:dimensiona)::cords

allocate (probei(n:n,nprobes))


probei(n:n,:)=0

kmaxe=xmpielrank(n)

	if (dimensiona.eq.3)then
	do inv=1,nprobes
	delta=tolbig
	do i=1,kmaxe

		
		
		
		call compute_centre3d(i,cords)
		vext(1,1:3)=cords(1:3)
		vext(2,1:3)=probec(inv,1:3)
		dist=distance3(n,vext)
		
		
		if (dist.lt.delta) then
			delta=dist
			l=i
		end if
	end do
	
	call mpi_barrier(mpi_comm_world,ierror)
	dumout=delta
		dumin=0.0d0
		call mpi_allreduce(dumout,dumin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
		call mpi_barrier(mpi_comm_world,ierror)
	if (abs(dumin-delta).le.1.0e-15) then
	probei(n,inv)=l
	end if
	end do

	else
      	do inv=1,nprobes
	delta=tolbig
	do i=1,kmaxe

		
		call compute_centre2d(i,cords)
		vext(1,1:2)=cords(1:2)
		vext(2,1:2)=probec(inv,1:2)
		dist=distance2(n,vext)
		
		
		if (dist.le.delta) then
			delta=dist
			l=i
		end if
	end do
	
	call mpi_barrier(mpi_comm_world,ierror)
	
	dumout=delta
		dumin=0.0d0
		call mpi_allreduce(dumout,dumin,1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
		call mpi_barrier(mpi_comm_world,ierror)
	if (abs(dumin-delta).le.1.0e-15) then
	probei(n,inv)=l
	
	end if
	end do
	end if

	
end subroutine probepos


subroutine  anglex(a_rot,b_rot,anglefacex)
 !> @brief
!> this subroutine computes the angles according to the quadrant sign
implicit none
#ifdef gpu
!$omp declare target
#endif
real,intent(inout)::anglefacex
real,intent(in)::a_rot,b_rot

if ((a_rot.ne.zero).and.(b_rot.ne.zero))then
	if ((a_rot.gt.zero).and.(b_rot.gt.zero))then
		anglefacex=atan(b_rot/a_rot)
	end if
	if ((a_rot.lt.zero).and.(b_rot.gt.zero))then
		anglefacex=pi+(atan(b_rot/a_rot))
	end if
	if ((a_rot.lt.zero).and.(b_rot.lt.zero))then
		anglefacex=pi+(atan(b_rot/a_rot))
	end if
	if ((a_rot.gt.zero).and.(b_rot.lt.zero))then
		anglefacex=(2.0d0*pi)+(atan(b_rot/a_rot))
	end if
end if
if ((a_rot.eq.zero).and.(b_rot.ne.zero))then
	if (b_rot.gt.zero) then
		anglefacex=(pi/2.0d0)
	end if
	if (b_rot.lt.zero) then
		anglefacex=3.0d0*(pi/2.0d0)
	end if
end if
if ((a_rot.ne.zero).and.(b_rot.eq.zero))then
	if (a_rot.gt.zero) then
		anglefacex=zero
	end if
	if (a_rot.lt.zero) then
		anglefacex=pi
	end if
end if
end subroutine anglex
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!function to calculate the angles based on the coordinates !!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!of the normal plane!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine angley(c_rot,root_rot,anglefacey)
 !> @brief
!> this subroutine computes the angles according to the quadrant sign
implicit none
#ifdef gpu
!$omp declare target
#endif
real,intent(inout)::anglefacey
real,intent(in)::c_rot,root_rot

if (c_rot.eq.zero)then
	anglefacey=acos(zero)
end if
if (c_rot.ne.zero)then
	if (c_rot.gt.zero)then
		anglefacey=acos(c_rot/root_rot)
	end if
	if (c_rot.lt.zero)then
		anglefacey=acos(c_rot/root_rot)
	end if
end if
end subroutine angley


subroutine angle2d(vext,anglefacex,anglefacey)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,intent(inout)::anglefacex,anglefacey
real,dimension(1:8,1:dimensiona),intent(in)::vext
real::length

length=distance2(n,vext)

anglefacex=(vext(2,2)-vext(1,2))/length
anglefacey=-(vext(2,1)-vext(1,1))/length



end subroutine angle2d










end module transform
