module fluxes
use library
use transform
use local
use riemann
use flow_operations
use declaration
use basis
use dg_functions
implicit none


 contains
 
 
subroutine solution_integ(i,solution_integ2)
 implicit none
#ifdef gpu
!$omp declare target
#endif
 real,dimension(1:nof_variables),intent(inout)::solution_integ2
 integer,intent(in)::i
 
 solution_integ2=dg_vol_integral2(n,i)
 
 
 end subroutine solution_integ

 subroutine solution_integ_s(i,solution_integ_strong)
 implicit none
#ifdef gpu
!$omp declare target
#endif
  real,dimension(1:nof_variables),intent(inout)::solution_integ_strong
 integer,intent(in)::i

 solution_integ_strong=dg_vol_integral_strong(n,i)


 end subroutine solution_integ_s


  subroutine solution_integ_w(i,solution_integ_weak)
 implicit none
#ifdef gpu
!$omp declare target
#endif
 real,dimension(1:nof_variables),intent(inout)::solution_integ_weak
 integer,intent(in)::i

 solution_integ_weak=dg_vol_integral_weak(n,i)


 end subroutine solution_integ_w




subroutine calculate_fluxeshi(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe
		real::godflux2,sum_detect
		integer::l,ngp,iqp
		real,dimension(1:numberofpoints2)::weights_temp
		real,dimension(1:8,1:dimensiona)::vext
		real::angle1,angle2,nx,ny,nz,normalvect
		integer::facex,pointx,iconsidered,nfx,lfx,rowfx
		real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,hllcflux,rhllcflux
		real,dimension(1:idegfree+1,1:nof_variables)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ
	kmaxe=xmpielrank(n)
	



#ifdef gpu
! !$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do i=1,kmaxe


        iconsidered=i
        if (dg.eq.1)then
        rhs_valdg(:,:,i) = 0.0d0
        dg_rhs = 0.0d0
        dg_rhs_surf_integ = 0.0d0
        dg_rhs_vol_integ = 0.0d0
        end if
	
        if (dg.eq.1) then
            dg_rhs_vol_integ = dg_vol_integral(n,iconsidered)
            
        end if
	
				  
		if (ielem_interior(i).eq.0)then
		    rhs_val(1,i)=zero
		    
		    do l=1,ielem_ifca(i)
				  godflux2=zero
				  angle1=ielem_faceanglex(l,i)
				  angle2=ielem_faceangley(l,i)
				  facex=l
				  normalvect=((cos(angle1)*sin(angle2))*lamx)+((sin(angle1)*sin(angle2))*lamy)+((cos(angle2))*lamz)
				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad
					weights_temp(1:iqp)=weights_q(1:iqp)
					
				  else
					iqp=qp_triangle
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  

				  
				  
				  do ngp=1,iqp
				  pointx=ngp
				  if (dg.eq.1) then
                        cleft = rec_uleft_dg(1:nof_variables, l, ngp,i)
                        cright = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                    else !fv
				  
				  
				  
				      cleft(1)=rec_uleft(1,l,ngp,i)
				      cright(1)=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                    end if
				      call exact_riemann_solver(n,cleft,cright,normalvect,hllcflux)
				      if (dg.eq.1)then
				      rhllcflux(1)=hllcflux(1)
				      
				      dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
				      
				       
				      
				      else
				      godflux2=godflux2+(hllcflux(1)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				  end do
				    rhs_val(1,i)=rhs_val(1,i)+godflux2

		    end do
		end if
		if (ielem_interior(i).eq.1)then
		    rhs_val(1,i)=zero
		    
		    do l=1,ielem_ifca(i)
				      facex=l
				  angle1=ielem_faceanglex(l,i)
				  angle2=ielem_faceangley(l,i)
				  normalvect=((cos(angle1)*sin(angle2))*lamx)+((sin(angle1)*sin(angle2))*lamy)+((cos(angle2))*lamz)
 				  
 				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad
					weights_temp(1:iqp)=weights_q(1:iqp)
				  else
					iqp=qp_triangle
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  
				  

				  godflux2=zero
				  do ngp=1,iqp
				  pointx = ngp
                    if (dg == 1) then
                        cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
                        else
				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
				      end if
					    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                 
								  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								 if (dg == 1) then
                                    cright = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                else  
                                cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
								end if
								  else
								  !not periodic ones in my cpu
								  cright(1:nof_variables)=cleft(1:nof_variables)
								    
								  end if
							else
							     if (dg == 1) then
                                cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))

                            else  
                            cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                            end if
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
									 if (dg == 1) then
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                                    nfx  = ielem_ineighn(l,i)
                                    lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                    rowfx = bound_offset(nfx) + lfx - 1
                                    cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                else 
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
								end if
								end if
							else 								
								 if (dg == 1) then
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                else
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                                end if
! 								   
							end if
					    end if
				      
				      call exact_riemann_solver(n,cleft,cright,normalvect,hllcflux)
				      
				      if (dg.eq.1)then
				      rhllcflux(1)=hllcflux(1)
				      dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
				      
				       
				     
				      else
				      
				      godflux2=godflux2+(hllcflux(1)*(weights_temp(ngp)*ielem_surf(l,i)))	
				      
				      end if
				  end do
				    rhs_val(1,i)=rhs_val(1,i)+godflux2
				    
		    end do
		end if
! 				
        if (dg.eq.1)then
        dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
        rhs_valdg(:,:,i) = rhs_valdg(:,:,i) + dg_rhs
        


        end if


	end do
#ifdef gpu
! !$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end subroutine calculate_fluxeshi
	
	
	
	
	
	
subroutine calculate_fluxeshi2d(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation in 2d
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe
	real::godflux2,sum_detect,lamxl,lamyl
	integer::l,k,ngp,iqp, neighbor_index, neighbor_face_index
	real,dimension(1:8,1:dimensiona)::vext
	 real::angle1,angle2,nx,ny,nz,normalvect
    integer::facex,pointx,iconsidered,nfx,lfx,rowfx
    real,dimension(1:numberofpoints2)::weights_temp
    real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,hllcflux,rhllcflux
    real,dimension(1:idegfree+1,1:nof_variables)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ


	kmaxe = xmpielrank(n)


#ifdef gpu
! !$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do i=1,kmaxe


	weights_temp(1:qp_line_n)=weights_l(1:qp_line_n)
	
        if (initcond.eq.3)then
            lamxl=-ielem_yyc(i)+0.5d0
            lamyl=ielem_xxc(i)-0.5
        else
			lamxl=lamx
			lamyl=lamy
		end if
        
        
        if (dg.eq.1)then
        rhs_valdg(:,:,i) = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
        else
        rhs_val(:,i)=zero
        end if
        
        
        iconsidered=i
        
        if (dg.eq.1) then
        
            dg_rhs_vol_integ = dg_vol_integral(n,iconsidered)
         
        end if
        
		if (ielem_interior(i).eq.0)then ! element is interior

		    do l=1,ielem_ifca(i)
                godflux2=zero
                nx=ielem_faceanglex(l,i)
                ny=ielem_faceangley(l,i)
                facex=l
                
                normalvect=(nx*lamxl)+(ny*lamyl)

                iqp=qp_line_n
                
                neighbor_index = ielem_ineigh(l,i)
                neighbor_face_index = ielem_ineighn(l,i)
                
                do ngp=1,iqp
                    pointx=ngp


                    if (dg.eq.1) then
                        cleft = rec_uleft_dg(1:nof_variables, l, ngp,i)
                        cright = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                    else !fv
                        cleft(1)=rec_uleft(1,l,ngp,i)
                        cright(1)=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                    end if
                    
                    call exact_riemann_solver(n,cleft,cright,normalvect,hllcflux)




                    if (dg.eq.1) then
                        ! riemann flux at interface quadrature points times basis
                        rhllcflux(1)=hllcflux(1)
                         
                        dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)



                    else !fv
                        godflux2=godflux2+(hllcflux(1)*(weights_temp(ngp)*ielem_surf(l,i)))
                    end if
                end do
                
                if (dg /= 1) rhs_val(1,i)=rhs_val(1,i)+godflux2





		    end do



		
        else if (ielem_interior(i).eq.1)then

            do l=1,ielem_ifca(i)
                facex = l
                nx=ielem_faceanglex(l,i)
                ny=ielem_faceangley(l,i)
                normalvect=(nx*lamxl)+(ny*lamyl)
                iqp=qp_line_n
                
                godflux2=zero
                do ngp=1,iqp
                    pointx = ngp

                    if (dg == 1) then
                        cleft = rec_uleft_dg(1:nof_variables, l, ngp,i)
                    else
                        cleft(1)=rec_uleft(1,l,ngp,i)
                    end if
                    
                    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                            if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                if (dg == 1) then
                                    cright = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                else !fv
                                    cright(1) = rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                end if
                            else !not periodic ones in my cpu
                                cright(1:nof_variables)=cleft(1:nof_variables)
                            end if
!                            
                        else
                            if (dg == 1) then
                                cright = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
!                                 
                            else !fv
                                cright(1) = rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                            end if
                        end if
                    else !in other cpus they can only be periodic or mpi neighbours
                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                            if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                            
                                if (dg == 1) then
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                                    nfx  = ielem_ineighn(l,i)
                                    lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                    rowfx = bound_offset(nfx) + lfx - 1
                                    cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                else
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                                    nfx  = ielem_ineighn(l,i)
                                    lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                    rowfx = bound_offset(nfx) + lfx - 1
                                    cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                                end if
                            end if
                        else
                            if (dg == 1) then
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                            else
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                            end if
                        end if
                        
                    end if
                    
                    call exact_riemann_solver(n,cleft,cright,normalvect,hllcflux)



                    
                    if (dg.eq.1) then
                        !riemann flux at interface quadrature points times basis
                        rhllcflux(1)=hllcflux(1)
                         dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
                         




                    else !fv
                        godflux2=godflux2+(hllcflux(1)*(weights_temp(ngp)*ielem_surf(l,i)))
                    end if
  
                end do
                
                if (dg /= 1) rhs_val(1,i)=rhs_val(1,i)+godflux2
		    end do



		end if
		
		
		if (dg == 1) dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ



		
        if (dg == 1) rhs_valdg(:,:,i) = rhs_valdg(:,:,i) + dg_rhs
        


	end do
#ifdef gpu
! !$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end subroutine calculate_fluxeshi2d

subroutine calculate_fluxeshi_convective(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe,ii
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::l,ngp,iqp,ikas,igoflux, icaseb,jx,jx2,b_code,srf
	real::sum_detect,norms,tempxx
	real,dimension(1:numberofpoints2)::weights_temp
	integer::iconsidered,facex,pointx,nfx,lfx,rowfx
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:nof_variables)::leftv,rightv,srf_speedrot
	real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
	real,dimension(1:dimensiona)::pox,poy,poz
	real,dimension(1:nof_variables)::srf_speed
	real,dimension(1:8,1:dimensiona)::vext
    real::x1,y1,z1
    real,dimension(1:idegfree+1,1:nof_variables)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ

	kmaxe=xmpielrank(n)


#ifdef gpu
! !$omp target teams distribute parallel do private(i,&
! !$omp godflux2, dg_vol_rec,rhllcflux,hllcflux,l,ngp,iqp,ikas,igoflux, icaseb,jx,jx2,b_code,&
! !$omp srf,sum_detect,norms,tempxx,weights_temp,iconsidered,facex,pointx,nfx,lfx,rowfx,&
! !$omp angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3,cleft,cright,cleft_rot,cright_rot,&
! !$omp leftv,rightv,srf_speedrot,cturbl,srf_speed,cturbr,pox,poy,poz,vext,x1,y1,z1,dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	i=el_int(ii)

	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero




	iconsidered=i
            mp_source3=zero        

        if(mrf.eq.1)then
            srf=rec_mrf(i)
        end if
		     if (dg.eq.1) then
		     rhs_valdg(:,:,i) = zero
            dg_rhs = zero
            dg_rhs_surf_integ = zero
            dg_rhs_vol_integ = zero
            
             dg_rhs_vol_integ = dg_vol_integral(n,iconsidered)

            end if
            
		    
		    do l=1,ielem_ifca(i) !for all their faces
				b_code=0

				  godflux2=zero
				  mp_source2=zero
 				  angle1=ielem_faceanglex(l,i)
 				  angle2=ielem_faceangley(l,i)
 				  nx=(cos(angle1)*sin(angle2))
				  ny=(sin(angle1)*sin(angle2))
				  nz=(cos(angle2))
				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad_n
					weights_temp(1:iqp)=weights_q(1:iqp)
				  else
					iqp=qp_triangle_n
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  
				  facex=l

				  
				  do ngp=1,iqp	!for all the gaussian quadrature points
				  pointx=ngp


		call get_states_interior(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)

				  

			
						  call rotatef(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  end if
						  
				      
				      select case(iriemann)
				      
				      case(1)			!hllc



				      
				      call hllc_riemann_solver(n,iconsidered, facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)



				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if

				      case(9)			!hll

				      call hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb(n,rhllcflux,hllcflux,angle1,angle2)






				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      case(3)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				     				      
				      case(4)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				      
				       
				       
				       
				       case(5)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				      
				       end select
				       
				      if (dg.eq.1)then
				      
				      dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
				      
				      
				      
				      else 
				       
				       
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      
				      end if
				      
				      
				      
				      
				      
				      
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				     
				     
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					  
					  
					  
					  norms=0.5*(cleft_rot(2)+cright_rot(2))


					if (rec_mrf(iconsidered).eq.1)then
					  norms=norms-srf_speedrot(2)
					end if
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
					  

					  end if
					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
					  
					 
					  
					  
				      end if
				      
				      
				  end do
				    
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				      if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if
				    

				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

				    end if

		    end do
                 if (multispecies.eq.1)then
                 rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
                 
                 end if
                 
                 if (dg.eq.1)then

                   dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
                    rhs_valdg(:,:,i) = rhs_valdg(:,:,i) + dg_rhs
                    
                    if (multispecies.eq.1)then


						rhs_valdg(1,8,i)=rhs_valdg(1,8,i)-(u_c_val(1,8,i)*mp_source3)

					end if

                   




                 end if

	end do
#ifdef gpu
! !$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	



#ifdef gpu
! ! !$omp target teams distribute parallel do private(i,&
! ! !$omp godflux2, dg_vol_rec,rhllcflux,hllcflux,l,ngp,iqp,ikas,igoflux, icaseb,jx,jx2,b_code,&
! ! !$omp srf,sum_detect,norms,tempxx,weights_temp,iconsidered,facex,pointx,nfx,lfx,rowfx,&
! ! !$omp angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3,cleft,cright,cleft_rot,cright_rot,&
! ! !$omp leftv,rightv,srf_speedrot,cturbl,srf_speed,cturbr,pox,poy,poz,vext,x1,y1,z1,dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)


	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero



	iconsidered=i	
		 mp_source3=zero		
				
	if(mrf.eq.1)then
        srf=rec_mrf(i)
    end if	   
		    if (dg.eq.1) then
		    rhs_valdg(:,:,i)=zero
            dg_rhs = zero
            dg_rhs_surf_integ = zero
            dg_rhs_vol_integ = zero
            
             dg_rhs_vol_integ = dg_vol_integral(n,iconsidered)
                
            end if

		    do l=1,ielem_ifca(i)
		    facex=l




		    
				      
				  angle1=ielem_faceanglex(l,i)
				  angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
 				  
 				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad_n
					weights_temp(1:iqp)=weights_q(1:iqp)
				  else
					iqp=qp_triangle_n
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  

				  
				  godflux2=zero
				  
				   mp_source2=zero
								  
				  
				  do ngp=1,iqp
				  pointx = ngp
				  
				      b_code=0




				      call get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
				      

				      
				      
			
						  call rotatef(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
! 						
						  
						  
								    
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  end if
						  
				      
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
								      
								      cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
								      cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
				      
				      
				      
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver(n,iconsidered, facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      



				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if

				      case(9)			!hllc

				      call hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if




				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb(n,rhllcflux,hllcflux,angle1,angle2)




				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(3)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				      
				      
				       case(4)		!roe
				       
				       
				       
				       if (b_code.gt.0)then
				       
				       call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				       
				       
				       else
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				     
				      
				      
				      call rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      
				       rhllcflux=hllcflux
				       
				       end if
				     			
                                    case(5)			!roe
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      
				     if (b_code.le.0)then
				      
				      call troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      else
				      
				      
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      end if
				      
				       rhllcflux=hllcflux
				     
				      
				       end select
				       
				     if (dg.eq.1)then
				      dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
				      
				       
				      
				      else
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      
				      end if
				      
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				       
				      
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					   
					  
					   
                     norms=0.5*(cleft_rot(2)+cright_rot(2))

				  if (rec_mrf(iconsidered).eq.1)then
                            norms=norms-srf_speedrot(2)
					end if
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
					  
					  
					  end if
					  
					  if ((b_code.eq.3).or.(b_code.eq.4))then
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
					  
					  end if
					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
					  
					  
					  
					  
					  
					  
				      end if
				      
				      
				  end do
				   
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				    if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

				    end if
! 				    end if
		    end do
		     if (multispecies.eq.1)then
                 rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
                 
                 end if
                 
                    if (dg.eq.1)then
                    dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
                    rhs_valdg(:,:,i) = rhs_valdg(:,:,i) + dg_rhs
                    
                    if (multispecies.eq.1)then

						rhs_valdg(1,8,i)=rhs_valdg(1,8,i)-(u_c_val(1,8,i)*mp_source3)


					end if



                    end if

	end do
#ifdef gpu
! !$omp end target teams distribute parallel do
#else
!$omp end do
#endif








end subroutine calculate_fluxeshi_convective
	
	
subroutine calculate_fluxeshi_convective2d(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe,ii
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::l,ngp,iqp,ikas,igoflux, icaseb,kxk,b_code
	real::sum_detect,norms
	real,dimension(1:numberofpoints2)::weights_temp
	integer::iconsidered,facex,pointx
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:nof_variables)::leftv,rightv,srf_speedrot
	real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
	real,dimension(1:dimensiona)::pox,poy,poz
	real::x1,y1,z1
	real,dimension(1:nof_variables)::srf_speed
	real,dimension(1:8,1:dimensiona)::vext
	real,dimension(1:idegfree+1,1:nof_variables)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ



	kmaxe=xmpielrank(n)
	


	







#ifdef gpu
! !$omp target teams distribute parallel do private(i)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	i=el_int(ii)


	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero

	iconsidered=i
            if (dg.eq.1) then
            rhs_valdg(:,:,i) = zero
            dg_rhs = zero
            dg_rhs_surf_integ = zero
            dg_rhs_vol_integ = zero
            
            
            dg_rhs_vol_integ = dg_vol_integral(n,iconsidered)
            
            
            end if
	

                    
                mp_source3=zero     
		    do l=1,ielem_ifca(i) !for all their faces
		    


				  godflux2=zero
				  mp_source2=zero
 				  nx=ielem_faceanglex(l,i)
 				  ny=ielem_faceangley(l,i)
 				  angle1=nx
 				  angle2=ny
				 b_code=0
					iqp=qp_line_n
				
				  do ngp=1,iqp	!for all the gaussian quadrature points
				  pointx=ngp
				   facex=l
		    pointx=ngp



				  call get_states_interior2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)



			
						  call rotatef2d(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef2d(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  
						  !realgas state
						  if (realgas.eq.1)then
						  call fix_conservative_state(cleft_rot)
						  call fix_conservative_state(cright_rot)

						  end if
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht2d(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  end if
						  
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if


				      case(9)			!hll

				      call hll_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if



				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
                      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      case(3)			!roe
				      
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny)
				      
				      rhllcflux(1:nof_variables)=hllcflux(1:nof_variables)
                                        
                                        
                                         case(4)			!roe
				      
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      call rroe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,b_code)
				      
				      rhllcflux=hllcflux
				     
				      
				       end select
				       
				       
				       
				       if (dg.eq.1) then
                        ! riemann flux at interface quadrature points times basis
                         

                         
                        dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
                       
                         

                    else !fv
				       
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      
				      end if
				      
				       if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					    
                     norms=0.5*(cleft_rot(2)+cright_rot(2))
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
					  
! 					  
					  end if
					  
 					
					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
					  
! 					 
					  
				      end if
				      
				      
				  end do
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				    if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if
				    
				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

! 				    
				    end if
! 				    end if
		    end do
		    if (multispecies.eq.1)then

                 rhs_val(7,i)=rhs_val(7,i)-(u_c_val(1,7,i)*mp_source3)!*ielem_totvolume(i))

                 end if
                 
               if (dg == 1) then

               
                  
               
               dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
                rhs_valdg(:,:,i) = rhs_valdg(:,:,i) + dg_rhs  
                
                if (multispecies.eq.1)then

					rhs_valdg(1,7,i)=rhs_valdg(1,7,i)-(u_c_val(1,7,i)*mp_source3)


                 
                 
                 
                 end if
                
              
                 end if
                 

	end do
#ifdef gpu
! !$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	

	
#ifdef gpu
! !$omp target teams distribute parallel do private(i)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)

	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero


	iconsidered=i	
	mp_source3=zero  
	
            if (dg.eq.1) then
            rhs_valdg(:,:,i)=zero
            dg_rhs = zero
            dg_rhs_surf_integ = zero
            dg_rhs_vol_integ = zero
            
            dg_rhs_vol_integ = dg_vol_integral(n,iconsidered)
               
            end if
				

		    
		    do l=1,ielem_ifca(i)
		    facex=l
				  igoflux=0


				    b_code=0  
				 nx=ielem_faceanglex(l,i)
 				  ny=ielem_faceangley(l,i)
 				  angle1=nx
 				  angle2=ny
				 
					iqp=qp_line_n
				  godflux2=zero
				  mp_source2=zero
				  do ngp=1,iqp
				  pointx=ngp
				  facex=l

				  call get_states_bounds2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
				  
				  
				  

				      
				      
			
						  call rotatef2d(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef2d(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
! 						   
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  
						  call lmacht2d(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  
						  
						  end if
						  
				      
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
								      
								      cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
								      cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
				      
				      
				      !realgas state
						  if (realgas.eq.1)then
						  call fix_conservative_state(cleft_rot)
						  call fix_conservative_state(cright_rot)

						  end if



				      
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if






				      case(9)			!hll

				      call hll_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(3)			!roe
				      
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny)
				      
				      rhllcflux=hllcflux
				     				      
				     
				        case(4)			!roe
				      
				      
				       
				      
				      if ((b_code.le.0))then
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      
				      
				      
				      call rroe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,b_code)
				      rhllcflux=hllcflux
				      else
				      

				     
				      
				      call rusanov_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				      
				      
				      end if
				      
				      
				       end select
				       
				      if (dg.eq.1) then
                        
                         dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
                         
                         else
				      
                         godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				       
				      
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  norms=0.5*(cleft_rot(2)+cright_rot(2))
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
! 					  
					  end if
					  
					  if ((b_code.eq.4).or.(b_code.eq.3))then
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
					  end if
					  
					  
					 
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
! 					  
					  
					  
				      end if
				      
				      
				  end do
				   
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				    if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

				    end if
		    end do
		     if (multispecies.eq.1)then
                 rhs_val(7,i)=rhs_val(7,i)-(u_c_val(1,7,i)*mp_source3)!*ielem_totvolume(i))
                 
                 end if
                 
                  if (dg == 1) then
                
               dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
                rhs_valdg(:,:,i) = rhs_valdg(:,:,i) + dg_rhs  
               
                if (multispecies.eq.1)then
                  rhs_valdg(1,7,i)=rhs_valdg(1,7,i)-(u_c_val(1,7,i)*mp_source3)

                 end if
                 
                 end if
                 
                 
                 

	end do
#ifdef gpu
! !!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine calculate_fluxeshi_convective2d

	
	

subroutine calculate_fluxeshi_diffusive(n)
!> @brief
!> this subroutine computes the diffusive fluxes for euler-navier-stokes equations
	implicit none
	integer,intent(in)::n
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ii,nvar,kc,iex,ittt,ikas,igoflux, icaseb,kk,b_code,srf,k
	real::sum_detect,norms
	integer::iconsidered,facex,pointx,rg_i,rg_j,idxy
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:nof_variables)::leftv,rightv,srf_speedrot,tempx_l,rtempx_l
	real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
	real,dimension(1:dimensiona)::pox,poy,poz,rg_sumfl,rg_sumfr
	real,dimension(1:nof_variables)::srf_speed
	real,dimension(1:numberofpoints2)::weights_temp
	real,dimension(1:8,1:dimensiona)::vext
	real::mp_pinfl,mp_pinfr,gammal,gammar,rho12l,rho12r
	real,dimension(1:4)::viscl,laml
	real,dimension(1:2)::turbmv
    real,dimension(1)::etvm
    real,dimension(1:20)::eddyfl,eddyfr
    real,dimension(1:nof_variables-1,1:dims)::lcvgrad,rcvgrad
	real,dimension(turbulenceequations+passivescalar,1:dims)::lcvgrad_t,rcvgrad_t
	real,dimension(1:nof_variables)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(1:nof_species)::rgs_htr_i, rgs_hvib_i,y_av
	real,dimension(3,3)::taul,taur,tau
	real,dimension(3)::q,nnn,nall,qvib
	real::ux,uy,uz,vx,vy,vz,wx,wy,wz,rho12,u12,v12,w12 ,damp,vdamp,tempxx 
	real,dimension(1:idegfree+1,1:nof_variables)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ
	real::mp_ttr,mp_tv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_laml,mp_lamr,y_face
	real,dimension(1:nof_species)::mp_d_eff,rg_difl,rg_difr,rg_enthl,rg_enthr,mp_htr,mp_hvib,rg_enthvbl,rg_enthvbr,mp_mu_i,mp_ktr_i
	real,dimension(1:dimensiona)::rg_sum_tr,rg_sumfr_tr,rg_sumfl_v,rg_sumfr_v,rg_sum_v,gradyl,gradyr,jl,jr,rg_sumfl_tr
	real,dimension(1:2)::qtr,qv
	real,dimension(1:dimensiona)::rg_sum_htr,rg_sum_hv,grady,rg_qv,rg_qtr
	real,dimension(1:nof_species)::rg_dif_av,rg_enth_av,rg_enthvb_av
	real::mp_lam_av,mp_ktr_mix_av,sum_y1,sum_y2
	real,dimension(1:dimensiona):: sumi
	real,dimension(1:dimensiona,1:nof_species) :: i_raw





	kmaxe=xmpielrank(n)
	

	


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	i=el_int(ii)
	iconsidered=i

		if (dg.eq.1) then
            dg_rhs = zero
            dg_rhs_surf_integ = zero
            dg_rhs_vol_integ = zero
            end if


                    
		   damp=lamx
            if( br2_yn == 2) damp = 0.0d0
		   if(mrf.eq.1)then
                srf=rec_mrf(i)
            end if
		    do l=1,ielem_ifca(i) !for all their faces

				  godflux2=zero
 				  angle1=ielem_faceanglex(l,i)
 				  angle2=ielem_faceangley(l,i)
 				  nx=(cos(angle1)*sin(angle2))
				  ny=(sin(angle1)*sin(angle2))
				  nz=(cos(angle2))
				  nnn(1)=nx;nnn(2)=ny;nnn(3)=nz
				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad_n
					weights_temp(1:iqp)=weights_q(1:iqp)
				  else
					iqp=qp_triangle_n
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  

				  
				  do ngp=1,iqp	!for all the gaussian quadrature points

					facex=l
					pointx=ngp

					call calculate_interior_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)





			
					
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  call rotatef(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  call rotateb(n,cleft,cleft_rot,angle1,angle2)
						  call rotateb(n,cright,cright_rot,angle1,angle2)
						  end if
						  
				      
				        leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)

				        if ((realgas.eq.1).or.(multispecies.eq.1))then

								if (multispecies.eq.1)then

								call get_visc_conduct(n,leftv,rightv,viscl,laml)

								end if


						else

								call get_visc_conduct(n,leftv,rightv,viscl,laml)

						end if
				     
							if (turbulence.eq.1)then
								if (turbulencemodel.eq.1)then
								turbmv(1)=cturbl(1);  turbmv(2)=cturbr(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
								call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
								end if
								if (turbulencemodel.eq.2)then
								eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
								eddyfl(4:6)= lcvgrad(1,1:3);eddyfl(7:9)=lcvgrad(2,1:3)
								eddyfl(10:12)=lcvgrad(3,1:3);eddyfl(13:15)=lcvgrad_t(1,1:3)
								eddyfl(16:18)=lcvgrad_t(2,1:3)


								eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
								eddyfr(4:6)= rcvgrad(1,1:3);eddyfr(7:9)=rcvgrad(2,1:3);eddyfr(10:12)=rcvgrad(3,1:3)
								eddyfr(13:15)=rcvgrad_t(1,1:3);eddyfl(16:18)=rcvgrad_t(2,1:3)
									call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
								end if
						end if
				       


						taul = zero;tau=zero;taur=zero;q=zero;ux=zero;uy=zero;uz=zero;vx=zero;vy=zero;vz=zero;wx=zero;wy=zero;wz=zero;
					    fxv=zero;fyv=zero;fzv=zero;rho12 =zero;
						u12=zero;v12=zero;w12=zero
					  

					  
					  vdamp=(4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
                                        nall(1)=nx;nall(2)=ny;nall(3)=nz


					leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)
					call cons2div(n,leftv,mp_pinfl,gammal)
					call cons2div(n,rightv,mp_pinfr,gammar)
					rho12 = oo2*(leftv(1)+rightv(1))
					rho12l=leftv(1);rho12r=rightv(1)
					u12   = oo2*(leftv(2)+rightv(2))
					  v12   = oo2*(leftv(3)+rightv(3))
					  w12   = oo2*(leftv(4)+rightv(4))

					 IF (nof_species.GT.0)THEN
					 do rg_i=1,nof_species
					  idxy = dimensiona+3 + rg_i
					  y_av(rg_i)=0.5d0*(leftv(idxy)+rightv(idxy))
					  end do
					  END IF

                    do k=1,nof_variables-1
					lcvgrad(k,1:3)=((lcvgrad(k,1:3)+rcvgrad(k,1:3))/(2.0d0))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(1:3)*(rightv(k+1)-leftv(k+1)))
					end do



				      !now compute all the temperature gradients +real gas
				      if ((realgas.eq.1).or.(multispecies.eq.1))then


								if (realgas .eq. 1) then

								!diffusion coefficient for species rg_diffl(1:nof_species),rg_diffr(1:nof_species)
								!viscosity for the mixture viscl(1)-left,viscl(2)-right
								!thermal conductivity mix laml(1)-left,laml(2)-right
								!vibrational thermal conductivity!mplaml-left,mp_lamr-right
								!species diffusion coefficients !rg_difl(1:nof_species),rg_difr(1:nof_species)
								!enthalpies for species	!rg_enthl(1:nof_species),rg_enthr(1:nof_species)
								!vibrational enthalpies for species !!rg_enthvbl(1:nof_species),rg_enthvbr(1:nof_species)



								leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)
								call multispecies_mixtures_rg(leftv,mp_mu_mix,mp_ktr_mix,mp_kve,rg_difl,rg_enthl,rg_enthvbl,gammal)

								mp_laml=mp_kve
								laml(1)=mp_ktr_mix
								viscl(1)=mp_mu_mix

								call multispecies_mixtures_rg(rightv,mp_mu_mix,mp_ktr_mix,mp_kve,rg_difr,rg_enthr,rg_enthvbr,gammar)

								mp_lamr=mp_kve
								laml(2)=mp_ktr_mix
								viscl(2)=mp_mu_mix

								mp_ktr_mix_av  = 0.5d0 * (laml(2)  + laml(1))     ! translational conductivity

								mp_lam_av  = 0.5d0 * (mp_lamr  + mp_laml)     ! ! vibrational conductivity

								rg_dif_av=0.5d0 * (rg_difr  + rg_difl)

								rg_enth_av=0.5d0 * (rg_enthr  + rg_enthl)
								rg_enthvb_av=0.5d0*(rg_enthvbr +rg_enthvbl)





											if (turbulence .eq. 1) then
											q(1:3)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:3))
											end if


													rg_sum_htr(1:dimensiona) = 0.0d0
													rg_sum_hv(1:dimensiona)  = 0.0d0



													if (rg_relax.ge.1)then

													! ---- first loop: raw fick fluxes i_k = -ρ d_k ∇y_k ----
													do rg_i = 1, nof_species

														idxy = dimensiona + 2 + rg_i     ! index of y_k at face

														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)




														! raw mixture-averaged diffusion flux
														i_raw(1:dimensiona,rg_i) = -rho12 * rg_dif_av(rg_i) * grady(1:dimensiona)





														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)
													end do


													! ---- second loop: mass-conserving flux j_k ----
													do rg_i = 1, nof_species

														! face-averaged mass fraction (must match your reconstruction)
														y_face = y_av(rg_i)

														! mass-conserving diffusion flux
														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - y_face * sumi(1:dimensiona)

														! add species diffusion fluxes into fxv / fyv
														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)

! 														! energy diffusion accumulation
 														rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
 																				+ rg_enth_av(rg_i)*jl(1:dimensiona)
 														rg_sum_hv(1:dimensiona)  = rg_sum_hv(1:dimensiona) &
 																				+ rg_enthvb_av(rg_i)*jl(1:dimensiona)
													end do
													end if

													! =============================================================
													! 4. conductive heat fluxes (fourier)
													! =============================================================
													rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona+1,1:dimensiona)
													rg_qv (1:dimensiona) = -mp_lam_av     * lcvgrad(dimensiona+2,1:dimensiona)

													! =============================================================
													! 5. total diffusive fluxes
													! =============================================================
													q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
																	+ rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

													qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

													! ------------------------------------------------
													! 6. add to conservative flux vectors
													! ------------------------------------------------
													! total energy equation (ρe)
													fxv(5) = fxv(5) + q(1)
													fyv(5) = fyv(5) + q(2)
													fzv(5) = fzv(5) + q(3)

													! vibrational energy equation (ρev)
													fxv(6) = fxv(6) + qvib(1)
													fyv(6) = fyv(6) + qvib(2)
													fzv(6) = fzv(6) + qvib(3)




								end if	!real gas ends





							if (multispecies.eq.1)then
							!viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

											if (turbulence .eq. 1) then
											q(1:3)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:3))
											else
											q(1:3) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:3))
											end if

							fxv(5) = fxv(5) - q(1);fyv(5) = fyv(5) - q(2);fzv(5) = fzv(5) - q(3)
							end if

				end if


				      if ((realgas.eq.0))then
							if (turbulence .eq. 1) then
							q(1:3)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:3))
							else
							q(1:3) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:3))
							end if

							fxv(5) = fxv(5) - q(1);fyv(5) = fyv(5) - q(2);fzv(5) = fzv(5) - q(3)

				      end if




							
					 
					  !left state derivatives
					 
					  ! determine taul!!
					   ux = lcvgrad(1,1); uy = lcvgrad(1,2); uz = lcvgrad(1,3);
					  vx = lcvgrad(2,1); vy = lcvgrad(2,2); vz = lcvgrad(2,3);
					  wx = lcvgrad(3,1); wy = lcvgrad(3,2); wz = lcvgrad(3,3);
					  
					  
					  
                                
					  
					  
					  ! tau_xx
					  taul(1,1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
					  ! tau_yy
					  taul(2,2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
					  ! tau_zz
					  taul(3,3) = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy

					  ! tau_xy
					  taul(1,2) = (uy + vx);taul(2,1) = taul(1,2)

					  ! tau_xz
					  taul(1,3) = (wx + uz);taul(3,1) = taul(1,3)

					  ! tau_yz
					  taul(2,3) = (vz + wy);taul(3,2) = taul(2,3)
					  !end determine taul
					  
					  
					  
					 ! average and multiplay by viscosity
					  if ( turbulence .eq. 1) then
					    tau = oo2*(( (viscl(1)+viscl(3))+(viscl(2)+viscl(4))))*taul
					  else
					    tau =oo2*(( (viscl(1))+(viscl(2))))*taul
					  end if

					    ! now addition into momentum fluxes
					  do kc=2,4
					      fxv(kc) = fxv(kc) + tau(1,kc-1)
					      fyv(kc) = fyv(kc) + tau(2,kc-1)
					      fzv(kc) = fzv(kc) + tau(3,kc-1)
					  enddo



					    ! compute interface velocities



					  
								   					  
					  
					  fxv(5) = fxv(5) + u12*tau(1,1) + v12*tau(1,2) + w12*tau(1,3)
					  fyv(5) = fyv(5) + u12*tau(2,1) + v12*tau(2,2) + w12*tau(2,3)
					  fzv(5) = fzv(5) + u12*tau(3,1) + v12*tau(3,2) + w12*tau(3,3)




		
		
					  hllcflux(1:nof_variables)=(nx*fxv+ny*fyv+nz*fzv)	


					  if (dg.eq.1)then
					  rhllcflux(1:nof_variables)=hllcflux(1:nof_variables)

					  dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)

					  else
					  
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))

				      end if
				      
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (turbulence.eq.1)then
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2)))+(oo2*(viscl(3)+viscl(4))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,3)+rcvgrad_t(1:turbulenceequations+passivescalar,3))*oo2*nz))
					  else
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,3)+rcvgrad_t(1:turbulenceequations+passivescalar,3))*oo2*nz))
					  
					  end if
						if (turbulencemodel.eq.1)then
						hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)/sigma
						end if					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				  end do
				  
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)-&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

 				    end if

		    end do

		    if (dg.eq.1)then




                    dg_rhs = dg_rhs_surf_integ

                    rhs_valdg(:,:,i) = rhs_valdg(:,:,i) - dg_rhs




                 end if





	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	

	
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)
	iconsidered=i	
			if (dg.eq.1) then
				dg_rhs = zero
				dg_rhs_surf_integ = zero
				dg_rhs_vol_integ = zero



				end if
				
		     if(mrf.eq.1)then
                srf=rec_mrf(i)
            end if
		     
		   
		    do l=1,ielem_ifca(i)



                                     igoflux=0




				   damp=lamx  
                    if( br2_yn == 2) damp = 0.0d0
				  angle1=ielem_faceanglex(l,i)
				  angle2=ielem_faceangley(l,i)
				   nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
 				   nnn(1)=nx;nnn(2)=ny;nnn(3)=nz
 				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad_n
					weights_temp(1:iqp)=weights_q(1:iqp)
				  else
					iqp=qp_triangle_n
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  

				  
				  godflux2=zero
				  
				  b_code=0
								  
				  
				  do ngp=1,iqp
				  facex=l
				  pointx=ngp

				  call calculate_bounded_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)




				      
				    
			
					
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  call rotatef(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  call rotateb(n,cleft,cleft_rot,angle1,angle2)
						  call rotateb(n,cright,cright_rot,angle1,angle2)
						  end if
						  
				      
				        leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)


				        if ((realgas.eq.1).or.(multispecies.eq.1))then

								if (multispecies.eq.1)then

								call get_visc_conduct(n,leftv,rightv,viscl,laml)

								end if


						else

					call get_visc_conduct(n,leftv,rightv,viscl,laml)

						end if



				     
					    if (turbulence.eq.1)then
						      if (turbulencemodel.eq.1)then
							  turbmv(1)=cturbl(1);  turbmv(2)=cturbr(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
						      if (turbulencemodel.eq.2)then
							  eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
							  eddyfl(4:6)= lcvgrad(1,1:3);eddyfl(7:9)=lcvgrad(2,1:3)
							  eddyfl(10:12)=lcvgrad(3,1:3);eddyfl(13:15)=lcvgrad_t(1,1:3)
							  eddyfl(16:18)=lcvgrad_t(2,1:3)
							    
							    
							  eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
							  eddyfr(4:6)= rcvgrad(1,1:3);eddyfr(7:9)=rcvgrad(2,1:3);eddyfr(10:12)=rcvgrad(3,1:3)
							  eddyfr(13:15)=rcvgrad_t(1,1:3);eddyfl(16:18)=rcvgrad_t(2,1:3)
							    call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
					  end if


                                            
					  taul = zero;tau=zero;taur=zero;q=zero;ux=zero;uy=zero;uz=zero;vx=zero;vy=zero;vz=zero;wx=zero;wy=zero;wz=zero;
					    fxv=zero;fyv=zero;fzv=zero;rho12 =zero;
					  u12=zero;v12=zero;w12=zero  
					  

					 if ((b_code.gt.0))then
					  damp=zero
 					  end if




					  
					  vdamp=(4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
                                        nall(1)=nx;nall(2)=ny;nall(3)=nz



                                              leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)
					call cons2div(n,leftv,mp_pinfl,gammal)
					call cons2div(n,rightv,mp_pinfr,gammar)
					rho12 = oo2*(leftv(1)+rightv(1))
					rho12l=leftv(1);rho12r=rightv(1)
					u12   = oo2*(leftv(2)+rightv(2))
					  v12   = oo2*(leftv(3)+rightv(3))



					 if (nof_species.gt.0)then
					  do rg_i=1,nof_species
					  idxy = dimensiona+3 + rg_i
					  y_av(rg_i)=0.5d0*(leftv(idxy)+rightv(idxy))
					  end do
					  end if



                                         do k=1,nof_variables-1
											lcvgrad(k,1:3)=((lcvgrad(k,1:3)+rcvgrad(k,1:3))/(2.0d0))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(1:3)*(rightv(k+1)-leftv(k+1)))
										 end do


					!now compute all the temperature gradients +real gas
				      if ((realgas.eq.1).or.(multispecies.eq.1))then

								if (realgas.eq.1)then

								!diffusion coefficient for species rg_diffl(1:nof_species),rg_diffr(1:nof_species)
								!viscosity for the mixture viscl(1)-left,viscl(2)-right
								!thermal conductivity mix laml(1)-left,laml(2)-right
								!vibrational thermal conductivity!mplaml-left,mp_lamr-right
								!species diffusion coefficients !rg_difl(1:nof_species),rg_difr(1:nof_species)
								!enthalpies for species	!rg_enthl(1:nof_species),rg_enthr(1:nof_species)
								!vibrational enthalpies for species !!rg_enthvbl(1:nof_species),rg_enthvbr(1:nof_species)



								leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)
								call multispecies_mixtures_rg(leftv,mp_mu_mix,mp_ktr_mix,mp_kve,rg_difl,rg_enthl,rg_enthvbl,gammal)

								mp_laml=mp_kve
								laml(1)=mp_ktr_mix
								viscl(1)=mp_mu_mix

								call multispecies_mixtures_rg(rightv,mp_mu_mix,mp_ktr_mix,mp_kve,rg_difr,rg_enthr,rg_enthvbr,gammar)

								mp_lamr=mp_kve
								laml(2)=mp_ktr_mix
								viscl(2)=mp_mu_mix

								mp_ktr_mix_av = 0.5d0*(laml(1) + laml(2))
													mp_lam_av     = 0.5d0*(mp_laml + mp_lamr)

													! mixture-averaged species diffusion coefficients
													rg_dif_av(:) = 0.5d0*(rg_difl(:) + rg_difr(:))

													! face-averaged species enthalpies
													rg_enth_av(:)   = 0.5d0*(rg_enthl(:)   + rg_enthr(:))
													rg_enthvb_av(:) = 0.5d0*(rg_enthvbl(:) + rg_enthvbr(:))





											if (turbulence .eq. 1) then
											q(1:3)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:3))
											end if


													rg_sum_htr(1:dimensiona) = 0.0d0
													rg_sum_hv(1:dimensiona)  = 0.0d0

													! =============================================================
													! 3. mixture-averaged diffusion fluxes (no pressure terms)
													! =============================================================
													rg_sum_htr(1:dimensiona) = 0.0d0
													rg_sum_hv(1:dimensiona)  = 0.0d0

													sumi(1:dimensiona) = 0.0d0    ! σ i_k


													if (rg_relax.ge.1)then

													! ---- first loop: raw fick fluxes i_k = -ρ d_k ∇y_k ----
													do rg_i = 1, nof_species

														idxy = dimensiona + 2 + rg_i     ! index of y_k at face

														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)

														! raw mixture-averaged diffusion flux
														i_raw(1:dimensiona,rg_i) = -rho12 * rg_dif_av(rg_i) * grady(1:dimensiona)



														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)
													end do





													! ---- second loop: mass-conserving flux j_k ----
													do rg_i = 1, nof_species

														! face-averaged mass fraction (must match your reconstruction)
														y_face = y_av(rg_i)

														! mass-conserving diffusion flux
														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - y_face * sumi(1:dimensiona)



														! add species diffusion fluxes into fxv / fyv
														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)




														! energy diffusion accumulation
														rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
																				+ rg_enth_av(rg_i)*jl(1:dimensiona)
														rg_sum_hv(1:dimensiona)  = rg_sum_hv(1:dimensiona) &
																				+ rg_enthvb_av(rg_i)*jl(1:dimensiona)
													end do



													end if



													! =============================================================
													! 4. conductive heat fluxes (fourier)
													! =============================================================
													rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona+1,1:dimensiona)
													rg_qv (1:dimensiona) = -mp_lam_av     * lcvgrad(dimensiona+2,1:dimensiona)

													! =============================================================
													! 5. total diffusive fluxes
													! =============================================================
													q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
																	+ rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

													qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

													! ------------------------------------------------
													! 6. add to conservative flux vectors
													! ------------------------------------------------
													! total energy equation (ρe)
													fxv(5) = fxv(5) - q(1)
													fyv(5) = fyv(5) - q(2)
													fzv(5) = fyv(5) - q(3)

													! vibrational energy equation (ρev)
													fxv(6) = fxv(6) + qvib(1)
													fyv(6) = fyv(6) + qvib(2)
													fzv(6) = fzv(6) + qvib(3)




								end if	!real gas ends





							if (multispecies.eq.1)then
							!viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

											if (turbulence .eq. 1) then
											q(1:3)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:3))
											else
											q(1:3) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:3))
											end if

							fxv(5) = fxv(5) - q(1);fyv(5) = fyv(5) - q(2);fzv(5) = fzv(5) - q(3)
							end if




				      end if
				       if (realgas.eq.0)then
							if (turbulence .eq. 1) then
							q(1:3)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:3))
							else
							q(1:3) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:3))
							end if

							fxv(5) = fxv(5) - q(1);fyv(5) = fyv(5) - q(2);fzv(5) = fzv(5) - q(3)

				      end if
				       				  

					  

							
					 
					  !left state derivatives
					 
					  ! determine taul!!
					   ux = lcvgrad(1,1); uy = lcvgrad(1,2); uz = lcvgrad(1,3);
					  vx = lcvgrad(2,1); vy = lcvgrad(2,2); vz = lcvgrad(2,3);
					  wx = lcvgrad(3,1); wy = lcvgrad(3,2); wz = lcvgrad(3,3);
					  
					  
					  
                                
					  
					  
					  ! tau_xx
					  taul(1,1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
					  ! tau_yy
					  taul(2,2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
					  ! tau_zz
					  taul(3,3) = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy

					  ! tau_xy
					  taul(1,2) = (uy + vx);taul(2,1) = taul(1,2)

					  ! tau_xz
					  taul(1,3) = (wx + uz);taul(3,1) = taul(1,3)

					  ! tau_yz
					  taul(2,3) = (vz + wy);taul(3,2) = taul(2,3)
					  !end determine taul
					  
					  
					  
					 ! average and multiplay by viscosity
					  if ( turbulence .eq. 1) then
					    tau = oo2*(( (viscl(1)+viscl(3))+(viscl(2)+viscl(4))))*taul
					  else
					    tau =oo2*(( (viscl(1))+(viscl(2))))*taul
					  end if

					    ! now addition into momentum fluxes
					  do kc=2,4
					      fxv(kc) = fxv(kc) + tau(1,kc-1)
					      fyv(kc) = fyv(kc) + tau(2,kc-1)
					      fzv(kc) = fzv(kc) + tau(3,kc-1)
					  enddo



					  
								   					  
					  
					  fxv(5) = fxv(5) + u12*tau(1,1) + v12*tau(1,2) + w12*tau(1,3)
					  fyv(5) = fyv(5) + u12*tau(2,1) + v12*tau(2,2) + w12*tau(2,3)
					  fzv(5) = fzv(5) + u12*tau(3,1) + v12*tau(3,2) + w12*tau(3,3)




		
		
					  hllcflux(1:nof_variables)=(nx*fxv+ny*fyv+nz*fzv)	


					  if (realgas.eq.1)then
						if ((b_code.eq.4).and.(catalytic_wall.eq.0))then
						hllcflux(7:nof_variables)=zero
						end if
						end if


					   if (dg.eq.1)then
					   rhllcflux(1:nof_variables)=hllcflux(1:nof_variables)

					  dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)


					  else
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					   if (turbulence.eq.1)then
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2)))+(oo2*(viscl(3)+viscl(4))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,3)+rcvgrad_t(1:turbulenceequations+passivescalar,3))*oo2*nz))
					  else
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,3)+rcvgrad_t(1:turbulenceequations+passivescalar,3))*oo2*nz))
					  
					  end if
						if (turbulencemodel.eq.1)then
						hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)/sigma
						end if					  

					if ((b_code.eq.3))then
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
					  
					  end if




					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				  end do
				  
				    
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)-&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

				    end if
		    end do

		     if (dg.eq.1)then



                    dg_rhs = dg_rhs_surf_integ

                    rhs_valdg(:,:,i) = rhs_valdg(:,:,i) - dg_rhs






                 end if



	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end subroutine calculate_fluxeshi_diffusive




subroutine calculate_fluxeshi_diffusive2d(n)
!> @brief
!> this subroutine computes the diffusive fluxes for euler-navier-stokes equations in 2d
	implicit none
	integer,intent(in)::n
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ii,nvar,kc,iex,ittt,ikas,igoflux, icaseb,kk,b_code
	real::sum_detect,norms
	integer::iconsidered,facex,pointx,k,rg_i,rg_j,idxy,icompute
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:nof_variables)::leftv,rightv,srf_speedrot,tempx_l,rtempx_l
	real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
	real,dimension(1:dimensiona)::pox,poy,poz,rg_sumfl,rg_sumfr
	real,dimension(1:nof_variables)::srf_speed
	real,dimension(1:numberofpoints2)::weights_temp
	real,dimension(1:8,1:dimensiona)::vext
	real::mp_pinfl,mp_pinfr,gammal,gammar,rho12l,rho12r
	real,dimension(1:4)::viscl,laml
	real,dimension(1:nof_species)::rgs_htr_i, rgs_hvib_i,y_av,x_av
	real,dimension(1:2)::turbmv
    real,dimension(1)::etvm
    real,dimension(1:20)::eddyfl,eddyfr
    real,dimension(1:nof_variables-1,1:dims)::lcvgrad,rcvgrad,jtmp
	real,dimension(turbulenceequations+passivescalar,1:dims)::lcvgrad_t,rcvgrad_t
	real,dimension(1:nof_variables)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(2,2)::taul,taur,tau
	real,dimension(2)::q,nall,qvib
	real::ux,uy,uz,vx,vy,vz,wx,wy,wz,rho12,u12,v12,w12,damp,vdamp,y_face
	real,dimension(1:idegfree+1,1:nof_variables)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ
	real::mp_ttr,mp_tv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_laml,mp_lamr
	real,dimension(1:nof_species)::mp_d_eff,rg_difl,rg_difr,rg_enthl,rg_enthr,mp_htr,mp_hvib,rg_enthvbl,rg_enthvbr,mp_mu_i,mp_ktr_i
	real,dimension(1:dimensiona)::rg_sum_tr,rg_sumfr_tr,rg_sumfl_tr,rg_sumfl_v,rg_sumfr_v,rg_sum_v,gradyl,gradyr,jl,jr
	real,dimension(1:2)::qtr,qv
	real,dimension(1:dimensiona)::rg_sum_htr,rg_sum_hv,grady,rg_qv,rg_qtr
	real,dimension(1:nof_species)::rg_dif_av,rg_enth_av,rg_enthvb_av
	real::mp_lam_av,mp_ktr_mix_av,sum_y1,sum_y2
	real,dimension(1:dimensiona):: sumi,gradx,grad_a
	real,dimension(1:dimensiona,1:nof_species) :: i_raw
	real :: molar_sum, mbar


    
	
	
	
	kmaxe=xmpielrank(n)
	

	
	weights_temp(1:qp_line_n)=weights_l(1:qp_line_n)

	


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	i=el_int(ii)
	iconsidered=i
		   
			if (dg.eq.1) then
				dg_rhs = zero
				dg_rhs_surf_integ = zero
				dg_rhs_vol_integ = zero



				end if



		    do l=1,ielem_ifca(i) !for all their faces
!                 if (ielem_reorient(l,i).eq.0)then
                    damp=lamx
                    if( br2_yn == 2) damp = 0.0d0
                    
				  godflux2=zero
 				  angle1=ielem_faceanglex(l,i)
 				  angle2=ielem_faceangley(l,i)
 				  nx=angle1
				  ny=angle2
				 
				  
					iqp=qp_line_n
				
				  
				  do ngp=1,iqp	!for all the gaussian quadrature points
				  facex=l
				  pointx=ngp

				  call calculate_interior_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)




			
					
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  call rotatef2d(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef2d(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht2d(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
						  call rotateb2d(n,cright,cright_rot,angle1,angle2)
						  end if
						  
				      
				        leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)


				        if ((realgas.eq.1).or.(multispecies.eq.1))then

								if (multispecies.eq.1)then

								call get_visc_conduct(n,leftv,rightv,viscl,laml)

								end if


				        else

								call get_visc_conduct(n,leftv,rightv,viscl,laml)
						end if
				     
					    if (turbulence.eq.1)then
						      if (turbulencemodel.eq.1)then
							  turbmv(1)=cturbl(1);  turbmv(2)=cturbr(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
						      if (turbulencemodel.eq.2)then
							  eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
							  eddyfl(4:5)= lcvgrad(1,1:2);eddyfl(6:7)=lcvgrad(2,1:2)
							  eddyfl(8:9)=lcvgrad_t(1,1:2)
							  eddyfl(10:11)=lcvgrad_t(2,1:2)
							    
							    
							  eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
							  eddyfr(4:5)= rcvgrad(1,1:2);eddyfr(6:7)=rcvgrad(2,1:2)
							  eddyfr(8:9)=rcvgrad_t(1,1:2);eddyfr(10:11)=rcvgrad_t(2,1:2)
							    call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
					  end if
				       
				       



                                     taul = zero;tau=zero;taur=zero;q=zero;ux=zero;uy=zero;uz=zero;vx=zero;vy=zero;vz=zero;wx=zero;wy=zero;wz=zero;
					    fxv=zero;fyv=zero;fzv=zero;rho12 =zero;
					  u12=zero;v12=zero;w12=zero 
				       
! 				      
					  
					  vdamp=(4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
                                        nall(1)=nx;nall(2)=ny



                               leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)
					call cons2div(n,leftv,mp_pinfl,gammal)
					call cons2div(n,rightv,mp_pinfr,gammar)
					rho12 = oo2*(leftv(1)+rightv(1))
					rho12l=leftv(1);rho12r=rightv(1)
					u12   = oo2*(leftv(2)+rightv(2))
					  v12   = oo2*(leftv(3)+rightv(3))


					  if (nof_species.GT.0)then
					  do rg_i=1,nof_species
					  idxy = dimensiona+3 + rg_i
					  y_av(rg_i)=0.5d0*(leftv(idxy)+rightv(idxy))




					  end do
					  end if





					do k=1,nof_variables-1
					lcvgrad(k,1:2)=((lcvgrad(k,1:2)+rcvgrad(k,1:2))/(2.0d0))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(1:2)*(rightv(k+1)-leftv(k+1)))
					end do



					  				

                           		    !now compute all the temperature gradients +real gas
				      if ((realgas.eq.1).or.(multispecies.eq.1))then


												if (realgas .eq. 1) then

													! =============================================================
													! 1. compute mixture transport properties (left & right)
													! =============================================================
													leftv(1:nof_variables)  = cleft(1:nof_variables)
													call multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, &
																				rg_difl, rg_enthl, rg_enthvbl, gammal)


													laml(1)   = mp_ktr_mix

													mp_laml   = mp_kve
													viscl(1)=mp_mu_mix


													rightv(1:nof_variables) = cright(1:nof_variables)
													call multispecies_mixtures_rg(rightv, mp_mu_mix, mp_ktr_mix, mp_kve, &
																				rg_difr, rg_enthr, rg_enthvbr, gammar)

													! mixture thermal conductivities (average)

													laml(2)   = mp_ktr_mix
													mp_lamr   = mp_kve
													viscl(2)=mp_mu_mix


													mp_ktr_mix_av = 0.5d0*(laml(1) + laml(2))
													mp_lam_av     = 0.5d0*(mp_laml + mp_lamr)

													! mixture-averaged species diffusion coefficients
													rg_dif_av(:) = 0.5d0*(rg_difl(:) + rg_difr(:))

													! face-averaged species enthalpies
													rg_enth_av(:)   = 0.5d0*(rg_enthl(:)   + rg_enthr(:))
													rg_enthvb_av(:) = 0.5d0*(rg_enthvbl(:) + rg_enthvbr(:))



													! =============================================================
													! 2. turbulence shortcut
													! =============================================================
													if (turbulence .eq. 1) then
														q(1:2) = -0.5d0*(laml(3)+laml(4))*lcvgrad(dimensiona+1,1:2)

													end if

													! =============================================================
													! 3. mixture-averaged diffusion fluxes (no pressure terms)
													! =============================================================


													rg_sum_htr(1:dimensiona) = 0.0d0
													rg_sum_hv(1:dimensiona)  = 0.0d0
													rg_qtr=0.0d0; rg_qv=0.0d0

													sumi(1:dimensiona) = 0.0d0    ! σ i_k

													if (rg_relax.ge.1)then

													! ---- first loop: raw fick fluxes i_k = -ρ d_k ∇y_k ----
! 													do rg_i = 1, nof_species
!
! 														idxy = dimensiona + 2 + rg_i     ! index of y_k at face
!
! 														!u,v,ttr,tvib,
!
! 														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)
!
!
!
!
! 														! raw mixture-averaged diffusion flux
! 														i_raw(1:dimensiona,rg_i) = -rho12 * rg_dif_av(rg_i) * grady(1:dimensiona)
!
!
!
!
!
! 														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)
! 													end do

													sumi(1:dimensiona) = 0.0d0
													i_raw(1:dimensiona,1:nof_species) = 0.0d0

													! compute face mole fractions from face mass fractions
													molar_sum = 0.0d0
													do rg_j = 1, nof_species
														molar_sum = molar_sum + max(y_av(rg_j),0.0d0) / rg_molm(rg_j)
													end do

													molar_sum = max(molar_sum,1.0d-30)
													mbar = 1.0d0 / molar_sum

													do rg_j = 1, nof_species
														x_av(rg_j) = (max(y_av(rg_j),0.0d0) / rg_molm(rg_j)) / molar_sum
													end do

													! grad_a = grad(sum_k Y_k/M_k)
													grad_a(1:dimensiona) = 0.0d0

													do rg_j = 1, nof_species
														idxy = dimensiona + 2 + rg_j
														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)

														grad_a(1:dimensiona) = grad_a(1:dimensiona) + &
															grady(1:dimensiona) / rg_molm(rg_j)
													end do

													! raw mixture-averaged flux:
													! J_i_raw = -rho * D_i * (M_i/Mmix) * grad(X_i)
													do rg_i = 1, nof_species

														idxy = dimensiona + 2 + rg_i
														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)

														gradx(1:dimensiona) = &
															( grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
															- (max(y_av(rg_i),0.0d0) / rg_molm(rg_i)) * grad_a(1:dimensiona) ) &
															/ (molar_sum * molar_sum)

														i_raw(1:dimensiona,rg_i) = &
															-rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)

														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)

													end do


													! ---- second loop: mass-conserving flux j_k ----
													do rg_i = 1, nof_species

														! face-averaged mass fraction (must match your reconstruction)
														y_face = y_av(rg_i)

														! mass-conserving diffusion flux
														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - y_face * sumi(1:dimensiona)

														! add species diffusion fluxes into fxv / fyv
														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)

! 														! energy diffusion accumulation
 														rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
 																				+ rg_enth_av(rg_i)*jl(1:dimensiona)
 														rg_sum_hv(1:dimensiona)  = rg_sum_hv(1:dimensiona) &
 																				+ rg_enthvb_av(rg_i)*jl(1:dimensiona)
													end do


													! =============================================================
													! 4. conductive heat fluxes (fourier)
													! =============================================================
													rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona+1,1:dimensiona)
													rg_qv (1:dimensiona) = -mp_lam_av     * lcvgrad(dimensiona+2,1:dimensiona)


													end if

													! =============================================================
													! 5. total diffusive fluxes
													! =============================================================
													q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
																	+ rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

													qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

													! =============================================================
													! 6. add to conservative flux vectors
													! =============================================================
													fxv(4) = fxv(4) - q(1)
													fyv(4) = fyv(4) - q(2)

													fxv(5) = fxv(5) - qvib(1)
													fyv(5) = fyv(5) - qvib(2)






!

											end if	!real gas ends





									if (multispecies.eq.1)then
									!viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

													if (turbulence .eq. 1) then
													q(1:2)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:2))
													else
													q(1:2) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:2))
													end if

									fxv(4) = fxv(4) - q(1);fyv(4) = fyv(4) - q(2)
									end if




				      end if

				      if (realgas.eq.0)then

											if (turbulence .eq. 1) then
											q(1:2)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:2))
											else
											q(1:2) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:2))
											end if
											!q(1:2) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:2))
											fxv(4) = fxv(4) - q(1);fyv(4) = fyv(4) - q(2)

				      end if





					  
					  
					  
					   
				       				  

					 
					  !left state derivatives
					  ux = lcvgrad(1,1); uy = lcvgrad(1,2)
					  vx = lcvgrad(2,1); vy = lcvgrad(2,2)
					  ! determine taul!!
					 

					  ! tau_xx
					  taul(1,1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy 
					  ! tau_yy
					  taul(2,2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux 
					  ! tau_zz
					 

					  ! tau_xy
					  taul(1,2) = (uy + vx);taul(2,1) = taul(1,2)

					  
					 ! average and multiplay by viscosity
					  if ( turbulence .eq. 1) then
					    tau = oo2*(( (viscl(1)+viscl(3)))+( (viscl(2)+viscl(4))))*taul
					  else
					    tau = oo2*((viscl(1))+(viscl(2)))*taul
					  end if

					    ! now addition into momentum fluxes
					  do kc=2,3
					      fxv(kc) = fxv(kc) + tau(1,kc-1)
					      fyv(kc) = fyv(kc) + tau(2,kc-1)
					      
					  enddo





					  fxv(4) = fxv(4) + u12*tau(1,1) + v12*tau(1,2) 
					  fyv(4) = fyv(4) + u12*tau(2,1) + v12*tau(2,2)








! 					 
					  hllcflux(1:nof_variables)=(nx*fxv+ny*fyv)	
		
					  if (dg.eq.1)then

					  rhllcflux(1:nof_variables)=hllcflux(1:nof_variables)
						
						dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
  
						else
					  
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
						end if
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (turbulence.eq.1)then
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2)))+(oo2*(viscl(3)+viscl(4))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny))
					  else
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny))
					  
					  end if
						if (turbulencemodel.eq.1)then
						hllcflux(5)=hllcflux(5)/sigma
						end if					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				  end do
				  
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)-&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

 				    end if
		    end do

			if (dg.eq.1)then



				dg_rhs = dg_rhs_surf_integ

				

				rhs_valdg(:,:,i) = rhs_valdg(:,:,i) - dg_rhs






			 end if



	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	


	
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)
	iconsidered=i	
			if (dg.eq.1) then
				dg_rhs = zero
				dg_rhs_surf_integ = zero
				dg_rhs_vol_integ = zero



				end if
		    
		    
		    do l=1,ielem_ifca(i)
                damp=lamx
                if( br2_yn == 2) damp = 0.0d0
				      b_code=0
				 godflux2=zero
 				  angle1=ielem_faceanglex(l,i)
 				  angle2=ielem_faceangley(l,i)
 				  nx=angle1
				  ny=angle2
				 
				  
					iqp=qp_line_n
				
				  
				 
								  
				  
				  do ngp=1,iqp



				   facex=l
				   pointx=ngp

				   call calculate_bounded_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)




				      
				    
			
					
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  call rotatef2d(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef2d(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht2d(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
						  call rotateb2d(n,cright,cright_rot,angle1,angle2)
						  end if
						  
				      
				        leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)

						if ((realgas.eq.1).or.(multispecies.eq.1))then

								if (multispecies.eq.1)then

								call get_visc_conduct(n,leftv,rightv,viscl,laml)

								end if


				        else

								call get_visc_conduct(n,leftv,rightv,viscl,laml)
						end if
				     
					    if (turbulence.eq.1)then
						      if (turbulencemodel.eq.1)then
							  turbmv(1)=cturbl(1);  turbmv(2)=cturbr(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      if (b_code.eq.4)then
						      viscl(3:4)=zero
						      laml(3:4)=zero
						      end if
						      end if
						      if (turbulencemodel.eq.2)then
							  eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
							  eddyfl(4:5)= lcvgrad(1,1:2);eddyfl(6:7)=lcvgrad(2,1:2)
							  eddyfl(8:9)=lcvgrad_t(1,1:2)
							  eddyfl(10:11)=lcvgrad_t(2,1:2)
							    
							    
							  eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
							  eddyfr(4:5)= rcvgrad(1,1:2);eddyfr(6:7)=rcvgrad(2,1:2)
							  eddyfr(8:9)=rcvgrad_t(1,1:2);eddyfr(10:11)=rcvgrad_t(2,1:2)
							    call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
					  end if
				       
				       
				       
				       



                                          taul = zero;tau=zero;taur=zero;q=zero;ux=zero;uy=zero;uz=zero;vx=zero;vy=zero;vz=zero;wx=zero;wy=zero;wz=zero;
					    fxv=zero;fyv=zero;fzv=zero;rho12 =zero;
					  u12=zero;v12=zero;w12=zero 
				       
				      if ((b_code.gt.0))then
					  damp=zero
 					  end if




				       
					   vdamp=(4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
                                        nall(1)=nx;nall(2)=ny


                                              leftv(1:nof_variables)=cleft(1:nof_variables);rightv(1:nof_variables)=cright(1:nof_variables)
					call cons2div(n,leftv,mp_pinfl,gammal)
					call cons2div(n,rightv,mp_pinfr,gammar)
					rho12 = oo2*(leftv(1)+rightv(1))
					rho12l=leftv(1);rho12r=rightv(1)
					u12   = oo2*(leftv(2)+rightv(2))
					  v12   = oo2*(leftv(3)+rightv(3))

					  if (nof_species.GT.0)then
					  do rg_i=1,nof_species
					  idxy = dimensiona+3 + rg_i
					  y_av(rg_i)=0.5d0*(leftv(idxy)+rightv(idxy))




					  end do
					  end if






					do k=1,nof_variables-1
					lcvgrad(k,1:2)=((lcvgrad(k,1:2)+rcvgrad(k,1:2))/(2.0d0))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(1:2)*(rightv(k+1)-leftv(k+1)))
					end do


					icompute=0


                           		    !now compute all the temperature gradients +real gas
				      if ((realgas.eq.1).or.(multispecies.eq.1))then


												if (realgas .eq. 1) then

													! =============================================================
													! 1. compute mixture transport properties (left & right)
													! =============================================================
													leftv(1:nof_variables)  = cleft(1:nof_variables)
													call multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, &
																				rg_difl, rg_enthl, rg_enthvbl, gammal)


													laml(1)   = mp_ktr_mix

													mp_laml   = mp_kve
													viscl(1)=mp_mu_mix


													rightv(1:nof_variables) = cright(1:nof_variables)
													call multispecies_mixtures_rg(rightv, mp_mu_mix, mp_ktr_mix, mp_kve, &
																				rg_difr, rg_enthr, rg_enthvbr, gammar)

													! mixture thermal conductivities (average)

													laml(2)   = mp_ktr_mix
													mp_lamr   = mp_kve
													viscl(2)=mp_mu_mix


													mp_ktr_mix_av = 0.5d0*(laml(1) + laml(2))
													mp_lam_av     = 0.5d0*(mp_laml + mp_lamr)

													! mixture-averaged species diffusion coefficients
													rg_dif_av(:) = 0.5d0*(rg_difl(:) + rg_difr(:))

													! face-averaged species enthalpies
													rg_enth_av(:)   = 0.5d0*(rg_enthl(:)   + rg_enthr(:))
													rg_enthvb_av(:) = 0.5d0*(rg_enthvbl(:) + rg_enthvbr(:))











													! =============================================================
													! 2. turbulence shortcut
													! =============================================================
													if (turbulence .eq. 1) then
														q(1:2) = -0.5d0*(laml(3)+laml(4))*lcvgrad(dimensiona+1,1:2)

													end if

													! =============================================================
													! 3. mixture-averaged diffusion fluxes (no pressure terms)
													! =============================================================
													rg_sum_htr(1:dimensiona) = 0.0d0
													rg_sum_hv(1:dimensiona)  = 0.0d0
													rg_qtr=0.0d0; rg_qv=0.0d0

													sumi(1:dimensiona) = 0.0d0    ! σ i_k

													if (rg_relax.ge.1)then
															if ((b_code.eq.4).and.(catalytic_wall.eq.0))then
																icompute=1
															end if
! 															if ((b_code.eq.1).or.(b_code.eq.2))then
! 																icompute=1
! 															end if

													if (icompute.eq.0)then

													! ---- first loop: raw fick fluxes i_k = -ρ d_k ∇y_k ----
! 													do rg_i = 1, nof_species
!
! 														idxy = dimensiona + 2 + rg_i     ! index of y_k at face
!
! 														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)
!
! 														! raw mixture-averaged diffusion flux
! 														i_raw(1:dimensiona,rg_i) = -rho12 * rg_dif_av(rg_i) * grady(1:dimensiona)
!
!
!
! 														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)
! 													end do

													sumi(1:dimensiona) = 0.0d0
													i_raw(1:dimensiona,1:nof_species) = 0.0d0

													! compute face mole fractions from face mass fractions
													molar_sum = 0.0d0
													do rg_j = 1, nof_species
														molar_sum = molar_sum + max(y_av(rg_j),0.0d0) / rg_molm(rg_j)
													end do

													molar_sum = max(molar_sum,1.0d-30)
													mbar = 1.0d0 / molar_sum

													do rg_j = 1, nof_species
														x_av(rg_j) = (max(y_av(rg_j),0.0d0) / rg_molm(rg_j)) / molar_sum
													end do

													! grad_a = grad(sum_k Y_k/M_k)
													grad_a(1:dimensiona) = 0.0d0

													do rg_j = 1, nof_species
														idxy = dimensiona + 2 + rg_j
														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)

														grad_a(1:dimensiona) = grad_a(1:dimensiona) + &
															grady(1:dimensiona) / rg_molm(rg_j)
													end do

													! raw mixture-averaged flux:
													! J_i_raw = -rho * D_i * (M_i/Mmix) * grad(X_i)
													do rg_i = 1, nof_species

														idxy = dimensiona + 2 + rg_i
														grady(1:dimensiona) = lcvgrad(idxy,1:dimensiona)

														gradx(1:dimensiona) = &
															( grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
															- (max(y_av(rg_i),0.0d0) / rg_molm(rg_i)) * grad_a(1:dimensiona) ) &
															/ (molar_sum * molar_sum)

														i_raw(1:dimensiona,rg_i) = &
															-rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)

														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)

													end do






													! ---- second loop: mass-conserving flux j_k ----
													do rg_i = 1, nof_species

														! face-averaged mass fraction (must match your reconstruction)
														y_face = y_av(rg_i)

														! mass-conserving diffusion flux
														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - y_face * sumi(1:dimensiona)



														! add species diffusion fluxes into fxv / fyv
														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)




														! energy diffusion accumulation
														rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
																				+ rg_enth_av(rg_i)*jl(1:dimensiona)
														rg_sum_hv(1:dimensiona)  = rg_sum_hv(1:dimensiona) &
																				+ rg_enthvb_av(rg_i)*jl(1:dimensiona)
													end do


													end if
! 													end if



													! =============================================================
													! 4. conductive heat fluxes (fourier)
													! =============================================================
													rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona+1,1:dimensiona)
													rg_qv (1:dimensiona) = -mp_lam_av     * lcvgrad(dimensiona+2,1:dimensiona)


													end if

													! =============================================================
													! 5. total diffusive fluxes
													! =============================================================
													q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
																	+ rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

													qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

													! =============================================================
													! 6. add to conservative flux vectors
													! =============================================================
													fxv(4) = fxv(4) - q(1)
													fyv(4) = fyv(4) - q(2)

													fxv(5) = fxv(5) - qvib(1)
													fyv(5) = fyv(5) - qvib(2)









!

											end if	!real gas ends






												if (multispecies.eq.1)then
												!viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

																if (turbulence .eq. 1) then
																q(1:2)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:2))
																else
																q(1:2) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:2))
																end if

												fxv(4) = fxv(4) - q(1);fyv(4) = fyv(4) - q(2)
												end if




				      end if


				      if (realgas.eq.0)then
							if (turbulence .eq. 1) then
											q(1:2)=  - oo2* ((laml(3)+ (laml(4)))*lcvgrad(dimensiona+1,1:2))
											else
											q(1:2) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:2))
											end if

							fxv(4) = fxv(4) - q(1);fyv(4) = fyv(4) - q(2)

				      end if




							
					 
					  !left state derivatives
					  ux = lcvgrad(1,1); uy = lcvgrad(1,2)
					  vx = lcvgrad(2,1); vy = lcvgrad(2,2)
					  ! determine taul!!
					 

					  ! tau_xx
					  taul(1,1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy 
					  ! tau_yy
					  taul(2,2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux 
					  ! tau_zz
					 

					  ! tau_xy
					  taul(1,2) = (uy + vx);taul(2,1) = taul(1,2)

					  
					 ! average and multiplay by viscosity
					  if ( turbulence .eq. 1) then
					    tau = oo2*(( (viscl(1)+viscl(3)))+( (viscl(2)+viscl(4))))*taul
					  else
					    tau = oo2*((viscl(1))+(viscl(2)))*taul
					  end if

					    ! now addition into momentum fluxes
					  do kc=2,3
					      fxv(kc) = fxv(kc) + tau(1,kc-1)
					      fyv(kc) = fyv(kc) + tau(2,kc-1)
					      
					  enddo


					  

					  fxv(4) = fxv(4) + u12*tau(1,1) + v12*tau(1,2) 
					  fyv(4) = fyv(4) + u12*tau(2,1) + v12*tau(2,2) 





! 					 
					  hllcflux(1:nof_variables)=(nx*fxv+ny*fyv)			

						if (realgas.eq.1)then
						if ((b_code.eq.4).and.(catalytic_wall.eq.0))then
						hllcflux(6:nof_variables)=zero
						end if
						end if










					  
					  if (dg.eq.1)then

					   rhllcflux(1:nof_variables)=hllcflux(1:nof_variables)

						dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)





  
						else
		




			      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
						end if
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (turbulence.eq.1)then
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2)))+(oo2*(viscl(3)+viscl(4))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny))
					  else
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar) =&
					  ((oo2*(viscl(1)+viscl(2))))*&
		    (((lcvgrad_t(1:turbulenceequations+passivescalar,1)+rcvgrad_t(1:turbulenceequations+passivescalar,1))*oo2*nx)+&
		    ((lcvgrad_t(1:turbulenceequations+passivescalar,2)+rcvgrad_t(1:turbulenceequations+passivescalar,2))*oo2*ny))
					  
					  end if
						if (turbulencemodel.eq.1)then
						hllcflux(5)=hllcflux(5)/sigma
						end if			
						
						if ((b_code.eq.3))then
					  hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
					  
					  end if
						
						
						
						
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				  end do
				  
				    
				     rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)
				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)-&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

				    end if
		    end do

			if (dg.eq.1)then



				dg_rhs = dg_rhs_surf_integ


				rhs_valdg(:,:,i) = rhs_valdg(:,:,i) - dg_rhs






			 end if



	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	



end subroutine calculate_fluxeshi_diffusive2d




subroutine calculate_fluxeshi_convective_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
	implicit none
	integer,intent(in)::n
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ii,ikas,igoflux, icaseb,jx,jx2,b_code
	real::sum_detect,norms,tempxx
	real,dimension(1:numberofpoints2)::weights_temp
	integer::iconsidered,facex,pointx,n_node
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:nof_variables)::leftv,rightv,srf_speedrot
	real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
	real,dimension(1:dimensiona)::pox,poy,poz
	real,dimension(1:nof_variables)::srf_speed
	real,dimension(1:8,1:dimensiona)::vext
	real,dimension(1:8,1:dimensiona)::nodes_list
	real,dimension(1:dimensiona)::cords
	integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf
	kmaxe=xmpielrank(n)
	
	kmaxe=xmpielrank(n)


        do i=1,xmpielrank(n)
            if (ielem_recalc(i).eq.1)then
            rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero 
            end if
        end do



#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	i=el_int(ii)
	iconsidered=i
            mp_source3=zero        
		    b_code=0
		    if (ielem_recalc(i).eq.1)then 
		    
		    do l=1,ielem_ifca(i) !for all their faces
!                                      if (ielem_reorient(l,i).eq.0)then
				  godflux2=zero
				  mp_source2=zero
 				  angle1=ielem_faceanglex(l,i)
 				  angle2=ielem_faceangley(l,i)
 				  nx=(cos(angle1)*sin(angle2))
				  ny=(sin(angle1)*sin(angle2))
				  nz=(cos(angle2))
				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad_n
					weights_temp(1:iqp)=weights_q(1:iqp)
				  else
					iqp=qp_triangle_n
					weights_temp(1:iqp)=weights_t(1:iqp)
				  end if
				  do ngp=1,iqp	!for all the gaussian quadrature points
                    
				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state
				      if ((cascade.eq.2).and.(ielem_mood(i).ge.1))then
				      cleft(1:nof_variables)=u_c_val(3,1:nof_variables,i)
				      end if
				      
				      
				      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))

				       if ((cascade.eq.2).and.(ielem_mood(ielem_ineighn(l,i)).ge.1))then
				      cright(1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
				      end if
				      
				      !right mean flow state
				      
					if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.1)then
                        
					    cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i) !left additional equations flow state
					    cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
					  else
					    cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
					    cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
					  end if
					  
					  
					  cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
					  cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
					end if
			
						  call rotatef(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  end if
						  
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver(n,iconsidered, facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if


				      case(9)			!hll

				      call hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if


				      case(3)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				     				      
				      case(4)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				      
				       
				       
				       
				       case(5)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				      
				       end select
				       
				       
				       
				       
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				     
				     
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  norms=(nx*(u_c_val(1,2,i)/u_c_val(1,1,i)))&
						+(ny*(u_c_val(1,3,i)/u_c_val(1,1,i)))&
						+(nz*(u_c_val(1,4,i)/u_c_val(1,1,i)))
					      if (norms.ge.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
					      end if
					      if (norms.lt.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
					      end if
					  
					  end if
					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
					  
					 
					  
					  
				      end if
				      
				      
				  end do
				    
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				      if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if
				    
!  				     rhs_val(1:nof_variables,ielem_ineigh(l,i))=rhs_val(1:nof_variables,ielem_ineigh(l,i))-godflux2(1:nof_variables)
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
!  				    rhst_val(1:turbulenceequations+passivescalar,ielem_ineigh(l,i))=rhst_val(1:turbulenceequations+passivescalar,ielem_ineigh(l,i))-&
! 				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				    end if
!  				    end if
		    end do
                 if (multispecies.eq.1)then
                 rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
                 
                 end if
            end if
	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)
	iconsidered=i	
		 mp_source3=zero		
		   if (ielem_recalc(i).eq.1)then
		    
		    do l=1,ielem_ifca(i)
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

		    

		    
				      
				  angle1=ielem_faceanglex(l,i)
				  angle2=ielem_faceangley(l,i)
				nx=(cos(angle1)*sin(angle2))
				ny=(sin(angle1)*sin(angle2))
				nz=(cos(angle2))
 				  
 				  if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad_n
					weights_temp(1:iqp)=weights_q(1:iqp)
					n_node=4
				  else
					iqp=qp_triangle_n
					weights_temp(1:iqp)=weights_t(1:iqp)
					n_node=3
				  end if
				  godflux2=zero
				  
				   mp_source2=zero
								  
				  
				  do ngp=1,iqp
				      b_code=0
				      
				      
				      
				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
				     if ((cascade.eq.2).and.(ielem_mood(i).ge.1))then
				     cleft(1:nof_variables)=u_c_val(3,1:nof_variables,i)
				     
				     end if
				      
				      
					 if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
						if (icoupleturb.eq.1)then
						
						
						
							cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i)
						
							
						else
							cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
						end if
					end if
				      
				      
					    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								  
								  cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
								  
								 
								  

								  
								  if ((cascade.eq.2).and.(ielem_mood(ielem_ineigh(l,i)).ge.1))then
								   cright(1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
								  
								  end if
								  
								  
								  
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
									if (cascade.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
									else
									cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									
									end if
									  
									  
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
								    end if
								  
								  
								  
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call  coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

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
								    
								    leftv(1:nof_variables)=cleft(1:nof_variables)
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    
								    
								    
								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    cright(1:nof_variables)=rightv(1:nof_variables)
								    			  				  
								    
								  end if
							else
							
							
                                  
							      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
							      if ((cascade.eq.2).and.(ielem_mood(ielem_ineigh(l,i)).ge.1))then
							      cright(1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
							      
							      end if
							      
							      
								  if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
                                       
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
									  
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
								    end if
							      
							      
							      
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
					     
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                    
! 									  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)

                                    if ((cascade.eq.2).and.(boundhirm(rowfx).ge.1))then


                                    nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								cright(1:nof_variables)=solhir(rowf,1:nof_variables)

!
!                                      cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!  					(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)
!
                                    end if
									   
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
									 
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                                   
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
								    end if
									  
									  

								end if
							else 			
							
                                     
! 								  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1

								  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)

								  
								  if ((cascade.eq.2).and.(boundhirm(rowfx).ge.1))then

								     nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								cright(1:nof_variables)=solhir(rowf,1:nof_variables)

! 								  cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)
								  
								  end if
								  
								  
! 								 
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
									
									
									
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									  
									   
									   
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
								    end if
								  
! 								   
							end if
					    end if
				      
				      
			
						  call rotatef(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
! 						
						  
						  
								    
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  end if
						  
				      
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
								      
								      cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
								      cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
				      
				      
				      
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver(n,iconsidered, facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(3)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      rhllcflux=hllcflux
				      
				      
				      case(9)			!hll

				      call hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				       case(4)			!roe
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				     
				      
				      
				      
				      
				     if (b_code.le.0)then
				      
				      call rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      else
				      
				      
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      end if
				      
				       rhllcflux=hllcflux
				     			
                                    case(5)			!roe
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      
				     if (b_code.le.0)then
				      
				      call troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      else
				      
				      
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      end if
				      
				       rhllcflux=hllcflux
				     
				      
				       end select
				       
				     
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				       
				      
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  norms=(nx*(u_c_val(1,2,i)/u_c_val(1,1,i)))&
						+(ny*(u_c_val(1,3,i)/u_c_val(1,1,i)))&
						+(nz*(u_c_val(1,4,i)/u_c_val(1,1,i)))
					      if (norms.ge.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
					      end if
					      if (norms.lt.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
					      end if
					  
					  end if
					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
					  
					  
					  
				      end if
				      
				      
				  end do
				   
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				    if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if
! 				    if ((igoflux.eq.1))then
! 				    rhs_val(1:nof_variables,ielem_ineigh(l,i))=rhs_val(1:nof_variables,ielem_ineigh(l,i))-godflux2(1:nof_variables)
! 				    end if
				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
! 				     if ((igoflux.eq.1))then
! 				     rhst_val(1:turbulenceequations+passivescalar,ielem_ineigh(l,i))=rhst_val(1:turbulenceequations+passivescalar,ielem_ineigh(l,i))-&
! 				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
! 				    end if
				    end if
! 				    end if
		    end do
		     if (multispecies.eq.1)then
                 rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
                 
                 end if
                end if
	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine calculate_fluxeshi_convective_mood



subroutine calculate_fluxeshi_convective2d_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
	implicit none
	integer,intent(in)::n
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ii,ikas,igoflux, icaseb,kxk,b_code
	real::sum_detect,norms
	real,dimension(1:numberofpoints2)::weights_temp
	integer::iconsidered,facex,pointx,n_node
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:nof_variables)::leftv,rightv,srf_speedrot
	real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
	real,dimension(1:dimensiona)::pox,poy,poz
	real,dimension(1:nof_variables)::srf_speed
	real,dimension(1:8,1:dimensiona)::vext,nodes_list
	real,dimension(1:dimensiona)::cords
	integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf
	kmaxe=xmpielrank(n)
	

	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)


	
	do i=1,kmaxe
	if (ielem_recalc(i).eq.1)then
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero 
	end if
	end do
	

#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	i=el_int(ii)
	iconsidered=i
! 		    rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero 
          if (ielem_recalc(i).eq.1)then             
                mp_source3=zero     
		    do l=1,ielem_ifca(i) !for all their faces

                                
                                
				  godflux2=zero
				  mp_source2=zero
 				  nx=ielem_faceanglex(l,i)
 				  ny=ielem_faceangley(l,i)
 				  angle1=nx
 				  angle2=ny
				 b_code=0
					iqp=qp_line_n
				
				  do ngp=1,iqp	!for all the gaussian quadrature points
				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state
				      
				       if ((cascade.eq.2).and.(ielem_mood(i).ge.1))then
				      cleft(1:nof_variables)=u_c_val(3,1:nof_variables,i)
				      end if
				      
				      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i)) !right mean flow state
				      
				      if ((cascade.eq.2).and.(ielem_mood(ielem_ineigh(l,i)).ge.1))then
				      cright(1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
				      end if
				      
				      
! 				      
					if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.1)then
					    cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i) !left additional equations flow state
					    cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
					  else
					    cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
					    cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
					  end if
					  
					  
					  cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
					  cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
					end if
			
						  call rotatef2d(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef2d(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  call lmacht2d(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  end if
						  
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      case(3)			!roe
				      
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny)
				      
				      rhllcflux(1:nof_variables)=hllcflux(1:nof_variables)
                                        


                      case(9)			!hll

				      call hll_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
                                        
                                         case(4)			!roe
				      
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      call rroe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,b_code)
				      
				      rhllcflux=hllcflux
				     
				      
				       end select
				       
				       
				       
				       
				       
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				       if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  norms=(nx*(u_c_val(1,2,i)/u_c_val(1,1,i)))&
						+(ny*(u_c_val(1,3,i)/u_c_val(1,1,i)))
					      if (norms.ge.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
					      end if
					      if (norms.lt.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
					      end if
					  
					  end if
					  
 					
					  
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
					  
! 					 
					  
				      end if
				      
				      
				  end do
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				    if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if
				    
				    
! 				     rhs_val(1:nof_variables,ielem_ineigh(l,i))=rhs_val(1:nof_variables,ielem_ineigh(l,i))-godflux2(1:nof_variables)
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
! 				    rhst_val(1:turbulenceequations+passivescalar,ielem_ineigh(l,i))=rhst_val(1:turbulenceequations+passivescalar,ielem_ineigh(l,i))-&
!  				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
! 				    
				    end if
! 				    end if
		    end do
		    if (multispecies.eq.1)then
                 rhs_val(7,i)=rhs_val(7,i)-(u_c_val(1,7,i)*mp_source3)!*ielem_totvolume(i))
                 
                 end if
                 
                 
            end if     
	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	
	
#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)
	iconsidered=i	
	mp_source3=zero  
		 if (ielem_recalc(i).eq.1)then	

		    
		    do l=1,ielem_ifca(i)
				  igoflux=0
				  if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                        if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                                icaseb=1        !periodic mine
!                                                                     
                                                        else
                                                                icaseb=3        !physical
!                                                                        
                                                        end if
                                                   
                                                    else
!                                                                  
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

				    b_code=0  
				 nx=ielem_faceanglex(l,i)
 				  ny=ielem_faceangley(l,i)
 				  angle1=nx
 				  angle2=ny
				 
					iqp=qp_line_n
				  godflux2=zero
				  mp_source2=zero
				  do ngp=1,iqp
				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
				      
				      if ((cascade.eq.2).and.(ielem_mood(i).ge.1))then
				     cleft(1:nof_variables)=u_c_val(3,1:nof_variables,i)
				     
				     end if
				      
				      
				      
					 if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
						if (icoupleturb.eq.1)then
							cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i)
						else
							cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
						end if
					end if
				      
				      
					    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								  cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
								  
								   if ((cascade.eq.2).and.(ielem_mood(ielem_ineigh(l,i)).ge.1))then
								   cright(1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
								  
								  end if
								  
								  
								  
								  
								  
								  
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
								    end if
								  
								  
								  ikas=1
								  
								  else
								  !not periodic ones in my cpu
								   
								  facex=l;iconsidered=i
								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
								    cords(1:2)=zero
								    n_node=2
								    cords(1:2)=cordinates2(n,nodes_list,n_node)
							    
								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    
								    
								    leftv(1:nof_variables)=cleft(1:nof_variables)
! 								    
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    cright(1:nof_variables)=rightv(1:nof_variables)
! 				  				   
				  				  	 ikas=2			  				  
								    
								  end if
							else
							      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
 							     
 							      if ((cascade.eq.2).and.(ielem_mood(ielem_ineigh(l,i)).ge.1))then
							      cright(1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
							      
							      end if
 							     
 							     
 							     
 							     
 							     
								  if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
								    end if
							      
							      
							       ikas=3
							      
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
					    
					    
					     
					    
					    
						
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
! 									  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1

									  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
									  
									  if ((cascade.eq.2).and.(boundhirm(rowfx).ge.1))then
!                                     cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

									   nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								cright(1:nof_variables)=solhir(rowf,1:nof_variables)

                                    
                                    end if
									  
									  
									  
									  
 									   
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
									
									
									
									
								    end if
									  
									  

								end if
							else 			
							
! 								  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
! 								  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
								  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
 								  
 								  
 								   if ((cascade.eq.2).and.(boundhirm(rowfx).ge.1))then
!                                     cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

									   nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								cright(1:nof_variables)=solhir(rowf,1:nof_variables)

                                    
                                    end if
 								  
 								  
 								  
 								  
 								  
 								  
 								  
 								  
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol&
! 									   (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
									
									
								    end if
								   ikas=4
! 								   
							end if
					    end if
				      
				      
			
						  call rotatef2d(n,cright_rot,cright,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
						  call rotatef2d(n,cleft_rot,cleft,angle1,angle2)	!rotate wrt to normalvector of face and solve 1d riemann problem
! 						   
						  
						  
						  if ((lmach.eq.1))then    !application of the low mach number correction
						  leftv(1:nof_variables)=cleft_rot(1:nof_variables); rightv(1:nof_variables)=cright_rot(1:nof_variables)
						  
						  call lmacht2d(n,leftv,rightv)
						  cleft_rot(1:nof_variables)=leftv(1:nof_variables);cright_rot(1:nof_variables)=rightv(1:nof_variables);
						  
						  
						  end if
						  
				      
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
								      
								      cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
								      cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
				      
				      
				      
				      
				      select case(iriemann)
				      
				      case(1)			!hllc
				      
				      call hllc_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if
				      
				      case(9)			!hll

				      call hll_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
				      end if



				      case(3)			!roe
				      
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny)
				      
				      rhllcflux=hllcflux
				     				      
				     
				        case(4)			!roe
				      
				      
				       
				      
				      if ((b_code.le.0))then
				      
				      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb2d(n,cright,cright_rot,angle1,angle2)
				      
				      
				      
				      call rroe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,b_code)
				      rhllcflux=hllcflux
				      else
				      
!
				     
				      
				      call rusanov_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot)
				      call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
				      
				      
				      
				      end if
				      
				      
				       end select
				       
				     
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				       
				      
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  norms=(nx*(u_c_val(1,2,i)/u_c_val(1,1,i)))&
						+(ny*(u_c_val(1,3,i)/u_c_val(1,1,i)))
						
					      if (norms.ge.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
					      end if
					      if (norms.lt.zero)then
						rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
					      end if
					  
					  end if
					  
					  if ((b_code.eq.4).or.(b_code.eq.3))then
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
					  end if
					  
					  
					 
					  godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
					  (rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)*(weights_temp(ngp)*ielem_surf(l,i)))
! 					  
					  
					  
				      end if
				      
				      
				  end do
				   
				    rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
				    if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    rhst_val(1:turbulenceequations+passivescalar,i)=rhst_val(1:turbulenceequations+passivescalar,i)+&
				    godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

				    end if
		    end do
		     if (multispecies.eq.1)then
                 rhs_val(7,i)=rhs_val(7,i)-(u_c_val(1,7,i)*mp_source3)!*ielem_totvolume(i))
                 
                 end if
                 
                 end if
	end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine calculate_fluxeshi_convective2d_mood


! !experimental subroutine---not for production
! subroutine vertex_neighbours_values(n)
! implicit none
! integer,intent(in)::n
! integer::i,j,k,l,m,nj
! real,dimension(4,30,nof_variables)::vertex_neighbours_vals	!maximum allocation not efficient
!
! i=iconsidered
!
! write(680+n,*),"element number",ielem_ihexgl(i)
!
! do nj=1,ielem_nonodes(i)
!
! 	if (rec_local(i).eq.1)then
! 		do j=1,ielem_nojecount(nj,i)
! 			vertex_neighbours_vals(nj,j,1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,ielem_nodes_neighbours(nj,j,i,i)))
! 		!rec_ihexl(1,l,i)	!local numbering in my cpu
! 		!rec_ihexg(1,l,i)	!global numbering in my cpu
! 			write(680+n,*),"node",nj,"node neighbour",j,"local number",ielem_nodes_neighbours(nj,j,i), "values",vertex_neighbours_vals(nj,j,1:nof_variables)
!
! 		end do
!
! 	else
! 			do j=1,ielem_nojecount(nj,i)
!
! 				if (rec_ihexb(1,ielem_nodes_neighbours(nj,j,i,i)).eq.n)then
! 					vertex_neighbours_vals(nj,j,1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,ielem_nodes_neighbours(nj,j,i,i)))
! 				else
! 					write(500+n,*)rec_ihexn(1,ielem_nodes_neighbours(nj,j,i,i)),ielem_nodes_neighbours(nj,j,i)
! 					vertex_neighbours_vals(nj,j,1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_nodes_neighbours(nj,j,i,i)))%sol(rec_ihexl(1,ielem_nodes_neighbours(nj,j,i,i)),1:nof_variables)
!
! 				end if
! 			write(680+n,*),"node",nj,"node neighbour",j,"local number",ielem_nodes_neighbours(nj,j,i), "values",vertex_neighbours_vals(nj,j,1:nof_variables)
! 			end do
! 	end if
!
! end do
!
! !now you have collected for every vertex the values
!
!
! end subroutine vertex_neighbours_values






end module fluxes
