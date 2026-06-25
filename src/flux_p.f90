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
 
 call dg_vol_integral2_fill(n,i,solution_integ2)
 
 
 end subroutine solution_integ

 subroutine solution_integ_s(i,solution_integ_strong)
 implicit none
#ifdef gpu
!$omp declare target
#endif
  real,dimension(1:nof_variables),intent(inout)::solution_integ_strong
 integer,intent(in)::i

 call dg_vol_integral_strong_fill(n,i,solution_integ_strong)


 end subroutine solution_integ_s


  subroutine solution_integ_w(i,solution_integ_weak)
 implicit none
#ifdef gpu
!$omp declare target
#endif
 real,dimension(1:nof_variables),intent(inout)::solution_integ_weak
 integer,intent(in)::i

 call dg_vol_integral_weak_fill(n,i,solution_integ_weak)


 end subroutine solution_integ_w




subroutine calculate_fluxeshi(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation
		implicit none
		integer,intent(in)::n
		integer::i,kmaxe

		kmaxe=xmpielrank(n)
	



#ifdef xpu
if (dg.eq.1)then
!$omp barrier
!$omp do
	do i=1,kmaxe
	call calculate_fluxeshi_cell(n,i)
	end do
!$omp end do
else
	call calculate_fluxeshi_cell_xpu_fv(n)
end if
#else
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(i)
#else
!$omp barrier
!$omp do
#endif
	do i=1,kmaxe
	call calculate_fluxeshi_cell(n,i)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
#endif





end subroutine calculate_fluxeshi

subroutine calculate_fluxeshi_cell_xpu_fv(n)
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe
	real::godflux2
	integer::l,ngp,iqp,nfx,lfx,rowfx
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real::angle1,angle2,normalvect
	real::cleft1,cright1,hllcflux1

	kmaxe=xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n,kmaxe) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, zero, lamx, lamy, lamz, qp_quad) &
!$omp& firstprivate(qp_triangle) &
!$omp& map(alloc: weights_q, weights_t, ielem_interior, ielem_ifca) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_types_faces, ielem_surf) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode) &
!$omp& map(alloc: ielem_qface, ielem_inter_id, ielem_indexf, bound_offset, boundhir, rec_uleft, rhs_val) &
!$omp& private(i,l,ngp,iqp,nfx,lfx,rowfx,weights_temp,angle1,angle2,normalvect,cleft1,cright1,hllcflux1,godflux2)
#endif
	do i=1,kmaxe
	rhs_val(1,i)=zero

	do l=1,ielem_ifca(i)
	  godflux2=zero
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

	  do ngp=1,iqp
	    cleft1=rec_uleft(1,l,ngp,i)

	    if (ielem_interior(i).eq.0)then
	      cright1=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
	    else
	      if (ielem_ineighb(l,i).eq.n)then
	        if (ielem_ibounds(l,i).gt.0)then
	          if (ibound_icode(ielem_ibounds(l,i)).eq.5)then
	            cright1=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
	          else
	            cright1=cleft1
	          end if
	        else
	          cright1=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
	        end if
	      else
	        if (ielem_ibounds(l,i).gt.0)then
	          if (ibound_icode(ielem_ibounds(l,i)).eq.5)then
	            nfx=ielem_ineighn(l,i)
	            lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
	            rowfx=bound_offset(nfx)+lfx-1
	            cright1=boundhir(rowfx,1)
	          else
	            cright1=cleft1
	          end if
	        else
	          nfx=ielem_ineighn(l,i)
	          lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
	          rowfx=bound_offset(nfx)+lfx-1
	          cright1=boundhir(rowfx,1)
	        end if
	      end if
	    end if

	    hllcflux1=zero
	    if (normalvect.gt.zero)then
	      hllcflux1=normalvect*cleft1
	    else
	      hllcflux1=normalvect*cright1
	    end if
	    godflux2=godflux2+(hllcflux1*(weights_temp(ngp)*ielem_surf(l,i)))
	  end do
	  rhs_val(1,i)=rhs_val(1,i)+godflux2
	end do
	end do

#ifdef xpu
!$omp end target teams distribute parallel do
#endif

end subroutine calculate_fluxeshi_cell_xpu_fv

subroutine calculate_fluxeshi_cell(n,i)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::i
	integer,intent(in)::n
	integer::kmaxe
		real::godflux2
		integer::l,ngp,iqp
		real,dimension(1:gpu_max_qp_face)::weights_temp
		real,dimension(1:8,1:gpu_max_dim)::vext
		real::angle1,angle2,nx,ny,nz,normalvect
		integer::facex,pointx,iconsidered,nfx,lfx,rowfx
		real,dimension(1:gpu_max_nvar_total)::cleft,cright,hllcflux,rhllcflux




        iconsidered=i


	        if (dg.eq.1)then
	        rhs_valdg(:,:,i) = 0.0d0
	        end if
		
	        if (dg.eq.1) then
	            call dg_vol_integral_accumulate_rhs(n,iconsidered,-1.0d0)
	            
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
				      
					      call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
				      
				       
				      
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
					      call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
				      
				       
				     
				      else
				      
				      godflux2=godflux2+(hllcflux(1)*(weights_temp(ngp)*ielem_surf(l,i)))	
				      
				      end if
				  end do
				    rhs_val(1,i)=rhs_val(1,i)+godflux2

				    
		    end do
		end if
! 				
	end subroutine calculate_fluxeshi_cell

	
	
	
	
	
	
subroutine calculate_fluxeshi2d(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation in 2d
		implicit none
		integer,intent(in)::n
		integer::i,kmaxe


		kmaxe = xmpielrank(n)


#ifdef xpu
if (dg.eq.1)then
!$omp barrier
!$omp do
	do i=1,kmaxe
	call calculate_fluxeshi2d_cell(n,i)
	end do
!$omp end do
else
	call calculate_fluxeshi2d_cell_xpu_fv(n)
end if
#else
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(i)
#else
!$omp barrier
!$omp do
#endif
	do i=1,kmaxe
	call calculate_fluxeshi2d_cell(n,i)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
#endif

end subroutine calculate_fluxeshi2d

subroutine calculate_fluxeshi2d_cell_xpu_fv(n)
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe
	real::godflux2,lamxl,lamyl
	integer::l,ngp,nfx,lfx,rowfx
	real::nx,ny,normalvect
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real::cleft1,cright1,hllcflux1

	kmaxe=xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n,kmaxe) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, zero, lamx, lamy, initcond, qp_line_n) &
!$omp& map(alloc: weights_l, ielem_interior, ielem_ifca, ielem_xxc, ielem_yyc, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_surf, ielem_ineighn, ielem_ineigh, ielem_ineighb) &
!$omp& map(alloc: ielem_ibounds, ibound_icode, ielem_qface, ielem_inter_id, ielem_indexf) &
!$omp& map(alloc: bound_offset, boundhir, rec_uleft, rhs_val) &
!$omp& private(i,l,ngp,nfx,lfx,rowfx,nx,ny,normalvect,weights_temp,cleft1,cright1,hllcflux1,godflux2,lamxl,lamyl)
#endif
	do i=1,kmaxe
	weights_temp(1:qp_line_n)=weights_l(1:qp_line_n)

	if (initcond.eq.3)then
	  lamxl=-ielem_yyc(i)+0.5d0
	  lamyl=ielem_xxc(i)-0.5d0
	else
	  lamxl=lamx
	  lamyl=lamy
	end if

	rhs_val(1,i)=zero

	do l=1,ielem_ifca(i)
	  godflux2=zero
	  nx=ielem_faceanglex(l,i)
	  ny=ielem_faceangley(l,i)
	  normalvect=(nx*lamxl)+(ny*lamyl)

	  do ngp=1,qp_line_n
	    cleft1=rec_uleft(1,l,ngp,i)

	    if (ielem_interior(i).eq.0)then
	      cright1=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
	    else
	      if (ielem_ineighb(l,i).eq.n)then
	        if (ielem_ibounds(l,i).gt.0)then
	          if (ibound_icode(ielem_ibounds(l,i)).eq.5)then
	            cright1=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
	          else
	            cright1=cleft1
	          end if
	        else
	          cright1=rec_uleft(1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
	        end if
	      else
	        if (ielem_ibounds(l,i).gt.0)then
	          if (ibound_icode(ielem_ibounds(l,i)).eq.5)then
	            nfx=ielem_ineighn(l,i)
	            lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
	            rowfx=bound_offset(nfx)+lfx-1
	            cright1=boundhir(rowfx,1)
	          else
	            cright1=cleft1
	          end if
	        else
	          nfx=ielem_ineighn(l,i)
	          lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
	          rowfx=bound_offset(nfx)+lfx-1
	          cright1=boundhir(rowfx,1)
	        end if
	      end if
	    end if

	    hllcflux1=zero
	    if (normalvect.gt.zero)then
	      hllcflux1=normalvect*cleft1
	    else
	      hllcflux1=normalvect*cright1
	    end if
	    godflux2=godflux2+(hllcflux1*(weights_temp(ngp)*ielem_surf(l,i)))
	  end do
	  rhs_val(1,i)=rhs_val(1,i)+godflux2
	end do
	end do

#ifdef xpu
!$omp end target teams distribute parallel do
#endif

end subroutine calculate_fluxeshi2d_cell_xpu_fv

subroutine calculate_fluxeshi2d_cell(n,i)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::i
	integer,intent(in)::n
	integer::kmaxe
	real::godflux2,sum_detect,lamxl,lamyl
	integer::l,k,ngp,iqp,neighbor_index,neighbor_face_index
	real,dimension(1:8,1:gpu_max_dim)::vext
	 real::angle1,angle2,nx,ny,nz,normalvect
    integer::facex,pointx,iconsidered,nfx,lfx,rowfx
    real,dimension(1:gpu_max_qp_face)::weights_temp
    real,dimension(1:gpu_max_nvar_total)::cleft,cright,hllcflux,rhllcflux





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
        else
        rhs_val(:,i)=zero
        end if
        
        
        iconsidered=i
        
        if (dg.eq.1) then
        
            call dg_vol_integral_accumulate_rhs(n,iconsidered,-1.0d0)
         
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
                         
                        call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)



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
                         call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
                         




                    else !fv
                        godflux2=godflux2+(hllcflux(1)*(weights_temp(ngp)*ielem_surf(l,i)))
                    end if
  
                end do
                
                if (dg /= 1) rhs_val(1,i)=rhs_val(1,i)+godflux2
		    end do



		end if
		
		
	end subroutine calculate_fluxeshi2d_cell


subroutine calculate_fluxeshi_convective(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
		implicit none
		integer,intent(in)::n
		integer::ii


	if ((dg.eq.0).and.(multispecies.eq.0).and.(realgas.eq.0).and. &
     & ((iriemann.eq.1).or.(iriemann.eq.2).or.(iriemann.eq.9)))then
#ifdef xpu
!$omp barrier
!$omp master
#endif
	call calculate_fluxeshi_convective_inner_cell_ideal(n)
	call calculate_fluxeshi_convective_bound_cell_fv_ideal(n)
#ifdef xpu
!$omp end master
!$omp barrier
#endif
	return
	end if


#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	call calculate_fluxeshi_convective_inner_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	


	if (dg.eq.1)then
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_convective_bound_cell_dg_volume(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_convective_bound_cell_dg_surface(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	else
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_convective_bound_cell_fv(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if



end subroutine calculate_fluxeshi_convective

subroutine calculate_fluxeshi_convective_inner_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,ii
	integer::i,l,ngp,iqp,ikas,igoflux,icaseb,jx,jx2,b_code,srf,iv,nvt
	real,dimension(1:gpu_max_nvar_total)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	real::sum_detect,norms,tempxx,wface
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx,nfx,lfx,rowfx
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext
	real::x1,y1,z1

	i=el_int(ii)

	nvt=turbulenceequations+passivescalar
	do iv=1,nof_variables
	rhs_val(iv,i)=zero
	end do
	if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	do iv=1,nvt
	rhst_val(iv,i)=zero
	end do
	end if




	iconsidered=i
            mp_source3=zero        

        if(mrf.eq.1)then
            srf=rec_mrf(i)
        end if
			     if (dg.eq.1) then
			     rhs_valdg(:,:,i) = zero
	            
	             call dg_vol_integral_accumulate_rhs(n,iconsidered,-1.0d0)

	            end if
            
		    
		    do l=1,ielem_ifca(i) !for all their faces
				b_code=0

				  do iv=1,nof_variables+nvt
				  godflux2(iv)=zero
				  end do
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
				      
				      do iv=1,nvt
				      rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
				      end do
				      end if

				      case(9)			!hll

				      call hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				      call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      do iv=1,nvt
				      rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
				      end do
				      end if
				      
				      case(2)			!rusanov
				      
				      call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
				       call rotateb(n,rhllcflux,hllcflux,angle1,angle2)






				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      do iv=1,nvt
				      rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
				      end do
				      end if
				      case(3)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      do iv=1,nof_variables+nvt
				      rhllcflux(iv)=hllcflux(iv)
				      end do
				     				      
				      case(4)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      do iv=1,nof_variables+nvt
				      rhllcflux(iv)=hllcflux(iv)
				      end do
				      
				       
				       
				       
				       case(5)			!roe
				      
				      
				      call rotateb(n,cleft,cleft_rot,angle1,angle2)
				      call rotateb(n,cright,cright_rot,angle1,angle2)
				      call troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
				      
				      do iv=1,nof_variables+nvt
				      rhllcflux(iv)=hllcflux(iv)
				      end do
				      
				       end select
				       
				      wface=weights_temp(ngp)*ielem_surf(l,i)
				      if (dg.eq.1)then
				      
					      call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
				      
				      
				      
				      else 
				       
				       
				      
				      do iv=1,nof_variables
				      godflux2(iv)=godflux2(iv)+(rhllcflux(iv)*wface)
				      end do
				      
				      end if
				      
				      
				      
				      
				      
				      
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*wface
                        end if
				     
				     
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					  
					  
					  
					  norms=0.5*(cleft_rot(2)+cright_rot(2))


					if (rec_mrf(iconsidered).eq.1)then
					  norms=norms-srf_speedrot(2)
					end if
					  do iv=1,nvt
					  rhllcflux(nof_variables+iv)=0.5d0*((norms*(cturbl(iv)+cturbr(iv)))+(abs(norms)*(cturbl(iv)-cturbr(iv))))
					  end do
					  
					  

					  end if
					  
					  do iv=1,nvt
					  godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+(rhllcflux(nof_variables+iv)*wface)
					  end do
					  
					 
					  
					  
				      end if
				      
				      
				  end do
				    
				    do iv=1,nof_variables
				    rhs_val(iv,i)=rhs_val(iv,i)+godflux2(iv)
				    end do
				      if (multispecies.eq.1)then
                        mp_source3=mp_source3+mp_source2
                        end if
				    

				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    do iv=1,nvt
				    rhst_val(iv,i)=rhst_val(iv,i)+godflux2(nof_variables+iv)
				    end do

				    end if

		    end do
                 if (multispecies.eq.1)then
                 rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
                 
                 end if
                 
	                 if (dg.eq.1)then
	                    if (multispecies.eq.1)then


						rhs_valdg(1,8,i)=rhs_valdg(1,8,i)-(u_c_val(1,8,i)*mp_source3)

					end if

                   




                 end if

end subroutine calculate_fluxeshi_convective_inner_cell


subroutine calculate_fluxeshi_convective_bound_point(n,iconsidered,facex,pointx,rhllcflux,mp_source1,b_code)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,iconsidered,facex,pointx
	integer,intent(out)::b_code
	integer::iv,nvt
	real,intent(out)::mp_source1
	real,dimension(1:gpu_max_nvar_total),intent(out)::rhllcflux
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot,hllcflux
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real::angle1,angle2,nx,ny,nz,norms

	nvt=turbulenceequations+passivescalar
	b_code=0
	mp_source1=zero
	rhllcflux=zero
	hllcflux=zero
	srf_speedrot=zero

	angle1=ielem_faceanglex(facex,iconsidered)
	angle2=ielem_faceangley(facex,iconsidered)
	nx=(cos(angle1)*sin(angle2))
	ny=(sin(angle1)*sin(angle2))
	nz=(cos(angle2))

	call get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
	     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)

	call rotatef(n,cright_rot,cright,angle1,angle2)
	call rotatef(n,cleft_rot,cleft,angle1,angle2)

	if ((lmach.eq.1))then
	  leftv(1:nof_variables)=cleft_rot(1:nof_variables)
	  rightv(1:nof_variables)=cright_rot(1:nof_variables)
	  call lmacht(n,leftv,rightv)
	  cleft_rot(1:nof_variables)=leftv(1:nof_variables)
	  cright_rot(1:nof_variables)=rightv(1:nof_variables)
	end if

	if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	  do iv=1,nvt
	    cleft_rot(nof_variables+iv)=cturbl(iv)
	    cright_rot(nof_variables+iv)=cturbr(iv)
	  end do
	end if

	select case(iriemann)
	case(1)
	  call hllc_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	  call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
	  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	    do iv=1,nvt
	      rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
	    end do
	  end if
	case(9)
	  call hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	  call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
	  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	    do iv=1,nvt
	      rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
	    end do
	  end if
	case(2)
	  call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	  call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
	  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	    do iv=1,nvt
	      rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
	    end do
	  end if
	case(3)
	  call rotateb(n,cleft,cleft_rot,angle1,angle2)
	  call rotateb(n,cright,cright_rot,angle1,angle2)
	  call roe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	  do iv=1,nof_variables+nvt
	    rhllcflux(iv)=hllcflux(iv)
	  end do
	case(4)
	  if (b_code.gt.0)then
	    call rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	    call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
	    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      do iv=1,nvt
	        rhllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)
	      end do
	    end if
	  else
	    call rotateb(n,cleft,cleft_rot,angle1,angle2)
	    call rotateb(n,cright,cright_rot,angle1,angle2)
	    call rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	    do iv=1,nof_variables+nvt
	      rhllcflux(iv)=hllcflux(iv)
	    end do
	  end if
	case(5)
	  call rotateb(n,cleft,cleft_rot,angle1,angle2)
	  call rotateb(n,cright,cright_rot,angle1,angle2)
	  if (b_code.le.0)then
	    call troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	  else
	    call roe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	  end if
	  do iv=1,nof_variables+nvt
	    rhllcflux(iv)=hllcflux(iv)
	  end do
	end select

	if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	  if (icoupleturb.eq.0)then
	    norms=0.5d0*(cleft_rot(2)+cright_rot(2))
	    if (rec_mrf(iconsidered).eq.1)then
	      norms=norms-srf_speedrot(2)
	    end if
	    do iv=1,nvt
	      rhllcflux(nof_variables+iv)=0.5d0*((norms*(cturbl(iv)+cturbr(iv)))+(abs(norms)*(cturbl(iv)-cturbr(iv))))
	    end do
	  end if
	  if ((b_code.eq.3).or.(b_code.eq.4))then
	    do iv=1,nvt
	      rhllcflux(nof_variables+iv)=zero
	    end do
	  end if
	end if
end subroutine calculate_fluxeshi_convective_bound_point


subroutine calculate_fluxeshi_convective_bound_cell_fv(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,ii
	integer::i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code
	real,dimension(1:gpu_max_nvar_total)::godflux2,rhllcflux
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real::mp_source1,mp_source2,mp_source3,wface

	i=el_bnd(ii)
	nvt=turbulenceequations+passivescalar
	do iv=1,nof_variables
	  rhs_val(iv,i)=zero
	end do
	if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	  do iv=1,nvt
	    rhst_val(iv,i)=zero
	  end do
	end if

	iconsidered=i
	mp_source3=zero

	do l=1,ielem_ifca(i)
	  facex=l
	  if (ielem_types_faces(l,i).eq.5)then
	    iqp=qp_quad_n
	    weights_temp(1:iqp)=weights_q(1:iqp)
	  else
	    iqp=qp_triangle_n
	    weights_temp(1:iqp)=weights_t(1:iqp)
	  end if

	  do iv=1,nof_variables+nvt
	    godflux2(iv)=zero
	  end do
	  mp_source2=zero

	  do ngp=1,iqp
	    pointx=ngp
	    call calculate_fluxeshi_convective_bound_point(n,iconsidered,facex,pointx,rhllcflux,mp_source1,b_code)
	    wface=weights_temp(ngp)*ielem_surf(l,i)
	    do iv=1,nof_variables
	      godflux2(iv)=godflux2(iv)+(rhllcflux(iv)*wface)
	    end do
	    if (multispecies.eq.1)then
	      mp_source2=mp_source2+mp_source1*wface
	    end if
	    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      do iv=1,nvt
	        godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+(rhllcflux(nof_variables+iv)*wface)
	      end do
	    end if
	  end do

	  do iv=1,nof_variables
	    rhs_val(iv,i)=rhs_val(iv,i)+godflux2(iv)
	  end do
	  if (multispecies.eq.1)then
	    mp_source3=mp_source3+mp_source2
	  end if
	  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	    do iv=1,nvt
	      rhst_val(iv,i)=rhst_val(iv,i)+godflux2(nof_variables+iv)
	    end do
	  end if
	end do

	if (multispecies.eq.1)then
	  rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
	end if
end subroutine calculate_fluxeshi_convective_bound_cell_fv


subroutine calculate_fluxeshi_convective_inner_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz
real::angle1,angle2,nx,ny,nz,wface,norms

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, turbulenceequations, passivescalar, turbulence, icoupleturb) &
!$omp& firstprivate(lmach, lmach_style, iriemann, mrf, per_rot, angle_per, zero, oo2, gamma, adda, qp_quad_n) &
!$omp& firstprivate(qp_triangle_n) &
!$omp& map(alloc: weights_q, weights_t, el_int, ielem_ifca, ielem_types_faces) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_ineighn, ielem_ineigh) &
!$omp& map(alloc: rec_mrf, rec_rotvel, rec_uleft, rec_uleftturb) &
!$omp& map(alloc: u_ct_val, rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#endif
do ii=1,nof_interior
i=el_int(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
rhs_val(1:nof_variables,i)=zero
if (nvt.gt.0) rhst_val(1:nvt,i)=zero
do l=1,ielem_ifca(i)
  facex=l
  b_code=0
  godflux2=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  nx=cos(angle1)*sin(angle2)
  ny=sin(angle1)*sin(angle2)
  nz=cos(angle2)
  if (ielem_types_faces(l,i).eq.5)then
    iqp=qp_quad_n
  else
    iqp=qp_triangle_n
  end if
  do ngp=1,iqp
    pointx=ngp
    hllcflux=zero
    rhllcflux=zero
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    cturbl=zero
    cturbr=zero
    call get_states_interior_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    call rotatef(n,cright_rot,cright,angle1,angle2)
    call rotatef(n,cleft_rot,cleft,angle1,angle2)
    if (lmach.eq.1)then
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
    end if
    call riemann_flux_ideal3d(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,srf_speedrot)
    call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
    if (nvt.gt.0) rhllcflux(nof_variables+1:nof_variables+nvt)=hllcflux(nof_variables+1:nof_variables+nvt)
    if (nvt.gt.0)then
      if (icoupleturb.eq.0)then
        norms=0.5d0*(cleft_rot(2)+cright_rot(2))
        if (rec_mrf(iconsidered).eq.1) norms=norms-srf_speedrot(2)
        do iv=1,nvt
          rhllcflux(nof_variables+iv)=0.5d0*((norms*(cturbl(iv)+cturbr(iv)))+(abs(norms)*(cturbl(iv)-cturbr(iv))))
        end do
      end if
    end if
    if (ielem_types_faces(l,i).eq.5)then
      wface=weights_q(ngp)*ielem_surf(l,i)
    else
      wface=weights_t(ngp)*ielem_surf(l,i)
    end if
    do iv=1,nof_variables+nvt
      godflux2(iv)=godflux2(iv)+rhllcflux(iv)*wface
    end do
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)+godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_convective_inner_cell_ideal


subroutine calculate_fluxeshi_convective_bound_cell_fv_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz
real::angle1,angle2,nx,ny,nz,wface,norms

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, turbulenceequations, passivescalar, turbulence, icoupleturb) &
!$omp& firstprivate(lmach, lmach_style, iriemann, mrf, per_rot, angle_per, zero, oo2, gamma, adda, qp_quad_n) &
!$omp& firstprivate(r_gas, pres, rres, boundtype, initcond, Mach_in, uvel, vvel, wvel, press_outlet, swirl, tolsmall, qp_triangle_n) &
!$omp& map(alloc: weights_q, weights_t, el_bnd, ielem_ifca, ielem_types_faces) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_facediss, ielem_ineighn, ielem_ineigh) &
!$omp& map(alloc: ielem_ineighb, ielem_ibounds, ibound_icode, ielem_qface, ielem_inter_id, ielem_indexf) &
!$omp& map(alloc: bound_offset) &
!$omp& map(alloc: boundhir, ielem_nodes_faces, inoder4_cord, rec_mrf, rec_rotvel, rec_uleft, rec_uleftturb) &
!$omp& map(alloc: u_ct_val, rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#endif
do ii=1,nof_bounded
i=el_bnd(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
rhs_val(1:nof_variables,i)=zero
if (nvt.gt.0) rhst_val(1:nvt,i)=zero
do l=1,ielem_ifca(i)
  facex=l
  b_code=0
  godflux2=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  nx=cos(angle1)*sin(angle2)
  ny=sin(angle1)*sin(angle2)
  nz=cos(angle2)
  if (ielem_types_faces(l,i).eq.5)then
    iqp=qp_quad_n
  else
    iqp=qp_triangle_n
  end if
  do ngp=1,iqp
    pointx=ngp
    hllcflux=zero
    rhllcflux=zero
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    cturbl=zero
    cturbr=zero
    call get_states_bounds_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    call rotatef(n,cright_rot,cright,angle1,angle2)
    call rotatef(n,cleft_rot,cleft,angle1,angle2)
    if (lmach.eq.1)then
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
    end if
    call riemann_flux_ideal3d(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,srf_speedrot)
    call rotateb(n,rhllcflux,hllcflux,angle1,angle2)
    if (nvt.gt.0) rhllcflux(nof_variables+1:nof_variables+nvt)=hllcflux(nof_variables+1:nof_variables+nvt)
    if (nvt.gt.0)then
      if (icoupleturb.eq.0)then
        norms=0.5d0*(cleft_rot(2)+cright_rot(2))
        if (rec_mrf(iconsidered).eq.1) norms=norms-srf_speedrot(2)
        do iv=1,nvt
          rhllcflux(nof_variables+iv)=0.5d0*((norms*(cturbl(iv)+cturbr(iv)))+(abs(norms)*(cturbl(iv)-cturbr(iv))))
        end do
      end if
      if ((b_code.eq.3).or.(b_code.eq.4)) rhllcflux(nof_variables+1:nof_variables+nvt)=zero
    end if
    if (ielem_types_faces(l,i).eq.5)then
      wface=weights_q(ngp)*ielem_surf(l,i)
    else
      wface=weights_t(ngp)*ielem_surf(l,i)
    end if
    do iv=1,nof_variables+nvt
      godflux2(iv)=godflux2(iv)+rhllcflux(iv)*wface
    end do
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)+godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_convective_bound_cell_fv_ideal


subroutine calculate_fluxeshi_convective_bound_cell_dg_volume(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,ii
	integer::i,iv,kbasis,nvt,iconsidered,number_of_dog

	i=el_bnd(ii)
	nvt=turbulenceequations+passivescalar
	do iv=1,nof_variables
	  rhs_val(iv,i)=zero
	end do
	if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	  do iv=1,nvt
	    rhst_val(iv,i)=zero
	  end do
	end if

	iconsidered=i
	number_of_dog=ielem_idegfree(iconsidered)
	do iv=1,nof_variables
	  do kbasis=1,number_of_dog+1
	    rhs_valdg(kbasis,iv,i)=zero
	  end do
	end do
	call dg_vol_integral_accumulate_rhs(n,iconsidered,-1.0d0)
end subroutine calculate_fluxeshi_convective_bound_cell_dg_volume


subroutine calculate_fluxeshi_convective_bound_cell_dg_surface(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,ii
	integer::i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code
	real,dimension(1:gpu_max_nvar_total)::rhllcflux
	real,dimension(1:gpu_max_extra_transport)::turbrhs
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real::mp_source1,mp_source2,mp_source3,wface

	i=el_bnd(ii)
	nvt=turbulenceequations+passivescalar
	iconsidered=i
	mp_source3=zero

	do l=1,ielem_ifca(i)
	  facex=l
	  if (ielem_types_faces(l,i).eq.5)then
	    iqp=qp_quad_n
	    weights_temp(1:iqp)=weights_q(1:iqp)
	  else
	    iqp=qp_triangle_n
	    weights_temp(1:iqp)=weights_t(1:iqp)
	  end if

	  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	    do iv=1,nvt
	      turbrhs(iv)=zero
	    end do
	  end if
	  mp_source2=zero

	  do ngp=1,iqp
	    pointx=ngp
	    call calculate_fluxeshi_convective_bound_point(n,iconsidered,facex,pointx,rhllcflux,mp_source1,b_code)
	    wface=weights_temp(ngp)*ielem_surf(l,i)
	    call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
	    if (multispecies.eq.1)then
	      mp_source2=mp_source2+mp_source1*wface
	    end if
	    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      do iv=1,nvt
	        turbrhs(iv)=turbrhs(iv)+(rhllcflux(nof_variables+iv)*wface)
	      end do
	    end if
	  end do

	  if (multispecies.eq.1)then
	    mp_source3=mp_source3+mp_source2
	  end if
	  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	    do iv=1,nvt
	      rhst_val(iv,i)=rhst_val(iv,i)+turbrhs(iv)
	    end do
	  end if
	end do

	if (multispecies.eq.1)then
	  rhs_val(8,i)=rhs_val(8,i)-(u_c_val(1,8,i)*mp_source3)
	end if

	if (multispecies.eq.1)then
	  rhs_valdg(1,8,i)=rhs_valdg(1,8,i)-(u_c_val(1,8,i)*mp_source3)
	end if
end subroutine calculate_fluxeshi_convective_bound_cell_dg_surface



	

	
	
subroutine calculate_fluxeshi_convective2d(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
		implicit none
		integer,intent(in)::n
		integer::ii
	

	if ((dg.eq.0).and.(multispecies.eq.0).and.(realgas.eq.0).and. &
     & ((iriemann.eq.1).or.(iriemann.eq.2).or.(iriemann.eq.9)))then
#ifdef xpu
!$omp barrier
!$omp master
#endif
	call calculate_fluxeshi_convective2d_inner_cell_ideal(n)
	call calculate_fluxeshi_convective2d_bound_cell_ideal(n)
#ifdef xpu
!$omp end master
!$omp barrier
#endif
	return
	end if


#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	call calculate_fluxeshi_convective2d_inner_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	

	
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_convective2d_bound_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine calculate_fluxeshi_convective2d

subroutine calculate_fluxeshi_convective2d_inner_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	integer::i,kmaxe
	real,dimension(1:gpu_max_nvar_total)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::l,ngp,iqp,ikas,igoflux,icaseb,kxk,b_code
	real::sum_detect,norms
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real::x1,y1,z1
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext




	i=el_int(ii)


	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero

	iconsidered=i
	            if (dg.eq.1) then
	            rhs_valdg(:,:,i) = zero
	            
	            
	            call dg_vol_integral_accumulate_rhs(n,iconsidered,-1.0d0)
            
            
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
! 						  if (realgas.eq.1)then
! 						  call fix_conservative_state(cleft_rot)
! 						  call fix_conservative_state(cright_rot)
!
! 						  end if
						  
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
                         

                         
	                        call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
                       
                         

                    else !fv
				       
				      
				      godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      
				      end if
				      
				       if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					    
                     norms=0.5*(cleft_rot(2)+cright_rot(2))
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(1:turbulenceequations+passivescalar)+cturbr(1:turbulenceequations+passivescalar)))+(abs(norms)*(cturbl(1:turbulenceequations+passivescalar)-(cturbr(1:turbulenceequations+passivescalar)))))
					  
					  
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

               
                  
               
	                if (multispecies.eq.1)then

					rhs_valdg(1,7,i)=rhs_valdg(1,7,i)-(u_c_val(1,7,i)*mp_source3)


                 
                 
                 
                 end if
                
              
                 end if
                 

end subroutine calculate_fluxeshi_convective2d_inner_cell


subroutine calculate_fluxeshi_convective2d_bound_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	integer::i,kmaxe
	real,dimension(1:gpu_max_nvar_total)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::l,ngp,iqp,ikas,igoflux,icaseb,kxk,b_code
	real::sum_detect,norms
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real::x1,y1,z1
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext




	i=el_bnd(ii)

	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero


	iconsidered=i	
	mp_source3=zero  
	
	            if (dg.eq.1) then
	            rhs_valdg(:,:,i)=zero
	            
	            call dg_vol_integral_accumulate_rhs(n,iconsidered,-1.0d0)
               
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
! 						  if (realgas.eq.1)then
! 						  call fix_conservative_state(cleft_rot)
! 						  call fix_conservative_state(cright_rot)
!
! 						  end if



				      
				      
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
                        
	                         call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,1.0d0)
                         
                         else
				      
                         godflux2(1:nof_variables)=godflux2(1:nof_variables)+(rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem_surf(l,i)))
				      end if
				      
				      
				      if (multispecies.eq.1)then
                        mp_source2=mp_source2+mp_source1*(weights_temp(ngp)*ielem_surf(l,i))
                        end if
				       
				      
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  norms=0.5*(cleft_rot(2)+cright_rot(2))
					  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(1:turbulenceequations+passivescalar)+cturbr(1:turbulenceequations+passivescalar)))+(abs(norms)*(cturbl(1:turbulenceequations+passivescalar)-(cturbr(1:turbulenceequations+passivescalar)))))
					  
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
	                if (multispecies.eq.1)then
                  rhs_valdg(1,7,i)=rhs_valdg(1,7,i)-(u_c_val(1,7,i)*mp_source3)

                 end if
                 
                 end if
                 
                 
                 

end subroutine calculate_fluxeshi_convective2d_bound_cell


subroutine calculate_fluxeshi_convective2d_inner_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz
real::angle1,angle2,nx,ny,nz,wface,norms

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, turbulenceequations, passivescalar, turbulence, icoupleturb) &
!$omp& firstprivate(lmach, lmach_style, iriemann, per_rot, angle_per, zero, oo2, gamma, qp_line_n) &
!$omp& map(alloc: weights_l, el_int, ielem_ifca, ielem_faceanglex, ielem_faceangley, ielem_surf) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh) &
!$omp& map(alloc: rec_uleft, rec_uleftturb, u_ct_val, rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#endif
do ii=1,nof_interior
i=el_int(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
rhs_val(1:nof_variables,i)=zero
if (nvt.gt.0) rhst_val(1:nvt,i)=zero
do l=1,ielem_ifca(i)
  facex=l
  b_code=0
  godflux2=zero
  nx=ielem_faceanglex(l,i)
  ny=ielem_faceangley(l,i)
  nz=zero
  angle1=nx
  angle2=ny
  iqp=qp_line_n
  do ngp=1,iqp
    pointx=ngp
    hllcflux=zero
    rhllcflux=zero
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    cturbl=zero
    cturbr=zero
    call get_states_interior2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    call rotatef2d(n,cright_rot,cright,angle1,angle2)
    call rotatef2d(n,cleft_rot,cleft,angle1,angle2)
    if (lmach.eq.1)then
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht2d_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
    end if
    call riemann_flux_ideal2d(n,cleft_rot,cright_rot,hllcflux,srf_speedrot)
    call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
    if (nvt.gt.0) rhllcflux(nof_variables+1:nof_variables+nvt)=hllcflux(nof_variables+1:nof_variables+nvt)
    if (nvt.gt.0)then
      if (icoupleturb.eq.0)then
        norms=0.5d0*(cleft_rot(2)+cright_rot(2))
        do iv=1,nvt
          rhllcflux(nof_variables+iv)=0.5d0*((norms*(cturbl(iv)+cturbr(iv)))+(abs(norms)*(cturbl(iv)-cturbr(iv))))
        end do
      end if
    end if
    wface=weights_l(ngp)*ielem_surf(l,i)
    do iv=1,nof_variables+nvt
      godflux2(iv)=godflux2(iv)+rhllcflux(iv)*wface
    end do
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)+godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_convective2d_inner_cell_ideal


subroutine calculate_fluxeshi_convective2d_bound_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz
real::angle1,angle2,nx,ny,nz,wface,norms

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, turbulenceequations, passivescalar, turbulence, icoupleturb) &
!$omp& firstprivate(lmach, lmach_style, iriemann, per_rot, angle_per, zero, oo2, gamma, qp_line_n) &
!$omp& firstprivate(r_gas, pres, rres, boundtype, initcond, Mach_in, uvel, vvel, wvel, press_outlet, swirl, tolsmall) &
!$omp& map(alloc: weights_l, el_bnd, ielem_ifca, ielem_faceanglex, ielem_faceangley, ielem_surf) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode, ielem_qface) &
!$omp& map(alloc: ielem_inter_id, ielem_indexf, bound_offset, boundhir, ielem_nodes_faces, inoder4_cord) &
!$omp& map(alloc: rec_uleft, rec_uleftturb, u_ct_val, rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,rhllcflux,hllcflux,cleft,cright,cleft_rot,cright_rot) &
!$omp& private(leftv,rightv,srf_speedrot,cturbl,cturbr,pox,poy,poz) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,norms)
#endif
do ii=1,nof_bounded
i=el_bnd(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
rhs_val(1:nof_variables,i)=zero
if (nvt.gt.0) rhst_val(1:nvt,i)=zero
do l=1,ielem_ifca(i)
  facex=l
  b_code=0
  godflux2=zero
  nx=ielem_faceanglex(l,i)
  ny=ielem_faceangley(l,i)
  nz=zero
  angle1=nx
  angle2=ny
  iqp=qp_line_n
  do ngp=1,iqp
    pointx=ngp
    hllcflux=zero
    rhllcflux=zero
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    cturbl=zero
    cturbr=zero
    call get_states_bounds2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    call rotatef2d(n,cright_rot,cright,angle1,angle2)
    call rotatef2d(n,cleft_rot,cleft,angle1,angle2)
    if (lmach.eq.1)then
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht2d_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
    end if
    call riemann_flux_ideal2d(n,cleft_rot,cright_rot,hllcflux,srf_speedrot)
    call rotateb2d(n,rhllcflux,hllcflux,angle1,angle2)
    if (nvt.gt.0) rhllcflux(nof_variables+1:nof_variables+nvt)=hllcflux(nof_variables+1:nof_variables+nvt)
    if (nvt.gt.0)then
      if (icoupleturb.eq.0)then
        norms=0.5d0*(cleft_rot(2)+cright_rot(2))
        do iv=1,nvt
          rhllcflux(nof_variables+iv)=0.5d0*((norms*(cturbl(iv)+cturbr(iv)))+(abs(norms)*(cturbl(iv)-cturbr(iv))))
        end do
      end if
      if ((b_code.eq.3).or.(b_code.eq.4)) rhllcflux(nof_variables+1:nof_variables+nvt)=zero
    end if
    wface=weights_l(ngp)*ielem_surf(l,i)
    do iv=1,nof_variables+nvt
      godflux2(iv)=godflux2(iv)+rhllcflux(iv)*wface
    end do
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)+godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_convective2d_bound_cell_ideal


	
	

subroutine calculate_fluxeshi_diffusive(n)
!> @brief
!> this subroutine computes the diffusive fluxes for euler-navier-stokes equations
		implicit none
		integer,intent(in)::n
		integer::ii
		logical::use_noturb_3d

		if ((dg.eq.0).and.(multispecies.eq.0).and.(realgas.eq.0))then
#ifdef xpu
!$omp barrier
!$omp master
#endif
		use_noturb_3d=.false.
#ifdef xpu
		if ((dimensiona.eq.3).and.(turbulence.eq.0).and.(passivescalar.eq.0).and.(mrf.eq.0))then
			use_noturb_3d=.true.
		end if
#endif
		if (use_noturb_3d)then
			call calculate_fluxeshi_diffusive_inner_cell_ideal_noturb_3d(n)
			call calculate_fluxeshi_diffusive_bound_cell_ideal_noturb_3d(n)
		else
			call calculate_fluxeshi_diffusive_inner_cell_ideal(n)
			call calculate_fluxeshi_diffusive_bound_cell_ideal(n)
		end if
#ifdef xpu
!$omp end master
!$omp barrier
#endif
		return
		end if
	
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	call calculate_fluxeshi_diffusive_inner_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	

	
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_diffusive_bound_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end subroutine calculate_fluxeshi_diffusive


subroutine calculate_fluxeshi_diffusive_inner_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,nvar,kc,iex,ittt,ikas,igoflux,icaseb,kk,b_code,srf,k,iv,nvt
	real::sum_detect,norms,wface,muturb
	integer::iconsidered,facex,pointx,rg_i,rg_j,idxy
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot,tempx_l,rtempx_l
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz,rg_sumfl,rg_sumfr
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real,dimension(1:8,1:gpu_max_dim)::vext
	real::mp_pinfl,mp_pinfr,gammal,gammar,rho12l,rho12r
	real,dimension(1:4)::viscl,laml
	real,dimension(1:2)::turbmv
    real,dimension(1)::etvm
    real,dimension(1:20)::eddyfl,eddyfr
    real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
	real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
	real,dimension(1:gpu_max_nvar)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(1:gpu_max_species)::rgs_htr_i, rgs_hvib_i,y_av,y_corr
	real,dimension(3,3)::taul,taur,tau
	real,dimension(3)::q,nnn,nall,qvib
	real::ux,uy,uz,vx,vy,vz,wx,wy,wz,rho12,u12,v12,w12 ,damp,vdamp,tempxx 
	real::mp_ttr,mp_tv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_laml,mp_lamr,y_face
	real,dimension(1:gpu_max_species)::mp_d_eff,rg_difl,rg_difr,rg_enthl,rg_enthr,mp_htr,mp_hvib,rg_enthvbl,rg_enthvbr,mp_mu_i,mp_ktr_i
	real,dimension(1:gpu_max_dim)::rg_sum_tr,rg_sumfr_tr,rg_sumfl_v,rg_sumfr_v,rg_sum_v,gradyl,gradyr,jl,jr,rg_sumfl_tr
	real,dimension(1:2)::qtr,qv
	real,dimension(1:gpu_max_dim)::rg_sum_htr,rg_sum_hv,grady,rg_qv,rg_qtr
	real,dimension(1:gpu_max_species)::rg_dif_av,rg_enth_av,rg_enthvb_av
	real::mp_lam_av,mp_ktr_mix_av,sum_y1,sum_y2,sumY_face,molar_sum,mbar
	real,dimension(1:gpu_max_dim):: sumi,gradx,grad_a,sumJ,sumGradY
	real,dimension(1:gpu_max_dim,1:gpu_max_species) :: i_raw,gradY_all


	i=el_int(ii)
	iconsidered=i
	nvt=turbulenceequations+passivescalar

	                    
			   damp=lamx
            if( br2_yn == 2) damp = 0.0d0
		   if(mrf.eq.1)then
                srf=rec_mrf(i)
            end if
		    do l=1,ielem_ifca(i) !for all their faces

				  do iv=1,nof_variables+nvt
				  godflux2(iv)=zero
				  end do
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
								call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
								end if
								if (turbulencemodel.eq.2)then
								eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
								eddyfl(4:6)= lcvgrad(1,1:3);eddyfl(7:9)=lcvgrad(2,1:3)
								eddyfl(10:12)=lcvgrad(3,1:3);eddyfl(13:15)=lcvgrad_t(1,1:3)
								eddyfl(16:18)=lcvgrad_t(2,1:3)


								eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
								eddyfr(4:6)= rcvgrad(1,1:3);eddyfr(7:9)=rcvgrad(2,1:3);eddyfr(10:12)=rcvgrad(3,1:3)
								eddyfr(13:15)=rcvgrad_t(1,1:3);eddyfl(16:18)=rcvgrad_t(2,1:3)
									call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
								end if
						end if
				       


						taul = zero;tau=zero;taur=zero;q=zero;ux=zero;uy=zero;uz=zero;vx=zero;vy=zero;vz=zero;wx=zero;wy=zero;wz=zero;
					    do iv=1,nof_variables
					    fxv(iv)=zero
					    fyv(iv)=zero
					    fzv(iv)=zero
					    end do
					    rho12 =zero;
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
					do iv=1,3
					lcvgrad(k,iv)=((lcvgrad(k,iv)+rcvgrad(k,iv))/(2.0d0))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(iv)*(rightv(k+1)-leftv(k+1)))
					end do
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

													sumi(1:dimensiona) = 0.0d0
													sumJ(1:dimensiona) = 0.0d0
													i_raw(1:dimensiona,1:nof_species) = 0.0d0

													! normalise face mass fractions used by correction velocity
													sumY_face = 0.0d0
													do rg_j = 1, nof_species
														y_corr(rg_j) = max(y_av(rg_j),0.0d0)
														sumY_face = sumY_face + y_corr(rg_j)
													end do

													if (sumY_face .gt. 1.0d-30) then
														do rg_j = 1, nof_species
															y_corr(rg_j) = y_corr(rg_j) / sumY_face
														end do
													else
														do rg_j = 1, nof_species
															y_corr(rg_j) = 1.0d0 / dble(nof_species)
														end do
													end if

													! mixture molar sum from corrected face mass fractions
													molar_sum = 0.0d0
													do rg_j = 1, nof_species
														molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
													end do
													molar_sum = max(molar_sum,1.0d-30)
													mbar = 1.0d0 / molar_sum

													! project species gradients so sum_k grad(Y_k)=0
													sumGradY(1:dimensiona) = 0.0d0
													do rg_j = 1, nof_species
														idxy = dimensiona + 2 + rg_j
														gradY_all(1:dimensiona,rg_j) = lcvgrad(idxy,1:dimensiona)
														sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona,rg_j)
													end do
													do rg_j = 1, nof_species
														gradY_all(1:dimensiona,rg_j) = gradY_all(1:dimensiona,rg_j) - &
															y_corr(rg_j) * sumGradY(1:dimensiona)
													end do

													grad_a(1:dimensiona) = 0.0d0
													do rg_j = 1, nof_species
														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_j)
														grad_a(1:dimensiona) = grad_a(1:dimensiona) + grady(1:dimensiona) / rg_molm(rg_j)
													end do

													! raw mixture-averaged flux from grad(X_i)
													do rg_i = 1, nof_species
														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_i)
														gradx(1:dimensiona) = &
															( grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
															- (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona) ) &
															/ (molar_sum * molar_sum)

														i_raw(1:dimensiona,rg_i) = &
															-rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)
														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)
													end do

													! corrected species flux, conservative because sum(y_corr)=1
													do rg_i = 1, nof_species
														y_face = y_corr(rg_i)
														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - y_face * sumi(1:dimensiona)
														sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)
														fzv(dimensiona+3+rg_i) = fzv(dimensiona+3+rg_i) + jl(3)

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




		
		
					  do iv=1,nof_variables
					  hllcflux(iv)=nx*fxv(iv)+ny*fyv(iv)+nz*fzv(iv)
					  end do


					  if (dg.eq.1)then
					  do iv=1,nof_variables
					  rhllcflux(iv)=hllcflux(iv)
					  end do

						  call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,-1.0d0)

					  else
					  
				      wface=weights_temp(ngp)*ielem_surf(l,i)
				      do iv=1,nof_variables
				      godflux2(iv)=godflux2(iv)+(hllcflux(iv)*wface)
				      end do

				      end if
				      
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					  if (turbulence.eq.1)then
					  muturb=(oo2*(viscl(1)+viscl(2)))+(oo2*(viscl(3)+viscl(4)))
					  else
					  muturb=oo2*(viscl(1)+viscl(2))
					  end if
					  do iv=1,nvt
					  hllcflux(nof_variables+iv)=muturb*oo2*((lcvgrad_t(iv,1)+rcvgrad_t(iv,1))*nx+&
					  (lcvgrad_t(iv,2)+rcvgrad_t(iv,2))*ny+(lcvgrad_t(iv,3)+rcvgrad_t(iv,3))*nz)
					  end do
						if (turbulencemodel.eq.1)then
						do iv=1,nvt
						hllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)/sigma
						end do
						end if					  
					  wface=weights_temp(ngp)*ielem_surf(l,i)
					  do iv=1,nvt
					  godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+(hllcflux(nof_variables+iv)*wface)
					  end do
				      end if
				      
				      
				  end do
				  
				    do iv=1,nof_variables
				    rhs_val(iv,i)=rhs_val(iv,i)-godflux2(iv)
				    end do

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    do iv=1,nvt
				    rhst_val(iv,i)=rhst_val(iv,i)-godflux2(nof_variables+iv)
				    end do

 				    end if

		    end do

	end subroutine calculate_fluxeshi_diffusive_inner_cell


subroutine calculate_fluxeshi_diffusive_bound_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,nvar,kc,iex,ittt,ikas,igoflux,icaseb,kk,b_code,srf,k,iv,nvt
	real::sum_detect,norms,wface,muturb
	integer::iconsidered,facex,pointx,rg_i,rg_j,idxy
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot,tempx_l,rtempx_l
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz,rg_sumfl,rg_sumfr
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real,dimension(1:8,1:gpu_max_dim)::vext
	real::mp_pinfl,mp_pinfr,gammal,gammar,rho12l,rho12r
	real,dimension(1:4)::viscl,laml
	real,dimension(1:2)::turbmv
    real,dimension(1)::etvm
    real,dimension(1:20)::eddyfl,eddyfr
    real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
	real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
	real,dimension(1:gpu_max_nvar)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(1:gpu_max_species)::rgs_htr_i, rgs_hvib_i,y_av,y_corr
	real,dimension(3,3)::taul,taur,tau
	real,dimension(3)::q,nnn,nall,qvib
	real::ux,uy,uz,vx,vy,vz,wx,wy,wz,rho12,u12,v12,w12 ,damp,vdamp,tempxx 
	real::mp_ttr,mp_tv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_laml,mp_lamr,y_face
	real,dimension(1:gpu_max_species)::mp_d_eff,rg_difl,rg_difr,rg_enthl,rg_enthr,mp_htr,mp_hvib,rg_enthvbl,rg_enthvbr,mp_mu_i,mp_ktr_i
	real,dimension(1:gpu_max_dim)::rg_sum_tr,rg_sumfr_tr,rg_sumfl_v,rg_sumfr_v,rg_sum_v,gradyl,gradyr,jl,jr,rg_sumfl_tr
	real,dimension(1:2)::qtr,qv
	real,dimension(1:gpu_max_dim)::rg_sum_htr,rg_sum_hv,grady,rg_qv,rg_qtr
	real,dimension(1:gpu_max_species)::rg_dif_av,rg_enth_av,rg_enthvb_av
	real::mp_lam_av,mp_ktr_mix_av,sum_y1,sum_y2,sumY_face,molar_sum,mbar
	real,dimension(1:gpu_max_dim):: sumi,gradx,grad_a,sumJ,sumGradY
	real,dimension(1:gpu_max_dim,1:gpu_max_species) :: i_raw,gradY_all


	i=el_bnd(ii)
	iconsidered=i	
	nvt=turbulenceequations+passivescalar
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
				  

				  
					  do iv=1,nof_variables+nvt
					  godflux2(iv)=zero
					  end do
				  
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
							  call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
						      if (turbulencemodel.eq.2)then
							  eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
							  eddyfl(4:6)= lcvgrad(1,1:3);eddyfl(7:9)=lcvgrad(2,1:3)
							  eddyfl(10:12)=lcvgrad(3,1:3);eddyfl(13:15)=lcvgrad_t(1,1:3)
							  eddyfl(16:18)=lcvgrad_t(2,1:3)
							    
							    
							  eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
							  eddyfr(4:6)= rcvgrad(1,1:3);eddyfr(7:9)=rcvgrad(2,1:3);eddyfr(10:12)=rcvgrad(3,1:3)
							  eddyfr(13:15)=rcvgrad_t(1,1:3);eddyfl(16:18)=rcvgrad_t(2,1:3)
							    call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
					  end if


                                            
					  taul = zero;tau=zero;taur=zero;q=zero;ux=zero;uy=zero;uz=zero;vx=zero;vy=zero;vz=zero;wx=zero;wy=zero;wz=zero;
					    do iv=1,nof_variables
					    fxv(iv)=zero
					    fyv(iv)=zero
					    fzv(iv)=zero
					    end do
					    rho12 =zero;
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
											do iv=1,3
											lcvgrad(k,iv)=((lcvgrad(k,iv)+rcvgrad(k,iv))/(2.0d0))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(iv)*(rightv(k+1)-leftv(k+1)))
											end do
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

													sumi(1:dimensiona) = 0.0d0
													sumJ(1:dimensiona) = 0.0d0
													i_raw(1:dimensiona,1:nof_species) = 0.0d0

													! normalise face mass fractions used by correction velocity
													sumY_face = 0.0d0
													do rg_j = 1, nof_species
														y_corr(rg_j) = max(y_av(rg_j),0.0d0)
														sumY_face = sumY_face + y_corr(rg_j)
													end do

													if (sumY_face .gt. 1.0d-30) then
														do rg_j = 1, nof_species
															y_corr(rg_j) = y_corr(rg_j) / sumY_face
														end do
													else
														do rg_j = 1, nof_species
															y_corr(rg_j) = 1.0d0 / dble(nof_species)
														end do
													end if

													! mixture molar sum from corrected face mass fractions
													molar_sum = 0.0d0
													do rg_j = 1, nof_species
														molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
													end do
													molar_sum = max(molar_sum,1.0d-30)
													mbar = 1.0d0 / molar_sum

													! project species gradients so sum_k grad(Y_k)=0
													sumGradY(1:dimensiona) = 0.0d0
													do rg_j = 1, nof_species
														idxy = dimensiona + 2 + rg_j
														gradY_all(1:dimensiona,rg_j) = lcvgrad(idxy,1:dimensiona)
														sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona,rg_j)
													end do
													do rg_j = 1, nof_species
														gradY_all(1:dimensiona,rg_j) = gradY_all(1:dimensiona,rg_j) - &
															y_corr(rg_j) * sumGradY(1:dimensiona)
													end do

													grad_a(1:dimensiona) = 0.0d0
													do rg_j = 1, nof_species
														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_j)
														grad_a(1:dimensiona) = grad_a(1:dimensiona) + grady(1:dimensiona) / rg_molm(rg_j)
													end do

													! raw mixture-averaged flux from grad(X_i)
													do rg_i = 1, nof_species
														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_i)
														gradx(1:dimensiona) = &
															( grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
															- (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona) ) &
															/ (molar_sum * molar_sum)

														i_raw(1:dimensiona,rg_i) = &
															-rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)
														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)
													end do

													! corrected species flux, conservative because sum(y_corr)=1
													do rg_i = 1, nof_species
														y_face = y_corr(rg_i)
														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - y_face * sumi(1:dimensiona)
														sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)
														fzv(dimensiona+3+rg_i) = fzv(dimensiona+3+rg_i) + jl(3)

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
													fzv(5) = fzv(5) - q(3)

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




		
		
					  do iv=1,nof_variables
					  hllcflux(iv)=nx*fxv(iv)+ny*fyv(iv)+nz*fzv(iv)
					  end do


					  if (realgas.eq.1)then
						if ((b_code.eq.4).and.(catalytic_wall.eq.0))then
						do iv=7,nof_variables
						hllcflux(iv)=zero
						end do
						end if
						end if


					   if (dg.eq.1)then
					   do iv=1,nof_variables
					   rhllcflux(iv)=hllcflux(iv)
					   end do

						  call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,-1.0d0)


					  else
				      wface=weights_temp(ngp)*ielem_surf(l,i)
				      do iv=1,nof_variables
				      godflux2(iv)=godflux2(iv)+(hllcflux(iv)*wface)
				      end do
				      end if
				      
				      
				     
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then 
					   if (turbulence.eq.1)then
					  muturb=(oo2*(viscl(1)+viscl(2)))+(oo2*(viscl(3)+viscl(4)))
					  else
					  muturb=oo2*(viscl(1)+viscl(2))
					  end if
					  do iv=1,nvt
					  hllcflux(nof_variables+iv)=muturb*oo2*((lcvgrad_t(iv,1)+rcvgrad_t(iv,1))*nx+&
					  (lcvgrad_t(iv,2)+rcvgrad_t(iv,2))*ny+(lcvgrad_t(iv,3)+rcvgrad_t(iv,3))*nz)
					  end do
						if (turbulencemodel.eq.1)then
						do iv=1,nvt
						hllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)/sigma
						end do
						end if					  

					if ((b_code.eq.3))then
					  do iv=1,nvt
					  hllcflux(nof_variables+iv)=zero
					  end do
					  
					  end if




					  wface=weights_temp(ngp)*ielem_surf(l,i)
					  do iv=1,nvt
					  godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+(hllcflux(nof_variables+iv)*wface)
					  end do
				      end if
				      
				      
				  end do
				  
				    
				    do iv=1,nof_variables
				    rhs_val(iv,i)=rhs_val(iv,i)-godflux2(iv)
				    end do

				    
				    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				    do iv=1,nvt
				    rhst_val(iv,i)=rhst_val(iv,i)-godflux2(nof_variables+iv)
				    end do

				    end if
		    end do

	end subroutine calculate_fluxeshi_diffusive_bound_cell





#if 1
subroutine calculate_fluxeshi_diffusive_inner_cell_ideal_noturb_3d(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,k
real,dimension(1:gpu_max_nvar)::cleft,cright,cleft_rot,cright_rot,leftv,rightv
real,dimension(1:4,1:3)::lgrad,rgrad,grad
real::angle1,angle2,nx,ny,nz,wface,damp,vdamp
real::god2,god3,god4,god5,flux2,flux3,flux4,flux5
real::rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r,t_l,t_r,t0,mu_l,mu_r,mu_av,lam_l,lam_r,lam_av
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12
real::tau11,tau22,tau33,tau12,tau13,tau23,q1,q2,q3,fx5,fy5,fz5
real::ca1,sa1,ca2,sa2,ssl,ssr,q2l,q2r,mach2,mach,dus,dvs,dws,duu,dvv,dww,diff
real::eel,eer,ppl,ppr

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, lmach, lmach_style, zero, oo2) &
!$omp& firstprivate(gamma, r_gas, pres, rres, visc, suther, prandtl, lamx, br2_yn) &
!$omp& firstprivate(qp_quad_n, qp_triangle_n) &
!$omp& map(alloc: weights_q, weights_t, el_int, ielem_ifca, ielem_types_faces) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_dih) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, rec_uleft, rec_uleftv, rhs_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,k,cleft,cright,cleft_rot,cright_rot,leftv,rightv) &
!$omp& private(lgrad,rgrad,grad,angle1,angle2,nx,ny,nz,wface,damp,vdamp) &
!$omp& private(god2,god3,god4,god5,flux2,flux3,flux4,flux5) &
!$omp& private(rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r,t_l,t_r,t0,mu_l,mu_r,mu_av,lam_l,lam_r,lam_av) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,tau11,tau22,tau33,tau12,tau13,tau23) &
!$omp& private(q1,q2,q3,fx5,fy5,fz5,ca1,sa1,ca2,sa2,ssl,ssr,q2l,q2r,mach2,mach,dus,dvs,dws) &
!$omp& private(duu,dvv,dww,diff,eel,eer,ppl,ppr)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,k,cleft,cright,cleft_rot,cright_rot,leftv,rightv) &
!$omp& private(lgrad,rgrad,grad,angle1,angle2,nx,ny,nz,wface,damp,vdamp) &
!$omp& private(god2,god3,god4,god5,flux2,flux3,flux4,flux5) &
!$omp& private(rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r,t_l,t_r,t0,mu_l,mu_r,mu_av,lam_l,lam_r,lam_av) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,tau11,tau22,tau33,tau12,tau13,tau23) &
!$omp& private(q1,q2,q3,fx5,fy5,fz5,ca1,sa1,ca2,sa2,ssl,ssr,q2l,q2r,mach2,mach,dus,dvs,dws) &
!$omp& private(duu,dvv,dww,diff,eel,eer,ppl,ppr)
#endif
do ii=1,nof_interior
i=el_int(ii)
damp=lamx
if (br2_yn.eq.2) damp=zero
do l=1,ielem_ifca(i)
  god2=zero
  god3=zero
  god4=zero
  god5=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  ca1=cos(angle1)
  sa1=sin(angle1)
  ca2=cos(angle2)
  sa2=sin(angle2)
  nx=ca1*sa2
  ny=sa1*sa2
  nz=ca2
  if (ielem_types_faces(l,i).eq.5)then
    iqp=qp_quad_n
  else
    iqp=qp_triangle_n
  end if
  do ngp=1,iqp
    cleft(1:5)=rec_uleft(1:5,l,ngp,i)
    cright(1:5)=rec_uleft(1:5,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    do iv=1,3
      do k=1,4
        lgrad(k,iv)=rec_uleftv(iv,k,l,ngp,i)
        rgrad(k,iv)=rec_uleftv(iv,k,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      end do
    end do
    if (lmach.eq.1)then
      cleft_rot(1)=cleft(1)
      cright_rot(1)=cright(1)
      cleft_rot(2)=nx*cleft(2)+ny*cleft(3)+nz*cleft(4)
      cright_rot(2)=nx*cright(2)+ny*cright(3)+nz*cright(4)
      cleft_rot(3)=ca1*ca2*cleft(2)+sa1*ca2*cleft(3)-sa2*cleft(4)
      cright_rot(3)=ca1*ca2*cright(2)+sa1*ca2*cright(3)-sa2*cright(4)
      cleft_rot(4)=-sa1*cleft(2)+ca1*cleft(3)
      cright_rot(4)=-sa1*cright(2)+ca1*cright(3)
      cleft_rot(5)=cleft(5)
      cright_rot(5)=cright(5)
      rho_l=cleft_rot(1)
      rho_r=cright_rot(1)
      u_l=cleft_rot(2)/rho_l
      u_r=cright_rot(2)/rho_r
      v_l=cleft_rot(3)/rho_l
      v_r=cright_rot(3)/rho_r
      w_l=cleft_rot(4)/rho_l
      w_r=cright_rot(4)/rho_r
      eel=cleft_rot(5)
      eer=cright_rot(5)
      ppl=(gamma-1.0d0)*(eel-oo2*rho_l*(u_l*u_l+v_l*v_l+w_l*w_l))
      ppr=(gamma-1.0d0)*(eer-oo2*rho_r*(u_r*u_r+v_r*v_r+w_r*w_r))
      ssl=(gamma*ppl)/rho_l
      ssr=(gamma*ppr)/rho_r
      q2l=u_l*u_l+v_l*v_l+w_l*w_l
      q2r=u_r*u_r+v_r*v_r+w_r*w_r
      mach2=max(q2l/ssl,q2r/ssr)
      mach=sqrt(mach2)
      mach=min(mach,1.0d0)
      dus=u_r+u_l
      dvs=v_r+v_l
      dws=w_r+w_l
      duu=mach*(u_r-u_l)
      dvv=v_r-v_l
      dww=w_r-w_l
      diff=oo2*duu
      u_l=oo2*dus-diff
      u_r=oo2*dus+diff
      if (lmach_style.eq.1)then
        diff=oo2*dvv
        v_l=oo2*dvs-diff
        v_r=oo2*dvs+diff
        diff=oo2*dww
        w_l=oo2*dws-diff
        w_r=oo2*dws+diff
      end if
      cleft_rot(2)=rho_l*u_l
      cright_rot(2)=rho_r*u_r
      cleft_rot(3)=rho_l*v_l
      cright_rot(3)=rho_r*v_r
      cleft_rot(4)=rho_l*w_l
      cright_rot(4)=rho_r*w_r
      cleft(2)=nx*cleft_rot(2)+ca1*ca2*cleft_rot(3)-sa1*cleft_rot(4)
      cright(2)=nx*cright_rot(2)+ca1*ca2*cright_rot(3)-sa1*cright_rot(4)
      cleft(3)=ny*cleft_rot(2)+sa1*ca2*cleft_rot(3)+ca1*cleft_rot(4)
      cright(3)=ny*cright_rot(2)+sa1*ca2*cright_rot(3)+ca1*cright_rot(4)
      cleft(4)=nz*cleft_rot(2)-sa2*cleft_rot(3)
      cright(4)=nz*cright_rot(2)-sa2*cright_rot(3)
    end if
    rho_l=cleft(1)
    rho_r=cright(1)
    u_l=cleft(2)/rho_l
    v_l=cleft(3)/rho_l
    w_l=cleft(4)/rho_l
    u_r=cright(2)/rho_r
    v_r=cright(3)/rho_r
    w_r=cright(4)/rho_r
    t_l=((gamma-1.0d0)*(cleft(5)-oo2*rho_l*(u_l*u_l+v_l*v_l+w_l*w_l)))/(rho_l*r_gas)
    t_r=((gamma-1.0d0)*(cright(5)-oo2*rho_r*(u_r*u_r+v_r*v_r+w_r*w_r)))/(rho_r*r_gas)
    t0=pres/(rres*r_gas)
    mu_l=visc*((t_l/t0)*sqrt(t_l/t0))*((t0+suther*t0)/(t_l+suther*t0))
    mu_r=visc*((t_r/t0)*sqrt(t_r/t0))*((t0+suther*t0)/(t_r+suther*t0))
    lam_l=mu_l*r_gas*gamma/(prandtl*(gamma-1.0d0))
    lam_r=mu_r*r_gas*gamma/(prandtl*(gamma-1.0d0))
    mu_av=oo2*(mu_l+mu_r)
    lam_av=oo2*(lam_l+lam_r)
    leftv(2)=u_l
    leftv(3)=v_l
    leftv(4)=w_l
    leftv(5)=t_l
    rightv(2)=u_r
    rightv(3)=v_r
    rightv(4)=w_r
    rightv(5)=t_r
    vdamp=4.0d0/3.0d0
    do k=1,4
      grad(k,1)=oo2*(lgrad(k,1)+rgrad(k,1))+damp*((vdamp/abs(ielem_dih(l,i)))*nx*(rightv(k+1)-leftv(k+1)))
      grad(k,2)=oo2*(lgrad(k,2)+rgrad(k,2))+damp*((vdamp/abs(ielem_dih(l,i)))*ny*(rightv(k+1)-leftv(k+1)))
      grad(k,3)=oo2*(lgrad(k,3)+rgrad(k,3))+damp*((vdamp/abs(ielem_dih(l,i)))*nz*(rightv(k+1)-leftv(k+1)))
    end do
    q1=-lam_av*grad(4,1)
    q2=-lam_av*grad(4,2)
    q3=-lam_av*grad(4,3)
    ux=grad(1,1)
    uy=grad(1,2)
    uz=grad(1,3)
    vx=grad(2,1)
    vy=grad(2,2)
    vz=grad(2,3)
    wx=grad(3,1)
    wy=grad(3,2)
    wz=grad(3,3)
    tau11=mu_av*((4.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy-(2.0d0/3.0d0)*wz)
    tau22=mu_av*((4.0d0/3.0d0)*vy-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*wz)
    tau33=mu_av*((4.0d0/3.0d0)*wz-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy)
    tau12=mu_av*(uy+vx)
    tau13=mu_av*(wx+uz)
    tau23=mu_av*(vz+wy)
    flux2=nx*tau11+ny*tau12+nz*tau13
    flux3=nx*tau12+ny*tau22+nz*tau23
    flux4=nx*tau13+ny*tau23+nz*tau33
    u12=oo2*(leftv(2)+rightv(2))
    v12=oo2*(leftv(3)+rightv(3))
    w12=oo2*(leftv(4)+rightv(4))
    fx5=-q1+u12*tau11+v12*tau12+w12*tau13
    fy5=-q2+u12*tau12+v12*tau22+w12*tau23
    fz5=-q3+u12*tau13+v12*tau23+w12*tau33
    flux5=nx*fx5+ny*fy5+nz*fz5
    if (ielem_types_faces(l,i).eq.5)then
      wface=weights_q(ngp)*ielem_surf(l,i)
    else
      wface=weights_t(ngp)*ielem_surf(l,i)
    end if
    god2=god2+flux2*wface
    god3=god3+flux3*wface
    god4=god4+flux4*wface
    god5=god5+flux5*wface
  end do
  rhs_val(2,i)=rhs_val(2,i)-god2
  rhs_val(3,i)=rhs_val(3,i)-god3
  rhs_val(4,i)=rhs_val(4,i)-god4
  rhs_val(5,i)=rhs_val(5,i)-god5
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_diffusive_inner_cell_ideal_noturb_3d


#if 1
subroutine calculate_fluxeshi_diffusive_bound_cell_ideal_noturb_3d(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,k,code_per,nfx,lfx,rowfx,ittt
integer::b_code,n_node,iex,d,jj
real,dimension(1:gpu_max_nvar)::cleft,cright,cleft_rot,cright_rot,leftv,rightv
real,dimension(1:4,1:3)::lgrad,rgrad,grad
real,dimension(1:3)::rot_tmp
real,dimension(1:gpu_max_nvar)::subson1,subson2,subson3
real,dimension(1:gpu_max_dim)::pox,poy,poz,cords,nnt,g
real,dimension(1:8,1:gpu_max_dim)::vext,nodes_list
real,dimension(1:3,1:3)::mat_a,mat_b,mat_c,mat_r
real::angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp
real::god2,god3,god4,god5,flux2,flux3,flux4,flux5
real::rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r,t_l,t_r,t0,mu_l,mu_r,mu_av,lam_l,lam_r,lam_av
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12
real::tau11,tau22,tau33,tau12,tau13,tau23,q1,q2,q3,fx5,fy5,fz5
real::ca1,sa1,ca2,sa2,ssl,ssr,q2l,q2r,mach2,mach,dus,dvs,dws,duu,dvv,dww,diff
real::eel,eer,ppl,ppr,g_n,mp_pinfl,mp_pinfr,gammal,gammar,sps,vel,vnb

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, lmach, lmach_style, per_rot, angle_per, zero, oo2) &
!$omp& firstprivate(gamma, r_gas, pres, rres, visc, suther, prandtl, lamx, br2_yn) &
!$omp& firstprivate(thermal, boundtype, initcond, Mach_in, uvel, vvel, wvel, press_outlet, swirl, tolsmall) &
!$omp& firstprivate(qp_quad_n, qp_triangle_n) &
!$omp& map(alloc: weights_q, weights_t, el_bnd, ielem_ifca, ielem_types_faces) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_dih) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode, ielem_qface) &
!$omp& map(alloc: ielem_inter_id, ielem_indexf, bound_offset, boundhir, rec_uleft, rec_uleftv, rhs_val) &
!$omp& map(alloc: ielem_nodes_faces, inoder4_cord) &
!$omp& private(ii,i,l,ngp,iqp,iv,k,code_per,nfx,lfx,rowfx,ittt,b_code,n_node,iex,d,jj) &
!$omp& private(cleft,cright,cleft_rot,cright_rot,leftv,rightv,subson1,subson2,subson3) &
!$omp& private(lgrad,rgrad,grad,rot_tmp,pox,poy,poz,cords,nnt,g,vext,nodes_list,mat_a,mat_b,mat_c,mat_r) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,g_n) &
!$omp& private(god2,god3,god4,god5,flux2,flux3,flux4,flux5) &
!$omp& private(rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r,t_l,t_r,t0,mu_l,mu_r,mu_av,lam_l,lam_r,lam_av) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,tau11,tau22,tau33,tau12,tau13,tau23) &
!$omp& private(q1,q2,q3,fx5,fy5,fz5,ca1,sa1,ca2,sa2,ssl,ssr,q2l,q2r,mach2,mach,dus,dvs,dws) &
!$omp& private(duu,dvv,dww,diff,eel,eer,ppl,ppr,mp_pinfl,mp_pinfr,gammal,gammar,sps,vel,vnb)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,k,code_per,nfx,lfx,rowfx,ittt,b_code,n_node,iex,d,jj) &
!$omp& private(cleft,cright,cleft_rot,cright_rot,leftv,rightv,subson1,subson2,subson3) &
!$omp& private(lgrad,rgrad,grad,rot_tmp,pox,poy,poz,cords,nnt,g,vext,nodes_list,mat_a,mat_b,mat_c,mat_r) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,g_n) &
!$omp& private(god2,god3,god4,god5,flux2,flux3,flux4,flux5) &
!$omp& private(rho_l,rho_r,u_l,u_r,v_l,v_r,w_l,w_r,t_l,t_r,t0,mu_l,mu_r,mu_av,lam_l,lam_r,lam_av) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,tau11,tau22,tau33,tau12,tau13,tau23) &
!$omp& private(q1,q2,q3,fx5,fy5,fz5,ca1,sa1,ca2,sa2,ssl,ssr,q2l,q2r,mach2,mach,dus,dvs,dws) &
!$omp& private(duu,dvv,dww,diff,eel,eer,ppl,ppr,mp_pinfl,mp_pinfr,gammal,gammar,sps,vel,vnb)
#endif
do ii=1,nof_bounded
i=el_bnd(ii)
do l=1,ielem_ifca(i)
  god2=zero
  god3=zero
  god4=zero
  god5=zero
  damp=lamx
  if (br2_yn.eq.2) damp=zero
  damp_loc=damp
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  ca1=cos(angle1)
  sa1=sin(angle1)
  ca2=cos(angle2)
  sa2=sin(angle2)
  nx=ca1*sa2
  ny=sa1*sa2
  nz=ca2
  if (ielem_types_faces(l,i).eq.5)then
    iqp=qp_quad_n
  else
    iqp=qp_triangle_n
  end if
  do ngp=1,iqp
    cleft(1:5)=rec_uleft(1:5,l,ngp,i)
    do iv=1,3
      do k=1,4
        lgrad(k,iv)=rec_uleftv(iv,k,l,ngp,i)
      end do
    end do
    b_code=0
    cleft_rot=zero
    cright_rot=zero
    subson1=zero
    subson2=zero
    subson3=zero
    code_per=0
    if (ielem_ibounds(l,i).gt.0) code_per=ibound_icode(ielem_ibounds(l,i))
    if (ielem_ineighb(l,i).eq.n)then
      if ((ielem_ibounds(l,i).gt.0).and. &
          & (.not.((code_per.eq.5).or.(code_per.eq.50))))then
        b_code=ibound_icode(ielem_ibounds(l,i))
        leftv(1:5)=cleft(1:5)
        rightv=zero
        call coordinates_face_innerx(n,i,l,vext,nodes_list)
        if (ielem_types_faces(l,i).eq.5)then
          n_node=4
        else
          n_node=3
        end if
        call cordinates3(n,nodes_list,n_node,cords(1:3))
        pox=zero
        poy=zero
        poz=zero
        pox(1)=cords(1)
        poy(1)=cords(2)
        poz(1)=cords(3)
	        select case(b_code)
	        case(1)
	          if (initcond.eq.4440)then
	            call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
	          else
	            call inflow_ideal3d(initcond,pox,poy,poz,rightv)
	            if (boundtype.ne.0)then
	            call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
	            subson1(1:5)=rightv(1:5)
	            subson2(1:5)=leftv(1:5)
            sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
            vel=sqrt((subson2(2)*subson2(2))+(subson2(3)*subson2(3))+(subson2(4)*subson2(4)))
            call prim2cons2_ideal(n,leftv,rightv)
            if (vel/(sps+tolsmall).le.1.0d0)then
              subson3(5)=0.5d0*(subson1(5)+subson2(5)-subson2(1)*sps* &
                   & (nx*(subson1(2)-subson2(2))+ny*(subson1(3)-subson2(3))+nz*(subson1(4)-subson2(4))))
              subson3(1)=subson1(1)+(subson3(5)-subson1(5))/(sps*sps+tolsmall)
              subson3(2)=subson1(2)-nx*(subson1(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
              subson3(3)=subson1(3)-ny*(subson1(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
              subson3(4)=subson1(4)-nz*(subson1(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
	              call set_state_ideal3d(subson3(1),subson3(2),subson3(3),subson3(4),subson3(5),rightv)
	            end if
	          end if
	          end if
        case(2)
          if (boundtype.eq.0)then
            rightv(1:5)=leftv(1:5)
          else
            call outflow_ideal3d(initcond,pox,poy,poz,rightv)
            call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
            subson1(1:5)=rightv(1:5)
            subson2(1:5)=leftv(1:5)
            sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
            vel=subson2(2)*nx+subson2(3)*ny+subson2(4)*nz
            call prim2cons2_ideal(n,leftv,rightv)
            if (vel.gt.0.0d0)then
              if (vel/(sps+tolsmall).gt.1.0d0)then
                rightv(1:5)=leftv(1:5)
              else
                subson3(1:5)=subson2(1:5)
                subson3(5)=subson1(5)
                subson3(1)=subson2(1)+(subson3(5)-subson2(5))/(sps*sps+tolsmall)
                subson3(2)=subson2(2)+nx*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
                subson3(3)=subson2(3)+ny*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
                subson3(4)=subson2(4)+nz*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
                call set_state_ideal3d(subson3(1),subson3(2),subson3(3),subson3(4),subson3(5),rightv)
              end if
            end if
          end if
        case(9)
          call outlet_pressure_fix_ideal3d(n,leftv,rightv,nx,ny,nz,press_outlet)
        case(3)
          call rotatef(n,cleft_rot,leftv,angle1,angle2)
          cright_rot(1)=cleft_rot(1)
          cright_rot(2)=-cleft_rot(2)
          cright_rot(3)=cleft_rot(3)
          cright_rot(4)=cleft_rot(4)
          cright_rot(5)=cleft_rot(5)
          call rotateb(n,rightv,cright_rot,angle1,angle2)
        case(4)
          rightv(1:5)=leftv(1:5)
          rightv(2)=-leftv(2)
          rightv(3)=-leftv(3)
          rightv(4)=-leftv(4)
        case(6)
          call rotatef(n,cleft_rot,leftv,angle1,angle2)
          vnb=cleft_rot(2)
          call outflow_ideal3d(initcond,pox,poy,poz,rightv)
          call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
          subson1(1:5)=rightv(1:5)
          subson2(1:5)=leftv(1:5)
          sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
          call prim2cons2_ideal(n,leftv,rightv)
          if (vnb.le.0.0d0)then
            if (initcond.eq.4440)then
              call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
            else
              call inflow_ideal3d(initcond,pox,poy,poz,rightv)
            end if
          else
            if ((abs(vnb)).ge.sps)then
              rightv(1:5)=leftv(1:5)
            else
              subson3(5)=subson1(5)
              subson3(1)=subson2(1)+(subson3(5)-subson2(5))/(sps*sps+tolsmall)
              subson3(2)=subson2(2)+nx*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
              subson3(3)=subson2(3)+ny*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
              subson3(4)=subson2(4)+nz*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
              call set_state_ideal3d(subson3(1),subson3(2),subson3(3),subson3(4),subson3(5),rightv)
            end if
          end if
        case default
          rightv(1:5)=leftv(1:5)
        end select
        cright(1:5)=rightv(1:5)
        rgrad(1:4,1:3)=lgrad(1:4,1:3)
        damp_loc=zero
        if (((b_code.eq.4).or.(b_code.eq.2)).and.(thermal.ne.1))then
          nnt(1)=nx
          nnt(2)=ny
          nnt(3)=nz
          do iex=4,4
            g_n=zero
            do d=1,3
              g(d)=lgrad(iex,d)
              g_n=g_n+g(d)*nnt(d)
            end do
            do d=1,3
              g(d)=g(d)-g_n*nnt(d)
              lgrad(iex,d)=g(d)
              rgrad(iex,d)=g(d)
            end do
          end do
        else if (b_code.eq.3)then
          nnt(1)=nx
          nnt(2)=ny
          nnt(3)=nz
          do iex=4,4
            g_n=zero
            do d=1,3
              g(d)=lgrad(iex,d)
              g_n=g_n+g(d)*nnt(d)
            end do
            do d=1,3
              rgrad(iex,d)=lgrad(iex,d)-2.0d0*g_n*nnt(d)
            end do
          end do
          do iv=1,3
            do jj=1,3
              mat_r(iv,jj)=zero
            end do
            mat_r(iv,iv)=1.0d0
          end do
          do iv=1,3
            do jj=1,3
              mat_r(iv,jj)=mat_r(iv,jj)-2.0d0*nnt(iv)*nnt(jj)
            end do
          end do
          mat_a(1,1)=lgrad(1,1); mat_a(1,2)=lgrad(1,2); mat_a(1,3)=lgrad(1,3)
          mat_a(2,1)=lgrad(2,1); mat_a(2,2)=lgrad(2,2); mat_a(2,3)=lgrad(2,3)
          mat_a(3,1)=lgrad(3,1); mat_a(3,2)=lgrad(3,2); mat_a(3,3)=lgrad(3,3)
          do iv=1,3
            do jj=1,3
              mat_c(iv,jj)=zero
              do k=1,3
                mat_c(iv,jj)=mat_c(iv,jj)+mat_r(iv,k)*mat_a(k,jj)
              end do
            end do
          end do
          do iv=1,3
            do jj=1,3
              mat_b(iv,jj)=zero
              do k=1,3
                mat_b(iv,jj)=mat_b(iv,jj)+mat_c(iv,k)*mat_r(k,jj)
              end do
            end do
          end do
          do d=1,3
            rgrad(1,d)=mat_b(1,d)
            rgrad(2,d)=mat_b(2,d)
            rgrad(3,d)=mat_b(3,d)
          end do
        end if
      else
        cright(1:5)=rec_uleft(1:5,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
        do iv=1,3
          do k=1,4
            rgrad(k,iv)=rec_uleftv(iv,k,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
          end do
        end do
        if ((ielem_ibounds(l,i).gt.0).and.(per_rot.eq.1))then
          code_per=ibound_icode(ielem_ibounds(l,i))
          if (code_per.eq.50)then
            cright(2:4)=rotate_per_1(cright(2:4),code_per,angle_per)
            do iv=1,3
              rot_tmp(1:3)=rgrad(1:3,iv)
              rgrad(1:3,iv)=rotate_per_1(rot_tmp(1:3),code_per,angle_per)
            end do
          end if
        end if
      end if
    else
      nfx=ielem_ineighn(l,i)
      lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
      rowfx=bound_offset(nfx)+lfx-1
      cright(1:5)=boundhir(rowfx,1:5)
      ittt=0
      do k=1,4
        do iv=1,3
          ittt=ittt+1
          rgrad(k,iv)=boundhir(rowfx,5+ittt)
        end do
      end do
      if ((ielem_ibounds(l,i).gt.0).and.(per_rot.eq.1))then
        code_per=ibound_icode(ielem_ibounds(l,i))
        if (code_per.eq.50)then
          cright(2:4)=rotate_per_1(cright(2:4),code_per,angle_per)
          do iv=1,3
            rot_tmp(1:3)=rgrad(1:3,iv)
            rgrad(1:3,iv)=rotate_per_1(rot_tmp(1:3),code_per,angle_per)
          end do
        end if
      end if
    end if
    if (lmach.eq.1)then
      cleft_rot(1)=cleft(1)
      cright_rot(1)=cright(1)
      cleft_rot(2)=nx*cleft(2)+ny*cleft(3)+nz*cleft(4)
      cright_rot(2)=nx*cright(2)+ny*cright(3)+nz*cright(4)
      cleft_rot(3)=ca1*ca2*cleft(2)+sa1*ca2*cleft(3)-sa2*cleft(4)
      cright_rot(3)=ca1*ca2*cright(2)+sa1*ca2*cright(3)-sa2*cright(4)
      cleft_rot(4)=-sa1*cleft(2)+ca1*cleft(3)
      cright_rot(4)=-sa1*cright(2)+ca1*cright(3)
      cleft_rot(5)=cleft(5)
      cright_rot(5)=cright(5)
      rho_l=cleft_rot(1)
      rho_r=cright_rot(1)
      u_l=cleft_rot(2)/rho_l
      u_r=cright_rot(2)/rho_r
      v_l=cleft_rot(3)/rho_l
      v_r=cright_rot(3)/rho_r
      w_l=cleft_rot(4)/rho_l
      w_r=cright_rot(4)/rho_r
      eel=cleft_rot(5)
      eer=cright_rot(5)
      ppl=(gamma-1.0d0)*(eel-oo2*rho_l*(u_l*u_l+v_l*v_l+w_l*w_l))
      ppr=(gamma-1.0d0)*(eer-oo2*rho_r*(u_r*u_r+v_r*v_r+w_r*w_r))
      ssl=(gamma*ppl)/rho_l
      ssr=(gamma*ppr)/rho_r
      q2l=u_l*u_l+v_l*v_l+w_l*w_l
      q2r=u_r*u_r+v_r*v_r+w_r*w_r
      mach2=max(q2l/ssl,q2r/ssr)
      mach=sqrt(mach2)
      mach=min(mach,1.0d0)
      dus=u_r+u_l
      dvs=v_r+v_l
      dws=w_r+w_l
      duu=mach*(u_r-u_l)
      dvv=v_r-v_l
      dww=w_r-w_l
      diff=oo2*duu
      u_l=oo2*dus-diff
      u_r=oo2*dus+diff
      if (lmach_style.eq.1)then
        diff=oo2*dvv
        v_l=oo2*dvs-diff
        v_r=oo2*dvs+diff
        diff=oo2*dww
        w_l=oo2*dws-diff
        w_r=oo2*dws+diff
      end if
      cleft_rot(2)=rho_l*u_l
      cright_rot(2)=rho_r*u_r
      cleft_rot(3)=rho_l*v_l
      cright_rot(3)=rho_r*v_r
      cleft_rot(4)=rho_l*w_l
      cright_rot(4)=rho_r*w_r
      cleft(2)=nx*cleft_rot(2)+ca1*ca2*cleft_rot(3)-sa1*cleft_rot(4)
      cright(2)=nx*cright_rot(2)+ca1*ca2*cright_rot(3)-sa1*cright_rot(4)
      cleft(3)=ny*cleft_rot(2)+sa1*ca2*cleft_rot(3)+ca1*cleft_rot(4)
      cright(3)=ny*cright_rot(2)+sa1*ca2*cright_rot(3)+ca1*cright_rot(4)
      cleft(4)=nz*cleft_rot(2)-sa2*cleft_rot(3)
      cright(4)=nz*cright_rot(2)-sa2*cright_rot(3)
    end if
    rho_l=cleft(1)
    rho_r=cright(1)
    u_l=cleft(2)/rho_l
    v_l=cleft(3)/rho_l
    w_l=cleft(4)/rho_l
    u_r=cright(2)/rho_r
    v_r=cright(3)/rho_r
    w_r=cright(4)/rho_r
    t_l=((gamma-1.0d0)*(cleft(5)-oo2*rho_l*(u_l*u_l+v_l*v_l+w_l*w_l)))/(rho_l*r_gas)
    t_r=((gamma-1.0d0)*(cright(5)-oo2*rho_r*(u_r*u_r+v_r*v_r+w_r*w_r)))/(rho_r*r_gas)
    t0=pres/(rres*r_gas)
    mu_l=visc*((t_l/t0)*sqrt(t_l/t0))*((t0+suther*t0)/(t_l+suther*t0))
    mu_r=visc*((t_r/t0)*sqrt(t_r/t0))*((t0+suther*t0)/(t_r+suther*t0))
    lam_l=mu_l*r_gas*gamma/(prandtl*(gamma-1.0d0))
    lam_r=mu_r*r_gas*gamma/(prandtl*(gamma-1.0d0))
    mu_av=oo2*(mu_l+mu_r)
    lam_av=oo2*(lam_l+lam_r)
    leftv(2)=u_l
    leftv(3)=v_l
    leftv(4)=w_l
    leftv(5)=t_l
    rightv(2)=u_r
    rightv(3)=v_r
    rightv(4)=w_r
    rightv(5)=t_r
    vdamp=4.0d0/3.0d0
    do k=1,4
      grad(k,1)=oo2*(lgrad(k,1)+rgrad(k,1))+damp_loc*((vdamp/abs(ielem_dih(l,i)))*nx*(rightv(k+1)-leftv(k+1)))
      grad(k,2)=oo2*(lgrad(k,2)+rgrad(k,2))+damp_loc*((vdamp/abs(ielem_dih(l,i)))*ny*(rightv(k+1)-leftv(k+1)))
      grad(k,3)=oo2*(lgrad(k,3)+rgrad(k,3))+damp_loc*((vdamp/abs(ielem_dih(l,i)))*nz*(rightv(k+1)-leftv(k+1)))
    end do
    q1=-lam_av*grad(4,1)
    q2=-lam_av*grad(4,2)
    q3=-lam_av*grad(4,3)
    ux=grad(1,1)
    uy=grad(1,2)
    uz=grad(1,3)
    vx=grad(2,1)
    vy=grad(2,2)
    vz=grad(2,3)
    wx=grad(3,1)
    wy=grad(3,2)
    wz=grad(3,3)
    tau11=mu_av*((4.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy-(2.0d0/3.0d0)*wz)
    tau22=mu_av*((4.0d0/3.0d0)*vy-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*wz)
    tau33=mu_av*((4.0d0/3.0d0)*wz-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy)
    tau12=mu_av*(uy+vx)
    tau13=mu_av*(wx+uz)
    tau23=mu_av*(vz+wy)
    flux2=nx*tau11+ny*tau12+nz*tau13
    flux3=nx*tau12+ny*tau22+nz*tau23
    flux4=nx*tau13+ny*tau23+nz*tau33
    u12=oo2*(leftv(2)+rightv(2))
    v12=oo2*(leftv(3)+rightv(3))
    w12=oo2*(leftv(4)+rightv(4))
    fx5=-q1+u12*tau11+v12*tau12+w12*tau13
    fy5=-q2+u12*tau12+v12*tau22+w12*tau23
    fz5=-q3+u12*tau13+v12*tau23+w12*tau33
    flux5=nx*fx5+ny*fy5+nz*fz5
    if (ielem_types_faces(l,i).eq.5)then
      wface=weights_q(ngp)*ielem_surf(l,i)
    else
      wface=weights_t(ngp)*ielem_surf(l,i)
    end if
    god2=god2+flux2*wface
    god3=god3+flux3*wface
    god4=god4+flux4*wface
    god5=god5+flux5*wface
  end do
  rhs_val(2,i)=rhs_val(2,i)-god2
  rhs_val(3,i)=rhs_val(3,i)-god3
  rhs_val(4,i)=rhs_val(4,i)-god4
  rhs_val(5,i)=rhs_val(5,i)-god5
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_diffusive_bound_cell_ideal_noturb_3d


#endif
#endif
subroutine calculate_fluxeshi_diffusive_inner_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz,nall
real,dimension(1:4)::viscl,laml
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:gpu_max_nvar)::fxv,fyv,fzv
real,dimension(3,3)::taul,tau
real,dimension(3)::q
real::angle1,angle2,nx,ny,nz,wface,damp,vdamp,mp_pinfl,mp_pinfr,gammal,gammar
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,muturb

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, turbulenceequations, passivescalar, turbulence, dimensiona) &
!$omp& firstprivate(turbulencemodel, icoupleturb, lmach, lmach_style, mrf, per_rot, angle_per, zero, oo2) &
!$omp& firstprivate(gamma, r_gas, pres, rres, visc, suther, prandtl, prtu, lamx, br2_yn, sigma, tolsmall) &
!$omp& firstprivate(cv1, ispal, charlength) &
!$omp& firstprivate(ufreestream, sigma_om2, sigma_om1, sigma_k1, sigma_k2, beta_i1, beta_i2, alpha_inf1) &
!$omp& firstprivate(alpha_inf2, alpha_star0, alpha_starinf, r_k_sst, aa_1, qp_quad_n, qp_triangle_n) &
!$omp& map(alloc: weights_q, weights_t, el_int, ielem_ifca, ielem_types_faces) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_dih, ielem_walldist) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh) &
!$omp& map(alloc: rec_mrf, rec_rotvel, rec_uleft, rec_uleftturb, rec_uleftv, rec_uleftturbv, u_ct_val) &
!$omp& map(alloc: rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,fzv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,muturb)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,fzv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,muturb)
#endif
do ii=1,nof_interior
i=el_int(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
damp=lamx
if (br2_yn.eq.2) damp=zero
do l=1,ielem_ifca(i)
  godflux2=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  nx=cos(angle1)*sin(angle2)
  ny=sin(angle1)*sin(angle2)
  nz=cos(angle2)
  nall(1)=nx
  nall(2)=ny
  nall(3)=nz
  if (ielem_types_faces(l,i).eq.5)then
    iqp=qp_quad_n
  else
    iqp=qp_triangle_n
  end if
  do ngp=1,iqp
    facex=l
    pointx=ngp
    b_code=0
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    call calculate_interior_viscous_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    if (lmach.eq.1)then
      call rotatef(n,cright_rot,cright,angle1,angle2)
      call rotatef(n,cleft_rot,cleft,angle1,angle2)
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
      call rotateb(n,cleft,cleft_rot,angle1,angle2)
      call rotateb(n,cright,cright_rot,angle1,angle2)
    end if
    leftv(1:nof_variables)=cleft(1:nof_variables)
    rightv(1:nof_variables)=cright(1:nof_variables)
    call get_visc_conduct_ideal(n,leftv,rightv,viscl,laml)
    if (turbulence.eq.1)then
      eddyfl=zero
      eddyfr=zero
      turbmv=zero
      etvm=zero
      if (turbulencemodel.eq.1)then
        turbmv(1)=cturbl(1)
        turbmv(2)=cturbr(1)
        eddyfl(2)=turbmv(1)
        eddyfr(2)=turbmv(2)
      else if (turbulencemodel.eq.2)then
        eddyfl(1)=ielem_walldist(i)
        eddyfl(2)=cturbl(1)
        eddyfl(3)=cturbl(2)
        eddyfl(4:6)=lcvgrad(1,1:3)
        eddyfl(7:9)=lcvgrad(2,1:3)
        eddyfl(10:12)=lcvgrad(3,1:3)
        eddyfl(13:15)=lcvgrad_t(1,1:3)
        eddyfl(16:18)=lcvgrad_t(2,1:3)
        eddyfr(1)=ielem_walldist(i)
        eddyfr(2)=cturbr(1)
        eddyfr(3)=cturbr(2)
        eddyfr(4:6)=rcvgrad(1,1:3)
        eddyfr(7:9)=rcvgrad(2,1:3)
        eddyfr(10:12)=rcvgrad(3,1:3)
        eddyfr(13:15)=rcvgrad_t(1,1:3)
        eddyfr(16:18)=rcvgrad_t(2,1:3)
      end if
      call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
    end if
    fxv=zero
    fyv=zero
    fzv=zero
    taul=zero
    tau=zero
    q=zero
    call cons2div_ideal(n,leftv,mp_pinfl,gammal)
    call cons2div_ideal(n,rightv,mp_pinfr,gammar)
    u12=oo2*(leftv(2)+rightv(2))
    v12=oo2*(leftv(3)+rightv(3))
    w12=oo2*(leftv(4)+rightv(4))
    vdamp=4.0d0/3.0d0
    do k=1,nof_variables-1
      do iv=1,3
        lcvgrad(k,iv)=oo2*(lcvgrad(k,iv)+rcvgrad(k,iv))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(iv)*(rightv(k+1)-leftv(k+1)))
      end do
    end do
    if (turbulence.eq.1)then
      q(1:3)=-oo2*(laml(3)+laml(4))*lcvgrad(dimensiona+1,1:3)
    else
      q(1:3)=-oo2*(laml(1)+laml(2))*lcvgrad(dimensiona+1,1:3)
    end if
    fxv(5)=fxv(5)-q(1)
    fyv(5)=fyv(5)-q(2)
    fzv(5)=fzv(5)-q(3)
    ux=lcvgrad(1,1)
    uy=lcvgrad(1,2)
    uz=lcvgrad(1,3)
    vx=lcvgrad(2,1)
    vy=lcvgrad(2,2)
    vz=lcvgrad(2,3)
    wx=lcvgrad(3,1)
    wy=lcvgrad(3,2)
    wz=lcvgrad(3,3)
    taul(1,1)=(4.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy-(2.0d0/3.0d0)*wz
    taul(2,2)=(4.0d0/3.0d0)*vy-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*wz
    taul(3,3)=(4.0d0/3.0d0)*wz-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy
    taul(1,2)=uy+vx
    taul(2,1)=taul(1,2)
    taul(1,3)=wx+uz
    taul(3,1)=taul(1,3)
    taul(2,3)=vz+wy
    taul(3,2)=taul(2,3)
    if (turbulence.eq.1)then
      tau=oo2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))*taul
    else
      tau=oo2*(viscl(1)+viscl(2))*taul
    end if
    do kc=2,4
      fxv(kc)=fxv(kc)+tau(1,kc-1)
      fyv(kc)=fyv(kc)+tau(2,kc-1)
      fzv(kc)=fzv(kc)+tau(3,kc-1)
    end do
    fxv(5)=fxv(5)+u12*tau(1,1)+v12*tau(1,2)+w12*tau(1,3)
    fyv(5)=fyv(5)+u12*tau(2,1)+v12*tau(2,2)+w12*tau(2,3)
    fzv(5)=fzv(5)+u12*tau(3,1)+v12*tau(3,2)+w12*tau(3,3)
    do iv=1,nof_variables
      hllcflux(iv)=nx*fxv(iv)+ny*fyv(iv)+nz*fzv(iv)
    end do
    if (ielem_types_faces(l,i).eq.5)then
      wface=weights_q(ngp)*ielem_surf(l,i)
    else
      wface=weights_t(ngp)*ielem_surf(l,i)
    end if
    do iv=1,nof_variables
      godflux2(iv)=godflux2(iv)+hllcflux(iv)*wface
    end do
    if (nvt.gt.0)then
      if (turbulence.eq.1)then
        muturb=oo2*(viscl(1)+viscl(2))+oo2*(viscl(3)+viscl(4))
      else
        muturb=oo2*(viscl(1)+viscl(2))
      end if
      do iv=1,nvt
        hllcflux(nof_variables+iv)=muturb*oo2*((lcvgrad_t(iv,1)+rcvgrad_t(iv,1))*nx+(lcvgrad_t(iv,2)+rcvgrad_t(iv,2))*ny+ &
             & (lcvgrad_t(iv,3)+rcvgrad_t(iv,3))*nz)
        if (turbulencemodel.eq.1) hllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)/sigma
        godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+hllcflux(nof_variables+iv)*wface
      end do
    end if
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)-godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_diffusive_inner_cell_ideal


subroutine calculate_fluxeshi_diffusive_bound_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz,nall
real,dimension(1:4)::viscl,laml
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:gpu_max_nvar)::fxv,fyv,fzv
real,dimension(3,3)::taul,tau
real,dimension(3)::q
real::angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,mp_pinfl,mp_pinfr,gammal,gammar
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,muturb

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, turbulenceequations, passivescalar, turbulence, dimensiona) &
!$omp& firstprivate(turbulencemodel, icoupleturb, lmach, lmach_style, mrf, per_rot, angle_per, thermal, zero, oo2) &
!$omp& firstprivate(gamma, r_gas, pres, rres, visc, suther, prandtl, prtu, lamx, br2_yn, sigma, tolsmall) &
!$omp& firstprivate(boundtype, initcond, Mach_in, uvel, vvel, wvel, press_outlet, swirl, cv1, ispal, charlength) &
!$omp& firstprivate(ufreestream, sigma_om2, sigma_om1, sigma_k1, sigma_k2, beta_i1, beta_i2, alpha_inf1) &
!$omp& firstprivate(alpha_inf2, alpha_star0, alpha_starinf, r_k_sst, aa_1, qp_quad_n, qp_triangle_n) &
!$omp& map(alloc: weights_q, weights_t, el_bnd, ielem_ifca, ielem_types_faces) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_dih, ielem_walldist) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode, ielem_qface) &
!$omp& map(alloc: ielem_inter_id, ielem_indexf, bound_offset, boundhir, ielem_nodes_faces, inoder4_cord) &
!$omp& map(alloc: rec_mrf, rec_rotvel, rec_uleft, rec_uleftturb, rec_uleftv, rec_uleftturbv, u_ct_val) &
!$omp& map(alloc: rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,fzv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,muturb)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,fzv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,uz,vx,vy,vz,wx,wy,wz,u12,v12,w12,muturb)
#endif
do ii=1,nof_bounded
i=el_bnd(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
do l=1,ielem_ifca(i)
  godflux2=zero
  damp=lamx
  if (br2_yn.eq.2) damp=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  nx=cos(angle1)*sin(angle2)
  ny=sin(angle1)*sin(angle2)
  nz=cos(angle2)
  nall(1)=nx
  nall(2)=ny
  nall(3)=nz
  if (ielem_types_faces(l,i).eq.5)then
    iqp=qp_quad_n
  else
    iqp=qp_triangle_n
  end if
  do ngp=1,iqp
    facex=l
    pointx=ngp
    b_code=0
    damp_loc=damp
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    call calculate_bounded_viscous_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    if (b_code.gt.0) damp_loc=zero
    if (lmach.eq.1)then
      call rotatef(n,cright_rot,cright,angle1,angle2)
      call rotatef(n,cleft_rot,cleft,angle1,angle2)
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
      call rotateb(n,cleft,cleft_rot,angle1,angle2)
      call rotateb(n,cright,cright_rot,angle1,angle2)
    end if
    leftv(1:nof_variables)=cleft(1:nof_variables)
    rightv(1:nof_variables)=cright(1:nof_variables)
    call get_visc_conduct_ideal(n,leftv,rightv,viscl,laml)
    if (turbulence.eq.1)then
      eddyfl=zero
      eddyfr=zero
      turbmv=zero
      etvm=zero
      if (turbulencemodel.eq.1)then
        turbmv(1)=cturbl(1)
        turbmv(2)=cturbr(1)
        eddyfl(2)=turbmv(1)
        eddyfr(2)=turbmv(2)
      else if (turbulencemodel.eq.2)then
        eddyfl(1)=ielem_walldist(i)
        eddyfl(2)=cturbl(1)
        eddyfl(3)=cturbl(2)
        eddyfl(4:6)=lcvgrad(1,1:3)
        eddyfl(7:9)=lcvgrad(2,1:3)
        eddyfl(10:12)=lcvgrad(3,1:3)
        eddyfl(13:15)=lcvgrad_t(1,1:3)
        eddyfl(16:18)=lcvgrad_t(2,1:3)
        eddyfr(1)=ielem_walldist(i)
        eddyfr(2)=cturbr(1)
        eddyfr(3)=cturbr(2)
        eddyfr(4:6)=rcvgrad(1,1:3)
        eddyfr(7:9)=rcvgrad(2,1:3)
        eddyfr(10:12)=rcvgrad(3,1:3)
        eddyfr(13:15)=rcvgrad_t(1,1:3)
        eddyfr(16:18)=rcvgrad_t(2,1:3)
      end if
      call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
    end if
    fxv=zero
    fyv=zero
    fzv=zero
    taul=zero
    tau=zero
    q=zero
    call cons2div_ideal(n,leftv,mp_pinfl,gammal)
    call cons2div_ideal(n,rightv,mp_pinfr,gammar)
    u12=oo2*(leftv(2)+rightv(2))
    v12=oo2*(leftv(3)+rightv(3))
    w12=oo2*(leftv(4)+rightv(4))
    vdamp=4.0d0/3.0d0
    do k=1,nof_variables-1
      do iv=1,3
        lcvgrad(k,iv)=oo2*(lcvgrad(k,iv)+rcvgrad(k,iv))+damp_loc*((vdamp/abs(ielem_dih(l,i)))*nall(iv)*(rightv(k+1)-leftv(k+1)))
      end do
    end do
    if (turbulence.eq.1)then
      q(1:3)=-oo2*(laml(3)+laml(4))*lcvgrad(dimensiona+1,1:3)
    else
      q(1:3)=-oo2*(laml(1)+laml(2))*lcvgrad(dimensiona+1,1:3)
    end if
    fxv(5)=fxv(5)-q(1)
    fyv(5)=fyv(5)-q(2)
    fzv(5)=fzv(5)-q(3)
    ux=lcvgrad(1,1)
    uy=lcvgrad(1,2)
    uz=lcvgrad(1,3)
    vx=lcvgrad(2,1)
    vy=lcvgrad(2,2)
    vz=lcvgrad(2,3)
    wx=lcvgrad(3,1)
    wy=lcvgrad(3,2)
    wz=lcvgrad(3,3)
    taul(1,1)=(4.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy-(2.0d0/3.0d0)*wz
    taul(2,2)=(4.0d0/3.0d0)*vy-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*wz
    taul(3,3)=(4.0d0/3.0d0)*wz-(2.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy
    taul(1,2)=uy+vx
    taul(2,1)=taul(1,2)
    taul(1,3)=wx+uz
    taul(3,1)=taul(1,3)
    taul(2,3)=vz+wy
    taul(3,2)=taul(2,3)
    if (turbulence.eq.1)then
      tau=oo2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))*taul
    else
      tau=oo2*(viscl(1)+viscl(2))*taul
    end if
    do kc=2,4
      fxv(kc)=fxv(kc)+tau(1,kc-1)
      fyv(kc)=fyv(kc)+tau(2,kc-1)
      fzv(kc)=fzv(kc)+tau(3,kc-1)
    end do
    fxv(5)=fxv(5)+u12*tau(1,1)+v12*tau(1,2)+w12*tau(1,3)
    fyv(5)=fyv(5)+u12*tau(2,1)+v12*tau(2,2)+w12*tau(2,3)
    fzv(5)=fzv(5)+u12*tau(3,1)+v12*tau(3,2)+w12*tau(3,3)
    do iv=1,nof_variables
      hllcflux(iv)=nx*fxv(iv)+ny*fyv(iv)+nz*fzv(iv)
    end do
    if (ielem_types_faces(l,i).eq.5)then
      wface=weights_q(ngp)*ielem_surf(l,i)
    else
      wface=weights_t(ngp)*ielem_surf(l,i)
    end if
    do iv=1,nof_variables
      godflux2(iv)=godflux2(iv)+hllcflux(iv)*wface
    end do
    if (nvt.gt.0)then
      if (turbulence.eq.1)then
        muturb=oo2*(viscl(1)+viscl(2))+oo2*(viscl(3)+viscl(4))
      else
        muturb=oo2*(viscl(1)+viscl(2))
      end if
      do iv=1,nvt
        hllcflux(nof_variables+iv)=muturb*oo2*((lcvgrad_t(iv,1)+rcvgrad_t(iv,1))*nx+(lcvgrad_t(iv,2)+rcvgrad_t(iv,2))*ny+ &
             & (lcvgrad_t(iv,3)+rcvgrad_t(iv,3))*nz)
        if (turbulencemodel.eq.1) hllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)/sigma
        if (b_code.eq.3) hllcflux(nof_variables+iv)=zero
        godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+hllcflux(nof_variables+iv)*wface
      end do
    end if
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)-godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_diffusive_bound_cell_ideal


subroutine calculate_fluxeshi_diffusive2d(n)
!> @brief
!> this subroutine computes the diffusive fluxes for euler-navier-stokes equations in 2d
		implicit none
		integer,intent(in)::n
		integer::ii


	if ((dg.eq.0).and.(multispecies.eq.0).and.(realgas.eq.0))then
#ifdef xpu
!$omp barrier
!$omp master
#endif
	call calculate_fluxeshi_diffusive2d_inner_cell_ideal(n)
	call calculate_fluxeshi_diffusive2d_bound_cell_ideal(n)
#ifdef xpu
!$omp end master
!$omp barrier
#endif
	return
	end if







#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	call calculate_fluxeshi_diffusive2d_inner_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_diffusive2d_bound_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine calculate_fluxeshi_diffusive2d

subroutine calculate_fluxeshi_diffusive2d_inner_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,nvar,kc,iex,ittt,ikas,igoflux,icaseb,kk,b_code
	real::sum_detect,norms
	integer::iconsidered,facex,pointx,k,rg_i,rg_j,idxy,icompute
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot,tempx_l,rtempx_l
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz,rg_sumfl,rg_sumfr
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real,dimension(1:8,1:gpu_max_dim)::vext
	real::mp_pinfl,mp_pinfr,gammal,gammar,rho12l,rho12r
	real,dimension(1:4)::viscl,laml
	real,dimension(1:gpu_max_species)::rgs_htr_i, rgs_hvib_i,y_av,x_av
	real,dimension(1:2)::turbmv
    real,dimension(1)::etvm
    real,dimension(1:20)::eddyfl,eddyfr
    real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad,jtmp
	real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
	real,dimension(1:gpu_max_nvar)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(2,2)::taul,taur,tau
	real,dimension(2)::q,nall,qvib
	real::ux,uy,uz,vx,vy,vz,wx,wy,wz,rho12,u12,v12,w12,damp,vdamp,y_face
	real::mp_ttr,mp_tv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_laml,mp_lamr
	real,dimension(1:gpu_max_species)::mp_d_eff,rg_difl,rg_difr,rg_enthl,rg_enthr,mp_htr,mp_hvib,rg_enthvbl,rg_enthvbr,mp_mu_i,mp_ktr_i
	real,dimension(1:gpu_max_dim)::rg_sum_tr,rg_sumfr_tr,rg_sumfl_tr,rg_sumfl_v,rg_sumfr_v,rg_sum_v,gradyl,gradyr,jl,jr
	real,dimension(1:2)::qtr,qv
	real,dimension(1:gpu_max_dim)::rg_sum_htr,rg_sum_hv,grady,rg_qv,rg_qtr
	real,dimension(1:gpu_max_species)::rg_dif_av,rg_enth_av,rg_enthvb_av
	real::mp_lam_av,mp_ktr_mix_av,sum_y1,sum_y2
	real,dimension(1:gpu_max_dim):: sumi,gradx,grad_a
	real,dimension(1:gpu_max_dim,1:gpu_max_species) :: i_raw
	 real :: sumY_face
	real :: y_corr(gpu_max_species)
	real :: sumJ(1:gpu_max_dim)
	real,dimension(1:gpu_max_dim)::sumGradY
	real,dimension(1:gpu_max_dim,1:gpu_max_species)::gradY_all
	real :: molar_sum, mbar


	i=el_int(ii)
	iconsidered=i

	weights_temp(1:qp_line_n)=weights_l(1:qp_line_n)

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
							  call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if
						      if (turbulencemodel.eq.2)then
							  eddyfl(1)=ielem_walldist(i);eddyfl(2)=cturbl(1);eddyfl(3)=cturbl(2)
							  eddyfl(4:5)= lcvgrad(1,1:2);eddyfl(6:7)=lcvgrad(2,1:2)
							  eddyfl(8:9)=lcvgrad_t(1,1:2)
							  eddyfl(10:11)=lcvgrad_t(2,1:2)


							  eddyfr(1)=ielem_walldist(i);eddyfr(2)=cturbr(1);eddyfr(3)=cturbr(2)
							  eddyfr(4:5)= rcvgrad(1,1:2);eddyfr(6:7)=rcvgrad(2,1:2)
							  eddyfr(8:9)=rcvgrad_t(1,1:2);eddyfr(10:11)=rcvgrad_t(2,1:2)
							    call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
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

													sumi(1:dimensiona) = 0.0d0

													if (rg_relax.ge.1)then


													sumi(1:dimensiona) = 0.0d0
													sumJ(1:dimensiona) = 0.0d0
													i_raw(1:dimensiona,1:nof_species) = 0.0d0

													! -------------------------------------------------------------
													! Normalise face mass fractions used in diffusion correction.
													! This is essential because:
													! sum(J_k) = sum(i_raw_k) - sum(Y_k)*sum(i_raw_k)
													! so conservation requires sum(Y_k)=1 at the face.
													! -------------------------------------------------------------
													sumY_face = 0.0d0

													do rg_j = 1, nof_species
														y_corr(rg_j) = max(y_av(rg_j), 0.0d0)
														sumY_face = sumY_face + y_corr(rg_j)
													end do

													if (sumY_face .gt. 1.0d-30) then
														do rg_j = 1, nof_species
															y_corr(rg_j) = y_corr(rg_j) / sumY_face
														end do
													else
														do rg_j = 1, nof_species
															y_corr(rg_j) = 1.0d0 / dble(nof_species)
														end do
													end if

													! -------------------------------------------------------------
													! Compute face mole fractions from corrected face mass fractions
													! -------------------------------------------------------------
													molar_sum = 0.0d0
													do rg_j = 1, nof_species
														molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
													end do

													molar_sum = max(molar_sum,1.0d-30)
													mbar = 1.0d0 / molar_sum

													do rg_j = 1, nof_species
														x_av(rg_j) = (y_corr(rg_j) / rg_molm(rg_j)) / molar_sum
													end do

													! -------------------------------------------------------------
													! Project species gradients so that sum_k grad(Y_k) = 0.
													! This removes a conservative but noisy pairwise diffusion mode
													! which can enter rhoEv through sum_k(hvib_k*J_k).
													! -------------------------------------------------------------
													sumGradY(1:dimensiona) = 0.0d0

													do rg_j = 1, nof_species
														idxy = dimensiona + 2 + rg_j
														gradY_all(1:dimensiona,rg_j) = lcvgrad(idxy,1:dimensiona)
														sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona,rg_j)
													end do

													do rg_j = 1, nof_species
														gradY_all(1:dimensiona,rg_j) = gradY_all(1:dimensiona,rg_j) - &
															y_corr(rg_j) * sumGradY(1:dimensiona)
													end do

													! -------------------------------------------------------------
													! grad_a = grad(sum_k Y_k/M_k), using projected grad(Y_k).
													! -------------------------------------------------------------
													grad_a(1:dimensiona) = 0.0d0

													do rg_j = 1, nof_species
														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_j)

														grad_a(1:dimensiona) = grad_a(1:dimensiona) + &
															grady(1:dimensiona) / rg_molm(rg_j)
													end do

													! -------------------------------------------------------------
													! Raw mixture-averaged flux:
													! i_raw_i = -rho * D_i * (M_i/Mmix) * grad(X_i)
													! -------------------------------------------------------------
													do rg_i = 1, nof_species

														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_i)

														gradx(1:dimensiona) = &
															( grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
															- (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona) ) &
															/ (molar_sum * molar_sum)

														i_raw(1:dimensiona,rg_i) = &
															-rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)

														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)

													end do

													! -------------------------------------------------------------
													! Mass-conserving corrected species diffusion flux:
													! J_i = i_raw_i - Y_i * sum_k(i_raw_k)
													! using corrected Y_i with sum(Y_i)=1.
													! -------------------------------------------------------------
													do rg_i = 1, nof_species

														y_face = y_corr(rg_i)

														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - &
																		  y_face * sumi(1:dimensiona)

														sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)

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

					call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,-1.0d0)

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

	end subroutine calculate_fluxeshi_diffusive2d_inner_cell


subroutine calculate_fluxeshi_diffusive2d_bound_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2,dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,nvar,kc,iex,ittt,ikas,igoflux,icaseb,kk,b_code
	real::sum_detect,norms
	integer::iconsidered,facex,pointx,k,rg_i,rg_j,idxy,icompute
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot,tempx_l,rtempx_l
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz,rg_sumfl,rg_sumfr
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:gpu_max_qp_face)::weights_temp
	real,dimension(1:8,1:gpu_max_dim)::vext
	real::mp_pinfl,mp_pinfr,gammal,gammar,rho12l,rho12r
	real,dimension(1:4)::viscl,laml
	real,dimension(1:gpu_max_species)::rgs_htr_i, rgs_hvib_i,y_av,x_av
	real,dimension(1:2)::turbmv
    real,dimension(1)::etvm
    real,dimension(1:20)::eddyfl,eddyfr
    real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad,jtmp
	real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
	real,dimension(1:gpu_max_nvar)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(2,2)::taul,taur,tau
	real,dimension(2)::q,nall,qvib
	real::ux,uy,uz,vx,vy,vz,wx,wy,wz,rho12,u12,v12,w12,damp,vdamp,y_face
	real::mp_ttr,mp_tv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_laml,mp_lamr
	real,dimension(1:gpu_max_species)::mp_d_eff,rg_difl,rg_difr,rg_enthl,rg_enthr,mp_htr,mp_hvib,rg_enthvbl,rg_enthvbr,mp_mu_i,mp_ktr_i
	real,dimension(1:gpu_max_dim)::rg_sum_tr,rg_sumfr_tr,rg_sumfl_tr,rg_sumfl_v,rg_sumfr_v,rg_sum_v,gradyl,gradyr,jl,jr
	real,dimension(1:2)::qtr,qv
	real,dimension(1:gpu_max_dim)::rg_sum_htr,rg_sum_hv,grady,rg_qv,rg_qtr
	real,dimension(1:gpu_max_species)::rg_dif_av,rg_enth_av,rg_enthvb_av
	real::mp_lam_av,mp_ktr_mix_av,sum_y1,sum_y2
	real,dimension(1:gpu_max_dim):: sumi,gradx,grad_a
	real,dimension(1:gpu_max_dim,1:gpu_max_species) :: i_raw
	 real :: sumY_face
	real :: y_corr(gpu_max_species)
	real :: sumJ(1:gpu_max_dim)
	real,dimension(1:gpu_max_dim)::sumGradY
	real,dimension(1:gpu_max_dim,1:gpu_max_species)::gradY_all
	real :: molar_sum, mbar


	i=el_bnd(ii)
	iconsidered=i

			weights_temp(1:qp_line_n)=weights_l(1:qp_line_n)





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
							  call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
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
							    call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
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

													sumi(1:dimensiona) = 0.0d0

													if (rg_relax.ge.1)then
															if ((b_code.eq.4).and.(catalytic_wall.eq.0))then
																icompute=1
															end if

													if (icompute.eq.0)then

													sumi(1:dimensiona) = 0.0d0
													sumJ(1:dimensiona) = 0.0d0
													i_raw(1:dimensiona,1:nof_species) = 0.0d0

													! -------------------------------------------------------------
													! Normalise face mass fractions used in diffusion correction.
													! This is essential because:
													! sum(J_k) = sum(i_raw_k) - sum(Y_k)*sum(i_raw_k)
													! so conservation requires sum(Y_k)=1 at the face.
													! -------------------------------------------------------------
													sumY_face = 0.0d0

													do rg_j = 1, nof_species
														y_corr(rg_j) = max(y_av(rg_j), 0.0d0)
														sumY_face = sumY_face + y_corr(rg_j)
													end do

													if (sumY_face .gt. 1.0d-30) then
														do rg_j = 1, nof_species
															y_corr(rg_j) = y_corr(rg_j) / sumY_face
														end do
													else
														do rg_j = 1, nof_species
															y_corr(rg_j) = 1.0d0 / dble(nof_species)
														end do
													end if

													! -------------------------------------------------------------
													! Compute face mole fractions from corrected face mass fractions
													! -------------------------------------------------------------
													molar_sum = 0.0d0
													do rg_j = 1, nof_species
														molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
													end do

													molar_sum = max(molar_sum,1.0d-30)
													mbar = 1.0d0 / molar_sum

													do rg_j = 1, nof_species
														x_av(rg_j) = (y_corr(rg_j) / rg_molm(rg_j)) / molar_sum
													end do

													! -------------------------------------------------------------
													! Project species gradients so that sum_k grad(Y_k) = 0.
													! This removes a conservative but noisy pairwise diffusion mode
													! which can enter rhoEv through sum_k(hvib_k*J_k).
													! -------------------------------------------------------------
													sumGradY(1:dimensiona) = 0.0d0

													do rg_j = 1, nof_species
														idxy = dimensiona + 2 + rg_j
														gradY_all(1:dimensiona,rg_j) = lcvgrad(idxy,1:dimensiona)
														sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona,rg_j)
													end do

													do rg_j = 1, nof_species
														gradY_all(1:dimensiona,rg_j) = gradY_all(1:dimensiona,rg_j) - &
															y_corr(rg_j) * sumGradY(1:dimensiona)
													end do

													! -------------------------------------------------------------
													! grad_a = grad(sum_k Y_k/M_k), using projected grad(Y_k).
													! -------------------------------------------------------------
													grad_a(1:dimensiona) = 0.0d0

													do rg_j = 1, nof_species
														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_j)

														grad_a(1:dimensiona) = grad_a(1:dimensiona) + &
															grady(1:dimensiona) / rg_molm(rg_j)
													end do

													! -------------------------------------------------------------
													! Raw mixture-averaged flux:
													! i_raw_i = -rho * D_i * (M_i/Mmix) * grad(X_i)
													! -------------------------------------------------------------
													do rg_i = 1, nof_species

														grady(1:dimensiona) = gradY_all(1:dimensiona,rg_i)

														gradx(1:dimensiona) = &
															( grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
															- (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona) ) &
															/ (molar_sum * molar_sum)

														i_raw(1:dimensiona,rg_i) = &
															-rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)

														sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona,rg_i)

													end do

													! -------------------------------------------------------------
													! Mass-conserving corrected species diffusion flux:
													! J_i = i_raw_i - Y_i * sum_k(i_raw_k)
													! using corrected Y_i with sum(Y_i)=1.
													! -------------------------------------------------------------
													do rg_i = 1, nof_species

														y_face = y_corr(rg_i)

														jl(1:dimensiona) = i_raw(1:dimensiona,rg_i) - &
																		  y_face * sumi(1:dimensiona)

														sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

														fxv(dimensiona+3+rg_i) = fxv(dimensiona+3+rg_i) + jl(1)
														fyv(dimensiona+3+rg_i) = fyv(dimensiona+3+rg_i) + jl(2)

														rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
																				+ rg_enth_av(rg_i)*jl(1:dimensiona)

														rg_sum_hv(1:dimensiona)  = rg_sum_hv(1:dimensiona) &
																				+ rg_enthvb_av(rg_i)*jl(1:dimensiona)

													end do

													! Optional diagnostic


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

				call dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,-1.0d0)






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

	end subroutine calculate_fluxeshi_diffusive2d_bound_cell





subroutine calculate_fluxeshi_diffusive2d_inner_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz,nall
real,dimension(1:4)::viscl,laml
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:gpu_max_nvar)::fxv,fyv
real,dimension(2,2)::taul,tau
real,dimension(2)::q
real::angle1,angle2,nx,ny,nz,wface,damp,vdamp,mp_pinfl,mp_pinfr,gammal,gammar
real::ux,uy,vx,vy,u12,v12,muturb

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, turbulenceequations, passivescalar, turbulence, dimensiona) &
!$omp& firstprivate(turbulencemodel, icoupleturb, lmach, lmach_style, per_rot, angle_per, zero, oo2, gamma) &
!$omp& firstprivate(r_gas, pres, rres, visc, suther, prandtl, prtu, lamx, br2_yn, sigma, tolsmall) &
!$omp& firstprivate(cv1, ispal, charlength, ufreestream, sigma_om2, sigma_om1, sigma_k1) &
!$omp& firstprivate(sigma_k2, beta_i1, beta_i2, alpha_inf1, alpha_inf2, alpha_star0, alpha_starinf, r_k_sst) &
!$omp& firstprivate(aa_1, qp_line_n) &
!$omp& map(alloc: weights_l, el_int, ielem_ifca, ielem_faceanglex, ielem_faceangley) &
!$omp& map(alloc: ielem_surf, ielem_dih, ielem_walldist) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh) &
!$omp& map(alloc: rec_uleft, rec_uleftturb, rec_uleftv, rec_uleftturbv, u_ct_val, rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,vx,vy,u12,v12,muturb)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,vx,vy,u12,v12,muturb)
#endif
do ii=1,nof_interior
i=el_int(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
damp=lamx
if (br2_yn.eq.2) damp=zero
do l=1,ielem_ifca(i)
  godflux2=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  nx=angle1
  ny=angle2
  nz=zero
  nall(1)=nx
  nall(2)=ny
  iqp=qp_line_n
  do ngp=1,iqp
    facex=l
    pointx=ngp
    b_code=0
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    call calculate_interior_viscous2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    if (lmach.eq.1)then
      call rotatef2d(n,cright_rot,cright,angle1,angle2)
      call rotatef2d(n,cleft_rot,cleft,angle1,angle2)
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht2d_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
      call rotateb2d(n,cright,cright_rot,angle1,angle2)
    end if
    leftv(1:nof_variables)=cleft(1:nof_variables)
    rightv(1:nof_variables)=cright(1:nof_variables)
    call get_visc_conduct_ideal(n,leftv,rightv,viscl,laml)
    if (turbulence.eq.1)then
      eddyfl=zero
      eddyfr=zero
      turbmv=zero
      etvm=zero
      if (turbulencemodel.eq.1)then
        turbmv(1)=cturbl(1)
        turbmv(2)=cturbr(1)
        eddyfl(2)=turbmv(1)
        eddyfr(2)=turbmv(2)
      else if (turbulencemodel.eq.2)then
        eddyfl(1)=ielem_walldist(i)
        eddyfl(2)=cturbl(1)
        eddyfl(3)=cturbl(2)
        eddyfl(4:5)=lcvgrad(1,1:2)
        eddyfl(6:7)=lcvgrad(2,1:2)
        eddyfl(8:9)=lcvgrad_t(1,1:2)
        eddyfl(10:11)=lcvgrad_t(2,1:2)
        eddyfr(1)=ielem_walldist(i)
        eddyfr(2)=cturbr(1)
        eddyfr(3)=cturbr(2)
        eddyfr(4:5)=rcvgrad(1,1:2)
        eddyfr(6:7)=rcvgrad(2,1:2)
        eddyfr(8:9)=rcvgrad_t(1,1:2)
        eddyfr(10:11)=rcvgrad_t(2,1:2)
      end if
      call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
    end if
    fxv=zero
    fyv=zero
    taul=zero
    tau=zero
    q=zero
    call cons2div_ideal(n,leftv,mp_pinfl,gammal)
    call cons2div_ideal(n,rightv,mp_pinfr,gammar)
    u12=oo2*(leftv(2)+rightv(2))
    v12=oo2*(leftv(3)+rightv(3))
    vdamp=4.0d0/3.0d0
    do k=1,nof_variables-1
      lcvgrad(k,1:2)=oo2*(lcvgrad(k,1:2)+rcvgrad(k,1:2))+damp*((vdamp/abs(ielem_dih(l,i)))*nall(1:2)*(rightv(k+1)-leftv(k+1)))
    end do
    if (turbulence.eq.1)then
      q(1:2)=-oo2*(laml(3)+laml(4))*lcvgrad(dimensiona+1,1:2)
    else
      q(1:2)=-oo2*(laml(1)+laml(2))*lcvgrad(dimensiona+1,1:2)
    end if
    fxv(4)=fxv(4)-q(1)
    fyv(4)=fyv(4)-q(2)
    ux=lcvgrad(1,1)
    uy=lcvgrad(1,2)
    vx=lcvgrad(2,1)
    vy=lcvgrad(2,2)
    taul(1,1)=(4.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy
    taul(2,2)=(4.0d0/3.0d0)*vy-(2.0d0/3.0d0)*ux
    taul(1,2)=uy+vx
    taul(2,1)=taul(1,2)
    if (turbulence.eq.1)then
      tau=oo2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))*taul
    else
      tau=oo2*(viscl(1)+viscl(2))*taul
    end if
    do kc=2,3
      fxv(kc)=fxv(kc)+tau(1,kc-1)
      fyv(kc)=fyv(kc)+tau(2,kc-1)
    end do
    fxv(4)=fxv(4)+u12*tau(1,1)+v12*tau(1,2)
    fyv(4)=fyv(4)+u12*tau(2,1)+v12*tau(2,2)
    do iv=1,nof_variables
      hllcflux(iv)=nx*fxv(iv)+ny*fyv(iv)
    end do
    wface=weights_l(ngp)*ielem_surf(l,i)
    do iv=1,nof_variables
      godflux2(iv)=godflux2(iv)+hllcflux(iv)*wface
    end do
    if (nvt.gt.0)then
      if (turbulence.eq.1)then
        muturb=oo2*(viscl(1)+viscl(2))+oo2*(viscl(3)+viscl(4))
      else
        muturb=oo2*(viscl(1)+viscl(2))
      end if
      do iv=1,nvt
        hllcflux(nof_variables+iv)=muturb*oo2*((lcvgrad_t(iv,1)+rcvgrad_t(iv,1))*nx+(lcvgrad_t(iv,2)+rcvgrad_t(iv,2))*ny)
        if (turbulencemodel.eq.1) hllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)/sigma
        godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+hllcflux(nof_variables+iv)*wface
      end do
    end if
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)-godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_diffusive2d_inner_cell_ideal


subroutine calculate_fluxeshi_diffusive2d_bound_cell_ideal(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code
real,dimension(1:gpu_max_nvar_total)::godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz,nall
real,dimension(1:4)::viscl,laml
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:gpu_max_nvar)::fxv,fyv
real,dimension(2,2)::taul,tau
real,dimension(2)::q
real::angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,mp_pinfl,mp_pinfr,gammal,gammar
real::ux,uy,vx,vy,u12,v12,muturb

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, turbulenceequations, passivescalar, turbulence, dimensiona) &
!$omp& firstprivate(turbulencemodel, icoupleturb, lmach, lmach_style, per_rot, angle_per, thermal, zero, oo2, gamma) &
!$omp& firstprivate(r_gas, pres, rres, visc, suther, prandtl, prtu, lamx, br2_yn, sigma, tolsmall, boundtype) &
!$omp& firstprivate(initcond, Mach_in, uvel, vvel, cv1, ispal, charlength, ufreestream, sigma_om2, sigma_om1, sigma_k1) &
!$omp& firstprivate(sigma_k2, beta_i1, beta_i2, alpha_inf1, alpha_inf2, alpha_star0, alpha_starinf, r_k_sst) &
!$omp& firstprivate(aa_1, qp_line_n) &
!$omp& map(alloc: weights_l, el_bnd, ielem_ifca, ielem_faceanglex, ielem_faceangley) &
!$omp& map(alloc: ielem_surf, ielem_dih, ielem_walldist) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode, ielem_qface) &
!$omp& map(alloc: ielem_inter_id, ielem_indexf, bound_offset, boundhir, ielem_nodes_faces, inoder4_cord) &
!$omp& map(alloc: rec_uleft, rec_uleftturb, rec_uleftv, rec_uleftturbv, u_ct_val, rhs_val, rhst_val) &
!$omp& private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,vx,vy,u12,v12,muturb)
#else
!$omp do private(ii,i,l,ngp,iqp,iv,k,kc,nvt,iconsidered,facex,pointx,b_code) &
!$omp& private(godflux2,hllcflux,cleft,cright,cleft_rot,cright_rot,leftv,rightv,srf_speedrot) &
!$omp& private(cturbl,cturbr,pox,poy,poz,nall,viscl,laml,turbmv,etvm,eddyfl,eddyfr) &
!$omp& private(lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t,fxv,fyv,taul,tau,q) &
!$omp& private(angle1,angle2,nx,ny,nz,wface,damp,damp_loc,vdamp,mp_pinfl,mp_pinfr,gammal,gammar) &
!$omp& private(ux,uy,vx,vy,u12,v12,muturb)
#endif
do ii=1,nof_bounded
i=el_bnd(ii)
iconsidered=i
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
do l=1,ielem_ifca(i)
  godflux2=zero
  damp=lamx
  if (br2_yn.eq.2) damp=zero
  angle1=ielem_faceanglex(l,i)
  angle2=ielem_faceangley(l,i)
  nx=angle1
  ny=angle2
  nz=zero
  nall(1)=nx
  nall(2)=ny
  iqp=qp_line_n
  do ngp=1,iqp
    facex=l
    pointx=ngp
    b_code=0
    damp_loc=damp
    srf_speedrot=zero
    cleft=zero
    cright=zero
    cleft_rot=zero
    cright_rot=zero
    call calculate_bounded_viscous2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
         & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
    if (b_code.gt.0) damp_loc=zero
    if (lmach.eq.1)then
      call rotatef2d(n,cright_rot,cright,angle1,angle2)
      call rotatef2d(n,cleft_rot,cleft,angle1,angle2)
      leftv(1:nof_variables)=cleft_rot(1:nof_variables)
      rightv(1:nof_variables)=cright_rot(1:nof_variables)
      call lmacht2d_ideal(n,leftv,rightv)
      cleft_rot(1:nof_variables)=leftv(1:nof_variables)
      cright_rot(1:nof_variables)=rightv(1:nof_variables)
      call rotateb2d(n,cleft,cleft_rot,angle1,angle2)
      call rotateb2d(n,cright,cright_rot,angle1,angle2)
    end if
    leftv(1:nof_variables)=cleft(1:nof_variables)
    rightv(1:nof_variables)=cright(1:nof_variables)
    call get_visc_conduct_ideal(n,leftv,rightv,viscl,laml)
    if (turbulence.eq.1)then
      eddyfl=zero
      eddyfr=zero
      turbmv=zero
      etvm=zero
      if (turbulencemodel.eq.1)then
        turbmv(1)=cturbl(1)
        turbmv(2)=cturbr(1)
        eddyfl(2)=turbmv(1)
        eddyfr(2)=turbmv(2)
      else if (turbulencemodel.eq.2)then
        eddyfl(1)=ielem_walldist(i)
        eddyfl(2)=cturbl(1)
        eddyfl(3)=cturbl(2)
        eddyfl(4:5)=lcvgrad(1,1:2)
        eddyfl(6:7)=lcvgrad(2,1:2)
        eddyfl(8:9)=lcvgrad_t(1,1:2)
        eddyfl(10:11)=lcvgrad_t(2,1:2)
        eddyfr(1)=ielem_walldist(i)
        eddyfr(2)=cturbr(1)
        eddyfr(3)=cturbr(2)
        eddyfr(4:5)=rcvgrad(1,1:2)
        eddyfr(6:7)=rcvgrad(2,1:2)
        eddyfr(8:9)=rcvgrad_t(1,1:2)
        eddyfr(10:11)=rcvgrad_t(2,1:2)
      end if
      call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
    end if
    fxv=zero
    fyv=zero
    taul=zero
    tau=zero
    q=zero
    call cons2div_ideal(n,leftv,mp_pinfl,gammal)
    call cons2div_ideal(n,rightv,mp_pinfr,gammar)
    u12=oo2*(leftv(2)+rightv(2))
    v12=oo2*(leftv(3)+rightv(3))
    vdamp=4.0d0/3.0d0
    do k=1,nof_variables-1
      lcvgrad(k,1:2)=oo2*(lcvgrad(k,1:2)+rcvgrad(k,1:2))+damp_loc*((vdamp/abs(ielem_dih(l,i)))*nall(1:2)*(rightv(k+1)-leftv(k+1)))
    end do
    if (turbulence.eq.1)then
      q(1:2)=-oo2*(laml(3)+laml(4))*lcvgrad(dimensiona+1,1:2)
    else
      q(1:2)=-oo2*(laml(1)+laml(2))*lcvgrad(dimensiona+1,1:2)
    end if
    fxv(4)=fxv(4)-q(1)
    fyv(4)=fyv(4)-q(2)
    ux=lcvgrad(1,1)
    uy=lcvgrad(1,2)
    vx=lcvgrad(2,1)
    vy=lcvgrad(2,2)
    taul(1,1)=(4.0d0/3.0d0)*ux-(2.0d0/3.0d0)*vy
    taul(2,2)=(4.0d0/3.0d0)*vy-(2.0d0/3.0d0)*ux
    taul(1,2)=uy+vx
    taul(2,1)=taul(1,2)
    if (turbulence.eq.1)then
      tau=oo2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))*taul
    else
      tau=oo2*(viscl(1)+viscl(2))*taul
    end if
    do kc=2,3
      fxv(kc)=fxv(kc)+tau(1,kc-1)
      fyv(kc)=fyv(kc)+tau(2,kc-1)
    end do
    fxv(4)=fxv(4)+u12*tau(1,1)+v12*tau(1,2)
    fyv(4)=fyv(4)+u12*tau(2,1)+v12*tau(2,2)
    do iv=1,nof_variables
      hllcflux(iv)=nx*fxv(iv)+ny*fyv(iv)
    end do
    wface=weights_l(ngp)*ielem_surf(l,i)
    do iv=1,nof_variables
      godflux2(iv)=godflux2(iv)+hllcflux(iv)*wface
    end do
    if (nvt.gt.0)then
      if (turbulence.eq.1)then
        muturb=oo2*(viscl(1)+viscl(2))+oo2*(viscl(3)+viscl(4))
      else
        muturb=oo2*(viscl(1)+viscl(2))
      end if
      do iv=1,nvt
        hllcflux(nof_variables+iv)=muturb*oo2*((lcvgrad_t(iv,1)+rcvgrad_t(iv,1))*nx+(lcvgrad_t(iv,2)+rcvgrad_t(iv,2))*ny)
        if (turbulencemodel.eq.1) hllcflux(nof_variables+iv)=hllcflux(nof_variables+iv)/sigma
        if (b_code.eq.3) hllcflux(nof_variables+iv)=zero
        godflux2(nof_variables+iv)=godflux2(nof_variables+iv)+hllcflux(nof_variables+iv)*wface
      end do
    end if
  end do
  rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-godflux2(1:nof_variables)
  if (nvt.gt.0) rhst_val(1:nvt,i)=rhst_val(1:nvt,i)-godflux2(nof_variables+1:nof_variables+nvt)
end do

end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine calculate_fluxeshi_diffusive2d_bound_cell_ideal


subroutine calculate_fluxeshi_convective_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
		implicit none
		integer,intent(in)::n
		integer::ii
	







#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	call calculate_fluxeshi_convective_mood_inner_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	


#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_convective_mood_bound_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine calculate_fluxeshi_convective_mood

subroutine calculate_fluxeshi_convective_mood_inner_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ikas,igoflux,icaseb,jx,jx2,b_code
	real::sum_detect,norms,tempxx
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx,n_node
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext
	real,dimension(1:8,1:gpu_max_dim)::nodes_list
	real,dimension(1:gpu_max_dim)::cords
	integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf

	i=el_int(ii)
	iconsidered=i
	if (ielem_recalc(i).eq.1)then
            rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero
            end if


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
end subroutine calculate_fluxeshi_convective_mood_inner_cell


subroutine calculate_fluxeshi_convective_mood_bound_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ikas,igoflux,icaseb,jx,jx2,b_code
	real::sum_detect,norms,tempxx
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx,n_node
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext
	real,dimension(1:8,1:gpu_max_dim)::nodes_list
	real,dimension(1:gpu_max_dim)::cords
	integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf

	i=el_bnd(ii)
	iconsidered=i	

			if (ielem_recalc(i).eq.1)then
            rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero
            end if



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
								    call cordinates3(n,nodes_list,n_node,cords(1:3))
							    
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
end subroutine calculate_fluxeshi_convective_mood_bound_cell




subroutine calculate_fluxeshi_convective2d_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
		implicit none
		integer,intent(in)::n
		integer::ii
	


#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_interior	!for all the interior elements
	call calculate_fluxeshi_convective2d_mood_inner_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	
	
#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n) private(ii)
#else
!$omp barrier
!$omp do
#endif
	do ii=1,nof_bounded
	call calculate_fluxeshi_convective2d_mood_bound_cell(n,ii)
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine calculate_fluxeshi_convective2d_mood

subroutine calculate_fluxeshi_convective2d_mood_inner_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ikas,igoflux,icaseb,kxk,b_code
	real::sum_detect,norms
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx,n_node
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext,nodes_list
	real,dimension(1:gpu_max_dim)::cords
	integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf

	i=el_int(ii)

	if (ielem_recalc(i).eq.1)then
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero
	end if
	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)


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
end subroutine calculate_fluxeshi_convective2d_mood_inner_cell


subroutine calculate_fluxeshi_convective2d_mood_bound_cell(n,ii)
	implicit none
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::ii
	integer,intent(in)::n
	real,dimension(1:gpu_max_nvar_total)::godflux2, dg_vol_rec,rhllcflux,hllcflux
	integer::i,l,ngp,kmaxe,iqp,ikas,igoflux,icaseb,kxk,b_code
	real::sum_detect,norms
	real,dimension(1:gpu_max_qp_face)::weights_temp
	integer::iconsidered,facex,pointx,n_node
	real::angle1,angle2,nx,ny,nz,mp_source1,mp_source2,mp_source3
	real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot
	real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speedrot
	real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
	real,dimension(1:gpu_max_dim)::pox,poy,poz
	real,dimension(1:gpu_max_nvar)::srf_speed
	real,dimension(1:8,1:gpu_max_dim)::vext,nodes_list
	real,dimension(1:gpu_max_dim)::cords
	integer::ibfc,nfx,lfx,rowfx,nf,lf,rowf

	i=el_bnd(ii)
	iconsidered=i	

	if (ielem_recalc(i).eq.1)then
	rhs_val(:,i)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(:,i)=zero
	end if
	weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)



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
								    call cordinates2(n,nodes_list,n_node,cords(1:2))
							    
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
end subroutine calculate_fluxeshi_convective2d_mood_bound_cell



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
