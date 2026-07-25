module gradients
use declaration, only: gpu_max_dim, gpu_max_dof, gpu_max_extra_transport, &
     gpu_max_faces, gpu_max_fnodes, gpu_max_neighbours, gpu_max_neighbours2, &
     gpu_max_nodes, gpu_max_nvar, gpu_max_nvar_total, gpu_max_passive, &
     gpu_max_qp_all, gpu_max_qp_face, gpu_max_qp_volume, gpu_max_species, &
     gpu_max_turbulence, gpu_max_typesten
use library
use transform
use flow_operations
implicit none

 contains


 subroutine compute_gradients_inner_mean_lsq_viscous_acc(n)!check_all
!> @brief
!> this subroutine computes the gradients of the primitve variables of each interior cell using the least-squares
implicit none
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::sols1,sols2,leftv
integer::i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local,iconsidered
real::diff,coef,mp_pinfl,gammal

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_variables, dimensiona, r_gas, gamma, zero) &
!$omp& map(alloc: xmpielrank) &
!$omp& map(alloc: u_c_val, solhir, rec_wall, ielem_ggs, ielem_idegfree, ielem_inumneighbours, rec_local, rec_ihexl) &
!$omp& map(alloc: rec_ihexb, rec_ihexn, halo_offset, rec_invmat_stencilt, rec_gradf) &
!$omp& private(i,iconsidered,iq,ll,imax,nf,lf,rowf,k,ideg_local,var2,diff,coef,mp_pinfl,gammal) &
!$omp& private(sols1,sols2,leftv)
#else
!$omp do
#endif


do i=1,xmpielrank(n)

if (rec_wall(i).eq.0)then
if (ielem_ggs(i).eq.0)then

imax=ielem_inumneighbours(i)-1
iconsidered=i
ll=1
ideg_local=ielem_idegfree(i)
	leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
call cons2div_ideal(n,leftv,mp_pinfl,gammal)
sols1(1:nof_variables-1)=leftv(2:nof_variables)

do var2=1,nof_variables-1
   do k=1,ideg_local
      rec_gradf(var2,k,iconsidered)=zero
   end do
end do

do iq=1,imax
   if (rec_local(i).eq.0)then
      leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
   else
      if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
         leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
      else
         nf=rec_ihexn(1,iq+1,rec_local(i))
         lf=rec_ihexl(1,iq+1,i)
         rowf=halo_offset(nf) + lf - 1
         leftv(1:nof_variables)=solhir(rowf,1:nof_variables)
      end if
   end if

   call cons2div_ideal(n,leftv,mp_pinfl,gammal)
   sols2(1:nof_variables-1)=leftv(2:nof_variables)

   do k=1,ideg_local
      coef=rec_invmat_stencilt(k,iq,ll,i)
      do var2=1,nof_variables-1
         diff=sols2(var2)-sols1(var2)
         rec_gradf(var2,k,iconsidered)=rec_gradf(var2,k,iconsidered)+coef*diff
      end do
   end do
end do


end if
end if


end do


#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end subroutine compute_gradients_inner_mean_lsq_viscous_acc



subroutine compute_gradients_wall_mean_lsq_viscous_acc(n)!check all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each non-interior cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(gpu_max_nvar)::sols1,sols2
real,dimension(1:gpu_max_nvar,1:gpu_max_neighbours)::matrix_1
real,dimension(1:gpu_max_nvar,1:gpu_max_dof)::matrix_2
real,dimension(1:gpu_max_dof,1:gpu_max_nvar)::sol_m
integer::i,var2,k0,g0,ttk,ivvm,iq,lq
real::attt,mp_pinfl,gammal
integer::ll,imax,nf,lf,rowf,iconsidered,number_of_dog
real,dimension(1:gpu_max_nvar)::leftv


#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(idegfree, nof_variables, nof_species, dimensiona, r_gas, thermal) &
!$omp& firstprivate(catalytic_wall, wall_temp, tolbig, zero) &
!$omp& map(alloc: xmpielrank) &
!$omp& map(alloc: u_c_val, solhir, rec_wall, ielem_ggs, ielem_idegfree, ielem_inumneighbours, rec_local, rec_ihexl) &
!$omp& map(alloc: rec_ihexb, rec_ihexn, halo_offset, rec_k0, rec_g0, rec_volume_w) &
!$omp& map(alloc: rec_weightl, rec_stencils, rec_wallcoeff, catalytic_con) &
!$omp& map(alloc: rec_vellsq, rec_tempsq, rec_velinvlsqmat, rec_tempsqmat, rec_wallcoefg, rec_gradf) &
!$omp& private(iconsidered,iq,ll,imax,nf,lf,rowf,k0,g0,ttk,ivvm,lq,number_of_dog) &
!$omp& private(var2,attt,mp_pinfl,gammal,sols1,sols2,leftv,matrix_1,matrix_2,sol_m)
#else
!$omp do
#endif


do i=1,xmpielrank(n)

if (rec_wall(i).gt.0)then
if (ielem_ggs(i).eq.0)then

imax=ielem_inumneighbours(i)-1
number_of_dog=ielem_idegfree(i)
ll=1
iconsidered=i

	    ll=1
	    k0=rec_k0(rec_wall(i))
	    g0=rec_g0(rec_wall(i))

	     matrix_1=zero;matrix_2=zero;sol_m=zero;
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
		call cons2div_ideal(n,leftv,mp_pinfl,gammal)

	       sols1(1:nof_variables-1)=leftv(2:nof_variables)
! 	       sols1(1)=leftv(5)/(leftv(1)*r_gas)

	      do iq=1,imax
			if (rec_local(i).eq.0)then
			leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
			else
				if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
				else
! 				leftv(1:nof_variables)=iexsolhir(rec_ihexn(1,iq+1,i))%sol(rec_ihexl(1,iq+1,i),1:nof_variables)
				    nf=rec_ihexn(1,iq+1,rec_local(i))
					lf=rec_ihexl(1,iq+1,i)
					rowf=halo_offset(nf) + lf - 1
					leftv(1:nof_variables)=solhir(rowf, 1:nof_variables)



				end if
			end if

		 call cons2div_ideal(n,leftv,mp_pinfl,gammal)
	       sols2(1:nof_variables-1)=leftv(2:nof_variables)
! 	       !sols2(1)=leftv(5)/(leftv(1)*r_gas)
  	        matrix_1(1:nof_variables-1,iq)= &
             rec_volume_w(1,iq+1,rec_wall(i))*rec_weightl(1,iq,rec_wall(i))* &
             (sols2(1:nof_variables-1)-sols1(1:nof_variables-1))

  	        !velocity gradients
  	        matrix_1(1:dimensiona,iq)=matrix_1(1:dimensiona,iq)+ &
             (sols1(1:dimensiona)*rec_stencils(ll,iq,k0,rec_wall(i)))/ &
             rec_wallcoeff(k0,rec_wall(i))

  	        !temperature now check

			if (thermal.eq.1)then
  	        matrix_1(dimensiona+1:nof_variables-nof_species-1,iq)= &
             matrix_1(dimensiona+1:nof_variables-nof_species-1,iq)+ &
             (sols1(dimensiona+1:nof_variables-nof_species-1)* &
             rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i))- &
             (wall_temp*rec_stencils(ll,iq,k0,rec_wall(i)))/ &
             rec_wallcoeff(k0,rec_wall(i))

  	        end if


   	        if (catalytic_wall.eq.1)then
   	        matrix_1(dimensiona+3:nof_variables-1,iq)= &
             matrix_1(dimensiona+3:nof_variables-1,iq)+ &
             (sols1(dimensiona+3:nof_variables-1)* &
             rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i))- &
             (catalytic_con(1:nof_species)*rec_stencils(ll,iq,k0,rec_wall(i)))/ &
             rec_wallcoeff(k0,rec_wall(i))



   	        end if



			end do
			do var2=1,nof_variables-1
		  matrix_2(var2,1:number_of_dog-1)=zero
		  if (var2.le.dimensiona)then		!velocity gradients
		  do iq=1,imax

		      do lq=1,number_of_dog-1
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
		      end do

		  end do
		  end if
		  if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then	!temperature gradients
		  do iq=1,imax
		     do lq=1,number_of_dog-1
		     if (thermal.eq.1)then
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
			else
			matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_tempsq(iq,lq,rec_wall(i))
			end if
		      end do
		  end do
		  end if
		   if (var2.gt.nof_variables-nof_species-1)then					!species
		   do iq=1,imax
		     do lq=1,number_of_dog-1

		     if (catalytic_wall.eq.1)then
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
		     else
			matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_tempsq(iq,lq,rec_wall(i))
			end if

		      end do
		  end do
		   end if

		if (var2.le.dimensiona)then
		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		end if
		if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then
		if (thermal.eq.1)then
		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		else
		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_tempsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		end if
		end if

		if (var2.gt.nof_variables-nof_species-1)then

		if (catalytic_wall.eq.1)then
			do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do

		else

		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_tempsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		end if
		end if

	     end do

		do var2=1,nof_variables-1

		 if (var2.le.dimensiona)then		!velocity gradients
		 rec_gradf(var2,1:idegfree,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		  end do
		  attt=zero
		  attt=-sols1(var2)
			  do ttk=1,number_of_dog
				    if (ttk.ne.k0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradf(var2,k0,iconsidered)=attt


		end if

		if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then	!temperature gradients
		if (thermal.eq.1)then
		rec_gradf(var2,1:idegfree,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		  end do
		  attt=zero
		  attt=wall_temp-sols1(var2)
			  do ttk=1,number_of_dog
				    if (ttk.ne.k0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradf(var2,k0,iconsidered)=attt

		else
		ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.g0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		    end do
		    attt=zero

			  do ttk=1,number_of_dog
				    if (ttk.ne.g0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoefg(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoefg(g0,rec_wall(i))
			    rec_gradf(var2,g0,iconsidered)=attt
		end if
		end if
		if (var2.gt.nof_variables-nof_species-1)then				!species gradients

			if (catalytic_wall.eq.1)then
		rec_gradf(var2,1:idegfree,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		  end do
		  attt=zero
		  attt=catalytic_con(var2-dimensiona-2)-sols1(var2)


			  do ttk=1,number_of_dog
				    if (ttk.ne.k0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradf(var2,k0,iconsidered)=attt

		else



			ivvm=0
				do ttk=1,number_of_dog
						if (ttk.eq.g0) cycle
						ivvm=ivvm+1
							rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
				end do
				attt=zero

				do ttk=1,number_of_dog
						if (ttk.ne.g0) &
					attt=attt-rec_gradf(var2,ttk,iconsidered)*&
								rec_wallcoefg(ttk,rec_wall(i))
				end do
					attt=attt/rec_wallcoefg(g0,rec_wall(i))
					rec_gradf(var2,g0,iconsidered)=attt

		end if

		end if


		end do



end if
end if

end do


#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif









end subroutine compute_gradients_wall_mean_lsq_viscous_acc


subroutine compute_gradients_inner_mean_ggs_viscous_acc(n)
!> @brief
!> Ideal-gas Green-Gauss primitive gradients for interior cells.
implicit none
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::sols1,sols2,leftv
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::sols_f
real,dimension(1:gpu_max_dim)::normal_all
real::oov2,mp_pinfl,gammal,angle1,angle2
integer::ii,i,j,k

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, dimensiona, zero, oo2) &
!$omp& map(alloc: el_int, rec_wall, ielem_ggs, u_c_val, ielem_totvolume) &
!$omp& map(alloc: ielem_ifca, ielem_faceanglex, ielem_faceangley, ielem_ineigh, ielem_surf, rec_grads) &
!$omp& private(ii,i,j,k,sols1,sols2,leftv,sols_f,normal_all,oov2,mp_pinfl,gammal,angle1,angle2)
#else
!$omp do private(ii,i,j,k,sols1,sols2,leftv,sols_f,normal_all,oov2,mp_pinfl,gammal,angle1,angle2)
#endif
do ii=1,nof_interior
  i=el_int(ii)
  if ((rec_wall(i).eq.0).and.(ielem_ggs(i).eq.1)) then
    sols_f=zero
    sols1=zero
    sols2=zero
    rec_grads(:,:,i)=zero
    oov2=1.0d0/ielem_totvolume(i)

    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
    call cons2div_ideal(n,leftv,mp_pinfl,gammal)
    sols1(1:nof_variables-1)=leftv(2:nof_variables)

    do j=1,ielem_ifca(i)
      angle1=ielem_faceanglex(j,i)
      angle2=ielem_faceangley(j,i)
      if (dimensiona.eq.3) then
        normal_all(1)=cos(angle1)*sin(angle2)
        normal_all(2)=sin(angle1)*sin(angle2)
        normal_all(3)=cos(angle2)
      else
        normal_all(1)=angle1
        normal_all(2)=angle2
      end if

      leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
      call cons2div_ideal(n,leftv,mp_pinfl,gammal)
      sols2(1:nof_variables-1)=leftv(2:nof_variables)

      do k=1,dimensiona
        sols_f(1:nof_variables-1,k)=sols_f(1:nof_variables-1,k)+ &
          (oo2*(sols2(1:nof_variables-1)+sols1(1:nof_variables-1))*normal_all(k)*ielem_surf(j,i)*oov2)
      end do
    end do

    do k=1,dimensiona
      rec_grads(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
    end do
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_inner_mean_ggs_viscous_acc


subroutine compute_gradients_inner_turb_ggs_viscous_acc(n)
!> @brief
!> Ideal-gas Green-Gauss turbulence/passive gradients for interior cells.
implicit none
integer,intent(in)::n
real,dimension(gpu_max_extra_transport)::sols1,sols2
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::sols_f
real,dimension(gpu_max_dim)::normal_all
real::oov2,angle1,angle2
integer::ii,i,j,k,var2,nvt

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_interior, nof_variables, turbulenceequations, passivescalar, dimensiona, zero, oo2) &
!$omp& map(alloc: el_int, rec_wall, u_c_val, u_ct_val, ielem_totvolume) &
!$omp& map(alloc: ielem_ifca, ielem_faceanglex, ielem_faceangley, ielem_ineigh, ielem_surf) &
!$omp& map(alloc: rec_gradientsturb, rec_grads) &
!$omp& private(ii,i,j,k,var2,nvt,sols1,sols2,sols_f,normal_all,oov2,angle1,angle2)
#else
!$omp do private(ii,i,j,k,var2,nvt,sols1,sols2,sols_f,normal_all,oov2,angle1,angle2)
#endif
do ii=1,nof_interior
  i=el_int(ii)
  nvt=turbulenceequations+passivescalar
  if ((rec_wall(i).eq.0).and.(nvt.gt.0)) then
    sols_f=zero
    sols1=zero
    sols2=zero
    oov2=1.0d0/ielem_totvolume(i)
    sols1(1:nvt)=u_ct_val(1,1:nvt,i)/u_c_val(1,1,i)

    do j=1,ielem_ifca(i)
      angle1=ielem_faceanglex(j,i)
      angle2=ielem_faceangley(j,i)
      if (dimensiona.eq.3) then
        normal_all(1)=cos(angle1)*sin(angle2)
        normal_all(2)=sin(angle1)*sin(angle2)
        normal_all(3)=cos(angle2)
      else
        normal_all(1)=angle1
        normal_all(2)=angle2
      end if

      sols2(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(j,i))/u_c_val(1,1,ielem_ineigh(j,i))
      do k=1,dimensiona
        sols_f(1:nvt,k)=sols_f(1:nvt,k)+(oo2*(sols2(1:nvt)+sols1(1:nvt))*normal_all(k)*ielem_surf(j,i)*oov2)
      end do
    end do

    do var2=1,nvt
      rec_gradientsturb(1,1:dimensiona,var2,i)=sols_f(var2,1:dimensiona)
      rec_grads(dimensiona+1+var2,1:dimensiona,i)=sols_f(var2,1:dimensiona)
    end do
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_inner_turb_ggs_viscous_acc


subroutine compute_gradients_turb_lsq_acc(n)
!> @brief
!> Least-squares turbulence/passive reconstruction gradients for all local cells.
implicit none
integer,intent(in)::n
real,dimension(gpu_max_extra_transport)::sols1,sols2
integer::i,var2,ll,iq,nf,lf,rowf,k,imax,ideg_local,nvt,number_of_dog,number_of_nei
real::diff,coef

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, ees, idegfree2, numneighbours2, zero) &
!$omp& map(alloc: xmpielrank, ielem_idegfree, ielem_inumneighbours, ielem_admis) &
!$omp& map(alloc: u_ct_val, solhir, rec_local, rec_ihexl, rec_ihexb, rec_ihexn, halo_offset) &
!$omp& map(alloc: rec_invmat_stencilt, rec_gradients2, rec_ihexlc, rec_invmat_stenciltc) &
!$omp& map(alloc: rec_gradientsc2, rec_ihexbc, rec_ihexnc) &
!$omp& private(i,var2,ll,iq,nf,lf,rowf,k,imax,ideg_local,nvt,number_of_dog,number_of_nei) &
!$omp& private(diff,coef,sols1,sols2)
#else
!$omp do private(i,var2,ll,iq,nf,lf,rowf,k,imax,ideg_local,nvt,number_of_dog,number_of_nei,diff,coef,sols1,sols2)
#endif
do i=1,xmpielrank(n)
  nvt=turbulenceequations+passivescalar
  if (nvt.gt.0) then
    number_of_dog=ielem_idegfree(i)
    number_of_nei=ielem_inumneighbours(i)
    imax=number_of_nei-1
    sols1=zero
    sols2=zero
    sols1(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,1,i))

    if (rec_local(i).eq.0)then
      do ll=1,ielem_admis(i)
        if ((ees.ne.5).or.(ll.eq.1))then
          ideg_local=number_of_dog
          do var2=1,nvt
            do k=1,ideg_local
              rec_gradients2(ll,k,var2,i)=zero
            end do
          end do
          do iq=1,imax
            sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(ll,iq+1,i))
            do k=1,ideg_local
              coef=rec_invmat_stencilt(k,iq,ll,i)
              do var2=1,nvt
                diff=sols2(var2)-sols1(var2)
                rec_gradients2(ll,k,var2,i)=rec_gradients2(ll,k,var2,i)+coef*diff
              end do
            end do
          end do
        else
          ideg_local=idegfree2
          do var2=1,nvt
            do k=1,ideg_local
              rec_gradientsc2(ll,k,var2,i)=zero
            end do
          end do
          do iq=1,numneighbours2-1
            sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexlc(ll,iq+1,i))
            do k=1,ideg_local
              coef=rec_invmat_stenciltc(k,iq,ll,i)
              do var2=1,nvt
                diff=sols2(var2)-sols1(var2)
                rec_gradientsc2(ll,k,var2,i)=rec_gradientsc2(ll,k,var2,i)+coef*diff
              end do
            end do
          end do
        end if
      end do
    else
      do ll=1,ielem_admis(i)
        if ((ees.ne.5).or.(ll.eq.1))then
          ideg_local=number_of_dog
          do var2=1,nvt
            do k=1,ideg_local
              rec_gradients2(ll,k,var2,i)=zero
            end do
          end do
          do iq=1,imax
            if (rec_ihexb(ll,iq+1,rec_local(i)).eq.n)then
              sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(ll,iq+1,i))
            else
              nf=rec_ihexn(ll,iq+1,rec_local(i))
              lf=rec_ihexl(ll,iq+1,i)
              rowf=halo_offset(nf)+lf-1
              sols2(1:nvt)=solhir(rowf,nof_variables+1:nof_variables+nvt)
            end if
            do k=1,ideg_local
              coef=rec_invmat_stencilt(k,iq,ll,i)
              do var2=1,nvt
                diff=sols2(var2)-sols1(var2)
                rec_gradients2(ll,k,var2,i)=rec_gradients2(ll,k,var2,i)+coef*diff
              end do
            end do
          end do
        else
          ideg_local=idegfree2
          do var2=1,nvt
            do k=1,ideg_local
              rec_gradientsc2(ll,k,var2,i)=zero
            end do
          end do
          do iq=1,numneighbours2-1
            if (rec_ihexbc(ll,iq+1,rec_local(i)).eq.n)then
              sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexlc(ll,iq+1,i))
            else
              nf=rec_ihexnc(ll,iq+1,rec_local(i))
              lf=rec_ihexlc(ll,iq+1,i)
              rowf=halo_offset(nf)+lf-1
              sols2(1:nvt)=solhir(rowf,nof_variables+1:nof_variables+nvt)
            end if
            do k=1,ideg_local
              coef=rec_invmat_stenciltc(k,iq,ll,i)
              do var2=1,nvt
                diff=sols2(var2)-sols1(var2)
                rec_gradientsc2(ll,k,var2,i)=rec_gradientsc2(ll,k,var2,i)+coef*diff
              end do
            end do
          end do
        end if
      end do
    end if
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_turb_lsq_acc


subroutine compute_gradients_turb_lsq_viscous_acc(n)
!> @brief
!> Least-squares viscous turbulence/passive gradients for non-wall cells.
implicit none
integer,intent(in)::n
real,dimension(gpu_max_extra_transport)::sols1,sols2
integer::i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local,nvt,number_of_dog,number_of_nei
real::diff,coef

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, zero) &
!$omp& map(alloc: xmpielrank, rec_wall, ielem_ggs, ielem_idegfree, ielem_inumneighbours) &
!$omp& map(alloc: u_c_val, u_ct_val, solhir, rec_local, rec_ihexl, rec_ihexb, rec_ihexn, halo_offset) &
!$omp& map(alloc: rec_invmat_stencilt, rec_gradientsturb) &
!$omp& private(i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local,nvt,number_of_dog,number_of_nei,diff,coef,sols1,sols2)
#else
!$omp do private(i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local,nvt,number_of_dog,number_of_nei,diff,coef,sols1,sols2)
#endif
do i=1,xmpielrank(n)
  nvt=turbulenceequations+passivescalar
  if ((nvt.gt.0).and.(rec_wall(i).eq.0).and.(ielem_ggs(i).eq.0)) then
    number_of_dog=ielem_idegfree(i)
    number_of_nei=ielem_inumneighbours(i)
    imax=number_of_nei-1
    ll=1
    ideg_local=number_of_dog
    sols1=zero
    sols2=zero
    sols1(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,1,i))/u_c_val(1,1,rec_ihexl(1,1,i))
    do var2=1,nvt
      do k=1,ideg_local
        rec_gradientsturb(1,k,var2,i)=zero
      end do
    end do

    do iq=1,imax
      if (rec_local(i).eq.0)then
        sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,iq+1,i))/u_c_val(1,1,rec_ihexl(1,iq+1,i))
      else
        if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
          sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,iq+1,i))/u_c_val(1,1,rec_ihexl(1,iq+1,i))
        else
          nf=rec_ihexn(ll,iq+1,rec_local(i))
          lf=rec_ihexl(ll,iq+1,i)
          rowf=halo_offset(nf)+lf-1
          sols2(1:nvt)=solhir(rowf,nof_variables+1:nof_variables+nvt)/solhir(rowf,1)
        end if
      end if
      do k=1,ideg_local
        coef=rec_invmat_stencilt(k,iq,ll,i)
        do var2=1,nvt
          diff=sols2(var2)-sols1(var2)
          rec_gradientsturb(1,k,var2,i)=rec_gradientsturb(1,k,var2,i)+coef*diff
        end do
      end do
    end do
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_turb_lsq_viscous_acc


subroutine compute_gradients_wall_turb_lsq_viscous_acc(n)
!> @brief
!> Least-squares viscous turbulence/passive gradients for wall cells.
implicit none
integer,intent(in)::n
real,dimension(1:gpu_max_extra_transport)::sols1,sols2
real,dimension(1:gpu_max_extra_transport,1:gpu_max_neighbours)::matrix_1
real,dimension(1:gpu_max_extra_transport,1:gpu_max_dof)::matrix_2
real,dimension(1:gpu_max_dof,1:gpu_max_extra_transport)::sol_m
real,dimension(1:gpu_max_extra_transport)::matrix_3
integer::i,var2,ii,k0,ttk,ivvm,iq,lq,imax,ll,nf,lf,rowf,number_of_dog,nvt
real::attt

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, turbulencemodel) &
!$omp& firstprivate(visc, beta_i1, tolbig, zero) &
!$omp& map(alloc: xmpielrank, rec_wall, ielem_ggs, ielem_idegfree, ielem_inumneighbours, rec_k0, ielem_walldist) &
!$omp& map(alloc: u_c_val, u_ct_val, solhir, rec_local, rec_ihexl, rec_ihexb, rec_ihexn, halo_offset) &
!$omp& map(alloc: rec_volume_w, rec_weightl, rec_stencils, rec_wallcoeff, rec_vellsq, rec_velinvlsqmat) &
!$omp& map(alloc: rec_gradientsturb) &
!$omp& private(i,var2,ii,k0,ttk,ivvm,iq,lq,imax,ll,nf,lf,rowf,number_of_dog,nvt,attt) &
!$omp& private(sols1,sols2,matrix_1,matrix_2,sol_m,matrix_3)
#else
!$omp do private(i,var2,ii,k0,ttk,ivvm,iq,lq,imax,ll,nf,lf,rowf,number_of_dog,nvt,attt,sols1,sols2,matrix_1,matrix_2,sol_m,matrix_3)
#endif
do i=1,xmpielrank(n)
  nvt=turbulenceequations+passivescalar
  if ((nvt.gt.0).and.(rec_wall(i).gt.0).and.(ielem_ggs(i).eq.0)) then
    imax=ielem_inumneighbours(i)-1
    ll=1
    number_of_dog=ielem_idegfree(i)
    k0=rec_k0(rec_wall(i))
    sols1=zero
    sols2=zero
    matrix_1=zero
    matrix_2=zero
    sol_m=zero
    matrix_3=zero
    sols1(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,1,i))/u_c_val(1,1,rec_ihexl(1,1,i))

    do iq=1,imax
      if (rec_local(i).eq.0)then
        sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,iq+1,i))/u_c_val(1,1,rec_ihexl(1,iq+1,i))
      else
        if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
          sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,iq+1,i))/u_c_val(1,1,rec_ihexl(1,iq+1,i))
        else
          nf=rec_ihexn(1,iq+1,rec_local(i))
          lf=rec_ihexl(1,iq+1,i)
          rowf=halo_offset(nf)+lf-1
          sols2(1:nvt)=solhir(rowf,nof_variables+1:nof_variables+nvt)/solhir(rowf,1)
        end if
      end if

      matrix_1(1:nvt,iq)=rec_volume_w(1,iq+1,rec_wall(i))*rec_weightl(1,iq,rec_wall(i))*(sols2(1:nvt)-sols1(1:nvt))
      matrix_1(1:nvt,iq)=matrix_1(1:nvt,iq)+((sols1(1:nvt)*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))
    end do

    matrix_3(1:nvt)=-sols1(1:nvt)
    if ((turbulencemodel.eq.2).and.(nvt.ge.2)) then
      matrix_3(2)=60.0d0*visc/(beta_i1*(ielem_walldist(i)**2))
    end if

    do var2=1,nvt
      matrix_2(var2,1:number_of_dog-1)=zero
      do iq=1,imax
        do lq=1,number_of_dog-1
          matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
        end do
      end do

      do ivvm=1,number_of_dog-1
        sol_m(ivvm,var2)=zero
        do lq=1,number_of_dog-1
          sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
        end do
      end do
    end do

    do var2=1,nvt
      rec_gradientsturb(1,1:number_of_dog,var2,i)=-tolbig
      ivvm=0
      do ttk=1,number_of_dog
        if (ttk.eq.k0) cycle
        ivvm=ivvm+1
        rec_gradientsturb(1,ttk,var2,i)=sol_m(ivvm,var2)
      end do
      attt=matrix_3(var2)
      do ttk=1,number_of_dog
        if (ttk.ne.k0) attt=attt-rec_gradientsturb(1,ttk,var2,i)*rec_wallcoeff(ttk,rec_wall(i))
      end do
      attt=attt/rec_wallcoeff(k0,rec_wall(i))
      rec_gradientsturb(1,k0,var2,i)=attt
    end do
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_wall_turb_lsq_viscous_acc


subroutine compute_gradients_mix_mean_ggs_viscous_acc(n)
!> @brief
!> Ideal-gas Green-Gauss primitive gradients for bounded/mixed cells.
implicit none
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::sols1,sols2,leftv,rightv,srf_speed,srf_speedrot,phi_f
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::sols_f
real,dimension(1:gpu_max_dim)::pox,poy,poz,normal_all
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_nvar_total)::cright_rot,cleft_rot
real::oov2,mp_pinfl,gammal,angle1,angle2,nx,ny,nz
integer::ii,i,j,k,b_code,facex,n_node,nf,lf,rowf,ibfc,nvt

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, turbulenceequations, passivescalar, dimensiona, zero, oo2) &
!$omp& firstprivate(per_rot, angle_per, thermal, wall_temp, initcond, boundtype, gamma, tolsmall, turbulence) &
!$omp& map(alloc: el_bnd, ielem_ggs, u_c_val, u_ct_val, solhir, ielem_totvolume, rec_grads) &
!$omp& map(alloc: ielem_ifca, ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_ineighb, ielem_ibounds) &
!$omp& map(alloc: ibound_icode, ielem_ineigh, rec_local, rec_ihexn, rec_ihexl, halo_offset, ielem_indexi) &
!$omp& map(alloc: rec_mrf, rec_rotvel, ielem_types_faces, ielem_nodes_faces, inoder4_cord) &
!$omp& private(ii,i,j,k,b_code,facex,n_node,nf,lf,rowf,ibfc,nvt) &
!$omp& private(sols1,sols2,leftv,rightv,srf_speed,srf_speedrot,phi_f,sols_f,pox,poy,poz,normal_all) &
!$omp& private(cturbl,cturbr,cright_rot,cleft_rot,oov2,mp_pinfl,gammal,angle1,angle2,nx,ny,nz)
#else
!$omp do private(ii,i,j,k,b_code,facex,n_node,nf,lf,rowf,ibfc,nvt,sols1,sols2,leftv,rightv,srf_speed,srf_speedrot,phi_f,sols_f,pox,poy,poz,normal_all,cturbl,cturbr,cright_rot,cleft_rot,oov2,mp_pinfl,gammal,angle1,angle2,nx,ny,nz)
#endif
do ii=1,nof_bounded
  i=el_bnd(ii)
  if (ielem_ggs(i).eq.1) then
    nvt=turbulenceequations+passivescalar
    sols_f=zero
    sols1=zero
    sols2=zero
    rec_grads(:,:,i)=zero
    oov2=1.0d0/ielem_totvolume(i)

    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
    call cons2div_ideal(n,leftv,mp_pinfl,gammal)
    sols1(1:nof_variables-1)=leftv(2:nof_variables)

    do j=1,ielem_ifca(i)
      facex=j
      b_code=0
      angle1=ielem_faceanglex(j,i)
      angle2=ielem_faceangley(j,i)
      if (dimensiona.eq.3) then
        normal_all(1)=cos(angle1)*sin(angle2)
        normal_all(2)=sin(angle1)*sin(angle2)
        normal_all(3)=cos(angle2)
      else
        normal_all(1)=angle1
        normal_all(2)=angle2
      end if
      nx=normal_all(1)
      ny=normal_all(2)
      nz=zero
      if (dimensiona.eq.3) nz=normal_all(3)
      srf_speed=zero
      srf_speedrot=zero
      if ((dimensiona.eq.3).and.(rec_mrf(i).eq.1))then
        srf_speed(2:4)=rec_rotvel(j,1,1:3,i)
        call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
      end if

      if (ielem_ineighb(j,i).eq.n)then
        if (ielem_ibounds(j,i).gt.0)then
          if ((ibound_icode(ielem_ibounds(j,i)).eq.5).or.(ibound_icode(ielem_ibounds(j,i)).eq.50))then
            sols2(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
            if ((dimensiona.eq.3).and.(per_rot.eq.1).and.(ibound_icode(ielem_ibounds(j,i)).eq.50)) &
              sols2(2:4)=rotate_per_1(sols2(2:4),ibound_icode(ielem_ibounds(j,i)),angle_per)
          else
            b_code=ibound_icode(ielem_ibounds(j,i))
            if (dimensiona.eq.3)then
              if (ielem_types_faces(facex,i).eq.5)then
                n_node=4
              else
                n_node=3
              end if
            else
              n_node=2
            end if
            pox=zero
            poy=zero
            poz=zero
            do k=1,n_node
              pox(1)=pox(1)+inoder4_cord(1,ielem_nodes_faces(facex,k,i))
              poy(1)=poy(1)+inoder4_cord(2,ielem_nodes_faces(facex,k,i))
              if (dimensiona.eq.3) poz(1)=poz(1)+inoder4_cord(3,ielem_nodes_faces(facex,k,i))
            end do
            pox(1)=pox(1)/real(n_node)
            poy(1)=poy(1)/real(n_node)
            if (dimensiona.eq.3) poz(1)=poz(1)/real(n_node)
            leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
            rightv=zero
            cturbl=zero
            cturbr=zero
            cright_rot=zero
            cleft_rot=zero
            if (nvt.gt.0) cturbl(1:nvt)=u_ct_val(1,1:nvt,i)
            ibfc=0
            if (dimensiona.eq.3)then
              call boundarys_ideal(n,b_code,i,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
                   & cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
            else
              call boundarys2d_ideal(n,b_code,i,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
                   & cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
            end if
            sols2(1:nof_variables)=rightv(1:nof_variables)
          end if
        else
          sols2(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
        end if
      else
        nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
        lf=rec_ihexl(1,ielem_indexi(j,i),i)
        rowf=halo_offset(nf)+lf-1
        sols2(1:nof_variables)=solhir(rowf,1:nof_variables)
        if ((dimensiona.eq.3).and.(ielem_ibounds(j,i).gt.0).and.(per_rot.eq.1))then
          if (ibound_icode(ielem_ibounds(j,i)).eq.50) &
            sols2(2:4)=rotate_per_1(sols2(2:4),ibound_icode(ielem_ibounds(j,i)),angle_per)
        end if
      end if

      leftv(1:nof_variables)=sols2(1:nof_variables)
      call cons2div_ideal(n,leftv,mp_pinfl,gammal)
      sols2(1:nof_variables-1)=leftv(2:nof_variables)
      if ((b_code.eq.4).and.(thermal.eq.1)) sols2(dimensiona+1:nof_variables-1)=wall_temp

      do k=1,dimensiona
        sols_f(1:nof_variables-1,k)=sols_f(1:nof_variables-1,k)+ &
          (oo2*(sols2(1:nof_variables-1)+sols1(1:nof_variables-1))*normal_all(k)*ielem_surf(j,i)*oov2)
      end do
    end do

    do k=1,dimensiona
      rec_grads(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
    end do
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_mix_mean_ggs_viscous_acc


subroutine compute_gradients_mix_turb_ggs_viscous_acc(n)
!> @brief
!> Ideal-gas Green-Gauss turbulence/passive gradients for bounded/mixed cells.
implicit none
integer,intent(in)::n
real,dimension(gpu_max_extra_transport)::sols1,sols2,cturbl,cturbr
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::sols_f
real,dimension(gpu_max_dim)::normal_all,pox,poy,poz
real,dimension(1:gpu_max_nvar)::leftv,rightv,srf_speed,srf_speedrot
real,dimension(1:gpu_max_nvar_total)::cright_rot,cleft_rot
real::oov2,angle1,angle2,nx,ny,nz
integer::ii,i,j,k,var2,b_code,facex,n_node,nf,lf,rowf,ibfc,nvt

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) &
!$omp& firstprivate(nof_bounded, nof_variables, turbulenceequations, passivescalar, dimensiona, zero, oo2) &
!$omp& firstprivate(initcond, boundtype, gamma, tolsmall, turbulence) &
!$omp& map(alloc: el_bnd, u_c_val, u_ct_val, solhir, ielem_totvolume, rec_gradientsturb, rec_grads) &
!$omp& map(alloc: ielem_ifca, ielem_faceanglex, ielem_faceangley, ielem_surf, ielem_ineighb, ielem_ibounds) &
!$omp& map(alloc: ibound_icode, ielem_ineigh, rec_local, rec_ihexn, rec_ihexl, halo_offset, ielem_indexi) &
!$omp& map(alloc: rec_mrf, ielem_types_faces, ielem_nodes_faces, inoder4_cord) &
!$omp& private(ii,i,j,k,var2,b_code,facex,n_node,nf,lf,rowf,ibfc,nvt) &
!$omp& private(sols1,sols2,cturbl,cturbr,sols_f,normal_all,pox,poy,poz,leftv,rightv,srf_speed,srf_speedrot) &
!$omp& private(cright_rot,cleft_rot,oov2,angle1,angle2,nx,ny,nz)
#else
!$omp do private(ii,i,j,k,var2,b_code,facex,n_node,nf,lf,rowf,ibfc,nvt,sols1,sols2,cturbl,cturbr,sols_f,normal_all,pox,poy,poz,leftv,rightv,srf_speed,srf_speedrot,cright_rot,cleft_rot,oov2,angle1,angle2,nx,ny,nz)
#endif
do ii=1,nof_bounded
  i=el_bnd(ii)
  nvt=turbulenceequations+passivescalar
  if (nvt.gt.0) then
    sols_f=zero
    sols1=zero
    sols2=zero
    oov2=1.0d0/ielem_totvolume(i)
    sols1(1:nvt)=u_ct_val(1,1:nvt,i)/u_c_val(1,1,i)

    do j=1,ielem_ifca(i)
      facex=j
      b_code=0
      angle1=ielem_faceanglex(j,i)
      angle2=ielem_faceangley(j,i)
      if (dimensiona.eq.3) then
        normal_all(1)=cos(angle1)*sin(angle2)
        normal_all(2)=sin(angle1)*sin(angle2)
        normal_all(3)=cos(angle2)
      else
        normal_all(1)=angle1
        normal_all(2)=angle2
      end if
      nx=normal_all(1)
      ny=normal_all(2)
      nz=zero
      if (dimensiona.eq.3) nz=normal_all(3)

      if (ielem_ineighb(j,i).eq.n)then
        if (ielem_ibounds(j,i).gt.0)then
          if ((ibound_icode(ielem_ibounds(j,i)).eq.5).or.(ibound_icode(ielem_ibounds(j,i)).eq.50))then
            sols2(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(j,i))/u_c_val(1,1,ielem_ineigh(j,i))
          else
            b_code=ibound_icode(ielem_ibounds(j,i))
            if (dimensiona.eq.3)then
              if (ielem_types_faces(facex,i).eq.5)then
                n_node=4
              else
                n_node=3
              end if
            else
              n_node=2
            end if
            pox=zero
            poy=zero
            poz=zero
            do k=1,n_node
              pox(1)=pox(1)+inoder4_cord(1,ielem_nodes_faces(facex,k,i))
              poy(1)=poy(1)+inoder4_cord(2,ielem_nodes_faces(facex,k,i))
              if (dimensiona.eq.3) poz(1)=poz(1)+inoder4_cord(3,ielem_nodes_faces(facex,k,i))
            end do
            pox(1)=pox(1)/real(n_node)
            poy(1)=poy(1)/real(n_node)
            if (dimensiona.eq.3) poz(1)=poz(1)/real(n_node)
            leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
            rightv=zero
            cturbl=zero
            cturbr=zero
            cturbl(1:nvt)=u_ct_val(1,1:nvt,i)
            cright_rot=zero
            cleft_rot=zero
            srf_speed=zero
            srf_speedrot=zero
            ibfc=0
            if (dimensiona.eq.3)then
              call boundarys_ideal(n,b_code,i,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
                   & cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
            else
              call boundarys2d_ideal(n,b_code,i,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
                   & cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
            end if
            sols2(1:nvt)=cturbr(1:nvt)/rightv(1)
          end if
        else
          sols2(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(j,i))/u_c_val(1,1,ielem_ineigh(j,i))
        end if
      else
        nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
        lf=rec_ihexl(1,ielem_indexi(j,i),i)
        rowf=halo_offset(nf)+lf-1
        sols2(1:nvt)=solhir(rowf,nof_variables+1:nof_variables+nvt)/solhir(rowf,1)
      end if

      do k=1,dimensiona
        sols_f(1:nvt,k)=sols_f(1:nvt,k)+(oo2*(sols2(1:nvt)+sols1(1:nvt))*normal_all(k)*ielem_surf(j,i)*oov2)
      end do
    end do

    do var2=1,nvt
      rec_gradientsturb(1,1:dimensiona,var2,i)=sols_f(var2,1:dimensiona)
      rec_grads(dimensiona+1+var2,1:dimensiona,i)=sols_f(var2,1:dimensiona)
    end do
  end if
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine compute_gradients_mix_turb_ggs_viscous_acc


 
subroutine allgrads_inner(n,iconsidered)
!> @brief
!> this subroutine calls the gradient approximation subroutines for every interior cell
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
integer::i
integer::number_of_dog,number_of_nei,imax
i=iconsidered
number_of_dog=ielem_idegfree(i)
number_of_nei=ielem_inumneighbours(i)
! imax=number_of_nei-1



if (dg.eq.1)then

call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)

else



select case(ielem_ggs(i))

case(0)	!least squares everything
    select case (itestcase)
 
    case(1,2,3)
      
    call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)


    if (initcond.eq.95)then
    call compute_gradients_inner_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
    end if

      
    case(4)
	call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)
    call compute_gradients_inner_mean_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)
		if (turbulence.eq.1)then
			call compute_gradients_inner_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
			call compute_gradients_turb_lsq(n,iconsidered,number_of_dog,number_of_nei)
			call compute_gradients_turb_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)
		end if





   end select
   
   
 case(1) !green gauss everything
      select case (itestcase)
 
 
       case(1,2,3)
      
      call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)


      if (initcond.eq.95)then
    call compute_gradients_inner_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
    end if



      case(4)
      call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)
      call compute_gradients_inner_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
      
      
      if (turbulence.eq.1)then
      call compute_gradients_turb_lsq(n,iconsidered,number_of_dog,number_of_nei)
      call compute_gradients_inner_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
      end if
      



      end select
      

 end select


 end if



end subroutine allgrads_inner

subroutine allgrads_mix(n,iconsidered)
!> @brief
!> this subroutine calls the gradient approximation subroutines for every non-interior cell
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
integer::i
integer::number_of_dog,number_of_nei,imax
i=iconsidered
number_of_dog=ielem_idegfree(i)
number_of_nei=ielem_inumneighbours(i)


if (dg.eq.1)then

		call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)


else







			if (fastest.ne.1)then !least squares
			select case(ielem_ggs(i))

				case(0)	!least squares everything
				select case (itestcase)





				case(1,2,3)

				call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)

				if (initcond.eq.95)then
				call compute_gradients_mix_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
				end if




				case(4)
				call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)

								if (rec_wall(i).eq.0)then
									call compute_gradients_inner_mean_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)
								else
									call compute_gradients_wall_mean_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)
								end if

								if (turbulence.eq.1)then
									call compute_gradients_mix_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
									call compute_gradients_turb_lsq(n,iconsidered,number_of_dog,number_of_nei)

									if (rec_wall(i).eq.0)then
										call compute_gradients_turb_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)
									else
										call compute_gradients_wall_turb_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)
									end if
								end if







			end select


			case(1)
				select case (itestcase)


				case(1,2,3)

				call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)


				if (initcond.eq.95)then
				call compute_gradients_mix_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
				end if

				case(4)
				call compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)
				call compute_gradients_mix_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)


				if (turbulence.eq.1)then
				call compute_gradients_turb_lsq(n,iconsidered,number_of_dog,number_of_nei)

				call compute_gradients_mix_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)
				end if


				end select
			end select


			end if

	end if
      


end subroutine allgrads_mix




subroutine allgrads_mix_av(n,iconsidered)
!> @brief
!> this subroutine calls the average gradient approximation subroutines for every non interior cell 
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
integer::i
integer::number_of_dog,number_of_nei,imax
i=iconsidered

number_of_dog=ielem_idegfree(i)
number_of_nei=ielem_inumneighbours(i)

 call compute_gradients_mix_mean_ggs_viscous_av(n,iconsidered,number_of_dog,number_of_nei)


end subroutine allgrads_mix_av


subroutine allgrads_inner_av(n,iconsidered)
!> @brief
!> this subroutine calls the average gradient approximation subroutines for every interior cell 
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
integer::i
integer::number_of_dog,number_of_nei,imax
i=iconsidered
number_of_dog=ielem_idegfree(i)
number_of_nei=ielem_inumneighbours(i)

 call compute_gradients_inner_mean_ggs_viscous_av(n,iconsidered,number_of_dog,number_of_nei)


end subroutine allgrads_inner_av



subroutine compute_gradients_mean_lsq(n,iconsidered,number_of_dog,number_of_nei)!check all
!> @brief
!> this subroutine computes the gradients of the conserved variables of each cell using the least-squares
   implicit none
#ifdef gpu
!$omp declare target
#endif
   integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
   real,dimension(1:gpu_max_nvar)::sols1,sols2,leftv
   real::mp_pinfl,gammal
   integer::i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local
   real::tempxx,diff,coef

   i=iconsidered
   imax=number_of_nei-1
   sols1=zero
   sols2=zero

   sols1(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
   if (wenwrt.eq.3)then
      leftv(1:nof_variables)=sols1(1:nof_variables)
      call cons2prim(n,leftv,mp_pinfl,gammal)
      sols1(1:nof_variables)=leftv(1:nof_variables)
   end if

   if (rec_local(i).eq.0)then
      do ll=1,ielem_admis(i)
         if ((ees.ne.5).or.(ll.eq.1))then
            ideg_local=number_of_dog
            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradients(ll,k,var2,iconsidered)=zero
               end do
            end do

            do iq=1,imax
               sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(ll,iq+1,i))
               if (wenwrt.eq.3)then
                  leftv(1:nof_variables)=sols2(1:nof_variables)
                  call cons2prim(n,leftv,mp_pinfl,gammal)
                  sols2(1:nof_variables)=leftv(1:nof_variables)
               end if

               if(per_rot.eq.1)then
                  if (rec_periodicflag(ll,iq+1,i).eq.2) then
                     tempxx=sols2(2)
                     sols2(2)=tempxx*cos(angle_per)-sols2(3)*sin(angle_per)
                     sols2(3)=tempxx*sin(angle_per)+sols2(3)*cos(angle_per)
                  end if
               end if

               do k=1,ideg_local
                  coef=rec_invmat_stencilt(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     rec_gradients(ll,k,var2,iconsidered)=rec_gradients(ll,k,var2,iconsidered)+coef*diff
                  end do
               end do
            end do
         else
            ideg_local=idegfree2
            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradientsc(ll,k,var2,iconsidered)=zero
               end do
            end do

            do iq=1,numneighbours2-1
               sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexlc(ll,iq+1,i))
               if (wenwrt.eq.3)then
                  leftv(1:nof_variables)=sols2(1:nof_variables)
                  call cons2prim(n,leftv,mp_pinfl,gammal)
                  sols2(1:nof_variables)=leftv(1:nof_variables)
               end if

               do k=1,ideg_local
                  coef=rec_invmat_stenciltc(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     rec_gradientsc(ll,k,var2,iconsidered)=rec_gradientsc(ll,k,var2,iconsidered)+coef*diff
                  end do
               end do
            end do
         end if
      end do
   else
      do ll=1,ielem_admis(i)
         if ((ees.ne.5).or.(ll.eq.1))then
            ideg_local=number_of_dog
            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradients(ll,k,var2,iconsidered)=zero
               end do
            end do

            do iq=1,imax
               if (rec_ihexb(ll,iq+1,rec_local(i)).eq.n)then
                  sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(ll,iq+1,i))
               else
                  nf=rec_ihexn(ll,iq+1,rec_local(i))
                  lf=rec_ihexl(ll,iq+1,i)
                  rowf=halo_offset(nf) + lf - 1
                  sols2(1:nof_variables)=solhir(rowf, 1:nof_variables)
               end if

               if (wenwrt.eq.3)then
                  leftv(1:nof_variables)=sols2(1:nof_variables)
                  call cons2prim(n,leftv,mp_pinfl,gammal)
                  sols2(1:nof_variables)=leftv(1:nof_variables)
               end if

               if(per_rot.eq.1)then
                  if (rec_periodicflag(ll,iq+1,i).eq.2) then
                     tempxx=sols2(2)
                     sols2(2)=tempxx*cos(angle_per)-sols2(3)*sin(angle_per)
                     sols2(3)=tempxx*sin(angle_per)+sols2(3)*cos(angle_per)
                  end if
               end if

               do k=1,ideg_local
                  coef=rec_invmat_stencilt(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     rec_gradients(ll,k,var2,iconsidered)=rec_gradients(ll,k,var2,iconsidered)+coef*diff
                  end do
               end do
            end do
         else
            ideg_local=idegfree2
            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradientsc(ll,k,var2,iconsidered)=zero
               end do
            end do

            do iq=1,numneighbours2-1
               if (rec_ihexbc(ll,iq+1,rec_local(i)).eq.n)then
                  sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexlc(ll,iq+1,i))
               else
                  nf=rec_ihexnc(ll,iq+1,rec_local(i))
                  lf=rec_ihexlc(ll,iq+1,i)
                  rowf=halo_offset(nf) + lf - 1
                  sols2(1:nof_variables)=solhir(rowf, 1:nof_variables)
               end if

               if (wenwrt.eq.3)then
                  leftv(1:nof_variables)=sols2(1:nof_variables)
                  call cons2prim(n,leftv,mp_pinfl,gammal)
                  sols2(1:nof_variables)=leftv(1:nof_variables)
               end if

               do k=1,ideg_local
                  coef=rec_invmat_stenciltc(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     rec_gradientsc(ll,k,var2,iconsidered)=rec_gradientsc(ll,k,var2,iconsidered)+coef*diff
                  end do
               end do
            end do
         end if
      end do
   end if

end subroutine compute_gradients_mean_lsq



subroutine compute_gradients_mean_lsq_acc(n)!check all
!> @brief
!> this subroutine computes the gradients of the conserved variables of each cell using the least-squares
   implicit none
   integer,intent(in)::n
   real,dimension(1:gpu_max_nvar)::sols1,sols2
   real,dimension(1:gpu_max_dof,1:gpu_max_nvar)::gradacc
   integer::i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local,number_of_dog,iconsidered
   real::tempxx,diff,coef
#ifdef xpu


   !$omp target teams distribute parallel do &
   !$omp& firstprivate(n) &
   !$omp& firstprivate(numneighbours2, idegfree2, nof_variables, ees, per_rot, angle_per, zero) &
   !$omp& map(alloc: xmpielrank, ielem_inumneighbours, ielem_idegfree, u_c_val) &
   !$omp& map(alloc: rec_ihexl, rec_local, ielem_admis, rec_periodicflag) &
   !$omp& map(alloc: rec_invmat_stencilt, rec_gradients, rec_ihexlc, rec_invmat_stenciltc) &
   !$omp& map(alloc: rec_gradientsc, rec_ihexb, rec_ihexn, halo_offset, solhir, ielem_xxc) &
   !$omp& map(alloc: rec_ihexbc, rec_ihexnc) &
   !$omp& private(i,var2,iconsidered,iq,ll,imax,nf,lf,rowf,k,ideg_local,number_of_dog,tempxx,diff,coef,sols1,sols2,gradacc)
#else
   !$omp do
#endif

   do i=1,xmpielrank(n)
   iconsidered=i
   imax=ielem_inumneighbours(i)-1
   number_of_dog=ielem_idegfree(i)
   sols1=zero
   sols2=zero

   sols1(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))


   if (rec_local(i).eq.0)then
      do ll=1,ielem_admis(i)
         if ((ees.ne.5).or.(ll.eq.1))then
            ideg_local=number_of_dog
            do var2=1,nof_variables
               do k=1,ideg_local
                  gradacc(k,var2)=zero
               end do
            end do

            do iq=1,imax
               sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(ll,iq+1,i))


               if(per_rot.eq.1)then
                  if (rec_periodicflag(ll,iq+1,i).eq.2) then
                     tempxx=sols2(2)
                     sols2(2)=tempxx*cos(angle_per)-sols2(3)*sin(angle_per)
                     sols2(3)=tempxx*sin(angle_per)+sols2(3)*cos(angle_per)
                  end if
               end if

               do k=1,ideg_local
                  coef=rec_invmat_stencilt(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     gradacc(k,var2)=gradacc(k,var2)+coef*diff
                  end do
               end do
            end do

            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradients(ll,k,var2,iconsidered)=gradacc(k,var2)
               end do
            end do
         else
            ideg_local=idegfree2
            do var2=1,nof_variables
               do k=1,ideg_local
                  gradacc(k,var2)=zero
               end do
            end do

            do iq=1,numneighbours2-1
               sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexlc(ll,iq+1,i))


               do k=1,ideg_local
                  coef=rec_invmat_stenciltc(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     gradacc(k,var2)=gradacc(k,var2)+coef*diff
                  end do
               end do
            end do

            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradientsc(ll,k,var2,iconsidered)=gradacc(k,var2)
               end do
            end do
         end if
      end do
   else
      do ll=1,ielem_admis(i)
         if ((ees.ne.5).or.(ll.eq.1))then
            ideg_local=number_of_dog
            do var2=1,nof_variables
               do k=1,ideg_local
                  gradacc(k,var2)=zero
               end do
            end do

            do iq=1,imax
               if (rec_ihexb(ll,iq+1,rec_local(i)).eq.n)then
                  sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(ll,iq+1,i))
               else
                  nf=rec_ihexn(ll,iq+1,rec_local(i))
                  lf=rec_ihexl(ll,iq+1,i)
                  rowf=halo_offset(nf) + lf - 1
                  sols2(1:nof_variables)=solhir(rowf, 1:nof_variables)
               end if



               if(per_rot.eq.1)then
                  if (rec_periodicflag(ll,iq+1,i).eq.2) then
                     tempxx=sols2(2)
                     sols2(2)=tempxx*cos(angle_per)-sols2(3)*sin(angle_per)
                     sols2(3)=tempxx*sin(angle_per)+sols2(3)*cos(angle_per)
                  end if
               end if

               do k=1,ideg_local
                  coef=rec_invmat_stencilt(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     gradacc(k,var2)=gradacc(k,var2)+coef*diff
                  end do
               end do
            end do

            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradients(ll,k,var2,iconsidered)=gradacc(k,var2)
               end do
            end do
         else
            ideg_local=idegfree2
            do var2=1,nof_variables
               do k=1,ideg_local
                  gradacc(k,var2)=zero
               end do
            end do

            do iq=1,numneighbours2-1
               if (rec_ihexbc(ll,iq+1,rec_local(i)).eq.n)then
                  sols2(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexlc(ll,iq+1,i))
               else
                  nf=rec_ihexnc(ll,iq+1,rec_local(i))
                  lf=rec_ihexlc(ll,iq+1,i)
                  rowf=halo_offset(nf) + lf - 1
                  sols2(1:nof_variables)=solhir(rowf, 1:nof_variables)
               end if



               do k=1,ideg_local
                  coef=rec_invmat_stenciltc(k,iq,ll,i)
                  do var2=1,nof_variables
                     diff=sols2(var2)-sols1(var2)
                     gradacc(k,var2)=gradacc(k,var2)+coef*diff
                  end do
               end do
            end do

            do var2=1,nof_variables
               do k=1,ideg_local
                  rec_gradientsc(ll,k,var2,iconsidered)=gradacc(k,var2)
               end do
            end do
         end if
      end do
   end if
   end do
#ifdef xpu
   !$omp end target teams distribute parallel do

!    if (allocated(rec_gradients)) then
!    end if
!    if (allocated(rec_gradientsc)) then
!    end if
#else
   !$omp end do
#endif

end subroutine compute_gradients_mean_lsq_acc








subroutine compute_gradients_inner_mean_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitve variables of each interior cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:gpu_max_nvar)::sols1,sols2,leftv
integer::i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local
real::mp_pinfl,gammal,diff,coef

imax=number_of_nei-1
i=iconsidered
ll=1
ideg_local=number_of_dog
sols1=zero
sols2=zero

leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
call cons2div(n,leftv,mp_pinfl,gammal)
sols1(1:nof_variables-1)=leftv(2:nof_variables)

do var2=1,nof_variables-1
   do k=1,ideg_local
      rec_gradf(var2,k,iconsidered)=zero
   end do
end do

do iq=1,imax
   if (rec_local(i).eq.0)then
      leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
   else
      if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
         leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
      else
         nf=rec_ihexn(1,iq+1,rec_local(i))
         lf=rec_ihexl(1,iq+1,i)
         rowf=halo_offset(nf) + lf - 1
         leftv(1:nof_variables)=solhir(rowf,1:nof_variables)
      end if
   end if

   call cons2div(n,leftv,mp_pinfl,gammal)
   sols2(1:nof_variables-1)=leftv(2:nof_variables)

   do k=1,ideg_local
      coef=rec_invmat_stencilt(k,iq,ll,i)
      do var2=1,nof_variables-1
         diff=sols2(var2)-sols1(var2)
         rec_gradf(var2,k,iconsidered)=rec_gradf(var2,k,iconsidered)+coef*diff
      end do
   end do
end do

end subroutine compute_gradients_inner_mean_lsq_viscous





subroutine compute_gradients_inner_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(gpu_max_extra_transport)::sols1,sols2
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::sols_f
real,dimension(gpu_max_dim)::normal_all
real::oov2,angle1,angle2
integer::i,j,k,l,var2
real,dimension(1:gpu_max_nvar)::leftv







i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)


if (dimensiona.eq.3)then


	  sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)/u_c_val(1,1,i)

do j=1,ielem_ifca(i)
			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=(cos(angle1)*sin(angle2))
				normal_all(2)=(sin(angle1)*sin(angle2))
				normal_all(3)=(cos(angle2))

			sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(j,i))/u_c_val(1,1,ielem_ineigh(j,i))


			do k=1,3
			sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do
end do

			  do var2=1,turbulenceequations+passivescalar
			    rec_gradientsturb(1,1:3,var2,iconsidered)=sols_f(var2,1:3)
			    rec_grads(4+var2,1:3,i)=sols_f(var2,1:3)
			 end do

else


 sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)/u_c_val(1,1,i)

do j=1,ielem_ifca(i)
			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=angle1
				normal_all(2)=angle2

			sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(j,i))/u_c_val(1,1,ielem_ineigh(j,i))


			do k=1,2
			sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do
end do

			  do var2=1,turbulenceequations+passivescalar
			    rec_gradientsturb(1,1:2,var2,iconsidered)=sols_f(var2,1:2)
			    rec_grads(3+var2,1:2,i)=sols_f(var2,1:2)
			 end do







end if




end subroutine compute_gradients_inner_turb_ggs_viscous





subroutine compute_gradients_turb_lsq(n,iconsidered,number_of_dog,number_of_nei) !check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(gpu_max_extra_transport)::sols1,sols2
integer::i,var2,ll,iq,nf,lf,rowf,k,imax,ideg_local,nvt
real::diff,coef

imax=number_of_nei-1
i=iconsidered
nvt=turbulenceequations+passivescalar
sols1=zero
sols2=zero

sols1(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,1,i))

if (rec_local(i).eq.0)then
   do ll=1,ielem_admis(i)
      if ((ees.ne.5).or.(ll.eq.1))then
         ideg_local=number_of_dog
         do var2=1,nvt
            do k=1,ideg_local
               rec_gradients2(ll,k,var2,iconsidered)=zero
            end do
         end do

         do iq=1,imax
            sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(ll,iq+1,i))
            do k=1,ideg_local
               coef=rec_invmat_stencilt(k,iq,ll,i)
               do var2=1,nvt
                  diff=sols2(var2)-sols1(var2)
                  rec_gradients2(ll,k,var2,iconsidered)=rec_gradients2(ll,k,var2,iconsidered)+coef*diff
               end do
            end do
         end do
      else
         ideg_local=idegfree2
         do var2=1,nvt
            do k=1,ideg_local
               rec_gradientsc2(ll,k,var2,iconsidered)=zero
            end do
         end do

         do iq=1,numneighbours2-1
            sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexlc(ll,iq+1,i))
            do k=1,ideg_local
               coef=rec_invmat_stenciltc(k,iq,ll,i)
               do var2=1,nvt
                  diff=sols2(var2)-sols1(var2)
                  rec_gradientsc2(ll,k,var2,iconsidered)=rec_gradientsc2(ll,k,var2,iconsidered)+coef*diff
               end do
            end do
         end do
      end if
   end do
else
   do ll=1,ielem_admis(i)
      if ((ees.ne.5).or.(ll.eq.1))then
         ideg_local=number_of_dog
         do var2=1,nvt
            do k=1,ideg_local
               rec_gradients2(ll,k,var2,iconsidered)=zero
            end do
         end do

         do iq=1,imax
            if (rec_ihexb(ll,iq+1,rec_local(i)).eq.n)then
               sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(ll,iq+1,i))
            else
               nf=rec_ihexn(ll,iq+1,rec_local(i))
               lf=rec_ihexl(ll,iq+1,i)
               rowf=halo_offset(nf) + lf - 1
               sols2(1:nvt)=solhir(rowf, nof_variables+1:nof_variables+nvt)
            end if

            do k=1,ideg_local
               coef=rec_invmat_stencilt(k,iq,ll,i)
               do var2=1,nvt
                  diff=sols2(var2)-sols1(var2)
                  rec_gradients2(ll,k,var2,iconsidered)=rec_gradients2(ll,k,var2,iconsidered)+coef*diff
               end do
            end do
         end do
      else
         ideg_local=idegfree2
         do var2=1,nvt
            do k=1,ideg_local
               rec_gradientsc2(ll,k,var2,iconsidered)=zero
            end do
         end do

         do iq=1,numneighbours2-1
            if (rec_ihexbc(ll,iq+1,rec_local(i)).eq.n)then
               sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexlc(ll,iq+1,i))
            else
               nf=rec_ihexnc(ll,iq+1,rec_local(i))
               lf=rec_ihexlc(ll,iq+1,i)
               rowf=halo_offset(nf) + lf - 1
               sols2(1:nvt)=solhir(rowf, nof_variables+1:nof_variables+nvt)
            end if

            do k=1,ideg_local
               coef=rec_invmat_stenciltc(k,iq,ll,i)
               do var2=1,nvt
                  diff=sols2(var2)-sols1(var2)
                  rec_gradientsc2(ll,k,var2,iconsidered)=rec_gradientsc2(ll,k,var2,iconsidered)+coef*diff
               end do
            end do
         end do
      end if
   end do
end if

end subroutine compute_gradients_turb_lsq



subroutine compute_gradients_turb_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each interior cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(gpu_max_extra_transport)::sols1,sols2
integer::i,var2,iq,ll,imax,nf,lf,rowf,k,ideg_local,nvt
real::diff,coef

imax=number_of_nei-1
i=iconsidered
ll=1
ideg_local=number_of_dog
nvt=turbulenceequations+passivescalar
sols1=zero
sols2=zero

sols1(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,1,i))/u_c_val(1,1,rec_ihexl(1,1,i))

do var2=1,nvt
   do k=1,ideg_local
      rec_gradientsturb(1,k,var2,iconsidered)=zero
   end do
end do

do iq=1,imax
   if (rec_local(i).eq.0)then
      sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,iq+1,i))/u_c_val(1,1,rec_ihexl(1,iq+1,i))
   else
      if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
         sols2(1:nvt)=u_ct_val(1,1:nvt,rec_ihexl(1,iq+1,i))/u_c_val(1,1,rec_ihexl(1,iq+1,i))
      else
         nf=rec_ihexn(ll,iq+1,rec_local(i))
         lf=rec_ihexl(ll,iq+1,i)
         rowf=halo_offset(nf) + lf - 1
         sols2(1:nvt)=solhir(rowf, nof_variables+1:nof_variables+nvt)/solhir(rowf, 1)
      end if
   end if

   do k=1,ideg_local
      coef=rec_invmat_stencilt(k,iq,ll,i)
      do var2=1,nvt
         diff=sols2(var2)-sols1(var2)
         rec_gradientsturb(1,k,var2,iconsidered)=rec_gradientsturb(1,k,var2,iconsidered)+coef*diff
      end do
   end do
end do

end subroutine compute_gradients_turb_lsq_viscous


subroutine compute_gradients_inner_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:gpu_max_nvar)::sols1,sols2,dudl,aver1,phi_f
real,dimension(1:gpu_max_nvar,gpu_max_dim)::sols_f
real,dimension(3)::normal_all,dih_vec,e_ih,sf
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,aorth,dih
integer::i,j,k,l
real,dimension(1:gpu_max_nvar)::leftv


i=iconsidered
sols_f=zero;sols1=zero;sols2=zero

rec_grads(:,:,i)=zero

oov2=1.0d0/ielem_totvolume(i)


if (dimensiona.eq.3)then


	    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
	    call cons2div(n,leftv,mp_pinfl,gammal)
	  sols1(1:nof_variables-1)=leftv(2:nof_variables)


do j=1,ielem_ifca(i)
			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=(cos(angle1)*sin(angle2))
				normal_all(2)=(sin(angle1)*sin(angle2))
				normal_all(3)=(cos(angle2))

			leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
			call cons2div(n,leftv,mp_pinfl,gammal)
			sols2(1:nof_variables-1)=leftv(2:nof_variables)

 			do k=1,dimensiona
 			sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem_surf(j,i)*oov2)
 			end do
end do

			do k=1,dimensiona
			rec_grads(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
			end do

else

				    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
	    call cons2div(n,leftv,mp_pinfl,gammal)
	  sols1(1:nof_variables-1)=leftv(2:nof_variables)


do j=1,ielem_ifca(i)
			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=angle1
				normal_all(2)=angle2


				dih_vec(1:dimensiona)=ielem_dih2(j,1:dimensiona,i)
				dih=ielem_dih(j,i)
				e_ih(1:dimensiona)=dih_vec(1:dimensiona)/dih

			leftv(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
			call cons2div(n,leftv,mp_pinfl,gammal)
			sols2(1:nof_variables-1)=leftv(2:nof_variables)

			do k=1,2
			sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do



! 					! build face area vector
! 					do k = 1, dimensiona
! 					sf(k) = normal_all(k) * ielem_surf(j,i)
! 					end do
!
! 					! orthogonal projected area along centroid-to-centroid line
! 					aorth = 0.0d0
! 					do k = 1, dimensiona
! 					aorth = aorth + sf(k) * e_ih(k)
! 					end do
!
! 					! face value: still simple average here
! 					phi_f(1:nof_variables) = oo2*(sols1(1:nof_variables) + sols2(1:nof_variables))
!
! 					! accumulate orthogonal gg contribution
! 					do k = 1, dimensiona
! 					sols_f(1:nof_variables,k) = sols_f(1:nof_variables,k) + &
! 						phi_f(1:nof_variables) * aorth * e_ih(k) * oov2
! 					end do





end do

			do k=1,dimensiona
			rec_grads(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
			end do





end if





end subroutine compute_gradients_inner_mean_ggs_viscous


subroutine compute_gradients_wall_mean_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)!check all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each non-interior cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(gpu_max_nvar)::sols1,sols2
real,dimension(1:gpu_max_nvar,1:gpu_max_neighbours)::matrix_1
real,dimension(1:gpu_max_nvar,1:gpu_max_dof)::matrix_2
real,dimension(1:gpu_max_dof,1:gpu_max_nvar)::sol_m
real,dimension(1:gpu_max_nvar)::matrix_3
integer::i,var2,ii,k0,g0,ttk,ivvm,iq,lq,irg
real::attt
integer::ll,imax,nf,lf,rowf
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::leftv






imax=number_of_nei-1








ll=1
i=iconsidered
sols1=zero;
sols2=zero

	    ll=1
	    k0=rec_k0(rec_wall(i))
	    g0=rec_g0(rec_wall(i))

	     matrix_1=zero;matrix_2=zero;sol_m=zero;
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
		call cons2div(n,leftv,mp_pinfl,gammal)

	       sols1(1:nof_variables-1)=leftv(2:nof_variables)
! 	       sols1(1)=leftv(5)/(leftv(1)*r_gas)

	      do iq=1,imax
			if (rec_local(i).eq.0)then
			leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
			else
				if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
				leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))
				else
! 				leftv(1:nof_variables)=iexsolhir(rec_ihexn(1,iq+1,i))%sol(rec_ihexl(1,iq+1,i),1:nof_variables)
				    nf=rec_ihexn(1,iq+1,rec_local(i))
					lf=rec_ihexl(1,iq+1,i)
					rowf=halo_offset(nf) + lf - 1
					leftv(1:nof_variables)=solhir(rowf, 1:nof_variables)



				end if
			end if

		 call cons2div(n,leftv,mp_pinfl,gammal)
	       sols2(1:nof_variables-1)=leftv(2:nof_variables)
! 	       !sols2(1)=leftv(5)/(leftv(1)*r_gas)
  	        matrix_1(1:nof_variables-1,iq)=(rec_volume_w(1,iq+1,rec_wall(i))*rec_weightl(1,iq,rec_wall(i))*(sols2(1:nof_variables-1)-sols1(1:nof_variables-1)))

  	        !velocity gradients
  	        matrix_1(1:dimensiona,iq)=matrix_1(1:dimensiona,iq)+((sols1(1:dimensiona)*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))

  	        !temperature now check

			if (thermal.eq.1)then
  	        matrix_1(dimensiona+1:nof_variables-nof_species-1,iq)=matrix_1(dimensiona+1:nof_variables-nof_species-1,iq)+((sols1(dimensiona+1:nof_variables-nof_species-1)*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))-(((wall_temp)*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))

  	        end if


   	        if (catalytic_wall.eq.1)then
   	        matrix_1(dimensiona+3:nof_variables-1,iq)=matrix_1(dimensiona+3:nof_variables-1,iq)+((sols1(dimensiona+3:nof_variables-1)*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))-(((catalytic_con(1:nof_species))*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))



   	        end if



		end do
		matrix_3(:)=zero
		matrix_3(1:dimensiona)=-sols1(1:dimensiona)


		do var2=1,nof_variables-1
		  matrix_2(var2,1:number_of_dog-1)=zero
		  if (var2.le.dimensiona)then		!velocity gradients
		  do iq=1,imax

		      do lq=1,number_of_dog-1
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
		      end do

		  end do
		  end if
		  if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then	!temperature gradients
		  do iq=1,imax
		     do lq=1,number_of_dog-1
		     if (thermal.eq.1)then
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
			else
			matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_tempsq(iq,lq,rec_wall(i))
			end if
		      end do
		  end do
		  end if
		   if (var2.gt.nof_variables-nof_species-1)then					!species
		   do iq=1,imax
		     do lq=1,number_of_dog-1

		     if (catalytic_wall.eq.1)then
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
		     else
			matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_tempsq(iq,lq,rec_wall(i))
			end if

		      end do
		  end do
		   end if

		if (var2.le.dimensiona)then
		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		end if
		if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then
		if (thermal.eq.1)then
		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		else
		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_tempsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		end if
		end if

		if (var2.gt.nof_variables-nof_species-1)then

		if (catalytic_wall.eq.1)then
			do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do

		else

		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_tempsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do
		end if
		end if

	     end do

		do var2=1,nof_variables-1

		 if (var2.le.dimensiona)then		!velocity gradients
		 rec_gradf(var2,1:idegfree,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		  end do
		  attt=zero
		  attt=-sols1(var2)
			  do ttk=1,number_of_dog
				    if (ttk.ne.k0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradf(var2,k0,iconsidered)=attt


		end if

		if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then	!temperature gradients
		if (thermal.eq.1)then
		rec_gradf(var2,1:idegfree,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		  end do
		  attt=zero
		  attt=wall_temp-sols1(var2)
			  do ttk=1,number_of_dog
				    if (ttk.ne.k0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradf(var2,k0,iconsidered)=attt

		else
		ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.g0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		    end do
		    attt=zero

			  do ttk=1,number_of_dog
				    if (ttk.ne.g0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoefg(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoefg(g0,rec_wall(i))
			    rec_gradf(var2,g0,iconsidered)=attt
		end if
		end if
		if (var2.gt.nof_variables-nof_species-1)then				!species gradients

			if (catalytic_wall.eq.1)then
		rec_gradf(var2,1:idegfree,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
		  end do
		  attt=zero
		  attt=catalytic_con(var2-dimensiona-2)-sols1(var2)


			  do ttk=1,number_of_dog
				    if (ttk.ne.k0) &
				  attt=attt-rec_gradf(var2,ttk,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradf(var2,k0,iconsidered)=attt

		else



			ivvm=0
				do ttk=1,number_of_dog
						if (ttk.eq.g0) cycle
						ivvm=ivvm+1
							rec_gradf(var2,ttk,iconsidered)=sol_m(ivvm,var2)
				end do
				attt=zero

				do ttk=1,number_of_dog
						if (ttk.ne.g0) &
					attt=attt-rec_gradf(var2,ttk,iconsidered)*&
								rec_wallcoefg(ttk,rec_wall(i))
				end do
					attt=attt/rec_wallcoefg(g0,rec_wall(i))
					rec_gradf(var2,g0,iconsidered)=attt

		end if

		end if


		end do













end subroutine compute_gradients_wall_mean_lsq_viscous





subroutine compute_gradients_mix_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulnece variables of each non-interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(gpu_max_extra_transport)::sols1,sols2
real,dimension(gpu_max_extra_transport,3)::sols_f
real,dimension(3)::normal_all,temp_vert
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,nx,ny,nz
integer::i,j,k,l,var2,b_code,facex,n_node,imax,nf,lf,rowf
real,dimension(1:gpu_max_nvar)::leftv,srf_speed,srf_speedrot,rightv
real,dimension(1:gpu_max_dim)::pox,poy,poz,cords
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_nvar_total)::cright_rot,cleft_rot
integer::ibfc



i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)



if (dimensiona.eq.3)then



	  sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)/u_c_val(1,1,i)


do j=1,ielem_ifca(i)
			 facex=j


			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=(cos(angle1)*sin(angle2))
				normal_all(2)=(sin(angle1)*sin(angle2))
				normal_all(3)=(cos(angle2))
				nx=normal_all(1);ny=normal_all(2);nz=normal_all(3)


			if (ielem_ineighb(j,i).eq.n)then	!my cpu only
			    if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in my cpu
				  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(j,i))/&
				  u_c_val(1,1,ielem_ineigh(j,i))
				  else
				  !not periodic ones in my cpu

				  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

				   if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if


				  call cordinates3(n,nodes_list,n_node,cords(1:3))
				  pox(1)=cords(1);poy(1)=cords(2);poz(1)=cords(3)





				  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
				  cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
				  b_code=ibound_icode(ielem_ibounds(j,i))
				  call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

				  sols2(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)/rightv(1)

				  end if
			    else
				    sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(j,i))/&
				    u_c_val(1,1,ielem_ineigh(j,i))




			    end if
			else	!in other cpus they can only be periodic or mpi neighbours

			      if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in other cpu

					!sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
					!(rec_ihexl(1,ielem_indexi(j,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
					!iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
					!(rec_ihexl(1,ielem_indexi(j,i)),1)

					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)/solhir(rowf,1)



				  end if
			      else


	!				sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
	!				(rec_ihexl(1,ielem_indexi(j,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
	!				iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
	!				(rec_ihexl(1,ielem_indexi(j,i)),1)


					nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)/solhir(rowf,1)




			     end if
			end if

			do k=1,3
			sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do
end do

					 do var2=1,turbulenceequations+passivescalar
			    rec_gradientsturb(1,1:3,var2,iconsidered)=sols_f(var2,1:3)
			    rec_grads(4+var2,1:3,i)=sols_f(var2,1:3)
			 end do


		else	!2d

		i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)

	  sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)/u_c_val(1,1,i)


do j=1,ielem_ifca(i)
			 facex=j


			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=angle1
				normal_all(2)=angle2
				nx=normal_all(1);ny=normal_all(2)


			if (ielem_ineighb(j,i).eq.n)then	!my cpu only
			    if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in my cpu
				  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(j,i))/&
				  u_c_val(1,1,ielem_ineigh(j,i))
				  else
				  !not periodic ones in my cpu

				  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
				  n_node=2
				  call cordinates2(n,nodes_list,n_node,cords(1:2))
				  pox(1)=cords(1);poy(1)=cords(2);

				    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)




				  cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
				  b_code=ibound_icode(ielem_ibounds(j,i))
				  call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

				  sols2(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)/rightv(1)

				  end if
			    else
				    sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(j,i))/&
				  u_c_val(1,1,ielem_ineigh(j,i))




			    end if
			else	!in other cpus they can only be periodic or mpi neighbours

			      if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in other cpu

! 					sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
! 					iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),1)

					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)/solhir(rowf,1)




				  end if
			      else


! 					sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
! 					iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),1)


					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)/solhir(rowf,1)


			     end if
			end if

			do k=1,2
			sols_f(1:turbulenceequations+passivescalar,k)=sols_f(1:turbulenceequations+passivescalar,k)+((oo2*(sols2(1:turbulenceequations+passivescalar)+sols1(1:turbulenceequations+passivescalar)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do
end do

					 do var2=1,turbulenceequations+passivescalar
			    rec_gradientsturb(1,1:2,var2,iconsidered)=sols_f(var2,1:2)
			    rec_grads(3+var2,1:2,i)=sols_f(var2,1:2)
			 end do



		end if















end subroutine compute_gradients_mix_turb_ggs_viscous





subroutine compute_gradients_wall_turb_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei) !check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each non-interior cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:gpu_max_extra_transport)::sols1,sols2
real,dimension(1:gpu_max_extra_transport,1:gpu_max_neighbours)::matrix_1
real,dimension(1:gpu_max_extra_transport,1:gpu_max_dof)::matrix_2
real,dimension(1:gpu_max_dof,1:gpu_max_extra_transport)::sol_m
real,dimension(1:gpu_max_extra_transport)::matrix_3
integer::i,var2,ii,k0,g0,ttk,ivvm,iq,lq,imax
real::attt
integer::ll,nf,lf,rowf






imax=number_of_nei-1

ll=1

i=iconsidered
sols1=zero;
sols2=zero

	    k0=rec_k0(rec_wall(i))



	     matrix_1=zero;matrix_2=zero;sol_m=zero;



		sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,1,i))/&
		u_c_val(1,1,rec_ihexl(1,1,i))








	      do iq=1,imax
	      if (rec_local(i).eq.0)then

	      sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,iq+1,i))/&
	      u_c_val(1,1,rec_ihexl(1,iq+1,i))
	      else

		 if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
		sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,iq+1,i))/&
		u_c_val(1,1,rec_ihexl(1,iq+1,i))
	    else

! 		sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,iq+1,i))%sol(rec_ihexl(1,iq+1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
! 		iexsolhir(rec_ihexn(1,iq+1,i))%sol(rec_ihexl(1,iq+1,i),1)


					nf=rec_ihexn(1,iq+1,rec_local(i))
					lf=rec_ihexl(1,iq+1,i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)/solhir(rowf,1)






	    end if
	      end if


  	        matrix_1(1:turbulenceequations+passivescalar,iq)=(rec_volume_w(1,iq+1,rec_wall(i))*rec_weightl(1,iq,rec_wall(i))*(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar)))
  	        matrix_1(1:turbulenceequations+passivescalar,iq)=matrix_1(1:turbulenceequations+passivescalar,iq)+((sols1(1:turbulenceequations+passivescalar)*rec_stencils(ll,iq,k0,rec_wall(i)))/rec_wallcoeff(k0,rec_wall(i)))
		end do
		matrix_3(1:turbulenceequations+passivescalar)=-sols1(1:turbulenceequations+passivescalar)

		if (turbulencemodel.eq.2)then
		matrix_3(2)=60.0d0*visc/(beta_i1*(ielem_walldist(iconsidered)**2))
		end if

		do var2=1,turbulenceequations+passivescalar
		  matrix_2(var2,1:number_of_dog-1)=zero

		  do iq=1,imax

		      do lq=1,number_of_dog-1
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
		      end do

		  end do


		do ivvm=1,number_of_dog-1
			sol_m(ivvm,var2)=zero
			do lq=1,number_of_dog-1
			sol_m(ivvm,var2)=sol_m(ivvm,var2)+rec_velinvlsqmat(ivvm,lq,rec_wall(i))*matrix_2(var2,lq)
			end do
			end do



	     end do

		do var2=1,turbulenceequations+passivescalar


		 rec_gradientsturb(1,1:number_of_dog,var2,iconsidered)=-tolbig
		    ivvm=0
		    do ttk=1,number_of_dog
				    if (ttk.eq.k0) cycle
					  ivvm=ivvm+1
					    rec_gradientsturb(1,ttk,var2,iconsidered)=sol_m(ivvm,var2)
		  end do
			  attt=zero
		  attt=matrix_3(var2)
				  do ttk=1,number_of_dog
					    if (ttk.ne.k0) &
					  attt=attt-rec_gradientsturb(1,ttk,var2,iconsidered)*&
						    rec_wallcoeff(ttk,rec_wall(i))
			  end do
			    attt=attt/rec_wallcoeff(k0,rec_wall(i))
			    rec_gradientsturb(1,k0,var2,iconsidered)=attt

		end do











end subroutine compute_gradients_wall_turb_lsq_viscous

subroutine compute_gradients_mix_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each non-interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:gpu_max_nvar)::sols1,sols2,dudl,aver1
real,dimension(1:gpu_max_nvar,3)::sols_f
real,dimension(3)::normal_all,dih_vec,e_ih,sf
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,nx,ny,nz,aorth,dih
integer::i,j,k,l,b_code,facex,n_node,imax,nf,lf,rowf
real,dimension(1:gpu_max_nvar)::leftv,srf_speed,srf_speedrot,rightv,phi_f
real,dimension(1:gpu_max_dim)::pox,poy,poz,cords
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_nvar_total)::cright_rot,cleft_rot
integer::ibfc



i=iconsidered
sols_f=zero;sols1=zero;sols2=zero
oov2=1.0d0/ielem_totvolume(i)

rec_grads(:,:,i)=zero


if (dimensiona.eq.3)then





	  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
	    call cons2div(n,leftv,mp_pinfl,gammal)
	  sols1(1:nof_variables-1)=leftv(2:nof_variables)


	  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)




do j=1,ielem_ifca(i)
			 facex=j
			 b_code=0

			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=(cos(angle1)*sin(angle2))
				normal_all(2)=(sin(angle1)*sin(angle2))
				normal_all(3)=(cos(angle2))
				nx=normal_all(1);ny=normal_all(2);nz=normal_all(3)

				dih_vec(1:dimensiona)=ielem_dih2(j,1:dimensiona,i)
				dih=ielem_dih(j,i)
				e_ih(1:dimensiona)=dih_vec(1:dimensiona)/dih


			if (rec_mrf(iconsidered).eq.1)then
			!retrieve rotational velocity in case of rotating reference frame to calculate
			!the correct value of the boundary condition
				srf_speed(2:4)=rec_rotvel(j,1,1:3,i)
				call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
			end	if


			if (ielem_ineighb(j,i).eq.n)then	!my cpu only
			    if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if  ((ibound_icode(ielem_ibounds(j,i)).eq.5).or.(ibound_icode(ielem_ibounds(j,i)).eq.50))then	!periodic in my cpu
				  sols2(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
					  if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(j,i)).eq.50))then
	                    sols2(2:4)=rotate_per_1(sols2(2:4),ibound_icode(ielem_ibounds(j,i)),angle_per)
					  end if
				  else
				  !not periodic ones in my cpu

				  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)
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

 				  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
				  b_code=ibound_icode(ielem_ibounds(j,i))
 				  call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

				  sols2(1:nof_variables)=rightv(1:nof_variables)

				  end if
			    else
				    sols2(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))




			    end if
			else	!in other cpus they can only be periodic or mpi neighbours

					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)

				      if (ielem_ibounds(j,i).gt.0)then	!check for periodic boundaries
					  if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(j,i)).eq.50))then
		                    sols2(2:4)=rotate_per_1(sols2(2:4),ibound_icode(ielem_ibounds(j,i)),angle_per)
					  end if
				      end if
			end if

			  leftv(1:nof_variables)=sols2(1:nof_variables)
			call cons2div(n,leftv,mp_pinfl,gammal)
			sols2(1:nof_variables-1)=leftv(2:nof_variables)

			if ((b_code.eq.4).and.(thermal.eq.1))then
				sols2(dimensiona+1:nof_variables-nof_species-1)=wall_temp
			end if

			if ((b_code.eq.4).and.(catalytic_wall.eq.1))then
				sols2(dimensiona+3:nof_variables-1)=catalytic_con(1:nof_species)
			end if



!
			do k=1,dimensiona
			sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do

! 			! build face area vector
! 					do k = 1, dimensiona
! 					sf(k) = normal_all(k) * ielem_surf(j,i)
! 					end do
!
! 					! orthogonal projected area along centroid-to-centroid line
! 					aorth = 0.0d0
! 					do k = 1, dimensiona
! 					aorth = aorth + sf(k) * e_ih(k)
! 					end do
!
! 					! face value: still simple average here
! 					phi_f(1:nof_variables) = oo2*(sols1(1:nof_variables) + sols2(1:nof_variables))
!
! 					! accumulate orthogonal gg contribution
! 					do k = 1, dimensiona
! 					sols_f(1:nof_variables,k) = sols_f(1:nof_variables,k) + &
! 						phi_f(1:nof_variables) * aorth * e_ih(k) * oov2
! 					end do
end do



			do k=1,dimensiona
			rec_grads(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
			end do



	else




	  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
	    call cons2div(n,leftv,mp_pinfl,gammal)
			sols1(1:nof_variables-1)=leftv(2:nof_variables)


! 	  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)




do j=1,ielem_ifca(i)




			 facex=j

			 b_code=0

			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=angle1
				normal_all(2)=angle2
				nx=normal_all(1);ny=normal_all(2)


				dih_vec(1:dimensiona)=ielem_dih2(j,1:dimensiona,i)
				dih=ielem_dih(j,i)
				e_ih(1:dimensiona)=dih_vec(1:dimensiona)/dih


			if (ielem_ineighb(j,i).eq.n)then	!my cpu only
			    if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in my cpu
				  sols2(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))
				  else
				  !not periodic ones in my cpu




				  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
				  n_node=2
				  call cordinates2(n,nodes_list,n_node,cords(1:2))
				  pox(1)=cords(1);poy(1)=cords(2)


				  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
				  b_code=ibound_icode(ielem_ibounds(j,i))

				  call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

				  sols2(1:nof_variables)=rightv(1:nof_variables)


				  end if
			    else
				    sols2(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(j,i))

			    end if
			else	!in other cpus they can only be periodic or mpi neighbours

					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)
			end if

			  leftv(1:nof_variables)=sols2(1:nof_variables)
			call cons2div(n,leftv,mp_pinfl,gammal)
			sols2(1:nof_variables-1)=leftv(2:nof_variables)

			if ((b_code.eq.4).and.(thermal.eq.1))then
				sols2(dimensiona+1:nof_variables-nof_species-1)=wall_temp
			end if

			if ((b_code.eq.4).and.(catalytic_wall.eq.1))then
				sols2(dimensiona+3:nof_variables-1)=catalytic_con(1:nof_species)
			end if


 			do k=1,dimensiona
 			sols_f(1:nof_variables-1,k)=sols_f(1:nof_variables-1,k)+((oo2*(sols2(1:nof_variables-1)+sols1(1:nof_variables-1)))*normal_all(k)*ielem_surf(j,i)*oov2)

 			end do


! 					! build face area vector
! 					do k = 1, dimensiona
! 					sf(k) = normal_all(k) * ielem_surf(j,i)
! 					end do
!
! 					! orthogonal projected area along centroid-to-centroid line
! 					aorth = 0.0d0
! 					do k = 1, dimensiona
! 					aorth = aorth + sf(k) * e_ih(k)
! 					end do
!
! 					! face value: still simple average here
! 					phi_f(1:nof_variables) = oo2*(sols1(1:nof_variables) + sols2(1:nof_variables))
!
! 					! accumulate orthogonal gg contribution
! 					do k = 1, dimensiona
! 					sols_f(1:nof_variables,k) = sols_f(1:nof_variables,k) + &
! 						phi_f(1:nof_variables) * aorth * e_ih(k) * oov2
! 					end do









end do



			do k=1,dimensiona
			rec_grads(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
			end do

	end if












end subroutine compute_gradients_mix_mean_ggs_viscous


subroutine compute_gradients_center(n,iconsidered)
!> @brief
!> this subroutine computes the gradients of the primitive variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,dimension(1:gpu_max_nvar)::sols1,sols2,dudl,aver1
real::oov2,titj
integer::i,j,k,l,iex
i=iconsidered


	do iex=1,nof_variables-1
				rec_grads(iex,1:dimensiona,i)=rec_uleftv(1:dimensiona,iex,1,1,i)



	end do




end subroutine compute_gradients_center



subroutine compute_gradients_mix_mean_ggs_viscous_av(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the averaged primitive variables of each non-interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:gpu_max_nvar)::sols1,sols2
real,dimension(1:gpu_max_nvar,3)::sols_f
real,dimension(3)::normal_all,temp_vert
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,nx,ny,nz
integer::i,j,k,l,var2,b_code,facex,n_node,nf,lf,rowf
real,dimension(1:gpu_max_nvar)::leftv,srf_speed,srf_speedrot,rightv
real,dimension(1:gpu_max_dim)::pox,poy,poz,cords
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_nvar_total)::cright_rot,cleft_rot
integer::ibfc








i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)




	  leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
	    call cons2div(n,leftv,mp_pinfl,gammal)
	  sols1(1:nof_variables-1)=leftv(2:nof_variables)



	  leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)




do j=1,ielem_ifca(i)
			 facex=j
			 b_code=0

			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=(cos(angle1)*sin(angle2))
				normal_all(2)=(sin(angle1)*sin(angle2))
				normal_all(3)=(cos(angle2))
				nx=normal_all(1);ny=normal_all(2);nz=normal_all(3)




			if (ielem_ineighb(j,i).eq.n)then	!my cpu only
			    if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in my cpu
				  sols2(1:nof_variables)=u_c_val(ind1,1:nof_variables,ielem_ineigh(j,i))
				  else
				  !not periodic ones in my cpu

				  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)

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

 				  leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
				  b_code=ibound_icode(ielem_ibounds(j,i))
 				  call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)


				  sols2(1:nof_variables)=rightv(1:nof_variables)



				  end if
			    else
				    sols2(1:nof_variables)=u_c_val(ind1,1:nof_variables,ielem_ineigh(j,i))




			    end if
			else	!in other cpus they can only be periodic or mpi neighbours

						nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
						lf=rec_ihexl(1,ielem_indexi(j,i),i)
						rowf=halo_offset(nf) + lf - 1
						sols2(1:nof_variables)=solhir(rowf,1:nof_variables)
			end if

			  leftv(1:nof_variables)=sols2(1:nof_variables)
			call cons2div(n,leftv,mp_pinfl,gammal)
			sols2(1:nof_variables-1)=leftv(2:nof_variables)


			if ((b_code.eq.4).and.(thermal.eq.1))then
				sols2(dimensiona+1:nof_variables-nof_species-1)=wall_temp
			end if
			if ((b_code.eq.4).and.(catalytic_wall.eq.1))then
				sols2(dimensiona+3:nof_variables-1)=catalytic_con(1:nof_species)
			end if







			do k=1,3
			sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do
end do



			do k=1,dimensiona
			rec_gradsav(1:nof_variables-1,k,i)=sols_f(1:nof_variables-1,k)
			end do








end subroutine compute_gradients_mix_mean_ggs_viscous_av


subroutine compute_gradients_inner_mean_ggs_viscous_av(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the averaged primitive variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(gpu_max_nvar)::sols1,sols2,leftv
real,dimension(gpu_max_nvar,3)::sols_f
real,dimension(3)::normal_all
real::oov2,mp_pinfl,gammal,angle1,angle2
integer::i,j,k,l


i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)


	    leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,i)
	    call cons2div(n,leftv,mp_pinfl,gammal)
	  sols1(1:nof_variables)=leftv(1:nof_variables)


do j=1,ielem_ifca(i)
			angle1=ielem_faceanglex(j,i)
			angle2=ielem_faceangley(j,i)
				normal_all(1)=(cos(angle1)*sin(angle2))
				normal_all(2)=(sin(angle1)*sin(angle2))
				normal_all(3)=(cos(angle2))

			leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,ielem_ineigh(j,i))
			call cons2div(n,leftv,mp_pinfl,gammal)
			sols2(1:nof_variables)=leftv(1:nof_variables)







			do k=1,3
			sols_f(1:nof_variables,k)=sols_f(1:nof_variables,k)+((oo2*(sols2(1:nof_variables)+sols1(1:nof_variables)))*normal_all(k)*ielem_surf(j,i)*oov2)

			end do
end do

			do k=1,dimensiona
			rec_gradsav(1:nof_variables-1,k,i)=sols_f(2:nof_variables,k)
			end do







end subroutine compute_gradients_inner_mean_ggs_viscous_av












end module gradients
