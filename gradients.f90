module gradients
use library
use flow_operations
implicit none


 contains

 
 
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
   real,dimension(1:nof_variables)::sols1
   real,dimension(1:nof_variables,1:typesten)::sols2
   real,dimension(1:numneighbours-1,1:nof_variables,1:typesten)::matrix_1
   real,dimension(1:idegfree,1:nof_variables,1:7)::sol_m
   real,dimension(1:nof_variables)::leftv,rightv
   real::mp_pinfl,gammal
   integer::i,var2,iq,ll,imax,nf,lf,rowf
   real::tempxx
	 i=iconsidered



   
   imax=numneighbours-1


   i=iconsidered
   sols1=zero; sols2=zero;  sol_m=zero; matrix_1=zero;


   sols1(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
   if (wenwrt.eq.3)then
   leftv(1:nof_variables)=sols1(1:nof_variables)
   call cons2prim(n,leftv,mp_pinfl,gammal)
   sols1(1:nof_variables)=leftv(1:nof_variables)
   end if



   if (rec_local(i).eq.0)then

      do ll=1,ielem_admis(i);


        if ((ees.ne.5).or.(ll.eq.1))then
             do iq=1,imax
            sols2(1:nof_variables,ll)=u_c_val(1,1:nof_variables,rec_ihexl(ll,iq+1,i))


            if (wenwrt.eq.3)then
            leftv(1:nof_variables)=sols2(1:nof_variables,ll)
            call cons2prim(n,leftv,mp_pinfl,gammal)
            sols2(1:nof_variables,ll)=leftv(1:nof_variables)
            end if

            if(per_rot.eq.1)then
                if (rec_periodicflag(ll,iq+1,i).eq.2) then
                    tempxx=sols2(2,ll)
                    sols2(2,ll)=tempxx*cos(angle_per)-sols2(3,ll)*sin(angle_per)
                    sols2(3,ll)=tempxx*sin(angle_per)+sols2(3,ll)*cos(angle_per)
                end if
                if (rec_periodicflag(ll,iq+1,i).eq.1) then
                    tempxx=sols2(2,ll)
                    sols2(2,ll)=tempxx*cos(-angle_per)-sols2(3,ll)*sin(-angle_per)
                    sols2(3,ll)=tempxx*sin(-angle_per)+sols2(3,ll)*cos(-angle_per)
                end if
            end if

            matrix_1(iq,1:nof_variables,ll)=(sols2(1:nof_variables,ll)-sols1(1:nof_variables))

            end do

        else

               do iq=1,numneighbours2-1
            sols2(1:nof_variables,ll)=u_c_val(1,1:nof_variables,rec_ihexlc(ll,iq+1,i))
            if (wenwrt.eq.3)then
            leftv(1:nof_variables)=sols2(1:nof_variables,ll)
            call cons2prim(n,leftv,mp_pinfl,gammal)
            sols2(1:nof_variables,ll)=leftv(1:nof_variables)
            end if
            matrix_1(iq,1:nof_variables,ll)=(sols2(1:nof_variables,ll)-sols1(1:nof_variables))


         end do

        end if




    end do
     do ll=1,ielem_admis(i);
        if ((ees.ne.5).or.(ll.eq.1))then

!          call dgemm('n','n',ielem_idegfree(i),nof_variables,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:nof_variables,ll),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:nof_variables,ll),ielem_idegfree(i))


		sol_m(1:ielem_idegfree(i),1:nof_variables,ll)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:nof_variables,ll))




         else
!          call dgemm('n','n',idegfree2,nof_variables,numneighbours2-1,&
!          alpha,rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:nof_variables,ll),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:nof_variables,ll),idegfree2)

		sol_m(1:idegfree2,1:nof_variables,ll)=matmul(rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),matrix_1(1:numneighbours2-1,1:nof_variables,ll))


         end if
      end do


       do ll=1,ielem_admis(i);
       if ((ees.ne.5).or.(ll.eq.1))then
      rec_gradients(ll,1:number_of_dog,1:nof_variables,iconsidered)=sol_m(1:number_of_dog,1:nof_variables,ll)
      else
      rec_gradientsc(ll,1:idegfree2,1:nof_variables,iconsidered)=sol_m(1:idegfree2,1:nof_variables,ll)
      end if
      end do

   else


      do ll=1,ielem_admis(i);
        if ((ees.ne.5).or.(ll.eq.1))then

                do iq=1,imax


                    if (rec_ihexb(ll,iq+1,rec_local(i)).eq.n)then

                    sols2(1:nof_variables,ll)=u_c_val(1,1:nof_variables,rec_ihexl(ll,iq+1,i))
                    else


!                     sols2(1:nof_variables,ll)=iexsolhir(rec_ihexn(ll,iq+1,rec_local(i)))%sol(rec_ihexl(ll,iq+1,i),1:nof_variables)
					nf=rec_ihexn(ll,iq+1,rec_local(i))
                    lf=rec_ihexl(ll,iq+1,i)
                    rowf=halo_offset(nf) + lf - 1
                    sols2(1:nof_variables,ll)=solhir(rowf, 1:nof_variables)
                    end if

                    if (wenwrt.eq.3)then
            leftv(1:nof_variables)=sols2(1:nof_variables,ll)
            call cons2prim(n,leftv,mp_pinfl,gammal)
            sols2(1:nof_variables,ll)=leftv(1:nof_variables)
            end if
                    if(per_rot.eq.1)then
                        if (rec_periodicflag(ll,iq+1,i).eq.1) then
!             write(2900+n,*),ielem_xxc(i),ielem_yyc(i),ielem_zzc(i)
                            if (ielem_xxc(i).gt.0.0d0) then
                                tempxx=sols2(2,ll)
                                sols2(2,ll)=sols2(3,ll)
                                sols2(3,ll)=-tempxx
                            end if
                            if (ielem_xxc(i).lt.0.0d0) then
                                tempxx=sols2(2,ll)
                                sols2(2,ll)=-sols2(3,ll)
                                sols2(3,ll)=tempxx
                            end if
                        end if
                    end if

                    matrix_1(iq,1:nof_variables,ll)=(sols2(1:nof_variables,ll)-sols1(1:nof_variables))

                end do





         else
                    do iq=1,numneighbours2-1

                    if (rec_ihexbc(ll,iq+1,rec_local(i)).eq.n)then

                    sols2(1:nof_variables,ll)=u_c_val(1,1:nof_variables,rec_ihexlc(ll,iq+1,i))
                    else
!                     sols2(1:nof_variables,ll)=iexsolhir(rec_ihexnc(ll,iq+1,rec_local(i)))%sol(rec_ihexlc(ll,iq+1,i),1:nof_variables)

                    nf=rec_ihexnc(ll,iq+1,rec_local(i))
                    lf=rec_ihexlc(ll,iq+1,i)
                    rowf=halo_offset(nf) + lf - 1
                    sols2(1:nof_variables,ll)=solhir(rowf, 1:nof_variables)




                    end if
                                if (wenwrt.eq.3)then
                        leftv(1:nof_variables)=sols2(1:nof_variables,ll)
                        call cons2prim(n,leftv,mp_pinfl,gammal)
                        sols2(1:nof_variables,ll)=leftv(1:nof_variables)
                        end if

                    matrix_1(iq,1:nof_variables,ll)=(sols2(1:nof_variables,ll)-sols1(1:nof_variables))

                end do




         end if


        end do

         do ll=1,ielem_admis(i);
        if ((ees.ne.5).or.(ll.eq.1))then


!          call dgemm('n','n',ielem_idegfree(i),nof_variables,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:nof_variables,ll),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:nof_variables,ll),ielem_idegfree(i))


		sol_m(1:ielem_idegfree(i),1:nof_variables,ll)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:nof_variables,ll))

         else


!          call dgemm('n','n',idegfree2,nof_variables,numneighbours2-1,&
!          alpha,rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:nof_variables,ll),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:nof_variables,ll),idegfree2)

		sol_m(1:idegfree2,1:nof_variables,ll)=matmul(rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),matrix_1(1:numneighbours2-1,1:nof_variables,ll))



         end if
      end do


       do ll=1,ielem_admis(i);
       if ((ees.ne.5).or.(ll.eq.1))then
      rec_gradients(ll,1:number_of_dog,1:nof_variables,iconsidered)=sol_m(1:number_of_dog,1:nof_variables,ll)

      else

      rec_gradientsc(ll,1:idegfree2,1:nof_variables,iconsidered)=sol_m(1:idegfree2,1:nof_variables,ll)

      end if
      end do

   end if







end subroutine compute_gradients_mean_lsq










subroutine compute_gradients_inner_mean_lsq_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitve variables of each interior cell using the least-squares
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:nof_variables)::sols1,sols2
real,dimension(1:numneighbours-1,1:nof_variables)::matrix_1
real,dimension(1:idegfree,1:nof_variables)::sol_m
integer::i,var2,iq,lq,ll,imax,nf,lf,rowf
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::leftv










imax=number_of_nei-1


i=iconsidered
sols1=zero;
sols2=zero







		ll=1


		if (rec_local(i).eq.0)then
		matrix_1=zero;
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
		call cons2div(n,leftv,mp_pinfl,gammal)

	       sols1(1:nof_variables-1)=leftv(2:nof_variables)

               do iq=1,imax
                leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))

		call cons2div(n,leftv,mp_pinfl,gammal)
	       sols2(1:nof_variables-1)=leftv(2:nof_variables)


  	        matrix_1(iq,1:nof_variables-1)=((sols2(1:nof_variables-1)-sols1(1:nof_variables-1)))

		end do



!             call dgemm('n','n',ielem_idegfree(i),nof_variables,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:nof_variables),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:nof_variables),ielem_idegfree(i))


		sol_m(1:ielem_idegfree(i),1:nof_variables-1)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:nof_variables-1))


		do var2=1,nof_variables-1
		rec_gradf(var2,1:number_of_dog,iconsidered)=sol_m(1:number_of_dog,var2)
		end do



		else

		matrix_1=zero;
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,1,i))
		call cons2div(n,leftv,mp_pinfl,gammal)

	       sols1(1:nof_variables-1)=leftv(2:nof_variables)

               do iq=1,imax
		  if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
		  leftv(1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,iq+1,i))

		  else
! 		  leftv(1:nof_variables)=iexsolhir(rec_ihexn(1,iq+1,rec_local(i)))%sol(rec_ihexl(1,iq+1,i),1:nof_variables)

		  nf=rec_ihexn(1,iq+1,rec_local(i))
          lf=rec_ihexl(1,iq+1,i)
          rowf=halo_offset(nf) + lf - 1
          leftv(1:nof_variables)=solhir(rowf, 1:nof_variables)



		  end if


		  call cons2div(n,leftv,mp_pinfl,gammal)
	       sols2(1:nof_variables-1)=leftv(2:nof_variables)


  	        matrix_1(iq,1:nof_variables-1)=((sols2(1:nof_variables-1)-sols1(1:nof_variables-1)))
		end do



!          call dgemm('n','n',ielem_idegfree(i),nof_variables,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:nof_variables),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:nof_variables),ielem_idegfree(i))

			sol_m(1:ielem_idegfree(i),1:nof_variables-1)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:nof_variables-1))




		do var2=1,nof_variables-1
		rec_gradf(var2,1:number_of_dog,iconsidered)=sol_m(1:number_of_dog,var2)
		end do

		end if



end subroutine compute_gradients_inner_mean_lsq_viscous





subroutine compute_gradients_inner_turb_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the turbulence variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(turbulenceequations+passivescalar)::sols1,sols2
real,dimension(turbulenceequations+passivescalar,nof_variables)::sols_f
real,dimension(dimensiona)::normal_all
real::oov2,angle1,angle2
integer::i,j,k,l,var2
real,dimension(1:nof_variables)::leftv







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
real,dimension(turbulenceequations+passivescalar)::sols1,sols2
real,dimension(1:numneighbours-1,1:turbulenceequations+passivescalar)::matrix_1
real,dimension(1:idegfree,1:turbulenceequations+passivescalar)::sol_m
integer::i,var2,ll,iq,il,ih,nf,lf,rowf
integer::imax











imax=number_of_nei-1
i=iconsidered
sols1=zero;
sols2=zero






	      sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,1,i))
	      if (rec_local(i).eq.0)then

	      do ll=1,ielem_admis(i);
		matrix_1=zero;
		if ((ees.ne.5).or.(ll.eq.1))then
		do iq=1,imax
		  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(ll,iq+1,i))
		  matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
		end do
		else
                do iq=1,numneighbours2-1
		  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexlc(ll,iq+1,i))
		  matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
		end do


                end if

             if ((ees.ne.5).or.(ll.eq.1))then



! c          call dgemm('n','n',ielem_idegfree(i),turbulenceequations+passivescalar,imax,&
! c          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
! c          ielem_idegfree(i),matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! c imax,beta,sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar),ielem_idegfree(i))

		sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:turbulenceequations+passivescalar))


            else


!         call dgemm('n','n',idegfree2,turbulenceequations+passivescalar,numneighbours2-1,&
!          alpha,rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:turbulenceequations+passivescalar),idegfree2)

			sol_m(1:idegfree2,1:turbulenceequations+passivescalar)=matmul(rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar))


            end if


                 if ((ees.ne.5).or.(ll.eq.1))then
		rec_gradients2(ll,1:number_of_dog,1:turbulenceequations+passivescalar,iconsidered)=sol_m(1:number_of_dog,1:turbulenceequations+passivescalar)

		else
		rec_gradientsc2(ll,1:idegfree2,1:turbulenceequations+passivescalar,iconsidered)=sol_m(1:idegfree2,1:turbulenceequations+passivescalar)


		end if
	        end do

	       else

	       do ll=1,ielem_admis(i);
		matrix_1=zero;
		 if ((ees.ne.5).or.(ll.eq.1))then
		do iq=1,imax

		  if (rec_ihexb(ll,iq+1,rec_local(i)).eq.n)then
		  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(ll,iq+1,i))

		  else
! 		  sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(ll,iq+1,i))%sol(rec_ihexl(ll,iq+1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)



		  nf=rec_ihexn(ll,iq+1,rec_local(i))
          lf=rec_ihexl(ll,iq+1,i)
          rowf=halo_offset(nf) + lf - 1
          sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)







		  end if

		   matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
		end do
		else
		do iq=1,numneighbours2-1

		  if (rec_ihexbc(ll,iq+1,rec_local(i)).eq.n)then
		  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexlc(ll,iq+1,i))

		  else
! 		  sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexnc(ll,iq+1,i))%sol(rec_ihexlc(ll,iq+1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

		  nf=rec_ihexnc(ll,iq+1,rec_local(i))
          lf=rec_ihexlc(ll,iq+1,i)
          rowf=halo_offset(nf) + lf - 1
          sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)



		  end if

		   matrix_1(iq,1:turbulenceequations+passivescalar)=(sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar))
		end do



		end if


                if ((ees.ne.5).or.(ll.eq.1))then




!          call dgemm('n','n',ielem_idegfree(i),turbulenceequations+passivescalar,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar),ielem_idegfree(i))

			sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:turbulenceequations+passivescalar))


                else


!             call dgemm('n','n',idegfree2,turbulenceequations+passivescalar,numneighbours2-1,&
!          alpha,rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),&
!          idegfree2,matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar),&
! numneighbours2-1,beta,sol_m(1:idegfree2,1:turbulenceequations+passivescalar),idegfree2)


			sol_m(1:idegfree2,1:turbulenceequations+passivescalar)=matmul(rec_invmat_stenciltc(1:idegfree2,1:numneighbours2-1,ll,i),matrix_1(1:numneighbours2-1,1:turbulenceequations+passivescalar))


                end if



		  if ((ees.ne.5).or.(ll.eq.1))then
		rec_gradients2(ll,1:number_of_dog,1:turbulenceequations+passivescalar,iconsidered)=sol_m(1:number_of_dog,1:turbulenceequations+passivescalar)

		else
		rec_gradientsc2(ll,1:idegfree2,1:turbulenceequations+passivescalar,iconsidered)=sol_m(1:idegfree2,1:turbulenceequations+passivescalar)


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
real,dimension(turbulenceequations+passivescalar)::sols1,sols2
real,dimension(1:numneighbours-1,1:turbulenceequations+passivescalar)::matrix_1
real,dimension(1:idegfree,1:turbulenceequations+passivescalar)::sol_m
integer::i,var2,iq,lq,ll,imax,ideg,nf,lf,rowf




imax=numneighbours-1

i=iconsidered
sols1=zero;
sols2=zero
ideg=ielem_idegfree(i)
		 ll=1







		if (rec_local(i).eq.0)then
		matrix_1=zero;
		sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,1,i))/&
		u_c_val(1,1,rec_ihexl(1,1,i))




               do iq=1,imax
                sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,iq+1,i))/&
                u_c_val(1,1,rec_ihexl(1,iq+1,i))

  	        matrix_1(iq,1:turbulenceequations+passivescalar)=((sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar)))
		end do




!         call dgemm('n','n',ielem_idegfree(i),turbulenceequations+passivescalar,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar),ielem_idegfree(i))



		sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:turbulenceequations+passivescalar))



		do var2=1,turbulenceequations+passivescalar
		rec_gradientsturb(1,1:number_of_dog,var2,iconsidered)=sol_m(1:number_of_dog,var2)
		end do


		else

		matrix_1=zero;
		sols1(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,1,i))&
		/u_c_val(1,1,rec_ihexl(1,1,i))


               do iq=1,imax
		  if (rec_ihexb(1,iq+1,rec_local(i)).eq.n)then
		  sols2(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,iq+1,i))&
		/u_c_val(1,1,rec_ihexl(1,iq+1,i))

		  else
! 		  sols2(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,iq+1,i))%sol(rec_ihexl(1,iq+1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
! 		  iexsolhir(rec_ihexn(1,iq+1,i))%sol(rec_ihexl(1,iq+1,i),1)

		  nf=rec_ihexn(ll,iq+1,rec_local(i))
          lf=rec_ihexl(ll,iq+1,i)
          rowf=halo_offset(nf) + lf - 1
          sols2(1:turbulenceequations+passivescalar)=solhir(rowf, nof_variables+1:nof_variables+turbulenceequations+passivescalar)/&
          solhir(rowf, 1)




		  end if

  	        matrix_1(iq,1:turbulenceequations+passivescalar)=((sols2(1:turbulenceequations+passivescalar)-sols1(1:turbulenceequations+passivescalar)))
		end do




!           call dgemm('n','n',ielem_idegfree(i),turbulenceequations+passivescalar,imax,&
!          alpha,rec_invmat_stencilt(1:idegfree,1:imax,ll),&
!          ielem_idegfree(i),matrix_1(1:imax,1:turbulenceequations+passivescalar),&
! imax,beta,sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar),ielem_idegfree(i))


			sol_m(1:ielem_idegfree(i),1:turbulenceequations+passivescalar)=matmul(rec_invmat_stencilt(1:idegfree,1:imax,ll,i),matrix_1(1:imax,1:turbulenceequations+passivescalar))




		do var2=1,turbulenceequations+passivescalar
		rec_gradientsturb(1,1:number_of_dog,var2,iconsidered)=sol_m(1:number_of_dog,var2)
		end do


		end if







end subroutine compute_gradients_turb_lsq_viscous


subroutine compute_gradients_inner_mean_ggs_viscous(n,iconsidered,number_of_dog,number_of_nei)!check_all
!> @brief
!> this subroutine computes the gradients of the primitive variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,number_of_dog,number_of_nei
real,dimension(1:nof_variables)::sols1,sols2,dudl,aver1,phi_f
real,dimension(1:nof_variables,dimensiona)::sols_f
real,dimension(3)::normal_all,dih_vec,e_ih,sf
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,aorth,dih
integer::i,j,k,l
real,dimension(1:nof_variables)::leftv


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
real,dimension(nof_variables-1)::sols1,sols2
real,dimension(1:nof_variables,1:numneighbours-1)::matrix_1
real,dimension(1:nof_variables,1:idegfree)::matrix_2
real,dimension(1:idegfree,1:nof_variables)::sol_m
real,dimension(1:nof_variables)::matrix_3
integer::i,var2,ii,k0,g0,ttk,ivvm,iq,lq,irg
real::attt
integer::ll,imax,nf,lf,rowf
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::leftv






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
		  matrix_2=zero
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
		sol_m(1:number_of_dog-1,var2)=matmul(rec_velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1,rec_wall(i)),matrix_2(var2,1:number_of_dog-1))
		end if
		if ((var2.gt.dimensiona).and.(var2.le.nof_variables-nof_species-1))then
		if (thermal.eq.1)then
		sol_m(1:number_of_dog-1,var2)=matmul(rec_velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1,rec_wall(i)),matrix_2(var2,1:number_of_dog-1))
		else
		sol_m(1:number_of_dog-1,var2)=matmul(rec_tempsqmat(1:number_of_dog-1,1:number_of_dog-1,rec_wall(i)),matrix_2(var2,1:number_of_dog-1))
		end if
		end if

		if (var2.gt.nof_variables-nof_species-1)then

		if (catalytic_wall.eq.1)then
			sol_m(1:number_of_dog-1,var2)=matmul(rec_velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1,rec_wall(i)),matrix_2(var2,1:number_of_dog-1))

		else

		sol_m(1:number_of_dog-1,var2)=matmul(rec_tempsqmat(1:number_of_dog-1,1:number_of_dog-1,rec_wall(i)),matrix_2(var2,1:number_of_dog-1))
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
real,dimension(turbulenceequations+passivescalar)::sols1,sols2
real,dimension(turbulenceequations+passivescalar,3)::sols_f
real,dimension(3)::normal_all,temp_vert
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,nx,ny,nz
integer::i,j,k,l,var2,b_code,facex,n_node,imax,nf,lf,rowf
real,dimension(1:nof_variables)::leftv,srf_speed,srf_speedrot,rightv
real,dimension(1:dimensiona)::pox,poy,poz,cords
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(turbulenceequations)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cright_rot,cleft_rot
integer::ibfc


if (dimensiona.eq.3)then

i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)

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


				  cords=cordinates3(n,nodes_list,n_node)
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
				  cords=cordinates2(n,nodes_list,n_node)
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
real,dimension(1:turbulenceequations+passivescalar)::sols1,sols2
real,dimension(1:turbulenceequations+passivescalar,1:numneighbours-1)::matrix_1
real,dimension(1:turbulenceequations+passivescalar,1:idegfree)::matrix_2
real,dimension(1:idegfree,1:turbulenceequations+passivescalar)::sol_m
real,dimension(1:turbulenceequations+passivescalar)::matrix_3
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
		  matrix_2=zero

		  do iq=1,imax

		      do lq=1,number_of_dog-1
		      matrix_2(var2,lq)=matrix_2(var2,lq)+matrix_1(var2,iq)*rec_vellsq(iq,lq,rec_wall(i))
		      end do

		  end do


		sol_m(1:number_of_dog-1,var2)=matmul(rec_velinvlsqmat(1:number_of_dog-1,1:number_of_dog-1,rec_wall(i)),matrix_2(var2,1:number_of_dog-1))



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
		  attt=-sols1(var2)
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
real,dimension(1:nof_variables)::sols1,sols2,dudl,aver1
real,dimension(1:nof_variables,3)::sols_f
real,dimension(3)::normal_all,dih_vec,e_ih,sf
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,nx,ny,nz,aorth,dih
integer::i,j,k,l,b_code,facex,n_node,imax,nf,lf,rowf
real,dimension(1:nof_variables)::leftv,srf_speed,srf_speedrot,rightv,phi_f
real,dimension(1:dimensiona)::pox,poy,poz,cords
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(turbulenceequations)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cright_rot,cleft_rot
integer::ibfc





if (dimensiona.eq.3)then

i=iconsidered
sols_f=zero;sols1=zero;sols2=zero

rec_grads(:,:,i)=zero

oov2=1.0d0/ielem_totvolume(i)




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
				  if(per_rot.eq.1)then
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
 				  cords(1:3)=cordinates3(n,nodes_list,n_node)

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

			      if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if ((ibound_icode(ielem_ibounds(j,i)).eq.5).or.(ibound_icode(ielem_ibounds(j,i)).eq.50))then	!periodic in other cpu

! 					sols2(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),1:nof_variables)


					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)




				      if(per_rot.eq.1)then
                        sols2(2:4)=rotate_per_1(sols2(2:4),ibound_icode(ielem_ibounds(j,i)),angle_per)
				      end if
				  end if
			      else


! 					sols2(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
!					(rec_ihexl(1,ielem_indexi(j,i)),1:nof_variables)


					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)



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

	i=iconsidered
sols_f=zero
oov2=1.0d0/ielem_totvolume(i)


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
				  cords=cordinates2(n,nodes_list,n_node)
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

			      if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in other cpu

! 					sols2(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),1:nof_variables)


					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)

				  end if
			      else


! 					sols2(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),1:nof_variables)


					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)


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


 			do k=1,dimensiona
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












end subroutine compute_gradients_mix_mean_ggs_viscous


subroutine compute_gradients_center(n,iconsidered)
!> @brief
!> this subroutine computes the gradients of the primitive variables of each interior cell using the green-gauss algorithm
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,dimension(1:nof_variables)::sols1,sols2,dudl,aver1
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
real,dimension(1:nof_variables)::sols1,sols2
real,dimension(1:nof_variables,3)::sols_f
real,dimension(3)::normal_all,temp_vert
real::oov2,titj,mp_pinfl,gammal,angle1,angle2,nx,ny,nz
integer::i,j,k,l,var2,b_code,facex,n_node,nf,lf,rowf
real,dimension(1:nof_variables)::leftv,srf_speed,srf_speedrot,rightv
real,dimension(1:dimensiona)::pox,poy,poz,cords
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(turbulenceequations)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cright_rot,cleft_rot
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
 				  cords(1:3)=cordinates3(n,nodes_list,n_node)

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

			      if (ielem_ibounds(j,i).gt.0)then	!check for boundaries
				  if (ibound_icode(ielem_ibounds(j,i)).eq.5)then	!periodic in other cpu

!					sols2(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
!					(rec_ihexl(1,ielem_indexi(j,i)),1:nof_variables)

					 nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)



				  end if
			      else


! 					sols2(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(j,i)))%sol&
! 					(rec_ihexl(1,ielem_indexi(j,i)),1:nof_variables)

					nf=rec_ihexn(1,ielem_indexi(j,i),rec_local(i))
					lf=rec_ihexl(1,ielem_indexi(j,i),i)
					rowf=halo_offset(nf) + lf - 1
					sols2(1:nof_variables)=solhir(rowf,1:nof_variables)

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
real,dimension(nof_variables)::sols1,sols2,leftv
real,dimension(nof_variables,3)::sols_f
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



subroutine example12(iconsidered)
!> @brief
!> this subroutine computes the gradients of the conserved variables of each cell using the least-squares
   implicit none
#ifdef gpu
!$omp declare target
#endif
   integer,intent(in)::iconsidered
   !real,dimension(1:5)::sols1
   !real,dimension(1:5)::sols2
   real,dimension(1:71,1:5)::matrix_1
   !real,dimension(1:35,1:5,1:7)::sol_m
   !real,dimension(1:5)::leftv,rightv
   !real::mp_pinfl,gammal
   integer::i,idegx,var2,iq,ll,imax,nf,lf,rowf,varx,eesxg,c,k,m
   !real::tempxx,sum
	 i=iconsidered
	varx=5
   imax=ielem_inumneighbours(i)-1
   idegx=ielem_idegfree(i)
   !sols1(1:varx)=u_c_val(1,1:varx,rec_ihexl(1,1,i))
   if (rec_local(i).eq.0)then

				do ll=1,ielem_admis(i);
						do iq=1,imax
                                                        
				                        matrix_1(iq,1:varx)=u_c_val(1,1:varx,rec_ihexl(ll,iq+1,i))-u_c_val(1,1:varx,i)
					!	matrix_1(iq,1:varx)=(sols2(1:varx)-sols1(1:varx))
		                                end do
                                        rec_gradients(ll,1:idegx,1:5,i)=matmul(rec_invmat_stencilt(1:idegx,1:imax,ll,i),matrix_1(1:imax,1:5))
				end do

                        end if
!					do ll=1,ielem_admis(i);
!					rec_gradients(ll,1:ielem_idegfree(i),1:5,iconsidered)=sol_m(1:ielem_idegfree(i),1:5,ll)

!					end do
!                        do ll = 1, ielem_admis(i)
!                                 do c = 1, 5
!                                        do m = 1, ielem_idegfree(i)
!                                             sum = 0.0d0
!                                            do k = 1, imax
!                                                                 sum = sum + rec_invmat_stencilt(m,k,ll,i) * matrix_1(k,c,ll)
!                                             end do
!                                                 rec_gradients(ll,m,c,iconsider







end subroutine example12


subroutine example12_row(iconsidered, m_in)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  integer, intent(in) :: iconsidered, m_in
  integer :: i, ll, c, k, imax, cellnbr
  real :: sols1(5)
  real :: sum

  i    = iconsidered
  imax = ielem_inumneighbours(i) - 1

  if (rec_local(i) .ne. 0) return

  ! Center cell values
  sols1(1:5) = u_c_val(1, 1:5, rec_ihexl(1,1,i))

  do ll = 1, ielem_admis(i)
    do c = 1, 5
      sum = 0.0d0
      do k = 1, imax
        cellnbr = rec_ihexl(ll, k+1, i)
        sum = sum + rec_invmat_stencilt(m_in, k, ll, i) * (u_c_val(1, c, cellnbr) - sols1(c))
      end do
      rec_gradients(ll, m_in, c, iconsidered) = sum
    end do
  end do

end subroutine example12_row


subroutine example1x(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, m

#ifdef gpu
  !$omp target teams distribute parallel do private(i)
#else
  !$omp parallel do private(i)
#endif
  do ii = 1, nof_interior
    do m = 1, ielem_idegfree(el_int(ii))   ! typically 35
      i = el_int(ii)
      call example12_row(i, m)
    end do
  end do
#ifdef gpu
  !$omp end target teams distribute parallel do
#else
  !$omp end parallel do
#endif

end subroutine example1x

subroutine example1xy(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, m
  integer :: ll, c, k, imax, cellnbr
  real:: sols1(5)
  real :: sum

#ifdef gpu
  !$omp target teams distribute parallel do schedule (static) &
  !$omp& private(i,ll,c,k,imax,cellnbr,sols1,sum)
#else
  !$omp parallel do  private(i,ll,c,k,imax,cellnbr,sols1,sum)
#endif
  do ii = 1, nof_interior
    do m = 1, ielem_idegfree(el_int(ii))   ! typically 35

      i = el_int(ii)

      if (rec_local(i) .ne. 0) cycle

      imax = ielem_inumneighbours(i) - 1

      ! Center cell values (5 vars)
      sols1(1:5) = u_c_val(1, 1:5, rec_ihexl(1,1,i))

      do ll = 1, ielem_admis(i)
        do c = 1, 5
          sum = 0.0d0
          do k = 1, imax
            cellnbr = rec_ihexl(ll, k+1, i)
            sum = sum + rec_invmat_stencilt(m, k, ll, i) * (u_c_val(1, c, cellnbr) - sols1(c))
          end do
          rec_gradients(ll, m, c, i) = sum
        end do
      end do

    end do
  end do
#ifdef gpu
  !$omp end target teams distribute parallel do
#else
  !$omp end parallel do
#endif

end subroutine example1xy

subroutine example1(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i

#ifdef gpu
!$omp target teams distribute parallel do private(i) schedule(static)
#else
!$omp do
#endif
  do ii = 1,nof_interior
    i = el_int(ii)
     call example12(i)
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
  !$omp end do
#endif



end subroutine example1








end module gradients
