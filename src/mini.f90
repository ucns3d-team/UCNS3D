program miniapp
implicit none
real,allocatable,dimension(:,:,:)::








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



end program
