module dg_functions

use basis
use declaration
use derivatives
use lapck
implicit none


contains




function dg_sol(n,iconsidered,x1,y1,z1)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> @brief
!> this function returns the dg solution at a given point (x_in, y_in)\n
!> requires: x_in, y_in: coordinates of the point where the solution is requested, num_variables: number of solution variables, num_dofs: number of basis terms
real,dimension(1:55)::basis_temp
integer,intent(in)::n,iconsidered
integer::i_dof, i_var,icompwrt,number_of_dog,number,k
real,dimension(1:nof_variables)::dg_sol
real::tempdp
real,intent(in)::x1,y1,z1


 icompwrt=-2
number=ielem_iorder(iconsidered)
number_of_dog=ielem_idegfree(iconsidered)





if (dimensiona.eq.2)then
basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
else
basis_temp = basis_rec(n, x1, y1,z1, number, iconsidered, number_of_dog,icompwrt)
end if






do i_var = 1, nof_variables
    tempdp=0.0d0
    do k=1,number_of_dog
      tempdp=tempdp+basis_temp(k)*u_c_valdg(1,i_var,k+1,iconsidered)
    end do
    dg_sol(i_var) = u_c_valdg(1,i_var,1,iconsidered)+tempdp

    !dg_sol(i_var) = u_c_valdg(1,i_var,1,iconsidered)+ dot_product(basis_temp(1:number_of_dog),u_c_valdg(1,i_var,2:number_of_dog+1,iconsidered))

end do

 icompwrt=0

end function dg_sol


function dg_solface(n,facex,pointx,iconsidered,number_of_dog)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> @brief
!> this function returns the dg solution at a given surface point (x_in, y_in)\n
!> requires: x_in, y_in: coordinates of the point where the solution is requested, num_variables: number of solution variables, num_dofs: number of basis terms
real,dimension(1:55)::basis_temp
integer::i_dof, i_var,icompwrt,k
integer,intent(in)::n,facex,pointx,iconsidered,number_of_dog
real,dimension(1:nof_variables)::dg_solface
real::x1,y1,z1,tempdp
integer::number

 icompwrt=-2

x1= rec_qpoints(facex,pointx,1,iconsidered)
 y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
 z1= rec_qpoints(facex,pointx,3,iconsidered)
end if







number=ielem_iorder(iconsidered)
if (dimensiona.eq.2)then
basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
else
basis_temp = basis_rec(n, x1, y1,z1, number, iconsidered, number_of_dog,icompwrt)
end if



do i_var = 1, nof_variables
    tempdp=0.0d0
    do k=1,number_of_dog
      tempdp=tempdp+basis_temp(k)*u_c_valdg(1,i_var,k+1,iconsidered)
    end do
    dg_solface(i_var) = u_c_valdg(1,i_var,1,iconsidered)+tempdp

!     dg_solface(i_var) = u_c_valdg(1,i_var,1,iconsidered) + dot_product(basis_temp(1:number_of_dog),u_c_valdg(1,i_var,2:number_of_dog+1,iconsidered))
end do





    icompwrt=0

end function dg_solface


function dg_sol_der(x1,y1,z1,number_of_dog,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> @brief
!> this function returns the derivative of the dg solution at a given point (x1, y1, [z1])
    integer::i_dof, i_var, i_dim,icompwrt,number,k
    integer,intent(in)::number_of_dog,iconsidered
    real::tempdp
    real,intent(in)::x1,y1,z1
    real,dimension(1:55,1:dimensiona)::basis_temp
    real,dimension(1:55)::basis_t
    real,dimension(1:nof_variables,1:dimensiona)::dg_sol_der

    icompwrt=-2
    number=ielem_iorder(iconsidered)


    if (dimensiona.eq.2)then
    basis_t = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
    else
    basis_t = basis_rec(n, x1, y1,z1, number, iconsidered, number_of_dog,icompwrt)
    end if



    if (br2_yn /= 1) then
        do i_dof = 1, number_of_dog
            if (dimensiona.eq.2) then
                if( poly == 1) then
                    basis_temp(i_dof,1) = df2dx(x1,y1,i_dof,iconsidered)
                    basis_temp(i_dof,2) = df2dy(x1,y1,i_dof,iconsidered)
                else if( poly == 4) then
                    basis_temp(i_dof,1) = tl2dx(x1,y1,i_dof,iconsidered)
                    basis_temp(i_dof,2) = tl2dy(x1,y1,i_dof,iconsidered)
                end if
            else
                basis_temp(i_dof,1) = tl3dx(x1,y1,z1,i_dof,iconsidered)
                basis_temp(i_dof,2) = tl3dy(x1,y1,z1,i_dof,iconsidered)
                basis_temp(i_dof,3) = tl3dz(x1,y1,z1,i_dof,iconsidered)
            end if
        end do
    end if

    icompwrt=-2
    do i_var = 1, nof_variables
        do i_dim = 1, dimensiona
            if (br2_yn == 1) then
                if(dimensiona.eq.2)then

! !                     dg_sol_der(i_var, i_dim) = u_c_br2_aux_var(1,i_var,i_dim,iconsidered) + dot_product(basis_rec2d(n,x1, y1, number, iconsidered, number_of_dog,icompwrt), u_c_br2_aux_var(2:num_dg_dofs,i_var,i_dim,iconsidered))


                    tempdp=0.0d0
                    do k=1,number_of_dog
                      tempdp=tempdp+basis_t(k)*u_c_br2_aux_var(k+1,i_var,i_dim,iconsidered)
                    end do
                    dg_sol_der(i_var, i_dim) = u_c_br2_aux_var(1,i_var,i_dim,iconsidered)+tempdp





                else
!
! !                     dg_sol_der(i_var, i_dim) = u_c_br2_aux_var(1,i_var,i_dim,iconsidered) + dot_product(basis_rec(n,x1, y1, z1, number, iconsidered, number_of_dog,icompwrt), u_c_br2_aux_var(2:num_dg_dofs,i_var,i_dim,iconsidered))


                     tempdp=0.0d0
                    do k=1,number_of_dog
                      tempdp=tempdp+basis_t(k)*u_c_br2_aux_var(k+1,i_var,i_dim,iconsidered)
                    end do
                    dg_sol_der(i_var, i_dim) = u_c_br2_aux_var(1,i_var,i_dim,iconsidered)+tempdp





                 end if
            else

!                 dg_sol_der(i_var, i_dim) = dot_product(basis_temp(1:number_of_dog,i_dim), u_c_valdg(1,i_var,2:num_dg_dofs,iconsidered))


               tempdp=0.0d0
                    do k=1,number_of_dog
                      tempdp=tempdp+basis_t(k)*u_c_valdg(1,i_var,k+1,iconsidered)
                    end do
                    dg_sol_der(i_var, i_dim) = tempdp




            end if
        end do
    end do
    icompwrt=0

end function dg_sol_der


function br2_local_lift(n, n_qp,facex,iconsidered,wequa2d)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> requires: iconsidered, facex, and weights_dg or wequa2d
integer,intent(in)::n, n_qp,facex,iconsidered
real,dimension(1:numberofpoints2),intent(in)::wequa2d
integer::i, l, ngp, iqp, i_dim,icompwrt,number_of_dog
real::angle1,angle2
real,dimension(dimensiona)::nnn
real,dimension(numberofpoints2)::weights_temp
real,dimension(nof_variables,dimensiona)::br2_local_lift
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv,srf_speedrot
real::mp_pinfr,gammar
real,dimension(1:55)::basis_temp
 real::x1,y1,z1
integer::b_code,pointx,number
real::nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar)::cturbl,cturbr
real,dimension(1:dimensiona)::pox,poy,poz



    i = iconsidered
    l = facex
    number_of_dog=ielem_idegfree(i)
    number=ielem_iorder(iconsidered)



	if (dimensiona == 3) then
        angle1=ielem_faceanglex(l,i)
        angle2=ielem_faceangley(l,i)
        nnn(1)=(cos(angle1)*sin(angle2))
        nnn(2)=(sin(angle1)*sin(angle2))
        nnn(3)=(cos(angle2))
    else
        nnn(1)=ielem_faceanglex(l,i)
        nnn(2)=ielem_faceangley(l,i)
    end if

    br2_local_lift = 0.0d0
    do ngp = 1, n_qp

        pointx = ngp

        icompwrt=-2
        x1 = rec_qpoints(facex,pointx,1,iconsidered)
        y1 = rec_qpoints(facex,pointx,2,iconsidered)
        if (dimensiona == 3) z1 = rec_qpoints(facex,pointx,3,iconsidered)

        if (dimensiona .eq. 2)then
            basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
        else
            basis_temp = basis_rec(n, x1, y1, z1, number, iconsidered, number_of_dog,icompwrt)
        end if







        call get_left_right_states(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)


        leftv = cleft
        rightv = cright
        call cons2div(n,leftv,mp_pinfl,gammal)
        call cons2div(n,rightv,mp_pinfr,gammar)


        do i_dim = 1, dimensiona
!             rhllcflux = nnn(i_dim) * (rightv - leftv)
            br2_local_lift(1:nof_variables,i_dim) = br2_local_lift(1:nof_variables,i_dim) + nnn(i_dim) * (rightv - leftv) * weights_temp(pointx) !sum(dg_surf_flux(n),dim=1)
!             if( n == 3 .and. iconsidered == 20) write(400+n,*) br2_local_lift, "br2_local_lift", dg_surf_flux(n), "dg_surf_flux", sum(dg_surf_flux(n),dim=1), "sum(dg_surf_flux(n),dim=1)"
        end do

    end do
    br2_local_lift = oo2 * br2_local_lift * ielem_surf(facex,iconsidered) / ielem_totvolume(iconsidered) ! / dg_vol_integral3(n)

    icompwrt=0


end function br2_local_lift



function dg_surf_flux(n,iconsidered,facex,pointx,weights_temp,rhllcflux)
implicit none
!> @brief
!> calculates the rhs flux term to be integrated in the dg formulation
#ifdef gpu
!$omp declare target
#endif
real,dimension(idegfree+1,nof_variables)::dg_surf_flux
integer::i,icompwrt
integer,intent(in)::n,iconsidered,facex,pointx
real,dimension(1:nof_variables),intent(in)::rhllcflux
real,dimension(1:numberofpoints2),intent(in)::weights_temp
real::x1,y1,z1
integer::number,number_of_dog

number_of_dog=ielem_idegfree(iconsidered)

number=ielem_iorder(iconsidered)

 icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
 y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
 z1= rec_qpoints(facex,pointx,3,iconsidered)
end if





number=ielem_iorder(iconsidered)
if (dimensiona.eq.2)then
do i = 1, nof_variables
    dg_surf_flux(1,i) = rhllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)

    dg_surf_flux(2:number_of_dog+1,i) = rhllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)*basis_rec2d(n,x1,y1,number,iconsidered,number_of_dog,icompwrt)

end do

else
do i = 1, nof_variables
    dg_surf_flux(1,i) = rhllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)

    dg_surf_flux(2:number_of_dog+1,i) = rhllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)*basis_rec(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt)
end do

end if



    icompwrt=0




end function dg_surf_flux


function dg_surf_fluxv(n,iconsidered,facex,pointx,weights_temp,hllcflux)
implicit none
!> @brief
!> calculates the rhs flux term to be integrated in the dg formulation
#ifdef gpu
!$omp declare target
#endif
real,dimension(idegfree+1,1:nof_variables)::dg_surf_fluxv
integer::i,icompwrt
integer,intent(in)::n,iconsidered,facex,pointx
real,dimension(1:nof_variables),intent(in)::hllcflux
real,dimension(1:nof_variables)::rhllcflux
real,dimension(1:numberofpoints2),intent(in)::weights_temp
real::x1,y1,z1
integer::number,number_of_dog
 icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
 y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
 z1= rec_qpoints(facex,pointx,3,iconsidered)
end if





number=ielem_iorder(iconsidered)
if (dimensiona.eq.2)then
do i = 1, nof_variables
    dg_surf_fluxv(1,i) = hllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)

    dg_surf_fluxv(2:number_of_dog+1,i) = hllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)*basis_rec2d(n,x1,y1,number,iconsidered,number_of_dog,icompwrt)

end do

else
do i = 1, nof_variables
    dg_surf_fluxv(1,i) = hllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)

    dg_surf_fluxv(2:number_of_dog+1,i) = hllcflux(i) * weights_temp(pointx)*ielem_surf(facex,iconsidered)*basis_rec(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt)
end do

end if





    icompwrt=0





end function dg_surf_fluxv




function dg_vol_integral(n,iconsidered)
implicit none
!> @brief
!> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
integer::i,nqp,i_qp,icompwrt,number,number_of_dog
real,dimension(idegfree+1,nof_variables)::dg_vol_integral
real,dimension(1:nof_variables)::leftv
real,dimension(1:nof_variables)::flux_term_x,flux_term_y,flux_term_z
real,dimension(1:nof_variables,1:dimensiona)::leftv_der
real::x1,y1,z1
real::mp_pinfl,gammal



    icompwrt=-2



    nqp = ielem_itotalpoints(iconsidered)


    number=ielem_iorder(iconsidered)
    number_of_dog = ielem_idegfree(iconsidered)


    dg_vol_integral(:,:) = 0.0d0




    do i=1,number_of_dog

        do i_qp = 1, nqp

            x1=qp_array_x(i_qp,iconsidered)
            y1=qp_array_y(i_qp,iconsidered)
            if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,iconsidered)
            end if

            if (itestcase.lt.3) then ! linear advection
                flux_term_x=dg_sol(n,iconsidered,x1,y1,z1)*lamx !flux in x-axis of linear advection in 2d (lamx*sol), &
     !& where lamx is the wave speed for x axis, and sol is the solution
                flux_term_y=dg_sol(n,iconsidered,x1,y1,z1)*lamy !flux in y-axis of linear advection in 2d (lamy*sol), &
    ! & where lamy is the wave speed for y axis, and sol is the solution

                if (dimensiona.eq.3)then
                flux_term_z = dg_sol(n,iconsidered,x1,y1,z1)*lamz
                end if

            end if

            if (itestcase.eq.3) then ! euler

                leftv=dg_sol(n,iconsidered,x1,y1,z1)




                    if (dimensiona.eq.2)then
                        call cons2prim(n,leftv,mp_pinfl,gammal)

                        call flux2dx(flux_term_x,leftv)
                        call flux2dy(flux_term_y,leftv)

                    else

                        call cons2prim(n,leftv,mp_pinfl,gammal)
                        call flux3dx(flux_term_x,leftv)
                        call flux3dy(flux_term_y,leftv)
                        call flux3dz(flux_term_z,leftv)

                    end if


            end if


           if (itestcase == 4) then ! ns
                leftv = dg_sol(n,iconsidered,x1,y1,z1)
                leftv_der = dg_sol_der(x1,y1,z1,number_of_dog,iconsidered)

                if( br2_yn == 2) then
                    call dcons2dprim(leftv_der,leftv)
                    leftv_der = leftv_der + sum(rec_br2_local_lift(:,:,:,iconsidered), dim=3)
                end if


                if (dimensiona.eq.2)then
                    call cons2prim(n,leftv,mp_pinfl,gammal)
                    call flux2dx(flux_term_x,leftv)
                    call flux2dy(flux_term_y,leftv)


                     call prim2cons(n,leftv)
                    call flux_visc2d(flux_term_x,flux_term_y,leftv,leftv_der)

                else
                    call cons2prim(n,leftv,mp_pinfl,gammal)
                    call flux3dx(flux_term_x,leftv)
                    call flux3dy(flux_term_y,leftv)
                    call flux3dz(flux_term_z,leftv)
                    call prim2cons(n,leftv)
                    call flux_visc3d(flux_term_x,flux_term_y,flux_term_z,leftv,leftv_der)




                end if
          end if

                    if (dimensiona.eq.2)then


                      if (poly.eq.1)then

                      dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array_qp_weight(i_qp,iconsidered) *(flux_term_x(:)*df2dx(x1,y1,i,iconsidered)+flux_term_y(:)*df2dy(x1,y1,i,iconsidered))

                      end if

                      if (poly.eq.4)then

                      dg_vol_integral(i+1,1:nof_variables) = dg_vol_integral(i+1,1:nof_variables)+ qp_array_qp_weight(i_qp,iconsidered)*(flux_term_x(1:nof_variables)*tl2dx(x1,y1,i,iconsidered)+flux_term_y(1:nof_variables)*tl2dy(x1,y1,i,iconsidered))

                      end if


                    else
                            if (poly.eq.1)then

                            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array_qp_weight(i_qp,iconsidered) *((flux_term_x(:)*dfx(x1,y1,z1,i,iconsidered)) +(flux_term_y(:)*dfy(x1,y1,z1,i,iconsidered))+(flux_term_z(:)*dfz(x1,y1,z1,i,iconsidered)))

                            end if

                            if (poly.eq.2)then

                            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array_qp_weight(i_qp,iconsidered) *((flux_term_x(:)*dlx(x1,y1,z1,i,iconsidered)) +(flux_term_y(:)*dly(x1,y1,z1,i,iconsidered))+(flux_term_z(:)*dlz(x1,y1,z1,i,iconsidered)))

                            end if



                            if (poly.eq.4)then

                            dg_vol_integral(i+1,:) = dg_vol_integral(i+1,:)+ qp_array_qp_weight(i_qp,iconsidered)*((flux_term_x(:)*tl3dx(x1,y1,z1,i,iconsidered)) +(flux_term_y(:)*tl3dy(x1,y1,z1,i,iconsidered))+(flux_term_z(:)*tl3dz(x1,y1,z1,i,iconsidered)))

                            end if

                    end if
         end do
     end do





end function dg_vol_integral








function dg_vol_integral2(n,i)
implicit none
!> @brief
!> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i
real,dimension(1:nof_variables)::dg_vol_integral2
integer::j,k,nqp,i_qp,i_var,icompwrt,iconsidered
real::ph,integ
real::x1,y1,z1,tempdp
integer::number,number_of_dog
real,dimension(1:nof_variables)::dg_sol2
real,dimension(1:55)::basis_temp



iconsidered=i
    icompwrt=-2


    nqp = ielem_itotalpoints(iconsidered)


    number=ielem_iorder(iconsidered)
    number_of_dog = ielem_idegfree(iconsidered)
         dg_vol_integral2(:) = 0.0d0





        do i_qp = 1, nqp
            x1=qp_array_x(i_qp,iconsidered)
            y1=qp_array_y(i_qp,iconsidered)
            if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,iconsidered)
            end if




            if (dimensiona.eq.2)then
            basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
            else
            basis_temp = basis_rec(n, x1, y1, z1,number, iconsidered, number_of_dog,icompwrt)
            end if









            do i_var = 1, nof_variables

                    tempdp=0.0d0
            do k=1,number_of_dog
              tempdp=tempdp+basis_temp(k)*u_c_valdg(1,i_var,k+1,iconsidered)
            end do



                dg_sol2(i_var) = tempdp!dot_product(basis_temp(1:number_of_dog), u_c_valdg(1,i_var,2:number_of_dog+1,iconsidered))
            end do



              dg_vol_integral2 = dg_vol_integral2+ (qp_array_qp_weight(i_qp,iconsidered)*dg_sol2/ielem_totvolume(iconsidered))


         end do


            dg_vol_integral2(:) =u_c_valdg(1,:,1,iconsidered)+dg_vol_integral2(:)

    icompwrt=0




end function dg_vol_integral2



function dg_vol_integral_strong(n,i)
implicit none
!> @brief
!> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i
real,dimension(1:nof_variables)::dg_vol_integral_strong
integer::j,k,nqp,i_qp,i_var,icompwrt,iconsidered
real::ph,integ
real::x1,y1,z1,tempdp
integer::number,number_of_dog
real,dimension(1:nof_variables)::dg_sol2
real,dimension(1:55)::basis_temp



iconsidered=i


    icompwrt=-2


    do j=1,nof_variables
    do k=1,idegfree
    u_cs_valdg(1,j,k+1,iconsidered)=u_c_valdg(1,j,k+1,iconsidered)*modal_filter_strong(k)
    end do
    end do



    nqp = ielem_itotalpoints(iconsidered)






    number=ielem_iorder(iconsidered)
    number_of_dog = ielem_idegfree(iconsidered)
         dg_vol_integral_strong(:) = 0.0d0





        do i_qp = 1, nqp
            x1=qp_array_x(i_qp,iconsidered)
            y1=qp_array_y(i_qp,iconsidered)
            if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,iconsidered)
            end if




            if (dimensiona.eq.2)then

            basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
            else
            basis_temp = basis_rec(n, x1, y1, z1,number, iconsidered, number_of_dog,icompwrt)

            end if



            do i_var = 1, nof_variables
                    tempdp=0.0d0
                do k=1,number_of_dog
                  tempdp=tempdp+basis_temp(k)*u_cs_valdg(1,i_var,k+1,iconsidered)
                end do


                dg_sol2(i_var) = tempdp!dot_product(basis_temp(1:number_of_dog), u_cs_valdg(1,i_var,2:number_of_dog+1,iconsidered))
            end do



              dg_vol_integral_strong= dg_vol_integral_strong+ (qp_array_qp_weight(i_qp,iconsidered)*dg_sol2/ielem_totvolume(iconsidered))


         end do


            dg_vol_integral_strong(:) =u_c_valdg(1,:,1,iconsidered)+dg_vol_integral_strong(:)

    icompwrt=0




end function dg_vol_integral_strong



function dg_vol_integral_weak(n,i)
implicit none
!> @brief
!> calculates the volume integral term in the dg rhs for scalar linear advection with speed = 1
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i
real,dimension(1:nof_variables)::dg_vol_integral_weak
integer::j,k,nqp,i_qp,i_var,icompwrt,iconsidered
real::ph,integ
real::x1,y1,z1,tempdp
integer::number,number_of_dog
real,dimension(1:nof_variables)::dg_sol2
real,dimension(1:55)::basis_temp



iconsidered=i


    icompwrt=-2



    do j=1,nof_variables
    do k=1,idegfree
    u_cw_valdg(1,j,k+1,iconsidered)=u_c_valdg(1,j,k+1,iconsidered)*modal_filter_weak(k)
    end do
    end do


    nqp = ielem_itotalpoints(iconsidered)






    number=ielem_iorder(iconsidered)
    number_of_dog = ielem_idegfree(iconsidered)
         dg_vol_integral_weak(:) = 0.0d0





        do i_qp = 1, nqp
            x1=qp_array_x(i_qp,iconsidered)
            y1=qp_array_y(i_qp,iconsidered)
            if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,iconsidered)
            end if




            if (dimensiona.eq.2)then

            basis_temp = basis_rec2d(n, x1, y1, number, iconsidered, number_of_dog,icompwrt)
            else
            basis_temp = basis_rec(n, x1, y1, z1,number, iconsidered, number_of_dog,icompwrt)

            end if



            do i_var = 1, nof_variables
             tempdp=0.0d0
                do k=1,number_of_dog
                  tempdp=tempdp+basis_temp(k)*u_cw_valdg(1,i_var,k+1,iconsidered)
                end do


                dg_sol2(i_var) = tempdp!dot_product(basis_temp(1:number_of_dog), u_cw_valdg(1,i_var,2:number_of_dog+1,iconsidered))
            end do



              dg_vol_integral_weak = dg_vol_integral_weak+ (qp_array_qp_weight(i_qp,iconsidered)*dg_sol2/ielem_totvolume(iconsidered))


         end do


            dg_vol_integral_weak(:) =u_c_valdg(1,:,1,iconsidered)+dg_vol_integral_weak(:)

    icompwrt=0

end function dg_vol_integral_weak



subroutine reconstruct_dg(n)
implicit none
integer,intent(in)::n
integer::i_face, i_elem, i_qp,iqp,icompwrt
integer::facex,pointx,iconsidered,number_of_dog


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i_elem = 1, xmpielrank(n)
    icompwrt=-2

    do i_face = 1, ielem_ifca(i_elem)
            !somewhere prestored the guassian quadrature points for your sides in another subroutine and you build only the cleft states and volume integral

            if (dimensiona.eq.2)then

            iqp=qp_line_n
            else
                if (ielem_types_faces(i_face,i_elem).eq.5)then
					iqp=qp_quad
				  else
					iqp=qp_triangle

				  end if
            end if

        do i_qp = 1,iqp! qp_line_n

            facex=i_face
            pointx=i_qp
            iconsidered=i_elem
            number_of_dog=ielem_idegfree(i_elem)

            rec_uleft_dg(:, facex, pointx,iconsidered) = dg_solface(n,facex,pointx,iconsidered,number_of_dog)



        end do
    end do

     icompwrt=0

end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine reconstruct_dg


subroutine reconstruct_br2_dg
implicit none
integer::i_face, i_elem, i_qp, iqp,j,k,iconsidered,facex,pointx,number_of_dog,number
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:numberofpoints2)::weights_temp
real,dimension(1:nof_variables,1:dimensiona)::leftv_der
real,dimension(1:nof_variables)::leftv
real::x1,y1,z1




! needed for dg_surf_flux, used in br2_local_lift


#ifdef gpu
!!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i_elem = 1, xmpielrank(n)
    iconsidered = i_elem
    do i_face = 1, ielem_ifca(i_elem)
        ! needed for br2_local_lift
        if (dimensiona.eq.2)then
            iqp=qp_line_n
            weights_temp(1:qp_line_n)=weights_l(1:qp_line_n)
        else
            if (ielem_types_faces(i_face,i_elem).eq.5)then
                iqp=qp_quad
                weights_temp(1:qp_quad)=weights_t(1:qp_quad)
            else
                iqp=qp_triangle
                weights_temp(1:qp_triangle)=weights_t(1:qp_triangle)
            end if
        end if
        number_of_dog=idegfree
        number=ielem_iorder(iconsidered)

        facex = i_face
        rec_br2_local_lift(:,:,i_face,iconsidered) = br2_local_lift(n, iqp,facex,iconsidered,weights_temp)

        do i_qp = 1, iqp
            pointx=i_qp
            if (br2_yn /= 1 .and. itestcase == 4) then
                x1 = rec_qpoints(facex, pointx,1,iconsidered)
                y1 = rec_qpoints(facex, pointx,2,iconsidered)
                if (dimensiona == 3) z1 = rec_qpoints(facex, pointx,3,iconsidered)

                leftv_der = dg_sol_der(x1,y1,z1,number_of_dog,iconsidered)
                call dcons2dprim(leftv_der,leftv)
                rec_br2_aux_var(:,:,i_face,i_qp,iconsidered) = leftv_der + rec_br2_local_lift(:,:,i_face,iconsidered) * br2_damping



                !now copy to other array
                rec_uleftv(1:dims,1:nof_variables-1,i_face,i_qp,iconsidered)=rec_br2_aux_var(1:nof_variables-1,1:dims,i_face,i_qp,iconsidered)


                !end copying
            end if
        end do
    end do








end do
#ifdef gpu
!!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine reconstruct_br2_dg

subroutine get_left_right_states(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer::i,l,ngp
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz

    i = iconsidered

    if (dimensiona == 3) then
        if (ielem_interior(i) == 0) then
            call get_states_interior(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
        else
            call get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
        end if
    else
        if (ielem_interior(i) == 0) then
            call get_states_interior2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
        else
            call get_states_bounds2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
        end if
    end if

end subroutine get_left_right_states


subroutine  get_states_interior(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer::i,l,ngp
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables)::srf_speed
l=facex
ngp=pointx
i=iconsidered


                      if (dg.eq.1) then
                        cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
                        cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                    else !fv

				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state
				      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i)) !right mean flow state

				      end if

					if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					  if (icoupleturb.eq.1)then
					    cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i) !left additional equations flow state
					    cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))!right additional equations flow state
					  else
					    cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
					    cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
					  end if

                      if (turbulence.eq.1)then
					  cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
 					  cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
 					  end if
					end if

					if (rec_mrf(i).eq.1)then

                        srf_speed(2:4)=rec_rotvel(l,ngp,1:3,i)
                        call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
					end if



end subroutine get_states_interior



subroutine  get_states_interior2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer::i,l,ngp
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz

l=facex
ngp=pointx
i=iconsidered


if (dg.eq.1) then
                        cleft(1:nof_variables)= rec_uleft_dg(1:nof_variables, l, ngp,i)
                        cright(1:nof_variables)= rec_uleft_dg(1:nof_variables,ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                        else !fv

                        cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state
                        cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i)) !right mean flow state

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

                      if (turbulence.eq.1)then
                        cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
                        cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)

                      end if
					end if


end subroutine get_states_interior2d



subroutine calculate_interior_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny, &
     & nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:nof_variables)::srf_speed
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables-1,1:dims)::lcvgrad,rcvgrad
real,dimension(turbulenceequations+passivescalar,1:dims)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:dimensiona)::cords
integer::i,l,ngp,ittt,nvar,iex,k
l=facex
ngp=pointx
i=iconsidered



					if (dg.eq.1) then
                        cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
                        cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
						else !fv

				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state
				      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i)) !right mean flow state

				      end if


                      do k=1,dimensiona
				      lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i);
				      rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
				      end do


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


					  do nvar=1,turbulenceequations+passivescalar
					  rcvgrad_t(nvar,1:3)=rec_uleftturbv(1:3,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
					  lcvgrad_t(nvar,1:3)=rec_uleftturbv(1:3,nvar,l,ngp,i)
					  end do


					end if



					if (rec_mrf(i).eq.1)then
						srf_speed(2:4)=rec_rotvel(l,ngp,1:3,i)
						call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
				end if



end subroutine calculate_interior_viscous



subroutine calculate_interior_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:nof_variables)::srf_speed
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables-1,1:dims)::lcvgrad,rcvgrad
real,dimension(turbulenceequations+passivescalar,1:dims)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:dimensiona)::cords
integer::i,l,ngp,ittt,nvar,iex,k
l=facex
ngp=pointx
i=iconsidered




					  if (dg.eq.1) then
                        cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
                        cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))



                        do k=1,dimensiona
				      lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i);
				      rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
				      end do





						else !fv
							cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state




				      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i)) !right mean flow state


				      do k=1,dimensiona
				      lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i);
				      rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
				      end do



						end if



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


					  do nvar=1,turbulenceequations+passivescalar
					  rcvgrad_t(nvar,1:2)=rec_uleftturbv(1:2,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
					  lcvgrad_t(nvar,1:2)=rec_uleftturbv(1:2,nvar,l,ngp,i)
					  end do


					end if





end subroutine calculate_interior_viscous2d


subroutine calculate_bounded_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny, &
     & nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:nof_variables)::srf_speed
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables-1,1:dims)::lcvgrad,rcvgrad
real,dimension(turbulenceequations+passivescalar,1:dims)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:dimensiona)::cords
real,dimension(1:3):: gl,nnt,g
real::g_n
integer :: d
integer::i,l,ngp,ittt,nvar,iex,kk,n_node,k,nfx,lfx,rowfx
integer::ibfc

l=facex
ngp=pointx
i=iconsidered




						if (dg == 1) then
									cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
									else
								cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
								end if



!   cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)	!left mean flow state



				      do k=1,dimensiona
				      lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i);
				      end do






					  if (rec_mrf(i).eq.1)then
						srf_speed(2:4)=rec_rotvel(l,ngp,1:3,i)
						call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
				end if

					 if ((turbulence.eq.1).or.(passivescalar.gt.0))then
						if (icoupleturb.eq.1)then
							cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i)
						else
							cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
						end if
						do nvar=1,turbulenceequations+passivescalar
						lcvgrad_t(nvar,1:3)=rec_uleftturbv(1:3,nvar,l,ngp,i)
						end do
					end if


					    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then

								  if (dg == 1) then
                                    cright(1:nof_variables)= rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                    else
                                    cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                    end if





								   do k=1,dimensiona
                                    rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                end do



								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if

									do nvar=1,turbulenceequations+passivescalar
									rcvgrad_t(nvar,1:3)=rec_uleftturbv(1:3,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									end do


								    end if
									if(per_rot.eq.1)then
                                        cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)

                                        do k=1,dimensiona
                                    rcvgrad(1:nof_variables-1,k)=rotate_per_1(rcvgrad(1:nof_variables-1,k),ibound_icode(ielem_ibounds(l,i)),angle_per)
                                end do


                                        rcvgrad_t(1,1:3)=rotate_per_1(rcvgrad_t(1,1:3),ibound_icode(ielem_ibounds(l,i)),angle_per)
                                    end if



								  else
								  !not periodic ones in my cpu


								  call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)
								    cords(1:3)=zero

								    if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if


								    cords(1:3)=cordinates3(n,nodes_list,n_node)

								    poy(1)=cords(2)
								    pox(1)=cords(1)
								    poz(1)=cords(3)

								    leftv(1:nof_variables)=cleft(1:nof_variables)
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    cright(1:nof_variables)=rightv(1:nof_variables)





				  				  	if (b_code.eq.4)then

                                            if (thermal.ne.1)then

                                                  do k=1,dimensiona
                                                  rcvgrad(1:dimensiona,k)=lcvgrad(1:dimensiona,k)
                                                  end do


                                                ! normal vector (global coords)
                                                nnt(1) = nx
                                                nnt(2) = ny
                                                if (dimensiona == 3) then
                                                  nnt(3) = nz
                                                else
                                                  nnt(3) = 0.0
                                                end if
                                                do iex=dimensiona+1,nof_variables-1
                                                ! load left gradient of temperature (global components)
                                                do d = 1, dimensiona
                                                  g(d) = lcvgrad(iex, d)
                                                end do
                                                if (dimensiona == 2) g(3) = 0.0

                                                ! normal component g_n = g · nnt
                                                g_n = 0.0
                                                do d = 1, dimensiona
                                                  g_n = g_n + g(d)*nnt(d)
                                                end do

                                                ! remove normal component: g_t = g - g_n * nnt
                                                do d = 1, dimensiona
                                                  g(d) = g(d) - g_n*nnt(d)
                                                end do

                                                ! write back into lcvgrad (this is the gradient used by the interior side)
                                                do d = 1, dimensiona
                                                  lcvgrad(iex, d) = g(d)
                                                end do

                                                ! for consistency with your averaging, set rcvgrad equal to the same wall gradient
                                                do d = 1, dimensiona
                                                  rcvgrad(iex, d) = g(d)
                                                end do
                                                end do





                                              if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                              rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                              end if

                                            else


                                                  do k=1,dimensiona
                                                  rcvgrad(1:nof_variables-1,k)=lcvgrad(1:nof_variables-1,k)
                                                  end do
                                                  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                                  rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                  end if

                                            end if




				  				    else

				  				     if ((b_code.ne.5).and.(b_code.gt.0))then

                                            ! load normal
                                            nnt(1) = nx
                                            nnt(2) = ny
                                            if (dimensiona == 3) then
                                              nnt(3) = nz
                                            else
                                              nnt(3) = 0.0d0
                                            end if

                                            do iex=1,nof_variables-1
                                            ! load left gradient of this variable

                                            do d = 1, dimensiona
                                              gl(d) = lcvgrad(iex, d)
                                            end do
                                            if (dimensiona == 2) gl(3) = 0.0d0

                                            ! compute normal component g_n = gl · nnt
                                            g_n = 0.0d0
                                            do d = 1, dimensiona
                                              g_n = g_n + gl(d)*nnt(d)
                                            end do

                                            ! reflect the normal component: gr = gl - 2*g_n*nnt
                                            do d = 1, dimensiona
                                              rcvgrad(iex, d) = gl(d) - 2.0*g_n*nnt(d)
                                            end do

                                            end do

                                            if ((turbulence.eq.1).or.(passivescalar.gt.0))then

                                            do iex=1,turbulenceequations+passivescalar
                                            ! load left gradient of this variable

                                            do d = 1, dimensiona
                                              gl(d) = lcvgrad_t(iex, d)
                                            end do
                                            if (dimensiona == 2) gl(3) = 0.0d0

                                            ! compute normal component g_n = gl · nnt
                                            g_n = 0.0d0
                                            do d = 1, dimensiona
                                              g_n = g_n + gl(d)*nnt(d)
                                            end do

                                            ! reflect the normal component: gr = gl - 2*g_n*nnt
                                            do d = 1, dimensiona
                                              rcvgrad_t(iex, d) = gl(d) - 2.0*g_n*nnt(d)
                                            end do

                                            end do

                                            end if
                                            end if

				  				    end if

!
!
								  end if
							else
							      if (dg == 1) then
                                cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                else
							      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
							      end if


							      do k=1,dimensiona
                                    rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                end do



								  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
									do nvar=1,turbulenceequations+passivescalar
									rcvgrad_t(nvar,1:3)=rec_uleftturbv(1:3,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									end do
								    end if




							end if
					    else	!in other cpus they can only be periodic or mpi neighbours






							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then	!periodic in other
									  if (dg == 1) then
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,&
!      & i),1:nof_variables)
                                    nfx  = ielem_ineighn(l,i)
                                    lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                    rowfx = bound_offset(nfx) + lfx - 1
                                    cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)




								else

! 									  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & 1:nof_variables)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                                end if

									  ittt=0
									  do iex=1,nof_variables-1
										do nvar=1,dims
										      ittt=ittt+1
										      if (dg.eq.1)then
!  										      rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
 										      nfx  = ielem_ineighn(l,i)
 										      lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
 										      rowfx = bound_offset(nfx) + lfx - 1
 										      rcvgrad(iex,nvar) = boundhir_dg(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
!
										      else

!                                                 rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l, &
!      & i))%facesol(ielem_qface(l,ngp,i),nof_variables+turbulenceequations+passivescalar+ittt)
                                                nfx  = ielem_ineighn(l,i)
                                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                                rowfx = bound_offset(nfx) + lfx - 1
                                                rcvgrad(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)


										      end if


										end do
									  end do






								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
								    end if



									  do iex=1,turbulenceequations+passivescalar
										do nvar=1,dims
										      ittt=ittt+1
! 									  rcvgrad_t(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad_t(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
										end do
									  end do

									  if(per_rot.eq.1)then


									   do k=1,dimensiona
                                    rcvgrad(1:nof_variables-1,k)=rotate_per_1(rcvgrad(1:nof_variables-1,k),ibound_icode(ielem_ibounds(l,i)),angle_per)
                                end do



										cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)

										  rcvgrad_t(1,1:3)=rotate_per_1(rcvgrad_t(1,1:3),ibound_icode(ielem_ibounds(l,i)),angle_per)
										end if



								end if
							else

								   if (dg == 1) then
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,&
!      & i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                else
! 								  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
								  end if
								  ittt=0
									  do iex=1,nof_variables-1
										do nvar=1,dims
										      ittt=ittt+1
										      if (dg.eq.1)then
!  										      rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
 										      nfx  = ielem_ineighn(l,i)
 										      lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
 										      rowfx = bound_offset(nfx) + lfx - 1
 										      rcvgrad(iex,nvar) = boundhir_dg(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
!
										      else




! 									  rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)


										      end if
										end do
									  end do



								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if

									 do iex=1,turbulenceequations+passivescalar
										do nvar=1,dims
										      ittt=ittt+1
! 									  rcvgrad_t(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad_t(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
										end do
									  end do


								    end if

!
							end if
					    end if



end subroutine calculate_bounded_viscous



subroutine calculate_bounded_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny, &
     & nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:nof_variables)::srf_speed,tempx_l,rtempx_l
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables-1,1:dims)::lcvgrad,rcvgrad
real,dimension(turbulenceequations+passivescalar,1:dims)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:dimensiona)::cords
real,dimension(1:3):: gl,nnt,g,grad
real::g_n
integer :: d
integer::i,l,ngp,ittt,nvar,iex,n_node,k,nfx,lfx,rowfx
integer::ibfc
l=facex
ngp=pointx
i=iconsidered

if (dg == 1) then
	cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
	else
cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
end if

                      do k=1,dimensiona
				      lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i)
				      end do


					 if ((turbulence.eq.1).or.(passivescalar.gt.0))then
						if (icoupleturb.eq.1)then
							cturbl(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,l,ngp,i)
						else
							cturbl(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
						end if
						do nvar=1,turbulenceequations+passivescalar
						lcvgrad_t(nvar,1:2)=rec_uleftturbv(1:2,nvar,l,ngp,i)
						end do
					end if


					    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								  !cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))

									if (dg == 1) then
										cright(1:nof_variables)= rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
										else
										cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
										end if





								  do k=1,dimensiona
                                  rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineigh(l,i),ngp,ielem_ineigh(l,i))
                                  end do



								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if

									do nvar=1,turbulenceequations+passivescalar
									rcvgrad_t(nvar,1:2)=rec_uleftturbv(1:2,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									end do


								    end if




								  else
								  !not periodic ones in my cpu


								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)


                                            n_node=2


								    cords(1:2)=zero
								    cords(1:2)=cordinates2(n,nodes_list,n_node)

								    poy(1)=cords(2)
								    pox(1)=cords(1)


								    leftv(1:nof_variables)=cleft(1:nof_variables)
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    cright(1:nof_variables)=rightv(1:nof_variables)

								   if (b_code.eq.4)then

                                                          if (thermal.ne.1)then

                                                                        do k=1,dimensiona
                                                                        rcvgrad(1:dimensiona,k)=lcvgrad(1:dimensiona,k)
                                                                        end do


                                                                      do iex=dimensiona+1,nof_variables-1


                                                                      tempx_l=0.0d0
                                                                      rtempx_l=0.0d0
                                                                      tempx_l(2)=lcvgrad(iex,1)
                                                                      tempx_l(3)=lcvgrad(iex,2)
                                                                      call rotatef2d(n,rtempx_l,tempx_l,angle1,angle2)
                                                                      rtempx_l(2)=-rtempx_l(2)
                                                                      call rotateb2d(n,tempx_l,rtempx_l,angle1,angle2)
                                                                      rcvgrad(iex,1)=tempx_l(2)
                                                                      rcvgrad(iex,2)=tempx_l(3)

                                                                      end do





                                                                    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                                                    rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                                    end if

                                                          else    !thermal


                                                                do k=1,dimensiona
                                                                rcvgrad(1:nof_variables-1,k)=lcvgrad(1:nof_variables-1,k)
                                                                end do
                                                                if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                                                rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                                end if


                                                                if (catalytic_wall.eq.0)then

                                                                     do iex = dimensiona+3,nof_variables-1     ! species indices
                                                                      tempx_l=0.0d0
                                                                      rtempx_l=0.0d0
                                                                      tempx_l(2)=lcvgrad(iex,1)
                                                                      tempx_l(3)=lcvgrad(iex,2)
                                                                      call rotatef2d(n,rtempx_l,tempx_l,angle1,angle2)
                                                                      rtempx_l(2)=-rtempx_l(2)
                                                                      call rotateb2d(n,tempx_l,rtempx_l,angle1,angle2)
                                                                      rcvgrad(iex,1)=tempx_l(2)
                                                                      rcvgrad(iex,2)=tempx_l(3)





                                                                    end do

                                                                end if





                                                  end if !thermal

				  				    else

                                                      if ((b_code.ne.5).and.(b_code.gt.0))then

                                                                    do k=1,dimensiona
                                                                    rcvgrad(1:nof_variables-1,k)=lcvgrad(1:nof_variables-1,k)
                                                                    end do
                                                                    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                                                    rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                                    end if






                                                      end if
                                    end if


                                end if
			else











								if (dg == 1) then
									cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
									else
									  cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									  end if






							      do k=1,dimensiona
                                  rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                  end do

								  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
									do nvar=1,turbulenceequations+passivescalar
									rcvgrad_t(nvar,1:2)=rec_uleftturbv(1:2,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									end do
								    end if




							end if
					    else	!in other cpus they can only be periodic or mpi neighbours






							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
									 ! cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)

									if (dg == 1) then
! 										cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i), &
!      & 1:nof_variables)
										nfx  = ielem_ineighn(l,i)
										lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
										rowfx = bound_offset(nfx) + lfx - 1
										cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
									else

! 										  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & 1:nof_variables)
										  nfx  = ielem_ineighn(l,i)
										  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
										  rowfx = bound_offset(nfx) + lfx - 1
										  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
									end if


									  ittt=0
									  do iex=1,nof_variables-1
										do nvar=1,dims
										      ittt=ittt+1
										      if (dg.eq.1)then
!  										      rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
 										      nfx  = ielem_ineighn(l,i)
 										      lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
 										      rowfx = bound_offset(nfx) + lfx - 1
 										      rcvgrad(iex,nvar) = boundhir_dg(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
!
										      else




! 									  rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)


										      end if
										end do
									  end do






								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if




									  do iex=1,turbulenceequations+passivescalar
										do nvar=1,dims
										      ittt=ittt+1
! 									  rcvgrad_t(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad_t(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
										end do
									  end do

								end if



								end if
							else

								  !cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
								if (dg == 1) then
! 									cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i), &
!      & 1:nof_variables)
									nfx  = ielem_ineighn(l,i)
									lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									rowfx = bound_offset(nfx) + lfx - 1
									cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
									else
! 									  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & 1:nof_variables)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
									  end if


								ittt=0
									  do iex=1,nof_variables-1
										do nvar=1,dims
										      ittt=ittt+1
										      if (dg.eq.1)then
!  										      rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
 										      nfx  = ielem_ineighn(l,i)
 										      lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
 										      rowfx = bound_offset(nfx) + lfx - 1
 										      rcvgrad(iex,nvar) = boundhir_dg(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
!
										      else



! 									  rcvgrad(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)


										      end if
										end do
									  end do



								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if

									 do iex=1,turbulenceequations+passivescalar
										do nvar=1,dims
										      ittt=ittt+1
! 									  rcvgrad_t(iex,nvar)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & nof_variables+turbulenceequations+passivescalar+ittt)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  rcvgrad_t(iex,nvar) = boundhir(rowfx,nof_variables+turbulenceequations+passivescalar+ittt)
										end do
									  end do


								    end if

!
							end if
					    end if






end subroutine calculate_bounded_viscous2d




subroutine  get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer::i,l,ngp,n_node
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,srf_speedrot
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables)::srf_speed
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:dimensiona)::cords
integer::ibfc,nfx,lfx,rowfx
l=facex
ngp=pointx
i=iconsidered



						if (dg == 1) then
									cleft(1:nof_variables) = rec_uleft_dg(1:nof_variables, l, ngp,i)
									else
								cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
								end if

								if (rec_mrf(i).eq.1)then
									srf_speed(2:4)=rec_rotvel(l,ngp,1:3,i)
									call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
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
								  !if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
									if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then

								 if (dg == 1) then
                                    cright(1:nof_variables)= rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                    else
                                    cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                    end if

									if(per_rot.eq.1)then
										cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
									  end if

								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
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

								    leftv(1:nof_variables)=cleft(1:nof_variables)
								    b_code=ibound_icode(ielem_ibounds(l,i))



								    call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    cright(1:nof_variables)=rightv(1:nof_variables)


								  end if
							else


                                   if (dg == 1) then
                                cright(1:nof_variables) = rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                else
							      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
							      end if

								  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
								    end if




							end if
					    else	!in other cpus they can only be periodic or mpi neighbours






							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								!if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
									if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
								 if (dg == 1) then
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,&
!      & i),1:nof_variables)
                                    nfx  = ielem_ineighn(l,i)
                                    lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                    rowfx = bound_offset(nfx) + lfx - 1
                                    cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
								else

! 									  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i), &
!      & 1:nof_variables)
									  nfx  = ielem_ineighn(l,i)
									  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									  rowfx = bound_offset(nfx) + lfx - 1
									  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                                end if
								if(per_rot.eq.1)then
									cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
								  end if

								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
								    end if



								end if
							else
                                   if (dg == 1) then
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol_dg(ielem_qface(l,ngp,&
!      & i),1:nof_variables)
                                nfx  = ielem_ineighn(l,i)
                                lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                rowfx = bound_offset(nfx) + lfx - 1
                                cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                else
! 								  cright(1:nof_variables)=iexboundhir(ielem_ineighn(l,i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
								  nfx  = ielem_ineighn(l,i)
								  lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
								  rowfx = bound_offset(nfx) + lfx - 1
								  cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
								  end if

!
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if
								    end if

!
							end if
					    end if

                         if (turbulence.eq.1)then
					  cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
 					  cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
 					  end if




end subroutine get_states_bounds



subroutine  get_states_bounds2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered, facex, pointx,n
integer::i,l,ngp,n_node
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot,srf_speedrot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv
real,dimension(1:nof_variables),intent(inout)::rightv
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables)::srf_speed
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:8,1:dimensiona)::nodes_list
real,dimension(1:dimensiona)::cords
real,dimension(1:6,1:4,1:dimensiona)::elem_listd
integer::ikas
integer::ibfc,nfx,lfx,rowfx
l=facex
ngp=pointx
i=iconsidered

if (dg == 1) then
                        cleft(1:nof_variables)= rec_uleft_dg(1:nof_variables, l, ngp,i)
                    else

				      cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)

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
                                        if (dg == 1) then
                                            cright(1:nof_variables)= rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
                                        else !fv
                                            cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                        end if
								    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                        if (icoupleturb.eq.1)then
                                        cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                        else
                                        cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
                                        end if
								    end if





								  ikas=1

								  else
								  !not periodic ones in my cpu


								 call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)



                                            n_node=2


								    cords(1:2)=zero
								    cords(1:2)=cordinates2(n,nodes_list,n_node)

								    poy(1)=cords(2)
								    pox(1)=cords(1)


								    leftv(1:nof_variables)=cleft(1:nof_variables)
!
								    b_code=ibound_icode(ielem_ibounds(l,i))
								    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
								    cright(1:nof_variables)=rightv(1:nof_variables)
                                    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
								    cturbr(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)

								    end if

!
				  				  	 ikas=2

								  end if
							else
                                    if (dg == 1) then
                                        cright(1:nof_variables)= rec_uleft_dg(1:nof_variables, ielem_ineighn(l,i), ngp,ielem_ineigh(l,i))
        !
                                    else !fv
                                        cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
                                    end if

								  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
									   cturbr(1:turbulenceequations+passivescalar)=rec_uleftturb(1:turbulenceequations+passivescalar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
									else
									 cturbr(1:turbulenceequations+passivescalar)=u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
									end if
								    end if


							       ikas=3

							end if
					    else	!in other cpus they can only be periodic or mpi neighbours






							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                        if (dg == 1) then
!                                             cright(1:nof_variables)=iexboundhir(ielem_ineighn(l, &
!      & i))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                                            nfx  = ielem_ineighn(l,i)
                                            lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                            rowfx = bound_offset(nfx) + lfx - 1
                                            cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                        else
!                                             cright(1:nof_variables)=iexboundhir(ielem_ineighn(l, &
!      & i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                                            nfx  = ielem_ineighn(l,i)
                                            lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                            rowfx = bound_offset(nfx) + lfx - 1
                                            cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
 									   end if
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									 nfx  = ielem_ineighn(l,i)
									 lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									 rowfx = bound_offset(nfx) + lfx - 1
									 cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									end if




								    end if



								end if
							else
                                    if (dg == 1) then
!                                         cright(1:nof_variables)=iexboundhir(ielem_ineighn(l, &
!      & i))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                                        nfx  = ielem_ineighn(l,i)
                                        lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                        rowfx = bound_offset(nfx) + lfx - 1
                                        cright(1:nof_variables) = boundhir_dg(rowfx,1:nof_variables)
                                    else
!                                         cright(1:nof_variables)=iexboundhir(ielem_ineighn(l, &
!      & i))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                                        nfx  = ielem_ineighn(l,i)
                                        lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
                                        rowfx = bound_offset(nfx) + lfx - 1
                                        cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                                    end if
								   if ((turbulence.eq.1).or.(passivescalar.gt.0))then
									if (icoupleturb.eq.1)then
! 									   cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									   nfx  = ielem_ineighn(l,i)
									   lfx = ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
									   rowfx = bound_offset(nfx) + lfx - 1
									   cturbr(1:turbulenceequations+passivescalar) = boundhir(rowfx,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
									else
! 									 cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(l,i))%facesol (ielem_qface(l,ngp,&
!      & i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
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


                          if (turbulence.eq.1)then
					  cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbl(1:turbulenceequations+passivescalar)
 					  cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
 					  end if





end subroutine get_states_bounds2d
















subroutine allocate_dg
!> @brief
!> allocates the gaussian quadrature volume points
    implicit none
    integer::i, k, i_qp, n_qp, i_face,nnd,iqp,idummy,loopc,kmaxe
    real,dimension(1:idegfree+1)::tempint
    real::tempf,tempx







kmaxe=xmpielrank(n)


if (dimensiona.eq.3)then
allocate(qp_array_x(numberofpoints,1:kmaxe))
allocate(qp_array_y(numberofpoints,1:kmaxe))
allocate(qp_array_z(numberofpoints,1:kmaxe))
allocate(qp_array_qp_weight(numberofpoints,1:kmaxe))
else
allocate(qp_array_x(numberofpoints,1:kmaxe))
allocate(qp_array_y(numberofpoints,1:kmaxe))
allocate(qp_array_z(numberofpoints,1:kmaxe))
allocate(qp_array_qp_weight(numberofpoints,1:kmaxe))
end if



      if (dg.eq.1)then
     if((itestcase == 4).or.((multispecies.eq.1).and.(br2_yn.eq.1))) then !ns
           allocate(rec_br2_aux_var(1:nof_variables, 1:dimensiona, max_faces, numberofpoints2,1:kmaxe))
          allocate(rec_br2_local_lift(1:nof_variables, 1:dimensiona, max_faces,1:kmaxe))
            rec_br2_aux_var= zero
            rec_br2_local_lift= zero

        end if
      end if





end subroutine allocate_dg







subroutine prestore_dg1(iconsidered)
!> @brief
!> prestores ielem_delta_xyz(i), qp_array, surf_qpoints, mass matrix
    implicit none
    integer,intent(in)::iconsidered
    integer::i, k, i_qp, n_qp, i_face,nnd,iqp,idummy,loopc,icompwrt,eltype,elem_dec,ixx,number
    integer::count_1
    real,dimension(1:idegfree+1)::tempint
    real,dimension(1:8,1:dimensiona)::nodes_list
    real,dimension(1:8,1:dimensiona)::vext
    real,dimension(1:dimensiona)::cords
    real,dimension(1:dimensiona,1:numberofpoints)::qpoints
    real,dimension(1:numberofpoints)::wequa3d
    real,dimension(1:6,1:4,1:dimensiona)::elem_listd
    real::voltemp,x1,y1,z1
    real::tempf,tempx


    number=ielem_iorder(iconsidered)

        if (dg.eq.1)then
        icompwrt=-2
        end if


        i=iconsidered
        !store volume quadrature points
        eltype=ielem_ishape(i)
        elem_dec=ielem_vdec(i)
        do k = 1,ielem_nonodes(i)
            nodes_list(k,1:dims)=dinoder(ielem_nodes(k,i))%cord(1:dims)
            vext(k,1:dims)=nodes_list(k,1:dims)
        end do


        if (dimensiona.eq.2)then

        call decompose2(n,eltype,nodes_list,elem_listd)
        else
        call decompose3(n,eltype,nodes_list,elem_listd)

        end if


        select case(ielem_ishape(i))


        case(1,2,3,4) !hexa
        count_1=0

             do k=1,elem_dec
                 vext(1:4,1:3)=elem_listd(k,1:4,1:3)

                call quadraturetetra(n,igqrules,vext,qpoints,wequa3d)

                voltemp=tetravolume(n,vext)

                do i_qp = 1, qp_tetra
                    count_1=count_1+1
                     qp_array_x(count_1,i) = (qpoints(1,i_qp)- ielem_xxc(i))
                    qp_array_y(count_1,i) = (qpoints(2,i_qp)- ielem_yyc(i))
                    qp_array_z(count_1,i) = (qpoints(3,i_qp)- ielem_zzc(i))

                    qp_array_qp_weight(count_1,i) = wequa3d(i_qp) * voltemp

                end do
            end do








        case(5)
            count_1=0

             do k=1,elem_dec
                 vext(1:3,1:2)=elem_listd(k,1:3,1:2)

                call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)

                voltemp=trianglevolume(n,vext)

                do i_qp = 1, qp_triangle
                    count_1=count_1+1

                    qp_array_x(count_1,i) = (qpoints(1,i_qp)- ielem_xxc(i))
                    qp_array_y(count_1,i) = (qpoints(2,i_qp)- ielem_yyc(i))



                    qp_array_qp_weight(count_1,i) = wequa3d(i_qp) * voltemp


                end do

            end do



        case(6)
            call quadraturetriangle(n,igqrules,vext,qpoints,wequa3d)
            n_qp = qp_triangle
            voltemp=trianglevolume(n,vext)

            do i_qp = 1, n_qp


                qp_array_x(i_qp,i) = (qpoints(1,i_qp) - ielem_xxc(i))
                    qp_array_y(i_qp,i) = (qpoints(2,i_qp) - ielem_yyc(i))



                    qp_array_qp_weight(i_qp,i) = wequa3d(i_qp) * voltemp

            end do


        end select





        if (dg.eq.1)then
        !store volume quadrature points
        integ_basis_dg_value(1:idegfree,iconsidered)=zero
        end if
        tempint=zero

            n_qp = ielem_itotalpoints(iconsidered)





                do i_qp = 1, n_qp
                    x1=qp_array_x(i_qp,i); y1=qp_array_y(i_qp,i)


            if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,i)

            end if
!
                    ixx=iconsidered
                    if (dg.eq.1)then
                    if (dimensiona.eq.2)then
                    number=ielem_iorder(iconsidered)
                    tempint(1:idegfree)=tempint(1:idegfree)+(basis_rec2d(n,x1,y1,number,ixx,idegfree,icompwrt)*qp_array_qp_weight(i_qp,i))

                    else
                    number=ielem_iorder(iconsidered)
                    tempint(1:idegfree)=tempint(1:idegfree)+(basis_rec(n,x1,y1,z1,number,ixx,idegfree,icompwrt)*qp_array_qp_weight(i_qp,i))

                    end if
                    end if
                end do
!
                if (dg.eq.1)then
                integ_basis_dg_value(1:idegfree,iconsidered)=tempint(1:idegfree)
                end if


        icompwrt=0



end subroutine










subroutine build_mass_matrix(n)
!> @brief
!> assembles the mass matrix
!> requires: globals: ielem, qp_quad, qp_triangle, mass_matrix
    implicit none
    integer,intent(in)::n
    integer::i_elem, i_qp, n_qp, i_dof, j_dof, kmaxe,icompwrt,number_of_dog,iconsidered,ixx,number
    real::integ_test,integ_sm1,integ_sm2,phx1,phx2,kron,maxs,mins,higher,phx,integ_mm,x1,y1,z1
    real,allocatable,dimension(:,:)::totalmm,invmm
    real,allocatable,dimension(:)::basis_vector
    kmaxe = xmpielrank(n)


    allocate(totalmm(1:num_dg_dofs,1:num_dg_dofs),invmm(1:num_dg_dofs,1:num_dg_dofs),basis_vector(1:idegfree))


    !$omp do
    do i_elem = 1, kmaxe
          totalmm(:,:) = zero
           invmm(:,:) = zero
            icompwrt=-2
           n_qp=ielem_itotalpoints(i_elem)




        iconsidered=i_elem
    number_of_dog = ielem_idegfree(i_elem)




    do i_dof = 1, num_dg_dofs
        do j_dof = 1, num_dg_dofs
            integ_mm = zero


                do i_qp = 1, n_qp
                    ixx = i_elem;


                    x1=qp_array_x(i_qp,i_elem);
                    y1=qp_array_y(i_qp,i_elem)
                    if (dimensiona.eq.3)then
                    z1=qp_array_z(i_qp,i_elem);
                    end if



             kron=1.0d0


             if (i_dof.ne.j_dof)then
             kron=1.0d0
             end if

                    if (dimensiona.eq.2)then
                    number=ielem_iorder(iconsidered)
                    basis_vector(1:idegfree) = basis_rec2d(n,x1,y1,number,ixx,number_of_dog,icompwrt)
                    else
                    number=ielem_iorder(iconsidered)
                    basis_vector(1:idegfree) = basis_rec(n,x1,y1,z1,number,ixx,number_of_dog,icompwrt)


                    end if

                        if (i_dof==1.and.j_dof==1)then
                            phx = 1.0d0


                        else if (i_dof==1.and.j_dof/=1)then
                            phx = basis_vector(j_dof-1)


                        else if (i_dof/=1.and.j_dof==1)then
                            phx = basis_vector(i_dof-1)



                        else if (i_dof/=1.and.j_dof/=1)then
                            phx = basis_vector(i_dof-1)*basis_vector(j_dof-1)

                        end if



                    integ_mm = integ_mm + phx*qp_array_qp_weight(i_qp,i_elem)*kron




                end do




        totalmm(i_dof, j_dof) = totalmm(i_dof, j_dof) + integ_mm


        end do
    end do



        call compmassinv(totalmm,invmm)
!
        m_1_val(:,:,i_elem)=invmm(:,:)
!



         icompwrt=0







end do
!$omp end do



    deallocate(totalmm,invmm,basis_vector)














end subroutine build_mass_matrix

subroutine compmassinv(totalmm,invmm)
!calculate the inverse of the input matrix with gauss-jordan elimination
implicit none
integer :: i,j,k,l,m,irow,num_dofs
real:: big,dum
real,allocatable,dimension(:,:)::a,b
real,allocatable,dimension(:,:),intent(in)::totalmm
real,allocatable,dimension(:,:),intent(inout)::invmm
allocate(a(1:num_dg_dofs,1:num_dg_dofs),b(1:num_dg_dofs,1:num_dg_dofs))
num_dofs=num_dg_dofs


a(:,:)=totalmm(:,:)
b(:,:)=zero




do i = 1,num_dofs
    do j = 1,num_dofs
        b(i,j) = 0.0d0
    end do
    b(i,i) = 1.0d0
end do

do i = 1,num_dofs
   big = a(i,i)
   do j = i,num_dofs
     if (a(j,i).gt.big) then
       big = a(j,i)
       irow = j
     end if
   end do
   ! interchange lines i with irow for both a() and b() matrices
   if (big.gt.a(i,i)) then
     do k = 1,num_dofs
       dum = a(i,k)                      ! matrix a()
       a(i,k) = a(irow,k)
       a(irow,k) = dum
       dum = b(i,k)                 ! matrix b()
       b(i,k) = b(irow,k)
       b(irow,k) = dum
     end do
   end if
   ! divide all entries in line i from a(i,j) by the value a(i,i);
   ! same operation for the identity matrix
   dum = a(i,i)
   do j = 1,num_dofs
     a(i,j) = a(i,j)/dum
     b(i,j) = b(i,j)/dum
   end do
   ! make zero all entries in the column a(j,i); same operation for indent()
   do j = i+1,num_dofs
     dum = a(j,i)
     do k = 1,num_dofs
       a(j,k) = a(j,k) - dum*a(i,k)
       b(j,k) = b(j,k) - dum*b(i,k)

     end do
   end do
end do

 do i = 1,num_dofs-1
   do j = i+1,num_dofs
     dum = a(i,j)
     do l = 1,num_dofs
       a(i,l) = a(i,l)-dum*a(j,l)
       b(i,l) = b(i,l)-dum*b(j,l)
     end do
   end do
 end do





 invmm(:,:)=b(:,:)
 deallocate(a,b)

 end subroutine compmassinv




end module dg_functions
