module dg_functions

use basis
use declaration
use derivatives
use lapck
use flow_operations, only: boundarys_ideal, boundarys2d_ideal
use transform
implicit none


contains


subroutine dg_sol_fill(n,iconsidered,x1,y1,z1,dg_sol_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> @brief
!> fills the dg solution at a given point without returning an array temporary
integer,intent(in)::n,iconsidered
real,intent(in)::x1,y1,z1
real,dimension(1:nof_variables),intent(out)::dg_sol_out
integer::i_var,icompwrt,number_of_dog,number,k
real::tempdp,basis_val

icompwrt=-2
number=ielem_iorder(iconsidered)
number_of_dog=ielem_idegfree(iconsidered)

do i_var = 1, nof_variables
    tempdp=0.0d0
    do k=1,number_of_dog
        if (dimensiona.eq.2)then
            basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
        else
            basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
        end if
        tempdp=tempdp+basis_val*u_c_valdg(1,i_var,k+1,iconsidered)
    end do
    dg_sol_out(i_var) = u_c_valdg(1,i_var,1,iconsidered)+tempdp
end do

icompwrt=0

end subroutine dg_sol_fill


subroutine dg_solface_fill(n,facex,pointx,iconsidered,number_of_dog,dg_solface_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> @brief
!> fills the dg solution at a surface point without returning an array temporary
integer,intent(in)::n,facex,pointx,iconsidered,number_of_dog
real,dimension(1:nof_variables),intent(out)::dg_solface_out
integer::i_var,icompwrt,k,number
real::x1,y1,z1,tempdp,basis_val

icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
    z1= rec_qpoints(facex,pointx,3,iconsidered)
end if
number=ielem_iorder(iconsidered)

do i_var = 1, nof_variables
    tempdp=0.0d0
    do k=1,number_of_dog
        if (dimensiona.eq.2)then
            basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
        else
            basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
        end if
        tempdp=tempdp+basis_val*u_c_valdg(1,i_var,k+1,iconsidered)
    end do
    dg_solface_out(i_var) = u_c_valdg(1,i_var,1,iconsidered)+tempdp
end do

icompwrt=0

end subroutine dg_solface_fill


subroutine dg_solface_store(n,facex,pointx,iconsidered,number_of_dog)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,facex,pointx,iconsidered,number_of_dog
integer::i_var,icompwrt,k
real::x1,y1,z1,tempdp,basis_val
integer::number

icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
    z1= rec_qpoints(facex,pointx,3,iconsidered)
end if
number=ielem_iorder(iconsidered)

do i_var = 1, nof_variables
    tempdp=0.0d0
    do k=1,number_of_dog
        if (dimensiona.eq.2)then
            basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
        else
            basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
        end if
        tempdp=tempdp+basis_val*u_c_valdg(1,i_var,k+1,iconsidered)
    end do
    rec_uleft_dg(i_var, facex, pointx, iconsidered) = u_c_valdg(1,i_var,1,iconsidered)+tempdp
end do

icompwrt=0

end subroutine dg_solface_store


subroutine dg_sol_der_fill(x1,y1,z1,number_of_dog,iconsidered,dg_sol_der_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> @brief
!> fills the derivative of the dg solution without returning an array temporary
integer,intent(in)::number_of_dog,iconsidered
real,intent(in)::x1,y1,z1
real,dimension(1:nof_variables,1:dimensiona),intent(out)::dg_sol_der_out
integer::i_var, i_dim,icompwrt,number,k
real::tempdp,basis_val

icompwrt=-2
number=ielem_iorder(iconsidered)

do i_var = 1, nof_variables
    do i_dim = 1, dimensiona
        tempdp=0.0d0
        do k=1,number_of_dog
            if (dimensiona.eq.2)then
                basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
            else
                basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
            end if
            if (br2_yn == 1) then
                tempdp=tempdp+basis_val*u_c_br2_aux_var(k+1,i_var,i_dim,iconsidered)
            else
                tempdp=tempdp+basis_val*u_c_valdg(1,i_var,k+1,iconsidered)
            end if
        end do
        if (br2_yn == 1) then
            dg_sol_der_out(i_var, i_dim) = u_c_br2_aux_var(1,i_var,i_dim,iconsidered)+tempdp
        else
            dg_sol_der_out(i_var, i_dim) = tempdp
        end if
    end do
end do

icompwrt=0

end subroutine dg_sol_der_fill


subroutine br2_local_lift_fill(n, n_qp,facex,iconsidered,wequa2d,br2_local_lift_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> requires: iconsidered, facex, and weights_dg or wequa2d
integer,intent(in)::n, n_qp,facex,iconsidered
real,dimension(1:numberofpoints2),intent(in)::wequa2d
real,dimension(1:nof_variables,1:dimensiona),intent(out)::br2_local_lift_out
integer::i, l, ngp, i_dim,number_of_dog
real::angle1,angle2
real,dimension(gpu_max_dim)::nnn
real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv,srf_speedrot
real::mp_pinfr,gammar
integer::b_code,pointx,number,i_var
real::nx,ny,nz
real,dimension(1:gpu_max_nvar_total)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_dim)::pox,poy,poz

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

br2_local_lift_out = 0.0d0
do ngp = 1, n_qp
    pointx = ngp
    call get_left_right_states(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
    leftv = cleft(1:nof_variables)
    rightv = cright(1:nof_variables)
    call cons2div(n,leftv,mp_pinfl,gammal)
    call cons2div(n,rightv,mp_pinfr,gammar)

    do i_dim = 1, dimensiona
        do i_var = 1, nof_variables
            br2_local_lift_out(i_var,i_dim) = br2_local_lift_out(i_var,i_dim) + &
                nnn(i_dim) * (rightv(i_var) - leftv(i_var)) * wequa2d(pointx)
        end do
    end do
end do

br2_local_lift_out = oo2 * br2_local_lift_out * ielem_surf(facex,iconsidered) / ielem_totvolume(iconsidered)

end subroutine br2_local_lift_fill


subroutine dg_surf_flux_accumulate(n,iconsidered,facex,pointx,weights_temp,rhllcflux,dg_surf_flux_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> accumulates the rhs surface flux term in the dg formulation without returning an array temporary
integer,intent(in)::n,iconsidered,facex,pointx
real,dimension(1:nof_variables),intent(in)::rhllcflux
real,dimension(1:numberofpoints2),intent(in)::weights_temp
real,dimension(1:idegfree+1,1:nof_variables),intent(inout)::dg_surf_flux_out
integer::i,icompwrt,kbasis
real::x1,y1,z1,wface,basis_val
integer::number,number_of_dog

number_of_dog=ielem_idegfree(iconsidered)
number=ielem_iorder(iconsidered)
icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
    z1= rec_qpoints(facex,pointx,3,iconsidered)
end if
wface = weights_temp(pointx)*ielem_surf(facex,iconsidered)

do i = 1, nof_variables
    dg_surf_flux_out(1,i) = dg_surf_flux_out(1,i) + rhllcflux(i) * wface
    do kbasis=1,number_of_dog
        if (dimensiona.eq.2)then
            basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,kbasis)
        else
            basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,kbasis)
        end if
        dg_surf_flux_out(kbasis+1,i) = dg_surf_flux_out(kbasis+1,i) + rhllcflux(i) * wface * basis_val
    end do
end do

icompwrt=0

end subroutine dg_surf_flux_accumulate


subroutine dg_surf_flux_accumulate_rhs(n,iconsidered,facex,pointx,weights_temp,rhllcflux,scale)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> accumulates the rhs surface flux term directly into rhs_valdg for one cell
integer,intent(in)::n,iconsidered,facex,pointx
real,dimension(1:nof_variables),intent(in)::rhllcflux
real,dimension(1:numberofpoints2),intent(in)::weights_temp
real,intent(in)::scale
integer::i,icompwrt,kbasis
real::x1,y1,z1,wface,basis_val,flux_weight
integer::number,number_of_dog

number_of_dog=ielem_idegfree(iconsidered)
number=ielem_iorder(iconsidered)
icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
    z1= rec_qpoints(facex,pointx,3,iconsidered)
end if
wface = scale*weights_temp(pointx)*ielem_surf(facex,iconsidered)

do i = 1, nof_variables
    flux_weight = rhllcflux(i) * wface
    rhs_valdg(1,i,iconsidered) = rhs_valdg(1,i,iconsidered) + flux_weight
    do kbasis=1,number_of_dog
        if (dimensiona.eq.2)then
            basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,kbasis)
        else
            basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,kbasis)
        end if
        rhs_valdg(kbasis+1,i,iconsidered) = rhs_valdg(kbasis+1,i,iconsidered) + flux_weight * basis_val
    end do
end do

icompwrt=0

end subroutine dg_surf_flux_accumulate_rhs


subroutine dg_surf_fluxv_fill(n,iconsidered,facex,pointx,weights_temp,hllcflux,dg_surf_fluxv_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> fills the viscous rhs surface flux term without returning an array temporary
integer,intent(in)::n,iconsidered,facex,pointx
real,dimension(1:nof_variables),intent(in)::hllcflux
real,dimension(1:numberofpoints2),intent(in)::weights_temp
real,dimension(1:idegfree+1,1:nof_variables),intent(out)::dg_surf_fluxv_out
integer::i,icompwrt,kbasis
real::x1,y1,z1,wface,basis_val
integer::number,number_of_dog

number_of_dog=ielem_idegfree(iconsidered)
number=ielem_iorder(iconsidered)
icompwrt=-2
x1= rec_qpoints(facex,pointx,1,iconsidered)
y1= rec_qpoints(facex,pointx,2,iconsidered)
if (dimensiona.eq.3)then
    z1= rec_qpoints(facex,pointx,3,iconsidered)
end if
wface = weights_temp(pointx)*ielem_surf(facex,iconsidered)
dg_surf_fluxv_out = 0.0d0

do i = 1, nof_variables
    dg_surf_fluxv_out(1,i) = hllcflux(i) * wface
    do kbasis=1,number_of_dog
        if (dimensiona.eq.2)then
            basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,kbasis)
        else
            basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,kbasis)
        end if
        dg_surf_fluxv_out(kbasis+1,i) = hllcflux(i) * wface * basis_val
    end do
end do

icompwrt=0

end subroutine dg_surf_fluxv_fill


subroutine dg_vol_integral_fill(n,iconsidered,dg_vol_integral_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> fills the volume integral term in the dg rhs without returning an array temporary
integer,intent(in)::n,iconsidered
real,dimension(1:idegfree+1,1:nof_variables),intent(out)::dg_vol_integral_out
integer::i,nqp,i_qp,number,number_of_dog,i_var,i_dim,i_face
real,dimension(1:gpu_max_nvar)::leftv
real,dimension(1:gpu_max_nvar)::flux_term_x,flux_term_y,flux_term_z
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::leftv_der
real::x1,y1,z1,lamxl,lamyl,lamzl
real::mp_pinfl,gammal

nqp = ielem_itotalpoints(iconsidered)
number=ielem_iorder(iconsidered)
number_of_dog = ielem_idegfree(iconsidered)
dg_vol_integral_out = 0.0d0

if (initcond.eq.3)then
    lamxl=-ielem_yyc(iconsidered)+0.5d0
    lamyl=ielem_xxc(iconsidered)-0.5d0
else
    lamxl=lamx
    lamyl=lamy
end if

if (dimensiona.eq.3) lamzl=lamz

do i=1,number_of_dog
    do i_qp = 1, nqp
        x1=qp_array_x(i_qp,iconsidered)
        y1=qp_array_y(i_qp,iconsidered)
        if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,iconsidered)
        end if

        if (itestcase.lt.3) then
            call dg_sol_fill(n,iconsidered,x1,y1,z1,leftv)
            do i_var=1,nof_variables
                flux_term_x(i_var)=leftv(i_var)*lamxl
                flux_term_y(i_var)=leftv(i_var)*lamyl
                if (dimensiona.eq.3) flux_term_z(i_var)=leftv(i_var)*lamzl
            end do
        end if

        if (itestcase.eq.3) then
            call dg_sol_fill(n,iconsidered,x1,y1,z1,leftv)
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

        if (itestcase == 4) then
            call dg_sol_fill(n,iconsidered,x1,y1,z1,leftv)
            call dg_sol_der_fill(x1,y1,z1,number_of_dog,iconsidered,leftv_der)
            if( br2_yn == 2) then
                call dcons2dprim(leftv_der,leftv)
                do i_var=1,nof_variables
                    do i_dim=1,dimensiona
                        do i_face=1,max_faces
                            leftv_der(i_var,i_dim)=leftv_der(i_var,i_dim)+rec_br2_local_lift(i_var,i_dim,i_face,iconsidered)
                        end do
                    end do
                end do
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
                do i_var=1,nof_variables
                    dg_vol_integral_out(i+1,i_var) = dg_vol_integral_out(i+1,i_var)+ qp_array_qp_weight(i_qp,iconsidered) * &
                        (flux_term_x(i_var)*df2dx(x1,y1,i,iconsidered)+flux_term_y(i_var)*df2dy(x1,y1,i,iconsidered))
                end do
            end if
            if (poly.eq.4)then
                do i_var=1,nof_variables
                    dg_vol_integral_out(i+1,i_var) = dg_vol_integral_out(i+1,i_var)+ qp_array_qp_weight(i_qp,iconsidered)* &
                        (flux_term_x(i_var)*tl2dx(x1,y1,i,iconsidered)+flux_term_y(i_var)*tl2dy(x1,y1,i,iconsidered))
                end do
            end if
        else
            if (poly.eq.1)then
                do i_var=1,nof_variables
                    dg_vol_integral_out(i+1,i_var) = dg_vol_integral_out(i+1,i_var)+ qp_array_qp_weight(i_qp,iconsidered) * &
                        ((flux_term_x(i_var)*dfx(x1,y1,z1,i,iconsidered)) +(flux_term_y(i_var)*dfy(x1,y1,z1,i,iconsidered))+(flux_term_z(i_var)*dfz(x1,y1,z1,i,iconsidered)))
                end do
            end if
            if (poly.eq.2)then
                do i_var=1,nof_variables
                    dg_vol_integral_out(i+1,i_var) = dg_vol_integral_out(i+1,i_var)+ qp_array_qp_weight(i_qp,iconsidered) * &
                        ((flux_term_x(i_var)*dlx(x1,y1,z1,i,iconsidered)) +(flux_term_y(i_var)*dly(x1,y1,z1,i,iconsidered))+(flux_term_z(i_var)*dlz(x1,y1,z1,i,iconsidered)))
                end do
            end if
            if (poly.eq.4)then
                do i_var=1,nof_variables
                    dg_vol_integral_out(i+1,i_var) = dg_vol_integral_out(i+1,i_var)+ qp_array_qp_weight(i_qp,iconsidered)* &
                        ((flux_term_x(i_var)*tl3dx(x1,y1,z1,i,iconsidered)) +(flux_term_y(i_var)*tl3dy(x1,y1,z1,i,iconsidered))+(flux_term_z(i_var)*tl3dz(x1,y1,z1,i,iconsidered)))
                end do
            end if
        end if
    end do
end do

end subroutine dg_vol_integral_fill


subroutine dg_vol_integral_accumulate_rhs(n,iconsidered,scale)
implicit none
#ifdef gpu
!$omp declare target
#endif
!> accumulates the volume integral term directly into rhs_valdg for one cell
integer,intent(in)::n,iconsidered
real,intent(in)::scale
integer::i,nqp,i_qp,number,number_of_dog,i_var,i_dim,i_face
real,dimension(1:gpu_max_nvar)::leftv
real,dimension(1:gpu_max_nvar)::flux_term_x,flux_term_y,flux_term_z
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::leftv_der
real::x1,y1,z1,lamxl,lamyl,lamzl,weighted_flux
real::mp_pinfl,gammal

nqp = ielem_itotalpoints(iconsidered)
number=ielem_iorder(iconsidered)
number_of_dog = ielem_idegfree(iconsidered)

if (initcond.eq.3)then
    lamxl=-ielem_yyc(iconsidered)+0.5d0
    lamyl=ielem_xxc(iconsidered)-0.5d0
else
    lamxl=lamx
    lamyl=lamy
end if

if (dimensiona.eq.3) lamzl=lamz

do i=1,number_of_dog
    do i_qp = 1, nqp
        x1=qp_array_x(i_qp,iconsidered)
        y1=qp_array_y(i_qp,iconsidered)
        if (dimensiona.eq.3)then
            z1=qp_array_z(i_qp,iconsidered)
        end if

        if (itestcase.lt.3) then
            call dg_sol_fill(n,iconsidered,x1,y1,z1,leftv)
            do i_var=1,nof_variables
                flux_term_x(i_var)=leftv(i_var)*lamxl
                flux_term_y(i_var)=leftv(i_var)*lamyl
                if (dimensiona.eq.3) flux_term_z(i_var)=leftv(i_var)*lamzl
            end do
        end if

        if (itestcase.eq.3) then
            call dg_sol_fill(n,iconsidered,x1,y1,z1,leftv)
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

        if (itestcase == 4) then
            call dg_sol_fill(n,iconsidered,x1,y1,z1,leftv)
            call dg_sol_der_fill(x1,y1,z1,number_of_dog,iconsidered,leftv_der)
            if( br2_yn == 2) then
                call dcons2dprim(leftv_der,leftv)
                do i_var=1,nof_variables
                    do i_dim=1,dimensiona
                        do i_face=1,max_faces
                            leftv_der(i_var,i_dim)=leftv_der(i_var,i_dim)+rec_br2_local_lift(i_var,i_dim,i_face,iconsidered)
                        end do
                    end do
                end do
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
                do i_var=1,nof_variables
                    weighted_flux = scale*qp_array_qp_weight(i_qp,iconsidered) * &
                        (flux_term_x(i_var)*df2dx(x1,y1,i,iconsidered)+flux_term_y(i_var)*df2dy(x1,y1,i,iconsidered))
                    rhs_valdg(i+1,i_var,iconsidered) = rhs_valdg(i+1,i_var,iconsidered) + weighted_flux
                end do
            end if
            if (poly.eq.4)then
                do i_var=1,nof_variables
                    weighted_flux = scale*qp_array_qp_weight(i_qp,iconsidered)* &
                        (flux_term_x(i_var)*tl2dx(x1,y1,i,iconsidered)+flux_term_y(i_var)*tl2dy(x1,y1,i,iconsidered))
                    rhs_valdg(i+1,i_var,iconsidered) = rhs_valdg(i+1,i_var,iconsidered) + weighted_flux
                end do
            end if
        else
            if (poly.eq.1)then
                do i_var=1,nof_variables
                    weighted_flux = scale*qp_array_qp_weight(i_qp,iconsidered) * &
                        ((flux_term_x(i_var)*dfx(x1,y1,z1,i,iconsidered)) + &
                        (flux_term_y(i_var)*dfy(x1,y1,z1,i,iconsidered)) + &
                        (flux_term_z(i_var)*dfz(x1,y1,z1,i,iconsidered)))
                    rhs_valdg(i+1,i_var,iconsidered) = rhs_valdg(i+1,i_var,iconsidered) + weighted_flux
                end do
            end if
            if (poly.eq.2)then
                do i_var=1,nof_variables
                    weighted_flux = scale*qp_array_qp_weight(i_qp,iconsidered) * &
                        ((flux_term_x(i_var)*dlx(x1,y1,z1,i,iconsidered)) + &
                        (flux_term_y(i_var)*dly(x1,y1,z1,i,iconsidered)) + &
                        (flux_term_z(i_var)*dlz(x1,y1,z1,i,iconsidered)))
                    rhs_valdg(i+1,i_var,iconsidered) = rhs_valdg(i+1,i_var,iconsidered) + weighted_flux
                end do
            end if
            if (poly.eq.4)then
                do i_var=1,nof_variables
                    weighted_flux = scale*qp_array_qp_weight(i_qp,iconsidered)* &
                        ((flux_term_x(i_var)*tl3dx(x1,y1,z1,i,iconsidered)) + &
                        (flux_term_y(i_var)*tl3dy(x1,y1,z1,i,iconsidered)) + &
                        (flux_term_z(i_var)*tl3dz(x1,y1,z1,i,iconsidered)))
                    rhs_valdg(i+1,i_var,iconsidered) = rhs_valdg(i+1,i_var,iconsidered) + weighted_flux
                end do
            end if
        end if
    end do
end do

end subroutine dg_vol_integral_accumulate_rhs


subroutine dg_vol_integral2_fill(n,i,dg_vol_integral2_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i
real,dimension(1:nof_variables),intent(out)::dg_vol_integral2_out
integer::k,nqp,i_qp,i_var,icompwrt,iconsidered
real::x1,y1,z1,tempdp,basis_val
integer::number,number_of_dog

iconsidered=i
icompwrt=-2
nqp = ielem_itotalpoints(iconsidered)
number=ielem_iorder(iconsidered)
number_of_dog = ielem_idegfree(iconsidered)
dg_vol_integral2_out = 0.0d0

do i_qp = 1, nqp
    x1=qp_array_x(i_qp,iconsidered)
    y1=qp_array_y(i_qp,iconsidered)
    if (dimensiona.eq.3)then
        z1=qp_array_z(i_qp,iconsidered)
    end if
    do i_var = 1, nof_variables
        tempdp=0.0d0
        do k=1,number_of_dog
            if (dimensiona.eq.2)then
                basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
            else
                basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
            end if
            tempdp=tempdp+basis_val*u_c_valdg(1,i_var,k+1,iconsidered)
        end do
        dg_vol_integral2_out(i_var) = dg_vol_integral2_out(i_var)+ (qp_array_qp_weight(i_qp,iconsidered)*tempdp/ielem_totvolume(iconsidered))
    end do
end do

do i_var = 1, nof_variables
    dg_vol_integral2_out(i_var) = u_c_valdg(1,i_var,1,iconsidered)+dg_vol_integral2_out(i_var)
end do
icompwrt=0

end subroutine dg_vol_integral2_fill


subroutine dg_vol_integral_strong_fill(n,i,dg_vol_integral_strong_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i
real,dimension(1:nof_variables),intent(out)::dg_vol_integral_strong_out
integer::j,k,nqp,i_qp,i_var,icompwrt,iconsidered
real::x1,y1,z1,tempdp,basis_val
integer::number,number_of_dog

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
dg_vol_integral_strong_out = 0.0d0

do i_qp = 1, nqp
    x1=qp_array_x(i_qp,iconsidered)
    y1=qp_array_y(i_qp,iconsidered)
    if (dimensiona.eq.3)then
        z1=qp_array_z(i_qp,iconsidered)
    end if
    do i_var = 1, nof_variables
        tempdp=0.0d0
        do k=1,number_of_dog
            if (dimensiona.eq.2)then
                basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
            else
                basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
            end if
            tempdp=tempdp+basis_val*u_cs_valdg(1,i_var,k+1,iconsidered)
        end do
        dg_vol_integral_strong_out(i_var)= dg_vol_integral_strong_out(i_var)+ (qp_array_qp_weight(i_qp,iconsidered)*tempdp/ielem_totvolume(iconsidered))
    end do
end do

do i_var = 1, nof_variables
    dg_vol_integral_strong_out(i_var) = u_c_valdg(1,i_var,1,iconsidered)+dg_vol_integral_strong_out(i_var)
end do
icompwrt=0

end subroutine dg_vol_integral_strong_fill


subroutine dg_vol_integral_weak_fill(n,i,dg_vol_integral_weak_out)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i
real,dimension(1:nof_variables),intent(out)::dg_vol_integral_weak_out
integer::j,k,nqp,i_qp,i_var,icompwrt,iconsidered
real::x1,y1,z1,tempdp,basis_val
integer::number,number_of_dog

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
dg_vol_integral_weak_out = 0.0d0

do i_qp = 1, nqp
    x1=qp_array_x(i_qp,iconsidered)
    y1=qp_array_y(i_qp,iconsidered)
    if (dimensiona.eq.3)then
        z1=qp_array_z(i_qp,iconsidered)
    end if
    do i_var = 1, nof_variables
        tempdp=0.0d0
        do k=1,number_of_dog
            if (dimensiona.eq.2)then
                basis_val = basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,k)
            else
                basis_val = basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,k)
            end if
            tempdp=tempdp+basis_val*u_cw_valdg(1,i_var,k+1,iconsidered)
        end do
        dg_vol_integral_weak_out(i_var)= dg_vol_integral_weak_out(i_var)+ (qp_array_qp_weight(i_qp,iconsidered)*tempdp/ielem_totvolume(iconsidered))
    end do
end do

do i_var = 1, nof_variables
    dg_vol_integral_weak_out(i_var) = u_c_valdg(1,i_var,1,iconsidered)+dg_vol_integral_weak_out(i_var)
end do
icompwrt=0

end subroutine dg_vol_integral_weak_fill


subroutine reconstruct_dg(n)
implicit none
integer,intent(in)::n
integer::i_face, i_elem, i_qp,iqp
integer::facex,pointx,iconsidered,number_of_dog


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(i_face, i_elem, i_qp, iqp, facex, pointx, iconsidered, number_of_dog)
#else
!$omp do
#endif
do i_elem = 1, xmpielrank(n)

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

            call dg_solface_store(n,facex,pointx,iconsidered,number_of_dog)



        end do
    end do

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine reconstruct_dg


subroutine reconstruct_br2_dg(n)
implicit none
integer,intent(in)::n
integer::i_elem




! needed for dg_surf_flux, used in br2_local_lift


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i_elem = 1, xmpielrank(n)
    call reconstruct_br2_dg_cell(n,i_elem)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine reconstruct_br2_dg

subroutine reconstruct_br2_dg_cell(n,i_elem)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,i_elem
integer::i_face, i_qp, iqp, iconsidered, facex, pointx, number_of_dog
integer::i_dim, i_var
real,dimension(1:gpu_max_qp_face)::weights_temp
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::leftv_der
real,dimension(1:gpu_max_nvar)::leftv
real::x1,y1,z1

iconsidered = i_elem
do i_face = 1, ielem_ifca(i_elem)
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

    facex = i_face
    call br2_local_lift_fill(n, iqp,facex,iconsidered,weights_temp,rec_br2_local_lift(:,:,i_face,iconsidered))

    do i_qp = 1, iqp
        pointx=i_qp
        if (br2_yn /= 1 .and. itestcase == 4) then
            x1 = rec_qpoints(facex, pointx,1,iconsidered)
            y1 = rec_qpoints(facex, pointx,2,iconsidered)
            if (dimensiona == 3) z1 = rec_qpoints(facex, pointx,3,iconsidered)

            call dg_sol_der_fill(x1,y1,z1,number_of_dog,iconsidered,leftv_der)
            call dcons2dprim(leftv_der,leftv)
            do i_var=1,nof_variables
                do i_dim=1,dimensiona
                    rec_br2_aux_var(i_var,i_dim,i_face,i_qp,iconsidered)=leftv_der(i_var,i_dim) &
                    +rec_br2_local_lift(i_var,i_dim,i_face,iconsidered)*br2_damping
                end do
            end do

            do i_dim=1,dims
                do i_var=1,nof_variables-1
                    rec_uleftv(i_dim,i_var,i_face,i_qp,iconsidered)=rec_br2_aux_var(i_var,i_dim,i_face,i_qp,iconsidered)
                end do
            end do
        end if
    end do
end do

end subroutine reconstruct_br2_dg_cell

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
real,dimension(1:gpu_max_nvar)::srf_speed
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


subroutine get_states_interior_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex,pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar)::srf_speed
integer::i,l,ngp,iv,nvt

i=iconsidered
l=facex
ngp=pointx
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
if (nvt.gt.0)then
  if (icoupleturb.eq.1)then
    do iv=1,nvt
      cturbl(iv)=rec_uleftturb(iv,l,ngp,i)
      cturbr(iv)=rec_uleftturb(iv,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    end do
  else
    do iv=1,nvt
      cturbl(iv)=u_ct_val(1,iv,i)
      cturbr(iv)=u_ct_val(1,iv,ielem_ineigh(l,i))
    end do
  end if
  do iv=1,nvt
    cleft_rot(nof_variables+iv)=cturbl(iv)
    cright_rot(nof_variables+iv)=cturbr(iv)
  end do
end if
if (rec_mrf(i).eq.1)then
  srf_speed=zero
  srf_speed(2:4)=rec_rotvel(l,ngp,1:3,i)
  call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
end if

end subroutine get_states_interior_ideal


subroutine get_states_interior2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex,pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
integer::i,l,ngp,iv,nvt

i=iconsidered
l=facex
ngp=pointx
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
if (nvt.gt.0)then
  if (icoupleturb.eq.1)then
    do iv=1,nvt
      cturbl(iv)=rec_uleftturb(iv,l,ngp,i)
      cturbr(iv)=rec_uleftturb(iv,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    end do
  else
    do iv=1,nvt
      cturbl(iv)=u_ct_val(1,iv,i)
      cturbr(iv)=u_ct_val(1,iv,ielem_ineigh(l,i))
    end do
  end if
  do iv=1,nvt
    cleft_rot(nof_variables+iv)=cturbl(iv)
    cright_rot(nof_variables+iv)=cturbr(iv)
  end do
end if

end subroutine get_states_interior2d_ideal


subroutine calculate_interior_viscous_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny, &
     & nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex,pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
integer::i,l,ngp,k,nvar

i=iconsidered
l=facex
ngp=pointx
call get_states_interior_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
do k=1,dimensiona
  lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i)
  rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
end do
if ((turbulence.eq.1).or.(passivescalar.gt.0))then
  do nvar=1,turbulenceequations+passivescalar
    lcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,l,ngp,i)
    rcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
  end do
end if

end subroutine calculate_interior_viscous_ideal


subroutine calculate_interior_viscous2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex,pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
integer::i,l,ngp,k,nvar

i=iconsidered
l=facex
ngp=pointx
call get_states_interior2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
do k=1,dimensiona
  lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i)
  rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
end do
if ((turbulence.eq.1).or.(passivescalar.gt.0))then
  do nvar=1,turbulenceequations+passivescalar
    lcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,l,ngp,i)
    rcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
  end do
end if

end subroutine calculate_interior_viscous2d_ideal



subroutine calculate_bounded_viscous_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny, &
     & nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,facex,pointx
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables,1:dimensiona),intent(inout)::lcvgrad,rcvgrad
real,dimension(1:turbulenceequations+passivescalar,1:dimensiona),intent(inout)::lcvgrad_t,rcvgrad_t
integer::i,l,ngp,k,nvar,iex,ittt,nfx,lfx,rowfx,nvt,d,ii,j
real,dimension(1:3)::nnt,g
real,dimension(1:3,1:3)::a,b,c,r
real::g_n

i=iconsidered
l=facex
ngp=pointx
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
lcvgrad=zero
rcvgrad=zero
if (nvt.gt.0)then
  lcvgrad_t=zero
  rcvgrad_t=zero
end if

call get_states_bounds_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)

leftv(1:nof_variables)=cleft(1:nof_variables)
rightv(1:nof_variables)=cright(1:nof_variables)
do k=1,dimensiona
  lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i)
end do
if (nvt.gt.0)then
  do nvar=1,nvt
    lcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,l,ngp,i)
  end do
end if

if (ielem_ineighb(l,i).eq.n)then
  if (ielem_ibounds(l,i).gt.0)then
    if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
      do k=1,dimensiona
        rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      end do
      if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
        do k=1,dimensiona
          rcvgrad(1:3,k)=rotate_per_1(rcvgrad(1:3,k),ibound_icode(ielem_ibounds(l,i)),angle_per)
        end do
      end if
      if (nvt.gt.0)then
        do nvar=1,nvt
          rcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
        end do
      end if
    else
      rcvgrad(1:nof_variables-1,1:dimensiona)=lcvgrad(1:nof_variables-1,1:dimensiona)
      if (nvt.gt.0) rcvgrad_t(1:nvt,1:dimensiona)=lcvgrad_t(1:nvt,1:dimensiona)
      if (((b_code.eq.4).or.(b_code.eq.2)).and.(thermal.ne.1))then
        nnt(1)=nx
        nnt(2)=ny
        nnt(3)=nz
        do iex=dimensiona+1,nof_variables-1
          g_n=zero
          do d=1,dimensiona
            g(d)=lcvgrad(iex,d)
            g_n=g_n+g(d)*nnt(d)
          end do
          do d=1,dimensiona
            g(d)=g(d)-g_n*nnt(d)
            lcvgrad(iex,d)=g(d)
            rcvgrad(iex,d)=g(d)
          end do
        end do
      else if (b_code.eq.3)then
        nnt(1)=nx
        nnt(2)=ny
        nnt(3)=nz
        do iex=dimensiona+1,nof_variables-1
          g_n=zero
          do d=1,dimensiona
            g(d)=lcvgrad(iex,d)
            g_n=g_n+g(d)*nnt(d)
          end do
          do d=1,dimensiona
            rcvgrad(iex,d)=lcvgrad(iex,d)-2.0d0*g_n*nnt(d)
          end do
        end do
        if (nvt.gt.0)then
          do iex=1,nvt
            g_n=zero
            do d=1,dimensiona
              g(d)=lcvgrad_t(iex,d)
              g_n=g_n+g(d)*nnt(d)
            end do
            do d=1,dimensiona
              rcvgrad_t(iex,d)=lcvgrad_t(iex,d)-2.0d0*g_n*nnt(d)
            end do
          end do
        end if
        if (dimensiona.eq.3)then
          do ii=1,3
            do j=1,3
              r(ii,j)=zero
            end do
            r(ii,ii)=1.0d0
          end do
          do ii=1,3
            do j=1,3
              r(ii,j)=r(ii,j)-2.0d0*nnt(ii)*nnt(j)
            end do
          end do
          a(1,1)=lcvgrad(1,1); a(1,2)=lcvgrad(1,2); a(1,3)=lcvgrad(1,3)
          a(2,1)=lcvgrad(2,1); a(2,2)=lcvgrad(2,2); a(2,3)=lcvgrad(2,3)
          a(3,1)=lcvgrad(3,1); a(3,2)=lcvgrad(3,2); a(3,3)=lcvgrad(3,3)
          do ii=1,3
            do j=1,3
              c(ii,j)=zero
              do k=1,3
                c(ii,j)=c(ii,j)+r(ii,k)*a(k,j)
              end do
            end do
          end do
          do ii=1,3
            do j=1,3
              b(ii,j)=zero
              do k=1,3
                b(ii,j)=b(ii,j)+c(ii,k)*r(k,j)
              end do
            end do
          end do
          do d=1,3
            rcvgrad(1,d)=b(1,d)
            rcvgrad(2,d)=b(2,d)
            rcvgrad(3,d)=b(3,d)
          end do
        end if
      end if
    end if
  else
    do k=1,dimensiona
      rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    end do
    if (nvt.gt.0)then
      do nvar=1,nvt
        rcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      end do
    end if
  end if
else
  nfx=ielem_ineighn(l,i)
  lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
  rowfx=bound_offset(nfx)+lfx-1
  ittt=0
  do iex=1,nof_variables-1
    do nvar=1,dimensiona
      ittt=ittt+1
      rcvgrad(iex,nvar)=boundhir(rowfx,nof_variables+nvt+ittt)
    end do
  end do
  if ((ielem_ibounds(l,i).gt.0).and.(per_rot.eq.1))then
    if (ibound_icode(ielem_ibounds(l,i)).eq.50)then
      do k=1,dimensiona
        rcvgrad(1:3,k)=rotate_per_1(rcvgrad(1:3,k),ibound_icode(ielem_ibounds(l,i)),angle_per)
      end do
    end if
  end if
  if (nvt.gt.0)then
    do iex=1,nvt
      do nvar=1,dimensiona
        ittt=ittt+1
        rcvgrad_t(iex,nvar)=boundhir(rowfx,nof_variables+nvt+ittt)
      end do
    end do
  end if
end if

end subroutine calculate_bounded_viscous_ideal


subroutine calculate_bounded_viscous2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny, &
     & nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n,iconsidered,facex,pointx
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:nof_variables,1:dimensiona),intent(inout)::lcvgrad,rcvgrad
real,dimension(1:turbulenceequations+passivescalar,1:dimensiona),intent(inout)::lcvgrad_t,rcvgrad_t
integer::i,l,ngp,k,nvar,iex,ittt,nfx,lfx,rowfx,nvt,d,ii,j
real,dimension(1:2)::nnt,g
real,dimension(1:2,1:2)::a,b,c,r
real::g_n

i=iconsidered
l=facex
ngp=pointx
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
lcvgrad=zero
rcvgrad=zero
if (nvt.gt.0)then
  lcvgrad_t=zero
  rcvgrad_t=zero
end if

call get_states_bounds2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)

leftv(1:nof_variables)=cleft(1:nof_variables)
rightv(1:nof_variables)=cright(1:nof_variables)
do k=1,dimensiona
  lcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,l,ngp,i)
end do
if (nvt.gt.0)then
  do nvar=1,nvt
    lcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,l,ngp,i)
  end do
end if

if (ielem_ineighb(l,i).eq.n)then
  if (ielem_ibounds(l,i).gt.0)then
    if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
      do k=1,dimensiona
        rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      end do
      if (nvt.gt.0)then
        do nvar=1,nvt
          rcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
        end do
      end if
    else
      rcvgrad(1:nof_variables-1,1:dimensiona)=lcvgrad(1:nof_variables-1,1:dimensiona)
      if (nvt.gt.0) rcvgrad_t(1:nvt,1:dimensiona)=lcvgrad_t(1:nvt,1:dimensiona)
      if (((b_code.eq.4).or.(b_code.eq.2)).and.(thermal.ne.1))then
        nnt(1)=nx
        nnt(2)=ny
        do iex=dimensiona+1,nof_variables-1
          g_n=zero
          do d=1,dimensiona
            g(d)=lcvgrad(iex,d)
            g_n=g_n+g(d)*nnt(d)
          end do
          do d=1,dimensiona
            g(d)=g(d)-g_n*nnt(d)
            lcvgrad(iex,d)=g(d)
            rcvgrad(iex,d)=g(d)
          end do
        end do
      else if (b_code.eq.3)then
        nnt(1)=nx
        nnt(2)=ny
        do iex=dimensiona+1,nof_variables-1
          g_n=zero
          do d=1,dimensiona
            g(d)=lcvgrad(iex,d)
            g_n=g_n+g(d)*nnt(d)
          end do
          do d=1,dimensiona
            rcvgrad(iex,d)=lcvgrad(iex,d)-2.0d0*g_n*nnt(d)
          end do
        end do
        if (nvt.gt.0)then
          do iex=1,nvt
            g_n=zero
            do d=1,dimensiona
              g(d)=lcvgrad_t(iex,d)
              g_n=g_n+g(d)*nnt(d)
            end do
            do d=1,dimensiona
              rcvgrad_t(iex,d)=lcvgrad_t(iex,d)-2.0d0*g_n*nnt(d)
            end do
          end do
        end if
        do ii=1,2
          do j=1,2
            r(ii,j)=zero
          end do
          r(ii,ii)=1.0d0
        end do
        do ii=1,2
          do j=1,2
            r(ii,j)=r(ii,j)-2.0d0*nnt(ii)*nnt(j)
          end do
        end do
        a(1,1)=lcvgrad(1,1); a(1,2)=lcvgrad(1,2)
        a(2,1)=lcvgrad(2,1); a(2,2)=lcvgrad(2,2)
        do ii=1,2
          do j=1,2
            c(ii,j)=zero
            do k=1,2
              c(ii,j)=c(ii,j)+r(ii,k)*a(k,j)
            end do
          end do
        end do
        do ii=1,2
          do j=1,2
            b(ii,j)=zero
            do k=1,2
              b(ii,j)=b(ii,j)+c(ii,k)*r(k,j)
            end do
          end do
        end do
        do d=1,2
          rcvgrad(1,d)=b(1,d)
          rcvgrad(2,d)=b(2,d)
        end do
      end if
    end if
  else
    do k=1,dimensiona
      rcvgrad(1:nof_variables-1,k)=rec_uleftv(k,1:nof_variables-1,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    end do
    if (nvt.gt.0)then
      do nvar=1,nvt
        rcvgrad_t(nvar,1:dimensiona)=rec_uleftturbv(1:dimensiona,nvar,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      end do
    end if
  end if
else
  nfx=ielem_ineighn(l,i)
  lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
  rowfx=bound_offset(nfx)+lfx-1
  ittt=0
  do iex=1,nof_variables-1
    do nvar=1,dimensiona
      ittt=ittt+1
      rcvgrad(iex,nvar)=boundhir(rowfx,nof_variables+nvt+ittt)
    end do
  end do
  if (nvt.gt.0)then
    do iex=1,nvt
      do nvar=1,dimensiona
        ittt=ittt+1
        rcvgrad_t(iex,nvar)=boundhir(rowfx,nof_variables+nvt+ittt)
      end do
    end do
  end if
end if

end subroutine calculate_bounded_viscous2d_ideal


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
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_dim)::cords
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
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_dim)::cords
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
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_dim)::cords
real,dimension(1:3):: gl,nnt,g
real,dimension(1:3,1:3)::a,b,c,r
real::g_n
integer :: d
integer::i,l,ngp,ittt,nvar,iex,kk,n_node,k,nfx,lfx,rowfx,j,ii
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
										if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
	                                        cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)

	                                        do k=1,dimensiona
	                                    rcvgrad(1:3,k)=rotate_per_1(rcvgrad(1:3,k),ibound_icode(ielem_ibounds(l,i)),angle_per)
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


								    call cordinates3(n,nodes_list,n_node,cords(1:3))

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

                                                  if (catalytic_wall.eq.0)then
                                                    nnt(1)=nx
                                                    nnt(2)=ny
                                                    nnt(3)=nz

                                                    do iex = dimensiona+3,nof_variables-1
                                                      g_n=0.0d0
                                                      do d=1,dimensiona
                                                        g(d)=lcvgrad(iex,d)
                                                        g_n=g_n+g(d)*nnt(d)
                                                      end do

                                                      do d=1,dimensiona
                                                        rcvgrad(iex,d)=lcvgrad(iex,d)-2.0d0*g_n*nnt(d)
                                                      end do
                                                    end do
                                                  end if

                                            end if




				  				    else

                                                  do k=1,dimensiona
                                                  rcvgrad(1:nof_variables-1,k)=lcvgrad(1:nof_variables-1,k)
                                                  end do
                                                  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                                  rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                  end if



				  				     if ((b_code.eq.3))then

                                         ! normal vector (global coords)
                                              nnt(1)=nx
                                              nnt(2)=ny
                                              if (dimensiona == 3) then
                                                nnt(3)=nz
                                              else
                                                nnt(3)=0.0
                                              end if

                                              do iex = dimensiona+1, nof_variables-1   ! temperature variable(s)
                                                ! g = grad(T)^- from left
                                                g_n = 0.0
                                                do d=1,dimensiona
                                                  g(d) = lcvgrad(iex,d)
                                                  g_n  = g_n + g(d)*nnt(d)
                                                end do
                                                if (dimensiona == 2) g(3)=0.0

                                                ! symmetry ghost gradient: grad^+ = grad^- - 2 (grad^-·n) n
                                                do d=1,dimensiona
                                                  rcvgrad(iex,d) = lcvgrad(iex,d) - 2.0*g_n*nnt(d)
                                                end do
                                              end do

                                              if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                              do iex = 1,turbulenceequations+passivescalar
                                                  !rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                  g_n = 0.0
                                                do d=1,dimensiona
                                                  g(d) = lcvgrad_t(iex,d)
                                                  g_n  = g_n + g(d)*nnt(d)
                                                end do
                                                if (dimensiona == 2) g(3)=0.0

                                                ! symmetry ghost gradient: grad^+ = grad^- - 2 (grad^-·n) n
                                                do d=1,dimensiona
                                                  rcvgrad_t(iex,d) = lcvgrad_t(iex,d) - 2.0*g_n*nnt(d)
                                                end do


                                                  end do
                                                  end if




                                              if (dimensiona == 3) then
                                                ! Build R = I - 2 n n^T (3x3)
                                                do ii=1,3
                                                  do j=1,3
                                                    R(ii,j)=0.0
                                                  end do
                                                  R(ii,ii)=1.0
                                                end do
                                                do ii=1,3
                                                  do j=1,3
                                                    R(ii,j) = R(ii,j) - 2.0*nnt(ii)*nnt(j)
                                                  end do
                                                end do

                                                ! A = grad(u)^- tensor from left: A(i,j)=d u_i / d x_j
                                                A(1,1)=lcvgrad(1,1); A(1,2)=lcvgrad(1,2); A(1,3)=lcvgrad(1,3)
                                                A(2,1)=lcvgrad(2,1); A(2,2)=lcvgrad(2,2); A(2,3)=lcvgrad(2,3)
                                                A(3,1)=lcvgrad(3,1); A(3,2)=lcvgrad(3,2); A(3,3)=lcvgrad(3,3)

                                                ! C = R*A
                                                do ii=1,3
                                                  do j=1,3
                                                    C(ii,j)=0.0
                                                    do k=1,3
                                                      C(ii,j)=C(ii,j)+R(ii,k)*A(k,j)
                                                    end do
                                                  end do
                                                end do

                                                ! B = C*R
                                                do ii=1,3
                                                  do j=1,3
                                                    B(ii,j)=0.0
                                                    do k=1,3
                                                      B(ii,j)=B(ii,j)+C(ii,k)*R(k,j)
                                                    end do
                                                  end do
                                                end do

                                                ! Store ghost gradients (rcvgrad) for velocity components iex=1..3
                                                do d=1,3
                                                  rcvgrad(1,d)=B(1,d)
                                                  rcvgrad(2,d)=B(2,d)
                                                  rcvgrad(3,d)=B(3,d)
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

										  if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then


										   do k=1,dimensiona
	                                    rcvgrad(1:3,k)=rotate_per_1(rcvgrad(1:3,k),ibound_icode(ielem_ibounds(l,i)),angle_per)
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
real,dimension(1:gpu_max_nvar)::srf_speed,tempx_l,rtempx_l
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar,1:gpu_max_dim)::lcvgrad,rcvgrad
real,dimension(gpu_max_extra_transport,1:gpu_max_dim)::lcvgrad_t,rcvgrad_t
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_dim)::cords
real,dimension(1:3):: gl,nnt,g,grad
real,dimension(1:2,1:2)::a,b,c,r
real::g_n
integer :: d
integer::i,l,ngp,ittt,nvar,iex,n_node,k,nfx,lfx,rowfx,j,ii
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




								  else
								  !not periodic ones in my cpu


								  call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)


                                            n_node=2


								    cords(1:2)=zero
								    call cordinates2(n,nodes_list,n_node,cords(1:2))

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

                                                     do k=1,dimensiona
                                                  rcvgrad(1:nof_variables-1,k)=lcvgrad(1:nof_variables-1,k)
                                                  end do
                                                  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                                  rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                  end if



                                                    if (b_code.eq.2) then


                                                  do k = 1, dimensiona
                                                    rcvgrad(1:dimensiona,k) = lcvgrad(1:dimensiona,k)
                                                  end do

                                                  ! Enforce zero normal gradient for scalar thermochemical variables:
                                                  ! T, Tv, species
                                                  do iex = dimensiona+1, nof_variables-1

                                                    tempx_l  = 0.0d0
                                                    rtempx_l = 0.0d0

                                                    ! gradient vector of scalar iex
                                                    tempx_l(2) = lcvgrad(iex,1)
                                                    tempx_l(3) = lcvgrad(iex,2)

                                                    ! rotate to local face coordinates
                                                    call rotatef2d(n,rtempx_l,tempx_l,angle1,angle2)

                                                    ! reflect normal component only
                                                    rtempx_l(2) = -rtempx_l(2)

                                                    ! keep tangential component unchanged:
                                                    ! rtempx_l(3) = rtempx_l(3)

                                                    ! rotate back to Cartesian
                                                    call rotateb2d(n,tempx_l,rtempx_l,angle1,angle2)

                                                    rcvgrad(iex,1) = tempx_l(2)
                                                    rcvgrad(iex,2) = tempx_l(3)

                                                  end do

                                                end if


                                                      if (b_code.eq.3)then
                                                       ! normal vector (global coords)
                                                      nnt(1)=nx
                                                      nnt(2)=ny
                                                      if (dimensiona == 3) then
                                                        nnt(3)=nz
                                                      else
                                                        nnt(3)=0.0
                                                      end if

                                                      do iex = dimensiona+1, nof_variables-1   ! temperature variable(s)
                                                        ! g = grad(T)^- from left
                                                        g_n = 0.0
                                                        do d=1,dimensiona
                                                          g(d) = lcvgrad(iex,d)
                                                          g_n  = g_n + g(d)*nnt(d)
                                                        end do
                                                        if (dimensiona == 2) g(3)=0.0

                                                        ! symmetry ghost gradient: grad^+ = grad^- - 2 (grad^-·n) n
                                                        do d=1,dimensiona
                                                          rcvgrad(iex,d) = lcvgrad(iex,d) - 2.0*g_n*nnt(d)
                                                        end do
                                                      end do

                                                      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
                                              do iex = 1,turbulenceequations+passivescalar
                                                  !rcvgrad_t(:,:)=lcvgrad_t(:,:)
                                                  g_n = 0.0
                                                do d=1,dimensiona
                                                  g(d) = lcvgrad_t(iex,d)
                                                  g_n  = g_n + g(d)*nnt(d)
                                                end do
                                                if (dimensiona == 2) g(3)=0.0

                                                ! symmetry ghost gradient: grad^+ = grad^- - 2 (grad^-·n) n
                                                do d=1,dimensiona
                                                  rcvgrad_t(iex,d) = lcvgrad_t(iex,d) - 2.0*g_n*nnt(d)
                                                end do


                                                  end do
                                                  end if









                                                          ! --- Build reflection matrix R = I - 2 n n^T ---

                                                          R(1,1) = 1.0 - 2.0*nx*nx
                                                          R(1,2) =      - 2.0*nx*ny
                                                          R(2,1) =      - 2.0*ny*nx
                                                          R(2,2) = 1.0 - 2.0*ny*ny

                                                          ! --- Build left gradient tensor A ---
                                                          ! A(i,j) = d u_i / d x_j

                                                          A(1,1) = lcvgrad(1,1)   ! du/dx
                                                          A(1,2) = lcvgrad(1,2)   ! du/dy
                                                          A(2,1) = lcvgrad(2,1)   ! dv/dx
                                                          A(2,2) = lcvgrad(2,2)   ! dv/dy

                                                          ! --- Compute C = R * A ---
                                                          do ii=1,2
                                                            do j=1,2
                                                              C(ii,j) = 0.0
                                                              do k=1,2
                                                                C(ii,j) = C(ii,j) + R(ii,k)*A(k,j)
                                                              end do
                                                            end do
                                                          end do

                                                          ! --- Compute B = C * R ---
                                                          do ii=1,2
                                                            do j=1,2
                                                              B(ii,j) = 0.0
                                                              do k=1,2
                                                                B(ii,j) = B(ii,j) + C(ii,k)*R(k,j)
                                                              end do
                                                            end do
                                                          end do

                                                          ! --- Store ghost gradients into rcvgrad ---
                                                          rcvgrad(1,1) = B(1,1)
                                                          rcvgrad(1,2) = B(1,2)
                                                          rcvgrad(2,1) = B(2,1)
                                                          rcvgrad(2,2) = B(2,2)



















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


subroutine get_states_bounds_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex,pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:8,1:gpu_max_dim)::vext,nodes_list
real,dimension(1:gpu_max_dim)::cords
integer::i,l,ngp,n_node,ibfc,nfx,lfx,rowfx,iv,nvt

i=iconsidered
l=facex
ngp=pointx
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
srf_speed=zero
cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
if (rec_mrf(i).eq.1)then
  srf_speed(2:4)=rec_rotvel(l,ngp,1:3,i)
  call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
end if
if (nvt.gt.0)then
  if (icoupleturb.eq.1)then
    cturbl(1:nvt)=rec_uleftturb(1:nvt,l,ngp,i)
  else
    cturbl(1:nvt)=u_ct_val(1,1:nvt,i)
  end if
end if

if (ielem_ineighb(l,i).eq.n)then
  if (ielem_ibounds(l,i).gt.0)then
    if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50)) &
        cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
      if (nvt.gt.0)then
        if (icoupleturb.eq.1)then
          cturbr(1:nvt)=rec_uleftturb(1:nvt,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
        else
          cturbr(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(l,i))
        end if
      end if
    else
      call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)
      if (ielem_types_faces(facex,iconsidered).eq.5)then
        n_node=4
      else
        n_node=3
      end if
      call cordinates3(n,nodes_list,n_node,cords(1:3))
      pox(1)=cords(1)
      poy(1)=cords(2)
      poz(1)=cords(3)
      leftv(1:nof_variables)=cleft(1:nof_variables)
      rightv=zero
      ibfc=0
      b_code=ibound_icode(ielem_ibounds(l,i))
      call boundarys_ideal(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
           & cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
      cright(1:nof_variables)=rightv(1:nof_variables)
    end if
  else
    cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    if (nvt.gt.0)then
      if (icoupleturb.eq.1)then
        cturbr(1:nvt)=rec_uleftturb(1:nvt,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      else
        cturbr(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(l,i))
      end if
    end if
  end if
else
  nfx=ielem_ineighn(l,i)
  lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
  rowfx=bound_offset(nfx)+lfx-1
  cright(1:nof_variables)=boundhir(rowfx,1:nof_variables)
  if ((ielem_ibounds(l,i).gt.0).and.(per_rot.eq.1))then
    if (ibound_icode(ielem_ibounds(l,i)).eq.50)then
      cright(2:4)=rotate_per_1(cright(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
    end if
  end if
  if (nvt.gt.0) cturbr(1:nvt)=boundhir(rowfx,nof_variables+1:nof_variables+nvt)
end if
if (nvt.gt.0)then
  do iv=1,nvt
    cleft_rot(nof_variables+iv)=cturbl(iv)
    cright_rot(nof_variables+iv)=cturbr(iv)
  end do
end if

end subroutine get_states_bounds_ideal


subroutine get_states_bounds2d_ideal(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
     & cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex,pointx,n
integer,intent(inout)::b_code
real,intent(inout)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cleft,cright,cright_rot,cleft_rot
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables),intent(inout)::leftv,rightv,srf_speedrot
real,dimension(1:dimensiona),intent(inout)::pox,poy,poz
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:8,1:gpu_max_dim)::vext,nodes_list
real,dimension(1:gpu_max_dim)::cords
integer::i,l,ngp,ibfc,nfx,lfx,rowfx,iv,nvt

i=iconsidered
l=facex
ngp=pointx
nvt=passivescalar
if (turbulence.eq.1) nvt=turbulenceequations+passivescalar
srf_speed=zero
cleft(1:nof_variables)=rec_uleft(1:nof_variables,l,ngp,i)
if (nvt.gt.0)then
  if (icoupleturb.eq.1)then
    cturbl(1:nvt)=rec_uleftturb(1:nvt,l,ngp,i)
  else
    cturbl(1:nvt)=u_ct_val(1,1:nvt,i)
  end if
end if

if (ielem_ineighb(l,i).eq.n)then
  if (ielem_ibounds(l,i).gt.0)then
    if (ibound_icode(ielem_ibounds(l,i)).eq.5)then
      cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      if (nvt.gt.0)then
        if (icoupleturb.eq.1)then
          cturbr(1:nvt)=rec_uleftturb(1:nvt,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
        else
          cturbr(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(l,i))
        end if
      end if
    else
      call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
      call cordinates2(n,nodes_list,2,cords(1:2))
      pox(1)=cords(1)
      poy(1)=cords(2)
      leftv(1:nof_variables)=cleft(1:nof_variables)
      rightv=zero
      ibfc=0
      b_code=ibound_icode(ielem_ibounds(l,i))
      call boundarys2d_ideal(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz, &
           & cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
      cright(1:nof_variables)=rightv(1:nof_variables)
    end if
  else
    cright(1:nof_variables)=rec_uleft(1:nof_variables,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
    if (nvt.gt.0)then
      if (icoupleturb.eq.1)then
        cturbr(1:nvt)=rec_uleftturb(1:nvt,ielem_ineighn(l,i),ngp,ielem_ineigh(l,i))
      else
        cturbr(1:nvt)=u_ct_val(1,1:nvt,ielem_ineigh(l,i))
      end if
    end if
  end if
else
  nfx=ielem_ineighn(l,i)
  lfx=ielem_qface(l,ngp,ielem_inter_id(ielem_indexf(i)))
  rowfx=bound_offset(nfx)+lfx-1
  cright(1:nof_variables)=boundhir(rowfx,1:nof_variables)
  if (nvt.gt.0) cturbr(1:nvt)=boundhir(rowfx,nof_variables+1:nof_variables+nvt)
end if
if (nvt.gt.0)then
  do iv=1,nvt
    cleft_rot(nof_variables+iv)=cturbl(iv)
    cright_rot(nof_variables+iv)=cturbr(iv)
  end do
end if

end subroutine get_states_bounds2d_ideal




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
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_dim)::cords
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

										if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
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
									if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
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
real,dimension(1:gpu_max_nvar)::srf_speed
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:8,1:gpu_max_dim)::nodes_list
real,dimension(1:gpu_max_dim)::cords
real,dimension(1:6,1:4,1:gpu_max_dim)::elem_listd
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
								    call cordinates2(n,nodes_list,n_node,cords(1:2))

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
integer :: kbasis
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
integer :: kbasis
    integer,intent(in)::n
    integer::i_elem, i_qp, n_qp, i_dof, j_dof, kmaxe,icompwrt,number_of_dog,iconsidered,ixx,number
    real::integ_test,integ_sm1,integ_sm2,phx1,phx2,kron,maxs,mins,higher,phx,integ_mm,x1,y1,z1
    real,dimension(1:num_dg_dofs,1:num_dg_dofs)::totalmm,invmm
    real,dimension(1:idegfree)::basis_vector
    kmaxe = xmpielrank(n)


    !$omp do
    do i_elem = 1, kmaxe
          do i_dof = 1, num_dg_dofs
              do j_dof = 1, num_dg_dofs
                  totalmm(i_dof,j_dof) = zero
                  invmm(i_dof,j_dof) = zero
              end do
          end do
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
                    basis_vector(1:number_of_dog)=basis_rec2d(n,x1,y1,number,ixx,number_of_dog,icompwrt)
                    else
                    number=ielem_iorder(iconsidered)
                    basis_vector(1:number_of_dog)=basis_rec(n,x1,y1,z1,number,ixx,number_of_dog,icompwrt)


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
        do i_dof = 1, num_dg_dofs
            do j_dof = 1, num_dg_dofs
                m_1_val(i_dof,j_dof,i_elem)=invmm(i_dof,j_dof)
            end do
        end do
!



         icompwrt=0







end do
!$omp end do















end subroutine build_mass_matrix

subroutine compmassinv(totalmm,invmm)
!calculate the inverse of the input matrix with gauss-jordan elimination
implicit none
integer :: i,j,k,l,m,irow,num_dofs
real:: big,dum
real,dimension(1:num_dg_dofs,1:num_dg_dofs)::a,b
real,dimension(1:num_dg_dofs,1:num_dg_dofs),intent(in)::totalmm
real,dimension(1:num_dg_dofs,1:num_dg_dofs),intent(inout)::invmm
num_dofs=num_dg_dofs


do i = 1,num_dofs
    do j = 1,num_dofs
        a(i,j) = totalmm(i,j)
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





 do i = 1,num_dofs
     do j = 1,num_dofs
         invmm(i,j)=b(i,j)
     end do
 end do

 end subroutine compmassinv




end module dg_functions
