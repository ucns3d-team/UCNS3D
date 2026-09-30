module profile
use declaration
use library
use basis

implicit none


 contains



 real function linear_init3d(n,pox,poy,poz)
 !> @brief
!> this function initialises the solution for linear advection in 3d,
!> various customisable profiles can be generated and assigned to each initcond code
implicit none
integer,intent(in)::n
real,dimension(1:dimensiona),intent(in)::pox,poy,poz





!coordinates=pox(1),poy(1),poz(1)

if (initcond.eq.0)then
if(((pox(1).ge.0.250).and.(pox(1).le.0.750)).and.((poz(1).ge.0.250).and.(poz(1).le.0.750)))then
	linear_init3d=1.00
else
	linear_init3d=1.00
end if
end if
if (initcond.eq.2)then
linear_init3d=(sin((2.00*pi)*(pox(1))))*&
(sin((2.00*pi)*(poy(1))))*(sin((2.00*pi)*(poz(1))))

end if

end function linear_init3d


 real function linear_init2d(n,pox,poy,poz)
  !> @brief
!> this function initialises the solution for linear advection in 2d,
!> various customisable profiles can be generated and assigned to each initcond code
implicit none
integer,intent(in)::n
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real::aadx,aady,sumf,rd
integer::ixg
!coordinates=pox(1),poy(1),poz(1)


sumf=zero
if (initcond.eq.1)then
 if(((pox(1).ge.0.250).and.(pox(1).le.0.750)).and.((poy(1).ge.0.250).and.(poy(1).le.0.750)))then
	linear_init2d=1.00
 else
	linear_init2d=0.00

end if
end if


if (initcond.eq.3)then

linear_init2d=0.00
if (sqrt(((pox(1)-0.250)**2)+((poy(1)-0.50)**2)).le.0.15)then
rd=(1.00/0.150)*sqrt(((pox(1)-0.250)**2)+((poy(1)-0.50)**2))

linear_init2d=0.250*(1.00+cos(pi*min(rd,1.00)))
end if

if (sqrt(((pox(1)-0.50)**2)+((poy(1)-0.250)**2)).le.0.15)then

rd=(1.00/0.150)*sqrt(((pox(1)-0.50)**2)+((poy(1)-0.250)**2))
linear_init2d=1.00-rd
end if

    if (sqrt(((pox(1)-0.50)**2)+((poy(1)-0.750)**2)).le.0.15)then

    rd=(1.00/0.150)*sqrt(((pox(1)-0.50)**2)+((poy(1)-0.750)**2))
	  if ((abs(pox(1)-0.5).ge.0.0250).or.(poy(1).gt.0.85))then

	  linear_init2d=1.00
	  else

	  linear_init2d=0.00

	  end if
    end if

end if



if (initcond.eq.5)then

linear_init2d=0.00
if ((sqrt(((pox(1)-0.50)**2)+((poy(1)-0.50)**2)).gt.0.25).and.(sqrt(((pox(1)-0.50)**2)+((poy(1)-0.50)**2)).lt.0.35))then

linear_init2d=1.0
else
linear_init2d=0.0


end if
end if


if (initcond.eq.2)then
linear_init2d=(sin((2.00*pi)*(pox(1))))*(sin((2.00*pi)*(poy(1))))

!linear_init2d=1.00
end if




if (initcond.eq.100001)then
linear_init2d=1.0
end if

if (initcond.eq.100002)then
linear_init2d=(sin((2.00*pi)*(pox(1))))*(sin((2.00*pi)*(poy(1))))
end if

if (initcond.eq.100003)then
linear_init2d=0.00
if (sqrt(((pox(1)-0.250)**2)+((poy(1)-0.50)**2)).le.0.15)then
rd=(1.00/0.150)*sqrt(((pox(1)-0.250)**2)+((poy(1)-0.50)**2))

linear_init2d=0.250*(1.00+cos(pi*min(rd,1.00)))
end if

if (sqrt(((pox(1)-0.50)**2)+((poy(1)-0.250)**2)).le.0.15)then

rd=(1.00/0.150)*sqrt(((pox(1)-0.50)**2)+((poy(1)-0.250)**2))
linear_init2d=1.00-rd
end if

    if (sqrt(((pox(1)-0.50)**2)+((poy(1)-0.750)**2)).le.0.15)then

    rd=(1.00/0.150)*sqrt(((pox(1)-0.50)**2)+((poy(1)-0.750)**2))
	  if ((abs(pox(1)-0.5).ge.0.0250).or.(poy(1).gt.0.85))then

	  linear_init2d=1.00
	  else

	  linear_init2d=0.00

	  end if
    end if
end if




end function linear_init2d


 subroutine initialise_euler3d(n,veccos,pox,poy,poz)
 implicit none
  !> @brief
!> this function initialises the solution for euler and navier-stokes equations in 3d,
!> various customisable profiles can be generated and assigned to each initcond code
integer,intent(in)::n
!coordinates=pox,poy,poz
!solution=veccos
!components from dat file gamma,uvel,wvel,vvel,pres,rres
!initcond= profile choice from data file
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::veccos
real,dimension(1:nof_species)::mp_r,mp_a,mp_ie
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx,vhx,amp,dvel,xin,yin,zin
real::rgg,tt1,khi_slope,khi_b,theeta,reeta,rg_ve,rg_tr,rg_chem,rg_density,rg_ev_total,rg_rmix
integer::u_cond1,u_cond2,u_cond3,u_cond4,ix,rg_i,rg_j
real::pr_radius,pr_beta,pr_machnumberfree,pr_pressurefree,pr_temperaturefree,pr_gammafree,pr_rgasfree,pr_xcenter,pr_ylength,pr_xlength,pr_ycenter,pr_densityfree,pr_cpconstant,pr_radiusvar,pr_velocityfree,pr_temperaturevar,drad,rg_tv,rg_ttr0,rg_tve0
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,dimension(1:nof_species)::rg_cvs
real,dimension(1:nof_species)::rg_tvsl,rg_ev
real,dimension(1:nof_variables)::leftv
real::acp,mscp,mvcp,vmcp,bcp,rcp,tcp,vfr,theta1,gammal,mp_pinfl




veccos(:)=zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read initial values from dat file
r1=rres
p1=pres
s1=sqrt((gamma*p1)/(r1))
u1=uvel
v1=vvel
w1=wvel

!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
!internal energy

ie1=((p1)/((gamma-1.00)*r1))

!total energy
e1=r1*(skin1+ie1)

!vector of conserved variables now
if (mrf.eq.1)then
veccos(1)=r1
veccos(2)=r1*u1+1.0e-15
veccos(3)=r1*v1+1.0e-15
veccos(4)=r1*w1+1.0e-15
veccos(5)=e1
else

veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1

end if


if (realgas.eq.1)then
u1=uvel
v1=vvel
w1=wvel
p1=pres
r1=rres


!first build mixture gas constant
rg_rmix=zero
do rg_i=1,nof_species
  rg_rmix=rg_rmix+rg_vf(rg_i)/rg_molm(rg_i)
end do

  rg_rmix=rgs_ru*rg_rmix


  rg_ttr0=rg_ttr    !initial temperature ttr
  rg_tve0=rg_tve    !initial temperature vib



  ! translational-rotational internal energy
rg_tr = 0.00
    do rg_i = 1, nof_species
      if (rg_i <= 3) then
        rg_cvs(rg_i) = (5.00 / 2.00) *  rgs_ru / rg_molm(rg_i)
      else
        rg_cvs(rg_i) = (3.00 / 2.00) *  rgs_ru / rg_molm(rg_i)
      end if
      rg_tr = rg_tr + rg_vf(rg_i) * rg_cvs(rg_i) * rg_ttr0
    end do

! vibrational energy
    rg_ev_total= 0.00
    do rg_i = 1, 3
      rg_ev_total = rg_ev_total + rg_vf(rg_i)  * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/rg_tve0) - 1.00))
    end do

rg_chem=zero

    ! chemical energy
 do rg_i=1,nof_species
        if (rg_hzero(rg_i).gt.1.0e-12)then
        rg_chem=rg_chem+(rg_vf(rg_i)*rg_hzero(rg_i)/rg_molm(rg_i))
        end if
end do

  !rg_chem=zero  !set it to zero for testing
! kinetic energy


skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))

  veccos(1)=r1
  veccos(2)=r1*u1
  veccos(3)=r1*v1
  veccos(4)=r1*w1
  veccos(5)=r1*(rg_ev_total+rg_tr+rg_chem+skin1)
  veccos(6)=r1*rg_ev_total


  do rg_i=1,nof_species
  veccos(6+rg_i)=r1 * rg_vf(rg_i)
  end do

!   leftv(1:nof_variables)=veccos(1:nof_variables)
!
!   write(280+n,*)"look here", rg_ev_total
!   write(280+n,*)"before conservative"
!   write(280+n,*)leftv(1:nof_variables)
!   call cons2prim(n,leftv,mp_pinfl,gammal)
!   write(280+n,*)"after primitive"
!   write(280+n,*)leftv(1:nof_variables)
!   call prim2cons(n,leftv)
!   write(280+n,*)"final  conservative"
!   write(280+n,*)leftv(1:nof_variables)



end if























if (turbulence.eq.1)then

  if (turbulencemodel.eq.1)then

  veccos(nof_variables+1)=visc*turbinit
  end if
  if (turbulencemodel.eq.2)then

   if (zero_turb_init .eq. 0) then
  veccos(nof_variables+1)=(1.50*(i_turb_inlet*ufreestream)**2)*r1
  veccos(nof_variables+2)=ufreestream/l_turb_inlet
 veccos(nof_variables+2)=(c_mu_inlet**(-0.250))*sqrt(veccos(5))&
			/l_turb_inlet*r1
  end if

  if (zero_turb_init .eq. 1) then
  veccos(nof_variables+1)=zero
  veccos(nof_variables+2)=ufreestream/l_turb_inlet
  end if

  end if


end if
if (passivescalar.gt.0)then

  veccos(nof_variables+turbulenceequations+1:nof_variables+turbulenceequations+passivescalar)=zero

end if




if (initcond.eq.10000)then	!shock density interaction


r1=0.50
p1=0.4127
u1=0.0
v1=0.0
w1=0.0




skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1





end if





if (initcond.eq.95)then	!taylor green initial profile
xin=pox(1)-pi
yin=poy(1)-pi
zin=poz(1)-pi

        if(boundtype.eq.1)then
r1=1.00
w1=0.00
p1=100.00+((r1/16.00)*((cos(2.00*zin))+2.00)*((cos(2.00*xin))+(cos(2.00*yin))))
u1=sin(xin)*cos(yin)*cos(zin)
v1=-cos(xin)*sin(yin)*cos(zin)



else

w1=0.00
p1=(1.00/(gamma*1.25*1.25))+((1.00/16.00)*((cos(2.00*zin))+2.00)*((cos(2.00*xin))+(cos(2.00*yin))))
r1=(p1*(gamma*1.25*1.25))
u1=sin(xin)*cos(yin)*cos(zin)
v1=-cos(xin)*sin(yin)*cos(zin)


end if
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
end if

if (initcond.eq.101)then	!shock density interaction
if (pox(1).lt.-4.00)then
r1=3.85710
u1=2.62940
v1=zero
w1=zero
p1=10.3330
else
r1=(1.00+0.20*sin(5.00*pox(1)))
u1=zero
v1=zero
w1=zero
p1=1
end if
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1





end if





if (initcond.eq.405)then
!test case 4.5 of coralic & colonius

if (pox(1).lt.-0.10)then
mp_r(1)=0.1663157890
mp_r(2)=1.6580
mp_a(1)=0.00
mp_a(2)=1.00
u1=114.490
v1= 0.00
w1=0.00
p1=159060.00


! skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))



r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now

else

!first within bubble region



if (sqrt(((pox(1)+0.050)**2)+((poy(1)-0.050)**2)+((poz(1)-0.050)**2)).le.0.0250)then
mp_r(1)=0.1663157890
mp_r(2)=1.2040
mp_a(1)=0.950
mp_a(2)=0.050
u1=0.00
v1=0.00
w1=0.00
p1=101325

! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now
else



mp_r(1)=0.1663157890
mp_r(2)=1.2040
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00
v1=0.00
w1=0.00
p1=101325


! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if






end if












veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
veccos(6)=mp_r(1)*mp_a(1)
veccos(7)=mp_r(2)*mp_a(2)
veccos(8)=mp_a(1)




end if


if (initcond.eq.470)then




        if (pox(1).le.1.0)then   !post shock concidions
        mp_r(2)=1.00 	    ! water density
        mp_r(1)=1.00 		! air density
        mp_a(2)=0.00 		! water volume fraction (everything is water here)
        mp_a(1)=1.00 		! air volume fraction
        u1=0.0	  	          ! m/s
        v1=0.0
        w1=0.0
        p1=1.0      		! pa


        r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
        mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
        mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
        ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
        skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
        e1=(r1*skin1)+ie1
        !vector of conserved variables now

        else

                if (sqrt(poy(1)**2+poz(1)**2).ge.1.50)then

                mp_r(2)=1.00 	! water density
                mp_r(1)=0.1250 		! air density
                mp_a(2)=0.00 		! water volume fraction (everything is water here)
                mp_a(1)=1.00 		! air volume fraction
                u1=0.0	  	          ! m/s
                v1=0.0
                w1=0.0
                p1=0.1      		! pa


                r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
                mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
                mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
                ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
                skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
                e1=(r1*skin1)+ie1
                !vector of conserved variables now

                else

                mp_r(2)=1.00 	! water density
                mp_r(1)=0.1250 		! air density
                mp_a(2)=1.00 		! water volume fraction (everything is water here)
                mp_a(1)=0.00 		! air volume fraction
                u1=0.0	  	          ! m/s
                v1=0.0
                w1=0.0
                p1=0.1      		! pa


                r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
                mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
                mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
                ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
                skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
                e1=(r1*skin1)+ie1
                !vector of conserved variables now



                end if

        end if



veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
veccos(6)=mp_r(1)*mp_a(1)
veccos(7)=mp_r(2)*mp_a(2)
veccos(8)=mp_a(1)


end if


if (initcond.eq.157)then
!test case 4.5 of coralic & colonius

if (sqrt(((pox(1)-200.0e-6)**2)+((poy(1)-150.0e-6)**2)+((poz(1)-150.0e-6)**2)).le.50.0e-6)then


mp_r(1)=1.225
mp_r(2)=1000.0
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.0
v1=0.00
w1=0.0
p1=100000.00


! skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))



r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now

else

!first within bubble region

if (pox(1).le.100.0e-6)then


p1=35e6
mp_r(1)=1.225
mp_r(2)=1000.0
mp_a(1)=0.00
mp_a(2)=1.00
u1=1647.0
v1=0.00
w1=0.00

else


mp_r(1)=1.225
mp_r(2)=1000.0
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.0
v1=0.00
w1=0.00
p1=100000.00


end if





! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if





veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
veccos(6)=mp_r(1)*mp_a(1)
veccos(7)=mp_r(2)*mp_a(2)
veccos(8)=mp_a(1)




end if

if (initcond.eq.411)then
!example vi paper5.pdf

!gamma_in(1) = 4.4 ! water
!gamma_in(2) = 1.4  ! air
!mp_pinf(1) = 6e8 !water from coralic and colonius or 2.218e8(abgrall203)
!mp_pinf(2) = 0 ! air

if (pox(1).le.0.00660)then
mp_r(2)=1323.650 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1=681.0580	  	! m/s
v1= 0.00
w1=0.00
p1=1.9e9      		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region



if (sqrt(((pox(1)-0.012)**2)+((poy(1)-0.012)**2)+((poz(1)-0.012)**2)).le.0.0030)then
mp_r(2)=1000.000 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=0.00 		! water volume fraction (everything is water here)
mp_a(1)=1.00 		! air volume fraction
u1= 0.00	  		! m/s
v1= 0.00
w1=0.00
p1= 100000    			! pa

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now


else

mp_r(2)=1000.00 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1= 0.00	  		! m/s
v1= 0.00
w1=0.00
p1= 100000    			! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if

end if

veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
veccos(6)=mp_r(1)*mp_a(1)
veccos(7)=mp_r(2)*mp_a(2)
veccos(8)=mp_a(1)

end if














if (initcond.eq.408)then
!test case 4.5 of coralic & colonius

if (pox(1).gt.0.100)then
mp_r(1)=6.03
mp_r(2)=1.6580
mp_a(1)=0.00
mp_a(2)=1.00
u1=-114.490
v1= 0.00
w1=0.00
p1=159060.00


! skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region



if (sqrt(((pox(1)-0.0790)**2)+((poy(1)-0.0350)**2)+((poz(1)-0.0350)**2)).le.(0.03250/2.00))then
mp_r(1)=6.03
mp_r(2)=1.2040
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.00
v1=0.00
w1=0.00
p1=101325

! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now
else



mp_r(1)=6.03
mp_r(2)=1.2040
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00
v1=0.00
w1=0.00
p1=101325


! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if






end if












veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
veccos(6)=mp_r(1)*mp_a(1)
veccos(7)=mp_r(2)*mp_a(2)
veccos(8)=mp_a(1)




end if



if (initcond.eq.103)then



r1=rres
p1=pres
s1=sqrt((gamma*p1)/(r1))
v1=vvel
w1=wvel
if (poy(1).gt.0.00)then
u1=uvel
else

u1=0.00
v1=-3.0
p1=press_outlet
end if








!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
!internal energy

ie1=((p1)/((gamma-1.00)*r1))

!total energy
e1=r1*(skin1+ie1)

!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=r1*w1
veccos(5)=e1
if (turbulence.eq.1)then

  if (turbulencemodel.eq.1)then

  veccos(6)=visc*turbinit
  end if
  if (turbulencemodel.eq.2)then

   if (zero_turb_init .eq. 0) then
  veccos(6)=(1.50*i_turb_inlet*(ufreestream**2))*r1
   veccos(7)=r1*veccos(6)/(10.0e-5*visc)
  end if

  if (zero_turb_init .eq. 1) then
  veccos(6)=(1.50*i_turb_inlet*(ufreestream**2))*r1
   veccos(7)=r1*veccos(6)/(10.0e-5*visc)
  end if

  end if


end if
if (passivescalar.gt.0)then

  veccos(nof_variables+turbulenceequations+1:nof_variables+turbulenceequations+passivescalar)=zero

end if


end if



!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine initialise_euler3d



 subroutine initialise_euler2d(n,veccos,pox,poy,poz)
 implicit none
   !> @brief
!> this function initialises the solution for euler and navier-stokes equations in 2d,
!> various customisable profiles can be generated and assigned to each initcond code
integer,intent(in)::n
!coordinates=pox,poy
!solution=veccos
!components from dat file gamma,uvel,wvel,vvel,pres,rres
!initcond= profile choice from data file
real::acp,mscp,mvcp,vmcp,bcp,rcp,tcp,vfr,theta1,theta405
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx,vhx,amp,dvel,rgg,tt1,khi_slope,khi_b,theeta,reeta,rg_ve,rg_tr,rg_chem,rg_density,rg_ev_total,rg_rmix
real::pr_radius,pr_beta,pr_machnumberfree,pr_pressurefree,pr_temperaturefree,pr_gammafree,pr_rgasfree,pr_xcenter,pr_ylength,pr_xlength,pr_ycenter,pr_densityfree,pr_cpconstant,pr_radiusvar,pr_velocityfree,pr_temperaturevar,drad,rg_tv,rg_ttr0,rg_tve0
integer::u_cond1,u_cond2,u_cond3,u_cond4,ix,rg_i,rg_j
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::veccos
real,dimension(1:nof_variables)::leftv
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,dimension(1:nof_species)::mp_r,mp_a,mp_ie
real,dimension(1:nof_species)::rg_cvs
real,dimension(1:nof_species)::rg_tvsl,rg_ev
real::mp_pinfl,gammal
veccos(:)=zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! read initial values from dat file
r1=rres
p1=pres
s1=sqrt((gamma*p1)/(r1))
u1=uvel
v1=vvel

!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy

ie1=((p1)/((gamma-1.00)*r1))

!total energy
e1=r1*(skin1+ie1)

!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1




if (realgas.eq.1)then
u1=uvel
v1=vvel
p1=pres
r1=rres


!first build mixture gas constant
rg_rmix=zero
do rg_i=1,nof_species
  rg_rmix=rg_rmix+rg_vf(rg_i)/rg_molm(rg_i)
end do

  rg_rmix=rgs_ru*rg_rmix


  rg_ttr0=rg_ttr    !initial temperature ttr
  rg_tve0=rg_tve    !initial temperature vib



  ! translational-rotational internal energy
rg_tr = 0.00
    do rg_i = 1, nof_species
      if (rg_i <= 3) then
        rg_cvs(rg_i) = (5.00 / 2.00) *  rgs_ru / rg_molm(rg_i)
      else
        rg_cvs(rg_i) = (3.00 / 2.00) *  rgs_ru / rg_molm(rg_i)
      end if
      rg_tr = rg_tr + rg_vf(rg_i) * rg_cvs(rg_i) * rg_ttr0
    end do

! vibrational energy
    rg_ev_total= 0.00
    do rg_i = 1, 3
      rg_ev_total = rg_ev_total + rg_vf(rg_i)  * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/rg_tve0) - 1.00))
    end do

rg_chem=zero

    ! chemical energy
 do rg_i=1,nof_species
        if (rg_hzero(rg_i).gt.1.0e-12)then
        rg_chem=rg_chem+(rg_vf(rg_i)*rg_hzero(rg_i)/rg_molm(rg_i))
        end if
end do

  !rg_chem=zero  !set it to zero for testing
! kinetic energy


skin1=(oo2)*((u1**2)+(v1**2))
r1=p1/(rg_rmix*rg_ttr0)

  veccos(1)=r1
  veccos(2)=r1*u1
  veccos(3)=r1*v1
  veccos(4)=r1*(rg_ev_total+rg_tr+rg_chem+skin1)
  veccos(5)=r1*rg_ev_total


  do rg_i=1,nof_species
  veccos(5+rg_i)=r1 * rg_vf(rg_i)
  end do

  leftv(1:nof_variables)=veccos(1:nof_variables)

!   write(280+n,*)"look here", rg_ev_total
!   write(280+n,*)"before conservative"
!   write(280+n,*)leftv(1:nof_variables)
!   call cons2prim(n,leftv,mp_pinfl,gammal)
!   write(280+n,*)"after primitive"
!   write(280+n,*)leftv(1:nof_variables)
!   call prim2cons(n,leftv)
!   write(280+n,*)"final  conservative"
!   write(280+n,*)leftv(1:nof_variables)



end if


if (turbulence.eq.1)then

  if (turbulencemodel.eq.1)then

  veccos(nof_variables+1)=visc*turbinit
  end if
  if (turbulencemodel.eq.2)then

   if (zero_turb_init .eq. 0) then
  veccos(nof_variables+1)=(1.50*(i_turb_inlet*ufreestream)**2)*r1
  veccos(nof_variables+2)=ufreestream/l_turb_inlet
 veccos(nof_variables+2)=(c_mu_inlet**(-0.250))*sqrt(veccos(5))&
			/l_turb_inlet*r1
  end if

  if (zero_turb_init .eq. 1) then
  veccos(nof_variables+1)=zero
  veccos(nof_variables+2)=ufreestream/l_turb_inlet
  end if

  end if


end if
if (passivescalar.gt.0)then

  veccos(4+turbulenceequations+1:nof_variables+turbulenceequations+passivescalar)=zero

end if












if (initcond.eq.95)then	!taylor green initial profile
if ((poy(1).ge.0.250).and.(poy(1).le.0.750))then

r1=2.00
u1=-0.50
v1=0.010*sin(2.00*pi*(pox(1)-0.5))
p1=2.5


end if

if ((poy(1).lt.0.250).or.(poy(1).gt.0.750))then
r1=1.00
u1=0.50
v1=0.010*sin(2.00*pi*(pox(1)-0.5))
p1=2.5
end if



!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
end if


if (initcond.eq.75)then	!taylor green initial profile

if (((pox(1).ge.(2.00)).and.(pox(1).le.(4.00))).and.((poy(1).ge.(2.00)).and.(poy(1).le.(4.00))))then
r1=1.00
v1=0.0
u1=0.0
p1=10.00
else
r1=0.20
v1=0.0
u1=0.0
p1=1.00


end if


!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
end if




if (initcond.eq.25)then	!taylor green initial profile

if (sqrt(((pox(1)-1.00)**2)+((poy(1)-1.00)**2)).le.0.40)then

r1=1.00
v1=0.0
u1=0.0
p1=1.00


else
r1=0.1250
v1=0.0
u1=0.0
p1=0.10



end if


!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
end if

if (initcond.eq.31)then	!taylor green initial profile

if (sqrt(((pox(1)-1.00)**2)+((poy(1)-1.00)**2)).le.0.40)then

r1=1.00
v1=0.0
u1=0.0
p1=1.00


else
r1=0.1250
v1=0.0
u1=0.0
p1=0.10



end if


!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
end if







IF (initcond.eq.65)THEN	!vortex evolution



! r1=1.0
! u1=uvel
! v1=vvel
! p1=1.00
! !rgg=((pox(1)**2)+(poy(1)**2))
! rgg=(((pox(1)-5.00)**2)+((poy(1)-5.00)**2))
! u1=u1+(((5.00/(2.00*pi)))*(exp(((1.00-rgg)/(2.00))))*(5.00-poy(1)))
! v1=v1+(((5.00/(2.00*pi)))*(exp(((1.00-rgg)/(2.00))))*(pox(1)-5.00))
! tt1=1.00-(((gamma-1)*(25.00/(8.00*pi**2*gamma)))*exp(1.00-rgg))
! r1=tt1**(1.00/(gamma-1.00))
! P1=R1*tt1
!
!
!
! !KINETIC ENERGY FIRST!
! SKIN1=(OO2)*((U1**2)+(V1**2))
! !INTERNAL ENERGY
! IE1=((P1)/((GAMMA-1.0)*R1))
! !TOTAL ENERGY
! E1=(P1/(GAMMA-1))+(R1*SKIN1)
! !VECTOR OF CONSERVED VARIABLES NOW
! VECCOS(1)=R1
! VECCOS(2)=R1*U1
! VECCOS(3)=R1*V1
! VECCOS(4)=E1

khi_slope=15.00
khi_b=tanh(khi_slope*(poy(1)-1)+7.50)-tanh(khi_slope*(poy(1)-1)-7.50)


R1=0.50+0.750*khi_b
u1=0.5*(khi_b-1.0)
v1=0.1*sin(2.00*pi*(pox(1)-1.00))
p1=1.00




!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY
IE1=((P1)/((GAMMA-1.0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF










if (initcond.eq.100)then	!taylor green initial profile
acp=0.0750
bcp=0.1750
mscp=1.50
mvcp=0.70
vmcp=sqrt(gamma)*mvcp

rcp=sqrt((pox(1)-0.250)**2+(poy(1)-0.50)**2)

if (rcp.le.acp)then
vfr=vmcp*rcp/acp
tcp=((rcp-0.175)*(gamma-1.00)*(vfr**2)/(rcp*gamma))+1.00
p1=tcp**(gamma/(gamma-1.00))
r1=tcp**(1.00/(gamma-1.00))
theta1=atan((poy(1)-0.50)/(pox(1)-0.250))
v1=vfr*cos(theta1)
u1=1.50*sqrt(gamma)-vfr*sin(theta1)



else

 if ((rcp.ge.acp).and.(rcp.le.bcp))then
vfr=vmcp*(acp/(acp**2-bcp**2))*(rcp-((bcp**2)/rcp))

vfr=vmcp*rcp/acp
tcp=((rcp-0.175)*(gamma-1.00)*(vfr**2)/(rcp*gamma))+1.00
p1=tcp**(gamma/(gamma-1.00))
r1=tcp**(1.00/(gamma-1.00))
theta1=atan((poy(1)-0.50)/(pox(1)-0.250))
v1=vfr*cos(theta1)
u1=1.50*sqrt(gamma)-vfr*sin(theta1)








else
vfr=0.00
r1=1.00
u1=1.50*sqrt(gamma)
v1=0.0
p1=1.00
end if
end if

if (pox(1).gt.0.5)then
r1=(9.00*gamma+9.00)/(9.00*gamma-1.00)
u1=sqrt(gamma)*((9.00*gamma-1.00)/(6.00*(gamma+1.00)))
p1=(7.00*gamma+2.00)/(2.00*gamma+2.00)
v1=0.00
end if



!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
end if


if (initcond.eq.101)then	!shock density interaction
if (pox(1).lt.-4.00)then
r1=3.85710
u1=2.62940
v1=zero

p1=10.3330
else
r1=(1.00+0.20*sin(5.00*pox(1)))
u1=zero
v1=zero

p1=1
end if

!kinetic energy first!
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
end if


if (initcond.eq.102)then	!shock density interaction
if (pox(1).lt.((1.00/6.00)+(poy(1)/(sqrt(3.00)))))then
r1=8.00
u1=8.25*cos(pi/6.00)
v1=-8.25*sin(pi/6.00)
p1=116.5
else
r1=1.40
u1=zero
v1=zero
p1=1.00
end if
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1





end if

if (initcond.eq.104)then	!dmr_domain1
if (pox(1).lt.(1.00/6.00))then
r1=8.00
u1=8.250
v1=0.00
p1=116.5
else
r1=1.40
u1=zero
v1=zero
p1=1.00
end if
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1





end if





if (initcond.eq.30)then	!shock density interaction
if (pox(1).le.zero)then
if (poy(1).le.zero)then
r1=0.138
u1=1.206
v1=1.206
p1=0.029
end if
if (poy(1).gt.zero)then
r1=0.5323
u1=1.206
v1=0.0
p1=0.3
end if
end if
if (pox(1).gt.zero)then
if (poy(1).le.zero)then
r1=0.5323
u1=0.0
v1=1.206
p1=0.3
end if
if (poy(1).gt.zero)then
r1=1.5
u1=0.0
v1=0.0
p1=1.5
end if
end if
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1





end if

if (initcond.eq.3)then	!shock density interaction
if (pox(1).ge.zero)then
r1=1.11750
u1=0.00
v1=0.00
p1=95000.00
else
r1=1.7522
u1=166.34345
v1=0.00
p1=180219.750
end if
skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1





end if





if (initcond.eq.123)then	!fast inviscid vortex evolution



pr_radius = 0.005
pr_beta = 0.2
pr_machnumberfree = 0.05
pr_pressurefree = 100000.0
pr_temperaturefree = 300.0
pr_gammafree = 1.4
pr_rgasfree = 287.15
pr_xcenter = 0.05 ; pr_ycenter = 0.05
pr_xlength = 0.1 ;  pr_ylength = 0.1

pr_densityfree = pr_pressurefree/ (pr_rgasfree*pr_temperaturefree)
 pr_cpconstant = pr_rgasfree*(pr_gammafree/(pr_gammafree-1.0))

pr_radiusvar = sqrt( (pox(1)-pr_xcenter)**2 + (poy(1)-pr_ycenter)**2 )/ pr_radius

pr_velocityfree = pr_machnumberfree * sqrt(pr_gammafree*pr_rgasfree*pr_temperaturefree)


u1 = pr_velocityfree * (1.0 - (pr_beta* ((poy(1)-pr_ycenter)/pr_radius )*exp((-pr_radiusvar**2)/2)  ) )
v1 = pr_velocityfree * ((pr_beta* ((poy(1)-pr_ycenter)/pr_radius )*exp((-pr_radiusvar**2)/2)  ) )
pr_temperaturevar = pr_temperaturefree - ( ((pr_velocityfree**2 *pr_beta**2.0)/(2.0 * pr_cpconstant))*exp(-pr_radiusvar**2) )

r1 = pr_densityfree * (pr_temperaturevar/pr_temperaturefree)**(1.0/(pr_gammafree-1.0))

 p1 = r1*pr_rgasfree*pr_temperaturevar
!
! w1=0.00
! p1=100.00+((r1/16.00)*((2.00*(cos(2.00*pox(1))))+(cos(2.00*poy(1)))-2.00))
! !uu=(sqrt((gamma*p1)/(r1)))*0.28
! !p1=1.0
! u1=sin(pox(1))*cos(poy(1))
! v1=-cos(pox(1))*sin(poy(1))
! !kinetic energy first!
 skin1=(oo2)*((u1**2)+(v1**2))
! !internal energy
 ie1=((p1)/((pr_gammafree-1.00)*r1))
! !total energy
 e1=(p1/(pr_gammafree-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1




end if



if (initcond.eq.124)then	!fast inviscid vortex evolution
pr_radius = 0.005
pr_beta = 0.2
pr_machnumberfree = 0.2
pr_pressurefree = 100000.0
pr_temperaturefree = 300.0
pr_gammafree = 1.4
pr_rgasfree = 287.15
pr_xcenter = 0.05 ; pr_ycenter = 0.05
pr_xlength = 0.1 ;  pr_ylength = 0.1

pr_densityfree = pr_pressurefree/ (pr_rgasfree*pr_temperaturefree)
 pr_cpconstant = pr_rgasfree*(pr_gammafree/(pr_gammafree-1.0))

pr_radiusvar = sqrt( (pox(1)-pr_xcenter)**2 + (poy(1)-pr_ycenter)**2 )/ pr_radius

pr_velocityfree = pr_machnumberfree * sqrt(pr_gammafree*pr_rgasfree*pr_temperaturefree)


u1 = pr_velocityfree * (1.0 - (pr_beta* ((poy(1)-pr_ycenter)/pr_radius )*exp((-pr_radiusvar**2)/2)  ) )
v1 = pr_velocityfree * ((pr_beta* ((poy(1)-pr_ycenter)/pr_radius )*exp((-pr_radiusvar**2)/2)  ) )
pr_temperaturevar = pr_temperaturefree - ( ((pr_velocityfree**2 *pr_beta**2.0)/(2.0 * pr_cpconstant))*exp(-pr_radiusvar**2) )

r1 = pr_densityfree * (pr_temperaturevar/pr_temperaturefree)**(1.0/(pr_gammafree-1.0))

 p1 = r1*pr_rgasfree*pr_temperaturevar
!
! w1=0.00
! p1=100.00+((r1/16.00)*((2.00*(cos(2.00*pox(1))))+(cos(2.00*poy(1)))-2.00))
! !uu=(sqrt((gamma*p1)/(r1)))*0.28
! !p1=1.0
! u1=sin(pox(1))*cos(poy(1))
! v1=-cos(pox(1))*sin(poy(1))
! !kinetic energy first!
 skin1=(oo2)*((u1**2)+(v1**2))
! !internal energy
 ie1=((p1)/((pr_gammafree-1.00)*r1))
! !total energy
 e1=(p1/(pr_gammafree-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1




end if






if (initcond.eq.401)then

if (sqrt(((pox(1)-1.00)**2)+((poy(1)-1.00)**2)).le.0.40)then
mp_r(1)=1.00
mp_r(2)=0.125
mp_a(1)=1.00
mp_a(2)=0.00
u1=zero
v1=zero
p1=1.00

else
mp_r(1)=1.00
mp_r(2)=0.125
mp_a(1)=0.00
mp_a(2)=1.00
u1=zero
v1=zero
p1=0.10

end if




r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now




veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)




end if


if (initcond.eq.402)then

if ((pox(1).ge.0.25).and.(pox(1).lt.0.75))then
mp_r(1)=10.00
mp_r(2)=1.00
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.50
v1=0.00
p1=1.00/1.40

else
mp_r(1)=10.00
mp_r(2)=1.00
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.50
v1=0.00
p1=1.00/1.40

end if




r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now




veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)




end if



if (initcond.eq.403)then
!test case 4.1 of wang, deiterding, pan and ren


mp_r(1)=7.00
mp_r(2)=1.00
mp_a(1)=0.50+0.250*(sin(pi*((pox(1)-1.00)+0.50)))
mp_a(2)=1.00-mp_a(1)
u1=1.00
v1=0.00
p1=1.00/1.40






r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now




veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)




end if




if (initcond.eq.157)then









if (sqrt(((pox(1)-200.0e-6)**2)+((poy(1)-150.0e-6)**2)).le.50.0e-6)then

mp_r(1)=1.225
mp_r(2)=1000.0
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.0
v1=0.00
p1=100000.00


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1

else




if (pox(1).le.100.0e-6)then


p1=35e6
mp_r(1)=1.225
mp_r(2)=1000.0
mp_a(1)=0.00
mp_a(2)=1.00
u1=1647
v1=0.00





else


mp_r(1)=1.225
mp_r(2)=1000.0
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.0
v1=0.00
p1=100000.00


end if

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1


end if






veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)






end if





if (initcond.eq.405)then
!test case 4.5 of coralic & colonius

drad=sqrt(((pox(1)+0.050)**2)+((poy(1)-0.050)**2))

theta405 = atan2(poy(1)-0.050, pox(1)+0.050)



if (pox(1).lt.-0.10)then
mp_r(1)=0.166315789
mp_r(2)=1.658
mp_a(1)=0.00
mp_a(2)=1.00
u1=114.490
v1= 0.00
p1=159060.00


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region

if (drad .le. (0.0250 + A405*cos(real(nof_perturbations405)*theta405 + 0.00))) then


mp_r(1)=0.1663157890
mp_r(2)=1.2040
mp_a(1)=0.950
mp_a(2)=0.050
u1=0.00
v1=0.00
p1=101325

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now
else



mp_r(1)=0.166315789
mp_r(2)=1.204
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00
v1=0.00
p1=101325


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if






end if












veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)




end if


if (initcond.eq.422)then
!test case 4.5 of coralic & colonius

if (pox(1).le.0.05)then

mp_r(1)=1000.00 	! water density
mp_r(2)=3.850 		! gas density
mp_a(1)=0.00
mp_a(2)=1.00
u1=567.30	  	! m/s
v1= 0.00
w1=0.00
p1=664000.00    		! pa

else
if ((sqrt(((pox(1)-0.0576)**2)+((poy(1)-0.0576)**2)).le.0.00480))then
! if ((sqrt(((pox(1)-0.0576)**2)+((poy(1)-0.0576)**2)+((poz(1)-0.0576)**2))).le.0.00480)then

mp_r(1)=1000.00 	! water density
mp_r(2)=1.2 		! gas density
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.00	  	! m/s
v1= 0.00
w1=0.00
p1=101000   		! pa


else

mp_r(1)=1000.00 	! water density
mp_r(2)=1.20 		! gas density
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00	  	! m/s
v1= 0.00
w1=0.00
p1=101000   		! pa



end if
end if




r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1







veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)




end if




if (initcond.eq.406)then
!test case 4.5 of coralic & colonius

if (pox(1).gt.0.10)then
mp_r(1)=0.1663157890
mp_r(2)=1.6580
mp_a(1)=0.00
mp_a(2)=1.00
u1=-114.490
v1= 0.00
p1=159060.00


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region
u_cond1=0;u_cond2=0;u_cond3=0; u_cond4=0

 if (((pox(1).ge.0.010).and.(pox(1).le.0.02)).or.((pox(1).ge.0.030).and.(pox(1).le.0.04)))then
    u_cond1=1
 end if
 if ((poy(1).ge.0.050).and.(poy(1).le.0.075))then
    u_cond2=1
 end if


 if ((poy(1).le.0.050).and.(poy(1).gt.0.02))then
    u_cond3=1
 end if
 if (u_cond3.eq.1)then
 if ((sqrt(((pox(1)-0.0250)**2)+((poy(1)-0.050)**2)).le.0.0150).and.(sqrt(((pox(1)-0.0250)**2)+((poy(1)-0.050)**2)).ge.0.0050)) then
    u_cond4=1
 end if
 end if


 if (((u_cond1.eq.1).and.(u_cond2.eq.1)).or.(u_cond4.eq.1))then

mp_r(1)=0.1663157890
mp_r(2)=1.2040
mp_a(1)=0.950
mp_a(2)=0.050
u1=0.00
v1=0.00
p1=101325

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now
else



mp_r(1)=0.1663157890
mp_r(2)=1.2040
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00
v1=0.00
p1=101325


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if






end if












veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)




end if





if (initcond.eq.408)then
!test case 4.5 of coralic & colonius

if (pox(1).gt.0.100)then
mp_r(1)=6.03
mp_r(2)=1.6580
mp_a(1)=0.00
mp_a(2)=1.00
u1=-114.490
v1= 0.00
w1=0.00
p1=159060.00


! skin1=(oo2)*((u1**2)+(v1**2)+(w1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region



if (sqrt(((pox(1)-0.0790)**2)+((poy(1)-0.0350)**2)).le.(0.03250/2.00))then
mp_r(1)=6.03
mp_r(2)=1.2040
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.00
v1=0.00
w1=0.00
p1=101325.00

! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now
else



mp_r(1)=6.03
mp_r(2)=1.2040
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00
v1=0.00
p1=101325.00


! skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if



end if



veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)


end if










if (initcond.eq.410)then
!test case 6.2 of the euler equations for multiphase compressible flow in conservation form simulation of shock–bubble interactions, robin k. s. hankin

!gamma_in(1) = 7.15 ! water
!gamma_in(2) = 1.4  ! air
!mp_pinf(1) = 3.43e8 !water from coralic and colonius or 2.218e8(abgrall203)
!mp_pinf(2) = 0 ! air

if (pox(1).le.-0.0070)then
mp_r(1)=1225.60 ! water density
mp_r(2)=1.20 ! air density
mp_a(1)=1.00 ! water volume fraction (everything is water here)
mp_a(2)=0.00 ! air volume fraction
u1=542.760	  ! m/s
v1= 0.00
p1=1.6e9      ! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region



if (sqrt(((pox(1))**2)+((poy(1))**2)).le.0.0030)then
mp_r(1)=1225.60 ! water density
mp_r(2)=1.20 ! air density
mp_a(1)=0.00 ! water volume fraction (everything is water here)
mp_a(2)=1.00 ! air volume fraction
u1= 0.00	  ! m/s
v1= 0.00
p1= 101325    ! pa

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now
else



mp_r(1)=1000.00 ! water density
mp_r(2)=1.20 ! air density
mp_a(1)=1.00 ! water volume fraction (everything is water here)
mp_a(2)=0.00 ! air volume fraction
u1= 0.00     ! m/s
v1= 0.00
p1=101325     ! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if

end if

veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if




if (initcond.eq.411)then
!example vi paper5.pdf

!gamma_in(1) = 4.4 ! water
!gamma_in(2) = 1.4  ! air
!mp_pinf(1) = 6e8 !water from coralic and colonius or 2.218e8(abgrall203)
!mp_pinf(2) = 0 ! air

if (pox(1).le.0.00660)then
mp_r(2)=1323.650 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1=681.0580	  	! m/s
v1= 0.00
p1=1.9e9      		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

!first within bubble region



if (sqrt(((pox(1)-0.012)**2)+((poy(1)-0.012)**2)).le.0.0030)then
mp_r(2)=1000.000 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=0.00 		! water volume fraction (everything is water here)
mp_a(1)=1.00 		! air volume fraction
u1= 0.00	  		! m/s
v1= 0.00
p1= 100000    			! pa

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now


else

mp_r(2)=1000.00 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1= 0.00	  		! m/s
v1= 0.00
p1= 100000    			! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1

!vector of conserved variables now
end if

end if

veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if


if (initcond.eq.412)then
!example 4.3 paper2.pdf

!gamma_in(1) = 1.4 ! air
!gamma_in(2) = 5.5  ! water
!mp_pinf(1) = 0 !water from coralic and colonius or 2.218e8(abgrall203)
!mp_pinf(2) = 1.505 ! air

if (pox(1).le.0.00)then
mp_r(1)=1.241 	! air density
mp_r(2)=0.991 		! water density
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.0	  	! m/s
v1= 0.00
p1=2.753     		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

mp_r(1)=1.241 	! air density
mp_r(2)=0.991 		! water density
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.0	  	! m/s
v1= 0.00
p1=3.059*10e-4     		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

end if

veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if

if (initcond.eq.420)then

mp_r(1)=7.0 	! air density
mp_r(2)=1.0 		! water density
mp_a(1)=0.50+0.25*sin(pi*(((pox(1)-0.50)*2.00)))
mp_a(2)=1.00-mp_a(1)
u1=1.00	  	! m/s
v1= 0.00
p1=1.00/1.40    		! pa




r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now



veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if


if (initcond.eq.421)then


if (pox(1).lt.0.00)then

mp_r(1)=1.241 	! air density
mp_r(2)=0.991 		! water density
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.00	  	! m/s
v1= 0.00
p1=2.7530    		! pa
else
mp_r(1)=1.241 	! air density
mp_r(2)=0.991 		! water density
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00	  	! m/s
v1= 0.00
p1=3.059*(10.0e-4)    		! pa


end if


skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
! mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)*mp_r(1)))
! mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)*mp_r(2)))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))

ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now



veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if

if (initcond.eq.430)then


if (poy(1).lt.0.3860)then

mp_r(1)=3.483 	! air density
mp_r(2)=867 		! water density
mp_a(1)=0.00
mp_a(2)=1.00
u1=0.00	  	! m/s
v1= 2
p1=300000    		! pa
else
mp_r(1)=3.483 	! air density
mp_r(2)=867		! water density
mp_a(1)=1.00
mp_a(2)=0.00
u1=0.00	  	! m/s
v1= 0.00
p1=300000    		! pa


end if


skin1=(oo2)*((u1**2)+(v1**2))
r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
! mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)*mp_r(1)))
! mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)*mp_r(2)))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))

ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now



veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if


if (initcond.eq.133)then

if (pox(1).lt.0.00)then


	p1=195557.25
	r1=p1/(350.50*287.0580)
	u1=168.62
	v1=0.00

	rhc1=r1
	rhc2=u1
	rhc3=v1
	rhc4=p1



	else


	u1=0.00
	v1=0.00
	p1=101325
	r1=p1/(288.150*287.0580)

	end if


skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1








end if


if (initcond.eq.266)then

if (poy(1).gt.1.00)then


	p1=20000
	r1=0.41
	u1=850
	v1=0.00





	else


	u1=0.00
	v1=0.00
	p1=100000
	r1=1.225

	end if


skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1








end if

if (initcond.eq.444)then




if (poy(1).le.0.00025)then   !post shock concidions
mp_r(2)=1323.650 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1=0.0	  	          ! m/s
v1=681.0580
p1=1.9e9      		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else


mp_r(2)=1000.000 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1= 0.00	  		! m/s
v1= 0.00
p1= 100000    			! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1




!first loop bubbles






do ix=1,nof_bubbles

if (sqrt(((pox(1)-bubble_centre(ix,1))**2)+((poy(1)-bubble_centre(ix,2))**2)).le.bubble_radius(ix))then
mp_r(2)=1000.000 	! water density
mp_r(1)=1.0 		! air density
mp_a(2)=0.00 		! water volume fraction (everything is water here)
mp_a(1)=1.00 		! air volume fraction
u1= 0.00	  		! m/s
v1= 0.00
p1= 100000    			! pa

r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
end if
end do


end if

veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)

end if




if (initcond.eq.222)then	!shock density interaction

r1=1.00
p1=1.0e-6
reeta=-1.00
theeta=atan(poy(1)/pox(1))

u1=reeta*cos(theeta)
v1=reeta*sin(theeta)




skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1





end if


if (initcond.eq.10000)then	!shock density interaction


r1=0.50
p1=0.4127
u1=0.0
v1=0.0




skin1=(oo2)*((u1**2)+(v1**2))
!internal energy
ie1=((p1)/((gamma-1.00)*r1))
!total energy
e1=(p1/(gamma-1))+(r1*skin1)
!vector of conserved variables now
veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1





end if

if (initcond.eq.470)then




if (pox(1).lt.1.0)then   !post shock concidions
mp_r(2)=1.00 	! water density
mp_r(1)=1.00 		! air density
mp_a(2)=0.00 		! water volume fraction (everything is water here)
mp_a(1)=1.00 		! air volume fraction
u1=0.0	  	          ! m/s
v1=0.0
p1=1.0      		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

if (poy(1).gt.1.5)then

mp_r(2)=1.00 	! water density
mp_r(1)=0.1250 		! air density
mp_a(2)=0.00 		! water volume fraction (everything is water here)
mp_a(1)=1.00 		! air volume fraction
u1=0.0	  	          ! m/s
v1=0.0
p1=0.1      		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now

else

mp_r(2)=1.00 	! water density
mp_r(1)=0.1250 		! air density
mp_a(2)=1.00 		! water volume fraction (everything is water here)
mp_a(1)=0.00 		! air volume fraction
u1=0.0	  	          ! m/s
v1=0.0
p1=0.1      		! pa


r1=(mp_r(1)*mp_a(1))+(mp_r(2)*mp_a(2))
mp_ie(1)=((p1+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.00)))
mp_ie(2)=((p1+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.00)))
ie1=(mp_ie(1)*mp_a(1))+(mp_ie(2)*mp_a(2))
skin1=(oo2)*((u1**2)+(v1**2))
e1=(r1*skin1)+ie1
!vector of conserved variables now



end if

end if




veccos(1)=r1
veccos(2)=r1*u1
veccos(3)=r1*v1
veccos(4)=e1
veccos(5)=mp_r(1)*mp_a(1)
veccos(6)=mp_r(2)*mp_a(2)
veccos(7)=mp_a(1)



end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine initialise_euler2d
























end module profile
