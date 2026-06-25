module flow_operations
use declaration
use mpiinfo
use transform


implicit none



 contains
 
 

 subroutine mrfswitch(n,iconsidered,facex,pointx,pox,poy)
implicit none

  !> @brief
!> this subroutine  check if the element is on the rotational/stationary reference frame and update the mrf_origin and srf_velocity srf accordingly

#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer, intent(in)::n,iconsidered,facex,pointx
real, dimension(3) ::mrf_origin, mrf_velocity,rotvel
integer:: rotframe_on
integer::ninv,k
real,dimension(1:dimensiona),intent(inout)::pox,poy
!internal variables
real, dimension(3) :: p1p2, pc, popc,po,pgp,poz !element coordinates, roation_axys, cylinder_center_coordinates, vector_element_center, rotational velocity at gaussian points, gausian points coordinates
real :: d1, d2, r1, theta, dpopc,tempdp

po=(pox(1:3))
pgp=(poy(1:3))
!body
do ninv=1,nrotors

pc(1:3)= (point1_gl(ninv,1:3)+point2_gl(ninv,1:3))/2  !center of cylinder
p1p2(1:3)=point2_gl(ninv,1:3)-point1_gl(ninv,1:3)          !axysvector
popc(1:3)=po(1:3)-pc(1:3)              ! vector elelement-centre
dpopc=sqrt((po(1)-pc(1))**2+(po(2)-pc(2))**2+(po(3)-pc(3))**2) !distance between element and center



tempdp=0.0d0

do k=1,3
  tempdp=tempdp+popc(k)*p1p2(k)
end do





theta= acos((tempdp)/(sqrt(popc(1)**2+popc(2)**2+popc(3)**2)*sqrt(p1p2(1)**2+p1p2(2)**2+p1p2(3)**2))) !angle between element vector and axys
d2=  dpopc*abs(cos(theta))
r1=dpopc*abs(sin(theta))
d1=sqrt((point1_gl(ninv,1)-pc(1))**2+(point1_gl(ninv,2)-pc(2))**2+(point1_gl(ninv,3)-pc(3))**2)

if ((d1.ge.d2).and.(r1.le.radius_gl(ninv))) then
   rotframe_on=1
    mrf_origin(1:3)=pc(1:3)
    pox(1:3)=pgp(1:3)-mrf_origin(1:3)
    mrf_velocity(1:3)=mrf_rot_gl(ninv)*(p1p2)/sqrt(p1p2(1)**2+p1p2(2)**2+p1p2(3)**2)
!     srf_velocity(1)=0.0
!     srf_velocity(2)=mrf_rot_gl
!     srf_velocity(3)=0.0
    poy(1:3)=mrf_velocity(1:3)
    call vect_function(pox,poy,rotvel)



    go to 606
else
    mrf_origin(1:3)=0.0
    rotframe_on=0
    mrf_velocity(1:3)=0.0
    rotvel(1:3)=0.0
end if

end do

606 continue

rec_mrf_origin(1:3,iconsidered)=mrf_origin(1:3)
rec_mrf_velocity(1:3,iconsidered)=mrf_velocity(1:3)
rec_rotvel(facex,pointx,1:3,iconsidered)=rotvel(1:3)
rec_mrf(iconsidered)=rotframe_on






end subroutine mrfswitch




subroutine gqp_sol(n,leftv,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(inout) :: leftv(1:nof_variables)
  real::coord_gqp(1:2),rd
  real::pi
  integer,intent(in)::iconsidered,n
  integer::i,j,k,facex,ngp

  pi=4.0d0*atan(1.0d0)
  i=iconsidered
  facex=1
  ngp=1
  coord_gqp(1:2)= rec_qpoints(facex,ngp,1:2,i)




    leftv(1)=0.0d0
if (sqrt(((coord_gqp(1)-0.25d0)**2)+((coord_gqp(2)-0.5d0)**2)).le.0.15)then
rd=(1.0d0/0.15d0)*sqrt(((coord_gqp(1)-0.25d0)**2)+((coord_gqp(2)-0.5d0)**2))

leftv(1)=0.25d0*(1.0d0+cos(pi*min(rd,1.0d0)))
end if

if (sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.25d0)**2)).le.0.15)then

rd=(1.0d0/0.15d0)*sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.25d0)**2))
leftv(1)=1.0d0-rd
end if

    if (sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.75d0)**2)).le.0.15)then

    rd=(1.0d0/0.15d0)*sqrt(((coord_gqp(1)-0.5d0)**2)+((coord_gqp(2)-0.75d0)**2))
	  if ((abs(coord_gqp(1)-0.5).ge.0.025d0).or.(coord_gqp(2).gt.0.85))then

	 leftv(1)=1.0d0
	  else

	  leftv(1)=0.0d0

	  end if
    end if








end subroutine gqp_sol


subroutine fix_conservative_state(leftv)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(inout) :: leftv(1:nof_variables)

  integer :: k, idxe, idxev
  real :: rho, rho_floor, tiny
  real :: rhoe, rhoev, dev
  real :: sumrhoy, scale
  real :: yk
  real,dimension(1:gpu_max_species)::y
  real :: u, v, w, ke
  real :: tmin, tv_max
  real :: echem, cv_mix, etr_min, rhoe_min
  real :: rhoev_cap, sum_viby
!
!   ------------------------
!   tunable constants
!   ------------------------
  rho_floor = 1.0d-14
  tiny      = 1.0d-30
  tmin      = 200.0d0      ! translational temperature floor
  tv_max    = 20000.0d0    ! cap to prevent tv blow-ups (try 12000-20000)

  idxe  = dimensiona + 2
  idxev = dimensiona + 3

  ! ------------------------
  ! 1) density floor
  ! ------------------------
  rho = leftv(1)
  if (rho < rho_floor) then
    rho = rho_floor
    leftv(1) = rho
  end if

  ! ------------------------
  ! 2) vib energy floor
  ! ------------------------
  rhoev = leftv(idxev)
  if (rhoev < 0.0d0) then
    rhoev = 0.0d0
    leftv(idxev) = rhoev
  end if

  ! ------------------------
  ! 3) species positivity + normalization: enforce sum rhoy = rho
  ! ------------------------
  sumrhoy = 0.0d0
  do k = 1, nof_species
    if (leftv(idxev + k) < 0.0d0) leftv(idxev + k) = 0.0d0
    sumrhoy = sumrhoy + leftv(idxev + k)
  end do

  if (sumrhoy > 0.0d0) then
    scale = rho / sumrhoy
    do k = 1, nof_species
      leftv(idxev + k) = leftv(idxev + k) * scale
      y(k) = leftv(idxev + k) / rho
    end do
  else
    ! pathological: assign all mass to species 1
    leftv(idxev + 1) = rho
    y(1) = 1.0d0
    do k = 2, nof_species
      leftv(idxev + k) = 0.0d0
      y(k) = 0.0d0
    end do
  end if

  ! ------------------------
  ! 4) cap vib energy to avoid tv spikes (composition-compatible)
  !     if rhoev too large for available vib species, reduce it and
  !     add the removed energy back into rhoe (conserve total).
  ! ------------------------
  rhoe  = leftv(idxe)
  rhoev = leftv(idxev)

  sum_viby = 0.0d0
  do k = 1, 3
    sum_viby = sum_viby + y(k)
  end do

  if (sum_viby < 1.0d-12) then
    ! no vib-capable mass: force rhoev -> 0, return energy to rhoe
    dev   = rhoev
    rhoev = 0.0d0
    rhoe  = rhoe + dev
  else
    ! cap corresponding to tv_max using your same oscillator formula
    rhoev_cap = 0.0d0
    do k = 1, 3
      if (y(k) > 0.0d0) then
        rhoev_cap = rhoev_cap + rho * y(k) * (rgs_ru / rg_molm(k)) * &
                    (rg_thetag(k) / (exp(rg_thetag(k)/tv_max) - 1.0d0))
      end if
    end do

    if (rhoev > rhoev_cap) then
      dev   = rhoev - rhoev_cap
      rhoev = rhoev_cap
      rhoe  = rhoe + dev   ! conserve energy by moving excess vib -> total
    end if
  end if

  leftv(idxev) = rhoev
  leftv(idxe)  = rhoe

  ! ------------------------
  ! 5) enforce tmin on translational energy via rhoe floor
  !     rhoe >= rho*(cv_mix*tmin + echem + ke) + rhoev
  ! ------------------------

  ! velocities + ke (note: ke here is per-mass)
  u = leftv(2) / rho
  v = leftv(3) / rho
  if (dimensiona == 3) then
    w = leftv(4) / rho
  else
    w = 0.0d0
  end if
  ke = 0.5d0 * (u*u + v*v + w*w)

  ! chemical energy per mass
  echem = 0.0d0
  do k = 1, nof_species
    if (rg_hzero(k) > 0.0d0) then
      echem = echem - y(k) * (rg_hzero(k) / rg_molm(k))
    end if
  end do

  ! cv_mix per mass
  cv_mix = 0.0d0
  do k = 1, nof_species
    if (k <= 3) then
      cv_mix = cv_mix + y(k) * (5.0d0/2.0d0) * (rgs_ru / rg_molm(k))
    else
      cv_mix = cv_mix + y(k) * (3.0d0/2.0d0) * (rgs_ru / rg_molm(k))
    end if
  end do
  if (cv_mix < tiny) cv_mix = tiny

  etr_min  = cv_mix * tmin
  rhoe_min = rho * (etr_min + echem + ke) + rhoev

  if (leftv(idxe) < rhoe_min) leftv(idxe) = rhoe_min



end subroutine fix_conservative_state









subroutine multispecies_mixtures_rg(leftv,mp_mu_mix,mp_ktr_mix,mp_kve,mp_d_eff,mp_htr,mp_hvib,gammal)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(in)::leftv
real,intent(inout)::mp_mu_mix,mp_ktr_mix,mp_kve,gammal
real,dimension(1:nof_species),intent(inout)::mp_htr,mp_hvib,mp_d_eff
real,dimension(1:gpu_max_species)::rgs_htr_i, rgs_hvib_i
real:: rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho,rgs_mu_mix,mp_pinfl
real, dimension(1:gpu_max_species):: rgs_x, rgs_y, rgs_deff, rgs_mu_i, rgs_ktr_i,tmp
real, dimension(1:gpu_max_species,1:gpu_max_species):: rgs_dij
real ::rgs_ktr_mix, rgs_kve,sumx,sumy
real,dimension(1:gpu_max_nvar)::temp_vect1,temp_vect2,temp_vect3,temp_vect4
integer::rg_i,k
real, parameter :: epsy = 1.0d-20   ! floor for tiny/negative noise
real, parameter :: epsx = 1.0d-20   ! floor for tiny/negative noise
real, parameter :: tiny = 1.0d-300  ! protect against division by zero

temp_vect1=zero  !copy active part of left vector
temp_vect2=zero
temp_vect1(1:nof_variables)=leftv(1:nof_variables)
temp_vect2(1:nof_variables)=leftv(1:nof_variables)

call cons2div(n,temp_vect1,mp_pinfl,gammal)

!now this vector holds the following variables rho, u,v,w,ttt,tvib,y_1,y_2,y_3,y_4,y_5 etc

rgs_ttr=temp_vect1(dimensiona+2)
rgs_tv=temp_vect1(dimensiona+3)

rgs_y(1:nof_species)=temp_vect1(dimensiona+4:nof_variables)



! convert to mole-fraction *proportional* values
do k = 1, nof_species
    tmp(k) = rgs_y(k) / rg_molm(k)
end do

! normalise to get actual mole fractions
sumx = 0.d0
do k = 1, nof_species
    sumx = sumx + tmp(k)
end do

if (sumx > tiny) then
    do k = 1, nof_species
        rgs_x(k) = tmp(k) / sumx
    end do
else
    ! fallback if something is seriously wrong
    do k = 1, nof_species
        rgs_x(k) = 1.d0 / nof_species
    end do
end if






call cons2prim(n,temp_vect2,mp_pinfl,gammal)

rgs_p_pa=temp_vect2(dimensiona+2)
rgs_rho=temp_vect2(1)



! rgs_ttr=rgs_ttr
! rgs_tv=rgs_tv


    call compute_real_gas_diffusion(rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho, rgs_x, rgs_y, &
                                     rgs_dij, rgs_deff, rgs_mu_i, rgs_mu_mix, rgs_ktr_i, rgs_ktr_mix, rgs_kve, &
                                     rgs_htr_i, rgs_hvib_i)

       mp_mu_mix=rgs_mu_mix
       mp_ktr_mix=rgs_ktr_mix
       mp_kve=rgs_kve
       mp_d_eff(1:nof_species)=rgs_deff(1:nof_species)
       mp_htr(1:nof_species)=rgs_htr_i(1:nof_species)
       mp_hvib(1:nof_species)=rgs_hvib_i(1:nof_species)



end subroutine multispecies_mixtures_rg






  subroutine compute_real_gas_diffusion( rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho, rgs_x, rgs_y, &
                                       rgs_dij, rgs_deff, rgs_mu_i, rgs_mu_mix, &
                                       rgs_ktr_i, rgs_ktr_mix, rgs_kve, &
                                       rgs_htr_i, rgs_hvib_i )
  implicit none

#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho
  real, intent(in)  :: rgs_x(1:gpu_max_species), rgs_y(1:gpu_max_species)
  real, intent(out) :: rgs_dij(1:gpu_max_species,1:gpu_max_species), rgs_deff(1:gpu_max_species)
  real, intent(out) :: rgs_mu_i(1:gpu_max_species), rgs_mu_mix
  real, intent(out) :: rgs_ktr_i(1:gpu_max_species), rgs_ktr_mix, rgs_kve
  real, intent(out) :: rgs_htr_i(1:gpu_max_species), rgs_hvib_i(1:gpu_max_species)

  ! ---- binary diffusion ----
  call compute_binary_diffusion( rgs_ttr, rgs_p_pa, rgs_dij )

  ! ---- mixture-average diffusion coefficient ----
  call compute_effective_diffusion( rgs_x, rgs_dij, rgs_deff )

  ! ---- species viscosity (blottner) ----
  call blottner_mu_species( rgs_ttr, rgs_mu_i )

  ! ---- mixture viscosity (wilke) ----
  call wilke_mixture_viscosity( rgs_x, rgs_mu_i, rgs_mu_mix )

  ! ---- species thermal conductivity (eucken) ----
  call eucken_ktr_species( rgs_mu_i, rgs_ktr_i )

  ! ---- mixture translational conductivity ----
  call mason_saxena_ktr_mixture( rgs_x, rgs_ktr_i, rgs_ktr_mix )

  ! ---- vibrational conductivity ----
  call vibrational_conductivity( rgs_rho, rgs_y, rgs_deff, rgs_tv, rgs_kve )

  ! ---- species translational and vibrational enthalpy ----
  call htr_air5 ( rgs_ttr, rgs_htr_i )
  call hvib_air5( rgs_tv , rgs_hvib_i )

end subroutine compute_real_gas_diffusion





  function omega11_neufeld(rgs_tstar) result(rgs_omega11)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in) :: rgs_tstar
  real             :: rgs_omega11
  real, parameter  :: rgs_a=1.06036d0, rgs_b=0.15610d0, rgs_c=0.19300d0, rgs_d=0.47635d0, &
                      rgs_e=1.03587d0, rgs_f=1.52996d0, rgs_g=1.76474d0, rgs_h=3.89411d0
  real :: tstar_eff

  tstar_eff = max(rgs_tstar, 1.0d-6)   ! avoid 0^(-b) and tiny t*

  rgs_omega11 = rgs_a / exp(rgs_b*log(tstar_eff))                        &
              + rgs_c*exp(-rgs_d*tstar_eff)                             &
              + rgs_e*exp(-rgs_f*tstar_eff)                             &
              + rgs_g*exp(-rgs_h*tstar_eff)

end function omega11_neufeld




subroutine compute_binary_diffusion(rgs_t, rgs_p_pa, rgs_dij)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_t, rgs_p_pa
  real, intent(out) :: rgs_dij(1:gpu_max_species,1:gpu_max_species)

  integer :: i, j
  real :: sig_ij, eps_ij, tstar, omega, p_atm, denom
  real, parameter :: tiny = 1d-30

  rgs_dij = 0.0d0

  ! pressure in atm
  p_atm = max(rgs_p_pa / rgs_pa_per_atm, 1d-12)

  do i = 1, nof_species
     do j = i+1, nof_species

        sig_ij = 0.5d0 * (rgs_sigmaa(i) + rgs_sigmaa(j))
        eps_ij = sqrt( rgs_eps_over_k(i) * rgs_eps_over_k(j) )

        ! non-dimensional temperature
        tstar  = max(rgs_t / eps_ij, 1d-6)
        omega  = omega11_neufeld(tstar)

        denom = p_atm * sig_ij**2 * omega * sqrt(1.d0/rgs_mg(i) + 1.d0/rgs_mg(j))
        denom = max(denom, tiny)

        rgs_dij(i,j) = (0.001858d0 * rgs_t*sqrt(rgs_t)) / denom * rgs_cm2s_to_m2s
        rgs_dij(j,i) = rgs_dij(i,j)

     end do
  end do

end subroutine compute_binary_diffusion










 subroutine compute_effective_diffusion(rgs_x, rgs_dij, rgs_deff)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_x(1:gpu_max_species)
  real, intent(in)  :: rgs_dij(1:gpu_max_species,1:gpu_max_species)
  real, intent(out) :: rgs_deff(1:gpu_max_species)

  integer :: i, j
  real :: sumj, one_minus_xi, sumx
  real :: xloc(1:gpu_max_species)
  real, parameter :: tiny = 1d-20
  real, parameter :: dmin = 1d-10   ! lower bound on diffusivity [m^2/s]

  ! work on a local copy to enforce positivity and normalisation


xloc = 0.0d0
xloc(1:nof_species) = rgs_x(1:nof_species)

! enforce positivity
do i = 1, nof_species
    if (xloc(i) < 0.d0) xloc(i) = 0.d0
end do

! enforce normalisation
sumx = sum(xloc(1:nof_species))
if (sumx > tiny) then
    xloc(1:nof_species) = xloc(1:nof_species) / sumx
else
    xloc(1:nof_species) = 1.d0 / max(1,nof_species)
end if

  do i = 1, nof_species
     sumj = 0.0d0

     do j = 1, nof_species
        if (j /= i) then
           sumj = sumj + xloc(j) / max(rgs_dij(i,j), tiny)
        end if
     end do

     one_minus_xi = max(1.0d0 - xloc(i), 0.0d0)

     if (one_minus_xi <= tiny .or. sumj <= tiny) then
        ! pure or nearly pure gas: use small floor (essentially no diffusive transport)
        rgs_deff(i) = dmin
     else
        rgs_deff(i) = one_minus_xi / sumj
        rgs_deff(i) = max(rgs_deff(i), dmin)
     end if
  end do

end subroutine compute_effective_diffusion



   subroutine blottner_mu_species(rgs_t, rgs_mu)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_t
  real, intent(out) :: rgs_mu(1:gpu_max_species)
  integer :: i
  real :: tlog

  tlog = log10(max(rgs_t, 50.d0))

  do i = 1, nof_species
     rgs_mu(i) = 1.d-7 * exp(2.302585092994045684d0 * (rgs_ab(i)*tlog*tlog + rgs_bb(i)*tlog + rgs_cb(i)))
     rgs_mu(i) = max(rgs_mu(i), 1d-12)
  end do
end subroutine blottner_mu_species

   subroutine wilke_mixture_viscosity(rgs_x, rgs_mu_i, rgs_mu_mix)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_x(1:gpu_max_species), rgs_mu_i(1:gpu_max_species)
  real, intent(out) :: rgs_mu_mix
  integer :: i, j
  real :: phi_ij, denom

  rgs_mu_mix = 0.d0

  do i = 1, nof_species
     denom = 0.d0

     do j = 1, nof_species
        if (i == j) then
            denom = denom + rgs_x(j)
        else
            phi_ij = ( 1.d0 + sqrt(rgs_mu_i(i)/rgs_mu_i(j)) * sqrt(sqrt(rg_molm(j)/rg_molm(i))) )**2 &
                     / ( sqrt(8.d0) * sqrt(1.d0 + rg_molm(i)/rg_molm(j)) )

            denom = denom + rgs_x(j) * phi_ij
        end if
     end do

     rgs_mu_mix = rgs_mu_mix + rgs_x(i) * rgs_mu_i(i) / max(denom, 1d-20)
  end do

end subroutine wilke_mixture_viscosity

   subroutine cp_tr_species(rgs_cp_tr)
  implicit none
#ifdef gpu
!$omp declare target
#endif
    real, intent(out) :: rgs_cp_tr(1:gpu_max_species)  ! j/kg-k
    integer :: rgs_i
    real :: rgs_rspec
    do rgs_i = 1, nof_species
      rgs_rspec = rgs_ru / rg_molm(rgs_i)
      if (rgs_i <= 3) then
        rgs_cp_tr(rgs_i) = 3.5d0 * rgs_rspec   ! diatomics: 7/2 r
      else
        rgs_cp_tr(rgs_i) = 2.5d0 * rgs_rspec   ! atoms:     5/2 r
      end if
    end do
  end subroutine cp_tr_species

   subroutine eucken_ktr_species(rgs_mu_i, rgs_ktr_i)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_mu_i(1:gpu_max_species)
  real, intent(out) :: rgs_ktr_i(1:gpu_max_species)

  integer :: i
  real :: rspec, cp

  do i = 1, nof_species
     rspec = rgs_ru / rg_molm(i)

     if (i <= 3) then
         cp = 3.5d0 * rspec          ! diatomic
         rgs_ktr_i(i) = rgs_mu_i(i) * (cp + 1.25d0*rspec)
     else
         cp = 2.5d0 * rspec          ! atomic
         rgs_ktr_i(i) = rgs_mu_i(i) * (cp + 1.50d0*rspec)
     end if

  end do

end subroutine eucken_ktr_species

  subroutine mason_saxena_ktr_mixture(rgs_x, rgs_ktr_i, rgs_ktr_mix)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_x(1:gpu_max_species), rgs_ktr_i(1:gpu_max_species)
  real, intent(out) :: rgs_ktr_mix

  integer :: i, j
  real :: psi_ij, denom

  rgs_ktr_mix = 0.d0

  do i = 1, nof_species
     denom = 0.d0

     do j = 1, nof_species
        if (i == j) then
            denom = denom + rgs_x(j)
        else
            psi_ij = (1.d0 + sqrt(rgs_ktr_i(i)/rgs_ktr_i(j)) * sqrt(sqrt(rg_molm(j)/rg_molm(i))) )**2 &
                      / ( sqrt(8.d0) * sqrt(1.d0 + rg_molm(i)/rg_molm(j)) )

            denom = denom + rgs_x(j) * psi_ij
        end if
     end do

     rgs_ktr_mix = rgs_ktr_mix + rgs_x(i) * rgs_ktr_i(i) / max(denom, 1d-20)

  end do

end subroutine mason_saxena_ktr_mixture

function rgs_cv_vibrational_diatomic(rgs_tv, rgs_theta) result(rgs_cv)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in) :: rgs_tv, rgs_theta
  real             :: rgs_cv, rgs_x, ex

  if (rgs_tv <= 1.0d0 .or. rgs_theta <= 0.0d0) then
    rgs_cv = 0.0d0
  else
    rgs_x = rgs_theta / rgs_tv
    if (rgs_x > 60.0d0) then
      ! asymptotic form for large x: cv ∝ x^2 e^(-x)
      rgs_cv = rgs_ru * (rgs_x*rgs_x) * exp(-rgs_x)
    else
      ex     = exp(rgs_x)
      rgs_cv = rgs_ru * (rgs_x*rgs_x) * ex / ( (ex - 1.0d0)**2 )
    end if
  end if
end function rgs_cv_vibrational_diatomic


function rgs_hvib_species(rgs_tv, rgs_theta, rgs_mi) result(rgs_hv)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in) :: rgs_tv, rgs_theta, rgs_mi
  real             :: rgs_hv, rgs_x, rgs_ri, ex

  if (rgs_theta <= 0.0d0 .or. rgs_tv <= 1.0d0) then
    rgs_hv = 0.0d0
    return
  end if

  rgs_ri = rgs_ru / rgs_mi
  rgs_x  = rgs_theta / rgs_tv

  if (rgs_x > 60.0d0) then
     ! asymptotic form: h_v ≈ r_i θ e^(-x)
     rgs_hv = rgs_ri * rgs_theta * exp(-rgs_x)
  else
     ex     = exp(rgs_x)
     rgs_hv = rgs_ri * rgs_theta / (ex - 1.0d0)
  end if

end function rgs_hvib_species





 subroutine hvib_air5(rgs_tv, rgs_hvib_i)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_tv
  real, intent(out) :: rgs_hvib_i(1:gpu_max_species)     ! j/kg
  integer :: rgs_i

  do rgs_i = 1, nof_species
    rgs_hvib_i(rgs_i) = rgs_hvib_species(rgs_tv, rg_thetag(rgs_i), rg_molm(rgs_i))
  end do

end subroutine hvib_air5






subroutine htr_air5(rgs_ttr, rgs_htr_i)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(in)  :: rgs_ttr
  real, intent(out) :: rgs_htr_i(1:gpu_max_species)
  real :: cp(1:gpu_max_species), h0_mass(1:gpu_max_species)
  integer :: i
  real, parameter :: tref = 298.15d0

  call cp_tr_species(cp)

  do i = 1, nof_species
     h0_mass(i)   = rg_hzero(i) / rg_molm(i)
     rgs_htr_i(i) = h0_mass(i) + cp(i) * (rgs_ttr - tref)
  end do

end subroutine htr_air5







subroutine vibrational_conductivity(rgs_rho, rgs_y, rgs_deff, rgs_tv, rgs_kve)
  implicit none
#ifdef gpu
!$omp declare target
#endif

  real, intent(in)  :: rgs_rho, rgs_y(1:gpu_max_species), rgs_deff(1:gpu_max_species), rgs_tv
  real, intent(out) :: rgs_kve

  real :: cv_vib_mass(1:gpu_max_species)
  integer :: i

  do i = 1, nof_species
     if (rg_thetag(i) > 0.0d0) then
        cv_vib_mass(i) = rgs_cv_vibrational_diatomic(rgs_tv, rg_thetag(i)) / rg_molm(i)
     else
        cv_vib_mass(i) = 0.0d0
     end if
  end do

  rgs_kve = 0.0d0
  do i = 1, nof_species
     rgs_kve = rgs_kve + rgs_rho * rgs_y(i) * rgs_deff(i) * cv_vib_mass(i)
  end do

  rgs_kve = max(rgs_kve, 0.0d0)

end subroutine vibrational_conductivity












subroutine multispecies_mixtures(leftv,mp_temp,mp_mu_mix,mp_k_mix,mp_cp_mix,gammal)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(inout)  :: mp_temp    ! temperature in kelvin
  real, intent(in),dimension(1:nof_variables)  :: leftv      ! vector of conserved variables
  real, intent(out) :: mp_mu_mix  ! mixture viscosity [pa.s]
  real, intent(out) :: mp_k_mix   ! mixture thermal conductivity [w/m.k]
  real, intent(out) :: mp_cp_mix  ! mixture cp [j/mol.k]
  real, intent(out) :: gammal     ! total
  real::molar_sum,mass_sum
  integer::i,j,k
  ! constants
  real,dimension(1:gpu_max_species) :: mp_viscl,mp_denom
  real,dimension(1:gpu_max_species) ::mp_mole_fraction,mp_mass_fraction
  real,dimension(1:gpu_max_species) :: mp_cp
  real,dimension(1:gpu_max_species) :: ml_laml,mp_denol
  real,dimension(1:gpu_max_species,1:gpu_max_species) :: mp_phi





  ! molar masses
  !mp_m(1)  = 2.016e-3    ! kg/mol   !defined in file
  !mp_m(2)  = 28.97e-3    ! kg/mol   !defined in file

  !mass sum
  mass_sum=zero
  do i=1,nof_species
      mp_mass_fraction(i)=leftv(nof_variables-nof_species-1+i)
      mass_sum=mass_sum+mp_mass_fraction(i)
  end do
  do i=1,nof_species
        mp_mass_fraction(i)=mp_mass_fraction(i)/mass_sum
  end do



  !mole fraction of each species
  molar_sum=zero
  do i=1,nof_species
      mp_mole_fraction(i)=mp_mass_fraction(i)/mp_m(i)
      molar_sum=molar_sum+mp_mole_fraction(i)
  end do
  do i=1,nof_species
        mp_mole_fraction(i)=mp_mole_fraction(i)/molar_sum
  end do



  !now get the temperature

  call get_temp_janaf(leftv,mp_temp,mp_mass_fraction,mp_mole_fraction)




  ! --- compute individual cp values [j/mol.k]
  do i=1,nof_species
   mp_cp(i)=rgs_ru*(mp_janaf(i,j,1) + mp_janaf(i,j,2)*mp_temp + mp_janaf(i,j,3)*mp_temp**2 + mp_janaf(i,j,4)*mp_temp**3 + mp_janaf(i,j,5)*mp_temp**4)
  end do
  ! --- compute individual viscosities [pa.s]
  do i=1,nof_species
    mp_viscl(i)=mp_brok_a(i) * exp(mp_brok_b(i)*log(mp_temp)) + mp_brok_c(i)
  end do

  ! --- compute interaction terms
  do i = 1, nof_species
    do j = 1, nof_species
      if (i == j) then
        mp_phi(i,j) = 1.0d0
      else
        mp_phi(i,j) = (1.0d0 + sqrt(mp_viscl(i)/mp_viscl(j)) * sqrt(sqrt(mp_m(j)/mp_m(i))))**2 / &
                   sqrt(8.0d0 * (1.0d0 + mp_m(i)/mp_m(j)))
      end if
    end do
  end do

  ! wilke’s mixture viscosity
  mp_mu_mix = 0.0d0
  do i = 1, nof_species
    mp_denom(i) = 0.0d0
    do j = 1, nof_species
      mp_denom(i) = mp_denom(i) + mp_mole_fraction(j) * mp_phi(i,j)
    end do
    mp_mu_mix = mp_mu_mix + mp_mole_fraction(i) * mp_viscl(i) / mp_denom(i)
  end do

  ! compute thermal conductivity via eucken
  do i = 1, nof_species
    ml_laml(i) = (mp_cp(i) + 1.25d0 * rgs_ru) * mp_viscl(i)
  end do


  mp_k_mix = 0.0d0
  do i = 1, nof_species
    mp_denol(i) = 0.0d0
    do j = 1, nof_species
      mp_denol(i) = mp_denol(i) + mp_mole_fraction(j) * mp_phi(i,j)
    end do
    mp_k_mix = mp_k_mix + mp_mole_fraction(i) * ml_laml(i) / mp_denol(i)
  end do

  ! mixture cp
  mp_cp_mix = 0.0d0
  do i = 1, nof_species
    mp_cp_mix = mp_cp_mix + mp_mole_fraction(i) * mp_cp(i)
  end do

  !  mixture gamma
  gammal = mp_cp_mix / (mp_cp_mix - rgs_ru)




end subroutine multispecies_mixtures













subroutine multispecies_temp(leftv,mp_temp)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  real, intent(inout)  :: mp_temp    ! temperature in kelvin
  real, intent(in),dimension(1:nof_variables)  :: leftv      ! vector of conserved variables
  real::molar_sum,mass_sum
  integer::i,j,k
  ! constants
  real,dimension(1:gpu_max_species) :: mp_viscl,mp_denom
  real,dimension(1:gpu_max_species) ::mp_mole_fraction,mp_mass_fraction
  real,dimension(1:gpu_max_species) :: mp_cp
  real,dimension(1:gpu_max_species) :: ml_laml,mp_denol
  real,dimension(1:gpu_max_species,1:gpu_max_species) :: mp_phi





  ! molar masses
  !mp_m(1)  = 2.016e-3    ! kg/mol   !please define in file
  !mp_m(2)  = 28.97e-3    ! kg/mol   !please define in file

  !mass sum
  mass_sum=zero
  do i=1,nof_species
      mp_mass_fraction(i)=leftv(nof_variables-nof_species-1+i)
      mass_sum=mass_sum+mp_mass_fraction(i)
  end do
  do i=1,nof_species
        mp_mass_fraction(i)=mp_mass_fraction(i)/mass_sum
  end do



  !mole fraction of each species
  molar_sum=zero
  do i=1,nof_species
      mp_mole_fraction(i)=mp_mass_fraction(i)/mp_m(i)
      molar_sum=molar_sum+mp_mole_fraction(i)
  end do
  do i=1,nof_species
        mp_mole_fraction(i)=mp_mole_fraction(i)/molar_sum
  end do



  !now get the temperature

  call get_temp_janaf(leftv,mp_temp,mp_mass_fraction,mp_mole_fraction)









end subroutine multispecies_temp




  subroutine normalize_mf(mpc_mp)
    implicit none
#ifdef gpu
!$omp declare target
#endif
    !! ensure mole fractions sum to 1 (robust to roundoff)
    real,dimension(1:nof_species),intent(inout) :: mpc_mp
    real:: s
    integer :: k
    s = zero
    do k=1,nof_species
    s = s + max(zero, mpc_mp(k))
    end do
    if (s > zero) then
      do k=1,nof_species
      mpc_mp(k) = max(zero, mpc_mp(k))/s
      end do
    else
      do k=1,nof_species
      mpc_mp(k) = zero
      end do
    end if
  end subroutine normalize_mf

  subroutine molef_to_massf(mpc_mp,mpc_mass)
    implicit none
#ifdef gpu
!$omp declare target
#endif
    !! convert mole fractions to mass fractions (for mass-specific energy)
    real,dimension(1:nof_species), intent(in)  :: mpc_mp
    real,dimension(1:nof_species), intent(out) :: mpc_mass
    real :: denom
    integer :: s
    denom = zero
    do s=1,nof_species
      denom = denom + mpc_mp(s)*mp_m(s)
    end do
    if (denom <= zero) then
      mpc_mass = zero
    else
      do s=1,nof_species
        mpc_mass(s) = mpc_mp(s)*mp_m(s) / denom
      end do
    end if
  end subroutine molef_to_massf

   subroutine species_cp_h_sensible(mp_temp, is, mp_cpl, mp_hl)
    implicit none
#ifdef gpu
!$omp declare target
#endif
    !! janaf 5-coeffs, sensible enthalpy with h(298.15k)=0 for each species
    !!   cp/r = a1 + a2*t + a3*t^2 + a4*t^3 + a5*t^-2
    !!    h/r = a1*t + a2*t^2/2 + a3*t^3/3 + a4*t^4/4 - a5/t  (then shifted so h(298.15)=0)
    real, intent(in)  :: mp_temp
    integer,  intent(in)  :: is
    real, intent(out) :: mp_cpl, mp_hl
    integer :: r
    real:: rs, mp_tt,mp_t2,mp_t3, mp_invt, mp_invt2, mp_a1,mp_a2,mp_a3,mp_a4,mp_a5, mp_href, mp_tref

    rs = rgs_ru / mp_m(is)
    r  = merge(1,2, mp_temp <= mp_tmid_in(is))
    mp_a1=mp_janaf(1,r,is); mp_a2=mp_janaf(2,r,is); mp_a3=mp_janaf(3,r,is); mp_a4=mp_janaf(4,r,is); mp_a5=mp_janaf(5,r,is)

     mp_tt=mp_temp; mp_t2=mp_tt*mp_tt; mp_t3=mp_t2*mp_tt; mp_invt=1.0d0/mp_tt; mp_invt2=mp_invt*mp_invt

    mp_cpl = rs * ( mp_a1 + mp_a2*mp_tt + mp_a3*mp_t2 + mp_a4*mp_t3 + mp_a5*mp_invt2 )

    ! raw enthalpy
    mp_hl  = rs * ( mp_a1*mp_tt + 0.5d0*mp_a2*mp_t2 + (mp_a3/3.0d0)*mp_t3 + 0.25d0*mp_a4*mp_t2*mp_t2 - mp_a5*mp_invt )

    !only for sensible enthalpies
    ! subtract reference to exclude formation/constant offsets: h(298.15 k) = 0
    mp_tref = 298.15d0
    mp_href = rs * ( mp_a1*mp_tref + 0.5d0*mp_a2*mp_tref**2 + (mp_a3/3.0d0)*mp_tref**3 + 0.25d0*mp_a4*mp_tref**4 - mp_a5/mp_tref )
    mp_hl = mp_hl - mp_href
  end subroutine species_cp_h_sensible

  subroutine mixture_props(mp_temp, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
  implicit none
#ifdef gpu
!$omp declare target
#endif
    !! mixture properties with mass-weighted energies; energies are sensible
    real, intent(in)  :: mp_temp
    real,dimension(1:nof_species), intent(in)  :: mp_mol_x
    real, intent(out) :: mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix
    real, dimension(1:gpu_max_species) :: mpc_mass
    real::mp_cp_s, mp_h_s
    integer :: s
    call molef_to_massf(mp_mol_x,mpc_mass)
    mp_rmix  = zero
    mp_cp_mix= zero
    mp_e_mix = zero
    do s=1,nof_species
      if (mpc_mass(s) <= zero) cycle
      mp_rmix = mp_rmix + mpc_mass(s) * (rgs_ru/mp_m(s))
      call species_cp_h_sensible(mp_temp, s, mp_cp_s, mp_h_s)
      mp_e_mix  = mp_e_mix  + mpc_mass(s) * (mp_h_s - (rgs_ru/mp_m(s))*mp_temp)  ! e = h - rs*t (sensible e)
      mp_cp_mix = mp_cp_mix + mpc_mass(s) * mp_cp_s
    end do
    mp_cv_mix = mp_cp_mix - mp_rmix
  end subroutine mixture_props


  subroutine normalize_x(x)
  implicit none
#ifdef gpu
!$omp declare target
#endif
    !! ensure mole fractions sum to 1 (robust to roundoff)
    real, intent(inout) :: x(nof_species)
    real :: s
    integer :: k
    s = 0.0d0
    do k=1,nof_species; s = s + max(0.0d0, x(k)); end do
    if (s > 0.0d0) then
      do k=1,nof_species; x(k) = max(0.0d0, x(k))/s; end do
    else
      do k=1,nof_species; x(k) = 0.0d0; end do
    end if
  end subroutine normalize_x



  subroutine get_temp_janaf(leftv, mp_temp,mp_mass_fraction,mp_mole_fraction)
   implicit none
#ifdef gpu
!$omp declare target
#endif
    !! recover temperature from conserved variables using mole fractions x
    !!
    !! inputs:
    !!   rho   [kg/m^3], u,v,w [m/s], rhoe [j/m^3], x_in(ns) mole fractions (sum≈1)
    !! optional:
    !!   t_low,t_high [k] bracket (default 150..20000 clipped to species ranges)
    !!   tol residual on |e_mix - e_target| [j/kg] (default ~1e-6*|e_target|+1e-3)
    !!   max_iter (default 60)
    !!
    !! outputs:
    !!   t [k], info: 0 ok; 1 bracket width tol; 2 iter limit; <0 bracketing issue
    real, dimension(1:nof_variables),intent(in):: leftv
    real, intent(out) :: mp_temp
    real,dimension(1:nof_species),intent(in) ::mp_mole_fraction,mp_mass_fraction
    integer:: mp_info,s
    real:: mp_tol,mp_sum,u,v,w
    integer:: mp_max_iter,i,j,k
    real,dimension(1:gpu_max_species)::mp_mol_x
    real:: mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix
    real:: mp_kin, mp_e_tgt, mp_tlo, mp_thi, mp_fa, mp_fb, mp_f, mp_df, mp_at, mp_bt, mp_tn, mp_fn, mp_tguess, mp_atol
    integer  :: mp_it, mp_itmax, mp_s
    logical::mp_bracket


    ! copy/normalize x
    mp_mol_x = zero





  mp_mol_x(1:nof_species)=mp_mole_fraction(1:nof_species)


    call normalize_x(mp_mol_x)

    u=leftv(2)/leftv(1)
    v=leftv(3)/leftv(1)

    if (dimensiona.eq.3)then
    w=leftv(4)/leftv(1)
    else
    w=zero
    end if

    mp_itmax = 50
    mp_kin   = 0.5d0*(u*u + v*v + w*w)
    mp_e_tgt = (leftv(nof_variables-nof_species-1)/leftv(1)) - mp_kin        ! target sensible internal energy [j/kg]

    ! temperature bracket honoring species ranges
    mp_tlo = 150.0d0
    mp_thi = 20000.0d0
    do s=1,nof_species
      if (mp_mol_x(s) <= zero) cycle
      mp_tlo = max(mp_tlo, mp_tlow_in(s))
      mp_thi = min(mp_thi, mp_thigh_in(s))
    end do
    if (mp_tlo >= mp_thi) then
      mp_info = -3; mp_temp = 300.0d0; return
    end if

    ! initial guess from cv around 300 k
    call mixture_props(300.0d0, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
    mp_tguess = 300.0d0 + (mp_e_tgt - mp_e_mix)/max(1.0e-8, mp_cv_mix)
    mp_temp = min(max(mp_tguess, mp_tlo), mp_thi)

    ! establish/verify bracket
    call mixture_props(mp_tlo, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix); mp_fa = mp_e_mix - mp_e_tgt
    call mixture_props(mp_thi, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix); mp_fb = mp_e_mix - mp_e_tgt
    mp_bracket = (mp_fa*mp_fb <= 0.0)
    mp_at = mp_tlo; mp_bt = mp_thi

    do mp_it = 1, mp_itmax
      call mixture_props(mp_temp, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
      mp_f  = mp_e_mix - mp_e_tgt
      mp_df = max(mp_cv_mix, 1.0e-20)

      mp_atol = 1.0e-6*max(1.0d0,abs(mp_e_tgt)) + 1.0e-3
      if (abs(mp_f) <= mp_atol) then
        mp_info = 0; return
      end if

      ! safeguarded newton step
      mp_tn = mp_temp - mp_f/mp_df
      if ((.not. mp_bracket) .or. (mp_tn <= mp_at) .or. (mp_tn >= mp_bt)) mp_tn = 0.5d0*(mp_at + mp_bt)

      call mixture_props(mp_tn, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
      mp_fn = mp_e_mix - mp_e_tgt

      if (mp_fa*mp_fn <= 0.0) then
        mp_bt = mp_tn; mp_fb = mp_fn; mp_bracket = .true.
      else
        mp_at = mp_tn; mp_fa = mp_fn
      end if
      mp_temp = mp_tn

      if (abs(mp_bt - mp_at) <= 1.0e-8*max(1.0,mp_temp)) then
        mp_info = 1; return
      end if
    end do

    mp_info = 2  ! iteration limit
  end subroutine get_temp_janaf
























 
 
 function fluxeval2d(leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar)::fluxeval2d
real,dimension(1:nof_variables),intent(in)::leftv
real::p,u,v,w,e,r,s,gm,skin,ien,pi
r=leftv(1)
u=leftv(2)
v=leftv(3)
p=leftv(4)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)

fluxeval2d(1)=r*u
fluxeval2d(2)=(r*(u**2))+p
fluxeval2d(3)=r*u*v
fluxeval2d(4)=u*(e+p)

end function

function fluxeval3d(leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar)::fluxeval3d
real,dimension(1:nof_variables),intent(in)::leftv
real::p,u,v,w,e,r,s,gm,skin,ien,pi
r=leftv(1)
u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)

fluxeval3d(1)=r*u
fluxeval3d(2)=(r*(u**2))+p
fluxeval3d(3)=r*u*v
fluxeval3d(4)=r*u*w
fluxeval3d(5)=u*(e+p)

end function
 
 
subroutine cons2prim(n,leftv,mp_pinfl,gammal)
implicit none
!> @brief
!> this subroutine transforms one vector of conservative variables to primitive variables
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::temps
real,dimension(1:nof_variables),intent(inout)::leftv
real,intent(inout)::mp_pinfl,gammal
real::oodensity,mp_density,mp_stiff,skinx
real::p_sat,p_tol, rho_g,rho_l, ss_g, ss_l, pp, p_gl, void_frac,p_temp,etr
real,dimension(1:gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:gpu_max_species)::rg_vftemp
real,dimension(1:gpu_max_species)::rg_cvs,rg_cps
real,dimension(1:gpu_max_species)::rg_tvsl,rg_ev
real::rg_ve,rg_tr,rg_chem,rg_density,rg_ev_total,rg_rmix,rg_kin,rg_cvt,rg_cpt
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,tve,ttr,rg_f,rg_df,sum1,sum2,sum3,u,v,w
real:: rhoe, rhoev
real,dimension(1:gpu_max_species):: rho_i,y
real:: ke, etot, echem, rmix,evib
real:: cv_mix, cp_mix,rho
integer :: i, idxe, idxev






rg_maxiter=20

if (nof_variables.gt.1)then     !only then

p_sat =2000
p_tol =10e-5

      if (multispecies.eq.1) then   !multispecies


            if (mp_modelc.eq.0)then
                  sum1=zero;
                  sum2=zero;
                  sum3=zero
                  do rg_i=1,nof_species
                      sum1=sum1+leftv(dimensiona+2+rg_i)
                  end do
                  mp_density=sum1

                  do rg_i=1,nof_species-1
                  sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
                  mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
                  mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
                  end do

                  mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
                  mp_vft(nof_species)=(1.0d0-sum2)
                  sum3=zero
                  do rg_i=1,nof_species
                    sum3=sum3+mp_ar(rg_i)
                    end do
                    gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
                    oodensity=1.0d0/mp_density

                  temps(1)=mp_density
                  temps(2)=leftv(2)*oodensity; u=temps(2)
                  temps(3)=leftv(3)*oodensity;v=temps(3)
                  w=zero
                    if (dimensiona.eq.3)then
                    temps(4)=leftv(4)*oodensity;w=temps(4)
                    end if
                  sum3=zero
                  do rg_i=1,nof_species
                    sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
                    end do

                  mp_stiff=sum3*(gammal-1.0d0)
                  sum2=zero
                  do rg_i=1,nof_species
                    sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
                    end do

                    mp_pinfl=sum2

                    skinx=(u*u)+(v*v)+(w*w)

                  temps(dimensiona+2)=(((gammal-1.0d0))*((leftv(dimensiona+2))-oo2*temps(1)*skinx))-mp_stiff
                  temps(dimensiona+3:nof_variables)=leftv(dimensiona+3:nof_variables)
                  leftv(1:nof_variables)=temps(1:nof_variables)
            end if  !mp modelc



      else

                        if (realgas.eq.1)then

                                temps(:)=0.0d0



!                                  call fix_conservative_state(leftv)

                                !first get total density-correct !
                                ! note:
                                ! total density rho is a primary conserved variable from continuity.
                                ! do not recompute rho as sum(rho_i) here.
                                ! species equations are not perfectly conservative numerically and
                                ! redefining rho severely degrades robustness in hypersonic flows.

                                rho=leftv(1)

                                u=leftv(2)/rho
                                v=leftv(3)/rho
                                if (dimensiona.eq.3)then
                                w=leftv(4)/rho
                                else
                                w=zero
                                end if

                               idxe  = dimensiona + 2          ! index of ρe
                                idxev = dimensiona + 3          ! index of ρe_vib

                                rhoe  = leftv(idxe)
                                rhoev = leftv(idxev)

                                do i = 1, nof_species
                                  rho_i(i) = leftv(idxev + i)   ! species densities ρ_i
                                end do


                                ! -------------------------------
                                ! 2. mass fractions
                                ! -------------------------------
                                do i = 1, nof_species
                                  y(i) = rho_i(i) / rho
                                end do

                                ! -------------------------------
                                ! 3. energies
                                ! -------------------------------
                                ke   = 0.5d0 * (u*u + v*v + w*w)
                                etot = rhoe / rho
                                evib = rhoev / rho

                                ! chemical energy: e_chem = - σ y_i * h°_i / m_i  [j/kg]
                                  echem = 0.0d0
                                  do i = 1, nof_species
                                    if (rg_hzero(i) > 0.0d0) then
                                      echem = echem - y(i) * (rg_hzero(i) / rg_molm(i))
                                    end if
                                  end do

                                  ! translational–rotational internal energy
                                  etr = etot - evib - echem - ke
                                  if (etr < 0.0d0) etr = 1.0d-12   ! safety

                                  ! -------------------------------
                                  ! 4. mixture gas constant rmix
                                  ! -------------------------------
                                  rmix = 0.0d0
                                  do i = 1, nof_species
                                    rmix = rmix + y(i) / rg_molm(i)
                                  end do
                                  rmix = rgs_ru * rmix    ! [j/kg/k]

                                  ! -------------------------------
                                  ! 5. mixture cv, cp, gamma
                                  !    here cv_i = 5/2 r for i<=3 (diatomic),
                                  !               3/2 r for i>3  (monatomic)
                                  ! -------------------------------
                                  cv_mix = 0.0d0
                                  do i = 1, nof_species
                                    if (i <= 3) then
                                      cv_mix = cv_mix + y(i) * (5.0d0/2.0d0) * (rgs_ru / rg_molm(i))
                                    else
                                      cv_mix = cv_mix + y(i) * (3.0d0/2.0d0) * (rgs_ru / rg_molm(i))
                                    end if
                                  end do

                                  cp_mix = cv_mix + rmix
                                  gammal  = cp_mix / cv_mix

                                  ! -------------------------------
                                  ! 6. translational temperature ttr and pressure
                                  !     etr = cv_mix * ttr  → ttr = etr/cv_mix
                                  !     p   = ρ rmix ttr
                                  ! -------------------------------
                                  if (cv_mix > 1.0d-20) then
                                    temps(dimensiona+2) = rho * rmix * (etr / cv_mix)
                                  else
                                    temps(dimensiona+2) = 0.0d0
                                  end if

                                mp_pinfl=zero


                                leftv(1)=rho
                                leftv(2)=u
                                leftv(3)=v
                                if (dimensiona.eq.3)then
                                leftv(4)=w
                                end if
                                leftv(dimensiona+2)=temps(dimensiona+2) !pressure
                                leftv(dimensiona+3)=evib
                                do i = 1, nof_species
                                  leftv(idxev + i)=y(i)
                                end do


                                else      !no real gas


                                oodensity=1.0d0/leftv(1)

                                temps(1)=leftv(1)


                                temps(2)=leftv(2)*oodensity; u=temps(2)
                                temps(3)=leftv(3)*oodensity; v=temps(3)
                                w=zero
                                if (dimensiona.eq.3)then
                                temps(4)=leftv(4)*oodensity; w=temps(4)
                                end if
                                skinx=(u**2 + v**2+ w**2)

                                temps(dimensiona+2)=((gamma-1.0d0))*((leftv(dimensiona+2))-oo2*leftv(1)*skinx)

                                leftv(1:nof_variables)=temps(1:nof_variables)
                                end if                    !

        end if              !multispecies if

end if


end subroutine cons2prim







subroutine cons2div(n,leftv,mp_pinfl,gammal)
implicit none
!> @brief
!> this subroutine transforms one vector of conservative variables to the variables required for the diffusion (r,u,v,w,t,y1,y2,y3..) or (r,u,v,w,ttr,tve,y1,y2,etc)
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::temps
real,dimension(1:nof_variables),intent(inout)::leftv
real,intent(inout)::mp_pinfl,gammal
real::oodensity,mp_density,mp_stiff,mp_temp
real::p_sat,p_tol, rho_g,rho_l, ss_g, ss_l, pp, p_gl, void_frac,p_temp,etr
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:gpu_max_species)::rg_vftemp
real,dimension(1:gpu_max_species)::rg_cvs
real,dimension(1:gpu_max_species)::rg_tvsl,rg_ev
real::rg_ve,rg_tr,rg_chem,rg_density,rg_ev_total,rg_rmix,rg_kin
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,tve,ttr,rg_f,rg_df,sum1,sum2,sum3,cv_i
real:: rhoe, rhoev, ke, etot, ev, u,v,w,cv_mix
real:: g, f, df, theta, ei, dei
real, dimension(1:gpu_max_species) :: rho_i,  cv_s
real :: echem, sumy
real:: rg_temp, rg_t_old,denom,tmpexp
integer :: i, iter
real, parameter :: rgtol = 1.0d-10
real, parameter :: rg_t_lo = 50.0d0, rg_t_hi = 20000.0d0
real:: rho
real:: rg_ttr2, rg_tv
real,dimension(1:gpu_max_species):: y







if (nof_variables.gt.1)then

if (dimensiona.eq.3)then

p_sat =2000
p_tol =10e-5

      if (multispecies.eq.1) then
            if (mp_modelc.eq.0)then
            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1

            do rg_i=1,nof_species-1
            sum2=sum2+leftv(5+nof_species+rg_i)
            mp_ar(rg_i)=leftv(5+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(5+nof_species+rg_i)
            end do

            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
              sum3=sum3+mp_ar(rg_i)
              end do
              gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
              oodensity=1.0d0/mp_density

            temps(1)=mp_density
            temps(2)=leftv(2)*oodensity
            temps(3)=leftv(3)*oodensity
            temps(4)=leftv(4)*oodensity
            sum3=zero
            do rg_i=1,nof_species
              sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
              end do

            mp_stiff=sum3*(gammal-1.0d0)
             sum2=zero
            do rg_i=1,nof_species
              sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
              end do

              mp_pinfl=sum2

            temps(5)=(((gammal-1.0d0))*((leftv(5))-oo2*temps(1)*(((temps(2))**2)+((temps(3))**2)+((temps(4))**2))))-mp_stiff

            temps(6:nof_variables)=leftv(6:nof_variables)



            call multispecies_temp(leftv,mp_temp)


            temps(5)=mp_temp


!             if(cavitation.eq.1)then
!             rho_g = leftv(6)/leftv(8)
!             rho_l = leftv(7)/leftv(8)
!             ss_g = sqrt(gamma_in(1)*(temps(5)+mp_pinf(1))/rho_g)
!             ss_l = sqrt(gamma_in(2)*(temps(5)+mp_pinf(2))/rho_l)
!
!             p_gl=rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_g-rho_l)/((rho_g*rho_g*ss_g*ss_g)-(rho_l*rho_l*ss_l*ss_l))
!             void_frac=(rho_g*ss_g*ss_g*rho_l*ss_l*ss_l*(rho_l+(leftv(8)*(rho_g-rho_l))))/(rho_l*((rho_g*ss_g*ss_g)-leftv(8)*((rho_g*ss_g*ss_g)-(rho_l*ss_l*ss_l))))
!
!             p_temp=temps(5)
!
!               if ((temps(5).gt.p_tol).and.(temps(5).lt.p_sat))then
!               p_temp=p_sat+p_gl*log(void_frac)
!               end if
!
!               if (temps(5).lt.p_tol)then
!               p_temp=p_tol
!               end if
!             temps(5)=p_temp
!
!             end if




      leftv(1:nof_variables)=temps(1:nof_variables)

            end if !(mp_modelc.eq.0)then

      else



      if (realgas.eq.1)then

              temps(:)=0.0d0

! !                call fix_conservative_state(leftv)
              !first get total density-correct

               !first get total density-correct !
               ! note:
               ! total density rho is a primary conserved variable from continuity.
               ! do not recompute rho as sum(rho_i) here.
               ! species equations are not perfectly conservative numerically and
               ! redefining rho severely degrades robustness in hypersonic flows.

              ! ============================================================
              ! 1. unpack conservative variables
              ! ============================================================

              rho  = leftv(1)
              u    = leftv(2) / rho
              v    = leftv(3) / rho
              if (dimensiona.eq.3)then
              w   =leftv(4) / rho
              else
              w=zero
              end if

              rhoe = leftv(dimensiona+2)
              rhoev = leftv(dimensiona+3)

              do i=1,nof_species
                rho_i(i) = leftv(dimensiona+3+i)
              end do

              ! ============================================================
              ! 2. recover mass fractions
              ! ============================================================
              sumy = 0.0d0
              do i=1,nof_species
                y(i) = rho_i(i) / rho
                sumy = sumy + y(i)
              end do

!               ! normalize (safety)
!               do i=1,nof_species
!                 y(i) = y(i) / sumy
!               end do

              ! ============================================================
              ! 3. compute kinetic energy and total specific energy
              ! ============================================================
              ke   = 0.5d0*(u*u + v*v + w*w)
              etot = rhoe / rho

              ! vibrational specific energy (already known from conservative vars)
              ev = rhoev / rho

              ! ============================================================
              ! 4. compute chemical energy
              ! ============================================================
              echem = 0.0d0
              do i=1,nof_species
                if (rg_hzero(i) > 0.0d0) then
                  ! hzero is j/mol → convert to j/kg
                  echem = echem - y(i) * (rg_hzero(i) / rg_molm(i))
                end if
              end do

              ! ============================================================
              ! 5. compute remaining translational energy
              !     etr = etot - ev - echem - ke
              ! ============================================================
              etr = etot - ev - echem - ke
              if (etr < 0.0d0) etr = 1d-12    ! safety

              ! ============================================================
              ! 6. compute ttr (no newton needed: etr = cv_mix * ttr)
              ! ============================================================
              cv_mix = 0.0d0
              do i = 1, nof_species
                if (i <= 3) then
                  cv_i = (5.0d0/2.0d0)*(rgs_ru/rg_molm(i))
                else
                  cv_i = (3.0d0/2.0d0)*(rgs_ru/rg_molm(i))
                end if
                cv_mix = cv_mix + y(i)*cv_i
              end do

              if (cv_mix > rgs_tiny) then
                rg_ttr2 = etr / cv_mix
              else
                rg_ttr2 = rg_t_lo
              end if

              rg_ttr2 = max(rg_t_lo, min(rg_t_hi, rg_ttr2))

              ! ============================================================
              ! 7. solve for tv using newton iteration
              !     ev = σ y_i * (r/m_i) * θ/(exp(θ/tv)-1)
              !     --> use ttr as initial guess
              ! ============================================================

              rg_tv = rg_ttr2   ! <-- use ttr as initial guess for tv
              rg_tv = max(rg_t_lo, min(rg_t_hi, rg_tv))

              do iter=1,40
                g  = 0.0d0
                df = 0.0d0

                do i=1,3   ! only n2,o2,no vibrate
                  if (y(i) < 1d-16) cycle
                  theta = rg_thetag(i)
                  if (theta <= 0.0d0) cycle

                  ei = theta / rg_tv

                  if (ei > 60.0d0) then
                    ! overflow-safe asymptotic form
                    tmpexp = exp(-ei)
                    g  = g  + y(i)*(rgs_ru/rg_molm(i))*theta*tmpexp
                    df = df + y(i)*(rgs_ru/rg_molm(i))*theta*(ei/rg_tv)*tmpexp
                  else
                    tmpexp = exp(ei)
                    denom  = tmpexp - 1.0d0
                    g  = g  + y(i)*(rgs_ru/rg_molm(i))*(theta/denom)
                    df = df + y(i)*(rgs_ru/rg_molm(i))*theta*ei*tmpexp / (denom*denom*rg_tv)
                  end if
                end do

                g = g - ev        ! target equation ev(tv) - ev_known = 0

                if (abs(g) < 1.0d-12 * max(1.0d0, abs(ev))) exit
                if (abs(df) < 1d-20) exit

                rg_tv = rg_tv - g/df

                ! enforce bounds
                if (rg_tv < rg_t_lo) rg_tv = rg_t_lo
                if (rg_tv > rg_t_hi) rg_tv = rg_t_hi
              end do



              leftv(1)=rho
              leftv(2)=u
              leftv(3)=v

              if (dimensiona.eq.3)then
              leftv(4)=w
              end if

              leftv(dimensiona+2)=rg_ttr2
              leftv(dimensiona+3)=rg_tv



              do i=1,nof_species
                 leftv(dimensiona+3+i)=y(i)
              end do



              else















      oodensity=1.0d0/leftv(1)

      temps(1)=leftv(1)
      temps(2)=leftv(2)*oodensity
      temps(3)=leftv(3)*oodensity
      temps(4)=leftv(4)*oodensity
      temps(5)=((gamma-1.0d0))*((leftv(5))-oo2*leftv(1)*(((temps(2))**2)+((temps(3))**2)+((temps(4))**2)))

      temps(5)=  leftv(5)/(leftv(1)*r_gas)  !temperature

      leftv(1:nof_variables)=temps(1:nof_variables)

      end if

end if




else

p_sat =2000
p_tol =10e-5


if ((multispecies.eq.1)) then
            if (mp_modelc.eq.0)then
            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1

            do rg_i=1,nof_species-1
            sum2=sum2+leftv(4+nof_species+rg_i)
            mp_ar(rg_i)=leftv(4+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(4+nof_species+rg_i)
            end do

             mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
              sum3=sum3+mp_ar(rg_i)
              end do
              gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
              oodensity=1.0d0/mp_density

              temps(1)=mp_density
              temps(2)=leftv(2)*oodensity
              temps(3)=leftv(3)*oodensity

            sum3=zero
            do rg_i=1,nof_species
              sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
              end do

            mp_stiff=sum3*(gammal-1.0d0)

             sum2=zero
            do rg_i=1,nof_species
              sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
              end do

              mp_pinfl=sum2

            temps(4)=(((gammal-1.0d0))*((leftv(4))-oo2*temps(1)*(((temps(2))**2)+((temps(3))**2))))-mp_stiff

            temps(5:nof_variables)=leftv(5:nof_variables)

            call multispecies_temp(leftv,mp_temp)


            temps(4)=mp_temp








 leftv(1:nof_variables)=temps(1:nof_variables)

           end  if !(mp_modelc.eq.0)then

else


if (realgas.eq.1)then

temps(:)=0.0d0

!            call fix_conservative_state(leftv)

              !first get total density-correct

!                do rg_i = 1, nof_species
!                                   if (leftv(dimensiona+3+rg_i).lt.0.0d0)then
!                                       leftv(dimensiona+3+rg_i)=0.0d0
!                                   end if
!                                 end do






!               do rg_i=1,nof_species
!               temps(1)=temps(1)+leftv(dimensiona+3+rg_i)
!               end do
!               rho  = temps(1)


              ! ============================================================
              ! 1. unpack conservative variables
              ! ============================================================

              rho  = leftv(1)
              u    = leftv(2) / rho
              v    = leftv(3) / rho
              if (dimensiona.eq.3)then
              w   =leftv(4) / rho
              else
              w=zero
              end if

              rhoe = leftv(dimensiona+2)
              rhoev = leftv(dimensiona+3)

              do i=1,nof_species
                rho_i(i) = leftv(dimensiona+3+i)
              end do

              ! ============================================================
              ! 2. recover mass fractions
              ! ============================================================
              sumy = 0.0d0
              do i=1,nof_species
                y(i) = rho_i(i) / rho
                sumy = sumy + y(i)
              end do

!               ! normalize (safety)
!               do i=1,nof_species
!                 y(i) = y(i) / sumy
!               end do

              ! ============================================================
              ! 3. compute kinetic energy and total specific energy
              ! ============================================================
              ke   = 0.5d0*(u*u + v*v + w*w)
              etot = rhoe / rho

              ! vibrational specific energy (already known from conservative vars)
              ev = rhoev / rho

              ! ============================================================
              ! 4. compute chemical energy
              ! ============================================================
              echem = 0.0d0
              do i=1,nof_species
                if (rg_hzero(i) > 0.0d0) then
                  ! hzero is j/mol → convert to j/kg
                  echem = echem - y(i) * (rg_hzero(i) / rg_molm(i))
                end if
              end do

              ! ============================================================
              ! 5. compute remaining translational energy
              !     etr = etot - ev - echem - ke
              ! ============================================================
              etr = etot - ev - echem - ke
              if (etr < 0.0d0) etr = 1d-12    ! safety

              ! ============================================================
              ! 6. compute ttr (no newton needed: etr = cv_mix * ttr)
              ! ============================================================
              cv_mix = 0.0d0
              do i = 1, nof_species
                if (i <= 3) then
                  cv_i = (5.0d0/2.0d0)*(rgs_ru/rg_molm(i))
                else
                  cv_i = (3.0d0/2.0d0)*(rgs_ru/rg_molm(i))
                end if
                cv_mix = cv_mix + y(i)*cv_i
              end do

              if (cv_mix > rgs_tiny) then
                rg_ttr2 = etr / cv_mix
              else
                rg_ttr2 = rg_t_lo
              end if

              rg_ttr2 = max(rg_t_lo, min(rg_t_hi, rg_ttr2))

              ! ============================================================
              ! 7. solve for tv using newton iteration
              !     ev = σ y_i * (r/m_i) * θ/(exp(θ/tv)-1)
              !     --> use ttr as initial guess
              ! ============================================================

              rg_tv = rg_ttr2   ! <-- use ttr as initial guess for tv
              rg_tv = max(rg_t_lo, min(rg_t_hi, rg_tv))

              do iter=1,40
                g  = 0.0d0
                df = 0.0d0

                do i=1,3   ! only n2,o2,no vibrate
                  if (y(i) < 1d-16) cycle
                  theta = rg_thetag(i)
                  if (theta <= 0.0d0) cycle

                  ei = theta / rg_tv

                  if (ei > 60.0d0) then
                    ! overflow-safe asymptotic form
                    tmpexp = exp(-ei)
                    g  = g  + y(i)*(rgs_ru/rg_molm(i))*theta*tmpexp
                    df = df + y(i)*(rgs_ru/rg_molm(i))*theta*(ei/rg_tv)*tmpexp
                  else
                    tmpexp = exp(ei)
                    denom  = tmpexp - 1.0d0
                    g  = g  + y(i)*(rgs_ru/rg_molm(i))*(theta/denom)
                    df = df + y(i)*(rgs_ru/rg_molm(i))*theta*ei*tmpexp / (denom*denom*rg_tv)
                  end if
                end do

                g = g - ev        ! target equation ev(tv) - ev_known = 0

                if (abs(g) < 1.0d-12 * max(1.0d0, abs(ev))) exit
                if (abs(df) < 1d-20) exit

                rg_tv = rg_tv - g/df

                ! enforce bounds
                if (rg_tv < rg_t_lo) rg_tv = rg_t_lo
                if (rg_tv > rg_t_hi) rg_tv = rg_t_hi
              end do



              leftv(1)=rho
              leftv(2)=u
              leftv(3)=v

              if (dimensiona.eq.3)then
              leftv(4)=w
              end if

              leftv(dimensiona+2)=rg_ttr2
              leftv(dimensiona+3)=rg_tv



              do i=1,nof_species
                 leftv(dimensiona+3+i)=y(i)
              end do









else















oodensity=1.0d0/leftv(1)

temps(1)=leftv(1)
temps(2)=leftv(2)*oodensity
temps(3)=leftv(3)*oodensity
temps(4)=((gamma-1.0d0))*((leftv(4))-oo2*leftv(1)*(((temps(2))**2)+((temps(3))**2)))

temps(4)=  leftv(4)/(leftv(1)*r_gas)  !temperature

leftv(1:nof_variables)=temps(1:nof_variables)



end if





end if




end if


end if


end subroutine cons2div














subroutine cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
implicit none
!> @brief
!> this subroutine transforms two vector of conservative variables to primitive variables
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv
real,intent(inout)::mp_pinfl,gammal
real,dimension(1:nof_variables),intent(inout)::rightv
real,intent(inout)::mp_pinfr,gammar

call cons2prim(n,leftv,mp_pinfl,gammal)
call cons2prim(n,rightv,mp_pinfr,gammar)


end subroutine cons2prim2


subroutine lmacht(n,leftv,rightv)
implicit none
!> @brief
!> this subroutine applies the low-mach number correction to two vectors of conserved variables
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables),intent(inout)::rightv
real::mp_pinfr,gammar
real::q2l,q2r,uul,uur,vvl,vvr,wwr,wwl,rhol,rhor,etal,etar,duu,dvv,dww
real::mach2,mach,cma,dus,dvs,dws,diff,c1o2,ssl,ssr,ppl,ppr,eel,eer,tole,mlm

tole=zero

 eel=leftv(5)
 eer=rightv(5)

 call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)


      c1o2=0.50d0
      rhol=leftv(1)
      uul=leftv(2)
      vvl=leftv(3)
      wwl=leftv(4)
       ppl=leftv(5)
     rhor=rightv(1)
      uur=rightv(2)
      vvr=rightv(3)
      wwr=rightv(4)
       ppr=rightv(5)
		
      q2l=(uul*uul)+(vvl*vvl)+(wwl*wwl)
      q2r=(uur*uur)+(vvr*vvr)+(wwr*wwr)

!	ssl=((gamma*ppl)/(rhol)); ssr=((gamma*ppr)/(rhor))
                                  
!	mlm=((abs(ppr-ppl))/(0.5*rres*((gamma*pres/rres))))




      if ((multispecies.eq.1).or.(realgas.eq.1))then
		 ssl=sqrt((leftv(5)+mp_pinfl)*gammal/leftv(1))
		  ssr=sqrt((rightv(5)+mp_pinfr)*gammar/rightv(1))
		else

      ssl=((gamma*ppl)/(rhol))
      ssr=((gamma*ppr)/(rhor))
      end if


      

      cma=1.0d0

      duu=uur-uul
      dvv=vvr-vvl
      dww=wwr-wwl

!       if(lmach.eq.1) then !standard proportional to du^2
         mach2=max(q2l/ssl,q2r/ssr)
         mach=sqrt(mach2)
         mach=min(cma*mach,1.0d0)
!       end if

      dus=uur+uul
      dvs=vvr+vvl
      dws=wwr+wwl

      duu=mach*duu
      dvv=mach*dvv
      dww=mach*dww


      diff=c1o2*duu

      uul=(dus*c1o2-diff)
      uur=(dus*c1o2+diff)

      if (lmach_style.eq.1)then
       diff=c1o2*dvv

      vvl=(dvs*c1o2-diff)
      vvr=(dvs*c1o2+diff)

      diff=c1o2*dww

      wwl=(dws*c1o2-diff)
      wwr=(dws*c1o2+diff)
      end if
	
    	
       leftv(1)=rhol
       leftv(2)=uul*rhol
       leftv(3)=vvl*rhol
       leftv(4)=wwl*rhol
       leftv(5)=eel
	
      rightv(1)=rhor 	
       rightv(2)=uur*rhor
       rightv(3)=vvr*rhor
       rightv(4)=wwr*rhor
       rightv(5)=eer




end subroutine lmacht



subroutine lmacht2d(n,leftv,rightv)
implicit none
!> @brief
!> this subroutine applies the low-mach number correction to two vectors of conserved variables 2d
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables),intent(inout)::rightv
real::mp_pinfr,gammar
real::q2l,q2r,uul,uur,vvl,vvr,wwr,wwl,rhol,rhor,etal,etar,duu,dvv,dww
real::mach2,mach,cma,dus,dvs,dws,diff,c1o2,ssl,ssr,ppl,ppr,eel,eer,tole,mlm

tole=tolsmall



 eel=leftv(4)
 eer=rightv(4)


 
 call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)



      c1o2=0.5d0
      rhol=leftv(1)
      uul=leftv(2)
      vvl=leftv(3)
       ppl=leftv(4)
       
     rhor=rightv(1)
      uur=rightv(2)
      vvr=rightv(3)
       ppr=rightv(4)
		
       q2l=(uul*uul)+(vvl*vvl)
       q2r=(uur*uur)+(vvr*vvr)
      
      
      
    if ((multispecies.eq.1).or.(realgas.eq.1))then
		 ssl=sqrt((leftv(4)+mp_pinfl)*gammal/leftv(1))
		  ssr=sqrt((rightv(4)+mp_pinfr)*gammar/rightv(1))
		else

      ssl=((gamma*ppl)/(rhol))
      ssr=((gamma*ppr)/(rhor))
      end if

      cma=1.0d0

      duu=uur-uul
      dvv=vvr-vvl
      
      

!       if(lmach.eq.1) then !standard proportional to du^2
         mach2=max(q2l/ssl,q2r/ssr)
         mach=sqrt(mach2)
         mach=min(cma*mach,1.0d0)
!       end if

      dus=uur+uul
      dvs=vvr+vvl
     

       !ul+zul+ur-zur=(ul+ur)-z(ur-ul))
      duu=mach*duu
      dvv=mach*dvv
      


      diff=c1o2*duu
	!if (uul*uur.gt.tole)then
	
      uul=(dus*c1o2)-diff
      uur=(dus*c1o2)+diff
      
if (lmach_style.eq.1)then
      diff=c1o2*dvv
      vvl=(dvs*c1o2-diff)
      vvr=(dvs*c1o2+diff)
end if
	
       leftv(1)=rhol
       leftv(2)=uul*rhol
       leftv(3)=vvl*rhol
       leftv(4)=eel
	
       rightv(1)=rhor 	
       rightv(2)=uur*rhor
       rightv(3)=vvr*rhor
       rightv(4)=eer




end subroutine lmacht2d































subroutine prim2cons(n,leftv)
implicit none
!> @brief
! !> this subroutine transforms one vector of primitive variables to conservative variables
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::temps
real::oodensity,skin1,ie1,mp_density,mp_stiff,rg_ve,rg_tr,rg_chem
real,dimension(1:nof_variables),intent(inout)::leftv
real::mp_pinfl,gammal,p
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:gpu_max_species)::rg_vftemp
real,dimension(1:gpu_max_species)::rg_cvs,rg_theta
real,dimension(1:gpu_max_species)::rg_tvsl,rg_ev,y
real::rg_density,rg_ev_total,rg_rmix,rg_kin
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,rg_tol,tve,ttr,rg_f,rg_df,sum1,sum2,sum3,u,v,w
real:: rmix
real:: cv_mix, e_tr, e_chem, ke, e_tot,evib,rho
integer :: i, idxe, idxev

rg_maxiter=20

if (nof_variables.gt.1)then

if (multispecies.eq.1) then

            if (mp_modelc.eq.0)then
            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1

            do rg_i=1,nof_species-1
            sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
            mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
            end do
            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)

             sum3=zero
            do rg_i=1,nof_species
              sum3=sum3+mp_ar(rg_i)
              end do
              gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
 

 
          temps(1)=mp_density
          temps(2)=leftv(2)*temps(1);u=leftv(2)
          temps(3)=leftv(3)*temps(1);v=leftv(3)
          w=zero
          if (dimensiona.eq.3)then
                    temps(4)=leftv(4)*temps(1);w=leftv(4)
                    end if
          skin1=(oo2)*((u*u)+(v*v)+(w*w))

            sum3=zero
            do rg_i=1,nof_species
              sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
              end do

            mp_stiff=sum3*(gammal-1.0d0)


! mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)

      ie1=((leftv(dimensiona+2)+mp_stiff)/((gammal-1.0d0)*temps(1)))
      temps(dimensiona+2)=temps(1)*(ie1+skin1)
      temps(dimensiona+3:nof_variables)=leftv(dimensiona+3:nof_variables)
      leftv(1:nof_variables)=temps(1:nof_variables)
 
            end if !(mp_modelc.eq.0)then
 
 else


            if (realgas.eq.1)then


            rho=leftv(1)
            u=leftv(2)
            v=leftv(3)

            if (dimensiona.eq.3)then
            w=leftv(4)
            else
            w=zero
            end if

            idxe  = dimensiona + 2
            idxev = dimensiona + 3

            p=leftv(idxe)
            evib=leftv(idxev)

                              do i = 1, nof_species
                                  y(i) = leftv(idxev + i)   ! species mass fractions
                                end do


              rmix = 0.0d0

             do i = 1, nof_species
              rmix = rmix + y(i) / rg_molm(i)
            end do




            ! -------------------------------
            ! 1. mixture gas constant rmix
            ! p = ρ rmix ttr → ttr = p / (ρ rmix)
            ! -------------------------------
            rmix = 0.0d0
            do i = 1, nof_species
              rmix = rmix + y(i) / rg_molm(i)
            end do
            rmix = rgs_ru * rmix

            ! -------------------------------
            ! 2. mixture cv and translational energy
            !    cv_i = 5/2 r for i<=3 (diatomic),
            !           3/2 r for i>3  (monatomic)
            ! -------------------------------
            cv_mix = 0.0d0
            do i = 1, nof_species
              if (i <= 3) then
                cv_mix = cv_mix + y(i) * (5.0d0/2.0d0) * (rgs_ru / rg_molm(i))
              else
                cv_mix = cv_mix + y(i) * (3.0d0/2.0d0) * (rgs_ru / rg_molm(i))
              end if
            end do

            ! translational temperature ttr from p = ρ rmix ttr
            ! then e_tr = cv_mix * ttr
            if (rmix > 1.0d-20 .and. cv_mix > 1.0d-20) then
              e_tr = cv_mix * (p / (rho * rmix))
            else
              e_tr = 0.0d0
            end if

            ! -------------------------------
            ! 3. chemical energy: e_chem = - σ y_i * h°_i/m_i
            ! -------------------------------
            e_chem = 0.0d0
!             do i = 1, nof_species
!               if (rg_hzero(i) > 0.0d0) then
!                 e_chem = e_chem - y(i) * (rg_hzero(i) / rg_molm(i))
!               end if
!             end do

            ! -------------------------------
            ! 4. kinetic and total specific energy
            ! -------------------------------
            ke    = 0.5d0 * (u*u + v*v + w*w)
            e_tot = e_tr + evib + e_chem + ke

            ! -------------------------------
            ! 5. fill conservative vector
            ! -------------------------------
            leftv(1) = rho
            leftv(2) = rho * u
            leftv(3) = rho * v
            if (dimensiona == 3) then
              leftv(4) = rho * w
            else
              leftv(4) = 0.0d0
            end if



            leftv(idxe)  = rho * e_tot      ! ρe_total
            leftv(idxev) = rho * evib       ! ρ e_vib

            do i = 1, nof_species
              leftv(idxev + i) = rho * y(i) ! ρ_i
            end do




            else



                      u=leftv(2)
                      v=leftv(3)
                      w=zero
                            if (dimensiona.eq.3)then
                                w=leftv(4)
                                end if
                      skin1=(oo2)*((u*u)+(v*v)+(w*w))
            ie1=((leftv(dimensiona+2))/((gamma-1.0d0)*leftv(1)))

            oodensity=1.0d0/leftv(1)

            temps(1)=leftv(1)
            temps(2)=leftv(2)*leftv(1)
            temps(3)=leftv(3)*leftv(1)
            if (dimensiona.eq.3)then
            temps(4)=leftv(4)*leftv(1)
            end if
            temps(dimensiona+2)=leftv(1)*(ie1+skin1)

            leftv(1:nof_variables)=temps(1:nof_variables)

            end if

end if

end if






end subroutine prim2cons







subroutine prim2cons2(n,leftv,rightv)
implicit none
!> @brief
!> this subroutine transforms two vectors of primitive variables to conservative variables
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:gpu_max_nvar)::temps
real::oodensity,skin1,ie1,mp_density,mp_stiff
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real::mp_pinfl,gammal,mp_pinfr,gammar
integer::rg_i,rg_j
real,dimension(1:gpu_max_species)::rg_vftemp
real,dimension(1:gpu_max_species)::rg_cvs,rg_theta
real,dimension(1:gpu_max_species)::rg_tvsl,rg_ev
real::rg_density,rg_ev_total,rg_rmix,rg_kin
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,rg_tol,tve,ttr,rg_f,rg_df,sum1,sum2,sum3

if (nof_variables.gt.1)then

call prim2cons(n,leftv)
call prim2cons(n,rightv)


end if


end subroutine prim2cons2


subroutine cons2prim_ideal(n,leftv,mp_pinfl,gammal)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv
real,intent(inout)::mp_pinfl,gammal
real::oodensity,u,v,w,skinx

mp_pinfl=zero
gammal=gamma
if (nof_variables.gt.1)then
  oodensity=1.0d0/leftv(1)
  u=leftv(2)*oodensity
  v=leftv(3)*oodensity
  w=zero
  if (dimensiona.eq.3) w=leftv(4)*oodensity
  skinx=(u*u)+(v*v)+(w*w)
  leftv(2)=u
  leftv(3)=v
  if (dimensiona.eq.3) leftv(4)=w
  leftv(dimensiona+2)=((gamma-1.0d0))*((leftv(dimensiona+2))-oo2*leftv(1)*skinx)
end if

end subroutine cons2prim_ideal


subroutine cons2div_ideal(n,leftv,mp_pinfl,gammal)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv
real,intent(inout)::mp_pinfl,gammal
real::oodensity

mp_pinfl=zero
gammal=gamma
if (nof_variables.gt.1)then
  oodensity=1.0d0/leftv(1)
  leftv(2)=leftv(2)*oodensity
  leftv(3)=leftv(3)*oodensity
  if (dimensiona.eq.3)then
    leftv(4)=leftv(4)*oodensity
    leftv(5)=leftv(5)/(leftv(1)*r_gas)
  else
    leftv(4)=leftv(4)/(leftv(1)*r_gas)
  end if
end if

end subroutine cons2div_ideal


subroutine cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,intent(inout)::mp_pinfl,mp_pinfr,gammal,gammar

call cons2prim_ideal(n,leftv,mp_pinfl,gammal)
call cons2prim_ideal(n,rightv,mp_pinfr,gammar)

end subroutine cons2prim2_ideal


subroutine prim2cons_ideal(n,leftv)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv
real::u,v,w,skin1,ie1

if (nof_variables.gt.1)then
  u=leftv(2)
  v=leftv(3)
  w=zero
  if (dimensiona.eq.3) w=leftv(4)
  skin1=oo2*((u*u)+(v*v)+(w*w))
  ie1=leftv(dimensiona+2)/((gamma-1.0d0)*leftv(1))
  leftv(2)=leftv(1)*u
  leftv(3)=leftv(1)*v
  if (dimensiona.eq.3) leftv(4)=leftv(1)*w
  leftv(dimensiona+2)=leftv(1)*(ie1+skin1)
end if

end subroutine prim2cons_ideal


subroutine prim2cons2_ideal(n,leftv,rightv)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv,rightv

call prim2cons_ideal(n,leftv)
call prim2cons_ideal(n,rightv)

end subroutine prim2cons2_ideal


subroutine lmacht_ideal(n,leftv,rightv)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real::mp_pinfl,mp_pinfr,gammal,gammar
real::q2l,q2r,uul,uur,vvl,vvr,wwl,wwr,rhol,rhor,duu,dvv,dww
real::mach2,mach,cma,dus,dvs,dws,diff,c1o2,ssl,ssr,ppl,ppr,eel,eer

eel=leftv(5)
eer=rightv(5)
call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
c1o2=0.5d0
rhol=leftv(1); uul=leftv(2); vvl=leftv(3); wwl=leftv(4); ppl=leftv(5)
rhor=rightv(1); uur=rightv(2); vvr=rightv(3); wwr=rightv(4); ppr=rightv(5)
q2l=(uul*uul)+(vvl*vvl)+(wwl*wwl)
q2r=(uur*uur)+(vvr*vvr)+(wwr*wwr)
ssl=(gamma*ppl)/rhol
ssr=(gamma*ppr)/rhor
cma=1.0d0
duu=uur-uul
dvv=vvr-vvl
dww=wwr-wwl
mach2=max(q2l/ssl,q2r/ssr)
mach=sqrt(mach2)
mach=min(cma*mach,1.0d0)
dus=uur+uul
dvs=vvr+vvl
dws=wwr+wwl
duu=mach*duu
dvv=mach*dvv
dww=mach*dww
diff=c1o2*duu
uul=dus*c1o2-diff
uur=dus*c1o2+diff
if (lmach_style.eq.1)then
  diff=c1o2*dvv
  vvl=dvs*c1o2-diff
  vvr=dvs*c1o2+diff
  diff=c1o2*dww
  wwl=dws*c1o2-diff
  wwr=dws*c1o2+diff
end if
leftv(1)=rhol
leftv(2)=uul*rhol
leftv(3)=vvl*rhol
leftv(4)=wwl*rhol
leftv(5)=eel
rightv(1)=rhor
rightv(2)=uur*rhor
rightv(3)=vvr*rhor
rightv(4)=wwr*rhor
rightv(5)=eer

end subroutine lmacht_ideal


subroutine lmacht2d_ideal(n,leftv,rightv)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real::mp_pinfl,mp_pinfr,gammal,gammar
real::q2l,q2r,uul,uur,vvl,vvr,rhol,rhor,duu,dvv
real::mach2,mach,cma,dus,dvs,diff,c1o2,ssl,ssr,ppl,ppr,eel,eer

eel=leftv(4)
eer=rightv(4)
call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
c1o2=0.5d0
rhol=leftv(1); uul=leftv(2); vvl=leftv(3); ppl=leftv(4)
rhor=rightv(1); uur=rightv(2); vvr=rightv(3); ppr=rightv(4)
q2l=(uul*uul)+(vvl*vvl)
q2r=(uur*uur)+(vvr*vvr)
ssl=(gamma*ppl)/rhol
ssr=(gamma*ppr)/rhor
cma=1.0d0
duu=uur-uul
dvv=vvr-vvl
mach2=max(q2l/ssl,q2r/ssr)
mach=sqrt(mach2)
mach=min(cma*mach,1.0d0)
dus=uur+uul
dvs=vvr+vvl
duu=mach*duu
dvv=mach*dvv
diff=c1o2*duu
uul=dus*c1o2-diff
uur=dus*c1o2+diff
if (lmach_style.eq.1)then
  diff=c1o2*dvv
  vvl=dvs*c1o2-diff
  vvr=dvs*c1o2+diff
end if
leftv(1)=rhol
leftv(2)=uul*rhol
leftv(3)=vvl*rhol
leftv(4)=eel
rightv(1)=rhor
rightv(2)=uur*rhor
rightv(3)=vvr*rhor
rightv(4)=eer

end subroutine lmacht2d_ideal












subroutine inflow(initcond,pox,poy,poz,inflow_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inflow in 3d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::inflow_out
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real::p,u,v,w,e,r,s,gm,skin,ien,pi
real::xf,yf,zf
real:: theta_0,vtang, vradial,gammar
real::mp_density,mp_stiff,sum1,sum2,sum3
real,dimension(1:gpu_max_nvar)::vect_in
real,dimension(gpu_max_species)::mp_ar,mp_ie
integer::rg_i,rg_j
real::khx,t1l,vhx,amp,dvel,rgg,tt1,khi_slope,khi_b,theeta,reeta,rg_ve,rg_tr,rg_chem,rg_density,rg_ev_total,rg_rmix
real,dimension(1:gpu_max_species)::rg_cvs
real,dimension(1:gpu_max_species)::rg_tvsl,rg_ev
real::rg_tv,rg_ttr0,rg_tve0

if (multispecies.eq.1) then



p=pres
u=uvel
v=vvel
w=wvel


! vect_in(1)=r
vect_in(2)=u
vect_in(3)=v
vect_in(4)=w
vect_in(5)=p
sum3=0.0d0
do rg_i=1,nof_species
vect_in(5+rg_i)=mp_r_in(rg_i)*mp_a_in(rg_i)
sum3=sum3+mp_r_in(rg_i)*mp_a_in(rg_i)
end do


vect_in(1)=sum3


do rg_i=1,nof_species-1
vect_in(5+nof_species+rg_i)=mp_a_in(rg_i)
end do

   call prim2cons(n,vect_in)


do rg_i=1,nof_variables
inflow_out(rg_i)=vect_in(rg_i)
end do


!
!
!
!
!
!
! mp_ar(1)=mp_a_in(1)/(gamma_in(1)-1.0d0)
! mp_ar(2)=mp_a_in(2)/(gamma_in(2)-1.0d0)
! gammar=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumption
!
! gm=gammar
!
! r=(mp_r_in(1)*mp_a_in(1))+(mp_r_in(2)*mp_a_in(2))
! mp_ie(1)=((p+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.0d0)))
! mp_ie(2)=((p+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.0d0)))
!
! ien=(mp_ie(1)*mp_a_in(1))+(mp_ie(2)*mp_a_in(2))
! ! !kinetic energy first!
! skin=(oo2)*((u**2)+(v**2)+(w**2))
! ! !total energy
! e=(r*skin)+ien
!
! !vector of conserved variables now
! inflow(1)=r
! inflow(2)=r*u
! inflow(3)=r*v
! inflow(4)=r*w
! inflow(5)=e
! inflow(6)=mp_r_in(1)*mp_a_in(1)
! inflow(7)=mp_r_in(2)*mp_a_in(2)
! inflow(8)=mp_a_in(1)



else

r=rres
gm=gamma
p=pres
u=uvel
v=vvel
w=wvel



if (initcond.eq.10000)then
if (sqrt(((poy(1)-0.0)**2)+((poz(1)-0.5)**2)).le.0.05)then
!if (((poy(1).ge.-0.05).and.(poy(1).le.0.05)).and.((poz(1).ge.-0.05).and.(poz(1).le.0.05)))then
p=0.4127
	r=5
	u=30.0
	v=0.0d0
	w=0.0d0

end if
end if







!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
!vector of conserved variables now
inflow_out(1)=r
inflow_out(2)=r*u
inflow_out(3)=r*v
inflow_out(4)=r*w
inflow_out(5)=e

end if


if (swirl.eq.1)then

if (pox(1).lt.-0.03)then
xf=pox(1)
yf=poy(1)
zf=poz(1)
theta_0=atan2(zf,yf)
vtang=18.0375d0
vradial=-12.63d0
u=0.0d0

v=-vtang*sin(theta_0)+vradial*cos(theta_0)
w=vtang*cos(theta_0)+vradial*sin(theta_0)



else

 u=70.06d0
  v=0.0d0
  w=0.0d0


end if

r=rres
gm=gamma
p=pres
s=sqrt((gm*p)/(r))
  
  
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
!vector of conserved variables now
inflow_out(1)=r
inflow_out(2)=r*u
inflow_out(3)=r*v
inflow_out(4)=r*w
inflow_out(5)=e  
  
  
end if


if (realgas.eq.1)then
u=uvel
v=vvel
w=wvel
p=pres
r=rres



!first build mixture gas constant
rg_rmix=zero
do rg_i=1,nof_species
  rg_rmix=rg_rmix+rg_vf(rg_i)/rg_molm(rg_i)
end do

  rg_rmix=rgs_ru*rg_rmix

  rg_ttr0=rg_ttr
  rg_tve0=rg_tve


! translational-rotational internal energy
rg_tr = 0.0d0
    do rg_i = 1, nof_species
      if (rg_i <= 3) then
        rg_cvs(rg_i) = (5.0d0 / 2.0d0) * rgs_ru / rg_molm(rg_i)
      else
        rg_cvs(rg_i) = (3.0d0 / 2.0d0) * rgs_ru / rg_molm(rg_i)
      end if
      rg_tr = rg_tr + rg_vf(rg_i) * rg_cvs(rg_i) * rg_ttr0
    end do


! vibrational energy
    rg_ev_total= 0.0d0
    do rg_i = 1, 3
      rg_ev_total = rg_ev_total + rg_vf(rg_i)  * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/rg_tve0) - 1.0d0))
    end do

rg_chem=zero

  ! chemical energy
 do rg_i=1,nof_species
        if (rg_hzero(rg_i).gt.1.0e-12)then
        rg_chem=rg_chem-(rg_vf(rg_i)*rg_hzero(rg_i)/rg_molm(rg_i))
        end if
end do
!rg_chem=zero





! kinetic energy

skin=(oo2)*((u**2)+(v**2)+(w**2))


inflow_out(1)=r
inflow_out(2)=r*u
inflow_out(3)=r*v
inflow_out(4)=r*w
inflow_out(5)=r*(rg_ev_total+rg_tr+rg_chem+skin)
inflow_out(6)=r*rg_ev_total

do rg_i=1,nof_species
inflow_out(6+rg_i)=r * rg_vf(rg_i)
end do



end if






end subroutine inflow



subroutine vect_function(pox,poy,vect_out)
implicit none
!> @brief
!> this makes a multipliciation between two vectors
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,dimension(3),intent(out)::vect_out
real,dimension(1:dimensiona),intent(in)::pox,poy

vect_out(1)=(poy(2)*pox(3))-(poy(3)*pox(2))
vect_out(2)=(poy(3)*pox(1))-(poy(1)*pox(3))
vect_out(3)=(poy(1)*pox(2))-(poy(2)*pox(1))





end subroutine vect_function

subroutine inflow2d(initcond,pox,poy,inflow2d_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inflow in 2d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::inflow2d_out
integer,intent(in)::initcond
real,dimension(1:2),intent(in)::pox,poy
real::p,u,v,w,e,r,s,gm,skin,ien,pi,ps
real::xf,yf,zf,lit_a,lit_o
real:: theta_0,vtang, vradial,gammar
real::mp_density,mp_stiff,sum1,sum2,sum3
real,dimension(1:gpu_max_nvar)::vect_in
real,dimension(gpu_max_species)::mp_ar,mp_ie
integer::rg_i,rg_j
real::khx,t1l,vhx,amp,dvel,rgg,tt1,khi_slope,khi_b,theeta,reeta,rg_ve,rg_tr,rg_chem,rg_density,rg_ev_total,rg_rmix
real,dimension(1:gpu_max_species)::rg_cvs
real,dimension(1:gpu_max_species)::rg_tvsl,rg_ev
real::rg_tv,rg_ttr0,rg_tve0





if (multispecies.eq.1) then


r=rres
p=pres
u=uvel
v=vvel

!time variable boundary condition
if (initcond.eq.430)then
v=(179299.375638680*(t**5)) - (82455.0868677361*(t**4)) + (14299.8472891299*(t**3)) - (1281.65548492021*(t**2)) + (62.2260666329356*t) - (0.0419033282181554)
end if


! if (initcond.eq.157)then
! ps=35*10e6
! lit_a=1.48*10e8
! lit_o=1.21*10e8
! p=pres+2.0d0*ps*exp(-lit_a*t)*cos((lit_o*t)+(pi/3.0))
! uvel=0.0
! vvel=0.0
! end if



! vect_in(1)=r
vect_in(2)=u
vect_in(3)=v
vect_in(4)=p
sum3=0.0d0
do rg_i=1,nof_species
vect_in(4+rg_i)=mp_r_in(rg_i)*mp_a_in(rg_i)
sum3=sum3+mp_r_in(rg_i)*mp_a_in(rg_i)
end do
vect_in(1)=sum3
do rg_i=1,nof_species-1
vect_in(4+nof_species+rg_i)=mp_a_in(rg_i)
end do

   call prim2cons(n,vect_in)


do rg_i=1,nof_variables
inflow2d_out(rg_i)=vect_in(rg_i)
end do

! !vector of conserved variables now
! inflow2d(1)=r
! inflow2d(2)=r*u
! inflow2d(3)=r*v
! inflow2d(4)=e
! inflow2d(5)=mp_r_in(1)*mp_a_in(1)
! inflow2d(6)=mp_r_in(2)*mp_a_in(2)
! inflow2d(7)=mp_a_in(1)



else




r=rres
gm=gamma
p=pres
u=uvel
v=vvel

if (initcond.eq.133)then
p=195557.25
	r=p/(350.5d0*287.058d0)
	u=168.62
	v=0.0d0

end if

if (initcond.eq.10000)then
if ((poy(1).ge.-0.05).and.(poy(1).le.0.05))then
p=0.4127
	r=5
	u=30.0
	v=0.0d0

end if
end if




if (initcond.eq.790)then
r=(2.4d0*6**2)/((0.4*6**2)+2)
u=(6*sqrt(1.4))*(70/(2.4*36))
v=0.0d0
p=(2.8*36-0.4)/(2.4)


end if




!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
!vector of conserved variables now
inflow2d_out(1)=r
inflow2d_out(2)=r*u
inflow2d_out(3)=r*v
inflow2d_out(4)=e


endif

if (realgas.eq.1)then
u=uvel
v=vvel
p=pres
r=rres



!first build mixture gas constant
rg_rmix=zero
do rg_i=1,nof_species
  rg_rmix=rg_rmix+rg_vf(rg_i)/rg_molm(rg_i)
end do

  rg_rmix=rgs_ru*rg_rmix

  rg_ttr0=rg_ttr
  rg_tve0=rg_tve


! translational-rotational internal energy
rg_tr = 0.0d0
    do rg_i = 1, nof_species
      if (rg_i <= 3) then
        rg_cvs(rg_i) = (5.0d0 / 2.0d0) * rgs_ru / rg_molm(rg_i)
      else
        rg_cvs(rg_i) = (3.0d0 / 2.0d0) * rgs_ru / rg_molm(rg_i)
      end if
      rg_tr = rg_tr + rg_vf(rg_i) * rg_cvs(rg_i) * rg_ttr0
    end do


! vibrational energy
    rg_ev_total= 0.0d0
    do rg_i = 1, 3
      rg_ev_total = rg_ev_total + rg_vf(rg_i)  * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/rg_tve0) - 1.0d0))
    end do

rg_chem=zero

  ! chemical energy
 do rg_i=1,nof_species
        if (rg_hzero(rg_i).gt.1.0e-12)then
        rg_chem=rg_chem-(rg_vf(rg_i)*rg_hzero(rg_i)/rg_molm(rg_i))
        end if
end do
!rg_chem=zero





! kinetic energy

skin=(oo2)*((u**2)+(v**2))


inflow2d_out(1)=r
inflow2d_out(2)=r*u
inflow2d_out(3)=r*v
inflow2d_out(4)=r*(rg_ev_total+rg_tr+rg_chem+skin)
inflow2d_out(5)=r*rg_ev_total

do rg_i=1,nof_species
inflow2d_out(5+rg_i)=r * rg_vf(rg_i)
end do



end if






end subroutine inflow2d

subroutine outflow2d(initcond,pox,poy,outflow2d_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 2d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::outflow2d_out
integer,intent(in)::initcond
real,dimension(1:2),intent(in)::pox,poy
real::p,u,v,w,e,r,s,gm,skin,ien,pi
real::xf,yf,zf,gammar
real,dimension(gpu_max_species)::mp_ar,mp_ie


if (multispecies.eq.1) then



p=pres
u=uvel
v=vvel
mp_ar(1)=mp_a_in(1)/(gamma_in(1)-1.0d0)  
mp_ar(2)=mp_a_in(2)/(gamma_in(2)-1.0d0)
gammar=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumption

gm=gammar

r=(mp_r_in(1)*mp_a_in(1))+(mp_r_in(2)*mp_a_in(2))
mp_ie(1)=((p+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.0d0)))
mp_ie(2)=((p+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.0d0)))
ien=(mp_ie(1)*mp_a_in(1))+(mp_ie(2)*mp_a_in(2))
! !kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
! !total energy
e=(r*skin)+ien

!vector of conserved variables now
outflow2d_out(1)=r
outflow2d_out(2)=r*u
outflow2d_out(3)=r*v
outflow2d_out(4)=e
outflow2d_out(5)=mp_r_in(1)*mp_a_in(1)
outflow2d_out(6)=mp_r_in(2)*mp_a_in(2)
outflow2d_out(7)=mp_a_in(1)

else

r=rres
gm=gamma
p=pres
u=uvel
v=vvel


!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
!vector of conserved variables now
outflow2d_out(1)=r
outflow2d_out(2)=r*u
outflow2d_out(3)=r*v
outflow2d_out(4)=e
















end if

end subroutine outflow2d
subroutine outflow(initcond,pox,poy,poz,outflow_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::outflow_out
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real::p,u,v,w,e,r,s,gm,skin,ien,pi
real::xf,yf,zf,gammar
real,dimension(gpu_max_species)::mp_ar,mp_ie

if (multispecies.eq.1) then



p=pres
u=uvel
v=vvel
w=wvel
mp_ar(1)=mp_a_in(1)/(gamma_in(1)-1.0d0)  
mp_ar(2)=mp_a_in(2)/(gamma_in(2)-1.0d0)
gammar=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumption

gm=gammar

r=(mp_r_in(1)*mp_a_in(1))+(mp_r_in(2)*mp_a_in(2))
mp_ie(1)=((p+(gamma_in(1)*mp_pinf(1)))/((gamma_in(1)-1.0d0)))
mp_ie(2)=((p+(gamma_in(2)*mp_pinf(2)))/((gamma_in(2)-1.0d0)))
ien=(mp_ie(1)*mp_a_in(1))+(mp_ie(2)*mp_a_in(2))
! !kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
! !total energy
e=(r*skin)+ien

!vector of conserved variables now
outflow_out(1)=r
outflow_out(2)=r*u
outflow_out(3)=r*v
outflow_out(4)=r*w
outflow_out(5)=e
outflow_out(6)=mp_r_in(1)*mp_a_in(1)
outflow_out(7)=mp_r_in(2)*mp_a_in(2)
outflow_out(8)=mp_a_in(1)

else

r=rres
gm=gamma
p=pres

if (initcond.eq.977)then
 p=101325
 end if 

u=uvel
v=vvel
w=wvel

!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
!vector of conserved variables now
outflow_out(1)=r
outflow_out(2)=r*u
outflow_out(3)=r*v
outflow_out(4)=r*w
outflow_out(5)=e

end if

end subroutine outflow

subroutine outflow2(initcond,pox,poy,poz,outflow2_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::outflow2_out
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real::p,u,v,w,e,r,s,gm,skin,ien,pi
real::xf,yf,zf
real,dimension(gpu_max_species)::mp_ar,mp_ie




r=rres
gm=gamma
p=press_outlet
u=uvel
v=vvel
w=wvel






!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
!vector of conserved variables now
outflow2_out(1)=r
outflow2_out(2)=r*u
outflow2_out(3)=r*v
outflow2_out(4)=r*w
outflow2_out(5)=e



end subroutine outflow2


subroutine bleed2d(iconsidered,facex,pox,poy,bleed2d_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::bleed2d_out
integer,intent(in)::iconsidered, facex
real,dimension(1:dimensiona),intent(in)::pox,poy
integer::ibleedn,rg_i
real::p,u,v,w,e,r,s,gm,skin,ien,pi
real::xf,yf,zf
real::cell_area
!------coefficients for bleed-----!
real::bleed_qsonic_s,bleed_mdotsonic_s,bleed_area,bleed_mdotsonic,bleed_region
real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2
real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot






cell_area=ielem_surf(facex,iconsidered)




!leftv is the input in conservative variables that we transform to primitive variables
call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
r=leftv(1)
u=leftv(2)
v=leftv(3)
p=leftv(4)







!now find the bleed number conditions for this face
ibleedn=ielem_bleedn(facex,iconsidered)



bleed_region=sqrt(((bleed_start(ibleedn,1)-bleed_end(ibleedn,1))**2)+((bleed_start(ibleedn,2)-bleed_end(ibleedn,2))**2))

!now compute the bleed area

bleed_area=bleed_porosity(ibleedn)*bleed_region

!now compute the bleed m dot sonic-s mass flow rate eq.17 , https://doi.org/10.2514/1.b37474

! bleed_mdotsonic_s=bleed_area*p*(sqrt(gamma*r/p))*(((gamma+1.0d0)/(2))**((gamma+1)/(2*(1-gamma))))

!now compute qsonic eq. 22
bleed_qsonic_s=0.598+0.0307*(bleed_plenum(ibleedn)/p)-0.5936*((bleed_plenum(ibleedn)/p)**2)

!equation 16
bleed_mdotsonic=bleed_qsonic_s*bleed_mdotsonic_s

!now that i have computed the bleed mass flow rate i have to distribute to the surface area of this tagged cell in the normal direction




call prim2cons2(n,leftv,rightv)

rightv(1:nof_variables)=leftv(1:nof_variables)


call rotatef2d(n,cright_rot,rightv,angle1,angle2)



cright_rot(2)=bleed_mdotsonic/cell_area



call rotateb2d(n,rightv,cright_rot,angle1,angle2)




do rg_i=1,nof_variables
bleed2d_out(rg_i)=rightv(rg_i)
end do


end subroutine bleed2d

subroutine bleed3d(iconsidered,facex,pox,poy,poz,bleed3d_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(out)::bleed3d_out
integer,intent(in)::iconsidered, facex
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
integer::ibleedn,rg_i
real::p,u,v,w,e,r,s,gm,skin,ien,pi
real::xf,yf,zf
real::cell_area
!------coefficients for bleed-----!
real::bleed_qsonic_s,bleed_mdotsonic_s,bleed_area,bleed_mdotsonic,bleed_region
real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2
real,dimension(1:gpu_max_nvar_total)::cleft,cright,cleft_rot,cright_rot


do rg_i=1,nof_variables
bleed3d_out(rg_i)=rightv(rg_i)
end do


end subroutine bleed3d



subroutine pass_inlet(initcond,pox,poy,poz,pass_inlet_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inlet for a passive scalar
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:passivescalar),intent(out)::pass_inlet_out
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
integer::rg_i

do rg_i=1,passivescalar
pass_inlet_out(rg_i)=1.0d0*rres
end do

end subroutine pass_inlet

subroutine pass_inlet2d(initcond,pox,poy,pass_inlet2d_out)
implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inlet for a passive scalar in 2d
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:passivescalar),intent(out)::pass_inlet2d_out
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy
integer::rg_i

do rg_i=1,passivescalar
pass_inlet2d_out(rg_i)=1.0d0
end do

end subroutine pass_inlet2d








subroutine shear_x(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml

 

 
vortet1(1:3,1:3) = rec_grads(1:3,1:3,iconsidered)


 ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
 vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
 wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=(cos(angle1)*sin(angle2))
 ny=(sin(angle1)*sin(angle2))
 nz=(cos(angle2))
 leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
tauyx=(uy + vx)
tauzx=(wx + uz)
tauzy=(vz + wy)
ssx=(viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))
ssy=(viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))
ssz=(viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))




if(rframe.eq.0)then
shear_temp=-ssx/(0.5*rres*ufreestream*ufreestream)
else
shear_temp=-ssx/(0.5*rres*v_ref*v_ref)
end if

end subroutine shear_x




subroutine shear_y(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in y-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 
vortet1(1:3,1:3) = rec_grads(1:3,1:3,iconsidered)


 ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
 vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
 wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=(cos(angle1)*sin(angle2))
 ny=(sin(angle1)*sin(angle2))
 nz=(cos(angle2))
 leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
tauyx=(uy + vx)
tauzx=(wx + uz)
tauzy=(vz + wy)
ssx=(viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))
ssy=(viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))
ssz=(viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))


if(rframe.eq.0)then
shear_temp=-ssy/(0.5*rres*ufreestream*ufreestream)
else
shear_temp=-ssy/(0.5*rres*v_ref*v_ref)
end if

end subroutine shear_y


subroutine shear_z(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in z-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 
vortet1(1:3,1:3) = rec_grads(1:3,1:3,iconsidered)


 ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
 vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
 wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=(cos(angle1)*sin(angle2))
 ny=(sin(angle1)*sin(angle2))
 nz=(cos(angle2))
 leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
tauyx=(uy + vx)
tauzx=(wx + uz)
tauzy=(vz + wy)
ssx=(viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))
ssy=(viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))
ssz=(viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))


if(rframe.eq.0)then
shear_temp=-ssz/(0.5*rres*ufreestream*ufreestream)
else
shear_temp=-ssz/(0.5*rres*v_ref*v_ref)
end if

end subroutine shear_z







subroutine heat_x(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp,lam_qflux
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 integer::i,k,j,kmaxe,gqi_points,nnd,im
 real,dimension(1:gpu_max_nvar)::leftv
 real,dimension(1:3)::temp_grad
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz,surface_temp
real,dimension(1:4)::viscl,laml
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:gpu_max_dim,1:gpu_max_qp_face)::qpoints2d
real,dimension(1:gpu_max_qp_face)::wequa2d


i=iconsidered
ssx=zero;ssy=zero;ssz=zero
j=facex
			      angle1=ielem_faceanglex(j,i)
			      angle2=ielem_faceangley(j,i)
			      nx=(cos(angle1)*sin(angle2))
			      ny=(sin(angle1)*sin(angle2))
			      nz=(cos(angle2))

                select case(ielem_types_faces(j,i))
				case (5)
					  gqi_points=qp_quad_n



					  if(reduce_comp.eq.1)then
					  wequa2d=1.0d0;
					  else
					    nnd=4
				      do k=1,nnd
					vext(k,1:dims)=inoder4_cord(1:dims,ielem_nodes_faces(j,k,i))
				      end do
					  call  quadraturequad3d(n,igqrules,vext,qpoints2d,wequa2d)
					  end if
					  surface_temp=ielem_surf(j,i)


				case(6)
					gqi_points=qp_triangle_n


					if(reduce_comp.eq.1)then
					  wequa2d=1.0d0;
					  else
					  nnd=3
					do k=1,nnd
					  vext(k,1:dims)=inoder4_cord(1:dims,ielem_nodes_faces(j,k,i))
					end do
					call quadraturetriang(n,igqrules,vext,qpoints2d,wequa2d)
					end if
 					    surface_temp=ielem_surf(j,i)



				end select




				do im=1,gqi_points
				temp_grad(1:3)=rec_uleftv(1:3,dimensiona+1,j,im,i)
				if (dg.eq.1)then
				  leftv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)
				  rightv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)


				  else
				  leftv(1:nof_variables)=rec_uleft(:,j,im,i)
				  rightv(1:nof_variables)=rec_uleft(:,j,im,i)
				  end if


					call get_visc_conduct(n,leftv,rightv,viscl,laml)

                              if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then

                              turbmv(1)=rec_uleftturb(1,j,im,i)

							  turbmv(2)=rec_uleftturb(1,j,im,i)
							  eddyfl(2)=turbmv(1);
							  eddyfr(2)=turbmv(2)
							  call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if

					if (turbulence .eq. 1) then
					lam_qflux=laml(3)
					else
					lam_qflux=laml(1)
					end if



				  ssx=ssx+lam_qflux*temp_grad(1)*wequa2d(im)*nx!*surface_temp
				  ssy=ssy+lam_qflux*temp_grad(2)*wequa2d(im)*ny!*surface_temp
				  ssz=ssz+lam_qflux*temp_grad(3)*wequa2d(im)*nz!*surface_temp
               end do








shear_temp=-(ssx+ssy+ssz)



end subroutine heat_x








subroutine heat_x2d(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 integer::i,k,j,kmaxe,gqi_points,nnd,im
 real,dimension(1:gpu_max_nvar)::leftv
 real,dimension(1:2)::temp_grad
real::mp_pinfl,gammal,lam_qflux
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz,surface_temp
real,dimension(1:4)::viscl,laml
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:gpu_max_dim,1:gpu_max_qp_face)::qpoints2d
real,dimension(1:gpu_max_qp_face)::wequa2d
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
real,dimension(1:20)::eddyfl,eddyfr


i=iconsidered
ssx=zero;ssy=zero
j=facex
			     nx=ielem_faceanglex(j,i)
			      ny=ielem_faceangley(j,i)

                gqi_points=qp_line_n
					   if(reduce_comp.eq.1)then
					  wequa2d=1.0d0;
					  else
					  nnd=2
				      do k=1,nnd
					vext(k,1:dims)=inoder4_cord(1:dims,ielem_nodes_faces(j,k,i))
				      end do

					  call  quadratureline(n,igqrules,vext,qpoints2d,wequa2d)
					  end if
					  surface_temp=ielem_surf(j,i)




				do im=1,gqi_points
				temp_grad(1:2)=rec_uleftv(1:2,dimensiona+1,j,im,i)

				if (dg.eq.1)then
				  leftv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)
				  rightv(1:nof_variables)=rec_uleft_dg(1:nof_variables, j,im,i)


				  else
				  leftv(1:nof_variables)=rec_uleft(:,j,im,i)
				  rightv(1:nof_variables)=rec_uleft(:,j,im,i)
				  end if

				  call get_visc_conduct(n,leftv,rightv,viscl,laml)

                              if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then

                              turbmv(1)=rec_uleftturb(1,j,im,i)

							  turbmv(2)=rec_uleftturb(1,j,im,i)
							  eddyfl(2)=turbmv(1);
							  eddyfr(2)=turbmv(2)
							  call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
						      end if

					if (turbulence .eq. 1) then
					lam_qflux=laml(3)
					else
					lam_qflux=laml(1)
					end if

                      ssx=ssx+lam_qflux*temp_grad(1)*wequa2d(im)*nx!*surface_temp
                      ssy=ssy+lam_qflux*temp_grad(2)*wequa2d(im)*ny!*surface_temp





               end do



shear_temp=-(ssx+ssy)



end subroutine heat_x2d


subroutine heat_y2d(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 integer::i,k,j,kmaxe,gqi_points,nnd,im
 real,dimension(1:gpu_max_nvar)::leftv
 real,dimension(1:3)::temp_grad
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz,surface_temp
real,dimension(1:4)::viscl,laml
real,dimension(1:8,1:gpu_max_dim)::vext
real,dimension(1:gpu_max_dim,1:gpu_max_qp_face)::qpoints2d
real,dimension(1:gpu_max_qp_face)::wequa2d


i=iconsidered
ssy=zero
j=facex
			     nx=ielem_faceanglex(j,i)
			      ny=ielem_faceangley(j,i)

                gqi_points=qp_line_n
					   if(reduce_comp.eq.1)then
					  wequa2d=1.0d0;
					  else
					  nnd=2
				      do k=1,nnd
					vext(k,1:dims)=inoder4_cord(1:dims,ielem_nodes_faces(j,k,i))
				      end do

					  call  quadratureline(n,igqrules,vext,qpoints2d,wequa2d)
					  end if
					  surface_temp=ielem_surf(j,i)




				do im=1,gqi_points
				temp_grad(1:2)=rec_uleftv(1:2,dimensiona+1,j,im,i)
				  ssy=ssy-0.026*temp_grad(2)*wequa2d(im)*surface_temp
               end do



shear_temp=ssy



end subroutine heat_y2d




















subroutine shear_x_av(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the average shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml

 
 
!  if (rungekutta.eq.4)then
! 	      ind1=7
! 	      else
! 	      ind1=5
! 	      end if

vortet1(1:3,1:3) = rec_gradsav(1:3,1:3,iconsidered)


 ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
 vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
 wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=(cos(angle1)*sin(angle2))
 ny=(sin(angle1)*sin(angle2))
 nz=(cos(angle2))
 leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
tauyx=(uy + vx)
tauzx=(wx + uz)
tauzy=(vz + wy)
ssx=(viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))
ssy=(viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))
ssz=(viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))


shear_temp=-ssx/(0.5*rres*ufreestream*ufreestream)



end subroutine shear_x_av



subroutine shear_y_av(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the average shear stresses in y-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml

 
vortet1(1:3,1:3) = rec_gradsav(1:3,1:3,iconsidered)


 ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
 vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
 wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=(cos(angle1)*sin(angle2))
 ny=(sin(angle1)*sin(angle2))
 nz=(cos(angle2))
 leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
tauyx=(uy + vx)
tauzx=(wx + uz)
tauzy=(vz + wy)
ssx=(viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))
ssy=(viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))
ssz=(viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))


shear_temp=-ssy/(0.5*rres*ufreestream*ufreestream)



end subroutine shear_y_av


subroutine shear_z_av(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the avergage shear stresses in z-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml

 
vortet1(1:3,1:3) = rec_gradsav(1:3,1:3,iconsidered)


 ux = vortet1(1,1);uy = vortet1(1,2);uz = vortet1(1,3)
 vx = vortet1(2,1);vy = vortet1(2,2);vz = vortet1(2,3)
 wx = vortet1(3,1);wy = vortet1(3,2);wz = vortet1(3,3)


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=(cos(angle1)*sin(angle2))
 ny=(sin(angle1)*sin(angle2))
 nz=(cos(angle2))
 leftv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(ind1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=(4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
tauyy=(4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
tauzz=(4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
tauyx=(uy + vx)
tauzx=(wx + uz)
tauzy=(vz + wy)
ssx=(viscl(1)*((nx*tauxx)+(ny*tauyx)+(nz*tauzx)))
ssy=(viscl(1)*((nx*tauyx)+(ny*tauyy)+(nz*tauzy)))
ssz=(viscl(1)*((nx*tauzx)+(ny*tauzy)+(nz*tauzz)))


shear_temp=-ssz/(0.5*rres*ufreestream*ufreestream)


end subroutine shear_z_av

subroutine shear_x2d(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml

integer::gqi_points,im
real,dimension(1:8,1:gpu_max_dim)::vext
	real,dimension(1:gpu_max_dim,1:gpu_max_qp_face)::qpoints2d
	real,dimension(1:gpu_max_qp_face)::wequa2d
 ssx=zero; ssp=zero; ssy=zero; 
 
 gqi_points=qp_line_n
call quadratureline(n,igqrules,vext,qpoints2d,wequa2d)
 
 do im=1,gqi_points
 if (ielem_ggs(iconsidered).eq.1)then
vortet1(1:2,1:2) = rec_grads(1:2,1:2,iconsidered)
else
 vortet1(1,1:2)=rec_uleftv(1:2,1,facex,im,iconsidered)
vortet1(2,1:2)=rec_uleftv(1:2,2,facex,im,iconsidered)

end if

 ux = vortet1(1,1);uy = vortet1(1,2)
 vx = vortet1(2,1);vy = vortet1(2,2)
 


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=angle1
 ny=angle2
 
 leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=2.0d0*ux
tauyy=2.0d0*vy

tauyx=(uy + vx)

ssx=ssx+((viscl(1)*((ny*tauyx)))*wequa2d(im))
ssy=ssy+((viscl(1)*((nx*tauyx)))*wequa2d(im))
end do

shear_temp=-ssx/(0.5*rres*ufreestream*ufreestream)




end subroutine shear_x2d



subroutine shear_y2d(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the shear stresses in y-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp
real::ux,uy,uz,vx,vy,vz,wx,wy,wz,tauxx,tauyy,tauzz,tauyx,tauzx,tauzy
 real::ssx,ssy,ssz,ssp
 real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1
 real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
integer::gqi_points,im
real,dimension(1:8,1:gpu_max_dim)::vext
	real,dimension(1:gpu_max_dim,1:gpu_max_qp_face)::qpoints2d
	real,dimension(1:gpu_max_qp_face)::wequa2d
 ssx=zero; ssp=zero; ssy=zero; 
 
 gqi_points=qp_line_n
call quadratureline(n,igqrules,vext,qpoints2d,wequa2d)
 
 do im=1,gqi_points
 if (ielem_ggs(iconsidered).eq.1)then
vortet1(1:2,1:2) = rec_grads(1:2,1:2,iconsidered)
else
 vortet1(1,1:2)=rec_uleftv(1:2,1,facex,im,iconsidered)
vortet1(2,1:2)=rec_uleftv(1:2,2,facex,im,iconsidered)

end if

 ux = vortet1(1,1);uy = vortet1(1,2)
 vx = vortet1(2,1);vy = vortet1(2,2)
 


 angle1=ielem_faceanglex(facex,iconsidered)
 angle2=ielem_faceangley(facex,iconsidered)
 nx=angle1
 ny=angle2
 
 leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
 rightv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)

 call get_visc_conduct(n,leftv,rightv,viscl,laml)

ssx=zero; ssp=zero; ssy=zero; ssz=zero

tauxx=2.0d0*ux
tauyy=2.0d0*vy

tauyx=(uy + vx)

ssx=ssx+((viscl(1)*((ny*tauyx)))*wequa2d(im))
ssy=ssy+((viscl(1)*((nx*tauyx)))*wequa2d(im))
end do

shear_temp=-ssx/(0.5*rres*ufreestream*ufreestream)


end subroutine shear_y2d




subroutine shear_x2d_av(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the average shear stresses in x-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp





shear_temp=0.0d0


end subroutine shear_x2d_av



subroutine shear_y2d_av(iconsidered,facex,shear_temp)
implicit none
!> @brief
!> this subroutine computes the average shear stresses in y-axis
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::iconsidered,facex
real,intent(inout)::shear_temp

shear_temp=0.0d0


end subroutine shear_y2d_av



subroutine get_visc_conduct(n,leftv,rightv,viscl,laml)
	implicit none
!> @brief
!> this subroutine computes the viscosity and thermal conductivity according to sutherland's law or others
#ifdef gpu
!$omp declare target
#endif
        real,dimension(1:nof_variables),intent(in)::leftv,rightv
        real,dimension(1:4),intent(inout)::viscl,laml
	integer,intent(in)::n
	real::kinetic,u,v,w,t0l,t1l,t0r,t1r
	real::mp_temp,mp_mu_mix,mp_k_mix,mp_cp_mix,gammal,gammar
	real,dimension(1:gpu_max_species)::mp_d_eff,mp_mu_i,mp_ktr_i
	real::mp_ktr_mix,mp_kve,mp_pinfl,mp_pinfr
	real,dimension(1:gpu_max_nvar)::left_temp,right_temp


          left_temp=leftv
          right_temp=rightv
	if ((multispecies.eq.1).or.(realgas.eq.1))then

          if (multispecies.eq.1)then

              call cons2prim(n,left_temp,mp_pinfl,gammal)
              call prim2cons(n,left_temp)

              call multispecies_mixtures(left_temp,mp_temp,mp_mu_mix,mp_k_mix,mp_cp_mix,gammal)

              viscl(1)=mp_mu_mix
              laml(1)=mp_k_mix

              call cons2prim(n,right_temp,mp_pinfr,gammar)
              call prim2cons(n,right_temp)
              call multispecies_mixtures(right_temp,mp_temp,mp_mu_mix,mp_k_mix,mp_cp_mix,gammar)

              viscl(2)=mp_mu_mix
              laml(2)=mp_k_mix

          end if


          if (realgas.eq.1)then
            call multispecies_mixtures_rg(left_temp,mp_mu_mix,mp_ktr_mix,mp_kve,mp_d_eff,mp_mu_i,mp_ktr_i,gammal)
								laml(1)=mp_ktr_mix
								viscl(1)=mp_mu_mix

            call multispecies_mixtures_rg(right_temp,mp_mu_mix,mp_ktr_mix,mp_kve,mp_d_eff,mp_mu_i,mp_ktr_i,gammar)
								laml(2)=mp_ktr_mix
								viscl(2)=mp_mu_mix


          end if



	else


        call cons2prim2(n,left_temp,right_temp,mp_pinfl,mp_pinfr,gammal,gammar)
		
		t1l=left_temp(dimensiona+2)/(left_temp(1)*r_gas)
		t0l=pres/(rres*r_gas)
		
		
		t1r=right_temp(dimensiona+2)/(right_temp(1)*r_gas)
		t0r=pres/(rres*r_gas)

             
	      	
              viscl(1)=visc*((t1l/t0l)*sqrt(t1l/t0l))*((t0l+(suther*t0l))/(t1l+(suther*t0l)))
              viscl(2)=visc*((t1r/t0r)*sqrt(t1r/t0r))*((t0r+(suther*t0r))/(t1r+(suther*t0r)))

	      
	      laml(1)=viscl(1)*r_gas*gamma/(prandtl*(gamma-1.d0))
	      laml(2)=viscl(2)*r_gas*gamma/(prandtl*(gamma-1.d0))
	  


	  end if
	
	     
	   

      

  end subroutine get_visc_conduct
  

subroutine get_visc_conduct_ideal(n,leftv,rightv,viscl,laml)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::leftv,rightv
real,dimension(1:4),intent(inout)::viscl,laml
real,dimension(1:gpu_max_nvar)::left_temp,right_temp
real::mp_pinfl,mp_pinfr,gammal,gammar,t0l,t1l,t0r,t1r

left_temp(1:nof_variables)=leftv(1:nof_variables)
right_temp(1:nof_variables)=rightv(1:nof_variables)
call cons2prim2_ideal(n,left_temp,right_temp,mp_pinfl,mp_pinfr,gammal,gammar)
t1l=left_temp(dimensiona+2)/(left_temp(1)*r_gas)
t1r=right_temp(dimensiona+2)/(right_temp(1)*r_gas)
t0l=pres/(rres*r_gas)
t0r=t0l
viscl(1)=visc*((t1l/t0l)*sqrt(t1l/t0l))*((t0l+(suther*t0l))/(t1l+(suther*t0l)))
viscl(2)=visc*((t1r/t0r)*sqrt(t1r/t0r))*((t0r+(suther*t0r))/(t1r+(suther*t0r)))
laml(1)=viscl(1)*r_gas*gamma/(prandtl*(gamma-1.0d0))
laml(2)=viscl(2)*r_gas*gamma/(prandtl*(gamma-1.0d0))
viscl(3)=zero
viscl(4)=zero
laml(3)=zero
laml(4)=zero

end subroutine get_visc_conduct_ideal




  
subroutine vortexcalc(n)
implicit none
!> @brief
!> this subroutine computes the q-criterion
integer,intent(in)::n
integer::kmaxe,i,ihgt,ihgj
real::snorm,onorm
real,dimension(3,3)::tvort,svort,ovort
real,dimension(1:3,1:3)::vortet1
 	 kmaxe=xmpielrank(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& firstprivate(n,kmaxe) &
!$omp& map(alloc: rec_grads, ielem_vortex) &
!$omp& private(i, ihgt, ihgj, snorm, onorm, tvort, svort, ovort, vortet1)
#else
!$omp do
#endif
do i=1,kmaxe     
	    
                vortet1(1:3,1:3)=rec_grads(1:3,1:3,i)    
	    
	    do ihgt=1,3; do ihgj=1,3
	    tvort(ihgt,ihgj)=vortet1(ihgj,ihgt)
	      end do; end do
	      svort=0.5d0*(vortet1+tvort)
	      ovort=0.5d0*(vortet1-tvort)
	      snorm=sqrt((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+&
(svort(1,3)*svort(1,3))+(svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))+(svort(2,3)*svort(2,3))&
+(svort(3,1)*svort(3,1))+(svort(3,2)*svort(3,2))+(svort(3,3)*svort(3,3)))
	      onorm=sqrt((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+(ovort(1,3)*ovort(1,3))+&
(ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2))+(ovort(2,3)*ovort(2,3))+(ovort(3,1)*ovort(3,1))+&
(ovort(3,2)*ovort(3,2))+(ovort(3,3)*ovort(3,3)))
	      
	      ielem_vortex(1,i)=(0.5d0*((onorm**2)-(snorm**2)))
		
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end subroutine vortexcalc



subroutine enstrophy_calc(n)
implicit none
!> @brief
!> this subroutine computes the q-criterion
integer,intent(in)::n
integer::kmaxe,i,ihgt,ihgj
real::snorm,onorm
real,dimension(3,3)::tvort,ovort
real::ux,uy,uz,vx,vy,vz,wx,wy,wz
real::rho_i,u_i,v_i,w_i,kinetic_i,p_i,t_i,t_ref,t_ratio,mu_i
real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1


 	 kmaxe=xmpielrank(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,boundtype,gamma,oo2,r_gas,pres,rres,visc,suther,zero) &
!$omp& map(alloc: rec_grads,u_c_val,ielem_vortex,ielem_totvolume) &
!$omp& private(i, ihgt, ihgj, snorm, onorm, tvort, ovort, vortet1) &
!$omp& private(ux, uy, uz, vx, vy, vz, wx, wy, wz) &
!$omp& private(rho_i, u_i, v_i, w_i, kinetic_i, p_i, t_i, t_ref, t_ratio, mu_i)
#else
!$omp do
#endif
do i=1,kmaxe

                vortet1(1:3,1:3)=rec_grads(1:3,1:3,i)

	    do ihgt=1,3; do ihgj=1,3
	    tvort(ihgt,ihgj)=vortet1(ihgj,ihgt)
	      end do; end do

	      ovort=(vortet1-tvort)
	      onorm=((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+(ovort(1,3)*ovort(1,3))+&
(ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2))+(ovort(2,3)*ovort(2,3))+(ovort(3,1)*ovort(3,1))+&
(ovort(3,2)*ovort(3,2))+(ovort(3,3)*ovort(3,3)))

           if(boundtype.eq.1)then

	      ielem_vortex(2,i)=(0.5d0*(onorm*u_c_val(1,1,i)))

	      else

                      ux = rec_grads(1,1,i); uy = rec_grads(1,2,i); uz = rec_grads(1,3,i);
					  vx = rec_grads(2,1,i); vy = rec_grads(2,2,i); vz = rec_grads(2,3,i);
					  wx = rec_grads(3,1,i); wy = rec_grads(3,2,i); wz = rec_grads(3,3,i);

          rho_i=u_c_val(1,1,i)
          u_i=u_c_val(1,2,i)/rho_i
          v_i=u_c_val(1,3,i)/rho_i
          w_i=zero
          if (dimensiona.eq.3) w_i=u_c_val(1,4,i)/rho_i
          kinetic_i=(u_i*u_i)+(v_i*v_i)+(w_i*w_i)
          p_i=(gamma-1.0d0)*(u_c_val(1,dimensiona+2,i)-oo2*rho_i*kinetic_i)
          t_i=p_i/(rho_i*r_gas)
          t_ref=pres/(rres*r_gas)
          t_ratio=t_i/t_ref
          mu_i=visc*(t_ratio*sqrt(t_ratio))*((t_ref+(suther*t_ref))/(t_i+(suther*t_ref)))

          snorm=((wy-vz)**2)+((uz-wx)**2)+((vx-uy)**2)
          ielem_vortex(2,i)=ielem_totvolume(i)*snorm*mu_i
          ielem_vortex(3,i)=(4.0/3.0)*mu_i*((ux+vy+wz)**2)*ielem_totvolume(i)





	      end if

end do

#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end subroutine enstrophy_calc


subroutine vortexcalc2d(n)
implicit none
!> @brief
!> this subroutine computes the q criterion for 2d
integer,intent(in)::n
integer::kmaxe,i,ihgt,ihgj
real::snorm,onorm
real,dimension(2,2)::tvort,svort,ovort
real,dimension(1:gpu_max_dim,1:gpu_max_dim)::vortet1

 	 kmaxe=xmpielrank(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& firstprivate(n,kmaxe) &
!$omp& map(alloc: rec_grads, ielem_vortex) &
!$omp& private(i, ihgt, ihgj, snorm, onorm, tvort, svort, ovort, vortet1)
#else
!$omp do
#endif
do i=1,kmaxe     
	    
                vortet1(1:2,1:2)=rec_grads(1:2,1:2,i)    
	    
	    do ihgt=1,2; do ihgj=1,2
	    tvort(ihgt,ihgj)=vortet1(ihgj,ihgt)
	      end do; end do
	      svort(1:2,1:2)=0.5d0*(vortet1(1:2,1:2)+tvort(1:2,1:2))
	      ovort(1:2,1:2)=0.5d0*(vortet1(1:2,1:2)-tvort(1:2,1:2))
	      snorm=sqrt((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+&
 (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2)))
	      onorm=sqrt((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+&
(ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2)))
	      
	      ielem_vortex(1,i)=(0.5d0*((onorm**2)-(snorm**2)))
		
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end subroutine vortexcalc2d





subroutine boundarys_realgas_outflow(n,leftv,rightv,nx,ny,nz)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::rightv
real,intent(in)::nx,ny,nz
real,dimension(1:gpu_max_nvar)::tempvect1,tempvect2
real,dimension(1:gpu_max_species)::y_s_g
real::mp_pinfl,gammal,agrt,u,v,w,un,mn,rho_g,rg_rmix,rmix,t_g,tv_g
real::rg_tr,rg_ev_total,rg_chem,skin1,rg_cv
integer::rg_i

tempvect1=zero
tempvect2=zero
tempvect1(1:nof_variables)=leftv(1:nof_variables)
tempvect2(1:nof_variables)=leftv(1:nof_variables)

call cons2prim(n,tempvect1,mp_pinfl,gammal)
agrt=sqrt((tempvect1(5)+mp_pinfl)*gammal/tempvect1(1))
u=tempvect1(2)
v=tempvect1(3)
w=tempvect1(4)
un=u*nx+v*ny+w*nz

call cons2div(n,tempvect2,mp_pinfl,gammal)

if (un.lt.0.0d0)then
  rightv=leftv
else
  mn=un/agrt

  if (mn.ge.1.0d0)then
    rightv(1:nof_variables)=leftv(1:nof_variables)
  else
    rho_g=tempvect1(1)
    y_s_g(1:nof_species)=tempvect1(dimensiona+4:nof_variables)

    rg_rmix=zero
    do rg_i=1,nof_species
      rg_rmix=rg_rmix+y_s_g(rg_i)/rg_molm(rg_i)
    end do

    rg_rmix=rgs_ru*rg_rmix
    rmix=rg_rmix
    t_g=pres/(rho_g*rmix)
    tv_g=tempvect2(6)

    rg_tr=zero
    do rg_i=1,nof_species
      if (rg_i.le.3)then
        rg_cv=(5.0d0/2.0d0)*rgs_ru/rg_molm(rg_i)
      else
        rg_cv=(3.0d0/2.0d0)*rgs_ru/rg_molm(rg_i)
      end if
      rg_tr=rg_tr+y_s_g(rg_i)*rg_cv*t_g
    end do

    rg_ev_total=zero
    do rg_i=1,3
      rg_ev_total=rg_ev_total+y_s_g(rg_i)*(rgs_ru/rg_molm(rg_i))*(rg_thetag(rg_i)/(exp(rg_thetag(rg_i)/tv_g)-1.0d0))
    end do

    rg_chem=zero
    do rg_i=1,nof_species
      if (rg_hzero(rg_i).gt.1.0e-12)then
        rg_chem=rg_chem-(y_s_g(rg_i)*rg_hzero(rg_i)/rg_molm(rg_i))
      end if
    end do

    skin1=oo2*((u**2)+(v**2)+(w*w))

    rightv(1)=rho_g
    rightv(2)=rho_g*u
    rightv(3)=rho_g*v
    rightv(4)=rho_g*w
    rightv(5)=rho_g*(rg_ev_total+rg_tr+rg_chem+skin1)
    rightv(6)=rho_g*rg_ev_total

    do rg_i=1,nof_species
      rightv(6+rg_i)=rho_g*y_s_g(rg_i)
    end do
  end if
end if

end subroutine boundarys_realgas_outflow


subroutine boundarys2d_realgas_outflow(n,leftv,rightv,nx,ny)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::rightv
real,intent(in)::nx,ny
real,dimension(1:gpu_max_nvar)::tempvect1,tempvect2
real,dimension(1:gpu_max_species)::y_s_g
real::mp_pinfl,gammal,agrt,u,v,un,mn,rho_g,rg_rmix,rmix,t_g,tv_g
real::rg_tr,rg_ev_total,rg_chem,skin1,rg_cv
integer::rg_i

tempvect1=zero
tempvect2=zero
tempvect1(1:nof_variables)=leftv(1:nof_variables)
tempvect2(1:nof_variables)=leftv(1:nof_variables)

call cons2prim(n,tempvect1,mp_pinfl,gammal)
agrt=sqrt((tempvect1(4)+mp_pinfl)*gammal/tempvect1(1))
u=tempvect1(2)
v=tempvect1(3)
un=u*nx+v*ny

call cons2div(n,tempvect2,mp_pinfl,gammal)

if (un.lt.0.0d0)then
  rightv=leftv
else
  mn=un/agrt

  if (mn.ge.1.0d0)then
    rightv(1:nof_variables)=leftv(1:nof_variables)
  else
    rho_g=tempvect1(1)
    y_s_g(1:nof_species)=tempvect1(dimensiona+4:nof_variables)

    rg_rmix=zero
    do rg_i=1,nof_species
      rg_rmix=rg_rmix+y_s_g(rg_i)/rg_molm(rg_i)
    end do

    rg_rmix=rgs_ru*rg_rmix
    rmix=rg_rmix
    t_g=pres/(rho_g*rmix)
    tv_g=tempvect2(5)

    rg_tr=zero
    do rg_i=1,nof_species
      if (rg_i.le.3)then
        rg_cv=(5.0d0/2.0d0)*rgs_ru/rg_molm(rg_i)
      else
        rg_cv=(3.0d0/2.0d0)*rgs_ru/rg_molm(rg_i)
      end if
      rg_tr=rg_tr+y_s_g(rg_i)*rg_cv*t_g
    end do

    rg_ev_total=zero
    do rg_i=1,3
      rg_ev_total=rg_ev_total+y_s_g(rg_i)*(rgs_ru/rg_molm(rg_i))*(rg_thetag(rg_i)/(exp(rg_thetag(rg_i)/tv_g)-1.0d0))
    end do

    rg_chem=zero
    do rg_i=1,nof_species
      if (rg_hzero(rg_i).gt.1.0e-12)then
        rg_chem=rg_chem-(y_s_g(rg_i)*rg_hzero(rg_i)/rg_molm(rg_i))
      end if
    end do

    skin1=oo2*((u**2)+(v**2))

    rightv(1)=rho_g
    rightv(2)=rho_g*u
    rightv(3)=rho_g*v
    rightv(4)=rho_g*(rg_ev_total+rg_tr+rg_chem+skin1)
    rightv(5)=rho_g*rg_ev_total

    do rg_i=1,nof_species
      rightv(5+rg_i)=rho_g*y_s_g(rg_i)
    end do
  end if
end if

end subroutine boundarys2d_realgas_outflow


subroutine set_state_ideal3d(r,u,v,w,p,state)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,intent(in)::r,u,v,w,p
real,dimension(1:nof_variables),intent(out)::state
real::skin,ien

skin=oo2*((u*u)+(v*v)+(w*w))
ien=p/((gamma-1.0d0)*r)
state(1)=r
state(2)=r*u
state(3)=r*v
state(4)=r*w
state(5)=r*(skin+ien)

end subroutine set_state_ideal3d


subroutine set_state_ideal2d(r,u,v,p,state)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,intent(in)::r,u,v,p
real,dimension(1:nof_variables),intent(out)::state
real::skin,ien

skin=oo2*((u*u)+(v*v))
ien=p/((gamma-1.0d0)*r)
state(1)=r
state(2)=r*u
state(3)=r*v
state(4)=r*(skin+ien)

end subroutine set_state_ideal2d


subroutine outlet_pressure_fix_ideal3d(n,leftv,rightv,nx,ny,nz,pout)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,intent(in)::nx,ny,nz,pout
real,dimension(1:gpu_max_nvar)::subson2,subson3
real::mp_pinfl,mp_pinfr,gammal,gammar,sps,vel,skins,ikins

rightv(1:nof_variables)=leftv(1:nof_variables)
call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
subson2(1:nof_variables)=leftv(1:nof_variables)
sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
vel=subson2(2)*nx+subson2(3)*ny+subson2(4)*nz
call prim2cons2_ideal(n,leftv,rightv)

if ((vel.gt.0.0d0).and.(vel/(sps+tolsmall).le.1.0d0))then
  subson3(1:nof_variables)=subson2(1:nof_variables)
  subson3(5)=max(pout,tolsmall)
  subson3(1)=subson2(1)+(subson3(5)-subson2(5))/(sps*sps+tolsmall)
  subson3(2)=subson2(2)+nx*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
  subson3(3)=subson2(3)+ny*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
  subson3(4)=subson2(4)+nz*(subson2(5)-subson3(5))/((sps+tolsmall)*(subson2(1)+tolsmall))
  subson3(1)=max(subson3(1),tolsmall)
  skins=oo2*((subson3(2)*subson3(2))+(subson3(3)*subson3(3))+(subson3(4)*subson3(4)))
  ikins=subson3(5)/((gamma-1.0d0)*subson3(1))
  rightv(1)=subson3(1)
  rightv(2)=subson3(1)*subson3(2)
  rightv(3)=subson3(1)*subson3(3)
  rightv(4)=subson3(1)*subson3(4)
  rightv(5)=subson3(1)*(ikins+skins)
end if

end subroutine outlet_pressure_fix_ideal3d


subroutine outlet_pressure_fix_ideal2d(n,leftv,rightv,nx,ny,pout)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,intent(in)::nx,ny,pout
real,dimension(1:gpu_max_nvar)::subson2,subson3
real::mp_pinfl,mp_pinfr,gammal,gammar,sps,vel,skins,ikins

rightv(1:nof_variables)=leftv(1:nof_variables)
call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
subson2(1:nof_variables)=leftv(1:nof_variables)
sps=sqrt((gamma*subson2(4))/(subson2(1)+tolsmall))
vel=subson2(2)*nx+subson2(3)*ny
call prim2cons2_ideal(n,leftv,rightv)

if ((vel.gt.0.0d0).and.(vel/(sps+tolsmall).le.1.0d0))then
  subson3(1:nof_variables)=subson2(1:nof_variables)
  subson3(4)=max(pout,tolsmall)
  subson3(1)=subson2(1)+(subson3(4)-subson2(4))/(sps*sps+tolsmall)
  subson3(2)=subson2(2)+nx*(subson2(4)-subson3(4))/((sps+tolsmall)*(subson2(1)+tolsmall))
  subson3(3)=subson2(3)+ny*(subson2(4)-subson3(4))/((sps+tolsmall)*(subson2(1)+tolsmall))
  subson3(1)=max(subson3(1),tolsmall)
  skins=oo2*((subson3(2)*subson3(2))+(subson3(3)*subson3(3)))
  ikins=subson3(4)/((gamma-1.0d0)*subson3(1))
  rightv(1)=subson3(1)
  rightv(2)=subson3(1)*subson3(2)
  rightv(3)=subson3(1)*subson3(3)
  rightv(4)=subson3(1)*(ikins+skins)
end if

end subroutine outlet_pressure_fix_ideal2d


subroutine update_outlet_mach_controller(n)
implicit none
integer,intent(in)::n
integer::ii,i,l,ibc
real::mach_sum_l,area_sum_l,p_sum_l,mach_sum_g,area_sum_g,p_sum_g
real::rho,uu,vv,ww,e,vel2,p,a,mach,area,err,limited_err,factor
logical::heref

if ((mach_outlet_target.le.0.0d0).or.(mach_outlet_update_freq.le.0)) return
if ((realgas.ne.0).or.(multispecies.ne.0)) return

mach_sum_l=0.0d0
area_sum_l=0.0d0
p_sum_l=0.0d0

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do reduction(+:mach_sum_l,area_sum_l,p_sum_l) &
!$omp& firstprivate(nof_bounded, dimensiona, gamma, tolsmall, oo2) &
!$omp& map(alloc: el_bnd, ielem_ifca, ielem_ibounds, ibound_icode, ielem_surf, u_c_val) &
!$omp& private(ii,i,l,ibc,rho,uu,vv,ww,e,vel2,p,a,mach,area)
#endif
do ii=1,nof_bounded
  i=el_bnd(ii)
  do l=1,ielem_ifca(i)
    ibc=ielem_ibounds(l,i)
    if (ibc.gt.0)then
      if (ibound_icode(ibc).eq.9)then
        rho=max(u_c_val(1,1,i),tolsmall)
        uu=u_c_val(1,2,i)/rho
        vv=u_c_val(1,3,i)/rho
        if (dimensiona.eq.3)then
          ww=u_c_val(1,4,i)/rho
          e=u_c_val(1,5,i)
        else
          ww=0.0d0
          e=u_c_val(1,4,i)
        end if
        vel2=(uu*uu)+(vv*vv)+(ww*ww)
        p=(gamma-1.0d0)*(e-rho*oo2*vel2)
        if (p.gt.tolsmall)then
          a=sqrt(gamma*p/rho)
          mach=sqrt(max(vel2,0.0d0))/(a+tolsmall)
          area=max(ielem_surf(l,i),0.0d0)
          mach_sum_l=mach_sum_l+(mach*area)
          area_sum_l=area_sum_l+area
          p_sum_l=p_sum_l+(p*area)
        end if
      end if
    end if
  end do
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#endif

mach_sum_g=0.0d0
area_sum_g=0.0d0
p_sum_g=0.0d0
call mpi_allreduce(mach_sum_l,mach_sum_g,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(area_sum_l,area_sum_g,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(p_sum_l,p_sum_g,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

if (area_sum_g.le.tolsmall) return

mach_outlet_average=mach_sum_g/area_sum_g
if (press_outlet.le.tolsmall) press_outlet=max(p_sum_g/(area_sum_g+tolsmall),tolsmall)

err=mach_outlet_average-mach_outlet_target
limited_err=max(min(mach_outlet_relax*err,0.2d0),-0.2d0)
factor=exp(limited_err)
press_outlet=max(press_outlet*factor,tolsmall)

#if defined(gpu) || defined(xpu)
!$omp target update to(press_outlet)
#endif

if (n.eq.0)then
  inquire(file='OUTLET_MACH.dat',exist=heref)
  open(91,file='OUTLET_MACH.dat',form='formatted',status='unknown',action='write',position='append')
  if (.not.heref) write(91,'(A)') '# it t mach_area target_mach press_outlet outlet_area'
  write(91,'(I10,5(1X,ES20.10))') it,t,mach_outlet_average,mach_outlet_target,press_outlet,area_sum_g
  flush(91)
  close(91)
end if

end subroutine update_outlet_mach_controller


subroutine inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,inflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::leftv
real,intent(in)::nx,ny,nz
real,dimension(1:nof_variables),intent(out)::inflow_out
real::rho_l,u_l,v_l,w_l,e_l,p_l,mach
real::pt0,t0,trat,temp,r,p,a,vel,rgas_local,pratio,nmag

rho_l=max(leftv(1),tolsmall)
u_l=leftv(2)/rho_l
v_l=leftv(3)/rho_l
w_l=leftv(4)/rho_l
e_l=leftv(5)
p_l=(gamma-1.0d0)*(e_l-rho_l*oo2*((u_l*u_l)+(v_l*v_l)+(w_l*w_l)))
if (p_l.le.tolsmall) p_l=max(pres,tolsmall)

pt0=max(total_pressure_inlet,tolsmall)
t0=max(total_temperature_inlet,tolsmall)
rgas_local=max(r_gas,tolsmall)
p=max(min(p_l,pt0),tolsmall)
pratio=max(pt0/(p+tolsmall),1.0d0)
mach=sqrt(max((pratio**((gamma-1.0d0)/gamma)-1.0d0)*(2.0d0/(gamma-1.0d0)),0.0d0))
mach=min(max(mach,0.0d0),0.999d0)
trat=1.0d0+0.5d0*(gamma-1.0d0)*mach*mach
temp=t0/trat
p=pt0/(trat**(gamma/(gamma-1.0d0)))
r=max(p/(rgas_local*temp),tolsmall)
a=sqrt(gamma*rgas_local*temp)
vel=mach*a
nmag=max(sqrt((nx*nx)+(ny*ny)+(nz*nz)),tolsmall)
call set_state_ideal3d(r,-vel*nx/nmag,-vel*ny/nmag,-vel*nz/nmag,p,inflow_out)

end subroutine inflow_4440_from_left_ideal3d


subroutine inflow_4440_from_left_ideal2d(n,leftv,nx,ny,inflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::leftv
real,intent(in)::nx,ny
real,dimension(1:nof_variables),intent(out)::inflow_out
real::rho_l,u_l,v_l,e_l,p_l,mach
real::pt0,t0,trat,temp,r,p,a,vel,rgas_local,pratio,nmag

rho_l=max(leftv(1),tolsmall)
u_l=leftv(2)/rho_l
v_l=leftv(3)/rho_l
e_l=leftv(4)
p_l=(gamma-1.0d0)*(e_l-rho_l*oo2*((u_l*u_l)+(v_l*v_l)))
if (p_l.le.tolsmall) p_l=max(pres,tolsmall)

pt0=max(total_pressure_inlet,tolsmall)
t0=max(total_temperature_inlet,tolsmall)
rgas_local=max(r_gas,tolsmall)
p=max(min(p_l,pt0),tolsmall)
pratio=max(pt0/(p+tolsmall),1.0d0)
mach=sqrt(max((pratio**((gamma-1.0d0)/gamma)-1.0d0)*(2.0d0/(gamma-1.0d0)),0.0d0))
mach=min(max(mach,0.0d0),0.999d0)
trat=1.0d0+0.5d0*(gamma-1.0d0)*mach*mach
temp=t0/trat
p=pt0/(trat**(gamma/(gamma-1.0d0)))
r=max(p/(rgas_local*temp),tolsmall)
a=sqrt(gamma*rgas_local*temp)
vel=mach*a
nmag=max(sqrt((nx*nx)+(ny*ny)),tolsmall)
call set_state_ideal2d(r,-vel*nx/nmag,-vel*ny/nmag,p,inflow_out)

end subroutine inflow_4440_from_left_ideal2d


subroutine inflow_4440_sa_value(n,inflow_state,sa_value)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::inflow_state
real,intent(out)::sa_value
real,dimension(1:4)::viscl,laml

call get_visc_conduct_ideal(n,inflow_state,inflow_state,viscl,laml)
sa_value=max(viscl(1)*turbinit,tolsmall)

end subroutine inflow_4440_sa_value


subroutine inflow_ideal3d(initcond,pox,poy,poz,inflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,dimension(1:nof_variables),intent(out)::inflow_out
real::r,u,v,w,p,xf,yf,zf,theta_0,vtang,vradial

r=rres
u=uvel
v=vvel
w=wvel
p=pres
if (initcond.eq.10000)then
  if (sqrt((poy(1)*poy(1))+((poz(1)-0.5d0)*(poz(1)-0.5d0))).le.0.05d0)then
    p=0.4127d0
    r=5.0d0
    u=30.0d0
    v=0.0d0
    w=0.0d0
  end if
end if
if (swirl.eq.1)then
  if (pox(1).lt.-0.03d0)then
    xf=pox(1)
    yf=poy(1)
    zf=poz(1)
    theta_0=atan2(zf,yf)
    vtang=18.0375d0
    vradial=-12.63d0
    u=0.0d0
    v=-vtang*sin(theta_0)+vradial*cos(theta_0)
    w=vtang*cos(theta_0)+vradial*sin(theta_0)
  else
    u=70.06d0
    v=0.0d0
    w=0.0d0
  end if
  r=rres
  p=pres
end if
call set_state_ideal3d(r,u,v,w,p,inflow_out)

end subroutine inflow_ideal3d


subroutine outflow_ideal3d(initcond,pox,poy,poz,outflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,dimension(1:nof_variables),intent(out)::outflow_out
real::p

p=pres
if (initcond.eq.977) p=101325.0d0
call set_state_ideal3d(rres,uvel,vvel,wvel,p,outflow_out)

end subroutine outflow_ideal3d


subroutine outflow2_ideal3d(initcond,pox,poy,poz,outflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,dimension(1:nof_variables),intent(out)::outflow_out

call set_state_ideal3d(rres,uvel,vvel,wvel,press_outlet,outflow_out)

end subroutine outflow2_ideal3d


subroutine inflow_ideal2d(initcond,pox,poy,inflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy
real,dimension(1:nof_variables),intent(out)::inflow_out
real::r,u,v,w,p

r=rres
u=uvel
v=vvel
w=0.0d0
p=pres
call set_state_ideal2d(r,u,v,p,inflow_out)

end subroutine inflow_ideal2d


subroutine outflow_ideal2d(initcond,pox,poy,outflow_out)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::initcond
real,dimension(1:dimensiona),intent(in)::pox,poy
real,dimension(1:nof_variables),intent(out)::outflow_out

call set_state_ideal2d(rres,uvel,vvel,pres,outflow_out)

end subroutine outflow_ideal2d


subroutine boundarys_ideal(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n,b_code,iconsidered,facex
integer,intent(inout)::ibfc
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,dimension(1:nof_variables),intent(in)::srf_speedrot,srf_speed
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,intent(in)::angle1,angle2,nx,ny,nz
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cright_rot,cleft_rot
real,dimension(1:gpu_max_nvar)::subson1,subson2,subson3
real::mp_pinfl,mp_pinfr,gammal,gammar,sps,skins,ikins,vel,vnb,sa_bc

select case(b_code)
case(1)
  if (initcond.eq.4440)then
    call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
  else
    call inflow_ideal3d(initcond,pox,poy,poz,rightv)
    if (boundtype.ne.0)then
	    call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
	    subson1(1:nof_variables)=rightv(1:nof_variables)
	    subson2(1:nof_variables)=leftv(1:nof_variables)
    sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
    vel=sqrt((subson2(2)*subson2(2))+(subson2(3)*subson2(3))+(subson2(4)*subson2(4)))
    call prim2cons2_ideal(n,leftv,rightv)
    if (vel/(sps+tolsmall).le.1.0d0)then
      subson3(5)=0.5d0*(subson1(5)+subson2(5)-subson2(1)*sps*(nx*(subson1(2)-subson2(2))+ny*(subson1(3)-subson2(3))+nz*(subson1(4)-subson2(4))))
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
    rightv(1:nof_variables)=leftv(1:nof_variables)
  else
    call outflow_ideal3d(initcond,pox,poy,poz,rightv)
    call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    subson1(1:nof_variables)=rightv(1:nof_variables)
    subson2(1:nof_variables)=leftv(1:nof_variables)
    sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
    vel=subson2(2)*nx+subson2(3)*ny+subson2(4)*nz
    call prim2cons2_ideal(n,leftv,rightv)
    if (vel.gt.0.0d0)then
      if (vel/(sps+tolsmall).gt.1.0d0)then
        rightv(1:nof_variables)=leftv(1:nof_variables)
      else
        subson3(1:nof_variables)=subson2(1:nof_variables)
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
  if (rec_mrf(iconsidered).eq.1)then
    rightv(1)=leftv(1)
    rightv(2)=-leftv(2)+2.0d0*leftv(1)*srf_speed(2)
    rightv(3)=-leftv(3)+2.0d0*leftv(1)*srf_speed(3)
    rightv(4)=-leftv(4)+2.0d0*leftv(1)*srf_speed(4)
    rightv(5)=leftv(5)+2.0d0*leftv(1)*((srf_speed(2)*srf_speed(2))+(srf_speed(3)*srf_speed(3))+(srf_speed(4)*srf_speed(4))) &
              -2.0d0*(leftv(2)*srf_speed(2)+leftv(3)*srf_speed(3)+leftv(4)*srf_speed(4))
  else
    rightv(1:nof_variables)=leftv(1:nof_variables)
    rightv(2)=-leftv(2)
    rightv(3)=-leftv(3)
    rightv(4)=-leftv(4)
  end if
case(6)
  call rotatef(n,cleft_rot,leftv,angle1,angle2)
  vnb=cleft_rot(2)
  call outflow_ideal3d(initcond,pox,poy,poz,rightv)
  call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
  subson1(1:nof_variables)=rightv(1:nof_variables)
  subson2(1:nof_variables)=leftv(1:nof_variables)
  sps=sqrt((gamma*subson2(5))/(subson2(1)+tolsmall))
  call prim2cons2_ideal(n,leftv,rightv)
  if (vnb.le.0.0d0)then
    ibfc=-1
    if (initcond.eq.4440)then
      call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
    else
      call inflow_ideal3d(initcond,pox,poy,poz,rightv)
    end if
  else
    ibfc=-2
    if ((abs(vnb)).ge.sps)then
      rightv(1:nof_variables)=leftv(1:nof_variables)
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
  rightv(1:nof_variables)=leftv(1:nof_variables)
end select
if ((turbulence.eq.1).or.(passivescalar.gt.0)) cturbr(:)=cturbl(:)
if ((b_code.eq.1).and.(initcond.eq.4440).and.(turbulence.eq.1).and.(turbulencemodel.eq.1))then
  call inflow_4440_sa_value(n,rightv,sa_bc)
  cturbr(1)=sa_bc
end if

end subroutine boundarys_ideal


subroutine boundarys2d_ideal(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n,b_code,iconsidered,facex
integer,intent(inout)::ibfc
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,dimension(1:nof_variables),intent(in)::srf_speedrot,srf_speed
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,intent(in)::angle1,angle2,nx,ny,nz
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cright_rot,cleft_rot
real,dimension(1:gpu_max_nvar)::subson1,subson2,subson3
real::mp_pinfl,mp_pinfr,gammal,gammar,sps,skins,ikins,vel,vnb,sa_bc

select case(b_code)
case(1)
  if (initcond.eq.4440)then
    call inflow_4440_from_left_ideal2d(n,leftv,nx,ny,rightv)
  else
    call inflow_ideal2d(initcond,pox,poy,rightv)
  end if
case(2)
  if (boundtype.eq.0)then
    rightv(1:nof_variables)=leftv(1:nof_variables)
  else
    call outflow_ideal2d(initcond,pox,poy,rightv)
    call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    subson1(1:nof_variables)=rightv(1:nof_variables)
    subson2(1:nof_variables)=leftv(1:nof_variables)
    sps=sqrt((gamma*subson2(4))/(subson2(1)+tolsmall))
    vel=subson2(2)*nx+subson2(3)*ny
    call prim2cons2_ideal(n,leftv,rightv)
    if ((vel.gt.0.0d0).and.(vel/(sps+tolsmall).le.1.0d0))then
      subson3(4)=subson1(4)
      subson3(1)=subson2(1)+(subson3(4)-subson2(4))/(sps*sps+tolsmall)
      subson3(2)=subson2(2)+nx*(subson2(4)-subson3(4))/((sps+tolsmall)*(subson2(1)+tolsmall))
      subson3(3)=subson2(3)+ny*(subson2(4)-subson3(4))/((sps+tolsmall)*(subson2(1)+tolsmall))
      call set_state_ideal2d(subson3(1),subson3(2),subson3(3),subson3(4),rightv)
    else if (vel.gt.0.0d0)then
      rightv(1:nof_variables)=leftv(1:nof_variables)
    end if
  end if
case(3)
  call rotatef2d(n,cleft_rot,leftv,angle1,angle2)
  cright_rot(1)=cleft_rot(1)
  cright_rot(2)=-cleft_rot(2)
  cright_rot(3)=cleft_rot(3)
  cright_rot(4)=cleft_rot(4)
  call rotateb2d(n,rightv,cright_rot,angle1,angle2)
case(4)
  rightv(1:nof_variables)=leftv(1:nof_variables)
  rightv(2)=-leftv(2)
  rightv(3)=-leftv(3)
case(6)
  call rotatef2d(n,cleft_rot,leftv,angle1,angle2)
  vnb=cleft_rot(2)/cleft_rot(1)
  if (vnb.le.0.0d0)then
    ibfc=-1
    if (initcond.eq.4440)then
      call inflow_4440_from_left_ideal2d(n,leftv,nx,ny,rightv)
    else
      call inflow_ideal2d(initcond,pox,poy,rightv)
    end if
	  else
	    ibfc=-2
	    call outflow_ideal2d(initcond,pox,poy,rightv)
	  end if
		case(9)
		  call outlet_pressure_fix_ideal2d(n,leftv,rightv,nx,ny,press_outlet)
	case default
	  rightv(1:nof_variables)=leftv(1:nof_variables)
end select
if ((turbulence.eq.1).or.(passivescalar.gt.0)) cturbr(:)=cturbl(:)
if ((b_code.eq.1).and.(initcond.eq.4440).and.(turbulence.eq.1).and.(turbulencemodel.eq.1))then
  call inflow_4440_sa_value(n,rightv,sa_bc)
  cturbr(1)=sa_bc
end if

end subroutine boundarys2d_ideal


subroutine boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
implicit none
!> @brief
!> this subroutine applies the boundary condition to each bounded cell
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,b_code,iconsidered,facex
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
integer,intent(inout)::ibfc
real,dimension(1:nof_variables),intent(in)::srf_speedrot,srf_speed
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,intent(in)::angle1,angle2,nx,ny,nz
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cright_rot,cleft_rot
real,dimension(1:gpu_max_nvar)::subson1,subson2,subson3
real::mp_pinfl,mp_pinfr,gammal,gammar,u,v,w,sa_bc
real::sps,skins,ikins,vel,vnb,theeta,reeta
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx,vhx,amp,dvel,rgg,tt1
real::mn,un,agrt,rho_g,rg_rmix,rmix,t_g,tv_g,rg_tr,rg_ev_total,rg_chem
real,dimension(1:gpu_max_passive)::pass_tmp
integer::rg_i





select case(b_code)

	    
	    case(1)!inflow subsonic or supersonic will be chosen based on mach number
	    if ((initcond.eq.4440).and.(realgas.eq.0).and.(multispecies.eq.0))then
	    call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
	    if ((turbulence.eq.1).or.(passivescalar.gt.0)) cturbr(:)=cturbl(:)
	    if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
	      call inflow_4440_sa_value(n,rightv,sa_bc)
	      cturbr(1)=sa_bc
	    end if
	    else
	    if (boundtype.eq.0)then	!supersonic
    
    call inflow(initcond,pox,poy,poz,rightv)
    
          if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit*rightv(1)
	      end if
    
        
    
    
    else		!subsonic
    call inflow(initcond,pox,poy,poz,rightv)
    
    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    
    subson1(1:nof_variables)=rightv(1:nof_variables)
    subson2(1:nof_variables)=leftv(1:nof_variables)
    sps=sqrt((gamma*subson2(5))/(subson2(1)))
    vel=sqrt(subson2(2)**2+subson2(3)**2+subson2(4)**2)
    call prim2cons2(n,leftv,rightv)
    
      if (vel/(sps+tolsmall).gt.1.0d0)then	!supersonic
      
      
      call inflow(initcond,pox,poy,poz,rightv)

       if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit*rightv(1)
	      end if
      
      
      
      else		!subsonic
      
      subson3(5)=0.5*((subson1(5))+(subson2(5))-(subson2(1)*sps*((nx*(subson1(2)-subson2(2)))+(ny*(subson1(3)-subson2(3)))&
+(nz*(subson1(4)-subson2(4))))))
      subson3(1)=subson1(1)+(subson3(5)-subson1(5))/(sps**2)
      subson3(2)=subson1(2)-(nx*(subson1(5)-subson3(5)))/(sps*subson2(1))
      subson3(3)=subson1(3)-(ny*(subson1(5)-subson3(5)))/(sps*subson2(1))
      subson3(4)=subson1(4)-(nz*(subson1(5)-subson3(5)))/(sps*subson2(1))
      
     
      
       rightv(1)=subson3(1)
    rightv(2)=subson3(2)*subson3(1)
    rightv(3)=subson3(3)*subson3(1)
    rightv(4)=subson3(4)*subson3(1)
    skins=oo2*((subson3(2)**2)+(subson3(3)**2)+(subson3(4)**2))
    ikins=subson3(5)/((gamma-1.0d0)*(subson3(1)))
    rightv(5)=(subson3(1)*(ikins))+(subson3(1)*skins)
      
      
      
      end if

      
      
      
      
      
      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      
      
      
      
	      if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit*rightv(1)
	      end if
	      if ((turbulence.eq.1).and.(turbulencemodel.eq.2))then	 
		cturbr(1)=(1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
		cturbr(2)=rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
        if (rec_mrf(iconsidered).eq.1)then
        cturbr(1)=(1.5d0*i_turb_inlet*(kinit_srf**2))*rightv(1)!k initialization
		cturbr(2)=rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
        end if
	      end if
 
	      if (passivescalar.gt.0)then
	      call pass_inlet(initcond,pox,poy,poz,pass_tmp)
	      do rg_i=1,passivescalar
	        cturbr(turbulenceequations+rg_i)=pass_tmp(rg_i)*rightv(1)
	      end do
	      end if
      end if
      
      
      
      
      
      
      
      
    
	    end if
	    end if
	    
	    
	    case(2)!outflow subsonic or supersonic will be chosen based on mach number
     if (boundtype.eq.0)then
      rightv(1:nof_variables)=leftv(1:nof_variables)
      
     else

      if (realgas.eq.1)then
                  call boundarys_realgas_outflow(n,leftv,rightv,nx,ny,nz)
                  else  !realgas

                                    ! ------------------------------------------------------------
                                    ! Build a reference/outlet state (used for prescribed pressure,
                                    ! and for backflow if it happens)
                                    ! ------------------------------------------------------------
                                    call outflow(initcond,pox,poy,poz,rightv)

                                    ! Convert interior/exterior to primitive for BC logic
                                    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

                                    subson2(1:nof_variables) = leftv(1:nof_variables)   ! interior primitive
                                    subson1(1:nof_variables) = rightv(1:nof_variables)  ! "target" primitive (from outflow)

                                    ! Speed of sound from interior
                                    sps = sqrt( (gamma*subson2(5)) / (subson2(1) + tolsmall) )

                                    ! Normal velocity from interior
                                    vel = subson2(2)*nx + subson2(3)*ny + subson2(4)*nz

                                    ! ------------------------------------------------------------
                                    ! Decide boundary type
                                    ! ------------------------------------------------------------

                                    ! 1) Backflow / inflow at outlet: treat as inflow (use reference state)
                                    if (vel .le. 0.0d0) then

                                      ! Keep the primitive outlet state from outflow() (subson1)
                                      subson3(1:nof_variables) = subson1(1:nof_variables)

                                    else

                                      ! 2) Outflow: check normal Mach number
                                      if ( vel/(sps + tolsmall) .gt. 1.0d0 ) then
                                          ! Supersonic outflow: extrapolate (copy interior primitive)
                                          subson3(1:nof_variables) = subson2(1:nof_variables)

                                      else
                                          ! Subsonic outflow: impose pressure (from reference state),
                                          ! correct rho and velocity using linearized characteristics
                                          subson3(1:nof_variables) = subson2(1:nof_variables)

                                          subson3(5) = subson1(5)   ! prescribed outlet static pressure (e.g. p_inf)

                                          subson3(1) = subson2(1) + (subson3(5) - subson2(5)) / (sps*sps + tolsmall)

                                          subson3(2) = subson2(2) + (nx*(subson2(5) - subson3(5))) / ((sps + tolsmall)*(subson2(1) + tolsmall))
                                          subson3(3) = subson2(3) + (ny*(subson2(5) - subson3(5))) / ((sps + tolsmall)*(subson2(1) + tolsmall))
                                          subson3(4) = subson2(4) + (nz*(subson2(5) - subson3(5))) / ((sps + tolsmall)*(subson2(1) + tolsmall))

                                      end if
                                    end if

                                    ! ------------------------------------------------------------
                                    ! Convert primitive subson3 -> conservative rightv
                                    ! ------------------------------------------------------------
                                    rightv(1) = subson3(1)
                                    rightv(2) = subson3(2)*subson3(1)
                                    rightv(3) = subson3(3)*subson3(1)
                                    rightv(4) = subson3(4)*subson3(1)

                                    skins = oo2 * ( (subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2) )
                                    ikins = subson3(5) / ( (gamma - 1.0d0) * (subson3(1) + tolsmall) )

                                    rightv(5) = (subson3(1)*ikins) + (subson3(1)*skins)





     









                end if
  end if
    
    
      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	     
		  cturbr(:)=cturbl(:)
	
      end if
    
    
    
    
    
    
    
    
    
    case(9)!outlets subsonic or supersonic will be chosen based on mach number
    ! if (boundtype.eq.0)then
     ! rightv(1:nof_variables)=leftv(1:nof_variables)

     !else

	     if ((realgas.eq.0).and.(multispecies.eq.0))then
	       call outlet_pressure_fix_ideal3d(n,leftv,rightv,nx,ny,nz,press_outlet)
	     else
       call outflow2(initcond,pox,poy,poz,rightv)
     end if



      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

		  cturbr(:)=cturbl(:)

      end if
    
    

     case(99)    !bleed boundary



     call bleed3d(iconsidered,facex,pox,poy,poz,rightv)

      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

		  cturbr(:)=cturbl(:)

      end if



    
    case(3)!symmetry
    
			      call rotatef(n,cleft_rot,leftv,angle1,angle2)
				if (rec_mrf(iconsidered).eq.1)then
                    cright_rot(1)=cleft_rot(1)
                    cright_rot(2)=-(cleft_rot(2))+2.0d0*cleft_rot(1)*srf_speedrot(2)
                    cright_rot(3)=cleft_rot(3)
                    cright_rot(4)=cleft_rot(4)
                    cright_rot(5)=cleft_rot(5)+2.0d0*cleft_rot(1)*(srf_speedrot(2)**2)-2.0d0*cleft_rot(2)*srf_speedrot(2)
                else
                    cright_rot(1)=cleft_rot(1)
                    cright_rot(2)=-cleft_rot(2)
                    cright_rot(3)=cleft_rot(3)
                    cright_rot(4)=cleft_rot(4)
                    cright_rot(5)=cleft_rot(5)

                if((multispecies.eq.1).or.(realgas.eq.1))then
                      cright_rot(6:nof_variables)=cleft_rot(6:nof_variables)
                    end if

				end if
				     if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					    cturbr(:)=cturbl(:)

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
				      end if
				
			      call rotateb(n,rightv,cright_rot,angle1,angle2)
    
    
    
    case(4)!wall
    
    
			      if (itestcase.eq.3)then
			      
			       call rotatef(n,cleft_rot,leftv,angle1,angle2)
			      if (rec_mrf(iconsidered).eq.1)then
                        cright_rot(1)=cleft_rot(1)
                        cright_rot(2)=-(cleft_rot(2))+2.0d0*cleft_rot(1)*srf_speedrot(2)
                        cright_rot(3)=cleft_rot(3)
                        cright_rot(4)=cleft_rot(4)
                        cright_rot(5)=cleft_rot(5)+cleft_rot(1)*(srf_speedrot(2)**2)*2.0d0-2.0d0*cleft_rot(2)*srf_speedrot(2)
			      else
         		      cright_rot(:)=cleft_rot(:)
			      cright_rot(2)=-cleft_rot(2)
                  end if
				     if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					    cturbr(:)=cturbl(:)

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
				      end if
				
				
				
			      call rotateb(n,rightv,cright_rot,angle1,angle2)
			      
			      
			      
			      else
                  if (rec_mrf(iconsidered).eq.1)then
                    rightv(1)=leftv(1)
                    rightv(2)=-leftv(2)+2.0d0*leftv(1)*srf_speed(2)
                    rightv(3)=-leftv(3)+2.0d0*leftv(1)*srf_speed(3)
                    rightv(4)=-leftv(4)+2.0d0*leftv(1)*srf_speed(4)
                    rightv(5)=leftv(5)+2.0d0*leftv(1)*(srf_speed(2)**2+srf_speed(3)**2+srf_speed(4)**2)&
                                            -2.0d0*(leftv(2)*srf_speed(2)+leftv(3)*srf_speed(3)+leftv(4)*srf_speed(4))
    
                  else
                    rightv(1:nof_variables)=leftv(1:nof_variables)
                    rightv(2)=-leftv(2)
                    rightv(3)=-leftv(3)
                    rightv(4)=-leftv(4)

                    if (realgas.eq.1)then
                    if (catalytic_wall.eq.1)then
                    rightv(dimensiona+4:nof_variables)=catalytic_con(1:nof_species)*leftv(1)


                    end if
                    end if








                  end if
    
    
    
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					if ((turbulence.ne.1).or.(turbulencemodel.ne.2))then
					    cturbr(:)=-cturbl(:)

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    -cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
					else
					     cturbr(1)=-cturbl(1)
					     cturbr(2)=60.0d0*visc/(beta_i1*(ielem_walldist(iconsidered)**2))

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    -cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
					  
					
					
					
					
					end if
				      end if
    
    
				end if
    
    
    
    
    
    case(6)!farfield inflow or outflow, subsonic or supersonic will be chosen based on mach number
	    call rotatef(n,cleft_rot,leftv,angle1,angle2)
	    vnb=cleft_rot(2)
	  
	    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    
	    subson1(1:nof_variables)=rightv(1:nof_variables)
	    subson2(1:nof_variables)=leftv(1:nof_variables)
	    sps=sqrt((gamma*subson2(5))/(subson2(1)))
	    vel = sqrt(subson2(2)**2+subson2(3)**2+subson2(4)**2)
	  
	    call prim2cons2(n,leftv,rightv)

	  if (vnb.le.0.0d0)then		!inflow
			ibfc=-1
	
		  if ((abs(vnb)).ge.sps)then
				!supersonic
				if ((initcond.eq.4440).and.(realgas.eq.0).and.(multispecies.eq.0))then
				call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
				else
				call inflow(initcond,pox,poy,poz,rightv)
				end if
					
		  else
				!subsonic
			
			
			if ((initcond.eq.4440).and.(realgas.eq.0).and.(multispecies.eq.0))then
			call inflow_4440_from_left_ideal3d(n,leftv,nx,ny,nz,rightv)
			else
			call inflow(initcond,pox,poy,poz,rightv)
			end if
	  	        
			  call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
			  
			
			subson1(1:nof_variables)=rightv(1:nof_variables)
			subson2(1:nof_variables)=leftv(1:nof_variables)
			sps=sqrt((gamma*subson2(5))/(subson2(1)))
			vel=sqrt(subson2(2)**2+subson2(3)**2+subson2(4)**2)
	             call prim2cons2(n,leftv,rightv)
				    
		    subson3(5)=0.5*((subson1(5))+(subson2(5))-(subson2(1)*sps*((nx*(subson1(2)-subson2(2)))+(ny*(subson1(3)-subson2(3)))&
	      +(nz*(subson1(4)-subson2(4))))))
		    subson3(1)=subson1(1)+(subson3(5)-subson1(5))/(sps**2)
		    subson3(2)=subson1(2)-(nx*(subson1(5)-subson3(5)))/(sps*subson2(1))
		    subson3(3)=subson1(3)-(ny*(subson1(5)-subson3(5)))/(sps*subson2(1))
		    subson3(4)=subson1(4)-(nz*(subson1(5)-subson3(5)))/(sps*subson2(1))
		    
		    
		      rightv(1)=subson3(1)
    rightv(2)=subson3(2)*subson3(1)
    rightv(3)=subson3(3)*subson3(1)
    rightv(4)=subson3(4)*subson3(1)
    skins=oo2*((subson3(2)**2)+(subson3(3)**2)+(subson3(4)**2))
    ikins=subson3(5)/((gamma-1.0d0)*(subson3(1)))
    rightv(5)=(subson3(1)*(ikins))+(subson3(1)*skins)
		  
		  end if
		  
		  
		  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      
      
      
      
	       if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit*rightv(1)
	      end if
	      if ((turbulence.eq.1).and.(turbulencemodel.eq.2))then	 
		cturbr(1)=(1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
		cturbr(2)=rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
		if (rec_mrf(iconsidered).eq.1)then
            cturbr(1)=(1.5d0*i_turb_inlet*(kinit_srf**2))*rightv(1)!k initialization
            cturbr(2)=rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
		end if
	      end if
 
	      if (passivescalar.gt.0)then
	      call pass_inlet(initcond,pox,poy,poz,pass_tmp)
	      do rg_i=1,passivescalar
	        cturbr(turbulenceequations+rg_i)=pass_tmp(rg_i)*rightv(1)
	      end do
	      end if
		    end if
		  
		  
      
	else
      
	  !outflow
	
	    ibfc=-2
		if ((abs(vnb)).ge.sps)then
		      
		      rightv(1:nof_variables)=leftv(1:nof_variables)
		
		else
		
		
		
		call outflow(initcond,pox,poy,poz,rightv)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
    
    subson1(1:nof_variables)=rightv(1:nof_variables)
    subson2(1:nof_variables)=leftv(1:nof_variables)
     
     
     call prim2cons2(n,leftv,rightv)
   
    
    
    subson3(5)=subson1(5)
    subson3(1)=subson2(1)+(subson3(5)-subson2(5))/(sps**2)
    subson3(2)=subson2(2)+(nx*(subson2(5)-subson3(5)))/(sps*subson2(1))
    subson3(3)=subson2(3)+(ny*(subson2(5)-subson3(5)))/(sps*subson2(1))
    subson3(4)=subson2(4)+(nz*(subson2(5)-subson3(5)))/(sps*subson2(1))
! 							
    rightv(1)=subson3(1)
    rightv(2)=subson3(2)*subson3(1)
    rightv(3)=subson3(3)*subson3(1)
    rightv(4)=subson3(4)*subson3(1)
    skins=oo2*((subson3(2)**2)+(subson3(3)**2)+(subson3(4)**2))
    ikins=subson3(5)/((gamma-1.0d0)*(subson3(1)))
    rightv(5)=(subson3(1)*(ikins))+(subson3(1)*skins)
	
	      
	      
	       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	     
		  cturbr(:)=cturbl(:)
	
		end if
	      end if
	      

	end if

	
! 	call rotatef(n,cright_rot,rightv,angle1,angle2)
    
    

			      


end select



end subroutine boundarys


subroutine boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
implicit none
!> @brief
!> this subroutine applies the boundary condition to each bounded cell
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,b_code,iconsidered,facex
integer,intent(inout)::ibfc
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,dimension(1:nof_variables),intent(in)::srf_speedrot,srf_speed
real,dimension(1:dimensiona),intent(in)::pox,poy,poz
real,intent(in)::angle1,angle2,nx,ny,nz
real,dimension(1:turbulenceequations+passivescalar),intent(inout)::cturbl,cturbr
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::cright_rot,cleft_rot
real,dimension(1:gpu_max_nvar)::subson1,subson2,subson3
real::mp_pinfl,mp_pinfr,gammal,gammar,u,v,sa_bc
real::sps,skins,ikins,vel,vnb,theeta,reeta
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx,vhx,amp,dvel,rgg,tt1
real::mn,un,agrt,rho_g,rg_rmix,rmix,t_g,tv_g,rg_tr,rg_ev_total,rg_chem
real,dimension(1:gpu_max_passive)::pass_tmp
integer::rg_i





select case(b_code)

	    
	    case(1)!inflow subsonic or supersonic will be chosen based on mach number
	    if ((initcond.eq.4440).and.(realgas.eq.0).and.(multispecies.eq.0))then
	    call inflow_4440_from_left_ideal2d(n,leftv,nx,ny,rightv)
	    if ((turbulence.eq.1).or.(passivescalar.gt.0)) cturbr(:)=cturbl(:)
	    if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
	      call inflow_4440_sa_value(n,rightv,sa_bc)
	      cturbr(1)=sa_bc
	    end if
	    else
	    if (boundtype.eq.0)then	!supersonic
    
    call inflow2d(initcond,pox,poy,rightv)
    
     if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit*rightv(1)
	      end if
    
        
    
    
    else		!subsonic
    call inflow2d(initcond,pox,poy,rightv)

    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    
    subson1(1:nof_variables)=rightv(1:nof_variables)
    subson2(1:nof_variables)=leftv(1:nof_variables)
    sps=sqrt((gamma*subson2(4))/(subson2(1)))
    vel=sqrt(subson2(2)**2+subson2(3)**2)
    

    call prim2cons2(n,leftv,rightv)
    

      if (vel/(sps+tolsmall).gt.1.0d0)then	!supersonic
      
      
      call inflow2d(initcond,pox,poy,rightv)
      
      
      
      else		!subsonic
      
      
      subson3(4)=0.5*((subson1(4))+(subson2(4))-(subson2(1)*sps*((nx*(subson1(2)-subson2(2)))+(ny*(subson1(3)-subson2(3)))&
)))
      subson3(1)=subson1(1)+(subson3(4)-subson1(4))/(sps**2)
      subson3(2)=subson1(2)-(nx*(subson1(4)-subson3(4)))/(sps*subson2(1))
      subson3(3)=subson1(3)-(ny*(subson1(4)-subson3(4)))/(sps*subson2(1))
      
      
      		    
      
      
    rightv(1)=subson3(1)
    rightv(2)=subson3(2)*subson3(1)
    rightv(3)=subson3(3)*subson3(1)
    
    skins=oo2*((subson3(2)**2)+(subson3(3)**2))
    ikins=subson3(4)/((gamma-1.0d0)*(subson3(1)))
    rightv(4)=(subson3(1)*(ikins))+(subson3(1)*skins)
      
      
       
      
      end if

      
       end if
      
      
      
      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      
      
      
      
	      if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit
	      end if
	     if ((turbulence.eq.1).and.(turbulencemodel.eq.2))then	 
		cturbr(1)=(1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
		cturbr(2)=rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
		
		
	      end if
 
	      if (passivescalar.gt.0)then
	      call pass_inlet2d(initcond,pox,poy,pass_tmp)
	      do rg_i=1,passivescalar
	        cturbr(turbulenceequations+rg_i)=pass_tmp(rg_i)*rightv(1)
	      end do
	      end if
      end if
      
      
      
      
      
      
      
      
    
	   
	    
	    
	    end if
	    case(2)!outflow subsonic or supersonic will be chosen based on mach number
     if (boundtype.eq.0)then
      rightv(1:nof_variables)=leftv(1:nof_variables)
      
     else
     




                if (realgas.eq.1)then
                  call boundarys2d_realgas_outflow(n,leftv,rightv,nx,ny)
                  else
                        call outflow2d(initcond,pox,poy,rightv)
                        call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

                        subson1(1:nof_variables)=rightv(1:nof_variables)
                        subson2(1:nof_variables)=leftv(1:nof_variables)

                        sps=sqrt((gamma*subson2(4))/(subson2(1)))
                        vel=subson2(2)*nx + subson2(3)*ny
!                         sqrt(subson2(2)**2+subson2(3)**2)

                        sps=sqrt((gamma*subson2(4))/(subson2(1)))

                        call prim2cons2(n,leftv,rightv)

                        if (vel/(sps+tolsmall).gt.1.0d0)then	!supersonic
                        rightv(1:nof_variables)=leftv(1:nof_variables)


                          else
                        subson3(4)=subson1(4)
                        subson3(1)=subson2(1)+(subson3(4)-subson2(4))/(sps**2)
                        subson3(2)=subson2(2)+(nx*(subson2(4)-subson3(4)))/(sps*subson2(1))
                        subson3(3)=subson2(3)+(ny*(subson2(4)-subson3(4)))/(sps*subson2(1))

                    !
                        rightv(1)=subson3(1)
                        rightv(2)=subson3(2)*subson3(1)
                        rightv(3)=subson3(3)*subson3(1)

                        skins=oo2*((subson3(2)**2)+(subson3(3)**2))
                        ikins=subson3(4)/((gamma-1.0d0)*(subson3(1)))
                        rightv(4)=(subson3(1)*(ikins))+(subson3(1)*skins)


                        end if

                  end if
    end if

    
    
      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	     
		  cturbr(:)=cturbl(:)
	
      end if
    
    
    
    
    
    
     case(99)  !bleed


     call bleed2d(iconsidered,facex,pox,poy,rightv)



      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

		  cturbr(:)=cturbl(:)

      end if
    
    
    
    
    
    
    case(3)!symmetry
    
    
    
                            if ((initcond.eq.102).or.(initcond.eq.30).or.(initcond.eq.222))then	!shock density interaction
                                   if ((initcond.eq.102))then	!shock density interaction
                            if (pox(1).lt.((1.0d0/6.0d0)+((1.0d0+20.0d0*t)/(sqrt(3.0d0)))))then
                            r1=8.0d0
                            u1=8.25*cos(pi/6.0d0)
                            v1=-8.25*sin(pi/6.0d0)
                            p1=116.5
                            else
                            r1=1.4d0
                            u1=zero
                            v1=zero
                            p1=1.0d0
                            end if
                            skin1=(oo2)*((u1**2)+(v1**2))
                            !internal energy 
                            ie1=((p1)/((gamma-1.0d0)*r1))
                            !total energy
                            e1=(p1/(gamma-1))+(r1*skin1)
                            !vector of conserved variables now
                            rightv(1)=r1
                            rightv(2)=r1*u1
                            rightv(3)=r1*v1
                            rightv(4)=e1
                            end if
    
                            
                            
                            
                            
                            
                             if ((initcond.eq.222))then
                             if (sqrt((pox(1)**2)+(poy(1)**2)).lt.(t/3.0d0))then
                            r1=16.0d0
                            u1=0.0
                            v1=0.0
                            p1=16.0d0/3.0d0
                            else
                            r1=1.0d0+(t/sqrt((pox(1)**2)+(poy(1)**2)))
                            reeta=-1
                            theeta=atan(poy(1)/pox(1))
                            u1=reeta*cos(theeta)
                            v1=reeta*sin(theeta)
                            p1=1.0e-6
                            end if
                            skin1=(oo2)*((u1**2)+(v1**2))
                            !internal energy 
                            ie1=((p1)/((gamma-1.0d0)*r1))
                            !total energy
                            e1=(p1/(gamma-1))+(r1*skin1)
                            !vector of conserved variables now
                            rightv(1)=r1
                            rightv(2)=r1*u1
                            rightv(3)=r1*v1
                            rightv(4)=e1
                             
                             
                             
                             end if
			     
			       if ((initcond.eq.30))then
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
                            ie1=((p1)/((gamma-1.0d0)*r1))
                            !total energy
                            e1=(p1/(gamma-1))+(r1*skin1)
                            !vector of conserved variables now
                           ! rightv(1)=leftv(1)
                            !rightv(2)=0.0
                            !rightv(3)=0.0
                            !rightv(4)=leftv(4)
                            
                             rightv(1)=r1
                            rightv(2)=r1*u1
                            rightv(3)=r1*v1
                            rightv(4)=e1
			      
			      end if
			      
			      
			      
			      
			      else
			      
			      
			       call rotatef2d(n,cleft_rot,leftv,angle1,angle2)
			      
                  cright_rot(1)=cleft_rot(1)
			      cright_rot(2)=-cleft_rot(2)
			      cright_rot(3)=cleft_rot(3)
			      cright_rot(4)=cleft_rot(4)
			      
			       if((multispecies.eq.1).or.(realgas.eq.1))then
                      cright_rot(5:nof_variables)=cleft_rot(5:nof_variables)
                    
                    end if
			     
					 
				     if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					    cturbr(:)=cturbl(:)

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
				      end if
				
			    	
				
			      call rotateb2d(n,rightv,cright_rot,angle1,angle2)
                            end if
    
    
    case(4)!wall
    
			     if (itestcase.eq.3)then
			      
			       call rotatef2d(n,cleft_rot,leftv,angle1,angle2)
			      
			      
                      if ((multispecies.eq.1).or.(realgas.eq.1))then
                          cright_rot(:)=cleft_rot(:)
                      cright_rot(2)=-cleft_rot(2)





                      else
                        cright_rot(:)=cleft_rot(:)
                        cright_rot(1)=cleft_rot(1)
                        cright_rot(2)=-cleft_rot(2)
                        cright_rot(3)=cleft_rot(3)
                        cright_rot(4)=cleft_rot(4)
                      end if
			     
			      
					 
				     if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					    cturbr(:)=cturbl(:)

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
				      end if
				
				
				
			      call rotateb2d(n,rightv,cright_rot,angle1,angle2)
			      
			    
			      
			      else


                  rightv(1:nof_variables)=leftv(1:nof_variables)


			      rightv(2)=-leftv(2)
			      rightv(3)=-leftv(3)
			      
                if (realgas.eq.1)then
                if (catalytic_wall.eq.1)then
                rightv(dimensiona+4:nof_variables)=catalytic_con(1:nof_species)*leftv(1)


                end if
                end if
    

    
    
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					if ((turbulence.ne.1).or.(turbulencemodel.ne.2))then
					    cturbr(:)=-cturbl(:)

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    -cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
					else
					     cturbr(1)=-cturbl(1)
					     cturbr(2)=60.0d0*visc/(beta_i1*(ielem_walldist(iconsidered)**2))

					  if (passivescalar.gt.0)then
					  cturbr(turbulenceequations+1:turbulenceequations+passivescalar)=&
						    -cturbl(turbulenceequations+1:turbulenceequations+passivescalar)

					  end if
					  
					
					
					
					
					end if
				      end if
    
    
				end if
    
    
    
    
    
    
    
    
    
    case(6)!farfield inflow or outflow, subsonic or supersonic will be chosen based on mach number
	 call rotatef2d(n,cleft_rot,leftv,angle1,angle2)
	    vnb=cleft_rot(2)/cleft_rot(1)
	    
	    
	  
	    
	    call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    
	    subson1(1:nof_variables)=rightv(1:nof_variables)
	    subson2(1:nof_variables)=leftv(1:nof_variables)
	    sps=sqrt((gamma*subson2(4))/(subson2(1)))
	    vel=sqrt(subson2(2)**2+subson2(3)**2)
	  
	  call prim2cons2(n,leftv,rightv)

	  if (vnb.le.0.0d0)then		!inflow
			ibfc=-1
	
		  if ((abs(vnb)).ge.sps)then
				!supersonic
				if ((initcond.eq.4440).and.(realgas.eq.0).and.(multispecies.eq.0))then
				call inflow_4440_from_left_ideal2d(n,leftv,nx,ny,rightv)
				else
				call inflow2d(initcond,pox,poy,rightv)
				end if
					
		  else
				!subsonic
				
			if ((initcond.eq.4440).and.(realgas.eq.0).and.(multispecies.eq.0))then
			call inflow_4440_from_left_ideal2d(n,leftv,nx,ny,rightv)
			else
			call inflow2d(initcond,pox,poy,rightv)
			end if
	  	  
			  
			  call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
			
			subson1(1:nof_variables)=rightv(1:nof_variables)
			subson2(1:nof_variables)=leftv(1:nof_variables)
			sps=sqrt((gamma*subson2(4))/(subson2(1)))
			vel=sqrt(subson2(2)**2+subson2(3)**2)
			  call prim2cons2(n,leftv,rightv)
				    
		    subson3(4)=0.5d0*((subson1(4))+(subson2(4))-(subson2(1)*sps*((nx*(subson1(2)-subson2(2)))+(ny*(subson1(3)-subson2(3))))))
		    subson3(1)=subson1(1)+(subson3(4)-subson1(4))/(sps**2)
		    subson3(2)=subson1(2)-(nx*(subson1(4)-subson3(4)))/(sps*subson2(1))
		    subson3(3)=subson1(3)-(ny*(subson1(4)-subson3(4)))/(sps*subson2(1))
		    
		    

		    rightv(1)=subson3(1)
    rightv(2)=subson3(2)*subson3(1)
    rightv(3)=subson3(3)*subson3(1)
    
    skins=oo2*((subson3(2)**2)+(subson3(3)**2))
    ikins=subson3(4)/((gamma-1.0d0)*(subson3(1)))
    rightv(4)=(subson3(1)*(ikins))+(subson3(1)*skins)
		  
		  end if
		  
		  
		  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      
      
      
      
	        if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
		  cturbr(1)=visc*turbinit*rightv(1)
	      end if
	     if ((turbulence.eq.1).and.(turbulencemodel.eq.2))then	 
		cturbr(1)=(1.5d0*i_turb_inlet*(ufreestream**2))*rightv(1)!k initialization
		cturbr(2)=rightv(1)*cturbr(1)/(10.0e-5*visc)!omega initialization
	      end if
 
	      if (passivescalar.gt.0)then
	      call pass_inlet2d(initcond,pox,poy,pass_tmp)
	      do rg_i=1,passivescalar
	        cturbr(turbulenceequations+rg_i)=pass_tmp(rg_i)*rightv(1)
	      end do
	      end if
		    end if
		  
		  
      
	else
      
	  !outflow
	
	    ibfc=-2
		if ((abs(vnb)).ge.sps)then
		
		      rightv(1:nof_variables)=leftv(1:nof_variables)
		
		else
		
		
		call outflow2d(initcond,pox,poy,rightv)
		 call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
    
    subson1(1:nof_variables)=rightv(1:nof_variables)
    subson2(1:nof_variables)=leftv(1:nof_variables)
     
     call prim2cons2(n,leftv,rightv)
     
   
    
    
    subson3(4)=subson1(4)
    subson3(1)=subson2(1)+(subson3(4)-subson2(4))/(sps**2)
    subson3(2)=subson2(2)+(nx*(subson2(4)-subson3(4)))/(sps*subson2(1))
    subson3(3)=subson2(3)+(ny*(subson2(4)-subson3(4)))/(sps*subson2(1))
    
! 							
    rightv(1)=subson3(1)
    rightv(2)=subson3(2)*subson3(1)
    rightv(3)=subson3(3)*subson3(1)
    
    skins=oo2*((subson3(2)**2)+(subson3(3)**2))
    ikins=subson3(4)/((gamma-1.0d0)*(subson3(1)))
    rightv(4)=(subson3(1)*(ikins))+(subson3(1)*skins)
	
	      
	      
	       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	     
		  cturbr(:)=cturbl(:)
	
		end if
	        end if
	      

		end if

		

	    
	    
	
				      
	
	    case(9)
	    if ((realgas.eq.0).and.(multispecies.eq.0))then
	      call outlet_pressure_fix_ideal2d(n,leftv,rightv,nx,ny,press_outlet)
	    else
	      call outflow2d(initcond,pox,poy,rightv)
	    end if
	    if ((turbulence.eq.1).or.(passivescalar.gt.0))then
	      cturbr(:)=cturbl(:)
	    end if
	

	end select



end subroutine boundarys2d


subroutine compute_eigenvectors(n,rveigl,rveigr,eigvl,eigvr,gamma)
implicit none
!> @brief
!> this subroutine computes the left and right eigenvectors 
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::rveigl,rveigr
real,intent(in)::gamma
real,dimension(1:nof_variables,1:nof_variables),intent(inout)::eigvl,eigvr
real::rs,us,vs,ws,es,ps,vvs,as,hs,gammam1,vsd,oor1,oor2
integer::ivgt

eigvr=zero
gammam1=gamma-1.0d0
oor1=1.0d0/rveigl(1)
oor2=1.0d0/rveigr(1)

rs=oo2*(rveigl(1)+rveigr(1))
us=oo2*((rveigl(2)*oor1)+(rveigr(2)*oor2))
vs=oo2*((rveigl(3)*oor1)+(rveigr(3)*oor2))
ws=oo2*((rveigl(4)*oor1)+(rveigr(4)*oor2))
es=oo2*(rveigl(5)+rveigr(5))
 
vvs=(us**2)+(vs**2)+(ws**2)
vsd=oo2*vvs
ps=(gamma-1.0d0)*(es - oo2*rs*vvs)
as=sqrt(gamma*ps/rs)
hs=(oo2*vvs) + ((as**2)/gammam1)
eigvr(1,1)=1.0d0		; eigvr(1,2)=1.0d0	; eigvr(1,3)=0.0d0	; eigvr(1,4)=0.0d0	; eigvr(1,5)=1.0d0
eigvr(2,1)=us-as	; eigvr(2,2)=us		; eigvr(2,3)=0.0d0	; eigvr(2,4)=0.0d0	; eigvr(2,5)=us+as
eigvr(3,1)=vs		; eigvr(3,2)=vs		; eigvr(3,3)=1.0d0	; eigvr(3,4)=0.0d0	; eigvr(3,5)=vs
eigvr(4,1)=ws		; eigvr(4,2)=ws		; eigvr(4,3)=0.0d0	; eigvr(4,4)=1.0d0	; eigvr(4,5)=ws
eigvr(5,1)=hs-(us*as)	; eigvr(5,2)=oo2*vvs	; eigvr(5,3)=vs		; eigvr(5,4)=ws		; eigvr(5,5)=hs+(us*as)

eigvl(1,1)=(hs*us + as*us**2 + as*vs**2 - as*vsd - us*vsd + as*ws**2)/ (2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(1,2)=(-hs - as*us + vsd)/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(1,3)= -((as*vs)/(2.0d0*as*hs - 2.0d0*as*vsd))
eigvl(1,4)=  -((as*ws)/(2.0d0*as*hs - 2.0d0*as*vsd))
eigvl(1,5)=as/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(2,1)=(2.0*as*hs - 2.0*as*us**2 - 2.0d0*as*vs**2 -2.0d0*as*ws**2)/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(2,2)=(2.0d0*as*us)/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(2,3)=(2.0d0*as*vs)/(2.0d0*as*hs - 2.0d0*as*vsd)
 eigvl(2,4)=(2.0d0*as*ws)/(2.0d0*as*hs - 2.0d0*as*vsd) 
 eigvl(2,5)=(-2.0d0*as)/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(3,1)=(-2.0*as*hs*vs + 2.0*as*vs*vsd)/(2.0*as*hs - 2.0*as*vsd)
eigvl(3,2)=0.0d0
eigvl(3,3)=1.0d0
eigvl(3,4)=0.0d0
 eigvl(3,5)=0.0d0
eigvl(4,1)=(-2.0d0*as*hs*ws + 2.0d0*as*vsd*ws)/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(4,2)=0.0d0
 eigvl(4,3)=0.0d0
eigvl(4,4)=1.0d0
eigvl(4,5)=0.0d0
eigvl(5,1)=(-(hs*us) + as*us**2 + as*vs**2 - as*vsd + us*vsd + as*ws**2)/ (2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(5,2)=(hs - as*us - vsd)/(2.0d0*as*hs - 2.0d0*as*vsd)
eigvl(5,3)=-((as*vs)/(2.0d0*as*hs - 2.0d0*as*vsd))
 eigvl(5,4)=-((as*ws)/(2.0d0*as*hs - 2.0d0*as*vsd))
eigvl(5,5)= as/(2.0d0*as*hs - 2.0d0*as*vsd)
 
 

end subroutine compute_eigenvectors








subroutine compute_jacobianse(n,iconsidered,eigvl,rveigl,gamma,angle1,angle2,srf_speedrot,nx,ny,nz)
implicit none
!> @brief
!> this subroutine computes the jacobians for the implicit time stepping
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,intent(in)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables),intent(in)::rveigl,srf_speedrot
real,intent(in)::gamma
real,dimension(1:nof_variables,1:nof_variables),intent(inout)::eigvl
real::rs,us,vs,ws,es,ps,vvs,as,hs,gammam1,vsd,phi,a1,a2,a3,oors,vn
integer::ivgt,i
real :: evs, etr
real :: temps,gammal,mp_pinfl
integer :: rg_i
real,dimension(1:gpu_max_nvar)::leftv

if (realgas.eq.0)then


a2=gamma-1.0d0
a3=gamma-2.0d0
oors=1.0d0/rveigl(1)
rs=(rveigl(1))
us=(rveigl(2)*oors)
vs=(rveigl(3)*oors)
ws=(rveigl(4)*oors)
es=(rveigl(5)*oors)
phi=oo2*(a2)*((us*us)+(vs*vs)+(ws*ws))
 a1=gamma*es-phi

 
 
vvs=nx*us+ny*vs+nz*ws

if (rec_mrf(iconsidered).eq.1)then
    eigvl(1,1)=0.0d0-srf_speedrot(2);eigvl(1,2)=nx	; 		eigvl(1,3)=ny	; 		eigvl(1,4)=nz	; 		eigvl(1,5)=0.0d0
    eigvl(2,1)=nx*phi-us*vvs	;eigvl(2,2)=vvs-a3*nx*us-srf_speedrot(2)	;eigvl(2,3)=ny*us-a2*nx*vs	; eigvl(2,4)=nz*us-a2*nx*ws; eigvl(2,5)=a2*nx
    eigvl(3,1)=ny*phi-vs*vvs	;eigvl(3,2)=nx*vs-a2*ny*us	; eigvl(3,3)=vvs-a3*ny*vs-srf_speedrot(2)	; eigvl(3,4)=nz*vs-a2*ny*ws; eigvl(3,5)=a2*ny
    eigvl(4,1)=nz*phi-ws*vvs	;eigvl(4,2)=nx*ws-a2*nz*us	; eigvl(4,3)=ny*ws-a2*nz*vs	; eigvl(4,4)=vvs-a3*nz*ws-srf_speedrot(2); eigvl(4,5)=a2*nz
    eigvl(5,1)=vvs*(phi-a1)         ; eigvl(5,2)=nx*a1-a2*us*vvs	; eigvl(5,3)=ny*a1-a3*vs*vvs; eigvl(5,4)=nz*a1-a2*ws*vvs; eigvl(5,5)=gamma*vvs-srf_speedrot(2)
else
    eigvl(1,1)=0.0d0		; 		eigvl(1,2)=nx	; 		eigvl(1,3)=ny	; 		eigvl(1,4)=nz	; 		eigvl(1,5)=0.0d0
    eigvl(2,1)=nx*phi-us*vvs	; 	eigvl(2,2)=vvs-a3*nx*us	; 	eigvl(2,3)=ny*us-a2*nx*vs	; eigvl(2,4)=nz*us-a2*nx*ws	; eigvl(2,5)=a2*nx
    eigvl(3,1)=ny*phi-vs*vvs		; eigvl(3,2)=nx*vs-a2*ny*us	; eigvl(3,3)=vvs-a3*ny*vs	; eigvl(3,4)=nz*vs-a2*ny*ws	; eigvl(3,5)=a2*ny
    eigvl(4,1)=nz*phi-ws*vvs		; eigvl(4,2)=nx*ws-a2*nz*us	; eigvl(4,3)=ny*ws-a2*nz*vs	; eigvl(4,4)=vvs-a3*nz*ws	; eigvl(4,5)=a2*nz
    eigvl(5,1)=vvs*(phi-a1)	;		 eigvl(5,2)=nx*a1-a2*us*vvs	; eigvl(5,3)=ny*a1-a2*vs*vvs	; eigvl(5,4)=nz*a1-a2*ws*vvs	; eigvl(5,5)=gamma*vvs
end if

 


 end if

 if (realgas.eq.1)then

 leftv = rveigl
call cons2prim(n, leftv, mp_pinfl, gammal)

! === extract =================================================
rs  = rveigl(1)
us  = rveigl(2) / rs
vs  = rveigl(3) / rs
ws  = rveigl(4) / rs
es  = rveigl(5) / rs
evs = rveigl(6) / rs

! effective translational energy for pressure
etr = es - evs

a2 = gammal - 1.0d0
a3 = gammal - 2.0d0

phi = 0.5d0*a2*(us*us + vs*vs + ws*ws)
a1  = gammal*etr - phi

vn = nx*us + ny*vs + nz*ws

eigvl(:,:) = 0.0d0

! === 5x5 euler block (ρ, ρu, ρv, ρw, ρe) =====================
eigvl(1,1)=0.0d0;   eigvl(1,2)=nx;   eigvl(1,3)=ny;   eigvl(1,4)=nz;   eigvl(1,5)=0.0d0

eigvl(2,1)=nx*phi - us*vn
eigvl(2,2)=vn - a3*nx*us
eigvl(2,3)=ny*us - a2*nx*vs
eigvl(2,4)=nz*us - a2*nx*ws
eigvl(2,5)=a2*nx

eigvl(3,1)=ny*phi - vs*vn
eigvl(3,2)=nx*vs - a2*ny*us
eigvl(3,3)=vn - a3*ny*vs
eigvl(3,4)=nz*vs - a2*ny*ws
eigvl(3,5)=a2*ny

eigvl(4,1)=nz*phi - ws*vn
eigvl(4,2)=nx*ws - a2*nz*us
eigvl(4,3)=ny*ws - a2*nz*vs
eigvl(4,4)=vn - a3*nz*ws
eigvl(4,5)=a2*nz

eigvl(5,1)=vn*(phi - a1)
eigvl(5,2)=nx*a1 - a2*us*vn
eigvl(5,3)=ny*a1 - a2*vs*vn
eigvl(5,4)=nz*a1 - a2*ws*vn
eigvl(5,5)=gammal*vn

! === vibrational scalar advection (ρev) ======================
eigvl(6,6) = vn

! === species scalar advection ================================
do i = 7, nof_variables
  eigvl(i,i) = vn
end do




 end if




 

end subroutine compute_jacobianse




subroutine compute_eigenvectors2d(n,rveigl,rveigr,eigvl,eigvr,gamma)
implicit none
!> @brief
!> this subroutine computes the left and right eigenvectors  in 2d
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,dimension(1:nof_variables),intent(in)::rveigl,rveigr
real,intent(in)::gamma
real,dimension(1:nof_variables,1:nof_variables),intent(inout)::eigvl,eigvr
real::rs,us,vs,ws,es,ps,vvs,as,hs,gammam1,vsd,g8,s1,s2,vtots,oor1,oor2
integer::ivgt,j,k

eigvr=zero
gammam1=gamma-1.0d0
oor1=1.0d0/rveigl(1)
oor2=1.0d0/rveigr(1)


eigvr=0.0d0
gammam1=gamma-1.0d0
rs=0.5*(rveigl(1)+rveigr(1))
us=0.5*((rveigl(2)*oor1)+(rveigr(2)*oor2))
vs=0.5*((rveigl(3)*oor1)+(rveigr(3)*oor2))
es=0.5*(rveigl(4)+rveigr(4))
g8 = gamma - 1.0d0
 
vtots = us**2 + vs**2
   ps = (gamma-1)*(es - oo2*rs*vtots)
   as = sqrt(gamma*ps/rs)
   hs = oo2*vtots + (as**2)/(g8)
  

   eigvr(1,1) =1.d0;        eigvr(1,2) = 1.d0;                eigvr(1,3) =0.d0;  eigvr(1,4) = 1.d0  
   eigvr(2,1) =us-as;     eigvr(2,2) = us;                eigvr(2,3) =0.d0;  eigvr(2,4) = us+as 
   eigvr(3,1) =vs;        eigvr(3,2) = vs;                eigvr(3,3) =1.d0;  eigvr(3,4) = vs   
   eigvr(4,1) =hs-us*as;  eigvr(4,2) = oo2*vtots; eigvr(4,3) =vs ; eigvr(4,4) = hs+us*as 

   s1 = as/(gamma-1d0)
   s2 = as**2/(gamma-1d0)

   eigvl(1,1) =   hs + s1*(us-as);   eigvl(1,2) = -(us+s1); eigvl(1,3) = -vs;            eigvl(1,4) = 1.d0
   eigvl(2,1) =-2d0*hs+4d0*s2;           eigvl(2,2) = 2d0*us;     eigvl(2,3) = 2d0*vs;           eigvl(2,4) = -2.d0
   eigvl(3,1) =-2d0*vs*s2;             eigvl(3,2) =   0.d0;     eigvl(3,3) = 2d0*s2;           eigvl(3,4) = 0.d0
   eigvl(4,1) = hs- s1*(us+as);      eigvl(4,2) = -us+s1;   eigvl(4,3) = - vs;           eigvl(4,4) = 1.d0 

do j=1,4
	do k=1,4

   eigvl(j,k) = eigvl(j,k)/(2d0*s2)
	end do
end do
 
 

end subroutine compute_eigenvectors2d


 subroutine compute_jacobianse2d(n,eigvl,rveigl,gamma,angle1,angle2,nx,ny,nz)
implicit none
 !> @brief
!> this subroutine computes the jacobians for the implicit time stepping in 2d
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n
real,intent(in)::angle1,angle2,nx,ny,nz
real,dimension(1:nof_variables),intent(in)::rveigl
real,intent(in)::gamma
real,dimension(1:nof_variables,1:nof_variables),intent(inout)::eigvl
real::rs,us,vs,es,ps,vvs,as,hs,gammam1,vsd,phi,a1,a2,a3,oors
integer::ivgt
real :: evs, etr,yk
real :: temps,gammal,mp_pinfl
integer :: i,rg_i,row,k
real,dimension(1:gpu_max_nvar)::leftv

if (realgas.eq.0)then
a2=gamma-1.0d0
a3=gamma-2.0d0
oors=1.0d0/rveigl(1)
rs=(rveigl(1))
us=(rveigl(2)*oors)
vs=(rveigl(3)*oors)
es=(rveigl(4)*oors)
phi=oo2*(a2)*((us*us)+(vs*vs))
 a1=gamma*es-phi

 
vvs=nx*us+ny*vs


eigvl(1,1)=0.0d0		; 		eigvl(1,2)=nx	; 		eigvl(1,3)=ny	; 		eigvl(1,4)=0.0d0
eigvl(2,1)=nx*phi-us*vvs	; 	eigvl(2,2)=vvs-a3*nx*us	; 	eigvl(2,3)=ny*us-a2*nx*vs	; eigvl(2,4)=a2*nx
eigvl(3,1)=ny*phi-vs*vvs		; eigvl(3,2)=nx*vs-a2*ny*us	; eigvl(3,3)=vvs-a3*ny*vs	; eigvl(3,4)=a2*ny

eigvl(4,1)=vvs*(phi-a1)	;		 eigvl(4,2)=nx*a1-a2*us*vvs	; eigvl(4,3)=ny*a1-a2*vs*vvs	; eigvl(4,4)=gamma*vvs



end if

if (realgas.eq.1)then
leftv=rveigl

call cons2prim(n,leftv,mp_pinfl,gammal)

! === extract ===============================================


rs   = rveigl(1)
us   = rveigl(2)/rs
vs   = rveigl(3)/rs
es   = rveigl(4)/rs
evs  = rveigl(5)/rs

! effective translational energy for pressure
etr = es - evs   ! (chemical echem handled in flux elsewhere)

a2 = gammal - 1.0d0
a3 = gammal - 2.0d0

phi = 0.5d0*a2*(us*us + vs*vs)
a1  = gammal*etr - phi

vvs = nx*us + ny*vs

eigvl(:,:) = 0.0d0

! === 4x4 euler block =======================================
eigvl(1,1)=0.0d0;      eigvl(1,2)=nx;       eigvl(1,3)=ny;        eigvl(1,4)=0.0d0

eigvl(2,1)=nx*phi-us*vvs
eigvl(2,2)=vvs-a3*nx*us
eigvl(2,3)=ny*us-a2*nx*vs
eigvl(2,4)=a2*nx

eigvl(3,1)=ny*phi-vs*vvs
eigvl(3,2)=nx*vs-a2*ny*us
eigvl(3,3)=vvs-a3*ny*vs
eigvl(3,4)=a2*ny

eigvl(4,1)=vvs*(phi-a1)
eigvl(4,2)=nx*a1-a2*us*vvs
eigvl(4,3)=ny*a1-a2*vs*vvs
eigvl(4,4)=gammal*vvs


! === vibrational scalar advection (ρev) ======================
eigvl(5,5) = vvs

! === species scalar advection ================================
do i = 6, nof_variables
  eigvl(i,i) = vvs
end do





end if







!
! ! === vibrational scalar advection ===========================
! eigvl(5,5) = vvs
!
! ! === species scalar advection ===============================
! do i=6,nof_variables
!   eigvl(i,i) = vvs
! end do
!
! end if
 

end subroutine compute_jacobianse2d









subroutine eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
	implicit none
!> @brief
!> compatibility wrapper for the ideal-gas eddy-viscosity implementation
#ifdef gpu
!$omp declare target
#endif
    real,dimension(1:2),intent(inout)::turbmv
    real,dimension(1:nof_variables),intent(in)::leftv,rightv
    real,dimension(1),intent(inout)::etvm
    real,dimension(1:4),intent(inout)::viscl,laml
    real,dimension(1:20),intent(inout)::eddyfl,eddyfr
	integer,intent(in)::n

	call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)

end subroutine eddyvisco


subroutine eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
	implicit none
!> @brief
!> this subroutine computes the tubulent eddy viscosity for turbulence models
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
    real,dimension(1:2),intent(inout)::turbmv
    real,dimension(1:nof_variables),intent(in)::leftv,rightv
    real,dimension(1),intent(inout)::etvm
    real,dimension(1:4),intent(inout)::viscl,laml
    real,dimension(1:20),intent(inout)::eddyfl,eddyfr
	integer,intent(in)::n
	real::ml,smmm,mtt,chi,tolepsma,chipow3,fv1
	integer:: ihgt, ihgj
	real,dimension(3,3)::vortet,tvort,svort
	real:: ux,uy,uz,vx,vy,vz,wx,wy,wz
		real:: wally, d_omplus, phi_2, phi_1, f_1, f_2, k_0, om_0,alpha_inf, alpha_star, alpha_star0_l, mu_turb
	real:: sigma_k_l, sigma_om_l, sigma_k_r, sigma_om_r
	real::snorm,dervk_dervom,re_t_sst,rho_0,beta_i




	tolepsma = 10e-16
	
	if (turbulence.eq.0)then
	viscl(4)=zero;viscl(3)=zero
	laml(4)=zero;laml(3)=zero
	else
	
	
	
!modified on 19/6/2013
  select case(turbulencemodel)
  
   case(1)
	  turbmv(1)=eddyfl(2)
	  turbmv(2)=eddyfr(2)  
         chi     = abs ( max(turbmv(1),tolepsma) / max(viscl(1),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
	viscl(3) = turbmv(1)*fv1
	chi     = abs ( max(turbmv(2),tolepsma) / max(viscl(2),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
		viscl(4) = turbmv(2)*fv1


  case(2)

  
  !for left cell
  
 vortet(1,1:3) = eddyfl(4:6)
 vortet(2,1:3) = eddyfl(7:9) 
 vortet(3,1:3) = eddyfl(10:12)

  ux = vortet(1,1);uy = vortet(1,2);uz = vortet(1,3)
  vx = vortet(2,1);vy = vortet(2,2);vz = vortet(2,3)
  wx = vortet(3,1);wy = vortet(3,2);wz = vortet(3,3)

  do ihgt=1,3
  do ihgj=1,3
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
  end do

svort=0.5*(vortet+tvort)
snorm=sqrt(2.0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+(svort(1,3)*svort(1,3))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))+(svort(2,3)*svort(2,3))+& 
	       (svort(3,1)*svort(3,1))+(svort(3,2)*svort(3,2))+(svort(3,3)*svort(3,3))))
		   
 wally=eddyfl(1)
 rho_0=leftv(1)
 k_0=max(tolepsma,eddyfl(2)/leftv(1))
 om_0=max(eddyfl(3)/leftv(1),ufreestream/charlength/10.0)


 dervk_dervom=(eddyfl(13)*eddyfl(16))+(eddyfl(14)*eddyfl(17))+(eddyfl(15)*eddyfl(18))
		   
d_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(1)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally)) 

f_1=tanh(phi_1**4)
f_2=tanh(phi_2**2)
re_t_sst=rho_0*k_0/(viscl(1)*om_0)
beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
alpha_star0_l=beta_i/3.0    
alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0_l+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)


viscl(3)=rho_0*k_0/om_0/max(1.0/alpha_star,snorm*f_2/(aa_1*om_0))


!added 20/6/2013
sigma_k_l=sigma_k1/f_1+sigma_k2/f_2
sigma_om_l=sigma_om1/f_1+sigma_om2/f_2



  if (eddyfr(1).gt.0.0)then
!for right
 vortet(1,1:3) = eddyfr(4:6)
 vortet(2,1:3) = eddyfr(7:9) 
 vortet(3,1:3) = eddyfr(10:12)

  ux = vortet(1,1);uy = vortet(1,2);uz = vortet(1,3)
  vx = vortet(2,1);vy = vortet(2,2);vz = vortet(2,3)
  wx = vortet(3,1);wy = vortet(3,2);wz = vortet(3,3)

  do ihgt=1,3
  do ihgj=1,3
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
  end do

svort=0.5*(vortet+tvort)
snorm=sqrt(2.0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+(svort(1,3)*svort(1,3))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))+(svort(2,3)*svort(2,3))+& 
	       (svort(3,1)*svort(3,1))+(svort(3,2)*svort(3,2))+(svort(3,3)*svort(3,3))))
		   
 wally=eddyfr(1)
 rho_0=rightv(1)
 k_0=eddyfr(2)/rightv(1)
 om_0=max(eddyfr(3)/rightv(1),1.0e-6)

 ! eddyfl(13:15)=rec_grads(4,1:3,k)
!eddyfl(16:18)=rec_grads(5,1:3,k)

 dervk_dervom=(eddyfr(13)*eddyfr(16))+(eddyfr(14)*eddyfr(17))+(eddyfr(15)*eddyfr(18))
		   
d_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    !i need derivative of k
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(2)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally)) 

f_1=tanh(phi_1**4)
f_2=tanh(phi_2**2)
re_t_sst=rho_0*k_0/(viscl(2)*om_0)
beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
alpha_star0_l=beta_i/3.0    
alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0_l+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)


viscl(4)=rho_0*k_0/om_0/max(1.0/alpha_star,snorm*f_2/(aa_1*om_0))

!added 20/6/2013
sigma_k_r=sigma_k1/f_1+sigma_k2/f_2
sigma_om_r=sigma_om1/f_1+sigma_om2/f_2


else
viscl(4)=-viscl(3)
sigma_k_r=sigma_k_l
sigma_om_r=sigma_om_l



end if


end select
		  
		  
    
		  viscl(3) = min(10000000*visc,viscl(3))  
		  viscl(4) = min(10000000*visc,viscl(4))		  
  

	 laml(3)=( viscl(3)*r_gas*gamma/(prtu*(gamma-1)) ) + ( viscl(1)*r_gas*gamma/(prandtl*(gamma-1)) )
	 laml(4)=( viscl(4)*r_gas*gamma/(prtu*(gamma-1)) ) + ( viscl(2)*r_gas*gamma/(prandtl*(gamma-1)) )
	 viscl(3)=max(0.0d0,viscl(3))
	 viscl(4)=max(0.0d0,viscl(4))
	 
	 if ((turbmv(1).lt.zero).or.(turbmv(2).lt.zero))then
	 viscl(3)=0.0d0
	 viscl(4)=0.0d0
	 end if
	 
	 etvm(1) = ( 0.5*(viscl(1)+viscl(2)) ) +  ( 0.5*(viscl(3)+viscl(4)) )




	  
!added on 20/6/2013---------------------------------------------------------------
!after limiting these variables, we compute the diffusion for the turbulent variables
if (turbulencemodel .eq. 2) then
!--------eddyfl/r(19)=gamma_k_l/r
!--------eddyfl/r(20)=gamma_om_l/r  

eddyfl(19)=viscl(1)+viscl(3)/sigma_k_l
eddyfr(19)=viscl(2)+viscl(4)/sigma_k_r

eddyfl(20)=viscl(1)+viscl(3)/sigma_om_l
eddyfr(20)=viscl(2)+viscl(4)/sigma_om_r
end if
end if




  end subroutine eddyvisco_ideal






subroutine eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
	implicit none
!> @brief
!> compatibility wrapper for the ideal-gas 2d eddy-viscosity implementation
#ifdef gpu
!$omp declare target
#endif
    real,dimension(1:2),intent(inout)::turbmv
    real,dimension(1:nof_variables),intent(in)::leftv,rightv
    real,dimension(1),intent(inout)::etvm
    real,dimension(1:4),intent(inout)::viscl,laml
    real,dimension(1:20),intent(inout)::eddyfl,eddyfr
	integer,intent(in)::n

	call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)

end subroutine eddyvisco2d


subroutine eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
	implicit none
!> @brief
!> this subroutine computes the tubulent eddy viscosity for turbulence models
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
    real,dimension(1:2),intent(inout)::turbmv
    real,dimension(1:nof_variables),intent(in)::leftv,rightv
    real,dimension(1),intent(inout)::etvm
    real,dimension(1:4),intent(inout)::viscl,laml
    real,dimension(1:20),intent(inout)::eddyfl,eddyfr
	integer,intent(in)::n
	real::ml,smmm,mtt,chi,tolepsma,chipow3,fv1
	integer:: ihgt, ihgj
	real,dimension(2,2)::vortet,tvort,svort
	real:: ux,uy,uz,vx,vy,vz,wx,wy,wz
		real:: wally, d_omplus, phi_2, phi_1, f_1, f_2, k_0, om_0,alpha_inf, alpha_star, alpha_star0_l, mu_turb
	real:: sigma_k_l, sigma_om_l, sigma_k_r, sigma_om_r
	real::snorm,dervk_dervom,re_t_sst,rho_0,beta_i




	tolepsma = tolsmall
	
	if (turbulence.eq.0)then
	viscl(4)=zero;viscl(3)=zero
	laml(4)=zero;laml(3)=zero
	else
	
	
	
!modified on 19/6/2013
  select case(turbulencemodel)
  
   case(1)
            if (ispal.eq.2)then
	  turbmv(1)=eddyfl(2)
	  turbmv(2)=eddyfr(2)  
	  
	  chi     = abs ((turbmv(1)) / (viscl(1)))
         !chi     = abs ( max(turbmv(1),tolepsma) / max(viscl(1),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
	viscl(3) = turbmv(1)*fv1
	 chi     = abs ((turbmv(2)) / (viscl(2)))
! 	chi     = abs ( max(turbmv(2),tolepsma) / max(viscl(2),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
		viscl(4) = turbmv(2)*fv1
        else
        turbmv(1)=eddyfl(2)
	  turbmv(2)=eddyfr(2)  
	  
	 
         chi     = abs ( max(turbmv(1),tolepsma) / max(viscl(1),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
	viscl(3) = turbmv(1)*fv1
! 	 chi     = abs ((turbmv(2)) / (viscl(2)))
 	chi     = abs ( max(turbmv(2),tolepsma) / max(viscl(2),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
		viscl(4) = turbmv(2)*fv1
        
        
        end if

  case(2)

  
  !for left cell
  
 vortet(1,1:2) = eddyfl(4:5)
 vortet(2,1:2) = eddyfl(6:7) 
 

  ux = vortet(1,1);uy = vortet(1,2)
  vx = vortet(2,1);vy = vortet(2,2)
 

  do ihgt=1,2
  do ihgj=1,2
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
  end do

svort=0.5*(vortet+tvort)
snorm=sqrt(2.0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))))
		   
 wally=eddyfl(1)
 rho_0=leftv(1)
 k_0=max(tolepsma,eddyfl(2)/leftv(1))
 om_0=max(eddyfl(3)/leftv(1),ufreestream/charlength/10.0)


 dervk_dervom=(eddyfl(8)*eddyfl(10))+(eddyfl(9)*eddyfl(11))
		   
d_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(1)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally)) 

f_1=tanh(phi_1**4)
f_2=tanh(phi_2**2)
re_t_sst=rho_0*k_0/(viscl(1)*om_0)
beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
alpha_star0_l=beta_i/3.0    
alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0_l+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)


viscl(3)=rho_0*k_0/om_0/max(1.0/alpha_star,snorm*f_2/(aa_1*om_0))


!added 20/6/2013
sigma_k_l=sigma_k1/f_1+sigma_k2/f_2
sigma_om_l=sigma_om1/f_1+sigma_om2/f_2



  if (eddyfr(1).gt.0.0)then
vortet(1,1:2) = eddyfl(4:5)
 vortet(2,1:2) = eddyfl(6:7) 
 

  ux = vortet(1,1);uy = vortet(1,2)
  vx = vortet(2,1);vy = vortet(2,2)
 

  do ihgt=1,2
  do ihgj=1,2
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
  end do

svort=0.5*(vortet+tvort)
snorm=sqrt(2.0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))))
		   
 wally=eddyfr(1)
 rho_0=rightv(1)
 k_0=eddyfr(2)/rightv(1)
 om_0=max(eddyfr(3)/rightv(1),1.0e-6)

 ! eddyfl(13:15)=rec_grads(4,1:3,k)
!eddyfl(16:18)=rec_grads(5,1:3,k)

 dervk_dervom=(eddyfl(8)*eddyfl(10))+(eddyfl(9)*eddyfl(11))
		   
d_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    !i need derivative of k
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(2)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*d_omplus*wally*wally)) 

f_1=tanh(phi_1**4)
f_2=tanh(phi_2**2)
re_t_sst=rho_0*k_0/(viscl(2)*om_0)
beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
alpha_star0_l=beta_i/3.0    
alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0_l+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)


viscl(4)=rho_0*k_0/om_0/max(1.0/alpha_star,snorm*f_2/(aa_1*om_0))

!added 20/6/2013
sigma_k_r=sigma_k1/f_1+sigma_k2/f_2
sigma_om_r=sigma_om1/f_1+sigma_om2/f_2


else
viscl(4)=-viscl(3)
sigma_k_r=sigma_k_l
sigma_om_r=sigma_om_l



end if


end select
		  
		  
    
		  viscl(3) = min(10000000*visc,viscl(3))  
		  viscl(4) = min(10000000*visc,viscl(4))		  
  

	 laml(3)=( viscl(3)*r_gas*gamma/(prtu*(gamma-1)) ) + ( viscl(1)*r_gas*gamma/(prandtl*(gamma-1)) )
	 laml(4)=( viscl(4)*r_gas*gamma/(prtu*(gamma-1)) ) + ( viscl(2)*r_gas*gamma/(prandtl*(gamma-1)) )
	 viscl(3)=max(0.0d0,viscl(3))
	 viscl(4)=max(0.0d0,viscl(4))
	 
	 if ((turbmv(1).lt.zero).or.(turbmv(2).lt.zero))then
	 viscl(3)=0.0d0
	 viscl(4)=0.0d0
	 end if
	 
	 etvm(1) = ( 0.5*(viscl(1)+viscl(2)) ) +  ( 0.5*(viscl(3)+viscl(4)) )




	  
!added on 20/6/2013---------------------------------------------------------------
!after limiting these variables, we compute the diffusion for the turbulent variables
if (turbulencemodel .eq. 2) then
!--------eddyfl/r(19)=gamma_k_l/r
!--------eddyfl/r(20)=gamma_om_l/r  

eddyfl(12)=viscl(1)+viscl(3)/sigma_k_l
eddyfr(13)=viscl(2)+viscl(4)/sigma_k_r

eddyfl(12)=viscl(1)+viscl(3)/sigma_om_l
eddyfr(13)=viscl(2)+viscl(4)/sigma_om_r
end if
end if




end subroutine eddyvisco2d_ideal



  

  

subroutine trajectories
implicit none
integer::i,j,k,traj1,traj2,traj3,traj4,kmaxe,writeid,writeconf,num_vg
real::win1,win2,win3,win4,post,post1,post2,post3,post4
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
kmaxe=xmpielrank(n)
post1=tolbig
post2=tolbig
post3=-tolbig
traj1=0
traj2=0
traj3=0

if (dimensiona.eq.3)then
num_vg=5
else
num_vg=4

end if


if (initcond.eq.157)then
pos_l(1:2)=zero
pos_g(1:2)=zero
do i=1,kmaxe
      leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      if (dimensiona.eq.3)then
      call cons2prim(n,leftv,mp_pinfl,gammal)
      else
      call cons2prim(n,leftv,mp_pinfl,gammal)
      end if

     if ((leftv(nof_variables)).gt.0.0d0)then  !total volume of gas evolution
            pos_l(2)=pos_l(2)+(leftv(nof_variables))*ielem_totvolume(i)
     end if
     if (leftv(num_vg).gt.pos_l(1))then  !total volume of gas evolution
            pos_l(1)=max(pos_l(1),leftv(num_vg))
     end if



end do
     !find position globally
call mpi_allreduce(pos_l(1),pos_g(1),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
call mpi_allreduce(pos_l(2),pos_g(2),1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)





if (n.eq.0)then

open(70,file='volumex.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),pos_g(2)
close(70)

end if

call mpi_barrier(mpi_comm_world,ierror)






end if




if (initcond.eq.430)then
pos_l(1:2)=zero
pos_g(1:2)=zero
do i=1,kmaxe
    if (u_c_val(1,7,i).gt.0.1d0)then   !volume fraction of gas to be used for lowest location tracking
        if (ielem_yyc(i).le.post3)then
            post3=ielem_yyc(i)
            traj1=i
        end if
     end if   
     if (u_c_val(1,7,i).gt.0.0d0)then  !total volume of gas evolution
            pos_l(2)=pos_l(2)+u_c_val(1,7,i)*ielem_totvolume(i)
     end if
end do

pos_l(1)=post3  !lowest position in my local  cpu


!find position globally
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(pos_l(1),pos_g(1),1,mpi_double_precision,mpi_min,mpi_comm_world,ierror)

!now select which cpu has the smallest
writeid=100000000
if (abs(pos_g(1)-pos_l(1)).le.tolsmall/1000)then
if (traj1.gt.0)then
writeid=n
end if
end if
!and if more than one, select the smallest id one
call mpi_barrier(mpi_comm_world,ierror)
call mpi_allreduce(writeid,writeconf,1,mpi_integer,mpi_min,mpi_comm_world,ierror)



if (writeid.eq.writeconf)then

leftv(1:nof_variables)=u_c_val(1,1:nof_variables,traj1)

call cons2prim(n,leftv,mp_pinfl,gammal)
open(70,file='position.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
close(70)
!the cpu that holds this cell will write its position
end if



call mpi_barrier(mpi_comm_world,ierror)
  
call mpi_allreduce(pos_l(2:2),pos_g(2:2),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

if (n.eq.0)then

open(70,file='volume.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7)')t,pos_g(2)
close(70)

end if

call mpi_barrier(mpi_comm_world,ierror)




end if


if (initcond.eq.405)then


do i=1,kmaxe
    if (dimensiona.eq.3)then
    if (u_c_val(1,8,i).gt.0.1d0)then
        if (ielem_xxc(i).le.post1)then
            post1=ielem_xxc(i)
            traj1=i
        end if
        if (((ielem_yyc(i).le.0.0515).and.(ielem_yyc(i).ge.0.0485)).and.((ielem_zzc(i).le.0.0515).and.(ielem_zzc(i).ge.0.0485)))then
            if(ielem_xxc(i).le.post2)then
                post2=ielem_xxc(i)
                traj2=i
            end if
            if (ielem_xxc(i).ge.post3)then
                post3=ielem_xxc(i)
                traj3=i
            end if
        end if
    end if
    else
    if (u_c_val(1,7,i).ge.0.4d0)then
        if (ielem_xxc(i).le.post1)then
            post1=ielem_xxc(i)
            traj1=i
        end if
        if (((ielem_yyc(i).le.0.053).and.(ielem_yyc(i).ge.0.048)))then
            if(ielem_xxc(i).le.post2)then
                post2=ielem_xxc(i)
                traj2=i
            end if
            if (ielem_xxc(i).ge.post3)then
                post3=ielem_xxc(i)
                traj3=i
            end if
        end if
    end if


    end if
end do

pos_l(1)=post1
pos_l(2)=post2
pos_l(3)=post3

call mpi_barrier(mpi_comm_world,ierror)
  
call mpi_allreduce(pos_l(1:2),pos_g(1:2),2,mpi_double_precision,mpi_min,mpi_comm_world,ierror)

!now select which cpu has the smallest
writeid=100000000
writeconf=-1
if (abs(pos_g(1)-pos_l(1)).le.tolsmall/1000)then
if (traj1.gt.0)then
writeid=n
end if
end if
!and if more than one, select the smallest id one

!if (traj1.gt.0)then
call mpi_allreduce(writeid,writeconf,1,mpi_integer,mpi_min,mpi_comm_world,ierror)
!end if



if (writeid.eq.writeconf)then

        if (traj1.gt.0)then
leftv(1:nof_variables)=u_c_val(1,1:nof_variables,traj1)


if (dimensiona.eq.3)then
call cons2prim(n,leftv,mp_pinfl,gammal)

open(70,file='pos1.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7),leftv(8)
close(70)
else
call cons2prim(n,leftv,mp_pinfl,gammal)

open(70,file='pos1.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
close(70)

end if




end if




end if



writeid=100000000
writeconf=-1
if (abs(pos_g(2)-pos_l(2)).le.tolsmall/1000)then
if (traj2.gt.0)then
writeid=n
end if
end if

!if (traj2.gt.0)then
call mpi_allreduce(writeid,writeconf,1,mpi_integer,mpi_min,mpi_comm_world,ierror)
!end if


if (writeid.eq.writeconf)then

        if (traj2.gt.0)then
leftv(1:nof_variables)=u_c_val(1,1:nof_variables,traj2)


if (dimensiona.eq.3)then

call cons2prim(n,leftv,mp_pinfl,gammal)

open(71,file='pos2.dat',form='formatted',action='write',position='append')
write(71,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(2),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7),leftv(8)
close(71)

else
call cons2prim(n,leftv,mp_pinfl,gammal)

open(70,file='pos2.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(2),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
close(70)


end if



end if
end if


call mpi_allreduce(pos_l(3),pos_g(3),1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)



writeid=100000000
writeconf=-1
if (abs(pos_g(3)-pos_l(3)).le.tolsmall/1000)then
if (traj3.gt.0)then
writeid=n
end if
end if


!if (traj3.gt.0)then
call mpi_allreduce(writeid,writeconf,1,mpi_integer,mpi_min,mpi_comm_world,ierror)
!end if



if (writeid.eq.writeconf)then
        if (traj3.gt.0)then
leftv(1:nof_variables)=u_c_val(1,1:nof_variables,traj3)


if (dimensiona.eq.3)then



call cons2prim(n,leftv,mp_pinfl,gammal)

open(71,file='pos3.dat',form='formatted',action='write',position='append')
write(71,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(3),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7),leftv(8)
close(71)


else
call cons2prim(n,leftv,mp_pinfl,gammal)

open(70,file='pos3.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(3),leftv(1),leftv(2),leftv(3),leftv(4),leftv(5),leftv(6),leftv(7)
close(70)

end if




end if






end if

call mpi_barrier(mpi_comm_world,ierror)


end if

if ((initcond.eq.408).or.(initcond.eq.422))then

pos_l(1:2)=zero
pos_g(1:2)=zero

do i=1,kmaxe
    if (u_c_val(1,8,i).gt.0.0d0)then
            pos_l(1)=pos_l(1)+ielem_totvolume(i)
            pos_l(2)=pos_l(2)+u_c_val(1,8,i)*ielem_totvolume(i)
    end if
end do


call mpi_barrier(mpi_comm_world,ierror)
  
call mpi_allreduce(pos_l(1:2),pos_g(1:2),2,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

if (n.eq.0)then

open(70,file='volume.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),pos_g(2)
close(70)

end if

call mpi_barrier(mpi_comm_world,ierror)


end if


if ((initcond.eq.411).or.(initcond.eq.444))then

pos_l(1:3)=zero
pos_g(1:3)=zero

do i=1,kmaxe
    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
    if (dimensiona.eq.2)then
    call cons2prim(n,leftv,mp_pinfl,gammal)
    pos_l(3)=max(abs(leftv(4)),pos_l(3))

    if (u_c_val(1,7,i).gt.0.0d0)then
            pos_l(1)=pos_l(1)+ielem_totvolume(i)
            pos_l(2)=pos_l(2)+u_c_val(1,7,i)*ielem_totvolume(i)
    end if
    else
    call cons2prim(n,leftv,mp_pinfl,gammal)
    pos_l(3)=max(abs(leftv(5)),pos_l(3))

    if (u_c_val(1,8,i).gt.0.0d0)then
            pos_l(1)=pos_l(1)+ielem_totvolume(i)
            pos_l(2)=pos_l(2)+u_c_val(1,7,i)*ielem_totvolume(i)
    end if


    end if
    
end do


call mpi_barrier(mpi_comm_world,ierror)
  
call mpi_allreduce(pos_l(1:2),pos_g(1:2),2,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)

call mpi_allreduce(pos_l(3),pos_g(3),1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)

if (n.eq.0)then

open(70,file='volume.dat',form='formatted',action='write',position='append')
write(70,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,pos_g(1),pos_g(2),pos_g(3)
close(70)

end if

call mpi_barrier(mpi_comm_world,ierror)


end if


end subroutine

subroutine flux2dx(flux_term_x,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real::p,u,v,w,e,r,s,gm,skin,ien,pi, ie1, mp_stiff, mp_density,gammal,gammar,sum1,sum2,sum3
real,dimension(gpu_max_species)::mp_ar,mp_ie
integer::rg_i,rg_j
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_x
real,dimension(gpu_max_species)::mp_vft

if(multispecies.eq.1)then


            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1
            do rg_i=1,nof_species-1
            sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
            mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
            end do
            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+mp_ar(rg_i)
            end do
            gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
             sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
            end do

            mp_stiff=sum3*(gammal-1.0d0)
             sum2=zero
            do rg_i=1,nof_species
            sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
            end do

 
r=mp_density
u=leftv(2)
v=leftv(3)
p=leftv(4)

!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy 
ie1=((p+mp_stiff)/((gammal-1.0d0)*r))
!total energy
e=r*(skin+ie1)
flux_term_x(1)=r*u
flux_term_x(2)=(r*(u**2))+p
flux_term_x(3)=r*u*v
flux_term_x(4)=u*(e+p)
flux_term_x(5:nof_variables)=leftv(5:nof_variables)*u
else

r=leftv(1)
u=leftv(2)
v=leftv(3)
p=leftv(4)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)
flux_term_x(1)=r*u
flux_term_x(2)=(r*(u**2))+p
flux_term_x(3)=r*u*v
flux_term_x(4)=u*(e+p)

end if

end subroutine


subroutine flux2dy(flux_term_y,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real::p,u,v,w,e,r,s,gm,skin,ien,pi, ie1, mp_stiff, mp_density,gammal,gammar,sum1,sum2,sum3
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_y
if(multispecies.eq.1)then


            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1
            do rg_i=1,nof_species-1
            sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
            mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
            end do
            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+mp_ar(rg_i)
            end do
            gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
             sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
            end do

            mp_stiff=sum3*(gammal-1.0d0)
             sum2=zero
            do rg_i=1,nof_species
            sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
            end do










! mp_ar(1)=leftv(7)/(gamma_in(1)-1.0d0)
! mp_ar(2)=(1.0d0-leftv(7))/(gamma_in(2)-1.0d0)
! gammal=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumptio
! mp_stiff=((leftv(7)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(7))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
! mp_density = leftv(5)+leftv(6)
 
r=mp_density
u=leftv(2)
v=leftv(3)
p=leftv(4)

!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy 
ie1=((p+mp_stiff)/((gammal-1.0d0)*r))
!total energy
e=r*(skin+ie1)
flux_term_y(1)=r*v
flux_term_y(2)=r*u*v
flux_term_y(3)=(r*(v**2))+p
flux_term_y(4)=v*(e+p)
flux_term_y(5:7)=leftv(5:7)*v
else


r=leftv(1)
u=leftv(2)
v=leftv(3)
p=leftv(4)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2))
!internal energy
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)

flux_term_y(1)=r*v
flux_term_y(2)=r*u*v
flux_term_y(3)=(r*v*v)+p
flux_term_y(4)=v*(e+p)

endif

end subroutine







subroutine flux3dx(flux_term_x,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real::p,u,v,w,e,r,s,gm,skin,ien,pi,ie1,mp_stiff,mp_density,gammal,gammar,sum1,sum2,sum3
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_x

if(multispecies.eq.1)then

 sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1
            do rg_i=1,nof_species-1
            sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
            mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
            end do
            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+mp_ar(rg_i)
            end do
            gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
             sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
            end do

            mp_stiff=sum3*(gammal-1.0d0)
             sum2=zero
            do rg_i=1,nof_species
            sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
            end do



! mp_ar(1)=leftv(8)/(gamma_in(1)-1.0d0)
! mp_ar(2)=(1.0d0-leftv(8))/(gamma_in(2)-1.0d0)
! gammal=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumptio
! mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)


! mp_density = leftv(6)+leftv(7)

r=mp_density

u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)

!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy
ie1=((p+mp_stiff)/((gammal-1.0d0)*r))
!total energy
e=r*(skin+ie1)
flux_term_x(1)=r*u
flux_term_x(2)=(r*(u**2))+p
flux_term_x(3)=r*u*v
flux_term_x(4)=r*u*w
flux_term_x(5)=u*(e+p)
flux_term_x(6:8)=leftv(6:8)*u
else
r=leftv(1)
u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)

flux_term_x(1)=r*u
flux_term_x(2)=(r*(u**2))+p
flux_term_x(3)=r*u*v
flux_term_x(4)=r*u*w
flux_term_x(5)=u*(e+p)

end if

end subroutine



subroutine flux3dy(flux_term_y,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real::p,u,v,w,e,r,s,gm,skin,ien,pi,ie1, mp_stiff,mp_density,gammal,gammar,sum1,sum2,sum3
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_y

if(multispecies.eq.1)then



 sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1
            do rg_i=1,nof_species-1
            sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
            mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
            end do
            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+mp_ar(rg_i)
            end do
            gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
             sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
            end do

            mp_stiff=sum3*(gammal-1.0d0)
             sum2=zero
            do rg_i=1,nof_species
            sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
            end do







! mp_ar(1)=leftv(8)/(gamma_in(1)-1.0d0)
! mp_ar(2)=(1.0d0-leftv(8))/(gamma_in(2)-1.0d0)
! gammal=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumptio
! mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
!
!
! mp_density = leftv(6)+leftv(7)

r=mp_density


u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)

!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy
ie1=((p+mp_stiff)/((gammal-1.0d0)*r))
!total energy
e=r*(skin+ie1)
flux_term_y(1)=r*v
flux_term_y(2)=r*u*v
flux_term_y(3)=(r*(v**2))+p
flux_term_y(4)=r*v*w
flux_term_y(5)=v*(e+p)
flux_term_y(6:8)=leftv(6:8)*v
else
r=leftv(1)
u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)

flux_term_y(1)=r*v
flux_term_y(2)=r*u*v
flux_term_y(3)=(r*(v*v))+p
flux_term_y(4)=r*v*w
flux_term_y(5)=v*(e+p)

end if
end subroutine

subroutine flux3dz(flux_term_z,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real::p,u,v,w,e,r,s,gm,skin,ien,pi,ie1, mp_stiff,mp_density,gammal,gammar,sum1,sum2,sum3
real,dimension(gpu_max_species)::mp_ar,mp_ie,mp_vft
integer::rg_i,rg_j
real,dimension(1:nof_variables),intent(in)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_z

if(multispecies.eq.1)then



 sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+leftv(dimensiona+2+rg_i)
            end do
            mp_density=sum1
            do rg_i=1,nof_species-1
            sum2=sum2+leftv(dimensiona+2+nof_species+rg_i)
            mp_ar(rg_i)=leftv(dimensiona+2+nof_species+rg_i)/(gamma_in(rg_i)-1.0d0)
            mp_vft(rg_i)=leftv(dimensiona+2+nof_species+rg_i)
            end do
            mp_ar(nof_species)=(1.0d0-sum2)/(gamma_in(nof_species)-1.0d0)
            mp_vft(nof_species)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+mp_ar(rg_i)
            end do
            gammal=(1.0d0/(sum3))+1.0d0    !mixture gamma isobaric assumption
             sum3=zero
            do rg_i=1,nof_species
            sum3=sum3+(mp_vft(rg_i)*(gamma_in(rg_i)/(gamma_in(rg_i)-1.0d0))*mp_pinf(rg_i))
            end do

            mp_stiff=sum3*(gammal-1.0d0)
             sum2=zero
            do rg_i=1,nof_species
            sum2=sum2+(mp_vft(rg_i)*mp_pinf(rg_i))
            end do










! mp_ar(1)=leftv(8)/(gamma_in(1)-1.0d0)
! mp_ar(2)=(1.0d0-leftv(8))/(gamma_in(2)-1.0d0)
! gammal=(1.0d0/(mp_ar(1)+mp_ar(2)))+1.0d0    !mixture gamma isobaric assumptio
! mp_stiff=((leftv(8)*(gamma_in(1)/(gamma_in(1)-1.0d0))*mp_pinf(1))+((1.0d0-leftv(8))*(gamma_in(2)/(gamma_in(2)-1.0d0))*mp_pinf(2)))*(gammal-1.0d0)
!
! mp_density = leftv(6)+leftv(7)
!
 r=mp_density
u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)

!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy
ie1=((p+mp_stiff)/((gammal-1.0d0)*r))
!total energy
e=r*(skin+ie1)
flux_term_z(1)=r*w
flux_term_z(2)=r*u*w
flux_term_z(3)=r*v*w
flux_term_z(4)=(r*(w*w))+p
flux_term_z(5)=w*(e+p)
flux_term_z(6:8)=leftv(6:8)*w
else
r=leftv(1)
u=leftv(2)
v=leftv(3)
w=leftv(4)
p=leftv(5)
gm=gamma
!kinetic energy first!
skin=(oo2)*((u**2)+(v**2)+(w**2))
!internal energy 
ien=((p)/((gm-1.0d0)*r))
!total energy
e=r*(skin+ien)

flux_term_z(1)=r*w
flux_term_z(2)=r*u*w
flux_term_z(3)=r*v*w
flux_term_z(4)=(r*(w**2))+p
flux_term_z(5)=w*(e+p)
end if
end subroutine






subroutine flux_visc2d(flux_term_x,flux_term_y,leftv,leftv_der)
implicit none
#ifdef gpu
!$omp declare target
#endif
real:: u, v, ux, uy, vx, vy, tx, ty, tauxx, tauxy, tauyy
real,dimension(1:nof_variables),intent(inout)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_x,flux_term_y
real,dimension(1:4)::viscl,laml
real,dimension(1:nof_variables,1:dimensiona),intent(in)::leftv_der

    call get_visc_conduct(n, leftv, leftv,viscl,laml)

       u = leftv(2)
       v = leftv(3)
      ux = leftv_der(2,1)
      uy = leftv_der(2,2)
      vx = leftv_der(3,1)
      vy = leftv_der(3,2)
      tx = leftv_der(4,1)
      ty = leftv_der(4,2)

    tauxx = 2.0d0 / 3.0d0 * viscl(1) * (2 * ux - vy)
    tauyy = 2.0d0 / 3.0d0 * viscl(1) * (2 * vy - ux)
    tauxy = viscl(1) * (uy + vx)

    flux_term_x(2) = flux_term_x(2) - tauxx
    flux_term_x(3) = flux_term_x(3) - tauxy
    flux_term_x(4) = flux_term_x(4) - (u * tauxx + v * tauxy + laml(1) * tx)
    flux_term_y(2) = flux_term_y(2) - tauxy
    flux_term_y(3) = flux_term_y(3) - tauyy
    flux_term_y(4) = flux_term_y(4) - (u * tauxy + v * tauyy + laml(1) * ty)


end subroutine


subroutine flux_visc3d(flux_term_x,flux_term_y,flux_term_z,leftv,leftv_der)
implicit none
#ifdef gpu
!$omp declare target
#endif
real:: u, v, ux, uy, vx, vy, tx, ty, tauxx, tauxy, tauyy, w, uz, vz, wx, wy, wz, tz, tauxz, tauyz, tauzz
real,dimension(1:nof_variables),intent(inout)::leftv
real,dimension(1:nof_variables),intent(inout)::flux_term_x,flux_term_y,flux_term_z
real,dimension(1:nof_variables,1:dimensiona),intent(in)::leftv_der
real,dimension(1:4)::viscl,laml


    call get_visc_conduct(n, leftv, leftv,viscl,laml)

! variables extrapolated at boundary
       u = leftv(2)
       v = leftv(3)
       w = leftv(4)

      ux = leftv_der(2,1)
      uy = leftv_der(2,2)
      uz = leftv_der(2,3)

      vx = leftv_der(3,1)
      vy = leftv_der(3,2)
      vz = leftv_der(3,3)

      wx = leftv_der(4,1)
      wy = leftv_der(4,2)
      wz = leftv_der(4,3)

      tx = leftv_der(5,1)
      ty = leftv_der(5,2)
      tz = leftv_der(5,3)

    tauxx = 2.0d0 / 3.0d0 * viscl(1) * (2 * ux - vy - wz)
    tauyy = 2.0d0 / 3.0d0 * viscl(1) * (2 * vy - ux - wz)
    tauxy = viscl(1) * (uy + vx)

    tauxz = viscl(1) * (uz + wx)
    tauyz = viscl(1) * (wy + vz)
    tauzz = 2.0d0 / 3.0d0 * viscl(1) * (2 * wz - ux - vy)

    flux_term_x(2) = flux_term_x(2) - tauxx
    flux_term_x(3) = flux_term_x(3) - tauxy
    flux_term_x(4) = flux_term_x(4) - tauxz
    flux_term_x(5) = flux_term_x(5) - (u * tauxx + v * tauxy + w * tauxz + laml(1) * tx)

    flux_term_y(2) = flux_term_y(2) - tauxy
    flux_term_y(3) = flux_term_y(3) - tauyy
    flux_term_y(4) = flux_term_y(4) - tauyz
    flux_term_y(5) = flux_term_y(5) - (u * tauxy + v * tauyy + w * tauyz + laml(1) * ty)

    flux_term_z(2) = flux_term_z(2) - tauxz
    flux_term_z(3) = flux_term_z(3) - tauyz
    flux_term_z(4) = flux_term_z(4) - tauzz
    flux_term_z(5) = flux_term_z(5) - (u * tauxz + v * tauyz + w * tauzz + laml(1) * tz)

end subroutine








subroutine dcons2dprim(leftv_der,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables,1:dimensiona),intent(inout)::leftv_der
real,dimension(1:nof_variables),intent(in)::leftv
integer:: i_dim, i_var
real:: vel_grad_dot


    do i_dim = 1, dimensiona
        do i_var = 2, nof_variables-1
            leftv_der(i_var,i_dim) = (leftv_der(i_var,i_dim) - leftv_der(1,i_dim) * leftv(i_var)) / leftv(1) ! ux = (rhoux - rhox * u) / rho
        end do

        vel_grad_dot = zero
        do i_var = 2, nof_variables-1
            vel_grad_dot = vel_grad_dot + leftv(i_var) * leftv_der(i_var,i_dim)
        end do

        leftv_der(nof_variables,i_dim) = (gamma - 1.0d0) * ((leftv(1) * leftv_der(nof_variables,i_dim) - leftv(nof_variables) * leftv_der(1,i_dim)) / leftv(1) ** 2 - vel_grad_dot)
        ! tx = (gamma-1) * ((rho * ex - e * rhox) / rho ** 2 - (u * ux + v * vx + w * wx))
    end do

end subroutine dcons2dprim





end module flow_operations
