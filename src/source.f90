module source
use library
use transform
use local
use riemann
use flow_operations
use declaration
implicit none

contains
subroutine sources_computation(n)
	implicit none
!> @brief
!> sources computation 
	integer,intent(in)::n
	integer::i,kmaxe,iconsidered
	real,dimension(turbulenceequations)::source_t

	
	
	kmaxe=xmpielrank(n)
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, iconsidered, source_t)
#else
!$omp do
#endif
	do i=1,kmaxe
		iconsidered=i

		if (turbulence.eq.1)then
		call sources(n,iconsidered,source_t)
		rhst_val(1:turbulenceequations,i)=rhst_val(1:turbulenceequations,i)-(source_t(1:turbulenceequations)*ielem_totvolume(i))

		end if


	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine sources_computation


subroutine sources_computation_rot(n)
	implicit none
	integer,intent(in)::n
	integer::i,kmaxe,srf
	real::oodensity
	real,dimension(5)::source_t2
	real,dimension(1:dimensiona)::pox,poy,poz

	
	
	kmaxe=xmpielrank(n)
	if(srfg.eq.1)then
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, source_t2,srf,oodensity,pox,poy,poz)
#else
!$omp do
#endif
	do i=1,kmaxe
		
		 oodensity=1.0d0/u_c_val(1,1,i)
        
                source_t2(1)=zero
                source_t2(2)=u_c_val(1,2,i)*oodensity-uvel
                source_t2(3)=u_c_val(1,3,i)*oodensity-vvel
                source_t2(4)=u_c_val(1,4,i)*oodensity-wvel
                source_t2(5)=zero
        
                pox(1:3)= source_t2(2:4)  
                poy(1:3)=srf_velocity(1:3)
                source_t2(2:4)=u_c_val(1,1,i)*vect_function(pox,poy)
		rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+(source_t2(1:nof_variables)*ielem_totvolume(i))
		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	if(mrf.eq.1)then
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, source_t2,srf,oodensity,pox,poy,poz)
#else
!$omp do
#endif
	do i=1,kmaxe
		srf=rec_mrf(i)
		if (rec_mrf(i).eq.1)then
            oodensity=1.0d0/u_c_val(1,1,i)
                source_t2(1)=zero
                source_t2(2)=u_c_val(1,2,i)*oodensity!-uvel
                source_t2(3)=u_c_val(1,3,i)*oodensity!-vvel
                source_t2(4)=u_c_val(1,4,i)*oodensity!-wvel
                source_t2(5)=zero
        
                pox(1:3)= source_t2(2:4)  
                poy(1:3)=rec_mrf_velocity(1:3,i)
                source_t2(2:4)=u_c_val(1,1,i)*vect_function(pox,poy)
		rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)+(source_t2(1:nof_variables)*ielem_totvolume(i))
		end if
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
end subroutine sources_computation_rot





subroutine sources_realgas(n,iconsidered,source_r,src_jac_diag)
  implicit none
  !> sources computation in 2d/3d for 5-species, 2-t real gas
#ifdef gpu
!$omp declare target
#endif

  integer,intent(in) :: n, iconsidered
  real,dimension(1:nof_variables),intent(inout)::source_r, src_jac_diag

  integer :: i, j, k, l, rg_i, rg_j
  real,dimension(nof_species) :: rg_r, rg_rm, rg_n      ! rg_r = ρ_s
  real :: rg_dq_rad, rg_qtv, rg_qw, velx, vely, velz
  real,dimension(nof_species) :: rg_ev, rg_esv
  real,dimension(3) :: rg_tvs
  real,dimension(nof_species) :: rg_dw_s
  real :: rg_ke1, rg_ke2, rg_ke3, rg_ke4, rg_ke5
  real :: rg_t, rg_ta, rg_tv
  real :: rg_z, rg_speed
  real :: rg_r1, rg_r2, rg_r3, rg_r4, rg_r5
  real :: rg_kf1, rg_kf2, rg_kf3, rg_kf4, rg_kf5
  real :: rg_kb1, rg_kb2, rg_kb3, rg_kb4, rg_kb5

  real :: rg_kf11,rg_kf12,rg_kf13,rg_kf14,rg_kf15
  real :: rg_kf21,rg_kf22,rg_kf23,rg_kf24,rg_kf25
  real :: rg_kf31,rg_kf32,rg_kf33,rg_kf34,rg_kf35
  real :: rg_kb11,rg_kb12,rg_kb13,rg_kb14,rg_kb15
  real :: rg_kb21,rg_kb22,rg_kb23,rg_kb24,rg_kb25
  real :: rg_kb31,rg_kb32,rg_kb33,rg_kb34,rg_kb35

  real :: kf11m,kf12m,kf13m,kf14m,kf15m
  real :: kf21m,kf22m,kf23m,kf24m,kf25m
  real :: kf31m,kf32m,kf33m,kf34m,kf35m
  real :: kb11m,kb12m,kb13m,kb14m,kb15m
  real :: kb21m,kb22m,kb23m,kb24m,kb25m
  real :: kb31m,kb32m,kb33m,kb34m,kb35m
  real :: kf4m,kf5m,kb4m,kb5m
  real :: kf1_eff,kb1_eff,kf2_eff,kb2_eff,kf3_eff,kb3_eff
  real :: rg_c(5)            ! molar concentrations [mol/m^3]
  real :: wdot(5)            ! reaction rates [mol/(m^3 s)]
  integer :: nu(5,5)         ! stoichiometric matrix ν(s,r)
  integer :: s, r

  real,dimension(nof_species,nof_species) :: rg_mu
  real :: temp1, temp2,taumw,taupk
  real,dimension(nof_species) :: rg_tsmw, rg_tsp, rg_tssum, rg_ablot, rg_ev_eq
  real,dimension(1:nof_variables) :: leftv, tempvect
  real :: mp_pinfl, gammal,q_chem
  real :: p_atm, t_13, rg_pressure, n_tot, exp_argtv, exp_argt, rs
  integer :: rg_molx
  real,dimension(1:nof_species) :: rg_mass_amu

  real, parameter :: kb = 1.380649d-23     ! boltzmann [j/k]
  real, parameter :: na = 6.02214076d23    ! avogadro [1/mol]

  real,dimension(1:nof_species) :: tau_mw    !! millikan-white [s]
  real,dimension(1:nof_species) :: tau_park  ! park [s]
  real,dimension(1:nof_species) :: tau_tot   ! combined [s]
  real :: num, den,rggtf
  real :: mu_sr, a_sr, b_sr
  real :: log_p_tau, p_tau, tau_sr
  real :: sigma0, sigma_v, v_th, m_s,maxloss
  real :: rg_alpha, rg_alpha_i, rg_rho_min,sum_sourcex
  integer :: idxe, idxev, idxy1,idxn2,idx
 real    :: rho_min, lambda_s, lambda_e
 real    :: tau_v_mix, rhoe_local, yi, sumy,q_ttv,tv_mix,t_eff
  real :: dt_loc
   real :: rhoe_old, rhoe_new
   real :: rhoev_old, rhoev_new, rhoev_eq
   real :: alpha_v
   real :: prod, dest, lambda_k, rhoy_old, rhoy_new
   real :: tau_v_eff,kvt,rhoy_n2,rho,sumrhoy
   real :: rho_vib, b_v, lambda_v


  ! ---------------------------------------------------------------------------
  ! start
  ! ---------------------------------------------------------------------------
  i = iconsidered

  source_r(1:nof_variables) = 0.0d0
  tau_mw(:)   = 1.0d30
  tau_park(:) = 1.0d30
  tau_tot(:)  = 1.0d30

  !-----------------------------------------------------------
  ! get conservative variables and primitive variables
  !-----------------------------------------------------------
  leftv(1:nof_variables) = u_c_val(1,1:nof_variables,i)

  ! species densities ρ_s = ρ y_s
  do rg_i = 1, nof_species
     rg_r(rg_i) = leftv(dimensiona+3+rg_i)
  end do

  !-----------------------------------------------------------
  ! 1) molar concentrations c_s = rho_s / m_s
  !-----------------------------------------------------------
  do s = 1, 5
     rg_c(s) = rg_r(s) / rg_molm(s)
  end do

  call cons2prim(n,leftv,mp_pinfl,gammal)

  rg_pressure = leftv(dimensiona+2)
  p_atm       = rg_pressure / rgs_pa_per_atm

  ! molar masses in kg/mol -> "amu-like" for mw reduced-mass formula
  rg_mass_amu = rg_molm * 1000.0d0

  tempvect(1:nof_variables) = u_c_val(1,1:nof_variables,i)
  call cons2div(n,tempvect,mp_pinfl,gammal)

  ! leftv:   (ρ, u, v, ρe, ρev, ρy1..ρy5)
  ! tempvect:(ρ, u, v, ttr, tv, y1..y5)

  ! temperatures
  rg_t  = tempvect(dimensiona+2)   ! ttr
  rg_tv = tempvect(dimensiona+3)   ! tv
  rg_z  = 1.0d0 / rg_t
  t_13  = rg_t**(-1.0d0/3.0d0)

  ! adjusted two-temperature model (candler/park style)
  rg_ta = rg_t**0.6d0 * rg_tv**0.4d0

  ! flow speed (not used in source terms yet)
  velx = leftv(2)
  vely = leftv(3)
  velz = 0.0d0
  if (dimensiona.eq.3) then
     velz = leftv(4)
  end if
  rg_speed = sqrt(velx**2 + vely**2 + velz**2)

  !-----------------------------------------------------------
  ! chemical source: reaction rates (candler 5-reaction model)
  ! using ρ-based rate coefficients cf [m^3/(kg s)]
  !-----------------------------------------------------------
  rg_dw_s(:) = 0.0d0



if (rg_nof_reactions == 5) then

  !
  ! candler 5-reaction model (your original code)
  !

  ! equilibrium constants k_e(t) (candler curve fits, t = ttr)
  rg_ke1 = exp(3.898d0 -12.611d0*rg_z +0.683d0*rg_z**2 -0.118d0*rg_z**3 +0.006d0*rg_z**4)
  rg_ke2 = exp(1.335d0 - 4.127d0*rg_z -0.616d0*rg_z**2 +0.093d0*rg_z**3 -0.005d0*rg_z**4)
  rg_ke3 = exp(1.549d0 - 7.784d0*rg_z +0.228d0*rg_z**2 -0.043d0*rg_z**3 +0.002d0*rg_z**4)
  rg_ke4 = exp(2.349d0 - 4.828d0*rg_z +0.455d0*rg_z**2 -0.075d0*rg_z**3 +0.004d0*rg_z**4)
  rg_ke5 = exp(0.215d0 - 3.652d0*rg_z +0.843d0*rg_z**2 -0.136d0*rg_z**3 +0.007d0*rg_z**4)

  ! --- reaction 1: n2 + m <-> 2n + m (in m^3/mol/s) ---
  kf11m = 3.78d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n2
  kf12m = 3.78d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o2
  kf13m = 3.78d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = no
  kf14m = 1.11d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n
  kf15m = 1.11d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o

  kb11m = kf11m / rg_ke1
  kb12m = kf12m / rg_ke1
  kb13m = kf13m / rg_ke1
  kb14m = kf14m / rg_ke1
  kb15m = kf15m / rg_ke1

  kf1_eff = kf11m*rg_c(1) + kf12m*rg_c(2) + kf13m*rg_c(3) &
          + kf14m*rg_c(4) + kf15m*rg_c(5)
  kb1_eff = kb11m*rg_c(1) + kb12m*rg_c(2) + kb13m*rg_c(3) &
          + kb14m*rg_c(4) + kb15m*rg_c(5)

  wdot(1) = kf1_eff*rg_c(1) - kb1_eff*rg_c(4)**2   ! n2 + m <-> 2n + m

  ! --- reaction 2: o2 + m <-> 2o + m ---
  kf21m = 2.75d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = n2
  kf22m = 2.75d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = o2
  kf23m = 2.75d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = no
  kf24m = 8.25d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = n
  kf25m = 8.25d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = o

  kb21m = kf21m / rg_ke2
  kb22m = kf22m / rg_ke2
  kb23m = kf23m / rg_ke2
  kb24m = kf24m / rg_ke2
  kb25m = kf25m / rg_ke2

  kf2_eff = kf21m*rg_c(1) + kf22m*rg_c(2) + kf23m*rg_c(3) &
          + kf24m*rg_c(4) + kf25m*rg_c(5)
  kb2_eff = kb21m*rg_c(1) + kb22m*rg_c(2) + kb23m*rg_c(3) &
          + kb24m*rg_c(4) + kb25m*rg_c(5)

  wdot(2) = kf2_eff*rg_c(2) - kb2_eff*rg_c(5)**2   ! o2 + m <-> 2o + m

  ! --- reaction 3: no + m <-> n + o + m ---
  kf31m = 2.30d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = n2
  kf32m = 2.30d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = o2
  kf33m = 2.30d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = no
  kf34m = 4.60d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = n
  kf35m = 4.60d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = o

  kb31m = kf31m / rg_ke3
  kb32m = kf32m / rg_ke3
  kb33m = kf33m / rg_ke3
  kb34m = kf34m / rg_ke3
  kb35m = kf35m / rg_ke3

  kf3_eff = kf31m*rg_c(1) + kf32m*rg_c(2) + kf33m*rg_c(3) &
          + kf34m*rg_c(4) + kf35m*rg_c(5)
  kb3_eff = kb31m*rg_c(1) + kb32m*rg_c(2) + kb33m*rg_c(3) &
          + kb34m*rg_c(4) + kb35m*rg_c(5)

  wdot(3) = kf3_eff*rg_c(3) - kb3_eff*rg_c(4)*rg_c(5) ! no + m <-> n + o + m

  ! --- reaction 4: n2 + o <-> no + n (bimolecular, 2-body) ---
  kf4m = 3.18d10 * 1.0d-6 * rg_ta**(0.1d0) * exp(-3.77d4/rg_ta)
  kb4m = kf4m / rg_ke4

  wdot(4) = kf4m*rg_c(1)*rg_c(5) - kb4m*rg_c(3)*rg_c(4)

  ! --- reaction 5: no + o <-> o2 + n (bimolecular) ---
  kf5m = 2.16d5 * 1.0d-6 * rg_ta**(1.29d0) * exp(-1.922d4/rg_ta)
  kb5m = kf5m / rg_ke5

  wdot(5) = kf5m*rg_c(3)*rg_c(5) - kb5m*rg_c(2)*rg_c(4)

  ! build stoichiometry & mass production (unchanged)
  do s = 1,5
     do r = 1,5
        nu(s,r) = 0
     end do
  end do

  ! r1: n2 -> 2n
  nu(1,1) = -1
  nu(4,1) = +2

  ! r2: o2 -> 2o
  nu(2,2) = -1
  nu(5,2) = +2

  ! r3: no -> n + o
  nu(3,3) = -1
  nu(4,3) = +1
  nu(5,3) = +1

  ! r4: n2 + o -> no + n
  nu(1,4) = -1
  nu(5,4) = -1
  nu(3,4) = +1
  nu(4,4) = +1

  ! r5: no + o -> o2 + n
  nu(3,5) = -1
  nu(5,5) = -1
  nu(2,5) = +1
  nu(4,5) = +1

  rg_dw_s(:) = 0.0d0
  do r = 1,5
     do s = 1,5
        rg_dw_s(s) = rg_dw_s(s) + nu(s,r) * rg_molm(s) * wdot(r)
     end do
  end do


end if




if (rg_nof_reactions == 6) then

  !
  ! park 2-t “6-reaction” model (option a, 5 neutral reactions)
  ! same topology as candler (r1..r5), park coefficients + 3rd bodies
  !

  ! equilibrium constants (keep your candler curve-fits)
  rg_ke1 = exp(3.898d0 -12.611d0*rg_z +0.683d0*rg_z**2 -0.118d0*rg_z**3 +0.006d0*rg_z**4)
  rg_ke2 = exp(1.335d0 - 4.127d0*rg_z -0.616d0*rg_z**2 +0.093d0*rg_z**3 -0.005d0*rg_z**4)
  rg_ke3 = exp(1.549d0 - 7.784d0*rg_z +0.228d0*rg_z**2 -0.043d0*rg_z**3 +0.002d0*rg_z**4)
  rg_ke4 = exp(2.349d0 - 4.828d0*rg_z +0.455d0*rg_z**2 -0.075d0*rg_z**3 +0.004d0*rg_z**4)
  rg_ke5 = exp(0.215d0 - 3.652d0*rg_z +0.843d0*rg_z**2 -0.136d0*rg_z**3 +0.007d0*rg_z**4)

  ! ---- r1: n2 + m <-> 2n + m (park 2t, table 1) ----
  !   kf = c_m t_a^{-1.6} exp(-113200/t_a), cm^3/(mol s)
  !   m = n2,o2,no : c = 7.0e21
  !   m = n,o      : c = 3.0e22
  kf11m = 7.0d21 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n2
  kf12m = 7.0d21 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o2
  kf13m = 7.0d21 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = no
  kf14m = 3.0d22 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n
  kf15m = 3.0d22 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o

  kb11m = kf11m / rg_ke1
  kb12m = kf12m / rg_ke1
  kb13m = kf13m / rg_ke1
  kb14m = kf14m / rg_ke1
  kb15m = kf15m / rg_ke1

  kf1_eff = kf11m*rg_c(1) + kf12m*rg_c(2) + kf13m*rg_c(3) &
          + kf14m*rg_c(4) + kf15m*rg_c(5)
  kb1_eff = kb11m*rg_c(1) + kb12m*rg_c(2) + kb13m*rg_c(3) &
          + kb14m*rg_c(4) + kb15m*rg_c(5)

  wdot(1) = kf1_eff*rg_c(1) - kb1_eff*rg_c(4)**2

  ! ---- r2: o2 + m <-> 2o + m (park 2t, table 1) ----
  !   kf = c_m t_a^{-1.5} exp(-59500/t_a), cm^3/(mol s)
  !   m = n2,o2,no : c = 2.0e21
  !   m = n,o      : c = 1.0e22
  kf21m = 2.0d21 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = n2
  kf22m = 2.0d21 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = o2
  kf23m = 2.0d21 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = no
  kf24m = 1.0d22 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = n
  kf25m = 1.0d22 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = o

  kb21m = kf21m / rg_ke2
  kb22m = kf22m / rg_ke2
  kb23m = kf23m / rg_ke2
  kb24m = kf24m / rg_ke2
  kb25m = kf25m / rg_ke2

  kf2_eff = kf21m*rg_c(1) + kf22m*rg_c(2) + kf23m*rg_c(3) &
          + kf24m*rg_c(4) + kf25m*rg_c(5)
  kb2_eff = kb21m*rg_c(1) + kb22m*rg_c(2) + kb23m*rg_c(3) &
          + kb24m*rg_c(4) + kb25m*rg_c(5)

  wdot(2) = kf2_eff*rg_c(2) - kb2_eff*rg_c(5)**2

  ! ---- r3: no + m <-> n + o + m (park/koshi) ----
  !   kf ≈ c_m t_a^{0} exp(-75500/t_a), cm^3/(mol s)
  !   m = no,o,n : c = 1.1e17
  !   m = n2,o2 : c = (1/22)*1.1e17  [koshi]
  kf31m = (1.1d17/22.0d0) * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = n2
  kf32m = (1.1d17/22.0d0) * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = o2
  kf33m =  1.1d17        * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = no
  kf34m =  1.1d17        * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = n
  kf35m =  1.1d17        * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = o

  kb31m = kf31m / rg_ke3
  kb32m = kf32m / rg_ke3
  kb33m = kf33m / rg_ke3
  kb34m = kf34m / rg_ke3
  kb35m = kf35m / rg_ke3

  kf3_eff = kf31m*rg_c(1) + kf32m*rg_c(2) + kf33m*rg_c(3) &
          + kf34m*rg_c(4) + kf35m*rg_c(5)
  kb3_eff = kb31m*rg_c(1) + kb32m*rg_c(2) + kb33m*rg_c(3) &
          + kb34m*rg_c(4) + kb35m*rg_c(5)

  wdot(3) = kf3_eff*rg_c(3) - kb3_eff*rg_c(4)*rg_c(5)

  ! ---- r4: n2 + o <-> no + n (park exchange, table 1) ----
  !   kf = 6.4e7 t^{-1} exp(-38370/t), cm^3/(mol s)
  kf4m = 6.4d7 * 1.0d-6 * rg_t**(-1.0d0) * exp(-3.837d4/rg_t)
  kb4m = kf4m / rg_ke4

  wdot(4) = kf4m*rg_c(1)*rg_c(5) - kb4m*rg_c(3)*rg_c(4)

  ! ---- r5: no + o <-> o2 + n (park exchange, table 1) ----
  !   kr_exo = 8.4e12 exp(-19450/t), cm^3/(mol s)
  !   take kf endothermic = kr_exo / k_eq
  kf5m = 8.4d12 * 1.0d-6 * exp(-1.945d4/rg_t)
  kb5m = kf5m / rg_ke5

  wdot(5) = kf5m*rg_c(3)*rg_c(5) - kb5m*rg_c(2)*rg_c(4)

  ! stoichiometry is the same as candler (same 5 reactions)
  do s = 1,5
     do r = 1,5
        nu(s,r) = 0
     end do
  end do

  ! r1: n2 -> 2n
  nu(1,1) = -1
  nu(4,1) = +2

  ! r2: o2 -> 2o
  nu(2,2) = -1
  nu(5,2) = +2

  ! r3: no -> n + o
  nu(3,3) = -1
  nu(4,3) = +1
  nu(5,3) = +1

  ! r4: n2 + o -> no + n
  nu(1,4) = -1
  nu(5,4) = -1
  nu(3,4) = +1
  nu(4,4) = +1

  ! r5: no + o -> o2 + n
  nu(3,5) = -1
  nu(5,5) = -1
  nu(2,5) = +1
  nu(4,5) = +1

  rg_dw_s(:) = 0.0d0
  do r = 1,5
     do s = 1,5
        rg_dw_s(s) = rg_dw_s(s) + nu(s,r) * rg_molm(s) * wdot(r)
     end do
  end do


end if










  rg_rho_min = 1.0d-20    ! or whatever you use as "zero" density
  rg_alpha   = 1.0d0

	! first: compute most restrictive scaling factor alpha
	do rg_i = 1, nof_species

		! if density is essentially zero, don't allow further destruction
		if (rg_r(rg_i) <= rg_rho_min .and. rg_dw_s(rg_i) < 0.0d0) then
			rg_dw_s(rg_i) = 0.0d0
			cycle
		end if

		! only destruction can cause negativity
! 		if (rg_dw_s(rg_i) < 0.0d0) then
! 			! need: rho_i + dt * alpha * dw_i >= 0
! 			!  => alpha <= rho_i / (-dt * dw_i)
! 			rg_alpha_i = rg_r(rg_i) / (-ielem_dtl(iconsidered) * rg_dw_s(rg_i))
!
! 			if (rg_alpha_i < 1.0d-2) rg_alpha = 1.0d-2
! 		end if
	end do

	! clamp alpha to [0,1]
	if (rg_alpha > 1.0d0) rg_alpha = 1.0d0
	if (rg_alpha < 0.0d0) rg_alpha = 0.0d0

	! second: apply scaling if needed
	if (rg_alpha < 1.0d0) then
		do rg_i = 1, nof_species
			rg_dw_s(rg_i) = rg_alpha * rg_dw_s(rg_i)
		end do
	end if

	! optional tiny cutoff (for cleanliness, won’t break stoichiometry)
	do rg_i = 1, nof_species
		if (abs(rg_dw_s(rg_i)) < 1.0d-50) rg_dw_s(rg_i) = 0.0d0
	end do






  !-----------------------------------------------------------
  ! 1) millikan-white vt relaxation times (species 1..3)
  !-----------------------------------------------------------
  do rg_i = 1, 3
     num = 0.0d0
     den = 0.0d0

     do rg_j = 1, nof_species
        if (rg_r(rg_j) <= 0.0d0) cycle

        ! reduced mass μ_sr in amu
        mu_sr = rg_mass_amu(rg_i)*rg_mass_amu(rg_j) / (rg_mass_amu(rg_i)+rg_mass_amu(rg_j))

        ! mw coefficients
        a_sr = 0.00116d0 * sqrt(mu_sr) * rg_thetag(rg_i)**(4.0d0/3.0d0)
        b_sr = 0.015d0   * mu_sr**0.25d0

        ! p * tau_sr in atm*s
        log_p_tau = a_sr*(t_13 - b_sr) - 18.42d0
        p_tau       = exp(log_p_tau)

        ! tau_sr [s]: p_tau has atm*s, divide by p in atm
        tau_sr = p_tau / p_atm

        ! mixture rule with mole fractions:
        ! 1/tau_s = σ_r (x_r / tau_sr)
        ! x_r ~ n_r / σ n_r, and n_r ~ ρ_r / m_r
        num = num + (rg_r(rg_j)/rg_molm(rg_j)) / tau_sr
        den = den + (rg_r(rg_j)/rg_molm(rg_j))
     end do

     if (num > 0.0d0 .and. den > 0.0d0) then
        tau_mw(rg_i) = den / num
     else
        tau_mw(rg_i) = 1.0d30
     end if



  end do




  sigma0 = 3.0d-21    ! [m^2]




! number density [1/m^3] based on translational/rotational t
n_tot = rg_pressure / (kb * rg_t)   ! rg_t = t_tr

! single vibrational temperature for the mixture
         ! <-- your global vibrational temperature

! avoid non-physical tv
tv_mix = max(rg_tv, 300.0d0)

! park-style weighting between t and tv (0.7 is a common choice)
q_ttv = 0.7d0

! effective two-temperature average (same for all species)
! geometric mean: t_eff = t^q * tv^(1-q)
t_eff = rg_t**q_ttv * tv_mix**(1.0d0 - q_ttv)
t_eff = max(t_eff, 300.0d0)

do rg_i = 1, 3
   ! particle mass [kg/particle]
   m_s = rg_molm(rg_i) / na

   ! effective cross-section using t_eff instead of rg_t
   sigma_v = sigma0 * (50000.0d0 / t_eff)**2

   ! thermal speed using t_eff
   v_th = sqrt( 8.0d0 * kb * t_eff / (pi * m_s) )

   ! tau_park [s]
   tau_park(rg_i) = 1.0d0 / ( n_tot * sigma_v * v_th )

end do
























  !-----------------------------------------------------------
  ! 3) combine millikan-white + park
  !    1/tau_tot = 1/tau_mw + 1/tau_park
  !-----------------------------------------------------------


 do rg_i = 1, 3

  ! local copies
  taumw  = tau_mw(rg_i)
  taupk  = tau_park(rg_i)

  ! replace non-positive with 'disabled'
  if (taumw <= 0.0d0) taumw = 1.0d30
  if (taupk <= 0.0d0) taupk = 1.0d30

  ! now combine them safely
  if (taumw >= 1.0d29 .and. taupk >= 1.0d29) then
     ! both basically 'off'
     tau_tot(rg_i) = 1.0d30

  elseif (taumw >= 1.0d29) then
     ! only park valid
     tau_tot(rg_i) = taupk

  elseif (taupk >= 1.0d29) then
     ! only mw valid
     tau_tot(rg_i) = taumw

  else
     ! both valid: use harmonic sum
      rggtf = taumw / taupk
       rggtf = max(1.0d0, rggtf)
       rggtf = min(10d0, rggtf)
       tau_tot(rg_i) = taumw / rggtf
!        tau_tot(rg_i) = 1.0d0 / ( 1.0d0/taumw + 1.0d0/taupk )

  end if


  rg_tssum(rg_i) = tau_tot(rg_i)

end do

  !-----------------------------------------------------------
  ! 4) vibrational energies e_v(tv) and e_v_eq(t)
  !-----------------------------------------------------------
  do rg_i = 1, 3
     rs         = rgs_ru / rg_molm(rg_i)
     exp_argtv  = rg_thetag(rg_i) / rg_tv
     rg_ev(rg_i) = rs * rg_thetag(rg_i) / (exp(exp_argtv) - 1.0d0)

     exp_argt      = rg_thetag(rg_i) / rg_t
     rg_ev_eq(rg_i)= rs * rg_thetag(rg_i) / (exp(exp_argt) - 1.0d0)
  end do

  !-----------------------------------------------------------
  ! 5) vt relaxation source: qtv
  !    qtv = sum_i rho_i (e_v_eq - e_v) / tau_i
  !-----------------------------------------------------------
  rg_qtv = 0.0d0
  do rg_i = 1, 3
     if (rg_tssum(rg_i) .gt. 0.0d0) then
        rg_qtv = rg_qtv + rg_r(rg_i) * (rg_ev_eq(rg_i) - rg_ev(rg_i)) / rg_tssum(rg_i)
     end if
  end do

  !-----------------------------------------------------------
  ! 6) reactive vibrational source: qw
  !    qw = sum_i (dw_i * e_v(tv))  (only vibrating species)
  !-----------------------------------------------------------
  rg_qw = 0.0d0
   do rg_i = 1, 3
      rg_qw = rg_qw + rg_dw_s(rg_i) * rg_ev(rg_i)
   end do

  !-----------------------------------------------------------
  ! 7) assemble source terms for conservative variables
  !
  !   ρe  = total energy (includes vibrational energy)
  !   ρev = vibrational energy only
  !
  !   -> vt exchange qtv is internal transfer between ρe and ρev,
  !      so it should not be removed from ρe (to conserve total energy).
  !   -> chemistry contributes q_chem via heats of formation.
  !-----------------------------------------------------------
  source_r(1:nof_variables) = 0.0d0

  ! chemical heat release (or absorption)
  q_chem = 0.0d0
  do rg_i = 1, nof_species
     ! rg_hzero in [j/mol] -> divide by molar mass [kg/mol] -> [j/kg]
     q_chem = q_chem + rg_dw_s(rg_i) * ( -rg_hzero(rg_i) / rg_molm(rg_i) )
  end do


  ! total energy equation: chemistry only
  source_r(dimensiona+2) = q_chem

  ! vibrational energy equation: vt + reactive vibrational source
  source_r(dimensiona+3) = rg_qtv + rg_qw

  ! species equations: mass production rates
  source_r(dimensiona+4:nof_variables) = rg_dw_s(1:nof_species)





	! ------------------------------------------------------------
  ! diagonal source jacobian for implicit time stepping
  ! ------------------------------------------------------------
  src_jac_diag(1:nof_variables) = 0.0d0

  idxe  = dimensiona + 2          ! ρe
  idxev = dimensiona + 3          ! ρev
  idxy1 = dimensiona + 4          ! first species equation

  rho_min = 1.0d-20

  ! --- species diagonal terms: ds_s / d(ρ_s) ≈ -(-ω_s / ρ_s) = ω_s / ρ_s ---
  do rg_i = 1, nof_species
     if (rg_r(rg_i) > rho_min .and. rg_dw_s(rg_i) < 0.0d0) then
        ! local destruction timescale τ_s = -ρ_s / ω_s  (ω_s < 0)
        lambda_s = -rg_dw_s(rg_i) / rg_r(rg_i)   ! [1/s]  > 0
        src_jac_diag(idxy1 + rg_i - 1) = -lambda_s
     else
        src_jac_diag(idxy1 + rg_i - 1) = 0.0d0
     end if
  end do

  ! --- vibrational energy diagonal: ds_ev / d(ρev) ≈ -1/τ_v_mix ---
  ! use a simple mixture vt timescale based on tau_tot (or tau_mw) for n2,o2,no
  tau_v_mix = 1.0d30
  sumy      = 0.0d0

  do rg_i = 1, 3   ! vibrating species: n2, o2, no
     if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0d29) then
        yi       = rg_r(rg_i) / leftv(1)          ! y_i ≈ ρ_i/ρ
        sumy     = sumy + yi
        tau_v_mix = min(tau_v_mix, tau_tot(rg_i)) ! use most restrictive vt time
     end if
  end do

  if (tau_v_mix < 1.0d29) then
     src_jac_diag(idxev) = -1.0d0 / tau_v_mix
  else
     src_jac_diag(idxev) = 0.0d0
  end if

  ! --- (optional) total energy diagonal: ds_e / d(ρe) ≈ -1/τ_echem ---
  rhoe_local = u_c_val(1,idxe,i)
  if (rhoe_local > 1.0d-12 .and. q_chem < 0.0d0) then
     lambda_e = -q_chem / rhoe_local   ! [1/s]
     src_jac_diag(idxe) = -lambda_e
  else
     src_jac_diag(idxe) = 0.0d0
  end if

  ! if chemistry/vt are turned off, also zero the jacobian
!    src_jac_diag(:) = 0.0d0


  !point implicit treatment of the

!   if (rungekutta.eq.3)then
!
!   dt_loc = dt
!
!   else
!
!    dt_loc = ielem_dtl(i)
!
!    end if
!
!
!    ! ----------------------------
!    ! (a) total energy: chemistry only (explicit, same as before)
!    ! ----------------------------
!    ! original pathway:
!    !   rhs_e -= q_chem * v
!    !   rhoe  -= dt * (rhs_e / v)
!    ! ==> rhoe += dt * q_chem
!    !
!    rhoe_old = u_c_val(1,dimensiona+2,i)
!    rhoe_new = rhoe_old + dt_loc * q_chem
!    u_c_val(1,dimensiona+2,i) = rhoe_new
!
!    ! ----------------------------
!    ! (b) vibrational energy: vt relaxation + reactive vib source
!    ! ----------------------------
!    ! use be-stable form:
!    !   d(rhoev)/dt = (rhoev_eq(t) - rhoev)/tau_v_eff + qw
!    ! relaxation part implicit, qw explicit.
!    !
!
!
!    rhoev_old = u_c_val(1,dimensiona+3,i)
!
! ! equilibrium mixture vib energy at t_tr:
! rhoev_eq = 0.0d0
! do rg_i = 1, 3
!   rhoev_eq = rhoev_eq + rg_r(rg_i) * rg_ev_eq(rg_i)
! end do
!
! ! stiffness kvt:
! kvt = 0.0d0
! do rg_i = 1, 3
!   if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0d29) then
!      kvt = kvt + rg_r(rg_i) / tau_tot(rg_i)
!   end if
! end do
!
! if (kvt > 0.0d0) then
!   alpha_v   = dt_loc * kvt
!   rhoev_new = (rhoev_old + dt_loc*(kvt*rhoev_eq + rg_qw)) / (1.0d0 + alpha_v)
! else
!   rhoev_new = rhoev_old + dt_loc * rg_qw   ! no vt exchange, only reactive vib
! end if
!
! if (rhoev_new < 0.0d0) rhoev_new = 0.0d0
! u_c_val(1,dimensiona+3,i) = rhoev_new
!
!
!
!
!
!
!
! 				! --- point-implicit chemistry step with n2 (k=1) as closure species ---
!
! 				! get local density (adjust index to your layout)
! 				rho = u_c_val(1,1,i)    ! <- check this index in your code
!
! 				! 1) update species k = 2..nof_species with patankar/be
! 				sumrhoy = 0.0d0
!
! 				do k = 2, nof_species
! 				idx = dimensiona + 3 + k
!
! 				rhoy_old = u_c_val(1,idx,i)
!
! 				prod = max(rg_dw_s(k), 0.0d0)
! 				dest = max(-rg_dw_s(k), 0.0d0)
!
! 				lambda_k = dest / max(rhoy_old, rho_min)
!
! 				rhoy_new = (rhoy_old + dt_loc * prod) / (1.0d0 + dt_loc * lambda_k)
!
! 				if (rhoy_new < 0.0d0) rhoy_new = 0.0d0
!
! 				u_c_val(1,idx,i) = rhoy_new
! 				sumrhoy = sumrhoy + rhoy_new
! 				end do
!
! 				! 2) closure: n2 (k=1) accumulates the difference so that sum rhoy = rho
! 				idxn2 = dimensiona + 3 + 1
!
! 				rhoy_n2 = rho - sumrhoy
!
! 				! optional safety: avoid tiny negative due to roundoff
! 				if (rhoy_n2 < 0.0d0) then
! 				! you can choose: clip to enforce a floor
! 				rhoy_n2 = max(rhoy_n2, 0.0d0)
! 				end if
!
! 				u_c_val(1,idxn2,i) = rhoy_n2




























end subroutine sources_realgas





subroutine sources_realgas_pi(n,iconsidered)
  implicit none
  !> sources computation in 2d/3d for 5-species, 2-t real gas
#ifdef gpu
!$omp declare target
#endif

  integer,intent(in) :: n, iconsidered
  real,dimension(1:nof_variables)::source_r, src_jac_diag

  integer :: i, j, k, l, rg_i, rg_j
  real,dimension(nof_species) :: rg_r, rg_rm, rg_n      ! rg_r = ρ_s
  real :: rg_dq_rad, rg_qtv, rg_qw, velx, vely, velz
  real,dimension(nof_species) :: rg_ev, rg_esv
  real,dimension(3) :: rg_tvs
  real,dimension(nof_species) :: rg_dw_s
  real :: rg_ke1, rg_ke2, rg_ke3, rg_ke4, rg_ke5
  real :: rg_t, rg_ta, rg_tv
  real :: rg_z, rg_speed
  real :: rg_r1, rg_r2, rg_r3, rg_r4, rg_r5
  real :: rg_kf1, rg_kf2, rg_kf3, rg_kf4, rg_kf5
  real :: rg_kb1, rg_kb2, rg_kb3, rg_kb4, rg_kb5

  real :: rg_kf11,rg_kf12,rg_kf13,rg_kf14,rg_kf15
  real :: rg_kf21,rg_kf22,rg_kf23,rg_kf24,rg_kf25
  real :: rg_kf31,rg_kf32,rg_kf33,rg_kf34,rg_kf35
  real :: rg_kb11,rg_kb12,rg_kb13,rg_kb14,rg_kb15
  real :: rg_kb21,rg_kb22,rg_kb23,rg_kb24,rg_kb25
  real :: rg_kb31,rg_kb32,rg_kb33,rg_kb34,rg_kb35

  real :: kf11m,kf12m,kf13m,kf14m,kf15m
  real :: kf21m,kf22m,kf23m,kf24m,kf25m
  real :: kf31m,kf32m,kf33m,kf34m,kf35m
  real :: kb11m,kb12m,kb13m,kb14m,kb15m
  real :: kb21m,kb22m,kb23m,kb24m,kb25m
  real :: kb31m,kb32m,kb33m,kb34m,kb35m
  real :: kf4m,kf5m,kb4m,kb5m
  real :: kf1_eff,kb1_eff,kf2_eff,kb2_eff,kf3_eff,kb3_eff
  real :: rg_c(5)            ! molar concentrations [mol/m^3]
  real :: wdot(5)            ! reaction rates [mol/(m^3 s)]
  integer :: nu(5,5)         ! stoichiometric matrix ν(s,r)
  integer :: s, r

  real,dimension(nof_species,nof_species) :: rg_mu
  real :: temp1, temp2,taumw,taupk
  real,dimension(nof_species) :: rg_tsmw, rg_tsp, rg_tssum, rg_ablot, rg_ev_eq
  real,dimension(1:nof_variables) :: leftv, tempvect,tempvx
  real :: mp_pinfl, gammal,q_chem
  real :: p_atm, t_13, rg_pressure, n_tot, exp_argtv, exp_argt, rs
  integer :: rg_molx
  real,dimension(1:nof_species) :: rg_mass_amu

  real, parameter :: kb = 1.380649d-23     ! boltzmann [j/k]
  real, parameter :: na = 6.02214076d23    ! avogadro [1/mol]

  real,dimension(1:nof_species) :: tau_mw    !! millikan-white [s]
  real,dimension(1:nof_species) :: tau_park  ! park [s]
  real,dimension(1:nof_species) :: tau_tot   ! combined [s]
  real :: num, den,rggtf
  real :: mu_sr, a_sr, b_sr
  real :: log_p_tau, p_tau, tau_sr
  real :: sigma0, sigma_v, v_th, m_s,maxloss
  real :: rg_alpha, rg_alpha_i, rg_rho_min,sum_sourcex
  integer :: idxe, idxev, idxy1,idxn2,idx
 real    :: rho_min, lambda_s, lambda_e
 real    :: tau_v_mix, rhoe_local, yi, sumy,q_ttv,tv_mix,t_eff
  real :: dt_loc
   real :: rhoe_old, rhoe_new
   real :: rhoev_old, rhoev_new, rhoev_eq
   real :: alpha_v
   real :: prod, dest, lambda_k, rhoy_old, rhoy_new
   real :: tau_v_eff,kvt,rhoy_n2,rho,sumrhoy
   real :: rho_vib, b_v, lambda_v


  ! ---------------------------------------------------------------------------
  ! start
  ! ---------------------------------------------------------------------------
  i = iconsidered

  source_r(1:nof_variables) = 0.0d0
  tau_mw(:)   = 1.0d30
  tau_park(:) = 1.0d30
  tau_tot(:)  = 1.0d30

  !firstly the base values








  !-----------------------------------------------------------
  ! get conservative variables and primitive variables
  !-----------------------------------------------------------
  leftv(1:nof_variables) = u_c_val(1,1:nof_variables,i)

  ! species densities ρ_s = ρ y_s
  do rg_i = 1, nof_species
     rg_r(rg_i) = leftv(dimensiona+3+rg_i)
  end do

  !-----------------------------------------------------------
  ! 1) molar concentrations c_s = rho_s / m_s
  !-----------------------------------------------------------
  do s = 1, 5
     rg_c(s) = rg_r(s) / rg_molm(s)
  end do

  call cons2prim(n,leftv,mp_pinfl,gammal)

  rg_pressure = leftv(dimensiona+2)
  p_atm       = rg_pressure / rgs_pa_per_atm

  ! molar masses in kg/mol -> "amu-like" for mw reduced-mass formula
  rg_mass_amu = rg_molm * 1000.0d0

  tempvect(1:nof_variables) = u_c_val(1,1:nof_variables,i)
  call cons2div(n,tempvect,mp_pinfl,gammal)

  ! leftv:   (ρ, u, v, ρe, ρev, ρy1..ρy5)
  ! tempvect:(ρ, u, v, ttr, tv, y1..y5)

  ! temperatures
  rg_t  = tempvect(dimensiona+2)   ! ttr
  rg_tv = tempvect(dimensiona+3)   ! tv
  rg_z  = 1.0d0 / rg_t
  t_13  = rg_t**(-1.0d0/3.0d0)

  ! adjusted two-temperature model (candler/park style)
  rg_ta = rg_t**0.6d0 * rg_tv**0.4d0

  ! flow speed (not used in source terms yet)
  velx = leftv(2)
  vely = leftv(3)
  velz = 0.0d0
  if (dimensiona.eq.3) then
     velz = leftv(4)
  end if
  rg_speed = sqrt(velx**2 + vely**2 + velz**2)

  !-----------------------------------------------------------
  ! chemical source: reaction rates (candler 5-reaction model)
  ! using ρ-based rate coefficients cf [m^3/(kg s)]
  !-----------------------------------------------------------
  rg_dw_s(:) = 0.0d0



if (rg_nof_reactions == 5) then

  !
  ! candler 5-reaction model (your original code)
  !

  ! equilibrium constants k_e(t) (candler curve fits, t = ttr)
  rg_ke1 = exp(3.898d0 -12.611d0*rg_z +0.683d0*rg_z**2 -0.118d0*rg_z**3 +0.006d0*rg_z**4)
  rg_ke2 = exp(1.335d0 - 4.127d0*rg_z -0.616d0*rg_z**2 +0.093d0*rg_z**3 -0.005d0*rg_z**4)
  rg_ke3 = exp(1.549d0 - 7.784d0*rg_z +0.228d0*rg_z**2 -0.043d0*rg_z**3 +0.002d0*rg_z**4)
  rg_ke4 = exp(2.349d0 - 4.828d0*rg_z +0.455d0*rg_z**2 -0.075d0*rg_z**3 +0.004d0*rg_z**4)
  rg_ke5 = exp(0.215d0 - 3.652d0*rg_z +0.843d0*rg_z**2 -0.136d0*rg_z**3 +0.007d0*rg_z**4)

  ! --- reaction 1: n2 + m <-> 2n + m (in m^3/mol/s) ---
  kf11m = 3.78d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n2
  kf12m = 3.78d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o2
  kf13m = 3.78d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = no
  kf14m = 1.11d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n
  kf15m = 1.11d18 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o

  kb11m = kf11m / rg_ke1
  kb12m = kf12m / rg_ke1
  kb13m = kf13m / rg_ke1
  kb14m = kf14m / rg_ke1
  kb15m = kf15m / rg_ke1

  kf1_eff = kf11m*rg_c(1) + kf12m*rg_c(2) + kf13m*rg_c(3) &
          + kf14m*rg_c(4) + kf15m*rg_c(5)
  kb1_eff = kb11m*rg_c(1) + kb12m*rg_c(2) + kb13m*rg_c(3) &
          + kb14m*rg_c(4) + kb15m*rg_c(5)

  wdot(1) = kf1_eff*rg_c(1) - kb1_eff*rg_c(4)**2   ! n2 + m <-> 2n + m

  ! --- reaction 2: o2 + m <-> 2o + m ---
  kf21m = 2.75d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = n2
  kf22m = 2.75d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = o2
  kf23m = 2.75d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = no
  kf24m = 8.25d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = n
  kf25m = 8.25d16 * 1.0d-6 * rg_ta**(-1.0d0) * exp(-5.95d4/rg_ta)   ! m = o

  kb21m = kf21m / rg_ke2
  kb22m = kf22m / rg_ke2
  kb23m = kf23m / rg_ke2
  kb24m = kf24m / rg_ke2
  kb25m = kf25m / rg_ke2

  kf2_eff = kf21m*rg_c(1) + kf22m*rg_c(2) + kf23m*rg_c(3) &
          + kf24m*rg_c(4) + kf25m*rg_c(5)
  kb2_eff = kb21m*rg_c(1) + kb22m*rg_c(2) + kb23m*rg_c(3) &
          + kb24m*rg_c(4) + kb25m*rg_c(5)

  wdot(2) = kf2_eff*rg_c(2) - kb2_eff*rg_c(5)**2   ! o2 + m <-> 2o + m

  ! --- reaction 3: no + m <-> n + o + m ---
  kf31m = 2.30d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = n2
  kf32m = 2.30d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = o2
  kf33m = 2.30d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = no
  kf34m = 4.60d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = n
  kf35m = 4.60d14 * 1.0d-6 * rg_ta**(-0.5d0) * exp(-7.55d4/rg_ta)   ! m = o

  kb31m = kf31m / rg_ke3
  kb32m = kf32m / rg_ke3
  kb33m = kf33m / rg_ke3
  kb34m = kf34m / rg_ke3
  kb35m = kf35m / rg_ke3

  kf3_eff = kf31m*rg_c(1) + kf32m*rg_c(2) + kf33m*rg_c(3) &
          + kf34m*rg_c(4) + kf35m*rg_c(5)
  kb3_eff = kb31m*rg_c(1) + kb32m*rg_c(2) + kb33m*rg_c(3) &
          + kb34m*rg_c(4) + kb35m*rg_c(5)

  wdot(3) = kf3_eff*rg_c(3) - kb3_eff*rg_c(4)*rg_c(5) ! no + m <-> n + o + m

  ! --- reaction 4: n2 + o <-> no + n (bimolecular, 2-body) ---
  kf4m = 3.18d10 * 1.0d-6 * rg_ta**(0.1d0) * exp(-3.77d4/rg_ta)
  kb4m = kf4m / rg_ke4

  wdot(4) = kf4m*rg_c(1)*rg_c(5) - kb4m*rg_c(3)*rg_c(4)

  ! --- reaction 5: no + o <-> o2 + n (bimolecular) ---
  kf5m = 2.16d5 * 1.0d-6 * rg_ta**(1.29d0) * exp(-1.922d4/rg_ta)
  kb5m = kf5m / rg_ke5

  wdot(5) = kf5m*rg_c(3)*rg_c(5) - kb5m*rg_c(2)*rg_c(4)

  ! build stoichiometry & mass production (unchanged)
  do s = 1,5
     do r = 1,5
        nu(s,r) = 0
     end do
  end do

  ! r1: n2 -> 2n
  nu(1,1) = -1
  nu(4,1) = +2

  ! r2: o2 -> 2o
  nu(2,2) = -1
  nu(5,2) = +2

  ! r3: no -> n + o
  nu(3,3) = -1
  nu(4,3) = +1
  nu(5,3) = +1

  ! r4: n2 + o -> no + n
  nu(1,4) = -1
  nu(5,4) = -1
  nu(3,4) = +1
  nu(4,4) = +1

  ! r5: no + o -> o2 + n
  nu(3,5) = -1
  nu(5,5) = -1
  nu(2,5) = +1
  nu(4,5) = +1

  rg_dw_s(:) = 0.0d0
  do r = 1,5
     do s = 1,5
        rg_dw_s(s) = rg_dw_s(s) + nu(s,r) * rg_molm(s) * wdot(r)
     end do
  end do


end if




if (rg_nof_reactions == 6) then

  !
  ! park 2-t “6-reaction” model (option a, 5 neutral reactions)
  ! same topology as candler (r1..r5), park coefficients + 3rd bodies
  !

  ! equilibrium constants (keep your candler curve-fits)
  rg_ke1 = exp(3.898d0 -12.611d0*rg_z +0.683d0*rg_z**2 -0.118d0*rg_z**3 +0.006d0*rg_z**4)
  rg_ke2 = exp(1.335d0 - 4.127d0*rg_z -0.616d0*rg_z**2 +0.093d0*rg_z**3 -0.005d0*rg_z**4)
  rg_ke3 = exp(1.549d0 - 7.784d0*rg_z +0.228d0*rg_z**2 -0.043d0*rg_z**3 +0.002d0*rg_z**4)
  rg_ke4 = exp(2.349d0 - 4.828d0*rg_z +0.455d0*rg_z**2 -0.075d0*rg_z**3 +0.004d0*rg_z**4)
  rg_ke5 = exp(0.215d0 - 3.652d0*rg_z +0.843d0*rg_z**2 -0.136d0*rg_z**3 +0.007d0*rg_z**4)

  ! ---- r1: n2 + m <-> 2n + m (park 2t, table 1) ----
  !   kf = c_m t_a^{-1.6} exp(-113200/t_a), cm^3/(mol s)
  !   m = n2,o2,no : c = 7.0e21
  !   m = n,o      : c = 3.0e22
  kf11m = 7.0d21 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n2
  kf12m = 7.0d21 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o2
  kf13m = 7.0d21 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = no
  kf14m = 3.0d22 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = n
  kf15m = 3.0d22 * 1.0d-6 * rg_ta**(-1.6d0) * exp(-1.132d5/rg_ta)  ! m = o

  kb11m = kf11m / rg_ke1
  kb12m = kf12m / rg_ke1
  kb13m = kf13m / rg_ke1
  kb14m = kf14m / rg_ke1
  kb15m = kf15m / rg_ke1

  kf1_eff = kf11m*rg_c(1) + kf12m*rg_c(2) + kf13m*rg_c(3) &
          + kf14m*rg_c(4) + kf15m*rg_c(5)
  kb1_eff = kb11m*rg_c(1) + kb12m*rg_c(2) + kb13m*rg_c(3) &
          + kb14m*rg_c(4) + kb15m*rg_c(5)

  wdot(1) = kf1_eff*rg_c(1) - kb1_eff*rg_c(4)**2

  ! ---- r2: o2 + m <-> 2o + m (park 2t, table 1) ----
  !   kf = c_m t_a^{-1.5} exp(-59500/t_a), cm^3/(mol s)
  !   m = n2,o2,no : c = 2.0e21
  !   m = n,o      : c = 1.0e22
  kf21m = 2.0d21 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = n2
  kf22m = 2.0d21 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = o2
  kf23m = 2.0d21 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = no
  kf24m = 1.0d22 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = n
  kf25m = 1.0d22 * 1.0d-6 * rg_ta**(-1.5d0) * exp(-5.95d4/rg_ta)   ! m = o

  kb21m = kf21m / rg_ke2
  kb22m = kf22m / rg_ke2
  kb23m = kf23m / rg_ke2
  kb24m = kf24m / rg_ke2
  kb25m = kf25m / rg_ke2

  kf2_eff = kf21m*rg_c(1) + kf22m*rg_c(2) + kf23m*rg_c(3) &
          + kf24m*rg_c(4) + kf25m*rg_c(5)
  kb2_eff = kb21m*rg_c(1) + kb22m*rg_c(2) + kb23m*rg_c(3) &
          + kb24m*rg_c(4) + kb25m*rg_c(5)

  wdot(2) = kf2_eff*rg_c(2) - kb2_eff*rg_c(5)**2

  ! ---- r3: no + m <-> n + o + m (park/koshi) ----
  !   kf ≈ c_m t_a^{0} exp(-75500/t_a), cm^3/(mol s)
  !   m = no,o,n : c = 1.1e17
  !   m = n2,o2 : c = (1/22)*1.1e17  [koshi]
  kf31m = (1.1d17/22.0d0) * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = n2
  kf32m = (1.1d17/22.0d0) * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = o2
  kf33m =  1.1d17        * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = no
  kf34m =  1.1d17        * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = n
  kf35m =  1.1d17        * 1.0d-6 * exp(-7.55d4/rg_ta)   ! m = o

  kb31m = kf31m / rg_ke3
  kb32m = kf32m / rg_ke3
  kb33m = kf33m / rg_ke3
  kb34m = kf34m / rg_ke3
  kb35m = kf35m / rg_ke3

  kf3_eff = kf31m*rg_c(1) + kf32m*rg_c(2) + kf33m*rg_c(3) &
          + kf34m*rg_c(4) + kf35m*rg_c(5)
  kb3_eff = kb31m*rg_c(1) + kb32m*rg_c(2) + kb33m*rg_c(3) &
          + kb34m*rg_c(4) + kb35m*rg_c(5)

  wdot(3) = kf3_eff*rg_c(3) - kb3_eff*rg_c(4)*rg_c(5)

  ! ---- r4: n2 + o <-> no + n (park exchange, table 1) ----
  !   kf = 6.4e7 t^{-1} exp(-38370/t), cm^3/(mol s)
  kf4m = 6.4d7 * 1.0d-6 * rg_t**(-1.0d0) * exp(-3.837d4/rg_t)
  kb4m = kf4m / rg_ke4

  wdot(4) = kf4m*rg_c(1)*rg_c(5) - kb4m*rg_c(3)*rg_c(4)

  ! ---- r5: no + o <-> o2 + n (park exchange, table 1) ----
  !   kr_exo = 8.4e12 exp(-19450/t), cm^3/(mol s)
  !   take kf endothermic = kr_exo / k_eq
  kf5m = 8.4d12 * 1.0d-6 * exp(-1.945d4/rg_t)
  kb5m = kf5m / rg_ke5

  wdot(5) = kf5m*rg_c(3)*rg_c(5) - kb5m*rg_c(2)*rg_c(4)

  ! stoichiometry is the same as candler (same 5 reactions)
  do s = 1,5
     do r = 1,5
        nu(s,r) = 0
     end do
  end do

  ! r1: n2 -> 2n
  nu(1,1) = -1
  nu(4,1) = +2

  ! r2: o2 -> 2o
  nu(2,2) = -1
  nu(5,2) = +2

  ! r3: no -> n + o
  nu(3,3) = -1
  nu(4,3) = +1
  nu(5,3) = +1

  ! r4: n2 + o -> no + n
  nu(1,4) = -1
  nu(5,4) = -1
  nu(3,4) = +1
  nu(4,4) = +1

  ! r5: no + o -> o2 + n
  nu(3,5) = -1
  nu(5,5) = -1
  nu(2,5) = +1
  nu(4,5) = +1

  rg_dw_s(:) = 0.0d0
  do r = 1,5
     do s = 1,5
        rg_dw_s(s) = rg_dw_s(s) + nu(s,r) * rg_molm(s) * wdot(r)
     end do
  end do


end if










  rg_rho_min = 1.0d-20    ! or whatever you use as "zero" density
  rg_alpha   = 1.0d0

	! first: compute most restrictive scaling factor alpha
	do rg_i = 1, nof_species

		! if density is essentially zero, don't allow further destruction
		if (rg_r(rg_i) <= rg_rho_min .and. rg_dw_s(rg_i) < 0.0d0) then
			rg_dw_s(rg_i) = 0.0d0
			cycle
		end if

		! only destruction can cause negativity
! 		if (rg_dw_s(rg_i) < 0.0d0) then
! 			! need: rho_i + dt * alpha * dw_i >= 0
! 			!  => alpha <= rho_i / (-dt * dw_i)
! 			rg_alpha_i = rg_r(rg_i) / (-ielem_dtl(iconsidered) * rg_dw_s(rg_i))
!
! 			if (rg_alpha_i < 1.0d-2) rg_alpha = 1.0d-2
! 		end if
	end do

	! clamp alpha to [0,1]
	if (rg_alpha > 1.0d0) rg_alpha = 1.0d0
	if (rg_alpha < 0.0d0) rg_alpha = 0.0d0

	! second: apply scaling if needed
	if (rg_alpha < 1.0d0) then
		do rg_i = 1, nof_species
			rg_dw_s(rg_i) = rg_alpha * rg_dw_s(rg_i)
		end do
	end if

	! optional tiny cutoff (for cleanliness, won’t break stoichiometry)
	do rg_i = 1, nof_species
		if (abs(rg_dw_s(rg_i)) < 1.0d-50) rg_dw_s(rg_i) = 0.0d0
	end do






  !-----------------------------------------------------------
  ! 1) millikan-white vt relaxation times (species 1..3)
  !-----------------------------------------------------------
  do rg_i = 1, 3
     num = 0.0d0
     den = 0.0d0

     do rg_j = 1, nof_species
        if (rg_r(rg_j) <= 0.0d0) cycle

        ! reduced mass μ_sr in amu
        mu_sr = rg_mass_amu(rg_i)*rg_mass_amu(rg_j) / (rg_mass_amu(rg_i)+rg_mass_amu(rg_j))

        ! mw coefficients
        a_sr = 0.00116d0 * sqrt(mu_sr) * rg_thetag(rg_i)**(4.0d0/3.0d0)
        b_sr = 0.015d0   * mu_sr**0.25d0

        ! p * tau_sr in atm*s
        log_p_tau = a_sr*(t_13 - b_sr) - 18.42d0
        p_tau       = exp(log_p_tau)

        ! tau_sr [s]: p_tau has atm*s, divide by p in atm
        tau_sr = p_tau / p_atm

        ! mixture rule with mole fractions:
        ! 1/tau_s = σ_r (x_r / tau_sr)
        ! x_r ~ n_r / σ n_r, and n_r ~ ρ_r / m_r
        num = num + (rg_r(rg_j)/rg_molm(rg_j)) / tau_sr
        den = den + (rg_r(rg_j)/rg_molm(rg_j))
     end do

     if (num > 0.0d0 .and. den > 0.0d0) then
        tau_mw(rg_i) = den / num
     else
        tau_mw(rg_i) = 1.0d30
     end if



  end do




  sigma0 = 3.0d-21    ! [m^2]




! number density [1/m^3] based on translational/rotational t
n_tot = rg_pressure / (kb * rg_t)   ! rg_t = t_tr

! single vibrational temperature for the mixture
         ! <-- your global vibrational temperature

! avoid non-physical tv
tv_mix = max(rg_tv, 300.0d0)

! park-style weighting between t and tv (0.7 is a common choice)
q_ttv = 0.7d0

! effective two-temperature average (same for all species)
! geometric mean: t_eff = t^q * tv^(1-q)
t_eff = rg_t**q_ttv * tv_mix**(1.0d0 - q_ttv)
t_eff = max(t_eff, 300.0d0)

do rg_i = 1, 3
   ! particle mass [kg/particle]
   m_s = rg_molm(rg_i) / na

   ! effective cross-section using t_eff instead of rg_t
   sigma_v = sigma0 * (50000.0d0 / t_eff)**2

   ! thermal speed using t_eff
   v_th = sqrt( 8.0d0 * kb * t_eff / (pi * m_s) )

   ! tau_park [s]
   tau_park(rg_i) = 1.0d0 / ( n_tot * sigma_v * v_th )

end do
























  !-----------------------------------------------------------
  ! 3) combine millikan-white + park
  !    1/tau_tot = 1/tau_mw + 1/tau_park
  !-----------------------------------------------------------


 do rg_i = 1, 3

  ! local copies
  taumw  = tau_mw(rg_i)
  taupk  = tau_park(rg_i)

  ! replace non-positive with 'disabled'
  if (taumw <= 0.0d0) taumw = 1.0d30
  if (taupk <= 0.0d0) taupk = 1.0d30

  ! now combine them safely
  if (taumw >= 1.0d29 .and. taupk >= 1.0d29) then
     ! both basically 'off'
     tau_tot(rg_i) = 1.0d30

  elseif (taumw >= 1.0d29) then
     ! only park valid
     tau_tot(rg_i) = taupk

  elseif (taupk >= 1.0d29) then
     ! only mw valid
     tau_tot(rg_i) = taumw

  else
     ! both valid: use harmonic sum
      rggtf = taumw / taupk
       rggtf = max(1.0d0, rggtf)
       rggtf = min(10d0, rggtf)
       tau_tot(rg_i) = taumw / rggtf
!         tau_tot(rg_i) = 1.0d0 / ( 1.0d0/taumw + 1.0d0/taupk )

  end if


  rg_tssum(rg_i) = tau_tot(rg_i)

end do

  !-----------------------------------------------------------
  ! 4) vibrational energies e_v(tv) and e_v_eq(t)
  !-----------------------------------------------------------
  do rg_i = 1, 3
     rs         = rgs_ru / rg_molm(rg_i)
     exp_argtv  = rg_thetag(rg_i) / rg_tv
     rg_ev(rg_i) = rs * rg_thetag(rg_i) / (exp(exp_argtv) - 1.0d0)

     exp_argt      = rg_thetag(rg_i) / rg_t
     rg_ev_eq(rg_i)= rs * rg_thetag(rg_i) / (exp(exp_argt) - 1.0d0)
  end do

  !-----------------------------------------------------------
  ! 5) vt relaxation source: qtv
  !    qtv = sum_i rho_i (e_v_eq - e_v) / tau_i
  !-----------------------------------------------------------
  rg_qtv = 0.0d0
  do rg_i = 1, 3
     if (rg_tssum(rg_i) .gt. 0.0d0) then
        rg_qtv = rg_qtv + rg_r(rg_i) * (rg_ev_eq(rg_i) - rg_ev(rg_i)) / rg_tssum(rg_i)
     end if
  end do

  !-----------------------------------------------------------
  ! 6) reactive vibrational source: qw
  !    qw = sum_i (dw_i * e_v(tv))  (only vibrating species)
  !-----------------------------------------------------------
  rg_qw = 0.0d0
   do rg_i = 1, 3
      rg_qw = rg_qw + rg_dw_s(rg_i) * rg_ev(rg_i)
   end do

  !-----------------------------------------------------------
  ! 7) assemble source terms for conservative variables
  !
  !   ρe  = total energy (includes vibrational energy)
  !   ρev = vibrational energy only
  !
  !   -> vt exchange qtv is internal transfer between ρe and ρev,
  !      so it should not be removed from ρe (to conserve total energy).
  !   -> chemistry contributes q_chem via heats of formation.
  !-----------------------------------------------------------
  source_r(1:nof_variables) = 0.0d0

  ! chemical heat release (or absorption)
  q_chem = 0.0d0
  do rg_i = 1, nof_species
     ! rg_hzero in [j/mol] -> divide by molar mass [kg/mol] -> [j/kg]
     q_chem = q_chem + rg_dw_s(rg_i) * ( -rg_hzero(rg_i) / rg_molm(rg_i) )
  end do


  ! total energy equation: chemistry only
  source_r(dimensiona+2) = q_chem

  ! vibrational energy equation: vt + reactive vibrational source
  source_r(dimensiona+3) = rg_qtv + rg_qw

  ! species equations: mass production rates
  source_r(dimensiona+4:nof_variables) = rg_dw_s(1:nof_species)





	! ------------------------------------------------------------
  ! diagonal source jacobian for implicit time stepping
  ! ------------------------------------------------------------
  src_jac_diag(1:nof_variables) = 0.0d0

  idxe  = dimensiona + 2          ! ρe
  idxev = dimensiona + 3          ! ρev
  idxy1 = dimensiona + 4          ! first species equation

  rho_min = 1.0d-20

  ! --- species diagonal terms: ds_s / d(ρ_s) ≈ -(-ω_s / ρ_s) = ω_s / ρ_s ---
  do rg_i = 1, nof_species
     if (rg_r(rg_i) > rho_min .and. rg_dw_s(rg_i) < 0.0d0) then
        ! local destruction timescale τ_s = -ρ_s / ω_s  (ω_s < 0)
        lambda_s = -rg_dw_s(rg_i) / rg_r(rg_i)   ! [1/s]  > 0
        src_jac_diag(idxy1 + rg_i - 1) = -lambda_s
     else
        src_jac_diag(idxy1 + rg_i - 1) = 0.0d0
     end if
  end do

  ! --- vibrational energy diagonal: ds_ev / d(ρev) ≈ -1/τ_v_mix ---
  ! use a simple mixture vt timescale based on tau_tot (or tau_mw) for n2,o2,no
  tau_v_mix = 1.0d30
  sumy      = 0.0d0

  do rg_i = 1, 3   ! vibrating species: n2, o2, no
     if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0d29) then
        yi       = rg_r(rg_i) / leftv(1)          ! y_i ≈ ρ_i/ρ
        sumy     = sumy + yi
        tau_v_mix = min(tau_v_mix, tau_tot(rg_i)) ! use most restrictive vt time
     end if
  end do

  if (tau_v_mix < 1.0d29) then
     src_jac_diag(idxev) = -1.0d0 / tau_v_mix
  else
     src_jac_diag(idxev) = 0.0d0
  end if

  ! --- (optional) total energy diagonal: ds_e / d(ρe) ≈ -1/τ_echem ---
  rhoe_local = u_c_val(1,idxe,i)
  if (rhoe_local > 1.0d-12 .and. q_chem < 0.0d0) then
     lambda_e = -q_chem / rhoe_local   ! [1/s]
     src_jac_diag(idxe) = -lambda_e
  else
     src_jac_diag(idxe) = 0.0d0
  end if

  ! if chemistry/vt are turned off, also zero the jacobian
!    src_jac_diag(:) = 0.0d0


  !point implicit treatment of the






  if (rungekutta.eq.3)then

  dt_loc = dt

  else

   dt_loc = ielem_dtl(i)

   end if










   ! ----------------------------
   ! (a) total energy: chemistry only (explicit, same as before)
   ! ----------------------------
   ! original pathway:
   !   rhs_e -= q_chem * v
   !   rhoe  -= dt * (rhs_e / v)
   ! ==> rhoe += dt * q_chem
   !
   rhoe_old = u_c_val(1,dimensiona+2,i)
   rhoe_new = rhoe_old + dt_loc * q_chem
   u_c_val(1,dimensiona+2,i) = rhoe_new

   ! ----------------------------
   ! (b) vibrational energy: vt relaxation + reactive vib source
   ! ----------------------------
   ! use be-stable form:
   !   d(rhoev)/dt = (rhoev_eq(t) - rhoev)/tau_v_eff + qw
   ! relaxation part implicit, qw explicit.
   !


   rhoev_old = u_c_val(1,dimensiona+3,i)

   ! Backward-Euler scalar approximation to
   !   d(rhoev)/dt = sum_s rho_s * (ev_eq_s - ev_s) / tau_s + rg_qw
   ! with the mixture vibrational specific energy approximated as
   !   ev_mix = rhoev / rho_vib,
   ! where rho_vib = sum_s rho_s over vibrating species.
   ! This gives
   !   d(rhoev)/dt = b_v - lambda_v * rhoev + rg_qw
   ! where
   !   b_v      = sum_s rho_s * ev_eq_s / tau_s
   !   lambda_v = (sum_s rho_s / tau_s) / rho_vib

   rho_vib = 0.0d0
   b_v     = 0.0d0
   kvt     = 0.0d0

   do rg_i = 1, 3
      if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0d29) then
         rho_vib = rho_vib + rg_r(rg_i)
         b_v     = b_v     + rg_r(rg_i) * rg_ev_eq(rg_i) / tau_tot(rg_i)
         kvt     = kvt     + rg_r(rg_i) / tau_tot(rg_i)
      end if
   end do

   if (rho_vib > rho_min .and. kvt > 0.0d0) then
      lambda_v  = kvt / rho_vib
      alpha_v   = dt_loc * lambda_v
      rhoev_new = (rhoev_old + dt_loc * (b_v + rg_qw)) / (1.0d0 + alpha_v)
   else
      rhoev_new = rhoev_old + dt_loc * rg_qw   ! no vt exchange, only reactive vib
   end if

   if (rhoev_new < 0.0d0) rhoev_new = 0.0d0
   u_c_val(1,dimensiona+3,i) = rhoev_new







				! --- point-implicit chemistry step with n2 (k=1) as closure species ---

				! get local density (adjust index to your layout)
				rho = u_c_val(1,1,i)    ! <- check this index in your code

				! 1) update species k = 2..nof_species with patankar/be
				sumrhoy = 0.0d0

				do k = 2, nof_species
				idx = dimensiona + 3 + k

				rhoy_old = u_c_val(1,idx,i)

				prod = max(rg_dw_s(k), 0.0d0)
				dest = max(-rg_dw_s(k), 0.0d0)

				lambda_k = dest / max(rhoy_old, rho_min)

				rhoy_new = (rhoy_old + dt_loc * prod) / (1.0d0 + dt_loc * lambda_k)

				if (rhoy_new < 0.0d0) rhoy_new = 0.0d0

				u_c_val(1,idx,i) = rhoy_new
				sumrhoy = sumrhoy + rhoy_new
				end do

				! 2) closure: n2 (k=1) accumulates the difference so that sum rhoy = rho
				idxn2 = dimensiona + 3 + 1

				rhoy_n2 = rho - sumrhoy

				! optional safety: avoid tiny negative due to roundoff
				if (rhoy_n2 < 0.0d0) then
				! you can choose: clip to enforce a floor
				rhoy_n2 = max(rhoy_n2, 0.0d0)
				end if

				u_c_val(1,idxn2,i) = rhoy_n2




























end subroutine sources_realgas_pi



! !---------------------------------------------
! ! u_base: state after fluxes (conv + visc)
! !---------------------------------------------
! rho_base   = u_c_val(1,1,i)
! rhoe_base  = u_c_val(1,dimensiona+2,i)
! rhoev_base = u_c_val(1,dimensiona+3,i)
!
! ! store base species
! do k = 1, nof_species
!   idx = dimensiona + 3 + k
!   rhoy_base(k) = u_c_val(1,idx,i)
! end do
!
! ! initial guess for post-chem state = base state
! rho    = rho_base
! rhoe   = rhoe_base
! rhoev  = rhoev_base
! do k = 1, nof_species
!   rhoy(k) = rhoy_base(k)
! end do
!
! !---------------------------------------------
! ! local picard iterations for chemistry
! !---------------------------------------------
! do it = 1, n_it   ! n_it = 2 or 3 is usually enough
!
!   ! (1) recompute thermo and source terms from current guess
!   !     - build mass fractions y_k = rhoy(k)/rho
!   !     - compute t, tv from (rho, rhoe, rhoev, rhoy)
!   !     - compute rg_r, rg_ev_eq, tau_tot, q_chem, rg_qw, rg_dw_s(k), etc.
!   call compute_thermo_and_sources_from_state( &
!        rho, rhoe, rhoev, rhoy,              &
!        rg_r, rg_ev_eq, tau_tot,            &
!        q_chem, rg_qw, rg_dw_s )
!
!   ! --- (a) total energy: explicit chemistry source, be in coefficients sense ---
!   ! we want: rhoe^{n+1} = rhoe_base + dt_loc * q_chem(u^it)
!   rhoe_new = rhoe_base + dt_loc * q_chem
!
!   ! --- (b) vibrational energy: be relax + explicit reactive part ---
!   rhoev_eq = 0.0d0
!   do rg_i = 1, 3
!     rhoev_eq = rhoev_eq + rg_r(rg_i) * rg_ev_eq(rg_i)
!   end do
!
!   kvt = 0.0d0
!   do rg_i = 1, 3
!     if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0d29) then
!        kvt = kvt + rg_r(rg_i) / tau_tot(rg_i)
!     end if
!   end do
!
!   if (kvt > 0.0d0) then
!     alpha_v   = dt_loc * kvt
!     rhoev_new = (rhoev_base + dt_loc*(kvt*rhoev_eq + rg_qw)) / (1.0d0 + alpha_v)
!   else
!     rhoev_new = rhoev_base + dt_loc * rg_qw
!   end if
!
!   if (rhoev_new < 0.0d0) rhoev_new = 0.0d0
!
!   ! --- (c) species: patankar be with n2 closure, from base state ---
!
!   sumrhoy = 0.0d0
!
!   do k = 2, nof_species
!     rhoy_old = rhoy_base(k)     ! base state in "old" term
!
!     prod = max(rg_dw_s(k), 0.0d0)
!     dest = max(-rg_dw_s(k), 0.0d0)
!
!     lambda_k = 0.0d0
!     if (rhoy_old > rho_min .and. dest > 0.0d0) then
!        lambda_k = dest / rhoy_old
!     end if
!
!     rhoy_new_k = (rhoy_old + dt_loc * prod) / (1.0d0 + dt_loc * lambda_k)
!
!     if (rhoy_new_k < 0.0d0) rhoy_new_k = 0.0d0
!
!     rhoy(k) = rhoy_new_k
!     sumrhoy = sumrhoy + rhoy_new_k
!   end do
!
!   ! closure species k=1 (n2)
!   rhoy(1) = rho_base - sumrhoy
!   if (rhoy(1) < 0.0d0) rhoy(1) = max(rhoy(1), 0.0d0)
!
!   ! (optional) renormalize rhoy(k) so that sum_k rhoy(k) = rho_base exactly
!   sumrhoytot = 0.0d0
!   do k = 1, nof_species
!     sumrhoytot = sumrhoytot + rhoy(k)
!   end do
!   if (sumrhoytot > 0.0d0) then
!     fac = rho_base / sumrhoytot
!     do k = 1, nof_species
!       rhoy(k) = rhoy(k) * fac
!     end do
!   end if
!
!   ! update iteration guess
!   rhoe  = rhoe_new
!   rhoev = rhoev_new
!
! end do   ! end chemistry iterations
!
! !---------------------------------------------
! ! write final state back to u_c
! !---------------------------------------------
! u_c_val(1,1,i)             = rho_base              ! density unchanged by chem
! u_c_val(1,dimensiona+2,i)  = rhoe
! u_c_val(1,dimensiona+3,i)  = rhoev
! do k = 1, nof_species
!   idx = dimensiona + 3 + k
!   u_c_val(1,idx,i) = rhoy(k)
! end do






subroutine sources_derivatives_computation(n)
	implicit none
!> @brief
!> sources derivative computation for implicit time stepping
	integer,intent(in)::n
	integer::i,kmaxe,iconsidered
	real,dimension(turbulenceequations)::source_t
	
	
	kmaxe=xmpielrank(n)
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, iconsidered, source_t)
#else
!$omp do
#endif
	do i=1,kmaxe
		iconsidered=i
		call sources_derivatives(n,iconsidered,source_t)
		sht(i,1:turbulenceequations)=(source_t(1:turbulenceequations)*ielem_totvolume(i))
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine sources_derivatives_computation



subroutine sources(n,iconsidered,source_t)
implicit none
!> @brief
!> sources computation procedure
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,dimension(turbulenceequations),intent(inout)::source_t
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx
real::vhx,amp,dvel,omega,squaret,tch_x,tch_x3,tch_fv1,tch_fv2
real::tch_rs,tch_r,tch_g,tch_glim,tch_fw,tch_dif,tch_dest,tch_prod
integer::i,k,j,l,ihgt,ihgj,iex, lowre
real::snorm,onorm,divnorm,ax,ay,az,tch_shh,tch_sav,verysmall,onesix,prodterm1,stild,rr
real::gg,fw,destterm,fodt,srcfull,dbpr,dbdi,dbde,dby,dbx,prodtermfinal
real:: r_des,f_des, ddw,f_des_sst,l_t_des
real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,prodterm2,sfac,sss,usss,ssss,s_bar,kron
real:: uz,vz,wx,wy,wz
real:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !for sas only
real,dimension(3,3)::vortet,tvort,svort,ovort
real,dimension(3)::vortem,vertf,vertex
real,dimension(turbulenceequations,1:3)::derivturb
!declarations for k-omega
real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
real:: sigma_k, sigma_om, f_1, f_2,phi_1,phi_2, d_omplus !-------------- those are for diffusion too!
real:: alpha_raw,alpha_star, re_t_sst,alpha_inf
real:: beta_stari, beta_i, beta_raw, beta_star
real:: k_0, om_0, wally
real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
real:: kx,ky,kz,omx,omy,omz,rho,nu,nu_tilde,difterm,prod
real,dimension(1:nof_variables)::leftv,rightv
real::mp_pinfl,gammal
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
real,dimension(1)::etvm

i=iconsidered

verysmall = 10e-16
	
vortet(1:3,1:3) = rec_grads(1:3,1:3,i)


ux = vortet(1,1);uy = vortet(1,2);uz = vortet(1,3)
vx = vortet(2,1);vy = vortet(2,2);vz = vortet(2,3)
wx = vortet(3,1);wy = vortet(3,2);wz = vortet(3,3)


do ihgt=1,3
  do ihgj=1,3
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
end do

svort=0.5*(vortet+tvort)
ovort=0.5*(vortet-tvort)





snorm=sqrt(2.0d0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+(svort(1,3)*svort(1,3))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))+(svort(2,3)*svort(2,3))+& 
	       (svort(3,1)*svort(3,1))+(svort(3,2)*svort(3,2))+(svort(3,3)*svort(3,3))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
onorm=sqrt(2.0d0*((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+(ovort(1,3)*ovort(1,3))+&
	       (ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2))+(ovort(2,3)*ovort(2,3))+&
	       (ovort(3,1)*ovort(3,1))+(ovort(3,2)*ovort(3,2))+(ovort(3,3)*ovort(3,3))))
omega=onorm


divnorm=ux+vy+wz  !careful with the sign. if it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)+(wz*wz)))&
	+((uy+vx)*(uy+vx)+(uz+wx)*(uz+wx)+(wy+vz)*(wy+vz))&
	-(2.0/3.0*(ux+vy+wz)*(ux+vy+wz)))


derivturb(1,1:3) = rec_grads(5,1:3,i)
squaret=(sqrt((derivturb(1,1)**2)+(derivturb(1,2)**2)+(derivturb(1,3)**2)))**2


leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
rightv(1:nof_variables)=leftv(1:nof_variables)

call get_visc_conduct(n,leftv,rightv,viscl,laml)




turbmv(1)=u_ct_val(1,1,i)
turbmv(2)=turbmv(1)



 select case(turbulencemodel)
 
 case(1)   !! spalart–almaras model (compressible ρν~ formulation)

    ! optional rotation/curvature correction
    if (rot_corr .eq. 1) then
        omega = omega + 2.0d0 * min(0.0d0, snorm - onorm)
    end if

    eddyfl(2) = turbmv(1)
    eddyfr(2) = turbmv(1)

    call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)

    rho = leftv(1)
    nu  = viscl(1) / rho
    nu_tilde = turbmv(1) / rho      ! convert φ = ρ ν~  →  ν~

     nu_tilde = max(nu_tilde, 1.0d-12 * visc)

    ! χ = ν~/ν
    tch_x = nu_tilde / nu
    tch_x3 = tch_x*tch_x*tch_x

    ! fv1
    tch_fv1 = tch_x3 / (tch_x3 + (cv1*cv1*cv1))

    ! fv2
    tch_fv2 = 1.0d0 - tch_x / (1.0d0 + tch_x * tch_fv1)

    ! physical wall distance
    ddw = ielem_walldist(i)
    ddw= max(ddw, 10e-12)

    ! --- des correction ---
    if (des_model .eq. 1) then
        cell_volume = ielem_totvolume(i)
        delta_cell  = cell_volume**0.333333333333333d0
        ddw = min(ddw, c_des_sa * delta_cell)
    end if

    if (des_model .eq. 2) then
        cell_volume = ielem_totvolume(i)
        delta_cell  = cell_volume**0.333333333333333d0

        r_des = min(10.0d0, (viscl(1)+viscl(3)) / (snorm*(kappa*ddw)**2 + 1.0d-16))
        f_des = 1.0d0 - tanh((8.0d0*r_des)**3)

        ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1.0d-16), 1.0d-16)
    end if

    ! --- stilde ---
    prodterm1 = (nu_tilde) / (kappa*kappa*ddw*ddw)
    stild = max( omega + tch_fv2*prodterm1 , 0.3d0 * omega)

    ! --- production term ---
    prodtermfinal = cb1 * rho * nu_tilde * stild

    ! --- destruction term ---
    rr = (nu_tilde) / (kappa*kappa*ddw*ddw*stild + 1.0d-20)
    rr = min(rr, 10.0d0)

    gg = rr + cw2 * (rr**6 - rr)
    fw = gg * (((1.0d0 + cw3**6) / (gg**6 + cw3**6))**(1.0d0/6.0d0))

    destterm = cw1 * fw * rho * (nu_tilde/ddw)**2

    ! --- diffusion model ---
    ! 1/sigma * [ (mu+ρν~) laplacian(nu~)  +  cb2 ρ |grad(nu~)|² ]
    ! here squaret = |grad(nu~)|²  (your notation)

    difterm = (cb2 * rho * squaret) / sigma

    ! total source:
    source_t(1) = prodtermfinal + difterm - destterm




  case(2)		!k omega sst

			      eddyfl(1)=ielem_walldist(i)
			      eddyfl(2)=u_ct_val(1,1,i)
			      eddyfl(3)=u_ct_val(1,2,i)
			      eddyfl(4:6)=rec_grads(1,1:3,i)
			      eddyfl(7:9)=rec_grads(2,1:3,i)
			      eddyfl(10:12)=rec_grads(3,1:3,i)

			      eddyfl(13:15)=rec_grads(4,1:3,i)
			      eddyfl(16:18)=rec_grads(5,1:3,i)
								    

			      eddyfr=eddyfl

			      call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
      k_0=(max(verysmall,u_ct_val(1,1,i)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0=max(1.0e-1*ufreestream/charlength,u_ct_val(1,2,i)/leftv(1))
      wally=ielem_walldist(i)
	

	      !calculate here k and omega gradients
		  omx=rec_grads(6,1,i)
		  omy=rec_grads(6,2,i)
		  omz=rec_grads(6,3,i)
		  kx=rec_grads(5,1,i)
		  ky=rec_grads(5,2,i)
		  kz=rec_grads(5,3,i)

		  dervk_dervom= kx*omx+ky*omy+kz*omz

			!parameters

			!-----needed in diffusion--------
			d_omplus=max(2*leftv(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(1)/(leftv(1)*wally*wally*om_0))
			phi_1=min(phi_2, 4.0*leftv(1)*k_0/(sigma_om2*d_omplus*wally*wally)) 

				f_1=tanh(phi_1**4)
				f_2=tanh(phi_2**2)
				!--------------------------------


				alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
				alpha_star=alpha_starinf
				alpha_raw=alpha_inf


				beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
				alpha_star0=beta_i/3.0
				beta_stari=beta_starinf




				lowre=0
				!low-re correction------------------------------------------------------------------

				if (lowre.eq.1) then
				re_t_sst=leftv(1)*k_0/(viscl(1)*om_0)  !limiters for this???

				alpha_star=alpha_starinf*(alpha_star0+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)
				alpha_raw=alpha_inf/alpha_star*(alpha_0+re_t_sst/r_om_sst)/(1.0+re_t_sst/r_om_sst)

				beta_stari=beta_starinf*(4.0/15.0+(re_t_sst/r_beta)**4)/(1.0+(re_t_sst/r_beta)**4)

				end if
				      !--------------------------------------------------------------------------

				      !no mach number corrections for the beta
				      beta_star=beta_stari
				      beta_raw=beta_i


				      !production terms
					      !production of k
					      if (vort_model.eq.1) then
						      prod_k=viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm
						      prod_om=alpha_raw*leftv(1)*snorm*onorm
					      else
						      !prod_k=viscl(3)*snorm**2
						      !exact formulation: 
						      prod_k=viscl(3)*(snorm**2)!-2.0/3.0*divnorm*divnorm)-2.0/3.0*leftv(1)*k_0*divnorm
						      prod_k=min(prod_k,10.0*leftv(1)*beta_star*k_0*om_0)

						      !intelligent way of limiting: 
						      !prod_k=min(prod_k,max(10.0*leftv(1)*beta_star*k_0*om_0,&
						      !	    viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm))
					      
					      
						      !production of omega  (menter does this before correcting prod_k, 
						      !but in  article of 2003 he applies the correction to both)	
						      prod_om=alpha_raw*leftv(1)*snorm**2
				      ! 		prod_om=min(prod_om,10.0*leftv(1)*beta_star*om_0*om_0)
					      end if

				      !destruction terms
					      !destruction of k
					      ydest_k=leftv(1)*beta_star*k_0*om_0
					      !destruction of omega
					      ydest_om=leftv(1)*beta_raw*om_0*om_0
					      
				      !crossed-diffusion term
					      !crossed diffusion of omega
					      diff_om=2.0*(1-f_1)*leftv(1)/(om_0*sigma_om2)*dervk_dervom
	 
	
				!qsas term: scale adaptive
					if (qsas_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
						!calculate here second derivative of u 
						!declare all variables
						      do iex=1,3
							vortet(iex,1:3)=rec_grads(3+turbulenceequations+iex,1:3,i)
							
						      end do
						
						uxx=vortet(1,1) ; uyy=vortet(1,2); uzz=vortet(1,3)				
						vxx=vortet(2,1) ; vyy=vortet(2,2); vzz=vortet(2,3)	
						wxx=vortet(3,1) ; wyy=vortet(3,2); wzz=vortet(3,3)	
						!compute here cell volume
						cell_volume=ielem_totvolume(i)
						
						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
						dervk2=kx*kx+ky*ky+kz*kz
						dervom2=omx*omx+omy*omy+omz*omz
						
						delta_cell=cell_volume**0.333333333333333
						
						l_sas=sqrt(k_0)/(beta_star**0.25*om_0)
						l_vk=max(kappa*snorm/u_lapl, &       !this switch provides high wave-number damping
							c_smg*delta_cell*sqrt(kappa*eta2_sas/(beta_raw/beta_star-alpha_raw)))
						
						q_sas1=leftv(1)*eta2_sas*kappa*snorm**2*(l_sas/l_vk)**2
						q_sas2= -c_sas*2*leftv(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)

						q_sas=max(q_sas1+q_sas2,0.0)
					else
					q_sas=zero
					end if

				    !final source terms


					    !des-sst model (if qsas_model=2)
					    if (qsas_model .eq.2) then
					    cell_volume=ielem_totvolume(i)
					    delta_cell=cell_volume**0.333333333333333
					    l_t_des=sqrt(k_0)/(beta_star*om_0)
					    f_des_sst=max(1.0, l_t_des/(c_des_sst*delta_cell)*(1-f_2))
					    !the (1-f_2) is meant to protect the boundary layer. will result in same 
					    !separation point that standard s-a
					    ydest_k=f_des_sst*ydest_k
					    end if
					    
				    srcfull_k=prod_k-ydest_k
				    srcfull_om=prod_om-ydest_om+diff_om+q_sas


				    !filling the output vector
				      source_t(1) = srcfull_k
				      source_t(2) = srcfull_om

				  
	
end select

end subroutine sources


subroutine sources_derivatives(n,iconsidered,source_t)
implicit none
!> @brief
!> sources derivative computation for implicit time stepping
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,dimension(turbulenceequations),intent(inout)::source_t
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx
real::vhx,amp,dvel,omega,squaret,tch_x,tch_x3,tch_fv1,tch_fv2
real::tch_rs,tch_r,tch_g,tch_glim,tch_fw,tch_dif,tch_dest,tch_prod
integer::i,k,j,l,ihgt,ihgj,iex, lowre
real::snorm,onorm,divnorm,ax,ay,az,tch_shh,tch_sav,verysmall,onesix,prodterm1,stild,rr
real::gg,fw,destterm,fodt,srcfull,dbpr,dbdi,dbde,dby,dbx,prodtermfinal
real:: r_des,f_des, ddw,f_des_sst,l_t_des
real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,prodterm2,sfac,sss,usss,ssss,s_bar,kron
real:: uz,vz,wx,wy,wz
real:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !for sas only
real,dimension(3,3)::vortet,tvort,svort,ovort
real,dimension(3)::vortem,vertf,vertex
real,dimension(turbulenceequations,1:3)::derivturb
!declarations for k-omega
real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
real:: sigma_k, sigma_om, f_1, f_2,phi_1,phi_2, d_omplus !-------------- those are for diffusion too!
real:: alpha_raw,alpha_star, re_t_sst,alpha_inf
real:: beta_stari, beta_i, beta_raw, beta_star
real:: k_0, om_0, wally
real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2,rho,nu,nu_tilde,difterm,prod,dstild_dnut,dest,dfw_dgg
real:: kx,ky,kz,omx,omy,omz,dprod_dnut,ddif_dnut,ddest_dnut,dif, drr_dnut ,dfw_drr,dgg_drr
real,dimension(1:nof_variables)::leftv,rightv
real::mp_pinfl,gammal
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
    real,dimension(1)::etvm


i=iconsidered

verysmall = 10e-16
	
vortet(1:3,1:3) = rec_grads(1:3,1:3,i)


ux = vortet(1,1);uy = vortet(1,2);uz = vortet(1,3)
vx = vortet(2,1);vy = vortet(2,2);vz = vortet(2,3)
wx = vortet(3,1);wy = vortet(3,2);wz = vortet(3,3)


do ihgt=1,3
  do ihgj=1,3
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
end do

svort=0.5*(vortet+tvort)
ovort=0.5*(vortet-tvort)





snorm=sqrt(2.0d0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+(svort(1,3)*svort(1,3))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))+(svort(2,3)*svort(2,3))+& 
	       (svort(3,1)*svort(3,1))+(svort(3,2)*svort(3,2))+(svort(3,3)*svort(3,3))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
onorm=sqrt(2.0d0*((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+(ovort(1,3)*ovort(1,3))+&
	       (ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2))+(ovort(2,3)*ovort(2,3))+&
	       (ovort(3,1)*ovort(3,1))+(ovort(3,2)*ovort(3,2))+(ovort(3,3)*ovort(3,3))))
omega=onorm


divnorm=ux+vy+wz  !careful with the sign. if it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)+(wz*wz)))&
	+((uy+vx)*(uy+vx)+(uz+wx)*(uz+wx)+(wy+vz)*(wy+vz))&
	-(2.0/3.0*(ux+vy+wz)*(ux+vy+wz)))


derivturb(1,1:3) = rec_grads(5,1:3,i)
squaret=(sqrt((derivturb(1,1)**2)+(derivturb(1,2)**2)+(derivturb(1,3)**2)))**2


leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
rightv(1:nof_variables)=leftv(1:nof_variables)
call get_visc_conduct(n,leftv,rightv,viscl,laml)




turbmv(1)=u_ct_val(1,1,i)
turbmv(2)=turbmv(1)



 select case(turbulencemodel)
 
 
 
 case(1) !!spalart almaras model

    ! (compressible ρν~ formulation)

    rho = leftv(1)
    nu = viscl(1) / rho
    nu_tilde = turbmv(1) / rho       ! convert φ=ρν~ → ν~
    nu_tilde = max(nu_tilde, 1.0d-12 * visc)
    ! χ
    tch_x = nu_tilde / nu
    tch_x3 = tch_x*tch_x*tch_x

    tch_fv1 = tch_x3 / (tch_x3 + cv1*cv1*cv1)
    tch_fv2 = 1.0d0 - tch_x / (1.0d0 + tch_x*tch_fv1)

    ddw = ielem_walldist(i)
    ddw= max(ddw, 10e-12)

    ! ----- des corrections -----
    if (des_model .eq. 1) then
        cell_volume = ielem_totvolume(i)
        delta_cell = cell_volume**0.333333333333333d0
        ddw = min(ddw, c_des_sa*delta_cell)
    end if

    if (des_model .eq. 2) then
        cell_volume = ielem_totvolume(i)
        delta_cell = cell_volume**0.333333333333333d0

        r_des = min(10.0d0, (viscl(1)+viscl(3)) /&
        (snorm*(kappa*ddw)**2 + 1.0d-16) )
        f_des = 1.0d0 - tanh((8.0d0*r_des)**3)

        ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell,1.0d-16),&
        1.0d-16)
    end if

    ! ----- stilde -----
    prodterm1 = nu_tilde / (kappa*kappa*ddw*ddw)

    stild = max(omega + tch_fv2*prodterm1, 0.3d0*omega)

    !--------------------------
    !       production
    !--------------------------
    prod = cb1 * rho * nu_tilde * stild

    ! d(prod)/d(ρν~)
    dstild_dnut =&
    tch_fv2 / (kappa*kappa*ddw*ddw)   ! only active if stild=ω+fv2*prodterm1

    if (stild .eq. 0.3d0*omega) dstild_dnut = 0.0d0

    dprod_dnut =&
    cb1 * ( stild + nu_tilde*dstild_dnut )

    !--------------------------
    !       destruction
    !--------------------------
    rr = nu_tilde / (kappa*kappa*ddw*ddw*stild + 1.0d-20)
    rr = min(rr, 10.0d0)

    gg = rr + cw2*(rr**6 - rr)
    fw = gg*((1.0d0+cw3**6)/(gg**6+cw3**6))**(1.0d0/6.0d0)

    dest = cw1 * fw * rho * (nu_tilde/ddw)**2

    ! d(fw)/d(rr)
    dgg_drr = 1.0d0 + cw2*(6.0d0*rr**5 - 1.0d0)

    dfw_dgg =&
    ((1.0d0+cw3**6)/(gg**6+cw3**6))**(1.0d0/6.0d0)&
    * (1.0d0 - (gg**6/(gg**6+cw3**6)))

    dfw_drr = dfw_dgg * dgg_drr




    ! drr/d(ρν~)
    drr_dnut =&
    (kappa*kappa*ddw*ddw*stild - nu_tilde*kappa*kappa*ddw*ddw*dstild_dnut)&
    /(kappa*kappa*ddw*ddw*stild + 1.0d-20)**2

    ! destruction derivative
    ddest_dnut =&
    cw1*( dfw_drr*drr_dnut * rho*(nu_tilde/ddw)**2&
    + fw*rho*2.0d0*(nu_tilde/ddw**2) / ddw )

    !--------------------------
    !       diffusion
    !--------------------------
    ! your model:  dif = cb2*rho*|∇ν~|² / σ
    dif = cb2 * rho * squaret / sigma


    ! jacobian:
    ! d(|grad(ν~)|²) / d(ρν~) = 0  (no dependence on ν~ itself)
    ddif_dnut = 0.0d0

    !--------------------------
    !     total jacobian
    !--------------------------
    source_t(1) = dprod_dnut + ddif_dnut - ddest_dnut


  case(2)		!k omega sst

			      eddyfl(1)=ielem_walldist(i)
			      eddyfl(2)=u_ct_val(1,1,i)
			      eddyfl(3)=u_ct_val(1,2,i)
			      eddyfl(4:6)=rec_grads(1,1:3,i)
			      eddyfl(7:9)=rec_grads(2,1:3,i)
			      eddyfl(10:12)=rec_grads(3,1:3,i)

			      eddyfl(13:15)=rec_grads(4,1:3,i)
			      eddyfl(16:18)=rec_grads(5,1:3,i)
								    

			      eddyfr=eddyfl

			      call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
      k_0=(max(verysmall,u_ct_val(1,1,i)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0=max(1.0e-1*ufreestream/charlength,u_ct_val(1,2,i)/leftv(1))
!       wally=ielem_walldist(i)
! 	
! 
! 	      !calculate here k and omega gradients
! 		  omx=rec_grads(5,1,i)
! 		  omy=rec_grads(5,2,i)
! 		  omz=rec_grads(5,3,i)
! 		  kx=rec_grads(4,1,i)
! 		  ky=rec_grads(4,2,i)
! 		  kz=rec_grads(4,3,i)
! 
! 		  dervk_dervom= kx*omx+ky*omy+kz*omz
! 
! 			!parameters
! 
! 			!-----needed in diffusion--------
! 			d_omplus=max(2*leftv(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
! 			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(1)/(leftv(1)*wally*wally*om_0))
! 			phi_1=min(phi_2, 4.0*leftv(1)*k_0/(sigma_om2*d_omplus*wally*wally)) 
! 
! 				f_1=tanh(phi_1**4)
! 				f_2=tanh(phi_2**2)
! 				!--------------------------------
! 
! 
! 				alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
! 				alpha_star=alpha_starinf
! 				alpha_raw=alpha_inf
! 
! 
! 				beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
! 				alpha_star0=beta_i/3.0
! 				beta_stari=beta_starinf
! 
! 
! 
! 
! 				lowre=0
! 				!low-re correction------------------------------------------------------------------
! 
! 				if (lowre.eq.1) then
! 				re_t_sst=leftv(1)*k_0/(viscl(1)*om_0)  !limiters for this???
! 
! 				alpha_star=alpha_starinf*(alpha_star0+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)
! 				alpha_raw=alpha_inf/alpha_star*(alpha_0+re_t_sst/r_om_sst)/(1.0+re_t_sst/r_om_sst)
! 
! 				beta_stari=beta_starinf*(4.0/15.0+(re_t_sst/r_beta)**4)/(1.0+(re_t_sst/r_beta)**4)
! 
! 				end if
! 				      !--------------------------------------------------------------------------
! 
! 				      !no mach number corrections for the beta
! 				      beta_star=beta_stari
! 				      beta_raw=beta_i
! 
! 
! 				      !production terms
! 					      !production of k
! 					      if (vort_model.eq.1) then
! 						      prod_k=viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm
! 						      prod_om=alpha_raw*leftv(1)*snorm*onorm
! 					      else
! 						      !prod_k=viscl(3)*snorm**2
! 						      !exact formulation: 
! 						      prod_k=viscl(3)*(snorm**2)!-2.0/3.0*divnorm*divnorm)-2.0/3.0*leftv(1)*k_0*divnorm
! 						      prod_k=min(prod_k,10.0*leftv(1)*beta_star*k_0*om_0)
! 
! 						      !intelligent way of limiting: 
! 						      !prod_k=min(prod_k,max(10.0*leftv(1)*beta_star*k_0*om_0,&
! 						      !	    viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm))
! 					      
! 					      
! 						      !production of omega  (menter does this before correcting prod_k, 
! 						      !but in  article of 2003 he applies the correction to both)	
! 						      prod_om=alpha_raw*leftv(1)*snorm**2
! 				      ! 		prod_om=min(prod_om,10.0*leftv(1)*beta_star*om_0*om_0)
! 					      end if
! 
! 				      !destruction terms
! 					      !destruction of k
! 					      ydest_k=leftv(1)*beta_star*k_0*om_0
! 					      !destruction of omega
! 					      ydest_om=leftv(1)*beta_raw*om_0*om_0
! 					      
! 				      !crossed-diffusion term
! 					      !crossed diffusion of omega
! 					      diff_om=2.0*(1-f_1)*leftv(1)/(om_0*sigma_om2)*dervk_dervom
! 	 
! 	
! 				!qsas term: scale adaptive
! 					if (qsas_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
! 						!calculate here second derivative of u 
! 						!declare all variables
! 						      do iex=1,3
! 							vortet(iex,1:3)=rec_grads(3+turbulenceequations+iex,1:3,i)
! 							
! 						      end do
! 						
! 						uxx=vortet(1,1) ; uyy=vortet(1,2); uzz=vortet(1,3)				
! 						vxx=vortet(2,1) ; vyy=vortet(2,2); vzz=vortet(2,3)	
! 						wxx=vortet(3,1) ; wyy=vortet(3,2); wzz=vortet(3,3)	
! 						!compute here cell volume
! 						cell_volume=ielem_totvolume(i)
! 						
! 						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
! 						dervk2=kx*kx+ky*ky+kz*kz
! 						dervom2=omx*omx+omy*omy+omz*omz
! 						
! 						delta_cell=cell_volume**0.333333333333333
! 						
! 						l_sas=sqrt(k_0)/(beta_star**0.25*om_0)
! 						l_vk=max(kappa*snorm/u_lapl, &       !this switch provides high wave-number damping
! 							c_smg*delta_cell*sqrt(kappa*eta2_sas/(beta_raw/beta_star-alpha_raw)))
! 						
! 						q_sas1=leftv(1)*eta2_sas*kappa*snorm**2*(l_sas/l_vk)**2
! 						q_sas2= -c_sas*2*leftv(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)
! 
! 						q_sas=max(q_sas1+q_sas2,0.0)
! 					else
! 					q_sas=zero
! 					end if
! 
! 				    !final source terms
! 
! 
! 					    !des-sst model (if qsas_model=2)
! 					    if (qsas_model .eq.2) then
! 					    cell_volume=ielem_totvolume(i)
! 					    delta_cell=cell_volume**0.333333333333333
! 					    l_t_des=sqrt(k_0)/(beta_star*om_0)
! 					    f_des_sst=max(1.0, l_t_des/(c_des_sst*delta_cell)*(1-f_2))
! 					    !the (1-f_2) is meant to protect the boundary layer. will result in same 
! 					    !separation point that standard s-a
! 					    ydest_k=f_des_sst*ydest_k
! 					    end if
! 					    
! 				    srcfull_k=prod_k-ydest_k
! 				    srcfull_om=prod_om-ydest_om+diff_om+q_sas


				    !filling the output vector
				      source_t(1) = beta_starinf*om_0
				      source_t(2) = beta_starinf*k_0

				  
	
end select

end subroutine sources_derivatives

subroutine sources_computation2d(n)
	implicit none
!> @brief
!> sources  computation in 2d
	integer,intent(in)::n
	integer::i,kmaxe,iconsidered
	real,dimension(turbulenceequations)::source_t
	real,dimension(1:nof_variables)::source_r,src_jac_diag

	
	
	kmaxe=xmpielrank(n)
	
	
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, iconsidered, source_t, source_r, src_jac_diag)
#else
!$omp do
#endif
	do i=1,kmaxe
		iconsidered=i

		if (turbulence.eq.1)then

		call sources2d(n,iconsidered,source_t)
		rhst_val(1:turbulenceequations,i)=rhst_val(1:turbulenceequations,i)-(source_t(1:turbulenceequations)*ielem_totvolume(i))
		end if
		if (realgas.eq.1)then
		if (rg_relax.eq.2)then
		call sources_realgas(n,iconsidered,source_r,src_jac_diag)
		rhs_val(1:nof_variables,i)=rhs_val(1:nof_variables,i)-(source_r(1:nof_variables)*ielem_totvolume(i))
		sht_rg(i,1:nof_variables)=src_jac_diag(1:nof_variables)*ielem_totvolume(i)
		end if
		end if



	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	
	
end subroutine sources_computation2d

subroutine sources_derivatives_computation2d(n)
	implicit none
!> @brief
!> sources  derivatives computation in 2d
	integer,intent(in)::n
	integer::i,kmaxe,iconsidered
	real,dimension(turbulenceequations)::source_t
	
	
	kmaxe=xmpielrank(n)
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, iconsidered, source_t)
#else
!$omp do
#endif
	do i=1,kmaxe
		iconsidered=i
		call sources_derivatives2d(n,iconsidered,source_t)
		sht(i,1:turbulenceequations)=(source_t(1:turbulenceequations)*ielem_totvolume(i))
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine sources_derivatives_computation2d



subroutine sources2d(n,iconsidered,source_t)
implicit none
!> @brief
!> sources  computation in 2d
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,dimension(turbulenceequations),intent(inout)::source_t
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx
real::vhx,amp,dvel,omega,squaret,tch_x,tch_x3,tch_fv1,tch_fv2
real::tch_rs,tch_r,tch_g,tch_glim,tch_fw,tch_dif,tch_dest,tch_prod
integer::i,k,j,l,ihgt,ihgj,iex, lowre
real::snorm,onorm,divnorm,ax,ay,az,tch_shh,tch_sav,verysmall,onesix,prodterm1,stild,rr
real::gg,fw,destterm,fodt,srcfull,dbpr,dbdi,dbde,dby,dbx,prodtermfinal
real:: r_des,f_des, ddw,f_des_sst,l_t_des
real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,prodterm2,sfac,sss,usss,ssss,s_bar,kron
real:: uz,vz,wx,wy,wz
real:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !for sas only
real,dimension(2,2)::vortet,tvort,svort,ovort
real,dimension(2)::vortem,vertf,vertex
real,dimension(turbulenceequations,1:2)::derivturb
!declarations for k-omega
real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
real:: sigma_k, sigma_om, f_1, f_2,phi_1,phi_2, d_omplus !-------------- those are for diffusion too!
real:: alpha_raw,alpha_star, re_t_sst,alpha_inf
real:: beta_stari, beta_i, beta_raw, beta_star
real:: k_0, om_0, wally
real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
real:: kx,ky,kz,omx,omy,omz,rho,nu,nu_tilde,difterm,prod
real,dimension(1:nof_variables)::leftv,rightv
real::mp_pinfl,gammal
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
real,dimension(1)::etvm

i=iconsidered

verysmall = 10e-16



	
vortet(1:2,1:2) = rec_grads(1:2,1:2,i)


ux = vortet(1,1);uy = vortet(1,2)
vx = vortet(2,1);vy = vortet(2,2)







do ihgt=1,2
  do ihgj=1,2
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
end do

svort=0.5*(vortet+tvort)
ovort=0.5*(vortet-tvort)





snorm=sqrt(2.0d0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
onorm=sqrt(2.0d0*((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+&
	       (ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2))))
omega=onorm


divnorm=ux+vy+wz  !careful with the sign. if it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)))&
	+((uy+vx)*(uy+vx))&
	-(2.0/3.0*(ux+vy)*(ux+vy)))


derivturb(1,1:2) = rec_grads(4,1:2,i)
squaret=(sqrt((derivturb(1,1)**2)+(derivturb(1,2)**2)))**2


leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
rightv(1:nof_variables)=leftv(1:nof_variables)

call get_visc_conduct(n,leftv,rightv,viscl,laml)




turbmv(1)=u_ct_val(1,1,i)
turbmv(2)=turbmv(1)



 select case(turbulencemodel)
 
 
 


    
      case(1)   !! spalart–almaras model (compressible ρν~ formulation)

    ! optional rotation/curvature correction
    if (rot_corr .eq. 1) then
        omega = omega + 2.0d0 * min(0.0d0, snorm - onorm)
    end if

    eddyfl(2) = turbmv(1)
    eddyfr(2) = turbmv(1)

    call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)

    rho = leftv(1)
    nu  = viscl(1) / rho
    nu_tilde = turbmv(1) / rho      ! convert φ = ρ ν~  →  ν~
    nu_tilde= max(nu_tilde, 1.0d-12 * visc)

    ! χ = ν~/ν
    tch_x = nu_tilde / nu
    tch_x3 = tch_x*tch_x*tch_x

    ! fv1
    tch_fv1 = tch_x3 / (tch_x3 + (cv1*cv1*cv1))

    ! fv2
    tch_fv2 = 1.0d0 - tch_x / (1.0d0 + tch_x * tch_fv1)

    ! physical wall distance
    ddw = ielem_walldist(i)
    ddw= max(ddw, 10e-12)

    ! --- des correction ---
    if (des_model .eq. 1) then
        cell_volume = ielem_totvolume(i)
        delta_cell  = cell_volume**0.333333333333333d0
        ddw = min(ddw, c_des_sa * delta_cell)
    end if

    if (des_model .eq. 2) then
        cell_volume = ielem_totvolume(i)
        delta_cell  = cell_volume**0.333333333333333d0

        r_des = min(10.0d0, (viscl(1)+viscl(3)) / (snorm*(kappa*ddw)**2 + 1.0d-16))
        f_des = 1.0d0 - tanh((8.0d0*r_des)**3)

        ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell, 1.0d-16), 1.0d-16)
    end if

    ! --- stilde ---
    prodterm1 = (nu_tilde) / (kappa*kappa*ddw*ddw)
    stild = max( omega + tch_fv2*prodterm1 , 0.3d0 * omega )

    ! --- production term ---
    prodtermfinal = cb1 * rho * nu_tilde * stild

    ! --- destruction term ---
    rr = (nu_tilde) / (kappa*kappa*ddw*ddw*stild + 1.0d-20)
    rr = min(rr, 10.0d0)

    gg = rr + cw2 * (rr**6 - rr)
    fw = gg * (((1.0d0 + cw3**6) / (gg**6 + cw3**6))**(1.0d0/6.0d0))

    destterm = cw1 * fw * rho * (nu_tilde/ddw)**2

    ! --- diffusion model ---
    ! 1/sigma * [ (mu+ρν~) laplacian(nu~)  +  cb2 ρ |grad(nu~)|² ]
    ! here squaret = |grad(nu~)|²  (your notation)

    difterm = (cb2 * rho * squaret) / sigma

    ! total source:
    source_t(1) = prodtermfinal + difterm - destterm


  case(2)		!k omega sst

			      eddyfl(1)=ielem_walldist(i)
			      eddyfl(2)=u_ct_val(1,1,i)
			      eddyfl(3)=u_ct_val(1,2,i)
			      eddyfl(4:5)=rec_grads(1,1:2,i)
			      eddyfl(6:7)=rec_grads(2,1:2,i)
			      eddyfl(8:9)=rec_grads(4,1:2,i)
			      eddyfl(10:11)=rec_grads(5,1:2,i)
								    

			      eddyfr=eddyfl

			      call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
      k_0=(max(verysmall,u_ct_val(1,1,i)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0=max(1.0e-1*ufreestream/charlength,u_ct_val(1,2,i)/leftv(1))
      wally=ielem_walldist(i)
	

	      !calculate here k and omega gradients
		  omx=rec_grads(5,1,i)
		  omy=rec_grads(5,2,i)
		  
		  kx=rec_grads(4,1,i)
		  ky=rec_grads(4,2,i)
		

		  dervk_dervom= kx*omx+ky*omy

			!parameters

			!-----needed in diffusion--------
			d_omplus=max(2*leftv(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(1)/(leftv(1)*wally*wally*om_0))
			phi_1=min(phi_2, 4.0*leftv(1)*k_0/(sigma_om2*d_omplus*wally*wally)) 

				f_1=tanh(phi_1**4)
				f_2=tanh(phi_2**2)
				!--------------------------------


				alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
				alpha_star=alpha_starinf
				alpha_raw=alpha_inf


				beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
				alpha_star0=beta_i/3.0
				beta_stari=beta_starinf




				lowre=0
				!low-re correction------------------------------------------------------------------

				if (lowre.eq.1) then
				re_t_sst=leftv(1)*k_0/(viscl(1)*om_0)  !limiters for this???

				alpha_star=alpha_starinf*(alpha_star0+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)
				alpha_raw=alpha_inf/alpha_star*(alpha_0+re_t_sst/r_om_sst)/(1.0+re_t_sst/r_om_sst)

				beta_stari=beta_starinf*(4.0/15.0+(re_t_sst/r_beta)**4)/(1.0+(re_t_sst/r_beta)**4)

				end if
				      !--------------------------------------------------------------------------

				      !no mach number corrections for the beta
				      beta_star=beta_stari
				      beta_raw=beta_i


				      !production terms
					      !production of k
					      if (vort_model.eq.1) then
						      prod_k=viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm
						      prod_om=alpha_raw*leftv(1)*snorm*onorm
					      else
						      !prod_k=viscl(3)*snorm**2
						      !exact formulation: 
						      prod_k=viscl(3)*(snorm**2)!-2.0/3.0*divnorm*divnorm)-2.0/3.0*leftv(1)*k_0*divnorm
						      prod_k=min(prod_k,10.0*leftv(1)*beta_star*k_0*om_0)

						      !intelligent way of limiting: 
						      !prod_k=min(prod_k,max(10.0*leftv(1)*beta_star*k_0*om_0,&
						      !	    viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm))
					      
					      
						      !production of omega  (menter does this before correcting prod_k, 
						      !but in  article of 2003 he applies the correction to both)	
						      prod_om=alpha_raw*leftv(1)*snorm**2
				      ! 		prod_om=min(prod_om,10.0*leftv(1)*beta_star*om_0*om_0)
					      end if

				      !destruction terms
					      !destruction of k
					      ydest_k=leftv(1)*beta_star*k_0*om_0
					      !destruction of omega
					      ydest_om=leftv(1)*beta_raw*om_0*om_0
					      
				      !crossed-diffusion term
					      !crossed diffusion of omega
					      diff_om=2.0*(1-f_1)*leftv(1)/(om_0*sigma_om2)*dervk_dervom
	 
	
				!qsas term: scale adaptive
					if (qsas_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
						!calculate here second derivative of u 
						!declare all variables
						      do iex=1,2
							vortet(iex,1:2)=rec_grads(3+turbulenceequations+iex,1:2,i)
							
						      end do
						
						uxx=vortet(1,1) ; uyy=vortet(1,2); 				
						vxx=vortet(2,1) ; vyy=vortet(2,2); 	
						!compute here cell volume
						cell_volume=ielem_totvolume(i)
						
						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
						dervk2=kx*kx+ky*ky+kz*kz
						dervom2=omx*omx+omy*omy+omz*omz
						
						delta_cell=cell_volume**0.333333333333333
						
						l_sas=sqrt(k_0)/(beta_star**0.25*om_0)
						l_vk=max(kappa*snorm/u_lapl, &       !this switch provides high wave-number damping
							c_smg*delta_cell*sqrt(kappa*eta2_sas/(beta_raw/beta_star-alpha_raw)))
						
						q_sas1=leftv(1)*eta2_sas*kappa*snorm**2*(l_sas/l_vk)**2
						q_sas2= -c_sas*2*leftv(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)

						q_sas=max(q_sas1+q_sas2,0.0)
					else
					q_sas=zero
					end if

				    !final source terms


					    !des-sst model (if qsas_model=2)
					    if (qsas_model .eq.2) then
					    cell_volume=ielem_totvolume(i)
					    delta_cell=cell_volume**0.333333333333333
					    l_t_des=sqrt(k_0)/(beta_star*om_0)
					    f_des_sst=max(1.0, l_t_des/(c_des_sst*delta_cell)*(1-f_2))
					    !the (1-f_2) is meant to protect the boundary layer. will result in same 
					    !separation point that standard s-a
					    ydest_k=f_des_sst*ydest_k
					    end if
					    
				    srcfull_k=prod_k-ydest_k
				    srcfull_om=prod_om-ydest_om+diff_om+q_sas


				    !filling the output vector
				      source_t(1) = srcfull_k
				      source_t(2) = srcfull_om

				  
	
end select





end subroutine sources2d


subroutine sources_derivatives2d(n,iconsidered,source_t)
implicit none
!> @brief
!> sources derivatives computation in 2d
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,dimension(turbulenceequations),intent(inout)::source_t
real::intenergy,r1,u1,v1,w1,et1,s1,ie1,p1,skin1,e1,rs,us,vs,ws,khx
real::vhx,amp,dvel,omega,squaret,tch_x,tch_x3,tch_fv1,tch_fv2
real::tch_rs,tch_r,tch_g,tch_glim,tch_fw,tch_dif,tch_dest,tch_prod
integer::i,k,j,l,ihgt,ihgj,iex, lowre
real::snorm,onorm,divnorm,ax,ay,az,tch_shh,tch_sav,verysmall,onesix,prodterm1,stild,rr
real::gg,fw,destterm,fodt,srcfull,dbpr,dbdi,dbde,dby,dbx,prodtermfinal
real:: r_des,f_des, ddw,f_des_sst,l_t_des
real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,prodterm2,sfac,sss,usss,ssss,s_bar,kron
real:: uz,vz,wx,wy,wz
real:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !for sas only
real,dimension(2,2)::vortet,tvort,svort,ovort
real,dimension(2)::vortem,vertf,vertex
real,dimension(turbulenceequations,1:3)::derivturb
!declarations for k-omega
real:: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
real:: sigma_k, sigma_om, f_1, f_2,phi_1,phi_2, d_omplus !-------------- those are for diffusion too!
real:: alpha_raw,alpha_star, re_t_sst,alpha_inf
real:: beta_stari, beta_i, beta_raw, beta_star
real:: k_0, om_0, wally
real:: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
real:: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2,rho,nu,nu_tilde,difterm,prod,dstild_dnut,dest,dfw_dgg
real:: kx,ky,kz,omx,omy,omz,dprod_dnut,ddif_dnut,ddest_dnut,dif, drr_dnut ,dfw_drr,dgg_drr
real,dimension(1:nof_variables)::leftv,rightv
real::mp_pinfl,gammal
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
    real,dimension(1)::etvm

i=iconsidered

verysmall = 10e-16
	
vortet(1:2,1:2) = rec_grads(1:2,1:2,i)


ux = vortet(1,1);uy = vortet(1,2)
vx = vortet(2,1);vy = vortet(2,2)


do ihgt=1,2
  do ihgj=1,2
  tvort(ihgt,ihgj)=vortet(ihgj,ihgt)
  end do
end do

svort=0.5*(vortet+tvort)
ovort=0.5*(vortet-tvort)





snorm=sqrt(2.0d0*((svort(1,1)*svort(1,1))+(svort(1,2)*svort(1,2))+&
	       (svort(2,1)*svort(2,1))+(svort(2,2)*svort(2,2))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
onorm=sqrt(2.0d0*((ovort(1,1)*ovort(1,1))+(ovort(1,2)*ovort(1,2))+&
	       (ovort(2,1)*ovort(2,1))+(ovort(2,2)*ovort(2,2))))
omega=onorm


divnorm=ux+vy !careful with the sign. if it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)))&
	+((uy+vx)*(uy+vx))&
	-(2.0/3.0*(ux+vy)*(ux+vy)))


derivturb(1,1:2) = rec_grads(4,1:2,i)
squaret=(sqrt((derivturb(1,1)**2)+(derivturb(1,2)**2)))**2


leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
rightv(1:nof_variables)=leftv(1:nof_variables)
call get_visc_conduct(n,leftv,rightv,viscl,laml)




turbmv(1)=u_ct_val(1,1,i)
turbmv(2)=turbmv(1)


 select case(turbulencemodel)
 
 
 
 case(1) !!spalart almaras model

    ! (compressible ρν~ formulation)

    rho = leftv(1)
    nu = viscl(1) / rho
    nu_tilde = turbmv(1) / rho       ! convert φ=ρν~ → ν~
    nu_tilde = max(nu_tilde, 1.0d-12 * visc)

    ! χ
    tch_x = nu_tilde / nu
    tch_x3 = tch_x*tch_x*tch_x

    tch_fv1 = tch_x3 / (tch_x3 + cv1*cv1*cv1)
    tch_fv2 = 1.0d0 - tch_x / (1.0d0 + tch_x*tch_fv1)

    ddw = ielem_walldist(i)
    ddw= max(ddw, 10e-12)

    ! ----- des corrections -----
    if (des_model .eq. 1) then
        cell_volume = ielem_totvolume(i)
        delta_cell = cell_volume**0.333333333333333d0
        ddw = min(ddw, c_des_sa*delta_cell)
    end if

    if (des_model .eq. 2) then
        cell_volume = ielem_totvolume(i)
        delta_cell = cell_volume**0.333333333333333d0

        r_des = min(10.0d0, (viscl(1)+viscl(3)) /&
        (snorm*(kappa*ddw)**2 + 1.0d-16) )
        f_des = 1.0d0 - tanh((8.0d0*r_des)**3)

        ddw = max(ddw - f_des*max(ddw - c_des_sa*delta_cell,1.0d-16),&
        1.0d-16)
    end if

    ! ----- stilde -----
    prodterm1 = nu_tilde / (kappa*kappa*ddw*ddw)

    stild = max(omega + tch_fv2*prodterm1, 0.3d0*omega)

    !--------------------------
    !       production
    !--------------------------
    prod = cb1 * rho * nu_tilde * stild

    ! d(prod)/d(ρν~)
    dstild_dnut =&
    tch_fv2 / (kappa*kappa*ddw*ddw)   ! only active if stild=ω+fv2*prodterm1

    if (stild .eq. 0.3d0*omega) dstild_dnut = 0.0d0

    dprod_dnut =&
    cb1 * ( stild + nu_tilde*dstild_dnut )

    !--------------------------
    !       destruction
    !--------------------------
    rr = nu_tilde / (kappa*kappa*ddw*ddw*stild + 1.0d-20)
    rr = min(rr, 10.0d0)

    gg = rr + cw2*(rr**6 - rr)
    fw = gg*((1.0d0+cw3**6)/(gg**6+cw3**6))**(1.0d0/6.0d0)

    dest = cw1 * fw * rho * (nu_tilde/ddw)**2

    ! d(fw)/d(rr)
    dgg_drr = 1.0d0 + cw2*(6.0d0*rr**5 - 1.0d0)

    dfw_dgg =&
    ((1.0d0+cw3**6)/(gg**6+cw3**6))**(1.0d0/6.0d0)&
    * (1.0d0 - (gg**6/(gg**6+cw3**6)))

    dfw_drr = dfw_dgg * dgg_drr




    ! drr/d(ρν~)
    drr_dnut =&
    (kappa*kappa*ddw*ddw*stild - nu_tilde*kappa*kappa*ddw*ddw*dstild_dnut)&
    /(kappa*kappa*ddw*ddw*stild + 1.0d-20)**2

    ! destruction derivative
    ddest_dnut =&
    cw1*( dfw_drr*drr_dnut * rho*(nu_tilde/ddw)**2&
    + fw*rho*2.0d0*(nu_tilde/ddw**2) / ddw )

    !--------------------------
    !       diffusion
    !--------------------------
    ! your model:  dif = cb2*rho*|∇ν~|² / σ
    dif = cb2 * rho * squaret / sigma


    ! jacobian:
    ! d(|grad(ν~)|²) / d(ρν~) = 0  (no dependence on ν~ itself)
    ddif_dnut = 0.0d0

    !--------------------------
    !     total jacobian
    !--------------------------
    source_t(1) = dprod_dnut + ddif_dnut - ddest_dnut


  case(2)		!k omega sst

			   eddyfl(1)=ielem_walldist(i)
			      eddyfl(2)=u_ct_val(1,1,i)
			      eddyfl(3)=u_ct_val(1,2,i)
			      eddyfl(4:5)=rec_grads(1,1:2,i)
			      eddyfl(6:7)=rec_grads(2,1:2,i)
			      eddyfl(8:9)=rec_grads(4,1:2,i)
			      eddyfl(10:11)=rec_grads(5,1:2,i)
								    

			      eddyfr=eddyfl

			      call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
      k_0=(max(verysmall,u_ct_val(1,1,i)/leftv(1))) !first subindex makes reference to the time-stepping
      om_0=max(1.0e-1*ufreestream/charlength,u_ct_val(1,2,i)/leftv(1))
!       wally=ielem_walldist(i)
! 	
! 
! 	      !calculate here k and omega gradients
! 		  omx=rec_grads(5,1,i)
! 		  omy=rec_grads(5,2,i)
! 		  omz=rec_grads(5,3,i)
! 		  kx=rec_grads(4,1,i)
! 		  ky=rec_grads(4,2,i)
! 		  kz=rec_grads(4,3,i)
! 
! 		  dervk_dervom= kx*omx+ky*omy+kz*omz
! 
! 			!parameters
! 
! 			!-----needed in diffusion--------
! 			d_omplus=max(2*leftv(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
! 			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*viscl(1)/(leftv(1)*wally*wally*om_0))
! 			phi_1=min(phi_2, 4.0*leftv(1)*k_0/(sigma_om2*d_omplus*wally*wally)) 
! 
! 				f_1=tanh(phi_1**4)
! 				f_2=tanh(phi_2**2)
! 				!--------------------------------
! 
! 
! 				alpha_inf=f_1*alpha_inf1+(1.0-f_1)*alpha_inf2
! 				alpha_star=alpha_starinf
! 				alpha_raw=alpha_inf
! 
! 
! 				beta_i=f_1*beta_i1+(1.0-f_1)*beta_i2
! 				alpha_star0=beta_i/3.0
! 				beta_stari=beta_starinf
! 
! 
! 
! 
! 				lowre=0
! 				!low-re correction------------------------------------------------------------------
! 
! 				if (lowre.eq.1) then
! 				re_t_sst=leftv(1)*k_0/(viscl(1)*om_0)  !limiters for this???
! 
! 				alpha_star=alpha_starinf*(alpha_star0+re_t_sst/r_k_sst)/(1.0+re_t_sst/r_k_sst)
! 				alpha_raw=alpha_inf/alpha_star*(alpha_0+re_t_sst/r_om_sst)/(1.0+re_t_sst/r_om_sst)
! 
! 				beta_stari=beta_starinf*(4.0/15.0+(re_t_sst/r_beta)**4)/(1.0+(re_t_sst/r_beta)**4)
! 
! 				end if
! 				      !--------------------------------------------------------------------------
! 
! 				      !no mach number corrections for the beta
! 				      beta_star=beta_stari
! 				      beta_raw=beta_i
! 
! 
! 				      !production terms
! 					      !production of k
! 					      if (vort_model.eq.1) then
! 						      prod_k=viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm
! 						      prod_om=alpha_raw*leftv(1)*snorm*onorm
! 					      else
! 						      !prod_k=viscl(3)*snorm**2
! 						      !exact formulation: 
! 						      prod_k=viscl(3)*(snorm**2)!-2.0/3.0*divnorm*divnorm)-2.0/3.0*leftv(1)*k_0*divnorm
! 						      prod_k=min(prod_k,10.0*leftv(1)*beta_star*k_0*om_0)
! 
! 						      !intelligent way of limiting: 
! 						      !prod_k=min(prod_k,max(10.0*leftv(1)*beta_star*k_0*om_0,&
! 						      !	    viscl(3)*onorm*snorm!-2.0/3.0*leftv(1)*k_0*divnorm))
! 					      
! 					      
! 						      !production of omega  (menter does this before correcting prod_k, 
! 						      !but in  article of 2003 he applies the correction to both)	
! 						      prod_om=alpha_raw*leftv(1)*snorm**2
! 				      ! 		prod_om=min(prod_om,10.0*leftv(1)*beta_star*om_0*om_0)
! 					      end if
! 
! 				      !destruction terms
! 					      !destruction of k
! 					      ydest_k=leftv(1)*beta_star*k_0*om_0
! 					      !destruction of omega
! 					      ydest_om=leftv(1)*beta_raw*om_0*om_0
! 					      
! 				      !crossed-diffusion term
! 					      !crossed diffusion of omega
! 					      diff_om=2.0*(1-f_1)*leftv(1)/(om_0*sigma_om2)*dervk_dervom
! 	 
! 	
! 				!qsas term: scale adaptive
! 					if (qsas_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
! 						!calculate here second derivative of u 
! 						!declare all variables
! 						      do iex=1,3
! 							vortet(iex,1:3)=rec_grads(3+turbulenceequations+iex,1:3,i)
! 							
! 						      end do
! 						
! 						uxx=vortet(1,1) ; uyy=vortet(1,2); uzz=vortet(1,3)				
! 						vxx=vortet(2,1) ; vyy=vortet(2,2); vzz=vortet(2,3)	
! 						wxx=vortet(3,1) ; wyy=vortet(3,2); wzz=vortet(3,3)	
! 						!compute here cell volume
! 						cell_volume=ielem_totvolume(i)
! 						
! 						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
! 						dervk2=kx*kx+ky*ky+kz*kz
! 						dervom2=omx*omx+omy*omy+omz*omz
! 						
! 						delta_cell=cell_volume**0.333333333333333
! 						
! 						l_sas=sqrt(k_0)/(beta_star**0.25*om_0)
! 						l_vk=max(kappa*snorm/u_lapl, &       !this switch provides high wave-number damping
! 							c_smg*delta_cell*sqrt(kappa*eta2_sas/(beta_raw/beta_star-alpha_raw)))
! 						
! 						q_sas1=leftv(1)*eta2_sas*kappa*snorm**2*(l_sas/l_vk)**2
! 						q_sas2= -c_sas*2*leftv(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)
! 
! 						q_sas=max(q_sas1+q_sas2,0.0)
! 					else
! 					q_sas=zero
! 					end if
! 
! 				    !final source terms
! 
! 
! 					    !des-sst model (if qsas_model=2)
! 					    if (qsas_model .eq.2) then
! 					    cell_volume=ielem_totvolume(i)
! 					    delta_cell=cell_volume**0.333333333333333
! 					    l_t_des=sqrt(k_0)/(beta_star*om_0)
! 					    f_des_sst=max(1.0, l_t_des/(c_des_sst*delta_cell)*(1-f_2))
! 					    !the (1-f_2) is meant to protect the boundary layer. will result in same 
! 					    !separation point that standard s-a
! 					    ydest_k=f_des_sst*ydest_k
! 					    end if
! 					    
! 				    srcfull_k=prod_k-ydest_k
! 				    srcfull_om=prod_om-ydest_om+diff_om+q_sas


				    !filling the output vector
				      source_t(1) = beta_starinf*om_0
				      source_t(2) = beta_starinf*k_0

				  
	
end select

end subroutine sources_derivatives2d

end module source
