SUBROUTINE MULTISPECIES_MIXTURES_RG(LEFTV, MP_mu_mix, MP_ktr_mix, MP_kve, &
                                    MP_D_eff, MP_HTR, MP_HVIB, GAMMAL)

  implicit none
  real, dimension(1:nof_Variables), intent(in)    :: LEFTV
  real, intent(inout)                             :: MP_mu_mix, MP_ktr_mix, MP_kve, GAMMAL
  real, dimension(1:nof_species), intent(inout)   :: MP_HTR, MP_HVIB, MP_D_eff

  ! Local
  real, dimension(1:nof_variables) :: prim
  real :: RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho, MP_PINFl
  real, dimension(1:nof_SPECIES) :: RGS_x, RGS_Y
  real, dimension(1:nof_SPECIES) :: RGS_mu_i, RGS_ktr_i
  real, dimension(1:nof_SPECIES) :: RGS_htr_i, RGS_hvib_i
  real, dimension(1:nof_SPECIES,1:nof_SPECIES) :: RGS_Dij
  real, dimension(1:nof_SPECIES) :: RGS_Deff
  real :: RGS_ktr_mix, RGS_kve
  real :: sumx
  integer :: rg_i

  ! Convert conserved → primitive state consistently
  prim = LEFTV
  call CONS2PRIM(N, prim, MP_PINFl, GAMMAL)

  ! Extract state variables
  RGS_rho = prim(1)
  RGS_Ttr = prim(dimensiona+2)
  RGS_Tv  = prim(dimensiona+3)
  RGS_P_pa = prim(dimensiona+4)     ! <-- ENSURE your CONS2PRIM places P here

  RGS_Y(:) = prim(nof_variables-nof_species+1 : nof_variables)

  ! --------- FIXED MOLE FRACTION COMPUTATION ----------
  do rg_i = 1, nof_species
     RGS_x(rg_i) = RGS_Y(rg_i) / RG_MOLM(rg_i)
  end do

  sumx = sum(RGS_x)
  sumx = max(sumx, 1d-300)
  RGS_x = RGS_x / sumx      ! normalize properly

  ! --------- CALL TRANSPORT / THERMO MODULE ----------
  call COMPUTE_REAL_GAS_DIFFUSION( RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho, RGS_x, RGS_Y, &
                                   RGS_Dij, RGS_Deff, RGS_mu_i, MP_mu_mix, &
                                   RGS_ktr_i, MP_ktr_mix, MP_kve, RGS_htr_i, RGS_hvib_i )

  MP_D_eff = RGS_Deff
  MP_HTR   = RGS_htr_i
  MP_HVIB  = RGS_hvib_i

END SUBROUTINE MULTISPECIES_MIXTURES_RG



subroutine compute_effective_diffusion(RGS_x, RGS_Dij, RGS_Deff)
  implicit none
  real, intent(in)  :: RGS_x(5)
  real, intent(in)  :: RGS_Dij(5,5)
  real, intent(out) :: RGS_Deff(5)

  integer :: i, j
  real :: sumj, one_minus_xi
  real, parameter :: tiny = 1d-20

  do i = 1, 5
     sumj = 0.0d0

     do j = 1, 5
        if (j /= i) then
           sumj = sumj + RGS_x(j) / max(RGS_Dij(i,j), tiny)
        end if
     end do

     one_minus_xi = max(1.0d0 - RGS_x(i), 0.0d0)

     ! Correct physical behavior:
     if (sumj <= tiny) then
         ! Nearly pure gas: D_i = typical self-diffusion limit
         RGS_Deff(i) = 1d-4   !! physically realistic fallback
     else
         RGS_Deff(i) = one_minus_xi / sumj
     end if

     ! More appropriate lower bound
     RGS_Deff(i) = max(RGS_Deff(i), 1d-10)
  end do

end subroutine compute_effective_diffusion




subroutine htr_air5(RGS_Ttr, RGS_htr_i)
  implicit none
  real, intent(in)  :: RGS_Ttr
  real, intent(out) :: RGS_htr_i(5)
  real :: Cp(5), h0_mass(5)
  integer :: i
  real, parameter :: Tref = 298.15d0

  call cp_tr_species(Cp)

  do i = 1, 5
     h0_mass(i)   = RG_hzero(i) / RG_MOLM(i)
     RGS_htr_i(i) = h0_mass(i) + Cp(i) * (RGS_Ttr - Tref)
  end do

end subroutine htr_air5


subroutine vibrational_conductivity(RGS_rho, RGS_Y, RGS_Deff, RGS_Tv, RGS_kve)
  implicit none

  real, intent(in)  :: RGS_rho, RGS_Y(5), RGS_Deff(5), RGS_Tv
  real, intent(out) :: RGS_kve

  real :: Cv_vib_mass(5)
  integer :: i

  do i = 1, 5
     if (RG_thetaG(i) > 0.0d0) then
        Cv_vib_mass(i) = RGS_cv_vibrational_diatomic(RGS_Tv, RG_thetaG(i)) / RG_MOLM(i)
     else
        Cv_vib_mass(i) = 0.0d0
     end if
  end do

  RGS_kve = 0.0d0
  do i = 1, 5
     RGS_kve = RGS_kve + RGS_rho * RGS_Y(i) * RGS_Deff(i) * Cv_vib_mass(i)
  end do

  RGS_kve = max(RGS_kve, 0.0d0)

end subroutine vibrational_conductivity



subroutine COMPUTE_REAL_GAS_DIFFUSION( RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho, RGS_x, RGS_Y, &
                                       RGS_Dij, RGS_Deff, RGS_mu_i, RGS_mu_mix, &
                                       RGS_ktr_i, RGS_ktr_mix, RGS_kve, &
                                       RGS_htr_i, RGS_hvib_i )

  implicit none
  real, intent(in)  :: RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho
  real, intent(in)  :: RGS_x(5), RGS_Y(5)
  real, intent(out) :: RGS_Dij(5,5), RGS_Deff(5)
  real, intent(out) :: RGS_mu_i(5), RGS_mu_mix
  real, intent(out) :: RGS_ktr_i(5), RGS_ktr_mix, RGS_kve
  real, intent(out) :: RGS_htr_i(5), RGS_hvib_i(5)

  ! ---- Binary diffusion ----
  call compute_binary_diffusion( RGS_Ttr, RGS_P_pa, RGS_Dij )

  ! ---- Mixture-average diffusion coefficient ----
  call compute_effective_diffusion( RGS_x, RGS_Dij, RGS_Deff )

  ! ---- Species viscosity (Blottner) ----
  call blottner_mu_species( RGS_Ttr, RGS_mu_i )

  ! ---- Mixture viscosity (Wilke) ----
  call wilke_mixture_viscosity( RGS_x, RGS_mu_i, RGS_mu_mix )

  ! ---- Species thermal conductivity (Eucken) ----
  call eucken_ktr_species( RGS_mu_i, RGS_ktr_i )

  ! ---- Mixture translational conductivity ----
  call mason_saxena_ktr_mixture( RGS_x, RGS_ktr_i, RGS_ktr_mix )

  ! ---- Vibrational conductivity ----
  call vibrational_conductivity( RGS_rho, RGS_Y, RGS_Deff, RGS_Tv, RGS_kve )

  ! ---- Species translational and vibrational enthalpy ----
  call htr_air5 ( RGS_Ttr, RGS_htr_i )
  call hvib_air5( RGS_Tv , RGS_hvib_i )

end subroutine COMPUTE_REAL_GAS_DIFFUSION



subroutine compute_binary_diffusion(RGS_T, RGS_P_pa, RGS_Dij)
  implicit none
  real, intent(in)  :: RGS_T, RGS_P_pa
  real, intent(out) :: RGS_Dij(5,5)

  integer :: i, j
  real :: sig_ij, eps_ij, Tstar, omega, P_atm, denom
  real, parameter :: tiny = 1d-30, Tmin = 200d0

  RGS_Dij = 0.0d0

  P_atm = max(RGS_P_pa / RGS_Pa_per_atm, 1d-12)
  Tstar = max(RGS_T, Tmin)

  do i = 1, 5
     do j = i+1, 5

        sig_ij = 0.5d0 * (RGS_sigmaA(i) + RGS_sigmaA(j))
        eps_ij = sqrt( RGS_eps_over_k(i) * RGS_eps_over_k(j) )

        Tstar = max(RGS_T / eps_ij, 1d-6)
        omega = omega11_neufeld(Tstar)

        denom = P_atm * sig_ij**2 * omega * sqrt(1.d0/RGS_Mg(i) + 1.d0/RGS_Mg(j))
        denom = max(denom, tiny)

        RGS_Dij(i,j) = (0.001858d0 * RGS_T**1.5d0) / denom * RGS_cm2s_to_m2s
        RGS_Dij(j,i) = RGS_Dij(i,j)

     end do
  end do

end subroutine compute_binary_diffusion


subroutine blottner_mu_species(RGS_T, RGS_mu)
  implicit none
  real, intent(in)  :: RGS_T
  real, intent(out) :: RGS_mu(5)
  integer :: i
  real :: Tlog

  Tlog = log10(max(RGS_T, 50.d0))

  do i = 1, 5
     RGS_mu(i) = 1.d-7 * 10.d0**( RGS_aB(i)*Tlog*Tlog + RGS_bB(i)*Tlog + RGS_cB(i) )
     RGS_mu(i) = max(RGS_mu(i), 1d-12)
  end do
end subroutine blottner_mu_species


subroutine wilke_mixture_viscosity(RGS_x, RGS_mu_i, RGS_mu_mix)
  implicit none
  real, intent(in)  :: RGS_x(5), RGS_mu_i(5)
  real, intent(out) :: RGS_mu_mix
  integer :: i, j
  real :: phi_ij, denom

  RGS_mu_mix = 0.d0

  do i = 1, 5
     denom = 0.d0

     do j = 1, 5
        if (i == j) then
            denom = denom + RGS_x(j)
        else
            phi_ij = ( 1.d0 + sqrt(RGS_mu_i(i)/RGS_mu_i(j)) * (RG_MOLM(j)/RG_MOLM(i))**0.25 )**2 &
                     / ( sqrt(8.d0) * sqrt(1.d0 + RG_MOLM(i)/RG_MOLM(j)) )

            denom = denom + RGS_x(j) * phi_ij
        end if
     end do

     RGS_mu_mix = RGS_mu_mix + RGS_x(i) * RGS_mu_i(i) / max(denom, 1d-20)
  end do

end subroutine wilke_mixture_viscosity


subroutine eucken_ktr_species(RGS_mu_i, RGS_ktr_i)
  implicit none
  real, intent(in)  :: RGS_mu_i(5)
  real, intent(out) :: RGS_ktr_i(5)

  integer :: i
  real :: Rspec, Cp

  do i = 1, 5
     Rspec = RGS_Ru / RG_MOLM(i)

     if (i <= 3) then
         Cp = 3.5d0 * Rspec          ! diatomic
         RGS_ktr_i(i) = RGS_mu_i(i) * (Cp + 1.25d0*Rspec)
     else
         Cp = 2.5d0 * Rspec          ! atomic
         RGS_ktr_i(i) = RGS_mu_i(i) * (Cp + 1.50d0*Rspec)
     end if

  end do

end subroutine eucken_ktr_species


subroutine mason_saxena_ktr_mixture(RGS_x, RGS_ktr_i, RGS_ktr_mix)
  implicit none
  real, intent(in)  :: RGS_x(5), RGS_ktr_i(5)
  real, intent(out) :: RGS_ktr_mix

  integer :: i, j
  real :: psi_ij, denom

  RGS_ktr_mix = 0.d0

  do i = 1, 5
     denom = 0.d0

     do j = 1, 5
        if (i == j) then
            denom = denom + RGS_x(j)
        else
            psi_ij = (1.d0 + sqrt(RGS_ktr_i(i)/RGS_ktr_i(j)) * (RG_MOLM(j)/RG_MOLM(i))**0.25 )**2 &
                      / ( sqrt(8.d0) * sqrt(1.d0 + RG_MOLM(i)/RG_MOLM(j)) )

            denom = denom + RGS_x(j) * psi_ij
        end if
     end do

     RGS_ktr_mix = RGS_ktr_mix + RGS_x(i) * RGS_ktr_i(i) / max(denom, 1d-20)

  end do

end subroutine mason_saxena_ktr_mixture






