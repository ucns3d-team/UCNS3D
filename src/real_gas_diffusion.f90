module transport_twoT_air5
  implicit none
  private
  public :: air5_compute_all, air5_species_names
  public :: air5_fill_blottner   ! helper to set your a,b,c at runtime if you prefer

  integer, parameter :: dp = selected_real_kind(15, 300)
  real(dp), parameter :: Ru = 8.31446261815324_dp          ! J/mol-K
  real(dp), parameter :: Pa_per_atm = 101325.0_dp
  real(dp), parameter :: cm2s_to_m2s = 1.0e-4_dp
  real(dp), parameter :: tiny = 1.0e-30_dp

  ! Species ordering
  character(len=*), parameter :: air5_species_names(5) = (/ 'N2', 'O2', 'NO', 'N ', 'O ' /)

  ! Molar masses [kg/mol]
  real(dp), parameter :: M(5) = (/ 28.0134e-3_dp, 31.9988e-3_dp, 30.0061e-3_dp, 14.0067e-3_dp, 15.9994e-3_dp /)
  ! Same in g/mol for the standard D_ij correlation
  real(dp), parameter :: Mg(5) = M * 1.0e3_dp

  ! Lennard–Jones parameters (typical Gupta–Yos / Hirschfelder set; adjust if your database differs)
  ! sigma [Angstrom], epsilon/k [K]
  real(dp), parameter :: sigmaA(5)     = (/ 3.621_dp, 3.458_dp, 3.491_dp, 3.298_dp, 3.050_dp /)
  real(dp), parameter :: eps_over_k(5) = (/ 97.53_dp, 106.70_dp, 116.70_dp, 71.40_dp, 80.00_dp /)

  ! Diatomic vibrational characteristic temperatures [K] (atoms have 0)
  ! N2: 3390, O2: 2270, NO: 2740 (common values)
  real(dp), parameter :: theta_v(5) = (/ 3390._dp, 2270._dp, 2740._dp, 0._dp, 0._dp /)

  ! ---- Blottner viscosity coefficients (YOU MUST CONFIRM / SUPPLY) ----
  ! Form used here:  mu_i [Pa·s] = 1e-7 * 10**( a_i*(log10 T)^2 + b_i*log10 T + c_i ), T in K
  ! Provide {a,b,c} for [N2, O2, NO, N, O]. Values vary by source; set yours below.
  real(dp) :: aB(5), bB(5), cB(5)
  logical  :: blottner_initialized = .false.
  ! ---------------------------------------------------------------------

contains

  pure function omega11_neufeld(Tstar) result(omega11)
    real(dp), intent(in) :: Tstar
    real(dp) :: omega11
    real(dp), parameter :: A=1.06036_dp, B=0.15610_dp, C=0.19300_dp, D=0.47635_dp, &
                           E=1.03587_dp, F=1.52996_dp, G=1.76474_dp, H=3.89411_dp
    omega11 = A / (Tstar**B) + C*exp(-D*Tstar) + E*exp(-F*Tstar) + G*exp(-H*Tstar)
  end function omega11_neufeld

  subroutine compute_binary_diffusion(T, P_pa, Dij)
    real(dp), intent(in)  :: T, P_pa
    real(dp), intent(out) :: Dij(5,5)
    integer :: i, j
    real(dp) :: sig_ij, epsk_ij, Tstar, omega, P_atm, denom
    Dij = 0.0_dp
    P_atm = max(P_pa / Pa_per_atm, 1.0e-12_dp)
    do i = 1, 5
      do j = i+1, 5
        sig_ij  = 0.5_dp*(sigmaA(i) + sigmaA(j))
        epsk_ij = sqrt(eps_over_k(i) * eps_over_k(j))
        Tstar   = max(T / epsk_ij, 1.0e-8_dp)
        omega   = omega11_neufeld(Tstar)
        denom   = P_atm * (sig_ij**2) * omega * sqrt( 1.0_dp/Mg(i) + 1.0_dp/Mg(j) )
        Dij(i,j) = (0.001858_dp * T**1.5_dp) / max(denom,tiny) * cm2s_to_m2s
        Dij(j,i) = Dij(i,j)
      end do
    end do
  end subroutine compute_binary_diffusion

  subroutine compute_effective_diffusion(x, Dij, Deff)
    real(dp), intent(in)  :: x(5)
    real(dp), intent(in)  :: Dij(5,5)
    real(dp), intent(out) :: Deff(5)
    integer :: i, j
    real(dp) :: sumj
    do i = 1, 5
      sumj = 0.0_dp
      do j = 1, 5
        if (j /= i) sumj = sumj + x(j) / max(Dij(i,j), tiny)
      end do
      Deff(i) = merge( 1.0_dp/sumj, 0.0_dp, sumj>0.0_dp )
    end do
  end subroutine compute_effective_diffusion

  pure subroutine blottner_mu_species(T, mu)
    real(dp), intent(in)  :: T
    real(dp), intent(out) :: mu(5)
    integer :: i
    real(dp) :: lt
    if (.not. blottner_initialized) call set_default_blottner_(&
         'WARNING: using placeholder Blottner coefficients — replace with your dataset.')
    lt = log10(max(T,1.0_dp))
    do i = 1, 5
      mu(i) = 1.0e-7_dp * 10.0_dp**( aB(i)*lt*lt + bB(i)*lt + cB(i) )
    end do
  end subroutine blottner_mu_species

  pure subroutine wilke_mixture_viscosity(x, mu_i, mu_mix)
    real(dp), intent(in)  :: x(5), mu_i(5)
    real(dp), intent(out) :: mu_mix
    integer :: i, j
    real(dp) :: phi_ij, denom
    mu_mix = 0.0_dp
    do i = 1, 5
      denom = 0.0_dp
      do j = 1, 5
        if (i == j) then
          denom = denom + x(j)
        else
          phi_ij = ( 1.0_dp + sqrt(mu_i(i)/mu_i(j)) * (M(j)/M(i))**0.25_dp )**2 &
                   / ( sqrt(8.0_dp) * sqrt(1.0_dp + M(i)/M(j)) )
          denom = denom + x(j) * phi_ij
        end if
      end do
      mu_mix = mu_mix + x(i) * mu_i(i) / max(denom, tiny)
    end do
  end subroutine wilke_mixture_viscosity

  pure subroutine eucken_ktr_species(mu_i, ktr_i)   !THIS IS ALREADY IMPLEMENTED WE NEED TO NOTICE THE I<3
    real(dp), intent(in)  :: mu_i(5)
    real(dp), intent(out) :: ktr_i(5)
    integer :: i
    real(dp) :: Rspec, Cp_tr
    do i = 1, 5
      ! Cp_tr: atoms 5/2 R, diatomics (trans+rot fully excited) 7/2 R
      Rspec = Ru / M(i)
      if (i <= 3) then
        Cp_tr = 3.5_dp * Rspec   ! N2, O2, NO
      else
        Cp_tr = 2.5_dp * Rspec   ! N, O
      end if
      ktr_i(i) = mu_i(i) * ( Cp_tr + 1.25_dp * Rspec )
    end do
  end subroutine eucken_ktr_species

  pure subroutine mason_saxena_ktr_mixture(x, ktr_i, ktr_mix)
    real(dp), intent(in)  :: x(5), ktr_i(5)
    real(dp), intent(out) :: ktr_mix
    integer :: i, j
    real(dp) :: psi_ij, denom
    ktr_mix = 0.0_dp
    do i = 1, 5
      denom = 0.0_dp
      do j = 1, 5
        if (i == j) then
          denom = denom + x(j)
        else
          psi_ij = ( 1.0_dp + sqrt(ktr_i(i)/ktr_i(j)) * (M(j)/M(i))**0.25_dp )**2 &
                   / ( sqrt(8.0_dp) * sqrt(1.0_dp + M(i)/M(j)) )
          denom = denom + x(j) * psi_ij
        end if
      end do
      ktr_mix = ktr_mix + x(i) * ktr_i(i) / max(denom, tiny)
    end do
  end subroutine mason_saxena_ktr_mixture

  pure function cv_vibrational_diatomic(Tv, theta) result(cv)
    ! Vibrational cv of a single diatomic mode (no anharmonicity/electronic)
    ! cv = R * (theta/T)^2 * exp(theta/T) / (exp(theta/T)-1)^2
    real(dp), intent(in) :: Tv, theta
    real(dp) :: cv, x
    if (Tv <= 1.0_dp .or. theta <= 0.0_dp) then
      cv = 0.0_dp
    else
      x = theta / Tv
      cv = Ru * (x*x) * exp(x) / ( (exp(x) - 1.0_dp)**2 )
    end if
  end function cv_vibrational_diatomic

  pure subroutine vibrational_conductivity(rho, Y, Deff, Tv, kve)
    real(dp), intent(in)  :: rho, Y(5), Deff(5), Tv
    real(dp), intent(out) :: kve
    real(dp) :: cv_ve(5), Rspec
    integer :: i
    ! Build cv_ve: diatomics get vibrational cv; atoms get 0 (electronic omitted here)
    do i = 1, 5
      if (theta_v(i) > 0.0_dp) then
        cv_ve(i) = cv_vibrational_diatomic(Tv, theta_v(i)) / M(i)   ! J/kg-K
      else
        cv_ve(i) = 0.0_dp
      end if
    end do
    kve = 0.0_dp
    do i = 1, 5
      kve = kve + rho * Y(i) * Deff(i) * cv_ve(i)
    end do
  end subroutine vibrational_conductivity

  ! ---------------- High-level one-shot API ----------------
  subroutine air5_compute_all(Ttr, Tv, P_pa, rho, x, Y, &
                              Dij, Deff, mu_i, mu_mix, ktr_i, ktr_mix, kve)
    real(dp), intent(in)  :: Ttr, Tv, P_pa, rho
    real(dp), intent(in)  :: x(5), Y(5)
    real(dp), intent(out) :: Dij(5,5), Deff(5), mu_i(5), mu_mix
    real(dp), intent(out) :: ktr_i(5), ktr_mix, kve
    call compute_binary_diffusion(Ttr, P_pa, Dij)
    call compute_effective_diffusion(x, Dij, Deff)
    call blottner_mu_species(Ttr, mu_i)
    call wilke_mixture_viscosity(x, mu_i, mu_mix)
    call eucken_ktr_species(mu_i, ktr_i)
    call mason_saxena_ktr_mixture(x, ktr_i, ktr_mix)
    call vibrational_conductivity(rho, Y, Deff, Tv, kve)
  end subroutine air5_compute_all

  ! ---------------- Convenience: set/get Blottner ----------------
  subroutine air5_fill_blottner(a_in, b_in, c_in)
    real(dp), intent(in) :: a_in(5), b_in(5), c_in(5)
    aB = a_in;  bB = b_in;  cB = c_in
    blottner_initialized = .true.
  end subroutine air5_fill_blottner

  subroutine set_default_blottner_(msg)
    character(len=*), intent(in), optional :: msg
    ! PLACEHOLDER values to keep code runnable. Replace with your vetted table!
    ! >>> BEGIN: replace these with your Blottner coefficients <<<
    ! Example-only numbers (NOT authoritative) for demonstration:
    aB = (/ 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp /)
    bB = (/ 0.9420_dp, 0.9850_dp, 0.9600_dp, 0.8500_dp, 0.9000_dp /)
    cB = (/ -7.2200_dp, -7.2500_dp, -7.2300_dp, -7.0000_dp, -7.0500_dp /)
    ! >>> END: replace <<<
    blottner_initialized = .true.
    if (present(msg)) then
      ! no I/O inside pure routines; this sub is non-pure on purpose—silence in library context
    end if
  end subroutine set_default_blottner_

end module transport_twoT_air5



call air5_fill_blottner( (/0.0617092805_dp, 0.1033860707_dp, 0.1003927101_dp, 0.0267099871_dp, 0.0467424774_dp/), &
                         (/0.3180000000_dp, -0.0826000000_dp, -0.0336000000_dp, 0.6030000000_dp, 0.4290000000_dp/), &
                         (/1.0924723545_dp,  2.0044907665_dp,  1.8394588634_dp, 0.6147484244_dp, 0.9621840099_dp/) )




program demo_air5
  use transport_twoT_air5
  implicit none
  integer, parameter :: dp = selected_real_kind(15, 300)
  real(dp) :: Ttr, Tv, P, rho
  real(dp) :: x(5), Y(5)
  real(dp) :: Dij(5,5), Deff(5), mu_i(5), mu_mix, ktr_i(5), ktr_mix, kve
  real(dp) :: a(5), b(5), c(5)

  ! Example flow state
  Ttr = 8000._dp
  Tv  = 6000._dp
  P   = 5000._dp      ! Pa
  rho = 0.02_dp       ! kg/m^3

  ! Example composition (must be consistent x & Y with M; set either from your solver)
  x = (/ 0.6_dp, 0.2_dp, 0.1_dp, 0.05_dp, 0.05_dp /)
  Y = x * ( (/28.0134_dp,31.9988_dp,30.0061_dp,14.0067_dp,15.9994_dp/)*1e-3_dp ) &
        / sum( x * ( (/28.0134_dp,31.9988_dp,30.0061_dp,14.0067_dp,15.9994_dp/)*1e-3_dp ) )

  ! Supply your Blottner table here (example placeholders below!)
  a = (/ 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp, 0.0000_dp /)
  b = (/ 0.9420_dp, 0.9850_dp, 0.9600_dp, 0.8500_dp, 0.9000_dp /)
  c = (/ -7.2200_dp,-7.2500_dp,-7.2300_dp,-7.0000_dp,-7.0500_dp /)
  call air5_fill_blottner(a,b,c)

  call air5_compute_all(Ttr, Tv, P, rho, x, Y, Dij, Deff, mu_i, mu_mix, ktr_i, ktr_mix, kve)

  ! ...use Dij/Deff/mu_i/mu_mix/ktr_i/ktr_mix/kve in your fluxes/closures ...
end program demo_air5









