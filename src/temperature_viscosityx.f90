module temp_from_energy_janaf5_X
  implicit none
  private
  public :: set_thermo_tables, T_from_conserved_X

  integer, parameter :: dp = kind(1.0d0)
  integer :: NS = 0
  real(dp), allocatable :: MP_M(:)                 ! kg/mol
  real(dp), allocatable :: A(:,:,:)              ! JANAF a(5,2,NS)
  real(dp), allocatable :: Tlow(:), Thigh(:), Tmid(:)
  logical :: include_hf0 = .false.               ! .false. => sensible energies (common)

contains

  subroutine set_thermo_tables(NOF_SPECIES, mw_in, a_in, Tlow_in, Thigh_in, Tmid_in, use_hf0)
    integer, intent(in) :: NOF_SPECIES
    real(dp), intent(in) :: mw_in(NOF_SPECIES)
    real(dp), intent(in) :: a_in(5,2,NOF_SPECIES)
    real(dp), intent(in) :: Tlow_in(NOF_SPECIES), Thigh_in(NOF_SPECIES), Tmid_in(NOF_SPECIES)
    logical,  intent(in), optional :: use_hf0
    NS = NOF_SPECIES
    allocate(MP_M(NS), A(5,2,NS), Tlow(NS), Thigh(NS), Tmid(NS))
    MW = mw_in; A = a_in; Tlow=Tlow_in; Thigh=Thigh_in; Tmid=Tmid_in
    if (present(use_hf0)) include_hf0 = use_hf0
  end subroutine set_thermo_tables

  pure subroutine X_to_Y(X, Y)
    real(dp), intent(in)  :: X(NS)
    real(dp), intent(out) :: Y(NS)
    integer :: s
    real(dp) :: denom
    denom = 0.0_dp
    do s=1,NS; denom = denom + X(s)*MP_M(s); end do
    if (denom <= 0.0_dp) then
      Y = 0.0_dp
    else
      do s=1,NS; Y(s) = X(s)*MP_M(s)/denom; end do
    end if
  end subroutine X_to_Y

  pure subroutine species_cp_h(T, is, cp, h)
    ! JANAF 5-coeffs:
    !   cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^-2
    !    h/R = a1*T + a2*T^2/2 + a3*T^3/3 + a4*T^4/4 - a5/T + C
    real(dp), intent(in)  :: T
    integer,  intent(in)  :: is
    real(dp), intent(out) :: cp, h
    integer :: r
    real(dp) :: Ru, Rs, t,t2,t3, invT, invT2, a1,a2,a3,a4,a5, href
    Ru = 8.31446261815324_dp; Rs = Ru / MP_M(is)
    r = merge(1,2, T <= Tmid(is))
    a1=A(1,r,is); a2=A(2,r,is); a3=A(3,r,is); a4=A(4,r,is); a5=A(5,r,is)
    t=T; t2=t*t; t3=t2*t; invT=1.0_dp/t; invT2=invT*invT
    cp = Rs * (a1 + a2*t + a3*t2 + a4*t3 + a5*invT2)
    h  = Rs * (a1*t + 0.5_dp*a2*t2 + (a3/3.0_dp)*t3 + 0.25_dp*a4*t2*t2 - a5*invT)
    if (.not. include_hf0) then
      href = Rs * (a1*298.15_dp + 0.5_dp*a2*298.15_dp**2 + (a3/3.0_dp)*298.15_dp**3 + &
                   0.25_dp*a4*298.15_dp**4 - a5/298.15_dp)
      h = h - href
    end if
  end subroutine species_cp_h

  pure subroutine mix_props_T_X(T, X, Rmix, cp_mix, cv_mix, e_mix)
    ! Mixture via mass-weighted energy with X->Y internally (energy is mass-specific).
    real(dp), intent(in)  :: T
    real(dp), intent(in)  :: X(NS)
    real(dp), intent(out) :: Rmix, cp_mix, cv_mix, e_mix
    real(dp) :: Y(NS), Ru, cp_s, h_s
    integer :: s
    call X_to_Y(X, Y)
    Ru = 8.31446261815324_dp
    Rmix=0.0_dp; cp_mix=0.0_dp; e_mix=0.0_dp
    do s=1,NS
      if (Y(s) <= 0.0_dp) cycle
      Rmix = Rmix + Y(s) * (Ru/MP_M(s))
      call species_cp_h(T, s, cp_s, h_s)
      e_mix  = e_mix  + Y(s) * (h_s - (Ru/MP_M(s))*T) ! e = h - Rs*T
      cp_mix = cp_mix + Y(s) * cp_s
    end do
    cv_mix = cp_mix - Rmix
  end subroutine mix_props_T_X

  subroutine T_from_conserved_X(rho, u, v, w, rhoE, X, T, info, T_low, T_high, tol, max_iter)
    ! Recover T from conservative state, using mole fractions X.
    real(dp), intent(in)  :: rho, u, v, w, rhoE
    real(dp), intent(in)  :: X(NS)    ! mole fractions (sum~1)
    real(dp), intent(out) :: T
    integer,  intent(out) :: info     ! 0 ok; 1 width tol; 2 iters; <0 bracket issue
    real(dp), intent(in), optional :: T_low, T_high, tol
    integer,  intent(in), optional :: max_iter

    real(dp) :: kin, e_tgt, Rmix, cp_mix, cv_mix, e_mix
    real(dp) :: Tlo, Thi, fa, fb, f, df, aT, bT, Tn, fn, Tguess, atol
    integer  :: it, itmax, s
    logical  :: bracket

    itmax = merge(max_iter, 60, present(max_iter))
    kin   = 0.5_dp*(u*u + v*v + w*w)
    e_tgt = (rhoE/rho) - kin                          ! sensible target [J/kg]

    ! Species-aware bracket
    Tlo = merge(T_low, 150.0_dp, present(T_low))
    Thi = merge(T_high, 20000.0_dp, present(T_high))
    do s=1,NS
      if (X(s) <= 0.0_dp) cycle
      Tlo = max(Tlo, Tlow(s));  Thi = min(Thi, Thigh(s))
    end do
    if (Tlo >= Thi) then
      info = -3; T = 300.0_dp; return
    end if

    ! Initial guess via cv at 300 K
    call mix_props_T_X(300.0_dp, X, Rmix, cp_mix, cv_mix, e_mix)
    Tguess = 300.0_dp + (e_tgt - e_mix)/max(1.0e-8_dp, cv_mix)
    T = min(max(Tguess, Tlo), Thi)

    ! Bracket residuals
    call mix_props_T_X(Tlo, X, Rmix, cp_mix, cv_mix, e_mix); fa = e_mix - e_tgt
    call mix_props_T_X(Thi, X, Rmix, cp_mix, cv_mix, e_mix); fb = e_mix - e_tgt
    bracket = (fa*fb <= 0.0_dp)
    aT = Tlo; bT = Thi

    do it = 1, itmax
      call mix_props_T_X(T, X, Rmix, cp_mix, cv_mix, e_mix)
      f  = e_mix - e_tgt
      df = cv_mix

      atol = merge(tol, 1.0e-6_dp*max(1.0_dp,abs(e_tgt)) + 1.0e-3_dp, present(tol))
      if (abs(f) <= atol) then
        info = 0; return
      end if

      ! Safeguarded Newton
      if (df > 1.0e-20_dp) then
        Tn = T - f/df
      else
        Tn = 0.5_dp*(aT + bT)
      end if
      if ((.not. bracket) .or. (Tn <= aT) .or. (Tn >= bT)) Tn = 0.5_dp*(aT + bT)

      call mix_props_T_X(Tn, X, Rmix, cp_mix, cv_mix, e_mix)
      fn = e_mix - e_tgt

      if (fa*fn <= 0.0_dp) then
        bT = Tn; fb = fn; bracket = .true.
      else
        aT = Tn; fa = fn
      end if
      T = Tn

      if (abs(bT - aT) <= 1.0e-8_dp*max(1.0_dp,T)) then
        info = 1; return
      end if
    end do
    info = 2
  end subroutine T_from_conserved_X

end module temp_from_energy_janaf5_X
