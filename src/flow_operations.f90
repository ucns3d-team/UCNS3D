module flow_operations
   use declaration
   use mpiinfo
   use transform

   implicit none

contains

   subroutine fill_turbulence_gradients(nd, i, packed)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: nd, i
      real, intent(inout) :: packed(20)
      integer :: v, d, base
      do v = 1, nd
         do d = 1, nd
            packed(3 + (v - 1) * nd + d) = rec_grads(i, v, d)
         end do
      end do
      base = 3 + nd * nd
      do v = 1, 2
         do d = 1, nd
            packed(base + (v - 1) * nd + d) = rec_grads(i, nof_variables - 1 + v, d)
         end do
      end do
   end subroutine fill_turbulence_gradients

   real function translational_cp(q) result(cp)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real :: y(gpu_max_species), temp, rmix, cv, energy
      integer :: k
      cp = gamma * r_gas / (gamma - 1.0e0)
      if (realgas == 1) then
         cp = 0.0e0
         do k = 1, nof_species
            cp = cp + (q(dimensiona + 3 + k) / q(1)) * merge(3.5e0, 2.5e0, k <= 3) * rgs_ru / rg_molm(k)
         end do
      elseif (multispecies == 1) then
         call multispecies_thermo(q, temp, y, rmix, cp, cv, energy)
      end if
   end function translational_cp

   subroutine eddy_transport_pair(nd, ql, qr, mu, conductivity, turbmv, etvm, el, er)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: nd
      real, intent(in) :: ql(nof_variables), qr(nof_variables)
      real, intent(inout) :: mu(4), conductivity(4), turbmv(2), etvm(1), el(20), er(20)
      real :: skl, swl, skr, swr
      mu(3:4) = 0.0e0; conductivity(3:4) = 0.0e0
      if (turbulence == 0) return
      turbmv = [el(2), er(2)]
      call eddy_viscosity_state(nd, ql(1), mu(1), el, .true., mu(3), skl, swl)
      if (turbulencemodel == 2 .and. er(1) <= 0.0e0) then
         mu(4) = 0.0e0; skr = skl; swr = swl
      else
         call eddy_viscosity_state(nd, qr(1), mu(2), er, .false., mu(4), skr, swr)
      end if
      if (turbulencemodel == 1 .and. min(turbmv(1), turbmv(2)) < 0.0e0) mu(3:4) = 0.0e0
      ! Slots 3:4 contain total conductivity, with the EOS-specific molecular
      ! values preserved. Vibrational conduction is handled by the real-gas flux.
      conductivity(3) = conductivity(1) + mu(3) * translational_cp(ql) / prtu
      conductivity(4) = conductivity(2) + mu(4) * translational_cp(qr) / prtu
      etvm(1) = 0.5e0 * sum(mu)
      if (turbulencemodel == 2) then
         el(19) = mu(1) + mu(3) / skl; er(19) = mu(2) + mu(4) / skr
         el(20) = mu(1) + mu(3) / swl; er(20) = mu(2) + mu(4) / swr
      end if
   end subroutine eddy_transport_pair

   subroutine eddy_viscosity_state(nd, rho, mu, e, left_side, mut, sigk, sigw)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: nd
      real, intent(in) :: rho, mu, e(20)
      logical, intent(in) :: left_side
      real, intent(out) :: mut, sigk, sigw
      real :: chi, chi3, k0, om, wall, strain, sij, cross, dplus, phi1, phi2, f1, f2, ret, beta0, astar
      integer :: i, j, kg, wg
      mut = 0.0e0; sigk = 1.0e0; sigw = 1.0e0
      select case(turbulencemodel)
      case(1)
         chi = abs(max(e(2), 1.0e-15) / max(mu, 1.0e-15))
         chi3 = chi * chi * chi
         mut = e(2) * chi3 / (chi3 + cv1**3)
      case(2)
         strain = 0.0e0
         do i = 1, nd
            do j = 1, nd
               sij = 0.5e0 * (e(3 + (i - 1) * nd + j) + e(3 + (j - 1) * nd + i))
               strain = strain + sij * sij
            end do
         end do
         strain = sqrt(2.0e0 * strain)
         k0 = max(e(2) / rho, 1.0e-15)
         om = max(e(3) / rho, 1.0e-6)
         if (left_side) om = max(e(3) / rho, ufreestream / charlength / 10.0e0, 1.0e-6)
         wall = max(e(1), 1.0e-15)
         kg = 3 + nd * nd; wg = kg + nd
         cross = 0.0e0
         do j = 1, nd
            cross = cross + e(kg + j) * e(wg + j)
         end do
         dplus = max(2.0e0 * rho / sigma_om2 / om * cross, 1.0e-10)
         phi2 = max(sqrt(k0) / (0.09e0 * om * wall), 500.0e0 * mu / (rho * wall * wall * om))
         phi1 = min(phi2, 4.0e0 * rho * k0 / (sigma_om2 * dplus * wall * wall))
         f1 = tanh(phi1**4); f2 = tanh(phi2**2)
         ret = rho * k0 / (max(mu, 1.0e-30) * om)
         beta0 = f1 * beta_i1 + (1.0e0 - f1) * beta_i2
         astar = alpha_starinf * (beta0 / 3.0e0 + ret / r_k_sst) / (1.0e0 + ret / r_k_sst)
         mut = rho * k0 / om / max(1.0e0 / astar, strain * f2 / (aa_1 * om))
         sigk = sigma_k1 / max(f1, 1.0e-30) + sigma_k2 / max(f2, 1.0e-30)
         sigw = sigma_om1 / max(f1, 1.0e-30) + sigma_om2 / max(f2, 1.0e-30)
      end select
      mut = max(0.0e0, min(10000000.0e0 * visc, mut))
   end subroutine eddy_viscosity_state

   subroutine allaire_properties(q, rho, gm, pinf, stiff)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: rho, gm, pinf, stiff
      real :: a, asum, den
      integer :: k
      rho = sum(q(dimensiona + 3:dimensiona + 2 + nof_species))
      asum = 0.0e0; den = 0.0e0; pinf = 0.0e0; stiff = 0.0e0
      do k = 1, nof_species
         if (k < nof_species) then
            a = q(dimensiona + 2 + nof_species + k)
            asum = asum + a
         else
            a = 1.0e0 - asum
         end if
         den = den + a / (gamma_in(k) - 1.0e0)
         pinf = pinf + a * mp_pinf(k)
         stiff = stiff + a * gamma_in(k) / (gamma_in(k) - 1.0e0) * mp_pinf(k)
      end do
      gm = 1.0e0 + 1.0e0 / den
      stiff = stiff * (gm - 1.0e0)
   end subroutine allaire_properties

   subroutine multispecies_thermo(q, temp, y, rmix, cp, cv, energy)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: temp, y(gpu_max_species), rmix, cp, cv, energy
      real :: x(gpu_max_species), total
      integer :: k
      total = sum(q(dimensiona + 3:dimensiona + 2 + nof_species))
      do k = 1, nof_species
         y(k) = q(dimensiona + 2 + k) / total
         x(k) = y(k) / mp_m(k)
      end do
      total = sum(x(1:nof_species))
      x(1:nof_species) = x(1:nof_species) / total
      call get_temp_janaf(q, temp, y, x)
      call mixture_props(temp, x, rmix, cp, cv, energy)
   end subroutine multispecies_thermo
   logical function conservative_state_valid(q) result(valid)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real :: rho, thermal, pressure, gm, pinf, stiff
      integer :: k
      valid = .false.
      if (any(q /= q) .or. any(abs(q) > 1.0e100)) return
      if (q(1) <= tolsmall) return
      if (multispecies == 1 .and. mp_modelc == 1) then
         call multispecies_pressure(q, pressure, valid)
         return
      end if
      rho = q(1)
      if (multispecies == 1) then
         if (any(q(dimensiona + 3:dimensiona + 2 + nof_species) < 0.0e0)) return
         if (any(q(dimensiona + 3 + nof_species:nof_variables) < 0.0e0)) return
         if (sum(q(dimensiona + 3 + nof_species:nof_variables)) > 1.0e0) return
         call allaire_properties(q, rho, gm, pinf, stiff)
         if (rho <= tolsmall) return
      end if
      thermal = q(dimensiona + 2) - 0.5e0 * sum(q(2:dimensiona + 1)**2) / rho
      if (realgas == 1) then
         if (q(dimensiona + 3) < 0.0e0) return
         thermal = thermal - q(dimensiona + 3)
         do k = 1, nof_species
            if (q(dimensiona + 3 + k) < 0.0e0) return
            thermal = thermal - q(dimensiona + 3 + k) * max(rg_hzero(k), 0.0e0) / rg_molm(k)
         end do
         if (sum(q(dimensiona + 4:nof_variables)) <= tolsmall) return
      elseif (multispecies == 1) then
         thermal = (gm - 1.0e0) * thermal - stiff
      end if
      valid = thermal > tolsmall
   end function conservative_state_valid

   subroutine multispecies_pressure(q, pressure, valid)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: pressure
      logical, intent(out) :: valid
      real :: y(gpu_max_species), temp, rmix, cp, cv, energy, target_energy
      valid = .false.; pressure = 0.0e0
      if (any(q /= q) .or. any(abs(q) > 1.0e100)) return
      if (q(1) <= 1.0e-30) return
      if (any(q(dimensiona + 3:nof_variables) < 0.0e0)) return
      if (sum(q(dimensiona + 3:nof_variables)) <= 1.0e-30) return
      call multispecies_thermo(q, temp, y, rmix, cp, cv, energy)
      target_energy = q(dimensiona + 2) / q(1) - 0.5e0 * sum((q(2:dimensiona + 1) / q(1))**2)
      ! A clamped temperature must not silently validate an out-of-range state.
      if (abs(energy - target_energy) > 2.0e-6 * max(abs(target_energy), 1.0e0) + 2.0e-3) return
      pressure = q(1) * rmix * temp
      valid = pressure > 0.0e0 .and. cv > 0.0e0
   end subroutine multispecies_pressure

   subroutine cons2div_realgas(q, pinf, gm)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: q(nof_variables)
      real, intent(out) :: pinf, gm
      real :: rho, ev, etr, rmix, cv, temp, tv, y(gpu_max_species)
      call realgas_thermo(q, rho, y, ev, etr, rmix, cv)
      temp = realgas_ttr(etr, cv)
      tv = realgas_vibrational_temperature(y, ev, temp)
      gm = (cv + rmix) / cv
      pinf = 0.0e0
      q(2:dimensiona + 1) = q(2:dimensiona + 1) / rho
      q(dimensiona + 2) = temp
      q(dimensiona + 3) = tv
      q(dimensiona + 4:nof_variables) = y(1:nof_species)
   end subroutine cons2div_realgas

   subroutine realgas_thermo(q, rho, y, ev, etr, rmix, cv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: rho, y(gpu_max_species), ev, etr, rmix, cv
      real :: ke, echem, rs, velocity
      integer :: k
      ! Density belongs to continuity; never replace it by sum(rho_species).
      rho = q(1)
      ke = 0.0e0
      do k = 2, dimensiona + 1
         velocity = q(k) / rho
         ke = ke + velocity * velocity
      end do
      ke = 0.5e0 * ke
      ev = q(dimensiona + 3) / rho
      rmix = 0.0e0; cv = 0.0e0; echem = 0.0e0
      do k = 1, nof_species
         y(k) = q(dimensiona + 3 + k) / rho
         rs = rgs_ru / rg_molm(k)
         rmix = rmix + y(k) / rg_molm(k)
         cv = cv + y(k) * merge(2.5e0, 1.5e0, k <= 3) * rs
         if (rg_hzero(k) > 0.0e0) echem = echem + y(k) * rg_hzero(k) / rg_molm(k)
      end do
      rmix = rgs_ru * rmix
      etr = q(dimensiona + 2) / rho - ev - echem - ke
      if (etr < 0.0e0) etr = 1.0e-12
   end subroutine realgas_thermo

   real function realgas_ttr(etr, cv) result(temp)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: etr, cv
      temp = 50.0e0
      if (cv > rgs_tiny) temp = etr / cv
      temp = max(50.0e0, min(20000.0e0, temp))
   end function realgas_ttr

   real function realgas_vibrational_temperature(y, ev, ttr) result(tv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: y(gpu_max_species), ev, ttr
      real :: g, df, theta, ei, tmpexp, denom
      integer :: iter, k
      tv = ttr
      do iter = 1, 40
         g = 0.0e0; df = 0.0e0
         do k = 1, 3
            if (y(k) < 1.0e-16) cycle
            theta = rg_thetag(k)
            if (theta <= 0.0e0) cycle
            ei = theta / tv
            if (ei > 60.0e0) then
               tmpexp = exp(-ei)
               g = g + y(k) * (rgs_ru / rg_molm(k)) * theta * tmpexp
               df = df + y(k) * (rgs_ru / rg_molm(k)) * theta * (ei / tv) * tmpexp
            else
               tmpexp = exp(ei)
               denom = tmpexp - 1.0e0
               g = g + y(k) * (rgs_ru / rg_molm(k)) * (theta / denom)
               df = df + y(k) * (rgs_ru / rg_molm(k)) * theta * ei * tmpexp / (denom * denom * tv)
            end if
         end do
         g = g - ev
         if (abs(g) < 1.0e-12 * max(1.0e0, abs(ev))) exit
         if (abs(df) < 1.0e-20) exit
         tv = max(50.0e0, min(20000.0e0, tv - g / df))
      end do
   end function realgas_vibrational_temperature

   subroutine realgas_mole_fractions(y, x)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: y(gpu_max_species)
      real, intent(out) :: x(gpu_max_species)
      real :: sx
      integer :: k
      sx = 0.0e0
      do k = 1, nof_species
         x(k) = y(k) / rg_molm(k)
         sx = sx + x(k)
      end do
      if (sx > 1.0e-300) then
         x(1:nof_species) = x(1:nof_species) / sx
      else
         x(1:nof_species) = 1.0e0 / nof_species
      end if
   end subroutine realgas_mole_fractions

   subroutine realgas_mu_ktr(q, mu, ktr)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: mu, ktr
      real :: rho, ev, etr, rmix, cv
      real :: y(gpu_max_species), x(gpu_max_species), species(gpu_max_species), conductivity(gpu_max_species)
      ! No vibrational inversion or diffusion is needed for these two outputs.
      call realgas_thermo(q, rho, y, ev, etr, rmix, cv)
      call realgas_mole_fractions(y, x)
      call blottner_mu_species(realgas_ttr(etr, cv), species)
      call wilke_mixture_viscosity(x, species, mu)
      call eucken_ktr_species(species, conductivity)
      call mason_saxena_ktr_mixture(x, conductivity, ktr)
   end subroutine realgas_mu_ktr

   subroutine realgas_transport(q, mu, ktr, kv, diff, htr, hv, gm, tv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: mu, ktr, kv, gm, tv
      real, intent(out) :: diff(nof_species), htr(nof_species), hv(nof_species)
      real :: rho, ev, etr, rmix, cv, ttr, pressure
      real :: y(gpu_max_species), x(gpu_max_species), species(gpu_max_species), conductivity(gpu_max_species)
      real :: deff(gpu_max_species), ht(gpu_max_species), hvib(gpu_max_species)
      call realgas_thermo(q, rho, y, ev, etr, rmix, cv)
      ttr = realgas_ttr(etr, cv)
      tv = realgas_vibrational_temperature(y, ev, ttr)
      gm = (cv + rmix) / cv
      pressure = 0.0e0
      if (cv > 1.0e-20) pressure = rho * rmix * (etr / cv)
      call realgas_mole_fractions(y, x)
      call effective_diffusion_direct(ttr, pressure, x, deff)
      call blottner_mu_species(ttr, species)
      call wilke_mixture_viscosity(x, species, mu)
      call eucken_ktr_species(species, conductivity)
      call mason_saxena_ktr_mixture(x, conductivity, ktr)
      call vibrational_conductivity(rho, y, deff, tv, kv)
      call htr_air5(ttr, ht)
      call hvib_air5(tv, hvib)
      diff = deff(1:nof_species); htr = ht(1:nof_species); hv = hvib(1:nof_species)
   end subroutine realgas_transport

   subroutine effective_diffusion_direct(temp, pressure, x, diff)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: temp, pressure, x(gpu_max_species)
      real, intent(out) :: diff(gpu_max_species)
      real :: xn(gpu_max_species), sx, denom, dij, sig, eps, omega, mass, patm
      integer :: i, j
      xn(1:nof_species) = max(x(1:nof_species), 0.0e0)
      sx = sum(xn(1:nof_species))
      if (sx > 1.0e-20) then
         xn(1:nof_species) = xn(1:nof_species) / sx
      else
         xn(1:nof_species) = 1.0e0 / max(1, nof_species)
      end if
      diff(1:nof_species) = 0.0e0
      patm = max(pressure / rgs_pa_per_atm, 1.0e-12)
      ! Evaluate each pair once; accumulate reciprocal diffusion sums in O(Ns) storage.
      do i = 1, nof_species
         do j = i + 1, nof_species
            sig = 0.5e0 * (rgs_sigmaa(i) + rgs_sigmaa(j))
            eps = sqrt(rgs_eps_over_k(i) * rgs_eps_over_k(j))
            omega = omega11_neufeld(max(temp / eps, 1.0e-6))
            mass = sqrt(1.0e0 / rgs_mg(i) + 1.0e0 / rgs_mg(j))
            denom = max(patm * sig**2 * omega, 1.0e-30)
            dij = (0.0018580e0 * temp * sqrt(temp) * mass) / denom * rgs_cm2s_to_m2s
            diff(i) = diff(i) + xn(j) / max(dij, 1.0e-20)
            diff(j) = diff(j) + xn(i) / max(dij, 1.0e-20)
         end do
      end do
      do i = 1, nof_species
         if (1.0e0 - xn(i) <= 1.0e-20 .or. diff(i) <= 1.0e-20) then
            diff(i) = 1.0e-10
         else
            diff(i) = max((1.0e0 - xn(i)) / diff(i), 1.0e-10)
         end if
      end do
   end subroutine effective_diffusion_direct

   subroutine realgas_log_equilibrium(t, logkc)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: t
      real, intent(out) :: logkc(5)
      real, parameter :: tref = 298.15e0, pref = 101325.0e0
! Standard entropies S/R at 298.15 K, 1 atm, in N2,O2,NO,N,O order.
! NASA Glenn TP-2002-211556, as distributed in Cantera airNASA9.yaml.
! Only the entropy integration constants are taken from NASA9. Enthalpy and
! heat capacity below follow OUR constant-cv + harmonic-oscillator EOS.
      real, parameter :: sref(5) = [23.0452199314975e0, 24.6736689675416e0, &
         & 25.3470090449906e0, 18.4378968249228e0, 19.3710101485200e0]
      real :: g(5), cp_r, theta, x, b, sv, svref, ev_rt, s_r
      integer :: k
      do k = 1, 5
         cp_r = 2.5e0
         if (k .le. 3) cp_r = 3.5e0
         theta = rg_thetag(k)
         sv = 0.0e0; svref = 0.0e0; ev_rt = 0.0e0
         if (theta .gt. 0.0e0) then
            x = theta / t; b = exp(-x)
            ev_rt = x * b / (1.0e0 - b)
            sv = ev_rt - log(1.0e0 - b)
            x = theta / tref; b = exp(-x)
            svref = x * b / (1.0e0 - b) - log(1.0e0 - b)
         end if
         s_r = sref(k) + cp_r * log(t / tref) + sv - svref
         g(k) = max(rg_hzero(k), 0.0e0) / (rgs_ru * t) + cp_r + ev_rt - s_r
      end do
! Kc uses mol/m^3. Dissociation has delta_nu=1; exchange has delta_nu=0.
      logkc(1) = g(1) - 2.0e0 * g(4) + log(pref / (rgs_ru * t))
      logkc(2) = g(2) - 2.0e0 * g(5) + log(pref / (rgs_ru * t))
      logkc(3) = g(3) - g(4) - g(5) + log(pref / (rgs_ru * t))
! Derive dependent reactions exactly: no inconsistent independently rounded fits.
      logkc(4) = logkc(1) - logkc(3)
      logkc(5) = logkc(3) - logkc(2)
   end subroutine realgas_log_equilibrium

   subroutine realgas_euler_jacobian(n, q, normal, amat)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: q(nof_variables), normal(3)
      real, intent(out) :: amat(gpu_max_nvar, gpu_max_nvar)
      real :: prim(gpu_max_nvar), dp(gpu_max_nvar), dvn(gpu_max_nvar)
      real :: rho, vel(3), vn, beta, rmix, cvmix, cv_s, rs_s, t, p, gamma_mix, pinf, h0
      integer :: k, j, idxe, idxev, idx
      amat = 0.0e0; dp = 0.0e0; dvn = 0.0e0; vel = 0.0e0
      idxe = dimensiona + 2; idxev = idxe + 1
      rho = q(1)
      prim = 0.0e0; prim(1:nof_variables) = q
      call cons2prim(n, prim, pinf, gamma_mix)
      p = prim(idxe); rmix = 0.0e0; cvmix = 0.0e0
      do k = 1, nof_species
         rs_s = rgs_ru / rg_molm(k)
         cv_s = 1.5e0 * rs_s
         if (k .le. 3) cv_s = 2.5e0 * rs_s
         rmix = rmix + q(idxev + k) * rs_s / rho
         cvmix = cvmix + q(idxev + k) * cv_s / rho
      end do
      beta = rmix / cvmix; t = p / (rho * rmix)
      vel(1:dimensiona) = q(2:dimensiona + 1) / rho
      vn = sum(vel(1:dimensiona) * normal(1:dimensiona))
      dp(1) = 0.5e0 * beta * sum(vel * vel)
      dp(2:dimensiona + 1) = -beta * vel(1:dimensiona)
      dp(idxe) = beta; dp(idxev) = -beta
      do k = 1, nof_species
         rs_s = rgs_ru / rg_molm(k)
         cv_s = 1.5e0 * rs_s
         if (k .le. 3) cv_s = 2.5e0 * rs_s
         h0 = max(rg_hzero(k), 0.0e0) / rg_molm(k)
         dp(idxev + k) = (rs_s - beta * cv_s) * t - beta * h0
      end do
      dvn(1) = -vn / rho
      dvn(2:dimensiona + 1) = normal(1:dimensiona) / rho
      amat(1, 2:dimensiona + 1) = normal(1:dimensiona)
      do k = 1, dimensiona
         amat(k + 1, 1:nof_variables) = q(k + 1) * dvn(1:nof_variables) + normal(k) * dp(1:nof_variables)
         amat(k + 1, k + 1) = amat(k + 1, k + 1) + vn
      end do
      amat(idxe, 1:nof_variables) = (q(idxe) + p) * dvn(1:nof_variables) + vn * dp(1:nof_variables)
      amat(idxe, idxe) = amat(idxe, idxe) + vn
      do idx = idxev, nof_variables
         amat(idx, 1:nof_variables) = q(idx) * dvn(1:nof_variables)
         amat(idx, idx) = amat(idx, idx) + vn
      end do
   end subroutine realgas_euler_jacobian

   real function realgas_diffusion_radius(q, mu_t, ktr_t) result(radius)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables), mu_t, ktr_t
      real :: mu, ktr, kv, gamma_mix, tv, rho, cvtr, cvv, rs_s, yi
      real :: diff(gpu_max_species), htr(gpu_max_species), hv(gpu_max_species)
      integer :: k
      call realgas_transport(q, mu, ktr, kv, diff, htr, hv, gamma_mix, tv)
      rho = max(q(1), tolsmall); cvtr = 0.0e0; cvv = 0.0e0
      do k = 1, nof_species
         yi = max(q(dimensiona + 3 + k) / rho, 0.0e0)
         rs_s = rgs_ru / rg_molm(k)
         if (k .le. 3) then
            cvtr = cvtr + yi * 2.5e0 * rs_s
         else
            cvtr = cvtr + yi * 1.5e0 * rs_s
         end if
         if (rg_thetag(k) .gt. 0.0e0) then
            cvv = cvv + yi * rgs_cv_vibrational_diatomic(tv, rg_thetag(k)) / rg_molm(k)
         end if
      end do
      radius = max((4.0e0 / 3.0e0) * (mu + max(mu_t, 0.0e0)) / rho, &
         & (ktr + max(ktr_t, 0.0e0)) / (rho * max(cvtr, tolsmall)))
      if (cvv .gt. 1.0e-30) radius = max(radius, kv / (rho * cvv))
      if (nof_species .gt. 0) radius = max(radius, maxval(diff(1:nof_species)))
   end function realgas_diffusion_radius

   real function realgas_face_diffusion_speed(ql, qr, distance, viscl, laml) result(speed)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: ql(nof_variables), qr(nof_variables), distance, viscl(4), laml(4)
      real :: mut_l, mut_r, kt_l, kt_r
      mut_l = 0.0e0; mut_r = 0.0e0; kt_l = 0.0e0; kt_r = 0.0e0
      if (turbulence .eq. 1) then
         mut_l = max(viscl(3), 0.0e0); mut_r = max(viscl(4), 0.0e0)
         ! eddyvisco returns TOTAL conductivity in slots 3:4, not an increment.
         kt_l = max(laml(3) - laml(1), 0.0e0); kt_r = max(laml(4) - laml(2), 0.0e0)
      end if
! Units m/s. Face measure is applied ONCE by the residual/matrix assembly.
      speed = max(realgas_diffusion_radius(ql, mut_l, kt_l), realgas_diffusion_radius(qr, mut_r, kt_r)) &
         & / max(abs(distance), tolsmall)
   end function realgas_face_diffusion_speed
   subroutine mrfswitch(n, iconsidered, facex, pointx, pox, poy)
      implicit none

      !> @brief
!> this subroutine  check if the element is on the rotational/stationary reference frame and update the mrf_origin and srf_velocity srf accordingly

#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered, facex, pointx
      real, dimension(3) :: mrf_origin, mrf_velocity, rotvel
      integer :: rotframe_on
      integer :: ninv, k
      real, dimension(1:dimensiona), intent(inout) :: pox, poy
!internal variables
      !element coordinates, roation_axys, cylinder_center_coordinates, vector_element_center, rotational velocity at gaussian points, gausian points coordinates
      real, dimension(3) :: p1p2, pc, popc, po, pgp, poz
      real :: d1, d2, r1, theta, dpopc, tempdp

      po = (pox(1:3))
      pgp = (poy(1:3))
!body
      do ninv = 1, nrotors

         pc(1:3) = (point1_gl(ninv, 1:3) + point2_gl(ninv, 1:3)) / 2 !center of cylinder
         p1p2(1:3) = point2_gl(ninv, 1:3) - point1_gl(ninv, 1:3) !axysvector
         popc(1:3) = po(1:3) - pc(1:3) ! vector elelement-centre
         dpopc = sqrt((po(1) - pc(1))**2 + (po(2) - pc(2))**2 + (po(3) - pc(3))**2) !distance between element and center

         tempdp = 0.00

         do k = 1, 3
            tempdp = tempdp + popc(k) * p1p2(k)
         end do

         !angle between element vector and axys
         theta = acos((tempdp) / (sqrt(popc(1)**2 + popc(2)**2 + popc(3)**2) * sqrt(p1p2(1)**2 + p1p2(2)**2 + p1p2(3)**2)))
         d2 = dpopc * abs(cos(theta))
         r1 = dpopc * abs(sin(theta))
         d1 = sqrt((point1_gl(ninv, 1) - pc(1))**2 + (point1_gl(ninv, 2) - pc(2))**2 + (point1_gl(ninv, 3) - pc(3))**2)

         if ((d1 .ge. d2) .and. (r1 .le. radius_gl(ninv))) then
            rotframe_on = 1
            mrf_origin(1:3) = pc(1:3)
            pox(1:3) = pgp(1:3) - mrf_origin(1:3)
            mrf_velocity(1:3) = mrf_rot_gl(ninv) * (p1p2) / sqrt(p1p2(1)**2 + p1p2(2)**2 + p1p2(3)**2)
!     srf_velocity(1)=0.0
!     srf_velocity(2)=mrf_rot_gl
!     srf_velocity(3)=0.0
            poy(1:3) = mrf_velocity(1:3)
            call vect_function(pox, poy, rotvel)

            go to 606
         else
            mrf_origin(1:3) = 0.0
            rotframe_on = 0
            mrf_velocity(1:3) = 0.0
            rotvel(1:3) = 0.0
         end if

      end do

606   continue

      rec_mrf_origin(iconsidered, 1:3) = mrf_origin(1:3)
      rec_mrf_velocity(iconsidered, 1:3) = mrf_velocity(1:3)
      rec_rotvel(iconsidered, facex, pointx, 1:3) = rotvel(1:3)
      rec_mrf(iconsidered) = rotframe_on

   end subroutine mrfswitch

   subroutine gqp_sol(n, leftv, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: leftv(1:nof_variables)
      real :: coord_gqp(1:2), rd
      real :: pi
      integer, intent(in) :: iconsidered, n
      integer :: i, j, k, facex, ngp

      pi = 4.00 * atan(1.00)
      i = iconsidered
      facex = 1
      ngp = 1
      coord_gqp(1:2) = rec_qpoints(i, facex, ngp, 1:2)

      leftv(1) = 0.00
      if (sqrt(((coord_gqp(1) - 0.250)**2) + ((coord_gqp(2) - 0.50)**2)) .le. 0.15) then
         rd = (1.00 / 0.150) * sqrt(((coord_gqp(1) - 0.250)**2) + ((coord_gqp(2) - 0.50)**2))

         leftv(1) = 0.250 * (1.00 + cos(pi * min(rd, 1.00)))
      end if

      if (sqrt(((coord_gqp(1) - 0.50)**2) + ((coord_gqp(2) - 0.250)**2)) .le. 0.15) then

         rd = (1.00 / 0.150) * sqrt(((coord_gqp(1) - 0.50)**2) + ((coord_gqp(2) - 0.250)**2))
         leftv(1) = 1.00 - rd
      end if

      if (sqrt(((coord_gqp(1) - 0.50)**2) + ((coord_gqp(2) - 0.750)**2)) .le. 0.15) then

         rd = (1.00 / 0.150) * sqrt(((coord_gqp(1) - 0.50)**2) + ((coord_gqp(2) - 0.750)**2))
         if ((abs(coord_gqp(1) - 0.5) .ge. 0.0250) .or. (coord_gqp(2) .gt. 0.85)) then

            leftv(1) = 1.00
         else

            leftv(1) = 0.00

         end if
      end if

   end subroutine gqp_sol

   subroutine fix_conservative_state(leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: leftv(1:nof_variables)

      integer :: k, idxe, idxev
      real :: rho, rho_floor, tiny
      real :: rhoe, rhoev, dev
      real :: sumrhoy, scale
      real, dimension(1:gpu_max_species) :: y
      real :: u, v, w, ke
      real :: tmin, tv_max
      real :: echem, cv_mix, etr_min, rhoe_min
      real :: rhoev_cap, sum_viby
!
!   ------------------------
!   tunable constants
!   ------------------------
      rho_floor = 1.0e-14
      tiny = 1.0e-30
      tmin = 200.00 ! translational temperature floor
      tv_max = 20000.00 ! cap to prevent tv blow-ups (try 12000-20000)

      idxe = dimensiona + 2
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
      if (rhoev < 0.00) then
         rhoev = 0.00
         leftv(idxev) = rhoev
      end if

      ! ------------------------
      ! 3) species positivity + normalization: enforce sum rhoy = rho
      ! ------------------------
      sumrhoy = 0.00
      do k = 1, nof_species
         if (leftv(idxev + k) < 0.00) leftv(idxev + k) = 0.00
         sumrhoy = sumrhoy + leftv(idxev + k)
      end do

      if (sumrhoy > 0.00) then
         scale = rho / sumrhoy
         do k = 1, nof_species
            leftv(idxev + k) = leftv(idxev + k) * scale
            y(k) = leftv(idxev + k) / rho
         end do
      else
         ! Pathological: restore the configured freestream composition.
         do k = 1, nof_species
            leftv(idxev + k) = rho * rg_vf(k)
            y(k) = rg_vf(k)
         end do
      end if

      ! ------------------------
      ! 4) cap vib energy to avoid tv spikes (composition-compatible)
      !     if rhoev too large for available vib species, reduce it and
      !     add the removed energy back into rhoe (conserve total).
      ! ------------------------
      rhoe = leftv(idxe)
      rhoev = leftv(idxev)

      sum_viby = 0.00
      do k = 1, 3
         sum_viby = sum_viby + y(k)
      end do

      if (sum_viby < 1.0e-12) then
         ! no vib-capable mass: force rhoev -> 0.  rhoe is total energy, so
         ! leaving it unchanged conserves total energy and returns this amount
         ! to the recovered translational mode.
         dev = rhoev
         rhoev = 0.00
      else
         ! cap corresponding to tv_max using your same oscillator formula
         rhoev_cap = 0.00
         do k = 1, 3
            if (y(k) > 0.00) then
               rhoev_cap = rhoev_cap + rho * y(k) * (rgs_ru / rg_molm(k)) * &
                  & (rg_thetag(k) / (exp(rg_thetag(k) / tv_max) - 1.00))
            end if
         end do

         if (rhoev > rhoev_cap) then
            dev = rhoev - rhoev_cap
            rhoev = rhoev_cap
            ! rhoe is total energy; do not add dev here.  Reducing rhoev while
            ! keeping rhoe fixed moves the excess into recovered translational energy.
         end if
      end if

      leftv(idxev) = rhoev
      leftv(idxe) = rhoe

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
         w = 0.00
      end if
      ke = 0.50 * (u * u + v * v + w * w)

      ! chemical formation energy per mass. Keep this sign consistent with
      ! prim2cons/cons2prim: atoms carry positive formation enthalpy.
      echem = 0.00
      do k = 1, nof_species
         if (rg_hzero(k) > 0.00) then
            echem = echem + y(k) * (rg_hzero(k) / rg_molm(k))
         end if
      end do

      ! cv_mix per mass
      cv_mix = 0.00
      do k = 1, nof_species
         if (k <= 3) then
            cv_mix = cv_mix + y(k) * (5.00 / 2.00) * (rgs_ru / rg_molm(k))
         else
            cv_mix = cv_mix + y(k) * (3.00 / 2.00) * (rgs_ru / rg_molm(k))
         end if
      end do
      if (cv_mix < tiny) cv_mix = tiny

      etr_min = cv_mix * tmin
      rhoe_min = rho * (etr_min + echem + ke) + rhoev

      if (leftv(idxe) < rhoe_min) leftv(idxe) = rhoe_min

   end subroutine fix_conservative_state

   subroutine multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, mp_d_eff, mp_htr, mp_hvib, gammal)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: mp_mu_mix, mp_ktr_mix, mp_kve, gammal
      real, intent(out) :: mp_d_eff(nof_species), mp_htr(nof_species), mp_hvib(nof_species)
      real :: tv
      call realgas_transport(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, mp_d_eff, mp_htr, mp_hvib, gammal, tv)
   end subroutine multispecies_mixtures_rg

   subroutine compute_real_gas_diffusion(rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho, rgs_x, rgs_y, &
      & rgs_dij, rgs_deff, rgs_mu_i, rgs_mu_mix, &
      & rgs_ktr_i, rgs_ktr_mix, rgs_kve, &
      & rgs_htr_i, rgs_hvib_i)
      implicit none

#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho
      real, intent(in) :: rgs_x(1:gpu_max_species), rgs_y(1:gpu_max_species)
      real, intent(out) :: rgs_dij(1:gpu_max_species, 1:gpu_max_species), rgs_deff(1:gpu_max_species)
      real, intent(out) :: rgs_mu_i(1:gpu_max_species), rgs_mu_mix
      real, intent(out) :: rgs_ktr_i(1:gpu_max_species), rgs_ktr_mix, rgs_kve
      real, intent(out) :: rgs_htr_i(1:gpu_max_species), rgs_hvib_i(1:gpu_max_species)

      ! ---- binary diffusion ----
      call compute_binary_diffusion(rgs_ttr, rgs_p_pa, rgs_dij)

      ! ---- mixture-average diffusion coefficient ----
      call compute_effective_diffusion(rgs_x, rgs_dij, rgs_deff)

      ! ---- species viscosity (blottner) ----
      call blottner_mu_species(rgs_ttr, rgs_mu_i)

      ! ---- mixture viscosity (wilke) ----
      call wilke_mixture_viscosity(rgs_x, rgs_mu_i, rgs_mu_mix)

      ! ---- species thermal conductivity (eucken) ----
      call eucken_ktr_species(rgs_mu_i, rgs_ktr_i)

      ! ---- mixture translational conductivity ----
      call mason_saxena_ktr_mixture(rgs_x, rgs_ktr_i, rgs_ktr_mix)

      ! ---- vibrational conductivity ----
      call vibrational_conductivity(rgs_rho, rgs_y, rgs_deff, rgs_tv, rgs_kve)

      ! ---- species translational and vibrational enthalpy ----
      call htr_air5(rgs_ttr, rgs_htr_i)
      call hvib_air5(rgs_tv, rgs_hvib_i)

   end subroutine compute_real_gas_diffusion

   function omega11_neufeld(rgs_tstar) result(rgs_omega11)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_tstar
      real :: rgs_omega11
      real, parameter :: rgs_a = 1.060360, rgs_b = 0.156100, rgs_c = 0.193000, rgs_d = 0.476350, &
         & rgs_e = 1.035870, rgs_f = 1.529960, rgs_g = 1.764740, rgs_h = 3.894110
      real :: tstar_eff

      tstar_eff = max(rgs_tstar, 1.0e-6) ! avoid 0^(-b) and tiny t*

      rgs_omega11 = rgs_a / exp(rgs_b * log(tstar_eff)) &
         & + rgs_c * exp(-rgs_d * tstar_eff) &
         & + rgs_e * exp(-rgs_f * tstar_eff) &
         & + rgs_g * exp(-rgs_h * tstar_eff)

   end function omega11_neufeld

   subroutine compute_binary_diffusion(rgs_t, rgs_p_pa, rgs_dij)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_t, rgs_p_pa
      real, intent(out) :: rgs_dij(1:gpu_max_species, 1:gpu_max_species)

      integer :: i, j
      real :: sig_ij, eps_ij, tstar, omega, p_atm, denom, mass_factor
      real, parameter :: tiny = 1e-30

      rgs_dij = 0.00

      ! pressure in atm
      p_atm = max(rgs_p_pa / rgs_pa_per_atm, 1e-12)

      do i = 1, nof_species
         do j = i + 1, nof_species

            sig_ij = 0.50 * (rgs_sigmaa(i) + rgs_sigmaa(j))
            eps_ij = sqrt(rgs_eps_over_k(i) * rgs_eps_over_k(j))

            ! non-dimensional temperature
            tstar = max(rgs_t / eps_ij, 1e-6)
            omega = omega11_neufeld(tstar)

            mass_factor = sqrt(1.0 / rgs_mg(i) + 1.0 / rgs_mg(j))
            denom = p_atm * sig_ij**2 * omega
            denom = max(denom, tiny)

            rgs_dij(i, j) = (0.0018580 * rgs_t * sqrt(rgs_t) * mass_factor) / denom * rgs_cm2s_to_m2s
            rgs_dij(j, i) = rgs_dij(i, j)

         end do
      end do

   end subroutine compute_binary_diffusion

   subroutine compute_effective_diffusion(rgs_x, rgs_dij, rgs_deff)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_x(1:gpu_max_species)
      real, intent(in) :: rgs_dij(1:gpu_max_species, 1:gpu_max_species)
      real, intent(out) :: rgs_deff(1:gpu_max_species)

      integer :: i, j
      real :: sumj, one_minus_xi, sumx
      real :: xloc(1:gpu_max_species)
      real, parameter :: tiny = 1e-20
      real, parameter :: dmin = 1e-10 ! lower bound on diffusivity [m^2/s]

      ! work on a local copy to enforce positivity and normalisation

      xloc = 0.00
      xloc(1:nof_species) = rgs_x(1:nof_species)

! enforce positivity
      do i = 1, nof_species
         if (xloc(i) < 0.0) xloc(i) = 0.0
      end do

! enforce normalisation
      sumx = sum(xloc(1:nof_species))
      if (sumx > tiny) then
         xloc(1:nof_species) = xloc(1:nof_species) / sumx
      else
         xloc(1:nof_species) = 1.0 / max(1, nof_species)
      end if

      do i = 1, nof_species
         sumj = 0.00

         do j = 1, nof_species
            if (j /= i) then
               sumj = sumj + xloc(j) / max(rgs_dij(i, j), tiny)
            end if
         end do

         one_minus_xi = max(1.00 - xloc(i), 0.00)

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
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_t
      real, intent(out) :: rgs_mu(1:gpu_max_species)
      integer :: i
      real :: tlog

      tlog = log10(max(rgs_t, 50.0))

      do i = 1, nof_species
         rgs_mu(i) = 1.e-7 * exp(2.3025850929940456840 * (rgs_ab(i) * tlog * tlog + rgs_bb(i) * tlog + rgs_cb(i)))
         rgs_mu(i) = max(rgs_mu(i), 1e-12)
      end do
   end subroutine blottner_mu_species

   subroutine wilke_mixture_viscosity(rgs_x, rgs_mu_i, rgs_mu_mix)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_x(1:gpu_max_species), rgs_mu_i(1:gpu_max_species)
      real, intent(out) :: rgs_mu_mix
      integer :: i, j
      real :: phi_ij, denom

      rgs_mu_mix = 0.0

      do i = 1, nof_species
         denom = 0.0

         do j = 1, nof_species
            if (i == j) then
               denom = denom + rgs_x(j)
            else
               phi_ij = (1.0 + sqrt(rgs_mu_i(i) / rgs_mu_i(j)) * sqrt(sqrt(rg_molm(j) / rg_molm(i))))**2 &
                  & / (sqrt(8.0) * sqrt(1.0 + rg_molm(i) / rg_molm(j)))

               denom = denom + rgs_x(j) * phi_ij
            end if
         end do

         rgs_mu_mix = rgs_mu_mix + rgs_x(i) * rgs_mu_i(i) / max(denom, 1e-20)
      end do

   end subroutine wilke_mixture_viscosity

   subroutine cp_tr_species(rgs_cp_tr)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(out) :: rgs_cp_tr(1:gpu_max_species) ! j/kg-k
      integer :: rgs_i
      real :: rgs_rspec
      do rgs_i = 1, nof_species
         rgs_rspec = rgs_ru / rg_molm(rgs_i)
         if (rgs_i <= 3) then
            rgs_cp_tr(rgs_i) = 3.50 * rgs_rspec ! diatomics: 7/2 r
         else
            rgs_cp_tr(rgs_i) = 2.50 * rgs_rspec ! atoms:     5/2 r
         end if
      end do
   end subroutine cp_tr_species

   subroutine eucken_ktr_species(rgs_mu_i, rgs_ktr_i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_mu_i(1:gpu_max_species)
      real, intent(out) :: rgs_ktr_i(1:gpu_max_species)

      integer :: i
      real :: rspec, cp

      do i = 1, nof_species
         rspec = rgs_ru / rg_molm(i)

         if (i <= 3) then
            cp = 3.50 * rspec ! diatomic
            rgs_ktr_i(i) = rgs_mu_i(i) * (cp + 1.250 * rspec)
         else
            cp = 2.50 * rspec ! atomic
            rgs_ktr_i(i) = rgs_mu_i(i) * (cp + 1.500 * rspec)
         end if

      end do

   end subroutine eucken_ktr_species

   subroutine mason_saxena_ktr_mixture(rgs_x, rgs_ktr_i, rgs_ktr_mix)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_x(1:gpu_max_species), rgs_ktr_i(1:gpu_max_species)
      real, intent(out) :: rgs_ktr_mix

      integer :: i, j
      real :: psi_ij, denom

      rgs_ktr_mix = 0.0

      do i = 1, nof_species
         denom = 0.0

         do j = 1, nof_species
            if (i == j) then
               denom = denom + rgs_x(j)
            else
               psi_ij = (1.0 + sqrt(rgs_ktr_i(i) / rgs_ktr_i(j)) * sqrt(sqrt(rg_molm(j) / rg_molm(i))))**2 &
                  & / (sqrt(8.0) * sqrt(1.0 + rg_molm(i) / rg_molm(j)))

               denom = denom + rgs_x(j) * psi_ij
            end if
         end do

         rgs_ktr_mix = rgs_ktr_mix + rgs_x(i) * rgs_ktr_i(i) / max(denom, 1e-20)

      end do

   end subroutine mason_saxena_ktr_mixture

   function rgs_cv_vibrational_diatomic(rgs_tv, rgs_theta) result(rgs_cv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_tv, rgs_theta
      real :: rgs_cv, rgs_x, ex

      if (rgs_tv <= 1.00 .or. rgs_theta <= 0.00) then
         rgs_cv = 0.00
      else
         rgs_x = rgs_theta / rgs_tv
         if (rgs_x > 60.00) then
            ! asymptotic form for large x: cv ∝ x^2 e^(-x)
            rgs_cv = rgs_ru * (rgs_x * rgs_x) * exp(-rgs_x)
         else
            ex = exp(rgs_x)
            rgs_cv = rgs_ru * (rgs_x * rgs_x) * ex / ((ex - 1.00)**2)
         end if
      end if
   end function rgs_cv_vibrational_diatomic

   function rgs_hvib_species(rgs_tv, rgs_theta, rgs_mi) result(rgs_hv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_tv, rgs_theta, rgs_mi
      real :: rgs_hv, rgs_x, rgs_ri, ex

      if (rgs_theta <= 0.00 .or. rgs_tv <= 1.00) then
         rgs_hv = 0.00
         return
      end if

      rgs_ri = rgs_ru / rgs_mi
      rgs_x = rgs_theta / rgs_tv

      if (rgs_x > 60.00) then
         ! asymptotic form: h_v ≈ r_i θ e^(-x)
         rgs_hv = rgs_ri * rgs_theta * exp(-rgs_x)
      else
         ex = exp(rgs_x)
         rgs_hv = rgs_ri * rgs_theta / (ex - 1.00)
      end if

   end function rgs_hvib_species

   subroutine hvib_air5(rgs_tv, rgs_hvib_i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_tv
      real, intent(out) :: rgs_hvib_i(1:gpu_max_species) ! j/kg
      integer :: rgs_i

      do rgs_i = 1, nof_species
         rgs_hvib_i(rgs_i) = rgs_hvib_species(rgs_tv, rg_thetag(rgs_i), rg_molm(rgs_i))
      end do

   end subroutine hvib_air5

   subroutine htr_air5(rgs_ttr, rgs_htr_i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rgs_ttr
      real, intent(out) :: rgs_htr_i(1:gpu_max_species)
      real :: cp(1:gpu_max_species), h0_mass(1:gpu_max_species)
      integer :: i

      call cp_tr_species(cp)

      do i = 1, nof_species
         h0_mass(i) = rg_hzero(i) / rg_molm(i)
         rgs_htr_i(i) = h0_mass(i) + cp(i) * rgs_ttr
      end do

   end subroutine htr_air5

   subroutine vibrational_conductivity(rgs_rho, rgs_y, rgs_deff, rgs_tv, rgs_kve)
      implicit none
#ifdef xpu
!$omp declare target
#endif

      real, intent(in) :: rgs_rho, rgs_y(1:gpu_max_species), rgs_deff(1:gpu_max_species), rgs_tv
      real, intent(out) :: rgs_kve

      real :: cv_vib_mass(1:gpu_max_species)
      integer :: i

      do i = 1, nof_species
         if (rg_thetag(i) > 0.00) then
            cv_vib_mass(i) = rgs_cv_vibrational_diatomic(rgs_tv, rg_thetag(i)) / rg_molm(i)
         else
            cv_vib_mass(i) = 0.00
         end if
      end do

      rgs_kve = 0.00
      do i = 1, nof_species
         rgs_kve = rgs_kve + rgs_rho * rgs_y(i) * rgs_deff(i) * cv_vib_mass(i)
      end do

      rgs_kve = max(rgs_kve, 0.00)

   end subroutine vibrational_conductivity

   subroutine multispecies_mixtures(leftv, mp_temp, mp_mu_mix, mp_k_mix, mp_cp_mix, gammal)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: mp_temp, mp_mu_mix, mp_k_mix, mp_cp_mix, gammal
      real :: y(gpu_max_species), x(gpu_max_species), mu(gpu_max_species)
      real :: rmix, cv, energy, total, phi, den, cp_i, h_i, k_i
      integer :: i, j
      call multispecies_thermo(leftv, mp_temp, y, rmix, mp_cp_mix, cv, energy)
      gammal = mp_cp_mix / cv
      total = 0.0e0
      do i = 1, nof_species
         x(i) = y(i) / mp_m(i)
         total = total + x(i)
         ! The supplied table contains Blottner log10 coefficients (Pa s).
         phi = log10(max(mp_temp, 1.0e0))
         mu(i) = max(1.0e-12, 1.0e-7 * exp(log(10.0e0) * &
            & (mp_brok_a(i) * phi * phi + mp_brok_b(i) * phi + mp_brok_c(i))))
      end do
      x(1:nof_species) = x(1:nof_species) / total
      mp_mu_mix = 0.0e0; mp_k_mix = 0.0e0
      do i = 1, nof_species
         den = 0.0e0
         do j = 1, nof_species
            phi = 1.0e0
            if (i /= j) phi = (1.0e0 + sqrt(mu(i) / mu(j)) * sqrt(sqrt(mp_m(j) / mp_m(i))))**2 &
               & / sqrt(8.0e0 * (1.0e0 + mp_m(i) / mp_m(j)))
            den = den + x(j) * phi
         end do
         call species_cp_h_sensible(mp_temp, i, cp_i, h_i)
         k_i = (cp_i + 1.25e0 * rgs_ru / mp_m(i)) * mu(i)
         mp_mu_mix = mp_mu_mix + x(i) * mu(i) / den
         mp_k_mix = mp_k_mix + x(i) * k_i / den
      end do
   end subroutine multispecies_mixtures

   subroutine multispecies_temp(leftv, mp_temp)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: mp_temp
      real :: y(gpu_max_species), rmix, cp, cv, energy
      call multispecies_thermo(leftv, mp_temp, y, rmix, cp, cv, energy)
   end subroutine multispecies_temp

   subroutine normalize_mf(mpc_mp)
      implicit none
#ifdef xpu
!$omp declare target
#endif
    !! ensure mole fractions sum to 1 (robust to roundoff)
      real, dimension(1:nof_species), intent(inout) :: mpc_mp
      real :: s
      integer :: k
      s = zero
      do k = 1, nof_species
         s = s + max(zero, mpc_mp(k))
      end do
      if (s > zero) then
         do k = 1, nof_species
            mpc_mp(k) = max(zero, mpc_mp(k)) / s
         end do
      else
         do k = 1, nof_species
            mpc_mp(k) = zero
         end do
      end if
   end subroutine normalize_mf

   subroutine molef_to_massf(mpc_mp, mpc_mass)
      implicit none
#ifdef xpu
!$omp declare target
#endif
    !! convert mole fractions to mass fractions (for mass-specific energy)
      real, dimension(1:nof_species), intent(in) :: mpc_mp
      real, dimension(1:nof_species), intent(out) :: mpc_mass
      real :: denom
      integer :: s
      denom = zero
      do s = 1, nof_species
         denom = denom + mpc_mp(s) * mp_m(s)
      end do
      if (denom <= zero) then
         mpc_mass = zero
      else
         do s = 1, nof_species
            mpc_mass(s) = mpc_mp(s) * mp_m(s) / denom
         end do
      end if
   end subroutine molef_to_massf

   subroutine species_cp_h_sensible(mp_temp, is, mp_cpl, mp_hl)
      implicit none
#ifdef xpu
!$omp declare target
#endif
    !! janaf 5-coeffs, sensible enthalpy with h(298.15k)=0 for each species
    !!   cp/r = a1 + a2*t + a3*t^2 + a4*t^3 + a5*t^4
    !!    h/r = a1*t + a2*t^2/2 + a3*t^3/3 + a4*t^4/4 + a5*t^5/5  (then shifted so h(298.15)=0)
      real, intent(in) :: mp_temp
      integer, intent(in) :: is
      real, intent(out) :: mp_cpl, mp_hl
      integer :: r
      real :: rs, tref, split
      rs = rgs_ru / mp_m(is)
      r = merge(1, 2, mp_temp <= mp_tmid_in(is))
      mp_cpl = rs * (mp_janaf(is, r, 1) + mp_temp * (mp_janaf(is, r, 2) + mp_temp * &
         & (mp_janaf(is, r, 3) + mp_temp * (mp_janaf(is, r, 4) + mp_temp * mp_janaf(is, r, 5)))))
      tref = 298.15e0; split = mp_tmid_in(is)
      if (mp_temp <= split .and. tref <= split) then
         mp_hl = rs * (janaf_cp_integral(mp_temp, is, 1) - janaf_cp_integral(tref, is, 1))
      elseif (mp_temp > split .and. tref > split) then
         mp_hl = rs * (janaf_cp_integral(mp_temp, is, 2) - janaf_cp_integral(tref, is, 2))
      elseif (mp_temp > split) then
         mp_hl = rs * (janaf_cp_integral(mp_temp, is, 2) - janaf_cp_integral(split, is, 2) &
            & + janaf_cp_integral(split, is, 1) - janaf_cp_integral(tref, is, 1))
      else
         mp_hl = rs * (janaf_cp_integral(mp_temp, is, 1) - janaf_cp_integral(split, is, 1) &
            & + janaf_cp_integral(split, is, 2) - janaf_cp_integral(tref, is, 2))
      end if
   end subroutine species_cp_h_sensible

   real function janaf_cp_integral(temp, species, range_id) result(value)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: temp
      integer, intent(in) :: species, range_id
      value = temp * (mp_janaf(species, range_id, 1) + temp * (mp_janaf(species, range_id, 2) / 2.0e0 &
         & + temp * (mp_janaf(species, range_id, 3) / 3.0e0 + temp * (mp_janaf(species, range_id, 4) / 4.0e0 &
         & + temp * mp_janaf(species, range_id, 5) / 5.0e0))))
   end function janaf_cp_integral

   subroutine mixture_props(mp_temp, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
      implicit none
#ifdef xpu
!$omp declare target
#endif
    !! mixture properties with mass-weighted energies; energies are sensible
      real, intent(in) :: mp_temp
      real, dimension(1:nof_species), intent(in) :: mp_mol_x
      real, intent(out) :: mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix
      real, dimension(1:gpu_max_species) :: mpc_mass
      real :: mp_cp_s, mp_h_s
      integer :: s
      call molef_to_massf(mp_mol_x, mpc_mass)
      mp_rmix = zero
      mp_cp_mix = zero
      mp_e_mix = zero
      do s = 1, nof_species
         if (mpc_mass(s) <= zero) cycle
         mp_rmix = mp_rmix + mpc_mass(s) * (rgs_ru / mp_m(s))
         call species_cp_h_sensible(mp_temp, s, mp_cp_s, mp_h_s)
         mp_e_mix = mp_e_mix + mpc_mass(s) * (mp_h_s - (rgs_ru / mp_m(s)) * mp_temp) ! e = h - rs*t (sensible e)
         mp_cp_mix = mp_cp_mix + mpc_mass(s) * mp_cp_s
      end do
      mp_cv_mix = mp_cp_mix - mp_rmix
   end subroutine mixture_props

   subroutine normalize_x(x)
      implicit none
#ifdef xpu
!$omp declare target
#endif
    !! ensure mole fractions sum to 1 (robust to roundoff)
      real, intent(inout) :: x(nof_species)
      real :: s
      integer :: k
      s = 0.00
      do k = 1, nof_species; s = s + max(0.00, x(k)); end do
      if (s > 0.00) then
         do k = 1, nof_species; x(k) = max(0.00, x(k)) / s; end do
      else
         do k = 1, nof_species; x(k) = 0.00; end do
      end if
   end subroutine normalize_x

   subroutine get_temp_janaf(leftv, mp_temp, mp_mass_fraction, mp_mole_fraction)
      implicit none
#ifdef xpu
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
      real, dimension(1:nof_variables), intent(in) :: leftv
      real, intent(out) :: mp_temp
      real, dimension(1:nof_species), intent(in) :: mp_mole_fraction, mp_mass_fraction
      integer :: mp_info, s
      real :: mp_tol, mp_sum, u, v, w
      integer :: mp_max_iter, i, j, k
      real, dimension(1:gpu_max_species) :: mp_mol_x
      real :: mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix
      real :: mp_kin, mp_e_tgt, mp_tlo, mp_thi, mp_fa, mp_fb, mp_f, mp_df, mp_at, mp_bt, mp_tn, mp_fn, mp_tguess, mp_atol
      integer :: mp_it, mp_itmax, mp_s
      logical :: mp_bracket

      ! copy/normalize x
      mp_mol_x = zero

      mp_mol_x(1:nof_species) = mp_mole_fraction(1:nof_species)

      call normalize_x(mp_mol_x)

      u = leftv(2) / leftv(1)
      v = leftv(3) / leftv(1)

      if (dimensiona .eq. 3) then
         w = leftv(4) / leftv(1)
      else
         w = zero
      end if

      mp_itmax = 50
      mp_kin = 0.50 * (u * u + v * v + w * w)
      mp_e_tgt = (leftv(dimensiona + 2) / leftv(1)) - mp_kin ! target sensible internal energy [j/kg]

      ! temperature bracket honoring species ranges
      mp_tlo = 150.00
      mp_thi = 20000.00
      do s = 1, nof_species
         if (mp_mol_x(s) <= zero) cycle
         mp_tlo = max(mp_tlo, mp_tlow_in(s))
         mp_thi = min(mp_thi, mp_thigh_in(s))
      end do
      if (mp_tlo >= mp_thi) then
         mp_info = -3; mp_temp = 300.00; return
      end if

      ! initial guess from cv around 300 k
      call mixture_props(300.00, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
      mp_tguess = 300.00 + (mp_e_tgt - mp_e_mix) / max(1.0e-8, mp_cv_mix)
      mp_temp = min(max(mp_tguess, mp_tlo), mp_thi)

      ! establish/verify bracket
      call mixture_props(mp_tlo, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix); mp_fa = mp_e_mix - mp_e_tgt
      call mixture_props(mp_thi, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix); mp_fb = mp_e_mix - mp_e_tgt
      mp_bracket = (mp_fa * mp_fb <= 0.0)
      mp_at = mp_tlo; mp_bt = mp_thi

      do mp_it = 1, mp_itmax
         call mixture_props(mp_temp, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
         mp_f = mp_e_mix - mp_e_tgt
         mp_df = max(mp_cv_mix, 1.0e-20)

         mp_atol = 1.0e-6 * max(1.00, abs(mp_e_tgt)) + 1.0e-3
         if (abs(mp_f) <= mp_atol) then
            mp_info = 0; return
         end if

         ! safeguarded newton step
         mp_tn = mp_temp - mp_f / mp_df
         if ((.not. mp_bracket) .or. (mp_tn <= mp_at) .or. (mp_tn >= mp_bt)) mp_tn = 0.50 * (mp_at + mp_bt)

         call mixture_props(mp_tn, mp_mol_x, mp_rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
         mp_fn = mp_e_mix - mp_e_tgt

         if (mp_fa * mp_fn <= 0.0) then
            mp_bt = mp_tn; mp_fb = mp_fn; mp_bracket = .true.
         else
            mp_at = mp_tn; mp_fa = mp_fn
         end if
         mp_temp = mp_tn

         if (abs(mp_bt - mp_at) <= 1.0e-8 * max(1.0, mp_temp)) then
            mp_info = 1; return
         end if
      end do

      mp_info = 2 ! iteration limit
   end subroutine get_temp_janaf

   function fluxeval2d(leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:gpu_max_nvar) :: fluxeval2d
      real, dimension(1:nof_variables), intent(in) :: leftv
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      p = leftv(4)
      gm = gamma
!kinetic energy first!
      skin = (oo2) * ((u**2) + (v**2))
!internal energy
      ien = ((p) / ((gm - 1.00) * r))
!total energy
      e = r * (skin + ien)

      fluxeval2d(1) = r * u
      fluxeval2d(2) = (r * (u**2)) + p
      fluxeval2d(3) = r * u * v
      fluxeval2d(4) = u * (e + p)

   end function

   function fluxeval3d(leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:gpu_max_nvar) :: fluxeval3d
      real, dimension(1:nof_variables), intent(in) :: leftv
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      w = leftv(4)
      p = leftv(5)
      gm = gamma
!kinetic energy first!
      skin = (oo2) * ((u**2) + (v**2) + (w**2))
!internal energy
      ien = ((p) / ((gm - 1.00) * r))
!total energy
      e = r * (skin + ien)

      fluxeval3d(1) = r * u
      fluxeval3d(2) = (r * (u**2)) + p
      fluxeval3d(3) = r * u * v
      fluxeval3d(4) = r * u * w
      fluxeval3d(5) = u * (e + p)

   end function

   subroutine cons2prim(n, q, pinf, gm)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(inout) :: q(nof_variables)
      real, intent(out) :: pinf, gm
      real :: rho, ev, etr, rmix, cv, y(gpu_max_species), temp, cp, energy, stiff, ke
      integer :: k
      pinf = 0.0e0; gm = gamma
      if (nof_variables <= 1) return
      if (realgas == 1) then
         call realgas_thermo(q, rho, y, ev, etr, rmix, cv)
         gm = (cv + rmix) / cv
         q(2:dimensiona + 1) = q(2:dimensiona + 1) / rho
         q(dimensiona + 2) = 0.0e0
         if (cv > 1.0e-20) q(dimensiona + 2) = rho * rmix * (etr / cv)
         q(dimensiona + 3) = ev
         q(dimensiona + 4:nof_variables) = y(1:nof_species)
      elseif (multispecies == 1) then
         if (mp_modelc == 0) then
            call allaire_properties(q, rho, gm, pinf, stiff)
            ke = 0.0e0
            do k = 2, dimensiona + 1
               q(k) = q(k) / rho
               ke = ke + q(k)**2
            end do
            q(1) = rho
            q(dimensiona + 2) = (gm - 1.0e0) * (q(dimensiona + 2) - 0.5e0 * rho * ke) - stiff
         else
            call multispecies_thermo(q, temp, y, rmix, cp, cv, energy)
            rho = q(1)
            gm = cp / cv
            q(2:dimensiona + 1) = q(2:dimensiona + 1) / rho
            q(dimensiona + 2) = rho * rmix * temp
            ! Preserve partial densities, as required by the multispecies fluxes.
         end if
      else
         rho = q(1)
         ke = 0.0e0
         do k = 2, dimensiona + 1
            q(k) = q(k) / rho
            ke = ke + q(k)**2
         end do
         q(dimensiona + 2) = (gamma - 1.0e0) * (q(dimensiona + 2) - oo2 * rho * ke)
      end if
   end subroutine cons2prim

   subroutine cons2div(n, q, pinf, gm)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(inout) :: q(nof_variables)
      real, intent(out) :: pinf, gm
      real :: temp, y(gpu_max_species), rmix, cp, cv, energy
      pinf = 0.0e0; gm = gamma
      if (nof_variables <= 1) return
      if (realgas == 1) then
         call cons2div_realgas(q, pinf, gm)
      elseif (multispecies == 1) then
         call multispecies_thermo(q, temp, y, rmix, cp, cv, energy)
         if (mp_modelc == 0) then
            call cons2prim(n, q, pinf, gm)
         else
            pinf = 0.0e0; gm = cp / cv
            q(2:dimensiona + 1) = q(2:dimensiona + 1) / q(1)
         end if
         q(dimensiona + 2) = temp
      else
         call cons2prim(n, q, pinf, gm)
         q(dimensiona + 2) = q(dimensiona + 2) / (q(1) * r_gas)
      end if
   end subroutine cons2div

   subroutine cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
      implicit none
!> @brief
!> this subroutine transforms two vector of conservative variables to primitive variables
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(inout) :: leftv
      real, intent(inout) :: mp_pinfl, gammal
      real, dimension(1:nof_variables), intent(inout) :: rightv
      real, intent(inout) :: mp_pinfr, gammar

      call cons2prim(n, leftv, mp_pinfl, gammal)
      call cons2prim(n, rightv, mp_pinfr, gammar)

   end subroutine cons2prim2

   subroutine lmacht(n, leftv, rightv)
      implicit none
!> @brief
!> this subroutine applies the low-mach number correction to two vectors of conserved variables
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(inout) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:nof_variables), intent(inout) :: rightv
      real :: mp_pinfr, gammar
      real :: q2l, q2r, uul, uur, vvl, vvr, wwr, wwl, rhol, rhor, etal, etar, duu, dvv, dww
      real :: mach2, mach, cma, dus, dvs, dws, diff, c1o2, ssl, ssr, ppl, ppr, eel, eer, tole, mlm

      tole = zero

      eel = leftv(5)
      eer = rightv(5)

      call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

      c1o2 = 0.500
      rhol = leftv(1)
      uul = leftv(2)
      vvl = leftv(3)
      wwl = leftv(4)
      ppl = leftv(5)
      rhor = rightv(1)
      uur = rightv(2)
      vvr = rightv(3)
      wwr = rightv(4)
      ppr = rightv(5)

      q2l = (uul * uul) + (vvl * vvl) + (wwl * wwl)
      q2r = (uur * uur) + (vvr * vvr) + (wwr * wwr)

!        ssl=((gamma*ppl)/(rhol)); ssr=((gamma*ppr)/(rhor))

!        mlm=((abs(ppr-ppl))/(0.5*rres*((gamma*pres/rres))))

      if ((multispecies .eq. 1) .or. (realgas .eq. 1)) then
         ssl = sqrt((leftv(5) + mp_pinfl) * gammal / leftv(1))
         ssr = sqrt((rightv(5) + mp_pinfr) * gammar / rightv(1))
      else

         ssl = ((gamma * ppl) / (rhol))
         ssr = ((gamma * ppr) / (rhor))
      end if

      cma = 1.00

      duu = uur - uul
      dvv = vvr - vvl
      dww = wwr - wwl

!       if(lmach.eq.1) then !standard proportional to du^2
      mach2 = max(q2l / ssl, q2r / ssr)
      mach = sqrt(mach2)
      mach = min(cma * mach, 1.00)
!       end if

      dus = uur + uul
      dvs = vvr + vvl
      dws = wwr + wwl

      duu = mach * duu
      dvv = mach * dvv
      dww = mach * dww

      diff = c1o2 * duu

      uul = (dus * c1o2 - diff)
      uur = (dus * c1o2 + diff)

      if (lmach_style .eq. 1) then
         diff = c1o2 * dvv

         vvl = (dvs * c1o2 - diff)
         vvr = (dvs * c1o2 + diff)

         diff = c1o2 * dww

         wwl = (dws * c1o2 - diff)
         wwr = (dws * c1o2 + diff)
      end if

      leftv(1) = rhol
      leftv(2) = uul * rhol
      leftv(3) = vvl * rhol
      leftv(4) = wwl * rhol
      leftv(5) = eel

      rightv(1) = rhor
      rightv(2) = uur * rhor
      rightv(3) = vvr * rhor
      rightv(4) = wwr * rhor
      rightv(5) = eer

   end subroutine lmacht

   subroutine lmacht2d(n, leftv, rightv)
      implicit none
!> @brief
!> this subroutine applies the low-mach number correction to two vectors of conserved variables 2d
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(inout) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:nof_variables), intent(inout) :: rightv
      real :: mp_pinfr, gammar
      real :: q2l, q2r, uul, uur, vvl, vvr, wwr, wwl, rhol, rhor, etal, etar, duu, dvv, dww
      real :: mach2, mach, cma, dus, dvs, dws, diff, c1o2, ssl, ssr, ppl, ppr, eel, eer, tole, mlm

      tole = tolsmall

      eel = leftv(4)
      eer = rightv(4)

      call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

      c1o2 = 0.50
      rhol = leftv(1)
      uul = leftv(2)
      vvl = leftv(3)
      ppl = leftv(4)

      rhor = rightv(1)
      uur = rightv(2)
      vvr = rightv(3)
      ppr = rightv(4)

      q2l = (uul * uul) + (vvl * vvl)
      q2r = (uur * uur) + (vvr * vvr)

      if ((multispecies .eq. 1) .or. (realgas .eq. 1)) then
         ssl = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1))
         ssr = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1))
      else

         ssl = ((gamma * ppl) / (rhol))
         ssr = ((gamma * ppr) / (rhor))
      end if

      cma = 1.00

      duu = uur - uul
      dvv = vvr - vvl

!       if(lmach.eq.1) then !standard proportional to du^2
      mach2 = max(q2l / ssl, q2r / ssr)
      mach = sqrt(mach2)
      mach = min(cma * mach, 1.00)
!       end if

      dus = uur + uul
      dvs = vvr + vvl

      !ul+zul+ur-zur=(ul+ur)-z(ur-ul))
      duu = mach * duu
      dvv = mach * dvv

      diff = c1o2 * duu
      !if (uul*uur.gt.tole)then

      uul = (dus * c1o2) - diff
      uur = (dus * c1o2) + diff

      if (lmach_style .eq. 1) then
         diff = c1o2 * dvv
         vvl = (dvs * c1o2 - diff)
         vvr = (dvs * c1o2 + diff)
      end if

      leftv(1) = rhol
      leftv(2) = uul * rhol
      leftv(3) = vvl * rhol
      leftv(4) = eel

      rightv(1) = rhor
      rightv(2) = uur * rhor
      rightv(3) = vvr * rhor
      rightv(4) = eer

   end subroutine lmacht2d

   subroutine prim2cons(n, q)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(inout) :: q(nof_variables)
      real :: rho, rmix, cv, echem, ke, rs, y, pressure, ev, gm, pinf, stiff, temp, cp, energy
      real :: x(gpu_max_species)
      integer :: k
      if (nof_variables <= 1) return
      if (realgas == 1) then
         rho = q(1); pressure = q(dimensiona + 2); ev = q(dimensiona + 3)
         rmix = 0.0e0; cv = 0.0e0; echem = 0.0e0
         do k = 1, nof_species
            y = q(dimensiona + 3 + k)
            rs = rgs_ru / rg_molm(k)
            rmix = rmix + y / rg_molm(k)
            cv = cv + y * merge(2.5e0, 1.5e0, k <= 3) * rs
            if (rg_hzero(k) > 0.0e0) echem = echem + y * rg_hzero(k) / rg_molm(k)
         end do
         rmix = rgs_ru * rmix
         ke = 0.0e0
         do k = 2, dimensiona + 1
            ke = ke + q(k)**2
            q(k) = rho * q(k)
         end do
         q(dimensiona + 2) = rho * (cv * (pressure / (rho * rmix)) + ev + echem + 0.5e0 * ke)
         q(dimensiona + 3:nof_variables) = rho * q(dimensiona + 3:nof_variables)
      elseif (multispecies == 1) then
         pressure = q(dimensiona + 2)
         if (mp_modelc == 0) then
            call allaire_properties(q, rho, gm, pinf, stiff)
            energy = (pressure + stiff) / ((gm - 1.0e0) * rho)
         else
            rho = q(1)
            rmix = 0.0e0
            do k = 1, nof_species
               x(k) = q(dimensiona + 2 + k) / (rho * mp_m(k))
               rmix = rmix + x(k)
            end do
            x(1:nof_species) = x(1:nof_species) / rmix
            temp = pressure / (rho * rgs_ru * rmix)
            call mixture_props(temp, x, rmix, cp, cv, energy)
         end if
         ke = 0.0e0
         do k = 2, dimensiona + 1
            ke = ke + q(k)**2
            q(k) = rho * q(k)
         end do
         q(1) = rho
         q(dimensiona + 2) = rho * (energy + 0.5e0 * ke)
      else
         rho = q(1)
         energy = q(dimensiona + 2) / ((gamma - 1.0e0) * rho)
         ke = 0.0e0
         do k = 2, dimensiona + 1
            ke = ke + q(k)**2
            q(k) = rho * q(k)
         end do
         q(dimensiona + 2) = rho * (energy + oo2 * ke)
      end if
   end subroutine prim2cons

   subroutine prim2cons2(n, leftv, rightv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(inout) :: leftv(nof_variables), rightv(nof_variables)
      call prim2cons(n, leftv)
      call prim2cons(n, rightv)
   end subroutine prim2cons2

   subroutine multispecies_reference_state(q)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(out) :: q(nof_variables)
      integer :: k
      q = 0.0e0
      q(2) = uvel; q(3) = vvel
      if (dimensiona == 3) q(4) = wvel
      q(dimensiona + 2) = pres
      do k = 1, nof_species
         q(dimensiona + 2 + k) = mp_r_in(k) * mp_a_in(k)
         q(1) = q(1) + q(dimensiona + 2 + k)
      end do
      if (mp_modelc == 0) then
         do k = 1, nof_species - 1
            q(dimensiona + 2 + nof_species + k) = mp_a_in(k)
         end do
      end if
      call prim2cons(0, q)
   end subroutine multispecies_reference_state

   subroutine primitive_euler_flux(axis, primitive, flux)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: axis
      real, intent(in) :: primitive(nof_variables)
      real, intent(out) :: flux(nof_variables)
      real :: conservative(gpu_max_nvar), velocity, pressure
      integer :: k
      conservative(1:nof_variables) = primitive
      call prim2cons(0, conservative)
      velocity = primitive(axis + 1); pressure = primitive(dimensiona + 2)
      do k = 1, nof_variables
         flux(k) = velocity * conservative(k)
      end do
      flux(axis + 1) = flux(axis + 1) + pressure
      flux(dimensiona + 2) = flux(dimensiona + 2) + velocity * pressure
   end subroutine primitive_euler_flux

   subroutine inflow(initcond, pox, poy, poz, inflow_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inflow in 3d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: inflow_out
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      real :: xf, yf, zf
      real :: theta_0, vtang, vradial, gammar
      real :: mp_density, mp_stiff, sum1, sum2, sum3
      real, dimension(1:gpu_max_nvar) :: vect_in
      real, dimension(gpu_max_species) :: mp_ar, mp_ie
      integer :: rg_i, rg_j
      real :: khx, t1l, vhx, amp, dvel, rgg, tt1, khi_slope, khi_b, theeta, reeta, rg_ve, rg_tr, rg_chem, rg_density, rg_ev_total, &
         & rg_rmix
      real, dimension(1:gpu_max_species) :: rg_cvs
      real, dimension(1:gpu_max_species) :: rg_tvsl, rg_ev
      real :: rg_tv, rg_ttr0, rg_tve0

      if (multispecies .eq. 1) then
         call multispecies_reference_state(inflow_out)
      else

         r = rres
         gm = gamma
         p = pres
         u = uvel
         v = vvel
         w = wvel

         if (initcond .eq. 10000) then
            if (sqrt(((poy(1) - 0.0)**2) + ((poz(1) - 0.5)**2)) .le. 0.05) then
!if (((poy(1).ge.-0.05).and.(poy(1).le.0.05)).and.((poz(1).ge.-0.05).and.(poz(1).le.0.05)))then
               p = 0.4127
               r = 5
               u = 30.0
               v = 0.00
               w = 0.00

            end if
         end if

!kinetic energy first!
         skin = (oo2) * ((u**2) + (v**2) + (w**2))
!internal energy
         ien = ((p) / ((gm - 1.00) * r))
!total energy
         e = r * (skin + ien)
!vector of conserved variables now
         inflow_out(1) = r
         inflow_out(2) = r * u
         inflow_out(3) = r * v
         inflow_out(4) = r * w
         inflow_out(5) = e

      end if

      if (swirl .eq. 1) then

         if (pox(1) .lt. -0.03) then
            xf = pox(1)
            yf = poy(1)
            zf = poz(1)
            theta_0 = atan2(zf, yf)
            vtang = 18.03750
            vradial = -12.630
            u = 0.00

            v = -vtang * sin(theta_0) + vradial * cos(theta_0)
            w = vtang * cos(theta_0) + vradial * sin(theta_0)

         else

            u = 70.060
            v = 0.00
            w = 0.00

         end if

         r = rres
         gm = gamma
         p = pres
         s = sqrt((gm * p) / (r))

!kinetic energy first!
         skin = (oo2) * ((u**2) + (v**2) + (w**2))
!internal energy
         ien = ((p) / ((gm - 1.00) * r))
!total energy
         e = r * (skin + ien)
!vector of conserved variables now
         inflow_out(1) = r
         inflow_out(2) = r * u
         inflow_out(3) = r * v
         inflow_out(4) = r * w
         inflow_out(5) = e

      end if

      if (realgas .eq. 1) then
         u = uvel
         v = vvel
         w = wvel
         p = pres
         r = rres

!first build mixture gas constant
         rg_rmix = zero
         do rg_i = 1, nof_species
            rg_rmix = rg_rmix + rg_vf(rg_i) / rg_molm(rg_i)
         end do

         rg_rmix = rgs_ru * rg_rmix

         rg_ttr0 = rg_ttr
         rg_tve0 = rg_tve

! translational-rotational internal energy
         rg_tr = 0.00
         do rg_i = 1, nof_species
            if (rg_i <= 3) then
               rg_cvs(rg_i) = (5.00 / 2.00) * rgs_ru / rg_molm(rg_i)
            else
               rg_cvs(rg_i) = (3.00 / 2.00) * rgs_ru / rg_molm(rg_i)
            end if
            rg_tr = rg_tr + rg_vf(rg_i) * rg_cvs(rg_i) * rg_ttr0
         end do

! vibrational energy
         rg_ev_total = 0.00
         do rg_i = 1, 3
            rg_ev_total = rg_ev_total + rg_vf(rg_i) * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i) / &
               & rg_tve0) - 1.00))
         end do

         rg_chem = zero

         ! chemical energy
         do rg_i = 1, nof_species
            if (rg_hzero(rg_i) .gt. 1.0e-12) then
               rg_chem = rg_chem + (rg_vf(rg_i) * rg_hzero(rg_i) / rg_molm(rg_i))
            end if
         end do
!rg_chem=zero

! kinetic energy

         skin = (oo2) * ((u**2) + (v**2) + (w**2))

         inflow_out(1) = r
         inflow_out(2) = r * u
         inflow_out(3) = r * v
         inflow_out(4) = r * w
         inflow_out(5) = r * (rg_ev_total + rg_tr + rg_chem + skin)
         inflow_out(6) = r * rg_ev_total

         do rg_i = 1, nof_species
            inflow_out(6 + rg_i) = r * rg_vf(rg_i)
         end do

      end if

   end subroutine inflow

   subroutine vect_function(pox, poy, vect_out)
      implicit none
!> @brief
!> this makes a multipliciation between two vectors
#ifdef xpu
!$omp declare target
#endif
      real, dimension(3), intent(out) :: vect_out
      real, dimension(1:dimensiona), intent(in) :: pox, poy

      vect_out(1) = (poy(2) * pox(3)) - (poy(3) * pox(2))
      vect_out(2) = (poy(3) * pox(1)) - (poy(1) * pox(3))
      vect_out(3) = (poy(1) * pox(2)) - (poy(2) * pox(1))

   end subroutine vect_function

   subroutine inflow2d(initcond, pox, poy, inflow2d_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inflow in 2d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: inflow2d_out
      integer, intent(in) :: initcond
      real, dimension(1:2), intent(in) :: pox, poy
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi, ps
      real :: xf, yf, zf, lit_a, lit_o
      real :: theta_0, vtang, vradial, gammar
      real :: mp_density, mp_stiff, sum1, sum2, sum3
      real, dimension(1:gpu_max_nvar) :: vect_in
      real, dimension(gpu_max_species) :: mp_ar, mp_ie
      integer :: rg_i, rg_j
      real :: khx, t1l, vhx, amp, dvel, rgg, tt1, khi_slope, khi_b, theeta, reeta, rg_ve, rg_tr, rg_chem, rg_density, rg_ev_total, &
         & rg_rmix
      real, dimension(1:gpu_max_species) :: rg_cvs
      real, dimension(1:gpu_max_species) :: rg_tvsl, rg_ev
      real :: rg_tv, rg_ttr0, rg_tve0

      if (multispecies .eq. 1) then
         call multispecies_reference_state(inflow2d_out)
      else

         r = rres
         gm = gamma
         p = pres
         u = uvel
         v = vvel

         if (initcond .eq. 133) then
            p = 195557.25
            r = p / (350.50 * 287.0580)
            u = 168.62
            v = 0.00

         end if

         if (initcond .eq. 10000) then
            if ((poy(1) .ge. -0.05) .and. (poy(1) .le. 0.05)) then
               p = 0.4127
               r = 5
               u = 30.0
               v = 0.00

            end if
         end if

         if (initcond .eq. 790) then
            r = (2.40 * 6**2) / ((0.4 * 6**2) + 2)
            u = (6 * sqrt(1.4)) * (70 / (2.4 * 36))
            v = 0.00
            p = (2.8 * 36 - 0.4) / (2.4)

         end if

!kinetic energy first!
         skin = (oo2) * ((u**2) + (v**2))
!internal energy
         ien = ((p) / ((gm - 1.00) * r))
!total energy
         e = r * (skin + ien)
!vector of conserved variables now
         inflow2d_out(1) = r
         inflow2d_out(2) = r * u
         inflow2d_out(3) = r * v
         inflow2d_out(4) = e

      end if

      if (realgas .eq. 1) then
         u = uvel
         v = vvel
         p = max(pres, tolsmall)

!first build mixture gas constant
         rg_rmix = zero
         do rg_i = 1, nof_species
            rg_rmix = rg_rmix + rg_vf(rg_i) / rg_molm(rg_i)
         end do

         rg_rmix = rgs_ru * rg_rmix

         rg_ttr0 = max(rg_ttr, tolsmall)
         rg_tve0 = max(rg_tve, tolsmall)
         if (rg_rmix .gt. tolsmall) then
            r = p / (rg_rmix * rg_ttr0)
         else
            r = rres
         end if

! translational-rotational internal energy
         rg_tr = 0.00
         do rg_i = 1, nof_species
            if (rg_i <= 3) then
               rg_cvs(rg_i) = (5.00 / 2.00) * rgs_ru / rg_molm(rg_i)
            else
               rg_cvs(rg_i) = (3.00 / 2.00) * rgs_ru / rg_molm(rg_i)
            end if
            rg_tr = rg_tr + rg_vf(rg_i) * rg_cvs(rg_i) * rg_ttr0
         end do

! vibrational energy
         rg_ev_total = 0.00
         do rg_i = 1, 3
            rg_ev_total = rg_ev_total + rg_vf(rg_i) * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i) / &
               & rg_tve0) - 1.00))
         end do

         rg_chem = zero

         ! chemical energy
         do rg_i = 1, nof_species
            if (rg_hzero(rg_i) .gt. 1.0e-12) then
               rg_chem = rg_chem + (rg_vf(rg_i) * rg_hzero(rg_i) / rg_molm(rg_i))
            end if
         end do
!rg_chem=zero

! kinetic energy

         skin = (oo2) * ((u**2) + (v**2))

         inflow2d_out(1) = r
         inflow2d_out(2) = r * u
         inflow2d_out(3) = r * v
         inflow2d_out(4) = r * (rg_ev_total + rg_tr + rg_chem + skin)
         inflow2d_out(5) = r * rg_ev_total

         do rg_i = 1, nof_species
            inflow2d_out(5 + rg_i) = r * rg_vf(rg_i)
         end do

      end if

   end subroutine inflow2d

   subroutine outflow2d(initcond, pox, poy, outflow2d_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 2d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: outflow2d_out
      integer, intent(in) :: initcond
      real, dimension(1:2), intent(in) :: pox, poy
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      real :: xf, yf, zf, gammar
      real, dimension(gpu_max_species) :: mp_ar, mp_ie

      if (realgas .eq. 1) then
         call inflow2d(initcond, pox, poy, outflow2d_out)
         return
      end if

      if (multispecies .eq. 1) then
         call multispecies_reference_state(outflow2d_out)
      else

         r = rres
         gm = gamma
         p = pres
         u = uvel
         v = vvel

!kinetic energy first!
         skin = (oo2) * ((u**2) + (v**2))
!internal energy
         ien = ((p) / ((gm - 1.00) * r))
!total energy
         e = r * (skin + ien)
!vector of conserved variables now
         outflow2d_out(1) = r
         outflow2d_out(2) = r * u
         outflow2d_out(3) = r * v
         outflow2d_out(4) = e

      end if

   end subroutine outflow2d
   subroutine outflow(initcond, pox, poy, poz, outflow_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: outflow_out
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      real :: xf, yf, zf, gammar
      real, dimension(gpu_max_species) :: mp_ar, mp_ie

      if (realgas == 1) then
         call inflow(initcond, pox, poy, poz, outflow_out)
         return
      end if
      if (multispecies .eq. 1) then
         call multispecies_reference_state(outflow_out)
      else

         r = rres
         gm = gamma
         p = pres

         if (initcond .eq. 977) then
            p = 101325
         end if

         u = uvel
         v = vvel
         w = wvel

!kinetic energy first!
         skin = (oo2) * ((u**2) + (v**2) + (w**2))
!internal energy
         ien = ((p) / ((gm - 1.00) * r))
!total energy
         e = r * (skin + ien)
!vector of conserved variables now
         outflow_out(1) = r
         outflow_out(2) = r * u
         outflow_out(3) = r * v
         outflow_out(4) = r * w
         outflow_out(5) = e

      end if

   end subroutine outflow

   subroutine outflow2(initcond, pox, poy, poz, outflow2_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: outflow2_out
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      real :: xf, yf, zf
      real, dimension(gpu_max_species) :: mp_ar, mp_ie

      r = rres
      gm = gamma
      p = press_outlet
      u = uvel
      v = vvel
      w = wvel

!kinetic energy first!
      skin = (oo2) * ((u**2) + (v**2) + (w**2))
!internal energy
      ien = ((p) / ((gm - 1.00) * r))
!total energy
      e = r * (skin + ien)
!vector of conserved variables now
      outflow2_out(1) = r
      outflow2_out(2) = r * u
      outflow2_out(3) = r * v
      outflow2_out(4) = r * w
      outflow2_out(5) = e

   end subroutine outflow2

   subroutine bleed2d(iconsidered, facex, pox, poy, bleed2d_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: bleed2d_out
      integer, intent(in) :: iconsidered, facex
      real, dimension(1:dimensiona), intent(in) :: pox, poy
      integer :: ibleedn, rg_i
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      real :: xf, yf, zf
      real :: cell_area
!------coefficients for bleed-----!
      real :: bleed_qsonic_s, bleed_mdotsonic_s, bleed_area, bleed_mdotsonic, bleed_region
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot

      cell_area = ielem_surf(iconsidered, facex)

!leftv is the input in conservative variables that we transform to primitive variables
      call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
      r = leftv(1)
      u = leftv(2)
      v = leftv(3)
      p = leftv(4)

!now find the bleed number conditions for this face
      ibleedn = ielem_bleedn(facex, iconsidered)

      bleed_region = sqrt(((bleed_start(ibleedn, 1) - bleed_end(ibleedn, 1))**2) + ((bleed_start(ibleedn, 2) - bleed_end(ibleedn, &
         & 2))**2))

!now compute the bleed area

      bleed_area = bleed_porosity(ibleedn) * bleed_region

!now compute the bleed m dot sonic-s mass flow rate eq.17 , https://doi.org/10.2514/1.b37474

! bleed_mdotsonic_s=bleed_area*p*(sqrt(gamma*r/p))*(((gamma+1.00)/(2))**((gamma+1)/(2*(1-gamma))))

!now compute qsonic eq. 22
      bleed_qsonic_s = 0.598 + 0.0307 * (bleed_plenum(ibleedn) / p) - 0.5936 * ((bleed_plenum(ibleedn) / p)**2)

!equation 16
      bleed_mdotsonic = bleed_qsonic_s * bleed_mdotsonic_s

!now that i have computed the bleed mass flow rate i have to distribute to the surface area of this tagged cell in the normal direction

      call prim2cons2(n, leftv, rightv)

      rightv(1:nof_variables) = leftv(1:nof_variables)

      call rotatef2d(n, cright_rot, rightv, angle1, angle2)

      cright_rot(2) = bleed_mdotsonic / cell_area

      call rotateb2d(n, rightv, cright_rot, angle1, angle2)

      do rg_i = 1, nof_variables
         bleed2d_out(rg_i) = rightv(rg_i)
      end do

   end subroutine bleed2d

   subroutine bleed3d(iconsidered, facex, pox, poy, poz, bleed3d_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the outflow in 3d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(out) :: bleed3d_out
      integer, intent(in) :: iconsidered, facex
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      integer :: ibleedn, rg_i
      real :: p, u, v, w, e, r, s, gm, skin, ien, pi
      real :: xf, yf, zf
      real :: cell_area
!------coefficients for bleed-----!
      real :: bleed_qsonic_s, bleed_mdotsonic_s, bleed_area, bleed_mdotsonic, bleed_region
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot

      do rg_i = 1, nof_variables
         bleed3d_out(rg_i) = rightv(rg_i)
      end do

   end subroutine bleed3d

   subroutine pass_inlet(initcond, pox, poy, poz, pass_inlet_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inlet for a passive scalar
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:passivescalar), intent(out) :: pass_inlet_out
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      integer :: rg_i

      do rg_i = 1, passivescalar
         pass_inlet_out(rg_i) = 1.00 * rres
      end do

   end subroutine pass_inlet

   subroutine pass_inlet2d(initcond, pox, poy, pass_inlet2d_out)
      implicit none
!> @brief
!> this function applies a prescribed boundary condition to  the inlet for a passive scalar in 2d
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:passivescalar), intent(out) :: pass_inlet2d_out
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy
      integer :: rg_i

      do rg_i = 1, passivescalar
         pass_inlet2d_out(rg_i) = 1.00
      end do

   end subroutine pass_inlet2d

   subroutine shear_x(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      ! This diagnostic uses the full 3D velocity-gradient tensor.
      real, dimension(3, 3) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

      vortet1(1:3, 1:3) = rec_grads(iconsidered, 1:3, 1:3)

      ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
      vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
      wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))
      leftv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)
      rightv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)

      call get_visc_conduct(n, leftv, rightv, viscl, laml)

      ssx = zero; ssp = zero; ssy = zero; ssz = zero

      tauxx = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
      tauyy = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
      tauzz = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy
      tauyx = (uy + vx)
      tauzx = (wx + uz)
      tauzy = (vz + wy)
      ssx = (viscl(1) * ((nx * tauxx) + (ny * tauyx) + (nz * tauzx)))
      ssy = (viscl(1) * ((nx * tauyx) + (ny * tauyy) + (nz * tauzy)))
      ssz = (viscl(1) * ((nx * tauzx) + (ny * tauzy) + (nz * tauzz)))

      if (rframe .eq. 0) then
         shear_temp = -ssx / (0.5 * rres * ufreestream * ufreestream)
      else
         shear_temp = -ssx / (0.5 * rres * v_ref * v_ref)
      end if

   end subroutine shear_x

   subroutine shear_y(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in y-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      ! This diagnostic uses the full 3D velocity-gradient tensor.
      real, dimension(3, 3) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

      vortet1(1:3, 1:3) = rec_grads(iconsidered, 1:3, 1:3)

      ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
      vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
      wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))
      leftv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)
      rightv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)

      call get_visc_conduct(n, leftv, rightv, viscl, laml)

      ssx = zero; ssp = zero; ssy = zero; ssz = zero

      tauxx = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
      tauyy = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
      tauzz = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy
      tauyx = (uy + vx)
      tauzx = (wx + uz)
      tauzy = (vz + wy)
      ssx = (viscl(1) * ((nx * tauxx) + (ny * tauyx) + (nz * tauzx)))
      ssy = (viscl(1) * ((nx * tauyx) + (ny * tauyy) + (nz * tauzy)))
      ssz = (viscl(1) * ((nx * tauzx) + (ny * tauzy) + (nz * tauzz)))

      if (rframe .eq. 0) then
         shear_temp = -ssy / (0.5 * rres * ufreestream * ufreestream)
      else
         shear_temp = -ssy / (0.5 * rres * v_ref * v_ref)
      end if

   end subroutine shear_y

   subroutine shear_z(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in z-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      ! This diagnostic uses the full 3D velocity-gradient tensor.
      real, dimension(3, 3) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

      vortet1(1:3, 1:3) = rec_grads(iconsidered, 1:3, 1:3)

      ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
      vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
      wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))
      leftv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)
      rightv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)

      call get_visc_conduct(n, leftv, rightv, viscl, laml)

      ssx = zero; ssp = zero; ssy = zero; ssz = zero

      tauxx = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
      tauyy = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
      tauzz = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy
      tauyx = (uy + vx)
      tauzx = (wx + uz)
      tauzy = (vz + wy)
      ssx = (viscl(1) * ((nx * tauxx) + (ny * tauyx) + (nz * tauzx)))
      ssy = (viscl(1) * ((nx * tauyx) + (ny * tauyy) + (nz * tauzy)))
      ssz = (viscl(1) * ((nx * tauzx) + (ny * tauzy) + (nz * tauzz)))

      if (rframe .eq. 0) then
         shear_temp = -ssz / (0.5 * rres * ufreestream * ufreestream)
      else
         shear_temp = -ssz / (0.5 * rres * v_ref * v_ref)
      end if

   end subroutine shear_z

   subroutine realgas_wall_heat_transport(n, leftv, ktr_qflux, kve_qflux)
      implicit none
!> @brief
!> Return real-gas translational and vibrational conductivities for wall heat flux output.
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(in) :: leftv
      real, intent(out) :: ktr_qflux, kve_qflux
      real, dimension(1:gpu_max_nvar) :: div_state, prim_state
      real, dimension(1:gpu_max_species) :: rgs_x, rgs_y, rgs_deff, rgs_mu_i, rgs_ktr_i
      real, dimension(1:gpu_max_species) :: rgs_htr_i, rgs_hvib_i
      real, dimension(1:gpu_max_species, 1:gpu_max_species) :: rgs_dij
      real :: rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho, rgs_mu_mix, rgs_ktr_mix, rgs_kve
      real :: sumy, sumx, mp_pinfl, gammal
      integer :: rg_i

      div_state = zero
      prim_state = zero
      rgs_x = zero
      rgs_y = zero
      div_state(1:nof_variables) = leftv(1:nof_variables)
      prim_state(1:nof_variables) = leftv(1:nof_variables)

      call cons2div(n, div_state, mp_pinfl, gammal)
      call cons2prim(n, prim_state, mp_pinfl, gammal)

      rgs_rho = max(div_state(1), tolsmall)
      rgs_p_pa = max(prim_state(dimensiona + 2), tolsmall)

      rgs_ttr = div_state(dimensiona + 2)
      rgs_tv = div_state(dimensiona + 3)

      rgs_ttr = max(rgs_ttr, 50.00)
      rgs_tv = max(rgs_tv, 50.00)

      sumy = zero
      do rg_i = 1, nof_species
         rgs_y(rg_i) = max(div_state(dimensiona + 3 + rg_i), zero)
         sumy = sumy + rgs_y(rg_i)
      end do

      if (sumy .gt. 1.0e-30) then
         do rg_i = 1, nof_species
            rgs_y(rg_i) = rgs_y(rg_i) / sumy
         end do
      else
         do rg_i = 1, nof_species
            rgs_y(rg_i) = 1.00 / real(nof_species)
         end do
      end if

      sumx = zero
      do rg_i = 1, nof_species
         rgs_x(rg_i) = rgs_y(rg_i) / rg_molm(rg_i)
         sumx = sumx + rgs_x(rg_i)
      end do

      if (sumx .gt. 1.0e-300) then
         do rg_i = 1, nof_species
            rgs_x(rg_i) = rgs_x(rg_i) / sumx
         end do
      else
         do rg_i = 1, nof_species
            rgs_x(rg_i) = 1.00 / real(nof_species)
         end do
      end if

      call compute_real_gas_diffusion(rgs_ttr, rgs_tv, rgs_p_pa, rgs_rho, rgs_x, rgs_y, &
         & rgs_dij, rgs_deff, rgs_mu_i, rgs_mu_mix, &
         & rgs_ktr_i, rgs_ktr_mix, rgs_kve, &
         & rgs_htr_i, rgs_hvib_i)

      ktr_qflux = rgs_ktr_mix
      kve_qflux = rgs_kve

   end subroutine realgas_wall_heat_transport

   subroutine heat_x(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp, lam_qflux, ktr_qflux, kve_qflux
      real, dimension(1:gpu_max_dim, 1:gpu_max_dim) :: vortet1
      integer :: i, k, j, kmaxe, gqi_points, nnd, im
      real, dimension(1:gpu_max_nvar) :: leftv
      real, dimension(1:3) :: temp_grad, temp_grad_v
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz, surface_temp
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:gpu_max_dim, 1:gpu_max_qp_face) :: qpoints2d
      real, dimension(1:gpu_max_qp_face) :: wequa2d

      i = iconsidered
      ssx = zero; ssy = zero; ssz = zero
      j = facex
      angle1 = ielem_faceanglex(i, j)
      angle2 = ielem_faceangley(i, j)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))

      select case (ielem_types_faces(i, j))
      case (5)
         gqi_points = qp_quad_n

         if (reduce_comp .eq. 1) then
            wequa2d = 1.00;
         else
            nnd = 4
            do k = 1, nnd
               vext(k, 1:dims) = inoder4_cord(1:dims, ielem_nodes_faces(i, j, k))
            end do
            call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
         end if
         surface_temp = ielem_surf(i, j)

      case (6)
         gqi_points = qp_triangle_n

         if (reduce_comp .eq. 1) then
            wequa2d = 1.00;
         else
            nnd = 3
            do k = 1, nnd
               vext(k, 1:dims) = inoder4_cord(1:dims, ielem_nodes_faces(i, j, k))
            end do
            call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
         end if
         surface_temp = ielem_surf(i, j)

      end select

      do im = 1, gqi_points
         temp_grad(1:3) = rec_uleftv(i, 1:3, dimensiona + 1, j, im)
         if (realgas .eq. 1) temp_grad_v(1:3) = rec_uleftv(i, 1:3, dimensiona + 2, j, im)
         if (dg .eq. 1) then
            leftv(1:nof_variables) = rec_uleft_dg(i, 1:nof_variables, j, im)
            rightv(1:nof_variables) = rec_uleft_dg(i, 1:nof_variables, j, im)

         else
            leftv(1:nof_variables) = rec_uleft(i, :, j, im)
            rightv(1:nof_variables) = rec_uleft(i, :, j, im)
         end if

         if (realgas .eq. 1) then
            call realgas_wall_heat_transport(n, leftv, ktr_qflux, kve_qflux)
            ssx = ssx + (ktr_qflux * temp_grad(1) + kve_qflux * temp_grad_v(1)) * wequa2d(im) * nx
            ssy = ssy + (ktr_qflux * temp_grad(2) + kve_qflux * temp_grad_v(2)) * wequa2d(im) * ny
            ssz = ssz + (ktr_qflux * temp_grad(3) + kve_qflux * temp_grad_v(3)) * wequa2d(im) * nz
         else
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then

               turbmv(1) = rec_uleftturb(i, 1, j, im)

               turbmv(2) = rec_uleftturb(i, 1, j, im)
               eddyfl(2) = turbmv(1);
               eddyfr(2) = turbmv(2)
               call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            if (turbulence .eq. 1) then
               lam_qflux = laml(3)
            else
               lam_qflux = laml(1)
            end if

            ssx = ssx + lam_qflux * temp_grad(1) * wequa2d(im) * nx !*surface_temp
            ssy = ssy + lam_qflux * temp_grad(2) * wequa2d(im) * ny !*surface_temp
            ssz = ssz + lam_qflux * temp_grad(3) * wequa2d(im) * nz !*surface_temp
         end if
      end do

      shear_temp = -(ssx + ssy + ssz)

   end subroutine heat_x

   subroutine heat_x2d(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ssx, ssy, ssz, ssp, ktr_qflux, kve_qflux
      real, dimension(1:gpu_max_dim, 1:gpu_max_dim) :: vortet1
      integer :: i, k, j, kmaxe, gqi_points, nnd, im
      real, dimension(1:gpu_max_nvar) :: leftv
      real, dimension(1:2) :: temp_grad, temp_grad_v
      real :: mp_pinfl, gammal, lam_qflux
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz, surface_temp
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:gpu_max_dim, 1:gpu_max_qp_face) :: qpoints2d
      real, dimension(1:gpu_max_qp_face) :: wequa2d
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      i = iconsidered
      ssx = zero; ssy = zero
      j = facex
      nx = ielem_faceanglex(i, j)
      ny = ielem_faceangley(i, j)

      gqi_points = qp_line_n
      if (reduce_comp .eq. 1) then
         wequa2d = 1.00;
      else
         nnd = 2
         do k = 1, nnd
            vext(k, 1:dims) = inoder4_cord(1:dims, ielem_nodes_faces(i, j, k))
         end do

         call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
      end if
      surface_temp = ielem_surf(i, j)

      do im = 1, gqi_points
         temp_grad(1:2) = rec_uleftv(i, 1:2, dimensiona + 1, j, im)
         if (realgas .eq. 1) temp_grad_v(1:2) = rec_uleftv(i, 1:2, dimensiona + 2, j, im)

         if (dg .eq. 1) then
            leftv(1:nof_variables) = rec_uleft_dg(i, 1:nof_variables, j, im)
            rightv(1:nof_variables) = rec_uleft_dg(i, 1:nof_variables, j, im)

         else
            leftv(1:nof_variables) = rec_uleft(i, :, j, im)
            rightv(1:nof_variables) = rec_uleft(i, :, j, im)
         end if

         if (realgas .eq. 1) then
            call realgas_wall_heat_transport(n, leftv, ktr_qflux, kve_qflux)
            ssx = ssx + (ktr_qflux * temp_grad(1) + kve_qflux * temp_grad_v(1)) * wequa2d(im) * nx
            ssy = ssy + (ktr_qflux * temp_grad(2) + kve_qflux * temp_grad_v(2)) * wequa2d(im) * ny
         else
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then

               turbmv(1) = rec_uleftturb(i, 1, j, im)

               turbmv(2) = rec_uleftturb(i, 1, j, im)
               eddyfl(2) = turbmv(1);
               eddyfr(2) = turbmv(2)
               call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            if (turbulence .eq. 1) then
               lam_qflux = laml(3)
            else
               lam_qflux = laml(1)
            end if

            ssx = ssx + lam_qflux * temp_grad(1) * wequa2d(im) * nx !*surface_temp
            ssy = ssy + lam_qflux * temp_grad(2) * wequa2d(im) * ny !*surface_temp

         end if

      end do

      shear_temp = -(ssx + ssy)

   end subroutine heat_x2d

   subroutine heat_y2d(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ssx, ssy, ssz, ssp
      real, dimension(1:gpu_max_dim, 1:gpu_max_dim) :: vortet1
      integer :: i, k, j, kmaxe, gqi_points, nnd, im
      real, dimension(1:gpu_max_nvar) :: leftv
      real, dimension(1:3) :: temp_grad
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz, surface_temp
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:gpu_max_dim, 1:gpu_max_qp_face) :: qpoints2d
      real, dimension(1:gpu_max_qp_face) :: wequa2d

      i = iconsidered
      ssy = zero
      j = facex
      nx = ielem_faceanglex(i, j)
      ny = ielem_faceangley(i, j)

      gqi_points = qp_line_n
      if (reduce_comp .eq. 1) then
         wequa2d = 1.00;
      else
         nnd = 2
         do k = 1, nnd
            vext(k, 1:dims) = inoder4_cord(1:dims, ielem_nodes_faces(i, j, k))
         end do

         call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
      end if
      surface_temp = ielem_surf(i, j)

      do im = 1, gqi_points
         temp_grad(1:2) = rec_uleftv(i, 1:2, dimensiona + 1, j, im)
         ssy = ssy - 0.026 * temp_grad(2) * wequa2d(im) * surface_temp
      end do

      shear_temp = ssy

   end subroutine heat_y2d

   subroutine shear_x_av(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the average shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      ! This diagnostic uses the full 3D velocity-gradient tensor.
      real, dimension(3, 3) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

!  if (rungekutta.eq.4)then
!               ind1=7
!               else
!               ind1=5
!               end if

      vortet1(1:3, 1:3) = rec_gradsav(iconsidered, 1:3, 1:3)

      ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
      vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
      wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))
      leftv(1:nof_variables) = u_c_val(iconsidered, ind1, 1:nof_variables)
      rightv(1:nof_variables) = u_c_val(iconsidered, ind1, 1:nof_variables)

      call get_visc_conduct(n, leftv, rightv, viscl, laml)

      ssx = zero; ssp = zero; ssy = zero; ssz = zero

      tauxx = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
      tauyy = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
      tauzz = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy
      tauyx = (uy + vx)
      tauzx = (wx + uz)
      tauzy = (vz + wy)
      ssx = (viscl(1) * ((nx * tauxx) + (ny * tauyx) + (nz * tauzx)))
      ssy = (viscl(1) * ((nx * tauyx) + (ny * tauyy) + (nz * tauzy)))
      ssz = (viscl(1) * ((nx * tauzx) + (ny * tauzy) + (nz * tauzz)))

      shear_temp = -ssx / (0.5 * rres * ufreestream * ufreestream)

   end subroutine shear_x_av

   subroutine shear_y_av(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the average shear stresses in y-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      ! This diagnostic uses the full 3D velocity-gradient tensor.
      real, dimension(3, 3) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

      vortet1(1:3, 1:3) = rec_gradsav(iconsidered, 1:3, 1:3)

      ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
      vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
      wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))
      leftv(1:nof_variables) = u_c_val(iconsidered, ind1, 1:nof_variables)
      rightv(1:nof_variables) = u_c_val(iconsidered, ind1, 1:nof_variables)

      call get_visc_conduct(n, leftv, rightv, viscl, laml)

      ssx = zero; ssp = zero; ssy = zero; ssz = zero

      tauxx = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
      tauyy = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
      tauzz = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy
      tauyx = (uy + vx)
      tauzx = (wx + uz)
      tauzy = (vz + wy)
      ssx = (viscl(1) * ((nx * tauxx) + (ny * tauyx) + (nz * tauzx)))
      ssy = (viscl(1) * ((nx * tauyx) + (ny * tauyy) + (nz * tauzy)))
      ssz = (viscl(1) * ((nx * tauzx) + (ny * tauzy) + (nz * tauzz)))

      shear_temp = -ssy / (0.5 * rres * ufreestream * ufreestream)

   end subroutine shear_y_av

   subroutine shear_z_av(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the avergage shear stresses in z-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      ! This diagnostic uses the full 3D velocity-gradient tensor.
      real, dimension(3, 3) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

      vortet1(1:3, 1:3) = rec_gradsav(iconsidered, 1:3, 1:3)

      ux = vortet1(1, 1); uy = vortet1(1, 2); uz = vortet1(1, 3)
      vx = vortet1(2, 1); vy = vortet1(2, 2); vz = vortet1(2, 3)
      wx = vortet1(3, 1); wy = vortet1(3, 2); wz = vortet1(3, 3)

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))
      leftv(1:nof_variables) = u_c_val(iconsidered, ind1, 1:nof_variables)
      rightv(1:nof_variables) = u_c_val(iconsidered, ind1, 1:nof_variables)

      call get_visc_conduct(n, leftv, rightv, viscl, laml)

      ssx = zero; ssp = zero; ssy = zero; ssz = zero

      tauxx = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
      tauyy = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
      tauzz = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy
      tauyx = (uy + vx)
      tauzx = (wx + uz)
      tauzy = (vz + wy)
      ssx = (viscl(1) * ((nx * tauxx) + (ny * tauyx) + (nz * tauzx)))
      ssy = (viscl(1) * ((nx * tauyx) + (ny * tauyy) + (nz * tauzy)))
      ssz = (viscl(1) * ((nx * tauzx) + (ny * tauzy) + (nz * tauzz)))

      shear_temp = -ssz / (0.5 * rres * ufreestream * ufreestream)

   end subroutine shear_z_av

   subroutine shear_x2d(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      real, dimension(1:gpu_max_dim, 1:gpu_max_dim) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml

      integer :: gqi_points, im
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:gpu_max_dim, 1:gpu_max_qp_face) :: qpoints2d
      real, dimension(1:gpu_max_qp_face) :: wequa2d
      ssx = zero; ssp = zero; ssy = zero;
      gqi_points = qp_line_n
      call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)

      do im = 1, gqi_points
         if (ielem_ggs(iconsidered) .eq. 1) then
            vortet1(1:2, 1:2) = rec_grads(iconsidered, 1:2, 1:2)
         else
            vortet1(1, 1:2) = rec_uleftv(iconsidered, 1:2, 1, facex, im)
            vortet1(2, 1:2) = rec_uleftv(iconsidered, 1:2, 2, facex, im)

         end if

         ux = vortet1(1, 1); uy = vortet1(1, 2)
         vx = vortet1(2, 1); vy = vortet1(2, 2)

         angle1 = ielem_faceanglex(iconsidered, facex)
         angle2 = ielem_faceangley(iconsidered, facex)
         nx = angle1
         ny = angle2

         leftv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)
         rightv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)

         call get_visc_conduct(n, leftv, rightv, viscl, laml)

         ssx = zero; ssp = zero; ssy = zero; ssz = zero

         tauxx = 2.00 * ux
         tauyy = 2.00 * vy

         tauyx = (uy + vx)

         ssx = ssx + ((viscl(1) * ((nx * tauxx) + (ny * tauyx))) * wequa2d(im))
         ssy = ssy + ((viscl(1) * ((nx * tauyx) + (ny * tauyy))) * wequa2d(im))
      end do

      shear_temp = -ssx / (0.5 * rres * ufreestream * ufreestream)

   end subroutine shear_x2d

   subroutine shear_y2d(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the shear stresses in y-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, tauxx, tauyy, tauzz, tauyx, tauzx, tauzy
      real :: ssx, ssy, ssz, ssp
      real, dimension(1:gpu_max_dim, 1:gpu_max_dim) :: vortet1
      real, dimension(1:gpu_max_nvar) :: leftv
      real :: mp_pinfl, gammal
      real, dimension(1:gpu_max_nvar) :: rightv
      real :: mp_pinfr, gammar
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:4) :: viscl, laml
      integer :: gqi_points, im
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:gpu_max_dim, 1:gpu_max_qp_face) :: qpoints2d
      real, dimension(1:gpu_max_qp_face) :: wequa2d
      ssx = zero; ssp = zero; ssy = zero;
      gqi_points = qp_line_n
      call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)

      do im = 1, gqi_points
         if (ielem_ggs(iconsidered) .eq. 1) then
            vortet1(1:2, 1:2) = rec_grads(iconsidered, 1:2, 1:2)
         else
            vortet1(1, 1:2) = rec_uleftv(iconsidered, 1:2, 1, facex, im)
            vortet1(2, 1:2) = rec_uleftv(iconsidered, 1:2, 2, facex, im)

         end if

         ux = vortet1(1, 1); uy = vortet1(1, 2)
         vx = vortet1(2, 1); vy = vortet1(2, 2)

         angle1 = ielem_faceanglex(iconsidered, facex)
         angle2 = ielem_faceangley(iconsidered, facex)
         nx = angle1
         ny = angle2

         leftv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)
         rightv(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)

         call get_visc_conduct(n, leftv, rightv, viscl, laml)

         ssx = zero; ssp = zero; ssy = zero; ssz = zero

         tauxx = 2.00 * ux
         tauyy = 2.00 * vy

         tauyx = (uy + vx)

         ssy = ssy + ((viscl(1) * ((nx * tauyx) + (ny * tauyy))) * wequa2d(im))
      end do

      shear_temp = -ssy / (0.5 * rres * ufreestream * ufreestream)

   end subroutine shear_y2d

   subroutine shear_x2d_av(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the average shear stresses in x-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp

      shear_temp = 0.00

   end subroutine shear_x2d_av

   subroutine shear_y2d_av(iconsidered, facex, shear_temp)
      implicit none
!> @brief
!> this subroutine computes the average shear stresses in y-axis
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered, facex
      real, intent(inout) :: shear_temp

      shear_temp = 0.00

   end subroutine shear_y2d_av

   subroutine get_visc_conduct(n, leftv, rightv, viscl, laml)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: leftv(nof_variables), rightv(nof_variables)
      real, intent(out) :: viscl(4), laml(4)
      call molecular_transport(n, leftv, viscl(1), laml(1))
      call molecular_transport(n, rightv, viscl(2), laml(2))
      viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
   end subroutine get_visc_conduct

   subroutine molecular_transport(n, q, mu, conductivity)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: q(nof_variables)
      real, intent(out) :: mu, conductivity
      real :: temp, cp, gm, tref, ke
      integer :: k
      if (realgas == 1) then
         call realgas_mu_ktr(q, mu, conductivity)
      elseif (multispecies == 1) then
         call multispecies_mixtures(q, temp, mu, conductivity, cp, gm)
      else
         ke = 0.0e0
         do k = 2, dimensiona + 1
            ke = ke + (q(k) / q(1))**2
         end do
         temp = (gamma - 1.0e0) * (q(dimensiona + 2) / q(1) - 0.5e0 * ke) / r_gas
         if (initcond .eq. 4440 .and. duct_viscosity_ref .gt. 0.0e0) then
            tref = total_temperature_inlet
         else
            tref = pres / (rres * r_gas)
         end if
         mu = visc * (temp / tref) * sqrt(temp / tref) * (tref + suther * tref) / (temp + suther * tref)
         conductivity = mu * r_gas * gamma / (prandtl * (gamma - 1.0e0))
      end if
   end subroutine molecular_transport

   subroutine vortexcalc(n)
      implicit none
!> @brief
!> this subroutine computes the q-criterion
      integer, intent(in) :: n
      integer :: kmaxe, i, ihgt, ihgj
      real :: snorm, onorm
      kmaxe = xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n, kmaxe) &
!$omp& map(alloc: rec_grads, ielem_vortex) &
!$omp& private(i, ihgt, ihgj, snorm, onorm)
#else
!$omp do &
!$omp& private(i, ihgj, ihgt, onorm, snorm)
#endif
      do i = 1, kmaxe
         block
            real, dimension(3, 3) :: tvort, svort, ovort
            real, dimension(1:3, 1:3) :: vortet1


            vortet1(1:3, 1:3) = rec_grads(i, 1:3, 1:3)

            do ihgt = 1, 3; do ihgj = 1, 3
                  tvort(ihgt, ihgj) = vortet1(ihgj, ihgt)
               end do; end do
            svort = 0.50 * (vortet1 + tvort)
            ovort = 0.50 * (vortet1 - tvort)
            snorm = sqrt((svort(1, 1) * svort(1, 1)) + (svort(1, 2) * svort(1, 2)) + &
               & (svort(1, 3) * svort(1, 3)) + (svort(2, 1) * svort(2, 1)) + (svort(2, 2) * svort(2, 2)) + (svort(2, 3) * svort(2, &
               & 3)) &
               & + (svort(3, 1) * svort(3, 1)) + (svort(3, 2) * svort(3, 2)) + (svort(3, 3) * svort(3, 3)))
            onorm = sqrt((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + (ovort(1, 3) * ovort(1, 3)) + &
               & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2)) + (ovort(2, 3) * ovort(2, 3)) + (ovort(3, 1) * ovort(3, &
               & 1)) + &
               & (ovort(3, 2) * ovort(3, 2)) + (ovort(3, 3) * ovort(3, 3)))

            ielem_vortex(i, 1) = (0.50 * ((onorm**2) - (snorm**2)))

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine vortexcalc

   subroutine enstrophy_calc(n)
      implicit none
!> @brief
!> this subroutine computes the q-criterion
      integer, intent(in) :: n
      integer :: kmaxe, i, ihgt, ihgj
      real :: snorm, onorm
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz
      real :: rho_i, u_i, v_i, w_i, kinetic_i, p_i, t_i, t_ref, t_ratio, mu_i
      ! This diagnostic uses the full 3D velocity-gradient tensor.

      kmaxe = xmpielrank(n)
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n, kmaxe, nof_variables, dimensiona, initcond, boundtype, gamma, oo2, r_gas, pres, rres, visc, suther, zero) &
!$omp& firstprivate(duct_viscosity_ref, total_temperature_inlet) &
!$omp& map(alloc: rec_grads, u_c_val, ielem_vortex, ielem_totvolume) &
!$omp& private(i, ihgt, ihgj, snorm, onorm) &
!$omp& private(ux, uy, uz, vx, vy, vz, wx, wy, wz) &
!$omp& private(rho_i, u_i, v_i, w_i, kinetic_i, p_i, t_i, t_ref, t_ratio, mu_i)
#else
!$omp do &
!$omp& private(i, ihgj, ihgt, kinetic_i, mu_i, onorm, p_i, rho_i, snorm, t_i, t_ratio, t_ref) &
!$omp& private(u_i, ux, uy, uz, v_i, vx, vy, vz, w_i, wx, wy, wz)
#endif
      do i = 1, kmaxe
         block
            real, dimension(3, 3) :: tvort, ovort
            real, dimension(3, 3) :: vortet1


            vortet1(1:3, 1:3) = rec_grads(i, 1:3, 1:3)

            do ihgt = 1, 3; do ihgj = 1, 3
                  tvort(ihgt, ihgj) = vortet1(ihgj, ihgt)
               end do; end do

            ovort = (vortet1 - tvort)
            onorm = ((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + (ovort(1, 3) * ovort(1, 3)) + &
               & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2)) + (ovort(2, 3) * ovort(2, 3)) + (ovort(3, 1) * ovort(3, &
               & 1)) + &
               & (ovort(3, 2) * ovort(3, 2)) + (ovort(3, 3) * ovort(3, 3)))

            if (boundtype .eq. 1) then

               ielem_vortex(i, 2) = (0.50 * (onorm * u_c_val(i, 1, 1)))

            else

               ux = rec_grads(i, 1, 1); uy = rec_grads(i, 1, 2); uz = rec_grads(i, 1, 3);
               vx = rec_grads(i, 2, 1); vy = rec_grads(i, 2, 2); vz = rec_grads(i, 2, 3);
               wx = rec_grads(i, 3, 1); wy = rec_grads(i, 3, 2); wz = rec_grads(i, 3, 3);
               rho_i = u_c_val(i, 1, 1)
               u_i = u_c_val(i, 1, 2) / rho_i
               v_i = u_c_val(i, 1, 3) / rho_i
               w_i = zero
               if (dimensiona .eq. 3) w_i = u_c_val(i, 1, 4) / rho_i
               kinetic_i = (u_i * u_i) + (v_i * v_i) + (w_i * w_i)
               p_i = (gamma - 1.00) * (u_c_val(i, 1, dimensiona + 2) - oo2 * rho_i * kinetic_i)
               t_i = p_i / (rho_i * r_gas)
               if (initcond .eq. 4440 .and. duct_viscosity_ref .gt. 0.0e0) then
                  t_ref = total_temperature_inlet
               else
                  t_ref = pres / (rres * r_gas)
                  if ((initcond .eq. 95) .and. (boundtype .ne. 1)) t_ref = 1.00 / (gamma * 1.250 * 1.250 * r_gas)
               end if
               t_ratio = t_i / t_ref
               mu_i = visc * (t_ratio * sqrt(t_ratio)) * ((t_ref + (suther * t_ref)) / (t_i + (suther * t_ref)))

               snorm = ((wy - vz)**2) + ((uz - wx)**2) + ((vx - uy)**2)
               ielem_vortex(i, 2) = ielem_totvolume(i) * snorm * mu_i
               ielem_vortex(i, 3) = (4.0 / 3.0) * mu_i * ((ux + vy + wz)**2) * ielem_totvolume(i)

            end if

         end block
      end do

#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine enstrophy_calc

   subroutine vortexcalc2d(n)
      implicit none
!> @brief
!> this subroutine computes the q criterion for 2d
      integer, intent(in) :: n
      integer :: kmaxe, i, ihgt, ihgj
      real :: snorm, onorm

      kmaxe = xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n, kmaxe) &
!$omp& map(alloc: rec_grads, ielem_vortex) &
!$omp& private(i, ihgt, ihgj, snorm, onorm)
#else
!$omp do &
!$omp& private(i, ihgj, ihgt, onorm, snorm)
#endif
      do i = 1, kmaxe
         block
            real, dimension(2, 2) :: tvort, svort, ovort
            real, dimension(1:gpu_max_dim, 1:gpu_max_dim) :: vortet1


            vortet1(1:2, 1:2) = rec_grads(i, 1:2, 1:2)

            do ihgt = 1, 2; do ihgj = 1, 2
                  tvort(ihgt, ihgj) = vortet1(ihgj, ihgt)
               end do; end do
            svort(1:2, 1:2) = 0.50 * (vortet1(1:2, 1:2) + tvort(1:2, 1:2))
            ovort(1:2, 1:2) = 0.50 * (vortet1(1:2, 1:2) - tvort(1:2, 1:2))
            snorm = sqrt((svort(1, 1) * svort(1, 1)) + (svort(1, 2) * svort(1, 2)) + &
               & (svort(2, 1) * svort(2, 1)) + (svort(2, 2) * svort(2, 2)))
            onorm = sqrt((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + &
               & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2)))

            ielem_vortex(i, 1) = (0.50 * ((onorm**2) - (snorm**2)))

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine vortexcalc2d

   subroutine boundarys_realgas_outflow(n, initcond, pox, poy, poz, leftv, rightv, nx, ny, nz)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      real, dimension(1:nof_variables), intent(in) :: leftv
      real, dimension(1:nof_variables), intent(inout) :: rightv
      real, intent(in) :: nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: tempvect1, tempvect2
      real, dimension(1:gpu_max_species) :: y_s_g
      real :: mp_pinfl, gammal, agrt, u, v, w, un, mn, rho_l, rho_g, p_l, p_g, u_g, v_g, w_g, rg_rmix, rmix, t_g, tv_g
      real :: rg_tr, rg_ev_total, rg_chem, skin1, rg_cv
      integer :: rg_i

      tempvect1 = zero
      tempvect2 = zero
      tempvect1(1:nof_variables) = leftv(1:nof_variables)
      tempvect2(1:nof_variables) = leftv(1:nof_variables)

      call cons2prim(n, tempvect1, mp_pinfl, gammal)
      agrt = sqrt((tempvect1(5) + mp_pinfl) * gammal / tempvect1(1))
      rho_l = tempvect1(1)
      u = tempvect1(2)
      v = tempvect1(3)
      w = tempvect1(4)
      p_l = tempvect1(5)
      un = u * nx + v * ny + w * nz

      call cons2div(n, tempvect2, mp_pinfl, gammal)

      if (un .le. 0.00) then
         call inflow(initcond, pox, poy, poz, rightv)
      else
         mn = un / agrt

         if (mn .ge. 1.00) then
            rightv(1:nof_variables) = leftv(1:nof_variables)
         else
            p_g = max(pres, tolsmall)
            rho_g = rho_l + (p_g - p_l) / ((agrt * agrt) + tolsmall)
            rho_g = max(rho_g, tolsmall)
            u_g = u + (nx * (p_l - p_g)) / ((agrt + tolsmall) * (rho_l + tolsmall))
            v_g = v + (ny * (p_l - p_g)) / ((agrt + tolsmall) * (rho_l + tolsmall))
            w_g = w + (nz * (p_l - p_g)) / ((agrt + tolsmall) * (rho_l + tolsmall))
            y_s_g(1:nof_species) = tempvect1(dimensiona + 4:nof_variables)

            rg_rmix = zero
            do rg_i = 1, nof_species
               rg_rmix = rg_rmix + y_s_g(rg_i) / rg_molm(rg_i)
            end do

            rg_rmix = rgs_ru * rg_rmix
            rmix = rg_rmix
            t_g = p_g / (rho_g * rmix)
            tv_g = tempvect2(6)

            rg_tr = zero
            do rg_i = 1, nof_species
               if (rg_i .le. 3) then
                  rg_cv = (5.00 / 2.00) * rgs_ru / rg_molm(rg_i)
               else
                  rg_cv = (3.00 / 2.00) * rgs_ru / rg_molm(rg_i)
               end if
               rg_tr = rg_tr + y_s_g(rg_i) * rg_cv * t_g
            end do

            rg_ev_total = zero
            do rg_i = 1, 3
               rg_ev_total = rg_ev_total + y_s_g(rg_i) * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i) / &
                  & tv_g) - 1.00))
            end do

            rg_chem = zero
            do rg_i = 1, nof_species
               if (rg_hzero(rg_i) .gt. 1.0e-12) then
                  rg_chem = rg_chem + (y_s_g(rg_i) * rg_hzero(rg_i) / rg_molm(rg_i))
               end if
            end do

            skin1 = oo2 * ((u_g**2) + (v_g**2) + (w_g * w_g))

            rightv(1) = rho_g
            rightv(2) = rho_g * u_g
            rightv(3) = rho_g * v_g
            rightv(4) = rho_g * w_g
            rightv(5) = rho_g * (rg_ev_total + rg_tr + rg_chem + skin1)
            rightv(6) = rho_g * rg_ev_total

            do rg_i = 1, nof_species
               rightv(6 + rg_i) = rho_g * y_s_g(rg_i)
            end do
         end if
      end if

   end subroutine boundarys_realgas_outflow

   subroutine boundarys2d_realgas_outflow(n, initcond, pox, poy, leftv, rightv, nx, ny)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      integer, intent(in) :: initcond
      real, dimension(1:dimensiona), intent(in) :: pox, poy
      real, dimension(1:nof_variables), intent(in) :: leftv
      real, dimension(1:nof_variables), intent(inout) :: rightv
      real, intent(in) :: nx, ny
      real, dimension(1:gpu_max_nvar) :: tempvect1, tempvect2
      real, dimension(1:gpu_max_species) :: y_s_g
      real :: mp_pinfl, gammal, agrt, u, v, un, mn, rho_l, rho_g, p_l, p_g, u_g, v_g, rg_rmix, rmix, t_g, tv_g
      real :: rg_tr, rg_ev_total, rg_chem, skin1, rg_cv
      integer :: rg_i

      tempvect1 = zero
      tempvect2 = zero
      tempvect1(1:nof_variables) = leftv(1:nof_variables)
      tempvect2(1:nof_variables) = leftv(1:nof_variables)

      call cons2prim(n, tempvect1, mp_pinfl, gammal)
      agrt = sqrt((tempvect1(4) + mp_pinfl) * gammal / tempvect1(1))
      rho_l = tempvect1(1)
      u = tempvect1(2)
      v = tempvect1(3)
      p_l = tempvect1(4)
      un = u * nx + v * ny

      call cons2div(n, tempvect2, mp_pinfl, gammal)

      if (un .le. 0.00) then
         call inflow2d(initcond, pox, poy, rightv)
      else
         mn = un / agrt

         if (mn .ge. 1.00) then
            rightv(1:nof_variables) = leftv(1:nof_variables)
         else
            p_g = max(pres, tolsmall)
            rho_g = rho_l + (p_g - p_l) / ((agrt * agrt) + tolsmall)
            rho_g = max(rho_g, tolsmall)
            u_g = u + (nx * (p_l - p_g)) / ((agrt + tolsmall) * (rho_l + tolsmall))
            v_g = v + (ny * (p_l - p_g)) / ((agrt + tolsmall) * (rho_l + tolsmall))
            y_s_g(1:nof_species) = tempvect1(dimensiona + 4:nof_variables)

            rg_rmix = zero
            do rg_i = 1, nof_species
               rg_rmix = rg_rmix + y_s_g(rg_i) / rg_molm(rg_i)
            end do

            rg_rmix = rgs_ru * rg_rmix
            rmix = rg_rmix
            t_g = p_g / (rho_g * rmix)
            tv_g = tempvect2(5)

            rg_tr = zero
            do rg_i = 1, nof_species
               if (rg_i .le. 3) then
                  rg_cv = (5.00 / 2.00) * rgs_ru / rg_molm(rg_i)
               else
                  rg_cv = (3.00 / 2.00) * rgs_ru / rg_molm(rg_i)
               end if
               rg_tr = rg_tr + y_s_g(rg_i) * rg_cv * t_g
            end do

            rg_ev_total = zero
            do rg_i = 1, 3
               rg_ev_total = rg_ev_total + y_s_g(rg_i) * (rgs_ru / rg_molm(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i) / &
                  & tv_g) - 1.00))
            end do

            rg_chem = zero
            do rg_i = 1, nof_species
               if (rg_hzero(rg_i) .gt. 1.0e-12) then
                  rg_chem = rg_chem + (y_s_g(rg_i) * rg_hzero(rg_i) / rg_molm(rg_i))
               end if
            end do

            skin1 = oo2 * ((u_g**2) + (v_g**2))

            rightv(1) = rho_g
            rightv(2) = rho_g * u_g
            rightv(3) = rho_g * v_g
            rightv(4) = rho_g * (rg_ev_total + rg_tr + rg_chem + skin1)
            rightv(5) = rho_g * rg_ev_total

            do rg_i = 1, nof_species
               rightv(5 + rg_i) = rho_g * y_s_g(rg_i)
            end do
         end if
      end if

   end subroutine boundarys2d_realgas_outflow

   logical function duct_outlet_face_flux_enabled()
! The stationary ideal-gas FV duct uses a characteristic boundary face state.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      duct_outlet_face_flux_enabled = initcond == 4440 .and. realgas == 0 .and. multispecies == 0 &
         & .and. dg == 0 .and. mrf == 0 .and. rframe == 0 .and. lmach == 0
   end function duct_outlet_face_flux_enabled

   subroutine pressure_outlet_face_flux(n, qface, tface, normal, flux, un)
! Flux of the prescribed face state: no Riemann contact diffusion into the reservoir.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: qface(nof_variables), tface(turbulenceequations + passivescalar), normal(3)
      real, intent(out) :: flux(gpu_max_nvar_total), un
      real :: prim(gpu_max_nvar), pinf, gm
      integer :: nd, ie, nt
      nd = dimensiona; ie = nd + 2; nt = turbulenceequations + passivescalar
      prim = 0.0e0; prim(1:nof_variables) = qface
      call cons2prim(n, prim, pinf, gm)
      un = sum(prim(2:nd + 1) * normal(1:nd))
      flux = 0.0e0
      flux(1:nof_variables) = un * qface
      flux(2:nd + 1) = flux(2:nd + 1) + prim(ie) * normal(1:nd)
      flux(ie) = un * (qface(ie) + prim(ie))
      flux(nof_variables + 1:nof_variables + nt) = un * tface
   end subroutine pressure_outlet_face_flux

   subroutine pressure_outlet_face_jacobian(n, qcell, normal, pout, jac, scalar_jac)
! Analytic d F(B(U))/dU for the ideal-gas duct; differentiate the active branch.
! scalar_jac holds d F_transport/d transport with mean flow frozen.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: qcell(nof_variables), normal(3), pout
      real, intent(out) :: jac(gpu_max_nvar, gpu_max_nvar), scalar_jac(gpu_max_extra_transport)
      real :: prim(gpu_max_nvar), qb(gpu_max_nvar), v(3), vb(3), dv(3), dmb(3)
      real :: rho, p, a, pg, pinf, gm, rhob, unb, drho, dp, drhob, dpb, deb, dunb, acoustic, dacoustic, c
      integer :: nd, ie, j, flag, nt
      logical :: backflow
      nd = dimensiona; ie = nd + 2; nt = turbulenceequations + passivescalar
      prim = 0.0e0; prim(1:nof_variables) = qcell
      call cons2prim(n, prim, pinf, gm)
      rho = prim(1); p = prim(ie); a = sqrt(gamma * p / rho)
      v = 0.0e0; v(1:nd) = prim(2:nd + 1)
      call pressure_outlet_state(n, qcell, qb(1:nof_variables), normal, pout, backflow, flag)
      rhob = qb(1); vb = 0.0e0; vb(1:nd) = qb(2:nd + 1) / rhob
      unb = sum(vb * normal); pg = max(pout, tolsmall)
      if (flag == -3) pg = p
      acoustic = (p - pg) / (rho * a)
      c = (gamma - 1.0e0) / gamma + pg / (gamma * p)
      jac = 0.0e0; scalar_jac = 0.0e0
      do j = 1, nof_variables
         drho = 0.0e0; dp = 0.0e0; dv = 0.0e0
         if (j == 1) then
            drho = 1.0e0; dv = -v / rho
            dp = 0.5e0 * (gamma - 1.0e0) * sum(v * v)
         elseif (j <= nd + 1) then
            dv(j - 1) = 1.0e0 / rho; dp = -(gamma - 1.0e0) * v(j - 1)
         elseif (j == ie) then
            dp = gamma - 1.0e0
         end if
         dpb = 0.0e0
         if (flag == -3) then
            drhob = drho; dpb = dp
         else
            dacoustic = (p + pg) * dp / (2.0e0 * p * rho * a) - acoustic * drho / (2.0e0 * rho)
            if (backflow) then
               drhob = 0.0e0; dv = normal * (sum(normal * dv) + dacoustic)
            else
               drhob = c * drho - rho * pg * dp / (gamma * p * p)
               if (rho + (pg - p) / (a * a) <= tolsmall) drhob = 0.0e0
               dv = dv + normal * dacoustic
            end if
         end if
         dmb = vb * drhob + rhob * dv
         deb = dpb / (gamma - 1.0e0) + 0.5e0 * sum(vb * vb) * drhob + rhob * sum(vb * dv)
         dunb = sum(normal * dv)
         jac(1, j) = unb * drhob + rhob * dunb
         jac(2:nd + 1, j) = unb * dmb(1:nd) + qb(2:nd + 1) * dunb + normal(1:nd) * dpb
         jac(ie, j) = unb * (deb + dpb) + (qb(ie) + pg) * dunb
      end do
      if (.not. backflow) then
         scalar_jac(1:nt) = unb
         if (turbulence == 1 .and. turbulencemodel == 1) scalar_jac(1) = unb * rhob / rho
      end if
   end subroutine pressure_outlet_face_jacobian

   subroutine pressure_outlet_state(n, qleft, qright, normal, pout, backflow, regime)
! Common CPU/XPU outlet entry for every equation of state.
! On backflow, incoming temperature/composition use the configured reference gas.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: qleft(nof_variables), normal(3), pout
      real, intent(out) :: qright(nof_variables)
      logical, intent(out) :: backflow
      integer, intent(out), optional :: regime
      real :: prim(gpu_max_nvar), hv(gpu_max_species), vel(3), refvel(3)
      real :: pinf, gm, rho, p, a, un, pg, rgmix, tref, ratio, un_g, refun
      integer :: nd, ie, k
      nd = dimensiona; ie = nd + 2
      prim = 0.0e0; prim(1:nof_variables) = qleft
      pinf = 0.0e0; gm = gamma
      call cons2prim(n, prim, pinf, gm)
      rho = prim(1); p = prim(ie)
      a = sqrt(max(gm * (p + pinf) / rho, tolsmall))
      vel = 0.0e0; vel(1:nd) = prim(2:nd + 1)
      un = sum(vel * normal); backflow = (un .lt. 0.0e0)
      qright = qleft
      if (present(regime)) regime = -2
      if (present(regime)) then
         if (backflow) regime = -1
      end if
      if (un .ge. a) then
         if (present(regime)) regime = -3
         return
      end if
      pg = max(pout, tolsmall)
      un_g = un + (p - pg) / (rho * a)
      if (duct_outlet_face_flux_enabled()) then
         backflow = un_g < 0.0e0
         if (present(regime)) regime = merge(-1, -2, backflow)
      end if
      if (backflow) then
         ! One outgoing acoustic characteristic; reference incoming thermal/composition data.
         refvel = [uvel, vvel, wvel]; if (nd .eq. 2) refvel(3) = 0.0e0
         refun = sum(refvel * normal)
         prim(2:nd + 1) = refvel(1:nd) + (un_g - refun) * normal(1:nd)
         if (realgas .eq. 1) then
            prim(nd + 4:nof_variables) = rg_vf(1:nof_species)
            rgmix = 0.0e0
            do k = 1, nof_species
               rgmix = rgmix + rg_vf(k) * rgs_ru / rg_molm(k)
            end do
            prim(1) = pg / (rgmix * max(rg_ttr, tolsmall))
            call hvib_air5(max(rg_tve, tolsmall), hv)
            prim(nd + 3) = sum(rg_vf(1:nof_species) * hv(1:nof_species))
         else
            tref = pres / max(rres * r_gas, tolsmall)
            if ((initcond .eq. 4440) .and. (total_temperature_inlet .gt. 0.0e0)) tref = total_temperature_inlet
            prim(1) = pg / max(r_gas * tref, tolsmall)
         end if
      else
         prim(1) = max(rho + (pg - p) / (a * a), tolsmall)
         prim(2:nd + 1) = vel(1:nd) + (un_g - un) * normal(1:nd)
      end if
      if (multispecies .eq. 1) then
         ratio = prim(1) / rho
         prim(nd + 3:nd + 2 + nof_species) = prim(nd + 3:nd + 2 + nof_species) * ratio
      end if
      prim(ie) = pg
      call prim2cons(n, prim)
      qright = prim(1:nof_variables)
   end subroutine pressure_outlet_state

   subroutine pressure_outlet_transport(n, qright, transport)
! Incoming transported quantities use the same reference data as the inlet.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: qright(nof_variables)
      real, intent(out) :: transport(turbulenceequations + passivescalar)
      real :: sa_bc
      transport = 0.0e0
      if (turbulence .eq. 1) then
         if (turbulencemodel .eq. 1) then
            transport(1) = visc * turbinit * qright(1)
            if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
               call inflow_4440_sa_value(n, qright, sa_bc)
               transport(1) = sa_bc
            end if
         else if (turbulencemodel .eq. 2) then
            transport(1) = 1.5e0 * i_turb_inlet * ufreestream**2 * qright(1)
            transport(2) = qright(1) * transport(1) / max(10.0e-5 * visc, tolsmall)
         end if
      end if
      if (passivescalar .gt. 0) transport(turbulenceequations + 1:) = qright(1) * rres
   end subroutine pressure_outlet_transport
   subroutine write_inlet_reynolds_4440(n)
! Area average of rho_face * |u_face| * L / mu_ref, with L = 1 metre.
! Evaluate the 4440 inlet BC from the current adjacent cell averages (code 1).
! Called by the master thread on every MPI rank, independently of the controller.
      implicit none
      integer, intent(in) :: n
      integer :: ii, i, l, ibc, output_unit, file_size
      real :: sums_local(4), sums_global(4), area, speed, reynolds_inlet
      real :: qleft(gpu_max_nvar), qin(gpu_max_nvar), normal(3), angle1, angle2
      logical :: exists

      if (initcond .ne. 4440 .or. duct_viscosity_ref .le. 0.0e0) return
      if (realgas .ne. 0 .or. multispecies .ne. 0) return
      if (inlet_reynolds_freq .le. 0) return
      if (mod(it, inlet_reynolds_freq) .ne. 0) return
      sums_local = 0.0e0
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& reduction(+:sums_local) firstprivate(n, nof_bounded, dimensiona, nof_variables) &
!$omp& map(alloc:el_bnd, ielem_ifca, ielem_ibounds, ibound_icode, ielem_surf, u_c_val) &
!$omp& map(alloc:ielem_faceanglex, ielem_faceangley) &
!$omp& private(ii, i, l, ibc, qleft, qin, normal, angle1, angle2, area, speed)
#endif
      do ii = 1, nof_bounded
         i = el_bnd(ii)
         do l = 1, ielem_ifca(i)
            ibc = ielem_ibounds(i, l)
            if (ibc .le. 0) cycle
            if (ibound_icode(ibc) .ne. 1) cycle
            area = ielem_surf(i, l)
            if (area .le. 0.0e0) cycle
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            if (dimensiona .eq. 3) then
               normal = [cos(angle1) * sin(angle2), sin(angle1) * sin(angle2), cos(angle2)]
            else
               normal = [angle1, angle2, 0.0e0]
            end if
            qleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            call inflow_4440_from_left(n, qleft(1:nof_variables), normal, qin(1:nof_variables))
            speed = sqrt(sum((qin(2:dimensiona + 1) / qin(1))**2))
            sums_local(1) = sums_local(1) + area
            sums_local(2) = sums_local(2) + area * qin(1)
            sums_local(3) = sums_local(3) + area * speed
            sums_local(4) = sums_local(4) + area * qin(1) * speed
         end do
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#endif
      call mpi_allreduce(sums_local, sums_global, 4, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      if (sums_global(1) .le. 0.0e0) return
      reynolds_inlet = sums_global(4) / sums_global(1) / duct_viscosity_ref
      if (n .eq. 0) then
         inquire (file='INLET_REYNOLDS.dat', exist=exists, size=file_size)
         open (newunit=output_unit, file='INLET_REYNOLDS.dat', form='formatted', status='unknown', &
            & action='write', position='append')
         if (.not. exists .or. file_size .eq. 0) then
            write (output_unit, '(A)') '# it t Re_inlet_area rho_inlet_area U_inlet_area mu_ref L_ref inlet_area'
         end if
         write (output_unit, '(I10,7(1X,ES20.10))') it, t, reynolds_inlet, sums_global(2) / sums_global(1), &
            & sums_global(3) / sums_global(1), duct_viscosity_ref, 1.0e0, sums_global(1)
         flush (output_unit)
         close (output_unit)
      end if
   end subroutine write_inlet_reynolds_4440

   subroutine update_outlet_mach_controller(n)
      implicit none
      integer, intent(in) :: n
      integer :: ii, i, l, ibc, eff_unit
      real :: mach_sum_l, area_sum_l, p_sum_l, pt_sum_l, pt_mdot_sum_l, mdot_sum_l, backflow_sum_l
      real :: mach_sum_g, area_sum_g, p_sum_g, pt_sum_g, pt_mdot_sum_g, mdot_sum_g, backflow_sum_g
      real :: rho, uu, vv, ww, e, vel2, p, a, mach, area, err, limited_err, factor
      real :: angle1, angle2, nx, ny, nz, un, mdot, pt, pt_area_out, pt_mass_out
      real :: eta_area, eta_mass, loss_area, loss_mass, eta_area_percent, eta_mass_percent
      logical :: heref, controller_active

      if ((mach_outlet_target .le. 0.00) .or. (mach_outlet_update_freq .le. 0)) return
      if ((realgas .ne. 0) .or. (multispecies .ne. 0)) return

      mach_sum_l = 0.00
      area_sum_l = 0.00
      p_sum_l = 0.00
      pt_sum_l = 0.00
      pt_mdot_sum_l = 0.00
      mdot_sum_l = 0.00
      backflow_sum_l = 0.00

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& reduction(+:mach_sum_l, area_sum_l, p_sum_l, pt_sum_l, pt_mdot_sum_l, mdot_sum_l, backflow_sum_l) &
!$omp& firstprivate(nof_bounded, dimensiona, gamma, tolsmall, oo2) &
!$omp& map(alloc: el_bnd, ielem_ifca, ielem_ibounds, ibound_icode, ielem_surf, u_c_val) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley) &
!$omp& private(ii, i, l, ibc, rho, uu, vv, ww, e, vel2, p, a, mach, area, angle1, angle2, nx, ny, nz, un, mdot, pt)
#endif
      do ii = 1, nof_bounded
         i = el_bnd(ii)
         do l = 1, ielem_ifca(i)
            ibc = ielem_ibounds(i, l)
            if (ibc .gt. 0) then
               if (ibound_icode(ibc) .eq. 9) then
                  rho = max(u_c_val(i, 1, 1), tolsmall)
                  uu = u_c_val(i, 1, 2) / rho
                  vv = u_c_val(i, 1, 3) / rho
                  if (dimensiona .eq. 3) then
                     ww = u_c_val(i, 1, 4) / rho
                     e = u_c_val(i, 1, 5)
                  else
                     ww = 0.00
                     e = u_c_val(i, 1, 4)
                  end if
                  vel2 = (uu * uu) + (vv * vv) + (ww * ww)
                  p = (gamma - 1.00) * (e - rho * oo2 * vel2)
                  if (p .gt. tolsmall) then
                     a = sqrt(gamma * p / rho)
                     mach = sqrt(max(vel2, 0.00)) / (a + tolsmall)
                     area = max(ielem_surf(i, l), 0.00)
                     pt = p * (1.00 + (0.50 * (gamma - 1.00) * mach * mach))**(gamma / (gamma - 1.00))
                     mach_sum_l = mach_sum_l + (mach * area)
                     area_sum_l = area_sum_l + area
                     p_sum_l = p_sum_l + (p * area)
                     pt_sum_l = pt_sum_l + (pt * area)
                     angle1 = ielem_faceanglex(i, l)
                     angle2 = ielem_faceangley(i, l)
                     if (dimensiona .eq. 3) then
                        nx = cos(angle1) * sin(angle2)
                        ny = sin(angle1) * sin(angle2)
                        nz = cos(angle2)
                     else
                        nx = angle1
                        ny = angle2
                        nz = 0.00
                     end if
                     un = (uu * nx) + (vv * ny) + (ww * nz)
                     mdot = rho * un * area
                     if (mdot .gt. 0.00) then
                        mdot_sum_l = mdot_sum_l + mdot
                        pt_mdot_sum_l = pt_mdot_sum_l + (pt * mdot)
                     else
                        backflow_sum_l = backflow_sum_l - mdot
                     end if
                  end if
               end if
            end if
         end do
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#endif

      mach_sum_g = 0.00
      area_sum_g = 0.00
      p_sum_g = 0.00
      pt_sum_g = 0.00
      pt_mdot_sum_g = 0.00
      mdot_sum_g = 0.00
      backflow_sum_g = 0.00
      call mpi_allreduce(mach_sum_l, mach_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(area_sum_l, area_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(p_sum_l, p_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(pt_sum_l, pt_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(pt_mdot_sum_l, pt_mdot_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(mdot_sum_l, mdot_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
      call mpi_allreduce(backflow_sum_l, backflow_sum_g, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)

      if (area_sum_g .le. tolsmall) return

      mach_outlet_average = mach_sum_g / area_sum_g
      if (mach_outlet_filter_ready .eq. 0) then
         mach_outlet_filtered = mach_outlet_average
         mach_outlet_filter_ready = 1
      else
         mach_outlet_filtered = ((1.00 - mach_outlet_filter_alpha) * mach_outlet_filtered) + (mach_outlet_filter_alpha * &
            & mach_outlet_average)
      end if
      controller_active = (it .ge. mach_outlet_start_iter)
      pt_area_out = pt_sum_g / (area_sum_g + tolsmall)
      pt_mass_out = -1.00
      if (mdot_sum_g .gt. tolsmall) pt_mass_out = pt_mdot_sum_g / mdot_sum_g
      eta_area = -1.00
      eta_mass = -1.00
      loss_area = -1.00
      loss_mass = -1.00
      eta_area_percent = -1.00
      eta_mass_percent = -1.00
      if (total_pressure_inlet .gt. tolsmall) then
         eta_area = pt_area_out / total_pressure_inlet
         loss_area = 1.00 - eta_area
         eta_area_percent = 100.00 * eta_area
         if (pt_mass_out .gt. tolsmall) then
            eta_mass = pt_mass_out / total_pressure_inlet
            loss_mass = 1.00 - eta_mass
            eta_mass_percent = 100.00 * eta_mass
         end if
      end if
      if (press_outlet .le. tolsmall) press_outlet = max(p_sum_g / (area_sum_g + tolsmall), tolsmall)

      if (controller_active) then
         err = mach_outlet_filtered - mach_outlet_target
         limited_err = max(min(mach_outlet_relax * err, 0.050), -0.050)
         factor = exp(limited_err)
         press_outlet = max(press_outlet * factor, tolsmall)
      end if

#ifdef xpu
!$omp target update to(press_outlet)
#endif

      if (n .eq. 0) then
         inquire (file='OUTLET_MACH.dat', exist=heref)
         open (91, file='OUTLET_MACH.dat', form='formatted', status='unknown', action='write', position='append')
         if (.not. heref) write (91, '(A)') '# it t mach_area target_mach press_outlet outlet_area mach_filtered controller_active'
         write (91, '(I10,6(1X,ES20.10),1X,I4)') it, t, mach_outlet_average, mach_outlet_target, press_outlet, area_sum_g, &
            & mach_outlet_filtered, merge(1, 0, controller_active)
         flush (91)
         close (91)
         inquire (file='DUCT_EFFICIENCY.dat', exist=heref)
         open (newunit=eff_unit, file='DUCT_EFFICIENCY.dat', form='formatted', status='unknown', action='write', position='append')
         if (.not. heref) write (eff_unit, '(A)') '# it t pt_in pt_out_area pt_out_mass eta_area eta_mass '// &
            & 'loss_area loss_mass mach_area press_outlet outlet_area mdot_out mdot_backflow eta_area_percent eta_mass_percent '// &
            & 'mach_filtered controller_active'
         write (eff_unit, '(I10,16(1X,ES20.10),1X,I4)') it, t, total_pressure_inlet, pt_area_out, pt_mass_out, &
            & eta_area, eta_mass, loss_area, loss_mass, mach_outlet_average, press_outlet, area_sum_g, mdot_sum_g, backflow_sum_g, &
            & eta_area_percent, eta_mass_percent, &
            & mach_outlet_filtered, merge(1, 0, controller_active)
         flush (eff_unit)
         close (eff_unit)
      end if

   end subroutine update_outlet_mach_controller

   subroutine inflow_4440_from_left(n, leftv, normal, inflow_out)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(in) :: leftv
      real, intent(in) :: normal(3)
      real, dimension(1:nof_variables), intent(out) :: inflow_out
      real :: rho_l, ke, p_l, mach
      integer :: k
      real :: pt0, t0, trat, temp, r, p, a, vel, rgas_local, pratio, nmag

      rho_l = max(leftv(1), tolsmall)
      ke = 0.0e0
      do k = 2, dimensiona + 1
         ke = ke + (leftv(k) / rho_l)**2
      end do
      p_l = (gamma - 1.00) * (leftv(dimensiona + 2) - rho_l * oo2 * ke)
      if (p_l .le. tolsmall) p_l = max(pres, tolsmall)

      pt0 = max(total_pressure_inlet, tolsmall)
      t0 = max(total_temperature_inlet, tolsmall)
      rgas_local = max(r_gas, tolsmall)
      p = max(min(p_l, pt0), tolsmall)
      pratio = max(pt0 / (p + tolsmall), 1.00)
      mach = sqrt(max((pratio**((gamma - 1.00) / gamma) - 1.00) * (2.00 / (gamma - 1.00)), 0.00))
      mach = min(max(mach, 0.00), 0.9990)
      trat = 1.00 + 0.50 * (gamma - 1.00) * mach * mach
      temp = t0 / trat
      p = pt0 / (trat**(gamma / (gamma - 1.00)))
      r = max(p / (rgas_local * temp), tolsmall)
      a = sqrt(gamma * rgas_local * temp)
      vel = mach * a
      nmag = max(sqrt(sum(normal(1:dimensiona)**2)), tolsmall)
      inflow_out = 0.0e0
      inflow_out(1) = r
      inflow_out(2:dimensiona + 1) = -vel * normal(1:dimensiona) / nmag
      inflow_out(dimensiona + 2) = p
      call prim2cons(n, inflow_out)

   end subroutine inflow_4440_from_left


   subroutine inflow_4440_sa_value(n, inflow_state, sa_value)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(in) :: inflow_state
      real, intent(out) :: sa_value
      real, dimension(1:4) :: viscl, laml

      call get_visc_conduct(n, inflow_state, inflow_state, viscl, laml)
      sa_value = max(viscl(1) * turbinit, tolsmall)

   end subroutine inflow_4440_sa_value

   subroutine boundarys(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, cturbl, cturbr, &
      & cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
      implicit none
!> @brief
!> this subroutine applies the boundary condition to each bounded cell
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, b_code, iconsidered, facex
      real, dimension(1:nof_variables), intent(inout) :: leftv, rightv
      integer, intent(inout) :: ibfc
      real, dimension(1:nof_variables), intent(in) :: srf_speedrot, srf_speed
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      real, intent(in) :: angle1, angle2, nx, ny, nz
      real, dimension(1:turbulenceequations + passivescalar), intent(inout) :: cturbl, cturbr
      real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_nvar) :: subson1, subson2, subson3
      real :: mp_pinfl, mp_pinfr, gammal, gammar, u, v, w, sa_bc
      real :: sps, skins, ikins, vel, vnb, theeta, reeta
      real :: intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx, vhx, amp, dvel, rgg, tt1
      real :: mn, un, agrt, rho_g, rg_rmix, rmix, t_g, tv_g, rg_tr, rg_ev_total, rg_chem
      real, dimension(1:gpu_max_passive) :: pass_tmp
      integer :: rg_i

      logical :: outlet_backflow

! Define the flag before selecting the boundary condition.
      ibfc = 0
      rightv(1:nof_variables) = leftv(1:nof_variables)
      cturbr = cturbl
      if (b_code .eq. 99) error stop 'Boundary code 99: legacy bleed model is incomplete'
      if (b_code .eq. 9) then
         call pressure_outlet_state(n, leftv, rightv, [nx, ny, nz], press_outlet, outlet_backflow, ibfc)
         cturbr = cturbl
         ! Preserve the kinematic SA variable when the pressure outlet changes density.
         if ((.not. outlet_backflow) .and. (turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
            cturbr(1) = cturbl(1) * rightv(1) / leftv(1)
         end if
         if (outlet_backflow) call pressure_outlet_transport(n, rightv, cturbr)
         return
      end if
      select case (b_code)

      case (1) !inflow subsonic or supersonic will be chosen based on mach number
         if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
            call inflow_4440_from_left(n, leftv, [nx, ny, nz], rightv)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) cturbr(:) = cturbl(:)
            if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
               call inflow_4440_sa_value(n, rightv, sa_bc)
               cturbr(1) = sa_bc
            end if
         else
            if (boundtype .eq. 0) then !supersonic

               call inflow(initcond, pox, poy, poz, rightv)

               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                  cturbr(1) = visc * turbinit * rightv(1)
               end if

            else !subsonic
               call inflow(initcond, pox, poy, poz, rightv)

               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)
               if (realgas .eq. 1) then
                  sps = sqrt((gammal * subson2(5)) / (subson2(1) + tolsmall))
               else
                  sps = sqrt((gamma * subson2(5)) / (subson2(1) + tolsmall))
               end if
               vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
               call prim2cons2(n, leftv, rightv)

               if (vel / (sps + tolsmall) .gt. 1.00) then !supersonic

                  call inflow(initcond, pox, poy, poz, rightv)

                  if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                     cturbr(1) = visc * turbinit * rightv(1)
                  end if

               else !subsonic

                  subson3(5) = 0.5 * ((subson1(5)) + (subson2(5)) - (subson2(1) * sps * ((nx * (subson1(2) - subson2(2))) + (ny * &
                     & (subson1(3) - subson2(3))) &
                     & + (nz * (subson1(4) - subson2(4))))))
                  subson3(1) = subson1(1) + (subson3(5) - subson1(5)) / (sps**2)
                  subson3(2) = subson1(2) - (nx * (subson1(5) - subson3(5))) / (sps * subson2(1))
                  subson3(3) = subson1(3) - (ny * (subson1(5) - subson3(5))) / (sps * subson2(1))
                  subson3(4) = subson1(4) - (nz * (subson1(5) - subson3(5))) / (sps * subson2(1))

                  if (realgas .eq. 1) then
                     rightv(1:nof_variables) = subson1(1:nof_variables)
                     rightv(1) = subson3(1)
                     rightv(2) = subson3(2)
                     rightv(3) = subson3(3)
                     rightv(4) = subson3(4)
                     rightv(5) = subson3(5)
                     call prim2cons(n, rightv)
                  else
                     rightv(1) = subson3(1)
                     rightv(2) = subson3(2) * subson3(1)
                     rightv(3) = subson3(3) * subson3(1)
                     rightv(4) = subson3(4) * subson3(1)
                     skins = oo2 * ((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
                     ikins = subson3(5) / ((gamma - 1.00) * (subson3(1)))
                     rightv(5) = (subson3(1) * (ikins)) + (subson3(1) * skins)
                  end if

               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                     cturbr(1) = visc * turbinit * rightv(1)
                  end if
                  if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 2)) then
                     cturbr(1) = (1.50 * i_turb_inlet * (ufreestream**2)) * rightv(1) !k initialization
                     cturbr(2) = rightv(1) * cturbr(1) / (10.0e-5 * visc) !omega initialization
                     if (rec_mrf(iconsidered) .eq. 1) then
                        cturbr(1) = (1.50 * i_turb_inlet * (kinit_srf**2)) * rightv(1) !k initialization
                        cturbr(2) = rightv(1) * cturbr(1) / (10.0e-5 * visc) !omega initialization
                     end if
                  end if

                  if (passivescalar .gt. 0) then
                     call pass_inlet(initcond, pox, poy, poz, pass_tmp)
                     do rg_i = 1, passivescalar
                        cturbr(turbulenceequations + rg_i) = pass_tmp(rg_i) * rightv(1)
                     end do
                  end if
               end if

            end if
         end if

      case (2) !outflow subsonic or supersonic will be chosen based on mach number
         if (boundtype .eq. 0) then
            rightv(1:nof_variables) = leftv(1:nof_variables)

         else

            if (realgas .eq. 1) then
               call boundarys_realgas_outflow(n, initcond, pox, poy, poz, leftv, rightv, nx, ny, nz)
            else !realgas

               ! ------------------------------------------------------------
               ! Build a reference/outlet state (used for prescribed pressure,
               ! and for backflow if it happens)
               ! ------------------------------------------------------------
               call outflow(initcond, pox, poy, poz, rightv)

               ! Convert interior/exterior to primitive for BC logic
               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson2(1:nof_variables) = leftv(1:nof_variables) ! interior primitive
               subson1(1:nof_variables) = rightv(1:nof_variables) ! "target" primitive (from outflow)

               ! Speed of sound from interior
               sps = sqrt((gamma * subson2(5)) / (subson2(1) + tolsmall))

               ! Normal velocity from interior
               vel = subson2(2) * nx + subson2(3) * ny + subson2(4) * nz

               ! ------------------------------------------------------------
               ! Decide boundary type
               ! ------------------------------------------------------------

               ! 1) Backflow / inflow at outlet: treat as inflow (use reference state)
               if (vel .le. 0.00) then

                  ! Keep the primitive outlet state from outflow() (subson1)
                  subson3(1:nof_variables) = subson1(1:nof_variables)

               else

                  ! 2) Outflow: check normal Mach number
                  if (vel / (sps + tolsmall) .gt. 1.00) then
                     ! Supersonic outflow: extrapolate (copy interior primitive)
                     subson3(1:nof_variables) = subson2(1:nof_variables)

                  else
                     ! Subsonic outflow: impose pressure (from reference state),
                     ! correct rho and velocity using linearized characteristics
                     subson3(1:nof_variables) = subson2(1:nof_variables)

                     subson3(5) = subson1(5) ! prescribed outlet static pressure (e.g. p_inf)

                     subson3(1) = subson2(1) + (subson3(5) - subson2(5)) / (sps * sps + tolsmall)

                     subson3(2) = subson2(2) + (nx * (subson2(5) - subson3(5))) / ((sps + tolsmall) * (subson2(1) + tolsmall))
                     subson3(3) = subson2(3) + (ny * (subson2(5) - subson3(5))) / ((sps + tolsmall) * (subson2(1) + tolsmall))
                     subson3(4) = subson2(4) + (nz * (subson2(5) - subson3(5))) / ((sps + tolsmall) * (subson2(1) + tolsmall))

                  end if
               end if

               ! ------------------------------------------------------------
               ! Convert primitive subson3 -> conservative rightv
               ! ------------------------------------------------------------
               rightv(1) = subson3(1)
               rightv(2) = subson3(2) * subson3(1)
               rightv(3) = subson3(3) * subson3(1)
               rightv(4) = subson3(4) * subson3(1)

               skins = oo2 * ((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
               ikins = subson3(5) / ((gamma - 1.00) * (subson3(1) + tolsmall))

               rightv(5) = (subson3(1) * ikins) + (subson3(1) * skins)

            end if
         end if

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbr(:) = cturbl(:)

         end if

      case (99) !bleed boundary

         call bleed3d(iconsidered, facex, pox, poy, poz, rightv)

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbr(:) = cturbl(:)

         end if

      case (3) !symmetry

         call rotatef(n, cleft_rot, leftv, angle1, angle2)
         if (rec_mrf(iconsidered) .eq. 1) then
            cright_rot(1) = cleft_rot(1)
            cright_rot(2) = -(cleft_rot(2)) + 2.00 * cleft_rot(1) * srf_speedrot(2)
            cright_rot(3) = cleft_rot(3)
            cright_rot(4) = cleft_rot(4)
            cright_rot(5) = cleft_rot(5) + 2.00 * cleft_rot(1) * (srf_speedrot(2)**2) - 2.00 * cleft_rot(2) * srf_speedrot(2)
         else
            cright_rot(1) = cleft_rot(1)
            cright_rot(2) = -cleft_rot(2)
            cright_rot(3) = cleft_rot(3)
            cright_rot(4) = cleft_rot(4)
            cright_rot(5) = cleft_rot(5)

            if ((multispecies .eq. 1) .or. (realgas .eq. 1)) then
               cright_rot(6:nof_variables) = cleft_rot(6:nof_variables)
            end if

         end if
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            cturbr(:) = cturbl(:)

            if (passivescalar .gt. 0) then
               cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                  & cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

            end if
         end if

         call rotateb(n, rightv, cright_rot, angle1, angle2)

      case (4) !wall

         if (itestcase .eq. 3) then

            call rotatef(n, cleft_rot, leftv, angle1, angle2)
            if (rec_mrf(iconsidered) .eq. 1) then
               cright_rot(1) = cleft_rot(1)
               cright_rot(2) = -(cleft_rot(2)) + 2.00 * cleft_rot(1) * srf_speedrot(2)
               cright_rot(3) = cleft_rot(3)
               cright_rot(4) = cleft_rot(4)
               cright_rot(5) = cleft_rot(5) + cleft_rot(1) * (srf_speedrot(2)**2) * 2.00 - 2.00 * cleft_rot(2) * srf_speedrot(2)
            else
               cright_rot(:) = cleft_rot(:)
               cright_rot(2) = -cleft_rot(2)
            end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               cturbr(:) = cturbl(:)

               if (passivescalar .gt. 0) then
                  cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                     & cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

               end if
            end if

            call rotateb(n, rightv, cright_rot, angle1, angle2)

         else
            if (rec_mrf(iconsidered) .eq. 1) then
               rightv(1) = leftv(1)
               rightv(2) = -leftv(2) + 2.00 * leftv(1) * srf_speed(2)
               rightv(3) = -leftv(3) + 2.00 * leftv(1) * srf_speed(3)
               rightv(4) = -leftv(4) + 2.00 * leftv(1) * srf_speed(4)
               rightv(5) = leftv(5) + 2.00 * leftv(1) * (srf_speed(2)**2 + srf_speed(3)**2 + srf_speed(4)**2) &
                  & - 2.00 * (leftv(2) * srf_speed(2) + leftv(3) * srf_speed(3) + leftv(4) * srf_speed(4))

            else
               rightv(1:nof_variables) = leftv(1:nof_variables)
               rightv(2) = -leftv(2)
               rightv(3) = -leftv(3)
               rightv(4) = -leftv(4)

               if (realgas .eq. 1) then
                  if (catalytic_wall .eq. 1) then
                     rightv(dimensiona + 4:nof_variables) = catalytic_con(1:nof_species) * leftv(1)

                  end if
               end if

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if ((turbulence .ne. 1) .or. (turbulencemodel .ne. 2)) then
                  cturbr(:) = -cturbl(:)

                  if (passivescalar .gt. 0) then
                     cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                        & -cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

                  end if
               else
                  cturbr(1) = -cturbl(1)
                  cturbr(2) = 60.00 * visc / (beta_i1 * (ielem_walldist(iconsidered)**2))

                  if (passivescalar .gt. 0) then
                     cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                        & -cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

                  end if

               end if
            end if

         end if

      case (6) !farfield inflow or outflow, subsonic or supersonic will be chosen based on mach number
         call rotatef(n, cleft_rot, leftv, angle1, angle2)
         vnb = cleft_rot(2) / cleft_rot(1)

         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         subson1(1:nof_variables) = rightv(1:nof_variables)
         subson2(1:nof_variables) = leftv(1:nof_variables)
         if (realgas .eq. 1) then
            sps = sqrt((gammal * subson2(5)) / (subson2(1) + tolsmall))
         else
            sps = sqrt((gamma * subson2(5)) / (subson2(1) + tolsmall))
         end if
         vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)

         call prim2cons2(n, leftv, rightv)

         if (vnb .le. 0.00) then !inflow
            ibfc = -1

            if ((abs(vnb)) .ge. sps) then
               !supersonic
               if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
                  call inflow_4440_from_left(n, leftv, [nx, ny, nz], rightv)
               else
                  call inflow(initcond, pox, poy, poz, rightv)
               end if

            else
               !subsonic

               if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
                  call inflow_4440_from_left(n, leftv, [nx, ny, nz], rightv)
               else
                  call inflow(initcond, pox, poy, poz, rightv)
               end if

               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)
               if (realgas .eq. 1) then
                  sps = sqrt((gammal * subson2(5)) / (subson2(1) + tolsmall))
               else
                  sps = sqrt((gamma * subson2(5)) / (subson2(1) + tolsmall))
               end if
               vel = sqrt(subson2(2)**2 + subson2(3)**2 + subson2(4)**2)
               call prim2cons2(n, leftv, rightv)

               subson3(5) = 0.5 * ((subson1(5)) + (subson2(5)) - (subson2(1) * sps * ((nx * (subson1(2) - subson2(2))) + (ny * &
                  & (subson1(3) - subson2(3))) &
                  & + (nz * (subson1(4) - subson2(4))))))
               subson3(1) = subson1(1) + (subson3(5) - subson1(5)) / (sps**2)
               subson3(2) = subson1(2) - (nx * (subson1(5) - subson3(5))) / (sps * subson2(1))
               subson3(3) = subson1(3) - (ny * (subson1(5) - subson3(5))) / (sps * subson2(1))
               subson3(4) = subson1(4) - (nz * (subson1(5) - subson3(5))) / (sps * subson2(1))

               if (realgas .eq. 1) then
                  rightv(1:nof_variables) = subson1(1:nof_variables)
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2)
                  rightv(3) = subson3(3)
                  rightv(4) = subson3(4)
                  rightv(5) = subson3(5)
                  call prim2cons(n, rightv)
               else
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2) * subson3(1)
                  rightv(3) = subson3(3) * subson3(1)
                  rightv(4) = subson3(4) * subson3(1)
                  skins = oo2 * ((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
                  ikins = subson3(5) / ((gamma - 1.00) * (subson3(1)))
                  rightv(5) = (subson3(1) * (ikins)) + (subson3(1) * skins)
               end if

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                  cturbr(1) = visc * turbinit * rightv(1)
               end if
               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 2)) then
                  cturbr(1) = (1.50 * i_turb_inlet * (ufreestream**2)) * rightv(1) !k initialization
                  cturbr(2) = rightv(1) * cturbr(1) / (10.0e-5 * visc) !omega initialization
                  if (rec_mrf(iconsidered) .eq. 1) then
                     cturbr(1) = (1.50 * i_turb_inlet * (kinit_srf**2)) * rightv(1) !k initialization
                     cturbr(2) = rightv(1) * cturbr(1) / (10.0e-5 * visc) !omega initialization
                  end if
               end if

               if (passivescalar .gt. 0) then
                  call pass_inlet(initcond, pox, poy, poz, pass_tmp)
                  do rg_i = 1, passivescalar
                     cturbr(turbulenceequations + rg_i) = pass_tmp(rg_i) * rightv(1)
                  end do
               end if
            end if

         else

            !outflow

            ibfc = -2
            if ((abs(vnb)) .ge. sps) then

               rightv(1:nof_variables) = leftv(1:nof_variables)

            else

               call outflow(initcond, pox, poy, poz, rightv)

               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)

               call prim2cons2(n, leftv, rightv)

               subson3(5) = subson1(5)
               subson3(1) = subson2(1) + (subson3(5) - subson2(5)) / (sps**2)
               subson3(2) = subson2(2) + (nx * (subson2(5) - subson3(5))) / (sps * subson2(1))
               subson3(3) = subson2(3) + (ny * (subson2(5) - subson3(5))) / (sps * subson2(1))
               subson3(4) = subson2(4) + (nz * (subson2(5) - subson3(5))) / (sps * subson2(1))
!
               if (realgas .eq. 1) then
                  rightv(1:nof_variables) = subson2(1:nof_variables)
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2)
                  rightv(3) = subson3(3)
                  rightv(4) = subson3(4)
                  rightv(5) = subson3(5)
                  call prim2cons(n, rightv)
               else
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2) * subson3(1)
                  rightv(3) = subson3(3) * subson3(1)
                  rightv(4) = subson3(4) * subson3(1)
                  skins = oo2 * ((subson3(2)**2) + (subson3(3)**2) + (subson3(4)**2))
                  ikins = subson3(5) / ((gamma - 1.00) * (subson3(1)))
                  rightv(5) = (subson3(1) * (ikins)) + (subson3(1) * skins)
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cturbr(:) = cturbl(:)

               end if
            end if

         end if

!         call rotatef(n,cright_rot,rightv,angle1,angle2)

      end select

   end subroutine boundarys

   subroutine boundarys2d(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, cturbl, cturbr, &
      & cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
      implicit none
!> @brief
!> this subroutine applies the boundary condition to each bounded cell
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, b_code, iconsidered, facex
      integer, intent(inout) :: ibfc
      real, dimension(1:nof_variables), intent(inout) :: leftv, rightv
      real, dimension(1:nof_variables), intent(in) :: srf_speedrot, srf_speed
      real, dimension(1:dimensiona), intent(in) :: pox, poy, poz
      real, intent(in) :: angle1, angle2, nx, ny, nz
      real, dimension(1:turbulenceequations + passivescalar), intent(inout) :: cturbl, cturbr
      real, dimension(1:nof_variables + turbulenceequations + passivescalar), intent(inout) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_nvar) :: subson1, subson2, subson3
      real :: mp_pinfl, mp_pinfr, gammal, gammar, u, v, sa_bc
      real :: sps, skins, ikins, vel, vnb, theeta, reeta
      real :: intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx, vhx, amp, dvel, rgg, tt1
      real :: mn, un, agrt, rho_g, rg_rmix, rmix, t_g, tv_g, rg_tr, rg_ev_total, rg_chem
      real, dimension(1:gpu_max_passive) :: pass_tmp
      integer :: rg_i

      logical :: outlet_backflow

! Define the flag before selecting the boundary condition.
      ibfc = 0
      rightv(1:nof_variables) = leftv(1:nof_variables)
      cturbr = cturbl
      if (b_code .eq. 99) error stop 'Boundary code 99: legacy bleed model is incomplete'
      if (b_code .eq. 9) then
         call pressure_outlet_state(n, leftv, rightv, [nx, ny, 0.0e0], press_outlet, outlet_backflow, ibfc)
         cturbr = cturbl
         ! Preserve the kinematic SA variable when the pressure outlet changes density.
         if ((.not. outlet_backflow) .and. (turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
            cturbr(1) = cturbl(1) * rightv(1) / leftv(1)
         end if
         if (outlet_backflow) call pressure_outlet_transport(n, rightv, cturbr)
         return
      end if
      select case (b_code)

      case (1) !inflow subsonic or supersonic will be chosen based on mach number
         if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
            call inflow_4440_from_left(n, leftv, [nx, ny, 0.0e0], rightv)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) cturbr(:) = cturbl(:)
            if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
               call inflow_4440_sa_value(n, rightv, sa_bc)
               cturbr(1) = sa_bc
            end if
         else
            if (boundtype .eq. 0) then !supersonic

               call inflow2d(initcond, pox, poy, rightv)

               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                  cturbr(1) = visc * turbinit * rightv(1)
               end if

            else !subsonic
               call inflow2d(initcond, pox, poy, rightv)

               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)
               if (realgas .eq. 1) then
                  sps = sqrt((gammal * subson2(4)) / (subson2(1) + tolsmall))
               else
                  sps = sqrt((gamma * subson2(4)) / (subson2(1) + tolsmall))
               end if
               vel = sqrt(subson2(2)**2 + subson2(3)**2)

               call prim2cons2(n, leftv, rightv)

               if (vel / (sps + tolsmall) .gt. 1.00) then !supersonic

                  call inflow2d(initcond, pox, poy, rightv)

               else !subsonic

                  subson3(4) = 0.5 * ((subson1(4)) + (subson2(4)) - (subson2(1) * sps * ((nx * (subson1(2) - subson2(2))) + (ny * &
                     & (subson1(3) - subson2(3))) &
                     & )))
                  subson3(1) = subson1(1) + (subson3(4) - subson1(4)) / (sps**2)
                  subson3(2) = subson1(2) - (nx * (subson1(4) - subson3(4))) / (sps * subson2(1))
                  subson3(3) = subson1(3) - (ny * (subson1(4) - subson3(4))) / (sps * subson2(1))

                  if (realgas .eq. 1) then
                     rightv(1:nof_variables) = subson1(1:nof_variables)
                     rightv(1) = subson3(1)
                     rightv(2) = subson3(2)
                     rightv(3) = subson3(3)
                     rightv(4) = subson3(4)
                     call prim2cons(n, rightv)
                  else
                     rightv(1) = subson3(1)
                     rightv(2) = subson3(2) * subson3(1)
                     rightv(3) = subson3(3) * subson3(1)

                     skins = oo2 * ((subson3(2)**2) + (subson3(3)**2))
                     ikins = subson3(4) / ((gamma - 1.00) * (subson3(1)))
                     rightv(4) = (subson3(1) * (ikins)) + (subson3(1) * skins)
                  end if

               end if

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                  cturbr(1) = visc * turbinit
               end if
               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 2)) then
                  cturbr(1) = (1.50 * i_turb_inlet * (ufreestream**2)) * rightv(1) !k initialization
                  cturbr(2) = rightv(1) * cturbr(1) / (10.0e-5 * visc) !omega initialization

               end if

               if (passivescalar .gt. 0) then
                  call pass_inlet2d(initcond, pox, poy, pass_tmp)
                  do rg_i = 1, passivescalar
                     cturbr(turbulenceequations + rg_i) = pass_tmp(rg_i) * rightv(1)
                  end do
               end if
            end if

         end if
      case (2) !outflow subsonic or supersonic will be chosen based on mach number
         if (boundtype .eq. 0) then
            rightv(1:nof_variables) = leftv(1:nof_variables)

         else

            if (realgas .eq. 1) then
               call boundarys2d_realgas_outflow(n, initcond, pox, poy, leftv, rightv, nx, ny)
            else
               call outflow2d(initcond, pox, poy, rightv)
               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)

               sps = sqrt((gamma * subson2(4)) / (subson2(1)))
               vel = subson2(2) * nx + subson2(3) * ny
!                         sqrt(subson2(2)**2+subson2(3)**2)

               sps = sqrt((gamma * subson2(4)) / (subson2(1)))

               call prim2cons2(n, leftv, rightv)

               if (vel / (sps + tolsmall) .gt. 1.00) then !supersonic
                  rightv(1:nof_variables) = leftv(1:nof_variables)

               else
                  subson3(4) = subson1(4)
                  subson3(1) = subson2(1) + (subson3(4) - subson2(4)) / (sps**2)
                  subson3(2) = subson2(2) + (nx * (subson2(4) - subson3(4))) / (sps * subson2(1))
                  subson3(3) = subson2(3) + (ny * (subson2(4) - subson3(4))) / (sps * subson2(1))

                  !
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2) * subson3(1)
                  rightv(3) = subson3(3) * subson3(1)

                  skins = oo2 * ((subson3(2)**2) + (subson3(3)**2))
                  ikins = subson3(4) / ((gamma - 1.00) * (subson3(1)))
                  rightv(4) = (subson3(1) * (ikins)) + (subson3(1) * skins)

               end if

            end if
         end if

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbr(:) = cturbl(:)

         end if

      case (99) !bleed

         call bleed2d(iconsidered, facex, pox, poy, rightv)

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbr(:) = cturbl(:)

         end if

      case (3) !symmetry

         if ((initcond .eq. 102) .or. (initcond .eq. 30) .or. (initcond .eq. 222)) then !shock density interaction
            if ((initcond .eq. 102)) then !shock density interaction
               if (pox(1) .lt. ((1.00 / 6.00) + ((1.00 + 20.00 * t) / (sqrt(3.00))))) then
                  r1 = 8.00
                  u1 = 8.25 * cos(pi / 6.00)
                  v1 = -8.25 * sin(pi / 6.00)
                  p1 = 116.5
               else
                  r1 = 1.40
                  u1 = zero
                  v1 = zero
                  p1 = 1.00
               end if
               skin1 = (oo2) * ((u1**2) + (v1**2))
               !internal energy
               ie1 = ((p1) / ((gamma - 1.00) * r1))
               !total energy
               e1 = (p1 / (gamma - 1)) + (r1 * skin1)
               !vector of conserved variables now
               rightv(1) = r1
               rightv(2) = r1 * u1
               rightv(3) = r1 * v1
               rightv(4) = e1
            end if

            if ((initcond .eq. 222)) then
               if (sqrt((pox(1)**2) + (poy(1)**2)) .lt. (t / 3.00)) then
                  r1 = 16.00
                  u1 = 0.0
                  v1 = 0.0
                  p1 = 16.00 / 3.00
               else
                  r1 = 1.00 + (t / sqrt((pox(1)**2) + (poy(1)**2)))
                  reeta = -1
                  theeta = atan(poy(1) / pox(1))
                  u1 = reeta * cos(theeta)
                  v1 = reeta * sin(theeta)
                  p1 = 1.0e-6
               end if
               skin1 = (oo2) * ((u1**2) + (v1**2))
               !internal energy
               ie1 = ((p1) / ((gamma - 1.00) * r1))
               !total energy
               e1 = (p1 / (gamma - 1)) + (r1 * skin1)
               !vector of conserved variables now
               rightv(1) = r1
               rightv(2) = r1 * u1
               rightv(3) = r1 * v1
               rightv(4) = e1

            end if

            if ((initcond .eq. 30)) then
               if (pox(1) .le. zero) then
                  if (poy(1) .le. zero) then
                     r1 = 0.138
                     u1 = 1.206
                     v1 = 1.206
                     p1 = 0.029
                  end if
                  if (poy(1) .gt. zero) then
                     r1 = 0.5323
                     u1 = 1.206
                     v1 = 0.0
                     p1 = 0.3
                  end if
               end if
               if (pox(1) .gt. zero) then
                  if (poy(1) .le. zero) then
                     r1 = 0.5323
                     u1 = 0.0
                     v1 = 1.206
                     p1 = 0.3
                  end if
                  if (poy(1) .gt. zero) then
                     r1 = 1.5
                     u1 = 0.0
                     v1 = 0.0
                     p1 = 1.5
                  end if
               end if
               skin1 = (oo2) * ((u1**2) + (v1**2))
               !internal energy
               ie1 = ((p1) / ((gamma - 1.00) * r1))
               !total energy
               e1 = (p1 / (gamma - 1)) + (r1 * skin1)
               !vector of conserved variables now
               ! rightv(1)=leftv(1)
               !rightv(2)=0.0
               !rightv(3)=0.0
               !rightv(4)=leftv(4)

               rightv(1) = r1
               rightv(2) = r1 * u1
               rightv(3) = r1 * v1
               rightv(4) = e1

            end if

         else

            call rotatef2d(n, cleft_rot, leftv, angle1, angle2)

            cright_rot(1) = cleft_rot(1)
            cright_rot(2) = -cleft_rot(2)
            cright_rot(3) = cleft_rot(3)
            cright_rot(4) = cleft_rot(4)

            if ((multispecies .eq. 1) .or. (realgas .eq. 1)) then
               cright_rot(5:nof_variables) = cleft_rot(5:nof_variables)

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               cturbr(:) = cturbl(:)

               if (passivescalar .gt. 0) then
                  cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                     & cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

               end if
            end if

            call rotateb2d(n, rightv, cright_rot, angle1, angle2)
         end if

      case (4) !wall

         if (itestcase .eq. 3) then

            call rotatef2d(n, cleft_rot, leftv, angle1, angle2)

            if ((multispecies .eq. 1) .or. (realgas .eq. 1)) then
               cright_rot(:) = cleft_rot(:)
               cright_rot(2) = -cleft_rot(2)

            else
               cright_rot(:) = cleft_rot(:)
               cright_rot(1) = cleft_rot(1)
               cright_rot(2) = -cleft_rot(2)
               cright_rot(3) = cleft_rot(3)
               cright_rot(4) = cleft_rot(4)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               cturbr(:) = cturbl(:)

               if (passivescalar .gt. 0) then
                  cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                     & cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

               end if
            end if

            call rotateb2d(n, rightv, cright_rot, angle1, angle2)

         else

            rightv(1:nof_variables) = leftv(1:nof_variables)

            rightv(2) = -leftv(2)
            rightv(3) = -leftv(3)

            if (realgas .eq. 1) then
               if (catalytic_wall .eq. 1) then
                  rightv(dimensiona + 4:nof_variables) = catalytic_con(1:nof_species) * leftv(1)

               end if
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if ((turbulence .ne. 1) .or. (turbulencemodel .ne. 2)) then
                  cturbr(:) = -cturbl(:)

                  if (passivescalar .gt. 0) then
                     cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                        & -cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

                  end if
               else
                  cturbr(1) = -cturbl(1)
                  cturbr(2) = 60.00 * visc / (beta_i1 * (ielem_walldist(iconsidered)**2))

                  if (passivescalar .gt. 0) then
                     cturbr(turbulenceequations + 1:turbulenceequations + passivescalar) = &
                        & -cturbl(turbulenceequations + 1:turbulenceequations + passivescalar)

                  end if

               end if
            end if

         end if

      case (6) !farfield inflow or outflow, subsonic or supersonic will be chosen based on mach number
         call rotatef2d(n, cleft_rot, leftv, angle1, angle2)
         vnb = cleft_rot(2) / cleft_rot(1)

         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         subson1(1:nof_variables) = rightv(1:nof_variables)
         subson2(1:nof_variables) = leftv(1:nof_variables)
         if (realgas .eq. 1) then
            sps = sqrt((gammal * subson2(4)) / (subson2(1) + tolsmall))
         else
            sps = sqrt((gamma * subson2(4)) / (subson2(1) + tolsmall))
         end if
         vel = sqrt(subson2(2)**2 + subson2(3)**2)

         call prim2cons2(n, leftv, rightv)

         if (vnb .le. 0.00) then !inflow
            ibfc = -1

            if ((abs(vnb)) .ge. sps) then
               !supersonic
               if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
                  call inflow_4440_from_left(n, leftv, [nx, ny, 0.0e0], rightv)
               else
                  call inflow2d(initcond, pox, poy, rightv)
               end if

            else
               !subsonic

               if ((initcond .eq. 4440) .and. (realgas .eq. 0) .and. (multispecies .eq. 0)) then
                  call inflow_4440_from_left(n, leftv, [nx, ny, 0.0e0], rightv)
               else
                  call inflow2d(initcond, pox, poy, rightv)
               end if

               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)
               if (realgas .eq. 1) then
                  sps = sqrt((gammal * subson2(4)) / (subson2(1) + tolsmall))
               else
                  sps = sqrt((gamma * subson2(4)) / (subson2(1) + tolsmall))
               end if
               vel = sqrt(subson2(2)**2 + subson2(3)**2)
               call prim2cons2(n, leftv, rightv)

               subson3(4) = 0.50 * ((subson1(4)) + (subson2(4)) - (subson2(1) * sps * ((nx * (subson1(2) - subson2(2))) + (ny * &
                  & (subson1(3) - subson2(3))))))
               subson3(1) = subson1(1) + (subson3(4) - subson1(4)) / (sps**2)
               subson3(2) = subson1(2) - (nx * (subson1(4) - subson3(4))) / (sps * subson2(1))
               subson3(3) = subson1(3) - (ny * (subson1(4) - subson3(4))) / (sps * subson2(1))

               if (realgas .eq. 1) then
                  rightv(1:nof_variables) = subson1(1:nof_variables)
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2)
                  rightv(3) = subson3(3)
                  rightv(4) = subson3(4)
                  call prim2cons(n, rightv)
               else
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2) * subson3(1)
                  rightv(3) = subson3(3) * subson3(1)

                  skins = oo2 * ((subson3(2)**2) + (subson3(3)**2))
                  ikins = subson3(4) / ((gamma - 1.00) * (subson3(1)))
                  rightv(4) = (subson3(1) * (ikins)) + (subson3(1) * skins)
               end if

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 1)) then
                  cturbr(1) = visc * turbinit * rightv(1)
               end if
               if ((turbulence .eq. 1) .and. (turbulencemodel .eq. 2)) then
                  cturbr(1) = (1.50 * i_turb_inlet * (ufreestream**2)) * rightv(1) !k initialization
                  cturbr(2) = rightv(1) * cturbr(1) / (10.0e-5 * visc) !omega initialization
               end if

               if (passivescalar .gt. 0) then
                  call pass_inlet2d(initcond, pox, poy, pass_tmp)
                  do rg_i = 1, passivescalar
                     cturbr(turbulenceequations + rg_i) = pass_tmp(rg_i) * rightv(1)
                  end do
               end if
            end if

         else

            !outflow

            ibfc = -2
            if ((abs(vnb)) .ge. sps) then

               rightv(1:nof_variables) = leftv(1:nof_variables)

            else

               call outflow2d(initcond, pox, poy, rightv)
               call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

               subson1(1:nof_variables) = rightv(1:nof_variables)
               subson2(1:nof_variables) = leftv(1:nof_variables)

               call prim2cons2(n, leftv, rightv)

               subson3(4) = subson1(4)
               subson3(1) = subson2(1) + (subson3(4) - subson2(4)) / (sps**2)
               subson3(2) = subson2(2) + (nx * (subson2(4) - subson3(4))) / (sps * subson2(1))
               subson3(3) = subson2(3) + (ny * (subson2(4) - subson3(4))) / (sps * subson2(1))

!
               if (realgas .eq. 1) then
                  rightv(1:nof_variables) = subson2(1:nof_variables)
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2)
                  rightv(3) = subson3(3)
                  rightv(4) = subson3(4)
                  call prim2cons(n, rightv)
               else
                  rightv(1) = subson3(1)
                  rightv(2) = subson3(2) * subson3(1)
                  rightv(3) = subson3(3) * subson3(1)

                  skins = oo2 * ((subson3(2)**2) + (subson3(3)**2))
                  ikins = subson3(4) / ((gamma - 1.00) * (subson3(1)))
                  rightv(4) = (subson3(1) * (ikins)) + (subson3(1) * skins)
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cturbr(:) = cturbl(:)

               end if
            end if

         end if

      end select

   end subroutine boundarys2d

   subroutine compute_eigenvectors(n, rveigl, rveigr, eigvl, eigvr, gamma)
      implicit none
!> @brief
!> this subroutine computes the left and right eigenvectors
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(in) :: rveigl, rveigr
      real, intent(in) :: gamma
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar), intent(inout) :: eigvl, eigvr
      real :: rs, us, vs, ws, es, ps, vvs, as, hs, gammam1, vsd, oor1, oor2
      integer :: ivgt

      eigvr = zero
      gammam1 = gamma - 1.00
      oor1 = 1.00 / rveigl(1)
      oor2 = 1.00 / rveigr(1)

      rs = oo2 * (rveigl(1) + rveigr(1))
      us = oo2 * ((rveigl(2) * oor1) + (rveigr(2) * oor2))
      vs = oo2 * ((rveigl(3) * oor1) + (rveigr(3) * oor2))
      ws = oo2 * ((rveigl(4) * oor1) + (rveigr(4) * oor2))
      es = oo2 * (rveigl(5) + rveigr(5))

      vvs = (us**2) + (vs**2) + (ws**2)
      vsd = oo2 * vvs
      ps = (gamma - 1.00) * (es - oo2 * rs * vvs)
      as = sqrt(gamma * ps / rs)
      hs = (oo2 * vvs) + ((as**2) / gammam1)
      eigvr(1, 1) = 1.00; eigvr(1, 2) = 1.00; eigvr(1, 3) = 0.00; eigvr(1, 4) = 0.00; eigvr(1, 5) = 1.00
      eigvr(2, 1) = us - as; eigvr(2, 2) = us; eigvr(2, 3) = 0.00; eigvr(2, 4) = 0.00; eigvr(2, 5) = us + as
      eigvr(3, 1) = vs; eigvr(3, 2) = vs; eigvr(3, 3) = 1.00; eigvr(3, 4) = 0.00; eigvr(3, 5) = vs
      eigvr(4, 1) = ws; eigvr(4, 2) = ws; eigvr(4, 3) = 0.00; eigvr(4, 4) = 1.00; eigvr(4, 5) = ws
      eigvr(5, 1) = hs - (us * as); eigvr(5, 2) = oo2 * vvs; eigvr(5, 3) = vs; eigvr(5, 4) = ws; eigvr(5, 5) = hs + (us * as)

      eigvl(1, 1) = (hs * us + as * us**2 + as * vs**2 - as * vsd - us * vsd + as * ws**2) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(1, 2) = (-hs - as * us + vsd) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(1, 3) = -((as * vs) / (2.00 * as * hs - 2.00 * as * vsd))
      eigvl(1, 4) = -((as * ws) / (2.00 * as * hs - 2.00 * as * vsd))
      eigvl(1, 5) = as / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(2, 1) = (2.0 * as * hs - 2.0 * as * us**2 - 2.00 * as * vs**2 - 2.00 * as * ws**2) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(2, 2) = (2.00 * as * us) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(2, 3) = (2.00 * as * vs) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(2, 4) = (2.00 * as * ws) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(2, 5) = (-2.00 * as) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(3, 1) = (-2.0 * as * hs * vs + 2.0 * as * vs * vsd) / (2.0 * as * hs - 2.0 * as * vsd)
      eigvl(3, 2) = 0.00
      eigvl(3, 3) = 1.00
      eigvl(3, 4) = 0.00
      eigvl(3, 5) = 0.00
      eigvl(4, 1) = (-2.00 * as * hs * ws + 2.00 * as * vsd * ws) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(4, 2) = 0.00
      eigvl(4, 3) = 0.00
      eigvl(4, 4) = 1.00
      eigvl(4, 5) = 0.00
      eigvl(5, 1) = (-(hs * us) + as * us**2 + as * vs**2 - as * vsd + us * vsd + as * ws**2) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(5, 2) = (hs - as * us - vsd) / (2.00 * as * hs - 2.00 * as * vsd)
      eigvl(5, 3) = -((as * vs) / (2.00 * as * hs - 2.00 * as * vsd))
      eigvl(5, 4) = -((as * ws) / (2.00 * as * hs - 2.00 * as * vsd))
      eigvl(5, 5) = as / (2.00 * as * hs - 2.00 * as * vsd)

   end subroutine compute_eigenvectors

   subroutine compute_jacobianse(n, iconsidered, eigvl, rveigl, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
      implicit none
!> @brief
!> this subroutine computes the jacobians for the implicit time stepping
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, intent(in) :: angle1, angle2, nx, ny, nz
      real, dimension(1:nof_variables), intent(in) :: rveigl, srf_speedrot
      real, intent(in) :: gamma
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar), intent(inout) :: eigvl
      real :: rs, us, vs, ws, es, ps, vvs, as, hs, gammam1, vsd, phi, a1, a2, a3, oors, vn
      integer :: ivgt, i, idxev
      real :: evs, etr, yk, echem, dpdrho
      real :: temps, gammal, mp_pinfl
      integer :: rg_i
      real, dimension(1:gpu_max_nvar) :: leftv

      if (realgas .eq. 0) then

         a2 = gamma - 1.00
         a3 = gamma - 2.00
         oors = 1.00 / rveigl(1)
         rs = (rveigl(1))
         us = (rveigl(2) * oors)
         vs = (rveigl(3) * oors)
         ws = (rveigl(4) * oors)
         es = (rveigl(5) * oors)
         phi = oo2 * (a2) * ((us * us) + (vs * vs) + (ws * ws))
         a1 = gamma * es - phi

         vvs = nx * us + ny * vs + nz * ws

         if (rec_mrf(iconsidered) .eq. 1) then
            eigvl(1, 1) = 0.00 - srf_speedrot(2); eigvl(1, 2) = nx; eigvl(1, 3) = ny; eigvl(1, 4) = nz; eigvl(1, 5) = 0.00
            eigvl(2, 1) = nx * phi - us * vvs; eigvl(2, 2) = vvs - a3 * nx * us - srf_speedrot(2); eigvl(2, &
               & 3) = ny * us - a2 * nx * vs; eigvl(2, 4) = nz * us - a2 * nx * ws; eigvl(2, 5) = a2 * nx
            eigvl(3, 1) = ny * phi - vs * vvs; eigvl(3, 2) = nx * vs - a2 * ny * us; eigvl(3, 3) = vvs - a3 * ny * vs - &
               & srf_speedrot(2); eigvl(3, 4) = nz * vs - a2 * ny * ws; eigvl(3, 5) = a2 * ny
            eigvl(4, 1) = nz * phi - ws * vvs; eigvl(4, 2) = nx * ws - a2 * nz * us; eigvl(4, 3) = ny * ws - a2 * nz * vs
            eigvl(4, 4) = vvs - a3 * nz * ws - srf_speedrot(2); eigvl(4, 5) = a2 * nz
            eigvl(5, 1) = vvs * (phi - a1); eigvl(5, 2) = nx * a1 - a2 * us * vvs; eigvl(5, 3) = ny * a1 - a3 * vs * vvs; eigvl(5, &
               & 4) = nz * a1 - a2 * ws * vvs; eigvl(5, 5) = gamma * vvs - srf_speedrot(2)
         else
            eigvl(1, 1) = 0.00; eigvl(1, 2) = nx; eigvl(1, 3) = ny; eigvl(1, 4) = nz; eigvl(1, 5) = 0.00
            eigvl(2, 1) = nx * phi - us * vvs; eigvl(2, 2) = vvs - a3 * nx * us; eigvl(2, 3) = ny * us - a2 * nx * vs; eigvl(2, &
               & 4) = nz * us - a2 * nx * ws; eigvl(2, 5) = a2 * nx
            eigvl(3, 1) = ny * phi - vs * vvs; eigvl(3, 2) = nx * vs - a2 * ny * us; eigvl(3, 3) = vvs - a3 * ny * vs; eigvl(3, &
               & 4) = nz * vs - a2 * ny * ws; eigvl(3, 5) = a2 * ny
            eigvl(4, 1) = nz * phi - ws * vvs; eigvl(4, 2) = nx * ws - a2 * nz * us; eigvl(4, 3) = ny * ws - a2 * nz * vs
            eigvl(4, 4) = vvs - a3 * nz * ws; eigvl(4, 5) = a2 * nz
            eigvl(5, 1) = vvs * (phi - a1); eigvl(5, 2) = nx * a1 - a2 * us * vvs; eigvl(5, 3) = ny * a1 - a2 * vs * vvs; eigvl(5, &
               & 4) = nz * a1 - a2 * ws * vvs; eigvl(5, 5) = gamma * vvs
         end if

      end if

      if (realgas .eq. 1) then

         call realgas_euler_jacobian(n, rveigl, [nx, ny, nz], eigvl)

      end if

   end subroutine compute_jacobianse

   subroutine compute_eigenvectors2d(n, rveigl, rveigr, eigvl, eigvr, gamma)
      implicit none
!> @brief
!> this subroutine computes the left and right eigenvectors  in 2d
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:nof_variables), intent(in) :: rveigl, rveigr
      real, intent(in) :: gamma
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar), intent(inout) :: eigvl, eigvr
      real :: rs, us, vs, ws, es, ps, vvs, as, hs, gammam1, vsd, g8, s1, s2, vtots, oor1, oor2
      integer :: ivgt, j, k

      eigvr = zero
      gammam1 = gamma - 1.00
      oor1 = 1.00 / rveigl(1)
      oor2 = 1.00 / rveigr(1)

      eigvr = 0.00
      gammam1 = gamma - 1.00
      rs = 0.5 * (rveigl(1) + rveigr(1))
      us = 0.5 * ((rveigl(2) * oor1) + (rveigr(2) * oor2))
      vs = 0.5 * ((rveigl(3) * oor1) + (rveigr(3) * oor2))
      es = 0.5 * (rveigl(4) + rveigr(4))
      g8 = gamma - 1.00

      vtots = us**2 + vs**2
      ps = (gamma - 1) * (es - oo2 * rs * vtots)
      as = sqrt(gamma * ps / rs)
      hs = oo2 * vtots + (as**2) / (g8)

      eigvr(1, 1) = 1.0; eigvr(1, 2) = 1.0; eigvr(1, 3) = 0.0; eigvr(1, 4) = 1.0
      eigvr(2, 1) = us - as; eigvr(2, 2) = us; eigvr(2, 3) = 0.0; eigvr(2, 4) = us + as
      eigvr(3, 1) = vs; eigvr(3, 2) = vs; eigvr(3, 3) = 1.0; eigvr(3, 4) = vs
      eigvr(4, 1) = hs - us * as; eigvr(4, 2) = oo2 * vtots; eigvr(4, 3) = vs; eigvr(4, 4) = hs + us * as

      s1 = as / (gamma - 1.0)
      s2 = as**2 / (gamma - 1.0)

      eigvl(1, 1) = hs + s1 * (us - as); eigvl(1, 2) = -(us + s1); eigvl(1, 3) = -vs; eigvl(1, 4) = 1.0
      eigvl(2, 1) = -2.0 * hs + 4.0 * s2; eigvl(2, 2) = 2.0 * us; eigvl(2, 3) = 2.0 * vs; eigvl(2, 4) = -2.0
      eigvl(3, 1) = -2.0 * vs * s2; eigvl(3, 2) = 0.0; eigvl(3, 3) = 2.0 * s2; eigvl(3, 4) = 0.0
      eigvl(4, 1) = hs - s1 * (us + as); eigvl(4, 2) = -us + s1; eigvl(4, 3) = -vs; eigvl(4, 4) = 1.0

      do j = 1, 4
         do k = 1, 4

            eigvl(j, k) = eigvl(j, k) / (2.0 * s2)
         end do
      end do

   end subroutine compute_eigenvectors2d

   subroutine compute_jacobianse2d(n, eigvl, rveigl, gamma, angle1, angle2, nx, ny, nz)
      implicit none
      !> @brief
!> this subroutine computes the jacobians for the implicit time stepping in 2d
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: angle1, angle2, nx, ny, nz
      real, dimension(1:nof_variables), intent(in) :: rveigl
      real, intent(in) :: gamma
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar), intent(inout) :: eigvl
      real :: rs, us, vs, es, ps, vvs, as, hs, gammam1, vsd, phi, a1, a2, a3, oors
      integer :: ivgt
      real :: evs, etr, yk, echem, dpdrho
      real :: temps, gammal, mp_pinfl
      integer :: i, rg_i, row, k, idxev
      real, dimension(1:gpu_max_nvar) :: leftv

      if (realgas .eq. 0) then
         a2 = gamma - 1.00
         a3 = gamma - 2.00
         oors = 1.00 / rveigl(1)
         rs = (rveigl(1))
         us = (rveigl(2) * oors)
         vs = (rveigl(3) * oors)
         es = (rveigl(4) * oors)
         phi = oo2 * (a2) * ((us * us) + (vs * vs))
         a1 = gamma * es - phi

         vvs = nx * us + ny * vs

         eigvl(1, 1) = 0.00; eigvl(1, 2) = nx; eigvl(1, 3) = ny; eigvl(1, 4) = 0.00
         eigvl(2, 1) = nx * phi - us * vvs; eigvl(2, 2) = vvs - a3 * nx * us; eigvl(2, 3) = ny * us - a2 * nx * vs; eigvl(2, &
            & 4) = a2 * nx
         eigvl(3, 1) = ny * phi - vs * vvs; eigvl(3, 2) = nx * vs - a2 * ny * us; eigvl(3, 3) = vvs - a3 * ny * vs; eigvl(3, &
            & 4) = a2 * ny

         eigvl(4, 1) = vvs * (phi - a1); eigvl(4, 2) = nx * a1 - a2 * us * vvs; eigvl(4, 3) = ny * a1 - a2 * vs * vvs; eigvl(4, &
            & 4) = gamma * vvs

      end if

      if (realgas .eq. 1) then
         call realgas_euler_jacobian(n, rveigl, [nx, ny, nz], eigvl)

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

   subroutine eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: leftv(nof_variables), rightv(nof_variables)
      real, intent(inout) :: viscl(4), laml(4), turbmv(2), etvm(1), eddyfl(20), eddyfr(20)
      call eddy_transport_pair(3, leftv, rightv, viscl, laml, turbmv, etvm, eddyfl, eddyfr)
   end subroutine eddyvisco

   subroutine eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: leftv(nof_variables), rightv(nof_variables)
      real, intent(inout) :: viscl(4), laml(4), turbmv(2), etvm(1), eddyfl(20), eddyfr(20)
      call eddy_transport_pair(2, leftv, rightv, viscl, laml, turbmv, etvm, eddyfl, eddyfr)
   end subroutine eddyvisco2d

   subroutine trajectories
      implicit none
      integer :: i, j, k, traj1, traj2, traj3, traj4, kmaxe, writeid, writeconf, num_vg
      real :: win1, win2, win3, win4, post, post1, post2, post3, post4
      real, dimension(1:nof_variables) :: leftv
      real :: mp_pinfl, gammal
      kmaxe = xmpielrank(n)
      post1 = tolbig
      post2 = tolbig
      post3 = -tolbig
      traj1 = 0
      traj2 = 0
      traj3 = 0

      if (dimensiona .eq. 3) then
         num_vg = 5
      else
         num_vg = 4

      end if

      if (initcond .eq. 157) then
         pos_l(1:2) = zero
         pos_g(1:2) = zero
         do i = 1, kmaxe
            leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            if (dimensiona .eq. 3) then
               call cons2prim(n, leftv, mp_pinfl, gammal)
            else
               call cons2prim(n, leftv, mp_pinfl, gammal)
            end if

            if ((leftv(nof_variables)) .gt. 0.00) then !total volume of gas evolution
               pos_l(2) = pos_l(2) + (leftv(nof_variables)) * ielem_totvolume(i)
            end if
            if (leftv(num_vg) .gt. pos_l(1)) then !total volume of gas evolution
               pos_l(1) = max(pos_l(1), leftv(num_vg))
            end if

         end do
         !find position globally
         call mpi_allreduce(pos_l(1), pos_g(1), 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)
         call mpi_allreduce(pos_l(2), pos_g(2), 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)

         if (n .eq. 0) then

            open (70, file='volumex.dat', form='formatted', action='write', position='append')
            write (70, '(e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), pos_g(2)
            close (70)

         end if

         call mpi_barrier(mpi_comm_world, ierror)

      end if

      if (initcond .eq. 430) then
         pos_l(1:2) = zero
         pos_g(1:2) = zero
         do i = 1, kmaxe
            if (u_c_val(i, 1, 7) .gt. 0.10) then !volume fraction of gas to be used for lowest location tracking
               if (ielem_yyc(i) .le. post3) then
                  post3 = ielem_yyc(i)
                  traj1 = i
               end if
            end if
            if (u_c_val(i, 1, 7) .gt. 0.00) then !total volume of gas evolution
               pos_l(2) = pos_l(2) + u_c_val(i, 1, 7) * ielem_totvolume(i)
            end if
         end do

         pos_l(1) = post3 !lowest position in my local  cpu

!find position globally
         call mpi_barrier(mpi_comm_world, ierror)
         call mpi_allreduce(pos_l(1), pos_g(1), 1, mpi_double_precision, mpi_min, mpi_comm_world, ierror)

!now select which cpu has the smallest
         writeid = 100000000
         if (abs(pos_g(1) - pos_l(1)) .le. tolsmall / 1000) then
            if (traj1 .gt. 0) then
               writeid = n
            end if
         end if
!and if more than one, select the smallest id one
         call mpi_barrier(mpi_comm_world, ierror)
         call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)

         if (writeid .eq. writeconf) then

            leftv(1:nof_variables) = u_c_val(traj1, 1, 1:nof_variables)

            call cons2prim(n, leftv, mp_pinfl, gammal)
            open (70, file='position.dat', form='formatted', action='write', position='append')
            write (70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), leftv(1), &
               & leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7)
            close (70)
!the cpu that holds this cell will write its position
         end if

         call mpi_barrier(mpi_comm_world, ierror)

         call mpi_allreduce(pos_l(2:2), pos_g(2:2), 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)

         if (n .eq. 0) then

            open (70, file='volume.dat', form='formatted', action='write', position='append')
            write (70, '(e14.7,1x,e14.7,1x,e14.7)') t, pos_g(2)
            close (70)

         end if

         call mpi_barrier(mpi_comm_world, ierror)

      end if

      if (initcond .eq. 405) then

         do i = 1, kmaxe
            if (dimensiona .eq. 3) then
               if (u_c_val(i, 1, 8) .gt. 0.10) then
                  if (ielem_xxc(i) .le. post1) then
                     post1 = ielem_xxc(i)
                     traj1 = i
                  end if
                  if (((ielem_yyc(i) .le. 0.0515) .and. (ielem_yyc(i) .ge. 0.0485)) .and. ((ielem_zzc(i) .le. 0.0515) .and. &
                     & (ielem_zzc(i) .ge. 0.0485))) then
                     if (ielem_xxc(i) .le. post2) then
                        post2 = ielem_xxc(i)
                        traj2 = i
                     end if
                     if (ielem_xxc(i) .ge. post3) then
                        post3 = ielem_xxc(i)
                        traj3 = i
                     end if
                  end if
               end if
            else
               if (u_c_val(i, 1, 7) .ge. 0.40) then
                  if (ielem_xxc(i) .le. post1) then
                     post1 = ielem_xxc(i)
                     traj1 = i
                  end if
                  if (((ielem_yyc(i) .le. 0.053) .and. (ielem_yyc(i) .ge. 0.048))) then
                     if (ielem_xxc(i) .le. post2) then
                        post2 = ielem_xxc(i)
                        traj2 = i
                     end if
                     if (ielem_xxc(i) .ge. post3) then
                        post3 = ielem_xxc(i)
                        traj3 = i
                     end if
                  end if
               end if

            end if
         end do

         pos_l(1) = post1
         pos_l(2) = post2
         pos_l(3) = post3

         call mpi_barrier(mpi_comm_world, ierror)

         call mpi_allreduce(pos_l(1:2), pos_g(1:2), 2, mpi_double_precision, mpi_min, mpi_comm_world, ierror)

!now select which cpu has the smallest
         writeid = 100000000
         writeconf = -1
         if (abs(pos_g(1) - pos_l(1)) .le. tolsmall / 1000) then
            if (traj1 .gt. 0) then
               writeid = n
            end if
         end if
!and if more than one, select the smallest id one

!if (traj1.gt.0)then
         call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
!end if

         if (writeid .eq. writeconf) then

            if (traj1 .gt. 0) then
               leftv(1:nof_variables) = u_c_val(traj1, 1, 1:nof_variables)

               if (dimensiona .eq. 3) then
                  call cons2prim(n, leftv, mp_pinfl, gammal)

                  open (70, file='pos1.dat', form='formatted', action='write', position='append')
                  write (70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, &
                     & pos_g(1), leftv(1), leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7), leftv(8)
                  close (70)
               else
                  call cons2prim(n, leftv, mp_pinfl, gammal)

                  open (70, file='pos1.dat', form='formatted', action='write', position='append')
                  write (70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), &
                     & leftv(1), leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7)
                  close (70)

               end if

            end if

         end if

         writeid = 100000000
         writeconf = -1
         if (abs(pos_g(2) - pos_l(2)) .le. tolsmall / 1000) then
            if (traj2 .gt. 0) then
               writeid = n
            end if
         end if

!if (traj2.gt.0)then
         call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
!end if

         if (writeid .eq. writeconf) then

            if (traj2 .gt. 0) then
               leftv(1:nof_variables) = u_c_val(traj2, 1, 1:nof_variables)

               if (dimensiona .eq. 3) then

                  call cons2prim(n, leftv, mp_pinfl, gammal)

                  open (71, file='pos2.dat', form='formatted', action='write', position='append')
                  write (71, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, &
                     & pos_g(2), leftv(1), leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7), leftv(8)
                  close (71)

               else
                  call cons2prim(n, leftv, mp_pinfl, gammal)

                  open (70, file='pos2.dat', form='formatted', action='write', position='append')
                  write (70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, pos_g(2), &
                     & leftv(1), leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7)
                  close (70)

               end if

            end if
         end if

         call mpi_allreduce(pos_l(3), pos_g(3), 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)

         writeid = 100000000
         writeconf = -1
         if (abs(pos_g(3) - pos_l(3)) .le. tolsmall / 1000) then
            if (traj3 .gt. 0) then
               writeid = n
            end if
         end if

!if (traj3.gt.0)then
         call mpi_allreduce(writeid, writeconf, 1, mpi_integer, mpi_min, mpi_comm_world, ierror)
!end if

         if (writeid .eq. writeconf) then
            if (traj3 .gt. 0) then
               leftv(1:nof_variables) = u_c_val(traj3, 1, 1:nof_variables)

               if (dimensiona .eq. 3) then

                  call cons2prim(n, leftv, mp_pinfl, gammal)

                  open (71, file='pos3.dat', form='formatted', action='write', position='append')
                  write (71, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, &
                     & pos_g(3), leftv(1), leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7), leftv(8)
                  close (71)

               else
                  call cons2prim(n, leftv, mp_pinfl, gammal)

                  open (70, file='pos3.dat', form='formatted', action='write', position='append')
                  write (70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, pos_g(3), &
                     & leftv(1), leftv(2), leftv(3), leftv(4), leftv(5), leftv(6), leftv(7)
                  close (70)

               end if

            end if

         end if

         call mpi_barrier(mpi_comm_world, ierror)

      end if

      if ((initcond .eq. 408) .or. (initcond .eq. 422)) then

         pos_l(1:2) = zero
         pos_g(1:2) = zero

         do i = 1, kmaxe
            if (u_c_val(i, 1, 8) .gt. 0.00) then
               pos_l(1) = pos_l(1) + ielem_totvolume(i)
               pos_l(2) = pos_l(2) + u_c_val(i, 1, 8) * ielem_totvolume(i)
            end if
         end do

         call mpi_barrier(mpi_comm_world, ierror)

         call mpi_allreduce(pos_l(1:2), pos_g(1:2), 2, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)

         if (n .eq. 0) then

            open (70, file='volume.dat', form='formatted', action='write', position='append')
            write (70, '(e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), pos_g(2)
            close (70)

         end if

         call mpi_barrier(mpi_comm_world, ierror)

      end if

      if ((initcond .eq. 411) .or. (initcond .eq. 444)) then

         pos_l(1:3) = zero
         pos_g(1:3) = zero

         do i = 1, kmaxe
            leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            if (dimensiona .eq. 2) then
               call cons2prim(n, leftv, mp_pinfl, gammal)
               pos_l(3) = max(abs(leftv(4)), pos_l(3))

               if (u_c_val(i, 1, 7) .gt. 0.00) then
                  pos_l(1) = pos_l(1) + ielem_totvolume(i)
                  pos_l(2) = pos_l(2) + u_c_val(i, 1, 7) * ielem_totvolume(i)
               end if
            else
               call cons2prim(n, leftv, mp_pinfl, gammal)
               pos_l(3) = max(abs(leftv(5)), pos_l(3))

               if (u_c_val(i, 1, 8) .gt. 0.00) then
                  pos_l(1) = pos_l(1) + ielem_totvolume(i)
                  pos_l(2) = pos_l(2) + u_c_val(i, 1, 7) * ielem_totvolume(i)
               end if

            end if

         end do

         call mpi_barrier(mpi_comm_world, ierror)

         call mpi_allreduce(pos_l(1:2), pos_g(1:2), 2, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)

         call mpi_allreduce(pos_l(3), pos_g(3), 1, mpi_double_precision, mpi_max, mpi_comm_world, ierror)

         if (n .eq. 0) then

            open (70, file='volume.dat', form='formatted', action='write', position='append')
            write (70, '(e14.7,1x,e14.7,1x,e14.7,1x,e14.7)') t, pos_g(1), pos_g(2), pos_g(3)
            close (70)

         end if

         call mpi_barrier(mpi_comm_world, ierror)

      end if

   end subroutine

   subroutine flux2dx(flux_term_x, leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: flux_term_x(nof_variables)
      call primitive_euler_flux(1, leftv, flux_term_x)
   end subroutine flux2dx

   subroutine flux2dy(flux_term_y, leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: flux_term_y(nof_variables)
      call primitive_euler_flux(2, leftv, flux_term_y)
   end subroutine flux2dy

   subroutine flux3dx(flux_term_x, leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: flux_term_x(nof_variables)
      call primitive_euler_flux(1, leftv, flux_term_x)
   end subroutine flux3dx

   subroutine flux3dy(flux_term_y, leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: flux_term_y(nof_variables)
      call primitive_euler_flux(2, leftv, flux_term_y)
   end subroutine flux3dy

   subroutine flux3dz(flux_term_z, leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: leftv(nof_variables)
      real, intent(out) :: flux_term_z(nof_variables)
      call primitive_euler_flux(3, leftv, flux_term_z)
   end subroutine flux3dz

   subroutine flux_visc2d(flux_term_x, flux_term_y, leftv, leftv_der)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real :: u, v, ux, uy, vx, vy, tx, ty, tauxx, tauxy, tauyy
      real, dimension(1:nof_variables), intent(inout) :: leftv
      real, dimension(1:nof_variables), intent(inout) :: flux_term_x, flux_term_y
      real, dimension(1:4) :: viscl, laml
      ! Match the active spatial extent; the leading workspace stride stays fixed.
      real, dimension(1:gpu_max_nvar, 1:dimensiona), intent(in) :: leftv_der

      call get_visc_conduct(n, leftv, leftv, viscl, laml)

      u = leftv(2) / leftv(1)
      v = leftv(3) / leftv(1)
      ux = leftv_der(2, 1)
      uy = leftv_der(2, 2)
      vx = leftv_der(3, 1)
      vy = leftv_der(3, 2)
      tx = leftv_der(4, 1)
      ty = leftv_der(4, 2)

      tauxx = 2.00 / 3.00 * viscl(1) * (2 * ux - vy)
      tauyy = 2.00 / 3.00 * viscl(1) * (2 * vy - ux)
      tauxy = viscl(1) * (uy + vx)

      flux_term_x(2) = flux_term_x(2) - tauxx
      flux_term_x(3) = flux_term_x(3) - tauxy
      flux_term_x(4) = flux_term_x(4) - (u * tauxx + v * tauxy + laml(1) * tx)
      flux_term_y(2) = flux_term_y(2) - tauxy
      flux_term_y(3) = flux_term_y(3) - tauyy
      flux_term_y(4) = flux_term_y(4) - (u * tauxy + v * tauyy + laml(1) * ty)

   end subroutine

   subroutine flux_visc3d(flux_term_x, flux_term_y, flux_term_z, leftv, leftv_der)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real :: u, v, ux, uy, vx, vy, tx, ty, tauxx, tauxy, tauyy, w, uz, vz, wx, wy, wz, tz, tauxz, tauyz, tauzz
      real, dimension(1:nof_variables), intent(inout) :: leftv
      real, dimension(1:nof_variables), intent(inout) :: flux_term_x, flux_term_y, flux_term_z
      ! Match the active spatial extent; the leading workspace stride stays fixed.
      real, dimension(1:gpu_max_nvar, 1:dimensiona), intent(in) :: leftv_der
      real, dimension(1:4) :: viscl, laml

      call get_visc_conduct(n, leftv, leftv, viscl, laml)

! variables extrapolated at boundary
      u = leftv(2) / leftv(1)
      v = leftv(3) / leftv(1)
      w = leftv(4) / leftv(1)

      ux = leftv_der(2, 1)
      uy = leftv_der(2, 2)
      uz = leftv_der(2, 3)

      vx = leftv_der(3, 1)
      vy = leftv_der(3, 2)
      vz = leftv_der(3, 3)

      wx = leftv_der(4, 1)
      wy = leftv_der(4, 2)
      wz = leftv_der(4, 3)

      tx = leftv_der(5, 1)
      ty = leftv_der(5, 2)
      tz = leftv_der(5, 3)

      tauxx = 2.00 / 3.00 * viscl(1) * (2 * ux - vy - wz)
      tauyy = 2.00 / 3.00 * viscl(1) * (2 * vy - ux - wz)
      tauxy = viscl(1) * (uy + vx)

      tauxz = viscl(1) * (uz + wx)
      tauyz = viscl(1) * (wy + vz)
      tauzz = 2.00 / 3.00 * viscl(1) * (2 * wz - ux - vy)

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

   subroutine dcons2dprim(leftv_der, leftv)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      ! Derivative of cons2div: density, velocities, temperature(s), species.
      ! The input state and gradients are conservative. Work one direction at
      ! a time; no thermodynamic Jacobian or coefficient-sized scratch is needed.
      real, intent(inout) :: leftv_der(gpu_max_nvar, gpu_max_dim)
      real, intent(in) :: leftv(nof_variables)
      integer :: d, k, ie, is, iv
      real :: rho, drho, velocity, du, energy, de, ev, dev, etr, rmix, cv, cp, temp, tv
      real :: y(gpu_max_species), dy, dcv, dchem, dvib, cvv, theta, arg, ex, es, cvs, rs
      real :: total, dtotal, h, correction, velocity_rho, velocity_drho
      ie = dimensiona + 2
      rho = leftv(1)
      if (nof_variables <= 1) return
      if (realgas == 1) then
         call realgas_thermo(leftv, rho, y, ev, etr, rmix, cv)
         temp = realgas_ttr(etr, cv)
         tv = realgas_vibrational_temperature(y, ev, temp)
      elseif (multispecies == 1) then
         call multispecies_thermo(leftv, temp, y, rmix, cp, cv, energy)
         total = sum(leftv(ie + 1:ie + nof_species))
      else
         cv = r_gas / (gamma - 1.0e0)
      end if
      do d = 1, dimensiona
         drho = leftv_der(1, d)
         de = (leftv_der(ie, d) - leftv(ie) * drho / rho) / rho
         do k = 2, dimensiona + 1
            velocity = leftv(k) / rho
            du = (leftv_der(k, d) - velocity * drho) / rho
            de = de - velocity * du
         end do
         if (realgas == 1) then
            dev = (leftv_der(ie + 1, d) - ev * drho) / rho
            dcv = zero; dchem = zero; dvib = zero; cvv = zero
            do k = 1, nof_species
               is = ie + 1 + k
               dy = (leftv_der(is, d) - y(k) * drho) / rho
               rs = rgs_ru / rg_molm(k)
               dcv = dcv + merge(2.5e0, 1.5e0, k <= 3) * rs * dy
               dchem = dchem + max(rg_hzero(k), zero) / rg_molm(k) * dy
               if (k <= 3) then
                  theta = rg_thetag(k)
                  if (theta > zero) then
                     arg = theta / tv
                     if (arg > 60.0e0) then
                        es = rs * theta * exp(-arg)
                        cvs = es * arg / tv
                     else
                        ex = exp(arg)
                        es = rs * theta / (ex - 1.0e0)
                        cvs = es * arg * ex / ((ex - 1.0e0) * tv)
                     end if
                     dvib = dvib + es * dy
                     cvv = cvv + y(k) * cvs
                  end if
               end if
               leftv_der(is, d) = dy
            end do
            leftv_der(ie, d) = (de - dev - dchem - temp * dcv) / cv
            leftv_der(ie + 1, d) = zero
            if (cvv > 1.0e-20) leftv_der(ie + 1, d) = (dev - dvib) / cvv
         elseif (multispecies == 1) then
            dtotal = sum(leftv_der(ie + 1:ie + nof_species, d))
            correction = zero
            do k = 1, nof_species
               dy = (leftv_der(ie + k, d) - y(k) * dtotal) / total
               call species_cp_h_sensible(temp, k, cp, h)
               es = h - (rgs_ru / mp_m(k)) * temp
               correction = correction + es * dy
            end do
            leftv_der(ie, d) = (de - correction) / cv
         else
            leftv_der(ie, d) = de / cv
         end if
         velocity_rho = rho; velocity_drho = drho
         if (multispecies == 1) then
            if (mp_modelc == 0) then
               velocity_rho = total; velocity_drho = dtotal
               leftv_der(1, d) = dtotal
            end if
         end if
         do iv = 2, dimensiona + 1
            velocity = leftv(iv) / velocity_rho
            leftv_der(iv, d) = (leftv_der(iv, d) - velocity * velocity_drho) / velocity_rho
         end do
      end do
   end subroutine dcons2dprim

end module flow_operations
