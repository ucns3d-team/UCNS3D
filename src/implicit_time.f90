module implicit_time
   use library
   use transform
   use flow_operations
   use implicit_fluxes
   use communications
   use declaration
   implicit none

#ifdef xpu
#define RELAX_NVAR gpu_max_nvar
#define RELAX_EXTRA gpu_max_extra_transport
#define RELAX_NTOTAL gpu_max_nvar_total
#define RELAX_DIM gpu_max_dim
#else
#define RELAX_NVAR nof_variables
#define RELAX_EXTRA turbulenceequations+passivescalar
#define RELAX_NTOTAL nof_variables+turbulenceequations+passivescalar
#define RELAX_DIM dimensiona
#endif

#ifndef XPU_IMPLICIT_BATCH
#define XPU_IMPLICIT_BATCH 16384
#endif
   integer, parameter :: implicit_lm_batch = XPU_IMPLICIT_BATCH
   real, allocatable :: lm_diag(:, :, :), lm_off(:, :, :, :), lm_diagt(:, :), lm_offt(:, :, :)
#ifdef xpu
!$omp declare target (lm_diag, lm_off, lm_diagt, lm_offt)
#endif
   integer, allocatable :: sweep_order(:), sweep_offset(:)
   integer :: sweep_levels = 0, sweep_cells = -1
#ifdef xpu
!$omp declare target (sweep_order)
#endif

contains

   subroutine invert_implicit_block(i, block)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      real, intent(inout) :: block(:, :)
      integer :: j
      ! Chemistry, rotation and the duct outlet couple equations; diagonal blocks need
      ! no private matrix copy, including additional species equations.
      if (realgas == 1 .or. (dimensiona == 3 .and. rframe /= 0) .or. duct_outlet_face_flux_enabled()) then
         call invert_realgas_block(i, block)
      else
         do j = 1, nof_variables
            block(j, j) = 1.0e0 / block(j, j)
         end do
      end if
   end subroutine invert_implicit_block


   subroutine prepare_implicit_workspace(ncells)
      implicit none
      integer, intent(in) :: ncells
      integer :: slots, extra
      slots = max(1, min(implicit_lm_batch, ncells))
      extra = max(1, turbulenceequations + passivescalar)
      if (allocated(lm_diag)) then
         if (size(lm_diag, 1) == slots .and. size(lm_diag, 2) == nof_variables .and. &
            & size(lm_diagt, 2) == extra) return
#ifdef xpu
!$omp target exit data map(delete: lm_diag, lm_off, lm_diagt, lm_offt)
#endif
         deallocate(lm_diag, lm_off, lm_diagt, lm_offt)
      end if
      allocate(lm_diag(slots, nof_variables, nof_variables), &
         & lm_off(slots, gpu_max_faces, nof_variables, nof_variables), &
         & lm_diagt(slots, extra), lm_offt(slots, gpu_max_faces, extra))
#ifdef xpu
!$omp target enter data map(alloc: lm_diag, lm_off, lm_diagt, lm_offt)
#endif
   end subroutine prepare_implicit_workspace

   subroutine implicit_sweep_range(n, first, last, direction, low_storage, ordered)
      implicit none
      integer, intent(in) :: n, first, last, direction
      logical, intent(in) :: low_storage, ordered
      integer :: chunk, chunk_end, stride, at, i, slot
      stride = max(1, last - first + 1)
      if (low_storage) stride = size(lm_diag, 1)
      do chunk = first, last, stride
         chunk_end = min(last, chunk + stride - 1)
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n, chunk, chunk_end, direction, low_storage, ordered) private(i, slot) &
!$omp& map(alloc: sweep_order) &
!$omp& private(at) &
!$omp& map(alloc: bleed_end, bleed_plenum, bleed_porosity, bleed_start, bound_offset, boundhiri, catalytic_con) &
!$omp& map(alloc: gamma_in, halo_offset, ibound_icode, ielem_bleedn, ielem_dih, ielem_dtl, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_ibounds, ielem_ifca, ielem_indexf, ielem_indexi, ielem_ineigh) &
!$omp& map(alloc: ielem_ineighb, ielem_ineighn, ielem_inter_id, ielem_interior, ielem_nodes_faces, ielem_qface) &
!$omp& map(alloc: ielem_reorient, ielem_surf, ielem_totvolume, ielem_types_faces, ielem_walldist, ielem_walltrans) &
!$omp& map(alloc: impdiag, impdiagt, impdu, impdu_old, impoff, impofft, inoder4_cord, lm_diag, lm_diagt, lm_off) &
!$omp& map(alloc: lm_offt, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_grads) &
!$omp& map(alloc: rec_ihexl, rec_ihexn, rec_local, rec_mrf, rec_rotvel, rg_hzero, rg_molm, rg_thetag, rg_vf, rgs_ab) &
!$omp& map(alloc: rgs_bb, rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, rhs_val, rhst_val, sht, solhir, u_c_val) &
!$omp& map(alloc: u_ct_val)
#else
!$omp do private(i, slot) &
!$omp& private(at)
#endif
         do at = chunk, chunk_end
            i = at
            if (ordered) i = sweep_order(at)
            slot = at - chunk + 1
            call implicit_sweep_cell(n, i, direction, low_storage, slot)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end do
   end subroutine implicit_sweep_range

   subroutine prepare_implicit_sweeps(n, ncells)
      implicit none
      integer, intent(in) :: n, ncells
      integer, allocatable :: degree(:), offset(:), edge(:), cursor(:), level(:), count(:)
      integer :: i, j, f, k, a, b, at
      if (sweep_cells == ncells) return
      if (allocated(sweep_order)) then
#ifdef xpu
!$omp target exit data map(delete: sweep_order)
#endif
         deallocate(sweep_order, sweep_offset)
      end if
      allocate(degree(ncells), offset(ncells + 1), cursor(ncells), level(ncells))
      degree = 0
      ! Orient every local adjacency from lower to higher cell index. Include
      ! both directions even if the mesh neighbour list is not reciprocal.
      do i = 1, ncells
         do f = 1, ielem_ifca(i)
            if (ielem_interior(i) /= 0) then
               if (ielem_ineighb(i, f) /= n) cycle
               if (ielem_ibounds(i, f) > 0) then
                  if (ibound_icode(ielem_ibounds(i, f)) /= 5 .and. ibound_icode(ielem_ibounds(i, f)) /= 50) cycle
               end if
            end if
            j = ielem_ineigh(i, f)
            if (j < 1 .or. j > ncells .or. j == i) cycle
            a = min(i, j)
            degree(a) = degree(a) + 1
         end do
      end do
      offset(1) = 1
      do i = 1, ncells
         offset(i + 1) = offset(i) + degree(i)
      end do
      allocate(edge(offset(ncells + 1) - 1))
      cursor = offset(1:ncells)
      do i = 1, ncells
         do f = 1, ielem_ifca(i)
            if (ielem_interior(i) /= 0) then
               if (ielem_ineighb(i, f) /= n) cycle
               if (ielem_ibounds(i, f) > 0) then
                  if (ibound_icode(ielem_ibounds(i, f)) /= 5 .and. ibound_icode(ielem_ibounds(i, f)) /= 50) cycle
               end if
            end if
            j = ielem_ineigh(i, f)
            if (j < 1 .or. j > ncells .or. j == i) cycle
            a = min(i, j); b = max(i, j)
            edge(cursor(a)) = b
            cursor(a) = cursor(a) + 1
         end do
      end do
      level = 1
      do i = 1, ncells
         do k = offset(i), offset(i + 1) - 1
            j = edge(k)
            level(j) = max(level(j), level(i) + 1)
         end do
      end do
      sweep_levels = 0
      if (ncells > 0) sweep_levels = maxval(level)
      allocate(count(sweep_levels), sweep_offset(sweep_levels + 1), sweep_order(ncells))
      count = 0
      do i = 1, ncells
         count(level(i)) = count(level(i)) + 1
      end do
      sweep_offset(1) = 1
      do k = 1, sweep_levels
         sweep_offset(k + 1) = sweep_offset(k) + count(k)
      end do
      count = 0
      do i = 1, ncells
         k = level(i); at = sweep_offset(k) + count(k)
         sweep_order(at) = i
         count(k) = count(k) + 1
      end do
      sweep_cells = ncells
#ifdef xpu
!$omp target enter data map(to: sweep_order)
#endif
   end subroutine prepare_implicit_sweeps

   subroutine relaxation(n, use_low_storage)
      ! Shared 2D/3D driver. Both backends execute the same cell routines.
      implicit none
      integer, intent(in) :: n
      logical, intent(in), optional :: use_low_storage
      logical :: low_storage
      integer :: ncells, i, sweep, lev, first, last
      low_storage = .false.
      if (present(use_low_storage)) low_storage = use_low_storage
      ncells = xmpielrank(n)
      if (low_storage) then
!$omp master
         call prepare_implicit_workspace(ncells)
!$omp end master
!$omp barrier
      end if
      call implicit_zero_impdu(ncells)
      if (.not. low_storage) then
         if (dimensiona == 3) then
            call calculate_jacobian(n)
         else
            call calculate_jacobian_2d(n)
         end if
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n, ncells) &
!$omp& private(i) &
!$omp& map(alloc: impdiag, impdiagt, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, ncells
            call implicit_prepare_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
      if (relax /= 1) then
!$omp master
         call prepare_implicit_sweeps(n, ncells)
!$omp end master
!$omp barrier
      end if
      do sweep = 1, 10
         if (relax == 1) then
            call exhboundhigher2(n)
            call implicit_snapshot_impdu(ncells)
            call implicit_sweep_range(n, 1, ncells, 0, low_storage, .false.)
         else
            ! Each level contains no neighbouring cells. Synchronous launches
            ! preserve ascending/descending Gauss-Seidel dependencies.
            do lev = 1, sweep_levels
               first = sweep_offset(lev); last = sweep_offset(lev + 1) - 1
               call implicit_sweep_range(n, first, last, 1, low_storage, .true.)
            end do
            call exhboundhigher2(n)
            do lev = sweep_levels, 1, -1
               first = sweep_offset(lev); last = sweep_offset(lev + 1) - 1
               call implicit_sweep_range(n, first, last, 2, low_storage, .true.)
            end do
            call exhboundhigher2(n)
         end if
      end do
   end subroutine relaxation

   subroutine implicit_prepare_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      if (dimensiona == 3) then
         if (rframe == 0) call relaxation_target_cell_1(n, i)
         if (srfg == 1) call relaxation_target_cell_2(n, i)
         if (mrf == 1) call relaxation_target_cell_3(n, i)
         if (turbulence > 0 .or. passivescalar > 0) call relaxation_target_cell_4(n, i)
      else
         call relaxation2d_target_cell_1(n, i)
         if (turbulence > 0 .or. passivescalar > 0) call relaxation2d_target_cell_2(n, i)
      end if
   end subroutine implicit_prepare_cell

   subroutine implicit_sweep_cell(n, i, direction, low_storage, slot)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i, direction, slot
      logical, intent(in) :: low_storage
      if (low_storage) then
         if (dimensiona == 3) then
            select case(direction)
            case(0)
               call relaxation_lm_target_cell_1(n, i, lm_diag(slot:slot, :, :), lm_diagt(slot:slot, :), &
                  & lm_off(slot:slot, :, :, :), lm_offt(slot:slot, :, :))
            case(1)
               call relaxation_lm_target_cell_2(n, i, lm_diag(slot:slot, :, :), lm_diagt(slot:slot, :), &
                  & lm_off(slot:slot, :, :, :), lm_offt(slot:slot, :, :))
            case(2)
               call relaxation_lm_target_cell_3(n, i, lm_diag(slot:slot, :, :), lm_diagt(slot:slot, :), &
                  & lm_off(slot:slot, :, :, :), lm_offt(slot:slot, :, :))
            end select
         else
            select case(direction)
            case(0)
               call relaxation_lm2d_target_cell_1(n, i, lm_diag(slot:slot, :, :), lm_diagt(slot:slot, :), &
                  & lm_off(slot:slot, :, :, :), lm_offt(slot:slot, :, :))
            case(1)
               call relaxation_lm2d_target_cell_2(n, i, lm_diag(slot:slot, :, :), lm_diagt(slot:slot, :), &
                  & lm_off(slot:slot, :, :, :), lm_offt(slot:slot, :, :))
            case(2)
               call relaxation_lm2d_target_cell_3(n, i, lm_diag(slot:slot, :, :), lm_diagt(slot:slot, :), &
                  & lm_off(slot:slot, :, :, :), lm_offt(slot:slot, :, :))
            end select
         end if
      elseif (dimensiona == 3) then
         select case(direction)
         case(0)
            call relaxation_target_cell_5(n, i)
         case(1)
            call relaxation_target_cell_6(n, i)
         case(2)
            call relaxation_target_cell_7(n, i)
         end select
      else
         select case(direction)
         case(0)
            call relaxation2d_target_cell_3(n, i)
         case(1)
            call relaxation2d_target_cell_4(n, i)
         case(2)
            call relaxation2d_target_cell_5(n, i)
         end select
      end if
   end subroutine implicit_sweep_cell

   logical function implicit_boundary_state_ok(q)
! Validate before applying a nonlinear EOS to a perturbed physical state.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: q(nof_variables)
      implicit_boundary_state_ok = conservative_state_valid(q)
   end function implicit_boundary_state_ok

   subroutine implicit_boundary_state(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, q, t, qg, tg, flag)
! Evaluate B(U) without modifying U or shared solver state.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, bcode, icell, iface
      real, intent(in) :: pox(dimensiona), poy(dimensiona), poz(dimensiona), angle1, angle2, nx, ny, nz
      real, intent(in) :: q(nof_variables), t(turbulenceequations + passivescalar)
      real, intent(out) :: qg(nof_variables), tg(turbulenceequations + passivescalar)
      integer, intent(out) :: flag
      real :: ql(gpu_max_nvar), tl(gpu_max_extra_transport)
      real :: cl(gpu_max_nvar_total), cr(gpu_max_nvar_total), speed(gpu_max_nvar), speedrot(gpu_max_nvar)
      ql = 0.0e0; ql(1:nof_variables) = q; tl = 0.0e0; tl(1:size(t)) = t
      qg = q; tg = t; cl = 0.0e0; cr = 0.0e0; speed = 0.0e0; speedrot = 0.0e0; flag = 0
      if (dimensiona .eq. 3) then
         if (rec_mrf(icell) .eq. 1) then
            speed(2:4) = rec_rotvel(icell, iface, 1, 1:3)
            call rotatef(n, speedrot, speed, angle1, angle2)
         end if
         call boundarys(n, bcode, icell, iface, ql, qg, pox, poy, poz, angle1, angle2, nx, ny, nz, tl, tg, cr, cl, speed, &
            & speedrot, flag)
      else
         call boundarys2d(n, bcode, icell, iface, ql, qg, pox, poy, poz, angle1, angle2, nx, ny, nz, tl, tg, cr, cl, speed, &
            & speedrot, flag)
      end if
   end subroutine implicit_boundary_state

   subroutine implicit_boundary_increment(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, dq, dg, dt, dgt, &
      & qghost, tghost)
! Directional linearization dUg = B'(U) dU, about the physical cell state.
! Never pass dU to a thermodynamic boundary routine as though it were U.
! Centered differences are used where admissible; at switches use the base branch.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, bcode, icell, iface
      real, intent(in) :: pox(dimensiona), poy(dimensiona), poz(dimensiona), angle1, angle2, nx, ny, nz
      real, intent(in) :: dq(nof_variables), dt(turbulenceequations + passivescalar)
      real, intent(out) :: dg(nof_variables), dgt(turbulenceequations + passivescalar)
      real, intent(out), optional :: qghost(nof_variables), tghost(turbulenceequations + passivescalar)
      real :: q0(gpu_max_nvar), qp(gpu_max_nvar), qm(gpu_max_nvar), g0(gpu_max_nvar), gp(gpu_max_nvar), gm(gpu_max_nvar)
      real :: t0(gpu_max_extra_transport), tp(gpu_max_extra_transport), tm(gpu_max_extra_transport)
      real :: gt0(gpu_max_extra_transport), gtp(gpu_max_extra_transport), gtm(gpu_max_extra_transport), scale(gpu_max_nvar)
      real :: h, normdir, rho, energy, coefficient
      integer :: nv, nt, flag0, flagp, flagm, attempt, iv
      logical :: plusok, minusok
      nv = nof_variables; nt = turbulenceequations + passivescalar
      dg = 0.0e0; dgt = 0.0e0
      q0 = 0.0e0; q0(1:nv) = u_c_val(icell, 1, 1:nv)
      t0 = 0.0e0
      if (nt .gt. 0) t0(1:nt) = u_ct_val(icell, 1, 1:nt)
      call implicit_boundary_state(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, q0, t0, g0, gt0, flag0)
      if (present(qghost)) qghost = g0(1:nv)
      if (present(tghost)) tghost = gt0(1:nt)
      if (all(dq .eq. 0.0e0) .and. all(dt .eq. 0.0e0)) return
      rho = max(abs(q0(1)), tolsmall); energy = max(abs(q0(dimensiona + 2)), rho)
      scale = rho; scale(2:dimensiona + 1) = sqrt(rho * energy); scale(dimensiona + 2) = energy
      if (realgas .eq. 1) scale(dimensiona + 3) = energy
      scale = max(scale, abs(q0))
      normdir = maxval(abs(dq) / scale(1:nv))
      if (nt .gt. 0) normdir = max(normdir, maxval(abs(dt) / max(abs(t0(1:nt)), rho)))
      if (normdir .lt. 1.0e-150) return
      h = 1.0e-6 / normdir
      plusok = .false.; minusok = .false.
      do attempt = 1, 12
         qp(1:nv) = q0(1:nv) + h * dq; qm(1:nv) = q0(1:nv) - h * dq
         tp(1:nt) = t0(1:nt) + h * dt; tm(1:nt) = t0(1:nt) - h * dt
         plusok = implicit_boundary_state_ok(qp(1:nv)); minusok = implicit_boundary_state_ok(qm(1:nv))
         if (plusok) then
            call implicit_boundary_state(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, qp, tp, gp, gtp, flagp)
            plusok = (flagp .eq. flag0)
         end if
         if (minusok) then
            call implicit_boundary_state(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, qm, tm, gm, gtm, flagm)
            minusok = (flagm .eq. flag0)
         end if
         if (plusok .or. minusok) exit
         h = 0.1e0 * h
      end do
      if (plusok .and. minusok) then
         dg = (gp(1:nv) - gm(1:nv)) / (2.0e0 * h); dgt = (gtp(1:nt) - gtm(1:nt)) / (2.0e0 * h)
      else if (plusok) then
         dg = (gp(1:nv) - g0(1:nv)) / h; dgt = (gtp(1:nt) - gt0(1:nt)) / h
      else if (minusok) then
         dg = (g0(1:nv) - gm(1:nv)) / h; dgt = (gt0(1:nt) - gtm(1:nt)) / h
      end if
! At a composition corner, a mixed direction can be inadmissible on BOTH sides.
! Differentiate each active coordinate one-sided, then form the linear combination.
      if (.not. plusok .and. .not. minusok) then
         do iv = 1, nv + nt
            if (iv .le. nv) then
               coefficient = dq(iv); h = 1.0e-6 * scale(iv)
            else
               coefficient = dt(iv - nv); h = 1.0e-6 * max(abs(t0(iv - nv)), rho)
            end if
            if (coefficient .eq. 0.0e0) cycle
            do attempt = 1, 12
               qp = q0; qm = q0; tp = t0; tm = t0
               if (iv .le. nv) then
                  qp(iv) = qp(iv) + h; qm(iv) = qm(iv) - h
               else
                  tp(iv - nv) = tp(iv - nv) + h; tm(iv - nv) = tm(iv - nv) - h
               end if
               plusok = implicit_boundary_state_ok(qp(1:nv)); minusok = implicit_boundary_state_ok(qm(1:nv))
               if (plusok) then
                  call implicit_boundary_state(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, qp, tp, gp, gtp, &
                     & flagp)
                  plusok = (flagp .eq. flag0)
               end if
               if (minusok) then
                  call implicit_boundary_state(n, bcode, icell, iface, pox, poy, poz, angle1, angle2, nx, ny, nz, qm, tm, gm, gtm, &
                     & flagm)
                  minusok = (flagm .eq. flag0)
               end if
               if (plusok .or. minusok) exit
               h = 0.1e0 * h
            end do
            if (plusok .and. minusok) then
               dg = dg + coefficient * (gp(1:nv) - gm(1:nv)) / (2.0e0 * h)
               dgt = dgt + coefficient * (gtp(1:nt) - gtm(1:nt)) / (2.0e0 * h)
            else if (plusok) then
               dg = dg + coefficient * (gp(1:nv) - g0(1:nv)) / h
               dgt = dgt + coefficient * (gtp(1:nt) - gt0(1:nt)) / h
            else if (minusok) then
               dg = dg + coefficient * (g0(1:nv) - gm(1:nv)) / h
               dgt = dgt + coefficient * (gt0(1:nt) - gtm(1:nt)) / h
            end if
         end do
      end if
   end subroutine implicit_boundary_increment

   subroutine invert_realgas_block(icell, block)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: icell
      real, intent(inout) :: block(:, :)
      real :: a(gpu_max_nvar, gpu_max_nvar)
      real :: scale(gpu_max_nvar), rowscale(gpu_max_nvar), row(gpu_max_nvar)
      real :: rho, energy, piv, factor
      integer :: i, j, k, p
      logical :: diagonal
! Rusanov contributes a scalar diagonal; active chemistry adds a coupled block.
! Retain scalar Rusanov inversion when chemistry adds no off-diagonal terms.
      diagonal = .true.
      do j = 1, nof_variables
         do i = 1, nof_variables
            if ((i .ne. j) .and. (block(i, j) .ne. 0.0e0)) diagonal = .false.
         end do
      end do
      if (diagonal) then
         do i = 1, nof_variables
            if (.not. (abs(block(i, i)) .gt. tiny(1.0e0))) error stop 'Singular real-gas implicit diagonal'
            block(i, i) = 1.0e0 / block(i, i)
         end do
         return
      end if
! Scale conservative variables before partial-pivoted Gauss-Jordan elimination.
      rho = max(u_c_val(icell, 1, 1), 1.0e-30)
      energy = max(abs(u_c_val(icell, 1, dimensiona + 2)), rho)
      scale = rho
      scale(2:dimensiona + 1) = sqrt(rho * energy)
      scale(dimensiona + 2:min(dimensiona + 3, nof_variables)) = energy
      a = 0.0e0
      do i = 1, nof_variables
         a(i, 1:nof_variables) = block(i, 1:nof_variables) * scale(1:nof_variables) / scale(i)
         rowscale(i) = maxval(abs(a(i, 1:nof_variables)))
         if (.not. (rowscale(i) .gt. tiny(1.0e0))) error stop 'Singular real-gas implicit block: zero/nonfinite row'
         a(i, 1:nof_variables) = a(i, 1:nof_variables) / rowscale(i)
         block(i, 1:nof_variables) = 0.0e0
         block(i, i) = 1.0e0
      end do
      do k = 1, nof_variables
         p = k
         do i = k + 1, nof_variables
            if (abs(a(i, k)) .gt. abs(a(p, k))) p = i
         end do
         if (.not. (abs(a(p, k)) .gt. tiny(1.0e0))) error stop 'Singular real-gas implicit block: zero/nonfinite pivot'
         if (p .ne. k) then
            row = a(k, :); a(k, :) = a(p, :); a(p, :) = row
            row(1:nof_variables) = block(k, 1:nof_variables)
            block(k, 1:nof_variables) = block(p, 1:nof_variables)
            block(p, 1:nof_variables) = row(1:nof_variables)
         end if
         piv = a(k, k)
         a(k, :) = a(k, :) / piv; block(k, 1:nof_variables) = block(k, 1:nof_variables) / piv
         do i = 1, nof_variables
            if (i .eq. k) cycle
            factor = a(i, k)
            a(i, :) = a(i, :) - factor * a(k, :)
            block(i, 1:nof_variables) = block(i, 1:nof_variables) - factor * block(k, 1:nof_variables)
         end do
      end do
      do j = 1, nof_variables
         block(1:nof_variables, j) = scale(1:nof_variables) * block(1:nof_variables, j) / (rowscale(j) * scale(j))
      end do
   end subroutine invert_realgas_block
   subroutine implicit_zero_impdu(kmaxe)
      implicit none
      integer, intent(in) :: kmaxe
      integer :: i
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(kmaxe) private(i) &
!$omp& firstprivate(nof_variables, passivescalar, turbulenceequations, zero) &
!$omp& map(alloc: impdu)
#else
!$omp do &
!$omp& private(i)
#endif
      do i = 1, kmaxe
         impdu(i, 1:nof_variables + turbulenceequations + passivescalar) = zero
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine implicit_zero_impdu

   subroutine implicit_snapshot_impdu(kmaxe)
      implicit none
      integer, intent(in) :: kmaxe
      integer :: i
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(kmaxe) private(i) &
!$omp& firstprivate(nof_variables, passivescalar, turbulenceequations) &
!$omp& map(alloc: impdu, impdu_old)
#else
!$omp do &
!$omp& private(i)
#endif
      do i = 1, kmaxe
         impdu_old(i, 1:nof_variables + turbulenceequations + passivescalar) = &
            & impdu(i, 1:nof_variables + turbulenceequations + passivescalar)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine implicit_snapshot_impdu

   subroutine relaxation_target_cell_1(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call invert_implicit_block(i, impdiag(i, :, :))
   end subroutine relaxation_target_cell_1

   subroutine relaxation_target_cell_2(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call invert_implicit_block(i, impdiag(i, :, :))
   end subroutine relaxation_target_cell_2

   subroutine relaxation_target_cell_3(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call invert_implicit_block(i, impdiag(i, :, :))
   end subroutine relaxation_target_cell_3

   subroutine relaxation_target_cell_4(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      impdiagt(i, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(i, 1:turbulenceequations + passivescalar), 1.0e-30)
   end subroutine relaxation_target_cell_4

   subroutine relaxation_target_cell_5(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: n_node
      integer :: ibfc
      integer :: srf
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00

      if (mrf .eq. 1) then
         srf = rec_mrf(i)
      end if
      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))
            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)
         end if
      end if
      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            du1(1:nof_variables) = zero
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               dut1(:) = zero
            end if

            du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                  & turbulenceequations + passivescalar)
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) &
                  & - impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = (cos(angle1) * sin(angle2))
            ny = (sin(angle1) * sin(angle2))
            nz = (cos(angle2))
            if (rec_mrf(i) .eq. 1) then
               !retrieve rotational velocity in case of rotating reference frame to calculate the correct value of the boundary condition
               srf_speed(2:4) = rec_rotvel(i, l, 1, 1:3)
               call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
            end if

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  !periodic in my cpu
                  if ((ibound_icode(ielem_ibounds(i, l)) .eq. 5) .or. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then
                     du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                     if ((per_rot .eq. 1) .and. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then
                        du1(2:4) = rotate_per_1(du1(2:4), ibound_icode(ielem_ibounds(i, l)), angle_per)
                     end if
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), &
                           & nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

                     end if

                  else
                     !not periodic ones in my cpu

                     facex = l; iconsidered = i
                     call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)

                     if (ielem_types_faces(iconsidered, facex) .eq. 5) then
                        n_node = 4
                     else
                        n_node = 3
                     end if

                     cords(1:3) = zero
                     call cordinates3(n, nodes_list, n_node, cords(1:3))

                     poy(1) = cords(2)
                     pox(1) = cords(1)
                     poz(1) = cords(3)

                     leftv(1:nof_variables) = impdu_old(i, 1:nof_variables)
                     cturbl = 0.0e0
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                        & cturbl(1:turbulenceequations + passivescalar) = impdu_old(i, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)
                     b_code = ibound_icode(ielem_ibounds(i, l))
                     call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & leftv, du1, cturbl, dut1)

                  end if
               else
                  du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                     dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else !in other cpus they can only be periodic or mpi neighbours

               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  !periodic in other cpu
                  if ((ibound_icode(ielem_ibounds(i, l)) .eq. 5) .or. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)
                     if ((per_rot .eq. 1) .and. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then
                        du1(2:4) = rotate_per_1(du1(2:4), ibound_icode(ielem_ibounds(i, l)), angle_per)
                     end if
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                  nfx = ielem_ineighn(i, l)
                  lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                  rowfx = bound_offset(nfx) + lfx - 1
                  du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

!
               end if
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) &
                  & - impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(i, mm_i, mm_j) * b1_imp(mm_j)
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(i, &
            & 1:turbulenceequations + passivescalar) * b1t(1:turbulenceequations + passivescalar)
      end if
   end subroutine relaxation_target_cell_5

   subroutine relaxation_target_cell_6(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: igoflux
      integer :: icaseb
      integer :: n_node
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      dummy12(:) = zero
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         dummy12t(:) = zero
      end if

      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))
! dummy12(:)=zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))

            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)
! dummy12(:)=zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)
! dummy12t(:)=zero
         end if
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 0) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            igoflux = 0
            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     icaseb = 1 !periodic mine
                  else
                     icaseb = 3 !physical
                  end if

               else

                  icaseb = 2 !no boundaries interior

               end if
            else
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu
                     icaseb = 4
                  end if
               else
                  icaseb = 5

               end if

            end if
            if (icaseb .le. 2) then
               if (ielem_reorient(i, l) .eq. 0) then
                  igoflux = 1
               else
                  igoflux = 0
               end if
            else
               igoflux = 2
            end if
            if (ielem_reorient(i, l) .eq. 0) then
               srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else
                        !not periodic ones in my cpu

                        facex = l; iconsidered = i
                        call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)

                        if (ielem_types_faces(iconsidered, facex) .eq. 5) then
                           n_node = 4
                        else
                           n_node = 3
                        end if

                        cords(1:3) = zero
                        call cordinates3(n, nodes_list, n_node, cords(1:3))

                        poy(1) = cords(2)
                        pox(1) = cords(1)
                        poz(1) = cords(3)

                        leftv(1:nof_variables) = impdu(i, 1:nof_variables)
                        cturbl = 0.0e0
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                           & cturbl(1:turbulenceequations + passivescalar) = impdu(i, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)
                        b_code = ibound_icode(ielem_ibounds(i, l))
                        call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & leftv, du1, cturbl, dut1)

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(i, mm_i, mm_j) * (b1_imp(mm_j) - dummy12(mm_j))
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(i, &
            & 1:turbulenceequations + passivescalar) * &
            & (b1t(1:turbulenceequations + passivescalar) - dummy12t(1:turbulenceequations + passivescalar))

      end if
   end subroutine relaxation_target_cell_6

   subroutine relaxation_target_cell_7(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: igoflux
      integer :: icaseb
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      dummy12(:) = zero
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         dummy12t(:) = zero
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 1) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3

            igoflux = 0
            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     icaseb = 1 !periodic mine
                  else
                     icaseb = 3 !physical
                  end if

               else

                  icaseb = 2 !no boundaries interior

               end if
            else
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu
                     icaseb = 4
                  end if
               else
                  icaseb = 5

               end if

            end if
            if (icaseb .le. 2) then
               if (ielem_reorient(i, l) .eq. 1) then
                  igoflux = 1
               else
                  igoflux = 0
               end if
            else
               igoflux = 2
            end if
            if (ielem_reorient(i, l) .eq. 1) then

               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if

      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(i, mm_i, mm_j) * dummy12(mm_j)
         end do
         impdu(i, mm_i) = impdu(i, mm_i) - mm_sum
      end do
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdu(i, &
            & nof_variables + 1:nof_variables + turbulenceequations + passivescalar) - &
            & (impdiagt(i, 1:turbulenceequations + passivescalar) * dummy12t(1:turbulenceequations + passivescalar))

      end if
   end subroutine relaxation_target_cell_7

   subroutine relaxation_lm(n)
      implicit none
      integer, intent(in) :: n
      call relaxation(n, .true.)
   end subroutine relaxation_lm

   subroutine relaxation_lm_target_cell_1(n, i, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: impdiag(:, :, :), impdiagt(:, :), impoff(:, :, :, :), impofft(:, :, :)
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: n_node
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      iconsidered = i
      call calculate_jacobianlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      call invert_implicit_block(i, impdiag(1, :, :))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdiagt(1, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(1, 1:turbulenceequations + passivescalar), &
            & 1.0e-30)
      end if

      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))

         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))

            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)

         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)

         end if
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            du1(1:nof_variables) = zero
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               dut1(:) = zero
            end if

            du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                  & turbulenceequations + passivescalar)
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) - &
                  & impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = (cos(angle1) * sin(angle2))
            ny = (sin(angle1) * sin(angle2))
            nz = (cos(angle2))

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), &
                           & nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

                     end if

                  else
                     !not periodic ones in my cpu

                     facex = l; iconsidered = i
                     call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)

                     if (ielem_types_faces(iconsidered, facex) .eq. 5) then
                        n_node = 4
                     else
                        n_node = 3
                     end if

                     cords(1:3) = zero
                     call cordinates3(n, nodes_list, n_node, cords(1:3))

                     poy(1) = cords(2)
                     pox(1) = cords(1)
                     poz(1) = cords(3)

                     leftv(1:nof_variables) = impdu_old(i, 1:nof_variables)
                     cturbl = 0.0e0
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                        & cturbl(1:turbulenceequations + passivescalar) = impdu_old(i, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)
                     b_code = ibound_icode(ielem_ibounds(i, l))
                     call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & leftv, du1, cturbl, dut1)

                  end if
               else
                  du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                     dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else !in other cpus they can only be periodic or mpi neighbours

               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                  nfx = ielem_ineighn(i, l)
                  lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                  rowfx = bound_offset(nfx) + lfx - 1
                  du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

!
               end if
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) - &
                  & impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f
      end if

      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(1, mm_i, mm_j) * b1_imp(mm_j)
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(1, &
            & 1:turbulenceequations + passivescalar) * b1t(1:turbulenceequations + passivescalar)
      end if
   end subroutine relaxation_lm_target_cell_1

   subroutine relaxation_lm_target_cell_2(n, i, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: impdiag(:, :, :), impdiagt(:, :), impoff(:, :, :, :), impofft(:, :, :)
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: n_node
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      iconsidered = i
      call calculate_jacobianlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      call invert_implicit_block(i, impdiag(1, :, :))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdiagt(1, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(1, 1:turbulenceequations + passivescalar), &
            & 1.0e-30)
      end if

      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))
         dummy12(:) = zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))
               dummy12t(:) = zero
            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)
         dummy12(:) = zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)
            dummy12t(:) = zero
         end if
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 0) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 0) then
               srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else
                        !not periodic ones in my cpu

                        facex = l; iconsidered = i
                        call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)

                        if (ielem_types_faces(iconsidered, facex) .eq. 5) then
                           n_node = 4
                        else
                           n_node = 3
                        end if

                        cords(1:3) = zero
                        call cordinates3(n, nodes_list, n_node, cords(1:3))

                        poy(1) = cords(2)
                        pox(1) = cords(1)
                        poz(1) = cords(3)

                        leftv(1:nof_variables) = impdu(i, 1:nof_variables)
                        cturbl = 0.0e0
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                           & cturbl(1:turbulenceequations + passivescalar) = impdu(i, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)
                        b_code = ibound_icode(ielem_ibounds(i, l))
                        call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & leftv, du1, cturbl, dut1)

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if

      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(1, mm_i, mm_j) * (b1_imp(mm_j) - dummy12(mm_j))
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(1, &
            & 1:turbulenceequations + passivescalar) * &
            & (b1t(1:turbulenceequations + passivescalar) - dummy12t(1:turbulenceequations + passivescalar))

      end if
   end subroutine relaxation_lm_target_cell_2

   subroutine relaxation_lm_target_cell_3(n, i, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: impdiag(:, :, :), impdiagt(:, :), impoff(:, :, :, :), impofft(:, :, :)
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00

      iconsidered = i
      call calculate_jacobianlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      call invert_implicit_block(i, impdiag(1, :, :))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdiagt(1, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(1, 1:turbulenceequations + passivescalar), &
            & 1.0e-30)
      end if

      dummy12(:) = zero
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         dummy12t(:) = zero
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 1) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 1) then
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(1, mm_i, mm_j) * dummy12(mm_j)
         end do
         impdu(i, mm_i) = impdu(i, mm_i) - mm_sum
      end do
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdu(i, &
            & nof_variables + 1:nof_variables + turbulenceequations + passivescalar) &
            & - ((impdiagt(1, 1:turbulenceequations + passivescalar) * dummy12t(1:turbulenceequations + passivescalar)))

      end if
   end subroutine relaxation_lm_target_cell_3

   subroutine relaxation2d(n)
      implicit none
      integer, intent(in) :: n
      call relaxation(n, .false.)
   end subroutine relaxation2d

   subroutine relaxation2d_target_cell_1(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call invert_implicit_block(i, impdiag(i, :, :))
   end subroutine relaxation2d_target_cell_1

   subroutine relaxation2d_target_cell_2(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      impdiagt(i, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(i, 1:turbulenceequations + passivescalar), 1.0e-30)
   end subroutine relaxation2d_target_cell_2

   subroutine relaxation2d_target_cell_3(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: n_node
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))
            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)
         end if
      end if
      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            du1(1:nof_variables) = zero
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               dut1(:) = zero
            end if

            du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                  & turbulenceequations + passivescalar)
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) &
                  & - impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = angle1
            ny = angle2

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), &
                           & nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

                     end if

                  else
                     !not periodic ones in my cpu

                     facex = l; iconsidered = i
                     call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                     n_node = 2
                     cords(1:2) = zero
                     call cordinates2(n, nodes_list, n_node, cords(1:2))

                     poy(1) = cords(2)
                     pox(1) = cords(1)

                     leftv(1:nof_variables) = impdu_old(i, 1:nof_variables)
                     cturbl = 0.0e0
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                        & cturbl(1:turbulenceequations + passivescalar) = impdu_old(i, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)
                     b_code = ibound_icode(ielem_ibounds(i, l))
                     call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & leftv, du1, cturbl, dut1)

                  end if
               else
                  du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                     dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else !in other cpus they can only be periodic or mpi neighbours

               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                  nfx = ielem_ineighn(i, l)
                  lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                  rowfx = bound_offset(nfx) + lfx - 1
                  du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

!
               end if
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) &
                  & - impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(i, mm_i, mm_j) * b1_imp(mm_j)
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(i, &
            & 1:turbulenceequations + passivescalar) * b1t(1:turbulenceequations + passivescalar)
      end if
      !loop elements

   end subroutine relaxation2d_target_cell_3

   subroutine relaxation2d_target_cell_4(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: igoflux
      integer :: icaseb
      integer :: n_node
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      dummy12(:) = zero
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         dummy12t(:) = zero
      end if

      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))

         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))
! dummy12t(:)=zero
            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)
! dummy12(:)=zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)
! dummy12t(:)=zero
         end if
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_ineigh(i, l) .lt. i) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            igoflux = 0
            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     icaseb = 1 !periodic mine
                  else
                     icaseb = 3 !physical
                  end if

               else

                  icaseb = 2 !no boundaries interior

               end if
            else
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu
                     icaseb = 4
                  end if
               else
                  icaseb = 5

               end if

            end if
            if (icaseb .le. 2) then
               if (ielem_ineigh(i, l) .lt. i) then
                  igoflux = 1
               else
                  igoflux = 0
               end if
            else
               igoflux = 2
            end if
            if (ielem_ineigh(i, l) .lt. i) then
               srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = angle1
               ny = angle2

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else
                        !not periodic ones in my cpu

                        facex = l; iconsidered = i
                        call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                        n_node = 2
                        cords(1:2) = zero
                        call cordinates2(n, nodes_list, n_node, cords(1:2))

                        poy(1) = cords(2)
                        pox(1) = cords(1)

                        leftv(1:nof_variables) = impdu(i, 1:nof_variables)
                        cturbl = 0.0e0
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                           & cturbl(1:turbulenceequations + passivescalar) = impdu(i, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)
                        b_code = ibound_icode(ielem_ibounds(i, l))
                        call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & leftv, du1, cturbl, dut1)

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(i, mm_i, mm_j) * (b1_imp(mm_j) - dummy12(mm_j))
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(i, &
            & 1:turbulenceequations + passivescalar) * &
            & (b1t(1:turbulenceequations + passivescalar) - dummy12t(1:turbulenceequations + passivescalar))

      end if
   end subroutine relaxation2d_target_cell_4

   subroutine relaxation2d_target_cell_5(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: igoflux
      integer :: icaseb
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      dummy12(:) = zero
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         dummy12t(:) = zero
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_ineigh(i, l) .gt. i) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            igoflux = 0

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     icaseb = 1 !periodic mine
                  else
                     icaseb = 3 !physical
                  end if

               else

                  icaseb = 2 !no boundaries interior

               end if
            else
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu
                     icaseb = 4
                  end if
               else
                  icaseb = 5

               end if

            end if
            if (icaseb .le. 2) then
               if (ielem_reorient(i, l) .eq. 1) then
                  igoflux = 1
               else
                  igoflux = 0
               end if
            else
               igoflux = 2
            end if
            if (ielem_ineigh(i, l) .gt. i) then
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = angle1
               ny = angle2

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                                                                 dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                                                                             du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                                                                 dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                     !
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(i, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(i, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if

! impdu(i,1:nof_variables)=explicit_matrix_product(impdiag(i,1:nof_variables,1:nof_variables),impdu(i,1:nof_variables)-dummy12(1:nof_variables))

      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(i, mm_i, mm_j) * dummy12(mm_j)
         end do
         impdu(i, mm_i) = impdu(i, mm_i) - mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdu(i, &
            & nof_variables + 1:nof_variables + turbulenceequations + passivescalar) - &
            & (impdiagt(i, 1:turbulenceequations + passivescalar) * dummy12t(1:turbulenceequations + passivescalar))

      end if
   end subroutine relaxation2d_target_cell_5

   subroutine relaxation_lm2d(n)
      implicit none
      integer, intent(in) :: n
      call relaxation(n, .true.)
   end subroutine relaxation_lm2d

   subroutine relaxation_lm2d_target_cell_1(n, i, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: impdiag(:, :, :), impdiagt(:, :), impoff(:, :, :, :), impofft(:, :, :)
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: n_node
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      iconsidered = i
      call calculate_jacobian_2dlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      call invert_implicit_block(i, impdiag(1, :, :))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdiagt(1, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(1, 1:turbulenceequations + passivescalar), &
            & 1.0e-30)
      end if

      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))

         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))

            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)

         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)

         end if
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            du1(1:nof_variables) = zero
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               dut1(:) = zero
            end if

            du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                  & turbulenceequations + passivescalar)
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) - &
                  & impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = angle1
            ny = angle2

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), &
                           & nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

                     end if

                  else
                     !not periodic ones in my cpu

                     facex = l; iconsidered = i
                     call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                     n_node = 2
                     cords(1:2) = zero
                     call cordinates2(n, nodes_list, n_node, cords(1:2))

                     poy(1) = cords(2)
                     pox(1) = cords(1)

                     leftv(1:nof_variables) = impdu_old(i, 1:nof_variables)
                     cturbl = 0.0e0
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                        & cturbl(1:turbulenceequations + passivescalar) = impdu_old(i, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)
                     b_code = ibound_icode(ielem_ibounds(i, l))
                     call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & leftv, du1, cturbl, dut1)

                  end if
               else
                  du1(1:nof_variables) = impdu_old(ielem_ineigh(i, l), 1:nof_variables)
                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                     dut1(1:turbulenceequations + passivescalar) = impdu_old(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else !in other cpus they can only be periodic or mpi neighbours

               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                  nfx = ielem_ineighn(i, l)
                  lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                  rowfx = bound_offset(nfx) + lfx - 1
                  du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                  if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

!
               end if
            end if

            do mm_i = 1, nof_variables
               mm_sum = zero
               do mm_j = 1, nof_variables
                  mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
               end do
               b1_imp(mm_i) = b1_imp(mm_i) - mm_sum
            end do
            if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
               b1t(1:turbulenceequations + passivescalar) = b1t(1:turbulenceequations + passivescalar) - &
                  & impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar)
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(1, mm_i, mm_j) * b1_imp(mm_j)
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(1, &
            & 1:turbulenceequations + passivescalar) * b1t(1:turbulenceequations + passivescalar)
      end if
   end subroutine relaxation_lm2d_target_cell_1

   subroutine relaxation_lm2d_target_cell_2(n, i, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: impdiag(:, :, :), impdiagt(:, :), impoff(:, :, :, :), impofft(:, :, :)
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      integer :: nvar
      integer :: ibfc
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: b1t(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: facex
      integer :: n_node
      integer :: b_code
      integer :: nfx
      integer :: lfx
      integer :: rowfx
      real :: angle1
      real :: angle2
      real :: nx
      real :: ny
      real :: nz
      real, dimension(1:RELAX_NTOTAL) :: cright_rot
      real, dimension(1:RELAX_NTOTAL) :: cleft_rot
      real, dimension(1:RELAX_EXTRA) :: cturbl
      real, dimension(1:RELAX_EXTRA) :: cturbr
      real, dimension(1:RELAX_NVAR) :: leftv
      real, dimension(1:RELAX_NVAR) :: srf_speedrot
      real, dimension(1:RELAX_NVAR) :: srf_speed
      real, dimension(1:RELAX_NVAR) :: rightv
      real, dimension(1:RELAX_DIM) :: pox
      real, dimension(1:RELAX_DIM) :: poy
      real, dimension(1:RELAX_DIM) :: poz
      real, dimension(1:8, 1:RELAX_DIM) :: vext
      real, dimension(1:8, 1:RELAX_DIM) :: nodes_list
      real, dimension(3) :: cords

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00
      iconsidered = i
      call calculate_jacobian_2dlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      call invert_implicit_block(i, impdiag(1, :, :))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdiagt(1, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(1, 1:turbulenceequations + passivescalar), &
            & 1.0e-30)
      end if

      if (iscoun .ne. 1) then
         b1_imp(1:nof_variables) = -(rhs_val(i, 1:nof_variables) + ((((1.5 * u_c_val(i, 1, 1:nof_variables)) - (2.00 * u_c_val(i, &
            & 2, 1:nof_variables)) + (0.50 * u_c_val(i, 3, 1:nof_variables))) / (dt)) * ielem_totvolume(i)))
         dummy12(:) = zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do nvar = 1, turbulenceequations + passivescalar
               b1t(nvar) = -(rhst_val(i, nvar) + ((((1.50 * u_ct_val(i, 1, nvar)) - (2.00 * u_ct_val(i, 2, &
                  & nvar)) + (0.50 * u_ct_val(i, 3, nvar))) / (dt)) * ielem_totvolume(i)))
               dummy12t(:) = zero
            end do
         end if
      else
         b1_imp(1:nof_variables) = -rhs_val(i, 1:nof_variables)
         dummy12(:) = zero
         if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            b1t(1:turbulenceequations + passivescalar) = -rhst_val(i, 1:turbulenceequations + passivescalar)
            dummy12t(:) = zero
         end if
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 0) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 0) then
               srf_speed = 0.0e0; srf_speedrot = 0.0e0; pox = 0.0e0; poy = 0.0e0; poz = 0.0e0
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = angle1
               ny = angle2

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else
                        !not periodic ones in my cpu

                        facex = l; iconsidered = i
                        call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                        n_node = 2
                        cords(1:2) = zero
                        call cordinates2(n, nodes_list, n_node, cords(1:2))

                        poy(1) = cords(2)
                        pox(1) = cords(1)

                        leftv(1:nof_variables) = impdu(i, 1:nof_variables)
                        cturbl = 0.0e0
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) &
                           & cturbl(1:turbulenceequations + passivescalar) = impdu(i, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)
                        b_code = ibound_icode(ielem_ibounds(i, l))
                        call implicit_boundary_increment(n, b_code, iconsidered, facex, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & leftv, du1, cturbl, dut1)

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) + &
                     & (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if

      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(1, mm_i, mm_j) * (b1_imp(mm_j) - dummy12(mm_j))
         end do
         impdu(i, mm_i) = mm_sum
      end do

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdiagt(1, &
            & 1:turbulenceequations + passivescalar) * &
            & (b1t(1:turbulenceequations + passivescalar) - dummy12t(1:turbulenceequations + passivescalar))

      end if
   end subroutine relaxation_lm2d_target_cell_2

   subroutine relaxation_lm2d_target_cell_3(n, i, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(inout) :: impdiag(:, :, :), impdiagt(:, :), impoff(:, :, :, :), impofft(:, :, :)
      integer, intent(in) :: n, i
      integer :: mm_i
      integer :: mm_j
      real :: mm_sum
      integer :: l
      real :: b1_imp(1:RELAX_NVAR)
      real :: du1(1:RELAX_NVAR)
      real :: dummy12(1:RELAX_NVAR)
      real :: dur(RELAX_NVAR)
      real :: dul(RELAX_NVAR)
      real :: durr(RELAX_NVAR)
      real :: dulr(RELAX_NVAR)
      real :: dut1(RELAX_EXTRA)
      real :: dummy12t(RELAX_EXTRA)
      integer :: iconsidered
      integer :: nfx
      integer :: lfx
      integer :: rowfx

      du1 = 0.00; b1_imp = 0.00; dur = 0.00; dul = 0.00; durr = 0.00; dulr = 0.00

      iconsidered = i
      call calculate_jacobian_2dlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      call invert_implicit_block(i, impdiag(1, :, :))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         impdiagt(1, 1:turbulenceequations + passivescalar) = 1.00 / max(impdiagt(1, 1:turbulenceequations + passivescalar), &
            & 1.0e-30)
      end if

      dummy12(:) = zero
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         dummy12t(:) = zero
      end if

      if (ielem_interior(i) .eq. 0) then
         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 1) then
               du1(1:nof_variables) = zero
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dut1(:) = zero
               end if

               du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f

      else

         do l = 1, ielem_ifca(i) !loop3
            if (ielem_reorient(i, l) .eq. 1) then

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                           dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables &
                              & + turbulenceequations + passivescalar)

                        end if

                     else

                     end if
                  else
                     du1(1:nof_variables) = impdu(ielem_ineigh(i, l), 1:nof_variables)
                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

                        dut1(1:turbulenceequations + passivescalar) = impdu(ielem_ineigh(i, l), nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                 du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                        if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                           rowfx = bound_offset(nfx) + lfx - 1
                           dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        end if

                     end if
                  else
!                                                               du1(1:nof_variables)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                     rowfx = bound_offset(nfx) + lfx - 1
                     du1(1:nof_variables) = boundhiri(rowfx, 1:nof_variables)

                     if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

!                                                                   dut1(1:turbulenceequations+passivescalar)=iexboundhiri(ielem_ineighn(i,l))%facesol(ielem_qface(l,1,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, 1)
                        rowfx = bound_offset(nfx) + lfx - 1
                        dut1(1:turbulenceequations + passivescalar) = boundhiri(rowfx, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

!
                  end if
               end if

               do mm_i = 1, nof_variables
                  mm_sum = zero
                  do mm_j = 1, nof_variables
                     mm_sum = mm_sum + impoff(1, l, mm_i, mm_j) * du1(mm_j)
                  end do
                  dummy12(mm_i) = dummy12(mm_i) + mm_sum
               end do
               if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
                  dummy12t(1:turbulenceequations + passivescalar) = dummy12t(1:turbulenceequations + passivescalar) &
                     & + (impofft(1, l, 1:turbulenceequations + passivescalar) * dut1(1:turbulenceequations + passivescalar))
               end if
            end if
         end do !loop f
      end if
      do mm_i = 1, nof_variables
         mm_sum = zero
         do mm_j = 1, nof_variables
            mm_sum = mm_sum + impdiag(1, mm_i, mm_j) * dummy12(mm_j)
         end do
         impdu(i, mm_i) = impdu(i, mm_i) - mm_sum
      end do
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         impdu(i, nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = impdu(i, &
            & nof_variables + 1:nof_variables + turbulenceequations + passivescalar) &
            & - ((impdiagt(1, 1:turbulenceequations + passivescalar) * dummy12t(1:turbulenceequations + passivescalar)))

      end if
   end subroutine relaxation_lm2d_target_cell_3


end module implicit_time
