module implicit_fluxes
   use declaration
   use library
   use transform
   use local
   use riemann
   use flow_operations
   use source
   implicit none

contains
   subroutine add_duct_outlet_jacobian(n, i, qcell, normal, radius, area, block, off, diagt, offt)
! Replace the Rusanov boundary approximation by dF(B(U))/dU in the cell block.
! The interior .5*A(U)*n terms were cancelled using the closed-cell area sum;
! restore the missing outlet term as -.5*A(U)*n. Keep existing viscous damping.
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      real, intent(in) :: qcell(nof_variables), normal(3), radius, area
      real, intent(inout) :: block(:, :), off(:, :)
      real, intent(inout), optional :: diagt(:), offt(:)
      real :: jac(gpu_max_nvar, gpu_max_nvar), euler(gpu_max_nvar, gpu_max_nvar)
      real :: transport(gpu_max_extra_transport), speed(gpu_max_nvar)
      integer :: j, nv, nt
      nv = nof_variables; nt = turbulenceequations + passivescalar
      call pressure_outlet_face_jacobian(n, qcell, normal, press_outlet, jac, transport)
      speed = 0.0e0; euler = 0.0e0
      if (dimensiona == 3) then
         call compute_jacobianse(n, i, euler, qcell, gamma, 0.0e0, 0.0e0, speed, normal(1), normal(2), normal(3))
      else
         call compute_jacobianse2d(n, euler, qcell, gamma, 0.0e0, 0.0e0, normal(1), normal(2), 0.0e0)
      end if
      block(1:nv, 1:nv) = block(1:nv, 1:nv) + area * (jac(1:nv, 1:nv) - 0.5e0 * euler(1:nv, 1:nv))
      do j = 1, nv
         block(j, j) = block(j, j) - 0.5e0 * radius * area
      end do
      off(1:nv, 1:nv) = 0.0e0
      if (nt > 0) then
         diagt(1:nt) = diagt(1:nt) + area * (transport(1:nt) - 0.5e0 * radius)
         offt(1:nt) = 0.0e0
      end if
   end subroutine add_duct_outlet_jacobian

   subroutine add_time_and_realgas_source(n, i, block)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      real, intent(inout) :: block(:, :)
      real :: jac(gpu_max_nvar, gpu_max_nvar), time_diagonal
      real :: qcell(gpu_max_nvar)
      integer :: j
      time_diagonal = ielem_totvolume(i) / ielem_dtl(i)
      if (rungekutta .ne. 10) time_diagonal = time_diagonal + 1.5e0 * ielem_totvolume(i) / dt
      do j = 1, nof_variables
         block(j, j) = block(j, j) + time_diagonal
      end do
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2)) then
         qcell(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         call realgas_source_jacobian(n, qcell, jac)
         block(1:nof_variables, 1:nof_variables) = block(1:nof_variables, 1:nof_variables) &
            & - ielem_totvolume(i) * jac(1:nof_variables, 1:nof_variables)
      end if
   end subroutine add_time_and_realgas_source
   subroutine calculate_jacobian(n)
      implicit none
      !> @brief
      !> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      integer :: i, kmaxe, ii

      kmaxe = xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_interior) &
!$omp& map(alloc: el_int, gamma_in, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_ineigh) &
!$omp& map(alloc: ielem_surf, ielem_walldist, impdiag, impdiagt, impoff, impofft, mp_janaf, mp_m, mp_pinf) &
!$omp& map(alloc: mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_grads, rec_mrf, rec_rotvel, rg_hzero, rg_molm, rg_thetag) &
!$omp& map(alloc: rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_interior !for all the interior elements
         call calculate_jacobian_inner_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_bounded) &
!$omp& map(alloc: catalytic_con, el_bnd, gamma_in, halo_offset, ibound_icode, ielem_dih, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_ibounds, ielem_ifca, ielem_indexi, ielem_ineigh, ielem_ineighb) &
!$omp& map(alloc: ielem_nodes_faces, ielem_surf, ielem_types_faces, ielem_walldist, impdiag, impdiagt, impoff) &
!$omp& map(alloc: impofft, inoder4_cord, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rec_grads, rec_ihexl, rec_ihexn, rec_local, rec_mrf, rec_rotvel, rg_hzero, rg_molm) &
!$omp& map(alloc: rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, solhir, u_c_val) &
!$omp& map(alloc: u_ct_val)
#else
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_bounded
         call calculate_jacobian_bound_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

      !add the contribution of the source term to the jacobian of the diagonal matrix
      if (srfg .eq. 1) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_totvolume, impdiag, srf_velocity)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call calculate_jacobian_loop1_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
      if (mrf .eq. 1) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_totvolume, impdiag, rec_mrf, rec_mrf_velocity)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call calculate_jacobian_loop2_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if

      if (rungekutta .eq. 10) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: gamma_in, ielem_dtl, ielem_totvolume, impdiag, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rg_hzero, rg_molm, rg_thetag, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call calculate_jacobian_loop3_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      else
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: gamma_in, ielem_dtl, ielem_totvolume, impdiag, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rg_hzero, rg_molm, rg_thetag, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
!
         do i = 1, kmaxe
            call calculate_jacobian_loop4_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         if (turbulence .eq. 1) call sources_derivatives_computation(n)
         if (rungekutta .eq. 10) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, impdiagt, sht)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_loop5_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
         else
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, impdiagt, sht)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_loop6_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

         end if
      end if

   end subroutine calculate_jacobian

   subroutine calculate_jacobian_loop1_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, srf, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      impdiag(i, 2, 3) = -srf_velocity(3) * ielem_totvolume(i)
      impdiag(i, 2, 4) = srf_velocity(2) * ielem_totvolume(i)
      impdiag(i, 3, 2) = srf_velocity(3) * ielem_totvolume(i)
      impdiag(i, 3, 4) = -srf_velocity(1) * ielem_totvolume(i)
      impdiag(i, 4, 2) = -srf_velocity(2) * ielem_totvolume(i)
      impdiag(i, 4, 3) = srf_velocity(1) * ielem_totvolume(i)
   end subroutine calculate_jacobian_loop1_cell

   subroutine calculate_jacobian_loop2_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, srf, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      srf = rec_mrf(i)
      if (rec_mrf(i) .eq. 1) then
         impdiag(i, 2, 3) = -rec_mrf_velocity(i, 3) * ielem_totvolume(i)
         impdiag(i, 2, 4) = rec_mrf_velocity(i, 2) * ielem_totvolume(i)
         impdiag(i, 3, 2) = rec_mrf_velocity(i, 3) * ielem_totvolume(i)
         impdiag(i, 3, 4) = -rec_mrf_velocity(i, 1) * ielem_totvolume(i)
         impdiag(i, 4, 2) = -rec_mrf_velocity(i, 2) * ielem_totvolume(i)
         impdiag(i, 4, 3) = rec_mrf_velocity(i, 1) * ielem_totvolume(i)
      end if
      srf = 0
   end subroutine calculate_jacobian_loop2_cell

   subroutine calculate_jacobian_loop3_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call add_time_and_realgas_source(n, i, impdiag(i, :, :))
   end subroutine calculate_jacobian_loop3_cell

   subroutine calculate_jacobian_loop4_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call add_time_and_realgas_source(n, i, impdiag(i, :, :))
   end subroutine calculate_jacobian_loop4_cell

   subroutine calculate_jacobian_loop5_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, srf, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & (ielem_dtl(i)))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / (ielem_dtl(i))) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / (ielem_dtl(i))), turb_diag_abs_floor))
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i))
         end do
      end if
   end subroutine calculate_jacobian_loop5_cell

   subroutine calculate_jacobian_loop6_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, srf, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations

            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + ielem_totvolume(i) * ((1.00 / &
               & ielem_dtl(i)) + (1.50 / dt))))
            impdiagt(i, nvar) = max((ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiagt(i, nvar)) - sht(i, &
               & nvar)), max(turb_diag_floor_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt))), &
               & turb_diag_abs_floor))
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiagt(i, 1))
         end do
      end if
   end subroutine calculate_jacobian_loop6_cell

   subroutine calculate_jacobian_inner_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, srf, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      identity1(:, :) = zero
      do l = 1, nof_variables
         identity1(l, l) = 1.00
      end do

      i = el_int(ii)
      iconsidered = i
      impdiag(i, :, :) = zero
      impoff(i, :, :, :) = zero
      if (turbulence .eq. 1) then
         impdiagt(i, :) = zero
         impofft(i, :, :) = zero
      end if

      if (mrf .eq. 1) then
         srf = rec_mrf(i)
      end if
      do l = 1, ielem_ifca(i) !for all their faces
         godflux2 = zero
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = (cos(angle1) * sin(angle2))
         ny = (sin(angle1) * sin(angle2))
         nz = (cos(angle2))
         mul1 = ielem_surf(i, l)
         b_code = 0
         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            !left additional equations flow state
            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
            cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, 1:turbulenceequations + passivescalar)

         end if

         call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call rotateb(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         if (realgas .eq. 0) then
            asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         else
            asound1 = sqrt((leftv(5) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt((rightv(5) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         end if
         if (rec_mrf(i) .eq. 1) then
            !retrieve the rotational velocity (at the gaussian point just for second order)
            srf_speed(2:4) = rec_rotvel(i, l, 1, 1:3)
            call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
            !calculate the new eigenvalue for rotating reference frame
            asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1) - srf_speedrot(2))
            asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1) - srf_speedrot(2))
         end if
         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2
            mul1 = ielem_surf(i, l)

            viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
               & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
               & * viscots / (ielem_dih(i, l) * prandtl)))
            if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
            vpp = max(asound1, asound2) + viscots

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)

                  call fill_turbulence_gradients(3, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if (turbulence .eq. 1) then
                  do nvar = 1, turbulenceequations

                     vpp = max(asound1, asound2) + viscots
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2 * vpp)) * mul1)
                  end do
               end if
               if (passivescalar .gt. 0) then
                  do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                     viscl(1) = viscl(1) / schmidt_lam
                     viscl(2) = viscl(2) / schmidt_lam
                     if (turbulence .eq. 1) then
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                     else
                        viscots = 0.5 * ((viscl(1) + viscl(2)))

                     end if

                     viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                     vpp = max(asound1, asound2) + viscots
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * ((vpp)) * mul1)
                  end do
               end if
            end if
         else

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

         end if

      end do
   end subroutine calculate_jacobian_inner_cell

   subroutine calculate_jacobian_bound_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, srf, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      identity1(:, :) = zero
      do l = 1, nof_variables
         identity1(l, l) = 1.00
      end do

      i = el_bnd(ii)
      iconsidered = i

      impdiag(i, :, :) = 0.0
      impoff(i, :, :, :) = 0.0
      if (turbulence .eq. 1) then
         impdiagt(i, :) = 0.0
         impofft(i, :, :) = 0.0
      end if
      if (mrf .eq. 1) then
         srf = rec_mrf(i)
      end if
      do l = 1, ielem_ifca(i)
         b_code = 0
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = (cos(angle1) * sin(angle2))
         ny = (sin(angle1) * sin(angle2))
         nz = (cos(angle2))
         mul1 = ielem_surf(i, l)
         b_code = 0
         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         if (rec_mrf(i) .eq. 1) then
            !retrieve rotational velocity in case of rotating reference frame to calculate the correct value of the boundary condition
            srf_speed(2:4) = rec_rotvel(i, l, 1, 1:3)
            call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
         end if
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)

         end if

         if (ielem_ineighb(i, l) .eq. n) then !my cpu only
            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               !periodic in my cpu
               if ((ibound_icode(ielem_ibounds(i, l)) .eq. 5) .or. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)
                  if ((per_rot .eq. 1) .and. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then
                     cright(2:4) = rotate_per_1(cright(2:4), ibound_icode(ielem_ibounds(i, l)), angle_per)
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)

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

                  leftv(1:nof_variables) = cleft(1:nof_variables)
                  b_code = ibound_icode(ielem_ibounds(i, l))

                  call boundarys(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, cturbl, &
                     & cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                  cright(1:nof_variables) = rightv(1:nof_variables)

               end if
            else
               cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                     & 1:turbulenceequations + passivescalar)

               end if

            end if
         else !in other cpus they can only be periodic or mpi neighbours

            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               !periodic in other cpu
               if ((ibound_icode(ielem_ibounds(i, l)) .eq. 5) .or. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then

                  ! cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                  ! (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                  if ((per_rot .eq. 1) .and. (ibound_icode(ielem_ibounds(i, l)) .eq. 50)) then
                     cright(2:4) = rotate_per_1(cright(2:4), ibound_icode(ielem_ibounds(i, l)), angle_per)
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     !                                                      cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                     !                                                      (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else

               ! cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
               ! (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

               nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
               lf = rec_ihexl(i, 1, ielem_indexi(i, L))
               rowf = halo_offset(nf) + lf - 1
               cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  !cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                  !(rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)

               end if

!
            end if
         end if

         call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef(n, cleft_rot, cleft, angle1, angle2)

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call rotateb(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
         if (realgas .eq. 0) then
            asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         else
            asound1 = sqrt((leftv(5) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt((rightv(5) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         end if
         if (rec_mrf(i) .eq. 1) then
            asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1) - srf_speedrot(2))
            asound2 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cright_rot(2) / cright_rot(1) - srf_speedrot(2))
         end if

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2
            mul1 = ielem_surf(i, l)

            viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
               & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
               & * viscots / (ielem_dih(i, l) * prandtl)))
            if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
            vpp = max(asound1, asound2) + viscots

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(3, i, eddyfl)


                  eddyfr = eddyfl

                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if (turbulence .eq. 1) then
                  do nvar = 1, turbulenceequations

                     vpp = max(asound1, asound2) + viscots
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2 * vpp)) * mul1)
                  end do
               end if
               if (passivescalar .gt. 0) then
                  do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                     viscl(1) = viscl(1) / schmidt_lam
                     viscl(2) = viscl(2) / schmidt_lam

                     if (turbulence .eq. 1) then
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                     else
                        viscots = 0.5 * ((viscl(1) + viscl(2)))

                     end if

                     viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                     vpp = max(asound1, asound2) + viscots
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * ((vpp)) * mul1)
                  end do
               end if
            end if
         else

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

         end if

         if (b_code == 9 .and. duct_outlet_face_flux_enabled()) then
            if (turbulenceequations + passivescalar > 0) then
               call add_duct_outlet_jacobian(n, i, cleft(1:nof_variables), [nx, ny, nz], max(asound1, asound2), mul1, &
                  & impdiag(i, :, :), impoff(i, l, :, :), impdiagt(i, :), impofft(i, l, :))
            else
               call add_duct_outlet_jacobian(n, i, cleft(1:nof_variables), [nx, ny, nz], max(asound1, asound2), mul1, &
                  & impdiag(i, :, :), impoff(i, l, :, :))
            end if
         end if

      end do
   end subroutine calculate_jacobian_bound_cell

   subroutine calculate_jacobian_2d(n)
      implicit none
      !> @brief
      !> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      integer :: i, kmaxe, ii

      kmaxe = xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_interior) &
!$omp& map(alloc: el_int, gamma_in, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_ineigh) &
!$omp& map(alloc: ielem_surf, ielem_walldist, impdiag, impdiagt, impoff, impofft, mp_janaf, mp_m, mp_pinf) &
!$omp& map(alloc: mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_grads, rg_hzero, rg_molm, rg_thetag, rgs_ab, rgs_bb) &
!$omp& map(alloc: rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_interior !for all the interior elements
         call calculate_jacobian_2d_inner_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_bounded) &
!$omp& map(alloc: bleed_end, bleed_plenum, bleed_porosity, bleed_start, catalytic_con, el_bnd, gamma_in, halo_offset) &
!$omp& map(alloc: ibound_icode, ielem_bleedn, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ibounds) &
!$omp& map(alloc: ielem_ifca, ielem_indexi, ielem_ineigh, ielem_ineighb, ielem_nodes_faces, ielem_surf) &
!$omp& map(alloc: ielem_walldist, impdiag, impdiagt, impoff, impofft, inoder4_cord, mp_a_in, mp_janaf, mp_m, mp_pinf) &
!$omp& map(alloc: mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_grads, rec_ihexl, rec_ihexn, rec_local, rg_hzero) &
!$omp& map(alloc: rg_molm, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, solhir) &
!$omp& map(alloc: u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_bounded
         call calculate_jacobian_2d_bound_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

      if (realgas .eq. 0) then
         if (rungekutta .eq. 10) then

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, impdiag)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_2d_loop1_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
         else
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, impdiag)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_2d_loop2_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
         end if

      end if

      if (realgas .eq. 1) then
         if (rungekutta .eq. 10) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: gamma_in, ielem_dtl, ielem_totvolume, impdiag, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rg_hzero, rg_molm, rg_thetag, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_2d_loop3_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
         else
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: gamma_in, ielem_dtl, ielem_totvolume, impdiag, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rg_hzero, rg_molm, rg_thetag, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_2d_loop4_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
         end if

      end if

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         if (turbulence .eq. 1) call sources_derivatives_computation2d(n)
         if (rungekutta .eq. 10) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, impdiagt, sht)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_2d_loop5_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
         else
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(i) &
!$omp& firstprivate(kmaxe) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, impdiagt, sht)
#else
!$omp do &
!$omp& private(i)
#endif
            do i = 1, kmaxe
               call calculate_jacobian_2d_loop6_cell(n, i)
            end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

         end if
      end if

   end subroutine calculate_jacobian_2d

   subroutine calculate_jacobian_2d_loop1_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, k, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, kas
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      do j = 1, nof_variables

         impdiag(i, j, j) = impdiag(i, j, j) + (ielem_totvolume(i) / ielem_dtl(i))

      end do
   end subroutine calculate_jacobian_2d_loop1_cell

   subroutine calculate_jacobian_2d_loop2_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, k, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, kas
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      do j = 1, nof_variables

         impdiag(i, j, j) = ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiag(i, j, j))

      end do
   end subroutine calculate_jacobian_2d_loop2_cell

   subroutine calculate_jacobian_2d_loop3_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call add_time_and_realgas_source(n, i, impdiag(i, :, :))
   end subroutine calculate_jacobian_2d_loop3_cell

   subroutine calculate_jacobian_2d_loop4_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      call add_time_and_realgas_source(n, i, impdiag(i, :, :))
   end subroutine calculate_jacobian_2d_loop4_cell

   subroutine calculate_jacobian_2d_loop5_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, k, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, kas
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
!
            impdiagt(i, nvar) = impdiagt(i, nvar) + ((ielem_totvolume(i) / (ielem_dtl(i)))) - sht(i, nvar)
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem_totvolume(i) / (ielem_dtl(i)))
         end do
      end if
   end subroutine calculate_jacobian_2d_loop5_cell

   subroutine calculate_jacobian_2d_loop6_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, k, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, kas
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
            impdiagt(i, nvar) = ielem_totvolume(i) * ((1.00 / (ielem_dtl(i))) + (1.50 / dt)) + (impdiagt(i, 1)) - sht(i, nvar)
!     impdiagt(i,nvar)=(ielem_totvolume(i)*((1.00/ielem_dtl(i))+((1.50/dt)*impdiagt(i,nvar))))-sht(i,nvar)
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = ielem_totvolume(i) * ((1.00 / (ielem_dtl(i))) + (1.50 / dt)) + (impdiagt(i, 1))
!     impdiagt(i,nvar)=(ielem_totvolume(i)*((1.00/ielem_dtl(i))+((1.50/dt)*impdiagt(i,nvar))))
         end do
      end if
   end subroutine calculate_jacobian_2d_loop6_cell

   subroutine calculate_jacobian_2d_inner_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc, k, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, kas
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      identity1(:, :) = zero
      do l = 1, nof_variables
         identity1(l, l) = 1.00
      end do

      i = el_int(ii)
      iconsidered = i
      impdiag(i, :, :) = zero
      impoff(i, :, :, :) = zero

      if (turbulence .eq. 1) then
         impdiagt(i, :) = zero
         impofft(i, :, :) = zero
      end if

      do l = 1, ielem_ifca(i) !for all their faces
         b_code = 0
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = angle1
         ny = angle2
         mul1 = ielem_surf(i, l)

         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            !left additional equations flow state
            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
            cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, 1:turbulenceequations + passivescalar)

         end if

         call rotatef2d(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb2d(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         if (realgas .eq. 0) then
            asound1 = sqrt(leftv(4) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(4) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         else
            asound1 = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         end if

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2
            mul1 = ielem_surf(i, l)

            if (realgas .eq. 0) then
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots
            else
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gammal / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots

            end if

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(2, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if (turbulence .eq. 1) then
                  do nvar = 1, turbulenceequations

                     vpp = max(asound1, asound2) + viscots

                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2 * vpp)) * mul1)
                  end do
               end if
               if (passivescalar .gt. 0) then

                  do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                     viscl(1) = viscl(1) / schmidt_lam
                     viscl(2) = viscl(2) / schmidt_lam

                     if (turbulence .eq. 1) then
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                     else
                        viscots = 0.5 * ((viscl(1) + viscl(2)))

                     end if

                     viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                     vpp = max(asound1, asound2) + viscots
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * ((vpp)) * mul1)
                  end do
               end if
            end if
         else

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

         end if
      end do
   end subroutine calculate_jacobian_2d_inner_cell

   subroutine calculate_jacobian_2d_bound_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc, k, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, kas
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      identity1(:, :) = zero
      do l = 1, nof_variables
         identity1(l, l) = 1.00
      end do

      i = el_bnd(ii)
      iconsidered = i

      impdiag(i, :, :) = zero
      impoff(i, :, :, :) = zero
      if (turbulence .eq. 1) then
         impdiagt(i, :) = zero
         impofft(i, :, :) = zero
      end if

      do l = 1, ielem_ifca(i)
         mul1 = ielem_surf(i, l)
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = angle1
         ny = angle2

         b_code = 0
         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)

         end if

         if (ielem_ineighb(i, l) .eq. n) then !my cpu only
            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)

                  end if

                  kas = 1

               else
                  !not periodic ones in my cpu

                  facex = l; iconsidered = i
                  call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                  n_node = 2
                  cords(1:2) = zero
                  call cordinates2(n, nodes_list, n_node, cords(1:2))

                  poy(1) = cords(2)
                  pox(1) = cords(1)

                  leftv(1:nof_variables) = cleft(1:nof_variables)
                  b_code = ibound_icode(ielem_ibounds(i, l))

                  if (realgas .eq. 0) then
                     call boundarys2d(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & cturbl, cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                     cright(1:nof_variables) = rightv(1:nof_variables)
                  else
                     cright(1:nof_variables) = cleft(1:nof_variables)
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                        cturbr(1:turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
                     end if
                  end if

                  kas = 2

               end if
            else
               cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                     & 1:turbulenceequations + passivescalar)

               end if

               kas = 3

            end if
         else !in other cpus they can only be periodic or mpi neighbours

            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

                  !cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                  !(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                  kas = 4

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     ! cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                     !(rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else
               kas = 5

               !  cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
               !  (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

               nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
               lf = rec_ihexl(i, 1, ielem_indexi(i, L))
               rowf = halo_offset(nf) + lf - 1
               cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  !cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                  !(rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)

               end if

!
            end if
         end if

         call rotatef2d(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb2d(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         if (realgas .eq. 0) then
            asound1 = sqrt(leftv(4) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(4) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         else
            asound1 = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         end if

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2
            mul1 = ielem_surf(i, l)
            if (realgas .eq. 0) then
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots
            else
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gammal / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots

            end if

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(2, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if (turbulence .eq. 1) then
                  do nvar = 1, turbulenceequations

                     vpp = max(asound1, asound2) + viscots
!
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2 * vpp)) * mul1)
                  end do
               end if
               if (passivescalar .gt. 0) then
                  do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                     viscl(1) = viscl(1) / schmidt_lam
                     viscl(2) = viscl(2) / schmidt_lam

                     if (turbulence .eq. 1) then
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                     else
                        viscots = 0.5 * ((viscl(1) + viscl(2)))

                     end if

                     viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                     vpp = max(asound1, asound2) + viscots
                     impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * ((vpp)) * mul1)
                     impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * ((vpp)) * mul1)
                  end do
               end if
            end if
         else

            impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
               & identity1(1:nof_variables, 1:nof_variables))) * mul1)
            call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
            convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
            impoff(i, l, 1:nof_variables, 1:nof_variables) = impoff(i, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
               & convj(1:nof_variables, 1:nof_variables)) &
               & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

         end if

         if (b_code == 9 .and. duct_outlet_face_flux_enabled()) then
            if (turbulenceequations + passivescalar > 0) then
               call add_duct_outlet_jacobian(n, i, cleft(1:nof_variables), [nx, ny, 0.0e0], max(asound1, asound2), mul1, &
                  & impdiag(i, :, :), impoff(i, l, :, :), impdiagt(i, :), impofft(i, l, :))
            else
               call add_duct_outlet_jacobian(n, i, cleft(1:nof_variables), [nx, ny, 0.0e0], max(asound1, asound2), mul1, &
                  & impdiag(i, :, :), impoff(i, l, :, :))
            end if
         end if

      end do
   end subroutine calculate_jacobian_2d_bound_cell

   subroutine calculate_jacobianlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 3d with low-memory footprint
      integer, intent(in) :: n, iconsidered
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, intent(inout) :: impdiagt(:, :)
      real, intent(inout) :: impdiag(:, :, :), impofft(:, :, :)
      real, intent(inout) :: impoff(:, :, :, :)
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl
      real, dimension(gpu_max_turbulence) :: source_t

      identity1(:, :) = zero

      identity1(:, :) = zero
      do j = 1, nof_variables
         identity1(j, j) = 1.00
      end do

      i = iconsidered
      if (ielem_interior(iconsidered) .eq. 0) then

         i = iconsidered
         impdiag(1, :, :) = zero
         impoff(1, :, :, :) = zero
         if (turbulence .eq. 1) then
            impdiagt(1, :) = zero
            impofft(1, :, :) = zero
         end if

         do l = 1, ielem_ifca(i) !for all their faces
            godflux2 = zero
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = (cos(angle1) * sin(angle2))
            ny = (sin(angle1) * sin(angle2))
            nz = (cos(angle2))
            mul1 = ielem_surf(i, l)

            cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               !left additional equations flow state
               cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
               cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                  & 1:turbulenceequations + passivescalar)

            end if

            call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then !application of the low mach number correction
               leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
               call lmacht(n, leftv, rightv)
               cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               call rotatef(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotatef(n, cleft, cleft_rot, angle1, angle2)

            end if

            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

            asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))

            vpp = max(asound1, asound2)

            if (itestcase .eq. 4) then
               leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
               viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
               call get_visc_conduct(n, leftv, rightv, viscl, laml)

               viscots = (viscl(1) + viscl(2)) * oo2
               mul1 = ielem_surf(i, l)

               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots

               if (turbulence .eq. 1) then
                  if (turbulencemodel .eq. 1) then
                     turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                     call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if
                  if (turbulencemodel .eq. 2) then

                     eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)

                     call fill_turbulence_gradients(3, i, eddyfl)


                     eddyfr = eddyfl
                     call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if

                  viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                     & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                     & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

                  if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

                  vpp = max(asound1, asound2) + viscots
               end if

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (turbulence .eq. 1) then
                     do nvar = 1, turbulenceequations

                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2 * vpp)) * mul1)
                     end do
                  end if
                  if (passivescalar .gt. 0) then
                     do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                        viscl(1) = viscl(1) / schmidt_lam
                        viscl(2) = viscl(2) / schmidt_lam
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                        viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2 * ((vpp)) * mul1)
                     end do
                  end if
               end if
            else

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            end if
         end do

      else
         i = iconsidered
         impdiag(1, :, :) = zero
         impoff(1, :, :, :) = zero
         if (turbulence .eq. 1) then
            impdiagt(1, :) = zero
            impofft(1, :, :) = zero
         end if

         do l = 1, ielem_ifca(i)
            mul1 = ielem_surf(i, l)
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = (cos(angle1) * sin(angle2))
            ny = (sin(angle1) * sin(angle2))
            nz = (cos(angle2))

            b_code = 0
            cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)

            end if

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                           & 1:turbulenceequations + passivescalar)

                     end if

                  else
                     !not periodic ones in my cpu

                     facex = l;
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

                     leftv(1:nof_variables) = cleft(1:nof_variables)
                     b_code = ibound_icode(ielem_ibounds(i, l))

                     call boundarys(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & cturbl, cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                     cright(1:nof_variables) = rightv(1:nof_variables)

                  end if
               else
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)

                  end if

               end if
            else !in other cpus they can only be periodic or mpi neighbours

               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

                     !cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                     !(rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        ! cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                        !(rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                        nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                        lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                        rowf = halo_offset(nf) + lf - 1
                        cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else

                  ! cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                  ! (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     ! cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
                     !(rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

!
               end if
            end if

            call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then !application of the low mach number correction
               leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
               call lmacht(n, leftv, rightv)
               cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               call rotatef(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotatef(n, cleft, cleft_rot, angle1, angle2)

            end if

            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

            asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))

            vpp = max(asound1, asound2)

            if (itestcase .eq. 4) then
               leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
               viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
               call get_visc_conduct(n, leftv, rightv, viscl, laml)

               viscots = (viscl(1) + viscl(2)) * oo2
               mul1 = ielem_surf(i, l)

               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots

               if (turbulence .eq. 1) then
                  if (turbulencemodel .eq. 1) then
                     turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                     call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if
                  if (turbulencemodel .eq. 2) then
                     eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)

                     call fill_turbulence_gradients(3, i, eddyfl)


                     eddyfr = eddyfl
                     call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if

                  viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                     & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                     & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

                  if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

                  vpp = max(asound1, asound2) + viscots
               end if

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (turbulence .eq. 1) then
                     do nvar = 1, turbulenceequations

                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2 * vpp)) * mul1)
                     end do
                  end if
                  if (passivescalar .gt. 0) then
                     do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                        viscl(1) = viscl(1) / schmidt_lam
                        viscl(2) = viscl(2) / schmidt_lam
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                        viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2 * ((vpp)) * mul1)
                     end do
                  end if
               end if
            else

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            end if

            if (b_code == 9 .and. duct_outlet_face_flux_enabled()) then
               call add_duct_outlet_jacobian(n, i, cleft(1:nof_variables), [nx, ny, nz], max(asound1, asound2), mul1, &
                  & impdiag(1, :, :), impoff(1, l, :, :), impdiagt(1, :), impofft(1, l, :))
            end if

         end do
      end if

      call add_time_and_realgas_source(n, i, impdiag(1, :, :))

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         call sources_derivatives(n, iconsidered, source_t)
         sht(i, 1:turbulenceequations) = (source_t(1:turbulenceequations) * ielem_totvolume(i))
         if (rungekutta .eq. 10) then

            if (turbulence .eq. 1) then
               do nvar = 1, turbulenceequations
                  sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(1, nvar) + (ielem_totvolume(i) / &
                     & ielem_dtl(i))))
                  impdiagt(1, nvar) = max((impdiagt(1, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
                     & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
               end do
            end if
            if (passivescalar .gt. 0) then
               do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                  impdiagt(1, nvar) = impdiagt(1, nvar) + (ielem_totvolume(i) / ielem_dtl(i))
               end do
            end if

         else

            if (turbulence .eq. 1) then
               do nvar = 1, turbulenceequations
                  sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) &
                     & + ((1.50 / dt) * impdiagt(1, nvar)))))
                  impdiagt(1, nvar) = max(((ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + ((1.50 / dt) * impdiagt(1, &
                     & nvar)))) - sht(i, nvar)), max(turb_diag_floor_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + ((1.50 &
                     & / dt) * impdiagt(1, nvar)))), turb_diag_abs_floor))
               end do
            end if
            if (passivescalar .gt. 0) then
               do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                  impdiagt(1, nvar) = (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + ((1.50 / dt) * impdiagt(1, nvar))))
               end do
            end if

         end if
      end if

   end subroutine calculate_jacobianlm

   subroutine calculate_jacobian_2dlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d with low-memory footprint
      integer, intent(in) :: n, iconsidered
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, intent(inout) :: impdiagt(:, :)
      real, intent(inout) :: impdiag(:, :, :), impofft(:, :, :)
      real, intent(inout) :: impoff(:, :, :, :)
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl
      real, dimension(gpu_max_turbulence) :: source_t

      identity1(:, :) = zero
      do j = 1, nof_variables
         identity1(j, j) = 1.00
      end do

      i = iconsidered
      if (ielem_interior(iconsidered) .eq. 0) then

         i = iconsidered
         impdiag(1, :, :) = zero
         impoff(1, :, :, :) = zero
         if (turbulence .eq. 1) then
            impdiagt(1, :) = zero
            impofft(1, :, :) = zero
         end if

         do l = 1, ielem_ifca(i) !for all their faces
            godflux2 = zero
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = angle1
            ny = angle2
            mul1 = ielem_surf(i, l)

            cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               !left additional equations flow state
               cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
               cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                  & 1:turbulenceequations + passivescalar)

            end if

            call rotatef2d(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then !application of the low mach number correction
               leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
               call lmacht2d(n, leftv, rightv)
               cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotateb2d(n, cright, cright_rot, angle1, angle2)
               call rotateb2d(n, cleft, cleft_rot, angle1, angle2) !rotate wrt to normalvector of face and so
            end if

            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

            if (realgas .eq. 0) then
               asound1 = sqrt(leftv(4) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
               asound2 = sqrt(rightv(4) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
            else
               asound1 = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
               asound2 = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
            end if

            vpp = max(asound1, asound2)

            if (itestcase .eq. 4) then
               leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
               viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
               call get_visc_conduct(n, leftv, rightv, viscl, laml)

               viscots = (viscl(1) + viscl(2)) * oo2
               mul1 = ielem_surf(i, l)

               if (realgas .eq. 0) then
                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                     & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                     & * viscots / (ielem_dih(i, l) * prandtl)))
               else
                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                     & / ielem_dih(i, l), ((gammal / (0.5 * (cleft(1) + cright(1)))) &
                     & * viscots / (ielem_dih(i, l) * prandtl)))
               end if
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots

               if (turbulence .eq. 1) then
                  if (turbulencemodel .eq. 1) then
                     turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                     call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if
                  if (turbulencemodel .eq. 2) then
                     eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                     call fill_turbulence_gradients(2, i, eddyfl)


                     eddyfr = eddyfl
                     call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if

                  viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                     & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                     & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

                  if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

                  vpp = max(asound1, asound2) + viscots
               end if

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (turbulence .eq. 1) then
                     do nvar = 1, turbulenceequations

                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2 * vpp)) * mul1)
                     end do
                  end if
                  if (passivescalar .gt. 0) then
                     do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                        viscl(1) = viscl(1) / schmidt_lam
                        viscl(2) = viscl(2) / schmidt_lam
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                        viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2 * ((vpp)) * mul1)
                     end do
                  end if
               end if
            else

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            end if
         end do
      else

         impdiag(1, :, :) = 0.0
         impoff(1, :, :, :) = 0.0
         if (turbulence .eq. 1) then
            impdiagt(1, :) = 0.0
            impofft(1, :, :) = 0.0
         end if

         do l = 1, ielem_ifca(i)
            mul1 = ielem_surf(i, l)
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = angle1
            ny = angle2

            b_code = 0
            cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

               cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)

            end if

            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                           & 1:turbulenceequations + passivescalar)

                     end if

                  else
                     !not periodic ones in my cpu

                     facex = l;
                     call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                     n_node = 2
                     cords(1:2) = zero
                     call cordinates2(n, nodes_list, n_node, cords(1:2))

                     poy(1) = cords(2)
                     pox(1) = cords(1)

                     leftv(1:nof_variables) = cleft(1:nof_variables)
                     b_code = ibound_icode(ielem_ibounds(i, l))

                     if (realgas .eq. 0) then
                        call boundarys2d(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & cturbl, cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                        cright(1:nof_variables) = rightv(1:nof_variables)
                     else
                        cright(1:nof_variables) = cleft(1:nof_variables)
                        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                           cturbr(1:turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
                        end if
                     end if

                  end if
               else
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)

                  end if

               end if
            else !in other cpus they can only be periodic or mpi neighbours

               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                               cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

!                                                               cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                        nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                        lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                        rowf = halo_offset(nf) + lf - 1
                        cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                           & turbulenceequations + passivescalar)

                     end if

                  end if
               else

!                                                               cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

!                                                               cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

!
               end if
            end if

            call rotatef2d(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then !application of the low mach number correction
               leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
               call lmacht2d(n, leftv, rightv)
               cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotateb2d(n, cright, cright_rot, angle1, angle2)
               call rotateb2d(n, cleft, cleft_rot, angle1, angle2) !rotate wrt to normalvector of face and so
            end if

            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

            if (realgas .eq. 0) then
               asound1 = sqrt(leftv(4) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
               asound2 = sqrt(rightv(4) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
            else
               asound1 = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
               asound2 = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
            end if

            vpp = max(asound1, asound2)

            if (itestcase .eq. 4) then
               leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
               viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
               call get_visc_conduct(n, leftv, rightv, viscl, laml)

               viscots = (viscl(1) + viscl(2)) * oo2
               mul1 = ielem_surf(i, l)

               if (realgas .eq. 0) then
                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                     & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                     & * viscots / (ielem_dih(i, l) * prandtl)))
               else
                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                     & / ielem_dih(i, l), ((gammal / (0.5 * (cleft(1) + cright(1)))) &
                     & * viscots / (ielem_dih(i, l) * prandtl)))
               end if
               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
               vpp = max(asound1, asound2) + viscots

               if (turbulence .eq. 1) then
                  if (turbulencemodel .eq. 1) then
                     turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                     call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if
                  if (turbulencemodel .eq. 2) then
                     eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                     call fill_turbulence_gradients(2, i, eddyfl)


                     eddyfr = eddyfl
                     call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                  end if

                  viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                  viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                     & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                     & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

                  if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

                  vpp = max(asound1, asound2) + viscots
               end if

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (turbulence .eq. 1) then
                     do nvar = 1, turbulenceequations

                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2 * vpp)) * mul1)
                     end do
                  end if
                  if (passivescalar .gt. 0) then
                     do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                        viscl(1) = viscl(1) / schmidt_lam
                        viscl(2) = viscl(2) / schmidt_lam
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

                        viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                        vpp = max(asound1, asound2) + viscots
                        impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2 * ((vpp)) * mul1)
                        impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2 * ((vpp)) * mul1)
                     end do
                  end if
               end if
            else

               impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2 * ((vpp * &
                  & identity1(1:nof_variables, 1:nof_variables))) * mul1)
               call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
               convj(1:nof_variables, 1:nof_variables) = eigvl(1:nof_variables, 1:nof_variables)
               impoff(1, l, 1:nof_variables, 1:nof_variables) = impoff(1, l, 1:nof_variables, 1:nof_variables) + (((oo2 * &
                  & convj(1:nof_variables, 1:nof_variables)) &
                  & - ((oo2 * vpp) * identity1(1:nof_variables, 1:nof_variables))) * mul1)

            end if

            if (b_code == 9 .and. duct_outlet_face_flux_enabled()) then
               call add_duct_outlet_jacobian(n, i, cleft(1:nof_variables), [nx, ny, 0.0e0], max(asound1, asound2), mul1, &
                  & impdiag(1, :, :), impoff(1, l, :, :), impdiagt(1, :), impofft(1, l, :))
            end if

         end do
      end if

      call add_time_and_realgas_source(n, i, impdiag(1, :, :))

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
         call sources_derivatives2d(n, iconsidered, source_t)
         sht(i, 1:turbulenceequations) = (source_t(1:turbulenceequations) * ielem_totvolume(i))
         if (rungekutta .eq. 10) then

            if (turbulence .eq. 1) then
               do nvar = 1, turbulenceequations
                  sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(1, nvar) + (ielem_totvolume(i) / &
                     & ielem_dtl(i))))
                  impdiagt(1, nvar) = max((impdiagt(1, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
                     & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
               end do
            end if
            if (passivescalar .gt. 0) then
               do nvar = turbulenceequations + 1, turbulenceequations + passivescalar

                  impdiagt(1, nvar) = impdiagt(1, nvar) + (ielem_totvolume(i) / ielem_dtl(i))
               end do
            end if

         else

            if (turbulence .eq. 1) then
               do nvar = 1, turbulenceequations
                  sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) &
                     & + (1.50 / dt)) + (impdiagt(1, 1))))
                  impdiagt(1, nvar) = max((ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiagt(1, 1)) - sht(i, &
                     & nvar)), max(turb_diag_floor_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt))), &
                     & turb_diag_abs_floor))

               end do
            end if
            if (passivescalar .gt. 0) then
               do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                  impdiagt(1, nvar) = ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiagt(1, 1))
!     impdiagt(1,nvar)=(ielem_totvolume(i)*((1.00/ielem_dtl(i))+((1.50/dt)*impdiagt(1,nvar))))
               end do
            end if

         end if
      end if

   end subroutine calculate_jacobian_2dlm

   subroutine calculate_jacobian_2d_mf(n)
      implicit none
      !> @brief
      !> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      integer :: i, kmaxe, ii
      kmaxe = xmpielrank(n)

!$omp do
      do ii = 1, nof_interior !for all the interior elements
         call calculate_jacobian_2d_mf_inner_cell(n, ii)
      end do
!$omp end do

!$omp do
      do ii = 1, nof_bounded
         call calculate_jacobian_2d_mf_bound_cell(n, ii)
      end do
!$omp end do

      if (rungekutta .eq. 10) then
!$omp do
         do i = 1, kmaxe
            call calculate_jacobian_2d_mf_loop1_cell(n, i)
         end do
!$omp end do
      else
!$omp do
         do i = 1, kmaxe
            call calculate_jacobian_2d_mf_loop2_cell(n, i)
         end do
!$omp end do
      end if

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         if (turbulence .eq. 1) call sources_derivatives_computation2d(n)
         if (rungekutta .eq. 10) then
!$omp do
            do i = 1, kmaxe
               call calculate_jacobian_2d_mf_loop3_cell(n, i)
            end do
!$omp end do
         else
!$omp do
            do i = 1, kmaxe
               call calculate_jacobian_2d_mf_loop4_cell(n, i)
            end do
!$omp end do

         end if
      end if

   end subroutine calculate_jacobian_2d_mf

   subroutine calculate_jacobian_2d_mf_loop1_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz

      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      impdiag_mf(i) = (impdiag_mf(i)) + (ielem_totvolume(i) / ielem_dtl(i))

      if ((realgas .eq. 1) .and. (rg_relax .eq. 2)) then
         ! Pack the strided state before evaluating the scalar source bound.
         leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         impdiag_mf(i) = impdiag_mf(i) + ielem_totvolume(i) * &
            & realgas_source_radius(n, leftv)
      end if
   end subroutine calculate_jacobian_2d_mf_loop1_cell

   subroutine calculate_jacobian_2d_mf_loop2_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz

      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      impdiag_mf(i) = ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiag_mf(i))
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2)) then
         ! Pack the strided state before evaluating the scalar source bound.
         leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         impdiag_mf(i) = impdiag_mf(i) + ielem_totvolume(i) * &
            & realgas_source_radius(n, leftv)
      end if
   end subroutine calculate_jacobian_2d_mf_loop2_cell

   subroutine calculate_jacobian_2d_mf_loop3_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz

      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & ielem_dtl(i))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & ielem_dtl(i))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
         end do
      end if
   end subroutine calculate_jacobian_2d_mf_loop3_cell

   subroutine calculate_jacobian_2d_mf_loop4_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz

      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
!
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + ielem_totvolume(i) * ((1.00 / &
               & ielem_dtl(i)) + (1.50 / dt))))
            impdiagt(i, nvar) = max((ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiagt(i, nvar)) - sht(i, &
               & nvar)), max(turb_diag_floor_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt))), &
               & turb_diag_abs_floor))
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & ielem_dtl(i))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
         end do
      end if

   end subroutine calculate_jacobian_2d_mf_loop4_cell

   subroutine calculate_jacobian_2d_mf_inner_cell(n, ii)
      implicit none
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc, kas, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz

      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      i = el_int(ii)
      iconsidered = i
      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
         impdiagt(i, 1:turbulenceequations + passivescalar) = 0.0
         impofft(i, :, :) = zero
      end if
      do l = 1, ielem_ifca(i) !for all their faces
         b_code = 0
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = angle1
         ny = angle2
         mul1 = ielem_surf(i, l)

         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

!                                      cleft(1:nof_variables)=rec_uleft(i,1:nof_variables,l,1)
!                                    cright(1:nof_variables)=rec_uleft(ielem_ineigh(i,l),1:nof_variables,ielem_ineighn(i,l),1)

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            !left additional equations flow state
            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
            cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, 1:turbulenceequations + passivescalar)

         end if

         call rotatef2d(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb2d(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         if (realgas .eq. 0) then
            asound1 = sqrt(leftv(4) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(4) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         else
            asound1 = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         end if

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
            if (realgas .eq. 0) then
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
            else
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gammal / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
            end if
            if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
            vpp = max(asound1, asound2) + viscots

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(2, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

               !                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if
         end if

         impdiag_mf(i) = impdiag_mf(i) + (oo2 * vpp * mul1)
         impoff_mf(i, l) = vpp

         if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations + passivescalar
               impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * vpp * mul1)
               impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * vpp * mul1)
            end do
         end if

      end do
   end subroutine calculate_jacobian_2d_mf_inner_cell

   subroutine calculate_jacobian_2d_mf_bound_cell(n, ii)
      implicit none
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc, kas, j
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux
      integer :: b_code, nf, lf, rowf
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz

      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      i = el_bnd(ii)
      iconsidered = i

      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
         impdiagt(i, 1:turbulenceequations + passivescalar) = zero
         impofft(i, :, :) = zero
      end if
      do l = 1, ielem_ifca(i)
         mul1 = ielem_surf(i, l)
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = angle1
         ny = angle2

         b_code = 0
         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
!                                       cleft(1:nof_variables)=rec_uleft(i,1:nof_variables,l,1)
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)

         end if

         if (ielem_ineighb(i, l) .eq. n) then !my cpu only
            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)
!                                                                   cright(1:nof_variables)=rec_uleft(ielem_ineigh(i,l),1:nof_variables,ielem_ineighn(i,l),1)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)

                  end if

                  kas = 1

               else
                  !not periodic ones in my cpu

                  facex = l; iconsidered = i
                  call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                  n_node = 2
                  cords(1:2) = zero
                  call cordinates2(n, nodes_list, n_node, cords(1:2))

                  poy(1) = cords(2)
                  pox(1) = cords(1)

                  leftv(1:nof_variables) = cleft(1:nof_variables)
                  b_code = ibound_icode(ielem_ibounds(i, l))

                  if (realgas .eq. 0) then
                     call boundarys2d(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                        & cturbl, cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                     cright(1:nof_variables) = rightv(1:nof_variables)
                  else
                     cright(1:nof_variables) = cleft(1:nof_variables)
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                        cturbr(1:turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
                     end if
                  end if

                  kas = 2

               end if
            else
               cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)
!                                                               cright(1:nof_variables)=rec_uleft(ielem_ineigh(i,l),1:nof_variables,ielem_ineighn(i,l),1)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                     & 1:turbulenceequations + passivescalar)

               end if

               kas = 3

            end if
         else !in other cpus they can only be periodic or mpi neighbours

            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                               cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                  kas = 4

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

!                                                               cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else
               kas = 5

!                                                               cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

               nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
               lf = rec_ihexl(i, 1, ielem_indexi(i, L))
               rowf = halo_offset(nf) + lf - 1
               cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

!
!                                                               cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)

               end if

!
            end if
         end if

         call rotatef2d(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb2d(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         if (realgas .eq. 0) then
            asound1 = sqrt(leftv(4) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt(rightv(4) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         else
            asound1 = sqrt((leftv(4) + mp_pinfl) * gammal / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
            asound2 = sqrt((rightv(4) + mp_pinfr) * gammar / rightv(1)) + abs(cright_rot(2) / cright_rot(1))
         end if

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
            if (realgas .eq. 0) then
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
            else
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
                  & / ielem_dih(i, l), ((gammal / (0.5 * (cleft(1) + cright(1)))) &
                  & * viscots / (ielem_dih(i, l) * prandtl)))
            end if
            if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
            vpp = max(asound1, asound2) + viscots

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(2, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if (turbulence .eq. 1) then
                  do nvar = 1, turbulenceequations

                     vpp = max(asound1, asound2) + viscots
!

                  end do
               end if
               if (passivescalar .gt. 0) then
                  do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                     viscl(1) = viscl(1) / schmidt_lam
                     viscl(2) = viscl(2) / schmidt_lam

                     if (turbulence .eq. 1) then
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                     else
                        viscots = 0.5 * ((viscl(1) + viscl(2)))

                     end if

                     viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                     vpp = max(asound1, asound2) + viscots

                  end do
               end if
            end if
         else

         end if

         impdiag_mf(i) = impdiag_mf(i) + (oo2 * vpp * mul1)
         impoff_mf(i, l) = vpp
         if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations + passivescalar
               impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * vpp * mul1)
               impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * vpp * mul1)
            end do
         end if

      end do
   end subroutine calculate_jacobian_2d_mf_bound_cell

   subroutine calculate_jacobian_3d_mf(n)
      implicit none
      !> @brief
      !> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      integer :: i, kmaxe, ii
      kmaxe = xmpielrank(n)

!$omp do
      do ii = 1, nof_interior !for all the interior elements
         call calculate_jacobian_3d_mf_inner_cell(n, ii)
      end do
!$omp end do

!$omp do
      do ii = 1, nof_bounded
         call calculate_jacobian_3d_mf_bound_cell(n, ii)
      end do
!$omp end do

      if (rungekutta .eq. 10) then
!$omp do
         do i = 1, kmaxe
            call calculate_jacobian_3d_mf_loop1_cell(n, i)
         end do
!$omp end do
      else
!$omp do
         do i = 1, kmaxe
            call calculate_jacobian_3d_mf_loop2_cell(n, i)
         end do
!$omp end do
      end if

      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

         if (turbulence .eq. 1) call sources_derivatives_computation(n)
         if (rungekutta .eq. 10) then
!$omp do
            do i = 1, kmaxe
               call calculate_jacobian_3d_mf_loop3_cell(n, i)
            end do
!$omp end do
         else
!$omp do
            do i = 1, kmaxe
               call calculate_jacobian_3d_mf_loop4_cell(n, i)
            end do
!$omp end do

         end if
      end if

   end subroutine calculate_jacobian_3d_mf

   subroutine calculate_jacobian_3d_mf_loop1_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list

      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      impdiag_mf(i) = (impdiag_mf(i)) + (ielem_totvolume(i) / ielem_dtl(i))

      if ((realgas .eq. 1) .and. (rg_relax .eq. 2)) then
         ! Pack the strided state before evaluating the scalar source bound.
         leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         impdiag_mf(i) = impdiag_mf(i) + ielem_totvolume(i) * &
            & realgas_source_radius(n, leftv)
      end if
   end subroutine calculate_jacobian_3d_mf_loop1_cell

   subroutine calculate_jacobian_3d_mf_loop2_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list

      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      impdiag_mf(i) = ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiag_mf(i))
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2)) then
         ! Pack the strided state before evaluating the scalar source bound.
         leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         impdiag_mf(i) = impdiag_mf(i) + ielem_totvolume(i) * &
            & realgas_source_radius(n, leftv)
      end if
   end subroutine calculate_jacobian_3d_mf_loop2_cell

   subroutine calculate_jacobian_3d_mf_loop3_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list

      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & ielem_dtl(i))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & ielem_dtl(i))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
         end do
      end if
   end subroutine calculate_jacobian_3d_mf_loop3_cell

   subroutine calculate_jacobian_3d_mf_loop4_cell(n, i)
      implicit none
      integer, intent(in) :: i
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list

      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      if (turbulence .eq. 1) then
         do nvar = 1, turbulenceequations
!
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + ielem_totvolume(i) * ((1.00 / &
               & ielem_dtl(i)) + (1.50 / dt))))
            impdiagt(i, nvar) = max((ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt)) + (impdiagt(i, nvar)) - sht(i, &
               & nvar)), max(turb_diag_floor_frac * (ielem_totvolume(i) * ((1.00 / ielem_dtl(i)) + (1.50 / dt))), &
               & turb_diag_abs_floor))
         end do
      end if
      if (passivescalar .gt. 0) then
         do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            sht(i, nvar) = min(max(sht(i, nvar), 0.00), turb_source_cap_frac * (impdiagt(i, nvar) + (ielem_totvolume(i) / &
               & ielem_dtl(i))))
            impdiagt(i, nvar) = max((impdiagt(i, nvar) + (ielem_totvolume(i) / ielem_dtl(i)) - sht(i, nvar)), &
               & max(turb_diag_floor_frac * (ielem_totvolume(i) / ielem_dtl(i)), turb_diag_abs_floor))
         end do
      end if

   end subroutine calculate_jacobian_3d_mf_loop4_cell

   subroutine calculate_jacobian_3d_mf_inner_cell(n, ii)
      implicit none
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc, kas
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list

      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      i = el_int(ii)
      iconsidered = i
      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
         impdiagt(i, 1:turbulenceequations + passivescalar) = 0.0
         impofft(i, :, :) = zero
      end if
      do l = 1, ielem_ifca(i) !for all their faces
         b_code = 0
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = angle1
         ny = angle2
         mul1 = ielem_surf(i, l)

         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
         cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)

!                                      cleft(1:nof_variables)=rec_uleft(i,1:nof_variables,l,1)
!                                    cright(1:nof_variables)=rec_uleft(ielem_ineigh(i,l),1:nof_variables,ielem_ineighn(i,l),1)

         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            !left additional equations flow state
            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
            cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, 1:turbulenceequations + passivescalar)

         end if

         call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotateb(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
         asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
            viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
               & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
               & * viscots / (ielem_dih(i, l) * prandtl)))
            if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
            vpp = max(asound1, asound2) + viscots

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(3, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

               !                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if
         end if

         impdiag_mf(i) = impdiag_mf(i) + (oo2 * vpp * mul1)
         impoff_mf(i, l) = vpp

         if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations + passivescalar
               impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * vpp * mul1)
               impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * vpp * mul1)
            end do
         end if

      end do
   end subroutine calculate_jacobian_3d_mf_inner_cell

   subroutine calculate_jacobian_3d_mf_bound_cell(n, ii)
      implicit none
      integer, intent(in) :: ii
      !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2
      integer :: i, l, ngp, kmaxe, iqp, nvar, n_node, ibfc, kas
      real :: sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: identity1
      real, dimension(gpu_max_nvar, gpu_max_nvar) :: convj, diffj
      integer :: iconsidered, facex, pointx, igoflux, nf, lf, rowf
      integer :: b_code
      real :: angle1, angle2, nx, ny, nz
      real, dimension(1:gpu_max_nvar) :: cleft, cright
      real, dimension(1:gpu_max_nvar_total) :: cright_rot, cleft_rot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_nvar) :: leftv, srf_speedrot, srf_speed
      real, dimension(1:gpu_max_nvar) :: rightv
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list

      real, dimension(1:4) :: viscl, laml
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm
      real, dimension(1:20) :: eddyfl, eddyfr

      real, dimension(3) :: cords
      real :: mp_pinfl, gammal
      real :: mp_pinfr, gammar
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: eigvl

      i = el_bnd(ii)
      iconsidered = i

      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
         impdiagt(i, 1:turbulenceequations + passivescalar) = zero
         impofft(i, :, :) = zero
      end if
      do l = 1, ielem_ifca(i)
         mul1 = ielem_surf(i, l)
         angle1 = ielem_faceanglex(i, l)
         angle2 = ielem_faceangley(i, l)
         nx = angle1
         ny = angle2

         b_code = 0
         cleft(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
!                                       cleft(1:nof_variables)=rec_uleft(i,1:nof_variables,l,1)
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

            cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)

         end if

         if (ielem_ineighb(i, l) .eq. n) then !my cpu only
            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)
!                                                                   cright(1:nof_variables)=rec_uleft(ielem_ineigh(i,l),1:nof_variables,ielem_ineighn(i,l),1)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)

                  end if

                  kas = 1

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

                  leftv(1:nof_variables) = cleft(1:nof_variables)
                  b_code = ibound_icode(ielem_ibounds(i, l))

                  call boundarys(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, cturbl, &
                     & cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                  cright(1:nof_variables) = rightv(1:nof_variables)

                  kas = 2

               end if
            else
               cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 1, 1:nof_variables)
!                                                               cright(1:nof_variables)=rec_uleft(ielem_ineigh(i,l),1:nof_variables,ielem_ineighn(i,l),1)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                     & 1:turbulenceequations + passivescalar)

               end if

               kas = 3

            end if
         else !in other cpus they can only be periodic or mpi neighbours

            if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
               if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                               cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

!
                  kas = 4

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

!                                                               cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                     nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                     lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                     rowf = halo_offset(nf) + lf - 1
                     cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                        & turbulenceequations + passivescalar)

                  end if

               end if
            else
               kas = 5

!
!                                                               cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

!
               nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
               lf = rec_ihexl(i, 1, ielem_indexi(i, L))
               rowf = halo_offset(nf) + lf - 1
               cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

!                                                               cturbr(1:turbulenceequations+passivescalar)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                                               (rec_ihexl(1,ielem_indexi(l,i)),nof_variables+1:nof_variables+turbulenceequations+passivescalar)

                  nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                  lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                  rowf = halo_offset(nf) + lf - 1
                  cturbr(1:turbulenceequations + passivescalar) = solhir(rowf, nof_variables + 1:nof_variables + &
                     & turbulenceequations + passivescalar)

               end if

!
            end if
         end if

         call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
         call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

         if ((lmach .eq. 1)) then !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb(n, cright, cright_rot, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotateb(n, cleft, cleft_rot, angle1, angle2)

         end if

         leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
         call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

         asound1 = sqrt(leftv(5) * gamma / leftv(1)) + abs(cleft_rot(2) / cleft_rot(1))
         asound2 = sqrt(rightv(5) * gamma / rightv(1)) + abs(cright_rot(2) / cright_rot(1))

         vpp = max(asound1, asound2)

         if (itestcase .eq. 4) then
            leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
            viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0
            call get_visc_conduct(n, leftv, rightv, viscl, laml)

            viscots = (viscl(1) + viscl(2)) * oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
            viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots &
               & / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) &
               & * viscots / (ielem_dih(i, l) * prandtl)))
            if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)
            vpp = max(asound1, asound2) + viscots

            if (turbulence .eq. 1) then
               if (turbulencemodel .eq. 1) then
                  turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if
               if (turbulencemodel .eq. 2) then
                  eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                  call fill_turbulence_gradients(3, i, eddyfl)


                  eddyfr = eddyfl
                  call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
               end if

               viscots = oo2 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem_dih(i,l))
               viscots = max((4.0 / (3.0 * 0.5 * (cleft(1) + cright(1)))) * viscots * &
                  & 1.0e0 / ielem_dih(i, l), ((gamma / (0.5 * (cleft(1) + cright(1)))) * &
                  & viscots / (ielem_dih(i, l) * (prandtl + prtu))))

               if (realgas .eq. 1) viscots = realgas_face_diffusion_speed(cleft, cright, ielem_dih(i, l), viscl, laml)

               vpp = max(asound1, asound2) + viscots
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               if (turbulence .eq. 1) then
                  do nvar = 1, turbulenceequations

                     vpp = max(asound1, asound2) + viscots
!

                  end do
               end if
               if (passivescalar .gt. 0) then
                  do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                     viscl(1) = viscl(1) / schmidt_lam
                     viscl(2) = viscl(2) / schmidt_lam

                     if (turbulence .eq. 1) then
                        viscl(3) = viscl(3) / schmidt_turb
                        viscl(4) = viscl(4) / schmidt_turb
                        viscots = 0.5 * ((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                     else
                        viscots = 0.5 * ((viscl(1) + viscl(2)))

                     end if

                     viscots = (2.0 * viscots) / ((cleft(1) + cright(1)) * ielem_dih(i, l))
                     vpp = max(asound1, asound2) + viscots

                  end do
               end if
            end if
         else

         end if

         impdiag_mf(i) = impdiag_mf(i) + (oo2 * vpp * mul1)
         impoff_mf(i, l) = vpp
         if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations + passivescalar
               impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2 * vpp * mul1)
               impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2 * vpp * mul1)
            end do
         end if

      end do
   end subroutine calculate_jacobian_3d_mf_bound_cell

end module implicit_fluxes
