module fluxes
   use library
   use transform
   use local
   use riemann
   use flow_operations
   use declaration
   use basis
   use dg_functions
   implicit none

contains

   subroutine subtract_volume_fraction_source(i, divergence, dg_storage)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      real, intent(in) :: divergence
      logical, intent(in) :: dg_storage
      integer :: v
      if (multispecies /= 1 .or. mp_modelc /= 0) return
      do v = dimensiona + 3 + nof_species, nof_variables
         if (dg_storage) then
            rhs_valdg(i, 1, v) = rhs_valdg(i, 1, v) - u_c_val(i, 1, v) * divergence
         else
            rhs_val(i, v) = rhs_val(i, v) - u_c_val(i, 1, v) * divergence
         end if
      end do
   end subroutine subtract_volume_fraction_source


   subroutine solution_integ(i, solution_integ2)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(inout) :: solution_integ2
      integer, intent(in) :: i

      call dg_vol_integral2_fill(n, i, solution_integ2)

   end subroutine solution_integ

   subroutine solution_integ_s(i, solution_integ_strong)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(inout) :: solution_integ_strong
      integer, intent(in) :: i

      call dg_vol_integral_strong_fill(n, i, solution_integ_strong)

   end subroutine solution_integ_s

   subroutine solution_integ_w(i, solution_integ_weak)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:nof_variables), intent(inout) :: solution_integ_weak
      integer, intent(in) :: i

      call dg_vol_integral_weak_fill(n, i, solution_integ_weak)

   end subroutine solution_integ_w

   subroutine calculate_fluxeshi(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation
      implicit none
      integer, intent(in) :: n
      integer :: i, kmaxe

      kmaxe = xmpielrank(n)

      if (dg .eq. 1) then
!$omp barrier
!$omp do
         do i = 1, kmaxe
            call calculate_fluxeshi_cell(n, i)
         end do
!$omp end do
      else
         call calculate_fluxeshi_cell_fv(n)
      end if

   end subroutine calculate_fluxeshi

   subroutine calculate_fluxeshi_cell_fv(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, kmaxe
      real :: godflux2
      integer :: l, ngp, iqp, nfx, lfx, rowfx
      real :: angle1, angle2, normalvect
      real :: cleft1, cright1, hllcflux1

      kmaxe = xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n, kmaxe) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, zero, lamx, lamy, lamz, qp_quad) &
!$omp& firstprivate(qp_triangle) &
!$omp& map(alloc: weights_q, weights_t, ielem_interior, ielem_ifca) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_types_faces, ielem_surf) &
!$omp& map(alloc: ielem_ineighn, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode) &
!$omp& map(alloc: ielem_qface, ielem_inter_id, ielem_indexf, bound_offset, boundhir, rec_uleft, rhs_val) &
!$omp& private(i, l, ngp, iqp, nfx, lfx, rowfx, angle1, angle2, normalvect, cleft1, cright1, hllcflux1, godflux2)
#else
!$omp do &
!$omp& private(i, l, ngp, iqp, nfx, lfx, rowfx, angle1, angle2, normalvect, cleft1, cright1, hllcflux1, godflux2)
#endif
      do i = 1, kmaxe
         block
            real, dimension(1:gpu_max_qp_face) :: weights_temp

            rhs_val(i, 1) = zero

            do l = 1, ielem_ifca(i)
               godflux2 = zero
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               normalvect = ((cos(angle1) * sin(angle2)) * lamx) + ((sin(angle1) * sin(angle2)) * lamy) + ((cos(angle2)) * lamz)
               if (ielem_types_faces(i, l) .eq. 5) then
                  iqp = qp_quad
                  weights_temp(1:iqp) = weights_q(1:iqp)
               else
                  iqp = qp_triangle
                  weights_temp(1:iqp) = weights_t(1:iqp)
               end if

               do ngp = 1, iqp
                  cleft1 = rec_uleft(i, 1, l, ngp)

                  if (ielem_interior(i) .eq. 0) then
                     cright1 = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                  else
                     if (ielem_ineighb(i, l) .eq. n) then
                        if (ielem_ibounds(i, l) .gt. 0) then
                           if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then
                              cright1 = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                           else
                              cright1 = cleft1
                           end if
                        else
                           cright1 = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                        end if
                     else
                        if (ielem_ibounds(i, l) .gt. 0) then
                           if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then
                              nfx = ielem_ineighn(i, l)
                              lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                              rowfx = bound_offset(nfx) + lfx - 1
                              cright1 = boundhir(rowfx, 1)
                           else
                              cright1 = cleft1
                           end if
                        else
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cright1 = boundhir(rowfx, 1)
                        end if
                     end if
                  end if

                  hllcflux1 = zero
                  if (normalvect .gt. zero) then
                     hllcflux1 = normalvect * cleft1
                  else
                     hllcflux1 = normalvect * cright1
                  end if
                  godflux2 = godflux2 + (hllcflux1 * (weights_temp(ngp) * ielem_surf(i, l)))
               end do
               rhs_val(i, 1) = rhs_val(i, 1) + godflux2
            end do
         end block
      end do

#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine calculate_fluxeshi_cell_fv

   subroutine calculate_fluxeshi_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      integer, intent(in) :: n
      integer :: kmaxe
      real :: godflux2
      integer :: l, ngp, iqp
      real, dimension(1:gpu_max_qp_face) :: weights_temp
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real :: angle1, angle2, nx, ny, nz, normalvect
      integer :: facex, pointx, iconsidered, nfx, lfx, rowfx
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, hllcflux, rhllcflux

      iconsidered = i

      if (dg .eq. 1) then
         rhs_valdg(i, :, :) = 0.00
      end if

      if (dg .eq. 1) then
         call dg_vol_integral_accumulate_rhs(n, iconsidered, -1.00)

      end if

      if (ielem_interior(i) .eq. 0) then
         rhs_val(i, 1) = zero

         do l = 1, ielem_ifca(i)
            godflux2 = zero
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            facex = l
            normalvect = ((cos(angle1) * sin(angle2)) * lamx) + ((sin(angle1) * sin(angle2)) * lamy) + ((cos(angle2)) * lamz)
            if (ielem_types_faces(i, l) .eq. 5) then
               iqp = qp_quad
               weights_temp(1:iqp) = weights_q(1:iqp)

            else
               iqp = qp_triangle
               weights_temp(1:iqp) = weights_t(1:iqp)
            end if

            do ngp = 1, iqp
               pointx = ngp
               if (dg .eq. 1) then
                  cleft = rec_uleft_dg(i, 1:nof_variables, l, ngp)
                  cright = rec_uleft_dg(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
               else !fv

                  cleft(1) = rec_uleft(i, 1, l, ngp)
                  cright(1) = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
               end if
               call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)

               if (dg .eq. 1) then
                  rhllcflux(1) = hllcflux(1)

                  call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

               else
                  godflux2 = godflux2 + (hllcflux(1) * (weights_temp(ngp) * ielem_surf(i, l)))
               end if
            end do
            rhs_val(i, 1) = rhs_val(i, 1) + godflux2

         end do
      end if
      if (ielem_interior(i) .eq. 1) then
         rhs_val(i, 1) = zero

         do l = 1, ielem_ifca(i)
            facex = l
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            normalvect = ((cos(angle1) * sin(angle2)) * lamx) + ((sin(angle1) * sin(angle2)) * lamy) + ((cos(angle2)) * lamz)

            if (ielem_types_faces(i, l) .eq. 5) then
               iqp = qp_quad
               weights_temp(1:iqp) = weights_q(1:iqp)
            else
               iqp = qp_triangle
               weights_temp(1:iqp) = weights_t(1:iqp)
            end if

            godflux2 = zero
            do ngp = 1, iqp
               pointx = ngp
               if (dg == 1) then
                  cleft(1:nof_variables) = rec_uleft_dg(i, 1:nof_variables, l, ngp)
               else
                  cleft(1:nof_variables) = rec_uleft(i, 1:nof_variables, l, ngp)
               end if
               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries

                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        if (dg == 1) then
                           cright = rec_uleft_dg(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
                        else
                           cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
                        end if
                     else
                        !not periodic ones in my cpu
                        cright(1:nof_variables) = cleft(1:nof_variables)

                     end if
                  else
                     if (dg == 1) then
                        cright(1:nof_variables) = rec_uleft_dg(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)

                     else
                        cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
                     end if
                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu
                        if (dg == 1) then
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cright(1:nof_variables) = boundhir_dg(rowfx, 1:nof_variables)
                        else
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)
                        end if
                     end if
                  else
                     if (dg == 1) then
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                        rowfx = bound_offset(nfx) + lfx - 1
                        cright(1:nof_variables) = boundhir_dg(rowfx, 1:nof_variables)
                     else
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                        rowfx = bound_offset(nfx) + lfx - 1
                        cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)
                     end if
!
                  end if
               end if

               call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)

               if (dg .eq. 1) then
                  rhllcflux(1) = hllcflux(1)
                  call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

               else

                  godflux2 = godflux2 + (hllcflux(1) * (weights_temp(ngp) * ielem_surf(i, l)))

               end if
            end do
            rhs_val(i, 1) = rhs_val(i, 1) + godflux2

         end do
      end if
!
   end subroutine calculate_fluxeshi_cell

   subroutine calculate_fluxeshi2d(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation in 2d
      implicit none
      integer, intent(in) :: n
      integer :: i, kmaxe

      kmaxe = xmpielrank(n)

      if (dg .eq. 1) then
!$omp barrier
!$omp do
         do i = 1, kmaxe
            call calculate_fluxeshi2d_cell(n, i)
         end do
!$omp end do
      else
         call calculate_fluxeshi2d_cell_fv(n)
      end if

   end subroutine calculate_fluxeshi2d

   subroutine calculate_fluxeshi2d_cell_fv(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, kmaxe
      real :: godflux2, lamxl, lamyl
      integer :: l, ngp, nfx, lfx, rowfx
      real :: nx, ny, normalvect
      real :: cleft1, cright1, hllcflux1

      kmaxe = xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n, kmaxe) &
!$omp& firstprivate(nof_variables, turbulenceequations, passivescalar, zero, lamx, lamy, initcond, qp_line_n) &
!$omp& map(alloc: weights_l, ielem_interior, ielem_ifca, ielem_xxc, ielem_yyc, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_surf, ielem_ineighn, ielem_ineigh, ielem_ineighb) &
!$omp& map(alloc: ielem_ibounds, ibound_icode, ielem_qface, ielem_inter_id, ielem_indexf) &
!$omp& map(alloc: bound_offset, boundhir, rec_uleft, rhs_val) &
!$omp& private(i, l, ngp, nfx, lfx, rowfx, nx, ny, normalvect, cleft1, cright1, hllcflux1, godflux2, lamxl, lamyl)
#else
!$omp do &
!$omp& private(i, l, ngp, nfx, lfx, rowfx, nx, ny, normalvect, cleft1, cright1, hllcflux1, godflux2, lamxl, lamyl)
#endif
      do i = 1, kmaxe
         block
            real, dimension(1:gpu_max_qp_face) :: weights_temp

            weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)

            if (initcond .eq. 3) then
               lamxl = -ielem_yyc(i) + 0.50
               lamyl = ielem_xxc(i) - 0.50
            else
               lamxl = lamx
               lamyl = lamy
            end if

            rhs_val(i, 1) = zero

            do l = 1, ielem_ifca(i)
               godflux2 = zero
               nx = ielem_faceanglex(i, l)
               ny = ielem_faceangley(i, l)
               normalvect = (nx * lamxl) + (ny * lamyl)

               do ngp = 1, qp_line_n
                  cleft1 = rec_uleft(i, 1, l, ngp)

                  if (ielem_interior(i) .eq. 0) then
                     cright1 = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                  else
                     if (ielem_ineighb(i, l) .eq. n) then
                        if (ielem_ibounds(i, l) .gt. 0) then
                           if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then
                              cright1 = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                           else
                              cright1 = cleft1
                           end if
                        else
                           cright1 = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                        end if
                     else
                        if (ielem_ibounds(i, l) .gt. 0) then
                           if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then
                              nfx = ielem_ineighn(i, l)
                              lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                              rowfx = bound_offset(nfx) + lfx - 1
                              cright1 = boundhir(rowfx, 1)
                           else
                              cright1 = cleft1
                           end if
                        else
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cright1 = boundhir(rowfx, 1)
                        end if
                     end if
                  end if

                  hllcflux1 = zero
                  if (normalvect .gt. zero) then
                     hllcflux1 = normalvect * cleft1
                  else
                     hllcflux1 = normalvect * cright1
                  end if
                  godflux2 = godflux2 + (hllcflux1 * (weights_temp(ngp) * ielem_surf(i, l)))
               end do
               rhs_val(i, 1) = rhs_val(i, 1) + godflux2
            end do
         end block
      end do

#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine calculate_fluxeshi2d_cell_fv

   subroutine calculate_fluxeshi2d_cell(n, i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      integer, intent(in) :: n
      integer :: kmaxe
      real :: godflux2, sum_detect, lamxl, lamyl
      integer :: l, k, ngp, iqp, neighbor_index, neighbor_face_index
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real :: angle1, angle2, nx, ny, nz, normalvect
      integer :: facex, pointx, iconsidered, nfx, lfx, rowfx
      real, dimension(1:gpu_max_qp_face) :: weights_temp
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, hllcflux, rhllcflux

      weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)

      if (initcond .eq. 3) then
         lamxl = -ielem_yyc(i) + 0.50
         lamyl = ielem_xxc(i) - 0.5
      else
         lamxl = lamx
         lamyl = lamy
      end if

      if (dg .eq. 1) then
         rhs_valdg(i, :, :) = zero
      else
         rhs_val(i, :) = zero
      end if

      iconsidered = i

      if (dg .eq. 1) then

         call dg_vol_integral_accumulate_rhs(n, iconsidered, -1.00)

      end if

      if (ielem_interior(i) .eq. 0) then ! element is interior

         do l = 1, ielem_ifca(i)
            godflux2 = zero
            nx = ielem_faceanglex(i, l)
            ny = ielem_faceangley(i, l)
            facex = l

            normalvect = (nx * lamxl) + (ny * lamyl)

            iqp = qp_line_n

            neighbor_index = ielem_ineigh(i, l)
            neighbor_face_index = ielem_ineighn(i, l)

            do ngp = 1, iqp
               pointx = ngp

               if (dg .eq. 1) then
                  cleft = rec_uleft_dg(i, 1:nof_variables, l, ngp)
                  cright = rec_uleft_dg(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
               else !fv
                  cleft(1) = rec_uleft(i, 1, l, ngp)
                  cright(1) = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
               end if

               call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)

               if (dg .eq. 1) then
                  ! riemann flux at interface quadrature points times basis
                  rhllcflux(1) = hllcflux(1)

                  call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

               else !fv
                  godflux2 = godflux2 + (hllcflux(1) * (weights_temp(ngp) * ielem_surf(i, l)))
               end if
            end do

            if (dg /= 1) rhs_val(i, 1) = rhs_val(i, 1) + godflux2

         end do

      else if (ielem_interior(i) .eq. 1) then

         do l = 1, ielem_ifca(i)
            facex = l
            nx = ielem_faceanglex(i, l)
            ny = ielem_faceangley(i, l)
            normalvect = (nx * lamxl) + (ny * lamyl)
            iqp = qp_line_n

            godflux2 = zero
            do ngp = 1, iqp
               pointx = ngp

               if (dg == 1) then
                  cleft = rec_uleft_dg(i, 1:nof_variables, l, ngp)
               else
                  cleft(1) = rec_uleft(i, 1, l, ngp)
               end if

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        if (dg == 1) then
                           cright = rec_uleft_dg(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
                        else !fv
                           cright(1) = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                        end if
                     else !not periodic ones in my cpu
                        cright(1:nof_variables) = cleft(1:nof_variables)
                     end if
!
                  else
                     if (dg == 1) then
                        cright = rec_uleft_dg(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
!
                     else !fv
                        cright(1) = rec_uleft(ielem_ineigh(i, l), 1, ielem_ineighn(i, l), ngp)
                     end if
                  end if
               else !in other cpus they can only be periodic or mpi neighbours
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

                        if (dg == 1) then
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cright(1:nof_variables) = boundhir_dg(rowfx, 1:nof_variables)
                        else
!                                     cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)
                        end if
                     end if
                  else
                     if (dg == 1) then
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol_dg(ielem_qface(l,ngp,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                        rowfx = bound_offset(nfx) + lfx - 1
                        cright(1:nof_variables) = boundhir_dg(rowfx, 1:nof_variables)
                     else
!                                 cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                        rowfx = bound_offset(nfx) + lfx - 1
                        cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)
                     end if
                  end if

               end if

               call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)

               if (dg .eq. 1) then
                  !riemann flux at interface quadrature points times basis
                  rhllcflux(1) = hllcflux(1)
                  call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

               else !fv
                  godflux2 = godflux2 + (hllcflux(1) * (weights_temp(ngp) * ielem_surf(i, l)))
               end if

            end do

            if (dg /= 1) rhs_val(i, 1) = rhs_val(i, 1) + godflux2
         end do

      end if

   end subroutine calculate_fluxeshi2d_cell

   subroutine calculate_fluxeshi_convective(n)
      implicit none
      integer, intent(in) :: n

      call calculate_fluxeshi_convective_inner_cell(n)
      if (dg .eq. 1) then
         call calculate_fluxeshi_convective_bound_cell_dg_volume(n)
         call calculate_fluxeshi_convective_bound_cell_dg_surface(n)
      else
         call calculate_fluxeshi_convective_bound_cell_fv(n)
      end if
   end subroutine calculate_fluxeshi_convective

   subroutine calculate_fluxeshi_convective_inner_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, b_code, srf, iv, nvt
      real :: norms, wface
      integer :: iconsidered, facex, pointx
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(dg, icoupleturb, iriemann, lmach, mrf, multispecies, n, nof_interior, nof_variables) &
!$omp& firstprivate(passivescalar, qp_quad_n, qp_triangle_n, turbulence, turbulenceequations, zero) &
!$omp& map(alloc: el_int, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, ielem_types_faces, rec_mrf) &
!$omp& map(alloc: rhs_val, rhs_valdg, rhst_val, weights_q, weights_t) &
!$omp& private(angle1, angle2, b_code, facex) &
!$omp& private(i, iconsidered, ii, iqp, iv, l, mp_source1, mp_source2, mp_source3, ngp, norms, nvt) &
!$omp& private(nx, ny, nz, pointx, srf, wface) &
!$omp& map(alloc: gamma_in, ielem_facediss, ielem_idegfree, ielem_ineigh, ielem_ineighn, ielem_iorder) &
!$omp& map(alloc: ielem_itotalpoints, ielem_totvolume, ielem_xxc, ielem_yyc, integ_basis_dg_value, integ_basis_value) &
!$omp& map(alloc: integ_basis_valuec, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in, mp_tmid_in) &
!$omp& map(alloc: qp_array_qp_weight, qp_array_x, qp_array_y, qp_array_z, rec_br2_local_lift, rec_qpoints, rec_rotvel) &
!$omp& map(alloc: rec_uleft, rec_uleft_dg, rec_uleftturb, rec_volume, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb) &
!$omp& map(alloc: u_c_br2_aux_var, u_c_val, u_c_valdg, u_ct_val) &
!$omp& map(alloc: rg_thetag)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, facex) &
!$omp& private(i, iconsidered, ii, iqp, iv, l, mp_source1, mp_source2, mp_source3, ngp, norms, nvt) &
!$omp& private(nx, ny, nz, pointx, srf, wface)
#endif
      do ii = 1, nof_interior
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:gpu_max_dim) :: pox, poy, poz

            i = el_int(ii)

            nvt = turbulenceequations + passivescalar
            do iv = 1, nof_variables
               rhs_val(i, iv) = zero
            end do
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               do iv = 1, nvt
                  rhst_val(i, iv) = zero
               end do
            end if

            iconsidered = i
            mp_source3 = zero

            if (mrf .eq. 1) then
               srf = rec_mrf(i)
            end if
            if (dg .eq. 1) then
               rhs_valdg(i, :, :) = zero

               call dg_vol_integral_accumulate_rhs(n, iconsidered, -1.00)

            end if

            do l = 1, ielem_ifca(i) !for all their faces
               b_code = 0

               do iv = 1, nof_variables + nvt
                  godflux2(iv) = zero
               end do
               mp_source2 = zero
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))
               if (ielem_types_faces(i, l) .eq. 5) then
                  iqp = qp_quad_n
                  weights_temp(1:iqp) = weights_q(1:iqp)
               else
                  iqp = qp_triangle_n
                  weights_temp(1:iqp) = weights_t(1:iqp)
               end if

               facex = l

               do ngp = 1, iqp !for all the gaussian quadrature points
                  pointx = ngp

                  call get_states_interior(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, angle1, &
               & angle2, nx, ny, &
                     & nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright)

                  !rotate wrt to normalvector of face and solve 1d riemann problem
                  call rotatef(n, cright_rot, cright, angle1, angle2)
                  call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                     call lmacht(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                  end if

                  select case (iriemann)

                  case (1) !hllc

                     call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb(n, rhllcflux, hllcflux, angle1, angle2)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        do iv = 1, nvt
                           rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
                        end do
                     end if

                  case (9) !hll

                     call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        do iv = 1, nvt
                           rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
                        end do
                     end if

                  case (2) !rusanov

                     call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb(n, rhllcflux, hllcflux, angle1, angle2)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        do iv = 1, nvt
                           rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
                        end do
                     end if
                  case (3) !roe

                     call rotateb(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb(n, cright, cright_rot, angle1, angle2)
                     call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                     do iv = 1, nof_variables + nvt
                        rhllcflux(iv) = hllcflux(iv)
                     end do

                  case (4) !roe

                     call rotateb(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb(n, cright, cright_rot, angle1, angle2)
                     call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                     do iv = 1, nof_variables + nvt
                        rhllcflux(iv) = hllcflux(iv)
                     end do

                  case (5) !roe

                     call rotateb(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb(n, cright, cright_rot, angle1, angle2)
                     call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                     do iv = 1, nof_variables + nvt
                        rhllcflux(iv) = hllcflux(iv)
                     end do

                  end select

                  wface = weights_temp(ngp) * ielem_surf(i, l)
                  if (dg .eq. 1) then

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

                  else

                     do iv = 1, nof_variables
                        godflux2(iv) = godflux2(iv) + (rhllcflux(iv) * wface)
                     end do

                  end if

                  if (multispecies .eq. 1) then
                     mp_source2 = mp_source2 + mp_source1 * wface
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (icoupleturb .eq. 0) then !first order upwind flux

                        norms = 0.50 * (cleft_rot(2) / cleft_rot(1) + cright_rot(2) / cright_rot(1))

                        if (rec_mrf(iconsidered) .eq. 1) then
                           norms = norms - srf_speedrot(2)
                        end if
                        do iv = 1, nvt
                           rhllcflux(nof_variables + iv) = 0.50 * ((norms * (cturbl(iv) + cturbr(iv))) + (abs(norms) * &
                              & (cturbl(iv) - &
                              & cturbr(iv))))
                        end do

                     end if

                     do iv = 1, nvt
                        godflux2(nof_variables + iv) = godflux2(nof_variables + iv) + (rhllcflux(nof_variables + iv) * wface)
                     end do

                  end if

               end do

               do iv = 1, nof_variables
                  rhs_val(i, iv) = rhs_val(i, iv) + godflux2(iv)
               end do
               if (multispecies .eq. 1) then
                  mp_source3 = mp_source3 + mp_source2
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  do iv = 1, nvt
                     rhst_val(i, iv) = rhst_val(i, iv) + godflux2(nof_variables + iv)
                  end do

               end if

            end do
            if (multispecies .eq. 1) then
               call subtract_volume_fraction_source(i, mp_source3, .false.)

            end if

            if (dg .eq. 1) then
               if (multispecies .eq. 1) then

                  call subtract_volume_fraction_source(i, mp_source3, .true.)

               end if

            end if

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_convective_inner_cell

   subroutine calculate_fluxeshi_convective_bound_point(n, iconsidered, facex, pointx, rhllcflux, mp_source1, b_code)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered, facex, pointx
      integer, intent(out) :: b_code
      integer :: iv, nvt
      real, intent(out) :: mp_source1
      real, dimension(1:gpu_max_nvar_total), intent(out) :: rhllcflux
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot, hllcflux
      real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real :: angle1, angle2, nx, ny, nz, norms

      nvt = turbulenceequations + passivescalar
      b_code = 0
      mp_source1 = zero
      rhllcflux = zero
      hllcflux = zero
      srf_speedrot = zero

      angle1 = ielem_faceanglex(iconsidered, facex)
      angle2 = ielem_faceangley(iconsidered, facex)
      nx = (cos(angle1) * sin(angle2))
      ny = (sin(angle1) * sin(angle2))
      nz = (cos(angle2))

      call get_states_bounds(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
         & cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright)

      if (b_code == 9 .and. duct_outlet_face_flux_enabled()) then
         call pressure_outlet_face_flux(n, cright(1:nof_variables), cturbr(1:nvt), [nx, ny, nz], rhllcflux, mp_source1)
         return
      end if

      call rotatef(n, cright_rot, cright, angle1, angle2)
      call rotatef(n, cleft_rot, cleft, angle1, angle2)

      if ((lmach .eq. 1)) then
         leftv(1:nof_variables) = cleft_rot(1:nof_variables)
         rightv(1:nof_variables) = cright_rot(1:nof_variables)
         call lmacht(n, leftv, rightv)
         cleft_rot(1:nof_variables) = leftv(1:nof_variables)
         cright_rot(1:nof_variables) = rightv(1:nof_variables)
      end if

      if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
         do iv = 1, nvt
            cleft_rot(nof_variables + iv) = cturbl(iv)
            cright_rot(nof_variables + iv) = cturbr(iv)
         end do
      end if

      select case (iriemann)
      case (1)
         call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
         call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            do iv = 1, nvt
               rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
            end do
         end if
      case (9)
         call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
         call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            do iv = 1, nvt
               rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
            end do
         end if
      case (2)
         call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
         call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            do iv = 1, nvt
               rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
            end do
         end if
      case (3)
         call rotateb(n, cleft, cleft_rot, angle1, angle2)
         call rotateb(n, cright, cright_rot, angle1, angle2)
         call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
         do iv = 1, nof_variables + nvt
            rhllcflux(iv) = hllcflux(iv)
         end do
      case (4)
         if (b_code .gt. 0) then
            call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               do iv = 1, nvt
                  rhllcflux(nof_variables + iv) = hllcflux(nof_variables + iv)
               end do
            end if
         else
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            do iv = 1, nof_variables + nvt
               rhllcflux(iv) = hllcflux(iv)
            end do
         end if
      case (5)
         call rotateb(n, cleft, cleft_rot, angle1, angle2)
         call rotateb(n, cright, cright_rot, angle1, angle2)
         if (b_code .le. 0) then
            call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
         else
            call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
         end if
         do iv = 1, nof_variables + nvt
            rhllcflux(iv) = hllcflux(iv)
         end do
      end select

      if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
         if (icoupleturb .eq. 0) then
            norms = 0.50 * (cleft_rot(2) / cleft_rot(1) + cright_rot(2) / cright_rot(1))
            if (rec_mrf(iconsidered) .eq. 1) then
               norms = norms - srf_speedrot(2)
            end if
            do iv = 1, nvt
               rhllcflux(nof_variables + iv) = 0.50 * ((norms * (cturbl(iv) + cturbr(iv))) + (abs(norms) * (cturbl(iv) - &
                  & cturbr(iv))))
            end do
         end if
         if ((b_code .eq. 3) .or. (b_code .eq. 4)) then
            do iv = 1, nvt
               rhllcflux(nof_variables + iv) = zero
            end do
         end if
      end if
   end subroutine calculate_fluxeshi_convective_bound_point

   subroutine calculate_fluxeshi_convective_bound_cell_fv(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, iv, nvt, iconsidered, facex, pointx, b_code
      real :: mp_source1, mp_source2, mp_source3, wface

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(multispecies, n, nof_bounded, nof_variables, passivescalar, qp_quad_n, qp_triangle_n) &
!$omp& firstprivate(turbulence, turbulenceequations, zero) &
!$omp& map(alloc: el_bnd, ielem_ifca, ielem_surf, ielem_types_faces, rhs_val, rhst_val, weights_q, weights_t) &
!$omp& private(b_code, facex, i, iconsidered, ii, iqp, iv, l, mp_source1, mp_source2, mp_source3, ngp) &
!$omp& private(nvt, pointx, wface) &
!$omp& map(alloc: bound_offset, boundhir, boundhir_dg, catalytic_con, gamma_in, ibound_icode, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_facediss, ielem_ibounds, ielem_indexf, ielem_ineigh, ielem_ineighb) &
!$omp& map(alloc: ielem_ineighn, ielem_inter_id, ielem_nodes_faces, ielem_qface, ielem_walldist, inoder4_cord) &
!$omp& map(alloc: mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_mrf, rec_rotvel) &
!$omp& map(alloc: rec_uleft, rec_uleft_dg, rec_uleftturb, rg_hzero, rg_molm, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb) &
!$omp& map(alloc: rhs_valdg, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(b_code, facex, i, iconsidered, ii, iqp, iv, l, mp_source1, mp_source2, mp_source3, ngp) &
!$omp& private(nvt, pointx, wface)
#endif
      do ii = 1, nof_bounded
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux
            real, dimension(1:gpu_max_qp_face) :: weights_temp

            i = el_bnd(ii)
            nvt = turbulenceequations + passivescalar
            do iv = 1, nof_variables
               rhs_val(i, iv) = zero
            end do
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               do iv = 1, nvt
                  rhst_val(i, iv) = zero
               end do
            end if

            iconsidered = i
            mp_source3 = zero

            do l = 1, ielem_ifca(i)
               facex = l
               if (ielem_types_faces(i, l) .eq. 5) then
                  iqp = qp_quad_n
                  weights_temp(1:iqp) = weights_q(1:iqp)
               else
                  iqp = qp_triangle_n
                  weights_temp(1:iqp) = weights_t(1:iqp)
               end if

               do iv = 1, nof_variables + nvt
                  godflux2(iv) = zero
               end do
               mp_source2 = zero

               do ngp = 1, iqp
                  pointx = ngp
                  call calculate_fluxeshi_convective_bound_point(n, iconsidered, facex, pointx, rhllcflux, mp_source1, b_code)
                  wface = weights_temp(ngp) * ielem_surf(i, l)
                  do iv = 1, nof_variables
                     godflux2(iv) = godflux2(iv) + (rhllcflux(iv) * wface)
                  end do
                  if (multispecies .eq. 1) then
                     mp_source2 = mp_source2 + mp_source1 * wface
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     do iv = 1, nvt
                        godflux2(nof_variables + iv) = godflux2(nof_variables + iv) + (rhllcflux(nof_variables + iv) * wface)
                     end do
                  end if
               end do

               do iv = 1, nof_variables
                  rhs_val(i, iv) = rhs_val(i, iv) + godflux2(iv)
               end do
               if (multispecies .eq. 1) then
                  mp_source3 = mp_source3 + mp_source2
               end if
               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  do iv = 1, nvt
                     rhst_val(i, iv) = rhst_val(i, iv) + godflux2(nof_variables + iv)
                  end do
               end if
            end do

            if (multispecies .eq. 1) then
               call subtract_volume_fraction_source(i, mp_source3, .false.)
            end if
         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_convective_bound_cell_fv



   subroutine calculate_fluxeshi_convective_bound_cell_dg_volume(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, iv, kbasis, nvt, iconsidered, number_of_dog

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n, nof_bounded, nof_variables, passivescalar, turbulence, turbulenceequations, zero) &
!$omp& map(alloc: el_bnd, ielem_idegfree, rhs_val, rhs_valdg, rhst_val) &
!$omp& private(i, iconsidered, ii, iv, kbasis, number_of_dog, nvt) &
!$omp& map(alloc: gamma_in, ielem_iorder, ielem_itotalpoints, ielem_totvolume, ielem_xxc, ielem_yyc) &
!$omp& map(alloc: integ_basis_dg_value, integ_basis_value, integ_basis_valuec, mp_janaf, mp_m, mp_pinf, mp_thigh_in) &
!$omp& map(alloc: mp_tlow_in, mp_tmid_in, qp_array_qp_weight, qp_array_x, qp_array_y, qp_array_z, rec_br2_local_lift) &
!$omp& map(alloc: rec_volume, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb, u_c_br2_aux_var, u_c_valdg) &
!$omp& map(alloc: ielem_ifca, rg_thetag)
#else
!$omp do &
!$omp& private(i, iconsidered, ii, iv, kbasis, number_of_dog, nvt)
#endif
      do ii = 1, nof_bounded
         i = el_bnd(ii)
         nvt = turbulenceequations + passivescalar
         do iv = 1, nof_variables
            rhs_val(i, iv) = zero
         end do
         if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            do iv = 1, nvt
               rhst_val(i, iv) = zero
            end do
         end if

         iconsidered = i
         number_of_dog = ielem_idegfree(iconsidered)
         do iv = 1, nof_variables
            do kbasis = 1, number_of_dog + 1
               rhs_valdg(i, kbasis, iv) = zero
            end do
         end do
         call dg_vol_integral_accumulate_rhs(n, iconsidered, -1.00)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_convective_bound_cell_dg_volume

   subroutine calculate_fluxeshi_convective_bound_cell_dg_surface(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, iv, nvt, iconsidered, facex, pointx, b_code
      real :: mp_source1, mp_source2, mp_source3, wface

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(multispecies, n, nof_bounded, nof_variables, passivescalar, qp_quad_n, qp_triangle_n) &
!$omp& firstprivate(turbulence, turbulenceequations, zero) &
!$omp& map(alloc: el_bnd, ielem_ifca, ielem_surf, ielem_types_faces, rhst_val, weights_q, weights_t) &
!$omp& private(b_code, facex, i, iconsidered, ii, iqp, iv, l, mp_source1, mp_source2, mp_source3, ngp, nvt, pointx) &
!$omp& private(wface) &
!$omp& map(alloc: bound_offset, boundhir, boundhir_dg, catalytic_con, gamma_in, ibound_icode, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_facediss, ielem_ibounds, ielem_idegfree, ielem_indexf, ielem_ineigh) &
!$omp& map(alloc: ielem_ineighb, ielem_ineighn, ielem_inter_id, ielem_iorder, ielem_nodes_faces, ielem_qface) &
!$omp& map(alloc: ielem_totvolume, ielem_walldist, inoder4_cord, integ_basis_dg_value, integ_basis_value) &
!$omp& map(alloc: integ_basis_valuec, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in) &
!$omp& map(alloc: rec_mrf, rec_qpoints, rec_rotvel, rec_uleft, rec_uleft_dg, rec_uleftturb, rec_volume, rg_hzero) &
!$omp& map(alloc: rg_molm, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb, rhs_val, rhs_valdg, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(b_code, facex, i, iconsidered, ii, iqp, iv, l, mp_source1, mp_source2, mp_source3, ngp, nvt, pointx) &
!$omp& private(wface)
#endif
      do ii = 1, nof_bounded
         block
            real, dimension(1:gpu_max_nvar_total) :: rhllcflux
            real, dimension(1:gpu_max_extra_transport) :: turbrhs
            real, dimension(1:gpu_max_qp_face) :: weights_temp

            i = el_bnd(ii)
            nvt = turbulenceequations + passivescalar
            iconsidered = i
            mp_source3 = zero

            do l = 1, ielem_ifca(i)
               facex = l
               if (ielem_types_faces(i, l) .eq. 5) then
                  iqp = qp_quad_n
                  weights_temp(1:iqp) = weights_q(1:iqp)
               else
                  iqp = qp_triangle_n
                  weights_temp(1:iqp) = weights_t(1:iqp)
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  do iv = 1, nvt
                     turbrhs(iv) = zero
                  end do
               end if
               mp_source2 = zero

               do ngp = 1, iqp
                  pointx = ngp
                  call calculate_fluxeshi_convective_bound_point(n, iconsidered, facex, pointx, rhllcflux, mp_source1, b_code)
                  wface = weights_temp(ngp) * ielem_surf(i, l)
                  call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)
                  if (multispecies .eq. 1) then
                     mp_source2 = mp_source2 + mp_source1 * wface
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     do iv = 1, nvt
                        turbrhs(iv) = turbrhs(iv) + (rhllcflux(nof_variables + iv) * wface)
                     end do
                  end if
               end do

               if (multispecies .eq. 1) then
                  mp_source3 = mp_source3 + mp_source2
               end if
               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  do iv = 1, nvt
                     rhst_val(i, iv) = rhst_val(i, iv) + turbrhs(iv)
                  end do
               end if
            end do

            if (multispecies .eq. 1) then
               call subtract_volume_fraction_source(i, mp_source3, .false.)
            end if

            if (multispecies .eq. 1) then
               call subtract_volume_fraction_source(i, mp_source3, .true.)
            end if
         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_convective_bound_cell_dg_surface

   subroutine calculate_fluxeshi_convective2d(n)
      implicit none
      integer, intent(in) :: n

      call calculate_fluxeshi_convective2d_inner_cell(n)
      call calculate_fluxeshi_convective2d_bound_cell(n)
   end subroutine calculate_fluxeshi_convective2d

   subroutine calculate_fluxeshi_convective2d_inner_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i
      integer :: l, ngp, iqp, b_code
      real :: norms
      integer :: iconsidered, facex, pointx
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(dg, icoupleturb, iriemann, lmach, multispecies, n, nof_interior, nof_variables, passivescalar) &
!$omp& firstprivate(qp_line_n, turbulence, turbulenceequations, zero) &
!$omp& map(alloc: el_int, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, rhs_val, rhs_valdg, rhst_val) &
!$omp& map(alloc: weights_l) &
!$omp& private(angle1, angle2, b_code, facex) &
!$omp& private(i, iconsidered, ii, iqp, l, mp_source1, mp_source2, mp_source3, ngp, norms, nx, ny) &
!$omp& private(nz, pointx) &
!$omp& map(alloc: gamma_in, ielem_idegfree, ielem_ineigh, ielem_ineighn, ielem_iorder, ielem_itotalpoints) &
!$omp& map(alloc: ielem_totvolume, ielem_xxc, ielem_yyc, integ_basis_dg_value, integ_basis_value, integ_basis_valuec) &
!$omp& map(alloc: mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in, mp_tmid_in, qp_array_qp_weight, qp_array_x) &
!$omp& map(alloc: qp_array_y, qp_array_z, rec_br2_local_lift, rec_mrf, rec_qpoints, rec_rotvel, rec_uleft) &
!$omp& map(alloc: rec_uleft_dg, rec_uleftturb, rec_volume, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb, u_c_br2_aux_var) &
!$omp& map(alloc: u_c_val, u_c_valdg, u_ct_val) &
!$omp& map(alloc: rg_thetag)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, facex) &
!$omp& private(i, iconsidered, ii, iqp, l, mp_source1, mp_source2, mp_source3, ngp, norms, nx, ny) &
!$omp& private(nz, pointx)
#endif
      do ii = 1, nof_interior
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:gpu_max_dim) :: pox, poy, poz

            i = el_int(ii)

            weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
            rhs_val(i, :) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst_val(i, :) = zero

            iconsidered = i
            if (dg .eq. 1) then
               rhs_valdg(i, :, :) = zero

               call dg_vol_integral_accumulate_rhs(n, iconsidered, -1.00)

            end if

            mp_source3 = zero
            do l = 1, ielem_ifca(i) !for all their faces

               godflux2 = zero
               mp_source2 = zero
               nx = ielem_faceanglex(i, l)
               ny = ielem_faceangley(i, l)
               angle1 = nx
               angle2 = ny
               b_code = 0
               iqp = qp_line_n

               do ngp = 1, iqp !for all the gaussian quadrature points
                  pointx = ngp
                  facex = l
                  pointx = ngp

                  call get_states_interior(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, angle1, &
                     & angle2, nx, &
                     & ny, nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright)

                  call rotatef2d(n, cright_rot, cright, angle1, angle2)
                  !rotate wrt to normalvector of face and solve 1d riemann problem
                  call rotatef2d(n, cleft_rot, cleft, angle1, angle2)

               !realgas state
   !                                                   if (realgas.eq.1)then
   !                                                   call fix_conservative_state(cleft_rot)
   !                                                   call fix_conservative_state(cright_rot)
   !
   !                                                   end if

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                     call lmacht2d(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                  end if

                  select case (iriemann)

                  case (1) !hllc

                     call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & hllcflux(nof_variables + &
                           & 1:nof_variables + turbulenceequations + passivescalar)
                     end if

                  case (9) !hll

                     call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & hllcflux(nof_variables + &
                           & 1:nof_variables + turbulenceequations + passivescalar)
                     end if

                  case (2) !rusanov

                     call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & hllcflux(nof_variables + &
                           & 1:nof_variables + turbulenceequations + passivescalar)
                     end if
                  case (3) !roe

                     call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb2d(n, cright, cright_rot, angle1, angle2)
                     call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)

                     rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)

                  case (4) !roe

                     call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb2d(n, cright, cright_rot, angle1, angle2)
                     call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)

                     rhllcflux = hllcflux

                  end select

                  if (dg .eq. 1) then
                  ! riemann flux at interface quadrature points times basis

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

                  else !fv

                     godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))

                  end if

                  if (multispecies .eq. 1) then
                     mp_source2 = mp_source2 + mp_source1 * (weights_temp(ngp) * ielem_surf(i, l))
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (icoupleturb .eq. 0) then !first order upwind flux

                        norms = 0.50 * (cleft_rot(2) / cleft_rot(1) + cright_rot(2) / cright_rot(1))
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = 0.5 * ((norms * &
                           & (cturbl(1:turbulenceequations + passivescalar) + cturbr(1:turbulenceequations + passivescalar))) + &
                           & (abs(norms) * (cturbl(1:turbulenceequations + passivescalar) - (cturbr(1:turbulenceequations + &
                           & passivescalar)))))

   !
                     end if

                     godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar) + &
                        & (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))

   !

                  end if

               end do
               rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) + godflux2(1:nof_variables)
               if (multispecies .eq. 1) then
                  mp_source3 = mp_source3 + mp_source2
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) + &
                     & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

   !
               end if
   !                                     end if
            end do
            if (multispecies .eq. 1) then

               call subtract_volume_fraction_source(i, mp_source3, .false.)

            end if

            if (dg == 1) then

               if (multispecies .eq. 1) then

                  call subtract_volume_fraction_source(i, mp_source3, .true.)

               end if

            end if

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_convective2d_inner_cell

   subroutine calculate_fluxeshi_convective2d_bound_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i
      integer :: l, ngp, iqp, igoflux, b_code
      real :: norms
      integer :: iconsidered, facex, pointx
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(dg, icoupleturb, iriemann, lmach, multispecies, n, nof_bounded, nof_variables, passivescalar) &
!$omp& firstprivate(qp_line_n, turbulence, turbulenceequations, zero) &
!$omp& map(alloc: el_bnd, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, rhs_val, rhs_valdg, rhst_val) &
!$omp& map(alloc: weights_l) &
!$omp& private(angle1, angle2, b_code, facex) &
!$omp& private(i, iconsidered, igoflux, ii, iqp, l, mp_source1, mp_source2, mp_source3, ngp, norms) &
!$omp& private(nx, ny, nz, pointx) &
!$omp& map(alloc: bleed_end, bleed_plenum, bleed_porosity, bleed_start, bound_offset, boundhir, boundhir_dg) &
!$omp& map(alloc: catalytic_con, gamma_in, ibound_icode, ielem_bleedn, ielem_ibounds, ielem_idegfree, ielem_indexf) &
!$omp& map(alloc: ielem_ineigh, ielem_ineighb, ielem_ineighn, ielem_inter_id, ielem_iorder, ielem_itotalpoints) &
!$omp& map(alloc: ielem_nodes_faces, ielem_qface, ielem_totvolume, ielem_walldist, ielem_xxc, ielem_yyc, inoder4_cord) &
!$omp& map(alloc: integ_basis_dg_value, integ_basis_value, integ_basis_valuec, mp_a_in, mp_janaf, mp_m, mp_pinf) &
!$omp& map(alloc: mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in, qp_array_qp_weight, qp_array_x, qp_array_y) &
!$omp& map(alloc: qp_array_z, rec_br2_local_lift, rec_qpoints, rec_uleft, rec_uleft_dg, rec_uleftturb, rec_volume) &
!$omp& map(alloc: rg_hzero, rg_molm, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb, u_c_br2_aux_var, u_c_val, u_c_valdg) &
!$omp& map(alloc: u_ct_val)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, facex) &
!$omp& private(i, iconsidered, igoflux, ii, iqp, l, mp_source1, mp_source2, mp_source3, ngp, norms) &
!$omp& private(nx, ny, nz, pointx)
#endif
      do ii = 1, nof_bounded
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:gpu_max_dim) :: pox, poy, poz

            i = el_bnd(ii)

            weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)
            rhs_val(i, :) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst_val(i, :) = zero

            iconsidered = i
            mp_source3 = zero

            if (dg .eq. 1) then
               rhs_valdg(i, :, :) = zero

               call dg_vol_integral_accumulate_rhs(n, iconsidered, -1.00)

            end if

            do l = 1, ielem_ifca(i)
               facex = l
               igoflux = 0

               b_code = 0
               nx = ielem_faceanglex(i, l)
               ny = ielem_faceangley(i, l)
               angle1 = nx
               angle2 = ny

               iqp = qp_line_n
               godflux2 = zero
               mp_source2 = zero
               do ngp = 1, iqp
                  pointx = ngp
                  facex = l

                  call get_states_bounds2d(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, angle1, &
               & angle2, nx, ny, &
                     & nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright)

                  call rotatef2d(n, cright_rot, cright, angle1, angle2)
                  !rotate wrt to normalvector of face and solve 1d riemann problem
                  call rotatef2d(n, cleft_rot, cleft, angle1, angle2)
   !

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)

                     call lmacht2d(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
               & cturbl(1:turbulenceequations + &
                        & passivescalar)
                     cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
               & cturbr(1:turbulenceequations + &
                        & passivescalar)
                  end if

               !realgas state
   !                                                   if (realgas.eq.1)then
   !                                                   call fix_conservative_state(cleft_rot)
   !                                                   call fix_conservative_state(cright_rot)
   !
   !                                                   end if

                  if (b_code == 9 .and. duct_outlet_face_flux_enabled()) then
                     call pressure_outlet_face_flux(n, cright(1:nof_variables), cturbr(1:turbulenceequations + passivescalar), &
                        & [nx, ny, 0.0e0], rhllcflux, mp_source1)
                  else
                  select case (iriemann)

                  case (1) !hllc

                     call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & hllcflux(nof_variables + &
                           & 1:nof_variables + turbulenceequations + passivescalar)
                     end if

                  case (9) !hll

                     call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & hllcflux(nof_variables + &
                           & 1:nof_variables + turbulenceequations + passivescalar)
                     end if

                  case (2) !rusanov

                     call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & hllcflux(nof_variables + &
                           & 1:nof_variables + turbulenceequations + passivescalar)
                     end if

                  case (3) !roe

                     call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb2d(n, cright, cright_rot, angle1, angle2)
                     call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)

                     rhllcflux = hllcflux

                  case (4) !roe

                     if ((b_code .le. 0)) then

                        call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                        call rotateb2d(n, cright, cright_rot, angle1, angle2)

                        call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)
                        rhllcflux = hllcflux
                     else

                        call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                        call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                     end if

                  end select

                  end if

                  if (dg .eq. 1) then

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, 1.00)

                  else

                     godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))
                  end if

                  if (multispecies .eq. 1) then
                     mp_source2 = mp_source2 + mp_source1 * (weights_temp(ngp) * ielem_surf(i, l))
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (icoupleturb == 0 .and. .not. (b_code == 9 .and. duct_outlet_face_flux_enabled())) then

                        norms = 0.50 * (cleft_rot(2) / cleft_rot(1) + cright_rot(2) / cright_rot(1))
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = 0.5 * ((norms * &
                           & (cturbl(1:turbulenceequations + passivescalar) + cturbr(1:turbulenceequations + passivescalar))) + &
                           & (abs(norms) * (cturbl(1:turbulenceequations + passivescalar) - (cturbr(1:turbulenceequations + &
                           & passivescalar)))))

   !
                     end if

                     if ((b_code .eq. 4) .or. (b_code .eq. 3)) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
                     end if

                     godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar) + &
                        & (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))
   !

                  end if

               end do

               rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) + godflux2(1:nof_variables)
               if (multispecies .eq. 1) then
                  mp_source3 = mp_source3 + mp_source2
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) + &
                     & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

               end if
            end do
            if (multispecies .eq. 1) then
               call subtract_volume_fraction_source(i, mp_source3, .false.)

            end if

            if (dg == 1) then
               if (multispecies .eq. 1) then
                  call subtract_volume_fraction_source(i, mp_source3, .true.)

               end if

            end if

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_convective2d_bound_cell



   subroutine calculate_fluxeshi_diffusive(n)
      implicit none
      integer, intent(in) :: n

      call calculate_fluxeshi_diffusive_inner_cell(n)
      call calculate_fluxeshi_diffusive_bound_cell(n)
   end subroutine calculate_fluxeshi_diffusive

   subroutine calculate_fluxeshi_diffusive_inner_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, kc, b_code, srf, k, iv, nvt
      real :: wface, muturb
      integer :: iconsidered, facex, pointx, rg_i, rg_j, idxy
      real :: angle1, angle2, nx, ny, nz
      real :: mp_pinfl, mp_pinfr, gammal, gammar, rho12l, rho12r
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp
      real :: mp_mu_mix, mp_ktr_mix, mp_kve, mp_laml, mp_lamr, y_face
      real :: mp_lam_av, mp_ktr_mix_av, sumY_face, molar_sum, mbar

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(br2_yn, dg, dimensiona, lamx, lmach, mrf, multispecies, n, nof_interior, nof_species) &
!$omp& firstprivate(nof_variables, oo2, passivescalar, qp_quad_n, qp_triangle_n, realgas, rg_relax, sigma) &
!$omp& firstprivate(turbulence, turbulenceequations, turbulencemodel, zero) &
!$omp& map(alloc: el_int, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, ielem_types_faces) &
!$omp& map(alloc: ielem_walldist, rec_mrf, rg_molm, rhs_val, rhst_val, weights_q, weights_t) &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(iconsidered, idxy, ii, iqp, iv, k, kc, l, mbar) &
!$omp& private(molar_sum, mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix, mp_pinfl) &
!$omp& private(mp_pinfr, muturb, ngp, nvt, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j, rho12, rho12l, rho12r) &
!$omp& private(srf, sumy_face, u12, ux, uy, uz) &
!$omp& private(v12, vdamp, vx, vy, vz, w12, wface, wx, wy, wz, y_face) &
!$omp& map(alloc: gamma_in, ielem_idegfree, ielem_ineigh, ielem_ineighn, ielem_iorder, ielem_totvolume) &
!$omp& map(alloc: integ_basis_dg_value, integ_basis_value, integ_basis_valuec, mp_janaf, mp_m, mp_pinf, mp_thigh_in) &
!$omp& map(alloc: mp_tlow_in, mp_tmid_in, rec_qpoints, rec_rotvel, rec_uleft, rec_uleft_dg, rec_uleftturb) &
!$omp& map(alloc: rec_uleftturbv, rec_uleftv, rec_volume, rg_hzero, rg_thetag, rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k) &
!$omp& map(alloc: rgs_mg, rgs_sigmaa, rhs_valdg, u_ct_val)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(iconsidered, idxy, ii, iqp, iv, k, kc, l, mbar) &
!$omp& private(molar_sum, mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix, mp_pinfl) &
!$omp& private(mp_pinfr, muturb, ngp, nvt, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j, rho12, rho12l, rho12r) &
!$omp& private(srf, sumy_face, u12, ux, uy, uz) &
!$omp& private(v12, vdamp, vx, vy, vz, w12, wface, wx, wy, wz, y_face)
#endif
      do ii = 1, nof_interior
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:3) :: pox, poy, poz
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:4) :: viscl, laml
            real, dimension(1:2) :: turbmv
            real, dimension(1) :: etvm
            real, dimension(1:20) :: eddyfl, eddyfr
            real, dimension(1:gpu_max_nvar, 1:3) :: lcvgrad, rcvgrad
            real, dimension(gpu_max_extra_transport, 1:3) :: lcvgrad_t, rcvgrad_t
            real, dimension(1:gpu_max_nvar) :: fxv, fyv, fzv
            real, dimension(1:gpu_max_species) :: y_av, y_corr
            real, dimension(3, 3) :: taul, taur, tau
            real, dimension(3) :: q, nnn, nall, qvib
            real, dimension(1:gpu_max_species) :: rg_difl, rg_difr, rg_enthl, rg_enthr, rg_enthvbl, rg_enthvbr
            real, dimension(1:3) :: jl
            real, dimension(1:3) :: rg_sum_htr, rg_sum_hv, grady, rg_qv, rg_qtr
            real, dimension(1:gpu_max_species) :: rg_dif_av, rg_enth_av, rg_enthvb_av
            real, dimension(1:3) :: sumi, gradx, grad_a, sumJ, sumGradY
            real, dimension(1:3, 1:gpu_max_species) :: i_raw, gradY_all

            i = el_int(ii)
            iconsidered = i
            nvt = turbulenceequations + passivescalar

            damp = lamx
            if (br2_yn == 2) damp = 0.00
            if (mrf .eq. 1) then
               srf = rec_mrf(i)
            end if
            do l = 1, ielem_ifca(i) !for all their faces

               do iv = 1, nof_variables + nvt
                  godflux2(iv) = zero
               end do
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))
               nnn(1) = nx; nnn(2) = ny; nnn(3) = nz
               if (ielem_types_faces(i, l) .eq. 5) then
                  iqp = qp_quad_n
                  weights_temp(1:iqp) = weights_q(1:iqp)
               else
                  iqp = qp_triangle_n
                  weights_temp(1:iqp) = weights_t(1:iqp)
               end if

               do ngp = 1, iqp !for all the gaussian quadrature points

                  facex = l
                  pointx = ngp

                  call calculate_interior_viscous(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, &
               & angle1, angle2, &
                     & nx, ny, nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright, lcvgrad, rcvgrad, &
                        & lcvgrad_t, &
                     & rcvgrad_t)

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                     call rotatef(n, cright_rot, cright, angle1, angle2)
                     !rotate wrt to normalvector of face and solve 1d riemann problem
                     call rotatef(n, cleft_rot, cleft, angle1, angle2)
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                     call lmacht(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                     call rotateb(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb(n, cright, cright_rot, angle1, angle2)
                  end if

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)

                  viscl = 0.0e0; laml = 0.0e0
                  if (realgas == 1) then
                        call multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, rg_difl, rg_enthl, rg_enthvbl, gammal)

                        mp_laml = mp_kve
                        laml(1) = mp_ktr_mix
                        viscl(1) = mp_mu_mix

                        call multispecies_mixtures_rg(rightv, mp_mu_mix, mp_ktr_mix, mp_kve, rg_difr, rg_enthr, rg_enthvbr, gammar)

                        mp_lamr = mp_kve
                        laml(2) = mp_ktr_mix
                        viscl(2) = mp_mu_mix
                  else
                     call get_visc_conduct(n, leftv, rightv, viscl, laml)
                  end if

                  if (turbulence .eq. 1) then
                     if (turbulencemodel .eq. 1) then
                        turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                     if (turbulencemodel .eq. 2) then
                        eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                        eddyfl(4:6) = lcvgrad(1, 1:3); eddyfl(7:9) = lcvgrad(2, 1:3)
                        eddyfl(10:12) = lcvgrad(3, 1:3); eddyfl(13:15) = lcvgrad_t(1, 1:3)
                        eddyfl(16:18) = lcvgrad_t(2, 1:3)

                        eddyfr(1) = ielem_walldist(i); eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
                        eddyfr(4:6) = rcvgrad(1, 1:3); eddyfr(7:9) = rcvgrad(2, 1:3); eddyfr(10:12) = rcvgrad(3, 1:3)
                        eddyfr(13:15) = rcvgrad_t(1, 1:3); eddyfr(16:18) = rcvgrad_t(2, 1:3)
                        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                  end if

                  taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz &
               & = zero; wx = &
                     & zero; wy = zero; wz = zero;
                  do iv = 1, nof_variables
                     fxv(iv) = zero
                     fyv(iv) = zero
                     fzv(iv) = zero
                  end do
                  rho12 = zero;
                  u12 = zero; v12 = zero; w12 = zero

                  vdamp = (4.0 / 3.0) !*(( (viscl(1))+(viscl(2)))))
                  nall(1) = nx; nall(2) = ny; nall(3) = nz

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
                  call cons2div(n, leftv, mp_pinfl, gammal)
                  call cons2div(n, rightv, mp_pinfr, gammar)
                  rho12 = oo2 * (leftv(1) + rightv(1))
                  rho12l = leftv(1); rho12r = rightv(1)
                  u12 = oo2 * (leftv(2) + rightv(2))
                  v12 = oo2 * (leftv(3) + rightv(3))
                  w12 = oo2 * (leftv(4) + rightv(4))

                  if ((realgas .eq. 1) .and. (nof_species .gt. 0)) then
                     do rg_i = 1, nof_species
                        idxy = dimensiona + 3 + rg_i
                        y_av(rg_i) = 0.50 * (leftv(idxy) + rightv(idxy))
                     end do
                  END IF

                  do k = 1, nof_variables - 1
                     do iv = 1, 3
                        lcvgrad(k, iv) = ((lcvgrad(k, iv) + rcvgrad(k, iv)) / (2.00)) + damp * ((vdamp / abs(ielem_dih(i, &
                           & l))) * nall(iv) * (rightv(k + 1) - leftv(k + 1)))
                     end do
                  end do

               !now compute all the temperature gradients +real gas
                  if ((realgas .eq. 1) .or. (multispecies .eq. 1)) then

                     if (realgas .eq. 1) then

                     !diffusion coefficient for species rg_diffl(1:nof_species),rg_diffr(1:nof_species)
                     !viscosity for the mixture viscl(1)-left,viscl(2)-right
                     !thermal conductivity mix laml(1)-left,laml(2)-right
                     !vibrational thermal conductivity!mplaml-left,mp_lamr-right
                     !species diffusion coefficients !rg_difl(1:nof_species),rg_difr(1:nof_species)
                     !enthalpies for species        !rg_enthl(1:nof_species),rg_enthr(1:nof_species)
                     !vibrational enthalpies for species !!rg_enthvbl(1:nof_species),rg_enthvbr(1:nof_species)

                        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)


                        mp_ktr_mix_av = 0.50 * (laml(2) + laml(1)) ! translational conductivity
                        if (turbulence == 1) mp_ktr_mix_av = 0.5e0 * (laml(3) + laml(4))

                        mp_lam_av = 0.50 * (mp_lamr + mp_laml) ! ! vibrational conductivity

                        rg_dif_av = 0.50 * (rg_difr + rg_difl)

                        rg_enth_av = 0.50 * (rg_enthr + rg_enthl)
                        rg_enthvb_av = 0.50 * (rg_enthvbr + rg_enthvbl)

                        if (turbulence .eq. 1) then
                           q(1:3) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:3))
                        end if

                        rg_sum_htr(1:dimensiona) = 0.00
                        rg_sum_hv(1:dimensiona) = 0.00

                        if (rg_relax .ge. 1) then

                           sumi(1:dimensiona) = 0.00
                           sumJ(1:dimensiona) = 0.00
                           i_raw(1:dimensiona, 1:nof_species) = 0.00

                        ! normalise face mass fractions used by correction velocity
                           sumY_face = 0.00
                           do rg_j = 1, nof_species
                              y_corr(rg_j) = max(y_av(rg_j), 0.00)
                              sumY_face = sumY_face + y_corr(rg_j)
                           end do

                           if (sumY_face .gt. 1.0e-30) then
                              do rg_j = 1, nof_species
                                 y_corr(rg_j) = y_corr(rg_j) / sumY_face
                              end do
                           else
                              do rg_j = 1, nof_species
                                 y_corr(rg_j) = 1.00 / real(nof_species)
                              end do
                           end if

                        ! mixture molar sum from corrected face mass fractions
                           molar_sum = 0.00
                           do rg_j = 1, nof_species
                              molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
                           end do
                           molar_sum = max(molar_sum, 1.0e-30)
                           mbar = 1.00 / molar_sum

                        ! project species gradients so sum_k grad(Y_k)=0
                           sumGradY(1:dimensiona) = 0.00
                           do rg_j = 1, nof_species
                              idxy = dimensiona + 2 + rg_j
                              gradY_all(1:dimensiona, rg_j) = lcvgrad(idxy, 1:dimensiona)
                              sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona, rg_j)
                           end do
                           do rg_j = 1, nof_species
                              gradY_all(1:dimensiona, rg_j) = gradY_all(1:dimensiona, rg_j) - &
                                 & y_corr(rg_j) * sumGradY(1:dimensiona)
                           end do

                           grad_a(1:dimensiona) = 0.00
                           do rg_j = 1, nof_species
                              grady(1:dimensiona) = gradY_all(1:dimensiona, rg_j)
                              grad_a(1:dimensiona) = grad_a(1:dimensiona) + grady(1:dimensiona) / rg_molm(rg_j)
                           end do

                        ! raw mixture-averaged flux from grad(X_i)
                           do rg_i = 1, nof_species
                              grady(1:dimensiona) = gradY_all(1:dimensiona, rg_i)
                              gradx(1:dimensiona) = &
                                 & (grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
                                 & - (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona)) &
                                 & / (molar_sum * molar_sum)

                              i_raw(1:dimensiona, rg_i) = &
                                 & -rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)
                              sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona, rg_i)
                           end do

                        ! corrected species flux, conservative because sum(y_corr)=1
                           do rg_i = 1, nof_species
                              y_face = y_corr(rg_i)
                              jl(1:dimensiona) = i_raw(1:dimensiona, rg_i) - y_face * sumi(1:dimensiona)
                              sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

                              fxv(dimensiona + 3 + rg_i) = fxv(dimensiona + 3 + rg_i) - jl(1)
                              fyv(dimensiona + 3 + rg_i) = fyv(dimensiona + 3 + rg_i) - jl(2)
                              fzv(dimensiona + 3 + rg_i) = fzv(dimensiona + 3 + rg_i) - jl(3)

                              rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
                                 & + rg_enth_av(rg_i) * jl(1:dimensiona)
                              rg_sum_hv(1:dimensiona) = rg_sum_hv(1:dimensiona) &
                                 & + rg_enthvb_av(rg_i) * jl(1:dimensiona)
                           end do

                        end if

                     ! =============================================================
                     ! 4. conductive heat fluxes (fourier)
                     ! =============================================================
                        rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona + 1, 1:dimensiona)
                        rg_qv(1:dimensiona) = -mp_lam_av * lcvgrad(dimensiona + 2, 1:dimensiona)

                     ! =============================================================
                     ! 5. total diffusive fluxes
                     ! =============================================================
                        q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
                           & + rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

                        qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

                     ! ------------------------------------------------
                     ! 6. add to conservative flux vectors
                     ! ------------------------------------------------
                     ! total energy equation (ρe)
                        fxv(5) = fxv(5) - q(1)
                        fyv(5) = fyv(5) - q(2)
                        fzv(5) = fzv(5) - q(3)

                     ! vibrational energy equation (ρev)
                        fxv(6) = fxv(6) - qvib(1)
                        fyv(6) = fyv(6) - qvib(2)
                        fzv(6) = fzv(6) - qvib(3)

                     end if !real gas ends

                     if (multispecies .eq. 1) then
                     !viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

                        if (turbulence .eq. 1) then
                           q(1:3) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:3))
                        else
                           q(1:3) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:3))
                        end if

                        fxv(5) = fxv(5) - q(1); fyv(5) = fyv(5) - q(2); fzv(5) = fzv(5) - q(3)
                     end if

                  end if

                  if (realgas == 0 .and. multispecies == 0) then
                     if (turbulence .eq. 1) then
                        q(1:3) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:3))
                     else
                        q(1:3) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:3))
                     end if

                     fxv(5) = fxv(5) - q(1); fyv(5) = fyv(5) - q(2); fzv(5) = fzv(5) - q(3)

                  end if

               !left state derivatives

               ! determine taul!!
                  ux = lcvgrad(1, 1); uy = lcvgrad(1, 2); uz = lcvgrad(1, 3);
                  vx = lcvgrad(2, 1); vy = lcvgrad(2, 2); vz = lcvgrad(2, 3);
                  wx = lcvgrad(3, 1); wy = lcvgrad(3, 2); wz = lcvgrad(3, 3);
               ! tau_xx
                  taul(1, 1) = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
               ! tau_yy
                  taul(2, 2) = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
               ! tau_zz
                  taul(3, 3) = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy

               ! tau_xy
                  taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

               ! tau_xz
                  taul(1, 3) = (wx + uz); taul(3, 1) = taul(1, 3)

               ! tau_yz
                  taul(2, 3) = (vz + wy); taul(3, 2) = taul(2, 3)
               !end determine taul

               ! average and multiplay by viscosity
                  if (turbulence .eq. 1) then
                     tau = oo2 * (((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))) * taul
                  else
                     tau = oo2 * (((viscl(1)) + (viscl(2)))) * taul
                  end if

               ! now addition into momentum fluxes
                  do kc = 2, 4
                     fxv(kc) = fxv(kc) + tau(1, kc - 1)
                     fyv(kc) = fyv(kc) + tau(2, kc - 1)
                     fzv(kc) = fzv(kc) + tau(3, kc - 1)
                  end do

               ! compute interface velocities

                  fxv(5) = fxv(5) + u12 * tau(1, 1) + v12 * tau(1, 2) + w12 * tau(1, 3)
                  fyv(5) = fyv(5) + u12 * tau(2, 1) + v12 * tau(2, 2) + w12 * tau(2, 3)
                  fzv(5) = fzv(5) + u12 * tau(3, 1) + v12 * tau(3, 2) + w12 * tau(3, 3)

                  do iv = 1, nof_variables
                     hllcflux(iv) = nx * fxv(iv) + ny * fyv(iv) + nz * fzv(iv)
                  end do

                  if (dg .eq. 1) then
                     do iv = 1, nof_variables
                        rhllcflux(iv) = hllcflux(iv)
                     end do

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, -1.00)

                  else

                     wface = weights_temp(ngp) * ielem_surf(i, l)
                     do iv = 1, nof_variables
                        godflux2(iv) = godflux2(iv) + (hllcflux(iv) * wface)
                     end do

                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (turbulence .eq. 1) then
                        muturb = (oo2 * (viscl(1) + viscl(2))) + (oo2 * (viscl(3) + viscl(4)))
                     else
                        muturb = oo2 * (viscl(1) + viscl(2))
                     end if
                     do iv = 1, nvt
                        hllcflux(nof_variables + iv) = muturb * oo2 * ((lcvgrad_t(iv, 1) + rcvgrad_t(iv, 1)) * nx + &
                           & (lcvgrad_t(iv, 2) + rcvgrad_t(iv, 2)) * ny + (lcvgrad_t(iv, 3) + rcvgrad_t(iv, 3)) * nz)
                     end do
                     if (turbulencemodel .eq. 1) then
                        do iv = 1, nvt
                           hllcflux(nof_variables + iv) = hllcflux(nof_variables + iv) / sigma
                        end do
                     end if
                     wface = weights_temp(ngp) * ielem_surf(i, l)
                     do iv = 1, nvt
                        godflux2(nof_variables + iv) = godflux2(nof_variables + iv) + (hllcflux(nof_variables + iv) * wface)
                     end do
                  end if

               end do

               do iv = 1, nof_variables
                  rhs_val(i, iv) = rhs_val(i, iv) - godflux2(iv)
               end do

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  do iv = 1, nvt
                     rhst_val(i, iv) = rhst_val(i, iv) - godflux2(nof_variables + iv)
                  end do

               end if

            end do

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_diffusive_inner_cell

   subroutine calculate_fluxeshi_diffusive_bound_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, kc, igoflux, b_code, srf, k, iv, nvt
      real :: wface, muturb
      integer :: iconsidered, facex, pointx, rg_i, rg_j, idxy, icompute
      real :: angle1, angle2, nx, ny, nz
      real :: mp_pinfl, mp_pinfr, gammal, gammar, rho12l, rho12r
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp
      real :: mp_mu_mix, mp_ktr_mix, mp_kve, mp_laml, mp_lamr, y_face
      real :: mp_lam_av, mp_ktr_mix_av, sumY_face, molar_sum, mbar

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(br2_yn, catalytic_wall, dg, dimensiona, lamx, lmach, mrf, multispecies, n, nof_bounded) &
!$omp& firstprivate(nof_species, nof_variables, oo2, passivescalar, qp_quad_n, qp_triangle_n, realgas, rg_relax) &
!$omp& firstprivate(sigma, turbulence, turbulenceequations, turbulencemodel, zero) &
!$omp& map(alloc: el_bnd, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, ielem_types_faces) &
!$omp& map(alloc: ielem_walldist, rec_mrf, rg_molm, rhs_val, rhst_val, weights_q, weights_t) &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(icompute, iconsidered, idxy, igoflux, ii, iqp, iv, k, kc, l) &
!$omp& private(mbar, molar_sum, mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix) &
!$omp& private(mp_pinfl, mp_pinfr, muturb, ngp, nvt, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j, rho12, rho12l, rho12r) &
!$omp& private(srf, sumy_face, u12, ux, uy, uz) &
!$omp& private(v12, vdamp, vx, vy, vz, w12, wface, wx, wy, wz, y_face) &
!$omp& map(alloc: bound_offset, boundhir, boundhir_dg, catalytic_con, gamma_in, ibound_icode, ielem_ibounds) &
!$omp& map(alloc: ielem_idegfree, ielem_indexf, ielem_ineigh, ielem_ineighb, ielem_ineighn, ielem_inter_id) &
!$omp& map(alloc: ielem_iorder, ielem_nodes_faces, ielem_qface, ielem_totvolume, inoder4_cord, integ_basis_dg_value) &
!$omp& map(alloc: integ_basis_value, integ_basis_valuec, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in) &
!$omp& map(alloc: mp_tlow_in, mp_tmid_in, rec_qpoints, rec_rotvel, rec_uleft, rec_uleft_dg, rec_uleftturb) &
!$omp& map(alloc: rec_uleftturbv, rec_uleftv, rec_volume, rg_hzero, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb) &
!$omp& map(alloc: rgs_eps_over_k, rgs_mg, rgs_sigmaa, rhs_valdg, u_ct_val)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(icompute, iconsidered, idxy, igoflux, ii, iqp, iv, k, kc, l) &
!$omp& private(mbar, molar_sum, mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix) &
!$omp& private(mp_pinfl, mp_pinfr, muturb, ngp, nvt, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j, rho12, rho12l, rho12r) &
!$omp& private(srf, sumy_face, u12, ux, uy, uz) &
!$omp& private(v12, vdamp, vx, vy, vz, w12, wface, wx, wy, wz, y_face)
#endif
      do ii = 1, nof_bounded
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:3) :: pox, poy, poz
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:4) :: viscl, laml
            real, dimension(1:2) :: turbmv
            real, dimension(1) :: etvm
            real, dimension(1:20) :: eddyfl, eddyfr
            real, dimension(1:gpu_max_nvar, 1:3) :: lcvgrad, rcvgrad
            real, dimension(gpu_max_extra_transport, 1:3) :: lcvgrad_t, rcvgrad_t
            real, dimension(1:gpu_max_nvar) :: fxv, fyv, fzv
            real, dimension(1:gpu_max_species) :: y_av, y_corr
            real, dimension(3, 3) :: taul, taur, tau
            real, dimension(3) :: q, nnn, nall, qvib
            real, dimension(1:gpu_max_species) :: rg_difl, rg_difr, rg_enthl, rg_enthr, rg_enthvbl, rg_enthvbr
            real, dimension(1:3) :: jl
            real, dimension(1:3) :: rg_sum_htr, rg_sum_hv, grady, rg_qv, rg_qtr
            real, dimension(1:gpu_max_species) :: rg_dif_av, rg_enth_av, rg_enthvb_av
            real, dimension(1:3) :: sumi, gradx, grad_a, sumJ, sumGradY
            real, dimension(1:3, 1:gpu_max_species) :: i_raw, gradY_all

            i = el_bnd(ii)
            iconsidered = i
            nvt = turbulenceequations + passivescalar
            if (mrf .eq. 1) then
               srf = rec_mrf(i)
            end if

            do l = 1, ielem_ifca(i)

               igoflux = 0

               damp = lamx
               if (br2_yn == 2) damp = 0.00
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = (cos(angle1) * sin(angle2))
               ny = (sin(angle1) * sin(angle2))
               nz = (cos(angle2))
               nnn(1) = nx; nnn(2) = ny; nnn(3) = nz
               if (ielem_types_faces(i, l) .eq. 5) then
                  iqp = qp_quad_n
                  weights_temp(1:iqp) = weights_q(1:iqp)
               else
                  iqp = qp_triangle_n
                  weights_temp(1:iqp) = weights_t(1:iqp)
               end if

               do iv = 1, nof_variables + nvt
                  godflux2(iv) = zero
               end do

               b_code = 0

               do ngp = 1, iqp
                  facex = l
                  pointx = ngp

                  call calculate_bounded_viscous(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, angle1, &
                     & angle2, &
                     & nx, ny, nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright, lcvgrad, rcvgrad, &
                        & lcvgrad_t, &
                     & rcvgrad_t)

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                     call rotatef(n, cright_rot, cright, angle1, angle2)
                     !rotate wrt to normalvector of face and solve 1d riemann problem
                     call rotatef(n, cleft_rot, cleft, angle1, angle2)
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                     call lmacht(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                     call rotateb(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb(n, cright, cright_rot, angle1, angle2)
                  end if

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)

                  viscl = 0.0e0; laml = 0.0e0
                  if (realgas == 1) then
                        call multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, rg_difl, rg_enthl, rg_enthvbl, gammal)

                        mp_laml = mp_kve
                        laml(1) = mp_ktr_mix
                        viscl(1) = mp_mu_mix

                        call multispecies_mixtures_rg(rightv, mp_mu_mix, mp_ktr_mix, mp_kve, rg_difr, rg_enthr, rg_enthvbr, gammar)

                        mp_lamr = mp_kve
                        laml(2) = mp_ktr_mix
                        viscl(2) = mp_mu_mix
                  else
                     call get_visc_conduct(n, leftv, rightv, viscl, laml)
                  end if

                  if (turbulence .eq. 1) then
                     if (turbulencemodel .eq. 1) then
                        turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                     if (turbulencemodel .eq. 2) then
                        eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                        eddyfl(4:6) = lcvgrad(1, 1:3); eddyfl(7:9) = lcvgrad(2, 1:3)
                        eddyfl(10:12) = lcvgrad(3, 1:3); eddyfl(13:15) = lcvgrad_t(1, 1:3)
                        eddyfl(16:18) = lcvgrad_t(2, 1:3)

                        eddyfr(1) = ielem_walldist(i); eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
                        eddyfr(4:6) = rcvgrad(1, 1:3); eddyfr(7:9) = rcvgrad(2, 1:3); eddyfr(10:12) = rcvgrad(3, 1:3)
                        eddyfr(13:15) = rcvgrad_t(1, 1:3); eddyfr(16:18) = rcvgrad_t(2, 1:3)
                        call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                  end if

                  taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz &
               & = zero; wx = &
                     & zero; wy = zero; wz = zero;
                  do iv = 1, nof_variables
                     fxv(iv) = zero
                     fyv(iv) = zero
                     fzv(iv) = zero
                  end do
                  rho12 = zero;
                  u12 = zero; v12 = zero; w12 = zero

                  if ((b_code .gt. 0)) then
                     damp = zero
                  end if

                  vdamp = (4.0 / 3.0) !*(( (viscl(1))+(viscl(2)))))
                  nall(1) = nx; nall(2) = ny; nall(3) = nz

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
                  call cons2div(n, leftv, mp_pinfl, gammal)
                  call cons2div(n, rightv, mp_pinfr, gammar)
                  rho12 = oo2 * (leftv(1) + rightv(1))
                  rho12l = leftv(1); rho12r = rightv(1)
                  u12 = oo2 * (leftv(2) + rightv(2))
                  v12 = oo2 * (leftv(3) + rightv(3))
                  w12 = oo2 * (leftv(4) + rightv(4))

                  if ((realgas .eq. 1) .and. (nof_species .gt. 0)) then
                     do rg_i = 1, nof_species
                        idxy = dimensiona + 3 + rg_i
                        y_av(rg_i) = 0.50 * (leftv(idxy) + rightv(idxy))
                     end do
                  end if

                  do k = 1, nof_variables - 1
                     do iv = 1, 3
                        lcvgrad(k, iv) = ((lcvgrad(k, iv) + rcvgrad(k, iv)) / (2.00)) + damp * ((vdamp / abs(ielem_dih(i, &
                           & l))) * nall(iv) * (rightv(k + 1) - leftv(k + 1)))
                     end do
                  end do

                  icompute = 0

               !now compute all the temperature gradients +real gas
                  if ((realgas .eq. 1) .or. (multispecies .eq. 1)) then

                     if (realgas .eq. 1) then

                     !diffusion coefficient for species rg_diffl(1:nof_species),rg_diffr(1:nof_species)
                     !viscosity for the mixture viscl(1)-left,viscl(2)-right
                     !thermal conductivity mix laml(1)-left,laml(2)-right
                     !vibrational thermal conductivity!mplaml-left,mp_lamr-right
                     !species diffusion coefficients !rg_difl(1:nof_species),rg_difr(1:nof_species)
                     !enthalpies for species        !rg_enthl(1:nof_species),rg_enthr(1:nof_species)
                     !vibrational enthalpies for species !!rg_enthvbl(1:nof_species),rg_enthvbr(1:nof_species)

                        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)


                        mp_ktr_mix_av = 0.50 * (laml(1) + laml(2))
                        if (turbulence == 1) mp_ktr_mix_av = 0.5e0 * (laml(3) + laml(4))
                        mp_lam_av = 0.50 * (mp_laml + mp_lamr)

                     ! mixture-averaged species diffusion coefficients
                        rg_dif_av(:) = 0.50 * (rg_difl(:) + rg_difr(:))

                     ! face-averaged species enthalpies
                        rg_enth_av(:) = 0.50 * (rg_enthl(:) + rg_enthr(:))
                        rg_enthvb_av(:) = 0.50 * (rg_enthvbl(:) + rg_enthvbr(:))

                        if (turbulence .eq. 1) then
                           q(1:3) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:3))
                        end if

                        rg_sum_htr(1:dimensiona) = 0.00
                        rg_sum_hv(1:dimensiona) = 0.00

                     ! =============================================================
                     ! 3. mixture-averaged diffusion fluxes (no pressure terms)
                     ! =============================================================
                        rg_sum_htr(1:dimensiona) = 0.00
                        rg_sum_hv(1:dimensiona) = 0.00

                        sumi(1:dimensiona) = 0.00 ! σ i_k

                        if (rg_relax .ge. 1) then
                           if ((b_code .eq. 4) .and. (catalytic_wall .eq. 0)) then
                              icompute = 1
                           end if

                           if (icompute .eq. 0) then

                              sumi(1:dimensiona) = 0.00
                              sumJ(1:dimensiona) = 0.00
                              i_raw(1:dimensiona, 1:nof_species) = 0.00

                           ! normalise face mass fractions used by correction velocity
                              sumY_face = 0.00
                              do rg_j = 1, nof_species
                                 y_corr(rg_j) = max(y_av(rg_j), 0.00)
                                 sumY_face = sumY_face + y_corr(rg_j)
                              end do

                              if (sumY_face .gt. 1.0e-30) then
                                 do rg_j = 1, nof_species
                                    y_corr(rg_j) = y_corr(rg_j) / sumY_face
                                 end do
                              else
                                 do rg_j = 1, nof_species
                                    y_corr(rg_j) = 1.00 / real(nof_species)
                                 end do
                              end if

                           ! mixture molar sum from corrected face mass fractions
                              molar_sum = 0.00
                              do rg_j = 1, nof_species
                                 molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
                              end do
                              molar_sum = max(molar_sum, 1.0e-30)
                              mbar = 1.00 / molar_sum

                           ! project species gradients so sum_k grad(Y_k)=0
                              sumGradY(1:dimensiona) = 0.00
                              do rg_j = 1, nof_species
                                 idxy = dimensiona + 2 + rg_j
                                 gradY_all(1:dimensiona, rg_j) = lcvgrad(idxy, 1:dimensiona)
                                 sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona, rg_j)
                              end do
                              do rg_j = 1, nof_species
                                 gradY_all(1:dimensiona, rg_j) = gradY_all(1:dimensiona, rg_j) - &
                                    & y_corr(rg_j) * sumGradY(1:dimensiona)
                              end do

                              grad_a(1:dimensiona) = 0.00
                              do rg_j = 1, nof_species
                                 grady(1:dimensiona) = gradY_all(1:dimensiona, rg_j)
                                 grad_a(1:dimensiona) = grad_a(1:dimensiona) + grady(1:dimensiona) / rg_molm(rg_j)
                              end do

                           ! raw mixture-averaged flux from grad(X_i)
                              do rg_i = 1, nof_species
                                 grady(1:dimensiona) = gradY_all(1:dimensiona, rg_i)
                                 gradx(1:dimensiona) = &
                                    & (grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
                                    & - (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona)) &
                                    & / (molar_sum * molar_sum)

                                 i_raw(1:dimensiona, rg_i) = &
                                    & -rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)
                                 sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona, rg_i)
                              end do

                           ! corrected species flux, conservative because sum(y_corr)=1
                              do rg_i = 1, nof_species
                                 y_face = y_corr(rg_i)
                                 jl(1:dimensiona) = i_raw(1:dimensiona, rg_i) - y_face * sumi(1:dimensiona)
                                 sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

                                 fxv(dimensiona + 3 + rg_i) = fxv(dimensiona + 3 + rg_i) - jl(1)
                                 fyv(dimensiona + 3 + rg_i) = fyv(dimensiona + 3 + rg_i) - jl(2)
                                 fzv(dimensiona + 3 + rg_i) = fzv(dimensiona + 3 + rg_i) - jl(3)

                                 rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
                                    & + rg_enth_av(rg_i) * jl(1:dimensiona)
                                 rg_sum_hv(1:dimensiona) = rg_sum_hv(1:dimensiona) &
                                    & + rg_enthvb_av(rg_i) * jl(1:dimensiona)
                              end do

                           end if
                        end if

                     ! =============================================================
                     ! 4. conductive heat fluxes (fourier)
                     ! =============================================================
                        rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona + 1, 1:dimensiona)
                        rg_qv(1:dimensiona) = -mp_lam_av * lcvgrad(dimensiona + 2, 1:dimensiona)

                     ! =============================================================
                     ! 5. total diffusive fluxes
                     ! =============================================================
                        q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
                           & + rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

                        qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

                     ! ------------------------------------------------
                     ! 6. add to conservative flux vectors
                     ! ------------------------------------------------
                     ! total energy equation (ρe)
                        fxv(5) = fxv(5) - q(1)
                        fyv(5) = fyv(5) - q(2)
                        fzv(5) = fzv(5) - q(3)

                     ! vibrational energy equation (ρev)
                        fxv(6) = fxv(6) - qvib(1)
                        fyv(6) = fyv(6) - qvib(2)
                        fzv(6) = fzv(6) - qvib(3)

                     end if !real gas ends

                     if (multispecies .eq. 1) then
                     !viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

                        if (turbulence .eq. 1) then
                           q(1:3) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:3))
                        else
                           q(1:3) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:3))
                        end if

                        fxv(5) = fxv(5) - q(1); fyv(5) = fyv(5) - q(2); fzv(5) = fzv(5) - q(3)
                     end if

                  end if
                  if (realgas == 0 .and. multispecies == 0) then
                     if (turbulence .eq. 1) then
                        q(1:3) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:3))
                     else
                        q(1:3) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:3))
                     end if

                     fxv(5) = fxv(5) - q(1); fyv(5) = fyv(5) - q(2); fzv(5) = fzv(5) - q(3)

                  end if

               !left state derivatives

               ! determine taul!!
                  ux = lcvgrad(1, 1); uy = lcvgrad(1, 2); uz = lcvgrad(1, 3);
                  vx = lcvgrad(2, 1); vy = lcvgrad(2, 2); vz = lcvgrad(2, 3);
                  wx = lcvgrad(3, 1); wy = lcvgrad(3, 2); wz = lcvgrad(3, 3);
               ! tau_xx
                  taul(1, 1) = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy - (2.00 / 3.00) * wz
               ! tau_yy
                  taul(2, 2) = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux - (2.00 / 3.00) * wz
               ! tau_zz
                  taul(3, 3) = (4.00 / 3.00) * wz - (2.00 / 3.00) * ux - (2.00 / 3.00) * vy

               ! tau_xy
                  taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

               ! tau_xz
                  taul(1, 3) = (wx + uz); taul(3, 1) = taul(1, 3)

               ! tau_yz
                  taul(2, 3) = (vz + wy); taul(3, 2) = taul(2, 3)
               !end determine taul

               ! average and multiplay by viscosity
                  if (turbulence .eq. 1) then
                     tau = oo2 * (((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))) * taul
                  else
                     tau = oo2 * (((viscl(1)) + (viscl(2)))) * taul
                  end if

               ! now addition into momentum fluxes
                  do kc = 2, 4
                     fxv(kc) = fxv(kc) + tau(1, kc - 1)
                     fyv(kc) = fyv(kc) + tau(2, kc - 1)
                     fzv(kc) = fzv(kc) + tau(3, kc - 1)
                  end do

                  fxv(5) = fxv(5) + u12 * tau(1, 1) + v12 * tau(1, 2) + w12 * tau(1, 3)
                  fyv(5) = fyv(5) + u12 * tau(2, 1) + v12 * tau(2, 2) + w12 * tau(2, 3)
                  fzv(5) = fzv(5) + u12 * tau(3, 1) + v12 * tau(3, 2) + w12 * tau(3, 3)

                  do iv = 1, nof_variables
                     hllcflux(iv) = nx * fxv(iv) + ny * fyv(iv) + nz * fzv(iv)
                  end do

                  if (realgas .eq. 1) then
                     if ((b_code .eq. 4) .and. (catalytic_wall .eq. 0)) then
                        do iv = 7, nof_variables
                           hllcflux(iv) = zero
                        end do
                     end if
                  end if

                  if (dg .eq. 1) then
                     do iv = 1, nof_variables
                        rhllcflux(iv) = hllcflux(iv)
                     end do

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, -1.00)

                  else
                     wface = weights_temp(ngp) * ielem_surf(i, l)
                     do iv = 1, nof_variables
                        godflux2(iv) = godflux2(iv) + (hllcflux(iv) * wface)
                     end do
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (turbulence .eq. 1) then
                        muturb = (oo2 * (viscl(1) + viscl(2))) + (oo2 * (viscl(3) + viscl(4)))
                     else
                        muturb = oo2 * (viscl(1) + viscl(2))
                     end if
                     do iv = 1, nvt
                        hllcflux(nof_variables + iv) = muturb * oo2 * ((lcvgrad_t(iv, 1) + rcvgrad_t(iv, 1)) * nx + &
                           & (lcvgrad_t(iv, 2) + rcvgrad_t(iv, 2)) * ny + (lcvgrad_t(iv, 3) + rcvgrad_t(iv, 3)) * nz)
                     end do
                     if (turbulencemodel .eq. 1) then
                        do iv = 1, nvt
                           hllcflux(nof_variables + iv) = hllcflux(nof_variables + iv) / sigma
                        end do
                     end if

                     if ((b_code .eq. 3)) then
                        do iv = 1, nvt
                           hllcflux(nof_variables + iv) = zero
                        end do

                     end if

                     wface = weights_temp(ngp) * ielem_surf(i, l)
                     do iv = 1, nvt
                        godflux2(nof_variables + iv) = godflux2(nof_variables + iv) + (hllcflux(nof_variables + iv) * wface)
                     end do
                  end if

               end do

               do iv = 1, nof_variables
                  rhs_val(i, iv) = rhs_val(i, iv) - godflux2(iv)
               end do

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  do iv = 1, nvt
                     rhst_val(i, iv) = rhst_val(i, iv) - godflux2(nof_variables + iv)
                  end do

               end if
            end do

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_diffusive_bound_cell





   subroutine calculate_fluxeshi_diffusive2d(n)
      implicit none
      integer, intent(in) :: n

      call calculate_fluxeshi_diffusive2d_inner_cell(n)
      call calculate_fluxeshi_diffusive2d_bound_cell(n)
   end subroutine calculate_fluxeshi_diffusive2d

   subroutine calculate_fluxeshi_diffusive2d_inner_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, kc, b_code
      integer :: iconsidered, facex, pointx, k, rg_i, rg_j, idxy
      real :: angle1, angle2, nx, ny, nz
      real :: mp_pinfl, mp_pinfr, gammal, gammar, rho12l, rho12r
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp, y_face
      real :: mp_mu_mix, mp_ktr_mix, mp_kve, mp_laml, mp_lamr
      real :: mp_lam_av, mp_ktr_mix_av
      real :: sumY_face
      real :: molar_sum, mbar

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(br2_yn, dg, dimensiona, lamx, lmach, multispecies, n, nof_interior, nof_species, nof_variables) &
!$omp& firstprivate(oo2, passivescalar, qp_line_n, realgas, rg_relax, sigma, turbulence, turbulenceequations) &
!$omp& firstprivate(turbulencemodel, zero) &
!$omp& map(alloc: el_int, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, ielem_walldist) &
!$omp& map(alloc: rg_molm, rhs_val, rhst_val, weights_l) &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(iconsidered, idxy, ii, iqp, k, kc, l, mbar, molar_sum) &
!$omp& private(mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix, mp_pinfl, mp_pinfr) &
!$omp& private(ngp, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j) &
!$omp& private(rho12, rho12l, rho12r) &
!$omp& private(sumy_face, u12, ux, uy, uz, v12, vdamp, vx, vy, vz, w12) &
!$omp& private(wx, wy, wz, y_face) &
!$omp& map(alloc: gamma_in, ielem_idegfree, ielem_ineigh, ielem_ineighn, ielem_iorder, ielem_totvolume) &
!$omp& map(alloc: integ_basis_dg_value, integ_basis_value, integ_basis_valuec, mp_janaf, mp_m, mp_pinf, mp_thigh_in) &
!$omp& map(alloc: mp_tlow_in, mp_tmid_in, rec_mrf, rec_qpoints, rec_rotvel, rec_uleft, rec_uleft_dg, rec_uleftturb) &
!$omp& map(alloc: rec_uleftturbv, rec_uleftv, rec_volume, rg_hzero, rg_thetag, rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k) &
!$omp& map(alloc: rgs_mg, rgs_sigmaa, rhs_valdg, u_ct_val)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(iconsidered, idxy, ii, iqp, k, kc, l, mbar, molar_sum) &
!$omp& private(mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix, mp_pinfl, mp_pinfr) &
!$omp& private(ngp, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j) &
!$omp& private(rho12, rho12l, rho12r) &
!$omp& private(sumy_face, u12, ux, uy, uz, v12, vdamp, vx, vy, vz, w12) &
!$omp& private(wx, wy, wz, y_face)
#endif
      do ii = 1, nof_interior
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:gpu_max_dim) :: pox, poy, poz
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:4) :: viscl, laml
            real, dimension(1:gpu_max_species) :: y_av, x_av
            real, dimension(1:2) :: turbmv
            real, dimension(1) :: etvm
            real, dimension(1:20) :: eddyfl, eddyfr
            real, dimension(1:gpu_max_nvar, 1:gpu_max_dim) :: lcvgrad, rcvgrad
            real, dimension(gpu_max_extra_transport, 1:gpu_max_dim) :: lcvgrad_t, rcvgrad_t
            real, dimension(1:gpu_max_nvar) :: fxv, fyv, fzv
            real, dimension(2, 2) :: taul, taur, tau
            real, dimension(2) :: q, nall, qvib
            real, dimension(1:gpu_max_species) :: rg_difl, rg_difr, rg_enthl, rg_enthr, rg_enthvbl, rg_enthvbr
            real, dimension(1:gpu_max_dim) :: jl
            real, dimension(1:gpu_max_dim) :: rg_sum_htr, rg_sum_hv, grady, rg_qv, rg_qtr
            real, dimension(1:gpu_max_species) :: rg_dif_av, rg_enth_av, rg_enthvb_av
            real, dimension(1:gpu_max_dim) :: sumi, gradx, grad_a
            real, dimension(1:gpu_max_dim, 1:gpu_max_species) :: i_raw
            real :: y_corr(gpu_max_species)
            real :: sumJ(1:gpu_max_dim)
            real, dimension(1:gpu_max_dim) :: sumGradY
            real, dimension(1:gpu_max_dim, 1:gpu_max_species) :: gradY_all

            i = el_int(ii)
            iconsidered = i

            weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)

            do l = 1, ielem_ifca(i) !for all their faces
   !                 if (ielem_reorient(l,i).eq.0)then
               damp = lamx
               if (br2_yn == 2) damp = 0.00

               godflux2 = zero
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = angle1
               ny = angle2

               iqp = qp_line_n

               do ngp = 1, iqp !for all the gaussian quadrature points
                  facex = l
                  pointx = ngp

                  call calculate_interior_viscous(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, &
               & angle1, angle2, &
                     & nx, ny, nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright, lcvgrad, rcvgrad, &
                        & lcvgrad_t, &
                     & rcvgrad_t)

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                  !rotate wrt to normalvector of face and solve 1d riemann problem
                     call rotatef2d(n, cright_rot, cright, angle1, angle2)
                     call rotatef2d(n, cleft_rot, cleft, angle1, angle2)
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                     call lmacht2d(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                     call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb2d(n, cright, cright_rot, angle1, angle2)
                  end if

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)

                  viscl = 0.0e0; laml = 0.0e0
                  if (realgas == 1) then
                        call multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, &
                           & rg_difl, rg_enthl, rg_enthvbl, gammal)

                        laml(1) = mp_ktr_mix

                        mp_laml = mp_kve
                        viscl(1) = mp_mu_mix

                        rightv(1:nof_variables) = cright(1:nof_variables)
                        call multispecies_mixtures_rg(rightv, mp_mu_mix, mp_ktr_mix, mp_kve, &
                           & rg_difr, rg_enthr, rg_enthvbr, gammar)

                     ! mixture thermal conductivities (average)

                        laml(2) = mp_ktr_mix
                        mp_lamr = mp_kve
                        viscl(2) = mp_mu_mix
                  else
                     call get_visc_conduct(n, leftv, rightv, viscl, laml)
                  end if

                  if (turbulence .eq. 1) then
                     if (turbulencemodel .eq. 1) then
                        turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                     if (turbulencemodel .eq. 2) then
                        eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                        eddyfl(4:5) = lcvgrad(1, 1:2); eddyfl(6:7) = lcvgrad(2, 1:2)
                        eddyfl(8:9) = lcvgrad_t(1, 1:2)
                        eddyfl(10:11) = lcvgrad_t(2, 1:2)

                        eddyfr(1) = ielem_walldist(i); eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
                        eddyfr(4:5) = rcvgrad(1, 1:2); eddyfr(6:7) = rcvgrad(2, 1:2)
                        eddyfr(8:9) = rcvgrad_t(1, 1:2); eddyfr(10:11) = rcvgrad_t(2, 1:2)
                        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                  end if

                  taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz &
               & = zero; wx = &
                     & zero; wy = zero; wz = zero;
                  fxv = zero; fyv = zero; fzv = zero; rho12 = zero;
                  u12 = zero; v12 = zero; w12 = zero

   !

                  vdamp = (4.0 / 3.0) !*(( (viscl(1))+(viscl(2)))))
                  nall(1) = nx; nall(2) = ny

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
                  call cons2div(n, leftv, mp_pinfl, gammal)
                  call cons2div(n, rightv, mp_pinfr, gammar)
                  rho12 = oo2 * (leftv(1) + rightv(1))
                  rho12l = leftv(1); rho12r = rightv(1)
                  u12 = oo2 * (leftv(2) + rightv(2))
                  v12 = oo2 * (leftv(3) + rightv(3))

                  if ((realgas .eq. 1) .and. (nof_species .gt. 0)) then
                     do rg_i = 1, nof_species
                        idxy = dimensiona + 3 + rg_i
                        y_av(rg_i) = 0.50 * (leftv(idxy) + rightv(idxy))

                     end do
                  end if

                  do k = 1, nof_variables - 1
                     lcvgrad(k, 1:2) = ((lcvgrad(k, 1:2) + rcvgrad(k, 1:2)) / (2.00)) + damp * ((vdamp / abs(ielem_dih(i, &
                        & l))) * nall(1:2) * (rightv(k + 1) - leftv(k + 1)))
                  end do

               !now compute all the temperature gradients +real gas
                  if ((realgas .eq. 1) .or. (multispecies .eq. 1)) then

                     if (realgas .eq. 1) then

                     ! =============================================================
                     ! 1. compute mixture transport properties (left & right)
                     ! =============================================================
                        leftv(1:nof_variables) = cleft(1:nof_variables)


                        mp_ktr_mix_av = 0.50 * (laml(1) + laml(2))
                        if (turbulence == 1) mp_ktr_mix_av = 0.5e0 * (laml(3) + laml(4))
                        mp_lam_av = 0.50 * (mp_laml + mp_lamr)

                     ! mixture-averaged species diffusion coefficients
                        rg_dif_av(:) = 0.50 * (rg_difl(:) + rg_difr(:))

                     ! face-averaged species enthalpies
                        rg_enth_av(:) = 0.50 * (rg_enthl(:) + rg_enthr(:))
                        rg_enthvb_av(:) = 0.50 * (rg_enthvbl(:) + rg_enthvbr(:))

                     ! =============================================================
                     ! 2. turbulence shortcut
                     ! =============================================================
                        if (turbulence .eq. 1) then
                           q(1:2) = -0.50 * (laml(3) + laml(4)) * lcvgrad(dimensiona + 1, 1:2)

                        end if

                     ! =============================================================
                     ! 3. mixture-averaged diffusion fluxes (no pressure terms)
                     ! =============================================================
                        rg_sum_htr(1:dimensiona) = 0.00
                        rg_sum_hv(1:dimensiona) = 0.00
                        rg_qtr = 0.00; rg_qv = 0.00

                        sumi(1:dimensiona) = 0.00

                        if (rg_relax .ge. 1) then

                           sumi(1:dimensiona) = 0.00
                           sumJ(1:dimensiona) = 0.00
                           i_raw(1:dimensiona, 1:nof_species) = 0.00

                        ! -------------------------------------------------------------
                        ! Normalise face mass fractions used in diffusion correction.
                        ! This is essential because:
                        ! sum(J_k) = sum(i_raw_k) - sum(Y_k)*sum(i_raw_k)
                        ! so conservation requires sum(Y_k)=1 at the face.
                        ! -------------------------------------------------------------
                           sumY_face = 0.00

                           do rg_j = 1, nof_species
                              y_corr(rg_j) = max(y_av(rg_j), 0.00)
                              sumY_face = sumY_face + y_corr(rg_j)
                           end do

                           if (sumY_face .gt. 1.0e-30) then
                              do rg_j = 1, nof_species
                                 y_corr(rg_j) = y_corr(rg_j) / sumY_face
                              end do
                           else
                              do rg_j = 1, nof_species
                                 y_corr(rg_j) = 1.00 / real(nof_species)
                              end do
                           end if

                        ! -------------------------------------------------------------
                        ! Compute face mole fractions from corrected face mass fractions
                        ! -------------------------------------------------------------
                           molar_sum = 0.00
                           do rg_j = 1, nof_species
                              molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
                           end do

                           molar_sum = max(molar_sum, 1.0e-30)
                           mbar = 1.00 / molar_sum

                           do rg_j = 1, nof_species
                              x_av(rg_j) = (y_corr(rg_j) / rg_molm(rg_j)) / molar_sum
                           end do

                        ! -------------------------------------------------------------
                        ! Project species gradients so that sum_k grad(Y_k) = 0.
                        ! This removes a conservative but noisy pairwise diffusion mode
                        ! which can enter rhoEv through sum_k(hvib_k*J_k).
                        ! -------------------------------------------------------------
                           sumGradY(1:dimensiona) = 0.00

                           do rg_j = 1, nof_species
                              idxy = dimensiona + 2 + rg_j
                              gradY_all(1:dimensiona, rg_j) = lcvgrad(idxy, 1:dimensiona)
                              sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona, rg_j)
                           end do

                           do rg_j = 1, nof_species
                              gradY_all(1:dimensiona, rg_j) = gradY_all(1:dimensiona, rg_j) - &
                                 & y_corr(rg_j) * sumGradY(1:dimensiona)
                           end do

                        ! -------------------------------------------------------------
                        ! grad_a = grad(sum_k Y_k/M_k), using projected grad(Y_k).
                        ! -------------------------------------------------------------
                           grad_a(1:dimensiona) = 0.00

                           do rg_j = 1, nof_species
                              grady(1:dimensiona) = gradY_all(1:dimensiona, rg_j)

                              grad_a(1:dimensiona) = grad_a(1:dimensiona) + &
                                 & grady(1:dimensiona) / rg_molm(rg_j)
                           end do

                        ! -------------------------------------------------------------
                        ! Raw mixture-averaged flux:
                        ! i_raw_i = -rho * D_i * (M_i/Mmix) * grad(X_i)
                        ! -------------------------------------------------------------
                           do rg_i = 1, nof_species

                              grady(1:dimensiona) = gradY_all(1:dimensiona, rg_i)

                              gradx(1:dimensiona) = &
                                 & (grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
                                 & - (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona)) &
                                 & / (molar_sum * molar_sum)

                              i_raw(1:dimensiona, rg_i) = &
                                 & -rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)

                              sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona, rg_i)

                           end do

                        ! -------------------------------------------------------------
                        ! Mass-conserving corrected species diffusion flux:
                        ! J_i = i_raw_i - Y_i * sum_k(i_raw_k)
                        ! using corrected Y_i with sum(Y_i)=1.
                        ! -------------------------------------------------------------
                           do rg_i = 1, nof_species

                              y_face = y_corr(rg_i)

                              jl(1:dimensiona) = i_raw(1:dimensiona, rg_i) - &
                                 & y_face * sumi(1:dimensiona)

                              sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

                              fxv(dimensiona + 3 + rg_i) = fxv(dimensiona + 3 + rg_i) - jl(1)
                              fyv(dimensiona + 3 + rg_i) = fyv(dimensiona + 3 + rg_i) - jl(2)

                              rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
                                 & + rg_enth_av(rg_i) * jl(1:dimensiona)

                              rg_sum_hv(1:dimensiona) = rg_sum_hv(1:dimensiona) &
                                 & + rg_enthvb_av(rg_i) * jl(1:dimensiona)

                           end do

                        ! =============================================================
                        ! 4. conductive heat fluxes (fourier)
                        ! =============================================================
                           rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona + 1, 1:dimensiona)
                           rg_qv(1:dimensiona) = -mp_lam_av * lcvgrad(dimensiona + 2, 1:dimensiona)

                        end if

                     ! =============================================================
                     ! 5. total diffusive fluxes
                     ! =============================================================
                        q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
                           & + rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

                        qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

                     ! =============================================================
                     ! 6. add to conservative flux vectors
                     ! =============================================================
                        fxv(4) = fxv(4) - q(1)
                        fyv(4) = fyv(4) - q(2)

                        fxv(5) = fxv(5) - qvib(1)
                        fyv(5) = fyv(5) - qvib(2)

   !

                     end if !real gas ends

                     if (multispecies .eq. 1) then
                     !viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

                        if (turbulence .eq. 1) then
                           q(1:2) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:2))
                        else
                           q(1:2) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:2))
                        end if

                        fxv(4) = fxv(4) - q(1); fyv(4) = fyv(4) - q(2)
                     end if

                  end if

                  if (realgas == 0 .and. multispecies == 0) then

                     if (turbulence .eq. 1) then
                        q(1:2) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:2))
                     else
                        q(1:2) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:2))
                     end if
                  !q(1:2) =  - oo2* ((laml(1)+ (laml(2)))*lcvgrad(dimensiona+1,1:2))
                     fxv(4) = fxv(4) - q(1); fyv(4) = fyv(4) - q(2)

                  end if

               !left state derivatives
                  ux = lcvgrad(1, 1); uy = lcvgrad(1, 2)
                  vx = lcvgrad(2, 1); vy = lcvgrad(2, 2)
               ! determine taul!!

               ! tau_xx
                  taul(1, 1) = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy
               ! tau_yy
                  taul(2, 2) = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux
               ! tau_zz

               ! tau_xy
                  taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

               ! average and multiplay by viscosity
                  if (turbulence .eq. 1) then
                     tau = oo2 * (((viscl(1) + viscl(3))) + ((viscl(2) + viscl(4)))) * taul
                  else
                     tau = oo2 * ((viscl(1)) + (viscl(2))) * taul
                  end if

               ! now addition into momentum fluxes
                  do kc = 2, 3
                     fxv(kc) = fxv(kc) + tau(1, kc - 1)
                     fyv(kc) = fyv(kc) + tau(2, kc - 1)

                  end do

                  fxv(4) = fxv(4) + u12 * tau(1, 1) + v12 * tau(1, 2)
                  fyv(4) = fyv(4) + u12 * tau(2, 1) + v12 * tau(2, 2)

   !
                  hllcflux(1:nof_variables) = nx * fxv(1:nof_variables) + ny * fyv(1:nof_variables)

                  if (dg .eq. 1) then

                     rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, -1.00)

                  else

                     godflux2(1:nof_variables) = godflux2(1:nof_variables) + (hllcflux(1:nof_variables) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (turbulence .eq. 1) then
                        hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & ((oo2 * (viscl(1) + viscl(2))) + (oo2 * (viscl(3) + viscl(4)))) * &
                           & (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + &
               & passivescalar, &
                           & 1)) * oo2 * nx) + &
                           & ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + &
                              & passivescalar, &
                           & 2)) * oo2 * ny))
                     else
                        hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & ((oo2 * (viscl(1) + viscl(2)))) * &
                           & (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + &
               & passivescalar, &
                           & 1)) * oo2 * nx) + &
                           & ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + &
                              & passivescalar, &
                           & 2)) * oo2 * ny))

                     end if
                     if (turbulencemodel .eq. 1) then
                        hllcflux(5) = hllcflux(5) / sigma
                     end if
                     godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar) + &
                        & (hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))
                  end if

               end do

               rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) - godflux2(1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) - &
                     & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

               end if
            end do

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_diffusive2d_inner_cell

   subroutine calculate_fluxeshi_diffusive2d_bound_cell(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, l, ngp, iqp, kc, b_code
      integer :: iconsidered, facex, pointx, k, rg_i, rg_j, idxy, icompute
      real :: angle1, angle2, nx, ny, nz
      real :: mp_pinfl, mp_pinfr, gammal, gammar, rho12l, rho12r
      real :: ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp, y_face
      real :: mp_mu_mix, mp_ktr_mix, mp_kve, mp_laml, mp_lamr
      real :: mp_lam_av, mp_ktr_mix_av
      real :: sumY_face
      real :: molar_sum, mbar

      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(br2_yn, catalytic_wall, dg, dimensiona, lamx, lmach, multispecies, n, nof_bounded, nof_species) &
!$omp& firstprivate(nof_variables, oo2, passivescalar, qp_line_n, realgas, rg_relax, sigma, turbulence) &
!$omp& firstprivate(turbulenceequations, turbulencemodel, zero) &
!$omp& map(alloc: el_bnd, ielem_dih, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_surf, ielem_walldist) &
!$omp& map(alloc: rg_molm, rhs_val, rhst_val, weights_l) &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(icompute, iconsidered, idxy, ii, iqp, k, kc, l, mbar) &
!$omp& private(molar_sum, mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix, mp_pinfl) &
!$omp& private(mp_pinfr, ngp, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j) &
!$omp& private(rho12, rho12l, rho12r) &
!$omp& private(sumy_face, u12, ux, uy, uz, v12, vdamp, vx, vy) &
!$omp& private(vz, w12, wx, wy, wz, y_face) &
!$omp& map(alloc: bleed_end, bleed_plenum, bleed_porosity, bleed_start, bound_offset, boundhir, boundhir_dg) &
!$omp& map(alloc: catalytic_con, gamma_in, ibound_icode, ielem_bleedn, ielem_ibounds, ielem_idegfree, ielem_indexf) &
!$omp& map(alloc: ielem_ineigh, ielem_ineighb, ielem_ineighn, ielem_inter_id, ielem_iorder, ielem_nodes_faces) &
!$omp& map(alloc: ielem_qface, ielem_totvolume, inoder4_cord, integ_basis_dg_value, integ_basis_value) &
!$omp& map(alloc: integ_basis_valuec, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in) &
!$omp& map(alloc: rec_qpoints, rec_uleft, rec_uleft_dg, rec_uleftturb, rec_uleftturbv, rec_uleftv, rec_volume) &
!$omp& map(alloc: rg_hzero, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, rhs_valdg) &
!$omp& map(alloc: u_ct_val)
#else
!$omp do &
!$omp& private(angle1, angle2, b_code, damp) &
!$omp& private(facex, gammal, gammar, i) &
!$omp& private(icompute, iconsidered, idxy, ii, iqp, k, kc, l, mbar) &
!$omp& private(molar_sum, mp_ktr_mix, mp_ktr_mix_av, mp_kve, mp_lam_av, mp_laml, mp_lamr, mp_mu_mix, mp_pinfl) &
!$omp& private(mp_pinfr, ngp, nx, ny, nz, pointx) &
!$omp& private(rg_i, rg_j) &
!$omp& private(rho12, rho12l, rho12r) &
!$omp& private(sumy_face, u12, ux, uy, uz, v12, vdamp, vx, vy) &
!$omp& private(vz, w12, wx, wy, wz, y_face)
#endif
      do ii = 1, nof_bounded
         block
            real, dimension(1:gpu_max_nvar_total) :: godflux2, rhllcflux, hllcflux
            real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
            real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
            real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
            real, dimension(1:gpu_max_dim) :: pox, poy, poz
            real, dimension(1:gpu_max_qp_face) :: weights_temp
            real, dimension(1:4) :: viscl, laml
            real, dimension(1:gpu_max_species) :: y_av, x_av
            real, dimension(1:2) :: turbmv
            real, dimension(1) :: etvm
            real, dimension(1:20) :: eddyfl, eddyfr
            real, dimension(1:gpu_max_nvar, 1:gpu_max_dim) :: lcvgrad, rcvgrad
            real, dimension(gpu_max_extra_transport, 1:gpu_max_dim) :: lcvgrad_t, rcvgrad_t
            real, dimension(1:gpu_max_nvar) :: fxv, fyv, fzv
            real, dimension(2, 2) :: taul, taur, tau
            real, dimension(2) :: q, nall, qvib
            real, dimension(1:gpu_max_species) :: rg_difl, rg_difr, rg_enthl, rg_enthr, rg_enthvbl, rg_enthvbr
            real, dimension(1:gpu_max_dim) :: jl
            real, dimension(1:gpu_max_dim) :: rg_sum_htr, rg_sum_hv, grady, rg_qv, rg_qtr
            real, dimension(1:gpu_max_species) :: rg_dif_av, rg_enth_av, rg_enthvb_av
            real, dimension(1:gpu_max_dim) :: sumi, gradx, grad_a
            real, dimension(1:gpu_max_dim, 1:gpu_max_species) :: i_raw
            real :: y_corr(gpu_max_species)
            real :: sumJ(1:gpu_max_dim)
            real, dimension(1:gpu_max_dim) :: sumGradY
            real, dimension(1:gpu_max_dim, 1:gpu_max_species) :: gradY_all

            i = el_bnd(ii)
            iconsidered = i

            weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)

            do l = 1, ielem_ifca(i)
               damp = lamx
               if (br2_yn == 2) damp = 0.00
               b_code = 0
               godflux2 = zero
               angle1 = ielem_faceanglex(i, l)
               angle2 = ielem_faceangley(i, l)
               nx = angle1
               ny = angle2

               iqp = qp_line_n

               do ngp = 1, iqp

                  facex = l
                  pointx = ngp

                  call calculate_bounded_viscous2d(n, b_code, iconsidered, facex, pointx, leftv, rightv, pox, poy, poz, &
               & angle1, angle2, &
                     & nx, ny, nz, cturbl, cturbr, cright_rot, cleft_rot, srf_speedrot, cleft, cright, lcvgrad, rcvgrad, &
                        & lcvgrad_t, &
                     & rcvgrad_t)

                  if ((lmach .eq. 1)) then !application of the low mach number correction
                  !rotate wrt to normalvector of face and solve 1d riemann problem
                     call rotatef2d(n, cright_rot, cright, angle1, angle2)
                     call rotatef2d(n, cleft_rot, cleft, angle1, angle2)
                     leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                     call lmacht2d(n, leftv, rightv)
                     cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
                     call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb2d(n, cright, cright_rot, angle1, angle2)
                  end if

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)

                  viscl = 0.0e0; laml = 0.0e0
                  if (realgas == 1) then
                        call multispecies_mixtures_rg(leftv, mp_mu_mix, mp_ktr_mix, mp_kve, &
                           & rg_difl, rg_enthl, rg_enthvbl, gammal)

                        laml(1) = mp_ktr_mix

                        mp_laml = mp_kve
                        viscl(1) = mp_mu_mix

                        rightv(1:nof_variables) = cright(1:nof_variables)
                        call multispecies_mixtures_rg(rightv, mp_mu_mix, mp_ktr_mix, mp_kve, &
                           & rg_difr, rg_enthr, rg_enthvbr, gammar)

                     ! mixture thermal conductivities (average)

                        laml(2) = mp_ktr_mix
                        mp_lamr = mp_kve
                        viscl(2) = mp_mu_mix
                  else
                     call get_visc_conduct(n, leftv, rightv, viscl, laml)
                  end if

                  if (turbulence .eq. 1) then
                     if (turbulencemodel .eq. 1) then
                        turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
                        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                        if (b_code .eq. 4) then
                           viscl(3:4) = zero
                           ! These slots hold total conductivity, so retain the molecular part.
                           laml(3:4) = laml(1:2)
                        end if
                     end if
                     if (turbulencemodel .eq. 2) then
                        eddyfl(1) = ielem_walldist(i); eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
                        eddyfl(4:5) = lcvgrad(1, 1:2); eddyfl(6:7) = lcvgrad(2, 1:2)
                        eddyfl(8:9) = lcvgrad_t(1, 1:2)
                        eddyfl(10:11) = lcvgrad_t(2, 1:2)

                        eddyfr(1) = ielem_walldist(i); eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
                        eddyfr(4:5) = rcvgrad(1, 1:2); eddyfr(6:7) = rcvgrad(2, 1:2)
                        eddyfr(8:9) = rcvgrad_t(1, 1:2); eddyfr(10:11) = rcvgrad_t(2, 1:2)
                        call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
                     end if
                  end if

                  taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz &
               & = zero; wx = &
                     & zero; wy = zero; wz = zero;
                  fxv = zero; fyv = zero; fzv = zero; rho12 = zero;
                  u12 = zero; v12 = zero; w12 = zero

                  if ((b_code .gt. 0)) then
                     damp = zero
                  end if

                  vdamp = (4.0 / 3.0) !*(( (viscl(1))+(viscl(2)))))
                  nall(1) = nx; nall(2) = ny

                  leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
                  call cons2div(n, leftv, mp_pinfl, gammal)
                  call cons2div(n, rightv, mp_pinfr, gammar)
                  rho12 = oo2 * (leftv(1) + rightv(1))
                  rho12l = leftv(1); rho12r = rightv(1)
                  u12 = oo2 * (leftv(2) + rightv(2))
                  v12 = oo2 * (leftv(3) + rightv(3))

                  if ((realgas .eq. 1) .and. (nof_species .gt. 0)) then
                     do rg_i = 1, nof_species
                        idxy = dimensiona + 3 + rg_i
                        y_av(rg_i) = 0.50 * (leftv(idxy) + rightv(idxy))

                     end do
                  end if

                  do k = 1, nof_variables - 1
                     lcvgrad(k, 1:2) = ((lcvgrad(k, 1:2) + rcvgrad(k, 1:2)) / (2.00)) + damp * ((vdamp / abs(ielem_dih(i, &
                        & l))) * nall(1:2) * (rightv(k + 1) - leftv(k + 1)))
                  end do

                  icompute = 0

               !now compute all the temperature gradients +real gas
                  if ((realgas .eq. 1) .or. (multispecies .eq. 1)) then

                     if (realgas .eq. 1) then

                     ! =============================================================
                     ! 1. compute mixture transport properties (left & right)
                     ! =============================================================
                        leftv(1:nof_variables) = cleft(1:nof_variables)


                        mp_ktr_mix_av = 0.50 * (laml(1) + laml(2))
                        if (turbulence == 1) mp_ktr_mix_av = 0.5e0 * (laml(3) + laml(4))
                        mp_lam_av = 0.50 * (mp_laml + mp_lamr)

                     ! mixture-averaged species diffusion coefficients
                        rg_dif_av(:) = 0.50 * (rg_difl(:) + rg_difr(:))

                     ! face-averaged species enthalpies
                        rg_enth_av(:) = 0.50 * (rg_enthl(:) + rg_enthr(:))
                        rg_enthvb_av(:) = 0.50 * (rg_enthvbl(:) + rg_enthvbr(:))

                     ! =============================================================
                     ! 2. turbulence shortcut
                     ! =============================================================
                        if (turbulence .eq. 1) then
                           q(1:2) = -0.50 * (laml(3) + laml(4)) * lcvgrad(dimensiona + 1, 1:2)

                        end if

                     ! =============================================================
                     ! 3. mixture-averaged diffusion fluxes (no pressure terms)
                     ! =============================================================
                        rg_sum_htr(1:dimensiona) = 0.00
                        rg_sum_hv(1:dimensiona) = 0.00
                        rg_qtr = 0.00; rg_qv = 0.00

                        sumi(1:dimensiona) = 0.00

                        if (rg_relax .ge. 1) then
                           if ((b_code .eq. 4) .and. (catalytic_wall .eq. 0)) then
                              icompute = 1
                           end if

                           if (icompute .eq. 0) then

                              sumi(1:dimensiona) = 0.00
                              sumJ(1:dimensiona) = 0.00
                              i_raw(1:dimensiona, 1:nof_species) = 0.00

                           ! -------------------------------------------------------------
                           ! Normalise face mass fractions used in diffusion correction.
                           ! This is essential because:
                           ! sum(J_k) = sum(i_raw_k) - sum(Y_k)*sum(i_raw_k)
                           ! so conservation requires sum(Y_k)=1 at the face.
                           ! -------------------------------------------------------------
                              sumY_face = 0.00

                              do rg_j = 1, nof_species
                                 y_corr(rg_j) = max(y_av(rg_j), 0.00)
                                 sumY_face = sumY_face + y_corr(rg_j)
                              end do

                              if (sumY_face .gt. 1.0e-30) then
                                 do rg_j = 1, nof_species
                                    y_corr(rg_j) = y_corr(rg_j) / sumY_face
                                 end do
                              else
                                 do rg_j = 1, nof_species
                                    y_corr(rg_j) = 1.00 / real(nof_species)
                                 end do
                              end if

                           ! -------------------------------------------------------------
                           ! Compute face mole fractions from corrected face mass fractions
                           ! -------------------------------------------------------------
                              molar_sum = 0.00
                              do rg_j = 1, nof_species
                                 molar_sum = molar_sum + y_corr(rg_j) / rg_molm(rg_j)
                              end do

                              molar_sum = max(molar_sum, 1.0e-30)
                              mbar = 1.00 / molar_sum

                              do rg_j = 1, nof_species
                                 x_av(rg_j) = (y_corr(rg_j) / rg_molm(rg_j)) / molar_sum
                              end do

                           ! -------------------------------------------------------------
                           ! Project species gradients so that sum_k grad(Y_k) = 0.
                           ! This removes a conservative but noisy pairwise diffusion mode
                           ! which can enter rhoEv through sum_k(hvib_k*J_k).
                           ! -------------------------------------------------------------
                              sumGradY(1:dimensiona) = 0.00

                              do rg_j = 1, nof_species
                                 idxy = dimensiona + 2 + rg_j
                                 gradY_all(1:dimensiona, rg_j) = lcvgrad(idxy, 1:dimensiona)
                                 sumGradY(1:dimensiona) = sumGradY(1:dimensiona) + gradY_all(1:dimensiona, rg_j)
                              end do

                              do rg_j = 1, nof_species
                                 gradY_all(1:dimensiona, rg_j) = gradY_all(1:dimensiona, rg_j) - &
                                    & y_corr(rg_j) * sumGradY(1:dimensiona)
                              end do

                           ! -------------------------------------------------------------
                           ! grad_a = grad(sum_k Y_k/M_k), using projected grad(Y_k).
                           ! -------------------------------------------------------------
                              grad_a(1:dimensiona) = 0.00

                              do rg_j = 1, nof_species
                                 grady(1:dimensiona) = gradY_all(1:dimensiona, rg_j)

                                 grad_a(1:dimensiona) = grad_a(1:dimensiona) + &
                                    & grady(1:dimensiona) / rg_molm(rg_j)
                              end do

                           ! -------------------------------------------------------------
                           ! Raw mixture-averaged flux:
                           ! i_raw_i = -rho * D_i * (M_i/Mmix) * grad(X_i)
                           ! -------------------------------------------------------------
                              do rg_i = 1, nof_species

                                 grady(1:dimensiona) = gradY_all(1:dimensiona, rg_i)

                                 gradx(1:dimensiona) = &
                                    & (grady(1:dimensiona) / rg_molm(rg_i) * molar_sum &
                                    & - (y_corr(rg_i) / rg_molm(rg_i)) * grad_a(1:dimensiona)) &
                                    & / (molar_sum * molar_sum)

                                 i_raw(1:dimensiona, rg_i) = &
                                    & -rho12 * rg_dif_av(rg_i) * (rg_molm(rg_i) / mbar) * gradx(1:dimensiona)

                                 sumi(1:dimensiona) = sumi(1:dimensiona) + i_raw(1:dimensiona, rg_i)

                              end do

                           ! -------------------------------------------------------------
                           ! Mass-conserving corrected species diffusion flux:
                           ! J_i = i_raw_i - Y_i * sum_k(i_raw_k)
                           ! using corrected Y_i with sum(Y_i)=1.
                           ! -------------------------------------------------------------
                              do rg_i = 1, nof_species

                                 y_face = y_corr(rg_i)

                                 jl(1:dimensiona) = i_raw(1:dimensiona, rg_i) - &
                                    & y_face * sumi(1:dimensiona)

                                 sumJ(1:dimensiona) = sumJ(1:dimensiona) + jl(1:dimensiona)

                                 fxv(dimensiona + 3 + rg_i) = fxv(dimensiona + 3 + rg_i) - jl(1)
                                 fyv(dimensiona + 3 + rg_i) = fyv(dimensiona + 3 + rg_i) - jl(2)

                                 rg_sum_htr(1:dimensiona) = rg_sum_htr(1:dimensiona) &
                                    & + rg_enth_av(rg_i) * jl(1:dimensiona)

                                 rg_sum_hv(1:dimensiona) = rg_sum_hv(1:dimensiona) &
                                    & + rg_enthvb_av(rg_i) * jl(1:dimensiona)

                              end do

                           ! Optional diagnostic

                           end if
   !                                                                                                         end if

                        ! =============================================================
                        ! 4. conductive heat fluxes (fourier)
                        ! =============================================================
                           rg_qtr(1:dimensiona) = -mp_ktr_mix_av * lcvgrad(dimensiona + 1, 1:dimensiona)
                           rg_qv(1:dimensiona) = -mp_lam_av * lcvgrad(dimensiona + 2, 1:dimensiona)

                        end if

                     ! =============================================================
                     ! 5. total diffusive fluxes
                     ! =============================================================
                        q(1:dimensiona) = rg_qtr(1:dimensiona) + rg_qv(1:dimensiona) &
                           & + rg_sum_htr(1:dimensiona) + rg_sum_hv(1:dimensiona)

                        qvib(1:dimensiona) = rg_qv(1:dimensiona) + rg_sum_hv(1:dimensiona)

                     ! =============================================================
                     ! 6. add to conservative flux vectors
                     ! =============================================================
                        fxv(4) = fxv(4) - q(1)
                        fyv(4) = fyv(4) - q(2)

                        fxv(5) = fxv(5) - qvib(1)
                        fyv(5) = fyv(5) - qvib(2)

   !

                     end if !real gas ends

                     if (multispecies .eq. 1) then
                     !viscous stress + heat conduction in momentum/energy equation only not on each individual component on allaire

                        if (turbulence .eq. 1) then
                           q(1:2) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:2))
                        else
                           q(1:2) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:2))
                        end if

                        fxv(4) = fxv(4) - q(1); fyv(4) = fyv(4) - q(2)
                     end if

                  end if

                  if (realgas == 0 .and. multispecies == 0) then
                     if (turbulence .eq. 1) then
                        q(1:2) = -oo2 * ((laml(3) + (laml(4))) * lcvgrad(dimensiona + 1, 1:2))
                     else
                        q(1:2) = -oo2 * ((laml(1) + (laml(2))) * lcvgrad(dimensiona + 1, 1:2))
                     end if

                     fxv(4) = fxv(4) - q(1); fyv(4) = fyv(4) - q(2)

                  end if

               !left state derivatives
                  ux = lcvgrad(1, 1); uy = lcvgrad(1, 2)
                  vx = lcvgrad(2, 1); vy = lcvgrad(2, 2)
               ! determine taul!!

               ! tau_xx
                  taul(1, 1) = (4.00 / 3.00) * ux - (2.00 / 3.00) * vy
               ! tau_yy
                  taul(2, 2) = (4.00 / 3.00) * vy - (2.00 / 3.00) * ux
               ! tau_zz

               ! tau_xy
                  taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

               ! average and multiplay by viscosity
                  if (turbulence .eq. 1) then
                     tau = oo2 * (((viscl(1) + viscl(3))) + ((viscl(2) + viscl(4)))) * taul
                  else
                     tau = oo2 * ((viscl(1)) + (viscl(2))) * taul
                  end if

               ! now addition into momentum fluxes
                  do kc = 2, 3
                     fxv(kc) = fxv(kc) + tau(1, kc - 1)
                     fyv(kc) = fyv(kc) + tau(2, kc - 1)

                  end do

                  fxv(4) = fxv(4) + u12 * tau(1, 1) + v12 * tau(1, 2)
                  fyv(4) = fyv(4) + u12 * tau(2, 1) + v12 * tau(2, 2)

   !
                  hllcflux(1:nof_variables) = nx * fxv(1:nof_variables) + ny * fyv(1:nof_variables)

                  if (realgas .eq. 1) then
                     if ((b_code .eq. 4) .and. (catalytic_wall .eq. 0)) then
                        hllcflux(6:nof_variables) = zero
                     end if
                  end if

                  if (dg .eq. 1) then

                     rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)

                     call dg_surf_flux_accumulate_rhs(n, iconsidered, facex, pointx, weights_temp, rhllcflux, -1.00)

                  else

                     godflux2(1:nof_variables) = godflux2(1:nof_variables) + (hllcflux(1:nof_variables) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                     if (turbulence .eq. 1) then
                        hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & ((oo2 * (viscl(1) + viscl(2))) + (oo2 * (viscl(3) + viscl(4)))) * &
                           & (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + &
               & passivescalar, &
                           & 1)) * oo2 * nx) + &
                           & ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + &
                              & passivescalar, &
                           & 2)) * oo2 * ny))
                     else
                        hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
                           & ((oo2 * (viscl(1) + viscl(2)))) * &
                           & (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + &
               & passivescalar, &
                           & 1)) * oo2 * nx) + &
                           & ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + &
                              & passivescalar, &
                           & 2)) * oo2 * ny))

                     end if
                     if (turbulencemodel .eq. 1) then
                        hllcflux(5) = hllcflux(5) / sigma
                     end if

                     if ((b_code .eq. 3)) then
                        hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero

                     end if

                     godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar) + &
                        & (hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                        & ielem_surf(i, l)))
                  end if

               end do

               rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) - godflux2(1:nof_variables)

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) - &
                     & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

               end if
            end do

         end block
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine calculate_fluxeshi_diffusive2d_bound_cell



   subroutine calculate_fluxeshi_convective_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
      implicit none
      integer, intent(in) :: n
      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_interior) &
!$omp& map(alloc: el_int, gamma_in, ielem_faceanglex, ielem_faceangley, ielem_facediss, ielem_ifca, ielem_ineigh) &
!$omp& map(alloc: ielem_ineighn, ielem_mood, ielem_recalc, ielem_surf, ielem_types_faces, mp_janaf, mp_m, mp_pinf) &
!$omp& map(alloc: mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_mrf, rec_uleft, rec_uleftturb, rg_hzero, rg_molm, rhs_val) &
!$omp& map(alloc: rhs_valdg, rhst_val, u_c_val, u_ct_val, weights_q, weights_t)
#else
!$omp barrier
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_interior !for all the interior elements
         call calculate_fluxeshi_convective_mood_inner_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_bounded) &
!$omp& map(alloc: bound_offset, boundhir, boundhirm, catalytic_con, el_bnd, gamma_in, halo_offset, ibound_icode) &
!$omp& map(alloc: ielem_faceanglex, ielem_faceangley, ielem_facediss, ielem_ibounds, ielem_ifca, ielem_indexf) &
!$omp& map(alloc: ielem_indexi, ielem_ineigh, ielem_ineighb, ielem_ineighn, ielem_inter_id, ielem_mood) &
!$omp& map(alloc: ielem_nodes_faces, ielem_qface, ielem_recalc, ielem_surf, ielem_types_faces, ielem_walldist) &
!$omp& map(alloc: inoder4_cord, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in) &
!$omp& map(alloc: rec_ihexl, rec_ihexn, rec_local, rec_mrf, rec_uleft, rec_uleftturb, rg_hzero, rg_molm, rg_thetag) &
!$omp& map(alloc: rg_vf, rgs_ab, rgs_bb, rgs_cb, rhs_val, rhs_valdg, rhst_val, solhir, u_c_val, u_ct_val, weights_q) &
!$omp& map(alloc: weights_t)
#else
!$omp barrier
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_bounded
         call calculate_fluxeshi_convective_mood_bound_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine calculate_fluxeshi_convective_mood

   subroutine calculate_fluxeshi_convective_mood_inner_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2, dg_vol_rec, rhllcflux, hllcflux
      integer :: i, l, ngp, kmaxe, iqp, ikas, igoflux, icaseb, jx, jx2, b_code
      real :: sum_detect, norms, tempxx
      real, dimension(1:gpu_max_qp_face) :: weights_temp
      integer :: iconsidered, facex, pointx, n_node
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
      real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:gpu_max_nvar) :: srf_speed
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      integer :: ibfc, nfx, lfx, rowfx, nf, lf, rowf

      i = el_int(ii)
      iconsidered = i
      if (ielem_recalc(i) .eq. 1) then
         rhs_val(i, :) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst_val(i, :) = zero
      end if

      mp_source3 = zero
      b_code = 0
      if (ielem_recalc(i) .eq. 1) then

         do l = 1, ielem_ifca(i) !for all their faces
!                                      if (ielem_reorient(l,i).eq.0)then
            godflux2 = zero
            mp_source2 = zero
            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = (cos(angle1) * sin(angle2))
            ny = (sin(angle1) * sin(angle2))
            nz = (cos(angle2))
            if (ielem_types_faces(i, l) .eq. 5) then
               iqp = qp_quad_n
               weights_temp(1:iqp) = weights_q(1:iqp)
            else
               iqp = qp_triangle_n
               weights_temp(1:iqp) = weights_t(1:iqp)
            end if
            do ngp = 1, iqp !for all the gaussian quadrature points

               cleft(1:nof_variables) = rec_uleft(i, 1:nof_variables, l, ngp) !left mean flow state
               if ((cascade .eq. 2) .and. (ielem_mood(i) .ge. 1)) then
                  cleft(1:nof_variables) = u_c_val(i, 3, 1:nof_variables)
               end if

               cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)

               if ((cascade .eq. 2) .and. (ielem_mood(ielem_ineighn(i, l)) .ge. 1)) then
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 3, 1:nof_variables)
               end if

               !right mean flow state

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then

                     !left additional equations flow state
                     cturbl(1:turbulenceequations + passivescalar) = rec_uleftturb(i, 1:turbulenceequations + passivescalar, l, ngp)
                     !right additional equations flow state
                     cturbr(1:turbulenceequations + passivescalar) = rec_uleftturb(ielem_ineigh(i, l), &
                        & 1:turbulenceequations + passivescalar, ielem_ineighn(i, l), ngp)
                  else
                     cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)
                  end if

                  cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations &
                     & + passivescalar)
                  cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations &
                     & + passivescalar)
               end if

               call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

               if ((lmach .eq. 1)) then !application of the low mach number correction
                  leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                  call lmacht(n, leftv, rightv)
                  cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               end if

               select case (iriemann)

               case (1) !hllc

                  call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (2) !rusanov

                  call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (9) !hll

                  call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (3) !roe

                  call rotateb(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb(n, cright, cright_rot, angle1, angle2)
                  call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  rhllcflux = hllcflux

               case (4) !roe

                  call rotateb(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb(n, cright, cright_rot, angle1, angle2)
                  call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  rhllcflux = hllcflux

               case (5) !roe

                  call rotateb(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb(n, cright, cright_rot, angle1, angle2)
                  call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  rhllcflux = hllcflux

               end select

               godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables) * (weights_temp(ngp) * &
                  & ielem_surf(i, l)))
               if (multispecies .eq. 1) then
                  mp_source2 = mp_source2 + mp_source1 * (weights_temp(ngp) * ielem_surf(i, l))
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 0) then !first order upwind flux

                     norms = (nx * (u_c_val(i, 1, 2) / u_c_val(i, 1, 1))) &
                        & + (ny * (u_c_val(i, 1, 3) / u_c_val(i, 1, 1))) &
                        & + (nz * (u_c_val(i, 1, 4) / u_c_val(i, 1, 1)))
                     if (norms .ge. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbl(1:turbulenceequations + passivescalar)
                     end if
                     if (norms .lt. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbr(1:turbulenceequations + passivescalar)
                     end if

                  end if

                  godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                     & 1:nof_variables + turbulenceequations + passivescalar) + &
                     & (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                     & ielem_surf(i, l)))

               end if

            end do

            rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) + godflux2(1:nof_variables)
            if (multispecies .eq. 1) then
               mp_source3 = mp_source3 + mp_source2
            end if

!                                       rhs_val(ielem_ineigh(i,l),1:nof_variables)=rhs_val(ielem_ineigh(i,l),1:nof_variables)-godflux2(1:nof_variables)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) + &
                  & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!                                      rhst_val(ielem_ineigh(i,l),1:turbulenceequations+passivescalar)=rhst_val(ielem_ineigh(i,l),1:turbulenceequations+passivescalar)-&
!                                     godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if
!                                      end if
         end do
         if (multispecies .eq. 1) then
            call subtract_volume_fraction_source(i, mp_source3, .false.)

         end if
      end if
   end subroutine calculate_fluxeshi_convective_mood_inner_cell

   subroutine calculate_fluxeshi_convective_mood_bound_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2, dg_vol_rec, rhllcflux, hllcflux
      integer :: i, l, ngp, kmaxe, iqp, ikas, igoflux, icaseb, jx, jx2, b_code
      real :: sum_detect, norms, tempxx
      real, dimension(1:gpu_max_qp_face) :: weights_temp
      integer :: iconsidered, facex, pointx, n_node
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
      real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:gpu_max_nvar) :: srf_speed
      real, dimension(1:8, 1:gpu_max_dim) :: vext
      real, dimension(1:8, 1:gpu_max_dim) :: nodes_list
      real, dimension(3) :: cords
      integer :: ibfc, nfx, lfx, rowfx, nf, lf, rowf

      i = el_bnd(ii)
      iconsidered = i

      if (ielem_recalc(i) .eq. 1) then
         rhs_val(i, :) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst_val(i, :) = zero
      end if

      mp_source3 = zero
      if (ielem_recalc(i) .eq. 1) then

         do l = 1, ielem_ifca(i)
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

            angle1 = ielem_faceanglex(i, l)
            angle2 = ielem_faceangley(i, l)
            nx = (cos(angle1) * sin(angle2))
            ny = (sin(angle1) * sin(angle2))
            nz = (cos(angle2))

            if (ielem_types_faces(i, l) .eq. 5) then
               iqp = qp_quad_n
               weights_temp(1:iqp) = weights_q(1:iqp)
               n_node = 4
            else
               iqp = qp_triangle_n
               weights_temp(1:iqp) = weights_t(1:iqp)
               n_node = 3
            end if
            godflux2 = zero

            mp_source2 = zero

            do ngp = 1, iqp
               b_code = 0

               cleft(1:nof_variables) = rec_uleft(i, 1:nof_variables, l, ngp)
               if ((cascade .eq. 2) .and. (ielem_mood(i) .ge. 1)) then
                  cleft(1:nof_variables) = u_c_val(i, 3, 1:nof_variables)

               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then

                     cturbl(1:turbulenceequations + passivescalar) = rec_uleftturb(i, 1:turbulenceequations + passivescalar, l, ngp)

                  else
                     cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
                  end if
               end if

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu

                        cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)

                        if ((cascade .eq. 2) .and. (ielem_mood(ielem_ineigh(i, l)) .ge. 1)) then
                           cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 3, 1:nof_variables)

                        end if

                        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                           if (icoupleturb .eq. 1) then
                              if (cascade .eq. 1) then
                                 !right additional equations flow state
                                 cturbr(1:turbulenceequations + passivescalar) = rec_uleftturb(ielem_ineigh(i, l), &
                                    & 1:turbulenceequations + passivescalar, ielem_ineighn(i, l), ngp)
                              else
                                 cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                                    & 1:turbulenceequations + passivescalar)

                              end if

                           else
                              cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                                 & 1:turbulenceequations + passivescalar)
                           end if
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

                        call boundarys(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & cturbl, cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                        cright(1:nof_variables) = rightv(1:nof_variables)

                     end if
                  else

                     cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)
                     if ((cascade .eq. 2) .and. (ielem_mood(ielem_ineigh(i, l)) .ge. 1)) then
                        cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 3, 1:nof_variables)

                     end if

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                        if (icoupleturb .eq. 1) then

                           !right additional equations flow state
                           cturbr(1:turbulenceequations + passivescalar) = rec_uleftturb(ielem_ineigh(i, l), &
                              & 1:turbulenceequations + passivescalar, ielem_ineighn(i, l), ngp)

                        else
                           cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                              & 1:turbulenceequations + passivescalar)
                        end if
                     end if

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu

!                                                                           cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                        rowfx = bound_offset(nfx) + lfx - 1
                        cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)

                        if ((cascade .eq. 2) .and. (boundhirm(rowfx) .ge. 1)) then

                           nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                           lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                           rowf = halo_offset(nf) + lf - 1
                           cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

!
!                                      cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                          (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)
!
                        end if

                        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                           if (icoupleturb .eq. 1) then

!                                                                            cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                              nfx = ielem_ineighn(i, l)
                              lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                              rowfx = bound_offset(nfx) + lfx - 1
                              cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                                 & turbulenceequations + passivescalar)

                           else
!                                                                          cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                              nfx = ielem_ineighn(i, l)
                              lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                              rowfx = bound_offset(nfx) + lfx - 1
                              cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                                 & turbulenceequations + passivescalar)
                           end if
                        end if

                     end if
                  else

!                                                                   cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                     rowfx = bound_offset(nfx) + lfx - 1

                     cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)

                     if ((cascade .eq. 2) .and. (boundhirm(rowfx) .ge. 1)) then

                        nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                        lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                        rowf = halo_offset(nf) + lf - 1
                        cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

!                                                                   cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                         (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                     end if

!
                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                        if (icoupleturb .eq. 1) then

!                                                                            cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)

                        else
!                                                                          cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)
                        end if
                     end if

!
                  end if
               end if

               call rotatef(n, cright_rot, cright, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotatef(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
!

               if ((lmach .eq. 1)) then !application of the low mach number correction
                  leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                  call lmacht(n, leftv, rightv)
                  cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations &
                     & + passivescalar)
                  cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations &
                     & + passivescalar)
               end if

               select case (iriemann)

               case (1) !hllc

                  call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb(n, rhllcflux, hllcflux, angle1, angle2)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if
               case (2) !rusanov

                  call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (3) !roe

                  call rotateb(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb(n, cright, cright_rot, angle1, angle2)
                  call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  rhllcflux = hllcflux

               case (9) !hll

                  call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if
               case (4) !roe

                  call rotateb(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb(n, cright, cright_rot, angle1, angle2)

                  if (b_code .le. 0) then

                     call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  else

                     call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  end if

                  rhllcflux = hllcflux

               case (5) !roe

                  call rotateb(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb(n, cright, cright_rot, angle1, angle2)

                  if (b_code .le. 0) then

                     call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  else

                     call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)

                  end if

                  rhllcflux = hllcflux

               end select

               godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables) * (weights_temp(ngp) * &
                  & ielem_surf(i, l)))
               if (multispecies .eq. 1) then
                  mp_source2 = mp_source2 + mp_source1 * (weights_temp(ngp) * ielem_surf(i, l))
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 0) then !first order upwind flux

                     norms = (nx * (u_c_val(i, 1, 2) / u_c_val(i, 1, 1))) &
                        & + (ny * (u_c_val(i, 1, 3) / u_c_val(i, 1, 1))) &
                        & + (nz * (u_c_val(i, 1, 4) / u_c_val(i, 1, 1)))
                     if (norms .ge. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbl(1:turbulenceequations + passivescalar)
                     end if
                     if (norms .lt. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbr(1:turbulenceequations + passivescalar)
                     end if

                  end if

                  godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                     & 1:nof_variables + turbulenceequations + passivescalar) + &
                     & (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                     & ielem_surf(i, l)))

               end if

            end do

            rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) + godflux2(1:nof_variables)
            if (multispecies .eq. 1) then
               mp_source3 = mp_source3 + mp_source2
            end if
!                                     if ((igoflux.eq.1))then
!                                     rhs_val(ielem_ineigh(i,l),1:nof_variables)=rhs_val(ielem_ineigh(i,l),1:nof_variables)-godflux2(1:nof_variables)
!                                     end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) + &
                  & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!                                      if ((igoflux.eq.1))then
!                                      rhst_val(ielem_ineigh(i,l),1:turbulenceequations+passivescalar)=rhst_val(ielem_ineigh(i,l),1:turbulenceequations+passivescalar)-&
!                                     godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
!                                     end if
            end if
!                                     end if
         end do
         if (multispecies .eq. 1) then
            call subtract_volume_fraction_source(i, mp_source3, .false.)

         end if
      end if
   end subroutine calculate_fluxeshi_convective_mood_bound_cell

   subroutine calculate_fluxeshi_convective2d_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
      implicit none
      integer, intent(in) :: n
      integer :: ii

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_interior) &
!$omp& map(alloc: el_int, gamma_in, ielem_faceanglex, ielem_faceangley, ielem_ifca, ielem_ineigh, ielem_ineighn) &
!$omp& map(alloc: ielem_mood, ielem_recalc, ielem_surf, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in, mp_tmid_in) &
!$omp& map(alloc: rec_uleft, rec_uleftturb, rg_hzero, rg_molm, rhs_val, rhs_valdg, rhst_val, u_c_val, u_ct_val) &
!$omp& map(alloc: weights_l)
#else
!$omp barrier
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_interior !for all the interior elements
         call calculate_fluxeshi_convective2d_mood_inner_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n) private(ii) &
!$omp& firstprivate(nof_bounded) &
!$omp& map(alloc: bleed_end, bleed_plenum, bleed_porosity, bleed_start, bound_offset, boundhir, boundhirm) &
!$omp& map(alloc: catalytic_con, el_bnd, gamma_in, halo_offset, ibound_icode, ielem_bleedn, ielem_faceanglex) &
!$omp& map(alloc: ielem_faceangley, ielem_ibounds, ielem_ifca, ielem_indexf, ielem_indexi, ielem_ineigh) &
!$omp& map(alloc: ielem_ineighb, ielem_ineighn, ielem_inter_id, ielem_mood, ielem_nodes_faces, ielem_qface) &
!$omp& map(alloc: ielem_recalc, ielem_surf, ielem_walldist, inoder4_cord, mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in) &
!$omp& map(alloc: mp_thigh_in, mp_tlow_in, mp_tmid_in, rec_ihexl, rec_ihexn, rec_local, rec_uleft, rec_uleftturb) &
!$omp& map(alloc: rg_hzero, rg_molm, rg_thetag, rg_vf, rgs_ab, rgs_bb, rgs_cb, rhs_val, rhs_valdg, rhst_val, solhir) &
!$omp& map(alloc: u_c_val, u_ct_val, weights_l)
#else
!$omp barrier
!$omp do &
!$omp& private(ii)
#endif
      do ii = 1, nof_bounded
         call calculate_fluxeshi_convective2d_mood_bound_cell(n, ii)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

   end subroutine calculate_fluxeshi_convective2d_mood

   subroutine calculate_fluxeshi_convective2d_mood_inner_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2, dg_vol_rec, rhllcflux, hllcflux
      integer :: i, l, ngp, kmaxe, iqp, ikas, igoflux, icaseb, kxk, b_code
      real :: sum_detect, norms
      real, dimension(1:gpu_max_qp_face) :: weights_temp
      integer :: iconsidered, facex, pointx, n_node
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
      real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:gpu_max_nvar) :: srf_speed
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      integer :: ibfc, nfx, lfx, rowfx, nf, lf, rowf

      i = el_int(ii)

      if (ielem_recalc(i) .eq. 1) then
         rhs_val(i, :) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst_val(i, :) = zero
      end if
      weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)

      iconsidered = i
!                     rhs_val(i,:)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst_val(i,:)=zero
      if (ielem_recalc(i) .eq. 1) then
         mp_source3 = zero
         do l = 1, ielem_ifca(i) !for all their faces

            godflux2 = zero
            mp_source2 = zero
            nx = ielem_faceanglex(i, l)
            ny = ielem_faceangley(i, l)
            angle1 = nx
            angle2 = ny
            b_code = 0
            iqp = qp_line_n

            do ngp = 1, iqp !for all the gaussian quadrature points
               cleft(1:nof_variables) = rec_uleft(i, 1:nof_variables, l, ngp) !left mean flow state

               if ((cascade .eq. 2) .and. (ielem_mood(i) .ge. 1)) then
                  cleft(1:nof_variables) = u_c_val(i, 3, 1:nof_variables)
               end if

               !right mean flow state
               cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)

               if ((cascade .eq. 2) .and. (ielem_mood(ielem_ineigh(i, l)) .ge. 1)) then
                  cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 3, 1:nof_variables)
               end if

!
               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then
                     !left additional equations flow state
                     cturbl(1:turbulenceequations + passivescalar) = rec_uleftturb(i, 1:turbulenceequations + passivescalar, l, ngp)
                     !right additional equations flow state
                     cturbr(1:turbulenceequations + passivescalar) = rec_uleftturb(ielem_ineigh(i, l), &
                        & 1:turbulenceequations + passivescalar, ielem_ineighn(i, l), ngp)
                  else
                     cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
                     cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                        & 1:turbulenceequations + passivescalar)
                  end if

                  cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations &
                     & + passivescalar)
                  cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations &
                     & + passivescalar)
               end if

               !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotatef2d(n, cright_rot, cright, angle1, angle2)
               call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem

               if ((lmach .eq. 1)) then !application of the low mach number correction
                  leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
                  call lmacht2d(n, leftv, rightv)
                  cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               end if

               select case (iriemann)

               case (1) !hllc

                  call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (2) !rusanov

                  call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if
               case (3) !roe

                  call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb2d(n, cright, cright_rot, angle1, angle2)
                  call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)

                  rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)

               case (9) !hll

                  call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (4) !roe

                  call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb2d(n, cright, cright_rot, angle1, angle2)
                  call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)

                  rhllcflux = hllcflux

               end select

               godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables) * (weights_temp(ngp) * &
                  & ielem_surf(i, l)))
               if (multispecies .eq. 1) then
                  mp_source2 = mp_source2 + mp_source1 * (weights_temp(ngp) * ielem_surf(i, l))
               end if
               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 0) then !first order upwind flux

                     norms = (nx * (u_c_val(i, 1, 2) / u_c_val(i, 1, 1))) &
                        & + (ny * (u_c_val(i, 1, 3) / u_c_val(i, 1, 1)))
                     if (norms .ge. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbl(1:turbulenceequations + passivescalar)
                     end if
                     if (norms .lt. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbr(1:turbulenceequations + passivescalar)
                     end if

                  end if

                  godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                     & 1:nof_variables + turbulenceequations + passivescalar) + &
                     & (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                     & ielem_surf(i, l)))

!

               end if

            end do
            rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) + godflux2(1:nof_variables)
            if (multispecies .eq. 1) then
               mp_source3 = mp_source3 + mp_source2
            end if

!                                      rhs_val(ielem_ineigh(i,l),1:nof_variables)=rhs_val(ielem_ineigh(i,l),1:nof_variables)-godflux2(1:nof_variables)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) + &
                  & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!                                     rhst_val(ielem_ineigh(i,l),1:turbulenceequations+passivescalar)=rhst_val(ielem_ineigh(i,l),1:turbulenceequations+passivescalar)-&
!                                      godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
!
            end if
!                                     end if
         end do
         if (multispecies .eq. 1) then
            call subtract_volume_fraction_source(i, mp_source3, .false.)

         end if

      end if
   end subroutine calculate_fluxeshi_convective2d_mood_inner_cell

   subroutine calculate_fluxeshi_convective2d_mood_bound_cell(n, ii)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: ii
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar_total) :: godflux2, dg_vol_rec, rhllcflux, hllcflux
      integer :: i, l, ngp, kmaxe, iqp, ikas, igoflux, icaseb, kxk, b_code
      real :: sum_detect, norms
      real, dimension(1:gpu_max_qp_face) :: weights_temp
      integer :: iconsidered, facex, pointx, n_node
      real :: angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
      real, dimension(1:gpu_max_nvar_total) :: cleft, cright, cleft_rot, cright_rot
      real, dimension(1:gpu_max_nvar) :: leftv, rightv, srf_speedrot
      real, dimension(1:gpu_max_extra_transport) :: cturbl, cturbr
      real, dimension(1:gpu_max_dim) :: pox, poy, poz
      real, dimension(1:gpu_max_nvar) :: srf_speed
      real, dimension(1:8, 1:gpu_max_dim) :: vext, nodes_list
      real, dimension(3) :: cords
      integer :: ibfc, nfx, lfx, rowfx, nf, lf, rowf

      i = el_bnd(ii)
      iconsidered = i

      if (ielem_recalc(i) .eq. 1) then
         rhs_val(i, :) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst_val(i, :) = zero
      end if
      weights_temp(1:qp_line_n) = weights_l(1:qp_line_n)

      mp_source3 = zero
      if (ielem_recalc(i) .eq. 1) then

         do l = 1, ielem_ifca(i)
            igoflux = 0
            if (ielem_ineighb(i, l) .eq. n) then !my cpu only
               if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                  if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                     icaseb = 1 !periodic mine
!
                  else
                     icaseb = 3 !physical
!
                  end if

               else
!
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

            b_code = 0
            nx = ielem_faceanglex(i, l)
            ny = ielem_faceangley(i, l)
            angle1 = nx
            angle2 = ny

            iqp = qp_line_n
            godflux2 = zero
            mp_source2 = zero
            do ngp = 1, iqp
               cleft(1:nof_variables) = rec_uleft(i, 1:nof_variables, l, ngp)

               if ((cascade .eq. 2) .and. (ielem_mood(i) .ge. 1)) then
                  cleft(1:nof_variables) = u_c_val(i, 3, 1:nof_variables)

               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then
                     cturbl(1:turbulenceequations + passivescalar) = rec_uleftturb(i, 1:turbulenceequations + passivescalar, l, ngp)
                  else
                     cturbl(1:turbulenceequations + passivescalar) = u_ct_val(i, 1, 1:turbulenceequations + passivescalar)
                  end if
               end if

               if (ielem_ineighb(i, l) .eq. n) then !my cpu only
                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in my cpu
                        cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)

                        if ((cascade .eq. 2) .and. (ielem_mood(ielem_ineigh(i, l)) .ge. 1)) then
                           cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 3, 1:nof_variables)

                        end if

                        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                           if (icoupleturb .eq. 1) then
                              !right additional equations flow state
                              cturbr(1:turbulenceequations + passivescalar) = rec_uleftturb(ielem_ineigh(i, l), &
                                 & 1:turbulenceequations + passivescalar, ielem_ineighn(i, l), ngp)
                           else
                              cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                                 & 1:turbulenceequations + passivescalar)
                           end if
                        end if

                        ikas = 1

                     else
                        !not periodic ones in my cpu

                        facex = l; iconsidered = i
                        call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                        cords(1:2) = zero
                        n_node = 2
                        call cordinates2(n, nodes_list, n_node, cords(1:2))

                        poy(1) = cords(2)
                        pox(1) = cords(1)

                        leftv(1:nof_variables) = cleft(1:nof_variables)
!
                        b_code = ibound_icode(ielem_ibounds(i, l))
                        call boundarys2d(n, b_code, iconsidered, facex, leftv, rightv, pox, poy, poz, angle1, angle2, nx, ny, nz, &
                           & cturbl, cturbr, cright_rot, cleft_rot, srf_speed, srf_speedrot, ibfc)
                        cright(1:nof_variables) = rightv(1:nof_variables)
!
                        ikas = 2

                     end if
                  else
                     cright(1:nof_variables) = rec_uleft(ielem_ineigh(i, l), 1:nof_variables, ielem_ineighn(i, l), ngp)

                     if ((cascade .eq. 2) .and. (ielem_mood(ielem_ineigh(i, l)) .ge. 1)) then
                        cright(1:nof_variables) = u_c_val(ielem_ineigh(i, l), 3, 1:nof_variables)

                     end if

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                        if (icoupleturb .eq. 1) then
                           !right additional equations flow state
                           cturbr(1:turbulenceequations + passivescalar) = rec_uleftturb(ielem_ineigh(i, l), &
                              & 1:turbulenceequations + passivescalar, ielem_ineighn(i, l), ngp)
                        else
                           cturbr(1:turbulenceequations + passivescalar) = u_ct_val(ielem_ineigh(i, l), 1, &
                              & 1:turbulenceequations + passivescalar)
                        end if
                     end if

                     ikas = 3

                  end if
               else !in other cpus they can only be periodic or mpi neighbours

                  if (ielem_ibounds(i, l) .gt. 0) then !check for boundaries
                     if (ibound_icode(ielem_ibounds(i, l)) .eq. 5) then !periodic in other cpu
!                                                                           cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                        nfx = ielem_ineighn(i, l)
                        lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                        rowfx = bound_offset(nfx) + lfx - 1

                        cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)

                        if ((cascade .eq. 2) .and. (boundhirm(rowfx) .ge. 1)) then
!                                     cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                         (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                           nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                           lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                           rowf = halo_offset(nf) + lf - 1
                           cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                        end if

                        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                           if (icoupleturb .eq. 1) then
!                                                                            cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                              nfx = ielem_ineighn(i, l)
                              lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                              rowfx = bound_offset(nfx) + lfx - 1
                              cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                                 & turbulenceequations + passivescalar)
                           else
!                                                                          cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                              nfx = ielem_ineighn(i, l)
                              lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                              rowfx = bound_offset(nfx) + lfx - 1
                              cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                                 & turbulenceequations + passivescalar)
                           end if

                        end if

                     end if
                  else

!                                                                   cright(1:nof_variables)=iexboundhir(ielem_ineighn(i,l))%facesol(ielem_qface(l,ngp,i),1:nof_variables)
                     nfx = ielem_ineighn(i, l)
                     lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                     rowfx = bound_offset(nfx) + lfx - 1
!                                                                   cright(1:nof_variables) = boundhir(rowfx,1:nof_variables)
                     cright(1:nof_variables) = boundhir(rowfx, 1:nof_variables)

                     if ((cascade .eq. 2) .and. (boundhirm(rowfx) .ge. 1)) then
!                                     cright(1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i)))%sol&
!                                         (rec_ihexl(1,ielem_indexi(l,i)),1:nof_variables)

                        nf = rec_ihexn(rec_local(i), 1, ielem_indexi(i, L))
                        lf = rec_ihexl(i, 1, ielem_indexi(i, L))
                        rowf = halo_offset(nf) + lf - 1
                        cright(1:nof_variables) = solhir(rowf, 1:nof_variables)

                     end if

                     if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                        if (icoupleturb .eq. 1) then
!                                                                            cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)
                        else
!                                                                          cturbr(1:turbulenceequations+passivescalar)=iexboundhir(ielem_ineighn(i,l))%facesol&
!                                                                            (ielem_qface(l,ngp,i),nof_variables+1:nof_variables+turbulenceequations+passivescalar)!right additional equations flow state
                           nfx = ielem_ineighn(i, l)
                           lfx = ielem_qface(ielem_inter_id(ielem_indexf(i)), l, ngp)
                           rowfx = bound_offset(nfx) + lfx - 1
                           cturbr(1:turbulenceequations + passivescalar) = boundhir(rowfx, nof_variables + 1:nof_variables + &
                              & turbulenceequations + passivescalar)
                        end if

                     end if
                     ikas = 4
!
                  end if
               end if

               !rotate wrt to normalvector of face and solve 1d riemann problem
               call rotatef2d(n, cright_rot, cright, angle1, angle2)
               call rotatef2d(n, cleft_rot, cleft, angle1, angle2) !rotate wrt to normalvector of face and solve 1d riemann problem
!

               if ((lmach .eq. 1)) then !application of the low mach number correction
                  leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)

                  call lmacht2d(n, leftv, rightv)
                  cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                  cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations &
                     & + passivescalar)
                  cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations &
                     & + passivescalar)
               end if

               select case (iriemann)

               case (1) !hllc

                  call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (2) !rusanov

                  call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (9) !hll

                  call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
                  call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = hllcflux(nof_variables + &
                        & 1:nof_variables + turbulenceequations + passivescalar)
                  end if

               case (3) !roe

                  call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                  call rotateb2d(n, cright, cright_rot, angle1, angle2)
                  call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)

                  rhllcflux = hllcflux

               case (4) !roe

                  if ((b_code .le. 0)) then

                     call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                     call rotateb2d(n, cright, cright_rot, angle1, angle2)

                     call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)
                     rhllcflux = hllcflux
                  else

!

                     call rusanov_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot)
                     call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)

                  end if

               end select

               godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables) * (weights_temp(ngp) * &
                  & ielem_surf(i, l)))
               if (multispecies .eq. 1) then
                  mp_source2 = mp_source2 + mp_source1 * (weights_temp(ngp) * ielem_surf(i, l))
               end if

               if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 0) then !first order upwind flux

                     norms = (nx * (u_c_val(i, 1, 2) / u_c_val(i, 1, 1))) &
                        & + (ny * (u_c_val(i, 1, 3) / u_c_val(i, 1, 1)))

                     if (norms .ge. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbl(1:turbulenceequations + passivescalar)
                     end if
                     if (norms .lt. zero) then
                        rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = (norms) * &
                           & cturbr(1:turbulenceequations + passivescalar)
                     end if

                  end if

                  if ((b_code .eq. 4) .or. (b_code .eq. 3)) then
                     rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
                  end if

                  godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = godflux2(nof_variables + &
                     & 1:nof_variables + turbulenceequations + passivescalar) + &
                     & (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) * (weights_temp(ngp) * &
                     & ielem_surf(i, l)))
!

               end if

            end do

            rhs_val(i, 1:nof_variables) = rhs_val(i, 1:nof_variables) + godflux2(1:nof_variables)
            if (multispecies .eq. 1) then
               mp_source3 = mp_source3 + mp_source2
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhst_val(i, 1:turbulenceequations + passivescalar) = rhst_val(i, 1:turbulenceequations + passivescalar) + &
                  & godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)

            end if
         end do
         if (multispecies .eq. 1) then
            call subtract_volume_fraction_source(i, mp_source3, .false.)

         end if

      end if
   end subroutine calculate_fluxeshi_convective2d_mood_bound_cell

! !experimental subroutine---not for production
! subroutine vertex_neighbours_values(n)
! implicit none
! integer,intent(in)::n
! integer::i,j,k,l,m,nj
! real,dimension(4,30,nof_variables)::vertex_neighbours_vals        !maximum allocation not efficient
!
! i=iconsidered
!
! write(680+n,*),"element number",ielem_ihexgl(i)
!
! do nj=1,ielem_nonodes(i)
!
!         if (rec_local(i).eq.1)then
!                 do j=1,ielem_nojecount(nj,i)
!                         vertex_neighbours_vals(nj,j,1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,ielem_nodes_neighbours(nj,j,i,i)))
!                 !rec_ihexl(1,l,i)        !local numbering in my cpu
!                 !rec_ihexg(1,l,i)        !global numbering in my cpu
!                         write(680+n,*),"node",nj,"node neighbour",j,"local number",ielem_nodes_neighbours(nj,j,i), "values",vertex_neighbours_vals(nj,j,1:nof_variables)
!
!                 end do
!
!         else
!                         do j=1,ielem_nojecount(nj,i)
!
!                                 if (rec_ihexb(1,ielem_nodes_neighbours(nj,j,i,i)).eq.n)then
!                                         vertex_neighbours_vals(nj,j,1:nof_variables)=u_c_val(1,1:nof_variables,rec_ihexl(1,ielem_nodes_neighbours(nj,j,i,i)))
!                                 else
!                                         write(500+n,*)rec_ihexn(1,ielem_nodes_neighbours(nj,j,i,i)),ielem_nodes_neighbours(nj,j,i)
!                                         vertex_neighbours_vals(nj,j,1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_nodes_neighbours(nj,j,i,i)))%sol(rec_ihexl(1,ielem_nodes_neighbours(nj,j,i,i)),1:nof_variables)
!
!                                 end if
!                         write(680+n,*),"node",nj,"node neighbour",j,"local number",ielem_nodes_neighbours(nj,j,i), "values",vertex_neighbours_vals(nj,j,1:nof_variables)
!                         end do
!         end if
!
! end do
!
! !now you have collected for every vertex the values
!
!
! end subroutine vertex_neighbours_values

! Bridges for the general host fluxes; ideal-gas XPU kernels return before them.


end module fluxes
