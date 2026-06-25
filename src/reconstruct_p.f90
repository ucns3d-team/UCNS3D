module recon
use declaration
use derivatives
use library
use transform
use local
use lapck
use gradients
use basis
implicit none

 contains


subroutine clear_rec_uleft_for_weno(n)
implicit none
integer,intent(in)::n
integer::i,l,ngp,iqp,iex,kmaxe,nextra

kmaxe=xmpielrank(n)

#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe,dimensiona,qp_quad,qp_triangle,qp_line,nof_variables,zero) &
!$omp& map(alloc: rec_uleft, ielem_ifca, ielem_types_faces) &
!$omp& private(i,l,ngp,iqp,iex)
#else
!$omp do private(i,l,ngp,iqp,iex)
#endif
do i=1,kmaxe
  do l=1,ielem_ifca(i)
    if (dimensiona.eq.3) then
      if (ielem_types_faces(l,i).eq.5) then
        iqp=qp_quad
      else
        iqp=qp_triangle
      end if
    else
      iqp=qp_line
    end if
    do ngp=1,iqp
      do iex=1,nof_variables
        rec_uleft(iex,l,ngp,i)=zero
      end do
    end do
  end do
end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

nextra=turbulenceequations+passivescalar
if ((nextra.gt.0).and.allocated(rec_uleftturb)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe,nextra) &
!$omp& firstprivate(dimensiona, qp_quad, qp_triangle, qp_line, zero) &
!$omp& map(alloc: rec_uleftturb, ielem_ifca, ielem_types_faces) &
!$omp& private(i,l,ngp,iqp,iex)
#else
!$omp do private(i,l,ngp,iqp,iex)
#endif
  do i=1,kmaxe
    do l=1,ielem_ifca(i)
      if (dimensiona.eq.3) then
        if (ielem_types_faces(l,i).eq.5) then
          iqp=qp_quad
        else
          iqp=qp_triangle
        end if
      else
        iqp=qp_line
      end if
      do ngp=1,iqp
        do iex=1,nextra
          rec_uleftturb(iex,l,ngp,i)=zero
        end do
      end do
    end do
  end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if

end subroutine clear_rec_uleft_for_weno








subroutine average_stresses(n)
implicit none
!> @brief
!> subroutine for calling the computation of the average shear stresses
integer,intent(in)::n
integer::ii,i,iconsidered
integer::number_of_dog,number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered,i,number_of_dog,number_of_nei)
#else
!$omp do private(iconsidered,i,number_of_dog,number_of_nei)
#endif
do ii=1,nof_interior
      i=el_int(ii)
      iconsidered=i
      number_of_dog=ielem_idegfree(i)
      number_of_nei=ielem_inumneighbours(i)
      call compute_gradients_inner_mean_ggs_viscous_av(n,iconsidered,number_of_dog,number_of_nei)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered,i,number_of_dog,number_of_nei)
#else
!$omp do private(iconsidered,i,number_of_dog,number_of_nei)
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)
	iconsidered=i
	number_of_dog=ielem_idegfree(i)
	number_of_nei=ielem_inumneighbours(i)
	call compute_gradients_mix_mean_ggs_viscous_av(n,iconsidered,number_of_dog,number_of_nei)
end do	
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine average_stresses
	
subroutine memory_fast(n)
  !> @brief
  !> subroutine for storing the gaussian quadrature points at the cell interfaces
  implicit none
integer :: mm_i
integer :: mm_j
real :: mm_old(1:3)
  integer, intent(in) :: n
  integer :: i, k, kmaxe, idummy, l, nnd, iqp, ngp, iex
  integer :: iconsidered, facex, pointx
  real, dimension(1:dimensiona,1:numberofpoints2) :: qpoints2d
  real, dimension(1:numberofpoints2) :: wequa2d
  real, dimension(1:8,1:dimensiona) :: vext, nodes_list
  real, dimension(1:dimensiona) :: pox, poy, poz
  real, dimension(3) :: rotvec

  kmaxe = xmpielrank(n)

  allocate(rec_qpoints(max_faces,numberofpoints2,1:dimensiona,1:kmaxe)); rec_qpoints = zero

!   if (initcond.gt.100000)allocate(rec_qpoints_p(max_faces,numberofpoints2,1:dimensiona,1:kmaxe)); rec_qpoints = zero

  allocate(rec_mrf(1:kmaxe)); rec_mrf = 0
  if (srfg == 1) then
    allocate(rec_rpoints(max_faces,numberofpoints2,1:dimensiona,1:kmaxe)); rec_rpoints = zero
    allocate(rec_rotvel(max_faces,numberofpoints2,1:dimensiona,1:kmaxe));  rec_rotvel  = zero
  end if
  if (mrf == 1) then
    ! keep both; your original code allocated them under mrf too
    if (.not. allocated(rec_rpoints)) then
      allocate(rec_rpoints(max_faces,numberofpoints2,1:dimensiona,1:kmaxe)); rec_rpoints = zero
    end if
    if (.not. allocated(rec_rotvel)) then
      allocate(rec_rotvel(max_faces,numberofpoints2,1:dimensiona,1:kmaxe));  rec_rotvel  = zero
    end if
    allocate(rec_mrf_origin(1:3,1:kmaxe));   rec_mrf_origin   = zero
    allocate(rec_mrf_velocity(1:3,1:kmaxe)); rec_mrf_velocity = zero
  end if

  if (dimensiona == 3) then

    do i = 1, kmaxe
      iconsidered = i

      do l = 1, ielem_ifca(i)
        idummy = 0

        if ((iperiodicity == 1) .and. (ielem_interior(i) == 1)) then
          if (ielem_ibounds(l,i) > 0) then
            if ((ibound_icode(ielem_ibounds(l,i)) == 5) .or. (ibound_icode(ielem_ibounds(l,i)) == 50)) then
              idummy = 1
            end if
          end if

          if (ielem_types_faces(l,i) == 5) then
            iqp = qp_quad
            nnd = 4
            if (idummy == 0) then
              do k = 1, nnd
                vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%cord(1:dims)
              end do
            else
              facex = l
              call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
            end if
            call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
          else
            iqp = qp_triangle
            nnd = 3
            if (idummy == 0) then
              do k = 1, nnd
                vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
              end do
            else
              facex = l
              call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
            end if
            call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
          end if

        else
          if (ielem_types_faces(l,i) == 5) then
            iqp = qp_quad
            nnd = 4
            do k = 1, nnd
              vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
            end do
            call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
          else
            iqp = qp_triangle
            nnd = 3
            do k = 1, nnd
              vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
            end do
            call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
          end if
        end if

        do ngp = 1, iqp
          if (srfg == 1) then
            rec_rpoints(l,ngp,1:3,i) = qpoints2d(1:3,ngp)
            pox(1:3) = rec_rpoints(l,ngp,1:3,i) - srf_origin(1:3)
            poy(1:3) = srf_velocity(1:3)
            call vect_function(pox,poy,rotvec)
            rec_rotvel(l,ngp,1,i) = rotvec(1)
            rec_rotvel(l,ngp,2,i) = rotvec(2)
            rec_rotvel(l,ngp,3,i) = rotvec(3)
          end if

          if (mrf == 1) then
            rec_rpoints(l,ngp,1:3,i) = qpoints2d(1:3,ngp)
            pox(1) = ielem_xxc(i); pox(2) = ielem_yyc(i); pox(3) = ielem_zzc(i)
            poy(1:3) = rec_rpoints(l,ngp,1:3,i)
            facex  = l
            pointx = ngp
            call mrfswitch(n, iconsidered, facex, pointx, pox, poy)
          end if
        end do
      end do

      ! --- store qpoints in reference space ---
      do l = 1, ielem_ifca(i)
        idummy = 0

        if ((iperiodicity == 1) .and. (ielem_interior(i) == 1)) then
          if (ielem_ibounds(l,i) > 0) then
            if ((ibound_icode(ielem_ibounds(l,i)) == 5) .or. (ibound_icode(ielem_ibounds(l,i)) == 50)) then
              idummy = 1
            end if
          end if

          if (ielem_types_faces(l,i) == 5) then
            iqp = qp_quad
            nnd = 4
            if (idummy == 0) then
              do k = 1, nnd
                vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
                  do mm_j=1,3
                    mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                  end do
                  do mm_i=1,3
                    vext(k,mm_i)=zero
                    do mm_j=1,3
                      vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                    end do
                  end do
              end do
            else
              facex = l
              call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
              do k = 1, nnd
                  do mm_j=1,3
                    mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                  end do
                  do mm_i=1,3
                    vext(k,mm_i)=zero
                    do mm_j=1,3
                      vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                    end do
                  end do
              end do
            end if
            call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
          else
            iqp = qp_triangle
            nnd = 3
            if (idummy == 0) then
              do k = 1, nnd
                vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
                  do mm_j=1,3
                    mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                  end do
                  do mm_i=1,3
                    vext(k,mm_i)=zero
                    do mm_j=1,3
                      vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                    end do
                  end do
              end do
            else
              facex = l
              call coordinates_face_period1(n, iconsidered, facex, vext, nodes_list)
              do k = 1, nnd
                  do mm_j=1,3
                    mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                  end do
                  do mm_i=1,3
                    vext(k,mm_i)=zero
                    do mm_j=1,3
                      vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                    end do
                  end do
              end do
            end if
            call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
          end if

        else
          if (ielem_types_faces(l,i) == 5) then
            iqp = qp_quad
            nnd = 4
            do k = 1, nnd
              vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
                do mm_j=1,3
                  mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                end do
                do mm_i=1,3
                  vext(k,mm_i)=zero
                  do mm_j=1,3
                    vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                  end do
                end do
            end do
            call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
          else
            iqp = qp_triangle
            nnd = 3
            do k = 1, nnd
              vext(k,1:3) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
                do mm_j=1,3
                  mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                end do
                do mm_i=1,3
                  vext(k,mm_i)=zero
                  do mm_j=1,3
                    vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                  end do
                end do
            end do
            call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
          end if
        end if

        do ngp = 1, iqp
          rec_qpoints(l,ngp,1:3,i) = qpoints2d(1:3,ngp)
        end do
      end do
    end do

  else
    ! ------------------- 2d -------------------
    do i = 1, kmaxe
      iconsidered = i


      do l = 1, ielem_ifca(i)



        idummy = 0

        if ((iperiodicity == 1) .and. (ielem_interior(i) == 1)) then
          if (ielem_ibounds(l,i) > 0) then
            if ((ibound_icode(ielem_ibounds(l,i)) == 5) .or. (ibound_icode(ielem_ibounds(l,i)) == 50)) then
              idummy = 1
            end if
          end if

          iqp = qp_line
          nnd = 2
          if (idummy == 0) then
                do k = 1, nnd
                  vext(k,1:2) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
                    do mm_j=1,2
                      mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                    end do
                    do mm_i=1,2
                      vext(k,mm_i)=zero
                      do mm_j=1,2
                        vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                      end do
                    end do
                end do
          else
                facex = l
                call coordinates_face_period2d1(n, iconsidered, facex, vext, nodes_list)
                    do k = 1, nnd
                        do mm_j=1,2
                          mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                        end do
                        do mm_i=1,2
                          vext(k,mm_i)=zero
                          do mm_j=1,2
                            vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                          end do
                        end do
                    end do
          end if
          call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)

        else
          iqp = qp_line
          nnd = 2

!               if (initcond.gt.100000)then
!
!
!               do k = 1, nnd
!                 vext(k,1:2) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
! !                 vext(k,1:2) = explicit_matrix_product(rec_invccjac(:,:,i), vext(k,1:2) - rec_vext_ref(1:2,i))
!               end do
!               call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
!
!
!               do ngp = 1, iqp
!               rec_qpoints_p(l,ngp,1:2,i) = qpoints2d(1:2,ngp)
!             end do





              do k = 1, nnd
                vext(k,1:2) = dinoder(ielem_nodes_faces(l,k,i))%CORD(1:dims)
                  do mm_j=1,2
                    mm_old(mm_j)=vext(k,mm_j)-rec_vext_ref(mm_j,i)
                  end do
                  do mm_i=1,2
                    vext(k,mm_i)=zero
                    do mm_j=1,2
                      vext(k,mm_i)=vext(k,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
                    end do
                  end do
              end do
              call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)








        end if

        do ngp = 1, iqp
          rec_qpoints(l,ngp,1:2,i) = qpoints2d(1:2,ngp)
        end do
      end do
    end do
  end if

end subroutine memory_fast





subroutine extrapolate_bound_linear(n,du,facex,pointx,iconsidered)
  implicit none
!> @brief
!> Pointwise extrapolation for the linear reconstruction using a single increment vector.
#ifdef gpu
!$omp declare target
#endif
  integer,intent(in) :: n,facex,pointx,iconsidered
  real,intent(in)    :: du(1:nof_variables+turbulenceequations+passivescalar)   ! size = NVTOT
  real :: mp_pinfl,gammal
  real,dimension(1:gpu_max_nvar) :: leftv

  if (wenwrt.eq.3) then
    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,iconsidered)
    call cons2prim(n,leftv,mp_pinfl,gammal)
    leftv(1:nof_variables)=leftv(1:nof_variables)+du(1:nof_variables)
    call prim2cons(n,leftv)
    rec_uleft(1:nof_variables,facex,pointx,iconsidered)=rec_uleft(1:nof_variables,facex,pointx,iconsidered)+leftv(1:nof_variables)
  else
    rec_uleft(1:nof_variables,facex,pointx,iconsidered)=(u_c_val(1,1:nof_variables,iconsidered) + du(1:nof_variables))
  end if

  if (turbulenceequations.ge.1) then
    rec_uleftturb(1:turbulenceequations+passivescalar,facex,pointx,iconsidered) = &
      (u_ct_val(1,1:turbulenceequations+passivescalar,iconsidered) + &
       du(nof_variables+1:nof_variables+turbulenceequations+passivescalar))
  end if

end subroutine extrapolate_bound_linear




subroutine extrapolate_bound_muscl(du,facex,pointx,iconsidered,slope)
  implicit none
!> @brief
!> Pointwise MUSCL extrapolation using a single increment vector (no usol/psi storage).
#ifdef gpu
!$omp declare target
#endif
  integer,intent(in) :: facex,pointx,iconsidered
  real,intent(in)    :: du(1:nof_variables+turbulenceequations+passivescalar)      ! size = nvtot (cons + turb + passive)
  real,intent(in)    :: slope(1:nof_variables+turbulenceequations+passivescalar)   ! same size as du
  real,dimension(1:gpu_max_nvar) :: leftv
  real :: mp_pinfl,gammal

  if (wenwrt.eq.3) then
    leftv(1:nof_variables) = u_c_val(1,1:nof_variables,iconsidered)
    call cons2prim(n,leftv,mp_pinfl,gammal)
    leftv(1:nof_variables) = leftv(1:nof_variables) + du(1:nof_variables) * slope(1:nof_variables)
    call prim2cons(n,leftv)
    rec_uleft(1:nof_variables,facex,pointx,iconsidered) =  leftv(1:nof_variables)
  else
    rec_uleft(1:nof_variables,facex,pointx,iconsidered) = &
      (u_c_val(1,1:nof_variables,iconsidered) + du(1:nof_variables) * slope(1:nof_variables))
  end if

  if (turbulenceequations.ge.1) then
    rec_uleftturb(1:turbulenceequations+passivescalar,facex,pointx,iconsidered) = &
      (u_ct_val(1,1:turbulenceequations+passivescalar,iconsidered) + &
       du(nof_variables+1:nof_variables+turbulenceequations+passivescalar) * &
       slope(nof_variables+1:nof_variables+turbulenceequations+passivescalar))
  end if

end subroutine extrapolate_bound_muscl















subroutine cp_reconstruction_cweno(iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif

integer :: kbasis
  integer, intent(in)  :: iconsidered
  real    :: divbyzero
  integer :: power

  integer :: i, l, ngp, iqp, ll, ll2, k, j, iex
  integer :: iadmis, n_faces, ideg_local, icompwrt
  integer :: ideg_max
  real    :: lwcx1, ax, ay, az
  real    :: tau_weno, sumomega
  real    :: mp_pinfl,gammal
  real    :: lamc(1:gpu_max_typesten), lambdaal(1:gpu_max_typesten), omegaatilde(1:gpu_max_typesten), omegaal(1:gpu_max_typesten)
  real    :: smooth(1:gpu_max_typesten)
  real    :: resvec(1:gpu_max_nvar)
  real    :: leftv(1:gpu_max_nvar)
  real    :: inv_lamc1

  real :: phi(1:gpu_max_dof)

  ! grad0 scratch
  real :: grad0(1:gpu_max_dof)

  ! WENO weights needed later
  real :: weno(1:gpu_max_nvar, 1:gpu_max_typesten)

  ! -------------------------
  ! Setup
  ! -------------------------
  i          = iconsidered
  iadmis     = ielem_admis(i)
  n_faces    = ielem_ifca(i)
  ideg_local = idegfree
  lwcx1      = ielem_linc(i)
  divbyzero  = 1.0e-6
  power      = 4

  ! Linear weights (lamc) depend only on element, not variable
  lamc(1) = (1.0d0 - (1.0d0 / lwcx1))
  if (iadmis > 1) then
    lamc(2:iadmis) = (1.0d0 - lamc(1)) / real(iadmis-1)
  end if
  inv_lamc1 = 1.0d0 / lamc(1)




  ! -------------------------
  ! Compute WENO weights per variable
  ! -------------------------
  do iex = 1, nof_variables

    smooth(1:iadmis)      = 0.0d0
    omegaatilde(1:iadmis) = 0.0d0
    omegaal(1:iadmis)     = 0.0d0

    ! ---- build grad0(:)  ----
    ! grad0 = inv_lamc1 * ( rec_gradients(1) - sum_{m=2..iadmis} lamc(m)*rec_gradientsc(m) )
    grad0(1:ideg_local) = rec_gradients(1, 1:ideg_local, iex, i)
    do ll = 2, iadmis
      do k = 1, idegfree2
        grad0(k) = grad0(k) - lamc(ll) * rec_gradientsc(ll, k, iex, i)
      end do
    end do
    grad0(1:ideg_local) = inv_lamc1 * grad0(1:ideg_local)

    ! ---- smoothness indicator for ll=1  ----
    smooth(1) = 0.0d0
    do j = 1, ideg_local
      do k = 1, ideg_local
        smooth(1) = smooth(1) + grad0(j) * rec_indicator(j,k,i) * grad0(k)
      end do
    end do

    ! ---- smoothness indicators for ll>=2  ----
    do ll = 2, iadmis
      smooth(ll) = 0.0d0
      do j = 1, idegfree2
        do k = 1, idegfree2
          smooth(ll) = smooth(ll) + rec_gradientsc(ll,j,iex,i) * rec_indicatorc(j,k,i) * rec_gradientsc(ll,k,iex,i)
        end do
      end do
    end do

    ! ---- lambda ----
    lambdaal(1:iadmis) = lamc(1:iadmis)

    ! ---- omega tilde ----
    if (wenoz .eq. 1) then
      tau_weno = 0.0d0
      do ll = 1, iadmis
        tau_weno = tau_weno + abs(smooth(1) - smooth(ll))
      end do
      if (iadmis > 1) tau_weno = tau_weno / real(iadmis-1)

      do ll = 1, iadmis
        omegaatilde(ll) = lambdaal(ll) * (1.0d0 + (tau_weno / (divbyzero + smooth(ll)))**4)
      end do
    else
      do ll = 1, iadmis
        omegaatilde(ll) = lambdaal(ll) / ((divbyzero + smooth(ll))**4)
      end do
    end if

    ! ---- normalize ----
    sumomega = 0.0d0
    do ll = 1, iadmis
      sumomega = sumomega + omegaatilde(ll)
    end do

    do ll = 1, iadmis
      omegaal(ll)   = omegaatilde(ll) / sumomega
      weno(iex,ll)  = omegaal(ll)
    end do

     if (iex .eq. 1)ielem_wcx(i) = weno(1,1)

  end do

  ! -------------------------
  ! Reconstruction
  ! -------------------------
  do l = 1, n_faces
    if (dimensiona .eq. 3) then
      if (ielem_types_faces(l,i) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if
    do ngp = 1, iqp
      do iex = 1, nof_variables
        rec_uleft(iex,l,ngp,iconsidered) = 0.0d0
      end do
    end do
  end do
   if (dg .eq. 1) dg2fv(1:ideg_local,:,i) = 0.0d0

  do ll = 1, iadmis

    do l = 1, n_faces

      if (dimensiona .eq. 3) then
        if (ielem_types_faces(l,i) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if

      do ngp = 1, iqp

        ax = rec_qpoints(l,ngp,1,i)
        ay = rec_qpoints(l,ngp,2,i)
         if (dimensiona .eq. 3)az = rec_qpoints(l,ngp,3,i)

        resvec(:) = 0.0d0

        if (ll .ne. 1) then
          ! ---- low-order basis (degree idegfree2) -> phi(1:idegfree2) ----

          if (dimensiona .eq. 3) then
              do kbasis=1,idegfree2
                phi(kbasis) = basis_rec_value(n,ax,ay,az,iorder2,i,idegfree2,1,kbasis)
              end do
          else
              do kbasis=1,idegfree2
                phi(kbasis) = basis_rec2d_value(n,ax,ay,iorder2,i,idegfree2,1,kbasis)
              end do
          end if

          do k = 1, idegfree2
            do iex = 1, nof_variables
              resvec(iex) = resvec(iex) + phi(k) * rec_gradientsc(ll, k, iex, i)
            end do
          end do

        else
          ! ---- full-order basis (degree ideg_local) -> phi(1:ideg_local) ----

          if (dimensiona .eq. 3) then
              do kbasis=1,ideg_local
                phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg_local,0,kbasis)
              end do
          else
              do kbasis=1,ideg_local
                phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg_local,0,kbasis)
              end do
          end if

          ! resvec(iex) = inv_lamc1 * ( dot(phi, rec_gradients(1,:,iex)) - sum_{m=2..iadmis} lamc(m)*dot(phi, rec_gradientsc(m,:,iex)) )
          do iex = 1, nof_variables
            ! dot(phi, rec_gradients(1))
            do k = 1, ideg_local
              resvec(iex) = resvec(iex) + phi(k) * rec_gradients(1, k, iex, i)
            end do

            ! subtract lamc(m)*dot(phi, rec_gradientsc(m))
            do ll2 = 2, iadmis
              do k = 1, idegfree2
                resvec(iex) = resvec(iex) - lamc(ll2) * phi(k) * rec_gradientsc(ll2, k, iex, i)
              end do
            end do

            resvec(iex) = inv_lamc1 * resvec(iex)
          end do

        end if

        if (wenwrt.eq.3) then
          leftv(1:nof_variables) = u_c_val(1,1:nof_variables,iconsidered)
          call cons2prim(n,leftv,mp_pinfl,gammal)
          rec_uleft(1:nof_variables,l,ngp,iconsidered) = &
            rec_uleft(1:nof_variables,l,ngp,iconsidered) + &
            ((leftv(1:nof_variables)+resvec(1:nof_variables)) * weno(1:nof_variables,ll))
        else
          rec_uleft(1:nof_variables,l,ngp,iconsidered) = &
            rec_uleft(1:nof_variables,l,ngp,iconsidered) + &
            (u_c_val(1,1:nof_variables,iconsidered) + resvec(1:nof_variables)) * weno(1:nof_variables,ll)
        end if

      end do
    end do

    ! -------------------------
    ! DG accumulation (on-the-fly; ll==1 uses grad0 scratch per variable)
    ! -------------------------
    if (dg .eq. 1) then
      if (ll .eq. 1) then
        do iex = 1, nof_variables
          grad0(1:ideg_local) = rec_gradients(1, 1:ideg_local, iex, i)
          do ll2 = 2, iadmis
            do k = 1, idegfree2
              grad0(k) = grad0(k) - lamc(ll2) * rec_gradientsc(ll2, k, iex, i)
            end do
          end do
          grad0(1:ideg_local) = inv_lamc1 * grad0(1:ideg_local)

          dg2fv(1:ideg_local,iex,i) = dg2fv(1:ideg_local,iex,i) + grad0(1:ideg_local) * weno(iex,1)
        end do
      else
        do iex = 1, nof_variables
          dg2fv(1:idegfree2,iex,i) = dg2fv(1:idegfree2,iex,i) + rec_gradientsc(ll,1:idegfree2,iex,i) * weno(iex,ll)
        end do
      end if
    end if

  end do

  ! -------------------------
  ! Convert to conservative if requested
  ! -------------------------
  if (wenwrt .eq. 3) then
    do l = 1, n_faces
      if (dimensiona .eq. 3) then
        if (ielem_types_faces(l,i) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if

      do ngp = 1, iqp
        leftv(1:nof_variables) = rec_uleft(1:nof_variables,l,ngp,i)
        call prim2cons(n,leftv)
        rec_uleft(1:nof_variables,l,ngp,i) = leftv(1:nof_variables)
      end do
    end do
  end if

#ifndef gpu

#endif


end subroutine cp_reconstruction_cweno










subroutine cp_reconstruction_cweno_acc25(n)
  implicit none
  integer,intent(in)::n
  integer :: i,kmaxe,l,ngp,iqp,ll,ll2,k,j,iex
  integer :: iadmis,n_faces,ideg_local,degree,zpow,ypow,xpow,kk,ktail,ip
  integer :: dbg_i,dbg_l,dbg_ngp,dbg_iqp,dbg_unit
  integer :: dbg_k,dbg_j,dbg_iex,dbg_ll
  integer :: dbg_degree,dbg_zpow,dbg_ypow,dbg_xpow,dbg_ip,dbg_kk,dbg_ktail
  real :: lwcx1,divbyzero,inv_lamc1,tau_weno,sumomega
  real :: denom,denom2,denom4,ratio,ratio2,ratio4
  real :: ax,ay,az,volinv
  real :: raw_basis,phik,xterm,yterm,zterm
  real :: lamc(1:gpu_max_typesten)
  real :: lambdaal(1:gpu_max_typesten)
  real :: omegaatilde(1:gpu_max_typesten)
  real :: omegaal(1:gpu_max_typesten)
  real :: smooth(1:gpu_max_typesten)
  real :: weno(1:gpu_max_nvar,1:gpu_max_typesten)
  real :: grad0(1:gpu_max_dof)
  real :: resvec(1:gpu_max_nvar)


  kmaxe=xmpielrank(n)
#ifdef xpu

  !$omp target teams distribute parallel do &
  !$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,dg,idegfree,idegfree2) &
  !$omp& firstprivate(iorder2,qp_quad,qp_triangle,qp_line,lwci1,wenoz,zero) &
  !$omp& map(alloc: xmpielrank, ielem_admis, ielem_ifca, ielem_idegfree, ielem_linc) &
  !$omp& map(alloc: rec_gradients, rec_gradientsc, rec_indicator, rec_indicatorc, ielem_wcx) &
  !$omp& map(alloc: ielem_types_faces, rec_uleft, dg2fv, rec_qpoints, u_c_val) &
  !$omp& map(alloc: ielem_iorder, rec_volume, integ_basis_value, integ_basis_valuec) &
  !$omp& private(i,l,ngp,iqp,ll,ll2,k,j,iex,iadmis,n_faces,ideg_local) &
  !$omp& private(degree,zpow,ypow,xpow,kk,ktail,ip,lwcx1,divbyzero,inv_lamc1) &
  !$omp& private(tau_weno,sumomega,denom,denom2,denom4,ratio,ratio2,ratio4) &
  !$omp& private(ax,ay,az,volinv,raw_basis,phik,xterm,yterm,zterm) &
  !$omp& private(lamc,lambdaal,omegaatilde,omegaal,smooth,weno,grad0,resvec)
#else
  !$omp do private(i,l,ngp,iqp,ll,ll2,k,j,iex,iadmis,n_faces,ideg_local) &
  !$omp& private(degree,zpow,ypow,xpow,kk,ktail,ip,lwcx1,divbyzero,inv_lamc1) &
  !$omp& private(tau_weno,sumomega,denom,denom2,denom4,ratio,ratio2,ratio4) &
  !$omp& private(ax,ay,az,volinv,raw_basis,phik,xterm,yterm,zterm) &
  !$omp& private(lamc,lambdaal,omegaatilde,omegaal,smooth,weno,grad0,resvec)
#endif
  do i=1,kmaxe
    iadmis=ielem_admis(i)
    n_faces=ielem_ifca(i)
    ideg_local=ielem_idegfree(i)
    lwcx1=lwci1
    ielem_linc(i)=lwcx1
    divbyzero=1.0d-6
    volinv=1.0d0/rec_volume(1,1,i)

    lamc(1)=1.0d0-(1.0d0/lwcx1)
    do ll=2,gpu_max_typesten
      lamc(ll)=zero
    end do
    if (iadmis.gt.1) then
      do ll=2,iadmis
        lamc(ll)=(1.0d0-lamc(1))/real(iadmis-1)
      end do
    end if
    inv_lamc1=1.0d0/lamc(1)

    do ll=1,gpu_max_typesten
      lambdaal(ll)=zero
      do iex=1,nof_variables
        weno(iex,ll)=zero
      end do
    end do
    do ll=1,iadmis
      lambdaal(ll)=lamc(ll)
    end do

    do iex=1,nof_variables
      do ll=1,gpu_max_typesten
        smooth(ll)=zero
        omegaatilde(ll)=zero
        omegaal(ll)=zero
      end do

      do k=1,ideg_local
        grad0(k)=rec_gradients(1,k,iex,i)
      end do
      do ll=2,iadmis
        do k=1,idegfree2
          grad0(k)=grad0(k)-lamc(ll)*rec_gradientsc(ll,k,iex,i)
        end do
      end do
      do k=1,ideg_local
        grad0(k)=grad0(k)*inv_lamc1
      end do

      smooth(1)=zero
      do j=1,ideg_local
        do k=1,ideg_local
          smooth(1)=smooth(1)+grad0(j)*rec_indicator(j,k,i)*grad0(k)
        end do
      end do

      do ll=2,iadmis
        smooth(ll)=zero
        do j=1,idegfree2
          do k=1,idegfree2
            smooth(ll)=smooth(ll)+rec_gradientsc(ll,j,iex,i)*rec_indicatorc(j,k,i)*rec_gradientsc(ll,k,iex,i)
          end do
        end do
      end do

      if (wenoz.eq.1) then
        tau_weno=zero
        do ll=1,iadmis
          tau_weno=tau_weno+abs(smooth(1)-smooth(ll))
        end do
        if (iadmis.gt.1) tau_weno=tau_weno/real(iadmis-1)
        do ll=1,iadmis
          denom=divbyzero+smooth(ll)
          ratio=tau_weno/denom
          ratio2=ratio*ratio
          ratio4=ratio2*ratio2
          omegaatilde(ll)=lambdaal(ll)*(1.0d0+ratio4)
        end do
      else
        do ll=1,iadmis
          denom=divbyzero+smooth(ll)
          denom2=denom*denom
          denom4=denom2*denom2
          omegaatilde(ll)=lambdaal(ll)/denom4
        end do
      end if

      sumomega=zero
      do ll=1,iadmis
        sumomega=sumomega+omegaatilde(ll)
      end do
      do ll=1,iadmis
        omegaal(ll)=omegaatilde(ll)/sumomega
        weno(iex,ll)=omegaal(ll)
      end do
      if (iex.eq.1) ielem_wcx(i)=weno(1,1)
    end do

    do l=1,n_faces
      if (dimensiona.eq.3) then
        if (ielem_types_faces(l,i).eq.5) then
          iqp=qp_quad
        else
          iqp=qp_triangle
        end if
      else
        iqp=qp_line
      end if
      do ngp=1,iqp
        do iex=1,nof_variables
          rec_uleft(iex,l,ngp,i)=zero
        end do
      end do
    end do

    if (dg.eq.1) then
      do k=1,ideg_local
        do iex=1,nof_variables
          dg2fv(k,iex,i)=zero
        end do
      end do
    end if

    do ll=1,iadmis
      do l=1,n_faces
        if (dimensiona.eq.3) then
          if (ielem_types_faces(l,i).eq.5) then
            iqp=qp_quad
          else
            iqp=qp_triangle
          end if
        else
          iqp=qp_line
        end if

        do ngp=1,iqp
          ax=rec_qpoints(l,ngp,1,i)
          ay=rec_qpoints(l,ngp,2,i)
          az=zero
          if (dimensiona.eq.3) az=rec_qpoints(l,ngp,3,i)

          do iex=1,nof_variables
            resvec(iex)=zero
          end do

          if (ll.ne.1) then
            kk=0
            if (dimensiona.eq.3) then
              do degree=1,iorder2
                do zpow=0,degree
                  zterm=1.0d0
                  do ip=1,zpow
                    zterm=zterm*az
                  end do
                  do ypow=0,degree-zpow
                    xpow=degree-zpow-ypow
                    kk=kk+1
                    if (kk.le.idegfree2) then
                      xterm=1.0d0
                      do ip=1,xpow
                        xterm=xterm*ax
                      end do
                      yterm=1.0d0
                      do ip=1,ypow
                        yterm=yterm*ay
                      end do
                      raw_basis=xterm*yterm*zterm
                      phik=raw_basis-integ_basis_valuec(kk,i)*volinv
                      do iex=1,nof_variables
                        resvec(iex)=resvec(iex)+phik*rec_gradientsc(ll,kk,iex,i)
                      end do
                    end if
                  end do
                end do
              end do
            else
              do degree=1,iorder2
                do ypow=0,degree
                  xpow=degree-ypow
                  kk=kk+1
                  if (kk.le.idegfree2) then
                    xterm=1.0d0
                    do ip=1,xpow
                      xterm=xterm*ax
                    end do
                    yterm=1.0d0
                    do ip=1,ypow
                      yterm=yterm*ay
                    end do
                    raw_basis=xterm*yterm
                    phik=raw_basis-integ_basis_valuec(kk,i)*volinv
                    do iex=1,nof_variables
                      resvec(iex)=resvec(iex)+phik*rec_gradientsc(ll,kk,iex,i)
                    end do
                  end if
                end do
              end do
            end if
            do ktail=kk+1,idegfree2
              phik=-integ_basis_valuec(ktail,i)*volinv
              do iex=1,nof_variables
                resvec(iex)=resvec(iex)+phik*rec_gradientsc(ll,ktail,iex,i)
              end do
            end do
          else
            kk=0
            if (dimensiona.eq.3) then
              do degree=1,ielem_iorder(i)
                do zpow=0,degree
                  zterm=1.0d0
                  do ip=1,zpow
                    zterm=zterm*az
                  end do
                  do ypow=0,degree-zpow
                    xpow=degree-zpow-ypow
                    kk=kk+1
                    if (kk.le.ideg_local) then
                      xterm=1.0d0
                      do ip=1,xpow
                        xterm=xterm*ax
                      end do
                      yterm=1.0d0
                      do ip=1,ypow
                        yterm=yterm*ay
                      end do
                      raw_basis=xterm*yterm*zterm
                      phik=raw_basis-integ_basis_value(kk,i)*volinv
                      do iex=1,nof_variables
                        resvec(iex)=resvec(iex)+phik*rec_gradients(1,kk,iex,i)
                      end do
                      if (kk.le.idegfree2) then
                        do ll2=2,iadmis
                          do iex=1,nof_variables
                            resvec(iex)=resvec(iex)-lamc(ll2)*phik*rec_gradientsc(ll2,kk,iex,i)
                          end do
                        end do
                      end if
                    end if
                  end do
                end do
              end do
            else
              do degree=1,ielem_iorder(i)
                do ypow=0,degree
                  xpow=degree-ypow
                  kk=kk+1
                  if (kk.le.ideg_local) then
                    xterm=1.0d0
                    do ip=1,xpow
                      xterm=xterm*ax
                    end do
                    yterm=1.0d0
                    do ip=1,ypow
                      yterm=yterm*ay
                      end do
                      raw_basis=xterm*yterm
                      phik=raw_basis-integ_basis_value(kk,i)*volinv
                      do iex=1,nof_variables
                        resvec(iex)=resvec(iex)+phik*rec_gradients(1,kk,iex,i)
                      end do
                      if (kk.le.idegfree2) then
                        do ll2=2,iadmis
                          do iex=1,nof_variables
                            resvec(iex)=resvec(iex)-lamc(ll2)*phik*rec_gradientsc(ll2,kk,iex,i)
                          end do
                        end do
                      end if
                    end if
                  end do
                end do
              end if
              do ktail=kk+1,ideg_local
                phik=-integ_basis_value(ktail,i)*volinv
                do iex=1,nof_variables
                  resvec(iex)=resvec(iex)+phik*rec_gradients(1,ktail,iex,i)
                end do
                if (ktail.le.idegfree2) then
                  do ll2=2,iadmis
                    do iex=1,nof_variables
                      resvec(iex)=resvec(iex)-lamc(ll2)*phik*rec_gradientsc(ll2,ktail,iex,i)
                    end do
                  end do
                end if
              end do

            do iex=1,nof_variables
              resvec(iex)=resvec(iex)*inv_lamc1
            end do
          end if

          do iex=1,nof_variables
            rec_uleft(iex,l,ngp,i)=rec_uleft(iex,l,ngp,i)+(u_c_val(1,iex,i)+resvec(iex))*weno(iex,ll)
          end do
        end do
      end do

      if (dg.eq.1) then
        if (ll.eq.1) then
          do iex=1,nof_variables
            do k=1,ideg_local
              grad0(k)=rec_gradients(1,k,iex,i)
            end do
            do ll2=2,iadmis
              do k=1,idegfree2
                grad0(k)=grad0(k)-lamc(ll2)*rec_gradientsc(ll2,k,iex,i)
              end do
            end do
            do k=1,ideg_local
              dg2fv(k,iex,i)=dg2fv(k,iex,i)+(grad0(k)*inv_lamc1*weno(iex,1))
            end do
          end do
        else
          do iex=1,nof_variables
            do k=1,idegfree2
              dg2fv(k,iex,i)=dg2fv(k,iex,i)+rec_gradientsc(ll,k,iex,i)*weno(iex,ll)
            end do
          end do
        end if
      end if
    end do
  end do
#ifdef xpu
  !$omp end target teams distribute parallel do
#else
  !$omp end do
#endif
end subroutine cp_reconstruction_cweno_acc25








subroutine cp_reconstruction_weno_acc2(n)
  implicit none
  integer,intent(in)::n
  integer :: i,kmaxe,l,ngp,iqp,ll,k,j,iex
  integer :: iadmis,n_faces,ideg_local,degree,zpow,ypow,xpow,kk,ktail,ip
  real :: lwcx1,divbyzero,sumomega,denom
  real :: ax,ay,az,volinv
  real :: raw_basis,phik,xterm,yterm,zterm
  real :: lambdaal(1:gpu_max_typesten)
  real :: omegaatilde(1:gpu_max_typesten)
  real :: omegaal(1:gpu_max_typesten)
  real :: smooth(1:gpu_max_typesten)
  real :: weno(1:gpu_max_nvar,1:gpu_max_typesten)
  real :: resvec(1:gpu_max_nvar)

  kmaxe=xmpielrank(n)
#ifdef xpu

  !$omp target teams distribute parallel do &
  !$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,dg,qp_quad,qp_triangle,qp_line) &
  !$omp& firstprivate(lwci1,zero) &
  !$omp& map(alloc: xmpielrank) &
  !$omp& map(alloc: u_c_val, rec_uleft, dg2fv, ielem_wcx) &
  !$omp& map(alloc: ielem_admis, ielem_ifca, ielem_linc, ielem_types_faces) &
  !$omp& map(alloc: ielem_iorder, ielem_idegfree) &
  !$omp& map(alloc: rec_qpoints, rec_gradients, rec_indicator) &
  !$omp& map(alloc: rec_volume, integ_basis_value) &
  !$omp& private(i,l,ngp,iqp,ll,k,j,iex,iadmis,n_faces,ideg_local) &
  !$omp& private(degree,zpow,ypow,xpow,kk,ktail,ip,lwcx1,divbyzero,sumomega,denom) &
  !$omp& private(ax,ay,az,volinv,raw_basis,phik,xterm,yterm,zterm) &
  !$omp& private(lambdaal,omegaatilde,omegaal,smooth,weno,resvec)
#else
  !$omp do private(i,l,ngp,iqp,ll,k,j,iex,iadmis,n_faces,ideg_local) &
  !$omp& private(degree,zpow,ypow,xpow,kk,ktail,ip,lwcx1,divbyzero,sumomega,denom) &
  !$omp& private(ax,ay,az,volinv,raw_basis,phik,xterm,yterm,zterm) &
  !$omp& private(lambdaal,omegaatilde,omegaal,smooth,weno,resvec)
#endif
  do i=1,kmaxe
    iadmis=ielem_admis(i)
    n_faces=ielem_ifca(i)
    ideg_local=ielem_idegfree(i)
    lwcx1=ielem_linc(i)
    if (lwcx1.eq.zero) lwcx1=lwci1
    ielem_linc(i)=lwcx1
    divbyzero=1.0d-6
    volinv=1.0d0/rec_volume(1,1,i)

    do iex=1,nof_variables
      do ll=1,gpu_max_typesten
        smooth(ll)=zero
        omegaatilde(ll)=zero
        omegaal(ll)=zero
        lambdaal(ll)=zero
        weno(iex,ll)=zero
      end do

      do ll=1,iadmis
        smooth(ll)=zero
        do j=1,ideg_local
          do k=1,ideg_local
            smooth(ll)=smooth(ll)+rec_gradients(ll,j,iex,i)*rec_indicator(j,k,i)*rec_gradients(ll,k,iex,i)
          end do
        end do
      end do

      do ll=1,iadmis
        lambdaal(ll)=1.0d0
      end do
      lambdaal(1)=lwcx1

      do ll=1,iadmis
        denom=divbyzero+smooth(ll)
        denom=denom*denom*denom*denom
        omegaatilde(ll)=lambdaal(ll)/denom
      end do

      sumomega=zero
      do ll=1,iadmis
        sumomega=sumomega+omegaatilde(ll)
      end do

      do ll=1,iadmis
        omegaal(ll)=omegaatilde(ll)/sumomega
        weno(iex,ll)=omegaal(ll)
      end do
      if (iex.eq.1) ielem_wcx(i)=weno(1,1)
    end do

    do l=1,n_faces
      if (dimensiona.eq.3) then
        if (ielem_types_faces(l,i).eq.5) then
          iqp=qp_quad
        else
          iqp=qp_triangle
        end if
      else
        iqp=qp_line
      end if
      do ngp=1,iqp
        do iex=1,nof_variables
          rec_uleft(iex,l,ngp,i)=zero
        end do
      end do
    end do

    if (dg.eq.1) then
      do k=1,ideg_local
        do iex=1,nof_variables
          dg2fv(k,iex,i)=zero
        end do
      end do
    end if

    do ll=1,iadmis
      do l=1,n_faces
        if (dimensiona.eq.3) then
          if (ielem_types_faces(l,i).eq.5) then
            iqp=qp_quad
          else
            iqp=qp_triangle
          end if
        else
          iqp=qp_line
        end if

        do ngp=1,iqp
          ax=rec_qpoints(l,ngp,1,i)
          ay=rec_qpoints(l,ngp,2,i)
          az=zero
          if (dimensiona.eq.3) az=rec_qpoints(l,ngp,3,i)

          do iex=1,nof_variables
            resvec(iex)=zero
          end do

          kk=0
          if (dimensiona.eq.3) then
            do degree=1,ielem_iorder(i)
              do zpow=0,degree
                zterm=1.0d0
                do ip=1,zpow
                  zterm=zterm*az
                end do
                do ypow=0,degree-zpow
                  xpow=degree-zpow-ypow
                  kk=kk+1
                  if (kk.le.ideg_local) then
                    xterm=1.0d0
                    do ip=1,xpow
                      xterm=xterm*ax
                    end do
                    yterm=1.0d0
                    do ip=1,ypow
                      yterm=yterm*ay
                    end do
                    raw_basis=xterm*yterm*zterm
                    phik=raw_basis-integ_basis_value(kk,i)*volinv
                    do iex=1,nof_variables
                      resvec(iex)=resvec(iex)+phik*rec_gradients(ll,kk,iex,i)
                    end do
                  end if
                end do
              end do
            end do
          else
            do degree=1,ielem_iorder(i)
              do ypow=0,degree
                xpow=degree-ypow
                kk=kk+1
                if (kk.le.ideg_local) then
                  xterm=1.0d0
                  do ip=1,xpow
                    xterm=xterm*ax
                  end do
                  yterm=1.0d0
                  do ip=1,ypow
                    yterm=yterm*ay
                  end do
                  raw_basis=xterm*yterm
                  phik=raw_basis-integ_basis_value(kk,i)*volinv
                  do iex=1,nof_variables
                    resvec(iex)=resvec(iex)+phik*rec_gradients(ll,kk,iex,i)
                  end do
                end if
              end do
            end do
          end if

          do ktail=kk+1,ideg_local
            phik=-integ_basis_value(ktail,i)*volinv
            do iex=1,nof_variables
              resvec(iex)=resvec(iex)+phik*rec_gradients(ll,ktail,iex,i)
            end do
          end do

          do iex=1,nof_variables
            rec_uleft(iex,l,ngp,i)=rec_uleft(iex,l,ngp,i)+(u_c_val(1,iex,i)+resvec(iex))*weno(iex,ll)
          end do
        end do
      end do

      if (dg.eq.1) then
        do iex=1,nof_variables
          do k=1,ideg_local
            dg2fv(k,iex,i)=dg2fv(k,iex,i)+rec_gradients(ll,k,iex,i)*weno(iex,ll)
          end do
        end do
      end if
    end do
  end do
#ifdef xpu
  !$omp end target teams distribute parallel do
#else
  !$omp end do
#endif
end subroutine cp_reconstruction_weno_acc2





















subroutine cp_reconstruction_cweno_turb(iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif

integer :: kbasis
  integer, intent(in)  :: iconsidered
  real    :: divbyzero
  integer :: power

  integer :: i, l, ngp, iqp, ll, ll2, k, j, iex
  integer :: iadmis, n_faces, ideg_local, icompwrt
  integer :: ideg_max
  real    :: lwcx1, ax, ay, az
  real    :: tau_weno, sumomega
  real    :: lamc(1:7), lambdaal(1:7), omegaatilde(1:7), omegaal(1:7)
  real    :: smooth(1:7)
  real    :: resvec(1:gpu_max_extra_transport)
  real    :: leftv(1:gpu_max_extra_transport)
  real    :: inv_lamc1

  real :: phi(1:gpu_max_dof)

  ! grad0 scratch
  real :: grad0(1:gpu_max_dof)

  ! WENO weights needed later
  real :: weno(1:gpu_max_extra_transport, 1:7)

  ! -------------------------
  ! Setup
  ! -------------------------
  i          = iconsidered
  iadmis     = ielem_admis(i)
  n_faces    = ielem_ifca(i)
  ideg_local = idegfree
  lwcx1      = ielem_linc(i)

  divbyzero  = 1.0e-6
  power      = 4

  ! Linear weights (lamc) depend only on element, not variable
  lamc(1) = (1.0d0 - (1.0d0 / lwcx1))
  if (iadmis > 1) then
    lamc(2:iadmis) = (1.0d0 - lamc(1)) / real(iadmis-1)
  end if
  inv_lamc1 = 1.0d0 / lamc(1)




  ! -------------------------
  ! Compute WENO weights per variable
  ! -------------------------
  do iex = 1, turbulenceequations+passivescalar

    smooth(1:iadmis)      = 0.0d0
    omegaatilde(1:iadmis) = 0.0d0
    omegaal(1:iadmis)     = 0.0d0

    ! ---- build grad0(:)  ----
    ! grad0 = inv_lamc1 * ( rec_gradients(1) - sum_{m=2..iadmis} lamc(m)*rec_gradientsc(m) )
    grad0(1:ideg_local) = rec_gradients2(1, 1:ideg_local, iex, i)
    do ll = 2, iadmis
      do k = 1, idegfree2
        grad0(k) = grad0(k) - lamc(ll) * rec_gradientsc2(ll, k, iex, i)
      end do
    end do
    grad0(1:ideg_local) = inv_lamc1 * grad0(1:ideg_local)

    ! ---- smoothness indicator for ll=1  ----
    smooth(1) = 0.0d0
    do j = 1, ideg_local
      do k = 1, ideg_local
        smooth(1) = smooth(1) + grad0(j) * rec_indicator(j,k,i) * grad0(k)
      end do
    end do

    ! ---- smoothness indicators for ll>=2  ----
    do ll = 2, iadmis
      smooth(ll) = 0.0d0
      do j = 1, idegfree2
        do k = 1, idegfree2
          smooth(ll) = smooth(ll) + rec_gradientsc2(ll,j,iex,i) * rec_indicatorc(j,k,i) * rec_gradientsc2(ll,k,iex,i)
        end do
      end do
    end do

    ! ---- lambda ----
    lambdaal(1:iadmis) = lamc(1:iadmis)

    ! ---- omega tilde ----
    if (wenoz .eq. 1) then
      tau_weno = 0.0d0
      do ll = 1, iadmis
        tau_weno = tau_weno + abs(smooth(1) - smooth(ll))
      end do
      if (iadmis > 1) tau_weno = tau_weno / real(iadmis-1)

      do ll = 1, iadmis
        omegaatilde(ll) = lambdaal(ll) * (1.0d0 + (tau_weno / (divbyzero + smooth(ll)))**4)
      end do
    else
      do ll = 1, iadmis
        omegaatilde(ll) = lambdaal(ll) / ((divbyzero + smooth(ll))**4)
      end do
    end if

    ! ---- normalize ----
    sumomega = 0.0d0
    do ll = 1, iadmis
      sumomega = sumomega + omegaatilde(ll)
    end do

    do ll = 1, iadmis
      omegaal(ll)   = omegaatilde(ll) / sumomega
      weno(iex,ll)  = omegaal(ll)
    end do

    if (iex .eq. 1) ielem_wcx(i) = weno(1,1)

  end do

  ! -------------------------
  ! Reconstruction
  ! -------------------------
  do l = 1, n_faces
    if (dimensiona .eq. 3) then
      if (ielem_types_faces(l,i) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if
    do ngp = 1, iqp
      do iex = 1, turbulenceequations+passivescalar
        rec_uleftturb(iex,l,ngp,i) = 0.0d0
      end do
    end do
  end do
!   if (dg .eq. 1) dg2fv(1:ideg_local,:,i) = 0.0d0

  do ll = 1, iadmis

    do l = 1, n_faces

      if (dimensiona .eq. 3) then
        if (ielem_types_faces(l,i) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if

      do ngp = 1, iqp

        ax = rec_qpoints(l,ngp,1,i)
        ay = rec_qpoints(l,ngp,2,i)
        if (dimensiona .eq. 3) az = rec_qpoints(l,ngp,3,i)

        resvec(:) = 0.0d0

        if (ll .ne. 1) then
          ! ---- low-order basis (degree idegfree2) -> phi(1:idegfree2) ----

          if (dimensiona .eq. 3) then
              do kbasis=1,idegfree2
                phi(kbasis) = basis_rec_value(n,ax,ay,az,iorder2,i,idegfree2,1,kbasis)
              end do
          else
              do kbasis=1,idegfree2
                phi(kbasis) = basis_rec2d_value(n,ax,ay,iorder2,i,idegfree2,1,kbasis)
              end do
          end if

          do k = 1, idegfree2
            do iex = 1, turbulenceequations+passivescalar
              resvec(iex) = resvec(iex) + phi(k) * rec_gradientsc2(ll, k, iex, i)
            end do
          end do

        else
          ! ---- full-order basis (degree ideg_local) -> phi(1:ideg_local) ----

          if (dimensiona .eq. 3) then
              do kbasis=1,ideg_local
                phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg_local,0,kbasis)
              end do
          else
              do kbasis=1,ideg_local
                phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg_local,0,kbasis)
              end do
          end if

          ! resvec(iex) = inv_lamc1 * ( dot(phi, rec_gradients(1,:,iex)) - sum_{m=2..iadmis} lamc(m)*dot(phi, rec_gradientsc(m,:,iex)) )
          do iex = 1, turbulenceequations+passivescalar
            ! dot(phi, rec_gradients(1))
            do k = 1, ideg_local
              resvec(iex) = resvec(iex) + phi(k) * rec_gradients2(1, k, iex, i)
            end do

            ! subtract lamc(m)*dot(phi, rec_gradientsc(m))
            do ll2 = 2, iadmis
              do k = 1, idegfree2
                resvec(iex) = resvec(iex) - lamc(ll2) * phi(k) * rec_gradientsc2(ll2, k, iex, i)
              end do
            end do

            resvec(iex) = inv_lamc1 * resvec(iex)
          end do

        end if

        call extrapolate_boundt(resvec, l, ngp, i, ll, weno)

      end do
    end do



  end do






end subroutine cp_reconstruction_cweno_turb



subroutine cp_reconstruction_weno_turb(iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif

integer :: kbasis
  integer, intent(in)  :: iconsidered
  real    :: divbyzero
  integer:: power

  integer :: i, l, ngp, iqp, ll, k, j, iex
  integer :: iadmis, n_faces, ideg_local, icompwrt
  real    :: lwcx1, ax, ay, az
  real    :: sumomega
  real    :: lambdaal(1:7), omegaatilde(1:7), omegaal(1:7)
  real    :: smooth(1:7)
  real    :: resvec(1:gpu_max_extra_transport)
  real    :: leftv(1:gpu_max_extra_transport)
  real    :: weno(1:gpu_max_extra_transport,1:7)

  real :: phi(1:gpu_max_dof)

  i          = iconsidered
  iadmis     = ielem_admis(i)
  n_faces    = ielem_ifca(i)
  ideg_local = idegfree
  lwcx1      = ielem_linc(i)

  divbyzero  = 1.0e-6
  power      = 4


  ! -------------------------
  ! WENO weights per variable (quadratic form on-the-fly)
  ! -------------------------
  do iex = 1, turbulenceequations+passivescalar

    smooth(1:iadmis)      = 0.0d0
    omegaatilde(1:iadmis) = 0.0d0
    omegaal(1:iadmis)     = 0.0d0

    do ll = 1, iadmis
      smooth(ll) = 0.0d0
      do j = 1, ideg_local
        do k = 1, ideg_local
          smooth(ll) = smooth(ll) + rec_gradients(ll,j,iex,i) * rec_indicator(j,k,i) * rec_gradients2(ll,k,iex,i)
        end do
      end do
    end do

    lambdaal(1:iadmis) = 1.0d0
    lambdaal(1)        = lwcx1

    do ll = 1, iadmis
      omegaatilde(ll) = lambdaal(ll) / ((divbyzero + smooth(ll))**4)
    end do

    sumomega = 0.0d0
    do ll = 1, iadmis
      sumomega = sumomega + omegaatilde(ll)
    end do

    do ll = 1, iadmis
      omegaal(ll)  = omegaatilde(ll) / sumomega
      weno(iex,ll) = omegaal(ll)
    end do

    if (iex .eq. 1) ielem_wcx(i) = weno(1,1)

  end do

  ! -------------------------
  ! Reconstruction (reuse phi)
  ! -------------------------
  do l = 1, n_faces
    if (dimensiona .eq. 3) then
      if (ielem_types_faces(l,i) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if
    do ngp = 1, iqp
      do iex = 1, turbulenceequations+passivescalar
        rec_uleftturb(iex,l,ngp,i) = 0.0d0
      end do
    end do
  end do


  do ll = 1, iadmis

    do l = 1, n_faces

      if (dimensiona .eq. 3) then
        if (ielem_types_faces(l,i) .eq. 5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if

      do ngp = 1, iqp

        ax = rec_qpoints(l,ngp,1,i)
        ay = rec_qpoints(l,ngp,2,i)
        if (dimensiona .eq. 3) az = rec_qpoints(l,ngp,3,i)

        icompwrt = 0
        if (dimensiona .eq. 3) then
            do kbasis=1,ideg_local
              phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg_local,icompwrt,kbasis)
            end do
        else
            do kbasis=1,ideg_local
              phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg_local,icompwrt,kbasis)
            end do
        end if

        resvec(:) = 0.0d0
        do k = 1, ideg_local
          do iex = 1, turbulenceequations+passivescalar
            resvec(iex) = resvec(iex) + phi(k) * rec_gradients2(ll,k,iex,i)
          end do
        end do

        call extrapolate_boundt(resvec, l, ngp, i, ll, weno)

      end do
    end do



  end do




end subroutine cp_reconstruction_weno_turb






subroutine cp_reconstruction_weno(iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif

integer :: kbasis
  integer, intent(in)  :: iconsidered
  real    :: divbyzero
  integer:: power

  integer :: i, l, ngp, iqp, ll, k, j, iex
  integer :: iadmis, n_faces, ideg_local, icompwrt
  real    :: lwcx1, ax, ay, az
  real    :: sumomega
  real    :: lambdaal(1:7), omegaatilde(1:7), omegaal(1:7)
  real    :: smooth(1:7)
  real    :: resvec(1:gpu_max_nvar)
  real    :: leftv(1:gpu_max_nvar)
  real    :: weno(1:gpu_max_nvar,1:7)

  real :: phi(1:gpu_max_dof)

  i          = iconsidered
  iadmis     = ielem_admis(i)
  n_faces    = ielem_ifca(i)
  ideg_local = ielem_idegfree(i)
  lwcx1      = ielem_linc(i)

  divbyzero  = 1.0e-6
  power      = 4


  ! -------------------------
  ! WENO weights per variable (quadratic form on-the-fly)
  ! -------------------------
  do iex = 1, nof_variables

    smooth(1:iadmis)      = 0.0d0
    omegaatilde(1:iadmis) = 0.0d0
    omegaal(1:iadmis)     = 0.0d0

    do ll = 1, iadmis
      smooth(ll) = 0.0d0
      do j = 1, ideg_local
        do k = 1, ideg_local
          smooth(ll) = smooth(ll) + rec_gradients(ll,j,iex,i) * rec_indicator(j,k,i) * rec_gradients(ll,k,iex,i)
        end do
      end do
    end do

    lambdaal(1:iadmis) = 1.0d0
    lambdaal(1)        = lwcx1

    do ll = 1, iadmis
      omegaatilde(ll) = lambdaal(ll) / ((divbyzero + smooth(ll))**4)
    end do

    sumomega = 0.0d0
    do ll = 1, iadmis
      sumomega = sumomega + omegaatilde(ll)
    end do

    do ll = 1, iadmis
      omegaal(ll)  = omegaatilde(ll) / sumomega
      weno(iex,ll) = omegaal(ll)
    end do

    if (iex .eq. 1) ielem_wcx(i) = weno(1,1)

  end do

  ! -------------------------
  ! Reconstruction (reuse phi)
  ! -------------------------
  do l = 1, n_faces
    if (dimensiona .eq. 3) then
      if (ielem_types_faces(l,i) .eq. 5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if
    do ngp = 1, iqp
      do iex = 1, nof_variables
        rec_uleft(iex,l,ngp,i) = 0.0d0
      end do
    end do
  end do
  if (dg .eq. 1) dg2fv(1:ideg_local,:,i) = 0.0d0

  do ll = 1, iadmis

    do l = 1, n_faces

      if (dimensiona .eq. 3) then
       if (ielem_types_faces(l,i) .eq. 5) then
          iqp = qp_quad
       else
         iqp = qp_triangle
       end if
     else
       iqp = qp_line
     end if

      do ngp = 1, iqp

        ax = rec_qpoints(l,ngp,1,i)
        ay = rec_qpoints(l,ngp,2,i)
        if (dimensiona .eq. 3) az = rec_qpoints(l,ngp,3,i)

        icompwrt = 0
        if (dimensiona .eq. 3) then
            do kbasis=1,ideg_local
              phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg_local,icompwrt,kbasis)
            end do
        else
            do kbasis=1,ideg_local
              phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg_local,icompwrt,kbasis)
            end do
        end if

        resvec(:) = 0.0d0
        do k = 1, ideg_local
          do iex = 1, nof_variables
            resvec(iex) = resvec(iex) + phi(k) * rec_gradients(ll,k,iex,i)
          end do
        end do

        call extrapolate_bound(resvec, l, ngp, i, ll, weno)

      end do
    end do

   if (dg .eq. 1) then
     do iex = 1, nof_variables
       dg2fv(1:ideg_local,iex,i) = dg2fv(1:ideg_local,iex,i) + rec_gradients(ll,1:ideg_local,iex,i) * weno(iex,ll)
     end do
   end if

  end do

  ! -------------------------
  ! Convert to conservative
  ! -------------------------
  if (wenwrt .eq. 3) then
   do l = 1, n_faces
     if (dimensiona .eq. 3) then
       if (ielem_types_faces(l,i) .eq. 5) then
         iqp = qp_quad
       else
         iqp = qp_triangle
       end if
     else
       iqp = qp_line
     end if

     do ngp = 1, iqp
       leftv(1:nof_variables) = rec_uleft(1:nof_variables,l,ngp,i)
       call prim2cons(n,leftv)
       rec_uleft(1:nof_variables,l,ngp,i) = leftv(1:nof_variables)
     end do
   end do
  end if



end subroutine cp_reconstruction_weno










subroutine wenoweights_cons(n)
implicit none
!> @brief
!> subroutine for weno type reconstruction in 3d
integer,intent(in)::n
integer::i,ii
integer::kmaxe
kmaxe=xmpielrank(n)

if (ees.eq.5)then


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
ielem_linc(i)=lwci1
            if (adda.eq.1) call adda_filter(n,i)

               call cp_reconstruction_cweno(i)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


else

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
ielem_linc(i)=lwci1
                 if (adda.eq.1) call adda_filter(n,i)

                call cp_reconstruction_weno(i)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end if








end subroutine wenoweights_cons



subroutine wenoweights_turb(n)
implicit none
!> @brief
!> subroutine for weno type reconstruction in 3d
integer,intent(in)::n
integer::i
integer::iconsidered,kmaxe
kmaxe=xmpielrank(n)

if (ees.eq.5)then

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered)
#else
!$omp do
#endif
do i=1,kmaxe
iconsidered=i
ielem_linc(i)=lwci1


                call cp_reconstruction_cweno_turb(iconsidered)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


else

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered)
#else
!$omp do
#endif
do i=1,kmaxe
iconsidered=i
ielem_linc(i)=lwci1


                call cp_reconstruction_weno_turb(iconsidered)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





end if



end subroutine wenoweights_turb



subroutine extrapolate_bound(resvec,facex,pointx,iconsidered,llx,weno)
implicit none
!> @brief
!> subroutine for extrapolating the reconstructed solution at the cell interfaces
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::facex,pointx,iconsidered,llx
real,intent(in) :: weno(1:nof_variables,1:7)
real,intent(in) :: resvec(1:nof_variables)
real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal




				    if (wenwrt.eq.3) then  ! primitive accumulation, but increments are in conservative space
                    leftv(1:nof_variables) = u_c_val(1,1:nof_variables,iconsidered)
                    call cons2prim(n,leftv,mp_pinfl,gammal)

                    rec_uleft(1:nof_variables,facex,pointx,iconsidered) = &
                      rec_uleft(1:nof_variables,facex,pointx,iconsidered) + ((leftv(1:nof_variables)+resvec(1:nof_variables)) * weno(1:nof_variables,llx))

                  else
                    rec_uleft(1:nof_variables,facex,pointx,iconsidered) = &
                      rec_uleft(1:nof_variables,facex,pointx,iconsidered) + &
                      (u_c_val(1,1:nof_variables,iconsidered) + resvec(1:nof_variables)) * weno(1:nof_variables,llx)
                  end if






end subroutine extrapolate_bound



subroutine extrapolate_boundt(resvec,facex,pointx,iconsidered,llx,weno)
implicit none
!> @brief
!> subroutine for extrapolating the reconstructed solution at the cell interfaces
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::facex,pointx,iconsidered,llx
real,intent(in) :: weno(turbulenceequations+passivescalar,typesten)
real,intent(in) :: resvec(turbulenceequations+passivescalar)
real::mp_pinfl,gammal


				     rec_uleftturb(1:turbulenceequations+passivescalar,facex,pointx,iconsidered)=rec_uleftturb(1:turbulenceequations+passivescalar,facex,pointx,iconsidered)&
				     +(u_ct_val(1,1:turbulenceequations+passivescalar,iconsidered)+resvec(1:turbulenceequations+passivescalar))*weno(1:turbulenceequations+passivescalar,llx)





end subroutine extrapolate_boundt


subroutine wenoweights_char(n)
implicit none
!> @brief
!> subroutine for weno type reconstruction in 3d
integer,intent(in)::n
integer :: ii, i, iconsidered,kmaxe
kmaxe=xmpielrank(n)

if (ees.eq.5)then



#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered)
#else
!$omp do
#endif
do i=1,kmaxe
iconsidered=i

ielem_linc(i)=lwci1

                call characteristic_reconstruction_cweno(iconsidered)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



else


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered)
#else
!$omp do
#endif
do i=1,kmaxe
iconsidered=i

ielem_linc(i)=lwci1



                call characteristic_reconstruction_weno(iconsidered)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end if




end subroutine wenoweights_char





subroutine characteristic_reconstruction_cweno(iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer :: kbasis
integer,intent(in)    :: iconsidered
integer               ::  power
real                  :: divbyzero

integer :: i, l, ngp, iqp, icompwrt, k, j, c
integer :: iadmis, n_faces, ideg
real    :: angle1, angle2, nx, ny, nz, lwcx1, ax, ay, az
real    :: veigl(gpu_max_nvar), veigr(gpu_max_nvar), rveigl(gpu_max_nvar), rveigr(gpu_max_nvar)
real    :: eigvl(gpu_max_nvar,gpu_max_nvar), eigvr(gpu_max_nvar,gpu_max_nvar)

real    :: lamc(gpu_max_typesten), inv_lamc1
real    :: limiteddw_c(0:gpu_max_dof)
real    :: phi(1:gpu_max_dof)
real    :: s



divbyzero  = 1.0e-6
power      = 4

i      = iconsidered
iadmis = ielem_admis(i)
n_faces= ielem_ifca(i)
ideg   = idegfree
lwcx1  = ielem_linc(i)

! linear weights for EES=5 (depend only on element)
lamc(:) = 0.0d0
lamc(1) = (1.0d0 - (1.0d0/lwcx1))
if (iadmis > 1) lamc(2:iadmis) = (1.0d0 - lamc(1)) / real(iadmis-1)
inv_lamc1 = 1.0d0 / lamc(1)

do l=1,n_faces

  angle1 = ielem_faceanglex(l,i)
  angle2 = ielem_faceangley(l,i)

  if (dimensiona.eq.3) then
    nx = cos(angle1)*sin(angle2)
    ny = sin(angle1)*sin(angle2)
    nz = cos(angle2)
  else
    nx = angle1
    ny = angle2
    nz = 0.0d0
  end if

  veigl(1:nof_variables) = u_c_val(1,1:nof_variables,i)

  if (dimensiona.eq.3) then
    call rotatef(n,rveigl,veigl,angle1,angle2)
  else
    call rotatef2d(n,rveigl,veigl,angle1,angle2)
  end if

  call weno_neighbour(i,l,veigl,veigr,nx,ny,nz,angle1,angle2)

  if (dimensiona.eq.3) then
    call rotatef(n,rveigr,veigr,angle1,angle2)
    call compute_eigenvectors(n,rveigl,rveigr,eigvl,eigvr,gamma)
  else
    call rotatef2d(n,rveigr,veigr,angle1,angle2)
    call compute_eigenvectors2d(n,rveigl,rveigr,eigvl,eigvr,gamma)
  end if

  ! quadrature points on this face
  if (dimensiona.eq.3) then
    if (ielem_types_faces(l,i).eq.5) then
      iqp = qp_quad
    else
      iqp = qp_triangle
    end if
  else
    iqp = qp_line
  end if

	  ! initialize output on this face for accumulating over characteristic components)
	  do ngp=1,iqp
	    do j=1,nof_variables
	      rec_uleft(j,l,ngp,i) = 0.0d0
	    end do
	  end do

  ! loop characteristic component-by-component
  do c = 1, nof_variables

    call char_build_limiteddw_cweno(i, iadmis, ideg, lamc, inv_lamc1, eigvl, c, divbyzero, power, limiteddw_c)

    do ngp=1,iqp

      ax = rec_qpoints(l,ngp,1,i)
      ay = rec_qpoints(l,ngp,2,i)
      if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)

      icompwrt = 0
      if (ideg > 0) then
        if (dimensiona.eq.3) then
            do kbasis=1,ideg
              phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,0,kbasis)
            end do
        else
            do kbasis=1,ideg
              phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg,0,kbasis)
            end do
        end if
      end if

      s = limiteddw_c(0)
      do k = 1, ideg
        s = s + phi(k) * limiteddw_c(k)
      end do

      do j = 1, nof_variables
        rec_uleft(j,l,ngp,i) = rec_uleft(j,l,ngp,i) + eigvr(j,c) * s
      end do

    end do ! ngp

  end do ! c

end do  ! faces

end subroutine characteristic_reconstruction_cweno



subroutine characteristic_reconstruction_weno(iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer :: kbasis
integer,intent(in)    :: iconsidered
integer               :: power
real                  :: divbyzero

integer :: i, l, ngp, iqp, icompwrt, k, j, c
integer :: iadmis, n_faces, ideg
real    :: angle1, angle2, nx, ny, nz, lwcx1, ax, ay, az
real    :: veigl(gpu_max_nvar), veigr(gpu_max_nvar), rveigl(gpu_max_nvar), rveigr(gpu_max_nvar)
real    :: eigvl(gpu_max_nvar,gpu_max_nvar), eigvr(gpu_max_nvar,gpu_max_nvar)

real    :: limiteddw_c(0:gpu_max_dof)
real    :: phi(1:gpu_max_dof)
real    :: s

i      = iconsidered
iadmis = ielem_admis(i)
n_faces= ielem_ifca(i)
ideg   = idegfree
lwcx1  = ielem_linc(i)

divbyzero  = 1.0e-6
power      = 4



do l=1,n_faces

  angle1 = ielem_faceanglex(l,i)
  angle2 = ielem_faceangley(l,i)

  if (dimensiona.eq.3) then
    nx = cos(angle1)*sin(angle2)
    ny = sin(angle1)*sin(angle2)
    nz = cos(angle2)
  else
    nx = angle1
    ny = angle2
    nz = 0.0d0
  end if

  veigl(1:nof_variables) = u_c_val(1,1:nof_variables,i)

  if (dimensiona.eq.3) then
    call rotatef(n,rveigl,veigl,angle1,angle2)
  else
    call rotatef2d(n,rveigl,veigl,angle1,angle2)
  end if

  call weno_neighbour(i,l,veigl,veigr,nx,ny,nz,angle1,angle2)

  if (dimensiona.eq.3) then
    call rotatef(n,rveigr,veigr,angle1,angle2)
    call compute_eigenvectors(n,rveigl,rveigr,eigvl,eigvr,gamma)
  else
    call rotatef2d(n,rveigr,veigr,angle1,angle2)
    call compute_eigenvectors2d(n,rveigl,rveigr,eigvl,eigvr,gamma)
  end if

  ! quadrature points on this face
  if (dimensiona.eq.3) then
    if (ielem_types_faces(l,i).eq.5) then
      iqp = qp_quad
    else
      iqp = qp_triangle
    end if
  else
    iqp = qp_line
  end if

	  ! initialize output on this face  to accumulate  characteristic components
	  do ngp=1,iqp
	    do j=1,nof_variables
	      rec_uleft(j,l,ngp,i) = 0.0d0
	    end do
	  end do

  do c = 1, nof_variables

    call char_build_limiteddw_weno(i, iadmis, ideg, lwcx1, eigvl, c, divbyzero, power, limiteddw_c)

    do ngp=1,iqp

      ax = rec_qpoints(l,ngp,1,i)
      ay = rec_qpoints(l,ngp,2,i)
      if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)

      icompwrt = 0
      if (ideg > 0) then
        if (dimensiona.eq.3) then
            do kbasis=1,ideg
              phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,0,kbasis)
            end do
        else
            do kbasis=1,ideg
              phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg,0,kbasis)
            end do
        end if
      end if

      s = limiteddw_c(0)
      do k = 1, ideg
        s = s + phi(k) * limiteddw_c(k)
      end do

      do j = 1, nof_variables
        rec_uleft(j,l,ngp,i) = rec_uleft(j,l,ngp,i) + eigvr(j,c) * s
      end do

    end do ! ngp

  end do ! c

end do  ! faces

end subroutine characteristic_reconstruction_weno



subroutine char_build_limiteddw_cweno(i, iadmis, ideg, lamc, inv_lamc1, eigvl, c, divbyzero, power, limiteddw_c)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer, intent(in) :: i, iadmis, ideg, c, power
real,    intent(in) :: lamc(typesten), inv_lamc1, divbyzero
real,    intent(in) :: eigvl(nof_variables,nof_variables)
real,    intent(out):: limiteddw_c(0:idegfree)

integer :: ll, ll2, k, j, itarget, kdot
real    :: smooth_c(gpu_max_typesten), omega_c(gpu_max_typesten), omegat(gpu_max_typesten)
real    :: tau_weno, sumomegat
real    :: a(gpu_max_dof), tmp(gpu_max_dof)
real    :: corr, gkj

limiteddw_c(0:ideg) = 0.0d0
smooth_c(1:iadmis)  = 0.0d0

!----------------------------
! PASS 1: smoothness indicators for this characteristic component
!----------------------------
do ll = 1, iadmis

  if (ll .eq. 1) then
    itarget = ideg

    do k = 1, itarget
      a(k) = 0.0d0
      do j = 1, nof_variables
        corr = 0.0d0
        if (k.le.idegfree2)then
        do ll2 = 2, iadmis
          corr = corr + lamc(ll2) * rec_gradientsc(ll2, k, j, i)
        end do
        end if

        gkj  = inv_lamc1 * (rec_gradients(1, k, j, i) - corr)
        a(k) = a(k) + eigvl(c,j) * gkj
      end do
    end do

    tmp(1:itarget) = 0.0d0
    do k = 1, itarget
      do j = 1, itarget
        tmp(k) = tmp(k) + rec_indicator(k,j,i) * a(j)
      end do
    end do

    smooth_c(ll) = 0.0d0
    do kdot = 1, itarget
      smooth_c(ll) = smooth_c(ll) + a(kdot) * tmp(kdot)
    end do

  else
    itarget = idegfree2

    do k = 1, itarget
      a(k) = 0.0d0
      do j = 1, nof_variables
        a(k) = a(k) + eigvl(c,j) * rec_gradientsc(ll, k, j, i)
      end do
    end do

    tmp(1:itarget) = 0.0d0
    do k = 1, itarget
      do j = 1, itarget
        tmp(k) = tmp(k) + rec_indicatorc(k,j,i) * a(j)
      end do
    end do

    smooth_c(ll) = 0.0d0
    do kdot = 1, itarget
      smooth_c(ll) = smooth_c(ll) + a(kdot) * tmp(kdot)
    end do
  end if

end do

!----------------------------
! nonlinear weights omega_c
!----------------------------
sumomegat = 0.0d0

if (wenoz .eq. 1) then
  tau_weno = 0.0d0
  do ll = 1, iadmis
    tau_weno = tau_weno + abs(smooth_c(1) - smooth_c(ll))
  end do
  if (iadmis > 1) tau_weno = tau_weno / real(iadmis-1)

  do ll = 1, iadmis
    omegat(ll) = lamc(ll) * (1.0d0 + (tau_weno/(divbyzero+smooth_c(ll)))**4)
    sumomegat  = sumomegat + omegat(ll)
  end do
else
  do ll = 1, iadmis
    omegat(ll) = lamc(ll) / ((divbyzero+smooth_c(ll))**4)
    sumomegat  = sumomegat + omegat(ll)
  end do
end if

do ll = 1, iadmis
  omega_c(ll) = omegat(ll) / sumomegat
end do

!----------------------------
! PASS 2: build limiteddw_c(k)
!----------------------------
do ll = 1, iadmis

  ! degree 0 in characteristic space (cell average)
  a(1) = 0.0d0
  do j = 1, nof_variables
    a(1) = a(1) + eigvl(c,j) * u_c_val(1,j,i)
  end do
  limiteddw_c(0) = limiteddw_c(0) + omega_c(ll) * a(1)

  if (ll .eq. 1) then
    itarget = ideg
    do k = 1, itarget
      a(k) = 0.0d0
      do j = 1, nof_variables
        corr = 0.0d0
         if (k.le.idegfree2)then
        do ll2 = 2, iadmis
          corr = corr + lamc(ll2) * rec_gradientsc(ll2, k, j, i)
        end do
        end if
        gkj  = inv_lamc1 * (rec_gradients(1, k, j, i) - corr)
        a(k) = a(k) + eigvl(c,j) * gkj
      end do
      limiteddw_c(k) = limiteddw_c(k) + omega_c(ll) * a(k)
    end do
  else
    itarget = idegfree2
    do k = 1, itarget
      a(k) = 0.0d0
      do j = 1, nof_variables
        a(k) = a(k) + eigvl(c,j) * rec_gradientsc(ll, k, j, i)
      end do
      limiteddw_c(k) = limiteddw_c(k) + omega_c(ll) * a(k)
    end do
  end if

end do

end subroutine char_build_limiteddw_cweno




subroutine char_build_limiteddw_weno(i, iadmis, ideg, lwcx1, eigvl, c, divbyzero, power, limiteddw_c)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer, intent(in) :: i, iadmis, ideg, c, power
real,    intent(in) :: lwcx1, divbyzero
real,    intent(in) :: eigvl(nof_variables,nof_variables)
real,    intent(out):: limiteddw_c(0:idegfree)

integer :: ll, k, j, kdot
real    :: smooth_c(gpu_max_typesten), omega_c(gpu_max_typesten), omegat(gpu_max_typesten)
real    :: sumomegat, lambda_ll
real    :: a(gpu_max_dof), tmp(gpu_max_dof)

limiteddw_c(0:ideg) = 0.0d0
smooth_c(1:iadmis)  = 0.0d0

!----------------------------
! PASS 1: smoothness indicators
!----------------------------
do ll = 1, iadmis

  do k = 1, ideg
    a(k) = 0.0d0
    do j = 1, nof_variables
      a(k) = a(k) + eigvl(c,j) * rec_gradients(ll, k, j, i)
    end do
  end do

  tmp(1:ideg) = 0.0d0
  do k = 1, ideg
    do j = 1, ideg
      tmp(k) = tmp(k) + rec_indicator(k,j,i) * a(j)
    end do
  end do

  smooth_c(ll) = 0.0d0
  do kdot = 1, ideg
    smooth_c(ll) = smooth_c(ll) + a(kdot) * tmp(kdot)
  end do

end do

!----------------------------
! nonlinear weights omega_c
!----------------------------
sumomegat = 0.0d0
do ll = 1, iadmis
  lambda_ll = 1.0d0
  if (ll .eq. 1) lambda_ll = lwcx1
  omegat(ll) = lambda_ll / ((divbyzero + smooth_c(ll))**4)
  sumomegat  = sumomegat + omegat(ll)
end do

do ll = 1, iadmis
  omega_c(ll) = omegat(ll) / sumomegat
end do

!----------------------------
! PASS 2: build limiteddw_c(k)
!----------------------------
do ll = 1, iadmis

  a(1) = 0.0d0
  do j = 1, nof_variables
    a(1) = a(1) + eigvl(c,j) * u_c_val(1,j,i)
  end do
  limiteddw_c(0) = limiteddw_c(0) + omega_c(ll) * a(1)

  do k = 1, ideg
    a(k) = 0.0d0
    do j = 1, nof_variables
      a(k) = a(k) + eigvl(c,j) * rec_gradients(ll, k, j, i)
    end do
    limiteddw_c(k) = limiteddw_c(k) + omega_c(ll) * a(k)
  end do

end do

end subroutine char_build_limiteddw_weno



subroutine weno_neighbour(iconsidered,facex,veigl,veigr,nx,ny,nz,angle1,angle2)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:nof_variables),intent(inout)::veigr
real,dimension(1:nof_variables),intent(inout)::veigl
real,intent(in)::nx,ny,nz,angle1,angle2
integer,intent(in)::iconsidered,facex
integer::idummy
real::mp_pinfl,gammal
integer::i,j,k,l,var2,b_code,n_node,nf,lf,rowf
real,dimension(1:gpu_max_nvar)::leftv,srf_speed,srf_speedrot,rightv
real,dimension(1:gpu_max_dim)::pox,poy,poz,cords
real,dimension(1:8,1:gpu_max_dim)::vext,nodes_list
real,dimension(1:gpu_max_extra_transport)::cturbl,cturbr
real,dimension(1:gpu_max_nvar_total)::cright_rot,cleft_rot
integer::ibfc



if (dimensiona.eq.3)then

l=facex
i=iconsidered


 if (ielem_interior(i).eq.0)then
     veigr(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i));

 else


                if (rec_mrf(i).eq.1)then
					srf_speed(2:4)=rec_rotvel(l,1,1:3,i)
					call rotatef(n,srf_speedrot,srf_speed,angle1,angle2)
                end if
				if (ielem_ineighb(l,i).eq.n)then	!my cpu only
				      if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
					  if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then	!periodic in my cpu
					  veigr(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
					  idummy=1
						  if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
	                        veigr(2:4)=rotate_per_1(veigr(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
						  end if
					  else
					  !not periodic ones in my cpu

					   call coordinates_face_innerx(n,iconsidered,facex,vext,nodes_list)


					    if (ielem_types_faces(facex,iconsidered).eq.5)then
                                            n_node=4
                                    else
                                            n_node=3
                                    end if

				  cords(1:3)=zero
 				  call cordinates3(n,nodes_list,n_node,cords(1:3))

				  poy(1)=cords(2)
				  pox(1)=cords(1)
				  poz(1)=cords(3)

 				  leftv(1:nof_variables)=veigl(1:nof_variables)
				  b_code=ibound_icode(ielem_ibounds(l,i))
 				 call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

				  veigr(1:nof_variables)=rightv(1:nof_variables)
				      	  end if
				      else
				      !fluid neighbour
				      veigr(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
				      end if
				else
			      !other my cpu
				    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
					  if ((ibound_icode(ielem_ibounds(l,i)).eq.5).or.(ibound_icode(ielem_ibounds(l,i)).eq.50))then	!periodic in other cpu
					!  veigr(1:nof_variables)=(iexsolhir(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
					!(rec_ihexl(1,ielem_indexi(l,i),i),1:nof_variables))

					 nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
                     lf=rec_ihexl(1,ielem_indexi(L,i),i)
                     rowf=halo_offset(nf) + lf - 1
                     veigr(1:nof_variables)=solhir(rowf,1:nof_variables)



					  idummy=1
						  if ((per_rot.eq.1).and.(ibound_icode(ielem_ibounds(l,i)).eq.50))then
	                        veigr(2:4)=rotate_per_1(veigr(2:4),ibound_icode(ielem_ibounds(l,i)),angle_per)
						  end if
					  end if
				    else

! 				      veigr(1:nof_variables)=(iexsolhir(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
! 					(rec_ihexl(1,ielem_indexi(l,i),i),1:nof_variables))


                     nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
                     lf=rec_ihexl(1,ielem_indexi(L,i),i)
                     rowf=halo_offset(nf) + lf - 1
                     veigr(1:nof_variables)=solhir(rowf,1:nof_variables)

				    end if

				end if



 end if

else

l=facex
i=iconsidered


 if (ielem_interior(i).eq.0)then
     veigr(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i));

 else

 if (ielem_ineighb(l,i).eq.n)then	!my cpu only
				      if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
					  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
					  veigr(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
					  idummy=1
					  else
					  !not periodic ones in my cpu

					   call coordinates_face_inner2dx(n,iconsidered,facex,vext,nodes_list)
					   n_node=2
				  cords(1:2)=zero
 				  call cordinates2(n,nodes_list,n_node,cords(1:2))


				  pox(1)=cords(1)
                  poy(1)=cords(2)

 				  leftv(1:nof_variables)=veigl(1:nof_variables)
				  b_code=ibound_icode(ielem_ibounds(l,i))
 				  call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)

				  veigr(1:nof_variables)=rightv(1:nof_variables)
				      	  end if
				      else
				      !fluid neighbour
				      veigr(1:nof_variables)=u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
				      end if
				else
			      !other my cpu
				    if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
					  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
! 					  veigr(1:nof_variables)=(iexsolhir(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
! 					(rec_ihexl(1,ielem_indexi(l,i),i),1:nof_variables))
                     nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
                     lf=rec_ihexl(1,ielem_indexi(L,i),i)
                     rowf=halo_offset(nf) + lf - 1
                     veigr(1:nof_variables)=solhir(rowf,1:nof_variables)


					  idummy=1
					  end if
				    else

! 				      veigr(1:nof_variables)=(iexsolhir(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
! 					(rec_ihexl(1,ielem_indexi(l,i),i),1:nof_variables))

                     nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
                     lf=rec_ihexl(1,ielem_indexi(L,i),i)
                     rowf=halo_offset(nf) + lf - 1
                     veigr(1:nof_variables)=solhir(rowf,1:nof_variables)

				    end if

				end if

 end if


end if

end subroutine weno_neighbour
















subroutine find_bounds(iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
implicit none
!> @brief
!> Compute local bounds and variation measures for MUSCL limiting WITHOUT storing neighbor states.
!> Fills:
!>   utmin, utmax  : component-wise min/max over stencil
!>   aver_vars     : component-wise average over stencil
!>   sumvars       : sum of abs differences to reference state (first stencil entry)
!>   maxvars       : max absolute value over stencil
#ifdef gpu
!$omp declare target
#endif
integer,intent(in) :: iconsidered
real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout) :: maxvars,aver_vars,sumvars,utmin,utmax

integer :: i,l,iq,k,j,nf,lf,rowf
integer:: nvtot
real :: mp_pinfl,gammal
real,dimension(1:gpu_max_nvar) :: leftv
real,dimension(1:gpu_max_nvar_total) :: uvec, uref

i = iconsidered
nvtot=nof_variables+turbulenceequations+passivescalar

! init outputs
aver_vars = 0.0d0
sumvars  = 0.0d0
maxvars  = 0.0d0
utmin    = huge(1.0)
utmax    = -huge(1.0)

k = 0

! -----------------------
! Build stencil on-the-fly
! -----------------------

uref(1:nof_variables)=u_c_val(1,1:nof_variables,i)
if ((turbulence.eq.1).or.(passivescalar.gt.0))then
uref(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,i)
end if

if (extended_bounds.eq.0) then
  ! strict bounds: center + (some) face neighbors

  uvec = zero
  uvec(1:nof_variables) = u_c_val(1,1:nof_variables,i)

  if (turbulenceequations.ge.1) then
    uvec(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,i)
  end if
  call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)

  if (ielem_interior(i).eq.0) then
    do l=1,ielem_ifca(i)
      uvec(1:nof_variables) = u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
      if (turbulenceequations.ge.1) then
        uvec(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
      end if
      call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
    end do

  else
    do l=1,ielem_ifca(i)

      if (ielem_ineighb(l,i).eq.n) then
        ! neighbor on my cpu
        if (ielem_ibounds(l,i).gt.0) then
          ! boundary: only periodic contributes
          if (ibound_icode(ielem_ibounds(l,i)).eq.5) then
            uvec(1:nof_variables) = u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
            if (turbulenceequations.ge.1) then
              uvec(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
            end if
            call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
          end if
        else
          ! fluid neighbor
          uvec(1:nof_variables) = u_c_val(1,1:nof_variables,ielem_ineigh(l,i))
          if (turbulenceequations.ge.1) then
            uvec(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,ielem_ineigh(l,i))
          end if
          call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
        end if

      else
        ! neighbor on other cpu: only periodic or mpi neighbor (halo)
        if (ielem_ibounds(l,i).gt.0) then
          if (ibound_icode(ielem_ibounds(l,i)).eq.5) then
            nf   = rec_ihexn(1,ielem_indexi(l,i),rec_local(i))
            lf   = rec_ihexl(1,ielem_indexi(l,i),i)
            rowf = halo_offset(nf) + lf - 1
            uvec(1:nof_variables) = solhir(rowf,1:nof_variables)
            if (turbulenceequations.ge.1) then
              uvec(nof_variables+1:nvtot) = solhir(rowf,nof_variables+1:nvtot)
            end if
            call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
          end if
        else
          nf   = rec_ihexn(1,ielem_indexi(l,i),rec_local(i))
          lf   = rec_ihexl(1,ielem_indexi(l,i),i)
          rowf = halo_offset(nf) + lf - 1
          uvec(1:nof_variables) = solhir(rowf,1:nof_variables)
          if (turbulenceequations.ge.1) then
            uvec(nof_variables+1:nvtot) = solhir(rowf,nof_variables+1:nvtot)
          end if
          call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
        end if
      end if

    end do
  end if

else
  ! extended bounds: use rec_ihexl / rec_ihexb neighbor lists
  do iq=1,ielem_inumneighbours(i)

    if (rec_local(i).eq.0) then
      uvec(1:nof_variables) = u_c_val(1,1:nof_variables,rec_ihexl(1,iq,i))
      if (turbulenceequations.ge.1) then
        uvec(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,iq,i))
      end if
      call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)

    else
      if (rec_ihexb(1,iq,rec_local(i)).eq.n) then
        uvec(1:nof_variables) = u_c_val(1,1:nof_variables,rec_ihexl(1,iq,i))
        if (turbulenceequations.ge.1) then
          uvec(nof_variables+1:nvtot) = u_ct_val(1,1:turbulenceequations+passivescalar,rec_ihexl(1,iq,i))
        end if
        call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
      else
        nf   = rec_ihexn(1,iq,i)
        lf   = rec_ihexl(1,iq,i)
        rowf = halo_offset(nf) + lf - 1
        uvec(1:nof_variables) = solhir(rowf,1:nof_variables)
        if (turbulenceequations.ge.1) then
          uvec(nof_variables+1:nvtot) = solhir(rowf,nof_variables+1:nvtot)
        end if
        call process_state(uvec,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
      end if
    end if

  end do
end if

if (k > 0) aver_vars = aver_vars / real(k)
end subroutine find_bounds




subroutine compute_muscl_reconstruction(iconsidered,utmin,utmax)
  implicit none
integer :: kbasis
!> @brief
!> subroutine for computing limited MUSCL reconstructed solution (on-the-fly per GP; minimal storage)
#ifdef gpu
!$omp declare target
#endif
  integer,intent(in)::iconsidered
  integer:: nvtot
  real,intent(in) :: utmin(1:nof_variables+turbulenceequations+passivescalar),utmax(1:nof_variables+turbulenceequations+passivescalar)
  integer :: i,l,ngp,iex,iqp,k,ideg
  real :: ax,ay,az,mp_pinfl,gammal,limvbg,rat,u0,slope_species
  real,dimension(1:gpu_max_nvar) :: leftv,slope_d,candid_lim,delta
  real :: phi(1:gpu_max_dof)
  real :: slope(1:gpu_max_nvar_total)
  real :: psi_min(1:gpu_max_nvar_total)
  real :: usol1(1:gpu_max_nvar_total)
  real :: psi1 (1:gpu_max_nvar_total)
  real :: du   (1:gpu_max_nvar_total)
  nvtot=nof_variables+turbulenceequations+passivescalar
  i = iconsidered
  ideg = idegfree

  do l=1,ielem_ifca(i)
    if (dimensiona.eq.3) then
      if (ielem_types_faces(l,i).eq.5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if
    do ngp=1,iqp
      do iex=1,nof_variables
        rec_uleft(iex,l,ngp,i) = zero
      end do
      if (turbulenceequations.ge.1) then
        do iex=1,turbulenceequations+passivescalar
          rec_uleftturb(iex,l,ngp,i) = zero
        end do
      end if
    end do
  end do

  psi_min(:) = tolbig

  ! --- compute psi_min directly (no storage of usol/psi over faces/GP) ---
  do l=1,ielem_ifca(i)

    if (dimensiona.eq.3) then
      if (ielem_types_faces(l,i).eq.5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if

    do ngp=1,iqp
      ax = rec_qpoints(l,ngp,1,i)
      ay = rec_qpoints(l,ngp,2,i)
      if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)

      if (dimensiona.eq.3) then
          do kbasis=1,ideg
            phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,0,kbasis)
          end do
      else
          do kbasis=1,ideg
            phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg,0,kbasis)
          end do
      end if

      ! unlimited reconstructed state at this GP (stored only in usol1(:))
      delta(1:nof_variables) = zero
      do k=1,ideg
        delta(1:nof_variables) = delta(1:nof_variables) + phi(k) * rec_gradients(1,k,1:nof_variables,i)
      end do

      if (wenwrt.eq.3) then
        leftv(1:nof_variables) = u_c_val(1,1:nof_variables,i)
        call cons2prim(n,leftv,mp_pinfl,gammal)
        usol1(1:nof_variables) = leftv(1:nof_variables) + delta(1:nof_variables)
      else
        usol1(1:nof_variables) = u_c_val(1,1:nof_variables,i) + delta(1:nof_variables)
      end if

      if (turbulenceequations.ge.1) then
        usol1(nof_variables+1:NVTOT) = u_ct_val(1,1:turbulenceequations+passivescalar,i)
        do k=1,ideg
          usol1(nof_variables+1:NVTOT) = usol1(nof_variables+1:NVTOT) + &
            phi(k) * rec_gradients2(1,k,1:turbulenceequations+passivescalar,i)
        end do
      end if
      ! compute psi(:) at this GP (new 1D limiter interface) and accumulate minimum
      call slope_limiters(n,i,utmin,utmax,usol1,psi1)
      psi_min(:) = min(psi_min(:), psi1(:))





    end do
  end do

  slope(:) = psi_min(:)

  if (realgas.eq.1)then
  slope_species = 1.0d0

        do l = 1, nof_species
          iex = dimensiona + 3 + l
          slope_species = min(slope_species, slope(iex))
        end do

        do l = 1, nof_species
          iex = dimensiona + 3 + l
          slope(iex) = slope(iex) * slope_species
        end do

  end if




  ielem_wcx(i) = slope(1)
  if (dg.eq.1) then
    do iex = 1, nvtot



      dg2fv(1:ideg,iex,i) = rec_gradients(1,ideg,iex,i) * slope(iex)
    end do
  end if

!   ! --- extra positivity limiter branch kept as in original ---
!   if (limiter.eq.1) then
!     if (realgas.eq.1) then
!       slope_d(:) = 1.0d0
!
!       do l=1,ielem_ifca(i)
!         if (dimensiona.eq.3) then
!           if (ielem_types_faces(l,i).eq.5) then
!             iqp = qp_quad
!           else
!             iqp = qp_triangle
!           end if
!         else
!           iqp = qp_line
!         end if
!
!         do ngp=1,iqp
!           ! recompute unlimited state at this GP (no stored usol)
!           ax = rec_qpoints(l,ngp,1,i)
!           ay = rec_qpoints(l,ngp,2,i)
!           if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)
!
!           if (dimensiona.eq.3) then
!             phi(1:ideg) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,0)
!           else
!             phi(1:ideg) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg,0)
!           end if
!
!           delta(:) = zero
!           do k=1,ideg
!             delta(:) = delta(:) + phi(k) * rec_gradients(1,k,1:nof_variables,i)
!           end do
!
!           candid_lim(1:nof_variables) = u_c_val(1,1:nof_variables,i) + (delta(:) * slope(1:nof_variables))
!
!           do iex=1,nof_variables
!             if ((iex.ge.2).and.(iex.lt.dimensiona+3)) cycle
!
!             if (u_c_val(1,iex,i) .gt. zero) then
!               if (candid_lim(iex) .lt. zero) then
!                 rat = u_c_val(1,iex,i) / (u_c_val(1,iex,i) - candid_lim(iex))
!                 rat = max(0.0d0, min(1.0d0, rat))
!                 slope_d(iex) = min(slope_d(iex), rat)
!               end if
!             else
!               slope_d(iex) = zero
!             end if
!           end do
!         end do
!       end do
!
!       do iex=1,nof_variables
!         if ((iex.ge.2).and.(iex.lt.dimensiona+3)) cycle
!         slope(iex) = slope(iex) * slope_d(iex)
!       end do
!     end if
!   end if

  ! --- apply slopes: convert usol to increments and extrapolate at faces ---
  do l=1,ielem_ifca(i)

    if (dimensiona.eq.3) then
      if (ielem_types_faces(l,i).eq.5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if

    do ngp=1,iqp
      ! recompute du at this GP (unlimited increment relative to cell average)
      ax = rec_qpoints(l,ngp,1,i)
      ay = rec_qpoints(l,ngp,2,i)
      if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)

      if (dimensiona.eq.3) then
          do kbasis=1,ideg
            phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,0,kbasis)
          end do
      else
          do kbasis=1,ideg
            phi(kbasis) = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg,0,kbasis)
          end do
      end if

      delta(1:nof_variables) = zero
      do k=1,ideg
        delta(1:nof_variables) = delta(1:nof_variables) + phi(k) * rec_gradients(1,k,1:nof_variables,i)
      end do

      if (wenwrt.eq.3) then
        ! du in primitive variables (consistent with extrapolate_bound_muscl_point)
        du(1:nof_variables) = delta(1:nof_variables)
      else
        du(1:nof_variables) = delta(1:nof_variables)
      end if

      if (turbulenceequations.ge.1) then
        du(nof_variables+1:NVTOT) = zero
        do k=1,ideg
          du(nof_variables+1:NVTOT) = du(nof_variables+1:NVTOT) + &
            phi(k) * rec_gradients2(1,k,1:turbulenceequations+passivescalar,i)
        end do
      end if

      call extrapolate_bound_muscl(du,l,ngp,i,slope)
    end do
  end do

end subroutine compute_muscl_reconstruction





subroutine muscl(n)
  implicit none
!> @brief
!> subroutine for muscl type reconstruction
  integer,intent(in)::n
  integer :: i, ii, iconsidered
  real,dimension(1:gpu_max_nvar_total) :: maxvars, aver_vars, sumvars, utmin, utmax


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i,iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
#else
!$omp do
#endif
  do ii=1,nof_interior
    i=el_int(ii)
    iconsidered=i
     if (((ielem_troubled(i).eq.1).and.(ielem_reduce(i).eq.1)).or.((ielem_full(i).eq.0).and.(ielem_troubled(i).eq.1)))then

        if (mood.eq.1)then
              if (ielem_mood(i).gt.0)then
                  if (adda.eq.1) call adda_filter(n,iconsidered)

                call find_bounds(iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
                call compute_muscl_reconstruction(iconsidered,utmin,utmax)
              end if


        else

            if (adda.eq.1) call adda_filter(n,iconsidered)

              call find_bounds(iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
              call compute_muscl_reconstruction(iconsidered,utmin,utmax)


        end if


    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i,iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
#else
!$omp do
#endif
  do ii=1,nof_bounded
    i=el_bnd(ii)
    iconsidered=i



     if (((ielem_troubled(i).eq.1).and.(ielem_reduce(i).eq.1)).or.((ielem_full(i).eq.0).and.(ielem_troubled(i).eq.1)))then


      if (mood.eq.1)then
              if (ielem_mood(i).gt.0)then
                  if (adda.eq.1) call adda_filter(n,iconsidered)

                call find_bounds(iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
                call compute_muscl_reconstruction(iconsidered,utmin,utmax)
              end if


        else

          if (adda.eq.1) call adda_filter(n,iconsidered)

            call find_bounds(iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)
            call compute_muscl_reconstruction(iconsidered,utmin,utmax)


        end if




     end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine muscl


subroutine muscl_acc(n)
  implicit none
!> @brief
!> Low-pressure ideal-gas MUSCL reconstruction for the XPU path.
  integer,intent(in)::n
  integer :: i,kmaxe,l,ngp,iqp,ideg,iex,iq,nf,lf,rowf
  integer :: degree,zpow,ypow,xpow,kk,ktail,ip
  integer :: active_cell
  logical :: include_neighbor
  real :: ax,ay,az,phik,raw_basis,volinv,xterm,yterm,zterm
  real :: rho,oodensity,uu,vv,ww,skinx,skin1,ie1
  real :: d2,sfd,dplus,dmin,epsi2,delu,delu2,edge,edge3
  real :: psi,psi2,sig_1,y_fun,y2,y3,s_y,pol_mog
  real :: du(1:gpu_max_nvar)
  real :: usol(1:gpu_max_nvar)
  real :: uwork(1:gpu_max_nvar)
  real :: utmin(1:gpu_max_nvar)
  real :: utmax(1:gpu_max_nvar)
  real :: psi_min(1:gpu_max_nvar)

  kmaxe=xmpielrank(n)
#ifdef xpu
  !$omp target teams distribute parallel do &
  !$omp& firstprivate(nof_variables, dimensiona, dg, ees, wenwrt, poly, limiter) &
  !$omp& firstprivate(idegfree, qp_quad, qp_triangle, qp_line, zero, gamma, oo2, tolsmall) &
  !$omp& firstprivate(mood, extended_bounds, n, kmaxe) &
  !$omp& map(alloc: xmpielrank, u_c_val, solhir, rec_uleft, rec_gradients, dg2fv) &
  !$omp& map(alloc: ielem_ifca, ielem_types_faces, ielem_iorder, rec_qpoints, rec_volume, integ_basis_value) &
  !$omp& map(alloc: ielem_troubled, ielem_reduce, ielem_full, ielem_hybrid, ielem_mood, ielem_wcx) &
  !$omp& map(alloc: ielem_interior, ielem_ineigh, ielem_ineighb, ielem_ibounds, ibound_icode) &
  !$omp& map(alloc: rec_local, rec_ihexl, rec_ihexb, rec_ihexn, ielem_indexi, ielem_inumneighbours, halo_offset) &
  !$omp& map(alloc: ielem_minedge, ielem_idegfree) &
  !$omp& private(i,l,ngp,iqp,ideg,iex,iq,nf,lf,rowf,degree,zpow,ypow,xpow,kk,ktail,ip) &
  !$omp& private(active_cell,ax,ay,az,phik,raw_basis,volinv,xterm,yterm,zterm) &
  !$omp& private(rho,oodensity,uu,vv,ww,skinx,skin1,ie1,d2,sfd,dplus,dmin,epsi2,delu,delu2) &
  !$omp& private(edge,edge3,psi,psi2,sig_1,y_fun,y2,y3,s_y,pol_mog,du,usol,uwork,utmin,utmax,psi_min)
#else
  !$omp do private(i,l,ngp,iqp,ideg,iex,iq,nf,lf,rowf,degree,zpow,ypow,xpow,kk,ktail,ip) &
  !$omp& private(active_cell,ax,ay,az,phik,raw_basis,volinv,xterm,yterm,zterm) &
  !$omp& private(rho,oodensity,uu,vv,ww,skinx,skin1,ie1,d2,sfd,dplus,dmin,epsi2,delu,delu2) &
  !$omp& private(edge,edge3,psi,psi2,sig_1,y_fun,y2,y3,s_y,pol_mog,du,usol,uwork,utmin,utmax,psi_min)
#endif
  do i=1,kmaxe
    active_cell=0
    if (((ielem_troubled(i).eq.1).and.(ielem_reduce(i).eq.1)).or. &
        ((ielem_full(i).eq.0).and.(ielem_troubled(i).eq.1)))then
      active_cell=1
      if (mood.eq.1)then
        if (ielem_mood(i).le.0) active_cell=0
      end if
    end if

    if (active_cell.eq.1)then
      ideg=idegfree

      do iex=1,nof_variables
        utmin(iex)=1.0d300
        utmax(iex)=-1.0d300
        psi_min(iex)=1.0d300
      end do

      do iex=1,nof_variables
        uwork(iex)=u_c_val(1,iex,i)
      end do
      if (wenwrt.eq.3)then
        oodensity=1.0d0/uwork(1)
        uu=uwork(2)*oodensity
        vv=uwork(3)*oodensity
        ww=zero
        if (dimensiona.eq.3) ww=uwork(4)*oodensity
        skinx=uu*uu+vv*vv+ww*ww
        uwork(2)=uu
        uwork(3)=vv
        if (dimensiona.eq.3) uwork(4)=ww
        uwork(dimensiona+2)=(gamma-1.0d0)*(uwork(dimensiona+2)-oo2*uwork(1)*skinx)
      end if
      do iex=1,nof_variables
        utmin(iex)=min(utmin(iex),uwork(iex))
        utmax(iex)=max(utmax(iex),uwork(iex))
      end do

      if (extended_bounds.eq.0)then
        if (ielem_interior(i).eq.0)then
          do l=1,ielem_ifca(i)
            do iex=1,nof_variables
              uwork(iex)=u_c_val(1,iex,ielem_ineigh(l,i))
            end do
            if (wenwrt.eq.3)then
              oodensity=1.0d0/uwork(1)
              uu=uwork(2)*oodensity
              vv=uwork(3)*oodensity
              ww=zero
              if (dimensiona.eq.3) ww=uwork(4)*oodensity
              skinx=uu*uu+vv*vv+ww*ww
              uwork(2)=uu
              uwork(3)=vv
              if (dimensiona.eq.3) uwork(4)=ww
              uwork(dimensiona+2)=(gamma-1.0d0)*(uwork(dimensiona+2)-oo2*uwork(1)*skinx)
            end if
            do iex=1,nof_variables
              utmin(iex)=min(utmin(iex),uwork(iex))
              utmax(iex)=max(utmax(iex),uwork(iex))
            end do
          end do
        else
          do l=1,ielem_ifca(i)
            if (ielem_ineighb(l,i).eq.n)then
              include_neighbor=.false.
              if (ielem_ibounds(l,i).le.0)then
                include_neighbor=.true.
              else
                if (ibound_icode(ielem_ibounds(l,i)).eq.5) include_neighbor=.true.
              end if
              if (include_neighbor)then
                do iex=1,nof_variables
                  uwork(iex)=u_c_val(1,iex,ielem_ineigh(l,i))
                end do
                if (wenwrt.eq.3)then
                  oodensity=1.0d0/uwork(1)
                  uu=uwork(2)*oodensity
                  vv=uwork(3)*oodensity
                  ww=zero
                  if (dimensiona.eq.3) ww=uwork(4)*oodensity
                  skinx=uu*uu+vv*vv+ww*ww
                  uwork(2)=uu
                  uwork(3)=vv
                  if (dimensiona.eq.3) uwork(4)=ww
                  uwork(dimensiona+2)=(gamma-1.0d0)*(uwork(dimensiona+2)-oo2*uwork(1)*skinx)
                end if
                do iex=1,nof_variables
                  utmin(iex)=min(utmin(iex),uwork(iex))
                  utmax(iex)=max(utmax(iex),uwork(iex))
                end do
              end if
            else
              include_neighbor=.false.
              if (ielem_ibounds(l,i).le.0)then
                include_neighbor=.true.
              else
                if (ibound_icode(ielem_ibounds(l,i)).eq.5) include_neighbor=.true.
              end if
              if (include_neighbor)then
                nf=rec_ihexn(1,ielem_indexi(l,i),rec_local(i))
                lf=rec_ihexl(1,ielem_indexi(l,i),i)
                rowf=halo_offset(nf)+lf-1
                do iex=1,nof_variables
                  uwork(iex)=solhir(rowf,iex)
                end do
                if (wenwrt.eq.3)then
                  oodensity=1.0d0/uwork(1)
                  uu=uwork(2)*oodensity
                  vv=uwork(3)*oodensity
                  ww=zero
                  if (dimensiona.eq.3) ww=uwork(4)*oodensity
                  skinx=uu*uu+vv*vv+ww*ww
                  uwork(2)=uu
                  uwork(3)=vv
                  if (dimensiona.eq.3) uwork(4)=ww
                  uwork(dimensiona+2)=(gamma-1.0d0)*(uwork(dimensiona+2)-oo2*uwork(1)*skinx)
                end if
                do iex=1,nof_variables
                  utmin(iex)=min(utmin(iex),uwork(iex))
                  utmax(iex)=max(utmax(iex),uwork(iex))
                end do
              end if
            end if
          end do
        end if
      else
        do iq=1,ielem_inumneighbours(i)
          if ((rec_local(i).eq.0).or.(rec_ihexb(1,iq,rec_local(i)).eq.n))then
            do iex=1,nof_variables
              uwork(iex)=u_c_val(1,iex,rec_ihexl(1,iq,i))
            end do
          else
            nf=rec_ihexn(1,iq,i)
            lf=rec_ihexl(1,iq,i)
            rowf=halo_offset(nf)+lf-1
            do iex=1,nof_variables
              uwork(iex)=solhir(rowf,iex)
            end do
          end if
          if (wenwrt.eq.3)then
            oodensity=1.0d0/uwork(1)
            uu=uwork(2)*oodensity
            vv=uwork(3)*oodensity
            ww=zero
            if (dimensiona.eq.3) ww=uwork(4)*oodensity
            skinx=uu*uu+vv*vv+ww*ww
            uwork(2)=uu
            uwork(3)=vv
            if (dimensiona.eq.3) uwork(4)=ww
            uwork(dimensiona+2)=(gamma-1.0d0)*(uwork(dimensiona+2)-oo2*uwork(1)*skinx)
          end if
          do iex=1,nof_variables
            utmin(iex)=min(utmin(iex),uwork(iex))
            utmax(iex)=max(utmax(iex),uwork(iex))
          end do
        end do
      end if

      do l=1,ielem_ifca(i)
        if (dimensiona.eq.3)then
          if (ielem_types_faces(l,i).eq.5)then
            iqp=qp_quad
          else
            iqp=qp_triangle
          end if
        else
          iqp=qp_line
        end if
        do ngp=1,iqp
          ax=rec_qpoints(l,ngp,1,i)
          ay=rec_qpoints(l,ngp,2,i)
          az=zero
          if (dimensiona.eq.3) az=rec_qpoints(l,ngp,3,i)
          volinv=1.0d0/rec_volume(1,1,i)

          do iex=1,nof_variables
            du(iex)=zero
          end do

          kk=0
          if (dimensiona.eq.3)then
            do degree=1,ielem_iorder(i)
              do zpow=0,degree
                zterm=1.0d0
                do ip=1,zpow
                  zterm=zterm*az
                end do
                do ypow=0,degree-zpow
                  xpow=degree-zpow-ypow
                  kk=kk+1
                  if (kk.le.ideg)then
                    xterm=1.0d0
                    do ip=1,xpow
                      xterm=xterm*ax
                    end do
                    yterm=1.0d0
                    do ip=1,ypow
                      yterm=yterm*ay
                    end do
                    raw_basis=xterm*yterm*zterm
                    phik=raw_basis-integ_basis_value(kk,i)*volinv
                    do iex=1,nof_variables
                      du(iex)=du(iex)+phik*rec_gradients(1,kk,iex,i)
                    end do
                  end if
                end do
              end do
            end do
          else
            do degree=1,ielem_iorder(i)
              do ypow=0,degree
                xpow=degree-ypow
                kk=kk+1
                if (kk.le.ideg)then
                  xterm=1.0d0
                  do ip=1,xpow
                    xterm=xterm*ax
                  end do
                  yterm=1.0d0
                  do ip=1,ypow
                    yterm=yterm*ay
                  end do
                  raw_basis=xterm*yterm
                  phik=raw_basis-integ_basis_value(kk,i)*volinv
                  do iex=1,nof_variables
                    du(iex)=du(iex)+phik*rec_gradients(1,kk,iex,i)
                  end do
                end if
              end do
            end do
          end if

          do ktail=kk+1,ideg
            phik=-integ_basis_value(ktail,i)*volinv
            do iex=1,nof_variables
              du(iex)=du(iex)+phik*rec_gradients(1,ktail,iex,i)
            end do
          end do

          if (wenwrt.eq.3)then
            do iex=1,nof_variables
              usol(iex)=u_c_val(1,iex,i)
            end do
            oodensity=1.0d0/usol(1)
            uu=usol(2)*oodensity
            vv=usol(3)*oodensity
            ww=zero
            if (dimensiona.eq.3) ww=usol(4)*oodensity
            skinx=uu*uu+vv*vv+ww*ww
            usol(2)=uu
            usol(3)=vv
            if (dimensiona.eq.3) usol(4)=ww
            usol(dimensiona+2)=(gamma-1.0d0)*(usol(dimensiona+2)-oo2*usol(1)*skinx)
            do iex=1,nof_variables
              usol(iex)=usol(iex)+du(iex)
            end do
          else
            do iex=1,nof_variables
              usol(iex)=u_c_val(1,iex,i)+du(iex)
            end do
          end if

          do iex=1,nof_variables
            if (wenwrt.eq.3)then
              uwork(iex)=u_c_val(1,iex,i)
              if (iex.eq.2) uwork(iex)=u_c_val(1,iex,i)/u_c_val(1,1,i)
              if (iex.eq.3) uwork(iex)=u_c_val(1,iex,i)/u_c_val(1,1,i)
              if ((dimensiona.eq.3).and.(iex.eq.4)) uwork(iex)=u_c_val(1,iex,i)/u_c_val(1,1,i)
              if (iex.eq.dimensiona+2)then
                rho=u_c_val(1,1,i)
                uu=u_c_val(1,2,i)/rho
                vv=u_c_val(1,3,i)/rho
                ww=zero
                if (dimensiona.eq.3) ww=u_c_val(1,4,i)/rho
                uwork(iex)=(gamma-1.0d0)*(u_c_val(1,iex,i)-oo2*rho*(uu*uu+vv*vv+ww*ww))
              end if
            else
              uwork(iex)=u_c_val(1,iex,i)
            end if

            d2=usol(iex)-uwork(iex)
            if (abs(d2).le.zero)then
              psi=1.0d0
            else
              if (d2.gt.zero)then
                sfd=(utmax(iex)-uwork(iex))/d2
                dplus=utmax(iex)-uwork(iex)
              else
                sfd=(utmin(iex)-uwork(iex))/d2
                dplus=utmin(iex)-uwork(iex)
              end if
              select case(limiter)
              case(1)
                psi=min(1.0d0,sfd)
              case(10)
                psi=min(0.3d0,sfd)
              case(2,3,9)
                pol_mog=-((4.0d0/27.0d0)*sfd*sfd*sfd)+sfd
                if (sfd.lt.1.5d0)then
                  psi=pol_mog
                else
                  psi=1.0d0
                end if
                edge=0.1d0*ielem_minedge(i)
                edge3=edge*edge*edge
                delu=utmax(iex)-utmin(iex)
                delu2=delu*delu
                if (delu2.le.edge3)then
                  sig_1=1.0d0
                else if ((edge3.lt.delu2).and.(delu2.lt.2.0d0*edge3))then
                  y_fun=(delu2-edge3)/edge3
                  y2=y_fun*y_fun
                  y3=y2*y_fun
                  s_y=2.0d0*y3-3.0d0*y2+1.0d0
                  sig_1=s_y
                else
                  sig_1=0.0d0
                end if
                psi2=sig_1+(1.0d0-sig_1)*psi
                psi=psi2
              case(4,8)
                dmin=usol(iex)-uwork(iex)
                dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
                edge=0.1d0*ielem_minedge(i)
                epsi2=edge*edge*edge
                psi=(1.0d0/dmin)*((((dplus*dplus)+epsi2)*dmin+(2.0d0*(dmin*dmin)*dplus))/ &
                  ((dplus*dplus)+(2.0d0*dmin*dmin)+(dmin*dplus)+epsi2))
                delu=utmax(iex)-utmin(iex)
                delu2=delu*delu
                if (delu2.le.epsi2)then
                  sig_1=1.0d0
                else if ((epsi2.lt.delu2).and.(delu2.lt.2.0d0*epsi2))then
                  y_fun=(delu2-epsi2)/epsi2
                  y2=y_fun*y_fun
                  y3=y2*y_fun
                  s_y=2.0d0*y3-3.0d0*y2+1.0d0
                  sig_1=s_y
                else
                  sig_1=0.0d0
                end if
                psi2=sig_1+(1.0d0-sig_1)*psi
                psi=psi2
              case(5)
                psi=(sfd*sfd+sfd)/(sfd*sfd+1.0d0)
              case(6)
                psi=2.0d0*sfd/(sfd+1.0d0)
              case(7)
                dmin=usol(iex)-uwork(iex)
                dmin=sign(1.0,dmin)*(abs(dmin)+tolsmall)
                edge=0.1d0*ielem_minedge(i)
                epsi2=edge*edge*edge
                psi=(1.0d0/dmin)*((((dplus*dplus)+epsi2)*dmin+(2.0d0*(dmin*dmin)*dplus))/ &
                  ((dplus*dplus)+(2.0d0*dmin*dmin)+(dmin*dplus)+epsi2))
              case default
                psi=min(1.0d0,sfd)
              end select
            end if
            psi_min(iex)=min(psi_min(iex),psi)
          end do
        end do
      end do

      do iex=1,nof_variables
        if (psi_min(iex).gt.1.0d200) psi_min(iex)=1.0d0
      end do
      ielem_wcx(i)=psi_min(1)

      if (dg.eq.1)then
        do iex=1,nof_variables
          do kk=1,ielem_idegfree(i)
            dg2fv(kk,iex,i)=rec_gradients(1,kk,iex,i)*psi_min(iex)
          end do
        end do
      end if

      do l=1,ielem_ifca(i)
        if (dimensiona.eq.3)then
          if (ielem_types_faces(l,i).eq.5)then
            iqp=qp_quad
          else
            iqp=qp_triangle
          end if
        else
          iqp=qp_line
        end if
        do ngp=1,iqp
          ax=rec_qpoints(l,ngp,1,i)
          ay=rec_qpoints(l,ngp,2,i)
          az=zero
          if (dimensiona.eq.3) az=rec_qpoints(l,ngp,3,i)
          volinv=1.0d0/rec_volume(1,1,i)
          do iex=1,nof_variables
            du(iex)=zero
          end do

          kk=0
          if (dimensiona.eq.3)then
            do degree=1,ielem_iorder(i)
              do zpow=0,degree
                zterm=1.0d0
                do ip=1,zpow
                  zterm=zterm*az
                end do
                do ypow=0,degree-zpow
                  xpow=degree-zpow-ypow
                  kk=kk+1
                  if (kk.le.ideg)then
                    xterm=1.0d0
                    do ip=1,xpow
                      xterm=xterm*ax
                    end do
                    yterm=1.0d0
                    do ip=1,ypow
                      yterm=yterm*ay
                    end do
                    raw_basis=xterm*yterm*zterm
                    phik=raw_basis-integ_basis_value(kk,i)*volinv
                    do iex=1,nof_variables
                      du(iex)=du(iex)+phik*rec_gradients(1,kk,iex,i)
                    end do
                  end if
                end do
              end do
            end do
          else
            do degree=1,ielem_iorder(i)
              do ypow=0,degree
                xpow=degree-ypow
                kk=kk+1
                if (kk.le.ideg)then
                  xterm=1.0d0
                  do ip=1,xpow
                    xterm=xterm*ax
                  end do
                  yterm=1.0d0
                  do ip=1,ypow
                    yterm=yterm*ay
                  end do
                  raw_basis=xterm*yterm
                  phik=raw_basis-integ_basis_value(kk,i)*volinv
                  do iex=1,nof_variables
                    du(iex)=du(iex)+phik*rec_gradients(1,kk,iex,i)
                  end do
                end if
              end do
            end do
          end if

          do ktail=kk+1,ideg
            phik=-integ_basis_value(ktail,i)*volinv
            do iex=1,nof_variables
              du(iex)=du(iex)+phik*rec_gradients(1,ktail,iex,i)
            end do
          end do

          if (wenwrt.eq.3)then
            do iex=1,nof_variables
              usol(iex)=u_c_val(1,iex,i)
            end do
            oodensity=1.0d0/usol(1)
            uu=usol(2)*oodensity
            vv=usol(3)*oodensity
            ww=zero
            if (dimensiona.eq.3) ww=usol(4)*oodensity
            skinx=uu*uu+vv*vv+ww*ww
            usol(2)=uu+du(2)*psi_min(2)
            usol(3)=vv+du(3)*psi_min(3)
            if (dimensiona.eq.3) usol(4)=ww+du(4)*psi_min(4)
            usol(1)=usol(1)+du(1)*psi_min(1)
            usol(dimensiona+2)=(gamma-1.0d0)*(u_c_val(1,dimensiona+2,i)-oo2*u_c_val(1,1,i)*skinx) + &
              du(dimensiona+2)*psi_min(dimensiona+2)
            uu=usol(2)
            vv=usol(3)
            ww=zero
            if (dimensiona.eq.3) ww=usol(4)
            skin1=oo2*(uu*uu+vv*vv+ww*ww)
            ie1=usol(dimensiona+2)/((gamma-1.0d0)*usol(1))
            rec_uleft(1,l,ngp,i)=usol(1)
            rec_uleft(2,l,ngp,i)=usol(1)*uu
            rec_uleft(3,l,ngp,i)=usol(1)*vv
            if (dimensiona.eq.3) rec_uleft(4,l,ngp,i)=usol(1)*ww
            rec_uleft(dimensiona+2,l,ngp,i)=usol(1)*(ie1+skin1)
          else
            do iex=1,nof_variables
              rec_uleft(iex,l,ngp,i)=u_c_val(1,iex,i)+du(iex)*psi_min(iex)
            end do
          end if
        end do
      end do
    end if
  end do
#ifdef xpu
  !$omp end target teams distribute parallel do
#else
  !$omp end do
#endif

end subroutine muscl_acc





subroutine solutiontriav2(n)
  implicit none
integer :: mm_i
integer :: mm_j
real :: mm_old(1:gpu_max_dim)
real :: mm_sum
!> @brief
!> subroutine for extrapolating the unlimited reconstructed values for diffusive fluxes in 3d
  integer,intent(in) :: n

  integer :: i,l,iex,nvar,ngp
  integer :: k,ideg,nfaces,ggs,iqp
  integer :: ihgt,ihgj
  real    :: ax,ay,az
  real,dimension(1:gpu_max_dim) :: ugradloc
  real,dimension(1:gpu_max_dim,1:gpu_max_dim) :: ainvjt
  real,dimension(1:gpu_max_dim,1:gpu_max_dim) :: vextc
  real :: dx,dy,dz,gk
  real :: mp_pinfl,gammal
  real,dimension(1:gpu_max_nvar) :: leftv

  integer :: kmaxe

  kmaxe=xmpielrank(n)



#ifdef xpu
  !$omp target teams distribute parallel do &
  !$omp& firstprivate(n) &
  !$omp& firstprivate(nof_variables, dimensiona, turbulence, passivescalar, turbulenceequations, icoupleturb) &
  !$omp& firstprivate(poly, zero, qp_quad, qp_triangle, qp_line) &
  !$omp& map(alloc: xmpielrank) &
  !$omp& map(alloc: ielem_ifca, ielem_types_faces, ielem_idegfree, ielem_ggs, ielem_totvolume) &
  !$omp& map(alloc: rec_qpoints, rec_invccjac, rec_vext_ref, ielem_xxc, ielem_yyc, ielem_zzc) &
  !$omp& map(alloc: rec_gradf, rec_grads, rec_uleftv) &
  !$omp& private(i,l,iex,nvar,ngp,k,ideg,nfaces,ggs,iqp,ihgt,ihgj) &
  !$omp& private(ax,ay,az,ugradloc,ainvjt,vextc,dx,dy,dz,gk,mp_pinfl,gammal,leftv) &
  !$omp& private(mm_i,mm_j,mm_old,mm_sum)
#elif defined(gpu)
  !$omp target teams distribute parallel do &
  !$omp& private(i,l,iex,nvar,ngp, k,ideg,nfaces,ggs,iqp, ihgt,ihgj, ax,ay,az,ugradloc,ainvjt,vextc, dx,dy,dz,gk,mp_pinfl, gammal, leftv) &
  !$omp& private(mm_i,mm_j,mm_old,mm_sum)
#else
  !$omp do
#endif
  do i = 1, xmpielrank(n)

    nfaces = ielem_ifca(i)
    do l = 1, nfaces
      if (dimensiona.eq.3) then
        if (ielem_types_faces(l,i).eq.5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if
      do ngp = 1, iqp
        do iex = 1, nof_variables-1
          do ihgt = 1, dimensiona
            rec_uleftv(ihgt,iex,l,ngp,i) = zero
          end do
        end do
#ifndef xpu
        if ((turbulence.gt.0).or.(passivescalar.gt.0)) then
          do nvar = 1, turbulenceequations+passivescalar
            do ihgt = 1, dimensiona
              rec_uleftturbv(ihgt,nvar,l,ngp,i) = zero
            end do
          end do
        end if
#endif
      end do
    end do

    ! inv Jacobian transpose (your original transpose fill)
    do ihgt=1,dimensiona
      do ihgj=1,dimensiona
        ainvjt(ihgt,ihgj) = rec_invccjac(ihgj,ihgt,i)
      end do
    end do

    ideg   = ielem_idegfree(i)
    nfaces = ielem_ifca(i)
    ggs    = ielem_ggs(i)

    do l = 1, nfaces

      if (dimensiona.eq.3) then
        if (ielem_types_faces(l,i).eq.5) then
          iqp = qp_quad
        else
          iqp = qp_triangle
        end if
      else
        iqp = qp_line
      end if

      select case (ggs)

      case (0)
        ! ------------------------------------------
        ! GGS=0 : compute gradients using basis derivs
        ! ------------------------------------------
        do ngp = 1, iqp
          ax = rec_qpoints(l,ngp,1,i)
          ay = rec_qpoints(l,ngp,2,i)
          if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)

          ! turbulence/passive first
#ifndef xpu
          if ((turbulence.gt.0).or.(passivescalar.gt.0)) then

            if (icoupleturb.eq.0) then
              do nvar=1,turbulenceequations+passivescalar
                rec_uleftturb(nvar,l,ngp,i) = u_ct_val(1,nvar,i)
              end do
            end if

            do nvar=1,turbulenceequations+passivescalar

              ugradloc = zero

              if (dimensiona.eq.3) then
                select case (poly)
                case (1)
                  do k=1,ideg
                    dx = dfx(ax,ay,az,k,i); dy = dfy(ax,ay,az,k,i); dz = dfz(ax,ay,az,k,i)
                    gk = rec_gradientsturb(1,k,nvar,i)
                    ugradloc(1)=ugradloc(1)+gk*dx
                    ugradloc(2)=ugradloc(2)+gk*dy
                    ugradloc(3)=ugradloc(3)+gk*dz
                  end do
                case (2)
                  do k=1,ideg
                    dx = dlx(ax,ay,az,k,i); dy = dly(ax,ay,az,k,i); dz = dlz(ax,ay,az,k,i)
                    gk = rec_gradientsturb(1,k,nvar,i)
                    ugradloc(1)=ugradloc(1)+gk*dx
                    ugradloc(2)=ugradloc(2)+gk*dy
                    ugradloc(3)=ugradloc(3)+gk*dz
                  end do
                case (4)
                  do k=1,ideg
                    dx = tl3dx(ax,ay,az,k,i); dy = tl3dy(ax,ay,az,k,i); dz = tl3dz(ax,ay,az,k,i)
                    gk = rec_gradientsturb(1,k,nvar,i)
                    ugradloc(1)=ugradloc(1)+gk*dx
                    ugradloc(2)=ugradloc(2)+gk*dy
                    ugradloc(3)=ugradloc(3)+gk*dz
                  end do
                end select
              else
                if (poly.eq.4) then
                  do k=1,ideg
                    dx = tl2dx(ax,ay,k,i); dy = tl2dy(ax,ay,k,i)
                    gk = rec_gradientsturb(1,k,nvar,i)
                    ugradloc(1)=ugradloc(1)+gk*dx
                    ugradloc(2)=ugradloc(2)+gk*dy
                  end do
                else
                  do k=1,ideg
                    dx = df2dx(ax,ay,k,i); dy = df2dy(ax,ay,k,i)
                    gk = rec_gradientsturb(1,k,nvar,i)
                    ugradloc(1)=ugradloc(1)+gk*dx
                    ugradloc(2)=ugradloc(2)+gk*dy
                  end do
                end if
              end if

                do mm_i=1,dimensiona
                  mm_sum=zero
                  do mm_j=1,dimensiona
                    mm_sum=mm_sum+ainvjt(mm_i,mm_j)*ugradloc(mm_j)
                  end do
                  rec_uleftturbv(mm_i,nvar,l,ngp,i)=mm_sum
                end do
            end do
          end if
#endif

          ! viscous gradients (mean flow)
          do iex = 1, nof_variables-1

            ugradloc = zero

            if (dimensiona.eq.3) then
              select case (poly)
              case (1)
                do k=1,ideg
                  dx = dfx(ax,ay,az,k,i); dy = dfy(ax,ay,az,k,i); dz = dfz(ax,ay,az,k,i)
                  gk = rec_gradf(iex,k,i)
                  ugradloc(1)=ugradloc(1)+gk*dx
                  ugradloc(2)=ugradloc(2)+gk*dy
                  ugradloc(3)=ugradloc(3)+gk*dz
                end do
              case (2)
                do k=1,ideg
                  dx = dlx(ax,ay,az,k,i); dy = dly(ax,ay,az,k,i); dz = dlz(ax,ay,az,k,i)
                  gk = rec_gradf(iex,k,i)
                  ugradloc(1)=ugradloc(1)+gk*dx
                  ugradloc(2)=ugradloc(2)+gk*dy
                  ugradloc(3)=ugradloc(3)+gk*dz
                end do
              case (4)
                do k=1,ideg
                  dx = tl3dx(ax,ay,az,k,i); dy = tl3dy(ax,ay,az,k,i); dz = tl3dz(ax,ay,az,k,i)
                  gk = rec_gradf(iex,k,i)
                  ugradloc(1)=ugradloc(1)+gk*dx
                  ugradloc(2)=ugradloc(2)+gk*dy
                  ugradloc(3)=ugradloc(3)+gk*dz
                end do
              end select
            else
              if (poly.eq.4) then
                do k=1,ideg
                  dx = tl2dx(ax,ay,k,i); dy = tl2dy(ax,ay,k,i)
                  gk = rec_gradf(iex,k,i)
                  ugradloc(1)=ugradloc(1)+gk*dx
                  ugradloc(2)=ugradloc(2)+gk*dy
                end do
              else
                do k=1,ideg
                  dx = df2dx(ax,ay,k,i); dy = df2dy(ax,ay,k,i)
                  gk = rec_gradf(iex,k,i)
                  ugradloc(1)=ugradloc(1)+gk*dx
                  ugradloc(2)=ugradloc(2)+gk*dy
                end do
              end if
            end if

              do mm_i=1,dimensiona
                mm_sum=zero
                do mm_j=1,dimensiona
                  mm_sum=mm_sum+ainvjt(mm_i,mm_j)*ugradloc(mm_j)
                end do
                rec_uleftv(mm_i,iex,l,ngp,i)=mm_sum
              end do




          end do

        end do  ! ngp

      case (1)
        ! ------------------------------------------
        ! GGS=1 : gradients already stored in rec_grads
        ! ------------------------------------------
        do ngp=1,iqp

#ifndef xpu
          if ((turbulence.gt.0).or.(passivescalar.gt.0)) then
            if (icoupleturb.eq.0) then
              do nvar=1,turbulenceequations+passivescalar
                rec_uleftturb(nvar,l,ngp,i)=u_ct_val(1,nvar,i)
              end do
            end if
            do nvar=1,turbulenceequations+passivescalar
              rec_uleftturbv(1:dimensiona,nvar,l,ngp,i) = rec_grads(dimensiona+1+nvar,1:dimensiona,i)
            end do
          end if
#endif

          do iex=1,nof_variables-1
            rec_uleftv(1:dimensiona,iex,l,ngp,i) = rec_grads(iex,1:dimensiona,i)
          end do

        end do

      end select

    end do  ! faces

    ! ----  final "element center gradient" part, compute derivs on-the-fly ----
    if (ggs.eq.0) then

      vextc(1,1) = ielem_xxc(i)
      vextc(1,2) = ielem_yyc(i)
      if (dimensiona.eq.3) vextc(1,3) = ielem_zzc(i)

        do mm_j=1,dimensiona
          mm_old(mm_j)=vextc(1,mm_j)-rec_vext_ref(mm_j,i)
        end do
        do mm_i=1,dimensiona
          vextc(1,mm_i)=zero
          do mm_j=1,dimensiona
            vextc(1,mm_i)=vextc(1,mm_i)+rec_invccjac(mm_i,mm_j,i)*mm_old(mm_j)
          end do
        end do

      ax = vextc(1,1)
      ay = vextc(1,2)
      if (dimensiona.eq.3) az = vextc(1,3)

      do iex=1,nof_variables-1

        ugradloc = zero

        if (dimensiona.eq.3) then
          select case (poly)
          case (1)
            do k=1,ideg
              dx = dfx(ax,ay,az,k,i); dy = dfy(ax,ay,az,k,i); dz = dfz(ax,ay,az,k,i)
              gk = rec_gradf(iex,k,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
              ugradloc(3)=ugradloc(3)+gk*dz
            end do
          case (2)
            do k=1,ideg
              dx = dlx(ax,ay,az,k,i); dy = dly(ax,ay,az,k,i); dz = dlz(ax,ay,az,k,i)
              gk = rec_gradf(iex,k,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
              ugradloc(3)=ugradloc(3)+gk*dz
            end do
          case (4)
            do k=1,ideg
              dx = tl3dx(ax,ay,az,k,i); dy = tl3dy(ax,ay,az,k,i); dz = tl3dz(ax,ay,az,k,i)
              gk = rec_gradf(iex,k,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
              ugradloc(3)=ugradloc(3)+gk*dz
            end do
          end select
        else
          if (poly.eq.4) then
            do k=1,ideg
              dx = tl2dx(ax,ay,k,i); dy = tl2dy(ax,ay,k,i)
              gk = rec_gradf(iex,k,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
            end do
          else
            do k=1,ideg
              dx = df2dx(ax,ay,k,i); dy = df2dy(ax,ay,k,i)
              gk = rec_gradf(iex,k,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
            end do
          end if
        end if

          do mm_i=1,dimensiona
            mm_sum=zero
            do mm_j=1,dimensiona
              mm_sum=mm_sum+ainvjt(mm_i,mm_j)*ugradloc(mm_j)
            end do
            rec_grads(iex,mm_i,i)=mm_sum*ielem_totvolume(i)
          end do
      end do

#ifndef xpu
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        do iex=1,turbulenceequations+passivescalar

        ugradloc = zero

        if (dimensiona.eq.3) then
          select case (poly)
          case (1)
            do k=1,ideg
              dx = dfx(ax,ay,az,k,i); dy = dfy(ax,ay,az,k,i); dz = dfz(ax,ay,az,k,i)
              gk = rec_gradientsturb(1,k,iex,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
              ugradloc(3)=ugradloc(3)+gk*dz
            end do
          case (2)
            do k=1,ideg
              dx = dlx(ax,ay,az,k,i); dy = dly(ax,ay,az,k,i); dz = dlz(ax,ay,az,k,i)
              gk = rec_gradientsturb(1,k,iex,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
              ugradloc(3)=ugradloc(3)+gk*dz
            end do
          case (4)
            do k=1,ideg
              dx = tl3dx(ax,ay,az,k,i); dy = tl3dy(ax,ay,az,k,i); dz = tl3dz(ax,ay,az,k,i)
              gk = rec_gradientsturb(1,k,iex,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
              ugradloc(3)=ugradloc(3)+gk*dz
            end do
          end select
        else
          if (poly.eq.4) then
            do k=1,ideg
              dx = tl2dx(ax,ay,k,i); dy = tl2dy(ax,ay,k,i)
              gk = rec_gradientsturb(1,k,iex,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
            end do
          else
            do k=1,ideg
              dx = df2dx(ax,ay,k,i); dy = df2dy(ax,ay,k,i)
             gk = rec_gradientsturb(1,k,iex,i)
              ugradloc(1)=ugradloc(1)+gk*dx
              ugradloc(2)=ugradloc(2)+gk*dy
            end do
          end if
        end if

          do mm_i=1,dimensiona
            mm_sum=zero
            do mm_j=1,dimensiona
              mm_sum=mm_sum+ainvjt(mm_i,mm_j)*ugradloc(mm_j)
            end do
            rec_grads(nof_variables-1+iex,mm_i,i)=mm_sum*ielem_totvolume(i)
          end do
      end do


      end if
#endif

    end if

  end do
#ifdef xpu
  !$omp end target teams distribute parallel do
#elif defined(gpu)
  !$omp end target teams distribute parallel do
#else
  !$omp end do
#endif

end subroutine solutiontriav2


subroutine compute_linear_reconstruction_acc2(n)
  implicit none
  integer,intent(in)::n
  integer :: i,kmaxe,l,ngp,iqp,ideg,iex
  integer :: degree,zpow,ypow,xpow,kk,ktail,ip
  real :: ax,ay,az,phik
  real :: raw_basis,volinv,xterm,yterm,zterm
  real :: du(1:gpu_max_nvar_total)

  kmaxe=xmpielrank(n)
#ifdef xpu

  !$omp target teams distribute parallel do &
  !$omp& firstprivate(nof_variables, dimensiona, turbulenceequations, passivescalar, dg, ees, wenwrt, poly) &
  !$omp& firstprivate(idegfree, qp_quad, qp_triangle, qp_line, zero, gamma, oo2, r_gas) &
  !$omp& map(alloc: xmpielrank) &
  !$omp& map(alloc: u_c_val, rec_uleft) &
  !$omp& map(alloc: ielem_ifca, ielem_types_faces, ielem_iorder, rec_qpoints, rec_gradients) &
  !$omp& map(alloc: rec_volume, ielem_totvolume, integ_basis_value) &
  !$omp& firstprivate(n,kmaxe) &
  !$omp& private(i,l,ngp,iqp,ideg,iex,ax,ay,az,phik,du) &
  !$omp& private(degree,zpow,ypow,xpow,kk,ktail,ip,raw_basis,volinv,xterm,yterm,zterm)
#else
  !$omp do private(i,l,ngp,iqp,ideg,iex,ax,ay,az,phik,du) &
  !$omp& private(degree,zpow,ypow,xpow,kk,ktail,ip,raw_basis,volinv,xterm,yterm,zterm)
#endif
  do i=1,kmaxe
       ideg = idegfree

       do l=1,ielem_ifca(i)
         if (dimensiona.eq.3) then
           if (ielem_types_faces(l,i).eq.5) then
             iqp = qp_quad
           else
             iqp = qp_triangle
           end if
         else
           iqp = qp_line
         end if
         do ngp=1,iqp
           do iex=1,nof_variables
             rec_uleft(iex,l,ngp,i)=zero
           end do
         end do
       end do

       do l=1,ielem_ifca(i)

         if (dimensiona.eq.3) then
           if (ielem_types_faces(l,i).eq.5) then
             iqp = qp_quad
           else
             iqp = qp_triangle
           end if
         else
           iqp = qp_line
         end if

         do ngp=1,iqp
           ax = rec_qpoints(l,ngp,1,i)
           ay = rec_qpoints(l,ngp,2,i)
           az = zero
           if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)
           volinv = 1.0d0/rec_volume(1,1,i)

           do iex=1,nof_variables
             du(iex) = zero
           end do


           kk = 0
           if (dimensiona.eq.3) then
             do degree=1,ielem_iorder(i)
               do zpow=0,degree
                 zterm = 1.0d0
                 do ip=1,zpow
                   zterm = zterm*az
                 end do
                 do ypow=0,degree-zpow
                   xpow = degree-zpow-ypow
                   kk = kk + 1
                   if (kk.le.ideg) then
                     xterm = 1.0d0
                     do ip=1,xpow
                       xterm = xterm*ax
                     end do
                     yterm = 1.0d0
                     do ip=1,ypow
                       yterm = yterm*ay
                     end do
                     raw_basis = xterm*yterm*zterm
                     phik = raw_basis - integ_basis_value(kk,i)*volinv
                     do iex=1,nof_variables
                       du(iex) = du(iex) + phik * rec_gradients(1,kk,iex,i)
                     end do

                   end if
                 end do
               end do
             end do
           else
             do degree=1,ielem_iorder(i)
               do ypow=0,degree
                 xpow = degree-ypow
                 kk = kk + 1
                 if (kk.le.ideg) then
                   xterm = 1.0d0
                   do ip=1,xpow
                     xterm = xterm*ax
                   end do
                   yterm = 1.0d0
                   do ip=1,ypow
                     yterm = yterm*ay
                   end do
                   raw_basis = xterm*yterm
                   phik = raw_basis - integ_basis_value(kk,i)*volinv
                   do iex=1,nof_variables
                     du(iex) = du(iex) + phik * rec_gradients(1,kk,iex,i)
                   end do

                 end if
               end do
             end do
           end if

           do ktail=kk+1,ideg
             phik = -integ_basis_value(ktail,i)*volinv
             do iex=1,nof_variables
               du(iex) = du(iex) + phik * rec_gradients(1,ktail,iex,i)
             end do
           end do

           do iex=1,nof_variables
             rec_uleft(iex,l,ngp,i) = u_c_val(1,iex,i) + du(iex)
           end do

         end do
       end do
  end do
#ifdef xpu
  !$omp end target teams distribute parallel do
#else
  !$omp end do
#endif
end subroutine compute_linear_reconstruction_acc2


subroutine least_squares(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i

#ifdef gpu
!$omp target teams distribute parallel do private(i)
#else
!$omp do private(i)
#endif
  do ii = 1, nof_interior
     i = el_int(ii)
     call allgrads_inner(n, i)
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

#ifdef gpu
!$omp target teams distribute parallel do private(i)
#else
!$omp do private(i)
#endif
  do ii = 1, nof_bounded
     i = el_bnd(ii)
     call allgrads_mix(n, i)
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine least_squares


! subroutine least_squares(n)
!   implicit none
!   integer, intent(in) :: n
!
!
!   !firstly all the mean
! #ifdef xpu
! !$omp barrier
! !$omp master
! call compute_gradients_mean_lsq_acc(n)
! !$omp end master
! !$omp barrier
! #else
! call compute_gradients_mean_lsq_acc(n)
! #endif
!
!
!
! if (initcond.eq.95) then
!
! call compute_gradients_inner_mean_lsq_viscous(n)  !mean viscous
!
!
!
!
! else
!
! call compute_gradients_inner_mean_lsq_viscous(n)  !mean viscous
! call compute_gradients_wall_mean_lsq_viscous      !wall viscous (only if i have some in my cpu)
!
! end if
!
!
!
!   if (nof_interior.gt.0) then
!     if (dg.eq.1) then
!       call lsq_int_mean_lsq_stage(n)
!     else
!       if ((itestcase.ge.1).and.(itestcase.le.4)) then
!         call lsq_int_mean_lsq_stage(n)
!       end if
!
!       if ((initcond.eq.95).and.((itestcase.ge.1).and.(itestcase.le.3))) then
!         call lsq_int_mean_ggs_init_stage(n)
!       end if
!
!       if (itestcase.eq.4) then
!         call lsq_int_mean_lsq_visc_stage(n)
!         call lsq_int_mean_ggs_case4_stage(n)
!
!         if (turbulence.eq.1) then
!           call lsq_int_turb_ggs_g0_stage(n)
!           call lsq_int_turb_lsq_g0_stage(n)
!           call lsq_int_turb_lsq_visc_stage(n)
!           call lsq_int_turb_lsq_g1_stage(n)
!           call lsq_int_turb_ggs_g1_stage(n)
!         end if
!       end if
!     end if
!   end if
!
!   if (nof_bounded.gt.0) then
!     if (dg.eq.1) then
!       call lsq_bnd_mean_lsq_stage(n)
!     else if (fastest.ne.1) then
!       if ((itestcase.ge.1).and.(itestcase.le.4)) then
!         call lsq_bnd_mean_lsq_stage(n)
!       end if
!
!       if ((initcond.eq.95).and.((itestcase.ge.1).and.(itestcase.le.3))) then
!         call lsq_bnd_mean_ggs_init_stage(n)
!       end if
!
!       if (itestcase.eq.4) then
!         call lsq_bnd_mean_lsq_visc_stage(n)
!         call lsq_bnd_wall_mean_lsq_stage(n)
!         call lsq_bnd_mean_ggs_case4_stage(n)
!
!         if (turbulence.eq.1) then
!           call lsq_bnd_turb_ggs_g0_stage(n)
!           call lsq_bnd_turb_lsq_g0_stage(n)
!           call lsq_bnd_turb_lsq_visc_stage(n)
!           call lsq_bnd_wall_turb_lsq_stage(n)
!           call lsq_bnd_turb_lsq_g1_stage(n)
!           call lsq_bnd_turb_ggs_g1_stage(n)
!         end if
!       end if
!     end if
!   end if
!
! end subroutine



! subroutine least_squares(n)
!   implicit none
!   integer, intent(in) :: n
!
!
!   !firstly all the mean
! #ifdef xpu
! !$omp barrier
! !$omp master
! call compute_gradients_mean_lsq_acc(n)
! !$omp end master
! !$omp barrier
! #else
! call compute_gradients_mean_lsq_acc(n)
! #endif
!
!
!
! if (initcond.eq.95) then
!
! call compute_gradients_inner_mean_lsq_viscous(n)
! call compute_gradients_inner_mean_lsq_viscous(n)
!
! else
!
! call compute_gradients_inner_mean_lsq_viscous(n)
!
!
! end if
!
!
!
!   if (nof_interior.gt.0) then
!     if (dg.eq.1) then
!       call lsq_int_mean_lsq_stage(n)
!     else
!       if ((itestcase.ge.1).and.(itestcase.le.4)) then
!         call lsq_int_mean_lsq_stage(n)
!       end if
!
!       if ((initcond.eq.95).and.((itestcase.ge.1).and.(itestcase.le.3))) then
!         call lsq_int_mean_ggs_init_stage(n)
!       end if
!
!       if (itestcase.eq.4) then
!         call lsq_int_mean_lsq_visc_stage(n)
!         call lsq_int_mean_ggs_case4_stage(n)
!
!         if (turbulence.eq.1) then
!           call lsq_int_turb_ggs_g0_stage(n)
!           call lsq_int_turb_lsq_g0_stage(n)
!           call lsq_int_turb_lsq_visc_stage(n)
!           call lsq_int_turb_lsq_g1_stage(n)
!           call lsq_int_turb_ggs_g1_stage(n)
!         end if
!       end if
!     end if
!   end if
!
!   if (nof_bounded.gt.0) then
!     if (dg.eq.1) then
!       call lsq_bnd_mean_lsq_stage(n)
!     else if (fastest.ne.1) then
!       if ((itestcase.ge.1).and.(itestcase.le.4)) then
!         call lsq_bnd_mean_lsq_stage(n)
!       end if
!
!       if ((initcond.eq.95).and.((itestcase.ge.1).and.(itestcase.le.3))) then
!         call lsq_bnd_mean_ggs_init_stage(n)
!       end if
!
!       if (itestcase.eq.4) then
!         call lsq_bnd_mean_lsq_visc_stage(n)
!         call lsq_bnd_wall_mean_lsq_stage(n)
!         call lsq_bnd_mean_ggs_case4_stage(n)
!
!         if (turbulence.eq.1) then
!           call lsq_bnd_turb_ggs_g0_stage(n)
!           call lsq_bnd_turb_lsq_g0_stage(n)
!           call lsq_bnd_turb_lsq_visc_stage(n)
!           call lsq_bnd_wall_turb_lsq_stage(n)
!           call lsq_bnd_turb_lsq_g1_stage(n)
!           call lsq_bnd_turb_ggs_g1_stage(n)
!         end if
!       end if
!     end if
!   end if
!
! end subroutine

subroutine lsq_int_mean_lsq_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
!     if ((dg.eq.1).or.(ielem_ggs(i).eq.0).or.(ielem_ggs(i).eq.1)) then
      call compute_gradients_mean_lsq(n,i,number_of_dog,number_of_nei)
!     end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_mean_lsq_stage

subroutine lsq_int_mean_ggs_init_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
!     if ((ielem_ggs(i).eq.0).or.(ielem_ggs(i).eq.1)) then
      call compute_gradients_inner_mean_ggs_viscous(n,i,number_of_dog,number_of_nei)
!     end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_mean_ggs_init_stage

subroutine lsq_int_mean_lsq_visc_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.0) then
      call compute_gradients_inner_mean_lsq_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_mean_lsq_visc_stage

subroutine lsq_int_mean_ggs_case4_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.1) then
      call compute_gradients_inner_mean_ggs_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_mean_ggs_case4_stage

subroutine lsq_int_turb_ggs_g0_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.0) then
      call compute_gradients_inner_turb_ggs_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_turb_ggs_g0_stage

subroutine lsq_int_turb_lsq_g0_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.0) then
      call compute_gradients_turb_lsq(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_turb_lsq_g0_stage

subroutine lsq_int_turb_lsq_visc_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.0) then
      call compute_gradients_turb_lsq_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_turb_lsq_visc_stage

subroutine lsq_int_turb_lsq_g1_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.1) then
      call compute_gradients_turb_lsq(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_turb_lsq_g1_stage

subroutine lsq_int_turb_ggs_g1_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_interior
    i = el_int(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.1) then
      call compute_gradients_inner_turb_ggs_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_int_turb_ggs_g1_stage

subroutine lsq_bnd_mean_lsq_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
!     if ((dg.eq.1).or.(ielem_ggs(i).eq.0).or.(ielem_ggs(i).eq.1)) then
      call compute_gradients_mean_lsq(n,i,number_of_dog,number_of_nei)
!     end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_mean_lsq_stage

subroutine lsq_bnd_mean_ggs_init_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
!     if ((ielem_ggs(i).eq.0).or.(ielem_ggs(i).eq.1)) then
      call compute_gradients_mix_mean_ggs_viscous(n,i,number_of_dog,number_of_nei)
!     end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_mean_ggs_init_stage

subroutine lsq_bnd_mean_lsq_visc_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if ((ielem_ggs(i).eq.0).and.(rec_wall(i).eq.0)) then
      call compute_gradients_inner_mean_lsq_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_mean_lsq_visc_stage

subroutine lsq_bnd_wall_mean_lsq_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if ((ielem_ggs(i).eq.0).and.(rec_wall(i).ne.0)) then
      call compute_gradients_wall_mean_lsq_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_wall_mean_lsq_stage

subroutine lsq_bnd_mean_ggs_case4_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.1) then
      call compute_gradients_mix_mean_ggs_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_mean_ggs_case4_stage

subroutine lsq_bnd_turb_ggs_g0_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.0) then
      call compute_gradients_mix_turb_ggs_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_turb_ggs_g0_stage

subroutine lsq_bnd_turb_lsq_g0_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.0) then
      call compute_gradients_turb_lsq(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_turb_lsq_g0_stage

subroutine lsq_bnd_turb_lsq_visc_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if ((ielem_ggs(i).eq.0).and.(rec_wall(i).eq.0)) then
      call compute_gradients_turb_lsq_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_turb_lsq_visc_stage

subroutine lsq_bnd_wall_turb_lsq_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if ((ielem_ggs(i).eq.0).and.(rec_wall(i).ne.0)) then
      call compute_gradients_wall_turb_lsq_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_wall_turb_lsq_stage

subroutine lsq_bnd_turb_lsq_g1_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.1) then
      call compute_gradients_turb_lsq(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_turb_lsq_g1_stage

subroutine lsq_bnd_turb_ggs_g1_stage(n)
  implicit none
  integer, intent(in) :: n
  integer :: ii, i, number_of_dog, number_of_nei
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n) private(ii,i,number_of_dog,number_of_nei)
#else
!$omp do private(ii,i,number_of_dog,number_of_nei)
#endif
  do ii = 1, nof_bounded
    i = el_bnd(ii)
    number_of_dog = ielem_idegfree(i)
    number_of_nei = ielem_inumneighbours(i)
    if (ielem_ggs(i).eq.1) then
      call compute_gradients_mix_turb_ggs_viscous(n,i,number_of_dog,number_of_nei)
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end subroutine lsq_bnd_turb_ggs_g1_stage







! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!subroutine called to employ the least sqares linear interpolation!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!for determining the slopes of each cell in each direction!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!in a weighted average way with the inverse distance !!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine piecewise_constant(n)
implicit none
!> @brief
!> subroutine for first-order scheme
integer,intent(in)::n
integer :: i, iex
integer::kmaxe
kmaxe=xmpielrank(n)




#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iex)
#else
!$omp do
#endif
	do i=1,kmaxe
    do iex=1,nof_variables
 	rec_uleft(iex,:,:,i)=u_c_val(1,iex,i)
    end do
	
	if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do iex=1,turbulenceequations+passivescalar
	rec_uleftturb(iex,:,:,i)=u_ct_val(1,iex,i)
    end do
	end if
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 
end subroutine piecewise_constant






subroutine linear_scheme(n)
 implicit none
!> @brief
!> subroutine for linear type reconstruction
 integer,intent(in)::n
integer::i
integer::kmaxe
kmaxe=xmpielrank(n)


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(n,kmaxe) private(i)
#else
!$omp do
#endif
do i=1,kmaxe
                    call compute_linear_reconstruction(n,i)
      end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



end subroutine linear_scheme







subroutine compute_linear_reconstruction(n,iconsidered)
  implicit none
!> @brief
!> subroutine for computing unlimited reconstructed solution (on-the-fly basis; stack only)
#ifdef gpu
!$omp declare target
#endif
  integer,intent(in)::n,iconsidered
	  integer :: i,l,ngp,iqp,k,ideg,iex
  real :: ax,ay,az,phik
  real :: du(1:gpu_max_nvar_total)

  i = iconsidered
  ideg = idegfree

  do l=1,ielem_ifca(i)
    if (dimensiona.eq.3) then
      if (ielem_types_faces(l,i).eq.5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if
    do ngp=1,iqp
      do iex=1,nof_variables
        rec_uleft(iex,l,ngp,i)=zero
      end do
      if (turbulenceequations.ge.1) then
        do iex=1,turbulenceequations+passivescalar
          rec_uleftturb(iex,l,ngp,i)=zero
        end do
      end if
    end do
  end do

  do l=1,ielem_ifca(i)

    if (dimensiona.eq.3) then
      if (ielem_types_faces(l,i).eq.5) then
        iqp = qp_quad
      else
        iqp = qp_triangle
      end if
    else
      iqp = qp_line
    end if

    do ngp=1,iqp
      ax = rec_qpoints(l,ngp,1,i)
      ay = rec_qpoints(l,ngp,2,i)
      if (dimensiona.eq.3) az = rec_qpoints(l,ngp,3,i)

      do iex=1,nof_variables
        du(iex) = zero
      end do
      if (turbulenceequations.ge.1) then
        do iex=1,turbulenceequations+passivescalar
          du(nof_variables+iex) = zero
        end do
      end if

      do k=1,ideg
        if (dimensiona.eq.3) then
          phik = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,0,k)
        else
          phik = basis_rec2d_value(n,ax,ay,ielem_iorder(i),i,ideg,0,k)
        end if
        do iex=1,nof_variables
          du(iex) = du(iex) + phik * rec_gradients(1,k,iex,i)
        end do
        if (turbulenceequations.ge.1) then
          do iex=1,turbulenceequations+passivescalar
            du(nof_variables+iex) = du(nof_variables+iex) + phik * rec_gradients2(1,k,iex,i)
          end do
        end if
      end do

      call extrapolate_bound_linear(n,du,l,ngp,i)
    end do
  end do

end subroutine compute_linear_reconstruction










logical function reconstruction_acc_supported()
implicit none
logical::supported_iweno

select case(iweno)
case(-1,0)
  supported_iweno=.true.
case(1)
  supported_iweno=(wenwrt.ne.2).and.(wenwrt.ne.3)
case default
  supported_iweno=.false.
end select

reconstruction_acc_supported=(dg.eq.0).and.(poly.eq.1).and.(firstorder.eq.0).and. &
& (iorder.le.gpu_max_iorder).and.supported_iweno
end function reconstruction_acc_supported


subroutine arbitrary_order(n)
implicit none
!> @brief
!> subroutine controlling the reconstruction in 3d
integer,intent(in)::n
integer::kmaxe,i
logical::use_common_acc


kmaxe=xmpielrank(n)
ielem_reduce(1:kmaxe)=0

use_common_acc=(multispecies.eq.0).and.(realgas.eq.0).and.(turbulence.eq.0).and. &
& (passivescalar.eq.0).and.reconstruction_acc_supported()

if (use_common_acc)then
#ifdef xpu
!$omp barrier
!$omp master
#endif
call compute_gradients_mean_lsq_acc(n)
call clear_rec_uleft_for_weno(n)
select case(iweno)
case(1)
  if ((wenwrt.ne.2).and.(wenwrt.ne.3))then
      if (ees.eq.5)then
        call cp_reconstruction_cweno_acc25(n)
      else
        call cp_reconstruction_weno_acc2(n)
      end if
  end if
  call checksol_acc(n)
  call muscl_acc(n)
  call checksolx_acc(n)
case(-1)
  call muscl_acc(n)
  call checksolx_acc(n)
case(0)
  call compute_linear_reconstruction_acc2(n)
  call checksolx_acc(n)
end select
if (itestcase.eq.4)then
  if (greengo.eq.0)then
    call compute_gradients_inner_mean_lsq_viscous_acc(n)
    if (allocated(rec_volume_w)) call compute_gradients_wall_mean_lsq_viscous_acc(n)
  else
    call compute_gradients_inner_mean_ggs_viscous_acc(n)
    call compute_gradients_mix_mean_ggs_viscous_acc(n)
  end if
  call solutiontriav2(n)
end if
#ifdef xpu
!$omp end master
!$omp barrier
#endif
return
end if


#ifdef xpu
!$omp barrier
!$omp master
call compute_gradients_mean_lsq_acc(n)
if ((iweno.eq.1).and.(wenwrt.ne.2).and.(wenwrt.ne.3)) then
		  if (ees.eq.5) then
		    call clear_rec_uleft_for_weno(n)
		    call cp_reconstruction_cweno_acc25(n)

		  else
		    call cp_reconstruction_weno_acc2(n)
		  end if
	  call checksolx_acc(n)
else
  call compute_linear_reconstruction_acc2(n)
  call checksolx_acc(n)
end if
if (itestcase.eq.4)then
  if (realgas.eq.0)then
    if (greengo.eq.0)then
      call compute_gradients_inner_mean_lsq_viscous_acc(n)
      if (allocated(rec_volume_w)) call compute_gradients_wall_mean_lsq_viscous_acc(n)
    else
      call compute_gradients_inner_mean_ggs_viscous_acc(n)
      call compute_gradients_mix_mean_ggs_viscous_acc(n)
    end if
    if (turbulence.eq.1)then
      call compute_gradients_turb_lsq_acc(n)
      call compute_gradients_inner_turb_ggs_viscous_acc(n)
      call compute_gradients_mix_turb_ggs_viscous_acc(n)
      if (greengo.eq.0)then
        call compute_gradients_turb_lsq_viscous_acc(n)
        if (allocated(rec_volume_w)) call compute_gradients_wall_turb_lsq_viscous_acc(n)
      end if
    end if
    call solutiontriav2(n)
  end if
end if
!$omp end master
!$omp barrier
#else
call least_squares(n)

select case(iweno)

  case(1)

  call clear_rec_uleft_for_weno(n)
  if (wenwrt.eq.2)then
  call wenoweights_char(n)
		  else
		  call wenoweights_cons(n)
		  end if
		  if (((turbulence.eq.1).or.(passivescalar.gt.0)) .and. (icoupleturb.eq.1)) then
	  call wenoweights_turb(n)
	  end if

  call checksol(n)
  call muscl(n)
  call checksolx(n)

  case(-1)

  call muscl(n)
  call checksolx(n)

  case(0)
  if (firstorder.eq.1)then
  call piecewise_constant(n)
  else
  call linear_scheme(n)
  call checksolx(n)
  end if

end select

if (itestcase.eq.4)then
  call solutiontriav2(n)
end if
#endif

end subroutine arbitrary_order

subroutine arbitrary_order_update_to_device()
implicit none
#ifdef gpu
!$omp target update to(dt,iscoun,it,t)

if (allocated(u_c_val)) then
!$omp target update to(u_c_val)
end if
if (allocated(u_ct_val)) then
!$omp target update to(u_ct_val)
end if
if (allocated(u_c_valdg)) then
!$omp target update to(u_c_valdg)
end if
if (allocated(dg2fv)) then
!$omp target update to(dg2fv)
end if
if (allocated(ielem_reduce)) then
!$omp target update to(ielem_reduce)
end if
if (allocated(ielem_troubled)) then
!$omp target update to(ielem_troubled)
end if
if (allocated(ielem_condition)) then
!$omp target update to(ielem_condition)
end if

if (allocated(solhir)) then
!$omp target update to(solhir)
end if
if (allocated(solhird)) then
!$omp target update to(solhird)
end if
if (allocated(solhir_flat)) then
!$omp target update to(solhir_flat)
end if

if (allocated(boundhir)) then
!$omp target update to(boundhir)
end if
if (allocated(boundhir_dg)) then
!$omp target update to(boundhir_dg)
end if
if (allocated(boundhiri)) then
!$omp target update to(boundhiri)
end if
if (allocated(boundhirm)) then
!$omp target update to(boundhirm)
end if
if (allocated(boundhir_flat)) then
!$omp target update to(boundhir_flat)
end if
if (allocated(boundhir_dgflat)) then
!$omp target update to(boundhir_dgflat)
end if
if (allocated(boundhiri_flat)) then
!$omp target update to(boundhiri_flat)
end if
#endif
end subroutine arbitrary_order_update_to_device

subroutine arbitrary_order_update_from_device()
implicit none
#ifdef gpu
#endif
end subroutine arbitrary_order_update_from_device



subroutine viscous_dg_ggs(n)
implicit none
!> @brief
!> subroutine controlling the reconstruction in 3d
integer,intent(in)::n
integer::i,iex,l,ngp,iqp,iconsidered
integer::kmaxe
kmaxe=xmpielrank(n)




if (dimensiona.eq.3)then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered,iex,l,iqp,ngp)
#else
!$omp do
#endif
do i=1,kmaxe
iconsidered=i


			do l=1,ielem_ifca(i)

						if (ielem_types_faces(l,i).eq.5)then
							iqp=qp_quad
						else
							iqp=qp_triangle
						end if

						do ngp=1,iqp
							rec_uleftv(1:3,1:nof_variables-1,l,ngp,i) = rec_grads(1:nof_variables-1,1:3,i)
						end do
			end do




end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



else
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(iconsidered,iex,l,iqp,ngp)
#else
!$omp do
#endif
do i=1,kmaxe
iconsidered=i



		do l=1,ielem_ifca(i)

                iqp=qp_line_n

                do ngp=1,iqp



			rec_uleftv(1:2,1:nof_variables-1,l,ngp,i) = rec_grads(1:nof_variables-1,1:2,i)

				end do
				end do



end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif












end if


end subroutine viscous_dg_ggs



subroutine checksol_acc(n)
implicit none
!> @brief
!> Ideal-gas accelerated variant of checksol for the XPU reconstruction path.
integer,intent(in)::n
integer :: i,l,ngp,iqp,iex
integer :: reduce1,kmaxe
real    :: jump_cond
real    :: mp_pinfl,mp_pinfr,gammal,gammar
real,dimension(1:gpu_max_nvar) :: leftv,rightv

kmaxe=xmpielrank(n)

if (itestcase.ge.3)then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(itestcase, nof_variables, dimensiona, qp_quad, qp_triangle, qp_line, zero, gamma, oo2) &
!$omp& firstprivate(r_gas) &
!$omp& map(alloc: xmpielrank) &
!$omp& map(alloc: u_c_val, rec_uleft, ielem_reduce, ielem_troubled, ielem_full, ielem_hybrid) &
!$omp& map(alloc: ielem_ifca, ielem_types_faces) &
!$omp& firstprivate(n,kmaxe) &
!$omp& private(i,l,ngp,iqp,iex,reduce1,jump_cond,mp_pinfl,mp_pinfr,gammal,gammar,leftv,rightv)
#else
!$omp do private(i,l,ngp,iqp,iex,reduce1,jump_cond,mp_pinfl,mp_pinfr,gammal,gammar,leftv,rightv)
#endif
  do i=1,kmaxe
    ielem_reduce(i)=0
    reduce1=0
    jump_cond=0.9d0

    if (ielem_troubled(i).eq.1)then
      if (ielem_full(i).eq.0)then
        ielem_reduce(i)=1
      end if

      if (ielem_full(i).eq.1)then
        do l=1,ielem_ifca(i)
          if (dimensiona.eq.3)then
            if (ielem_types_faces(l,i).eq.5)then
              iqp=qp_quad
            else
              iqp=qp_triangle
            end if
          else
            iqp=qp_line
          end if

          do ngp=1,iqp
            do iex=1,nof_variables
              leftv(iex)=rec_uleft(iex,l,ngp,i)
              rightv(iex)=u_c_val(1,iex,i)
            end do
            call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

            do iex=1,nof_variables
              if ((iex.ge.2).and.(iex.le.dimensiona+1)) cycle
              if ((abs(leftv(iex)-rightv(iex))).ge.(jump_cond*rightv(iex)))then
                reduce1=1
                ielem_reduce(i)=1
              end if
            end do
          end do
        end do

        if (ielem_hybrid(i).eq.1)then
          reduce1=1
          ielem_reduce(i)=1
        end if
      end if
    end if
  end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if

end subroutine checksol_acc


subroutine checksolx_acc(n)
implicit none
!> @brief
!> Ideal-gas accelerated variant of checksolx for the XPU reconstruction path.
integer,intent(in)::n
integer :: i,l,ngp,iqp,iex,k
integer :: reduce1,kmaxe
real    :: jump_cond
real    :: mp_pinfl,mp_pinfr,gammal,gammar
real,dimension(1:gpu_max_nvar) :: leftv,rightv

kmaxe=xmpielrank(n)

if (itestcase.ge.3)then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(itestcase, nof_variables, dimensiona, turbulence, turbulenceequations, passivescalar) &
!$omp& firstprivate(icoupleturb, dg, qp_quad, qp_triangle, qp_line, zero, gamma, oo2, r_gas) &
!$omp& map(alloc: xmpielrank) &
!$omp& map(alloc: u_c_val, u_ct_val, rec_uleft, rec_uleftturb, dg2fv) &
!$omp& map(alloc: ielem_reduce, ielem_troubled, ielem_hybrid, ielem_ifca, ielem_types_faces, ielem_idegfree) &
!$omp& firstprivate(n,kmaxe) &
!$omp& private(i,l,ngp,iqp,iex,k,reduce1,jump_cond,mp_pinfl,mp_pinfr,gammal,gammar,leftv,rightv)
#else
!$omp do private(i,l,ngp,iqp,iex,k,reduce1,jump_cond,mp_pinfl,mp_pinfr,gammal,gammar,leftv,rightv)
#endif
  do i=1,kmaxe
    jump_cond=0.8d0
    reduce1=0

    if (ielem_troubled(i).eq.1)then
      do l=1,ielem_ifca(i)
        if (dimensiona.eq.3)then
          if (ielem_types_faces(l,i).eq.5)then
            iqp=qp_quad
          else
            iqp=qp_triangle
          end if
        else
          iqp=qp_line
        end if

        do ngp=1,iqp
          do iex=1,nof_variables
            leftv(iex)=rec_uleft(iex,l,ngp,i)
            rightv(iex)=u_c_val(1,iex,i)
          end do
          call cons2prim2_ideal(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

          do iex=1,nof_variables
            if ((iex.ge.2).and.(iex.le.dimensiona+1)) cycle
            if ((abs(leftv(iex)-rightv(iex))).ge.(jump_cond*rightv(iex)))then
              reduce1=1
              ielem_reduce(i)=2
            end if
          end do
        end do
      end do

      if (ielem_hybrid(i).eq.1)then
        reduce1=1
        ielem_reduce(i)=1
      end if

      if (reduce1.ge.1)then
        do l=1,ielem_ifca(i)
          if (dimensiona.eq.3)then
            if (ielem_types_faces(l,i).eq.5)then
              iqp=qp_quad
            else
              iqp=qp_triangle
            end if
          else
            iqp=qp_line
          end if

          do ngp=1,iqp
            do iex=1,nof_variables
              rec_uleft(iex,l,ngp,i)=u_c_val(1,iex,i)
            end do
          end do
        end do

        if (dg.eq.1)then
          do iex=1,nof_variables
            do k=1,ielem_idegfree(i)
              dg2fv(k,iex,i)=zero
            end do
          end do
        end if
      end if

      if (turbulence.eq.1)then
        if (icoupleturb.eq.1)then
          reduce1=0
          do l=1,ielem_ifca(i)
            if (dimensiona.eq.3)then
              if (ielem_types_faces(l,i).eq.5)then
                iqp=qp_quad
              else
                iqp=qp_triangle
              end if
            else
              iqp=qp_line
            end if

            do ngp=1,iqp
              leftv(1)=rec_uleftturb(1,l,ngp,i)
              rightv(1)=u_ct_val(1,1,i)
              if (((abs(leftv(1)-rightv(1))).ge.(0.6d0*rightv(1))).or.(leftv(1).le.zero)) then
                reduce1=1
              end if
            end do
          end do

          if (ielem_hybrid(i).eq.1)then
            reduce1=1
            ielem_reduce(i)=1
          end if

          if (reduce1.eq.1)then
            do l=1,ielem_ifca(i)
              if (dimensiona.eq.3)then
                if (ielem_types_faces(l,i).eq.5)then
                  iqp=qp_quad
                else
                  iqp=qp_triangle
                end if
              else
                iqp=qp_line
              end if

              do ngp=1,iqp
                rec_uleftturb(1,l,ngp,i)=u_ct_val(1,1,i)
              end do
            end do
          end if
        end if
      end if
    end if
  end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if

end subroutine checksolx_acc
























subroutine checksol(n)
implicit none
!> @brief
!> subroutine for checking the reconstructed solution
integer,intent(in)::n
integer :: i,l,ngp,iqp,iex
integer :: reduce1,kmaxe
real    :: jump_cond,jump,rhol,rscale,templ
real    :: mp_pinfl,gammal,sumx
real    :: mp_pinfr,gammar,temp_scale
real,dimension(1:gpu_max_nvar) :: leftv,rightv
kmaxe=xmpielrank(n)











if (itestcase.ge.3)then

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(l,ngp,iqp,iex,reduce1,jump,rhol,rscale,templ, mp_pinfl,mp_pinfr,gammal,gammar,sumx,temp_scale, leftv,rightv,jump_cond)
#else
!$omp do
#endif
	do i=1,kmaxe	!all elements

        ielem_reduce(i)=0;reduce1=0
        jump_cond=0.9

        if (realgas==1) jump_cond=0.3


        if (ielem_troubled(i).eq.1)then

        if (ielem_full(i).eq.0)then
			ielem_reduce(i)=1
        end if


        if (ielem_full(i).eq.1)then
		
		  do l=1,ielem_ifca(i)	!faces2
				  if (dimensiona.eq.3)then

			      if (ielem_types_faces(l,i).eq.5)then
				    iqp=qp_quad
			      else
				    iqp=qp_triangle
			      end if
			      else

					 iqp=qp_line
			      end if

				  do ngp=1,iqp
					      
							do iex=1,nof_variables
							leftv(iex)=rec_uleft(iex,l,ngp,i)
							rightv(iex)=u_c_val(1,iex,i)
							end do
							call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)




                                                    if (realgas.eq.0)then
                                                    do iex=1,nof_variables
                                                           if ((iex.ge.2).and.(iex.le.dimensiona+1)) cycle



                                                            if (((abs(leftv(iex)-rightv(iex))).ge.(jump_cond*rightv(iex))))then
																	reduce1=1
																ielem_reduce(i)=1
															end if

                                                    end do
                                                    end if
                                                      if (realgas.eq.1)then
                                                              do iex=1,nof_variables       !loop rho,u,v,w,e,p
                                                                    if ((iex.ge.2).and.(iex.le.dimensiona+1)) cycle


                                                                      jump  = abs(leftv(iex) - rightv(iex))
                                                                      temp_scale = rightv(iex)

                                                                      if (abs(temp_scale).gt.10e-300)then
                                                                        if ((jump .ge. jump_cond*(temp_scale)).or.(rightv(iex).lt.0.0d0))then
                                                                                reduce1=1
                                                                                ielem_reduce(i)=1    !jump from not species
                                                                                exit

                                                                        end if
                                                                      end if

                                                              end do
                                                    end if











				
					
				  end do
		end do	
		
		
		if (ielem_hybrid(i).eq.1)then
		reduce1=1
		ielem_reduce(i)=1
		end if
		


		
		
		end if
		
		
		end if




	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
		

		
end if		
		
		

end subroutine checksol

 
subroutine checksolx(n)
implicit none
!> @brief
!> subroutine for checking the reconstructed solution
integer,intent(in)::n
integer::i,l,ngp,iqp,iex,k
integer::reduce1,kmaxe
real::jump_cond,jump,rhol,rscale
real,dimension(1:gpu_max_nvar)::leftv
real::mp_pinfl,gammal,sumx
real,dimension(1:gpu_max_nvar)::rightv
real::mp_pinfr,gammar,temp_scale

kmaxe=xmpielrank(n)






if (itestcase.ge.3)then

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(l,ngp,iqp,iex,k,reduce1,jump,rhol,rscale,mp_pinfl,mp_pinfr,gammal,gammar,sumx,temp_scale, leftv,rightv,jump_cond)
#else
!$omp do
#endif
	do i=1,kmaxe
            jump_cond=0.8

            if (Realgas.eq.1)jump_cond=0.2

			reduce1=0

        if (ielem_troubled(i).eq.1)then
        

						do l=1,ielem_ifca(i)	!faces2
								if (dimensiona.eq.3)then

									if (ielem_types_faces(l,i).eq.5)then
										iqp=qp_quad
									else
										iqp=qp_triangle
									end if
									else

										iqp=qp_line
									end if
										do ngp=1,iqp


													do iex=1,nof_variables
													leftv(iex)=rec_uleft(iex,l,ngp,i)
													rightv(iex)=u_c_val(1,iex,i)
													end do

													call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

                                                            if (realgas.eq.0)then
                                                    do iex=1,nof_variables
                                                           if ((iex.ge.2).and.(iex.le.dimensiona+1)) cycle

                                                            if (((abs(leftv(iex)-rightv(iex))).ge.(jump_cond*rightv(iex))))then
																	reduce1=1
																ielem_reduce(i)=2
															end if

                                                    end do
                                                    end if
                                                     if (realgas.eq.1)then
                                                              do iex=1,nof_variables       !loop rho,u,v,w,e,p
                                                                    if ((iex.ge.2).and.(iex.le.dimensiona+1)) cycle


                                                                      jump  = abs(leftv(iex) - rightv(iex))
                                                                      temp_scale = rightv(iex)

                                                                      if (abs(temp_scale).gt.10e-300)then
                                                                        if ((jump .ge. jump_cond*(temp_scale)).or.(rightv(iex).lt.0.0d0))then
                                                                                reduce1=1
                                                                                ielem_reduce(i)=iex    !jump from not species
                                                                                exit

                                                                        end if
                                                                      end if

                                                              end do
                                                    end if


										end do
						end do


					if (ielem_hybrid(i).eq.1)then
					reduce1=1
					ielem_reduce(i)=1
					end if





						if (reduce1.ge.1)then
							do l=1,ielem_ifca(i)
								if (dimensiona.eq.3)then
									if (ielem_types_faces(l,i).eq.5)then
										iqp=qp_quad
									else
										iqp=qp_triangle
									end if
								else
									iqp=qp_line
								end if
								do ngp=1,iqp
									do iex=1,nof_variables
										rec_uleft(iex,l,ngp,i)=u_c_val(1,iex,i)
									end do
								end do
							end do





							if (dg.eq.1)then
							do iex=1,nof_variables
								do k=1,ielem_idegfree(i)
									dg2fv(k,iex,i)=zero
								end do
							end do

							end if

					end if


					if (turbulence.eq.1)then
							if (icoupleturb.eq.1)then
							reduce1=0
								do l=1,ielem_ifca(i)	!faces2

										if (dimensiona.eq.3)then

										if (ielem_types_faces(l,i).eq.5)then
											iqp=qp_quad
										else
											iqp=qp_triangle
										end if
										else

											iqp=qp_line
										end if

										do ngp=1,iqp
											leftv(1)=rec_uleftturb(1,l,ngp,i)
											rightv(1)=u_ct_val(1,1,i)
												if (((abs(leftv(1)-rightv(1))).ge.(0.6*rightv(1))).or.(leftv(1).le.zero)) then
													reduce1=1
												end if
										end do
								end do

								if (ielem_hybrid(i).eq.1)then
								reduce1=1
								ielem_reduce(i)=1
								end if

									if (reduce1.eq.1)then
									do l=1,ielem_ifca(i)
										if (dimensiona.eq.3)then
											if (ielem_types_faces(l,i).eq.5)then
												iqp=qp_quad
											else
												iqp=qp_triangle
											end if
										else
											iqp=qp_line
										end if
										do ngp=1,iqp
											rec_uleftturb(1,l,ngp,i)=u_ct_val(1,1,i)
										end do
									end do
									end if
							end if
					end if





		end if



	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end if



end subroutine checksolx
 
 
subroutine slope_limiters(n,iconsidered,utmin,utmax,usol,psi)
implicit none
!> @brief
!> Pointwise MUSCL slope limiter evaluated at a single quadrature point.
!> Inputs:
!>   usol(:)     : unlimited reconstructed values at the current quadrature point
!>   utmin/utmax : min/max bounds for element iconsidered (per variable)
!> Output:
!>   psi(:)      : limiter factors in [0,1] (or as defined by selected limiter), per variable
#ifdef gpu
!$omp declare target
#endif
integer, intent(in) :: n, iconsidered
real,    intent(in) :: utmin(1:nof_variables+turbulenceequations+passivescalar), utmax(1:nof_variables+turbulenceequations+passivescalar)
real,    intent(in) :: usol(1:nof_variables+turbulenceequations+passivescalar)
real,    intent(out):: psi(1:nof_variables+turbulenceequations+passivescalar)
integer :: i, iex, nvtot
real :: u0, d2, sfd, epsi2, dmin, dplus, kappa_ven, psi2, sig_1, delu, y_fun, s_y, pol_mog

i = iconsidered
nvtot = nof_variables+turbulenceequations+passivescalar

kappa_ven = 0.1

do iex = 1, nvtot

  ! cell-centered value for this variable (no utemp storage)
  if (iex <= nof_variables) then
    u0 = u_c_val(1,iex,i)
  else
    u0 = u_ct_val(1,iex-nof_variables,i)
  end if

  psi2 = zero
  d2 = usol(iex) - u0

  if (abs(d2) <= zero) then
    psi(iex) = 1.0d0

  else if (d2 > zero) then
    sfd = (utmax(iex) - u0) / d2

    select case (limiter)

    case (1)
      psi(iex) = min(1.0d0, sfd)                 ! Barth-Jespersen

    case (10)
      psi(iex) = min(0.3d0, sfd)                 ! very restrictive BJ

    case (2)
      pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
      if (sfd < 1.5d0) then
        psi(iex) = pol_mog
      else
        psi(iex) = 1.0d0
      end if

      ! BJ Michalak blending
      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = (2.0d0*y_fun**3) - (3.0d0*y_fun**2) + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (9)
      pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
      if (sfd < 1.5d0) then
        psi(iex) = pol_mog
      else
        psi(iex) = 1.0d0
      end if

      ! same blending as case(2)
      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = (2.0d0*y_fun**3) - (3.0d0*y_fun**2) + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (3)
      pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
      if (sfd < 1.5d0) then
        psi(iex) = pol_mog
      else
        psi(iex) = 1.0d0
      end if

      ! extended-stencil blending
      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (4)    ! VKM + blending
      dmin = usol(iex) - u0
      dmin = sign(1.0, dmin) * (abs(dmin) + tolsmall)
      dplus = utmax(iex) - u0
      epsi2 = (kappa_ven*ielem_minedge(i))**3
      psi(iex) = (1.0d0/dmin) * ( (((dplus**2)+epsi2)*dmin + (2.0d0*(dmin**2)*dplus)) / ((dplus**2) + (2.0d0*dmin**2) + (dmin*dplus) + epsi2) )

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (5)
      psi(iex) = (sfd**2 + sfd) / (sfd**2 + 1.0d0)   ! Van Albada

    case (6)
      psi(iex) = 2.0d0*sfd / (sfd + 1.0d0)           ! Van Leer

    case (7)    ! Venkatakrishnan
      dmin  = usol(iex) - u0
      dmin  = sign(1.0, dmin) * (abs(dmin) + tolsmall)
      dplus = utmax(iex) - u0
      epsi2 = (kappa_ven*ielem_minedge(i))**3
      psi(iex) = (1.0d0/dmin) * ( (((dplus**2)+epsi2)*dmin + (2.0d0*(dmin**2)*dplus)) / ((dplus**2) + (2.0d0*dmin**2) + (dmin*dplus) + epsi2) )

    case (8)    ! Venkatakrishnan + blending
      dmin  = usol(iex) - u0
      dmin  = sign(1.0, dmin) * (abs(dmin) + tolsmall)
      dplus = utmax(iex) - u0
      epsi2 = (kappa_ven*ielem_minedge(i))**3
      psi(iex) = (1.0d0/dmin) * ( (((dplus**2)+epsi2)*dmin + (2.0d0*(dmin**2)*dplus)) / ((dplus**2) + (2.0d0*dmin**2) + (dmin*dplus) + epsi2) )

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    end select

  else
    sfd = (utmin(iex) - u0) / d2

    select case (limiter)

    case (1)
      psi(iex) = min(1.0d0, sfd)

    case (10)
      psi(iex) = min(0.3d0, sfd)

    case (2)
      pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
      if (sfd < 1.5d0) then
        psi(iex) = pol_mog
      else
        psi(iex) = 1.0d0
      end if

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (9)
      pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
      if (sfd < 1.5d0) then
        psi(iex) = pol_mog
      else
        psi(iex) = 1.0d0
      end if

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (3)
      pol_mog = -((4.0d0/27.0d0)*sfd**3) + sfd
      if (sfd < 1.5d0) then
        psi(iex) = pol_mog
      else
        psi(iex) = 1.0d0
      end if

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (4)
      dmin  = usol(iex) - u0
      dmin  = sign(1.0, dmin) * (abs(dmin) + tolsmall)
      dplus = utmin(iex) - u0
      epsi2 = (kappa_ven*ielem_minedge(i))**3
      psi(iex) = (1.0d0/dmin) * ( (((dplus**2)+epsi2)*dmin + (2.0d0*(dmin**2)*dplus)) / ((dplus**2) + (2.0d0*dmin**2) + (dmin*dplus) + epsi2) )

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    case (5)
      psi(iex) = (sfd**2 + sfd) / (sfd**2 + 1.0d0)

    case (6)
      psi(iex) = 2.0d0*sfd / (sfd + 1.0d0)

    case (7)
      dmin  = usol(iex) - u0
      dmin  = sign(1.0, dmin) * (abs(dmin) + tolsmall)
      dplus = utmin(iex) - u0
      epsi2 = (kappa_ven*ielem_minedge(i))**3
      psi(iex) = (1.0d0/dmin) * ( (((dplus**2)+epsi2)*dmin + (2.0d0*(dmin**2)*dplus)) / ((dplus**2) + (2.0d0*dmin**2) + (dmin*dplus) + epsi2) )

    case (8)
      dmin  = usol(iex) - u0
      dmin  = sign(1.0, dmin) * (abs(dmin) + tolsmall)
      dplus = utmin(iex) - u0
      epsi2 = (kappa_ven*ielem_minedge(i))**3
      psi(iex) = (1.0d0/dmin) * ( (((dplus**2)+epsi2)*dmin + (2.0d0*(dmin**2)*dplus)) / ((dplus**2) + (2.0d0*dmin**2) + (dmin*dplus) + epsi2) )

      delu = utmax(iex) - utmin(iex)
      if (delu**2 <= ((kappa_ven*ielem_minedge(i))**3)) then
        sig_1 = 1.0d0
      else if ( ((kappa_ven*ielem_minedge(i))**3) < delu**2 .and. delu**2 < 2.0d0*((kappa_ven*ielem_minedge(i))**3) ) then
        y_fun = (delu**2 - (kappa_ven*ielem_minedge(i))**3) / ((kappa_ven*ielem_minedge(i))**3)
        s_y   = 2.0d0*y_fun**3 - 3.0d0*y_fun**2 + 1.0d0
        sig_1 = s_y
      else
        sig_1 = 0.0d0
      end if

      psi2 = sig_1 + (1.0d0 - sig_1)*psi(iex)
      psi(iex) = psi2

    end select

  end if

end do

end subroutine slope_limiters





subroutine process_state(uin,uref,k,utmin,utmax,sumvars,aver_vars,maxvars)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  integer::nvtot
  real,    intent(in)    :: uin(nof_variables+turbulenceequations+passivescalar),uref(nof_variables+turbulenceequations+passivescalar)
  integer, intent(inout) :: k
  real,    intent(inout) :: utmin(nof_variables+turbulenceequations+passivescalar), utmax(nof_variables+turbulenceequations+passivescalar)
  real,    intent(inout) :: sumvars(nof_variables+turbulenceequations+passivescalar), aver_vars(nof_variables+turbulenceequations+passivescalar), maxvars(nof_variables+turbulenceequations+passivescalar)
  ! locals
  integer :: j
  real :: mp_pinfl, gammal
  real :: uwork(gpu_max_nvar_total)
  real :: leftv(1:gpu_max_nvar)
  nvtot=nof_variables+turbulenceequations+passivescalar
  uwork(1:nvtot) = uin(1:nvtot)

  ! Convert conservative->primitive for flow variables if requested
  if (wenwrt == 3) then
    leftv(1:nof_variables) = uwork(1:nof_variables)
    call cons2prim(n, leftv, mp_pinfl, gammal)
    uwork(1:nof_variables) = leftv(1:nof_variables)
  end if

  k = k + 1


    do j = 1, nvtot
      sumvars(j) = sumvars(j) + abs(uwork(j) - uref(j))
    end do


  do j = 1, nvtot
    utmin(j)     = min(utmin(j), uwork(j))
    utmax(j)     = max(utmax(j), uwork(j))
    aver_vars(j) = aver_vars(j) + uwork(j)
    maxvars(j)   = max(maxvars(j), abs(uwork(j)))
  end do

end subroutine process_state




subroutine trouble_indicator1
implicit none
integer::i,l,j,k,kmaxe,iqp,ngp,iex
integer::trouble
integer::iconsidered,facex,pointx
real,dimension(1:gpu_max_nvar)::leftv,rightv
real,dimension(1:gpu_max_nvar)::maxvars,aver_vars,sumvars,utmin,utmax
real :: usol(1:gpu_max_nvar)




if (code_profile.ne.102)then

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i, l, j, k, iqp, ngp, iex, trouble, iconsidered, facex, pointx, leftv, rightv, maxvars, aver_vars, sumvars, utmin, utmax, usol)
#else
!$omp do
#endif
do i = 1, xmpielrank(n)
iconsidered=i





    call find_bounds(iconsidered,maxvars,aver_vars,sumvars,utmin,utmax)

    do l = 1, ielem_ifca(i)

            if (dimensiona.eq.2)then

            iqp=qp_line_n
            else
                if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad
				  else
					iqp=qp_triangle
				  end if
            end if

                do ngp = 1,iqp!
                    facex=l
                    pointx=ngp
                    usol(1:nof_variables)=rec_uleft_dg(1:nof_variables,facex,pointx,iconsidered)
                    leftv(1:nof_variables)=usol(1:nof_variables)
                    call pad_dg(iconsidered,leftv)
                    call nad_dg(iconsidered,facex,pointx,leftv,rightv,usol,maxvars,aver_vars,sumvars,utmin,utmax)
                end do

    end do







end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if




end subroutine



subroutine trouble_indicator2
implicit none
integer::i,l,j,k,kmaxe,iqp,ngp,iex,ndof
integer::trouble, ifree,i_deg
integer::iconsidered,facex,pointx

kmaxe=xmpielrank(n)

if (code_profile.ne.102)then


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i, l, j, k, iqp, ngp, iex, ndof, trouble, ifree, i_deg, iconsidered, facex, pointx)
#else
!$omp do
#endif
do i = 1, xmpielrank(n)
iconsidered=i



    if (ielem_troubled(i).eq.1)then

    do l = 1, ielem_ifca(i)

            if (dimensiona.eq.2)then

            iqp=qp_line_n
            else
                if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad
				  else
					iqp=qp_triangle
				  end if
            end if

                do ngp = 1,iqp!
                    facex=l
                    pointx=ngp
                    rec_uleft_dg(:, facex,pointx,iconsidered)=rec_uleft(:, facex,pointx,iconsidered)
                end do

    end do

      do iex=1,nof_variables

!
      u_c_valdg(1,iex,2:idegfree+1,iconsidered)=dg2fv(1:idegfree,iex,iconsidered)

      end do




 end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if



end subroutine












    
    
subroutine pad_dg(iconsidered,leftv)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer::i,l,j,k,kmaxe,iqp,ngp,iex
integer::trouble
integer,intent(in)::iconsidered
real,dimension(1:nof_variables),intent(inout)::leftv
real::mp_pinfl,gammal

    i=iconsidered
                                                if (itestcase.ge.3)then
														if (dimensiona.eq.3)then

                                                    call cons2prim(n,leftv,mp_pinfl,gammal)


                                                    if(multispecies.eq.1)then

															if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then
																ielem_troubled(i)=1;ielem_condition(i)=1
															end if
															if ((leftv(5).le.zero).or.(leftv(5).ne.leftv(5)))then
																ielem_troubled(i)=1;ielem_condition(i)=1
															end if
	! 														if ((leftv(5).le.-mp_pinf(1)).or.(leftv(5).le.-mp_pinf(2)))then
	! 															ielem_troubled(i)=1;ielem_condition(i)=1
	!                                                         end if


														!	if((leftv(nof_variables).lt.-0.05).or.leftv(nof_variables).ne.leftv(8)) then
														!		ielem_troubled(i) =1; ielem_condition(i)=1
														!	end if
													else

						!
															if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then
															ielem_troubled(i)=1;ielem_condition(i)=1


															end if
															if ((leftv(5).le.zero).or.(leftv(5).ne.leftv(5)))then
															ielem_troubled(i)=1;ielem_condition(i)=1
															end if
													end if
                                                else
												call cons2prim(n,leftv,mp_pinfl,gammal)
                                                if (multispecies.eq.1)then

													if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then
														ielem_troubled(i)=1;ielem_condition(i)=1
													end if
													if ((leftv(4).le.zero))then
														ielem_troubled(i)=1;ielem_condition(i)=1
													end if
													if ((leftv(4).ne.leftv(4)))then
														ielem_troubled(i)=1;ielem_condition(i)=1
													end if

! if ((leftv(4).le.zero).or.(leftv(4).ne.leftv(4)))then
! 														ielem_troubled(i)=1;ielem_condition(i)=1
! 													end if

													if((leftv(nof_variables).lt.zero).or.leftv(nof_variables).gt.1.0d0) then
														ielem_troubled(i) =1; ielem_condition(i)=1
													end if





                                                else

                                                        if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then
                                                        ielem_troubled(i)=1;ielem_condition(i)=1
                                                        end if
                                                        if ((leftv(4).le.zero).or.(leftv(4).ne.leftv(4)))then
                                                        ielem_troubled(i)=1;ielem_condition(i)=1
                                                        end if

                                                end if

                                                end if
                                                end if



                                                if (ielem_hybrid(i).eq.1)then
                                                   ielem_troubled(i)=1
                                                end if





end subroutine


subroutine nad_dg(iconsidered,facex,pointx,leftv,rightv,usol,maxvars,aver_vars,sumvars,utmin,utmax)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer::i,l,j,k,kmaxe,iqp,ngp,iex
integer::trouble,img
integer,intent(in)::iconsidered,facex,pointx
real::par1,par2,d2,minb,maxb
real,dimension(1:gpu_max_nvar)::nad_dg_el
real,dimension(1:nof_variables),intent(inout)::leftv,rightv
real,dimension(1:nof_variables),intent(in)::maxvars,aver_vars,sumvars,utmin,utmax
real,dimension(1:nof_variables),intent(in) :: usol
real::mp_pinfl,gammal




! par1=1e-4
! par2=4e-1

        if (dimensiona.eq.3)then

        img=5
    else

        img=4
    end if



       select case(indicator_type)


             case(1)       !mood indicator

            do iex=1,nof_variables
			nad_dg_el(iex)=max(indicator_par1,(indicator_par2)*(utmax(iex)-utmin(iex)))
			end do
			leftv(1:nof_variables)=usol(1:nof_variables)


			if (dimensiona.eq.2)then
			call cons2prim(n,leftv,mp_pinfl,gammal)
			else
			call cons2prim(n,leftv,mp_pinfl,gammal)
			end if

            do iex=1,nof_variables
                if ((leftv(iex).lt.(utmin(iex)-nad_dg_el(iex))).or.(leftv(iex).gt.(utmax(iex)+nad_dg_el(iex))))then
                    ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
                end if
            end do



             case(11)       !mood indicator

            do iex=1,nof_variables
			nad_dg_el(iex)=max(indicator_par1,(indicator_par2)*(utmax(iex)-utmin(iex)))
			end do


            do iex=1,nof_variables
                if ((usol(iex).lt.(utmin(iex)-nad_dg_el(iex))).or.(usol(iex).gt.(utmax(iex)+nad_dg_el(iex))))then
                    ielem_condition(iconsidered)=1
                     ielem_troubled(iconsidered)=1;

                end if
            end do


            case(2)         !shu indicator

                        do iex=1,nof_variables


                         if ((sumvars(iex)/maxvars(iex)).gt.indicator_par1)then
                            ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
                        end if

                        end do

             case(22)         !shock detector, indicator (only density and energy)

                        do iex=1,nof_variables

                        if ((iex.eq.1).or.(iex.eq.img))then
                         if ((sumvars(iex)/maxvars(iex)).gt.indicator_par1)then
                            ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
                        end if
                        end if

                        end do

             case(3)         !dmp

                        do iex=1,nof_variables


                         if ((usol(iex).gt.(utmax(iex))).or.(usol(iex).lt.(utmin(iex))))then
                            ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
                        end if

                        end do



             case(4)    !minmod
                        do iex=1,nof_variables



                         if ((iex.eq.1).or.(iex.eq.img))then
                        if (abs(usol(iex)-u_c_val(1,iex,iconsidered)).gt.(indicator_par1*u_c_val(1,iex,iconsidered)))then
                        ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
                        end if
                        end if
                        end do


             case(5)    !all troubled


                        ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1


              case(6)       !mood indicator only density & energy

            do iex=1,nof_variables
            if ((iex.eq.1).or.(iex.eq.img))then
			nad_dg_el(iex)=max(indicator_par1,(indicator_par2)*(utmax(iex)-utmin(iex)))
			end if
			end do
			leftv(1:nof_variables)=usol(1:nof_variables)

! 			if (dimensiona.eq.2)then
! 			call cons2prim(n,leftv,mp_pinfl,gammal)
! 			else
! 			call cons2prim(n,leftv,mp_pinfl,gammal)
! 			end if

            do iex=1,nof_variables
            if ((iex.eq.1).or.(iex.eq.img))then

                if ((leftv(iex).lt.(utmin(iex)-nad_dg_el(iex))).or.(leftv(iex).gt.(utmax(iex)+nad_dg_el(iex))))then
                    ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
                end if
                end if
            end do


			case(9)       !mood indicator only density & energy

			 if (ielem_filtered(iconsidered).eq.1)then
				ielem_troubled(iconsidered)=1;ielem_condition(iconsidered)=1
			end if








            end select


end subroutine nad_dg









subroutine apply_filter(n)
    implicit none
    integer,intent(in) :: n
    integer :: i, j, k, kmaxe

!     kmaxe = xmpielrank(n)

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,k)
#else
!$omp do
#endif
    do i = 1, xmpielrank(n)
      if (ielem_filtered(i) == 1) then
        do j = 1, nof_variables
          do k = 1, idegfree
            rhs_valdg(k+1, j, i) = rhs_valdg(k+1, j, i) * modal_filter_weak(k)
          end do
        end do
      end if
    end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  end subroutine apply_filter


subroutine apply_filter_dg(n)
    implicit none
    integer,intent(in) :: n
    integer :: i, kmaxe
    real    :: ex1, ex2, energy_ratio
    real    :: u2, u3, u4, ws2, ws3, ws4, ww2, ww3, ww4



#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(ex1,ex2,energy_ratio,u2,u3,u4,ws2,ws3,ws4,ww2,ww3,ww4)
#else
!$omp do
#endif
    do i = 1, xmpielrank(n)
      ielem_filtered(i) = 0

      u2  = u_c_val (1,2,i);  u3  = u_c_val (1,3,i);  u4  = u_c_val (1,4,i)
      ws2 = u_cs_val(1,2,i);  ws3 = u_cs_val(1,3,i);  ws4 = u_cs_val(1,4,i)
      ww2 = u_cw_val(1,2,i);  ww3 = u_cw_val(1,3,i);  ww4 = u_cw_val(1,4,i)

      ex2 = (u2-ww2)*(u2-ww2) + (u3-ww3)*(u3-ww3) + (u4-ww4)*(u4-ww4)
      ex1 = (u2-ws2)*(u2-ws2) + (u3-ws3)*(u3-ws3) + (u4-ws4)*(u4-ws4)

      energy_ratio = (ex2 + 1.0d-31) / (ex1 + 1.0d-31)

      ielem_er1dt(i) = (ielem_er1(i) - ex1) / dt
      ielem_er2dt(i) = (ielem_er2(i) - ex2) / dt

      ielem_er (i) = energy_ratio
      ielem_er1(i) = ex1
      ielem_er2(i) = ex2

      if (ielem_er(i) > 1.15d0) then
        ielem_filtered(i) = 1
      end if
    end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  end subroutine apply_filter_dg



  subroutine apply_filter1(n)
    implicit none
    integer,intent(in) :: n
    integer :: i, j, k, kmaxe



#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,k)
#else
!$omp do
#endif
    do i = 1, xmpielrank(n)
      do j = 1, nof_variables
        do k = 1, idegfree
          u_c_valdg(1, j, k+1, i) = u_c_valdg(1, j, k+1, i) * modal_filter(k)
        end do
      end do
    end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  end subroutine apply_filter1


  subroutine apply_filter2(n)
    implicit none
    integer,intent(in) :: n
    integer :: i, j, k, kmaxe

    ! Equivalent to your original apply_filter2 but with NO temp array.
!     kmaxe = xmpielrank(n)

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,k)
#else
!$omp do
#endif
    do i = 1, xmpielrank(n)
      do j = 1, nof_variables
        do k = 1, idegfree
          rhs_valdg(k+1, j, i) = rhs_valdg(k+1, j, i) * modal_filter(k)
        end do
      end do
    end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  end subroutine apply_filter2





subroutine filter_init(n)
    implicit none
    integer,intent(in) :: n
    integer :: fil_i, j
    integer :: last, curr
    real    :: xorder, rfil_alpha, rfil_nc, rfil_s, rfil_i
    real    :: filx, coeff

    ! Defensive init
    modal_filter(1:idegfree)        = 1.0d0
    modal_filter_weak(1:idegfree)   = 1.0d0
    modal_filter_strong(1:idegfree) = 1.0d0

    last = 0
    xorder = real(iorder, kind=8)

      ! Weak/strong step filters:
      ! Weak
      last = 0
      do fil_i = 1, iorder
        curr = (((fil_i+1)*(fil_i+2)*(fil_i+3))/6) - 1
        if (fil_i < iorder) then
          coeff = 1.0d0
        else
          coeff = 0.0d0
        end if
        do j = last+1, curr
          if (j >= 1 .and. j <= idegfree) modal_filter_weak(j) = coeff
        end do
        last = curr
      end do

      ! Strong
      last = 0
      do fil_i = 1, iorder
        curr = (((fil_i+1)*(fil_i+2)*(fil_i+3))/6) - 1
        if (fil_i < iorder-1) then
          coeff = 1.0d0
        else
          coeff = 0.0d0
        end if
        do j = last+1, curr
          if (j >= 1 .and. j <= idegfree) modal_filter_strong(j) = coeff
        end do
        last = curr
      end do



  end subroutine filter_init







subroutine adda_filter(n,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer :: kbasis
  integer,intent(in) :: n, iconsidered
  integer :: i, j, k, ngp, ideg, icompwrt
  real    :: ax, ay, az, wgt
  real    :: ex1, ex2
  real    :: s_unf, s_str, s_weak
  real    :: val_unf, val_str, val_weak
  real    :: mp_pinfl, mp_pinfr, gammal, gammar
  real    :: energy_ratio
  real    :: phi(1:gpu_max_dof)

  i    = iconsidered
  ideg = ielem_idegfree(i)

  if (dg.eq.1) return

  ex1 = 0.0d0
  ex2 = 0.0d0

  icompwrt = 0

  if (adda_type.eq.1) then
    ! single evaluation at element center
    ax = 0.0d0; ay = 0.0d0; az = 0.0d0

      do kbasis=1,ideg
        phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,icompwrt,kbasis)
      end do

    do k = 2, 4
      s_unf  = 0.0d0
      s_str  = 0.0d0
      s_weak = 0.0d0

      do j = 1, ideg
        s_unf  = s_unf  + phi(j) * rec_gradients(1,j,k,i)
        s_str  = s_str  + phi(j) * rec_gradients(1,j,k,i) * adda_filter_strong(j)
        s_weak = s_weak + phi(j) * rec_gradients(1,j,k,i) * adda_filter_weak(j)
      end do

      val_unf  = u_c_val(1,k,i) + s_unf
      val_str  = u_c_val(1,k,i) + s_str
      val_weak = u_c_val(1,k,i) + s_weak

      ex1 = ex1 + (val_unf - val_str )*(val_unf - val_str )
      ex2 = ex2 + (val_unf - val_weak)*(val_unf - val_weak)
    end do

  else if (adda_type.eq.2) then
    ! quadrature average
    do ngp = 1, ielem_itotalpoints(i)
      ax  = qp_array_x(ngp,i)
      ay  = qp_array_y(ngp,i)
      az  = qp_array_z(ngp,i)
      wgt = qp_array_qp_weight(ngp,i)

        do kbasis=1,ideg
          phi(kbasis) = basis_rec_value(n,ax,ay,az,ielem_iorder(i),i,ideg,icompwrt,kbasis)
        end do

      do k = 2, 4
        s_unf  = 0.0d0
        s_str  = 0.0d0
        s_weak = 0.0d0

        do j = 1, ideg
          s_unf  = s_unf  + phi(j) * rec_gradients(1,j,k,i)
          s_str  = s_str  + phi(j) * rec_gradients(1,j,k,i) * adda_filter_strong(j)
          s_weak = s_weak + phi(j) * rec_gradients(1,j,k,i) * adda_filter_weak(j)
        end do

        val_unf  = u_c_val(1,k,i) + s_unf
        val_str  = u_c_val(1,k,i) + s_str
        val_weak = u_c_val(1,k,i) + s_weak

        ex1 = ex1 + (val_unf - val_str )*(val_unf - val_str ) * wgt
        ex2 = ex2 + (val_unf - val_weak)*(val_unf - val_weak) * wgt
      end do
    end do

  else
    ex1 = 0.0d0
    ex2 = 0.0d0
  end if

  energy_ratio = (ex2 + 1.0d-31) / (ex1 + 1.0d-31)

  ielem_er1dt(i) = (ielem_er1(i) - ex1) / dt
  ielem_er2dt(i) = (ielem_er2(i) - ex2) / dt

  ielem_er (i) = energy_ratio
  ielem_er1(i) = ex1
  ielem_er2(i) = ex2

  ielem_er1er2(i) = abs(ielem_er1dt(i)) / ielem_er2dt(i)

  call apply_adda_filter(n,i)

end subroutine adda_filter





subroutine apply_adda_filter(n,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
  integer,intent(in) :: n, iconsidered
  integer :: i
  real    :: lwcx1

  i = iconsidered
  lwcx1 = lwci1

  if (rungekutta.eq.11) then
    if (iscoun.eq.1) then
      lwcx1 = lwci1
      if (ielem_er(i).gt.1.2)  lwcx1 = 10.0d0
      if (ielem_er(i).le.0.95) lwcx1 = 1000.0d0
      ielem_lwcx2(i) = lwcx1
    else
      lwcx1 = ielem_lwcx2(i)
    end if
  else
    if (ielem_er(i).gt.1.2)  lwcx1 = 10.0d0
    if (ielem_er(i).le.0.95) lwcx1 = 1000.0d0
    ielem_lwcx2(i) = lwcx1
  end if

  if (ielem_full(i).eq.0) then
    ielem_lwcx2(i) = -10.0d0
  end if

  ielem_linc(i) = lwcx1

end subroutine apply_adda_filter




subroutine fix_dissipation(n)
implicit none
real::check1
real::check
integer::i,j,k,l
integer,intent(in)::n
integer::kmaxe
kmaxe=xmpielrank(n)


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe

	if (ielem_full(i).eq.1)then
	!1)reduce dissipation
	if (ielem_lwcx2(i).gt.10)then
		if (ielem_wcx(i).ge.0.999)then
		ielem_diss(i)=max(ielem_diss(i)-0.1,0.5d0)	!reduce dissipation even more
		else
		ielem_diss(i)=1.0d0							!increase dissipation if shock
		end if
	else
	!2) increase dissipation
		if (ielem_wcx(i).ge.0.999)then
		ielem_diss(i)=min(ielem_diss(i)+0.1,1.0d0)	!increase dissipation even more
		else
		ielem_diss(i)=1.0d0							!increase dissipation if shock
		end if
	end if
	end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end subroutine fix_dissipation




subroutine fix_dissipation2(n)
implicit none
real::check1
real::check
integer::i,j,k,l,iconsidered
integer,intent(in)::n
integer::kmaxe
kmaxe=xmpielrank(n)



#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i, iconsidered)
#else
!$omp do
#endif
do i=1,kmaxe
	iconsidered=i

	call find_bounds_diss(iconsidered)

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end subroutine fix_dissipation2






subroutine find_bounds_diss(iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer::i,l,j,k,kmaxe,iqp,ngp,iex,ik,iq,nf,lf,rowf
integer,intent(in)::iconsidered
real,dimension(1:6,1)::utemp

i=iconsidered


utemp(:,:)=1.0d0



		ielem_facediss(:,i)=1.0d0


            if (ielem_interior(i).eq.0)then
                do l = 1, ielem_ifca(i)

                    utemp(l,1)=ielem_diss(ielem_ineigh(l,i))
                end do
            end if

            if (ielem_interior(i).eq.1)then
			    do l=1,ielem_ifca(i)
                    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                            if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu

                                utemp(l,1)=ielem_diss(ielem_ineigh(l,i))
                                else
                                !not periodic ones in my cpu
                                end if
                            else

                                utemp(l,1)=ielem_diss(ielem_ineigh(l,i))
                            end if
                    else	!in other cpus they can only be periodic or mpi neighbours

                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                            if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu

!                             utemp(l,1)=iexsolhird(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
!                             (rec_ihexl(1,ielem_indexi(l,i),i),1)

                            nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								utemp(l,1)=solhird(rowf)


                            end if
                        else




!                         utemp(l,1)=iexsolhird(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
!                         (rec_ihexl(1,ielem_indexi(l,i),i),1)


                        nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
								lf=rec_ihexl(1,ielem_indexi(L,i),i)
								rowf=halo_offset(nf) + lf - 1
								utemp(l,1)=solhird(rowf)




                        end if

                    end if

			  end do
         end if





            if (ielem_interior(i).eq.0)then
                do l = 1, ielem_ifca(i)
                    ielem_facediss(l,i)=max(utemp(l,1),ielem_diss(i))

                end do
            end if

            if (ielem_interior(i).eq.1)then
			    do l=1,ielem_ifca(i)
                    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                            if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                ielem_facediss(l,i)=max(utemp(l,1),ielem_diss(i))
!                                 utemp(l,1)=ielem_diss(ielem_ineigh(l,i))
                                else
                                !not periodic ones in my cpu
                                end if
                            else
                                ielem_facediss(l,i)=max(utemp(l,1),ielem_diss(i))

                            end if
                    else	!in other cpus they can only be periodic or mpi neighbours

                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                            if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu

                            ielem_facediss(l,i)=max(utemp(l,1),ielem_diss(i))
                            end if
                        else


                        ielem_facediss(l,i)=max(utemp(l,1),ielem_diss(i))
                        end if

                    end if

			  end do
         end if








end subroutine find_bounds_diss



subroutine vfbp_limiter
implicit none
real::vf_qpsol, vf_avsol, lthresh, hthresh,scaling,scaling1,scaling2, pd1_qpsol, pd1_avsol,pd2_qpsol, pd2_avsol
integer::i,ii,icd,k,number_of_nei,idummy,l,nnd,ngp, iqp,i_elem,i_face, img1,img2
real, dimension(1:gpu_max_qp_face)::scaling_f,scaling_f1,scaling_f2
integer::facex,pointx,iconsidered,number_of_dog


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i, vf_qpsol, vf_avsol, lthresh, hthresh, scaling, scaling1, scaling2, pd1_qpsol, pd1_avsol, pd2_qpsol, pd2_avsol, ii, icd, k, number_of_nei, idummy, l, nnd, ngp, iqp, i_elem, i_face, img1, img2, facex, pointx, iconsidered, number_of_dog, scaling_f, scaling_f1, scaling_f2)
#else
!$omp do
#endif
do i = 1, xmpielrank(n)
    lthresh = 1.0e-16
    hthresh = 1.0d0-lthresh


    do l = 1, ielem_ifca(i)



    if (dimensiona.eq.2)then

            iqp=qp_line_n
            img1=5
            img2=6
            else
            img1=6
            img2=7
                if (ielem_types_faces(l,i).eq.5)then
					iqp=qp_quad
				  else
					iqp=qp_triangle

				  end if
            end if

 !compute scaling factors at each face and each qp

     do ngp = 1,iqp! qp_line_n

                facex=l
                pointx=ngp
                iconsidered=i
                number_of_dog=ielem_idegfree(i)

                vf_qpsol = rec_uleft_dg(nof_variables, facex, pointx,iconsidered)
                vf_avsol = u_c_valdg(1,nof_variables,1,i)

                pd1_qpsol = rec_uleft_dg(img1, facex, pointx,iconsidered)
                pd1_avsol = u_c_valdg(1,img1,1,i)

                pd2_qpsol = rec_uleft_dg(img2, facex, pointx,iconsidered)
                pd2_avsol = u_c_valdg(1,img2,1,i)



        !volume fraction
                if (vf_qpsol.lt.lthresh) then

                    scaling_f(ngp) = (vf_avsol - lthresh)/(vf_avsol - vf_qpsol)

                else if ((vf_qpsol.gt.lthresh).and.(vf_qpsol.lt.hthresh)) then

                    scaling_f(ngp) = 1.0d0

                else if (vf_qpsol.gt.hthresh) then

                    scaling_f(ngp) = (hthresh - vf_avsol)/(vf_qpsol - vf_avsol)

                end if




        !partial densities
                if (pd1_qpsol.lt.lthresh) then

                    scaling_f1(ngp) = (pd1_avsol - lthresh)/(pd1_avsol - pd1_qpsol)
                else

                    scaling_f1(ngp) = 1.0d0

                end if



                if (pd2_qpsol.lt.lthresh) then

                    scaling_f2(ngp) = (pd2_avsol - lthresh)/(pd2_avsol - pd2_qpsol)

                else

                    scaling_f2(ngp) = 1.0d0

                end if



            end do



            scaling = minval(scaling_f(:))
            scaling1 = minval(scaling_f1(:))
            scaling2 = minval(scaling_f2(:))





            do ngp = 1,iqp! qp_line_n

                facex=l
                pointx=ngp
                iconsidered=i
                number_of_dog=ielem_idegfree(i)

                vf_qpsol = rec_uleft_dg(nof_variables, facex, pointx,iconsidered)
                vf_avsol = u_c_valdg(1,nof_variables,1,i)

                pd1_qpsol = rec_uleft_dg(img1, facex, pointx,iconsidered)
                pd1_avsol = u_c_valdg(1,img1,1,i)

                pd2_qpsol = rec_uleft_dg(img2, facex, pointx,iconsidered)
                pd2_avsol = u_c_valdg(1,img2,1,i)


            !volume fraction scaling at qps
                if((scaling.gt.0.0d0).and.(scaling.lt.1.0d0)) then
                    rec_uleft_dg(nof_variables, facex, pointx,iconsidered) = vf_avsol + scaling*(vf_qpsol - vf_avsol)

                end if

            !partial densities scaling at qps

                if((scaling1.gt.0.0d0).and.(scaling1.lt.1.0d0)) then
                    rec_uleft_dg(img1, facex, pointx,iconsidered) = pd1_avsol + scaling1*(pd1_qpsol - pd1_avsol)

                end if

                if((scaling2.gt.0.0d0).and.(scaling2.lt.1.0d0)) then
                    rec_uleft_dg(img2, facex, pointx,iconsidered) = pd2_avsol + scaling2*(pd2_qpsol - pd2_avsol)

                end if




            end do





    end do

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end subroutine vfbp_limiter



subroutine adda_filter_init(n)
  implicit none
  integer,intent(in)::n
  integer :: fil_i, j
  integer :: last, curr
  real    :: f2, f3

  ! Strong/weak filters depend only on iorder and idegfree

  adda_filter_strong(1:idegfree) = 0.0d0
  adda_filter_weak  (1:idegfree) = 0.0d0

  last = 0
  do fil_i = 1, iorder
    curr = (((fil_i+1)*(fil_i+2)*(fil_i+3))/6) - 1

    ! ---- strong filter coefficient for this "order band"
    if (iorder == 2) then
      f2 = 0.0d0
    else
      if (fil_i < 2) then
        f2 = 1.0d0
      else
        f2 = 0.0d0
      end if
    end if

    ! ---- weak filter coefficient for this "order band"
    if (fil_i <= 2) then
      f3 = 1.0d0
    else
      f3 = 0.0d0
    end if

    do j = last+1, curr
      if (j >= 1 .and. j <= idegfree) then
        adda_filter_strong(j) = f2
        adda_filter_weak  (j) = f3
      end if
    end do

    last = curr
  end do

end subroutine adda_filter_init







! ! ---------------------------------------------------------------------------------------------!

! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module recon
