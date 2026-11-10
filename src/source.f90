module source
   use mpi
   use library
   use transform
   use local
   use riemann
   use flow_operations
   use declaration
   implicit none

   integer, parameter :: realgas_source_max_subcycles = 64
   integer, parameter :: realgas_source_newton_maxiter = 12
   integer, parameter :: realgas_source_newton_linesearch = 10

   ! Steady SA-RC state: local cells and existing reconstruction halos.
   real, allocatable :: rc_xyz(:,:), rc_strain(:,:), rc_inverse(:,:,:), rc_scale(:,:)
   real, allocatable :: rc_send(:), rc_recv(:), sa_rc_fr1(:)
   integer, allocatable :: rc_requests(:)
   integer :: rc_cells = 0
#ifdef xpu
!$omp declare target(rc_xyz, rc_strain, rc_inverse, rc_scale, rc_send, rc_recv, sa_rc_fr1, rc_cells)
#endif
   private :: rc_xyz, rc_strain, rc_inverse, rc_scale, rc_send, rc_recv, rc_requests
   private :: rc_cells, rc_fail, rc_neighbour, rc_displacement, rc_exchange

contains
   subroutine axisymmetric_source_state(n,q,qt,g,gt,radius,wall,src,st)
      ! S in dQ/dt + dFx/dx + dFr/dr = S. Planar FV areas/fluxes are unchanged.
      ! Cell-centre quadrature is second order on smooth, regular solutions.
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: n
      real,intent(in) :: q(gpu_max_nvar),qt(gpu_max_extra_transport)
      real,intent(in) :: g(gpu_max_nvar,2),gt(gpu_max_extra_transport,2),radius,wall
      real,intent(out) :: src(gpu_max_nvar),st(gpu_max_extra_transport)
      real :: prim(gpu_max_nvar),mu,ktr,kve,cvv,gm,pinf,rho,u,v,hoop,mut,sk,sw
      real :: packed(20),txx,txr,trr,ttt,coef,rmix,grada,gsum,jsum,ysum
      real :: dif(gpu_max_species),htr(gpu_max_species),hv(gpu_max_species)
      real :: yy(gpu_max_species),gy(gpu_max_species),jflux(gpu_max_species)
      integer :: k,nt
      src=0.0; st=0.0
      prim=q; call cons2prim(n,prim(1:nof_variables),pinf,gm)
      rho=q(1); u=q(2)/rho; v=q(3)/rho; hoop=v/radius
      src(1:nof_variables)=-hoop*q(1:nof_variables)
      src(4)=src(4)-hoop*prim(4)
      nt=turbulenceequations+passivescalar
      st(1:nt)=-hoop*qt(1:nt)
      if(itestcase/=4) return

      kve=0.0; cvv=0.0
      if(realgas==1) then
         call multispecies_mixtures_rg(q(1:nof_variables),mu,ktr,kve,dif,htr,hv,gm,cvv)
      else
         call molecular_transport(n,q(1:nof_variables),mu,ktr)
      end if
      mut=0.0; sk=1.0; sw=1.0
      if(turbulence==1) then
         packed=0.0; packed(1)=wall; packed(2:3)=qt(1:2)
         packed(4:5)=g(1,:); packed(6:7)=g(2,:)
         packed(8:9)=gt(1,:); packed(10:11)=gt(2,:); packed(12)=hoop
         call eddy_viscosity_state(2,rho,mu,packed,.true.,mut,sk,sw)
         ktr=ktr+mut*translational_cp(q(1:nof_variables))/prtu
         kve=kve+mut*cvv/prtu
      end if
      call axis_stress(mu+mut,g(1,1),g(1,2),g(2,1),g(2,2),hoop,txx,txr,trr,ttt)
      src(2)=src(2)+txr/radius
      src(3)=src(3)+(trr-ttt)/radius
      src(4)=src(4)+(u*txr+v*trr)/radius
      if(realgas==1) then
         ! Match the existing RG_RELAX transport selection in the face flux.
         if(rg_relax>=1) then
            yy(1:nof_species)=max(q(6:5+nof_species)/rho,0.0)
            ysum=sum(yy(1:nof_species))
            if(ysum>1.0e-30) then
               yy(1:nof_species)=yy(1:nof_species)/ysum
            else
               yy(1:nof_species)=1.0/real(nof_species)
            end if
            gy(1:nof_species)=g(5:4+nof_species,2)
            gsum=sum(gy(1:nof_species))
            gy(1:nof_species)=gy(1:nof_species)-yy(1:nof_species)*gsum
            rmix=sum(yy(1:nof_species)/rg_molm(1:nof_species))
            grada=sum(gy(1:nof_species)/rg_molm(1:nof_species))
            do k=1,nof_species
               jflux(k)=-rho*dif(k)*(gy(k)-yy(k)*grada/rmix)-mut/schmidt_turb*gy(k)
            end do
            jsum=sum(jflux(1:nof_species))
            jflux(1:nof_species)=jflux(1:nof_species)-yy(1:nof_species)*jsum
            src(6:5+nof_species)=src(6:5+nof_species)-jflux(1:nof_species)/radius
            src(4)=src(4)+(ktr*g(3,2)+kve*g(4,2)- &
               sum((htr(1:nof_species)+hv(1:nof_species))*jflux(1:nof_species)))/radius
            src(5)=src(5)+(kve*g(4,2)-sum(hv(1:nof_species)*jflux(1:nof_species)))/radius
         end if
      else
         src(4)=src(4)+ktr*g(3,2)/radius
      end if
      do k=1,nt
         coef=mu+mut
         if(turbulence==1.and.k<=turbulenceequations) then
            if(turbulencemodel==1) coef=(mu+max(qt(1),0.0))/sigma
            if(turbulencemodel==2) then
               if(k==1) coef=mu+mut/sk
               if(k==2) coef=mu+mut/sw
            end if
         end if
         st(k)=st(k)+coef*gt(k,2)/radius
      end do
   end subroutine axisymmetric_source_state

   subroutine axisymmetric_cell_data(i,q,qt,g,gt,wall)
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: i
      real,intent(out) :: q(gpu_max_nvar),qt(gpu_max_extra_transport)
      real,intent(out) :: g(gpu_max_nvar,2),gt(gpu_max_extra_transport,2),wall
      integer :: nt
      q=0.0; qt=0.0; g=0.0; gt=0.0; wall=0.0
      q(1:nof_variables)=u_c_val(i,1,1:nof_variables)
      nt=turbulenceequations+passivescalar
      if(nt>0) qt(1:nt)=u_ct_val(i,1,1:nt)
      if(itestcase==4) then
         g(1:nof_variables-1,:)=rec_grads(i,1:nof_variables-1,1:2)
         if(nt>0) gt(1:nt,:)=rec_grads(i,nof_variables:nof_variables+nt-1,1:2)
         if(turbulence==1) wall=ielem_walldist(i)
      end if
   end subroutine axisymmetric_cell_data

   subroutine sources_axisymmetric_cell(n,i)
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: n,i
      real :: q(gpu_max_nvar),qt(gpu_max_extra_transport),g(gpu_max_nvar,2),gt(gpu_max_extra_transport,2)
      real :: src(gpu_max_nvar),st(gpu_max_extra_transport),wall
      integer :: nt
      call axisymmetric_cell_data(i,q,qt,g,gt,wall)
      call axisymmetric_source_state(n,q,qt,g,gt,axis_radius(i),wall,src,st)
      rhs_val(i,1:nof_variables)=rhs_val(i,1:nof_variables)-ielem_totvolume(i)*src(1:nof_variables)
      nt=turbulenceequations+passivescalar
      if(nt>0) rhst_val(i,1:nt)=rhst_val(i,1:nt)-ielem_totvolume(i)*st(1:nt)
   end subroutine sources_axisymmetric_cell

   logical function axisymmetric_source_admissible(q) result(ok)
      ! A Jacobian probe must be in the domain of the transport laws. In
      ! particular, positive density alone does not make Sutherland sqrt(T)
      ! valid in a cold, high-Mach cell. Do not clip the state or temperature.
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real,intent(in) :: q(gpu_max_nvar)
      real :: rho,ev,etr,rmix,cv,y(gpu_max_species),vel2,temp
      ok=.false.
      if(.not.all(ieee_is_finite(q(1:nof_variables)))) return
      if(q(1)<=tiny(1.0)) return
      if(realgas==1) then
         if(any(q(6:5+nof_species)<0.0).or.q(5)<0.0) return
         call realgas_thermo(q(1:nof_variables),rho,y,ev,etr,rmix,cv)
         if(.not.all(ieee_is_finite([rho,ev,etr,rmix,cv]))) return
         if(etr<=0.0.or.cv<=0.0.or.rmix<=0.0) return
      else
         rho=q(1)
         vel2=(q(2)/rho)**2+(q(3)/rho)**2
         ! Use the same arithmetic as molecular_transport, not a pressure floor.
         temp=(gamma-1.0)*(q(4)/rho-0.5*vel2)/r_gas
         if(.not.ieee_is_finite(temp)) return
         if(temp<=0.0) return
      end if
      ok=.true.
   end function axisymmetric_source_admissible

   subroutine axisymmetric_source_jacobian(n,i,jac,tdiag)
      ! Local block derivative with reconstructed gradients frozen, as in the existing
      ! approximate implicit operator. Chemistry is assembled separately, once.
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: n,i
      real,intent(out) :: jac(gpu_max_nvar,gpu_max_nvar),tdiag(gpu_max_extra_transport)
      real :: q(gpu_max_nvar),qp(gpu_max_nvar),qt(gpu_max_extra_transport),tp(gpu_max_extra_transport)
      real :: g(gpu_max_nvar,2),gt(gpu_max_extra_transport,2),base(gpu_max_nvar),sp(gpu_max_nvar)
      real :: bt(gpu_max_extra_transport),st(gpu_max_extra_transport),wall,h,scale(gpu_max_nvar),rho,energy,delta
      integer :: k,nt,attempt,direction
      logical :: evaluated
      jac=0.0; tdiag=0.0
      if(axisymmetric_enabled==0) return
      call axisymmetric_cell_data(i,q,qt,g,gt,wall)
      if(.not.axisymmetric_source_admissible(q)) error stop 'Invalid base state in axisymmetric source Jacobian'
      call axisymmetric_source_state(n,q,qt,g,gt,axis_radius(i),wall,base,bt)
      rho=q(1); energy=max(abs(q(4)),rho)
      scale=rho; scale(2:3)=sqrt(rho*energy); scale(4)=energy
      if(realgas==1) scale(5)=max(abs(q(5)),rho)
      scale(1:nof_variables)=max(scale(1:nof_variables),abs(q(1:nof_variables)))
      do k=1,nof_variables
         h=sqrt(epsilon(1.0))*scale(k)
         evaluated=.false.
         do attempt=1,32
            ! Keep the existing forward difference when admissible; otherwise
            ! use the backward side before reducing the probe. This matters
            ! when increasing momentum would consume the last internal energy.
            do direction=1,-1,-2
               qp=q; qp(k)=q(k)+real(direction)*h
               delta=qp(k)-q(k)
               if(delta==0.0) cycle
               if(.not.axisymmetric_source_admissible(qp)) cycle
               call axisymmetric_source_state(n,qp,qt,g,gt,axis_radius(i),wall,sp,st)
               jac(1:nof_variables,k)=(sp(1:nof_variables)-base(1:nof_variables))/delta
               evaluated=.true.
               exit
            end do
            if(evaluated) exit
            h=0.5*h
         end do
         if(.not.evaluated) error stop 'No admissible axisymmetric Jacobian perturbation'
      end do
      nt=turbulenceequations+passivescalar
      do k=1,nt
         h=sqrt(epsilon(1.0))*max(abs(qt(k)),1.0e-8*rho)
         tp=qt; tp(k)=tp(k)+h
         call axisymmetric_source_state(n,q,tp,g,gt,axis_radius(i),wall,sp,st)
         tdiag(k)=(st(k)-bt(k))/h
      end do
   end subroutine axisymmetric_source_jacobian

   subroutine add_axisymmetric_implicit(n,i,block)
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: n,i
      real,intent(inout) :: block(:,:)
      real :: jac(gpu_max_nvar,gpu_max_nvar),td(gpu_max_extra_transport)
      if(axisymmetric_enabled==0) return
      call axisymmetric_source_jacobian(n,i,jac,td)
      block(1:nof_variables,1:nof_variables)=block(1:nof_variables,1:nof_variables)- &
         ielem_totvolume(i)*jac(1:nof_variables,1:nof_variables)
   end subroutine add_axisymmetric_implicit

   subroutine add_axisymmetric_transport_implicit(n,i,diag)
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: n,i
      real,intent(inout) :: diag(:)
      real :: jac(gpu_max_nvar,gpu_max_nvar),td(gpu_max_extra_transport),floor_value
      integer :: k
      if(axisymmetric_enabled==0) return
      call axisymmetric_source_jacobian(n,i,jac,td)
      floor_value=max(turb_diag_abs_floor,turb_diag_floor_frac*ielem_totvolume(i)/ielem_dtl(i))
      do k=1,turbulenceequations+passivescalar
         diag(k)=max(floor_value,diag(k)-ielem_totvolume(i)*td(k))
      end do
   end subroutine add_axisymmetric_transport_implicit

   real function axisymmetric_source_radius(n,i) result(rate)
#ifdef xpu
!$omp declare target
#endif
      integer,intent(in) :: n,i
      real :: jac(gpu_max_nvar,gpu_max_nvar),td(gpu_max_extra_transport)
      real :: scale(gpu_max_nvar),rho,energy
      integer :: j,k
      rate=0.0
      if(axisymmetric_enabled==0) return
      call axisymmetric_source_jacobian(n,i,jac,td)
      ! Scale the block before taking a row-sum bound: its raw entries have different units.
      rho=u_c_val(i,1,1); energy=max(abs(u_c_val(i,1,4)),rho)
      scale=rho; scale(2:3)=sqrt(rho*energy); scale(4)=energy
      if(realgas==1) scale(5)=max(abs(u_c_val(i,1,5)),rho)
      do j=1,nof_variables
         rate=max(rate,sum(abs(jac(j,1:nof_variables))*scale(1:nof_variables))/scale(j))
      end do
      if(turbulenceequations+passivescalar>0) rate=max(rate,maxval(abs(td(1:turbulenceequations+passivescalar))))
   end function axisymmetric_source_radius

   real function realgas_arrhenius(a, b, theta, t, logkc) result(rate)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: a, b, theta, t, logkc
      real :: lograte
! a is already in SI molar units. Subtract ln(Kc) BEFORE exponentiation:
! the forward coefficient can underflow at a cold wall while kb stays finite.
      lograte = log(a) + b * log(t) - theta / t - logkc
      if (lograte .lt. log(tiny(1.0e0))) then
         rate = 0.0e0
      else
         rate = exp(min(lograte, log(huge(1.0e0)) - 4.0e0))
      end if
   end function realgas_arrhenius

   subroutine realgas_air5_rates(t, tv, c, wdot)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: t, tv, c(5)
      real, intent(out) :: wdot(5)
      real :: logkc(5), a(5, 3), b(3), theta(3), kf(3), kb(3), ta, af4, bf4, af5, bf5, theta4, theta5, tf
      integer :: r, k
      wdot = 0.0e0
      if ((rg_nof_reactions .ne. 5) .and. (rg_nof_reactions .ne. 6)) return
      call realgas_log_equilibrium(t, logkc)
      ta = exp(0.6e0 * log(t) + 0.4e0 * log(tv))
      theta = [113200.0e0, 59500.0e0, 75500.0e0]
      if (rg_nof_reactions .eq. 6) then
         a(:, 1) = [7.0e21, 7.0e21, 7.0e21, 3.0e22, 3.0e22]
         a(:, 2) = [2.0e21, 2.0e21, 2.0e21, 1.0e22, 1.0e22]
         a(:, 3) = [5.0e15, 5.0e15, 1.1e17, 1.1e17, 1.1e17]
         b = [-1.6e0, -1.5e0, 0.0e0]
         ! Park N2+O -> NO+N: 6.4e17 cm^3/(mol s), NOT 6.4e7.
         af4 = 6.4e17; bf4 = -1.0e0; theta4 = 38370.0e0
         af5 = 8.4e12; bf5 = 0.0e0; theta5 = 19450.0e0
         tf = t
      else
         ! Retain the legacy option-5 forward coefficients; use the same SI
         ! equilibrium convention and translational reverse temperature as option 6.
         a(:, 1) = [3.78e18, 3.78e18, 3.78e18, 1.11e18, 1.11e18]
         a(:, 2) = [2.75e16, 2.75e16, 2.75e16, 8.25e16, 8.25e16]
         a(:, 3) = [2.30e14, 2.30e14, 2.30e14, 4.60e14, 4.60e14]
         b = [-1.6e0, -1.0e0, -0.5e0]
         af4 = 3.18e10; bf4 = 0.1e0; theta4 = 37700.0e0
         af5 = 2.16e5; bf5 = 1.29e0; theta5 = 19220.0e0
         tf = ta
      end if
      a = a * 1.0e-6
      kf = 0.0e0; kb = 0.0e0
      do r = 1, 3
         do k = 1, 5
            kf(r) = kf(r) + c(k) * realgas_arrhenius(a(k, r), b(r), theta(r), ta, 0.0e0)
            kb(r) = kb(r) + c(k) * realgas_arrhenius(a(k, r), b(r), theta(r), t, logkc(r))
         end do
      end do
      wdot(1) = kf(1) * c(1) - kb(1) * c(4) * c(4)
      wdot(2) = kf(2) * c(2) - kb(2) * c(5) * c(5)
      wdot(3) = kf(3) * c(3) - kb(3) * c(4) * c(5)
      wdot(4) = realgas_arrhenius(af4 * 1.0e-6, bf4, theta4, tf, 0.0e0) * c(1) * c(5) &
         & - realgas_arrhenius(af4 * 1.0e-6, bf4, theta4, t, logkc(4)) * c(3) * c(4)
      wdot(5) = realgas_arrhenius(af5 * 1.0e-6, bf5, theta5, tf, 0.0e0) * c(3) * c(5) &
         & - realgas_arrhenius(af5 * 1.0e-6, bf5, theta5, t, logkc(5)) * c(2) * c(4)
   end subroutine realgas_air5_rates

   subroutine realgas_source_jacobian(n, qstate, jac)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: qstate(nof_variables)
      real, intent(out) :: jac(gpu_max_nvar, gpu_max_nvar)
      real :: qp(gpu_max_nvar), qm(gpu_max_nvar), s0(gpu_max_nvar), sp(gpu_max_nvar), sm(gpu_max_nvar)
      real :: diag(gpu_max_nvar), scale(gpu_max_nvar), h, rho, energy, momentum
      integer :: j, attempt
      logical :: plus_ok, minus_ok
      jac = 0.0e0
      if ((realgas .ne. 1) .or. (rg_relax .ne. 2)) return
      rho = max(qstate(1), 1.0e-30)
      energy = max(abs(qstate(dimensiona + 2)), rho)
      momentum = sqrt(rho * energy)
      scale = rho
      scale(2:dimensiona + 1) = momentum
      scale(dimensiona + 2:dimensiona + 3) = energy
      call sources_realgas_state(n, qstate, s0, diag)
      do j = 1, nof_variables
         h = 1.0e-6 * max(abs(qstate(j)), scale(j))
         do attempt = 1, 20
            qp = 0.0e0; qm = 0.0e0
            qp(1:nof_variables) = qstate; qm(1:nof_variables) = qstate
            qp(j) = qstate(j) + h; qm(j) = qstate(j) - h
            plus_ok = realgas_source_state_admissible(qp)
            minus_ok = realgas_source_state_admissible(qm)
            if (plus_ok .or. minus_ok) exit
            h = 0.5e0 * h
         end do
         if (plus_ok) call sources_realgas_state(n, qp, sp, diag)
         if (minus_ok) call sources_realgas_state(n, qm, sm, diag)
         if (plus_ok .and. minus_ok) then
            jac(1:nof_variables, j) = (sp(1:nof_variables) - sm(1:nof_variables)) / (2.0e0 * h)
         elseif (plus_ok) then
            jac(1:nof_variables, j) = (sp(1:nof_variables) - s0(1:nof_variables)) / h
         elseif (minus_ok) then
            jac(1:nof_variables, j) = (s0(1:nof_variables) - sm(1:nof_variables)) / h
         end if
      end do
   end subroutine realgas_source_jacobian

   real function realgas_source_radius(n, q) result(radius)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, intent(in) :: q(nof_variables)
      real :: jac(gpu_max_nvar, gpu_max_nvar), scale(gpu_max_nvar), rho, energy
      integer :: i
      call realgas_source_jacobian(n, q, jac)
      rho = max(q(1), 1.0e-30); energy = max(abs(q(dimensiona + 2)), rho)
      scale = rho; scale(2:dimensiona + 1) = sqrt(rho * energy)
      scale(dimensiona + 2:dimensiona + 3) = energy
      radius = 0.0e0
      do i = 1, nof_variables
         radius = max(radius, sum(abs(jac(i, 1:nof_variables)) * scale(1:nof_variables)) / scale(i))
      end do
   end function realgas_source_radius
   subroutine sources_computation(n)
      implicit none
!> @brief
!> sources computation
      integer, intent(in) :: n
      integer :: i, kmaxe

      call sa_rc_prepare(n)
      kmaxe = xmpielrank(n)
      if ((turbulence .eq. 1) .and. (kmaxe .gt. 0)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe, n) &
!$omp& private(i) &
!$omp& map(alloc: ielem_totvolume, ielem_walldist, ielem_walltrans, mp_janaf, mp_m, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rec_grads, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb, rhst_val, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_computation_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2) .and. (rungekutta .ge. 10) .and. (kmaxe .gt. 0)) then
#ifdef xpu
         !$omp target teams distribute parallel do &
!$omp& firstprivate(nof_variables, kmaxe) &
!$omp& map(alloc: sht_rg) &
!$omp& private(i)
#else
         !$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            ! Clear the implicit source Jacobian before assembling this step.
            sht_rg(i, 1:nof_variables) = 0.00
         end do
#ifdef xpu
         !$omp end target teams distribute parallel do
#else
         !$omp end do
#endif
      end if
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2) .and. (kmaxe .gt. 0)) then
#ifdef xpu
         !$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe, n) &
!$omp& private(i) &
!$omp& map(alloc: gamma_in, ielem_totvolume, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in, mp_tmid_in, rg_hzero) &
!$omp& map(alloc: rg_molm, rg_thetag, rhs_val, sht_rg, u_c_val)
#else
         !$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_computation_realgas_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
   end subroutine sources_computation

   subroutine sources_computation_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      integer :: iv
      real, dimension(gpu_max_turbulence) :: source_t

      if (turbulence .eq. 1) then
         call sources(n, iconsidered, source_t)
         do iv = 1, turbulenceequations
            rhst_val(iconsidered, iv) = rhst_val(iconsidered, iv) - source_t(iv) * ielem_totvolume(iconsidered)
         end do
      end if

   end subroutine sources_computation_cell

   subroutine sources_computation_realgas_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      integer :: iv
      real :: source_r(gpu_max_nvar), src_jac_diag(gpu_max_nvar)
! Pseudo-time changes the update, never the physical source residual.
      call sources_realgas(n, iconsidered, source_r, src_jac_diag, 0.0e0)

      do iv = 1, nof_variables
         rhs_val(iconsidered, iv) = rhs_val(iconsidered, iv) - source_r(iv) * ielem_totvolume(iconsidered)
         if (rungekutta .ge. 10) sht_rg(iconsidered, iv) = src_jac_diag(iv) * ielem_totvolume(iconsidered)
      end do

   end subroutine sources_computation_realgas_cell

   subroutine realgas_integrate_source_subcycles(n, iconsidered, dt_loc, source_r, src_jac_diag, source_ok)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, intent(in) :: dt_loc
      real, dimension(1:gpu_max_nvar), intent(out) :: source_r, src_jac_diag
      logical, intent(out) :: source_ok
      integer :: iv, isub, nsub
      real :: dt_sub, onsub
      logical :: attempt_ok, step_ok
      real, dimension(1:gpu_max_nvar) :: qold, qwork, qtrial, jac_sub, jac_accum

      source_ok = .false.
      source_r(1:nof_variables) = zero
      src_jac_diag(1:nof_variables) = zero

      if (dt_loc .le. 1.0e-300) return

      qold(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)

      nsub = 1
      do while (nsub .le. realgas_source_max_subcycles)
         onsub = 1.00 / real(nsub)
         dt_sub = dt_loc * onsub
         qwork(1:nof_variables) = qold(1:nof_variables)
         jac_accum(1:nof_variables) = zero
         attempt_ok = .true.

         do isub = 1, nsub
            call realgas_coupled_source_newton(n, iconsidered, qwork, dt_sub, qtrial, jac_sub, step_ok)

            if (.not. step_ok) then
               attempt_ok = .false.
               exit
            end if

            qwork(1:nof_variables) = qtrial(1:nof_variables)
            jac_accum(1:nof_variables) = jac_accum(1:nof_variables) + jac_sub(1:nof_variables)
         end do

         if (attempt_ok) then
            do iv = 1, nof_variables
               source_r(iv) = (qwork(iv) - qold(iv)) / dt_loc
               src_jac_diag(iv) = jac_accum(iv) * onsub
            end do
            source_ok = .true.
            exit
         end if

         nsub = 2 * nsub
      end do

      u_c_val(iconsidered, 1, 1:nof_variables) = qold(1:nof_variables)

   end subroutine realgas_integrate_source_subcycles

   subroutine realgas_coupled_source_newton(n, iconsidered, qold, dt_sub, qnew, jac_diag_eff, step_ok)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, intent(in) :: dt_sub
      real, dimension(1:gpu_max_nvar), intent(in) :: qold
      real, dimension(1:gpu_max_nvar), intent(out) :: qnew, jac_diag_eff
      logical, intent(out) :: step_ok
      integer :: iter, ls, ia, ib, idxa, idxb, nloc
      real :: resnorm, testnorm, eps, alpha, scale
      logical :: linear_ok, accepted
      real, dimension(1:gpu_max_nvar) :: qcurr, qtest, qpert
      real, dimension(1:gpu_max_nvar) :: res, res_test, rhs, delta
      real, dimension(1:gpu_max_nvar) :: source_eval, jac_eval, source_pert, jac_pert
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar) :: amat

      step_ok = .false.
      qnew(1:nof_variables) = qold(1:nof_variables)
      jac_diag_eff(1:nof_variables) = zero

      if ((realgas .ne. 1) .or. (nof_species .le. 0)) then
         step_ok = .true.
         return
      end if

      nloc = nof_species + 1
      if ((dimensiona + 3 + nof_species) .gt. nof_variables) return
      if (nloc .gt. gpu_max_nvar) return

      qcurr(1:nof_variables) = qold(1:nof_variables)
      if (.not. realgas_source_state_admissible(qcurr)) return

      do iter = 1, realgas_source_newton_maxiter
         call realgas_source_residual(n, iconsidered, qold, qcurr, dt_sub, res, source_eval, jac_eval, resnorm)
         if ((resnorm .ne. resnorm) .or. (abs(resnorm) .gt. 1.0e300)) return

         if (resnorm .le. 1.0e-9) then
            qnew(1:nof_variables) = qcurr(1:nof_variables)
            jac_diag_eff(1:nof_variables) = jac_eval(1:nof_variables)
            step_ok = .true.
            return
         end if

         do ib = 1, nloc
            if (ib .eq. 1) then
               idxb = dimensiona + 3
            else
               idxb = dimensiona + 2 + ib
            end if

            eps = 1.0e-6 * max(abs(qcurr(idxb)), 1.0e-12)
            qpert(1:nof_variables) = qcurr(1:nof_variables)
            qpert(idxb) = qpert(idxb) + eps

            if (.not. realgas_source_state_admissible(qpert)) then
               eps = 0.50 * eps
               qpert(1:nof_variables) = qcurr(1:nof_variables)
               qpert(idxb) = qpert(idxb) + eps
               if (.not. realgas_source_state_admissible(qpert)) return
            end if

            call realgas_evaluate_source_state(n, iconsidered, qpert, dt_sub, source_pert, jac_pert)

            do ia = 1, nloc
               if (ia .eq. 1) then
                  idxa = dimensiona + 3
               else
                  idxa = dimensiona + 2 + ia
               end if

               amat(ia, ib) = -dt_sub * (source_pert(idxa) - source_eval(idxa)) / eps
               if (ia .eq. ib) amat(ia, ib) = amat(ia, ib) + 1.00
            end do
         end do

         do ia = 1, nloc
            rhs(ia) = -res(ia)
         end do

         call realgas_solve_dense(nloc, amat, rhs, delta, linear_ok)
         if (.not. linear_ok) return

         alpha = 1.00
         accepted = .false.
         do ls = 1, realgas_source_newton_linesearch
            qtest(1:nof_variables) = qcurr(1:nof_variables)
            do ia = 1, nloc
               if (ia .eq. 1) then
                  idxa = dimensiona + 3
               else
                  idxa = dimensiona + 2 + ia
               end if
               qtest(idxa) = qcurr(idxa) + alpha * delta(ia)
            end do

            if (realgas_source_state_admissible(qtest)) then
               call realgas_source_residual(n, iconsidered, qold, qtest, dt_sub, res_test, source_eval, jac_eval, testnorm)
               if ((testnorm .eq. testnorm) .and. (abs(testnorm) .lt. 1.0e300)) then
                  if ((testnorm .le. resnorm) .or. (alpha .le. 1.0e-3)) then
                     accepted = .true.
                     exit
                  end if
               end if
            end if
            alpha = 0.50 * alpha
         end do

         if (.not. accepted) return

         qcurr(1:nof_variables) = qtest(1:nof_variables)
      end do

      call realgas_source_residual(n, iconsidered, qold, qcurr, dt_sub, res, source_eval, jac_eval, resnorm)
      if (resnorm .le. 1.0e-7) then
         qnew(1:nof_variables) = qcurr(1:nof_variables)
         jac_diag_eff(1:nof_variables) = jac_eval(1:nof_variables)
         step_ok = .true.
      end if

   end subroutine realgas_coupled_source_newton

   subroutine realgas_source_residual(n, iconsidered, qold, qstate, dt_sub, res, source_eval, jac_eval, resnorm)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, intent(in) :: dt_sub
      real, dimension(1:gpu_max_nvar), intent(in) :: qold, qstate
      real, dimension(1:gpu_max_nvar), intent(out) :: res, source_eval, jac_eval
      real, intent(out) :: resnorm
      integer :: ia, idxa, nloc
      real :: scale

      res(1:gpu_max_nvar) = zero
      resnorm = zero
      nloc = nof_species + 1

      call realgas_evaluate_source_state(n, iconsidered, qstate, dt_sub, source_eval, jac_eval)

      do ia = 1, nloc
         if (ia .eq. 1) then
            idxa = dimensiona + 3
         else
            idxa = dimensiona + 2 + ia
         end if
         res(ia) = qstate(idxa) - qold(idxa) - dt_sub * source_eval(idxa)
         scale = max(max(abs(qold(idxa)), abs(qstate(idxa))), 1.0e-30)
         resnorm = max(resnorm, abs(res(ia)) / scale)
      end do

   end subroutine realgas_source_residual

   subroutine realgas_evaluate_source_state(n, iconsidered, qstate, dt_sub, source_eval, jac_eval)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, intent(in) :: dt_sub
      real, dimension(1:gpu_max_nvar), intent(in) :: qstate
      real, dimension(1:gpu_max_nvar), intent(out) :: source_eval, jac_eval
      call sources_realgas_state(n, qstate, source_eval, jac_eval)

   end subroutine realgas_evaluate_source_state

   subroutine realgas_solve_dense(n, a, b, x, ok)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n
      real, dimension(1:gpu_max_nvar, 1:gpu_max_nvar), intent(inout) :: a
      real, dimension(1:gpu_max_nvar), intent(inout) :: b
      real, dimension(1:gpu_max_nvar), intent(out) :: x
      logical, intent(out) :: ok
      integer :: i, j, k, pivot
      real :: pivabs, tmp, factor, sumv

      ok = .false.
      x(1:gpu_max_nvar) = zero

      do k = 1, n - 1
         pivot = k
         pivabs = abs(a(k, k))
         do i = k + 1, n
            if (abs(a(i, k)) .gt. pivabs) then
               pivot = i
               pivabs = abs(a(i, k))
            end if
         end do

         if ((pivabs .ne. pivabs) .or. (pivabs .le. 1.0e-300)) return

         if (pivot .ne. k) then
            do j = k, n
               tmp = a(k, j)
               a(k, j) = a(pivot, j)
               a(pivot, j) = tmp
            end do
            tmp = b(k)
            b(k) = b(pivot)
            b(pivot) = tmp
         end if

         do i = k + 1, n
            factor = a(i, k) / a(k, k)
            a(i, k) = zero
            do j = k + 1, n
               a(i, j) = a(i, j) - factor * a(k, j)
            end do
            b(i) = b(i) - factor * b(k)
         end do
      end do

      if ((abs(a(n, n)) .ne. abs(a(n, n))) .or. (abs(a(n, n)) .le. 1.0e-300)) return

      do i = n, 1, -1
         sumv = b(i)
         do j = i + 1, n
            sumv = sumv - a(i, j) * x(j)
         end do
         if ((abs(a(i, i)) .ne. abs(a(i, i))) .or. (abs(a(i, i)) .le. 1.0e-300)) return
         x(i) = sumv / a(i, i)
         if ((x(i) .ne. x(i)) .or. (abs(x(i)) .gt. 1.0e300)) return
      end do

      ok = .true.

   end subroutine realgas_solve_dense

   logical function realgas_source_state_admissible(qtrial)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real, dimension(1:gpu_max_nvar), intent(in) :: qtrial
      integer :: iv, k, idxe, idxev
      real :: rho, vel2, etot, evib, echem, etr, rmix, cv_mix, pressure
      real :: rhoy, species_sum, y_i
      real :: state_min_abs, state_max_abs

      realgas_source_state_admissible = .false.
      state_min_abs = 1.0e-300
      state_max_abs = 1.0e300

      if ((realgas .ne. 1) .or. (nof_species .le. 0)) then
         realgas_source_state_admissible = .true.
         return
      end if

      idxe = dimensiona + 2
      idxev = dimensiona + 3
      if (idxev + nof_species .gt. nof_variables) return

      do iv = 1, nof_variables
         if ((qtrial(iv) .ne. qtrial(iv)) .or. (abs(qtrial(iv)) .gt. state_max_abs)) return
      end do

      rho = qtrial(1)
      if ((rho .ne. rho) .or. (abs(rho) .gt. state_max_abs)) return
      if (rho .le. state_min_abs) return

      species_sum = zero
      do k = 1, nof_species
         rhoy = qtrial(idxev + k)
         if ((rhoy .ne. rhoy) .or. (abs(rhoy) .gt. state_max_abs)) return
         if (rhoy .lt. zero) return
         species_sum = species_sum + max(rhoy, zero)
      end do
      if (species_sum .le. state_min_abs) return

      if ((qtrial(idxev) .ne. qtrial(idxev)) .or. (abs(qtrial(idxev)) .gt. state_max_abs)) return
      if (qtrial(idxev) .lt. -state_min_abs) return

      vel2 = zero
      if (nof_variables .ge. 2) vel2 = vel2 + (qtrial(2) / rho) * (qtrial(2) / rho)
      if ((dimensiona .ge. 2) .and. (nof_variables .ge. 3)) vel2 = vel2 + (qtrial(3) / rho) * (qtrial(3) / rho)
      if ((dimensiona .eq. 3) .and. (nof_variables .ge. 4)) vel2 = vel2 + (qtrial(4) / rho) * (qtrial(4) / rho)
      if ((vel2 .ne. vel2) .or. (abs(vel2) .gt. state_max_abs)) return

      etot = qtrial(idxe) / rho
      evib = max(qtrial(idxev), zero) / rho
      if ((etot .ne. etot) .or. (abs(etot) .gt. state_max_abs)) return
      if ((evib .ne. evib) .or. (abs(evib) .gt. state_max_abs)) return

      echem = zero
      rmix = zero
      cv_mix = zero
      do k = 1, nof_species
         if (abs(rg_molm(k)) .le. tolsmall) return
         y_i = max(qtrial(idxev + k), zero) / rho
         if ((y_i .ne. y_i) .or. (abs(y_i) .gt. state_max_abs)) return
         if (rg_hzero(k) .gt. zero) echem = echem + (y_i * (rg_hzero(k) / rg_molm(k)))
         rmix = rmix + (y_i / rg_molm(k))
         if (k .le. 3) then
            cv_mix = cv_mix + (y_i * 2.50 * (rgs_ru / rg_molm(k)))
         else
            cv_mix = cv_mix + (y_i * 1.50 * (rgs_ru / rg_molm(k)))
         end if
      end do

      rmix = rgs_ru * rmix
      if ((rmix .le. tolsmall) .or. (cv_mix .le. tolsmall)) return

      etr = etot - evib - echem - (oo2 * vel2)
      if ((etr .ne. etr) .or. (abs(etr) .gt. state_max_abs)) return
      if (etr .le. state_min_abs) return

      pressure = rho * rmix * (etr / cv_mix)
      if ((pressure .ne. pressure) .or. (abs(pressure) .gt. state_max_abs)) return
      if (pressure .le. state_min_abs) return

      realgas_source_state_admissible = .true.

   end function realgas_source_state_admissible

   subroutine sources_computation_rot(n)
      implicit none
      integer, intent(in) :: n
      integer :: i, kmaxe

      kmaxe = xmpielrank(n)
      if ((srfg .eq. 1) .and. (kmaxe .gt. 0)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe) &
!$omp& private(i) &
!$omp& map(alloc: ielem_totvolume, rhs_val, srf_velocity, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_computation_rot_srf_cell(i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
      if ((mrf .eq. 1) .and. (kmaxe .gt. 0)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe) &
!$omp& private(i) &
!$omp& map(alloc: ielem_totvolume, rec_mrf, rec_mrf_velocity, rhs_val, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_computation_rot_mrf_cell(i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
   end subroutine sources_computation_rot

   subroutine sources_computation_rot_srf_cell(i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      integer :: iv
      real :: oodensity
      real, dimension(5) :: source_t2
      real, dimension(3) :: pox, poy
      real, dimension(3) :: rotvec

      oodensity = 1.00 / u_c_val(i, 1, 1)
      source_t2(1) = zero
      source_t2(2) = u_c_val(i, 1, 2) * oodensity - uvel
      source_t2(3) = u_c_val(i, 1, 3) * oodensity - vvel
      source_t2(4) = u_c_val(i, 1, 4) * oodensity - wvel
      source_t2(5) = zero
      pox(1:3) = source_t2(2:4)
      poy(1:3) = srf_velocity(1:3)
      call vect_function(pox, poy, rotvec)
      source_t2(2) = u_c_val(i, 1, 1) * rotvec(1)
      source_t2(3) = u_c_val(i, 1, 1) * rotvec(2)
      source_t2(4) = u_c_val(i, 1, 1) * rotvec(3)

      do iv = 1, nof_variables
         rhs_val(i, iv) = rhs_val(i, iv) + source_t2(iv) * ielem_totvolume(i)
      end do

   end subroutine sources_computation_rot_srf_cell

   subroutine sources_computation_rot_mrf_cell(i)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i
      integer :: iv
      real :: oodensity
      real, dimension(5) :: source_t2
      real, dimension(3) :: pox, poy
      real, dimension(3) :: rotvec

      if (rec_mrf(i) .eq. 1) then
         oodensity = 1.00 / u_c_val(i, 1, 1)
         source_t2(1) = zero
         source_t2(2) = u_c_val(i, 1, 2) * oodensity
         source_t2(3) = u_c_val(i, 1, 3) * oodensity
         source_t2(4) = u_c_val(i, 1, 4) * oodensity
         source_t2(5) = zero
         pox(1:3) = source_t2(2:4)
         poy(1:3) = rec_mrf_velocity(i, 1:3)
         call vect_function(pox, poy, rotvec)
         source_t2(2) = u_c_val(i, 1, 1) * rotvec(1)
         source_t2(3) = u_c_val(i, 1, 1) * rotvec(2)
         source_t2(4) = u_c_val(i, 1, 1) * rotvec(3)

         do iv = 1, nof_variables
            rhs_val(i, iv) = rhs_val(i, iv) + source_t2(iv) * ielem_totvolume(i)
         end do
      end if

   end subroutine sources_computation_rot_mrf_cell

   subroutine sources_realgas(n, iconsidered, source_r, src_jac_diag, source_dt)
      implicit none
      !> sources computation in 2d/3d for 5-species, 2-t real gas
#ifdef xpu
!$omp declare target
#endif

      integer, intent(in) :: n, iconsidered
      real, intent(in) :: source_dt
      real, intent(out) :: source_r(nof_variables), src_jac_diag(nof_variables)
      real :: qcell(gpu_max_nvar)
! Pack the strided cell-first state into a contiguous local workspace.
      qcell(1:nof_variables) = u_c_val(iconsidered, 1, 1:nof_variables)
! source_dt retained for API compatibility; physical rates do not depend on it.
      call sources_realgas_state(n, qcell, source_r, src_jac_diag)
   end subroutine sources_realgas

   subroutine sources_realgas_state(n, qstate, source_r, src_jac_diag)
      implicit none
      !> sources computation in 2d/3d for 5-species, 2-t real gas
#ifdef xpu
!$omp declare target
#endif

      integer, intent(in) :: n
      real, intent(in) :: qstate(nof_variables)
      real, dimension(1:nof_variables), intent(inout) :: source_r, src_jac_diag

      integer :: i, j, k, l, rg_i, rg_j
      real, dimension(gpu_max_species) :: rg_r, rg_rm, rg_n ! rg_r = ρ_s
      real :: rg_dq_rad, rg_qtv, rg_qw, velx, vely, velz
      real, dimension(gpu_max_species) :: rg_ev, rg_esv
      real, dimension(3) :: rg_tvs
      real, dimension(gpu_max_species) :: rg_dw_s
      real :: rg_ke1, rg_ke2, rg_ke3, rg_ke4, rg_ke5
      real :: ke1_si, ke2_si, ke3_si
      real :: rg_t, rg_ta, rg_tv
      real :: rg_z, rg_speed
      real :: rg_r1, rg_r2, rg_r3, rg_r4, rg_r5
      real :: rg_kf1, rg_kf2, rg_kf3, rg_kf4, rg_kf5
      real :: rg_kb1, rg_kb2, rg_kb3, rg_kb4, rg_kb5

      real :: rg_kf11, rg_kf12, rg_kf13, rg_kf14, rg_kf15
      real :: rg_kf21, rg_kf22, rg_kf23, rg_kf24, rg_kf25
      real :: rg_kf31, rg_kf32, rg_kf33, rg_kf34, rg_kf35
      real :: rg_kb11, rg_kb12, rg_kb13, rg_kb14, rg_kb15
      real :: rg_kb21, rg_kb22, rg_kb23, rg_kb24, rg_kb25
      real :: rg_kb31, rg_kb32, rg_kb33, rg_kb34, rg_kb35

      real :: kf11m, kf12m, kf13m, kf14m, kf15m
      real :: kf21m, kf22m, kf23m, kf24m, kf25m
      real :: kf31m, kf32m, kf33m, kf34m, kf35m
      real :: kb11m, kb12m, kb13m, kb14m, kb15m
      real :: kb21m, kb22m, kb23m, kb24m, kb25m
      real :: kb31m, kb32m, kb33m, kb34m, kb35m
      real :: kf4m, kf5m, kb4m, kb5m
      real :: kf1_eff, kb1_eff, kf2_eff, kb2_eff, kf3_eff, kb3_eff
      real :: rg_c(5) ! molar concentrations [mol/m^3]
      real :: wdot(5) ! reaction rates [mol/(m^3 s)]
      integer :: nu(5, 5) ! stoichiometric matrix ν(s,r)
      integer :: s, r

      real, dimension(gpu_max_species, gpu_max_species) :: rg_mu
      real :: temp1, temp2, taumw, taupk
      real, dimension(gpu_max_species) :: rg_tsmw, rg_tsp, rg_tssum, rg_ablot, rg_ev_eq
      real, dimension(1:gpu_max_nvar) :: leftv, tempvect
      real :: mp_pinfl, gammal, q_chem
      real :: p_atm, t_13, rg_pressure, n_tot, exp_argtv, exp_argt, rs
      integer :: rg_molx
      real, dimension(1:gpu_max_species) :: rg_mass_amu

      real, parameter :: kb = 1.380649e-23 ! boltzmann [j/k]
      real, parameter :: na = 6.02214076e23 ! avogadro [1/mol]
      real, parameter :: rg_exp_min = -700.00
      real, parameter :: rg_exp_max = 700.00

      real, dimension(1:gpu_max_species) :: tau_mw !! millikan-white [s]
      real, dimension(1:gpu_max_species) :: tau_park ! park [s]
      real, dimension(1:gpu_max_species) :: tau_tot ! combined [s]
      real :: num, den
      real :: mu_sr, a_sr, b_sr
      real :: log_p_tau, p_tau, tau_sr
      real :: sigma0, sigma_v, v_th, m_s, maxloss
      real :: rg_alpha, rg_alpha_i, rg_rho_min, sum_sourcex
      integer :: idxe, idxev, idxy1, idxn2, idx
      real :: rho_min, lambda_s, lambda_e
      real :: tau_v_mix, rhoe_local, yi, sumy, q_ttv, tv_mix, t_eff
      real :: dt_loc
      real :: rhoe_old, rhoe_new
      real :: rhoev_old, rhoev_new, rhoev_eq
      real :: alpha_v
      real :: prod, dest, lambda_k, rhoy_old, rhoy_new
      real :: tau_v_eff, kvt, rhoy_n2, rho, sumrhoy
      real :: rho_vib, b_v, lambda_v

      ! ---------------------------------------------------------------------------
      ! start
      ! ---------------------------------------------------------------------------
      i = 1

      source_r(1:nof_variables) = 0.00
      src_jac_diag(1:nof_variables) = 0.00
      tau_mw(:) = 1.0e30
      tau_park(:) = 1.0e30
      tau_tot(:) = 1.0e30

      if (nof_species .ne. 5) return

      !-----------------------------------------------------------
      ! get conservative variables and primitive variables
      !-----------------------------------------------------------
      leftv(1:nof_variables) = qstate(1:nof_variables)

      ! species densities ρ_s = ρ y_s
      do rg_i = 1, nof_species
         rg_r(rg_i) = max(leftv(dimensiona + 3 + rg_i), 0.00)
      end do

      !-----------------------------------------------------------
      ! 1) molar concentrations c_s = rho_s / m_s
      !-----------------------------------------------------------
      do s = 1, 5
         rg_c(s) = rg_r(s) / max(rg_molm(s), 1.0e-30)
      end do

      call cons2prim(n, leftv, mp_pinfl, gammal)

      rg_pressure = max(leftv(dimensiona + 2), 1.0e-12)
      p_atm = max(rg_pressure / rgs_pa_per_atm, 1.0e-30)

      ! molar masses in kg/mol -> "amu-like" for mw reduced-mass formula
      rg_mass_amu(:) = 1.0e30
      do rg_i = 1, nof_species
         rg_mass_amu(rg_i) = max(rg_molm(rg_i), 1.0e-30) * 1000.00
      end do

      tempvect(1:nof_variables) = qstate(1:nof_variables)
      call cons2div(n, tempvect, mp_pinfl, gammal)

      ! leftv:   (ρ, u, v, ρe, ρev, ρy1..ρy5)
      ! tempvect:(ρ, u, v, ttr, tv, y1..y5)

      ! temperatures
      rg_t = max(tempvect(dimensiona + 2), 50.00) ! ttr
      rg_tv = max(tempvect(dimensiona + 3), 50.00) ! tv
      rg_z = 10000.00 / rg_t
      t_13 = exp((-1.00 / 3.00) * log(rg_t))

      ! adjusted two-temperature model (candler/park style)
      rg_ta = exp(0.60 * log(rg_t) + 0.40 * log(rg_tv))

      ! flow speed (not used in source terms yet)
      velx = leftv(2)
      vely = leftv(3)
      velz = 0.00
      if (dimensiona .eq. 3) then
         velz = leftv(4)
      end if
      rg_speed = sqrt(velx**2 + vely**2 + velz**2)

      !-----------------------------------------------------------
      ! chemical source: reaction rates (candler 5-reaction model)
      ! using ρ-based rate coefficients cf [m^3/(kg s)]
      !-----------------------------------------------------------
      rg_dw_s(:) = 0.00

      call realgas_air5_rates(rg_t, rg_tv, rg_c, wdot)
      rg_dw_s(1) = rg_molm(1) * (-wdot(1) - wdot(4))
      rg_dw_s(2) = rg_molm(2) * (-wdot(2) + wdot(5))
      rg_dw_s(3) = rg_molm(3) * (-wdot(3) + wdot(4) - wdot(5))
      rg_dw_s(4) = rg_molm(4) * (2.0e0 * wdot(1) + wdot(3) + wdot(4) + wdot(5))
      rg_dw_s(5) = rg_molm(5) * (2.0e0 * wdot(2) + wdot(3) - wdot(4) - wdot(5))

      ! The active real-gas source path integrates chemistry with a local
      ! implicit Newton/subcycle solve and checks accepted states for
      ! nonnegative species. Do not additionally clip the source rates here:
      ! the old cell-local alpha limiter switches on/off discontinuously and
      ! can seed checkerboard chemistry on refined meshes.

      ! optional tiny cutoff (for cleanliness, will not break stoichiometry)
      do rg_i = 1, nof_species
         if (abs(rg_dw_s(rg_i)) < 1.0e-50) rg_dw_s(rg_i) = 0.00
      end do

      !-----------------------------------------------------------
      ! 1) millikan-white vt relaxation times (species 1..3)
      !-----------------------------------------------------------
      do rg_i = 1, 3
         num = 0.00
         den = 0.00

         do rg_j = 1, nof_species
            if (rg_r(rg_j) <= 0.00) cycle

            ! reduced mass μ_sr in amu
            mu_sr = rg_mass_amu(rg_i) * rg_mass_amu(rg_j) / (rg_mass_amu(rg_i) + rg_mass_amu(rg_j))

            ! mw coefficients
            a_sr = 0.001160 * sqrt(mu_sr) * exp((4.00 / 3.00) * log(rg_thetag(rg_i)))
            b_sr = 0.0150 * sqrt(sqrt(mu_sr))

            ! p * tau_sr in atm*s
            log_p_tau = a_sr * (t_13 - b_sr) - 18.420
            p_tau = exp(log_p_tau)

            ! tau_sr [s]: p_tau has atm*s, divide by p in atm
            tau_sr = p_tau / p_atm

            ! mixture rule with mole fractions:
            ! 1/tau_s = σ_r (x_r / tau_sr)
            ! x_r ~ n_r / σ n_r, and n_r ~ ρ_r / m_r
            num = num + (rg_r(rg_j) / rg_molm(rg_j)) / tau_sr
            den = den + (rg_r(rg_j) / rg_molm(rg_j))
         end do

         if (num > 0.00 .and. den > 0.00) then
            tau_mw(rg_i) = den / num
         else
            tau_mw(rg_i) = 1.0e30
         end if

      end do

      sigma0 = 3.0e-21 ! [m^2]

! number density [1/m^3] based on translational/rotational t
      n_tot = max(rg_pressure / (kb * rg_t), 1.0e-300) ! rg_t = t_tr

! single vibrational temperature for the mixture
      ! <-- your global vibrational temperature

! avoid non-physical tv
      tv_mix = max(rg_tv, 300.00)

! park-style weighting between t and tv (0.7 is a common choice)
      q_ttv = 0.70

! effective two-temperature average (same for all species)
! geometric mean: t_eff = t^q * tv^(1-q)
      t_eff = tv_mix * exp(q_ttv * log(rg_t / tv_mix))
      t_eff = max(t_eff, 300.00)

      do rg_i = 1, 3
         ! particle mass [kg/particle]
         m_s = max(rg_molm(rg_i), 1.0e-30) / na

         ! effective cross-section using t_eff instead of rg_t
         sigma_v = sigma0 * (50000.00 / t_eff)**2

         ! thermal speed using t_eff
         v_th = sqrt(8.00 * kb * t_eff / (pi * m_s))

         ! tau_park [s]
         tau_park(rg_i) = 1.00 / max(n_tot * sigma_v * v_th, 1.0e-300)

      end do

      !-----------------------------------------------------------
      ! 3) combine millikan-white + park collision-limited correction.
      !    Park's high-temperature correction is additive:
      !      tau_tot = tau_mw + tau_park
      !    Using a harmonic sum makes relaxation faster than both limits and
      !    drives rhoev/tvb into the numerical caps.
      !-----------------------------------------------------------

      do rg_i = 1, 3

         ! local copies
         taumw = tau_mw(rg_i)
         taupk = tau_park(rg_i)

         ! replace non-positive with 'disabled'
         if (taumw <= 0.00) taumw = 1.0e30
         if (taupk <= 0.00) taupk = 1.0e30

         ! now combine them safely
         if (taumw >= 1.0e29 .and. taupk >= 1.0e29) then
            ! both basically 'off'
            tau_tot(rg_i) = 1.0e30

         elseif (taumw >= 1.0e29) then
            ! only park valid
            tau_tot(rg_i) = taupk

         elseif (taupk >= 1.0e29) then
            ! only mw valid
            tau_tot(rg_i) = taumw

         else
            ! both valid: collision-limited additive correction
            tau_tot(rg_i) = taumw + taupk

         end if

         rg_tssum(rg_i) = tau_tot(rg_i)

      end do

      !-----------------------------------------------------------
      ! 4) vibrational energies e_v(tv) and e_v_eq(t)
      !-----------------------------------------------------------
      do rg_i = 1, 3
         rs = rgs_ru / rg_molm(rg_i)
         exp_argtv = rg_thetag(rg_i) / rg_tv
         rg_ev(rg_i) = rs * rg_thetag(rg_i) / (exp(exp_argtv) - 1.00)

         exp_argt = rg_thetag(rg_i) / rg_t
         rg_ev_eq(rg_i) = rs * rg_thetag(rg_i) / (exp(exp_argt) - 1.00)
      end do

      !-----------------------------------------------------------
      ! 5) vt relaxation source: qtv
      !    qtv = sum_i rho_i (e_v_eq - e_v) / tau_i
      !-----------------------------------------------------------
      rg_qtv = 0.00
      do rg_i = 1, 3
         if (rg_tssum(rg_i) .gt. 0.00) then
            rg_qtv = rg_qtv + rg_r(rg_i) * (rg_ev_eq(rg_i) - rg_ev(rg_i)) / rg_tssum(rg_i)
         end if
      end do

      !-----------------------------------------------------------
      ! 6) reactive vibrational source: qw
      !    qw = sum_i (dw_i * e_v(tv))  (only vibrating species)
      !-----------------------------------------------------------
      rg_qw = 0.00
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
      source_r(1:nof_variables) = 0.00

      ! Chemical formation-energy rate. The conservative total energy already
      ! includes formation enthalpy through prim2cons/cons2prim, so this is not
      ! also added to the total-energy equation.
      q_chem = 0.00
      do rg_i = 1, nof_species
         q_chem = q_chem + rg_dw_s(rg_i) * (rg_hzero(rg_i) / rg_molm(rg_i))
      end do

      idxe = dimensiona + 2 ! ρe
      idxev = dimensiona + 3 ! ρev
      idxy1 = dimensiona + 4 ! first species equation

      ! total energy equation: no chemistry source when formation enthalpy is
      ! carried in the conservative energy state.
      source_r(idxe) = 0.00

      ! vibrational energy equation: vt + reactive vibrational source
      source_r(idxev) = rg_qtv + rg_qw

      ! species equations: mass production rates
      source_r(idxy1:nof_variables) = rg_dw_s(1:nof_species)

      ! ------------------------------------------------------------
      ! diagonal source jacobian for implicit time stepping
      ! ------------------------------------------------------------
      src_jac_diag(1:nof_variables) = 0.00

      rho_min = 1.0e-20

      ! --- species diagonal terms: ds_s / d(ρ_s) ≈ -(-ω_s / ρ_s) = ω_s / ρ_s ---
      do rg_i = 1, nof_species
         if (rg_r(rg_i) > rho_min .and. rg_dw_s(rg_i) < 0.00) then
            ! local destruction timescale τ_s = -ρ_s / ω_s  (ω_s < 0)
            lambda_s = -rg_dw_s(rg_i) / rg_r(rg_i) ! [1/s]  > 0
            src_jac_diag(idxy1 + rg_i - 1) = -lambda_s
         else
            src_jac_diag(idxy1 + rg_i - 1) = 0.00
         end if
      end do

      ! --- vibrational energy diagonal: ds_ev / d(ρev) ≈ -1/τ_v_mix ---
      ! use a simple mixture vt timescale based on tau_tot (or tau_mw) for n2,o2,no
      tau_v_mix = 1.0e30
      sumy = 0.00

      do rg_i = 1, 3 ! vibrating species: n2, o2, no
         if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0e29) then
            yi = rg_r(rg_i) / leftv(1) ! y_i ≈ ρ_i/ρ
            sumy = sumy + yi
            tau_v_mix = min(tau_v_mix, tau_tot(rg_i)) ! use most restrictive vt time
         end if
      end do

      if (tau_v_mix < 1.0e29) then
         src_jac_diag(idxev) = -1.00 / tau_v_mix
      else
         src_jac_diag(idxev) = 0.00
      end if

      ! No direct total-energy chemistry source diagonal; chemistry stiffness is
      ! carried by species and vibrational-energy source diagonals.
      src_jac_diag(idxe) = 0.00

   end subroutine sources_realgas_state

   subroutine sources_realgas_pi(n, iconsidered)
      implicit none
      !> sources computation in 2d/3d for 5-species, 2-t real gas
#ifdef xpu
!$omp declare target
#endif

      integer, intent(in) :: n, iconsidered
      real, dimension(1:gpu_max_nvar) :: source_r, src_jac_diag

      integer :: i, j, k, l, rg_i, rg_j
      real, dimension(gpu_max_species) :: rg_r, rg_rm, rg_n ! rg_r = ρ_s
      real :: rg_dq_rad, rg_qtv, rg_qw, velx, vely, velz
      real, dimension(gpu_max_species) :: rg_ev, rg_esv
      real, dimension(3) :: rg_tvs
      real, dimension(gpu_max_species) :: rg_dw_s
      real :: rg_ke1, rg_ke2, rg_ke3, rg_ke4, rg_ke5
      real :: ke1_si, ke2_si, ke3_si
      real :: rg_t, rg_ta, rg_tv
      real :: rg_z, rg_speed
      real :: rg_r1, rg_r2, rg_r3, rg_r4, rg_r5
      real :: rg_kf1, rg_kf2, rg_kf3, rg_kf4, rg_kf5
      real :: rg_kb1, rg_kb2, rg_kb3, rg_kb4, rg_kb5

      real :: rg_kf11, rg_kf12, rg_kf13, rg_kf14, rg_kf15
      real :: rg_kf21, rg_kf22, rg_kf23, rg_kf24, rg_kf25
      real :: rg_kf31, rg_kf32, rg_kf33, rg_kf34, rg_kf35
      real :: rg_kb11, rg_kb12, rg_kb13, rg_kb14, rg_kb15
      real :: rg_kb21, rg_kb22, rg_kb23, rg_kb24, rg_kb25
      real :: rg_kb31, rg_kb32, rg_kb33, rg_kb34, rg_kb35

      real :: kf11m, kf12m, kf13m, kf14m, kf15m
      real :: kf21m, kf22m, kf23m, kf24m, kf25m
      real :: kf31m, kf32m, kf33m, kf34m, kf35m
      real :: kb11m, kb12m, kb13m, kb14m, kb15m
      real :: kb21m, kb22m, kb23m, kb24m, kb25m
      real :: kb31m, kb32m, kb33m, kb34m, kb35m
      real :: kf4m, kf5m, kb4m, kb5m
      real :: kf1_eff, kb1_eff, kf2_eff, kb2_eff, kf3_eff, kb3_eff
      real :: rg_c(5) ! molar concentrations [mol/m^3]
      real :: wdot(5) ! reaction rates [mol/(m^3 s)]
      integer :: nu(5, 5) ! stoichiometric matrix ν(s,r)
      integer :: s, r

      real, dimension(gpu_max_species, gpu_max_species) :: rg_mu
      real :: temp1, temp2, taumw, taupk
      real, dimension(gpu_max_species) :: rg_tsmw, rg_tsp, rg_tssum, rg_ablot, rg_ev_eq
      real, dimension(1:gpu_max_nvar) :: leftv, tempvect, tempvx
      real :: mp_pinfl, gammal, q_chem
      real :: p_atm, t_13, rg_pressure, n_tot, exp_argtv, exp_argt, rs
      integer :: rg_molx
      real, dimension(1:gpu_max_species) :: rg_mass_amu

      real, parameter :: kb = 1.380649e-23 ! boltzmann [j/k]
      real, parameter :: na = 6.02214076e23 ! avogadro [1/mol]
      real, parameter :: rg_exp_min = -700.00
      real, parameter :: rg_exp_max = 700.00

      real, dimension(1:gpu_max_species) :: tau_mw !! millikan-white [s]
      real, dimension(1:gpu_max_species) :: tau_park ! park [s]
      real, dimension(1:gpu_max_species) :: tau_tot ! combined [s]
      real :: num, den
      real :: mu_sr, a_sr, b_sr
      real :: log_p_tau, p_tau, tau_sr
      real :: sigma0, sigma_v, v_th, m_s, maxloss
      real :: rg_alpha, rg_alpha_i, rg_rho_min, sum_sourcex
      integer :: idxe, idxev, idxy1, idxn2, idx
      real :: rho_min, lambda_s, lambda_e
      real :: tau_v_mix, rhoe_local, yi, sumy, q_ttv, tv_mix, t_eff
      real :: dt_loc
      real :: rhoe_old, rhoe_new
      real :: rhoev_old, rhoev_new, rhoev_eq
      real :: alpha_v
      real :: prod, dest, lambda_k, rhoy_old, rhoy_new
      real :: tau_v_eff, kvt, rhoy_n2, rho, sumrhoy
      real :: rho_vib, b_v, lambda_v

      ! ---------------------------------------------------------------------------
      ! start
      ! ---------------------------------------------------------------------------
      i = iconsidered

      source_r(1:nof_variables) = 0.00
      src_jac_diag(1:nof_variables) = 0.00
      tau_mw(:) = 1.0e30
      tau_park(:) = 1.0e30
      tau_tot(:) = 1.0e30
      if (nof_species .ne. 5) return

      !firstly the base values

      !-----------------------------------------------------------
      ! get conservative variables and primitive variables
      !-----------------------------------------------------------
      leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)

      ! species densities ρ_s = ρ y_s
      do rg_i = 1, nof_species
         rg_r(rg_i) = max(leftv(dimensiona + 3 + rg_i), 0.00)
      end do

      !-----------------------------------------------------------
      ! 1) molar concentrations c_s = rho_s / m_s
      !-----------------------------------------------------------
      do s = 1, 5
         rg_c(s) = rg_r(s) / max(rg_molm(s), 1.0e-30)
      end do

      call cons2prim(n, leftv, mp_pinfl, gammal)

      rg_pressure = max(leftv(dimensiona + 2), 1.0e-12)
      p_atm = max(rg_pressure / rgs_pa_per_atm, 1.0e-30)

      ! molar masses in kg/mol -> "amu-like" for mw reduced-mass formula
      rg_mass_amu(:) = 1.0e30
      do rg_i = 1, nof_species
         rg_mass_amu(rg_i) = max(rg_molm(rg_i), 1.0e-30) * 1000.00
      end do

      tempvect(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
      call cons2div(n, tempvect, mp_pinfl, gammal)

      ! leftv:   (ρ, u, v, ρe, ρev, ρy1..ρy5)
      ! tempvect:(ρ, u, v, ttr, tv, y1..y5)

      ! temperatures
      rg_t = max(tempvect(dimensiona + 2), 50.00) ! ttr
      rg_tv = max(tempvect(dimensiona + 3), 50.00) ! tv
      rg_z = 10000.00 / rg_t
      t_13 = exp((-1.00 / 3.00) * log(rg_t))

      ! adjusted two-temperature model (candler/park style)
      rg_ta = exp(0.60 * log(rg_t) + 0.40 * log(rg_tv))

      ! flow speed (not used in source terms yet)
      velx = leftv(2)
      vely = leftv(3)
      velz = 0.00
      if (dimensiona .eq. 3) then
         velz = leftv(4)
      end if
      rg_speed = sqrt(velx**2 + vely**2 + velz**2)

      !-----------------------------------------------------------
      ! chemical source: reaction rates (candler 5-reaction model)
      ! using ρ-based rate coefficients cf [m^3/(kg s)]
      !-----------------------------------------------------------
      rg_dw_s(:) = 0.00

      if (rg_nof_reactions == 5) then

         !
         ! candler 5-reaction model (your original code)
         !

         ! equilibrium constants k_e(t) (candler curve fits, t = ttr)
         rg_ke1 = exp(min(max(3.8980 - 12.6110 * rg_z + 0.6830 * rg_z**2 - 0.1180 * rg_z**3 + 0.0060 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke2 = exp(min(max(1.3350 - 4.1270 * rg_z - 0.6160 * rg_z**2 + 0.0930 * rg_z**3 - 0.0050 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke3 = exp(min(max(1.5490 - 7.7840 * rg_z + 0.2280 * rg_z**2 - 0.0430 * rg_z**3 + 0.0020 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke4 = exp(min(max(2.3490 - 4.8280 * rg_z + 0.4550 * rg_z**2 - 0.0750 * rg_z**3 + 0.0040 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke5 = exp(min(max(0.2150 - 3.6520 * rg_z + 0.8430 * rg_z**2 - 0.1360 * rg_z**3 + 0.0070 * rg_z**4, rg_exp_min), &
            & rg_exp_max))

         ! --- reaction 1: n2 + m <-> 2n + m (in m^3/mol/s) ---
         kf11m = 3.78e18 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = n2
         kf12m = 3.78e18 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = o2
         kf13m = 3.78e18 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = no
         kf14m = 1.11e18 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = n
         kf15m = 1.11e18 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = o

         kb11m = kf11m / rg_ke1
         kb12m = kf12m / rg_ke1
         kb13m = kf13m / rg_ke1
         kb14m = kf14m / rg_ke1
         kb15m = kf15m / rg_ke1

         kf1_eff = kf11m * rg_c(1) + kf12m * rg_c(2) + kf13m * rg_c(3) &
            & + kf14m * rg_c(4) + kf15m * rg_c(5)
         kb1_eff = kb11m * rg_c(1) + kb12m * rg_c(2) + kb13m * rg_c(3) &
            & + kb14m * rg_c(4) + kb15m * rg_c(5)

         wdot(1) = kf1_eff * rg_c(1) - kb1_eff * rg_c(4)**2 ! n2 + m <-> 2n + m

         ! --- reaction 2: o2 + m <-> 2o + m ---
         kf21m = 2.75e16 * 1.0e-6 * (1.00 / rg_ta) * exp(-5.95e4 / rg_ta) ! m = n2
         kf22m = 2.75e16 * 1.0e-6 * (1.00 / rg_ta) * exp(-5.95e4 / rg_ta) ! m = o2
         kf23m = 2.75e16 * 1.0e-6 * (1.00 / rg_ta) * exp(-5.95e4 / rg_ta) ! m = no
         kf24m = 8.25e16 * 1.0e-6 * (1.00 / rg_ta) * exp(-5.95e4 / rg_ta) ! m = n
         kf25m = 8.25e16 * 1.0e-6 * (1.00 / rg_ta) * exp(-5.95e4 / rg_ta) ! m = o

         kb21m = kf21m / rg_ke2
         kb22m = kf22m / rg_ke2
         kb23m = kf23m / rg_ke2
         kb24m = kf24m / rg_ke2
         kb25m = kf25m / rg_ke2

         kf2_eff = kf21m * rg_c(1) + kf22m * rg_c(2) + kf23m * rg_c(3) &
            & + kf24m * rg_c(4) + kf25m * rg_c(5)
         kb2_eff = kb21m * rg_c(1) + kb22m * rg_c(2) + kb23m * rg_c(3) &
            & + kb24m * rg_c(4) + kb25m * rg_c(5)

         wdot(2) = kf2_eff * rg_c(2) - kb2_eff * rg_c(5)**2 ! o2 + m <-> 2o + m

         ! --- reaction 3: no + m <-> n + o + m ---
         kf31m = 2.30e14 * 1.0e-6 * (1.00 / sqrt(rg_ta)) * exp(-7.55e4 / rg_ta) ! m = n2
         kf32m = 2.30e14 * 1.0e-6 * (1.00 / sqrt(rg_ta)) * exp(-7.55e4 / rg_ta) ! m = o2
         kf33m = 2.30e14 * 1.0e-6 * (1.00 / sqrt(rg_ta)) * exp(-7.55e4 / rg_ta) ! m = no
         kf34m = 4.60e14 * 1.0e-6 * (1.00 / sqrt(rg_ta)) * exp(-7.55e4 / rg_ta) ! m = n
         kf35m = 4.60e14 * 1.0e-6 * (1.00 / sqrt(rg_ta)) * exp(-7.55e4 / rg_ta) ! m = o

         kb31m = kf31m / rg_ke3
         kb32m = kf32m / rg_ke3
         kb33m = kf33m / rg_ke3
         kb34m = kf34m / rg_ke3
         kb35m = kf35m / rg_ke3

         kf3_eff = kf31m * rg_c(1) + kf32m * rg_c(2) + kf33m * rg_c(3) &
            & + kf34m * rg_c(4) + kf35m * rg_c(5)
         kb3_eff = kb31m * rg_c(1) + kb32m * rg_c(2) + kb33m * rg_c(3) &
            & + kb34m * rg_c(4) + kb35m * rg_c(5)

         wdot(3) = kf3_eff * rg_c(3) - kb3_eff * rg_c(4) * rg_c(5) ! no + m <-> n + o + m

         ! --- reaction 4: n2 + o <-> no + n (bimolecular, 2-body) ---
         kf4m = 3.18e10 * 1.0e-6 * exp(0.10 * log(rg_ta)) * exp(-3.77e4 / rg_ta)
         kb4m = kf4m / rg_ke4

         wdot(4) = kf4m * rg_c(1) * rg_c(5) - kb4m * rg_c(3) * rg_c(4)

         ! --- reaction 5: no + o <-> o2 + n (bimolecular) ---
         kf5m = 2.16e5 * 1.0e-6 * exp(1.290 * log(rg_ta)) * exp(-1.922e4 / rg_ta)
         kb5m = kf5m / rg_ke5

         wdot(5) = kf5m * rg_c(3) * rg_c(5) - kb5m * rg_c(2) * rg_c(4)

         ! build stoichiometry & mass production (unchanged)
         do s = 1, 5
            do r = 1, 5
               nu(s, r) = 0
            end do
         end do

         ! r1: n2 -> 2n
         nu(1, 1) = -1
         nu(4, 1) = +2

         ! r2: o2 -> 2o
         nu(2, 2) = -1
         nu(5, 2) = +2

         ! r3: no -> n + o
         nu(3, 3) = -1
         nu(4, 3) = +1
         nu(5, 3) = +1

         ! r4: n2 + o -> no + n
         nu(1, 4) = -1
         nu(5, 4) = -1
         nu(3, 4) = +1
         nu(4, 4) = +1

         ! r5: no + o -> o2 + n
         nu(3, 5) = -1
         nu(5, 5) = -1
         nu(2, 5) = +1
         nu(4, 5) = +1

         rg_dw_s(:) = 0.00
         do r = 1, 5
            do s = 1, 5
               rg_dw_s(s) = rg_dw_s(s) + nu(s, r) * rg_molm(s) * wdot(r)
            end do
         end do

      end if

      if (rg_nof_reactions == 6) then

         !
         ! park 2-t “6-reaction” model (option a, 5 neutral reactions)
         ! same topology as candler (r1..r5), park coefficients + 3rd bodies
         !

         ! equilibrium constants (keep your candler curve-fits)
         rg_ke1 = exp(min(max(3.8980 - 12.6110 * rg_z + 0.6830 * rg_z**2 - 0.1180 * rg_z**3 + 0.0060 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke2 = exp(min(max(1.3350 - 4.1270 * rg_z - 0.6160 * rg_z**2 + 0.0930 * rg_z**3 - 0.0050 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke3 = exp(min(max(1.5490 - 7.7840 * rg_z + 0.2280 * rg_z**2 - 0.0430 * rg_z**3 + 0.0020 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke4 = exp(min(max(2.3490 - 4.8280 * rg_z + 0.4550 * rg_z**2 - 0.0750 * rg_z**3 + 0.0040 * rg_z**4, rg_exp_min), &
            & rg_exp_max))
         rg_ke5 = exp(min(max(0.2150 - 3.6520 * rg_z + 0.8430 * rg_z**2 - 0.1360 * rg_z**3 + 0.0070 * rg_z**4, rg_exp_min), &
            & rg_exp_max))

         ! Dissociation Keq curve fits are cgs concentration based.  The code uses mol/m^3.
         ke1_si = rg_ke1 * 1.0e6
         ke2_si = rg_ke2 * 1.0e6
         ke3_si = rg_ke3 * 1.0e6

         ! ---- r1: n2 + m <-> 2n + m (park 2t, table 1) ----
         !   kf = c_m t_a^{-1.6} exp(-113200/t_a), cm^3/(mol s)
         !   m = n2,o2,no : c = 7.0e21
         !   m = n,o      : c = 3.0e22
         kf11m = 7.0e21 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = n2
         kf12m = 7.0e21 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = o2
         kf13m = 7.0e21 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = no
         kf14m = 3.0e22 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = n
         kf15m = 3.0e22 * 1.0e-6 * exp(-1.60 * log(rg_ta)) * exp(-1.132e5 / rg_ta) ! m = o

         kb11m = 7.0e21 * 1.0e-6 * exp(-1.60 * log(rg_t)) * exp(-1.132e5 / rg_t) / ke1_si
         kb12m = 7.0e21 * 1.0e-6 * exp(-1.60 * log(rg_t)) * exp(-1.132e5 / rg_t) / ke1_si
         kb13m = 7.0e21 * 1.0e-6 * exp(-1.60 * log(rg_t)) * exp(-1.132e5 / rg_t) / ke1_si
         kb14m = 3.0e22 * 1.0e-6 * exp(-1.60 * log(rg_t)) * exp(-1.132e5 / rg_t) / ke1_si
         kb15m = 3.0e22 * 1.0e-6 * exp(-1.60 * log(rg_t)) * exp(-1.132e5 / rg_t) / ke1_si

         kf1_eff = kf11m * rg_c(1) + kf12m * rg_c(2) + kf13m * rg_c(3) &
            & + kf14m * rg_c(4) + kf15m * rg_c(5)
         kb1_eff = kb11m * rg_c(1) + kb12m * rg_c(2) + kb13m * rg_c(3) &
            & + kb14m * rg_c(4) + kb15m * rg_c(5)

         wdot(1) = kf1_eff * rg_c(1) - kb1_eff * rg_c(4)**2

         ! ---- r2: o2 + m <-> 2o + m (park 2t, table 1) ----
         !   kf = c_m t_a^{-1.5} exp(-59500/t_a), cm^3/(mol s)
         !   m = n2,o2,no : c = 2.0e21
         !   m = n,o      : c = 1.0e22
         kf21m = 2.0e21 * 1.0e-6 * (1.00 / (rg_ta * sqrt(rg_ta))) * exp(-5.95e4 / rg_ta) ! m = n2
         kf22m = 2.0e21 * 1.0e-6 * (1.00 / (rg_ta * sqrt(rg_ta))) * exp(-5.95e4 / rg_ta) ! m = o2
         kf23m = 2.0e21 * 1.0e-6 * (1.00 / (rg_ta * sqrt(rg_ta))) * exp(-5.95e4 / rg_ta) ! m = no
         kf24m = 1.0e22 * 1.0e-6 * (1.00 / (rg_ta * sqrt(rg_ta))) * exp(-5.95e4 / rg_ta) ! m = n
         kf25m = 1.0e22 * 1.0e-6 * (1.00 / (rg_ta * sqrt(rg_ta))) * exp(-5.95e4 / rg_ta) ! m = o

         kb21m = 2.0e21 * 1.0e-6 * (1.00 / (rg_t * sqrt(rg_t))) * exp(-5.95e4 / rg_t) / ke2_si
         kb22m = 2.0e21 * 1.0e-6 * (1.00 / (rg_t * sqrt(rg_t))) * exp(-5.95e4 / rg_t) / ke2_si
         kb23m = 2.0e21 * 1.0e-6 * (1.00 / (rg_t * sqrt(rg_t))) * exp(-5.95e4 / rg_t) / ke2_si
         kb24m = 1.0e22 * 1.0e-6 * (1.00 / (rg_t * sqrt(rg_t))) * exp(-5.95e4 / rg_t) / ke2_si
         kb25m = 1.0e22 * 1.0e-6 * (1.00 / (rg_t * sqrt(rg_t))) * exp(-5.95e4 / rg_t) / ke2_si

         kf2_eff = kf21m * rg_c(1) + kf22m * rg_c(2) + kf23m * rg_c(3) &
            & + kf24m * rg_c(4) + kf25m * rg_c(5)
         kb2_eff = kb21m * rg_c(1) + kb22m * rg_c(2) + kb23m * rg_c(3) &
            & + kb24m * rg_c(4) + kb25m * rg_c(5)

         wdot(2) = kf2_eff * rg_c(2) - kb2_eff * rg_c(5)**2

         ! ---- r3: no + m <-> n + o + m (park/koshi) ----
         !   kf ≈ c_m t_a^{0} exp(-75500/t_a), cm^3/(mol s)
         !   m = no,o,n : c = 1.1e17
         !   m = n2,o2 : c = (1/22)*1.1e17  [koshi]
         kf31m = (1.1e17 / 22.00) * 1.0e-6 * exp(-7.55e4 / rg_ta) ! m = n2
         kf32m = (1.1e17 / 22.00) * 1.0e-6 * exp(-7.55e4 / rg_ta) ! m = o2
         kf33m = 1.1e17 * 1.0e-6 * exp(-7.55e4 / rg_ta) ! m = no
         kf34m = 1.1e17 * 1.0e-6 * exp(-7.55e4 / rg_ta) ! m = n
         kf35m = 1.1e17 * 1.0e-6 * exp(-7.55e4 / rg_ta) ! m = o

         kb31m = (1.1e17 / 22.00) * 1.0e-6 * exp(-7.55e4 / rg_t) / ke3_si
         kb32m = (1.1e17 / 22.00) * 1.0e-6 * exp(-7.55e4 / rg_t) / ke3_si
         kb33m = 1.1e17 * 1.0e-6 * exp(-7.55e4 / rg_t) / ke3_si
         kb34m = 1.1e17 * 1.0e-6 * exp(-7.55e4 / rg_t) / ke3_si
         kb35m = 1.1e17 * 1.0e-6 * exp(-7.55e4 / rg_t) / ke3_si

         kf3_eff = kf31m * rg_c(1) + kf32m * rg_c(2) + kf33m * rg_c(3) &
            & + kf34m * rg_c(4) + kf35m * rg_c(5)
         kb3_eff = kb31m * rg_c(1) + kb32m * rg_c(2) + kb33m * rg_c(3) &
            & + kb34m * rg_c(4) + kb35m * rg_c(5)

         wdot(3) = kf3_eff * rg_c(3) - kb3_eff * rg_c(4) * rg_c(5)

         ! ---- r4: n2 + o <-> no + n (park exchange, table 1) ----
         !   kf = 6.4e7 t^{-1} exp(-38370/t), cm^3/(mol s)
         kf4m = 6.4e7 * 1.0e-6 * (1.00 / rg_t) * exp(-3.837e4 / rg_t)
         kb4m = kf4m / rg_ke4

         wdot(4) = kf4m * rg_c(1) * rg_c(5) - kb4m * rg_c(3) * rg_c(4)

         ! ---- r5: no + o <-> o2 + n (park exchange, table 1) ----
         !   kf = 8.4e12 exp(-19450/t), cm^3/(mol s)
         kf5m = 8.4e12 * 1.0e-6 * exp(-1.945e4 / rg_t)
         kb5m = kf5m / rg_ke5

         wdot(5) = kf5m * rg_c(3) * rg_c(5) - kb5m * rg_c(2) * rg_c(4)

         ! stoichiometry is the same as candler (same 5 reactions)
         do s = 1, 5
            do r = 1, 5
               nu(s, r) = 0
            end do
         end do

         ! r1: n2 -> 2n
         nu(1, 1) = -1
         nu(4, 1) = +2

         ! r2: o2 -> 2o
         nu(2, 2) = -1
         nu(5, 2) = +2

         ! r3: no -> n + o
         nu(3, 3) = -1
         nu(4, 3) = +1
         nu(5, 3) = +1

         ! r4: n2 + o -> no + n
         nu(1, 4) = -1
         nu(5, 4) = -1
         nu(3, 4) = +1
         nu(4, 4) = +1

         ! r5: no + o -> o2 + n
         nu(3, 5) = -1
         nu(5, 5) = -1
         nu(2, 5) = +1
         nu(4, 5) = +1

         rg_dw_s(:) = 0.00
         do r = 1, 5
            do s = 1, 5
               rg_dw_s(s) = rg_dw_s(s) + nu(s, r) * rg_molm(s) * wdot(r)
            end do
         end do

      end if

      rg_rho_min = 1.0e-20 ! or whatever you use as "zero" density
      rg_alpha = 1.00
      dt_loc = max(ielem_dtl(iconsidered), 1.0e-30)

      ! first: compute most restrictive scaling factor alpha
      do rg_i = 1, nof_species

         ! if density is essentially zero, don't allow further destruction
         if (rg_r(rg_i) <= rg_rho_min .and. rg_dw_s(rg_i) < 0.00) then
            rg_alpha = 0.00
            cycle
         end if

         ! only destruction can cause negativity
         if (rg_dw_s(rg_i) < 0.00) then
            ! Keep the point-implicit source update positive while preserving
            ! stoichiometry by scaling the full production vector.
            rg_alpha_i = 0.90 * rg_r(rg_i) / max(-dt_loc * rg_dw_s(rg_i), 1.0e-300)
            rg_alpha = min(rg_alpha, rg_alpha_i)
         end if
      end do

      ! clamp alpha to [0,1]
      if (rg_alpha > 1.00) rg_alpha = 1.00
      if (rg_alpha < 0.00) rg_alpha = 0.00

      ! second: apply scaling if needed
      if (rg_alpha < 1.00) then
         do rg_i = 1, nof_species
            rg_dw_s(rg_i) = rg_alpha * rg_dw_s(rg_i)
         end do
      end if

      ! optional tiny cutoff (for cleanliness, will not break stoichiometry)
      do rg_i = 1, nof_species
         if (abs(rg_dw_s(rg_i)) < 1.0e-50) rg_dw_s(rg_i) = 0.00
      end do

      !-----------------------------------------------------------
      ! 1) millikan-white vt relaxation times (species 1..3)
      !-----------------------------------------------------------
      do rg_i = 1, 3
         num = 0.00
         den = 0.00

         do rg_j = 1, nof_species
            if (rg_r(rg_j) <= 0.00) cycle

            ! reduced mass μ_sr in amu
            mu_sr = rg_mass_amu(rg_i) * rg_mass_amu(rg_j) / (rg_mass_amu(rg_i) + rg_mass_amu(rg_j))

            ! mw coefficients
            a_sr = 0.001160 * sqrt(mu_sr) * exp((4.00 / 3.00) * log(rg_thetag(rg_i)))
            b_sr = 0.0150 * sqrt(sqrt(mu_sr))

            ! p * tau_sr in atm*s
            log_p_tau = a_sr * (t_13 - b_sr) - 18.420
            p_tau = exp(log_p_tau)

            ! tau_sr [s]: p_tau has atm*s, divide by p in atm
            tau_sr = p_tau / p_atm

            ! mixture rule with mole fractions:
            ! 1/tau_s = σ_r (x_r / tau_sr)
            ! x_r ~ n_r / σ n_r, and n_r ~ ρ_r / m_r
            num = num + (rg_r(rg_j) / rg_molm(rg_j)) / tau_sr
            den = den + (rg_r(rg_j) / rg_molm(rg_j))
         end do

         if (num > 0.00 .and. den > 0.00) then
            tau_mw(rg_i) = den / num
         else
            tau_mw(rg_i) = 1.0e30
         end if

      end do

      sigma0 = 3.0e-21 ! [m^2]

! number density [1/m^3] based on translational/rotational t
      n_tot = max(rg_pressure / (kb * rg_t), 1.0e-300) ! rg_t = t_tr

! single vibrational temperature for the mixture
      ! <-- your global vibrational temperature

! avoid non-physical tv
      tv_mix = max(rg_tv, 300.00)

! park-style weighting between t and tv (0.7 is a common choice)
      q_ttv = 0.70

! effective two-temperature average (same for all species)
! geometric mean: t_eff = t^q * tv^(1-q)
      t_eff = tv_mix * exp(q_ttv * log(rg_t / tv_mix))
      t_eff = max(t_eff, 300.00)

      do rg_i = 1, 3
         ! particle mass [kg/particle]
         m_s = max(rg_molm(rg_i), 1.0e-30) / na

         ! effective cross-section using t_eff instead of rg_t
         sigma_v = sigma0 * (50000.00 / t_eff)**2

         ! thermal speed using t_eff
         v_th = sqrt(8.00 * kb * t_eff / (pi * m_s))

         ! tau_park [s]
         tau_park(rg_i) = 1.00 / max(n_tot * sigma_v * v_th, 1.0e-300)

      end do

      !-----------------------------------------------------------
      ! 3) combine millikan-white + park collision-limited correction.
      !    Park's high-temperature correction is additive:
      !      tau_tot = tau_mw + tau_park
      !    Using a harmonic sum makes relaxation faster than both limits and
      !    drives rhoev/tvb into the numerical caps.
      !-----------------------------------------------------------

      do rg_i = 1, 3

         ! local copies
         taumw = tau_mw(rg_i)
         taupk = tau_park(rg_i)

         ! replace non-positive with 'disabled'
         if (taumw <= 0.00) taumw = 1.0e30
         if (taupk <= 0.00) taupk = 1.0e30

         ! now combine them safely
         if (taumw >= 1.0e29 .and. taupk >= 1.0e29) then
            ! both basically 'off'
            tau_tot(rg_i) = 1.0e30

         elseif (taumw >= 1.0e29) then
            ! only park valid
            tau_tot(rg_i) = taupk

         elseif (taupk >= 1.0e29) then
            ! only mw valid
            tau_tot(rg_i) = taumw

         else
            ! both valid: collision-limited additive correction
            tau_tot(rg_i) = taumw + taupk

         end if

         rg_tssum(rg_i) = tau_tot(rg_i)

      end do

      !-----------------------------------------------------------
      ! 4) vibrational energies e_v(tv) and e_v_eq(t)
      !-----------------------------------------------------------
      do rg_i = 1, 3
         rs = rgs_ru / rg_molm(rg_i)
         exp_argtv = rg_thetag(rg_i) / rg_tv
         rg_ev(rg_i) = rs * rg_thetag(rg_i) / (exp(exp_argtv) - 1.00)

         exp_argt = rg_thetag(rg_i) / rg_t
         rg_ev_eq(rg_i) = rs * rg_thetag(rg_i) / (exp(exp_argt) - 1.00)
      end do

      !-----------------------------------------------------------
      ! 5) vt relaxation source: qtv
      !    qtv = sum_i rho_i (e_v_eq - e_v) / tau_i
      !-----------------------------------------------------------
      rg_qtv = 0.00
      do rg_i = 1, 3
         if (rg_tssum(rg_i) .gt. 0.00) then
            rg_qtv = rg_qtv + rg_r(rg_i) * (rg_ev_eq(rg_i) - rg_ev(rg_i)) / rg_tssum(rg_i)
         end if
      end do

      !-----------------------------------------------------------
      ! 6) reactive vibrational source: qw
      !    qw = sum_i (dw_i * e_v(tv))  (only vibrating species)
      !-----------------------------------------------------------
      rg_qw = 0.00
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
      source_r(1:nof_variables) = 0.00

      ! Chemical formation-energy rate. The conservative total energy already
      ! includes formation enthalpy through prim2cons/cons2prim, so this is not
      ! also added to the total-energy equation.
      q_chem = 0.00
      do rg_i = 1, nof_species
         q_chem = q_chem + rg_dw_s(rg_i) * (rg_hzero(rg_i) / rg_molm(rg_i))
      end do

      idxe = dimensiona + 2 ! ρe
      idxev = dimensiona + 3 ! ρev
      idxy1 = dimensiona + 4 ! first species equation

      ! total energy equation: no chemistry source when formation enthalpy is
      ! carried in the conservative energy state.
      source_r(idxe) = 0.00

      ! vibrational energy equation: vt + reactive vibrational source
      source_r(idxev) = rg_qtv + rg_qw

      ! species equations: mass production rates
      source_r(idxy1:nof_variables) = rg_dw_s(1:nof_species)

      ! ------------------------------------------------------------
      ! diagonal source jacobian for implicit time stepping
      ! ------------------------------------------------------------
      src_jac_diag(1:nof_variables) = 0.00

      rho_min = 1.0e-20

      ! --- species diagonal terms: ds_s / d(ρ_s) ≈ -(-ω_s / ρ_s) = ω_s / ρ_s ---
      do rg_i = 1, nof_species
         if (rg_r(rg_i) > rho_min .and. rg_dw_s(rg_i) < 0.00) then
            ! local destruction timescale τ_s = -ρ_s / ω_s  (ω_s < 0)
            lambda_s = -rg_dw_s(rg_i) / rg_r(rg_i) ! [1/s]  > 0
            src_jac_diag(idxy1 + rg_i - 1) = -lambda_s
         else
            src_jac_diag(idxy1 + rg_i - 1) = 0.00
         end if
      end do

      ! --- vibrational energy diagonal: ds_ev / d(ρev) ≈ -1/τ_v_mix ---
      ! use a simple mixture vt timescale based on tau_tot (or tau_mw) for n2,o2,no
      tau_v_mix = 1.0e30
      sumy = 0.00

      do rg_i = 1, 3 ! vibrating species: n2, o2, no
         if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0e29) then
            yi = rg_r(rg_i) / leftv(1) ! y_i ≈ ρ_i/ρ
            sumy = sumy + yi
            tau_v_mix = min(tau_v_mix, tau_tot(rg_i)) ! use most restrictive vt time
         end if
      end do

      if (tau_v_mix < 1.0e29) then
         src_jac_diag(idxev) = -1.00 / tau_v_mix
      else
         src_jac_diag(idxev) = 0.00
      end if

      ! No direct total-energy chemistry source diagonal; chemistry stiffness is
      ! carried by species and vibrational-energy source diagonals.
      src_jac_diag(idxe) = 0.00

      if ((rungekutta .eq. 5) .or. (rungekutta .ge. 10)) then
         dt_loc = max(ielem_dtl(i), 1.0e-30)
      else
         dt_loc = max(dt, 1.0e-30)
      end if

      ! ----------------------------
      ! (a) total energy: chemistry only (explicit, same as before)
      ! ----------------------------
      ! original pathway:
      !   rhs_e -= q_chem * v
      !   rhoe  -= dt * (rhs_e / v)
      ! ==> rhoe += dt * q_chem
      !
      rhoe_old = u_c_val(i, 1, dimensiona + 2)
      rhoe_new = rhoe_old + dt_loc * q_chem
      u_c_val(i, 1, dimensiona + 2) = rhoe_new

      ! ----------------------------
      ! (b) vibrational energy: vt relaxation + reactive vib source
      ! ----------------------------
      ! use be-stable form:
      !   d(rhoev)/dt = (rhoev_eq(t) - rhoev)/tau_v_eff + qw
      ! relaxation part implicit, qw explicit.
      !

      rhoev_old = u_c_val(i, 1, dimensiona + 3)

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

      rho_vib = 0.00
      b_v = 0.00
      kvt = 0.00

      do rg_i = 1, 3
         if (rg_r(rg_i) > rho_min .and. tau_tot(rg_i) < 1.0e29) then
            rho_vib = rho_vib + rg_r(rg_i)
            b_v = b_v + rg_r(rg_i) * rg_ev_eq(rg_i) / tau_tot(rg_i)
            kvt = kvt + rg_r(rg_i) / tau_tot(rg_i)
         end if
      end do

      if (rho_vib > rho_min .and. kvt > 0.00) then
         lambda_v = kvt / rho_vib
         alpha_v = dt_loc * lambda_v
         rhoev_new = (rhoev_old + dt_loc * (b_v + rg_qw)) / (1.00 + alpha_v)
      else
         rhoev_new = rhoev_old + dt_loc * rg_qw ! no vt exchange, only reactive vib
      end if

      if (rhoev_new < 0.00) rhoev_new = 0.00
      u_c_val(i, 1, dimensiona + 3) = rhoev_new

      ! ----------------------------
      ! (c) species chemistry
      ! ----------------------------
      ! The source vector has already been scaled by one cell-local alpha,
      ! so applying the same increment to every species keeps the reaction
      ! stoichiometry intact.  Do not use one species as a closure reservoir:
      ! that can preserve sum(Y)=1 while creating pure-species cells.
      rho = u_c_val(i, 1, 1)
      if (rho <= rho_min) then
         rho = rho_min
         u_c_val(i, 1, 1) = rho
      end if

      sumrhoy = 0.00
      do k = 1, nof_species
         idx = dimensiona + 3 + k
         rhoy_new = u_c_val(i, 1, idx) + dt_loc * rg_dw_s(k)
         if (rhoy_new < 0.00) rhoy_new = 0.00
         u_c_val(i, 1, idx) = rhoy_new
         sumrhoy = sumrhoy + rhoy_new
      end do

      if (sumrhoy > rho_min) then
         do k = 1, nof_species
            idx = dimensiona + 3 + k
            u_c_val(i, 1, idx) = rho * u_c_val(i, 1, idx) / sumrhoy
         end do
      else
         do k = 1, nof_species
            idx = dimensiona + 3 + k
            u_c_val(i, 1, idx) = rho * rg_vf(k)
         end do
      end if

   end subroutine sources_realgas_pi

   subroutine sources_derivatives_computation(n)
      implicit none
!> @brief
!> sources derivative computation for implicit time stepping
      integer, intent(in) :: n
      integer :: i, kmaxe

      kmaxe = xmpielrank(n)
      if ((turbulence .eq. 1) .and. (kmaxe .gt. 0)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe, n) &
!$omp& private(i) &
!$omp& map(alloc: ielem_totvolume, ielem_walldist, ielem_walltrans, mp_janaf, mp_m, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rec_grads, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb, sht, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_derivatives_computation_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
   end subroutine sources_derivatives_computation

   subroutine sources_derivatives_computation_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      integer :: iv
      real, dimension(gpu_max_turbulence) :: source_t

      call sources_derivatives(n, iconsidered, source_t)
      do iv = 1, turbulenceequations
         sht(iconsidered, iv) = source_t(iv) * ielem_totvolume(iconsidered)
      end do

   end subroutine sources_derivatives_computation_cell

   real function transition_source_factor(iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: iconsidered
      real :: xtr

      transition_source_factor = 1.00
      if (transition_model .ne. 1) return

      xtr = real(transition_direction) * ielem_walltrans(iconsidered)

      if (transition_ramp_length .le. 0.00) then
         if (xtr .lt. 0.00) transition_source_factor = 0.00
      else
         xtr = xtr / transition_ramp_length
         if (xtr .le. 0.00) then
            transition_source_factor = 0.00
         else if (xtr .ge. 1.00) then
            transition_source_factor = 1.00
         else if (transition_ramp_type .eq. 2) then
            transition_source_factor = xtr * xtr * (3.00 - 2.00 * xtr)
         else
            transition_source_factor = xtr
         end if
      end if

   end function transition_source_factor

   subroutine sources_sa_rc(n, i, src, jac)
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i
      real, intent(out) :: src, jac
      real :: q(gpu_max_nvar), mu, conductivity, gu(3,3), w(3,3), om, gsq, sf
      real :: strain, distance, delta, rd, fd, fr1, rotation_delta, mut, sigk, sigw, packed(20)
      ! Shared SA-noft2 algebra: rot_corr=0 (SA), 1 (legacy SA-R), 2 (SA-RC).
      q(1:nof_variables) = u_c_val(i,1,1:nof_variables)
      call molecular_transport(n, q, mu, conductivity)
      gu = 0.0; gu(1:dimensiona,1:dimensiona) = rec_grads(i,1:dimensiona,1:dimensiona)
      if (dimensiona==2 .and. axisymmetric_enabled==1) &
         gu(3,3)=q(3)/(q(1)*axis_radius(i))
      w = 0.5*(gu-transpose(gu)); om = sqrt(2.0*sum(w*w))
      strain = sqrt(0.5*sum((gu+transpose(gu))**2))
      gsq = sum(rec_grads(i,nof_variables,1:dimensiona)**2)
      fr1 = 1.0; rotation_delta = 0.0
      if (rot_corr == 2) fr1 = sa_rc_fr1(i)
      if (rot_corr == 1) rotation_delta = 2.0*min(0.0,strain-om)
      distance = max(ielem_walldist(i),1.0e-11)
      if (des_model /= 0) then
         delta = ielem_totvolume(i)**(1.0/real(dimensiona))
         if (des_model == 1) distance = min(distance,c_des_sa*delta)
         if (des_model == 2) then
            packed = 0.0; packed(2) = u_ct_val(i,1,1)
            call eddy_viscosity_state(dimensiona,q(1),mu,packed,.true.,mut,sigk,sigw)
            ! DDES shielding uses kinematic viscosity and the full velocity-gradient norm.
            rd = min(10.0,(mu+mut)/(q(1)*max(sqrt(sum(gu*gu)),1.0e-15)*(kappa*distance)**2))
            fd = 1.0-tanh((8.0*rd)**3)
            distance = max(distance-fd*max(distance-c_des_sa*delta,0.0),1.0e-11)
         end if
      end if
      ! DES distance and reconstructed gradients are frozen in this local Jacobian.
      call sa_noft2_rc_source(q(1), mu, u_ct_val(i,1,1), om, distance, gsq, fr1, &
                             1.0e-12*mu/q(1), cb1, cb2, cv1, cw1, cw2, cw3, kappa, sigma, src, jac, rotation_delta)
      sf = transition_source_factor(i)
      src = sf*src; jac = sf*jac
   end subroutine sources_sa_rc
   subroutine sources(n, iconsidered, source_t)
      implicit none
!> @brief
!> sources computation procedure
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, dimension(turbulenceequations), intent(inout) :: source_t
      real :: alpha_star0_local
      real :: intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
      real :: vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
      real :: tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
      integer :: i, k, j, l, ihgt, ihgj, iex, lowre
      real :: snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
      real :: gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
      real :: r_des, f_des, ddw, f_des_sst, l_t_des
      real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
      real :: uz, vz, wx, wy, wz
      real :: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz !for sas only
      real, dimension(3, 3) :: vortet, tvort, svort, ovort
!declarations for k-omega
      real :: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
      real :: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
      real :: alpha_raw, alpha_star, re_t_sst, alpha_inf
      real :: beta_stari, beta_i, beta_raw, beta_star
      real :: k_0, om_0, wally
      real :: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
      real :: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
      real :: kx, ky, kz, omx, omy, omz, rho, nu, nu_tilde, difterm, prod
      real, dimension(1:gpu_max_nvar) :: leftv, rightv
      real :: mp_pinfl, gammal
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm

      real :: rc_src, rc_jac

      i = iconsidered
      if (turbulencemodel == 1) then
         call sources_sa_rc(n, i, rc_src, rc_jac)
         source_t(1) = rc_src
         return
      end if


      verysmall = 10e-16

      vortet(1:3, 1:3) = rec_grads(i, 1:3, 1:3)

      ux = vortet(1, 1); uy = vortet(1, 2); uz = vortet(1, 3)
      vx = vortet(2, 1); vy = vortet(2, 2); vz = vortet(2, 3)
      wx = vortet(3, 1); wy = vortet(3, 2); wz = vortet(3, 3)

      do ihgt = 1, 3
         do ihgj = 1, 3
            tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
         end do
      end do

      svort = 0.5 * (vortet + tvort)
      ovort = 0.5 * (vortet - tvort)

      snorm = sqrt(2.00 * ((svort(1, 1) * svort(1, 1)) + (svort(1, 2) * svort(1, 2)) + (svort(1, 3) * svort(1, 3)) + &
         & (svort(2, 1) * svort(2, 1)) + (svort(2, 2) * svort(2, 2)) + (svort(2, 3) * svort(2, 3)) + &
         & (svort(3, 1) * svort(3, 1)) + (svort(3, 2) * svort(3, 2)) + (svort(3, 3) * svort(3, 3))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
      onorm = sqrt(2.00 * ((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + (ovort(1, 3) * ovort(1, 3)) + &
         & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2)) + (ovort(2, 3) * ovort(2, 3)) + &
         & (ovort(3, 1) * ovort(3, 1)) + (ovort(3, 2) * ovort(3, 2)) + (ovort(3, 3) * ovort(3, 3))))
      omega = onorm

      divnorm = ux + vy + wz !careful with the sign. if it becomes very big, it can produce negative production

      usss = sqrt((2.0 * ((ux * ux) + (vy * vy) + (wz * wz))) &
         & + ((uy + vx) * (uy + vx) + (uz + wx) * (uz + wx) + (wy + vz) * (wy + vz)) &
         & - (2.0 / 3.0 * (ux + vy + wz) * (ux + vy + wz)))

      squaret = (sqrt((rec_grads(i, nof_variables, 1)**2) + (rec_grads(i, nof_variables, 2)**2) &
         & + (rec_grads(i, nof_variables, 3)**2)))**2

      leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
      rightv(1:nof_variables) = leftv(1:nof_variables)

      call molecular_transport(n, leftv, viscl(1), laml(1))
      viscl(2) = viscl(1); laml(2) = laml(1)
      viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0

      turbmv(1) = u_ct_val(i, 1, 1)
      turbmv(2) = turbmv(1)

      select case (turbulencemodel)

      case (2) !k omega sst

         eddyfl(1) = ielem_walldist(i)
         eddyfl(2) = u_ct_val(i, 1, 1)
         eddyfl(3) = u_ct_val(i, 1, 2)
         call fill_turbulence_gradients(3, i, eddyfl)


         eddyfr = eddyfl

         call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
         k_0 = (max(verysmall, u_ct_val(i, 1, 1) / leftv(1))) !first subindex makes reference to the time-stepping
         om_0 = max(1.0e-6, u_ct_val(i, 1, 2) / leftv(1))
         wally = max(ielem_walldist(i),1.0e-15)

         !calculate here k and omega gradients
         omx = rec_grads(i, nof_variables + 1, 1)
         omy = rec_grads(i, nof_variables + 1, 2)
         omz = rec_grads(i, nof_variables + 1, 3)
         kx = rec_grads(i, nof_variables, 1)
         ky = rec_grads(i, nof_variables, 2)
         kz = rec_grads(i, nof_variables, 3)

         dervk_dervom = kx * omx + ky * omy + kz * omz

         call sst_blending(leftv(1), viscl(1), k_0, om_0, wally, dervk_dervom, f_1, f_2)
         alpha_raw = f_1 * alpha_inf1 + (1.0-f_1) * alpha_inf2
         beta_star = beta_starinf
         beta_raw = f_1 * beta_i1 + (1.0-f_1) * beta_i2

         ! SST-2003m: use the same limited production in both equations.
         ! vort_model=1 remains the explicitly selected Kato-Launder variant.
         if (vort_model == 1) then
            prod_k = viscl(3) * snorm * onorm
         else
            prod_k = viscl(3) * snorm**2
         end if
         prod_k = min(prod_k, 10.0 * leftv(1) * beta_star * k_0 * om_0)
         prod_om = alpha_raw * prod_k * leftv(1) / max(viscl(3),1.0e-30)

         !destruction terms
         !destruction of k
         ydest_k = leftv(1) * beta_star * k_0 * om_0
         !destruction of omega
         ydest_om = leftv(1) * beta_raw * om_0 * om_0

         !crossed-diffusion term
         !crossed diffusion of omega
         diff_om = 2.0 * (1 - f_1) * leftv(1) / (om_0 * sigma_om2) * dervk_dervom

         !qsas term: scale adaptive
         if (qsas_model .eq. 1) then !<------------!!!!!!!!!!!!!!!!!!!
            !calculate here second derivative of u
            !declare all variables
            do iex = 1, 3
               vortet(iex, 1:3) = rec_grads(i, 3 + turbulenceequations + iex, 1:3)

            end do

            uxx = vortet(1, 1); uyy = vortet(1, 2); uzz = vortet(1, 3)
            vxx = vortet(2, 1); vyy = vortet(2, 2); vzz = vortet(2, 3)
            wxx = vortet(3, 1); wyy = vortet(3, 2); wzz = vortet(3, 3)
            !compute here cell volume
            cell_volume = ielem_totvolume(i)

            u_lapl = sqrt((uxx + uyy + uzz)**2 + (vxx + vyy + vzz)**2 + (wxx + wyy + wzz)**2)
            dervk2 = kx * kx + ky * ky + kz * kz
            dervom2 = omx * omx + omy * omy + omz * omz

            delta_cell = exp((1.00 / 3.00) * log(cell_volume))

            l_sas = sqrt(k_0) / (sqrt(sqrt(beta_star)) * om_0)
            l_vk = max(kappa * snorm / u_lapl, & !this switch provides high wave-number damping
               & c_smg * delta_cell * sqrt(kappa * eta2_sas / (beta_raw / beta_star - alpha_raw)))

            q_sas1 = leftv(1) * eta2_sas * kappa * snorm**2 * (l_sas / l_vk)**2
            q_sas2 = -c_sas * 2 * leftv(1) * k_0 / sigma_phi * max(dervk2 / k_0**2, dervom2 / om_0**2)

            q_sas = max(q_sas1 + q_sas2, 0.0)
         else
            q_sas = zero
         end if

         !final source terms

         !des-sst model (if qsas_model=2)
         if (qsas_model .eq. 2) then
            cell_volume = ielem_totvolume(i)
            delta_cell = exp((1.00 / 3.00) * log(cell_volume))
            l_t_des = sqrt(k_0) / (beta_star * om_0)
            f_des_sst = max(1.0, l_t_des / (c_des_sst * delta_cell) * (1 - f_2))
            !the (1-f_2) is meant to protect the boundary layer. will result in same
            !separation point that standard s-a
            ydest_k = f_des_sst * ydest_k
         end if

         srcfull_k = prod_k - ydest_k
         srcfull_om = prod_om - ydest_om + diff_om + q_sas

         !filling the output vector
         source_t(1) = srcfull_k
         source_t(2) = srcfull_om

      end select

   end subroutine sources

   subroutine sst_destruction_diagonal(n, i, nd, diagonal)
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, i, nd
      real, intent(out) :: diagonal(2)
      real :: q(gpu_max_nvar), mu, conductivity, rho, kin, omega, cross, f1, f2, beta_sst, factor, delta
      integer :: d
      q(1:nof_variables) = u_c_val(i,1,1:nof_variables)
      rho = max(q(1),1.0e-30)
      kin = max(u_ct_val(i,1,1)/rho,1.0e-15)
      omega = max(u_ct_val(i,1,2)/rho,1.0e-6)
      call molecular_transport(n, q, mu, conductivity)
      cross = 0.0
      do d=1,nd
         cross = cross + rec_grads(i,nof_variables,d)*rec_grads(i,nof_variables+1,d)
      end do
      call sst_blending(rho, mu, kin, omega, ielem_walldist(i), cross, f1, f2)
      beta_sst = f1*beta_i1+(1.0-f1)*beta_i2
      factor = 1.0
      if (qsas_model == 2) then
         delta = ielem_totvolume(i)**(1.0/3.0)
         factor = max(1.0, sqrt(kin)/(beta_starinf*omega*c_des_sst*delta)*(1.0-f2))
      end if
      ! Freeze blending, production, cross diffusion and optional DES factor.
      ! These negative terms add damping when subtracted from the implicit matrix.
      diagonal(1) = -factor*beta_starinf*omega
      diagonal(2) = -2.0*beta_sst*omega
   end subroutine sst_destruction_diagonal

   subroutine sources_derivatives(n, iconsidered, source_t)
      implicit none
!> @brief
!> sources derivative computation for implicit time stepping
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, dimension(turbulenceequations), intent(inout) :: source_t
      real :: alpha_star0_local
      real :: intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
      real :: vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
      real :: tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
      integer :: i, k, j, l, ihgt, ihgj, iex, lowre
      real :: snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
      real :: gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
      real :: r_des, f_des, ddw, f_des_sst, l_t_des
      real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
      real :: uz, vz, wx, wy, wz
      real :: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz !for sas only
      real, dimension(3, 3) :: vortet, tvort, svort, ovort
!declarations for k-omega
      real :: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
      real :: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
      real :: alpha_raw, alpha_star, re_t_sst, alpha_inf
      real :: beta_stari, beta_i, beta_raw, beta_star
      real :: k_0, om_0, wally
      real :: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
      real :: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2, rho, nu, nu_tilde, difterm, prod, dstild_dnut, dest, dfw_dgg
      real :: kx, ky, kz, omx, omy, omz, dprod_dnut, ddif_dnut, ddest_dnut, dif, drr_dnut, dfw_drr, dgg_drr
      real, dimension(1:gpu_max_nvar) :: leftv, rightv
      real :: mp_pinfl, gammal
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm

      real :: rc_src, rc_jac

      i = iconsidered
      if (turbulencemodel == 1) then
         call sources_sa_rc(n, i, rc_src, rc_jac)
         source_t(1) = rc_jac
         return
      end if


      verysmall = 10e-16

      vortet(1:3, 1:3) = rec_grads(i, 1:3, 1:3)

      ux = vortet(1, 1); uy = vortet(1, 2); uz = vortet(1, 3)
      vx = vortet(2, 1); vy = vortet(2, 2); vz = vortet(2, 3)
      wx = vortet(3, 1); wy = vortet(3, 2); wz = vortet(3, 3)

      do ihgt = 1, 3
         do ihgj = 1, 3
            tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
         end do
      end do

      svort = 0.5 * (vortet + tvort)
      ovort = 0.5 * (vortet - tvort)

      snorm = sqrt(2.00 * ((svort(1, 1) * svort(1, 1)) + (svort(1, 2) * svort(1, 2)) + (svort(1, 3) * svort(1, 3)) + &
         & (svort(2, 1) * svort(2, 1)) + (svort(2, 2) * svort(2, 2)) + (svort(2, 3) * svort(2, 3)) + &
         & (svort(3, 1) * svort(3, 1)) + (svort(3, 2) * svort(3, 2)) + (svort(3, 3) * svort(3, 3))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
      onorm = sqrt(2.00 * ((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + (ovort(1, 3) * ovort(1, 3)) + &
         & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2)) + (ovort(2, 3) * ovort(2, 3)) + &
         & (ovort(3, 1) * ovort(3, 1)) + (ovort(3, 2) * ovort(3, 2)) + (ovort(3, 3) * ovort(3, 3))))
      omega = onorm

      divnorm = ux + vy + wz !careful with the sign. if it becomes very big, it can produce negative production

      usss = sqrt((2.0 * ((ux * ux) + (vy * vy) + (wz * wz))) &
         & + ((uy + vx) * (uy + vx) + (uz + wx) * (uz + wx) + (wy + vz) * (wy + vz)) &
         & - (2.0 / 3.0 * (ux + vy + wz) * (ux + vy + wz)))

      squaret = (sqrt((rec_grads(i, nof_variables, 1)**2) + (rec_grads(i, nof_variables, 2)**2) &
         & + (rec_grads(i, nof_variables, 3)**2)))**2

      leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
      rightv(1:nof_variables) = leftv(1:nof_variables)
      call molecular_transport(n, leftv, viscl(1), laml(1))
      viscl(2) = viscl(1); laml(2) = laml(1)
      viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0

      turbmv(1) = u_ct_val(i, 1, 1)
      turbmv(2) = turbmv(1)

      select case (turbulencemodel)

      case (2) ! SST frozen-coefficient destruction Jacobian, dS/d(rho*k,rho*omega)
         call sst_destruction_diagonal(n, i, 3, source_t)

      end select

   end subroutine sources_derivatives

   subroutine sources_computation2d(n)
      implicit none
!> @brief
!> sources  computation in 2d
      integer, intent(in) :: n
      integer :: i, kmaxe

      if (.not.sources_required_2d()) return
      call sa_rc_prepare(n)
      kmaxe = xmpielrank(n)
      if(axisymmetric_enabled==1) then
#ifdef xpu
!$omp target teams distribute parallel do firstprivate(n,kmaxe) private(i) &
!$omp& map(alloc: axis_radius,ielem_totvolume,ielem_walldist,u_c_val,u_ct_val,rec_grads,rhs_val,rhst_val) &
!$omp& map(alloc: gamma_in,mp_janaf,mp_m,mp_pinf,mp_thigh_in,mp_tlow_in,mp_tmid_in,rg_hzero,rg_molm) &
!$omp& map(alloc: rg_thetag,rgs_ab,rgs_bb,rgs_cb,rgs_eps_over_k,rgs_mg,rgs_sigmaa)
#else
!$omp do private(i)
#endif
         do i=1,kmaxe
            call sources_axisymmetric_cell(n,i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
      if ((turbulence .eq. 1) .and. (kmaxe .gt. 0)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe, n) &
!$omp& private(i) &
!$omp& map(alloc: ielem_totvolume, ielem_walldist, ielem_walltrans, mp_janaf, mp_m, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rec_grads, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb, rhst_val, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_computation2d_turb_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2) .and. (rungekutta .ge. 10) .and. (kmaxe .gt. 0)) then
#ifdef xpu
         !$omp target teams distribute parallel do &
!$omp& firstprivate(nof_variables, kmaxe) &
!$omp& map(alloc: sht_rg) &
!$omp& private(i)
#else
         !$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            ! Clear the implicit source Jacobian before assembling this step.
            sht_rg(i, 1:nof_variables) = 0.00
         end do
#ifdef xpu
         !$omp end target teams distribute parallel do
#else
         !$omp end do
#endif
      end if
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2) .and. (kmaxe .gt. 0)) then
#ifdef xpu
         !$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe, n) &
!$omp& private(i) &
!$omp& map(alloc: gamma_in, ielem_totvolume, mp_janaf, mp_m, mp_pinf, mp_thigh_in, mp_tlow_in, mp_tmid_in, rg_hzero) &
!$omp& map(alloc: rg_molm, rg_thetag, rhs_val, sht_rg, u_c_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_computation2d_realgas_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if

   end subroutine sources_computation2d

   subroutine sources_computation2d_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      if (axisymmetric_enabled==1) call sources_axisymmetric_cell(n,iconsidered)
      if (turbulence .eq. 1) call sources_computation2d_turb_cell(n, iconsidered)
      if ((realgas .eq. 1) .and. (rg_relax .eq. 2)) call sources_computation2d_realgas_cell(n, iconsidered)

   end subroutine sources_computation2d_cell

   subroutine sources_computation2d_turb_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      integer :: iv
      real, dimension(gpu_max_turbulence) :: source_t

      call sources2d(n, iconsidered, source_t)
      do iv = 1, turbulenceequations
         rhst_val(iconsidered, iv) = rhst_val(iconsidered, iv) - source_t(iv) * ielem_totvolume(iconsidered)
      end do

   end subroutine sources_computation2d_turb_cell

   subroutine sources_computation2d_realgas_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      integer :: iv
      real :: source_r(gpu_max_nvar), src_jac_diag(gpu_max_nvar)
! Pseudo-time changes the update, never the physical source residual.
      call sources_realgas(n, iconsidered, source_r, src_jac_diag, 0.0e0)

      do iv = 1, nof_variables
         rhs_val(iconsidered, iv) = rhs_val(iconsidered, iv) - source_r(iv) * ielem_totvolume(iconsidered)
         if (rungekutta .ge. 10) sht_rg(iconsidered, iv) = src_jac_diag(iv) * ielem_totvolume(iconsidered)
      end do

   end subroutine sources_computation2d_realgas_cell

   subroutine sources_derivatives_computation2d(n)
      implicit none
!> @brief
!> sources  derivatives computation in 2d
      integer, intent(in) :: n
      integer :: i, kmaxe

      kmaxe = xmpielrank(n)
      if ((turbulence .eq. 1) .and. (kmaxe .gt. 0)) then
#ifdef xpu
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe, n) &
!$omp& private(i) &
!$omp& map(alloc: ielem_totvolume, ielem_walldist, ielem_walltrans, mp_janaf, mp_m, mp_thigh_in, mp_tlow_in) &
!$omp& map(alloc: mp_tmid_in, rec_grads, rg_hzero, rg_molm, rgs_ab, rgs_bb, rgs_cb, sht, u_c_val, u_ct_val)
#else
!$omp do &
!$omp& private(i)
#endif
         do i = 1, kmaxe
            call sources_derivatives_computation2d_cell(n, i)
         end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
      end if
   end subroutine sources_derivatives_computation2d

   subroutine sources_derivatives_computation2d_cell(n, iconsidered)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      integer :: iv
      real, dimension(gpu_max_turbulence) :: source_t

      call sources_derivatives2d(n, iconsidered, source_t)
      do iv = 1, turbulenceequations
         sht(iconsidered, iv) = source_t(iv) * ielem_totvolume(iconsidered)
      end do

   end subroutine sources_derivatives_computation2d_cell

   subroutine sources2d(n, iconsidered, source_t)
      implicit none
!> @brief
!> sources  computation in 2d
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, dimension(turbulenceequations), intent(inout) :: source_t
      real :: alpha_star0_local
      real :: intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
      real :: vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
      real :: tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
      integer :: i, k, j, l, ihgt, ihgj, iex, lowre
      real :: snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
      real :: gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
      real :: r_des, f_des, ddw, f_des_sst, l_t_des
      real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
      real :: uz, vz, wx, wy, wz
      real :: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz !for sas only
      real, dimension(2, 2) :: vortet, tvort, svort, ovort
!declarations for k-omega
      real :: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
      real :: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
      real :: alpha_raw, alpha_star, re_t_sst, alpha_inf
      real :: beta_stari, beta_i, beta_raw, beta_star
      real :: k_0, om_0, wally
      real :: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
      real :: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2
      real :: kx, ky, kz, omx, omy, omz, rho, nu, nu_tilde, difterm, prod
      real, dimension(1:gpu_max_nvar) :: leftv, rightv
      real :: mp_pinfl, gammal
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm

      real :: rc_src, rc_jac

      i = iconsidered
      if (turbulencemodel == 1) then
         call sources_sa_rc(n, i, rc_src, rc_jac)
         source_t(1) = rc_src
         return
      end if


      verysmall = 10e-16

      vortet(1:2, 1:2) = rec_grads(i, 1:2, 1:2)

      ux = vortet(1, 1); uy = vortet(1, 2)
      vx = vortet(2, 1); vy = vortet(2, 2)

      do ihgt = 1, 2
         do ihgj = 1, 2
            tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
         end do
      end do

      svort = 0.5 * (vortet + tvort)
      ovort = 0.5 * (vortet - tvort)

      snorm = sqrt(2.00 * ((svort(1, 1) * svort(1, 1)) + (svort(1, 2) * svort(1, 2)) + &
         & (svort(2, 1) * svort(2, 1)) + (svort(2, 2) * svort(2, 2))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
      onorm = sqrt(2.00 * ((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + &
         & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2))))
      omega = onorm

      wz=0.0
      if (axisymmetric_enabled==1) wz=u_c_val(i,1,3)/(u_c_val(i,1,1)*axis_radius(i))
      if (axisymmetric_enabled==1) snorm=sqrt(snorm**2+2.0*wz**2)
      divnorm = ux + vy + wz !careful with the sign. if it becomes very big, it can produce negative production

      usss = sqrt((2.0 * ((ux * ux) + (vy * vy))) &
         & + ((uy + vx) * (uy + vx)) &
         & - (2.0 / 3.0 * (ux + vy) * (ux + vy)))

      squaret = (sqrt((rec_grads(i, nof_variables, 1)**2) + (rec_grads(i, nof_variables, 2)**2)))**2

      leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
      rightv(1:nof_variables) = leftv(1:nof_variables)

      call molecular_transport(n, leftv, viscl(1), laml(1))
      viscl(2) = viscl(1); laml(2) = laml(1)
      viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0

      turbmv(1) = u_ct_val(i, 1, 1)
      turbmv(2) = turbmv(1)

      select case (turbulencemodel)

      case (2) !k omega sst

         eddyfl(1) = ielem_walldist(i)
         eddyfl(2) = u_ct_val(i, 1, 1)
         eddyfl(3) = u_ct_val(i, 1, 2)
         call fill_turbulence_gradients(2, i, eddyfl)


         eddyfr = eddyfl

         call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
         k_0 = (max(verysmall, u_ct_val(i, 1, 1) / leftv(1))) !first subindex makes reference to the time-stepping
         om_0 = max(1.0e-6, u_ct_val(i, 1, 2) / leftv(1))
         wally = max(ielem_walldist(i),1.0e-15)

         !calculate here k and omega gradients
         omx = rec_grads(i, nof_variables + 1, 1)
         omy = rec_grads(i, nof_variables + 1, 2)

         kx = rec_grads(i, nof_variables, 1)
         ky = rec_grads(i, nof_variables, 2)

         dervk_dervom = kx * omx + ky * omy

         call sst_blending(leftv(1), viscl(1), k_0, om_0, wally, dervk_dervom, f_1, f_2)
         alpha_raw = f_1 * alpha_inf1 + (1.0-f_1) * alpha_inf2
         beta_star = beta_starinf
         beta_raw = f_1 * beta_i1 + (1.0-f_1) * beta_i2

         ! SST-2003m: use the same limited production in both equations.
         ! vort_model=1 remains the explicitly selected Kato-Launder variant.
         if (vort_model == 1) then
            prod_k = viscl(3) * snorm * onorm
         else
            prod_k = viscl(3) * snorm**2
         end if
         prod_k = min(prod_k, 10.0 * leftv(1) * beta_star * k_0 * om_0)
         prod_om = alpha_raw * prod_k * leftv(1) / max(viscl(3),1.0e-30)

         !destruction terms
         !destruction of k
         ydest_k = leftv(1) * beta_star * k_0 * om_0
         !destruction of omega
         ydest_om = leftv(1) * beta_raw * om_0 * om_0

         !crossed-diffusion term
         !crossed diffusion of omega
         diff_om = 2.0 * (1 - f_1) * leftv(1) / (om_0 * sigma_om2) * dervk_dervom

         !qsas term: scale adaptive
         if (qsas_model .eq. 1) then !<------------!!!!!!!!!!!!!!!!!!!
            !calculate here second derivative of u
            !declare all variables
            do iex = 1, 2
               vortet(iex, 1:2) = rec_grads(i, 3 + turbulenceequations + iex, 1:2)

            end do

            uxx = vortet(1, 1); uyy = vortet(1, 2);
            vxx = vortet(2, 1); vyy = vortet(2, 2);
            !compute here cell volume
            cell_volume = ielem_totvolume(i)

            u_lapl = sqrt((uxx + uyy + uzz)**2 + (vxx + vyy + vzz)**2 + (wxx + wyy + wzz)**2)
            dervk2 = kx * kx + ky * ky + kz * kz
            dervom2 = omx * omx + omy * omy + omz * omz

            delta_cell = exp((1.00 / 3.00) * log(cell_volume))

            l_sas = sqrt(k_0) / (sqrt(sqrt(beta_star)) * om_0)
            l_vk = max(kappa * snorm / u_lapl, & !this switch provides high wave-number damping
               & c_smg * delta_cell * sqrt(kappa * eta2_sas / (beta_raw / beta_star - alpha_raw)))

            q_sas1 = leftv(1) * eta2_sas * kappa * snorm**2 * (l_sas / l_vk)**2
            q_sas2 = -c_sas * 2 * leftv(1) * k_0 / sigma_phi * max(dervk2 / k_0**2, dervom2 / om_0**2)

            q_sas = max(q_sas1 + q_sas2, 0.0)
         else
            q_sas = zero
         end if

         !final source terms

         !des-sst model (if qsas_model=2)
         if (qsas_model .eq. 2) then
            cell_volume = ielem_totvolume(i)
            delta_cell = exp((1.00 / 3.00) * log(cell_volume))
            l_t_des = sqrt(k_0) / (beta_star * om_0)
            f_des_sst = max(1.0, l_t_des / (c_des_sst * delta_cell) * (1 - f_2))
            !the (1-f_2) is meant to protect the boundary layer. will result in same
            !separation point that standard s-a
            ydest_k = f_des_sst * ydest_k
         end if

         srcfull_k = prod_k - ydest_k
         srcfull_om = prod_om - ydest_om + diff_om + q_sas

         !filling the output vector
         source_t(1) = srcfull_k
         source_t(2) = srcfull_om

      end select

   end subroutine sources2d

   subroutine sources_derivatives2d(n, iconsidered, source_t)
      implicit none
!> @brief
!> sources derivatives computation in 2d
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: n, iconsidered
      real, dimension(turbulenceequations), intent(inout) :: source_t
      real :: alpha_star0_local
      real :: intenergy, r1, u1, v1, w1, et1, s1, ie1, p1, skin1, e1, rs, us, vs, ws, khx
      real :: vhx, amp, dvel, omega, squaret, tch_x, tch_x3, tch_fv1, tch_fv2
      real :: tch_rs, tch_r, tch_g, tch_glim, tch_fw, tch_dif, tch_dest, tch_prod
      integer :: i, k, j, l, ihgt, ihgj, iex, lowre
      real :: snorm, onorm, divnorm, ax, ay, az, tch_shh, tch_sav, verysmall, onesix, prodterm1, stild, rr
      real :: gg, fw, destterm, fodt, srcfull, dbpr, dbdi, dbde, dby, dbx, prodtermfinal
      real :: r_des, f_des, ddw, f_des_sst, l_t_des
      real :: ux, uy, vx, vy, shear, sratio, prodmod, cvor, stildmod, prodterm2, sfac, sss, usss, ssss, s_bar, kron
      real :: uz, vz, wx, wy, wz
      real :: uxx, uyy, uzz, vxx, vyy, vzz, wxx, wyy, wzz !for sas only
      real, dimension(2, 2) :: vortet, tvort, svort, ovort
!declarations for k-omega
      real :: srcfull_k, srcfull_om, prod_k, prod_om, ydest_k, ydest_om, diff_om, q_sas
      real :: sigma_k, sigma_om, f_1, f_2, phi_1, phi_2, d_omplus !-------------- those are for diffusion too!
      real :: alpha_raw, alpha_star, re_t_sst, alpha_inf
      real :: beta_stari, beta_i, beta_raw, beta_star
      real :: k_0, om_0, wally
      real :: dervk_dervom, dervom2, dervk2, u_lapl !generalization of the velocity laplacian
      real :: l_sas, l_vk, delta_cell, cell_volume, q_sas1, q_sas2, rho, nu, nu_tilde, difterm, prod, dstild_dnut, dest, dfw_dgg
      real :: kx, ky, kz, omx, omy, omz, dprod_dnut, ddif_dnut, ddest_dnut, dif, drr_dnut, dfw_drr, dgg_drr
      real, dimension(1:gpu_max_nvar) :: leftv, rightv
      real :: mp_pinfl, gammal
      real, dimension(1:4) :: viscl, laml
      real, dimension(1:20) :: eddyfl, eddyfr
      real, dimension(1:2) :: turbmv
      real, dimension(1) :: etvm

      real :: rc_src, rc_jac

      i = iconsidered
      if (turbulencemodel == 1) then
         call sources_sa_rc(n, i, rc_src, rc_jac)
         source_t(1) = rc_jac
         return
      end if


      verysmall = 10e-16

      vortet(1:2, 1:2) = rec_grads(i, 1:2, 1:2)

      ux = vortet(1, 1); uy = vortet(1, 2)
      vx = vortet(2, 1); vy = vortet(2, 2)

      do ihgt = 1, 2
         do ihgj = 1, 2
            tvort(ihgt, ihgj) = vortet(ihgj, ihgt)
         end do
      end do

      svort = 0.5 * (vortet + tvort)
      ovort = 0.5 * (vortet - tvort)

      snorm = sqrt(2.00 * ((svort(1, 1) * svort(1, 1)) + (svort(1, 2) * svort(1, 2)) + &
         & (svort(2, 1) * svort(2, 1)) + (svort(2, 2) * svort(2, 2))))
!quadratic mean of the strain tensor (defined as svort). also needed in sst
      onorm = sqrt(2.00 * ((ovort(1, 1) * ovort(1, 1)) + (ovort(1, 2) * ovort(1, 2)) + &
         & (ovort(2, 1) * ovort(2, 1)) + (ovort(2, 2) * ovort(2, 2))))
      omega = onorm

      wz=0.0
      if (axisymmetric_enabled==1) wz=u_c_val(i,1,3)/(u_c_val(i,1,1)*axis_radius(i))
      if (axisymmetric_enabled==1) snorm=sqrt(snorm**2+2.0*wz**2)
      divnorm = ux + vy + wz !careful with the sign. if it becomes very big, it can produce negative production

      usss = sqrt((2.0 * ((ux * ux) + (vy * vy))) &
         & + ((uy + vx) * (uy + vx)) &
         & - (2.0 / 3.0 * (ux + vy) * (ux + vy)))

      squaret = (sqrt((rec_grads(i, nof_variables, 1)**2) + (rec_grads(i, nof_variables, 2)**2)))**2

      leftv(1:nof_variables) = u_c_val(i, 1, 1:nof_variables)
      rightv(1:nof_variables) = leftv(1:nof_variables)
      call molecular_transport(n, leftv, viscl(1), laml(1))
      viscl(2) = viscl(1); laml(2) = laml(1)
      viscl(3:4) = 0.0e0; laml(3:4) = 0.0e0

      turbmv(1) = u_ct_val(i, 1, 1)
      turbmv(2) = turbmv(1)

      select case (turbulencemodel)

      case (2) ! SST frozen-coefficient destruction Jacobian, dS/d(rho*k,rho*omega)
         call sst_destruction_diagonal(n, i, 2, source_t)

      end select

   end subroutine sources_derivatives2d

! Synchronize only the source terms selected for this run. Geometry and material
! tables are unchanged after startup and already have valid host copies.


   ! Spalart--Shur SA-noft2-RC algebra, shared by CPU and OpenMP target callers.
   ! Shur et al., AIAA J. 38(5), 2000. Geometry preparation assumes a fixed mesh.
   real function sa_rc_factor(grad_u, material_strain) result(fr1)
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: grad_u(3,3), material_strain(3,3)
      real :: s(3,3), w(3,3), a(3,3), ds(3,3), scale, sn, wn, d2, rt
      integer :: i, j, k
      scale = maxval(abs(grad_u))
      fr1 = 1.0
      if (scale <= tiny(1.0)**0.25) return
      a = grad_u / scale
      s = 0.5*(a + transpose(a)); w = 0.5*(a - transpose(a))
      ! Scaling avoids fourth powers of dimensional strain/vorticity.
      ds = (material_strain / scale) / scale
      sn = sqrt(2.0*sum(s*s)); wn = sqrt(2.0*sum(w*w))
      d2 = 0.5*(sn*sn + wn*wn)
      rt = 0.0
      do k = 1, 3
         do j = 1, 3
            do i = 1, 3
               rt = rt + 2.0*w(i,k)*s(j,k)*ds(i,j)
            end do
         end do
      end do
      rt = rt/(d2*d2)
      ! cr1=1, cr2=12, cr3=1. No Fluent-style [0,1.25] clipping.
      ! S/(S+omega) is r*/(1+r*) without division by zero at omega=0.
      fr1 = 4.0*sn/(sn+wn)*(1.0-atan(12.0*rt))-1.0
   end function sa_rc_factor

   subroutine sa_noft2_rc_source(rho, mu, phi, omega, distance, grad_nut_sq, fr1, nut_floor, &
                               cb1, cb2, cv1, cw1, cw2, cw3, kap, sig, src, jac, rotation_delta)
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: rho, mu, phi, omega, distance, grad_nut_sq, fr1, nut_floor
      real, intent(in) :: cb1, cb2, cv1, cw1, cw2, cw3, kap, sig
      real, intent(out) :: src, jac
      real, intent(in), optional :: rotation_delta
      real :: rotprod
      real :: nut, nu, chi, chi3, cv13, fv1, fv2, dfv1, dfv2, den, kd2
      real :: st, dst, rr, dr, gg, dg, fw, dfw, ff, d, prod, dest, dif
      nut = max(phi/rho, nut_floor); nu = mu/rho
      d = max(distance, 1.0e-11); kd2 = (kap*d)**2
      chi = nut/nu; chi3 = chi**3; cv13 = cv1**3
      fv1 = chi3/(chi3+cv13)
      dfv1 = 3.0*chi**2*cv13/((chi3+cv13)**2*nu)
      den = 1.0+chi*fv1
      fv2 = 1.0-chi/den
      dfv2 = -(1.0/nu-chi*chi*dfv1)/(den*den)
      st = omega + fv2*nut/kd2
      dst = (fv2+nut*dfv2)/kd2
      if (st <= 0.30*omega) then
         st = 0.30*omega; dst = 0.0
      end if
      ! The destruction model uses UNCORRECTED, nonnegative stilde.
      rr = 10.0; dr = 0.0
      if (st > 0.0) then
         rr = nut/(kd2*st)
         if (rr >= 10.0) then
            rr = 10.0
         else
            dr = (1.0 - nut*dst/st)/(kd2*st)
         end if
      end if
      gg = rr + cw2*(rr**6-rr)
      dg = (1.0+cw2*(6.0*rr**5-1.0))*dr
      ff = ((1.0+cw3**6)/(gg**6+cw3**6))**(1.0/6.0)
      fw = gg*ff
      dfw = ff*cw3**6/(gg**6+cw3**6)*dg
      rotprod = 0.0
      if (present(rotation_delta)) rotprod = rotation_delta
      prod = cb1*rho*nut*(fr1*st+rotprod)
      dest = cw1*rho*fw*(nut/d)**2
      dif = cb2*rho*grad_nut_sq/sig
      src = prod + dif - dest
      ! d(src)/d(phi), phi=rho*nut, with rho, mu, gradients and fr1 frozen.
      jac = cb1*(fr1*(st+nut*dst)+rotprod) - cw1*(dfw*(nut/d)**2+2.0*fw*nut/d**2)
      if (phi/rho <= nut_floor) jac = 0.0
   end subroutine sa_noft2_rc_source

   subroutine sa_rc_inverse(a, ndim, ainv, ok)
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: a(3,3)
      integer, intent(in) :: ndim
      real, intent(out) :: ainv(3,3)
      logical, intent(out) :: ok
      real :: work(3,6), row(6), pivot
      integer :: i, j, k, p
      work = 0.0; ainv = 0.0; ok = .false.
      work(:,1:3) = a
      do i = 1, ndim
         work(i,3+i) = 1.0
      end do
      do k = 1, ndim
         p = k
         do i = k+1, ndim
            if (abs(work(i,k)) > abs(work(p,k))) p = i
         end do
         if (abs(work(p,k)) <= 100.0*epsilon(1.0)*maxval(abs(a))) return
         row = work(k,:); work(k,:) = work(p,:); work(p,:) = row
         pivot = work(k,k); work(k,:) = work(k,:)/pivot
         do j = 1, ndim
            if (j == k) cycle
            pivot = work(j,k); work(j,:) = work(j,:) - pivot*work(k,:)
         end do
      end do
      ainv(1:ndim,1:ndim) = work(1:ndim,4:3+ndim)
      ok = .true.
   end subroutine sa_rc_inverse

   real function sa_source_jacobian_for_diagonal(value, correction, model) result(value_used)
#ifdef xpu
!$omp declare target
#endif
      real, intent(in) :: value
      integer, intent(in) :: correction, model
      value_used = max(value,0.0)
      ! Retain implicit destruction for every SA variant and SST. Callers
      ! cap positive derivatives and retain their existing diagonal floor.
      if (model == 1 .or. model == 2) value_used = value
   end function sa_source_jacobian_for_diagonal

   subroutine rc_fail(message)
      character(*), intent(in) :: message
      integer :: ierr
      write(*,*) 'SA-RC error: ', trim(message)
      call mpi_abort(mpi_comm_world, 72, ierr)
   end subroutine

   integer function rc_neighbour(i, q, rank) result(j)
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i, q, rank
      integer :: group
      j = rec_ihexl(i,1,q)
      if (rec_local(i) /= 0) then
         if (rec_ihexb(rec_local(i),1,q) /= rank) then
            group = rec_ihexn(rec_local(i),1,q)
            j = rc_cells + halo_offset(group) + j - 1
         end if
      end if
   end function

   function rc_displacement(i,j) result(dx)
#ifdef xpu
!$omp declare target
#endif
      integer, intent(in) :: i,j
      real :: dx(3)
      if (j <= rc_cells) then
         dx = [ielem_xxc(j),ielem_yyc(j),0.0]
         if (dimensiona == 3) dx(3) = ielem_zzc(j)
      else
         dx = rc_xyz(j-rc_cells,:)
      end if
      dx(1) = dx(1)-ielem_xxc(i); dx(2) = dx(2)-ielem_yyc(i)
      if (dimensiona == 3) dx(3) = dx(3)-ielem_zzc(i)
   end function rc_displacement

   subroutine rc_exchange()
      ! Called by the master thread only. Buffers are packed by halo cell,
      ! six values per cell; no noncontiguous MPI temporaries or device MPI.
      integer :: g, nr, offset, count, ierr
      nr = 0
      do g = 1, ineedhalo
         count = 6*halo_len(g)
         if (count == 0) cycle
         offset = 6*(halo_offset(g)-1)+1; nr = nr+1
         call mpi_irecv(rc_recv(offset), count, mpi_double_precision, halo_proc(g), 4371, &
                        mpi_comm_world, rc_requests(nr), ierr)
      end do
      do g = 1, ineedhalos
         count = 6*halos_len(g)
         if (count == 0) cycle
         offset = 6*(halos_offset(g)-1)+1; nr = nr+1
         call mpi_isend(rc_send(offset), count, mpi_double_precision, halos_proc(g), 4371, &
                        mpi_comm_world, rc_requests(nr), ierr)
      end do
      call mpi_waitall(nr, rc_requests, mpi_statuses_ignore, ierr)
   end subroutine

   subroutine sa_rc_initialize(rank)
      integer, intent(in) :: rank
      integer :: i, j, q, a, b, bad, all_bad, ierr
      real :: dx(3), v(3), aa(3,3), ai(3,3), scale(3), radius
      logical :: ok
      if (rot_corr /= 2 .or. turbulence /= 1 .or. turbulencemodel /= 1) return
      if (allocated(sa_rc_fr1)) call rc_fail('initialization called twice')
      if (.not. allocated(halo_len)) call rc_fail('higher-order halo layout is required')
      if (halo_total > huge(0)/6 .or. halos_total > huge(0)/6) call rc_fail('local halo count overflow')
      rc_cells = xmpielrank(rank)
      allocate(rc_xyz(halo_total,3), rc_strain(rc_cells+halo_total,6))
      allocate(rc_inverse(rc_cells,3,3), rc_scale(rc_cells,3), sa_rc_fr1(rc_cells))
      allocate(rc_send(max(1,6*halos_total)), rc_recv(max(1,6*halo_total)))
      allocate(rc_requests(max(1,ineedhalo+ineedhalos)))
      rc_xyz = 0.0; rc_strain = 0.0; rc_inverse = 0.0; rc_scale = 1.0; sa_rc_fr1 = 1.0
      rc_send = 0.0; rc_recv = 0.0
      do i = 1, halos_total
         j = solhi_loc(i)
         rc_send(6*i-5:6*i-3) = [ielem_xxc(j),ielem_yyc(j),0.0]
         if (dimensiona == 3) rc_send(6*i-3) = ielem_zzc(j)
      end do
      call rc_exchange()
      do i = 1, halo_total
         rc_xyz(i,:) = rc_recv(6*i-5:6*i-3)
      end do
      ! Inverse-distance weighted point LS, independently of the MUSCL/WENO
      ! polynomial order. Column equilibration handles stretched wall cells.
      bad = 0
!$omp parallel do private(i,j,q,a,b,dx,v,aa,ai,scale,radius,ok) reduction(+:bad)
      do i = 1, rc_cells
         scale = 0.0
         do q = 2, ielem_inumneighbours(i)
            j = rc_neighbour(i,q,rank)
            dx = rc_displacement(i,j)
            radius = sqrt(sum(dx*dx))
            if (radius <= tiny(1.0)**0.25) cycle
            scale = max(scale,abs(dx/radius))
         end do
         scale = max(scale, tiny(1.0)**0.25)
         aa = 0.0
         do q = 2, ielem_inumneighbours(i)
            j = rc_neighbour(i,q,rank)
            dx = rc_displacement(i,j)
            radius = sqrt(sum(dx*dx))
            if (radius <= tiny(1.0)**0.25) cycle
            v = dx/radius/scale
            do b = 1, dimensiona
               do a = 1, dimensiona
                  aa(a,b) = aa(a,b) + v(a)*v(b)
               end do
            end do
         end do
         call sa_rc_inverse(aa, dimensiona, ai, ok)
         if (.not. ok) bad = bad+1
         rc_inverse(i,:,:) = ai; rc_scale(i,:) = scale
      end do
!$omp end parallel do
      call mpi_allreduce(bad, all_bad, 1, mpi_integer, mpi_sum, mpi_comm_world, ierr)
      if (all_bad > 0) then
         if (rank == 0) write(*,*) 'SA-RC rank-deficient central stencils: ', all_bad
         call rc_fail('strain-gradient stencil cannot resolve all spatial directions')
      end if
#ifdef xpu
!$omp target enter data map(to: rc_xyz, rc_strain, rc_inverse, rc_scale, rc_send, rc_recv, sa_rc_fr1, rc_cells)
#endif
      if (rank == 0) write(*,*) 'SA-noft2-RC: steady Spalart-Shur; cr=(1,12,1); unclipped fr1; strain LS + MPI halos'
   end subroutine

   subroutine sa_rc_prepare(rank)
      integer, intent(in) :: rank
      integer :: i, j, a, b, k, q, c
      real :: gu(3,3), strain(3,3), ds(3,3), material(6), rhs(3,6), grad(3), dx(3), v(3), vel(3), radius
      if (rot_corr /= 2 .or. turbulence /= 1 .or. turbulencemodel /= 1) return
      ! Every thread enters this routine; all MPI calls are funneled to master.
#ifdef xpu
!$omp target teams distribute parallel do private(i,gu,strain) firstprivate(dimensiona,rc_cells) &
!$omp& map(alloc: rec_grads,rc_strain)
#else
!$omp do private(i,gu,strain)
#endif
      do i = 1, rc_cells
         gu = 0.0; gu(1:dimensiona,1:dimensiona) = rec_grads(i,1:dimensiona,1:dimensiona)
         strain = 0.5*(gu+transpose(gu))
         rc_strain(i,:) = [strain(1,1),strain(2,2),strain(3,3),strain(1,2),strain(1,3),strain(2,3)]
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do private(i,c) firstprivate(halos_total) &
!$omp& map(alloc: rc_send,rc_strain,solhi_loc)
#else
!$omp end do
!$omp do private(i,c)
#endif
      do i = 1, halos_total
         do c = 1, 6
            rc_send(6*(i-1)+c) = rc_strain(solhi_loc(i),c)
         end do
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
!$omp target update from(rc_send)
#else
!$omp end do
#endif
!$omp master
      call rc_exchange()
!$omp end master
!$omp barrier
#ifdef xpu
!$omp target update to(rc_recv)
!$omp target teams distribute parallel do private(i,c) firstprivate(halo_total,rc_cells) &
!$omp& map(alloc: rc_recv,rc_strain)
#else
!$omp do private(i,c)
#endif
      do i = 1, halo_total
         do c = 1, 6
            rc_strain(rc_cells+i,c) = rc_recv(6*(i-1)+c)
         end do
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
!$omp target teams distribute parallel do &
!$omp& private(i,j,a,b,k,q,c,gu,ds,material,rhs,grad,dx,v,vel,radius) &
!$omp& firstprivate(rank,dimensiona,rc_cells) &
!$omp& map(alloc: rc_xyz,rc_strain,rc_inverse,rc_scale,sa_rc_fr1,rec_grads,u_c_val) &
!$omp& map(alloc: ielem_xxc,ielem_yyc,ielem_zzc) &
!$omp& map(alloc: rec_ihexl,rec_local,rec_ihexb,rec_ihexn,halo_offset,ielem_inumneighbours)
#else
!$omp end do
!$omp do private(i,j,a,b,k,q,c,gu,ds,material,rhs,grad,dx,v,vel,radius)
#endif
      do i = 1, rc_cells
         rhs = 0.0; vel = 0.0
         vel(1:dimensiona) = u_c_val(i,1,2:dimensiona+1)/u_c_val(i,1,1)
         do q = 2, ielem_inumneighbours(i)
            j = rc_neighbour(i,q,rank)
            dx = rc_displacement(i,j)
            radius = sqrt(sum(dx*dx))
            if (radius <= tiny(1.0)**0.25) cycle
            v = dx/radius/rc_scale(i,:)
            do c = 1, 6
               rhs(:,c) = rhs(:,c) + v*(rc_strain(j,c)-rc_strain(i,c))/radius
            end do
         end do
         do c = 1, 6
            grad = 0.0
            do a = 1, dimensiona
               do b = 1, dimensiona
                  grad(a) = grad(a) + rc_inverse(i,a,b)*rhs(b,c)
               end do
            end do
            grad = grad/rc_scale(i,:)
            material(c) = sum(vel*grad)
         end do
         ds(1,:) = [material(1),material(4),material(5)]
         ds(2,:) = [material(4),material(2),material(6)]
         ds(3,:) = [material(5),material(6),material(3)]
         gu = 0.0; gu(1:dimensiona,1:dimensiona) = rec_grads(i,1:dimensiona,1:dimensiona)
         sa_rc_fr1(i) = sa_rc_factor(gu,ds)
      end do
#ifdef xpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
   end subroutine


   subroutine sa_rc_release()
      if (.not. allocated(sa_rc_fr1)) return
#ifdef xpu
!$omp target exit data map(delete: rc_xyz,rc_strain,rc_inverse,rc_scale,rc_send,rc_recv,sa_rc_fr1,rc_cells)
#endif
      deallocate(rc_xyz,rc_strain,rc_inverse,rc_scale,rc_send,rc_recv,sa_rc_fr1,rc_requests)
      rc_cells = 0
   end subroutine
end module source
