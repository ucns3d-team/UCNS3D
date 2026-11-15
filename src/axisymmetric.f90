! Cylindrical geometry for the planar finite-volume equations, without swirl.
! Volumes in the residual remain planar areas. Do NOT radius-weight their fluxes.
module axisymmetric_support
   use declaration
   use mpi
   use mpiinfo, only: n
   use transform, only: quadratureline
   use mesh_local_nodes, only: node_coords, node_xyz, node_slot, local_nodes_enabled
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: read_axisymmetric_input, validate_axisymmetric_options, prepare_axisymmetric_geometry
   public :: snap_axisymmetric_nodes
   public :: sources_required_2d, axis_ratio, axis_stress, axisymmetric_restart_metadata
contains
   logical function sources_required_2d() result(required_source)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      required_source = (turbulence == 1) .or. (axisymmetric_enabled == 1)
      ! Keep this separate: RG_RELAX need not have been initialised for a perfect gas.
      if (realgas == 1) required_source = required_source .or. (rg_relax == 2)
   end function sources_required_2d

   subroutine axis_error(message)
      implicit none
      character(*), intent(in) :: message
      integer :: rank, ierr
      call mpi_comm_rank(mpi_comm_world, rank, ierr)
      if (rank == 0) write(*,'(A)') 'AXISSYMMETRIC.DAT: '//message
      call mpi_abort(mpi_comm_world, 1, ierr)
   end subroutine axis_error

   subroutine read_axisymmetric_input
      implicit none
      integer :: unit, ios
      logical :: exists
      character(512) :: message
      namelist /axisymmetric/ axis_y, axis_tolerance_relative
      axisymmetric_enabled = 0
      axis_y = 0.0
      axis_tolerance_relative = 1.0e-10
      axis_tolerance = 0.0
      inquire(file='AXISSYMMETRIC.DAT', exist=exists)
      if (.not. exists) return
      open(newunit=unit, file='AXISSYMMETRIC.DAT', status='old', action='read', iostat=ios, iomsg=message)
      if (ios /= 0) call axis_error(trim(message))
      read(unit, nml=axisymmetric, iostat=ios, iomsg=message)
      close(unit)
      if (ios /= 0) call axis_error('expected &axisymmetric namelist: '//trim(message))
      if (.not. ieee_is_finite(axis_y)) call axis_error('axis_y must be finite')
      if (.not. ieee_is_finite(axis_tolerance_relative)) call axis_error('axis tolerance must be finite')
      if (axis_tolerance_relative <= 0.0 .or. axis_tolerance_relative > 1.0e-4) &
         call axis_error('axis_tolerance_relative must be in (0, 1e-4]')
      axisymmetric_enabled = 1
   end subroutine read_axisymmetric_input

   subroutine validate_axisymmetric_options
      implicit none
      if (axisymmetric_enabled == 0) return
      if (dimensiona /= 2) call axis_error('requires a 2D meridional mesh, x axial and r=y-axis_y >= 0')
      if (itestcase /= 3 .and. itestcase /= 4) call axis_error('only Euler/Navier-Stokes equations are supported')
      if (dg /= 0 .or. mood /= 0) call axis_error('DG and MOOD are not supported by this finite-volume implementation')
      ! Higher-order reconstruction is allowed with the existing cell-centred sources.
      ! Their quadrature remains second order for smooth, regular solutions and
      ! can limit the spatial accuracy of the complete axisymmetric discretisation.
      if (fastest_q /= 1) call axis_error('requires stored face quadrature (fastest_q=1)')
      if (reduce_comp /= 0) call axis_error('quadrature-free fluxes are unsupported in axisymmetric mode')
      if (multispecies /= 0) call axis_error('MULTISPECIES multiphase models are unsupported; REALGAS is supported')
      if (rframe /= 0 .or. movement /= 0 .or. swirl /= 0) call axis_error('requires a stationary mesh without swirl or rotation')
      if (abs(aoa) > 1.0e-12 .or. abs(vvel) > 1.0e-12 .or. abs(wvel) > 1.0e-12) &
         call axis_error('freestream must be axial: zero incidence and zero radial/azimuthal velocity')
      if (turbulence == 1) then
         if (itestcase /= 4) call axis_error('turbulence requires Navier-Stokes equations')
         if (turbulencemodel /= 1 .and. turbulencemodel /= 2) call axis_error('only SA and SST RANS are supported')
         if (des_model /= 0 .or. qsas_model /= 0 .or. ihybrid /= 0) call axis_error('DES/DDES/SAS/hybrid modes require 3D')
         if (rot_corr /= 0) call axis_error('rotation/curvature turbulence corrections are not supported in axisymmetric mode')
      end if
      if (n == 0) write(*,'(A,ES16.8)') 'Axisymmetric finite volumes (no swirl), axis y = ', axis_y
   end subroutine validate_axisymmetric_options

   subroutine snap_axisymmetric_nodes
      implicit none
      ! Remove coordinate roundoff at r=0 BEFORE constructing any cell metrics.
      ! This changes only the in-memory mesh; the user's mesh file is never edited.
      integer :: i,j,id,slot,ierr,nc,count,global_count
      real :: bounds(4),global_bounds(4),coord(2),length,shift,maxshift,global_shift
      if(axisymmetric_enabled==0) return
      nc=xmpielrank(n)
      bounds=-huge(1.0)
      do i=1,nc
         do j=1,ielem_nonodes(i)
            coord=node_coords(ielem_nodes(j,i))
            bounds=max(bounds,[coord(1),-coord(1),coord(2),-coord(2)])
         end do
      end do
      call mpi_allreduce(bounds,global_bounds,4,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
      length=max(global_bounds(1)+global_bounds(2),global_bounds(3)+global_bounds(4))
      if(.not.ieee_is_finite(length).or.length<=0.0) call axis_error('invalid mesh extent')
      axis_tolerance=max(axis_tolerance_relative*length, &
         64.0*epsilon(1.0)*max(length,abs(axis_y),maxval(abs(global_bounds))))
      count=0;maxshift=0.0
      do i=1,nc
         do j=1,ielem_nonodes(i)
            id=ielem_nodes(j,i);coord=node_coords(id);shift=coord(2)-axis_y
            if(abs(shift)<=axis_tolerance.and.shift/=0.0) then
               if(local_nodes_enabled) then
                  slot=node_slot(id);node_xyz(2,slot)=axis_y
               else
                  dinoder(id)%cord(2)=axis_y
               end if
               count=count+1;maxshift=max(maxshift,abs(shift))
            end if
         end do
      end do
      call mpi_allreduce(count,global_count,1,mpi_integer,mpi_sum,mpi_comm_world,ierr)
      call mpi_allreduce(maxshift,global_shift,1,mpi_double_precision,mpi_max,mpi_comm_world,ierr)
      if(n==0.and.global_count>0) write(*,'(A,I0,A,ES12.4)') &
         'Axis coordinate roundoff removed from ',global_count,' local vertex copies; maximum shift = ',global_shift
   end subroutine snap_axisymmetric_nodes

   subroutine prepare_axisymmetric_geometry
      implicit none
      ! Called before vertex coordinates are released, outside OpenMP parallel regions.
      integer :: i, f, j, ng, ierr, bad, bad_global, nc, bc
      real :: bounds(4), global_bounds(4), coord(2), a(2), b(2), length, centre(2), det
      real :: roundoff
      real :: vertices(8,2), points(2,numberofpoints2), weights(numberofpoints2)
      if (axisymmetric_enabled == 0) return
      nc = xmpielrank(n)
      bounds = [-huge(1.0), -huge(1.0), -huge(1.0), -huge(1.0)]
      do i=1,nc
         do j=1,ielem_nonodes(i)
            coord = node_coords(ielem_nodes(j,i))
            bounds = max(bounds, [coord(1), -coord(1), coord(2), -coord(2)])
         end do
      end do
      call mpi_allreduce(bounds, global_bounds, 4, mpi_double_precision, mpi_max, mpi_comm_world, ierr)
      length = max(global_bounds(1)+global_bounds(2), global_bounds(3)+global_bounds(4))
      if (.not. ieee_is_finite(length) .or. length <= 0.0) call axis_error('invalid mesh extent')
      axis_tolerance = axis_tolerance_relative*length
      roundoff = 64.0*epsilon(1.0)*max(length,abs(axis_y),maxval(abs(global_bounds)))
      axis_tolerance = max(axis_tolerance,roundoff)
      allocate(axis_radius(nc),axis_face_radius(nc,max_faces,numberofpoints2),axis_face(nc,max_faces))
      axis_radius = ielem_yyc(1:nc)-axis_y
      axis_face_radius = 0.0; axis_face = .false.
      bad = 0
      do i=1,nc
         if (.not. ieee_is_finite(axis_radius(i)) .or. axis_radius(i) <= roundoff) bad=1
         do j=1,ielem_nonodes(i)
            coord = node_coords(ielem_nodes(j,i))
            if (.not. all(ieee_is_finite(coord))) bad=1
            ! Tolerance identifies boundaries; it does not legalise negative-radius cells.
            if (coord(2)-axis_y < -roundoff) bad=1
         end do
         do f=1,ielem_ifca(i)
            a = node_coords(ielem_nodes_faces(i,f,1))
            b = node_coords(ielem_nodes_faces(i,f,2))
            axis_face(i,f) = max(abs(a(2)-axis_y),abs(b(2)-axis_y)) <= axis_tolerance
            bc=0
            if(ielem_ibounds(i,f)>0) bc=ibound_icode(ielem_ibounds(i,f))
            if (axis_face(i,f)) then
               if(bc/=3) bad=1
               if(max(abs(a(2)-axis_y),abs(b(2)-axis_y))>roundoff) bad=1
            end if
            vertices=0.0; vertices(1,:)=a; vertices(2,:)=b
            call quadratureline(n,igqrules,vertices,points,weights)
            axis_face_radius(i,f,1:qp_line)=max(points(2,1:qp_line)-axis_y,0.0)
            if(axis_face(i,f)) axis_face_radius(i,f,:)=0.0
         end do
      end do
      call mpi_allreduce(bad,bad_global,1,mpi_integer,mpi_max,mpi_comm_world,ierr)
      if(bad_global/=0) call axis_error('mesh must have positive-radius cell centres, nonnegative nodes, and BC 3 on exact axis faces')
   end subroutine prepare_axisymmetric_geometry

   real function axis_ratio(radial_velocity, radius, radial_derivative) result(ratio)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real,intent(in) :: radial_velocity, radius, radial_derivative
      ! Use a regular limit only AT the axis. Small positive radii are never clipped.
      if(radius>0.0) then
         ratio=radial_velocity/radius
      else
         ratio=radial_derivative
      end if
   end function axis_ratio

   subroutine axis_stress(mu, ux, ur, vx, vr, hoop, txx, txr, trr, ttt)
      implicit none
#ifdef xpu
!$omp declare target
#endif
      real,intent(in) :: mu,ux,ur,vx,vr,hoop
      real,intent(out) :: txx,txr,trr,ttt
      real :: trace
      trace=(2.0/3.0)*mu*(ux+vr+hoop)
      txx=2.0*mu*ux-trace; txr=mu*(ur+vx)
      trr=2.0*mu*vr-trace; ttt=2.0*mu*hoop-trace
   end subroutine axis_stress

   subroutine axisymmetric_restart_metadata(reading)
      implicit none
      logical,intent(in) :: reading
      logical :: exists
      integer :: unit,ios,mode,ierr
      real :: saved_axis
      if (reading) then
         inquire(file='AXISSYMMETRIC.restart',exist=exists)
         if (.not.exists) then
            if(axisymmetric_enabled==1) call axis_error('restart lacks AXISSYMMETRIC.restart geometry metadata; start fresh')
            return
         end if
         open(newunit=unit,file='AXISSYMMETRIC.restart',status='old',action='read',iostat=ios)
         if(ios/=0) call axis_error('cannot read restart geometry metadata')
         read(unit,*,iostat=ios) mode,saved_axis
         close(unit)
         if(ios/=0) call axis_error('invalid restart geometry metadata')
         if(mode/=axisymmetric_enabled) call axis_error('restart planar/axisymmetric mode mismatch')
         if(mode==1) then
            if(.not.ieee_is_finite(saved_axis)) call axis_error('invalid restart axis')
            if(saved_axis/=axis_y) call axis_error('restart axis location mismatch')
         end if
      else
         if(n/=0) return
         ! Preserve legacy planar directories unless geometry metadata already exists.
         inquire(file='AXISSYMMETRIC.restart',exist=exists)
         if(axisymmetric_enabled==0.and..not.exists) return
         open(newunit=unit,file='AXISSYMMETRIC.restart',status='replace',action='write',iostat=ios)
         if(ios/=0) call axis_error('cannot write restart geometry metadata')
         write(unit,'(I1,1X,ES25.17)',iostat=ios) axisymmetric_enabled,axis_y
         close(unit)
         if(ios/=0) call axis_error('failed writing restart geometry metadata')
      end if
   end subroutine axisymmetric_restart_metadata
end module axisymmetric_support
