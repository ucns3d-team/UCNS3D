MODULE MESHQULITY_module
    USE LIBRARY
    USE MPIINFO
    USE OMP_LIB
    USE DECLARATION
    ! USE FLUXES
    ! USE TRANSFORM
    ! USE TRANSFORM_MovingMesh
    ! USE PRESTORE_MOVINGMESH
    ! USE COMMUNICATIONS
    ! USE RECON

    IMPLICIT NONE

CONTAINS

function NodeMeshQuality(p_minus, p, p_plus, result)
    implicit none
    real,dimension(1:2),intent(in)::p_minus, p, p_plus
    real,intent(out)::result
    logical::NodeMeshQuality

    if (dimensiona.eq.2) then
        select case (relaxation_centre_type)
        case (5)
            NodeMeshQuality = JacobiCondNumber2D(p_minus, p, p_plus, result)
        case (7)
            NodeMeshQuality = OddyMetric2D(p_minus, p, p_plus, result)
        case DEFAULT
            ! print*,"invalid node relaxation type in NodeMeshQuality"
            NodeMeshQuality = OddyMetric2D(p_minus, p, p_plus, result)
        end select
    else
        print*,"invalid node relaxation type in NodeMeshQuality"
    end if 
end function NodeMeshQuality





function JacobiCondNumber2D(p_minus, p, p_plus, result)
    implicit none
    real,dimension(1:2),intent(in)::p_minus, p, p_plus
    real,intent(out)::result
    logical::JacobiCondNumber2D

    real::area, len_minus2, len_plus2
    integer::i

    len_minus2 = zero
    len_plus2 = zero

    do i=1, 2
        len_minus2 = len_minus2 + ((p_minus(i) - p(i))**2)
        len_plus2  = len_plus2  + ((p_plus(i)  - p(i))**2)
    end do

    area = ((p(1)-p_plus(1))*(p(2)-p_minus(2))) - ((p(2)-p_plus(2))*(p(1)-p_minus(1)))
    if (area.le.zero) Then
        JacobiCondNumber2D = .false.
    else
        JacobiCondNumber2D = .true.
    end if

    result = (len_minus2 + len_plus2)/(2.0*area)

end function





function OddyMetric2D(p_minus, p, p_plus, result)
    implicit none
    real,dimension(1:2),intent(in)::p_minus, p, p_plus
    real,intent(out)::result
    logical::OddyMetric2D

    real::area, len_minus2, len_plus2
    integer::i

    len_minus2 = zero
    len_plus2 = zero

    do i=1, 2
        len_minus2 = len_minus2 + ((p_minus(i) - p(i))**2)
        len_plus2  = len_plus2  + ((p_plus(i)  - p(i))**2)
    end do

    area = ((p(1)-p_plus(1))*(p(2)-p_minus(2))) - ((p(2)-p_plus(2))*(p(1)-p_minus(1)))
    if (area.le.zero) Then
        OddyMetric2D = .false.
    else
        OddyMetric2D = .true.
    end if

    result = 0.5*(((len_minus2 + len_plus2)**2)/(area**2) - 2.0)

end function





! function JacobiCondNumberAverage(num, values)
!     implicit none
!     integer::num
!     real,dimension(1:num),intent(in)::values
!     real::JacobiCondNumberAverage

!     integer::p
!     integer::i
!     p = 2

!     JacobiCondNumberAverage = zero

!     do i=1, num
!         JacobiCondNumberAverage = JacobiCondNumberAverage + (values(i)**p)
!     end do

!     JacobiCondNumberAverage = JacobiCondNumberAverage / real(num)

!     if (p.eq.2) then
!         JacobiCondNumberAverage = sqrt(JacobiCondNumberAverage)
!     else if (p.gt.2) then
!         JacobiCondNumberAverage = JacobiCondNumberAverage**(1.0/real(p))
!     end if

! end function





! subroutine find_node_JacobiCondNumber(moved, position_index, d_t, N)
!     implicit none
!     integer,intent(in)::moved, position_index, N
!     real,intent(in)::d_t

!     ! integer M
!     integer::node_index, node_plus_index, node_minus_index, cell_index, index, cell_num_nodes
!     integer::cpu_index, cpu
!     integer::iter, i, j, k, i_plus, i_minus
!     real,dimension(1:dimensiona)::p, p_minus, p_plus
!     real::val
!     integer::power

!     integer,dimension(2*isize)::requests
!     integer::num_requests, count

!     ! M = omp_get_thread_num()
!     if (num_values_to_send_per_node.lt.(2*dimensiona)) then
!         print *,"something went wrong sorry :("
!         call abort
!     end if

!     power = 2

!     !$omp do
!     do iter = 1,my_num_interface_nodes 
!         node_index = local_interface_nodes(iter)
!         ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
!         do cpu_index = 1,local_nodes(node_index)%num_cpus
!             cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
!             index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (2*dimensiona)
!             do j = 1, local_nodes(node_index)%num_local_neighbours
!                 cell_index = local_nodes(node_index)%local_neighbours(j)
!                 cell_num_nodes = ielem(N, cell_index)%nonodes
!                 do i = 1, cell_num_nodes
!                     if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
!                         exit
!                     end if
!                 end do
!                 i_plus = i+1
!                 if (i_plus.gt.cell_num_nodes) then
!                     i_plus = 1
!                 end if
!                 i_minus = i-1
!                 if (i_minus.lt.1) then
!                     i_minus = cell_num_nodes
!                 end if
!                 node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
!                 node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
!                 p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
!                 p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
!                 if (moved.ne.0) then
!                     p_plus(:)  = p_plus   + (local_nodes(node_plus_index)%velocity(:)*d_t)
!                     p_minus(:) = p_minus  + (local_nodes(node_minus_index)%velocity(:)*d_t)
!                 end if

!                 do k = 1,dimensiona
!                     index = index +1
!                     if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
!                         print *, "copying too much data to send buffer from", N, "to", cpu 
!                     end if
!                     node_snd_buffer(cpu)%data(index) = p_minus(k)
!                 end do
!                 do k = 1,dimensiona
!                     index = index +1
!                     if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
!                         print *, "copying too much data to send buffer from", N, "to", cpu 
!                     end if
!                     node_snd_buffer(cpu)%data(index) = p_plus(k)
!                 end do
!             end do
!         end do
!     end do
!     !$omp end do

!     !$omp barrier

!     num_requests = 0
!     !$omp master
!         ! print *, "inside send_rcv part on CPU", N
!         do cpu_index = 0, isize-1
!             ! if (cpu_index.ne.N) then
!                 if (node_snd_count(cpu_index).gt.0) then
!                     count = node_snd_count(cpu_index) * (2*dimensiona)
!                     num_requests = num_requests + 1
!                     CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
!                     ! print *, "sending from", N, "to", cpu_index, count, "values"
!                 end if
!                 if (node_rcv_count(cpu_index).gt.0) then
!                     count = node_rcv_count(cpu_index) * (2*dimensiona)
!                     num_requests = num_requests + 1
!                     CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
!                     ! print *, N, "waiting to receive", count, "values from", cpu_index
!                 end if
!             ! else
!             !     if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
!             !         print *,"send receive count missmatch on CPU", n
!             !         call abort
!             !     end if
!             !     count = node_snd_count(cpu_index) * (2*dimensiona)
!             !     do i = 1, count
!             !         node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
!             !     end do
!             ! end if
!         end do

!         ! print*,"CPU", n, "witing on", num_requests, "requests"
!         CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

!     !$omp end master

!     !$omp do
!     do node_index = 1, kmaxn 
!         local_nodes(node_index)%JacobiCondNumber = zero

!         do iter = 1, local_nodes(node_index)%num_local_neighbours
!             cell_index = local_nodes(node_index)%local_neighbours(iter)
!             cell_num_nodes = ielem(N, cell_index)%nonodes
!             do i = 1, cell_num_nodes
!                 if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
!                     exit
!                 end if
!             end do
!             i_plus = i+1
!             if (i_plus.gt.cell_num_nodes) then
!                 i_plus = 1
!             end if
!             i_minus = i-1
!             if (i_minus.lt.1) thenO2D(p_minus, p, 
!                 i_minus = cell_num_nodes
!             end if
!             node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
!             node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
!             p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
!             p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
!             p(:)       = local_nodes(node_index)%positions(position_index,:)
!             if (moved.ne.0) then
!                 p_plus(:)  = p_plus(:)  + (local_nodes(node_plus_index)%velocity(:)*d_t)
!                 p_minus(:) = p_minus(:) + (local_nodes(node_minus_index)%velocity(:)*d_t)
!                 p(:)       = p(:)       + (local_nodes(node_index)%velocity(:)*d_t)
!             end if

!             if (dimensiona.eq.2) then
!                 val = JacobiCondNumber2D(p_minus(:), p(:), p_plus(:))
!             else
!                 print*,"not implemented yet"
!                 call abort
!             end if

!             local_nodes(node_index)%JacobiCondNumber = local_nodes(node_index)%JacobiCondNumber + (val**power)
!         end do

!         do iter = 1, local_nodes(node_index)%num_cpus
!             cpu_index = local_nodes(node_index)%rcv_offsets(iter)%cpu
!             index = (local_nodes(node_index)%rcv_offsets(iter)%lower - 1)*(2*dimensiona)
!             do j = local_nodes(node_index)%rcv_offsets(iter)%lower, local_nodes(node_index)%rcv_offsets(iter)%upper
!                 if (dimensiona.eq.2) then
!                     p_minus(:) = node_rcv_buffer(cpu_index)%data(index+1              : index+dimensiona)
!                     p_plus(:)  = node_rcv_buffer(cpu_index)%data(index+(dimensiona+1) : index+(2*dimensiona))
!                     p(:)       = local_nodes(node_index)%positions(position_index,:)
!                     if (moved.ne.0) then
!                         p(:) = p(:) + (local_nodes(node_index)%velocity(:)*d_t)
!                     end if

!                     val = JacobiCondNumber2D(p_minus(:), p(:), p_plus(:))
!                 else
!                     print*,"not implemented yet"
!                     call abort
!                 end if
!                 local_nodes(node_index)%JacobiCondNumber = local_nodes(node_index)%JacobiCondNumber + (val**power)
                
!                 index = index + (2*dimensiona)
!             end do
!         end do

!         local_nodes(node_index)%JacobiCondNumber = local_nodes(node_index)%JacobiCondNumber / real(local_nodes(node_index)%num_neighbours)
!         if (power.eq.2) then
!             local_nodes(node_index)%JacobiCondNumber = sqrt(local_nodes(node_index)%JacobiCondNumber)
!         else if (power.gt.2) Then
!             local_nodes(node_index)%JacobiCondNumber = local_nodes(node_index)%JacobiCondNumber**(1.0/real(power))
!         end if

!     end do
!     !$omp end do
    
!     !$omp barrier

!     !$omp master
!         CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!     !$omp end master

!     !$omp barrier

! end subroutine find_node_JacobiCondNumber




function compute_relaxation_gradient(num_point_pairs, points, p, result)
    implicit none
    integer,intent(in)::num_point_pairs
    real,dimension(1:dimensiona,1:2*num_point_pairs),intent(in)::points
    real,dimension(1:dimensiona),intent(in)::p
    real,dimension(1:dimensiona),intent(out)::result
    logical::compute_relaxation_gradient

    if (relaxation_centre_type.eq.5) then
        compute_relaxation_gradient = compute_JacobiCondNumber_gradient(num_point_pairs, points, p, result)
    else if (relaxation_centre_type.eq.7) then
        compute_relaxation_gradient = compute_OddyMetric_gradient(num_point_pairs, points, p, result)
    else
        print*,"invalid relaxation_centre_type"
    end if

end function compute_relaxation_gradient





function compute_JacobiCondNumber_gradient(num_point_pairs, points, p, result)
    implicit none
    integer,intent(in)::num_point_pairs
    real,dimension(1:dimensiona,1:2*num_point_pairs),intent(in)::points
    real,dimension(1:dimensiona),intent(in)::p
    real,dimension(1:dimensiona),intent(out)::result
    logical::compute_JacobiCondNumber_gradient

    real::area, len2_plus, len2_minus
    integer::pair, minus_index, plus_index, i, other_i
    real,dimension(1:dimensiona)::p_plus, p_minus

    result(:) = zero
    compute_JacobiCondNumber_gradient = .true.

    if (dimensiona.eq.2) then
        do pair = 1, num_point_pairs
            plus_index = 2*pair
            minus_index = plus_index-1 
            p_plus(:) = points(:,plus_index)
            p_minus(:) = points(:,minus_index)
            len2_plus = zero
            len2_minus  = zero

            area = ((p(1)-p_plus(1))*(p(2)-p_minus(2))) - ((p(2)-p_plus(2))*(p(1)-p_minus(1)))
            if (area.lt.zero) then
                ! print*,"negative area in compute_JacobiCondNumber_gradient"
                compute_JacobiCondNumber_gradient = .false.
            else if (area.eq.zero) then
                print*,"area = 0 in compute_JacobiCondNumber_gradient"
                compute_JacobiCondNumber_gradient = .false.
            end if
            do i = 1, dimensiona
                len2_minus = len2_minus + ((p(i)-p_minus(i))**2)
                len2_plus = len2_plus + ((p(i)-p_plus(i))**2)
            end do

            result(1) = result(1) + (((2.0*p(1)) - p_minus(1) - p_plus(1))/area) + (((len2_plus+len2_minus)/(2.0*(area**2)))*(p_minus(2) - p_plus(2)))
            result(2) = result(2) + (((2.0*p(2)) - p_minus(2) - p_plus(2))/area) + (((len2_plus+len2_minus)/(2.0*(area**2)))*(p_plus(1) - p_minus(1)))
        end do
    else
        print*,"3D is not supported in compute_JacobiCondNumber_gradient in meshquality module"
    end if

end function compute_JacobiCondNumber_gradient





function compute_OddyMetric_gradient(num_point_pairs, points, p, result)
    implicit none
    integer,intent(in)::num_point_pairs
    real,dimension(1:dimensiona,1:2*num_point_pairs),intent(in)::points
    real,dimension(1:dimensiona),intent(in)::p
    real,dimension(1:dimensiona),intent(out)::result
    logical::compute_OddyMetric_gradient

    real::area, len2_plus, len2_minus
    integer::pair, minus_index, plus_index, i, other_i
    real,dimension(1:dimensiona)::p_plus, p_minus

    result(:) = zero
    compute_OddyMetric_gradient = .true.

    if (dimensiona.eq.2) then
        do pair = 1, num_point_pairs
            plus_index = 2*pair
            minus_index = plus_index-1 
            p_plus(:) = points(:,plus_index)
            p_minus(:) = points(:,minus_index)
            len2_plus = zero
            len2_minus  = zero

            area = ((p(1)-p_plus(1))*(p(2)-p_minus(2))) - ((p(2)-p_plus(2))*(p(1)-p_minus(1)))
            if (area.lt.zero) then
                ! print*,"negative area in compute_OddyMetric_gradient"
                compute_OddyMetric_gradient = .false.
            else if (area.eq.zero) then
                print*,"area = 0 in compute_OddyMetric_gradient"
                compute_OddyMetric_gradient = .false.
            end if
            do i = 1, dimensiona
                len2_minus = len2_minus + ((p(i)-p_minus(i))**2)
                len2_plus = len2_plus + ((p(i)-p_plus(i))**2)
            end do

            result(1) = result(1) + ((2.0*((2.0*p(1)) - p_minus(1) - p_plus(1))*(len2_plus+len2_minus))/(area**2)) + ((((len2_plus+len2_minus)**2)/(area**3))*(p_minus(2) - p_plus(2)))
            result(2) = result(2) + ((2.0*((2.0*p(2)) - p_minus(2) - p_plus(2))*(len2_plus+len2_minus))/(area**2)) + ((((len2_plus+len2_minus)**2)/(area**3))*(p_plus(1) - p_minus(1)))
        end do
    else
        print*,"3D is not supported in compute_OddyMetric_gradient in meshquality module"
    end if

end function compute_OddyMetric_gradient




function compute_relaxation_hessian(num_point_pairs, points, p, result)
    implicit none
    integer,intent(in)::num_point_pairs
    real,dimension(1:dimensiona,1:2*num_point_pairs),intent(in)::points
    real,dimension(1:dimensiona),intent(in)::p
    real,dimension(1:dimensiona,1:dimensiona),intent(out)::result
    logical::compute_relaxation_hessian

    if (relaxation_centre_type.eq.5) then
        compute_relaxation_hessian = compute_JacobiCondNumber_hessian(num_point_pairs, points, p, result)
    else if (relaxation_centre_type.eq.7) then
        compute_relaxation_hessian = compute_OddyMetric_hessian(num_point_pairs, points, p, result)
    else
        print*,"invalid relaxation_centre_type"
    end if

end function compute_relaxation_hessian





function compute_JacobiCondNumber_hessian(num_point_pairs, points, p, result)
    implicit none
    integer,intent(in)::num_point_pairs
    real,dimension(1:dimensiona,1:2*num_point_pairs),intent(in)::points
    real,dimension(1:dimensiona),intent(in)::p
    real,dimension(1:dimensiona,1:dimensiona),intent(out)::result
    logical::compute_JacobiCondNumber_hessian

    real::area, len2_plus, len2_minus
    integer::pair, minus_index, plus_index, i, other_i
    real,dimension(1:dimensiona)::p_plus, p_minus

    result(:,:) = zero
    compute_JacobiCondNumber_hessian = .true.

    if (dimensiona.eq.2) then
        do pair = 1, num_point_pairs
            plus_index = 2*pair
            minus_index = plus_index-1 
            p_plus(:) = points(:,plus_index)
            p_minus(:) = points(:,minus_index)
            len2_plus = zero
            len2_minus  = zero

            area = ((p(1)-p_plus(1))*(p(2)-p_minus(2))) - ((p(2)-p_plus(2))*(p(1)-p_minus(1)))
            if (area.lt.zero) then
                ! print*,"negative area in compute_JacobiCondNumber_hessian"
                compute_JacobiCondNumber_hessian = .false.
            else if (area.eq.zero) then
                print*,"area = 0 in compute_JacobiCondNumber_hessian"
                compute_JacobiCondNumber_hessian = .false.
            end if
            do i = 1, dimensiona
                len2_minus = len2_minus + ((p(i)-p_minus(i))**2)
                len2_plus = len2_plus + ((p(i)-p_plus(i))**2)
            end do

            result(1,1) = result(1,1) + (2.0/Area)
            result(1,1) = result(1,1) + (2.0*((2.0*p(1)) - p_minus(1) - p_plus(1))*(p_minus(2) - p_plus(2))/(area**2))
            result(1,1) = result(1,1) + ((len2_plus+len2_minus)*((p_minus(2) - p_plus(2))**2)/(area**3))
            result(2,2) = result(2,2) + (2.0/Area)
            result(2,2) = result(2,2) + (2.0*((2.0*p(2)) - p_minus(2) - p_plus(2))*(p_plus(1) - p_minus(1))/(area**2))
            result(2,2) = result(2,2) + ((len2_plus+len2_minus)*((p_plus(1) - p_minus(1))**2)/(area**3))
            
            result(1,2) = result(1,2) + (((2.0*p(1)) - p_minus(1) - p_plus(1))*(p_plus(1) - p_minus(1))/(area**2))
            result(1,2) = result(1,2) + ((2.0*p(2) - p_minus(2) - p_plus(2))*(p_minus(2) - p_plus(2))/(area**2))
            result(1,2) = result(1,2) + ((len2_plus+len2_minus)*(p_plus(1) - p_minus(1))*(p_minus(2) - p_plus(2)))
        end do
        result(2,1) = result(1,2)
    else
        print*,"3D is not supported in compute_JacobiCondNumber_gradient in meshquality module"
    end if

end function compute_JacobiCondNumber_hessian





function compute_OddyMetric_hessian(num_point_pairs, points, p, result)
    implicit none
    integer,intent(in)::num_point_pairs
    real,dimension(1:dimensiona,1:2*num_point_pairs),intent(in)::points
    real,dimension(1:dimensiona),intent(in)::p
    real,dimension(1:dimensiona,1:dimensiona),intent(out)::result
    logical::compute_OddyMetric_hessian

    real::area, len2_plus, len2_minus
    integer::pair, minus_index, plus_index, i, other_i
    real,dimension(1:dimensiona)::p_plus, p_minus

    result(:,:) = zero
    compute_OddyMetric_hessian = .true.

    if (dimensiona.eq.2) then
        do pair = 1, num_point_pairs
            plus_index = 2*pair
            minus_index = plus_index-1 
            p_plus(:) = points(:,plus_index)
            p_minus(:) = points(:,minus_index)
            len2_plus = zero
            len2_minus  = zero

            area = ((p(1)-p_plus(1))*(p(2)-p_minus(2))) - ((p(2)-p_plus(2))*(p(1)-p_minus(1)))
            if (area.lt.zero) then
                ! print*,"negative area in compute_OddyMetric_hessian"
                compute_OddyMetric_hessian = .false.
            else if (area.eq.zero) then
                print*,"area = 0 in compute_OddyMetric_hessian"
                compute_OddyMetric_hessian = .false.
            end if
            do i = 1, dimensiona
                len2_minus = len2_minus + ((p(i)-p_minus(i))**2)
                len2_plus = len2_plus + ((p(i)-p_plus(i))**2)
            end do

            result(1,1) = result(1,1) + ((4.0*(((2.0*p(1)) - p_minus(1) - p_plus(1))**2))/(area**2))
            result(1,1) = result(1,1) + ((4.0*(len2_plus+len2_minus))/(area**2))
            result(1,1) = result(1,1) - ((4.0*(len2_plus+len2_minus)*((2.0*p(1))-p_minus(1)-p_plus(1))*(p_plus(2)-p_minus(2)))/(area**3))
            result(1,1) = result(1,1) + ((3.0*((len2_plus+len2_minus)**2)*((p_plus(2)-p_minus(2))**2))/(area**4))

            result(2,2) = result(2,2) + ((4.0*(((2.0*p(2)) - p_minus(2) - p_plus(2))**2))/(area**2))
            result(2,2) = result(2,2) + ((4.0*(len2_plus+len2_minus))/(area**2))
            result(2,2) = result(2,2) - ((4.0*(len2_plus+len2_minus)*((2.0*p(2))-p_minus(2)-p_plus(2))*(p_minus(1)-p_plus(1)))/(area**3))
            result(2,2) = result(2,2) + ((3.0*((len2_plus+len2_minus)**2)*((p_minus(1)-p_plus(1))**2))/(area**4))
            
            result(1,2) = result(1,2) + ((4.0*((2.0*p(1)) - p_minus(1) - p_plus(1))*((2.0*p(2)) - p_minus(2) - p_plus(2)))/(area**2))
            result(1,2) = result(1,2) - ((2.0*(len2_plus+len2_minus)*((2.0*p(1))-p_minus(1)-p_plus(1))*(p_minus(1)-p_plus(1)))/(area**3))
            result(1,2) = result(1,2) - ((2.0*(len2_plus+len2_minus)*((2.0*p(2))-p_minus(2)-p_plus(2))*(p_plus(2)-p_minus(2)))/(area**3))
            result(1,2) = result(1,2) + ((3.0*((len2_plus+len2_minus)**2)*(p_minus(1)-p_plus(1))*(p_plus(2)-p_minus(2)))/(area**4))
        end do
        result(2,1) = result(1,2)
    else
        print*,"3D is not supported in compute_OddyMetric_gradient in meshquality module"
    end if

end function compute_OddyMetric_hessian





function invert_mat(size, input, output)
    implicit none
    integer,intent(in)::size
    real,dimension(1:size,1:size),intent(in)::input
    real,dimension(1:size,1:size),intent(out)::output
    logical::invert_mat

    real::det

    invert_mat = .true.

    if (size.eq.1) then
        if (input(1,1).eq.zero) then
            invert_mat = .false.
        end if
        output(1,1) = 1.0/input(1,1)
    else if (size.eq.2) then
        det = (input(1,1)*input(2,2)) - (input(1,2)*input(2,1))
        if (det.eq.zero) then
            invert_mat = .false.
            print*,"dividing by zero in invert_mat"
        else if (det.lt.zero) Then
            print*,"negative determinant in invert_mat"
        end if
        output(1,1) = input(2,2)/det
        output(2,2) = input(1,1)/det
        output(1,2) = (-1.0)*input(1,2)/det
        output(2,1) = (-1.0)*input(2,1)/det
    else
        print*,"unsupported size in invert_mat in meshquality module"
    end if

end function invert_mat





subroutine find_node_Jacobi_relaxation_velocity(moved, position_index, d_t, N)
    implicit none
    integer,intent(in)::moved, position_index, N
    real,intent(in)::d_t

    ! integer M
    integer::node_index, node_plus_index, node_minus_index, cell_index, index, node_num_neighbours, cell_num_nodes
    integer::cpu_index, cpu
    integer::iter, i, j, k, i_plus, i_minus, counter, pair
    real,dimension(1:dimensiona)::p, p_minus, p_plus
    real::val, change_len2, treshold2, multiple, direction_len, first_derivative, second_derivative, change_dot, helper 
    real,dimension(1:dimensiona,1:(2*max_num_node_neighbours))::points
    real,dimension(1:dimensiona)::relaxed_point, new_relaxed, change, relaxation_gradient, direction
    real,dimension(1:dimensiona, 1:dimensiona)::relaxation_hessian, relaxation_hessian_inverse
    logical::valid, valid_helper, valid1, valid2, valid3, done

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()
    if (num_values_to_send_per_node.lt.(2*dimensiona)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    ! treshold = 0.0001
    treshold2 = ((xmax(n)-xmin(n))*(ymax(n)-ymin(n))/imaxe)*0.0001
    ! print*,"treshold2 =", treshold2

    !$omp do
    do iter = 1,my_num_interface_nodes 
        node_index = local_interface_nodes(iter)
        ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
        do cpu_index = 1,local_nodes(node_index)%num_cpus
            cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
            index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (2*dimensiona)
            do j = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(j)
                cell_num_nodes = ielem(N, cell_index)%nonodes
                do i = 1, cell_num_nodes
                    if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                        exit
                    end if
                end do
                i_plus = i+1
                if (i_plus.gt.cell_num_nodes) then
                    i_plus = 1
                end if
                i_minus = i-1
                if (i_minus.lt.1) then
                    i_minus = cell_num_nodes
                end if
                node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
                node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
                p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
                p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
                if (moved.ne.0) then
                    p_plus(1:dimensiona)  = p_plus(1:dimensiona)   + (local_nodes(node_plus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                    p_minus(1:dimensiona) = p_minus(1:dimensiona)  + (local_nodes(node_minus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                end if

                do k = 1, dimensiona
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_Jacobi_relaxation_velocity"
                    end if
                    node_snd_buffer(cpu)%data(index) = p_minus(k)
                end do
                do k = 1, dimensiona
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_Jacobi_relaxation_velocity"
                    end if
                    node_snd_buffer(cpu)%data(index) = p_plus(k)
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp barrier

    num_requests = 0
    !$omp master
        ! print *, "inside send_rcv part on CPU", N
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                if (node_snd_count(cpu_index).gt.0) then
                    count = node_snd_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 19+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 19+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * (2*dimensiona)
                do i = 1, count
                    node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
                end do
            end if
        end do

        ! print*,"CPU", n, "witing on", num_requests, "requests"
        CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

    !$omp end master

    !$omp barrier

    !$omp do
    do node_index = 1, kmaxn 
        local_nodes(node_index)%mesh_quality_before = zero
        node_num_neighbours = local_nodes(node_index)%num_neighbours

        local_nodes(node_index)%relaxation_velocity(:) = zero

        p(:) = local_nodes(node_index)%positions(position_index,:)
        if (moved.ne.0) then
            p(1:dimensiona) = p(1:dimensiona) + (local_nodes(node_index)%lagrangian_velocity(1:dimensiona)*d_t)
        end if

        counter = 0
        valid = .true.
        do iter = 1, local_nodes(node_index)%num_local_neighbours
            cell_index = local_nodes(node_index)%local_neighbours(iter)
            cell_num_nodes = ielem(N, cell_index)%nonodes
            do i = 1, cell_num_nodes
                if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                    exit
                end if
            end do
            i_plus = i+1
            if (i_plus.gt.cell_num_nodes) then
                i_plus = 1
            end if
            i_minus = i-1
            if (i_minus.lt.1) then
                i_minus = cell_num_nodes
            end if
            node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
            node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
            p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
            p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
            if (moved.ne.0) then
                p_plus(1:dimensiona)  = p_plus(1:dimensiona)  + (local_nodes(node_plus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                p_minus(1:dimensiona) = p_minus(1:dimensiona) + (local_nodes(node_minus_index)%lagrangian_velocity(1:dimensiona)*d_t)
            end if

            valid = valid.and.NodeMeshQuality(p_minus(:), p(:), p_plus(:), helper)
            if (valid) then
                local_nodes(node_index)%mesh_quality_before = local_nodes(node_index)%mesh_quality_before + helper
            else
                local_nodes(node_index)%mesh_quality_before = -1.0*(abs(local_nodes(node_index)%mesh_quality_before) + abs(helper))
            end if
            ! if (dimensiona.eq.2) then
            !     val = JacobiCondNumber2D(p_minus(:), p(:), p_plus(:), valid_helper)
            !     valid = valid.and.valid_helper
            ! else
            !     print*,"not implemented yet"
            !     call abort
            ! end if

            counter = counter+1
            points(:,counter) = p_minus(:)
            counter = counter+1
            points(:,counter) = p_plus(:)
        end do

        do iter = 1, local_nodes(node_index)%num_cpus
            cpu_index = local_nodes(node_index)%rcv_offsets(iter)%cpu
            index = (local_nodes(node_index)%rcv_offsets(iter)%lower - 1)*(2*dimensiona)
            do j = local_nodes(node_index)%rcv_offsets(iter)%lower, local_nodes(node_index)%rcv_offsets(iter)%upper
                if (dimensiona.eq.2) then
                    p_minus(:) = node_rcv_buffer(cpu_index)%data(index+1              : index+dimensiona)
                    p_plus(:)  = node_rcv_buffer(cpu_index)%data(index+(dimensiona+1) : index+(2*dimensiona))
                else
                    print*,"not implemented yet"
                    call abort
                end if
                valid = valid.and.NodeMeshQuality(p_minus(:), p(:), p_plus(:), helper)
                if (valid) then
                    local_nodes(node_index)%mesh_quality_before = local_nodes(node_index)%mesh_quality_before + helper
                else
                    local_nodes(node_index)%mesh_quality_before = -1.0*(abs(local_nodes(node_index)%mesh_quality_before) + abs(helper))
                end if

                counter = counter+1
                points(:,counter) = p_minus(:)
                counter = counter+1
                points(:,counter) = p_plus(:)

                index = index + (2*dimensiona)
            end do
        end do

        if (counter.ne.(2*node_num_neighbours)) Then
            print*,"wrong number of points in find_node_Jacobi_relaxation_velocity"
        end if

        !if (valid) then
        !     relaxed_point(:) = p(:)
        ! else
        !     print*,"alternative starting point in node", node_index
            relaxed_point(:) = zero
            do i = 1, 2*node_num_neighbours
                relaxed_point(:) = relaxed_point(:) + (points(:,i)/real(2*node_num_neighbours))
            end do
        ! end if
        change_len2 = 10000.0

        ! change_len2 = 1000.0
        ! ! if (node_num_neighbours.gt.1) then
        ! if (local_nodes(node_index)%boundary.eq.0) then  
        !     counter = 0
        !     do while ((change_len2.gt.(treshold*treshold)).and.(counter.lt.3))
        !     ! do while (change_len2.gt.(treshold*treshold))
        !         call compute_relaxation_gradient(node_num_neighbours, points(:,:), relaxed_point(:), relaxation_gradient(:))
        !         call compute_relaxation_hessian(node_num_neighbours, points(:,:), relaxed_point(:), relaxation_hessian(:,:))
        !         call invert_mat(dimensiona, relaxation_hessian(:,:), relaxation_hessian_inverse(:,:))

        !         change(:) = zero
        !         change_len2 = zero
        !         do i = 1, dimensiona
        !             do j = 1, dimensiona
        !                 change(i) = change(i) + (0.5*(relaxation_hessian_inverse(i,j)*relaxation_gradient(j)))
        !             end do
        !             change_len2 = change_len2 + (change(i)**2)
        !         end do

        !         relaxed_point(:) = relaxed_point(:) - change(:)
        !         counter = counter+1
        !     end do
        ! end if
        
        if (local_nodes(node_index)%boundary.eq.0) then  
            valid1 = compute_relaxation_gradient(node_num_neighbours, points(:,:), relaxed_point(:), relaxation_gradient(:))
            valid2 = compute_relaxation_hessian(node_num_neighbours, points(:,:), relaxed_point(:), relaxation_hessian(:,:))
            valid3 = invert_mat(dimensiona, relaxation_hessian(:,:), relaxation_hessian_inverse(:,:))
            valid = valid1.and.valid2.and.valid3

            if (.not.valid) then
                print*,"something went wrong and the initial guess in subroutine find_node_Jacobi_relaxation_velocity did not provide valid gradient and hessian inverse"
            end if

            counter = 0
            do while ((change_len2.gt.treshold2).and.(valid).and.(counter.lt.20))
                multiple = 1.0

                107 continue
                change_len2 = zero
                change(:) = zero
                do i = 1, dimensiona
                    do j = 1, dimensiona
                        change(i) = change(i) + (multiple*(relaxation_hessian_inverse(i,j)*relaxation_gradient(j)))
                    end do
                    change_len2 = change_len2 + (change(i)**2)
                end do

                new_relaxed(:) = relaxed_point(:) - change(:)

                valid1 = compute_relaxation_gradient(node_num_neighbours, points(:,:), new_relaxed(:), relaxation_gradient(:))
                valid2 = compute_relaxation_hessian(node_num_neighbours, points(:,:), new_relaxed(:), relaxation_hessian(:,:))
                valid3 = invert_mat(dimensiona, relaxation_hessian(:,:), relaxation_hessian_inverse(:,:))
                valid = valid1.and.valid2.and.valid3
                if (valid) then
                    relaxed_point(:) = new_relaxed(:)
                    counter = counter+1
                    ! print*,"internal node", node_index, "counter", counter
                else
                    if (multiple.gt.0.06125) then
                        multiple = 0.5*multiple
                        print*,"decrementing multiple internal cell"
                        goto 107
                    else
                        print*,"internal cell fail"
                        done = .true.
                    end if
                end if
            end do

        else if (node_num_neighbours.gt.1) then ! not a corrner cell
            direction(:) = zero
            do i = 1, node_num_neighbours
                direction(:) = direction(:) + points(:,2*i) - points(:,(2*i)-1)
            end do
            direction_len = zero
            do i = 1, dimensiona
                direction_len = direction_len + (direction(i)**2)
            end do
            if (direction_len.le.0.0000000001) then
                print*,"something went wrong and the direction vector in subroutine find_node_Jacobi_relaxation_velocity has zero length"
            end if
            direction_len = sqrt(direction_len)
            direction(:) = direction(:)/direction_len
            ! print*,"p(:) =", p(1), p(2), "direction(:) =", direction(1), direction(2)

            ! if (.not.valid) then
                change(:) = relaxed_point(:) - p(:)
                change_dot = zero
                do i = 1, dimensiona
                    change_dot = change_dot + (change(i)*direction(i))
                end do
                change(:) = change_dot * direction(:)
                relaxed_point(:) = p(:) + change(:)
            ! end if

            valid1 = compute_relaxation_gradient(node_num_neighbours, points(:,:), relaxed_point(:), relaxation_gradient(:))
            valid2 = compute_relaxation_hessian(node_num_neighbours, points(:,:), relaxed_point(:), relaxation_hessian(:,:))
            first_derivative = zero
            second_derivative = zero
            do i = 1, dimensiona
                first_derivative = first_derivative + (relaxation_gradient(i)*direction(i))
            end do
            do i = 1, dimensiona
                do j = 1, dimensiona
                    second_derivative = second_derivative + (direction(i)*relaxation_hessian(i,j)*direction(j))
                end do    
            end do
            valid = valid1.and.valid2.and.(second_derivative.ne.zero)

            if (.not.valid) then
                print*,"something went wrong and the initial guess in subroutine find_node_Jacobi_relaxation_velocity did not provide valid gradient and hessian inverse"
            end if

            counter = 0
            do while ((change_len2.gt.treshold2).and.(valid).and.(counter.lt.20))
                multiple = 1.0

                207 continue
                change(:) = (multiple*(first_derivative/second_derivative))*direction(:)
                change_len2 = (multiple*(first_derivative/second_derivative))**2

                new_relaxed(:) = relaxed_point(:) - change(:)

                valid1 = compute_relaxation_gradient(node_num_neighbours, points(:,:), new_relaxed(:), relaxation_gradient(:))
                valid2 = compute_relaxation_hessian(node_num_neighbours, points(:,:), new_relaxed(:), relaxation_hessian(:,:))
                first_derivative = zero
                second_derivative = zero
                do i = 1, dimensiona
                    first_derivative = first_derivative + (relaxation_gradient(i)*direction(i))
                end do
                do i = 1, dimensiona
                    do j = 1, dimensiona
                        second_derivative = second_derivative + (direction(i)*relaxation_hessian(i,j)*direction(j))
                    end do    
                end do
                valid = valid1.and.valid2.and.(second_derivative.ne.zero)
                if (valid) then
                    relaxed_point(:) = new_relaxed(:)
                    counter = counter+1
                    ! print*,"boundary node", node_index, "counter", counter
                else
                    if (multiple.gt.0.06125) then
                        print*,"decrementing multiple boundary cell"
                        multiple = 0.5*multiple
                        goto 207
                    else
                        print*,"boundary cell fail"
                        done = .true.
                    end if
                end if
            end do
        end if
        if (relaxed_point(1).lt.xmin(n)) then
            relaxed_point(1) = 0.5*(p(1) + xmin(n))
        else if (relaxed_point(1).gt.xmax(n)) then
            relaxed_point(1) = 0.5*(p(1) + xmax(n))
        end if
        if (relaxed_point(2).lt.ymin(n)) then
            relaxed_point(2) = 0.5*(p(2) + ymin(n))
        else if (relaxed_point(2).gt.ymax(n)) then
            relaxed_point(2) = 0.5*(p(2) + ymax(n))
        end if
        local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (relaxed_point(1:dimensiona) - p(1:dimensiona))/d_t

        local_nodes(node_index)%mesh_quality_after = zero
        valid = .true.
        do pair = 1, node_num_neighbours
            i_plus = 2*pair
            i_minus = i_plus-1 
            p_plus(:) = points(:,i_plus)
            p_minus(:) = points(:,i_minus)
            
            valid = valid.or.NodeMeshQuality(p_minus(:), relaxed_point(:), p_plus(:), helper)
            if (valid) then
                local_nodes(node_index)%mesh_quality_after = local_nodes(node_index)%mesh_quality_after + helper
            else
                local_nodes(node_index)%mesh_quality_after = -1.0*(abs(local_nodes(node_index)%mesh_quality_after) + abs(helper))
            end if
        end do
    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

end subroutine find_node_Jacobi_relaxation_velocity





subroutine find_mesh_quality_before_relaxation(moved, position_index, d_t, N)
    implicit none
    integer,intent(in)::moved, position_index, N
    real,intent(in)::d_t

    ! integer M
    integer::node_index, node_plus_index, node_minus_index, cell_index, index, node_num_neighbours, cell_num_nodes
    integer::cpu_index, cpu
    integer::iter, i, j, k, i_plus, i_minus, counter
    real::helper
    logical::valid
    real,dimension(1:dimensiona)::p, p_minus, p_plus
    real,dimension(1:dimensiona,1:(2*max_num_node_neighbours))::points
    
    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()
    if (num_values_to_send_per_node.lt.(2*dimensiona)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
    do iter = 1,my_num_interface_nodes 
        node_index = local_interface_nodes(iter)
        ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
        do cpu_index = 1,local_nodes(node_index)%num_cpus
            cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
            index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (2*dimensiona)
            do j = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(j)
                cell_num_nodes = ielem(N, cell_index)%nonodes
                do i = 1, cell_num_nodes
                    if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                        exit
                    end if
                end do
                i_plus = i+1
                if (i_plus.gt.cell_num_nodes) then
                    i_plus = 1
                end if
                i_minus = i-1
                if (i_minus.lt.1) then
                    i_minus = cell_num_nodes
                end if
                node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
                node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
                p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
                p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
                if (moved.ne.0) then
                    p_plus(1:dimensiona)  = p_plus(1:dimensiona)   + (local_nodes(node_plus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                    p_minus(1:dimensiona) = p_minus(1:dimensiona)  + (local_nodes(node_minus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                end if

                do k = 1, dimensiona
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_mesh_quality_before_relaxation"
                    end if
                    node_snd_buffer(cpu)%data(index) = p_minus(k)
                end do
                do k = 1, dimensiona
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_mesh_quality_before_relaxation"
                    end if
                    node_snd_buffer(cpu)%data(index) = p_plus(k)
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp barrier

    num_requests = 0
    !$omp master
        ! print *, "inside send_rcv part on CPU", N
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                if (node_snd_count(cpu_index).gt.0) then
                    count = node_snd_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 19+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 19+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * (2*dimensiona)
                do i = 1, count
                    node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
                end do
            end if
        end do

        ! print*,"CPU", n, "witing on", num_requests, "requests"
        CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

    !$omp end master

    !$omp barrier

    !$omp do
    do node_index = 1, kmaxn 
        local_nodes(node_index)%mesh_quality_before = zero
        node_num_neighbours = local_nodes(node_index)%num_neighbours

        p(:) = local_nodes(node_index)%positions(position_index,:)
        if (moved.ne.0) then
            p(1:dimensiona) = p(1:dimensiona) + (local_nodes(node_index)%lagrangian_velocity(1:dimensiona)*d_t)
        end if

        counter = 0
        valid = .true.
        do iter = 1, local_nodes(node_index)%num_local_neighbours
            cell_index = local_nodes(node_index)%local_neighbours(iter)
            cell_num_nodes = ielem(N, cell_index)%nonodes
            do i = 1, cell_num_nodes
                if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                    exit
                end if
            end do
            i_plus = i+1
            if (i_plus.gt.cell_num_nodes) then
                i_plus = 1
            end if
            i_minus = i-1
            if (i_minus.lt.1) then
                i_minus = cell_num_nodes
            end if
            node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
            node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
            p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
            p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
            if (moved.ne.0) then
                p_plus(1:dimensiona)  = p_plus(1:dimensiona)  + (local_nodes(node_plus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                p_minus(1:dimensiona) = p_minus(1:dimensiona) + (local_nodes(node_minus_index)%lagrangian_velocity(1:dimensiona)*d_t)
            end if

            valid = valid.and.NodeMeshQuality(p_minus(:), p(:), p_plus(:), helper)
            if (valid) then
                local_nodes(node_index)%mesh_quality_before = local_nodes(node_index)%mesh_quality_before + helper
            else
                local_nodes(node_index)%mesh_quality_before = -1.0*(abs(local_nodes(node_index)%mesh_quality_before) + abs(helper))
            end if
            ! if (dimensiona.eq.2) then
            !     val = JacobiCondNumber2D(p_minus(:), p(:), p_plus(:), valid_helper)
            !     valid = valid.and.valid_helper
            ! else
            !     print*,"not implemented yet"
            !     call abort
            ! end if

            counter = counter+1
            points(:,counter) = p_minus(:)
            counter = counter+1
            points(:,counter) = p_plus(:)
        end do

        do iter = 1, local_nodes(node_index)%num_cpus
            cpu_index = local_nodes(node_index)%rcv_offsets(iter)%cpu
            index = (local_nodes(node_index)%rcv_offsets(iter)%lower - 1)*(2*dimensiona)
            do j = local_nodes(node_index)%rcv_offsets(iter)%lower, local_nodes(node_index)%rcv_offsets(iter)%upper
                if (dimensiona.eq.2) then
                    p_minus(:) = node_rcv_buffer(cpu_index)%data(index+1              : index+dimensiona)
                    p_plus(:)  = node_rcv_buffer(cpu_index)%data(index+(dimensiona+1) : index+(2*dimensiona))
                else
                    print*,"not implemented yet"
                    call abort
                end if
                valid = valid.and.NodeMeshQuality(p_minus(:), p(:), p_plus(:), helper)
                if (valid) then
                    local_nodes(node_index)%mesh_quality_before = local_nodes(node_index)%mesh_quality_before + helper
                else
                    local_nodes(node_index)%mesh_quality_before = -1.0*(abs(local_nodes(node_index)%mesh_quality_before) + abs(helper))
                end if

                counter = counter+1
                points(:,counter) = p_minus(:)
                counter = counter+1
                points(:,counter) = p_plus(:)

                index = index + (2*dimensiona)
            end do
        end do

        if (counter.ne.(2*node_num_neighbours)) Then
            print*,"wrong number of points in find_mesh_quality_before_relaxation"
        end if
    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

end subroutine find_mesh_quality_before_relaxation





subroutine find_mesh_quality_after_relaxation(moved, position_index, d_t, relaxation_timescale, N)
    implicit none
    integer,intent(in)::moved, position_index, N
    real,intent(in)::d_t, relaxation_timescale

    ! integer M
    integer::node_index, node_plus_index, node_minus_index, cell_index, index, node_num_neighbours, cell_num_nodes
    integer::cpu_index, cpu
    integer::iter, i, j, k, i_plus, i_minus, counter
    real::helper
    logical::valid
    real,dimension(1:dimensiona)::p, p_minus, p_plus
    real,dimension(1:dimensiona,1:(2*max_num_node_neighbours))::points
    
    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()
    if (num_values_to_send_per_node.lt.(2*dimensiona)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
    do iter = 1,my_num_interface_nodes 
        node_index = local_interface_nodes(iter)
        ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
        do cpu_index = 1,local_nodes(node_index)%num_cpus
            cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
            index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (2*dimensiona)
            do j = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(j)
                cell_num_nodes = ielem(N, cell_index)%nonodes
                do i = 1, cell_num_nodes
                    if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                        exit
                    end if
                end do
                i_plus = i+1
                if (i_plus.gt.cell_num_nodes) then
                    i_plus = 1
                end if
                i_minus = i-1
                if (i_minus.lt.1) then
                    i_minus = cell_num_nodes
                end if
                node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
                node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
                p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
                p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
                if (moved.ne.0) then
                    p_plus(1:dimensiona)  = p_plus(1:dimensiona)   + (local_nodes(node_plus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                    p_minus(1:dimensiona) = p_minus(1:dimensiona)  + (local_nodes(node_minus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                end if
                p_plus(1:dimensiona)  = p_plus(1:dimensiona)   + (local_nodes(node_plus_index)%relaxation_velocity(1:dimensiona)*relaxation_timescale)
                p_minus(1:dimensiona) = p_minus(1:dimensiona)  + (local_nodes(node_minus_index)%relaxation_velocity(1:dimensiona)*relaxation_timescale)

                do k = 1, dimensiona
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_mesh_quality_before_relaxation"
                    end if
                    node_snd_buffer(cpu)%data(index) = p_minus(k)
                end do
                do k = 1, dimensiona
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_mesh_quality_before_relaxation"
                    end if
                    node_snd_buffer(cpu)%data(index) = p_plus(k)
                end do
            end do
        end do
    end do
    !$omp end do

    !$omp barrier

    num_requests = 0
    !$omp master
        ! print *, "inside send_rcv part on CPU", N
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                if (node_snd_count(cpu_index).gt.0) then
                    count = node_snd_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 19+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 19+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * (2*dimensiona)
                do i = 1, count
                    node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
                end do
            end if
        end do

        ! print*,"CPU", n, "witing on", num_requests, "requests"
        CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

    !$omp end master

    !$omp barrier

    !$omp do
    do node_index = 1, kmaxn 
        local_nodes(node_index)%mesh_quality_after = zero
        node_num_neighbours = local_nodes(node_index)%num_neighbours

        p(:) = local_nodes(node_index)%positions(position_index,:)
        if (moved.ne.0) then
            p(1:dimensiona) = p(1:dimensiona) + (local_nodes(node_index)%lagrangian_velocity(1:dimensiona)*d_t)
        end if
        p(1:dimensiona) = p(1:dimensiona) + (local_nodes(node_index)%relaxation_velocity(1:dimensiona)*relaxation_timescale) 

        counter = 0
        valid = .true.
        do iter = 1, local_nodes(node_index)%num_local_neighbours
            cell_index = local_nodes(node_index)%local_neighbours(iter)
            cell_num_nodes = ielem(N, cell_index)%nonodes
            do i = 1, cell_num_nodes
                if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                    exit
                end if
            end do
            i_plus = i+1
            if (i_plus.gt.cell_num_nodes) then
                i_plus = 1
            end if
            i_minus = i-1
            if (i_minus.lt.1) then
                i_minus = cell_num_nodes
            end if
            node_plus_index  = ielem(N, cell_index)%nodes_counterclockwise(i_plus)
            node_minus_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
            p_plus(:)  = local_nodes(node_plus_index )%positions(position_index, :)
            p_minus(:) = local_nodes(node_minus_index)%positions(position_index, :)
            if (moved.ne.0) then
                p_plus(1:dimensiona)  = p_plus(1:dimensiona)  + (local_nodes(node_plus_index)%lagrangian_velocity(1:dimensiona)*d_t)
                p_minus(1:dimensiona) = p_minus(1:dimensiona) + (local_nodes(node_minus_index)%lagrangian_velocity(1:dimensiona)*d_t)
            end if

            valid = valid.and.NodeMeshQuality(p_minus(:), p(:), p_plus(:), helper)
            if (valid) then
                local_nodes(node_index)%mesh_quality_after = local_nodes(node_index)%mesh_quality_after + helper
            else
                local_nodes(node_index)%mesh_quality_after = -1.0*(abs(local_nodes(node_index)%mesh_quality_after) + abs(helper))
            end if

            counter = counter+1
            points(:,counter) = p_minus(:)
            counter = counter+1
            points(:,counter) = p_plus(:)
        end do

        do iter = 1, local_nodes(node_index)%num_cpus
            cpu_index = local_nodes(node_index)%rcv_offsets(iter)%cpu
            index = (local_nodes(node_index)%rcv_offsets(iter)%lower - 1)*(2*dimensiona)
            do j = local_nodes(node_index)%rcv_offsets(iter)%lower, local_nodes(node_index)%rcv_offsets(iter)%upper
                if (dimensiona.eq.2) then
                    p_minus(:) = node_rcv_buffer(cpu_index)%data(index+1              : index+dimensiona)
                    p_plus(:)  = node_rcv_buffer(cpu_index)%data(index+(dimensiona+1) : index+(2*dimensiona))
                else
                    print*,"not implemented yet"
                    call abort
                end if
                valid = valid.and.NodeMeshQuality(p_minus(:), p(:), p_plus(:), helper)
                if (valid) then
                    local_nodes(node_index)%mesh_quality_after = local_nodes(node_index)%mesh_quality_after + helper
                else
                    local_nodes(node_index)%mesh_quality_after = -1.0*(abs(local_nodes(node_index)%mesh_quality_after) + abs(helper))
                end if

                counter = counter+1
                points(:,counter) = p_minus(:)
                counter = counter+1
                points(:,counter) = p_plus(:)

                index = index + (2*dimensiona)
            end do
        end do

        if (counter.ne.(2*node_num_neighbours)) Then
            print*,"wrong number of points in find_mesh_quality_before_relaxation"
        end if
    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

end subroutine find_mesh_quality_after_relaxation





end module MESHQULITY_module

