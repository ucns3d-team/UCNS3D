MODULE MESHMOVEMENT_module
    USE LIBRARY
    USE MPIINFO
    USE OMP_LIB
    USE DECLARATION
    USE FLUXES
    USE TRANSFORM
    USE TRANSFORM_MovingMesh
    IMPLICIT NONE

CONTAINS





function min_abs(v1, v2) ! function returns the shorter of two vectors
    implicit none
    real,dimension(1:dimensiona)::v1, v2
    real,dimension(1:dimensiona)::min_abs
    real::sum1,sum2
    integer::i

    sum1 = zero
    sum2 = zero

    do i = 1,dimensiona
        sum1 = sum1 + (v1(i)*v1(i))
        sum2 = sum2 + (v2(i)*v2(i))
    end do
    if (sum1.lt.sum2) then
        min_abs = v1
    else
        min_abs = v2
    end if
end function





function trinagle_area(a, b, c)
    implicit none
    real,dimension(1:dimensiona)::a, b, c
    real,dimension(1:3)::v1, v2, cross
    real::trinagle_area
    integer::i

    v1 = zero
    v2 = Zero

    do i = 1, dimensiona
        v1(i) = a(i) - b(i)
        v2(i) = a(i) - c(i)
    end do

    cross(1) = (v1(2)*v2(3)) - (v1(3)*v2(2))
    cross(2) = (v1(3)*v2(1)) - (v1(1)*v2(3))
    cross(3) = (v1(1)*v2(2)) - (v1(2)*v2(1))

    trinagle_area = (cross(1)*cross(1)) + (cross(2)*cross(2)) + (cross(3)*cross(3))

    trinagle_area = sqrt(trinagle_area / 4.0)

end function trinagle_area





subroutine polygon_centre(vert_num, vertices, centre)
    ! this function assumes that the polygon is convex
    implicit none
    integer::vert_num
    real,intent(in),dimension(1:vert_num,1:dimensiona)::vertices
    real,intent(out),dimension(1:dimensiona)::centre
    real,dimension(1:dimensiona)::temp_centre, helper, v1, v2
    real::area, area_sum
    integer::i, j

    if (vert_num.le.0) then
        print *, "trying to find a centre of 0-gon"
        call abort()
    end if
    
    temp_centre = zero
    do i = 1,vert_num
        temp_centre(:) = temp_centre(:) + vertices(i,:)
    end do
    do i = 1,dimensiona
        temp_centre(i) = temp_centre(i) / real(vert_num)
    end do

    if (vert_num.le.3) then
        centre(1:dimensiona) = temp_centre(1:dimensiona)
    else
        centre = zero
        area_sum = zero
        do i = 1, vert_num
            j = i+1
            if (i.eq.vert_num) then
                j = 1
            end if

            v1 = vertices(i,:)
            v2 = vertices(j,:)

            helper(:) = ((v1(:) + v2(:) + temp_centre(:)) / 3.0)
            area = trinagle_area(v1(:), v2(:), temp_centre(:))

            area_sum = area_sum + area
            centre(:) = centre(:) + (helper(:) * area)
        end do

        centre(:) = centre(:) / area_sum
    end if

end subroutine polygon_centre





subroutine MOVE_NODES(time_step, index_from, index_to)
    implicit NONE
    integer,intent(in)::index_from,index_to
    real,intent(in)::time_step
    integer::node_index

    !$omp do
    do node_index=1,kmaxn 
        local_nodes(node_index)%positions(index_to,1:dimensiona) = local_nodes(node_index)%positions(index_from,1:dimensiona) &
                                                                 + time_step * local_nodes(node_index)%velocity(1:dimensiona)
    end do
    !$omp end do

end subroutine





subroutine COPY_BACK_LOCAL_NODES(index_from, index_to)
    implicit NONE
    integer,intent(in)::index_from,index_to
    integer::node_index

    !$omp do
    do node_index=1,kmaxn 
        local_nodes(node_index)%positions(index_to,1:dimensiona) = local_nodes(node_index)%positions(index_from,1:dimensiona)
    end do
    !$omp end do

end subroutine





function point_distance(point1, point2)
    implicit NONE
    real,dimension(dimensiona),intent(inout)::point1, point2
    real::point_distance
    real::x_coord_diff, y_coord_diff, z_coord_diff
    real::x1,x2,y1,y2,z1,z2
    real::distance2

    x_coord_diff = abs(point1(1) - point2(1))
    y_coord_diff = abs(point1(2) - point2(2))
    z_coord_diff = 0.0
    if (dimensiona.eq.3) then
        z_coord_diff = abs(point1(3) - point2(3))
    end if
    if (my_xper.gt.0.0) then
        x_coord_diff = modulo(x_coord_diff, my_xper)
        x_coord_diff = min(x_coord_diff, abs(my_xper - x_coord_diff))
    end if
    if (my_yper.gt.0.0) then
        y_coord_diff = modulo(y_coord_diff, my_yper)
        y_coord_diff = min(y_coord_diff, abs(my_yper - y_coord_diff))
    end if
    if ((dimensiona.eq.3).and.(my_zper.gt.0.0)) then
        z_coord_diff = modulo(z_coord_diff, my_zper)
        z_coord_diff = min(z_coord_diff, abs(my_zper - z_coord_diff))
    end if
    distance2 = x_coord_diff*x_coord_diff
    distance2 = distance2 + y_coord_diff*y_coord_diff
    if (dimensiona.eq.3) then 
        distance2 = distance2 + z_coord_diff*z_coord_diff
    end if
    point_distance = sqrt(distance2)

end function point_distance





subroutine establish_node_neighbours(N)
    implicit NONE
    integer,intent(in)::n
    integer::i, j, k, l, ii, iter
    integer::cell_index, node_index, node_index_2, edge_index, cpu_index, neighbour_cpu
    integer::num_requests, counter
    integer::kmaxe
    ! integer,dimension(isize)::num_nodes_to_send
    integer,dimension(0:isize-1)::num_interface_nodes
    integer,dimension(2*isize)::requests
    real::distance, distance2, coord_diff, x_coord_diff, y_coord_diff, z_coord_diff
    real::tolerance
    logical::close_enough, duplicate
    integer::helper_size, to_copy
    integer,allocatable,dimension(:)::helper

    type exchange_type
        real,allocatable,dimension(:)::coords
        integer,allocatable,dimension(:)::count
    end type

    type(exchange_type),dimension(0:isize-1)::buff

    tolerance = 0.000000001
    KMAXE=XMPIELRANK(N)

    my_num_interface_nodes = 0
    helper_size = 2
    max_num_node_neighbours = 2
    !$omp master
        allocate(helper(helper_size))
        allocate(node_rcv_count(0:isize-1))
        allocate(node_rcv_buffer(0:isize-1))
        allocate(node_snd_count(0:isize-1))
        allocate(node_snd_buffer(0:isize-1))

        do node_index=1,kmaxn 
            local_nodes(node_index)%num_local_neighbours = 0
            local_nodes(node_index)%num_neighbours       = 0
            local_nodes(node_index)%boundary             = 0 
            local_nodes(node_index)%num_cpus             = 0
        end do
    
        do cell_index = 1,kmaxe
            do i = 1,IELEM(N,cell_index)%nonodes
                node_index = IELEM(N,cell_index)%nodes(i)
                if (node_index.lt.1) then
                    print *,"invalid node index in call", cell_index, "CPU", N
                end if
                local_nodes(node_index)%Num_Local_Neighbours = local_nodes(node_index)%Num_Local_Neighbours + 1
            
                ! if (IELEM(N,cell_index)%INEIGHB(i).ne.N) then
                !     neighbour_cpu = IELEM(N,cell_index)%INEIGHB(i)
                !     ! num_nodes_to_send(neighbour_cpu) = num_nodes_to_send(neighbour_cpu) + 2
                ! end if
            end do
        end do

        do node_index=1,kmaxn 
            if (local_nodes(node_index)%Num_Local_Neighbours.eq.0) then
                print *, "0 local neighbours in node", Node_index, local_nodes(node_index)%positions(1,1), local_nodes(node_index)%positions(1,2), "CPU", N
            end if
            allocate(local_nodes(node_index)%Local_neighbours(local_nodes(node_index)%Num_Local_Neighbours))
            max_num_node_neighbours = max(max_num_node_neighbours, local_nodes(node_index)%Num_Local_Neighbours)
        end do

        do cell_index = 1,kmaxe
            do i = 1,IELEM(N,cell_index)%nonodes
                node_index = IELEM(N,cell_index)%nodes(i)
                local_nodes(node_index)%Num_Neighbours = local_nodes(node_index)%Num_Neighbours + 1
                local_nodes(node_index)%Local_neighbours(local_nodes(node_index)%Num_Neighbours) = cell_index
            end do
        end do

        do node_index=1,kmaxn 
            if (local_nodes(node_index)%Num_Neighbours.ne.local_nodes(node_index)%Num_Local_Neighbours) then
                print *, "missmatch in number of neighbours in node", node_index, "on cpu", n
            end if
        end do

        DO II=1,NOF_BOUNDED
            cell_index=EL_BND(II)
            do edge_index = 1,IELEM(N,cell_index)%nonodes
                node_index = IELEM(N,cell_index)%nodes(edge_index)
                
                if (edge_index.lt.IELEM(N,cell_index)%nonodes) then
                    node_index_2 = (IELEM(N,cell_index)%nodes(edge_index + 1))
                else
                    node_index_2 = (IELEM(N,cell_index)%nodes(1))
                end if
                if (IELEM(N,cell_index)%INEIGHB(edge_index).ne.N) then
                    ! local_nodes(node_index)%num_cpus = 2
                    local_nodes(node_index)%boundary = max(local_nodes(node_index)%boundary, 1)
                    ! local_nodes(node_index_2)%num_cpus = 2
                    local_nodes(node_index_2)%boundary = max(local_nodes(node_index_2)%boundary, 1)
                end if
                if (ielem(n,cell_index)%ibounds(edge_index).gt.0) then
                    if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.5) then
                        local_nodes(node_index)%boundary = 2
                        local_nodes(node_index_2)%boundary = 2
                    end if
                end if
            end do
        end do
  
        my_num_interface_nodes = 0

        do node_index=1,kmaxn 
            if ((local_nodes(node_index)%boundary.gt.0).or.(local_nodes(node_index)%num_local_neighbours.eq.0)) then
                my_num_interface_nodes = my_num_interface_nodes + 1
            end if
        end do
        allocate(local_interface_nodes(my_num_interface_nodes))
        i = 0
        do node_index=1,kmaxn
            if ((local_nodes(node_index)%boundary.gt.0).or.(local_nodes(node_index)%num_local_neighbours.eq.0)) then
                i = i+1
                local_interface_nodes(i) = node_index
            end if
        end do
        if (i.ne.my_num_interface_nodes) then
            print *, "missmatch in the number of interface nodes on CPU", N
        ! else
        !     print *, my_num_interface_nodes, "= num interface nodes on CPU", N
        end if

        Call MPI_ALLGATHER(my_num_interface_nodes, 1, MPI_INT, num_interface_nodes, 1, MPI_INT, MPI_COMM_WORLD, IERROR)

        my_num_interface_nodes = num_interface_nodes(N)

        do cpu_index = 0, isize-1
            allocate(buff(cpu_index)%coords(dimensiona*num_interface_nodes(cpu_index)))
            allocate(buff(cpu_index)%count(num_interface_nodes(cpu_index)))
        end do

        do i = 1, my_num_interface_nodes 
            node_index = local_interface_nodes(i)

            do j = 1, dimensiona
                buff(N)%coords(dimensiona*(i-1)+j) = local_nodes(node_index)%positions(1,j)
            end do
            buff(N)%count(i) = local_nodes(node_index)%Num_Local_Neighbours

            ! print *, "on CPU", N, i, "th interface node has coordinates", local_nodes(node_index)%positions(1,:), "and", local_nodes(node_index)%Num_Local_Neighbours, "local cell neighbours"
        end do

        num_requests = 0
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                num_requests = num_requests + 1
                CALL MPI_ISENDRECV(buff(N)%coords(1:my_num_interface_nodes*dimensiona), my_num_interface_nodes*dimensiona, MPI_DOUBLE_PRECISION, cpu_index, 789,&
                        buff(cpu_index)%coords(1:num_interface_nodes(cpu_index)*dimensiona), num_interface_nodes(cpu_index)*dimensiona, MPI_DOUBLE_PRECISION, cpu_index, 789,&
                        MPI_COMM_WORLD, requests(num_requests), IERROR)

                num_requests = num_requests + 1
                CALL MPI_ISENDRECV(buff(N)%count(1:my_num_interface_nodes), my_num_interface_nodes, MPI_INTEGER, cpu_index, 987,&
                        buff(cpu_index)%count(1:num_interface_nodes(cpu_index)), num_interface_nodes(cpu_index), MPI_INTEGER, cpu_index, 987,&
                        MPI_COMM_WORLD, requests(num_requests), IERROR)
            end if
        end do

        CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

        do i = 1, my_num_interface_nodes 
            node_index = local_interface_nodes(i)

            do j = 1, dimensiona
                buff(N)%coords(dimensiona*(i-1)+j) = local_nodes(node_index)%positions(1,j)
            end do
            buff(N)%count(i) = local_nodes(node_index)%Num_Local_Neighbours
        end do
       
        do cpu_index = 0, isize-1
            do i = 1, num_interface_nodes(cpu_index)
                do iter = 1, my_num_interface_nodes
                    node_index = local_interface_nodes(iter)
                    if (.not.((cpu_index.eq.N).and.(iter.eq.i))) then ! skip myself in my own buffer

                        distance = point_distance(local_nodes(node_index)%positions(1,1:dimensiona), buff(cpu_index)%coords((dimensiona*(i-1))+1: (dimensiona*(i-1))+dimensiona))

                        close_enough = .false.
                        if (distance.lt.tolerance) then
                            close_enough = .true.
                        end if

                        if (close_enough) then
                            local_nodes(node_index)%num_neighbours = local_nodes(node_index)%num_neighbours + buff(cpu_index)%count(i)
                            local_nodes(node_index)%num_cpus = local_nodes(node_index)%num_cpus + 1
                            ! print *, "on CPU", N, node_index, "matched with", i, "th node from", cpu_index 
                        end if
                    end if
                end do
            end do
        end do

        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "local_nodes(", node_index, ")%num_cpus = ", local_nodes(node_index)%num_cpus
            allocate(local_nodes(node_index)%rcv_offsets(local_nodes(node_index)%num_cpus))
            allocate(local_nodes(node_index)%snd_offsets(local_nodes(node_index)%num_cpus))
            local_nodes(node_index)%num_cpus = 0
        end do

        do cpu_index = 0, isize-1
            node_rcv_count(cpu_index) = 0
            node_snd_count(cpu_index) = 0

            !counter = 0
            do i = 1, num_interface_nodes(cpu_index)
                ! counter = 0
                do iter = 1, my_num_interface_nodes 
                    if (.not.((cpu_index.eq.n).and.(i.eq.iter))) then ! skip myself in my buffer
                        node_index = local_interface_nodes(iter)

                        distance = point_distance(local_nodes(node_index)%positions(1,1:dimensiona), buff(cpu_index)%coords(dimensiona*(i-1)+1: dimensiona*(i-1)+dimensiona))

                        close_enough = .false.
                        if (distance.lt.tolerance) then
                            close_enough = .true.
                        end if

                        if (close_enough) then
                                
                            local_nodes(node_index)%num_cpus = local_nodes(node_index)%num_cpus + 1
                            k = local_nodes(node_index)%num_cpus
                                
                            local_nodes(node_index)%rcv_offsets(k)%cpu = cpu_index
                            local_nodes(node_index)%rcv_offsets(k)%lower = node_rcv_count(cpu_index) + 1
                            local_nodes(node_index)%rcv_offsets(k)%upper = node_rcv_count(cpu_index) + buff(cpu_index)%count(i)
                            node_rcv_count(cpu_index) = node_rcv_count(cpu_index) + buff(cpu_index)%count(i)

                            local_nodes(node_index)%snd_offsets(k)%cpu = cpu_index
                            local_nodes(node_index)%snd_offsets(k)%lower = node_snd_count(cpu_index) + 1
                            local_nodes(node_index)%snd_offsets(k)%upper = node_snd_count(cpu_index) + local_nodes(node_index)%num_local_neighbours
                            node_snd_count(cpu_index) = node_snd_count(cpu_index) + local_nodes(node_index)%num_local_neighbours
                            ! print *, "on CPU", N, node_index, "rematched with", i, "th node from", cpu_index 
                        end if
                    end if
                end do
            end do
        end do

        ! node_snd_count(:) = 0

        ! do iter = 1,my_num_interface_nodes 
        !     node_index = local_interface_nodes(iter)
        !     do k = 1, local_nodes(node_index)%num_cpus-1
        !         cpu_index = local_nodes(node_index)%rcv_offsets(k)%cpu
        !         local_nodes(node_index)%snd_offsets(k)%cpu = cpu_index
        !         local_nodes(node_index)%snd_offsets(k)%lower = node_snd_count(cpu_index) + 1
        !         local_nodes(node_index)%snd_offsets(k)%upper = node_snd_count(cpu_index) + local_nodes(node_index)%Num_Local_Neighbours
        !         node_snd_count(cpu_index) = node_snd_count(cpu_index) + local_nodes(node_index)%Num_Local_Neighbours
        !     end do 
        ! end do

        do cpu_index = 0, isize-1
            ! if (cpu_index.ne.N) then
                allocate(node_snd_buffer(cpu_index)%data(node_snd_count(cpu_index)*num_values_to_send_per_node))
                ! print *, "on CPU", N, "node_snd_count(", cpu_index, ") = ", node_snd_count(cpu_index)
                allocate(node_rcv_buffer(cpu_index)%data(node_rcv_count(cpu_index)*num_values_to_send_per_node))
                ! print *, "on CPU", N, "node_rcv_count(", cpu_index, ") = ", node_rcv_count(cpu_index)
            ! end if
        end do
    !$omp end master

    deallocate(helper)
    do cpu_index = 0, isize-1
        deallocate(buff(cpu_index)%coords)
        deallocate(buff(cpu_index)%count)
    end do

end subroutine establish_node_neighbours





subroutine find_node_velocities(position_index, d_t, N)
    implicit none
    integer,intent(in)::position_index, N
    real,intent(in)::d_t
    real::x,y
    integer::node_index, cell_index, kmaxe

    KMAXE=XMPIELRANK(N)

    !$omp do
    do node_index = 1, kmaxn
        local_nodes(node_index)%velocity(:) = zero
    end do
    !$omp end do

    !$omp barrier

    if ((governingequations.eq.3).and.(INITCOND.eq.3)) then
        !$omp do
        do node_index = 1, kmaxn
            x = local_nodes(node_index)%positions(position_index,1)
            y = local_nodes(node_index)%positions(position_index,2)
            local_nodes(node_index)%velocity(1) = 0.5 - Y
            local_nodes(node_index)%velocity(2) = x - 0.5
        end do
        !$omp end do
    else
        ! !$omp do
        ! !$omp master
        ! do cell_index = 1, kmaxe
        !     call FirstOrderNodeSolverArithmetic(cell_index, N)
        ! end do
        ! !$omp end master
        ! !$omp end do

        ! !$omp do
        ! do node_index = 1, kmaxn
        !     x = local_nodes(node_index)%positions(position_index,1)
        !     y = local_nodes(node_index)%positions(position_index,2)
        !     local_nodes(node_index)%velocity(1) = 0.5 - (x / 9.0)
        !     local_nodes(node_index)%velocity(2) = 0.0
        ! end do
        ! !$omp end do
        if ((moving_mesh_mode.eq.1).or.(moving_mesh_mode.eq.2)) then
            call FirstOrderNodeAverage(1, N)
        else if (moving_mesh_mode.eq.3) then
            call FirstOrderNodeAverage_withVF(1, N)
        else if (moving_mesh_mode.eq.4) then
            call directReALE_node_velocity(1, position_index, d_t, N)
        else
            print *, "invalid moving mesh mode"
            call abort()
        end if
        call ModifyLagrangianVelocities(N)
    endif

    !$omp barrier

    call enforce_node_velocity_BC(position_index, N)

    !$omp barrier

end subroutine





! SUBROUTINE FirstOrderNodeSolverArithmetic(cell_index, N)
!     IMPLICIT NONE
!     INTEGER,INTENT(IN)::cell_index, N
!     INTEGER::J, K, L, M, num_neighbours, rho_index,d, node_index, node ! , I, NJ
!     REAL,DIMENSION(4,30,NOF_VARIABLES)::VERTEX_NEIGHBOURS_VALS	!MAXIMUM ALLOCATION NOT EFFICIENT
!     real,dimension(1:dimensiona)::old_node_velocity
!     real::epsilon
!     rho_index = 1
!     num_neighbours = 0
!     epsilon = 0.000001
!     ! I=ICONSIDERED
!     ! WRITE(680+N,*) "ELEMENT NUMBER", IELEM(N,cell_index)%IHEXGL
!     ! print *, "N", n, "cell_index", cell_index, "nonodes", IELEM(N,cell_index)%nonodes
!     ! call flush

!     DO node_index = 1, IELEM(N,cell_index)%nonodes
!         IF (ILOCAL_RECON3(cell_index)%LOCAL.eq.1)then
!             DO J= 1, IELEM(N,cell_index)%NOJECOUNT(node_index)
!                 VERTEX_NEIGHBOURS_VALS(node_index, J, 1:NOF_VARIABLES) = U_C(ILOCAL_RECON3(cell_index)%IHEXL(1,IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)))%val(1,1:nof_variables)
!                 ! ILOCAL_RECON3(I)%IHEXL(1,L)	!LOCAL NUMBERING IN MY CPU
!                 ! ILOCAL_RECON3(I)%IHEXG(1,L)	!GLOBAL NUMBERING IN MY CPU
!                 WRITE(680+N,*) "NODE", node_index, "NODE NEIGHBOUR", J, "LOCAL NUMBER",IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J), "VALUES", VERTEX_NEIGHBOURS_VALS(node_index,J,1:NOF_VARIABLES)
!             END DO
!         ELSE
! 			DO J=1,IELEM(N,cell_index)%NOJECOUNT(node_index)
! 				IF (ILOCAL_RECON3(cell_index)%IHEXB(1,IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)).EQ.N)THEN
! 					VERTEX_NEIGHBOURS_VALS(node_index,J,1:NOF_VARIABLES) = U_C(ILOCAL_RECON3(cell_index)%IHEXL(1,IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)))%val(1,1:nof_variables)
! 				ELSE
! 					! write(500+n,*) ILOCAL_RECON3(cell_index)%IHEXN(1,IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)),IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)
! 					VERTEX_NEIGHBOURS_VALS(node_index,J,1:NOF_VARIABLES) = IEXSOLHIR(ILOCAL_RECON3(cell_index)%IHEXN(1,IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)))%SOL(ILOCAL_RECON3(cell_index)%IHEXL(1,IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J)),1:nof_variables)
! 				END IF
! 			    WRITE(680+N,*) "NODE", node_index, "NODE NEIGHBOUR", J, "LOCAL NUMBER", IELEM(N,cell_index)%NODES_NEIGHBOURS(node_index,J), "VALUES", VERTEX_NEIGHBOURS_VALS(node_index,J,1:NOF_VARIABLES)
! 			END DO
! 	    END IF
!     END DO
!     ! NOW YOU HAVE COLLECTED FOR EVERY VERTEX THE VALUES

!     DO node_index = 1, IELEM(N,cell_index)%nonodes
!         node = ielem(N,cell_index)%nodes(node_index)
!         old_node_velocity(1:dimensiona) = local_nodes(node)%velocity(1:dimensiona)
!         local_nodes(node)%velocity(:) = zero
!         if (IELEM(N,cell_index)%NOJECOUNT(node_index).lt.3) then
!             print *, "This node has too few (", IELEM(N,cell_index)%NOJECOUNT(node_index), ") neighbouring cells:", node_index, "in cell", cell_index, N 
!         end if

!         DO J = 1, IELEM(N,cell_index)%NOJECOUNT(node_index)
!             if (VERTEX_NEIGHBOURS_VALS(node_index, J, rho_index).eq.Zero) then
!                 print *, "dividing by zero in cell", cell_index, "node", node_index, N 
!             end if
!             do d = 1, dimensiona
!                 local_nodes(node)%velocity(d) = local_nodes(node)%velocity(d) &
!                         + (VERTEX_NEIGHBOURS_VALS(node_index, J, d+rho_index) / VERTEX_NEIGHBOURS_VALS(node_index, J, rho_index))
!             end do
!         END DO
!         ! do d = 1, dimensiona
!         !     local_nodes(node)%velocity(d) = local_nodes(node)%velocity(d) &
!         !             + (U_C(cell_index)%val(1,d+rho_index) / U_C(cell_index)%val(1,rho_index))
!         ! end do
!         do d = 1, dimensiona
!             ! local_nodes(node)%velocity(d) = local_nodes(node)%velocity(d) / (real(IELEM(N,cell_index)%NOJECOUNT(node_index) + 1))
!             local_nodes(node)%velocity(d) = local_nodes(node)%velocity(d) / real(IELEM(N,cell_index)%NOJECOUNT(node_index))
!             local_nodes(node)%velocity(d) = local_nodes(node)%velocity(d) * mesh_volocity_multiple
            
!             if (local_nodes(node)%velocity(d).ne.local_nodes(node)%velocity(d)) then
!                 print *, "NaN node velocity in cell", cell_index, "node", node, N
!                 ! call abort
!             end if
!             if (old_node_velocity(d).ne.zero) then
!                 if (abs(local_nodes(node)%velocity(d)-old_node_velocity(d)).gt.epsilon) then
!                     ! print *, "node velocity changed at node", node, "coordinate", d, "old", old_node_velocity(d), "new", local_nodes(node)%velocity(d)  
!                     ! call abort
!                 end if
!             end if
!         end do
!         write(120+n,*) IELEM(N,cell_index)%ihexgl,node_index,local_nodes(node)%velocity(1:dimensiona)
!     END DO

! END SUBROUTINE FirstOrderNodeSolverArithmetic





SUBROUTINE FirstOrderNodeAverage(stage, N)
    implicit none
    integer,intent(in)::stage, N
    integer::i, j, k, iter, node_index, cell_index, cpu_index, cpu, index
    integer:: M
    real::rho,v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, gammal
    integer::rho_index

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()

    if (num_values_to_send_per_node.ne.dimensiona) then
        print *,"something went wrong sorry :("
        call abort
    end if

    rho_index = 1

    ! do cpu_index = 0, isize-1
    !     index = 0
    !     do i = 1, node_snd_count(cpu_index)
    !         do j = 1, dimensiona
    !             index = index +1
    !             node_snd_buffer(cpu_index)%data(index) = 1000.0
    !         end do
    !     end do
    !     if (index.ne.dimensiona*node_snd_count(cpu_index)) then
    !         print *, "something went wrong with rezeroing the send buffer"
    !     else
    !         print *, "snd_buffer on CPU", N, "thread", M, "to", cpu_index, "reset with", node_snd_count(cpu_index)*dimensiona, "values"
    !     end if
    ! end do
    ! !$omp barrier 

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * num_values_to_send_per_node
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index) = copy(rho_index + k)
                    end do
                end do
            end do
        end do
    !$omp end do

    ! do cpu_index = 0, isize-1
    !     index = 0
    !     do i = 1, node_rcv_count(cpu_index)
    !         do j = 1, dimensiona
    !             index = index +1
    !             node_rcv_buffer(cpu_index)%data(index) = 1000000.0
    !         end do
    !     end do
    !     if (index.ne.dimensiona*node_rcv_count(cpu_index)) then
    !         print *, "something went wrong with rezeroing the receive buffer"
    !     else
    !         print *, "rcv_buffer on CPU", N, "thread", M, "from", cpu_index, "reset with", node_rcv_count(cpu_index)*dimensiona, "values"
    !     end if
    ! end do

    !$omp barrier

    num_requests = 0
    !$omp master
        ! print *, "inside send_rcv part on CPU", N
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                if (node_snd_count(cpu_index).gt.0) then
                    count = node_snd_count(cpu_index) * num_values_to_send_per_node
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * num_values_to_send_per_node
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * num_values_to_send_per_node
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

            local_nodes(node_index)%velocity(:) = 0.0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)
                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                do j = 1, dimensiona
                    v = copy(rho_index+j)
                    local_nodes(node_index)%velocity(j) = local_nodes(node_index)%velocity(j) + v
                end do
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        v = node_rcv_buffer(cpu_index)%data(index)
                        local_nodes(node_index)%velocity(k) = local_nodes(node_index)%velocity(k) + v
                    end do
                end do
            end do

            do j = 1, dimensiona
                if (real(local_nodes(node_index)%num_neighbours).eq.0.0) then
                    print*,"dividing by real zero num_neighbours in node", node_index, "CPU", n
                end if
                if (local_nodes(node_index)%num_neighbours.eq.0) then
                    print*,"dividing by integer zero num_neighbours in node", node_index, "CPU", n
                end if
                local_nodes(node_index)%velocity(j) = local_nodes(node_index)%velocity(j) / real(local_nodes(node_index)%num_neighbours)
                if (moving_mesh_mode.eq.1) then
                    local_nodes(node_index)%velocity(j) = local_nodes(node_index)%velocity(j) * mesh_volocity_multiple
                end if
            end do

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE FirstOrderNodeAverage




SUBROUTINE FirstOrderNodeAverage_withVF(stage, N)
    implicit none
    integer,intent(in)::stage, N
    integer::i, j, k, iter, node_index, cell_index, cpu_index, cpu, index
    integer:: M
    real::rho, v, vf, vf_closest_to_half, local_mesh_velocity_multiple
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, gammal
    integer::rho_index, vf_index

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()

    if (num_values_to_send_per_node.ne.(dimensiona+1)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    rho_index = 1
    vf_index = dimensiona+5

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * num_values_to_send_per_node
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index) = copy(rho_index + k)
                    end do
                    index = index +1
                    if (index.gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                        print *, "copying too much data to send buffer from", N, "to", cpu 
                    end if
                    node_snd_buffer(cpu)%data(index) = copy(vf_index)
                end do
            end do
        end do
    !$omp end do

    ! do cpu_index = 0, isize-1
    !     index = 0
    !     do i = 1, node_rcv_count(cpu_index)
    !         do j = 1, dimensiona
    !             index = index +1
    !             node_rcv_buffer(cpu_index)%data(index) = 1000000.0
    !         end do
    !     end do
    !     if (index.ne.dimensiona*node_rcv_count(cpu_index)) then
    !         print *, "something went wrong with rezeroing the receive buffer"
    !     else
    !         print *, "rcv_buffer on CPU", N, "thread", M, "from", cpu_index, "reset with", node_rcv_count(cpu_index)*dimensiona, "values"
    !     end if
    ! end do

    !$omp barrier

    num_requests = 0
    !$omp master
        ! print *, "inside send_rcv part on CPU", N
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                if (node_snd_count(cpu_index).gt.0) then
                    count = node_snd_count(cpu_index) * num_values_to_send_per_node
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * num_values_to_send_per_node
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * num_values_to_send_per_node
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

            local_nodes(node_index)%velocity(:) = 0.0
            vf_closest_to_half = zero

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)
                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                do j = 1, dimensiona
                    v = copy(rho_index+j)
                    local_nodes(node_index)%velocity(j) = local_nodes(node_index)%velocity(j) + v
                end do
                vf = copy(vf_index)
                if (abs(vf-0.5).lt.abs(vf_closest_to_half-0.5)) then
                    vf_closest_to_half = vf
                end if
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        v = node_rcv_buffer(cpu_index)%data(index)
                        local_nodes(node_index)%velocity(k) = local_nodes(node_index)%velocity(k) + v
                    end do
                    index = index+1
                    vf = node_rcv_buffer(cpu_index)%data(index)
                    if (abs(vf-0.5).lt.abs(vf_closest_to_half-0.5)) then
                        vf_closest_to_half = vf
                    end if
                end do
            end do

            do j = 1, dimensiona
                if (real(local_nodes(node_index)%num_neighbours).eq.0.0) then
                    print*,"dividing by real zero num_neighbours in node", node_index, "CPU", n
                end if
                if (local_nodes(node_index)%num_neighbours.eq.0) then
                    print*,"dividing by integer zero num_neighbours in node", node_index, "CPU", n
                end if
                local_nodes(node_index)%velocity(j) = local_nodes(node_index)%velocity(j) / real(local_nodes(node_index)%num_neighbours)

                if (vf_closest_to_half.le.0.5) then
                    local_mesh_velocity_multiple = 2.0*vf_closest_to_half
                else
                    local_mesh_velocity_multiple = 2.0*(1.0-vf_closest_to_half)
                end if
                ! local_mesh_velocity_multiple = sin(pi*vf_closest_to_half)
                local_nodes(node_index)%velocity(j) = local_nodes(node_index)%velocity(j) * local_mesh_velocity_multiple
            end do

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE FirstOrderNodeAverage_withVF





SUBROUTINE directReALE_node_velocity(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, counter, node_index, cell_index, cpu_index, cpu, index, rho_index
    integer:: M
    real::rho, v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, dummy_gammal
    real,dimension(1:dimensiona)::helper_centre_position
    real,dimension(1:max_num_node_neighbours,1:dimensiona)::centre_positions


    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()
    
    rho_index = 1

    if (num_values_to_send_per_node.ne.dimensiona) then
        print *,"something went wrong sorry :("
        call abort
    end if

    ! do cpu_index = 0, isize-1
    !     index = 0
    !     do i = 1, node_snd_count(cpu_index)
    !         do j = 1, dimensiona
    !             index = index +1
    !             node_snd_buffer(cpu_index)%data(index) = 1000.0
    !         end do
    !     end do
    !     if (index.ne.dimensiona*node_snd_count(cpu_index)) then
    !         print *, "something went wrong with rezeroing the send buffer"
    !     else
    !         print *, "snd_buffer on CPU", N, "thread", M, "to", cpu_index, "reset with", node_snd_count(cpu_index)*dimensiona, "values"
    !     end if
    ! end do
    ! !$omp barrier 

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * num_values_to_send_per_node
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    call cons2prim(N, copy, dummuy_MP_PINFl, dummy_gammal)
                    helper_centre_position(1) = ielem(N, cell_index)%xxc + (copy(rho_index+1)*d_t)
                    helper_centre_position(2) = ielem(N, cell_index)%yyc + (copy(rho_index+2)*d_t)
                    if (dimensiona.eq.3) then
                        helper_centre_position(3) = ielem(N, cell_index)%zzc + (copy(rho_index+3)*d_t)
                    end if
                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index) = helper_centre_position(k)
                    end do
                end do
            end do
        end do
    !$omp end do

    ! do cpu_index = 0, isize-1
    !     index = 0
    !     do i = 1, node_rcv_count(cpu_index)
    !         do j = 1, dimensiona
    !             index = index +1
    !             node_rcv_buffer(cpu_index)%data(index) = 1000000.0
    !         end do
    !     end do
    !     if (index.ne.dimensiona*node_rcv_count(cpu_index)) then
    !         print *, "something went wrong with rezeroing the receive buffer"
    !     else
    !         print *, "rcv_buffer on CPU", N, "thread", M, "from", cpu_index, "reset with", node_rcv_count(cpu_index)*dimensiona, "values"
    !     end if
    ! end do

    !$omp barrier

    num_requests = 0
    !$omp master
        ! print *, "inside send_rcv part on CPU", N
        do cpu_index = 0, isize-1
            if (cpu_index.ne.N) then
                if (node_snd_count(cpu_index).gt.0) then
                    count = node_snd_count(cpu_index) * num_values_to_send_per_node
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * num_values_to_send_per_node
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * num_values_to_send_per_node
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

            local_nodes(node_index)%velocity(:) = 0.0
            counter = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1
                cell_index = local_nodes(node_index)%local_neighbours(i)

                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummuy_MP_PINFl, dummy_gammal)
                helper_centre_position(1) = ielem(N, cell_index)%xxc + (copy(rho_index+1)*d_t)
                helper_centre_position(2) = ielem(N, cell_index)%yyc + (copy(rho_index+2)*d_t)
                if (dimensiona.eq.3) then
                    helper_centre_position(3) = ielem(N, cell_index)%zzc + (copy(rho_index+3)*d_t)
                end if

                centre_positions(counter, :) = helper_centre_position(:)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    counter = counter+1
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index)
                    end do
                end do
            end do

            if (counter.ne.local_nodes(node_index)%num_neighbours) then
                print *, "something went wrong counter =/= local_nodes(node_index)%num_neighbours"
            end if
            call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%velocity(1:dimensiona))

            local_nodes(node_index)%velocity(1:dimensiona) = (local_nodes(node_index)%velocity(1:dimensiona) - local_nodes(node_index)%positions(position_index,1:dimensiona)) / d_t

            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) * mesh_volocity_multiple
        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE directReALE_node_velocity





subroutine ModifyLagrangianVelocities(N)
    implicit none
    integer,intent(in)::N
    integer::node_index
    integer::i
    real,dimension(1:dimensiona)::initial_velocity1, initial_velocity2
    real,dimension(1:dimensiona)::option1, option2

    initial_velocity1 = zero
    initial_velocity2 = zero
    
    if (moving_mesh_mode.eq.2) then
        if ((governingequations.ne.3).and.(initcond.eq.102)) then
            initial_velocity2(1) =  8.25*cos(pi/6.0d0)
            ! initial_velocity2(2) = -8.25*sin(pi/6.0d0)
        end if

        if ((governingequations.ne.3).and.(initcond.eq.101)) then
            initial_velocity2(1) = 2.6294d0
            initial_velocity2(2) = zero
        end if

        !$omp do
        do node_index = 1, kmaxn
            if ((governingequations.ne.3).and.(initcond.eq.102)) then
                local_nodes(node_index)%velocity(2) = zero
                ! print *, "zeroing y node velocity"
            end if
            option1 = local_nodes(node_index)%velocity(1:dimensiona) - initial_velocity1
            option2 = initial_velocity2 - local_nodes(node_index)%velocity(1:dimensiona)
            local_nodes(node_index)%velocity(1:dimensiona) = min_abs(option1, option2)
        end do
        !$omp end do
    end if

end subroutine ModifyLagrangianVelocities





subroutine enforce_node_velocity_BC(position_index, N)
    implicit none
    integer,intent(in)::position_index, n
    real::x, y, epsilon, v2
    integer:: node_index
    epsilon = 0.000001

    if (initcond.eq.101) then
        !$omp do
        do node_index = 1, kmaxn
            x = local_nodes(node_index)%positions(position_index,1)
            y = local_nodes(node_index)%positions(position_index,2)
            v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
            if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
                print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
            end if
            if (v2.gt.((5.0*mesh_volocity_multiple)*(5.0*mesh_volocity_multiple))) then
                print *, "too high node speed^2", v2, "in node", node_index, N
            end if
            if ((y.lt.epsilon).or.(y.gt.(1.0-epsilon))) then
                local_nodes(node_index)%velocity(2) = 0.0
            end if
            if (y.lt.0.0) then
                local_nodes(node_index)%velocity(2) = (-1.0)*y/dt
            end if
            if (y.gt.1.0) then
                local_nodes(node_index)%velocity(2) = (y - 1.0)/dt
            end if
            if (moving_mesh_mode.eq.4) then
                if ((x.le.-4.5+epsilon).or.(x.ge.4.5-epsilon)) then
                    local_nodes(node_index)%velocity(1) = 0.0
                end if
            end if

        end do
        !$omp end do
    
    else if (initcond.eq.102) then
        !$omp do
        do node_index = 1, kmaxn
            x = local_nodes(node_index)%positions(position_index,1)
            y = local_nodes(node_index)%positions(position_index,2)
            v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
            if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
                print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
            end if
            ! if (v2.gt.((10.0*mesh_volocity_multiple)*(10.0*mesh_volocity_multiple))) then
            !     print *, "too high node speed^2", v2, "in node", node_index, N
            ! end if
            if (moving_mesh_mode.eq.4) then
                if ((x.le.epsilon).or.(x.ge.4.0-epsilon)) then
                    local_nodes(node_index)%velocity(1) = 0.0
                end if
                if ((y.le.epsilon).or.(y.ge.1.0-epsilon)) then
                    local_nodes(node_index)%velocity(2) = 0.0
                end if
            end if

        end do
        !$omp end do

    else if (initcond.eq.470) then
        !$omp do
        do node_index = 1, kmaxn
            x = local_nodes(node_index)%positions(position_index,1)
            y = local_nodes(node_index)%positions(position_index,2)
            v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
            if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
                print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
            end if
            ! if (v2.gt.((10.0*mesh_volocity_multiple)*(10.0*mesh_volocity_multiple))) then
            !     print *, "too high node speed^2", v2, "in node", node_index, N
            ! end if
            if (moving_mesh_mode.eq.4) then
                if ((x.le.epsilon).or.(x.ge.7.0-epsilon)) then
                    local_nodes(node_index)%velocity(1) = 0.0
                end if
                if ((y.le.epsilon).or.(y.ge.3.0-epsilon)) then
                    local_nodes(node_index)%velocity(2) = 0.0
                end if
            end if

        end do
        !$omp end do
    
    else if (initcond.eq.405) then
        !$omp do
        do node_index = 1, kmaxn
            x = local_nodes(node_index)%positions(position_index,1)
            y = local_nodes(node_index)%positions(position_index,2)
            v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
            if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
                print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
            end if
            ! if (v2.gt.((10.0*mesh_volocity_multiple)*(10.0*mesh_volocity_multiple))) then
            !     print *, "too high node speed^2", v2, "in node", node_index, N
            ! end if
            if (moving_mesh_mode.eq.4) then
                if ((x.le.-0.25+epsilon).or.(x.ge.0.25-epsilon)) then
                    local_nodes(node_index)%velocity(1) = 0.0
                end if
                if ((y.le.epsilon).or.(y.ge.0.1-epsilon)) then
                    local_nodes(node_index)%velocity(2) = 0.0
                end if
            end if

        end do
        !$omp end do
    end if

end subroutine enforce_node_velocity_BC





SUBROUTINE QuadraturePoint_velocity(N,IGQRULES,NODE_velocities,QPOINTS2D_velocity)
	!> @brief
	!> This subroutine computes the velocity quadrature points for a line and returns it in QPOINTS2D_velocity(DIM,QP)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,IGQRULES
	REAL,DIMENSION(1:2,1:DIMENSIONA),INTENT(IN)::NODE_velocities
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2),INTENT(INOUT)::QPOINTS2D_velocity
	real,dimension(1:2,1:2)::VVA,VVA1
	REAL,dimension(1)::DETA
	REAL,DIMENSION(1:4)::VVNXI
	real,dimension(1:ALLS)::VVwg
	real,dimension(1:ALLS)::VVR1
	real,dimension(1:igqrules)::vvwpox,vvnpox,vvwpoy,vvnpoy,vvwpoz,vvnpoz
	REAL::R,S,TX,a,b,c,d,e,f,G,H,K
	INTEGER::Kk,J,ii,ij,ik,count1

	SELECT CASE(IGQRULES)
	
	  case(1)

		VVR1(1)=0.0d0

	  case(2)

		a=-0.5773502691896257
		b=0.5773502691896257
		vvnpox(1)=a	;  vvnpox(2)=b		
  
  		count1=0
 		do ii=1,2
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
  	
	  CASE(3)

		a=0.0d0
		b=-0.7745966692414834
		c=0.7745966692414834
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c
		
		count1=0
		do ii=1,3
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
 		
	  CASE(4)
		a=-0.3399810435848563
		b=0.3399810435848563
		c=-0.8611363115940526
		d=0.8611363115940526
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d
		
		count1=0
		do ii=1,4
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
			
	  CASE(5)

		a=0.0d0
		b=-0.5384693101056831
		c=0.5384693101056831
		d=-0.9061798459386640
		e=0.9061798459386640
		vvnpox(1)=a	;  vvnpox(2)=b ;   vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e
		
		count1=0
		do ii=1,5
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do
		
	  CASE(6)

		a=0.6612093864662645
		b=-0.6612093864662645
		c=-0.2386191860831969
		d=0.2386191860831969
		e=-0.9324695142031521
		F=0.9324695142031521
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e ;  vvnpox(6)=f
		
		count1=0
		do ii=1,6
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do

	  CASE(7,8,9)

		a=0.0000000000000000
		b=-0.8360311073266358
		C=0.8360311073266358
		D=-0.9681602395076261
		E=0.9681602395076261
		F=-0.3242534234038089
		G=0.3242534234038089
		H=-0.6133714327005904
		K=0.6133714327005904
		vvnpox(1)=a	;  vvnpox(2)=b ;  vvnpox(3)=c ;  vvnpox(4)=d ;  vvnpox(5)=e ;  vvnpox(6)=f ;  vvnpox(7)=G ;  vvnpox(8)=H ;  vvnpox(9)=K
	
		count1=0
 		do ii=1,9
			count1=count1+1
			VVR1(count1)=vvnpox(ii)
		end do	
			
	END SELECT

	QPOINTS2D_velocity(:,:)=0.0d0

	do kk=1,IGQRULES
		R=VVR1(kk)		
		QPOINTS2D_velocity(:,kk)=((NODE_velocities(1,1:2)+NODE_velocities(2,1:2))/2.0D0)+(R*((NODE_velocities(2,1:2)-NODE_velocities(1,1:2))/2.0D0))
	END DO

END SUBROUTINE QuadraturePoint_velocity





SUBROUTINE Find_QP_velocities(N)
    !> @brief
    !> Subroutine for computing and storing velocities of gaussian quadrature points at cell interfaces
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX
    INTEGER::ICONSIDERED,FACEX,POINTX
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D_velocity
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    ! REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
    REAL,DIMENSION(1:8,1:DIMENSIONA)::Node_velocities
    REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
    
    KMAXE=XMPIELRANK(N)
    
    if (dimensiona.eq.3)then
        ! DO I=1,KMAXE
        !     IF (IELEM(N,I)%ISHAPE.EQ.2)THEN
        !         ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !         IF (SRFG.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !         END IF
        !         IF (MRF.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_TRIANGLE,3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_ORIGIN(1:3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_VELOCITY(1:3))
        !         END IF
        !     ELSE
        !         ALLOCATE(ILOCAL_RECON3(I)%QPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
        !         IF (SRFG.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_QUAD,3))
        !         END IF
        !         IF (MRF.EQ.1)THEN
        !             ALLOCATE(ILOCAL_RECON3(I)%RPOINTS(IELEM(N,I)%IFCA,QP_QUAD,3))
        !             ALLOCATE(ILOCAL_RECON3(I)%ROTVEL(IELEM(N,I)%IFCA,QP_QUAD,3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_ORIGIN(1:3))
        !             ALLOCATE (ILOCAL_RECON3(I)%MRF_VELOCITY(1:3))
        !         END IF
        !     END IF
        !     ICONSIDERED=I
        !     DO L=1,IELEM(N,I)%IFCA
        !         IDUMMY=0
        !         if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
        !             IF (IELEM(N,I)%IBOUNDS(l).GT.0) THEN	!CHECK FOR BOUNDARIES
        !                 if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50)) then	!PERIODIC IN OTHER CPU
        !                     IDUMMY=1
        !                 END IF
        !             END IF	
        !             if (ielem(n,i)%types_faces(L).eq.5) then
        !                 iqp=qp_quad
        !                 NND=4
        !                 IF (IDUMMY.EQ.0)THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 ELSE
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     ! do K=1,nnd
        !                     !     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     ! END DO
        !                 END IF
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3
        !                 IF (IDUMMY.EQ.0)THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 ELSE
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     ! do K=1,nnd
        !                     !     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     ! END DO
        !                 END IF
                            
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         else
        !             if (ielem(n,i)%types_faces(L).eq.5)then
        !                 iqp=qp_quad
        !                 NND=4
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO 
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3 
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO  
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         end if
                
        !         do NGP=1,iqp			!for gqp
        !             ! ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !             IF (SRFG.EQ.1) THEN
        !                 ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !                 POX(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)-SRF_ORIGIN(1:3)
        !                 POY(1:3)=SRF_VELOCITY(1:3)
        !                 ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)=VECT_FUNCTION(POX,POY)
        !             END IF
                        
        !             IF (MRF.EQ.1) THEN
        !                 ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !                 POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
        !                 POY(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)
        !                 ICONSIDERED=I
        !                 FACEX=L
        !                 POINTX=NGP
        !                 CALL MRFSWITCH(N,ICONSIDERED,FACEX,POINTX,pox,poy)
                        
        !             END IF 
        !         END DO	!NGP
        !     END DO
    
        !     DO L=1,IELEM(N,I)%IFCA
        !         IDUMMY=0
        !         if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
        !             IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
        !                 if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
        !                     IDUMMY=1
        !                 END IF
        !             END IF	
        !             if (ielem(n,i)%types_faces(L).eq.5) then
        !                 iqp=qp_quad
        !                 NND=4
        !                 IF (IDUMMY.EQ.0) THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO	!NGP
        !                 ELSE
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 END IF
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3
        !                 IF (IDUMMY.EQ.0)THEN
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        
        !                 ELSE ! 2 dimensions (Michael's edit: ?????)
        !                     facex=l;
        !                     CALL coordinates_face_PERIOD1(n,iconsidered,facex,VEXT,NODES_LIST)
        !                     do K=1,nnd
        !                         VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                     END DO
        !                 END IF
                                                
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         else
        !             if (ielem(n,i)%types_faces(L).eq.5)then
        !                 iqp=qp_quad
        !                 NND=4
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO 
        !                 call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             else
        !                 iqp=QP_TRIANGLE
        !                 NND=3 
        !                 do K=1,nnd
        !                     VEXT(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
        !                     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
        !                 END DO  
        !                 call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
        !             end if
        !         end if
            
        !         do NGP=1,iqp			!for gqp
        !             ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
        !         END DO	!NGP
        !     END DO
        ! END DO
          
    else ! 2 dimensions
        
        DO I=1,KMAXE
            
            ICONSIDERED=I
        
            DO L=1,IELEM(N,I)%IFCA
                IDUMMY=0
            
                ! if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                !     IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
                !         if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
                !             IDUMMY=1
                !         END IF
                !     END IF
            
                !     IQP=QP_LINE
                !     NND=2
                !     IF (IDUMMY.EQ.0)THEN
                !         DO K=1,NND
                !             ! VEXT(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(L,K))%CORD(1:dims)
                !             Node_velocities(k,1:2) = local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(1:2)
                !             !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                !                 ! VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                !             !END IF
                !         END DO
                !     ELSE
                !         ! facex=l;
                !         ! CALL coordinates_face_PERIOD2D1(n,iconsidered,facex,VEXT,NODES_LIST)
                !         ! DO K=1,NND
                !         !     !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                !         !         VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                !         !     !END IF
                !         ! END DO
                !     END IF
                !     CALL QUADRATURELINE(N,IGQRULES,Node_velocities,QPOINTS2D_velocity,WEQUA2D)
                ! ELSE
                    IQP=QP_LINE
                    NND=2
                    DO K=1,NND
                        Node_velocities(k,1:2) = local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(1:2)
                        ! print*, "velocity", local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(1), local_nodes(IELEM(N,I)%NODES_FACES(L,K))%VELOCITY(2), "expected", (0.5 - local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(1,2)), (local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(1,1) - 0.5)
                        ! !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                        !     VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                        ! !END IF
                    END DO
                    CALL QUADRATURELINE(N,IGQRULES,Node_velocities,QPOINTS2D_velocity,WEQUA2D)
                ! END IF
            
                DO NGP=1,iqp !for gqp
                    ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2)=QPOINTS2D_velocity(1:2,NGP) ! Storing surface quadrature points
                END DO !NGP
            END DO
        END DO
    END IF
END SUBROUTINE Find_QP_velocities





SUBROUTINE Find_QP_positions(N, node_position_index)
    !> @brief
    !> Subroutine for storing the gaussian quadrature points at the cell interfaces
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N, node_position_index
    INTEGER::I,K,KMAXE,IDUMMY,L,NND,IQP,NGP,IEX
    INTEGER::ICONSIDERED,FACEX,POINTX
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
    REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ

    !@TT: ensure allocation outside of this one so that you can call it as many times as possible

    KMAXE=XMPIELRANK(N)

    if (dimensiona.eq.3)then
        DO I=1,KMAXE

            ICONSIDERED=I
            DO L=1,IELEM(N,I)%IFCA
                IDUMMY=0
                if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                    IF (IELEM(N,I)%IBOUNDS(l).GT.0) THEN	!CHECK FOR BOUNDARIES
                        if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50)) then	!PERIODIC IN OTHER CPU
                            IDUMMY=1
                        END IF
                    END IF	
                    if (ielem(n,i)%types_faces(L).eq.5) then
                        iqp=qp_quad
                        NND=4
                        IF (IDUMMY.EQ.0)THEN
                            do K=1,nnd
                                VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                                ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            END DO
                        ELSE
                            facex=l;
                            CALL coordinates_face_PERIOD1_MovingMesh(n,iconsidered,facex,VEXT,NODES_LIST, node_position_index)
                            ! do K=1,nnd
                            !     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            ! END DO
                        END IF
                        call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    else
                        iqp=QP_TRIANGLE
                        NND=3
                        IF (IDUMMY.EQ.0)THEN
                            do K=1,nnd
                                VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                                ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            END DO
                        ELSE
                            facex=l;
                            CALL coordinates_face_PERIOD1_MovingMesh(n,iconsidered,facex,VEXT,NODES_LIST, node_position_index)
                            ! do K=1,nnd
                            !     VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            ! END DO
                        END IF
                            
                        call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    end if
                else
                    if (ielem(n,i)%types_faces(L).eq.5)then
                        iqp=qp_quad
                        NND=4
                        do K=1,nnd
                            VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                            ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                        END DO 
                        call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    else
                        iqp=QP_TRIANGLE
                        NND=3 
                        do K=1,nnd
                            VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                            ! VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                        END DO  
                        call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    end if
                end if
                
                do NGP=1,iqp			!for gqp
                    ! ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
                    IF (SRFG.EQ.1) THEN
                        ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
                        POX(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)-SRF_ORIGIN(1:3)
                        POY(1:3)=SRF_VELOCITY(1:3)
                        ILOCAL_RECON3(I)%ROTVEL(L,NGP,1:3)=VECT_FUNCTION(POX,POY)
                    END IF
                        
                    IF (MRF.EQ.1) THEN
                        ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
                        POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
                        POY(1:3)=ILOCAL_RECON3(I)%RPOINTS(L,NGP,1:3)
                        ICONSIDERED=I
                        FACEX=L
                        POINTX=NGP
                        CALL MRFSWITCH(N,ICONSIDERED,FACEX,POINTX,pox,poy)
                        
                    END IF 
                END DO	!NGP
            END DO

            DO L=1,IELEM(N,I)%IFCA
                IDUMMY=0
                if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                    IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
                        if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
                            IDUMMY=1
                        END IF
                    END IF	
                    if (ielem(n,i)%types_faces(L).eq.5) then
                        iqp=qp_quad
                        NND=4
                        IF (IDUMMY.EQ.0) THEN
                            do K=1,nnd
                                VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            END DO	!NGP
                        ELSE
                            facex=l;
                            CALL coordinates_face_PERIOD1_MovingMesh(n,iconsidered,facex,VEXT,NODES_LIST, node_position_index)
                            do K=1,nnd
                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:), VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            END DO
                        END IF
                        call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    else
                        iqp=QP_TRIANGLE
                        NND=3
                        IF (IDUMMY.EQ.0)THEN
                            do K=1,nnd
                                VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            END DO
        
                        ELSE ! 2 dimensions (Michael's edit: ?????)
                            facex=l;
                            CALL coordinates_face_PERIOD1_MovingMesh(n,iconsidered,facex,VEXT,NODES_LIST, node_position_index)
                            do K=1,nnd
                                VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                            END DO
                        END IF
                                                
                        call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    end if
                else
                    if (ielem(n,i)%types_faces(L).eq.5)then
                        iqp=qp_quad
                        NND=4
                        do K=1,nnd
                            VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                            VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                        END DO 
                        call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    else
                        iqp=QP_TRIANGLE
                        NND=3 
                        do K=1,nnd
                            VEXT(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                            VEXT(k,1:3)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:3)-ILOCAL_RECON3(I)%VEXT_REF(1:3))
                        END DO  
                        call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                    end if
                end if
            
                do NGP=1,iqp			!for gqp
                    ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:3)=QPOINTS2D(1:3,NGP)
                END DO	!NGP
            END DO
        END DO
        
    else ! 2 dimensions
        
        DO I=1,KMAXE
            ICONSIDERED=I
        
            DO L=1,IELEM(N,I)%IFCA
                IDUMMY=0
            
                if ((iperiodicity.eq.1).and.(ielem(n,i)%interior.eq.1))then	
                    IF (IELEM(N,I)%IBOUNDS(l).GT.0)THEN	!CHECK FOR BOUNDARIES
                        if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
                            IDUMMY=1
                        END IF
                    END IF
            
                    IQP=QP_LINE
                    NND=2
                    IF (IDUMMY.EQ.0)THEN
                        DO K=1,NND
                            VEXT(k,1:2)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                            !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                                VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                            !END IF
                        END DO
                    ELSE
                        facex=l;
                        CALL coordinates_face_PERIOD2D1_MovingMesh(n, iconsidered, facex, VEXT, NODES_LIST, node_position_index)
                        DO K=1,NND
                            !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                                VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                            !END IF
                        END DO
                    END IF
                    CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                ELSE
                    IQP=QP_LINE
                    NND=2
                    DO K=1,NND
                        VEXT(k,1:2)=local_nodes(IELEM(N,I)%NODES_FACES(L,K))%positions(node_position_index, 1:dims)
                        !IF (DG /= 1) THEN ! Only transforming to reference space if not DG
                            VEXT(k,1:2)=MATMUL(ILOCAL_RECON3(I)%INVCCJAC(:,:),VEXT(K,1:2)-ILOCAL_RECON3(I)%VEXT_REF(1:2))
                        !END IF
                    END DO
                    CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
                END IF
            
                DO NGP=1,iqp !for gqp
                    ILOCAL_RECON3(I)%QPOINTS(L,NGP,1:2)=QPOINTS2D(1:2,NGP) ! Storing surface quadrature points
                    ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2) = 0.0
                END DO !NGP
            END DO
        END DO
    END IF
    
END SUBROUTINE Find_QP_positions

  


subroutine coordinates_face_PERIOD1_MovingMesh(n, iconsidered, facex, VEXT, NODES_LIST, node_position_index)
	!> @brief
	!> This subroutine retrieve the nodes of periodic faces of elements in 3D
	IMPLICIT NONE
	integer,intent(in)::n, iconsidered, facex, node_position_index
	REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(INOUT)::VEXT
	REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(INOUT)::NODES_LIST
	integer::nnd
	integer::i,k
	real::tempxx
	i=iconsidered

	select case (ielem(n,i)%types_faces(facex))
	  case(5)
		nnd=4
	  case(6)
		nnd=3
	end select
	
	VEXT(1,1)=IELEM(N,I)%XXC
	VEXT(1,2)=IELEM(N,I)%YYC
	VEXT(1,3)=IELEM(N,I)%ZZC
	      
	do K=1,nnd
		NODES_LIST(k,1:3)=local_nodes(IELEM(N,I)%NODES_FACES(FACEX,K))%positions(node_position_index,1:dims)
	END DO
	do K=1,nnd
		IF(PER_ROT.EQ.0)THEN
			IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
				NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
				NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,3)-vext(1,3)).GT.zPER*oo2)THEN
				NODES_LIST(k,3)=NODES_LIST(k,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
			end if
		ELSE
			if (IELEM(n,i)%reorient(facex).eq.1) then
				if (ibound(n,ielem(n,i)%ibounds(facex))%icode.eq.5) then
					tempxx=NODES_LIST(k,1)
					NODES_LIST(k,1)=tempxx*cos(-angle_per)-sin(-angle_per)*NODES_LIST(k,2)
					NODES_LIST(k,2)=tempxx*sin(-angle_per)+cos(-angle_per)*NODES_LIST(k,2)
				else
					tempxx=NODES_LIST(k,1)
					NODES_LIST(k,1)=tempxx*cos(angle_per)-sin(angle_per)*NODES_LIST(k,2)
					NODES_LIST(k,2)=tempxx*sin(angle_per)+cos(angle_per)*NODES_LIST(k,2)
                end if
             end if
		END IF
	END DO
	      
	do K=1,nnd
		VEXT(K,1:3)=NODES_LIST(k,1:3)
	END DO
	      
end subroutine coordinates_face_PERIOD1_MovingMesh





subroutine coordinates_face_PERIOD2d1_MovingMesh(n, iconsidered, facex, VEXT, NODES_LIST, node_position_index)
	!> @brief
	!> This subroutine retrieve the nodes of periodic edges of elements in 2D
	IMPLICIT NONE
	integer,intent(in)::n, iconsidered, facex, node_position_index
	REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(INOUT)::VEXT
	REAL,DIMENSION(1:8,1:DIMENSIONA),INTENT(INOUT)::NODES_LIST
	integer::nnd
	integer::i,k
	i=iconsidered

	nnd=2
	      
	VEXT(1,1)=IELEM(N,I)%XXC
	VEXT(1,2)=IELEM(N,I)%YYC
	
	do K=1,nnd
		NODES_LIST(k,1:2)=local_nodes(IELEM(N,I)%NODES_FACES(FACEX,K))%positions(node_position_index, 1:dims)
	END DO
	do K=1,nnd
		IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
			NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
		end if
		IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
			NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
		end if
	END DO
	
	do K=1,nnd
		VEXT(K,1:2)=NODES_LIST(k,1:2)
	END DO
	      
end subroutine coordinates_face_PERIOD2d1_MovingMesh





SUBROUTINE CALCULATE_FLUXESHI_MovingMesh_2D(N)
    !> @brief
    !> This subroutine computes the fluxes for linear-advection equation in 2D
    IMPLICIT NONE
    INTEGER,INTENT(IN)::N
    REAL::GODFLUX2,sum_detect,lamxl,lamyl
    INTEGER::I,L,K,NGP,KMAXE,IQP, NEIGHBOR_INDEX, NEIGHBOR_FACE_INDEX
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_TEMP,WEIGHTS_DG !Quadrature weights for interfaces
    REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    REAL::ANGLE1,ANGLE2,NX,NY,NZ,NORMALVECT, normal_qp_velocity, net_normal_velocity
    INTEGER::facex,POINTX,ICONSIDERED
    REAL,DIMENSION(1:NOF_VARIABLES)::CLEFT,CRIGHT,HLLCFLUX,RHLLCFLUX
    REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG
    real,dimension(1:dimensiona)::vleft,vright, qp_velocity
  
    IF (DG.EQ.1)THEN
        print *, "Discontinous Galerkin not supported in moving mesh mode"
        call abort
        allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
    END IF
    
    KMAXE = XMPIELRANK(N)
    
    CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
    WEIGHTS_TEMP(1:QP_LINE_N) = WEQUA2D(1:QP_LINE_N)
      
    !$OMP DO
    DO I=1,KMAXE
    
        if (initcond.eq.3)then
            lamxl =-ielem(n,i)%yyc + 0.5d0
            lamyl = ielem(n,i)%xxc - 0.5d0
        else
            lamxl=lamx
            lamyl=lamy
        end if
        
        if (dg.eq.1)then
            RHS(I)%VALDG = ZERO
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
        else
            RHS(I)%VAL=ZERO
        end if
        
        ICONSIDERED=I
        
        IF (DG.EQ.1) THEN
            DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
        END IF
          
        IF (IELEM(N,I)%INTERIOR.EQ.0) THEN ! Element is interior
  
            DO L=1,IELEM(N,I)%IFCA
                GODFLUX2=ZERO
                NX=IELEM(N,I)%FACEANGLEX(L)
                NY=IELEM(N,I)%FACEANGLEY(L)
                facex=l
                
                NORMALVECT=(NX*LAMXl)+(NY*LAMYl)

                IQP=QP_LINE_N
                
                NEIGHBOR_INDEX = IELEM(N,I)%INEIGH(L)
                NEIGHBOR_FACE_INDEX = IELEM(N,I)%INEIGHN(L)
                  
                DO NGP=1,IQP
                    POINTX=NGP
  
                    IF (DG.EQ.1) THEN
                        CLEFT  = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                    ELSE !FV
                        CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
                        CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                        ! CLEFT(1)=U_C(I)%VAL(1,1)
                        ! CRIGHT(1)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1)
                    END IF
                    if (CLEFT(1).ne.CLEFT(1)) print*, "CLEFT(1) NaN"
                    if (CRIGHT(1).ne.CRIGHT(1)) print*, "CRIGHT(1) NaN"
                    vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:2)
                    ! vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(IELEM(N,I)%INEIGHN(L),NGP,1:2)
                    if ((vleft(1) .ne.vleft(1) ).or.(vleft(2) .ne.vleft(2) )) print*, "vleft NaN"
                    !if ((vright(1).ne.vright(1)).or.(vright(2).ne.vright(2))) print*, "vright NaN"
                    ! qp_velocity = 0.5*(vleft + vright)
                    qp_velocity = vleft
                    normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                    net_normal_velocity = NORMALVECT - normal_qp_velocity
                    ! print*,"net normal velocity", net_normal_velocity
                    CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,net_normal_velocity,HLLCFLUX)
                    if (HLLCFLUX(1).ne.HLLCFLUX(1)) print*, "NaN HLLC flux"
                    IF (DG.EQ.1) THEN
                        ! Riemann flux at interface quadrature points times basis
                        RHLLCFLUX(1)=HLLCFLUX(1)
                        DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
                    ELSE !FV
                        GODFLUX2 = GODFLUX2 + (HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                    END IF
                END DO
                  
                ! if (RHS(I)%VAL(1).ne.RHS(I)%VAL(1)) print*, "NaN flux"
                IF (DG /= 1) RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2
  
            END DO
  
        ELSE IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
  
            DO L=1,IELEM(N,I)%IFCA
                FACEX = L
                NX=IELEM(N,I)%FACEANGLEX(L)
                NY=IELEM(N,I)%FACEANGLEY(L)
                NORMALVECT=(NX*LAMXl)+(NY*LAMYl)
                if (nx.ne.nx) print *, "nx NaN at an interface"
                if (ny.ne.ny) print *, "ny NaN at an interface"
                if (NORMALVECT.ne.NORMALVECT) print *, "NORMALVECT NaN at an interface"
                IQP=QP_LINE_N
                
                GODFLUX2=ZERO
                DO NGP=1,IQP
                    POINTX = NGP

                    IF (DG == 1) THEN
                        CLEFT = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                    ELSE
                        CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
                        ! CLEFT(1)=U_C(I)%VAL(1,1)
                    END IF
                    
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                IF (DG == 1) THEN
                                    CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE !FV
                                    CRIGHT(1) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                                    ! CRIGHT(1)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1)
                                END IF
                            ELSE !NOT PERIODIC ONES IN MY CPU
                                CRIGHT(1:nof_variables)=CLEFT(1:nof_variables)
                            END IF
                        ELSE
                            IF (DG == 1) THEN
                                CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                            ELSE !FV
                                CRIGHT(1) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                                ! CRIGHT(1)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1)
                            END IF
                        END IF
                        
                    ELSE !IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                            
                                IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                else
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                end if
                            END IF
                        ELSE
                            IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                            ELSE
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                            end if
                        END IF
                    END IF

                    if (CLEFT(1).ne.CLEFT(1)) print*, "CLEFT(1) NaN at an interface"
                    if (CRIGHT(1).ne.CRIGHT(1)) print*, "CRIGHT(1) NaN at an interface"
                    vleft  = ILOCAL_RECON3(I)%QPOINTS_velocity(L,NGP,1:dimensiona)
                    if ((vleft(1) .ne.vleft(1) ).or.(vleft(2) .ne.vleft(2) )) print*, "vleft NaN at an interface"
                    ! vright = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%QPOINTS_velocity(L,NGP,1:dimensiona)
                    ! qp_velocity = 0.5*(vleft+vright)
                    qp_velocity = vleft
                    normal_qp_velocity = (nx * qp_velocity(1)) + (ny * qp_velocity(2))
                    CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT - normal_qp_velocity,HLLCFLUX)
                    if (HLLCFLUX(1).ne.HLLCFLUX(1)) print*, "NaN HLLC flux at an interface"
                    IF (DG.EQ.1) THEN
                        !Riemann flux at interface quadrature points times basis
                        RHLLCFLUX(1)=HLLCFLUX(1)
                        DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX) 
                    ELSE !FV
                        GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                    END IF
    
                END DO

                ! if (RHS(I)%VAL(1).ne.RHS(I)%VAL(1)) print*, "NaN flux at an interface"
                IF (DG /= 1) RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2

            END DO
  
        END IF
          
        IF (DG == 1) DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
        IF (DG == 1) RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS
          
    END DO
    !$OMP END DO 
  
    IF (DG.EQ.1)THEN
        deallocate(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
    END IF
  
END SUBROUTINE CALCULATE_FLUXESHI_MovingMesh_2D




	

  





end module