MODULE MESHMOVEMENT_module
    USE LIBRARY
    USE MPIINFO
    USE OMP_LIB
    USE DECLARATION
    USE FLUXES
    USE TRANSFORM
    USE TRANSFORM_MovingMesh
    USE PRESTORE_MOVINGMESH
    USE MESHQULITY_module
    USE COMMUNICATIONS
    USE RECON

    IMPLICIT NONE

    real::max_gradient_magnitude

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





subroutine clamp(x, bottom, top)
    implicit none
    real, intent(in)::bottom, top
    real, intent(inout)::x

    if (top.lt.bottom) then
        print*,"invalid clamp", bottom, ">", top
    end if
    if (x.lt.bottom) then
        x = bottom
    end if
    if (x.gt.top) then
        x = top
    end if

end subroutine clamp





function trinagle_area(a, b, c)
    implicit none
    real,dimension(1:dimensiona)::a, b, c
    real,dimension(1:3)::v1, v2, cross
    real::trinagle_area
    integer::i

    v1 = zero
    v2 = Zero

    do i = 1, dimensiona
        v1(i) = b(i) - a(i)
        v2(i) = c(i) - b(i)
    end do

    cross(1) = (v1(2)*v2(3)) - (v1(3)*v2(2))
    cross(2) = (v1(3)*v2(1)) - (v1(1)*v2(3))
    cross(3) = (v1(1)*v2(2)) - (v1(2)*v2(1))

    trinagle_area = (cross(1)*cross(1)) + (cross(2)*cross(2)) + (cross(3)*cross(3))

    trinagle_area = sqrt(trinagle_area / 4.0)

end function trinagle_area





function trinagle_area_2D(a, b, c)
    implicit none
    real,dimension(1:dimensiona)::a, b, c
    real::trinagle_area_2D
    real,dimension(1:dimensiona)::v1, v2
    
    real::cross
    integer::i

    v1(:) = b(:) - a(:)
    v2(:) = c(:) - b(:)

    cross = (v1(1)*v2(2)) - (v1(2)*v2(1))

    trinagle_area_2D = cross*0.5

end function trinagle_area_2D





subroutine polygon_centre(vert_num, vertices, centre)
    ! this function assumes that the polygon is convex
    implicit none
    integer::vert_num
    real,intent(in),dimension(1:vert_num,1:dimensiona)::vertices
    real,intent(out),dimension(1:dimensiona)::centre
    integer::i

    if (vert_num.le.0) then
        print *, "trying to find a centre of 0-gon"
        call abort()
    end if
    
    centre(:) = zero
    do i = 1, vert_num
        centre(:) = centre(:) + (vertices(i,:)/real(vert_num))
    end do


end subroutine polygon_centre





subroutine pseudoVoronoi_centre(vert_num, vertices, centre)
    implicit none
    integer::vert_num
    real,intent(in),dimension(1:vert_num,1:dimensiona)::vertices
    real,intent(out),dimension(1:dimensiona)::centre
    real,dimension(1:dimensiona)::temp_centre, helper, v1, v2
    real,dimension(1:vert_num)::l, r
    real,dimension(1:vert_num, 1:dimensiona)::J
    ! real,dimension(1:vert_num, 1:dimensiona)::JtJinvJtt
    ! real,dimension(1:dimensiona, 1:vert_num)::Jt
    real,dimension(1:dimensiona, 1:vert_num)::JtJinvJt
    real,dimension(1:dimensiona, 1:dimensiona)::JtJ
    real,dimension(1:dimensiona, 1:dimensiona)::JtJinv
    real,dimension(1:dimensiona)::dl_av
    real,dimension(1:dimensiona)::delta
    real::l_av, det
    real::coord_diff
    integer::a, b, c, step

    if (dimensiona.ne.2) Then
        print*, "pseudo-Voronoi mesh is not supported in 3D"
        call abort()
    end if

    if (vert_num.le.0) then
        print*, "trying to find a centre of 0-gon"
        call abort()
    end if
    
    centre(:) = zero
    do a = 1, vert_num
        centre(:) = centre(:) + (vertices(a,:)/real(vert_num))
    end do

    do step = 1, vert_num - 1
        l_av = 0
        dl_av(:) = zero
        do a = 1, vert_num
            l(a) = zero
            do b = 1, dimensiona
                coord_diff = centre(b) - vertices(a, b)
                l(a) = l(a) + (coord_diff*coord_diff)
            end do
            l(a) = sqrt(l(a))
            l_av = l_av + (l(a)/real(vert_num))

            do b = 1, dimensiona
                J(a, b) = (centre(b) - vertices(a, b))/l(a)
                dl_av(b) = dl_av(b) + (J(a, b)/real(vert_num))
            end do
        end do
        
        do a = 1, vert_num
            r(a) = l(a) - l_av
            J(a, :) = J(a, :) - dl_av(:)
            ! do b = 1, dimensiona
                ! J(a, b) = J(a, b) - dl_av(b)
                ! Jt(b, a) = J(a, b)
            ! end do
        end do

        JtJ(:,:) = zero
        do a = 1, dimensiona
            do b = 1, dimensiona
                do c = 1, vert_num
                    ! JtJ(a, b) = JtJ(a, b) + (Jt(a, c) * J(c, b))
                    JtJ(a, b) = JtJ(a, b) + (J(c, a) * J(c, b))
                end do
            end do
        end do

        det = (JtJ(1,1)*JtJ(2,2)) - (JtJ(1,2)*JtJ(2,1))

        if (abs(det).lt.0.000000001) then
            goto 10001
        end if

        JtJinv(1,1) = JtJ(2,2)/det
        JtJinv(2,2) = JtJ(1,1)/det
        JtJinv(1,2) = (-1.0)*JtJ(1,2)/det
        JtJinv(2,1) = (-1.0)*JtJ(2,1)/det

        JtJinvJt(:,:) = zero 
        ! JtJinvJtt(:,:) = zero 
        do a = 1, dimensiona
            do b = 1, vert_num
                do c = 1, dimensiona
                    ! JtJinvJt(a,b) = JtJinvJt(a,b) + (JtJinv(a,c)*Jt(c,b))
                    JtJinvJt(a,b) = JtJinvJt(a,b) + (JtJinv(a,c)*J(b,c))
                    ! JtJinvJtt(b,a) = JtJinvJt(b,a) + (JtJinv(a,c)*J(b,c))
                end do
            end do
        end do

        delta(:) = zero
        do a = 1, dimensiona
            do b = 1, vert_num
                delta(a) = delta(a) + (JtJinvJt(a,b)*r(b))
            end do
        end do 
        
        centre(:) = centre(:) - delta(:)
    end do

    10001 continue

end subroutine pseudoVoronoi_centre





subroutine cell_centre(vert_num, vertices, centre)
    implicit none
    integer::vert_num
    real,intent(in),dimension(1:dimensiona, 1:vert_num)::vertices
    real,intent(out),dimension(1:dimensiona)::centre
    real,dimension(1:dimensiona)::temp_centre, helper, v1, v2
    real::area, area_sum
    integer::i, j

    if (.not.((vert_num.eq.3).or.(vert_num.eq.4))) then
        print *, "invalid cell size"
        call abort()
    end if
    
    temp_centre = zero
    do i = 1,vert_num
        temp_centre(:) = temp_centre(:) + (vertices(:,i)/real(vert_num))
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

            v1 = vertices(:,i)
            v2 = vertices(:,j)

            helper(:) = ((v1(:) + v2(:) + temp_centre(:)) / 3.0)
            area = trinagle_area_2D(v1(:), v2(:), temp_centre(:))

            area_sum = area_sum + area
            centre(:) = centre(:) + (helper(:) * area)
        end do

        centre(:) = centre(:) / area_sum
    end if

end subroutine cell_centre





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
    integer,intent(in)::index_from, index_to
    integer::node_index

    !$omp do
    do node_index=1, kmaxn 
        local_nodes(node_index)%positions(index_to,1:dimensiona) = local_nodes(node_index)%positions(index_from,1:dimensiona)
    end do
    !$omp end do

end subroutine





subroutine COPY_BACK_LOCAL_NODES_and_Volumes(N, index_from, index_to)
    implicit NONE
    integer,intent(in)::N, index_from, index_to
    integer::node_index, cell_index, kmaxe
    kmaxe = xmpielrank(N)

    !$omp do
    do node_index=1, kmaxn 
        local_nodes(node_index)%positions(index_to,1:dimensiona) = local_nodes(node_index)%positions(index_from,1:dimensiona)
    end do
    !$omp end do

    !$omp do
    do cell_index=1, kmaxe
        ielem(N, cell_index)%moving_volume(index_to) = ielem(N, cell_index)%moving_volume(index_from)
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
    integer::index1, index2, index3
    integer::cell_index, node_index, node_index_1, node_index_2, edge_index, cpu_index, neighbour_cpu
    integer::num_requests, counter
    integer::kmaxe
    ! integer,dimension(isize)::num_nodes_to_send
    integer,dimension(0:isize-1)::num_interface_nodes
    integer,dimension(2*2*isize)::requests
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
            local_nodes(node_index)%communication        = 0 
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
            ! do edge_index = 1,IELEM(N,cell_index)%nonodes
            !     node_index = IELEM(N,cell_index)%nodes(edge_index)
                
            !     if (edge_index.lt.IELEM(N,cell_index)%nonodes) then
            !         node_index_2 = (IELEM(N,cell_index)%nodes(edge_index + 1))
            !     else
            !         node_index_2 = (IELEM(N,cell_index)%nodes(1))
            !     end if
            !     if (IELEM(N,cell_index)%INEIGHB(edge_index).ne.N) then
            !         ! local_nodes(node_index)%num_cpus = 2
            !         local_nodes(node_index)%boundary = max(local_nodes(node_index)%boundary, 1)
            !         ! local_nodes(node_index_2)%num_cpus = 2
            !         local_nodes(node_index_2)%boundary = max(local_nodes(node_index_2)%boundary, 1)
            !     end if
            !     if (ielem(n,cell_index)%ibounds(edge_index).gt.0) then
            !         if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.5) then
            !             local_nodes(node_index)%boundary = 2
            !             local_nodes(node_index_2)%boundary = 2
            !         end if
            !     end if
            ! end do
            do edge_index = 1,IELEM(N,cell_index)%ifca
                node_index_1 = IELEM(N,cell_index)%nodes_faces(edge_index, 1)
                node_index_2 = IELEM(N,cell_index)%nodes_faces(edge_index, 2)
                if (IELEM(N,cell_index)%INEIGHB(edge_index).ne.N) then
                    ! local_nodes(node_index_1)%num_cpus = 2
                    ! local_nodes(node_index_1)%boundary = max(local_nodes(node_index_1)%boundary, 1)
                    ! local_nodes(node_index_2)%num_cpus = 2
                    ! local_nodes(node_index_2)%boundary = max(local_nodes(node_index_2)%boundary, 1)
                    local_nodes(node_index_1)%communication = 1
                    local_nodes(node_index_2)%communication = 1
                end if
                if (ielem(n,cell_index)%ibounds(edge_index).gt.0) then
                    if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.5) then
                        local_nodes(node_index_1)%communication = 1
                        local_nodes(node_index_2)%communication = 1
                    else
                        local_nodes(node_index_1)%boundary = 1
                        local_nodes(node_index_2)%boundary = 1
                    end if
                    if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.gt.100) then
                        local_nodes(node_index_1)%boundary = ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode
                        local_nodes(node_index_2)%boundary = ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode
                    end if
                end if
            end do
        end do
  
        my_num_interface_nodes = 0
        my_num_boundary_nodes = 0
        my_num_moving_nodes = 0
        do node_index=1,kmaxn 
            if ((local_nodes(node_index)%communication.gt.0).or.(local_nodes(node_index)%num_local_neighbours.eq.0)) then
                my_num_interface_nodes = my_num_interface_nodes + 1
            end if
            if (local_nodes(node_index)%boundary.gt.0) then
                my_num_boundary_nodes = my_num_boundary_nodes + 1
            end if
            if (local_nodes(node_index)%boundary.gt.100) then
                my_num_moving_nodes = my_num_moving_nodes + 1
            end if
        end do
        allocate(local_interface_nodes(my_num_interface_nodes))
        allocate(local_boundary_nodes(my_num_boundary_nodes))
        allocate(local_moving_nodes(my_num_moving_nodes))
        ! print*,N,"my_num_interface_nodes =", my_num_interface_nodes, "\n", N, "my_num_boundary_nodes =", my_num_boundary_nodes,"\n", N, "my_num_moving_nodes =", my_num_moving_nodes
        index1 = 0
        index2 = 0
        index3 = 0
        do node_index=1,kmaxn
            if ((local_nodes(node_index)%communication.gt.0).or.(local_nodes(node_index)%num_local_neighbours.eq.0)) then
                index1 = index1 + 1
                local_interface_nodes(index1) = node_index
            end if
            if (local_nodes(node_index)%boundary.gt.0) then
                index2 = index2 + 1
                local_boundary_nodes(index2) = node_index
            end if
            if (local_nodes(node_index)%boundary.gt.100) then
                index3 = index3 + 1
                local_moving_nodes(index3) = node_index
            end if
        end do
        if (index1.ne.my_num_interface_nodes) then
            print *, "missmatch in the number of interface nodes on CPU", N
        ! else
        !     print *, my_num_interface_nodes, "= num interface nodes on CPU", N
        end if
        if (index2.ne.my_num_boundary_nodes) then
            print *, "missmatch in the number of boundary nodes on CPU", N
        ! else
        !     print *, my_num_boundary_nodes, "= num boundary nodes on CPU", N
        end if
        if (index3.ne.my_num_moving_nodes) then
            print *, "missmatch in the number of moving nodes on CPU", N
        ! else
        !     print *, my_num_moving_nodes, "= num moving nodes on CPU", N
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
                ! num_requests = num_requests + 1
                ! CALL MPI_ISENDRECV(buff(N)%coords(1:my_num_interface_nodes*dimensiona), my_num_interface_nodes*dimensiona, MPI_DOUBLE_PRECISION, cpu_index, 789,&
                !         buff(cpu_index)%coords(1:num_interface_nodes(cpu_index)*dimensiona), num_interface_nodes(cpu_index)*dimensiona, MPI_DOUBLE_PRECISION, cpu_index, 789,&
                !         MPI_COMM_WORLD, requests(num_requests), IERROR)
                num_requests = num_requests + 1
                CALL MPI_ISEND(buff(N)%coords(1:my_num_interface_nodes*dimensiona), my_num_interface_nodes*dimensiona, MPI_DOUBLE_PRECISION, cpu_index, 789,&
                        MPI_COMM_WORLD, requests(num_requests), IERROR)
                num_requests = num_requests + 1
                CALL MPI_IRECV(buff(cpu_index)%coords(1:num_interface_nodes(cpu_index)*dimensiona), num_interface_nodes(cpu_index)*dimensiona, MPI_DOUBLE_PRECISION, cpu_index, 789,&
                        MPI_COMM_WORLD, requests(num_requests), IERROR)

                ! num_requests = num_requests + 1
                ! CALL MPI_ISENDRECV(buff(N)%count(1:my_num_interface_nodes), my_num_interface_nodes, MPI_INTEGER, cpu_index, 987,&
                !         buff(cpu_index)%count(1:num_interface_nodes(cpu_index)), num_interface_nodes(cpu_index), MPI_INTEGER, cpu_index, 987,&
                !         MPI_COMM_WORLD, requests(num_requests), IERROR)
                num_requests = num_requests + 1
                CALL MPI_ISEND(buff(N)%count(1:my_num_interface_nodes), my_num_interface_nodes, MPI_INTEGER, cpu_index, 987,&
                        MPI_COMM_WORLD, requests(num_requests), IERROR)
                num_requests = num_requests + 1
                CALL MPI_IRECV(buff(cpu_index)%count(1:num_interface_nodes(cpu_index)), num_interface_nodes(cpu_index), MPI_INTEGER, cpu_index, 987,&
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
                            if (local_nodes(node_index)%num_neighbours.gt.max_num_node_neighbours) then
                                max_num_node_neighbours = local_nodes(node_index)%num_neighbours
                            end if
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


subroutine find_node_lagrangian_velocity(stage, position_index, N)
    implicit none
    integer,intent(in)::stage, position_index, N

    select case(node_solver_type)
      case (1)
        call FirstOrderNodeAverage(stage, N)
      case(2)
        call FirstOrderNodeMassWeightedAverage(stage, position_index, N)
      case(3)
        call HighOrderNodeAverage(stage, position_index, N)
      case(4)
        call HighOrderNodeMassWeightedAverage(stage, position_index, N)
      case(5)
        call HighOrderUpstreamNodeAverage(stage, position_index, N)
      case default
        print*, "invalid node solver"
    end select

end subroutine find_node_lagrangian_velocity





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
        select case (moving_mesh_mode)
          case(1, 2)
            call find_node_lagrangian_velocity(1, position_index, N)
          case(3)
            call FirstOrderNodeAverage_withVF(1, N)
          case(4)
            call directReALE_node_velocity(1, position_index, d_t, N)
          case(5, 7)
            call find_node_lagrangian_velocity(1, position_index, N)
            call find_node_relaxation_velocity(1, position_index, d_t, N)
          case(6)
            call FirstOrderNodeAverage_withVf(1, N)
            call find_node_relaxation_velocity(1, position_index, d_t, N)
          case(8)
            call find_node_lagrangian_velocity(1, position_index, N)
            ! call find_node_relaxation_velocity_and_normalized_density_gradient(1, position_index, d_t, N)
            call find_node_normalized_density_gradient(1, position_index, d_t, N)
            call find_node_relaxation_velocity(1, position_index, d_t, N)
          case(9, 11)
            call find_node_lagrangian_velocity(1, position_index, N)
            call find_node_normalized_density_gradient(1, position_index, d_t, N)
            if (relaxation_centre_type.le.3) then
                call find_node_relaxation_velocity(1, position_index, d_t, N)
            else if (relaxation_centre_type.eq.4) then
                call find_node_centre_relaxation_velocity(1, position_index, d_t, N)
            else if ((relaxation_centre_type.eq.5).or.(relaxation_centre_type.eq.7)) then
                call find_node_Jacobi_relaxation_velocity(0, position_index, d_t, N)
            else
                print*,"invalid node relaxation algorithm"
            end if

          case(10, 13)
            call find_node_lagrangian_velocity(1, position_index, N)
            call find_node_normalized_density_gradient(1, position_index, d_t, N)
            if (moving_mesh_mode.eq.13) then
                call find_node_normalized_vf_gradient(1, position_index, d_t, N)
            end if
            !$omp barrier
            call enforce_node_lagrangian_velocity_BC(position_index, d_t, N)
            !$omp barrier
            if (relaxation_centre_type.le.3) then
                call find_moved_node_relaxation_velocity(1, position_index, d_t, N)
            else if (relaxation_centre_type.eq.4) then
                call find_moved_node_centre_relaxation_velocity(1, position_index, d_t, N)
            else if ((relaxation_centre_type.eq.5).or.(relaxation_centre_type.eq.7)) then
                call find_node_Jacobi_relaxation_velocity(1, position_index, d_t, N)
            else
                print*,"invalid node relaxation algorithm"
            end if

          case(14)
            call find_node_lagrangian_velocity(1, position_index, N)
            call find_node_normalized_density_gradient(1, position_index, d_t, N)
            if (governingequations.eq.-1) then
                call find_node_normalized_vf_gradient(1, position_index, d_t, N)
            end if
            !$omp barrier
            call enforce_node_lagrangian_velocity_BC(position_index, d_t, N)
            !$omp barrier
            if (relaxation_centre_type.le.3) then
                call find_moved_node_relaxation_velocity(1, position_index, d_t, N)
            else if (relaxation_centre_type.eq.4) then
                call find_moved_node_centre_relaxation_velocity(1, position_index, d_t, N)
            else if ((relaxation_centre_type.eq.5).or.(relaxation_centre_type.eq.7)) then
                call find_node_Jacobi_relaxation_velocity(1, position_index, d_t, N)
            else
                print*,"invalid node relaxation algorithm"
            end if

          case(12)
            call find_node_lagrangian_velocity(1, position_index, N)
            call find_node_volume_ratio(1, position_index, N)
            call find_node_normalized_density_gradient(1, position_index, d_t, N)
            call find_node_relaxation_velocity(1, position_index, d_t, N)

         case DEFAULT
            print *, "invalid moving mesh mode"
            call abort()
        end select
        !$omp barrier
        call CombineNodeVelocities(1, position_index, d_t, N)
    endif

    !$omp barrier

    call enforce_node_velocity_BC(position_index, d_t, N)

    !$omp barrier

    ! call clamp_node_velocity_to_CFL(position_index, d_t, N)
    call clamp_node_velocity_to_fraction(position_index, d_t, N)

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
!             local_nodes(node)%velocity(d) = local_nodes(node)%velocity(d) * mesh_velocity_multiple
            
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
    ! integer:: M
    real::rho,v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, gammal
    integer::rho_index

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.dimensiona) then
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
                    ! index = index +1
                    if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                        print *, "copying too much data to send buffer from", N, "to", cpu 
                    end if
                    node_snd_buffer(cpu)%data(index + k) = copy(rho_index + k)
                end do
                index = index + num_values_to_send_per_node
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

        local_nodes(node_index)%lagrangian_velocity(:) = 0.0

        do i = 1, local_nodes(node_index)%num_local_neighbours
            cell_index = local_nodes(node_index)%local_neighbours(i)
            copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
            call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
            do j = 1, dimensiona
                v = copy(rho_index+j)
                local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + v
            end do
        end do

        do i = 1, local_nodes(node_index)%num_cpus
            cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
            index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
            do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                do k = 1, dimensiona
                    ! index = index+1
                    if ((index+k).gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
                        print *, "copying too much data from receive buffer from", N, "to", cpu 
                    end if
                    v = node_rcv_buffer(cpu_index)%data(index+k)
                    local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + v
                end do
                index = index + num_values_to_send_per_node
            end do
        end do

        do j = 1, dimensiona
            if (local_nodes(node_index)%num_neighbours.eq.0) then
                print*,"dividing by zero num_neighbours in node", node_index, "CPU", n
            end if
            local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) / real(local_nodes(node_index)%num_neighbours)
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

    if (num_values_to_send_per_node.lt.(dimensiona+1)) then
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
                        if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index+k) = copy(rho_index + k)
                    end do
                    if ((index+dimensiona+1).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                        print *, "copying too much data to send buffer from", N, "to", cpu 
                    end if
                    node_snd_buffer(cpu)%data(index+dimensiona+1) = copy(vf_index)
                    index = index + num_values_to_send_per_node
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

            local_nodes(node_index)%lagrangian_velocity(:) = 0.0
            vf_closest_to_half = zero

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)
                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                do j = 1, dimensiona
                    v = copy(rho_index+j)
                    local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + v
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
                        local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + v
                    end do
                    index = index+1
                    vf = node_rcv_buffer(cpu_index)%data(index)
                    if (abs(vf-0.5).lt.abs(vf_closest_to_half-0.5)) then
                        vf_closest_to_half = vf
                    end if
                end do
            end do

            do j = 1, dimensiona
                if (local_nodes(node_index)%num_neighbours.eq.0) then
                    print*,"dividing by zero num_neighbours in node", node_index, "CPU", n
                end if
                local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) / real(local_nodes(node_index)%num_neighbours)

                if ((vf_closest_to_half.ge.0.0).and.(vf_closest_to_half.le.1.0)) then
                    if (lagrangian_mesh_velocity_multiple_function_type.eq.1) then
                        if (vf_closest_to_half.le.0.5) then
                            local_mesh_velocity_multiple = 2.0*vf_closest_to_half
                        else
                            local_mesh_velocity_multiple = 2.0*(1.0-vf_closest_to_half)
                        end if
                    else if (lagrangian_mesh_velocity_multiple_function_type.eq.2) then
                        local_mesh_velocity_multiple = sin(pi*vf_closest_to_half)
                    else if (lagrangian_mesh_velocity_multiple_function_type.eq.3) then
                        if ((vf_closest_to_half.lt.0.125).or.(vf_closest_to_half.gt.0.875)) then
                            local_mesh_velocity_multiple = zero
                        else if (vf_closest_to_half.lt.0.25) then
                            local_mesh_velocity_multiple = (vf_closest_to_half-0.125)*8.0
                        else if (vf_closest_to_half.gt.0.75) then
                            local_mesh_velocity_multiple = (0.875-vf_closest_to_half)*8.0
                        else
                            local_mesh_velocity_multiple = 1.0
                        end if
                    else
                        print*,"invalid lagrangian mesh velocity multiple function type"
                        call abort()
                    end if
                else
                    local_mesh_velocity_multiple = 0.0
                end if
                local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) * local_mesh_velocity_multiple
            end do

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE FirstOrderNodeAverage_withVF





SUBROUTINE FirstOrderNodeMassWeightedAverage(stage, node_position_index, N)
    implicit none
    integer,intent(in)::stage, node_position_index, N
    integer::i, j, k, iter, node_index, cell_index, cpu_index, cpu, index
    integer:: M
    real::rho, v, mass, mass_sum
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, gammal
    integer::rho_index

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.(dimensiona+1)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    rho_index = 1

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (dimensiona+1)
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                    mass = copy(rho_index) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)

                    do k = 1,dimensiona
                        if ((index+k).gt.node_snd_count(cpu) * (dimensiona+1)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index + k) = copy(rho_index + k)
                    end do
                    if ((index + dimensiona+1).gt.node_snd_count(cpu) * (dimensiona+1)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu 
                    end if
                    node_snd_buffer(cpu)%data(index + dimensiona+1) = mass

                    index = index + (dimensiona+1)
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
                ! if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                !     print *,"send receive count missmatch on CPU", n
                !     call abort
                ! end if
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

            local_nodes(node_index)%lagrangian_velocity(:) = 0.0
            mass_sum = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)
                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                mass = copy(rho_index) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)
                do j = 1, dimensiona
                    v = copy(rho_index+j)
                    local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + (v*mass)
                end do
                mass_sum = mass_sum + mass
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*(dimensiona+1)
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    ! if ((index+dimensiona+1).gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                    !     print *, "copying too much data from receive buffer from", N, "to", cpu 
                    ! end if
                    mass = node_rcv_buffer(cpu_index)%data(index + dimensiona+1)
                    do k = 1, dimensiona
                        ! if ((index+k).gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                        !     print *, "copying too much data from receive buffer from", N, "to", cpu 
                        ! end if
                        v = node_rcv_buffer(cpu_index)%data(index+k)
                        local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + (v*mass)
                    end do
                    mass_sum = mass_sum + mass

                    index = index + dimensiona+1
                end do
            end do

            do j = 1, dimensiona
                local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) / mass_sum
            end do

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE FirstOrderNodeMassWeightedAverage





SUBROUTINE HighOrderNodeAverage(stage, node_position_index, N)
    implicit none
    integer,intent(in)::stage, node_position_index, N
    integer::i, j, k, iter, node_index, face_index, face_node_index, cell_index, cpu_index, cpu, index, read_index, cell_num_faces
    integer::cell_num_nodes
    integer:: M
    real::rho, v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, gammal
    integer::rho_index
    real::delta
    integer::first_time

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.(dimensiona)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    rho_index = 1

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * dimensiona
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)

                    cell_num_faces = ielem(n,cell_index)%ifca
                    first_time = 0
                    do face_index = 1, cell_num_faces
                        do face_node_index = 1, 2
                            if (ielem(n,cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                                if (first_time.eq.0) then
                                    copy(:) = ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index)
                                    do k = 1, nof_variables
                                        if (copy(k).ne.copy(k)) then
                                            print*,"NaN node reconstructed value before send in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    first_time = first_time + 1
                                else
                                    do k = 1, nof_variables
                                        delta = abs(copy(k) - ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index))
                                        if (delta.gt.(0.001*abs(copy(k)))) then
                                            print*,"suspeciously large change", copy(k), "vs", ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index), "in reconstruced value", k, "before send in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    copy(:) = (copy(:) + ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index))*0.5
                                    do k = 1, nof_variables
                                        if (copy(k).ne.copy(k)) then
                                            print*,"NaN node reconstructed value before send in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    first_time = first_time + 1
                                end if
                            end if
                        end do
                    end do
                    if (first_time.lt.2) print*,"not enough reconstructed nodes"
                    if (first_time.gt.2) print*,"too many reconstructed nodes"

                    call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                    do k = 1, nof_variables
                        if (copy(k).ne.copy(k)) then
                            print*,"NaN node reconstructed value before send after cons2prim in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                        end if
                    end do

                    do k = 1,dimensiona
                        index = index+1
                        if (index.gt.node_snd_count(cpu) * dimensiona) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "HighOrderNodeAverage"
                        end if
                        node_snd_buffer(cpu)%data(index) = copy(rho_index + k)
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
                    count = node_snd_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                ! if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                !     print *,"send receive count missmatch on CPU", n
                !     call abort
                ! end if
                count = node_snd_count(cpu_index) * dimensiona
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

            local_nodes(node_index)%lagrangian_velocity(:) = 0.0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)

                cell_num_faces = ielem(n,cell_index)%ifca
                first_time = 0
                do face_index = 1, cell_num_faces
                    do face_node_index = 1, 2
                        if (ielem(n,cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                            if (first_time.eq.0) then
                                copy(:) = ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index)
                                do k = 1, nof_variables
                                    if (copy(k).ne.copy(k)) then
                                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                first_time = first_time + 1
                            else
                                do k = 1, nof_variables
                                    delta = abs(copy(k) - ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index))
                                    if (delta.gt.(0.001*abs(copy(k)))) then
                                        print*,"suspeciously large change", copy(k), "vs", ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index), "in reconstruced value", k, "in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                copy(:) = (copy(:) + ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index))*0.5
                                do k = 1, nof_variables
                                    if (copy(k).ne.copy(k)) then
                                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                first_time = first_time + 1
                            end if
                        end if
                    end do
                end do
                if (first_time.lt.2) print*,"not enough reconstructed nodes"
                if (first_time.gt.2) print*,"too many reconstructed nodes"

                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                do k = 1, nof_variables
                    if (copy(k).ne.copy(k)) then
                        print*,"NaN node reconstructed value after cons2prim in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                    end if
                end do

                do j = 1, dimensiona
                    v = copy(rho_index+j)
                    local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + (v/real(local_nodes(node_index)%num_neighbours))
                end do
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*dimensiona
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    do k = 1, dimensiona
                        ! if ((index+k).gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                        !     print *, "copying too much data from receive buffer from", N, "to", cpu 
                        ! end if
                        v = node_rcv_buffer(cpu_index)%data(index+k)
                        local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + (v/real(local_nodes(node_index)%num_neighbours))
                    end do
                    index = index + dimensiona
                end do
            end do

            do j = 1, dimensiona
                if (local_nodes(node_index)%lagrangian_velocity(j).ne.local_nodes(node_index)%lagrangian_velocity(j)) then
                    print*,"NaN lagrangian velocity in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") on CPU", N
                end if
            end do

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE HighOrderNodeAverage





SUBROUTINE HighOrderNodeMassWeightedAverage(stage, node_position_index, N)
    implicit none
    integer,intent(in)::stage, node_position_index, N
    integer::i, j, k, iter, node_index, face_index, face_node_index, cell_index, cpu_index, cpu, index, read_index, cell_num_faces
    integer::cell_num_nodes
    ! integer:: M
    real::rho, v, mass, v_mass, mass_sum
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, gammal
    integer,parameter::rho_index=1
    real::delta
    integer::first_time

    integer,dimension(2*isize)::requests
    integer::num_requests, count

   ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.(dimensiona+1)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (dimensiona+1)
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)

                    cell_num_faces = ielem(n,cell_index)%ifca
                    first_time = 0
                    do face_index = 1, cell_num_faces
                        do face_node_index = 1, 2
                            if (ielem(n,cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                                if (first_time.eq.0) then
                                    copy(:) = ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index)
                                    do k = 1, nof_variables
                                        if (copy(k).ne.copy(k)) then
                                            print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    first_time = first_time + 1
                                else
                                    do k = 1, nof_variables
                                        delta = abs(copy(k) - ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index))
                                        if (delta.gt.(0.001*abs(copy(k)))) then
                                            print*,"suspeciously large change in reconstruced value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    copy(:) = (copy(:) + ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index))*0.5
                                    do k = 1, nof_variables
                                        if (copy(k).ne.copy(k)) then
                                            print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    first_time = first_time + 1
                                end if
                            end if
                        end do
                    end do
                    if (first_time.lt.2) print*,"not enough reconstructed nodes"
                    if (first_time.gt.2) print*,"too many reconstructed nodes"

                    call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                    do k = 1, nof_variables
                        if (copy(k).ne.copy(k)) then
                            print*,"NaN node reconstructed value after cons2prim in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                        end if
                    end do
                    ! mass = copy(rho_index) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)
                    mass = u_c(cell_index)%val(stage, rho_index) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)

                    do k = 1,dimensiona
                        index = index+1
                        if (index.gt.node_snd_count(cpu) * (dimensiona+1)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index) = copy(rho_index + k)
                    end do
                    index = index + 1
                    if (index.gt.node_snd_count(cpu) * (dimensiona+1)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "HighOrderNodeMassWeightedAverage"
                    end if
                    node_snd_buffer(cpu)%data(index) = mass
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
                    count = node_snd_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                ! if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                !     print *,"send receive count missmatch on CPU", n
                !     call abort
                ! end if
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

            local_nodes(node_index)%lagrangian_velocity(:) = 0.0
            mass_sum = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)

                cell_num_faces = ielem(n,cell_index)%ifca
                first_time = 0
                do face_index = 1, cell_num_faces
                    do face_node_index = 1, 2
                        if (ielem(n,cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                            if (first_time.eq.0) then
                                copy(:) = ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index)
                                do k = 1, nof_variables
                                    if (copy(k).ne.copy(k)) then
                                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                first_time = first_time + 1
                            else
                                do k = 1, nof_variables
                                    delta = abs(copy(k) - ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index))
                                    if (delta.gt.(0.001*abs(copy(k)))) then
                                        print*,"suspeciously large change in reconstruced value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                copy(:) = (copy(:) + ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index))*0.5
                                do k = 1, nof_variables
                                    if (copy(k).ne.copy(k)) then
                                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                first_time = first_time + 1
                            end if
                        end if
                    end do
                end do
                if (first_time.lt.2) print*,"not enough reconstructed nodes"
                if (first_time.gt.2) print*,"too many reconstructed nodes"

                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                do k = 1, nof_variables
                    if (copy(k).ne.copy(k)) then
                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                    end if
                end do
                ! mass = copy(rho_index) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)
                mass = u_c(cell_index)%val(stage, rho_index) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)
                do j = 1, dimensiona
                    v = copy(rho_index+j)
                    local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + (v*mass)
                    ! v_mass = copy(rho_index+j) * IELEM(N, cell_index)%moving_VOLUME(node_position_index)
                    ! local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + v_mass
                end do
                mass_sum = mass_sum + mass
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*(dimensiona+1)
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    ! if ((index+dimensiona+1).gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                    !     print *, "copying too much data from receive buffer from", N, "to", cpu 
                    ! end if
                    mass = node_rcv_buffer(cpu_index)%data(index + dimensiona+1)
                    do k = 1, dimensiona
                        ! if ((index+k).gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                        !     print *, "copying too much data from receive buffer from", N, "to", cpu 
                        ! end if
                        v = node_rcv_buffer(cpu_index)%data(index+k)
                        local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + (v*mass)
                        ! v_mass = node_rcv_buffer(cpu_index)%data(index+k)
                        ! local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + v_mass
                    end do
                    mass_sum = mass_sum + mass

                    index = index + dimensiona+1
                end do
            end do

            do j = 1, dimensiona
                local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) / mass_sum
                if (local_nodes(node_index)%lagrangian_velocity(j).ne.local_nodes(node_index)%lagrangian_velocity(j)) then
                    print*,"NaN lagrangian velocity in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") on CPU", N
                end if
            end do

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE HighOrderNodeMassWeightedAverage





SUBROUTINE HighOrderUpstreamNodeAverage(stage, node_position_index, N)
    implicit none
    integer,intent(in)::stage, node_position_index, N
    integer::i, j, k, iter, node_index, face_index, face_node_index, cell_index, cpu_index, cpu, index, read_index, cell_num_faces
    integer::cell_num_nodes
    ! integer:: M
    real::rho, dot
    real,dimension(1:nof_variables)::copy
    real,dimension(1:dimensiona)::v, CellCentre, NodeCentreVector
    real::dummuy_MP_PINFl, gammal
    integer::rho_index
    real::delta
    integer::first_time, num_contributions
    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.(dimensiona+1)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    rho_index = 1

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (dimensiona+1)
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)

                    cell_num_faces = ielem(n,cell_index)%ifca
                    first_time = 0
                    do face_index = 1, cell_num_faces
                        do face_node_index = 1, 2
                            if (ielem(n,cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                                if (first_time.eq.0) then
                                    copy(:) = ilocal_recon3(cell_index)%node_values(:, face_index, face_node_index)
                                    do k = 1, nof_variables
                                        if (copy(k).ne.copy(k)) then
                                            print*,"NaN node reconstructed valuein node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    first_time = first_time + 1
                                else
                                    ! do k = 1, nof_variables
                                    !     delta = abs(copy(k) - ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index))
                                    !     if (delta.gt.(0.01*abs(copy(k)))) then
                                    !         print*,"suspeciously large change in reconstruced value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    !     end if
                                    ! end do
                                    copy(:) = (copy(:) + ilocal_recon3(cell_index)%node_values(:, face_index, face_node_index))*0.5
                                    do k = 1, nof_variables
                                        if (copy(k).ne.copy(k)) then
                                            print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                        end if
                                    end do
                                    first_time = first_time + 1
                                end if
                            end if
                        end do
                    end do
                    if (first_time.lt.2) print*,"not enough reconstructed nodes"
                    if (first_time.gt.2) print*,"too many reconstructed nodes"

                    call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                    do k = 1, nof_variables
                        if (copy(k).ne.copy(k)) then
                            print*,"NaN node reconstructed value after cons2prim in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                        end if
                    end do

                    do k = 1,dimensiona
                        if ((index+k).gt.node_snd_count(cpu) * (dimensiona+1)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        v(k) = copy(rho_index + k)
                        node_snd_buffer(cpu)%data(index + k) = v(k)
                    end do

                    CellCentre(1) = ielem(n,cell_index)%xxc
                    CellCentre(2) = ielem(n,cell_index)%yyc
                    if (dimensiona.eq.3) then
                        CellCentre(3) = ielem(n,cell_index)%zzc
                    end if
                    NodeCentreVector(1:dimensiona) = CellCentre(1:dimensiona) - local_nodes(node_index)%positions(node_position_index, 1:dimensiona)
                    dot = dot_product(NodeCentreVector, v)

                    node_snd_buffer(cpu)%data(index + dimensiona+1) = dot

                    index = index + (dimensiona+1)
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
                    count = node_snd_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                ! if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                !     print *,"send receive count missmatch on CPU", n
                !     call abort
                ! end if
                count = node_snd_count(cpu_index) * (dimensiona+1)
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

            local_nodes(node_index)%lagrangian_velocity(:) = 0.0
            num_contributions = 0
            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)

                cell_num_faces = ielem(n,cell_index)%ifca
                first_time = 0
                do face_index = 1, cell_num_faces
                    do face_node_index = 1, 2
                        if (ielem(n,cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                            if (first_time.eq.0) then
                                copy(:) = ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index)
                                do k = 1, nof_variables
                                    if (copy(k).ne.copy(k)) then
                                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                first_time = first_time + 1
                            else
                                ! do k = 1, nof_variables
                                !     delta = abs(copy(k) - ilocal_recon3(cell_index)%node_values(k, face_index, face_node_index))
                                !     if (delta.gt.(0.01*abs(copy(k)))) then
                                !         print*,"suspeciously large change in reconstruced value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                !     end if
                                ! end do
                                copy(:) = (copy(:) + ilocal_recon3(cell_index)%node_values(1:nof_variables, face_index, face_node_index))*0.5
                                do k = 1, nof_variables
                                    if (copy(k).ne.copy(k)) then
                                        print*,"NaN node reconstructed value in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                                    end if
                                end do
                                first_time = first_time + 1
                            end if
                        end if
                    end do
                end do
                if (first_time.lt.2) print*,"not enough reconstructed nodes"
                if (first_time.gt.2) print*,"too many reconstructed nodes"

                call cons2prim(N, copy, dummuy_MP_PINFl, gammal)
                do k = 1, nof_variables
                    if (copy(k).ne.copy(k)) then
                        print*,"NaN node reconstructed value after cons2prim in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") in cell", cell_index, "on CPU", N
                    end if
                end do
                v(1:dimensiona) = copy(2:dimensiona+1)

                CellCentre(1) = ielem(n,cell_index)%xxc
                CellCentre(2) = ielem(n,cell_index)%yyc
                if (dimensiona.eq.3) then
                    CellCentre(3) = ielem(n,cell_index)%zzc
                end if
                NodeCentreVector(1:dimensiona) = CellCentre(1:dimensiona) - local_nodes(node_index)%positions(node_position_index, 1:dimensiona)
                dot = dot_product(NodeCentreVector, v)

                if (dot.ge.zero) then
                    num_contributions = num_contributions+1
                    do j = 1, dimensiona
                        local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) + copy(rho_index+j)
                    end do
                end if
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*(dimensiona+1)
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    dot = node_rcv_buffer(cpu_index)%data(index+(dimensiona+1))
                    if (dot.gt.zero) Then
                        num_contributions = num_contributions+1
                        do k = 1, dimensiona
                            ! if ((index+k).gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                            !     print *, "copying too much data from receive buffer from", N, "to", cpu 
                            ! end if
                            local_nodes(node_index)%lagrangian_velocity(k) = local_nodes(node_index)%lagrangian_velocity(k) + node_rcv_buffer(cpu_index)%data(index+k)
                        end do
                    end if
                    index = index + (dimensiona + 1)
                end do
            end do

            if (num_contributions.gt.0) then
                do j = 1, dimensiona
                    local_nodes(node_index)%lagrangian_velocity(j) = local_nodes(node_index)%lagrangian_velocity(j) / real(num_contributions)
                    if (local_nodes(node_index)%lagrangian_velocity(j).ne.local_nodes(node_index)%lagrangian_velocity(j)) then
                        print*,"NaN lagrangian velocity in node", node_index, "(", local_nodes(node_index)%positions(node_position_index,1), local_nodes(node_index)%positions(node_position_index,2), ") on CPU", N
                    end if
                end do
            end if

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE HighOrderUpstreamNodeAverage





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
            if (relaxation_centre_type.eq.1) then
                call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
            else if (relaxation_centre_type.eq.2) then
                call pseudoVoronoi_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
                if ((local_nodes(node_index)%velocity(1).le.xmin(n)).or.&
                        (local_nodes(node_index)%velocity(1).ge.xmax(n)).or.&
                        (local_nodes(node_index)%velocity(2).le.ymin(n)).or.&
                        (local_nodes(node_index)%velocity(2).ge.ymax(n))) then
                    call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
                end if
            else
                print*,"invalid centre algorithm in node relaxation velocity"
            end if
            local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (local_nodes(node_index)%relaxation_velocity(1:dimensiona) - local_nodes(node_index)%positions(position_index,1:dimensiona)) / d_t

            ! local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%relaxation_velocity(1:dimensiona) * mesh_velocity_multiple
        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE directReALE_node_velocity





SUBROUTINE find_node_relaxation_velocity(stage, position_index, d_t, N)
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
    ! real,dimension(1:max_num_node_neighbours)::weights


    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()
    
    rho_index = 1

    if (num_values_to_send_per_node.lt.dimensiona) then
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
                    ! copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    ! call cons2prim(N, copy, dummuy_MP_PINFl, dummy_gammal)
                    helper_centre_position(1) = ielem(N, cell_index)%xxc ! + (copy(rho_index+1)*d_t)
                    helper_centre_position(2) = ielem(N, cell_index)%yyc ! + (copy(rho_index+2)*d_t)
                    if (dimensiona.eq.3) then
                        helper_centre_position(3) = ielem(N, cell_index)%zzc ! + (copy(rho_index+3)*d_t)
                    end if
                    do k = 1,dimensiona
                        if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index+k) = helper_centre_position(k)
                    end do
                    index = index + num_values_to_send_per_node
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

            local_nodes(node_index)%relaxation_velocity(:) = 0.0
            counter = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1
                cell_index = local_nodes(node_index)%local_neighbours(i)

                helper_centre_position(1) = ielem(N, cell_index)%xxc
                helper_centre_position(2) = ielem(N, cell_index)%yyc
                if (dimensiona.eq.3) then
                    helper_centre_position(3) = ielem(N, cell_index)%zzc
                end if
                centre_positions(counter, :) = helper_centre_position(:)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    counter = counter+1
                    do k = 1, dimensiona
                        if ((index+k).gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index+k)
                    end do
                    index = index + num_values_to_send_per_node
                end do
            end do

            if (counter.ne.local_nodes(node_index)%num_neighbours) then
                print *, "something went wrong counter =/= local_nodes(node_index)%num_neighbours"
            end if
            if (relaxation_centre_type.eq.1) then
                call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
            else if (relaxation_centre_type.eq.2) then
                call pseudoVoronoi_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
                if ((local_nodes(node_index)%velocity(1).le.xmin(n)).or.&
                        (local_nodes(node_index)%velocity(1).ge.xmax(n)).or.&
                        (local_nodes(node_index)%velocity(2).le.ymin(n)).or.&
                        (local_nodes(node_index)%velocity(2).ge.ymax(n))) then
                    call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
                end if
            else
                print*,"invalid centre algorithm in node relaxation velocity"
            end if

            local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (local_nodes(node_index)%relaxation_velocity(1:dimensiona) - local_nodes(node_index)%positions(position_index,1:dimensiona)) / d_t
        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_node_relaxation_velocity





SUBROUTINE find_node_centre_relaxation_velocity(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, counter, node_index, face_index, cell_index, cpu_index, cpu, index, other_node_index, face_node_index, other_face_node_index
    ! integer:: M
    real::rho, v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, dummy_gammal
    real,dimension(1:dimensiona)::helper,centre
    ! real,dimension(1:max_num_node_neighbours)::weights

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.dimensiona) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * dimensiona
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    helper(:) = zero
                    do face_index = 1, ielem(n, cell_index)%ifca
                        do face_node_index = 1, 2
                            if (ielem(n, cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                                other_face_node_index = 3 - face_node_index ! swaps 1 with 2 and vice versa
                                other_node_index = ielem(n, cell_index)%nodes_faces(face_index, other_face_node_index)
                                helper(:) = helper(:) + local_nodes(other_node_index)%positions(position_index,:)
                            end if
                        end do
                    end do
                    helper = helper(:)/(2.0*real(local_nodes(node_index)%num_neighbours))
                    
                    do k = 1,dimensiona
                        index = index+1
                        if (index.gt.node_snd_count(cpu) * dimensiona) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index) = helper(k)
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
                    count = node_snd_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * dimensiona
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

            local_nodes(node_index)%relaxation_velocity(:) = 0.0
            counter = 0

            centre(:) = zero

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1
                cell_index = local_nodes(node_index)%local_neighbours(i)

                helper(:) = zero
                do face_index = 1, ielem(n, cell_index)%ifca
                    do face_node_index = 1, 2
                        if (ielem(n, cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                            other_face_node_index = 3 - face_node_index ! swaps 1 with 2 and vice versa
                            other_node_index = ielem(n, cell_index)%nodes_faces(face_index, other_face_node_index)
                            helper(:) = helper(:) + local_nodes(other_node_index)%positions(position_index,:)
                        end if
                    end do
                end do
                helper = helper(:)/(2.0*real(local_nodes(node_index)%num_neighbours))

                centre(:) = centre(:) + helper(:)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*dimensiona
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * dimensiona) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        helper(k) = node_rcv_buffer(cpu_index)%data(index)
                    end do
                    centre(:) = centre(:) + helper(:)
                end do
            end do

            local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (centre(1:dimensiona) - local_nodes(node_index)%positions(position_index,1:dimensiona)) / d_t
        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_node_centre_relaxation_velocity





SUBROUTINE find_moved_node_centre_relaxation_velocity(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, counter, node_index, face_index, cell_index, cpu_index, cpu, index, other_node_index, face_node_index, other_face_node_index
    ! integer:: M
    real::rho, v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, dummy_gammal
    real,dimension(1:dimensiona)::helper,centre,new_position
    ! real,dimension(1:max_num_node_neighbours)::weights

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.dimensiona) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * dimensiona
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    helper(:) = zero
                    do face_index = 1, ielem(n, cell_index)%ifca
                        do face_node_index = 1, 2
                            if (ielem(n, cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                                other_face_node_index = 3 - face_node_index ! swaps 1 with 2 and vice versa
                                other_node_index = ielem(n, cell_index)%nodes_faces(face_index, other_face_node_index)
                                new_position(1:dimensiona) = local_nodes(other_node_index)%positions(position_index,:) &
                                                            + (d_t*local_nodes(other_node_index)%velocity(1:dimensiona))
                                helper(:) = helper(:) + new_position(:)
                            end if
                        end do
                    end do
                    helper = helper(:)/(2.0*real(local_nodes(node_index)%num_neighbours))
                    
                    do k = 1,dimensiona
                        index = index+1
                        if (index.gt.node_snd_count(cpu) * dimensiona) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index) = helper(k)
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
                    count = node_snd_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * dimensiona
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

            local_nodes(node_index)%relaxation_velocity(:) = 0.0
            counter = 0

            centre(:) = zero

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1
                cell_index = local_nodes(node_index)%local_neighbours(i)

                helper(:) = zero
                do face_index = 1, ielem(n, cell_index)%ifca
                    do face_node_index = 1, 2
                        if (ielem(n, cell_index)%nodes_faces(face_index, face_node_index).eq.node_index) then
                            other_face_node_index = 3 - face_node_index ! swaps 1 with 2 and vice versa
                            other_node_index = ielem(n, cell_index)%nodes_faces(face_index, other_face_node_index)
                            new_position(1:dimensiona) = local_nodes(other_node_index)%positions(position_index,:) &
                                                        + (d_t*local_nodes(other_node_index)%velocity(1:dimensiona))
                            helper(:) = helper(:) + new_position(:)
                        end if
                    end do
                end do
                helper = helper(:)/(2.0*real(local_nodes(node_index)%num_neighbours))

                centre(:) = centre(:) + helper(:)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*dimensiona
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * dimensiona) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        helper(k) = node_rcv_buffer(cpu_index)%data(index)
                    end do
                    centre(:) = centre(:) + helper(:)
                end do
            end do

            local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (centre(1:dimensiona) - local_nodes(node_index)%positions(position_index,1:dimensiona)) / d_t
        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_moved_node_centre_relaxation_velocity





SUBROUTINE find_moved_node_relaxation_velocity(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, counter, node_index, cell_node_index, cell_index, cpu_index, cpu, index, rho_index
    integer:: M
    real::rho, v
    real,dimension(1:nof_variables)::copy
    real::dummuy_MP_PINFl, dummy_gammal
    real,dimension(1:dimensiona)::helper_centre_position
    real,dimension(1:max_num_node_neighbours,1:dimensiona)::centre_positions
    real,dimension(1:dimensiona,1:4)::cell_nodes
    ! real,dimension(1:max_num_node_neighbours)::weights


    integer,dimension(2*isize)::requests
    integer::num_requests, count

    M = omp_get_thread_num()
    
    rho_index = 1

    if (num_values_to_send_per_node.lt.dimensiona) then
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
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * dimensiona
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    do k = 1, ielem(N, cell_index)%nonodes
                        cell_node_index = ielem(N, cell_index)%nodes_counterclockwise(k)
                        cell_nodes(:,k) = local_nodes(cell_node_index)%positions(position_index,:)
                        cell_nodes(:,k) = cell_nodes(:,k) + (local_nodes(cell_node_index)%velocity(1:dimensiona)*d_t)
                    end do

                    call cell_centre(ielem(N, cell_index)%nonodes, cell_nodes(:,1:ielem(N, cell_index)%nonodes), helper_centre_position(:))

                    do k = 1,dimensiona
                        if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
                            print *, "copying too much data to send buffer from", N, "to", cpu 
                        end if
                        node_snd_buffer(cpu)%data(index+k) = helper_centre_position(k)
                    end do
                    index = index + dimensiona
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
                    count = node_snd_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * dimensiona
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

            local_nodes(node_index)%relaxation_velocity(:) = 0.0
            counter = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1

                cell_index = local_nodes(node_index)%local_neighbours(i)
                do k = 1, ielem(N, cell_index)%nonodes
                    cell_node_index = ielem(N, cell_index)%nodes_counterclockwise(k)
                    cell_nodes(:,k) = local_nodes(cell_node_index)%positions(position_index,:)
                    cell_nodes(:,k) = cell_nodes(:,k) + (local_nodes(cell_node_index)%lagrangian_velocity(1:dimensiona)*d_t)
                end do

                call cell_centre(ielem(N, cell_index)%nonodes, cell_nodes(:,1:ielem(N, cell_index)%nonodes), helper_centre_position(:))

                centre_positions(counter, :) = helper_centre_position(:)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*dimensiona
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    counter = counter+1
                    do k = 1, dimensiona
                        if ((index+k).gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index+k)
                    end do
                    index = index + dimensiona
                end do
            end do

            if (counter.ne.local_nodes(node_index)%num_neighbours) then
                print *, "something went wrong counter =/= local_nodes(node_index)%num_neighbours"
            end if
            if (relaxation_centre_type.eq.1) then
                call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
            else if (relaxation_centre_type.eq.2) then
                call pseudoVoronoi_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
                if ((local_nodes(node_index)%velocity(1).le.xmin(n)).or.&
                        (local_nodes(node_index)%velocity(1).ge.xmax(n)).or.&
                        (local_nodes(node_index)%velocity(2).le.ymin(n)).or.&
                        (local_nodes(node_index)%velocity(2).ge.ymax(n))) then
                    call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
                end if
            else
                print*,"invalid centre algorithm in node relaxation velocity"
            end if

            local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (local_nodes(node_index)%relaxation_velocity(1:dimensiona) - (local_nodes(node_index)%positions(position_index,1:dimensiona)+(local_nodes(node_index)%lagrangian_velocity(1:dimensiona)*d_t))) / d_t
        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_moved_node_relaxation_velocity






! SUBROUTINE find_node_relaxation_velocity_and_normalized_density_gradient(stage, position_index, d_t, N)
!     implicit none
!     integer,intent(in)::stage, position_index, N
!     real,intent(in)::d_t
!     integer::i, j, k, iter, counter, node_index, cell_index, cpu_index, cpu, index, rho_index
!     integer::num_neighbours
!     integer:: M
!     real::rho, v
!     real,dimension(1:nof_variables)::copy
!     real::dummy_MP_PINFl, dummy_gammal
!     real,dimension(1:dimensiona)::helper_centre_position, node_center
!     real,dimension(1:max_num_node_neighbours,1:dimensiona)::centre_positions
!     real,dimension(1:max_num_node_neighbours)::rho_vector
!     ! real,dimension(1:(dimensiona+1),1:max_num_node_neighbours)::At
!     real,dimension(1:max_num_node_neighbours,1:(dimensiona+1))::A
!     real,dimension(1:(dimensiona+1),1:(dimensiona+1))::AtA
!     real,dimension(1:(dimensiona+1))::At_rho_vector, x
!     real,dimension(1:dimensiona)::coord_diff
!     real::rho_difference, distance2
!     real::swap_helper, factor
!     real::gradient_magnitude, my_max_gradient_magnitude
!     logical::solvable

!     integer,dimension(2*isize)::requests
!     integer::num_requests, count

!     M = omp_get_thread_num()
    
!     rho_index = 1

!     max_gradient_magnitude = 0.0

!     if (num_values_to_send_per_node.lt.(dimensiona+1)) then
!         print *,"something went wrong sorry :("
!         call abort
!     end if

!     ! do cpu_index = 0, isize-1
!     !     index = 0
!     !     do i = 1, node_snd_count(cpu_index)
!     !         do j = 1, dimensiona
!     !             index = index +1
!     !             node_snd_buffer(cpu_index)%data(index) = 1000.0
!     !         end do
!     !     end do
!     !     if (index.ne.dimensiona*node_snd_count(cpu_index)) then
!     !         print *, "something went wrong with rezeroing the send buffer"
!     !     else
!     !         print *, "snd_buffer on CPU", N, "thread", M, "to", cpu_index, "reset with", node_snd_count(cpu_index)*dimensiona, "values"
!     !     end if
!     ! end do
!     ! !$omp barrier 

!     !$omp do
!         do iter = 1,my_num_interface_nodes 
!             node_index = local_interface_nodes(iter)
!             ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
!             do cpu_index = 1,local_nodes(node_index)%num_cpus
!                 cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
!                 index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * num_values_to_send_per_node
!                 do j = 1, local_nodes(node_index)%num_local_neighbours
!                     cell_index = local_nodes(node_index)%local_neighbours(j)
!                     copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
!                     call cons2prim(N, copy, dummy_MP_PINFl, dummy_gammal)
!                     rho = copy(rho_index)

!                     helper_centre_position(1) = ielem(N, cell_index)%xxc
!                     helper_centre_position(2) = ielem(N, cell_index)%yyc
!                     if (dimensiona.eq.3) then
!                         helper_centre_position(3) = ielem(N, cell_index)%zzc
!                     end if
!                     do k = 1, dimensiona
!                         if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
!                             print *, "copying too much data to send buffer from", N, "to", cpu 
!                         end if
!                         node_snd_buffer(cpu)%data(index+k) = helper_centre_position(k)
!                     end do
!                     if ((index+dimensiona+1).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
!                         print *, "copying too much data to send buffer from", N, "to", cpu 
!                     end if
!                     node_snd_buffer(cpu)%data(index+dimensiona+1) = rho
!                     index = index + num_values_to_send_per_node
!                 end do
!             end do
!         end do
!     !$omp end do

!     ! do cpu_index = 0, isize-1
!     !     index = 0
!     !     do i = 1, node_rcv_count(cpu_index)
!     !         do j = 1, dimensiona
!     !             index = index +1
!     !             node_rcv_buffer(cpu_index)%data(index) = 1000000.0
!     !         end do
!     !     end do
!     !     if (index.ne.dimensiona*node_rcv_count(cpu_index)) then
!     !         print *, "something went wrong with rezeroing the receive buffer"
!     !     else
!     !         print *, "rcv_buffer on CPU", N, "thread", M, "from", cpu_index, "reset with", node_rcv_count(cpu_index)*dimensiona, "values"
!     !     end if
!     ! end do

!     !$omp barrier

!     num_requests = 0
!     !$omp master
!         ! print *, "inside send_rcv part on CPU", N
!         do cpu_index = 0, isize-1
!             if (cpu_index.ne.N) then
!                 if (node_snd_count(cpu_index).gt.0) then
!                     count = node_snd_count(cpu_index) * num_values_to_send_per_node
!                     num_requests = num_requests + 1
!                     CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
!                     ! print *, "sending from", N, "to", cpu_index, count, "values"
!                 end if
!                 if (node_rcv_count(cpu_index).gt.0) then
!                     count = node_rcv_count(cpu_index) * num_values_to_send_per_node
!                     num_requests = num_requests + 1
!                     CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
!                     ! print *, N, "waiting to receive", count, "values from", cpu_index
!                 end if
!             else
!                 if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
!                     print *,"send receive count missmatch on CPU", n
!                     call abort
!                 end if
!                 count = node_snd_count(cpu_index) * num_values_to_send_per_node
!                 do i = 1, count
!                     node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
!                 end do
!             end if
!         end do

!         ! print*,"CPU", n, "witing on", num_requests, "requests"
!         CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

!     !$omp end master

!     !$omp barrier

!     ! !$omp parallel do reduction(max: max_gradient_magnitude)
!     ! !$omp do reduction(max: max_gradient_magnitude)
!     !$omp do
!         do node_index = 1, kmaxn 

!             local_nodes(node_index)%relaxation_velocity(:) = 0.0
!             counter = 0

!             do i = 1, local_nodes(node_index)%num_local_neighbours
!                 counter = counter + 1
!                 cell_index = local_nodes(node_index)%local_neighbours(i)
!                 copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
!                 call cons2prim(N, copy, dummy_MP_PINFl, dummy_gammal)

!                 helper_centre_position(1) = ielem(N, cell_index)%xxc
!                 helper_centre_position(2) = ielem(N, cell_index)%yyc
!                 if (dimensiona.eq.3) then
!                     helper_centre_position(3) = ielem(N, cell_index)%zzc
!                 end if
!                 centre_positions(counter, :) = helper_centre_position(:)

!                 rho_vector(counter) = copy(rho_index)
!             end do

!             do i = 1, local_nodes(node_index)%num_cpus
!                 cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
!                 index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
!                 do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
!                     counter = counter+1
!                     do k = 1, dimensiona
!                         if ((index+k).gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
!                             print *, "copying too much data from receive buffer from", N, "to", cpu 
!                         end if
!                         centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index+k)
!                     end do
!                     rho_vector(counter) = node_rcv_buffer(cpu_index)%data(index+dimensiona+1)
!                     index = index + num_values_to_send_per_node
!                 end do
!             end do

!             if (counter.ne.local_nodes(node_index)%num_neighbours) then
!                 print *, "something went wrong counter =/= local_nodes(node_index)%num_neighbours"
!             end if
!             if (relaxation_centre_type.eq.1) then
!                 call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
!             else if (relaxation_centre_type.eq.2) then
!                 call pseudoVoronoi_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
!                 if ((local_nodes(node_index)%velocity(1).le.xmin(n)).or.&
!                         (local_nodes(node_index)%velocity(1).ge.xmax(n)).or.&
!                         (local_nodes(node_index)%velocity(2).le.ymin(n)).or.&
!                         (local_nodes(node_index)%velocity(2).ge.ymax(n))) then
!                     call polygon_centre(local_nodes(node_index)%num_neighbours, centre_positions(1:counter,1:dimensiona), local_nodes(node_index)%relaxation_velocity(1:dimensiona))
!                 end if
!             else
!                 print*,"invalid centre algorithm in node relaxation velocity"
!             end if
!             local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (local_nodes(node_index)%relaxation_velocity(1:dimensiona) - local_nodes(node_index)%positions(position_index,1:dimensiona)) / d_t

!             local_nodes(node_index)%density_gradient(:) = zero 
!             local_nodes(node_index)%normalized_density_gradient_magnitude = zero

!             num_neighbours = local_nodes(node_index)%num_neighbours
!             if (num_neighbours.gt.max_num_node_neighbours) then
!                 print*,"overfilled matrix"
!             end if
!             solvable = .false.
!             if (num_neighbours.ge.3) then
!                 do i = 1, num_neighbours 
!                     A(i, 1) = 1.0
!                     do j = 1, dimensiona
!                         A(i, j+1) = centre_positions(i, j) - local_nodes(node_index)%positions(position_index,j)
!                     end do
!                 end do

!                 At_rho_vector(:) = zero
!                 do i = 1, (dimensiona+1)
!                     do j = 1, num_neighbours
!                         ! At(i,j) = A(j,i)
!                         ! At_rho_vector(i) = At_rho_vector(i) + At(i,j)*rho_vector(j)
!                         At_rho_vector(i) = At_rho_vector(i) + A(j, i)*rho_vector(j)
!                     end do
!                 end do

!                 AtA(:,:) = zero
!                 do i = 1, (dimensiona+1)
!                     do j = 1, (dimensiona+1)
!                         do  k = 1, num_neighbours
!                             ! AtA(i,j) = At(i,j) + (At(i,k)*A(k,j))
!                             AtA(i,j) = AtA(i,j) + (A(k, i)*A(k, j))
!                         end do
!                     end do
!                 end do

!                 if (AtA(1,1).ne.real(num_neighbours)) then
!                     print *, "linear solver failuer\nAta(1,1)=/=num_neighbours"
!                     call abort()
!                 end if

!                 solvable = linear_solve(1+dimensiona, AtA, At_rho_vector, x)

!                 if (solvable) then
!                     gradient_magnitude = zero
!                     do i = 1, dimensiona
!                         local_nodes(node_index)%density_gradient(i) = x(i+1)
!                         gradient_magnitude = gradient_magnitude + (x(i+1)*x(i+1))
!                     end do
!                     gradient_magnitude = sqrt(gradient_magnitude)
!                     local_nodes(node_index)%normalized_density_gradient_magnitude = gradient_magnitude

!                     ! if (gradient_magnitude.gt.max_gradient_magnitude) then
!                     !     max_gradient_magnitude = gradient_magnitude
!                     ! end if
!                 end if
!             end if

!             if ((.not.solvable).and.(num_neighbours.ge.2)) then

!                 rho_difference = rho_vector(2) - rho_vector(1)
!                 distance2 = zero
!                 coord_diff(:) = centre_positions(2,:) - centre_positions(1,:)
!                 do i = 1, dimensiona
!                     distance2 = distance2 + (coord_diff(i)*coord_diff(i))
!                 end do
!                 gradient_magnitude = abs(rho_difference) / sqrt(distance2)

!                 local_nodes(node_index)%density_gradient(1:dimensiona) = (coord_diff(1:dimensiona) / sqrt(distance2)) * gradient_magnitude
!                 local_nodes(node_index)%normalized_density_gradient_magnitude = gradient_magnitude

!                 ! if (gradient_magnitude.gt.my_max_gradient_magnitude) then
!                 !     my_max_gradient_magnitude = gradient_magnitude
!                 ! end if

!             else ! num_neighbours.le.1
!                 ! print*,"cannot approximate density gradient" ! it should be a corner cell which normally does not move so it shouldn't be an issue
!             end if
                
!         end do
!     !$omp end do
!     ! !$omp end parallel do

!     max_gradient_magnitude = 0.0
!     !$omp barrier
    
!     !$omp do reduction(max: max_gradient_magnitude)
!         do node_index = 1, kmaxn 
!             if (local_nodes(node_index)%normalized_density_gradient_magnitude.gt.max_gradient_magnitude) then
!                 max_gradient_magnitude = local_nodes(node_index)%normalized_density_gradient_magnitude
!             end if
!         end do
!     !$omp end do
!     ! print*, "before", N, M, max_gradient_magnitude

!     !$omp barrier

!     !$omp master
!         my_max_gradient_magnitude = max_gradient_magnitude
!         CALL MPI_ALLREDUCE(my_max_gradient_magnitude, max_gradient_magnitude, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
!     !$omp end master

!     !$omp barrier

!     my_max_gradient_magnitude = max_gradient_magnitude
!     ! print*, "after", N, M, my_max_gradient_magnitude

!     !$omp do
!     do node_index = 1, kmaxn 
!         local_nodes(node_index)%density_gradient(:) = local_nodes(node_index)%density_gradient(:) / my_max_gradient_magnitude
!         local_nodes(node_index)%normalized_density_gradient_magnitude = local_nodes(node_index)%normalized_density_gradient_magnitude / my_max_gradient_magnitude
!     end do
!     !$omp end do
    
!     !$omp barrier

!     !$omp master
!         call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!     !$omp end master

! END SUBROUTINE find_node_relaxation_velocity_and_normalized_density_gradient





SUBROUTINE find_node_normalized_density_gradient(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, counter, node_index, cell_index, cpu_index, cpu, index, rho_index
    integer::num_neighbours
    integer:: M
    real::rho, v
    real,dimension(1:nof_variables)::copy
    real::dummy_MP_PINFl, dummy_gammal
    real,dimension(1:dimensiona)::helper_centre_position, node_center
    real,dimension(1:max_num_node_neighbours,1:dimensiona)::centre_positions
    real,dimension(1:max_num_node_neighbours)::rho_vector
    ! real,dimension(1:(dimensiona+1),1:max_num_node_neighbours)::At
    real,dimension(1:max_num_node_neighbours,1:(dimensiona+1))::A
    real,dimension(1:(dimensiona+1),1:(dimensiona+1))::AtA
    real,dimension(1:(dimensiona+1))::At_rho_vector, x
    real,dimension(1:dimensiona)::coord_diff
    real::rho_difference, distance2
    real::swap_helper, factor
    real::gradient_magnitude, my_max_gradient_magnitude
    logical::solvable

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()
    rho_index = 1

    max_gradient_magnitude = 0.0

    if (num_values_to_send_per_node.lt.(dimensiona+1)) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (dimensiona+1)
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    call cons2prim(N, copy, dummy_MP_PINFl, dummy_gammal)
                    rho = copy(rho_index)

                    helper_centre_position(1) = ielem(N, cell_index)%xxc
                    helper_centre_position(2) = ielem(N, cell_index)%yyc
                    if (dimensiona.eq.3) then
                        helper_centre_position(3) = ielem(N, cell_index)%zzc
                    end if
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_snd_count(cpu) * (dimensiona+1)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_normalized_density_gradient"
                        end if
                        node_snd_buffer(cpu)%data(index) = helper_centre_position(k)
                    end do
                    index = index+1
                    if ((index).gt.node_snd_count(cpu) * (dimensiona+1)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_normalized_density_gradient"
                    end if
                    node_snd_buffer(cpu)%data(index) = rho
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
                    count = node_snd_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * (dimensiona+1)
                do i = 1, count
                    node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
                end do
            end if
        end do

        ! print*,"CPU", n, "witing on", num_requests, "requests"
        CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

    !$omp end master

    !$omp barrier

    ! !$omp parallel do reduction(max: max_gradient_magnitude)
    ! !$omp do reduction(max: max_gradient_magnitude)
    !$omp do
        do node_index = 1, kmaxn 

            counter = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1
                cell_index = local_nodes(node_index)%local_neighbours(i)
                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummy_MP_PINFl, dummy_gammal)

                helper_centre_position(1) = ielem(N, cell_index)%xxc
                helper_centre_position(2) = ielem(N, cell_index)%yyc
                if (dimensiona.eq.3) then
                    helper_centre_position(3) = ielem(N, cell_index)%zzc
                end if
                centre_positions(counter, :) = helper_centre_position(:)

                rho_vector(counter) = copy(rho_index)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*(dimensiona+1)
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    counter = counter+1
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index)
                    end do
                    index = index+1
                    if (index.gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                        print *, "copying too much data from receive buffer from", N, "to", cpu 
                    end if
                    rho_vector(counter) = node_rcv_buffer(cpu_index)%data(index)
                end do
            end do

            if (counter.ne.local_nodes(node_index)%num_neighbours) then
                print *, "something went wrong counter =/= local_nodes(node_index)%num_neighbours"
            end if

            ! local_nodes(node_index)%density_gradient(:) = zero 
            local_nodes(node_index)%normalized_density_gradient_magnitude = zero

            num_neighbours = local_nodes(node_index)%num_neighbours
            if (num_neighbours.gt.max_num_node_neighbours) then
                print*,"overfilled matrix"
            end if
            solvable = .false.
            if (num_neighbours.ge.3) then
                do i = 1, num_neighbours 
                    A(i, 1) = 1.0
                    do j = 1, dimensiona
                        A(i, j+1) = centre_positions(i, j) - local_nodes(node_index)%positions(position_index,j)
                    end do
                end do

                At_rho_vector(:) = zero
                do i = 1, (dimensiona+1)
                    do j = 1, num_neighbours
                        ! At(i,j) = A(j,i)
                        ! At_rho_vector(i) = At_rho_vector(i) + At(i,j)*rho_vector(j)
                        At_rho_vector(i) = At_rho_vector(i) + A(j, i)*rho_vector(j)
                    end do
                end do

                AtA(:,:) = zero
                do i = 1, (dimensiona+1)
                    do j = 1, (dimensiona+1)
                        do  k = 1, num_neighbours
                            ! AtA(i,j) = At(i,j) + (At(i,k)*A(k,j))
                            AtA(i,j) = AtA(i,j) + (A(k, i)*A(k, j))
                        end do
                    end do
                end do

                if (AtA(1,1).ne.real(num_neighbours)) then
                    print *, "linear solver failuer\nAta(1,1)=/=num_neighbours"
                    call abort()
                end if

                solvable = linear_solve(1+dimensiona, AtA, At_rho_vector, x)

                if (solvable) then
                    gradient_magnitude = zero
                    do i = 1, dimensiona
                        ! local_nodes(node_index)%density_gradient(i) = x(i+1)
                        gradient_magnitude = gradient_magnitude + (x(i+1)*x(i+1))
                    end do
                    gradient_magnitude = sqrt(gradient_magnitude)
                    local_nodes(node_index)%normalized_density_gradient_magnitude = gradient_magnitude

                    ! if (gradient_magnitude.gt.max_gradient_magnitude) then
                    !     max_gradient_magnitude = gradient_magnitude
                    ! end if
                else
                    print*, "linear solver failure at", local_nodes(node_index)%positions(position_index,:)
                end if
            end if

            if ((.not.solvable).and.(num_neighbours.ge.2)) then
                rho_difference = rho_vector(2) - rho_vector(1)
                distance2 = zero
                coord_diff(:) = centre_positions(2,:) - centre_positions(1,:)
                do i = 1, dimensiona
                    distance2 = distance2 + (coord_diff(i)*coord_diff(i))
                end do
                gradient_magnitude = abs(rho_difference) / sqrt(distance2)

                ! local_nodes(node_index)%density_gradient(1:dimensiona) = (coord_diff(1:dimensiona) / sqrt(distance2)) * gradient_magnitude
                local_nodes(node_index)%normalized_density_gradient_magnitude = gradient_magnitude

                ! if (gradient_magnitude.gt.my_max_gradient_magnitude) then
                !     my_max_gradient_magnitude = gradient_magnitude
                ! end if

            else ! num_neighbours.le.1
                ! print*,"cannot approximate density gradient" ! it should be a corner cell which normally does not move so it shouldn't be an issue
            end if
                
        end do
    !$omp end do
    ! !$omp end parallel do

    max_gradient_magnitude = 0.0
    !$omp barrier
    
    !$omp do reduction(max: max_gradient_magnitude)
        do node_index = 1, kmaxn
            if (initcond.eq.102) then
                if (.not.((local_nodes(node_index)%positions(position_index,1).lt.(9.5*t+0.1)).and.(local_nodes(node_index)%positions(position_index,2).lt.0.1))) then
                    if (local_nodes(node_index)%normalized_density_gradient_magnitude.gt.max_gradient_magnitude) then
                        max_gradient_magnitude = local_nodes(node_index)%normalized_density_gradient_magnitude
                    end if
                end if
            else
                if (local_nodes(node_index)%normalized_density_gradient_magnitude.gt.max_gradient_magnitude) then
                    max_gradient_magnitude = local_nodes(node_index)%normalized_density_gradient_magnitude
                end if
            end if
        end do
    !$omp end do
    ! print*, "before", N, M, max_gradient_magnitude

    !$omp barrier

    !$omp master
        my_max_gradient_magnitude = max_gradient_magnitude
        CALL MPI_ALLREDUCE(my_max_gradient_magnitude, max_gradient_magnitude, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

    my_max_gradient_magnitude = max_gradient_magnitude
    ! print*, "after", N, M, my_max_gradient_magnitude

    !$omp do
    do node_index = 1, kmaxn 
        ! local_nodes(node_index)%density_gradient(:) = local_nodes(node_index)%density_gradient(:) / my_max_gradient_magnitude
        local_nodes(node_index)%normalized_density_gradient_magnitude = local_nodes(node_index)%normalized_density_gradient_magnitude / my_max_gradient_magnitude
    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_node_normalized_density_gradient





SUBROUTINE find_node_normalized_vf_gradient(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, counter, node_index, cell_index, cpu_index, cpu, index, vf_index
    integer::num_neighbours
    ! integer:: M
    real::vf
    real,dimension(1:nof_variables)::copy
    real::dummy_MP_PINFl, dummy_gammal
    real,dimension(1:dimensiona)::helper_centre_position, node_center
    real,dimension(1:max_num_node_neighbours,1:dimensiona)::centre_positions
    real,dimension(1:max_num_node_neighbours)::vf_vector
    ! real,dimension(1:(dimensiona+1),1:max_num_node_neighbours)::At
    real,dimension(1:max_num_node_neighbours,1:(dimensiona+1))::A
    real,dimension(1:(dimensiona+1),1:(dimensiona+1))::AtA
    real,dimension(1:(dimensiona+1))::At_vf_vector, x
    real,dimension(1:dimensiona)::coord_diff
    real::vf_difference, distance2
    real::swap_helper, factor
    real::gradient_magnitude, my_max_gradient_magnitude
    logical::solvable

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    max_gradient_magnitude = 0.0

    if (num_values_to_send_per_node.lt.(dimensiona+1)) then
        print *,"something went wrong sorry :("
        call abort
    end if
    if (governingequations.ne.-1) then
        print*,"volume fraction is only defined for diffuse-interface multi-phase flows"
        call abort
    end if
    vf_index = 7

    !$omp do
        do iter = 1,my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
            do cpu_index = 1,local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * (dimensiona+1)
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                    call cons2prim(N, copy, dummy_MP_PINFl, dummy_gammal)
                    vf = copy(vf_index)

                    helper_centre_position(1) = ielem(N, cell_index)%xxc
                    helper_centre_position(2) = ielem(N, cell_index)%yyc
                    if (dimensiona.eq.3) then
                        helper_centre_position(3) = ielem(N, cell_index)%zzc
                    end if
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_snd_count(cpu) * (dimensiona+1)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_normalized_density_gradient"
                        end if
                        node_snd_buffer(cpu)%data(index) = helper_centre_position(k)
                    end do
                    index = index+1
                    if ((index).gt.node_snd_count(cpu) * (dimensiona+1)) then
                        print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_normalized_density_gradient"
                    end if
                    node_snd_buffer(cpu)%data(index) = vf
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
                    count = node_snd_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (dimensiona+1)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * (dimensiona+1)
                do i = 1, count
                    node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
                end do
            end if
        end do

        ! print*,"CPU", n, "witing on", num_requests, "requests"
        CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

    !$omp end master

    !$omp barrier

    ! !$omp parallel do reduction(max: max_gradient_magnitude)
    ! !$omp do reduction(max: max_gradient_magnitude)
    !$omp do
        do node_index = 1, kmaxn 

            counter = 0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                counter = counter + 1
                cell_index = local_nodes(node_index)%local_neighbours(i)
                copy(:) = u_c(cell_index)%val(stage,1:nof_variables)
                call cons2prim(N, copy, dummy_MP_PINFl, dummy_gammal)

                helper_centre_position(1) = ielem(N, cell_index)%xxc
                helper_centre_position(2) = ielem(N, cell_index)%yyc
                if (dimensiona.eq.3) then
                    helper_centre_position(3) = ielem(N, cell_index)%zzc
                end if
                centre_positions(counter, :) = helper_centre_position(:)

                vf_vector(counter) = copy(vf_index)
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*(dimensiona+1)
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    counter = counter+1
                    do k = 1, dimensiona
                        index = index+1
                        if (index.gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                            print *, "copying too much data from receive buffer from", N, "to", cpu 
                        end if
                        centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index)
                    end do
                    index = index+1
                    if (index.gt.node_rcv_count(cpu_index) * (dimensiona+1)) then
                        print *, "copying too much data from receive buffer from", N, "to", cpu 
                    end if
                    vf_vector(counter) = node_rcv_buffer(cpu_index)%data(index)
                end do
            end do

            if (counter.ne.local_nodes(node_index)%num_neighbours) then
                print *, "something went wrong counter =/= local_nodes(node_index)%num_neighbours"
            end if

            ! local_nodes(node_index)%density_gradient(:) = zero 
            local_nodes(node_index)%normalized_vf_gradient_magnitude = zero

            num_neighbours = local_nodes(node_index)%num_neighbours
            if (num_neighbours.gt.max_num_node_neighbours) then
                print*,"overfilled matrix"
            end if
            solvable = .false.
            if (num_neighbours.ge.3) then
                do i = 1, num_neighbours 
                    A(i, 1) = 1.0
                    do j = 1, dimensiona
                        A(i, j+1) = centre_positions(i, j) - local_nodes(node_index)%positions(position_index,j)
                    end do
                end do

                At_vf_vector(:) = zero
                do i = 1, (dimensiona+1)
                    do j = 1, num_neighbours
                        At_vf_vector(i) = At_vf_vector(i) + A(j, i)*vf_vector(j)
                    end do
                end do

                AtA(:,:) = zero
                do i = 1, (dimensiona+1)
                    do j = 1, (dimensiona+1)
                        do  k = 1, num_neighbours
                            AtA(i,j) = AtA(i,j) + (A(k, i)*A(k, j))
                        end do
                    end do
                end do

                if (AtA(1,1).ne.real(num_neighbours)) then
                    print *, "linear solver failuer\nAta(1,1)=/=num_neighbours"
                    call abort()
                end if

                solvable = linear_solve(1+dimensiona, AtA, At_vf_vector, x)

                if (solvable) then
                    gradient_magnitude = zero
                    do i = 1, dimensiona
                        ! local_nodes(node_index)%density_gradient(i) = x(i+1)
                        gradient_magnitude = gradient_magnitude + (x(i+1)*x(i+1))
                    end do
                    gradient_magnitude = sqrt(gradient_magnitude)
                    local_nodes(node_index)%normalized_vf_gradient_magnitude = gradient_magnitude

                    ! if (gradient_magnitude.gt.max_gradient_magnitude) then
                    !     max_gradient_magnitude = gradient_magnitude
                    ! end if
                else
                    print*, "linear solver failure at", local_nodes(node_index)%positions(position_index,:)
                end if
            end if

            if ((.not.solvable).and.(num_neighbours.ge.2)) then
                vf_difference = vf_vector(2) - vf_vector(1)
                distance2 = zero
                coord_diff(:) = centre_positions(2,:) - centre_positions(1,:)
                do i = 1, dimensiona
                    distance2 = distance2 + (coord_diff(i)*coord_diff(i))
                end do
                gradient_magnitude = abs(vf_difference) / sqrt(distance2)

                ! local_nodes(node_index)%density_gradient(1:dimensiona) = (coord_diff(1:dimensiona) / sqrt(distance2)) * gradient_magnitude
                local_nodes(node_index)%normalized_vf_gradient_magnitude = gradient_magnitude

                ! if (gradient_magnitude.gt.my_max_gradient_magnitude) then
                !     my_max_gradient_magnitude = gradient_magnitude
                ! end if

            else ! num_neighbours.le.1
                ! print*,"cannot approximate density gradient" ! it should be a corner cell which normally does not move so it shouldn't be an issue
            end if
                
        end do
    !$omp end do
    ! !$omp end parallel do

    max_gradient_magnitude = 0.0
    !$omp barrier
    
    !$omp do reduction(max: max_gradient_magnitude)
    do node_index = 1, kmaxn
        if (local_nodes(node_index)%normalized_vf_gradient_magnitude.gt.max_gradient_magnitude) then
            max_gradient_magnitude = local_nodes(node_index)%normalized_vf_gradient_magnitude
        end if
    end do
    !$omp end do
    ! print*, "before", N, M, max_gradient_magnitude

    !$omp barrier

    !$omp master
        my_max_gradient_magnitude = max_gradient_magnitude
        CALL MPI_ALLREDUCE(my_max_gradient_magnitude, max_gradient_magnitude, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

    my_max_gradient_magnitude = max_gradient_magnitude
    ! print*, "after", N, M, my_max_gradient_magnitude

    !$omp do
    do node_index = 1, kmaxn 
        ! local_nodes(node_index)%density_gradient(:) = local_nodes(node_index)%density_gradient(:) / my_max_gradient_magnitude
        local_nodes(node_index)%normalized_vf_gradient_magnitude = local_nodes(node_index)%normalized_vf_gradient_magnitude / my_max_gradient_magnitude
    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_node_normalized_vf_gradient





SUBROUTINE find_node_volume_ratio(stage, node_position_index, N)
    implicit none
    integer,intent(in)::stage, node_position_index, N
    integer::i, j, k, iter, node_index, cell_index, cpu_index, cpu, index
    ! integer:: M
    real::volume, min_volume, max_volume

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.1) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1, my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            do cpu_index = 1, local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * 1
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    volume = ielem(n, cell_index)%moving_volume(node_position_index)
                    index = index + 1
                    node_snd_buffer(cpu)%data(index) = volume
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
                    count = node_snd_count(cpu_index) * 1
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * 1
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * 1
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

            local_nodes(node_index)%volume_ratio = 1.0
            min_volume = 10000000000.0
            max_volume = 0.0

            do i = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(i)

                volume = ielem(n, cell_index)%moving_volume(node_position_index)
                if (volume.lt.min_volume) then
                    min_volume = volume
                end if
                if (volume.gt.max_volume) then
                    max_volume = volume
                end if
            end do

            do i = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
                index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1) * 1
                do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
                    index = index+1
                    volume = node_rcv_buffer(cpu_index)%data(index)
                    if (volume.lt.min_volume) then
                        min_volume = volume
                    end if
                    if (volume.gt.max_volume) then
                        max_volume = volume
                    end if
                end do
            end do

            if (max_volume.lt.min_volume) then
                print*,"max_volume < min_volume", max_volume, min_volume
            end if
            if (min_volume.eq.zero) then
                print*,"min_volume = 0"
            end if

            local_nodes(node_index)%volume_ratio = max_volume / min_volume

        end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

END SUBROUTINE find_node_volume_ratio





subroutine CombineNodeVelocities(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::stage, position_index, N
    real::d_t
    integer::node_index
    integer::i
    real,dimension(1:dimensiona)::initial_velocity1, initial_velocity2
    real,dimension(1:dimensiona)::option1, option2
    real::gradient_copy, fraction, free_stream, val1, val2, gradient_term
    real::local_lagrangian_velocity_multiple, local_relaxation_velocity_multiple
    real::local_lagrangian_velocity_magnitude, local_relaxation_velocity_magnitude


    initial_velocity1 = zero
    initial_velocity2 = zero
    
    select case(moving_mesh_mode)
      case(1)
        !$omp do
        do node_index = 1, kmaxn
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) * lagrangian_mesh_velocity_multiple
        end do
        !$omp end do

      case(2)
        !$omp do
        do node_index = 1, kmaxn
            option1(:) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            option2(:) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            option2(1) = option2(1) - uvel
            option2(2) = option2(2) - vvel
            if (dimensiona.eq.3) then
                option2(3) = option2(3) - wvel
            end if
            local_nodes(node_index)%velocity(1:dimensiona) = min_abs(option1, option2)
        end do
        !$omp end do

      case(3)
        !$omp do
        do node_index = 1, kmaxn
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
        end do
        !$omp end do

      case(4)
        !$omp do
        do node_index = 1, kmaxn
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%relaxation_velocity(1:dimensiona)*mesh_velocity_multiple
        end do
        !$omp end do

      case(5)
        !$omp do
        do node_index = 1, kmaxn
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) * lagrangian_mesh_velocity_multiple &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * relaxation_mesh_velocity_multiple
        end do
        !$omp end do

      case(6)
        !$omp do
        do node_index = 1, kmaxn
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * relaxation_mesh_velocity_multiple
        end do
        !$omp end do

      case(7)
        !$omp do
        do node_index = 1, kmaxn
            option1(:) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            option2(:) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            option2(1) = option2(1) - uvel
            option2(2) = option2(2) - vvel
            if (dimensiona.eq.3) then
                option2(3) = option2(3) - wvel
            end if
            local_nodes(node_index)%velocity(1:dimensiona) = min_abs(option1, option2) &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * relaxation_mesh_velocity_multiple
        end do
        !$omp end do

      case(8)
        !$omp do
        do node_index = 1, kmaxn
            ! if (local_nodes(node_index)%num_neighbours.lt.3) then
            if (local_nodes(node_index)%num_neighbours.lt.2) then
                local_lagrangian_velocity_multiple = 0.0
            else
                if (gradient_copy.lt.zero) then
                    print*, "negative desnity gradient magnitude"
                else if (gradient_copy.lt.gradient_treshold) then
                    local_lagrangian_velocity_multiple = 0.0
                else if (gradient_copy.lt.(2.0*gradient_treshold)) then
                    local_lagrangian_velocity_multiple = (gradient_copy - gradient_treshold) / gradient_treshold
                else if (gradient_copy.le.1.0) then
                    local_lagrangian_velocity_multiple = 1.0
                else
                    print*, "normalized desnity gradient magnitude above 1.0"
                end if
            end if

            if ((local_lagrangian_velocity_multiple.lt.zero).or.(local_lagrangian_velocity_multiple.gt.1.0)) then
                print *, "invalid local lagrangian velocity multiple"
            end if
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) * local_lagrangian_velocity_multiple &
                                                            + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * relaxation_mesh_velocity_multiple
        end do
        !$omp end do

      case(9, 10)
        !$omp do
        do node_index = 1, kmaxn
            
            if ((initcond.eq.102).and.(local_nodes(node_index)%positions(position_index,1).lt.0.5).and.(local_nodes(node_index)%positions(position_index,2).lt.0.25))then
                local_relaxation_velocity_multiple = 1.0
            else
                ! if (local_nodes(node_index)%num_neighbours.lt.3) then
                if (local_nodes(node_index)%num_neighbours.lt.2) then
                    local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
                else
                    gradient_copy = local_nodes(node_index)%normalized_density_gradient_magnitude
                    if (gradient_copy.lt.zero) then
                        print*, "negative desnity gradient magnitude"
                    else if (gradient_copy.lt.lower_gradient_treshold) then
                        local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
                    else if (gradient_copy.le.upper_gradient_treshold) then
                        fraction = (gradient_copy - lower_gradient_treshold) / (upper_gradient_treshold - lower_gradient_treshold)
                        local_relaxation_velocity_multiple =  (fraction * (lower_relaxation_mesh_velocity_multiple - upper_relaxation_mesh_velocity_multiple)) + upper_relaxation_mesh_velocity_multiple
                    else if (gradient_copy.le.1.0) then
                        local_relaxation_velocity_multiple = lower_relaxation_mesh_velocity_multiple
                    else
                        local_relaxation_velocity_multiple = 1.0
                        print*, "normalized desnity gradient magnitude above 1.0"
                    end if
                end if
            end if

            if ((local_lagrangian_velocity_multiple.lt.zero).or.(local_lagrangian_velocity_multiple.gt.1.0)) then
                print *, "invalid local lagrangian velocity multiple"
            end if

            ! local_lagrangian_velocity_magnitude = zero
            ! local_relaxation_velocity_magnitude = zero
            ! do i = 1, dimensiona
            !     local_lagrangian_velocity_magnitude = local_lagrangian_velocity_magnitude + (local_nodes(node_index)%lagrangian_velocity(i)**2)
            !     local_relaxation_velocity_magnitude = local_relaxation_velocity_magnitude + (local_nodes(node_index)%relaxation_velocity(i)**2)
            ! end do
            ! if (dimensiona.eq.2) then
            !     local_lagrangian_velocity_magnitude = max(local_lagrangian_velocity_magnitude, ((uvel*uvel)+(vvel*vvel)))
            ! else
            !     local_lagrangian_velocity_magnitude = max(local_lagrangian_velocity_magnitude, ((uvel*uvel)+(vvel*vvel)+(wvel*wvel)))
            ! end if
            ! if (local_relaxation_velocity_magnitude.gt.(4.0*local_lagrangian_velocity_magnitude)) then
            !     local_lagrangian_velocity_magnitude = sqrt(local_lagrangian_velocity_magnitude) 
            !     local_relaxation_velocity_magnitude = sqrt(local_relaxation_velocity_magnitude)
            !     local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (2.0*local_lagrangian_velocity_magnitude/local_relaxation_velocity_magnitude)*local_nodes(node_index)%relaxation_velocity(1:dimensiona)
            ! end if
            
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * local_relaxation_velocity_multiple
        end do
        !$omp end do

      case(11)
        !$omp do
        do node_index = 1, kmaxn
            option1(:) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            option2(:) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            option2(1) = option2(1) - uvel
            option2(2) = option2(2) - vvel
            if (dimensiona.eq.3) then
                option2(3) = option2(3) - wvel
            end if
            val1 = zero
            val2 = zero
            do i = 1, dimensiona
                val1 = val1 + (option1(i)*option1(i))
                val2 = val2 + (option2(i)*option2(i))
            end do
            val1 = sqrt(val1)
            val2 = sqrt(val2)
            free_stream = sqrt((uvel*uvel)+(vvel*vvel))

            if (val1.le.(0.5*free_stream)) then
                local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona)
            else
                local_nodes(node_index)%velocity(1:dimensiona) = min(val1, val2) * (local_nodes(node_index)%lagrangian_velocity(1:dimensiona) / val1)
            end if
                                                           
            if (local_nodes(node_index)%num_neighbours.lt.2) then
                local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
            else
                gradient_copy = local_nodes(node_index)%normalized_density_gradient_magnitude
                if (gradient_copy.lt.zero) then
                    print*, "negative desnity gradient magnitude"
                else if (gradient_copy.lt.lower_gradient_treshold) then
                    local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
                else if (gradient_copy.le.upper_gradient_treshold) then
                    fraction = (gradient_copy - lower_gradient_treshold) / (upper_gradient_treshold - lower_gradient_treshold)
                    local_relaxation_velocity_multiple =  (fraction * (lower_relaxation_mesh_velocity_multiple - upper_relaxation_mesh_velocity_multiple)) + upper_relaxation_mesh_velocity_multiple
                else if (gradient_copy.le.1.0) then
                    local_relaxation_velocity_multiple = lower_relaxation_mesh_velocity_multiple
                else
                    local_relaxation_velocity_multiple = 1.0
                    print*, "normalized desnity gradient magnitude above 1.0"
                end if
            end if

            if ((local_lagrangian_velocity_multiple.lt.zero).or.(local_lagrangian_velocity_multiple.gt.1.0)) then
                print *, "invalid local lagrangian velocity multiple"
            end if
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * local_relaxation_velocity_multiple
        end do
        !$omp end do

      case(12)
        !$omp do
        do node_index = 1, kmaxn
            ! if (local_nodes(node_index)%num_neighbours.lt.3) then
            if (local_nodes(node_index)%volume_ratio.lt.1.0) then
                print*,"local_nodes(node_index)%volume_ratio < 1.0"
            end if
            if (local_nodes(node_index)%num_neighbours.lt.2) then
                local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
            else
                gradient_copy = local_nodes(node_index)%normalized_density_gradient_magnitude
                local_relaxation_velocity_multiple = (1.0 - gradient_copy)*(1.0 - gradient_copy) * local_nodes(node_index)%volume_ratio * upper_relaxation_mesh_velocity_multiple
            end if
            if (local_relaxation_velocity_multiple.gt.1.0) then
                local_relaxation_velocity_multiple = 1.0
            end if

            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * local_relaxation_velocity_multiple
        end do
        !$omp end do

      case(13)
        !$omp do
        do node_index = 1, kmaxn
            
            if ((initcond.eq.102).and.(local_nodes(node_index)%positions(position_index,1).lt.0.5).and.(local_nodes(node_index)%positions(position_index,2).lt.0.25))then
                local_relaxation_velocity_multiple = 1.0
            else
                ! if (local_nodes(node_index)%num_neighbours.lt.3) then
                if (local_nodes(node_index)%num_neighbours.lt.2) then
                    local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
                else
                    gradient_copy = 0.5*(local_nodes(node_index)%normalized_density_gradient_magnitude + local_nodes(node_index)%normalized_vf_gradient_magnitude)
                    if (gradient_copy.lt.zero) then
                        print*, "negative desnity gradient magnitude"
                    else if (gradient_copy.lt.lower_gradient_treshold) then
                        local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
                    else if (gradient_copy.le.upper_gradient_treshold) then
                        fraction = (gradient_copy - lower_gradient_treshold) / (upper_gradient_treshold - lower_gradient_treshold)
                        local_relaxation_velocity_multiple =  (fraction * (lower_relaxation_mesh_velocity_multiple - upper_relaxation_mesh_velocity_multiple)) + upper_relaxation_mesh_velocity_multiple
                    else if (gradient_copy.le.1.0) then
                        local_relaxation_velocity_multiple = lower_relaxation_mesh_velocity_multiple
                    else
                        local_relaxation_velocity_multiple = 1.0
                        print*, "normalized desnity gradient magnitude above 1.0"
                    end if
                end if
            end if

            if ((local_lagrangian_velocity_multiple.lt.zero).or.(local_lagrangian_velocity_multiple.gt.1.0)) then
                print *, "invalid local lagrangian velocity multiple"
            end if
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * local_relaxation_velocity_multiple
        end do
        !$omp end do

      case(14)
        !$omp do
        do node_index = 1, kmaxn
            if (local_nodes(node_index)%mesh_quality_after.lt.zero) then
                local_relaxation_velocity_multiple = zero
            else
                if (local_nodes(node_index)%mesh_quality_before.lt.zero) then
                    local_relaxation_velocity_multiple = upper_relaxation_mesh_velocity_multiple
                else
                    if (governingequations.eq.-1) then
                        gradient_copy = 0.5*(local_nodes(node_index)%normalized_density_gradient_magnitude + local_nodes(node_index)%normalized_vf_gradient_magnitude)
                    else
                        gradient_copy = local_nodes(node_index)%normalized_density_gradient_magnitude
                    end if

                    if ((gradient_copy.gt.1.0).or.(local_nodes(node_index)%num_neighbours.lt.3)) then
                        gradient_term = zero
                    else if (gradient_copy.ge.upper_gradient_treshold) then
                        gradient_term = 1.0
                    else
                        gradient_term = gradient_copy/upper_gradient_treshold
                    end if

                    local_relaxation_velocity_multiple = (local_nodes(node_index)%mesh_quality_before / real(local_nodes(node_index)%num_neighbours)) - 1.0
                    local_relaxation_velocity_multiple = local_relaxation_velocity_multiple - ((quality_treshold - 1.0)*gradient_term)
                    local_lagrangian_velocity_multiple = local_lagrangian_velocity_multiple * scaling
                    
                end if
            end if
            call clamp(local_relaxation_velocity_multiple, 0.0, upper_relaxation_mesh_velocity_multiple)
            ! if ((local_lagrangian_velocity_multiple.lt.zero).or.(local_lagrangian_velocity_multiple.gt.1.0)) then
            !     print *, "invalid local lagrangian velocity multiple"
            ! end if
            local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%lagrangian_velocity(1:dimensiona) &
                                                           + local_nodes(node_index)%relaxation_velocity(1:dimensiona) * local_relaxation_velocity_multiple
        end do
        !$omp end do

      case default
        print *, "invalid moving mesh mode"
    end select

    !$omp barrier

    ! call fix_concave_cells(stage, position_index, d_t, N)
    if (relaxation_centre_type.lt.5) then
        call fix_moved_concave_cells(stage, position_index, d_t, N)
    end if

    !$omp barrier

end subroutine CombineNodeVelocities





subroutine fix_concave_cells(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::N, position_index, stage
    real,intent(in)::d_t
    integer::node_index, previous_node_index, next_node_index
    integer::cell_index, cell_num_nodes
    integer::cpu, cpu_index
    integer::i, i_minus, i_plus, j, k, index
    real,dimension(1:dimensiona)::v1, v2
    real,dimension(1:dimensiona)::current_node_position, previous_node_position, next_node_position, desired_node_position
    real::len1, len2, v1_cross_v2

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    !$omp do
    do node_index = 1, kmaxn
        local_nodes(node_index)%relaxation_velocity(:) = zero
        current_node_position = local_nodes(node_index)%positions(position_index, 1:dimensiona)

        do k = 1, local_nodes(node_index)%num_local_neighbours
            cell_index = local_nodes(node_index)%local_neighbours(k)
            cell_num_nodes = ielem(N, cell_index)%nonodes
            do i = 1, cell_num_nodes
                if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                    exit
                end if
            end do

            i_minus = i-1
            if (i_minus.eq.0) then
                i_minus = cell_num_nodes
            end if
            i_plus = i+1
            if (i_plus.eq.cell_num_nodes+1) then
                i_plus = 1
            end if
            previous_node_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
            next_node_index = ielem(N, cell_index)%nodes_counterclockwise(i_plus)

            previous_node_position = local_nodes(previous_node_index)%positions(position_index, 1:dimensiona)
            next_node_position = local_nodes(next_node_index)%positions(position_index, 1:dimensiona)

            v1(:) = current_node_position(:) - previous_node_position(:)
		    v2(:) = next_node_position(:) - current_node_position(:)

            len1 = zero
            len2 = zero
            do j = 1, dimensiona
                len1 = len1 + (v1(j)*v1(j))
                len2 = len2 + (v2(j)*v2(j))
            end do
            len1 = sqrt(len1)
            len2 = sqrt(len2)

            v1_cross_v2 = cross_product_2D(v1, v2)

            ! if (v1_cross_v2.lt.zero) then ! concave angle
            if (v1_cross_v2.lt.-0.175*len1*len2) then ! too concave angle

                if ((local_nodes(node_index)%relaxation_velocity(1).ne.zero).or.(local_nodes(node_index)%relaxation_velocity(2).ne.zero)) then
                    print*,"two concave angles before communication"
                end if

                if (cell_num_nodes.eq.3) then
                    print*,"inverted triangle not trying to fix"
                else
                    if (dimensiona.eq.3) then
                        print*,"trying to fix a concave angle in a quadrilateral in cell", cell_index, ielem(n,cell_index)%xxc, ielem(n,cell_index)%yyc, ielem(n,cell_index)%zzc
                    else
                        print*,"trying to fix a concave angle in a quadrilateral in cell", cell_index, ielem(n,cell_index)%xxc, ielem(n,cell_index)%yyc
                    end if
                    desired_node_position(:) = previous_node_position(:) + (len2/(len1+len2))*(next_node_position(:)-previous_node_position(:))

                    local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (desired_node_position(:)-current_node_position(:))/d_t
                end if
            end if
        end do
    end do
    !$omp end do

    !$omp barrier

    !$omp do
    do i = 1, my_num_interface_nodes 
        node_index = local_interface_nodes(i)
        ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
        do cpu_index = 1,local_nodes(node_index)%num_cpus
            cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
            index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * dimensiona
            do j = 1, dimensiona
                node_snd_buffer(cpu)%data(index+j) = local_nodes(node_index)%relaxation_velocity(j)
            end do
            do j = dimensiona+1, local_nodes(node_index)%num_local_neighbours*dimensiona
                node_snd_buffer(cpu)%data(index+j) = zero
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
                    count = node_snd_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * dimensiona
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

        do i = 1, local_nodes(node_index)%num_cpus
            cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
            index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*dimensiona

            if ((local_nodes(node_index)%relaxation_velocity(1).eq.zero).and.(local_nodes(node_index)%relaxation_velocity(2).eq.zero)) then
                do j = 1, dimensiona
                    local_nodes(node_index)%relaxation_velocity(j) = node_rcv_buffer(cpu_index)%data(index+j)
                end do
            else
                if ((node_rcv_buffer(cpu_index)%data(index+1).ne.zero).or.(node_rcv_buffer(cpu_index)%data(index+2).ne.zero)) then
                    print*,"multiple convex angles after communication"
                end if
            end if
        end do

        local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) + local_nodes(node_index)%relaxation_velocity(1:dimensiona)
    end do
    !$omp end do

end subroutine fix_concave_cells





subroutine fix_moved_concave_cells(stage, position_index, d_t, N)
    implicit none
    integer,intent(in)::N, position_index, stage
    real,intent(in)::d_t
    integer::node_index, previous_node_index, next_node_index
    integer::cell_index, cell_num_nodes
    integer::cpu, cpu_index
    integer::i, i_minus, i_plus, j, k, index
    real,dimension(1:dimensiona)::v1, v2
    real,dimension(1:dimensiona)::current_node_position, previous_node_position, next_node_position, desired_node_position
    real::len1, len2, v1_cross_v2

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    !$omp do
    do node_index = 1, kmaxn
        local_nodes(node_index)%relaxation_velocity(:) = zero
        current_node_position = local_nodes(node_index)%positions(position_index, 1:dimensiona) &
                              + (local_nodes(node_index)%velocity(1:dimensiona) * d_t)

        do k = 1, local_nodes(node_index)%num_local_neighbours
            cell_index = local_nodes(node_index)%local_neighbours(k)
            cell_num_nodes = ielem(N, cell_index)%nonodes
            do i = 1, cell_num_nodes
                if (ielem(N, cell_index)%nodes_counterclockwise(i).eq.node_index) then
                    exit
                end if
            end do

            i_minus = i-1
            if (i_minus.eq.0) then
                i_minus = cell_num_nodes
            end if
            i_plus = i+1
            if (i_plus.eq.cell_num_nodes+1) then
                i_plus = 1
            end if
            previous_node_index = ielem(N, cell_index)%nodes_counterclockwise(i_minus)
            next_node_index = ielem(N, cell_index)%nodes_counterclockwise(i_plus)

            previous_node_position = local_nodes(previous_node_index)%positions(position_index, 1:dimensiona) &
                                   + (local_nodes(previous_node_index)%velocity(1:dimensiona) * d_t)
            next_node_position = local_nodes(next_node_index)%positions(position_index, 1:dimensiona) &
                                   + (local_nodes(next_node_index)%velocity(1:dimensiona) * d_t)

            v1(:) = current_node_position(:) - previous_node_position(:)
		    v2(:) = next_node_position(:) - current_node_position(:)

            len1 = zero
            len2 = zero
            do j = 1, dimensiona
                len1 = len1 + (v1(j)*v1(j))
                len2 = len2 + (v2(j)*v2(j))
            end do
            len1 = sqrt(len1)
            len2 = sqrt(len2)

            v1_cross_v2 = cross_product_2D(v1, v2)

            if (v1_cross_v2.lt.zero) then ! concave angle
            ! if (v1_cross_v2.lt.-0.175*len1*len2) then ! too concave angle

                if ((local_nodes(node_index)%relaxation_velocity(1).ne.zero).or.(local_nodes(node_index)%relaxation_velocity(2).ne.zero)) then
                    print*,"two concave angles before communication"
                end if

                if (cell_num_nodes.eq.3) then
                    print*,"inverted triangle not trying to fix"
                else
                    if (dimensiona.eq.3) then
                        print*,"trying to fix a concave angle in a quadrilateral in cell", cell_index, ielem(n,cell_index)%xxc, ielem(n,cell_index)%yyc, ielem(n,cell_index)%zzc
                    else
                        print*,"trying to fix a concave angle in a quadrilateral in cell", cell_index, ielem(n,cell_index)%xxc, ielem(n,cell_index)%yyc
                    end if
                    desired_node_position(:) = previous_node_position(:) + (len2/(len1+len2))*(next_node_position(:)-previous_node_position(:))

                    local_nodes(node_index)%relaxation_velocity(1:dimensiona) = (desired_node_position(:)-current_node_position(:))/d_t
                end if
            end if
        end do
    end do
    !$omp end do

    !$omp barrier

    !$omp do
    do i = 1, my_num_interface_nodes 
        node_index = local_interface_nodes(i)
        ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
        do cpu_index = 1,local_nodes(node_index)%num_cpus
            cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
            index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * dimensiona
            do j = 1, dimensiona
                node_snd_buffer(cpu)%data(index+j) = local_nodes(node_index)%relaxation_velocity(j)
            end do
            do j = dimensiona+1, local_nodes(node_index)%num_local_neighbours*dimensiona
                node_snd_buffer(cpu)%data(index+j) = zero
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
                    count = node_snd_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * dimensiona
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * dimensiona
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

        do i = 1, local_nodes(node_index)%num_cpus
            cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
            index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*dimensiona

            if ((local_nodes(node_index)%relaxation_velocity(1).eq.zero).and.(local_nodes(node_index)%relaxation_velocity(2).eq.zero)) then
                do j = 1, dimensiona
                    local_nodes(node_index)%relaxation_velocity(j) = node_rcv_buffer(cpu_index)%data(index+j)
                end do
            else
                if ((node_rcv_buffer(cpu_index)%data(index+1).ne.zero).or.(node_rcv_buffer(cpu_index)%data(index+2).ne.zero)) then
                    print*,"multiple convex angles after communication"
                end if
            end if
        end do

        local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) + local_nodes(node_index)%relaxation_velocity(1:dimensiona)
    end do
    !$omp end do

end subroutine fix_moved_concave_cells





! subroutine enforce_node_velocity_BC(position_index, N)
!     implicit none
!     integer,intent(in)::position_index, n
!     real::x, y, epsilon, v2
!     integer:: node_index
!     epsilon = 0.000001

!     if (initcond.eq.101) then
!         !$omp do
!         do node_index = 1, kmaxn
!             x = local_nodes(node_index)%positions(position_index,1)
!             y = local_nodes(node_index)%positions(position_index,2)
!             v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
!             if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
!                 print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
!             end if
!             ! if (v2.gt.((5.0*mesh_velocity_multiple)*(5.0*mesh_velocity_multiple))) then
!             !     print *, "too high node speed^2", v2, "in node", node_index, N
!             ! end if
!             if ((y.lt.epsilon).or.(y.gt.(1.0-epsilon))) then
!                 local_nodes(node_index)%velocity(2) = 0.0
!             end if
!             if (y.lt.0.0) then
!                 local_nodes(node_index)%velocity(2) = (-1.0)*y/dt
!             end if
!             if (y.gt.1.0) then
!                 local_nodes(node_index)%velocity(2) = (y - 1.0)/dt
!             end if
!             if (moving_mesh_mode.ge.4) then
!                 if ((x.le.-4.5+epsilon).or.(x.ge.4.5-epsilon)) then
!                     local_nodes(node_index)%velocity(1) = 0.0
!                 end if
!             end if

!         end do
!         !$omp end do
    
!     else if (initcond.eq.102) then
!         !$omp do
!         do node_index = 1, kmaxn
!             x = local_nodes(node_index)%positions(position_index,1)
!             y = local_nodes(node_index)%positions(position_index,2)
!             v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
!             if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
!                 print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
!             end if
!             ! if (v2.gt.((10.0*mesh_velocity_multiple)*(10.0*mesh_velocity_multiple))) then
!             !     print *, "too high node speed^2", v2, "in node", node_index, N
!             ! end if
!             if (moving_mesh_mode.ge.4) then
!                 if ((x.le.epsilon).or.(x.ge.4.0-epsilon)) then
!                     local_nodes(node_index)%velocity(1) = 0.0
!                 end if
!                 if ((y.le.epsilon).or.(y.ge.1.0-epsilon)) then
!                     local_nodes(node_index)%velocity(2) = 0.0
!                 end if
!             end if

!         end do
!         !$omp end do

!     else if (initcond.eq.470) then
!         !$omp do
!         do node_index = 1, kmaxn
!             x = local_nodes(node_index)%positions(position_index,1)
!             y = local_nodes(node_index)%positions(position_index,2)
!             v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
!             if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
!                 print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
!             end if
!             ! if (v2.gt.((10.0*mesh_velocity_multiple)*(10.0*mesh_velocity_multiple))) then
!             !     print *, "too high node speed^2", v2, "in node", node_index, N
!             ! end if
!             if (moving_mesh_mode.ge.4) then
!                 if ((x.le.epsilon).or.(x.ge.7.0-epsilon)) then
!                     local_nodes(node_index)%velocity(1) = 0.0
!                 end if
!                 if ((y.le.epsilon).or.(y.ge.3.0-epsilon)) then
!                     local_nodes(node_index)%velocity(2) = 0.0
!                 end if
!             end if

!         end do
!         !$omp end do
    
!     else if (initcond.eq.405) then
!         !$omp do
!         do node_index = 1, kmaxn
!             x = local_nodes(node_index)%positions(position_index,1)
!             y = local_nodes(node_index)%positions(position_index,2)
!             v2 = (local_nodes(node_index)%velocity(1)*local_nodes(node_index)%velocity(1)) + (local_nodes(node_index)%velocity(2)*local_nodes(node_index)%velocity(2))
!             if ((local_nodes(node_index)%velocity(1) .ne. local_nodes(node_index)%velocity(1)) .or. (local_nodes(node_index)%velocity(2) .ne. local_nodes(node_index)%velocity(2))) then
!                 print *, "NaN node velocity during BC check in node", node_index, local_nodes(node_index)%positions(position_index,:), N
!             end if
!             ! if (v2.gt.((10.0*mesh_velocity_multiple)*(10.0*mesh_velocity_multiple))) then
!             !     print *, "too high node speed^2", v2, "in node", node_index, N
!             ! end if
!             if (moving_mesh_mode.ge.4) then
!                 if ((x.le.-0.25+epsilon).or.(x.ge.0.25-epsilon)) then
!                     local_nodes(node_index)%velocity(1) = 0.0
!                 end if
!                 if ((y.le.epsilon).or.(y.ge.0.1-epsilon)) then
!                     local_nodes(node_index)%velocity(2) = 0.0
!                 end if
!             end if

!         end do
!         !$omp end do
!     end if

! end subroutine enforce_node_velocity_BC






! subroutine enforce_node_velocity_BC(position_index, N)
!     implicit none
!     integer,intent(in)::position_index, n
!     integer::i, ii, j, cell_index, edge_index, boundary_index, node_index, node_index_1, node_index_2
!     integer::apply, done
!     real,dimension(1:dimensiona)::edge, normalized
!     real::edge_len, dot, speed, y1, y2

!     !$omp master
!         DO II = 1, NOF_BOUNDED
!             cell_index = EL_BND(II)
!             do edge_index = 1, IELEM(N,cell_index)%ifca
                
!                 if (IELEM(N,cell_index)%INEIGHB(edge_index).ne.N) then
!                     cycle ! it is an interface between CPUs not a domain boundary
!                 end if
!                 if (ielem(n,cell_index)%ibounds(edge_index).eq.0) then
!                     cycle ! not a domian boundary
!                 end if
!                 node_index_1 = IELEM(N,cell_index)%nodes_faces(edge_index, 1)
!                 node_index_2 = IELEM(N,cell_index)%nodes_faces(edge_index, 2)
                
!                 apply = 0
!                 if ((moving_mesh_mode.gt.2).and.(relaxation_mesh_velocity_multiple.ne.zero)) then
!                     apply = 1
!                 else
!                     if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.3) then ! symmetry BC
!                         apply = 1
!                     end if
!                     if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.4) then ! wall BC
!                         apply = 1
!                     end if
!                 end if
!                 ! if (initcond.eq.102) then
!                 !     apply = 1
!                 ! end if
!                 if (apply.gt.0) then
!                     if (dimensiona.eq.3) Then
!                         print*, "moving mesh currently does not support 3D boundary conditions"
!                     else
!                         edge(:) = local_nodes(node_index_2)%positions(position_index,:) - local_nodes(node_index_1)%positions(position_index,:)

!                         done = 0
!                         if (initcond.eq.105) Then
!                             local_nodes(node_index_1)%velocity(:) = zero
!                             local_nodes(node_index_2)%velocity(:) = zero
!                         else
!                             if (edge(2).ne.zero) then
!                                 local_nodes(node_index_1)%velocity(1) = zero
!                                 local_nodes(node_index_2)%velocity(1) = zero
!                                 done = 1
!                             end if
!                             if (edge(1).ne.zero) then
!                                 local_nodes(node_index_1)%velocity(2) = zero
!                                 local_nodes(node_index_2)%velocity(2) = zero
!                                 done = 1
!                             end if
!                             if (done.eq.0) then
!                                 edge_len = sqrt((edge(1)*edge(1)) + (edge(2)*edge(2)))
!                                 ! speed = local_nodes(node_index_1)%velocity(1)*local_nodes(node_index_1)%velocity(1)
!                                 ! speed = speed + (local_nodes(node_index_1)%velocity(2)*local_nodes(node_index_1)%velocity(2))
!                                 ! speed = sqrt(speed)
!                                 dot = local_nodes(node_index_1)%velocity(1) * edge(1)
!                                 dot = dot + local_nodes(node_index_1)%velocity(2) * edge(2)
!                                 dot = dot / edge_len
!                                 local_nodes(node_index_1)%velocity(1) = edge(1) * dot / edge_len
!                                 local_nodes(node_index_1)%velocity(2) = edge(2) * dot / edge_len
!                                 local_nodes(node_index_1)%velocity(3) = zero

!                                 ! speed = local_nodes(node_index_2)%velocity(1)*local_nodes(node_index_2)%velocity(1)
!                                 ! speed = speed + (local_nodes(node_index_2)%velocity(2)*local_nodes(node_index_2)%velocity(2))
!                                 ! speed = sqrt(speed)
!                                 dot = local_nodes(node_index_2)%velocity(1) * edge(1)
!                                 dot = dot + local_nodes(node_index_2)%velocity(2) * edge(2)
!                                 dot = dot / edge_len
!                                 local_nodes(node_index_2)%velocity(1) = edge(1) * dot / edge_len
!                                 local_nodes(node_index_2)%velocity(2) = edge(2) * dot / edge_len
!                                 local_nodes(node_index_2)%velocity(3) = zero
!                             end if
!                         end if
!                     end if
!                 end if
!                 if (initcond.eq.102) then
!                     if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.1) then ! inflow
!                         y1 = local_nodes(node_index_1)%positions(position_index,2)
!                         y2 = local_nodes(node_index_2)%positions(position_index,2)
!                         if ((y1.eq.zero).and.(y2.eq.zero)) then
!                             local_nodes(node_index_1)%velocity(1) = zero
!                             local_nodes(node_index_2)%velocity(1) = zero
!                         end if
!                     end if
!                 end if
!                 if (BOUNDARY_MOVEMENT) then
!                     if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.gt.100) then ! moving wall
!                         boundary_index = ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode - 100
!                         local_nodes(node_index_1)%velocity(:) = zero
!                         local_nodes(node_index_2)%velocity(:) = zero
!                         if (((initcond.eq.105).and.(t.ge.0.01)).or.(initcond.ne.105)) then
!                             local_nodes(node_index_1)%velocity(1:dimensiona) = boundary_velocity(boundary_index, 1:dimensiona)
!                             local_nodes(node_index_2)%velocity(1:dimensiona) = boundary_velocity(boundary_index, 1:dimensiona)
!                         end if
!                     end if
!                 end if

!             end do 
!         end do
!     !$omp end master

! end subroutine enforce_node_velocity_BC





! subroutine enforce_node_velocity_BC(position_index, d_t, N)
!     implicit none
!     integer,intent(in)::position_index, n
!     real,intent(in)::d_t
!     integer::i, ii, j, cell_index, edge_index, boundary_index, node_index, node_index_1, node_index_2
!     integer::apply, done
!     real,dimension(1:dimensiona)::edge, normalized
!     real::edge_len, edge_len2, dot, speed, x, y

!     !$omp do
!     do ii = 1, my_num_boundary_nodes
!         node_index = local_boundary_nodes(ii)

!         ! if (local_nodes(node_index)%communication.gt.0) then
!         !     ! to ensure consistency between CPUs set (non-moving) boundary velocity to zero
!         !     local_nodes(node_index)%velocity(:) = zero
!         ! else
!         if ((moving_mesh_mode.gt.2).and.(local_nodes(node_index)%num_neighbours.eq.1)) then
!             local_nodes(node_index)%velocity(:) = zero
!         else
!             do i = 1, local_nodes(node_index)%num_local_neighbours
!                 cell_index = local_nodes(node_index)%local_neighbours(i)

!                 if (IELEM(N,cell_index)%INTERIOR.EQ.0) then
!                     ! no boundary conditions in this cell
!                     cycle
!                 end if

!                 do edge_index = 1, ielem(N, cell_index)%ifca
!                     if (IELEM(N, cell_index)%INEIGHB(edge_index).ne.N) then
!                         cycle ! it is an interface between CPUs not a domain boundary
!                     end if
!                     if (ielem(n,cell_index)%ibounds(edge_index).eq.0) then
!                         cycle ! not a domian boundary
!                     end if

!                     if ((IELEM(N,cell_index)%nodes_faces(edge_index, 1).eq.node_index).or.(IELEM(N,cell_index)%nodes_faces(edge_index, 2).eq.node_index)) then
!                         ! this node belongs to this edge
!                         apply = 0
!                         if (moving_mesh_mode.gt.2) then
!                             apply = 1
!                         else
!                             if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.3) then ! symmetry BC
!                                 apply = 1
!                             end if
!                             if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.4) then ! wall BC
!                                 apply = 1
!                             end if
!                             if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.gt.100) then ! moving boundary
!                                 apply = 1
!                             end if
!                         end if

!                         if (apply.gt.0) then
!                             if (dimensiona.eq.3) Then
!                                 print*, "moving mesh currently does not support 3D boundary conditions"
!                             else
!                                 node_index_1 = IELEM(N,cell_index)%nodes_faces(edge_index, 1)
!                                 node_index_2 = IELEM(N,cell_index)%nodes_faces(edge_index, 2)
!                                 edge(:) = local_nodes(node_index_2)%positions(position_index,:) - local_nodes(node_index_1)%positions(position_index,:)

!                                 if (edge(1).eq.zero) then
!                                     local_nodes(node_index)%velocity(1) = zero
!                                 else if (edge(2).eq.zero) then
!                                     local_nodes(node_index)%velocity(2) = zero
!                                 else
!                                     if (initcond.eq.105) then
!                                         local_nodes(node_index)%velocity(:) = zero
!                                     else
!                                         edge_len2 = (edge(1)*edge(1)) + (edge(2)*edge(2))
!                                         edge_len = sqrt(edge_len2)
!                                         edge(:) = edge(:)/edge_len

!                                         dot = local_nodes(node_index)%velocity(1) * edge(1)
!                                         dot = dot + (local_nodes(node_index)%velocity(2) * edge(2))

!                                         if (dot.gt.zero) then
!                                             dot = min(dot, 0.25*(edge_len/d_t))
!                                         else
!                                             dot = max(dot, -0.25*(edge_len/d_t))
!                                         end if
!                                         ! dot = zero
                                        
!                                         local_nodes(node_index)%velocity(1) = edge(1) * dot
!                                         local_nodes(node_index)%velocity(2) = edge(2) * dot
!                                         local_nodes(node_index)%velocity(3) = zero
!                                     end if
!                                 end if
!                             end if
!                         end if
!                         if (initcond.eq.102) then
!                             if (ibound(n,ielem(n,cell_index)%ibounds(edge_index))%icode.eq.1) then ! inflow
!                                 y = local_nodes(node_index)%positions(position_index, 2)
!                                 if (y.eq.zero) then
!                                     local_nodes(node_index)%velocity(1) = zero
!                                 end if
!                             end if
!                         end if
!                     end if
!                 end do
!             end do
!         end if
!     end do
!     !$omp end do

!     !$omp barrier

!     if (BOUNDARY_MOVEMENT) then
!         !$omp do
!         do i = 1, my_num_moving_nodes
!             node_index = local_moving_nodes(i)

!             if (local_nodes(node_index)%boundary.lt.100) then
!                 print*,"something went wrong with local_nodes(node_index)%boundary"
!             end if
!             boundary_index = local_nodes(node_index)%boundary-100
!             if (.not.((initcond.eq.105).and.(t.ge.(2.0*0.7/uvel)))) then
!                 local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) + boundary_velocity(boundary_index, 1:dimensiona)
!             end if
!         end do
!         !$omp end do
!     end if

!     !$omp barrier

! end subroutine enforce_node_velocity_BC





subroutine enforce_node_velocity_BC(position_index, d_t, N)
    implicit none
    integer,intent(in)::position_index, N
    real,intent(in)::d_t

    ! integer M
    integer::node_index, node_plus_index, node_minus_index, cell_index, index, node_num_neighbours, cell_num_nodes
    integer::cpu_index, cpu
    integer::iter, i, ii, j, k, i_plus, i_minus, counter
    real,dimension(1:dimensiona)::p, p_minus, p_plus
    real::total_direction_len2, len_plus, len_minus, dot
    real,dimension(1:dimensiona,1:(2*max_num_node_neighbours))::points
    real,dimension(1:dimensiona)::direction_minus, direction_plus, total_direction

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
        ! if (local_nodes(node_index)%boundary.gt.zero) then
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

                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_Jacobi_relaxation_velocity"
                        end if
                        node_snd_buffer(cpu)%data(index) = p_minus(k)
                    end do
                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_Jacobi_relaxation_velocity"
                        end if
                        node_snd_buffer(cpu)%data(index) = p_plus(k)
                    end do
                end do
            end do
        ! end if
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
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
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
    do ii = 1, my_num_boundary_nodes
        node_index = local_boundary_nodes(ii)

        node_num_neighbours = local_nodes(node_index)%num_neighbours

        p(:) = local_nodes(node_index)%positions(position_index,:)

        counter = 0
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

                counter = counter+1
                points(:,counter) = p_minus(:)
                counter = counter+1
                points(:,counter) = p_plus(:)

                index = index + (2*dimensiona)
            end do
        end do

        if (counter.ne.(2*node_num_neighbours)) Then
            print*,"wrong number of points in enforce_node_velocity_BC"
        end if

        total_direction(:) = zero
        do i = 1, node_num_neighbours
            direction_minus = points(:,(2*i)-1) - p(:)
            direction_plus = p(:) - points(:,2*i)
            
            len_minus= zero
            len_plus = zero
            do j = 1, dimensiona
                len_minus= len_minus+ (direction_minus(j)**2)
                len_plus = len_plus + (direction_plus(j)**2)
            end do
            len_minus= sqrt(len_minus)
            len_plus = sqrt(len_plus)

            total_direction(:) = total_direction(:) + (direction_minus(:)/len_minus)
            total_direction(:) = total_direction(:) + (direction_plus(:)/len_plus)
        end do

        total_direction_len2 = zero
        do j = 1, dimensiona
            total_direction_len2 = total_direction_len2 + (total_direction(j)**2)
        end do
        
        if (total_direction_len2.lt.3.0) then ! max angle 30 degrees
            local_nodes(node_index)%velocity(:) = zero
        else
            dot = zero
            do j = 1, dimensiona
                dot = dot + (local_nodes(node_index)%velocity(j) * total_direction(j))
            end do
                            
            local_nodes(node_index)%velocity(1) = dot * total_direction(1) / total_direction_len2
            local_nodes(node_index)%velocity(2) = dot * total_direction(2) / total_direction_len2
            local_nodes(node_index)%velocity(3) = zero
        end if

        if (initcond.eq.102) then
            if (((local_nodes(node_index)%positions(position_index,1).le.0.2).or.(local_nodes(node_index)%positions(1,1).le.0.2)).and. &
                ((local_nodes(node_index)%positions(position_index,2).le.0.001).or.(local_nodes(node_index)%positions(1,2).le.0.001))) then
                local_nodes(node_index)%velocity(:) = zero
            end if
        end if

    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

end subroutine enforce_node_velocity_BC





subroutine enforce_node_lagrangian_velocity_BC(position_index, d_t, N)
    implicit none
    integer,intent(in)::position_index, N
    real,intent(in)::d_t

    ! integer M
    integer::node_index, node_plus_index, node_minus_index, cell_index, index, node_num_neighbours, cell_num_nodes
    integer::cpu_index, cpu
    integer::iter, i, ii, j, k, i_plus, i_minus, counter
    real,dimension(1:dimensiona)::p, p_minus, p_plus
    real::total_direction_len2, len_plus, len_minus, dot
    real,dimension(1:dimensiona,1:(2*max_num_node_neighbours))::points
    real,dimension(1:dimensiona)::direction_minus, direction_plus, total_direction

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
        ! if (local_nodes(node_index)%boundary.gt.zero) then
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

                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_Jacobi_relaxation_velocity"
                        end if
                        node_snd_buffer(cpu)%data(index) = p_minus(k)
                    end do
                    do k = 1,dimensiona
                        index = index +1
                        if (index.gt.node_snd_count(cpu) * (2*dimensiona)) then
                            print *, "copying too much data to send buffer from", N, "to", cpu, "find_node_Jacobi_relaxation_velocity"
                        end if
                        node_snd_buffer(cpu)%data(index) = p_plus(k)
                    end do
                end do
            end do
        ! end if
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
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * (2*dimensiona)
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
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
    do ii = 1, my_num_boundary_nodes
        node_index = local_boundary_nodes(ii)

        node_num_neighbours = local_nodes(node_index)%num_neighbours

        p(:) = local_nodes(node_index)%positions(position_index,:)

        counter = 0
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

                counter = counter+1
                points(:,counter) = p_minus(:)
                counter = counter+1
                points(:,counter) = p_plus(:)

                index = index + (2*dimensiona)
            end do
        end do

        if (counter.ne.(2*node_num_neighbours)) Then
            print*,"wrong number of points in enforce_node_lagrangian_velocity_BC"
        end if

        total_direction(:) = zero
        do i = 1, node_num_neighbours
            direction_minus = points(:,(2*i)-1) - p(:)
            direction_plus = p(:) - points(:,2*i)
            
            len_minus= zero
            len_plus = zero
            do j = 1, dimensiona
                len_minus= len_minus+ (direction_minus(j)**2)
                len_plus = len_plus + (direction_plus(j)**2)
            end do
            len_minus= sqrt(len_minus)
            len_plus = sqrt(len_plus)

            total_direction(:) = total_direction(:) + (direction_minus(:)/len_minus)
            total_direction(:) = total_direction(:) + (direction_plus(:)/len_plus)
        end do

        total_direction_len2 = zero
        do j = 1, dimensiona
            total_direction_len2 = total_direction_len2 + (total_direction(j)**2)
        end do
        
        if (total_direction_len2.lt.3.0) then ! max angle 30 degrees
            local_nodes(node_index)%lagrangian_velocity(:) = zero
        else
            dot = zero
            do j = 1, dimensiona
                dot = dot + (local_nodes(node_index)%lagrangian_velocity(j) * total_direction(j))
            end do
                            
            local_nodes(node_index)%lagrangian_velocity(1) = dot * total_direction(1) / total_direction_len2
            local_nodes(node_index)%lagrangian_velocity(2) = dot * total_direction(2) / total_direction_len2
            local_nodes(node_index)%lagrangian_velocity(3) = zero
        end if

        if (initcond.eq.102) then
            if (((local_nodes(node_index)%positions(position_index,1).le.0.2).or.(local_nodes(node_index)%positions(1,1).le.0.2)).and. &
                ((local_nodes(node_index)%positions(position_index,2).le.0.001).or.(local_nodes(node_index)%positions(1,2).le.0.001))) then
                local_nodes(node_index)%lagrangian_velocity(:) = zero
            end if
        end if

    end do
    !$omp end do
    
    !$omp barrier

    !$omp master
        CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master

    !$omp barrier

end subroutine enforce_node_lagrangian_velocity_BC





subroutine clamp_node_velocity_to_CFL(position_index, d_t, N)
    implicit none
    integer,intent(in)::position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, node_index, cell_index, cpu_index, cpu, index, num_faces
    ! integer:: M
    real::volume, area, lengthscale, min_lengthscale, velocity_magnitude, max_velocity

    real, parameter::fraction = 0.5

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.1) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1, my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            do cpu_index = 1, local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * 1
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    volume = ielem(n, cell_index)%moving_volume(position_index)
                    if (dimensiona.eq.3) then
                        print*,"3D not implemented yet"
                    else
                        if (ielem(n, cell_index)%ishape.eq.5) then
                            num_faces = 4
                        else
                            num_faces = 3
                        end if
                    end if
                    area = zero
                    do i = 1, num_faces
                        area = area + ielem(n, cell_index)%surf(i)
                    end do
                    if (dimensiona.eq.3) then
                        lengthscale = 6.0 * volume / area
                    else
                        lengthscale = 4.0 * volume / area
                    end if

                    index = index + 1
                    node_snd_buffer(cpu)%data(index) = lengthscale
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
                    count = node_snd_count(cpu_index) * 1
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * 1
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * 1
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

            min_lengthscale = 10000000000.0

            do iter = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(iter)

                volume = ielem(n, cell_index)%moving_volume(position_index)
                if (dimensiona.eq.3) then
                    print*,"3D not implemented yet"
                else
                    if (ielem(n, cell_index)%ishape.eq.5) then
                        num_faces = 4
                    else
                        num_faces = 3
                    end if
                end if
                area = zero
                do i = 1, num_faces
                    area = area + ielem(n, cell_index)%surf(i)
                end do
                if (dimensiona.eq.3) then
                    lengthscale = 6.0 * volume / area
                else
                    lengthscale = 4.0 * volume / area
                end if

                if (lengthscale.lt.min_lengthscale) then
                    min_lengthscale = lengthscale
                end if

            end do

            do iter = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(iter)%cpu
                index = (local_nodes(node_index)%rcv_offsets(iter)%lower - 1) * 1
                do j = local_nodes(node_index)%rcv_offsets(iter)%lower, local_nodes(node_index)%rcv_offsets(iter)%upper
                    index = index+1
                    lengthscale = node_rcv_buffer(cpu_index)%data(index)
                    if (lengthscale.lt.min_lengthscale) then
                        min_lengthscale = lengthscale
                    end if
                end do
            end do

            max_velocity = fraction * CFL * min_lengthscale/d_t

            velocity_magnitude = zero
            do i = 1, dimensiona
                velocity_magnitude = velocity_magnitude + (local_nodes(node_index)%velocity(i)**2)
            end do
            velocity_magnitude = sqrt(velocity_magnitude)

            if (velocity_magnitude.gt.max_velocity) then
                local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) * (max_velocity/velocity_magnitude)
                if (dimensiona.eq.3) then
                    print*,"clamping velocity in node", node_index, local_nodes(node_index)%positions(position_index,1), local_nodes(node_index)%positions(position_index,2), local_nodes(node_index)%positions(position_index,3)
                else
                    print*,"clamping velocity in node", node_index, local_nodes(node_index)%positions(position_index,1), local_nodes(node_index)%positions(position_index,2)
                end if    
            end if

        end do
    !$omp end do
    
    !$omp barrier
    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master
    !$omp barrier

end subroutine clamp_node_velocity_to_CFL





subroutine clamp_node_velocity_to_fraction(position_index, d_t, N)
    implicit none
    integer,intent(in)::position_index, N
    real,intent(in)::d_t
    integer::i, j, k, iter, node_index, cell_index, cpu_index, cpu, index, num_faces
    ! integer:: M
    real::volume, area, lengthscale, min_lengthscale, velocity_magnitude, max_velocity

    real, parameter::fraction = 0.75

    integer,dimension(2*isize)::requests
    integer::num_requests, count

    ! M = omp_get_thread_num()

    if (num_values_to_send_per_node.lt.1) then
        print *,"something went wrong sorry :("
        call abort
    end if

    !$omp do
        do iter = 1, my_num_interface_nodes 
            node_index = local_interface_nodes(iter)
            do cpu_index = 1, local_nodes(node_index)%num_cpus
                cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
                index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * 1
                do j = 1, local_nodes(node_index)%num_local_neighbours
                    cell_index = local_nodes(node_index)%local_neighbours(j)
                    volume = ielem(n, cell_index)%moving_volume(position_index)
                    if (dimensiona.eq.3) then
                        print*,"3D not implemented yet"
                    else
                        if (ielem(n, cell_index)%ishape.eq.5) then
                            num_faces = 4
                        else
                            num_faces = 3
                        end if
                    end if
                    area = zero
                    do i = 1, num_faces
                        area = area + ielem(n, cell_index)%surf(i)
                    end do
                    if (dimensiona.eq.3) then
                        lengthscale = 6.0 * volume / area
                    else
                        lengthscale = 4.0 * volume / area
                    end if

                    index = index + 1
                    node_snd_buffer(cpu)%data(index) = lengthscale
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
                    count = node_snd_count(cpu_index) * 1
                    num_requests = num_requests + 1
                    CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, "sending from", N, "to", cpu_index, count, "values"
                end if
                if (node_rcv_count(cpu_index).gt.0) then
                    count = node_rcv_count(cpu_index) * 1
                    num_requests = num_requests + 1
                    CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
                    ! print *, N, "waiting to receive", count, "values from", cpu_index
                end if
            else
                if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
                    print *,"send receive count missmatch on CPU", n
                    call abort
                end if
                count = node_snd_count(cpu_index) * 1
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

            min_lengthscale = 10000000000.0

            do iter = 1, local_nodes(node_index)%num_local_neighbours
                cell_index = local_nodes(node_index)%local_neighbours(iter)

                volume = ielem(n, cell_index)%moving_volume(position_index)
                if (dimensiona.eq.3) then
                    print*,"3D not implemented yet"
                else
                    if (ielem(n, cell_index)%ishape.eq.5) then
                        num_faces = 4
                    else
                        num_faces = 3
                    end if
                end if
                area = zero
                do i = 1, num_faces
                    area = area + ielem(n, cell_index)%surf(i)
                end do
                if (dimensiona.eq.3) then
                    lengthscale = 6.0 * volume / area
                else
                    lengthscale = 4.0 * volume / area
                end if

                if (lengthscale.lt.min_lengthscale) then
                    min_lengthscale = lengthscale
                end if

            end do

            do iter = 1, local_nodes(node_index)%num_cpus
                cpu_index = local_nodes(node_index)%rcv_offsets(iter)%cpu
                index = (local_nodes(node_index)%rcv_offsets(iter)%lower - 1) * 1
                do j = local_nodes(node_index)%rcv_offsets(iter)%lower, local_nodes(node_index)%rcv_offsets(iter)%upper
                    index = index+1
                    lengthscale = node_rcv_buffer(cpu_index)%data(index)
                    if (lengthscale.lt.min_lengthscale) then
                        min_lengthscale = lengthscale
                    end if
                end do
            end do

            max_velocity = fraction * min_lengthscale/d_t

            velocity_magnitude = zero
            do i = 1, dimensiona
                velocity_magnitude = velocity_magnitude + (local_nodes(node_index)%velocity(i)**2)
            end do
            velocity_magnitude = sqrt(velocity_magnitude)

            if (velocity_magnitude.gt.max_velocity) then
                local_nodes(node_index)%velocity(1:dimensiona) = local_nodes(node_index)%velocity(1:dimensiona) * (max_velocity/velocity_magnitude)
                if (dimensiona.eq.3) then
                    print*,"clamping velocity in node", node_index, local_nodes(node_index)%positions(position_index,1), local_nodes(node_index)%positions(position_index,2), local_nodes(node_index)%positions(position_index,3)
                else
                    print*,"clamping velocity in node", node_index, local_nodes(node_index)%positions(position_index,1), local_nodes(node_index)%positions(position_index,2)
                end if    
            end if

        end do
    !$omp end do
    
    !$omp barrier
    !$omp master
        call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
    !$omp end master
    !$omp barrier

end subroutine clamp_node_velocity_to_fraction





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





! SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_2D(stage, position_index, N)
!     IMPLICIT NONE
!     integer::stage, position_index, N
!     integer::cell_index, iter
!     integer::my_deg_free
!     real,dimension(1:((IORDER+1)*(IORDER+2))/2 - 1)::gradient
!     real::gradient_magnitude, my_max_gradient_magnitude

!     my_deg_free = ((IORDER+1)*(IORDER+2))/2 - 1

!     !$omp master
!         CALL EXCH_CORDS_MovingMesh(N, 1)
!     !$omp end master
!     IF (FASTEST.EQ.1) THEN
!         CALL EXCHANGE_LOWER(N)
!     ELSE
!         CALL EXCHANGE_HIGHER(N)
!     END IF

!         ! CALL EXCHANGE_LOWER(N)
!         ! CALL EXCHANGE_HIGHER(N)
    
!     CALL LEAST_SQUARES(N)

!     my_max_gradient_magnitude = zero

!     !$omp parallel do reduction(max: max_gradient_magnitude)
!     DO cell_index = 1, XMPIELRANK(N)
!         gradient(1:my_deg_free) = ILOCAL_rECON5(cell_index)%GRADIENTS(1,1:my_deg_free,1)
!         gradient_magnitude = zero
!         do iter = 1,dimensiona
!             gradient_magnitude = gradient_magnitude + (gradient(iter)*gradient(iter))
!         end do
!         gradient_magnitude = sqrt(gradient_magnitude)

!         u_c(cell_index)%normalized_gradient = gradient_magnitude

!         if (my_max_gradient_magnitude.lt.gradient_magnitude) then
!             my_max_gradient_magnitude = gradient_magnitude
!         end if
!     END DO
!     !$omp end parallel do

!     ! print *, "first loop done", N, max_gradient_magnitude, my_max_gradient_magnitude

!     !$omp barrier

!     !$omp master
!        ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!         CALL MPI_ALLREDUCE(my_max_gradient_magnitude, max_gradient_magnitude, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
!         ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!     !$omp end master

!     !$omp barrier
!     my_max_gradient_magnitude = max_gradient_magnitude

!     ! print *, "reduction done", N, max_gradient_magnitude

!     !$omp do
!     DO cell_index = 1, XMPIELRANK(N)
!         u_c(cell_index)%normalized_gradient = u_c(cell_index)%normalized_gradient / my_max_gradient_magnitude
!     END DO
!     !$omp end do

! END SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_2D





! function linear_solve(n, matrix, rhs, solution)
!     implicit none
!     integer,intent(in)::N
!     real,intent(inout), dimension(1:n,1:n)::matrix
!     real,intent(inout), dimension(1:n)::rhs
!     real,intent(out),dimension(1:n)::solution
!     logical::linear_solve
!     real::swap_helper,factor
!     integer::i,j,k

!     linear_solve = .true.

!     do i = 1, N-1
!         do j = (i+1), N
!             if (matrix(j,i).ne.zero) then
!                 if (matrix(i,i).eq.zero) then
!                     do k = i, N
!                         swap_helper = matrix(i,k)
!                         matrix(i,k) = matrix(j,k)
!                         matrix(j,k) = swap_helper
!                     end do

!                     swap_helper = rhs(i)
!                     rhs(i) = rhs(j)
!                     rhs(j) = swap_helper
!                 else
!                     factor = matrix(j,i)/matrix(i,i) 
!                     matrix(j,:) = matrix(j,:) - (factor*matrix(i,:))
!                     rhs(j) = rhs(j) - (factor*rhs(i))
!                 end if
!             end if
!         end do
!     end do

!     do i = N, 1, -1
!         solution(i) = rhs(i)
!         do j = i+1, N
!             solution(i) = solution(i) - (matrix(i,j)*solution(j))
!         end do
!         if (matrix(i,i).eq.zero) then
!             if (solution(i).ne.zero) then
!                 print*,"linear solver failure - dividing by zero at (", i, i, ")"
!                 ! call abort()
!                 linear_solve = .false.
!             end if
!         else
!             solution(i) = solution(i)/matrix(i,i)
!         end if
!     end do

! end function linear_solve 





function linear_solve(n, matrix, rhs, solution)
    implicit none
    integer,intent(in)::N
    real,intent(inout), dimension(1:n,1:n)::matrix
    real,intent(inout), dimension(1:n)::rhs
    real,intent(out),dimension(1:n)::solution
    real, dimension(1:n,1:n)::inverse
    logical::linear_solve
    real::swap_helper,factor
    real::det
    integer::i,j,k

    linear_solve = .false.

    solution(:) = zero

    if (n.ne.3) then
        print *, "lineat solver can only handle 3x3 matrices"
        call abort
    end if

    det = (matrix(1,1)*matrix(2,2)*matrix(3,3)) &
        + (matrix(1,2)*matrix(2,3)*matrix(3,1)) &
        + (matrix(1,3)*matrix(2,1)*matrix(3,2)) &
        - (matrix(1,3)*matrix(2,2)*matrix(3,1)) &
        - (matrix(1,2)*matrix(2,1)*matrix(3,3)) &
        - (matrix(1,1)*matrix(2,3)*matrix(3,2))
    
    if (det.ne.zero) then
        linear_solve = .true.
        inverse(1,1) = ((matrix(2,2)*matrix(3,3))-(matrix(2,3)*matrix(3,2)))/det
        inverse(2,2) = ((matrix(1,1)*matrix(3,3))-(matrix(1,3)*matrix(3,1)))/det
        inverse(3,3) = ((matrix(1,1)*matrix(2,2))-(matrix(1,2)*matrix(2,1)))/det
        inverse(1,2) = (-1.0)*((matrix(1,2)*matrix(3,3))-(matrix(1,3)*matrix(3,2)))/det
        inverse(2,1) = (-1.0)*((matrix(2,1)*matrix(3,3))-(matrix(2,3)*matrix(3,1)))/det
        inverse(2,3) = (-1.0)*((matrix(1,1)*matrix(2,3))-(matrix(1,3)*matrix(2,1)))/det
        inverse(3,2) = (-1.0)*((matrix(1,1)*matrix(3,2))-(matrix(1,2)*matrix(3,1)))/det
        inverse(1,3) = ((matrix(1,2)*matrix(2,3))-(matrix(1,3)*matrix(2,2)))/det
        inverse(3,1) = ((matrix(2,1)*matrix(3,2))-(matrix(2,2)*matrix(3,1)))/det

        do i = 1, n
            do j = 1, n
                solution(i) = solution(i) + inverse(i,j)*rhs(j)
            end do
        end do
    else
        ! print*,"linear solver failure - zero determinant"
    end if

end function linear_solve 





! SUBROUTINE find_node_density_gradient(stage, position_index, N)
!     implicit none
!     integer,intent(in)::stage, position_index, N
!     ! real,intent(in)::d_t
!     integer::i, j, k, iter, counter, node_index, cell_index, cpu_index, cpu, index, rho_index
!     integer::num_neighbours
!     ! integer:: M
!     real::rho, v
!     real,dimension(1:nof_variables)::copy
!     real::dummuy_MP_PINFl, dummy_gammal
!     real,dimension(1:dimensiona)::helper_centre_position, node_center
!     real,dimension(1:max_num_node_neighbours,1:dimensiona)::centre_positions
!     real,dimension(1:max_num_node_neighbours)::rho_vector
!     real,dimension(1:(dimensiona+1),1:max_num_node_neighbours)::At
!     real,dimension(1:max_num_node_neighbours,1:(dimensiona+1))::A
!     real,dimension(1:(dimensiona+1),1:(dimensiona+1))::AtA
!     real,dimension(1:(dimensiona+1))::At_rho_vector, x
!     real::swap_helper, factor
!     logical::solvable

!     integer,dimension(2*isize)::requests
!     integer::num_requests, count

!     ! M = omp_get_thread_num()
    
!     rho_index = 1

!     if (num_values_to_send_per_node.lt.(dimensiona+1)) then
!         print *,"something went wrong sorry :("
!         call abort
!     end if

!     ! do cpu_index = 0, isize-1
!     !     index = 0
!     !     do i = 1, node_snd_count(cpu_index)
!     !         do j = 1, dimensiona
!     !             index = index +1
!     !             node_snd_buffer(cpu_index)%data(index) = 1000.0
!     !         end do
!     !     end do
!     !     if (index.ne.dimensiona*node_snd_count(cpu_index)) then
!     !         print *, "something went wrong with rezeroing the send buffer"
!     !     else
!     !         print *, "snd_buffer on CPU", N, "thread", M, "to", cpu_index, "reset with", node_snd_count(cpu_index)*dimensiona, "values"
!     !     end if
!     ! end do
!     ! !$omp barrier 

!     !$omp do
!         do iter = 1,my_num_interface_nodes 
!             node_index = local_interface_nodes(iter)
!             ! print *, "on CPU", N, "thread", M, "coping data of", node_index, "(", iter, ") to send buffer" 
!             do cpu_index = 1,local_nodes(node_index)%num_cpus
!                 cpu = local_nodes(node_index)%snd_offsets(cpu_index)%cpu
!                 index = (local_nodes(node_index)%snd_offsets(cpu_index)%lower -1) * num_values_to_send_per_node
!                 do j = 1, local_nodes(node_index)%num_local_neighbours
!                     cell_index = local_nodes(node_index)%local_neighbours(j)
!                     rho = u_c(cell_index)%val(stage,rho_index)
!                     helper_centre_position(1) = ielem(N, cell_index)%xxc
!                     helper_centre_position(2) = ielem(N, cell_index)%yyc
!                     if (dimensiona.eq.3) then
!                         helper_centre_position(3) = ielem(N, cell_index)%zzc
!                     end if
!                     do k = 1,dimensiona
!                         if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
!                             print *, "copying too much data to send buffer from", N, "to", cpu 
!                         end if
!                         node_snd_buffer(cpu)%data(index+k) = helper_centre_position(k)
!                     end do
!                     k = dimensiona+1
!                     if ((index+k).gt.node_snd_count(cpu) * num_values_to_send_per_node) then
!                         print *, "copying too much data to send buffer from", N, "to", cpu 
!                     end if
!                     node_snd_buffer(cpu)%data(index+k) = rho
!                     index = index + num_values_to_send_per_node
!                 end do
!             end do
!         end do
!     !$omp end do

!     ! do cpu_index = 0, isize-1
!     !     index = 0
!     !     do i = 1, node_rcv_count(cpu_index)
!     !         do j = 1, dimensiona
!     !             index = index +1
!     !             node_rcv_buffer(cpu_index)%data(index) = 1000000.0
!     !         end do
!     !     end do
!     !     if (index.ne.dimensiona*node_rcv_count(cpu_index)) then
!     !         print *, "something went wrong with rezeroing the receive buffer"
!     !     else
!     !         print *, "rcv_buffer on CPU", N, "thread", M, "from", cpu_index, "reset with", node_rcv_count(cpu_index)*dimensiona, "values"
!     !     end if
!     ! end do

!     !$omp barrier

!     num_requests = 0
!     !$omp master
!         ! print *, "inside send_rcv part on CPU", N
!         do cpu_index = 0, isize-1
!             if (cpu_index.ne.N) then
!                 if (node_snd_count(cpu_index).gt.0) then
!                     count = node_snd_count(cpu_index) * num_values_to_send_per_node
!                     num_requests = num_requests + 1
!                     CALL MPI_ISEND(node_snd_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
!                     ! print *, "sending from", N, "to", cpu_index, count, "values"
!                 end if
!                 if (node_rcv_count(cpu_index).gt.0) then
!                     count = node_rcv_count(cpu_index) * num_values_to_send_per_node
!                     num_requests = num_requests + 1
!                     CALL MPI_IRECV(node_rcv_buffer(cpu_index)%data(1:count), count, MPI_DOUBLE_PRECISION, cpu_index, 9+n+cpu_index, MPI_COMM_WORLD, requests(num_requests), IERROR)
!                     ! print *, N, "waiting to receive", count, "values from", cpu_index
!                 end if
!             else
!                 if (node_snd_count(cpu_index).ne.node_rcv_count(cpu_index)) then
!                     print *,"send receive count missmatch on CPU", n
!                     call abort
!                 end if
!                 count = node_snd_count(cpu_index) * num_values_to_send_per_node
!                 do i = 1, count
!                     node_rcv_buffer(cpu_index)%data(i) = node_snd_buffer(cpu_index)%data(i)
!                 end do
!             end if
!         end do

!         ! print*,"CPU", n, "witing on", num_requests, "requests"
!         CALL MPI_WAITALL(num_requests, requests, MPI_STATUSES_IGNORE, IERROR)

!     !$omp end master

!     !$omp barrier

!     !$omp do
!         do node_index = 1, kmaxn 

!             local_nodes(node_index)%relaxation_velocity(:) = 0.0
!             node_center(1:dimensiona) = local_nodes(node_index)%positions(position_index,1:dimensiona)
!             num_neighbours = local_nodes(node_index)%num_neighbours
!             counter = 0

!             do i = 1, local_nodes(node_index)%num_local_neighbours
!                 counter = counter + 1
!                 cell_index = local_nodes(node_index)%local_neighbours(i)

!                 helper_centre_position(1) = ielem(N, cell_index)%xxc
!                 helper_centre_position(2) = ielem(N, cell_index)%yyc
!                 if (dimensiona.eq.3) then
!                     helper_centre_position(3) = ielem(N, cell_index)%zzc
!                 end if

!                 centre_positions(counter, :) = helper_centre_position(:)
                
!                 A(counter,1) = 1.0
!                 A(counter,2) = helper_centre_position(1) - node_center(1)
!                 A(counter,3) = helper_centre_position(2) - node_center(2)
!                 if (dimensiona.eq.3) then
!                     A(counter,4) = helper_centre_position(3) - node_center(3)
!                 end if

!                 rho_vector(counter) = u_c(cell_index)%val(stage, rho_index)
!             end do

!             do i = 1, local_nodes(node_index)%num_cpus
!                 cpu_index = local_nodes(node_index)%rcv_offsets(i)%cpu
!                 index = (local_nodes(node_index)%rcv_offsets(i)%lower - 1)*num_values_to_send_per_node
!                 do j = local_nodes(node_index)%rcv_offsets(i)%lower, local_nodes(node_index)%rcv_offsets(i)%upper
!                     counter = counter+1
!                     A(counter,1) = 1.0
!                     do k = 1, dimensiona
!                         if ((index+k).gt.node_rcv_count(cpu_index) * num_values_to_send_per_node) then
!                             print *, "copying too much data from receive buffer from", N, "to", cpu 
!                         end if
!                         centre_positions(counter, k) = node_rcv_buffer(cpu_index)%data(index+k)
!                         A(counter,k+1) = centre_positions(counter, k) - node_center(k)
!                     end do
!                     k = dimensiona+1
!                     rho_vector(counter) = node_rcv_buffer(cpu_index)%data(index+k)
!                     index = index + num_values_to_send_per_node
!                 end do
!             end do

!             local_nodes(node_index)%density_gradient(:) = zero 

!             solvable = .false.
!             if (num_neighbours.ge.3) then
!                 if (counter.ne.num_neighbours) then
!                     print*,"something went wrong when filling the A matrix"
!                 end if

!                 At_rho_vector = zero
!                 do i = 1, (dimensiona+1)
!                     do j = 1, num_neighbours
!                         At(i,j) = A(j,i)
!                         ! At_rho_vector(i) = At_rho_vector(i) + At(i,j)*rho_vector(j)
!                         At_rho_vector(i) = At_rho_vector(i) + A(j, i)*rho_vector(j)
!                     end do
!                 end do

!                 AtA = zero
!                 do i = 1, (dimensiona+1)
!                     do j = 1, (dimensiona+1)
!                         do  k = 1, num_neighbours
!                             ! AtA(i,j) = At(i,j) + (At(i,k)*A(k,j))
!                             AtA(i,j) = AtA(i,j) + (A(k, i)*A(k, j))
!                         end do
!                     end do
!                 end do

!                 if (AtA(1,1).ne.real(num_neighbours)) then
!                     print *, "linear solver failuer\nAta(1,1)=/=num_neighbours"
!                     call abort()
!                 end if

!                 solvable = linear_solve(1+dimensiona, AtA, At_rho_vector, x)

!                 do i = 1, dimensiona
!                     local_nodes(node_index)%density_gradient(i) = x(i+1)
!                 end do
!             end if
                
!         end do
!     !$omp end do
    
!     !$omp barrier

!     !$omp master
!         call MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!     !$omp end master

! END SUBROUTINE find_node_density_gradient





! SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_2D(stage, position_index, N)
!     IMPLICIT NONE
!     integer::stage, position_index, N
!     integer::cell_index, node_index
!     integer::i, iter
!     integer::counter
!     real,dimension(1:dimensiona)::gradient
!     real::gradient_magnitude, my_max_gradient_magnitude

!     my_max_gradient_magnitude = zero

!     call find_node_density_gradient(stage, position_index, N)

!     !$omp parallel do reduction(max: max_gradient_magnitude)
!     DO cell_index = 1, XMPIELRANK(N)
!         counter = 0
!         gradient = zero
!         do i = 1, ielem(n, cell_index)%nonodes
!             node_index = ielem(n, cell_index)%nodes(i)
!             if (local_nodes(node_index)%num_neighbours.ge.3) then 
!                 counter = counter+1
!                 gradient(1:dimensiona) = gradient(1:dimensiona) + local_nodes(node_index)%density_gradient(1:dimensiona)
!             end if
!         end do

!         if (counter.eq.0) then
!             print *, "no valid neighbour for gradient in cell", cell_index
!         end if

!         gradient(1:dimensiona) = gradient(1:dimensiona)/counter

!         gradient_magnitude = zero
!         do iter = 1,dimensiona
!             gradient_magnitude = gradient_magnitude + (gradient(iter)*gradient(iter))
!         end do
!         gradient_magnitude = sqrt(gradient_magnitude)

!         u_c(cell_index)%normalized_gradient = gradient_magnitude

!         if (my_max_gradient_magnitude.lt.gradient_magnitude) then
!             my_max_gradient_magnitude = gradient_magnitude
!         end if
!     END DO
!     !$omp end parallel do

!     !$omp barrier

!     !$omp master
!        ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!         CALL MPI_ALLREDUCE(my_max_gradient_magnitude, max_gradient_magnitude, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, IERROR)
!         ! CALL MPI_BARRIER(MPI_COMM_WORLD, IERROR)
!     !$omp end master

!     !$omp barrier
!     my_max_gradient_magnitude = max_gradient_magnitude

!     ! print *, "reduction done", N, max_gradient_magnitude

!     !$omp do
!     DO cell_index = 1, XMPIELRANK(N)
!         u_c(cell_index)%normalized_gradient = u_c(cell_index)%normalized_gradient / my_max_gradient_magnitude
!     END DO
!     !$omp end do

! END SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_2D





! SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_from_precomputed(N)
!     IMPLICIT NONE
!     integer::N
!     integer::cell_index, node_index
!     integer::i, j, iter
!     integer::counter
!     real,dimension(1:dimensiona)::gradient
!     real::gradient_magnitude, factor

!     !$omp do
!     DO cell_index = 1, XMPIELRANK(N)
!         counter = 0
!         gradient = zero
!         do i = 1, ielem(n, cell_index)%nonodes
!             node_index = ielem(n, cell_index)%nodes(i)
!             ! if (local_nodes(node_index)%num_neighbours.ge.3) then 
!             if (local_nodes(node_index)%num_neighbours.ge.2) then 
!                 ! if ((local_nodes(node_index)%normalized_density_gradient_magnitude.gt.1.0).or.(local_nodes(node_index)%normalized_density_gradient_magnitude.lt.0.0)) then
!                 !     print*, "invalid normalized_density gradient in cell", cell_index, n
!                 ! end if 
!                 ! if (local_nodes(node_index)%normalized_density_gradient_magnitude.ne.zero) then
!                     counter = counter+1
!                     gradient(1:dimensiona) = gradient(1:dimensiona) + local_nodes(node_index)%density_gradient(1:dimensiona)
!                 ! end if
!                 ! gradient_magnitude = zero
!                 ! do j = 1, dimensiona
!                 !     gradient_magnitude = gradient_magnitude + (local_nodes(node_index)%density_gradient(j)*local_nodes(node_index)%density_gradient(j))
!                 ! end do
!                 ! gradient_magnitude = sqrt(gradient_magnitude)
!                 ! factor = local_nodes(node_index)%normalized_density_gradient_magnitude / gradient_magnitude
!                 ! gradient(1:dimensiona) = gradient(1:dimensiona) + (local_nodes(node_index)%density_gradient(1:dimensiona)*factor)   
!             end if
!         end do

!         if (counter.eq.0) then
!             print *, "no valid neighbour for gradient in cell", cell_index
!         end if

!         gradient(1:dimensiona) = gradient(1:dimensiona)/counter

!         gradient_magnitude = zero
!         do iter = 1,dimensiona
!             gradient_magnitude = gradient_magnitude + (gradient(iter)*gradient(iter))
!         end do
!         gradient_magnitude = sqrt(gradient_magnitude)

!         u_c(cell_index)%normalized_gradient = gradient_magnitude
!     END DO
!     !$omp end do

! END SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_from_precomputed





SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_from_precomputed(N)
    IMPLICIT NONE
    integer::N
    integer::cell_index, node_index
    integer::i, j, iter
    integer::counter
    real::gradient_magnitude

    !$omp do
    DO cell_index = 1, XMPIELRANK(N)
        counter = 0
        gradient_magnitude = zero
        do i = 1, ielem(n, cell_index)%nonodes
            node_index = ielem(n, cell_index)%nodes(i)
            ! if (local_nodes(node_index)%num_neighbours.ge.3) then 
            if (local_nodes(node_index)%num_neighbours.ge.2) then 
                counter = counter+1
                if (moving_mesh_mode.eq.13) then
                    gradient_magnitude = gradient_magnitude + (0.5*(local_nodes(node_index)%normalized_density_gradient_magnitude + local_nodes(node_index)%normalized_vf_gradient_magnitude))
                else
                    gradient_magnitude = gradient_magnitude + local_nodes(node_index)%normalized_density_gradient_magnitude
                end if
            end if
        end do

        if (counter.eq.0) then
            print *, "no valid neighbour for gradient in cell", cell_index
        end if

        gradient_magnitude = gradient_magnitude/real(counter)

        u_c(cell_index)%normalized_gradient = gradient_magnitude
    END DO
    !$omp end do

END SUBROUTINE FIND_NORMALIZED_DENSITY_GRADIENT_from_precomputed





end module