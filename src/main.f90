PROGRAM UCNS3D
!> @author
!> Panagiotis Tsoutsanis & Antonis Foivos Antoniadis
!> copyright: Panagiotis Tsoutsanis & Antonis Foivos Antoniadis
!> version: 3.9
!DESCRIPTION
!> @brief:
!> Main Driver of UCNS3D code
use mpiinfo
use translate
use declaration
use memory
use communications
use io
use partition
use library
use transform
use fluxes
use initialisation
use boundary
use advance
use recon
use local
use profile
use flow_operations
use gradients
use basis
use prestore
use riemann
use source
use implicit_time
use implicit_fluxes
use moodr
use omp_lib
use parameters



implicit none
integer::intxgt


external metis_partmeshdual
external parmetis_v3_partmeshkway
!call mpi_init(ierror)


call mpi_init_thread(mpi_thread_funneled,provided,ierror)
call mpi_comm_size(mpi_comm_world,isize,ierror)
call mpi_comm_rank(mpi_comm_world,n,ierror)


call open_input1(n,itt) !> open the input files

call tolerances  !> setup the tolerances values
call read_ucns3d !> read all the parameter files

 call close_input1(n,itt) !> close the input files


if (n.eq.0)then
  call translate_mesh !> translate the mesh from fluent msh format to native format

 end if
 call mpi_barrier(mpi_comm_world, ierror)


call timing(n,cpux1,cpux2,cpux3,cpux4,cpux5,cpux6,timex1,timex2,timex3,timex4,timex5,timex6) !> start the timers
 cpux1(1)=mpi_wtime()




 call open_arbitrary(n,imaxe,imaxn,imaxb) !> open the grid files


call shallocation(ieshape,imaxe) !> allocate arrays for shape of each element

  call mpi_barrier(mpi_comm_world, ierror)

call find_shape(n,imaxe,ieshape)	!> find the shape each element



call checkres  !> check the existence of restart/checkpoint files
call restore_outlet_mach_restart_pressure !> keep Mach outlet pressure continuous on restart

call xmpiallocate(xmpie,xmpil,xmpin,xmpielrank,xmpinrank,imaxe,imaxn,nproc) !> allocate memory for local and global numbering of elements and nodes




call mpi_barrier(mpi_comm_world, ierror)

    if (emetis.lt.6)then    !> choose a grid partitioning property if emetis<6 use serial metis
    if (n.eq.0) then

	if (emetis.eq.1)then
          call partitioner1(n,imaxe,imaxn,xmpie,ieshape)
        end if
        if (emetis .eq.2)then

	  call partitioner2(n,imaxe,imaxn,xmpie,ieshape)

	end if
	if (emetis .eq. 3) then
	  call partitioner3(n,imaxe,imaxn,xmpie,ieshape)
	end if
	if (emetis .eq. 4) then
	  call partitioner4(n,imaxe,imaxn,xmpie,ieshape)
	end if
	if (emetis.eq.5)then
	  call partitioner5(n,imaxe,imaxn,xmpie,ieshape)
	end if


    end if

      call mpi_bcast(xmpie,imaxe,mpi_integer,0,mpi_comm_world,ierror)

   else

    if (emetis.eq.6)then !> when emetis=6 use parmetis (default)
	  call partitioner6(n,imaxe,imaxn,xmpie,ieshape)
	end if
    end if



call mpi_barrier(mpi_comm_world, ierror)


call xmpifind(xmpie,xmpin,xmpielrank,xmpinrank,imaxe,imaxn,nproc) !> determine the number of elements in this process

call elallocation(n,xmpie,xmpielrank,imaxe,ieshape,itestcase,imaxb,xmin,xmax,ymin,ymax,zmin,zmax) !> allocate the appropriate memory for each elements




call read_input(n,xmpielrank,xmpinrank,xmpie,xmpin,dinode,imaxn,imaxe,imaxb,xmpinnumber,scaler,dinoder) !> read the grid files and populate the allocated memory values for vertex coordinates and numbering


call fix_offsets_local(n)

 call determine_size(n,iorder,iselem,iselemt,ioverst,ioverto,ilx,numneighbours,idegfree,imaxdegfree,iextend) !> determing the stencil sizes, number of polynomial coefficients etc.



 call gaussianpoints(igqrules,numberofpoints,numberofpoints2)   !> establish the number of gausian quadrature points for each element
 call quadalloc(numberofpoints,numberofpoints2) !>allocate the memory required for the gaussian integration rules


 call shdeallocation(ieshape,imaxe)                     !> deallocate memory for the shape allocation









call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0) then
   cpux3(1) = mpi_wtime()
   write(120+n,*)"timei_1",cpux3(1)-cpux1(1)     !> write in file the total wall clock time taken so far
end if

call allocate2
call neighbourss(n,imaxe,imaxn,xmpie,xmpin,xmpielrank,restart,dinoder) !> find neighbours of each element
 call allocate3
call mpi_barrier(mpi_comm_world,ierror)







 if (n.eq.0) then
   cpux3(1) = mpi_wtime()
   write(120+n,*)"timei_2",cpux3(1)-cpux1(1)   !> write in file the total wall clock time taken so far
end if

!$omp barrier
!!$omp parallel default(shared)

  call geometry_calc


!!$omp end parallel
!$omp barrier


 call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0)  then
   cpux3(1) = mpi_wtime()
   write(120+n,*)"timei_3",cpux3(1)-cpux1(1)
end if

!---start here
call read_bound(n,imaxb,xmpielrank)

   if (dg.eq.1)then
   if (filtering.eq.1)then
    allocate(modal_filter(1:idegfree),modal_filter_strong(1:idegfree),modal_filter_weak(1:idegfree))
    call filter_init(n)
  end if
	end if


	if (adda.eq.1)then
	allocate(adda_filter_strong(1:idegfree))
	allocate(adda_filter_weak(1:idegfree))
	call adda_filter_init(n)
	end if


!$omp barrier
!!$omp parallel default(shared)
   call apply_boundary(n,xper,yper,zper,iperiodicity,xmpielrank)
!$omp barrier
!!$omp end parallel




				call mpi_barrier(mpi_comm_world,ierror)
				!$omp barrier
				!$omp master
				if (n.eq.0)then

				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished applying boundary conditions"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time1=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier
				call mpi_barrier(mpi_comm_world,ierror)
  call xmpilocal

  !if (tecplot.lt.5)then
  call count_walls
  !end if

  call mpi_barrier(mpi_comm_world,ierror)


  !$omp barrier
!$omp master
if (lowmem.eq.0)call globalistx(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh)
if (lowmem.eq.1)call globalist(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh,glneighper)
 !$omp end master
!$omp barrier
call mpi_barrier(mpi_comm_world,ierror)

				!$omp barrier
				!$omp master
				if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished obtaining neighbours within my cpu"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time2=",cpux3(1)-cpux1(1)
				  close(63)
				  end if



 if (dimensiona.eq.3)then
 call cons(n,diconr,dicons,iperiodicity,xmpielrank,isize,diconrpa,diconrpm,diconspo,xper,yper,zper,diconrpf,numneighbours,typesten)
 else
call cons2d(n,diconr,dicons,iperiodicity,xmpielrank,isize,diconrpa,diconrpm,diconspo,xper,yper,zper,diconrpf,numneighbours,typesten)
 end if


				if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished obtaining neighbours within my cpu"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time22=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier

  if (lowmem.eq.0)call globalistx2(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh)

  if (lowmem.eq.1)call globalist2(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh,glneighper)
				  !$omp barrier
				  !$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished obtaining neighbours across all cpu"
				  cpux3(1) = mpi_wtime()
! ! ! 				  write(63,*)"time3=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier


if (ischeme.gt.1)then
 call mpi_barrier(mpi_comm_world,ierror)
  cpux3(1)=mpi_wtime()
 call allocate5
  !!$omp parallel default(shared)
 if (lowmem.eq.0)call detstenx(n)
 if (lowmem.eq.1)call detsten(n)
 !!$omp end parallel



				  !$omp barrier
				  !$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished obtaining the central stencils across all cpu"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time4=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier


  call localstallocation(n,xmpielrank,ilocalstencil,ilocalstencilper,typesten,numneighbours)
! call mpi_barrier(mpi_comm_world,ierror)



call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0)  then
   cpux3(1) = mpi_wtime()
   write(120+n,*)"time_i4",cpux3(1)-cpux1(1)
end if


!!$omp parallel default(shared)

  if ((ees.eq.0).or.(ees.ge.4))then
  if (lowmem.eq.0)call stenciilsx(n)
  if (lowmem.eq.1)call stenciils(n)
  end if
  if ((ees.gt.0).and.(ees.le.2))then
  if (lowmem.eq.0)call stenciils_eesx(n)
  if (lowmem.eq.1)call stenciils_ees(n)

  end if
!!$omp end parallel

   if (ees.eq.3)then
      call stencils3(n)

   end if



  deallocate(ilocalallelg)
  deallocate(ilocalallelgper)
  call globaldea




  call stencils(n,imaxe,xmpie,xmpielrank,ilocalstencil,typesten,numneighbours,restart)
    if (iadapt.eq.1)then
  call adapt_criterion
  end if
end if
call mpi_barrier(mpi_comm_world,ierror)

				  !$omp barrier
				  !$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finished obtaining the directional stencils across all cpu"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time5=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier





! call mpi_barrier(mpi_comm_world,ierror)
  call estabexhange(n,imaxe,xmpie,xmpin,xmpielrank,ilocalstencil,diexchanger,diexchanges,direcexr,direcexs,&
numneighbours,ischeme,isize,iperiodicity,typesten,xmpil)
!

 call renumber_neighbours(n,xmpie,xmpielrank,diexchanger,diexchanges)
!

!


call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0)  then
   cpux3(1) = mpi_wtime()
  write(120+n,*)"timei_5",cpux3(1)-cpux1(1)
end if



call deallocatempi1(n)

if (nprobes.gt.0)then
call probepos(n,probei)
end if

 !memory allocation for transformation to computational domain!

				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time6=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier






call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0) then
   cpux3(1) = mpi_wtime()
write(120+n,*)"timei_6",cpux3(1)-cpux1(1)
end if





call local_reconallocation3(n)

				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time66=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier


 call mpi_barrier(mpi_comm_world,ierror)




call exch_cords(n)
				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time67=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier


call mpi_barrier(mpi_comm_world,ierror)

if ((rungekutta.ge.10).and.(rungekutta.lt.12))then
call exch_cords2(n,isize,diexboundhiri,diexboundhisi,itestcase,numberofpoints2,diexchanger,diexchanges)
end if

					!$omp barrier
					!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time68=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier


if ((fastest.ne.1).and.(ischeme.ge.2))then
call allocate_basis_function(n,xmpielrank,idegfree)
end if

				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time69=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier

! !$omp parallel default(shared)
!  call memory1
! !$omp end parallel


call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0)  then
   cpux3(1) = mpi_wtime()
write(120+n,*)"timei_7",cpux3(1)-cpux1(1)
end if

					!$omp barrier
					!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time70=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier



  if (stencil_io.eq.1)then
   call stenprint(n)
  end if



call solex_alloc(n)


				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time71=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier



if (rungekutta.ge.2)then
!this is parallel


!!$omp parallel default(shared)
if (dimensiona.eq.3)then
call direct_side(n)
else
call direct_side2d(n)
end if
!!$omp end parallel
end if

					!$omp barrier
					!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time72=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier



call mpi_barrier(mpi_comm_world,ierror)



call exch_cord3(n)
if (iperiodicity.eq.1)then
call read_input_period(n,xmpielrank,xmpinrank,xmpie,xmpin,dinode,imaxn,imaxe,imaxb,xmpinnumber,scaler)
end if
 deallocate(dinoder2)

call mpi_barrier(mpi_comm_world,ierror)





call mpi_barrier(mpi_comm_world,ierror)
   cpux2(1) = mpi_wtime()

  if ((dg.eq.1).or.(adda_type.eq.2))then

 call allocate_dg

 end if

  call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0)  then
   cpux3(1) = mpi_wtime()
 write(120+n,*)"timei_8",cpux3(1)-cpux1(1)
end if

					!$omp barrier
					!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time7=",cpux3(1)-cpux1(1)

				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier

			if ((fastest.ne.1).and.(ischeme.ge.2)) call prestore_1(n)



					!$omp barrier
					!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"finishded prestoring reconstruction matrices"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time8=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier




call dealcordinates2	!for ale do not deallocate
call localsdeallocation(n,xmpielrank,ilocalstencil,ilocalstencilper,typesten,numneighbours)	!for ale do not deallocate
call deallocatempi2(n)	!for ale do not deallocate
call dealcordinates1(n,diexcordr,diexcords) !for ale do not deallocate


call mpi_barrier(mpi_comm_world,ierror)








  if (turbulence.eq.1)then
    if (dimensiona.eq.3)then
    call walldistancex(n,imaxe,xmpielrank)
    else
    call walldistance2d(n,imaxe,xmpielrank)
    end if
  end if




call mpi_barrier(mpi_comm_world,ierror)


				  !$omp barrier
				  !$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"started prestoring geometry information"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time9=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier

call mpi_barrier(mpi_comm_world,ierror)
 if (n.eq.0)then
   cpux3(1) = mpi_wtime()
  write(120+n,*)"timei_9",cpux3(1)-cpux1(1)
end if

!for ale this needs to be repeated
!!$omp parallel default(shared)
call grads_assign(n)
call find_angles(n)
!!$omp end parallel



!--------------thursday-------------!




				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"allocating solution  and flux variables"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time10=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier

call u_c_allocation(n,xmpielrank,itestcase)


if (dg.eq.1)then

!!$omp parallel default(shared)
call build_mass_matrix(n)
!!$omp end parallel

end if

!$omp barrier
!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"initialising"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time10=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier
call mpi_barrier(mpi_comm_world,ierror)

    !!$omp parallel default(shared)
    call initialise (n)
   !!$omp end parallel





if (restart.gt.0)then
   call rest_read(n)
end if




 call globaldea2(xmpil,xmpie)  !for ale do not deallocate


				call mpi_barrier(mpi_comm_world,ierror)
				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"flux allocation"
				  cpux3(1) = mpi_wtime()
				  write(63,*)"time10=",cpux3(1)-cpux1(1)
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier

 call sumflux_allocation(n)






 if (rungekutta.ge.10)then

  call impallocate(n)


 end if

!
!
!
! ! !----------------------------------------------------------------!
! ! !		advancement of solution in time			 !
! ! !           different options available depending on scheme	 !
! ! !----------------------------------------------------------------!

 call mpi_barrier(mpi_comm_world,ierror)
   cpux6(1)=mpi_wtime()

  call timers(n,cpux1,cpux2,cpux3,cpux4,cpux5,cpux6,timex1,timex2,timex3,timex4,timex5,timex6)


! !$omp parallel default(shared)
!   call memory2
! !$omp end parallel





  if (statistics.eq.1)then
  if (n.eq.0)then

    write(st_n_cpu,fmt='(i10)') isize
    thread_n=omp_get_max_threads()
    write(st_n_threads,fmt='(i10)') thread_n
    statfile="stats_mpi_"//trim(adjustl(st_n_cpu))//"_threads_"//trim(adjustl(st_n_threads))//".txt"
    open(133,file=statfile,form='formatted',status='replace',action='write')
    write(133,'(5x,a2,1x,a11,a11,a11,a11,a11,a11,a11,a11,a11,a11)')"#it","t_time","t_comm","t_comp","t_dgint","t_halo","t_recon","t_bound","t_adda","t_flux","t_update"
    close(133)

  end if
  end if


!------------friday--------------!


if (fastest_q.eq.1)then
    call memory_fast(n)




end if





call new_arrays(n)







call exch_cords_opt(n)



if (fastest.ne.1)then
    call exchange_higher_pre(n)
end if



call fix_nodes_local	!this has to be investigated for ale








call specify_write_variables(n)


call mpi_barrier(mpi_comm_world,ierror)

select case(tecplot)

case(5)



call partition_preparation(n)

if (outsurf.eq.1)then
call prepare_surfaces_v(n)
call partition_preparation_wallv(n)
end if

case(6)

call partition_preparation_p(n)
if (outsurf.eq.1)then
call partition_preparation_p_wall(n)
end if




end select


call mpi_barrier(mpi_comm_world,ierror)
				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"before flattening"
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier
				  call mpi_barrier(mpi_comm_world,ierror)




!this creates a flat array for receiving data from mpi neighbours!
call flat_comm(n)
call alloc_weights(n)


				call mpi_barrier(mpi_comm_world,ierror)
				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,*)"before mapping"
				  close(63)
				  end if
				  !$omp end master
				  !$omp barrier
				  call mpi_barrier(mpi_comm_world,ierror)

call omp_map_first(n)


		if (n.eq.0)print*,"ucns3d running"


!
			call mpi_barrier(mpi_comm_world,ierror)
				!$omp barrier
				!$omp master
				  if (n.eq.0)then
				  cpux3(1) = mpi_wtime()
				  end if
				  !$omp end master
				  !$omp barrier
				  call mpi_barrier(mpi_comm_world,ierror)






#ifdef xpu
!$omp target data map(to: n, aa_1, adda, adda_1, adda_1_s, adda_2, adda_2_s, adda_alpha_1, adda_alpha_2, &
!$omp& adda_type, allnodesgloball, allres, allresdt, alls, alpha, alpha_0, alpha_inf1, alpha_inf2, alpha_star0, alpha_starinf, &
!$omp& angle_per, aoa, average_restart, averaging, beta, beta_i1, beta_i2, beta_starinf, beta_t, betaas, binio, bleed, &
!$omp& bleed_number, bleed_type, boundtype, br2_damping, br2_yn, bubble_centre, bubble_radius, c_des_sa, c_des_sst, &
!$omp& c_mu_inlet, c_sas, c_smg, cascade, catalytic_wall, cavitation, cb1, cb2, cfl, cflmax, cflramp, cfw, charlength, &
!$omp& chunk_n, code_profile, ct1, ct2, ct3, ct4, cv1, cw1, cw2, cw3, d_corr, datatypeint, datatypex, datatypexx, &
!$omp& datatypey, datatypeyy, datatypez, des_model, dg, dimensiona, dims, dt, ees, ek_time, emetis, eta2_sas, every_time, &
!$omp& extended_bounds, extf, fastest, fastest_q, fastmovie, fil_alpha, fil_nc, fil_s, filter_type, filtering, firstorder, &
!$omp& firstrese, firstresk, firstresomega, firstrespass, firstresr, firstrest, firstresu, firstresv, firstresw, forcex, &
!$omp& forcey, forcez, gamma, governingequations, greengo, gridar1, gridar2, guassianquadra, hybridist, i_turb_inlet, iadapt, &
!$omp& ibcode, iboundary, ibside, icarlos1, icarlos2, icompact, icong, iconimp, iconsgvq, iconsr, icoupleturb, idegfree, &
!$omp& idegfree2, idegfree3, ievery, ievery2, ieveryav, iforce, igianagraps, igqrules, ihax1, ihybrid, iloop, iloopx, ilx, &
!$omp& imaxb, imaxdegfree, imaxdegfree2, imaxe, imaxn, in, ind1, indicator_type, init_mu_ratio, initcond, initialres, inum2, &
!$omp& inwhichel, iorder, iorder2, ioverst, ioverto, iperiodicity, ires_turb, ires_unsteady, iriemann, irs, ischeme, &
!$omp& iscoun, iselem, ispal, isplit, issf, istn, it, itestcase, itold, itotalb, itt, ivortex, iweightlsqr, iweno, iwmaxe, &
!$omp& jk, jtotal, jtotal1, jtotal2, jtotal3, jump_cond1, jump_cond2, jump_cond3, kappa, kappa_sst, kdum1, kdum2, kdum3, &
!$omp& kill, kinit_srf, kloopx, kmaxn, l0norm, l1norm, l2norm, l_turb_inlet, lam, lamps, ccfl,lamx, lamy, lamz, limiter, lmach, &
!$omp& lmach_style, lowmem, lowmemory, lwci1, m_t0, Mach_in, modeio, momentx, momenty, momentz, max_faces, max_fnodes, max_nodes, mood, mood_mode, mood_var1, &
!$omp& mood_var2, mood_var3, mood_var4, movement, mp_modelc, mrf, multispecies, n_boundaries, nderivative, nodes_i, &
!$omp& nodes_part, nof_bounded, nof_bubbles, nof_interior, nof_species, nof_variables, nprobes, nproc, nrotors, ntmax, &
!$omp& num_dg_dofs, num_dg_reconstruct_dofs, num_hexas, num_prisms, num_pyramids, num_tetras, numberofpoints, numberofpoints2, numneighbours, numneighbours2, oo2, out_time, &
!$omp& output_freq, outsurf, part1_end, part2_end, part3_end, part4_end, part5_end, passivescalar, per_rot, pi, poly, pr_t1, &
!$omp& pr_t2, pr_t3, pr_t4, pr_t5, pr_t6, pr_t7, pr_t8, prace_t1, prace_t2, prace_t3, prace_t4, prace_t5, prace_t6, &
!$omp& prace_t7, prace_t8, prace_t9, prace_tx1, prace_tx2, prace_tx3, prandtl, pres, press_outlet, mach_outlet_target, prev_turbmodel, prevres, &
!$omp& prtu, qp_hexa, qp_line, qp_line_n, qp_prism, qp_pyra, qp_quad, qp_quad_n, qp_tetra, qp_triangle, qp_triangle_n, &
!$omp& qrde, qsas_model, r_beta, r_gas, r_k_sst, r_om_sst, realgas, reduce_comp, relax, required, res_time, rescounter, &
!$omp& rescountert, residualfreq, reslimit, resmax, resmaxt, restart, reynolds, rframe, rg_kf_type, rg_nof_reactions, &
!$omp& rg_nof_tv_coef, rg_relax, rg_t_inf, rg_t_ref, rg_t_wall_init, rg_ttr, rg_tve, rhc1, rhc2, rhc3, rhc4, rot_corr, rres, &
!$omp& rungekutta, scaler, schmidt_lam, total_pressure_inlet, total_temperature_inlet, density_inlet,schmidt_turb, sigma, sigma_k1, sigma_k2, sigma_om1, sigma_om2, sigma_phi,a405,nof_perturbations405, &
!$omp& source_active, transition_model, transition_axis, transition_direction, transition_ramp_type, transition_location, transition_ramp_length, spatialorder, spatiladiscret, spkin, spos, srf_origin, srf_velocity, origin, srfg, &
!$omp& statistics, stencil_io, stennorm, subdiv, surfshear, suther, swirl, t, taylor, taylor_ens, taylor_ensx, &
!$omp& tecplot, temp_model, temporder, thermal, thread_n, timestep, tol_per, tolbig, tolsmall, totalvolume, totiw, totwalls, &
!$omp& totwallsc, turbinit, turbulence, turbulenceequations, turbulencemodel, twall, typ_countn, typ_countn_global, &
!$omp& typ_countn_global_w, typ_countn_w, typesten, tz1, ufreestream, unwou, upperlimit, upturblimit, uvel, v_ref, &
!$omp& vectorx, vectory, vectorz, visc, viscous_s, &
!$omp& voll, vorder, vort_model, vvel, wall_temp, wallc, wdatatypeint, wdatatypex, wdatatypexx, wdatatypey, wdatatypeyy, &
!$omp& wdatatypez, weight_lsqr, wenocentralweight, wenocnschar, wenoz, wenwrt, wkdum1, wkdum2, wkdum3, wnodes_part, &
!$omp& wpart1_end, wpart2_end, wpart3_end, wpart4_end, wpart5_end, write_variables, write_variables_av, write_variables_av_w, &
!$omp& write_variables_w, wvel, xper, pos_l,pos_g,ipos_l,ipos_g, yper, zero, zero_turb_init, zeta_star, zper, indicator_par1, indicator_par2, indicator_par3, jtot, adda_filter_strong, adda_filter_weak, bleed_end, bleed_plenum, bleed_porosity, bleed_start, bound_len, bounds_len, bound_offset, bounds_offset, boundhir_dg, boundhis_dg, catalytic_con, el_bnd, el_int, gamma_in, &
!$omp& halo_len, halos_len, halo_offset, halos_offset, halo_proc, halos_proc, ibound_cpun, ibound_face, ibound_ibid, ibound_ibl, ibound_icode, ibound_inum, ibound_ishape, ibound_localn, ibound_nibl, ibound_nlocal, ibound_t, ibound_t2, ibound_which, totk,totens,totensx,kill_nan, &
!$omp& ielem_admis, ielem_avars, ielem_bleedn, ielem_condition, ielem_condx, ielem_dih, ielem_dih2, ielem_diss, ielem_dtl, ielem_er, ielem_er1, ielem_er1dt, ielem_er1er2, ielem_er2, ielem_er2dt, ielem_erx, ielem_faceanglex, ielem_faceangley, ielem_facediss, &
!$omp& ielem_filtered, ielem_full, ielem_ggs, ielem_hybrid, ielem_ibounds, ielem_idegfree, ielem_ifca, ielem_ihex, ielem_ihexgl, ielem_indexf, ielem_indexi, ielem_ineigh, ielem_ineighb, ielem_ineighg, ielem_ineighn, ielem_inter_id, ielem_interior, ielem_inumneighbours, ielem_iorder, ielem_ishape, ielem_itotalpoints, ielem_linc, ielem_lwcx2, ielem_minedge, ielem_mode, ielem_mood, ielem_mood_o, &
!$omp& ielem_nodes, ielem_nodes_faces, ielem_nodes_faces_v, ielem_nodes_neighbours, ielem_nodes_v, ielem_nofbc, ielem_nojecount, ielem_nonodes, ielem_q_face_q_mapl, ielem_qface, ielem_recalc, ielem_reduce, ielem_reorient, ielem_stencil_dist, ielem_surf, ielem_totvolume, ielem_troubled, ielem_types_faces, ielem_vdec, ielem_viscx, ielem_vortex, ielem_walldist, ielem_walltrans, ielem_walls, ielem_wcx, ielem_xxc, ielem_yyc, ielem_zzc, &
!$omp& impdiag_mf, impoff_mf, impdiag, impdiagt, impdu, impoff, impofft, inoder4_bct, inoder4_cord, inoder4_itor, integ_basis_dg_value, integ_basis_value, integ_basis_valuec, m_1_val, modal_filter, modal_filter_strong, modal_filter_weak, &
!$omp& mp_a_in, mp_janaf, mp_m, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in, mrf_rot_gl, point1_gl, point2_gl, qp_array_qp_weight, qp_array_x, qp_array_y, qp_array_z, radius_gl, &
!$omp& rec_br2_aux_var, rec_br2_local_lift, rec_cgradientstemp, rec_cond, rec_findw, rec_g0, rec_gradf, rec_gradients, rec_gradients2, rec_gradientsc, rec_gradientsc2, rec_gradientstemp, rec_gradientstemp_wall, rec_gradientsturb, rec_gradientsturb_wall, rec_grads, rec_gradsav, &
!$omp& rec_ihexb, rec_ihexbc, rec_ihexg, rec_ihexgc, rec_ihexl, rec_ihexlc, rec_ihexn, rec_ihexnc, rec_indicator, rec_indicatorc, rec_invccjac, rec_invctjac, rec_invmat_stencilt, rec_invmat_stenciltc, rec_k0, rec_local, rec_mrf, rec_mrf_origin, rec_mrf_velocity, rec_periodicflag, rec_qpoints, rec_qpoints_p, rec_rotvel, rec_rpoints, rec_stencils, rec_stencilsc, rec_surf_qpoints, rec_tempsq, rec_tempsqmat, rec_uleft, rec_uleft_dg, rec_uleftturb, rec_uleftturbv, rec_uleftv, rec_uleftx, rec_velinvlsqmat, rec_vellsq, rec_velocitydof_wall, rec_vext_ref, rec_volume, rec_volume_w, rec_volumec, rec_wall, rec_wallcoeff, rec_wallcoefg, rec_weightl, rec_weno, rec_weno2, rec_wenos, &
!$omp& rg_hzero, rg_molm, rg_thetag, rg_tv_coef, rg_vf, rgs_ab, rgs_bb, rgs_cb, rgs_eps_over_k, rgs_mg, rgs_sigmaa, rhs_sol_mm_dg, rhs_val, rhs_valdg, rhst_val, sht, sht_rg, nodelist, xmpielrank, dg2fv, &
!$omp& u_c_br2_aux_var, u_c_rms, u_c_val, u_c_valdg, u_cs_val, u_cs_valdg, u_ct_val, u_cw_val, u_cw_valdg, u_e_val, &
!$omp& solhir, solhis, solhird, solhisd, solhi_loc, solhir_flat, solhis_flat, boundhiri, boundhisi, boundhir, boundhis, boundhirm, boundhism, boundhir_flat, boundhis_flat, boundhir_dgflat, boundhis_dgflat, boundhiri_flat, boundhisi_flat, &
!$omp& ineedhalo, ineedhalos, ineedbound, ineedbounds, bound_total, bounds_total, halo_total,ilength1,ilength2, halos_total, need_side, need_q, need_loc, bound_proc, bounds_proc,weights_t,weights_l,weights_q)

#endif













if (dimensiona.eq.3)then

#if defined(gpu) || defined(xpu)
    call time_marching(n)
#else
    !$omp parallel default(shared)
    call time_marching(n)
    !$omp end parallel
#endif

else

#if defined(gpu) || defined(xpu)
    call time_marching2(n)
#else
    !$omp parallel default(shared)
    call time_marching2(n)
    !$omp end parallel
#endif

end if

#ifdef xpu
!$omp end target data
#endif

call mpi_barrier(mpi_comm_world,ierror)


cpux3(1) = mpi_wtime()


if (n.eq.0)  then
write(120+n,*)"total time taken=",cpux3(1)-cpux2(1),"seconds"
end if



call mpi_barrier(mpi_comm_world,ierror)


call mpi_finalize(ierror)

if (n.eq.0) print*,"ucns3d finished running"




!time 3d to time 2d fix


END PROGRAM UCNS3D
