module declaration
!--------------------------------------------------------------------------------------------------------------------------!
!hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh!
!_____________________________________________ start header of module______________________________________________________!
!function: is to collect all the global variables, and data types of the code in one module
!**************************************************************************************************************************!
!_____________________________________________ end header of module________________________________________________________!
!--------------------------------------------------------------------------------------------------------------------------!

implicit none


! Compile-time upper bounds for fixed-size local work arrays.
! GPU builds should pass case-specific GPU_MAX_* values from gpu_max_flags.py.
! CPU builds do not need UCNS3D.DAT at compile time, so their fallback values
! are deliberately larger and can be overridden from the Makefile if needed.
#ifndef GPU_MAX_IORDER
#ifdef gpu
#define GPU_MAX_IORDER 4
#else
#define GPU_MAX_IORDER 6
#endif
#endif
#ifndef GPU_MAX_DIM
#define GPU_MAX_DIM 3
#endif
#ifndef GPU_MAX_SPECIES
#ifdef gpu
#define GPU_MAX_SPECIES 1
#else
#define GPU_MAX_SPECIES 5
#endif
#endif
#ifndef GPU_MAX_TURBULENCE
#ifdef gpu
#define GPU_MAX_TURBULENCE 1
#else
#define GPU_MAX_TURBULENCE 2
#endif
#endif
#ifndef GPU_MAX_PASSIVE
#ifdef gpu
#define GPU_MAX_PASSIVE 1
#else
#define GPU_MAX_PASSIVE 8
#endif
#endif
#ifndef GPU_MAX_EXTF
#define GPU_MAX_EXTF 3
#endif
#ifndef GPU_MAX_QP_FACE
#ifdef gpu
#define GPU_MAX_QP_FACE 25
#else
#define GPU_MAX_QP_FACE 36
#endif
#endif
#ifndef GPU_MAX_QP_VOLUME
#ifdef gpu
#define GPU_MAX_QP_VOLUME 210
#else
#define GPU_MAX_QP_VOLUME 336
#endif
#endif
#ifndef GPU_MAX_QP_ALL
#ifdef gpu
#define GPU_MAX_QP_ALL 125
#else
#define GPU_MAX_QP_ALL 216
#endif
#endif
#ifndef GPU_MAX_TYPESTEN
#ifdef gpu
#define GPU_MAX_TYPESTEN 7
#else
#define GPU_MAX_TYPESTEN 15
#endif
#endif
#ifndef GPU_MAX_NVAR
#ifdef gpu
#define GPU_MAX_NVAR 5
#else
#define GPU_MAX_NVAR 11
#endif
#endif
#ifndef GPU_MAX_NODES
#define GPU_MAX_NODES 8
#endif
#ifndef GPU_MAX_FACES
#define GPU_MAX_FACES 6
#endif
#ifndef GPU_MAX_FNODES
#define GPU_MAX_FNODES 4
#endif
#ifndef GPU_MAX_IDEGFREE
#if GPU_MAX_DIM == 2
#define GPU_MAX_IDEGFREE ((((GPU_MAX_IORDER + 1) * (GPU_MAX_IORDER + 2)) / 2) - 1)
#else
#define GPU_MAX_IDEGFREE ((((GPU_MAX_IORDER + 1) * (GPU_MAX_IORDER + 2) * (GPU_MAX_IORDER + 3)) / 6) - 1)
#endif
#endif
#ifndef GPU_MAX_DOF
#define GPU_MAX_DOF (GPU_MAX_IDEGFREE + 1)
#endif
#ifndef GPU_MAX_NEIGHBOURS
#if GPU_MAX_DIM == 2
#if GPU_MAX_IORDER == 1
#define GPU_MAX_NEIGHBOURS 5
#else
#define GPU_MAX_NEIGHBOURS (GPU_MAX_DOF * GPU_MAX_EXTF)
#endif
#else
#if GPU_MAX_IORDER == 1
#define GPU_MAX_NEIGHBOURS 9
#else
#define GPU_MAX_NEIGHBOURS (GPU_MAX_DOF * GPU_MAX_EXTF)
#endif
#endif
#endif
#ifndef GPU_MAX_NEIGHBOURS2
#if GPU_MAX_DIM == 2
#if GPU_MAX_IORDER == 1
#define GPU_MAX_NEIGHBOURS2 6
#elif GPU_MAX_IORDER <= 3
#define GPU_MAX_NEIGHBOURS2 7
#else
#define GPU_MAX_NEIGHBOURS2 11
#endif
#else
#if GPU_MAX_IORDER == 1
#define GPU_MAX_NEIGHBOURS2 7
#elif GPU_MAX_IORDER <= 3
#define GPU_MAX_NEIGHBOURS2 12
#else
#define GPU_MAX_NEIGHBOURS2 20
#endif
#endif
#endif

integer, parameter :: gpu_max_iorder = GPU_MAX_IORDER
integer, parameter :: gpu_max_dim = GPU_MAX_DIM
integer, parameter :: gpu_max_species = GPU_MAX_SPECIES
integer, parameter :: gpu_max_turbulence = GPU_MAX_TURBULENCE
integer, parameter :: gpu_max_passive = GPU_MAX_PASSIVE
integer, parameter :: gpu_max_extra_transport = gpu_max_turbulence + gpu_max_passive
integer, parameter :: gpu_max_qp_face = GPU_MAX_QP_FACE
integer, parameter :: gpu_max_qp_volume = GPU_MAX_QP_VOLUME
integer, parameter :: gpu_max_qp_all = GPU_MAX_QP_ALL
integer, parameter :: gpu_max_typesten = GPU_MAX_TYPESTEN
integer, parameter :: gpu_max_nodes = GPU_MAX_NODES
integer, parameter :: gpu_max_faces = GPU_MAX_FACES
integer, parameter :: gpu_max_fnodes = GPU_MAX_FNODES
integer, parameter :: gpu_max_idegfree = GPU_MAX_IDEGFREE
integer, parameter :: gpu_max_dof = GPU_MAX_DOF
integer, parameter :: gpu_max_neighbours = GPU_MAX_NEIGHBOURS
integer, parameter :: gpu_max_neighbours2 = GPU_MAX_NEIGHBOURS2
integer, parameter :: gpu_max_nvar = GPU_MAX_NVAR
integer, parameter :: gpu_max_nvar_total = gpu_max_nvar + gpu_max_extra_transport

!--------------------------------------------------------------------------------------------------------------------------!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! s.1.   integer  variables here                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!--------------------------------------------------------------------------------------------------------------------------!

integer::cascade,mood,mood_mode,kmaxn,multispecies,nof_species,dimensiona,lowmem,binio,nof_variables,chunk_n,dims,ires_turb,ind1,code_profile,ires_unsteady,lamps,itotalb,totiw,icompact,ees,iscoun,itold,lmach_style,weight_lsqr	!dimensions of problem
integer::governingequations,filtering,guassianquadra,temporder,iboundary,wenocnschar,required,nodes_i,swirl,iadapt,tecplot,stencil_io,surfshear,issf,n_boundaries,fastest_q,statistics,adda,alls
integer,allocatable,dimension(:,:)::jtot1,jtot2,jtot3,jtot,el_connect
integer::jtotal,jtotal1,jtotal2,jtotal3,fastmovie,movement,typ_countn_global,typ_countn,typ_countn_global_w,typ_countn_w
integer::filter_type,fil_nc,fil_s,fil_alpha		!filter values
integer::adda_type,adda_alpha_1,adda_alpha_2,adda_1_s,adda_2_s,adda_1,adda_2		!filter values
integer::kill			 ! flag for killing a simulation
integer::nof_perturbations405	!number of perturbations
integer::nderivative		  	! index of the numbering of the component of the polynomials
real::output_freq				!output frequency in simulation time
integer::extended_bounds		!bounds strict or relaxed
integer::iconsr					! index for identifying the considered element
integer::inwhichel				! index for identifying the considered element within a stencil
integer::subdiv					! index for the number of subdivisions of a face
integer::iconsgvq				! index for considered cell
integer::igianagraps				!index for considered cell
integer::it					! number of iterations
integer::iconimp				!index for considered cell
integer::ibcode					! index for boundary condition code
integer::ibside					! index for which side of cell is bounded
integer::cfw					!index for determining from the which section the boundary subroutine is called
integer::dg, br2_yn                     ! flag for dg discretisation
real:: br2_damping
real::ccfl
real::totk,totens,totensx		!totk,totens,totensx,kill_nan
real::res_sum,pos_l1,pos_l2,pos_l3,pos_l4	!res_sum,pos_l1,pos_l2,pos_l3,pos_l4,ipos_l1,ipos_l2
integer::ipos_l1,ipos_l2
integer::kill_nan
real,dimension(5)::pos_l,pos_g
integer,dimension(3)::ipos_l,ipos_g
real::r_gas						!specific gas constant
real,allocatable,dimension(:):: weights_q,weights_t,weights_l
integer::lowmemory				! memory usage flag
integer::fastest				! fastest mode flag for the code
integer::emetis					! type of metis partitioning
integer::allnodesgloball			! total number of nodes in the domain
integer::spatiladiscret				! spatial discretisation
integer::spatialorder    			! spatial order
integer::jk					!index
integer::numneighbours,numneighbours2				! number of neighbours in the stencils
integer::restart				! flag for determining if a restart file is present
integer::ihybrid				!hybrid mode where some part of the domain is solved with lower order and another with higher order schemes
integer::upperlimit !
integer::residualfreq				!frequency for writing the residuals
integer::irs					!implicit residual smoothing version
integer::greengo				!gradient approximation techniques
integer::lmach					!low mach preconditioning
integer::ispal					!type of spalart allmaras modification
integer::relax					!relaxation schemes gauss seidel, jacobian etc
integer::icong					!index for considered cell
integer::iloop					!index for considered cell
integer::cflramp				!cfl ramping activated for implicit time stepping
integer::averaging				!flag for activating averaging of transient data
integer::average_restart			!flag for having a restart file with the averaged solution
integer::prev_turbmodel				!flag of previous turbulence model in case of a restart file
integer::iforce					!index on how often to compute the forces
integer::boundtype				!type of boundary condition supersonic,  subsonic
integer::nproc					!number of cpus
integer::initcond				!initial condition type
integer::wenwrt					!weno wrt to which variable conserved, characteristics etc
integer::poly					!interpolating polynomials 1 generic, 2 legendre
integer::ivortex				!q criterion com putation
integer::icarlos1				!index for considered element passive scalars, complicated inflow condition
integer::icarlos2				!index for considered face passive scalars, complicated inflow condition
integer::imaxe					!total number of elements of the grid
integer::ntmax					!maximum number of iterations
integer::ihax1					!temporary integer used as pointer for operations regarding hexahedral elements
integer::limiter				!slope limiter
integer::iextend				!extend parameter for central big stencil only for weno cases
integer::idegfree,idegfree2,inum2,idegfree3		!degrees of freedom of polynomial it will become loca	l
integer::outsurf				!outsurf stands for computing forces
integer::unwou					!type of output files to write all partitions to one file etc
integer::igqrules				!gaussian quadrature rule
integer::imaxdegfree,imaxdegfree2				!maximum number of neighbours
integer::qp_hexa,qp_tetra,qp_pyra,qp_prism,qp_quad,qp_triangle,qp_line,qp_line_n,reduce_comp,qp_quad_n,qp_triangle_n !volume gaussian quadrature points
!> number of points allocated to qpoints, wequa (coords and weights for quadrature points)
integer::numberofpoints,numberofpoints2
integer::rescounter				!counter for residuals
integer::rescountert				!counter for turbulence residuals
integer::ilx					!degrees of freedom for required order of accuracy
integer::imaxb					!total number of boundary conditions
integer::imaxn					!total number of nodes present in the grid
integer::in					!index for considered cell
integer::iorder,iorder2					!order of accuracy (spatial)
integer::ioverst				! number of elements for the stencil overdetermined matrix required for avoiding ill conditions
integer::ioverto				! number of elements for the estab stencils overdetermined matrix required for avoiding ill conditions
integer::qrde					!type of solution for the least squares problem
integer::rungekutta     			!time stepping scheme type
integer::ischeme 				!type of scheme to be used
integer::isplit					!unsplit finit volume scheme
integer::iselem					!maximum number of elements required for the stencil dependent on the order of accuracy
integer::istn					!variable that takes each time the value of the considered element  to find neighbours in the stencil
integer::itt,wenoz					!integer variable used for opening files
integer::itestcase              !types of equations to be solved
integer::iweightlsqr    			!weighted least squres option of inverted distance
integer::typesten				!number of stencils for weno
integer::iperiodicity				!periodicity detected in domain
integer::iweno					!type of schemes weno, muscl, unlimited
integer::turbulence				!turbulence equations
integer::turbulencemodel			!which model
integer::turbulenceequations			!how many equations to solve
integer::icoupleturb				!coupled turbulence model
integer::iriemann				!riemann solver
integer::firstorder				!first order scheme active
integer::passivescalar				!additional transport equations
integer::qsas_model				!qsas model
integer::vort_model				!vorticity model
integer::zero_turb_init				!flag for type of initialization for k-omega
integer::des_model 				!integer switches for sst
integer:: nprobes,totwalls,nof_interior,nof_bounded,mrf			!number of probes for transient data
integer:: rot_corr,d_corr   !integer for turbulence corrections
integer:: transition_model,transition_axis,transition_direction,transition_ramp_type
integer::max_faces
integer::max_fnodes
integer::max_nodes
integer::num_hexas
integer::num_tetras
integer::num_pyramids
integer::num_prisms
!--------------------- variables for parallel partitioned output-------!
integer,allocatable,dimension(:)::dispart1,dispart2,dispart3,dispart4,dispart5,typ_nodesn,typ_nodesn_w
integer,allocatable,dimension(:)::iarray_part1,iarray_part2,iarray_part3,iarray_part4,iarray_part5,i_array_part2x
real,allocatable,dimension(:)::rarray_part4,rarray_part2,rarray_part3,rarray_part5
integer,allocatable,dimension(:)::wdispart1,wdispart2,wdispart3,wdispart4,wdispart5,wallcount_cpu_l,wallcount_cpu_g
integer,allocatable,dimension(:)::wiarray_part1,wiarray_part2,wiarray_part3,wiarray_part4,wiarray_part5
integer,allocatable,dimension(:)::offset_vtu,connect_vtu,type_vtu,offset_vtu_w,connect_vtu_w,type_vtu_w
real,allocatable,dimension(:,:)::sol_vtu,sol_vtu_w
real,allocatable,dimension(:)::nodes_vtu,nodes_vtu_w
real,allocatable,dimension(:)::wrarray_part4,wrarray_part2,wrarray_part3,wrarray_part5
real,allocatable,dimension(:,:)::rarray_part1,wrarray_part1
integer::part1_end,part2_end,part3_end,part4_end,part5_end
integer::wpart1_end,wpart2_end,wpart3_end,wpart4_end,wpart5_end
integer::kdum1,kdum2,write_variables,write_variables_av,nodes_part
integer::wkdum1,wkdum2,write_variables_w,write_variables_av_w,wnodes_part
integer::datatypex,datatypey,datatypez,datatypexx,datatypeyy,datatypeint
integer,dimension(1)::kdum3
character(len=25)::variable_names(20),variable_names_av(20)
integer::wdatatypex,wdatatypey,wdatatypez,wdatatypexx,wdatatypeyy,wdatatypeint
integer,dimension(1)::wkdum3
character(len=25)::variable_names_w(20),variable_names_av_w(20)
integer::kloopx,iloopx,totwallsc,iwmaxe
integer,allocatable,dimension(:)::wallit,offsetwall,wall_nodes,wallcx,offsetwc,offsetwc_g,wallcx_g,wallshape,wallshape_g,wallshape_g2
integer,allocatable,dimension(:,:)::wall_l
!--------------------- end of variables for parallel partitioned output-------!
real,dimension(3)::srf_origin,srf_velocity,origin
real::press_outlet,mach_outlet_target,mach_outlet_average,mach_outlet_relax
real::mach_outlet_filtered,mach_outlet_filter_alpha
integer::mach_outlet_update_freq,mach_outlet_start_iter,mach_outlet_filter_ready
integer::rframe,source_active
real::per_rot,angle_per,v_ref,kinit_srf,srfg,tol_per
real,allocatable,dimension(:,:)::point1_gl,point2_gl
real,allocatable,dimension(:)::radius_gl,mrf_rot_gl
integer:: nrotors
integer::num_dg_dofs !number of degrees of freedom for dg polynomial approximation
integer::num_dg_reconstruct_dofs !number of degrees of freedom for dg order + 1, not including num_dg_dofs
integer::indicator_type                !troubled indicator type
integer::viscous_s,jump_cond1,jump_cond2,jump_cond3
integer::cavitation
real::indicator_par1,indicator_par2,indicator_par3, bound_lim  !troubled indicator parameters
real::rhc1,rhc2,rhc3,rhc4
real::a405   !perturbations amplitude
real::total_pressure_inlet,total_temperature_inlet,density_inlet
real::transition_location,transition_ramp_length
real::prace_t1,prace_t2,prace_t3,prace_t4,prace_t5,prace_t6,prace_t7,prace_t8,prace_t9,pr_t1,pr_t2,pr_t3,pr_t4,pr_t5,pr_t6,pr_t7,pr_t8,prace_tx1,prace_tx2,prace_tx3
!------------------start bleed parameters-------------------!
integer::bleed_number,bleed,bleed_type
real,allocatable,dimension(:,:)::bleed_start,bleed_end
real,allocatable,dimension(:)::bleed_plenum,bleed_porosity
!------------------end bleed parameters-------------------!
!------------------------real gas effects section variables-----------------!
!species order 'n2', 'o2', 'no', 'n ', 'o '
integer::realgas			!flag for real gas effects
integer::rg_nof_reactions	!number of reactions
integer::rg_kf_type		    !type of reaction rate equations: 1. gupta, 2. candler
integer::rg_relax		    !type of expression for vibrational relaxation time: 1. low t mw, high t park; 2. sum of mw and park
real,allocatable,dimension(:)::rg_vf	!mass fraction of each species ([n2, o2, no, n, o])
real,allocatable,dimension(:)::rg_molm  !molar mass of each species ([n2, o2, no, n, o])
real,allocatable,dimension(:)::rg_hzero !enthalpy of formation: n2, o2, no, n, o
real,allocatable,dimension(:)::rg_thetag	!diatomic vibrational characteristic temperatures [k]; atoms 0
real,allocatable,dimension(:)::catalytic_con	!species concentrations at the wall
integer::catalytic_wall						!flag for catalytic wall
real::rg_t_inf								!temperature infinity (k)
real::rg_t_wall_init										!temperature_wall_initial (k)
real::rg_t_ref								!temperature_reference (k)
integer::rg_nof_tv_coef						!number of tv-t equation coefficients
real,allocatable,dimension(:,:)::rg_tv_coef	!tv-t equation coefficients ([n2,o2,no],[p1, p2, p3, p4]) (number of tv-t equation coefficients is not 0)
real,parameter :: rgs_ru = 8.31446261815324	! j/kmol-k	!universal gas constast
real::rg_ttr								!translational-rotational temperature
real::rg_tve								!vibrational temperature
real,allocatable,dimension(:)::rgs_mg 		!rg_molm * 1.0e3   ! g/mol for d_ij correlation
real, parameter :: rgs_pa_per_atm = 101325.0	!pascals per atmosphere	101325.0
real, parameter :: rgs_cm2s_to_m2s = 1.0e-4	!cm2s to m2s 1.0e-4
real, parameter :: rgs_tiny = 1.0e-30 		!tiny number for real gas 1.0e-30
! blottner viscosity coefficients (base-10 form)
real,allocatable,dimension(:):: rgs_ab, rgs_bb,rgs_cb
! classic lennard–jones parameters (σ [å], ε/k [k])
real,allocatable,dimension(:):: rgs_sigmaa     	!(/ 3.667, 3.467, 3.492, 3.298, 3.050 /)
real,allocatable,dimension(:):: rgs_eps_over_k !(/ 99.8 ,106.7 ,116.7 , 71.4, 80.0  /)
!-------------------------------------------------------------------------!
!----------------multiphysics-----------------------!
!constants read from file
integer::mp_modelc	!multispecies mode; 0=allaire (no diffusion of species), 1=multispecies
real,allocatable,dimension(:)::gamma_in	!gamma for each species
real,allocatable,dimension(:)::mp_a_in !volume fraction for each species in the inlet
real,allocatable,dimension(:)::mp_r_in !density for each species in the inflow
real,allocatable,dimension(:)::mp_pinf !p infinity for each species for the stiffened gas eos
real,allocatable,dimension(:)::mp_m	!molecular weight kg/mol
real,parameter,dimension(5)::mp_brok_a = (/ &
  0.06170928, 0.10338607, 0.10039271, 0.02670999, 0.04674248 /) ! A-Blottner coefficients for each species in base-10 form
real,parameter,dimension(5)::mp_brok_b = (/ &
  0.31800000, -0.08260000, -0.03360000, 0.60300000, 0.42900000 /) ! B-Blottner coefficients for each species in base-10 form
real,parameter,dimension(5)::mp_brok_c = (/ &
  1.09247235, 2.00449077, 1.83945886, 0.61474842, 0.96218401 /) ! C-Blottner coefficients for each species in base-10 form
real,allocatable,dimension(:)::mp_tlow_in, mp_thigh_in,mp_tmid_in	!mp_tlo_in
real,allocatable,dimension(:,:,:) :: mp_janaf !janaf coefficients for cp/r: cp = r*(a1 + a2*t + a3*t^2 + a4*t^3 + a5*t^4)
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! s.2.   integer allocatable variables here        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!--------------------------------------------------------------------------------------------------------------------------!

integer,allocatable,dimension(:,:)::probei	!probe positions
integer,allocatable,dimension(:)::nodes_offset_local,nodes_offset_local2,nodes_offset_localw,nodes_offset_local2w	!offsets for the nodes for each cell in my cpu from the global list
integer,allocatable,dimension(:)::nodes_offset,nodes_offset2,nodes_offsetw,nodes_offsetw2	!offsets for the nodes for each cell  globally
integer,allocatable,dimension(:,:)::glneigh	!globalist of all the neighbours in the grid
integer,allocatable,dimension(:,:)::glneighper	!globalist of all the rotational periodic neighbours in the grid
integer,allocatable,dimension(:)::xsize,xgo,xgo_v,chunk_size	!allocatable index for describing which stencil is considered within the directional stencil construction routine
integer,allocatable,dimension(:)::my_nodesl,my_nodesg	!allocatable arrays for local and global indexing of nodes belonging to each cpu (no duplications)
integer,allocatable,dimension(:)::ibound_t,ibound_t2	!allocatable index for describing the shyape of the considered cell
integer,allocatable,dimension(:)::iselemt	!total number of elements in each stencil (extended one)
integer,allocatable,dimension(:)::xmpiall,offset,xmpiwall,woffset,xmpi_re,xmpi_wre,offset_v,xmpiall_v		!allocatable index for stencil routine
integer,allocatable,dimension(:)::el_int,el_bnd
integer,allocatable,dimension(:)::pare,dose,pareel,pares
integer,allocatable,dimension(:,:)::doseel,soseel
integer,allocatable,dimension(:,:,:,:)::ilocalstencil	!4-d array for stencils
integer,allocatable,dimension(:,:,:,:)::ilocalstencilper	!4-d array for rotational periodic stencils
integer,allocatable,dimension(:,:,:,:)::ilocalallelg	!4-d array for stencils
integer,allocatable,dimension(:,:,:,:)::ilocalallelgper	    !4-d array for rotational periodic stencils
integer,allocatable,dimension(:)::xmpie			!global list of number of elements that belong to each cpu
integer,allocatable,dimension(:)::xmpil			!local list of number of elements that belong to each cpu
integer,allocatable,dimension(:)::xmpin			!not used can get rid of
integer,allocatable,dimension(:)::xmpielrank		!number of elements in this cpu
integer,allocatable,dimension(:)::xmpinrank		!not used can get rid of
integer,allocatable,dimension(:,:,:)::xmpinnumber	!not used can get rid of
integer,allocatable,dimension(:)::stcon			!dummy variable for recursive subroutine of stencils
integer,allocatable,dimension(:)::stconc		!dummy variable for recursive subroutine of stencils
integer,allocatable,dimension(:)::stcons		!dummy variable for recursive subroutine of stencils
integer,allocatable,dimension(:)::list,ineb,iperb,nodelist		!dummy variable for recursive subroutine of stencils
integer,allocatable,dimension(:,:,:,:)::ilocalalls      !dummy variable for recursive subroutine of stencils
integer::thermal,temp_model
!--------------------------------------------------------------------------------------------------------------------------!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! s.3.   real variables here        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!--------------------------------------------------------------------------------------------------------------------------!
real::wenocentralweight,timestep,oo2,zero,forcex,forcey,forcez,extf,vorder,mood_var1,mood_var2,mood_var3,mood_var4
real::pi					!pi trigonometri
real::taylor,taylor_ens,taylor_ensx					!only to be used for taylor green vortex
real::voll					!total volume of the domain
real::wall_temp					!wall temperature model
real::upturblimit				!upper turbulence viscosity ratio
real::hybridist					!upper turbulence viscosity ratio
real,parameter ::to4=3.0d0/4.0d0
real,parameter ::oo4=1.0d0/4.0d0
real,parameter ::to3=2.0d0/3.0d0
real,parameter ::oo3=1.0d0/3.0d0
real,parameter :: turb_diag_floor_frac=5.0d-2
real,parameter :: turb_source_cap_frac=5.0d-1
real,parameter :: turb_diag_abs_floor=1.0d-30
real::reynolds					!reynolds number
real::cflmax					!maximum number of allowable cfl to be used only with implicit
real::prevres					!previous residual in order to determine ramping strategy
real::vectorx					!setting up with respect to which plane the aoa is defined x
real::vectory					!setting up with respect to which plane the aoa is defined y
real::vectorz					!setting up with respect to which plane the aoa is defined z
real::gridar1					!aspect ratio of grid
real::gridar2					!maximum volume aspect ratio of stencils
real::t						!time
real::totalvolume				!total volume of domain
real::reslimit       				!limit to stop the simulation

real::res_time					!restart_time
real::ievery2					!how often to write output
real::wallc					!wall clock limit in seconds
real::betaas					!betaas of sutherland exponent
real::suther					!sutherland non dimensional constant of temperature
real::prandtl					!prandtl number
real::scaler					!for scaling the mesh
real::lwci1					!linear weno weight of central stencils				!
real::resmax					!residual dimensional mean flow
real::resmaxt					!residual dimensional turbulence
real::cfl					!cfl number
real::ievery,modeio					!how often to write output
real::ieveryav					!how often to write output of average solution
real:: firstresu				!residual dimensional mean flow
real::firstresv					!residual dimensional mean flow
real::firstresw					!residual dimensional mean flow
real::firstrese					!residual dimensional mean flow
real::firstresr					!residual dimensional mean flow
real::firstrest					!residual dimensional mean flow
real::firstresk					!residual dimensional mean flow
real::firstresomega				!residual dimensional mean flow
real::firstrespass				!residual dimensional mean flow
real::gamma
real::uvel					!u velocity
real::ufreestream				!u velocity
real::Mach_in					!prescribed inlet Mach number
real::vvel					!v-velocity
real::wvel					!w-velocity
real::pres					!pressure
real::rres					!density
real::spos					!speed of sound
real::spkin					!specific kinetic energy
real::lam					!heat conductivity (watt/(mk))
real::visc					!viscosity
real::lamx					!flags to be used for restarting
real::lamy					!flags to be used for restarting
real::lamz					!flags to be used for restarting
real::alpha,beta		!flags to be used for restarting
real::out_time					!final time to write output for unsteady simulations
real::every_time,ek_time					!final time to write output for unsteady simulations
real::xper					!periodicity in x axis
real::yper					!periodicity in y axis
real::zper					!periodicity in z axis
real::aoa					!angle of attack
real::charlength				!characteristic length if undefined, and reynolds undefined the free stream will be used for air
						!turbulence model constants
real::turbinit				!turbulence initial value vt/v
real::cb1
real::cb2
real::sigma
real::kappa
real::cw1
real::cw2
real::cw3
real::cv1
real::ct1
real::ct2
real::ct3
real::ct4
real::prtu
real::twall
character(len=30)::statfile,st_n_cpu,st_n_threads
integer::thread_n
real::sigma_k1
real::sigma_k2
real::sigma_om1
real::sigma_om2
real::aa_1
real::beta_i1
real::beta_i2
real::alpha_starinf
real::alpha_0
real::beta_starinf
real::r_beta
real::r_k_sst
real::beta_t
real::kappa_sst
real::r_om_sst
real::zeta_star
real::m_t0
real::c_smg
real::eta2_sas
real::sigma_phi
real::c_sas
real::alpha_inf1
real::alpha_inf2
real::alpha_star0
real::c_mu_inlet   						!functions for sst
real::l_turb_inlet
real::i_turb_inlet
real::init_mu_ratio
real::c_des_sa
real::c_des_sst
real:: schmidt_lam
real::schmidt_turb
real::tolsmall,tolbig,allresdt,tz1
real,dimension(1:gpu_max_nvar_total)::allres,initialres
real,dimension(100,3)::bubble_centre
real,dimension(100)::bubble_radius
integer::nof_bubbles

real:: momentx,momenty,momentz
!--------------------------------------------------------------------------------------------------------------------------!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! s.4.   real allocatable variables here        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!--------------------------------------------------------------------------------------------------------------------------!
real,allocatable,dimension(:)::modal_filter,adda_filter_weak,adda_filter_strong,modal_filter_strong,modal_filter_weak
real::l1norm		!l1 norm of solution for grid convergence studies of euler and linear advection equations
real::l2norm		!l2 norm of solution for grid convergence studies of euler and linear advection equations
real::l0norm,stennorm		!l0 norm of solution for grid convergence studies of euler and linear advection equations
real::dt 					!real time step size
real,allocatable,dimension(:)::avrg		!temporary solution averages
real,allocatable,dimension(:,:)::jac		!jacobians
real,allocatable,dimension(:,:)::inversejac	!inverse jacobians
real,allocatable,dimension(:)::deterjac		!determinant jacobians
real,allocatable,dimension(:,:)::probec		!probe positions
real,allocatable,dimension(:,:)::impdu		!implicit only change of solution
real,allocatable,dimension(:,:)::impdu_old	!implicit jacobi snapshot of change of solution
real,allocatable,dimension(:,:,:)::impdiag	!implicit only diagonal matrix d
real,allocatable,dimension(:,:,:,:)::impoff	!implicit only off diagonal matrix d
real,allocatable,dimension(:)::impdiag_mf	!implicit only diagonal matrix d
real,allocatable,dimension(:,:)::impoff_mf	!implicit only off diagonal matrix d
real,allocatable,dimension(:,:,:)::impofft	!implicit only off diagonal matrix d for turbulence
real,allocatable,dimension(:,:)::impdiagt	!implicit only  diagonal matrix d for turbulence
real,allocatable,dimension(:,:)::sht		!implicit only  source term jacobian
real,allocatable,dimension(:,:)::sht_rg		!implicit only  soruce term jacobian
real,allocatable,dimension(:)::maxdiff		!se
real,allocatable,dimension(:)::mindiff		!se
real,allocatable,dimension(:)::tmaxdiff		!se
real,allocatable,dimension(:)::tmindiff		!se
real,allocatable,dimension(:)::xmin		!se
real,allocatable,dimension(:)::xmax		!se
real,allocatable,dimension(:)::ymin		!se
real,allocatable,dimension(:)::ymax		!se
real,allocatable,dimension(:)::zmin		!se
real,allocatable,dimension(:)::zmax		!se
real,allocatable,dimension(:)::cpux1		!timer
real,allocatable,dimension(:)::cpux2		!timer
real,allocatable,dimension(:)::cpux3		!timer
real,allocatable,dimension(:)::cpux4		!timer
real,allocatable,dimension(:)::cpux5		!timer
real,allocatable,dimension(:)::cpux6		!timer
real,allocatable,dimension(:)::cpux7		!timer
real,allocatable,dimension(:)::timex1		!timer
real,allocatable,dimension(:)::timex2		!timer
real,allocatable,dimension(:)::timex3		!timer
real,allocatable,dimension(:)::timex4		!timer
real,allocatable,dimension(:)::timex5		!timer
real,allocatable,dimension(:)::timex6		!timer
real,allocatable,dimension(:,:)::ifin		!temporary pointer for sync output
real,allocatable,dimension(:,:)::tfin		!temporary pointer for sync output
real,allocatable,dimension(:,:)::centerr	!for directional stencils
integer, allocatable,dimension(:)::cand,cands,candr
integer, allocatable,dimension(:,:)::candxr,candxs,cand2s,cand2rt
integer, allocatable,dimension(:,:,:)::cand2r
real, allocatable,dimension(:,:)::xand2s,xand2rt
real, allocatable,dimension(:,:,:)::xand2r
real,allocatable,dimension(:)::flux_term_left_z,flux_term_left_x,flux_term_left_y
real,allocatable,dimension(:)::flux_term_right_z,flux_term_right_x,flux_term_right_y
real,allocatable,dimension(:,:)::sind1,sind2,sind3,sind4,sind5,sind6
integer::ineedhalo,ineedhalos				       !global flat arrays for the number of processors that i need halo cells from
integer,allocatable::ineedbound,ineedbounds		   !global flat arrays for the number of processors that  for boundaries
integer,allocatable::halo_len(:),halos_len(:)	   !how many elements from each processor for the halo cells
integer,allocatable::halo_offset(:),halos_offset(:)	   !their offset
integer,allocatable::bound_len(:),bounds_len(:)	   !how many elements from each processor for the boundary cells
integer,allocatable::need_side(:),need_q(:),need_loc(:)
integer, allocatable::halo_proc(:),halos_proc(:)
integer,allocatable::solhi_loc(:)
integer,allocatable::bound_proc(:),bounds_proc(:)
integer,allocatable::bound_offset(:),bounds_offset(:)   !their offset
real,allocatable::solhir(:,:),solhis(:,:)		   !receiving flat array for halo cells mean values
real,allocatable::solhird(:),solhisd(:)			!receiving flat array for halo cells adda terms
real,allocatable::boundhiri(:,:),boundhisi(:,:)		!receiving flat array for boundary cells implicit time stepping
real,allocatable::boundhir(:,:),boundhis(:,:)			!receiving flat array for boundary cells gradients and solutions
real,allocatable::boundhir_dg(:,:),boundhis_dg(:,:)		!receiving flat array for boundary cells gradients and solutions for dg
integer,allocatable::boundhirm(:),boundhism(:)		!receiving flat array for mood states
real,allocatable::solhir_flat(:),solhis_flat(:),boundhir_flat(:),boundhis_flat(:)
real,allocatable::boundhir_dgflat(:),boundhis_dgflat(:), boundhiri_flat(:), boundhisi_flat(:)
integer::bound_total,bounds_total, halo_total, halos_total,ilength1,ilength2

!--------------------------------------------------------------------------------------------------------------------------!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! s.5.   data type variables here        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo!
!--------------------------------------------------------------------------------------------------------------------------!
type::aneixx
      integer::cpu                              !cpu index
      integer,allocatable,dimension(:,:)::elem1 !elements indexing
      real,allocatable,dimension(:)::centers    !centers of cells for weno directionals stencils and compact stencils
end type aneixx

type(aneixx),allocatable,dimension(:)::dneix1

type aexchange_cord
	integer::procid        !id of processor
	integer::howmany       ! how many elements are needed
	real,allocatable,dimension(:,:,:)::nodecord    !the coordinates of the element
end type aexchange_cord

type(aexchange_cord),allocatable,dimension(:)::diexcordr  !receiving data type
type(aexchange_cord),allocatable,dimension(:)::diexcords  !sending data type


type aexchange_solhi
	integer::procid    !id of processor
	integer::iavt      !number of processors that we need to receive/or send data
	integer::iavc      !number of processors that we need to receive/or send data
	integer::fast
	integer::howmany   !how many elements are needed
	real,allocatable,dimension(:,:)::sol   !array to hold the values to be send/received
end type aexchange_solhi

type(aexchange_solhi),allocatable,dimension(:)::diexsolhir,diexsolhird    !receiving data type for halo cells of stencils
type(aexchange_solhi),allocatable,dimension(:)::diexsolhis,diexsolhisd    !sending data type for halo cells of stencils



type aexchange_boundhi
	integer::procid   !id of processor
	integer::iavt     !number of processors that we need to receive/or send data
	integer::fast
	integer::howmany  !how many elements are needed
	real,allocatable,dimension(:,:)::facesol,facesol_m, facesol_dg  !array holding the number of variables to be sent/received
	integer,allocatable,dimension(:,:)::vertpp !vertex mapping between different processes
end type aexchange_boundhi

type(aexchange_boundhi),allocatable,dimension(:)::diexboundhir !receiving data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries
type(aexchange_boundhi),allocatable,dimension(:)::diexboundhis !sending data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries
type(aexchange_boundhi),allocatable,dimension(:)::diexboundhiri !receiving data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping
type(aexchange_boundhi),allocatable,dimension(:)::diexboundhisi !sending data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping
type(aexchange_boundhi),allocatable,dimension(:)::diexboundhirr  !receiving data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping
type(aexchange_boundhi),allocatable,dimension(:)::diexboundhiss  !sending data type for boundary extrapolated values at gaussian quadrature points of inter-processor boundaries for implicit time stepping



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
type::anode_ne	!name of type for the set of nodes
	integer::itor  !index declaring if the coordinates of this node are needed for each cpu
	real,allocatable,dimension(:)::cord    !coordinates of the node
	integer,allocatable,dimension(:)::bct
	integer::itorm	!index refering to unique list index in cells
end type anode_ne
type::anode_lit	!name of type for the set of nodes
	integer::numberofneib      !number of elements that share this node
	integer,allocatable,dimension(:)::neibids,xne,xneib !ids of the elements, and counters
end type anode_lit

type(anode_ne),allocatable,dimension(:)::dinoder
type(anode_lit),allocatable,dimension(:)::dinoder2


integer,allocatable,dimension(:)::ieshape	!1-d array for shape of element

type::aexchange
	integer::procid					!processor id
	integer::tot                    ! total number of elements
	integer,allocatable,dimension(:)::whattheyneed   !element numbers
	integer,allocatable,dimension(:)::muchtheyneed  !number of elements they need
	integer,allocatable,dimension(:)::muchineed	!number of elements i need
	integer,allocatable,dimension(:)::whatineed	!element number global
	integer,allocatable,dimension(:)::sideineed !side that i need
	integer,allocatable,dimension(:)::sideineedn,qineed    !from which cpu, which gaussian quadrature point order
	integer,allocatable,dimension(:)::sidetheyneed,qtheyneed !side they need, which gaussian quadrature point order
	integer,allocatable,dimension(:)::sidetheyneedn
	integer,allocatable,dimension(:)::localref !local reference
	integer,allocatable,dimension(:,:)::nodex
end type aexchange


type(aexchange),allocatable,dimension(:)::diexchanger !receive
type(aexchange),allocatable,dimension(:)::diexchanges !send
type(aexchange),allocatable,dimension(:)::diexchanger1
type(aexchange),allocatable,dimension(:)::diexchanges1

type::asolexchange
	integer::procid        !processor id
	real,allocatable,dimension(:,:)::centres   !cell centres
	integer,allocatable,dimension(:,:)::nodes  !nodes
	real,allocatable,dimension(:,:)::sol,sol_dg       !solution
end type asolexchange

type(asolexchange),allocatable,dimension(:)::dsolchanger  !receives

type(asolexchange),allocatable,dimension(:)::dsolchanges  !sends

type::arecex
	integer::procid					!processor id
	integer::tot
	integer,allocatable,dimension(:)::whattheyneed   !element numbers
	integer,allocatable,dimension(:)::muchtheyneed  !number of elements they need
	integer,allocatable,dimension(:)::muchineed	!number of elements i need
	integer,allocatable,dimension(:)::whatineed	!element number global
	integer,allocatable,dimension(:)::ishape   !shape of elements
	integer,allocatable,dimension(:)::localref !local referencing
	real,allocatable,dimension(:,:)::centers   !barycentres
end type arecex

type(arecex),allocatable,dimension(:)::direcexr		!receive elements due to stencils
type(arecex),allocatable,dimension(:)::direcexr1		!receive elements due to stencils
type(arecex),allocatable,dimension(:)::direcexrg		!receive elements due to stencils
type(arecex),allocatable,dimension(:)::direcexsg		!send elements due to stencils
type(arecex),allocatable,dimension(:)::direcexs1		!send elements due to stencils
type(arecex),allocatable,dimension(:)::direcexs		!send elements due to stencil




type::aconnx	!name of type for the set of elements
	integer::procid	!cpu id
	integer,allocatable,dimension(:)::howmanyi	!how many element i need
	integer,allocatable,dimension(:)::retm     !how many entries to receive
	integer,allocatable,dimension(:)::ret      !how many entries to send
	integer,allocatable,dimension(:)::howmanythey !how many element they need
	integer,allocatable,dimension(:,:)::whichi !which element i need
	integer,allocatable,dimension(:,:)::whichthey  !which element they need
	real,allocatable,dimension(:,:)::facx      !fac
end type aconnx

type(aconnx),allocatable,dimension(:)::diconr
type(aconnx),allocatable,dimension(:)::dicons
type(aconnx),allocatable,dimension(:)::diconrpa
type(aconnx),allocatable,dimension(:)::diconrpm
type(aconnx),allocatable,dimension(:)::diconspo
type(aconnx),allocatable,dimension(:)::diconrpf



type::anode_number	!name of type for the set of nodes
	integer::noden	!identification number that can be used as a pointer inside an array
	integer::blockn !block number or cpu number of this element
	integer::nodegl
	real::x		!coordinates in x axis
	real::y		!coordinates in y axis
	real::z		!coordinates in z axis
end type anode_number


type(anode_number),allocatable,dimension(:,:)::dinode	  !1-d array for pointer type for nodes





real, allocatable :: u_c_val(:,:,:)				!mean flow variables
real, allocatable :: u_ct_val(:,:,:)			!turbulence + passive scalars
real, allocatable :: u_e_val(:,:,:)				!exact solutions
real, allocatable :: u_cs_val(:,:,:)			!mean flow variables with strong filter
real, allocatable :: u_cw_val(:,:,:)			!mean flow variables with weak filter
real, allocatable :: u_c_valdg(:,:,:,:)			!mean flow variables for each dof of DG
real, allocatable :: u_cs_valdg(:,:,:,:)		!mean flow variables for each dof of DG with strong filter
real, allocatable :: u_cw_valdg(:,:,:,:)		!mean flow variables for each dof of DG with weak filter
real, allocatable :: u_c_rms(:,:)				!time averaged rms
real, allocatable :: u_c_br2_aux_var(:,:,:,:)	!BR2 aux variables
real, allocatable :: m_1_val(:,:,:)				!dg mass matrices
real, allocatable :: rhs_val(:,:)				!rhs of mean flow variables
real, allocatable :: rhs_valdg(:,:,:)			!rhs of mean flow variables for each dof of DG
real, allocatable :: rhs_sol_mm_dg(:,:,:)		!rhs of mean flow variables for the matrix of DG
real, allocatable :: rhst_val(:,:)				!rhs of turbulence variables
real, allocatable :: integ_basis_value(:,:)
real, allocatable :: integ_basis_valuec(:,:)
real, allocatable :: integ_basis_dg_value(:,:)
real, allocatable :: dg2fv(:,:,:)
real, allocatable :: qp_array_x(:,:)
real, allocatable :: qp_array_y(:,:)
real, allocatable :: qp_array_z(:,:)
real, allocatable :: qp_array_qp_weight(:,:)
integer, allocatable :: inoder4_bct(:,:)
integer, allocatable :: inoder4_itor(:)
real, allocatable :: inoder4_cord(:,:)


!==================== ielem_ ====================

! integers
integer, allocatable :: ielem_admis(:)  ! (kmaxe)  number of admissible stencils
integer, allocatable :: ielem_bleedn(:,:)  ! (nof_faces,kmaxe)  bleed boundary conditions
integer, allocatable ::ielem_inter_id(:)
integer, allocatable ::ielem_indexf(:)				   ! (face index for interpocessor boundary)
integer, allocatable :: ielem_condx(:)                 ! (kmaxe)
integer, allocatable :: ielem_nofbc(:)                 ! (kmaxe)
integer, allocatable :: ielem_filtered(:)              ! (kmaxe)
integer, allocatable :: ielem_full(:)  ! (kmaxe)  specifies if sufficient number of stencils are found to proceed with weno for this cell
integer, allocatable :: ielem_ggs(:)  ! (kmaxe)  specifies with what algorithm to compute the gradients (green gauss, least squares or blend of them)
integer, allocatable :: ielem_hybrid(:)  ! (kmaxe)  flag for switching to lower order discretisation as a function of wall distance
integer, allocatable :: ielem_ibounds(:,:)  ! (nof_faces,kmaxe)  bounded codes for each bounded face
integer, allocatable :: ielem_idegfree(:)  ! (kmaxe)  degrees of freedom for polynomial selected
integer, allocatable :: ielem_ifca(:)  ! (kmaxe)  number of sides
integer, allocatable :: ielem_ihex(:)  ! (kmaxe) local index of each cell
integer, allocatable :: ielem_ihexgl(:)  ! (kmaxe)  global index of each cell
integer, allocatable :: ielem_indexi(:,:)  ! (nof_faces,kmaxe)  indexing for the cells
integer, allocatable :: ielem_ineigh(:,:)  ! (nof_faces,kmaxe)  neighbours local numbering
integer, allocatable :: ielem_ineighg(:,:)  ! (nof_faces,kmaxe)  neighbours global numbering
integer, allocatable :: ielem_ineighb(:,:)  ! (nof_faces,kmaxe)  neighbours cpu index
integer, allocatable :: ielem_ineighn(:,:)  ! (nof_faces,kmaxe)  neighbours numbering in other cpus
integer, allocatable :: ielem_interior(:)  ! (kmaxe)  specifies if this cell has any side bounded (interior cells get a value of 0, non interior ones get a value of 1)
integer, allocatable :: ielem_inumneighbours(:)  ! (kmaxe)  number of neighbours in each stencil
integer, allocatable :: ielem_iorder(:)  ! (kmaxe)  order of polynomials
integer, allocatable :: ielem_ishape(:)  ! (kmaxe)  > shape of element 1: hex 2: tet 3: pyramid 4: prism 5: quad 6: tri
integer, allocatable :: ielem_itotalpoints(:)          ! (kmaxe)
integer, allocatable :: ielem_mode(:)  ! (kmaxe)  specifies if decomposition canbe avoided for straight sided element (not implemented yet)
integer, allocatable :: ielem_mood(:)  ! (kmaxe)  mood flag for every element
integer, allocatable :: ielem_mood_o(:)  ! (kmaxe)  mood flag for every element
integer, allocatable :: ielem_nodes(:,:)  ! (nof_nodes,kmaxe)  nodes index
integer, allocatable :: ielem_nodes_faces(:,:,:)  ! (nof_faces,face_nodes,kmaxe)  counterclockwise numbering of the nodes for each face
integer, allocatable :: ielem_nodes_faces_v(:,:,:)  ! (nof_faces,face_nodes,kmaxe)  counterclockwise numbering of the nodes for each face
integer, allocatable :: ielem_nodes_neighbours(:,:,:)  ! not used  nodes index
integer, allocatable :: ielem_nodes_v(:,:)  ! (nof_nodes,kmaxe)  nodes_rearranged for output
integer, allocatable :: ielem_nojecount(:,:)  ! (nof_nodes,kmaxe)  number of nodes
integer, allocatable :: ielem_nonodes(:)  ! (kmaxe)  number of nodes
integer, allocatable :: ielem_recalc(:)  ! (kmaxe)  flag for recalculating the solution post-eriori
integer, allocatable :: ielem_reduce(:)                ! (kmaxe)
integer, allocatable :: ielem_reorient(:,:)  ! (nof_faces,kmaxe)  consistency across interface in terms of ordering of quadrature points
integer, allocatable :: ielem_troubled(:)              ! (kmaxe)
integer, allocatable :: ielem_types_faces(:,:)  ! (nof_faces,kmaxe)  type of each face (quadrilateral, triangle)
integer, allocatable :: ielem_vdec(:)  ! (kmaxe)  number of volume decompositions for each element
integer, allocatable :: ielem_walls(:)  ! (kmaxe)  flag to declare if this is cell bounded by a wall
integer, allocatable :: ielem_q_face_q_mapl(:,:,:)     ! (max_qp_face,nof_faces,kmaxe)

! reals
real, allocatable :: ielem_avars(:,:)  ! (nof_variables,kmaxe)  q criterion
real, allocatable :: ielem_condition(:)                ! (kmaxe)
real, allocatable :: ielem_dih(:,:)  ! (nof_faces,kmaxe)  distance across cell centres at each face
real, allocatable :: ielem_dih2(:,:,:)  ! (nof_faces,1:3,kmaxe)  distance across cell centres at each face
real, allocatable :: ielem_diss(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_dtl(:)  ! (kmaxe)  local time step size
real, allocatable :: ielem_er(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_er1(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_er1dt(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_er1er2(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_er2(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_er2dt(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_erx(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_faceanglex(:,:)  ! (nof_faces,kmaxe)  faceangle
real, allocatable :: ielem_faceangley(:,:)  ! (nof_faces,kmaxe)  faceangles y
real, allocatable :: ielem_facediss(:,:)  ! (nof_faces,kmaxe)  dissipation for mood
real, allocatable :: ielem_linc(:)  ! (kmaxe)  central stencil linear weight
real, allocatable :: ielem_lwcx2(:)  ! (kmaxe)  adda components
real, allocatable :: ielem_minedge(:)  ! (kmaxe)  inscribed sphere radius
real, allocatable :: ielem_stencil_dist(:)  ! (kmaxe)  stencil distance factor
real, allocatable :: ielem_surf(:,:)  ! (nof_faces,kmaxe)  surface area
real, allocatable :: ielem_totvolume(:)  ! (kmaxe)  volume of element
real, allocatable :: ielem_viscx(:)  ! (kmaxe)  local time step size
real, allocatable :: ielem_vortex(:,:)  ! (1:3,kmaxe)  q criterion
real, allocatable :: ielem_walldist(:)  ! (kmaxe)  wall distance
real, allocatable :: ielem_walltrans(:)  ! (kmaxe)  signed nearest-wall distance from transition location
real, allocatable :: ielem_wcx(:)                    ! (kmaxe)
real, allocatable :: ielem_xxc(:)  ! (kmaxe)  cell centre coordinates in x
real, allocatable :: ielem_yyc(:)  ! (kmaxe)  cell centre coordinates in y
real, allocatable :: ielem_zzc(:)  ! (kmaxe)  cell centre coordinates in z
integer, allocatable :: ielem_qface(:,:,:)



!==================== rec_ ====================

! integers
integer, allocatable :: rec_g0(:)    !constrained least squares gaussian elimination component
integer, allocatable :: rec_ihexb(:,:,:)    !cpu that that each cell belongs to
integer, allocatable :: rec_ihexbc(:,:,:)    !cpu that that each cell belongs to
integer, allocatable :: rec_ihexg(:,:,:)    !global index of cells
integer, allocatable :: rec_ihexgc(:,:,:)    !global index of cells
integer, allocatable :: rec_ihexl(:,:,:)    !local index of cells
integer, allocatable :: rec_ihexlc(:,:,:)    !local index of cells
integer, allocatable :: rec_ihexn(:,:,:)    !internal index from where to take the values from communicated messages
integer, allocatable :: rec_ihexnc(:,:,:)    !internal index from where to take the values from communicated messages
integer, allocatable :: rec_k0(:)    !constrained least squares gaussian elimination component
integer, allocatable :: rec_local(:)
integer, allocatable :: rec_wall(:)
integer, allocatable :: rec_mrf(:)
integer, allocatable :: rec_periodicflag(:,:,:)

! reals
real, allocatable :: rec_br2_aux_var(:,:,:,:,:)    !(var, dim, i_face, i_qp)
real, allocatable :: rec_br2_local_lift(:,:,:,:)    !(var, dim, i_face)
real, allocatable :: rec_cgradientstemp(:,:,:,:)    !unlimited gradients for temperature for wall cells
real, allocatable :: rec_cond(:,:)    !dummy variable used for gradient approximation estimation
real, allocatable :: rec_findw(:,:,:,:,:)    !weno weights for main equations with respect to characteristic variables
real, allocatable :: rec_gradf(:,:,:)    !unlimited gradients for velocity
real, allocatable :: rec_gradients(:,:,:,:)    !reconstructed gradients for main variables
real, allocatable :: rec_gradients2(:,:,:,:)    !reconstructed gradients for turbulent variables for wall cells
real, allocatable :: rec_gradientsc(:,:,:,:)    !reconstructed gradients for main variables
real, allocatable :: rec_gradientsc2(:,:,:,:)    !reconstructed gradients for turbulent variables for wall cells
real, allocatable :: rec_gradientstemp(:,:)    !unlimited gradients for temperature
real, allocatable :: rec_gradientstemp_wall(:,:,:)    !unlimited gradients for temperature for wall cells
real, allocatable :: rec_gradientsturb(:,:,:,:)    !unlimited gradients for turbulent variables for
real, allocatable :: rec_gradientsturb_wall(:,:,:,:)    !unlimited gradients for turbulent variables for wall cells
real, allocatable :: rec_grads(:,:,:)    !gradients obtained from green gauss approximation
real, allocatable :: rec_gradsav(:,:,:)    !gradients obtained from green gauss approximation
real, allocatable :: rec_indicator(:,:,:)    !precomputed smoothness indicators
real, allocatable :: rec_indicatorc(:,:,:)    !precomputed smoothness indicators
real, allocatable :: rec_invccjac(:,:,:)    !inverse jacobian
real, allocatable :: rec_invctjac(:,:,:)    !inverse jacobian transposed
real, allocatable :: rec_invmat_stencilt(:,:,:,:)    !pseudo inverse matrix for least squares reconstruction
real, allocatable :: rec_invmat_stenciltc(:,:,:,:)    !pseudo inverse matrix for least squares reconstruction
real, allocatable :: rec_mrf_origin(:,:)
real, allocatable :: rec_mrf_velocity(:,:)
real, allocatable :: rec_qpoints(:,:,:,:)    !quadrature points
real, allocatable :: rec_qpoints_p(:,:,:,:)    !quadrature points physical only for AI training
real, allocatable :: rec_rotvel(:,:,:,:)    !radius of qpoints, rotational velocity
real, allocatable :: rec_rpoints(:,:,:,:)    !radius of qpoints, rotational velocity
real, allocatable :: rec_stencils(:,:,:,:)    !stencils entries for matrix a (usually stored only for wall bounded cells)
real, allocatable :: rec_stencilsc(:,:,:,:)    !stencils entries for matrix a (usually stored only for wall bounded cells)
real, allocatable :: rec_surf_qpoints(:,:,:,:)    !physical space surface quadrature points (i_face, i_qp, xy) relative to cell center
real, allocatable :: rec_tempsq(:,:,:)    !constrained least squares reconstruction matrix for temperature gradient
real, allocatable :: rec_tempsqmat(:,:,:)    !constrained least squares reconstruction matrix for temperature gradient
real, allocatable :: rec_uleft(:,:,:,:)    !boundary extrapolated value for main equations variables for considered cell
real, allocatable :: rec_uleft_dg(:,:,:,:)    !boundary extrapolated solution value (var, i_face, i_qp)
real, allocatable :: rec_uleftturb(:,:,:,:)    !boundary extrapolated value for turbulent equations variables for considered cell
real, allocatable :: rec_uleftturbv(:,:,:,:,:)    !boundary extrapolated values for turbulent equations gradients for considered cell
real, allocatable :: rec_uleftv(:,:,:,:,:)    !boundary extrapolated values for main equations gradients for considered cell
real, allocatable :: rec_uleftx(:,:,:,:)    !boundary extrapolated value for main equations variables for considered cell
real, allocatable :: rec_velinvlsqmat(:,:,:)    !constrained least squares reconstruction matrix for velocity gradient
real, allocatable :: rec_vellsq(:,:,:)    !constrained least squares reconstruction matrix for velocity gradient
real, allocatable :: rec_velocitydof_wall(:,:,:,:)    !unlimited gradients for velocity for wall cells
real, allocatable :: rec_vext_ref(:,:)    !reference coordinates of the vertex by which the transformation has been based upon
real, allocatable :: rec_volume(:,:,:)    !volume of elements in the stencil
real, allocatable :: rec_volume_w(:,:,:)    !volume of elements in the stencil
real, allocatable :: rec_volumec(:,:,:)    !volume of elements in the stencil
real, allocatable :: rec_wallcoeff(:,:)    !constrained least squares gaussian elimination component
real, allocatable :: rec_wallcoefg(:,:)    !constrained least squares gaussian elimination component
real, allocatable :: rec_weightl(:,:,:)
real, allocatable :: rec_weno(:,:,:)    !weno weights for main equations variables
real, allocatable :: rec_weno2(:,:,:)    !weno weights for turbulent equations variables
real, allocatable :: rec_wenos(:,:,:,:,:)    !weno weights for main equations with respect to characteristic variables


!==================== ibound_ ====================

integer, allocatable :: ibound_ishape(:)    !indices for the shape, bc code, number, and which face is bounded
integer, allocatable :: ibound_icode(:)    !indices for the shape, bc code, number, and which face is bounded
integer, allocatable :: ibound_inum(:)    !indices for the shape, bc code, number, and which face is bounded
integer, allocatable :: ibound_which(:)    !indices for the shape, bc code, number, and which face is bounded
integer, allocatable :: ibound_face(:)    !indices for the shape, bc code, number, and which face is bounded
integer, allocatable :: ibound_ibid(:)    !indices for the shape, bc code, number, and which face is bounded
integer, allocatable :: ibound_nibl(:)      ! length of ibl(:) for each i
integer, allocatable :: ibound_nlocal(:)    ! length of localn(:)/cpun(:) for each i
integer, allocatable :: ibound_ibl(:,:)    !bc flag; (max_ibl,   nbound)
integer, allocatable :: ibound_localn(:,:)    !local number and cpu for each bounded surface; (max_local, nbound)
integer, allocatable :: ibound_cpun(:,:)    !local number and cpu for each bounded surface; (max_local, nbound)



#ifdef xpu
!$omp declare target (aa_1, adda_1, adda_1_s, adda_2, adda_2_s, adda_alpha_1, adda_alpha_2, adda_filter_strong, adda_filter_weak, adda_type)
!$omp declare target (alpha_0, alpha_inf1, alpha_inf2, alpha_star0, alpha_starinf, angle_per, average_restart, beta_i1, beta_i2, beta_starinf)
!$omp declare target (beta_t, bleed_end, bleed_number, bleed_plenum, bleed_porosity, bleed_start, bleed_type, bound_len, bound_offset, boundhir_dg)
!$omp declare target (br2_damping, br2_yn, bubble_centre, bubble_radius, c_des_sa, c_des_sst, c_mu_inlet, c_sas, c_smg, catalytic_con)
!$omp declare target (catalytic_wall, chunk_n, code_profile, d_corr, des_model, ek_time, el_bnd, el_int, eta2_sas, every_time)
!$omp declare target (extended_bounds, fastest_q, fil_alpha, fil_nc, fil_s, filter_type, gamma_in, halo_offset, i_turb_inlet, ibound_cpun)
!$omp declare target (ibound_face, ibound_ibid, ibound_ibl, ibound_icode, ibound_inum, ibound_ishape, ibound_localn, ibound_nibl, ibound_nlocal, ibound_t)
!$omp declare target (ibound_t2, ibound_which, ielem_admis, ielem_avars, ielem_bleedn, ielem_condition, ielem_condx, ielem_dih, ielem_dih2, ielem_diss)
!$omp declare target (ielem_dtl, ielem_er, ielem_er1, ielem_er1dt, ielem_er1er2, ielem_er2, ielem_er2dt, ielem_erx, ielem_faceanglex, ielem_faceangley)
!$omp declare target (ielem_facediss, ielem_filtered, ielem_full, ielem_ggs, ielem_hybrid, ielem_ibounds, ielem_idegfree, ielem_ifca, ielem_ihex, ielem_ihexgl)
!$omp declare target (ielem_indexf, ielem_indexi, ielem_ineigh, ielem_ineighb, ielem_ineighg, ielem_ineighn, ielem_inter_id, ielem_interior, ielem_inumneighbours, ielem_iorder)
!$omp declare target (ielem_ishape, ielem_itotalpoints, ielem_linc, ielem_lwcx2, ielem_minedge, ielem_mode, ielem_mood, ielem_mood_o, ielem_nodes, ielem_nodes_faces)
!$omp declare target (ielem_nodes_faces_v, ielem_nodes_neighbours, ielem_nodes_v, ielem_nofbc, ielem_nojecount, ielem_nonodes, ielem_q_face_q_mapl, ielem_qface, ielem_recalc, ielem_reduce)
!$omp declare target (ielem_reorient, ielem_stencil_dist, ielem_surf, ielem_totvolume, ielem_troubled, ielem_types_faces, ielem_vdec, ielem_viscx, ielem_vortex, ielem_walldist, ielem_walltrans)
!$omp declare target (ielem_walls, a405, nof_perturbations405, ielem_wcx, ielem_xxc, ielem_yyc, ielem_zzc, impdiag_mf, impoff_mf, indicator_type, init_mu_ratio, inoder4_bct)
!$omp declare target (inoder4_cord, inoder4_itor, integ_basis_dg_value, integ_basis_value, integ_basis_valuec, ires_turb, ires_unsteady, jump_cond1, jump_cond2, jump_cond3)
!$omp declare target (kappa_sst, kinit_srf, l_turb_inlet, lmach_style, m_1_val, m_t0, max_faces, max_fnodes, max_nodes, modal_filter)
!$omp declare target (modal_filter_strong, modal_filter_weak, mood_mode, mood_var1, mood_var2, mood_var3, mood_var4, mp_a_in)
!$omp declare target (mp_janaf, mp_m, mp_modelc, mp_pinf, mp_r_in, mp_thigh_in, mp_tlow_in, mp_tmid_in, mrf_rot_gl)
!$omp declare target (n_boundaries, nodes_i, nodes_part, nof_bounded, nof_bubbles, nof_interior, nof_species, nof_variables, num_dg_dofs, num_dg_reconstruct_dofs)
!$omp declare target (num_hexas, num_prisms, num_pyramids, num_tetras, out_time, output_freq, part1_end, part2_end, part3_end, part4_end)
!$omp declare target (part5_end, per_rot, point1_gl, point2_gl, pr_t1, pr_t2, pr_t3, pr_t4, pr_t5, pr_t6)
!$omp declare target (pr_t7, pr_t8, prace_t1, prace_t2, prace_t3, prace_t4, prace_t5, prace_t6, prace_t7, prace_t8)
!$omp declare target (prace_t9, prace_tx1, prace_tx2, prace_tx3, press_outlet, mach_outlet_target, prev_turbmodel, qp_array_qp_weight, qp_array_x, qp_array_y, qp_array_z)
!$omp declare target (qp_hexa, qp_line, qp_line_n, qp_prism, qp_pyra, qp_quad, qp_quad_n, qp_tetra, qp_triangle, qp_triangle_n)
!$omp declare target (qsas_model, r_beta, r_gas, r_k_sst, r_om_sst, radius_gl, rec_br2_aux_var, rec_br2_local_lift, rec_cgradientstemp, rec_cond)
!$omp declare target (rec_findw, rec_g0, rec_gradf, rec_gradients, rec_gradients2, rec_gradientsc, rec_gradientsc2, rec_gradientstemp, rec_gradientstemp_wall, rec_gradientsturb)
!$omp declare target (rec_gradientsturb_wall, rec_grads, rec_gradsav, rec_ihexb, rec_ihexbc, rec_ihexg, rec_ihexgc, rec_ihexl, rec_ihexlc, rec_ihexn)
!$omp declare target (rec_ihexnc, rec_indicator, rec_indicatorc, rec_invccjac, rec_invctjac, rec_invmat_stencilt, rec_invmat_stenciltc, rec_k0, rec_local, rec_mrf)
!$omp declare target (rec_mrf_origin, rec_mrf_velocity, rec_periodicflag, rec_qpoints, rec_qpoints_p, rec_rotvel, rec_rpoints, rec_stencils, rec_stencilsc, rec_surf_qpoints, rec_tempsq)
!$omp declare target (rec_tempsqmat, rec_uleft, rec_uleft_dg, rec_uleftturb, rec_uleftturbv, rec_uleftv, rec_uleftx, rec_velinvlsqmat, rec_vellsq, rec_velocitydof_wall)
!$omp declare target (rec_vext_ref, res_sum, pos_l1, pos_l2, pos_l3, pos_l4, total_pressure_inlet, total_temperature_inlet, density_inlet, ipos_l1, ipos_l2, rec_volume, rec_volume_w, rec_volumec, rec_wall, rec_wallcoeff, rec_wallcoefg, rec_weightl, rec_weno, rec_weno2)
!$omp declare target (rec_wenos, reduce_comp, res_time, rg_hzero, rg_kf_type, rg_molm, rg_nof_reactions, rg_nof_tv_coef, rg_relax, rg_t_inf)
!$omp declare target (rg_t_ref, rg_t_wall_init, rg_thetag, rg_ttr, rg_tv_coef, rg_tve, rg_vf, rgs_ab, rgs_bb, rgs_cb)
!$omp declare target (rgs_eps_over_k, rgs_mg, rgs_sigmaa, rhs_sol_mm_dg, rhs_val, rhs_valdg, rhst_val, rot_corr, schmidt_lam, schmidt_turb)
!$omp declare target (sht_rg, sigma_k1, sigma_k2, sigma_om1, sigma_om2, sigma_phi, source_active, srf_origin, srf_velocity)
!$omp declare target (stencil_io, taylor_ens, taylor_ensx, temp_model, thread_n, tol_per, typ_countn, typ_countn_global, typ_countn_global_w)
!$omp declare target (typ_countn_w, u_c_br2_aux_var, u_c_rms, u_c_val, u_c_valdg, u_cs_val, u_cs_valdg, u_ct_val, u_cw_val, u_cw_valdg)
!$omp declare target (u_e_val, v_ref, viscous_s, vort_model, wall_temp, weight_lsqr)
!$omp declare target (wnodes_part, wpart1_end, wpart2_end, wpart3_end, wpart4_end, wpart5_end, write_variables, write_variables_av, write_variables_av_w, write_variables_w)
!$omp declare target (zero_turb_init, zeta_star, adda, allnodesgloball, allres, allresdt, alls, alpha, aoa, averaging)
!$omp declare target (transition_model, transition_axis, transition_direction, transition_ramp_type, transition_location, transition_ramp_length)
!$omp declare target (beta, betaas, binio, bleed, boundtype, cascade, cavitation, origin)
!----mpi comm
!$omp declare target (halo_len, halos_len, halo_offset, halos_offset, ineedbound, ineedbounds)
!$omp declare target (bound_len, bounds_len, need_side, need_q, need_loc, halo_proc, halos_proc)
!$omp declare target (bound_proc, bounds_proc, bound_offset, bounds_offset, solhir, solhis, solhird)
!$omp declare target (solhisd, boundhiri, boundhisi, boundhir, boundhis, boundhir_dg, boundhis_dg, boundhirm)
!$omp declare target (boundhism, bound_total, bounds_total, halo_total, halos_total, solhi_loc)
!---mpi comm
!$omp declare target (cb1, cb2, cfl, cflmax, cflramp, cfw, charlength, ct1, ct2, ct3)
!$omp declare target (ct4, cv1, cw1, cw2, cw3, datatypeint, datatypex, datatypexx, datatypey, datatypeyy)
!$omp declare target (datatypez, dg, dg2fv, dimensiona, dims, dt, ees, emetis, extf, fastest)
!$omp declare target (fastmovie, filtering, firstorder, firstrese, firstresk, firstresomega, firstrespass, firstresr, firstrest, firstresu)
!$omp declare target (firstresv, firstresw, forcex, forcey, forcez, gamma, governingequations, greengo, gridar1)
!$omp declare target (gridar2, guassianquadra, hybridist, iadapt, ibcode, iboundary, ibside, icarlos1, icarlos2, icompact)
!$omp declare target (icong, iconimp, iconsgvq, iconsr, icoupleturb, idegfree, idegfree2, idegfree3, ievery, ievery2)
!$omp declare target (ieveryav, iforce, igianagraps, igqrules, ihax1, ihybrid, iloop, iloopx, ilx, imaxb)
!$omp declare target (imaxdegfree, imaxdegfree2, imaxe, imaxn, impdiag, impdiagt, impdu, impdu_old, impoff, impofft, in)
!$omp declare target (ineedbound, ineedhalo, initcond, initialres, inum2, inwhichel, iorder, iorder2, ioverst, ioverto)
!$omp declare target (iperiodicity, iriemann, irs, ischeme, iscoun, iselem, ispal, isplit, issf, istn)
!$omp declare target (it, itestcase, itold, itotalb, itt, ivortex, iweightlsqr, iweno, iwmaxe, jk)
!$omp declare target (jtotal, jtotal1, jtotal2, jtotal3, kappa, kdum1, kdum2, kdum3, kill, kloopx, jtot, bound_total, bounds_total, halo_total, halos_total)
!$omp declare target (kmaxn, l0norm, l1norm, l2norm, lam, lamps, lamx, lamy, lamz, limiter)
!$omp declare target (lmach, lowmem, lowmemory, lwci1, Mach_in, modeio, momentx, momenty, momentz, mood, movement)
!$omp declare target (mrf, multispecies, nderivative, nodelist, nprobes, nproc, nrotors, ntmax, numberofpoints, numberofpoints2)
!$omp declare target (numneighbours, numneighbours2, oo2, outsurf, passivescalar, pi, poly, prandtl, pres, prevres)
!$omp declare target (prtu, qrde, realgas, relax, required, rescounter, rescountert, residualfreq, reslimit, resmax)
!$omp declare target (resmaxt, restart, reynolds, rframe, rhc1, rhc2, rhc3, rhc4, rres, rungekutta)
!$omp declare target (scaler, sht, sigma, spatialorder, spatiladiscret, spkin, spos, srfg)
!$omp declare target (statistics, stennorm, subdiv, surfshear, suther, swirl, t, taylor, tecplot)
!$omp declare target (temporder, thermal, timestep, tolbig, tolsmall, totalvolume, totiw, totwalls, totwallsc, turbinit)
!$omp declare target (turbulence, turbulenceequations, turbulencemodel, twall, typesten, tz1, ufreestream, unwou, upperlimit, upturblimit)
!$omp declare target (uvel, vectorx, pos_l, pos_g, ipos_l, ipos_g, vectory, vectorz, visc, voll, vorder, vvel, wallc, wdatatypeint)
!$omp declare target (wdatatypex, wdatatypexx, wdatatypey, wdatatypeyy, wdatatypez, wenocentralweight, wenocnschar, wenoz, wenwrt, wkdum1)
!$omp declare target (wkdum2, wkdum3, wvel, xmpielrank, xper, yper, zero, zper)
!$omp declare target (indicator_par1, indicator_par2, indicator_par3, totk, totens, totensx, kill_nan)
!$omp declare target (solhir_flat, solhis_flat, boundhir_flat, boundhis_flat, boundhir_dgflat, boundhis_dgflat, boundhiri_flat, boundhisi_flat, weights_t, weights_q, weights_l, ilength1, ilength2, ccfl, ind1)

#endif






end module declaration
