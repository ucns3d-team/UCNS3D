module parameters
use mpiinfo
use declaration

implicit none

 contains


subroutine read_ucns3d
!> @brief
!> this subroutine reads the parameter file

	implicit none

 	integer :: inv,ix,ibleed,i,j
 	integer :: inv1
 	real :: angledum
	character(48)::stamp1,frame
	logical::here1,here2,here3,here5,here,here4,here7,here8,here9,bleedio,here10
	character(len=8)  :: date_w
	character(len=10) :: time_w

	mp_modelc=0	!allaire by default


 	
	call mpi_barrier(mpi_comm_world,ierror)
	inquire (file='RESTARTAV.dat',exist=here1)
	if (here1) then
	average_restart=1
	else
	average_restart=0
	end if
	source_active=0


	movement=0
	
 	frame='ROTFRAME.dat'
	inquire (file=frame,exist=here)
	if (here) then
	open(16,file=frame,form='formatted',status='old',action='read')
	read(16,*)!1
	read(16,*)!2
	read(16,*)!3
	read(16,*)!4
	read(16,*)rframe
	read(16,*)!6
	read(16,*)srf_origin(1),srf_origin(2),srf_origin(3)
	read(16,*)!8
	read(16,*)srf_velocity(1),srf_velocity(2),srf_velocity(3)
    read(16,*)!10    
	read(16,*)per_rot,angle_per,v_ref	
    read(16,*)!12
    read(16,*)nrotors
    allocate(point1_gl(nrotors,3),point2_gl(nrotors,3),radius_gl(nrotors),mrf_rot_gl(nrotors))
    do inv1=1,nrotors !storing multiple rotating frame coordinates
        read(16,*)point1_gl(inv1,1),point1_gl(inv1,2),point1_gl(inv1,3)
        read(16,*)point2_gl(inv1,1),point2_gl(inv1,2),point2_gl(inv1,3)
        read(16,*)radius_gl(inv1), mrf_rot_gl(inv1)
	end do
	close(16)
        if(per_rot.eq.1)then
	        tol_per=1.0e-8
            lowmem=1
            iperiodicity = 1 
        end if
        if(rframe.eq.2)then
            mrf=1
            srfg=0
        end if
        if (rframe.eq.1)then
            srfg=1
            mrf=0
        end if
	else
        rframe=0
        srfg=0
        mrf=0
	end if
	
	if ((mrf.eq.1).or.(srfg.eq.1))then
	 source_active=1
    kinit_srf=0.00001
    if (mrf.eq.1)then
            rot_corr=1
            d_corr=1
    end if
    if (srfg.eq.1)then
            rot_corr=0
            d_corr=1
    end if
	
	end if




	
	inquire (file='MULTISPECIES.DAT',exist=here2)
	if (here2) then
	multispecies=1
	open(14,file='MULTISPECIES.DAT',form='formatted',status='old',action='read')
	read(14,*)
	read(14,*)
	read(14,*)nof_species
    allocate(gamma_in(1:nof_species),mp_a_in(1:nof_species),mp_r_in(1:nof_species),mp_pinf(1:nof_species))
    read(14,*)gamma_in(1:nof_species)
    read(14,*)mp_a_in(1:nof_species)
    read(14,*)mp_r_in(1:nof_species)
    read(14,*)mp_pinf(1:nof_species)
    close(14)
	else
	multispecies=0
	end if

	inquire (file='MULTISPECIES_DIFF.DAT',exist=here2)
	if (here2) then
	open(14,file='MULTISPECIES_DIFF.DAT',form='formatted',status='old',action='read')
	read(14,*)
	read(14,*)mp_modelc			!mp model   0=allaire, 1=mp_species
    allocate(mp_m(1:nof_species),mp_brok_a(1:nof_species), mp_brok_b(1:nof_species), mp_brok_c(1:nof_species),mp_tlow_in(1:nof_species), mp_thigh_in(1:nof_species),mp_tmid_in(1:nof_species),mp_janaf(1:nof_species,1:2,1:5))
    read(14,*)mp_m(1:nof_species)
    read(14,*)mp_brok_a(1:nof_species)
    read(14,*)mp_brok_b(1:nof_species)
    read(14,*)mp_brok_c(1:nof_species)
    read(14,*)mp_tlow_in(1:nof_species)
    read(14,*)mp_thigh_in(1:nof_species)
    read(14,*)mp_tmid_in(1:nof_species)
    do i=1,5
    do j=1,2
    read(14,*)mp_janaf(1:nof_species,j,i)
    end do
    end do
    close(14)
	else

	end if




	inquire (file='REALGAS.DAT',exist=here10)
	if (here10) then
	realgas=1;
	open(14,file='REALGAS.DAT',form='formatted',status='old',action='read')
	read(14,*)
	read(14,*)
	read(14,*)nof_species
	read(14,*)rg_nof_reactions
	read(14,*)rg_kf_type
	read(14,*)rg_relax
    allocate(rg_vf(1:nof_species), rg_molm(1:nof_species),rg_hzero(1:nof_species),rg_thetag(1:nof_species))
    read(14,*)rg_vf(1:nof_species)
	read(14,*)rg_molm(1:nof_species)
	read(14,*)rg_hzero(1:nof_species)
	read(14,*)rg_thetag(1:nof_species)
    read(14,*)rg_t_inf
    read(14,*)rg_t_wall_init
    read(14,*)rg_t_ref
	read(14,*)rg_nof_tv_coef
	rg_ttr=rg_t_inf;	!need to specify initial values for those
	rg_tve=rg_t_inf;	!need to specify initial values for these as well




	if (rg_nof_tv_coef.gt.0)then
		allocate(rg_tv_coef(1:3,1:rg_nof_tv_coef))
		read(14,*)rg_tv_coef(1,1:rg_nof_tv_coef)
		read(14,*)rg_tv_coef(2,1:rg_nof_tv_coef)
		read(14,*)rg_tv_coef(3,1:rg_nof_tv_coef)
	end if




	allocate(rgs_mg(1:nof_species))
	rgs_mg=rg_molm * 1.0e3   ! g/mol for d_ij correlation





    allocate(rgs_sigmaa(1:nof_species),rgs_eps_over_k(1:nof_species),rgs_ab(1:nof_species),rgs_bb(1:nof_species),rgs_cb(1:nof_species))

    read(14,*)rgs_ab(1:nof_species)
    read(14,*)rgs_bb(1:nof_species)
	read(14,*)rgs_cb(1:nof_species)
	read(14,*)rgs_sigmaa(1:nof_species)
	read(14,*)rgs_eps_over_k(1:nof_species)



	read(14,*) catalytic_wall

	if (catalytic_wall.eq.1)then
	allocate(catalytic_con(1:nof_species))
	read(14,*)catalytic_con(1:nof_species)
	end if

    close(14)
	else
	realgas=0
	end if









	inquire (file='THERMAL.DAT',exist=here9)
	if (here9) then
	thermal=1
	open(31,file='THERMAL.DAT',form='formatted',status='old',action='read')
	read(31,*)temp_model !1=adiabatic,2=isothermal 3
    read(31,*)wall_temp
    close(31)
	else
	thermal=0
	end if



	inquire (file='MOOD.DAT',exist=here3)
	if (here3) then
	mood=1
	open(17,file='MOOD.DAT',form='formatted',status='old',action='read')
	read(17,*)
	read(17,*)
	read(17,*)mood_mode        !type of mood mode (1=relaxed, 0=original)
    read(17,*)mood_var1,mood_var2
    read(17,*)mood_var3,mood_var4
    if (n.eq.0)then
		print*,"mood active"
    end if
    close(17)
	else
	mood=0
	end if
	
	
	
	inquire (file='INDICATOR.DAT',exist=here4)
	if (here4) then
	open(18,file='INDICATOR.DAT',form='formatted',status='old',action='read')
	read(18,*)
	read(18,*)
	read(18,*)indicator_type        !type of indicator (1=mood,2=shu,3=...)
    read(18,*)indicator_par1        !parameter 1    
    read(18,*)indicator_par2        !parameter 2    
    read(18,*)indicator_par3        !parameter 3
    close(18)
	else
	indicator_type=11
	indicator_par1=10e-12;indicator_par2=10e-10;indicator_par3=zero;
	
	end if
			
	
	inquire (file='OUTLETS.DAT',exist=here9)
	if (here9) then
	open(29,file='OUTLETS.DAT',form='formatted',status='old',action='read')
	read(29,*)
	read(29,*)press_outlet
	print*,"i am reading the pressure for the outlet vents from the file"
    close(29)	
	end if



	inquire (file='FILTER.DAT',exist=here7)
	if (here7) then
	filtering=1
	open(19,file='FILTER.DAT',form='formatted',status='old',action='read')
	read(19,*)
	read(19,*)
	read(19,*)filter_type        !type of filter (1=exponential)
    read(19,*)fil_alpha
    read(19,*)fil_s
    read(19,*)fil_nc
    close(19)
	else
	filtering=0
	end if



	inquire (file='ADDA.DAT',exist=here8)
	if (here8) then
	adda=1
	open(19,file='ADDA.DAT',form='formatted',status='old',action='read')
	read(19,*)
	read(19,*)
	read(19,*)adda_type        !type of filter (1=exponential)
    read(19,*)adda_alpha_1,adda_alpha_2
    read(19,*)adda_1_s,adda_2_s
    read(19,*)adda_1,adda_2
    close(19)
	else
	adda=0
	adda_type=0
	end if


	
	
	call mpi_barrier(mpi_comm_world,ierror)
	
	
	
	
	
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)dimensiona,statistics,code_profile
	read(15,*)
	read(15,*)
	read(15,*)governingequations,initcond
	read(15,*)
	read(15,*)
	read(15,*)turbulence,icoupleturb,passivescalar
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)rres,ufreestream,vvel,wvel,pres
	read(15,*)
	read(15,*)
	read(15,*)aoa,vectorx,vectory,vectorz
	read(15,*)
	read(15,*)
	read(15,*)gamma,prandtl,reynolds,charlength
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)spatiladiscret,iriemann,spatialorder,limiter,poly
	read(15,*)
	read(15,*)
	read(15,*)wenocnschar,ees,wenoz,wenocentralweight
	read(15,*)
	read(15,*)
	read(15,*)temporder,cfl,timestep,upperlimit,reslimit
	read(15,*)
	read(15,*)
	read(15,*)iboundary,boundtype,scaler
	read(15,*)
	read(15,*)
	read(15,*)greengo,lmach
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)out_time,ntmax,wallc
	read(15,*)
	read(15,*)
	if (code_profile.lt.0)then
	read(15,*)tecplot,output_freq,ievery2,ieveryav,stencil_io
		ievery=10e15
	else					!steady state
		read(15,*)tecplot,ievery,ievery2,ieveryav,stencil_io
		output_freq=10e15
	end if
	read(15,*)
	read(15,*)
	read(15,*)averaging
	read(15,*)
	read(15,*)
	read(15,*)outsurf,iforce,surfshear
	read(15,*)
	read(15,*)
	read(15,*)ires_turb,ires_unsteady,lamps,prev_turbmodel  
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)
	read(15,*)nprobes
	read(15,*)

	    
	if (dimensiona.eq.3)then
		ccfl=(cfl/3.0d0)
	else
		ccfl=(cfl/2.0d0)
	end if
	if (rungekutta.eq.4)then
 	      ind1=7
	else
 	      ind1=5
	end if
	origin(1:3)=0.0d0


	fastmovie=0
	dg=0;
	!turbulence default values
	
	
			! 	turbulence model parameters:
			! ---options spalart almaras---
			ispal=1 !    ||spalart allmaras version:| 1:original |2: negative modification
			!des_model=0! 0              			|| 1-detached eddy simulation 2-delayed des
			! ---constants-------
			cb1=0.1355	! 0.1355				|| cb1 
			cb2=0.622! 0.622   			|| cb2
			sigma=0.666666667! 0.666666667      		|| sigma
			kappa=0.41  ! 0.41     			|| kappa
			cw1=3.23886781677! 3.23886781677      		|| cw1 (0 for (cb1/kappa**2)+1/sigma(1+cb2) or predifined)
			cw2=0.3! 0.3      			|| cw2
			cw3=2.0 ! 2.0      			|| cw3
			cv1=7.1! 7.1      			|| cv1
			ct1=1.0 ! 1.0      			|| ct1
			ct2=2.0! 2.0      			|| ct2
			ct3=1.1 ! 1.1      			|| ct3
			ct4=2.0 ! 2.0      			|| ct4
			prtu=0.9! 0.9				|| turbulent prandtl number
			twall=0! 0				|| wall temperature (kelvin) leave 0 for adiabatic (q_wall =0 <=> dt/dn=0)
			turbinit=3.0 ! 3.0	  			|| initial value for turbulence parameter (multiplyied by the freestream viscosity from given re)
			upturblimit=1000000! 1000000				|| upper limit for turbulence
			residualfreq=100! 10				|| residual compute every
			irs=0! 0				||implicit residual smoothing (doubles cfl)
			c_des_sa=0.61	! 0.61				||c_des_sa
			! =============================================================================
			! k-omega sst: 
			! ---options---
			vort_model=0! 0					||0:default strain-production for k, 1:vorticity-production for k
			qsas_model=0! 0					||hybrid models: 1-sas-sst  //  2-des-sst 
			zero_turb_init=0! 0					||zero turbulence option (if 1, initialization for k/w is done with zero turbulence)
			! --------------k-omega constants--------------------------
			sigma_k1=1.176470588! 1.176470588				||sigma_k1
			sigma_k2=1.0! 1.0					||sigma_k2
			sigma_om1=2.0	! 2.0					||sigma_om1
			sigma_om2=1.168! 1.168					||sigma_om2
			aa_1=0.31! 0.31					||aa_1
			beta_i1=0.075	! 0.075					||beta_i1
			beta_i2=0.0828! 0.0828					||beta_i2
			alpha_starinf=1.0! 1.0					||alpha_starinf
			alpha_0=0.111111111! 0.111111111				||alpha_0
			beta_starinf=0.09! 0.09					||beta_starinf
			r_beta=8.0! 8.0					||r_beta
			r_k_sst=6.0! 6.0					||r_k_sst
			beta_t=0.072	! 0.072					||beta_t
			kappa_sst=0.41! 0.41					||kappa   (same in spalart-allmaras)
			r_om_sst=2.95	! 2.95					||r_om_sst
			zeta_star= 1.5	! 1.5					||zeta_star  (only for mach corrections)
			m_t0=0.25! 0.25					||m_t0		(only for mach corrections)
			c_mu_inlet=0.09! 0.09					||c_mu_inlet
			c_smg=0.11! 0.11					||c_smg  (sas)
			eta2_sas=3.51! 3.51					||eta2_sas (sas)
			sigma_phi=0.666666667	! 0.666666667				||sigma_phi (sas)
			c_sas=2.0! 2.0					||c_sas   (sas)
			l_turb_inlet=0.0026! 0.0026					||l_turb_inlet (default turbulence lengthscale)
			i_turb_inlet=0.1! 0.1					||i_turb_inlet (default turbulene intensity)
			init_mu_ratio=0.01! 0.01					||init_mu_ratio	(default initial ratio for turbulent viscosity)
			c_des_sst=0.61! 0.61					||c_des_sst
			! =============================================================================
			! passive scalar transport
			schmidt_lam=10.0! 10.0					||laminar schmidt number
			schmidt_turb=0.7! 0.7					||turbulent schmidt number
	
	cavitation=0
	    extended_bounds=0
	    
	    select case(code_profile)
	
	
	case (0,888)       
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if ((initcond.eq.405).or.(initcond.eq.405))then
    iadapt=1
    end if
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10000	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=10000	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=0.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	


	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 
	 
	 des_model=0


	 case (199)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if ((initcond.eq.405).or.(initcond.eq.405))then
    iadapt=1
    end if
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=2	!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=300000000	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=300000000	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
        d_corr=1
        turbinit=0.1
	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0




	 case (501)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=20.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=50.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)


	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0
	 



	case (30)  !vertex based multidimensional-experimental option not for production     
	
		lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
		binio=1	    	!i/o (ascii=0, binary=1) 
		lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
		reduce_comp=0	!quadrature free flux=0 not true,1 true
		turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
		ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
		hybridist=0.0d0 !hybrid distance
		swirl=0		!swirling flow:0 deactivated, 1 activated
		iadapt=0	!adaptive numerical scheme (0 not true,1 true)
		if (initcond.eq.405)iadapt=1
		icompact=0	!compact stencil mode(0 not true,1 true)
		extf=3		!stencils stability values from 1.2 to 3 (default 2)
		weight_lsqr=0	!weighted least squares(0 not true,1 true)
		guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
		fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
			relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
		cflmax=30	!cflmax:to be used with ramping
		cflramp=0	!cfl ramping: |0: deactivated |1:activated
		emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
		itold=10000	!tolerance=n_iterations
		gridar1=20.0	! 0	  5.0    7.0  limit aspect ratio cells,
		gridar2=50.0	! limit volume cells
		fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
		lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
		lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
		
		
		if (iboundary.eq.1)then
		 lowmem=1
		 end if
		 
		 
		 des_model=0





	 case (8)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=2		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=100.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=200.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)


	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0
	 
	 
	  case (88)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=40.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=40.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)

	fastmovie=0
		
	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0

	 case (18)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=1	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=40000.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=40000.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=0.0d0;lamy=0.0d0;lamz=0.0d0	!linear advection coefficients (lamx, lamy,lamz)

	fastmovie=0

	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0


	 
	 case (98)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=2		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=1	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=20.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=40.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)

	fastmovie=0

	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0
	 
	 
	 case (9)  !robust
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3	!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=7.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	fastmovie=0
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 
	 
	 des_model=0





         case (777)  !robust

        lowmemory=0     !memory usage: |0: high(faster) |1:low (slower)||
        binio=1         !i/o (ascii=0, binary=1)
        lowmem=0        !global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
        reduce_comp=0   !quadrature free flux=0 not true,1 true
        turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
!       icoupleturb=0   !coupling turbulence model: |1:coupled | 0: decoupled
        ihybrid=0       !hybrid turbulence : |1:enabled|0:disabled
        hybridist=0.0d0 !hybrid distance
        swirl=0         !swirling flow:0 deactivated, 1 activated
        iadapt=1        !adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
        icompact=0      !compact stencil mode(0 not true,1 true)
        extf=3  !stencils stability values from 1.2 to 3 (default 2)
        weight_lsqr=0   !weighted least squares(0 not true,1 true)
        guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
        fastest_q=1     !store gqp points (1 =yes faster, 0= slower)
        relax=1         !relaxation parameter : |1:block jacobi |2: lu-sgs
        cflmax=30       !cflmax:to be used with ramping
        cflramp=0       !cfl ramping: |0: deactivated |1:activated
        emetis=6        !metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
        itold=10000     !tolerance=n_iterations
        gridar1=1000000.0    ! 0       5.0    7.0  limit aspect ratio cells,
        gridar2=7000000.0     ! limit volume cells
        fastest=0       ! 0                             ||fastest, no coordinate mapping (1: engaged,0:with transformation)
        lmach_style=0   !0                      ||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
        lamx=0.0d0;lamy=1.0d0;lamz=1.0d0        !linear advection coefficients (lamx, lamy,lamz)
        fastmovie=0
        if (iboundary.eq.1)then
         lowmem=1
         end if


         des_model=0
	 
	 
	 case (91)  !robust with matrix free lu-sgs
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3	!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=1	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=3	!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10000.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=70000.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 
	 
	 des_model=0
	 
	 
	 
	 
	case (1)           !for ddes
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	!icoupleturb=1	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=7.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	cavitation=0
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 turbinit=0.1
	 des_model=2
	 
	 
	 
	 case (11)           !for ddes with matrix free lu-sgs
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	!icoupleturb=1	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=10.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 des_model=2
	 
	 
	 
    case (100)           !for hybrid dg-fv method
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	!icoupleturb=1	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=3		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
    relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=10.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	dg=1   !0=deactivated fv only, 1=activated dg only, 2=hybrid
	
	bound_lim=0
    cavitation=0        
	
	
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 des_model=2
	 

	 br2_damping=3.0
	 br2_yn=2
	 
	 case (101)           !for hybrid dg-fv method
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	!icoupleturb=1	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=2		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=6!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=10.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	dg=1   !0=deactivated fv only, 1=activated dg only, 2=hybrid
	bound_lim=0
	cavitation=0
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 des_model=2
	 

	 br2_damping=3.0
	 br2_yn=2
	 
	 
	  case (102)           !for pure dg method
	
	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	!i/o (ascii=0, binary=1) 
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	!icoupleturb=1	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=2		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
    relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	!tolerance=n_iterations
	gridar1=10.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=10.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	!linear advection coefficients (lamx, lamy,lamz)
	dg=1   !0=deactivated fv only, 1=activated dg only, 2=hybrid
	
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	 des_model=2


	 br2_damping=3.0
	 br2_yn=2

	 
	
	
	case default
	
	lowmemory=0 	 !memory usage: |0: high(faster) |1:low (slower)|| 
	binio=1	    	 !i/o (ascii=0, binary=1) 
	lowmem=0    	 !global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	 !quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst	
	icoupleturb=0	 !coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	 !hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0  !hybrid distance
	swirl=0		 !swirling flow:0 deactivated, 1 activated
	iadapt=0	 !adaptive numerical scheme (0 not true,1 true)
	icompact=0	 !compact stencil mode(0 not true,1 true)
	extf=3		 !stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=0	 !weighted least squares(0 not true,1 true)
	guassianquadra=0 !gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	 !store gqp points (1 =yes faster, 0= slower)
        relax=2		 !relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	 !cflmax:to be used with ramping
	cflramp=0	 !cfl ramping: |0: deactivated |1:activated
	emetis=6    	 !metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets 
	itold=10000	 !tolerance=n_iterations
	gridar1=1000000.0	 ! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=1000000.0	 ! limit volume cells
	fastest=0	 ! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	 !0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=1.0d0;lamy=1.0d0;lamz=1.0d0	 !linear advection coefficients (lamx, lamy,lamz)
	ispal=1! 1				||spalart allmaras version:| 1:original |2: negative modification
	
	
	
	if (iboundary.eq.1)then
	 lowmem=1
	 end if
	
	
	des_model=0






	case (17)

	lowmemory=0 	!memory usage: |0: high(faster) |1:low (slower)||
	binio=1	    	!i/o (ascii=0, binary=1)
	lowmem=0    	!global arrays setting (0=without better suited for non periodic bound,1=with (large memory footprint))
	reduce_comp=0	!quadrature free flux=0 not true,1 true
	turbulencemodel=1 !turbulence model selection: |1:spalart-allmaras |2:k-w sst
! 	icoupleturb=0	!coupling turbulence model: |1:coupled | 0: decoupled
	ihybrid=0	!hybrid turbulence : |1:enabled|0:disabled
	hybridist=0.0d0 !hybrid distance
	swirl=0		!swirling flow:0 deactivated, 1 activated
	iadapt=0	!adaptive numerical scheme (0 not true,1 true)
    if (initcond.eq.405)iadapt=1
	icompact=0	!compact stencil mode(0 not true,1 true)
	extf=2		!stencils stability values from 1.2 to 3 (default 2)
	weight_lsqr=1	!weighted least squares(0 not true,1 true)
	guassianquadra=0!gaussian quadrature rule (1,2,5,6), default 0 will use the appropriate number
	fastest_q=1	!store gqp points (1 =yes faster, 0= slower)
        relax=1		!relaxation parameter : |1:block jacobi |2: lu-sgs
	cflmax=30	!cflmax:to be used with ramping
	cflramp=0	!cfl ramping: |0: deactivated |1:activated
	emetis=6    	!metis partitioner : 1: hybrid metis, 2:adaptive weights for hybrid grids, 3: uniform metis partionioner,4:nodal,6=parmets
	itold=10000	!tolerance=n_iterations
	gridar1=200000.0	! 0	  5.0    7.0  limit aspect ratio cells,
	gridar2=500000.0	! limit volume cells
	fastest=0	! 0		       		||fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||low mach treatment (1 activate, 0 disable),lmach_style(0=only normal component,1=all components)
	lamx=0.0d0;lamy=0.0d0;lamz=0.0d0	!linear advection coefficients (lamx, lamy,lamz)
	iweno=5

	if (iboundary.eq.1)then
	 lowmem=1
	 end if


	 des_model=0
	
	
	
	
	end select
	    
	    
	
	 if (n.eq.0)then
      open(63,file='history.txt',form='formatted',action='write',position='append')



		call date_and_time(date_w, time_w)

		write(63,'(A,1X,A2,"/",A2,"/",A4,1X,A2,":",A2)') &
			"ucns3d parallel computation start", &
				date_w(7:8), date_w(5:6), date_w(1:4), &
				time_w(1:2), time_w(3:4)







      end if

      
      
	    nproc=isize-1;vorder=iorder
      
	    !--------------------------1------------------------------!
	    !probe position allocation and reading of coordinates
	    !$omp master
	    if (nprobes.gt.0)then
	    if (dimensiona.eq.3)then
		  allocate(probec(1:nprobes,1:3))
		    do inv=1,nprobes
		      read(15,*)probec(inv,1),probec(inv,2),probec(inv,3)
		    end do
	    else
		  allocate(probec(1:nprobes,1:2))
	    
		    do inv=1,nprobes
		      read(15,*)probec(inv,1),probec(inv,2)
		    end do
	    end if
	    end if
	    !$omp end master
	     !--------------------------end 1-------------------------!
	     
	     
	     call mpi_barrier(mpi_comm_world,ierror)
	     
	     
	     !-------------------------2---------------------------------!
	     !number of variables setup
	     if (governingequations.le.2)then ! ns or euler
		  if (dimensiona.eq.3)then
		      nof_variables=5;dims=3
		  else
		      nof_variables=4;dims=2
		  end if
		  
		  
	      else ! linear advection
		  nof_variables=1
		  
		  if (dimensiona.eq.3)then
		    dims=3
		    else
		    dims=2
		    end if
	      end if
	      
	      !multiphase modification starts
	      if (multispecies.eq.1)then
			if (mp_modelc.eq.0)then
					if (dimensiona.eq.3)then
					nof_variables=5+nof_species+(nof_species-1);dims=3
					else
					nof_variables=4+nof_species+(nof_species-1);dims=2
					end if
			else
					if (dimensiona.eq.3)then
					nof_variables=5+nof_species;dims=3
					else
					nof_variables=4+nof_species;dims=2
					end if


			end if

		  end if
		  !multiphase modification ends

		  if (realgas.eq.1)then
					if (dimensiona.eq.3)then
					nof_variables=dimensiona+3+nof_species;dims=3
					else
					nof_variables=dimensiona+3+nof_species;dims=2
					end if

		  end if



		  
	    !--------------------------end 2-------------------------!


	    !-------------------------3---------------------------------!
	     !turbulence
	      if (turbulence .eq.1) then
		      if (turbulencemodel .eq. 1) then
			      turbulenceequations=1
		      else if (turbulencemodel .eq. 2) then
			      turbulenceequations=2
		      end if 
	      else
		      turbulenceequations=0
		      turbulencemodel=0
	      end if
	    !------------------------- end 3---------------------------------!
	   
	   
! 	   !-------------------------4---------------------------------!
	   !equations type
	   
	   select case(governingequations)
	   case(-1)
	   !multi-phase inviscid euler equations
	    if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'multi-phase inviscid euler solver engaged'
	      close(63)
	    end if
	    itestcase = 3;ivortex = 0
	   
	   case(1)
	   !navier stokes equations
	    if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'navier-stokes solver engaged'
	      close(63)
	    end if
	    itestcase = 4;ivortex = 1
	    
	  case(2)
	   !euler equations
	    if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'euler solver engaged'
	      close(63)
	    end if
	    itestcase = 3;ivortex = 0
	    
	    
	    
	  case(3)
	   !linear advection equation
	    if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'linear solver engaged sinewave'
	      close(63)
	    end if
	    itestcase = 1;ivortex = 0
	    
	 case(4)
	   !linear advection equation
	    if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'linear solver engaged step function'
	      close(63)
	    end if
	    itestcase = 2;ivortex = 0
	    
	  end select
	  
	  
	  
	  !-------------------------end equations 4---------------------------------!
	  
	  !-------------------------5---------------------------------!
	  !flow parameters
	  unwou = 3
	  betaas=1.5d0
	  suther=0.412158681d0
	  uvel=ufreestream
          if (uvel.lt.10.0)then
          r_gas=1.0
          else
                  r_gas=287.052874d0
                  suther=110.0/(pres/(rres*r_gas))

          end if
	  
	  if (initcond.eq.977)then
	  uvel=0.0d0;vvel=0.0d0;
	  ufreestream=wvel
	  end if
	  
	  
	   ! set pressure
	  if ( pres .lt. 0 ) pres = rres/gamma	
	  ! set dynamic free-stream viscosity
	  
	  if (rframe.eq.0) then
        visc = (rres*ufreestream*charlength)/reynolds
    else
        visc = (rres*v_ref*charlength)/reynolds
    end if
	  
	  if (swirl.eq.1)then
	  uvel=zero
	  end if
	  if (aoa .ne. 0.0d0) then
	  angledum=(aoa*pi)/180.0d0
	  uvel = cos(angledum)*ufreestream*vectorx
	  wvel = sin(angledum)*ufreestream*vectorz
	  vvel = sin(angledum)*ufreestream*vectory
	  end if
	  
	  if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      if (itestcase .eq. 4) then
	      if (rframe.eq.0)then
		write(63,*)'----reynolds number:',(rres*ufreestream*charlength)/visc
            else
            write(63,*)'----reynolds number:',(rres*v_ref*charlength)/visc
            
            end if
		end if
	      close(63)
	   end if
	  
	  
	  
	  !-------------------------end flow parameters 5---------------------------------!
	  
	  !-------------------------6---------------------------------!
	   !spatial & temporal discretisation
	   
	   iorder = max(1,spatialorder-1)
	   firstorder = 0
! 	   if (ees.eq.5)then
! 	   icompact=1
! 	   end if
	   iextend =2; if (icompact.eq.1)iextend = 10
	   iweno =0
	   wenwrt = 1
	   ischeme=3
	   typesten=1
	   lwci1=wenocentralweight


!  	   if (ees.eq.5)then
!   	    		lwci1=10**((wenocentralweight*14)+2)
!   	   end if



	  if (fastest.eq.1)ischeme=1
	    
	   
	   select case(spatiladiscret)



	   
	   case(1)	!no limiter
	      
	   if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'----1st order in space|linear scheme'
	      close(63)
	   end if
	   
	      if (spatialorder.eq.1) then 
	      firstorder = 1;iorder=1
	      else
	      firstorder = 0
	      end if

	     case(2)	!muscl type
	      
	   if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'----muscl scheme engaged----'
	      close(63)
	   end if
	  
	  iweno = -1
	  
	  
	   case(3)	!weno type
	      
	   if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'----weno scheme engaged----'
	      close(63)
	   end if
	  
	  iweno = 1
	  iweno = 1
	  iextend = 12	;
	  
	  
	  if (ees.eq.5)then
	   if (icompact.eq.1)then
	  iextend = 10
	  
	  
	  
	  else
	  iextend = 5
	  
	  end if
	  end if
	  wenwrt=wenocnschar
	  
	  if (dimensiona.eq.2)then
	    if ((ees.eq.0))then;typesten = 5;end if
	    if ((ees.eq.5))then;typesten = 5;end if
	    if (ees.eq.1)then;typesten = 9;end if
	    if (ees.eq.2)then;typesten = 9;end if
	    if (ees.eq.3)then;typesten = 5;end if
	  else
	    if ((ees.eq.0))then;typesten = 7;end if
	    if ((ees.ge.4))then;typesten = 7;end if
	    if (ees.eq.1)then ;typesten = 15;end if
	    if (ees.eq.2)then;typesten = 15;end if
	    if (ees.eq.3)then;typesten = 7;end if
	  end if
	  
	  
	  
	  end select


	  if (code_profile.eq.17)iweno=5
	  
        
	   
        ! temporal order
        rungekutta = temporder 
	  
       if ( iboundary .eq. 0 ) then 
	    iperiodicity = -3 
	    end if
	    if ( iboundary .eq. 1 ) then 
	    iperiodicity = 1 
	    end if

	    
	    
            if (dg.eq.1)then
                
                igqrules=min(iorder+1,6)
                else
                igqrules=min(iorder,6)
                end if
            

             
		alls=igqrules**(dimensiona)
             
        
	    
	    
	  !-------------------------end discretisation 6---------------------------------!




	call mpi_barrier(mpi_comm_world,ierror)

	    if (n.eq.0)then
	      open(63,file='history.txt',form='formatted',action='write',position='append')
	      write(63,*)'total number of processes:',isize,igqrules
	      close(63)
	  end if
        if(srfg.eq.1)then
            if (n.eq.0)then
                open(63,file='history.txt',form='formatted',action='write',position='append')
                write(63,*)'single reference frame engaged:'
                close(63)
            end if
        end if
	  		if(mrf.eq.1)then
            if (n.eq.0)then
                open(63,file='history.txt',form='formatted',action='write',position='append')
                write(63,*)'number of rotating frames engaged:',nrotors
                close(63)
            end if
        end if
		if(per_rot.eq.1)then
            iboundary=1
            lowmem=1
            if (n.eq.0)then
                open(63,file='history.txt',form='formatted',action='write',position='append')
                write(63,*)'rotational  periodicity engaged'
                close(63)
            end if
        end if
	  
	   
	  
	  if (initcond.eq.444)then
	inquire (file='BUBBLES.DAT',exist=here5)
	if (here5) then
	open(18,file='BUBBLES.DAT',form='formatted',status='old',action='read')
	read(18,*)nof_bubbles
	do ix=1,nof_bubbles
	read(18,*)bubble_centre(ix,1),bubble_centre(ix,2),bubble_radius(ix)
	end do
	close(18)
	end if
    end if
    
    if (initcond.eq.445)then

		inquire (file = 'BUBBLES.DAT',exist=here5)
		if (here5) then
		open(18,file='BUBBLES.DAT',form='formatted',status='old',action='read')
		read(18,*)nof_bubbles
		do ix=1,nof_bubbles
		read(18,*)bubble_centre(ix,1),bubble_centre(ix,2),bubble_centre(ix,3),bubble_radius(ix)
		end do
		close(18)

		end if
		end if
		

		inquire (file='BLEED.DAT',exist=bleedio)
	if (bleedio) then
	bleed=1
	open(27,file='BLEED.DAT',form='formatted',status='old',action='read')
	read(27,*)
	read(27,*)
	read(27,*)bleed_number !bleed_number=number of bleed boundary conditions
	read(27,*)bleed_type	!bleed_type=bleed type (slater=1,other=2,porous=3,etc)
	allocate(bleed_start(1:bleed_number,1:dims),bleed_end(1:bleed_number,1:dims),bleed_plenum(1:bleed_number), bleed_porosity(1:bleed_number))
	do ibleed=1,bleed_number
    read(27,*)bleed_start(ibleed,1:dims),bleed_end(ibleed,1:dims)
    end do
    do ibleed=1,bleed_number
    read(27,*)bleed_plenum(ibleed), bleed_porosity(ibleed)
    end do
    close(27)
	else
	bleed=0
	end if


	

	end subroutine read_ucns3d
	
	
	
end module parameters

