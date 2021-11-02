module parameters
USE MPIINFO
use DECLARATION

IMPLICIT NONE

 CONTAINS


SUBROUTINE READ_UCNS3D
!> @brief
!> This subroutine reads the parameter file

	IMPLICIT NONE

 	Integer :: INV
 	Real :: angledum
	CHARACTER(48)::STAMP1
	LOGICAL::HERE1,HERE2,HERE3


 	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	INQUIRE (FILE='RESTARTav.dat',EXIST=HERE1)
	IF (HERE1) THEN
	Average_restart=1
	Else
	Average_restart=0
	end if
	
	
	
	
	
	INQUIRE (FILE='MULTISPECIES.DAT',EXIST=HERE2)
	IF (HERE2) THEN
	MULTISPECIES=1
	OPEN(14,FILE='MULTISPECIES.DAT',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(14,*)
	READ(14,*)
	READ(14,*)NOF_SPECIES
    ALLOCATE(GAMMA_IN(1:NOF_SPECIES),MP_A_IN(1:NOF_SPECIES),MP_R_IN(1:NOF_SPECIES),MP_PINF(1:NOF_SPECIES))
    READ(14,*)GAMMA_IN(1:NOF_SPECIES)
    READ(14,*)MP_A_IN(1:NOF_SPECIES)
    READ(14,*)MP_R_IN(1:NOF_SPECIES)
    READ(14,*)MP_PINF(1:NOF_SPECIES)
    CLOSE(14)
	ELSE
	MULTISPECIES=0
	END IF


	INQUIRE (FILE='MOOD.DAT',EXIST=HERE3)
	IF (HERE3) THEN
	MOOD=1
	OPEN(17,FILE='MOOD.DAT',FORM='FORMATTED',STATUS='OLD',ACTION='READ')
	READ(17,*)
	READ(17,*)
	READ(17,*)MOOD_MODE        !TYPE OF MOOD MODE (1=RELAXED, 0=ORIGINAL)
    READ(17,*)MOOD_VAR1,MOOD_VAR2
    READ(17,*)MOOD_VAR3,MOOD_VAR4
    CLOSE(17)
	ELSE
	MOOD=0
	END IF
	
	
	
	
			
	
	
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	
	
	
	
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)DIMENSIONA,STATISTICS,CODE_PROFILE, RECONSTRUCT_HIGHER_ORDER_DG_DOFS_BOOLEAN
	READ(15,*)
	READ(15,*)
	READ(15,*)governingequations,INITCOND
	READ(15,*)
	READ(15,*)
	READ(15,*)turbulence,icoupleturb,PASSIVESCALAR
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)RRES,ufreestream,VVEL,WVEL,PRES
	READ(15,*)
	READ(15,*)
	READ(15,*)AOA,vectorx,vectory,vectorz
	READ(15,*)
	READ(15,*)
	READ(15,*)GAMMA,PRANDTL,Reynolds,CharLength
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)spatiladiscret,iRiemann,spatialorder,LIMITER,POLY
	READ(15,*)
	READ(15,*)
	READ(15,*)wenocnschar,EES,wenoz,wenocentralweight
	READ(15,*)
	READ(15,*)
	READ(15,*)temporder,CFL,timestep,upperlimit,reSLIMIT
	READ(15,*)
	READ(15,*)
	READ(15,*)iboundary,boundtype,SCALER
	READ(15,*)
	READ(15,*)
	READ(15,*)GREENGO,LMACH
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)OUT_TIME,NTMAX,WALLC
	READ(15,*)
	READ(15,*)
	READ(15,*)TECPLOT,IEVERY,IEVERY2,IEVERYAV,STENCIL_IO
	READ(15,*)
	READ(15,*)
	READ(15,*)Averaging
	READ(15,*)
	READ(15,*)
	READ(15,*)OUTSURF,IFORCE,surfshear
	READ(15,*)
	READ(15,*)
	READ(15,*)IRES_TURB,IRES_UNSTEADY,LAMPS,Prev_turbmodel  
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)
	READ(15,*)NPROBES
	READ(15,*)
	    
	
	DG=0
	!TURBULENCE DEFAULT VALUES
	
	
			! 	TURBULENCE MODEL PARAMETERS:
			! ---OPTIONS Spalart Almaras---
			ISPAL=1 !    ||SPALART ALLMARAS VERSION:| 1:ORIGINAL |2: NEGATIVE MODIFICATION
			!DES_model=0! 0              			|| 1-Detached Eddy Simulation 2-Delayed DES
			! ---CONSTANTS-------
			CB1=0.1355	! 0.1355				|| Cb1 
			CB2=0.622! 0.622   			|| Cb2
			SIGMA=0.666666667! 0.666666667      		|| sigma
			KAPPA=0.41  ! 0.41     			|| kappa
			CW1=3.23886781677! 3.23886781677      		|| CW1 (0 for (CB1/kappa**2)+1/sigma(1+CB2) or predifined)
			CW2=0.3! 0.3      			|| CW2
			CW3=2.0 ! 2.0      			|| CW3
			CV1=7.1! 7.1      			|| CV1
			CT1=1.0 ! 1.0      			|| CT1
			CT2=2.0! 2.0      			|| CT2
			CT3=1.1 ! 1.1      			|| CT3
			CT4=2.0 ! 2.0      			|| CT4
			PRTU=0.9! 0.9				|| Turbulent Prandtl Number
			TWALL=0! 0				|| Wall temperature (Kelvin) leave 0 for adiabatic (q_wall =0 <=> dT/dn=0)
			TURBINIT=3.0 ! 3.0	  			|| Initial value for turbulence parameter (multiplyied by the freestream viscosity from given Re)
			Upturblimit=1000000! 1000000				|| Upper limit for turbulence
			residualfreq=10! 10				|| Residual compute every
			IRS=0! 0				||IMPLICIT RESIDUAL SMOOTHING (DOUBLES CFL)
			C_DES_SA=0.61	! 0.61				||C_DES_SA
			! =============================================================================
			! K-OMEGA SST: 
			! ---OPTIONS---
			VORT_MODEL=0! 0					||0:Default strain-production for k, 1:Vorticity-production for k
			QSAS_MODEL=0! 0					||Hybrid models: 1-SAS-SST  //  2-DES-SST 
			ZERO_TURB_INIT=0! 0					||Zero turbulence option (if 1, initialization for k/w is done with zero turbulence)
			! --------------K-OMEGA CONSTANTS--------------------------
			SIGMA_K1=1.176470588! 1.176470588				||sigma_k1
			SIGMA_K2=1.0! 1.0					||sigma_k2
			SIGMA_OM1=2.0	! 2.0					||sigma_om1
			SIGMA_OM2=1.168! 1.168					||sigma_om2
			AA_1=0.31! 0.31					||aa_1
			BETA_I1=0.075	! 0.075					||beta_i1
			BETA_I2=0.0828! 0.0828					||beta_i2
			ALPHA_STARINF=1.0! 1.0					||alpha_starinf
			ALPHA_0=0.111111111! 0.111111111				||alpha_0
			BETA_STARINF=0.09! 0.09					||beta_starinf
			R_BETA=8.0! 8.0					||R_beta
			R_K_SST=6.0! 6.0					||R_k_SST
			BETA_T=0.072	! 0.072					||beta_t
			KAPPA_sst=0.41! 0.41					||kappa   (Same in Spalart-Allmaras)
			R_OM_SST=2.95	! 2.95					||R_om_SST
			ZETA_STAR= 1.5	! 1.5					||zeta_star  (Only for Mach corrections)
			M_T0=0.25! 0.25					||M_t0		(Only for Mach corrections)
			C_MU_INLET=0.09! 0.09					||C_mu_inlet
			C_SMG=0.11! 0.11					||C_smg  (SAS)
			ETA2_SAS=3.51! 3.51					||eta2_SAS (SAS)
			SIGMA_PHI=0.666666667	! 0.666666667				||sigma_phi (SAS)
			C_SAS=2.0! 2.0					||C_SAS   (SAS)
			L_turb_inlet=0.0026! 0.0026					||L_turb_inlet (Default turbulence lengthscale)
			I_turb_inlet=0.1! 0.1					||I_turb_inlet (Default turbulene intensity)
			Init_mu_ratio=0.01! 0.01					||Init_mu_ratio	(Default initial ratio for turbulent viscosity)
			C_DES_SST=0.61! 0.61					||C_DES_SST
			! =============================================================================
			! PASSIVE SCALAR TRANSPORT
			SCHMIDT_LAM=10.0! 10.0					||Laminar Schmidt number
			SCHMIDT_TURB=0.7! 0.7					||Turbulent Schmidt number
	
	
	    
	    
	    SELECT CASE(CODE_PROFILE)
	
	
	CASE (0)       
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
! 	icoupleturb=0	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
    if (initcond.eq.405)iadapt=1
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=2		!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 
	 
	 DES_model=0
	 
	 
	 
	 
	 CASE (9)  !robust
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
! 	icoupleturb=0	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
    if (initcond.eq.405)iadapt=1
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=3	!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 
	 
	 DES_model=0
	 
	 
	 CASE (91)  !robust WITH Matrix free LU-SGS
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
! 	icoupleturb=0	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
    if (initcond.eq.405)iadapt=1
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=3	!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=3		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 
	 
	 DES_model=0
	 
	 
	 
	 
	CASE (1)           !FOR DDES
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	!icoupleturb=1	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=3		!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 DES_model=2
	 
	 
	 
	 CASE (11)           !FOR DDES with matrix free LU-SGS
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	!icoupleturb=1	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	swirl=0		!swirling flow:0 deactivated, 1 activated
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=3		!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	!Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	!TOLERANCE=n_iterations
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	fastest=0	! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 DES_model=2
	 
	 
	 
    CASE (100)           !FOR DG method
	
	LOWMEMORY=0 	!MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	BINIO=1	    	!I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	!GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	REDUCE_COMP=0	!QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	TURBULENCEMODEL=1 !TURBULENCE MODEL SELECTION: |1:SPALART-ALLMARAS |2:K-W SST	
	!ICOUPLETURB=1	!COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	IHYBRID=0	!HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0 !HYBRID DISTANCE
	SWIRL=0		!SWIRLING FLOW:0 DEACTIVATED, 1 ACTIVATED
	IADAPT=0	!ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	!COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	EXTF=2		!STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	!WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	GUASSIANQUADRA=0!GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	!STORE GQP POINTS (1 =YES FASTER, 0= SLOWER)
    RELAX=1		!RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	!CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	!CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	EMETIS=6    	!METIS PARTITIONER : 1: HYBRID METIS, 2:ADAPTIVE WEIGHTS FOR HYBRID GRIDS, 3: UNIFORM METIS PARTIONIONER,4:NODAL,6=PARMETS 
	ITOLD=10000	!TOLERANCE=N_ITERATIONS
	GRIDAR1=5.0	! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	! LIMIT VOLUME CELLS
	FASTEST=0	! 0		       		||FASTEST, NO COORDINATE MAPPING (1: ENGAGED,0:WITH TRANSFORMATION)
	LMACH_STYLE=0	!0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),LMACH_STYLE(0=ONLY NORMAL COMPONENT,1=ALL COMPONENTS)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	!LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	DG=1    !0=DEACTIVATED FV ONLY, 1=ACTIVATED DG ONLY, 2=HYBRID
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	 DES_model=2
	 
	 
	 
	 
	
	
	CASE DEFAULT
	
	LOWMEMORY=0 	 !MEMORY USAGE: |0: HIGH(FASTER) |1:LOW (SLOWER)|| 
	binio=1	    	 !I/O (ASCII=0, BINARY=1) 
	LOWMEM=0    	 !GLOBAL ARRAYS SETTING (0=WITHOUT BETTER SUITED FOR NON PERIODIC BOUND,1=WITH (LARGE MEMORY FOOTPRINT))
	reduce_comp=0	 !QUADRATURE FREE FLUX=0 NOT TRUE,1 TRUE
	turbulencemodel=1 !TURBULENCE MODEL SELECTION: |1:Spalart-Allmaras |2:k-w SST	
	icoupleturb=0	 !COUPLING TURBULENCE MODEL: |1:COUPLED | 0: DECOUPLED
	ihybrid=0	 !HYBRID TURBULENCE : |1:ENABLED|0:DISABLED
	HYBRIDIST=0.0D0  !HYBRID DISTANCE
	swirl=0		 !swirling flow:0 deactivated, 1 activated
	IADAPT=0	 !ADAPTIVE NUMERICAL SCHEME (0 NOT TRUE,1 TRUE)
	ICOMPACT=0	 !COMPACT STENCIL MODE(0 NOT TRUE,1 TRUE)
	extf=2		 !STENCILS STABILITY VALUES FROM 1.2 TO 3 (DEFAULT 2)
	WEIGHT_LSQR=0	 !WEIGHTED LEAST SQUARES(0 NOT TRUE,1 TRUE)
	guassianquadra=0 !GAUSSIAN QUADRATURE RULE (1,2,5,6), DEFAULT 0 WILL USE THE APPROPRIATE NUMBER
	FASTEST_Q=1	 !STORE gqp POINTS (1 =YES FASTER, 0= SLOWER)
        relax=1		 !RELAXATION PARAMETER : |1:BLOCK JACOBI |2: LU-SGS
	CFLMAX=30	 !CFLMAX:TO BE USED WITH RAMPING
	CFLRAMP=0	 !CFL RAMPING: |0: DEACTIVATED |1:ACTIVATED
	emetis=6    	 !Metis partitioner : 1: Hybrid metis, 2:adaptive weights for hybrid grids, 3: Uniform metis partionioner,4:NODAL,6=PARMETS 
	itold=10000	 !TOLERANCE=n_iterations
	GRIDAR1=5.0	 ! 0	  5.0    7.0  LIMIT ASPECT RATIO CELLS,
	GRIDAR2=7.0	 ! LIMIT VOLUME CELLS
	fastest=0	 ! 0		       		||Fastest, no coordinate mapping (1: engaged,0:with transformation)
	lmach_style=0	 !0			||LOW MACH TREATMENT (1 ACTIVATE, 0 DISABLE),lmach_style(0=only normal component,1=all components)
	LAMX=1.0D0;LAMY=1.0D0;LAMZ=1.0D0	 !LINEAR ADVECTION COEFFICIENTS (LAMX, LAMY,LAMZ)
	ISPAL=1! 1				||SPALART ALLMARAS VERSION:| 1:ORIGINAL |2: NEGATIVE MODIFICATION
	
	
	
	if (iboundary.eq.1)then
	 LOWMEM=1
	 end if
	
	
	DES_model=0
	
	
	
	
	END SELECT
	    
	    
	
	 IF (N.EQ.0)THEN
      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
      CALL DATE_and_time(STAMP1)
      WRITE(63,*)"UCNS3D Parallel Computation Start",STAMP1
      END IF

      
      
	    NPROC=isize-1;vorder=iorder
      
	    !--------------------------1------------------------------!
	    !PROBE POSITION ALLOCATION AND READING OF COORDINATES
	    !$OMP MASTER
	    IF (NPROBES.GT.0)THEN
	    IF (DIMENSIONA.EQ.3)THEN
		  ALLOCATE(PROBEC(1:NPROBES,1:3))
		    DO INV=1,NPROBES
		      READ(15,*)PROBEC(INV,1),PROBEC(INV,2),PROBEC(INV,3)
		    END DO
	    ELSE
		  ALLOCATE(PROBEC(1:NPROBES,1:2))
	    
		    DO INV=1,NPROBES
		      READ(15,*)PROBEC(INV,1),PROBEC(INV,2)
		    END DO
	    END IF
	    END IF
	    !$OMP END MASTER
	     !--------------------------END 1-------------------------!
	     
	     
	     CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	     
	     
	     !-------------------------2---------------------------------!
	     !NUMBER OF VARIABLES SETUP
	     if (governingequations.le.2)then ! NS or Euler
		  if (dimensiona.eq.3)then
		      nof_variables=5;dims=3
		  else
		      nof_variables=4;dims=2
		  end if
		  
		  
	      else ! Linear advection
		  nof_variables=1
		  
		  if (dimensiona.eq.3)then
		    dims=3
		    else
		    dims=2
		    end if
	      end if
	      
	      !multiphase modification starts
	      IF (governingequations.EQ.-1)THEN
	      if (dimensiona.eq.3)then
		      nof_variables=5+NOF_SPECIES+(NOF_SPECIES-1);dims=3
		  else
		      nof_variables=4+NOF_SPECIES+(NOF_SPECIES-1);dims=2
		  end if
		  END IF
		  !multiphase modification ends
		  
	    !--------------------------END 2-------------------------!


	    !-------------------------3---------------------------------!
	     !TURBULENCE
	      if (turbulence .eq.1) then
		      if (turbulencemodel .eq. 1) then
			      turbulenceequations=1
		      else if (turbulencemodel .eq. 2) then
			      turbulenceequations=2
		      end if 
	      eLSE
		      turbulenceequations=0
		      turbulencemodel=0
	      END IF
	    !------------------------- END 3---------------------------------!
	   
	   
! 	   !-------------------------4---------------------------------!
	   !EQUATIONS TYPE
	   
	   SELECT CASE(governingequations)
	   CASE(-1)
	   !MULTI-PHASE INVISCID EULER EQUATIONS
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Multi-Phase inviscid Euler Solver Engaged'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 3;IVORTEX = 0
	   
	   CASE(1)
	   !NAVIER STOKES EQUATIONS
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Navier-Stokes Solver Engaged'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 4;IVORTEX = 1
	    
	  CASE(2)
	   !EULER EQUATIONS
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Euler Solver Engaged'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 3;IVORTEX = 0
	    
	    
	    
	  CASE(3)
	   !LINEAR ADVECTION EQUATION
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Linear Solver Engaged Sinewave'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 1;IVORTEX = 0
	    
	 CASE(4)
	   !LINEAR ADVECTION EQUATION
	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Linear Solver Engaged STEP FUNCTION'
	      CLOSE(63)
	    END IF
	    ITESTCASE = 2;IVORTEX = 0
	    
	  END SELECT
	  
	  
	  
	  !-------------------------END EQUATIONS 4---------------------------------!
	  
	  !-------------------------5---------------------------------!
	  !FLOW PARAMETERS
	  UNWOU = 3
	  BETAAS=1.5d0
	  SUTHER=0.412158681d0
	  uvel=ufreestream
	   ! Set pressure
	  if ( PRES .lt. 0 ) PRES = RRES/GAMMA	
	  ! Set dynamic free-stream viscosity
	  VISC = (RRES*ufreestream*CharLength)/Reynolds
	  if (swirl.eq.1)then
	  uvel=ZERO
	  end if
	  If (AOA .NE. 0.0D0) Then
	  angledum=(AOA*PI)/180.0d0
	  UVEL = COS(angledum)*ufreestream*VECTORX
	  WVEL = SIN(angledum)*ufreestream*VECTORZ
	  VVEL = SIN(angledum)*ufreestream*VECTORY
	  end if
	  
	  IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      if (ITESTCASE .eq. 4) Then
		write(63,*)'----Reynolds Number:',(RRES*ufreestream*CharLength)/VISC
		end if
	      CLOSE(63)
	   END IF
	  
	  
	  
	  !-------------------------END FLOW PARAMETERS 5---------------------------------!
	  
	  !-------------------------6---------------------------------!
	   !SPATIAL & TEMPORAL DISCRETISATION
	   
	   IORDER = max(1,spatialorder-1)
	   firstorder = 0
! 	   if (ees.eq.5)then
! 	   icompact=1
! 	   end if
	   IEXTEND =2; IF (ICOMPACT.EQ.1)IEXTEND = 8
	   IWENO =0
	   WENWRT = 1
	   ISCHEME=3
	   TYPESTEN=1
	   LWCI1=wenocentralweight
	  if (fastest.eq.1)ischeme=1
	    
	   
	   SELECT CASE(spatiladiscret)
	   
	   CASE(1)	!NO LIMITER
	      
	   IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'----1st Order in Space|LINEAR Scheme'
	      CLOSE(63)
	   END IF
	   
	      if (spatialorder.eq.1) then 
	      firstorder = 1;iorder=1
	      else
	      firstorder = 0
	      end if

	     CASE(2)	!MUSCL TYPE
	      
	   IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'----MUSCL Scheme Engaged----'
	      CLOSE(63)
	   END IF
	  
	  IWENO = -1
	  
	  
	   CASE(3)	!WENO TYPE
	      
	   IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'----WENO Scheme Engaged----'
	      CLOSE(63)
	   END IF
	  
	  IWENO = 1
	  IWENO = 1
	  IEXTEND = 12	;
	  
	  
	  if (EES.EQ.5)THEN
	   IF (ICOMPACT.EQ.1)THEN
	  IEXTEND = 8
	  
	  
	  
	  ELSE
	  IEXTEND = 5
	  
	  END IF
	  END IF
	  WENWRT=wenocnschar
	  
	  IF (DIMENSIONA.EQ.2)THEN
	    if ((EES.EQ.0).OR.(EES.GE.4))THEN;TYPESTEN = 5;End if
	    if (EES.EQ.1)THEN;TYPESTEN = 9;end if
	    if (EES.EQ.2)THEN;TYPESTEN = 9;End if
	    if (EES.EQ.3)THEN;TYPESTEN = 5;end if
	  ELSE
	    if ((EES.EQ.0).OR.(EES.GE.4))THEN;TYPESTEN = 7;End if
	    if (EES.EQ.1)THEN ;TYPESTEN = 15;End if
	    if (EES.EQ.2)THEN;TYPESTEN = 15;End if
	    if (EES.EQ.3)THEN;TYPESTEN = 7;End if
	  END IF
	  
	  
	  
	  END SELECT
	  
        IF (N.EQ.0)THEN
            OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
            write(63,*)'Order of Accuracy in space:',spatialorder
            CLOSE(63)
        END IF
	   
        ! Temporal order
        RUNGEKUTTA = temporder 
	  
        if ( iboundary .eq. 0 ) then 
            IPERIODICITY = -3 
        else if ( iboundary .eq. 1 ) then 
            IPERIODICITY = 1 
        end if

        IF (DG == 1) THEN
            IGQRULES = MIN(IORDER+1,6)
	    ELSE if (guassianquadra.eq.0) then
            IGQRULES=min(iorder,6)
        else 
            if (guassianquadra.gt.1)then
                IGQRULES = guassianquadra-1
            else
                IGQRULES =guassianquadra
            end if
        End if
	    
	    
	  !-------------------------END DISCRETISATION 6---------------------------------!




	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

	    IF (N.EQ.0)THEN
	      OPEN(63,FILE='history.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
	      write(63,*)'Total Number of Processes:',isize
	      CLOSE(63)
	  END IF

	  
	   
	  
	  

		
	

	END SUBROUTINE READ_UCNS3D
	
	
	
END MODULE PARAMETERS
