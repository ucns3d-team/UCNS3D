MODULE SOURCE
USE LIBRARY
USE TRANSFORM
USE LOCAL
USE RIEMANN
USE FLOW_OPERATIONS
USE DECLARATION
IMPLICIT NONE

contains
SUBROUTINE SOURCES_COMPUTATION(N)
!> @brief
!> Sources computation 
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE,ICONSIDERED
	REAL,DIMENSION(TURBULENCEEQUATIONS)::SOURCE_T

	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO
	DO I=1,KMAXE
		ICONSIDERED=I

		if (turbulence.eq.1)then
		CALL SOURCES(N,ICONSIDERED,SOURCE_T)
		RHST(I)%VAL(1:turbulenceequations)=RHST(I)%VAL(1:turbulenceequations)-(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)

		end if


	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_COMPUTATION


SUBROUTINE SOURCES_COMPUTATION_ROT(N)
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE,SRF
	REAL::OODENSITY
	real,dimension(5)::source_t2
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ

	
	
	KMAXE=XMPIELRANK(N)
	IF(SRFG.EQ.1)THEN
	!$OMP DO
	DO I=1,KMAXE
		
		 OODENSITY=1.0D0/U_C(I)%VAL(1,1)
        
                SOURCE_T2(1)=ZERO
                SOURCE_T2(2)=U_C(I)%VAL(1,2)*OODENSITY-UVEL
                SOURCE_T2(3)=U_C(I)%VAL(1,3)*OODENSITY-VVEL
                SOURCE_T2(4)=U_C(I)%VAL(1,4)*OODENSITY-WVEL
                SOURCE_T2(5)=ZERO
        
                POX(1:3)= SOURCE_T2(2:4)  
                POY(1:3)=SRF_VELOCITY(1:3)
                SOURCE_T2(2:4)=U_C(I)%VAL(1,1)*VECT_FUNCTION(POX,POY)
		RHS(I)%VAL(1:NOF_VARIABLES)=RHS(I)%VAL(1:NOF_VARIABLES)+(SOURCE_T2(1:NOF_VARIABLES)*ielem(n,I)%totvolume)
		
	END DO
	!$OMP END DO
	END IF
	IF(MRF.EQ.1)THEN
	!$OMP DO
	DO I=1,KMAXE
		SRF=ILOCAL_RECON3(I)%MRF
		IF (ILOCAL_RECON3(I)%MRF.EQ.1)THEN
            OODENSITY=1.0D0/U_C(I)%VAL(1,1)
                SOURCE_T2(1)=ZERO
                SOURCE_T2(2)=U_C(I)%VAL(1,2)*OODENSITY!-UVEL
                SOURCE_T2(3)=U_C(I)%VAL(1,3)*OODENSITY!-VVEL
                SOURCE_T2(4)=U_C(I)%VAL(1,4)*OODENSITY!-WVEL
                SOURCE_T2(5)=ZERO
        
                POX(1:3)= SOURCE_T2(2:4)  
                POY(1:3)=ILOCAL_RECON3(I)%MRF_VELOCITY(1:3)
                SOURCE_T2(2:4)=U_C(I)%VAL(1,1)*VECT_FUNCTION(POX,POY)
		RHS(I)%VAL(1:NOF_VARIABLES)=RHS(I)%VAL(1:NOF_VARIABLES)+(SOURCE_T2(1:NOF_VARIABLES)*ielem(n,I)%totvolume)
		END IF
	END DO
	!$OMP END DO
	END IF
END SUBROUTINE SOURCES_COMPUTATION_ROT





SUBROUTINE SOURCES_realgas(N,ICONSIDERED)
  !> Sources computation in 2D/3D for 5-species, 2-T real gas
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: N, ICONSIDERED
  real,dimension(1:nof_Variables)::source_r, SRC_JAC_DIAG

  INTEGER :: I, J, K, L, rg_i, RG_J
  REAL,DIMENSION(nof_species) :: rg_R, RG_RM, RG_N      ! RG_R = ρ_s
  REAL :: RG_DQ_RAD, RG_QTV, RG_QW, velx, vely, velz
  REAL,DIMENSION(nof_species) :: RG_EV, RG_ESV
  REAL,DIMENSION(3) :: RG_TVS
  REAL,DIMENSION(nof_species) :: RG_DW_S
  REAL :: RG_KE1, RG_KE2, RG_KE3, RG_KE4, RG_KE5
  REAL :: RG_T, RG_TA, RG_TV
  REAL :: RG_Z, rg_speed
  REAL :: RG_R1, RG_R2, RG_R3, RG_R4, RG_R5
  REAL :: RG_KF1, RG_KF2, RG_KF3, RG_KF4, RG_KF5
  REAL :: RG_KB1, RG_KB2, RG_KB3, RG_KB4, RG_KB5

  REAL :: RG_KF11,RG_KF12,RG_KF13,RG_KF14,RG_KF15
  REAL :: RG_KF21,RG_KF22,RG_KF23,RG_KF24,RG_KF25
  REAL :: RG_KF31,RG_KF32,RG_KF33,RG_KF34,RG_KF35
  REAL :: RG_KB11,RG_KB12,RG_KB13,RG_KB14,RG_KB15
  REAL :: RG_KB21,RG_KB22,RG_KB23,RG_KB24,RG_KB25
  REAL :: RG_KB31,RG_KB32,RG_KB33,RG_KB34,RG_KB35

  REAL :: Kf11m,Kf12m,Kf13m,Kf14m,Kf15m
  REAL :: Kf21m,Kf22m,Kf23m,Kf24m,Kf25m
  REAL :: Kf31m,Kf32m,Kf33m,Kf34m,Kf35m
  REAL :: Kb11m,Kb12m,Kb13m,Kb14m,Kb15m
  REAL :: Kb21m,Kb22m,Kb23m,Kb24m,Kb25m
  REAL :: Kb31m,Kb32m,Kb33m,Kb34m,Kb35m
  REAL :: Kf4m,Kf5m,Kb4m,Kb5m
  REAL :: Kf1_eff,Kb1_eff,Kf2_eff,Kb2_eff,Kf3_eff,Kb3_eff
  REAL :: RG_C(5)            ! molar concentrations [mol/m^3]
  REAL :: wdot(5)            ! reaction rates [mol/(m^3 s)]
  INTEGER :: nu(5,5)         ! stoichiometric matrix ν(s,r)
  INTEGER :: s, r

  REAL,DIMENSION(nof_species,nof_species) :: RG_MU
  REAL :: TEMP1, TEMP2,tauMW,tauPk
  REAL,DIMENSION(nof_species) :: RG_TSMW, RG_TSP, RG_TSSUM, RG_Ablot, RG_EV_EQ
  REAL,DIMENSION(1:NOF_VARIABLES) :: LEFTV, tempvect
  REAL :: MP_PINFl, GAMMAL,Q_chem
  REAL :: P_atm, t_13, rg_pressure, n_tot, exp_argTV, exp_argT, Rs
  INTEGER :: rg_molx
  REAL,DIMENSION(1:nof_species) :: rg_mass_amu

  REAL, PARAMETER :: kB = 1.380649D-23     ! Boltzmann [J/K]
  REAL, PARAMETER :: NA = 6.02214076D23    ! Avogadro [1/mol]

  REAL,DIMENSION(1:nof_species) :: tau_MW    !! Millikan-White [s]
  REAL,DIMENSION(1:nof_species) :: tau_Park  ! Park [s]
  REAL,DIMENSION(1:nof_species) :: tau_tot   ! combined [s]
  REAL :: num, den,rggtF
  REAL :: mu_sr, a_sr, b_sr
  REAL :: log_p_tau, p_tau, tau_sr
  REAL :: sigma0, sigma_v, v_th, m_s,maxLoss
  rEAL :: rg_alpha, rg_alpha_i, rg_rho_min,sum_sourcex
  INTEGER :: idxE, idxEv, idxY1
 REAL    :: rho_min, lambda_s, lambda_E
 REAL    :: tau_v_mix, rhoE_local, Yi, sumY,q_TTv,Tv_mix,T_eff
  REAL :: dt_loc
   REAL :: rhoE_old, rhoE_new
   REAL :: rhoEv_old, rhoEv_new, rhoEv_eq
   REAL :: alpha_v
   REAL :: prod, dest, lambda_k, rhoY_old, rhoY_new
   REAL :: tau_v_eff


  ! ---------------------------------------------------------------------------
  ! Start
  ! ---------------------------------------------------------------------------
  I = ICONSIDERED

  SOURCE_R(1:NOF_VARIABLES) = 0.0D0
  tau_MW(:)   = 1.0D30
  tau_Park(:) = 1.0D30
  tau_tot(:)  = 1.0D30

  !-----------------------------------------------------------
  ! Get conservative variables and primitive variables
  !-----------------------------------------------------------
  LEFTV(1:nof_variables) = U_C(I)%VAL(1,1:nof_variables)

  ! Species densities ρ_s = ρ Y_s
  DO RG_I = 1, nof_species
     RG_R(RG_I) = LEFTV(dimensiona+3+RG_I)
  END DO

  !-----------------------------------------------------------
  ! 1) Molar concentrations C_s = rho_s / M_s
  !-----------------------------------------------------------
  DO s = 1, 5
     RG_C(s) = RG_R(s) / RG_MOLM(s)
  END DO

  CALL CONS2PRIM(N,LEFTV,MP_PINFl,GAMMAL)

  rg_pressure = LEFTV(dimensiona+2)
  P_atm       = rg_pressure / RGS_Pa_per_atm

  ! molar masses in kg/mol -> "amu-like" for MW reduced-mass formula
  rg_mass_amu = RG_MOLM * 1000.0D0

  tempvect(1:nof_variables) = U_C(I)%VAL(1,1:nof_variables)
  CALL CONS2div(N,tempvect,MP_PINFl,GAMMAL)

  ! LEFTV:   (ρ, u, v, ρE, ρEv, ρY1..ρY5)
  ! tempvect:(ρ, u, v, Ttr, Tv, Y1..Y5)

  ! Temperatures
  RG_T  = tempvect(dimensiona+2)   ! Ttr
  RG_TV = tempvect(dimensiona+3)   ! Tv
  RG_Z  = 1.0D0 / RG_T
  t_13  = RG_T**(-1.0D0/3.0D0)

  ! Adjusted two-temperature model (Candler/Park style)
  RG_TA = RG_T**0.6D0 * RG_TV**0.4D0

  ! Flow speed (not used in source terms yet)
  velx = LEFTV(2)
  vely = LEFTV(3)
  velz = 0.0D0
  IF (dimensiona.EQ.3) THEN
     velz = LEFTV(4)
  END IF
  rg_speed = SQRT(velx**2 + vely**2 + velz**2)

  !-----------------------------------------------------------
  ! Chemical source: reaction rates (Candler 5-reaction model)
  ! Using ρ-based rate coefficients Cf [m^3/(kg s)]
  !-----------------------------------------------------------
  RG_DW_S(:) = 0.0D0

  IF (RG_NOF_REACTIONS == 5) THEN

!
!
!
!
!
!
!      ! Equilibrium constants K_E(T) (Candler curve fits, T = Ttr)
      RG_KE1 = EXP(3.898D0 -12.611D0*RG_Z +0.683D0*RG_Z**2 -0.118D0*RG_Z**3 +0.006D0*RG_Z**4)
      RG_KE2 = EXP(1.335D0 - 4.127D0*RG_Z -0.616D0*RG_Z**2 +0.093D0*RG_Z**3 -0.005D0*RG_Z**4)
      RG_KE3 = EXP(1.549D0 - 7.784D0*RG_Z +0.228D0*RG_Z**2 -0.043D0*RG_Z**3 +0.002D0*RG_Z**4)
      RG_KE4 = EXP(2.349D0 - 4.828D0*RG_Z +0.455D0*RG_Z**2 -0.075D0*RG_Z**3 +0.004D0*RG_Z**4)
      RG_KE5 = EXP(0.215D0 - 3.652D0*RG_Z +0.843D0*RG_Z**2 -0.136D0*RG_Z**3 +0.007D0*RG_Z**4)
!




! --- Reaction 1: N2 + M <-> 2N + M (in m^3/mol/s) ---
  Kf11m = 3.78D18 * 1.0D-6 * RG_TA**(-1.6D0) * EXP(-1.132D5/RG_TA)  ! M = N2
  Kf12m = 3.78D18 * 1.0D-6 * RG_TA**(-1.6D0) * EXP(-1.132D5/RG_TA)  ! M = O2
  Kf13m = 3.78D18 * 1.0D-6 * RG_TA**(-1.6D0) * EXP(-1.132D5/RG_TA)  ! M = NO
  Kf14m = 1.11D18 * 1.0D-6 * RG_TA**(-1.6D0) * EXP(-1.132D5/RG_TA)  ! M = N
  Kf15m = 1.11D18 * 1.0D-6 * RG_TA**(-1.6D0) * EXP(-1.132D5/RG_TA)  ! M = O

  Kb11m = Kf11m / RG_KE1
  Kb12m = Kf12m / RG_KE1
  Kb13m = Kf13m / RG_KE1
  Kb14m = Kf14m / RG_KE1
  Kb15m = Kf15m / RG_KE1

  ! mixture-averaged 3rd-body rate for reaction 1 [1/s]
  Kf1_eff = Kf11m*RG_C(1) + Kf12m*RG_C(2) + Kf13m*RG_C(3) &
          + Kf14m*RG_C(4) + Kf15m*RG_C(5)
  Kb1_eff = Kb11m*RG_C(1) + Kb12m*RG_C(2) + Kb13m*RG_C(3) &
          + Kb14m*RG_C(4) + Kb15m*RG_C(5)

  ! molar rate for reaction 1: [mol/(m^3 s)]
  wdot(1) = Kf1_eff*RG_C(1) - Kb1_eff*RG_C(4)**2   ! N2 + M <-> 2N + M


  ! --- Reaction 2: O2 + M <-> 2O + M ---
  Kf21m = 2.75D16 * 1.0D-6 * RG_TA**(-1.0D0) * EXP(-5.95D4/RG_TA)   ! M = N2
  Kf22m = 2.75D16 * 1.0D-6 * RG_TA**(-1.0D0) * EXP(-5.95D4/RG_TA)   ! M = O2
  Kf23m = 2.75D16 * 1.0D-6 * RG_TA**(-1.0D0) * EXP(-5.95D4/RG_TA)   ! M = NO
  Kf24m = 8.25D16 * 1.0D-6 * RG_TA**(-1.0D0) * EXP(-5.95D4/RG_TA)   ! M = N
  Kf25m = 8.25D16 * 1.0D-6 * RG_TA**(-1.0D0) * EXP(-5.95D4/RG_TA)   ! M = O

  Kb21m = Kf21m / RG_KE2
  Kb22m = Kf22m / RG_KE2
  Kb23m = Kf23m / RG_KE2
  Kb24m = Kf24m / RG_KE2
  Kb25m = Kf25m / RG_KE2

  Kf2_eff = Kf21m*RG_C(1) + Kf22m*RG_C(2) + Kf23m*RG_C(3) &
          + Kf24m*RG_C(4) + Kf25m*RG_C(5)
  Kb2_eff = Kb21m*RG_C(1) + Kb22m*RG_C(2) + Kb23m*RG_C(3) &
          + Kb24m*RG_C(4) + Kb25m*RG_C(5)

  wdot(2) = Kf2_eff*RG_C(2) - Kb2_eff*RG_C(5)**2   ! O2 + M <-> 2O + M


  ! --- Reaction 3: NO + M <-> N + O + M ---
  Kf31m = 2.30D14 * 1.0D-6 * RG_TA**(-0.5D0) * EXP(-7.55D4/RG_TA)   ! M = N2
  Kf32m = 2.30D14 * 1.0D-6 * RG_TA**(-0.5D0) * EXP(-7.55D4/RG_TA)   ! M = O2
  Kf33m = 2.30D14 * 1.0D-6 * RG_TA**(-0.5D0) * EXP(-7.55D4/RG_TA)   ! M = NO
  Kf34m = 4.60D14 * 1.0D-6 * RG_TA**(-0.5D0) * EXP(-7.55D4/RG_TA)   ! M = N
  Kf35m = 4.60D14 * 1.0D-6 * RG_TA**(-0.5D0) * EXP(-7.55D4/RG_TA)   ! M = O

  Kb31m = Kf31m / RG_KE3
  Kb32m = Kf32m / RG_KE3
  Kb33m = Kf33m / RG_KE3
  Kb34m = Kf34m / RG_KE3
  Kb35m = Kf35m / RG_KE3

  Kf3_eff = Kf31m*RG_C(1) + Kf32m*RG_C(2) + Kf33m*RG_C(3) &
          + Kf34m*RG_C(4) + Kf35m*RG_C(5)
  Kb3_eff = Kb31m*RG_C(1) + Kb32m*RG_C(2) + Kb33m*RG_C(3) &
          + Kb34m*RG_C(4) + Kb35m*RG_C(5)

  wdot(3) = Kf3_eff*RG_C(3) - Kb3_eff*RG_C(4)*RG_C(5) ! NO + M <-> N + O + M


  ! --- Reaction 4: N2 + O <-> NO + N (bimolecular, 2-body) ---
  ! You originally had RG_KF4 in m^3/(kg s); we reinterpret it here
  ! as cm^3/(mol s) with the same A,n,Ea, convert to m^3/(mol s).
  Kf4m = 3.18D10 * 1.0D-6 * RG_TA**(0.1D0) * EXP(-3.77D4/RG_TA)
  Kb4m = Kf4m / RG_KE4

  wdot(4) = Kf4m*RG_C(1)*RG_C(5) - Kb4m*RG_C(3)*RG_C(4)

  ! --- Reaction 5: NO + O <-> O2 + N (bimolecular) ---
  Kf5m = 2.16D5 * 1.0D-6 * RG_TA**(1.29D0) * EXP(-1.922D4/RG_TA)
  Kb5m = Kf5m / RG_KE5

  wdot(5) = Kf5m*RG_C(3)*RG_C(5) - Kb5m*RG_C(2)*RG_C(4)


  ! Optionally keep RG_R1..RG_R5 as "mass rates" based on wdot,
  ! if you want to inspect them (not used elsewhere in source terms):
  RG_R1 = wdot(1) * RG_MOLM(1)   ! ~mass rate scale for reaction 1
  RG_R2 = wdot(2) * RG_MOLM(2)
  RG_R3 = wdot(3) * RG_MOLM(3)
  RG_R4 = wdot(4) * RG_MOLM(1)   ! arbitrary scaling; for diagnostics only
  RG_R5 = wdot(5) * RG_MOLM(3)

  !-----------------------------------------------------------
  ! 3) Stoichiometric update: RG_DW_S(s) [kg/(m^3 s)]
  !    ν(s,r) mol-based coefficients:
  !      r1: N2 + M <-> 2N + M
  !      r2: O2 + M <-> 2O + M
  !      r3: NO + M <-> N + O + M
  !      r4: N2 + O <-> NO + N
  !      r5: NO + O <-> O2 + N
  !-----------------------------------------------------------
  DO s = 1,5
     DO r = 1,5
        nu(s,r) = 0
     END DO
  END DO

  ! r1: N2 -> 2N
  nu(1,1) = -1
  nu(4,1) = +2

  ! r2: O2 -> 2O
  nu(2,2) = -1
  nu(5,2) = +2

  ! r3: NO -> N + O
  nu(3,3) = -1
  nu(4,3) = +1
  nu(5,3) = +1

  ! r4: N2 + O -> NO + N
  nu(1,4) = -1
  nu(5,4) = -1
  nu(3,4) = +1
  nu(4,4) = +1

  ! r5: NO + O -> O2 + N
  nu(3,5) = -1
  nu(5,5) = -1
  nu(2,5) = +1
  nu(4,5) = +1

  ! Build mass production rates
  RG_DW_S(:) = 0.0D0
  DO r = 1,5
     DO s = 1,5
        RG_DW_S(s) = RG_DW_S(s) + nu(s,r) * RG_MOLM(s) * wdot(r)
     END DO
  END DO


!   sum_sourcex=zero
! 	DO rg_i = 1, nof_species
! 		sum_sourcex=sum_sourcex+ RG_DW_S(rg_i)
! 	END DO





  rg_rho_min = 1.0D-20    ! or whatever you use as "zero" density
  rg_alpha   = 1.0D0

	! First: compute most restrictive scaling factor alpha
	DO rg_i = 1, nof_species

		! If density is essentially zero, don't allow further destruction
		IF (RG_R(rg_i) <= rg_rho_min .AND. RG_DW_S(rg_i) < 0.0D0) THEN
			RG_DW_S(rg_i) = 0.0D0
			CYCLE
		END IF

		! Only destruction can cause negativity
! 		IF (RG_DW_S(rg_i) < 0.0D0) THEN
! 			! Need: rho_i + dt * alpha * DW_i >= 0
! 			!  => alpha <= rho_i / (-dt * DW_i)
! 			rg_alpha_i = RG_R(rg_i) / (-ielem(n,ICONSIDERED)%dtl * RG_DW_S(rg_i))
!
! 			IF (rg_alpha_i < 1.0d-2) rg_alpha = 1.0d-2
! 		END IF
	END DO

	! Clamp alpha to [0,1]
	IF (rg_alpha > 1.0D0) rg_alpha = 1.0D0
	IF (rg_alpha < 0.0D0) rg_alpha = 0.0D0

	! Second: apply scaling if needed
	IF (rg_alpha < 1.0D0) THEN
		DO rg_i = 1, nof_species
			RG_DW_S(rg_i) = rg_alpha * RG_DW_S(rg_i)
		END DO
	END IF

	! Optional tiny cutoff (for cleanliness, won’t break stoichiometry)
	DO rg_i = 1, nof_species
		IF (ABS(RG_DW_S(rg_i)) < 1.0D-50) RG_DW_S(rg_i) = 0.0D0
	END DO


! 	sum_sourcex=zero
! 	DO rg_i = 1, nof_species
! 		sum_sourcex=sum_sourcex+ RG_DW_S(rg_i)
! 	END DO
!
! 	write(110+n,*)sum_sourcex



  ELSE
     RG_DW_S(:) = 0.0D0
  END IF

  !-----------------------------------------------------------
  ! 1) Millikan-White VT relaxation times (species 1..3)
  !-----------------------------------------------------------
  DO RG_I = 1, 3
     num = 0.0D0
     den = 0.0D0

     DO RG_J = 1, nof_species
        IF (RG_R(RG_J) <= 0.0D0) CYCLE

        ! Reduced mass μ_sr in amu
        mu_sr = rg_mass_amu(RG_I)*rg_mass_amu(RG_J) / (rg_mass_amu(RG_I)+rg_mass_amu(RG_J))

        ! MW coefficients
        a_sr = 0.00116D0 * SQRT(mu_sr) * rg_thetag(RG_I)**(4.0D0/3.0D0)
        b_sr = 0.015D0   * mu_sr**0.25D0

        ! p * tau_sr in atm*s
        log_p_tau = a_sr*(t_13 - b_sr) - 18.42D0
        p_tau       = exp(log_p_tau)

        ! tau_sr [s]: p_tau has atm*s, divide by p in atm
        tau_sr = p_tau / P_atm

        ! mixture rule with mole fractions:
        ! 1/tau_s = Σ_r (x_r / tau_sr)
        ! x_r ~ n_r / Σ n_r, and n_r ~ ρ_r / M_r
        num = num + (RG_R(RG_J)/RG_MOLM(RG_J)) / tau_sr
        den = den + (RG_R(RG_J)/RG_MOLM(RG_J))
     END DO

     IF (num > 0.0D0 .AND. den > 0.0D0) THEN
        tau_MW(RG_I) = den / num
     ELSE
        tau_MW(RG_I) = 1.0D30
     END IF



  END DO

!   !-----------------------------------------------------------
!   ! 2) Park (1989) high-temperature correction (number density)
!   !    tau_Park_s = 1 / ( n_tot * sigma_v * v_th )
!   !    n_tot = p / (kB T)
!   !-----------------------------------------------------------
!   sigma0 = 3.0D-21    ! [m^2]
!
!   ! Number density [1/m^3]
!   n_tot = rg_pressure / (kB * RG_T)
!
!   DO RG_I = 1, 3
!      ! particle mass [kg/particle]
!      m_s = RG_MOLM(RG_I) / NA
!
!      ! effective cross-section (kept as in your original implementation)
!      sigma_v = sigma0 * (50000/RG_T)**2
!
!      ! thermal speed
!      v_th = SQRT( 8.0D0 * kB * RG_T / (pi * m_s) )
!
!      ! tau_Park [s]
!      tau_Park(RG_I) = 1.0D0 / ( n_tot * sigma_v * v_th )
!   END DO


  sigma0 = 3.0D-21    ! [m^2]




! Number density [1/m^3] based on translational/rotational T
n_tot = rg_pressure / (kB * RG_T)   ! RG_T = T_tr

! Single vibrational temperature for the mixture
         ! <-- your global vibrational temperature

! Avoid non-physical Tv
Tv_mix = MAX(RG_TV, 300.0D0)

! Park-style weighting between T and Tv (0.7 is a common choice)
q_TTv = 0.7D0

! Effective two-temperature average (same for all species)
! geometric mean: T_eff = T^q * Tv^(1-q)
T_eff = RG_T**q_TTv * Tv_mix**(1.0D0 - q_TTv)
T_eff = MAX(T_eff, 300.0D0)

DO RG_I = 1, 3
   ! particle mass [kg/particle]
   m_s = RG_MOLM(RG_I) / NA

   ! effective cross-section using T_eff instead of RG_T
   sigma_v = sigma0 * (50000.0D0 / T_eff)**2

   ! thermal speed using T_eff
   v_th = SQRT( 8.0D0 * kB * T_eff / (pi * m_s) )

   ! tau_Park [s]
   tau_Park(RG_I) = 1.0D0 / ( n_tot * sigma_v * v_th )

END DO
























  !-----------------------------------------------------------
  ! 3) Combine Millikan-White + Park
  !    1/tau_tot = 1/tau_MW + 1/tau_Park
  !-----------------------------------------------------------


 DO RG_I = 1, 3

  ! Local copies
  tauMW  = tau_MW(RG_I)
  tauPk  = tau_Park(RG_I)

  ! Replace non-positive with 'disabled'
  IF (tauMW <= 0.0D0) tauMW = 1.0D30
  IF (tauPk <= 0.0D0) tauPk = 1.0D30

  ! Now combine them safely
  IF (tauMW >= 1.0D29 .AND. tauPk >= 1.0D29) THEN
     ! both basically 'off'
     tau_tot(RG_I) = 1.0D30

  ELSEIF (tauMW >= 1.0D29) THEN
     ! only Park valid
     tau_tot(RG_I) = max(taupk,5.0e-9)

  ELSEIF (tauPk >= 1.0D29) THEN
     ! only MW valid
     tau_tot(RG_I) = tauMW

  ELSE
     ! both valid: use harmonic sum
     rggtF = tauMW / tauPk
      rggtF = MAX(1.0D0, rggtF)
      rggtF = MIN(20.0D0, rggtF)
      tau_tot(RG_I) = tauMW / rggtF
!      tau_tot(RG_I) = 1.0D0 / ( 1.0D0/tauMW + 1.0D0/tauPk )
!      ! or: tau_tot(RG_I) = MIN(tauMW, tauPk)
  END IF


  RG_TSSUM(RG_I) = tau_tot(RG_I)

END Do

  !-----------------------------------------------------------
  ! 4) Vibrational energies e_v(Tv) and e_v_eq(T)
  !-----------------------------------------------------------
  DO RG_I = 1, 3
     Rs         = RGS_Ru / RG_MOLM(RG_I)
     exp_argTV  = rg_thetag(RG_I) / RG_TV
     RG_EV(RG_I) = Rs * rg_thetag(RG_I) / (EXP(exp_argTV) - 1.0D0)

     exp_argT      = rg_thetag(RG_I) / RG_T
     RG_EV_EQ(RG_I)= Rs * rg_thetag(RG_I) / (EXP(exp_argT) - 1.0D0)
  END DO

  !-----------------------------------------------------------
  ! 5) VT relaxation source: QTV
  !    QTV = sum_i rho_i (e_v_eq - e_v) / tau_i
  !-----------------------------------------------------------
  RG_QTV = 0.0D0
  DO RG_I = 1, 3
     IF (RG_TSSUM(RG_I) .GT. 0.0D0) THEN
        RG_QTV = RG_QTV + RG_R(RG_I) * (RG_EV_EQ(RG_I) - RG_EV(RG_I)) / RG_TSSUM(RG_I)
     END IF
  END DO

  !-----------------------------------------------------------
  ! 6) Reactive vibrational source: QW
  !    QW = sum_i (dW_i * e_v(Tv))  (only vibrating species)
  !-----------------------------------------------------------
  RG_QW = 0.0D0
   DO RG_I = 1, 3
      RG_QW = RG_QW + RG_DW_S(RG_I) * RG_EV(RG_I)
   END DO

  !-----------------------------------------------------------
  ! 7) Assemble source terms for conservative variables
  !
  !   ρE  = total energy (includes vibrational energy)
  !   ρEv = vibrational energy only
  !
  !   -> VT exchange QTV is internal transfer between ρE and ρEv,
  !      so it should NOT be removed from ρE (to conserve total energy).
  !   -> Chemistry contributes Q_chem via heats of formation.
  !-----------------------------------------------------------
  SOURCE_R(1:nof_variables) = 0.0D0

  ! Chemical heat release (or absorption)
  Q_chem = 0.0D0
  DO rg_i = 1, nof_species
     ! rg_hzero in [J/mol] -> divide by molar mass [kg/mol] -> [J/kg]
     Q_chem = Q_chem + RG_DW_S(rg_i) * ( -rg_hzero(rg_i) / RG_MOLM(rg_i) )
  END DO


  ! Total energy equation: chemistry only
  SOURCE_R(dimensiona+2) = Q_chem

  ! Vibrational energy equation: VT + reactive vibrational source
  SOURCE_R(dimensiona+3) = RG_QTV + RG_QW

  ! Species equations: mass production rates
  SOURCE_R(dimensiona+4:nof_variables) = RG_DW_S(1:nof_species)

	SOURCE_R(:)=ZERO



	! ------------------------------------------------------------
  ! Diagonal source Jacobian for implicit time stepping
  ! ------------------------------------------------------------
  SRC_JAC_DIAG(1:NOF_VARIABLES) = 0.0D0

  idxE  = dimensiona + 2          ! ρE
  idxEv = dimensiona + 3          ! ρEv
  idxY1 = dimensiona + 4          ! first species equation

  rho_min = 1.0D-20

  ! --- Species diagonal terms: dS_s / d(ρ_s) ≈ -(-ω_s / ρ_s) = ω_s / ρ_s ---
  DO rg_i = 1, nof_species
     IF (RG_R(rg_i) > rho_min .AND. RG_DW_S(rg_i) < 0.0D0) THEN
        ! local destruction timescale τ_s = -ρ_s / ω_s  (ω_s < 0)
        lambda_s = -RG_DW_S(rg_i) / RG_R(rg_i)   ! [1/s]  > 0
        SRC_JAC_DIAG(idxY1 + rg_i - 1) = -lambda_s
     ELSE
        SRC_JAC_DIAG(idxY1 + rg_i - 1) = 0.0D0
     END IF
  END DO

  ! --- Vibrational energy diagonal: dS_Ev / d(ρEv) ≈ -1/τ_v_mix ---
  ! Use a simple mixture VT timescale based on tau_tot (or tau_MW) for N2,O2,NO
  tau_v_mix = 1.0D30
  sumY      = 0.0D0

  DO rg_i = 1, 3   ! vibrating species: N2, O2, NO
     IF (RG_R(rg_i) > rho_min .AND. tau_tot(rg_i) < 1.0D29) THEN
        Yi       = RG_R(rg_i) / LEFTV(1)          ! Y_i ≈ ρ_i/ρ
        sumY     = sumY + Yi
        tau_v_mix = MIN(tau_v_mix, tau_tot(rg_i)) ! use most restrictive VT time
     END IF
  END DO

  IF (tau_v_mix < 1.0D29) THEN
     SRC_JAC_DIAG(idxEv) = -1.0D0 / tau_v_mix
  ELSE
     SRC_JAC_DIAG(idxEv) = 0.0D0
  END IF

  ! --- (Optional) Total energy diagonal: dS_E / d(ρE) ≈ -1/τ_Echem ---
  rhoE_local = U_C(I)%VAL(1,idxE)
  IF (rhoE_local > 1.0D-12 .AND. Q_chem < 0.0D0) THEN
     lambda_E = -Q_chem / rhoE_local   ! [1/s]
     SRC_JAC_DIAG(idxE) = -lambda_E
  ELSE
     SRC_JAC_DIAG(idxE) = 0.0D0
  END IF

  ! If chemistry/VT are turned off, also zero the Jacobian
   SRC_JAC_DIAG(:) = 0.0D0


  !point implicit treatment of the


   dt_loc = IELEM(N,I)%DTL

   ! ----------------------------
   ! (A) Total energy: chemistry only (EXPLICIT, same as before)
   ! ----------------------------
   ! Original pathway:
   !   RHS_E -= Q_chem * V
   !   rhoE  -= dt * (RHS_E / V)
   ! ==> rhoE += dt * Q_chem
   !
   rhoE_old = U_C(I)%VAL(1, dimensiona+2)
   rhoE_new = rhoE_old + dt_loc * Q_chem
   U_C(I)%VAL(1, dimensiona+2) = rhoE_new

   ! ----------------------------
   ! (B) Vibrational energy: VT relaxation + reactive vib source
   ! ----------------------------
   ! Use BE-stable form:
   !   d(rhoEv)/dt = (rhoEv_eq(T) - rhoEv)/tau_v_eff + QW
   ! Relaxation part implicit, QW explicit.
   !
   rhoEv_old = U_C(I)%VAL(1, dimensiona+3)

   rhoEv_eq = 0.0D0
   DO rg_i = 1, 3
      rhoEv_eq = rhoEv_eq + RG_R(rg_i) * RG_EV_EQ(rg_i)
   END DO

   tau_v_eff = 1.0D30
   DO rg_i = 1, 3
      IF (RG_R(rg_i) > rho_min .AND. tau_tot(rg_i) < 1.0D29) THEN
         tau_v_eff = MIN(tau_v_eff, tau_tot(rg_i))
      END IF
   END DO

   IF (tau_v_eff < 1.0D29) THEN
      alpha_v = dt_loc / tau_v_eff
      rhoEv_new = (rhoEv_old + alpha_v * rhoEv_eq + dt_loc * RG_QW) / (1.0D0 + alpha_v)
   ELSE
      rhoEv_new = rhoEv_old + dt_loc * (RG_QTV + RG_QW)
   END IF

   IF (rhoEv_new < 0.0D0) rhoEv_new = 0.0D0
   U_C(I)%VAL(1, dimensiona+3) = rhoEv_new

   ! ----------------------------
   ! (C) Species: point-implicit, diagonal Patankar/BE style
   ! ----------------------------
   DO k = 1, nof_species
      rhoY_old = U_C(I)%VAL(1, dimensiona+3+k)

      prod = MAX(RG_DW_S(k), 0.0D0)
      dest = MAX(-RG_DW_S(k), 0.0D0)

      lambda_k = dest / MAX(rhoY_old, rho_min)

      rhoY_new = (rhoY_old + dt_loc * prod) / (1.0D0 + dt_loc * lambda_k)

      IF (rhoY_new < 0.0D0) rhoY_new = 0.0D0
      U_C(I)%VAL(1, dimensiona+3+k) = rhoY_new
   END DO





























END SUBROUTINE SOURCES_realgas















SUBROUTINE SOURCES_derivatives_COMPUTATION(N)
!> @brief
!> Sources derivative computation for implicit time stepping
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE,ICONSIDERED
	REAL,DIMENSION(TURBULENCEEQUATIONS)::SOURCE_T
	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO
	DO I=1,KMAXE
		ICONSIDERED=I
		CALL SOURCES_derivatives(N,ICONSIDERED,SOURCE_T)
		sht(I,1:turbulenceequations)=(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_derivatives_COMPUTATION



SUBROUTINE SOURCES(N,ICONSIDERED,SOURCE_T)
!> @brief
!> Sources computation procedure
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
REAL,DIMENSION(TURBULENCEEQUATIONS),INTENT(INOUT)::SOURCE_T
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX
REAL::VHX,AMP,DVEL,OMEGA,SQUARET,TCH_X,TCH_X3,TCH_FV1,TCH_FV2
REAL::TCH_RS,TCH_R,TCH_G,TCH_GLIM,TCH_FW,TCH_DIF,TCH_DEST,TCH_PROD
INTEGER::I,K,J,L,IHGT,IHGJ,IEX, LOWRE
REAL::SNORM,ONORM,DIVNORM,ax,ay,az,TCH_SHH,TCH_SAV,Verysmall,onesix,ProdTerm1,stild,rr
REAL::gg,FW,destterm,fodt,srcfull,DBPR,DBDI,DBDE,DBY,DBX,ProdTermfinal
REAL:: r_DES,f_DES, ddw,F_DES_SST,L_t_DES
Real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,ProdTerm2,sfac,sss,usss,ssss,S_bar,KRON
real:: uz,vz,wx,wy,wz
REAL:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !For SAS only
REAL,DIMENSION(3,3)::VORTET,TVORT,SVORT,OVORT
REAL,DIMENSION(3)::VORTEM,VERTF,vertex
real,dimension(TURBULENCEEQUATIONS,1:3)::DERIVTURB
!Declarations for k-omega
REAL:: srcfull_k, srcfull_om, Prod_k, Prod_om, Ydest_k, Ydest_om, Diff_om, Q_sas
REAL:: sigma_k, sigma_om, F_1, F_2,Phi_1,Phi_2, D_omplus !-------------- Those are for diffusion too!
REAL:: alpha_raw,alpha_star, Re_t_SST,alpha_inf
REAL:: beta_stari, beta_i, beta_raw, beta_star
REAL:: k_0, om_0, wally
REAL:: dervk_dervom, dervom2, dervk2, u_lapl !Generalization of the velocity Laplacian
REAL:: L_sas, L_vk, Delta_cell, Cell_volume, Q_sas1, Q_sas2
REAL:: kx,ky,kz,omx,omy,omz
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV
REAL::MP_PINFl,GAMMAL
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM

I=ICONSIDERED

Verysmall = 10e-16
	
VORTET(1:3,1:3) = ILOCAL_RECON3(I)%GRADS(1:3,1:3)


ux = Vortet(1,1);uy = Vortet(1,2);uz = Vortet(1,3)
vx = Vortet(2,1);vy = Vortet(2,2);vz = Vortet(2,3)
wx = Vortet(3,1);wy = Vortet(3,2);wz = Vortet(3,3)


DO IHGT=1,3
  DO IHGJ=1,3
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
END DO

sVORT=0.5*(VORTET+TVORT)
OVORT=0.5*(VORTET-TVORT)





SNORM=SQRT(2.0D0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+(SVORT(1,3)*SVORT(1,3))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))+(SVORT(2,3)*SVORT(2,3))+& 
	       (SVORT(3,1)*SVORT(3,1))+(SVORT(3,2)*SVORT(3,2))+(SVORT(3,3)*SVORT(3,3))))
!Quadratic mean of the strain tensor (defined as Svort). Also needed in SST
ONORM=SQRT(2.0D0*((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+(OVORT(1,3)*OVORT(1,3))+&
	       (OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2))+(OVORT(2,3)*OVORT(2,3))+&
	       (OVORT(3,1)*OVORT(3,1))+(OVORT(3,2)*OVORT(3,2))+(OVORT(3,3)*OVORT(3,3))))
OMEGA=ONORM


DIVNORM=ux+uy+uz  !Careful with the sign. If it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)+(wz*wz)))&
	+((uy+vx)*(uy+vx)+(uz+wx)*(uz+wx)+(wy+vz)*(wy+vz))&
	-(2.0/3.0*(ux+vy+wz)*(ux+vy+wz)))


DERIVTURB(1,1:3) = ILOCAL_RECON3(I)%GRADS(5,1:3)
SQUARET=(sqrt((DERIVTURB(1,1)**2)+(DERIVTURB(1,2)**2)+(DERIVTURB(1,3)**2)))**2


LEFTV(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
RIGHTV(1:nof_Variables)=LEFTV(1:nof_Variables)

CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)



 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	
        if (ROT_CORR.EQ.1) then
        OMEGA = OMEGA + 2.0d0*min(0.0d0,SNORM-ONORM)
    end if
    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
	      cw1 = cb1 / (kappa*kappa)
	      cw1 = cw1 + (1.0 + cb2) /sigma
	      TCH_X=   TURBMV(1) / VISCL(1) 
	      
			  if (TCH_X .lt. Verysmall) then
						SOURCE_T(1) = ZERO

			  else
			      TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
			      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
			      TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 


			ddw=IELEM(N,I)%WallDist
			
			if (DES_model .eq. 1) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			ddw=min(ddw,C_DES_SA*Delta_cell)
			end if
			
			if (DES_model .eq. 2) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			r_DES=min(10.0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			!Previous limiter is just for numerical reasons regarding tanh
			f_DES=1-tanh((8.0*r_DES)**3)
			ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			end if

			ProdTerm1 = (TURBMV(1))/(leftv(1)* KAPPA * KAPPA * ddw * ddw)
			Stild = max ( OMEGA + (TCH_fv2*ProdTerm1), 0.3*OMEGA)
			Prodtermfinal=Stild*turbmv(1)*cb1

			RR=MIN((TURBMV(1)/(((LEFTV(1)*STILD*KAPPA * KAPPA * (ddw) * (ddw)))+0.000000001)),10.0)

			!cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))

			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) )/(leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET/ SIGMA
			      
			      if (D_CORR.EQ.0) then
				SOURCE_T(1) = ProdTermfinal + fodt - destterm
            else 
				SOURCE_T(1) = ProdTermfinal + (fodt - destterm)*leftv(1)
            end if
			END IF
	 eLSE
		    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
			  TCH_X=   TURBMV(1) / VISCL(1) 
			    if (tch_x.gt.10.0D0)then
						
			      tch_x=tch_x
						
						
			    ELSE
						
						
				tch_x=0.05D0*log(1.0D0+exp(20.0D0*(tch_x)))
						
						
			      end if
			ddw=IELEM(N,I)%WallDist
			  

			  if (DES_model .eq. 1) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  ddw=min(ddw,C_DES_SA*Delta_cell)
			  end if
			  
			  if (DES_model .eq. 2) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  r_DES=min(10.0D0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			  !Previous limiter is just for numerical reasons regarding tanh
			  f_DES=1-tanh((8.0D0*r_DES)**3)
			  ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			  end if

		      


		  TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
		      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
		  TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 
		  ProdTerm1 = (TCH_fv2*tch_x*(TURBMV(1)))/( leftv(1)*KAPPA * KAPPA * (ddw) * (ddw))
		    
		  IF (PRODTERM1.GE.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ProdTerm1
		  END IF
		  IF (PRODTERM1.LT.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ((OMEGA*(((0.7D0*0.7D0)*(OMEGA))+(0.9*PRODTERM1)))/(((0.9-1.4)*OMEGA)-PRODTERM1))
		  END IF



		  Prodtermfinal=Stild*turbmv(1)*cb1*tch_x
		  ! 
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*Stild*KAPPA * KAPPA * (ddw) * (ddw))))

			cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))
		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
		  ! ! 				!  First order diffusion term
		  fodt  =   LEFTV(1)*cb2 * SQUARET / SIGMA
		  SOURCE_T(1)= ProdTermfinal + fodt - destterm

		  
	 
	 END IF


  CASE(2)		!K OMEGA SST

			      EDDYFL(1)=IELEM(N,I)%WALLDIST
			      EDDYFL(2)=U_CT(I)%VAL(1,1)
			      EDDYFL(3)=U_CT(I)%VAL(1,2)
			      EDDYFL(4:6)=ILOCAL_RECON3(I)%GRADS(1,1:3)
			      EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADS(2,1:3)
			      EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADS(3,1:3)

			      EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADS(4,1:3)
			      EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADS(5,1:3)
								    

			      eddyfr=eddyfl

			      CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
      k_0=(MAX(Verysmall,U_CT(I)%VAL(1,1)/LEFTV(1))) !First subindex makes reference to the time-stepping
      om_0=MAX(1.0e-1*ufreestream/CharLength,U_CT(I)%VAL(1,2)/LEFTV(1))
      wally=IELEM(N,I)%WallDist
	

	      !Calculate here k and omega gradients
		  OMX=ILOCAL_RECON3(I)%GRADS(6,1)
		  OMY=ILOCAL_RECON3(I)%GRADS(6,2)
		  OMZ=ILOCAL_RECON3(I)%GRADS(6,3)
		  KX=ILOCAL_RECON3(I)%GRADS(5,1)
		  KY=ILOCAL_RECON3(I)%GRADS(5,2)
		  KZ=ILOCAL_RECON3(I)%GRADS(5,3)

		  dervk_dervom= kx*omx+ky*omy+kz*omz

			!Parameters

			!-----NEEDED IN DIFFUSION--------
			D_omplus=max(2*LEFTV(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(1)/(LEFTV(1)*wally*wally*om_0))
			phi_1=min(phi_2, 4.0*LEFTV(1)*k_0/(sigma_om2*D_omplus*wally*wally)) 

				F_1=tanh(phi_1**4)
				F_2=tanh(phi_2**2)
				!--------------------------------


				alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
				alpha_star=alpha_starinf
				alpha_raw=alpha_inf


				beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
				alpha_star0=beta_i/3.0
				beta_stari=beta_starinf




				LOWRE=0
				!LOW-RE correction------------------------------------------------------------------

				if (lowre.eq.1) then
				Re_t_SST=LEFTV(1)*k_0/(VISCL(1)*om_0)  !Limiters for this???

				alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)
				alpha_raw=alpha_inf/alpha_star*(alpha_0+Re_t_SST/R_om_SST)/(1.0+Re_t_SST/R_om_SST)

				beta_stari=beta_starinf*(4.0/15.0+(Re_t_SST/R_beta)**4)/(1.0+(Re_t_SST/R_beta)**4)

				end if
				      !--------------------------------------------------------------------------

				      !No Mach number corrections for the beta
				      beta_star=beta_stari
				      beta_raw=beta_i


				      !PRODUCTION TERMS
					      !Production of k
					      if (vort_model.eq.1) then
						      Prod_k=VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM
						      Prod_om=alpha_raw*LEFTV(1)*SNORM*ONORM
					      else
						      !Prod_k=VISCL(3)*SNORM**2
						      !Exact formulation: 
						      Prod_k=VISCL(3)*(SNORM**2)!-2.0/3.0*DIVNORM*DIVNORM)-2.0/3.0*LEFTV(1)*k_0*DIVNORM
						      Prod_k=min(Prod_k,10.0*LEFTV(1)*beta_star*k_0*om_0)

						      !INTELLIGENT WAY OF LIMITING: 
						      !Prod_k=min(Prod_k,max(10.0*LEFTV(1)*beta_star*k_0*om_0,&
						      !	    VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM))
					      
					      
						      !Production of omega  (Menter does this before correcting Prod_k, 
						      !but in  article of 2003 he applies the correction to both)	
						      Prod_om=alpha_raw*LEFTV(1)*SNORM**2
				      ! 		Prod_om=min(Prod_om,10.0*LEFTV(1)*beta_star*om_0*om_0)
					      end if

				      !DESTRUCTION TERMS
					      !Destruction of k
					      Ydest_k=LEFTV(1)*beta_star*k_0*om_0
					      !Destruction of omega
					      Ydest_om=LEFTV(1)*beta_raw*om_0*om_0
					      
				      !CROSSED-DIFFUSION TERM
					      !Crossed diffusion of omega
					      Diff_om=2.0*(1-F_1)*LEFTV(1)/(om_0*sigma_om2)*dervk_dervom
	 
	
				!QSAS TERM: SCALE ADAPTIVE
					if (QSAS_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
						!Calculate here second derivative of u 
						!Declare all variables
						      DO IEX=1,3
							VORTET(IEX,1:3)=ILOCAL_RECON3(I)%GRADS(3+TURBULENCEEQUATIONS+IEX,1:3)
							
						      END DO
						
						uxx=VORTET(1,1) ; uyy=VORTET(1,2); uzz=VORTET(1,3)				
						vxx=VORTET(2,1) ; vyy=VORTET(2,2); vzz=VORTET(2,3)	
						wxx=VORTET(3,1) ; wyy=VORTET(3,2); wzz=VORTET(3,3)	
						!Compute here cell volume
						CELL_VOLUME=IELEM(N,I)%TOTVOLUME
						
						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
						dervk2=kx*kx+ky*ky+kz*kz
						dervom2=omx*omx+omy*omy+omz*omz
						
						Delta_cell=Cell_volume**0.333333333333333
						
						L_sas=sqrt(k_0)/(beta_star**0.25*om_0)
						L_vk=max(kappa*SNORM/u_lapl, &       !This switch provides high wave-number damping
							C_smg*Delta_cell*sqrt(kappa*eta2_SAS/(beta_raw/beta_star-alpha_raw)))
						
						Q_sas1=LEFTV(1)*eta2_SAS*kappa*SNORM**2*(L_sas/L_vk)**2
						Q_sas2= -C_SAS*2*LEFTV(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)

						Q_SAS=max(Q_sas1+Q_sas2,0.0)
					else
					Q_SAS=ZERO
					end if

				    !FINAL SOURCE TERMS


					    !DES-SST MODEL (if QSAS_model=2)
					    if (QSAS_model .eq.2) then
					    CELL_VOLUME=IELEM(N,I)%TOTVOLUME
					    Delta_cell=Cell_volume**0.333333333333333
					    L_t_DES=sqrt(k_0)/(beta_star*om_0)
					    F_DES_SST=max(1.0, L_t_DES/(C_DES_SST*Delta_cell)*(1-F_2))
					    !The (1-F_2) is meant to protect the boundary layer. Will result in same 
					    !separation point that standard S-A
					    Ydest_k=F_DES_SST*Ydest_k
					    end if
					    
				    srcfull_k=Prod_k-Ydest_k
				    srcfull_om=Prod_om-Ydest_om+Diff_om+Q_SAS


				    !Filling the output vector
				      SOURCE_T(1) = srcfull_k
				      SOURCE_T(2) = srcfull_om

				  
	
END SELECT

END SUBROUTINE SOURCES


SUBROUTINE SOURCES_DERIVATIVES(N,ICONSIDERED,SOURCE_t)
!> @brief
!> Sources derivative computation for implicit time stepping
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
REAL,DIMENSION(TURBULENCEEQUATIONS),INTENT(INOUT)::SOURCE_T
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX
REAL::VHX,AMP,DVEL,OMEGA,SQUARET,TCH_X,TCH_X3,TCH_FV1,TCH_FV2
REAL::TCH_RS,TCH_R,TCH_G,TCH_GLIM,TCH_FW,TCH_DIF,TCH_DEST,TCH_PROD
INTEGER::I,K,J,L,IHGT,IHGJ,IEX, LOWRE
REAL::SNORM,ONORM,DIVNORM,ax,ay,az,TCH_SHH,TCH_SAV,Verysmall,onesix,ProdTerm1,stild,rr
REAL::gg,FW,destterm,fodt,srcfull,DBPR,DBDI,DBDE,DBY,DBX,ProdTermfinal
REAL:: r_DES,f_DES, ddw,F_DES_SST,L_t_DES
Real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,ProdTerm2,sfac,sss,usss,ssss,S_bar,KRON
real:: uz,vz,wx,wy,wz
REAL:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !For SAS only
REAL,DIMENSION(3,3)::VORTET,TVORT,SVORT,OVORT
REAL,DIMENSION(3)::VORTEM,VERTF,vertex
real,dimension(TURBULENCEEQUATIONS,1:3)::DERIVTURB
!Declarations for k-omega
REAL:: srcfull_k, srcfull_om, Prod_k, Prod_om, Ydest_k, Ydest_om, Diff_om, Q_sas
REAL:: sigma_k, sigma_om, F_1, F_2,Phi_1,Phi_2, D_omplus !-------------- Those are for diffusion too!
REAL:: alpha_raw,alpha_star, Re_t_SST,alpha_inf
REAL:: beta_stari, beta_i, beta_raw, beta_star
REAL:: k_0, om_0, wally
REAL:: dervk_dervom, dervom2, dervk2, u_lapl !Generalization of the velocity Laplacian
REAL:: L_sas, L_vk, Delta_cell, Cell_volume, Q_sas1, Q_sas2
REAL:: kx,ky,kz,omx,omy,omz
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV
REAL::MP_PINFl,GAMMAL
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM


I=ICONSIDERED

Verysmall = 10e-16
	
VORTET(1:3,1:3) = ILOCAL_RECON3(I)%GRADS(1:3,1:3)


ux = Vortet(1,1);uy = Vortet(1,2);uz = Vortet(1,3)
vx = Vortet(2,1);vy = Vortet(2,2);vz = Vortet(2,3)
wx = Vortet(3,1);wy = Vortet(3,2);wz = Vortet(3,3)


DO IHGT=1,3
  DO IHGJ=1,3
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
END DO

sVORT=0.5*(VORTET+TVORT)
OVORT=0.5*(VORTET-TVORT)





SNORM=SQRT(2.0D0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+(SVORT(1,3)*SVORT(1,3))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))+(SVORT(2,3)*SVORT(2,3))+& 
	       (SVORT(3,1)*SVORT(3,1))+(SVORT(3,2)*SVORT(3,2))+(SVORT(3,3)*SVORT(3,3))))
!Quadratic mean of the strain tensor (defined as Svort). Also needed in SST
ONORM=SQRT(2.0D0*((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+(OVORT(1,3)*OVORT(1,3))+&
	       (OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2))+(OVORT(2,3)*OVORT(2,3))+&
	       (OVORT(3,1)*OVORT(3,1))+(OVORT(3,2)*OVORT(3,2))+(OVORT(3,3)*OVORT(3,3))))
OMEGA=ONORM


DIVNORM=ux+uy+uz  !Careful with the sign. If it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)+(wz*wz)))&
	+((uy+vx)*(uy+vx)+(uz+wx)*(uz+wx)+(wy+vz)*(wy+vz))&
	-(2.0/3.0*(ux+vy+wz)*(ux+vy+wz)))


DERIVTURB(1,1:3) = ILOCAL_RECON3(I)%GRADS(5,1:3)
SQUARET=(sqrt((DERIVTURB(1,1)**2)+(DERIVTURB(1,2)**2)+(DERIVTURB(1,3)**2)))**2


LEFTV(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
RIGHTV(1:nof_Variables)=LEFTV(1:nof_Variables)
CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)



 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
	      cw1 = cb1 / (kappa*kappa)
	      cw1 = cw1 + (1.0 + cb2) /sigma
	      TCH_X=   TURBMV(1) / VISCL(1) 
	      
			  if (TCH_X .lt. Verysmall) then
						SOURCE_T(1) = ZERO

			  else
			      TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
			      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
			      TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 


			ddw=IELEM(N,I)%WallDist
			
			if (DES_model .eq. 1) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			ddw=min(ddw,C_DES_SA*Delta_cell)
			end if
			
			if (DES_model .eq. 2) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			r_DES=min(10.0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			!Previous limiter is just for numerical reasons regarding tanh
			f_DES=1-tanh((8.0*r_DES)**3)
			ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			end if

			ProdTerm1 = (TURBMV(1))/(leftv(1)* KAPPA * KAPPA * ddw * ddw)
			Stild = max ( OMEGA + (TCH_fv2*ProdTerm1), 0.3*OMEGA)
			Prodtermfinal=Stild*turbmv(1)*cb1/leftv(1)

			RR=MIN((TURBMV(1)/(((LEFTV(1)*STILD*KAPPA * KAPPA * (ddw) * (ddw)))+0.000000001)),10.0)
			!cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))

			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) ) /( leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET / SIGMA
			      
			      
			        Prodtermfinal=Stild*cb1
			         
			      destterm  = 2.0* cw1 * fw * (TURBMV(1)/(LEFTV(1)*DDW**2))
			      fodt  =   2.0* cb2 * (SQRT(SQUARET)) / SIGMA
			      		      
			      
			      SOURCE_T(1)=  min(ProdTermfinal + fodt - destterm,ZERO)
			END IF
	 eLSE
		    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
			  TCH_X=   TURBMV(1) / VISCL(1) 
			    if (tch_x.gt.10.0D0)then
						
			      tch_x=tch_x
						
						
			    ELSE
						
						
				tch_x=0.05D0*log(1.0D0+exp(20.0D0*(tch_x)))
						
						
			      end if
			ddw=IELEM(N,I)%WallDist
			  

			  if (DES_model .eq. 1) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  ddw=min(ddw,C_DES_SA*Delta_cell)
			  end if
			  
			  if (DES_model .eq. 2) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  r_DES=min(10.0D0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			  !Previous limiter is just for numerical reasons regarding tanh
			  f_DES=1-tanh((8.0D0*r_DES)**3)
			  ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			  end if

		      


		  TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
		      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
		  TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 
		  ProdTerm1 = (TCH_fv2*tch_x*(TURBMV(1)))/( leftv(1)*KAPPA * KAPPA * (ddw) * (ddw))
		    
		  IF (PRODTERM1.GE.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ProdTerm1
		  END IF
		  IF (PRODTERM1.LT.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ((OMEGA*(((0.7D0*0.7D0)*(OMEGA))+(0.9*PRODTERM1)))/(((0.9-1.4)*OMEGA)-PRODTERM1))
		  END IF



		  Prodtermfinal=Stild*turbmv(1)*cb1*tch_x
		  ! 
		  ! 
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*Stild*KAPPA * KAPPA * (ddw) * (ddw))))
		!	cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))

		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
		  ! ! 				!  First order diffusion term
		  fodt  =   LEFTV(1)*cb2 * SQUARET / SIGMA
		  
			      			      		      
			   Prodtermfinal=Stild*cb1*tch_x
			         
			      destterm  = 2.0* cw1 * fw * (TURBMV(1)/(LEFTV(1)*DDW**2))
			      fodt  =   2.0* cb2 * (SQRT(SQUARET)) / SIGMA
			      		      
			      
			      SOURCE_T(1)=  min(ProdTermfinal + fodt - destterm,ZERO)
		  
	 
	 END IF


  CASE(2)		!K OMEGA SST

			      EDDYFL(1)=IELEM(N,I)%WALLDIST
			      EDDYFL(2)=U_CT(I)%VAL(1,1)
			      EDDYFL(3)=U_CT(I)%VAL(1,2)
			      EDDYFL(4:6)=ILOCAL_RECON3(I)%GRADS(1,1:3)
			      EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADS(2,1:3)
			      EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADS(3,1:3)

			      EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADS(4,1:3)
			      EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADS(5,1:3)
								    

			      eddyfr=eddyfl

			      CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
      k_0=(MAX(Verysmall,U_CT(I)%VAL(1,1)/LEFTV(1))) !First subindex makes reference to the time-stepping
      om_0=MAX(1.0e-1*ufreestream/CharLength,U_CT(I)%VAL(1,2)/LEFTV(1))
!       wally=IELEM(N,I)%WallDist
! 	
! 
! 	      !Calculate here k and omega gradients
! 		  OMX=ILOCAL_RECON3(I)%GRADS(5,1)
! 		  OMY=ILOCAL_RECON3(I)%GRADS(5,2)
! 		  OMZ=ILOCAL_RECON3(I)%GRADS(5,3)
! 		  KX=ILOCAL_RECON3(I)%GRADS(4,1)
! 		  KY=ILOCAL_RECON3(I)%GRADS(4,2)
! 		  KZ=ILOCAL_RECON3(I)%GRADS(4,3)
! 
! 		  dervk_dervom= kx*omx+ky*omy+kz*omz
! 
! 			!Parameters
! 
! 			!-----NEEDED IN DIFFUSION--------
! 			D_omplus=max(2*LEFTV(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
! 			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(1)/(LEFTV(1)*wally*wally*om_0))
! 			phi_1=min(phi_2, 4.0*LEFTV(1)*k_0/(sigma_om2*D_omplus*wally*wally)) 
! 
! 				F_1=tanh(phi_1**4)
! 				F_2=tanh(phi_2**2)
! 				!--------------------------------
! 
! 
! 				alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
! 				alpha_star=alpha_starinf
! 				alpha_raw=alpha_inf
! 
! 
! 				beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
! 				alpha_star0=beta_i/3.0
! 				beta_stari=beta_starinf
! 
! 
! 
! 
! 				LOWRE=0
! 				!LOW-RE correction------------------------------------------------------------------
! 
! 				if (lowre.eq.1) then
! 				Re_t_SST=LEFTV(1)*k_0/(VISCL(1)*om_0)  !Limiters for this???
! 
! 				alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)
! 				alpha_raw=alpha_inf/alpha_star*(alpha_0+Re_t_SST/R_om_SST)/(1.0+Re_t_SST/R_om_SST)
! 
! 				beta_stari=beta_starinf*(4.0/15.0+(Re_t_SST/R_beta)**4)/(1.0+(Re_t_SST/R_beta)**4)
! 
! 				end if
! 				      !--------------------------------------------------------------------------
! 
! 				      !No Mach number corrections for the beta
! 				      beta_star=beta_stari
! 				      beta_raw=beta_i
! 
! 
! 				      !PRODUCTION TERMS
! 					      !Production of k
! 					      if (vort_model.eq.1) then
! 						      Prod_k=VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM
! 						      Prod_om=alpha_raw*LEFTV(1)*SNORM*ONORM
! 					      else
! 						      !Prod_k=VISCL(3)*SNORM**2
! 						      !Exact formulation: 
! 						      Prod_k=VISCL(3)*(SNORM**2)!-2.0/3.0*DIVNORM*DIVNORM)-2.0/3.0*LEFTV(1)*k_0*DIVNORM
! 						      Prod_k=min(Prod_k,10.0*LEFTV(1)*beta_star*k_0*om_0)
! 
! 						      !INTELLIGENT WAY OF LIMITING: 
! 						      !Prod_k=min(Prod_k,max(10.0*LEFTV(1)*beta_star*k_0*om_0,&
! 						      !	    VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM))
! 					      
! 					      
! 						      !Production of omega  (Menter does this before correcting Prod_k, 
! 						      !but in  article of 2003 he applies the correction to both)	
! 						      Prod_om=alpha_raw*LEFTV(1)*SNORM**2
! 				      ! 		Prod_om=min(Prod_om,10.0*LEFTV(1)*beta_star*om_0*om_0)
! 					      end if
! 
! 				      !DESTRUCTION TERMS
! 					      !Destruction of k
! 					      Ydest_k=LEFTV(1)*beta_star*k_0*om_0
! 					      !Destruction of omega
! 					      Ydest_om=LEFTV(1)*beta_raw*om_0*om_0
! 					      
! 				      !CROSSED-DIFFUSION TERM
! 					      !Crossed diffusion of omega
! 					      Diff_om=2.0*(1-F_1)*LEFTV(1)/(om_0*sigma_om2)*dervk_dervom
! 	 
! 	
! 				!QSAS TERM: SCALE ADAPTIVE
! 					if (QSAS_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
! 						!Calculate here second derivative of u 
! 						!Declare all variables
! 						      DO IEX=1,3
! 							VORTET(IEX,1:3)=ILOCAL_RECON3(I)%GRADS(3+TURBULENCEEQUATIONS+IEX,1:3)
! 							
! 						      END DO
! 						
! 						uxx=VORTET(1,1) ; uyy=VORTET(1,2); uzz=VORTET(1,3)				
! 						vxx=VORTET(2,1) ; vyy=VORTET(2,2); vzz=VORTET(2,3)	
! 						wxx=VORTET(3,1) ; wyy=VORTET(3,2); wzz=VORTET(3,3)	
! 						!Compute here cell volume
! 						CELL_VOLUME=IELEM(N,I)%TOTVOLUME
! 						
! 						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
! 						dervk2=kx*kx+ky*ky+kz*kz
! 						dervom2=omx*omx+omy*omy+omz*omz
! 						
! 						Delta_cell=Cell_volume**0.333333333333333
! 						
! 						L_sas=sqrt(k_0)/(beta_star**0.25*om_0)
! 						L_vk=max(kappa*SNORM/u_lapl, &       !This switch provides high wave-number damping
! 							C_smg*Delta_cell*sqrt(kappa*eta2_SAS/(beta_raw/beta_star-alpha_raw)))
! 						
! 						Q_sas1=LEFTV(1)*eta2_SAS*kappa*SNORM**2*(L_sas/L_vk)**2
! 						Q_sas2= -C_SAS*2*LEFTV(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)
! 
! 						Q_SAS=max(Q_sas1+Q_sas2,0.0)
! 					else
! 					Q_SAS=ZERO
! 					end if
! 
! 				    !FINAL SOURCE TERMS
! 
! 
! 					    !DES-SST MODEL (if QSAS_model=2)
! 					    if (QSAS_model .eq.2) then
! 					    CELL_VOLUME=IELEM(N,I)%TOTVOLUME
! 					    Delta_cell=Cell_volume**0.333333333333333
! 					    L_t_DES=sqrt(k_0)/(beta_star*om_0)
! 					    F_DES_SST=max(1.0, L_t_DES/(C_DES_SST*Delta_cell)*(1-F_2))
! 					    !The (1-F_2) is meant to protect the boundary layer. Will result in same 
! 					    !separation point that standard S-A
! 					    Ydest_k=F_DES_SST*Ydest_k
! 					    end if
! 					    
! 				    srcfull_k=Prod_k-Ydest_k
! 				    srcfull_om=Prod_om-Ydest_om+Diff_om+Q_SAS


				    !Filling the output vector
				      SOURCE_T(1) = BETA_STARINF*om_0
				      SOURCE_T(2) = BETA_STARINF*k_0

				  
	
END SELECT

END SUBROUTINE SOURCES_DERIVATIVES

SUBROUTINE SOURCES_COMPUTATION2d(N)
!> @brief
!> Sources  computation in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE,ICONSIDERED
	REAL,DIMENSION(TURBULENCEEQUATIONS)::SOURCE_T

	
	
	KMAXE=XMPIELRANK(N)
	
	
	!$OMP DO
	DO I=1,KMAXE
		ICONSIDERED=I

		if (turbulence.eq.1)then

		CALL SOURCES2d(N,ICONSIDERED,SOURCE_T)
		RHST(I)%VAL(1:turbulenceequations)=RHST(I)%VAL(1:turbulenceequations)-(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)

		end if



	END DO
	!$OMP END DO 
	
	
END SUBROUTINE SOURCES_COMPUTATION2d

SUBROUTINE SOURCES_derivatives_COMPUTATION2d(N)
!> @brief
!> Sources  derivatives computation in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE,ICONSIDERED
	REAL,DIMENSION(TURBULENCEEQUATIONS)::SOURCE_T
	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO
	DO I=1,KMAXE
		ICONSIDERED=I
		CALL SOURCES_derivatives2d(N,ICONSIDERED,SOURCE_T)
		sht(I,1:turbulenceequations)=(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_derivatives_COMPUTATION2d



SUBROUTINE SOURCES2d(N,ICONSIDERED,SOURCE_T)
!> @brief
!> Sources  computation in 2D
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
REAL,DIMENSION(TURBULENCEEQUATIONS),INTENT(INOUT)::SOURCE_T
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX
REAL::VHX,AMP,DVEL,OMEGA,SQUARET,TCH_X,TCH_X3,TCH_FV1,TCH_FV2
REAL::TCH_RS,TCH_R,TCH_G,TCH_GLIM,TCH_FW,TCH_DIF,TCH_DEST,TCH_PROD
INTEGER::I,K,J,L,IHGT,IHGJ,IEX, LOWRE
REAL::SNORM,ONORM,DIVNORM,ax,ay,az,TCH_SHH,TCH_SAV,Verysmall,onesix,ProdTerm1,stild,rr
REAL::gg,FW,destterm,fodt,srcfull,DBPR,DBDI,DBDE,DBY,DBX,ProdTermfinal
REAL:: r_DES,f_DES, ddw,F_DES_SST,L_t_DES
Real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,ProdTerm2,sfac,sss,usss,ssss,S_bar,KRON
real:: uz,vz,wx,wy,wz
REAL:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !For SAS only
REAL,DIMENSION(2,2)::VORTET,TVORT,SVORT,OVORT
REAL,DIMENSION(2)::VORTEM,VERTF,vertex
real,dimension(TURBULENCEEQUATIONS,1:2)::DERIVTURB
!Declarations for k-omega
REAL:: srcfull_k, srcfull_om, Prod_k, Prod_om, Ydest_k, Ydest_om, Diff_om, Q_sas
REAL:: sigma_k, sigma_om, F_1, F_2,Phi_1,Phi_2, D_omplus !-------------- Those are for diffusion too!
REAL:: alpha_raw,alpha_star, Re_t_SST,alpha_inf
REAL:: beta_stari, beta_i, beta_raw, beta_star
REAL:: k_0, om_0, wally
REAL:: dervk_dervom, dervom2, dervk2, u_lapl !Generalization of the velocity Laplacian
REAL:: L_sas, L_vk, Delta_cell, Cell_volume, Q_sas1, Q_sas2
REAL:: kx,ky,kz,omx,omy,omz
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV
REAL::MP_PINFl,GAMMAL,MP_PINFR,GAMMAR
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM

I=ICONSIDERED

Verysmall = 10e-16



	
VORTET(1:2,1:2) = ILOCAL_RECON3(I)%GRADS(1:2,1:2)


ux = Vortet(1,1);uy = Vortet(1,2)
vx = Vortet(2,1);vy = Vortet(2,2)







DO IHGT=1,2
  DO IHGJ=1,2
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
END DO

sVORT=0.5*(VORTET+TVORT)
OVORT=0.5*(VORTET-TVORT)





SNORM=SQRT(2.0D0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))))
!Quadratic mean of the strain tensor (defined as Svort). Also needed in SST
ONORM=SQRT(2.0D0*((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+&
	       (OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2))))
OMEGA=ONORM


DIVNORM=ux+uy+uz  !Careful with the sign. If it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)))&
	+((uy+vx)*(uy+vx))&
	-(2.0/3.0*(ux+vy)*(ux+vy)))


DERIVTURB(1,1:2) = ILOCAL_RECON3(I)%GRADS(4,1:2)
SQUARET=(sqrt((DERIVTURB(1,1)**2)+(DERIVTURB(1,2)**2)))**2


LEFTV(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
RIGHTV(1:nof_Variables)=LEFTV(1:nof_Variables)

CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)



 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
	      cw1 = cb1 / (kappa*kappa)
	      cw1 = cw1 + (1.0 + cb2) /sigma
	      TCH_X=   TURBMV(1) / VISCL(1) 
	      
			  if (TCH_X .lt. Verysmall) then
						SOURCE_T(1) = ZERO

			  else
			      TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
			      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
			      TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 


			ddw=IELEM(N,I)%WallDist
			
			if (DES_model .eq. 1) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			ddw=min(ddw,C_DES_SA*Delta_cell)
			end if
			
			if (DES_model .eq. 2) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			r_DES=min(10.0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			!Previous limiter is just for numerical reasons regarding tanh
			f_DES=1-tanh((8.0*r_DES)**3)
			ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			end if

			ProdTerm1 = (TURBMV(1))/(leftv(1)* KAPPA * KAPPA * ddw * ddw)
			Stild = max ( OMEGA + (TCH_fv2*ProdTerm1), 0.3*OMEGA)
			Prodtermfinal=Stild*turbmv(1)*cb1

			RR=MIN((TURBMV(1)/(((LEFTV(1)*STILD*KAPPA * KAPPA * (ddw) * (ddw)))+0.000000001)),10.0)

			!cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))
			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) )/(leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET/ SIGMA
			      SOURCE_T(1) = ProdTermfinal + fodt - destterm
			     
			     
			END IF
	 eLSE
		    CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
			  TCH_X=   TURBMV(1) / VISCL(1) 
			    if (tch_x.gt.10.0D0)then
						
			      tch_x=tch_x
						
						
			    ELSE
						
						
				tch_x=0.05D0*log(1.0D0+exp(20.0D0*(tch_x)))
						
						
			      end if
			ddw=IELEM(N,I)%WallDist
			  

			  if (DES_model .eq. 1) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  ddw=min(ddw,C_DES_SA*Delta_cell)
			  end if
			  
			  if (DES_model .eq. 2) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  r_DES=min(10.0D0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			  !Previous limiter is just for numerical reasons regarding tanh
			  f_DES=1-tanh((8.0D0*r_DES)**3)
			  ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			  end if

		      


		  TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
		      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
		  TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 
		  ProdTerm1 = (TCH_fv2*tch_x*(TURBMV(1)))/( leftv(1)*KAPPA * KAPPA * (ddw) * (ddw))
		    
		  IF (PRODTERM1.GE.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ProdTerm1
		  END IF
		  IF (PRODTERM1.LT.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ((OMEGA*(((0.7D0*0.7D0)*(OMEGA))+(0.9*PRODTERM1)))/(((0.9-1.4)*OMEGA)-PRODTERM1))
		  END IF



		  Prodtermfinal=Stild*turbmv(1)*cb1*tch_x
		  ! 
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*Stild*KAPPA * KAPPA * (ddw) * (ddw))))

			cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))
		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
		  ! ! 				!  First order diffusion term
		  fodt  =   LEFTV(1)*cb2 * SQUARET / SIGMA
		  SOURCE_T(1)= ProdTermfinal + fodt - destterm

		  
	 
	 END IF


  CASE(2)		!K OMEGA SST

			      EDDYFL(1)=IELEM(N,I)%WALLDIST
			      EDDYFL(2)=U_CT(I)%VAL(1,1)
			      EDDYFL(3)=U_CT(I)%VAL(1,2)
			      EDDYFL(4:5)=ILOCAL_RECON3(I)%GRADS(1,1:2)
			      EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADS(2,1:2)
			      EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADS(4,1:2)
			      EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADS(5,1:2)
								    

			      eddyfr=eddyfl

			      CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
      k_0=(MAX(Verysmall,U_CT(I)%VAL(1,1)/LEFTV(1))) !First subindex makes reference to the time-stepping
      om_0=MAX(1.0e-1*ufreestream/CharLength,U_CT(I)%VAL(1,2)/LEFTV(1))
      wally=IELEM(N,I)%WallDist
	

	      !Calculate here k and omega gradients
		  OMX=ILOCAL_RECON3(I)%GRADS(5,1)
		  OMY=ILOCAL_RECON3(I)%GRADS(5,2)
		  
		  KX=ILOCAL_RECON3(I)%GRADS(4,1)
		  KY=ILOCAL_RECON3(I)%GRADS(4,2)
		

		  dervk_dervom= kx*omx+ky*omy

			!Parameters

			!-----NEEDED IN DIFFUSION--------
			D_omplus=max(2*LEFTV(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(1)/(LEFTV(1)*wally*wally*om_0))
			phi_1=min(phi_2, 4.0*LEFTV(1)*k_0/(sigma_om2*D_omplus*wally*wally)) 

				F_1=tanh(phi_1**4)
				F_2=tanh(phi_2**2)
				!--------------------------------


				alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
				alpha_star=alpha_starinf
				alpha_raw=alpha_inf


				beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
				alpha_star0=beta_i/3.0
				beta_stari=beta_starinf




				LOWRE=0
				!LOW-RE correction------------------------------------------------------------------

				if (lowre.eq.1) then
				Re_t_SST=LEFTV(1)*k_0/(VISCL(1)*om_0)  !Limiters for this???

				alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)
				alpha_raw=alpha_inf/alpha_star*(alpha_0+Re_t_SST/R_om_SST)/(1.0+Re_t_SST/R_om_SST)

				beta_stari=beta_starinf*(4.0/15.0+(Re_t_SST/R_beta)**4)/(1.0+(Re_t_SST/R_beta)**4)

				end if
				      !--------------------------------------------------------------------------

				      !No Mach number corrections for the beta
				      beta_star=beta_stari
				      beta_raw=beta_i


				      !PRODUCTION TERMS
					      !Production of k
					      if (vort_model.eq.1) then
						      Prod_k=VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM
						      Prod_om=alpha_raw*LEFTV(1)*SNORM*ONORM
					      else
						      !Prod_k=VISCL(3)*SNORM**2
						      !Exact formulation: 
						      Prod_k=VISCL(3)*(SNORM**2)!-2.0/3.0*DIVNORM*DIVNORM)-2.0/3.0*LEFTV(1)*k_0*DIVNORM
						      Prod_k=min(Prod_k,10.0*LEFTV(1)*beta_star*k_0*om_0)

						      !INTELLIGENT WAY OF LIMITING: 
						      !Prod_k=min(Prod_k,max(10.0*LEFTV(1)*beta_star*k_0*om_0,&
						      !	    VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM))
					      
					      
						      !Production of omega  (Menter does this before correcting Prod_k, 
						      !but in  article of 2003 he applies the correction to both)	
						      Prod_om=alpha_raw*LEFTV(1)*SNORM**2
				      ! 		Prod_om=min(Prod_om,10.0*LEFTV(1)*beta_star*om_0*om_0)
					      end if

				      !DESTRUCTION TERMS
					      !Destruction of k
					      Ydest_k=LEFTV(1)*beta_star*k_0*om_0
					      !Destruction of omega
					      Ydest_om=LEFTV(1)*beta_raw*om_0*om_0
					      
				      !CROSSED-DIFFUSION TERM
					      !Crossed diffusion of omega
					      Diff_om=2.0*(1-F_1)*LEFTV(1)/(om_0*sigma_om2)*dervk_dervom
	 
	
				!QSAS TERM: SCALE ADAPTIVE
					if (QSAS_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
						!Calculate here second derivative of u 
						!Declare all variables
						      DO IEX=1,2
							VORTET(IEX,1:2)=ILOCAL_RECON3(I)%GRADS(3+TURBULENCEEQUATIONS+IEX,1:2)
							
						      END DO
						
						uxx=VORTET(1,1) ; uyy=VORTET(1,2); 				
						vxx=VORTET(2,1) ; vyy=VORTET(2,2); 	
						!Compute here cell volume
						CELL_VOLUME=IELEM(N,I)%TOTVOLUME
						
						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
						dervk2=kx*kx+ky*ky+kz*kz
						dervom2=omx*omx+omy*omy+omz*omz
						
						Delta_cell=Cell_volume**0.333333333333333
						
						L_sas=sqrt(k_0)/(beta_star**0.25*om_0)
						L_vk=max(kappa*SNORM/u_lapl, &       !This switch provides high wave-number damping
							C_smg*Delta_cell*sqrt(kappa*eta2_SAS/(beta_raw/beta_star-alpha_raw)))
						
						Q_sas1=LEFTV(1)*eta2_SAS*kappa*SNORM**2*(L_sas/L_vk)**2
						Q_sas2= -C_SAS*2*LEFTV(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)

						Q_SAS=max(Q_sas1+Q_sas2,0.0)
					else
					Q_SAS=ZERO
					end if

				    !FINAL SOURCE TERMS


					    !DES-SST MODEL (if QSAS_model=2)
					    if (QSAS_model .eq.2) then
					    CELL_VOLUME=IELEM(N,I)%TOTVOLUME
					    Delta_cell=Cell_volume**0.333333333333333
					    L_t_DES=sqrt(k_0)/(beta_star*om_0)
					    F_DES_SST=max(1.0, L_t_DES/(C_DES_SST*Delta_cell)*(1-F_2))
					    !The (1-F_2) is meant to protect the boundary layer. Will result in same 
					    !separation point that standard S-A
					    Ydest_k=F_DES_SST*Ydest_k
					    end if
					    
				    srcfull_k=Prod_k-Ydest_k
				    srcfull_om=Prod_om-Ydest_om+Diff_om+Q_SAS


				    !Filling the output vector
				      SOURCE_T(1) = srcfull_k
				      SOURCE_T(2) = srcfull_om

				  
	
END SELECT





END SUBROUTINE SOURCES2d


SUBROUTINE SOURCES_DERIVATIVES2d(N,ICONSIDERED,SOURCE_t)
!> @brief
!> Sources derivatives computation in 2D
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
REAL,DIMENSION(TURBULENCEEQUATIONS),INTENT(INOUT)::SOURCE_T
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX
REAL::VHX,AMP,DVEL,OMEGA,SQUARET,TCH_X,TCH_X3,TCH_FV1,TCH_FV2
REAL::TCH_RS,TCH_R,TCH_G,TCH_GLIM,TCH_FW,TCH_DIF,TCH_DEST,TCH_PROD
INTEGER::I,K,J,L,IHGT,IHGJ,IEX, LOWRE
REAL::SNORM,ONORM,DIVNORM,ax,ay,az,TCH_SHH,TCH_SAV,Verysmall,onesix,ProdTerm1,stild,rr
REAL::gg,FW,destterm,fodt,srcfull,DBPR,DBDI,DBDE,DBY,DBX,ProdTermfinal
REAL:: r_DES,f_DES, ddw,F_DES_SST,L_t_DES
Real :: ux,uy,vx,vy,shear,sratio,prodmod,cvor,stildmod,ProdTerm2,sfac,sss,usss,ssss,S_bar,KRON
real:: uz,vz,wx,wy,wz
REAL:: uxx,uyy,uzz,vxx,vyy,vzz,wxx,wyy,wzz  !For SAS only
REAL,DIMENSION(2,2)::VORTET,TVORT,SVORT,OVORT
REAL,DIMENSION(2)::VORTEM,VERTF,vertex
real,dimension(TURBULENCEEQUATIONS,1:2)::DERIVTURB
!Declarations for k-omega
REAL:: srcfull_k, srcfull_om, Prod_k, Prod_om, Ydest_k, Ydest_om, Diff_om, Q_sas
REAL:: sigma_k, sigma_om, F_1, F_2,Phi_1,Phi_2, D_omplus !-------------- Those are for diffusion too!
REAL:: alpha_raw,alpha_star, Re_t_SST,alpha_inf
REAL:: beta_stari, beta_i, beta_raw, beta_star
REAL:: k_0, om_0, wally
REAL:: dervk_dervom, dervom2, dervk2, u_lapl !Generalization of the velocity Laplacian
REAL:: L_sas, L_vk, Delta_cell, Cell_volume, Q_sas1, Q_sas2
REAL:: kx,ky,kz,omx,omy,omz
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV
REAL::MP_PINFl,GAMMAL
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM

I=ICONSIDERED

Verysmall = 10e-16
	
VORTET(1:2,1:2) = ILOCAL_RECON3(I)%GRADS(1:2,1:2)


ux = Vortet(1,1);uy = Vortet(1,2)
vx = Vortet(2,1);vy = Vortet(2,2)


DO IHGT=1,2
  DO IHGJ=1,2
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
END DO

sVORT=0.5*(VORTET+TVORT)
OVORT=0.5*(VORTET-TVORT)





SNORM=SQRT(2.0D0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))))
!Quadratic mean of the strain tensor (defined as Svort). Also needed in SST
ONORM=SQRT(2.0D0*((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+&
	       (OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2))))
OMEGA=ONORM


DIVNORM=ux+uy !Careful with the sign. If it becomes very big, it can produce negative production

usss=sqrt((2.0*((ux*ux)+(vy*vy)))&
	+((uy+vx)*(uy+vx))&
	-(2.0/3.0*(ux+vy)*(ux+vy)))


DERIVTURB(1,1:2) = ILOCAL_RECON3(I)%GRADS(4,1:2)
SQUARET=(sqrt((DERIVTURB(1,1)**2)+(DERIVTURB(1,2)**2)))**2


LEFTV(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
RIGHTV(1:nof_Variables)=LEFTV(1:nof_Variables)
CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)


 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
	      cw1 = cb1 / (kappa*kappa)
	      cw1 = cw1 + (1.0 + cb2) /sigma
	      TCH_X=   TURBMV(1) / VISCL(1) 
	      
			  if (TCH_X .lt. Verysmall) then
						SOURCE_T(1) = ZERO

			  else
			      TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
			      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
			      TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 


			ddw=IELEM(N,I)%WallDist
			
			if (DES_model .eq. 1) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			ddw=min(ddw,C_DES_SA*Delta_cell)
			end if
			
			if (DES_model .eq. 2) then
			CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			Delta_cell=Cell_volume**0.333333333333333
			r_DES=min(10.0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			!Previous limiter is just for numerical reasons regarding tanh
			f_DES=1-tanh((8.0*r_DES)**3)
			ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			end if

			ProdTerm1 = (TURBMV(1))/(leftv(1)* KAPPA * KAPPA * ddw * ddw)
			Stild = max ( OMEGA + (TCH_fv2*ProdTerm1), 0.3*OMEGA)
			Prodtermfinal=Stild*turbmv(1)*cb1/leftv(1)

			RR=MIN((TURBMV(1)/(((LEFTV(1)*STILD*KAPPA * KAPPA * (ddw) * (ddw)))+0.000000001)),10.0)
			cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))
			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) ) /( leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET / SIGMA
			      
			      
			        Prodtermfinal=Stild*cb1
			         
			      destterm  = 2.0* cw1 * fw * (TURBMV(1)/(LEFTV(1)*(kappa**2)*(DDW**2)))
			      fodt  =   -2.0* cb2 * (SQRT(SQUARET))
			      		      
			      
			      
			      
! 			      IF (LAMZ.GT.0)THEN
!  			      SOURCE_T(1)=  min(ProdTermfinal + fodt - destterm,ZERO)
! 			      ELSE
			      SOURCE_T(1)=  min(zero,- destterm +min(fodt,ZERO))
! 			      END IF
			      
			      
			END IF
	 eLSE
		    CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
			  TCH_X=   TURBMV(1) / VISCL(1) 
			    if (tch_x.gt.10.0D0)then
						
			      tch_x=tch_x
						
						
			    ELSE
						
						
				tch_x=0.05D0*log(1.0D0+exp(20.0D0*(tch_x)))
						
						
			      end if
			ddw=IELEM(N,I)%WallDist
			  

			  if (DES_model .eq. 1) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  ddw=min(ddw,C_DES_SA*Delta_cell)
			  end if
			  
			  if (DES_model .eq. 2) then
			  CELL_VOLUME=IELEM(N,I)%TOTVOLUME
			  Delta_cell=Cell_volume**0.333333333333333
			  r_DES=min(10.0D0, (viscl(1)+viscl(3))/(SNORM*(KAPPA*ddw)**2+1e-16)) 
			  !Previous limiter is just for numerical reasons regarding tanh
			  f_DES=1-tanh((8.0D0*r_DES)**3)
			  ddw=max(ddw-f_DES*max(ddw-C_DES_SA*Delta_cell,1e-16),10.0e-16)
			  end if

		      


		  TCH_X3 = (TCH_X)*(TCH_X)*(TCH_X)
		      TCH_FV1 = TCH_X3/(TCH_X3+(CV1*CV1*CV1))
		  TCH_fv2   = 1.0D0 - (TCH_X/(1.0D0 + TCH_X*TCH_fv1)) 
		  ProdTerm1 = (TCH_fv2*tch_x*(TURBMV(1)))/( leftv(1)*KAPPA * KAPPA * (ddw) * (ddw))
		    
		  IF (PRODTERM1.GE.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ProdTerm1
		  END IF
		  IF (PRODTERM1.LT.(-0.7D0*OMEGA))THEN
		  Stild =  OMEGA + ((OMEGA*(((0.7D0*0.7D0)*(OMEGA))+(0.9*PRODTERM1)))/(((0.9-1.4)*OMEGA)-PRODTERM1))
		  END IF



		  Prodtermfinal=Stild*turbmv(1)*cb1*tch_x
		  ! 
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*Stild*KAPPA * KAPPA * (ddw) * (ddw))))

			cw2=0.21+(1.5/((1.0+(TCH_X/40))**2))
		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))/(ddw))**2)
		  ! ! 				!  First order diffusion term
		  fodt  =   LEFTV(1)*cb2 * SQUARET / SIGMA
		  
			      			      		      
			   Prodtermfinal=Stild*cb1*tch_x
			         
			      destterm  = 2.0* cw1 * fw * (TURBMV(1)/(LEFTV(1)*DDW**2))
			      fodt  =   2.0* cb2 * (SQRT(SQUARET)) / SIGMA
			      		      
			      
			      SOURCE_T(1)=  min(ProdTermfinal + fodt - destterm,ZERO)
		  
	 
	 END IF


  CASE(2)		!K OMEGA SST

			   EDDYFL(1)=IELEM(N,I)%WALLDIST
			      EDDYFL(2)=U_CT(I)%VAL(1,1)
			      EDDYFL(3)=U_CT(I)%VAL(1,2)
			      EDDYFL(4:5)=ILOCAL_RECON3(I)%GRADS(1,1:2)
			      EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADS(2,1:2)
			      EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADS(4,1:2)
			      EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADS(5,1:2)
								    

			      eddyfr=eddyfl

			      CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
      k_0=(MAX(Verysmall,U_CT(I)%VAL(1,1)/LEFTV(1))) !First subindex makes reference to the time-stepping
      om_0=MAX(1.0e-1*ufreestream/CharLength,U_CT(I)%VAL(1,2)/LEFTV(1))
!       wally=IELEM(N,I)%WallDist
! 	
! 
! 	      !Calculate here k and omega gradients
! 		  OMX=ILOCAL_RECON3(I)%GRADS(5,1)
! 		  OMY=ILOCAL_RECON3(I)%GRADS(5,2)
! 		  OMZ=ILOCAL_RECON3(I)%GRADS(5,3)
! 		  KX=ILOCAL_RECON3(I)%GRADS(4,1)
! 		  KY=ILOCAL_RECON3(I)%GRADS(4,2)
! 		  KZ=ILOCAL_RECON3(I)%GRADS(4,3)
! 
! 		  dervk_dervom= kx*omx+ky*omy+kz*omz
! 
! 			!Parameters
! 
! 			!-----NEEDED IN DIFFUSION--------
! 			D_omplus=max(2*LEFTV(1)/sigma_om2/om_0*dervk_dervom, 1e-10)
! 			phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(1)/(LEFTV(1)*wally*wally*om_0))
! 			phi_1=min(phi_2, 4.0*LEFTV(1)*k_0/(sigma_om2*D_omplus*wally*wally)) 
! 
! 				F_1=tanh(phi_1**4)
! 				F_2=tanh(phi_2**2)
! 				!--------------------------------
! 
! 
! 				alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
! 				alpha_star=alpha_starinf
! 				alpha_raw=alpha_inf
! 
! 
! 				beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
! 				alpha_star0=beta_i/3.0
! 				beta_stari=beta_starinf
! 
! 
! 
! 
! 				LOWRE=0
! 				!LOW-RE correction------------------------------------------------------------------
! 
! 				if (lowre.eq.1) then
! 				Re_t_SST=LEFTV(1)*k_0/(VISCL(1)*om_0)  !Limiters for this???
! 
! 				alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)
! 				alpha_raw=alpha_inf/alpha_star*(alpha_0+Re_t_SST/R_om_SST)/(1.0+Re_t_SST/R_om_SST)
! 
! 				beta_stari=beta_starinf*(4.0/15.0+(Re_t_SST/R_beta)**4)/(1.0+(Re_t_SST/R_beta)**4)
! 
! 				end if
! 				      !--------------------------------------------------------------------------
! 
! 				      !No Mach number corrections for the beta
! 				      beta_star=beta_stari
! 				      beta_raw=beta_i
! 
! 
! 				      !PRODUCTION TERMS
! 					      !Production of k
! 					      if (vort_model.eq.1) then
! 						      Prod_k=VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM
! 						      Prod_om=alpha_raw*LEFTV(1)*SNORM*ONORM
! 					      else
! 						      !Prod_k=VISCL(3)*SNORM**2
! 						      !Exact formulation: 
! 						      Prod_k=VISCL(3)*(SNORM**2)!-2.0/3.0*DIVNORM*DIVNORM)-2.0/3.0*LEFTV(1)*k_0*DIVNORM
! 						      Prod_k=min(Prod_k,10.0*LEFTV(1)*beta_star*k_0*om_0)
! 
! 						      !INTELLIGENT WAY OF LIMITING: 
! 						      !Prod_k=min(Prod_k,max(10.0*LEFTV(1)*beta_star*k_0*om_0,&
! 						      !	    VISCL(3)*ONORM*SNORM!-2.0/3.0*LEFTV(1)*k_0*DIVNORM))
! 					      
! 					      
! 						      !Production of omega  (Menter does this before correcting Prod_k, 
! 						      !but in  article of 2003 he applies the correction to both)	
! 						      Prod_om=alpha_raw*LEFTV(1)*SNORM**2
! 				      ! 		Prod_om=min(Prod_om,10.0*LEFTV(1)*beta_star*om_0*om_0)
! 					      end if
! 
! 				      !DESTRUCTION TERMS
! 					      !Destruction of k
! 					      Ydest_k=LEFTV(1)*beta_star*k_0*om_0
! 					      !Destruction of omega
! 					      Ydest_om=LEFTV(1)*beta_raw*om_0*om_0
! 					      
! 				      !CROSSED-DIFFUSION TERM
! 					      !Crossed diffusion of omega
! 					      Diff_om=2.0*(1-F_1)*LEFTV(1)/(om_0*sigma_om2)*dervk_dervom
! 	 
! 	
! 				!QSAS TERM: SCALE ADAPTIVE
! 					if (QSAS_model.eq.1) then  !<------------!!!!!!!!!!!!!!!!!!!
! 						!Calculate here second derivative of u 
! 						!Declare all variables
! 						      DO IEX=1,3
! 							VORTET(IEX,1:3)=ILOCAL_RECON3(I)%GRADS(3+TURBULENCEEQUATIONS+IEX,1:3)
! 							
! 						      END DO
! 						
! 						uxx=VORTET(1,1) ; uyy=VORTET(1,2); uzz=VORTET(1,3)				
! 						vxx=VORTET(2,1) ; vyy=VORTET(2,2); vzz=VORTET(2,3)	
! 						wxx=VORTET(3,1) ; wyy=VORTET(3,2); wzz=VORTET(3,3)	
! 						!Compute here cell volume
! 						CELL_VOLUME=IELEM(N,I)%TOTVOLUME
! 						
! 						u_lapl=sqrt((uxx+uyy+uzz)**2+(vxx+vyy+vzz)**2+(wxx+wyy+wzz)**2)
! 						dervk2=kx*kx+ky*ky+kz*kz
! 						dervom2=omx*omx+omy*omy+omz*omz
! 						
! 						Delta_cell=Cell_volume**0.333333333333333
! 						
! 						L_sas=sqrt(k_0)/(beta_star**0.25*om_0)
! 						L_vk=max(kappa*SNORM/u_lapl, &       !This switch provides high wave-number damping
! 							C_smg*Delta_cell*sqrt(kappa*eta2_SAS/(beta_raw/beta_star-alpha_raw)))
! 						
! 						Q_sas1=LEFTV(1)*eta2_SAS*kappa*SNORM**2*(L_sas/L_vk)**2
! 						Q_sas2= -C_SAS*2*LEFTV(1)*k_0/sigma_phi*max(dervk2/k_0**2,dervom2/om_0**2)
! 
! 						Q_SAS=max(Q_sas1+Q_sas2,0.0)
! 					else
! 					Q_SAS=ZERO
! 					end if
! 
! 				    !FINAL SOURCE TERMS
! 
! 
! 					    !DES-SST MODEL (if QSAS_model=2)
! 					    if (QSAS_model .eq.2) then
! 					    CELL_VOLUME=IELEM(N,I)%TOTVOLUME
! 					    Delta_cell=Cell_volume**0.333333333333333
! 					    L_t_DES=sqrt(k_0)/(beta_star*om_0)
! 					    F_DES_SST=max(1.0, L_t_DES/(C_DES_SST*Delta_cell)*(1-F_2))
! 					    !The (1-F_2) is meant to protect the boundary layer. Will result in same 
! 					    !separation point that standard S-A
! 					    Ydest_k=F_DES_SST*Ydest_k
! 					    end if
! 					    
! 				    srcfull_k=Prod_k-Ydest_k
! 				    srcfull_om=Prod_om-Ydest_om+Diff_om+Q_SAS


				    !Filling the output vector
				      SOURCE_T(1) = BETA_STARINF*om_0
				      SOURCE_T(2) = BETA_STARINF*k_0

				  
	
END SELECT

END SUBROUTINE SOURCES_DERIVATIVES2d

end module source
