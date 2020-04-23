MODULE RIEMANN
USE DECLARATION
USE LIBRARY
USE FLOW_OPERATIONS



IMPLICIT NONE

 CONTAINS

SUBROUTINE EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)
!> @brief
!> Subroutine for linear advection IVP
	IMPLICIT NONE
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	INTEGER,INTENT(IN)::N
	REAL,INTENT(IN)::NORMALVECT
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)::CLEFT
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(IN)::CRIGHT
			HLLCFLUX(1)=ZERO
			IF (NORMALVECT.GT.ZERO)THEN
				HLLCFLUX(1)=(NORMALVECT)*CLEFT(1)
			ELSE
				HLLCFLUX(1)=(NORMALVECT)*CRIGHT(1)
			END IF
END SUBROUTINE EXACT_RIEMANN_SOLVER






Subroutine HLLC_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
!> @brief
!> HLLC Riemann solver in 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT_ROT,ROTVL
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT_ROT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,SPED
	REAL::MUL,MUR,LASTL,LASTR
	
	
	
	      
	HLLCFLUX=ZERO
	ROTVL=ZERO
	ROTVR=ZERO
		TEMPFL=ZERO
		TEMPFR=ZERO
		FLSTAR=ZERO
		FRSTAR=ZERO
		ULSTAR=ZERO
		URSTAR=ZERO
		TEMPUL=ZERO
		TEMPUR=ZERO
		FL=ZERO
		FR=ZERO
		!CONSERVATIVE VARIABLES TO PRIMITIVE
		LEFTV(1:5)=CLEFT_ROT(1:5)
		RIGHTV(1:5)=CRIGHT_ROT(1:5)
		
		CALL CONS2PRIM2(N)
		
		ROTVL(1:5)=LEFTV(1:5)
		ROTVR(1:5)=RIGHTV(1:5)
		
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
		
		ROTVL(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=CLEFT_ROT(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
		ROTVR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=CRIGHT_ROT(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		END IF
		
		CALL ESTIMATE_WAVES(N,ROTVL,ROTVR,SL,SM,SR,GAMMA)
		
		
		
		
			!NOW CONDITIONS BASED ON WAVE SPEEDS!
			RL=ROTVL(1);UL=ROTVL(2);VL=ROTVL(3);WL=ROTVL(4);PL=ROTVL(5);EL=CLEFT_ROT(5)
			RR=ROTVR(1);UR=ROTVR(2);VR=ROTVR(3);WR=ROTVR(4);PR=ROTVR(5);ER=CRIGHT_ROT(5)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			
			RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVL(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
			RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
			

! 			IF (TURBULENCEMODEL.EQ.2)THEN
! 			PL=PL+((2.0D0/3.0D0)*EDDYFL(2))
! 			PR=PR+((2.0D0/3.0D0)*EDDYFR(2))  
! 
! 			END IF
			END IF


			FL(1)=RL*UL
			FL(2)=(RL*(UL**2))+PL
			FL(3)=RL*UL*VL
			FL(4)=RL*UL*WL
			FL(5)=UL*(EL+PL)
			
			
			FR(1)=RR*UR
			FR(2)=(RR*(UR**2))+PR
			FR(3)=RR*UR*VR
			FR(4)=RR*UR*WR
			FR(5)=UR*(ER+PR)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			FL(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UL
			FR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UR
			END IF
			
			MUL=RL*((SL(1)-UL)/(SL(1)-SM(1)))
			MUR=RR*((SR(1)-UR)/(SR(1)-SM(1)))
			LASTL=(EL/RL)+((SM(1)-UL)*(SM(1)+((PL)/(RL*(SL(1)-UL)))))
			LASTR=(ER/RR)+((SM(1)-UR)*(SM(1)+((PR)/(RR*(SR(1)-UR)))))
			ULSTAR(1)=MUL
			ULSTAR(2)=MUL*SM(1)
			ULSTAR(3)=MUL*VL
			ULSTAR(4)=MUL*WL
			ULSTAR(5)=MUL*LASTL

			
			URSTAR(1)=MUR
			URSTAR(2)=MUR*SM(1)
			URSTAR(3)=MUR*VR
			URSTAR(4)=MUR*WR
			URSTAR(5)=MUR*LASTR
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			ULSTAR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=MUL*RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)/rl
			URSTAR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=MUR*RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)/rr
			END IF
			
			FLSTAR(:)=FL(:)+SL(1)*(ULSTAR(:)-CLEFT_ROT(:))
			FRSTAR(:)=FR(:)+SR(1)*(URSTAR(:)-CRIGHT_ROT(:))
			
			IF (SL(1).GE.ZERO)THEN
				HLLCFLUX(:)=FL(:)
			END IF
			IF (SR(1).LE.ZERO)THEN
				HLLCFLUX(:)=FR(:)
			END IF
			IF ((SL(1).LE.ZERO).AND.(SM(1).GE.ZERO))THEN
				HLLCFLUX(:)=FLSTAR(:)
			END IF
			IF ((SR(1).GE.ZERO).AND.(SM(1).LE.ZERO))THEN
				HLLCFLUX(:)=FRSTAR(:)
			END IF

END SUBROUTINE HLLC_RIEMANN_SOLVER


SUBROUTINE ROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
!> @brief
!> ROE Riemann solver in 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT,ROTVL,sl,sr,sm
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,INTENT(IN)::GAMMA
	REAL,DIMENSION(5)::TEMPFL,TEMPFR,FLSTAR,FRSTAR,ULSTAR,URSTAR,TEMPUL,TEMPUR
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,RML,RMR,SPED
	REAL::MUL,MUR,LASTL,LASTR
	Real :: sqrtrhoL,sqrtrhoR,HL,HR,utilde,vtilde,wtilde,htilde,atilde,VelTilde
	
!ORIGINALLY OBTAINED FROM Katate Masatsuka, February 2009. http://www.cfdbooks.com

!Input
 real:: primL(5), primR(5) ! Input: primitive variables
 real:: njk(3)             ! Input: face normal vector

!Output

!Some constants
 real::   one = 1.0d0
 real::   two = 2.0d0
 real::  half = 0.5d0
 real:: fifth = 0.2d0

!Local variables
  real:: eig(4)                         ! Eigenvalues
 real:: rhoL, rhoR           ! Primitive variables.
 real:: qnL, qnR                     ! Normal velocities
 real:: aL, aR             ! Speed of sound, Total enthalpy
 real:: RT,rho,u,v,w,H,a,qn          ! Roe-averages
 real:: drho,dqn,dp,LdU(4)           ! Wave strengths
 real:: du, dv, dw                   ! Velocity differences
 real:: ws(4), R(5,4)                ! Wave speeds and right-eigenvectors
 real:: dws(4)                       ! Width of a parabolic fit for entropy fix
 real:: fL(5), fR(5), diss(5)        ! Fluxes ad dissipation term
 real:: SRp,SLm                        ! Wave speeds for the HLL part
 real:: nx1, ny1, nz1                  ! Vector along which HLL is applied
 real:: nx2, ny2, nz2                  ! Vector along which Roe is applied
 real:: alpha1, alpha2                 ! Projections of the new normals
 real:: abs_dq                         ! Magnitude of the velocity difference
 real:: temp, tempx, tempy, tempz      ! Temporary variables
 

! Face normal vector (unit vector)

                LEFTV(1:5)=CLEFT(1:5)
		RIGHTV(1:5)=CRIGHT(1:5)
		
		CALL CONS2PRIM2(N)
		
		priml(1:5)=LEFTV(1:5)
		primr(1:5)=RIGHTV(1:5)

!Primitive and other variables.

!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL                                        !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity

!Wave Strengths

   drho = rhoR - rhoL !Density difference
     dp =   pR - pL   !Pressure difference
    dqn =  qnR - qnL  !Normal velocity difference

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho                         !Shear wave strength (not really, just a factor)

!Absolute values of the wave Speeds

  ws(1) = abs(qn-a) !Left-moving acoustic wave
  ws(2) = abs(qn)   !Entropy wave
  ws(3) = abs(qn+a) !Right-moving acoustic wave
  ws(4) = abs(qn)   !Shear waves

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Right Eigenvectors
!Note: Two shear wave components are combined into one, so that tangent vectors
!      are not required. And that's why there are only 4 vectors here.
!      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx
  R(3,1) = v - a*ny
  R(4,1) = w - a*nz
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx
  R(3,3) = v + a*ny
  R(4,3) = w + a*nz
  R(5,3) = H + a*qn

! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  R(1,4) = zero
  R(2,4) = du - dqn*nx
  R(3,4) = dv - dqn*ny
  R(4,4) = dw - dqn*nz
  R(5,4) = u*du + v*dv + w*dw - qn*dqn

!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
         + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

    HLLCFLUX(1:5)= half * (fL(1:5) + fR(1:5) - diss(1:5))

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a



END SUBROUTINE ROE_RIEMANN_SOLVER

SUBROUTINE tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
!> @brief
!> tROE Riemann solver in 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT,ROTVL,sl,sr,sm
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,INTENT(IN)::GAMMA
	REAL,DIMENSION(5)::TEMPFL,TEMPFR,FLSTAR,FRSTAR,ULSTAR,URSTAR,TEMPUL,TEMPUR
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,RML,RMR,SPED
	REAL::MUL,MUR,LASTL,LASTR
	Real :: sqrtrhoL,sqrtrhoR,HL,HR,utilde,vtilde,wtilde,htilde,atilde,VelTilde
	
!ORIGINALLY OBTAINED FROM Katate Masatsuka, February 2009. http://www.cfdbooks.com

!Input
 real:: primL(5), primR(5) ! Input: primitive variables
 real:: njk(3)             ! Input: face normal vector

!Output

!Some constants
 real::   one = 1.0d0
 real::   two = 2.0d0
 real::  half = 0.5d0
 real:: fifth = 0.2d0
real:: eig(4)                         ! Eigenvalues
!Local variables
  
 real:: rhoL, rhoR           ! Primitive variables.
 real:: qnL, qnR                     ! Normal velocities
 real:: aL, aR             ! Speed of sound, Total enthalpy
 real:: RT,rho,u,v,w,H,a,qn          ! Roe-averages
 real:: drho,dqn,dp,LdU(5)           ! Wave strengths
 real:: du, dv, dw                   ! Velocity differences
 real:: ws(5), R(5,5)                ! Wave speeds and right-eigenvectors
 real:: dws(5)                       ! Width of a parabolic fit for entropy fix
 real:: fL(5), fR(5), diss(5)        ! Fluxes ad dissipation term
 real:: SRp,SLm                        ! Wave speeds for the HLL part
 real:: nx1, ny1, nz1                  ! Vector along which HLL is applied
 real:: nx2, ny2, nz2                  ! Vector along which Roe is applied
 real:: alpha1, alpha2                 ! Projections of the new normals
 real:: abs_dq                         ! Magnitude of the velocity difference
 real:: temp, tempx, tempy, tempz,lx,ly,lz,mx,my,mz, abs_n_cross_l, ql,qm,qll,qml,qlr,qmr,dql,dqm   ! Temporary variables

! Face normal vector (unit vector)

                LEFTV(1:5)=CLEFT(1:5)
		RIGHTV(1:5)=CRIGHT(1:5)
		
		CALL CONS2PRIM2(N)
		
		priml(1:5)=LEFTV(1:5)
		primr(1:5)=RIGHTV(1:5)

!Primitive and other variables.
   
!  Left state
    tempx = ny*ny + nz*nz
     tempy = nz*nz + nx*nx
     tempz = nx*nx + ny*ny

     if     ( tempx >= tempy .and. tempx >= tempz ) then
       lx =  zero
       ly = -nz
       lz =  ny
     elseif ( tempy >= tempx .and. tempy >= tempz ) then
       lx = -nz
       ly =  zero
       lz =  nx
     elseif ( tempz >= tempx .and. tempz >= tempy ) then
       lx = -ny
       ly =  nx
       lz =  zero
     else
      ! Impossible to happen
      write(*,*) "subroutine inviscid_roe: Impossible to happen. Please report the problem."
      stop
     endif

!     Make it the unit vector.
      temp = sqrt( lx*lx + ly*ly + lz*lz )
       lx = lx/temp
       ly = ly/temp
       lz = lz/temp

! m = (mx,my,mz)
!
! The other one, m = (mx,my,mz), is chosen as a vector orthogonal to both n and l
! defined by the vector product: m = n x l / |n x l|

  mx = ny*lz - nz*ly
  my = nz*lx - nx*lz
  mz = nx*ly - ny*lx

  abs_n_cross_l = sqrt(mx**2 + my**2 + mz**2)
  mx = mx / abs_n_cross_l
  my = my / abs_n_cross_l
  mz = mz / abs_n_cross_l

!(Do you like such ambiguous tangent vectors? Actually, the Roe flux can
! be implemented without any tangent vector. See "I do like CFD, VOL.1",
! or the subroutine "inviscid_roe_n" or "inviscid_rotated_rhll" for details.
! In fact, the resulting flux is independent of the choice of these tangent vectors
! as it should be.)

!Primitive and other variables.

!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
     qlL = uL*lx + vL*ly + wL*lz
     qmL = uL*mx + vL*my + wL*mz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
     qlR = uR*lx + vR*ly + wR*lz
     qmR = uR*mx + vR*my + wR*mz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL                                        !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity
    ql = u*lx + v*ly + w*lz                             !Roe-averaged face-tangent velocity
    qm = u*mx + v*my + w*mz                             !Roe-averaged face-tangent velocity

!Wave Strengths

   drho = rhoR - rhoL !Density difference
     dp =   pR - pL   !Pressure difference
    dqn =  qnR - qnL  !Normal velocity difference
    dql =  qlR - qlL  !Tangent velocity difference in l
    dqm =  qmR - qmL  !Tangent velocity difference in m

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho*dql                     !Shear wave strength
  LdU(5) = rho*dqm                     !Shear wave strength

!Absolute values of the wave speeds

  ws(1) = abs(qn-a) !Left-moving acoustic wave speed
  ws(2) = abs(qn)   !Entropy wave speed
  ws(3) = abs(qn+a) !Right-moving acoustic wave speed
  ws(4) = abs(qn)   !Shear wave speed
  ws(5) = abs(qn)   !Shear wave speed

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Right Eigenvectors

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx
  R(3,1) = v - a*ny
  R(4,1) = w - a*nz
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx
  R(3,3) = v + a*ny
  R(4,3) = w + a*nz
  R(5,3) = H + a*qn

! Shear wave
  R(1,4) = zero
  R(2,4) = lx
  R(3,4) = ly
  R(4,4) = lz
  R(5,4) = ql

! Shear wave
  R(1,5) = zero
  R(2,5) = mx
  R(3,5) = my
  R(4,5) = mz
  R(5,5) = qm

!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) + ws(3)*LdU(3)*R(:,3) &
         + ws(4)*LdU(4)*R(:,4) + ws(5)*LdU(5)*R(:,5)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

!   num_flux = half * (fL + fR - diss)
! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

    HLLCFLUX(1:5)= half * (fL + fR - diss)

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a



END SUBROUTINE tROE_RIEMANN_SOLVER

SUBROUTINE rROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
!> @brief
!> rROE Riemann solver in 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT,ROTVL,sl,sr,sm
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,INTENT(IN)::GAMMA
	REAL,DIMENSION(5)::TEMPFL,TEMPFR,FLSTAR,FRSTAR,ULSTAR,URSTAR,TEMPUL,TEMPUR
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,RML,RMR,SPED
	REAL::MUL,MUR,LASTL,LASTR
	Real :: sqrtrhoL,sqrtrhoR,HL,HR,utilde,vtilde,wtilde,htilde,atilde,VelTilde
	
!ORIGINALLY OBTAINED FROM Katate Masatsuka, February 2009. http://www.cfdbooks.com
!Input
 real:: primL(5), primR(5) ! Input: primitive variables
 real:: njk(3)             ! Input: face normal vector

!Output

!Some constants
 real::   one = 1.0d0
 real::   two = 2.0d0
 real::  half = 0.5d0
 real:: fifth = 0.2d0
real:: eig(4)                         ! Eigenvalues
!Local variables
  
 real:: rhoL, rhoR           ! Primitive variables.
 real:: qnL, qnR                     ! Normal velocities
 real:: aL, aR             ! Speed of sound, Total enthalpy
 real:: RT,rho,u,v,w,H,a,qn          ! Roe-averages
 real:: drho,dqn,dp,LdU(4)           ! Wave strengths
 real:: du, dv, dw                   ! Velocity differences
 real:: ws(4), R(5,4)                ! Wave speeds and right-eigenvectors
 real:: dws(4)                       ! Width of a parabolic fit for entropy fix
 real:: fL(5), fR(5), diss(5)        ! Fluxes ad dissipation term
 real:: SRp,SLm                        ! Wave speeds for the HLL part
 real:: nx1, ny1, nz1                  ! Vector along which HLL is applied
 real:: nx2, ny2, nz2                  ! Vector along which Roe is applied
 real:: alpha1, alpha2                 ! Projections of the new normals
 real:: abs_dq                         ! Magnitude of the velocity difference
 real:: temp, tempx, tempy, tempz      ! Temporary variables

! Face normal vector (unit vector)

                LEFTV(1:5)=CLEFT(1:5)
		RIGHTV(1:5)=CRIGHT(1:5)
		
		CALL CONS2PRIM2(N)
		
		priml(1:5)=LEFTV(1:5)
		primr(1:5)=RIGHTV(1:5)

!Primitive and other variables.
   
!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)

!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

!--------------------------------------------------------------------------------
!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
! Note: n1 and n2 may need to be frozen at some point during 
!       a steady calculation to fully make it converge. For time-accurate 
!       calculation, it is not necessary.
! Note: For a boundary face, you may want to set (nx2,ny2,nz2)=(nx,ny,nz),
!       (nx1,ny1)=tangent vector, i.e., use the Roe flux.

    abs_dq = sqrt( (uR-uL)**2 + (vR-vL)**2 + (wR-wL)**2 )

!n1 = Velocity difference vector: normal to shock or tangent to shear
  if ( abs_dq > tolsmall) then

       nx1 = (uR-uL)/abs_dq
       ny1 = (vR-vL)/abs_dq
       nz1 = (wR-wL)/abs_dq

!n1 = Face tangent vector if abs_dq is too small.
!     Note: There are infinitely many choices for the tangent vector.
!           The best choice may be discovered in future.
! Here, we choose a vector in a plane (essentially 2D vector).
! Note that we must be careful and make sure that the vector is not a zero vector.
  else

!  E.g., if n = (1,0,0), then, (0,-nz,ny) = (0,0,0). This is a problem.
!        To avoid such a situation, we choose the 2D tangential vector
!        having the maximum length.

     tempx = ny*ny + nz*nz
     tempy = nz*nz + nx*nx
     tempz = nx*nx + ny*ny

     if     ( tempx >= tempy .and. tempx >= tempz ) then
       nx1 =  zero
       ny1 = -nz
       nz1 =  ny
     elseif ( tempy >= tempx .and. tempy >= tempz ) then
       nx1 = -nz
       ny1 =  zero
       nz1 =  nx
     elseif ( tempz >= tempx .and. tempz >= tempy ) then
       nx1 = -ny
       ny1 =  nx
       nz1 =  zero
     else
      ! Impossible to happen
      write(*,*) "inviscid_rotated_rhll: Impossible to happen. Please report the problem."
      stop
     endif

!     Make it the unit vector.
      temp = sqrt( nx1*nx1 + ny1*ny1 + nz1*nz1 )
       nx1 = nx1/temp
       ny1 = ny1/temp
       nz1 = nz1/temp

  endif

    alpha1 = nx*nx1 + ny*ny1 + nz*nz1

! Make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
       nz1 = temp * nz1
    alpha1 = temp * alpha1

!n2 = direction perpendicular to n1.
!     Note: There are infinitely many choices for this vector.
!           The best choice may be discovered in future.
! Here, we employ the formula (4.4) in the paper:
!     (nx2,ny2,nz2) = (n1xn)xn1 / |(n1xn)xn1|    ('x' is the vector product.)

!  (tempx,tempy,tempz) = n1xn
     tempx = ny1*nz - nz1*ny
     tempy = nz1*nx - nx1*nz
     tempz = nx1*ny - ny1*nx

!  (nx2,ny2,nz2) = (n1xn)xn1
     nx2 = tempy*nz1 - tempz*ny1
     ny2 = tempz*nx1 - tempx*nz1
     nz2 = tempx*ny1 - tempy*nx1

!  Make n2 the unit vector
     temp = sqrt( nx2*nx2 + ny2*ny2 + nz2*nz2 )
       nx2 = nx2/temp
       ny2 = ny2/temp
       nz2 = nz2/temp

    alpha2 = nx*nx2 + ny*ny2 + nz*nz2

!  Make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
       nz2 = temp * nz2
    alpha2 = temp * alpha2

!--------------------------------------------------------------------------------
!Now we are going to compute the Roe flux with n2 as the normal with modified 
!wave speeds (5.12). NOTE: the Roe flux here is computed without tangent vectors.
!See "I do like CFD, VOL.1" for details: page 57, Equation (3.6.31).

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

      RT = sqrt(rhoR/rhoL)
     rho = RT*rhoL                                        !Roe-averaged density.
       u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
       v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
       w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
       H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
       a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound

!----------------------------------------------------
!Compute the wave speed estimates for the HLL part,
!following Einfeldt:
!
! B. Einfeldt, On Godunov-type methods for gas dynamics,
! SIAM Journal on Numerical Analysis 25 (2) (1988) 294–318.
!
! Note: HLL is actually applied to n1, but this is
!       all we need to incorporate HLL. See JCP2008 paper.

     qn  = u *nx1 + v *ny1 + w *nz1
     qnL = uL*nx1 + vL*ny1 + wL*nz1
     qnR = uR*nx1 + vR*ny1 + wR*nz1
     SLm = min( zero, qn - a, qnL - aL ) !Minimum wave speed estimate
     SRp = max( zero, qn + a, qnR + aR ) !Maximum wave speed estimate

! This is the only place where n1=(nx1,ny1,nz1) is used.
! n1=(nx1,ny1,nz1) is never used below.
!----------------------------------------------------

!Wave Strengths

     qn  = u *nx2 + v *ny2 + w *nz2
     qnL = uL*nx2 + vL*ny2 + wL*nz2
     qnR = uR*nx2 + vR*ny2 + wR*nz2

    drho = rhoR - rhoL  !Density difference
      dp =   pR - pL    !Pressure difference
     dqn =  qnR - qnL   !Normal velocity difference

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho                         !Shear wave strength (not really, just a factor)

!Wave Speed (Eigenvalues)

  eig(1) = qn-a !Left-moving acoustic wave velocity
  eig(2) = qn   !Entropy wave velocity
  eig(3) = qn+a !Right-moving acoustic wave velocity
  eig(4) = qn   !Shear wave velocity

!Absolute values of the wave speeds (Eigenvalues)

   ws(1) = abs(qn-a) !Left-moving acoustic wave speed
   ws(2) = abs(qn)   !Entropy wave speed
   ws(3) = abs(qn+a) !Right-moving acoustic wave speed
   ws(4) = abs(qn)   !Shear wave speed

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Combine the wave speeds for Rotated-RHLL: Eq.(5.12) in the original JCP2008 paper.

      ws = alpha2*ws - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig)/(SRp-SLm)

!Below, we compute the Roe dissipation term in the direction n2
!with the above modified wave speeds. HLL wave speeds act something like
!the entropy fix or eigenvalue limiting; they contribute only by the amount
!given by the fraction, alpha1 (less than or equal to 1.0). See JCP2008 paper.

!Right Eigenvectors:
!Note: Two shear wave components are combined into one, so that tangent vectors
!      are not required. And that's why there are only 4 vectors here.

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx2
  R(3,1) = v - a*ny2
  R(4,1) = w - a*nz2
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx2
  R(3,3) = v + a*ny2
  R(4,3) = w + a*nz2
  R(5,3) = H + a*qn

! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  R(1,4) = zero
  R(2,4) = du - dqn*nx2
  R(3,4) = dv - dqn*ny2
  R(4,4) = dw - dqn*nz2
  R(5,4) = u*du + v*dv + w*dw - qn*dqn

!Dissipation Term: Roe dissipation with the modified wave speeds.
! |An|dU = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ], where n=n2.

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
         + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)

!Compute the Rotated-RHLL flux. (It looks like the HLL flux with Roe dissipation.)

!   num_flux = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss
! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

    HLLCFLUX(1:5)= (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a



END SUBROUTINE rROE_RIEMANN_SOLVER


Subroutine RUSANOV_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
!> @brief
!> Rusanov Riemann solver in 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT_ROT,ROTVL
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT_ROT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,SPED
	REAL::MUL,MUR,LASTL,LASTR
	
	
	
	      
	HLLCFLUX=ZERO
	ROTVL=ZERO
	ROTVR=ZERO
		TEMPFL=ZERO
		TEMPFR=ZERO
		FLSTAR=ZERO
		FRSTAR=ZERO
		ULSTAR=ZERO
		URSTAR=ZERO
		TEMPUL=ZERO
		TEMPUR=ZERO
		FL=ZERO
		FR=ZERO
		!CONSERVATIVE VARIABLES TO PRIMITIVE
		LEFTV(1:5)=CLEFT_ROT(1:5)
		RIGHTV(1:5)=CRIGHT_ROT(1:5)
		
		CALL CONS2PRIM2(N)
		
		ROTVL(1:5)=LEFTV(1:5)
		ROTVR(1:5)=RIGHTV(1:5)
		
		
		
		
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
		
		ROTVL(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=CLEFT_ROT(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
		ROTVR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=CRIGHT_ROT(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		END IF
		
		CALL ESTIMATE_WAVES(N,ROTVL,ROTVR,SL,SM,SR,GAMMA)
			!NOW CONDITIONS BASED ON WAVE SPEEDS!
			RL=ROTVL(1);UL=ROTVL(2);VL=ROTVL(3);WL=ROTVL(4);PL=ROTVL(5);EL=CLEFT_ROT(5)
			RR=ROTVR(1);UR=ROTVR(2);VR=ROTVR(3);WR=ROTVR(4);PR=ROTVR(5);ER=CRIGHT_ROT(5)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			
			RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVL(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
			RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
			

! 			IF (TURBULENCEMODEL.EQ.2)THEN
! 			PL=PL+((2.0D0/3.0D0)*EDDYFL(2))
! 			PR=PR+((2.0D0/3.0D0)*EDDYFR(2))  
! 
! 			END IF
			END IF


			FL(1)=RL*UL
			FL(2)=(RL*(UL**2))+PL
			FL(3)=RL*UL*VL
			FL(4)=RL*UL*WL
			FL(5)=UL*(EL+PL)
			
			
			FR(1)=RR*UR
			FR(2)=(RR*(UR**2))+PR
			FR(3)=RR*UR*VR
			FR(4)=RR*UR*WR
			FR(5)=UR*(ER+PR)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			FL(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UL
			FR(6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UR
			END IF
			
			
			
			

			
			
			
			
			
			
			HLLCFLUX(:)=0.5d0*(FL(:)+FR(:))-0.5d0*MAX(ABS(SL(1)),ABS(SR(1)))*(Cright_ROT(:)-Cleft_ROT(:))
			
			

END SUBROUTINE RUSANOV_RIEMANN_SOLVER




SUBROUTINE ESTIMATE_WAVES(N,ROTVL,ROTVR,SL,SM,SR,GAMMA)
!> @brief
!> Wave speed estimates for 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,K
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVL
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	REAL::CL,CR,PR,PL,UL,UR,VL,VR,WR,WL,RL,RR
	REAL::CUP,PPV,PMIN,PMAX,QMAX,QUSER,BL,BR,COV,PM,UM
	REAL::G1,G2,G3,G4,G5,G6,G7,G8,GEL,GER,PQ,PTL,PTR
	G1 = (GAMMA - 1.0d0)/(2.0d0*GAMMA)
    	G2 = (GAMMA + 1.0d0)/(2.0d0*GAMMA)
   	G3 = 2.0d0*GAMMA/(GAMMA - 1.0d0)
    	G4 = 2.0d0/(GAMMA - 1.0d0)
    	G5 = 2.0d0/(GAMMA + 1.0d0)
    	G6 = (GAMMA - 1.0d0)/(GAMMA + 1.0d0)
   	G7 = (GAMMA - 1.0d0)/2.0d0
   	G8 = GAMMA - 1.0d0
	
	SL=0.0d0
	SR=0.0d0
	SM=0.0d0
	COV=0.0d0
	!BUILD LEFT STATE VARIABLES
	RL=ROTVL(1)
	UL=ROTVL(2)
	VL=ROTVL(3)
	WL=ROTVL(4)
	PL=ROTVL(5)
	CL=SQRT((PL*GAMMA)/(RL))
	!BUILD RIGHT  STATE VARIABLES
	RR=ROTVR(1)
	UR=ROTVR(2)
	VR=ROTVR(3)
	WR=ROTVR(4)
	PR=ROTVR(5)
	CR=SQRT((PR*GAMMA)/(RR))

	CUP=0.25d0*(RL+RR)*(CL+CR)
	PPV=0.5d0*(PL + PR) + 0.5d0*(UL - UR)*CUP
	PPV=MAX(0.0d0,PPV)
	PMIN=MIN(PL,PR)
	PMAX=MAX(PL,PR)
	QMAX=PMAX/PMIN
	QUSER=2.0d0

	 IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
  
!        Select PRVS Riemann solver
 
         PM = PPV
         UM = 0.5d0*(UL + UR) + 0.5d0*(PL - PR)/CUP 
  
      	ELSE
 
         BL = 1.0d0 - COV*RL
         BR = 1.0d0 - COV*RR
 
         IF(PPV.LT.PMIN)THEN
 
!           Select Two-Rarefaction Riemann solver
 
        
            PQ  = (PL/PR)**G1
            UM  = (PQ*UL/CL/BL + UR/CR/BR + G4*(PQ - 1.0d0)) 
            UM  = UM/(PQ/CL/BL + 1.0d0/CR/BR)
            PTL = 1.0d0 + G7*(UL - UM)/CL/BL
            PTR = 1.0d0 + G7*(UM - UR)/CR/BR
            PM  = 0.5d0*(PL*PTL**G3 + PR*PTR**G3)
         ELSE

!           Use Two-Shock Riemann solver with PVRS as estimate
 
!           introduce iterations with PVRS as initial guess
 
            DO K=1,4
             GEL = SQRT((G5*BL/RL)/(G6*PL + PPV))
             GER = SQRT((G5*BR/RR)/(G6*PR + PPV))
             PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
             UM  = 0.5d0*(UL + UR) + 0.5d0*(GER*(PM - PR) - GEL*(PM - PL))	     
             IF ( ABS((PM-PPV)/PM) .LE. 1D-8) GOTO 101
                PPV = PM
     	    END DO
         ENDIF
      ENDIF
      
      101 continue

!     Find speeds
 
      IF(PM.LE.PL)THEN
         SL(1) = UL - CL
      ELSE
         SL(1) = UL - CL*SQRT(1.0d0 + G2*(PM/PL - 1.0d0))
      ENDIF
 
      SM(1)= UM
 
      IF(PM.LE.PR)THEN
         SR(1)= UR + CR
      ELSE
         SR(1) = UR + CR*SQRT(1.0d0 + G2*(PM/PR - 1.0d0))
      ENDIF

END SUBROUTINE ESTIMATE_WAVES



Subroutine HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
!> @brief
!> HLLC Riemann solver in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k,riem_adapt
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT_ROT,ROTVL
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT_ROT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,SPED
	REAL::MUL,MUR,LASTL,LASTR,pgrad,om_P
	
	
	
	      
	HLLCFLUX=ZERO
	ROTVL=ZERO
	ROTVR=ZERO
		TEMPFL=ZERO
		TEMPFR=ZERO
		FLSTAR=ZERO
		FRSTAR=ZERO
		ULSTAR=ZERO
		URSTAR=ZERO
		TEMPUL=ZERO
		TEMPUR=ZERO
		FL=ZERO
		FR=ZERO
		!CONSERVATIVE VARIABLES TO PRIMITIVE
		LEFTV(1:4)=CLEFT_ROT(1:4)
		RIGHTV(1:4)=CRIGHT_ROT(1:4)
		
		CALL CONS2PRIM2d2(N)
		
		ROTVL(1:4)=LEFTV(1:4)
		ROTVR(1:4)=RIGHTV(1:4)
		
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
		
		ROTVL(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=CLEFT_ROT(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
		ROTVR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=CRIGHT_ROT(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		END IF
		
! 		CALL ESTIMATE_WAVES2d(N,ROTVL,ROTVR,SL,SM,SR,GAMMA)
		
		
		
		
		
		
		
		
		
		
		
			!NOW CONDITIONS BASED ON WAVE SPEEDS!
			RL=ROTVL(1);UL=ROTVL(2);VL=ROTVL(3);PL=ROTVL(4);EL=CLEFT_ROT(4)
			RR=ROTVR(1);UR=ROTVR(2);VR=ROTVR(3);PR=ROTVR(4);ER=CRIGHT_ROT(4)
			
			
			
			
			sl(1)=min((ul)-sqrt(gamma*pl/rl),(ur)-sqrt(gamma*pr/rr))
		sr(1)=max((ul)+sqrt(gamma*pl/rl),(ur)+sqrt(gamma*pr/rr))	
		
		sm(1)=(pr-pl+(rl*ul*(sl(1)-ul))-(rr*ur*(sr(1)-ur)))/((rl*(sl(1)-ul))-(rr*(sr(1)-ur)))
			
			
			
			
			
			
			
			
			
			
			
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			
			RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVL(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
			RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
			

! 			IF (TURBULENCEMODEL.EQ.2)THEN
! 			PL=PL+((2.0D0/3.0D0)*EDDYFL(2))
! 			PR=PR+((2.0D0/3.0D0)*EDDYFR(2))  
! 
! 			END IF
			END IF


			FL(1)=RL*UL
			FL(2)=(RL*(UL**2))+PL
			FL(3)=RL*UL*VL
			
			FL(4)=UL*(EL+PL)
			
			
			FR(1)=RR*UR
			FR(2)=(RR*(UR**2))+PR
			FR(3)=RR*UR*VR
			
			FR(4)=UR*(ER+PR)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			FL(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UL
			FR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UR
			END IF
			
			MUL=RL*((SL(1)-UL)/(SL(1)-SM(1)))
			MUR=RR*((SR(1)-UR)/(SR(1)-SM(1)))
			LASTL=(EL/RL)+((SM(1)-UL)*(SM(1)+((PL)/(RL*(SL(1)-UL)))))
			LASTR=(ER/RR)+((SM(1)-UR)*(SM(1)+((PR)/(RR*(SR(1)-UR)))))
			ULSTAR(1)=MUL
			ULSTAR(2)=MUL*SM(1)
			ULSTAR(3)=MUL*VL
			
			ULSTAR(4)=MUL*LASTL

			
			URSTAR(1)=MUR
			URSTAR(2)=MUR*SM(1)
			URSTAR(3)=MUR*VR
			
			URSTAR(4)=MUR*LASTR
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			ULSTAR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=MUL*RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)/rl
			URSTAR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=MUR*RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)/rr
			END IF
			
			FLSTAR(:)=FL(:)+SL(1)*(ULSTAR(:)-CLEFT_ROT(:))
			FRSTAR(:)=FR(:)+SR(1)*(URSTAR(:)-CRIGHT_ROT(:))
			
			IF (SL(1).GE.ZERO)THEN
				HLLCFLUX(:)=FL(:)
			END IF
			IF (SR(1).LE.ZERO)THEN
				HLLCFLUX(:)=FR(:)
			END IF
			IF ((SL(1).LE.ZERO).AND.(SM(1).GE.ZERO))THEN
				HLLCFLUX(:)=FLSTAR(:)
			END IF
			IF ((SR(1).GE.ZERO).AND.(SM(1).LE.ZERO))THEN
				HLLCFLUX(:)=FRSTAR(:)
			END IF
			
			
			
			!pgrad=abs(pl-pr)/min(pl,pr)
			!om_p=0.5d0-0.5d0*sign(pgrad-0.2,1.0d0)*(1.0-exp(-100.0d0*abs(pgrad-0.2)))
			
			!if(om_p.lt.0.9d0)then
			
			!else
			
			!sl(1)=abs(ul)+sqrt(gamma*pl/rl)
			!sr(1)=abs(ur)+sqrt(gamma*pr/rr)
			!HLLCFLUX(:)=0.5d0*(FL(:)+FR(:))-0.5d0*MAX(ABS(SL(1)),ABS(SR(1)))*(Cright_ROT(:)-Cleft_ROT(:))
			!end if
			
			
			
			
			
			
			

END SUBROUTINE HLLC_RIEMANN_SOLVER2d



SUBROUTINE ROE_RIEMANN_SOLVER2d(N,Cleft,Cright,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
!> @brief
!> Roe Riemann solver in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVL,CLEFT,CRIGHT
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	real :: uL(4), uR(4) !  Input: conservative variables rho*[1, u, v, E]
 
 real :: Roe(4)       ! Output: Roe flux function (upwind)
!Local constants
                      ! Ratio of specific heat.
 real ::  fifth, half, one, two    ! Numbers
!Local variables
 real :: tx, ty       ! Tangent vector (perpendicular to the face normal)
 real :: vxL, vxR, vyL, vyR             ! Velocity components.
 real :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real :: vnL, vnR, vtL, vtR             ! Normal and tangent velocities
 real :: aL, aR, HL, HR                 ! Speeds of sound.
 real :: RT,rho,vx,vy,H,a,vn, vt        ! Roe-averages
 real :: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 real :: ws(4),dws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 integer ::  j
!ORIGINALLY OBTAINED FROM Katate Masatsuka, February 2009. http://www.cfdbooks.com
!Constants.

      
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0

       ul(1:4)=cleft(1:4)
       ur(1:4)=cright(1:4)
       
       
!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  tx = -ny
  ty = nx

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
     vnL = vxL*nx+vyL*ny
     vtL = vxL*tx+vyL*ty
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
     vnR = vxR*nx+vyR*ny
     vtR = vxR*tx+vyR*ty
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
    vx = (vxL+RT*vxR)/(one+RT)
    vy = (vyL+RT*vyR)/(one+RT)
     H = ( HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty

!Wave Strengths
   drho = rhoR - rhoL 
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) = rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx
  Rv(3,1) = vy - a*ny
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = vt*a

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx
  Rv(3,4) = vy + a*ny
  Rv(4,4) =  H + vn*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*vnL
  fL(2) = rhoL*vnL * vxL + pL*nx
  fL(3) = rhoL*vnL * vyL + pL*ny
  fL(4) = rhoL*vnL *  HL

  fR(1) = rhoR*vnR
  fR(2) = rhoR*vnR * vxR + pR*nx
  fR(3) = rhoR*vnR * vyR + pR*ny
  fR(4) = rhoR*vnR *  HR

  hllcflux(1:4)= half * (fL(1:4) + fR(1:4) - diss(1:4))





END SUBROUTINE ROE_RIEMANN_SOLVER2d


SUBROUTINE RROE_RIEMANN_SOLVER2d(N,Cleft,Cright,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
!> @brief
!> rROE Riemann solver in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVL,CLEFT,CRIGHT
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
real :: uL(4), uR(4)    !  Input: conservative variables rho*[1, u, v, E]
 
 real :: Rotated_RHLL(4) ! Output: Rotated_RHLL flux function.
!Local constants
 
 real :: fifth, half, one, two    ! Numbers
 real :: eps                            ! 
!Local variables
 real :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
 real :: tx, ty                         ! Tangent vector (taken as n1)
 real :: alpha1, alpha2                 ! Projections of the new normals
 real :: vxL, vxR, vyL, vyR             ! Velocity components.
 real :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real :: vnL, vnR, vtL, vtR             ! Normal and tagent velocities
 real :: aL, aR, HL, HR                 ! Speeds of sound and total enthalpy
 real :: RT,rho,vx,vy,H,a               ! Roe-averages
 real :: vn, vt                         ! Normal and tagent velocities(Roe-average)
 real :: drho,dvx,dvy,dvn,dvt,dp,dV(4)  ! Wave strenghs
 real :: abs_dq                         ! Magnitude of the velocity difference
 real :: abs_ws(4),ws(4),dws(4), Rv(4,4)! Wave speeds and right-eigevectors
 real :: SRp,SLm                        ! Wave speeds for the HLL part
 real :: fL(4), fR(4), diss(4)          ! Fluxes ad dissipation term
 real :: temp
 integer ::  j

!Constants.
!ORIGINALLY OBTAINED FROM Katate Masatsuka, February 2009. http://www.cfdbooks.com
    
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0
       eps = tolsmall ! 1.0e-12 in the original paper (double precision)
       
       UL(1:4)=CLEFT(1:4)
       UR(1:4)=CRIGHT(1:4)

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
     vxL = uL(2)/uL(1)
     vyL = uL(3)/uL(1)
      pL = (gamma-one)*( uL(4) - half*rhoL*(vxL*vxL+vyL*vyL) )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(4) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
     vxR = uR(2)/uR(1)
     vyR = uR(3)/uR(1)
      pR = (gamma-one)*( uR(4) - half*rhoR*(vxR*vxR+vyR*vyR) )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(4) + pR ) / rhoR

     vnL = vxL*nx + vyL*ny
     vnR = vxR*nx + vyR*ny

!Compute the flux.
   fL(1) = rhoL*vnL
   fL(2) = rhoL*vnL * vxL + pL*nx
   fL(3) = rhoL*vnL * vyL + pL*ny
   fL(4) = rhoL*vnL *  HL

   fR(1) = rhoR*vnR
   fR(2) = rhoR*vnR * vxR + pR*nx
   fR(3) = rhoR*vnR * vyR + pR*ny
   fR(4) = rhoR*vnR *  HR

!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. For time-accurate 
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxR-vxL)**2+(vyR-vyL)**2 )
  if ( abs_dq > eps) then
       nx1 = (vxR-vxL)/abs_dq
       ny1 = (vyR-vyL)/abs_dq
  else
    nx1 = -ny 
    ny1 =  nx
  endif
  
  
  
  
  
    alpha1 = nx * nx1 + ny * ny1 
!   To make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
    alpha1 = temp * alpha1

! Take n2 as perpendicular to n1.
       nx2 = -ny1
       ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
!   To make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
    alpha2 = temp * alpha2
    
    
    if(b_code.gt.0)then
  nx2=nx
  ny2=ny
  nx1=-ny
  ny1=nx
  end if

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tagent vector, with modified wave speeds (5.12)

!Compute the Roe Averages
     RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
     vx = (vxL+RT*vxR)/(one+RT)
     vy = (vyL+RT*vyR)/(one+RT)
      H = ( HL+RT* HR)/(one+RT)
      a = sqrt( (gamma-one)*(H-half*(vx*vx+vy*vy)) )
     vn = vx*nx2+vy*ny2
     vt = vx*nx1+vy*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    vnL = vxL*nx2 + vyL*ny2
    vnR = vxR*nx2 + vyR*ny2
    vtL = vxL*nx1 + vyL*ny1
    vtR = vxR*nx1 + vyR*ny1

   drho = rhoR - rhoL 
     dp =   pR - pL
    dvn =  vnR - vnL
    dvt =  vtR - vtL

  dV(1) = (dp - rho*a*dvn )/(two*a*a)
  dV(2) =  rho*dvt/a
  dV(3) =  drho - dp/(a*a)
  dV(4) = (dp + rho*a*dvn )/(two*a*a)

!Wave Speeds for Roe flux part.
    ws(1) = vn-a
    ws(2) = vn
    ws(3) = vn
    ws(4) = vn+a
  abs_ws  = abs(ws)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
!only for the nonlinear fields.
  dws(1) = fifth
   if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  dws(4) = fifth
   if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

!HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( zero, vtR + aR, vt + a)
   SLm = min( zero, vtL - aL, vt - a)

!Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)

!Right Eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1
  ty = ny1

  Rv(1,1) = one    
  Rv(2,1) = vx - a*nx2
  Rv(3,1) = vy - a*ny2
  Rv(4,1) =  H - vn*a

  Rv(1,2) = zero
  Rv(2,2) = a*tx
  Rv(3,2) = a*ty
  Rv(4,2) = a*vt

  Rv(1,3) = one
  Rv(2,3) = vx
  Rv(3,3) = vy 
  Rv(4,3) = half*(vx*vx+vy*vy)

  Rv(1,4) = one
  Rv(2,4) = vx + a*nx2
  Rv(3,4) = vy + a*ny2
  Rv(4,4) =  H + vn*a

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dV(j)*Rv(i,j)
   end do
  end do

!Compute the Rotated-RHLL flux.
  hllcflux = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

 





END SUBROUTINE RROE_RIEMANN_SOLVER2d





Subroutine RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
!> @brief
!> Rusanov Riemann solver in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,k
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CLEFT_ROT,ROTVL
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::CRIGHT_ROT,ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::HLLCFLUX
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	REAL::RL,RR,PL,PR,EL,ER,UL,UR,VL,VR,WL,WR,SPED
	REAL::MUL,MUR,LASTL,LASTR
	
	
	
	      
	HLLCFLUX=ZERO
	ROTVL=ZERO
	ROTVR=ZERO
		TEMPFL=ZERO
		TEMPFR=ZERO
		FLSTAR=ZERO
		FRSTAR=ZERO
		ULSTAR=ZERO
		URSTAR=ZERO
		TEMPUL=ZERO
		TEMPUR=ZERO
		FL=ZERO
		FR=ZERO
		!CONSERVATIVE VARIABLES TO PRIMITIVE
		LEFTV(1:4)=CLEFT_ROT(1:4)
		RIGHTV(1:4)=CRIGHT_ROT(1:4)
		
		CALL CONS2PRIM2d2(N)
		
		ROTVL(1:4)=LEFTV(1:4)
		ROTVR(1:4)=RIGHTV(1:4)
		
		
		
		
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
		
		ROTVL(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=CLEFT_ROT(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
		ROTVR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=CRIGHT_ROT(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
		
		END IF
		
		CALL ESTIMATE_WAVES2D(N,ROTVL,ROTVR,SL,SM,SR,GAMMA)
			!NOW CONDITIONS BASED ON WAVE SPEEDS!
			RL=ROTVL(1);UL=ROTVL(2);VL=ROTVL(3);PL=ROTVL(4);EL=CLEFT_ROT(4)
			RR=ROTVR(1);UR=ROTVR(2);VR=ROTVR(3);PR=ROTVR(4);ER=CRIGHT_ROT(4)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			
			RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVL(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
			RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)=ROTVR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
			

! 			IF (TURBULENCEMODEL.EQ.2)THEN
! 			PL=PL+((2.0D0/3.0D0)*EDDYFL(2))
! 			PR=PR+((2.0D0/3.0D0)*EDDYFR(2))  
! 
! 			END IF
			END IF


			FL(1)=RL*UL
			FL(2)=(RL*(UL**2))+PL
			FL(3)=RL*UL*VL
			
			FL(4)=UL*(EL+PL)
			
			
			FR(1)=RR*UR
			FR(2)=(RR*(UR**2))+PR
			FR(3)=RR*UR*VR
			
			FR(4)=UR*(ER+PR)
			
			IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
			FL(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=RML(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UL
			FR(5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=RMR(1:0+TURBULENCEEQUATIONS+PASSIVESCALAR)*UR
			END IF
			
			
			
			
! 			
! 			FLSTAR(:)=FL(:)+SL(1)*(ULSTAR(:)-CLEFT_ROT(:))
! 			FRSTAR(:)=FR(:)+SR(1)*(URSTAR(:)-CRIGHT_ROT(:))
			
			
			
			sl(1)=abs(ul)+sqrt(gamma*pl/rl)
			sr(1)=abs(ur)+sqrt(gamma*pr/rr)
			
			
			
			HLLCFLUX(:)=0.5d0*(FL(:)+FR(:))-0.5d0*MAX(ABS(SL(1)),ABS(SR(1)))*(Cright_ROT(:)-Cleft_ROT(:))
			


END SUBROUTINE RUSANOV_RIEMANN_SOLVER2d




SUBROUTINE ESTIMATE_WAVES2d(N,ROTVL,ROTVR,SL,SM,SR,GAMMA)
!> @brief
!> Waves speed estimates in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,K
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVL
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVR
	REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::SL,SR,SM
	REAL,INTENT(IN)::GAMMA
	REAL::CL,CR,PR,PL,UL,UR,VL,VR,WR,WL,RL,RR
	REAL::CUP,PPV,PMIN,PMAX,QMAX,QUSER,BL,BR,COV,PM,UM
	REAL::G1,G2,G3,G4,G5,G6,G7,G8,GEL,GER,PQ,PTL,PTR
	G1 = (GAMMA - 1.0d0)/(2.0d0*GAMMA)
    	G2 = (GAMMA + 1.0d0)/(2.0d0*GAMMA)
   	G3 = 2.0d0*GAMMA/(GAMMA - 1.0d0)
    	G4 = 2.0d0/(GAMMA - 1.0d0)
    	G5 = 2.0d0/(GAMMA + 1.0d0)
    	G6 = (GAMMA - 1.0d0)/(GAMMA + 1.0d0)
   	G7 = (GAMMA - 1.0d0)/2.0d0
   	G8 = GAMMA - 1.0d0
	
	SL=0.0d0
	SR=0.0d0
	SM=0.0d0
	COV=0.0d0
	!BUILD LEFT STATE VARIABLES
	RL=ROTVL(1)
	UL=ROTVL(2)
	VL=ROTVL(3)
	
	PL=ROTVL(4)
	CL=SQRT((PL*GAMMA)/(RL))
	!BUILD RIGHT  STATE VARIABLES
	RR=ROTVR(1)
	UR=ROTVR(2)
	VR=ROTVR(3)
	
	PR=ROTVR(4)
	CR=SQRT((PR*GAMMA)/(RR))

	CUP=0.25d0*(RL+RR)*(CL+CR)
	PPV=0.5d0*(PL + PR) + 0.5d0*(UL - UR)*CUP
	PPV=MAX(0.0d0,PPV)
	PMIN=MIN(PL,PR)
	PMAX=MAX(PL,PR)
	QMAX=PMAX/PMIN
	QUSER=2.0d0

	 IF(QMAX.LE.QUSER.AND.(PMIN.LE.PPV.AND.PPV.LE.PMAX))THEN
  
!        Select PRVS Riemann solver
 
         PM = PPV
         UM = 0.5d0*(UL + UR) + 0.5d0*(PL - PR)/CUP 
  
      	ELSE
 
         BL = 1.0d0 - COV*RL
         BR = 1.0d0 - COV*RR
 
         IF(PPV.LT.PMIN)THEN
 
!           Select Two-Rarefaction Riemann solver
 
        
            PQ  = (PL/PR)**G1
            UM  = (PQ*UL/CL/BL + UR/CR/BR + G4*(PQ - 1.0d0)) 
            UM  = UM/(PQ/CL/BL + 1.0d0/CR/BR)
            PTL = 1.0d0 + G7*(UL - UM)/CL/BL
            PTR = 1.0d0 + G7*(UM - UR)/CR/BR
            PM  = 0.5d0*(PL*PTL**G3 + PR*PTR**G3)
         ELSE

!           Use Two-Shock Riemann solver with PVRS as estimate
 
!          introduce iterations with PVRS as initial guess
 
            DO K=1,20
             GEL = SQRT((G5*BL/RL)/(G6*PL + PPV))
             GER = SQRT((G5*BR/RR)/(G6*PR + PPV))
             PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER)
             UM  = 0.5d0*(UL + UR) + 0.5d0*(GER*(PM - PR) - GEL*(PM - PL))	     
             IF ( ABS((PM-PPV)/PM) .LE. 1D-8) GOTO 101
                PPV = PM
     	    END DO
         ENDIF
      ENDIF
      
      101 continue

!     Find speeds
 
      IF(PM.LE.PL)THEN
         SL(1) = UL - CL
      ELSE
         SL(1) = UL - CL*SQRT(1.0d0 + G2*(PM/PL - 1.0d0))
      ENDIF
 
      SM(1)= UM
 
      IF(PM.LE.PR)THEN
         SR(1)= UR + CR
      ELSE
         SR(1) = UR + CR*SQRT(1.0d0 + G2*(PM/PR - 1.0d0))
      ENDIF

END SUBROUTINE ESTIMATE_WAVES2d








END MODULE RIEMANN
