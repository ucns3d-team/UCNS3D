MODULE ADVANCE
USE LIBRARY
USE TRANSFORM
USE FLUXES
USE SOURCE
USE INITIALISATION
USE BOUNDARY
USE RECON
USE LOCAL
USE FLOW_OPERATIONS
USE COMMUNICATIONS
USE implicit_time
USE implicit_FLUXES
USE DECLARATION
USE IO
USE MOODR
IMPLICIT NONE

 CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!SUBROUTINE CALLED TO CALCULATE CFL NUMBER DEPENDING ON THE SCHEME!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!EMPLOYED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE CALCULATE_CFL(N)
!> @brief
!> subroutine for computing the global time step size in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,L,KMAXE,J,INGTMAX,INGTMIN,WHGU,WHGL,SRF
REAL::SUVI,SUV3,maxU,MINU
REAL::CCFL,VELN,AGRT
real,dimension(1:nof_Variables)::leftv,rightv
real,dimension(1:nof_Variables)::SRF_SPEED
real::MP_PINFL,gammal
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
REAL,DIMENSION(1)::ETVM

KMAXE=XMPIELRANK(N)
       
  CCFL=(CFL/3.0d0)
        
  DT=tolbig
	IF (ITESTCASE.LT.3)THEN
    !$OMP DO REDUCTION (MIN:DT)
    DO I=1,KMAXE
      VELN=MAX(ABS(LAMx),ABS(LAMy),ABS(LAMz))
      
      if (dg.eq.1)then
        DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))*(1.0D0/(2*IORDER+1)))
      else
        DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN))))
      end if
      
    END DO
    !$OMP END DO
	END IF
	
	IF (ITESTCASE.EQ.3) THEN
	  !$OMP DO REDUCTION (MIN:DT)
    DO I=1,KMAXE
        
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  IF (multispecies.eq.1) THEN
		    AGRT=SQRT((LEFTV(5)+MP_PINFL)*GAMMAl/LEFTV(1))
		  ELSE
		    AGRT=SQRT(LEFTV(5)*GAMMA/LEFTV(1))
		  END IF
      IF (RFRAME.eq.0) THEN
          VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
      END IF
      IF (SRFG.EQ.1) THEN
          POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
          POY(1:3)=SRF_VELOCITY
          SRF_SPEED=ZERO
          SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
          VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
      END IF
      IF(MRF.EQ.1)THEN
        SRF=ILOCAL_RECON3(I)%MRF
        IF (ILOCAL_RECON3(I)%MRF.EQ.0) THEN
          VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
        ELSE
          POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
          POX=POX-ILOCAL_RECON3(I)%MRF_ORIGIN
          POY(1:3)=ILOCAL_RECON3(I)%MRF_VELOCITY
          SRF_SPEED=ZERO
          SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
          VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
        END IF
      END IF
      if (dg.eq.1)then
        DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))*(1.0D0/(2*IORDER+1)))
      else
        DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN))))
      END IF
		
    END DO
    !$OMP END DO
	END IF
	
	
	
	IF (ITESTCASE.EQ.4)THEN
	  !$OMP DO REDUCTION (MIN:DT)
    DO I=1,KMAXE
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  RIGHTV(1:NOF_vARIABLES)=LEFTV(1:NOF_vARIABLES)
		  CALL SUTHERLAND(N,leftv,rightv,VISCL,LAML)
		  AGRT=SQRT(LEFTV(5)*GAMMA/LEFTV(1))
                
      IF (RFRAME.EQ.0) THEN
        VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
      END IF
      IF(SRFG.EQ.1)THEN
        POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
        POY(1:3)=SRF_VELOCITY
        SRF_SPEED=ZERO
        SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
        VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
      END IF          
      IF(MRF.EQ.1)THEN
        SRF=ILOCAL_RECON3(I)%MRF
        IF (ILOCAL_RECON3(I)%MRF.EQ.0) THEN
          VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
        ELSE
          POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
          POX(1:3)=POX(1:3)-ILOCAL_RECON3(I)%MRF_ORIGIN(1:3)
          POY(1:3)=ILOCAL_RECON3(I)%MRF_VELOCITY(1:3)
          SRF_SPEED=ZERO
          SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
          VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
        END IF
      END IF
      IF (TURBULENCE.EQ.1) THEN
        IF (TURBULENCEMODEL.EQ.1) THEN
          TURBMV(1)=U_CT(I)%VAL(1,1);  TURBMV(2)=U_CT(I)%VAL(1,1);
          eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
          CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
          LAML(1)=LAML(1)+LAML(3)
          VISCL(1)=VISCL(1)+VISCL(3)
        END IF
      END IF

      if (dg.eq.1)then
        DT=MIN(DT,(CCFL/(2*IORDER+1))*(IELEM(N,I)%MINEDGE/((ABS(VELN))+(2.0D0*MAX(((4.0/3.0)*VISCL(1)/LEFTV(1)),GAMMA*LAML(1)/(PRANDTL*LEFTV(1)))*((2*IORDER+1)/IELEM(N,I)%MINEDGE)))))
      else
        IELEM(N,I)%VISCX=VISCL(1)/VISC
        DT=MIN(DT,CCFL*(1.0D0/((ABS(VELN)/((IELEM(N,I)%MINEDGE))) + (0.5D0*(LAML(1)+VISCL(1))/((IELEM(N,I)%MINEDGE))**2))))
        ! DT=MIN(DT,(CCFL)*(IELEM(N,I)%MINEDGE/((ABS(VELN))+(2.0D0*MAX(((4.0/3.0)*VISCL(1)/LEFTV(1)),GAMMA*LAML(1)/(PRANDTL*LEFTV(1)))*(1.0D0/IELEM(N,I)%MINEDGE)))))
      end if
                             
    END DO
    !$OMP END DO
	END IF
	
  RETURN
        
END SUBROUTINE CALCULATE_CFL

SUBROUTINE CALCULATE_CFLL(N)
!> @brief
!> subroutine for computing the time step size for each cell in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,L,KMAXE,J,INGTMAX,INGTMIN,WHGU,WHGL,SRF
REAL::SUVI,SUV3,maxU,MINU
REAL::CCFL,VELN,AGRT
real,dimension(1:nof_Variables)::leftv,rightv
real,dimension(1:nof_Variables)::SRF_SPEED
real::MP_PINFL,gammal
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
REAL,DIMENSION(1)::ETVM
KMAXE=XMPIELRANK(N)
       
  CCFL=(CFL/3.0d0)
        
	IF (ITESTCASE.LT.3)THEN
	  !$OMP DO
    DO I=1,KMAXE
		  VELN=MAX(ABS(LAMx),ABS(LAMy),ABS(LAMz))
		  IELEM(N,I)%DTL=CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))
	  END DO
	  !$OMP END DO
	END IF
	
	IF (ITESTCASE.EQ.3)THEN
	  !$OMP DO
    DO I=1,KMAXE
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)

      IF (multispecies.eq.1)THEN
		    AGRT=SQRT((LEFTV(5)+MP_PINFL)*GAMMAl/LEFTV(1))
		  ELSE
		    AGRT=SQRT(LEFTV(5)*GAMMA/LEFTV(1))
		  END IF
      IF (RFRAME.eq.0) THEN
        VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
      END IF
      IF (SRFG.EQ.1) THEN
        POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
        POY(1:3)=SRF_VELOCITY
        SRF_SPEED=ZERO
        SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
        VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
      END IF
      IF(MRF.EQ.1)THEN
        SRF=ILOCAL_RECON3(I)%MRF
        IF (ILOCAL_RECON3(I)%MRF.EQ.1)THEN
          VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
        ELSE
          POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
          POX=POX-ILOCAL_RECON3(I)%MRF_ORIGIN
          POY(1:3)=ILOCAL_RECON3(I)%MRF_VELOCITY
          SRF_SPEED=ZERO
          SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
          VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
        END IF
      END IF

		  if (dg.eq.1)then
		
		    IELEM(N,I)%DTL=CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))*(1.0D0/(2*IORDER+1))
		  else
      
		    IELEM(N,I)%DTL=CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))
		  end if
	  END DO
	  !$OMP END DO
	END IF
	

	IF (ITESTCASE.EQ.4)THEN
	  !$OMP DO
    DO I=1,KMAXE
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		  CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
		  RIGHTV(1:NOF_vARIABLES)=LEFTV(1:NOF_vARIABLES)
		  CALL SUTHERLAND(N,leftv,rightv,VISCL,LAML)
      AGRT=SQRT(LEFTV(5)*GAMMA/LEFTV(1))
      IF (SRFG.EQ.0.AND.MRF.EQ.0) THEN
        VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
      END IF
      IF(SRFG.EQ.1)THEN
        POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
        POY(1:3)=SRF_VELOCITY
        SRF_SPEED=ZERO
        SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
        VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
      END IF
      IF(MRF.EQ.1)THEN
        SRF=ILOCAL_RECON3(I)%MRF
        IF (ILOCAL_RECON3(I)%MRF.EQ.0)THEN
          VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)),ABS(LEFTV(4)))+AGRT
        ELSE
          POX(1)=IELEM(N,I)%XXC;POX(2)=IELEM(N,I)%YYC;POX(3)=IELEM(N,I)%ZZC
          POX(1:3)=POX(1:3)-ILOCAL_RECON3(I)%MRF_ORIGIN(1:3)
          POY(1:3)=ILOCAL_RECON3(I)%MRF_VELOCITY
          SRF_SPEED=ZERO
          SRF_SPEED(2:4)=VECT_FUNCTION(POX,POY)
          VELN=MAX(ABS(LEFTV(2)-SRF_SPEED(2)),ABS(LEFTV(3)-SRF_SPEED(3)),ABS(LEFTV(4)-SRF_SPEED(4)))+AGRT
        END IF
      END IF
      IF (TURBULENCE.EQ.1)THEN
		    IF (TURBULENCEMODEL.EQ.1)THEN
		      TURBMV(1)=U_CT(I)%VAL(1,1);  TURBMV(2)=U_CT(I)%VAL(1,1);
		      eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		      CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
		      LAML(1)=LAML(1)+LAML(3)
		      VISCL(1)=VISCL(1)+VISCL(3)
		    END IF
		  END IF
		
		  if (dg.eq.1)then
		    IELEM(N,I)%DTL=(CCFL/(2*IORDER+1))*(IELEM(N,I)%MINEDGE/((ABS(VELN))+(2.0D0*MAX(((4.0/3.0)*VISCL(1)/LEFTV(1)),GAMMA*LAML(1)/(PRANDTL*LEFTV(1)))*((2*IORDER+1)/IELEM(N,I)%MINEDGE))))
		  else
		    IELEM(N,I)%DTL=CCFL*(1.0D0/((ABS(VELN)/((IELEM(N,I)%MINEDGE))) + (0.5D0*(LAML(1)+VISCL(1))/((IELEM(N,I)%MINEDGE))**2)))
		  end if
		
	  END DO
	  !$OMP END DO
	END IF
	

  RETURN
        
END SUBROUTINE CALCULATE_CFLL



SUBROUTINE CALCULATE_CFL2D(N)
!> @brief
!> subroutine for computing the global time step size in 2D
  IMPLICIT NONE
  INTEGER,INTENT(IN)::N
  INTEGER::I,K,L,KMAXE,J,INGTMAX,INGTMIN,WHGU,WHGL
  REAL::SUVI,SUV3,maxU,MINU,sum_DT1,sum_dt2
  REAL::CCFL,VELN,AGRT,lamxl,lamyl
  real,dimension(1:nof_Variables)::leftv,rightv
  real,dimension(1:nof_Variables)::SRF_SPEED
  real::MP_PINFL,gammal
  REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
  REAL,DIMENSION(1:4)::VISCL,LAML
  REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
  REAL,DIMENSION(1:2)::TURBMV
  REAL,DIMENSION(1)::ETVM
  KMAXE=XMPIELRANK(N)
       
  CCFL=(CFL/2.0d0)
  
  DT=tolbig
  
  
  !$OMP BARRIER
        
	IF (ITESTCASE.LT.3)THEN
	  !$OMP DO REDUCTION (MIN:DT)
    DO I=1,KMAXE
        
      IF (initcond.eq.3)THEN
        lamxl=-ielem(n,i)%yyc+0.5d0
        lamyl=ielem(n,i)%xxc-0.5
      Else
        lamxl=lamx
        lamyl=lamy
      end if
        
		  VELN=MAX(ABS(LAMxl),ABS(LAMyl))
		
      if (dg.eq.1)then
        DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))*(1.0D0/(2*IORDER+1)))
      else

        DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN))))
      end if
		
    END DO
    !$OMP END DO
	END IF
	
	IF (ITESTCASE.EQ.3)THEN
	  !$OMP DO REDUCTION (MIN:DT)
    DO I=1,KMAXE
        
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		
		  CALL cons2prim(N,leftv,MP_PINFl,gammal)
      IF (multispecies.eq.1)THEN
        AGRT=SQRT((LEFTV(4)+MP_PINFL)*GAMMAl/LEFTV(1))
      ELSE
        AGRT=SQRT(LEFTV(4)*GAMMA/LEFTV(1))
      END IF
      VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)))+AGRT
		
		  if (dg.eq.1)then
		    DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))*(1.0D0/(2*IORDER+1)))
		  else
		    DT=MIN(DT,CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN))))
		  end if

    END DO
    !$OMP END DO
	END IF
	
	
	IF (ITESTCASE.EQ.4)THEN
	  !$OMP DO REDUCTION (MIN:DT)
    DO I=1,KMAXE
      LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
      CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
      RIGHTV(1:NOF_vARIABLES)=LEFTV(1:NOF_vARIABLES)
      CALL SUTHERLAND2D(N,leftv,rightv,VISCL,LAML)
      AGRT=SQRT(LEFTV(4)*GAMMA/LEFTV(1))
      VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)))+AGRT
      IF (TURBULENCE.EQ.1) THEN
        IF (TURBULENCEMODEL.EQ.1) THEN
          TURBMV(1)=U_CT(I)%VAL(1,1);  TURBMV(2)=U_CT(I)%VAL(1,1);
          eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
          CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
          LAML(1)=LAML(1)+LAML(3)
          VISCL(1)=VISCL(1)+VISCL(3)
        END IF
      END IF

      if (dg.eq.1)then
        DT=MIN(DT,(CCFL/(2*IORDER+1))*(IELEM(N,I)%MINEDGE/((ABS(VELN))+(2.0D0*MAX(((4.0/3.0)*VISCL(1)/LEFTV(1)),GAMMA*LAML(1)/(PRANDTL*LEFTV(1)))*((2*IORDER+1)/IELEM(N,I)%MINEDGE)))))
      else
		    DT=MIN(DT,CCFL*(1.0D0/((ABS(VELN)/((IELEM(N,I)%MINEDGE))) + (0.5D0*(LAML(1)+VISCL(1))/((IELEM(N,I)%MINEDGE))**2))))
      END IF
  
    END DO
    !$OMP END DO
	END IF
	
  RETURN
        
END SUBROUTINE CALCULATE_CFL2D


SUBROUTINE CALCULATE_CFLL2D(N)
!> @brief
!> subroutine for computing the time step size for each cell in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,K,L,KMAXE,J,INGTMAX,INGTMIN,WHGU,WHGL
REAL::SUVI,SUV3,maxU,MINU
REAL::CCFL,VELN,AGRT
real,dimension(1:nof_Variables)::leftv,rightv
real,dimension(1:nof_Variables)::SRF_SPEED
reaL::MP_PINFL,gammal,MP_PINFr,gammar
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:2)::TURBMV
REAL,DIMENSION(1)::ETVM
KMAXE=XMPIELRANK(N)
       
  CCFL=(CFL/2.0d0)
           
	IF (ITESTCASE.LT.3)THEN
	  !$OMP DO
    DO I=1,KMAXE
		  VELN=MAX(ABS(LAMx),ABS(LAMy))
		  IELEM(N,I)%DTL=CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))
	  END DO
	  !$OMP END DO
	END IF
	
	IF (ITESTCASE.EQ.3)THEN
	  !$OMP DO
    DO I=1,KMAXE
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		
		  CALL cons2prim(N,leftv,MP_PINFl,gammal)
		  AGRT=SQRT(LEFTV(4)*GAMMA/LEFTV(1))
		  VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)))+AGRT
		  if (dg.eq.1)then
		
		    IELEM(N,I)%DTL=CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))*(1.0D0/(2*IORDER+1))
		  else
		
		    IELEM(N,I)%DTL=CCFL*((IELEM(N,I)%MINEDGE)/(ABS(VELN)))
		  end if
		
	  END DO
	  !$OMP END DO
	END IF
	
	
	
	IF (ITESTCASE.EQ.4)THEN
	  !$OMP DO
    DO I=1,KMAXE
		  LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		  CALL cons2prim(N,leftv,MP_PINFl,gammal)
		  RIGHTV(1:NOF_vARIABLES)=LEFTV(1:NOF_vARIABLES)
		  CALL SUTHERLAND2D(N,leftv,rightv,VISCL,LAML)
		  AGRT=SQRT(LEFTV(4)*GAMMA/LEFTV(1))
		  VELN=MAX(ABS(LEFTV(2)),ABS(LEFTV(3)))+AGRT
		  IF (TURBULENCE.EQ.1)THEN
		    IF (TURBULENCEMODEL.EQ.1)THEN
		      TURBMV(1)=U_CT(I)%VAL(1,1);  TURBMV(2)=U_CT(I)%VAL(1,1);
		      eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		      CALL EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
		      LAML(1)=LAML(1)+LAML(3)
		      VISCL(1)=VISCL(1)+VISCL(3)
		    END IF
		  END IF

		  if (dg.eq.1)then
        IELEM(N,I)%DTL=(CCFL/(2*IORDER+1))*(IELEM(N,I)%MINEDGE/((ABS(VELN))+(2.0D0*MAX(((4.0/3.0)*VISCL(1)/LEFTV(1)),GAMMA*LAML(1)/(PRANDTL*LEFTV(1)))*((2*IORDER+1)/IELEM(N,I)%MINEDGE))))
      else
		    IELEM(N,I)%DTL=CCFL*(1.0D0/((ABS(VELN)/((IELEM(N,I)%MINEDGE))) + (0.5D0*(LAML(1)+VISCL(1))/((IELEM(N,I)%MINEDGE))**2)))
		  end if
	  END DO
	  !$OMP END DO
	END IF
	
  RETURN
        
END SUBROUTINE CALCULATE_CFLL2D


! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !---------------------------------------------------------------------------------------------!
! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED TO ADVANCE SOLUTION BY ONE TIME STEP SIZE!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE RUNGE_KUTTA3_MOOD(N)
!> @brief
!> SSP RUNGE KUTTA 3RD-ORDER SCHEME
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::IAVR,nvar,I,KMAXE,inds
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0	


IF (MOOD.EQ.1)THEN
  INDS=4
ELSE
  INDS=1
END IF

IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_DIFFUSIVE(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_DIFFUSIVE(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
END IF


!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  IF (MOOD.EQ.1) THEN
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  END IF
  U_C(I)%VAL(INDS,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF
 
IF (MOOD.EQ.1)THEN

  if (MOOD_MODE.eq.3) then
    Call EXCHANGE_HIGHER_MOOD(N, INDS)
  endif
  CALL MOOD_OPERATOR_2(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1) THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(INDS,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
      IELEM(N,I)%MOOD_O=2
    END IF
  END DO
  !$OMP END DO

  if (MOOD_MODE.eq.3) then
    Call EXCHANGE_HIGHER_MOOD(N, INDS)
  endif
  CALL MOOD_OPERATOR_1(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1) THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
      IELEM(N,I)%MOOD_O=1
    ELSE
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
    END IF
  END DO
  !$OMP END DO
END IF

 
IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_DIFFUSIVE(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
    CALL CALCULATE_FLUXESHI(N)
      CASE(3)
    CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_DIFFUSIVE(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
END IF

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  U_C(I)%VAL(inds,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0)) THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(TO4*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(OO4*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((OO4))*((DT)*&
  ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

IF (MOOD.EQ.1)THEN
 
  if (MOOD_MODE.eq.3) then
    Call EXCHANGE_HIGHER_MOOD(N, INDS)
  endif
  CALL MOOD_OPERATOR_2(N)
  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1) THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(INDS,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
      ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    IELEM(N,I)%MOOD_O=2
    END IF
  END DO
  !$OMP END DO

  if (MOOD_MODE.eq.3) then
    Call EXCHANGE_HIGHER_MOOD(N, INDS)
  endif
  CALL MOOD_OPERATOR_1(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(1,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
      ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=1
    ELSE
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
    END IF
  END DO
  !$OMP END DO

END IF
 

IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_DIFFUSIVE(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_DIFFUSIVE(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
        CALL VORTEXCALC(N)
    END SELECT
END IF
!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (MOOD.EQ.1)THEN
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  END IF
  U_C(I)%VAL(INDS,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
  ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
END DO
!$OMP END DO


IF (MOOD.EQ.1)THEN
 
  if (MOOD_MODE.eq.3) then
    Call EXCHANGE_HIGHER_MOOD(N, INDS)
  endif
  CALL MOOD_OPERATOR_2(N)
  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(INDS,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
      ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=2
    END IF
  END DO
  !$OMP END DO

  if (MOOD_MODE.eq.3) then
    Call EXCHANGE_HIGHER_MOOD(N, INDS)
  endif
  CALL MOOD_OPERATOR_1(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(1,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
      ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=1
    ELSE
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
    END IF
  END DO
  !$OMP END DO
END IF


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=((OO3)*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+((TO3)*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(((TO3))*&
    ((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN

 CALL AVERAGING_T(N)
 
END IF
                    
END SUBROUTINE RUNGE_KUTTA3_MOOD


SUBROUTINE RUNGE_KUTTA3(N)
!> @brief
!> SSP RUNGE KUTTA 3RD-ORDER SCHEME
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0	

CALL CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
    IF (DG == 1) THEN
        U_C(I)%VALDG(2,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(2,1:NOF_VARIABLES,:) - DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
    ELSE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
    END IF
END DO
!$OMP END DO

! IF ((DG.EQ.1).AND.(FILTERING.EQ.1))THEN
!     CALL SOL_INTEG_DGx(N)
!     CALL APPLY_FILTER_DG(N)
! END IF

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
        U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
    END DO
    !$OMP END DO
END IF
 
CALL CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
    IF (DG == 1) THEN
        U_C(I)%VALDG(3,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
         U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=TO4*U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + OO4*U_C(I)%VALDG(3,1:NOF_VARIABLES,:) - OO4*DT* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))

    ELSE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(1,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
          ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
        U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(TO4*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(OO4*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((OO4))*((DT)*&
          ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
    END DO
    !$OMP END DO
END IF

CALL CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
    IF (DG == 1) THEN
        ! U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=OO3*U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + TO3*U_C(I)%VALDG(1,1:NOF_VARIABLES,:) - TO3*DT*TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
        ! CALL DGEMM('N','N',NUM_DG_DOFS,nof_variables,NUM_DG_DOFS,ALPHA,m_1(i)%val(1:NUM_DG_DOFS,1:NUM_DG_DOFS),NUM_DG_DOFS,&
        ! RHS(I)%VALDG(1:NUM_DG_DOFS,1:NOF_VARIABLES),&
        ! NUM_DG_DOFS,BETA,RHS(I)%SOL_MM_DG,NUM_DG_DOFS)

        RHS(I)%SOL_MM_DG(1:NUM_DG_DOFS,1:nof_variables)=matmul(m_1(i)%val(1:NUM_DG_DOFS,1:NUM_DG_DOFS),RHS(I)%VALDG(1:NUM_DG_DOFS,1:NOF_VARIABLES))

        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=OO3*U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + TO3*U_C(I)%VALDG(1,1:NOF_VARIABLES,:) - TO3*DT * TRANSPOSE(RHS(I)%SOL_MM_DG)

    ELSE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(1,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
          ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=((OO3)*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+((TO3)*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(((TO3))*&
          ((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
    END DO
    !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN

    CALL AVERAGING_T(N)
 
END IF
                       
END SUBROUTINE RUNGE_KUTTA3




SUBROUTINE RUNGE_KUTTA1(N)
!> @brief
!> SSP FORWARD EULER SCHEME
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,kx
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
reaL::t1,t2,t3
KMAXE=XMPIELRANK(N)




CALL CALL_FLUX_SUBROUTINES_3D

DO I=1,KMAXE
    IF (DG == 1)then
        if((U_C(i)%VALDG(1,1,1).ne. U_C(i)%VALDG(1,1,1))) THEN
            IF (N == 0) PRINT*, 'STOPPING BECAUSE NaNs1'
            STOP ! Stop if NaNs
        END IF
    end if
end do

!$OMP DO
DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  
    IF (DG == 1) THEN

        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:) - DT* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))!*OOVOLUME   
    else

        U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)-(dt*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
    end if
  
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    !$OMP DO
    DO I=1,KMAXE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
    END DO
    !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN

    CALL AVERAGING_T(N)
 
END IF

END SUBROUTINE RUNGE_KUTTA1



SUBROUTINE RUNGE_KUTTA2(N)
!> @brief
!> SSP RUNGE KUTTA 2ND-ORDER SCHEME
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0	

CALL CALL_FLUX_SUBROUTINES_3D


!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  
END DO
!$OMP END DO
IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    !$OMP DO
    DO I=1,KMAXE
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
      U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
    END DO
    !$OMP END DO
END IF

CALL CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  U_C(I)%VAL(1,1:NOF_VARIABLES)=(oo2*U_C(I)%VAL(2,1:NOF_VARIABLES))+(oo2*U_C(I)%VAL(1,1:NOF_VARIABLES))-(dt*oo2*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(oo2*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(oo2*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(dt*oo2*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
 END IF

IF (AVERAGING.EQ.1)THEN

  CALL AVERAGING_T(N)
 
END IF

                        
END SUBROUTINE RUNGE_KUTTA2


SUBROUTINE RUNGE_KUTTA5(N)
!> @brief
!> SSP RUNGE KUTTA 2ND-ORDER SCHEME FOR LOCAL TIME STEPPING
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)


CALL CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  
  IF (DG == 1) THEN
      U_C(I)%VALDG(2,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        
      U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(2,1:NOF_VARIABLES,:) - ielem(n,i)%dtl * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
         
    ELSE
    
      U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)-(ielem(n,i)%dtl*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)-(ielem(n,i)%dtl*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
 END IF


CALL CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  
  IF (DG == 1) THEN
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=(oo2*U_C(I)%VALDG(2,1:NOF_VARIABLES,:)) +(oo2*U_C(I)%VALDG(1,1:NOF_VARIABLES,:))- (ielem(n,i)%dtl *oo2* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES))))    
  ELSE
    U_C(I)%VAL(1,1:NOF_VARIABLES)=(oo2*U_C(I)%VAL(2,1:NOF_VARIABLES))+(oo2*U_C(I)%VAL(1,1:NOF_VARIABLES))-(ielem(n,i)%dtl*oo2*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  end if
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(oo2*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(oo2*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(ielem(n,i)%dtl*oo2*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF


IF (AVERAGING.EQ.1)THEN

  CALL AVERAGING_T(N)
 
END IF

                        
END SUBROUTINE RUNGE_KUTTA5


SUBROUTINE RUNGE_KUTTA5_2D(N)
!> @brief
!> SSP RUNGE KUTTA 2ND-ORDER SCHEME FOR LOCAL TIME STEPPING IN 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)


CALL CALL_FLUX_SUBROUTINES_2D


!$OMP DO
DO I=1,KMAXE

  IF (DG == 1) THEN
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VALDG(2,:,:)=U_C(I)%VALDG(1,:,:)
    U_C(I)%VALDG(1,:,:)=U_C(I)%VALDG(2,:,:) - (ielem(n,i)%dtl* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,:))))!*OOVOLUME
  ELSE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(ielem(n,i)%dtl*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  END IF

END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)-(ielem(n,i)%dtl*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
 END IF


CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE

  IF (DG == 1) THEN
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VALDG(1,:,:)=(oo2*U_C(I)%VALdg(2,:,:))+(oo2*U_C(I)%VALDG(1,:,:))-(ielem(n,i)%dtl* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,:))))
  ELSE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VAL(1,1:NOF_VARIABLES)=(oo2*U_C(I)%VAL(2,1:NOF_VARIABLES))+(oo2*U_C(I)%VAL(1,1:NOF_VARIABLES))-(ielem(n,i)%dtl*oo2*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(oo2*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(oo2*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(ielem(n,i)%dtl*oo2*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN

  CALL AVERAGING_T(N)
 
END IF
                     
END SUBROUTINE RUNGE_KUTTA5_2D


SUBROUTINE RUNGE_KUTTA2_2D(N)
!> @brief
!> SSP RUNGE KUTTA 2ND-ORDER SCHEME IN 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)


CALL CALL_FLUX_SUBROUTINES_2D


!$OMP DO
DO I=1,KMAXE
  IF (DG == 1) THEN
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VALDG(2,:,:)=U_C(I)%VALDG(1,:,:)
    U_C(I)%VALDG(1,:,:)=U_C(I)%VALDG(2,:,:) - (DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,:))))!*OOVOLUME
  ELSE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(dt*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  END IF
END DO
!$OMP END DO


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(dt*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF


CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
  IF (DG == 1) THEN
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VALDG(1,:,:)=OO2*U_C(I)%VALDG(2,:,:) + OO2*U_C(I)%VALDG(1,:,:) - (OO2*DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,:))))!*OOVOLUME
  ELSE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_C(I)%VAL(1,1:NOF_VARIABLES)=(oo2*U_C(I)%VAL(2,1:NOF_VARIABLES))+(oo2*U_C(I)%VAL(1,1:NOF_VARIABLES))-(dt*oo2*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))-(dt*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN

  CALL AVERAGING_T(N)
 
END IF
                        
END SUBROUTINE RUNGE_KUTTA2_2D



SUBROUTINE SOL_INTEG_DG(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL,DIMENSION(1:NOF_VARIABLES)::SOLUTION_INTEG2

  
KMAXE=XMPIELRANK(N)
!$OMP DO
DO I=1,KMAXE
  call SOLUTION_INTEG(I,SOLUTION_INTEG2)
  U_C(I)%VAL(1,1:NOF_VARIABLES)=SOLUTION_INTEG2(1:NOF_VARIABLES)
END DO
!$OMP END DO

 


END SUBROUTINE SOL_INTEG_DG


SUBROUTINE SOL_INTEG_DGx(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,iconsidered
REAL,DIMENSION(1:NOF_VARIABLES)::SOLUTION_INTEG2
REAL,DIMENSION(1:NOF_VARIABLES)::SOLUTION_INTEG_WEAK
REAL,DIMENSION(1:NOF_VARIABLES)::SOLUTION_INTEG_STRONG


KMAXE=XMPIELRANK(N)
!$OMP DO
DO I=1,KMAXE
  ICONSIDERED=I

  call SOLUTION_INTEG(I,SOLUTION_INTEG2)
  U_C(I)%VAL(1,:)=SOLUTION_INTEG2(1:NOF_VARIABLES)

  if (FILTERING.EQ.1)THEN
    call SOLUTION_INTEG_S(I,SOLUTION_INTEG_STRONG)
    call SOLUTION_INTEG_W(I,SOLUTION_INTEG_WEAK)

    U_CS(I)%VAL(1,:)=SOLUTION_INTEG_STRONG(1:NOF_VARIABLES)
    U_CW(I)%VAL(1,:)=SOLUTION_INTEG_WEAK(1:NOF_VARIABLES)
  END IF
END DO
!$OMP END DO




END SUBROUTINE SOL_INTEG_DGx

SUBROUTINE SOL_INTEG_DG_init(N)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL,DIMENSION(1:NOF_VARIABLES)::SOLUTION_INTEG2
 
KMAXE=XMPIELRANK(N)
!$OMP DO
DO I=1,KMAXE

  call SOLUTION_INTEG(I,SOLUTION_INTEG2)
  U_C(I)%VAL(1,:)=SOLUTION_INTEG2(1:NOF_VARIABLES)
  U_e(I)%VAL(1,:)=U_C(I)%VAL(1,:)
END DO
!$OMP END DO

END SUBROUTINE SOL_INTEG_DG_init







SUBROUTINE RUNGE_KUTTA1_2D(N)
!> @brief
!> SSP FORWARD EULER SCHEME IN 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,kx
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
REAL,DIMENSION(NOF_VARIABLES)::TEMPSOL
KMAXE=XMPIELRANK(N)


CALL CALL_FLUX_SUBROUTINES_2D

DO I=1,KMAXE
  IF (DG == 1)then
    if((U_C(i)%VALDG(1,1,1).ne. U_C(i)%VALDG(1,1,1))) THEN
      IF (N == 0) PRINT*, 'STOPPING BECAUSE NaNs1'
        STOP ! Stop if NaNs
    END IF

  end if
end do


!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  
  IF (DG == 1) THEN    
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:) - DT* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))  
  else
    U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)-(dt*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  end if
  ! WRITE(200+N,*)I,TEMPSOL,U_C(I)%VALDG(1,1:NOF_VARIABLES,1)
END DO
!$OMP END DO


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)-(dt*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF


IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

                      
END SUBROUTINE RUNGE_KUTTA1_2D





SUBROUTINE RUNGE_KUTTA3_2D(N)
!> @brief
!> SSP RUNGE KUTTA 3RD-ORDER SCHEME IN 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0	

CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
    IF (DG == 1) THEN
        U_C(I)%VALDG(2,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(2,1:NOF_VARIABLES,:) - DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
         
    ELSE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
    END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF


CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
    IF (DG == 1) THEN
        U_C(I)%VALDG(3,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=TO4*U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + OO4*U_C(I)%VALDG(3,1:NOF_VARIABLES,:) - OO4*DT* TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
    ELSE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(1,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
          ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(TO4*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(OO4*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((OO4))*((DT)*&
      ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
    IF (DG == 1) THEN
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=OO3*U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + TO3*U_C(I)%VALDG(1,1:NOF_VARIABLES,:) - TO3*DT*TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
    ELSE
        OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
        U_C(I)%VAL(1,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
          ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
    END IF
END DO
!$OMP END DO


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=((OO3)*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+((TO3)*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(((TO3))*&
    ((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

END SUBROUTINE RUNGE_KUTTA3_2D



SUBROUTINE RUNGE_KUTTA3_2D_MOOD(N)
!> @brief
!> SSP RUNGE KUTTA 3RD-ORDER SCHEME
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,INDS
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0	


IF (MOOD.EQ.1)THEN
  INDS=4
ELSE
  INDS=1
END IF

IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
    
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
END IF


!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  IF (MOOD.EQ.1)THEN
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  END IF
  U_C(I)%VAL(INDS,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF
 
IF (MOOD.EQ.1)THEN
 
  CALL MOOD_OPERATOR_2(N)
  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(INDS,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
      IELEM(N,I)%MOOD_O=2
    END IF
  END DO
  !$OMP END DO

  CALL MOOD_OPERATOR_1(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1) THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
      IELEM(N,I)%MOOD_O=1
    ELSE
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
    END IF
  END DO
  !$OMP END DO
END IF

 

IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
          CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
          CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
END IF

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  U_C(I)%VAL(inds,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
    ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
END DO
!$OMP END DO


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(TO4*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(OO4*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((OO4))*((DT)*&
      ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

IF (MOOD.EQ.1)THEN
 
  CALL MOOD_OPERATOR_2(N)
  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(INDS,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
        ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=2
    END IF
  END DO
  !$OMP END DO

  CALL MOOD_OPERATOR_1(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(1,1:NOF_VARIABLES)=(TO4*U_C(I)%VAL(2,1:NOF_VARIABLES))+(OO4*U_C(I)%VAL(3,1:NOF_VARIABLES))-(((OO4))*((DT)*&
        ((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=1
    ELSE
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
    END IF
  END DO
  !$OMP END DO
END IF
 

IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
        CALL VORTEXCALC2D(N)
    END SELECT
END IF
!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (MOOD.EQ.1)THEN
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
  END IF
  U_C(I)%VAL(INDS,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
    ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
END DO
!$OMP END DO


IF (MOOD.EQ.1)THEN

  CALL MOOD_OPERATOR_2(N)
  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1) THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(INDS,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
        ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=2
    END IF
  END DO
  !$OMP END DO

  CALL MOOD_OPERATOR_1(N)

  !$OMP DO
  DO I=1,KMAXE
    IF (IELEM(N,I)%RECALC.EQ.1)THEN
      OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
      U_C(I)%VAL(1,1:NOF_VARIABLES)=((OO3)*U_C(I)%VAL(2,1:NOF_VARIABLES))+((TO3)*U_C(I)%VAL(1,1:NOF_VARIABLES))-(((TO3))*&
        ((DT)*((RHS(I)%VAL(1:NOF_VARIABLES))*(OOVOLUME))))
      IELEM(N,I)%MOOD_O=1
    ELSE
      U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(4,1:NOF_VARIABLES)
    END IF
  END DO
  !$OMP END DO
END IF


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=((OO3)*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+((TO3)*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))-(((TO3))*&
      ((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

END SUBROUTINE RUNGE_KUTTA3_2D_MOOD



SUBROUTINE RUNGE_KUTTA4(N)
!> @brief
!> SSP RUNGE KUTTA 4TH-ORDER SCHEME
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
REAL::DUMPRACEIN,DUMPRACEOUT,flops_count
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0	



CALL CALL_FLUX_SUBROUTINES_3D


!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (DG == 1) THEN
    U_C(I)%VALDG(2,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(2,1:NOF_VARIABLES,:) - DT * 0.391752226571890 * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))!*OOVOLUME
  ELSE
    U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*0.391752226571890*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
  END IF
  
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*0.391752226571890*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF
 
 
 
 
IF (statistics.eq.1)THEN
    !$OMP BARRIER 
    !$OMP MASTER
    pr_t8=MPI_Wtime()
    prace_t7=pr_t8-pr_t7
    
    DUMPRACEIN=PRACE_t1
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t1=DUMPRACEOUT
    DUMPRACEIN=PRACE_t2
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t2=DUMPRACEOUT
    DUMPRACEIN=PRACE_t3
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t3=DUMPRACEOUT
    DUMPRACEIN=PRACE_t4
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t4=DUMPRACEOUT
    DUMPRACEIN=PRACE_t5
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t5=DUMPRACEOUT

    DUMPRACEIN=PRACE_t6
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t6=DUMPRACEOUT

    DUMPRACEIN=PRACE_t7
    CALL MPI_ALLREDUCE(DUMPRACEIN,DUMPRACEOUT,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)
    PRACE_t7=DUMPRACEOUT

    PRACE_TX1=PRACE_T2+PRACE_T4
    PRACE_TX2=PRACE_T1+PRACE_T3+PRACE_T5+PRACE_T6+PRACE_T7
    PRACE_TX3=PRACE_TX1+PRACE_TX2

   
    IF (N.EQ.0)THEN
      OPEN(133,FILE=STATFILE,FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
      WRITE(133,'(I6,1X,E11.4,E11.4,E11.4,E11.4,E11.4,E11.4,E11.4,E11.4,E11.4,E11.4)')it,PRACE_TX3,PRACE_TX1,PRACE_TX2,PRACE_T1,PRACE_T2,prace_t3,prace_t4,prace_t5,prace_t6,prace_t7
      CLOSE(133)
    END IF
    
    
    !$OMP END MASTER
    !$OMP BARRIER
    
END IF

call CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (DG == 1) THEN
      U_C(I)%VALDG(3,1:NOF_VARIABLES,:) = U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
      U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = 0.444370493651235 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + 0.555629506348765 * U_C(I)%VALDG(3,1:NOF_VARIABLES,:) - 0.368410593050371 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
  ELSE
      U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
      U_C(I)%VAL(1,1:NOF_VARIABLES)=(0.444370493651235 * U_C(I)%VAL(2,1:NOF_VARIABLES)) + (0.555629506348765 * U_C(I)%VAL(3,1:NOF_VARIABLES)) - 0.368410593050371 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
  END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.444370493651235*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.555629506348765*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((0.368410593050371))*((DT)*&
    ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

call CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (DG == 1) THEN
    U_C(I)%VALDG(4,1:NOF_VARIABLES,:) = U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = 0.620101851488403 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + 0.379898148511597 * U_C(I)%VALDG(4,1:NOF_VARIABLES,:) - 0.251891774271694 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
  ELSE
    U_C(I)%VAL(4,1:NOF_VARIABLES) = U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(1,1:NOF_VARIABLES) = 0.620101851488403 * U_C(I)%VAL(2,1:NOF_VARIABLES) + 0.379898148511597 * U_C(I)%VAL(4,1:NOF_VARIABLES) - 0.251891774271694 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
  END IF
END DO
!$OMP END DO


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(4,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.620101851488403*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.379898148511597*U_Ct(I)%VAL(4,1:turbulenceequations+passivescalar))-(((0.251891774271694))*((DT)*&
    ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

call CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (DG == 1) THEN
    U_C(I)%VALDG(5,1:NOF_VARIABLES,:) = U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
    U_C(I)%VALDG(6,1:NOF_VARIABLES,:) = - DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = 0.178079954393132 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + 0.821920045606868 * U_C(I)%VALDG(5,1:NOF_VARIABLES,:) - 0.544974750228521 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
  ELSE
    U_C(I)%VAL(5,1:NOF_VARIABLES) = U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(6,1:NOF_VARIABLES) = - DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
    U_C(I)%VAL(1,1:NOF_VARIABLES) = 0.178079954393132 * U_C(I)%VAL(2,1:NOF_VARIABLES) + 0.821920045606868 * U_C(I)%VAL(5,1:NOF_VARIABLES) - 0.544974750228521 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
  END IF
END DO
!$OMP END DO

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(5,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(6,1:turbulenceequations+passivescalar)=-((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME)))
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.178079954393132*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.821920045606868*U_Ct(I)%VAL(5,1:turbulenceequations+passivescalar))-(((0.544974750228521))*((DT)*&
    ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF


call CALL_FLUX_SUBROUTINES_3D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  IF (DG == 1) THEN
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = (0.00683325884039 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:)) + (0.517231671970585 * U_C(I)%VALDG(4,1:NOF_VARIABLES,:)) + (0.12759831133288 * U_C(I)%VALDG(5,1:NOF_VARIABLES,:)) + (0.34833675773694 * U_C(I)%VALDG(1,1:NOF_VARIABLES,:)) + (0.08460416338212 * U_C(I)%VALDG(6,1:NOF_VARIABLES,:)) - 0.22600748319395 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
  ELSE
    U_C(I)%VAL(1,1:NOF_VARIABLES) = (0.00683325884039 * U_C(I)%VAL(2,1:NOF_VARIABLES)) + (0.517231671970585 * U_C(I)%VAL(4,1:NOF_VARIABLES)) +  (0.12759831133288 * U_C(I)%VAL(5,1:NOF_VARIABLES)) + (0.34833675773694 * U_C(I)%VAL(1,1:NOF_VARIABLES)) + (0.08460416338212 * U_C(I)%VAL(6,1:NOF_VARIABLES)) - (0.22600748319395 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME)
  END IF
END DO
!$OMP END DO


IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.00683325884039*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.517231671970585*U_Ct(I)%VAL(4,1:turbulenceequations+passivescalar))+&
          (0.12759831133288*U_Ct(I)%VAL(5,1:turbulenceequations+passivescalar))+(0.34833675773694*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))+&
          (0.08460416338212*U_Ct(I)%VAL(6,1:turbulenceequations+passivescalar))-(0.22600748319395*(DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME)))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF


END SUBROUTINE RUNGE_KUTTA4


SUBROUTINE CALL_FLUX_SUBROUTINES_3D
IMPLICIT NONE
REAL::DUMPRACEIN,DUMPRACEOUT
INTEGER::KMAXE,i
KMAXE=XMPIELRANK(N)

  IF (statistics.eq.1)THEN
    !$OMP BARRIER
    !$OMP MASTER
    pr_T1=MPI_Wtime()
    !$OMP END MASTER
    !$OMP BARRIER
  END IF

  if (dg.eq.1)then
      CALL SOL_INTEG_DG(N) ! Calculates cell average of DG solution for FV
  END IF

  IF (DG.EQ.1)THEN
    !$OMP DO
    DO I=1,KMAXE
        ielem(n,i)%filtered=0
    END DO
    !$OMP END DO
  END IF

	IF ((DG.EQ.1).AND.(FILTERING.EQ.1))THEN
    CALL SOL_INTEG_DGx(N)
    CALL APPLY_FILTER_DG(N)
  END IF

  IF (statistics.eq.1)THEN
    !$OMP BARRIER
    !$OMP MASTER
    pr_t2=MPI_Wtime()
    prace_t1=pr_t2-pr_t1
    !$OMP END MASTER
    !$OMP BARRIER
  END IF

  IF (FASTEST.EQ.1) THEN
      CALL EXCHANGE_LOWER(N)
  ELSE
      CALL EXCHANGE_HIGHER(N)
  END IF
        
  IF (statistics.eq.1)THEN
    !$OMP BARRIER
    !$OMP MASTER
    pr_t3=MPI_Wtime()
    prace_t2=pr_t3-pr_t2
    !$OMP END MASTER
    !$OMP BARRIER
  END IF

  IF (DG == 1) THEN

        CALL RECONSTRUCT_DG(N) ! Extrapolates solution to faces
        IF(MULTISPECIES.EQ.1)THEN
          IF (BOUND_LIM == 1) THEN
            CALL VFBP_LIMITER
          END IF
        END IF

        CALL TROUBLE_INDICATOR1 ! Checks for troubled cells
        
    END IF
    
    CALL ARBITRARY_ORDER(N)
    
    IF (DG == 1) THEN
        CALL TROUBLE_INDICATOR2 ! Changes DG to FV
    end if

    CALL EXHBOUNDHIGHER(N)

    if (dg.eq.1)then

        CALL EXHBOUNDHIGHER_DG(N)

        IF (ITESTCASE.EQ.4)THEN

          IF( BR2_YN.eq.2) then
            CALL RECONSTRUCT_BR2_DG
            CALL EXHBOUNDHIGHER_DG2(N)
          END IF

          IF( BR2_YN.eq.0) then
            CALL VISCOUS_DG_GGS(N)
          END IF
        END IF
    end if
    
    IF (statistics.eq.1)THEN
      !$OMP BARRIER
      !$OMP MASTER
      pr_t4=MPI_Wtime()
      prace_t3=pr_t4-pr_t3
      !$OMP END MASTER
      !$OMP BARRIER
    END IF

    IF (statistics.eq.1)THEN
      !$OMP BARRIER
      !$OMP MASTER
      pr_t5=MPI_Wtime()
      prace_t4=pr_t5-pr_t4
      !$OMP END MASTER
      !$OMP BARRIER
    END IF

    if (adda.eq.1)then
      IF (RUNGEKUTTA.EQ.11)THEN
        IF (ISCOUN.EQ.1)THEN
          call fix_dissipation(n)
          call EXCHANGE_ADDA_DISS(N)
          call fix_dissipation2(n)
        END IF
      ELSE
        call fix_dissipation(n)
        call EXCHANGE_ADDA_DISS(N)
        call fix_dissipation2(n)
      END IF
    end if

    IF (statistics.eq.1)THEN
      !$OMP BARRIER
      !$OMP MASTER
      pr_t6=MPI_Wtime()
      prace_t5=pr_t6-pr_t5
      !$OMP END MASTER
      !$OMP BARRIER
    END IF

    !Modifies RHS
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
    
        if ((SOURCE_ACTIVE.EQ.1))then
          call SOURCES_COMPUTATION_ROT(N)
        end if
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_dIFfusive(N)
        if ((SOURCE_ACTIVE.EQ.1))then
          call SOURCES_COMPUTATION_ROT(N)
        end if
        IF (turbulence.eq.1)THEN
          CALL SOURCES_COMPUTATION(N)
        END IF
        CALL VORTEXCALC(N)

    END SELECT
	
    IF (INITCOND.EQ.95)THEN
      CALL ENSTROPHY_CALC(N)
    END IF

    IF (statistics.eq.1)THEN
      !$OMP BARRIER
      !$OMP MASTER
      pr_t7=MPI_Wtime()
      prace_t6=pr_t7-pr_t6
      !$OMP END MASTER
      !$OMP BARRIER
    END IF

    !SOL INTEGRATION TIME=PRACE_T1
    !COMMUNICATION TIME OF HALO CELLS=prace_t2
    !RECONSTRUCTION TIME=prace_t3
    !COMMUNICATION TIME OF EXBOUNDHIGHER=prace_t4
    !ADDA=prace_t5
    !FLUXES=PRACE_T6

    !UPDATE OF SOLUTION=prace_t7
    !TOTAL COMMUNICATION TIME=PRACE_T2+PRACE_T4=PRACE_TX1
    !TOTAL COMPUTATIONALS TIME=PRACE_T1+PRACE_T3+PRACE_T5+PRACE_T6+PRACE_T7=PRACE_TX2
    !TOTAL TIME=TOTAL COMMUNICATION TIME+TOTAL COMPUTATIONALS TIME=PRACE_TX3

END SUBROUTINE CALL_FLUX_SUBROUTINES_3D




SUBROUTINE CALL_FLUX_SUBROUTINES_2D
IMPLICIT NONE
INTEGER::I,ICONSIDERED

    if (dg.eq.1)then
        CALL SOL_INTEG_DG(N)
    END IF

    IF (FASTEST.EQ.1) THEN
        CALL EXCHANGE_LOWER(N)
    ELSE
        CALL EXCHANGE_HIGHER(N)
    END IF
    
    IF (DG == 1) THEN 
        CALL RECONSTRUCT_DG(N)
        CALL TROUBLE_INDICATOR1
    END IF
     
    CALL ARBITRARY_ORDER(N)
    
    IF (DG == 1) THEN
        CALL TROUBLE_INDICATOR2
    end if
    
    IF (BOUND_LIM == 1) THEN
        CALL VFBP_LIMITER
    END IF
    
    CALL EXHBOUNDHIGHER(N)
    
    if (dg.eq.1)then
      call EXHBOUNDHIGHER_dg(N)

      IF (ITESTCASE.EQ.4)THEN
        IF( BR2_YN.eq.2) then
          CALL RECONSTRUCT_BR2_DG
          CALL EXHBOUNDHIGHER_DG2(N)
        END IF

        IF( BR2_YN.eq.0) then
          CALL VISCOUS_DG_GGS(N)
        END IF
      END IF
    end if
    
    !Modifies RHS
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2D(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT


    !FOR TEST ONLY EXPERIMENTAL
!     IF (CODE_PROFILE.EQ.30)THEN
!
!         IF (PROBEI(N,1).gt.0) THEN
!           DO I=1,XMPIELRANK(N)
!             IF (IELEM(N,I)%IHEXGL.EQ.PROBEI(N,1))THEN
!             ICONSIDERED=i
!             END IF
!           END DO
!         CALL VERTEX_NEIGHBOURS_VALUES(N)
!         END IF
!     END IF
    !END TEST

END SUBROUTINE CALL_FLUX_SUBROUTINES_2D





SUBROUTINE RUNGE_KUTTA4_2D(N)
!> @brief
!> SSP RUNGE KUTTA 4TH-ORDER SCHEME IN 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,RK_STAGE
REAL::AVRGS,OOVOLUME,TO4,OO4,TO3,OO3
KMAXE=XMPIELRANK(N)
TO4=3.0D0/4.0D0
OO4=1.0D0/4.0D0
TO3=2.0D0/3.0D0
OO3=1.0D0/3.0D0
RK_STAGE = 0

CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME

    IF (DG == 1) THEN
        U_C(I)%VALDG(2,1:NOF_VARIABLES,:)=U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:)=U_C(I)%VALDG(2,1:NOF_VARIABLES,:) - DT * 0.391752226571890 * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))!*OOVOLUME
    ELSE
        U_C(I)%VAL(2,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(1,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)-(DT*0.391752226571890*(RHS(I)%VAL(1:NOF_VARIABLES)*OOVOLUME))
    END IF
END DO
!$OMP END DO

RK_STAGE = RK_STAGE + 1

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar)-(DT*0.391752226571890*(RHSt(I)%VAL(1:turbulenceequations+passivescalar)*OOVOLUME))
  END DO
  !$OMP END DO
END IF

CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
  OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
  
  IF (DG == 1) THEN
    U_C(I)%VALDG(3,1:NOF_VARIABLES,:) = U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
    U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = 0.444370493651235 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + 0.555629506348765 * U_C(I)%VALDG(3,1:NOF_VARIABLES,:) - 0.368410593050371 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
  ELSE
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(1,1:NOF_VARIABLES)=(0.444370493651235 * U_C(I)%VAL(2,1:NOF_VARIABLES)) + (0.555629506348765 * U_C(I)%VAL(3,1:NOF_VARIABLES)) - 0.368410593050371 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
  END IF
END DO
!$OMP END DO

RK_STAGE = RK_STAGE + 1

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.444370493651235*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.555629506348765*U_Ct(I)%VAL(3,1:turbulenceequations+passivescalar))-(((0.368410593050371))*((DT)*&
    ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    
    IF (DG == 1) THEN
        U_C(I)%VALDG(4,1:NOF_VARIABLES,:) = U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = 0.620101851488403 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + 0.379898148511597 * U_C(I)%VALDG(4,1:NOF_VARIABLES,:) - 0.251891774271694 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
    ELSE
        U_C(I)%VAL(4,1:NOF_VARIABLES) = U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(1,1:NOF_VARIABLES) = 0.620101851488403 * U_C(I)%VAL(2,1:NOF_VARIABLES) + 0.379898148511597 * U_C(I)%VAL(4,1:NOF_VARIABLES) - 0.251891774271694 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
    END IF
END DO
!$OMP END DO

RK_STAGE = RK_STAGE + 1

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(4,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.620101851488403*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.379898148511597*U_Ct(I)%VAL(4,1:turbulenceequations+passivescalar))-(((0.251891774271694))*((DT)*&
    ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF
 
CALL CALL_FLUX_SUBROUTINES_2D

!$OMP DO
DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    
    IF (DG == 1) THEN
        U_C(I)%VALDG(5,1:NOF_VARIABLES,:) = U_C(I)%VALDG(1,1:NOF_VARIABLES,:)
        U_C(I)%VALDG(6,1:NOF_VARIABLES,:) = - DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = 0.178079954393132 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:) + 0.821920045606868 * U_C(I)%VALDG(5,1:NOF_VARIABLES,:) - 0.544974750228521 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))
    ELSE
        U_C(I)%VAL(5,1:NOF_VARIABLES) = U_C(I)%VAL(1,1:NOF_VARIABLES)
        U_C(I)%VAL(6,1:NOF_VARIABLES) = - DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
        U_C(I)%VAL(1,1:NOF_VARIABLES) = 0.178079954393132 * U_C(I)%VAL(2,1:NOF_VARIABLES) + 0.821920045606868 * U_C(I)%VAL(5,1:NOF_VARIABLES) - 0.544974750228521 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME
    END IF
END DO
!$OMP END DO

RK_STAGE = RK_STAGE + 1

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(5,1:turbulenceequations+passivescalar)=U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)
    U_Ct(I)%VAL(6,1:turbulenceequations+passivescalar)=-((DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME)))
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.178079954393132*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.821920045606868*U_Ct(I)%VAL(5,1:turbulenceequations+passivescalar))-(((0.544974750228521))*((DT)*&
    ((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME))))
  END DO
  !$OMP END DO
END IF

CALL CALL_FLUX_SUBROUTINES_2D
IF (FASTEST /= 1)THEN
    IF (ITESTCASE.EQ.4)THEN
        CALL VORTEXCALC2D(N)
    END IF
END IF

!$OMP DO
DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    
    IF (DG == 1) THEN
        U_C(I)%VALDG(1,1:NOF_VARIABLES,:) = (0.00683325884039 * U_C(I)%VALDG(2,1:NOF_VARIABLES,:)) + (0.517231671970585 * U_C(I)%VALDG(4,1:NOF_VARIABLES,:)) + (0.12759831133288 * U_C(I)%VALDG(5,1:NOF_VARIABLES,:)) + (0.34833675773694 * U_C(I)%VALDG(1,1:NOF_VARIABLES,:)) + (0.08460416338212 * U_C(I)%VALDG(6,1:NOF_VARIABLES,:)) - 0.22600748319395 * DT * TRANSPOSE(MATMUL(m_1(i)%val(:,:), RHS(I)%VALDG(:,1:NOF_VARIABLES)))! * OOVOLUME
    ELSE
        U_C(I)%VAL(1,1:NOF_VARIABLES) = (0.00683325884039 * U_C(I)%VAL(2,1:NOF_VARIABLES)) + (0.517231671970585 * U_C(I)%VAL(4,1:NOF_VARIABLES)) +  (0.12759831133288 * U_C(I)%VAL(5,1:NOF_VARIABLES)) + (0.34833675773694 * U_C(I)%VAL(1,1:NOF_VARIABLES)) + (0.08460416338212 * U_C(I)%VAL(6,1:NOF_VARIABLES)) - (0.22600748319395 * DT * RHS(I)%VAL(1:NOF_VARIABLES) * OOVOLUME)
    END IF
END DO
!$OMP END DO

RK_STAGE = RK_STAGE + 1

IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    OOVOLUME=1.0D0/IELEM(N,I)%TOTVOLUME
    U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar)=(0.00683325884039*U_Ct(I)%VAL(2,1:turbulenceequations+passivescalar))+(0.517231671970585*U_Ct(I)%VAL(4,1:turbulenceequations+passivescalar))+&
          (0.12759831133288*U_Ct(I)%VAL(5,1:turbulenceequations+passivescalar))+(0.34833675773694*U_Ct(I)%VAL(1,1:turbulenceequations+passivescalar))+&
          (0.08460416338212*U_Ct(I)%VAL(6,1:turbulenceequations+passivescalar))-(0.22600748319395*(DT)*((RHSt(I)%VAL(1:turbulenceequations+passivescalar))*(OOVOLUME)))
  END DO
  !$OMP END DO
END IF

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

IF (DG == 1) then
  if (ALL(U_C(1)%VALDG(1,:,:) /= U_C(1)%VALDG(1,:,:))) THEN
    IF (N == 0) PRINT*, 'STOPPING BECAUSE NaNs'
      STOP ! Stop if NaNs
  end if
END IF

END SUBROUTINE RUNGE_KUTTA4_2D





SUBROUTINE IMPLICIT_TIMEs(N)
!> @brief
!> IMPLICIT APPROXIMATELY FACTORED TIME STEPPING SCHEME
IMPLICIT NONE
INTEGER::I,K,KMAXE,kill_nan
INTEGER,INTENT(IN)::N
reaL::verysmall
verysmall = tolsmall

KMAXE=XMPIELRANK(N)
IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        if ((SOURCE_ACTIVE.EQ.1))then
            call SOURCES_COMPUTATION_ROT(N)
        end if
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_diffusive(N)
        if ((SOURCE_ACTIVE.EQ.1))then
            call SOURCES_COMPUTATION_ROT(N)
        end if
        call VORTEXCALC(N)
        if (turbulence.eq.1)then
            call SOURCES_COMPUTATION(N)
        end if
    END SELECT
    
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        if ((SOURCE_ACTIVE.EQ.1))then
            call SOURCES_COMPUTATION_ROT(N)
        end if
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_diffusive(N)
        if ((SOURCE_ACTIVE.EQ.1))then
            call SOURCES_COMPUTATION_ROT(N)
        end if
        call VORTEXCALC(N)
        if (turbulence.eq.1)then
            call SOURCES_COMPUTATION(N)
        end if
    END SELECT
END IF

IF (RELAX.EQ.3)THEN
  CALL RELAXATION_LUMFREE(N)
ELSE
  IF (lowmemory.eq.0)THEN
    CALL RELAXATION(N)
  ELSE
    CALL RELAXATION_lm(N)
  END IF
END IF

kill_nan=0
 
!$OMP DO
DO I=1,KMAXE
    IF ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)).or.(impdu(i,5).ne.impdu(i,5)))THEN
        write(600+n,*)"nan present",ielem(n,i)%ihexgl,ielem(n,i)%ishape,ielem(n,i)%xxc, ielem(n,i)%yyc,ielem(n,i)%zzc
        write(600+n,*)ielem(n,i)%dih(:)
        write(500+n,'(3es14.6)')ielem(n,i)%xxc, ielem(n,i)%yyc,ielem(n,i)%zzc
        if(MRF.EQ.1)then
            write(700+n,*)'SRF -diagonal', ILOCAL_RECON3(I)%MRF ,I
            write(700+n,'(3es14.6)'),ielem(n,i)%xxc, ielem(n,i)%yyc,ielem(n,i)%zzc
            write(700+n,*),impdu(i,1),impdu(i,2),impdu(i,3),impdu(i,4),impdu(i,5)
        end if
        kill_nan=1
    END IF
    U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
END DO
!$OMP END DO
IF (kill_nan.eq.1)THEN
    stop
END IF

IF ((PASSIVESCALAR.GT.0).OR.(TURBULENCE.GT.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    do k=1,turbulenceequations+passivescalar
      IF (U_CT(I)%VAL(1,k)+IMPDU(I,5+k).ge.zero)THEN
        U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+0.4*IMPDU(i,5+k)
      END IF
    END do
  END DO
  !$OMP END DO
END IF

  
END SUBROUTINE IMPLICIT_TIMEs




SUBROUTINE IMPLICIT_TIMEs_2d(N) 
!> @brief
!> IMPLICIT APPROXIMATELY FACTORED TIME STEPPING SCHEME 2D
IMPLICIT NONE
INTEGER::I,K,KMAXE,kill_nan
INTEGER,INTENT(IN)::N
reaL::verysmall
verysmall = tolsmall

KMAXE=XMPIELRANK(N)
IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2d(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        ! CALL VORTEXCALC2D(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
    
ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2d(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        !CALL VORTEXCALC2D(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
END IF

IF (RELAX.EQ.3)THEN
 
    CALL RELAXATION_LUMFREE(N)
ELSE
    IF (lowmemory.eq.0)THEN
        CALL RELAXATION2d(N)
    ELSE
        CALL RELAXATION_lm2d(N)
    END IF
END IF

kill_nan=0
!$OMP DO
DO I=1,KMAXE

    IF ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)))THEN
        write(600+n,*)"nan present",ielem(n,i)%ihexgl,ielem(n,i)%ishape,ielem(n,i)%xxc, ielem(n,i)%yyc
        write(600+n,*)ielem(n,i)%dih(:)
        kill_nan=1
    END IF

    U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
END DO
!$OMP END DO

IF (kill_nan.eq.1)THEN
    stop
END IF

IF ((PASSIVESCALAR.GT.0).OR.(TURBULENCE.GT.0))THEN
  !$OMP DO
  DO I=1,KMAXE
    do k=1,turbulenceequations+passivescalar
      IF (ispal.eq.1)THEN
        IF (U_CT(I)%VAL(1,k)+IMPDU(I,4+k).ge.zero)THEN
          U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+0.4*IMPDU(i,4+k)
        END IF
      ELSE
        U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+0.4*IMPDU(i,4+k)
      END IF
    END do
  END DO
  !$OMP END DO

  ! IF (kill_nan.eq.2)THEN
  !   stop
  ! END IF

END IF

END SUBROUTINE IMPLICIT_TIMEs_2d





SUBROUTINE DUAL_TIME(N)
!> @brief
!> DUAL TIME STEPPING
IMPLICIT NONE
INTEGER::I,K,KMAXE,jj,kill_nan
INTEGER,INTENT(IN)::N
reaL::verysmall
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
verysmall = tolsmall

inner_tol=reslimit

KMAXE=XMPIELRANK(N)

IF (IT.EQ.RESTART)THEN
  !$OMP DO
  DO I=1,KMAXE 
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
    IF ((turbulence.gt.0).or.(passivescalar.gt.0)) THEN
      U_CT(I)%VAL(3,:)=U_CT(I)%VAL(1,:)
      U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    END IF
  END DO
  !$OMP END DO
END IF

firsti=0.0d0
DO JJ=1,upperlimit
  rsumfacei=zero;allresdt=zero;dummy3i=zero; 
  if (jj.eq.1)then
    iscoun=1
  else
    iscoun=2
  end if    
      
  CALL CALL_FLUX_SUBROUTINES_3D

  IF (relax.eq.3)THEN
    CALL RELAXATION_LUMFREE(N)
  ELSE
    IF (lowmemory.eq.0)THEN
      CALL RELAXATION(N)
    ELSE
      CALL RELAXATION_lm(N)
    END IF
  END IF

  kill_nan=0
 
  !$OMP BARRIER 
  !$OMP DO  REDUCTION(+:allresdt)
  DO I=1,KMAXE
      rsumfacei=sqrt(((IMPDU(I,1))**2)+((IMPDU(I,2))**2)+((IMPDU(I,3))**2)+((IMPDU(I,4))**2)+((IMPDU(I,5))**2))
      allresdt=allresdt+(rsumfacei*ielem(n,i)%totvolume)
      
      IF ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)).or.(impdu(i,5).ne.impdu(i,5)))THEN
          kill_nan=1
      END IF
  END do
  !$OMP END DO

  !$OMP BARRIER 
  IF (kill_nan.eq.1)THEN
      stop
  END IF

  !$OMP MASTER
  DUMMY3I=zero

  CALL MPI_ALLREDUCE(allresdt,DUMMY3i,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
  allresdt=dummy3i/TOTALVOLUME

  IF (allresdt.gt.firsti)THEN
      firsti=allresdt
  END IF

  allresdt=allresdt/firsti

  IF (n.eq.0)THEN
      OPEN(77,FILE='res1.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
      WRITE(77,*)allresdt,jj,it
      CLOSE(77)
  END IF

  !$OMP END MASTER
  !$OMP BARRIER 

  IF ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))THEN
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          IF (U_CT(I)%VAL(1,k)+IMPDU(I,nof_Variables+k).ge.zero)THEN
            U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
          END IF
        END do
      END IF
    END DO
    !$OMP END DO
    exit
  ELSE
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          IF (U_CT(I)%VAL(1,k)+IMPDU(I,nof_Variables+k).ge.zero)THEN
            U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
          END IF
        END do
      END IF
    END DO
    !$OMP END DO
  END IF

END DO

!$OMP DO
DO I=1,KMAXE 
  U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)
  U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
  !U_C(I)%VAL(1,1:nof_variables)=2.0*U_C(I)%VAL(2,1:nof_variables)-U_C(I)%VAL(3,1:nof_variables)
  
  IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    U_CT(I)%VAL(3,:)=U_CT(I)%VAL(2,:)
    U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    !U_CT(I)%VAL(1,:)=2.0*U_CT(I)%VAL(2,:)-U_CT(I)%VAL(3,:)
  end if
END DO
!$OMP END DO

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

END SUBROUTINE dual_time





SUBROUTINE DUAL_TIME_EX(N)
!> @brief
!> DUAL TIME STEPPING 2D
IMPLICIT NONE
INTEGER::I,K,KMAXE,jj
INTEGER,INTENT(IN)::N
reaL::verysmall
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
verysmall = tolsmall

inner_tol=reslimit

KMAXE=XMPIELRANK(N)

IF (IT.EQ.RESTART)THEN
  !$OMP DO
  DO I=1,KMAXE 
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
    IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
      U_CT(I)%VAL(3,:)=U_CT(I)%VAL(1,:)
      U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    END IF
  END DO
  !$OMP END DO
END IF

firsti=0.0d0
DO JJ=1,upperlimit
  rsumfacei=zero;allresdt=zero;dummy3i=zero; 
      
  IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_dIFfusive(N)
        CALL VORTEXCALC(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
    
  ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE(N)
        CALL CALCULATE_FLUXESHI_dIFfusive(N)
        CALL VORTEXCALC(N)
        IF (turbulence.eq.1)THEN
          CALL SOURCES_COMPUTATION(N)
        END IF
    END SELECT
  END IF

  CALL RELAXATION_EX(N)

  !$OMP BARRIER 
  !$OMP DO  REDUCTION(+:allresdt)
  DO I=1,KMAXE
      rsumfacei=sqrt(((IMPDU(I,1))**2)+((IMPDU(I,2))**2)+((IMPDU(I,3))**2)+((IMPDU(I,4))**2)+((IMPDU(I,5))**2))
        allresdt=allresdt+(rsumfacei*ielem(n,i)%totvolume)
  END do
  !$OMP END DO

  !$OMP BARRIER 
  !$OMP MASTER
  DUMMY3I=zero

  CALL MPI_ALLREDUCE(allresdt,DUMMY3i,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
  allresdt=dummy3i/TOTALVOLUME

  IF (allresdt.gt.firsti)THEN
    firsti=allresdt
  END IF

  allresdt=allresdt/firsti

  IF (n.eq.0)THEN
    write(777,*)allresdt,jj,it
  END IF

  !$OMP END MASTER
  !$OMP BARRIER 

  IF ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))THEN
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
        END do
      END IF
    END DO
    !$OMP END DO
    exit

  ELSE
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
        END do
      END IF
    END DO
    !$OMP END DO
  END IF

END DO

!$OMP DO
DO I=1,KMAXE 
  U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)
  U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
  U_C(I)%VAL(1,1:nof_variables)=2.0*U_C(I)%VAL(2,1:nof_variables)-U_C(I)%VAL(3,1:nof_variables)
  
  IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    U_CT(I)%VAL(3,:)=U_CT(I)%VAL(2,:)
    U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    U_CT(I)%VAL(1,:)=2.0*U_CT(I)%VAL(2,:)-U_CT(I)%VAL(3,:)
  END IF
END DO
!$OMP END DO

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

END SUBROUTINE DUAL_TIME_EX





SUBROUTINE DUAL_TIME_EX_2D(N)
!> @brief
!> DUAL TIME STEPPING 2D
IMPLICIT NONE
INTEGER::I,K,KMAXE,jj
INTEGER,INTENT(IN)::N
reaL::verysmall
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
verysmall = tolsmall

inner_tol=reslimit

KMAXE=XMPIELRANK(N)

IF (IT.EQ.RESTART)THEN
  !$OMP DO
  DO I=1,KMAXE 
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
    IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
      U_CT(I)%VAL(3,:)=U_CT(I)%VAL(1,:)
      U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    END IF
  END DO
  !$OMP END DO
END IF

firsti=0.0d0
DO JJ=1,upperlimit
  rsumfacei=zero;allresdt=zero;dummy3i=zero; 
      
  IF (FASTEST.EQ.1)THEN
    CALL EXCHANGE_LOWER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2d(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        CALL VORTEXCALC2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
    
  ELSE
    CALL EXCHANGE_HIGHER(N)
    CALL ARBITRARY_ORDER(N)
    CALL EXHBOUNDHIGHER(N)
    SELECT CASE(ITESTCASE)
      CASE(1,2)
        CALL CALCULATE_FLUXESHI2d(N)
      CASE(3)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
      CASE(4)
        CALL CALCULATE_FLUXESHI_CONVECTIVE2d(N)
        CALL CALCULATE_FLUXESHI_dIFfusive2d(N)
        CALL VORTEXCALC2d(N)
        IF (turbulence.eq.1)THEN
            CALL SOURCES_COMPUTATION2d(N)
        END IF
    END SELECT
  END IF

  CALL RELAXATION_EX(N)

  !$OMP BARRIER 
  !$OMP DO  REDUCTION(+:allresdt)
  DO I=1,KMAXE
        rsumfacei=sqrt(((IMPDU(I,1))**2)+((IMPDU(I,2))**2)+((IMPDU(I,3))**2)+((IMPDU(I,4))**2))
        allresdt=allresdt+(rsumfacei*ielem(n,i)%totvolume)
  END do
  !$OMP END DO

  !$OMP BARRIER 
  !$OMP MASTER
  DUMMY3I=zero

  CALL MPI_ALLREDUCE(allresdt,DUMMY3i,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
  allresdt=dummy3i/TOTALVOLUME

  IF (allresdt.gt.firsti)THEN
    firsti=allresdt
  END IF

  allresdt=allresdt/firsti

  IF (n.eq.0)THEN
    write(777,*)allresdt,jj,it
  END IF

  !$OMP END MASTER
  !$OMP BARRIER 

  IF ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))THEN
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
        END do
      END IF
    END DO
    !$OMP END DO
    exit

  ELSE
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
        END do
      END IF
    END DO
    !$OMP END DO
  END IF

END DO

!$OMP DO
DO I=1,KMAXE 
  U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)
  U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
  U_C(I)%VAL(1,1:nof_variables)=2.0*U_C(I)%VAL(2,1:nof_variables)-U_C(I)%VAL(3,1:nof_variables)
  
  IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    U_CT(I)%VAL(3,:)=U_CT(I)%VAL(2,:)
    U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    U_CT(I)%VAL(1,:)=2.0*U_CT(I)%VAL(2,:)-U_CT(I)%VAL(3,:)
  END IF
END DO
!$OMP END DO

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

END SUBROUTINE dual_time_ex_2d





SUBROUTINE DUAL_TIME_2d(N)
!> @brief
!> DUAL TIME STEPPING 2D
IMPLICIT NONE
INTEGER::I,K,KMAXE,nvar,jj,kill_nan
INTEGER,INTENT(IN)::N
reaL::verysmall
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
verysmall = tolsmall

inner_tol=reslimit

KMAXE=XMPIELRANK(N)

IF (IT.EQ.RESTART)THEN
  !$OMP DO
  DO I=1,KMAXE 
    U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
    U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
    
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
      U_CT(I)%VAL(3,:)=U_CT(I)%VAL(1,:)
      U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    end if
  END DO
  !$OMP END DO
END IF

firsti=0.0d0
DO JJ=1,upperlimit
  rsumfacei=zero;allresdt=zero;dummy3i=zero; 
  if (jj.eq.1)then
    iscoun=1
  else
    iscoun=2
  end if
      
  CALL CALL_FLUX_SUBROUTINES_2D

  IF (relax.eq.3)THEN
    CALL RELAXATION_LUMFREE(N)
  ELSE
    IF (lowmemory.eq.0)THEN
      CALL RELAXATION2d(N)
    ELSE
      CALL RELAXATION_lm2d(N)
    END IF
  END IF

  !$OMP BARRIER 
  !$OMP DO  REDUCTION(+:allresdt)
  DO I=1,KMAXE
      rsumfacei=sqrt(((IMPDU(I,1))**2)+((IMPDU(I,2))**2)+((IMPDU(I,3))**2)+((IMPDU(I,4))**2))
      allresdt=allresdt+(rsumfacei*ielem(n,i)%totvolume)

      IF ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)))THEN
        kill_nan=1
      END IF
  END do
  !$OMP END DO

  !$OMP BARRIER 

  IF (kill_nan.eq.1) THEN
    stop
  END IF

  !$OMP MASTER
  DUMMY3I=zero

  CALL MPI_ALLREDUCE(allresdt,DUMMY3i,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
  allresdt=dummy3i/TOTALVOLUME

  IF (allresdt.gt.firsti)THEN
    firsti=allresdt
  END IF

  allresdt=allresdt/firsti

  IF (n.eq.0)THEN
				  OPEN(77,FILE='res1.txt',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
				  WRITE(77,*)allresdt,jj,it
				  CLOSE(77)
  END IF

  !$OMP END MASTER
  !$OMP BARRIER 

  IF ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))THEN
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          IF (U_CT(I)%VAL(1,k)+IMPDU(I,nof_Variables+k).ge.zero)THEN
            U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
          END IF
        END do
      END IF
    END DO
    !$OMP END DO
    exit

  ELSE
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%VAL(1,1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)+IMPDU(I,1:nof_Variables)
      IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
        do k=1,turbulenceequations+passivescalar
          IF (U_CT(I)%VAL(1,k)+IMPDU(I,nof_Variables+k).ge.zero)THEN
            U_CT(I)%VAL(1,k)=U_CT(I)%VAL(1,k)+IMPDU(i,nof_Variables+k)
          END IF
        END do
      END IF
    END DO
    !$OMP END DO
  END IF
END DO

!$OMP DO
DO I=1,KMAXE 
  U_C(I)%VAL(3,1:NOF_VARIABLES)=U_C(I)%VAL(2,1:NOF_VARIABLES)
  U_C(I)%VAL(2,1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
  !U_C(I)%VAL(1,1:nof_variables)=(2.0*U_C(I)%VAL(2,1:nof_variables))-U_C(I)%VAL(3,1:nof_variables)
  
  IF ((turbulence.gt.0).or.(passivescalar.gt.0))THEN
    U_CT(I)%VAL(3,:)=U_CT(I)%VAL(2,:)
    U_CT(I)%VAL(2,:)=U_CT(I)%VAL(1,:)
    !U_CT(I)%VAL(1,:)=2.0*U_CT(I)%VAL(2,:)-U_CT(I)%VAL(3,:)
  end if
END DO
!$OMP END DO

IF (AVERAGING.EQ.1)THEN
  CALL AVERAGING_T(N)
END IF

END SUBROUTINE dual_time_2d





SUBROUTINE RELAXATION_EX(N)
IMPLICIT NONE
!> @brief
!> This subroutine solves the linear system for implicit time stepping either through MATRIX FREE LU-SGS low memory footprint
INTEGER,INTENT(IN)::N
INTEGER::I,L,K,II,SWEEPS,kmaxe,nvar,igoflux,icaseb,INDT1,INDT2,INDT3,IJK
real::dt1,dtau

kmaxe=xmpielrank(n)

impdu(:,:)=zero

INDT1=NOF_VARIABLES+1
INDT2=NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR
INDT3=TURBULENCEEQUATIONS+PASSIVESCALAR

!$OMP DO
DO I=1,KMAXE
  dt1=ielem(n,i)%totvolume/dt
  IMPDU(I,1:nof_variables)=DT1*(1.5D0*U_C(I)%VAL(1,1:nof_Variables)-2.0D0*U_C(I)%VAL(2,1:nof_Variables)+0.5D0*U_C(I)%VAL(3,1:nof_Variables))+RHS(I)%VAL(1:nof_variables)

  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
    impdu(I,INDT1:INDT2)=DT1*(1.5D0*U_CT(I)%VAL(1,1:INDT3)-2.0D0*U_CT(I)%VAL(2,1:INDT3)+0.5D0*U_CT(I)%VAL(3,1:INDT3))+RHST(I)%VAL(1:INDT3)
  END IF
END DO
!$OMP END DO 

!$OMP DO
DO I=1,KMAXE
  dtau=(ielem(n,i)%dtl/IELEM(N,I)%TOTVOLUME)*(1.0d0/(1.0D0+1.5D0*(ielem(n,i)%dtl/DT)))
  IMPDU(I,1:nof_variables)=-IMPDU(I,1:nof_variables)*dtau
  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
    impdu(I,INDT1:INDT2)=-impdu(I,INDT1:INDT2)*dtau
  END IF
END DO
!$OMP END DO 
                                       
END SUBROUTINE RELAXATION_EX





SUBROUTINE AVERAGING_T(N)
!> @brief
!> TEMPORAL AVERAGING FOR UNSTEADY SIMULATIONS
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,nvar
integer::ind1
KMAXE=XMPIELRANK(N)

IF (DIMENSIONA.EQ.3)THEN
  IF (rungekutta.eq.4)THEN
    ind1=7
  ELSE
    ind1=5
  END IF

  IF (tz1.gt.zero) THEN
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%val(ind1,:)=(((Tz1-dt)/(Tz1))*U_C(I)%val(ind1,:))+((dt*U_C(I)%VAL(1,:))/Tz1)
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) THEN
        DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
          U_CT(I)%val(ind1,NVAR)=(((Tz1-dt)/(Tz1))*U_CT(I)%val(ind1,NVAR))+((dt*U_CT(I)%VAL(1,NVAR))/(Tz1*U_C(I)%val(ind1,1)))
        END DO
      END IF
      !U,V,W,UV,UW,WV,PS
      U_C(I)%RMS(1)=SQRT(abs(((U_C(I)%RMS(1)**2)*((Tz1-dt)/(Tz1)))+(((U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1)-U_C(I)%val(ind1,2)/U_C(I)%val(ind1,1))**2)*dt/Tz1)))
      U_C(I)%RMS(2)=SQRT(abs(((U_C(I)%RMS(2)**2)*((Tz1-dt)/(Tz1)))+(((U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1)-U_C(I)%val(ind1,3)/U_C(I)%val(ind1,1))**2)*dt/Tz1)))
      U_C(I)%RMS(3)=SQRT(abs(((U_C(I)%RMS(3)**2)*((Tz1-dt)/(Tz1)))+(((U_C(I)%VAL(1,4)/U_C(I)%VAL(1,1)-U_C(I)%val(ind1,4)/U_C(I)%val(ind1,1))**2)*dt/Tz1)))
      
      U_C(I)%RMS(4)=(((U_C(I)%RMS(4))*((Tz1-DT)/(Tz1)))+&
        (((((U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1))-(U_C(I)%val(ind1,2)/U_C(I)%val(ind1,1)))*((U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1))-(U_C(I)%val(ind1,3)/U_C(I)%val(ind1,1)))))*DT/Tz1))
      U_C(I)%RMS(5)=(((U_C(I)%RMS(5))*((Tz1-DT)/(Tz1)))+&
        (((((U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1))-(U_C(I)%val(ind1,2)/U_C(I)%val(ind1,1)))*((U_C(I)%VAL(1,4)/U_C(I)%VAL(1,1))-(U_C(I)%val(ind1,4)/U_C(I)%val(ind1,1)))))*DT/Tz1))
      U_C(I)%RMS(6)=(((U_C(I)%RMS(6))*((Tz1-DT)/(Tz1)))+&
        (((((U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1))-(U_C(I)%val(ind1,3)/U_C(I)%val(ind1,1)))*((U_C(I)%VAL(1,4)/U_C(I)%VAL(1,1))-(U_C(I)%val(ind1,4)/U_C(I)%val(ind1,1)))))*DT/Tz1))
      IF ((PASSIVESCALAR.GT.0))THEN
        U_C(I)%RMS(7)=SQRT(abs(((U_C(I)%RMS(7)**2)*((Tz1-DT)/(Tz1)))+(((U_CT(I)%VAL(1,TURBULENCEEQUATIONS+1)&
          -U_CT(I)%val(ind1,TURBULENCEEQUATIONS+1))**2)*DT/Tz1)))
      END IF
   
    END DO
    !$OMP END DO
  ELSE
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%val(ind1,:)=ZERO;U_C(I)%RMS=zero
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
        U_CT(I)%val(ind1,:)=ZERO
      END IF
    END DO
    !$OMP END DO
  END IF

ELSE
  IF (t.gt.0.0)THEN
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%val(ind1,:)=(((Tz1-dt)/(Tz1))*U_C(I)%val(ind1,:))+((dt*U_C(I)%VAL(1,:))/Tz1)
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
        DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
          U_CT(I)%val(ind1,NVAR)=(((Tz1-dt)/(Tz1))*U_CT(I)%val(ind1,NVAR))+((dt*U_CT(I)%VAL(1,NVAR))/(Tz1*U_C(I)%val(ind1,1)))
        END DO
      END IF
      !U,V,UV,PS
      U_C(I)%RMS(1)=SQRT(abs(((U_C(I)%RMS(1)**2)*((Tz1-dt)/(Tz1)))+(((U_C(I)%VAL(1,2)-U_C(I)%val(ind1,2))**2)*dt/Tz1)))/U_C(I)%val(ind1,1)
      U_C(I)%RMS(2)=SQRT(abs(((U_C(I)%RMS(2)**2)*((Tz1-dt)/(Tz1)))+(((U_C(I)%VAL(1,3)-U_C(I)%val(ind1,3))**2)*dt/Tz1)))/U_C(I)%val(ind1,1)
    
      U_C(I)%RMS(3)=(((U_C(I)%RMS(4))*((Tz1-dt)/(Tz1)))+&
      ((((U_C(I)%VAL(1,2)-U_C(I)%val(ind1,2))*(U_C(I)%VAL(1,3)-U_C(I)%val(ind1,3))))*dt/Tz1))/U_C(I)%val(ind1,1)
    
      IF ((PASSIVESCALAR.GT.0))THEN
        U_C(I)%RMS(4)=SQRT(abs(((U_C(I)%RMS(4)**2)*((Tz1-dt)/(Tz1)))+(((U_CT(I)%VAL(1,TURBULENCEEQUATIONS+1)&
        -U_CT(I)%val(ind1,TURBULENCEEQUATIONS+1))**2)*dt/Tz1)))/U_C(I)%val(ind1,1)
      END IF
    
    END DO
    !$OMP END DO
  ELSE
    !$OMP DO
    DO I=1,KMAXE
      U_C(I)%val(ind1,:)=ZERO
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
        U_CT(I)%val(ind1,:)=ZERO
      END IF
    END DO
    !$OMP END DO
  END IF
END IF

IF (OUTSURF.EQ.1)THEN
  CALL EXCHANGE_HIGHER_AV(N)
  CALL AVERAGE_STRESSES(N)
END IF
   
END SUBROUTINE AVERAGING_T


! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!SUBROUTINE CALLED TO ADVANCE SOLUTION BY ONE TIME STEP SIZE!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE TIME_MARCHING(N)
!> @brief
!> TIME MARCHING SUBROUTINE 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(1:5)::DUMMYOUT,DUMMYIN
INTEGER::I,KMAXE,TTIME
real::dtiv
REAL::CPUT1,CPUT2,CPUT3,CPUT4,CPUT5,CPUT6,CPUT8,timec3,TIMEC1,TIMEC4,TIMEC8,TOTV1,TOTV2,DUMEtg1,DUMEtg2,TOTK,TZX1,TZX2,resolx,totens,totens1,totens2,totensx,totensx1,totensx2
      kill=0
      T=RES_TIME
      resolx=0.01
      iscoun=1
      kmaxe=XMPIELRANK(n)
      EVERY_TIME=((IDNINT(T/output_freq)) * output_freq)+output_freq

      TOTV1=0.0

!$OMP BARRIER
!$OMP MASTER 
    IF (INITCOND.eq.95)THEN                    
      CALL CHECKPOINTv3(N)
    end if
    CPUT1=CPUX1(1)
    CPUT4=CPUX1(1)
    CPUT5=CPUX1(1)
    CPUT8=CPUX1(1)
!$OMP END MASTER 
!$OMP BARRIER
      	      			
IT=RESTART
if (dg.eq.1) call SOL_INTEG_DG_init(N)
      
!$OMP BARRIER
!$OMP MASTER 
  if (tecplot.lt.5)then
      CALL GRID_WRITE
      IF (outsurf.eq.1)THEN
          CALL SURF_WRITE
      END IF
  end if
  IF ((Average_restart.eq.0).and.(averaging.eq.1)) THEN
      Tz1=0.0
  ELSE
      tz1=t
  END IF 
!$OMP END MASTER 
!$OMP BARRIER
      
!$OMP BARRIER
!$OMP MASTER
	CALL VOLUME_SOLUTION_WRITE
	IF (OUTSURF.EQ.1)THEN
	    CALL surface_SOLUTION_WRITE
	END IF
!$OMP END MASTER
!$OMP BARRIER
      
if ((it.eq.0).and.(initcond.eq.95))then
    call EXCHANGE_HIGHER(N)
    call ARBITRARY_ORDER(N)
    call ENSTROPHY_CALC(N)
end if
      
DO 
            
		CALL CALCULATE_CFL(N)
		    
		IF (RUNGEKUTTA.GE.5) CALL CALCULATE_CFLL(N)
		     
    IF (DG.EQ.1)THEN
        DO I=1,KMAXE
          ielem(n,i)%condition=0
          IELEM(N,I)%TROUBLED=0
        END DO
    END IF

		!$OMP BARRIER
		!$OMP MASTER
    DUMMYOUT(1)=DT
    CPUT2=MPI_WTIME()
    TIMEC8=CPUT2-CPUT8
    TIMEC1=CPUT2-CPUT1
    dummyout(2)=TIMEC1
    DUMMYIN=0.0
    timec3=cput2-cput4
    dummyout(3)=TIMEC3
    TIMEC4=CPUT2-CPUT5
    dummyout(4)=TIMEC4
    DUMMYOUT(5)=TIMEC8
			
    CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,5,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    dtiv=DUMMYIN(1)
    DT=DUMMYIN(1)
    TIMEC1=DUMMYIN(2)
    TIMEC3=DUMMYIN(3)
    TIMEC4=DUMMYIN(4)
    TIMEC8=DUMMYIN(5)
    IF (N.EQ.0)THEN
        OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
        WRITE(63,*)DT,it,"TIME STEP SIZE",T
        CLOSE(63)
    END IF
					
    IF (INITCOND.eq.95)THEN
        TOTK=0;TOTENS=0;totensx=0.0d0
        DO I=1,xmpielrank(n)

            TOTK=TOTK+IELEM(N,I)%TOTVOLUME*U_C(I)%VAL(1,1)*(1.0/2.0)*&
                (((U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1))**2)+((U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1))**2)+((U_C(I)%VAL(1,4)/U_C(I)%VAL(1,1))**2))

            if (BOUNDTYPE.eq.1)then
                TOTENS=TOTENS+(IELEM(N,I)%TOTVOLUME*U_C(I)%VAL(1,1)*(1.0/2.0)*&
                IELEM(N,I)%VORTEX(2))
            else
                TOTENS=TOTENS+(IELEM(N,I)%VORTEX(2))
                TOTENSx=TOTENSx+(IELEM(N,I)%VORTEX(3))
            end if
        END DO

        DUMEtg1=TOTK
        DUMEtg2=0.0
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        TOTK=DUMEtg2
        DUMEtg1=TOTENS
        DUMEtg2=0.0
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        TOTENS=DUMEtg2
        DUMEtg1=TOTENSx
        DUMEtg2=0.0
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        TOTENSx=DUMEtg2
        IF (N.EQ.0)THEN
            TOTV1=TOTK/((2.0*PI)**3)
            TOTENS1=TOTENS/(((2.0*PI)**3))
            TOTENSx1=TOTENSx/((2.0*PI)**3)
            IF (it.eq.0)THEN
                TAYLOR=TOTK
                TAYLOR_ENS=TOTENS
                TAYLOR_ENSx=TOTENSx
            END IF
        END IF

		END IF
			
    IF (rungekutta.GE.11)THEN
        dt=timestep
        IF (INITCOND.eq.95)THEN 
            DT=MIN(DT,OUT_TIME-T,EVERY_TIME-T)
        ELSE
            DT=MIN(DT,OUT_TIME-T,EVERY_TIME-T)
        END IF
    else
        IF (INITCOND.eq.95)THEN
            DT=MIN(DT,OUT_TIME-T,EVERY_TIME-T)
        ELSE
            DT=MIN(DT,OUT_TIME-T,EVERY_TIME-T)
        END IF
    end if

    IF (DG.EQ.1)then
        if (filtering.gt.0)then
            call filter(n)
        end if
    end if
			
		!$OMP END MASTER 
		!$OMP BARRIER	
			
		SELECT CASE(RUNGEKUTTA)
			
			CASE(1)
			  CALL RUNGE_KUTTA1(N)
			
			CASE(2)
			  CALL RUNGE_KUTTA2(N)
			
			CASE(3)
			
			  IF (MOOD.EQ.1)THEN
			      CALL RUNGE_KUTTA3_MOOD(N)
			  ELSE
			      CALL RUNGE_KUTTA3(N)
			  END IF
			
			CASE(4)
			  CALL RUNGE_KUTTA4(N)
			
			CASE(5)
			  CALL RUNGE_KUTTA5(N)
			
			case(10)
			  CALL IMPLICIT_TIMEs(N)
			
			case(11)
			  CALL dual_TIME(N)
			
			case(12)
			  CALL dual_TIME_EX(N)
			
		END SELECT
			
		if (dg.eq.1) call SOL_INTEG_DG(N)

		!$OMP BARRIER
		!$OMP MASTER
			
		IF (rungekutta.GE.11)THEN
 			  T=T+(DT)
 			  Tz1=Tz1+(DT)

		ELSE
        T=T+DT
			  tz1=tz1+DT
		END IF
			  
		IF (DG.EQ.1)THEN
        IF (CODE_PROFILE.ne.102)THEN
            IF ( mod(it, 100) .eq. 0) THEN
                CALL TROUBLED_HISTORY
            END IF
        end if
        IF ( filtering .eq. 1) THEN
            CALL FILTERED_HISTORY
        END IF
    end if

    ! IF ( mod(it, 100) .eq. 0) THEN
    !     CALL REDUCED_HISTORY
    ! END IF

		IF (INITCOND.eq.95)THEN                    
 				TOTK=0; TOTENS=0.0; totensx=0.0d0
 				DO I=1,xmpielrank(n)
 				       
            TOTK=TOTK+IELEM(N,I)%TOTVOLUME*U_C(I)%VAL(1,1)*(1.0/2.0)*&
                (((U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1))**2)+((U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1))**2)+((U_C(I)%VAL(1,4)/U_C(I)%VAL(1,1))**2))

            if (BOUNDTYPE.eq.1)then
                TOTENS=TOTENS+(IELEM(N,I)%TOTVOLUME*U_C(I)%VAL(1,1)*(1.0/2.0)*&
                IELEM(N,I)%VORTEX(2))
            else
                TOTENS=TOTENS+(IELEM(N,I)%VORTEX(2))
                TOTENSx=TOTENSx+(IELEM(N,I)%VORTEX(3))
            end if
				END DO
 				
 				DUMEtg1=TOTK
 				DUMEtg2=0.0
 				CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 				CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
 				TOTK=DUMEtg2

 				DUMEtg1=TOTENS
 				DUMEtg2=0.0
 				CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 				CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
 				TOTENS=DUMEtg2

 				DUMEtg1=TOTENSx
 				DUMEtg2=0.0
 				CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 				CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
 				TOTENSx=DUMEtg2

 				IF (N.EQ.0)THEN
            TOTV2=TOTK/((2.0*PI)**3)
            TOTENS2=TOTENS/((2.0*PI)**3)
            TOTENSx2=TOTENSx/((2.0*PI)**3)
            IF (it.eq.0)THEN
                TAYLOR=TOTK
                TAYLOR_ENS=TOTENS
                TAYLOR_ENSx=TOTENSx
            END IF

            IF (IT.EQ.0)THEN
                OPEN(73,FILE='ENERGY.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',POSITION='APPEND')
            ELSE
                OPEN(73,FILE='ENERGY.dat',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
            END IF
            IF (DG.EQ.1)THEN
                WRITE(73,'(E14.7,1X,E14.7,1X,E14.7)')T,TOTK/TAYLOR,-(TOTV2-TOTV1)/DT
            ELSE
                if (boundtype.eq.1)then
                    WRITE(73,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,TOTK/TAYLOR,-(TOTV2-TOTV1)/DT,TOTENS/TAYLOR_ENS
                else
                    WRITE(73,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,TOTK/TAYLOR,-(TOTV2-TOTV1)/DT,TOTENS,TOTENSx
                end if
            END IF
            CLOSE(73)
				END IF
 				
 				CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
 				
 				IF (ADDA.EQ.1)THEN
            TOTK=0
            DO I=1,xmpielrank(n)
                TOTK=TOTK+IELEM(N,I)%ER
            END DO
            DUMEtg1=TOTK
            DUMEtg2=0.0
            CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
            CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

            IF (N.EQ.0)THEN
                TOTK=DUMEtg2/IMAXE
                IF (IT.EQ.0)THEN
                    OPEN(123,FILE='ER.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',POSITION='APPEND')
                ELSE
                    OPEN(123,FILE='ER.dat',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
                END IF
                WRITE(123,*)T,TOTK
                CLOSE(123)
            end if

        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)

    ! END IF
    END IF
           
    ! IF ((initcond.eq.405).or.(initcond.eq.422).or.(initcond.eq.411).or.(initcond.eq.157))THEN
    !     IF ( mod(it, 1) .eq. 0)THEN
    !         CALL TRAJECTORIES
    ! 	  END IF
    !	END IF
    
    
    !$OMP END MASTER 
    !$OMP BARRIER
			
    IF ( mod(it, IForce) .eq. 0) THEN
        IF (OUTSURF.EQ.1) THEN   
            CALL forces
        END IF
    END IF
			
    IF ((rungekutta.ge.5).and.(rungekutta.lt.11))THEN
        IF ( mod(it, residualfreq) .eq. 0) THEN
            CALL RESIDUAL_COMPUTE
        END IF
    END IF

    !$OMP MASTER
    IF (NPROBES.GT.0) CALL PROBING
        
    IF (TIMEC1.GE.IEVERY)THEN
        CALL VOLUME_SOLUTION_WRITE
        IF (outsurf.eq.1)THEN
            CALL surface_SOLUTION_WRITE
        END IF
        CPUT1=MPI_WTIME()
    END IF
    
    IF (INITCOND.eq.95)THEN           
        if (abs(T - ((IDNINT(T/output_freq)) * output_freq)).le.tolsmall) then
            CALL VOLUME_SOLUTION_WRITE
            if (outsurf.eq.1)then
                call surface_SOLUTION_WRITE
            end if
                IF (INITCOND.eq.95)THEN                    
            CALL CHECKPOINTv4(N)
            END IF
            EVERY_TIME=EVERY_TIME+output_freq
        END IF
    ELSE
        IF (CODE_PROFILE.EQ.-1)THEN
            if (abs(T - ((IDNINT(T/output_freq)) * output_freq)).le.tolsmall) then
                CALL VOLUME_SOLUTION_WRITE
                if (outsurf.eq.1)then
                    call surface_SOLUTION_WRITE
                end if
                EVERY_TIME=EVERY_TIME+output_freq
            END IF
        END IF
    END IF
		
    IF (TIMEC8.GE.IEVERYAV)THEN
        IF (AVERAGING.EQ.1)THEN
            CALL VOLUME_SOLUTION_WRITE_av
            IF (outsurf.eq.1)THEN
                CALL surface_SOLUTION_WRITE_av
            END IF
        END IF
        CPUT8=MPI_WTIME()
    END IF
			
    IF (TIMEC4.GE.IEVERY2)THEN
        CALL CHECKPOINTING
        IF (AVERAGING.EQ.1)THEN
            CALL CHECKPOINTING_av
        END IF
        CPUT5=MPI_WTIME()
    END IF
			  
    !$OMP END MASTER 
    !$OMP BARRIER
			
		!$OMP MASTER
		IT=IT+1
			
		IF ((IT.EQ.NTMAX).OR.(TIMEC3.GE.WALLC).OR.(DTiv.GT.OUT_TIME))THEN
			  KILL=1
		END IF
			
		IF ((rungekutta.lt.5).or.(rungekutta.GE.11))THEN
        IF ((T.GE.OUT_TIME).OR.(DTiv.GT.OUT_TIME))THEN
            KILL=1
        END IF
		END IF
    !$OMP END MASTER 
    !$OMP BARRIER

		!$OMP MASTER
		IF (kill.eq.1)THEN
			
			  CALL VOLUME_SOLUTION_WRITE
			  IF (outsurf.eq.1)THEN
			      CALL surface_SOLUTION_WRITE
			  END IF
			  CALL CHECKPOINTING
			  IF (AVERAGING.EQ.1)THEN
			      CALL VOLUME_SOLUTION_WRITE_av
				    IF (outsurf.eq.1)THEN
				        CALL surface_SOLUTION_WRITE_av
				    END IF	    
				    CALL CHECKPOINTING_av
			  END IF
		END IF
			
    !$OMP END MASTER
    !$OMP BARRIER
			
		IF (kill.eq.1)THEN
		  	IF (itestcase.le.3)THEN  
			      CALL CALCULATE_ERROR(n)
			  END IF			  
			  return
		END IF
END do
		      
END SUBROUTINE TIME_MARCHING





SUBROUTINE TIME_MARCHING2(N)
!> @brief
!> TIME MARCHING SUBROUTINE 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(1:5)::DUMMYOUT,DUMMYIN
INTEGER::I,KMAXE
REAL::CPUT1,CPUT2,CPUT3,CPUT4,CPUT5,CPUT6,CPUT8,timec3,TIMEC1,TIMEC4,TIMEC8,TOTV1,TOTV2,DUMEtg1,DUMEtg2,TOTK
real::dtiv,flort
kmaxe=XMPIELRANK(n)
kill=0
T=res_time
iscoun=1

EVERY_TIME=((IDNINT(T/output_freq)) * output_freq)+output_freq


!$OMP MASTER
CPUT1=CPUX1(1)
CPUT4=CPUX1(1)
CPUT5=CPUX1(1)
CPUT8=CPUX1(1)
!$OMP END MASTER
!$OMP BARRIER


IT=RESTART
if (dg.eq.1) call SOL_INTEG_DG_init(N)
!$OMP BARRIER
!$OMP MASTER
if (tecplot.lt.5)then
    CALL GRID_WRITE
end if


CALL VOLUME_SOLUTION_WRITE
IF (outsurf.eq.1)THEN
    CALL SURF_WRITE
END IF

IF ((Average_restart.eq.0).and.(averaging.eq.1)) THEN
    Tz1=0.0
ELSE
    tz1=t
END IF
!$OMP END MASTER
!$OMP BARRIER

DO
    CALL CALCULATE_CFL2D(N)
    IF (RUNGEKUTTA.GE.5) CALL CALCULATE_CFLL2d(N)

    IF (DG.EQ.1)THEN
        DO I=1,KMAXE
            ielem(n,i)%condition=0
            IELEM(N,I)%TROUBLED=0
        END DO
    END IF


    !$OMP MASTER
    DUMMYOUT(1)=DT
    CPUT2=MPI_WTIME()
    TIMEC8=CPUT2-CPUT8
    TIMEC1=CPUT2-CPUT1
    DUMMYOUT(2)=TIMEC1
    DUMMYIN=0.0d0
    TIMEC3=CPUT2-CPUT4
    DUMMYOUT(3)=TIMEC3
    TIMEC4=CPUT2-CPUT5
    DUMMYOUT(4)=TIMEC4
    DUMMYOUT(5)=TIMEC8

    CALL MPI_ALLREDUCE(DUMMYOUT,DUMMYIN,5,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
    DTIV=DUMMYIN(1)
    DT=DUMMYIN(1)
    TIMEC1=DUMMYIN(2)
    TIMEC3=DUMMYIN(3)
    TIMEC4=DUMMYIN(4)
    TIMEC8=DUMMYIN(5)
    IF (N.EQ.0)THEN
        OPEN(63,FILE='history.txt',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
        WRITE(63,*)DT,it,"TIME STEP SIZE",T
        CLOSE(63)
    END IF



    IF (INITCOND.eq.95)THEN
        TOTK=0
        DO I=1,KMAXE
            TOTK=TOTK+IELEM(N,I)%TOTVOLUME*(1.0/2.0)*&
                (((U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1))**2)+((U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1))**2))
        END DO

        DUMEtg1=TOTK
        DUMEtg2=0.0
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
        CALL MPI_ALLREDUCE(DUMEtg1,DUMEtg2,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
        TOTK=DUMEtg2
        IF (N.EQ.0)THEN
            ! TOTV2=TOTK/((2.0*PI)**3)
            ! IF (it.eq.0)THEN
            ! 		TAYLOR=TOTK
            ! END IF
            IF (IT.EQ.0)THEN
                OPEN(73,FILE='ENERGY.dat',FORM='FORMATTED',STATUS='NEW',ACTION='WRITE',POSITION='APPEND')
            ELSE
                OPEN(73,FILE='ENERGY.dat',FORM='FORMATTED',STATUS='old',ACTION='WRITE',POSITION='APPEND')
            END IF
            WRITE(73,*)T,TOTK
            CLOSE(73)
        END IF

        CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
    END IF

    ! IF ((MULTISPECIES.EQ.1))THEN
    !     IF ((initcond.eq.405).or.(initcond.eq.411))THEN
    !         IF ( mod(it, 20) .eq. 0)THEN
    !             CALL TRAJECTORIES
    !         END IF
    !     END IF
    ! END IF

    IF (rungekutta.GE.11)THEN
        dt=timestep
        DT=MIN(DT,OUT_TIME-T,EVERY_TIME-T)
    ELSE
        DT=MIN(DT,OUT_TIME-T,EVERY_TIME-T)
    END IF


    !$OMP END MASTER
    !$OMP BARRIER

    SELECT CASE(RUNGEKUTTA)

      CASE(1)
        CALL RUNGE_KUTTA1_2d(N)

      CASE(2)
        CALL RUNGE_KUTTA2_2d(N)

      CASE(3)
        IF (MOOD.EQ.1)THEN
            CALL RUNGE_KUTTA3_2D_MOOD(N)
        ELSE
            CALL RUNGE_KUTTA3_2D(N)
        END IF

      CASE(4)
        CALL RUNGE_KUTTA4_2D(N)

      CASE(5)
        CALL RUNGE_KUTTA5_2D(N)

      CASE(10)
        CALL IMPLICIT_TIMEs_2d(N)

      CASE(11)
        CALL dual_TIME_2d(N)

      CASE(12)
        CALL dual_TIME_EX_2D(N)

    END SELECT

    if (dg.eq.1)call SOL_INTEG_DG(N)

! Increment time

    !$OMP BARRIER
    !$OMP MASTER
    IF (rungekutta.GE.11)THEN
        T=T+(DT)
        Tz1=Tz1+(DT)
    ELSE
        T=T+DT
        tz1=tz1+DT
    END IF

    IF (DG.EQ.1)THEN
        IF (CODE_PROFILE.ne.102)THEN
            IF ( mod(it, 100) .eq. 0) THEN
                CALL TROUBLED_HISTORY
            END IF
        END IF
    END IF

    IF (mood.gt.0)THEN
        CALL TROUBLED_HISTORY
    end if

! Write output

    !$OMP END MASTER
    !$OMP BARRIER
    IF ( mod(it, IForce) .eq. 0) THEN
        IF (OUTSURF.EQ.1) THEN
            CALL forces
        END IF
    END IF

    IF ((rungekutta.ge.5).and.(rungekutta.lt.11))THEN
        IF ( mod(it, residualfreq) .eq. 0) THEN
            CALL RESIDUAL_COMPUTE
        END IF
    END IF

    !$OMP MASTER
    IF (NPROBES.GT.0) CALL PROBING2D

    IF (TIMEC1.GE.IEVERY)THEN
        CALL VOLUME_SOLUTION_WRITE
        IF (outsurf.eq.1)THEN
            CALL surface_SOLUTION_WRITE
        END IF
        CPUT1=MPI_WTIME()
    END IF

    IF (TIMEC8.GE.IEVERYAV)THEN
        IF (AVERAGING.EQ.1)THEN
            CALL VOLUME_SOLUTION_WRITE_av
            IF (outsurf.eq.1) THEN
                CALL surface_SOLUTION_WRITE_av
            END IF
        END IF
        CPUT8=MPI_WTIME()
    END IF

    IF (CODE_PROFILE.EQ.-1)THEN
			  if (abs(T - ((IDNINT(T/output_freq)) * output_freq)).le.tolsmall) then

            CALL VOLUME_SOLUTION_WRITE
            if (outsurf.eq.1)then
                call surface_SOLUTION_WRITE
            end if
            EVERY_TIME=EVERY_TIME+output_freq
        END IF
    END IF

! Check end condition

    !$OMP END MASTER
    !$OMP BARRIER

    !$OMP MASTER
    IT=IT+1

    IF ((IT.EQ.NTMAX).OR.(TIMEC3.GE.WALLC).OR.(DTiv.GT.OUT_TIME))THEN
        KILL=1
    END IF

    IF ((rungekutta.lt.5).or.(rungekutta.GE.11))THEN
        IF ((T.GE.OUT_TIME).OR.(DTiv.GT.OUT_TIME))THEN
            KILL=1
        END IF
    END IF
    !$OMP END MASTER
    !$OMP BARRIER

    !$OMP MASTER
    IF (kill.eq.1)THEN

        CALL VOLUME_SOLUTION_WRITE
        IF (outsurf.eq.1)THEN
            CALL surface_SOLUTION_WRITE
        END IF
        CALL CHECKPOINTING
        IF (AVERAGING.EQ.1)THEN
            CALL VOLUME_SOLUTION_WRITE_av
            IF (outsurf.eq.1)THEN
                CALL surface_SOLUTION_WRITE_av
            END IF
            CALL CHECKPOINTING_av
        END IF
    END IF

    !$OMP END MASTER
    !$OMP BARRIER

    IF (kill.eq.1)THEN
        IF (itestcase.le.3)THEN
            CALL CALCULATE_ERROR(n)
        END IF

        return
    END IF

END DO

END SUBROUTINE TIME_MARCHING2





!---------------------------------------------------------------------------------------------!
END MODULE ADVANCE
