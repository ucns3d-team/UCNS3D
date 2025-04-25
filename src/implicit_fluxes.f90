module implicit_fluxes
use declaration
USE LIBRARY
USE TRANSFORM
USE LOCAL
USE RIEMANN
USE FLOW_OPERATIONS
use source
IMPLICIT NONE

contains
SUBROUTINE CALCULATE_JACOBIAN(N)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 3D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,N_NODE,IBFC
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB,tempxx,VISCOTS
	REAL,DIMENSION(NOF_variables,NOF_variables)::IDENTITY1
	real,dimension(NOF_variables,NOF_variables)::convj,diffj
	INTEGER::ICONSIDERED, FACEX, POINTX,igoflux
	INTEGER::B_CODE,srf
	REAL::ANGLE1,ANGLE2,NX,NY,NZ
	real,dimension(1:nof_variables)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
	real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
	real,dimension(1:nof_Variables)::RIGHTv
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:8,1:DIMENSIONA)::NODES_LIST
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	REAL,DIMENSION(1:4)::viscl,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
	real::MP_PINFL,gammal
    real::MP_PINFR,gammaR
   REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::EIGVL

	KMAXE=XMPIELRANK(N)
	IDENTITY1(:,:)=ZERO
	DO L=1,NOF_VARIABLES
	IDENTITY1(L,L)=1.0D0
	END DO
		

	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
		IMPDIAG(i,:,:)=zero
		IMPOFF(i,:,:,:)=zero
		if (turbulence.eq.1)then
		impdiagt(i,:)=zero
		IMPOFFt(i,:,:)=zero
		end if
	
		    
		    
        IF (MRF.EQ.1)THEN
            SRF=ILOCAL_RECON3(I)%MRF
        END IF 
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=(COS(ANGLE1)*SIN(ANGLE2))
				  NY=(SIN(ANGLE1)*SIN(ANGLE2))
				  NZ=(COS(ANGLE2))
				  mul1=IELEM(N,I)%SURF(L)
				  B_CODE=0
				  CLEFT(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				  CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				     				      
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)!left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)

					END IF
			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  
						  CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						   CALL ROTATEB(N,ClefT,Cleft_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						  ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  IF (ILOCAL_RECON3(i)%MRF.EQ.1)THEN
                                !RETRIEVE THE ROTATIONAL VELOCITY (AT THE GAUSSIAN POINT JUST FOR SECOND ORDER)
                                SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,1,1:3)
                                CALL ROTATEF(N,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
                                !CALCULATE THE NEW EIGENVALUE FOR ROTATING REFERENCE FRAME
                                ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1)-SRF_SPEEDROT(2))
                                ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1)-SRF_SPEEDROT(2))
                            END IF
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 
							  
							  
							  
							  EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      if (turbulence.eq.1)then
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
							      else
							      viscots=0.5*((viscl(1)+viscl(2)))
							      
							      end if
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
                                                
						  
						  
						  END IF
						  

		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		   IMPDIAG(i,:,:)=0.0
		IMPOFF(i,:,:,:)=0.0
		if (turbulence.eq.1)then
		impdiagt(i,:)=0.0
		IMPOFFt(i,:,:)=0.0
		end if
		IF (MRF.EQ.1)THEN
            SRF=ILOCAL_RECON3(I)%MRF
        END IF     
		    DO L=1,IELEM(N,I)%IFCA
				      B_CODE=0
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
 				  mul1=IELEM(N,I)%SURF(L)
				      B_CODE=0
				      CLEFT(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
                 IF (ILOCAL_RECON3(i)%MRF.EQ.1)THEN
                    !RETRIEVE ROTATIONAL VELOCITY IN CASE OF ROTATING REFERENCE FRAME TO CALCULATE THE CORRECT VALUE OF THE BOUNDARY CONDITION
                    SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,1,1:3)
                    CALL ROTATEF(N,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
                END IF
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN MY CPU
								  CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
                                    IF(PER_ROT.EQ.1)THEN
                                        CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                    END IF
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_innerx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)

								   if (ielem(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
                                            N_NODE=4
                                    else
                                            N_NODE=3
                                    end if

								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
				  				   
				  				  				  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if  ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								IF(PER_ROT.EQ.1)THEN
                                    CRIGHT(2:4)=ROTATE_PER_1(CRIGHT(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
								END IF 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							
								  IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				       CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
				      
				      
				      IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						   CALL ROTATEB(N,ClefT,Cleft_ROT,ANGLE1,ANGLE2)
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
                            ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
                            ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
                        IF (ILOCAL_RECON3(i)%MRF.EQ.1)THEN
                            ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1)-SRF_SPEEDROT(2))
                            ASOUND2=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(Cright_ROT(2)/Cright_ROT(1)-SRF_SPEEDROT(2))
                        END IF

						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							    EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    
							    
							  
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      
							  
							  
							  
							  if (turbulence.eq.1)then
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
							      else
							      viscots=0.5*((viscl(1)+viscl(2)))
							      
							      end if
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF

						
				   
				  
		    END DO
	END DO
	!$OMP END DO

    !ADD THE CONTRIBUTION OF THE SOURCE TERM TO THE JACOBIAN OF THE DIAGONAL MATRIX
        IF (SRFg.EQ.1) THEN
        !$OMP DO
            DO I=1,KMAXE
                IMPDIAG(i,2,3)=-SRF_VELOCITY(3)*ielem(n,I)%totvolume
                IMPDIAG(i,2,4)=SRF_VELOCITY(2)*ielem(n,I)%totvolume
                IMPDIAG(i,3,2)=SRF_VELOCITY(3)*ielem(n,I)%totvolume
                IMPDIAG(i,3,4)=-SRF_VELOCITY(1)*ielem(n,I)%totvolume
                IMPDIAG(i,4,2)=-SRF_VELOCITY(2)*ielem(n,I)%totvolume
                IMPDIAG(i,4,3)=SRF_VELOCITY(1)*ielem(n,I)%totvolume
            END DO
        !$OMP END DO
        END IF	
        IF (MRF.EQ.1) THEN
        !$OMP DO
            DO I=1,KMAXE
				SRF=ILOCAL_RECON3(I)%MRF
                IF (ILOCAL_RECON3(i)%MRF.EQ.1)THEN
                    IMPDIAG(i,2,3)=-ILOCAL_RECON3(I)%MRF_VELOCITY(3)*ielem(n,I)%totvolume
                    IMPDIAG(i,2,4)=ILOCAL_RECON3(I)%MRF_VELOCITY(2)*ielem(n,I)%totvolume
                    IMPDIAG(i,3,2)=ILOCAL_RECON3(I)%MRF_VELOCITY(3)*ielem(n,I)%totvolume
                    IMPDIAG(i,3,4)=-ILOCAL_RECON3(I)%MRF_VELOCITY(1)*ielem(n,I)%totvolume
                    IMPDIAG(i,4,2)=-ILOCAL_RECON3(I)%MRF_VELOCITY(2)*ielem(n,I)%totvolume
                    IMPDIAG(i,4,3)=ILOCAL_RECON3(I)%MRF_VELOCITY(1)*ielem(n,I)%totvolume
                END IF
				SRF=0
            END DO
        !$OMP END DO
        END IF
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		!$OMP DO
		do i=1,kmaxe
            DO L=1,NOF_VARIABLES
		    IMPDIAG(i,L,L)=IMPDIAG(i,L,L)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    END DO
		end do
		!$OMP END DO
	  ELSE
	!$OMP DO
! 	
	  do i=1,kmaxe
        DO L=1,NOF_VARIABLES
	    IMPDIAG(I,L,L)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,L,L))
	    END DO
	end do
	!$OMP END DO
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 if (turbulence.eq.1)CALL SOURCES_derivatives_COMPUTATION(N)
if (rungekutta.eq.10)then
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+(ielem(n,I)%totvolume/(ielem(n,I)%dtl))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
    end do
    end if
end do
!$OMP END DO
else
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations

    IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(i,nvar))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(i,1))
 end do
    end if
end do
!$OMP END DO

end if
end if
	

END SUBROUTINE CALCULATE_JACOBIAN
	
	

SUBROUTINE CALCULATE_JACOBIAN_2D(N)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,N_NODE,IBFC
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB,tempxx,VISCOTS
	REAL,DIMENSION(NOF_variables,NOF_variables)::IDENTITY1
	real,dimension(NOF_variables,NOF_variables)::convj,diffj
	INTEGER::ICONSIDERED, FACEX, POINTX,igoflux,kas
	INTEGER::B_CODE
	REAL::ANGLE1,ANGLE2,NX,NY,NZ
	real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
	real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
	real,dimension(1:nof_Variables)::RIGHTv
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	REAL,DIMENSION(1:4)::viscl,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR

	real::MP_PINFL,gammal
    real::MP_PINFR,gammaR
     REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::EIGVL



	KMAXE=XMPIELRANK(N)
	IDENTITY1(:,:)=ZERO
	IDENTITY1(1,1)=1.0D0
	IDENTITY1(2,2)=1.0D0
	IDENTITY1(3,3)=1.0D0
	IDENTITY1(4,4)=1.0D0
	
		

	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
		IMPDIAG(i,:,:)=zero
		IMPOFF(i,:,:,:)=zero
		if (turbulence.eq.1)then
		impdiagt(i,:)=zero
		IMPOFFt(i,:,:)=zero
		end if
	
		    
		    
		    
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				  B_CODE=0
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=ANGLE1
				  NY=ANGLE2
				  mul1=IELEM(N,I)%SURF(L)
				  
				  
				    CLEFT(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				   CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				     				      
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)!left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)

					END IF
			
						  CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEb2D(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						 ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
! 						  viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					 
						      
						     VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2D(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      
							  
							  
							  
							  if (turbulence.eq.1)then
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
							      else
							      viscots=0.5*((viscl(1)+viscl(2)))
							      
							      end if
							  
							  
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		   IMPDIAG(i,:,:)=zero
		IMPOFF(i,:,:,:)=zero
		if (turbulence.eq.1)then
		impdiagt(i,:)=zero
		IMPOFFt(i,:,:)=zero
		end if
		    
		    DO L=1,IELEM(N,I)%IFCA
				      mul1=IELEM(N,I)%SURF(L)
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
 				  
				      B_CODE=0
				      CLEFT(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  KAS=1
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2Dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								     CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
				  				    
				  				  	KAS=2			  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      KAS=3
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								KAS=4
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							KAS=5
								  IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				    
						  
						  
						 CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEb2D(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
! 						  viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
! 							   
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							        
							  if (turbulence.eq.1)then
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
							      else
							      viscots=0.5*((viscl(1)+viscl(2)))
							      
							      end if
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(i,1:nof_Variables,1:nof_Variables)=IMPDIAG(i,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(i,l,1:nof_Variables,1:nof_Variables)=IMPOFF(i,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
			
						
				   
				  
		    END DO
	END DO
	!$OMP END DO

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		!$OMP DO
		do i=1,kmaxe
				  
		    IMPDIAG(i,1,1)=IMPDIAG(i,1,1)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,2,2)=IMPDIAG(i,2,2)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,3,3)=IMPDIAG(i,3,3)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,4,4)=IMPDIAG(i,4,4)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    
		end do
		!$OMP END DO
	  ELSE
	!$OMP DO
	  do i=1,kmaxe
	      IMPDIAG(I,1,1)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,1,1))
	      IMPDIAG(I,2,2)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,2,2))
	      IMPDIAG(I,3,3)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,3,3))
	      IMPDIAG(I,4,4)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,4,4))
	    
	end do
	!$OMP END DO
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 
 if (turbulence.eq.1)CALL SOURCES_derivatives_COMPUTATION2D(N)
if (rungekutta.eq.10)then
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
!    
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+((ielem(n,I)%totvolume/(ielem(n,I)%dtl)))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+(ielem(n,I)%totvolume/(ielem(n,I)%dtl))
    end do
    end if
end do
!$OMP END DO
else
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*((1.0D0/(ielem(n,I)%dtl))+(1.5D0/DT))+(IMPDIAGt(i,1))-sht(i,nvar)
!     IMPDIAGt(i,nvar)=(ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+((1.5D0/dt)*IMPDIAGt(i,nvar))))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*((1.0D0/(ielem(n,I)%dtl))+(1.5D0/DT))+(IMPDIAGt(i,1))
!     IMPDIAGt(i,nvar)=(ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+((1.5D0/dt)*IMPDIAGt(i,nvar))))
    end do
    end if
end do
!$OMP END DO

end if
end if
	

END SUBROUTINE CALCULATE_JACOBIAN_2D
	
	
SUBROUTINE CALCULATE_JACOBIANLM(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 3D with low-memory footprint
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,ICONSIDERED
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,N_NODE,IBFC
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB,tempxx,VISCOTS
	REAL,DIMENSION(NOF_variables,NOF_variables)::IDENTITY1
	real,dimension(NOF_variables,NOF_variables)::convj,diffj
	INTEGER::FACEX, POINTX,igoflux
	INTEGER::B_CODE
	REAL::ANGLE1,ANGLE2,NX,NY,NZ
	real,dimension(1:nof_variables)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
	real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
	real,dimension(1:nof_Variables)::RIGHTv
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ

	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
	REAL,DIMENSION(1:4)::viscl,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR

	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	real::MP_PINFL,gammal
    real::MP_PINFR,gammaR
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IMPDIAGT
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::IMPDIAG,IMPOFFt
	REAL,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::IMPOFF
    REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::EIGVL
   REAL,DIMENSION(TURBULENCEEQUATIONS)::SOURCE_T



	IDENTITY1(:,:)=ZERO
	IDENTITY1(1,1)=1.0D0
	IDENTITY1(2,2)=1.0D0
	IDENTITY1(3,3)=1.0D0
	IDENTITY1(4,4)=1.0D0
	IDENTITY1(5,5)=1.0D0
		
	IF (IELEM(N,ICONSIDERED)%INTERIOR.EQ.0)THEN
	
	I=ICONSIDERED
		IMPDIAG(1,:,:)=ZERO
		IMPOFF(1,:,:,:)=ZERO
		if (turbulence.eq.1)then
		impdiagt(1,:)=ZERO
		IMPOFFt(1,:,:)=ZERO
		end if
	
		    
		    
		    
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=(COS(ANGLE1)*SIN(ANGLE2))
				  NY=(SIN(ANGLE1)*SIN(ANGLE2))
				  NZ=(COS(ANGLE2))
				  mul1=IELEM(N,I)%SURF(L)
				  
				    CLEFT(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				   CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				     				      
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)!left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)

					END IF
			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEF(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  
							 EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 
							  
							  
							  
							  EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
		    END DO
	
	ELSE
	I=ICONSIDERED			
		   IMPDIAG(1,:,:)=ZERO
		IMPOFF(1,:,:,:)=ZERO
		if (turbulence.eq.1)then
		impdiagt(1,:)=ZERO
		IMPOFFt(1,:,:)=ZERO
		end if
		    
		    DO L=1,IELEM(N,I)%IFCA
				     mul1=IELEM(N,I)%SURF(L) 
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
 				  
				      B_CODE=0
				      CLEFT(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								 facex=l;
								  CALL coordinates_face_innerx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)

								   if (ielem(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
                                            N_NODE=4
                                    else
                                            N_NODE=3
                                    end if

								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    
								    

								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
				  				   
				  				  				  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							
								  IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				       CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEF(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						  ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 
							  
							  
							  
							  EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
			
						
				   
				  
		    END DO
	END IF
	

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		
		
		    IMPDIAG(1,1,1)=IMPDIAG(1,1,1)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,2,2)=IMPDIAG(1,2,2)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,3,3)=IMPDIAG(1,3,3)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,4,4)=IMPDIAG(1,4,4)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,5,5)=IMPDIAG(1,5,5)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		
	  ELSE
	
! 	    IMPDIAG(1,1,1)=ielem(n,I)%totvolume*( ((dt+1.5d0*ielem(n,I)%dtl)/dt) +(IMPDIAG(1,1,1)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl
! 	    IMPDIAG(1,2,2)=ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt)+(IMPDIAG(1,2,2)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl
! 	    IMPDIAG(1,3,3)=ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt)+(IMPDIAG(1,3,3)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl
! 	    IMPDIAG(1,4,4)=ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt)+(IMPDIAG(1,4,4)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl
! 	    IMPDIAG(1,5,5)=ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt)+(IMPDIAG(1,5,5)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl

	IMPDIAG(1,1,1)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,1,1))
	      IMPDIAG(1,2,2)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,2,2))
	      IMPDIAG(1,3,3)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,3,3))
	      IMPDIAG(1,4,4)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,4,4))
	    IMPDIAG(1,5,5)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,5,5))
	
	
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 CALL SOURCES_derivatives(N,ICONSIDERED,SOURCE_T)
 sht(I,1:turbulenceequations)=(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
if (rungekutta.eq.10)then

    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    IMPDIAGt(1,nvar)=IMPDIAGt(1,nvar)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(1,nvar)=IMPDIAGt(1,nvar)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
    end do
    end if

else

    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    IMPDIAGt(1,nvar)=(ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+((1.5D0/dt)*IMPDIAGt(1,nvar))))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(1,nvar)=(ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+((1.5D0/dt)*IMPDIAGt(1,nvar))))
    end do
    end if


end if
end if
	

END SUBROUTINE CALCULATE_JACOBIANLM
	
	

SUBROUTINE CALCULATE_JACOBIAN_2DLM(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 2D with low-memory footprint
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N,ICONSIDERED
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,N_NODE,IBFC
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB,tempxx,VISCOTS
	REAL,DIMENSION(NOF_variables,NOF_variables)::IDENTITY1
	real,dimension(NOF_variables,NOF_variables)::convj,diffj
	INTEGER::FACEX, POINTX,igoflux
	INTEGER::B_CODE
	REAL::ANGLE1,ANGLE2,NX,NY,NZ
	real,dimension(1:nof_variables)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
	real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
	real,dimension(1:nof_Variables)::RIGHTv
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
	REAL,DIMENSION(1:4)::viscl,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
    REAL,DIMENSION(1:DIMENSIONA)::CORDS
	real::MP_PINFL,gammal
    real::MP_PINFR,gammaR
	REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::IMPDIAGT
	REAL,ALLOCATABLE,DIMENSION(:,:,:),INTENT(INOUT)::IMPDIAG,IMPOFFt
	REAL,ALLOCATABLE,DIMENSION(:,:,:,:),INTENT(INOUT)::IMPOFF
    REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::EIGVL
   REAL,DIMENSION(TURBULENCEEQUATIONS)::SOURCE_T



	IDENTITY1(:,:)=ZERO
	IDENTITY1(1,1)=1.0D0
	IDENTITY1(2,2)=1.0D0
	IDENTITY1(3,3)=1.0D0
	IDENTITY1(4,4)=1.0D0
	
		

	IF (IELEM(N,ICONSIDERED)%INTERIOR.EQ.0)THEN
	
	I=ICONSIDERED
		IMPDIAG(1,:,:)=zero
		IMPOFF(1,:,:,:)=zero
		if (turbulence.eq.1)then
		impdiagt(1,:)=zero
		IMPOFFt(1,:,:)=zero
		end if
	
		    
		    
		    
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=ANGLE1
				  NY=ANGLE2
				  mul1=IELEM(N,I)%SURF(L)
				  
				  
				    CLEFT(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				   CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				     				      
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)!left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)

					END IF

						  CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						   CALL ROTATEb2D(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and so
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							   EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
		    END DO
	else
	

				
		   IMPDIAG(1,:,:)=0.0
		IMPOFF(1,:,:,:)=0.0
		if (turbulence.eq.1)then
		impdiagt(1,:)=0.0
		IMPOFFt(1,:,:)=0.0
		end if
		    
		    DO L=1,IELEM(N,I)%IFCA
				      mul1=IELEM(N,I)%SURF(L)
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
 				  
				      B_CODE=0
				      CLEFT(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;
								  CALL coordinates_face_inner2Dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
				  				   
				  				  				  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							
								  IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				    CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						   CALL ROTATEb2D(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and so
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						 ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(i,nvar)=impdiagt(i,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(i,l,nvar)=impofft(i,l,nvar)-(((OO2*vpp))*mul1)
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  impdiagt(1,nvar)=impdiagt(1,nvar)+(OO2*((vpp))*MUL1)
							  IMPOFFt(1,l,nvar)=impofft(1,l,nvar)-(OO2*((vpp))*MUL1)
							  end do
							  end if
						  END IF
						  ELSE
						  
						  IMPDIAG(1,1:nof_Variables,1:nof_Variables)=IMPDIAG(1,1:nof_Variables,1:nof_Variables)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
						  convj=eigvl
						  IMPOFF(1,l,1:nof_Variables,1:nof_Variables)=IMPOFF(1,l,1:nof_Variables,1:nof_Variables)+(((OO2*CONVJ(1:nof_Variables,1:nof_Variables))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
			
						
				   
				  
		    END DO
	end if

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		
		    IMPDIAG(1,1,1)=IMPDIAG(1,1,1)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,2,2)=IMPDIAG(1,2,2)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,3,3)=IMPDIAG(1,3,3)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(1,4,4)=IMPDIAG(1,4,4)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    
		
	  ELSE
	    IMPDIAG(1,1,1)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,1,1))
	      IMPDIAG(1,2,2)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,2,2))
	      IMPDIAG(1,3,3)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,3,3))
	      IMPDIAG(1,4,4)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(1,4,4))
	    
	    

	
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 CALL SOURCES_derivatives2d(N,ICONSIDERED,SOURCE_T)
 sht(I,1:turbulenceequations)=(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
if (rungekutta.eq.10)then

    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    IMPDIAGt(1,nvar)=IMPDIAGt(1,nvar)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    
    IMPDIAGt(1,nvar)=IMPDIAGt(1,nvar)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
    end do
    end if

else

    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    IMPDIAGt(1,nvar)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(1,1))-sht(i,nvar)
    
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(1,nvar)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(1,1))
!     IMPDIAGt(1,nvar)=(ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+((1.5D0/dt)*IMPDIAGt(1,nvar))))
    end do
    end if


end if
end if
	

END SUBROUTINE CALCULATE_JACOBIAN_2DLM

SUBROUTINE CALCULATE_JACOBIAN_2D_MF(N)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,N_NODE,IBFC,KAS
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB,tempxx,VISCOTS
	REAL,DIMENSION(NOF_variables,NOF_variables)::IDENTITY1
	real,dimension(NOF_variables,NOF_variables)::convj,diffj
	INTEGER::ICONSIDERED, FACEX, POINTX,igoflux
	INTEGER::B_CODE
	REAL::ANGLE1,ANGLE2,NX,NY,NZ
	real,dimension(1:nof_variables)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
	real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
	real,dimension(1:nof_Variables)::RIGHTv
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ

	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	REAL,DIMENSION(1:4)::viscl,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR

	real::MP_PINFL,gammal
    real::MP_PINFR,gammaR
   REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::EIGVL
	KMAXE=XMPIELRANK(N)
	

	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
		IMPDIAG_MF(i)=zero
		IMPOFF_MF(i,:)=zero
        if (turbulence.eq.1)then
		impdiagt(I,1:turbulenceequations+PASSIVESCALAR)=0.0
		IMPOFFt(i,:,:)=zero
		end if
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				  B_CODE=0
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=ANGLE1
				  NY=ANGLE2
				  mul1=IELEM(N,I)%SURF(L)
				  
				  
				    CLEFT(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				   CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				    
! 				     CLEFT(1:nof_variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_variables,L,1)
! 				   CRIGHT(1:nof_variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),1)
				    
				    
				     		
				     		
				     		
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)!left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)

					END IF
			
						  CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEb2D(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						 ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  
! 						  viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
 
                                        IF (TURBULENCE.EQ.1)THEN
                                            IF (TURBULENCEMODEL.EQ.1)THEN
                                            TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
                                            Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
                                            END IF
                                            IF (TURBULENCEMODEL.EQ.2)THEN
                                            EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
                                            EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
                                            EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
                                            EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
                                                
                                                
                                            eddyfr=eddyfl
                                                Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
                                            END IF
                                    
                                            
                                            VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
                        
                ! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
                                            viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
                                            mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
                                            viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
                                            
                                            VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
                                        END IF
						  END IF
						  
						  
						   IMPDIAG_MF(i)=IMPDIAG_MF(i)+(OO2*vpp*MUL1)
						  IMPOFF_MF(i,l)=vpp
						  
						  if (turbulence.eq.1)then
                            IMPDIAGT(i,:)=IMPDIAGT(i,:)+(OO2*vpp*MUL1)
                            IMPOFFt(i,l,:)=impofft(i,l,:)-(OO2*((vpp))*MUL1)
                        end if
						  
						  
						  
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		IMPDIAG_MF(i)=zero
		IMPOFF_MF(i,:)=zero
		  if (turbulence.eq.1)then
		impdiagt(I,1:turbulenceequations+PASSIVESCALAR)=zero
		IMPOFFt(i,:,:)=zero
		end if   
		    DO L=1,IELEM(N,I)%IFCA
				      mul1=IELEM(N,I)%SURF(L)
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
 				  
				      B_CODE=0
				      CLEFT(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
! 				      CLEFT(1:nof_variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_variables,L,1)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
! 								  cRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),1)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  KAS=1
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2Dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
				  				    
				  				  	KAS=2			  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
! 							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),1)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      KAS=3
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							      
! 							      CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
							      
							      
							      
							      
							    END IF
								KAS=4
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							KAS=5
								  IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							      
! 							      CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
							      
							      
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				    
						  
						  
						 CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEb2D(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  
! 						  viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
! 							   
							 
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							        
							  if (turbulence.eq.1)then
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
							      else
							      viscots=0.5*((viscl(1)+viscl(2)))
							      
							      end if
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  
							  end do
							  end if
						  END IF
						  ELSE
						  
						  END IF
						  
						  
						  
						   IMPDIAG_MF(i)=IMPDIAG_MF(i)+(OO2*vpp*MUL1)
						  IMPOFF_MF(i,l)=vpp
                         if (turbulence.eq.1)then
                            IMPDIAGT(i,:)=IMPDIAGT(i,:)+(OO2*vpp*MUL1)
                            
                            IMPOFFt(i,l,:)=impofft(i,l,:)-(OO2*((vpp))*MUL1)
                        end if
						
				   
				  
		    END DO
	END DO
	!$OMP END DO

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		!$OMP DO
		do i=1,kmaxe
				  
		    IMPDIAG_MF(i)=(IMPDIAG_MF(i))+(IELEM(N,I)%TOTVOLUME/ielem(n,I)%dtl)
		    
		end do
		!$OMP END DO
	  ELSE
	!$OMP DO
	  do i=1,kmaxe
            IMPDIAG_MF(i)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG_MF(i))
	end do
	!$OMP END DO
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 
 if (turbulence.eq.1)CALL SOURCES_derivatives_COMPUTATION2D(N)
if (rungekutta.eq.10)then
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
   IMPDIAGT(i,NVAR)=IMPDIAGT(i,NVAR)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar) 
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
     IMPDIAGT(i,NVAR)=IMPDIAGT(i,NVAR)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar) 
    end do
    end if
end do
!$OMP END DO
else
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
!    
     IMPDIAGT(i,NVAR)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(i,nvar))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
     IMPDIAGT(i,NVAR)=IMPDIAGT(i,NVAR)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar) 
    end do
    end if

end do
!$OMP END DO

end if
end if
	

END SUBROUTINE CALCULATE_JACOBIAN_2D_MF



SUBROUTINE CALCULATE_JACOBIAN_3D_MF(N)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,N_NODE,IBFC,KAS
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB,tempxx,VISCOTS
	REAL,DIMENSION(NOF_variables,NOF_variables)::IDENTITY1
	real,dimension(NOF_variables,NOF_variables)::convj,diffj
	INTEGER::ICONSIDERED, FACEX, POINTX,igoflux
	INTEGER::B_CODE
	REAL::ANGLE1,ANGLE2,NX,NY,NZ
	real,dimension(1:nof_variables)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
	real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
	real,dimension(1:nof_Variables)::RIGHTv
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST

	REAL,DIMENSION(1:4)::viscl,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR

	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	real::MP_PINFL,gammal
    real::MP_PINFR,gammaR
   REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::EIGVL
	KMAXE=XMPIELRANK(N)
	

	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
		IMPDIAG_MF(i)=zero
		IMPOFF_MF(i,:)=zero
        if (turbulence.eq.1)then
		impdiagt(I,1:turbulenceequations+PASSIVESCALAR)=0.0
		IMPOFFt(i,:,:)=zero
		end if
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				  B_CODE=0
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=ANGLE1
				  NY=ANGLE2
				  mul1=IELEM(N,I)%SURF(L)
				  
				  
				    CLEFT(1:nof_variables)=U_C(I)%VAL(1,1:nof_variables)
				   CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
				    
! 				     CLEFT(1:nof_variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_variables,L,1)
! 				   CRIGHT(1:nof_variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),1)
				    
				    
				     		
				     		
				     		
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)!left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)

					END IF
			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  
! 						  viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
 
                                        IF (TURBULENCE.EQ.1)THEN
                                            IF (TURBULENCEMODEL.EQ.1)THEN
                                            TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
                                            Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
                                            END IF
                                            IF (TURBULENCEMODEL.EQ.2)THEN
                                            EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
                                            EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
                                                
                                                
                                            eddyfr=eddyfl
                                                Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
                                            END IF
                                    
                                            
                                            VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
                        
                ! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
                                            viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
                                            mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
                                            viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
                                            
                                            VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
                                        END IF
						  END IF
						  
						  
						   IMPDIAG_MF(i)=IMPDIAG_MF(i)+(OO2*vpp*MUL1)
						  IMPOFF_MF(i,l)=vpp
						  
						  if (turbulence.eq.1)then
                            IMPDIAGT(i,:)=IMPDIAGT(i,:)+(OO2*vpp*MUL1)
                            IMPOFFt(i,l,:)=impofft(i,l,:)-(OO2*((vpp))*MUL1)
                        end if
						  
						  
						  
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		IMPDIAG_MF(i)=zero
		IMPOFF_MF(i,:)=zero
		  if (turbulence.eq.1)then
		impdiagt(I,1:turbulenceequations+PASSIVESCALAR)=zero
		IMPOFFt(i,:,:)=zero
		end if   
		    DO L=1,IELEM(N,I)%IFCA
				      mul1=IELEM(N,I)%SURF(L)
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
 				  
				      B_CODE=0
				      CLEFT(1:nof_Variables)=U_C(I)%VAL(1,1:nof_Variables)
! 				      CLEFT(1:nof_variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_variables,L,1)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
! 								  cRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),1)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  KAS=1
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_innerx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)

								   if (ielem(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
                                            N_NODE=4
                                    else
                                            N_NODE=3
                                    end if

								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)
								   
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
				  				    
				  				  	KAS=2			  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
! 							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),1)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      KAS=3
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							      
! 							      CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
							      
							      
							      
							      
							    END IF
								KAS=4
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							KAS=5
								  IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:nof_Variables)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:nof_Variables)
							    ELSE
							     
							      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)
							      
! 							      CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
							      
							      
							    END IF
								
								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									IF (FASTEST.EQ.1)THEN
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL&
							      (IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    ELSE
							     
							      CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
							    END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				    
						  
						  
						 CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:nof_Variables)=CLEFT(1:nof_Variables);RIGHTV(1:nof_Variables)=CRIGHT(1:nof_Variables)						  
						  CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
						  
						ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  
! 						  viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							    EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  
						  
						  
						  if ((turbulence.eq.1).or.(passivescalar.gt.0))then
							  if (turbulence.eq.1)then
							  do nvar=1,turbulenceequations
							  
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
! 							   
							 
							  end do
							  end if
							  if (passivescalar.gt.0)then
							  do nvar=turbulenceequations+1,turbulenceequations+passivescalar
							  VISCL(1)=VISCL(1)/SCHMIDT_LAM
							      VISCL(2)=VISCL(2)/SCHMIDT_LAM
							        
							  if (turbulence.eq.1)then
							      VISCL(3)=VISCL(3)/SCHMIDT_TURB
							      VISCL(4)=VISCL(4)/SCHMIDT_TURB
							  viscots=0.5*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
							      else
							      viscots=0.5*((viscl(1)+viscl(2)))
							      
							      end if
												  
							  VISCOTS=(2.0*VISCOTS)/((cleft(1)+cRIGHT(1))*IELEM(N,I)%DIH(L))
							  vpp=MAX(ASOUND1,ASOUND2)+viscots
							  
							  end do
							  end if
						  END IF
						  ELSE
						  
						  END IF
						  
						  
						  
						   IMPDIAG_MF(i)=IMPDIAG_MF(i)+(OO2*vpp*MUL1)
						  IMPOFF_MF(i,l)=vpp
                         if (turbulence.eq.1)then
                            IMPDIAGT(i,:)=IMPDIAGT(i,:)+(OO2*vpp*MUL1)
                            
                            IMPOFFt(i,l,:)=impofft(i,l,:)-(OO2*((vpp))*MUL1)
                        end if
						
				   
				  
		    END DO
	END DO
	!$OMP END DO

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		!$OMP DO
		do i=1,kmaxe
				  
		    IMPDIAG_MF(i)=(IMPDIAG_MF(i))+(IELEM(N,I)%TOTVOLUME/ielem(n,I)%dtl)
		    
		end do
		!$OMP END DO
	  ELSE
	!$OMP DO
	  do i=1,kmaxe
            IMPDIAG_MF(i)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG_MF(i))
	end do
	!$OMP END DO
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 
 if (turbulence.eq.1)CALL SOURCES_derivatives_COMPUTATION(N)
if (rungekutta.eq.10)then
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
   IMPDIAGT(i,NVAR)=IMPDIAGT(i,NVAR)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar) 
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
     IMPDIAGT(i,NVAR)=IMPDIAGT(i,NVAR)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar) 
    end do
    end if
end do
!$OMP END DO
else
!$OMP DO
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
!    
     IMPDIAGT(i,NVAR)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(i,nvar))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
     IMPDIAGT(i,NVAR)=IMPDIAGT(i,NVAR)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar) 
    end do
    end if

end do
!$OMP END DO

end if
end if
	

END SUBROUTINE CALCULATE_JACOBIAN_3D_MF



END module implicit_fluxes
