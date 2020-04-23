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
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1,DXB
	REAL,DIMENSION(5,5)::IDENTITY1
	real,dimension(5,5)::convj,diffj
	KMAXE=XMPIELRANK(N)
	IDENTITY1(:,:)=ZERO
	IDENTITY1(1,1)=1.0D0
	IDENTITY1(2,2)=1.0D0
	IDENTITY1(3,3)=1.0D0
	IDENTITY1(4,4)=1.0D0
	IDENTITY1(5,5)=1.0D0
		

	!$OMP DO SCHEDULE (STATIC)
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
			
						  CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  
						  CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						   CALL ROTATEB(N,INVTRI,ClefT,Cleft_ROT,ANGLE1,ANGLE2)
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:5)=CLEFT(1:5);RIGHTV(1:5)=CRIGHT(1:5)						  
						  CALL CONS2PRIM2(N)
						  
						  ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 
							  
							  
							  
							  EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:5,1:5)=IMPDIAG(i,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:5,1:5)=IMPOFF(i,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
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
						  
						  IMPDIAG(i,1:5,1:5)=IMPDIAG(i,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:5,1:5)=IMPOFF(i,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
                                                
						  
						  
						  END IF
						  

		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		   IMPDIAG(i,:,:)=0.0
		IMPOFF(i,:,:,:)=0.0
		if (turbulence.eq.1)then
		impdiagt(i,:)=0.0
		IMPOFFt(i,:,:)=0.0
		end if
		    
		    DO L=1,IELEM(N,I)%IFCA
				      B_CODE=0
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
 				  mul1=IELEM(N,I)%SURF(L)
				      B_CODE=0
				      CLEFT(1:5)=U_C(I)%VAL(1,1:5)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:5)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:5)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    
								    
								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
								    cright(1:5)=rightv(1:5)
				  				   
				  				  				  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:5)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:5)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:5)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:5)
							    ELSE
							     
							      CRIGHT(1:5)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:5)
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
							      CRIGHT(1:5)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:5)
							    ELSE
							     
							      CRIGHT(1:5)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:5)
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
				      
				       CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	
				      
				      
				      IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						   CALL ROTATEB(N,INVTRI,ClefT,Cleft_ROT,ANGLE1,ANGLE2)
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:5)=CLEFT(1:5);RIGHTV(1:5)=CRIGHT(1:5)						  
						  CALL CONS2PRIM2(N)
						  
						 ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							    EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    
							    
							  
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:5,1:5)=IMPDIAG(i,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:5,1:5)=IMPOFF(i,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
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
						  
						  IMPDIAG(i,1:5,1:5)=IMPDIAG(i,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:5,1:5)=IMPOFF(i,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF

						
				   
				  
		    END DO
	END DO
	!$OMP END DO

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		!$OMP DO SCHEDULE(GUIDED)	
		do i=1,kmaxe
		    IMPDIAG(i,1,1)=IMPDIAG(i,1,1)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,2,2)=IMPDIAG(i,2,2)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,3,3)=IMPDIAG(i,3,3)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,4,4)=IMPDIAG(i,4,4)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,5,5)=IMPDIAG(i,5,5)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		end do
		!$OMP END DO
	  ELSE
	!$OMP DO SCHEDULE(GUIDED)
! 	
	  do i=1,kmaxe

	    IMPDIAG(I,1,1)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,1,1))
	      IMPDIAG(I,2,2)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,2,2))
	      IMPDIAG(I,3,3)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,3,3))
	      IMPDIAG(I,4,4)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,4,4))
	    IMPDIAG(I,5,5)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAG(i,5,5))
	end do
	!$OMP END DO
      END IF






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
 if (turbulence.eq.1)CALL SOURCES_derivatives_COMPUTATION(N)
if (rungekutta.eq.10)then
!$OMP DO SCHEDULE(GUIDED)
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)-sht(i,nvar)
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
!$OMP DO SCHEDULE(GUIDED)
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
    
    !IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt) +(IMPDIAGt(i,nvar)*ielem(n,i)%dtl/ielem(n,I)%totvolume)-(sht(i,nvar)/ielem(n,I)%totvolume))/ielem(n,i)%dtl
!     IMPDIAGt(i,nvar)=(ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt) +(IMPDIAGt(i,nvar)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl)-(sht(i,nvar))
    IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(i,1))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(i,nvar)=ielem(n,I)%totvolume*((1.0D0/ielem(n,I)%dtl)+(1.5D0/DT))+(IMPDIAGt(i,1))
!     IMPDIAGt(i,nvar)=(ielem(n,I)%totvolume*(((dt+1.5d0*ielem(n,I)%dtl)/dt)+(IMPDIAGt(i,nvar)*ielem(n,i)%dtl/ielem(n,I)%totvolume))/ielem(n,i)%dtl)
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
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar,KAS
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1
	REAL,DIMENSION(4,4)::IDENTITY1
	real,dimension(4,4)::convj,diffj
	KMAXE=XMPIELRANK(N)
	IDENTITY1(:,:)=ZERO
	IDENTITY1(1,1)=1.0D0
	IDENTITY1(2,2)=1.0D0
	IDENTITY1(3,3)=1.0D0
	IDENTITY1(4,4)=1.0D0
	
		

	!$OMP DO SCHEDULE (STATIC)
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
			
						  CALL ROTATEF2D(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						  CALL ROTATEb2D(N,invTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,invTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:4)=CLEFT(1:4);RIGHTV(1:4)=CRIGHT(1:4)						  
						  CALL CONS2PRIM2D2(N)
						  
						 ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV)
						  
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
							  Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					 
						      
						     VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:4,1:4)=IMPDIAG(i,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2D(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:4,1:4)=IMPOFF(i,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
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
						  
						  IMPDIAG(i,1:4,1:4)=IMPDIAG(i,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:4,1:4)=IMPOFF(i,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO SCHEDULE (STATIC) 
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
				      CLEFT(1:4)=U_C(I)%VAL(1,1:4)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:4)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:4)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  KAS=1
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2D(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    
								    
								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:4)=rightv(1:4)
				  				    
				  				  	KAS=2			  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:4)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:4)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      KAS=3
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:4)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:4)
							    ELSE
							     
							      CRIGHT(1:4)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:4)
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
							      CRIGHT(1:4)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:4)
							    ELSE
							     
							      CRIGHT(1:4)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:4)
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
				      
				    
						  
						  
						 CALL ROTATEF2D(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						  CALL ROTATEb2D(N,invTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,invTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	
						  
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:4)=CLEFT(1:4);RIGHTV(1:4)=CRIGHT(1:4)						  
						  CALL CONS2PRIM2D2(N)
						  
						ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV)
						  
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
							  Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
! 						      viscots=viscots/((0.5*(cleft(1)+cRIGHT(1)))*ielem(n,i)%dih(l))
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(i,1:4,1:4)=IMPDIAG(i,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2D(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:4,1:4)=IMPOFF(i,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
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
						  
						  IMPDIAG(i,1:4,1:4)=IMPDIAG(i,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(i,l,1:4,1:4)=IMPOFF(i,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
			
						
				   
				  
		    END DO
	END DO
	!$OMP END DO

	
	
	
	
	
	IF (RUNGEKUTTA.EQ.10)THEN
		!$OMP DO SCHEDULE(GUIDED)	
		do i=1,kmaxe
				  
		    IMPDIAG(i,1,1)=IMPDIAG(i,1,1)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,2,2)=IMPDIAG(i,2,2)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,3,3)=IMPDIAG(i,3,3)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    IMPDIAG(i,4,4)=IMPDIAG(i,4,4)+(ielem(n,I)%totvolume/ielem(n,I)%dtl)
		    
		end do
		!$OMP END DO
	  ELSE
	!$OMP DO SCHEDULE(GUIDED)	  
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
!$OMP DO SCHEDULE(GUIDED)
do i=1,kmaxe
    if (turbulence.eq.1)then
    do nvar=1,turbulenceequations
!    
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+(ielem(n,I)%totvolume/(ielem(n,I)%dtl*(5/CFL)))-sht(i,nvar)
    end do
    end if
    if (passivescalar.gt.0)then
    do nvar=turbulenceequations+1,turbulenceequations+passivescalar
    IMPDIAGt(i,nvar)=IMPDIAGt(i,nvar)+(ielem(n,I)%totvolume/(ielem(n,I)%dtl*(5/CFL)))
    end do
    end if
end do
!$OMP END DO
else
!$OMP DO SCHEDULE(GUIDED)
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
	
	
SUBROUTINE CALCULATE_JACOBIANLM(N)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 3D with low-memory footprint
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1
	REAL,DIMENSION(5,5)::IDENTITY1
	real,dimension(5,5)::convj,diffj
	KMAXE=XMPIELRANK(N)
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
			
						  CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  CALL ROTATEF(N,invTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,invTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:5)=CLEFT(1:5);RIGHTV(1:5)=CRIGHT(1:5)						  
						  CALL CONS2PRIM2(N)
						  
						ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  
							 EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 
							  
							  
							  
							  EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:5,1:5)=IMPDIAG(1,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:5,1:5)=IMPOFF(1,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
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
						  
						  IMPDIAG(1,1:5,1:5)=IMPDIAG(1,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:5,1:5)=IMPOFF(1,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
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
				      CLEFT(1:5)=U_C(I)%VAL(1,1:5)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:5)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:5)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner(N,Iconsidered,facex)
								    CORDS(1:3)=zero
								    CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    poz(1)=cords(3)
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    
								    
								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED)
								    cright(1:5)=rightv(1:5)
				  				   
				  				  				  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:5)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:5)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:5)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:5)
							    ELSE
							     
							      CRIGHT(1:5)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:5)
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
							      CRIGHT(1:5)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:5)
							    ELSE
							     
							      CRIGHT(1:5)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:5)
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
				      
				       CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  CALL ROTATEF(N,invTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,invTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	
						  
						  
						  
						  END IF
						  
						  				  
						  
						  LEFTV(1:5)=CLEFT(1:5);RIGHTV(1:5)=CRIGHT(1:5)						  
						  CALL CONS2PRIM2(N)
						  
						  ASOUND1=SQRT(LEFTV(5)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(5)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND(N,LEFTV,RIGHTV)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							 
							  
							  
							  
							  EDDYFL(4:6)= ILOCAL_RECON3(I)%GRADs(1,1:3);EDDYFL(7:9)=ILOCAL_RECON3(I)%GRADs(2,1:3)
							  EDDYFL(10:12)=ILOCAL_RECON3(I)%GRADs(3,1:3);EDDYFL(13:15)=ILOCAL_RECON3(I)%GRADs(5,1:3)
							  EDDYFL(16:18)=ILOCAL_RECON3(I)%GRADs(6,1:3)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:5,1:5)=IMPDIAG(1,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:5,1:5)=IMPOFF(1,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
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
						  
						  IMPDIAG(1,1:5,1:5)=IMPDIAG(1,1:5,1:5)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:5,1:5)=IMPOFF(1,l,1:5,1:5)+(((OO2*CONVJ(1:5,1:5))&
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
 CALL SOURCES_derivatives2D(N,ICONSIDERED)
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
	
	

SUBROUTINE CALCULATE_JACOBIAN_2DLM(N)
 !> @brief
!> This subroutine computes the approximate jacobian for implicit time stepping in 2D with low-memory footprint
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,nvar
	REAL::sum_detect,NORMS,VPP,ASOUND1,ASOUND2,MUL1
	REAL,DIMENSION(4,4)::IDENTITY1
	real,dimension(4,4)::convj,diffj
	KMAXE=XMPIELRANK(N)
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
			
						  CALL ROTATEF2D(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						   CALL ROTATEb2D(N,invTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,invTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and so
						  END IF
						  
						  				  
						  
						  LEFTV(1:4)=CLEFT(1:4);RIGHTV(1:4)=CRIGHT(1:4)						  
						  CALL CONS2PRIM2D2(N)
						  
						ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							   EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:4,1:4)=IMPDIAG(1,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2D(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:4,1:4)=IMPOFF(1,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
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
						  
						  IMPDIAG(1,1:4,1:4)=IMPDIAG(1,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:4,1:4)=IMPOFF(1,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
						  -((OO2*vpp)*IDENTITY1))*MUL1)
						  
						 
						  
						  
						  END IF
		    END DO
	else
	
	ICONSIDERED=I	
				
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
				      CLEFT(1:4)=U_C(I)%VAL(1,1:4)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:4)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:4)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2D(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    
								    
								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:4)=rightv(1:4)
				  				   
				  				  				  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:4)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:4)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								IF (FASTEST.EQ.1)THEN
							      CRIGHT(1:4)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:4)
							    ELSE
							     
							      CRIGHT(1:4)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:4)
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
							      CRIGHT(1:4)=SOLCHANGER(IELEM(N,I)%INEIGHN(l))%SOL(IELEM(N,i)%Q_FACE(l)%Q_MAPL(1),1:4)
							    ELSE
							     
							      CRIGHT(1:4)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:4)
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
				      
				    CALL ROTATEF2D(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2D(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						   CALL ROTATEb2D(N,invTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEb2D(N,invTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and so
						  END IF
						  
						  				  
						  
						  LEFTV(1:4)=CLEFT(1:4);RIGHTV(1:4)=CRIGHT(1:4)						  
						  CALL CONS2PRIM2D2(N)
						  
						 ASOUND1=SQRT(LEFTV(4)*GAMMA/LEFTV(1))+abs(CLEFT_ROT(2)/CLEFT_ROT(1))
						  ASOUND2=SQRT(RIGHTV(4)*GAMMA/RIGHTV(1))+abs(Cright_ROT(2)/Cright_ROT(1))
						  
						  VPP=MAX(ASOUND1,ASOUND2)
						  
						  IF (ITESTCASE.EQ.4)THEN
						  CALL SUTHERLAND2D(N,LEFTV,RIGHTV)
						  
						  viscots=(VISCL(1)+VISCL(2))*OO2
						  mul1=IELEM(N,I)%SURF(L)
						  
						  viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*mul1&
						  /ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))&
						  *viscots*mul1/(ielem(n,i)%dih(l)*prandtl)))
						  VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  
						  
						  
						  
						  IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= ILOCAL_RECON3(I)%GRADs(1,1:2);EDDYFL(6:7)=ILOCAL_RECON3(I)%GRADs(2,1:2)
							  EDDYFL(8:9)=ILOCAL_RECON3(I)%GRADs(4,1:2)
							  EDDYFL(10:11)=ILOCAL_RECON3(I)%GRADs(5,1:2)
							    
							    
							  eddyfr=eddyfl
							    Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					 
						      
						      VISCOTS=OO2*((viscl(1)+viscl(3))+(viscl(2)+viscl(4)))
		
						      viscots=max((4.0/(3.0*0.5*(cleft(1)+cRIGHT(1))))*viscots*&
						      mul1/ielem(n,i)%dih(l),((gamma/(0.5*(cleft(1)+cRIGHT(1))))*&
						      viscots*mul1/(ielem(n,i)%dih(l)*(prandtl+prtu))))
						      
						      VPP=MAX(ASOUND1,ASOUND2)+VISCOTS
						  END IF
						  
						  
						  IMPDIAG(1,1:4,1:4)=IMPDIAG(1,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2D(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:4,1:4)=IMPOFF(1,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
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
						  
						  IMPDIAG(1,1:4,1:4)=IMPDIAG(1,1:4,1:4)+(OO2*((vpp*identity1))*MUL1)
						  CALL COMPUTE_JACOBIANSE2d(N,EIGVL,Cright,GAMMA,ANGLE1,ANGLE2)
						  convj=eigvl
						  IMPOFF(1,l,1:4,1:4)=IMPOFF(1,l,1:4,1:4)+(((OO2*CONVJ(1:4,1:4))&
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
 CALL SOURCES_derivatives2d(N,ICONSIDERED)
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

END module implicit_fluxes
