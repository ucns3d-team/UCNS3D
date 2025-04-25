module implicit_time
USE LIBRARY
USE TRANSFORM
USE FLOW_OPERATIONS
use implicit_fluxes
use COMMUNICATIONS
USE DECLARATION
IMPLICIT NONE

 contains



subroutine RELAXATION(N)
 !> @brief
!> This subroutine solves the linear system for implicit time stepping either through jacobian or LU-SGS in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,K,II,SWEEPS,kmaxe,nvar,igoflux, icaseb,N_NODE,IBFC,SRF
real::impres1,impres2,impres3,TEMPXX
real:: w1,w2,w3,denx
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::LSCQM1
REAL::B1_imp(1:nof_Variables),DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::ICONSIDERED, FACEX, POINTX
INTEGER::NGP
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS

SWEEPS=10
kmaxe=xmpielrank(n)

impdu(:,:)=zero


du1=zero
b1_imp=zero
LSCQM1=zero
DUR=zero; dul=zero
DURR=zero; DULR=zero



call CALCULATE_JACOBIAN(N)
IF (RFRAME.EQ.0) THEN
!$OMP DO
do i=1,kmaxe
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(i,1:nof_Variables,1:nof_Variables)
impdiag(i,1,1)=1.0d0/lscqm1(1,1)
impdiag(i,2,2)=1.0d0/lscqm1(2,2)
impdiag(i,3,3)=1.0d0/lscqm1(3,3)
impdiag(i,4,4)=1.0d0/lscqm1(4,4)
impdiag(i,5,5)=1.0d0/lscqm1(5,5)
end do
!$OMP END DO
end if
if (SRFG.EQ.1)then
!$OMP DO
do i=1,kmaxe
  lscqm1(1:5,1:5)=impdiag(i,1:5,1:5)
    w1=LSCQM1(4,3)
    w2=LSCQM1(2,4)
    w3=LSCQM1(3,2)
    !INVERSE OF THE MATRIX IN CASE OF SOURCE TERM
    DENX=(LSCQM1(2,2)*LSCQM1(3,3)*LSCQM1(4,4)+LSCQM1(2,2)*w1**2+LSCQM1(3,3)*w2**2+LSCQM1(4,4)*w3**2)
    IMPDIAG(I,1,1)=1.0D0/LSCQM1(1,1)
    IMPDIAG(I,2,2)=(LSCQM1(3,3)*LSCQM1(4,4)+w1**2)/DENX            
    IMPDIAG(I,2,3)=(LSCQM1(4,4)*w3+w1*w2)/DENX
    IMPDIAG(I,2,4)=(-LSCQM1(3,3)*w2+w1*w3)/DENX            
    IMPDIAG(I,3,2)=(-LSCQM1(4,4)*w3+w1*w2)/DENX
    IMPDIAG(I,3,3)=(LSCQM1(2,2)*LSCQM1(4,4)+w2**2)/DENX
    IMPDIAG(I,3,4)=(LSCQM1(2,2)*w1+w2*w3)/DENX
    IMPDIAG(I,4,2)=(LSCQM1(3,3)*w2+w1*w3)/DENX            
    IMPDIAG(I,4,3)=(-LSCQM1(2,2)*w1+w2*w3)/DENX
    IMPDIAG(I,4,4)=(LSCQM1(2,2)*LSCQM1(3,3)+w3**2)/DENX
    IMPDIAG(I,5,5)=1.0D0/LSCQM1(5,5)
end do
!$OMP END DO
END IF
IF(MRF.EQ.1)THEN
 SRF=0
!$OMP DO
do i=1,kmaxe
	SRF=ILOCAL_RECON3(I)%MRF
	IF (ILOCAL_RECON3(i)%MRF.EQ.0)THEN
        lscqm1(1:5,1:5)=impdiag(i,1:5,1:5)
        impdiag(i,1,1)=1.0d0/lscqm1(1,1)
        impdiag(i,2,2)=1.0d0/lscqm1(2,2)
        impdiag(i,3,3)=1.0d0/lscqm1(3,3)
        impdiag(i,4,4)=1.0d0/lscqm1(4,4)
        impdiag(i,5,5)=1.0d0/lscqm1(5,5)
    ELSE
        lscqm1(1:5,1:5)=impdiag(i,1:5,1:5)
        w1=LSCQM1(4,3)
        w2=LSCQM1(2,4)
        w3=LSCQM1(3,2)
        !INVERSE OF THE MATRIX IN CASE OF SOURCE TERM
        DENX=(LSCQM1(2,2)*LSCQM1(3,3)*LSCQM1(4,4)+LSCQM1(2,2)*w1**2+LSCQM1(3,3)*w2**2+LSCQM1(4,4)*w3**2)
        IMPDIAG(I,1,1)=1.0D0/LSCQM1(1,1)
        IMPDIAG(I,2,2)=(LSCQM1(3,3)*LSCQM1(4,4)+w1**2)/DENX            
        IMPDIAG(I,2,3)=(LSCQM1(4,4)*w3+w1*w2)/DENX
        IMPDIAG(I,2,4)=(-LSCQM1(3,3)*w2+w1*w3)/DENX            
        IMPDIAG(I,3,2)=(-LSCQM1(4,4)*w3+w1*w2)/DENX
        IMPDIAG(I,3,3)=(LSCQM1(2,2)*LSCQM1(4,4)+w2**2)/DENX
        IMPDIAG(I,3,4)=(LSCQM1(2,2)*w1+w2*w3)/DENX
        IMPDIAG(I,4,2)=(LSCQM1(3,3)*w2+w1*w3)/DENX            
        IMPDIAG(I,4,3)=(-LSCQM1(2,2)*w1+w2*w3)/DENX
        IMPDIAG(I,4,4)=(LSCQM1(2,2)*LSCQM1(3,3)+w3**2)/DENX
        IMPDIAG(I,5,5)=1.0D0/LSCQM1(5,5)
    END IF
end do
!$OMP END DO
END IF

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
!$OMP DO
do i=1,kmaxe
impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(I,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END DO
!$OMP end DO
END IF


IF (RELAX.EQ.1)THEN
DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
    IF(MRF.EQ.1)THEN
        SRF=ILOCAL_RECON3(I)%MRF
    END IF
if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))
end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
end if
end if
if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(i,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
		IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		b1T(:)=b1T(:)-(IMPofft(i,L,:)*DUt1(:))
		end if
END DO	!loop f

			
else

DO L=1,IELEM(N,I)%IFCA	!loop3			
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
               IF (ILOCAL_RECON3(i)%MRF.EQ.1)THEN
                    !RETRIEVE ROTATIONAL VELOCITY IN CASE OF ROTATING REFERENCE FRAME TO CALCULATE THE CORRECT VALUE OF THE BOUNDARY CONDITION
                    SRF_SPEED(2:4)=ILOCAL_RECON3(I)%ROTVEL(L,1,1:3)
                    CALL ROTATEF(N,SRF_SPEEDROT,SRF_SPEED,ANGLE1,ANGLE2)
                END IF
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
								    IF(PER_ROT.EQ.1)THEN
                                        DU1(2:4)=ROTATE_PER_1(DU1(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
								    END IF
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
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
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightv(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if ((ibound(n,ielem(n,i)%ibounds(l))%icode.eq.5).or.(ibound(n,ielem(n,i)%ibounds(l))%icode.eq.50))then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
                                IF(PER_ROT.EQ.1)THEN
                                    DU1(2:4)=ROTATE_PER_1(DU1(2:4),ibound(n,ielem(n,i)%ibounds(l))%icode,angle_per)
                                END IF
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(i,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1T(:)=b1T(:)-(IMPofft(i,L,:)*DUt1(:))
end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=MATMUL(impdiag(i,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)


END DO!sweeps


else


DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DUMMY12T(:)=zero
end if

if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))
! dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))

end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)
! dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
! DUMMY12T(:)=zero
end if
end if



if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%reorient(l).eq.0)then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(i,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
                                 igoflux=0
		     IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                                                icaseb=1        !periodic mine
                                                        else
                                                                icaseb=3        !physical
                                                        end if
                                                   
                                                    ELSE
                                                    
                                                                icaseb=2!no boundaries interior
                                                    
                                                    end if
                                else
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                                                icaseb=4
                                                                end if
                                                    else
                                                                icaseb=5
                                                    
                                                    end if
                                
                                
                                end if
				  IF (icaseb.le.2)THEN
                                        IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
                                            igoflux=1
                                        else
                                            igoflux=0
                                        end if
                                else
                                        igoflux=2
                                end if
                                 if (IELEM(N,I)%REORIENT(l).EQ.0)then
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
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
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(i,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*&
(b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR))


end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)
 


!$OMP DO
do i=1,kmaxe,-1	!loop2
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DUMMY12T(:)=zero
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%reorient(l).eq.1)then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
+(impofft(i,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
end if
		end if
END DO	!loop f
			
else

DO L=1,IELEM(N,I)%IFCA	!loop3	

                          igoflux=0
		     IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                                                icaseb=1        !periodic mine
                                                        else
                                                                icaseb=3        !physical
                                                        end if
                                                   
                                                    ELSE
                                                    
                                                                icaseb=2!no boundaries interior
                                                    
                                                    end if
                                else
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                                                icaseb=4
                                                                end if
                                                    else
                                                                icaseb=5
                                                    
                                                    end if
                                
                                
                                end if
				  IF (icaseb.le.2)THEN
                                        IF (IELEM(N,I)%REORIENT(l).EQ.1)THEN
                                            igoflux=1
                                        else
                                            igoflux=0
                                        end if
                                else
                                        igoflux=2
                                end if
                                 if (IELEM(N,I)%REORIENT(l).EQ.1)then

				
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

 dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
+(impofft(i,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
end if
		end if
END DO	!loop f
end if

IMPDU(I,1:nof_variables)=IMPDU(I,1:nof_variables)-matmul(impdiag(i,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN

impdu(i,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdu(i,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)-&
(impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR))

end if
END DO	!loop elements
!$OMP end DO 
 
 call EXHBOUNDHIGHER2(N)
 

END DO!sweeps

end if







END SUBROUTINE RELAXATION


subroutine RELAXATION_lm(N)
 !> @brief
!> This subroutine solves the linear system for implicit time stepping either through jacobian or LU-SGS in 3D with low memory footprint
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,K,II,SWEEPS,kmaxe,nvar,IBFC,igoflux
real::impres1,impres2,impres3
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::LSCQM1
REAL::B1_imp(1:nof_Variables),DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::ICONSIDERED, FACEX, POINTX
INTEGER::NGP,n_node
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS
REAL,ALLOCATABLE,DIMENSION(:,:)::IMPDIAGT
REAL,ALLOCATABLE,DIMENSION(:,:,:)::IMPDIAG,IMPOFFt
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)::IMPOFF


ALLOCATE (IMPDIAG(1,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPOFF(1,6,1:nof_Variables,1:nof_Variables))
IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(1,6,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(1,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF






SWEEPS=4
kmaxe=xmpielrank(n)

impdu(iconsidered,:)=zero


du1=zero
b1_imp=zero
LSCQM1=zero
DUR=zero; dul=zero
DURR=zero; DULR=zero







IF (RELAX.EQ.1)THEN
DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
iconsidered=i
call CALCULATE_JACOBIANlm(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(1,1:nof_Variables,1:nof_Variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)
impdiag(1,5,5)=1.0d0/lscqm1(5,5)



IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END IF

if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))

end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)

end if
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(1,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
		IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		b1T(:)=b1T(:)-(IMPofft(1,L,:)*DUt1(:))
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3			
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_innerx(N,Iconsidered,facex,VEXT,NODES_LIST)

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
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1,6)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(1,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1T(:)=b1T(:)-(IMPofft(1,L,:)*DUt1(:))
end if
END DO	!loop f
end if

IMPDU(I,1:nof_variables)=MATMUL(impdiag(1,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)


END DO!sweeps


else


DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
iconsidered=i
call CALCULATE_JACOBIANlm(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(1,1:nof_Variables,1:nof_Variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)
impdiag(1,5,5)=1.0d0/lscqm1(5,5)



IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END IF




if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))
DUMMY12T(:)=zero
end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
DUMMY12T(:)=zero
end if
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%reorient(l).eq.0)then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(1,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
				if (ielem(n,i)%reorient(l).eq.0)then
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
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
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								 select case(b_code)
								    case(1)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(1,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f
end if

IMPDU(I,1:nof_variables)=matmul(impdiag(1,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*&
(b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR))


end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)
 


!$OMP DO
do i=1,kmaxe,-1	!loop2

iconsidered=i
call CALCULATE_JACOBIANlm(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(1,1:nof_Variables,1:nof_Variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)
impdiag(1,5,5)=1.0d0/lscqm1(5,5)



IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END IF







dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DUMMY12T(:)=zero
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%reorient(l).eq.1)then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
+(impofft(1,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
end if
		end if
END DO	!loop f

				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
				if (ielem(n,i)%reorient(l).eq.1)then
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

 dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
+(impofft(1,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
end if
		end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=IMPDU(I,1:nof_variables)-matmul(impdiag(1,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN

impdu(1,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(I,6:5+TURBULENCEEQUATIONS+PASSIVESCALAR)&
-((impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)))

end if
END DO	!loop elements
!$OMP end DO 
 
 call EXHBOUNDHIGHER2(N)
 

END DO!sweeps

end if



DEALLOCATE(IMPDIAG,IMPOFF,IMPOFFT,IMPDIAGT)



END SUBROUTINE RELAXATION_lm


subroutine RELAXATION2d(N)
 !> @brief
!> This subroutine solves the linear system for implicit time stepping either through jacobian or LU-SGS in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,K,II,SWEEPS,kmaxe,nvar,igoflux, icaseb,N_NODE,IBFC
real::impres1,impres2,impres3
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::LSCQM1
REAL::B1_imp(1:nof_Variables),DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::ICONSIDERED, FACEX, POINTX
INTEGER::NGP
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,srf_speed
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS

SWEEPS=10
kmaxe=xmpielrank(n)

impdu(:,:)=zero


du1=zero
b1_imp=zero
LSCQM1=zero
DUR=zero; dul=zero
DURR=zero; DULR=zero



call CALCULATE_JACOBIAN_2d(N)

!$OMP DO
do i=1,kmaxe
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(i,1:nof_Variables,1:nof_Variables)
impdiag(i,1,1)=1.0d0/lscqm1(1,1)
impdiag(i,2,2)=1.0d0/lscqm1(2,2)
impdiag(i,3,3)=1.0d0/lscqm1(3,3)
impdiag(i,4,4)=1.0d0/lscqm1(4,4)
end do
!$OMP END DO


IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
!$OMP DO
do i=1,kmaxe
impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(I,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END DO
!$OMP end DO
END IF


IF (RELAX.EQ.1)THEN
DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2

if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))
end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
end if
end if
if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(i,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
		IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		b1T(:)=b1T(:)-(IMPofft(i,L,:)*DUt1(:))
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3			
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=angle1
				NY=angle2
				
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    



								    
								     select case(b_code)
								    case(1)
								    du1(:)=zero
								    
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								   
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if


								

								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	

B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(i,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1T(:)=b1T(:)-(IMPofft(i,L,:)*DUt1(:))
end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=MATMUL(impdiag(i,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
end if
	!loop elements

end do
!$OMP end DO

 call EXHBOUNDHIGHER2(N)

 

END DO!sweeps


else


DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DUMMY12T(:)=zero
end if

if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))
! DUMMY12T(:)=zero
end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)
! dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
! DUMMY12T(:)=zero
end if
end if



if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%INEIGH(L).lt.i) then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(i,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
				  igoflux=0
		     IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                                                icaseb=1        !periodic mine
                                                        else
                                                                icaseb=3        !physical
                                                        end if
                                                   
                                                    ELSE
                                                    
                                                                icaseb=2!no boundaries interior
                                                    
                                                    end if
                                else
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                                                icaseb=4
                                                                end if
                                                    else
                                                                icaseb=5
                                                    
                                                    end if
                                
                                
                                end if
				  IF (icaseb.le.2)THEN
                                        if (ielem(n,i)%INEIGH(L).lt.i)then
                                            igoflux=1
                                        else
                                            igoflux=0
                                        end if
                                else
                                        igoflux=2
                                end if
                                 if (ielem(n,i)%INEIGH(L).lt.i)then
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=ANGLE1
				NY=ANGLE2
				
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    

								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6,9,99)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(i,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*&
(b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR))


end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)
 


!$OMP DO
do i=1,kmaxe	!loop2
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DUMMY12T(:)=zero
end if

if (ielem(n,i)%interior.eq.0)then
                        DO L=1,IELEM(N,I)%IFCA	!loop3
                                    if (ielem(n,i)%INEIGH(L).gt.i)then
                                                    DU1(1:nof_Variables)=zero
                                                    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                    DUT1(:)=zero
                                                    end if
                                                        
                                                        
                                
                                        du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
                                                                                            
                                                    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
                                                    DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)   
                                                    END IF	
                                        
                                                                    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1))
                                                IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
                                                +(impofft(i,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
                                                end if
                                        end if
                        END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
					  igoflux=0
                                    
                                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                                                                icaseb=1        !periodic mine
                                                                        else
                                                                                icaseb=3        !physical
                                                                        end if
                                                                
                                                                    ELSE
                                                                    
                                                                                icaseb=2!no boundaries interior
                                                                    
                                                                    end if
                                                else
                                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                                                if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                                                                icaseb=4
                                                                                end if
                                                                    else
                                                                                icaseb=5
                                                                    
                                                                    end if
                                                
                                                
                                                end if
                                                IF (icaseb.le.2)THEN
                                                        IF (IELEM(N,I)%REORIENT(l).EQ.1)THEN
                                                            igoflux=1
                                                        else
                                                            igoflux=0
                                                        end if
                                                else
                                                        igoflux=2
                                                end if
                                                                    if (ielem(n,i)%INEIGH(L).gt.i)then
                                                                    ANGLE1=IELEM(N,I)%FACEANGLEX(L)
                                                                    ANGLE2=IELEM(N,I)%FACEANGLEY(L)
                                                                    NX=angle1
                                                                    NY=angle2
                                                                    
                                            
                                                                                        IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                                                                                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                                                                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                                                                                                            DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
                                                                                                                                IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                                                                                            
                                                                                                                            DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
                                                                                                                                
                                                                                                                                end if
                                                                                                                                                                
                                                                                                                        
                                                                                                                        
                                                                                                                        ELSE
                                                                                                                        
                                                                                                                            
                                                                                                                            
                                                                                                                                                                                        
                                                                                                                            
                                                                                                                        END IF
                                                                                                                ELSE
                                                                                                                    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
                                                                                                                                IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                                                                                            
                                                                                                                            DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
                                                                                                                                
                                                                                                                                end if
                                                                                                                    
                                                                                                                    
                                                                                                                    
                                                                                                                    
                                                                                                                END IF
                                                                                            ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                                                                                            
                                                                                            
                                                                                                
                                                                                                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                                                                                if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                                                                                                                
                                                                                                                DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
                                                                    
                                                                                                                IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                                                                                
                                                                                                                DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
                                                                                                                
                                                                                                                end if
                                                                                        
                                                                                                                END IF
                                                                                                        ELSE 			
                                                                                                            DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
                                                                    
                                                                                                                IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                                                                                
                                                                                                                DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
                                                                                                                
                                                                                                                end if
                                                                                                        
                                                                                                                
                                                ! 								   
                                                                                                                END IF
                                                                                                    END IF
                                                    
                                                    

                                            dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(i,l,1:nof_variables,1:nof_variables),du1))
                                                                                            IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                                                                                            dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
                                                                                            +(impofft(i,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
                                                                                            end if
                                                            end if
END DO	!loop f
end if

! IMPDU(I,1:nof_variables)=matmul(impdiag(i,1:nof_variables,1:nof_variables),IMPDU(I,1:nof_variables)-dummy12(1:nof_variables))


IMPDU(I,1:nof_variables)=IMPDU(I,1:nof_variables)-matmul(impdiag(i,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))


IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN

impdu(i,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdu(i,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)-&
(impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR))

end if
END DO	!loop elements
!$OMP end DO 
 
 call EXHBOUNDHIGHER2(N)
 

END DO!sweeps

end if







END SUBROUTINE RELAXATION2d


subroutine RELAXATION_lm2d(N)
 !> @brief
!> This subroutine solves the linear system for implicit time stepping either through jacobian or LU-SGS in 2D with low memory footprint
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,L,K,II,SWEEPS,kmaxe,nvar,IBFC,igoflux
real::impres1,impres2,impres3
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::LSCQM1
REAL::B1_imp(1:nof_Variables),DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::ICONSIDERED, FACEX, POINTX,N_NODE
INTEGER::NGP
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,srf_Speed
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS
REAL,ALLOCATABLE,DIMENSION(:,:)::IMPDIAGT
REAL,ALLOCATABLE,DIMENSION(:,:,:)::IMPDIAG,IMPOFFt
REAL,ALLOCATABLE,DIMENSION(:,:,:,:)::IMPOFF


ALLOCATE (IMPDIAG(1,1:nof_Variables,1:nof_Variables))
ALLOCATE (IMPOFF(1,4,1:nof_Variables,1:nof_Variables))
IF ((ITESTCASE.EQ.4).AND.((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0)))THEN
ALLOCATE(IMPOFFt(1,4,TURBULENCEEQUATIONS+PASSIVESCALAR))
ALLOCATE(IMPDIAGT(1,TURBULENCEEQUATIONS+PASSIVESCALAR))
END IF

SWEEPS=4
kmaxe=xmpielrank(n)

impdu(iconsidered,:)=zero


du1=zero
b1_imp=zero
LSCQM1=zero
DUR=zero; dul=zero
DURR=zero; DULR=zero







IF (RELAX.EQ.1)THEN
DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
iconsidered=i
call CALCULATE_JACOBIAN_2dlm(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(1,1:nof_Variables,1:nof_Variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)



IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END IF

if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))

end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)

end if
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(1,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
		IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		b1T(:)=b1T(:)-(IMPofft(1,L,:)*DUt1(:))
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3			
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=angle1
				NY=angle2
				
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2dx(N,Iconsidered,facex,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    
								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    
								     select case(b_code)
								    case(1)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

B1_imp(1:nof_variables)=B1_imp(1:nof_variables)-MATMUL(IMPoff(1,L,1:nof_variables,1:nof_variables),DU1(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1T(:)=b1T(:)-(IMPofft(1,L,:)*DUt1(:))
end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=MATMUL(impdiag(1,1:nof_variables,1:nof_variables),b1_imp(1:nof_variables))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)


END DO!sweeps


else


DO II=1,SWEEPS	!loop1
!$OMP DO
do i=1,kmaxe	!loop2
iconsidered=i
call CALCULATE_JACOBIAN_2dlm(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(1,1:nof_Variables,1:nof_Variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)




IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END IF




if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DO NVAR=1,TURBULENCEEQUATIONS+PASSIVESCALAR
B1T(NVAR)=-(RHST(I)%VAL(NVAR)+((((1.5d0*U_CT(I)%VAL(1,NVAR))-(2.0d0*U_CT(I)%VAL(2,NVAR))+(0.5d0*U_CT(I)%VAL(3,NVAR)))/(dt))*IELEM(N,I)%TOTVOLUME))
DUMMY12T(:)=zero
end do
end if
else
B1_imp(1:nof_variables)=-RHS(I)%VAL(1:nof_variables)
dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHST(i)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
DUMMY12T(:)=zero
end if
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%reorient(l).eq.0)then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(1,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
				if (ielem(n,i)%reorient(l).eq.0)then
				ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=angle1
				NY=angle2
				
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)
								  N_NODE=2
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								   
								    
								    LEFTV(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								   

								    
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    
								    DU1(1:nof_variables)=rightV(1:nof_variables)
				  				    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
								    end if
								    
								    select case(b_code)
								    case(1)
								    du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    case(6)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1(1:nof_variables)))
		  IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
		  dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
		  (IMPofft(1,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
		  end if
		end if
END DO	!loop f
end if

IMPDU(I,1:nof_variables)=matmul(impdiag(1,1:nof_variables,1:nof_variables),(b1_imp(1:nof_variables)-dummy12(1:nof_variables)))

IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*&
(b1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR))


end if
END DO	!loop elements
!$OMP end DO

 call EXHBOUNDHIGHER2(N)
 


!$OMP DO
do i=1,kmaxe,-1	!loop2

iconsidered=i
call CALCULATE_JACOBIAN_2dlm(N,ICONSIDERED,impdiag,IMPDIAGT,IMPOFF,IMPOFFT)
  lscqm1(1:nof_Variables,1:nof_Variables)=impdiag(1,1:nof_Variables,1:nof_Variables)
impdiag(1,1,1)=1.0d0/lscqm1(1,1)
impdiag(1,2,2)=1.0d0/lscqm1(2,2)
impdiag(1,3,3)=1.0d0/lscqm1(3,3)
impdiag(1,4,4)=1.0d0/lscqm1(4,4)




IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)=1.0d0/IMPDIAGT(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
END IF







dummy12(:)=zero
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
DUMMY12T(:)=zero
end if

if (ielem(n,i)%interior.eq.0)then
DO L=1,IELEM(N,I)%IFCA	!loop3
	      if (ielem(n,i)%reorient(l).eq.1)then
			    DU1(1:nof_Variables)=zero
			    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			    DUT1(:)=zero
			    end if
				
				
	 
		du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
				     				      
		IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
		DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)   
		END IF	
		
		    dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
+(impofft(1,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
end if
		end if
END DO	!loop f


				
else

DO L=1,IELEM(N,I)%IFCA	!loop3	
				if (ielem(n,i)%reorient(l).eq.1)then
				
	
					IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
						IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								    DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
								  					  
								  
								  
								  ELSE
								  
								    
								    
				  				  				  				  
								    
								  END IF
							ELSE
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
									
									end if
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
								
								DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 			
							      DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
								  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
							
								  
! 								   
							END IF
					    END IF
	
	

 dummy12(1:nof_variables)=dummy12(1:nof_variables)+(matmul(IMPoff(1,l,1:nof_variables,1:nof_variables),du1))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)&
+(impofft(1,l,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dut1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))
end if
		end if
END DO	!loop f
end if
IMPDU(I,1:nof_variables)=IMPDU(I,1:nof_variables)-matmul(impdiag(1,1:nof_variables,1:nof_variables),dummy12(1:nof_variables))
IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN

impdu(1,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(I,5:4+TURBULENCEEQUATIONS+PASSIVESCALAR)&
-((impdiagt(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*dummy12t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)))

end if
END DO	!loop elements
!$OMP end DO 
 
 call EXHBOUNDHIGHER2(N)
 

END DO!sweeps

end if



DEALLOCATE(IMPDIAG,IMPOFF,IMPOFFT,IMPDIAGT)



END SUBROUTINE RELAXATION_lm2d






SUBROUTINE RELAXATION_LUMFREE(N)
IMPLICIT NONE
!> @brief
!> This subroutine solves the linear system for implicit time stepping either through MATRIX FREE LU-SGS low memory footprint
INTEGER,INTENT(IN)::N
INTEGER::I,L,K,II,SWEEPS,kmaxe,nvar,igoflux,icaseb,inds,MODEU
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES)::LSCQM1
REAL::B1_imp(1:nof_Variables),DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::ICONSIDERED, FACEX, POINTX
INTEGER::NGP
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS

SWEEPS=1
kmaxe=xmpielrank(n)
impdu(:,:)=zero

if (TURBULENCE.gt.0)then
dut1=zero
end if


if (dimensiona.eq.3)then
inds=5
CALL CALCULATE_JACOBIAN_3D_MF(n)

else
inds=4
CALL CALCULATE_JACOBIAN_2D_MF(n)
end if




DO II=1,SWEEPS	!loop1
du1=zero
call EXHBOUNDHIGHER2(N)

!$OMP DO
DO I=1,kmaxe   !FORWARD SWEEP !loop 1
iconsidered=i
                    b1_imp=zero
                    if (TURBULENCE.gt.0)b1T=zero
                    MODEU=1
                     if (ielem(n,i)%interior.eq.0)then  
                    CALL LUSGS_INTERIOR(N,ICONSIDERED,MODEU,B1_IMP,B1T)
                     ELSE
                    
                    CALL LUSGS_OUT(N,ICONSIDERED,MODEU,B1_IMP,B1T)
                     END IF
                     
                    if (iscoun.ne.1)then
B1_imp(1:nof_variables)=-(RHS(I)%VAL(1:nof_variables)+((((1.5*U_C(I)%VAL(1,1:nof_Variables))-(2.0d0*U_C(I)%VAL(2,1:nof_Variables))+(0.5d0*U_C(I)%VAL(3,1:nof_Variables)))/(dt))*IELEM(N,I)%TOTVOLUME))-B1_IMP(1:NOF_VARIABLES)
                    else
			
                    B1_IMP(1:NOF_VARIABLES)=(-RHS(I)%VAL(1:nof_variables))-B1_IMP(1:NOF_VARIABLES)
                    end if
                    
                    
                    IMPDU(I,1:nof_variables)=(B1_IMP(1:NOF_VARIABLES)/IMPDIAG_MF(i))
                    
                    
                    
                    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                    
                    if (iscoun.ne.1)then
                    B1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-(RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+((((1.5d0*U_CT(I)%VAL(1,1:TURBULENCEEQUATIONS+PASSIVESCALAR))-(2.0d0*U_CT(I)%VAL(2,1:TURBULENCEEQUATIONS+PASSIVESCALAR))+(0.5d0*U_CT(I)%VAL(3,1:TURBULENCEEQUATIONS+PASSIVESCALAR)))/(dt))*IELEM(N,I)%TOTVOLUME))-B1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                    impdu(i,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)=b1T(:)/impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                    
                    
                    else
                    b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=-RHSt(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-B1t(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                    impdu(i,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)=b1T(:)/impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                    end if
                    
                    
                        end if
                    
                    
                  
END DO
!$OMP END DO 

 call EXHBOUNDHIGHERlu(N)

!$OMP DO
DO I=kmaxe,1,-1     !FORWARD SWEEP !loop 1
iconsidered=i
                     b1_imp=zero
                    MODEU=0
                    if (TURBULENCE.gt.0)b1T=zero
                     if (ielem(n,i)%interior.eq.0)then  !if 1
                        CALL LUSGS_INTERIOR(N,ICONSIDERED,MODEU,B1_IMP,B1T)
                     ELSE
                    
                    CALL LUSGS_OUT(N,ICONSIDERED,MODEU,B1_IMP,B1T)
                     END IF

                   
                    
!                     IMPDU(I,1:nof_variables)=IMPDU(I,1:nof_variables)-(B1_IMP(1:NOF_VARIABLES)/IMPDIAG_MF(i))
                    
                    
                     IMPDU(I,1:nof_variables)=IMPDU(I,1:nof_variables)-(B1_IMP(1:NOF_VARIABLES)/IMPDIAG_MF(i))
                    
                    
                     IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                    impdu(i,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)=impdu(i,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)-b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)/impdiagt(i,1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                        end if
                    
                    
                 
END DO
!$OMP END DO




END DO!sweeps
                                       
END SUBROUTINE RELAXATION_LUMFREE                  
                    
                    
SUBROUTINE LUSGS_INTERIOR(N,ICONSIDERED,MODEU,B1_IMP,B1T)
IMPLICIT NONE
!> @brief
!> This subroutine solves the linear system for implicit time stepping either through MATRIX FREE LU-SGS low memory footprint
INTEGER,INTENT(IN)::N,ICONSIDERED,MODEU
INTEGER::I,L,K,II,inds
real::impres1,impres2,impres3
REAL,INTENT(INOUT)::B1_imp(1:nof_Variables),B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::FACEX, POINTX,igoflux
INTEGER::NGP
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ,VELNORMAL
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS
REAl::DELTAF(1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR)




I=ICONSIDERED

!     MODEU=1!LOWER
!     MODEU=0!UPPER
if (dimensiona.eq.3)then
inds=5
else
inds=4
end if

DO L=1,IELEM(N,I)%IFCA	!loop2     
                                    if (MODEU.EQ.0)THEN
                                    
							       if (IELEM(N,I)%INEIGH(L).lt.ICONSIDERED) CYCLE
							       
							       else
							       
							        if (IELEM(N,I)%INEIGH(L).gt.ICONSIDERED) CYCLE
							        
							       end if
        
        
        
         du1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
         CRIGHT(1:nof_variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_variables)
              IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
                  DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)   
                  CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
                  CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
               END IF	
        ANGLE1=IELEM(N,I)%FACEANGLEX(L)
        ANGLE2=IELEM(N,I)%FACEANGLEY(L)
        NX=ANGLE1
        NY=ANGLE2
        facex=l;
        CALL DELTAFX(ICONSIDERED,FACEX,DELTAF,VELNORMAL)
        
        
        
        
         b1_imp(1:NOF_VARIABLES)=b1_imp(1:NOF_VARIABLES)+(0.5D0*(DELTAF(1:NOF_VARIABLES)-VELNORMAL*du1(1:nof_variables))*IELEM(N,I)%SURF(L))
        
!         
       
        
        
        IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
        
!         b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+IMPofft(i,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
        
        b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+0.5D0*(DELTAF(inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)-VELNORMAL*duT1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))*IELEM(N,I)%SURF(L)
        END IF
        
        

        
END DO
        
        

        
        
        
END SUBROUTINE LUSGS_INTERIOR

                    
SUBROUTINE LUSGS_OUT(N,ICONSIDERED,MODEU,B1_IMP,B1T)
IMPLICIT NONE
!> @brief
!> This subroutine solves the linear system for implicit time stepping either through MATRIX FREE LU-SGS low memory footprint
INTEGER,INTENT(IN)::N,ICONSIDERED,MODEU
INTEGER::I,L,K,II,inds,interfaces,N_NODE
REAL,INTENT(INOUT)::B1_imp(1:nof_Variables),B1T(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DU1(1:nof_Variables),DU2(1:nof_Variables),DUMMY12(1:nof_Variables),C1_imp(1:nof_Variables)
REAL::DUR(nof_Variables),DUL(nof_Variables)
REAL::DURR(nof_Variables),DULR(nof_Variables)
REAL::DUT1(TURBULENCEEQUATIONS+PASSIVESCALAR)
REAL::DUMMY12T(TURBULENCEEQUATIONS+PASSIVESCALAR)
INTEGER::FACEX, POINTX,igoflux
INTEGER::NGP
INTEGER::B_CODE
REAL::ANGLE1,ANGLE2,NX,NY,NZ,VELNORMAL
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr
real,dimension(1:nof_Variables)::leftv,SRF_SPEEDROT,SRF_SPEED
real,dimension(1:nof_Variables)::RIGHTv
REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
REAL,DIMENSION(1:DIMENSIONA)::CORDS
INTEGER::IBFC
REAl::DELTAF(1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR)
I=ICONSIDERED

if (dimensiona.eq.3)then
inds=NOF_VARIABLES
else
inds=NOF_VARIABLES
end if

DO L=1,IELEM(N,I)%IFCA	!COND1   
        B_CODE=0
        
        
        
        IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY !COND2
        
		IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES !COND3
		if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU !COND4
		
		
                                    if (MODEU.EQ.0)THEN
							       if (IELEM(N,I)%INEIGH(L).LT.ICONSIDERED)CYCLE
							       else
							        if (IELEM(N,I)%INEIGH(L).gt.ICONSIDERED)CYCLE
							       end if
! 		if (ielem(n,i)%IELEM(N,I)%INEIGH(L).LT.ICONSIDERED)CYCLE
		
		
		
		
		 DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
		 
		 CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
		 
		 
		 
		 
		IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN  
        DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
        CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
        end if      
        ELSE                                !COND4
		!NOT PERIODIC ONES IN MY CPU
        facex=l;
        
        
        
        if (dimensiona.eq.2)then
		CALL coordinates_face_inner2dx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)

		N_NODE=2
		CORDS(1:2)=zero
		CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE);Pox(1)=cords(1);Poy(1)=cords(2)
		else
		CALL coordinates_face_innerx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)

		 if (ielem(n,ICONSIDERED)%types_faces(FACEX).eq.5)then
                                            N_NODE=4
                                    else
                                            N_NODE=3
                                    end if

		CORDS(1:3)=zero
		CORDS(1:3)=CORDINATES3(N,NODES_LIST,N_NODE);Pox(1)=cords(1);Poy(1)=cords(2);poz(1)=cords(3);
		
		end if
		
		
        LEFTV(1:nof_variables)=IMPDU(i,1:nof_variables)!U_C(ICONSIDERED)%VAL(1,1:NOF_VARIABLES)						           
        IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN  !COND5
		 cturbl(1:turbulenceequations+passivescalar)=IMPDU(I,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
		end if                                                !END COND5
		B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
		
		if (dimensiona.eq.2)then
			


		CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
		else
		
		CALL BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
		end if
		
		DU1(1:nof_variables)=rightV(1:nof_variables)
		
        IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN !COND6
		 dut1(1:turbulenceequations+passivescalar)=cturbr(1:turbulenceequations+passivescalar)
		end if                                              !END COND6  
                                    select case(b_code)     !COND 7
                                    case(1)
                                    du1(:)=zero
								    
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								   
								    end if
								    
								    case(2)
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    case(4)
        
								    
								    
								    case(6)
								    
								    if (ibfc.eq.-1)then
								   du1(:)=zero
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=zero
								    end if
								    
								    
								    else
								    du1(1:nof_variables)=IMPDU(I,1:nof_variables)
								    IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								    dut1(1:turbulenceequations+passivescalar)=IMPDU(I,inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
								    end if
								    
								    
								    
								    end if
								    
								    
								    end select   !END COND7  
								    
								    
				  				  				  				  
								    
								  END IF            !END COND4
							ELSE     !COND3
                                    
                                    
                                     if (MODEU.EQ.0)THEN
							       if (IELEM(N,I)%INEIGH(L).LT.ICONSIDERED)CYCLE
							       else
							        if (IELEM(N,I)%INEIGH(L).gt.ICONSIDERED)CYCLE
							       end if
                                    
                                    
							       DU1(1:nof_variables)=IMPDU(IELEM(N,I)%INEIGH(L),1:nof_variables)
							       CRIGHT(1:nof_Variables)=U_C(IELEM(N,I)%INEIGH(L))%VAL(1,1:nof_Variables)
							       
							       if (MODEU.EQ.0)THEN
							       if (IELEM(N,I)%INEIGH(L).LT.ICONSIDERED)CYCLE
							       else
							        if (IELEM(N,I)%INEIGH(L).gt.ICONSIDERED)CYCLE
							       end if
							       
							       
							       
									IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								      
								      DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IMPDU(IELEM(N,I)%INEIGH(L),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
									CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									end if
							      
							      
							      
							      
							END IF   !END COND3
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS !COND2
					    
					    
						
	IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
	
			if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                
                                    
                
			CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
			(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables)					
			DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      
			IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
			CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)					  
			DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
								  
								 end if
					
								END IF
							ELSE 	
							
                                    
                                    
				DU1(1:nof_variables)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),1:nof_variables)
	      	      CRIGHT(1:nof_Variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),1:nof_Variables) 
								IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
								  
				DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=IEXBOUNDHIRi(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(1),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)
					CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(l)))%SOL&
							      (ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(l)),inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)			  
								 end if
							
								  
! 								   
							END IF
					    END IF       
               
               
               
               
               
               
               
               
               
        ANGLE1=IELEM(N,I)%FACEANGLEX(L)
        ANGLE2=IELEM(N,I)%FACEANGLEY(L)
        NX=ANGLE1
        NY=ANGLE2
        facex=l;
        CALL DELTAFX(ICONSIDERED,FACEX,DELTAF,VELNORMAL)
        
       
        if (b_code.gt.0)then
        if (b_code.ne.5)then
        deltaf=0.0;VELNORMAL=0.0
        end if
        end if

        

        
        
         b1_imp(1:NOF_VARIABLES)=b1_imp(1:NOF_VARIABLES)+(0.5D0*(DELTAF(1:NOF_VARIABLES)-VELNORMAL*du1(1:nof_variables))*IELEM(N,I)%SURF(L))
        

        
        
        
       
        
        
                            
        IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
!          b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+IMPofft(i,L,1:TURBULENCEEQUATIONS+PASSIVESCALAR)*DUt1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
        b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=b1T(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+0.5D0*(DELTAF(inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)-VELNORMAL*duT1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))*IELEM(N,I)%SURF(L)
        END IF
        
END DO
        

        
END SUBROUTINE LUSGS_OUT  



SUBROUTINE DELTAFX(ICONSIDERED,FACEX,DELTAF,VELNORMAL)
IMPLICIT NONE
integer,intent(in)::ICONSIDERED,facex
real::uul,uur,angle1,angle2
real,dimension(1:nof_variables)::CDX1,cdx2
REAl,INTENT(INOUT)::DELTAF(1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR),VELNORMAL
real,dimension(1:nof_variables+turbulenceequations+PASSIVESCALAR)::cleft,cright,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT
real,dimension(1:turbulenceequations+PASSIVESCALAR)::cturbl,cturbr,dut1
real,dimension(1:nof_Variables)::leftv,du1
real,dimension(1:nof_Variables)::RIGHTv
REAL::MP_PINFL,MP_PINFR,GAMMAL,GAMMAR
integer::inds


                            if (dimensiona.eq.2)then
                            inds=NOF_VARIABLES
                            !2 rotate the solution of the neighbour along the face normal and store in cright_rot
                            CALL ROTATEF2D(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)
                            !3 copy to cleft the solution of the neighbour plus the change du
                            CLEFT(1:NOF_VARIABLES)=CRIGHT(1:NOF_VARIABLES)+DU1(1:NOF_VARIABLES)
                            !4 rotate the cleft along the face normal and store in cleft_rot
                            CALL ROTATEF2D(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
                            !5 Compute the spectral radius
                            
                            VELNORMAL=IMPOFF_MF(ICONSIDERED,FACEX)
                            !copy to temporary variables
                            RIGHTV(1:NOF_VARIABLES)=CRIGHT_ROT(1:NOF_VARIABLES)
                            LEFTV(1:NOF_VARIABLES)=CLEFT_ROT(1:NOF_VARIABLES)
                            
                            !transform leftv and rightv from cons to prim variables
                            call cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
                            
                            
                            uul=leftv(2);uur=rightv(2)
                            !evaluate the fluxes
                           
                            CDX2(1:NOF_VARIABLES)=FLUXEVAL2D(LEFTV)
                            
                            
                            leftv=rightv
                             CDX1(1:NOF_VARIABLES)=FLUXEVAL2D(LEFTV)
                            
                           
                            
                            DELTAf(1:NOF_VARIABLES)=CDX2(1:NOF_VARIABLES)-CDX1(1:NOF_VARIABLES)!-VELNORMAL*DU1(1:NOF_VARIABLES)
                           
                           
                            
                            
                            CLEFT_ROT(1:NOF_VARIABLES)=DELTAF(1:NOF_VARIABLES)
                            CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
                            DELTAF(1:NOF_VARIABLES)=CLEFT(1:NOF_VARIABLES)-VELNORMAL*DU1(1:NOF_VARIABLES)
                            
                            IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                            
                            DELTAF(inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)=(((CTURBr(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+duT1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))*uul)-(CTURBR(1:TURBULENCEEQUATIONS+PASSIVESCALAR)*uur))-&
                             VELNORMAL*DUT1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                            end if
                            
                            else
                            
                            
                            
                            inds=NOF_VARIABLES
                            !2 rotate the solution of the neighbour along the face normal and store in cright_rot
                            CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)
                            !3 copy to cleft the solution of the neighbour plus the change du
                            CLEFT(1:NOF_VARIABLES)=CRIGHT(1:NOF_VARIABLES)+DU1(1:NOF_VARIABLES)
                            !4 rotate the cleft along the face normal and store in cleft_rot
                            CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
                            !5 Compute the spectral radius
                            
                            VELNORMAL=IMPOFF_MF(ICONSIDERED,FACEX)
                            !copy to temporary variables
                            RIGHTV(1:NOF_VARIABLES)=CRIGHT_ROT(1:NOF_VARIABLES)
                            LEFTV(1:NOF_VARIABLES)=CLEFT_ROT(1:NOF_VARIABLES)
                            
                            !transform leftv and rightv from cons to prim variables
                            call CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
                            
                            
                            uul=leftv(2);uur=rightv(2)
                            !evaluate the fluxes
                           
                            CDX2(1:NOF_VARIABLES)=FLUXEVAL3D(LEFTV)
                            
                            
                            leftv=rightv
                             CDX1(1:NOF_VARIABLES)=FLUXEVAL3D(LEFTV)
                            
                           
                            
                            DELTAf(1:NOF_VARIABLES)=CDX2(1:NOF_VARIABLES)-CDX1(1:NOF_VARIABLES)!-VELNORMAL*DU1(1:NOF_VARIABLES)
                           
                           
                            
                            
                            CLEFT_ROT(1:NOF_VARIABLES)=DELTAF(1:NOF_VARIABLES)
                            CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
                            DELTAF(1:NOF_VARIABLES)=CLEFT(1:NOF_VARIABLES)-VELNORMAL*DU1(1:NOF_VARIABLES)
                            
                            IF ((TURBULENCE.GT.0).OR.(PASSIVESCALAR.GT.0))THEN
                            
                            DELTAF(inds+1:inds+TURBULENCEEQUATIONS+PASSIVESCALAR)=(((CTURBr(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+duT1(1:TURBULENCEEQUATIONS+PASSIVESCALAR))*uul)-(CTURBR(1:TURBULENCEEQUATIONS+PASSIVESCALAR)*uur))-&
                             VELNORMAL*DUT1(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
                            end if
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            end if
                            
                            
                            

END SUBROUTINE DELTAFX



















end module implicit_time
