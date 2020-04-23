MODULE FLUXES
USE LIBRARY
USE TRANSFORM
USE LOCAL
USE RIEMANN
USE FLOW_OPERATIONS
IMPLICIT NONE
! !**************************DEVELOPED BY PANAGIOTIS TSOUTSANIS**************************!
! !*****************************FMACS RESEARCH GROUP CRANFIELD **************************!
! !*****************************___CRANFIELD_____UNIVERSITY____**************************!
! !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
! !**************************************************************************************!
! !**************************************************************************************!
! !****************************__________DECLARATIONS_______*****************************!
! !**************************************************************************************!
! !************************_____SHARED___GLOBAL_____VARIABLES____************************!
! !**************************************************************************************!
! !**************************************************************************************!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!VARIABLES IN ALPHABETICAL ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 CONTAINS

SUBROUTINE CALCULATE_FLUXESHI(N)
!> @brief
!> This subroutine computes the fluxes for linear-advection equation
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL::GODFLUX2,sum_detect
	INTEGER::I,L,NGP,KMAXE,IQP
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREQUAD3D(N,IGQRULES)
	
	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
	
	call QUADRATURETRIANG(N,IGQRULES)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)

	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE
				  
		IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
		    RHS(I)%VAL(1)=ZERO
		    
		    DO L=1,IELEM(N,I)%IFCA
				  GODFLUX2=ZERO
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				  
				  NORMALVECT=((COS(ANGLE1)*SIN(ANGLE2))*LAMX)+((SIN(ANGLE1)*SIN(ANGLE2))*LAMY)+((COS(ANGLE2))*LAMZ)
				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
					iqp=QP_TRIANGLE
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  do NGP=1,iqp
				      CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
				      CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)

				      CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)
				      GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))	 								      
				  END DO
				    RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2

		    END DO
		END IF
		IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		    RHS(I)%VAL(1)=ZERO
		    
		    DO L=1,IELEM(N,I)%IFCA
				      
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				  NORMALVECT=((COS(ANGLE1)*SIN(ANGLE2))*LAMX)+((SIN(ANGLE1)*SIN(ANGLE2))*LAMY)+((COS(ANGLE2))*LAMZ)
 				  
 				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
					iqp=QP_TRIANGLE
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  GODFLUX2=ZERO
				  do NGP=1,iqp
				      CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)

								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								  CRIGHT(1:nof_variables)=CLEFT(1:nof_variables)
								    
								  END IF
							ELSE
							      CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)

								END IF
							ELSE 								
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
! 								   
							END IF
					    END IF
				      
				      CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)
				      GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))	
				  END DO
				    RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2
				    
		    END DO
		END IF
! 				
	END DO
	!$OMP END DO 
END SUBROUTINE CALCULATE_FLUXESHI
	
	
	
	
	
	
SUBROUTINE CALCULATE_FLUXESHI2D(N)
!> @brief
!> This subroutine computes the fluxes for linear-advection equation in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL::GODFLUX2,sum_detect
	INTEGER::I,L,NGP,KMAXE,IQP
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_TEMP
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATURELINE(N,IGQRULES)
	
! 	WEIGHTS_Q(1:QP_LINE)=WEQUA2D(1:QP_LINE)
	
	WEIGHTS_TEMP(1:QP_LINE_n)=WEQUA2D(1:QP_LINE_n)
	
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE
				if (initcond.eq.3)then
				  lamx=-ielem(n,i)%yyc+0.5d0
				  lamy=ielem(n,i)%xxc-0.5
				  end if
		IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
		    RHS(I)%VAL(1)=ZERO
		    
		    DO L=1,IELEM(N,I)%IFCA
				  GODFLUX2=ZERO
				  NX=IELEM(N,I)%FACEANGLEX(L)
				  NY=IELEM(N,I)%FACEANGLEY(L)
				  
! 				  
				  
				  NORMALVECT=(NX*LAMX)+(NY*LAMY)
				  
				  IQP=QP_LINE_n
				  do NGP=1,iqp
				      CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
				      CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)

				     
				      CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)

				      GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))	 
				  END DO
				    RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2

		    END DO
		END IF
		IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		    RHS(I)%VAL(1)=ZERO
		    
		    DO L=1,IELEM(N,I)%IFCA
				      
				  NX=IELEM(N,I)%FACEANGLEX(L)
				  NY=IELEM(N,I)%FACEANGLEY(L)
				  NORMALVECT=(NX*LAMX)+(NY*LAMY)
				  IQP=QP_LINE_n
 				  
!  				  
				  GODFLUX2=ZERO
				  do NGP=1,iqp
				      CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
! 								
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								  CRIGHT(1:nof_variables)=CLEFT(1:nof_variables)
! 								 
								    
								  END IF
							ELSE
							      CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
!  							       
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
!    									  
								END IF
							ELSE 								
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
!  									
! 								
							END IF
					    END IF
				      
				      CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)
!  				      
				      GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))	
!  				      
				  END DO
				    RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2
				    
		    END DO
		END IF
! 				
	END DO
	!$OMP END DO 
	END SUBROUTINE CALCULATE_FLUXESHI2D	
	
	

	
	
SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE(N)
!> @brief
!> This subroutine computes the convective fluxes for hyperbolic conservation laws
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,IKAS,igoflux, icaseb
	REAL::sum_detect,NORMS
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREQUAD3D(N,IGQRULES)
	
	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
	
	call QUADRATURETRIANG(N,IGQRULES)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)
	
	if(reduce_comp.eq.1)then
	WEIGHTS_T=1.0d0;WEIGHTS_Q=1.0d0
	end if
	
	
	
        do i=1,xmpielrank(n)
            RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
        end do
        !$OMP BARRIER
	!$OMP DO SCHEDULE (STATIC)
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
                    
		    B_CODE=0
		    
		    
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
!                                      IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=(COS(ANGLE1)*SIN(ANGLE2))
				  NY=(SIN(ANGLE1)*SIN(ANGLE2))
				  NZ=(COS(ANGLE2))
				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
					iqp=QP_TRIANGLE_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  do NGP=1,iqp	!for all the gaussian quadrature points
				      CLEFT(1:5)=ILOCAL_RECON3(I)%ULEFT(1:5,L,NGP)	!left mean flow state
				      CRIGHT(1:5)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:5,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
				      
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF
					  
					  
					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					END IF
			
						  CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  END IF
						  
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				      CALL ROTATEB(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				       CALL ROTATEB(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL ROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX=HLLCFLUX
				     				      
				      CASE(4)			!roe
				      
				      
				      CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL rROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				       
				       
				       
				       CASE(5)			!roe
				      
				      
				      CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				       END SELECT
				       
				       
				       
				       
				      
				      GODFLUX2(1:5)=GODFLUX2(1:5)+(RHLLCFLUX(1:5)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  NORMs=(nx*(U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1)))&
						+(nY*(U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1)))&
						+(nz*(U_C(I)%val(1,4)/U_C(I)%val(1,1)))
					      IF (NORMs.GE.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					      IF (NORMs.LT.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBR(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					  
					  END IF
					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
					  
					 
					  
					  
				      END IF
				      
				      
				  END DO
				    
				    RHS(I)%VAL(1:5)=RHS(I)%VAL(1:5)+GODFLUX2(1:5)
!  				     RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)-godflux2(1:5)
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
!  				    RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
! 				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				    end if
!  				    end if
		    END DO
	END DO
	!$OMP END DO
	
	 !$OMP BARRIER
	!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		   
		    
		    DO L=1,IELEM(N,I)%IFCA
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
! 				  IF (icaseb.le.2)THEN
!                                         IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
!                                             igoflux=1
!                                         else
!                                             igoflux=0
!                                         end if
!                                 else
!                                         igoflux=2
!                                 end if
!                                   if (igoflux.ge.1)then
		    
		    
		    
				      
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
 				  
 				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
					iqp=QP_TRIANGLE_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  GODFLUX2=ZERO
				  
				  
								  
				  
				  do NGP=1,iqp
				      B_CODE=0
				      CLEFT(1:5)=ILOCAL_RECON3(I)%ULEFT(1:5,L,NGP)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:5)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:5,IELEM(N,I)%INEIGHN(L),NGP)
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
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
							      CRIGHT(1:5)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:5,IELEM(N,I)%INEIGHN(L),NGP)
							      
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									   
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF
									  
									  

								END IF
							ELSE 			
							
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  
								  
! 								 
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				      
			
						  CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						
						  
						  
								    
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  END IF
						  
				      
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
								      
								      cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
								      cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
								    END IF
				      
				      
				      
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				      CALL ROTATEB(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				       CALL ROTATEB(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL ROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				      
				       CASE(4)			!roe
				      
				      CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				     
				      
				      
				      
				      
				     IF (B_CODE.LE.0)THEN
				      
				      CALL rROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      ELSE
				      
				      
				      CALL ROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      END IF
				      
				       RHLLCFLUX=HLLCFLUX
				     			
                                    CASE(5)			!roe
				      
				      CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      
				     IF (B_CODE.LE.0)THEN
				      
				      CALL tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      ELSE
				      
				      
				      CALL ROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      END IF
				      
				       RHLLCFLUX=HLLCFLUX
				     
				      
				       END SELECT
				       
				     
				      
				      GODFLUX2(1:5)=GODFLUX2(1:5)+(RHLLCFLUX(1:5)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				       
				      
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  NORMs=(nx*(U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1)))&
						+(nY*(U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1)))&
						+(nz*(U_C(I)%val(1,4)/U_C(I)%val(1,1)))
					      IF (NORMs.GE.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					      IF (NORMs.LT.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBR(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					  
					  END IF
					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
					  
					  
					  
				      END IF
				      
				      
				  END DO
				   
				    RHS(I)%VAL(1:5)=RHS(I)%VAL(1:5)+GODFLUX2(1:5)
! 				    if ((igoflux.eq.1))then
! 				    RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)-GODFLUX2(1:5)
! 				    end if
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				     if ((igoflux.eq.1))then
! 				     RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
! 				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    end if
				    end if
! 				    end if
		    END DO
	END DO
	!$OMP END DO

END SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE
	
	
SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE2d(N)
!> @brief
!> This subroutine computes the convective fluxes for hyperbolic conservation laws in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,IKAS,igoflux, icaseb
	REAL::sum_detect,NORMS
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUReline(N,IGQRULES)
	
	WEIGHTS_TEMP(1:qp_line_n)=WEQUA2D(1:qp_line_n)
	
	if(reduce_comp.eq.1)then
	WEIGHTS_TEMP=1.0d0
	end if
	
	do i=1,kmaxe
	RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
	end do
	
	!$OMP BARRIER
	!$OMP DO SCHEDULE (STATIC)
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
! 		    RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
!                     
                     
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
		    if (IELEM(N,I)%REORIENT(l).eq.IELEM(N,IELEM(N,I)%INEIGH(L))%REORIENT(IELEM(N,I)%INEIGHN(L)))then
		    
		     
		     end if
!                               IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
                                
                                
				  GODFLUX2=ZERO
 				  nx=IELEM(N,I)%FACEANGLEX(L)
 				  NY=IELEM(N,I)%FACEANGLEY(L)
 				  angle1=nx
 				  angle2=ny
				 b_code=0
					iqp=qp_line_n
				
				  do NGP=1,iqp	!for all the gaussian quadrature points
				      CLEFT(1:4)=ILOCAL_RECON3(I)%ULEFT(1:4,L,NGP)	!left mean flow state
				      CRIGHT(1:4)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:4,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
! 				      
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF
					  
					  
					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					END IF
			
						  CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						  END IF
						  
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				      CALL ROTATEB2d(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				       CALL ROTATEB2d(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB2d(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL ROE_RIEMANN_SOLVER2d(N,Cleft,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX(1:4)=HLLCFLUX(1:4)
                                        
                                        
                                         CASE(4)			!roe
				      
				      
				      CALL ROTATEB2d(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL rROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX=HLLCFLUX
				     
				      
				       END SELECT
				       
				       
				       
				       
				       
				      
				      GODFLUX2(1:4)=GODFLUX2(1:4)+(RHLLCFLUX(1:4)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  NORMs=(nx*(U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1)))&
						+(nY*(U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1)))
					      IF (NORMs.GE.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					      IF (NORMs.LT.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBR(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					  
					  END IF
					  
 					
					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
					  
! 					 
					  
				      END IF
				      
				      
				  END DO
				    RHS(I)%VAL(1:4)=RHS(I)%VAL(1:4)+GODFLUX2(1:4)
! 				     RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)-GODFLUX2(1:4)
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
!  				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    
				    END IF
! 				    end if
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
! 		    RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
		    
		    DO L=1,IELEM(N,I)%IFCA
				  igoflux=0
				  IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                                                    IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                                        if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                                                icaseb=1        !periodic mine
!                                                                     
                                                        else
                                                                icaseb=3        !physical
!                                                                        
                                                        end if
                                                   
                                                    ELSE
!                                                                  
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

				    b_code=0  
				 nx=IELEM(N,I)%FACEANGLEX(L)
 				  NY=IELEM(N,I)%FACEANGLEY(L)
 				  angle1=nx
 				  angle2=ny
				 
					iqp=qp_line_n
				  GODFLUX2=ZERO
				  do NGP=1,iqp
				      CLEFT(1:4)=ILOCAL_RECON3(I)%ULEFT(1:4,L,NGP)
				      
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:4)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:4,IELEM(N,I)%INEIGHN(L),NGP)
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF
								  
								  
								  IKAS=1
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2d(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
! 								    
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:4)=rightv(1:4)
! 				  				   
				  				  	 IKAS=2			  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:4)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:4,IELEM(N,I)%INEIGHN(L),NGP)
 							     
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF
							      
							      
							       IKAS=3
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
 									   
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
									
									
									
									
								    END IF
									  
									  

								END IF
							ELSE 			
							
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
 								  
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
									
									
								    END IF
								   IKAS=4
! 								   
							END IF
					    END IF
				      
				      
			
						  CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						   
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						  
						  
						  END IF
						  
				      
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
								      
								      cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
								      cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
								    END IF
				      
				      
				      
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				      CALL ROTATEB2d(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				       CALL ROTATEB2d(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB2d(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)   
				      CALL ROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      
				      RHLLCFLUX=HLLCFLUX
				     				      
				     
				        CASE(4)			!roe
				      
				      
				       
				      
				      IF ((B_CODE.le.0))THEN
				      
				      CALL ROTATEB2d(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)  
				      
				      
				      
				      CALL rROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,ROTVL,ROTVR,GAMMA,sl,sr,sm)
				      RHLLCFLUX=HLLCFLUX
				      ELSE
				      
! 				      CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
!                                       CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
				     
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,ROTVL,ROTVR,GAMMA,SL,SR,SM)
				      CALL ROTATEB2d(N,INVTRI,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				      
				      
				      END IF
				      
				      
				       END SELECT
				       
				     
				      
				      GODFLUX2(1:4)=GODFLUX2(1:4)+(RHLLCFLUX(1:4)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				       
				      
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  NORMs=(nx*(U_C(I)%VAL(1,2)/U_C(I)%VAL(1,1)))&
						+(nY*(U_C(I)%VAL(1,3)/U_C(I)%VAL(1,1)))
						
					      IF (NORMs.GE.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					      IF (NORMs.LT.ZERO)THEN
						rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=(NORMs)*CTURBR(1:TURBULENCEEQUATIONS+PASSIVESCALAR)
					      END IF
					  
					  END IF
					  
					  if ((b_code.eq.4).or.(b_code.eq.3))then
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=zero
					  end if
					  
					  
					 
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
! 					  
					  
					  
				      END IF
				      
				      
				  END DO
				   
				    RHS(I)%VAL(1:4)=RHS(I)%VAL(1:4)+GODFLUX2(1:4)
! 				    if ((igoflux.eq.1))then
! 				    RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)-GODFLUX2(1:4)
! 				    end if
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				     if ((igoflux.eq.1))then
! 				     RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
! 				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    end if
! 				    end if
				    end if
		    END DO
	END DO
	!$OMP END DO

END SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE2d

	
	

SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE(N)
!> @brief
!> This subroutine computes the diffusive fluxes for Euler-Navier-Stokes Equations
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,NVAR,KC,IEX,ITTT,IKAS,igoflux, icaseb
	REAL::sum_detect,NORMS
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP
	REAL,DIMENSION(4,3)::LCVGRAD,RCVGRAD
	REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR,3)::LCVGRAD_T,RCVGRAD_T
	real,dimension(5)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(3,3)::taul,taur,TAU
	REAL,DIMENSION(3)::Q,NNN,nall
	REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,RHO12,U12,V12,W12 ,damp,vdamp 
	
	
	
	
	
	
	
	
	
	
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREQUAD3D(N,IGQRULES)
	
	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
	
	call QUADRATURETRIANG(N,IGQRULES)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)

	if(reduce_comp.eq.1)then
	WEIGHTS_T=1.0d0;WEIGHTS_Q=1.0d0
	end if
	
	
	!$OMP DO SCHEDULE (STATIC)
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
                    
		   damp=lamx
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
!                         IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=(COS(ANGLE1)*SIN(ANGLE2))
				  NY=(SIN(ANGLE1)*SIN(ANGLE2))
				  NZ=(COS(ANGLE2))
				  NNN(1)=NX;NNN(2)=NY;NNN(3)=NZ
				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
					iqp=QP_TRIANGLE_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  do NGP=1,iqp	!for all the gaussian quadrature points
				      CLEFT(1:5)=ILOCAL_RECON3(I)%ULEFT(1:5,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,2,L,NGP);LCVGRAD(2,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,3,L,NGP);
				      LCVGRAD(3,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,4,L,NGP);LCVGRAD(4,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP);
				      CRIGHT(1:5)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:5,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
				      RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,IELEM(N,I)%INEIGHN(L),NGP);
			
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF
					  
					  
					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					  
					  
					  do nvar=1,turbulenceequations+passivescalar			    
					  RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
					  LCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(I)%ULEFTTURBV(1:3,nvar,L,NGP)
					  end do
					   				  
					  
					END IF
			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2) 
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL CONS2PRIM2(N)
					CALL SUTHERLAND(N,LEFTV,RIGHTV)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:6)= LCVGRAD(1,1:3);EDDYFL(7:9)=LCVGRAD(2,1:3)
							  EDDYFL(10:12)=LCVGRAD(3,1:3);EDDYFL(13:15)=LCVGRAD_T(1,1:3)
							  EDDYFL(16:18)=LCVGRAD_T(2,1:3)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:6)= RCVGRAD(1,1:3);EDDYFR(7:9)=RCVGRAD(2,1:3);EDDYFR(10:12)=RCVGRAD(3,1:3)
							  EDDYFR(13:15)=RCVGRAD_T(1,1:3);EDDYFL(16:18)=RCVGRAD_T(2,1:3)
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					  END IF
				       
! 				        TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
! 					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
! 					  U12=ZERO;V12=ZERO;W12=ZERO  
! 					  
! 					   
! 				       				  
! 					  if (turbulence .eq. 1) then
! 					  Q(1:3)=  - OO2* ((LAML(3)*LCVGRAD(4,1:3)) + (LAML(4)*RCVGRAD(4,1:3)))
! 					  else
! 					  Q(1:3) =  - OO2* ((LAML(1)*LCVGRAD(4,1:3)) + (LAML(2)*RCVGRAD(4,1:3)))
! 					  end if
! 					  
! 					  FXV(5) = FXV(5) - Q(1);FYV(5) = FYV(5) - Q(2);FZV(5) = FZV(5) - Q(3)
! 							
! 					 
! 					  !LEFT STATE DERIVATIVES
! 					  UX = LCVGRAD(1,1); UY = LCVGRAD(1,2); UZ = LCVGRAD(1,3);
! 					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2); VZ = LCVGRAD(2,3);
! 					  WX = LCVGRAD(3,1); WY = LCVGRAD(3,2); WZ = LCVGRAD(3,3);
! 					  ! DETERMINE TAUL!!
! 					 
! 
! 					  ! TAU_XX
! 					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
! 					  ! TAU_YY
! 					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
! 					  ! TAU_ZZ
! 					  TAUL(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
! 
! 					  ! tau_xy
! 					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)
! 
! 					  ! TAU_XZ
! 					  TAUL(1,3) = (WX + UZ);TAUL(3,1) = TAUL(1,3)
! 
! 					  ! TAU_YZ
! 					  TAUL(2,3) = (VZ + WY);TAUL(3,2) = TAUL(2,3)
! 					  !END DETERMINE TAUL
! 					  
! 					  
! 					  !RIGHT STATE DERIVATIVES
! 					  UX = RCVGRAD(1,1); UY = RCVGRAD(1,2); UZ = RCVGRAD(1,3);
! 					  VX = RCVGRAD(2,1); VY = RCVGRAD(2,2); VZ = RCVGRAD(2,3);
! 					  WX = RCVGRAD(3,1); WY = RCVGRAD(3,2); WZ = RCVGRAD(3,3);
! 					  ! DETERMINE TAUL!!
! 					  TAUR= ZERO
! 
! 					  ! TAU_XX
! 					  TAUR(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
! 					  ! TAU_YY
! 					  TAUR(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
! 					  ! TAU_ZZ
! 					  TAUR(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
! 
! 					  ! tau_xy
! 					  TAUR(1,2) = (UY + VX);TAUR(2,1) = TAUR(1,2)
! 
! 					  ! TAU_XZ
! 					  TAUR(1,3) = (WX + UZ);TAUR(3,1) = TAUR(1,3)
! 
! 					  ! TAU_YZ
! 					  TAUR(2,3) = (VZ + WY);TAUR(3,2) = TAUR(2,3)
! 					  !END DETERMINE TAUL
! 					  
! 					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
! 					  if ( turbulence .eq. 1) then
! 					    TAU = OO2*(( (VISCL(1)+VISCL(3))*TAUL)+( (VISCL(2)+VISCL(4)) *TAUR))
! 					  else
! 					    TAU = OO2*((VISCL(1)*TAUL)+(VISCL(2)*TAUR))
! 					  end if
! 
! 					    ! NOW ADDITION INTO MOMENTUM FLUXES
! 					  DO KC=2,4
! 					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
! 					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
! 					      FZV(KC) = FZV(KC) + TAU(3,KC-1)
! 					  ENDDO
! 
! 					    ! COMPUTE INTERFACE VELOCITIES
! 					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
! 					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
! 					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
! 					  W12   = OO2*(CLEFT(4)+CRIGHT(4))/RHO12
! 
! 					  FXV(5) = FXV(5) + U12*TAU(1,1) + V12*TAU(1,2) + W12*TAU(1,3)
! 					  FYV(5) = FYV(5) + U12*TAU(2,1) + V12*TAU(2,2) + W12*TAU(2,3)
! 					  FZV(5) = FZV(5) + U12*TAU(3,1) + V12*TAU(3,2) + W12*TAU(3,3)
! 		
! 				       
!                                                     CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
! 		
! !                                             if (b_Code.eq.4)then
! ! 					  damp=zero
! !  					  end if
!                                             
!                                              vdamp=((2*iorder+1)/(ielem(n,i)%dih(L)*1.25))*((4.0/3.0)*OO2*(( (VISCL(1))+(VISCL(2)))))
! 					  HLLCFLUX(1:5)=(NX*FXV+NY*FYV+NZ*FZV)+damp*vdamp*(cright_rot(1:5)-cleft_rot(1:5)) 

                                     TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
					  U12=ZERO;V12=ZERO;W12=ZERO  
					  

! 					  if (b_Code.eq.4)then
! 					  damp=zero
!  					  end if
					  
					  vdamp=(4.0/3.0)!*(( (VISCL(1))+(VISCL(2)))))
                                        nall(1)=nx;nall(2)=ny;nall(3)=nz
                                      LCVGRAD(1,1:3)=((LCVGRAD(1,1:3)+rCVGRAD(1,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(2)-leftv(2)))
					   LCVGRAD(2,1:3)=((LCVGRAD(2,1:3)+rCVGRAD(2,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(3)-leftv(3)))
					  LCVGRAD(3,1:3)=((LCVGRAD(3,1:3)+rCVGRAD(3,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(4)-leftv(4)))
                                        LCVGRAD(4,1:3)=((LCVGRAD(4,1:3)+rCVGRAD(4,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*((rightv(5)/rightv(1))-(leftv(5)/leftv(1))))
				       				  
					  if (turbulence .eq. 1) then
					  Q(1:3)=  - OO2* ((LAML(3)+ (LAML(4)))*LCVGRAD(4,1:3))
					  else
					  Q(1:3) =  - OO2* ((LAML(1)+ (LAML(2)))*LCVGRAD(4,1:3))
					  end if
					  
					  FXV(5) = FXV(5) - Q(1);FYV(5) = FYV(5) - Q(2);FZV(5) = FZV(5) - Q(3)
							
					 
					  !LEFT STATE DERIVATIVES
					 
					  ! DETERMINE TAUL!!
					   UX = LCVGRAD(1,1); UY = LCVGRAD(1,2); UZ = LCVGRAD(1,3);
					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2); VZ = LCVGRAD(2,3);
					  WX = LCVGRAD(3,1); WY = LCVGRAD(3,2); WZ = LCVGRAD(3,3);
					  
					  
					  
                                
					  
					  
					  ! TAU_XX
					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
					  ! TAU_YY
					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
					  ! TAU_ZZ
					  TAUL(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY

					  ! tau_xy
					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)

					  ! TAU_XZ
					  TAUL(1,3) = (WX + UZ);TAUL(3,1) = TAUL(1,3)

					  ! TAU_YZ
					  TAUL(2,3) = (VZ + WY);TAUL(3,2) = TAUL(2,3)
					  !END DETERMINE TAUL
					  
					  
					  
					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
					  if ( turbulence .eq. 1) then
					    TAU = OO2*(( (VISCL(1)+VISCL(3))+(VISCL(2)+VISCL(4))))*TAUL
					  else
					    TAU =OO2*(( (VISCL(1))+(VISCL(2))))*TAUL
					  end if

					    ! NOW ADDITION INTO MOMENTUM FLUXES
					  DO KC=2,4
					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
					      FZV(KC) = FZV(KC) + TAU(3,KC-1)
					  ENDDO

					    ! COMPUTE INTERFACE VELOCITIES
					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
					  W12   = OO2*(CLEFT(4)+CRIGHT(4))/RHO12

					  
								   					  
					  
					  FXV(5) = FXV(5) + U12*TAU(1,1) + V12*TAU(1,2) + W12*TAU(1,3)
					  FYV(5) = FYV(5) + U12*TAU(2,1) + V12*TAU(2,2) + W12*TAU(2,3)
					  FZV(5) = FZV(5) + U12*TAU(3,1) + V12*TAU(3,2) + W12*TAU(3,3)
		
		
					  HLLCFLUX(1:5)=(NX*FXV+NY*FYV+NZ*FZV)	

					  
				      GODFLUX2(1:5)=GODFLUX2(1:5)+(HLLCFLUX(1:5)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  IF (TURBULENCE.EQ.1)THEN
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2)))+(OO2*(VISCL(3)+VISCL(4))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    +((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  ELSE
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    +((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  
					  END IF
						IF (TURBULENCEMODEL.EQ.1)THEN
						HLLCFLUX(6)=HLLCFLUX(6)/SIGMA
						END IF					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    RHS(I)%VAL(1:5)=RHS(I)%VAL(1:5)-GODFLUX2(1:5)
!  				     RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)+GODFLUX2(1:5)
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
!  				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
 				    end if
!  				    end if
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		     
		   
		    DO L=1,IELEM(N,I)%IFCA
                        damp=lamx
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
! 				  IF (icaseb.le.2)THEN
!                                         IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
!                                             igoflux=1
!                                         else
!                                             igoflux=0
!                                         end if
!                                 else
!                                         igoflux=2
!                                 end if
!                                  if (igoflux.ge.1)then
				   damp=lamx  
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				   NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
 				   NNN(1)=NX;NNN(2)=NY;NNN(3)=NZ
 				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
				  else
					iqp=QP_TRIANGLE_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  GODFLUX2=ZERO
				  
				  b_code=0
								  
				  
				  do NGP=1,iqp
				      CLEFT(1:5)=ILOCAL_RECON3(I)%ULEFT(1:5,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,2,L,NGP);LCVGRAD(2,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,3,L,NGP);
				      LCVGRAD(3,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,4,L,NGP);LCVGRAD(4,1:3)=ILOCAL_RECON3(I)%ULEFTV(1:3,1,L,NGP);
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
						do nvar=1,turbulenceequations+passivescalar			    
						LCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(I)%ULEFTTURBV(1:3,nvar,L,NGP)
						end do
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:5)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:5,IELEM(N,I)%INEIGHN(L),NGP)
								  RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,L,NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,L,NGP);
								  RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,L,NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,L,NGP);
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									
									do nvar=1,turbulenceequations+passivescalar			    
									RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
									
									
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
								    
								    
								    RCVGRAD(:,:)=LCVGRAD(:,:)
				  				    RCVGRAD_T(:,:)=LCVGRAD_T(:,:)
				  				     if (B_CODE.eq.4)then	
				  				    rightv=zero
				  				    rightv(2:4)=LCVGRAD(4,1:3)
				  				    leftv=zero
				  				    
				  				    CALL ROTATEF(N,TRI,leftv,rightv,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
								    leftv(2)=-leftv(2)
								    CALL ROTATEB(N,INVTRI,rightv,leftv,ANGLE1,ANGLE2)
				  				    RCVGRAD(4,1:3)=rightv(2:4)
				  				  	end if
! 				  				 
! 				  				  
								  END IF
							ELSE
							      CRIGHT(1:5)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:5,IELEM(N,I)%INEIGHN(L),NGP)
							      RCVGRAD(1,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,3,IELEM(N,I)%INEIGHN(L),NGP);
							      RCVGRAD(3,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,4,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(4,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:3,1,IELEM(N,I)%INEIGHN(L),NGP);
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:NOF_VARIABLES)
									END IF
									do nvar=1,turbulenceequations+passivescalar			    
									RCVGRAD_T(nvar,1:3)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:3,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									  
									  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (IEX.EQ.1)THEN
									  RCVGRAD(4,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)	      
										      END IF
										END DO  
									  END DO
									  
									 

									  
									  
									  
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								    END IF
									  
									  
									  
									  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1 
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT) 
										END DO  
									  END DO
									  
									  
									  
									  

								END IF
							ELSE 			
							
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (IEX.EQ.1)THEN
									  RCVGRAD(4,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)	      
										      END IF
										END DO  
									  END DO
								  
								  
								  
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
									
									 DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1 
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT) 
										END DO  
									  END DO
									
									
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				    
			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:5)=CLEFT_ROT(1:5); RIGHTV(1:5)=CRIGHT_ROT(1:5)
						  CALL LMACHT(N)
						  CLEFT_ROT(1:5)=LEFTV(1:5);CRIGHT_ROT(1:5)=RIGHTV(1:5);
						  CALL ROTATEB(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2) 
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL CONS2PRIM2(N)
					CALL SUTHERLAND(N,LEFTV,RIGHTV)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:6)= LCVGRAD(1,1:3);EDDYFL(7:9)=LCVGRAD(2,1:3)
							  EDDYFL(10:12)=LCVGRAD(3,1:3);EDDYFL(13:15)=LCVGRAD_T(1,1:3)
							  EDDYFL(16:18)=LCVGRAD_T(2,1:3)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:6)= RCVGRAD(1,1:3);EDDYFR(7:9)=RCVGRAD(2,1:3);EDDYFR(10:12)=RCVGRAD(3,1:3)
							  EDDYFR(13:15)=RCVGRAD_T(1,1:3);EDDYFL(16:18)=RCVGRAD_T(2,1:3)
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					  END IF
				       
				       
				       
				       
! 					 TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
! 					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
! 					  U12=ZERO;V12=ZERO;W12=ZERO  
! 					  
! 					   
! 				       				  
! 					  if (turbulence .eq. 1) then
! 					  Q(1:3)=  - OO2* ((LAML(3)*LCVGRAD(4,1:3)) + (LAML(4)*RCVGRAD(4,1:3)))
! 					  else
! 					  Q(1:3) =  - OO2* ((LAML(1)*LCVGRAD(4,1:3)) + (LAML(2)*RCVGRAD(4,1:3)))
! 					  end if
! 					  
! 					  FXV(5) = FXV(5) - Q(1);FYV(5) = FYV(5) - Q(2);FZV(5) = FZV(5) - Q(3)
! 							
! 					 
! 					  !LEFT STATE DERIVATIVES
! 					  UX = LCVGRAD(1,1); UY = LCVGRAD(1,2); UZ = LCVGRAD(1,3);
! 					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2); VZ = LCVGRAD(2,3);
! 					  WX = LCVGRAD(3,1); WY = LCVGRAD(3,2); WZ = LCVGRAD(3,3);
! 					  ! DETERMINE TAUL!!
! 					 
! 
! 					  ! TAU_XX
! 					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
! 					  ! TAU_YY
! 					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
! 					  ! TAU_ZZ
! 					  TAUL(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
! 
! 					  ! tau_xy
! 					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)
! 
! 					  ! TAU_XZ
! 					  TAUL(1,3) = (WX + UZ);TAUL(3,1) = TAUL(1,3)
! 
! 					  ! TAU_YZ
! 					  TAUL(2,3) = (VZ + WY);TAUL(3,2) = TAUL(2,3)
! 					  !END DETERMINE TAUL
! 					  
! 					  
! 					  !RIGHT STATE DERIVATIVES
! 					  UX = RCVGRAD(1,1); UY = RCVGRAD(1,2); UZ = RCVGRAD(1,3);
! 					  VX = RCVGRAD(2,1); VY = RCVGRAD(2,2); VZ = RCVGRAD(2,3);
! 					  WX = RCVGRAD(3,1); WY = RCVGRAD(3,2); WZ = RCVGRAD(3,3);
! 					  ! DETERMINE TAUL!!
! 					  TAUR= ZERO
! 
! 					  ! TAU_XX
! 					  TAUR(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
! 					  ! TAU_YY
! 					  TAUR(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
! 					  ! TAU_ZZ
! 					  TAUR(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
! 
! 					  ! tau_xy
! 					  TAUR(1,2) = (UY + VX);TAUR(2,1) = TAUR(1,2)
! 
! 					  ! TAU_XZ
! 					  TAUR(1,3) = (WX + UZ);TAUR(3,1) = TAUR(1,3)
! 
! 					  ! TAU_YZ
! 					  TAUR(2,3) = (VZ + WY);TAUR(3,2) = TAUR(2,3)
! 					  !END DETERMINE TAUL
! 					  
! 					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
! 					  if ( turbulence .eq. 1) then
! 					    TAU = OO2*(( (VISCL(1)+VISCL(3))*TAUL)+( (VISCL(2)+VISCL(4)) *TAUR))
! 					  else
! 					    TAU = OO2*((VISCL(1)*TAUL)+(VISCL(2)*TAUR))
! 					  end if
! 
! 					    ! NOW ADDITION INTO MOMENTUM FLUXES
! 					  DO KC=2,4
! 					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
! 					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
! 					      FZV(KC) = FZV(KC) + TAU(3,KC-1)
! 					  ENDDO
! 
! 					    ! COMPUTE INTERFACE VELOCITIES
! 					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
! 					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
! 					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
! 					  W12   = OO2*(CLEFT(4)+CRIGHT(4))/RHO12
! 
! 					  FXV(5) = FXV(5) + U12*TAU(1,1) + V12*TAU(1,2) + W12*TAU(1,3)
! 					  FYV(5) = FYV(5) + U12*TAU(2,1) + V12*TAU(2,2) + W12*TAU(2,3)
! 					  FZV(5) = FZV(5) + U12*TAU(3,1) + V12*TAU(3,2) + W12*TAU(3,3)
! 		
! 		
!                                             CALL ROTATEF(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						  CALL ROTATEF(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
! 		
!                                             if (b_Code.eq.4)then
! 					  damp=zero
!  					  end if
!                                             
!                                              vdamp=((2*iorder+1)/(ielem(n,i)%dih(L)*1.25))*((4.0/3.0)*OO2*(( (VISCL(1))+(VISCL(2)))))
! 					  HLLCFLUX(1:5)=(NX*FXV+NY*FYV+NZ*FZV)+damp*vdamp*(cright_rot(1:5)-cleft_rot(1:5))		

                                            
					  TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
					  U12=ZERO;V12=ZERO;W12=ZERO  
					  

					 if ((b_Code.lt.5).and.(b_Code.gt.0))then
					  damp=zero
 					  end if
					  
					  vdamp=(4.0/3.0)!*(( (VISCL(1))+(VISCL(2)))))
                                        nall(1)=nx;nall(2)=ny;nall(3)=nz
                                      LCVGRAD(1,1:3)=((LCVGRAD(1,1:3)+rCVGRAD(1,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(2)-leftv(2)))
					   LCVGRAD(2,1:3)=((LCVGRAD(2,1:3)+rCVGRAD(2,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(3)-leftv(3)))
					  LCVGRAD(3,1:3)=((LCVGRAD(3,1:3)+rCVGRAD(3,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(4)-leftv(4)))
                                        LCVGRAD(4,1:3)=((LCVGRAD(4,1:3)+rCVGRAD(4,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*((rightv(5)/rightv(1))-(leftv(5)/leftv(1))))
				       				  
					  if (turbulence .eq. 1) then
					  Q(1:3)=  - OO2* ((LAML(3)+ (LAML(4)))*LCVGRAD(4,1:3))
					  else
					  Q(1:3) =  - OO2* ((LAML(1)+ (LAML(2)))*LCVGRAD(4,1:3))
					  end if
					  
					  FXV(5) = FXV(5) - Q(1);FYV(5) = FYV(5) - Q(2);FZV(5) = FZV(5) - Q(3)
							
					 
					  !LEFT STATE DERIVATIVES
					 
					  ! DETERMINE TAUL!!
					   UX = LCVGRAD(1,1); UY = LCVGRAD(1,2); UZ = LCVGRAD(1,3);
					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2); VZ = LCVGRAD(2,3);
					  WX = LCVGRAD(3,1); WY = LCVGRAD(3,2); WZ = LCVGRAD(3,3);
					  
					  
					  
                                
					  
					  
					  ! TAU_XX
					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
					  ! TAU_YY
					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
					  ! TAU_ZZ
					  TAUL(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY

					  ! tau_xy
					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)

					  ! TAU_XZ
					  TAUL(1,3) = (WX + UZ);TAUL(3,1) = TAUL(1,3)

					  ! TAU_YZ
					  TAUL(2,3) = (VZ + WY);TAUL(3,2) = TAUL(2,3)
					  !END DETERMINE TAUL
					  
					  
					  
					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
					  if ( turbulence .eq. 1) then
					    TAU = OO2*(( (VISCL(1)+VISCL(3))+(VISCL(2)+VISCL(4))))*TAUL
					  else
					    TAU =OO2*(( (VISCL(1))+(VISCL(2))))*TAUL
					  end if

					    ! NOW ADDITION INTO MOMENTUM FLUXES
					  DO KC=2,4
					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
					      FZV(KC) = FZV(KC) + TAU(3,KC-1)
					  ENDDO

					    ! COMPUTE INTERFACE VELOCITIES
					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
					  W12   = OO2*(CLEFT(4)+CRIGHT(4))/RHO12

					  
								   					  
					  
					  FXV(5) = FXV(5) + U12*TAU(1,1) + V12*TAU(1,2) + W12*TAU(1,3)
					  FYV(5) = FYV(5) + U12*TAU(2,1) + V12*TAU(2,2) + W12*TAU(2,3)
					  FZV(5) = FZV(5) + U12*TAU(3,1) + V12*TAU(3,2) + W12*TAU(3,3)
		
		
					  HLLCFLUX(1:5)=(NX*FXV+NY*FYV+NZ*FZV)	


				      GODFLUX2(1:5)=GODFLUX2(1:5)+(HLLCFLUX(1:5)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				      
				      
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					   IF (TURBULENCE.EQ.1)THEN
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2)))+(OO2*(VISCL(3)+VISCL(4))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    +((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  ELSE
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    +((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  
					  END IF
						IF (TURBULENCEMODEL.EQ.1)THEN
						HLLCFLUX(6)=HLLCFLUX(6)/SIGMA
						END IF					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    
				    RHS(I)%VAL(1:5)=RHS(I)%VAL(1:5)-GODFLUX2(1:5)
! 				    if ((igoflux.eq.1))then
! 				    RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:5)+GODFLUX2(1:5)
! 				    end if
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				     if ((igoflux.eq.1))then
! 				   RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
!  				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    end if
				    
! 				    end if
				    end if
		    END DO
	END DO
	!$OMP END DO

END SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE




SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE2d(N)
!> @brief
!> This subroutine computes the diffusive fluxes for Euler-Navier-Stokes Equations in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2
	INTEGER::I,L,NGP,KMAXE,IQP,ii,NVAR,KC,IEX,ITTT,IKAS,igoflux, icaseb
	REAL::sum_detect,NORMS
	REAL,DIMENSION(NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP
	REAL,DIMENSION(3,2)::LCVGRAD,RCVGRAD
	REAL,DIMENSION(TURBULENCEEQUATIONS+PASSIVESCALAR,2)::LCVGRAD_T,RCVGRAD_T
	real,dimension(4)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(2,2)::taul,taur,TAU
	REAL,DIMENSION(2)::Q,NALL
	REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,RHO12,U12,V12,W12,damp,vdamp  
	
	
	
	
	
	
	
	
	
	
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREline(N,IGQRULES)
	
	WEIGHTS_temp(1:QP_line_n)=WEQUA2D(1:QP_line_n)
	if (Reduce_comp.eq.1)then
	WEIGHTS_temp=1.0d0
	end if
	
	!$OMP DO SCHEDULE (STATIC)
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
		   
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
!                              IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
                                damp=lamx
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=angle1
				  NY=angle2
				 
				  
					iqp=qp_line_n
				
				  
				  do NGP=1,iqp	!for all the gaussian quadrature points
				      CLEFT(1:4)=ILOCAL_RECON3(I)%ULEFT(1:4,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
				      CRIGHT(1:4)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:4,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
				      RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
				      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP)
			
					IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.1)then
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp) !left additional equations flow state
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
					  ELSE
					    CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					    CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
					  END IF
					  
					  
					  cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
					  cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
					  
					  
					  do nvar=1,turbulenceequations+passivescalar			    
					  RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
					  LCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,nvar,L,NGP)
					  end do
					   				  
					   
					END IF
			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						  CALL ROTATEB2d(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB2d(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2) 
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL CONS2PRIM2d2(N)
					CALL SUTHERLAND2D(N,LEFTV,RIGHTV)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= LCVGRAD(1,1:2);EDDYFL(6:7)=LCVGRAD(2,1:2)
							  EDDYFL(8:9)=LCVGRAD_T(1,1:2)
							  EDDYFL(10:11)=LCVGRAD_T(2,1:2)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:5)= RCVGRAD(1,1:2);EDDYFR(6:7)=RCVGRAD(2,1:2)
							  EDDYFR(8:9)=RCVGRAD_T(1,1:2);EDDYFL(10:11)=RCVGRAD_T(2,1:2)
							    Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					  END IF
				       
				       
				       
				       
! 					   TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
! 					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
! 					  U12=ZERO;V12=ZERO;W12=ZERO  
! 					  
! 					   
! 				       				  
! 					  if (turbulence .eq. 1) then
! 					  Q(1:2)=  - OO2* ((LAML(3)*LCVGRAD(3,1:2)) + (LAML(4)*RCVGRAD(3,1:2)))
! 					  else
! 					  Q(1:2) =  - OO2* ((LAML(1)*LCVGRAD(3,1:2)) + (LAML(2)*RCVGRAD(3,1:2)))
! 					  end if
! 					  
! 					  FXV(4) = FXV(4) - Q(1);FYV(4) = FYV(4) - Q(2)
! 							
! 					 
! 					  !LEFT STATE DERIVATIVES
! 					  UX = LCVGRAD(1,1); UY = LCVGRAD(1,2)
! 					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2)
! 					  ! DETERMINE TAUL!!
! 					 
! 
! 					  ! TAU_XX
! 					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY 
! 					  ! TAU_YY
! 					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX 
! 					  ! TAU_ZZ
! 					 
! 
! 					  ! tau_xy
! 					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)
! 
! 					  
! 					  !END DETERMINE TAUL
! 					  
! 					  
! 					  !RIGHT STATE DERIVATIVES
! 					  UX = RCVGRAD(1,1); UY = RCVGRAD(1,2)
! 					  VX = RCVGRAD(2,1); VY = RCVGRAD(2,2)
! 					  
! 					  ! DETERMINE TAUL!!
! 					  TAUR= ZERO
! 
! 					  ! TAU_XX
! 					  TAUR(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY 
! 					  ! TAU_YY
! 					  TAUR(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX 
! 					  
! 
! 					  ! tau_xy
! 					  TAUR(1,2) = (UY + VX);TAUR(2,1) = TAUR(1,2)
! 
! 					  
! 					  !END DETERMINE TAUL
! 					  
! 					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
! 					  if ( turbulence .eq. 1) then
! 					    TAU = OO2*(( (VISCL(1)+VISCL(3))*TAUL)+( (VISCL(2)+VISCL(4)) *TAUR))
! 					  else
! 					    TAU = OO2*((VISCL(1)*TAUL)+(VISCL(2)*TAUR))
! 					  end if
! 
! 					    ! NOW ADDITION INTO MOMENTUM FLUXES
! 					  DO KC=2,3
! 					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
! 					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
! 					      
! 					  ENDDO
! 
! 					    ! COMPUTE INTERFACE VELOCITIES
! 					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
! 					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
! 					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
! 					  
! 
! 					  FXV(4) = FXV(4) + U12*TAU(1,1) + V12*TAU(1,2) 
! 					  FYV(4) = FYV(4) + U12*TAU(2,1) + V12*TAU(2,2) 
! 					 
! 		
! 		
!                                CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						  CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
! ! 		
! !                                             if (b_Code.eq.4)then
! ! 					  damp=zero
! !  					  end if
!                                             
!                                              vdamp=((2*iorder+1)/(ielem(n,i)%dih(L)*1.25))*((4.0/3.0)*OO2*(( (VISCL(1))+(VISCL(2)))))
! ! 					  HLLCFLUX(1:5)=(NX*FXV+NY*FYV+NZ*FZV)+damp*vdamp*(cright_rot(1:5)-cleft_rot(1:5)) 
! 					  HLLCFLUX(1:4)=(NX*FXV+NY*FYV)+damp*vdamp*(cright_rot(1:4)-cleft_rot(1:4))	


                                     TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
					  U12=ZERO;V12=ZERO;W12=ZERO 
				       
! 				      
					  
					  vdamp=(4.0/3.0)!*(( (VISCL(1))+(VISCL(2)))))
                                        nall(1)=nx;nall(2)=ny
                                      LCVGRAD(1,1:2)=((LCVGRAD(1,1:2)+rCVGRAD(1,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*(rightv(2)-leftv(2)))
					   LCVGRAD(2,1:2)=((LCVGRAD(2,1:2)+rCVGRAD(2,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*(rightv(3)-leftv(3)))
					  LCVGRAD(3,1:2)=((LCVGRAD(3,1:2)+rCVGRAD(3,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*((rightv(4)/rightv(1))-(leftv(4)/leftv(1))))
                                        
				       
					  
					  
					  
					  
					   
				       				  
					  if (turbulence .eq. 1) then
					  Q(1:2)=  - OO2* ((LAML(3) +(LAML(4)))*lCVGRAD(3,1:2))
					  else
					  Q(1:2)=  - OO2* ((LAML(1) +(LAML(2)))*lCVGRAD(3,1:2))
					  end if
					  
					  FXV(4) = FXV(4) - Q(1);FYV(4) = FYV(4) - Q(2)
							
					 
					  !LEFT STATE DERIVATIVES
					  UX = LCVGRAD(1,1); UY = LCVGRAD(1,2)
					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2)
					  ! DETERMINE TAUL!!
					 

					  ! TAU_XX
					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY 
					  ! TAU_YY
					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX 
					  ! TAU_ZZ
					 

					  ! tau_xy
					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)

					  
					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
					  if ( turbulence .eq. 1) then
					    TAU = OO2*(( (VISCL(1)+VISCL(3)))+( (VISCL(2)+VISCL(4))))*taul
					  else
					    TAU = OO2*((VISCL(1))+(VISCL(2)))*taul
					  end if

					    ! NOW ADDITION INTO MOMENTUM FLUXES
					  DO KC=2,3
					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
					      
					  ENDDO

					    ! COMPUTE INTERFACE VELOCITIES
					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
					  

					  FXV(4) = FXV(4) + U12*TAU(1,1) + V12*TAU(1,2) 
					  FYV(4) = FYV(4) + U12*TAU(2,1) + V12*TAU(2,2) 
! 					 
					  HLLCFLUX(1:4)=(NX*FXV+NY*FYV)	
		
					  			      
				      GODFLUX2(1:4)=GODFLUX2(1:4)+(HLLCFLUX(1:4)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  IF (TURBULENCE.EQ.1)THEN
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2)))+(OO2*(VISCL(3)+VISCL(4))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY))
					  ELSE
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY))
					  
					  END IF
						IF (TURBULENCEMODEL.EQ.1)THEN
						HLLCFLUX(5)=HLLCFLUX(5)/SIGMA
						END IF					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    RHS(I)%VAL(1:4)=RHS(I)%VAL(1:4)-GODFLUX2(1:4)
!  				     RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)+GODFLUX2(1:4)
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
!  				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
!  				    end if
 				    end if
		    END DO
	END DO
	!$OMP END DO
	
	
	!$OMP DO SCHEDULE (STATIC) 
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
				
		    
		    
		    DO L=1,IELEM(N,I)%IFCA
                            damp=lamx
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
! 				  IF (icaseb.le.2)THEN
!                                         IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
!                                             igoflux=1
!                                         else
!                                             igoflux=0
!                                         end if
!                                 else
!                                         igoflux=2
!                                 end if
!                                  if (igoflux.ge.1)then
                                    damp=lamx
				      b_code=0
				 GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=angle1
				  NY=angle2
				 
				  
					iqp=qp_line_n
				
				  
				 
								  
				  
				  do NGP=1,iqp
				   damp=lamx
				      CLEFT(1:4)=ILOCAL_RECON3(I)%ULEFT(1:4,L,NGP)	!left mean flow state
				      LCVGRAD(1,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,2,L,NGP);LCVGRAD(2,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,3,L,NGP);
				      LCVGRAD(3,1:2)=ILOCAL_RECON3(I)%ULEFTV(1:2,1,L,NGP)
					 IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
						if (icoupleturb.eq.1)then
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(I)%ULEFTTURB(1:turbulenceequations+PASSIVESCALAR,L,ngp)
						ELSE
							CTURBL(1:turbulenceequations+PASSIVESCALAR)=U_CT(I)%VAL(1,1:turbulenceequations+PASSIVESCALAR)
						end if
						do nvar=1,turbulenceequations+passivescalar			    
						LCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(I)%ULEFTTURBV(1:2,nvar,L,NGP)
						end do
					end if
				      
				      
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								  CRIGHT(1:4)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:4,IELEM(N,I)%INEIGHN(L),NGP)
								  RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,L,NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,L,NGP);
								  RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,L,NGP);
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									
									do nvar=1,turbulenceequations+passivescalar			    
									RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
									
									
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL coordinates_face_inner2d(N,Iconsidered,facex)
								    CORDS(1:2)=zero
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED)
								    cright(1:4)=rightv(1:4)
								    
								    RCVGRAD(:,:)=LCVGRAD(:,:)
				  				    RCVGRAD_T(:,:)=LCVGRAD_T(:,:)
				  				     if (B_CODE.eq.4)then	
				  				    rightv=zero
				  				    rightv(2:3)=LCVGRAD(3,1:2)
				  				    leftv=zero
				  				    
				  				    CALL ROTATEF2d(N,TRI,leftv,rightv,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
								    leftv(2)=-leftv(2)
								    CALL ROTATEB2d(N,INVTRI,rightv,leftv,ANGLE1,ANGLE2)
				  				    RCVGRAD(3,1:2)=rightv(2:3)
				  				  	end if
! 				  				 
! 				  				  
								  END IF
							ELSE
							      CRIGHT(1:4)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:4,IELEM(N,I)%INEIGHN(L),NGP)
							      RCVGRAD(1,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,2,IELEM(N,I)%INEIGHN(L),NGP);RCVGRAD(2,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,3,IELEM(N,I)%INEIGHN(L),NGP);
							      RCVGRAD(3,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTV(1:2,1,IELEM(N,I)%INEIGHN(L),NGP);
								  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
									do nvar=1,turbulenceequations+passivescalar			    
									RCVGRAD_T(nvar,1:2)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURBV(1:2,nvar,IELEM(N,I)%INEIGHN(L),NGP)
									end do
								    END IF
							      
							      
							      
							      
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
					    
					    
					     
					    
					    
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
									  
									  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (IEX.EQ.1)THEN
									  RCVGRAD(3,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)	      
										      END IF
										END DO  
									  END DO
									  
									 

									  
									  
									  
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
								   
									  
									  
									  
									  DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1 
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT) 
										END DO  
									  END DO
									  
								END IF
									  
									  

								END IF
							ELSE 			
							
								  CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								  ITTT=0
									  DO IEX=1,NOF_VARIABLES-1
										DO nvar=1,DIMS
										      ITTT=ITTT+1
										      IF (IEX.EQ.1)THEN
									  RCVGRAD(3,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)
										      ELSE
									  RCVGRAD(IEX-1,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT)	      
										      END IF
										END DO  
									  END DO
								  
								  
								  
								   IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL&
									   (IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),nof_variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)!right additional equations flow state
									END IF
									
									 DO IEX=1,TURBULENCEEQUATIONS+PASSIVESCALAR
										DO nvar=1,DIMS
										      ITTT=ITTT+1 
									  RCVGRAD_T(IEX,NVAR)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR+ITTT) 
										END DO  
									  END DO
									
									
								    END IF
								  
! 								   
							END IF
					    END IF
				      
				    
			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:4)=CLEFT_ROT(1:4); RIGHTV(1:4)=CRIGHT_ROT(1:4)
						  CALL LMACHT2d(N)
						  CLEFT_ROT(1:4)=LEFTV(1:4);CRIGHT_ROT(1:4)=RIGHTV(1:4);
						  CALL ROTATEB2d(N,INVTRI,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB2d(N,INVTRI,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2) 
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL CONS2PRIM2d2(N)
					CALL SUTHERLAND2d(N,LEFTV,RIGHTV)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      IF (B_CODE.EQ.4)THEN
						      VISCL(3:4)=ZERo
						      LAML(3:4)=ZERO
						      END IF
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= LCVGRAD(1,1:2);EDDYFL(6:7)=LCVGRAD(2,1:2)
							  EDDYFL(8:9)=LCVGRAD_T(1,1:2)
							  EDDYFL(10:11)=LCVGRAD_T(2,1:2)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:5)= RCVGRAD(1,1:2);EDDYFR(6:7)=RCVGRAD(2,1:2)
							  EDDYFR(8:9)=RCVGRAD_T(1,1:2);EDDYFL(10:11)=RCVGRAD_T(2,1:2)
							    Call EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
						      END IF
					  END IF
				       
				       
				       
				       
! 					   TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
! 					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
! 					  U12=ZERO;V12=ZERO;W12=ZERO  
! 					  
! 					   
! 				       				  
! 					 if (turbulence .eq. 1) then
! 					  Q(1:2)=  - OO2* ((LAML(3)*LCVGRAD(3,1:2)) + (LAML(4)*RCVGRAD(3,1:2)))
! 					  else
! 					  Q(1:2) =  - OO2* ((LAML(1)*LCVGRAD(3,1:2)) + (LAML(2)*RCVGRAD(3,1:2)))
! 					  end if
! 					  
! 					  FXV(4) = FXV(4) - Q(1);FYV(4) = FYV(4) - Q(2)
! 							
! 					 
! 					  !LEFT STATE DERIVATIVES
! 					  UX = LCVGRAD(1,1); UY = LCVGRAD(1,2)
! 					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2)
! 					  ! DETERMINE TAUL!!
! 					 
! 
! 					  ! TAU_XX
! 					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY 
! 					  ! TAU_YY
! 					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX 
! 					  ! TAU_ZZ
! 					 
! 
! 					  ! tau_xy
! 					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)
! 
! 					  
! 					  !END DETERMINE TAUL
! 					  
! 					  
! 					  !RIGHT STATE DERIVATIVES
! 					  UX = RCVGRAD(1,1); UY = RCVGRAD(1,2)
! 					  VX = RCVGRAD(2,1); VY = RCVGRAD(2,2)
! 					  
! 					  ! DETERMINE TAUL!!
! 					  TAUR= ZERO
! 
! 					  ! TAU_XX
! 					  TAUR(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY 
! 					  ! TAU_YY
! 					  TAUR(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX 
! 					  
! 
! 					  ! tau_xy
! 					  TAUR(1,2) = (UY + VX);TAUR(2,1) = TAUR(1,2)
! 
! 					  
! 					  !END DETERMINE TAUL
! 					  
! 					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
! 					  if ( turbulence .eq. 1) then
! 					    TAU = OO2*(( (VISCL(1)+VISCL(3))*TAUL)+( (VISCL(2)+VISCL(4)) *TAUR))
! 					  else
! 					    TAU = OO2*((VISCL(1)*TAUL)+(VISCL(2)*TAUR))
! 					  end if
! 
! 					    ! NOW ADDITION INTO MOMENTUM FLUXES
! 					  DO KC=2,3
! 					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
! 					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
! 					      
! 					  ENDDO
! 
! 					    ! COMPUTE INTERFACE VELOCITIES
! 					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
! 					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
! 					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
! 					  
! 
! 					  FXV(4) = FXV(4) + U12*TAU(1,1) + V12*TAU(1,2) 
! 					  FYV(4) = FYV(4) + U12*TAU(2,1) + V12*TAU(2,2) 
! 					  
! 					 
!                                                     CALL ROTATEF2d(N,TRI,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						  CALL ROTATEF2d(N,TRI,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)
! 		
!                                             if (b_Code.eq.4)then
! 					  damp=zero
!  					  end if
!                                             
!                                              vdamp=((2*iorder+1)/(ielem(n,i)%dih(L)*1.25))*((4.0/3.0)*OO2*(( (VISCL(1))+(VISCL(2)))))


                                          TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
					  U12=ZERO;V12=ZERO;W12=ZERO 
				       
				      if ((b_Code.lt.5).and.(b_Code.gt.0))then
					  damp=zero
 					  end if
				       
					   vdamp=(4.0/3.0)!*(( (VISCL(1))+(VISCL(2)))))
                                        nall(1)=nx;nall(2)=ny
                                      LCVGRAD(1,1:2)=((LCVGRAD(1,1:2)+rCVGRAD(1,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*(rightv(2)-leftv(2)))
					   LCVGRAD(2,1:2)=((LCVGRAD(2,1:2)+rCVGRAD(2,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*(rightv(3)-leftv(3)))
					  LCVGRAD(3,1:2)=((LCVGRAD(3,1:2)+rCVGRAD(3,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*((rightv(4)/rightv(1))-(leftv(4)/leftv(1))))
					  
					  
					  
					  
					   
				       				  
					  if (turbulence .eq. 1) then
					  Q(1:2)=  - OO2* ((LAML(3) +(LAML(4)))*lCVGRAD(3,1:2))
					  else
					  Q(1:2)=  - OO2* ((LAML(1) +(LAML(2)))*lCVGRAD(3,1:2))
					  end if
					  
					  FXV(4) = FXV(4) - Q(1);FYV(4) = FYV(4) - Q(2)
							
					 
					  !LEFT STATE DERIVATIVES
					  UX = LCVGRAD(1,1); UY = LCVGRAD(1,2)
					  VX = LCVGRAD(2,1); VY = LCVGRAD(2,2)
					  ! DETERMINE TAUL!!
					 

					  ! TAU_XX
					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY 
					  ! TAU_YY
					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX 
					  ! TAU_ZZ
					 

					  ! tau_xy
					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)

					  
					 ! AVERAGE AND MULTIPLAY BY VISCOSITY
					  if ( turbulence .eq. 1) then
					    TAU = OO2*(( (VISCL(1)+VISCL(3)))+( (VISCL(2)+VISCL(4))))*taul
					  else
					    TAU = OO2*((VISCL(1))+(VISCL(2)))*taul
					  end if

					    ! NOW ADDITION INTO MOMENTUM FLUXES
					  DO KC=2,3
					      FXV(KC) = FXV(KC) + TAU(1,KC-1)
					      FYV(KC) = FYV(KC) + TAU(2,KC-1)
					      
					  ENDDO

					    ! COMPUTE INTERFACE VELOCITIES
					  RHO12 = OO2*(CLEFT(1)+CRIGHT(1))
					  U12   = OO2*(CLEFT(2)+CRIGHT(2))/RHO12
					  V12   = OO2*(CLEFT(3)+CRIGHT(3))/RHO12
					  

					  FXV(4) = FXV(4) + U12*TAU(1,1) + V12*TAU(1,2) 
					  FYV(4) = FYV(4) + U12*TAU(2,1) + V12*TAU(2,2) 
! 					 
					  HLLCFLUX(1:4)=(NX*FXV+NY*FYV)			
					  
		
		
			      
				      GODFLUX2(1:4)=GODFLUX2(1:4)+(HLLCFLUX(1:4)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  IF (TURBULENCE.EQ.1)THEN
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2)))+(OO2*(VISCL(3)+VISCL(4))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY))
					  ELSE
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY))
					  
					  END IF
						IF (TURBULENCEMODEL.EQ.1)THEN
						HLLCFLUX(5)=HLLCFLUX(5)/SIGMA
						END IF					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    
				     RHS(I)%VAL(1:4)=RHS(I)%VAL(1:4)-GODFLUX2(1:4)
! 				    if ((igoflux.eq.1))then
! 				    RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:4)+GODFLUX2(1:4)
! 				    end if
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				     if ((igoflux.eq.1))then
! 				     RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
!  				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    end if
				    
! 				    end if
				    end if
		    END DO
	END DO
	!$OMP END DO

END SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE2d











END MODULE FLUXES
