MODULE FLUXES
USE LIBRARY
USE TRANSFORM
USE LOCAL
USE RIEMANN
USE FLOW_OPERATIONS
USE DECLARATION
USE BASIS
USE DG_FUNCTIONS
IMPLICIT NONE

 CONTAINS
 
 
SUBROUTINE SOLUTION_INTEG(I,SOLUTION_INTEG2)
 IMPLICIT NONE
 REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::SOLUTION_INTEG2
 INTEGER,INTENT(IN)::I
 
 SOLUTION_INTEG2=DG_VOL_INTEGRAL2(N,I)
 
 
 END SUBROUTINE SOLUTION_INTEG

 SUBROUTINE SOLUTION_INTEG_S(I,SOLUTION_INTEG_STRONG)
 IMPLICIT NONE
  REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::SOLUTION_INTEG_STRONG
 INTEGER,INTENT(IN)::I

 SOLUTION_INTEG_STRONG=DG_VOL_INTEGRAL_STRONG(N,I)


 END SUBROUTINE SOLUTION_INTEG_S


  SUBROUTINE SOLUTION_INTEG_W(I,SOLUTION_INTEG_WEAK)
 IMPLICIT NONE
 REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::SOLUTION_INTEG_WEAK
 INTEGER,INTENT(IN)::I

 SOLUTION_INTEG_WEAK=DG_VOL_INTEGRAL_WEAK(N,I)


 END SUBROUTINE SOLUTION_INTEG_W




SUBROUTINE CALCULATE_FLUXESHI(N)
!> @brief
!> This subroutine computes the fluxes for linear-advection equation
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL::GODFLUX2,sum_detect
	INTEGER::I,L,NGP,KMAXE,IQP
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP,WEIGHTS_DG
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    REAL::ANGLE1,ANGLE2,NX,NY,NZ,NORMALVECT
    INTEGER::facex,POINTX,ICONSIDERED
    REAL,DIMENSION(1:NOF_VARIABLES)::CLEFT,CRIGHT,HLLCFLUX,RHLLCFLUX
    REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG
	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	
	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)
	
	call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)

	IF (DG.EQ.1)THEN
	allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))

	END IF



	!$OMP DO
	DO I=1,KMAXE
	
        iconsidered=i
        if (dg.eq.1)then
        RHS(I)%VALDG = ZERO
        DG_RHS = ZERO
        DG_RHS_SURF_INTEG = ZERO
        DG_RHS_VOL_INTEG = ZERO
        end if
	
        IF (dg.EQ.1) THEN
            DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
            
        END IF
	
				  
		IF (IELEM(N,I)%INTERIOR.EQ.0)THEN
		    RHS(I)%VAL(1)=ZERO
		    
		    DO L=1,IELEM(N,I)%IFCA
				  GODFLUX2=ZERO
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				  facex=l
				  NORMALVECT=((COS(ANGLE1)*SIN(ANGLE2))*LAMX)+((SIN(ANGLE1)*SIN(ANGLE2))*LAMY)+((COS(ANGLE2))*LAMZ)
				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
					
				  else
					iqp=QP_TRIANGLE
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
				  end if
				  
				  if (dg.eq.1)WEIGHTS_dg(1:iqp)=WEIGHTS_TEMP(1:IQP)
				  
				  
				  
				  do NGP=1,iqp
				  pointx=ngp
				  IF (dg.EQ.1) THEN
                        CLEFT = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                    ELSE !FV
				  
				  
				  
				      CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
				      CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                    end if
				      CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)
				      if (dg.eq.1)then
				      RHLLCFLUX(1)=HLLCFLUX(1)
				      
				      DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
				      
				       
				      
				      else
				      GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      end if
				  END DO
				    RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2

		    END DO
		END IF
		IF (IELEM(N,I)%INTERIOR.EQ.1)THEN
		    RHS(I)%VAL(1)=ZERO
		    
		    DO L=1,IELEM(N,I)%IFCA
				      facex=l
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
				  
				  
				   if (dg.eq.1)WEIGHTS_dg(1:iqp)=WEIGHTS_TEMP(1:IQP)
				  GODFLUX2=ZERO
				  do NGP=1,iqp
				  POINTX = NGP
                    IF (dg == 1) THEN
                        CLEFT(1:nof_variables) = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        else
				      CLEFT(1:nof_variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_variables,L,NGP)
				      end if
					    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                                 
								  if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
								 IF (dg == 1) THEN
                                    CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE  
                                CRIGHT(1:nof_variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_variables,IELEM(N,I)%INEIGHN(L),NGP)
								end if
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								  CRIGHT(1:nof_variables)=CLEFT(1:nof_variables)
								    
								  END IF
							ELSE
							     IF (dg == 1) THEN
                                CRIGHT(1:NOF_VARIABLES) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)

                            ELSE  
                            CRIGHT(1:nof_variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_variables,IELEM(N,I)%INEIGHN(L),NGP)
                            end if
							END IF
					    ELSE	!IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
						
							IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
								if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
									 IF (dg == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                else 
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
								end if
								END IF
							ELSE 								
								 IF (dg == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                else
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                end if
! 								   
							END IF
					    END IF
				      
				      CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)
				      
				      if (dg.eq.1)then
				      RHLLCFLUX(1)=HLLCFLUX(1)
				      DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
				      
				       
				     
				      else
				      
				      GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))	
				      
				      end if
				  END DO
				    RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2
				    
		    END DO
		END IF
! 				
        IF (dg.eq.1)then
        DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
        RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS
        


        end if


	END DO
	!$OMP END DO 
	IF (DG.EQ.1)THEN
	DEALLOCATE(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
	END IF


END SUBROUTINE CALCULATE_FLUXESHI
	
	
	
	
	
	
SUBROUTINE CALCULATE_FLUXESHI2D(N)
!> @brief
!> This subroutine computes the fluxes for linear-advection equation in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL::GODFLUX2,sum_detect,lamxl,lamyl
	INTEGER::I,L,K,NGP,KMAXE,IQP, NEIGHBOR_INDEX, NEIGHBOR_FACE_INDEX
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_TEMP,WEIGHTS_DG !Quadrature weights for interfaces
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
	 REAL::ANGLE1,ANGLE2,NX,NY,NZ,NORMALVECT
    INTEGER::facex,POINTX,ICONSIDERED
    REAL,DIMENSION(1:NOF_VARIABLES)::CLEFT,CRIGHT,HLLCFLUX,RHLLCFLUX
    REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG


    IF (DG.EQ.1)THEN
	allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
	END IF
	
	
	KMAXE = XMPIELRANK(N)
	
	CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_TEMP(1:QP_LINE_N) = WEQUA2D(1:QP_LINE_N)
	
	!$OMP DO
	DO I=1,KMAXE
	
        if (initcond.eq.3)then
            lamxl=-ielem(n,i)%yyc+0.5d0
            lamyl=ielem(n,i)%xxc-0.5
        else
			lamxl=lamx
			lamyl=lamy
		end if
        
        
        if (dg.eq.1)then
        RHS(I)%VALDG = ZERO
        DG_RHS = ZERO
        DG_RHS_SURF_INTEG = ZERO
        DG_RHS_VOL_INTEG = ZERO
        else
        RHS(I)%VAL=ZERO
        end if
        
        
        ICONSIDERED=I
        
        IF (DG.EQ.1) THEN
        
            DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
         
        END IF
        
		IF (IELEM(N,I)%INTERIOR.EQ.0)THEN ! Element is interior

		    DO L=1,IELEM(N,I)%IFCA
                GODFLUX2=ZERO
                NX=IELEM(N,I)%FACEANGLEX(L)
                NY=IELEM(N,I)%FACEANGLEY(L)
                facex=l
                
                NORMALVECT=(NX*LAMXl)+(NY*LAMYl)

                IQP=QP_LINE_N
                
                NEIGHBOR_INDEX = IELEM(N,I)%INEIGH(L)
                NEIGHBOR_FACE_INDEX = IELEM(N,I)%INEIGHN(L)
                
                DO NGP=1,IQP
                    POINTX=NGP


                    IF (DG.EQ.1) THEN
                        CLEFT = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                        CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                    ELSE !FV
                        CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
                        CRIGHT(1)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                    END IF
                    
                    CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)




                    IF (DG.EQ.1) THEN
                        ! Riemann flux at interface quadrature points times basis
                        RHLLCFLUX(1)=HLLCFLUX(1)
                         
                        DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)



                    ELSE !FV
                        GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                    END IF
                END DO
                
                IF (DG /= 1) RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2





		    END DO



		
        ELSE IF (IELEM(N,I)%INTERIOR.EQ.1)THEN

            DO L=1,IELEM(N,I)%IFCA
                FACEX = L
                NX=IELEM(N,I)%FACEANGLEX(L)
                NY=IELEM(N,I)%FACEANGLEY(L)
                NORMALVECT=(NX*LAMXl)+(NY*LAMYl)
                IQP=QP_LINE_N
                
                GODFLUX2=ZERO
                DO NGP=1,IQP
                    POINTX = NGP

                    IF (DG == 1) THEN
                        CLEFT = ILOCAL_RECON3(I)%ULEFT_DG(1:NOF_VARIABLES, L, NGP)
                    ELSE
                        CLEFT(1)=ILOCAL_RECON3(I)%ULEFT(1,L,NGP)
                    END IF
                    
                    IF (IELEM(N,I)%INEIGHB(L).EQ.N)THEN	!MY CPU ONLY
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN MY CPU
                                IF (DG == 1) THEN
                                    CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
                                ELSE !FV
                                    CRIGHT(1) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                                END IF
                            ELSE !NOT PERIODIC ONES IN MY CPU
                                CRIGHT(1:nof_variables)=CLEFT(1:nof_variables)
                            END IF
!                            
                        ELSE
                            IF (DG == 1) THEN
                                CRIGHT = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT_DG(1:NOF_VARIABLES, IELEM(N,I)%INEIGHN(L), NGP)
!                                 
                            ELSE !FV
                                CRIGHT(1) = ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1,IELEM(N,I)%INEIGHN(L),NGP)
                            END IF
                        END IF
                    ELSE !IN OTHER CPUS THEY CAN ONLY BE PERIODIC OR MPI NEIGHBOURS
                        IF (IELEM(N,I)%IBOUNDS(L).GT.0)THEN	!CHECK FOR BOUNDARIES
                            if (ibound(n,ielem(n,i)%ibounds(L))%icode.eq.5)then	!PERIODIC IN OTHER CPU
                            
                                IF (DG == 1) THEN
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                else
                                    CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                                end if
                            END IF
                        ELSE
                            IF (DG == 1) THEN
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_DG(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                            ELSE
                                CRIGHT(1:nof_variables)=IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1:nof_variables)
                            end if
                        END IF
                        
                    END IF
                    
                    CALL EXACT_RIEMANN_SOLVER(N,CLEFT,CRIGHT,NORMALVECT,HLLCFLUX)



                    
                    IF (DG.EQ.1) THEN
                        !Riemann flux at interface quadrature points times basis
                        RHLLCFLUX(1)=HLLCFLUX(1)
                         DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
                         




                    ELSE !FV
                        GODFLUX2=GODFLUX2+(HLLCFLUX(1)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
                    END IF
  
                END DO
                
                IF (DG /= 1) RHS(I)%VAL(1)=RHS(I)%VAL(1)+GODFLUX2
		    END DO



		END IF
		
		
		IF (DG == 1) DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG



		
        IF (DG == 1) RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS
        

		
	END DO
	!$OMP END DO 

	IF (DG.EQ.1)THEN
	DEALLOCATE(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
	END IF



END SUBROUTINE CALCULATE_FLUXESHI2D

SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE(N)
!> @brief
!> This subroutine computes the convective fluxes for hyperbolic conservation laws
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2, DG_VOL_REC,RHLLCFLUX,HLLCFLUX
	INTEGER::I,L,NGP,KMAXE,IQP,ii,IKAS,igoflux, icaseb,jx,jx2,B_CODE,SRF
	REAL::sum_detect,NORMS,TEMPXX
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP,WEIGHTS_DG
	INTEGER::ICONSIDERED,FACEX,POINTX
	REAL::ANGLE1,ANGLE2,NX,NY,NZ,MP_SOURCE1,MP_SOURCE2,MP_SOURCE3
	REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV,SRF_SPEEDROT
	REAL,DIMENSION(1:TURBULENCEEQUATIONS+PASSIVESCALAR)::CTURBL,CTURBR
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:NOF_VARIABLES)::SRF_SPEED
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
    real::x1,y1,z1
    REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG


    IF (DG.EQ.1)THEN
	allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
	END IF



	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)

	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)

	call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)
	

	
	
	
        do i=1,xmpielrank(n)
            RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
          
        end do
        !$OMP BARRIER
	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
            MP_SOURCE3=ZERO        

        IF(MRF.EQ.1)THEN
            SRF=ILOCAL_RECON3(I)%MRF
        END IF
		     IF (DG.EQ.1) THEN
		     RHS(I)%VALDG = ZERO
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
            
             DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)

            END IF
            
		    
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
				B_CODE=0

				  GODFLUX2=ZERO
				  MP_SOURCE2=ZERO
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
				  
				  facex=l
				  if (dg.eq.1)WEIGHTS_dg(1:iqp)=WEIGHTS_TEMP(1:IQP)
				  
				  do NGP=1,iqp	!for all the gaussian quadrature points
				  pointx=ngp


		CALL GET_STATES_INTERIOR(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT)

				  

			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  END IF
						  
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc



				      
				      CALL HLLC_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)



				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if

				      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)






				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				     				      
				      CASE(4)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL rROE_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				       
				       
				       
				       CASE(5)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				       END SELECT
				       
				      if (dg.eq.1)then
				      
				      DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
				      
				      
				      
				      else 
				       
				       
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				      end if
				      
				      
				      
				      
				      
				      
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				     
				     
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					  
					  
					  
					  NORMS=0.5*(CLEFT_ROT(2)+CRIGHT_ROT(2))


					IF (ILOCAL_RECON3(Iconsidered)%MRF.EQ.1)THEN
					  NORMS=NORMS-SRF_SPEEDROT(2)
					END IF
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
					  

					  END IF
					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
					  
					 
					  
					  
				      END IF
				      
				      
				  END DO
				    
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF
				    

				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

				    end if

		    END DO
                 IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(8)=RHS(I)%VAL(8)-(U_C(I)%VAL(1,8)*MP_SOURCE3)
                 
                 END IF
                 
                 IF (DG.eq.1)then

                   DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
                    RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS
                    
                    IF (MULTISPECIES.EQ.1)THEN


						RHS(I)%VALdg(1,8)=RHS(I)%VALdg(1,8)-(U_C(I)%VAL(1,8)*MP_SOURCE3)

					end if

                   




                 end if
                 
	END DO
	!$OMP END DO
	
	 !$OMP BARRIER
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
		 MP_SOURCE3=ZERO		
				
	IF(MRF.EQ.1)THEN
        SRF=ILOCAL_RECON3(I)%MRF
    END IF	   
		    IF (DG.EQ.1) THEN
		    rhs(i)%valdg=zero
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
            
             DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
                
            END IF

		    DO L=1,IELEM(N,I)%IFCA
		    FACEX=L




		    
				      
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
				  
				  if (dg.eq.1)WEIGHTS_dg(1:iqp)=WEIGHTS_TEMP(1:IQP)
				  
				  GODFLUX2=ZERO
				  
				   MP_SOURCE2=ZERO
								  
				  
				  do NGP=1,iqp
				  POINTX = NGP
				  
				      B_CODE=0




				      CALL GET_STATES_BOUNDS(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT)
				      

				      
				      
			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						
						  
						  
								    
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  END IF
						  
				      
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
								      
								      cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
								      cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
								    END IF
				      
				      
				      
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      



				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if

				      CASE(9)			!hllc

				      CALL HLL_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if




				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)




				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				      
				       CASE(4)		!roe
				       
				       
				       
				       IF (B_CODE.gt.0)THEN
				       
				       CALL RUSANOV_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       
				       
				       ELSE
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				     
				      
				      
				      CALL rROE_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      
				       RHLLCFLUX=HLLCFLUX
				       
				       END IF
				     			
                                    CASE(5)			!roe
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      
				     IF (B_CODE.LE.0)THEN
				      
				      CALL tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      ELSE
				      
				      
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      END IF
				      
				       RHLLCFLUX=HLLCFLUX
				     
				      
				       END SELECT
				       
				     if (dg.eq.1)then
				      DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
				      
				       
				      
				      else
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				      END IF
				      
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				       
				      
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					   
					  
					   
                     NORMS=0.5*(CLEFT_ROT(2)+CRIGHT_ROT(2))

				  IF (ILOCAL_RECON3(Iconsidered)%MRF.EQ.1)THEN
                            NORMS=NORMS-SRF_SPEEDROT(2)
					END IF
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
					  
					  
					  END IF
					  
					  if ((b_code.eq.3).or.(b_code.eq.4))then
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=zero
					  
					  end if
					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
					  
					  
					  
					  
					  
					  
				      END IF
				      
				      
				  END DO
				   
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				    IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF

				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

				    end if
! 				    end if
		    END DO
		     IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(8)=RHS(I)%VAL(8)-(U_C(I)%VAL(1,8)*MP_SOURCE3)
                 
                 END IF
                 
                    IF (DG.eq.1)then
                    DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
                    RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS
                    
                    IF (MULTISPECIES.EQ.1)THEN

						RHS(I)%VALdg(1,8)=RHS(I)%VALdg(1,8)-(U_C(I)%VAL(1,8)*MP_SOURCE3)


					END IF



                    end if
                
	END DO
	!$OMP END DO

	IF (DG.EQ.1)THEN
	DEALLOCATE(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
	END IF






END SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE
	
	
SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE2d(N)
!> @brief
!> This subroutine computes the convective fluxes for hyperbolic conservation laws in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2, DG_VOL_REC,RHLLCFLUX,HLLCFLUX
	INTEGER::I,L,NGP,KMAXE,IQP,ii,IKAS,igoflux, icaseb,kxk,B_CODE
	REAL::sum_detect,NORMS
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_TEMP,WEIGHTS_DG
	INTEGER::ICONSIDERED,FACEX,POINTX
	REAL::ANGLE1,ANGLE2,NX,NY,NZ,MP_SOURCE1,MP_SOURCE2,MP_SOURCE3
	REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV,SRF_SPEEDROT
	REAL,DIMENSION(1:TURBULENCEEQUATIONS+PASSIVESCALAR)::CTURBL,CTURBR
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	real::x1,y1,z1
	REAL,DIMENSION(1:NOF_VARIABLES)::SRF_SPEED
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
	REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG


    IF (DG.EQ.1)THEN
	allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
	END IF
	KMAXE=XMPIELRANK(N)
	
	CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_TEMP = WEQUA2D(1:QP_LINE_N)
	
	if(reduce_comp.eq.1)then
	WEIGHTS_TEMP=1.0d0
	end if
	
	do i=1,kmaxe
	RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
	
    
	
	ICONSIDERED=I
        
        
	
	
	end do




		
	!$OMP BARRIER
	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
            IF (DG.EQ.1) THEN
            RHS(I)%VALDG = ZERO
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
            
            
            DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
            
            
            END IF
	

                    
                MP_SOURCE3=ZERO     
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
		    


				  GODFLUX2=ZERO
				  MP_SOURCE2=ZERO
 				  nx=IELEM(N,I)%FACEANGLEX(L)
 				  NY=IELEM(N,I)%FACEANGLEY(L)
 				  angle1=nx
 				  angle2=ny
				 b_code=0
					iqp=qp_line_n
				
				  do NGP=1,iqp	!for all the gaussian quadrature points
				  POINTX=NGP
				   FACEX=L
		    POINTX=NGP



				  CALL GET_STATES_INTERIOR2D(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT)



			
						  CALL ROTATEF2d(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  END IF
						  
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if


				      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if



				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
                      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY)
				      
				      RHLLCFLUX(1:nof_Variables)=HLLCFLUX(1:nof_Variables)
                                        
                                        
                                         CASE(4)			!roe
				      
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL RROE_RIEMANN_SOLVER2d(N,Cleft,Cright,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,b_code)
				      
				      RHLLCFLUX=HLLCFLUX
				     
				      
				       END SELECT
				       
				       
				       
				       IF (DG.EQ.1) THEN
                        ! Riemann flux at interface quadrature points times basis
                         

                         
                        DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
                       
                         

                    ELSE !FV
				       
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      
				      END IF
				      
				       IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					    
                     NORMS=0.5*(CLEFT_ROT(2)+CRIGHT_ROT(2))
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
					  
! 					  
					  END IF
					  
 					
					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
					  
! 					 
					  
				      END IF
				      
				      
				  END DO
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				    IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF
				    
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

! 				    
				    END IF
! 				    end if
		    END DO
		    IF (MULTISPECIES.EQ.1)THEN

                 RHS(I)%VAL(7)=RHS(I)%VAL(7)-(U_C(I)%VAL(1,7)*MP_SOURCE3)!*ielem(n,I)%totvolume)

                 END IF
                 
               IF (DG == 1) THEN

               
                  
               
               DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
                RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS  
                
                IF (MULTISPECIES.EQ.1)THEN

					RHS(I)%VALdg(1,7)=RHS(I)%VALdg(1,7)-(U_C(I)%VAL(1,7)*MP_SOURCE3)


                 
                 
                 
                 END IF
                
              
                 END IF
                 
                 
	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
	MP_SOURCE3=ZERO  
	
            IF (DG.EQ.1) THEN
            rhs(i)%valdg=zero
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
            
            DG_RHS_VOL_INTEG = DG_VOL_INTEGRAL(N,ICONSIDERED)
               
            END IF
				

		    
		    DO L=1,IELEM(N,I)%IFCA
		    FACEX=L
				  igoflux=0


				    b_code=0  
				 nx=IELEM(N,I)%FACEANGLEX(L)
 				  NY=IELEM(N,I)%FACEANGLEY(L)
 				  angle1=nx
 				  angle2=ny
				 
					iqp=qp_line_n
				  GODFLUX2=ZERO
				  MP_SOURCE2=ZERO
				  do NGP=1,iqp
				  POINTX=NGP
				  FACEX=L

				  CALL GET_STATES_BOUNDS2D(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT)
				  
				  
				  

				      
				      
			
						  CALL ROTATEF2d(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						   
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  
						  
						  END IF
						  
				      
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
								      
								      cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
								      cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
								    END IF
				      
				      
				      
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if






				      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY)
				      
				      RHLLCFLUX=HLLCFLUX
				     				      
				     
				        CASE(4)			!roe
				      
				      
				       
				      
				      IF ((B_CODE.le.0))THEN
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      
				      
				      
				      CALL RROE_RIEMANN_SOLVER2d(N,Cleft,Cright,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,b_code)
				      RHLLCFLUX=HLLCFLUX
				      ELSE
				      

				     
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				      
				      
				      END IF
				      
				      
				       END SELECT
				       
				      IF (DG.EQ.1) THEN
                        
                         DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
                         
                         ELSE
				      
                         GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				       
				      
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  if (icoupleturb.eq.0)then	!first order upwind flux
					    
					  NORMS=0.5*(CLEFT_ROT(2)+CRIGHT_ROT(2))
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
					  
! 					  
					  END IF
					  
					  if ((b_code.eq.4).or.(b_code.eq.3))then
					  rHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=zero
					  end if
					  
					  
					 
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
! 					  
					  
					  
				      END IF
				      
				      
				  END DO
				   
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				    IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF

				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

				    end if
		    END DO
		     IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(7)=RHS(I)%VAL(7)-(U_C(I)%VAL(1,7)*MP_SOURCE3)!*ielem(n,I)%totvolume)
                 
                 END IF
                 
                  IF (DG == 1) THEN
                
               DG_RHS = DG_RHS_SURF_INTEG - DG_RHS_VOL_INTEG
                RHS(I)%VALDG = RHS(I)%VALDG + DG_RHS  
               
                IF (MULTISPECIES.EQ.1)THEN
                  RHS(I)%VALdg(1,7)=RHS(I)%VALdg(1,7)-(U_C(I)%VAL(1,7)*MP_SOURCE3)

                 END IF
                 
                 END IF
                 
                 
                 
                 
	END DO
	!$OMP END DO
	IF (DG.EQ.1)THEN
	DEALLOCATE(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
	END IF



END SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE2d

	
	

SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE(N)
!> @brief
!> This subroutine computes the diffusive fluxes for Euler-Navier-Stokes Equations
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2,DG_VOL_REC,RHLLCFLUX,HLLCFLUX
	INTEGER::I,L,NGP,KMAXE,IQP,ii,NVAR,KC,IEX,ITTT,IKAS,igoflux, icaseb,KK,B_CODE,SRF
	REAL::sum_detect,NORMS
	INTEGER::ICONSIDERED,FACEX,POINTX
	REAL::ANGLE1,ANGLE2,NX,NY,NZ,MP_SOURCE1,MP_SOURCE2,MP_SOURCE3
	REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV,SRF_SPEEDROT
	REAL,DIMENSION(1:TURBULENCEEQUATIONS+PASSIVESCALAR)::CTURBL,CTURBR
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:NOF_VARIABLES)::SRF_SPEED
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP,WEIGHTS_DG
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
	REAl::MP_PINFL,MP_PINFR,GAMMAL,GAMMAR
	REAL,DIMENSION(1:4)::VISCL,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
    REAL,DIMENSION(1:nof_Variables-1,1:dims)::LCVGRAD,RCVGRAD
	REAL,DIMENSION(turbulenceequations+passivescalar,1:dims)::LCVGRAD_T,RCVGRAD_T
	real,dimension(5)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(3,3)::taul,taur,TAU
	REAL,DIMENSION(3)::Q,NNN,nall
	REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,RHO12,U12,V12,W12 ,damp,vdamp,TEMPXX 
	REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG


    IF (DG.EQ.1)THEN
	allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
	END IF


	KMAXE=XMPIELRANK(N)
	
	call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)

	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)

	call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)

	if(reduce_comp.eq.1)then
	WEIGHTS_T=1.0d0;WEIGHTS_Q=1.0d0
	end if
	
	
	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I

		IF (DG.EQ.1) THEN
            DG_RHS = ZERO
            DG_RHS_SURF_INTEG = ZERO
            DG_RHS_VOL_INTEG = ZERO
            END IF


                    
		   damp=lamx
            IF( BR2_YN == 2) DAMP = 0.0D0
		   IF(MRF.EQ.1)THEN
                SRF=ILOCAL_RECON3(I)%MRF
            END IF
		    DO L=1,IELEM(N,I)%IFCA !for all their faces

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
				  
				  IF( DG == 1) WEIGHTS_DG = WEIGHTS_TEMP
				  
				  do NGP=1,iqp	!for all the gaussian quadrature points

					FACEX=L
					POINTX=NGP

					CALL CALCULATE_INTERIOR_VISCOUS(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT,LCVGRAD,RCVGRAD,LCVGRAD_T,RCVGRAD_T)





			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
					CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:6)= LCVGRAD(1,1:3);EDDYFL(7:9)=LCVGRAD(2,1:3)
							  EDDYFL(10:12)=LCVGRAD(3,1:3);EDDYFL(13:15)=LCVGRAD_T(1,1:3)
							  EDDYFL(16:18)=LCVGRAD_T(2,1:3)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:6)= RCVGRAD(1,1:3);EDDYFR(7:9)=RCVGRAD(2,1:3);EDDYFR(10:12)=RCVGRAD(3,1:3)
							  EDDYFR(13:15)=RCVGRAD_T(1,1:3);EDDYFL(16:18)=RCVGRAD_T(2,1:3)
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					  END IF
				       


						TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
					  U12=ZERO;V12=ZERO;W12=ZERO  
					  

					  
					  vdamp=(4.0/3.0)!*(( (VISCL(1))+(VISCL(2)))))
                                        nall(1)=nx;nall(2)=ny;nall(3)=nz
                                      LCVGRAD(1,1:3)=((LCVGRAD(1,1:3)+rCVGRAD(1,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(2)-leftv(2)))
					   LCVGRAD(2,1:3)=((LCVGRAD(2,1:3)+rCVGRAD(2,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(3)-leftv(3)))
					  LCVGRAD(3,1:3)=((LCVGRAD(3,1:3)+rCVGRAD(3,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*(rightv(4)-leftv(4)))
                                        LCVGRAD(4,1:3)=((LCVGRAD(4,1:3)+rCVGRAD(4,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*((rightv(5)/(rightv(1)*R_gas))-(leftv(5)/(leftv(1)*R_gas))))
				       				  
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
		
		
					  HLLCFLUX(1:nof_Variables)=(NX*FXV+NY*FYV+NZ*FZV)	


					  if (dg.eq.1)then
					  RHLLCFLUX(1:nof_Variables)=HLLCFLUX(1:nof_Variables)

					  DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)

					  else
					  
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(HLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))

				      end if
				      
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					  IF (TURBULENCE.EQ.1)THEN
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2)))+(OO2*(VISCL(3)+VISCL(4))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  ELSE
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  
					  END IF
						IF (TURBULENCEMODEL.EQ.1)THEN
						HLLCFLUX(6)=HLLCFLUX(6)/SIGMA
						END IF					  
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)

				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

 				    end if

		    END DO

		    IF (DG.eq.1)then




                    DG_RHS = DG_RHS_SURF_INTEG

                    RHS(I)%VALDG = RHS(I)%VALDG - DG_RHS




                 end if





	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
			IF (DG.EQ.1) THEN
				DG_RHS = ZERO
				DG_RHS_SURF_INTEG = ZERO
				DG_RHS_VOL_INTEG = ZERO



				END IF
				
		     IF(MRF.EQ.1)THEN
                SRF=ILOCAL_RECON3(I)%MRF
            END IF
		     
		   
		    DO L=1,IELEM(N,I)%IFCA



                                     igoflux=0




				   damp=lamx  
                    IF( BR2_YN == 2) DAMP = 0.0D0
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
				  
				  IF( DG == 1) WEIGHTS_DG = WEIGHTS_TEMP
				  
				  GODFLUX2=ZERO
				  
				  b_code=0
								  
				  
				  do NGP=1,iqp
				  FACEX=l
				  pointx=ngp

				  CALL CALCULATE_BOUNDED_VISCOUS(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT,LCVGRAD,RCVGRAD,LCVGRAD_T,RCVGRAD_T)




				      
				    
			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
					CALL SUTHERLAND(N,LEFTV,RIGHTV,VISCL,LAML)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:6)= LCVGRAD(1,1:3);EDDYFL(7:9)=LCVGRAD(2,1:3)
							  EDDYFL(10:12)=LCVGRAD(3,1:3);EDDYFL(13:15)=LCVGRAD_T(1,1:3)
							  EDDYFL(16:18)=LCVGRAD_T(2,1:3)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:6)= RCVGRAD(1,1:3);EDDYFR(7:9)=RCVGRAD(2,1:3);EDDYFR(10:12)=RCVGRAD(3,1:3)
							  EDDYFR(13:15)=RCVGRAD_T(1,1:3);EDDYFL(16:18)=RCVGRAD_T(2,1:3)
							    Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					  END IF


                                            
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
                                        LCVGRAD(4,1:3)=((LCVGRAD(4,1:3)+rCVGRAD(4,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:3)*((rightv(5)/(rightv(1)*R_gas))-(leftv(5)/(leftv(1)*R_gas))))
				       				  
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
		
		
					  HLLCFLUX(1:nof_Variables)=(NX*FXV+NY*FYV+NZ*FZV)	


					   if (dg.eq.1)then
					   RHLLCFLUX(1:nof_Variables)=HLLCFLUX(1:nof_Variables)

					  DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)


					  else
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(HLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      end if
				      
				      
				     
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
					   IF (TURBULENCE.EQ.1)THEN
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2)))+(OO2*(VISCL(3)+VISCL(4))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  ELSE
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR) =&
					  ((OO2*(VISCL(1)+VISCL(2))))*&
		    (((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,1))*OO2*NX)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,2))*OO2*NY)+&
		    ((LCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3)+RCVGRAD_T(1:TURBULENCEEQUATIONS+PASSIVESCALAR,3))*OO2*NZ))
					  
					  END IF
						IF (TURBULENCEMODEL.EQ.1)THEN
						HLLCFLUX(6)=HLLCFLUX(6)/SIGMA
						END IF					  

					if ((b_code.eq.3))then
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=zero
					  
					  end if




					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)

				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

				    end if
		    END DO

		     IF (DG.eq.1)then



                    DG_RHS = DG_RHS_SURF_INTEG

                    RHS(I)%VALDG = RHS(I)%VALDG - DG_RHS






                 end if



	END DO
	!$OMP END DO


	IF (DG.EQ.1)THEN
	DEALLOCATE(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
	END IF


END SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE




SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE2d(N)
!> @brief
!> This subroutine computes the diffusive fluxes for Euler-Navier-Stokes Equations in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2,DG_VOL_REC,RHLLCFLUX,HLLCFLUX
	INTEGER::I,L,NGP,KMAXE,IQP,ii,NVAR,KC,IEX,ITTT,IKAS,igoflux, icaseb,KK,B_CODE
	REAL::sum_detect,NORMS
	INTEGER::ICONSIDERED,FACEX,POINTX
	REAL::ANGLE1,ANGLE2,NX,NY,NZ,MP_SOURCE1,MP_SOURCE2,MP_SOURCE3
	REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV,SRF_SPEEDROT
	REAL,DIMENSION(1:TURBULENCEEQUATIONS+PASSIVESCALAR)::CTURBL,CTURBR
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:NOF_VARIABLES)::SRF_SPEED
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP,WEIGHTS_DG
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
	REAl::MP_PINFL,MP_PINFR,GAMMAL,GAMMAR
	REAL,DIMENSION(1:4)::VISCL,LAML
	REAL,DIMENSION(1:2)::TURBMV
    REAL,DIMENSION(1)::ETVM
    REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
    REAL,DIMENSION(1:nof_Variables-1,1:dims)::LCVGRAD,RCVGRAD
	REAL,DIMENSION(turbulenceequations+passivescalar,1:dims)::LCVGRAD_T,RCVGRAD_T
	real,dimension(4)::fxv,fyv,fzv,tem_pn,rtem_pn
	real,dimension(2,2)::taul,taur,TAU
	REAL,DIMENSION(2)::Q,NALL
	REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,RHO12,U12,V12,W12,damp,vdamp  
	REAL,allocatable,dimension(:,:)::DG_RHS, DG_RHS_VOL_INTEG, DG_RHS_SURF_INTEG


    IF (DG.EQ.1)THEN
	allocate(DG_RHS(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_VOL_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES), DG_RHS_SURF_INTEG(1:NUM_DG_DOFS,1:NOF_VARIABLES))
	END IF
	
	
	
	KMAXE=XMPIELRANK(N)
	
	 CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	
	WEIGHTS_temp(1:QP_line_n)=WEQUA2D(1:QP_line_n)
	if (Reduce_comp.eq.1)then
	WEIGHTS_temp=1.0d0
	end if
	
	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
		   
			IF (DG.EQ.1) THEN
				DG_RHS = ZERO
				DG_RHS_SURF_INTEG = ZERO
				DG_RHS_VOL_INTEG = ZERO



				END IF



		    DO L=1,IELEM(N,I)%IFCA !for all their faces
!                 IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
                    damp=LAMX
                    IF( BR2_YN == 2) DAMP = 0.0D0
                    
				  GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=angle1
				  NY=angle2
				 
				  
					iqp=qp_line_n
				
				  
				  do NGP=1,iqp	!for all the gaussian quadrature points
				  facex=L
				  POINTX=NGP

				  CALL CALCULATE_INTERIOR_VISCOUS2D(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT,LCVGRAD,RCVGRAD,LCVGRAD_T,RCVGRAD_T)




			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF2d(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,leftv,rightv)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
					CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
						      IF (TURBULENCEMODEL.EQ.2)THEN
							  EDDYFL(1)=IELEM(N,I)%WALLDIST;EDDYFL(2)=CTURBL(1);EDDYFL(3)=CTURBL(2)
							  EDDYFL(4:5)= LCVGRAD(1,1:2);EDDYFL(6:7)=LCVGRAD(2,1:2)
							  EDDYFL(8:9)=LCVGRAD_T(1,1:2)
							  EDDYFL(10:11)=LCVGRAD_T(2,1:2)
							    
							    
							  EDDYFR(1)=IELEM(N,I)%WALLDIST;EDDYFR(2)=CTURBR(1);EDDYFR(3)=CTURBR(2)
							  EDDYFR(4:5)= RCVGRAD(1,1:2);EDDYFR(6:7)=RCVGRAD(2,1:2)
							  EDDYFR(8:9)=RCVGRAD_T(1,1:2);EDDYFL(10:11)=RCVGRAD_T(2,1:2)
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					  END IF
				       
				       



                                     TAUL = ZERO;TAU=ZERO;TAUR=ZERO;Q=ZERO;UX=ZERO;UY=ZERO;UZ=ZERO;VX=ZERO;VY=ZERO;VZ=ZERO;WX=ZERO;WY=ZERO;WZ=ZERO;
					    FXV=ZERO;FYV=ZERO;FZV=ZERO;RHO12 =ZERO;
					  U12=ZERO;V12=ZERO;W12=ZERO 
				       
! 				      
					  
					  vdamp=(4.0/3.0)!*(( (VISCL(1))+(VISCL(2)))))
                                        nall(1)=nx;nall(2)=ny
					  				
                                      LCVGRAD(1,1:2)=((LCVGRAD(1,1:2)+rCVGRAD(1,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*(rightv(2)-leftv(2)))
					   LCVGRAD(2,1:2)=((LCVGRAD(2,1:2)+rCVGRAD(2,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*(rightv(3)-leftv(3)))
					  LCVGRAD(3,1:2)=((LCVGRAD(3,1:2)+rCVGRAD(3,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*((rightv(4)/(rightv(1)*R_gas))-(leftv(4)/(leftv(1)*R_gas))))
                                        
				       
					  
					  
					  
					  
					   
				       				  
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
					  HLLCFLUX(1:nof_Variables)=(NX*FXV+NY*FYV)	
		
					  if (dg.eq.1)then

					  RHLLCFLUX(1:nof_Variables)=HLLCFLUX(1:nof_Variables)
						
						DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)
  
						else
					  
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(HLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
						end if
				     
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
				  
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)

				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

 				    end if
		    END DO

			IF (DG.eq.1)then



				DG_RHS = DG_RHS_SURF_INTEG

				

				RHS(I)%VALDG = RHS(I)%VALDG - DG_RHS






			 end if



	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
			IF (DG.EQ.1) THEN
				DG_RHS = ZERO
				DG_RHS_SURF_INTEG = ZERO
				DG_RHS_VOL_INTEG = ZERO



				END IF
		    
		    
		    DO L=1,IELEM(N,I)%IFCA
                damp=lamx
                IF( BR2_YN == 2) DAMP = 0.0D0
				      b_code=0
				 GODFLUX2=ZERO
 				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
 				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
 				  NX=angle1
				  NY=angle2
				 
				  
					iqp=qp_line_n
				
				  
				 
								  
				  
				  do NGP=1,iqp



				   FACEX=L
				   POINTX=NGP

				   CALL CALCULATE_BOUNDED_VISCOUS2D(N,B_CODE,ICONSIDERED,FACEX,POINTX,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEEDROT,CLEFT,CRIGHT,LCVGRAD,RCVGRAD,LCVGRAD_T,RCVGRAD_T)




				      
				    
			
					
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  CALL ROTATEF2d(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,leftv,rightv)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
						  CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
						  END IF
						  
				      
				        LEFTV(1:NOF_vARIABLES)=CLEFT(1:NOF_vARIABLES);RIGHTV(1:NOF_vARIABLES)=CRIGHT(1:NOF_vARIABLES)
					CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
					CALL SUTHERLAND2D(N,LEFTV,RIGHTV,VISCL,LAML)
				     
					    IF (TURBULENCE.EQ.1)THEN
						      IF (TURBULENCEMODEL.EQ.1)THEN
							  TURBMV(1)=CTURBL(1);  TURBMV(2)=CTURBR(1);eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
							  Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
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
							    Call EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF
					  END IF
				       
				       
				       
				       



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
					  LCVGRAD(3,1:2)=((LCVGRAD(3,1:2)+rCVGRAD(3,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(L)))*nall(1:2)*((rightv(4)/(rightv(1)*R_gas))-(leftv(4)/(leftv(1)*R_gas))))
					  
					  
					  
					  
					   
				       				  
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
					  HLLCFLUX(1:nof_Variables)=(NX*FXV+NY*FYV)			
					  
					  if (dg.eq.1)then

					   RHLLCFLUX(1:nof_Variables)=HLLCFLUX(1:nof_Variables)

						DG_RHS_SURF_INTEG = DG_RHS_SURF_INTEG + DG_SURF_FLUX(N,ICONSIDERED,FACEX,POINTX,WEIGHTS_TEMP,RHLLCFLUX)





  
						else
		



			      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(HLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
						end if
				     
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
						
						if ((b_code.eq.3))then
					  HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=zero
					  
					  end if
						
						
						
						
					  GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)+&
					  (HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      END IF
				      
				      
				  END DO
				  
				    
				     RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)
				    
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)

				    end if
		    END DO

			IF (DG.eq.1)then



				DG_RHS = DG_RHS_SURF_INTEG


				RHS(I)%VALDG = RHS(I)%VALDG - DG_RHS






			 end if



	END DO
	!$OMP END DO
	IF (DG.EQ.1)THEN
	DEALLOCATE(DG_RHS,DG_RHS_VOL_INTEG,DG_RHS_SURF_INTEG)
	END IF



END SUBROUTINE CALCULATE_FLUXESHI_DIFFUSIVE2d




SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE_MOOD(N)
!> @brief
!> This subroutine computes the convective fluxes for hyperbolic conservation laws
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2, DG_VOL_REC,RHLLCFLUX,HLLCFLUX
	INTEGER::I,L,NGP,KMAXE,IQP,ii,IKAS,igoflux, icaseb,jx,jx2,B_CODE
	REAL::sum_detect,NORMS,TEMPXX
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_Q,WEIGHTS_T,WEIGHTS_TEMP,WEIGHTS_DG
	INTEGER::ICONSIDERED,FACEX,POINTX,N_NODE
	REAL::ANGLE1,ANGLE2,NX,NY,NZ,MP_SOURCE1,MP_SOURCE2,MP_SOURCE3
	REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV,SRF_SPEEDROT
	REAL,DIMENSION(1:TURBULENCEEQUATIONS+PASSIVESCALAR)::CTURBL,CTURBR
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:NOF_VARIABLES)::SRF_SPEED
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
    REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
    REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
	REAL,DIMENSION(1:8,1:DIMENSIONA)::NODES_LIST
	REAL,DIMENSION(1:DIMENSIONA)::CORDS
	INTEGER::IBFC
	KMAXE=XMPIELRANK(N)
	
	KMAXE=XMPIELRANK(N)

	call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)

	WEIGHTS_Q(1:QP_QUAD)=WEQUA2D(1:QP_QUAD)

	call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_T(1:QP_TRIANGLE)=WEQUA2D(1:QP_TRIANGLE)
	
	if(reduce_comp.eq.1)then
	WEIGHTS_T=1.0d0;WEIGHTS_Q=1.0d0
	end if
        do i=1,xmpielrank(n)
            IF (IELEM(N,I)%RECALC.EQ.1)THEN
            RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
            END IF
        end do
        !$OMP BARRIER
	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
            MP_SOURCE3=ZERO        
		    B_CODE=0
		    IF (IELEM(N,I)%RECALC.EQ.1)THEN 
		    
		    DO L=1,IELEM(N,I)%IFCA !for all their faces
!                                      IF (IELEM(N,I)%REORIENT(l).EQ.0)THEN
				  GODFLUX2=ZERO
				  MP_SOURCE2=ZERO
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
                    
				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      IF ((CASCADE.EQ.2).AND.(IELEM(N,I)%MOOD.EQ.1))THEN
				      CLEFT(1:nof_Variables)=U_c(I)%VAL(3,1:nof_variables)
				      END IF
				      
				      
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
				       IF ((CASCADE.EQ.2).AND.(IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD.EQ.1))THEN
				      CRIGHT(1:nof_Variables)=U_c(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
				      END IF
				      
				      !right mean flow state
				      
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
			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  END IF
						  
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if


				      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if


				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				     				      
				      CASE(4)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL rROE_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				       
				       
				       
				       CASE(5)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				       END SELECT
				       
				       
				       
				       
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				     
				     
				     
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
				    
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF
				    
!  				     RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)-godflux2(1:nof_Variables)
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))then
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
!  				    RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
! 				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				    end if
!  				    end if
		    END DO
                 IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(8)=RHS(I)%VAL(8)-(U_C(I)%VAL(1,8)*MP_SOURCE3)
                 
                 END IF
            END IF
	END DO
	!$OMP END DO
	
	 !$OMP BARRIER
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
		 MP_SOURCE3=ZERO		
		   IF (IELEM(N,I)%RECALC.EQ.1)THEN
		    
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

		    

		    
				      
				  ANGLE1=IELEM(N,I)%FACEANGLEX(L)
				  ANGLE2=IELEM(N,I)%FACEANGLEY(L)
				NX=(COS(ANGLE1)*SIN(ANGLE2))
				NY=(SIN(ANGLE1)*SIN(ANGLE2))
				NZ=(COS(ANGLE2))
 				  
 				  if (ielem(n,i)%types_faces(L).eq.5)then
					iqp=qp_quad_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_Q(1:IQP)
					N_NODE=4
				  else
					iqp=QP_TRIANGLE_n
					WEIGHTS_TEMP(1:IQP)=WEIGHTS_T(1:IQP)
					N_NODE=3
				  end if
				  GODFLUX2=ZERO
				  
				   MP_SOURCE2=ZERO
								  
				  
				  do NGP=1,iqp
				      B_CODE=0
				      
				      
				      
				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
				     IF ((CASCADE.EQ.2).AND.(IELEM(N,I)%MOOD.EQ.1))THEN
				     CLEFT(1:nof_Variables)=U_c(I)%VAL(3,1:nof_variables)
				     
				     END IF
				      
				      
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
								  
								  CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
								  
								 
								  
								  
								  
								  IF ((CASCADE.EQ.2).AND.(IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD.EQ.1))THEN
								   CRIGHT(1:nof_variables)=U_c(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
								  
								  END IF
								  
								  
								  
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
									if (icoupleturb.eq.1)then
									IF (CASCADE.EQ.1)THEN
									   CTURBR(1:turbulenceequations+PASSIVESCALAR)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFTTURB&
									  (1:turbulenceequations+PASSIVESCALAR,IELEM(N,I)%INEIGHN(L),ngp)!right additional equations flow state
									ELSE
									CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									
									END IF
									  
									  
									ELSE
									 CTURBR(1:turbulenceequations+PASSIVESCALAR)=U_CT(IELEM(N,I)%INEIGH(L))%VAL(1,1:turbulenceequations+PASSIVESCALAR)
									END IF
								    END IF
								  
								  
								  
								  
								  ELSE
								  !NOT PERIODIC ONES IN MY CPU
								   
								  facex=l;iconsidered=i
								  CALL  coordinates_face_innerx(N,ICONSIDERED,FACEX,VEXT,NODES_LIST)

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
							
							
                                  
							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
							      IF ((CASCADE.EQ.2).AND.(IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD.EQ.1))THEN
							      CRIGHT(1:nof_Variables)=U_c(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
							      
							      END IF
							      
							      
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
                                    IF ((CASCADE.EQ.2).AND.(IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_M(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1).GT.0.5))THEN
                                    CRIGHT(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                                    
                                    END IF
									   
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
								  
								  IF ((CASCADE.EQ.2).AND.(IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_M(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1).GT.0.5))THEN
								  CRIGHT(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
								  
								  END IF
								  
								  
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
				      
				      
			
						  CALL ROTATEF(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						
						  
						  
								    
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  END IF
						  
				      
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
								      
								      cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
								      cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
								    END IF
				      
				      
				      
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      RHLLCFLUX=HLLCFLUX
				      
				      
				      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				       CASE(4)			!roe
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				     
				      
				      
				      
				      
				     IF (B_CODE.LE.0)THEN
				      
				      CALL rROE_RIEMANN_SOLVER(N,iconsidered,facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      ELSE
				      
				      
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      END IF
				      
				       RHLLCFLUX=HLLCFLUX
				     			
                                    CASE(5)			!roe
				      
				      CALL ROTATEB(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      
				     IF (B_CODE.LE.0)THEN
				      
				      CALL tROE_RIEMANN_SOLVER(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      ELSE
				      
				      
				      CALL ROE_RIEMANN_SOLVER(N,iconsidered, facex,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,NZ)
				      
				      END IF
				      
				       RHLLCFLUX=HLLCFLUX
				     
				      
				       END SELECT
				       
				     
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				       
				      
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
				   
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				    IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF
! 				    if ((igoflux.eq.1))then
! 				    RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)
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
		     IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(8)=RHS(I)%VAL(8)-(U_C(I)%VAL(1,8)*MP_SOURCE3)
                 
                 END IF
                END IF
	END DO
	!$OMP END DO

END SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE_MOOD



SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE2d_MOOD(N)
!> @brief
!> This subroutine computes the convective fluxes for hyperbolic conservation laws in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	REAL,DIMENSION(1:NOF_variables+TURBULENCEEQUATIONS+PASSIVESCALAR)::GODFLUX2, DG_VOL_REC,RHLLCFLUX,HLLCFLUX
	INTEGER::I,L,NGP,KMAXE,IQP,ii,IKAS,igoflux, icaseb,kxk,B_CODE
	REAL::sum_detect,NORMS
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEIGHTS_TEMP,WEIGHTS_DG
	INTEGER::ICONSIDERED,FACEX,POINTX,n_node
	REAL::ANGLE1,ANGLE2,NX,NY,NZ,MP_SOURCE1,MP_SOURCE2,MP_SOURCE3
	REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV,RIGHTV,SRF_SPEEDROT
	REAL,DIMENSION(1:TURBULENCEEQUATIONS+PASSIVESCALAR)::CTURBL,CTURBR
	REAL,DIMENSION(1:DIMENSIONA)::POX,POY,POZ
	REAL,DIMENSION(1:NOF_VARIABLES)::SRF_SPEED
	REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT,NODES_LIST
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
	real,dimension(1:dimensiona)::cords
	INTEGER::IBFC
	KMAXE=XMPIELRANK(N)
	
	CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
	WEIGHTS_TEMP = WEQUA2D(1:QP_LINE_N)

	if(reduce_comp.eq.1)then
	WEIGHTS_TEMP=1.0d0
	end if
	
	do i=1,kmaxe
	IF (IELEM(N,I)%RECALC.EQ.1)THEN
	RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
	END IF
	end do
	
	!$OMP BARRIER
	!$OMP DO
	DO II=1,NOF_INTERIOR	!for all the interior elements
	I=EL_INT(II)
	ICONSIDERED=I
! 		    RHS(I)%VAL(:)=ZERO;IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0)) RHST(I)%VAL(:)=ZERO 
          IF (IELEM(N,I)%RECALC.EQ.1)THEN             
                MP_SOURCE3=ZERO     
		    DO L=1,IELEM(N,I)%IFCA !for all their faces

                                
                                
				  GODFLUX2=ZERO
				  MP_SOURCE2=ZERO
 				  nx=IELEM(N,I)%FACEANGLEX(L)
 				  NY=IELEM(N,I)%FACEANGLEY(L)
 				  angle1=nx
 				  angle2=ny
				 b_code=0
					iqp=qp_line_n
				
				  do NGP=1,iqp	!for all the gaussian quadrature points
				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)	!left mean flow state
				      
				       IF ((CASCADE.EQ.2).AND.(IELEM(N,I)%MOOD.EQ.1))THEN
				      CLEFT(1:nof_Variables)=U_c(I)%VAL(3,1:nof_variables)
				      END IF
				      
				      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP) !right mean flow state
				      
				      IF ((CASCADE.EQ.2).AND.(IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD.EQ.1))THEN
				      CRIGHT(1:nof_Variables)=U_c(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
				      END IF
				      
				      
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
			
						  CALL ROTATEF2d(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  CALL LMACHT2d(N,leftv,rightv)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  END IF
						  
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY)
				      
				      RHLLCFLUX(1:nof_Variables)=HLLCFLUX(1:nof_Variables)
                                        


                      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
                                        
                                         CASE(4)			!roe
				      
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL RROE_RIEMANN_SOLVER2d(N,Cleft,Cright,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,b_code)
				      
				      RHLLCFLUX=HLLCFLUX
				     
				      
				       END SELECT
				       
				       
				       
				       
				       
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				       IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
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
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				    IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF
				    
				    
! 				     RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)
				    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
				    RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(I)%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)+&
				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)=RHST(IELEM(N,I)%INEIGH(L))%VAL(1:TURBULENCEEQUATIONS+PASSIVESCALAR)-&
!  				    GODFLUX2(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
! 				    
				    END IF
! 				    end if
		    END DO
		    IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(7)=RHS(I)%VAL(7)-(U_C(I)%VAL(1,7)*MP_SOURCE3)!*ielem(n,I)%totvolume)
                 
                 END IF
                 
                 
            END IF     
	END DO
	!$OMP END DO
	
	
	!$OMP DO
	DO II=1,NOF_BOUNDED
	I=EL_BND(II)
	ICONSIDERED=I	
	MP_SOURCE3=ZERO  
		 IF (IELEM(N,I)%RECALC.EQ.1)THEN	

		    
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
				  MP_SOURCE2=ZERO
				  do NGP=1,iqp
				      CLEFT(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(1:nof_Variables,L,NGP)
				      
				      IF ((CASCADE.EQ.2).AND.(IELEM(N,I)%MOOD.EQ.1))THEN
				     CLEFT(1:nof_Variables)=U_c(I)%VAL(3,1:nof_variables)
				     
				     END IF
				      
				      
				      
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
								  CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
								  
								   IF ((CASCADE.EQ.2).AND.(IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD.EQ.1))THEN
								   CRIGHT(1:nof_variables)=U_c(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
								  
								  END IF
								  
								  
								  
								  
								  
								  
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
								  CALL coordinates_face_inner2dx(N,Iconsidered,facex,vext,nodes_list)
								    CORDS(1:2)=zero
								    N_NODE=2
								    CORDS(1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
							    
								    Poy(1)=cords(2)
								    Pox(1)=cords(1)
								    
								    
								    LEFTV(1:nof_variables)=CLEFT(1:nof_variables)
! 								    
								    B_CODE=ibound(n,ielem(n,i)%ibounds(l))%icode
								    CALL BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
								    cright(1:nof_Variables)=rightv(1:nof_Variables)
! 				  				   
				  				  	 IKAS=2			  				  
								    
								  END IF
							ELSE
							      CRIGHT(1:nof_Variables)=ILOCAL_RECON3(IELEM(N,I)%INEIGH(L))%ULEFT(1:nof_Variables,IELEM(N,I)%INEIGHN(L),NGP)
 							     
 							      IF ((CASCADE.EQ.2).AND.(IELEM(N,(IELEM(N,I)%INEIGH(L)))%MOOD.EQ.1))THEN
							      CRIGHT(1:nof_Variables)=U_c(IELEM(N,I)%INEIGH(L))%VAL(3,1:nof_variables)
							      
							      END IF
 							     
 							     
 							     
 							     
 							     
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
									  
									  
									  IF ((CASCADE.EQ.2).AND.(IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_M(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1).GT.0.5))THEN
                                    CRIGHT(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                                    
                                    END IF
									  
									  
									  
									  
 									   
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
 								  
 								  
 								   IF ((CASCADE.EQ.2).AND.(IEXBOUNDHIR(IELEM(N,I)%INEIGHN(L))%FACESOL_M(IELEM(N,I)%Q_FACE(L)%Q_MAPL(NGP),1).GT.0.5))THEN
                                    CRIGHT(1:nof_variables)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%INDEXI(L)))%SOL&
					(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%INDEXI(L)),1:nof_variables)
                                    
                                    END IF
 								  
 								  
 								  
 								  
 								  
 								  
 								  
 								  
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
				      
				      
			
						  CALL ROTATEF2d(N,CRIGHT_ROT,CRIGHT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
						  CALL ROTATEF2d(N,CLEFT_ROT,CLEFT,ANGLE1,ANGLE2)	!rotate wrt to normalvector of face and solve 1D Riemann problem
! 						   
						  
						  
						  IF ((LMACH.EQ.1))THEN    !application of the low mach number correction
						  LEFTV(1:nof_Variables)=CLEFT_ROT(1:nof_Variables); RIGHTV(1:nof_Variables)=CRIGHT_ROT(1:nof_Variables)
						  
						  CALL LMACHT2d(N,LEFTV,RIGHTV)
						  CLEFT_ROT(1:nof_Variables)=LEFTV(1:nof_Variables);CRIGHT_ROT(1:nof_Variables)=RIGHTV(1:nof_Variables);
						  
						  
						  END IF
						  
				      
								    IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN 
								      
								      cleft_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBL(1:turbulenceequations+PASSIVESCALAR)
								      cright_rot(nof_Variables+1:nof_variables+turbulenceequations+PASSIVESCALAR)=CTURBr(1:turbulenceequations+PASSIVESCALAR)
								    END IF
				      
				      
				      
				      
				      SELECT CASE(iRiemann)
				      
				      CASE(1)			!hllc
				      
				      CALL HLLC_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(2)			!rusanov
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				       CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then
				      
				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if
				      
				      CASE(9)			!hll

				      CALL HLL_RIEMANN_SOLVER2d(N,CLEFT_ROT,CRIGHT_ROT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				       if ((turbulence.eq.1).or.(passivescalar.gt.0))then

				      RHLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)=HLLCFLUX(NOF_VARIABLES+1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)
				      end if



				      CASE(3)			!roe
				      
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      CALL ROE_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY)
				      
				      RHLLCFLUX=HLLCFLUX
				     				      
				     
				        CASE(4)			!roe
				      
				      
				       
				      
				      IF ((B_CODE.le.0))THEN
				      
				      CALL ROTATEB2d(N,CLEFT,CLEFT_ROT,ANGLE1,ANGLE2)
				      CALL ROTATEB2d(N,CRIGHT,CRIGHT_ROT,ANGLE1,ANGLE2)
				      
				      
				      
				      CALL RROE_RIEMANN_SOLVER2d(N,Cleft,Cright,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT,NX,NY,b_code)
				      RHLLCFLUX=HLLCFLUX
				      ELSE
				      
!
				     
				      
				      CALL RUSANOV_RIEMANN_SOLVER2d(N,CLEFT,CRIGHT,HLLCFLUX,MP_SOURCE1,SRF_SPEEDROT)
				      CALL ROTATEB2d(N,RHLLCFLUX,HLLCFLUX,ANGLE1,ANGLE2)
				      
				      
				      
				      END IF
				      
				      
				       END SELECT
				       
				     
				      
				      GODFLUX2(1:nof_Variables)=GODFLUX2(1:nof_Variables)+(RHLLCFLUX(1:nof_Variables)*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L)))
				      IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE2=MP_SOURCE2+MP_SOURCE1*(WEIGHTS_TEMP(NGP)*IELEM(N,I)%SURF(L))
                        END IF
				       
				      
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
				   
				    RHS(I)%VAL(1:nof_Variables)=RHS(I)%VAL(1:nof_Variables)+GODFLUX2(1:nof_Variables)
				    IF (MULTISPECIES.EQ.1)THEN
                        MP_SOURCE3=MP_SOURCE3+MP_SOURCE2
                        END IF
! 				    if ((igoflux.eq.1))then
! 				    RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)=RHS(IELEM(N,I)%INEIGH(L))%VAL(1:nof_Variables)-GODFLUX2(1:nof_Variables)
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
		     IF (MULTISPECIES.EQ.1)THEN
                 RHS(I)%VAL(7)=RHS(I)%VAL(7)-(U_C(I)%VAL(1,7)*MP_SOURCE3)!*ielem(n,I)%totvolume)
                 
                 END IF
                 
                 END IF
	END DO
	!$OMP END DO

END SUBROUTINE CALCULATE_FLUXESHI_CONVECTIVE2d_MOOD


! !EXPERIMENTAL SUBROUTINE---NOT FOR PRODUCTION
! SUBROUTINE VERTEX_NEIGHBOURS_VALUES(N)
! IMPLICIT NONE
! INTEGER,INTENT(IN)::N
! INTEGER::I,J,K,L,M,NJ
! REAL,DIMENSION(4,30,NOF_VARIABLES)::VERTEX_NEIGHBOURS_VALS	!MAXIMUM ALLOCATION NOT EFFICIENT
!
! I=ICONSIDERED
!
! WRITE(680+N,*),"ELEMENT NUMBER",IELEM(N,I)%IHEXGL
!
! DO NJ=1,IELEM(N,I)%nonodes
!
! 	IF (ILOCAL_RECON3(I)%LOCAL.eq.1)then
! 		DO J=1,IELEM(N,I)%NOJECOUNT(NJ)
! 			VERTEX_NEIGHBOURS_VALS(NJ,J,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)))%val(1,1:nof_variables)
! 		!ILOCAL_RECON3(I)%IHEXL(1,L)	!LOCAL NUMBERING IN MY CPU
! 		!ILOCAL_RECON3(I)%IHEXG(1,L)	!GLOBAL NUMBERING IN MY CPU
! 			WRITE(680+N,*),"NODE",NJ,"NODE NEIGHBOUR",J,"LOCAL NUMBER",IELEM(N,I)%NODES_NEIGHBOURS(NJ,J), "VALUES",VERTEX_NEIGHBOURS_VALS(NJ,J,1:NOF_VARIABLES)
!
! 		END DO
!
! 	ELSE
! 			DO J=1,IELEM(N,I)%NOJECOUNT(NJ)
!
! 				IF (ILOCAL_RECON3(I)%IHEXB(1,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)).EQ.N)THEN
! 					VERTEX_NEIGHBOURS_VALS(NJ,J,1:NOF_VARIABLES)=U_C(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)))%val(1,1:nof_variables)
! 				ELSE
! 					write(500+n,*)ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)),IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)
! 					VERTEX_NEIGHBOURS_VALS(NJ,J,1:NOF_VARIABLES)=IEXSOLHIR(ILOCAL_RECON3(I)%IHEXN(1,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)))%SOL(ILOCAL_RECON3(I)%IHEXL(1,IELEM(N,I)%NODES_NEIGHBOURS(NJ,J)),1:nof_variables)
!
! 				END IF
! 			WRITE(680+N,*),"NODE",NJ,"NODE NEIGHBOUR",J,"LOCAL NUMBER",IELEM(N,I)%NODES_NEIGHBOURS(NJ,J), "VALUES",VERTEX_NEIGHBOURS_VALS(NJ,J,1:NOF_VARIABLES)
! 			END DO
! 	END IF
!
! END DO
!
! !NOW YOU HAVE COLLECTED FOR EVERY VERTEX THE VALUES
!
!
! END SUBROUTINE VERTEX_NEIGHBOURS_VALUES






END MODULE FLUXES
