MODULE SOURCE
USE LIBRARY
USE TRANSFORM
USE LOCAL
USE RIEMANN
USE FLOW_OPERATIONS
IMPLICIT NONE

contains
SUBROUTINE SOURCES_COMPUTATION(N)
!> @brief
!> Sources computation 
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE
	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE
		ICONSIDERED=I
		CALL SOURCES(N,ICONSIDERED)
		RHST(I)%VAL(1:turbulenceequations)=RHST(I)%VAL(1:turbulenceequations)-(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_COMPUTATION

SUBROUTINE SOURCES_derivatives_COMPUTATION(N)
!> @brief
!> Sources derivative computation for implicit time stepping
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE
	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE
		ICONSIDERED=I
		CALL SOURCES_derivatives(N,ICONSIDERED)
		sht(I,1:turbulenceequations)=(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_derivatives_COMPUTATION



SUBROUTINE SOURCES(N,ICONSIDERED)
!> @brief
!> Sources computation procedure
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
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


I=ICONSIDERED

Verysmall = TOLSMALL
	
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


LEFTV(1:5)=U_C(I)%VAL(1,1:5)

CALL CONS2PRIM(N)
RIGHTV(1:5)=LEFTV(1:5)
CALL SUTHERLAND(N,LEFTV,RIGHTV)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)



 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) )/(leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET/ SIGMA
			      SOURCE_T(1) = ProdTermfinal + fodt - destterm
			END IF
	 eLSE
		    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*KAPPA * KAPPA * (ddw) * (ddw))))


		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))**2.0)/( (ddw) * (ddw)))
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

			      CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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


SUBROUTINE SOURCES_DERIVATIVES(N,ICONSIDERED)
!> @brief
!> Sources derivative computation for implicit time stepping
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
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


I=ICONSIDERED

Verysmall = TOLSMALL
	
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


LEFTV(1:5)=U_C(I)%VAL(1,1:5)

CALL CONS2PRIM(N)
RIGHTV(1:5)=LEFTV(1:5)
CALL SUTHERLAND(N,LEFTV,RIGHTV)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)



 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
		    CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*KAPPA * KAPPA * (ddw) * (ddw))))


		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))**2.0)/( (ddw) * (ddw)))
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

			      CALL EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
	INTEGER::I,KMAXE
	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE
		ICONSIDERED=I
		CALL SOURCES2d(N,ICONSIDERED)
		RHST(I)%VAL(1:turbulenceequations)=RHST(I)%VAL(1:turbulenceequations)-(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_COMPUTATION2d

SUBROUTINE SOURCES_derivatives_COMPUTATION2d(N)
!> @brief
!> Sources  derivatives computation in 2D
	IMPLICIT NONE
	INTEGER,INTENT(IN)::N
	INTEGER::I,KMAXE
	
	
	KMAXE=XMPIELRANK(N)
	!$OMP DO SCHEDULE (STATIC)
	DO I=1,KMAXE
		ICONSIDERED=I
		CALL SOURCES_derivatives2d(N,ICONSIDERED)
		sht(I,1:turbulenceequations)=(SOURCE_T(1:turbulenceequations)*ielem(n,I)%totvolume)
	END DO
	!$OMP END DO 
END SUBROUTINE SOURCES_derivatives_COMPUTATION2d



SUBROUTINE SOURCES2d(N,ICONSIDERED)
!> @brief
!> Sources  computation in 2D
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
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


I=ICONSIDERED

Verysmall = TOLSMALL
	
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


LEFTV(1:4)=U_C(I)%VAL(1,1:4)
RIGHTV(1:4)=LEFTV(1:4)
CALL CONS2PRIM2d2(N)

CALL SUTHERLAND2d(N,LEFTV,RIGHTV)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)



 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) )/(leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET/ SIGMA
			      SOURCE_T(1) = ProdTermfinal + fodt - destterm
			     
			     
			END IF
	 eLSE
		    CALL EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*KAPPA * KAPPA * (ddw) * (ddw))))


		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))**2.0)/( (ddw) * (ddw)))
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

			      CALL EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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


SUBROUTINE SOURCES_DERIVATIVES2d(N,ICONSIDERED)
!> @brief
!> Sources derivatives computation in 2D
implicit none
INTEGER,INTENT(IN)::N,ICONSIDERED
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


I=ICONSIDERED

Verysmall = TOLSMALL
	
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


LEFTV(1:4)=U_C(I)%VAL(1,1:4)

CALL CONS2PRIM2d(N)
RIGHTV(1:4)=LEFTV(1:4)
CALL SUTHERLAND2d(N,LEFTV,RIGHTV)




TURBMV(1)=U_CT(I)%VAL(1,1)
TURBMV(2)=TURBMV(1)


 SELECT CASE(TURBULENCEMODEL)
 
 
 
 CASE(1) !!SPALART ALMARAS MODEL	

    
      if (ISPAL .eq.1) then
		eddyfl(2)=turbmv(1)
		eddyfr(2)=turbmv(2)
			
		onesix = 1.0D0/6.0D0

	    CALL EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
			gg	= rr +  ( CW2 * (rr**6 - rr) )
			Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
						!  Destruction term
			      destterm  = cw1 * fw * ((( TURBMV(1) ) /( leftv(1)*(ddw)))**2)
			      ! ! 	!  First order diffusion term
			      fodt  =   cb2 * SQUARET / SIGMA
			      
			      
			        Prodtermfinal=Stild*cb1
			         
			      destterm  = 2.0* cw1 * fw * (TURBMV(1)/(LEFTV(1)*(kappa**2)*(DDW**2)))
			      fodt  =   -2.0* cb2 * (SQRT(SQUARET))
			      		      
			      
			      
			      
			      IF (LAMZ.GT.0)THEN
 			      SOURCE_T(1)=  min(ProdTermfinal + fodt - destterm,ZERO)
			      ELSE
			      SOURCE_T(1)=  min(zero,- destterm +min(fodt,ZERO))
			      END IF
			      
			      
			END IF
	 eLSE
		    CALL EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
		  RR=((TURBMV(1)*TCH_X/(LEFTV(1)*KAPPA * KAPPA * (ddw) * (ddw))))


		  gg	= rr +  ( CW2 * (rr**6 - rr) )
		  Fw    = gg * (((1.0 + cw3**6) / (gg**6 + cw3**6))**onesix)
		  ! ! 				!  Destruction term

		  destterm  = cw1 * fw *leftv(1)* (((turbmv(1)*tch_x/leftv(1))**2.0)/( (ddw) * (ddw)))
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

			      CALL EDDYVISCO2d(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR)
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
