module riemann
use declaration
use library
use flow_operations



implicit none

 contains

subroutine exact_riemann_solver(n,cleft,cright,normalvect,hllcflux)
	implicit none
!> @brief
!> subroutine for linear advection ivp
#ifdef gpu
!$omp declare target
#endif
	real,dimension(1:nof_variables),intent(inout)::hllcflux
	integer,intent(in)::n
	real,intent(in)::normalvect
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cright
			hllcflux(1)=zero
			if (normalvect.gt.zero)then
				hllcflux(1)=(normalvect)*cleft(1)
			else
				hllcflux(1)=(normalvect)*cright(1)
			end if
end subroutine exact_riemann_solver






subroutine hll_riemann_solver(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	implicit none
!> @brief
!> hllc riemann solver in 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft_rot,cright_rot
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fl,fr
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fhll
	
	
	
	      
	hllcflux=zero
	rotvl=zero
	rotvr=zero
		tempfl=zero
		tempfr=zero
		flstar=zero
		frstar=zero
		ulstar=zero
		urstar=zero
		tempul=zero
		tempur=zero
		fl=zero
		fr=zero
		!conservative variables to primitive
		leftv(1:nof_variables)=cleft_rot(1:nof_variables)
		rightv(1:nof_variables)=cright_rot(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		rotvl(1:nof_variables)=leftv(1:nof_variables)
		rotvr(1:nof_variables)=rightv(1:nof_variables)
		
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then
		
		rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		
		end if
		
		!call estimate_waves(n,rotvl,rotvr,sl,sm,sr,gamma)
		
		
		
		
			!now conditions based on wave speeds!
			rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);wl=rotvl(4);pl=rotvl(5);el=cleft_rot(5)
			rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);wr=rotvr(4);pr=rotvr(5);er=cright_rot(5)
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			
			rml(1:0+turbulenceequations+passivescalar)=rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			rmr(1:0+turbulenceequations+passivescalar)=rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			

! 			if (turbulencemodel.eq.2)then
! 			pl=pl+((2.0d0/3.0d0)*eddyfl(2))
! 			pr=pr+((2.0d0/3.0d0)*eddyfr(2))  
! 
! 			end if
			end if
			
			
			
			
			
			if ((multispecies.eq.1).or.(realgas.eq.1))then
		ccl=sqrt(gammal*(pl+mp_pinfl)/rl)
		ccr=sqrt(gammar*(pr+mp_pinfr)/rr)
		else
		ccl=sqrt(gamma*pl/rl)
		ccr=sqrt(gamma*pr/rr)
		end if
		
		!einfeldt approximations
		cc2=sqrt(((((ccl**2)*sqrt(rl))+((ccr**2)*sqrt(rr)))/(sqrt(rl)+sqrt(rr)))+(0.5d0*((sqrt(rl)*sqrt(rr))/((sqrt(rl)+sqrt(rr))**2))*((ur-ul)**2)))
		uu2=(((ul*sqrt(rl))+(ur*sqrt(rr)))/(sqrt(rl)+sqrt(rr)))
		sl(1)=min(ul-ccl,uu2-cc2); sr(1)=max(ur+ccr,uu2+cc2)
 		sl(1)=min(sl(1),0.0d0); sr(1)=max(sr(1),0.0d0)
		sm(1)=(pr-pl+(rl*ul*(sl(1)-ul))-(rr*ur*(sr(1)-ur)))/((rl*(sl(1)-ul))-(rr*(sr(1)-ur)))
        




		

			fl(1)=rl*ul
			fl(2)=(rl*(ul**2))+pl
			fl(3)=rl*ul*vl
			fl(4)=rl*ul*wl
			fl(5)=ul*(el+pl)
			
			
			fr(1)=rr*ur
			fr(2)=(rr*(ur**2))+pr
			fr(3)=rr*ur*vr
			fr(4)=rr*ur*wr
			fr(5)=ur*(er+pr)



			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
			fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
			end if


			 if ((multispecies.eq.1).or.(realgas.eq.1))then
			 if (realgas.eq.1)then
			fl(6:nof_variables)=rotvl(6:nof_variables)*ul*rl
			fr(6:nof_variables)=rotvr(6:nof_variables)*ur*rr
			end if
			if (multispecies.eq.1)then
			fl(6:nof_variables)=rotvl(6:nof_variables)*ul
			fr(6:nof_variables)=rotvr(6:nof_variables)*ur
			end if

			end if

			

            sl(1)=min(ul-ccl,ur-ccr)
        sr(1)=max(ul+ccl,ur+ccr)

			
			!fhll=(sr(1)*fl(:)-sl(1)*fr(:)+(sl(1)*sr(1)*(cright_rot(:)-celft_rot(:))))/(sr(1)-sl(1))))
			
			
			fhll(:)=(sr(1)*fl(:)-sl(1)*fr(:)+(sl(1)*sr(1)*(cright_rot(:)-cleft_rot(:))))/(sr(1)-sl(1))

			

			
			if (sl(1).ge.zero)then
				hllcflux(:)=fl(:)

				if (multispecies.eq.1)then

                mp_source1=ul

                end if

			end if
			if (sr(1).le.zero)then
				hllcflux(:)=fr(:)


				if (multispecies.eq.1)then

                mp_source1=ur

                end if

			end if
			if ((sl(1).le.zero).and.(sr(1).ge.zero))then
				hllcflux(:)=fhll(:)


				if (multispecies.eq.1)then

                mp_source1=(ul+ur)/2



                end if

			end if

			
			


end subroutine hll_riemann_solver



subroutine hllc_riemann_solver(n,iconsidered, facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	implicit none
!> @brief
!> hllc riemann solver in 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,iconsidered, facex
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft_rot,cright_rot
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fl,fr
	real,dimension(turbulenceequations+passivescalar)::rml,rmr




	hllcflux=zero
	rotvl=zero
	rotvr=zero
		tempfl=zero
		tempfr=zero
		flstar=zero
		frstar=zero
		ulstar=zero
		urstar=zero
		tempul=zero
		tempur=zero
		fl=zero
		fr=zero
		!conservative variables to primitive
		leftv(1:nof_variables)=cleft_rot(1:nof_variables)
		rightv(1:nof_variables)=cright_rot(1:nof_variables)

		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

		rotvl(1:nof_variables)=leftv(1:nof_variables)
		rotvr(1:nof_variables)=rightv(1:nof_variables)

		if ((turbulence.eq.1).or.(passivescalar.gt.0))then

		rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

		end if

		!call estimate_waves(n,rotvl,rotvr,sl,sm,sr,gamma)




			!now conditions based on wave speeds!
			rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);wl=rotvl(4);pl=rotvl(5);el=cleft_rot(5)
			rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);wr=rotvr(4);pr=rotvr(5);er=cright_rot(5)

			if ((turbulence.eq.1).or.(passivescalar.gt.0))then

			rml(1:0+turbulenceequations+passivescalar)=rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			rmr(1:0+turbulenceequations+passivescalar)=rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)


! 			if (turbulencemodel.eq.2)then
! 			pl=pl+((2.0d0/3.0d0)*eddyfl(2))
! 			pr=pr+((2.0d0/3.0d0)*eddyfr(2))
!
! 			end if
			end if





			if ((multispecies.eq.1).or.(realgas.eq.1))then
		ccl=sqrt(gammal*(pl+mp_pinfl)/rl)
		ccr=sqrt(gammar*(pr+mp_pinfr)/rr)
		else
		ccl=sqrt(gamma*pl/rl)
		ccr=sqrt(gamma*pr/rr)
		end if

		!einfeldt approximations
		cc2=sqrt(((((ccl**2)*sqrt(rl))+((ccr**2)*sqrt(rr)))/(sqrt(rl)+sqrt(rr)))+(0.5d0*((sqrt(rl)*sqrt(rr))/((sqrt(rl)+sqrt(rr))**2))*((ur-ul)**2)))
		uu2=(((ul*sqrt(rl))+(ur*sqrt(rr)))/(sqrt(rl)+sqrt(rr)))
	sl(1)=min(ul-ccl,uu2-cc2); sr(1)=max(ur+ccr,uu2+cc2)
 		sl(1)=min(sl(1),0.0d0); sr(1)=max(sr(1),0.0d0)
		sm(1)=(pr-pl+(rl*ul*(sl(1)-ul))-(rr*ur*(sr(1)-ur)))/((rl*(sl(1)-ul))-(rr*(sr(1)-ur)))






			fl(1)=rl*ul
			fl(2)=(rl*(ul**2))+pl
			fl(3)=rl*ul*vl
			fl(4)=rl*ul*wl
			fl(5)=ul*(el+pl)


			fr(1)=rr*ur
			fr(2)=(rr*(ur**2))+pr
			fr(3)=rr*ur*vr
			fr(4)=rr*ur*wr
			fr(5)=ur*(er+pr)

			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
			fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
			end if

			mul=rl*((sl(1)-ul)/(sl(1)-sm(1)))
			mur=rr*((sr(1)-ur)/(sr(1)-sm(1)))
			lastl=(el/rl)+((sm(1)-ul)*(sm(1)+((pl)/(rl*(sl(1)-ul)))))
			lastr=(er/rr)+((sm(1)-ur)*(sm(1)+((pr)/(rr*(sr(1)-ur)))))
			ulstar(1)=mul
			ulstar(2)=mul*sm(1)
			ulstar(3)=mul*vl
			ulstar(4)=mul*wl
			ulstar(5)=mul*lastl

			if ((multispecies.eq.1).or.(realgas.eq.1))then

			 if (realgas.eq.1)then
			fl(6:nof_variables)=rotvl(6:nof_variables)*ul*rl
			ulstar(6:nof_variables)=mul*rotvl(6:nof_variables)
			end if
			if (multispecies.eq.1)then
			fl(6:nof_variables)=rotvl(6:nof_variables)*ul
            ulstar(6:nof_variables)=mul*rotvl(6:nof_variables)/rl
			end if
			end if

			urstar(1)=mur
			urstar(2)=mur*sm(1)
			urstar(3)=mur*vr
			urstar(4)=mur*wr
			urstar(5)=mur*lastr

			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			ulstar(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=mul*rml(1:0+turbulenceequations+passivescalar)/rl
			urstar(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=mur*rmr(1:0+turbulenceequations+passivescalar)/rr
			end if




			if ((multispecies.eq.1).or.(realgas.eq.1))then

			 if (realgas.eq.1)then
			fr(6:nof_variables)=rotvr(6:nof_variables)*ur*rr
			urstar(6:nof_variables)=mur*rotvr(6:nof_variables)
			end if
			if (multispecies.eq.1)then
			fr(6:nof_variables)=rotvr(6:nof_variables)*ur
            urstar(6:nof_variables)=mur*rotvr(6:nof_variables)/rr
			end if
			end if



			flstar(:)=fl(:)+sl(1)*(ulstar(:)-cleft_rot(:))
			frstar(:)=fr(:)+sr(1)*(urstar(:)-cright_rot(:))

			if (rec_mrf(iconsidered).eq.0)then

			if (sl(1).ge.zero)then
				hllcflux(:)=fl(:)
				if (multispecies.eq.1)then
                mp_source1=ul
                end if
			end if
			if (sr(1).le.zero)then
				hllcflux(:)=fr(:)
				if (multispecies.eq.1)then
                mp_source1=ur
                end if
			end if
			if ((sl(1).le.zero).and.(sm(1).ge.zero))then
				hllcflux(:)=flstar(:)
				if (multispecies.eq.1)then
                mp_source1=ul+sl(1)*(((sl(1)-ul)/(sl(1)-sm(1)))-1.0d0)
                end if
			end if
			if ((sr(1).ge.zero).and.(sm(1).le.zero))then
				hllcflux(:)=frstar(:)
				if (multispecies.eq.1)then
                mp_source1=ur+sr(1)*(((sr(1)-ur)/(sr(1)-sm(1)))-1.0d0)
                end if
			end if


			hllcflux(:)=(((1.0d0+sign(1.0d0,sm(1)))/2.0d0)*(fl(:)+sl(1)*(ulstar(:)-cleft_rot(:))))+&
			(((1.0d0-sign(1.0d0,sm(1)))/2.0d0)*(fr(:)+sr(1)*(urstar(:)-cright_rot(:))))


            else
                if ((sl(1)-srf_speedrot(2)).ge.zero)then
                        hllcflux(:)=fl(:)-srf_speedrot(2)*cleft_rot(:)
                end if
                if ((sr(1)-srf_speedrot(2)).le.zero)then
                        hllcflux(:)=fr(:)-srf_speedrot(2)*cright_rot(:)
                end if
                if (((sl(1)-srf_speedrot(2)).le.zero).and.((sm(1)-srf_speedrot(2)).ge.zero))then
                        hllcflux(:)=flstar(:)-srf_speedrot(2)*ulstar(:)
                end if
                if (((sr(1)-srf_speedrot(2)).ge.zero).and.((sm(1)-srf_speedrot(2)).le.zero))then
                        hllcflux(:)=frstar(:)-srf_speedrot(2)*urstar(:)
                end if
            end if

end subroutine hllc_riemann_solver






subroutine roe_riemann_solver(n,iconsidered, facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	implicit none
!> @brief
!> roe riemann solver in 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,iconsidered,facex
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft,cright
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real,intent(in)::nx,ny,nz
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	real :: sqrtrhol,sqrtrhor,hl,hr,utilde,vtilde,wtilde,htilde,atilde,veltilde
	
!originally obtained from katate masatsuka, february 2009. http://www.cfdbooks.com

!input
 real:: priml(5), primr(5) ! input: primitive variables
 real:: njk(3)             ! input: face normal vector

!output

!some constants
 real::   one = 1.0d0
 real::   two = 2.0d0
 real::  half = 0.5d0
 real:: fifth = 0.2d0

!local variables
  real:: eig(4)                         ! eigenvalues
 real:: rhol, rhor           ! primitive variables.
 real:: qnl, qnr                     ! normal velocities
 real:: al, ar             ! speed of sound, total enthalpy
 real:: rt,rho,u,v,w,h,a,qn          ! roe-averages
 real:: drho,dqn,dp,ldu(4)           ! wave strengths
 real:: du, dv, dw                   ! velocity differences
 real:: ws(4), r(5,4)                ! wave speeds and right-eigenvectors
 real:: dws(4)                       ! width of a parabolic fit for entropy fix
 real:: fl(5), fr(5), diss(5)        ! fluxes ad dissipation term
 real:: srp,slm                        ! wave speeds for the hll part
 real:: nx1, ny1, nz1                  ! vector along which hll is applied
 real:: nx2, ny2, nz2                  ! vector along which roe is applied
 real:: alpha1, alpha2                 ! projections of the new normals
 real:: abs_dq                         ! magnitude of the velocity difference
 real:: temp, tempx, tempy, tempz      ! temporary variables
 

! face normal vector (unit vector)

		leftv(1:nof_variables)=cleft(1:nof_variables)
		rightv(1:nof_variables)=cright(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		priml(1:nof_variables)=leftv(1:nof_variables)
		primr(1:nof_variables)=rightv(1:nof_variables)

!primitive and other variables.

!  left state
    rhol = priml(1)
      ul = priml(2)
      vl = priml(3)
      wl = priml(4)
     qnl = ul*nx + vl*ny + wl*nz
      pl = priml(5)
      al = sqrt(gamma*pl/rhol)
      hl = al*al/(gamma-one) + half*(ul*ul+vl*vl+wl*wl)
!  right state
    rhor = primr(1)
      ur = primr(2)
      vr = primr(3)
      wr = primr(4)
     qnr = ur*nx + vr*ny + wr*nz
      pr = primr(5)
      ar = sqrt(gamma*pr/rhor)
      hr = ar*ar/(gamma-one) + half*(ur*ur+vr*vr+wr*wr)

!first compute the roe-averaged quantities

!  note: see http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the roe-averaged density.

    rt = sqrt(rhor/rhol)
   rho = rt*rhol                                        !roe-averaged density
     u = (ul + rt*ur)/(one + rt)                        !roe-averaged x-velocity
     v = (vl + rt*vr)/(one + rt)                        !roe-averaged y-velocity
     w = (wl + rt*wr)/(one + rt)                        !roe-averaged z-velocity
     h = (hl + rt*hr)/(one + rt)                        !roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(h-half*(u*u + v*v + w*w)) ) !roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !roe-averaged face-normal velocity

!wave strengths

   drho = rhor - rhol !density difference
     dp =   pr - pl   !pressure difference
    dqn =  qnr - qnl  !normal velocity difference

  ldu(1) = (dp - rho*a*dqn )/(two*a*a) !left-moving acoustic wave strength
  ldu(2) =  drho - dp/(a*a)            !entropy wave strength
  ldu(3) = (dp + rho*a*dqn )/(two*a*a) !right-moving acoustic wave strength
  ldu(4) = rho                         !shear wave strength (not really, just a factor)

!absolute values of the wave speeds

  ws(1) = abs(qn-a) !left-moving acoustic wave
  ws(2) = abs(qn)   !entropy wave
  ws(3) = abs(qn+a) !right-moving acoustic wave
  ws(4) = abs(qn)   !shear waves

!harten's entropy fix jcp(1983), 49, pp357-393: only for the nonlinear fields.
!note: it avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!right eigenvectors
!note: two shear wave components are combined into one, so that tangent vectors
!      are not required. and that's why there are only 4 vectors here.
!      see "i do like cfd, vol.1" about how tangent vectors are eliminated.

! left-moving acoustic wave
  r(1,1) = one    
  r(2,1) = u - a*nx
  r(3,1) = v - a*ny
  r(4,1) = w - a*nz
  r(5,1) = h - a*qn

! entropy wave
  r(1,2) = one
  r(2,2) = u
  r(3,2) = v 
  r(4,2) = w
  r(5,2) = half*(u*u + v*v + w*w)

! right-moving acoustic wave
  r(1,3) = one
  r(2,3) = u + a*nx
  r(3,3) = v + a*ny
  r(4,3) = w + a*nz
  r(5,3) = h + a*qn

! two shear wave components combined into one (wave strength incorporated).
  du = ur - ul
  dv = vr - vl
  dw = wr - wl
  r(1,4) = zero
  r(2,4) = du - dqn*nx
  r(3,4) = dv - dqn*ny
  r(4,4) = dw - dqn*nz
  r(5,4) = u*du + v*dv + w*dw - qn*dqn

!dissipation term: |an|(ur-ul) = r|lambda|l*du = sum_k of [ ws(k) * r(:,k) * l*du(k) ]

 diss(:) = ws(1)*ldu(1)*r(:,1) + ws(2)*ldu(2)*r(:,2) &
         + ws(3)*ldu(3)*r(:,3) + ws(4)*ldu(4)*r(:,4)

!compute the physical flux: fl = fn(ul) and fr = fn(ur)

  fl(1) = rhol*qnl
  fl(2) = rhol*qnl * ul + pl*nx
  fl(3) = rhol*qnl * vl + pl*ny
  fl(4) = rhol*qnl * wl + pl*nz
  fl(5) = rhol*qnl * hl

  fr(1) = rhor*qnr
  fr(2) = rhor*qnr * ur + pr*nx
  fr(3) = rhor*qnr * vr + pr*ny
  fr(4) = rhor*qnr * wr + pr*nz
  fr(5) = rhor*qnr * hr

! this is the numerical flux: roe flux = 1/2 *[  fn(ul)+fn(ur) - |an|(ur-ul) ]

if (adda.eq.1)then

    hllcflux(1:nof_variables)= half * (fl(1:nof_variables) + fr(1:nof_variables) - diss(1:nof_variables)*ielem_facediss(facex,iconsidered))

else

    hllcflux(1:nof_variables)= half * (fl(1:nof_variables) + fr(1:nof_variables) - diss(1:nof_variables))

end if
!normal max wave speed in the normal direction.
!  wsn = abs(qn) + a



end subroutine roe_riemann_solver

subroutine troe_riemann_solver(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	implicit none
!> @brief
!> troe riemann solver in 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft,cright
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	real :: sqrtrhol,sqrtrhor,hl,hr,utilde,vtilde,wtilde,htilde,atilde,veltilde
	real,intent(in)::nx,ny,nz
	
!originally obtained from katate masatsuka, february 2009. http://www.cfdbooks.com

!input
 real:: priml(5), primr(5) ! input: primitive variables
 real:: njk(3)             ! input: face normal vector

!output

!some constants
 real::   one = 1.0d0
 real::   two = 2.0d0
 real::  half = 0.5d0
 real:: fifth = 0.2d0
real:: eig(4)                         ! eigenvalues
!local variables
  
 real:: rhol, rhor           ! primitive variables.
 real:: qnl, qnr                     ! normal velocities
 real:: al, ar             ! speed of sound, total enthalpy
 real:: rt,rho,u,v,w,h,a,qn          ! roe-averages
 real:: drho,dqn,dp,ldu(5)           ! wave strengths
 real:: du, dv, dw                   ! velocity differences
 real:: ws(5), r(5,5)                ! wave speeds and right-eigenvectors
 real:: dws(5)                       ! width of a parabolic fit for entropy fix
 real:: fl(5), fr(5), diss(5)        ! fluxes ad dissipation term
 real:: srp,slm                        ! wave speeds for the hll part
 real:: nx1, ny1, nz1                  ! vector along which hll is applied
 real:: nx2, ny2, nz2                  ! vector along which roe is applied
 real:: alpha1, alpha2                 ! projections of the new normals
 real:: abs_dq                         ! magnitude of the velocity difference
 real:: temp, tempx, tempy, tempz,lx,ly,lz,mx,my,mz, abs_n_cross_l, ql,qm,qll,qml,qlr,qmr,dql,dqm   ! temporary variables

! face normal vector (unit vector)

                leftv(1:nof_variables)=cleft(1:nof_variables)
		rightv(1:nof_variables)=cright(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		priml(1:nof_variables)=leftv(1:nof_variables)
		primr(1:nof_variables)=rightv(1:nof_variables)

!primitive and other variables.
   
!  left state
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
      ! impossible to happen
      stop
     endif

!     make it the unit vector.
      temp = sqrt( lx*lx + ly*ly + lz*lz )
       lx = lx/temp
       ly = ly/temp
       lz = lz/temp

! m = (mx,my,mz)
!
! the other one, m = (mx,my,mz), is chosen as a vector orthogonal to both n and l
! defined by the vector product: m = n x l / |n x l|

  mx = ny*lz - nz*ly
  my = nz*lx - nx*lz
  mz = nx*ly - ny*lx

  abs_n_cross_l = sqrt(mx**2 + my**2 + mz**2)
  mx = mx / abs_n_cross_l
  my = my / abs_n_cross_l
  mz = mz / abs_n_cross_l

!(do you like such ambiguous tangent vectors? actually, the roe flux can
! be implemented without any tangent vector. see "i do like cfd, vol.1",
! or the subroutine "inviscid_roe_n" or "inviscid_rotated_rhll" for details.
! in fact, the resulting flux is independent of the choice of these tangent vectors
! as it should be.)

!primitive and other variables.

!  left state
    rhol = priml(1)
      ul = priml(2)
      vl = priml(3)
      wl = priml(4)
     qnl = ul*nx + vl*ny + wl*nz
     qll = ul*lx + vl*ly + wl*lz
     qml = ul*mx + vl*my + wl*mz
      pl = priml(5)
      al = sqrt(gamma*pl/rhol)
      hl = al*al/(gamma-one) + half*(ul*ul+vl*vl+wl*wl)
!  right state
    rhor = primr(1)
      ur = primr(2)
      vr = primr(3)
      wr = primr(4)
     qnr = ur*nx + vr*ny + wr*nz
     qlr = ur*lx + vr*ly + wr*lz
     qmr = ur*mx + vr*my + wr*mz
      pr = primr(5)
      ar = sqrt(gamma*pr/rhor)
      hr = ar*ar/(gamma-one) + half*(ur*ur+vr*vr+wr*wr)

!first compute the roe-averaged quantities

!  note: see http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the roe-averaged density.

    rt = sqrt(rhor/rhol)
   rho = rt*rhol                                        !roe-averaged density
     u = (ul + rt*ur)/(one + rt)                        !roe-averaged x-velocity
     v = (vl + rt*vr)/(one + rt)                        !roe-averaged y-velocity
     w = (wl + rt*wr)/(one + rt)                        !roe-averaged z-velocity
     h = (hl + rt*hr)/(one + rt)                        !roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(h-half*(u*u + v*v + w*w)) ) !roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !roe-averaged face-normal velocity
    ql = u*lx + v*ly + w*lz                             !roe-averaged face-tangent velocity
    qm = u*mx + v*my + w*mz                             !roe-averaged face-tangent velocity

!wave strengths

   drho = rhor - rhol !density difference
     dp =   pr - pl   !pressure difference
    dqn =  qnr - qnl  !normal velocity difference
    dql =  qlr - qll  !tangent velocity difference in l
    dqm =  qmr - qml  !tangent velocity difference in m

  ldu(1) = (dp - rho*a*dqn )/(two*a*a) !left-moving acoustic wave strength
  ldu(2) =  drho - dp/(a*a)            !entropy wave strength
  ldu(3) = (dp + rho*a*dqn )/(two*a*a) !right-moving acoustic wave strength
  ldu(4) = rho*dql                     !shear wave strength
  ldu(5) = rho*dqm                     !shear wave strength

!absolute values of the wave speeds

  ws(1) = abs(qn-a) !left-moving acoustic wave speed
  ws(2) = abs(qn)   !entropy wave speed
  ws(3) = abs(qn+a) !right-moving acoustic wave speed
  ws(4) = abs(qn)   !shear wave speed
  ws(5) = abs(qn)   !shear wave speed

!harten's entropy fix jcp(1983), 49, pp357-393: only for the nonlinear fields.
!note: it avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!right eigenvectors

! left-moving acoustic wave
  r(1,1) = one    
  r(2,1) = u - a*nx
  r(3,1) = v - a*ny
  r(4,1) = w - a*nz
  r(5,1) = h - a*qn

! entropy wave
  r(1,2) = one
  r(2,2) = u
  r(3,2) = v 
  r(4,2) = w
  r(5,2) = half*(u*u + v*v + w*w)

! right-moving acoustic wave
  r(1,3) = one
  r(2,3) = u + a*nx
  r(3,3) = v + a*ny
  r(4,3) = w + a*nz
  r(5,3) = h + a*qn

! shear wave
  r(1,4) = zero
  r(2,4) = lx
  r(3,4) = ly
  r(4,4) = lz
  r(5,4) = ql

! shear wave
  r(1,5) = zero
  r(2,5) = mx
  r(3,5) = my
  r(4,5) = mz
  r(5,5) = qm

!dissipation term: |an|(ur-ul) = r|lambda|l*du = sum_k of [ ws(k) * r(:,k) * l*du(k) ]

 diss(:) = ws(1)*ldu(1)*r(:,1) + ws(2)*ldu(2)*r(:,2) + ws(3)*ldu(3)*r(:,3) &
         + ws(4)*ldu(4)*r(:,4) + ws(5)*ldu(5)*r(:,5)

!compute the physical flux: fl = fn(ul) and fr = fn(ur)

  fl(1) = rhol*qnl
  fl(2) = rhol*qnl * ul + pl*nx
  fl(3) = rhol*qnl * vl + pl*ny
  fl(4) = rhol*qnl * wl + pl*nz
  fl(5) = rhol*qnl * hl

  fr(1) = rhor*qnr
  fr(2) = rhor*qnr * ur + pr*nx
  fr(3) = rhor*qnr * vr + pr*ny
  fr(4) = rhor*qnr * wr + pr*nz
  fr(5) = rhor*qnr * hr

! this is the numerical flux: roe flux = 1/2 *[  fn(ul)+fn(ur) - |an|(ur-ul) ]

!   num_flux = half * (fl + fr - diss)
! this is the numerical flux: roe flux = 1/2 *[  fn(ul)+fn(ur) - |an|(ur-ul) ]

    hllcflux(1:nof_variables)= half * (fl + fr - diss)

!normal max wave speed in the normal direction.
!  wsn = abs(qn) + a



end subroutine troe_riemann_solver

subroutine rroe_riemann_solver(n,iconsidered,facex,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,nz)
	implicit none
!> @brief
!> rroe riemann solver in 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,iconsidered,facex
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft,cright
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	real :: sqrtrhol,sqrtrhor,hl,hr,utilde,vtilde,wtilde,htilde,atilde,veltilde
	real,intent(in)::nx,ny,nz
	
!originally obtained from katate masatsuka, february 2009. http://www.cfdbooks.com
!input
 real:: priml(5), primr(5) ! input: primitive variables
 real:: njk(3)             ! input: face normal vector

!output

!some constants
 real::   one = 1.0d0
 real::   two = 2.0d0
 real::  half = 0.5d0
 real:: fifth = 0.2d0
real:: eig(4)                         ! eigenvalues
!local variables
  
 real:: rhol, rhor           ! primitive variables.
 real:: qnl, qnr                     ! normal velocities
 real:: al, ar             ! speed of sound, total enthalpy
 real:: rt,rho,u,v,w,h,a,qn          ! roe-averages
 real:: drho,dqn,dp,ldu(4)           ! wave strengths
 real:: du, dv, dw                   ! velocity differences
 real:: ws(4), r(5,4)                ! wave speeds and right-eigenvectors
 real:: dws(4)                       ! width of a parabolic fit for entropy fix
 real:: fl(5), fr(5), diss(5)        ! fluxes ad dissipation term
 real:: srp,slm                        ! wave speeds for the hll part
 real:: nx1, ny1, nz1                  ! vector along which hll is applied
 real:: nx2, ny2, nz2                  ! vector along which roe is applied
 real:: alpha1, alpha2                 ! projections of the new normals
 real:: abs_dq                         ! magnitude of the velocity difference
 real:: temp, tempx, tempy, tempz      ! temporary variables

! face normal vector (unit vector)

                leftv(1:nof_variables)=cleft(1:nof_variables)
		rightv(1:nof_variables)=cright(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		priml(1:nof_variables)=leftv(1:nof_variables)
		primr(1:nof_variables)=rightv(1:nof_variables)

!primitive and other variables.
   
!  left state
    rhol = priml(1)
      ul = priml(2)
      vl = priml(3)
      wl = priml(4)
     qnl = ul*nx + vl*ny + wl*nz
      pl = priml(5)
      al = sqrt(gamma*pl/rhol)
      hl = al*al/(gamma-one) + half*(ul*ul+vl*vl+wl*wl)

!  right state
    rhor = primr(1)
      ur = primr(2)
      vr = primr(3)
      wr = primr(4)
     qnr = ur*nx + vr*ny + wr*nz
      pr = primr(5)
      ar = sqrt(gamma*pr/rhor)
      hr = ar*ar/(gamma-one) + half*(ur*ur+vr*vr+wr*wr)

!compute the physical flux: fl = fn(ul) and fr = fn(ur)

  fl(1) = rhol*qnl
  fl(2) = rhol*qnl * ul + pl*nx
  fl(3) = rhol*qnl * vl + pl*ny
  fl(4) = rhol*qnl * wl + pl*nz
  fl(5) = rhol*qnl * hl

  fr(1) = rhor*qnr
  fr(2) = rhor*qnr * ur + pr*nx
  fr(3) = rhor*qnr * vr + pr*ny
  fr(4) = rhor*qnr * wr + pr*nz
  fr(5) = rhor*qnr * hr

!--------------------------------------------------------------------------------
!define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
! note: n1 and n2 may need to be frozen at some point during 
!       a steady calculation to fully make it converge. for time-accurate 
!       calculation, it is not necessary.
! note: for a boundary face, you may want to set (nx2,ny2,nz2)=(nx,ny,nz),
!       (nx1,ny1)=tangent vector, i.e., use the roe flux.

    abs_dq = sqrt( (ur-ul)**2 + (vr-vl)**2 + (wr-wl)**2 )

!n1 = velocity difference vector: normal to shock or tangent to shear
  if ( abs_dq > tolsmall) then

       nx1 = (ur-ul)/abs_dq
       ny1 = (vr-vl)/abs_dq
       nz1 = (wr-wl)/abs_dq

!n1 = face tangent vector if abs_dq is too small.
!     note: there are infinitely many choices for the tangent vector.
!           the best choice may be discovered in future.
! here, we choose a vector in a plane (essentially 2d vector).
! note that we must be careful and make sure that the vector is not a zero vector.
  else

!  e.g., if n = (1,0,0), then, (0,-nz,ny) = (0,0,0). this is a problem.
!        to avoid such a situation, we choose the 2d tangential vector
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
      ! impossible to happen
      !write(*,*) "inviscid_rotated_rhll: impossible to happen. please report the problem."
!       stop
     endif

!     make it the unit vector.
      temp = sqrt( nx1*nx1 + ny1*ny1 + nz1*nz1 )
       nx1 = nx1/temp
       ny1 = ny1/temp
       nz1 = nz1/temp

  endif

    alpha1 = nx*nx1 + ny*ny1 + nz*nz1

! make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
       nz1 = temp * nz1
    alpha1 = temp * alpha1

!n2 = direction perpendicular to n1.
!     note: there are infinitely many choices for this vector.
!           the best choice may be discovered in future.
! here, we employ the formula (4.4) in the paper:
!     (nx2,ny2,nz2) = (n1xn)xn1 / |(n1xn)xn1|    ('x' is the vector product.)

!  (tempx,tempy,tempz) = n1xn
     tempx = ny1*nz - nz1*ny
     tempy = nz1*nx - nx1*nz
     tempz = nx1*ny - ny1*nx

!  (nx2,ny2,nz2) = (n1xn)xn1
     nx2 = tempy*nz1 - tempz*ny1
     ny2 = tempz*nx1 - tempx*nz1
     nz2 = tempx*ny1 - tempy*nx1

!  make n2 the unit vector
     temp = sqrt( nx2*nx2 + ny2*ny2 + nz2*nz2 )
       nx2 = nx2/temp
       ny2 = ny2/temp
       nz2 = nz2/temp

    alpha2 = nx*nx2 + ny*ny2 + nz*nz2

!  make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
       nz2 = temp * nz2
    alpha2 = temp * alpha2

!--------------------------------------------------------------------------------
!now we are going to compute the roe flux with n2 as the normal with modified 
!wave speeds (5.12). note: the roe flux here is computed without tangent vectors.
!see "i do like cfd, vol.1" for details: page 57, equation (3.6.31).

!first compute the roe-averaged quantities

!  note: see http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the roe-averaged density.

      rt = sqrt(rhor/rhol)
     rho = rt*rhol                                        !roe-averaged density.
       u = (ul + rt*ur)/(one + rt)                        !roe-averaged x-velocity
       v = (vl + rt*vr)/(one + rt)                        !roe-averaged y-velocity
       w = (wl + rt*wr)/(one + rt)                        !roe-averaged z-velocity
       h = (hl + rt*hr)/(one + rt)                        !roe-averaged total enthalpy
       a = sqrt( (gamma-one)*(h-half*(u*u + v*v + w*w)) ) !roe-averaged speed of sound

!----------------------------------------------------
!compute the wave speed estimates for the hll part,
!following einfeldt:
!
! b. einfeldt, on godunov-type methods for gas dynamics,
! siam journal on numerical analysis 25 (2) (1988) 294–318.
!
! note: hll is actually applied to n1, but this is
!       all we need to incorporate hll. see jcp2008 paper.

     qn  = u *nx1 + v *ny1 + w *nz1
     qnl = ul*nx1 + vl*ny1 + wl*nz1
     qnr = ur*nx1 + vr*ny1 + wr*nz1
     slm = min( zero, qn - a, qnl - al ) !minimum wave speed estimate
     srp = max( zero, qn + a, qnr + ar ) !maximum wave speed estimate

! this is the only place where n1=(nx1,ny1,nz1) is used.
! n1=(nx1,ny1,nz1) is never used below.
!----------------------------------------------------

!wave strengths

     qn  = u *nx2 + v *ny2 + w *nz2
     qnl = ul*nx2 + vl*ny2 + wl*nz2
     qnr = ur*nx2 + vr*ny2 + wr*nz2

    drho = rhor - rhol  !density difference
      dp =   pr - pl    !pressure difference
     dqn =  qnr - qnl   !normal velocity difference

  ldu(1) = (dp - rho*a*dqn )/(two*a*a) !left-moving acoustic wave strength
  ldu(2) =  drho - dp/(a*a)            !entropy wave strength
  ldu(3) = (dp + rho*a*dqn )/(two*a*a) !right-moving acoustic wave strength
  ldu(4) = rho                         !shear wave strength (not really, just a factor)

!wave speed (eigenvalues)

  eig(1) = qn-a !left-moving acoustic wave velocity
  eig(2) = qn   !entropy wave velocity
  eig(3) = qn+a !right-moving acoustic wave velocity
  eig(4) = qn   !shear wave velocity

!absolute values of the wave speeds (eigenvalues)

   ws(1) = abs(qn-a) !left-moving acoustic wave speed
   ws(2) = abs(qn)   !entropy wave speed
   ws(3) = abs(qn+a) !right-moving acoustic wave speed
   ws(4) = abs(qn)   !shear wave speed

!harten's entropy fix jcp(1983), 49, pp357-393: only for the nonlinear fields.
!note: it avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!combine the wave speeds for rotated-rhll: eq.(5.12) in the original jcp2008 paper.

      ws = alpha2*ws - (alpha1*two*srp*slm + alpha2*(srp+slm)*eig)/(srp-slm)

!below, we compute the roe dissipation term in the direction n2
!with the above modified wave speeds. hll wave speeds act something like
!the entropy fix or eigenvalue limiting; they contribute only by the amount
!given by the fraction, alpha1 (less than or equal to 1.0). see jcp2008 paper.

!right eigenvectors:
!note: two shear wave components are combined into one, so that tangent vectors
!      are not required. and that's why there are only 4 vectors here.

! left-moving acoustic wave
  r(1,1) = one    
  r(2,1) = u - a*nx2
  r(3,1) = v - a*ny2
  r(4,1) = w - a*nz2
  r(5,1) = h - a*qn

! entropy wave
  r(1,2) = one
  r(2,2) = u
  r(3,2) = v 
  r(4,2) = w
  r(5,2) = half*(u*u + v*v + w*w)

! right-moving acoustic wave
  r(1,3) = one
  r(2,3) = u + a*nx2
  r(3,3) = v + a*ny2
  r(4,3) = w + a*nz2
  r(5,3) = h + a*qn

! two shear wave components combined into one (wave strength incorporated).
  du = ur - ul
  dv = vr - vl
  dw = wr - wl
  r(1,4) = zero
  r(2,4) = du - dqn*nx2
  r(3,4) = dv - dqn*ny2
  r(4,4) = dw - dqn*nz2
  r(5,4) = u*du + v*dv + w*dw - qn*dqn

!dissipation term: roe dissipation with the modified wave speeds.
! |an|du = r|lambda|l*du = sum_k of [ ws(k) * r(:,k) * l*du(k) ], where n=n2.

 diss(:) = ws(1)*ldu(1)*r(:,1) + ws(2)*ldu(2)*r(:,2) &
         + ws(3)*ldu(3)*r(:,3) + ws(4)*ldu(4)*r(:,4)

!compute the rotated-rhll flux. (it looks like the hll flux with roe dissipation.)

!   num_flux = (srp*fl - slm*fr)/(srp-slm) - half*diss
! this is the numerical flux: roe flux = 1/2 *[  fn(ul)+fn(ur) - |an|(ur-ul) ]

    if (adda.eq.0)then

    hllcflux(1:nof_variables)= (srp*fl - slm*fr)/(srp-slm) - half*diss

    else
     hllcflux(1:nof_variables)= (srp*fl - slm*fr)/(srp-slm) - half*diss*ielem_facediss(facex,iconsidered)

    end if
!normal max wave speed in the normal direction.
!  wsn = abs(qn) + a



end subroutine rroe_riemann_solver


subroutine rusanov_riemann_solver(n,iconsidered,facex,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	implicit none
!> @brief
!> rusanov riemann solver in 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,iconsidered,facex
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft_rot,cright_rot
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fl,fr
	real,dimension(turbulenceequations+passivescalar)::rml,rmr

	
	
	
	      
	hllcflux=zero
	rotvl=zero
	rotvr=zero
		tempfl=zero
		tempfr=zero
		flstar=zero
		frstar=zero
		ulstar=zero
		urstar=zero
		tempul=zero
		tempur=zero
		fl=zero
		fr=zero
		!conservative variables to primitive
		leftv(1:nof_variables)=cleft_rot(1:nof_variables)
		rightv(1:nof_variables)=cright_rot(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		rotvl(1:nof_variables)=leftv(1:nof_variables)
		rotvr(1:nof_variables)=rightv(1:nof_variables)
		
		
		
		
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then
		
		rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		
		end if
		
		!call estimate_waves(n,rotvl,rotvr,sl,sm,sr,gamma)
        if (rec_mrf(iconsidered).eq.0)then



					!now conditions based on wave speeds!
					rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);wl=rotvl(4);pl=rotvl(5);el=cleft_rot(5)
					rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);wr=rotvr(4);pr=rotvr(5);er=cright_rot(5)
			
					if ((turbulence.eq.1).or.(passivescalar.gt.0))then

					rml(1:0+turbulenceequations+passivescalar)=rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
					rmr(1:0+turbulenceequations+passivescalar)=rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)


		! 			if (turbulencemodel.eq.2)then
		! 			pl=pl+((2.0d0/3.0d0)*eddyfl(2))
		! 			pr=pr+((2.0d0/3.0d0)*eddyfr(2))
		!
		! 			end if
					end if


					fl(1)=rl*ul
					fl(2)=(rl*(ul**2))+pl
					fl(3)=rl*ul*vl
					fl(4)=rl*ul*wl
					fl(5)=ul*(el+pl)
					if ((multispecies.eq.1).or.(realgas.eq.1))then
					if (realgas.eq.1)then
					fl(6:nof_variables)=rotvl(6:nof_variables)*ul*rl
					end if
					if (multispecies.eq.1)then
					fl(6:nof_variables)=rotvl(6:nof_variables)*ul
					end if
					end if

					fr(1)=rr*ur
					fr(2)=(rr*(ur**2))+pr
					fr(3)=rr*ur*vr
					fr(4)=rr*ur*wr
					fr(5)=ur*(er+pr)

					if ((multispecies.eq.1).or.(realgas.eq.1))then
					if (realgas.eq.1)then
					fr(6:nof_variables)=rotvr(6:nof_variables)*ur*rr
					end if
					if (multispecies.eq.1)then
					fr(6:nof_variables)=rotvr(6:nof_variables)*ur
					end if
					end if

					if ((turbulence.eq.1).or.(passivescalar.gt.0))then
					fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
					fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
					end if
			
			
			
			

			
					if ((multispecies.eq.1).or.(realgas.eq.1))then
					sl(1)=abs(ul)+sqrt(gammal*(pl+mp_pinfl)/rl)
					sr(1)=abs(ur)+sqrt(gammar*(pr+mp_pinfr)/rr)
					mp_source1=0.5d0*(ul+ur)!(max(abs(sl(1)),abs(sr(1))))
					else
					sl(1)=abs(ul)+sqrt(gamma*pl/rl)
					sr(1)=abs(ur)+sqrt(gamma*pr/rr)
					endif
			
					if (adda.eq.0)then
					hllcflux(:)=0.5d0*(fl(:)+fr(:))-0.5d0*max(abs(sl(1)),abs(sr(1)))*(cright_rot(:)-cleft_rot(:))

					else
					hllcflux(:)=0.5d0*(fl(:)+fr(:))-0.5d0*max(abs(sl(1)),abs(sr(1)))*ielem_facediss(facex,iconsidered)*(cright_rot(:)-cleft_rot(:))
					end if
			

			else


			!now conditions based on wave speeds!
			rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);wl=rotvl(4);pl=rotvl(5);el=cleft_rot(5)
			rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);wr=rotvr(4);pr=rotvr(5);er=cright_rot(5)
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			
			rml(1:0+turbulenceequations+passivescalar)=rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			rmr(1:0+turbulenceequations+passivescalar)=rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			
			end if

			fl(1)=rl*ul
			fl(2)=(rl*(ul**2))+pl
			fl(3)=rl*ul*vl
			fl(4)=rl*ul*wl
			fl(5)=ul*(el+pl)
			fr(1)=rr*ur
			fr(2)=(rr*(ur**2))+pr
			fr(3)=rr*ur*vr
			fr(4)=rr*ur*wr
			fr(5)=ur*(er+pr)
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
			fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
			end if
 			fl(:)=fl(:)-srf_speedrot(2)*(cleft_rot(:))
 			fr(:)=fr(:)-srf_speedrot(2)*(cright_rot(:))
			hllcflux(:)=0.5d0*(fl(:)+fr(:))-0.5d0*max(abs(sl(1)-srf_speedrot(2)),abs(sr(1)-srf_speedrot(2)))*(cright_rot(:)-cleft_rot(:))
			
			
        end if


end subroutine rusanov_riemann_solver




subroutine estimate_waves(n,rotvl,rotvr,sl,sm,sr)
	implicit none
!> @brief
!> wave speed estimates for 3d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::rotvl,rotvr
	real,dimension(1),intent(inout)::sl,sr,sm
	real::cl,cr,pr,pl,ul,ur,vl,vr,wr,wl,rl,rr
	real::cup,ppv,pmin,pmax,qmax,quser,bl,br,cov,pm,um
	real::g1,g2,g3,g4,g5,g6,g7,g8,gel,ger,pq,ptl,ptr
	g1 = (gamma - 1.0d0)/(2.0d0*gamma)
    	g2 = (gamma + 1.0d0)/(2.0d0*gamma)
   	g3 = 2.0d0*gamma/(gamma - 1.0d0)
    	g4 = 2.0d0/(gamma - 1.0d0)
    	g5 = 2.0d0/(gamma + 1.0d0)
    	g6 = (gamma - 1.0d0)/(gamma + 1.0d0)
   	g7 = (gamma - 1.0d0)/2.0d0
   	g8 = gamma - 1.0d0
	
	sl=0.0d0
	sr=0.0d0
	sm=0.0d0
	cov=0.0d0
	!build left state variables
	rl=rotvl(1)
	ul=rotvl(2)
	vl=rotvl(3)
	wl=rotvl(4)
	pl=rotvl(5)
	cl=sqrt((pl*gamma)/(rl))
	!build right  state variables
	rr=rotvr(1)
	ur=rotvr(2)
	vr=rotvr(3)
	wr=rotvr(4)
	pr=rotvr(5)
	cr=sqrt((pr*gamma)/(rr))

	cup=0.25d0*(rl+rr)*(cl+cr)
	ppv=0.5d0*(pl + pr) + 0.5d0*(ul - ur)*cup
	ppv=max(0.0d0,ppv)
	pmin=min(pl,pr)
	pmax=max(pl,pr)
	qmax=pmax/pmin
	quser=2.0d0

	 if(qmax.le.quser.and.(pmin.le.ppv.and.ppv.le.pmax))then
  
!        select prvs riemann solver
 
         pm = ppv
         um = 0.5d0*(ul + ur) + 0.5d0*(pl - pr)/cup 
  
      	else
 
         bl = 1.0d0 - cov*rl
         br = 1.0d0 - cov*rr
 
         if(ppv.lt.pmin)then
 
!           select two-rarefaction riemann solver
 
        
            pq  = (pl/pr)**g1
            um  = (pq*ul/cl/bl + ur/cr/br + g4*(pq - 1.0d0)) 
            um  = um/(pq/cl/bl + 1.0d0/cr/br)
            ptl = 1.0d0 + g7*(ul - um)/cl/bl
            ptr = 1.0d0 + g7*(um - ur)/cr/br
            pm  = 0.5d0*(pl*ptl**g3 + pr*ptr**g3)
         else

!           use two-shock riemann solver with pvrs as estimate
 
!           introduce iterations with pvrs as initial guess
 
            do k=1,4
             gel = sqrt((g5*bl/rl)/(g6*pl + ppv))
             ger = sqrt((g5*br/rr)/(g6*pr + ppv))
             pm  = (gel*pl + ger*pr - (ur - ul))/(gel + ger)
             um  = 0.5d0*(ul + ur) + 0.5d0*(ger*(pm - pr) - gel*(pm - pl))	     
             if ( abs((pm-ppv)/pm) .le. 1d-8) goto 101
                ppv = pm
     	    end do
         endif
      endif
      
      101 continue

!     find speeds
 
      if(pm.le.pl)then
         sl(1) = ul - cl
      else
         sl(1) = ul - cl*sqrt(1.0d0 + g2*(pm/pl - 1.0d0))
      endif
 
      sm(1)= um
 
      if(pm.le.pr)then
         sr(1)= ur + cr
      else
         sr(1) = ur + cr*sqrt(1.0d0 + g2*(pm/pr - 1.0d0))
      endif

end subroutine estimate_waves



subroutine hllc_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	implicit none
!> @brief
!> hllc riemann solver in 2d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k,riem_adapt
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft_rot,cright_rot
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fl,fr
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	
	
	
	      
	hllcflux=zero
	rotvl=zero
	rotvr=zero
		tempfl=zero
		tempfr=zero
		flstar=zero
		frstar=zero
		ulstar=zero
		urstar=zero
		tempul=zero
		tempur=zero
		fl=zero
		fr=zero
		!conservative variables to primitive
		leftv(1:nof_variables)=cleft_rot(1:nof_variables)
		rightv(1:nof_variables)=cright_rot(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		rotvl(1:nof_variables)=leftv(1:nof_variables)
		rotvr(1:nof_variables)=rightv(1:nof_variables)
		
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then
		
		rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		
		end if
		
! 		call estimate_waves2d(n,rotvl,rotvr,sl,sm,sr,gamma)
		
		
		
		
		
		
		
		
		
		
		
			!now conditions based on wave speeds!
			rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);pl=rotvl(4);el=cleft_rot(4)
			rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);pr=rotvr(4);er=cright_rot(4)
			
			
			
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		ccl=sqrt(gammal*(pl+mp_pinfl)/rl)
		ccr=sqrt(gammar*(pr+mp_pinfr)/rr)
		else
		ccl=sqrt(gamma*pl/rl)
		ccr=sqrt(gamma*pr/rr)
		end if	
			
! 		if (multispecies.eq.1)then
! 		sl(1)=min((ul)-sqrt(gammal*(pl+mp_pinfl)/rl),(ur)-sqrt(gammar*(pr+mp_pinfr)/rr))
! 		sr(1)=max((ul)+sqrt(gammal*(pl+mp_pinfl)/rl),(ur)+sqrt(gammar*(pr+mp_pinfr)/rr))	
! 		else
!         sl(1)=min((ul)-sqrt(gamma*pl/rl),(ur)-sqrt(gamma*pr/rr))
! 		sr(1)=max((ul)+sqrt(gamma*pl/rl),(ur)+sqrt(gamma*pr/rr))	
! 		end if
! 		sm(1)=(pr-pl+(rl*ul*(sl(1)-ul))-(rr*ur*(sr(1)-ur)))/((rl*(sl(1)-ul))-(rr*(sr(1)-ur)))
			
			
		cc2=sqrt(((((ccl**2)*sqrt(rl))+((ccr**2)*sqrt(rr)))/(sqrt(rl)+sqrt(rr)))+(0.5d0*((sqrt(rl)*sqrt(rr))/((sqrt(rl)+sqrt(rr))**2))*((ur-ul)**2)))
		uu2=(((ul*sqrt(rl))+(ur*sqrt(rr)))/(sqrt(rl)+sqrt(rr)))
		sl(1)=min(ul-ccl,uu2-cc2);
		sr(1)=max(ur+ccr,uu2+cc2)
  		sl(1)=min(sl(1),0.0d0); sr(1)=max(sr(1),0.0d0)
		sm(1)=(pr-pl+(rl*ul*(sl(1)-ul))-(rr*ur*(sr(1)-ur)))/((rl*(sl(1)-ul))-(rr*(sr(1)-ur)))	
			
			
			
			
			
			
			
			
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			
			rml(1:0+turbulenceequations+passivescalar)=rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			rmr(1:0+turbulenceequations+passivescalar)=rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			

! 			if (turbulencemodel.eq.2)then
! 			pl=pl+((2.0d0/3.0d0)*eddyfl(2))
! 			pr=pr+((2.0d0/3.0d0)*eddyfr(2))  
! 
! 			end if
			end if


			fl(1)=rl*ul
			fl(2)=(rl*(ul**2))+pl
			fl(3)=rl*ul*vl
			
			fl(4)=ul*(el+pl)
			
			
			fr(1)=rr*ur
			fr(2)=(rr*(ur**2))+pr
			fr(3)=rr*ur*vr
			
			fr(4)=ur*(er+pr)
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
			fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
			end if
			
			mul=rl*((sl(1)-ul)/(sl(1)-sm(1)))
			mur=rr*((sr(1)-ur)/(sr(1)-sm(1)))
			lastl=(el/rl)+((sm(1)-ul)*(sm(1)+((pl)/(rl*(sl(1)-ul)))))
			lastr=(er/rr)+((sm(1)-ur)*(sm(1)+((pr)/(rr*(sr(1)-ur)))))
			ulstar(1)=mul
			ulstar(2)=mul*sm(1)
			ulstar(3)=mul*vl
			
			ulstar(4)=mul*lastl
			
			if ((multispecies.eq.1).or.(realgas.eq.1))then

			 if (realgas.eq.1)then
			fl(5:nof_variables)=rotvl(5:nof_variables)*ul*rl
			ulstar(5:nof_variables)=mul*rotvl(5:nof_variables)
			end if
			if (multispecies.eq.1)then
			fl(5:nof_variables)=rotvl(5:nof_variables)*ul
            ulstar(5:nof_variables)=mul*rotvl(5:nof_variables)/rl
			end if
			end if
			
			

			
			urstar(1)=mur
			urstar(2)=mur*sm(1)
			urstar(3)=mur*vr
			
			urstar(4)=mur*lastr
			
			if ((multispecies.eq.1).or.(realgas.eq.1))then

			 if (realgas.eq.1)then
			fr(5:nof_variables)=rotvr(5:nof_variables)*ur*rr
			urstar(5:nof_variables)=mur*rotvr(5:nof_variables)
			end if
			if (multispecies.eq.1)then
			fr(5:nof_variables)=rotvr(5:nof_variables)*ur
            urstar(5:nof_variables)=mur*rotvr(5:nof_variables)/rr
			end if
			end if
			
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			ulstar(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=mul*rml(1:0+turbulenceequations+passivescalar)/rl
			urstar(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=mur*rmr(1:0+turbulenceequations+passivescalar)/rr
			end if
			
			flstar(:)=fl(:)+sl(1)*(ulstar(:)-cleft_rot(:))
			frstar(:)=fr(:)+sr(1)*(urstar(:)-cright_rot(:))
			
			if (sl(1).ge.zero)then
				hllcflux(:)=fl(:)
				if (multispecies.eq.1)then

                mp_source1=ul

                end if
				
			end if
			if (sr(1).le.zero)then
				hllcflux(:)=fr(:)
				if (multispecies.eq.1)then
                mp_source1=ur
                end if
			end if
			if ((sl(1).le.zero).and.(sm(1).ge.zero))then
				hllcflux(:)=flstar(:)
				if (multispecies.eq.1)then


                mp_source1=ul+sl(1)*(((sl(1)-ul)/(sl(1)-sm(1)))-1.0d0)

                end if
			end if
			if ((sr(1).ge.zero).and.(sm(1).le.zero))then
				hllcflux(:)=frstar(:)
				if (multispecies.eq.1)then

                mp_source1=ur+sr(1)*(((sr(1)-ur)/(sr(1)-sm(1)))-1.0d0)

                end if
			end if
			
			
			
			hllcflux(:)=(((1.0d0+sign(1.0d0,sm(1)))/2.0d0)*(fl(:)+sl(1)*(ulstar(:)-cleft_rot(:))))+&
			(((1.0d0-sign(1.0d0,sm(1)))/2.0d0)*(fr(:)+sr(1)*(urstar(:)-cright_rot(:))))
			
			
			
			!pgrad=abs(pl-pr)/min(pl,pr)
			!om_p=0.5d0-0.5d0*sign(pgrad-0.2,1.0d0)*(1.0-exp(-100.0d0*abs(pgrad-0.2)))
			
			!if(om_p.lt.0.9d0)then
			
			!else
			
			!sl(1)=abs(ul)+sqrt(gamma*pl/rl)
			!sr(1)=abs(ur)+sqrt(gamma*pr/rr)
			!hllcflux(:)=0.5d0*(fl(:)+fr(:))-0.5d0*max(abs(sl(1)),abs(sr(1)))*(cright_rot(:)-cleft_rot(:))
			!end if
			
			
			
			
			
			
			

end subroutine hllc_riemann_solver2d

subroutine hll_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	implicit none
!> @brief
!> hllc riemann solver in 2d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k,riem_adapt
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft_rot,cright_rot
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fl,fr
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fhll





	hllcflux=zero
	rotvl=zero
	rotvr=zero
		tempfl=zero
		tempfr=zero
		flstar=zero
		frstar=zero
		ulstar=zero
		urstar=zero
		tempul=zero
		tempur=zero
		fl=zero
		fr=zero
		!conservative variables to primitive
		leftv(1:nof_variables)=cleft_rot(1:nof_variables)
		rightv(1:nof_variables)=cright_rot(1:nof_variables)

		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)

		rotvl(1:nof_variables)=leftv(1:nof_variables)
		rotvr(1:nof_variables)=rightv(1:nof_variables)

		if ((turbulence.eq.1).or.(passivescalar.gt.0))then

		rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)

		end if

! 		call estimate_waves2d(n,rotvl,rotvr,sl,sm,sr,gamma)











			!now conditions based on wave speeds!
			rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);pl=rotvl(4);el=cleft_rot(4)
			rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);pr=rotvr(4);er=cright_rot(4)



		if ((multispecies.eq.1).or.(realgas.eq.1))then
		ccl=sqrt(gammal*(pl+mp_pinfl)/rl)
		ccr=sqrt(gammar*(pr+mp_pinfr)/rr)
		else
		ccl=sqrt(gamma*pl/rl)
		ccr=sqrt(gamma*pr/rr)
		end if

! 		if (multispecies.eq.1)then
! 		sl(1)=min((ul)-sqrt(gammal*(pl+mp_pinfl)/rl),(ur)-sqrt(gammar*(pr+mp_pinfr)/rr))
! 		sr(1)=max((ul)+sqrt(gammal*(pl+mp_pinfl)/rl),(ur)+sqrt(gammar*(pr+mp_pinfr)/rr))
! 		else
!         sl(1)=min((ul)-sqrt(gamma*pl/rl),(ur)-sqrt(gamma*pr/rr))
! 		sr(1)=max((ul)+sqrt(gamma*pl/rl),(ur)+sqrt(gamma*pr/rr))
! 		end if
! 		sm(1)=(pr-pl+(rl*ul*(sl(1)-ul))-(rr*ur*(sr(1)-ur)))/((rl*(sl(1)-ul))-(rr*(sr(1)-ur)))






			fl(1)=rl*ul
			fl(2)=(rl*(ul**2))+pl
			fl(3)=rl*ul*vl

			fl(4)=ul*(el+pl)


			fr(1)=rr*ur
			fr(2)=(rr*(ur**2))+pr
			fr(3)=rr*ur*vr

			fr(4)=ur*(er+pr)

			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
			fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
			end if



			     sl(1)=min(ul-ccl,ur-ccr)
        sr(1)=max(ul+ccl,ur+ccr)


        if ((multispecies.eq.1).or.(realgas.eq.1))then
			 if (realgas.eq.1)then
			fl(5:nof_variables)=rotvl(5:nof_variables)*ul*rl
			fr(5:nof_variables)=rotvr(5:nof_variables)*ur*rr
			end if
			if (multispecies.eq.1)then
			fl(5:nof_variables)=rotvl(5:nof_variables)*ul
			fr(5:nof_variables)=rotvr(5:nof_variables)*ur
			end if

			end if


			fhll(:)=(sr(1)*fl(:)-sl(1)*fr(:)+(sl(1)*sr(1)*(cright_rot(:)-cleft_rot(:))))/(sr(1)-sl(1))







			if (sl(1).ge.zero)then
				hllcflux(:)=fl(:)

				if (multispecies.eq.1)then

                mp_source1=ul

                end if

			end if
			if (sr(1).le.zero)then
				hllcflux(:)=fr(:)

				if (multispecies.eq.1)then

                mp_source1=ur

                end if

			end if
			if ((sl(1).le.zero).and.(sr(1).ge.zero))then
				hllcflux(:)=fhll(:)

				if (multispecies.eq.1)then

                mp_source1=(ul+ur)/2



                end if

			end if





			!pgrad=abs(pl-pr)/min(pl,pr)
			!om_p=0.5d0-0.5d0*sign(pgrad-0.2,1.0d0)*(1.0-exp(-100.0d0*abs(pgrad-0.2)))

			!if(om_p.lt.0.9d0)then

			!else

			!sl(1)=abs(ul)+sqrt(gamma*pl/rl)
			!sr(1)=abs(ur)+sqrt(gamma*pr/rr)
			!hllcflux(:)=0.5d0*(fl(:)+fr(:))-0.5d0*max(abs(sl(1)),abs(sr(1)))*(cright_rot(:)-cleft_rot(:))
			!end if








end subroutine hll_riemann_solver2d

subroutine roe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny)
	implicit none
!> @brief
!> roe riemann solver in 2d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft,cright
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real,intent(in)::nx,ny
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,el,er,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	real :: sqrtrhol,sqrtrhor,utilde,vtilde,wtilde,htilde,atilde,veltilde
	real :: ul(4), ur(4) !  input: conservative variables rho*[1, u, v, e]
 
 real :: roe(4)       ! output: roe flux function (upwind)
!local constants
                      ! ratio of specific heat.
 real ::  fifth, half, one, two    ! numbers
!local variables
 real :: tx, ty       ! tangent vector (perpendicular to the face normal)
 real :: vxl, vxr, vyl, vyr             ! velocity components.
 real :: rhol, rhor, pl, pr             ! primitive variables.
 real :: vnl, vnr, vtl, vtr             ! normal and tangent velocities
 real :: al, ar, hl, hr                 ! speeds of sound.
 real :: rt,rho,vx,vy,h,a,vn, vt        ! roe-averages
 real :: drho,dvx,dvy,dvn,dvt,dp,dv(4)  ! wave strenghs
 real :: ws(4),dws(4), rv(4,4)          ! wave speeds and right-eigevectors
 real :: fl(4), fr(4), diss(4)          ! fluxes ad dissipation term
 integer ::  j
!originally obtained from katate masatsuka, february 2009. http://www.cfdbooks.com
!constants.

      
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0

       ul(1:nof_variables)=cleft(1:nof_variables)
       ur(1:nof_variables)=cright(1:nof_variables)
       
       
!tangent vector (do you like it? actually, roe flux can be implemented 
! without any tangent vector. see "i do like cfd, vol.1" for details.)
  tx = -ny
  ty = nx

!primitive and other variables.
!  left state
    rhol = ul(1)
     vxl = ul(2)/ul(1)
     vyl = ul(3)/ul(1)
     vnl = vxl*nx+vyl*ny
     vtl = vxl*tx+vyl*ty
      pl = (gamma-one)*( ul(4) - half*rhol*(vxl*vxl+vyl*vyl) )
      al = sqrt(gamma*pl/rhol)
      hl = ( ul(4) + pl ) / rhol
!  right state
    rhor = ur(1)
     vxr = ur(2)/ur(1)
     vyr = ur(3)/ur(1)
     vnr = vxr*nx+vyr*ny
     vtr = vxr*tx+vyr*ty
      pr = (gamma-one)*( ur(4) - half*rhor*(vxr*vxr+vyr*vyr) )
      ar = sqrt(gamma*pr/rhor)
      hr = ( ur(4) + pr ) / rhor

!first compute the roe averages
    rt = sqrt(rhor/rhol)
   rho = rt*rhol
    vx = (vxl+rt*vxr)/(one+rt)
    vy = (vyl+rt*vyr)/(one+rt)
     h = ( hl+rt* hr)/(one+rt)
     a = sqrt( (gamma-one)*(h-half*(vx*vx+vy*vy)) )
    vn = vx*nx+vy*ny
    vt = vx*tx+vy*ty

!wave strengths
   drho = rhor - rhol 
     dp =   pr - pl
    dvn =  vnr - vnl
    dvt =  vtr - vtl

  dv(1) = (dp - rho*a*dvn )/(two*a*a)
  dv(2) = rho*dvt/a
  dv(3) =  drho - dp/(a*a)
  dv(4) = (dp + rho*a*dvn )/(two*a*a)

!wave speed
  ws(1) = abs(vn-a)
  ws(2) = abs(vn)
  ws(3) = abs(vn)
  ws(4) = abs(vn+a)

!harten's entropy fix jcp(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!right eigenvectors
  rv(1,1) = one    
  rv(2,1) = vx - a*nx
  rv(3,1) = vy - a*ny
  rv(4,1) =  h - vn*a

  rv(1,2) = zero
  rv(2,2) = a*tx
  rv(3,2) = a*ty
  rv(4,2) = vt*a

  rv(1,3) = one
  rv(2,3) = vx
  rv(3,3) = vy 
  rv(4,3) = half*(vx*vx+vy*vy)

  rv(1,4) = one
  rv(2,4) = vx + a*nx
  rv(3,4) = vy + a*ny
  rv(4,4) =  h + vn*a

!dissipation term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dv(j)*rv(i,j)
   end do
  end do

!compute the flux.
  fl(1) = rhol*vnl
  fl(2) = rhol*vnl * vxl + pl*nx
  fl(3) = rhol*vnl * vyl + pl*ny
  fl(4) = rhol*vnl *  hl

  fr(1) = rhor*vnr
  fr(2) = rhor*vnr * vxr + pr*nx
  fr(3) = rhor*vnr * vyr + pr*ny
  fr(4) = rhor*vnr *  hr

  hllcflux(1:nof_variables)= half * (fl(1:nof_variables) + fr(1:nof_variables) - diss(1:nof_variables))





end subroutine roe_riemann_solver2d


subroutine rroe_riemann_solver2d(n,cleft,cright,hllcflux,mp_source1,srf_speedrot,nx,ny,b_code)
	implicit none
!> @brief
!> rroe riemann solver in 2d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n,b_code
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft,cright
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real,intent(in)::nx,ny
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,el,er,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
real :: ul(4), ur(4)    !  input: conservative variables rho*[1, u, v, e]
 
 real :: rotated_rhll(4) ! output: rotated_rhll flux function.
!local constants
 
 real :: fifth, half, one, two    ! numbers
 real :: eps                            ! 
!local variables
 real :: nx1, ny1, nx2, ny2             ! rotated normals, n1 and n2
 real :: tx, ty                         ! tangent vector (taken as n1)
 real :: alpha1, alpha2                 ! projections of the new normals
 real :: vxl, vxr, vyl, vyr             ! velocity components.
 real :: rhol, rhor, pl, pr             ! primitive variables.
 real :: vnl, vnr, vtl, vtr             ! normal and tagent velocities
 real :: al, ar, hl, hr                 ! speeds of sound and total enthalpy
 real :: rt,rho,vx,vy,h,a               ! roe-averages
 real :: vn, vt                         ! normal and tagent velocities(roe-average)
 real :: drho,dvx,dvy,dvn,dvt,dp,dv(4)  ! wave strenghs
 real :: abs_dq                         ! magnitude of the velocity difference
 real :: abs_ws(4),ws(4),dws(4), rv(4,4)! wave speeds and right-eigevectors
 real :: srp,slm                        ! wave speeds for the hll part
 real :: fl(4), fr(4), diss(4)          ! fluxes ad dissipation term
 real :: temp
 integer ::  j

!constants.
!originally obtained from katate masatsuka, february 2009. http://www.cfdbooks.com
    
     fifth = 0.2
      half = 0.5
       one = 1.0
       two = 2.0
       eps = tolsmall ! 1.0e-12 in the original paper (double precision)
       
       ul(1:nof_variables)=cleft(1:nof_variables)
       ur(1:nof_variables)=cright(1:nof_variables)

!primitive and other variables.
!  left state
    rhol = ul(1)
     vxl = ul(2)/ul(1)
     vyl = ul(3)/ul(1)
      pl = (gamma-one)*( ul(4) - half*rhol*(vxl*vxl+vyl*vyl) )
      al = sqrt(gamma*pl/rhol)
      hl = ( ul(4) + pl ) / rhol
!  right state
    rhor = ur(1)
     vxr = ur(2)/ur(1)
     vyr = ur(3)/ur(1)
      pr = (gamma-one)*( ur(4) - half*rhor*(vxr*vxr+vyr*vyr) )
      ar = sqrt(gamma*pr/rhor)
      hr = ( ur(4) + pr ) / rhor

     vnl = vxl*nx + vyl*ny
     vnr = vxr*nx + vyr*ny

!compute the flux.
   fl(1) = rhol*vnl
   fl(2) = rhol*vnl * vxl + pl*nx
   fl(3) = rhol*vnl * vyl + pl*ny
   fl(4) = rhol*vnl *  hl

   fr(1) = rhor*vnr
   fr(2) = rhor*vnr * vxr + pr*nx
   fr(3) = rhor*vnr * vyr + pr*ny
   fr(4) = rhor*vnr *  hr

!define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(nb: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. for time-accurate 
!     calculation, this is fine.)
! nb: for a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    abs_dq = sqrt( (vxr-vxl)**2+(vyr-vyl)**2 )
  if ( abs_dq > eps) then
       nx1 = (vxr-vxl)/abs_dq
       ny1 = (vyr-vyl)/abs_dq
  else
    nx1 = -ny 
    ny1 =  nx
  endif
  
  
  
  
  
    alpha1 = nx * nx1 + ny * ny1 
!   to make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
    alpha1 = temp * alpha1

! take n2 as perpendicular to n1.
       nx2 = -ny1
       ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
!   to make alpha2 always positive.
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

!now we are going to compute the roe flux with n2 as the normal
!and n1 as the tagent vector, with modified wave speeds (5.12)

!compute the roe averages
     rt = sqrt(rhor/rhol)
    rho = rt*rhol
     vx = (vxl+rt*vxr)/(one+rt)
     vy = (vyl+rt*vyr)/(one+rt)
      h = ( hl+rt* hr)/(one+rt)
      a = sqrt( (gamma-one)*(h-half*(vx*vx+vy*vy)) )
     vn = vx*nx2+vy*ny2
     vt = vx*nx1+vy*ny1

!wave strengths (remember that n2 is the normal and n1 is the tangent.)
    vnl = vxl*nx2 + vyl*ny2
    vnr = vxr*nx2 + vyr*ny2
    vtl = vxl*nx1 + vyl*ny1
    vtr = vxr*nx1 + vyr*ny1

   drho = rhor - rhol 
     dp =   pr - pl
    dvn =  vnr - vnl
    dvt =  vtr - vtl

  dv(1) = (dp - rho*a*dvn )/(two*a*a)
  dv(2) =  rho*dvt/a
  dv(3) =  drho - dp/(a*a)
  dv(4) = (dp + rho*a*dvn )/(two*a*a)

!wave speeds for roe flux part.
    ws(1) = vn-a
    ws(2) = vn
    ws(3) = vn
    ws(4) = vn+a
  abs_ws  = abs(ws)

!harten's entropy fix jcp(1983), 49, pp357-393:
!only for the nonlinear fields.
  dws(1) = fifth
   if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  dws(4) = fifth
   if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

!hll wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   srp = max( zero, vtr + ar, vt + a)
   slm = min( zero, vtl - al, vt - a)

!modified wave speeds for the rotated-rhll flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(srp+slm)*ws + two*alpha1*srp*slm )/ (srp-slm)

!right eigenvectors: with n2 as normal and n1 as tangent.
  tx = nx1
  ty = ny1

  rv(1,1) = one    
  rv(2,1) = vx - a*nx2
  rv(3,1) = vy - a*ny2
  rv(4,1) =  h - vn*a

  rv(1,2) = zero
  rv(2,2) = a*tx
  rv(3,2) = a*ty
  rv(4,2) = a*vt

  rv(1,3) = one
  rv(2,3) = vx
  rv(3,3) = vy 
  rv(4,3) = half*(vx*vx+vy*vy)

  rv(1,4) = one
  rv(2,4) = vx + a*nx2
  rv(3,4) = vy + a*ny2
  rv(4,4) =  h + vn*a

!dissipation term: roe dissipation with the modified wave speeds.
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*dv(j)*rv(i,j)
   end do
  end do

!compute the rotated-rhll flux.
  hllcflux = (srp*fl - slm*fr)/(srp-slm) - half*diss

 





end subroutine rroe_riemann_solver2d





subroutine rusanov_riemann_solver2d(n,cleft_rot,cright_rot,hllcflux,mp_source1,srf_speedrot)
	implicit none
!> @brief
!> rusanov riemann solver in 2d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::cleft_rot,cright_rot
	real,dimension(1:nof_variables),intent(in)::srf_speedrot
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::rotvl,rotvr
	real,dimension(1:nof_variables)::leftv,rightv
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(inout)::hllcflux
	real,dimension(1)::sl,sr,sm
	real,intent(inout)::mp_source1
	real::mp_pinfl,mp_pinfr,gammal,gammar
	real::rl,rr,pl,pr,el,er,ul,ur,vl,vr,wl,wr,sped
	real::mul,mur,lastl,lastr,cc2,uu2,ccl,ccr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempfl,tempfr
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::flstar,frstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::ulstar,urstar
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::tempul,tempur
	real,dimension(1:nof_variables+turbulenceequations+passivescalar)::fl,fr
	real,dimension(turbulenceequations+passivescalar)::rml,rmr
	
	
	
	      
	hllcflux=zero
	rotvl=zero
	rotvr=zero
		tempfl=zero
		tempfr=zero
		flstar=zero
		frstar=zero
		ulstar=zero
		urstar=zero
		tempul=zero
		tempur=zero
		fl=zero
		fr=zero
		!conservative variables to primitive
		leftv(1:nof_variables)=cleft_rot(1:nof_variables)
		rightv(1:nof_variables)=cright_rot(1:nof_variables)
		
		call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
		
		rotvl(1:nof_variables)=leftv(1:nof_variables)
		rotvr(1:nof_variables)=rightv(1:nof_variables)
		
		
		
		
		if ((turbulence.eq.1).or.(passivescalar.gt.0))then
		
		rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cleft_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=cright_rot(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
		
		end if
		
		!call estimate_waves2d(n,rotvl,rotvr,sl,sm,sr,gamma)
			!now conditions based on wave speeds!
			rl=rotvl(1);ul=rotvl(2);vl=rotvl(3);pl=rotvl(4);el=cleft_rot(4)
			rr=rotvr(1);ur=rotvr(2);vr=rotvr(3);pr=rotvr(4);er=cright_rot(4)
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			
			rml(1:0+turbulenceequations+passivescalar)=rotvl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			rmr(1:0+turbulenceequations+passivescalar)=rotvr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
			

! 			if (turbulencemodel.eq.2)then
! 			pl=pl+((2.0d0/3.0d0)*eddyfl(2))
! 			pr=pr+((2.0d0/3.0d0)*eddyfr(2))  
! 
! 			end if
			end if


			fl(1)=rl*ul
			fl(2)=(rl*(ul**2))+pl
			fl(3)=rl*ul*vl
			
			fl(4)=ul*(el+pl)
			
			if ((multispecies.eq.1).or.(realgas.eq.1))then
			fl(dimensiona+3:nof_variables)=rotvl(dimensiona+3:nof_variables)*ul*rl
    
			end if
			
			
			fr(1)=rr*ur
			fr(2)=(rr*(ur**2))+pr
			fr(3)=rr*ur*vr
			
			fr(4)=ur*(er+pr)
			if ((multispecies.eq.1).or.(realgas.eq.1))then
			fr(dimensiona+3:nof_variables)=rotvr(dimensiona+3:nof_variables)*ur*rr
			end if
			
			if ((turbulence.eq.1).or.(passivescalar.gt.0))then
			fl(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rml(1:0+turbulenceequations+passivescalar)*ul
			fr(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=rmr(1:0+turbulenceequations+passivescalar)*ur
			end if
			
			
			
			
! 			
! 			flstar(:)=fl(:)+sl(1)*(ulstar(:)-cleft_rot(:))
! 			frstar(:)=fr(:)+sr(1)*(urstar(:)-cright_rot(:))
			
			
			if ((multispecies.eq.1).or.(realgas.eq.1))then
			sl(1)=abs(ul)+sqrt(gammal*(pl+mp_pinfl)/rl)
			sr(1)=abs(ur)+sqrt(gammar*(pr+mp_pinfr)/rr)
			mp_source1=0.5d0*(ul+ur)!-0.5d0*(max(abs(sl(1)),abs(sr(1))))*(ur-ul)
			else
			sl(1)=abs(ul)+sqrt(gamma*pl/rl)
			sr(1)=abs(ur)+sqrt(gamma*pr/rr)
			endif
! 			write(190+n,*)"why here"
! 			write(190+n,*)sl(1),sr(1)
! 			write(190+n,*)fl
! 			write(190+n,*)fr
			
			
			hllcflux(:)=0.5d0*(fl(:)+fr(:))-0.5d0*max(abs(sl(1)),abs(sr(1)))*(cright_rot(:)-cleft_rot(:))
			


end subroutine rusanov_riemann_solver2d




subroutine estimate_waves2d(n,rotvl,rotvr,sl,sm,sr,gamma)
	implicit none
!> @brief
!> waves speed estimates in 2d
#ifdef gpu
!$omp declare target
#endif
	integer,intent(in)::n
	integer::i,k
	real,dimension(1:nof_variables+turbulenceequations+passivescalar),intent(in)::rotvl,rotvr
	real,dimension(1),intent(inout)::sl,sr,sm
	real,intent(in)::gamma
	real::cl,cr,pr,pl,ul,ur,vl,vr,wr,wl,rl,rr
	real::cup,ppv,pmin,pmax,qmax,quser,bl,br,cov,pm,um
	real::g1,g2,g3,g4,g5,g6,g7,g8,gel,ger,pq,ptl,ptr
	g1 = (gamma - 1.0d0)/(2.0d0*gamma)
    	g2 = (gamma + 1.0d0)/(2.0d0*gamma)
   	g3 = 2.0d0*gamma/(gamma - 1.0d0)
    	g4 = 2.0d0/(gamma - 1.0d0)
    	g5 = 2.0d0/(gamma + 1.0d0)
    	g6 = (gamma - 1.0d0)/(gamma + 1.0d0)
   	g7 = (gamma - 1.0d0)/2.0d0
   	g8 = gamma - 1.0d0
	
	sl=0.0d0
	sr=0.0d0
	sm=0.0d0
	cov=0.0d0
	!build left state variables
	rl=rotvl(1)
	ul=rotvl(2)
	vl=rotvl(3)
	
	pl=rotvl(4)
	cl=sqrt((pl*gamma)/(rl))
	!build right  state variables
	rr=rotvr(1)
	ur=rotvr(2)
	vr=rotvr(3)
	
	pr=rotvr(4)
	cr=sqrt((pr*gamma)/(rr))

	cup=0.25d0*(rl+rr)*(cl+cr)
	ppv=0.5d0*(pl + pr) + 0.5d0*(ul - ur)*cup
	ppv=max(0.0d0,ppv)
	pmin=min(pl,pr)
	pmax=max(pl,pr)
	qmax=pmax/pmin
	quser=2.0d0

	 if(qmax.le.quser.and.(pmin.le.ppv.and.ppv.le.pmax))then
  
!        select prvs riemann solver
 
         pm = ppv
         um = 0.5d0*(ul + ur) + 0.5d0*(pl - pr)/cup 
  
      	else
 
         bl = 1.0d0 - cov*rl
         br = 1.0d0 - cov*rr
 
         if(ppv.lt.pmin)then
 
!           select two-rarefaction riemann solver
 
        
            pq  = (pl/pr)**g1
            um  = (pq*ul/cl/bl + ur/cr/br + g4*(pq - 1.0d0)) 
            um  = um/(pq/cl/bl + 1.0d0/cr/br)
            ptl = 1.0d0 + g7*(ul - um)/cl/bl
            ptr = 1.0d0 + g7*(um - ur)/cr/br
            pm  = 0.5d0*(pl*ptl**g3 + pr*ptr**g3)
         else

!           use two-shock riemann solver with pvrs as estimate
 
!          introduce iterations with pvrs as initial guess
 
            do k=1,20
             gel = sqrt((g5*bl/rl)/(g6*pl + ppv))
             ger = sqrt((g5*br/rr)/(g6*pr + ppv))
             pm  = (gel*pl + ger*pr - (ur - ul))/(gel + ger)
             um  = 0.5d0*(ul + ur) + 0.5d0*(ger*(pm - pr) - gel*(pm - pl))	     
             if ( abs((pm-ppv)/pm) .le. 1d-8) goto 101
                ppv = pm
     	    end do
         endif
      endif
      
      101 continue

!     find speeds
 
      if(pm.le.pl)then
         sl(1) = ul - cl
      else
         sl(1) = ul - cl*sqrt(1.0d0 + g2*(pm/pl - 1.0d0))
      endif
 
      sm(1)= um
 
      if(pm.le.pr)then
         sr(1)= ur + cr
      else
         sr(1) = ur + cr*sqrt(1.0d0 + g2*(pm/pr - 1.0d0))
      endif

end subroutine estimate_waves2d








end module riemann
