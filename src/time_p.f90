module advance
use library
use transform
use fluxes
use source
use initialisation
use boundary
use recon
use local
use flow_operations
use communications
use implicit_time
use implicit_fluxes
use declaration
use io
use moodr
implicit none

real,parameter :: fv_update_floor_fraction=1.0d-14
real,parameter :: fv_update_max_ref_ratio=1.0d8
real,parameter :: fv_update_min_ratio=1.0d-3
real,parameter :: fv_update_max_ratio=1.0d3
real,parameter :: fv_update_min_abs=1.0d-300
real,parameter :: fv_update_max_abs=1.0d300
real,parameter :: fv_update_species_closure_tol=1.0d-6
! Keep cutting the local FV update until an admissible nonzero step is found.
integer,parameter :: fv_update_max_backtracks=40
real,parameter :: realgas_source_dtl_cfl=2.0d-1
real,parameter :: realgas_source_species_floor_fraction=1.0d-10
real,parameter :: realgas_source_energy_floor_fraction=1.0d-10
real,parameter :: realgas_source_dtl_min_abs=1.0d-30

 contains

real function timestep_viscous_diffusivity(qcons,qprim,viscl,laml)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar),intent(in)::qcons,qprim
real,dimension(1:4),intent(in)::viscl,laml
integer::rg_i,idx
real::rho,cv_mix,y_i,rspec
real::mom_diff,thermal_diff,species_diff
real::mp_mu_mix,mp_ktr_mix,mp_kve,gammal
real,dimension(1:gpu_max_species)::mp_d_eff,mp_htr,mp_hvib

rho=max(qprim(1),tolsmall)
mom_diff=(4.0d0/3.0d0)*max(viscl(1),0.0d0)/rho
species_diff=0.0d0

if ((realgas.eq.1).and.(nof_species.gt.0))then
  cv_mix=0.0d0
  do rg_i=1,nof_species
    idx=dimensiona+3+rg_i
    if (idx.le.nof_variables)then
      y_i=max(qprim(idx),0.0d0)
      rspec=rgs_ru/max(abs(rg_molm(rg_i)),tolsmall)
      if (rg_i.le.3)then
        cv_mix=cv_mix+y_i*2.5d0*rspec
      else
        cv_mix=cv_mix+y_i*1.5d0*rspec
      end if
    end if
  end do
  thermal_diff=max(laml(1),0.0d0)/(rho*max(cv_mix,tolsmall))

  mp_mu_mix=0.0d0
  mp_ktr_mix=0.0d0
  mp_kve=0.0d0
  gammal=gamma
  mp_d_eff=0.0d0
  mp_htr=0.0d0
  mp_hvib=0.0d0
  call multispecies_mixtures_rg(qcons,mp_mu_mix,mp_ktr_mix,mp_kve,mp_d_eff,mp_htr,mp_hvib,gammal)
  do rg_i=1,nof_species
    if ((mp_d_eff(rg_i).eq.mp_d_eff(rg_i)).and.(mp_d_eff(rg_i).gt.0.0d0))then
      species_diff=max(species_diff,mp_d_eff(rg_i))
    end if
  end do
else
  thermal_diff=max(laml(1),0.0d0)*(gamma-1.0d0)/(rho*max(r_gas,tolsmall))
end if

timestep_viscous_diffusivity=max(mom_diff,thermal_diff,species_diff)

end function timestep_viscous_diffusivity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!subroutine called to calculate cfl number depending on the scheme!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!employed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calculate_cfl(n)
!> @brief
!> subroutine for computing the global time step size in 3d
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::veln,agrt,diff_scale
real,dimension(1:gpu_max_nvar)::leftv,rightv
real,dimension(4)::srf_speed
real,dimension(3)::rotvec
real::mp_pinfl,gammal
real,dimension(3)::pox,poy
real,dimension(1:4)::viscl,laml
	real,dimension(1:20)::eddyfl,eddyfr
	real,dimension(1:2)::turbmv
	real,dimension(1)::etvm
#ifdef xpu
	real::cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_w,cfl_q2,cfl_p,cfl_temp,cfl_tref,cfl_tratio
	real::cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1
	real::cfl_rx,cfl_ry,cfl_rz,cfl_sx,cfl_sy,cfl_sz,cfl_vdiff
#endif
#if defined(gpu) || defined(xpu)
	real::dt_reduce
#define CFL_DT dt_reduce
#else
#define CFL_DT dt
#endif

kmaxe=xmpielrank(n)

#if defined(gpu) || defined(xpu)
  dt_reduce = tolbig
#else
  dt = tolbig
#endif
	if (itestcase.lt.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) private(veln) &
!$omp& map(alloc: ielem_minedge) &
!$omp& firstprivate(ccfl, dg, iorder, lamx, lamy, lamz)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe
		veln=max(abs(lamx),abs(lamy),abs(lamz))

		if (dg.eq.1)then

		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))

		else

		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln))))
		end if

	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if

	if (itestcase.eq.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) private(veln,leftv,mp_pinfl,gammal,agrt,pox,poy,rotvec,srf_speed) &
!$omp& map(alloc: ielem_minedge, ielem_xxc, ielem_yyc, ielem_zzc, rec_mrf, rec_mrf_origin, rec_mrf_velocity) &
!$omp& map(alloc: srf_velocity, u_c_val) &
!$omp& firstprivate(ccfl, dg, gamma, iorder, mrf, multispecies, nof_variables, realgas, rframe, srfg, zero)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe



		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
#ifdef xpu
		call cons2prim_ideal(n,leftv,mp_pinfl,gammal)
#else
		call cons2prim(n,leftv,mp_pinfl,gammal)
#endif
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(5)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(5)*gamma/leftv(1))
		end if
        if (rframe.eq.0) then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
        end if
        if (srfg.eq.1) then
            pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
            poy(1:3)=srf_velocity
            srf_speed=zero
            call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
            if (rec_mrf(i).eq.0)then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox=pox-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
		if (dg.eq.1)then
		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))

		else
		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln))))

		end if


	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if



	if (itestcase.eq.4)then
#ifdef xpu
	cfl_tref=pres/(rres*r_gas)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) &
!$omp& firstprivate(kmaxe,dg,iorder,ccfl,gamma,prandtl,r_gas,visc,suther,turbulence,turbulencemodel,cv1,prtu,cfl_tref,srfg,mrf) &
!$omp& map(alloc: u_c_val,ielem_minedge,ielem_xxc,ielem_yyc,ielem_zzc,srf_velocity) &
!$omp& private(veln,agrt,cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_w,cfl_q2,cfl_p,cfl_temp,cfl_tratio) &
!$omp& private(cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1,cfl_rx,cfl_ry,cfl_rz,cfl_sx,cfl_sy,cfl_sz,cfl_vdiff)
#elif defined(gpu)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) &
!$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,dg,iorder,ccfl,gamma,prandtl) &
!$omp& firstprivate(r_gas,pres,rres,visc,suther,zero,rframe,srfg,mrf,turbulence,turbulencemodel) &
!$omp& map(alloc: xmpielrank,u_c_val,ielem_minedge,ielem_xxc,ielem_yyc,ielem_zzc) &
!$omp& map(alloc: srf_velocity) &
!$omp& private(veln,leftv,rightv,mp_pinfl,gammal,agrt,pox,poy,rotvec,srf_speed,viscl,laml,turbmv,etvm,eddyfl,eddyfr,diff_scale)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe
#ifdef xpu
		cfl_rho=u_c_val(1,1,i)
		cfl_orho=1.0d0/cfl_rho
		cfl_u=u_c_val(1,2,i)*cfl_orho
		cfl_v=u_c_val(1,3,i)*cfl_orho
		cfl_w=u_c_val(1,4,i)*cfl_orho
		cfl_q2=(cfl_u*cfl_u)+(cfl_v*cfl_v)+(cfl_w*cfl_w)
		cfl_p=(gamma-1.0d0)*(u_c_val(1,5,i)-0.5d0*cfl_rho*cfl_q2)
		agrt=sqrt(cfl_p*gamma*cfl_orho)
		cfl_temp=cfl_p*cfl_orho/r_gas
		cfl_tratio=cfl_temp/cfl_tref
		cfl_mu=visc*(cfl_tratio*sqrt(cfl_tratio))*((cfl_tref+suther*cfl_tref)/(cfl_temp+suther*cfl_tref))
		cfl_lam=cfl_mu*r_gas*gamma/(prandtl*(gamma-1.0d0))
		veln=max(abs(cfl_u),abs(cfl_v),abs(cfl_w))+agrt
		if (srfg.eq.1)then
			cfl_sx=(srf_velocity(2)*ielem_zzc(i))-(srf_velocity(3)*ielem_yyc(i))
			cfl_sy=(srf_velocity(3)*ielem_xxc(i))-(srf_velocity(1)*ielem_zzc(i))
			cfl_sz=(srf_velocity(1)*ielem_yyc(i))-(srf_velocity(2)*ielem_xxc(i))
			veln=max(abs(cfl_u-cfl_sx),abs(cfl_v-cfl_sy),abs(cfl_w-cfl_sz))+agrt
		end if
		if (mrf.eq.1)then
			if (rec_mrf(i).eq.0)then
				veln=max(abs(cfl_u),abs(cfl_v),abs(cfl_w))+agrt
			else
				cfl_rx=ielem_xxc(i)-rec_mrf_origin(1,i)
				cfl_ry=ielem_yyc(i)-rec_mrf_origin(2,i)
				cfl_rz=ielem_zzc(i)-rec_mrf_origin(3,i)
				cfl_sx=(rec_mrf_velocity(2,i)*cfl_rz)-(rec_mrf_velocity(3,i)*cfl_ry)
				cfl_sy=(rec_mrf_velocity(3,i)*cfl_rx)-(rec_mrf_velocity(1,i)*cfl_rz)
				cfl_sz=(rec_mrf_velocity(1,i)*cfl_ry)-(rec_mrf_velocity(2,i)*cfl_rx)
				veln=max(abs(cfl_u-cfl_sx),abs(cfl_v-cfl_sy),abs(cfl_w-cfl_sz))+agrt
			end if
		end if
		if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
			cfl_mut=u_ct_val(1,1,i)
			cfl_chi=abs(max(cfl_mut,1.0d-15)/max(cfl_mu,1.0d-15))
			cfl_chi3=cfl_chi*cfl_chi*cfl_chi
			cfl_fv1=cfl_chi3/(cfl_chi3+(cv1*cv1*cv1))
			cfl_mut=cfl_mut*cfl_fv1
			cfl_mut=min(10000000.0d0*visc,cfl_mut)
			cfl_lam=cfl_lam+(cfl_mut*r_gas*gamma/(prtu*(gamma-1.0d0)))
			cfl_mut=max(0.0d0,cfl_mut)
			if (u_ct_val(1,1,i).lt.0.0d0)cfl_mut=0.0d0
			cfl_mu=cfl_mu+cfl_mut
		end if
		cfl_vdiff=max(((4.0d0/3.0d0)*cfl_mu*cfl_orho),cfl_lam*(gamma-1.0d0)*cfl_orho/r_gas)
		if (dg.eq.1)then
			CFL_DT=min(CFL_DT,(ccfl/(2*iorder+1))*(ielem_minedge(i)/(abs(veln)+(2.0d0*cfl_vdiff*((2*iorder+1)/ielem_minedge(i))))))
		else
			CFL_DT=min(CFL_DT,ccfl*(1.0d0/((abs(veln)/ielem_minedge(i))+(0.5d0*cfl_vdiff/(ielem_minedge(i)**2)))))
		end if
#else
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		rightv(1:nof_variables)=leftv(1:nof_variables)
		call get_visc_conduct(n,leftv,rightv,viscl,laml)

        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		call cons2prim(n,leftv,mp_pinfl,gammal)
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(5)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(5)*gamma/leftv(1))
		end if




        if (rframe.eq.0) then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
        end if
        if(srfg.eq.1)then
            pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
            poy(1:3)=srf_velocity
            srf_speed=zero
            call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
            if (rec_mrf(i).eq.0)then
                veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox(1:3)=pox(1:3)-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
		if (turbulence.eq.1)then
		if (turbulencemodel.eq.1)then
		turbmv(1)=u_ct_val(1,1,i);  turbmv(2)=u_ct_val(1,1,i);
		eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if

		diff_scale=timestep_viscous_diffusivity(rightv,leftv,viscl,laml)

		if (dg.eq.1)then


		CFL_DT=min(CFL_DT,(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*diff_scale*((2*iorder+1)/ielem_minedge(i))))))


		else

         CFL_DT=min(CFL_DT,ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*diff_scale/((ielem_minedge(i)))**2))))



         end if
#endif

	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if


#if defined(gpu) || defined(xpu)
        dt = dt_reduce
#endif
#undef CFL_DT
        return

end subroutine calculate_cfl

subroutine calculate_cfll(n)
!> @brief
!> subroutine for computing the time step size for each cell in 3d
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::veln,agrt,diff_scale
real,dimension(1:gpu_max_nvar)::leftv,rightv
real,dimension(4)::srf_speed
real,dimension(3)::rotvec
real::mp_pinfl,gammal
real,dimension(3)::pox,poy
real,dimension(1:4)::viscl,laml
	real,dimension(1:20)::eddyfl,eddyfr
	real,dimension(1:2)::turbmv
	real,dimension(1)::etvm
#ifdef xpu
	real::cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_w,cfl_q2,cfl_p,cfl_temp,cfl_tref,cfl_tratio
	real::cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1
	real::cfl_rx,cfl_ry,cfl_rz,cfl_sx,cfl_sy,cfl_sz,cfl_vdiff
#endif
	kmaxe=xmpielrank(n)




	if (itestcase.lt.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private (veln) &
!$omp& map(alloc: ielem_dtl, ielem_minedge) &
!$omp& firstprivate(ccfl, lamx, lamy, lamz)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		veln=max(abs(lamx),abs(lamy),abs(lamz))
		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))
	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if

	if (itestcase.eq.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(veln,leftv,mp_pinfl,gammal,agrt,pox,poy,rotvec,srf_speed) &
!$omp& map(alloc: ielem_dtl, ielem_minedge, ielem_xxc, ielem_yyc, ielem_zzc, rec_mrf, rec_mrf_origin) &
!$omp& map(alloc: rec_mrf_velocity, srf_velocity, u_c_val) &
!$omp& firstprivate(ccfl, dg, gamma, iorder, mrf, multispecies, nof_variables, realgas, rframe, srfg, zero)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
#ifdef xpu
		call cons2prim_ideal(n,leftv,mp_pinfl,gammal)
#else
		call cons2prim(n,leftv,mp_pinfl,gammal)
#endif

       if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(5)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(5)*gamma/leftv(1))
		end if
        if (rframe.eq.0) then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
        end if
        if (srfg.eq.1) then
            pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
            poy(1:3)=srf_velocity
            srf_speed=zero
            call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
            if (rec_mrf(i).eq.1)then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox=pox-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
		if (dg.eq.1)then

		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1))

		else

		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))

		end if
	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if



	if (itestcase.eq.4)then
#ifdef xpu
	cfl_tref=pres/(rres*r_gas)
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe,dg,iorder,ccfl,gamma,prandtl,r_gas,visc,suther,turbulence,turbulencemodel,cv1,prtu,cfl_tref,srfg,mrf) &
!$omp& map(alloc: u_c_val,ielem_dtl,ielem_minedge,ielem_xxc,ielem_yyc,ielem_zzc,srf_velocity) &
!$omp& private(veln,agrt,cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_w,cfl_q2,cfl_p,cfl_temp,cfl_tratio) &
!$omp& private(cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1,cfl_rx,cfl_ry,cfl_rz,cfl_sx,cfl_sy,cfl_sz,cfl_vdiff)
#elif defined(gpu)
!$omp target teams distribute parallel do &
!$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,dg,iorder,ccfl,gamma,prandtl) &
!$omp& firstprivate(r_gas,pres,rres,visc,suther,zero,rframe,srfg,mrf,turbulence,turbulencemodel) &
!$omp& map(alloc: xmpielrank,u_c_val,ielem_dtl,ielem_minedge,ielem_xxc,ielem_yyc,ielem_zzc) &
!$omp& map(alloc: srf_velocity) &
!$omp& private(veln,leftv,rightv,mp_pinfl,gammal,agrt,pox,poy,rotvec,srf_speed,viscl,laml,turbmv,etvm,eddyfl,eddyfr,diff_scale)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe

#ifdef xpu
		cfl_rho=u_c_val(1,1,i)
		cfl_orho=1.0d0/cfl_rho
		cfl_u=u_c_val(1,2,i)*cfl_orho
		cfl_v=u_c_val(1,3,i)*cfl_orho
		cfl_w=u_c_val(1,4,i)*cfl_orho
		cfl_q2=(cfl_u*cfl_u)+(cfl_v*cfl_v)+(cfl_w*cfl_w)
		cfl_p=(gamma-1.0d0)*(u_c_val(1,5,i)-0.5d0*cfl_rho*cfl_q2)
		agrt=sqrt(cfl_p*gamma*cfl_orho)
		cfl_temp=cfl_p*cfl_orho/r_gas
		cfl_tratio=cfl_temp/cfl_tref
		cfl_mu=visc*(cfl_tratio*sqrt(cfl_tratio))*((cfl_tref+suther*cfl_tref)/(cfl_temp+suther*cfl_tref))
		cfl_lam=cfl_mu*r_gas*gamma/(prandtl*(gamma-1.0d0))
		veln=max(abs(cfl_u),abs(cfl_v),abs(cfl_w))+agrt
		if (srfg.eq.1)then
			cfl_sx=(srf_velocity(2)*ielem_zzc(i))-(srf_velocity(3)*ielem_yyc(i))
			cfl_sy=(srf_velocity(3)*ielem_xxc(i))-(srf_velocity(1)*ielem_zzc(i))
			cfl_sz=(srf_velocity(1)*ielem_yyc(i))-(srf_velocity(2)*ielem_xxc(i))
			veln=max(abs(cfl_u-cfl_sx),abs(cfl_v-cfl_sy),abs(cfl_w-cfl_sz))+agrt
		end if
		if (mrf.eq.1)then
			if (rec_mrf(i).eq.0)then
				veln=max(abs(cfl_u),abs(cfl_v),abs(cfl_w))+agrt
			else
				cfl_rx=ielem_xxc(i)-rec_mrf_origin(1,i)
				cfl_ry=ielem_yyc(i)-rec_mrf_origin(2,i)
				cfl_rz=ielem_zzc(i)-rec_mrf_origin(3,i)
				cfl_sx=(rec_mrf_velocity(2,i)*cfl_rz)-(rec_mrf_velocity(3,i)*cfl_ry)
				cfl_sy=(rec_mrf_velocity(3,i)*cfl_rx)-(rec_mrf_velocity(1,i)*cfl_rz)
				cfl_sz=(rec_mrf_velocity(1,i)*cfl_ry)-(rec_mrf_velocity(2,i)*cfl_rx)
				veln=max(abs(cfl_u-cfl_sx),abs(cfl_v-cfl_sy),abs(cfl_w-cfl_sz))+agrt
			end if
		end if
		if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
			cfl_mut=u_ct_val(1,1,i)
			cfl_chi=abs(max(cfl_mut,1.0d-15)/max(cfl_mu,1.0d-15))
			cfl_chi3=cfl_chi*cfl_chi*cfl_chi
			cfl_fv1=cfl_chi3/(cfl_chi3+(cv1*cv1*cv1))
			cfl_mut=cfl_mut*cfl_fv1
			cfl_mut=min(10000000.0d0*visc,cfl_mut)
			cfl_lam=cfl_lam+(cfl_mut*r_gas*gamma/(prtu*(gamma-1.0d0)))
			cfl_mut=max(0.0d0,cfl_mut)
			if (u_ct_val(1,1,i).lt.0.0d0)cfl_mut=0.0d0
			cfl_mu=cfl_mu+cfl_mut
		end if
		cfl_vdiff=max(((4.0d0/3.0d0)*cfl_mu*cfl_orho),cfl_lam*(gamma-1.0d0)*cfl_orho/r_gas)
		if (dg.eq.1)then
			ielem_dtl(i)=(ccfl/(2*iorder+1))*(ielem_minedge(i)/(abs(veln)+(2.0d0*cfl_vdiff*((2*iorder+1)/ielem_minedge(i)))))
		else
			ielem_dtl(i)=ccfl*(1.0d0/((abs(veln)/ielem_minedge(i))+(0.5d0*cfl_vdiff/(ielem_minedge(i)**2))))
		end if
#else
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		rightv(1:nof_variables)=leftv(1:nof_variables)
		call get_visc_conduct(n,leftv,rightv,viscl,laml)

        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		call cons2prim(n,leftv,mp_pinfl,gammal)
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(5)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(5)*gamma/leftv(1))
		end if








        if (srfg.eq.0.and.mrf.eq.0) then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
        end if
        if(srfg.eq.1)then
            pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
            poy(1:3)=srf_velocity
            srf_speed=zero
            call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
           if (rec_mrf(i).eq.0)then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox(1:3)=pox(1:3)-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                call vect_function(pox,poy,rotvec)
            srf_speed(2)=rotvec(1)
            srf_speed(3)=rotvec(2)
            srf_speed(4)=rotvec(3)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
        if (turbulence.eq.1)then
		if (turbulencemodel.eq.1)then
		turbmv(1)=u_ct_val(1,1,i);  turbmv(2)=u_ct_val(1,1,i);
		eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		call eddyvisco_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if


		diff_scale=timestep_viscous_diffusivity(rightv,leftv,viscl,laml)

		if (dg.eq.1)then


		ielem_dtl(i)=(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*diff_scale*((2*iorder+1)/ielem_minedge(i)))))


		else

		ielem_dtl(i)=ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*diff_scale/((ielem_minedge(i)))**2)))
		end if
#endif

	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if


        call limit_realgas_source_timestep(n)

        return

end subroutine calculate_cfll



subroutine calculate_cfl2d(n)
!> @brief
!> subroutine for computing the global time step size in 2d
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::veln,agrt,lamxl,lamyl,diff_scale
real,dimension(1:gpu_max_nvar)::leftv,rightv
real::mp_pinfl,gammal
real,dimension(1:4)::viscl,laml
	real,dimension(1:20)::eddyfl,eddyfr
	real,dimension(1:2)::turbmv
	real,dimension(1)::etvm
#ifdef xpu
	real::cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_q2,cfl_p,cfl_temp,cfl_tref,cfl_tratio
	real::cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1,cfl_vdiff
#endif
#if defined(gpu) || defined(xpu)
	real::dt_reduce
#define CFL_DT dt_reduce
#else
#define CFL_DT dt
#endif
kmaxe=xmpielrank(n)



#if defined(gpu) || defined(xpu)
  dt_reduce = tolbig
#else
  dt = tolbig
#endif




	if (itestcase.lt.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) private(veln,lamxl,lamyl) &
!$omp& map(alloc: ielem_minedge, ielem_xxc, ielem_yyc) &
!$omp& firstprivate(ccfl, dg, initcond, iorder, lamx, lamy)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe

				if (initcond.eq.3)then
				  lamxl=-ielem_yyc(i)+0.5d0
				  lamyl=ielem_xxc(i)-0.5
				  else
				  lamxl=lamx
				  lamyl=lamy
                end if


		veln=max(abs(lamxl),abs(lamyl))

		if (dg.eq.1)then

		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))

		else

		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln))))
		end if

	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if

	if (itestcase.eq.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) private(veln,leftv,mp_pinfl,gammal,agrt) &
!$omp& map(alloc: ielem_minedge, u_c_val) &
!$omp& firstprivate(ccfl, dg, gamma, iorder, multispecies, nof_variables, realgas)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe


		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
#ifdef xpu
		call cons2prim_ideal(n,leftv,mp_pinfl,gammal)
#else
		call cons2prim(n,leftv,mp_pinfl,gammal)
#endif
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(4)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(4)*gamma/leftv(1))
		end if
		veln=max(abs(leftv(2)),abs(leftv(3)))+agrt



		if (dg.eq.1)then

		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))


		else

		CFL_DT=min(CFL_DT,ccfl*((ielem_minedge(i))/(abs(veln))))
		end if



	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if



	if (itestcase.eq.4)then
#ifdef xpu
	cfl_tref=pres/(rres*r_gas)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) &
!$omp& firstprivate(kmaxe,dg,iorder,ccfl,gamma,prandtl,r_gas,visc,suther,turbulence,turbulencemodel,cv1,prtu,cfl_tref) &
!$omp& map(alloc: u_c_val,ielem_minedge) &
!$omp& private(veln,agrt,cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_q2,cfl_p,cfl_temp,cfl_tratio) &
!$omp& private(cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1,cfl_vdiff)
#elif defined(gpu)
!$omp target teams distribute parallel do &
!$omp& map(tofrom: dt_reduce) &
!$omp& reduction(min: dt_reduce) &
!$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,dg,iorder,ccfl,gamma,prandtl) &
!$omp& firstprivate(r_gas,pres,rres,visc,suther,zero,turbulence,turbulencemodel) &
!$omp& map(alloc: xmpielrank,u_c_val,ielem_minedge) &
!$omp& private(veln,leftv,rightv,mp_pinfl,gammal,agrt,viscl,laml,turbmv,etvm,eddyfl,eddyfr,diff_scale)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe

#ifdef xpu
		cfl_rho=u_c_val(1,1,i)
		cfl_orho=1.0d0/cfl_rho
		cfl_u=u_c_val(1,2,i)*cfl_orho
		cfl_v=u_c_val(1,3,i)*cfl_orho
		cfl_q2=(cfl_u*cfl_u)+(cfl_v*cfl_v)
		cfl_p=(gamma-1.0d0)*(u_c_val(1,4,i)-0.5d0*cfl_rho*cfl_q2)
		agrt=sqrt(cfl_p*gamma*cfl_orho)
		cfl_temp=cfl_p*cfl_orho/r_gas
		cfl_tratio=cfl_temp/cfl_tref
		cfl_mu=visc*(cfl_tratio*sqrt(cfl_tratio))*((cfl_tref+suther*cfl_tref)/(cfl_temp+suther*cfl_tref))
		cfl_lam=cfl_mu*r_gas*gamma/(prandtl*(gamma-1.0d0))
		veln=max(abs(cfl_u),abs(cfl_v))+agrt
		if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
			cfl_mut=u_ct_val(1,1,i)
			cfl_chi=abs(max(cfl_mut,1.0d-15)/max(cfl_mu,1.0d-15))
			cfl_chi3=cfl_chi*cfl_chi*cfl_chi
			cfl_fv1=cfl_chi3/(cfl_chi3+(cv1*cv1*cv1))
			cfl_mut=cfl_mut*cfl_fv1
			cfl_mut=min(10000000.0d0*visc,cfl_mut)
			cfl_lam=cfl_lam+(cfl_mut*r_gas*gamma/(prtu*(gamma-1.0d0)))
			cfl_mut=max(0.0d0,cfl_mut)
			if (u_ct_val(1,1,i).lt.0.0d0)cfl_mut=0.0d0
			cfl_mu=cfl_mu+cfl_mut
		end if
		cfl_vdiff=max(((4.0d0/3.0d0)*cfl_mu*cfl_orho),cfl_lam*(gamma-1.0d0)*cfl_orho/r_gas)
		if (dg.eq.1)then
			CFL_DT=min(CFL_DT,(ccfl/(2*iorder+1))*(ielem_minedge(i)/(abs(veln)+(2.0d0*cfl_vdiff*((2*iorder+1)/ielem_minedge(i))))))
		else
			CFL_DT=min(CFL_DT,ccfl*(1.0d0/((abs(veln)/ielem_minedge(i))+(0.5d0*cfl_vdiff/(ielem_minedge(i)**2)))))
		end if
#else
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		rightv(1:nof_variables)=leftv(1:nof_variables)
		call get_visc_conduct(n,leftv,rightv,viscl,laml)

        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		call cons2prim(n,leftv,mp_pinfl,gammal)
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(4)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(4)*gamma/leftv(1))
		end if


		veln=max(abs(leftv(2)),abs(leftv(3)))+agrt
		if (turbulence.eq.1)then
		if (turbulencemodel.eq.1)then
		turbmv(1)=u_ct_val(1,1,i);  turbmv(2)=u_ct_val(1,1,i);
		eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if


    diff_scale=timestep_viscous_diffusivity(rightv,leftv,viscl,laml)

    if (dg.eq.1)then


      CFL_DT=min(CFL_DT,(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*diff_scale*((2*iorder+1)/ielem_minedge(i))))))


      else



			CFL_DT=min(CFL_DT,ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*diff_scale/((ielem_minedge(i)))**2))))
	      end if
#endif

	  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if


#if defined(gpu) || defined(xpu)
        dt = dt_reduce
#endif
#undef CFL_DT
        return

end subroutine calculate_cfl2d


subroutine calculate_cfll2d(n)
!> @brief
!> subroutine for computing the time step size for each cell in 2d
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::veln,agrt,diff_scale
real,dimension(1:gpu_max_nvar)::leftv,rightv
real::mp_pinfl,gammal
real,dimension(1:4)::viscl,laml
	real,dimension(1:20)::eddyfl,eddyfr
	real,dimension(1:2)::turbmv
	real,dimension(1)::etvm
#ifdef xpu
	real::cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_q2,cfl_p,cfl_temp,cfl_tref,cfl_tratio
	real::cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1,cfl_vdiff
#endif
	kmaxe=xmpielrank(n)




	if (itestcase.lt.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private (veln) &
!$omp& map(alloc: ielem_dtl, ielem_minedge) &
!$omp& firstprivate(ccfl, lamx, lamy)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		veln=max(abs(lamx),abs(lamy))
		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))
	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if

	if (itestcase.eq.3)then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(veln,leftv,mp_pinfl,gammal,agrt) &
!$omp& map(alloc: ielem_dtl, ielem_minedge, u_c_val) &
!$omp& firstprivate(ccfl, dg, gamma, iorder, multispecies, nof_variables, realgas)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
#ifdef xpu
		call cons2prim_ideal(n,leftv,mp_pinfl,gammal)
#else
		call cons2prim(n,leftv,mp_pinfl,gammal)
#endif
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(4)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(4)*gamma/leftv(1))
		end if



	        !	agrt=sqrt(leftv(4)*gamma/leftv(1))
		veln=max(abs(leftv(2)),abs(leftv(3)))+agrt


		if (dg.eq.1)then

		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1))

		else

		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))

		end if




	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if



	if (itestcase.eq.4)then
#ifdef xpu
	cfl_tref=pres/(rres*r_gas)
!$omp target teams distribute parallel do &
!$omp& firstprivate(kmaxe,dg,iorder,ccfl,gamma,prandtl,r_gas,visc,suther,turbulence,turbulencemodel,cv1,prtu,cfl_tref) &
!$omp& map(alloc: u_c_val,ielem_dtl,ielem_minedge) &
!$omp& private(veln,agrt,cfl_rho,cfl_orho,cfl_u,cfl_v,cfl_q2,cfl_p,cfl_temp,cfl_tratio) &
!$omp& private(cfl_mu,cfl_lam,cfl_mut,cfl_chi,cfl_chi3,cfl_fv1,cfl_vdiff)
#elif defined(gpu)
!$omp target teams distribute parallel do &
!$omp& firstprivate(n,kmaxe,nof_variables,dimensiona,dg,iorder,ccfl,gamma,prandtl) &
!$omp& firstprivate(r_gas,pres,rres,visc,suther,zero,turbulence,turbulencemodel) &
!$omp& map(alloc: xmpielrank,u_c_val,ielem_dtl,ielem_minedge) &
!$omp& private(veln,leftv,rightv,mp_pinfl,gammal,agrt,viscl,laml,turbmv,etvm,eddyfl,eddyfr,diff_scale)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
#ifdef xpu
		cfl_rho=u_c_val(1,1,i)
		cfl_orho=1.0d0/cfl_rho
		cfl_u=u_c_val(1,2,i)*cfl_orho
		cfl_v=u_c_val(1,3,i)*cfl_orho
		cfl_q2=(cfl_u*cfl_u)+(cfl_v*cfl_v)
		cfl_p=(gamma-1.0d0)*(u_c_val(1,4,i)-0.5d0*cfl_rho*cfl_q2)
		agrt=sqrt(cfl_p*gamma*cfl_orho)
		cfl_temp=cfl_p*cfl_orho/r_gas
		cfl_tratio=cfl_temp/cfl_tref
		cfl_mu=visc*(cfl_tratio*sqrt(cfl_tratio))*((cfl_tref+suther*cfl_tref)/(cfl_temp+suther*cfl_tref))
		cfl_lam=cfl_mu*r_gas*gamma/(prandtl*(gamma-1.0d0))
		veln=max(abs(cfl_u),abs(cfl_v))+agrt
		if ((turbulence.eq.1).and.(turbulencemodel.eq.1))then
			cfl_mut=u_ct_val(1,1,i)
			cfl_chi=abs(max(cfl_mut,1.0d-15)/max(cfl_mu,1.0d-15))
			cfl_chi3=cfl_chi*cfl_chi*cfl_chi
			cfl_fv1=cfl_chi3/(cfl_chi3+(cv1*cv1*cv1))
			cfl_mut=cfl_mut*cfl_fv1
			cfl_mut=min(10000000.0d0*visc,cfl_mut)
			cfl_lam=cfl_lam+(cfl_mut*r_gas*gamma/(prtu*(gamma-1.0d0)))
			cfl_mut=max(0.0d0,cfl_mut)
			if (u_ct_val(1,1,i).lt.0.0d0)cfl_mut=0.0d0
			cfl_mu=cfl_mu+cfl_mut
		end if
		cfl_vdiff=max(((4.0d0/3.0d0)*cfl_mu*cfl_orho),cfl_lam*(gamma-1.0d0)*cfl_orho/r_gas)
		if (dg.eq.1)then
			ielem_dtl(i)=(ccfl/(2*iorder+1))*(ielem_minedge(i)/(abs(veln)+(2.0d0*cfl_vdiff*((2*iorder+1)/ielem_minedge(i)))))
		else
			ielem_dtl(i)=ccfl*(1.0d0/((abs(veln)/ielem_minedge(i))+(0.5d0*cfl_vdiff/(ielem_minedge(i)**2))))
		end if
#else
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		rightv(1:nof_variables)=leftv(1:nof_variables)
		call get_visc_conduct(n,leftv,rightv,viscl,laml)

        leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		call cons2prim(n,leftv,mp_pinfl,gammal)
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(4)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(4)*gamma/leftv(1))
		end if










		veln=max(abs(leftv(2)),abs(leftv(3)))+agrt
		if (turbulence.eq.1)then
		if (turbulencemodel.eq.1)then
		turbmv(1)=u_ct_val(1,1,i);  turbmv(2)=u_ct_val(1,1,i);
		eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		call eddyvisco2d_ideal(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if






		diff_scale=timestep_viscous_diffusivity(rightv,leftv,viscl,laml)

		if (dg.eq.1)then


      ielem_dtl(i)=(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*diff_scale*((2*iorder+1)/ielem_minedge(i)))))


      else

			ielem_dtl(i)=ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*diff_scale/((ielem_minedge(i)))**2)))


			end if
#endif
	end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if


        call limit_realgas_source_timestep(n)

        return

end subroutine calculate_cfll2d


! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !---------------------------------------------------------------------------------------------!
! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!subroutine called to advance solution by one time step size!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine runge_kutta3_mood(n)
!> @brief
!> ssp runge kutta 3rd-order scheme
implicit none
integer,intent(in)::n
integer::iavr,nvar,i,kmaxe,inds
real::avrgs,oovolume
kmaxe=xmpielrank(n)



! if (mood.eq.1)then
! inds=4
! else
! inds=1
! end if

if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if

    end select
else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if

    end select
end if


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

 if (mood.eq.1)then

 call mood_operator_2(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(4,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
    ielem_mood_o(i)=2
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
    ielem_mood_o(i)=1
   else
   u_c_val(1,1:nof_variables,i)=u_c_val(4,1:nof_variables,i)
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if
    end select
else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if
    end select
end if

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo4, to4)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo4, passivescalar, to4, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

 if (mood.eq.1)then

 call mood_operator_2(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo4, to4)
#else
	!$omp do
#endif
do i=1,kmaxe
     if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(4,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=2
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo4, to4)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=1
   else
   u_c_val(1,1:nof_variables,i)=u_c_val(4,1:nof_variables,i)
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if (turbulence.eq.1)then
    call sources_computation(n)
    end if
    end select
else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if (turbulence.eq.1)then
    call sources_computation(n)
    end if
    call vortexcalc(n)
    end select
end if
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo3, to3)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if (mood.eq.1)then

 call mood_operator_2(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo3, to3)
#else
	!$omp do
#endif
do i=1,kmaxe
     if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(4,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(3,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=2
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

call mood_operator_1(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo3, to3)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=1
   else
   u_c_val(1,1:nof_variables,i)=u_c_val(4,1:nof_variables,i)
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo3, passivescalar, to3, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta3_mood


subroutine runge_kutta3(n)
!> @brief
!> ssp runge kutta 3rd-order scheme
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
logical::accepted_update
real::avrgs,oovolume,update_alpha
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
kmaxe=xmpielrank(n)




call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)

           call dg_mass_rhs_product_cell(i)
           do mm_var=1,nof_variables
             do mm_idof=1,num_dg_dofs
               u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
             end do
           end do




    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        old_state=zero
        candidate=zero
        old_state(1:nof_variables)=u_c_val(2,1:nof_variables,i)
        candidate(1:nof_variables)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
        if (accepted_update) then
          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        else
          u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
        end if
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif







if ((turbulence.gt.0).or.(passivescalar.gt.0))then

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 end if


call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs, oo4, to4)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
  if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
           call dg_mass_rhs_product_cell(i)
           do mm_var=1,nof_variables
             do mm_idof=1,num_dg_dofs
               u_c_valdg(1,mm_var,mm_idof,i) = to4*u_c_valdg(2,mm_var,mm_idof,i) + oo4*u_c_valdg(3,mm_var,mm_idof,i) - oo4*dt* rhs_sol_mm_dg(mm_idof,mm_var,i)
             end do
           end do



    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        old_state=zero
        candidate=zero
        old_state(1:nof_variables)=u_c_val(3,1:nof_variables,i)
        candidate(1:nof_variables)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
        ((rhs_val(1:nof_variables,i))*(oovolume))))
        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
        if (accepted_update) then
          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        else
          u_c_val(1,1:nof_variables,i)=u_c_val(3,1:nof_variables,i)
        end if
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




if ((turbulence.gt.0).or.(passivescalar.gt.0))then

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo4, passivescalar, to4, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 end if







call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs, oo3, to3)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
  if (dg == 1) then



        call dg_mass_rhs_product_cell(i)


          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = oo3*u_c_valdg(2,mm_var,mm_idof,i) + to3*u_c_valdg(1,mm_var,mm_idof,i) - to3*dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do


    else
        oovolume=1.0d0/ielem_totvolume(i)
        old_state=zero
        candidate=zero
        old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
        candidate(1:nof_variables)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
        ((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
        if (accepted_update) then
          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        end if
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




if ((turbulence.gt.0).or.(passivescalar.gt.0))then

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo3, passivescalar, to3, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 end if







if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta3




subroutine runge_kutta1(n)
!> @brief
!> ssp forward euler scheme
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe,kx
real::avrgs,oovolume
real::t1,t2,t3
kmaxe=xmpielrank(n)




call call_flux_subroutines_3d



#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)

 if (dg == 1) then


          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(1,mm_var,mm_idof,i) - dt* rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do


  else



  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  end if

end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if




if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta1



subroutine runge_kutta2(n)
!> @brief
!> ssp runge kutta 2nd-order scheme
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_3d


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))

end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo2)
#else
	!$omp do
#endif
do i=1,kmaxe
 oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-(dt*oo2*(rhs_val(1:nof_variables,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo2, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(oo2*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo2*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(dt*oo2*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if







if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta2


subroutine runge_kutta5(n)
!> @brief
!> ssp runge kutta 2nd-order scheme for local time stepping
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
logical::accepted_update
real::avrgs,oovolume,update_alpha
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
kmaxe=xmpielrank(n)


call call_flux_subroutines_3d







#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, nof_variables, num_dg_dofs)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)


  if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)

          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - ielem_dtl(i) * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do

    else

  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  old_state=zero
  candidate=zero
  old_state(1:nof_variables)=u_c_val(2,1:nof_variables,i)
  candidate(1:nof_variables)=u_c_val(2,1:nof_variables,i)-(ielem_dtl(i)*(rhs_val(1:nof_variables,i)*oovolume))
  call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
  if (accepted_update) then
    u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
  else
    u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
  end if
  end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(ielem_dtl(i)*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, nof_variables, num_dg_dofs, oo2)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)

  if (dg == 1) then


          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = (oo2*u_c_valdg(2,mm_var,mm_idof,i)) +(oo2*u_c_valdg(1,mm_var,mm_idof,i))- (ielem_dtl(i) *oo2* rhs_sol_mm_dg(mm_idof,mm_var,i))
            end do
          end do

    else



  old_state=zero
  candidate=zero
  old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
  candidate(1:nof_variables)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-&
                             (ielem_dtl(i)*oo2*(rhs_val(1:nof_variables,i)*oovolume))
  call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
  if (accepted_update) then
    u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
  end if
  end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(oo2, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(oo2*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo2*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(ielem_dtl(i)*oo2*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta5

subroutine runge_kutta5_2d(n)
!> @brief
!> ssp runge kutta 2nd-order scheme for local time stepping in 2d
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
logical::accepted_update
real::avrgs,oovolume,update_alpha
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
kmaxe=xmpielrank(n)


call call_flux_subroutines_2d


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, nof_variables, num_dg_dofs)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe

  if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
        u_c_valdg(2,:,:,i)=u_c_valdg(1,:,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - (ielem_dtl(i)* rhs_sol_mm_dg(mm_idof,mm_var,i))
            end do
          end do
    else
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  old_state=zero
  candidate=zero
  old_state(1:nof_variables)=u_c_val(2,1:nof_variables,i)
  candidate(1:nof_variables)=u_c_val(2,1:nof_variables,i)-(ielem_dtl(i)*(rhs_val(1:nof_variables,i)*oovolume))
  call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
  if (accepted_update) then
    u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
  else
    u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
  end if
  end if




end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(ielem_dtl(i)*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, nof_variables, num_dg_dofs, oo2)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe

   if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = (oo2*u_c_valdg(2,mm_var,mm_idof,i))+(oo2*u_c_valdg(1,mm_var,mm_idof,i))-(ielem_dtl(i)* rhs_sol_mm_dg(mm_idof,mm_var,i))
            end do
          end do
    else
  oovolume=1.0d0/ielem_totvolume(i)
  old_state=zero
  candidate=zero
  old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
  candidate(1:nof_variables)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-&
                             (ielem_dtl(i)*oo2*(rhs_val(1:nof_variables,i)*oovolume))
  call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
  if (accepted_update) then
    u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
  end if

  end if

end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_dtl, ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(oo2, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(oo2*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo2*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(ielem_dtl(i)*oo2*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta5_2d


subroutine runge_kutta2_2d(n)
!> @brief
!> ssp runge kutta 2nd-order scheme in 2d
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_2d


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
        u_c_valdg(2,:,:,i)=u_c_valdg(1,:,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - (dt * rhs_sol_mm_dg(mm_idof,mm_var,i))
            end do
          end do
    else
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs, oo2)
#else
	!$omp do
#endif
do i=1,kmaxe
if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = oo2*u_c_valdg(2,mm_var,mm_idof,i) + oo2*u_c_valdg(1,mm_var,mm_idof,i) - (oo2*dt * rhs_sol_mm_dg(mm_idof,mm_var,i))
            end do
          end do
    else
 oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-(dt*oo2*(rhs_val(1:nof_variables,i)*oovolume))
  end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(u_ct_val(2,1:turbulenceequations+passivescalar,i))-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta2_2d



subroutine sol_integ_dg(n)
implicit none
integer,intent(in)::n
integer::i,kmaxe
real,dimension(1:gpu_max_nvar)::solution_integ2


kmaxe=xmpielrank(n)
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(solution_integ2)
#else
	!$omp do
#endif
do i=1,kmaxe
  call solution_integ(i,solution_integ2)
  u_c_val(1,1:nof_variables,i)=solution_integ2(1:nof_variables)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine sol_integ_dg


subroutine sol_integ_dgx(n)
implicit none
integer,intent(in)::n
integer::i,kmaxe
real,dimension(1:gpu_max_nvar)::solution_integ2
real,dimension(1:gpu_max_nvar)::solution_integ_weak
real,dimension(1:gpu_max_nvar)::solution_integ_strong


kmaxe=xmpielrank(n)
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(solution_integ2,solution_integ_weak,solution_integ_strong)
#else
	!$omp do
#endif
do i=1,kmaxe




  call solution_integ(i,solution_integ2)
  u_c_val(1,:,i)=solution_integ2(1:nof_variables)



  if (filtering.eq.1)then
  call solution_integ_s(i,solution_integ_strong)
  call solution_integ_w(i,solution_integ_weak)

  u_cs_val(1,:,i)=solution_integ_strong(1:nof_variables)
  u_cw_val(1,:,i)=solution_integ_weak(1:nof_variables)
  end if


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif




end subroutine sol_integ_dgx

subroutine sol_integ_dg_init(n)
implicit none
integer,intent(in)::n
integer::i,kmaxe
real,dimension(1:gpu_max_nvar)::solution_integ2

kmaxe=xmpielrank(n)
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(solution_integ2)
#else
	!$omp do
#endif
do i=1,kmaxe

  call solution_integ(i,solution_integ2)
  u_c_val(1,:,i)=solution_integ2(1:nof_variables)
  u_e_val(1,:,i)=u_c_val(1,:,i)


end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif






end subroutine sol_integ_dg_init







subroutine runge_kutta1_2d(n)
!> @brief
!> ssp forward euler scheme in 2d
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe,kx
real::avrgs,oovolume
real,dimension(gpu_max_nvar)::tempsol
kmaxe=xmpielrank(n)



call call_flux_subroutines_2d





#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)


  if (dg == 1) then

          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(1,mm_var,mm_idof,i) - dt* rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do

  else

  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  end if


end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



if (averaging.eq.1)then

 call averaging_t(n)

end if






end subroutine runge_kutta1_2d







subroutine runge_kutta3_2d(n)
!> @brief
!> ssp runge kutta 3rd-order scheme in 2d
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
logical::accepted_update
real::avrgs,oovolume,update_alpha
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
kmaxe=xmpielrank(n)


call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
    if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)



            call dg_mass_rhs_product_cell(i)
            do mm_var=1,nof_variables
              do mm_idof=1,num_dg_dofs
                u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
              end do
            end do

    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        old_state=zero
        candidate=zero
        old_state(1:nof_variables)=u_c_val(2,1:nof_variables,i)
        candidate(1:nof_variables)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
        if (accepted_update) then
          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        else
          u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
        end if
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs, oo4, to4)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
    if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)


            call dg_mass_rhs_product_cell(i)
            do mm_var=1,nof_variables
              do mm_idof=1,num_dg_dofs
                u_c_valdg(1,mm_var,mm_idof,i) = to4*u_c_valdg(2,mm_var,mm_idof,i) + oo4*u_c_valdg(3,mm_var,mm_idof,i) - oo4*dt* rhs_sol_mm_dg(mm_idof,mm_var,i)
              end do
            end do
    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        old_state=zero
        candidate=zero
        old_state(1:nof_variables)=u_c_val(3,1:nof_variables,i)
        candidate(1:nof_variables)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
        ((rhs_val(1:nof_variables,i))*(oovolume))))
        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
        if (accepted_update) then
          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        else
          u_c_val(1,1:nof_variables,i)=u_c_val(3,1:nof_variables,i)
        end if
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo4, passivescalar, to4, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs, oo3, to3)
#else
	!$omp do private(oovolume,candidate,old_state,limited,update_alpha,accepted_update)
#endif
do i=1,kmaxe
    if (dg == 1) then


          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = oo3*u_c_valdg(2,mm_var,mm_idof,i) + to3*u_c_valdg(1,mm_var,mm_idof,i) - to3*dt*rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        oovolume=1.0d0/ielem_totvolume(i)
        old_state=zero
        candidate=zero
        old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
        candidate(1:nof_variables)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
        ((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
        if (accepted_update) then
          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        end if
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo3, passivescalar, to3, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if




if (averaging.eq.1)then

 call averaging_t(n)

end if


end subroutine runge_kutta3_2d

subroutine runge_kutta3_2d_mood(n)
!> @brief
!> ssp runge kutta 3rd-order scheme
implicit none
integer,intent(in)::n
integer::i,kmaxe,inds
real::avrgs,oovolume
kmaxe=xmpielrank(n)



if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)

    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select

else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select
end if


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

 if (mood.eq.1)then

 call mood_operator_2(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(4,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
    ielem_mood_o(i)=2
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

call mood_operator_1(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
    ielem_mood_o(i)=1
   else
   u_c_val(1,1:nof_variables,i)=u_c_val(4,1:nof_variables,i)
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
 end if



if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select
else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select
end if

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo4, to4)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo4, passivescalar, to4, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

 if (mood.eq.1)then

 call mood_operator_2(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo4, to4)
#else
	!$omp do
#endif
do i=1,kmaxe
     if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(4,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=2
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo4, to4)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=1
   else
   u_c_val(1,1:nof_variables,i)=u_c_val(4,1:nof_variables,i)
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select
else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    call vortexcalc2d(n)
    end select
end if
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo3, to3)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(3,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if (mood.eq.1)then

 call mood_operator_2(n)
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo3, to3)
#else
	!$omp do
#endif
do i=1,kmaxe
     if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(4,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(3,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=2
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_mood_o, ielem_recalc, ielem_totvolume, rhs_val, u_c_val) &
!$omp& firstprivate(dt, nof_variables, oo3, to3)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (ielem_recalc(i).eq.1)then
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(3,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
    ielem_mood_o(i)=1
   else
   u_c_val(1,1:nof_variables,i)=u_c_val(4,1:nof_variables,i)
   end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, oo3, passivescalar, to3, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



if (averaging.eq.1)then

 call averaging_t(n)

end if

end subroutine runge_kutta3_2d_mood



subroutine runge_kutta4(n)
!> @brief
!> ssp runge kutta 4th-order scheme
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
real::dumpracein,dumpraceout,flops_count
kmaxe=xmpielrank(n)




call call_flux_subroutines_3d










#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - dt * 0.391752226571890 * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*0.391752226571890*(rhs_val(1:nof_variables,i)*oovolume))
    end if

end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*0.391752226571890*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if




    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t8=mpi_wtime()
    prace_t7=pr_t8-pr_t7




   dumpracein=prace_t1
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t1=dumpraceout
   dumpracein=prace_t2
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t2=dumpraceout
   dumpracein=prace_t3
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t3=dumpraceout
   dumpracein=prace_t4
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t4=dumpraceout
   dumpracein=prace_t5
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t5=dumpraceout

   dumpracein=prace_t6
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t6=dumpraceout

   dumpracein=prace_t7
 call mpi_allreduce(dumpracein,dumpraceout,1,mpi_double_precision,mpi_max,mpi_comm_world,ierror)
   prace_t7=dumpraceout

   prace_tx1=prace_t2+prace_t4
   prace_tx2=prace_t1+prace_t3+prace_t5+prace_t6+prace_t7
   prace_tx3=prace_tx1+prace_tx2





    if (n.eq.0)then
    open(133,file=statfile,form='formatted',status='old',action='write',position='append')
    write(133,'(i6,1x,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4,e11.4)')it,prace_tx3,prace_tx1,prace_tx2,prace_t1,prace_t2,prace_t3,prace_t4,prace_t5,prace_t6,prace_t7
    close(133)
    end if


    !$omp end master
    !$omp barrier


    end if

                call call_flux_subroutines_3d


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = 0.444370493651235 * u_c_valdg(2,mm_var,mm_idof,i) + 0.555629506348765 * u_c_valdg(3,mm_var,mm_idof,i) - 0.368410593050371 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=(0.444370493651235 * u_c_val(2,1:nof_variables,i)) + (0.555629506348765 * u_c_val(3,1:nof_variables,i)) - 0.368410593050371 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
    u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.444370493651235*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.555629506348765*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((0.368410593050371))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


            call call_flux_subroutines_3d


#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
   if (dg == 1) then
        u_c_valdg(4,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = 0.620101851488403 * u_c_valdg(2,mm_var,mm_idof,i) + 0.379898148511597 * u_c_valdg(4,mm_var,mm_idof,i) - 0.251891774271694 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(4,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i) = 0.620101851488403 * u_c_val(2,1:nof_variables,i) + 0.379898148511597 * u_c_val(4,1:nof_variables,i) - 0.251891774271694 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(4,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.620101851488403*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.379898148511597*u_ct_val(4,1:turbulenceequations+passivescalar,i))-(((0.251891774271694))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

            call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(5,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(6,mm_var,mm_idof,i) = - dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = 0.178079954393132 * u_c_valdg(2,mm_var,mm_idof,i) + 0.821920045606868 * u_c_valdg(5,mm_var,mm_idof,i) - 0.544974750228521 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(5,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(6,1:nof_variables,i) = - dt * rhs_val(1:nof_variables,i) * oovolume
        u_c_val(1,1:nof_variables,i) = 0.178079954393132 * u_c_val(2,1:nof_variables,i) + 0.821920045606868 * u_c_val(5,1:nof_variables,i) - 0.544974750228521 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
    u_ct_val(5,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(6,1:turbulenceequations+passivescalar,i)=-((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume)))
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.178079954393132*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.821920045606868*u_ct_val(5,1:turbulenceequations+passivescalar,i))-(((0.544974750228521))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


    call call_flux_subroutines_3d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = (0.00683325884039 * u_c_valdg(2,mm_var,mm_idof,i)) + (0.517231671970585 * u_c_valdg(4,mm_var,mm_idof,i)) + (0.12759831133288 * u_c_valdg(5,mm_var,mm_idof,i)) + (0.34833675773694 * u_c_valdg(1,mm_var,mm_idof,i)) + (0.08460416338212 * u_c_valdg(6,mm_var,mm_idof,i)) - 0.22600748319395 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(1,1:nof_variables,i) = (0.00683325884039 * u_c_val(2,1:nof_variables,i)) + (0.517231671970585 * u_c_val(4,1:nof_variables,i)) +  (0.12759831133288 * u_c_val(5,1:nof_variables,i)) + (0.34833675773694 * u_c_val(1,1:nof_variables,i)) + (0.08460416338212 * u_c_val(6,1:nof_variables,i)) - (0.22600748319395 * dt * rhs_val(1:nof_variables,i) * oovolume)
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
   u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.00683325884039*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.517231671970585*u_ct_val(4,1:turbulenceequations+passivescalar,i))+&
				(0.12759831133288*u_ct_val(5,1:turbulenceequations+passivescalar,i))+(0.34833675773694*u_ct_val(1,1:turbulenceequations+passivescalar,i))+&
				(0.08460416338212*u_ct_val(6,1:turbulenceequations+passivescalar,i))-(0.22600748319395*(dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume)))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


if (averaging.eq.1)then

 call averaging_t(n)

end if


end subroutine runge_kutta4


subroutine call_flux_subroutines_3d
implicit none
real::dumpracein,dumpraceout
integer::kmaxe,i,l,iqp,ngp
kmaxe=xmpielrank(n)



    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t1=mpi_wtime()
     !$omp end master
     !$omp barrier

    end if


    if (realgas.eq.1)then

        call normalise_species(n)

    end if

    if (dg.eq.1)then


    call sol_integ_dg(n) ! calculates cell average of dg solution for fv


    end if

      if ((dg.eq.1).and.(filtering.eq.1))then

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& map(alloc: ielem_filtered) &
!$omp& firstprivate(kmaxe)
#else
	!$omp do
#endif
        do i=1,kmaxe
        ielem_filtered(i)=0
        end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

    end if



	if ((dg.eq.1).and.(filtering.eq.1))then


            call sol_integ_dgx(n)

            call apply_filter_dg(n)

          end if








    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t2=mpi_wtime()
    prace_t1=pr_t2-pr_t1
    !$omp end master
    !$omp barrier

    end if

    if (fastest.eq.1) then

        call exchange_lower(n)

    else

        call exchange_higher(n)

    end if


    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t3=mpi_wtime()
    prace_t2=pr_t3-pr_t2
    !$omp end master
    !$omp barrier

    end if


    if (dg == 1) then


        call reconstruct_dg(n) ! extrapolates solution to faces

        if(multispecies.eq.1)then
          if (bound_lim == 1) then

            call vfbp_limiter

          end if
        end if


        call trouble_indicator1(n) ! checks for troubled cells


    end if


    call arbitrary_order(n)
!
! #ifdef gpu
!     !$omp target update from(rec_uleft)
!     !$omp barrier
!     !$omp master
!     do i=1,kmaxe
!           write(300+n,*) "element",ielem_ihexgl(i)
!           write(300+n,*) "sol",u_c_val(1,1,i)
!         do l=1,ielem_ifca(i)
!                   if (ielem_types_faces(l,i).eq.5)then
! 					iqp=qp_quad
! 				  else
! 					iqp=qp_triangle
! 				  end if
! 				  do ngp=1,iqp
! 				  write(300+n,*)l,ngp,rec_uleft(1,l,ngp,i)
!                   end do
!         end do
!     end do
!
!     !$omp end master
! #endif
!
!

     if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t4=mpi_wtime()
    prace_t3=pr_t4-pr_t3
    !$omp end master
    !$omp barrier
    end if

        if (dg == 1) then


        call trouble_indicator2(n) ! changes dg to fv


    end if


    call exhboundhigher(n)





    if (dg.eq.1)then


        call exhboundhigher_dg(n)


        if (itestcase.eq.4)then

          if( br2_yn.eq.2) then


          call reconstruct_br2_dg(n)

          call exhboundhigher_dg2(n)


          end if

          if( br2_yn.eq.0) then

          call viscous_dg_ggs(n)

          end if



        end if

    end if

!     if (statistics.eq.1)then
!    !$omp barrier
!    !$omp master
!    pr_t4=mpi_wtime()
!    prace_t3=pr_t4-pr_t3
!    !$omp end master
!    !$omp barrier
!    end if






    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t5=mpi_wtime()
    prace_t4=pr_t5-pr_t4
    !$omp end master
    !$omp barrier
    end if


    if (adda.eq.1)then

    if (rungekutta.eq.11)then

    if (iscoun.eq.1)then

    call fix_dissipation(n)

    call exchange_adda_diss(n)

    call fix_dissipation2(n)

    end if
    else

    call fix_dissipation(n)

    call exchange_adda_diss(n)

    call fix_dissipation2(n)

    end if


    end if



    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t6=mpi_wtime()
    prace_t5=pr_t6-pr_t5
    !$omp end master
    !$omp barrier
    end if



    !modifies rhs
    select case(itestcase)
    case(1,2)

    call calculate_fluxeshi(n)




    case(3)

    call calculate_fluxeshi_convective(n)





    if ((source_active.eq.1))then

    call sources_computation_rot(n)

    end if
    if ((realgas.eq.1))then

    call sources_computation(n)

    end if
    case(4)

    call calculate_fluxeshi_convective(n)

    call calculate_fluxeshi_diffusive(n)

    if ((source_active.eq.1))then

    call sources_computation_rot(n)

    end if
    if ((turbulence.eq.1).or.(realgas.eq.1))then

    call sources_computation(n)

    end if


      call vortexcalc(n)


    end select



    if (initcond.eq.95)then


    call  enstrophy_calc(n)

    end if


    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t7=mpi_wtime()
    prace_t6=pr_t7-pr_t6
    !$omp end master
    !$omp barrier
    end if














end subroutine call_flux_subroutines_3d





subroutine normalise_species(n)
  implicit none

  integer, intent(in) :: n
  integer :: i, k, kmaxe

  real :: rho, sumrhoy, scale, rel_err
!   real, parameter :: epsrho = 1.0d-18
!   real, parameter :: epssum = 1.0d-300
!
!   ! Much larger than roundoff, but still numerically tiny.
!   real, parameter :: tol_mass_closure = 1.0d-12

  real, parameter :: epsrho = 1.0d-14
real, parameter :: epssum = 1.0d-300
real, parameter :: tol_mass_closure = 1.0d-12

  real, dimension(1:gpu_max_species)   :: rhoy
  real, dimension(1:gpu_max_nvar) :: leftv
  logical :: bad, need_repair

  kmaxe = xmpielrank(n)

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(i,k,rho,sumrhoy,scale,rel_err,rhoy,leftv,bad,need_repair)
#else
  !$omp do private(i,k,rho,sumrhoy,scale,rel_err,rhoy,leftv,bad,need_repair)
#endif
  do i = 1, kmaxe

    ! First make the full conservative state sane.
    leftv(1:nof_variables) = u_c_val(1,1:nof_variables,i)
!     call fix_conservative_state(leftv)
    u_c_val(1,1:nof_variables,i) = leftv(1:nof_variables)

    rho = u_c_val(1,1,i)

    if (rho <= epsrho .or. rho /= rho .or. abs(rho) > huge(rho)*0.5d0) then

      rho = epsrho
      u_c_val(1,1,i) = rho

      do k = 1, nof_species
        u_c_val(1,dimensiona+3+k,i) = rho * rg_vf(k)
      end do

      leftv(1:nof_variables) = u_c_val(1,1:nof_variables,i)
!       call fix_conservative_state(leftv)
      u_c_val(1,1:nof_variables,i) = leftv(1:nof_variables)

      cycle

    end if

    need_repair = .false.
    sumrhoy = 0.0d0

    do k = 1, nof_species

      rhoy(k) = u_c_val(1,dimensiona+3+k,i)

      bad = (rhoy(k) /= rhoy(k))
      if (.not. bad) bad = (abs(rhoy(k)) > huge(rhoy(k))*0.5d0)

      if (bad) then
        rhoy(k) = 0.0d0
        need_repair = .true.
      elseif (rhoy(k) < 0.0d0) then
        rhoy(k) = 0.0d0
        need_repair = .true.
      end if

      sumrhoy = sumrhoy + rhoy(k)

    end do

    rel_err = abs(sumrhoy - rho) / max(abs(rho), epsrho)

    if (rel_err > tol_mass_closure) need_repair = .true.

    if (need_repair) then

      sumrhoy = 0.0d0
      do k = 1, nof_species
        rhoy(k) = max(rhoy(k), 0.0d0)
        sumrhoy = sumrhoy + rhoy(k)
      end do

      if (sumrhoy > epssum) then

        ! Smoothest conservative repair: preserve composition ratios.
        scale = rho / sumrhoy

        do k = 1, nof_species
          rhoy(k) = rhoy(k) * scale
        end do

	      else

	        ! Completely broken species state: reset to reference mass fractions.
	        do k = 1, nof_species
	          rhoy(k) = rho * rg_vf(k)
	        end do

	      end if

    end if

    do k = 1, nof_species
      u_c_val(1,dimensiona+3+k,i) = rhoy(k)
    end do

    leftv(1:nof_variables) = u_c_val(1,1:nof_variables,i)
    call fix_conservative_state(leftv)
    u_c_val(1,1:nof_variables,i) = leftv(1:nof_variables)

  end do

#ifdef gpu
!$omp end target teams distribute parallel do
#else
  !$omp end do
#endif


end subroutine normalise_species

subroutine call_flux_subroutines_2d
implicit none
integer::i,iconsidered




     if (realgas.eq.1)call normalise_species(n)



    if (dg.eq.1)then
    call sol_integ_dg(n)
    end if

    if (fastest.eq.1) then
        call exchange_lower(n)
    else
        call exchange_higher(n)
    end if


    if (dg == 1) then

        call reconstruct_dg(n)
        call trouble_indicator1(n)

    end if



        call arbitrary_order(n)


    if (dg == 1) then

        call trouble_indicator2(n)

    end if

    if (bound_lim == 1) then
            call vfbp_limiter
          end if


    call exhboundhigher(n)

    if (dg.eq.1)then
    call exhboundhigher_dg(n)

    if (itestcase.eq.4)then

          if( br2_yn.eq.2) then

          call reconstruct_br2_dg(n)

          call exhboundhigher_dg2(n)

          end if

          if( br2_yn.eq.0) then
          call viscous_dg_ggs(n)
          end if



        end if



    end if

    !modifies rhs
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    if (realgas.eq.1)then
    if (rg_relax.eq.2)then
    call sources_computation2d(n)
    end if
    end if
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
     if ((turbulence.eq.1).or.(realgas.eq.1))then
     call sources_computation2d(n)
     end if



    end select






end subroutine call_flux_subroutines_2d


subroutine runge_kutta4_2d(n)
!> @brief
!> ssp runge kutta 4th-order scheme in 2d
implicit none
integer :: mm_idof
integer :: mm_var
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)



call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)

    if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = u_c_valdg(2,mm_var,mm_idof,i) - dt * 0.391752226571890 * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*0.391752226571890*(rhs_val(1:nof_variables,i)*oovolume))
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*0.391752226571890*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)

    if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = 0.444370493651235 * u_c_valdg(2,mm_var,mm_idof,i) + 0.555629506348765 * u_c_valdg(3,mm_var,mm_idof,i) - 0.368410593050371 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=(0.444370493651235 * u_c_val(2,1:nof_variables,i)) + (0.555629506348765 * u_c_val(3,1:nof_variables,i)) - 0.368410593050371 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
    u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.444370493651235*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.555629506348765*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((0.368410593050371))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)

    if (dg == 1) then
        u_c_valdg(4,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = 0.620101851488403 * u_c_valdg(2,mm_var,mm_idof,i) + 0.379898148511597 * u_c_valdg(4,mm_var,mm_idof,i) - 0.251891774271694 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(4,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i) = 0.620101851488403 * u_c_val(2,1:nof_variables,i) + 0.379898148511597 * u_c_val(4,1:nof_variables,i) - 0.251891774271694 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(4,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.620101851488403*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.379898148511597*u_ct_val(4,1:turbulenceequations+passivescalar,i))-(((0.251891774271694))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)

    if (dg == 1) then
        u_c_valdg(5,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(6,mm_var,mm_idof,i) = - dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = 0.178079954393132 * u_c_valdg(2,mm_var,mm_idof,i) + 0.821920045606868 * u_c_valdg(5,mm_var,mm_idof,i) - 0.544974750228521 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(5,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(6,1:nof_variables,i) = - dt * rhs_val(1:nof_variables,i) * oovolume
        u_c_val(1,1:nof_variables,i) = 0.178079954393132 * u_c_val(2,1:nof_variables,i) + 0.821920045606868 * u_c_val(5,1:nof_variables,i) - 0.544974750228521 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
    u_ct_val(5,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(6,1:turbulenceequations+passivescalar,i)=-((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume)))
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.178079954393132*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.821920045606868*u_ct_val(5,1:turbulenceequations+passivescalar,i))-(((0.544974750228521))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if

call call_flux_subroutines_2d
if (fastest /= 1)then
    if (itestcase.eq.4)then
        call vortexcalc2d(n)
    end if
end if

#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhs_sol_mm_dg, rhs_val, u_c_val, u_c_valdg) &
!$omp& firstprivate(dg, dt, nof_variables, num_dg_dofs)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)

    if (dg == 1) then
          call dg_mass_rhs_product_cell(i)
          do mm_var=1,nof_variables
            do mm_idof=1,num_dg_dofs
              u_c_valdg(1,mm_var,mm_idof,i) = (0.00683325884039 * u_c_valdg(2,mm_var,mm_idof,i)) + (0.517231671970585 * u_c_valdg(4,mm_var,mm_idof,i)) + (0.12759831133288 * u_c_valdg(5,mm_var,mm_idof,i)) + (0.34833675773694 * u_c_valdg(1,mm_var,mm_idof,i)) + (0.08460416338212 * u_c_valdg(6,mm_var,mm_idof,i)) - 0.22600748319395 * dt * rhs_sol_mm_dg(mm_idof,mm_var,i)
            end do
          end do
    else
        u_c_val(1,1:nof_variables,i) = (0.00683325884039 * u_c_val(2,1:nof_variables,i)) + (0.517231671970585 * u_c_val(4,1:nof_variables,i)) +  (0.12759831133288 * u_c_val(5,1:nof_variables,i)) + (0.34833675773694 * u_c_val(1,1:nof_variables,i)) + (0.08460416338212 * u_c_val(6,1:nof_variables,i)) - (0.22600748319395 * dt * rhs_val(1:nof_variables,i) * oovolume)
    end if
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#if defined(gpu) || defined(xpu)
!$omp target teams distribute parallel do &
!$omp& private(oovolume) &
!$omp& map(alloc: ielem_totvolume, rhst_val, u_ct_val) &
!$omp& firstprivate(dt, passivescalar, turbulenceequations)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
   u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.00683325884039*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.517231671970585*u_ct_val(4,1:turbulenceequations+passivescalar,i))+&
				(0.12759831133288*u_ct_val(5,1:turbulenceequations+passivescalar,i))+(0.34833675773694*u_ct_val(1,1:turbulenceequations+passivescalar,i))+&
				(0.08460416338212*u_ct_val(6,1:turbulenceequations+passivescalar,i))-(0.22600748319395*(dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume)))
end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

if (averaging.eq.1)then

 call averaging_t(n)

end if



end subroutine runge_kutta4_2d




subroutine implicit_times(n)
!> @brief
!> implicit approximately factored time stepping scheme
implicit none
integer::i,j,k,kmaxe,kill_nan_global
integer,intent(in)::n
#ifdef gpu
logical::kill_nan_reduce
#endif
logical::bad_impdu,accepted_update
real::verysmall,du,uold, umin, dumax, unew,update_alpha
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
verysmall = tolsmall


if (realgas.eq.1)call normalise_species(n)

kmaxe=xmpielrank(n)
if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)

    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    if ((source_active.eq.1))then
    call sources_computation_rot(n)
    end if
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if ((source_active.eq.1))then
    call sources_computation_rot(n)
    end if
    call vortexcalc(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if
    end select

else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    if ((source_active.eq.1))then
    call sources_computation_rot(n)
    end if
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    if ((source_active.eq.1))then
    call sources_computation_rot(n)
    end if
    call vortexcalc(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if
    end select
end if

if (relax.eq.3)then

 !call relaxation_lumfree(n)

 else

if (lowmemory.eq.0)then

 call relaxation(n)

 else

 call relaxation_lm(n)

 end if

  end if

 kill_nan=0
#ifdef gpu
 kill_nan_reduce=.false.
#endif

#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& reduction(.or.:kill_nan_reduce)
#else
	!$omp do private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update) reduction(max:kill_nan)
#endif
	do i=1,kmaxe
	    bad_impdu=.false.
	    do j=1,nof_variables
	      if ((impdu(i,j).ne.impdu(i,j)).or.(abs(impdu(i,j)).gt.fv_update_max_abs)) bad_impdu=.true.
	    end do
	    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
	        if ((impdu(i,j).ne.impdu(i,j)).or.(abs(impdu(i,j)).gt.fv_update_max_abs)) bad_impdu=.true.
	      end do
	    end if

	    if (bad_impdu)then
	      if (realgas.eq.1)then
	        old_state=zero
	        old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
	        call fv_explicit_local_candidate(i,old_state,candidate)
	        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
	        if (accepted_update) then
	          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
	          impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
	          if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	            impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
	            -update_alpha*ielem_dtl(i)*rhst_val(1:turbulenceequations+passivescalar,i)/ielem_totvolume(i)
	          end if
	        else
	          if (.not.fv_state_update_admissible(old_state)) then
#ifdef gpu
	            kill_nan_reduce=.true.
#else
	            kill_nan=1
#endif
	          end if
	          impdu(i,1:nof_variables)=zero
	          if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	            impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	          end if
	        end if
	      else
#ifdef gpu
	        kill_nan_reduce=.true.
#else
	        kill_nan=1
#endif
	        impdu(i,1:nof_variables)=zero
	        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	        end if
	      end if
    else
      old_state=zero
      candidate=zero
      old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      candidate(1:nof_variables)=old_state(1:nof_variables)+impdu(i,1:nof_variables)
      call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
      if (accepted_update) then
        u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
          update_alpha*impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
        end if
	      else
	        if (realgas.eq.1)then
	          call fv_explicit_local_candidate(i,old_state,candidate)
	          call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
	          if (accepted_update) then
	            u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
	            impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
	            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	              impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
	              -update_alpha*ielem_dtl(i)*rhst_val(1:turbulenceequations+passivescalar,i)/ielem_totvolume(i)
	            end if
	          else
	            impdu(i,1:nof_variables)=zero
	            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	              impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	            end if
	          end if
	        else
	          impdu(i,1:nof_variables)=zero
	          if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	            impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	          end if
	        end if
	      end if
	    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
kill_nan=0
if (kill_nan_reduce) kill_nan=1
#else
!$omp end do
#endif


    !$omp barrier
    !$omp master
    kill_nan_global = 0
    call mpi_allreduce(kill_nan,kill_nan_global,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
    kill_nan = kill_nan_global
    if (kill_nan.eq.1)then
        kill=1
        if (n.eq.0)then
      print*,"killed due to divergence in implicit time stepping"
    end if
    end if

    !$omp end master
    !$omp barrier

    if (kill_nan.eq.1)then
        return
    end if

if (realgas.eq.1)then
if (rg_relax.eq.2)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
	!$omp do
#endif
do i=1,kmaxe
! Real-gas sources are assembled through rhs_val in sources_computation.
! Do not mutate u_c_val here with a local pseudo-time source step.
continue
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if
if (realgas.eq.1)call normalise_species(n)




if ((passivescalar.gt.0).or.(turbulence.gt.0))then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private (du,uold,umin,dumax,unew)
#else
	!$omp do
#endif
!   do i=1,kmaxe
!   do k=1,turbulenceequations+passivescalar
!   if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
!   u_ct_val(1,k,i)=u_ct_val(1,k,i)+0.4*impdu(i,nof_variables+k)
!   end if
!   end do
! end do

do i=1,kmaxe
  do k=1,turbulenceequations+passivescalar

    du   = 0.2d0*impdu(i,nof_variables+k)
    uold = u_ct_val(1,k,i)
    umin = 1.0d-12   !  floor

    ! optional relative limiter
    dumax = 0.25d0*max(abs(uold), umin)
    if (du.gt. dumax) du =  dumax
    if (du.lt.-dumax) du = -dumax

    unew = uold + du

    if (unew.lt.umin) then
      unew = umin
    end if

    u_ct_val(1,k,i) = unew

  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if










end subroutine implicit_times


subroutine implicit_times_2d(n)
!> @brief
!> implicit approximately factored time stepping scheme 2d
implicit none
integer::i,k,kmaxe,j,kill_nan_global
integer,intent(in)::n
#ifdef gpu
logical::kill_nan_reduce
#endif
logical::bad_impdu,accepted_update
real::verysmall,du,uold, umin, dumax, unew,update_alpha
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
verysmall = tolsmall

  if (realgas.eq.1)call normalise_species(n)
kmaxe=xmpielrank(n)
if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)

    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    if (realgas.eq.1)then
    call sources_computation2d(n)

    end if
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)

   ! call vortexcalc2d(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation2d(n)
    end if
    end select

else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    if (realgas.eq.1)then
    call sources_computation2d(n)
    end if
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)


    !call vortexcalc2d(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation2d(n)
    end if
    end select
end if



 if (relax.eq.3)then

 !call relaxation_lumfree(n)

 else
 if (lowmemory.eq.0)then

 call relaxation2d(n)

 else

 call relaxation_lm2d(n)

 end if
  end if

 kill_nan=0
#ifdef gpu
 kill_nan_reduce=.false.
#endif
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update) &
!$omp& reduction(.or.:kill_nan_reduce)
#else
!$omp do private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update) reduction(max:kill_nan)
#endif
	do i=1,kmaxe
	    bad_impdu=.false.
	    do j=1,nof_variables
	    if ((impdu(i,j).ne.impdu(i,j)).or.(abs(impdu(i,j)).gt.fv_update_max_abs)) bad_impdu=.true.
	    end do
	    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
	        if ((impdu(i,j).ne.impdu(i,j)).or.(abs(impdu(i,j)).gt.fv_update_max_abs)) bad_impdu=.true.
	      end do
	    end if

	    if (bad_impdu)then

	      if (realgas.eq.1)then
	        old_state=zero
	        old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
	        call fv_explicit_local_candidate(i,old_state,candidate)
	        call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
	        if (accepted_update) then
	          u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
	          impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
	          if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	            impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
	            -update_alpha*ielem_dtl(i)*rhst_val(1:turbulenceequations+passivescalar,i)/ielem_totvolume(i)
	          end if
	        else
	          if (.not.fv_state_update_admissible(old_state)) then
#ifdef gpu
	            kill_nan_reduce=.true.
#else
	            kill_nan=1
#endif
	          end if
	          impdu(i,1:nof_variables)=zero
	          if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	            impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	          end if
	        end if
	      else
#ifdef gpu
	        kill_nan_reduce=.true.
#else
	        kill_nan=1
#endif
	        impdu(i,1:nof_variables)=zero
	        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	        end if
	      end if
    else
      old_state=zero
      candidate=zero
      old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      candidate(1:nof_variables)=old_state(1:nof_variables)+impdu(i,1:nof_variables)
      call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
      if (accepted_update) then
        u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
          update_alpha*impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
        end if
	      else
	        if (realgas.eq.1)then
	          call fv_explicit_local_candidate(i,old_state,candidate)
	          call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
	          if (accepted_update) then
	            u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
	            impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
	            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	              impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
	              -update_alpha*ielem_dtl(i)*rhst_val(1:turbulenceequations+passivescalar,i)/ielem_totvolume(i)
	            end if
	          else
	            impdu(i,1:nof_variables)=zero
	            if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	              impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	            end if
	          end if
	        else
	          impdu(i,1:nof_variables)=zero
	          if ((turbulence.gt.0).or.(passivescalar.gt.0))then
	            impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
	          end if
	        end if
	      end if
	    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
kill_nan=0
if (kill_nan_reduce) kill_nan=1
#else
!$omp end do
#endif

    !$omp barrier
    !$omp master
    kill_nan_global = 0
    call mpi_allreduce(kill_nan,kill_nan_global,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
    kill_nan = kill_nan_global
    if (kill_nan.eq.1)then
        kill=1
        if (n.eq.0)then
      print*,"killed due to divergence in implicit time stepping"
    end if
    end if
    !$omp end master
    !$omp barrier

    if (kill_nan.eq.1)then
        return
    end if


 if (realgas.eq.1)call normalise_species(n)


if (realgas.eq.1)then
if (rg_relax.eq.2)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
! Real-gas sources are assembled through rhs_val in sources_computation.
! Do not mutate u_c_val here with a local pseudo-time source step.
continue
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if
if (realgas.eq.1)call normalise_species(n)

if ((passivescalar.gt.0).or.(turbulence.gt.0))then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private (du,uold,umin,dumax,unew)
#else
!$omp do
#endif
do i=1,kmaxe
  do k=1,turbulenceequations+passivescalar

    du   = 0.2d0*impdu(i,nof_variables+k)
    uold = u_ct_val(1,k,i)
    umin = 1.0d-12   !  floor

    ! optional relative limiter
    dumax = 0.25d0*max(abs(uold), umin)
    if (du.gt. dumax) du =  dumax
    if (du.lt.-dumax) du = -dumax

    unew = uold + du

    if (unew.lt.umin) then
      unew = umin
    end if

    u_ct_val(1,k,i) = unew

  end do
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end if










end subroutine implicit_times_2d


subroutine dual_time(n)
!> @brief
!> dual time stepping
implicit none
integer::i,j,k,kmaxe,jj,kill_nan_global
integer,intent(in)::n
#ifdef gpu
logical::kill_nan_reduce
real::allresdt_reduce
#endif
logical::bad_impdu,accepted_update
real::verysmall,update_alpha
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
verysmall = tolsmall

inner_tol=reslimit

kmaxe=xmpielrank(n)


if (it.eq.restart)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(1,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if





firsti=0.0d0
do jj=1,upperlimit
      rsumfacei=zero;allresdt=zero;dummy3i=zero;
   if (jj.eq.1)then
    iscoun=1
      else
      iscoun=2
      end if


#ifdef gpu
!$omp target update to(iscoun)
#endif

call call_flux_subroutines_3d


if (relax.eq.3)then



!call relaxation_lumfree(n)


else


if (lowmemory.eq.0)then

 call relaxation(n)

 else

 call relaxation_lm(n)

 end if
 end if


kill_nan=0
allresdt = 0.0d0
#ifdef gpu
kill_nan_reduce=.false.
allresdt_reduce = 0.0d0
#endif


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(rsumfacei,j) reduction(+:allresdt_reduce) reduction(.or.:kill_nan_reduce)
#else
!$omp do private(rsumfacei,j) reduction(+:allresdt) reduction(max:kill_nan)
#endif
do i = 1, kmaxe
    rsumfacei = sqrt( impdu(i,1)**2 + impdu(i,2)**2 + impdu(i,3)**2 + &
                      impdu(i,4)**2 + impdu(i,5)**2 )

#ifdef gpu
    allresdt_reduce = allresdt_reduce + rsumfacei * ielem_totvolume(i)
#else
    allresdt = allresdt + rsumfacei * ielem_totvolume(i)
#endif

    if ( (impdu(i,1) /= impdu(i,1)) .or. &
         (impdu(i,2) /= impdu(i,2)) .or. &
         (impdu(i,3) /= impdu(i,3)) .or. &
         (impdu(i,4) /= impdu(i,4)) .or. &
         (impdu(i,5) /= impdu(i,5)) ) then
#ifdef gpu
        kill_nan_reduce = .true.
#else
        kill_nan = 1
#endif
    end if
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
        if (impdu(i,j).ne.impdu(i,j)) then
#ifdef gpu
          kill_nan_reduce = .true.
#else
          kill_nan = 1
#endif
        end if
      end do
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
allresdt=allresdt_reduce
kill_nan=0
if (kill_nan_reduce) kill_nan=1
#else
!$omp end do
#endif

    !$omp barrier
    !$omp master
    kill_nan_global = 0
    call mpi_allreduce(kill_nan,kill_nan_global,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
    kill_nan = kill_nan_global
    if (kill_nan.eq.1)then
        kill=1
        if (n.eq.0)then
      print*,"killed due to divergence in implicit time stepping"
    end if
    end if
    !$omp end master
    !$omp barrier

    if (kill_nan.eq.1)then
        return
    end if


!$omp master
dummy3i=zero







call mpi_allreduce(allresdt,dummy3i,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allresdt=dummy3i/totalvolume

if (allresdt.gt.firsti)then
firsti=allresdt
end if

allresdt=allresdt/firsti


if (n.eq.0)then
				  open(77,file='res1.txt',form='formatted',action='write',position='append')
				  write(77,*)allresdt,jj,it
				  close(77)

end if

!$omp end master
!$omp barrier




  if ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#else
!$omp do private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#endif
  do i=1,kmaxe
    bad_impdu=.false.
    do j=1,nof_variables
      if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
    end do
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
        if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
      end do
    end if
    if (bad_impdu)then
      impdu(i,1:nof_variables)=zero
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
      end if
    else
      old_state=zero
      candidate=zero
      old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      candidate(1:nof_variables)=old_state(1:nof_variables)+impdu(i,1:nof_variables)
      call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
      if (accepted_update) then
        u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
          update_alpha*impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
          do k=1,turbulenceequations+passivescalar
            if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
              u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
            end if
          end do
        end if
      else
        impdu(i,1:nof_variables)=zero
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
        end if
      end if
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if (realgas.eq.1)then
if (rg_relax.eq.2)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
! Real-gas sources are assembled through rhs_val in sources_computation.
! Do not mutate u_c_val here with a local pseudo-time source step.
continue
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if
if (realgas.eq.1)call normalise_species(n)










  exit

else
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#else
!$omp do private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#endif
 do i=1,kmaxe
    bad_impdu=.false.
    do j=1,nof_variables
      if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
    end do
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
        if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
      end do
    end if
    if (bad_impdu)then
      impdu(i,1:nof_variables)=zero
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
      end if
    else
      old_state=zero
      candidate=zero
      old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      candidate(1:nof_variables)=old_state(1:nof_variables)+impdu(i,1:nof_variables)
      call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
      if (accepted_update) then
        u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
          update_alpha*impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
          do k=1,turbulenceequations+passivescalar
            if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
              u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
            end if
          end do
        end if
      else
        impdu(i,1:nof_variables)=zero
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
        end if
      end if
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if (realgas.eq.1)then
if (rg_relax.eq.2)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
! Real-gas sources are assembled through rhs_val in sources_computation.
! Do not mutate u_c_val here with a local pseudo-time source step.
continue
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if
if (realgas.eq.1)call normalise_species(n)

end if


end do

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  !u_c_val(1,1:nof_variables,i)=2.0*u_c_val(2,1:nof_variables,i)-u_c_val(3,1:nof_variables,i)


  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(2,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  !u_ct_val(1,:,i)=2.0*u_ct_val(2,:,i)-u_ct_val(3,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if (averaging.eq.1)then

 call averaging_t(n)

end if











end subroutine dual_time




















subroutine dual_time_ex(n)
!> @brief
!> dual time stepping 2d
implicit none
integer::i,k,kmaxe,jj
integer,intent(in)::n
#ifdef gpu
real::allresdt_reduce
#endif
real::verysmall
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
verysmall = tolsmall

inner_tol=reslimit

kmaxe=xmpielrank(n)


if (it.eq.restart)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(1,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if





firsti=0.0d0
do jj=1,upperlimit
      rsumfacei=zero;allresdt=zero;dummy3i=zero;






if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)

    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    call vortexcalc(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if
    end select

else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi(n)
    case(3)
    call calculate_fluxeshi_convective(n)
    case(4)
    call calculate_fluxeshi_convective(n)
    call calculate_fluxeshi_diffusive(n)
    call vortexcalc(n)
    if ((turbulence.eq.1).or.(realgas.eq.1))then
    call sources_computation(n)
    end if
    end select
end if


call relaxation_ex(n)





#ifdef gpu
allresdt_reduce=0.0d0
!$omp target teams distribute parallel do &
!$omp& private(rsumfacei) reduction(+:allresdt_reduce)
#else
!$omp do private(rsumfacei) reduction(+:allresdt)
#endif
do i=1,kmaxe
      rsumfacei=sqrt(((impdu(i,1))**2)+((impdu(i,2))**2)+((impdu(i,3))**2)+((impdu(i,4))**2)+((impdu(i,5))**2))
#ifdef gpu
      allresdt_reduce=allresdt_reduce+(rsumfacei*ielem_totvolume(i))
#else
      allresdt=allresdt+(rsumfacei*ielem_totvolume(i))
#endif
end do
#ifdef gpu
!$omp end target teams distribute parallel do
allresdt=allresdt_reduce
#else
!$omp end do
#endif

dummy3i=zero

call mpi_allreduce(allresdt,dummy3i,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allresdt=dummy3i/totalvolume



if (allresdt.gt.firsti)then
firsti=allresdt
end if

allresdt=allresdt/firsti

!$omp master
    if (n.eq.0)then
    write(777,*)allresdt,jj,it
    end if


!$omp end master
!$omp barrier


  if ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
  do i=1,kmaxe
 u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do k=1,turbulenceequations+passivescalar
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
    end do
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  exit

else
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
 do i=1,kmaxe
 u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do k=1,turbulenceequations+passivescalar
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
    end do
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if






end do






#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=2.0*u_c_val(2,1:nof_variables,i)-u_c_val(3,1:nof_variables,i)


  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(2,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  u_ct_val(1,:,i)=2.0*u_ct_val(2,:,i)-u_ct_val(3,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif












if (averaging.eq.1)then

 call averaging_t(n)

end if



end subroutine dual_time_ex





subroutine dual_time_ex_2d(n)
!> @brief
!> dual time stepping 2d
implicit none
integer::i,k,kmaxe,jj
integer,intent(in)::n
#ifdef gpu
real::allresdt_reduce
#endif
real::verysmall
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
verysmall = tolsmall

inner_tol=reslimit

kmaxe=xmpielrank(n)


if (it.eq.restart)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(1,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if





firsti=0.0d0
do jj=1,upperlimit
      rsumfacei=zero;allresdt=zero;dummy3i=zero;






if (fastest.eq.1)then
    call exchange_lower(n)
    call arbitrary_order(n)
    call exhboundhigher(n)

    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    call vortexcalc2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select

else
    call exchange_higher(n)
    call arbitrary_order(n)
    call exhboundhigher(n)
    select case(itestcase)
    case(1,2)
    call calculate_fluxeshi2d(n)
    case(3)
    call calculate_fluxeshi_convective2d(n)
    case(4)
    call calculate_fluxeshi_convective2d(n)
    call calculate_fluxeshi_diffusive2d(n)
    call vortexcalc2d(n)
    if (turbulence.eq.1)then
    call sources_computation2d(n)
    end if
    end select
end if


call relaxation_ex(n)





#ifdef gpu
allresdt_reduce=0.0d0
!$omp target teams distribute parallel do &
!$omp& private(rsumfacei) reduction(+:allresdt_reduce)
#else
!$omp do private(rsumfacei) reduction(+:allresdt)
#endif
do i=1,kmaxe
      rsumfacei=sqrt(((impdu(i,1))**2)+((impdu(i,2))**2)+((impdu(i,3))**2)+((impdu(i,4))**2))
#ifdef gpu
      allresdt_reduce=allresdt_reduce+(rsumfacei*ielem_totvolume(i))
#else
      allresdt=allresdt+(rsumfacei*ielem_totvolume(i))
#endif
end do
#ifdef gpu
!$omp end target teams distribute parallel do
allresdt=allresdt_reduce
#else
!$omp end do
#endif

dummy3i=zero

call mpi_allreduce(allresdt,dummy3i,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allresdt=dummy3i/totalvolume



if (allresdt.gt.firsti)then
firsti=allresdt
end if

allresdt=allresdt/firsti

!$omp master
    if (n.eq.0)then
    write(777,*)allresdt,jj,it
    end if


!$omp end master
!$omp barrier


  if ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
  do i=1,kmaxe
 u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do k=1,turbulenceequations+passivescalar
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
    end do
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  exit

else
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
 do i=1,kmaxe
 u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do k=1,turbulenceequations+passivescalar
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
    end do
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if






end do






#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=2.0*u_c_val(2,1:nof_variables,i)-u_c_val(3,1:nof_variables,i)


  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(2,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  u_ct_val(1,:,i)=2.0*u_ct_val(2,:,i)-u_ct_val(3,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif












if (averaging.eq.1)then

 call averaging_t(n)

end if





end subroutine dual_time_ex_2d



















subroutine dual_time_2d(n)
!> @brief
!> dual time stepping 2d
implicit none
integer::i,j,k,kmaxe,nvar,jj,kill_nan_global
integer,intent(in)::n
#ifdef gpu
logical::kill_nan_reduce
real::allresdt_reduce
#endif
logical::bad_impdu,accepted_update
real::verysmall,update_alpha
real::firsti,resmaxi,rsumfacei,suml2ri,dummy3i,inner_tol
real,dimension(1:gpu_max_nvar)::candidate,old_state,limited
verysmall = tolsmall

inner_tol=reslimit

kmaxe=xmpielrank(n)


if (it.eq.restart)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)

  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(1,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)

  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if



firsti=0.0d0
do jj=1,upperlimit
      rsumfacei=zero;allresdt=zero;dummy3i=zero;
      if (jj.eq.1)then
    iscoun=1
      else
      iscoun=2
      end if

#ifdef gpu
!$omp target update to(iscoun)
#endif

call call_flux_subroutines_2d


if (relax.eq.3)then



!call relaxation_lumfree(n)


else



if (lowmemory.eq.0)then

 call relaxation2d(n)

 else

 call relaxation_lm2d(n)

 end if

end if

kill_nan=0
allresdt = 0.0d0
#ifdef gpu
kill_nan_reduce=.false.
allresdt_reduce = 0.0d0
#endif


#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(rsumfacei,j) reduction(+:allresdt_reduce) reduction(.or.:kill_nan_reduce)
#else
!$omp do private(rsumfacei,j) reduction(+:allresdt) reduction(max:kill_nan)
#endif
do i=1,kmaxe
      rsumfacei=sqrt(((impdu(i,1))**2)+((impdu(i,2))**2)+((impdu(i,3))**2)+((impdu(i,4))**2))
#ifdef gpu
      allresdt_reduce=allresdt_reduce+(rsumfacei*ielem_totvolume(i))
#else
      allresdt=allresdt+(rsumfacei*ielem_totvolume(i))
#endif


      if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)))then
#ifdef gpu
      kill_nan_reduce=.true.
#else
      kill_nan=1
#endif
      end if
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
          if (impdu(i,j).ne.impdu(i,j)) then
#ifdef gpu
            kill_nan_reduce=.true.
#else
            kill_nan=1
#endif
          end if
        end do
      end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
allresdt=allresdt_reduce
kill_nan=0
if (kill_nan_reduce) kill_nan=1
#else
!$omp end do
#endif


  !$omp barrier
    !$omp master
    kill_nan_global = 0
    call mpi_allreduce(kill_nan,kill_nan_global,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
    kill_nan = kill_nan_global
    if (kill_nan.eq.1)then
        kill=1
        if (n.eq.0)then
      print*,"killed due to divergence in implicit time stepping"
    end if
    end if
    !$omp end master
    !$omp barrier

    if (kill_nan.eq.1)then
        return
    end if



!$omp master
dummy3i=zero

call mpi_allreduce(allresdt,dummy3i,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
allresdt=dummy3i/totalvolume



if (allresdt.gt.firsti)then
firsti=allresdt
end if

allresdt=allresdt/firsti

if (n.eq.0)then
				  open(77,file='res1.txt',form='formatted',action='write',position='append')
				  write(77,*)allresdt,jj,it
				  close(77)

end if


!$omp end master
!$omp barrier




  if ((allresdt.le.inner_tol).or.(jj.eq.upperlimit))then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#else
!$omp do private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#endif
  do i=1,kmaxe
    bad_impdu=.false.
    do j=1,nof_variables
      if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
    end do
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
        if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
      end do
    end if
    if (bad_impdu)then
      impdu(i,1:nof_variables)=zero
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
      end if
    else
      old_state=zero
      candidate=zero
      old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      candidate(1:nof_variables)=old_state(1:nof_variables)+impdu(i,1:nof_variables)
      call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
      if (accepted_update) then
        u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
          update_alpha*impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
          do k=1,turbulenceequations+passivescalar
            if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
              u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
            end if
          end do
        end if
      else
        impdu(i,1:nof_variables)=zero
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
        end if
      end if
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
if (realgas.eq.1)then
if (rg_relax.eq.2)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
! Real-gas sources are assembled through rhs_val in sources_computation.
! Do not mutate u_c_val here with a local pseudo-time source step.
continue
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if
if (realgas.eq.1)call normalise_species(n)


  exit

else
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#else
!$omp do private(j,bad_impdu,candidate,old_state,limited,update_alpha,accepted_update)
#endif
 do i=1,kmaxe
    bad_impdu=.false.
    do j=1,nof_variables
      if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
    end do
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
      do j=nof_variables+1,nof_variables+turbulenceequations+passivescalar
        if (impdu(i,j).ne.impdu(i,j)) bad_impdu=.true.
      end do
    end if
    if (bad_impdu)then
      impdu(i,1:nof_variables)=zero
      if ((turbulence.gt.0).or.(passivescalar.gt.0))then
        impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
      end if
    else
      old_state=zero
      candidate=zero
      old_state(1:nof_variables)=u_c_val(1,1:nof_variables,i)
      candidate(1:nof_variables)=old_state(1:nof_variables)+impdu(i,1:nof_variables)
      call fv_limit_candidate_update(old_state,candidate,limited,update_alpha,accepted_update)
      if (accepted_update) then
        u_c_val(1,1:nof_variables,i)=limited(1:nof_variables)
        impdu(i,1:nof_variables)=limited(1:nof_variables)-old_state(1:nof_variables)
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=&
          update_alpha*impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)
          do k=1,turbulenceequations+passivescalar
            if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
              u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
            end if
          end do
        end if
      else
        impdu(i,1:nof_variables)=zero
        if ((turbulence.gt.0).or.(passivescalar.gt.0))then
          impdu(i,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=zero
        end if
      end if
    end if
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if (realgas.eq.1)then
if (rg_relax.eq.2)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
! Real-gas sources are assembled through rhs_val in sources_computation.
! Do not mutate u_c_val here with a local pseudo-time source step.
continue
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if
if (realgas.eq.1)call normalise_species(n)

end if


end do





#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  u_c_val(3,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  !u_c_val(1,1:nof_variables,i)=(2.0*u_c_val(2,1:nof_variables,i))-u_c_val(3,1:nof_variables,i)


  if ((turbulence.gt.0).or.(passivescalar.gt.0))then
  u_ct_val(3,:,i)=u_ct_val(2,:,i)
  u_ct_val(2,:,i)=u_ct_val(1,:,i)
  !u_ct_val(1,:,i)=2.0*u_ct_val(2,:,i)-u_ct_val(3,:,i)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





if (averaging.eq.1)then

 call averaging_t(n)

end if

end subroutine dual_time_2d




subroutine relaxation_ex(n)
implicit none
!> @brief
!> this subroutine solves the linear system for implicit time stepping either through matrix free lu-sgs low memory footprint
integer,intent(in)::n
integer::i,l,k,ii,sweeps,kmaxe,nvar,igoflux,icaseb,indt1,indt2,indt3,ijk,j
real::dt1,dtau



kmaxe=xmpielrank(n)



indt1=nof_variables+1
indt2=nof_variables+turbulenceequations+passivescalar
indt3=turbulenceequations+passivescalar




!
!

!===========================================================
! Loop 1: build impdu
!===========================================================
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(dt1,j)
#else
!$omp do private(dt1,j)
#endif
do i = 1, kmaxe

    dt1 = ielem_totvolume(i) / dt

    impdu(i,:) = zero

    do j = 1, nof_variables
        impdu(i,j) = dt1 * ( &
             1.5d0 * u_c_val(1,j,i) &
           - 2.0d0 * u_c_val(2,j,i) &
           + 0.5d0 * u_c_val(3,j,i) ) &
           + rhs_val(j,i)
    end do

    if ( (turbulence.gt.0) .or. (passivescalar.gt.0) ) then
        do j = 1, turbulenceequations+passivescalar
            impdu(i,nof_variables+j) = dt1 * ( &
                 1.5d0 * u_ct_val(1,j,i) &
               - 2.0d0 * u_ct_val(2,j,i) &
               + 0.5d0 * u_ct_val(3,j,i) ) &
               + rhst_val(j,i)
        end do
    end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


!===========================================================
! Loop 2: scale impdu by dtau
!===========================================================
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(dtau,j)
#else
!$omp do private(dtau,j)
#endif
do i = 1, kmaxe

    dtau = (ielem_dtl(i) / ielem_totvolume(i)) * &
           (1.0d0 / (1.0d0 + 1.5d0 * (ielem_dtl(i) / dt)))

    do j = 1, nof_variables
        impdu(i,j) = -impdu(i,j) * dtau
    end do

    if ( (turbulence.gt.0) .or. (passivescalar.gt.0) ) then
        do j = 1, turbulenceequations+passivescalar
            impdu(i,nof_variables+j) = -impdu(i,nof_variables+j) * dtau
        end do
    end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end subroutine relaxation_ex





subroutine averaging_t(n)
!> @brief
!> temporal averaging for unsteady simulations
implicit none
integer,intent(in)::n
integer::i,kmaxe,nvar

kmaxe=xmpielrank(n)



if (dimensiona.eq.3)then





  if (tz1.gt.zero)then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(nvar)
#else
!$omp do
#endif
  do i=1,kmaxe
    u_c_val(ind1,:,i)=(((tz1-dt)/(tz1))*u_c_val(ind1,:,i))+((dt*u_c_val(1,:,i))/tz1)
      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
    do nvar=1,turbulenceequations+passivescalar

    u_ct_val(ind1,nvar,i)=(((tz1-dt)/(tz1))*u_ct_val(ind1,nvar,i))+((dt*u_ct_val(1,nvar,i))/(tz1*u_c_val(ind1,1,i)))

    end do
    end if
    !u,v,w,uv,uw,wv,ps
    u_c_rms(1,i)=sqrt(abs(((u_c_rms(1,i)**2)*((tz1-dt)/(tz1)))+(((u_c_val(1,2,i)/u_c_val(1,1,i)-u_c_val(ind1,2,i)/u_c_val(ind1,1,i))**2)*dt/tz1)))
    u_c_rms(2,i)=sqrt(abs(((u_c_rms(2,i)**2)*((tz1-dt)/(tz1)))+(((u_c_val(1,3,i)/u_c_val(1,1,i)-u_c_val(ind1,3,i)/u_c_val(ind1,1,i))**2)*dt/tz1)))
    u_c_rms(3,i)=sqrt(abs(((u_c_rms(3,i)**2)*((tz1-dt)/(tz1)))+(((u_c_val(1,4,i)/u_c_val(1,1,i)-u_c_val(ind1,4,i)/u_c_val(ind1,1,i))**2)*dt/tz1)))





    u_c_rms(4,i)=(((u_c_rms(4,i))*((tz1-dt)/(tz1)))+&
(((((u_c_val(1,2,i)/u_c_val(1,1,i))-(u_c_val(ind1,2,i)/u_c_val(ind1,1,i)))*((u_c_val(1,3,i)/u_c_val(1,1,i))-(u_c_val(ind1,3,i)/u_c_val(ind1,1,i)))))*dt/tz1))
    u_c_rms(5,i)=(((u_c_rms(5,i))*((tz1-dt)/(tz1)))+&
(((((u_c_val(1,2,i)/u_c_val(1,1,i))-(u_c_val(ind1,2,i)/u_c_val(ind1,1,i)))*((u_c_val(1,4,i)/u_c_val(1,1,i))-(u_c_val(ind1,4,i)/u_c_val(ind1,1,i)))))*dt/tz1))
    u_c_rms(6,i)=(((u_c_rms(6,i))*((tz1-dt)/(tz1)))+&
(((((u_c_val(1,3,i)/u_c_val(1,1,i))-(u_c_val(ind1,3,i)/u_c_val(ind1,1,i)))*((u_c_val(1,4,i)/u_c_val(1,1,i))-(u_c_val(ind1,4,i)/u_c_val(ind1,1,i)))))*dt/tz1))
    if ((passivescalar.gt.0))then
    u_c_rms(7,i)=sqrt(abs(((u_c_rms(7,i)**2)*((tz1-dt)/(tz1)))+(((u_ct_val(1,turbulenceequations+1,i)&
-u_ct_val(ind1,turbulenceequations+1,i))**2)*dt/tz1)))
    end if




  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  else
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
  do i=1,kmaxe
    u_c_val(ind1,:,i)=zero;u_c_rms(:,i)=zero
    if ((turbulence.eq.1).or.(passivescalar.gt.0))then


    u_ct_val(ind1,:,i)=zero


    end if

  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  end if

else
if (t.gt.0.0)then
#ifdef gpu
!$omp target teams distribute parallel do &
!$omp& private(nvar)
#else
!$omp  do
#endif
  do i=1,kmaxe
    u_c_val(ind1,:,i)=(((tz1-dt)/(tz1))*u_c_val(ind1,:,i))+((dt*u_c_val(1,:,i))/tz1)
      if ((turbulence.eq.1).or.(passivescalar.gt.0))then
    do nvar=1,turbulenceequations+passivescalar

    u_ct_val(ind1,nvar,i)=(((tz1-dt)/(tz1))*u_ct_val(ind1,nvar,i))+((dt*u_ct_val(1,nvar,i))/(tz1*u_c_val(ind1,1,i)))

    end do
    end if
    !u,v,uv,ps
    u_c_rms(1,i)=sqrt(abs(((u_c_rms(1,i)**2)*((tz1-dt)/(tz1)))+(((u_c_val(1,2,i)-u_c_val(ind1,2,i))**2)*dt/tz1)))/u_c_val(ind1,1,i)
    u_c_rms(2,i)=sqrt(abs(((u_c_rms(2,i)**2)*((tz1-dt)/(tz1)))+(((u_c_val(1,3,i)-u_c_val(ind1,3,i))**2)*dt/tz1)))/u_c_val(ind1,1,i)


    u_c_rms(3,i)=(((u_c_rms(4,i))*((tz1-dt)/(tz1)))+&
((((u_c_val(1,2,i)-u_c_val(ind1,2,i))*(u_c_val(1,3,i)-u_c_val(ind1,3,i))))*dt/tz1))/u_c_val(ind1,1,i)

    if ((passivescalar.gt.0))then
    u_c_rms(4,i)=sqrt(abs(((u_c_rms(4,i)**2)*((tz1-dt)/(tz1)))+(((u_ct_val(1,turbulenceequations+1,i)&
-u_ct_val(ind1,turbulenceequations+1,i))**2)*dt/tz1)))/u_c_val(ind1,1,i)
    end if

  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  else
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp  do
#endif
  do i=1,kmaxe
    u_c_val(ind1,:,i)=zero
    if ((turbulence.eq.1).or.(passivescalar.gt.0))then


    u_ct_val(ind1,:,i)=zero


    end if

  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
  end if
 end if



   if (outsurf.eq.1)then
   call exchange_higher_av(n)
   call average_stresses(n)
   end if

end subroutine averaging_t


! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !---------------------------------------------------------------------------------------------!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!subroutine called to advance solution by one time step size!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine time_marching(n)
!> @brief
!> time marching subroutine 3d
implicit none
integer,intent(in)::n
real,dimension(1:5)::dummyout,dummyin
integer::i,kmaxe,ttime,kill_local
logical::tgv_energy_prev_set
real::dtiv
real::cput1,cput2,cput3,cput4,cput5,cput6,cput8,timec3,timec1,timec4,timec8,totv1,totv2,dumetg1,dumetg2,tzx1,tzx2,resolx,totens1,totens2,totensx1,totensx2
#if defined(gpu) || defined(xpu)
real::totk_reduce,totens_reduce,totensx_reduce
#endif
      kill=0
      t=res_time
      resolx=0.01
      iscoun=1
      kmaxe=xmpielrank(n)
      every_time=((idnint(t/output_freq)) * output_freq)+output_freq

      totv1=0.0
      tgv_energy_prev_set=.false.



#if defined(gpu) || defined(xpu)
!$omp target update to(kill,iscoun)
#endif

!$omp barrier
!$omp master
    if (initcond.eq.95)then
    call checkpointv3(n)
    end if
	cput1=cpux1(1)
	cput4=cpux1(1)
	cput5=cpux1(1)
	cput8=cpux1(1)
!$omp end master
!$omp barrier



	it=restart
	if (dg.eq.1)call sol_integ_dg_init(n)

!$omp barrier
!$omp master
      if (tecplot.lt.5)then
        call grid_write
        if (outsurf.eq.1)then
        call surf_write
        end if
      end if
      if ((average_restart.eq.0).and.(averaging.eq.1)) then
				tz1=0.0d0
				else
				tz1=t
		end if


!$omp end master
!$omp barrier



!$omp barrier
!$omp master
	call volume_solution_write
	if (outsurf.eq.1)then
	call surface_solution_write
	end if
!$omp end master
!$omp barrier




      if ((it.eq.0).and.(initcond.eq.95))then

        call exchange_higher(n)

        call arbitrary_order(n)

#if defined(gpu) || defined(xpu)
!$omp barrier
!$omp master
#endif

        call enstrophy_calc(n)

#if defined(gpu) || defined(xpu)
!$omp end master
!$omp barrier
#endif


      end if

     do

		    call calculate_cfl(n)

		    if (rungekutta.ge.5)call calculate_cfll(n)


		    if (dg.eq.1)then
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp  do
#endif
              do i=1,kmaxe
              ielem_condition(i)=0
              ielem_troubled(i)=0
              end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
            end if




			!$omp barrier
			!$omp master
			dummyout(1)=dt
			cput2=mpi_wtime()
			timec8=cput2-cput8
			timec1=cput2-cput1
		        dummyout(2)=timec1
			dummyin=0.0
			timec3=cput2-cput4
			dummyout(3)=timec3
			timec4=cput2-cput5
			dummyout(4)=timec4
			dummyout(5)=timec8

			call mpi_allreduce(dummyout,dummyin,5,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
			dtiv=dummyin(1)
			dt=dummyin(1)
			timec1=dummyin(2)
			timec3=dummyin(3)
			timec4=dummyin(4)
			timec8=dummyin(5)
				!$omp end master
				!$omp barrier

                              if (initcond.eq.95)then
#if defined(gpu) || defined(xpu)
                          totk_reduce=0.0d0
                          totens_reduce=0.0d0
                          totensx_reduce=0.0d0
!$omp target teams distribute parallel do &
!$omp& reduction(+:totk_reduce,totens_reduce,totensx_reduce) &
!$omp& map(alloc: ielem_totvolume, ielem_vortex, u_c_val, xmpielrank) &
!$omp& firstprivate(boundtype)
#else
                          !$omp master
                          totk=0.0d0
                          totens=0.0d0
                          totensx=0.0d0
                          !$omp end master
                          !$omp barrier
                          !$omp do reduction(+:totk,totens,totensx)
#endif

                          do i=1,xmpielrank(n)

#if defined(gpu) || defined(xpu)
                              totk_reduce=totk_reduce+ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
          (((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2)+((u_c_val(1,4,i)/u_c_val(1,1,i))**2))
#else
                              totk=totk+ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
          (((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2)+((u_c_val(1,4,i)/u_c_val(1,1,i))**2))
#endif

                              if (boundtype.eq.1)then

#if defined(gpu) || defined(xpu)
                              totens_reduce=totens_reduce+(ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
                              ielem_vortex(2,i))
#else
                              totens=totens+(ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
                              ielem_vortex(2,i))
#endif
                              else

#if defined(gpu) || defined(xpu)
                              totens_reduce=totens_reduce+(ielem_vortex(2,i))
                              totensx_reduce=totensx_reduce+(ielem_vortex(3,i))
#else
                              totens=totens+(ielem_vortex(2,i))
                              totensx=totensx+(ielem_vortex(3,i))
#endif
                              end if

                          end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
                          totk=totk_reduce
                          totens=totens_reduce
                          totensx=totensx_reduce
#else
                          !$omp end do
#endif

                         !$omp barrier
                          !$omp master
                          dumetg1=totk
                          dumetg2=0.0
                          call mpi_barrier(mpi_comm_world,ierror)
                          call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
                          totk=dumetg2
                          dumetg1=totens
                          dumetg2=0.0
                          call mpi_barrier(mpi_comm_world,ierror)
                          call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
                          totens=dumetg2
                          dumetg1=totensx
                          dumetg2=0.0
                          call mpi_barrier(mpi_comm_world,ierror)
                          call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
                          totensx=dumetg2
                          if (n.eq.0)then
                          if (.not.tgv_energy_prev_set)then
                          totv1=totk/((2.0*pi)**3)
                          tgv_energy_prev_set=.true.
                          end if
                          totens1=totens/(((2.0*pi)**3))
                          totensx1=totensx/(((2.0*pi)**3))
                              if (it.eq.0)then
                              taylor=totk
                              taylor_ens=totens
                              taylor_ensx=totensx
                              end if



                          end if
                        !$omp end master
			!$omp barrier





			    end if
               !$omp barrier
			!$omp master
			if (rungekutta.ge.11)then
			dt=timestep
			if (initcond.eq.95)then
			dt=min(dt,out_time-t,every_time-t)
			else
			dt=min(dt,out_time-t,every_time-t)
			end if
			else
			if (initcond.eq.95)then
			dt=min(dt,out_time-t,every_time-t)
			else
			dt=min(dt,out_time-t,every_time-t)
			end if
			end if

#if defined(gpu) || defined(xpu)
!$omp target update to(dt)
#endif
				  if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='unknown',action='write',position='append')
				  write(63,'(A,I10,A,ES20.10,A,ES20.10,A,ES20.10)') 'time_step it=', it, ' dt_cfl=', dtiv, ' dt_used=', dt, ' t=', t
				  flush(63)
				  close(63)
				  end if
				if (dg.eq.1)then
			if (filtering.gt.0)then
			call apply_filter_dg(n)
			end if
			end if



			!$omp end master
			!$omp barrier
			select case(rungekutta)

			case(1)
			call runge_kutta1(n)

			case(2)
			call runge_kutta2(n)

			case(3)

			if (mood.eq.1)then
			call runge_kutta3_mood(n)
			else
			call runge_kutta3(n)
			end if

			case(4)
			call runge_kutta4(n)

			case(5)
			call runge_kutta5(n)

			case(10)
			call implicit_times(n)


			case(11)
			call dual_time(n)


			case(12)
			call dual_time_ex(n)


			end select


			if (dg.eq.1)call sol_integ_dg(n)
            if (realgas.eq.1) call normalise_species(n)
			!$omp barrier
			!$omp master


			if (rungekutta.ge.11)then
			 t=t+(dt)
			  tz1=tz1+(dt)

				else
	                       t=t+dt
				tz1=tz1+dt
				  end if
				if ((mach_outlet_target.gt.0.0d0).and.(mach_outlet_update_freq.gt.0))then
				  if (mod(it,mach_outlet_update_freq).eq.0) call update_outlet_mach_controller(n)
				end if
				 !$omp end master
				!$omp barrier

#if defined(gpu) || defined(xpu)
!$omp target update to(tz1,t)
#endif


			!$omp barrier

				if (dg.eq.1)then
                      if (code_profile.ne.102)then
                          if ( mod(it, 100) .eq. 0) then
                            call troubled_history
                          end if
                      end if
                      if ( filtering .eq. 1) then
                        if ( mod(it, 100) .eq. 0) then
                          call filtered_history
                        end if
                      end if
                end if

          if ( mod(it, 100) .eq. 0) then
            call reduced_history
          end if


			!$omp barrier


			if (initcond.eq.95)then

			    call exchange_higher(n)
			    call arbitrary_order(n)

#if defined(gpu) || defined(xpu)
!$omp barrier
!$omp master
#endif
			    call enstrophy_calc(n)
#if defined(gpu) || defined(xpu)
!$omp end master
!$omp barrier
#endif

#if defined(gpu) || defined(xpu)
                          totk_reduce=0.0d0
                          totens_reduce=0.0d0
                          totensx_reduce=0.0d0
!$omp target teams distribute parallel do &
!$omp& reduction(+:totk_reduce,totens_reduce,totensx_reduce) &
!$omp& map(alloc: ielem_totvolume, ielem_vortex, u_c_val, xmpielrank) &
!$omp& firstprivate(boundtype)
#else
                          !$omp master
                          totk=0.0d0
                          totens=0.0d0
                          totensx=0.0d0
                          !$omp end master
                          !$omp barrier
                          !$omp  do reduction(+:totk,totens,totensx)
#endif
				do i=1,xmpielrank(n)


#if defined(gpu) || defined(xpu)
                    totk_reduce=totk_reduce+ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
(((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2)+((u_c_val(1,4,i)/u_c_val(1,1,i))**2))
#else
                    totk=totk+ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
(((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2)+((u_c_val(1,4,i)/u_c_val(1,1,i))**2))
#endif






                              if (boundtype.eq.1)then

#if defined(gpu) || defined(xpu)
                              totens_reduce=totens_reduce+(ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
                              ielem_vortex(2,i))
#else
                              totens=totens+(ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
                              ielem_vortex(2,i))
#endif
                              else

#if defined(gpu) || defined(xpu)
                              totens_reduce=totens_reduce+(ielem_vortex(2,i))
                               totensx_reduce=totensx_reduce+(ielem_vortex(3,i))
#else
                              totens=totens+(ielem_vortex(2,i))
                               totensx=totensx+(ielem_vortex(3,i))
#endif
                              end if





				end do
#if defined(gpu) || defined(xpu)
!$omp end target teams distribute parallel do
                          totk=totk_reduce
                          totens=totens_reduce
                          totensx=totensx_reduce
#else
                          !$omp end do
#endif

             !$omp barrier
			!$omp master
				dumetg1=totk
				dumetg2=0.0
				call mpi_barrier(mpi_comm_world,ierror)
				call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
				totk=dumetg2

				dumetg1=totens
				dumetg2=0.0
				call mpi_barrier(mpi_comm_world,ierror)
				call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
				totens=dumetg2



				dumetg1=totensx
				dumetg2=0.0
				call mpi_barrier(mpi_comm_world,ierror)
				call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
				totensx=dumetg2

				if (n.eq.0)then
                          totv2=totk/((2.0*pi)**3)
                          totens2=totens/(((2.0*pi)**3))
                          totensx2=totensx/(((2.0*pi)**3))
                              if (it.eq.0)then
                              taylor=totk
                              taylor_ens=totens
                              taylor_ensx=totensx
                              end if

                            if (it.eq.0)then
                            open(73,file='energy.dat',form='formatted',status='new',action='write',position='append')
                            else
                            open(73,file='energy.dat',form='formatted',status='old',action='write',position='append')
                            end if
                            if (dg.eq.1)then
                            write(73,'(e14.7,1x,e14.7,1x,e14.7)')t,totk/taylor,-(totv2-totv1)/dt
                            else
                              if (boundtype.eq.1)then
                              write(73,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,totk/taylor,-(totv2-totv1)/dt,totens/taylor_ens
                              else
                              write(73,'(e14.7,1x,e14.7,1x,e14.7,1x,e14.7,1x,e14.7)')t,totv2,-(totv2-totv1)/dt,totens2,totensx2
                              end if
                            end if
                          close(73)
                          totv1=totv2
                          tgv_energy_prev_set=.true.
				end if



			call mpi_barrier(mpi_comm_world,ierror)
			!$omp end master
			!$omp barrier


          if (adda.eq.1) then

          totk = 0.0d0

#if defined(gpu) || defined(xpu)
        totk_reduce = 0.0d0
!$omp target teams distribute parallel do &
!$omp& reduction(+:totk_reduce) &
!$omp& map(alloc: ielem_er, xmpielrank) &
!$omp& firstprivate(n)
          do i = 1, xmpielrank(n)
              totk_reduce = totk_reduce + ielem_er(i)
          end do
!$omp end target teams distribute parallel do
        totk = totk_reduce
#else
        !$omp do reduction(+:totk)
          do i = 1, xmpielrank(n)
              totk = totk + ielem_er(i)
          end do
        !$omp end do
#endif

         !$omp barrier
		!$omp master
          dumetg1 = totk
          dumetg2 = 0.0d0
          call mpi_allreduce(dumetg1, dumetg2, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierror)

          if (n.eq.0) then
              totk = dumetg2 / imaxe

              if (it.eq.0) then
                open(123,file='er.dat',form='formatted',status='new',action='write',position='append')
              else
                open(123,file='er.dat',form='formatted',status='old',action='write',position='append')
              end if

              write(123,*) t, totk
              close(123)
          end if
	!$omp end master
			!$omp barrier

        end if




! 				end if
			    end if




	!$omp barrier
			!$omp master
			    if ((initcond.eq.405).or.(initcond.eq.422).or.(initcond.eq.411).or.(initcond.eq.157))then
			    if ( mod(it, 100) .eq. 0)then
#if defined(gpu) || defined(xpu)
!$omp target update from(u_c_val)
#endif
                                 call trajectories
			    end if
			    end if



			!$omp end master
			!$omp barrier




			if ( mod(it, iforce) .eq. 0) then
			if (outsurf.eq.1) then

				  call forces
			end if
			end if



			if ((rungekutta.ge.5).and.(rungekutta.lt.11))then
			if ( mod(it, residualfreq) .eq. 0) then

                               call residual_compute
                               call update_cfl_ramp(n)
			end if
			end if





			!$omp master
			if (nprobes.gt.0) then
			if ( mod(it, 100) .eq. 0) then
#if defined(gpu) || defined(xpu)


!$omp target update from(u_c_val)
#endif
			call probing

			end if
			end if



			if (timec1.ge.ievery)then

			    call volume_solution_write
			     if (outsurf.eq.1)then

			    call surface_solution_write
			    end if
			cput1=mpi_wtime()
			end if





			if (initcond.eq.95)then
                  if (abs(t - ((idnint(t/output_freq)) * output_freq)).le.tolsmall) then

                      call volume_solution_write
                      if (outsurf.eq.1)then
                      call surface_solution_write
                      end if
                      if (initcond.eq.95)then
                      call checkpointv4(n)
                      end if
                  every_time=every_time+output_freq
                  end if

            else


                  if ((code_profile.lt.0).or.(code_profile.eq.100).or.(code_profile.eq.101).or.(code_profile.eq.102))then
                        if (abs(t - ((idnint(t/output_freq)) * output_freq)).le.tolsmall) then
                            call volume_solution_write
                            if (outsurf.eq.1)then
                            call surface_solution_write
                            end if
                        every_time=every_time+output_freq
                        end if
                  end if

            end if






                if (timec8.ge.ieveryav)then
                      if (averaging.eq.1)then
                      call volume_solution_write_av
                        if (outsurf.eq.1)then
                        call surface_solution_write_av
                        end if
                      end if
                cput8=mpi_wtime()
                end if



                  if (timec4.ge.ievery2)then
                      call checkpointing
                        if (averaging.eq.1)then
                      call checkpointing_av
                        end if

                    cput5=mpi_wtime()

                  end if


			!$omp end master
			!$omp barrier


            !$omp barrier
			!$omp master

			it=it+1

			if ((it.eq.ntmax).or.(timec3.ge.wallc).or.(dtiv.gt.out_time))then
			 kill=1
			end if

			if ((rungekutta.lt.5).or.(rungekutta.ge.11))then
			if ((t.ge.out_time).or.(dtiv.gt.out_time))then
			kill=1
			end if
			end if
			kill_local=kill
			call mpi_allreduce(kill_local,kill,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
#if defined(gpu) || defined(xpu)
!$omp target update to(kill)
#endif

			!$omp end master
			!$omp barrier





			!$omp master
			if (kill.eq.1)then

			    call volume_solution_write
			      if (outsurf.eq.1)then
			      call surface_solution_write
			      end if
			    call checkpointing
			      if (averaging.eq.1)then
			      call volume_solution_write_av
				if (outsurf.eq.1)then
				  call surface_solution_write_av
				end if
				call checkpointing_av
			      end if
			end if

			!$omp end master
			!$omp barrier

			if (kill.eq.1)then
			if (itestcase.le.3)then

			call calculate_error(n)
			end if

			return
			end if



end do

end subroutine time_marching



subroutine time_marching2(n)
!> @brief
!> time marching subroutine 2d
implicit none
integer,intent(in)::n
real,dimension(1:5)::dummyout,dummyin
integer::i,kmaxe,kill_local
real::cput1,cput2,cput3,cput4,cput5,cput6,cput8,timec3,timec1,timec4,timec8,totv1,totv2,dumetg1,dumetg2
real::dtiv,flort
kmaxe=xmpielrank(n)
kill=0
t=res_time
iscoun=1

every_time=((idnint(t/output_freq)) * output_freq)+output_freq


#if defined(gpu) || defined(xpu)
!$omp target update to(kill,iscoun)
#endif

!$omp master
cput1=cpux1(1)
cput4=cpux1(1)
cput5=cpux1(1)
cput8=cpux1(1)
!$omp end master
!$omp barrier


it=restart
if (dg.eq.1)call sol_integ_dg_init(n)
!$omp barrier
!$omp master
if (tecplot.lt.5)then
  call grid_write
end if


call volume_solution_write
if (outsurf.eq.1)then
    call surf_write
end if

if ((average_restart.eq.0).and.(averaging.eq.1)) then
    tz1=0.0
else
    tz1=t
end if
!$omp end master
!$omp barrier



if (initcond.gt.100000)then
it=100000
call exchange_higher(n)
call arbitrary_order(n)
call writesols(n)
call volume_solution_write


else




do




    call calculate_cfl2d(n)

    if (rungekutta.ge.5) call calculate_cfll2d(n)

    if (dg.eq.1)then
        do i=1,kmaxe
        ielem_condition(i)=0
        ielem_troubled(i)=0
        end do
    end if




    !$omp master
    dummyout(1)=dt
    cput2=mpi_wtime()
    timec8=cput2-cput8
    timec1=cput2-cput1
    dummyout(2)=timec1
    dummyin=0.0d0
    timec3=cput2-cput4
    dummyout(3)=timec3
    timec4=cput2-cput5
    dummyout(4)=timec4
    dummyout(5)=timec8

    call mpi_allreduce(dummyout,dummyin,5,mpi_double_precision,mpi_min,mpi_comm_world,ierror)
    dtiv=dummyin(1)
    dt=dummyin(1)
    timec1=dummyin(2)
    timec3=dummyin(3)
    timec4=dummyin(4)
    timec8=dummyin(5)
                  if ( mod(it, 10) .eq. 0) then
				   if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				 write(63,'(A,I10,A,ES20.10,A,ES20.10)') 'it=', it, ' dt=', dt, ' t=', t
				  close(63)
				  end if
				  end if




     if ((multispecies.eq.1))then
         if((initcond.eq.405).or.(initcond.eq.411))then
#if defined(gpu) || defined(xpu)
!$omp target update from(u_c_val)
#endif
                 call trajectories
            ! end if
         end if
     end if

    if (rungekutta.ge.11)then
        dt=timestep
        dt=min(dt,out_time-t,every_time-t)
    else
        dt=min(dt,out_time-t,every_time-t)
    end if


    !$omp end master
    !$omp barrier

#if defined(gpu) || defined(xpu)
!$omp target update to(dt)
#endif

                select case(rungekutta)

                    case(1)
                        call runge_kutta1_2d(n)

                    case(2)
                        call runge_kutta2_2d(n)

                    case(3)
                        if (mood.eq.1)then
                            call runge_kutta3_2d_mood(n)
                        else
                            call runge_kutta3_2d(n)
                        end if

                    case(4)
                        call runge_kutta4_2d(n)

                    case(5)
                        call runge_kutta5_2d(n)

                    case(10)
                        call implicit_times_2d(n)

                    case(11)
                        call dual_time_2d(n)

                    case(12)
                        call dual_time_ex_2d(n)

                end select

          if (dg.eq.1)call sol_integ_dg(n)

	          ! Real-gas sources are handled by the residual assembly for all
	          ! temporal schemes; avoid a separate in-place source update here.


           if (realgas.eq.1) then
!            call diagnose_species_before_norm(IT)
           call normalise_species(n)
           end if


! increment time
    !$omp barrier
    !$omp master
    if (rungekutta.ge.11)then
        t=t+(dt)
        tz1=tz1+(dt)
	    else
	        t=t+dt
	        tz1=tz1+dt

	end if
	    if ((mach_outlet_target.gt.0.0d0).and.(mach_outlet_update_freq.gt.0))then
	        if (mod(it,mach_outlet_update_freq).eq.0) call update_outlet_mach_controller(n)
	    end if


	      !$omp end master
    !$omp barrier

#if defined(gpu) || defined(xpu)
!$omp target update to(tz1,t)
#endif


          !$omp barrier
          if (dg.eq.1)then
          if (code_profile.ne.102)then
          if ( mod(it, 100) .eq. 0) then
            call troubled_history
          end if
          end if
          end if

           if ( mod(it, 20) .eq. 0) then
            call reduced_history
          end if

          if (mood.gt.0)then
          call troubled_history
          end if

        !$omp barrier


! write output


    if ( mod(it, iforce) .eq. 0) then
        if (outsurf.eq.1) then

            call forces
        end if
    end if

    if ((rungekutta.ge.5).and.(rungekutta.lt.11))then
        if ( mod(it, residualfreq) .eq. 0) then

            call residual_compute
            call update_cfl_ramp(n)
        end if
    end if


    !$omp master
    if (nprobes.gt.0) call probing2d

    if (timec1.ge.ievery)then
!
        call volume_solution_write
        if (outsurf.eq.1)then
            call surface_solution_write
        end if
    cput1=mpi_wtime()
    end if


    if (timec8.ge.ieveryav)then
        if (averaging.eq.1)then
            call volume_solution_write_av
            if (outsurf.eq.1)then
                call surface_solution_write_av
            end if
        end if
        cput8=mpi_wtime()
    end if

           if ((code_profile.lt.0).or.(code_profile.eq.100).or.(code_profile.eq.101).or.(code_profile.eq.102))then



			if (abs(t - ((idnint(t/output_freq)) * output_freq)).le.tolsmall) then

                call volume_solution_write
			     if (outsurf.eq.1)then
			    call surface_solution_write
			    end if
			every_time=every_time+output_freq
            end if
            end if

        if (timec4.ge.ievery2)then
			    call checkpointing
			      if (averaging.eq.1)then
				call checkpointing_av
			      end if

			  cput5=mpi_wtime()

			end if


! check end condition

    !$omp end master
    !$omp barrier



          !$omp master
			it=it+1

			if ((it.eq.ntmax).or.(timec3.ge.wallc).or.(dtiv.gt.out_time))then
			 kill=1
			end if

			if ((rungekutta.lt.5).or.(rungekutta.ge.11))then
			if ((t.ge.out_time).or.(dtiv.gt.out_time))then
			kill=1
			end if
			end if
			kill_local=kill
			call mpi_allreduce(kill_local,kill,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
#if defined(gpu) || defined(xpu)
!$omp target update to(kill)
#endif
			!$omp end master
			!$omp barrier


        !$omp master
			if (kill.eq.1)then

			    call volume_solution_write
			      if (outsurf.eq.1)then
			      call surface_solution_write
			      end if
			    call checkpointing
			      if (averaging.eq.1)then
			      call volume_solution_write_av
				if (outsurf.eq.1)then
				  call surface_solution_write_av
				end if
				call checkpointing_av
			      end if
			end if

			!$omp end master
			!$omp barrier

			if (kill.eq.1)then
			if (itestcase.le.3)then

			call calculate_error(n)
			end if

			return
			end if


            !$omp barrier



end do


end if

end subroutine time_marching2


subroutine limit_realgas_source_timestep(n)
!> Restrict local pseudo-time by real-gas source stiffness.
implicit none
integer,intent(in)::n
integer::i,kmaxe
real::dt_source

if ((realgas.ne.1).or.(rg_relax.ne.2).or.(nof_species.le.0)) return

#if defined(xpu) && !defined(gpu)
return
#endif

kmaxe=xmpielrank(n)
if (kmaxe.le.0) return

#ifdef gpu
!$omp target teams distribute parallel do firstprivate(n,kmaxe) private(i,dt_source)
#else
!$omp barrier
!$omp do private(i,dt_source)
#endif
do i=1,kmaxe
  call realgas_source_timestep_cell(n,i,dt_source)
  if ((dt_source.eq.dt_source).and.(dt_source.gt.zero)) then
    ielem_dtl(i)=min(ielem_dtl(i),dt_source)
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine limit_realgas_source_timestep


subroutine realgas_source_timestep_cell(n,iconsidered,dt_source)
!> Estimate a local source timescale from chemistry and vibrational relaxation.
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::n,iconsidered
real,intent(out)::dt_source
integer::idxe,idxev,idx,k
real::rho,rhoe,rhoev,qscale,rate,sum_species_rate,dt_candidate
real,dimension(1:gpu_max_nvar)::source_r,src_jac_diag

dt_source=max(ielem_dtl(iconsidered),realgas_source_dtl_min_abs)

if ((realgas.ne.1).or.(rg_relax.ne.2).or.(nof_species.le.0)) return
if ((dimensiona+3+nof_species).gt.nof_variables) return

source_r=zero
src_jac_diag=zero
call sources_realgas(n,iconsidered,source_r,src_jac_diag,realgas_source_dtl_min_abs)

idxe=dimensiona+2
idxev=dimensiona+3
rho=max(abs(u_c_val(1,1,iconsidered)),realgas_source_dtl_min_abs)
rhoe=max(abs(u_c_val(1,idxe,iconsidered)),realgas_source_dtl_min_abs)
rhoev=abs(u_c_val(1,idxev,iconsidered))
sum_species_rate=zero

do k=1,nof_species
  idx=idxev+k
  rate=abs(source_r(idx))
  if ((rate.eq.rate).and.(rate.gt.1.0d-300)) then
    sum_species_rate=sum_species_rate+rate
  end if

  if ((source_r(idx).eq.source_r(idx)).and.(source_r(idx).lt.zero)) then
    qscale=max(abs(u_c_val(1,idx,iconsidered)), &
               realgas_source_species_floor_fraction*rho, &
               realgas_source_dtl_min_abs)
    dt_candidate=realgas_source_dtl_cfl*qscale/max(-source_r(idx),1.0d-300)
    if ((dt_candidate.eq.dt_candidate).and.(dt_candidate.gt.zero)) then
      dt_source=min(dt_source,dt_candidate)
    end if
  end if

  if ((src_jac_diag(idx).eq.src_jac_diag(idx)).and.(abs(src_jac_diag(idx)).gt.1.0d-300)) then
    dt_candidate=realgas_source_dtl_cfl/abs(src_jac_diag(idx))
    if ((dt_candidate.eq.dt_candidate).and.(dt_candidate.gt.zero)) then
      dt_source=min(dt_source,dt_candidate)
    end if
  end if
end do

if (sum_species_rate.gt.1.0d-300) then
  dt_candidate=realgas_source_dtl_cfl*rho/sum_species_rate
  if ((dt_candidate.eq.dt_candidate).and.(dt_candidate.gt.zero)) then
    dt_source=min(dt_source,dt_candidate)
  end if
end if

rate=abs(source_r(idxev))
if ((rate.eq.rate).and.(rate.gt.1.0d-300)) then
  qscale=max(rhoev, &
             realgas_source_energy_floor_fraction*rhoe, &
             realgas_source_dtl_min_abs)
  dt_candidate=realgas_source_dtl_cfl*qscale/rate
  if ((dt_candidate.eq.dt_candidate).and.(dt_candidate.gt.zero)) then
    dt_source=min(dt_source,dt_candidate)
  end if
end if

if ((src_jac_diag(idxev).eq.src_jac_diag(idxev)).and.(abs(src_jac_diag(idxev)).gt.1.0d-300)) then
  dt_candidate=realgas_source_dtl_cfl/abs(src_jac_diag(idxev))
  if ((dt_candidate.eq.dt_candidate).and.(dt_candidate.gt.zero)) then
    dt_source=min(dt_source,dt_candidate)
  end if
end if

dt_source=max(dt_source,realgas_source_dtl_min_abs)

end subroutine realgas_source_timestep_cell


subroutine fv_explicit_local_candidate(iconsidered,qold,qcandidate)
!> Fallback local explicit FV update when an implicit correction is unusable.
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::iconsidered
real,dimension(1:gpu_max_nvar),intent(in)::qold
real,dimension(1:gpu_max_nvar),intent(out)::qcandidate
integer::iv
real::oovolume

qcandidate=zero
qcandidate(1:nof_variables)=qold(1:nof_variables)

if (ielem_totvolume(iconsidered).le.fv_update_min_abs) return
oovolume=1.0d0/ielem_totvolume(iconsidered)

do iv=1,nof_variables
  qcandidate(iv)=qold(iv)-ielem_dtl(iconsidered)*rhs_val(iv,iconsidered)*oovolume
end do

end subroutine fv_explicit_local_candidate


subroutine fv_prepare_realgas_update_state(qin,qout,projected_ok)
!> Repair the species closure of a trial FV state before admissibility checks.
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar),intent(in)::qin
real,dimension(1:gpu_max_nvar),intent(out)::qout
logical,intent(out)::projected_ok
integer::k,idx
real::rho,species_sum,rhoy,scale

projected_ok=.true.
qout=zero
qout(1:nof_variables)=qin(1:nof_variables)

if ((realgas.ne.1).or.(nof_species.le.0)) return

if ((dimensiona+3+nof_species).gt.nof_variables) then
  projected_ok=.false.
  return
end if

rho=qout(1)
if ((rho.ne.rho).or.(abs(rho).gt.fv_update_max_abs)) then
  projected_ok=.false.
  return
end if
if (rho.le.fv_update_min_abs) then
  projected_ok=.false.
  return
end if

species_sum=zero
do k=1,nof_species
  idx=dimensiona+3+k
  rhoy=qout(idx)
  if ((rhoy.ne.rhoy).or.(abs(rhoy).gt.fv_update_max_abs)) then
    projected_ok=.false.
    return
  end if
  if (rhoy.lt.0.0d0) then
    projected_ok=.false.
    return
  end if
  species_sum=species_sum+rhoy
end do

if (species_sum.le.fv_update_min_abs) then
  projected_ok=.false.
  return
end if

scale=rho/species_sum
do k=1,nof_species
  idx=dimensiona+3+k
  qout(idx)=qout(idx)*scale
end do

end subroutine fv_prepare_realgas_update_state


subroutine fv_state_density_pressure(q,rho,pressure,state_ok)
!> Extract density and pressure from an FV conservative state.
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar),intent(in)::q
real,intent(out)::rho,pressure
logical,intent(out)::state_ok
integer::k,idxe,idxev
real::vel2
real::rho_mix,gamma_mix,ar_sum,vf,vf_sum,denom,stiff_sum
real::etot,evib,echem,etr,rmix,cv_mix,y_i

rho=zero
pressure=zero
state_ok=.false.

if (nof_variables.lt.1) return

rho=q(1)
if ((rho.ne.rho).or.(abs(rho).gt.fv_update_max_abs)) return
if (rho.le.fv_update_min_abs) return

if (nof_variables.lt.dimensiona+2) then
  state_ok=.true.
  return
end if

if ((multispecies.eq.1).and.(mp_modelc.eq.0).and.(nof_species.gt.0)) then
  if ((dimensiona+1+(2*nof_species)).le.nof_variables) then
    rho_mix=zero
    do k=1,nof_species
      rho_mix=rho_mix+q(dimensiona+2+k)
    end do

    if ((rho_mix.le.zero).or.(rho_mix.ne.rho_mix).or. &
        (abs(rho_mix).gt.fv_update_max_abs)) return

    vel2=zero
    if (nof_variables.ge.2) vel2=vel2+(q(2)/rho_mix)*(q(2)/rho_mix)
    if ((dimensiona.ge.2).and.(nof_variables.ge.3)) vel2=vel2+(q(3)/rho_mix)*(q(3)/rho_mix)
    if ((dimensiona.eq.3).and.(nof_variables.ge.4)) vel2=vel2+(q(4)/rho_mix)*(q(4)/rho_mix)

    vf_sum=zero
    ar_sum=zero
    stiff_sum=zero
    do k=1,nof_species-1
      vf=q(dimensiona+2+nof_species+k)
      vf_sum=vf_sum+vf
      denom=gamma_in(k)-1.0d0
      if (abs(denom).le.tolsmall) return
      ar_sum=ar_sum+(vf/denom)
      stiff_sum=stiff_sum+(vf*(gamma_in(k)/denom)*mp_pinf(k))
    end do

    vf=1.0d0-vf_sum
    denom=gamma_in(nof_species)-1.0d0
    if (abs(denom).le.tolsmall) return
    ar_sum=ar_sum+(vf/denom)
    stiff_sum=stiff_sum+(vf*(gamma_in(nof_species)/denom)*mp_pinf(nof_species))
    if (ar_sum.le.tolsmall) return

    gamma_mix=(1.0d0/ar_sum)+1.0d0
    pressure=((gamma_mix-1.0d0)*(q(dimensiona+2)-oo2*rho_mix*vel2))-((gamma_mix-1.0d0)*stiff_sum)
    if ((pressure.ne.pressure).or.(abs(pressure).gt.fv_update_max_abs)) return
    state_ok=.true.
    return
  end if
end if

vel2=zero
if (nof_variables.ge.2) vel2=vel2+(q(2)/rho)*(q(2)/rho)
if ((dimensiona.ge.2).and.(nof_variables.ge.3)) vel2=vel2+(q(3)/rho)*(q(3)/rho)
if ((dimensiona.eq.3).and.(nof_variables.ge.4)) vel2=vel2+(q(4)/rho)*(q(4)/rho)

if ((realgas.eq.1).and.(nof_species.gt.0)) then
  idxe=dimensiona+2
  idxev=dimensiona+3
  if ((idxev+nof_species).le.nof_variables) then
    etot=q(idxe)/rho
    evib=q(idxev)/rho
    echem=zero
    rmix=zero
    cv_mix=zero

    do k=1,nof_species
      if (abs(rg_molm(k)).le.tolsmall) return
      y_i=q(idxev+k)/rho
      if (rg_hzero(k).gt.zero) echem=echem+(y_i*(rg_hzero(k)/rg_molm(k)))
      rmix=rmix+(y_i/rg_molm(k))
      if (k.le.3) then
        cv_mix=cv_mix+(y_i*2.5d0*(rgs_ru/rg_molm(k)))
      else
        cv_mix=cv_mix+(y_i*1.5d0*(rgs_ru/rg_molm(k)))
      end if
    end do

    rmix=rgs_ru*rmix
    if ((rmix.le.tolsmall).or.(cv_mix.le.tolsmall)) return
    etr=etot-evib-echem-(oo2*vel2)
    pressure=rho*rmix*(etr/cv_mix)
    if ((pressure.ne.pressure).or.(abs(pressure).gt.fv_update_max_abs)) return
    state_ok=.true.
    return
  end if
end if

if (gamma.le.1.0d0) return
pressure=(gamma-1.0d0)*(q(dimensiona+2)-oo2*rho*vel2)
if ((pressure.ne.pressure).or.(abs(pressure).gt.fv_update_max_abs)) return
state_ok=.true.

end subroutine fv_state_density_pressure

logical function fv_state_pair_admissible(qold,qnew)
!> Check dimensional reference bounds and local density/pressure update ratios.
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar),intent(in)::qold,qnew
logical::old_ok,new_ok
integer::k
real::rho_old,rho_new,p_old,p_new
real::rho_scale,p_scale,rho_min,p_min,rho_max,p_max,rho_ratio,p_ratio
real::species_sum,rhoy_new

fv_state_pair_admissible=.false.

call fv_state_density_pressure(qold,rho_old,p_old,old_ok)
call fv_state_density_pressure(qnew,rho_new,p_new,new_ok)

if (.not.new_ok) return

if ((realgas.eq.1).and.(nof_species.gt.0)) then
  if ((dimensiona+3+nof_species).gt.nof_variables) return
  species_sum=zero
  do k=1,nof_species
    rhoy_new=qnew(dimensiona+3+k)
    if ((rhoy_new.ne.rhoy_new).or.(abs(rhoy_new).gt.fv_update_max_abs)) return
    if (rhoy_new.lt.0.0d0) return
    species_sum=species_sum+rhoy_new
  end do
  if (species_sum.le.fv_update_min_abs) return
  if (abs(species_sum-rho_new).gt.fv_update_species_closure_tol*max(abs(rho_new),fv_update_min_abs)) return
end if

rho_scale=max(abs(rres),abs(rho_old),1.0d0)
p_scale=max(abs(pres),abs(p_old),1.0d0)
rho_min=max(fv_update_min_abs,fv_update_floor_fraction*rho_scale)
p_min=max(fv_update_min_abs,fv_update_floor_fraction*p_scale)
rho_max=min(fv_update_max_abs,fv_update_max_ref_ratio*rho_scale)
p_max=min(fv_update_max_abs,fv_update_max_ref_ratio*p_scale)

if (rho_new.le.rho_min) return
if (p_new.le.p_min) return
if (rho_new.gt.rho_max) return
if (p_new.gt.p_max) return

if (old_ok) then
  if ((rho_old.gt.rho_min).and.(p_old.gt.p_min)) then
    rho_ratio=rho_new/rho_old
    p_ratio=p_new/p_old
    if (rho_ratio.lt.fv_update_min_ratio) return
    if (rho_ratio.gt.fv_update_max_ratio) return
    if (p_ratio.lt.fv_update_min_ratio) return
    if (p_ratio.gt.fv_update_max_ratio) return
  end if
end if

fv_state_pair_admissible=.true.

end function fv_state_pair_admissible

subroutine fv_limit_candidate_update(qold,qcandidate,qaccepted,alpha,accepted)
!> For real gas, profile 999, or implicit FV updates, bisect the
!> update until the state is admissible.
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar),intent(in)::qold,qcandidate
real,dimension(1:gpu_max_nvar),intent(out)::qaccepted
real,intent(out)::alpha
logical,intent(out)::accepted
integer::tries,iv
real::trial_alpha,low_alpha,high_alpha,best_alpha
logical::trial_ok,found_alpha,needs_admissibility
real,dimension(1:gpu_max_nvar)::qtrial,qprojected

accepted=.false.
alpha=zero
qaccepted=zero
qtrial=zero
qprojected=zero
qaccepted(1:nof_variables)=qold(1:nof_variables)

needs_admissibility=(realgas.eq.1).or.(code_profile.eq.999).or.&
                    (rungekutta.eq.10).or.(rungekutta.eq.11)

if (.not.needs_admissibility) then
  qaccepted(1:nof_variables)=qcandidate(1:nof_variables)
  alpha=1.0d0
  accepted=.true.
  return
end if

trial_alpha=1.0d0
do iv=1,nof_variables
  qtrial(iv)=qold(iv)+trial_alpha*(qcandidate(iv)-qold(iv))
end do

call fv_prepare_realgas_update_state(qtrial,qprojected,trial_ok)
if (trial_ok) then
  if (fv_state_pair_admissible(qold,qprojected)) then
    qaccepted(1:nof_variables)=qprojected(1:nof_variables)
    alpha=trial_alpha
    accepted=.true.
    return
  end if
end if

low_alpha=0.0d0
high_alpha=1.0d0
best_alpha=0.0d0
found_alpha=.false.

do tries=1,fv_update_max_backtracks
  trial_alpha=0.5d0*(low_alpha+high_alpha)
  do iv=1,nof_variables
    qtrial(iv)=qold(iv)+trial_alpha*(qcandidate(iv)-qold(iv))
  end do

  call fv_prepare_realgas_update_state(qtrial,qprojected,trial_ok)
  if (trial_ok) then
    if (fv_state_pair_admissible(qold,qprojected)) then
      best_alpha=trial_alpha
      qaccepted(1:nof_variables)=qprojected(1:nof_variables)
      found_alpha=.true.
      low_alpha=trial_alpha
    else
      high_alpha=trial_alpha
    end if
  else
    high_alpha=trial_alpha
  end if
end do

if (found_alpha) then
  alpha=best_alpha
  accepted=.true.
  return
end if

end subroutine fv_limit_candidate_update

logical function fv_state_update_admissible(q)
!> Backward-compatible candidate-only FV state check.
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,dimension(1:gpu_max_nvar),intent(in)::q
logical::state_ok
real::rho,pressure

call fv_state_density_pressure(q,rho,pressure,state_ok)
fv_state_update_admissible=state_ok
if (.not.state_ok) return
if (rho.le.max(fv_update_min_abs,fv_update_floor_fraction*max(abs(rres),1.0d0))) fv_state_update_admissible=.false.
if (pressure.le.max(fv_update_min_abs,fv_update_floor_fraction*max(abs(pres),1.0d0))) fv_state_update_admissible=.false.

end function fv_state_update_admissible

subroutine update_cfl_ramp(n)
!> Residual-based CFL ramp for steady/local implicit FV runs.
implicit none
integer,intent(in)::n
integer::nres
real::res_now,cfl_old,growth

!$omp barrier
!$omp master

if (cflramp.eq.1) then

  if (cflmax.le.zero) cflmax=max(cfl,1.0d0)
  cfl_old=cfl

  if (dimensiona.eq.3) then
    nres=min(5,nof_variables)
  else
    nres=min(4,nof_variables)
  end if
  if (turbulence.gt.0) nres=min(15,nof_variables+turbulenceequations)

  if (nres.gt.0) then
    res_now=maxval(allres(1:nres))

    if ((res_now.eq.res_now).and.(abs(res_now).lt.tolbig)) then
      cfl=min(cfl,cflmax)

      if ((prevres.le.zero).or.(prevres.ne.prevres).or.(abs(prevres).gt.tolbig)) then
        prevres=res_now
      else
        if (res_now.le.(0.98d0*prevres)) then
          growth=1.05d0
          if (res_now.le.(0.70d0*prevres)) growth=1.10d0
          if (res_now.le.(0.40d0*prevres)) growth=1.20d0
          cfl=min(cflmax,max(cfl,1.0d-6)*growth)
        else if (res_now.gt.(1.25d0*prevres)) then
          cfl=max(1.0d-3,0.70d0*cfl)
        end if
        prevres=res_now
      end if

      if (dimensiona.eq.3) then
        ccfl=cfl/3.0d0
      else
        ccfl=cfl/2.0d0
      end if

      if ((n.eq.0).and.(abs(cfl-cfl_old).gt.(1.0d-12*max(1.0d0,abs(cfl_old))))) then
        open(63,file='history.txt',form='formatted',status='unknown',action='write',position='append')
        write(63,'(A,I10,A,ES14.7,A,ES14.7,A,ES14.7)') 'cfl_ramp it=', it, ' residual=', res_now, &
                                                           ' cfl=', cfl, ' cflmax=', cflmax
        close(63)
      end if
    end if
  end if

end if

#if defined(gpu) || defined(xpu)
!$omp target update to(cfl,ccfl)
#endif

!$omp end master
!$omp barrier

end subroutine update_cfl_ramp

subroutine arbitrary_order_master(n)
!> AO-only GPU mode launches reconstruction from one host thread.
implicit none
integer,intent(in)::n
#if defined(gpu) && defined(GPU_AO_ONLY)
!$omp barrier
!$omp master
#endif
call arbitrary_order(n)
#if defined(gpu) && defined(GPU_AO_ONLY)
!$omp end master
!$omp barrier
#endif
end subroutine arbitrary_order_master



subroutine dg_mass_rhs_product_cell(i)
implicit none
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::i
integer::mm_idof,mm_jdof,mm_var
real::mm_sum

do mm_var=1,nof_variables
  do mm_idof=1,num_dg_dofs
    mm_sum=zero
    do mm_jdof=1,num_dg_dofs
      mm_sum=mm_sum + m_1_val(mm_idof,mm_jdof,i)*rhs_valdg(mm_jdof,mm_var,i)
    end do
    rhs_sol_mm_dg(mm_idof,mm_var,i)=mm_sum
  end do
end do

end subroutine dg_mass_rhs_product_cell

!---------------------------------------------------------------------------------------------!
end module advance
