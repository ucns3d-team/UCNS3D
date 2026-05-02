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

 contains
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
integer::i,k,l,kmaxe,j,ingtmax,ingtmin,whgu,whgl,srf
real::suvi,suv3,maxu,minu
real::veln,agrt
real,dimension(1:nof_variables)::leftv,rightv
real,dimension(1:nof_variables)::srf_speed
real::mp_pinfl,gammal
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
real,dimension(1)::etvm

kmaxe=xmpielrank(n)
       
#ifdef gpu
!$omp target map(present: dt)
#endif
  dt = tolbig
#ifdef gpu
!$omp end target
#endif
	if (itestcase.lt.3)then
#ifdef gpu
	!$omp target teams distribute parallel do map(present: dt) &
    !$omp& reduction(min: dt) private(veln)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe
		veln=max(abs(lamx),abs(lamy),abs(lamz))
		
		if (dg.eq.1)then
		
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))
		
		else
		
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln))))
		end if
		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	if (itestcase.eq.3)then
#ifdef gpu
	!$omp target teams distribute parallel do map(present: dt) &
    !$omp& reduction(min: dt) private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz,srf_speed,srf)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe
        
        
        
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
        if (srfg.eq.1) then
            pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
            poy(1:3)=srf_velocity
            srf_speed=zero
            srf_speed(2:4)=vect_function(pox,poy)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
            srf=rec_mrf(i)
            if (rec_mrf(i).eq.0)then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox=pox-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                srf_speed(2:4)=vect_function(pox,poy)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
		if (dg.eq.1)then
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))
		
		else
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln))))
		
		end if
		
		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
	
	if (itestcase.eq.4)then
#ifdef gpu
	!$omp target teams distribute parallel do map(present: dt) &
    !$omp& reduction(min: dt) private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz,srf_speed,srf,&
    !$omp&  viscl,laml,turbmv,etvm,eddyfl,eddyfr)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe
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
            srf_speed(2:4)=vect_function(pox,poy)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if          
        if(mrf.eq.1)then
            srf=rec_mrf(i)
            if (rec_mrf(i).eq.0)then
                veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox(1:3)=pox(1:3)-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                srf_speed(2:4)=vect_function(pox,poy)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
		if (turbulence.eq.1)then
		if (turbulencemodel.eq.1)then
		turbmv(1)=u_ct_val(1,1,i);  turbmv(2)=u_ct_val(1,1,i);
		eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(1)+laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if

		if (dg.eq.1)then
		

		dt=min(dt,(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem_minedge(i))))))
		
		
		else

         dt=min(dt,ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*(laml(1)+viscl(1))/((ielem_minedge(i)))**2))))



         end if
                             
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
        return
        
end subroutine calculate_cfl

subroutine calculate_cfll(n)
!> @brief
!> subroutine for computing the time step size for each cell in 3d
implicit none
integer,intent(in)::n
integer::i,k,l,kmaxe,j,ingtmax,ingtmin,whgu,whgl,srf
real::suvi,suv3,maxu,minu
real::veln,agrt
real,dimension(1:nof_variables)::leftv,rightv
real,dimension(1:nof_variables)::srf_speed
real::mp_pinfl,gammal
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
kmaxe=xmpielrank(n)
       

        
        
	if (itestcase.lt.3)then
#ifdef gpu
	!$omp target teams distribute parallel do private (veln)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		veln=max(abs(lamx),abs(lamy),abs(lamz))
		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	if (itestcase.eq.3)then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz,srf_speed,srf)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
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
        if (srfg.eq.1) then
            pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
            poy(1:3)=srf_velocity
            srf_speed=zero
            srf_speed(2:4)=vect_function(pox,poy)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
            srf=rec_mrf(i)
            if (rec_mrf(i).eq.1)then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox=pox-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                srf_speed(2:4)=vect_function(pox,poy)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
		if (dg.eq.1)then
		
		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1))
		
		else
		
		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))
		
		end if
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
	
	if (itestcase.eq.4)then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz,srf_speed,srf,&
    !$omp&  viscl,laml,turbmv,etvm,eddyfl,eddyfr)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe

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
            srf_speed(2:4)=vect_function(pox,poy)
            veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
        end if
        if(mrf.eq.1)then
            srf=rec_mrf(i)
           if (rec_mrf(i).eq.0)then
            veln=max(abs(leftv(2)),abs(leftv(3)),abs(leftv(4)))+agrt
            else
                pox(1)=ielem_xxc(i);pox(2)=ielem_yyc(i);pox(3)=ielem_zzc(i)
                pox(1:3)=pox(1:3)-rec_mrf_origin(1:3,i)
                poy(1:3)=rec_mrf_velocity(1:3,i)
                srf_speed=zero
                srf_speed(2:4)=vect_function(pox,poy)
                veln=max(abs(leftv(2)-srf_speed(2)),abs(leftv(3)-srf_speed(3)),abs(leftv(4)-srf_speed(4)))+agrt
            end if
        end if
        if (turbulence.eq.1)then
		if (turbulencemodel.eq.1)then
		turbmv(1)=u_ct_val(1,1,i);  turbmv(2)=u_ct_val(1,1,i);
		eddyfl(2)=turbmv(1); eddyfr(2)=turbmv(2)
		call eddyvisco(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(1)+laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if
		

		if (dg.eq.1)then


		ielem_dtl(i)=(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem_minedge(i)))))

		
		else
		
		ielem_dtl(i)=ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*(laml(1)+viscl(1))/((ielem_minedge(i)))**2)))
		end if
		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
        return
        
end subroutine calculate_cfll



subroutine calculate_cfl2d(n)
!> @brief
!> subroutine for computing the global time step size in 2d
implicit none
integer,intent(in)::n
integer::i,k,l,kmaxe,j,ingtmax,ingtmin,whgu,whgl
real::suvi,suv3,maxu,minu,sum_dt1,sum_dt2
real::veln,agrt,lamxl,lamyl
real,dimension(1:nof_variables)::leftv,rightv
real,dimension(1:nof_variables)::srf_speed
real::mp_pinfl,gammal
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
kmaxe=xmpielrank(n)
       

        
#ifdef gpu
!$omp target map(present: dt)
#endif
  dt = tolbig
#ifdef gpu
!$omp end target
#endif
        
        

        
	if (itestcase.lt.3)then
#ifdef gpu
	!$omp target teams distribute parallel do map(present: dt) &
    !$omp& reduction(min: dt) private(veln,lamxl,lamyl)
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
		
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))
		
		else
		
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln))))
		end if
		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	if (itestcase.eq.3)then
#ifdef gpu
	!$omp target teams distribute parallel do map(present: dt) &
    !$omp& reduction(min: dt) private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe
        
        
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		
        
		call cons2prim(n,leftv,mp_pinfl,gammal)
		if ((multispecies.eq.1).or.(realgas.eq.1))then
		agrt=sqrt((leftv(4)+mp_pinfl)*gammal/leftv(1))
		else
		agrt=sqrt(leftv(4)*gamma/leftv(1))
		end if
		veln=max(abs(leftv(2)),abs(leftv(3)))+agrt
		
		
		
		if (dg.eq.1)then
		
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln)))*(1.0d0/(2*iorder+1)))
		
		
		else
		
		dt=min(dt,ccfl*((ielem_minedge(i))/(abs(veln))))
		end if
		
		
		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
	
	if (itestcase.eq.4)then
#ifdef gpu
	!$omp target teams distribute parallel do map(present: dt) &
    !$omp& reduction(min: dt) private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz,&
    !$omp&  viscl,laml,turbmv,etvm,eddyfl,eddyfr)
#else
    !$omp barrier
	!$omp do reduction (min:dt)
#endif
        do i=1,kmaxe

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
		call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(1)+laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if


    if (dg.eq.1)then
		

      dt=min(dt,(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem_minedge(i))))))
      
      
      else


		
		dt=min(dt,ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*(laml(1)+viscl(1))/((ielem_minedge(i)))**2))))
      end if
  
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
        return
        
end subroutine calculate_cfl2d


subroutine calculate_cfll2d(n)
!> @brief
!> subroutine for computing the time step size for each cell in 2d
implicit none
integer,intent(in)::n
integer::i,k,l,kmaxe,j,ingtmax,ingtmin,whgu,whgl
real::suvi,suv3,maxu,minu
real::veln,agrt
real,dimension(1:nof_variables)::leftv,rightv
real,dimension(1:nof_variables)::srf_speed
real::mp_pinfl,gammal,mp_pinfr,gammar
real,dimension(1:dimensiona)::pox,poy,poz
real,dimension(1:4)::viscl,laml
real,dimension(1:20)::eddyfl,eddyfr
real,dimension(1:2)::turbmv
real,dimension(1)::etvm
kmaxe=xmpielrank(n)
       

        
        
	if (itestcase.lt.3)then
#ifdef gpu
	!$omp target teams distribute parallel do private (veln)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		veln=max(abs(lamx),abs(lamy))
		ielem_dtl(i)=ccfl*((ielem_minedge(i))/(abs(veln)))
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	if (itestcase.eq.3)then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
		leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
		call cons2prim(n,leftv,mp_pinfl,gammal)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
	
	if (itestcase.eq.4)then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(veln,leftv,mp_pinfl,gammal,agrt,veln,pox,poy,poz,&
    !$omp&  viscl,laml,turbmv,etvm,eddyfl,eddyfr)
#else
    !$omp barrier
	!$omp do
#endif
        do i=1,kmaxe
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
		call eddyvisco2d(n,viscl,laml,turbmv,etvm,eddyfl,eddyfr,leftv,rightv)
		laml(1)=laml(1)+laml(3)
		viscl(1)=viscl(1)+viscl(3)
		end if
		end if






		if (dg.eq.1)then


      ielem_dtl(i)=(ccfl/(2*iorder+1))*(ielem_minedge(i)/((abs(veln))+(2.0d0*max(((4.0/3.0)*viscl(1)/leftv(1)),gamma*laml(1)/(prandtl*leftv(1)))*((2*iorder+1)/ielem_minedge(i)))))


      else
		
		ielem_dtl(i)=ccfl*(1.0d0/((abs(veln)/((ielem_minedge(i)))) + (0.5d0*(laml(1)+viscl(1))/((ielem_minedge(i)))**2)))


		end if
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end if
	
	
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
    
    end select
end if


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if
 
 if (mood.eq.1)then
 
 call mood_operator_2(n)
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
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
    end select
end if

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

 if (mood.eq.1)then
 
 call mood_operator_2(n)
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
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
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if (mood.eq.1)then
 
 call mood_operator_2(n)
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

call mood_operator_1(n)

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_3d


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
        
         u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(2,1:nof_variables,:,i) - dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))



         
    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif





if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if
 
call call_flux_subroutines_3d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
         u_c_valdg(1,1:nof_variables,:,i)=to4*u_c_valdg(2,1:nof_variables,:,i) + oo4*u_c_valdg(3,1:nof_variables,:,i) - oo4*dt* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))



    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
        ((rhs_val(1:nof_variables,i))*(oovolume))))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if






call call_flux_subroutines_3d


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  if (dg == 1) then



        rhs_sol_mm_dg(1:num_dg_dofs,1:nof_variables,i)=matmul(m_1_val(1:num_dg_dofs,1:num_dg_dofs,i),rhs_valdg(1:num_dg_dofs,1:nof_variables,i))

        
        u_c_valdg(1,1:nof_variables,:,i)=oo3*u_c_valdg(2,1:nof_variables,:,i) + to3*u_c_valdg(1,1:nof_variables,:,i) - to3*dt * transpose(rhs_sol_mm_dg(:,:,i))


    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(1,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
        ((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe,kx
real::avrgs,oovolume
real::t1,t2,t3
kmaxe=xmpielrank(n)




call call_flux_subroutines_3d



#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  
 if (dg == 1) then
        
        
        u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i) - dt* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))!*oovolume
        
        
  else
  
  
  
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  end if
  
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#ifdef gpu
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


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



call call_flux_subroutines_3d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
 oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-(dt*oo2*(rhs_val(1:nof_variables,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(oo2*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo2*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(dt*oo2*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_3d







#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  
  
  if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
        
        u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(2,1:nof_variables,:,i) - ielem_dtl(i) * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
         
    else
    
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(ielem_dtl(i)*(rhs_val(1:nof_variables,i)*oovolume))
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(ielem_dtl(i)*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


call call_flux_subroutines_3d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  
  if (dg == 1) then
        
        
        u_c_valdg(1,1:nof_variables,:,i)=(oo2*u_c_valdg(2,1:nof_variables,:,i)) +(oo2*u_c_valdg(1,1:nof_variables,:,i))- (ielem_dtl(i) *oo2* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i))))
         
    else
  
  
  
  u_c_val(1,1:nof_variables,i)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-(ielem_dtl(i)*oo2*(rhs_val(1:nof_variables,i)*oovolume))
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(oo2*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo2*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(ielem_dtl(i)*oo2*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_2d


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe

  if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
        u_c_valdg(2,:,:,i)=u_c_valdg(1,:,:,i)
        u_c_valdg(1,:,:,i)=u_c_valdg(2,:,:,i) - (ielem_dtl(i)* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,:,i))))!*oovolume
    else
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(ielem_dtl(i)*(rhs_val(1:nof_variables,i)*oovolume))
  end if



  
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(ielem_dtl(i)*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if



call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe

   if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
        u_c_valdg(1,:,:,i)=(oo2*u_c_valdg(2,:,:,i))+(oo2*u_c_valdg(1,:,:,i))-(ielem_dtl(i)* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,:,i))))
    else
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-(ielem_dtl(i)*oo2*(rhs_val(1:nof_variables,i)*oovolume))

  end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(oo2*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo2*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(ielem_dtl(i)*oo2*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_2d


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
        u_c_valdg(2,:,:,i)=u_c_valdg(1,:,:,i)
        u_c_valdg(1,:,:,i)=u_c_valdg(2,:,:,i) - (dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,:,i))))!*oovolume
    else
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
if (dg == 1) then
oovolume=1.0d0/ielem_totvolume(i)
        u_c_valdg(1,:,:,i)=oo2*u_c_valdg(2,:,:,i) + oo2*u_c_valdg(1,:,:,i) - (oo2*dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,:,i))))!*oovolume
    else
 oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(1,1:nof_variables,i)=(oo2*u_c_val(2,1:nof_variables,i))+(oo2*u_c_val(1,1:nof_variables,i))-(dt*oo2*(rhs_val(1:nof_variables,i)*oovolume))
  end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(u_ct_val(2,1:turbulenceequations+passivescalar,i))-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
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
real,dimension(1:nof_variables)::solution_integ2

  
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
real,dimension(1:nof_variables)::solution_integ2
real,dimension(1:nof_variables)::solution_integ_weak
real,dimension(1:nof_variables)::solution_integ_strong


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
real,dimension(1:nof_variables)::solution_integ2
 
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
integer,intent(in)::n
integer::i,kmaxe,kx
real::avrgs,oovolume
real,dimension(nof_variables)::tempsol
kmaxe=xmpielrank(n)



call call_flux_subroutines_2d





#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  
  
  if (dg == 1) then
        
        u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i) - dt* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
        
  else
  
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
  end if
  

end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)


call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
        

        
          u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(2,1:nof_variables,:,i) - dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
         
    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
        
        
          u_c_valdg(1,1:nof_variables,:,i)=to4*u_c_valdg(2,1:nof_variables,:,i) + oo4*u_c_valdg(3,1:nof_variables,:,i) - oo4*dt* transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
        ((rhs_val(1:nof_variables,i))*(oovolume))))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    if (dg == 1) then
    

        u_c_valdg(1,1:nof_variables,:,i)=oo3*u_c_valdg(2,1:nof_variables,:,i) + to3*u_c_valdg(1,1:nof_variables,:,i) - to3*dt*transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
    else
        oovolume=1.0d0/ielem_totvolume(i)
        u_c_val(1,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(1,1:nof_variables,i))-(((to3))*&
        ((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
  end do
#ifdef gpu
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


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*(rhs_val(1:nof_variables,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if
 
 if (mood.eq.1)then
 
 call mood_operator_2(n)
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

call mood_operator_1(n)

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
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

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=(to4*u_c_val(2,1:nof_variables,i))+(oo4*u_c_val(3,1:nof_variables,i))-(((oo4))*((dt)*&
((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif






if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(to4*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(oo4*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((oo4))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

 if (mood.eq.1)then
 
 call mood_operator_2(n)
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
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
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
  u_c_val(4,1:nof_variables,i)=((oo3)*u_c_val(2,1:nof_variables,i))+((to3)*u_c_val(3,1:nof_variables,i))-(((to3))*&
((dt)*((rhs_val(1:nof_variables,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if (mood.eq.1)then
 
 call mood_operator_2(n)
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!
call mood_operator_1(n)

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
u_ct_val(1,1:turbulenceequations+passivescalar,i)=((oo3)*u_ct_val(2,1:turbulenceequations+passivescalar,i))+((to3)*u_ct_val(1,1:turbulenceequations+passivescalar,i))-(((to3))*&
((dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
real::dumpracein,dumpraceout,flops_count
kmaxe=xmpielrank(n)




call call_flux_subroutines_3d




    
    




#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(2,1:nof_variables,:,i) - dt * 0.391752226571890 * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))!*oovolume
    else
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*0.391752226571890*(rhs_val(1:nof_variables,i)*oovolume))
    end if
  
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*0.391752226571890*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
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


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(1,1:nof_variables,:,i) = 0.444370493651235 * u_c_valdg(2,1:nof_variables,:,i) + 0.555629506348765 * u_c_valdg(3,1:nof_variables,:,i) - 0.368410593050371 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
    else
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=(0.444370493651235 * u_c_val(2,1:nof_variables,i)) + (0.555629506348765 * u_c_val(3,1:nof_variables,i)) - 0.368410593050371 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
    u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.444370493651235*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.555629506348765*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((0.368410593050371))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


            call call_flux_subroutines_3d


#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
   if (dg == 1) then
        u_c_valdg(4,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(1,1:nof_variables,:,i) = 0.620101851488403 * u_c_valdg(2,1:nof_variables,:,i) + 0.379898148511597 * u_c_valdg(4,1:nof_variables,:,i) - 0.251891774271694 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
    else
        u_c_val(4,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i) = 0.620101851488403 * u_c_val(2,1:nof_variables,i) + 0.379898148511597 * u_c_val(4,1:nof_variables,i) - 0.251891774271694 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(4,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.620101851488403*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.379898148511597*u_ct_val(4,1:turbulenceequations+passivescalar,i))-(((0.251891774271694))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

            call call_flux_subroutines_3d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(5,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(6,1:nof_variables,:,i) = - dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
        u_c_valdg(1,1:nof_variables,:,i) = 0.178079954393132 * u_c_valdg(2,1:nof_variables,:,i) + 0.821920045606868 * u_c_valdg(5,1:nof_variables,:,i) - 0.544974750228521 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
    else
        u_c_val(5,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(6,1:nof_variables,i) = - dt * rhs_val(1:nof_variables,i) * oovolume
        u_c_val(1,1:nof_variables,i) = 0.178079954393132 * u_c_val(2,1:nof_variables,i) + 0.821920045606868 * u_c_val(5,1:nof_variables,i) - 0.544974750228521 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
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
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if


    call call_flux_subroutines_3d

#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  if (dg == 1) then
        u_c_valdg(1,1:nof_variables,:,i) = (0.00683325884039 * u_c_valdg(2,1:nof_variables,:,i)) + (0.517231671970585 * u_c_valdg(4,1:nof_variables,:,i)) + (0.12759831133288 * u_c_valdg(5,1:nof_variables,:,i)) + (0.34833675773694 * u_c_valdg(1,1:nof_variables,:,i)) + (0.08460416338212 * u_c_valdg(6,1:nof_variables,:,i)) - 0.22600748319395 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
    else
        u_c_val(1,1:nof_variables,i) = (0.00683325884039 * u_c_val(2,1:nof_variables,i)) + (0.517231671970585 * u_c_val(4,1:nof_variables,i)) +  (0.12759831133288 * u_c_val(5,1:nof_variables,i)) + (0.34833675773694 * u_c_val(1,1:nof_variables,i)) + (0.08460416338212 * u_c_val(6,1:nof_variables,i)) - (0.22600748319395 * dt * rhs_val(1:nof_variables,i) * oovolume)
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do &
    !$omp& private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
   u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.00683325884039*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.517231671970585*u_ct_val(4,1:turbulenceequations+passivescalar,i))+&
				(0.12759831133288*u_ct_val(5,1:turbulenceequations+passivescalar,i))+(0.34833675773694*u_ct_val(1,1:turbulenceequations+passivescalar,i))+&
				(0.08460416338212*u_ct_val(6,1:turbulenceequations+passivescalar,i))-(0.22600748319395*(dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume)))
end do
#ifdef gpu
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
integer::kmaxe,i
kmaxe=xmpielrank(n)


    if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t1=mpi_wtime()
     !$omp end master
     !$omp barrier

    end if


    if (realgas.eq.1)call normalise_species(n)

    if (dg.eq.1)then

    call sol_integ_dg(n) ! calculates cell average of dg solution for fv

    end if

      if ((dg.eq.1).and.(filtering.eq.1))then
#ifdef gpu
	!$omp target teams distribute parallel do
#else
	!$omp do
#endif
        do i=1,kmaxe
        ielem_filtered(i)=0
        end do
#ifdef gpu
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

        call trouble_indicator1 ! checks for troubled cells
        
    end if
    
    call arbitrary_order(n)
    
    
    if (dg == 1) then
        
        call trouble_indicator2 ! changes dg to fv

    end if

    call exhboundhigher(n)



    if (dg.eq.1)then

        call exhboundhigher_dg(n)

        if (itestcase.eq.4)then

          if( br2_yn.eq.2) then

          call reconstruct_br2_dg

          call exhboundhigher_dg2(n)

          end if

          if( br2_yn.eq.0) then
          call viscous_dg_ggs(n)
          end if



        end if

    end if
    
     if (statistics.eq.1)then
    !$omp barrier
    !$omp master
    pr_t4=mpi_wtime()
    prace_t3=pr_t4-pr_t3
    !$omp end master
    !$omp barrier
    end if






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
  real :: rho, sumrhoy, scale,mismatch,rel_err,rhoy_old_n2
  real, parameter :: epsrho = 1.0d-14
  real, parameter :: epssum = 1.0d-300
  real,parameter::tol_mass_closure = 1.0d-12
  real, dimension(1:nof_species) :: rhoy
  real,dimension(1:nof_variables)::leftv
  logical :: bad

  kmaxe = xmpielrank(n)

#ifdef gpu
	!$omp target teams distribute parallel do private(rho,sumrhoy,scale,mismatch,rel_err,rhoy_old_n2,rhoy,leftv,k)
#else
	!$omp do
#endif
  do i = 1, kmaxe


    leftv(1:nof_variables)=u_c_val(1,1:nof_variables,i)
    call fix_conservative_state(leftv)
    u_c_val(1,1:nof_variables,i)=leftv(1:nof_variables)

    rho = u_c_val(1,1,i)


    if (rho <= epsrho) then
      rho = max(rho, epsrho)
      do k = 1, nof_species
        u_c_val(1,dimensiona+3+k,i) = rho * rg_vf(k)
      end do
      cycle
    end if

    ! 1) read, nan/inf->0, clip negatives, accumulate
    sumrhoy = 0.0d0
    do k = 1, nof_species
      rhoy(k) = u_c_val(1,dimensiona+3+k,i)

      ! nan check: (x /= x) is true only for nan
      bad = (rhoy(k) /= rhoy(k))
      ! inf/huge check (portable enough): treat absurdly large as bad
      if (.not. bad) bad = (abs(rhoy(k)) > huge(rhoy(k))*0.5d0)

      if (bad) rhoy(k) = 0.0d0
      if (rhoy(k) < 0.0d0) rhoy(k) = 0.0d0

      sumrhoy = sumrhoy + rhoy(k)


    end do

    ! 2) enforce sum(rhoy)=rho without changing rho
    if (sumrhoy > epssum) then


      scale = rho / sumrhoy
      do k = 1, nof_species
        rhoy(k) = rhoy(k) * scale
      end do
      do k = 1, nof_species
      u_c_val(1,dimensiona+3+k,i) = rhoy(k)
    end do

    else
      ! everything got wiped out -> reset to reference mixture

      do k = 1, nof_species
        rhoy(k) = rho * rg_vf(k)
      end do
      do k = 1, nof_species
      u_c_val(1,dimensiona+3+k,i) = rhoy(k)
    end do

    end if

    ! 3) write back


  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
!








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
        call trouble_indicator1
        
    end if
    
    
        
        call arbitrary_order(n)
    
    
    if (dg == 1) then
        
        call trouble_indicator2

    end if
    
    if (bound_lim == 1) then
            call vfbp_limiter
          end if
    
    
    call exhboundhigher(n)
    
    if (dg.eq.1)then
    call exhboundhigher_dg(n)

    if (itestcase.eq.4)then

          if( br2_yn.eq.2) then

          call reconstruct_br2_dg

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
integer,intent(in)::n
integer::i,kmaxe
real::avrgs,oovolume
kmaxe=xmpielrank(n)



call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)

    if (dg == 1) then
        u_c_valdg(2,1:nof_variables,:,i)=u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(1,1:nof_variables,:,i)=u_c_valdg(2,1:nof_variables,:,i) - dt * 0.391752226571890 * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))!*oovolume
    else
        u_c_val(2,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=u_c_val(2,1:nof_variables,i)-(dt*0.391752226571890*(rhs_val(1:nof_variables,i)*oovolume))
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(2,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=u_ct_val(2,1:turbulenceequations+passivescalar,i)-(dt*0.391752226571890*(rhst_val(1:turbulenceequations+passivescalar,i)*oovolume))
  end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  
    if (dg == 1) then
        u_c_valdg(3,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(1,1:nof_variables,:,i) = 0.444370493651235 * u_c_valdg(2,1:nof_variables,:,i) + 0.555629506348765 * u_c_valdg(3,1:nof_variables,:,i) - 0.368410593050371 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
    else
        u_c_val(3,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i)=(0.444370493651235 * u_c_val(2,1:nof_variables,i)) + (0.555629506348765 * u_c_val(3,1:nof_variables,i)) - 0.368410593050371 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
    u_ct_val(3,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.444370493651235*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.555629506348765*u_ct_val(3,1:turbulenceequations+passivescalar,i))-(((0.368410593050371))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if

call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)
    
    if (dg == 1) then
        u_c_valdg(4,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(1,1:nof_variables,:,i) = 0.620101851488403 * u_c_valdg(2,1:nof_variables,:,i) + 0.379898148511597 * u_c_valdg(4,1:nof_variables,:,i) - 0.251891774271694 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
    else
        u_c_val(4,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(1,1:nof_variables,i) = 0.620101851488403 * u_c_val(2,1:nof_variables,i) + 0.379898148511597 * u_c_val(4,1:nof_variables,i) - 0.251891774271694 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
  u_ct_val(4,1:turbulenceequations+passivescalar,i)=u_ct_val(1,1:turbulenceequations+passivescalar,i)
  u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.620101851488403*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.379898148511597*u_ct_val(4,1:turbulenceequations+passivescalar,i))-(((0.251891774271694))*((dt)*&
((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume))))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
 end if
 
call call_flux_subroutines_2d

#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)
    
    if (dg == 1) then
        u_c_valdg(5,1:nof_variables,:,i) = u_c_valdg(1,1:nof_variables,:,i)
        u_c_valdg(6,1:nof_variables,:,i) = - dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
        u_c_valdg(1,1:nof_variables,:,i) = 0.178079954393132 * u_c_valdg(2,1:nof_variables,:,i) + 0.821920045606868 * u_c_valdg(5,1:nof_variables,:,i) - 0.544974750228521 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))
    else
        u_c_val(5,1:nof_variables,i) = u_c_val(1,1:nof_variables,i)
        u_c_val(6,1:nof_variables,i) = - dt * rhs_val(1:nof_variables,i) * oovolume
        u_c_val(1,1:nof_variables,i) = 0.178079954393132 * u_c_val(2,1:nof_variables,i) + 0.821920045606868 * u_c_val(5,1:nof_variables,i) - 0.544974750228521 * dt * rhs_val(1:nof_variables,i) * oovolume
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
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
#ifdef gpu
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

#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
    oovolume=1.0d0/ielem_totvolume(i)
    
    if (dg == 1) then
        u_c_valdg(1,1:nof_variables,:,i) = (0.00683325884039 * u_c_valdg(2,1:nof_variables,:,i)) + (0.517231671970585 * u_c_valdg(4,1:nof_variables,:,i)) + (0.12759831133288 * u_c_valdg(5,1:nof_variables,:,i)) + (0.34833675773694 * u_c_valdg(1,1:nof_variables,:,i)) + (0.08460416338212 * u_c_valdg(6,1:nof_variables,:,i)) - 0.22600748319395 * dt * transpose(matmul(m_1_val(:,:,i), rhs_valdg(:,1:nof_variables,i)))! * oovolume
    else
        u_c_val(1,1:nof_variables,i) = (0.00683325884039 * u_c_val(2,1:nof_variables,i)) + (0.517231671970585 * u_c_val(4,1:nof_variables,i)) +  (0.12759831133288 * u_c_val(5,1:nof_variables,i)) + (0.34833675773694 * u_c_val(1,1:nof_variables,i)) + (0.08460416338212 * u_c_val(6,1:nof_variables,i)) - (0.22600748319395 * dt * rhs_val(1:nof_variables,i) * oovolume)
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif



if ((turbulence.gt.0).or.(passivescalar.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do private(oovolume)
#else
	!$omp do
#endif
do i=1,kmaxe
  oovolume=1.0d0/ielem_totvolume(i)
   u_ct_val(1,1:turbulenceequations+passivescalar,i)=(0.00683325884039*u_ct_val(2,1:turbulenceequations+passivescalar,i))+(0.517231671970585*u_ct_val(4,1:turbulenceequations+passivescalar,i))+&
				(0.12759831133288*u_ct_val(5,1:turbulenceequations+passivescalar,i))+(0.34833675773694*u_ct_val(1,1:turbulenceequations+passivescalar,i))+&
				(0.08460416338212*u_ct_val(6,1:turbulenceequations+passivescalar,i))-(0.22600748319395*(dt)*((rhst_val(1:turbulenceequations+passivescalar,i))*(oovolume)))
end do
#ifdef gpu
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
integer::i,k,kmaxe,kill_nan_global
integer,intent(in)::n
real::verysmall,du,uold, umin, dumax, unew
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
    if (turbulence.eq.1)then
    call sources_computation(n)
    end if
    end select
end if

if (relax.eq.3)then
 
 call relaxation_lumfree(n)
 
 else
  
if (lowmemory.eq.0)then
   
 call relaxation(n)
   
 else
 
 call relaxation_lm(n)
 
 end if
   
  end if

 kill_nan=0
 
#ifdef gpu
	!$omp target teams distribute parallel do map(tofrom:kill_nan) reduction(max:kill_nan)
#else
	!$omp do reduction(max:kill_nan)
#endif
do i=1,kmaxe
    if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)).or.(impdu(i,5).ne.impdu(i,5)))then
        kill_nan=1
    end if
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
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
call sources_realgas_pi(n,i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if




if ((passivescalar.gt.0).or.(turbulence.gt.0))then
#ifdef gpu
	!$omp target teams distribute parallel do private (du,uold,umin,dumax,unew)
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
real::verysmall,du,uold, umin, dumax, unew
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
 
 call relaxation_lumfree(n)
 
 else
 if (lowmemory.eq.0)then
   
 call relaxation2d(n)
   
 else
 
 call relaxation_lm2d(n)
 
 end if
  end if

 kill_nan=0
#ifdef gpu
!$omp target teams distribute parallel do map(tofrom:kill_nan) reduction(max:kill_nan)
#else
!$omp do reduction(max:kill_nan)
#endif
do i=1,kmaxe
    do j=1,nof_variables
    if ((impdu(i,j).ne.impdu(i,j)))then

        kill_nan=1
    end if
    end do
  u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
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
call sources_realgas_pi(n,i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if

if ((passivescalar.gt.0).or.(turbulence.gt.0))then
#ifdef gpu
!$omp target teams distribute parallel do private (du,uold,umin,dumax,unew)
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
integer::i,k,kmaxe,jj,kill_nan_global
integer,intent(in)::n
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
   if (jj.eq.1)then
    iscoun=1
      else
      iscoun=2
      end if    
    
      
!$omp target update to(iscoun)
      
call call_flux_subroutines_3d


if (relax.eq.3)then



call relaxation_lumfree(n)


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
!$omp target teams distribute parallel do          &
!$omp& private(rsumfacei)                          &
!$omp& reduction(+:allresdt)                       &
!$omp& reduction(max:kill_nan)                     &
!$omp& map(tofrom: allresdt, kill_nan)             &
#else
!$omp do private(rsumfacei) reduction(+:allresdt) reduction(max:kill_nan)
#endif
do i = 1, kmaxe
    rsumfacei = sqrt( impdu(i,1)**2 + impdu(i,2)**2 + impdu(i,3)**2 + &
                      impdu(i,4)**2 + impdu(i,5)**2 )

    allresdt = allresdt + rsumfacei * ielem_totvolume(i)

    if ( (impdu(i,1) /= impdu(i,1)) .or. &
         (impdu(i,2) /= impdu(i,2)) .or. &
         (impdu(i,3) /= impdu(i,3)) .or. &
         (impdu(i,4) /= impdu(i,4)) .or. &
         (impdu(i,5) /= impdu(i,5)) ) then
        kill_nan = 1
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
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
!$omp target teams distribute parallel do
#else
!$omp do
#endif
  do i=1,kmaxe
 u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do k=1,turbulenceequations+passivescalar
	if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
	end if
    end do
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
call sources_realgas_pi(n,i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if










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
	if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
	end if
    end do
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
call sources_realgas_pi(n,i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if

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
    call vortexcalc(n)
    if (turbulence.eq.1)then
    call sources_computation(n)
    end if
    end select
end if


call relaxation_ex(n)





#ifdef gpu
!$omp target teams distribute parallel do          &
!$omp& private(rsumfacei)                          &
!$omp& reduction(+:allresdt)                       &
!$omp& map(tofrom: allresdt)             &
#else
!$omp do private(rsumfacei) reduction(+:allresdt)
#endif
do i=1,kmaxe
      rsumfacei=sqrt(((impdu(i,1))**2)+((impdu(i,2))**2)+((impdu(i,3))**2)+((impdu(i,4))**2)+((impdu(i,5))**2))
      allresdt=allresdt+(rsumfacei*ielem_totvolume(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
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
!$omp target teams distribute parallel do          &
!$omp& private(rsumfacei)                          &
!$omp& reduction(+:allresdt)                       &
!$omp& map(tofrom: allresdt)             &
#else
!$omp do private(rsumfacei) reduction(+:allresdt)
#endif
do i=1,kmaxe
      rsumfacei=sqrt(((impdu(i,1))**2)+((impdu(i,2))**2)+((impdu(i,3))**2)+((impdu(i,4))**2))
      allresdt=allresdt+(rsumfacei*ielem_totvolume(i))
end do
#ifdef gpu
!$omp end target teams distribute parallel do
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
integer::i,k,kmaxe,nvar,jj,kill_nan_global
integer,intent(in)::n
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
      if (jj.eq.1)then
    iscoun=1
      else
      iscoun=2
      end if
      
!$omp target update to(iscoun)
      
call call_flux_subroutines_2d


if (relax.eq.3)then



call relaxation_lumfree(n)


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
!$omp target teams distribute parallel do          &
!$omp& private(rsumfacei)                          &
!$omp& reduction(+:allresdt)                       &
!$omp& reduction(max:kill_nan)                     &
!$omp& map(tofrom: allresdt, kill_nan)             &
#else
!$omp do private(rsumfacei) reduction(+:allresdt) reduction(max:kill_nan)
#endif
do i=1,kmaxe
      rsumfacei=sqrt(((impdu(i,1))**2)+((impdu(i,2))**2)+((impdu(i,3))**2)+((impdu(i,4))**2))
      allresdt=allresdt+(rsumfacei*ielem_totvolume(i))


      if ((impdu(i,1).ne.impdu(i,1)).or.(impdu(i,2).ne.impdu(i,2)).or.(impdu(i,3).ne.impdu(i,3)).or.(impdu(i,4).ne.impdu(i,4)))then
      kill_nan=1
      end if

end do
#ifdef gpu
!$omp end target teams distribute parallel do
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
!$omp target teams distribute parallel do
#else
!$omp do
#endif
  do i=1,kmaxe
 u_c_val(1,1:nof_variables,i)=u_c_val(1,1:nof_variables,i)+impdu(i,1:nof_variables)
    if ((turbulence.gt.0).or.(passivescalar.gt.0))then
    do k=1,turbulenceequations+passivescalar
	if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
	end if
    end do
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
call sources_realgas_pi(n,i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if


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
	if (u_ct_val(1,k,i)+impdu(i,nof_variables+k).ge.zero)then
  u_ct_val(1,k,i)=u_ct_val(1,k,i)+impdu(i,nof_variables+k)
	end if
    end do
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
call sources_realgas_pi(n,i)
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
end if
end if

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
!$omp target teams distribute parallel do private(dt1,j)
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
!$omp target teams distribute parallel do private(dtau,j)
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
!$omp target teams distribute parallel do
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
!$omp  target teams distribute parallel do private(nvar)
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
integer::i,kmaxe,ttime
real::dtiv
real::cput1,cput2,cput3,cput4,cput5,cput6,cput8,timec3,timec1,timec4,timec8,totv1,totv2,dumetg1,dumetg2,tzx1,tzx2,resolx,totens1,totens2,totensx1,totensx2
      kill=0
      t=res_time
      resolx=0.01
      iscoun=1
      kmaxe=xmpielrank(n)
      every_time=((idnint(t/output_freq)) * output_freq)+output_freq

      totv1=0.0



#ifdef gpu
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
        call enstrophy_calc(n)


      end if
      
     do 
            
		    call calculate_cfl(n)
#ifdef gpu
		    !$omp target update from(dt)
#endif
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
                   if ( mod(it, 10) .eq. 0) then
				   if (n.eq.0)then
				  open(63,file='history.txt',form='formatted',status='old',action='write',position='append')
				  write(63,'(A,I10,A,ES20.10,A,ES20.10)') 'it=', it, ' dt=', dt, ' t=', t
				  close(63)
				  end if
				  end if
			!$omp end master 
			!$omp barrier
					
                              if (initcond.eq.95)then
                          totk=0;totens=0;totensx=0.0d0


#ifdef GPU
                          !$omp target teams distribute parallel do          &
                          !$omp& map(tofrom: totk, totens, totensx)          &
                          !$omp& reduction(+:totk,totens,totensx)
#else
                          !$omp do reduction(+:totk,totens,totensx)
#endif

                          do i=1,xmpielrank(n)

                              totk=totk+ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
          (((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2)+((u_c_val(1,4,i)/u_c_val(1,1,i))**2))

                              if (boundtype.eq.1)then

                              totens=totens+(ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
                              ielem_vortex(2,i))
                              else

                              totens=totens+(ielem_vortex(2,i))
                              totensx=totensx+(ielem_vortex(3,i))
                              end if

                          end do
#ifdef gpu
                          !$omp end target teams distribute parallel do
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
                          totv1=totk/((2.0*pi)**3)
                          totens1=totens/(((2.0*pi)**3))
                          totensx1=4.0*totensx/(3.0*((2.0*pi)**3))
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
            		
#ifdef gpu
			!$omp target update to(dt)
#endif
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
			 !$omp end master 
			!$omp barrier

#ifdef gpu
			!$omp target update to(tz1,t)
#endif


			!$omp barrier
			!$omp master
				if (dg.eq.1)then
                      if (code_profile.ne.102)then
                          if ( mod(it, 100) .eq. 0) then
#ifdef gpu
                      !$omp target update from(ielem_condition,ielem_mood_o)
#endif
                            call troubled_history
                          end if
                      end if
                      if ( filtering .eq. 1) then
                        if ( mod(it, 100) .eq. 0) then
#ifdef gpu
                    !$omp target update from(ielem_filtered,ielem_er,ielem_er1,ielem_er2)
#endif
                          call filtered_history
                        end if
                      end if
                end if

          if ( mod(it, 100) .eq. 0) then
#ifdef gpu
          !$omp target update from(ielem_reduce)
#endif
            call reduced_history
          end if

            !$omp end master
			!$omp barrier


			
			
			if (initcond.eq.95)then                    
 				totk=0; totens=0.0; totensx=0.0d0


#ifdef GPU
                          !$omp target teams distribute parallel do          &
                          !$omp& map(tofrom: totk, totens, totensx)          &
                          !$omp& reduction(+:totk,totens,totensx)
#else
                          !$omp  do reduction(+:totk,totens,totensx)
#endif

 				do i=1,xmpielrank(n)
 				
                   
                    totk=totk+ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
(((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2)+((u_c_val(1,4,i)/u_c_val(1,1,i))**2))
                    





                              if (boundtype.eq.1)then

                              totens=totens+(ielem_totvolume(i)*u_c_val(1,1,i)*(1.0/2.0)*&
                              ielem_vortex(2,i))
                              else

                              totens=totens+(ielem_vortex(2,i))
                               totensx=totensx+(ielem_vortex(3,i))
                              end if





				end do
#ifdef gpu
                          !$omp end target teams distribute parallel do
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
                          totens2=totens/(reynolds*((2.0*pi)**3))
                          totensx2=4.0*totensx/(3.0*reynolds*((2.0*pi)**3))
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
				end if

 				

			call mpi_barrier(mpi_comm_world,ierror)
			!$omp end master
			!$omp barrier
 				

          if (adda.eq.1) then

          totk = 0.0d0

#ifdef GPU
        !$omp target teams distribute parallel do     &
        !$omp& map(to: ielem_er)               &
        !$omp& map(tofrom: totk)                &
        !$omp& reduction(+:totk)
          do i = 1, xmpielrank(n)
              totk = totk + ielem_er(i)
          end do
        !$omp end target teams distribute parallel do
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
#ifdef gpu
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
			end if
			end if
			
			
			
		
			
			
			!$omp master
			if (nprobes.gt.0) then
			if ( mod(it, 100) .eq. 0) then
#ifdef gpu
                  !$omp target update from(u_c_val)
#endif
			call probing

			end if
			end if
					
			
			if (timec1.ge.ievery)then
#ifdef gpu
                  !$omp target update from(u_c_val)
                   if (allocated(u_ct_val)) then
                  !$omp target update from(u_ct_val)
                  end if
                  if (allocated(ielem_vortex)) then
                  !$omp target update from(ielem_vortex)
                  end if
                  if (allocated(ielem_troubled)) then
                  !$omp target update from(ielem_troubled)
                  end if
                   if (allocated(ielem_mood_o)) then
                  !$omp target update from(ielem_mood_o)
                  end if
#endif
			    call volume_solution_write
			     if (outsurf.eq.1)then
#ifdef gpu
			     !$omp target update from(rec_grads)
#endif

			    call surface_solution_write
			    end if
			cput1=mpi_wtime()
			end if
			
			if (initcond.eq.95)then           
			if (abs(t - ((idnint(t/output_freq)) * output_freq)).le.tolsmall) then
#ifdef gpu
                  !$omp target update from(u_c_val)
                   if (allocated(u_ct_val)) then
                  !$omp target update from(u_ct_val)
                  end if
                  if (allocated(ielem_vortex)) then
                  !$omp target update from(ielem_vortex)
                  end if
                  if (allocated(ielem_troubled)) then
                  !$omp target update from(ielem_troubled)
                  end if
                   if (allocated(ielem_mood_o)) then
                  !$omp target update from(ielem_mood_o)
                  end if
#endif
                call volume_solution_write
			     if (outsurf.eq.1)then
#ifdef gpu
			     !$omp target update from(rec_grads)
#endif
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
#ifdef gpu
                  !$omp target update from(u_c_val)
                   if (allocated(u_ct_val)) then
                  !$omp target update from(u_ct_val)
                  end if
                  if (allocated(ielem_vortex)) then
                  !$omp target update from(ielem_vortex)
                  end if
                  if (allocated(ielem_troubled)) then
                  !$omp target update from(ielem_troubled)
                  end if
                   if (allocated(ielem_mood_o)) then
                  !$omp target update from(ielem_mood_o)
                  end if
#endif
                call volume_solution_write
			     if (outsurf.eq.1)then
#ifdef gpu
			     !$omp target update from(rec_grads)
#endif
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
integer::i,kmaxe
real::cput1,cput2,cput3,cput4,cput5,cput6,cput8,timec3,timec1,timec4,timec8,totv1,totv2,dumetg1,dumetg2
real::dtiv,flort
kmaxe=xmpielrank(n)
kill=0
t=res_time
iscoun=1

every_time=((idnint(t/output_freq)) * output_freq)+output_freq




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



    if (initcond.eq.95)then
        totk=0
        do i=1,kmaxe
            totk=totk+ielem_totvolume(i)*(1.0/2.0)*&
                (((u_c_val(1,2,i)/u_c_val(1,1,i))**2)+((u_c_val(1,3,i)/u_c_val(1,1,i))**2))
        end do
!
        dumetg1=totk
        dumetg2=0.0
        call mpi_barrier(mpi_comm_world,ierror)
        call mpi_allreduce(dumetg1,dumetg2,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierror)
        totk=dumetg2
        if (n.eq.0)then
!  				totv2=totk/((2.0*pi)**3)
! 					if (it.eq.0)then
! 					taylor=totk
! 					end if
            if (it.eq.0)then
                open(73,file='energy.dat',form='formatted',status='new',action='write',position='append')
            else
                open(73,file='energy.dat',form='formatted',status='old',action='write',position='append')
            end if
            write(73,*)t,totk
            close(73)
        end if

        call mpi_barrier(mpi_comm_world,ierror)
    end if

     if ((multispecies.eq.1))then
         if((initcond.eq.405).or.(initcond.eq.411))then
            ! if ( mod(it, 20) .eq. 0)then
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

          if (rungekutta.eq.3)then
          if  (realgas.eq.1)then
          if (rg_relax.eq.1)then
          !$omp do
          do i=1,kmaxe
          call sources_realgas_pi(n,i)

          end do
          !$omp end do
          end if
          end if
          end if

!           if (realgas.eq.1) call normalise_species(n)
          


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




! write output

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

end subroutine time_marching2





!---------------------------------------------------------------------------------------------!
end module advance
