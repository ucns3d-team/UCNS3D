module moodr
use mpiinfo
use translate
use declaration
use memory
use communications
use io
use partition
use library
use transform
use fluxes
use initialisation
use boundary
use recon
use local
use profile
use flow_operations
use gradients
use basis
use prestore
use riemann
use source
use implicit_time
use implicit_fluxes
use omp_lib




implicit none

 contains


subroutine pad_nad(n)
implicit none
integer,intent(in)::n
integer::i,l,ngp,iqp,iex,kmaxe,k,ii,iconsidered
integer::reduce1,nf,lf,rowf
real,dimension(nof_variables)::nad_delta1
integer::pad_true,nad_true
real::relax_mood1,relax_mood2
real,dimension(1:nof_variables)::leftv
real::mp_pinfl,gammal
real,dimension(1:nof_variables)::rightv
real::mp_pinfr,gammar
real,dimension(imaxdegfree+1,1:nof_variables+turbulenceequations+passivescalar)::utemp
real,dimension(1:nof_variables)::utmin,utmax





kmaxe=xmpielrank(n)


if (cascade.eq.1)then
relax_mood1=mood_var1
relax_mood2=mood_var2
end if
if (cascade.eq.2)then
relax_mood1=mood_var3
relax_mood2=mood_var4
end if



if (itestcase.ge.3)then

#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(ii, i, l, ngp, iqp, iex, k, iconsidered, reduce1, nf, &
!$omp&         lf, rowf, nad_delta1, pad_true, nad_true, leftv, &
!$omp&         mp_pinfl, gammal, rightv, mp_pinfr, gammar, utemp, &
!$omp&         utmin, utmax)
#else
!$omp do
#endif
	do ii=1,nof_interior
	i=el_int(ii)
	iconsidered=i
	ielem_mood(i)=0
    reduce1=0
    pad_true=0
    nad_true=0
	
     !1 copy candidate solution at temp variable   
     leftv(1:nof_variables)=u_c_val(4,1:nof_variables,i)

     !2 transform conservative  to primitive and check if pressure and density are physically admissible if not pad_true=1
                                                if (dimensiona.eq.3)then
                                                
                                                call cons2prim(n,leftv,mp_pinfl,gammal)
						
						!
                                                    if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then						
                                                    pad_true=1
                                                    end if
                                                    if ((leftv(5).le.zero).or.(leftv(5).ne.leftv(5)))then						
                                                    pad_true=1
                                                    end if
                                                else
                                                    call cons2prim(n,leftv,mp_pinfl,gammal)
                                                    if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then						
                                                    pad_true=1
                                                    end if
                                                    if ((leftv(4).le.zero).or.(leftv(4).ne.leftv(4)))then						
                                                    pad_true=1
                                                    end if
                                                
                                                
                                                end if
 
 
        !3 the ones with a physically admissible solution need to get checked
 
		
		if (pad_true.eq.0)then
		
		
		utemp=zero
                
                
                !4 now establish a temporary array with the current solution from the direct side neighbours of considered cell
                
                k=0
			    utemp(1,1:nof_variables)=u_c_val(3,1:nof_variables,i)


			    leftv(1:nof_variables)=utemp(1,1:nof_variables)
			    call cons2prim(n,leftv,mp_pinfl,gammal)
			    utemp(1,1:nof_variables)=leftv(1:nof_variables)

			    k=1
			    do l=1,ielem_ifca(i)
                k=k+1
                utemp(k,1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))

                leftv(1:nof_variables)=utemp(k,1:nof_variables)
			    call cons2prim(n,leftv,mp_pinfl,gammal)
			    utemp(k,1:nof_variables)=leftv(1:nof_variables)


                end do
			    

                
                            
                !5 now establish the min and max bounds
                        do iex=1,nof_variables
                            utmin(iex)=minval(utemp(1:k,iex))
                            utmax(iex)=maxval(utemp(1:k,iex))
                        end do
			
			
			
			! mood_var1=0.0001;mood_var2=0.001; mood_mode=relaxed
			
			
			 do iex=1,nof_variables
			nad_delta1(iex)=max(relax_mood1,(relax_mood2)*(utmax(iex)-utmin(iex)))
			end do

		 
                        nad_true=0
                        !6 specify relaxed or original mood pattern
                        if (mood_mode.gt.0)then

                        leftv(1:nof_variables)=u_c_val(4,1:nof_variables,i)
                        call cons2prim(n,leftv,mp_pinfl,gammal)

                        do iex=1,nof_variables
                        
!                          if ((iex.eq.1).or.(iex.eq.nof_variables))then
                         if ((leftv(iex).lt.(utmin(iex)-nad_delta1(iex))).or.(leftv(iex).gt.(utmax(iex)+nad_delta1(iex))))then
                            nad_true=1
                         end if
!                         endif
                        
                        end do
                        else

                        leftv(1:nof_variables)=u_c_val(4,1:nof_variables,i)
                        call cons2prim(n,leftv,mp_pinfl,gammal)

                        do iex=1,nof_variables



!                          if ((iex.eq.1).or.(iex.eq.nof_variables))then
                        if ((leftv(iex).lt.(utmin(iex))).or.(leftv(iex).gt.(utmax(iex))))then
                            nad_true=1
!                         end if
                        end if
                        end do
                        
                        
                        
                        end if
                        
                        
		 
		 
		 end if

		 !7 now set the mood flag for each element
		 ielem_mood(i)=nad_true+pad_true
		 
        end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(ii, i, l, ngp, iqp, iex, k, iconsidered, reduce1, nf, &
!$omp&         lf, rowf, nad_delta1, pad_true, nad_true, leftv, &
!$omp&         mp_pinfl, gammal, rightv, mp_pinfr, gammar, utemp, &
!$omp&         utmin, utmax)
#else
!$omp do
#endif
	do ii=1,nof_bounded
	i=el_bnd(ii)
	iconsidered=i
            reduce1=0
		pad_true=0
		nad_true=0
		ielem_mood(i)=0
                                leftv(1:nof_variables)=u_c_val(4,1:nof_variables,i)
						if (dimensiona.eq.3)then
						call cons2prim2(n,leftv,rightv,mp_pinfl,mp_pinfr,gammal,gammar)
                                                    if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then						
                                                    pad_true=1
                                                    end if
                                                    if ((leftv(5).le.zero).or.(leftv(5).ne.leftv(5)))then						
                                                    pad_true=1
                                                    end if
                                                else
                                                    call cons2prim(n,leftv,mp_pinfl,gammal)
                                                    if ((leftv(1).le.zero).or.(leftv(1).ne.leftv(1)))then						
                                                    pad_true=1
                                                    end if
                                                    if ((leftv(4).le.zero).or.(leftv(4).ne.leftv(4)))then						
                                                    pad_true=1
                                                    end if
                                                
                                                
                                                end if
                                                
                                                
                                                
                                    
		
		
		
		
		if (pad_true.eq.0)then
		
		
		utemp=zero
                            k=0
			    utemp(1,1:nof_variables)=u_c_val(3,1:nof_variables,i)
			     leftv(1:nof_variables)=utemp(1,1:nof_variables)
			    call cons2prim(n,leftv,mp_pinfl,gammal)
			    utemp(1,1:nof_variables)=leftv(1:nof_variables)


			    k=1
                            
			    
			    
			    
			    
			    
                                do l=1,ielem_ifca(i)	!faces2
                                            if (ielem_ineighb(l,i).eq.n)then	!my cpu only
                                                        if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                            if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
                                                            k=k+1
                                                            utemp(k,1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
                                                            else
                                                            !not periodic ones in my cpu			  				  
                                                            end if
                                                        else
                                                                k=k+1
                                                                utemp(k,1:nof_variables)=u_c_val(3,1:nof_variables,ielem_ineigh(l,i))
                                                        end if
                                            else	!in other cpus they can only be periodic or mpi neighbours
                                            
                                                            if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                                                if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
                                                                    k=k+1
!                                                                     utemp(k,1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
! !                                                                     (rec_ihexl(1,ielem_indexi(l,i),i),1:nof_variables)
! !
! !
                                                                    nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
                                                                    lf=rec_ihexl(1,ielem_indexi(L,i),i)
                                                                    rowf=halo_offset(nf) + lf - 1
                                                                    utemp(k,1:nof_variables)=solhir(rowf,1:nof_variables)




                                                                end if
                                                            else
                                                            
                                                                    
                                                                    
                                                                    k=k+1
!                                                                     utemp(k,1:nof_variables)=iexsolhir(rec_ihexn(1,ielem_indexi(l,i),i))%sol&
!                                                                     (rec_ihexl(1,ielem_indexi(l,i),i),1:nof_variables)
!
                                                                    nf=rec_ihexn(1,ielem_indexi(L,i),rec_local(i))
                                                                    lf=rec_ihexl(1,ielem_indexi(L,i),i)
                                                                    rowf=halo_offset(nf) + lf - 1
                                                                    utemp(k,1:nof_variables)=solhir(rowf,1:nof_variables)


                                                                    
                                                                
                                                            end if
                                        end if


                                        leftv(1:nof_variables)=utemp(k,1:nof_variables)
                                        call cons2prim(n,leftv,mp_pinfl,gammal)
                                        utemp(k,1:nof_variables)=leftv(1:nof_variables)



                                end do
                
                            
		 
                        do iex=1,nof_variables
			utmin(iex)=minval(utemp(1:k,iex))
			utmax(iex)=maxval(utemp(1:k,iex))
			end do
			
			 do iex=1,nof_variables
			nad_delta1(iex)=max(relax_mood1,(relax_mood2)*(utmax(iex)-utmin(iex)))
			end do
		 
		 
                        nad_true=0
                        
                        
                        if (mood_mode.gt.0)then
                        leftv(1:nof_variables)=u_c_val(4,1:nof_variables,i)
                        call cons2prim(n,leftv,mp_pinfl,gammal)


                        do iex=1,nof_variables
                        
!                           if ((iex.eq.1).or.(iex.eq.nof_variables))then
                         if ((leftv(iex).lt.(utmin(iex)-nad_delta1(iex))).or.(leftv(iex).gt.(utmax(iex)+nad_delta1(iex))))then
                            nad_true=1
!                          end if
                        end if
                        
                        end do
                        else
                        leftv(1:nof_variables)=u_c_val(4,1:nof_variables,i)
                        call cons2prim(n,leftv,mp_pinfl,gammal)

                        do iex=1,nof_variables
!                          if ((iex.eq.1).or.(iex.eq.nof_variables))then
                        if ((leftv(iex).lt.(utmin(iex))).or.(leftv(iex).gt.(utmax(iex))))then
                            nad_true=1
!                         end if
                        end if
                        end do
                        
                        
                        
                        end if
                        

		 
		 
		 
		 end if
		 
		 
		 
		 
		 
		 ielem_mood(i)=nad_true+pad_true
		 

		 
		 
		 
        end do
        
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end if		
		

		

end subroutine pad_nad






subroutine mood_operator_2(n)
implicit none
integer::i,kmaxe
integer,intent(in)::n

kmaxe=xmpielrank(n)


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  ielem_recalc(i)=0
  ielem_mood(i)=0
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

  cascade=1
 call pad_nad(n)

 call exhboundhigher_mood(n)
 
 call fix_list(n)
 
 call muscl(n)
 
 call exhboundhigher(n)
 
 
 
 if (dimensiona.eq.2)then
 call calculate_fluxeshi_convective2d_mood(n)

 else
 
 call calculate_fluxeshi_convective_mood(n)
 
 end if
 
 
 
 
#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
    if (ielem_mood(i).eq.0)then
    ielem_mood_o(i)=iorder+1  !target polynomial/scheme admissible
    else
    ielem_mood_o(i)=2         !second order muscl admissible
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

end subroutine mood_operator_2





subroutine mood_operator_1(n)
implicit none
integer::i,kmaxe
integer,intent(in)::n
kmaxe=xmpielrank(n)


#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
  ielem_recalc(i)=0
  ielem_mood(i)=0
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif

 cascade=2
 call pad_nad(n)

 call exhboundhigher_mood(n)
 
 call fix_list(n)
 
 
 if (dimensiona.eq.2)then
 call calculate_fluxeshi_convective2d_mood(n)
 else
 call calculate_fluxeshi_convective_mood(n)
 
 
 end if
 
 
 
 
 

#ifdef gpu
!$omp target teams distribute parallel do
#else
!$omp do
#endif
do i=1,kmaxe
    if (ielem_mood(i).ge.1)then
    ielem_mood_o(i)=1     ! only first order solution admissible
    end if
end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif


end subroutine mood_operator_1


subroutine fix_list(n)
	implicit none
	integer,intent(in)::n
	real::godflux2,sum_detect
	integer::i,l,ngp,kmaxe,iqp,nfx,lfx,rowfx
	real,dimension(numberofpoints2)::weights_temp
	integer,dimension(1)::mright
	kmaxe=xmpielrank(n)
	
#ifdef gpu
!$omp target teams distribute parallel do  &
!$omp& private(i, godflux2, sum_detect, l, ngp, iqp, nfx, lfx, rowfx, &
!$omp&         weights_temp, mright)
#else
!$omp do
#endif
	do i=1,kmaxe
                ielem_recalc(i)=0
		if (ielem_mood(i).ge.1)then
		
                    ielem_recalc(i)=1
                else
		if (ielem_interior(i).eq.0)then
		    mright(1)=zero
		    
		    do l=1,ielem_ifca(i)
				 
				 
				      mright(1)=ielem_mood(ielem_ineigh(l,i))
				      
				  if (mright(1).ge.1)then
                                    ielem_recalc(i)=1
                                    end if	    

		    end do
		end if
		
		
		if (ielem_interior(i).eq.1)then
		    mright(1)=zero
		    
		    do l=1,ielem_ifca(i)
				      
				 
			 
				      
					    if (ielem_ineighb(l,i).eq.n)then	!my cpu only
							if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
								  if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in my cpu
								 mright(1)=ielem_mood(ielem_ineigh(l,i))
 								     
								  else
                                                                        
								  end if
							else
							      mright(1)=ielem_mood(ielem_ineigh(l,i))
!  							       
							end if
					    else	!in other cpus they can only be periodic or mpi neighbours
						
                                if (ielem_ibounds(l,i).gt.0)then	!check for boundaries
                                    if (ibound_icode(ielem_ibounds(l,i)).eq.5)then	!periodic in other cpu
    ! 									  cright(1)=iexboundhir(ielem_ineighn(l,i))%facesol_m(ielem_qface(l,1,i),1)
                                        nfx  = ielem_ineighn(l,i)
                                        lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
                                        rowfx = bound_offset(nfx) + lfx - 1
                                        mright(1) = boundhirm(rowfx)

                                    end if
                                else
    ! 								  cright(1)=iexboundhir(ielem_ineighn(l,i))%facesol_m(ielem_qface(l,1,i),1)
                                    nfx  = ielem_ineighn(l,i)
                                    lfx = ielem_qface(l,1,ielem_inter_id(ielem_indexf(i)))
                                    rowfx = bound_offset(nfx) + lfx - 1
                                    mright(1) = boundhirm(rowfx)

    !
                                end if
					    end if
				      
				     	  if (mright(1).ge.1)then
                                    ielem_recalc(i)=1
                                    end if	    
!  				      
				 
				    
		    end do
		end if
 		end if		
	end do
#ifdef gpu
!$omp end target teams distribute parallel do
#else
!$omp end do
#endif
	end subroutine fix_list








end module moodr
