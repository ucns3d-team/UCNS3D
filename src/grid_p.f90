module library
!> @brief
!> this module includes all the subroutines related to stencils and neighbours establishing
use mpiinfo
use declaration
use io
use transform
implicit none

contains

subroutine tolerances
!> @brief
!> this subroutine specifies the tolerances and other fixed numbers used throughout the code
implicit none

tolsmall=1.0e-13
tolbig=1.0e+13
oo2=1.0d0/2.0d0
zero=0.0d0
! pi=(acos(zero))*2
pi=4.0d0*atan(1.0d0)
alpha=1.0d0
beta=zero

end subroutine tolerances



function determina(eigvl)
!> @brief
!> this function computes the determinant of a matrix
implicit none
real,allocatable,dimension(:,:),intent(in)::eigvl
real,dimension(5,5)::matrix
real::determina
integer::i,j,k,l,nn
real::m,temp
logical :: detexists = .true.
matrix=eigvl
nn=5
    l = 1
    !convert to upper triangular form
    do k = 1, nn-1
        if (matrix(k,k) == 0) then
            detexists = .false.
            do i = k+1, nn
                if (matrix(i,k) /= 0) then
                    do j = 1, nn
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    end do
                    detexists = .true.
                    l=-l
                    exit
                endif
            end do
            if (detexists .eqv. .false.) then
                determina = 0
                return
            end if
        endif
        do j = k+1, nn
            m = matrix(j,k)/matrix(k,k)
            do i = k+1, nn
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            end do
        end do
    end do
   
    !calculate determinant by finding product of diagonal elements
    determina= l
    do i = 1, nn
        determina = determina * matrix(i,i)
    end do








end function determina


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !---------------------------------------------------------------------------------------------!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !---------------------------------------------------------------------------------------------!


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!

subroutine timers(n,cpux1,cpux2,cpux3,cpux4,cpux5,cpux6,timex1,timex2,timex3,timex4,timex5,timex6)
!> @brief
!> this subroutine establishes the timers
real,allocatable,dimension(:),intent(in)::cpux1,cpux2,cpux3,cpux4,cpux5,cpux6
real,allocatable,dimension(:),intent(inout)::timex1,timex2,timex3,timex4,timex5,timex6
integer,intent(in)::n

timex6(1)=cpux6(1)-cpux1(1)
timex5(1)=cpux6(1)-cpux5(1)
timex4(1)=cpux4(1)-cpux3(1)
timex3(1)=cpux3(1)-cpux2(1)
timex2(1)=cpux2(1)-cpux1(1)


end subroutine timers
 
! !---------------------------------------------------------------------------------------------!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




















subroutine xmpifind(xmpie,xmpin,xmpielrank,xmpinrank,imaxe,imaxn,nproc)
!> @brief
!> this subroutine finds the number of elements in each process
implicit none
integer,allocatable,dimension(:),intent(inout)::xmpie
integer,allocatable,dimension(:),intent(inout)::xmpin
integer,allocatable,dimension(:),intent(inout)::xmpielrank
integer,allocatable,dimension(:),intent(inout)::xmpinrank
integer::i,j,k,l,m,itemr
integer,intent(in)::nproc,imaxe,imaxn



k=0
	do i=1,imaxe
		if (xmpie(i).eq.n)then
		k=k+1
		
		
		end if
	end do
	
	xmpielrank(n)=k



end subroutine xmpifind




subroutine globalist(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh,glneighper)
!> @brief
!> this subroutine establishes the connectivity within each process
implicit none
integer,allocatable,dimension(:),intent(in)::xmpie
integer,allocatable,dimension(:),intent(inout)::xmpil
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer,intent(in)::imaxe,n,isize
real,allocatable,dimension(:,:),intent(inout)::centerr
integer,allocatable,dimension(:,:),intent(inout)::glneigh,glneighper
real,allocatable,dimension(:,:)::centerx
integer,allocatable,dimension(:,:)::glneighx
integer,allocatable,dimension(:,:)::glneights,glneightr
integer,allocatable,dimension(:)::glneightot
real,allocatable,dimension(:,:)::centerts,centertr
integer::i,j,k,kmaxe,icpuid,kj,tempi,tempt
real,dimension(1:dimensiona)::cords
real::xv,yc,zc
	if (n.eq.0) then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)"global number of elements",imaxe
	close(63)
	end if

if (dimensiona.eq.3)then

if ((typesten.gt.1).or.(icompact.ge.1))then
 allocate(centerr(1:imaxe,1:3))
end if
allocate(glneigh(1:imaxe,1:6))
allocate(glneighper(1:imaxe,1:6))

if ((typesten.gt.1).or.(icompact.ge.1))then
  centerr(:,:)=-tolbig
end if
 glneigh(:,:)=0
 glneighper(:,:)=0


 call mpi_barrier(mpi_comm_world,ierror)

 kmaxe=xmpielrank(n)
 

 do k=1,kmaxe
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	call compute_centre3d(k,cords)
	
 	centerr(ielem_ihexgl(k),1)=cords(1)
 	centerr(ielem_ihexgl(k),2)=cords(2)
 	centerr(ielem_ihexgl(k),3)=cords(3)
	end if
 	do j=1,ielem_ifca(k)
		

 		glneigh(ielem_ihexgl(k),j)=ielem_ineighg(j,k)

		
	end do
end do

	kj=0
	do k=1,kmaxe
		do j=1,ielem_ifca(k)
			if (ielem_interior(k).eq.0)then
			if ((glneigh(ielem_ihexgl(k),j).eq.0))then
				kj=kj+1
			end if
			else
			if ((glneigh(ielem_ihexgl(k),j).eq.0).and.(ielem_ibounds(j,k).eq.0))then
				kj=kj+1
			end if



			end if
			
		end do
	end do
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)kj,"=number of unknown neighbours before communication"
	close(63)
	end if
	
	
	icpuid=n

 call mpi_barrier(mpi_comm_world,ierror)



  allocate(glneightot(2))
  allocate(glneights(1:kmaxe,1:7))
  
  if ((typesten.gt.1).or.(icompact.ge.1))then
   allocate(centerts(1:kmaxe,1:3))
  end if
  glneights(:,:)=0
  

  do k=1,kmaxe
	  glneights(k,1)=ielem_ihexgl(k)
	  do j=1,ielem_ifca(k)
	  glneights(k,1+j)=glneigh(ielem_ihexgl(k),j)
	  end do
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:3)=centerr(ielem_ihexgl(k),1:3)
	  end if
  end do




!  !first

	do i=0,isize-1
		if (i.ne.n)then
			
			tempi=kmaxe
			call mpi_sendrecv(tempi,1,mpi_integer,i,icpuid,&
			tempt,1,mpi_integer,i,i,mpi_comm_world,status,ierror)
			allocate(glneightr(tempt,1:7))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  allocate(centertr(tempt,1:3))
			  end if
			glneightr(:,:)=0
			 if ((typesten.gt.1).or.(icompact.ge.1))then
			centertr(:,:)=0.0
			end if
			tempi=kmaxe
			
			call mpi_sendrecv(glneights(1:tempi,1:7),tempi*7,mpi_integer,i,icpuid,&
			glneightr(1:tempt,1:7),tempt*7,mpi_integer,i,i,mpi_comm_world,status,ierror)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			call mpi_sendrecv(centerts(1:tempi,1:3),tempi*3,mpi_double_precision,i,icpuid,&
			centertr(1:tempt,1:3),tempt*3,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
			end if
			
			do k=1,tempt
				      do j=1,6
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:3)=centertr(k,1:3)
					    end if
					  end if
				      end do
				      
			end do
			deallocate(glneightr)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deallocate(centertr)
			  end if


			
		end if
	  end do

    deallocate(glneightot)
  deallocate(glneights)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deallocate(centerts)
      end if


end if
if (dimensiona.eq.2)then

if ((typesten.gt.1).or.(icompact.ge.1))then
 allocate(centerr(1:imaxe,1:2))
end if
allocate(glneigh(1:imaxe,1:4))
if ((typesten.gt.1).or.(icompact.ge.1))then
  centerr(:,:)=-tolbig
end if
 glneigh(:,:)=0
 

 call mpi_barrier(mpi_comm_world,ierror)

 kmaxe=xmpielrank(n)
 do k=1,kmaxe
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	call compute_centre2d(k,cords)
	
 	centerr(ielem_ihexgl(k),1)=cords(1)
 	centerr(ielem_ihexgl(k),2)=cords(2)
 	
	end if
 	do j=1,ielem_ifca(k)
		glneigh(ielem_ihexgl(k),j)=ielem_ineighg(j,k)
	end do	

end do
	kj=0
	do k=1,kmaxe
		do j=1,ielem_ifca(k)
			if (ielem_interior(k).eq.0)then
			if ((glneigh(ielem_ihexgl(k),j).eq.0))then
				kj=kj+1
			end if
			else
			if ((glneigh(ielem_ihexgl(k),j).eq.0).and.(ielem_ibounds(j,k).eq.0))then
				kj=kj+1
			end if



			end if
			
		end do
	end do
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)kj,"=number of unknown neighbours before communication"
	close(63)
	end if
	
	
	icpuid=n

 call mpi_barrier(mpi_comm_world,ierror)



  allocate(glneightot(2))
  allocate(glneights(1:kmaxe,1:5))
  if ((typesten.gt.1).or.(icompact.ge.1))then
   allocate(centerts(1:kmaxe,1:2))
  end if
  glneights(:,:)=0
  do k=1,kmaxe
	  glneights(k,1)=ielem_ihexgl(k)
	  do j=1,ielem_ifca(k)
	  glneights(k,1+j)=glneigh(ielem_ihexgl(k),j)
	  end do
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:2)=centerr(ielem_ihexgl(k),1:2)
	  end if
  end do




!  !first

	do i=0,isize-1
		if (i.ne.n)then
			
			tempi=kmaxe
			call mpi_sendrecv(tempi,1,mpi_integer,i,icpuid,&
			tempt,1,mpi_integer,i,i,mpi_comm_world,status,ierror)
			allocate(glneightr(tempt,1:5))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  allocate(centertr(tempt,1:2))
			  end if
			glneightr(:,:)=0
			 if ((typesten.gt.1).or.(icompact.ge.1))then
			centertr(:,:)=0.0
			end if
			tempi=kmaxe
			
			call mpi_sendrecv(glneights(1:tempi,1:5),tempi*5,mpi_integer,i,icpuid,&
			glneightr(1:tempt,1:5),tempt*5,mpi_integer,i,i,mpi_comm_world,status,ierror)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			call mpi_sendrecv(centerts(1:tempi,1:2),tempi*2,mpi_double_precision,i,icpuid,&
			centertr(1:tempt,1:2),tempt*2,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
			end if
			
			do k=1,tempt
				      do j=1,4
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:2)=centertr(k,1:2)
					    end if
					  end if
				      end do
				      
			end do
			deallocate(glneightr)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deallocate(centertr)
			  end if


			
		end if
	  end do

    deallocate(glneightot)
  deallocate(glneights)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deallocate(centerts)
      end if


end if


call mpi_barrier(mpi_comm_world,ierror)

end subroutine globalist

subroutine globalist2(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh,glneighper)
!> @brief
!> this subroutine establishes the connectivity with elements from other processes
implicit none
integer,allocatable,dimension(:),intent(in)::xmpie
integer,allocatable,dimension(:),intent(inout)::xmpil
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer,intent(in)::imaxe,n,isize
real,allocatable,dimension(:,:),intent(inout)::centerr
integer,allocatable,dimension(:,:),intent(inout)::glneigh,glneighper
real,allocatable,dimension(:,:)::centerx
integer,allocatable,dimension(:,:)::glneighx
integer,allocatable,dimension(:,:)::glneights,glneightr,glneightsper,glneightrper
integer,allocatable,dimension(:)::glneightot
real,allocatable,dimension(:,:)::centerts,centertr
integer::i,j,k,kmaxe,icpuid,kj,tempi,tempt
real::xv,yc,zc
real,dimension(1:dimensiona)::cords


 call mpi_barrier(mpi_comm_world,ierror)

  if (dimensiona.eq.3)then

 kmaxe=xmpielrank(n)
		 do k=1,kmaxe
			if ((typesten.gt.1).or.(icompact.ge.1))then
			
			call compute_centre3d(k,cords)
			
		 	centerr(ielem_ihexgl(k),1)=cords(1)
		 	centerr(ielem_ihexgl(k),2)=cords(2)
		 	centerr(ielem_ihexgl(k),3)=cords(3)
			end if
		 	do j=1,ielem_ifca(k)
				

		 		glneigh(ielem_ihexgl(k),j)=ielem_ineighg(j,k)
			if (ielem_interior(k).eq.1) then
				if (ielem_ibounds(j,k).gt.0) then
				if (ibound_icode(ielem_ibounds(j,k)).eq.5) then
				glneighper(ielem_ihexgl(k),j)=1
				else if (ibound_icode(ielem_ibounds(j,k)).eq.50) then
				glneighper(ielem_ihexgl(k),j)=2
				end if
				end if
				end if
				
			end do
		end do

	kj=0
				do k=1,kmaxe
					do j=1,ielem_ifca(k)
						if (ielem_interior(k).eq.0)then
						if ((glneigh(ielem_ihexgl(k),j).eq.0))then
							kj=kj+1
						end if
						else
						if ((glneigh(ielem_ihexgl(k),j).eq.0).and.(ielem_ibounds(j,k).eq.0))then
							kj=kj+1
						end if



						end if
						
					end do
				end do
				if (n.eq.0)then
				open(63,file='history.txt',form='formatted',action='write',position='append')
				write(63,*)kj,"number of unknown neighbours after communication"
			      
				close(63)
				end if
	
	
	icpuid=n
 call mpi_barrier(mpi_comm_world,ierror)

 





	    
 allocate(glneightot(2))
  allocate(glneights(1:kmaxe,1:7))
allocate(glneightsper(1:kmaxe,1:7))

      if ((typesten.gt.1).or.(icompact.ge.1))then
   allocate(centerts(1:kmaxe,1:3))
      end if
  glneights(:,:)=0
glneightsper(:,:)=0
  do k=1,kmaxe
	  glneights(k,1)=ielem_ihexgl(k)
	  glneightsper(k,1)=ielem_ihexgl(k)
	  do j=1,ielem_ifca(k)
	  glneights(k,1+j)=glneigh(ielem_ihexgl(k),j)
	  if (ielem_interior(k).eq.1) then
        if  (ielem_ibounds(j,k).gt.0)then
            if((ibound_icode(ielem_ibounds(j,k)).eq.5).or.(ibound_icode(ielem_ibounds(j,k)).eq.50)) then
            if (ibound_icode(ielem_ibounds(j,k)).eq.5) then
            glneightsper(k,1+j)=1
            else if (ibound_icode(ielem_ibounds(j,k)).eq.50) then
            glneightsper(k,1+j)=2
            end if
        end if
	  end if
	  end if
	  end do
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:3)=centerr(ielem_ihexgl(k),1:3)
	  end if
  end do




 
	do i=0,isize-1
		if (i.ne.n)then


			tempi=kmaxe
			call mpi_sendrecv(tempi,1,mpi_integer,i,icpuid,&
			tempt,1,mpi_integer,i,i,mpi_comm_world,status,ierror)
			allocate(glneightr(tempt,1:7))
			allocate(glneightrper(tempt,1:7))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
 			  allocate(centertr(tempt,1:3))
			  centertr(:,:)=0.0
			  end if
			glneightr(:,:)=0
			glneightrper(:,:)=0
! 			
			
			call mpi_sendrecv(glneights(1:tempi,1:7),tempi*7,mpi_integer,i,icpuid,&
			glneightr(1:tempt,1:7),tempt*7,mpi_integer,i,i,mpi_comm_world,status,ierror)
			call mpi_sendrecv(glneightsper(1:tempi,1:7),tempi*7,mpi_integer,i,icpuid,&
			glneightrper(1:tempt,1:7),tempt*7,mpi_integer,i,i,mpi_comm_world,status,ierror)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			call mpi_sendrecv(centerts(1:tempi,1:3),tempi*3,mpi_double_precision,i,icpuid,&
			centertr(1:tempt,1:3),tempt*3,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
			end if

			do k=1,tempt
				      do j=1,6
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    glneighper(glneightr(k,1),j)=glneightrper(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:3)=centertr(k,1:3)
					    end if
					  end if
				      end do
				      
			end do
			deallocate(glneightr)
			deallocate(glneightrper)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deallocate(centertr)
			  end if



			
		end if
	  end do
    deallocate(glneightot)
  deallocate(glneights)
  deallocate(glneightsper)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deallocate(centerts)
      end if

else



 kmaxe=xmpielrank(n)
 do k=1,kmaxe
	if ((typesten.gt.1).or.(icompact.ge.1))then
	
	call compute_centre2d(k,cords)
	
 	centerr(ielem_ihexgl(k),1)=cords(1)
 	centerr(ielem_ihexgl(k),2)=cords(2)
 	
	end if
 	do j=1,ielem_ifca(k)
		

 		glneigh(ielem_ihexgl(k),j)=ielem_ineighg(j,k)

		
	end do
end do

	kj=0
	do k=1,kmaxe
		do j=1,ielem_ifca(k)
			if (ielem_interior(k).eq.0)then
			if ((glneigh(ielem_ihexgl(k),j).eq.0))then
				kj=kj+1
			end if
			else
			if ((glneigh(ielem_ihexgl(k),j).eq.0).and.(ielem_ibounds(j,k).eq.0))then
				kj=kj+1
			end if



			end if
			
		end do
	end do
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)kj,"number of unknown neighbours after communication"
	close(63)
	end if
	
	
	icpuid=n
 call mpi_barrier(mpi_comm_world,ierror)

 





	    
 allocate(glneightot(2))
  allocate(glneights(1:kmaxe,1:5))
      if ((typesten.gt.1).or.(icompact.ge.1))then
   allocate(centerts(1:kmaxe,1:2))
      end if
  glneights(:,:)=0
  do k=1,kmaxe
	  glneights(k,1)=ielem_ihexgl(k)
	  do j=1,ielem_ifca(k)
	  glneights(k,1+j)=glneigh(ielem_ihexgl(k),j)
	  end do
	  if ((typesten.gt.1).or.(icompact.ge.1))then
	  centerts(k,1:2)=centerr(ielem_ihexgl(k),1:2)
	  end if
  end do




 
	do i=0,isize-1
		if (i.ne.n)then


			tempi=kmaxe
			call mpi_sendrecv(tempi,1,mpi_integer,i,icpuid,&
			tempt,1,mpi_integer,i,i,mpi_comm_world,status,ierror)
			allocate(glneightr(tempt,1:5))
			  if ((typesten.gt.1).or.(icompact.ge.1))then
 			  allocate(centertr(tempt,1:2))
			  centertr(:,:)=0.0
			  end if
			glneightr(:,:)=0
! 			
			
			call mpi_sendrecv(glneights(1:tempi,1:5),tempi*5,mpi_integer,i,icpuid,&
			glneightr(1:tempt,1:5),tempt*5,mpi_integer,i,i,mpi_comm_world,status,ierror)

			if ((typesten.gt.1).or.(icompact.ge.1))then
			call mpi_sendrecv(centerts(1:tempi,1:2),tempi*2,mpi_double_precision,i,icpuid,&
			centertr(1:tempt,1:2),tempt*2,mpi_double_precision,i,i,mpi_comm_world,status,ierror)
			end if

			do k=1,tempt
				      do j=1,4
					  if (glneightr(k,j+1).gt.0)then
					    glneigh(glneightr(k,1),j)=glneightr(k,j+1)
					    if ((typesten.gt.1).or.(icompact.ge.1))then
					    centerr(glneightr(k,1),1:2)=centertr(k,1:2)
					    end if
					  end if
				      end do
				      
			end do
			deallocate(glneightr)
			  if ((typesten.gt.1).or.(icompact.ge.1))then
			  deallocate(centertr)
			  end if



			
		end if
	  end do
    deallocate(glneightot)
  deallocate(glneights)
      if ((typesten.gt.1).or.(icompact.ge.1))then
   deallocate(centerts)
      end if

end if


end subroutine globalist2




subroutine globalistx(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh)
!> @brief
!> this subroutine establishes the connectivity within each process without using a global list
implicit none
integer,allocatable,dimension(:),intent(in)::xmpie
integer,allocatable,dimension(:),intent(inout)::xmpil
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer,intent(in)::imaxe,n,isize
real,allocatable,dimension(:,:),intent(inout)::centerr
integer,allocatable,dimension(:,:),intent(inout)::glneigh
real,allocatable,dimension(:,:)::centerx
integer,allocatable,dimension(:,:)::glneighx
integer,allocatable,dimension(:,:)::glneights,glneightr
integer,allocatable,dimension(:)::glneightot
real,allocatable,dimension(:,:)::centerts,centertr
integer::i,j,k,kmaxe,icpuid,kj,tempi,tempt
real::xv,yc,zc
kmaxe=xmpielrank(n)
kj=0
	do k=1,kmaxe
		do j=1,ielem_ifca(k)
			if (ielem_interior(k).eq.0)then
			  if (ielem_ineighg(j,k).eq.0)then
				  kj=kj+1
			  end if
			else
			  if ((ielem_ineighg(j,k).eq.0).and.(ielem_ibounds(j,k).eq.0))then
				  kj=kj+1
			  end if
			end if
		end do
	end do
	
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)kj,"=number of unknown neighbours before communication"
	close(63)
	end if



end subroutine globalistx


subroutine globalistx2(n,xmpie,xmpil,xmpielrank,imaxe,isize,centerr,glneigh)
!> @brief
!> this subroutine establishes the connectivity with elements from other processes without a global list
implicit none
integer,allocatable,dimension(:),intent(in)::xmpie
integer,allocatable,dimension(:),intent(inout)::xmpil
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer,intent(in)::imaxe,n,isize
real,allocatable,dimension(:,:),intent(inout)::centerr
integer,allocatable,dimension(:,:),intent(inout)::glneigh
real,allocatable,dimension(:,:)::centerx
integer,allocatable,dimension(:,:)::glneighx
integer,allocatable,dimension(:,:)::glneights,glneightr
integer,allocatable,dimension(:)::glneightot
real,allocatable,dimension(:,:)::centerts,centertr
integer::i,j,k,kmaxe,icpuid,kj,tempi,tempt,ifst2,jk
integer:: n_requests
integer,allocatable, dimension(:) :: requests
integer::ifst,imax_cpu, imax_cput,sumcentral,sumcentral_t,approx
real::xv,yc,zc
real,dimension(1:dimensiona)::cords
kmaxe=xmpielrank(n)
kj=0
	do k=1,kmaxe
		do j=1,ielem_ifca(k)
			if (ielem_interior(k).eq.0)then
			  if (ielem_ineighg(j,k).eq.0)then
				  kj=kj+1
			  end if
			else
			  if ((ielem_ineighg(j,k).eq.0).and.(ielem_ibounds(j,k).eq.0))then
				  kj=kj+1
			  end if
			end if
		end do
	end do
	
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)kj,"=number of unknown neighbours after communication"
	close(63)
	end if
	
	allocate(cand(0:isize-1));cand=0		!1 allocate(cand)


	do k=1,kmaxe
		
		do j=1,ielem_ifca(k)
			if (ielem_ineighg(j,k).gt.0)then
! 			
			  if(xmpie(ielem_ineighg(j,k)).ne.n) then
				 cand(xmpie(ielem_ineighg(j,k)))=cand(xmpie(ielem_ineighg(j,k)))+1   
			  end if
			end if
			 
		end do
		
	end do
! 	
	!first establish the cpus that are needed!
	
	
	
	if (iperiodicity.eq.1)then
	  
	sumcentral_t=min(imaxe,iselem*isize)
	approx=min(isize,(sumcentral_t/kmaxe)*iextend)
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)approx,"=number of cpus required for periodic boundary conditions and stencils"
	close(63)
	end if
	end if
	
	
	
	
	ifst=0
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  ifst=ifst+1
! 		      
		  end if
	end do
	
	
	
	
	!first allocate memory for these additional ones
	allocate(cands(0:isize-1),candr(0:isize-1));cands=0;candr=0.0   !2 allocate(cands,candr,neix)
	allocate(dneix1(ifst))
	ifst=0
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  ifst=ifst+1
		      dneix1(ifst)%cpu=cand(i)
		      cands(i)=ifst
		  end if
	end do
	
	!now send to cpus of cpus!
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  cands(i)=ifst
		  end if
	end do
	call mpi_barrier(mpi_comm_world,ierror)
	
	
	
	n_requests = 0
      allocate(requests(2*ifst))					!3 allocate(requests)
	requests(:)=0
	do i=0,isize-1
	
		  if (cand(i).gt.0)then
		      n_requests = n_requests + 1
			      call mpi_isend(                                                     &
			    cands(i), 								& !sendbuf
			    1, mpi_integer,       						& !sendcount, sendtype
			    i, 0,                                       			 & !destination, tag
			    mpi_comm_world, requests(n_requests), ierror                       & !communicator, request handle, error
			)
	
		      n_requests = n_requests + 1
			      call mpi_irecv(                                                     &
				  candr(i),    & !recvbuf
				  1, mpi_integer,          & !recvcount, recvtype
				  i, 0,                                        & !source, tag
				  mpi_comm_world, requests(n_requests), ierror                     & !communicator, request handle, error
			      )
	
	
		  end if
	
	end do


	call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
	
	call mpi_barrier(mpi_comm_world,ierror)
	
	deallocate (requests)							!33 dallocate(requests)
	
	
	
	  imax_cpu=0
	do i=0,isize-1
		if (cand(i).gt.0)then
		imax_cpu=max(imax_cpu,cands(i),candr(i))
		end if
	end do
	
	!imax_cpu=imax_cpu*1.3
	
	call mpi_barrier(mpi_comm_world,ierror)
	
	call mpi_allreduce(imax_cpu,imax_cput,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
	
	
	allocate(candxs(1:ifst,1:imax_cput),candxr(1:ifst,1:imax_cput))		!4 allocate(candxs,candxr)
	candxs=-1
	candxr=-1
! 	
	
	
	ifst=0
	do i=0,isize-1
		  if (cand(i).gt.0)then
		  ifst=ifst+1
		  candxs(:,ifst)=i
		  end if
	end do
	
	
	
	
     
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      call mpi_barrier(mpi_comm_world,ierror)
      
      
      
 	ifst=0
	do i=0,isize-1
	
		  if (cand(i).gt.0)then
		      
		      
		      ifst=ifst+1
 		      call mpi_sendrecv(candxs(ifst:ifst,1:imax_cput)&
,imax_cput,mpi_integer,i,&
      n,candxr(ifst:ifst,1:imax_cput),&
imax_cput,mpi_integer,i,&
      i,mpi_comm_world,status,ierror)
	
	
		  end if
	
	end do


	
	
 	do i=1,ifst
	    do j=1,imax_cput
	    if ((candxr(i,j).ge.0).and.(candxr(i,j).ne.n))then
	     cand(candxr(i,j))=1
	    end if
	    end do
 	end do
	
	ifst2=0
	do i=0,isize-1
	    if (cand(i).gt.0)then
	    ifst2=ifst2+1
! 	    
	    end if
	end do
	
	
	
	
	
	
	
	
	
	
	ifst2=0
	do i=0,isize-1
	    if (cand(i).gt.0)then
	    ifst2=ifst2+1
	    cand(i)=ifst2
	    end if
	end do
	
	
	
	
!
	!now gather from all cpus the required ones
	
	imax_cpu=kmaxe
	imax_cput=0
	
	call mpi_allreduce(imax_cpu,imax_cput,1,mpi_integer,mpi_max,mpi_comm_world,ierror)
	
	
	
	
	if (dimensiona.eq.3)then
	allocate(cand2s(1:imax_cput,6))						!5 allocate(cand2s,cand2r)
	
	allocate(cand2r(1:ifst2,1:imax_cput,1:6))
	allocate(cand2rt(1:imax_cput,1:6))
	cand2s(:,:)=-100
	cand2r(:,:,:)=-200
	
	else

	allocate(cand2s(imax_cput,4))
	allocate(cand2r(ifst2,imax_cput,4))
	allocate(cand2rt(1:imax_cput,1:6))
	cand2s=-100
	cand2r=-200
	
	
	end if
	do i=1,kmaxe
	      do j=1,ielem_ifca(i)
			  cand2s(i,j)=ielem_ineighg(j,i)
	      end do
	end do
	
	
	
    
	
	
	if (dimensiona.eq.3)then
	ifst2=0
	
	do i=0,isize-1
	
		  if (cand(i).gt.0)then
		 	

		  ifst2=ifst2+1
	

		
		
		   call mpi_sendrecv(cand2s(1:imax_cput,1:6),imax_cput*6,mpi_integer,i,n,cand2rt(1:imax_cput,1:6),&
			imax_cput*6,mpi_integer,i,i,mpi_comm_world,status,ierror)

		do jk=1,imax_cput		
		cand2r(ifst2,jk,1:6)=cand2rt(jk,1:6)
		end do

		  end if
	
	end do
	
	else
	ifst2=0
	
	do i=0,isize-1
	
		  if (cand(i).gt.0)then
		  ifst2=ifst2+1
		   call mpi_sendrecv(cand2s(1:imax_cput,1:4)&
,imax_cput*4,mpi_integer,i,&
      n,cand2rt(1:imax_cput,1:4),&
imax_cput*4,mpi_integer,i,&
      i,mpi_comm_world,status,ierror)

		do jk=1,imax_cput		
		cand2r(ifst2,jk,1:4)=cand2rt(jk,1:4)
		end do
		
		  end if
	
	end do
	
	
	end if
	
! 	
	
	
	
	 if ((typesten.gt.1).or.(icompact.ge.1))then
	 
	if (dimensiona.eq.3)then
	allocate(xand2s(imax_cput,3))						!6 allocate(xand2s,xand2r)
	allocate(xand2r(ifst2,imax_cput,3))
	allocate(xand2rt(imax_cput,3))
	do k=1,kmaxe
	
	call compute_centre3d(k,cords)
	
 	xand2s(k,1:3)=cords(1:3)
 	
	end do
	
	else
	allocate(xand2s(imax_cput,2))
	allocate(xand2r(ifst2,imax_cput,2))
	allocate(xand2rt(imax_cput,2))
	do k=1,kmaxe
	
	call compute_centre2d(k,cords)
	
 	xand2s(k,1:2)=cords(1:2)
 	
	end do

	end if

	ifst2=0
	
	do i=0,isize-1
	
		  if (cand(i).gt.0)then
		  ifst2=ifst2+1
		   call mpi_sendrecv(xand2s(1:imax_cput,1:dims)&
,imax_cput*dims,mpi_double_precision,i,&
      n,xand2rt(1:imax_cput,1:dims),&
imax_cput*dims,mpi_double_precision,i,&
      i,mpi_comm_world,status,ierror)

		
		do jk=1,imax_cput		
		xand2r(ifst2,jk,1:dims)=xand2rt(jk,1:dims)
		end do

		  end if
	
	end do
	
	
	sumcentral=imax_cput*ifst2
	
	if (n.eq.0)then
	open(63,file='history.txt',form='formatted',action='write',position='append')
	write(63,*)iselem,"=number of total neighbours needed"
	write(63,*)imax_cput,"=number of neighbours per adjacent cpus"
	write(63,*)sumcentral,"=total number of neighbours per all adjacent cpus"
	close(63)
	end if
	
	
	end if
	call mpi_barrier(mpi_comm_world,ierror)
	
	
	
								    !now deallocate what is not needed
	    deallocate(cands,candr,dneix1,candxs,candxr)
								    
	
	
! 	
	    
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	



end subroutine globalistx2



subroutine globaldea
!> @brief
!> this subroutine deallocates the memory used for establishing the connectivity
implicit none
 call mpi_barrier(mpi_comm_world,ierror)
    if ((typesten.gt.1).or.(icompact.ge.1))then
	    if (lowmem.eq.0)deallocate(xand2s,xand2r,xand2rt)
	    if (lowmem.eq.1)deallocate(centerr)
    end if
	    if (lowmem.eq.0)deallocate(cand2s,cand2r,cand2rt,cand)
	    if (lowmem.eq.1)then
            deallocate(glneigh)

            if (dimensiona.eq.3)deallocate(glneighper)
        end if
 call mpi_barrier(mpi_comm_world,ierror)
 end subroutine globaldea



subroutine xmpilocal
!> @brief
!> this subroutine establishes the local numbering of the cells
implicit none
integer::i,kmaxe,count_block
integer,allocatable,dimension(:)::xmpic,bin,val
allocate(xmpic(0:isize-1))

xmpic=0
kmaxe=xmpielrank(n)

allocate(val(1:kmaxe))
allocate(xgo(1:kmaxe))

do i=1,kmaxe
xgo(i)=ielem_ihexgl(i)
end do

if (n.eq.0)then
allocate(xmpi_re(1:imaxe))
xmpi_re=zero
else
allocate(xmpi_re(1))

end if

do i=1,kmaxe
	xmpil(ielem_ihexgl(i))=i
end do

	
do i=1,imaxe
      xmpic(xmpie(i))=xmpic(xmpie(i))+1
      xmpil(i)=xmpic(xmpie(i))
end do




deallocate(xmpic)


allocate(xmpiall(0:isize-1),offset(0:isize-1))

xmpiall=0

xmpiall(n)=kmaxe



call mpi_allgather(kmaxe,1,mpi_integer,xmpiall,1,mpi_integer,mpi_comm_world,ierror)



offset(0)=0
do i=1,isize-1
offset(i)=offset(i-1)+xmpiall(i-1)
end do




do i=1,kmaxe
  val(i)=ielem_ihexgl(i)
end do





call mpi_gatherv(val(1:kmaxe),kmaxe,mpi_integer,xmpi_re,xmpiall,offset,mpi_integer,0,mpi_comm_world,ierror)






deallocate(val)


!  chunk_n=0
! do i=2,kmaxe
!     if (xgo(i).gt.xgo(i-1)+1)then
! 	chunk_n=chunk_n+1
!     end if
! end do
!   
!   write(100+n,*)chunk_n,n,kmaxe
!   
! allocate(chunk_size(chunk_n))
! 
!  chunk_n=0
!  count_block=0
! do i=2,kmaxe
!     if (xgo(i).gt.xgo(i-1)+1)then
! 	chunk_n=chunk_n+1
! 	chunk_size(chunk_n)=xgo(i)
!     end if
! end do
! 
! do i=1,chunk_n
! write(100+n,*)i,chunk_size(i)
! end do

!  chunk_n=0
! do i=1,kmaxe-1
!       if (xgo(i+1).eq.xgo(i)+1)then
!       
!       else
!       chunk_n=chunk_n+1
!       
!       end if
! end do
! 
! 
! allocate(chunk_size(chunk_n))
!  count_block=0
! do i=1,kmaxe-1
!       if (xgo(i+1).eq.xgo(i)+1)then
!       count_block=count_block+1
!       else
!       chunk_n=chunk_n+1
!       
!       count_block=count_block+1
!       end if
! end do
!     





end subroutine xmpilocal



subroutine count_walls
!> @brief
!> this subroutine allocates the appropriate memory for bounded walls indexing for writing files
implicit none
integer::i,kmaxe,iloop,j
integer,allocatable,dimension(:)::bin,val
kmaxe=xmpielrank(n)
iloop=0
do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		  iloop=iloop+1
	      end if
	  end if
	end do
   end if
end do



allocate(xmpiwall(0:isize-1),woffset(0:isize-1))

xmpiwall=0

xmpiwall(n)=iloop
call mpi_allgather(iloop,1,mpi_integer,xmpiwall,1,mpi_integer,mpi_comm_world,ierror)


woffset(0)=0
do i=1,isize-1
woffset(i)=woffset(i-1)+xmpiwall(i-1)
end do

allocate(val(iloop))
if (n.eq.0)then
allocate(xmpi_wre(1:totwalls))
else
allocate(xmpi_wre(1))

end if


iloop=0

do i=1,kmaxe
  if (ielem_interior(i).eq.1)then
	do j=1,ielem_ifca(i)
	  if (ielem_ibounds(j,i).gt.0)then
	      if ((ibound_icode(ielem_ibounds(j,i)).eq.4).or.(ibound_icode(ielem_ibounds(j,i)).eq.99))then
		  iloop=iloop+1
		  val(iloop)=ibound_inum(ielem_ibounds(j,i))
	      end if
	  end if
	end do
   end if
end do



call mpi_gatherv(val(1:iloop),iloop,mpi_integer,xmpi_wre,xmpiwall,woffset,mpi_integer,0,mpi_comm_world,ierror)




deallocate(val)





end subroutine count_walls




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!








subroutine inverf(r,invr,ivgt)
!> @brief
!> this subroutine inverts some special matrices
  implicit none
  integer,intent(in)::ivgt
  real,intent(in) ::r(1:ivgt-1,1:ivgt-1)
  real,intent(out)::invr(1:ivgt-1,1:ivgt-1)
  real::invvvr2(1:ivgt-1,1:ivgt-1)
 integer i,j,k,gt
 
  invr = 0.d0
 gt=ivgt-1
  do i=gt,1,-1
    invr(i,i) = 1./r(i,i)
	do j=i+1,gt
	  invr(i,j) = 0.
	  do k= 1,j-1
	   invr(i,j) = invr(i,j) - r(k,j)*invr(i,k)
	  enddo
 	  invr(i,j) =invr(i,j) /r(j,j)
	enddo
  enddo

 end subroutine inverf

subroutine allocatevectors(n,tri,invtri,rotvect,vectco,veigl,veigr,rveigl,rveigr,eigvl,eigvr)
!> @brief
!> this subroutine allocates vectors and matrices frequently used
implicit none
integer,intent(in)::n
real,allocatable,dimension(:,:),intent(inout)::tri
real,allocatable,dimension(:,:),intent(inout)::invtri
real,allocatable,dimension(:),intent(inout)::rotvect
real,allocatable,dimension(:),intent(inout)::vectco
real,allocatable,dimension(:),intent(inout)::veigl,veigr,rveigl,rveigr
real,allocatable,dimension(:,:),intent(inout)::eigvl,eigvr
allocate(tri(5,5))
allocate(invtri(5,5))
allocate(eigvl(5,5))
allocate(eigvr(5,5))
allocate(vectco(5+turbulenceequations+passivescalar))
allocate(rotvect(5+turbulenceequations+passivescalar))
allocate(veigl(5))
allocate(veigr(5))
allocate(rveigl(5))
allocate(rveigr(5))
tri=0.0
invtri=0.0
eigvl=0.0
eigvr=0.0
vectco=0.0
rotvect=0.0
veigl=0.0
veigr=0.0
rveigl=0.0
rveigr=0.0
end subroutine allocatevectors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for stencil!!!!!!!!!!!!!!!!!
! subroutine timecpu(tempcp)
! real,intent(inout)::tempcp
! call cpu_time(tempcp)
! end subroutine timecpu
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!subroutine called initially to allocate memory for stencil!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine localstallocation(n,xmpielrank,ilocalstencil,ilocalstencilper,typesten,numneighbours)
!> @brief
!> this subroutine memory for stencil allocation for all cells
	implicit none
	integer,allocatable,dimension(:,:,:,:),intent(inout)::ilocalstencil,ilocalstencilper
	integer,intent(in)::numneighbours
	integer,intent(in)::n
	integer,allocatable,dimension(:),intent(in)::xmpielrank
	integer,intent(in)::typesten
	integer::kmaxe,i
	kmaxe=xmpielrank(n)
	
	allocate (ilocalstencil(n:n,kmaxe,typesten,imaxdegfree+1))
	
	ilocalstencil(n:n,:,:,:)=0
	allocate (ilocalstencilper(n:n,kmaxe,typesten,imaxdegfree+1))
	
	ilocalstencilper(n:n,:,:,:)=0
	
	end subroutine localstallocation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!subroutine called for determining!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!the neighbouring cells at all directions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!kit is initially used for first level!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!neighbours and it is going to be modified!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine find_shape(n,imaxe,ieshape)
!> @brief
!> this subroutine finds the shape of each cell
	implicit none
	integer,intent(inout)::imaxe
	integer,allocatable,dimension(:),intent(inout)::ieshape
	integer::i1,i2,i3,i4,i5,i6,i7,i8,i,j,k,l,m,kk,info,fh,ll,mm,kn
	character(len=12)::celfile,proc
	integer,intent(in)::n
	info=1
	fh=1
	
	if (binio.eq.0)then
	 if (dimensiona.eq.3)then
		write(proc,fmt='(i10)') n
		celfile='GRID.cel'
		
		open(8,file=celfile,form='formatted',status='old',action='read')
	do j=1,imaxe	
		read(8,*)i,i1,i2,i3,i4,i5,i6,i7,i8
		if ((i5.ne.i6).and.(i6.ne.i7).and.(i7.ne.i8).and.(i3.ne.i4))then
		ieshape(j)=1	!hexahedral element
		end if
		if ((i3.eq.i4).and.(i5.eq.i6).and.(i6.eq.i7).and.(i7.eq.i8))then
		ieshape(j)=2	!tetrahedral element
		end if
		if ((i3.ne.i4).and.(i5.eq.i6).and.(i6.eq.i7).and.(i7.eq.i8))then
		ieshape(j)=3	!pyramidal element
		end if
		if ((i3.eq.i4).and.(i5.ne.i6).and.(i7.eq.i8))then
		ieshape(j)=4	!prism element
		end if
	end do	
	else
	write(proc,fmt='(i10)') n
		celfile='GRID.cel'
		open(8,file=celfile,form='formatted',status='old',action='read')
	do j=1,imaxe	
		read(8,*)i,i1,i2,i3,i4
		if ((i3.ne.i4))then
		ieshape(j)=5	!quadrilateral element
		else
		ieshape(j)=6	!triangle el element
		end if
		
	end do	



	end if
	else
	 if (dimensiona.eq.3)then
		write(proc,fmt='(i10)') n
		celfile='GRID.cel'
		
		open(8,file=celfile,form='unformatted',status='old',action='read')
	do j=1,imaxe	
		read(8)i,i1,i2,i3,i4,i5,i6,i7,i8
		if ((i5.ne.i6).and.(i6.ne.i7).and.(i7.ne.i8).and.(i3.ne.i4))then
		ieshape(j)=1	!hexahedral element
		end if
		if ((i3.eq.i4).and.(i5.eq.i6).and.(i6.eq.i7).and.(i7.eq.i8))then
		ieshape(j)=2	!tetrahedral element
		end if
		if ((i3.ne.i4).and.(i5.eq.i6).and.(i6.eq.i7).and.(i7.eq.i8))then
		ieshape(j)=3	!pyramidal element
		end if
		if ((i3.eq.i4).and.(i5.ne.i6).and.(i7.eq.i8))then
		ieshape(j)=4	!prism element
		end if
	end do	
	else
	write(proc,fmt='(i10)') n
		celfile='GRID.cel'
		open(8,file=celfile,form='unformatted',status='old',action='read')
	do j=1,imaxe	
		read(8)i,i1,i2,i3,i4
		if ((i3.ne.i4))then
		ieshape(j)=5	!quadrilateral element
		else
		ieshape(j)=6	!triangle el element
		end if
		
	end do	



	end if
	
	
	
	end if
	
	kn=0
	k=0
	l=0
	m=0
	kk=0
	ll=0
	mm=0
	do j=1,imaxe
	nodes_offset(j)=kn
	if (ieshape(j).eq.1)then	!hexa
		k=k+1
		kn=kn+8
	end if
	if (ieshape(j).eq.2)then	!tetra
		l=l+1
		kn=kn+4
	end if
	if (ieshape(j).eq.3)then	!pyramid
		m=m+1
		kn=kn+5
	end if
	if (ieshape(j).eq.4)then	!prism
		kk=kk+1
		kn=kn+6
	end if
	if (ieshape(j).eq.5)then	!quad
		ll=ll+1
		kn=kn+4
	end if
	if (ieshape(j).eq.6)then	!tri
		mm=mm+1
		kn=kn+3
	end if
	nodes_offset2(j)=kn
	end do
	

close(8)


if (n.eq.0)then
open(63,file='history.txt',form='formatted',action='write',position='append')
	
!-------------------for debugging only -----------------------------------------!
	write(63,*)"------------unstructured fluid dynamics solver ------------"
	write(63,*)"-----------------------grid  type--------------------------"
	if (dimensiona.eq.3 )then
	write(63,*)"-----------------------3d mode--------------------------"
	else
	write(63,*)"-----------------------2d mode--------------------------"
	end if
	!if ((l.gt.0).and.((k.gt.0).or.(m.gt.0).or.(n.gt.0)))then
	write(63,*)k,"hexahedral elements"
	write(63,*)l,"tetrahedral elements"
	write(63,*)m,"pyramidal elements"
	write(63,*)kk,"prismatic elements"
	write(63,*)ll,"quadrilateral elements"
	write(63,*)mm,"triangular elements"

	num_hexas=k
	num_tetras=l
	num_pyramids=m
	num_prisms=kk
	!end if
	

	

close(63)

if ((emetis .eq. 1) .or. (emetis .eq.2)) then
      allnodesgloball=(k*8)+(l*4)+(m*5)+(kk*6)
      else
      allnodesgloball=8*imaxe
    end if
    




end if

	call mpi_bcast(num_hexas, 1, mpi_integer, 0 ,mpi_comm_world,ierror )
	call mpi_bcast(num_tetras, 1, mpi_integer, 0 ,mpi_comm_world,ierror )
	call mpi_bcast(num_pyramids, 1, mpi_integer, 0 ,mpi_comm_world,ierror )
	call mpi_bcast(num_prisms, 1, mpi_integer, 0 ,mpi_comm_world,ierror )


	!-------------------for debugging only -----------------------------------------!
end subroutine find_shape

subroutine fix_offsets_local(n)
!> @brief
!> this assigns the offset for the nodes for each cell in my cpu from the global list
integer::i,j,k,l,m,kk,kmaxe
integer,intent(in)::n


kmaxe=xmpielrank(n)

if (tecplot.eq.5)then
do i=1,kmaxe
	nodes_offset_local(i)=nodes_offset(ielem_ihexgl(i))
	nodes_offset_local2(i)=nodes_offset2(ielem_ihexgl(i))

end do
end if




end subroutine fix_offsets_local






! ---------------------------------------------------------------------------------------------!
subroutine neighbourss(n,imaxe,imaxn,xmpie,xmpin,xmpielrank,restart,dinoder)
!> @brief
!> this subroutine finds the neighbours
implicit none
type(anode_ne),allocatable,dimension(:),intent(inout)::dinoder
integer,intent(in)::n,imaxe,restart,imaxn
integer,allocatable,dimension(:),intent(in)::xmpie,xmpin
integer,allocatable,dimension(:),intent(in)::xmpielrank
integer::i,j,ji,k,lm,iex,kmaxn,kk,kmaxe,ihax1,numnod,l,icn,l2,p,m,q,ij,c_n1,c_n2,c_n3,c_n4,d_n1,d_n2,d_n3,d_n4
integer::it1,it2,it3,it4,it5,it6,it7,it8,itx,inx,itf,ifg,print_i,i1,i2,i3,i4,i5,i6,i7,i8,imbg,icface
real::delta,cpuer
   character(len=20)::proc,neibfile,proc3
   kmaxe=xmpielrank(n)
	
	
	

	
		dinoder2(:)%numberofneib=0
		
	

	
        
	do i = 1, kmaxe

  select case (ielem_ishape(i))

  case (4)  ! prism

    ielem_types_faces(1:2, i) = 6
    ielem_types_faces(3:5, i) = 5
    ielem_nodes_faces(1:5, 1:4, i) = 0

    ielem_ineighg(1:5, i) = 0
    ielem_ifca(i) = 5
    ielem_vdec(i) = 3

    ielem_totvolume(i) = 0.0
    ielem_minedge(i)   = 0.0
    ielem_walldist(i)  = 0.0

    ! first face
    ielem_nodes_faces(1,1,i) = ielem_nodes(5,i)
    ielem_nodes_faces(1,2,i) = ielem_nodes(6,i)
    ielem_nodes_faces(1,3,i) = ielem_nodes(4,i)

    ! second face
    ielem_nodes_faces(2,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(2,2,i) = ielem_nodes(1,i)
    ielem_nodes_faces(2,3,i) = ielem_nodes(3,i)

    ! third face
    ielem_nodes_faces(3,1,i) = ielem_nodes(1,i)   ! 125,154
    ielem_nodes_faces(3,2,i) = ielem_nodes(2,i)
    ielem_nodes_faces(3,3,i) = ielem_nodes(5,i)
    ielem_nodes_faces(3,4,i) = ielem_nodes(4,i)

    ! fourth face
    ielem_nodes_faces(4,1,i) = ielem_nodes(2,i)   ! 236,265
    ielem_nodes_faces(4,2,i) = ielem_nodes(3,i)
    ielem_nodes_faces(4,3,i) = ielem_nodes(6,i)
    ielem_nodes_faces(4,4,i) = ielem_nodes(5,i)

    ! fifth face
    ielem_nodes_faces(5,1,i) = ielem_nodes(3,i)   ! 314,346
    ielem_nodes_faces(5,2,i) = ielem_nodes(1,i)
    ielem_nodes_faces(5,3,i) = ielem_nodes(4,i)
    ielem_nodes_faces(5,4,i) = ielem_nodes(6,i)


  case (1)  ! hexahedral

    ielem_types_faces(1:6, i) = 5
    ielem_vdec(i) = 6

    ielem_totvolume(i) = 0.0
    ielem_minedge(i)   = 0.0
    ielem_walldist(i)  = 0.0

    ielem_nodes_faces(:,:,i) = 0

    ielem_ineighg(:, i) = 0
    ielem_ifca(i) = 6

    ! first face
    ielem_nodes_faces(1,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(1,2,i) = ielem_nodes(1,i)   ! dec214,243
    ielem_nodes_faces(1,3,i) = ielem_nodes(4,i)
    ielem_nodes_faces(1,4,i) = ielem_nodes(3,i)

    ! second face
    ielem_nodes_faces(2,1,i) = ielem_nodes(5,i)
    ielem_nodes_faces(2,2,i) = ielem_nodes(6,i)   ! 567,578
    ielem_nodes_faces(2,3,i) = ielem_nodes(7,i)
    ielem_nodes_faces(2,4,i) = ielem_nodes(8,i)

    ! third face
    ielem_nodes_faces(3,1,i) = ielem_nodes(6,i)
    ielem_nodes_faces(3,2,i) = ielem_nodes(2,i)   ! 623,637
    ielem_nodes_faces(3,3,i) = ielem_nodes(3,i)
    ielem_nodes_faces(3,4,i) = ielem_nodes(7,i)

    ! fourth face
    ielem_nodes_faces(4,1,i) = ielem_nodes(8,i)
    ielem_nodes_faces(4,2,i) = ielem_nodes(4,i)   ! 841,815
    ielem_nodes_faces(4,3,i) = ielem_nodes(1,i)
    ielem_nodes_faces(4,4,i) = ielem_nodes(5,i)

    ! fifth face
    ielem_nodes_faces(5,1,i) = ielem_nodes(7,i)   ! 734,748
    ielem_nodes_faces(5,2,i) = ielem_nodes(3,i)
    ielem_nodes_faces(5,3,i) = ielem_nodes(4,i)
    ielem_nodes_faces(5,4,i) = ielem_nodes(8,i)

    ! sixth face
    ielem_nodes_faces(6,1,i) = ielem_nodes(5,i)   ! 512,526
    ielem_nodes_faces(6,2,i) = ielem_nodes(1,i)
    ielem_nodes_faces(6,3,i) = ielem_nodes(2,i)
    ielem_nodes_faces(6,4,i) = ielem_nodes(6,i)


  case (2)  ! tetrahedral

    ielem_types_faces(1:4, i) = 6
    ielem_nodes_faces(:,:,i) = 0

    ielem_ineighg(:, i) = 0
    ielem_vdec(i) = 1

    ielem_totvolume(i) = 0.0
    ielem_minedge(i)   = 0.0
    ielem_walldist(i)  = 0.0

    ielem_ifca(i) = 4

    ! first face
    ielem_nodes_faces(1,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(1,2,i) = ielem_nodes(4,i)
    ielem_nodes_faces(1,3,i) = ielem_nodes(1,i)

    ! second face
    ielem_nodes_faces(2,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(2,2,i) = ielem_nodes(1,i)
    ielem_nodes_faces(2,3,i) = ielem_nodes(3,i)

    ! third face
    ielem_nodes_faces(3,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(3,2,i) = ielem_nodes(3,i)
    ielem_nodes_faces(3,3,i) = ielem_nodes(4,i)

    ! fourth face
    ielem_nodes_faces(4,1,i) = ielem_nodes(3,i)
    ielem_nodes_faces(4,2,i) = ielem_nodes(1,i)
    ielem_nodes_faces(4,3,i) = ielem_nodes(4,i)


  case (3)  ! pyramidal

    ielem_types_faces(1, i)   = 5
    ielem_types_faces(2:5, i) = 6
    ielem_nodes_faces(:,:,i) = 0

    ielem_ineighg(:, i) = 0
    ielem_vdec(i) = 2

    ielem_totvolume(i) = 0.0
    ielem_minedge(i)   = 0.0
    ielem_walldist(i)  = 0.0

    ielem_ifca(i) = 5

    ! first face
    ielem_nodes_faces(1,1,i) = ielem_nodes(4,i)   ! 432,421
    ielem_nodes_faces(1,2,i) = ielem_nodes(3,i)
    ielem_nodes_faces(1,3,i) = ielem_nodes(2,i)
    ielem_nodes_faces(1,4,i) = ielem_nodes(1,i)

    ! second face
    ielem_nodes_faces(2,1,i) = ielem_nodes(1,i)
    ielem_nodes_faces(2,2,i) = ielem_nodes(2,i)
    ielem_nodes_faces(2,3,i) = ielem_nodes(5,i)

    ! third face
    ielem_nodes_faces(3,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(3,2,i) = ielem_nodes(3,i)
    ielem_nodes_faces(3,3,i) = ielem_nodes(5,i)

    ! fourth face
    ielem_nodes_faces(4,1,i) = ielem_nodes(3,i)
    ielem_nodes_faces(4,2,i) = ielem_nodes(4,i)
    ielem_nodes_faces(4,3,i) = ielem_nodes(5,i)

    ! fifth face
    ielem_nodes_faces(5,1,i) = ielem_nodes(4,i)
    ielem_nodes_faces(5,2,i) = ielem_nodes(1,i)
    ielem_nodes_faces(5,3,i) = ielem_nodes(5,i)


  case (5)  ! quadrilateral

    ielem_ifca(i) = 4

    ielem_nodes_faces(:,:,i) = 0
    ielem_vdec(i) = 2

    ielem_totvolume(i) = 0.0
    ielem_minedge(i)   = 0.0
    ielem_walldist(i)  = 0.0

    ielem_ineighg(:, i) = 0

    ! first face
    ielem_nodes_faces(1,1,i) = ielem_nodes(1,i)
    ielem_nodes_faces(1,2,i) = ielem_nodes(2,i)

    ! second face
    ielem_nodes_faces(2,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(2,2,i) = ielem_nodes(3,i)

    ! third face
    ielem_nodes_faces(3,1,i) = ielem_nodes(3,i)
    ielem_nodes_faces(3,2,i) = ielem_nodes(4,i)

    ! fourth face
    ielem_nodes_faces(4,1,i) = ielem_nodes(4,i)
    ielem_nodes_faces(4,2,i) = ielem_nodes(1,i)


  case (6)  ! triangular

    ielem_ifca(i) = 3

    ielem_nodes_faces(:,:,i) = 0
    ielem_vdec(i) = 1

    ielem_totvolume(i) = 0.0
    ielem_minedge(i)   = 0.0
    ielem_walldist(i)  = 0.0

    ielem_ineighg(:, i) = 0

    ! first face
    ielem_nodes_faces(1,1,i) = ielem_nodes(1,i)
    ielem_nodes_faces(1,2,i) = ielem_nodes(2,i)

    ! second face
    ielem_nodes_faces(2,1,i) = ielem_nodes(2,i)
    ielem_nodes_faces(2,2,i) = ielem_nodes(3,i)

    ! third face
    ielem_nodes_faces(3,1,i) = ielem_nodes(3,i)
    ielem_nodes_faces(3,2,i) = ielem_nodes(1,i)

  end select
end do

	
! 	
       

	do i=1,kmaxe
		do j=1,ielem_nonodes(i)
		

		k=ielem_nodes(j,i)
		dinoder2(k)%numberofneib=dinoder2(k)%numberofneib+1
		
		end do
		
		
	end do
	
	

	
	

	
	do i=1,imaxn
		if (dinoder2(i)%numberofneib.gt.0)then
		allocate(dinoder2(i)%neibids(dinoder2(i)%numberofneib))
		dinoder2(i)%neibids(:)=0
		dinoder2(i)%numberofneib=0
		end if
	end do
	
	
	
! 	
	
	do i=1,kmaxe
		do ji=1,ielem_nonodes(i)
		k=ielem_nodes(ji,i)
		dinoder2(k)%numberofneib=dinoder2(k)%numberofneib+1
		j=dinoder2(k)%numberofneib
		dinoder2(k)%neibids(j)=i
		end do
		
		
	end do
	
! 	


	print_i=kmaxe/5
open(63,file='history.txt',form='formatted',action='write',position='append')	
	if (n.eq. 0) then
	
	write(63,*)'find neigbours'
	end if
close(63)

 

	do ihax1=1,kmaxe


	
	
	
			numnod=ielem_nonodes(ihax1)
	 		nodelist(1:numnod) = ielem_nodes(1:numnod,ihax1)
	 		
			
	

		
				
				list(:)=0
				m=0
				do k=1,numnod
					p=nodelist(k)
						if (dinoder2(p)%numberofneib.gt.0)then
						do q=1,dinoder2(p)%numberofneib
							j=dinoder2(p)%neibids(q)
							icn=0
							do ij=1,m
								if (j.eq.list(ij))then
									icn=1
								end if
							end do
							if ((icn.eq.0).and.(ihax1.ne.j))then
								m=m+1
								list(m)=j
							end if
						end do
						end if
				end do
				
			      
    do l=1,ielem_ifca(ihax1)
  if (dimensiona.eq.3)then

  do p=1,m
  
  j=list(p)
	  if (j.eq.ihax1) cycle
	do l2=1,ielem_ifca(j)
		  
		  if (ielem_types_faces(l,ihax1).eq.ielem_types_faces(l2,j))then
		  if (ielem_types_faces(l,ihax1).eq.6)then
		  c_n4=0
		  c_n1=ielem_nodes_faces(l,1,ihax1)
		  c_n2=ielem_nodes_faces(l,2,ihax1)
		  c_n3=ielem_nodes_faces(l,3,ihax1)
		  d_n4=0
		  d_n1=ielem_nodes_faces(l2,1,j)
		  d_n2=ielem_nodes_faces(l2,2,j)
		  d_n3=ielem_nodes_faces(l2,3,j)
		  if (((c_n1.eq.d_n1).or.(c_n1.eq.d_n2).or.(c_n1.eq.d_n3)).and.&
		      ((c_n2.eq.d_n1).or.(c_n2.eq.d_n2).or.(c_n2.eq.d_n3)).and.&
		      ((c_n3.eq.d_n1).or.(c_n3.eq.d_n2).or.(c_n3.eq.d_n3)))then
			  ielem_ineighg(l,ihax1)=ielem_ihexgl(j)
			 
			  
		  goto 101
		  end if
		  else
		  c_n4=ielem_nodes_faces(l,4,ihax1)
		  c_n1=ielem_nodes_faces(l,1,ihax1)
		  c_n2=ielem_nodes_faces(l,2,ihax1)
		  c_n3=ielem_nodes_faces(l,3,ihax1)
		  d_n4=ielem_nodes_faces(l2,4,j)
		  d_n1=ielem_nodes_faces(l2,1,j)
		  d_n2=ielem_nodes_faces(l2,2,j)
		  d_n3=ielem_nodes_faces(l2,3,j)
		  if (((c_n1.eq.d_n1).or.(c_n1.eq.d_n2).or.(c_n1.eq.d_n3).or.(c_n1.eq.d_n4)).and.&
		      ((c_n2.eq.d_n1).or.(c_n2.eq.d_n2).or.(c_n2.eq.d_n3).or.(c_n2.eq.d_n4)).and.&
		      ((c_n3.eq.d_n1).or.(c_n3.eq.d_n2).or.(c_n3.eq.d_n3).or.(c_n3.eq.d_n4)).and.&
		      ((c_n4.eq.d_n1).or.(c_n4.eq.d_n2).or.(c_n4.eq.d_n3).or.(c_n4.eq.d_n4)))then
			  ielem_ineighg(l,ihax1)=ielem_ihexgl(j)
			  
			  
		  goto 101
		  end if
		  end if

		  
		 
		 
		  
		  
		  end if
	end do
  end do
  else
  do p=1,m
  j=list(p)
	  if (j.eq.ihax1) cycle
	do l2=1,ielem_ifca(ihax1)
		  c_n1=ielem_nodes_faces(l,1,ihax1)
		  c_n2=ielem_nodes_faces(l,2,ihax1)
		  

		  d_n1=ielem_nodes_faces(l2,1,j)
		  d_n2=ielem_nodes_faces(l2,2,j)
		 
		 
		  if (((c_n1.eq.d_n1).or.(c_n1.eq.d_n2)).and.&
		      ((c_n2.eq.d_n1).or.(c_n2.eq.d_n2)))then
			  ielem_ineighg(l,ihax1)=ielem_ihexgl(j)
			  
		  goto 101
		  end if
	end do
  end do
  end if





  
  101 continue
			end do
		end do	
     
	
		
		
		

	ji=0
      do i=1,kmaxe
	    ielem_interior(i)=0
	    do l=1,ielem_ifca(i)
	      if (ielem_ineighg(l,i).eq.0)then
	      ielem_interior(i)=1
	      ji=ji+1
	       end if
	    end do
      end do
     
      
		
end subroutine neighbourss
! !---------------------------------------------------------------------------------------------!







subroutine cons(n,diconr,dicons,iperiodicity,xmpielrank,isize,diconrpa,diconrpm,diconspo,xper,&
yper,zper,diconrpf,numneighbours,typesten)
!> @brief
!> this subroutine establishes the connectivity across different cpus
implicit none
integer,intent(in)::iperiodicity,n,isize,numneighbours,typesten
real::small,dist,temp_x
real,intent(in)::xper,yper,zper
integer::i,k,jjj,kj,j,l,im,io,ir,kmaxe,icpuid,kkj,countxsize
integer,dimension(6)::jx
real,dimension(3)::dumface1,dumface2
integer::faceper1,faceper2,p,kkj2,kkj3,kkj4,kkj5,jj1,p1,p2,p3,p4,q1,q2,q3,q4,j1,j2,j3,j4,code_per1,facex,iconsi,ixxff
integer,allocatable,dimension(:),intent(in)::xmpielrank
type(aconnx),allocatable,dimension(:),intent(inout)::diconr,diconrpa,diconrpm,diconrpf
type(aconnx),allocatable,dimension(:),intent(inout)::dicons,diconspo
real,dimension(1:dimensiona)::cords
real,dimension(1:8,1:dimensiona)::vext


kmaxe=xmpielrank(n)

small=tolsmall					! a small number!
!call mpi_barrier(mpi_comm_world,ierror)

           
kj=0; kkj=0
	do i=1,kmaxe
	if (ielem_interior(i).eq.1)then
		do k=1,ielem_ifca(i)
		
		  if (((ielem_ineighg(k,i).eq.0)))then
		      if(ielem_ibounds(k,i).eq.0)then
						kj=kj+1
		      else
		      if ((ibound_icode(ielem_ibounds(k,i)).eq.5).or.(ibound_icode(ielem_ibounds(k,i)).eq.50))then
						  kj=kj+1	
						
		      end if
		      end if
		end if
		
		

		end do
	end if
	end do
	
	




! allocate(dICONs(1:countxsize))   1st


allocate(diconr(n:n))
allocate(diconr(n)%howmanyi(1))


  diconr(n)%procid=n
  diconr(n)%howmanyi(1)=kj			!this is the total number of neighbours i need

	
	
	
		allocate(diconr(n)%whichi(diconr(n)%howmanyi(1),6))
		diconr(n)%whichi(:,:)=0
! 		
		kj=0
		do i=1,kmaxe
		if (ielem_interior(i).eq.1)then
		io=ielem_ifca(i)
		do k=1,io
		 if (((ielem_ineighg(k,i).eq.0)))then
		    if ((ielem_ibounds(k,i).eq.0))then
						kj=kj+1
				diconr(n)%whichi(kj,1)=i
				diconr(n)%whichi(kj,2)=k
		    
						
		      else
		      
		      if ((ibound_icode(ielem_ibounds(k,i)).eq.5).or.(ibound_icode(ielem_ibounds(k,i)).eq.50))then
						  kj=kj+1
				diconr(n)%whichi(kj,1)=i
				diconr(n)%whichi(kj,2)=k
		      end if
		      end if
		end if
		end do
		end if
		end do

		 




		do kj=1,diconr(n)%howmanyi(1)
			      
			      select case(ielem_types_faces(diconr(n)%whichi(kj,2),diconr(n)%whichi(kj,1)))
			      case(5)
			      diconr(n)%whichi(kj,3:6)=ielem_nodes_faces(diconr(n)%whichi(kj,2),1:4,diconr(n)%whichi(kj,1))


			      case(6)
			      diconr(n)%whichi(kj,3:5)=ielem_nodes_faces(diconr(n)%whichi(kj,2),1:3,diconr(n)%whichi(kj,1))
			       diconr(n)%whichi(kj,6)=0
				end select

			
		end do





do ihax1=1,diconr(n)%howmanyi(1)

		      select case(ielem_types_faces(diconr(n)%whichi(ihax1,2),diconr(n)%whichi(ihax1,1)))
		      
			      case(5)
			      diconr(n)%whichi(ihax1,3:6)=ielem_nodes_faces(diconr(n)%whichi(ihax1,2),1:4,diconr(n)%whichi(ihax1,1))
	 		       nodelist(1:4) = diconr(n)%whichi(ihax1,3:6)
	 		       
			
			p1=nodelist(1)
			    do q1=1,dinoder2(p1)%xne(1)
				  j1=dinoder2(p1)%xneib(q1)
				      p2=nodelist(2)
				      do q2=1,dinoder2(p2)%xne(1)
					  j2=dinoder2(p2)%xneib(q2)
					    if ((j2.eq.j1).and.(j2.ne.ielem_ihexgl(diconr(n)%whichi(ihax1,1))))then
					      p3=nodelist(3)
						  do q3=1,dinoder2(p3)%xne(1)
						      j3=dinoder2(p3)%xneib(q3)
							  if ((j3.eq.j2).and.(j3.ne.ielem_ihexgl(diconr(n)%whichi(ihax1,1))))then
							      p4=nodelist(4)
								  do q4=1,dinoder2(p4)%xne(1)
								    j4=dinoder2(p4)%xneib(q4)
									if ((j4.eq.j3).and.(j4.ne.ielem_ihexgl(diconr(n)%whichi(ihax1,1))))then
									ielem_ineighg(diconr(n)%whichi(ihax1,2),diconr(n)%whichi(ihax1,1))=j4
									go to 331
									end if
								  end do
							  end if
						  end do
					    end if
				     end do
			    end do
		case(6)
			      diconr(n)%whichi(ihax1,3:5)=ielem_nodes_faces(diconr(n)%whichi(ihax1,2),1:3,diconr(n)%whichi(ihax1,1))
	 		       nodelist(1:3) = diconr(n)%whichi(ihax1,3:5)
	 		       
			
			p1=nodelist(1)
			    do q1=1,dinoder2(p1)%xne(1)
				  j1=dinoder2(p1)%xneib(q1)
				      p2=nodelist(2)
				      do q2=1,dinoder2(p2)%xne(1)
					  j2=dinoder2(p2)%xneib(q2)
					    if ((j2.eq.j1).and.(j2.ne.ielem_ihexgl(diconr(n)%whichi(ihax1,1))))then
					      p3=nodelist(3)
						  do q3=1,dinoder2(p3)%xne(1)
						      j3=dinoder2(p3)%xneib(q3)
							  if ((j3.eq.j2).and.(j3.ne.ielem_ihexgl(diconr(n)%whichi(ihax1,1))))then
									ielem_ineighg(diconr(n)%whichi(ihax1,2),diconr(n)%whichi(ihax1,1))=j3
									go to 331
							  end if
						  end do
					    end if
				     end do
			    end do	
	      end select
	      331 continue
				
end do


            kj=0

	      do i=1,diconr(n)%howmanyi(1)
! 		    

		       if ((ielem_ineighg(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1)).eq.0))then
		      if (iperiodicity.eq.1)then
		      if ((((ibound_icode(ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))))).eq.5).or.&
		      (((ibound_icode(ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))))).eq.50))then
						  kj=kj+1
			
		      end if
		      else

					    kj=kj+1
		      end if
		      end if
		end do

!-------------------for debugging only -----------------------------------------!
		if (iperiodicity.eq.1)then
		    ! 		
				    
				    allocate(diconrpm(n:n))
				    allocate(diconrpm(n)%howmanyi(1))
				    diconrpm(n)%howmanyi(1)=kj
				    allocate(diconrpm(n)%whichi(kj,3))
				    diconrpm(n)%whichi(:,:)=0
				    kj=0
		    ! 		

				    


				    kj=0
				    do i=1,diconr(n)%howmanyi(1)


					  if ((ielem_ineighg(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1)).eq.0))then
					  if ((((ibound_icode(ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))))).eq.5).or.&
		      (((ibound_icode(ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))))).eq.50))then
								      kj=kj+1
					    diconrpm(n)%whichi(kj,1)=diconr(n)%whichi(i,1)
					    diconrpm(n)%whichi(kj,3)=ielem_ihexgl(diconr(n)%whichi(i,1))
					    diconrpm(n)%whichi(kj,2)=(diconr(n)%whichi(i,2))
					  end if
					  end if
				    end do
		    ! 		
				    
				    jjj=0
		    !-------------------for debugging only -----------------------------------------!
				    allocate(dICONRPA(n:n))
				    allocate(dICONRPf(1:isize-1))
		    
				    do i=1,isize-1
				    allocate(dICONRPf(i)%howmanyi(1))
				    dICONRPf(i)%howmanyi(1)=0
				    dICONRPf(i)%howmanyi(1)=diconrpm(n)%howmanyi(1)-jjj
				    allocate(dICONRPf(i)%whichi(dICONRPf(i)%howmanyi(1),1))
				    dICONRPf(i)%whichi(:,:)=0
				    end do
				    allocate(dICONRPA(n)%howmanyi(1))
				    dICONRPA(n)%howmanyi(1)=0
				    dICONRPA(n)%howmanyi(1)=diconrpm(n)%howmanyi(1)-jjj
				    allocate(dICONRPA(n)%whichi(diconrpm(n)%howmanyi(1)-jjj,3))
				    allocate(dICONRPA(n)%facx(diconrpm(n)%howmanyi(1)-jjj,3))
				    dICONRPA(n)%whichi(:,:)=0
				    dICONRPA(n)%facx(:,:)=0.d0
				    kj=0
				    do i=1,diconrpm(n)%howmanyi(1)
		   
					    if ((ielem_ineighg(diconrpm(n)%whichi(i,2),diconrpm(n)%whichi(i,1)).eq.0))then
					  if ((((ibound_icode(ielem_ibounds(diconrpm(n)%whichi(i,2),diconrpm(n)%whichi(i,1))))).eq.5).or.&
		      (((ibound_icode(ielem_ibounds(diconrpm(n)%whichi(i,2),diconrpm(n)%whichi(i,1))))).eq.50))then
		    ! 						  kj=kj+1
		    ! 			
					    kj=kj+1
					    dICONRPA(n)%whichi(kj,1)=diconrpm(n)%whichi(i,1)
					    dICONRPA(n)%whichi(kj,2)=diconrpm(n)%whichi(i,2)
					    dICONRPA(n)%whichi(kj,3)=diconrpm(n)%whichi(i,3)
					      if (ielem_types_faces(dICONRPA(n)%whichi(kj,2),dICONRPA(n)%whichi(kj,1)).eq.5)then

					      ixxff=4
					      else
					      ixxff=3
					      end if
						iconsi=dICONRPA(n)%whichi(kj,1); facex=dICONRPA(n)%whichi(kj,2)

					      call compute_centre3df(n,iconsi,facex,ixxff,cords)
					    

					    dICONRPA(n)%facx(kj,1)=cords(1)
					    dICONRPA(n)%facx(kj,2)=cords(2)
					    dICONRPA(n)%facx(kj,3)=cords(3)
					    end if
					    end if
		    ! 			end if
				    end do
				    
				    dumface1(:)=0.0
				    dumface2(:)=0.0
				    faceper1=0
				    faceper2=0
				    icpuid=n
				    call mpi_barrier(mpi_comm_world,ierror)
		   
		    !-------------------for debugging only -----------------------------------------!
				    allocate(diconspo(1:isize-1))
		    
				    do i=1,isize-1
				    allocate(diconspo(i)%howmanythey(1))
				    diconspo(i)%howmanythey(1)=0
				    end do
		    !-------------------for debugging only -----------------------------------------!
		    
		    !-------------------for debugging only -----------------------------------------!
				    k=0
				    do i=0,isize-1
				    if (i.ne.n)then
				    k=k+1
				    diconspo(k)%procid=i
				    end if
				    end do
		    !-------------------for debugging only -----------------------------------------!
		    
		    !-------------------for debugging only -----------------------------------------!
				    call mpi_barrier(mpi_comm_world,ierror)
		    !-------------------for debugging only -----------------------------------------!
		    
		    !-------------------for debugging only -----------------------------------------!
					    icpuid=n
				    k=0
				    do i=0,isize-1
					    if (i.ne.n)then
							    k=k+1
			call mpi_sendrecv(dICONRPA(n)%howmanyi(1),1,mpi_integer,i,icpuid,&
			diconspo(k)%howmanythey(1),1,mpi_integer,diconspo(k)%procid,diconspo(k)%procid,mpi_comm_world,status,ierror)
					    end if
				    end do


				    call mpi_barrier(mpi_comm_world,ierror)
				    do i=1,isize-1
				    allocate(diconspo(i)%whichthey(diconspo(i)%howmanythey(1),1))
				    allocate(diconspo(i)%facx(diconspo(i)%howmanythey(1),3))
				    diconspo(i)%whichthey(:,:)=0
				    diconspo(i)%facx(:,:)=0.d0
				    end do
				    call mpi_barrier(mpi_comm_world,ierror)
					    icpuid=n
				    k=0
		    do i=0,isize-1
		    if (i.ne.n)then
		      k=k+1
		      call mpi_sendrecv(dICONRPA(n)%facx(1:dICONRPA(n)%howmanyi(1),1:3),&
		    dICONRPA(n)%howmanyi(1)*3,mpi_double_precision,i,icpuid,&
		      diconspo(k)%facx(1:diconspo(k)%howmanythey(1),1:3),diconspo(k)%howmanythey(1)*3,&
		    mpi_double_precision,diconspo(k)%procid,diconspo(k)%procid,mpi_comm_world,status,ierror)
		    end if
		    end do

			      call mpi_barrier(mpi_comm_world,ierror)
			      k=0
			      jj1=0
			      do j=0,isize-1
				      if (j.ne.n)then
				      k=k+1
				      do i=1,diconspo(k)%howmanythey(1)
					vext(1,1)=diconspo(k)%facx(i,1); vext(1,2)=diconspo(k)%facx(i,2); vext(1,3)=diconspo(k)%facx(i,3)


				      do kj=1,dICONRPA(n)%howmanyi(1)
					vext(2,1)=dICONRPA(n)%facx(kj,1); vext(2,2)=dICONRPA(n)%facx(kj,2);vext(2,3)=dICONRPA(n)%facx(kj,3)

		dist=distance3(n,vext)

                        if(per_rot.eq.0)then
						    if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
						    ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
						    if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
						    jj1=jj1+1
	      ! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall).or.((abs(dist-zper)).lt.tolsmall))then

						    diconspo(k)%whichthey(i,1)=dICONRPA(n)%whichi(kj,3);go to 401
						    end if
						    end if
						    if (((abs(vext(2,2)-yper).lt.tolsmall).or.(abs((abs(vext(2,2)-yper))-yper).lt.tolsmall)).and.&
						    ((abs(vext(1,2)-yper).lt.tolsmall).or.(abs((abs(vext(1,2)-yper))-yper).lt.tolsmall)))then
						    if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall).and.(abs(vext(2,3)-vext(1,3)).lt.tolsmall))then
						    jj1=jj1+1
	      ! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall).or.((abs(dist-zper)).lt.tolsmall))then

						    diconspo(k)%whichthey(i,1)=dICONRPA(n)%whichi(kj,3);go to 401
						    end if
						    end if
						    
						    
						    if (((abs(vext(2,3)-zper).lt.tolsmall).or.(abs((abs(vext(2,3)-zper))-zper).lt.tolsmall)).and.&
						    ((abs(vext(1,3)-zper).lt.tolsmall).or.(abs((abs(vext(1,3)-zper))-zper).lt.tolsmall)))then
						    if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall).and.(abs(vext(2,1)-vext(1,1)).lt.tolsmall))then
						    jj1=jj1+1
	      ! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall).or.((abs(dist-zper)).lt.tolsmall))then

						  diconspo(k)%whichthey(i,1)=dICONRPA(n)%whichi(kj,3);go to 401
						  
						    end if
						    end if
                        else
                        
                            code_per1=ibound_icode(ielem_ibounds(dICONRPA(n)%whichi(kj,2),dICONRPA(n)%whichi(kj,1)))
                            vext(2,:)=rotate_per(vext(2,:),code_per1,angle_per)
                            if ((abs(vext(1,1)-vext(2,1)).lt.tol_per).and.(abs(vext(1,2)-vext(2,2)).lt.tol_per).and.&
                                      (abs(vext(1,3)-vext(2,3)).lt.tol_per)) then
!                                      write(9000+n,'(9es24.15)'),vext(1,1:3),vext(2,1:3),abs(vext(1,1:3)-vext(2,1:3))
                                jj1=jj1+1
                                diconspo(k)%whichthey(i,1)=dICONRPA(n)%whichi(kj,3);go to 401
                            end if
                        end if

					  end do
					  401 continue
					  end do
				      end if
			      end do
! 			
			!-------------------for debugging only -----------------------------------------!
		
		      !-------------------for debugging only -----------------------------------------!
			k=0
		!-------------------for debugging only -----------------------------------------!
! 		
		!-------------------for debugging only -----------------------------------------!
		call mpi_barrier(mpi_comm_world,ierror)
			icpuid=n
		k=0
do i=0,isize-1
if (i.ne.n)then
  k=k+1
  call mpi_sendrecv(diconspo(k)%whichthey(1:diconspo(k)%howmanythey(1),1:1),&
diconspo(k)%howmanythey(1),mpi_integer,i,icpuid,&
  dICONRPf(k)%whichi(1:dICONRPf(k)%howmanyi(1),1:1),dICONRPf(k)%howmanyi(1),&
mpi_integer,diconspo(k)%procid,diconspo(k)%procid,mpi_comm_world,status,ierror)
end if
end do

k=0
! 		
		call mpi_barrier(mpi_comm_world,ierror)
		do i=0,isize-1
			if (i.ne.n)then
			k=k+1
			do kj=1,dICONRPf(k)%howmanyi(1)
			if (dICONRPf(k)%whichi(kj,1).gt.0)then
! 			
 			ielem_ineighg(dICONRPA(n)%whichi(kj,2),dICONRPA(n)%whichi(kj,1))=dICONRPf(k)%whichi(kj,1)
			end if
			end do
			end if
		end do
		kj=0
		do i=1,dICONRPA(n)%howmanyi(1)
			
			if (ielem_ineighg(dICONRPA(n)%whichi(i,2),dICONRPA(n)%whichi(i,1)).eq.0)then
			  if ((((ibound_icode(ielem_ibounds(dICONRPA(n)%whichi(i,2),dICONRPA(n)%whichi(i,1))))).eq.5).or.&
			  (((ibound_icode(ielem_ibounds(dICONRPA(n)%whichi(i,2),dICONRPA(n)%whichi(i,1))))).eq.50))then
			kj=kj+1
				  
			      
			end if
			end if
		end do

		if (kj.gt.0)then
		  print*,"error detected in boundary conditions matching rules",kj,dICONRPA(n)%howmanyi(1)
		end if
		call mpi_barrier(mpi_comm_world,ierror)

			deallocate(diconrpa)
			deallocate(diconrpf)
			deallocate(diconrpm)
			deallocate(diconspo)
		end if
		
		call mpi_barrier(mpi_comm_world,ierror)
		
! 		deallocate(icons)
		deallocate(diconr)
		
! 		

	
			

end subroutine cons



subroutine cons2d(n,diconr,dicons,iperiodicity,xmpielrank,isize,diconrpa,diconrpm,diconspo,xper,&
yper,zper,diconrpf,numneighbours,typesten)
!> @brief
!> this subroutine establishes the connectivity across different cpus in 2d
implicit none
integer,intent(in)::iperiodicity,n,isize,numneighbours,typesten
real::small,dist
real,intent(in)::xper,yper,zper
integer::i,k,jjj,kj,j,l,im,io,ir,kmaxe,icpuid,kkj
integer,dimension(6)::jx
real,dimension(3)::dumface1,dumface2
integer::faceper1,faceper2,p,kkj2,kkj3,kkj4,kkj5,ixxff,iconsi,facex
real,dimension(1:dimensiona)::cords
integer,allocatable,dimension(:),intent(in)::xmpielrank
type(aconnx),allocatable,dimension(:),intent(inout)::diconr,diconrpa,diconrpm,diconrpf
type(aconnx),allocatable,dimension(:),intent(inout)::dicons,diconspo
real,dimension(1:8,1:dimensiona)::vext

kmaxe=xmpielrank(n)

small=tolsmall


           
kj=0; kkj=0
	do i=1,kmaxe
	if (ielem_interior(i).eq.1)then
		do k=1,ielem_ifca(i)
		
		  if (((ielem_ineighg(k,i).eq.0)))then
		      if(ielem_ibounds(k,i).eq.0)then
						kj=kj+1
		      else
		      if (((ibound_icode(ielem_ibounds(k,i)))).eq.5)then
						  kj=kj+1	
						
		      end if
		      end if
		end if
		
		

		end do
	end if
	end do


allocate(diconr(n:n))
allocate(dICONs(1:isize-1))
allocate(diconr(n)%howmanyi(1))


do i=1,isize-1
allocate(dICONs(i)%howmanythey(1))
dICONs(i)%howmanythey(1)=0
end do
diconr(n)%procid=n
diconr(n)%howmanyi(1)=kj
k=0
do i=0,isize-1
	 if (i.ne.n)then
	 k=k+1
	 dICONs(k)%procid=i
	 end if
end do
call mpi_barrier(mpi_comm_world,ierror)
icpuid=n
k=0

do i=0,isize-1
	if (i.ne.n)then
			k=k+1
			call mpi_sendrecv(diconr(n)%howmanyi(1),1,mpi_integer,i,icpuid,&
			dICONs(k)%howmanythey(1),1,mpi_integer,dICONs(k)%procid,dICONs(k)%procid,mpi_comm_world,status,ierror)
	end if
end do

		allocate(diconr(n)%whichi(diconr(n)%howmanyi(1),4))
		diconr(n)%whichi(:,:)=0
		do i=1,isize-1
		allocate(dICONs(i)%whichthey(dICONs(i)%howmanythey(1),4))
		allocate(dICONs(i)%ret(dICONs(i)%howmanythey(1)))
		allocate(dICONs(i)%retm(diconr(n)%howmanyi(1)))
		dICONs(i)%whichthey(:,:)=0
		dICONs(i)%ret(:)=0
		dICONs(i)%retm(:)=0
		end do		
! 	
		kj=0
		do i=1,kmaxe
		if (ielem_interior(i).eq.1)then
		io=ielem_ifca(i)
		do k=1,io
		 if (((ielem_ineighg(k,i).eq.0)))then
		    if ((ielem_ibounds(k,i).eq.0))then
						kj=kj+1
				diconr(n)%whichi(kj,1)=i
				diconr(n)%whichi(kj,2)=k
		    
						
		      else
		      
		      if (((ibound_icode(ielem_ibounds(k,i)))).eq.5)then
						  kj=kj+1
				diconr(n)%whichi(kj,1)=i
				diconr(n)%whichi(kj,2)=k
		      end if
		      end if
		end if
		end do
		end if
		end do

		 



		do kj=1,diconr(n)%howmanyi(1)
			      
! 			     
			      diconr(n)%whichi(kj,3:4)=ielem_nodes_faces(diconr(n)%whichi(kj,2),1:2,diconr(n)%whichi(kj,1))


			
		end do
		call mpi_barrier(mpi_comm_world,ierror)
		icpuid=n
		k=0
		
		
do i=0,isize-1
if (i.ne.n)then
    k=k+1
    call mpi_sendrecv(diconr(n)%whichi(1:diconr(n)%howmanyi(1),1:4),diconr(n)%howmanyi(1)*4,mpi_integer,i,icpuid,&
    dICONs(k)%whichthey(1:dICONs(k)%howmanythey(1),1:4),dICONs(k)%howmanythey(1)*4,&
mpi_integer,dICONs(k)%procid,dICONs(k)%procid,mpi_comm_world,status,ierror)
end if
end do





do kj=1,diconr(n)%howmanyi(1)
do i=1,isize-1
      do k=1,dICONs(i)%howmanythey(1)
 	      if (dinoder(dICONs(i)%whichthey(k,3))%itor.gt.0)then
	      
	      if ((((diconr(n)%whichi(kj,3)).eq.(dICONs(i)%whichthey(k,3))).or.((diconr(n)%whichi(kj,3))&
.eq.(dICONs(i)%whichthey(k,4)))).and.(((diconr(n)%whichi(kj,4)).eq.(dICONs(i)%whichthey(k,3))).or.((diconr(n)%whichi(kj,4))&
.eq.(dICONs(i)%whichthey(k,4)))))then
	      dICONs(i)%ret(k)=ielem_ihexgl(diconr(n)%whichi(kj,1))
	      
	      go to 801
	      end if
	      end if
      end do
end do
801 continue
end do

		call mpi_barrier(mpi_comm_world,ierror)

		k=0
		icpuid=n
		k=0
do i=0,isize-1
if (i.ne.n)then
  k=k+1
  call mpi_sendrecv(dICONs(k)%ret(1:dICONs(k)%howmanythey(1)),dICONs(k)%howmanythey(1),mpi_integer,dICONs(k)%procid,icpuid,&
  dICONs(k)%retm(1:diconr(n)%howmanyi(1)),diconr(n)%howmanyi(1),mpi_integer,i,i,mpi_comm_world,status,ierror)
end if
end do

			
			do i=1,diconr(n)%howmanyi(1)
				k=0
				do j=0,isize-1
					if (j.ne.n)then
					k=k+1
					if (dICONs(k)%retm(i).gt.0)then
					ielem_ineighg(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))=dICONs(k)%retm(i)
					end if
					end if
				end do
			end do
		  
		  kj=0
		do i=1,diconr(n)%howmanyi(1)


		       if ((ielem_ineighg(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1)).eq.0))then
		      if (iperiodicity.eq.1)then
		      
! 		     
		      if (ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1)).ne.0)then
		      if (((ibound_icode(ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))))).eq.5)then
						  kj=kj+1
! 			
		      end if
		      end if
		      else

					    kj=kj+1
		      end if
		      end if
		end do

	     
		
!-------------------for debugging only -----------------------------------------!



!-------------------for debugging only -----------------------------------------!
		if (iperiodicity.eq.1)then
! 		
		
		allocate(diconrpm(n:n))
		allocate(diconrpm(n)%howmanyi(1))
		diconrpm(n)%howmanyi(1)=kj
		allocate(diconrpm(n)%whichi(kj,3))
		diconrpm(n)%whichi(:,:)=0
		kj=0
! 		
		


		 kj=0
		do i=1,diconr(n)%howmanyi(1)


		       if ((ielem_ineighg(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1)).eq.0))then
		      if (((ibound_icode(ielem_ibounds(diconr(n)%whichi(i,2),diconr(n)%whichi(i,1))))).eq.5)then
						  kj=kj+1
			diconrpm(n)%whichi(kj,1)=diconr(n)%whichi(i,1)
			diconrpm(n)%whichi(kj,3)=ielem_ihexgl(diconr(n)%whichi(i,1))
			diconrpm(n)%whichi(kj,2)=(diconr(n)%whichi(i,2))
		      end if
		      end if
		end do
! 		
		
		jjj=0
!-------------------for debugging only -----------------------------------------!
		allocate(dICONRPA(n:n))
		allocate(dICONRPf(1:isize-1))
! 		 
! 		dICONRPA(n:n)=0
! 		dICONRPf(1:isize-1)=0
		do i=1,isize-1
		allocate(dICONRPf(i)%howmanyi(1))
		dICONRPf(i)%howmanyi(1)=0
		dICONRPf(i)%howmanyi(1)=diconrpm(n)%howmanyi(1)-jjj
		allocate(dICONRPf(i)%whichi(dICONRPf(i)%howmanyi(1),1))
		dICONRPf(i)%whichi(:,:)=0
		end do
		allocate(dICONRPA(n)%howmanyi(1))
		dICONRPA(n)%howmanyi(1)=0
		dICONRPA(n)%howmanyi(1)=diconrpm(n)%howmanyi(1)-jjj
		allocate(dICONRPA(n)%whichi(diconrpm(n)%howmanyi(1)-jjj,3))
		allocate(dICONRPA(n)%facx(diconrpm(n)%howmanyi(1)-jjj,2))
		dICONRPA(n)%whichi(:,:)=0
		dICONRPA(n)%facx(:,:)=0.d0
		kj=0
		do i=1,diconrpm(n)%howmanyi(1)
! 			
			 if ((ielem_ineighg(diconrpm(n)%whichi(i,2),diconrpm(n)%whichi(i,1)).eq.0))then
		      if (((ibound_icode(ielem_ibounds(diconrpm(n)%whichi(i,2),diconrpm(n)%whichi(i,1))))).eq.5)then
! 						  kj=kj+1
! 			
			kj=kj+1
			dICONRPA(n)%whichi(kj,1)=diconrpm(n)%whichi(i,1)
			dICONRPA(n)%whichi(kj,2)=diconrpm(n)%whichi(i,2)
			dICONRPA(n)%whichi(kj,3)=diconrpm(n)%whichi(i,3)
! 			  if (ielem_types_faces(dICONRPA(n)%whichi(kj,2),dICONRPA(n)%whichi(kj,1)).eq.5)then

			  ixxff=2
! 			  else
! 			  ixxff=3
! 			  end if
			    iconsi=dICONRPA(n)%whichi(kj,1); facex=dICONRPA(n)%whichi(kj,2)

			  call compute_centre2df(n,iconsi,facex,ixxff,cords)
			

			dICONRPA(n)%facx(kj,1)=cords(1)
			dICONRPA(n)%facx(kj,2)=cords(2)
			
			end if
			end if
! 			end if
		end do
		
		dumface1(:)=0.0
		dumface2(:)=0.0
		faceper1=0
		faceper2=0
		icpuid=n
		call mpi_barrier(mpi_comm_world,ierror)
! 		
!-------------------for debugging only -----------------------------------------!
		allocate(diconspo(1:isize-1))
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
		do i=1,isize-1
		allocate(diconspo(i)%howmanythey(1))
		diconspo(i)%howmanythey(1)=0
		end do
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
		k=0
		do i=0,isize-1
		if (i.ne.n)then
		k=k+1
		diconspo(k)%procid=i
	 	end if
		end do
!-------------------for debugging only -----------------------------------------!
!-----------------------------!
		call mpi_barrier(mpi_comm_world,ierror)
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
			icpuid=n
		k=0
		do i=0,isize-1
			if (i.ne.n)then
					k=k+1
    call mpi_sendrecv(dICONRPA(n)%howmanyi(1),1,mpi_integer,i,icpuid,&
    diconspo(k)%howmanythey(1),1,mpi_integer,diconspo(k)%procid,diconspo(k)%procid,mpi_comm_world,status,ierror)
			end if
		end do


		call mpi_barrier(mpi_comm_world,ierror)
		do i=1,isize-1
		allocate(diconspo(i)%whichthey(diconspo(i)%howmanythey(1),1))
		allocate(diconspo(i)%facx(diconspo(i)%howmanythey(1),3))
		diconspo(i)%whichthey(:,:)=0
		diconspo(i)%facx(:,:)=0.d0
		end do
		call mpi_barrier(mpi_comm_world,ierror)
			icpuid=n
		k=0
do i=0,isize-1
if (i.ne.n)then
  k=k+1
  call mpi_sendrecv(dICONRPA(n)%facx(1:dICONRPA(n)%howmanyi(1),1:2),&
dICONRPA(n)%howmanyi(1)*2,mpi_double_precision,i,icpuid,&
  diconspo(k)%facx(1:diconspo(k)%howmanythey(1),1:2),diconspo(k)%howmanythey(1)*2,&
mpi_double_precision,diconspo(k)%procid,diconspo(k)%procid,mpi_comm_world,status,ierror)
end if
end do
!-------------------for debugging only -----------------------------------------!

		call mpi_barrier(mpi_comm_world,ierror)
		k=0
		do j=0,isize-1
			if (j.ne.n)then
			k=k+1
			do i=1,diconspo(k)%howmanythey(1)
			  vext(1,1)=diconspo(k)%facx(i,1); vext(1,2)=diconspo(k)%facx(i,2)


			do kj=1,dICONRPA(n)%howmanyi(1)
			  vext(2,1)=dICONRPA(n)%facx(kj,1); vext(2,2)=dICONRPA(n)%facx(kj,2)

  dist=distance2(n,vext)

 if (((abs(vext(2,1)-xper).lt.tolsmall).or.(abs((abs(vext(2,1)-xper))-xper).lt.tolsmall)).and.&
				      ((abs(vext(1,1)-xper).lt.tolsmall).or.(abs((abs(vext(1,1)-xper))-xper).lt.tolsmall)))then
				      if ((abs(vext(2,2)-vext(1,2)).lt.tolsmall))then


! 				      if (((abs(dist-xper)).lt.tolsmall).or.((abs(dist-yper)).lt.tolsmall))then
						
				     diconspo(k)%whichthey(i,1)=dICONRPA(n)%whichi(kj,3);go to 901
				      end if
				      end if
				      if (((abs(vext(2,2)-yper).lt.tolsmall).or.(abs((abs(vext(2,2)-yper))-yper).lt.tolsmall)).and.&
				      ((abs(vext(1,2)-yper).lt.tolsmall).or.(abs((abs(vext(1,2)-yper))-yper).lt.tolsmall)))then
				      if ((abs(vext(2,1)-vext(1,1)).lt.tolsmall))then



						
				      diconspo(k)%whichthey(i,1)=dICONRPA(n)%whichi(kj,3);go to 901
				      end if
				      end if




			    end do
			    901 continue
			    end do
			end if
		end do
				
			!-------------------for debugging only -----------------------------------------!
		
		      !-------------------for debugging only -----------------------------------------!
		k=0
		!-------------------for debugging only -----------------------------------------!
! 		
		!-------------------for debugging only -----------------------------------------!
		call mpi_barrier(mpi_comm_world,ierror)
			icpuid=n
		k=0
do i=0,isize-1
if (i.ne.n)then
  k=k+1
  call mpi_sendrecv(diconspo(k)%whichthey(1:diconspo(k)%howmanythey(1),1:1),&
diconspo(k)%howmanythey(1),mpi_integer,i,icpuid,&
  dICONRPf(k)%whichi(1:dICONRPf(k)%howmanyi(1),1:1),dICONRPf(k)%howmanyi(1),&
mpi_integer,diconspo(k)%procid,diconspo(k)%procid,mpi_comm_world,status,ierror)
end if
end do
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
k=0
! 		
		call mpi_barrier(mpi_comm_world,ierror)
		do i=0,isize-1
			if (i.ne.n)then
			k=k+1
			do kj=1,dICONRPf(k)%howmanyi(1)
			if (dICONRPf(k)%whichi(kj,1).gt.0)then
! 			
 			ielem_ineighg(dICONRPA(n)%whichi(kj,2),dICONRPA(n)%whichi(kj,1))=dICONRPf(k)%whichi(kj,1)
			end if
			end do
			end if
		end do
		kj=0
		do i=1,dICONRPA(n)%howmanyi(1)
			
			if (ielem_ineighg(dICONRPA(n)%whichi(i,2),dICONRPA(n)%whichi(i,1)).eq.0)then
			  if (((ibound_icode(ielem_ibounds(dICONRPA(n)%whichi(i,2),dICONRPA(n)%whichi(i,1))))).eq.5)then
			kj=kj+1
				  
			      
			end if
			end if
		end do
! 		
		if (kj.gt.0)then
		  print*,"error detected in boundary conditions matching rules"
		end if
		call mpi_barrier(mpi_comm_world,ierror)
! 		
		deallocate(diconrpa)
		deallocate(diconrpf)
		deallocate(diconrpm)
		deallocate(diconspo)
		end if
		
		call mpi_barrier(mpi_comm_world,ierror)
		
		deallocate(dicons)
		deallocate(diconr)

! 		

	
			

end subroutine cons2d

!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine detstenx(n)
!> @brief
!> this subroutine builds the large stencils
implicit none
integer,intent(in)::n
integer::i,k,jjj,kj,j,l,im,io,ir,kmaxe,icpuid,itrue,ic,print_i,kx,kxk,kk,candid,candid2
integer,dimension(1)::stcon,stconc,stcons,stcong,isosa,iistart,ifsat,ix
integer::pik,jk,ibn,igd1,icomp_set,i2comp_set,extf2,test6
real::distf,x_ste,y_ste,z_ste,max_sten,min_sten,max_sten2,min_sten2,rcomp_set,testdist,testdiv
real,dimension(1,1:iselem)::ilocalallelg3,ilocalallelgd


kmaxe=xmpielrank(n)






icpuid=n




print_i=kmaxe/20

!$omp do
do k=1,kmaxe

	ilocalallelg(n,k,1,1)=ielem_ihexgl(k)

	stcon(1)=k
	stconc(1)=k
	stcons(1)=0
	stcong(1)=ielem_ihexgl(k)
	isosa(1)=1
	iistart(1)=1
	ifsat(1)=0
	ix(1)=0

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!


call allsx(stcon,stconc,stcons,stcong,isosa,ifsat,iistart,ix)

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!

        
       
        
! 	if (iweno.eq.1)then
	
      if (icompact.ge.1)then      !if one
	
                    if (dimensiona.eq.3)then        !if dimensiona
                                        ibn=ielem_ifca(k)
                                        ilocalallelg3(1,1:iselem)=zero
                                        
                                        ilocalallelg3(1,1:iselemt(n)-ielem_ifca(k))=ilocalallelg(n,k,1,1:iselemt(n)-ibn)
                                        ilocalallelgd(:,:)=zero
                                                        do j=2,iselemt(n)-ielem_ifca(k)
							    candid2=ilocalallelg(n,k,1,j)
							    if (cand(xmpie(candid2)).gt.0)then
                                                            
                                                            x_ste=xand2r(cand(xmpie(candid2)),xmpil(candid2),1)
                                                            y_ste=xand2r(cand(xmpie(candid2)),xmpil(candid2),2)
                                                            z_ste=xand2r(cand(xmpie(candid2)),xmpil(candid2),3)
                                                            if (iperiodicity.eq.1)then      !if three
                                                                    
                                                                        if(abs(x_ste-ielem_xxc(k)).gt.xper*oo2)then
                                                                        x_ste=x_ste+(xper*sign(1.0,ielem_xxc(k)-xper*oo2))
                                                                        end if
                                                                        if(abs(y_ste-ielem_yyc(k)).gt.yper*oo2)then
                                                                        y_ste=y_ste+(yper*sign(1.0,ielem_yyc(k)-yper*oo2))
                                                                        end if
                                                                        if(abs(z_ste-ielem_zzc(k)).gt.zper*oo2)then
                                                                        z_ste=z_ste+(zper*sign(1.0,ielem_zzc(k)-zper*oo2))
                                                                        end if
                                                            
                                                            end if          !if three
                                                        ilocalallelgd(1,j)=sqrt(((x_ste-ielem_xxc(k))**2)+&
                                                ((y_ste-ielem_yyc(k))**2)+((z_ste-ielem_zzc(k))**2)+(tolsmall*j))
							    else
							    write(99+n,*)"killed here due to stencils requiring more cpus"
							    stop
							    
							    
							    end if
                                                        end do
                                        
                                        max_sten2=tolsmall
                                        min_sten2=tolbig
                                        
                                        do kk=1,ilx
                                        
                                                min_sten2=min(min_sten2,ilocalallelgd(1,kk))
                                        
                                                max_sten2=max(max_sten2,ilocalallelgd(1,kk))
                                        end do
	
	
      
                                    igd1=1
                                    testdist=zero
                                    max_sten=tolsmall
                                    min_sten=tolbig
                                                                            do j=1,ielem_ifca(k)
                                                    if (ielem_ineighg(j,k).gt.0)then
                                                    igd1=igd1+1
                                                    testdist=testdist+ilocalallelgd(1,igd1)
                                                        if (ilocalallelgd(1,igd1).ge.max_sten)then
                                                        max_sten=ilocalallelgd(1,igd1)
                                                        
                                                        end if
                                                        if (ilocalallelgd(1,igd1).le.min_sten)then
                                                        min_sten=ilocalallelgd(1,igd1)
                                                        end if
                                                    end if
                                                end do
                                      testdiv=igd1-1
	testdist=(testdist)/(testdiv*2)
!
                                            if (icompact.eq.2)then   !if icompact 
                                ! 	       ielem_inumneighbours(k)=max(ielem_inumneighbours(k))
                                        imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
                                        
                                        kxk=0
                                        
                                                            do j=ielem_ifca(k)+2,iselemt(n)-ibn
                                                            do jk=iselemt(n)-ibn,j+1,-1
                                                                        if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                                        pik= ilocalallelg3(1,jk-1)
                                                                        distf=ilocalallelgd(1,jk-1)
                                                                            ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                                        ilocalallelg3(1,jk) = pik
                                                                            ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                                        ilocalallelgd(1,jk)=distf
                                                                        endif
                                                            enddo ! j
                                                            enddo  ! i
                                                                            
                                                                            
                                                                            
                                            ilocalallelg(n,k,1,ielem_ifca(k)+2:iselemt(n))=ilocalallelg3(1,ielem_ifca(k)+2:iselemt(n))
                                            
                                            
                                            
                                            
                                            end if           !if three
	      
	      
	      
                                                if (icompact.eq.3)then   !if icompact 3
                                        kxk=0
                                        icomp_set=igd1+1
                                !  	  
                                                                    if ((testdist.le.(1.4*ielem_minedge(k))).and.(max_sten/min_sten.le.(1.4)))then
                                                                ielem_vortex(1,k)=100.0
                                                                
                                                                                                
                                                                                                
                                                                                                            do j=icomp_set,iselemt(n)-ibn
                                                                                                            do jk=iselemt(n)-icomp_set,j+1,-1
                                                                                                        if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                                                                        pik= ilocalallelg3(1,jk-1)
                                                                                                        distf=ilocalallelgd(1,jk-1)
                                                                                                            ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                                                                        ilocalallelg3(1,jk) = pik
                                                                                                            ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                                                                        ilocalallelgd(1,jk)=distf
                                                                                                        endif
                                                                                                    enddo ! j
                                                                                                    enddo  ! i
                                                                                                    
                                                                                                    
                                                                        max_sten2=((testdist)*(iorder+1))  
                                                                        
                                                                            do kk=icomp_set,iselemt(n)-ibn
                                                                            if  ((ilocalallelgd(1,kk)).le.max_sten2)then
                                                                                
                                                                                kxk=kxk+1
                                                                            
                                                                            end if
                                                                            end do
                                                                            
                                                                    
                                                                    end if     !if three
                                            test6=ielem_inumneighbours(k)*1.2
                                            
                                            ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,test6))
                                        imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
                                        ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))
                                        
                                        
                                        
                                        
                                        end if      !if three
                                        if ((icompact.eq.1).and.((max_sten/min_sten).le.1.2))then      !if icompact 1
                        
                                            extf2=extf      
                                            icomp_set=(ilx*extf2/2)+1
                                            
                                            kxk=0
                                        
                                                                do j=icomp_set,iselemt(n)-ibn
                                                                do jk=iselemt(n)-icomp_set-ibn,j+1,-1
                                                                        if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                                        pik= ilocalallelg3(1,jk-1)
                                                                        distf=ilocalallelgd(1,jk-1)
                                                                            ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                                        ilocalallelg3(1,jk) = pik
                                                                            ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                                        ilocalallelgd(1,jk)=distf
                                                                        endif
                                                        enddo ! j
                                                        enddo  ! i
                                            
                                            
                                                            do kk=icomp_set,iselemt(n)-ibn
                                                                if  ((ilocalallelgd(1,kk)).le.max_sten2)then
                                            !                          
                                                                    kxk=kxk+1
                                                                
                                                                end if
                                                                end do
                                !  	

                                        ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,ilx*extf2))
                                        imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
                                        ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))

                                        
                                !  	 end if
                                        
                                            
                                        end if
	
	   
	   
	   
	  else
     
      
         
	 
	   ibn=ielem_ifca(k)
	ilocalallelg3(1,1:iselem)=zero
	
	ilocalallelg3(1,1:iselemt(n)-ielem_ifca(k))=ilocalallelg(n,k,1,1:iselemt(n)-ibn)
! 	
	ilocalallelgd(:,:)=zero
	do j=2,iselemt(n)-ielem_ifca(k)
	
	    candid2=ilocalallelg(n,k,1,j)
		if (cand(xmpie(candid2)).gt.0)then
                 
                 x_ste=xand2r(cand(xmpie(candid2)),xmpil(candid2),1)
                 y_ste=xand2r(cand(xmpie(candid2)),xmpil(candid2),2)
                                                            
	
	
	
	   
	    if (iperiodicity.eq.1)then
		     
			if(abs(x_ste-ielem_xxc(k)).gt.xper*oo2)then
			x_ste=x_ste+(xper*sign(1.0,ielem_xxc(k)-xper*oo2))
			end if
			if(abs(y_ste-ielem_yyc(k)).gt.yper*oo2)then
			y_ste=y_ste+(yper*sign(1.0,ielem_yyc(k)-yper*oo2))
			end if
	    
	    
	    end if
! 	    if (icompact.eq.1)then
	  ilocalallelgd(1,j)=(sqrt(((x_ste-ielem_xxc(k))**2)+&
                   ((y_ste-ielem_yyc(k))**2)))+(tolsmall)
!             else
!             ilocalallelgd(1,j)=max(((sqrt(((x_ste-ielem_xxc(k))**2)+&
!                    ((y_ste-ielem_yyc(k))**2)+(tolsmall)))/(sqrt((x_ste-ielem_xxc(k))**2)+(tolsmall))),((sqrt(((x_ste-ielem_xxc(k))**2)+&
!                    ((y_ste-ielem_yyc(k))**2)+(tolsmall)))/(sqrt((y_ste-ielem_yyc(k))**2)+(tolsmall))))
!             
!             end if
	  else
	  write(99+n,*)"killed here due to stencils requiring more cpus"
							    stop
	  
	  
	  end if
	end do
	
	
	max_sten2=tolsmall
	min_sten2=tolbig
	
	
        
        
          do kk=1,ilx
                max_sten2=max(max_sten2,ilocalallelgd(1,kk))
                min_sten2=min(min_sten2,ilocalallelgd(1,kk))
          end do
!           
            
      
	igd1=1
	testdist=zero
	max_sten=tolsmall
	min_sten=tolbig
	do j=1,ielem_ifca(k)
	      if (ielem_ineighg(j,k).gt.0)then
	      igd1=igd1+1
	      testdist=testdist+ilocalallelgd(1,igd1)
		if (ilocalallelgd(1,igd1).ge.max_sten)then
		  max_sten=ilocalallelgd(1,igd1)
		  
		end if
		if (ilocalallelgd(1,igd1).le.min_sten)then
		  min_sten=ilocalallelgd(1,igd1)
		end if
	      end if
	end do
	testdiv=igd1-1
	testdist=(testdist)/(testdiv*2)
! 	    min_sten=2.0*ielem_minedge(k)
! 	        ielem_vortex(1,k)=max_sten/min_sten
	      if (icompact.eq.2)then
! 	       ielem_inumneighbours(k)=max(ielem_inumneighbours(k))
 	  imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
 	  
 	   kxk=0
	
                    do j=ielem_ifca(k)+2,iselemt(n)-ibn
                    do jk=iselemt(n)-ibn,j+1,-1
                if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                pik= ilocalallelg3(1,jk-1)
                distf=ilocalallelgd(1,jk-1)
                    ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                ilocalallelg3(1,jk) = pik
                    ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                ilocalallelgd(1,jk)=distf
                endif
            enddo ! j
            enddo  ! i
 	  
 	  
	   ilocalallelg(n,k,1,ielem_ifca(k)+2:iselemt(n))=ilocalallelg3(1,ielem_ifca(k)+2:iselemt(n))
	      
	      
	      
	      
	      end if
	      
	        if (icompact.eq.3)then
	   kxk=0
	   icomp_set=igd1+1

            if ((testdist.le.(1.4*ielem_minedge(k))).and.(max_sten/min_sten.le.(1.4)))then
	   ielem_vortex(1,k)=100.0
	   
                                        
                                        
                                                    do j=icomp_set,iselemt(n)-ibn
                                                    do jk=iselemt(n)-icomp_set,j+1,-1
                                                if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                pik= ilocalallelg3(1,jk-1)
                                                distf=ilocalallelgd(1,jk-1)
                                                    ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                ilocalallelg3(1,jk) = pik
                                                    ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                ilocalallelgd(1,jk)=distf
                                                endif
                                            enddo ! j
                                            enddo  ! i
                                            
                                            
                max_sten2=((testdist)*(iorder+1))  
                
	     do kk=icomp_set,iselemt(n)-ibn
                    if  ((ilocalallelgd(1,kk)).le.max_sten2)then
                         
                        kxk=kxk+1
                       
                    end if
                    end do
                    
 	    
	    end if
	    test6=ielem_inumneighbours(k)*2.6
	    
	    ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,test6))
 	  imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
	   ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))
	   
	   
	   
	   
	   end if
	      
	      
	     if ((icompact.eq.1).and.((max_sten/min_sten).le.1.2))then   

            extf2=extf          !extension
	    icomp_set=(ilx*extf2/2)+3      
	   
	    kxk=0
	
                    do j=icomp_set,iselemt(n)-ibn
                    do jk=iselemt(n)-icomp_set-ibn,j+1,-1
                if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                pik= ilocalallelg3(1,jk-1)
                distf=ilocalallelgd(1,jk-1)
                    ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                ilocalallelg3(1,jk) = pik
                    ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                ilocalallelgd(1,jk)=distf
                endif
            enddo ! j
            enddo  ! i
            
            
                   do kk=icomp_set,iselemt(n)-ibn
                    if  ((ilocalallelgd(1,kk)).le.max_sten2)then
!                          
                        kxk=kxk+1
                       
                    end if
                    end do


         ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,ilx*extf2))
 	  imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
	   ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))

 	 
!  	 end if
 	 
	    
	end if
	end if
        
       
    end if


	
end do
!$omp end do
    
 
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!

end subroutine detstenx





! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine allsx(stcon,stconc,stcons,stcong,isosa,ifsat,iistart,ix)
!> @brief
!> this recursive subroutine builds the large stencils until the prescribed number of elements is reached
implicit none
integer::i,k,jjj,kj,j,l,im,io,ir,kmaxe,icpuid,itrue,jloop,candid
integer,dimension(1),intent(inout)::stcon,stconc,stcons,stcong,isosa,iistart,ifsat,ix
kmaxe=xmpielrank(n)
      if (dimensiona.eq.3)then
			  jloop=6
		else
			  jloop=4
		end if
		icpuid=n
		iistart(1)=iistart(1)+1
		if (xmpie(stcong(1)).eq.n)then		!global array checking if element belongs to my cpu
			l=xmpil(stcong(1))		!if yes then set l=to the number of this cpu
			stconc(1)=l			!use this stconc(n)=l for referencing after
		
		do j=1,ielem_ifca(stconc(1))		!loop all the sides of this element
			if (ielem_ineighg(j,stconc(1)).gt.0)then  !if the neighbour is gt.0 then
			ix(1)=ielem_ineighg(j,stconc(1))		!set the ix(n) as the global index of this element
			call check(n,stcon,ix,ifsat)	!check if this element is already in the list
			if (ifsat(1).eq.1)then					!if not then include
			if (isosa(1).le.iselemt(n)-1)then
			isosa(1)=isosa(1)+1
			ilocalallelg(n,stcon(1),1,isosa(1))=ielem_ineighg(j,stconc(1))
			end if
			end if
			end if
! 			
		end do
		end if
		
		if (xmpie(stcong(1)).ne.n)then	!if this element belongs to another cpu then
		      if (cand(xmpie(stcong(1))).gt.0)then
		      do j=1,jloop			!loop all the sides of this particular element (i need to find the cpu that it belongs) therefore indexing to parent cpu
			      candid=cand2r(cand(xmpie(stcong(1))),xmpil(stcong(1)),j)
			
			      if(candid.gt.0) then
			      
			      ix(1)=candid
			      call check(n,stcon,ix,ifsat)
			      if (ifsat(1).eq.1)then
			      if (isosa(1).le.iselemt(n)-1)then
			      isosa(1)=isosa(1)+1
			      ilocalallelg(n,stcon(1),1,isosa(1))=candid
			      end if
			      end if
			      end if
		      end do
		      else
		
		      
		      
		      write(99+n,*)"killed due to stencils requiring more cpus"
		      stop
		      end if
		
! 		
		end if
		if (isosa(1).lt.iselemt(n))then
			!-------------------for debugging only -----------------------------------------!
! 			
			!-------------------for debugging only -----------------------------------------!
			stcong(1)=ilocalallelg(n,stcon(1),1,iistart(1))
			!-------------------for debugging only -----------------------------------------!
! 			
			!-------------------for debugging only -----------------------------------------!
			call allsx(stcon,stconc,stcons,stcong,isosa,ifsat,iistart,ix)
		end if
		if (isosa(1).eq.iselemt(n))then
			return
		end if

			
end subroutine allsx



subroutine detsten(n)
!> @brief
!> this subroutine builds the large stencils suitable for periodic boundaries
implicit none
integer,intent(in)::n
integer::i,k,jjj,kj,j,l,im,io,ir,kmaxe,icpuid,itrue,ic,print_i,kx,kxk,kk
integer,dimension(1)::stcon,stconc,stcons,stcong,isosa,iistart,ifsat,ix
integer::pik,jk,ibn,igd1,icomp_set,i2comp_set,extf2,test6
real::distf,x_ste,y_ste,z_ste,max_sten,min_sten,max_sten2,min_sten2,rcomp_set,testdist,testdiv
real,dimension(1,1:iselem)::ilocalallelg3,ilocalallelgd
kmaxe=xmpielrank(n)






icpuid=n


print_i=kmaxe/20

!$omp do
do k=1,kmaxe

	stcon(1)=k
	stconc(1)=k
	stcons(1)=0
	stcong(1)=ielem_ihexgl(k)
	isosa(1)=1
	iistart(1)=1
	ifsat(1)=0
	ix(1)=0

! 	ilocalallel(n,k,1,1)=ielem_ihex(k)
	ilocalallelg(n,k,1,1)=ielem_ihexgl(k)
	ilocalallelgper(n,k,1,1)=0
! 	ilocalalls(n,k,1,1)=ielem_ifca(k)

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
! 	call alls(ic,isize,ielem,stcon,stconc,stcons,stcong,iselemt,ilocalalls,&
! ilocalallelg,ilocalallelgper,xmpie,isosa,ifsat,iistart,ix,xmpielrank,pare,dose,pareel,doseel,pares,soseel,ifin,tfin,xmpil,glneigh,glneighper)
   call allsf(stcon,stconc,stcons,stcong,isosa,ifsat,iistart,ix)


        
        
        
! 	if (iweno.eq.1)then
	
      if (icompact.ge.1)then      !if one
	
                    if (dimensiona.eq.3)then        !if dimensiona
                                        ibn=ielem_ifca(k)
                                        ilocalallelg3(1,1:iselem)=zero
                                        
                                        ilocalallelg3(1,1:iselemt(n)-ielem_ifca(k))=ilocalallelg(n,k,1,1:iselemt(n)-ibn)
                                        ilocalallelgd(:,:)=zero
                                                        do j=2,iselemt(n)-ielem_ifca(k)
                                                            
                                                            x_ste=centerr(ilocalallelg(n,k,1,j),1);y_ste=centerr(ilocalallelg(n,k,1,j),2);z_ste=centerr(ilocalallelg(n,k,1,j),3)
                                                            if (iperiodicity.eq.1)then      !if three
                                                                    
                                                                        if(abs(x_ste-ielem_xxc(k)).gt.xper*oo2)then
                                                                        x_ste=x_ste+(xper*sign(1.0,ielem_xxc(k)-xper*oo2))
                                                                        end if
                                                                        if(abs(y_ste-ielem_yyc(k)).gt.yper*oo2)then
                                                                        y_ste=y_ste+(yper*sign(1.0,ielem_yyc(k)-yper*oo2))
                                                                        end if
                                                                        if(abs(z_ste-ielem_zzc(k)).gt.zper*oo2)then
                                                                        z_ste=z_ste+(zper*sign(1.0,ielem_zzc(k)-zper*oo2))
                                                                        end if
                                                            
                                                            end if          !if three
                                                        ilocalallelgd(1,j)=sqrt(((x_ste-ielem_xxc(k))**2)+&
                                                ((y_ste-ielem_yyc(k))**2)+((z_ste-ielem_zzc(k))**2)+(tolsmall*j))
                                                        end do
                                        
                                        max_sten2=tolsmall
                                        min_sten2=tolbig
                                        
                                        do kk=1,ilx
                                        
                                                min_sten2=min(min_sten2,ilocalallelgd(1,kk))
                                        
                                                max_sten2=max(max_sten2,ilocalallelgd(1,kk))
                                        end do
	
	
      
                                    igd1=1
                                    testdist=zero
                                    max_sten=tolsmall
                                    min_sten=tolbig
                                                                            do j=1,ielem_ifca(k)
                                                    if (ielem_ineighg(j,k).gt.0)then
                                                    igd1=igd1+1
                                                    testdist=testdist+ilocalallelgd(1,igd1)
                                                        if (ilocalallelgd(1,igd1).ge.max_sten)then
                                                        max_sten=ilocalallelgd(1,igd1)
                                                        
                                                        end if
                                                        if (ilocalallelgd(1,igd1).le.min_sten)then
                                                        min_sten=ilocalallelgd(1,igd1)
                                                        end if
                                                    end if
                                                end do
                                      testdiv=igd1-1
	testdist=(testdist)/(testdiv*2)
! 	    min_sten=2.0*ielem_minedge(k)
! 	        ielem_vortex(1,k)=max_sten/min_sten
                                            if (icompact.eq.2)then   !if icompact 
                                ! 	       ielem_inumneighbours(k)=max(ielem_inumneighbours(k))
                                        imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
                                        
                                        kxk=0
                                        
                                                            do j=ielem_ifca(k)+2,iselemt(n)-ibn
                                                            do jk=iselemt(n)-ibn,j+1,-1
                                                                        if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                                        pik= ilocalallelg3(1,jk-1)
                                                                        distf=ilocalallelgd(1,jk-1)
                                                                            ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                                        ilocalallelg3(1,jk) = pik
                                                                            ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                                        ilocalallelgd(1,jk)=distf
                                                                        endif
                                                            enddo ! j
                                                            enddo  ! i
                                                                            
                                                                            
                                                                            
                                            ilocalallelg(n,k,1,ielem_ifca(k)+2:iselemt(n))=ilocalallelg3(1,ielem_ifca(k)+2:iselemt(n))
                                            
                                            
                                            
                                            
                                            end if           !if three
	      
	      
	      
                                                if (icompact.eq.3)then   !if icompact 3
                                        kxk=0
                                        icomp_set=igd1+1
                                !  	   
                                                                    if ((testdist.le.(1.4*ielem_minedge(k))).and.(max_sten/min_sten.le.(1.4)))then
                                                                ielem_vortex(1,k)=100.0
                                                                
                                                                                                
                                                                                                
                                                                                                            do j=icomp_set,iselemt(n)-ibn
                                                                                                            do jk=iselemt(n)-icomp_set,j+1,-1
                                                                                                        if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                                                                        pik= ilocalallelg3(1,jk-1)
                                                                                                        distf=ilocalallelgd(1,jk-1)
                                                                                                            ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                                                                        ilocalallelg3(1,jk) = pik
                                                                                                            ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                                                                        ilocalallelgd(1,jk)=distf
                                                                                                        endif
                                                                                                    enddo ! j
                                                                                                    enddo  ! i
                                                                                                    
                                                                                                    
                                                                        max_sten2=((testdist)*(iorder+1))  
                                                                        
                                                                            do kk=icomp_set,iselemt(n)-ibn
                                                                            if  ((ilocalallelgd(1,kk)).le.max_sten2)then
                                                                                
                                                                                kxk=kxk+1
                                                                            
                                                                            end if
                                                                            end do
                                                                            
                                                                    
                                                                    end if     !if three
                                            test6=ielem_inumneighbours(k)*1.2
                                            
                                            ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,test6))
                                        imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
                                        ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))
                                        
                                        
                                        
                                        
                                        end if      !if three
                                        if ((icompact.eq.1).and.((max_sten/min_sten).le.1.2))then      !if icompact 1
                        
                                            extf2=extf      
                                            icomp_set=(ilx*extf2/2)+1
                                            
                                            kxk=0
                                        
                                                                do j=icomp_set,iselemt(n)-ibn
                                                                do jk=iselemt(n)-icomp_set-ibn,j+1,-1
                                                                        if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                                        pik= ilocalallelg3(1,jk-1)
                                                                        distf=ilocalallelgd(1,jk-1)
                                                                            ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                                        ilocalallelg3(1,jk) = pik
                                                                            ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                                        ilocalallelgd(1,jk)=distf
                                                                        endif
                                                        enddo ! j
                                                        enddo  ! i
                                            
                                            
                                                            do kk=icomp_set,iselemt(n)-ibn
                                                                if  ((ilocalallelgd(1,kk)).le.max_sten2)then
                                            !                          
                                                                    kxk=kxk+1
                                                                
                                                                end if
                                                                end do
                                !  	 

                                        ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,ilx*extf2))
                                        imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
                                        ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))

                                        
                                !  	 end if
                                        
                                            
                                        end if
	
	   
	   
	   
	  else
     
      
         
	 
	   ibn=ielem_ifca(k)
	ilocalallelg3(1,1:iselem)=zero
	
	ilocalallelg3(1,1:iselemt(n)-ielem_ifca(k))=ilocalallelg(n,k,1,1:iselemt(n)-ibn)
! 	
	ilocalallelgd(:,:)=zero
	do j=2,iselemt(n)-ielem_ifca(k)
	    x_ste=centerr(ilocalallelg(n,k,1,j),1);y_ste=centerr(ilocalallelg(n,k,1,j),2)
	    if (iperiodicity.eq.1)then
		     
			if(abs(x_ste-ielem_xxc(k)).gt.xper*oo2)then
			x_ste=x_ste+(xper*sign(1.0,ielem_xxc(k)-xper*oo2))
			end if
			if(abs(y_ste-ielem_yyc(k)).gt.yper*oo2)then
			y_ste=y_ste+(yper*sign(1.0,ielem_yyc(k)-yper*oo2))
			end if
	    
	    
	    end if
! 	    if (icompact.eq.1)then
	  ilocalallelgd(1,j)=(sqrt(((x_ste-ielem_xxc(k))**2)+&
                   ((y_ste-ielem_yyc(k))**2)))+(tolsmall)
!             else
!             ilocalallelgd(1,j)=max(((sqrt(((x_ste-ielem_xxc(k))**2)+&
!                    ((y_ste-ielem_yyc(k))**2)+(tolsmall)))/(sqrt((x_ste-ielem_xxc(k))**2)+(tolsmall))),((sqrt(((x_ste-ielem_xxc(k))**2)+&
!                    ((y_ste-ielem_yyc(k))**2)+(tolsmall)))/(sqrt((y_ste-ielem_yyc(k))**2)+(tolsmall))))
!             
!             end if
	end do
	
	
	max_sten2=tolsmall
	min_sten2=tolbig
	
	
        
        
          do kk=1,ilx
                max_sten2=max(max_sten2,ilocalallelgd(1,kk))
                min_sten2=min(min_sten2,ilocalallelgd(1,kk))
          end do
!            
            
            
      
	igd1=1
	testdist=zero
	max_sten=tolsmall
	min_sten=tolbig
	do j=1,ielem_ifca(k)
	      if (ielem_ineighg(j,k).gt.0)then
	      igd1=igd1+1
	      testdist=testdist+ilocalallelgd(1,igd1)
		if (ilocalallelgd(1,igd1).ge.max_sten)then
		  max_sten=ilocalallelgd(1,igd1)
		  
		end if
		if (ilocalallelgd(1,igd1).le.min_sten)then
		  min_sten=ilocalallelgd(1,igd1)
		end if
	      end if
	end do
	testdiv=igd1-1
	testdist=(testdist)/(testdiv*2)
! 	    min_sten=2.0*ielem_minedge(k)
! 	        ielem_vortex(1,k)=max_sten/min_sten
	      if (icompact.eq.2)then
! 	       ielem_inumneighbours(k)=max(ielem_inumneighbours(k))
 	  imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
 	  
 	   kxk=0
	
                    do j=ielem_ifca(k)+2,iselemt(n)-ibn
                    do jk=iselemt(n)-ibn,j+1,-1
                if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                pik= ilocalallelg3(1,jk-1)
                distf=ilocalallelgd(1,jk-1)
                    ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                ilocalallelg3(1,jk) = pik
                    ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                ilocalallelgd(1,jk)=distf
                endif
            enddo ! j
            enddo  ! i
 	  
 	  
	   ilocalallelg(n,k,1,ielem_ifca(k)+2:iselemt(n))=ilocalallelg3(1,ielem_ifca(k)+2:iselemt(n))
	      
	      
	      
	      
	      end if
	      
	        if (icompact.eq.3)then
	   kxk=0
	   icomp_set=igd1+1
!  	   
            if ((testdist.le.(1.4*ielem_minedge(k))).and.(max_sten/min_sten.le.(1.4)))then
	   ielem_vortex(1,k)=100.0
	   
                                        
                                        
                                                    do j=icomp_set,iselemt(n)-ibn
                                                    do jk=iselemt(n)-icomp_set,j+1,-1
                                                if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                                                pik= ilocalallelg3(1,jk-1)
                                                distf=ilocalallelgd(1,jk-1)
                                                    ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                                                ilocalallelg3(1,jk) = pik
                                                    ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                                                ilocalallelgd(1,jk)=distf
                                                endif
                                            enddo ! j
                                            enddo  ! i
                                            
                                            
                max_sten2=((testdist)*(iorder+1))  
                
	     do kk=icomp_set,iselemt(n)-ibn
                    if  ((ilocalallelgd(1,kk)).le.max_sten2)then
                         
                        kxk=kxk+1
                       
                    end if
                    end do
                    
 	    
	    end if
	    test6=ielem_inumneighbours(k)*2.6
	    
	    ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,test6))
 	  imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
	   ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))
	   
	   
	   
	   
	   end if
	      
	      
	     if ((icompact.eq.1).and.((max_sten/min_sten).le.1.2))then   

            extf2=extf          !extension
	    icomp_set=(ilx*extf2/2)+3      
	   
	    kxk=0
	
                    do j=icomp_set,iselemt(n)-ibn
                    do jk=iselemt(n)-icomp_set-ibn,j+1,-1
                if (ilocalallelgd(1,jk-1) .gt. ilocalallelgd(1,jk)) then
                pik= ilocalallelg3(1,jk-1)
                distf=ilocalallelgd(1,jk-1)
                    ilocalallelg3(1,jk-1) = ilocalallelg3(1,jk)
                ilocalallelg3(1,jk) = pik
                    ilocalallelgd(1,jk-1)=ilocalallelgd(1,jk)
                ilocalallelgd(1,jk)=distf
                endif
            enddo ! j
            enddo  ! i
            
            
                   do kk=icomp_set,iselemt(n)-ibn
                    if  ((ilocalallelgd(1,kk)).le.max_sten2)then
!                         
                        kxk=kxk+1
                       
                    end if
                    end do
!  	

         ielem_inumneighbours(k)=max(ielem_inumneighbours(k),min(kxk+icomp_set,ilx*extf2))
 	  imaxdegfree=max(imaxdegfree,ielem_inumneighbours(k)-1)
	   ilocalallelg(n,k,1,icomp_set:iselemt(n))=ilocalallelg3(1,icomp_set:iselemt(n))

 	 
!  	 end if
 	 
	    
	end if
	end if
        
       
    end if


	
end do
!$omp end do
    
 
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!

end subroutine detsten
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine allsf(stcon,stconc,stcons,stcong,isosa,ifsat,iistart,ix)
!> @brief
!> this recursive subroutine builds the large stencils suitable for periodic boundaries
implicit none
integer::i,k,jjj,kj,j,l,im,io,ir,kmaxe,icpuid,itrue,jloop,flag_per
integer,dimension(1),intent(inout)::stcon,stconc,stcons,stcong,isosa,iistart,ifsat,ix
kmaxe=xmpielrank(n)
      if (dimensiona.eq.3)then
			  jloop=6
		else
			  jloop=4
		end if
		icpuid=n
		flag_per=ilocalallelgper(n,stcon(1),1,iistart(1))
		iistart(1)=iistart(1)+1
		if (xmpie(stcong(1)).eq.n)then
			l=xmpil(stcong(1))
			stconc(1)=l
		
		do j=1,ielem_ifca(stconc(1))
			if (ielem_ineighg(j,stconc(1)).gt.0)then
			ix(1)=ielem_ineighg(j,stconc(1))
			call check(n,stcon,ix,ifsat)
			if (ifsat(1).eq.1)then
			if (isosa(1).le.iselemt(n)-1)then
			isosa(1)=isosa(1)+1
			ilocalallelg(n,stcon(1),1,isosa(1))=ielem_ineighg(j,stconc(1))
			if (dimensiona.eq.3)then
			if (flag_per.eq.0) then
                if (ielem_interior(stconc(1)).eq.1)then
                if(ielem_ibounds(j,stconc(1)).gt.0)then
                if((ibound_icode(ielem_ibounds(j,stconc(1))).eq.5)&
                            .or.(ibound_icode(ielem_ibounds(j,stconc(1))).eq.50)) then
                    if (ibound_icode(ielem_ibounds(j,stconc(1))).eq.5) then
                        ilocalallelgper(n,stcon(1),1,isosa(1))=1
                    else
                        ilocalallelgper(n,stcon(1),1,isosa(1))=2
                    end if
                else
                    ilocalallelgper(n,stcon(1),1,isosa(1))=0
                end if
            end if
            end if
			else if (flag_per.eq.1) then
                if (ielem_interior(stconc(1)).eq.1)then
                if (ielem_ibounds(j,stconc(1)).gt.0)then
                if((ibound_icode(ielem_ibounds(j,stconc(1))).eq.5)&
                            .or.(ibound_icode(ielem_ibounds(j,stconc(1))).eq.50)) then
                    if (ibound_icode(ielem_ibounds(j,stconc(1))).eq.5) then
                        ilocalallelgper(n,stcon(1),1,isosa(1))=1
                    else
                        ilocalallelgper(n,stcon(1),1,isosa(1))=0
                    end if
                else
                    ilocalallelgper(n,stcon(1),1,isosa(1))=1
                end if
                end if
                end if
            else
                if (ielem_interior(stconc(1)).eq.1)then
                if(ielem_ibounds(j,stconc(1)).gt.0)then
                if((ibound_icode(ielem_ibounds(j,stconc(1))).eq.5)&
                                .or.(ibound_icode(ielem_ibounds(j,stconc(1))).eq.50)) then
                    if (ibound_icode(ielem_ibounds(j,stconc(1))).eq.5) then
                        ilocalallelgper(n,stcon(1),1,isosa(1))=0
                    else
                        ilocalallelgper(n,stcon(1),1,isosa(1))=2
                    end if
                else
                    ilocalallelgper(n,stcon(1),1,isosa(1))=2
                end if
			end if
			end if
			end if
			end if
			end if
			end if
			end if
! 			
		end do
		end if
		if (xmpie(stcong(1)).ne.n)then
		
		do j=1,jloop
			if(glneigh(stcong(1),j).gt.0) then
			
			ix(1)=glneigh(stcong(1),j)
			call check(n,stcon,ix,ifsat)

			if (ifsat(1).eq.1)then
			if (isosa(1).le.iselemt(n)-1)then
			isosa(1)=isosa(1)+1
			ilocalallelg(n,stcon(1),1,isosa(1))=glneigh(stcong(1),j)
			if (dimensiona.eq.3)then
			if (flag_per.eq.0) then
                if (glneighper(stcong(1),j).eq.1) then
                    ilocalallelgper(n,stcon(1),1,isosa(1))=1
                else if (glneighper(stcong(1),j).eq.2) then
                    ilocalallelgper(n,stcon(1),1,isosa(1))=2
                else
                    ilocalallelgper(n,stcon(1),1,isosa(1))=0
                end if
			else if (flag_per.eq.1) then
                if (glneighper(stcong(1),j).eq.1) then
                    ilocalallelgper(n,stcon(1),1,isosa(1))=1
                else if (glneighper(stcong(1),j).eq.2) then
                    ilocalallelgper(n,stcon(1),1,isosa(1))=0
                else
                    ilocalallelgper(n,stcon(1),1,isosa(1))=1
                end if
            else 
                if (glneighper(stcong(1),j).eq.1) then
                    ilocalallelgper(n,stcon(1),1,isosa(1))=0
                else if (glneighper(stcong(1),j).eq.2) then
                    ilocalallelgper(n,stcon(1),1,isosa(1))=2
                else
                    ilocalallelgper(n,stcon(1),1,isosa(1))=2
                end if
			end if
			end if
			end if
			end if
			end if
		end do
		end if
		if (isosa(1).lt.iselemt(n))then
			!-------------------for debugging only -----------------------------------------!
! 			
			!-------------------for debugging only -----------------------------------------!
			stcong(1)=ilocalallelg(n,stcon(1),1,iistart(1))
			!-------------------for debugging only -----------------------------------------!
! 			
			!-------------------for debugging only -----------------------------------------!
			call allsf (stcon,stconc,stcons,stcong,isosa,ifsat,iistart,ix)
		end if
		if (isosa(1).eq.iselemt(n))then
			return
		end if

			
end subroutine allsf



subroutine check(n,stcon,ix,ifsat)
!> @brief
!> this subroutine checks if some candidate elements already belong to the an existing list of neighbours
implicit none
integer,intent(in)::n
integer,dimension(1),intent(inout)::stcon,ix,ifsat
integer::i,j,k
ifsat(1)=1
do i=1,iselemt(n)
if (ix(1).eq.ilocalallelg(n,stcon(1),1,i))then
	ifsat(1)=0
end if
end do
end subroutine check













subroutine check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
!> @brief
!> this subroutine checks which candidate cells satisfy the directionality condition for directional stencils
implicit none
integer,intent(in)::n,n_node,is_periodic
integer,intent(inout)::isatisfied
real,dimension(1:dimensiona),intent(in)::vg,bc
real,dimension(1:dimensiona)::xcc,vgg
real,dimension(1:dimensiona,1:dimensiona)::vva1
real,dimension(1:8,1:dimensiona),intent(inout)::vext
real,dimension(1)::deta
real::small,tempxx
integer::isat1,isat2,isat4,isat3




small=tolsmall

isatisfied=0

select case(n_node)

case (4)
        
	 vext(1,:)=bc(:)
	 vext(2,1:3)=vext(6,1:3); vext(3,1:3)=vext(5,1:3);vext(4,1:3)=vext(7,1:3);  
	vgg(:)=vg(:)
	call computejacobians(n,vext,vva1,deta)
	 if (iperiodicity.eq.1) then
	 if(per_rot.eq.0)then
	 if (abs(vg(1) - bc(1)) .ge. xper/2.0d0)    vgg(1) = vg(1) + xper*sign(1.d0,bc(1) - xper/2.0d0)
          if (abs(vg(2) - bc(2)) .ge. yper/2.0d0)    vgg(2) = vg(2) + yper*sign(1.d0,bc(2) - yper/2.0d0)
	   if (abs(vg(3) - bc(3)) .ge. zper/2.0d0)    vgg(3) = vg(3) + zper*sign(1.d0,bc(3) - zper/2.0d0)
     else
        if (is_periodic.eq.2) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(angle_per)-vgg(2)*sin(angle_per)
            vgg(2)=tempxx*sin(angle_per)+vgg(2)*cos(angle_per)
        end if 
        if (is_periodic.eq.1) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(-angle_per)-vgg(2)*sin(-angle_per)
            vgg(2)=tempxx*sin(-angle_per)+vgg(2)*cos(-angle_per)
        end if 
     end if
    end if
          xcc(1:3) =  matmul(vva1(1:3,1:3),vgg(1:3) -vext(1,1:3))
          
           if ((xcc(1).ge.zero).and.(xcc(2).ge.zero).and.(xcc(3).ge.zero))then
	    isat1=1
	  else
	      isat1=0
	  end if
	  

	  
	  vext(1,:)=bc(:)
	 vext(2,1:3)=vext(6,1:3); vext(3,1:3)=vext(8,1:3);vext(4,1:3)=vext(7,1:3);  
	 call computejacobians(n,vext,vva1,deta)
	 if (iperiodicity.eq.1) then
	 if(per_rot.eq.0)then
	 if (abs(vg(1) - bc(1)) .ge. xper/2.0d0)    vgg(1) = vg(1) + xper*sign(1.d0,bc(1) - xper/2.0d0)
          if (abs(vg(2) - bc(2)) .ge. yper/2.0d0)    vgg(2) = vg(2) + yper*sign(1.d0,bc(2) - yper/2.0d0)
	   if (abs(vg(3) - bc(3)) .ge. zper/2.0d0)    vgg(3) = vg(3) + zper*sign(1.d0,bc(3) - zper/2.0d0)
     else
        if (is_periodic.eq.2) then
	    tempxx=vg(1)
	    vgg(1)=tempxx*cos(angle_per)-vgg(2)*sin(angle_per)
	    vgg(2)=tempxx*sin(angle_per)+vgg(2)*cos(angle_per)
        end if 
	if (is_periodic.eq.1) then
	    tempxx=vg(1)
	    vgg(1)=tempxx*cos(-angle_per)-vgg(2)*sin(-angle_per)
	    vgg(2)=tempxx*sin(-angle_per)+vgg(2)*cos(-angle_per)
        end if 
     end if
        end if
          xcc(1:3) =  matmul(vva1(1:3,1:3),vgg(1:3) -vext(1,1:3))
          if ((xcc(1).ge.zero).and.(xcc(2).ge.zero).and.(xcc(3).ge.zero))then
	    isat2=1
	  else
	      isat2=0
	  end if
	  
	  
	  
	   vext(1,:)=bc(:)
	 vext(2,1:3)=vext(5,1:3); vext(3,1:3)=vext(8,1:3);vext(4,1:3)=vext(7,1:3);  
	 call computejacobians(n,vext,vva1,deta)
	 if (iperiodicity.eq.1) then
	 if(per_rot.eq.0)then
	 if (abs(vg(1) - bc(1)) .ge. xper/2.0d0)    vgg(1) = vg(1) + xper*sign(1.d0,bc(1) - xper/2.0d0)
          if (abs(vg(2) - bc(2)) .ge. yper/2.0d0)    vgg(2) = vg(2) + yper*sign(1.d0,bc(2) - yper/2.0d0)
	   if (abs(vg(3) - bc(3)) .ge. zper/2.0d0)    vgg(3) = vg(3) + zper*sign(1.d0,bc(3) - zper/2.0d0)
     else
        if (is_periodic.eq.2) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(angle_per)-vgg(2)*sin(angle_per)
            vgg(2)=tempxx*sin(angle_per)+vgg(2)*cos(angle_per)
        end if 
        if (is_periodic.eq.1) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(-angle_per)-vgg(2)*sin(-angle_per)
            vgg(2)=tempxx*sin(-angle_per)+vgg(2)*cos(-angle_per)
        end if 
     end if
        end if
          xcc(1:3) =  matmul(vva1(1:3,1:3),vgg(1:3) -vext(1,1:3))
          if ((xcc(1).ge.zero).and.(xcc(2).ge.zero).and.(xcc(3).ge.zero))then
	    isat3=1
	  else
	      isat3=0
	  end if
	  
	  
          
          vext(1,:)=bc(:)
	 vext(2,1:3)=vext(5,1:3); vext(3,1:3)=vext(8,1:3);vext(4,1:3)=vext(6,1:3);  
	 call computejacobians(n,vext,vva1,deta)
	 if (iperiodicity.eq.1) then
	 if(per_rot.eq.0)then
	 if (abs(vg(1) - bc(1)) .ge. xper/2.0d0)    vgg(1) = vg(1) + xper*sign(1.d0,bc(1) - xper/2.0d0)
          if (abs(vg(2) - bc(2)) .ge. yper/2.0d0)    vgg(2) = vg(2) + yper*sign(1.d0,bc(2) - yper/2.0d0)
	   if (abs(vg(3) - bc(3)) .ge. zper/2.0d0)    vgg(3) = vg(3) + zper*sign(1.d0,bc(3) - zper/2.0d0)
     else
        if (is_periodic.eq.2) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(angle_per)-vgg(2)*sin(angle_per)
            vgg(2)=tempxx*sin(angle_per)+vgg(2)*cos(angle_per)
        end if 
        if (is_periodic.eq.1) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(-angle_per)-vgg(2)*sin(-angle_per)
            vgg(2)=tempxx*sin(-angle_per)+vgg(2)*cos(-angle_per)
        end if 
     end if
        end if
          xcc(1:3) =  matmul(vva1(1:3,1:3),vgg(1:3) -vext(1,1:3))
          if ((xcc(1).ge.zero).and.(xcc(2).ge.zero).and.(xcc(3).ge.zero))then
	    isat4=1
	  else
	      isat4=0
	  end if
          
          
	 if ((isat1.eq.1).or.(isat2.eq.1).or.(isat3.eq.1).or.(isat4.eq.1))then
	    isatisfied=1
	  else
	      isatisfied=0
	  end if
           


case (3)
vext(1,:)=bc(:)
	vgg(:)=vg(:)
	call computejacobians(n,vext,vva1,deta)
	 if (iperiodicity.eq.1) then
	 if(per_rot.eq.0)then
	 if (abs(vg(1) - bc(1)) .ge. xper/2.d0)    vgg(1) = vg(1) + xper*sign(1.d0,bc(1) - xper/2.d0)
          if (abs(vg(2) - bc(2)) .ge. yper/2.d0)    vgg(2) = vg(2) + yper*sign(1.d0,bc(2) - yper/2.d0)
	   if (abs(vg(3) - bc(3)) .ge. zper/2.d0)    vgg(3) = vg(3) + zper*sign(1.d0,bc(3) - zper/2.d0)
      else
        if (is_periodic.eq.2) then
	    
            tempxx=vg(1)
            vgg(1)=tempxx*cos(angle_per)-vgg(2)*sin(angle_per)
            vgg(2)=tempxx*sin(angle_per)+vgg(2)*cos(angle_per)
        end if 
        if (is_periodic.eq.1) then
            tempxx=vg(1)
            vgg(1)=tempxx*cos(-angle_per)-vgg(2)*sin(-angle_per)
            vgg(2)=tempxx*sin(-angle_per)+vgg(2)*cos(-angle_per)
        end if
      end if
        end if
          xcc =  matmul(vva1(1:3,1:3),vgg(1:3) -vext(1,1:3))
          if ((xcc(1).ge.zero).and.(xcc(2).ge.zero).and.(xcc(3).ge.zero))then
	    isat1=1
	  else
	      isat1=0
	  end if

 if ((isat1.eq.1))then
	    isatisfied=1
	  else
	      isatisfied=0
	  end if


case(2)


       vext(1,1:2)=bc(1:2)
	vgg(1:2)=vg(1:2)
	call computejacobians2(n,vext,vva1,deta)
	 if (iperiodicity.eq.1) then
	 if (abs(vg(1) - bc(1)) .ge. xper/2.0d0)    vgg(1) = vg(1) + xper*sign(1.0d0,bc(1) - xper/2.0d0)
          if (abs(vg(2) - bc(2)) .ge. yper/2.0d0)    vgg(2) = vg(2) + yper*sign(1.0d0,bc(2) - yper/2.0d0)
	 
         
        end if
          xcc(1:2) =  matmul(vva1(1:2,1:2),vgg(1:2) -vext(1,1:2))
          if ((xcc(1).ge.zero).and.(xcc(2).ge.zero))then
	    isat1=1
	  else
	      isat1=0
	  end if

 if ((isat1.eq.1))then
	    isatisfied=1
	  else
	      isatisfied=0
	  end if


end select





! 
 end subroutine check_condition





subroutine sortstencils(n)
!> @brief
!> this subroutine sorts the stencils with respect to their distance from the cell-centre
implicit none
integer,intent(in)::n
real,dimension(iselemt(n))::rdistl
integer,dimension(iselemt(n))::idistl
real::minds,xc,yc,zc,xl,yl,zl
integer::i,j,k,l,m,kmaxe
kmaxe=xmpielrank(n)
do i=1,kmaxe
	minds=0.0
	  xc=ielem_xxc(i)
	  yc=ielem_yyc(i)
	  zc=ielem_zzc(i)
	do j=2,iselemt(n)

	    if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	      xl=centerr((ilocalallelg(n,i,1,j)),1)
	      yl=centerr((ilocalallelg(n,i,1,j)),2)
	      zl=centerr((ilocalallelg(n,i,1,j)),3)
	    end if
	    if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then
	      xl=ielem_xxc(xmpil((ilocalallelg(n,i,1,j))))
	      yl=ielem_yyc(xmpil((ilocalallelg(n,i,1,j))))
	      zl=ielem_zzc(xmpil((ilocalallelg(n,i,1,j))))
	    end if
	      rdistl(j)=sqrt(((xl-xc)**2.0)+((yl-yc)**2.0)+((zl-zc)**2.0))
		
	end do
	do k=2,iselemt(n)
	do j=2,iselemt(n)-1
	if (rdistl(j)>rdistl(j+1))then
	  m=ilocalallelg(n,i,1,j)
	  ilocalallelg(n,i,1,j)=ilocalallelg(n,i,1,j+1)
	  ilocalallelg(n,i,1,j+1)=m
	end if
   enddo 
  enddo  

end do








end subroutine sortstencils




subroutine stenciils_eesx(n)
!> @brief
!> this subroutine builds all the directional stencils from the large stencil
implicit none
integer,intent(in)::n
real,dimension(1:dimensiona)::vg,bc
integer::iwhichsten,ishyape,isatisfied,iconsi,n_node
integer::i,j,k,l,m,o,p,kmaxe,icount,icpuid,iatrue,ig,il,ifg,stnsha,itgh,kk,kxk,ix,ifvs,ixcz,iadd,iadd2,iadd3,itarget,iaddx,iaddx1,ifno,igvd
integer::candid,candid2
real,dimension(1:8,1:dimensiona)::vext
integer::is_periodic


kmaxe=xmpielrank(n)

if (dimensiona.eq.3)then
!$omp do
do i=1,kmaxe
! 
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
		
		
			stnsha=ielem_ifca(i)
			    iconsi=i
		     bc(1)=ielem_xxc(i);	bc(2)=ielem_yyc(i);	bc(3)=ielem_zzc(i);
	    
		ishyape=ielem_ishape(i)
		
		 il=0
			do iaddx=1,stnsha	!for all stencils
			
			
			if (ielem_types_faces(iaddx,i).eq.5)then
			
			igvd=2
			
			
			else
			
			igvd=2
			
			
			
			end if
			
			do iaddx1=1,igvd
			il=il+1
			ifno=3
			
			if (ielem_types_faces(iaddx,i).eq.5)then
			
			
			if (iaddx1.eq.1)then
			
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
				
						
			else
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,4,i))%cord(1:3)
			end if
		
		
			else
			
			if (iaddx1.eq.1)then
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3)
			vext(4,1:3)=(dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)+dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3))/2.0d0
			
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(4,1:3))/2.0d0
			
			
			
			end if
			if (iaddx1.eq.2)then
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=(dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)+dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3))/2.0d0
			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
			
			
			
			end if
			
! 			if (iaddx1.eq.3)then
! 			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
! 			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
! 			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3)
! 			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
! 			
! 			
! 			
! 			end if
			
			
			
			
			
			
			end if
					




			if ((ielem_ineighg(iaddx,i).gt.0))then
			iwhichsten=il
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then
  !do ifg=1,kmaxe
	  ifg=xmpil((ilocalallelg(n,i,1,j)))
	  iconsi=ifg
	   vg(1)=ielem_xxc(iconsi) ;vg(2)=ielem_yyc(i);vg(3)=ielem_zzc(i)
		     n_node=ifno

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.ielem_inumneighbours(i))then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.ielem_inumneighbours(i))then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	  
	  candid2=ilocalallelg(n,i,1,j)
	  if (cand(xmpie(candid2)).gt.0)then
                 
                 vg(1:3)=xand2r(cand(xmpie(candid2)),xmpil(candid2),1:3)
	  
	  


	    call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
	    if (isatisfied.eq.1)then
	    if (itgh+1.le.ielem_inumneighbours(i))then
	    itgh=itgh+1
	    ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
	    end if
	    end if
	    if (itgh.eq.ielem_inumneighbours(i))then
	    exit
	    end if
	    
	    
	    else
	     write(99+n,*)"killed here 3 due to stencils requiring more cpus"
							    stop
	    
	    
	    end if
  
  end if
end if
					end do		!elements in stencil
			      end if
			end do		!directional stencils
			end do
	
	end do
	!$omp end do
end if



end if

if (dimensiona.eq.2)then

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe
! 	
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
		
		
			stnsha=ielem_ifca(i)
			    iconsi=i
		    iconsi=i
		 bc(1)=ielem_xxc(i);	bc(2)=ielem_yyc(i);
	    
		ishyape=ielem_ishape(i)
		  il=0
			do iaddx=1,stnsha	!for all stencils
			
			do iaddx1=1,2
			il=il+1
			ifno=2
			
			if (iaddx1.eq.1)then
			vext(2,1:2)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:2)
			vext(3,1:2)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:2)
			vext(3,1:2)=(vext(2,1:2)+vext(3,1:2))*oo2
			
			else
			vext(2,1:2)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:2)
			vext(3,1:2)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:2)
			vext(2,1:2)=(vext(2,1:2)+vext(3,1:2))*oo2
			end if
			
			
			
			if ((ielem_ineighg(iaddx,i).gt.0))then
			iwhichsten=il
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then

	  ifg=xmpil((ilocalallelg(n,i,1,j)))
	  iconsi=ifg
	    vg(1)=ielem_xxc(iconsi) ;vg(2)=ielem_yyc(i)
		     n_node=ifno

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.ielem_inumneighbours(i))then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.ielem_inumneighbours(i))then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	 candid2=ilocalallelg(n,i,1,j)
	  if (cand(xmpie(candid2)).gt.0)then
                 
                 vg(1:2)=xand2r(cand(xmpie(candid2)),xmpil(candid2),1:2)
! 	  vg(n,3)=centerr((ilocalallelg(n,i,1,j)),3)

	      call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
	      if (isatisfied.eq.1)then
	      if (itgh+1.le.ielem_inumneighbours(i))then
	      itgh=itgh+1
	      ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
	      end if
	      end if
	      if (itgh.eq.ielem_inumneighbours(i))then
	      exit
	      end if
	      else
	      
	      write(99+n,*)"killed here 2 due to stencils requiring more cpus"
							    stop
	      
	      end if
  
  end if
end if
					end do		!elements in stencil
			      end if
			      
			end do		!directional stencils
	
	
	
	end do
	end do
	!$omp end do
end if


















end if






end subroutine stenciils_eesx


subroutine stenciils_ees(n)
!> @brief
!> this subroutine builds all the directional stencils from the large stencil (suitable for periodic boundaries)
implicit none
integer,intent(in)::n
real,dimension(1:dimensiona)::vg,bc
integer::iwhichsten,ishyape,isatisfied,iconsi,n_node
integer::i,j,k,l,m,o,p,kmaxe,icount,icpuid,iatrue,ig,il,ifg,stnsha,itgh,kk,kxk,ix,ifvs,ixcz,iadd,iadd2,iadd3,itarget,iaddx,iaddx1,ifno,igvd
integer::candid,candid2
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dimensiona)::cords
integer::is_periodic
kmaxe=xmpielrank(n)
is_periodic=0
if (dimensiona.eq.3)then

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe
! 	
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
		
		
			stnsha=ielem_ifca(i)

		      call compute_centre3d(i,cords)
			
		bc(1)=cords(1)   ;   bc(2)=cords(2);    bc(3)=cords(3)
	    
		ishyape=ielem_ishape(i)
		
		 il=0
			do iaddx=1,stnsha	!for all stencils
			
			
			if (ielem_types_faces(iaddx,i).eq.5)then
			
			igvd=2
			
			
			else
			
			igvd=2
			
			
			
			end if
			
			do iaddx1=1,igvd
			il=il+1
			ifno=3
			
			if (ielem_types_faces(iaddx,i).eq.5)then
			
			
			if (iaddx1.eq.1)then
			
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
				
						
			else
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,4,i))%cord(1:3)
			end if
		
		
			else
			
			if (iaddx1.eq.1)then
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3)
			vext(4,1:3)=(dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)+dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3))/2.0d0
			
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(4,1:3))/2.0d0
			
			
			
			end if
			if (iaddx1.eq.2)then
			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
			vext(3,1:3)=(dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)+dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3))/2.0d0
			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
			
			
			
			end if
			
! 			if (iaddx1.eq.3)then
! 			vext(2,1:3)=dinoder(ielem_nodes_faces(iaddx,3,i))%cord(1:3)
! 			vext(3,1:3)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:3)
! 			vext(4,1:3)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:3)
! 			
! 			vext(4,1:3)=(vext(2,1:3)+vext(3,1:3)+vext(4,1:3))/3.0d0
! 			
! 			
! 			
! 			end if
			
			
			
			
			
			
			end if
					

! 			if (ifno.eq.3)then
! 			vext(2,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
! 			vext(3,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
! 			vext(4,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
! 			else
! 			vext(2,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
! 			vext(3,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
! 			vext(4,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
! 			vext(5,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
! 			vext(6,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
! 			vext(7,1:3)=dinoder(ielem_nodes_faces(il,4,i))%cord(1:3)
! 
! 			end if


			if ((ielem_ineighg(iaddx,i).gt.0))then
			iwhichsten=il
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then
  !do ifg=1,kmaxe
	  ifg=xmpil((ilocalallelg(n,i,1,j)))

	  call compute_centre3d(ifg,cords)
	  vg(1)=cords(1) ;vg(2)=cords(2) ; vg(3)=cords(3)
		     n_node=ifno

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.ielem_inumneighbours(i))then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.ielem_inumneighbours(i))then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	    vg(1)=centerr((ilocalallelg(n,i,1,j)),1)
	  vg(2)=centerr((ilocalallelg(n,i,1,j)),2)
	  vg(3)=centerr((ilocalallelg(n,i,1,j)),3)

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.ielem_inumneighbours(i))then
  itgh=itgh+1
  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.ielem_inumneighbours(i))then
  exit
  end if
  
  end if
end if
					end do		!elements in stencil
			      end if
			end do		!directional stencils
			end do
	
	end do
	!$omp end do
end if



end if

if (dimensiona.eq.2)then

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe

	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
		
		
			stnsha=ielem_ifca(i)

		      call compute_centre2d(i,cords)
			
		bc(1)=cords(1)   ;   bc(2)=cords(2);    !bc(n,3)=cords(3)
	    
		ishyape=ielem_ishape(i)
		  il=0
			do iaddx=1,stnsha	!for all stencils
			
			do iaddx1=1,2
			il=il+1
			ifno=2
			
			if (iaddx1.eq.1)then
			vext(2,1:2)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:2)
			vext(3,1:2)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:2)
			vext(3,1:2)=(vext(2,1:2)+vext(3,1:2))*oo2
			
			else
			vext(2,1:2)=dinoder(ielem_nodes_faces(iaddx,1,i))%cord(1:2)
			vext(3,1:2)=dinoder(ielem_nodes_faces(iaddx,2,i))%cord(1:2)
			vext(2,1:2)=(vext(2,1:2)+vext(3,1:2))*oo2
			end if
			
			
			
			if ((ielem_ineighg(iaddx,i).gt.0))then
			iwhichsten=il
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then

	  ifg=xmpil((ilocalallelg(n,i,1,j)))

	  call compute_centre2d(ifg,cords)
	  vg(1)=cords(1) ;vg(2)=cords(2)
		     n_node=ifno

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.ielem_inumneighbours(i))then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.ielem_inumneighbours(i))then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	    vg(1)=centerr((ilocalallelg(n,i,1,j)),1)
	  vg(2)=centerr((ilocalallelg(n,i,1,j)),2)
! 	  vg(n,3)=centerr((ilocalallelg(n,i,1,j)),3)

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.ielem_inumneighbours(i))then
  itgh=itgh+1
  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.ielem_inumneighbours(i))then
  exit
  end if
  
  end if
end if
					end do		!elements in stencil
			      end if
			      
			end do		!directional stencils
	
	
	
	end do
	end do
	!$omp end do
end if


















end if






end subroutine stenciils_ees





subroutine stenciilsx(n)
!> @brief
!> this subroutine builds all the directional stencils from the large stencil based on various aglorithms
implicit none
integer,intent(in)::n
real,dimension(1:dimensiona)::vg,bc
integer::iwhichsten,ishyape,isatisfied,iconsi,n_node
integer::i,j,k,l,m,o,p,kmaxe,icount,icpuid,iatrue,ig,il,ifg,stnsha,itgh,kk,kxk,ix,ifvs,ixcz,iadd,iadd2,iadd3,itarget,iaddx,iaddx1,ifno
integer::candid,candid2
real,dimension(1:8,1:dimensiona)::vext
integer::is_periodic
is_periodic=0
kmaxe=xmpielrank(n)

if (dimensiona.eq.3)then

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe
! 	
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
                 if (ees.eq.5)then
                   itarget=numneighbours2
                else 
                   itarget=ielem_inumneighbours(i)
                end if
                 
                    
		
			stnsha=ielem_ifca(i)
			    iconsi=i

		bc(1)=ielem_xxc(i);	bc(2)=ielem_yyc(i);	bc(3)=ielem_zzc(i);

	    
		ishyape=ielem_ishape(i)
			do il=1,stnsha	!for all stencils
			   if (ielem_types_faces(il,i).eq.5)then
			ifno=4
			else
			ifno=3
			end if
			

			if (ifno.eq.3)then
			vext(2,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
			else
			vext(2,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
			vext(5,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
			vext(6,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
			vext(7,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
            vext(8,1:3)=dinoder(ielem_nodes_faces(il,4,i))%cord(1:3)
			end if


			if ((ielem_ineighg(il,i).gt.0))then
			iwhichsten=il
			
! 			
			
			
			
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						
						
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then

	  ifg=xmpil((ilocalallelg(n,i,1,j)))
	  iconsi=ifg
	  vg(1)=ielem_xxc(iconsi);	vg(2)=ielem_yyc(iconsi);	vg(3)=ielem_zzc(iconsi);

		     n_node=ifno

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.itarget)then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.itarget)then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	  candid2=ilocalallelg(n,i,1,j)
		if (cand(xmpie(candid2)).gt.0)then
                 
                 vg(1:3)=xand2r(cand(xmpie(candid2)),xmpil(candid2),1:3)
                

		call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
		      if (isatisfied.eq.1)then
		      if (itgh+1.le.itarget)then
		      itgh=itgh+1
		      ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
		      end if
		      end if
		      if (itgh.eq.itarget)then
		      exit
		      end if
		else
		
		write(99+n,*)"killed here2 due to stencils requiring more cpus"
							    stop
		
		
		end if
		end if
end if
                        
					end do		!elements in stencil
			      end if
			end do		!directional stencils
	
	end do
	!$omp end do
end if



end if

if (dimensiona.eq.2)then

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe
! 	
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
		
                         if (ees.eq.5)then
                   itarget=numneighbours2
                else 
                   itarget=ielem_inumneighbours(i)
                end if
			stnsha=ielem_ifca(i)
			    iconsi=i
		 bc(1)=ielem_xxc(i);	bc(2)=ielem_yyc(i);
	    

			do il=1,stnsha	!for all stencils
			ifno=2
			vext(2,1:2)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:2)
			vext(3,1:2)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:2)

! 			


			if ((ielem_ineighg(il,i).gt.0))then
			iwhichsten=il
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
			
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then

	  ifg=xmpil((ilocalallelg(n,i,1,j)))
	  iconsi=ifg
	   vg(1)=ielem_xxc(iconsi) ;vg(2)=ielem_yyc(iconsi)

		     n_node=ifno

  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.itarget)then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.itarget)then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	  candid2=ilocalallelg(n,i,1,j)
	  if (cand(xmpie(candid2)).gt.0)then
                 
                 vg(1:2)=xand2r(cand(xmpie(candid2)),xmpil(candid2),1:2)
	  

	  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
	  if (isatisfied.eq.1)then
	  if (itgh+1.le.itarget)then
	  itgh=itgh+1
	  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
	  end if
	  end if
	  if (itgh.eq.itarget)then
	  exit
	  end if
  
	  else
	  
	  write(99+n,*)"killed here 2 due to stencils requiring more cpus"
							    stop
	  
	  end if
  end if
end if
					end do		!elements in stencil
			      end if
			end do		!directional stencils
! 	end do
	end do
	!$omp end do
end if

















end if






end subroutine stenciilsx


subroutine stenciils(n)
!> @brief
!> this subroutine builds all the directional stencils from the large stencil based on various algorithms (suitable for period boundaries)
implicit none
integer,intent(in)::n
real,dimension(1:dimensiona)::vg,bc
integer::iwhichsten,ishyape,isatisfied,iconsi,n_node
integer::i,j,k,l,m,o,p,kmaxe,icount,icpuid,iatrue,ig,il,ifg,stnsha,itgh,kk,kxk,ix,ifvs,ixcz,iadd,iadd2,iadd3,itarget,iaddx,iaddx1,ifno
integer::candid,candid2
real,dimension(1:8,1:dimensiona)::vext
real,dimension(1:dimensiona)::cords
integer::is_periodic
is_periodic=0
kmaxe=xmpielrank(n)

if (dimensiona.eq.3)then

!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe
! 	
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
        ilocalstencilper(n,i,1,j)=ilocalallelgper(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
                 if (ees.eq.5)then
                   itarget=numneighbours2
                else 
                   itarget=ielem_inumneighbours(i)
                end if
                 
                    
		
			stnsha=ielem_ifca(i)

		      call compute_centre3d(i,cords)
			
		bc(1)=cords(1)   ;   bc(2)=cords(2);    bc(3)=cords(3)
	    
		ishyape=ielem_ishape(i)
			do il=1,stnsha	!for all stencils
			   if (ielem_types_faces(il,i).eq.5)then
			ifno=4
			else
			ifno=3
			end if
			

			if (ifno.eq.3)then
			vext(2,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
			else
			vext(2,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
			vext(3,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
			vext(4,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
			vext(5,1:3)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:3)
			vext(6,1:3)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:3)
			vext(7,1:3)=dinoder(ielem_nodes_faces(il,3,i))%cord(1:3)
                        vext(8,1:3)=dinoder(ielem_nodes_faces(il,4,i))%cord(1:3)
			end if


			if ((ielem_ineighg(il,i).gt.0))then
			iwhichsten=il
			
! 	
			
			
			
			
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						
						
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then
  !do ifg=1,kmaxe
	  ifg=xmpil((ilocalallelg(n,i,1,j)))

	  call compute_centre3d(ifg,cords)
	  vg(1)=cords(1) ;vg(2)=cords(2) ; vg(3)=cords(3)
		     n_node=ifno
  is_periodic=ilocalallelgper(n,i,1,j)
  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.itarget)then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  ilocalstencilper(n,i,il+1,itgh)=ilocalallelgper(n,i,1,j)
  end if
  end if
  if (itgh.eq.itarget)then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	    vg(1)=centerr((ilocalallelg(n,i,1,j)),1)
	  vg(2)=centerr((ilocalallelg(n,i,1,j)),2)
	  vg(3)=centerr((ilocalallelg(n,i,1,j)),3)
  is_periodic=ilocalallelgper(n,i,1,j) 
  call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.itarget)then
  itgh=itgh+1
  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  ilocalstencilper(n,i,il+1,itgh)=ilocalallelgper(n,i,1,j)
  end if
  end if
  if (itgh.eq.itarget)then
  exit
  end if
  
  end if
end if
                        
					end do		!elements in stencil
			      end if
			end do		!directional stencils
	
	end do
	!$omp end do
end if



end if

if (dimensiona.eq.2)then

!-------------------for debugging only -----------------------------------------!
	
!-------------------for debugging only -----------------------------------------!
!$omp do
do i=1,kmaxe
! 	
	do j=1,ielem_inumneighbours(i)
	      
		ilocalstencil(n,i,1,j)=ilocalallelg(n,i,1,j)
! 		
	end do
	
end do
!$omp end do
if (typesten.gt.1)then
	icpuid=n
	!$omp do
	do i=1,kmaxe	!for all elements
		
                         if (ees.eq.5)then
                   itarget=numneighbours2
                else 
                   itarget=ielem_inumneighbours(i)
                end if
			stnsha=ielem_ifca(i)
			    iconsi=i
		      call compute_centre2d(i,cords)
			
		bc(1)=cords(1)   ;   bc(2)=cords(2);    !bc(n,3)=cords(3)
	    
		ishyape=ielem_ishape(i)
			do il=1,stnsha	!for all stencils
			ifno=2

			vext(2,1:2)=dinoder(ielem_nodes_faces(il,1,i))%cord(1:2)
			vext(3,1:2)=dinoder(ielem_nodes_faces(il,2,i))%cord(1:2)




 			if ((ielem_ineighg(il,i).gt.0))then
			iwhichsten=il
			
			ilocalstencil(n,i,il+1,1)=ilocalallelg(n,i,1,1)
			
					itgh=1
					do j=2,iselemt(n)!for all stencil elements
						if ((ilocalallelg(n,i,1,j)).gt.0) then
						isatisfied=0
if (xmpie(ilocalallelg(n,i,1,j)).eq.n)then

	  ifg=xmpil((ilocalallelg(n,i,1,j)))

	  call compute_centre2d(ifg,cords)
	  vg(1)=cords(1) ;vg(2)=cords(2)
		     n_node=ifno

   call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.itarget)then
  itgh=itgh+1

  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.itarget)then
  exit
  end if
  
  end if
  if (xmpie(ilocalallelg(n,i,1,j)).ne.n)then
	  n_node=ifno
	    vg(1)=centerr((ilocalallelg(n,i,1,j)),1)
	  vg(2)=centerr((ilocalallelg(n,i,1,j)),2)
! 	  vg(n,3)=centerr((ilocalallelg(n,i,1,j)),3)

   call check_condition(n,n_node,bc,vg,vext,is_periodic,isatisfied)
  if (isatisfied.eq.1)then
  if (itgh+1.le.itarget)then
  itgh=itgh+1
  ilocalstencil(n,i,il+1,itgh)=ilocalallelg(n,i,1,j)
  end if
  end if
  if (itgh.eq.itarget)then
  exit
  end if
  
  end if
end if
					end do		!elements in stencil
			      end if
			end do		!directional stencils
	
	end do
! 	end do
	!$omp end do
end if

















end if






end subroutine stenciils



! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!subroutine for determining the number of elements!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!required for each stencil!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!for various order of accuracy!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine determine_size(n,iorder,iselem,iselemt,ioverst,ioverto,ilx,numneighbours,idegfree,imaxdegfree,iextend)
!> @brief
!> this subroutine determines the degress of freedom, neighbours and polynomial order for each stencil of each cell
	implicit none
	integer,intent(inout)::iorder
	integer,allocatable,dimension(:),intent(inout)::iselemt
	integer,intent(out)::ioverst,iselem
	integer,intent(out)::ioverto
	integer,intent(out)::idegfree
	integer,intent(out)::imaxdegfree
	integer,intent(out)::ilx     !number of degrees of freedom
	integer,intent(inout)::numneighbours
	integer,intent(in)::iextend,n
	integer::i,itemd,ilxx
	
	allocate(iselemt(n:n))
	
	if (dimensiona.eq.3) then
        if (iorder.ge.2) then
            ilx=((iorder+1)*(iorder+2)*(iorder+3))/6

            idegfree=ilx-1
            numneighbours=ilx*extf
            imaxdegfree=numneighbours-1
            iselem=(ilx*extf)*iextend
            iselemt(n:n)=iselem
            ioverst=iselem
            ioverto=iselem
            
           if (dg.eq.1)then
                numneighbours=ilx*extf
                imaxdegfree=numneighbours-1
                itemd=(ilx*extf)*iextend
                iselem=min(itemd,imaxe-1)
                iselemt(n:n)=iselem
                ioverst=iselem
                ioverto=iselem
                idegfree3=(((iorder)*(iorder+1)*(iorder+2))/6)-1
            else
                numneighbours=ilx*extf
                imaxdegfree=numneighbours-1
                itemd=(ilx*extf)*iextend
                iselem=min(itemd,imaxe-1)
                iselemt(n:n)=iselem
                ioverst=iselem
                ioverto=iselem
            end if
            
            
            
            
            
            
        
            select case(iorder)
            
            case (1,2,3)
                idegfree2=3
                iorder2=1
                numneighbours2=12
            
            case(4,5,6,7)
                idegfree2=9
                iorder2=2
                numneighbours2=20
            
            end select
        end if
        
        if (iorder.eq.1)then
            ilx=((iorder+1)*(iorder+2)*(iorder+3))/6
            idegfree=ilx-1
            numneighbours=9
            imaxdegfree=numneighbours-1
            iselem=(ilx*extf)*iextend
            iselemt(n:n)=iselem
            ioverst=iselem
            ioverto=iselem
            
            idegfree2=3
            iorder2=1
            numneighbours2=7
            
            
            if (dg.eq.1)then
                idegfree=ilx-1
            numneighbours=9
            imaxdegfree=numneighbours-1
            iselem=(ilx*extf)*iextend
            iselemt(n:n)=iselem
            ioverst=iselem
            ioverto=iselem
            idegfree3=(((iorder)*(iorder+1)*(iorder+2))/6)-1
            end if
            
            
            
            
            
            
        end if
        
	else ! 2 dimensions
        if (iorder.ge.2)then
            ilx=((iorder+1)*(iorder+2))/2
            idegfree=ilx-1
        
            if (dg.eq.1)then
                numneighbours=ilx*extf
                imaxdegfree=numneighbours-1
                itemd=(ilx*extf)*iextend
                iselem=min(itemd,imaxe-1)
                iselemt(n:n)=iselem
                ioverst=iselem
                ioverto=iselem
                idegfree3=(((iorder)*(iorder+1))/2)-1
            else
                numneighbours=ilx*extf
                imaxdegfree=numneighbours-1
                itemd=(ilx*extf)*iextend
                iselem=min(itemd,imaxe-1)
                iselemt(n:n)=iselem
                ioverst=iselem
                ioverto=iselem
            end if
        
        
            select case(iorder)
            
            case (1,2,3)
                idegfree2=2
                iorder2=1
                numneighbours2=7
            
            
            case(4,5,6,7)
                idegfree2=5
                iorder2=2
                numneighbours2=11
            
            end select
        
        end if
        
        if (iorder.eq.1)then
            ilx=((iorder+1)*(iorder+2))/2
            idegfree=ilx-1
            
            
            if (dg.eq.1)then
                numneighbours=5
                imaxdegfree=numneighbours-1
                itemd=20
                iselem=min(itemd,imaxe-1)
                iselemt(n:n)=iselem
                ioverst=iselem
                ioverto=iselem
                idegfree3=(((iorder)*(iorder+1))/2)-1
            else
                numneighbours=ilx*extf
                imaxdegfree=numneighbours-1
                itemd=(ilx*extf)*iextend
                iselem=min(itemd,imaxe-1)
                iselemt(n:n)=iselem
                ioverst=iselem
                ioverto=iselem
            end if
            
            iorder2=1
            idegfree2=2
            numneighbours2=6
        end if
	end if
	
    do i=1,xmpielrank(n)
        ielem_inumneighbours(i)=numneighbours
        ielem_idegfree(i)=idegfree
        ielem_iorder(i)=iorder
    end do
        
        
    if (dg == 1) then
        num_dg_dofs = ilx
        if (dimensiona == 2) then 
            num_dg_reconstruct_dofs = (iorder + 2) * (iorder + 3) / 2 - num_dg_dofs
        else 
            num_dg_reconstruct_dofs = (iorder+2)*(iorder+3)*(iorder+4)/6 - num_dg_dofs
        end if
    end if


	end subroutine determine_size
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!subroutine called for reading dat file!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!for the simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!with information such as!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!order of accuracy methods etc!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------------------------------!	
subroutine gaussianpoints(igqrules,numberofpoints,numberofpoints2)
!> @brief
!> this subroutine determines quadrature points required for each spatial order of accuracy and for each cell type
implicit none
integer,intent(inout)::igqrules,numberofpoints,numberofpoints2
if (dimensiona.eq.3)then
if (igqrules.eq.1)then
	qp_hexa=1;qp_tetra=1;qp_pyra=1;qp_prism=1;
	qp_quad=1;qp_triangle=1
end if
if (igqrules.eq.2)then

	qp_hexa=8;qp_tetra=4;qp_pyra=5;qp_prism=6;
	qp_quad=4;qp_triangle=3
	
	
end if
if (igqrules.eq.3)then
	qp_hexa=27;qp_tetra=14;qp_pyra=15;qp_prism=18;
	qp_quad=9;qp_triangle=6
	
end if
if (igqrules.eq.4)then
	qp_hexa=64;qp_tetra=24;qp_pyra=15;qp_prism=24;
	qp_quad=16;qp_triangle=10
	
end if
if (igqrules.eq.5)then
	qp_hexa=125;qp_tetra=35;qp_pyra=15;qp_prism=50;
	qp_quad=25;qp_triangle=15
	
end if
if (igqrules.eq.6)then
	qp_hexa=216;qp_tetra=56;qp_pyra=15;qp_prism=60;
	qp_quad=36;qp_triangle=21
	
end if
if (igqrules.ge.7)then
	qp_hexa=216;qp_tetra=56;qp_pyra=15;qp_prism=60;
	qp_quad=36;qp_triangle=36
	
end if
if (reduce_comp.eq.1)then
qp_quad_n=1;qp_triangle_n=1
else
qp_quad_n=qp_quad;qp_triangle_n=qp_triangle
end if

else


if (igqrules.eq.1)then
	qp_quad=1;qp_triangle=1;qp_line=1
	
end if
if (igqrules.eq.2)then
	qp_quad=4;qp_triangle=3;qp_line=2
end if
if (igqrules.eq.3)then
	qp_quad=9;qp_triangle=6;qp_line=3
end if
if (igqrules.eq.4)then
	qp_quad=16;qp_triangle=10;qp_line=4
end if
if (igqrules.eq.5)then
	qp_quad=25;qp_triangle=15;qp_line=5
end if

if (igqrules.eq.6)then
	qp_quad=36;qp_triangle=21;qp_line=6
end if

if (igqrules.ge.7)then
	qp_quad=36;qp_triangle=36;qp_line=9
end if


if (reduce_comp.eq.1)then
qp_line_n=1
else
qp_line_n=qp_line
end if


end if
end subroutine gaussianpoints








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine writeout
!> @brief
!> this subroutine is solely for debugging purposes
integer::ijd,ic2,kmaxe,i
real,allocatable,dimension(:)::xx,yy,zz
kmaxe=xmpielrank(n)

if (n.eq.0)then
call open_input(n,itt)
	write(97,*)'title="final solution"'
 	write(97,*)'variables="x","y","z"'
 	
 	
write(97,*) 'zone n=',4,',e=',1,',zonetype = fequadrilateral,','datapacking = block'
allocate(xx(imaxn),yy(imaxn),zz(imaxn))
do i=1,kmaxe
    if (ielem_ishape(i).eq.6)then
	      do ijd=1,ielem_nonodes(i)
		write(97,*)dinoder(ielem_nodes(ijd,i))%cord(1)
	    end do
	    do ijd=1,ielem_nonodes(i)
		write(97,*)dinoder(ielem_nodes(ijd,i))%cord(2)
	    end do
! 	      do ijd=1,ielem_nonodes(i)
! 		write(97,*)dinoder(ielem_nodes(ijd,i))%cord(3)
! 	    end do
	    
    exit

  end if

end do

do i=1,kmaxe
if (ielem_ishape(i).eq.6)then

write(97,*)1,2,3,3
exit
end if

end do

call close_input(n,itt)
end if



call mpi_barrier(mpi_comm_world,ierror)
end subroutine





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine stencils(n,imaxe,xmpie,xmpielrank,ilocalstencil,typesten,numneighbours,restart)
!> @brief
!> this subroutine is establishing which of the stencils are admissible and which cells can use the weno algorithms
	implicit none
	integer,intent(in)::n,imaxe
	integer,allocatable,dimension(:),intent(in)::xmpie
	integer,allocatable,dimension(:),intent(in)::xmpielrank
	integer,allocatable,dimension(:,:,:,:),intent(inout)::ilocalstencil
	integer,intent(in)::numneighbours
	integer,intent(in)::typesten,restart
	integer::i,j,ji,k,lm,kmaxn,kk,kmaxe,iaa,kx,l,itrr,itrx,itry,itarget
	integer::it1,it2,it3,it4,it5,it6,it7,it8,itx,inx,itf,ifg,kmh,jng,ichg
	integer,dimension(12,1000)::items

	
	kmaxe=xmpielrank(n)
	call mpi_barrier(mpi_comm_world,ierror)
	
	do i=1,kmaxe
                 ielem_full(i)=0
                   
                  ielem_troubled(i)=1
                  
                  if(dg.eq.1)then
                  ielem_troubled(i)=0
                 
                  end if
			kx=0
		do lm=1,typesten
                                if ((lm.eq.1).or.(ees.ne.5))then
                                itarget=ielem_inumneighbours(i)
                                else
                                itarget=numneighbours2
                                end if
                                
				kk=0
				do ifg=1,itarget
					if (ilocalstencil(n,i,lm,ifg).gt.0)then
					kk=kk+1
					end if
				end do
				if (kk.eq.itarget)then
				kx=kx+1
				end if

			
		end do
			ielem_admis(i)=kx
			if ((ees.eq.0))then
                            
                            if ((initcond.eq.101).or.(initcond.eq.102))then
                            
                            if (ielem_admis(i).eq.ielem_ifca(i)+1)then
                            ielem_full(i)=1
!                             
                            end if
                            
                            else
                            if (ielem_admis(i).gt.3)then
                            ielem_full(i)=1
                            end if
                            
                            
                            
                            end if
			end if
			
			if ((ees.ge.4))then
                            
                            if ((initcond.eq.101).or.(initcond.eq.405))then
                            
                            if (ielem_admis(i).eq.ielem_ifca(i)+1)then
                            ielem_full(i)=1
!                             
                            end if
                            
                            else
                            if (ielem_admis(i).eq.ielem_ifca(i)+1)then
                            ielem_full(i)=1
                            end if
                            
                            
                            
                            end if
			end if
			if (ees.eq.1)then
                            if (ielem_admis(i).gt.(ielem_ifca(i)+1))then
                            ielem_full(i)=1
                            end if
			end if
			if (ees.eq.2)then
                            if (ielem_admis(i).gt.(ielem_ifca(i)+1))then
                            ielem_full(i)=1
                            end if
			end if
			
			
			
                        
			
! 			
	end do	
	
	if (ees.eq.2)then
	  do i=1,kmaxe
            itrr=0

				
	    
	    if (ielem_admis(i).eq.((ielem_ifca(i)*2)+1))then
	      itrr=1    
	    
	    
	    end if
	    if (itrr.eq.1)then
	      ielem_admis(i)=ielem_ifca(i)+1
	      itrx=1;itry=1
	      do ifg=2,ielem_inumneighbours(i)
			if (mod(ifg,2).eq.0)then
			itrx=itrx+1
				if (ielem_ishape(i).eq.5)then !quad
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
				items(3,ifg)=ilocalstencil(n,i,4,itrx)
				items(4,ifg)=ilocalstencil(n,i,5,itrx)
                                end if
				if (ielem_ishape(i).eq.6)then !tri
				
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
				items(3,ifg)=ilocalstencil(n,i,4,itrx)
! 				
				
				
				end if
                                
                                if (ielem_ishape(i).eq.1)then !hexa
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
				items(3,ifg)=ilocalstencil(n,i,6,itrx)
				items(4,ifg)=ilocalstencil(n,i,7,itrx)
				items(5,ifg)=ilocalstencil(n,i,10,itrx)
				items(6,ifg)=ilocalstencil(n,i,11,itrx)
                                end if
                                 if (ielem_ishape(i).eq.2)then !tetra
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
				items(3,ifg)=ilocalstencil(n,i,6,itrx)
				items(4,ifg)=ilocalstencil(n,i,8,itrx)
                                end if
                                 if (ielem_ishape(i).eq.3)then !pyramidal
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
				items(3,ifg)=ilocalstencil(n,i,4,itrx)
				items(4,ifg)=ilocalstencil(n,i,5,itrx)
				items(5,ifg)=ilocalstencil(n,i,7,itrx)
                                end if
                                if (ielem_ishape(i).eq.4)then !prism
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
				items(3,ifg)=ilocalstencil(n,i,6,itrx)
				items(4,ifg)=ilocalstencil(n,i,7,itrx)
				items(5,ifg)=ilocalstencil(n,i,8,itrx)
                                end if
			else
			itry=itry+1
			    if (ielem_ishape(i).eq.5)then !quad
			    items(1,ifg)=ilocalstencil(n,i,6,itry)
			    items(2,ifg)=ilocalstencil(n,i,7,itry)
			    items(3,ifg)=ilocalstencil(n,i,8,itry)
			    items(4,ifg)=ilocalstencil(n,i,9,itry)
 			    end if
			    if (ielem_ishape(i).eq.6)then !tri
			    
			    items(1,ifg)=ilocalstencil(n,i,5,itry)
			    items(2,ifg)=ilocalstencil(n,i,6,itry)
			    items(3,ifg)=ilocalstencil(n,i,7,itry)
			    end if
			    if (ielem_ishape(i).eq.1)then !hexa
			    items(1,ifg)=ilocalstencil(n,i,4,itry)
			    items(2,ifg)=ilocalstencil(n,i,5,itry)
			    items(3,ifg)=ilocalstencil(n,i,8,itry)
			    items(4,ifg)=ilocalstencil(n,i,9,itry)
			    items(5,ifg)=ilocalstencil(n,i,12,itry)
			    items(6,ifg)=ilocalstencil(n,i,13,itry)
 			    end if
 			    if (ielem_ishape(i).eq.2)then !tetra
			    items(1,ifg)=ilocalstencil(n,i,5,itry)
			    items(2,ifg)=ilocalstencil(n,i,6,itry)
			    items(3,ifg)=ilocalstencil(n,i,9,itry)
			    items(4,ifg)=ilocalstencil(n,i,7,itry)
 			    end if
 			     if (ielem_ishape(i).eq.3)then !pyra
			    items(1,ifg)=ilocalstencil(n,i,10,itry)
			    items(2,ifg)=ilocalstencil(n,i,6,itry)
			    items(3,ifg)=ilocalstencil(n,i,8,itry)
			    items(4,ifg)=ilocalstencil(n,i,9,itry)
			    items(5,ifg)=ilocalstencil(n,i,11,itry)
			    
 			    end if
 			    if (ielem_ishape(i).eq.4)then !prism
				items(1,ifg)=ilocalstencil(n,i,4,itry)
				items(2,ifg)=ilocalstencil(n,i,5,itry)
				items(3,ifg)=ilocalstencil(n,i,9,itry)
				items(4,ifg)=ilocalstencil(n,i,11,itry)
				items(5,ifg)=ilocalstencil(n,i,10,itry)
                                end if
			end if
	    end do
	    if (ielem_ishape(i).eq.5)then
	    ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4,2:ielem_inumneighbours(i))=items(3,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,5,2:ielem_inumneighbours(i))=items(4,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,6:typesten,:)=0
	    
	    
	    end if
	    if (ielem_ishape(i).eq.6)then
	      ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4,2:ielem_inumneighbours(i))=items(3,2:ielem_inumneighbours(i))
	     ilocalstencil(n,i,5:typesten,:)=0
	    
	    end if
	    if (ielem_ishape(i).eq.1)then
	    ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4,2:ielem_inumneighbours(i))=items(3,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,5,2:ielem_inumneighbours(i))=items(4,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,6,2:ielem_inumneighbours(i))=items(5,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,7,2:ielem_inumneighbours(i))=items(6,2:ielem_inumneighbours(i))
	     ilocalstencil(n,i,8:typesten,:)=0
	    
	    end if
	    if (ielem_ishape(i).eq.2)then
	    ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4,2:ielem_inumneighbours(i))=items(3,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,5,2:ielem_inumneighbours(i))=items(4,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,6:typesten,:)=0
	    
	    
	    end if
	      if ((ielem_ishape(i).eq.3).or.(ielem_ishape(i).eq.4))then
	    ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4,2:ielem_inumneighbours(i))=items(3,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,5,2:ielem_inumneighbours(i))=items(4,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,6,2:ielem_inumneighbours(i))=items(5,2:ielem_inumneighbours(i))
	     ilocalstencil(n,i,7:typesten,:)=0
	    
	    end if
	    
	    end if
	  
	  
	  
	  end do
	
	
	
	end if
	
	
	
	
	if (ees.eq.4)then
	  do i=1,kmaxe
            itrr=0

				
	    
	    if (ielem_admis(i).eq.((ielem_ifca(i))+1))then
	      itrr=1    
	    
	    
	    end if
	    if ((itrr.eq.1).and.((ielem_ishape(i).eq.5).or.(ielem_ishape(i).eq.1)))then
              if (ielem_ishape(i).eq.5)then
	      ielem_admis(i)=3
	      end if
	      if (ielem_ishape(i).eq.1)then
	      ielem_admis(i)=4
	      end if
	      
	      itrx=1;itry=1
	      do ifg=2,ielem_inumneighbours(i)
			if (mod(ifg,2).eq.0)then
			itrx=itrx+1
				if (ielem_ishape(i).eq.5)then !quad
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,3,itrx)
                                end if
				
                                
                                if (ielem_ishape(i).eq.1)then !hexa
				items(1,ifg)=ilocalstencil(n,i,2,itrx)
				items(2,ifg)=ilocalstencil(n,i,4,itrx)
				items(3,ifg)=ilocalstencil(n,i,6,itrx)
                                end if
                                 
			else
			itry=itry+1
			    if (ielem_ishape(i).eq.5)then !quad
			    items(1,ifg)=ilocalstencil(n,i,4,itry)
			    items(2,ifg)=ilocalstencil(n,i,5,itry)
 			    end if
			   
			    if (ielem_ishape(i).eq.1)then !hexa
			    items(1,ifg)=ilocalstencil(n,i,3,itry)
			    items(2,ifg)=ilocalstencil(n,i,5,itry)
			    items(3,ifg)=ilocalstencil(n,i,7,itry)
			    
 			    end if
 			   
			end if
	    end do
	    if (ielem_ishape(i).eq.5)then
	    ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	     ilocalstencil(n,i,4:typesten,:)=0
	    
	    end if
	   
	    if (ielem_ishape(i).eq.1)then
	    ilocalstencil(n,i,2,2:ielem_inumneighbours(i))=items(1,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,3,2:ielem_inumneighbours(i))=items(2,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4,2:ielem_inumneighbours(i))=items(3,2:ielem_inumneighbours(i))
	    ilocalstencil(n,i,4:typesten,:)=0
	    end if
	   
	    
	    end if
	  
	  
	  
	  end do
	
	
	
	end if
	
	
	
	
	
	
	
	
	
	
	
	
	
      call mpi_barrier(mpi_comm_world,ierror)
	end subroutine stencils

	
	
	subroutine stencils3(n)
	!> @brief
!> this subroutine is establishing which of the stencils are admissible under different set of rules and which cells can use the weno algorithms
	implicit none
	integer,intent(in)::n
	integer::i,j,ji,k,lm,kmaxn,kk,kmaxe,iaa,kx,l,itrr,itrx,itry
	integer::it1,it2,it3,it4,it5,it6,it7,it8,itx,inx,itf,ifg,kmh,jng,ichg

	
	kmaxe=xmpielrank(n)
	call mpi_barrier(mpi_comm_world,ierror)
	
	do i=1,kmaxe
		
			kx=1
		do lm=1,typesten
				do ifg=2,ielem_inumneighbours(i)
					kx=kx+1
					ilocalstencil(n,i,lm,ifg)=ilocalallelg(n,i,1,kx)
				end do
		end do
		
			ilocalstencil(n,i,:,1)=ilocalallelg(n,i,1,1)
! 			
	end do	
	
	
      call mpi_barrier(mpi_comm_world,ierror)
	end subroutine stencils3


subroutine adapt_criterion
!> @brief
!> this subroutine is establishing a region for which to use a very high-order discretisation and a lower one outside this region
	implicit none
	integer::kmaxe,i,fc
	real::xmin_ad,xmax_ad
	  kmaxe=xmpielrank(n)

          xmax_ad=-tolbig
          xmin_ad=tolbig

	  if (initcond.eq.405)then
	  xmin_ad=-0.2d0
	  xmax_ad=0.2d0

	  end if
	  if (initcond.eq.422)then
	  xmin_ad=0.01
	  xmax_ad=0.15

	  end if

	  do i=1,kmaxe
	      fc=0


	      if (ielem_xxc(i).lt.xmin_ad)then
	      fc=1

	      end if
	       if (ielem_xxc(i).gt.xmax_ad)then

		fc=1
	      end if

 	if (fc.eq.1)then
 		ielem_hybrid(i)=1

 	end if


	if (fc.eq.1)then
		ielem_full(i)=0

	end if


	  end do




	end subroutine adapt_criterion



!---------------------------------------------------------------------------------------------!	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_blocks(n,xmpielrank,xmpinrank,xmpie,xmpin,dinode,imaxn,imaxe,imaxb,xmpinnumber)
!> @brief
!> this subroutine is solely for debugging purposes
implicit none
type(anode_number),allocatable,dimension(:,:),intent(inout)::dinode
integer,allocatable,dimension(:,:,:),intent(inout)::xmpinnumber
integer,intent(in)::n,imaxn,imaxe
integer,allocatable,dimension(:),intent(in)::xmpie,xmpin
integer,allocatable,dimension(:),intent(in)::xmpielrank,xmpinrank
integer,intent(inout)::imaxb
integer::i,j,ji,k,lm,iex,kmaxn,kk,kmaxe
integer::it1,it2,it3,it4,it5,it6,it7,it8,itx,inx
real::x,y,z
kmaxe=xmpielrank(n)

call open_input(n,itt)
do i=1,imaxn
	 read(9,*)inx,x,y,z
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
end do
call close_input(n,itt)
call open_input(n,itt)
do i=1,imaxn
	read(9,*)inx,x,y,z
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
end do
call close_input(n,itt)
call open_input(n,itt)
do i=1,imaxn
	read(9,*)inx,x,y,z
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
end do
call close_input(n,itt)
!-------------------for debugging only -----------------------------------------!


end subroutine write_blocks	




















!---------------------------------------------------------------------------------------------!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!subroutine called for opening grid files!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!and every other related files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!subroutine called for opening grid files!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!and every other related files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function rotate_per(vect_per,code_per,angle_period)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(3),intent(in)::vect_per
real,intent(in)::angle_period
real,dimension(3)::rotate_per
integer,intent(in)::code_per
real::angle_temp


if (code_per.eq.50) then
    angle_temp=-angle_period
    else
    angle_temp=angle_period
end if

rotate_per(1)=vect_per(1)*cos(angle_temp)-vect_per(2)*sin(angle_temp)
rotate_per(2)=vect_per(1)*sin(angle_temp)+vect_per(2)*cos(angle_temp)
rotate_per(3)=vect_per(3)


end function
function rotate_per_1(vect_per,code_per,angle_period)
implicit none
#ifdef gpu
!$omp declare target
#endif
real,dimension(3),intent(in)::vect_per
real,intent(in)::angle_period
real,dimension(3)::rotate_per_1
integer,intent(in)::code_per
real::angle_temp


if (code_per.eq.50) then
    angle_temp=angle_period
    else
    angle_temp=-angle_period
end if

rotate_per_1(1)=vect_per(1)*cos(angle_temp)-vect_per(2)*sin(angle_temp)
rotate_per_1(2)=vect_per(1)*sin(angle_temp)+vect_per(2)*cos(angle_temp)
rotate_per_1(3)=vect_per(3)


end function



subroutine new_arrays(n)
implicit none
integer,intent(in)::n
integer::i,j,k,kmaxe

j=0
k=0

kmaxe=xmpielrank(n)

do i=1,kmaxe
    if (ielem_interior(i).eq.1)then
    
    j=j+1
    else
    
    k=k+1
    
    end if
end do
nof_interior=k
nof_bounded=j

 allocate(el_int(k));el_int(:)=0
 allocate(el_bnd(j));el_int(:)=0

j=0
k=0

do i=1,kmaxe
    if (ielem_interior(i).eq.1)then
    
    j=j+1
    el_bnd(j)=i
    else
    
    k=k+1
    el_int(k)=i
    end if
end do

end subroutine new_arrays





!---------------------------------------------------------------------------------------------!
end module library
