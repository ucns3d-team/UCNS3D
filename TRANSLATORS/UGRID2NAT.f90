PROGRAM UGRID_TRANSLATE
IMPLICIT NONE
	INTEGER::I,J,K,L,N,I1,I2,I3,I4,I5,I6,I7,I8,IOS,IOX,IOY,IMAXE,IMAXB,IMAXN,ICG,KX,NBOUND,DIP
	INTEGER::afnnodesg ! = number of nodes
        INTEGER::afntface  ! = number of boundary triangles
        INTEGER::afnqface  ! = number of boundary quads
        INTEGER::afntet    ! = number of volume TETRA_4 elements
        INTEGER::afnpyr    ! = number of volume PYRA_5 elements
        INTEGER::afnprz    ! = number of volume PENTA_6 elements
        INTEGER::afnhex    ! = number of volume HEXA_8 elements
        INTEGER,ALLOCATABLE,DIMENSION(:)::IBID,IBX,IBXX,ifacetag
	INTEGER,ALLOCATABLE,DIMENSION(:,:)::if2nt,if2nq,ic2nt,ic2np,ic2nz,ic2nh
	real,allocatable,dimension(:)::x,y,z  
                          
                          
                          
OPEN(180,FILE="grid.ugrid",FORM='UNFORMATTED',STATUS='OLD',ACCESS='STREAM',CONVERT="BIG_ENDIAN")
read(180)afnnodesg,afntface,afnqface, afntet, afnpyr, afnprz, afnhex
print*,afnnodesg,afntface, afnqface, afntet, afnpyr, afnprz, afnhex

allocate(x(afnnodesg),y(afnnodesg),z(afnnodesg))
allocate(if2nt(3,afntface),if2nq(4,afnqface),ifacetag(afntface+afnqface),ic2nt(4,afntet),ic2np(5,afnpyr),ic2nz(6,afnprz),ic2nh(8,afnhex))

imaxe=afntet+afnpyr+afnprz+afnhex
imaxb=afntface+afnqface
imaxn=afnnodesg

			 do i=1,afnnodesg;read(180)x(i),y(i),z(i);end do


			 do i=1,afntface;do j=1,3;read(180)if2nt(j,i); end do;end do
print*,"2"
			 do i=1,afnqface;do j=1,4;read(180)if2nq(j,i); end do;end do
print*,"3"
			 do i=1,afntface+afnqface;read(180)ifacetag(i);end do
print*,"4"
			 do i=1,afntet; do j=1,4; read(180)ic2nt(j,i);end do; end do
print*,"5"
			 do i=1,afnpyr; do j=1,5; read(180)ic2np(j,i);end do; end do
print*,"6"
			 do i=1,afnprz; do j=1,6; read(180)ic2nz(j,i);end do; end do
print*,"7"
			 do i=1,afnhex; do j=1,8; read(180)ic2nh(j,i);end do; end do
print*,"8"

close(180)
OPEN(120,FILE="grid.mapbc",STATUS='OLD',FORM='FORMATTED')
READ (120,*) NBOUND

ALLOCATE (IBID(NBOUND),IBX(NBOUND),IBXX(NBOUND))
DO I=1,NBOUND
  READ(120,*)IBID(I),IBX(I)
END DO




DO I=1,NBOUND

	select case(ibx(i))
	
	case (5000,5050)	!farfield
	IBXX(i)=6
	case(6662)	!symmetry	
	IBXX(i)=3
	case(4000)	!wall
	IBXX(i)=4
	case(7031)	!outflow
	IBXX(i)=2
	
	case(7036)	!inflow
	IBXX(i)=1
	
	case(6100)	!periodicity
	IBXX(i)=5
	
	end select
	
end do

	    
	    
	    
	    
	    
			
			
			  !write nodes first
			  OPEN(12,FILE="GRID.vrt",FORM='unformatted',ACTION='WRITE') 
			 do i=1,afnnodesg
			write(12)i,x(i),y(i),z(i)
			! write(12,"(5X,I8,2X,ES21.14,2X,ES21.14,2X,ES21.14)")i,x(i),y(i),z(i)
			 end do
			 close(12)
			 
			 !end nodes writing
			 
			 
			 !write elements now
			 kx=0
			 OPEN(11,FILE="GRID.cel",FORM='unformatted',ACTION='WRITE')
			 
			 !tetra: 1 2 3 3 4 4 4 4
			 do i=1,afntet
			    kx=kx+1
			    !write(11,"(9I10)")kx,ic2nt(1,i),ic2nt(2,i),ic2nt(3,i),ic2nt(3,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i)
				write(11)kx,ic2nt(1,i),ic2nt(2,i),ic2nt(3,i),ic2nt(3,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i),ic2nt(4,i)

			 end do
			 !pyramid: 1 2 3 4 5 5 5 5
			 
			 
			 do i=1,afnpyr
			    kx=kx+1
! 			    write(150,*)kx,ic2np(1,i),ic2np(2,i),ic2np(5,i),ic2np(4,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
			    	!write(11,"(9I10)")kx,ic2np(1,i),ic2np(4,i),ic2np(5,i),ic2np(2,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
				write(11)kx,ic2np(1,i),ic2np(4,i),ic2np(5,i),ic2np(2,i),ic2np(3,i),ic2np(3,i),ic2np(3,i),ic2np(3,i)
			 end do
			 !prism: 1 2 3 3 4 5 6 6
			 
			 
			 
			 
			 
			 do i=1,afnprz
			    kx=kx+1
 			    !write(11,"(9I10)")kx,ic2nz(1,i),ic2nz(2,i),ic2nz(3,i),ic2nz(3,i),ic2nz(4,i),ic2nz(5,i),ic2nz(6,i),ic2nz(6,i)
				write(11)kx,ic2nz(1,i),ic2nz(2,i),ic2nz(3,i),ic2nz(3,i),ic2nz(4,i),ic2nz(5,i),ic2nz(6,i),ic2nz(6,i)
			 end do
			 !hexa: 1 2 3 4 5 6 7 8
			  do i=1,afnhex
			    kx=kx+1
			    !write(11,"(9I10)")kx,ic2nh(1,i),ic2nh(2,i),ic2nh(3,i),ic2nh(4,i),ic2nh(5,i),ic2nh(6,i),ic2nh(7,i),ic2nh(8,i)
				write(11)kx,ic2nh(1,i),ic2nh(2,i),ic2nh(3,i),ic2nh(4,i),ic2nh(5,i),ic2nh(6,i),ic2nh(7,i),ic2nh(8,i)
			 end do
			 close(11)
			 !end writing elements
			 
			 !now write the boundary file
			 OPEN(10,FILE="GRID.bnd",FORM='unformatted',ACTION='WRITE')
			 kx=0
			 !triangle: 1 2 3 3
			  do i=1,afntface
			    kx=kx+1
			    !write(10,"(6I12)")kx,if2nt(1,i),if2nt(2,i),if2nt(3,i),if2nt(3,i),ibxx(ifacetag(kx))
			   write(10)kx,if2nt(1,i),if2nt(2,i),if2nt(3,i),if2nt(3,i),ibxx(ifacetag(kx))
			 end do
			 !quad: 1 2 3 4
			  do i=1,afnqface
			    kx=kx+1
			    !write(10,"(6I12)")kx,if2nq(1,i),if2nq(2,i),if2nq(3,i),if2nq(4,i),ibxx(ifacetag(kx))
				write(10)kx,if2nq(1,i),if2nq(2,i),if2nq(3,i),if2nq(4,i),ibxx(ifacetag(kx))
			 end do
			 


END PROGRAM













