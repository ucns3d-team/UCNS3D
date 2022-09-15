MODULE TRANSFORM
USE MPIINFO
USE DECLARATION


IMPLICIT NONE

 CONTAINS
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 REAL FUNCTION TRIANGLEAREA(N)
!> @brief
!> This function computes the area of a triangle in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
	VVB(1:3)=VEXT(1,1:3)
	VVC(1:3)=VEXT(2,1:3)
	VVD(1:3)=VEXT(3,1:3)
	
	VVA(1,1)=VVB(2);VVA(2,1)=VVC(2);VVA(3,1)=VVD(2)
	VVA(1,2)=VVB(3);VVA(2,2)=VVC(3);VVA(3,2)=VVD(3)
	VVA(1,3)=1.0d0;VVA(2,3)=1.0d0;VVA(3,3)=1.0D0
		VVJACOBSURF(1)=VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3)))
	VVA(1,1)=VVB(3);VVA(2,1)=VVC(3);VVA(3,1)=VVD(3)
	VVA(1,2)=VVB(1);VVA(2,2)=VVC(1);VVA(3,2)=VVD(1)
	VVA(1,3)=1.0d0;VVA(2,3)=1.0d0;VVA(3,3)=1.0D0
		VVJACOBSURF(2)=VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3)))
	VVA(1,1)=VVB(1);VVA(2,1)=VVC(1);VVA(3,1)=VVD(1)
	VVA(1,2)=VVB(2);VVA(2,2)=VVC(2);VVA(3,2)=VVD(2)
	VVA(1,3)=1.0d0; VVA(2,3)=1.0d0;	VVA(3,3)=1.0D0
		VVJACOBSURF(3)=VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3)))
		TRIANGLEAREA=((OO2)*(SQRT((VVJACOBSURF(1)**2)+(VVJACOBSURF(2)**2)+(VVJACOBSURF(3)**2))))

END FUNCTION TRIANGLEAREA



 REAL FUNCTION QUADAREA(N)
 !> @brief
!> This function computes the area of a quadrilateral in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
	
	
	VVE(1:3)=VEXT(4,1:3)-VEXT(2,1:3)
	VVD(1:3)=VEXT(3,1:3)-VEXT(1,1:3)
	
	
	QUADAREA=OO2*sqrt((((VVE(2)*VVD(3))-(VVE(3)*VVD(2)))**2)+(((VVE(3)*VVD(1))-(VVE(1)*VVD(3)))**2)+(((VVE(1)*VVD(2))-(VVE(2)*VVD(1)))**2))
	
	

END FUNCTION QUADAREA


 REAL FUNCTION LINEAREA(N)
 !> @brief
!> This function computes the length of an edge in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
	
	
	
	VVE(1:2)=VEXT(2,1:2)-VEXT(1,1:2)
	
	
	
	linearea=sqrt((vve(1)**2)+(vve(2)**2))
	
	

END FUNCTION LINEAREA

REAL FUNCTION QUADVOLUME(N)
!> @brief
!> This function computes the area of quad in 2D
IMPLICIT NONE
!!$OMP THREADPRIVATE(QUADVOLUME)
INTEGER,INTENT(IN)::N
real::s,t,r,vol
integer::kK,II
	
VVXI(1)=-1.0d0; VVeta(1)=-1.0d0;
VVXI(2)=1.0d0; VVeta(2)=-1.0d0;
VVXI(3)=1.0d0; VVeta(3)=1.0d0; 
VVXI(4)=-1.0d0; VVeta(4)=1.0d0; 



VVnallx(:)=0.0d0;VVnally(:)=0.0


do Kk=1,qp_QUAD
r=QPOINTS(1,Kk)
s=QPOINTS(2,Kk)


do iI=1,4
VVNXI(1)=-(0.25D0)*(1.0D0-s); VVNETA(1)=-(0.25D0)*(1.D0-r);
VVNXI(2)=(0.25D0)*(1.0D0-s); VVNETA(2)=-(0.25D0)*(1.D0+r);
VVNXI(3)=(0.25D0)*(1.0D0+s); VVNETA(3)=(0.25D0)*(1.D0+r); 
VVNXI(4)=-(0.25D0)*(1.0D0+s); VVNETA(4)=(0.25D0)*(1.D0-r); 


VVnallx(ii)=VVnallx(ii)+(VVNXI(ii)*WEQUA3D(Kk))
VVnally(ii)=VVnally(ii)+(VVNETA(ii)*WEQUA3D(Kk))

end do
end do



VVa=0.0d0
VVa1=0.0d0

do iI=1,4

VVnxi(ii)=VVnallx(ii)
VVneta(ii)=VVnally(ii)





    VVa(1,1)=VVa(1,1)+VVnxi(ii)*vext(ii,1); VVa(1,2)=VVa(1,2)+VVnxi(ii)*vext(ii,2)
    VVa(2,1)=VVa(2,1)+VVneta(ii)*vext(ii,1); VVa(2,2)=VVa(2,2)+VVneta(ii)*vext(ii,2)
    
end do

DETA(1)=(VVA(1,1)*VVA(2,2))-(VVA(1,2)*VVA(2,1))

VVA1(1,1)=(VVA(2,2))
VVA1(1,2)=-(VVA(1,2))
VVA1(2,1)=-(VVA(2,1))
VVA1(2,2)=(VVA(1,1))

vol=DETA(1)*4.0d0
VVa1=VVa1/DETA(1)
! Deta(1)=vol

quadvolume=VOL



END FUNCTION QUADVOLUME



REAL FUNCTION TRIANGLEVOLUME(N)
!> @brief
!> This function computes the area of triangle in 2D
IMPLICIT NONE
!!$OMP THREADPRIVATE(TRIANGLEVOLUME)
INTEGER,INTENT(IN)::N
real::s,t,r,vol


VVA(1,1)=VEXT(1,1)-VEXT(3,1)
VVA(1,2)=VEXT(1,2)-VEXT(3,2)
VVA(2,1)=VEXT(2,1)-VEXT(3,1)
VVA(2,2)=VEXT(2,2)-VEXT(3,2)	
VOL=(VVA(1,1)*VVA(2,2))-(VVA(2,1)*VVA(1,2))


VVA1(1,1)=(VVA(2,2))
VVA1(1,2)=-(VVA(1,2))
VVA1(2,1)=-(VVA(2,1))
VVA1(2,2)=(VVA(1,1))
vol=vol*0.50d0
VVA1=VVA1/VOL
Deta(1)=vol
TRIANGLEVOLUME=VOL



END FUNCTION TRIANGLEVOLUME


! ! ! ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION TETRAVOLUME(N)
!> @brief
!> This function computes the volume of a tetrahedrals 

IMPLICIT NONE
!!$OMP THREADPRIVATE(TETRAVOLUME)
INTEGER,INTENT(IN)::N


	VVB(1:3)=VEXT(1,1:3)
	VVC(1:3)=VEXT(2,1:3)
	VVD(1:3)=VEXT(3,1:3)
	VVE(1:3)=vext(4,1:3)
	
	

	VVA(1,1)=VVC(2);VVA(2,1)=VVD(2);VVA(3,1)=VVE(2)
	VVA(1,2)=VVC(3);VVA(2,2)=VVD(3);VVA(3,2)=VVE(3)
	VVA(1,3)=1.0d0;VVA(2,3)=1.0d0;VVA(3,3)=1.0d0
	VVJACOBVOLUME(1)=VVB(1)*(VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3))))
	VVA(1,1)=VVC(1);VVA(2,1)=VVD(1);VVA(3,1)=VVE(1)
	VVA(1,2)=VVC(3);VVA(2,2)=VVD(3);VVA(3,2)=VVE(3)
	VVA(1,3)=1.0d0;VVA(2,3)=1.0d0;VVA(3,3)=1.0d0
	VVJACOBVOLUME(2)=(-VVB(2))*(VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3))))
	VVA(1,1)=VVC(1);VVA(2,1)=VVD(1);VVA(3,1)=VVE(1)
	VVA(1,2)=VVC(2);VVA(2,2)=VVD(2);VVA(3,2)=VVE(2)
	VVA(1,3)=1.0d0;VVA(2,3)=1.0d0;VVA(3,3)=1.0
	VVJACOBVOLUME(3)=VVB(3)*(VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3))))
	VVA(1,1)=VVC(1);VVA(2,1)=VVD(1);VVA(3,1)=VVE(1)
	VVA(1,2)=VVC(2);VVA(2,2)=VVD(2);VVA(3,2)=VVE(2)
	VVA(1,3)=VVC(3);VVA(2,3)=VVD(3);VVA(3,3)=VVE(3)
	VVJACOBVOLUME(4)=(-1.0)*(VVA(1,1)*((VVA(3,3)*VVA(2,2))-(VVA(3,2)*VVA(2,3)))-VVA(2,1)*&
((VVA(3,3)*VVA(1,2))-(VVA(3,2)*VVA(1,3)))+VVA(3,1)*((VVA(2,3)*VVA(1,2))-(VVA(2,2)*VVA(1,3))))

	TETRAVOLUME=abs((0.166666666666666)*(VVJACOBVOLUME(1)+VVJACOBVOLUME(2)+VVJACOBVOLUME(3)+VVJACOBVOLUME(4)))

	
	
END FUNCTION TETRAVOLUME

SUBROUTINE COMPUTEJACOBIANS
!> @brief
!> This function computes the jacobian of a tetrahedral
IMPLICIT NONE


     vva(:,1) = vext(2,:)-vext(1,:)
     vva(:,2) = vext(3,:)-vext(1,:)
     vva(:,3) = vext(4,:)-vext(1,:)
	 Deta(1) =  vva(1,1)*((vva(3,3)*vva(2,2))-(vva(3,2)*vva(2,3))) & 
            - vva(2,1)*((vva(3,3)*vva(1,2))-(vva(3,2)*vva(1,3))) &
			+ vva(3,1)*((vva(2,3)*vva(1,2))-(vva(2,2)*vva(1,3)))


	 vva1(1,1) =  VVA(3,3)*VVA(2,2) - VVA(3,2)*VVA(2,3)
	 vva1(1,2) = -(VVA(3,3)*VVA(1,2) - VVA(3,2)*VVA(1,3))
	 vva1(1,3) =  VVA(2,3)*VVA(1,2) - VVA(2,2)*VVA(1,3)

	 vva1(2,1) = -(VVA(3,3)*VVA(2,1)-VVA(3,1)*VVA(2,3) )   
	 vva1(2,2) = VVA(3,3)*VVA(1,1) -VVA(3,1)*VVA(1,3)
	 vva1(2,3) = -(VVA(2,3)*VVA(1,1)-VVA(2,1)*VVA(1,3)) 

	 vva1(3,1) = VVA(3,2)*VVA(2,1)-VVA(3,1)*VVA(2,2)
	 vva1(3,2) =  -(VVA(3,2)*VVA(1,1)-VVA(3,1)*VVA(1,2))
	 vva1(3,3) =   VVA(2,2)*VVA(1,1)-VVA(2,1)*VVA(1,2)

	 vva1 = vva1/Deta(1)

	
	
END SUBROUTINE COMPUTEJACOBIANs

SUBROUTINE COMPUTeJACOBIANS2
!> @brief
!> This function computes the volume of a triangle
implicit none
VVA(1,1) = VEXT(2,1) - VEXT(1,1); 	 VVA(1,2) = VEXT(3,1) - VEXT(1,1)
	  VVA(2,1) = VEXT(2,2) - VEXT(1,2); 	 VVA(2,2) = VEXT(3,2) - VEXT(1,2)
      DeTA(1) = VVA(1,1)*VVA(2,2) - VVA(1,2)*VVA(2,1)
	  VVA1(1,1) = VEXT(3,2) - Vext(1,2);  VVA1(1,2) = -(VEXT(3,1) - VEXT(1,1))
	  VVA1(2,1) = -(VEXT(2,2) - Vext(1,2));   VVA1(2,2) = VEXT(2,1) - VEXT(1,1)
      VVA1(:,:) = VVA1(:,:)/DETA(1)




END SUBROUTINE COMPUTEJACOBIANS2



REAL FUNCTION hexaVOLUME(N)
!> @brief
!> This function computes the volume of a hexahedral
IMPLICIT NONE
!!$OMP THREADPRIVATE(hexaVOLUME)
INTEGER,INTENT(IN)::N
real::s,t,r,vol
integer::kk,Ii


	
VVXI(1)=-1.0d0; VVeta(1)=-1.0d0; VVzeta(1)=-1.0d0
VVXI(2)=1.0d0; VVeta(2)=-1.0d0; VVzeta(2)=-1.0d0
VVXI(3)=1.0d0; VVeta(3)=1.0d0; VVzeta(3)=-1.0d0
VVXI(4)=-1.0d0; VVeta(4)=1.0d0; VVzeta(4)=-1.0d0
VVXI(5)=-1.0d0; VVeta(5)=-1.0d0; VVzeta(5)=1.0d0
VVXI(6)=1.0d0; VVeta(6)=-1.0d0; VVzeta(6)=1.0d0
VVXI(7)=1.0d0; VVeta(7)=1.0d0; VVzeta(7)=1.0d0
VVXI(8)=-1.0d0; VVeta(8)=1.0d0; VVzeta(8)=1.0d0


VVnallx(:)=0.0d0;VVnally(:)=0.0d0;VVnallz(:)=0.0d0


do kk=1,qp_hexa
r=QPOINTS(1,Kk)
s=QPOINTS(2,Kk)
t=QPOINTS(3,Kk)

do ii=1,8
VVNXI(1)=-(1.0d0/8.0d0)*(1.0-s)*(1.0d0-t); VVNETA(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-T); VVNZETA(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-s);
VVNXI(2)=(1.0d0/8.0d0)*(1.0-s)*(1.0d0-t); VVNETA(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-T); VVNZETA(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-s);
VVNXI(3)=(1.0d0/8.0d0)*(1.0+s)*(1.0d0-t); VVNETA(3)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-T); VVNZETA(3)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+s);
VVNXI(4)=-(1.0d0/8.0d0)*(1.0+s)*(1.0d0-t); VVNETA(4)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-T); VVNZETA(4)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+s);
VVNXI(5)=-(1.0d0/8.0d0)*(1.0-s)*(1.0d0+t); VVNETA(5)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+t); VVNZETA(5)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-s);
VVNXI(6)=(1.0d0/8.0d0)*(1.0-s)*(1.0d0+t); VVNETA(6)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+t); VVNZETA(6)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-s);
VVNXI(7)=(1.0d0/8.0d0)*(1.0+s)*(1.0d0+t); VVNETA(7)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+t); VVNZETA(7)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+s);
VVNXI(8)=-(1.0d0/8.0d0)*(1.0+s)*(1.0d0+t); VVNETA(8)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+t); VVNZETA(8)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+s);

VVnallx(ii)=VVnallx(ii)+(VVNXI(ii)*WEQUA3D(kk))
VVnally(ii)=VVnally(ii)+(VVNETA(ii)*WEQUA3D(kk))
VVnallz(ii)=VVnallz(ii)+(VVNZETA(ii)*WEQUA3D(kk))
end do
end do



VVa=0.0d0
VVa1=0.0d0

do ii=1,8

!  r=xi(ii)
!  s=eta(ii)
!  t=zeta(ii)
VVnxi(ii)=VVnallx(ii)
VVneta(ii)=VVnally(ii)
VVnzeta(ii)=VVnallz(ii)





    VVa(1,1)=VVa(1,1)+VVnxi(ii)*vext(ii,1); VVa(1,2)=VVa(1,2)+VVnxi(ii)*vext(ii,2); VVa(1,3)=VVa(1,3)+VVnxi(ii)*vext(ii,3)
    VVa(2,1)=VVa(2,1)+VVneta(ii)*vext(ii,1); VVa(2,2)=VVa(2,2)+VVneta(ii)*vext(ii,2); VVa(2,3)=VVa(2,3)+VVneta(ii)*vext(ii,3)
    VVa(3,1)=VVa(3,1)+VVnzeta(ii)*vext(ii,1); VVa(3,2)=VVa(3,2)+VVnzeta(ii)*vext(ii,2); VVa(3,3)=VVa(3,3)+VVnzeta(ii)*vext(ii,3)
end do

vol=(VVA(1,1)*VVA(2,2)*VVA(3,3))-(VVA(1,1)*VVA(2,3)*VVA(3,2))-(VVA(1,2)*VVA(2,1)*VVA(3,3))+&
      (VVA(1,2)*VVA(2,3)*VVA(3,1))+(VVA(1,3)*VVA(2,1)*VVA(3,2))-(VVA(1,3)*VVA(2,2)*VVA(3,1))

! vol=vol*8.0d0




 VVA1(1,1)=(VVA(2,2)*VVA(3,3))-(VVA(2,3)*VVA(3,2));VVA1(1,2)=((VVA(1,3)*VVA(3,2))-(VVA(3,3)*VVA(1,2)));VVA1(1,3)=(VVA(1,2)*VVA(2,3))-(VVA(2,2)*VVA(1,3));
VVA1(2,1)=((VVA(2,3)*VVA(3,1))-(VVA(3,3)*VVA(2,1)));VVA1(2,2)=((VVA(1,1)*VVA(3,3))-(VVA(3,1)*VVA(1,3)));VVA1(2,3)=((VVA(1,3)*VVA(2,1))-(VVA(2,3)*VVA(1,1)));
VVA1(3,1)=(VVA(2,1)*VVA(3,2))-(VVA(3,1)*VVA(2,2));VVA1(3,2)=((VVA(1,2)*VVA(3,1))-(VVA(3,2)*VVA(1,1)));VVA1(3,3)=(VVA(1,1)*VVA(2,2))-(VVA(2,1)*VVA(1,2));





deta(1)=(VVA(1,1)*VVA1(1,1))+(VVA(1,2)*VVA1(2,1))+(VVA(1,3)*VVA1(3,1))


VOL=DETA(1)*8.0D0
deta(1)=deta(1)

VVa1=VVa1/DETA(1)




HEXAVOLUME=VOL




END FUNCTION hexaVOLUME



REAL FUNCTION PYRAVOLUME(N)
!> @brief
!> This function computes the volume of a pyramid 
IMPLICIT NONE
!!$OMP THREADPRIVATE(PYRAVOLUME)
INTEGER,INTENT(IN)::N
real::s,t,r,vol
integer::kk,ii

	
vvxi(1)=-1.0d0; vveta(1)=-1.0d0; vvzeta(1)=-1.0d0
vvxi(2)=1.0d0; vveta(2)=-1.0d0; vvzeta(2)=-1.0d0
vvxi(3)=1.0d0; vveta(3)=1.0d0; vvzeta(3)=-1.0d0
vvxi(4)=-1.0d0; vveta(4)=1.0d0; vvzeta(4)=-1.0d0
vvxi(5)=0.0d0; vveta(5)=0.0d0; vvzeta(5)=1.0d0



vvnallx(:)=0.0d0;vvnally(:)=0.0d0;vvnallz(:)=0.0d0


do kk=1,qp_PYRA
r=QPOINTS(1,kk)
s=QPOINTS(2,kk)
t=QPOINTS(3,kk)

do ii=1,5
vvnxi(1)=-(1.0d0/8.0d0)*(1.0-s)*(1.0d0-t); vvneta(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-t); vvnzeta(1)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-s);
vvnxi(2)=(1.0d0/8.0d0)*(1.0-s)*(1.0d0-t); vvneta(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-t); vvnzeta(2)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-s);
vvnxi(3)=(1.0d0/8.0d0)*(1.0+s)*(1.0d0-t); vvneta(3)=(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0-t); vvnzeta(3)=-(1.0d0/8.0d0)*(1.0d0+r)*(1.0d0+s);
vvnxi(4)=-(1.0d0/8.0d0)*(1.0+s)*(1.0d0-t); vvneta(4)=(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0-t); vvnzeta(4)=-(1.0d0/8.0d0)*(1.0d0-r)*(1.0d0+s);
vvnxi(5)=0.0d0; vvneta(5)=0.0d0; vvnzeta(5)=0.5d0;


vvnallx(ii)=vvnallx(ii)+(vvnxi(ii)*WEQUA3D(kk))
vvnally(ii)=vvnally(ii)+(vvneta(ii)*WEQUA3D(kk))
vvnallz(ii)=vvnallz(ii)+(vvnzeta(ii)*WEQUA3D(kk))
end do
end do



vva=0.0d0
vva1=0.0d0

do ii=1,5


vvnxi(ii)=vvnallx(ii)
vvneta(ii)=vvnally(ii)
vvnzeta(ii)=vvnallz(ii)



    vva(1,1)=vva(1,1)+vvnxi(ii)*vext(ii,1); vva(1,2)=vva(1,2)+vvnxi(ii)*vext(ii,2); vva(1,3)=vva(1,3)+vvnxi(ii)*vext(ii,3)
    vva(2,1)=vva(2,1)+vvneta(ii)*vext(ii,1); vva(2,2)=vva(2,2)+vvneta(ii)*vext(ii,2); vva(2,3)=vva(2,3)+vvneta(ii)*vext(ii,3)
    vva(3,1)=vva(3,1)+vvnzeta(ii)*vext(ii,1);vva(3,2)=vva(3,2)+vvnzeta(ii)*vext(ii,2); vva(3,3)=vva(3,3)+vvnzeta(ii)*vext(ii,3)
end do


 VVA1(1,1)=(VVA(2,2)*VVA(3,3))-(VVA(2,3)*VVA(3,2));VVA1(1,2)=((VVA(1,3)*VVA(3,2))-(VVA(3,3)*VVA(1,2)));VVA1(1,3)=(VVA(1,2)*VVA(2,3))-(VVA(2,2)*VVA(1,3));
VVA1(2,1)=((VVA(2,3)*VVA(3,1))-(VVA(3,3)*VVA(2,1)));VVA1(2,2)=((VVA(1,1)*VVA(3,3))-(VVA(3,1)*VVA(1,3)));VVA1(2,3)=((VVA(1,3)*VVA(2,1))-(VVA(2,3)*VVA(1,1)));
VVA1(3,1)=(VVA(2,1)*VVA(3,2))-(VVA(3,1)*VVA(2,2));VVA1(3,2)=((VVA(1,2)*VVA(3,1))-(VVA(3,2)*VVA(1,1)));VVA1(3,3)=(VVA(1,1)*VVA(2,2))-(VVA(2,1)*VVA(1,2));





deta(1)=(VVA(1,1)*VVA1(1,1))+(VVA(1,2)*VVA1(2,1))+(VVA(1,3)*VVA1(3,1))

DETA(1)=DETA(1)

VVa1=VVa1/DETA(1)

VOL=DETA(1)

PYRAVOLUME=VOL










END FUNCTION PYRAVOLUME


REAL FUNCTION PRISMVOLUME(N)
!> @brief
!> This function computes the volume of a prism 
IMPLICIT NONE
! !$OMP THREADPRIVATE(PRISMVOLUME)
INTEGER,INTENT(IN)::N
real::s,t,r,vol
integer::kk,ii




	
vvxi(1)=1.0d0; vveta(1)=0.0d0; vvzeta(1)=-1.0d0
vvxi(2)=0.0d0; vveta(2)=1.0d0; vvzeta(2)=-1.0d0
vvxi(3)=0.0d0; vveta(3)=0.0d0; vvzeta(3)=-1.0d0
vvxi(4)=1.0d0; vveta(4)=0.0d0; vvzeta(4)=0.0d0
vvxi(5)=0.0d0; vveta(5)=1.0d0; vvzeta(5)=0.0d0
vvxi(6)=0.0d0; vveta(6)=0.0d0; vvzeta(6)=0.0d0

 
vvnallx(:)=0.0d0;vvnally(:)=0.0d0;vvnallz(:)=0.0d0
do ii=1,6

do kk=1,qp_prism
r=qpoints(1,kk)
s=qpoints(2,kk)
t=qpoints(3,kk)



vvnxi(1)=0.5d0*(1.0d0-t); vvneta(1)=0.0d0; vvnzeta(1)=-0.5d0*r;
vvnxi(2)=0.0d0; vvneta(2)=0.5d0*(1.0d0-t); vvnzeta(2)=-0.5d0*s;
vvnxi(3)=-0.5d0*(1.0d0-t); vvneta(3)=-0.5d0*(1.0d0-t); vvnzeta(3)=-0.5d0*(1.0-r-s);
vvnxi(4)=0.5d0*(1.0d0+t); vvneta(4)=0.0d0; vvnzeta(4)=0.5d0*r;
vvnxi(5)=0.0d0; vvneta(5)=0.5d0*(1.0d0+t); vvnzeta(5)=0.5d0*s;
vvnxi(6)=-0.5d0*(1.0d0+t); vvneta(6)=-0.5d0*(1.0d0+t); vvnzeta(6)=0.5d0*(1.0-r-s);


vvnallx(ii)=vvnallx(ii)+(vvnxi(ii)*WEQUA3D(kk))
vvnally(ii)=vvnally(ii)+(vvneta(ii)*WEQUA3D(kk))
vvnallz(ii)=vvnallz(ii)+(vvnzeta(ii)*WEQUA3D(kk))
end do
end do



vva=0.0d0
vva1=0.0d0

do ii=1,6

!  r=xi(ii)
!  s=eta(ii)
!  t=zeta(ii)
vvnxi(ii)=vvnallx(ii)
vvneta(ii)=vvnally(ii)
vvnzeta(ii)=vvnallz(ii)





    VVA(1,1)=VVA(1,1)+vvnxi(ii)*vext(ii,1); VVA(1,2)=VVA(1,2)+vvnxi(ii)*vext(ii,2); VVA(1,3)=VVA(1,3)+vvnxi(ii)*vext(ii,3)
    VVA(2,1)=VVA(2,1)+vvneta(ii)*vext(ii,1); VVA(2,2)=VVA(2,2)+vvneta(ii)*vext(ii,2); VVA(2,3)=VVA(2,3)+vvneta(ii)*vext(ii,3)
    VVA(3,1)=VVA(3,1)+vvnzeta(ii)*vext(ii,1); VVA(3,2)=VVA(3,2)+vvnzeta(ii)*vext(ii,2); VVA(3,3)=VVA(3,3)+vvnzeta(ii)*vext(ii,3)
end do


 VVA1(1,1)=(VVA(2,2)*VVA(3,3))-(VVA(2,3)*VVA(3,2));VVA1(1,2)=((VVA(1,3)*VVA(3,2))-(VVA(3,3)*VVA(1,2)));VVA1(1,3)=(VVA(1,2)*VVA(2,3))-(VVA(2,2)*VVA(1,3));
VVA1(2,1)=((VVA(2,3)*VVA(3,1))-(VVA(3,3)*VVA(2,1)));VVA1(2,2)=((VVA(1,1)*VVA(3,3))-(VVA(3,1)*VVA(1,3)));VVA1(2,3)=((VVA(1,3)*VVA(2,1))-(VVA(2,3)*VVA(1,1)));
VVA1(3,1)=(VVA(2,1)*VVA(3,2))-(VVA(3,1)*VVA(2,2));VVA1(3,2)=((VVA(1,2)*VVA(3,1))-(VVA(3,2)*VVA(1,1)));VVA1(3,3)=(VVA(1,1)*VVA(2,2))-(VVA(2,1)*VVA(1,2));





deta(1)=(VVA(1,1)*VVA1(1,1))+(VVA(1,2)*VVA1(2,1))+(VVA(1,3)*VVA1(3,1))

DETA(1)=DETA(1)

VVa1=VVa1/DETA(1)

VOL=DETA(1)




PRISMVOLUME=VOL



END FUNCTION PRISMVOLUME


 FUNCTION CORDINATES3(N,NODES_LIST,N_NODE)
 !> @brief
!> This function computes the centre of 3d element 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,N_NODE
 REAL,DIMENSION(3)::CORDINATES3
real::rnode
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(in)::NODES_LIST
  rnode=n_node
 CORDINATES3(1)=sum(nodes_list(1:n_node,1))/rnode
 CORDINATES3(2)=sum(nodes_list(1:n_node,2))/rnode
 CORDINATES3(3)=sum(nodes_list(1:n_node,3))/rnode


end function CORDINATES3


 FUNCTION distance3(N)
 !> @brief
!> This function computes the distance between two points in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
 REAL::distance3
real::rnode

  distance3=sqrt(((vext(1,1)-vext(2,1))**2)+((vext(1,2)-vext(2,2))**2)+((vext(1,3)-vext(2,3))**2))
  

end function distance3

 FUNCTION distance2(N)
 !> @brief
!> This function computes the distance between two points in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
 REAL::distance2
real::rnode

  distance2=sqrt(((vext(1,1)-vext(2,1))**2)+((vext(1,2)-vext(2,2))**2))
  

end function distance2

 FUNCTION CORDINATES2(N,NODES_LIST,N_NODE)
  !> @brief
!> This function computes the centre of 2d element 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,N_NODE
 REAL,DIMENSION(2)::CORDINATES2
real::rnode
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(in)::NODES_LIST
  rnode=n_node
 CORDINATES2(1)=sum(nodes_list(1:n_node,1))/rnode
 CORDINATES2(2)=sum(nodes_list(1:n_node,2))/rnode
 


end function CORDINATES2



REAL FUNCTION CELL_CENTRE_CORD2(N,CORDS,NODES_LIST,N_NODE)
 !> @brief
!> This function computes the centre of 2d element 
IMPLICIT NONE
INTEGER,INTENT(IN)::N,N_NODE
REAL,ALLOCATABLE,DIMENSION(:),INTENT(OUT)::CORDS
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(in)::NODES_LIST

 cords(1)=sum(nodes_list(1:n_node,1))/n_node
 cords(2)=sum(nodes_list(1:n_node,2))/n_node
 
 CELL_CENTRE_CORD2=cords(1)

end function CELL_CENTRE_CORD2



FUNCTION comp_max_diff(N,NODES_LIST,N_NODE)
!> @brief
!> This function computes the maximum coordinates value given the nodes location
INTEGER,INTENT(IN)::N,N_NODE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(in)::NODES_LIST
REAL,DIMENSION(1:2)::comp_max_diff
INTEGER::Idex
REAL,DIMENSION(1:2)::tempDIFF
tempDiff=0.0d0
comp_max_diff(1:2)=0.0d0


DO Idex=2,N_NODE
    tempDIFF(1)=abs(nodes_list(1,1)-nodes_list(Idex,1))
    tempDIFF(2)=abs(nodes_list(1,2)-nodes_list(Idex,2))
    if (tempDiff(1).gt.comp_max_diff(1)) then
        comp_max_diff(1)=tempDiff(1)
    end if
    if (tempDiff(2).gt.comp_max_diff(2)) then
        comp_max_diff(2)=tempDiff(2)
    end if
END DO

END FUNCTION COMP_MAX_DIFF



FUNCTION comp_min_diff(N,NODES_LIST,N_NODE)
!> @brief
!> This function computes the minimum coordinates value given the nodes location
INTEGER,INTENT(IN)::N,N_NODE
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(in)::NODES_LIST
REAL,DIMENSION(1:2)::comp_min_diff
INTEGER::Idex
REAL,DIMENSION(1:2)::tempDIFF
tempDiff=0.0d0
comp_min_diff(1:2)=0.0d0

DO Idex=2,N_NODE
    tempDIFF(1)=abs(nodes_list(1,1)-nodes_list(Idex,1))
    tempDIFF(2)=abs(nodes_list(1,2)-nodes_list(Idex,2))
    if (tempDiff(1).lt.comp_min_diff(1)) then
        comp_min_diff(1)=tempDiff(1)
    end if
    if (tempDiff(2).lt.comp_min_diff(2)) then
        comp_min_diff(2)=tempDiff(2)
    end if
END DO

END FUNCTION COMP_MIN_DIFF




SUBROUTINE DECOMPOSE3
 !> @brief
!> This function decomposes element into tetrahedrals (counterclockwise numbering)
IMPLICIT NONE

ELEM_LISTD(:,:,:)=zero


SELECT CASE(ELTYPE)
    
    CASE(1)
      ELEM_LISTD(1,1,:)=NODES_LIST(1,:);ELEM_LISTD(1,2,:)=NODES_LIST(6,:);ELEM_LISTD(1,3,:)=NODES_LIST(8,:);ELEM_LISTD(1,4,:)=NODES_LIST(5,:)
      ELEM_LISTD(2,1,:)=NODES_LIST(1,:);ELEM_LISTD(2,2,:)=NODES_LIST(2,:);ELEM_LISTD(2,3,:)=NODES_LIST(8,:);ELEM_LISTD(2,4,:)=NODES_LIST(6,:)
      ELEM_LISTD(3,1,:)=NODES_LIST(2,:);ELEM_LISTD(3,2,:)=NODES_LIST(7,:);ELEM_LISTD(3,3,:)=NODES_LIST(8,:);ELEM_LISTD(3,4,:)=NODES_LIST(6,:)
      ELEM_LISTD(4,1,:)=NODES_LIST(1,:);ELEM_LISTD(4,2,:)=NODES_LIST(8,:);ELEM_LISTD(4,3,:)=NODES_LIST(3,:);ELEM_LISTD(4,4,:)=NODES_LIST(4,:)
      ELEM_LISTD(5,1,:)=NODES_LIST(1,:);ELEM_LISTD(5,2,:)=NODES_LIST(8,:);ELEM_LISTD(5,3,:)=NODES_LIST(2,:);ELEM_LISTD(5,4,:)=NODES_LIST(3,:)
      ELEM_LISTD(6,1,:)=NODES_LIST(2,:);ELEM_LISTD(6,2,:)=NODES_LIST(8,:);ELEM_LISTD(6,3,:)=NODES_LIST(7,:);ELEM_LISTD(6,4,:)=NODES_LIST(3,:)

    CASE(2)
      ELEM_LISTD(1,1,1:3)=NODES_LIST(1,1:3);ELEM_LISTD(1,2,1:3)=NODES_LIST(2,1:3);ELEM_LISTD(1,3,1:3)=NODES_LIST(3,1:3);ELEM_LISTD(1,4,1:3)=NODES_LIST(4,1:3)
    CASE(3)
      ELEM_LISTD(1,1,:)=NODES_LIST(1,:);ELEM_LISTD(1,2,:)=NODES_LIST(2,:);ELEM_LISTD(1,3,:)=NODES_LIST(3,:);ELEM_LISTD(1,4,:)=NODES_LIST(5,:)
      ELEM_LISTD(2,1,:)=NODES_LIST(1,:);ELEM_LISTD(2,2,:)=NODES_LIST(3,:);ELEM_LISTD(2,3,:)=NODES_LIST(4,:);ELEM_LISTD(2,4,:)=NODES_LIST(5,:)
     



    CASE(4)
    
    
    
!      ELEM_LISTD(1,1,:)=NODES_LIST(1,:);ELEM_LISTD(1,2,:)=NODES_LIST(3,:);ELEM_LISTD(1,3,:)=NODES_LIST(2,:);ELEM_LISTD(1,4,:)=NODES_LIST(4,:)
!       ELEM_LISTD(2,1,:)=NODES_LIST(4,:);ELEM_LISTD(2,2,:)=NODES_LIST(3,:);ELEM_LISTD(2,3,:)=NODES_LIST(5,:);ELEM_LISTD(2,4,:)=NODES_LIST(6,:)
!       ELEM_LISTD(3,1,:)=NODES_LIST(3,:);ELEM_LISTD(3,2,:)=NODES_LIST(5,:);ELEM_LISTD(3,3,:)=NODES_LIST(2,:);ELEM_LISTD(3,4,:)=NODES_LIST(4,:)
      ELEM_LISTD(1,1,:)=NODES_LIST(1,:);ELEM_LISTD(1,2,:)=NODES_LIST(2,:);ELEM_LISTD(1,3,:)=NODES_LIST(3,:);ELEM_LISTD(1,4,:)=NODES_LIST(6,:)
      ELEM_LISTD(2,1,:)=NODES_LIST(1,:);ELEM_LISTD(2,2,:)=NODES_LIST(2,:);ELEM_LISTD(2,3,:)=NODES_LIST(6,:);ELEM_LISTD(2,4,:)=NODES_LIST(5,:)
      ELEM_LISTD(3,1,:)=NODES_LIST(1,:);ELEM_LISTD(3,2,:)=NODES_LIST(5,:);ELEM_LISTD(3,3,:)=NODES_LIST(6,:);ELEM_LISTD(3,4,:)=NODES_LIST(4,:)
     
 


  END SELECT





!  SELECT CASE(ELTYPE)
!     
!     CASE(1)
!       ELEM_LISTD(1,1,:)=NODES_LIST(5,:);ELEM_LISTD(1,2,:)=NODES_LIST(7,:);ELEM_LISTD(1,3,:)=NODES_LIST(8,:);ELEM_LISTD(1,4,:)=NODES_LIST(4,:)
!       ELEM_LISTD(2,1,:)=NODES_LIST(5,:);ELEM_LISTD(2,2,:)=NODES_LIST(6,:);ELEM_LISTD(2,3,:)=NODES_LIST(7,:);ELEM_LISTD(2,4,:)=NODES_LIST(4,:)
!       ELEM_LISTD(3,1,:)=NODES_LIST(6,:);ELEM_LISTD(3,2,:)=NODES_LIST(3,:);ELEM_LISTD(3,3,:)=NODES_LIST(7,:);ELEM_LISTD(3,4,:)=NODES_LIST(4,:)
!       ELEM_LISTD(4,1,:)=NODES_LIST(1,:);ELEM_LISTD(4,2,:)=NODES_LIST(5,:);ELEM_LISTD(4,3,:)=NODES_LIST(4,:);ELEM_LISTD(4,4,:)=NODES_LIST(2,:)
!       ELEM_LISTD(5,1,:)=NODES_LIST(5,:);ELEM_LISTD(5,2,:)=NODES_LIST(2,:);ELEM_LISTD(5,3,:)=NODES_LIST(6,:);ELEM_LISTD(5,4,:)=NODES_LIST(4,:)
!       ELEM_LISTD(6,1,:)=NODES_LIST(4,:);ELEM_LISTD(6,2,:)=NODES_LIST(3,:);ELEM_LISTD(6,3,:)=NODES_LIST(2,:);ELEM_LISTD(6,4,:)=NODES_LIST(6,:)
!       
! 
!     CASE(2)
!       ELEM_LISTD(1,1,1:3)=NODES_LIST(1,1:3);ELEM_LISTD(1,2,1:3)=NODES_LIST(3,1:3);ELEM_LISTD(1,3,1:3)=NODES_LIST(2,1:3);ELEM_LISTD(1,4,1:3)=NODES_LIST(4,1:3)
!     CASE(3)
!       ELEM_LISTD(1,1,:)=NODES_LIST(1,:);ELEM_LISTD(1,2,:)=NODES_LIST(4,:);ELEM_LISTD(1,3,:)=NODES_LIST(2,:);ELEM_LISTD(1,4,:)=NODES_LIST(5,:)
!       ELEM_LISTD(2,1,:)=NODES_LIST(2,:);ELEM_LISTD(2,2,:)=NODES_LIST(4,:);ELEM_LISTD(2,3,:)=NODES_LIST(3,:);ELEM_LISTD(2,4,:)=NODES_LIST(5,:)
!      
! 
! 
! 
!     CASE(4)
!     !1324  !4356  !3524
!     
!     
!      
!      ELEM_LISTD(1,1,1:3)=NODES_LIST(1,1:3);ELEM_LISTD(1,2,1:3)=NODES_LIST(3,1:3);ELEM_LISTD(1,3,1:3)=NODES_LIST(2,1:3);ELEM_LISTD(1,4,1:3)=NODES_LIST(4,1:3)
!       ELEM_LISTD(2,1,1:3)=NODES_LIST(4,1:3);ELEM_LISTD(2,2,1:3)=NODES_LIST(3,1:3);ELEM_LISTD(2,3,1:3)=NODES_LIST(5,1:3);ELEM_LISTD(2,4,1:3)=NODES_LIST(6,1:3)
!       ELEM_LISTD(3,1,1:3)=NODES_LIST(3,1:3);ELEM_LISTD(3,2,1:3)=NODES_LIST(5,1:3);ELEM_LISTD(3,3,1:3)=NODES_LIST(2,1:3);ELEM_LISTD(3,4,1:3)=NODES_LIST(4,1:3)
! 
!     
!  
! 
! 
!   END SELECT

 

end SUBROUTINE DECOMPOSE3



subroutine DECOMPOSE2
 !> @brief
!> This function writes decomposed triangle element nodes into ELEM_LISTD from NODES_LIST (counterclockwise numbering)
implicit none
!!$OMP THREADPRIVATE(DECOMPOSE2)

 SELECT CASE(ELTYPE)
    
    CASE(5)
      ELEM_LISTD(1,1,1:2)=NODES_LIST(1,1:2);ELEM_LISTD(1,2,1:2)=NODES_LIST(2,1:2);ELEM_LISTD(1,3,1:2)=NODES_LIST(3,1:2)
      ELEM_LISTD(2,1,1:2)=NODES_LIST(1,1:2);ELEM_LISTD(2,2,1:2)=NODES_LIST(3,1:2);ELEM_LISTD(2,3,1:2)=NODES_LIST(4,1:2)
     

    CASE(6)
      ELEM_LISTD(1,1,1:2)=NODES_LIST(1,1:2);ELEM_LISTD(1,2,1:2)=NODES_LIST(2,1:2);ELEM_LISTD(1,3,1:2)=NODES_LIST(3,1:2)
   
     
 


  END SELECT


end subroutine DECOMPOSE2


SUBROUTINE EDGE_CALCULATOR(N)
 !> @brief
!> This subroutine computes the radius of inscribed sphere or circle
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::I,KMAXE,L
REAL::EDGEL,DIST
KMAXE=XMPIELRANK(N)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(I) 
 IF (DIMENSIONA.EQ.3)THEN
!$OMP DO SCHEDULE (STATIC) 
    DO I=1,KMAXE
    ICONSIDERED=I
    
 	IELEM(N,I)%MINEDGE=(3.0D0*IELEM(N,I)%TOTVOLUME)/(SUM(IELEM(N,I)%SURF(1:IELEM(N,I)%IFCA)))
	
	DO L=1,IELEM(N,I)%IFCA
	FACEX=L
				  CALL coordinates_face_inner(N,Iconsidered,facex)
				  
 				  VEXT(2,1:3)=CORDINATES3(N,NODES_LIST,N_NODE)
				  VEXT(1,1)=IELEM(N,I)%XXC;VEXT(1,2)=IELEM(N,I)%YYC; VEXT(1,3)=IELEM(N,I)%ZZC
				  DIST=DISTANCE3(N)
				  
				
!   				  IELEM(N,I)%MINEDGE=DIST
  				   IELEM(N,I)%MINEDGE=MIN(DIST,IELEM(N,I)%MINEDGE)
! 				  
! !                                      IELEM(N,I)%MINEDGE=dist*2
 	END DO
	
    END DO
!$OMP END DO 

ELSE

!$OMP DO SCHEDULE (STATIC) 
    DO I=1,KMAXE
    ICONSIDERED=I
    
	IELEM(N,I)%MINEDGE=(2.0D0*IELEM(N,I)%TOTVOLUME)/(SUM(IELEM(N,I)%SURF(1:IELEM(N,I)%IFCA)))
	
	DO L=1,IELEM(N,I)%IFCA
	FACEX=L
				  CALL coordinates_face_inner2D(N,Iconsidered,facex)
				  
 				  VEXT(2,1:2)=CORDINATES2(N,NODES_LIST,N_NODE)
				  VEXT(1,1)=IELEM(N,I)%XXC;VEXT(1,2)=IELEM(N,I)%YYC; 
				  DIST=DISTANCE2(N)
				  
				  
				  IELEM(N,I)%MINEDGE=MIN(DIST,IELEM(N,I)%MINEDGE)
				  
	
	END DO
	
    END DO
!$OMP END DO 
END IF
!$OMP END PARALLEL




END SUBROUTINE EDGE_CALCULATOR

SUBROUTINE VOLUME_CALCULATOR3(N)
 !> @brief
!> This subroutine computes the volume of elements
IMPLICIT NONE
INTEGER,INTENT(IN)::N
!$ integer::OMP_IN_PARALLEL,OMP_GET_THREAD_NUM
INTEGER::I,K,KMAXE,jx,JX2
real::DUMV1,DUMV2,dumv3,DUMV5
 KMAXE=XMPIELRANK(N)
DUMV5=ZERO
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(DUMV1,DUMV2,I,JX,jx2,K) 
!$OMP DO 
    DO I=1,KMAXE
    
    VEXT=0.0d0
    NODES_LIST=0.0d0
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
    ELEM_LISTD=0.0d0
     IELEM(N,I)%TOTVOLUME=0.0d0
      
      jx=IELEM(N,I)%NONODES
      

	  do K=1,jx
	    JX2=IELEM(N,I)%NODES(k)
	    NODES_LIST(k,:)=inoder(JX2)%CORD(:)
	    VEXT(K,:)=NODES_LIST(k,:)
	  END DO
	  CALL DECOMPOSE3!(N,ELTYPE,ELEM_DEC)
    
      SELECT CASE(ielem(n,i)%ishape)

      CASE(1)
      CALL QUADRATUREHEXA(N,IGQRULES)
      
      DUMV1=HEXAVOLUME(N)
      DUMV2=0.0d0
       do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)

	  DUMV2=DUMV2+TETRAVOLUME(N)
	
	END DO
	
	IF (ABS(DUMV2-DUMV1).LE.(0.001d0*ABS(DUMV2)))THEN
	IELEM(N,I)%TOTVOLUME=DUMV1
	IELEM(N,I)%MODE=0
	ELSE
	IELEM(N,I)%TOTVOLUME=DUMV2
	IELEM(N,I)%MODE=1
	END IF
	
	IF (DUMV1.LE.ZERO)THEN
	IELEM(N,I)%MODE=1
	IELEM(N,I)%TOTVOLUME=DUMV2
	END IF
	
	ielem(n,i)%mode=1
	IELEM(N,I)%TOTVOLUME=DUMV2

      CASE(2)
      VEXT(1:4,1:3)=ELEM_LISTD(1,1:4,1:3)
	
       
      
	IELEM(N,I)%TOTVOLUME=TETRAVOLUME(N)
	IELEM(N,I)%MODE=1
	
	
	
	
	
	
    
      CASE(3)
      
      CALL QUADRATUREPYRA(N,IGQRULES)
      DUMV1=PYRAVOLUME(N)
      
      
      DUMV2=0.0d0
       do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
            dumv3=TETRAVOLUME(N)

            
	  DUMV2=DUMV2+TETRAVOLUME(N)
    
	END DO

	IELEM(N,I)%TOTVOLUME=DUMV2
	IELEM(N,I)%MODE=1

	
	
	
	
       CASE(4)
      CALL QUADRATUREPRISM(N,IGQRULES)
      
      DUMV1=PRISMVOLUME(N)

      
      
      
      
      
      DUMV2=0.0d0
       do K=1,ELEM_DEC
	VEXT(1:4,1:3)=ELEM_LISTD(k,1:4,1:3)
             dumv3=TETRAVOLUME(N)

	  DUMV2=DUMV2+TETRAVOLUME(N)
    
	END DO
	IF (ABS(DUMV2-DUMV1).LE.(0.001d0*ABS(DUMV2)))THEN
	IELEM(N,I)%TOTVOLUME=DUMV1
	IELEM(N,I)%MODE=0
	ELSE
	IELEM(N,I)%TOTVOLUME=DUMV2
	IELEM(N,I)%MODE=1
	END IF
	IF (DUMV1.LE.ZERO)THEN
	IELEM(N,I)%MODE=1
	END IF
	
 	
 	IELEM(N,I)%MODE=1
        IELEM(N,I)%TOTVOLUME=DUMV2
      END SELECT
   
   
    
   
    END DO
!$OMP END DO 
!$OMP END PARALLEL 
    
!$OMP BARRIER 
!$OMP MASTER
DUMV5=ZERO
DO I=1,KMAXE
    DUMV5=DUMV5+IELEM(N,I)%TOTVOLUME
END DO
 CALL MPI_ALLREDUCE(DUMV5,TOTALVOLUME,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
!$OMP END MASTER
!$OMP BARRIER 

END SUBROUTINE VOLUME_CALCULATOR3



SUBROUTINE VOLUME_CALCULATOR2(N)
 !> @brief
!> This subroutine computes the volume of elements in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
!$ integer::OMP_IN_PARALLEL,OMP_GET_THREAD_NUM
INTEGER::I,K,KMAXE,jx,JX2
real::DUMV1,DUMV2,dumr,DUMV5
 KMAXE=XMPIELRANK(N)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(DUMV1,DUMV2,dumr,I,JX,jx2,K) 
!$OMP DO SCHEDULE (STATIC) 
 DO I=1,KMAXE
    ELTYPE=IELEM(N,I)%ISHAPE
    ELEM_DEC=IELEM(N,I)%VDEC
     IELEM(N,I)%TOTVOLUME=0.0d0
	  do K=1,IELEM(N,I)%NONODES
	    NODES_LIST(k,1:2)=inoder(IELEM(N,I)%NODES(K))%CORD(1:2)
	    vext(k,1:2)=NODES_LIST(k,1:2)
	  END DO
	 call DECOMPOSE2

	    SELECT CASE(ielem(n,i)%ishape)

      CASE(5)

      CALL QUADRATUREQUAD(N,IGQRULES)
      DUMV1=QUADVOLUME(N)
!       do K=1,IELEM(N,I)%IFCA
! 	VEXT(1,1:2)=inoder(IELEM(N,I)%NODES_FACES(K,1))%CORD(1:2)
! 	VEXT(2,1:2)=inoder(IELEM(N,I)%NODES_FACES(K,2))%CORD(1:2)
! 	CALL QUADRATURELINE(N,IGQRULES)
!       END DO
      

      
      DUMV2=0.0d0
      

      
      
       do K=1,ELEM_DEC
	VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
	  DUMV2=DUMV2+TRIANGLEVOLUME(N)
    
	END DO
	IF (ABS(DUMV2-DUMV1).LE.(0.001*DUMV2))THEN
	IELEM(N,I)%TOTVOLUME=DUMV1
	IELEM(N,I)%MODE=0
	ELSE
	IELEM(N,I)%TOTVOLUME=DUMV2
	IELEM(N,I)%MODE=1
	END IF
 	IELEM(N,I)%MODE=1
 	IELEM(N,I)%TOTVOLUME=DUMV2
     
      CASE(6)

      DUMV1=TRIANGLEVOLUME(N)

      

      DUMV2=0.0d0
       do K=1,ELEM_DEC
	VEXT(1:3,1:2)=ELEM_LISTD(k,1:3,1:2)
	  DUMV2=DUMV2+TRIANGLEVOLUME(N)
    
	END DO
	IF (ABS(DUMV2-DUMV1).LE.(0.001d0*DUMV2))THEN
	IELEM(N,I)%TOTVOLUME=DUMV1
	IELEM(N,I)%MODE=0
	ELSE
	IELEM(N,I)%TOTVOLUME=DUMV2
	IELEM(N,I)%MODE=1
	END IF
        IELEM(N,I)%MODE=1
 	IELEM(N,I)%TOTVOLUME=DUMV2
     
    END SELECT
   
    
    END DO
!$OMP END DO 
!$OMP END PARALLEL

 
    
!$OMP BARRIER 
!$OMP MASTER
DUMV5=ZERO
DO I=1,KMAXE
    DUMV5=DUMV5+IELEM(N,I)%TOTVOLUME
END DO
 CALL MPI_ALLREDUCE(DUMV5,TOTALVOLUME,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
!$OMP END MASTER
!$OMP BARRIER 

END SUBROUTINE VOLUME_CALCULATOR2




SUBROUTINE SURFACE_CALCULATOR3(N)
 !> @brief
!> This subroutine computes the surface area of elements in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
!$ integer::OMP_IN_PARALLEL,OMP_GET_THREAD_NUM
INTEGER::I,K,KMAXE,jx,JX2,J,nnd
real::DUMV1,DUMV2
 KMAXE=XMPIELRANK(N)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(DUMV1,DUMV2,I,JX,jx2,K,J,nnd) 
!$OMP DO SCHEDULE (STATIC) 
    DO I=1,KMAXE
    
    
    
    DO J=1,IELEM(N,I)%IFCA
    select case(ielem(n,i)%types_faces(j))
				case (5)
					 
					  NND=4
				      do K=1,nnd
					VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO
					  
					  
					  IELEM(N,I)%surf(J)=QUADarea(N)
					  
				  
				case(6)
					
					NND=3
					do K=1,nnd
					  VEXT(k,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
					END DO
					    
					
 					    IELEM(N,I)%surf(J)=TRIANGLEAREA(N)
 					    
 					    
					  
				end select
    
    
    
    END DO
    
    
    
    DO J=1,IELEM(N,I)%IFCA
                        select case(ielem(n,i)%types_faces(j))
                        case (5)
                                DUMV2=ZERO
				
					 
					
				     
					VEXT(1,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,1))%CORD(1:dims)
					VEXT(2,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,2))%CORD(1:dims)
					VEXT(3,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,3))%CORD(1:dims)
				     
					  
					  
					  DUMV2=DUMV2+TRIANGLEAREA(N)
					  
					  VEXT(1,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,1))%CORD(1:dims)
					VEXT(2,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,3))%CORD(1:dims)
					VEXT(3,1:dims)=inoder(IELEM(N,I)%NODES_FACES(J,4))%CORD(1:dims)
					  DUMV2=DUMV2+TRIANGLEAREA(N)
					  
					  if (abs((DUMV2-IELEM(N,I)%surf(J))/IELEM(N,I)%surf(J))*100.0d0.gt.10.0d0)then
! 					  
					  
					  
					  IELEM(N,I)%surf(J)=dumv2
                                          end if
                                          
                                          
                                          
                                          
				
					
    
    end select
    END DO
    
    
    
    
    
    
    
    END DO
!$OMP END DO 
!$OMP END PARALLEL 
    


END SUBROUTINE SURFACE_CALCULATOR3



SUBROUTINE SURFACE_CALCULATOR2(N)
 !> @brief
!> This subroutine computes the length of edges of elements in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
!$ INTEGER::OMP_IN_PARALLEL,OMP_GET_THREAD_NUM
INTEGER::I,K,KMAXE,JX,JX2,NND,J
REAL::DUMV1,DUMV2,DUMR

KMAXE=XMPIELRANK(N)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(DUMV1,DUMV2,DUMR,I,JX,JX2,K,NND,J) 
!$OMP DO SCHEDULE (STATIC) 

DO I=1,KMAXE
    DO J=1,IELEM(N,I)%IFCA
        NND=2
        
        DO K=1,NND
            VEXT(K,1:DIMS)=INODER(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:DIMS)
        END DO
        
        IELEM(N,I)%SURF(J)=LINEAREA(N)
    END DO
END DO
!$OMP END DO 
!$OMP END PARALLEL 

END SUBROUTINE SURFACE_CALCULATOR2










SUBROUTINE COMPUTE_CENTRE3d(N,Iconsi)
 !> @brief
!> This subroutine computes the cell centre of elements in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,Iconsi
integer::k
    N_NODE=IELEM(N,Iconsi)%NONODES
    do K=1,IELEM(N,Iconsi)%NONODES
      NODES_LIST(k,1:3)=inoder(IELEM(N,Iconsi)%NODES(K))%CORD(1:3)
    END DO
    CORDS=CORDINATES3(N,NODES_LIST,N_NODE)
   
   
   ielem(n,iconsi)%dxx=0.5*(maxval(nodes_list(1:n_node,1))-minval(nodes_list(1:n_node,1)))
   ielem(n,iconsi)%dyy=0.5*(maxval(nodes_list(1:n_node,2))-minval(nodes_list(1:n_node,2)))
   ielem(n,iconsi)%dzz=0.5*(maxval(nodes_list(1:n_node,3))-minval(nodes_list(1:n_node,3)))


END SUBROUTINE

SUBROUTINE COMPUTE_CENTRE3dF(N,Iconsi,facex,IXXFF)
 !> @brief
!> This subroutine retrieve the nodes of faces of elements in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::N,Iconsi,facex,IXXFF
integer::k
    N_NODE=ixxff
  

    do K=1,ixxff

      NODES_LIST(k,1:3)=inoder(IELEM(N,Iconsi)%NODES_FACES(facex,K))%CORD(1:3)
    END DO
   


    CORDS=CORDINATES3(N,NODES_LIST,N_NODE)
   


END SUBROUTINE


subroutine coordinates_face_inner(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieve the nodes of interior faces of elements in 3D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	      select case (ielem(n,iconsidered)%types_faces(facex))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      
	      do K=1,nnd
		  NODES_LIST(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(facex,K))%CORD(1:3)
		  VEXT(K,1:3)=NODES_LIST(k,1:3)
	      END DO
	      
	      
	      N_NODE=NND
	      
	     
end subroutine coordinates_face_inner


subroutine coordinates_face_inner2d(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieves the nodes of edges of elements in 2D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	      nnd=2
	      
	      
	      do K=1,nnd
		  NODES_LIST(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(facex,K))%CORD(1:2)
		  VEXT(K,1:2)=NODES_LIST(k,1:2)
	      END DO
	      
	      
	      N_NODE=NND
	      
	     
end subroutine coordinates_face_inner2d


subroutine coordinates_face_innerx(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieve the nodes of interior faces of elements in 3D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	      select case (ielem(n,iconsidered)%types_faces(facex))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      
	      do K=1,nnd
		  NODES_LIST(k,1:3)=inoder4(IELEM(N,I)%NODES_FACES(facex,K))%CORD(1:3)
		  VEXT(K,1:3)=NODES_LIST(k,1:3)
	      END DO
	      
	      
	      N_NODE=NND
	      
	     
end subroutine coordinates_face_innerx


subroutine coordinates_face_inner2dx(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieves the nodes of edges of elements in 2D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	      nnd=2
	      
	      
	      do K=1,nnd
		  NODES_LIST(k,1:2)=inoder4(IELEM(N,I)%NODES_FACES(facex,K))%CORD(1:2)
		  VEXT(K,1:2)=NODES_LIST(k,1:2)
	      END DO
	      
	      
	      N_NODE=NND
	      
	     
end subroutine coordinates_face_inner2dx


subroutine coordinates_face_PERIOD1(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieve the nodes of periodic faces of elements in 3D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	      select case (ielem(n,iconsidered)%types_faces(facex))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      VEXT(1,1)=IELEM(N,I)%XXC
	      VEXT(1,2)=IELEM(N,I)%YYC
	      VEXT(1,3)=IELEM(N,I)%ZZC
	      
			do K=1,nnd
			  NODES_LIST(k,1:3)=inoder(IELEM(N,I)%NODES_FACES(FACEX,K))%CORD(1:dims)
			END DO
			do K=1,nnd
			IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
			NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
			NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,3)-vext(1,3)).GT.zPER*oo2)THEN
			NODES_LIST(k,3)=NODES_LIST(k,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
			end if
			END DO
	      
	      
	      do K=1,nnd
		  VEXT(K,1:3)=NODES_LIST(k,1:3)
	      END DO
	      
	     N_NODE=NND
end subroutine coordinates_face_PERIOD1


subroutine coordinates_face_PERIOD2d1(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieve the nodes of periodic edges of elements in 2D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	     nnd=2
	      
	      VEXT(1,1)=IELEM(N,I)%XXC
	      VEXT(1,2)=IELEM(N,I)%YYC
	      
	      
			do K=1,nnd
			  NODES_LIST(k,1:2)=inoder(IELEM(N,I)%NODES_FACES(FACEX,K))%CORD(1:dims)
			END DO
			do K=1,nnd
			IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
			NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
			NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
			end if
			END DO
	      
	      
	      do K=1,nnd
		  VEXT(K,1:2)=NODES_LIST(k,1:2)
	      END DO
	      
	     N_NODE=NND
end subroutine coordinates_face_PERIOD2d1



subroutine coordinates_face_PERIOD(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieve the nodes of periodic faces of elements in 3D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	      select case (ielem(n,iconsidered)%types_faces(facex))
	      case(5)
	      nnd=4
	      case(6)
	      nnd=3
	      end select
	      
	      VEXT(1,1)=IELEM(N,I)%XXC
	      VEXT(1,2)=IELEM(N,I)%YYC
	      VEXT(1,3)=IELEM(N,I)%ZZC
	      
			do K=1,nnd
			  NODES_LIST(k,1:3)=inoder4(IELEM(N,I)%NODES_FACES(FACEX,K))%CORD(1:dims)
			END DO
			do K=1,nnd
			IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
			NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
			NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,3)-vext(1,3)).GT.zPER*oo2)THEN
			NODES_LIST(k,3)=NODES_LIST(k,3)+(zPER*SIGN(1.0,vext(1,3)-zPER*oo2))
			end if
			END DO
	      
	      
	      do K=1,nnd
		  VEXT(K,1:3)=NODES_LIST(k,1:3)
	      END DO
	      
	     N_NODE=NND
end subroutine coordinates_face_PERIOD


subroutine coordinates_face_PERIOD2d(n,iconsidered,facex)
 !> @brief
!> This subroutine retrieve the nodes of periodic edges of elements in 2D
IMPLICIT NONE
integer,intent(in)::n,iconsidered,facex
integer::nnd
integer::i,k
i=iconsidered


	     nnd=2
	      
	      VEXT(1,1)=IELEM(N,I)%XXC
	      VEXT(1,2)=IELEM(N,I)%YYC
	      
	      
			do K=1,nnd
			  NODES_LIST(k,1:2)=inoder4(IELEM(N,I)%NODES_FACES(FACEX,K))%CORD(1:dims)
			END DO
			do K=1,nnd
			IF(ABS(NODES_LIST(k,1)-vext(1,1)).GT.XPER*oo2)THEN
			NODES_LIST(k,1)=NODES_LIST(k,1)+(XPER*SIGN(1.0,vext(1,1)-XPER*oo2))
			end if
			IF(ABS(NODES_LIST(k,2)-vext(1,2)).GT.yPER*oo2)THEN
			NODES_LIST(k,2)=NODES_LIST(k,2)+(yPER*SIGN(1.0,vext(1,2)-yPER*oo2))
			end if
			END DO
	      
	      
	      do K=1,nnd
		  VEXT(K,1:2)=NODES_LIST(k,1:2)
	      END DO
	      
	     N_NODE=NND
end subroutine coordinates_face_PERIOD2d



SUBROUTINE COMPUTE_CENTRE2dF(N,Iconsi,facex,IXXFF)
 !> @brief
!> This subroutine retrieves the nodes of the vertices of edges of 2D elements
IMPLICIT NONE
INTEGER,INTENT(IN)::N,Iconsi,facex,IXXFF
integer::k
    N_NODE=ixxff
    
    do K=1,N_NODE
      NODES_LIST(k,1:2)=inoder(IELEM(N,Iconsi)%NODES_FACES(facex,K))%CORD(1:2)
    END DO
    CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
   


END SUBROUTINE


SUBROUTINE COMPUTE_CENTRE2d(N,Iconsi)
 !> @brief
!> This subroutine retrieves the nodes of the vertices of 2D elements
IMPLICIT NONE
INTEGER,INTENT(IN)::N,Iconsi
integer::k
    N_NODE=IELEM(N,Iconsi)%NONODES
    do K=1,IELEM(N,Iconsi)%NONODES
      NODES_LIST(k,1:2)=inoder(IELEM(N,Iconsi)%NODES(K))%CORD(1:2)
    END DO
    CORDS=CORDINATES2(N,NODES_LIST,N_NODE)
   


END SUBROUTINE


SUBROUTINE CENTRE(N)
 !> @brief
!> This subroutine computes the cell centres
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I
KMAXE=XMPIELRANK(N)

!$OMP PARALLEL DEFAULT(SHARED) private(i)
IF (DIMENSIONA.EQ.3)THEN


!$OMP DO SCHEDULE (STATIC) 
DO I=1,KMAXE
    CALL COMPUTE_CENTRE3d(N,I)
    IELEM(N,I)%XXC=CORDS(1);IELEM(N,I)%YYC=CORDS(2);IELEM(N,I)%ZZC=CORDS(3);
    ielem(n,i)%xxc=ielem(n,i)%xxc
    ielem(n,i)%yyc=ielem(n,i)%yyc
    ielem(n,i)%zzc=ielem(n,i)%zzc
END DO
!$OMP END DO


ELSE


!$OMP DO SCHEDULE (STATIC)  
DO I=1,KMAXE
    CALL COMPUTE_CENTRE2d(N,I)
    IELEM(N,I)%XXC=CORDS(1);IELEM(N,I)%YYC=CORDS(2)
END DO
!$OMP END DO



END IF
!$OMP END PARALLEL



END SUBROUTINE CENTRE

SUBROUTINE QUADRATURETRIANG(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for triangle in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
INTEGER::Kk

WEQUA2D=0.0d0
QPOINTS2D=0.0d0


SELECT CASE(IGQRULES)

case(1)
		
	    vvwg(1) = 1.0d0
	    VVR1(1)=1.0d0/3.0d0;	VVR2(1)=1.0d0/3.0d0;	VVR3(1)=1.0d0/3.0d0
case(2)
		vvwg(1)=0.33333333333333333333
  		vvwg(2)=0.33333333333333333333
  		vvwg(3)=0.33333333333333333333

		VVR1(1)=0.666666666666667 ;VVR2(1)=0.166666666666667 ;VVR3(1)=0.166666666666667 
		VVR1(2)=0.166666666666667 ;VVR2(2)=0.666666666666667 ;VVR3(2)=0.166666666666667 
		VVR1(3)=0.166666666666667 ;VVR2(3)=0.166666666666667 ;VVR3(3)=0.666666666666667 


case(3)
		VVR1(1)=0.816847572980440 ;VVR2(1)=0.091576213509780 ;VVR3(1)=0.091576213509780 ;vvwg(1)=0.109951743655333
		VVR1(2)=0.091576213509780 ;VVR2(2)=0.816847572980440 ;VVR3(2)=0.091576213509780 ;vvwg(2)=0.109951743655333
		VVR1(3)=0.091576213509780 ;VVR2(3)=0.091576213509780 ;VVR3(3)=0.816847572980440 ;vvwg(3)=0.109951743655333
		VVR1(4)=0.445948490915964 ;VVR2(4)=0.445948490915964 ;VVR3(4)=0.108103018168071 ;vvwg(4)=0.223381589678000
		VVR1(5)=0.445948490915964 ;VVR2(5)=0.108103018168071 ;VVR3(5)=0.445948490915964 ;vvwg(5)=0.223381589678000
		VVR1(6)=0.108103018168071 ;VVR2(6)=0.445948490915964 ;VVR3(6)=0.445948490915964 ;vvwg(6)=0.223381589678000

case(4)
		
VVR1(1)=0.888871894660413 ;VVR2(1)=0.055564052669793 ;VVR3(1)=0.055564052669793 ;vvwg(1)=0.041955512996649
VVR1(2)=0.055564052669793 ;VVR2(2)=0.888871894660413 ;VVR3(2)=0.055564052669793 ;vvwg(2)=0.041955512996649
VVR1(3)=0.055564052669793 ;VVR2(3)=0.055564052669793 ;VVR3(3)=0.888871894660413 ;vvwg(3)=0.041955512996649
VVR1(4)=0.295533711735893 ;VVR2(4)=0.634210747745723 ;VVR3(4)=0.070255540518384 ;vvwg(4)=0.112098412070887
VVR1(5)=0.295533711735893 ;VVR2(5)=0.070255540518384 ;VVR3(5)=0.634210747745723 ;vvwg(5)=0.112098412070887
VVR1(6)=0.070255540518384 ;VVR2(6)=0.295533711735893 ;VVR3(6)=0.634210747745723 ;vvwg(6)=0.112098412070887
VVR1(7)=0.634210747745723 ;VVR2(7)=0.295533711735893 ;VVR3(7)=0.070255540518384 ;vvwg(7)=0.112098412070887
VVR1(8)=0.634210747745723 ;VVR2(8)=0.070255540518384 ;VVR3(8)=0.295533711735893 ;vvwg(8)=0.112098412070887
VVR1(9)=0.070255540518384 ;VVR2(9)=0.634210747745723 ;VVR3(9)=0.295533711735893 ;vvwg(9)=0.112098412070887
VVR1(10)=0.333333333333333 ;VVR2(10)=0.333333333333333 ;VVR3(10)=0.333333333333333 ;vvwg(10)=0.201542988584730


case(5)

VVR1(1)=0.928258244608533;VVR2(1)= 0.035870877695734 ;VVR3(1)=0.035870877695734 ;vvwg(1)=0.017915455012303
VVR1(2)=0.035870877695734;VVR2(2)= 0.928258244608533 ;VVR3(2)=0.035870877695734 ;vvwg(2)=0.017915455012303
VVR1(3)=0.035870877695734;VVR2(3)= 0.035870877695734 ;VVR3(3)=0.928258244608533 ;vvwg(3)=0.017915455012303
VVR1(4)=0.516541208464066;VVR2(4)= 0.241729395767967 ;VVR3(4)=0.241729395767967 ;vvwg(4)=0.127712195881265
VVR1(5)=0.241729395767967;VVR2(5)= 0.516541208464066 ;VVR3(5)=0.241729395767967 ;vvwg(5)=0.127712195881265
VVR1(6)=0.241729395767967;VVR2(6)= 0.241729395767967 ;VVR3(6)=0.516541208464066 ;vvwg(6)=0.127712195881265
VVR1(7)=0.474308787777079;VVR2(7)= 0.474308787777079 ;VVR3(7)=0.051382424445843 ;vvwg(7)=0.076206062385535
VVR1(8)=0.474308787777079;VVR2(8)= 0.051382424445843 ;VVR3(8)=0.474308787777079 ;vvwg(8)=0.076206062385535
VVR1(9)=0.051382424445843;VVR2(9)= 0.474308787777079 ;VVR3(9)=0.474308787777079 ;vvwg(9)=0.076206062385535
VVR1(10)=0.201503881881800;VVR2(10)= 0.751183631106484 ;VVR3(10)=0.047312487011716 ;vvwg(10)=0.055749810027115
VVR1(11)=0.201503881881800;VVR2(11)= 0.047312487011716 ;VVR3(11)=0.751183631106484 ;vvwg(11)=0.055749810027115
VVR1(12)=0.047312487011716;VVR2(12)= 0.201503881881800 ;VVR3(12)=0.751183631106484 ;vvwg(12)=0.055749810027115
VVR1(13)=0.751183631106484;VVR2(13)= 0.201503881881800 ;VVR3(13)=0.047312487011716 ;vvwg(13)=0.055749810027115
VVR1(14)=0.751183631106484;VVR2(14)= 0.047312487011716 ;VVR3(14)=0.201503881881800 ;vvwg(14)=0.055749810027115
VVR1(15)=0.047312487011716;VVR2(15)= 0.751183631106484 ;VVR3(15)=0.201503881881800 ;vvwg(15)=0.055749810027115








case(6)

VVR1(1)=0.943774095634672   ;VVR2(1)=0.028112952182664  ;VVR3(1)=0.028112952182664 ;vvwg(1)=0.010359374696538
VVR1(2)=0.028112952182664   ;VVR2(2)=0.943774095634672  ;VVR3(2)=0.028112952182664 ;vvwg(2)=0.010359374696538
VVR1(3)=0.028112952182664   ;VVR2(3)=0.028112952182664  ;VVR3(3)=0.943774095634672 ;vvwg(3)=0.010359374696538
VVR1(4)=0.645721803061365   ;VVR2(4)=0.177139098469317  ;VVR3(4)=0.177139098469317 ;vvwg(4)=0.075394884326738
VVR1(5)=0.177139098469317   ;VVR2(5)=0.645721803061365  ;VVR3(5)=0.177139098469317 ;vvwg(5)=0.075394884326738
VVR1(6)=0.177139098469317   ;VVR2(6)=0.177139098469317  ;VVR3(6)=0.645721803061365 ;vvwg(6)=0.075394884326738
VVR1(7)=0.405508595867433   ;VVR2(7)=0.405508595867433  ;VVR3(7)=0.188982808265134 ;vvwg(7)=0.097547802373242
VVR1(8)=0.405508595867433   ;VVR2(8)=0.188982808265134  ;VVR3(8)=0.405508595867433 ;vvwg(8)=0.097547802373242
VVR1(9)=0.188982808265134   ;VVR2(9)=0.405508595867433  ;VVR3(9)=0.405508595867433 ;vvwg(9)=0.097547802373242
VVR1(10)=0.148565812270887 ;VVR2(10)=0.817900980028499 ;VVR3(10)=0.033533207700614 ;vvwg(10)=0.028969269372473
VVR1(11)=0.148565812270887 ;VVR2(11)=0.033533207700614 ;VVR3(11)=0.817900980028499 ;vvwg(11)=0.028969269372473
VVR1(12)=0.033533207700614 ;VVR2(12)=0.148565812270887 ;VVR3(12)=0.817900980028499 ;vvwg(12)=0.028969269372473
VVR1(13)=0.817900980028499 ;VVR2(13)=0.148565812270887 ;VVR3(13)=0.033533207700614 ;vvwg(13)=0.028969269372473
VVR1(14)=0.817900980028499 ;VVR2(14)=0.033533207700614 ;VVR3(14)=0.148565812270887 ;vvwg(14)=0.028969269372473
VVR1(15)=0.033533207700614 ;VVR2(15)=0.817900980028499 ;VVR3(15)=0.148565812270887 ;vvwg(15)=0.028969269372473
VVR1(16)=0.357196298615681 ;VVR2(16)=0.604978911775132 ;VVR3(16)=0.037824789609186 ;vvwg(16)=0.046046366595935
VVR1(17)=0.357196298615681 ;VVR2(17)=0.037824789609186 ;VVR3(17)=0.604978911775132 ;vvwg(17)=0.046046366595935
VVR1(18)=0.037824789609186 ;VVR2(18)=0.357196298615681 ;VVR3(18)=0.604978911775132 ;vvwg(18)=0.046046366595935
VVR1(19)=0.604978911775132 ;VVR2(19)=0.357196298615681 ;VVR3(19)=0.037824789609186 ;vvwg(19)=0.046046366595935
VVR1(20)=0.604978911775132 ;VVR2(20)=0.037824789609186 ;VVR3(20)=0.357196298615681 ;vvwg(20)=0.046046366595935
VVR1(21)=0.037824789609186 ;VVR2(21)=0.604978911775132 ;VVR3(21)=0.357196298615681 ;vvwg(21)=0.046046366595935




case(7,8,9)

VVR1(1)=0.957657154441070
VVR1(2)=0.021171422779465
VVR1(3)=0.021171422779465
VVR1(4)=0.798831205208225
VVR1(5)=0.100584397395888
VVR1(6)=0.100584397395888
VVR1(7)=0.457923384576135
VVR1(8)=0.271038307711932
VVR1(9)=0.271038307711932
VVR1(10)=0.440191258403832
VVR1(11)=0.440191258403832
VVR1(12)=0.119617483192335
VVR1(13)=0.101763679498021
VVR1(14)=0.101763679498021
VVR1(15)=0.018256679074748
VVR1(16)=0.879979641427232
VVR1(17)=0.879979641427232
VVR1(18)=0.018256679074748
VVR1(19)=0.394033271669987
VVR1(20)=0.394033271669987
VVR1(21)=0.023404705466341
VVR1(22)=0.582562022863673
VVR1(23)=0.582562022863673
VVR1(24)=0.023404705466341
VVR1(25)=0.226245530909229
VVR1(26)=0.226245530909229
VVR1(27)=0.022223854547989
VVR1(28)=0.751530614542782
VVR1(29)=0.751530614542782
VVR1(30)=0.022223854547989
VVR1(31)=0.635737183263105
VVR1(32)=0.635737183263105
VVR1(33)=0.115183589115563
VVR1(34)=0.249079227621332
VVR1(35)=0.249079227621332
VVR1(36)=0.115183589115563



VVR2(1)=0.021171422779465
VVR2(2)=0.957657154441070
VVR2(3)=0.021171422779465
VVR2(4)=0.100584397395888
VVR2(5)=0.798831205208225
VVR2(6)=0.100584397395888
VVR2(7)=0.271038307711932
VVR2(8)=0.457923384576135
VVR2(9)=0.271038307711932
VVR2(10)=0.440191258403832
VVR2(11)=0.119617483192335
VVR2(12)=0.440191258403832
VVR2(13)=0.879979641427232
VVR2(14)=0.018256679074748
VVR2(15)=0.101763679498021
VVR2(16)=0.101763679498021
VVR2(17)=0.018256679074748
VVR2(18)=0.879979641427232
VVR2(19)=0.582562022863673
VVR2(20)=0.023404705466341
VVR2(21)=0.394033271669987
VVR2(22)=0.394033271669987
VVR2(23)=0.023404705466341
VVR2(24)=0.582562022863673
VVR2(25)=0.751530614542782
VVR2(26)=0.022223854547989
VVR2(27)=0.226245530909229
VVR2(28)=0.226245530909229
VVR2(29)=0.022223854547989
VVR2(30)=0.751530614542782
VVR2(31)=0.249079227621332
VVR2(32)=0.115183589115563
VVR2(33)=0.635737183263105
VVR2(34)=0.635737183263105
VVR2(35)=0.115183589115563
VVR2(36)=0.249079227621332




VVR3(1)=0.021171422779465
VVR3(2)=0.021171422779465
VVR3(3)=0.957657154441070
VVR3(4)=0.100584397395888
VVR3(5)=0.100584397395888
VVR3(6)=0.798831205208225
VVR3(7)=0.271038307711932
VVR3(8)=0.271038307711932
VVR3(9)=0.457923384576135
VVR3(10)=0.119617483192335
VVR3(11)=0.440191258403832
VVR3(12)=0.440191258403832
VVR3(13)=0.018256679074748
VVR3(14)=0.879979641427232
VVR3(15)=0.879979641427232
VVR3(16)=0.018256679074748
VVR3(17)=0.101763679498021
VVR3(18)=0.101763679498021
VVR3(19)=0.023404705466341
VVR3(20)=0.582562022863673
VVR3(21)=0.582562022863673
VVR3(22)=0.023404705466341
VVR3(23)=0.394033271669987
VVR3(24)=0.394033271669987
VVR3(25)=0.022223854547989
VVR3(26)=0.751530614542782
VVR3(27)=0.751530614542782
VVR3(28)=0.022223854547989
VVR3(29)=0.226245530909229
VVR3(30)=0.226245530909229
VVR3(31)=0.115183589115563
VVR3(32)=0.249079227621332
VVR3(33)=0.249079227621332
VVR3(34)=0.115183589115563
VVR3(35)=0.635737183263105
VVR3(36)=0.635737183263105



vvwg(1)=0.005639123786910
vvwg(2)=0.005639123786910
vvwg(3)=0.005639123786910
vvwg(4)=0.027148968192278
vvwg(5)=0.027148968192278
vvwg(6)=0.027148968192278
vvwg(7)=0.063100912533359
vvwg(8)=0.063100912533359
vvwg(9)=0.063100912533359
vvwg(10)=0.051752795679899
vvwg(11)=0.051752795679899
vvwg(12)=0.051752795679899
vvwg(13)=0.009866753574646
vvwg(14)=0.009866753574646
vvwg(15)=0.009866753574646
vvwg(16)=0.009866753574646
vvwg(17)=0.009866753574646
vvwg(18)=0.009866753574646
vvwg(19)=0.022008204800147
vvwg(20)=0.022008204800147
vvwg(21)=0.022008204800147
vvwg(22)=0.022008204800147
vvwg(23)=0.022008204800147
vvwg(24)=0.022008204800147
vvwg(25)=0.016644570076736
vvwg(26)=0.016644570076736
vvwg(27)=0.016644570076736
vvwg(28)=0.016644570076736
vvwg(29)=0.016644570076736
vvwg(30)=0.016644570076736
vvwg(31)=0.044326238118914
vvwg(32)=0.044326238118914
vvwg(33)=0.044326238118914
vvwg(34)=0.044326238118914
vvwg(35)=0.044326238118914
vvwg(36)=0.044326238118914















		
END select



		DO Kk=1,qp_triangle
			WEQUA2D(kK)=vvwg(Kk)
			QPOINTS2D(:,kk)=(VVR1(kk)*VEXT(1,:))+(VVR2(kk)*VEXT(2,:))+(VVR3(kk)*VEXT(3,:))
			
		END DO




END SUBROUTINE QUADRATURETRIANG

SUBROUTINE QUADRATURETRIANGLE(N,IGQRULES)
!> @brief
!> This subroutine computes the quadrature points and weights for triangle in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
INTEGER::Kk

WEQUA3D=0.0d0
QPOINTS=0.0d0

select case(IGQRULES)


case(1)
		
	    vvwg(1) = 1.0d0
	    VVR1(1)=1.0d0/3.0d0;	VVR2(1)=1.0d0/3.0d0;	VVR3(1)=1.0d0/3.0d0

case(2)
		vvwg(1)=0.33333333333333333333
  		vvwg(2)=0.33333333333333333333
  		vvwg(3)=0.33333333333333333333

		VVR1(1)=0.666666666666667 ;VVR2(1)=0.166666666666667 ;VVR3(1)=0.166666666666667 
		VVR1(2)=0.166666666666667 ;VVR2(2)=0.666666666666667 ;VVR3(2)=0.166666666666667 
		VVR1(3)=0.166666666666667 ;VVR2(3)=0.166666666666667 ;VVR3(3)=0.666666666666667 

case(3)
		VVR1(1)=0.816847572980440 ;VVR2(1)=0.091576213509780 ;VVR3(1)=0.091576213509780 ;vvwg(1)=0.109951743655333
		VVR1(2)=0.091576213509780 ;VVR2(2)=0.816847572980440 ;VVR3(2)=0.091576213509780 ;vvwg(2)=0.109951743655333
		VVR1(3)=0.091576213509780 ;VVR2(3)=0.091576213509780 ;VVR3(3)=0.816847572980440 ;vvwg(3)=0.109951743655333
		VVR1(4)=0.445948490915964 ;VVR2(4)=0.445948490915964 ;VVR3(4)=0.108103018168071 ;vvwg(4)=0.223381589678000
		VVR1(5)=0.445948490915964 ;VVR2(5)=0.108103018168071 ;VVR3(5)=0.445948490915964 ;vvwg(5)=0.223381589678000
		VVR1(6)=0.108103018168071 ;VVR2(6)=0.445948490915964 ;VVR3(6)=0.445948490915964 ;vvwg(6)=0.223381589678000

case(4)
		
VVR1(1)=0.888871894660413 ;VVR2(1)=0.055564052669793 ;VVR3(1)=0.055564052669793 ;vvwg(1)=0.041955512996649
VVR1(2)=0.055564052669793 ;VVR2(2)=0.888871894660413 ;VVR3(2)=0.055564052669793 ;vvwg(2)=0.041955512996649
VVR1(3)=0.055564052669793 ;VVR2(3)=0.055564052669793 ;VVR3(3)=0.888871894660413 ;vvwg(3)=0.041955512996649
VVR1(4)=0.295533711735893 ;VVR2(4)=0.634210747745723 ;VVR3(4)=0.070255540518384 ;vvwg(4)=0.112098412070887
VVR1(5)=0.295533711735893 ;VVR2(5)=0.070255540518384 ;VVR3(5)=0.634210747745723 ;vvwg(5)=0.112098412070887
VVR1(6)=0.070255540518384 ;VVR2(6)=0.295533711735893 ;VVR3(6)=0.634210747745723 ;vvwg(6)=0.112098412070887
VVR1(7)=0.634210747745723 ;VVR2(7)=0.295533711735893 ;VVR3(7)=0.070255540518384 ;vvwg(7)=0.112098412070887
VVR1(8)=0.634210747745723 ;VVR2(8)=0.070255540518384 ;VVR3(8)=0.295533711735893 ;vvwg(8)=0.112098412070887
VVR1(9)=0.070255540518384 ;VVR2(9)=0.634210747745723 ;VVR3(9)=0.295533711735893 ;vvwg(9)=0.112098412070887
VVR1(10)=0.333333333333333 ;VVR2(10)=0.333333333333333 ;VVR3(10)=0.333333333333333 ;vvwg(10)=0.201542988584730


case(5)

VVR1(1)=0.928258244608533;VVR2(1)= 0.035870877695734 ;VVR3(1)=0.035870877695734 ;vvwg(1)=0.017915455012303
VVR1(2)=0.035870877695734;VVR2(2)= 0.928258244608533 ;VVR3(2)=0.035870877695734 ;vvwg(2)=0.017915455012303
VVR1(3)=0.035870877695734;VVR2(3)= 0.035870877695734 ;VVR3(3)=0.928258244608533 ;vvwg(3)=0.017915455012303
VVR1(4)=0.516541208464066;VVR2(4)= 0.241729395767967 ;VVR3(4)=0.241729395767967 ;vvwg(4)=0.127712195881265
VVR1(5)=0.241729395767967;VVR2(5)= 0.516541208464066 ;VVR3(5)=0.241729395767967 ;vvwg(5)=0.127712195881265
VVR1(6)=0.241729395767967;VVR2(6)= 0.241729395767967 ;VVR3(6)=0.516541208464066 ;vvwg(6)=0.127712195881265
VVR1(7)=0.474308787777079;VVR2(7)= 0.474308787777079 ;VVR3(7)=0.051382424445843 ;vvwg(7)=0.076206062385535
VVR1(8)=0.474308787777079;VVR2(8)= 0.051382424445843 ;VVR3(8)=0.474308787777079 ;vvwg(8)=0.076206062385535
VVR1(9)=0.051382424445843;VVR2(9)= 0.474308787777079 ;VVR3(9)=0.474308787777079 ;vvwg(9)=0.076206062385535
VVR1(10)=0.201503881881800;VVR2(10)= 0.751183631106484 ;VVR3(10)=0.047312487011716 ;vvwg(10)=0.055749810027115
VVR1(11)=0.201503881881800;VVR2(11)= 0.047312487011716 ;VVR3(11)=0.751183631106484 ;vvwg(11)=0.055749810027115
VVR1(12)=0.047312487011716;VVR2(12)= 0.201503881881800 ;VVR3(12)=0.751183631106484 ;vvwg(12)=0.055749810027115
VVR1(13)=0.751183631106484;VVR2(13)= 0.201503881881800 ;VVR3(13)=0.047312487011716 ;vvwg(13)=0.055749810027115
VVR1(14)=0.751183631106484;VVR2(14)= 0.047312487011716 ;VVR3(14)=0.201503881881800 ;vvwg(14)=0.055749810027115
VVR1(15)=0.047312487011716;VVR2(15)= 0.751183631106484 ;VVR3(15)=0.201503881881800 ;vvwg(15)=0.055749810027115








case(6)

VVR1(1)=0.943774095634672 ;VVR2(1)=0.028112952182664 ;VVR3(1)=0.028112952182664 ;vvwg(1)=0.010359374696538
VVR1(2)=0.028112952182664 ;VVR2(2)=0.943774095634672 ;VVR3(2)=0.028112952182664 ;vvwg(2)=0.010359374696538
VVR1(3)=0.028112952182664 ;VVR2(3)=0.028112952182664 ;VVR3(3)=0.943774095634672 ;vvwg(3)=0.010359374696538
VVR1(4)=0.645721803061365 ;VVR2(4)=0.177139098469317 ;VVR3(4)=0.177139098469317 ;vvwg(4)=0.075394884326738
VVR1(5)=0.177139098469317 ;VVR2(5)=0.645721803061365 ;VVR3(5)=0.177139098469317 ;vvwg(5)=0.075394884326738
VVR1(6)=0.177139098469317 ;VVR2(6)=0.177139098469317 ;VVR3(6)=0.645721803061365 ;vvwg(6)=0.075394884326738
VVR1(7)=0.405508595867433 ;VVR2(7)=0.405508595867433 ;VVR3(7)=0.188982808265134 ;vvwg(7)=0.097547802373242
VVR1(8)=0.405508595867433 ;VVR2(8)=0.188982808265134 ;VVR3(8)=0.405508595867433 ;vvwg(8)=0.097547802373242
VVR1(9)=0.188982808265134 ;VVR2(9)=0.405508595867433 ;VVR3(9)=0.405508595867433 ;vvwg(9)=0.097547802373242
VVR1(10)=0.148565812270887 ;VVR2(10)=0.817900980028499 ;VVR3(10)=0.033533207700614 ;vvwg(10)=0.028969269372473
VVR1(11)=0.148565812270887 ;VVR2(11)=0.033533207700614 ;VVR3(11)=0.817900980028499 ;vvwg(11)=0.028969269372473
VVR1(12)=0.033533207700614 ;VVR2(12)=0.148565812270887 ;VVR3(12)=0.817900980028499 ;vvwg(12)=0.028969269372473
VVR1(13)=0.817900980028499 ;VVR2(13)=0.148565812270887 ;VVR3(13)=0.033533207700614 ;vvwg(13)=0.028969269372473
VVR1(14)=0.817900980028499 ;VVR2(14)=0.033533207700614 ;VVR3(14)=0.148565812270887 ;vvwg(14)=0.028969269372473
VVR1(15)=0.033533207700614 ;VVR2(15)=0.817900980028499 ;VVR3(15)=0.148565812270887 ;vvwg(15)=0.028969269372473
VVR1(16)=0.357196298615681 ;VVR2(16)=0.604978911775132 ;VVR3(16)=0.037824789609186 ;vvwg(16)=0.046046366595935
VVR1(17)=0.357196298615681 ;VVR2(17)=0.037824789609186 ;VVR3(17)=0.604978911775132 ;vvwg(17)=0.046046366595935
VVR1(18)=0.037824789609186 ;VVR2(18)=0.357196298615681 ;VVR3(18)=0.604978911775132 ;vvwg(18)=0.046046366595935
VVR1(19)=0.604978911775132 ;VVR2(19)=0.357196298615681 ;VVR3(19)=0.037824789609186 ;vvwg(19)=0.046046366595935
VVR1(20)=0.604978911775132 ;VVR2(20)=0.037824789609186 ;VVR3(20)=0.357196298615681 ;vvwg(20)=0.046046366595935
VVR1(21)=0.037824789609186 ;VVR2(21)=0.604978911775132 ;VVR3(21)=0.357196298615681 ;vvwg(21)=0.046046366595935

case(7,8,9)

VVR1(1)=0.957657154441070
VVR1(2)=0.021171422779465
VVR1(3)=0.021171422779465
VVR1(4)=0.798831205208225
VVR1(5)=0.100584397395888
VVR1(6)=0.100584397395888
VVR1(7)=0.457923384576135
VVR1(8)=0.271038307711932
VVR1(9)=0.271038307711932
VVR1(10)=0.440191258403832
VVR1(11)=0.440191258403832
VVR1(12)=0.119617483192335
VVR1(13)=0.101763679498021
VVR1(14)=0.101763679498021
VVR1(15)=0.018256679074748
VVR1(16)=0.879979641427232
VVR1(17)=0.879979641427232
VVR1(18)=0.018256679074748
VVR1(19)=0.394033271669987
VVR1(20)=0.394033271669987
VVR1(21)=0.023404705466341
VVR1(22)=0.582562022863673
VVR1(23)=0.582562022863673
VVR1(24)=0.023404705466341
VVR1(25)=0.226245530909229
VVR1(26)=0.226245530909229
VVR1(27)=0.022223854547989
VVR1(28)=0.751530614542782
VVR1(29)=0.751530614542782
VVR1(30)=0.022223854547989
VVR1(31)=0.635737183263105
VVR1(32)=0.635737183263105
VVR1(33)=0.115183589115563
VVR1(34)=0.249079227621332
VVR1(35)=0.249079227621332
VVR1(36)=0.115183589115563



VVR2(1)=0.021171422779465
VVR2(2)=0.957657154441070
VVR2(3)=0.021171422779465
VVR2(4)=0.100584397395888
VVR2(5)=0.798831205208225
VVR2(6)=0.100584397395888
VVR2(7)=0.271038307711932
VVR2(8)=0.457923384576135
VVR2(9)=0.271038307711932
VVR2(10)=0.440191258403832
VVR2(11)=0.119617483192335
VVR2(12)=0.440191258403832
VVR2(13)=0.879979641427232
VVR2(14)=0.018256679074748
VVR2(15)=0.101763679498021
VVR2(16)=0.101763679498021
VVR2(17)=0.018256679074748
VVR2(18)=0.879979641427232
VVR2(19)=0.582562022863673
VVR2(20)=0.023404705466341
VVR2(21)=0.394033271669987
VVR2(22)=0.394033271669987
VVR2(23)=0.023404705466341
VVR2(24)=0.582562022863673
VVR2(25)=0.751530614542782
VVR2(26)=0.022223854547989
VVR2(27)=0.226245530909229
VVR2(28)=0.226245530909229
VVR2(29)=0.022223854547989
VVR2(30)=0.751530614542782
VVR2(31)=0.249079227621332
VVR2(32)=0.115183589115563
VVR2(33)=0.635737183263105
VVR2(34)=0.635737183263105
VVR2(35)=0.115183589115563
VVR2(36)=0.249079227621332




VVR3(1)=0.021171422779465
VVR3(2)=0.021171422779465
VVR3(3)=0.957657154441070
VVR3(4)=0.100584397395888
VVR3(5)=0.100584397395888
VVR3(6)=0.798831205208225
VVR3(7)=0.271038307711932
VVR3(8)=0.271038307711932
VVR3(9)=0.457923384576135
VVR3(10)=0.119617483192335
VVR3(11)=0.440191258403832
VVR3(12)=0.440191258403832
VVR3(13)=0.018256679074748
VVR3(14)=0.879979641427232
VVR3(15)=0.879979641427232
VVR3(16)=0.018256679074748
VVR3(17)=0.101763679498021
VVR3(18)=0.101763679498021
VVR3(19)=0.023404705466341
VVR3(20)=0.582562022863673
VVR3(21)=0.582562022863673
VVR3(22)=0.023404705466341
VVR3(23)=0.394033271669987
VVR3(24)=0.394033271669987
VVR3(25)=0.022223854547989
VVR3(26)=0.751530614542782
VVR3(27)=0.751530614542782
VVR3(28)=0.022223854547989
VVR3(29)=0.226245530909229
VVR3(30)=0.226245530909229
VVR3(31)=0.115183589115563
VVR3(32)=0.249079227621332
VVR3(33)=0.249079227621332
VVR3(34)=0.115183589115563
VVR3(35)=0.635737183263105
VVR3(36)=0.635737183263105



vvwg(1)=0.005639123786910
vvwg(2)=0.005639123786910
vvwg(3)=0.005639123786910
vvwg(4)=0.027148968192278
vvwg(5)=0.027148968192278
vvwg(6)=0.027148968192278
vvwg(7)=0.063100912533359
vvwg(8)=0.063100912533359
vvwg(9)=0.063100912533359
vvwg(10)=0.051752795679899
vvwg(11)=0.051752795679899
vvwg(12)=0.051752795679899
vvwg(13)=0.009866753574646
vvwg(14)=0.009866753574646
vvwg(15)=0.009866753574646
vvwg(16)=0.009866753574646
vvwg(17)=0.009866753574646
vvwg(18)=0.009866753574646
vvwg(19)=0.022008204800147
vvwg(20)=0.022008204800147
vvwg(21)=0.022008204800147
vvwg(22)=0.022008204800147
vvwg(23)=0.022008204800147
vvwg(24)=0.022008204800147
vvwg(25)=0.016644570076736
vvwg(26)=0.016644570076736
vvwg(27)=0.016644570076736
vvwg(28)=0.016644570076736
vvwg(29)=0.016644570076736
vvwg(30)=0.016644570076736
vvwg(31)=0.044326238118914
vvwg(32)=0.044326238118914
vvwg(33)=0.044326238118914
vvwg(34)=0.044326238118914
vvwg(35)=0.044326238118914
vvwg(36)=0.044326238118914
		
		
END select



		DO Kk=1,qp_triangle
			WEQUA3D(kk)=vvwg(kk)
			QPOINTS(:,kk)=(VVR1(kk)*VEXT(1,:))+(VVR2(kk)*VEXT(2,:))+(VVR3(kk)*VEXT(3,:))
			
		END DO




END SUBROUTINE QUADRATURETRIANGLE


SUBROUTINE QUADRATUREQUAD(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for quadrilateral in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
REAL::R,S,T,a,b,c,d,e,f
REAL::a1,b1,c1,d1,e1,f1
INTEGER::Kk,J,ii,ij,ik,count1,alls


 WEQUA3D=0.0d0
  QPOINTS=0.0d0

SELECT CASE(IGQRULES)
 

 case(1)

		vvwg(1) = 4.0d0
	    VVR1(1)=0.0d0	;VVR2(1)=0.0d0	


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
  vvnpoy(1)= a	;vvnpoy(2)=b	
  vvnpoz(1)= a	;vvnpoz(2)=b

  vvwpox(1)=a1	;vvwpox(2)=b1	
  vvwpoy(1)= a1	;vvwpoy(2)=b1	
  vvwpoz(1)= a1	;vvwpoz(2)=b1
  
  count1=0
 do ii=1,2
    do ij=1,2
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
  	


  	
 CASE(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1
  count1=0
 do ii=1,3
    do ij=1,3
     
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
 

			
 CASE(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1
  count1=0
 do ii=1,4
    do ij=1,4
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
     
    end do
end do


			
CASE(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  E1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1
  count1=0
 do ii=1,5
    do ij=1,5
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=(vvwpox(ii)*vvwpox(ij))
      
    end do
end do
	
			
CASE(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  F=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  F1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e;vvnpoy(6)=f
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e;vvnpoz(6)=f
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1 ;vvwpoy(6)=f1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1 ;vvwpoz(6)=f1
  count1=0
 do ii=1,6
    do ij=1,6
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
	
			
END SELECT
		QPOINTS(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.25d0
! 		  WEQUA3D(:)=vvwg(:)
		do kk=1,qp_quad
			 WEQUA3D(kk)=vvwg(kk)
			R=VVR1(kk); S=VVR2(kk);
			VVnxi(1)=(0.25d0)*(1.0d0-R)*(1.0d0-s)
			VVnxi(2)=(0.25d0)*(1.0d0+R)*(1.0d0-s)
			VVnxi(3)=(0.25d0)*(1.0d0+R)*(1.0d0+s)
			VVnxi(4)=(0.25d0)*(1.0d0-R)*(1.0d0+s)
			
			DO J=1,4
			QPOINTS(1:2,kk)=QPOINTS(1:2,kk)+(VVNXI(j)*VEXT(j,1:2))
			END DO
! 			

		END DO
		




END SUBROUTINE QUADRATUREQUAD




SUBROUTINE QUADRATUREQUAD3D(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for quadrilateral in 3D
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
REAL::R,S,T,a,b,c,d,e,f
REAL::a1,b1,c1,d1,e1,f1
INTEGER::Kk,J,ii,ij,ik,count1,alls


 WEQUA2D=0.0d0
  QPOINTS2D=0.0d0

SELECT CASE(IGQRULES)
 

 case(1)

		vvwg(1) = 4.0d0
	    VVR1(1)=0.0d0	;VVR2(1)=0.0d0	


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
  vvnpoy(1)= a	;vvnpoy(2)=b	
  vvnpoz(1)= a	;vvnpoz(2)=b

  vvwpox(1)=a1	;vvwpox(2)=b1	
  vvwpoy(1)= a1	;vvwpoy(2)=b1	
  vvwpoz(1)= a1	;vvwpoz(2)=b1
  
  count1=0
 do ii=1,2
    do ij=1,2
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
  	


  	
 CASE(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1
  count1=0
 do ii=1,3
    do ij=1,3
     
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
 

			
 CASE(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1
  count1=0
 do ii=1,4
    do ij=1,4
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
     
    end do
end do


			
CASE(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  E1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1
  count1=0
 do ii=1,5
    do ij=1,5
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=(vvwpox(ii)*vvwpox(ij))
      
    end do
end do
	
			
CASE(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  F=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  F1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e;vvnpoy(6)=f
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e;vvnpoz(6)=f
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1 ;vvwpoy(6)=f1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1 ;vvwpoz(6)=f1
  count1=0
 do ii=1,6
    do ij=1,6
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)
      
    end do
end do
	
			
END SELECT
		QPOINTS2D(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.25d0
! 		  WEQUA2D(:)=vvwg(:)
		do kk=1,qp_quad
			WEQUA2D(kk)=vvwg(kk)
			R=VVR1(kk); S=VVR2(kk);
			VVnxi(1)=(0.25d0)*(1.0d0-R)*(1.0d0-s)
			VVnxi(2)=(0.25d0)*(1.0d0+R)*(1.0d0-s)
			VVnxi(3)=(0.25d0)*(1.0d0+R)*(1.0d0+s)
			VVnxi(4)=(0.25d0)*(1.0d0-R)*(1.0d0+s)
			
			DO J=1,4
			QPOINTS2D(1:3,kk)=QPOINTS2D(1:3,kk)+(VVNXI(j)*VEXT(j,1:3))
			END DO
! 			

		END DO
		




END SUBROUTINE QUADRATUREQUAD3D


SUBROUTINE QUADRATURELINE(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points for a line and returns it in QPOINTS2D(DIM,QP)
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
REAL::R,S,T,a,b,c,d,e,f,G,H,K
REAL::a1,b1,c1,d1,e1,f1,G1,H1,K1
INTEGER::Kk,J,ii,ij,ik,count1,alls


 WEQUA2D=0.0d0
  QPOINTS2D=0.0d0

SELECT CASE(IGQRULES)
 

 case(1)

		vvwg(1) = 2.0d0
	    VVR1(1)=0.0d0	;VVR2(1)=0.0d0	


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
 

  vvwpox(1)=a1	;vvwpox(2)=b1	
  
  
  count1=0
 do ii=1,2
         
	count1=count1+1
	VVR1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
    
end do
  	


  	
 CASE(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  
  count1=0
 do ii=1,3
    
     
	count1=count1+1
	VVR1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
    
end do
 

			
 CASE(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  
  count1=0
 do ii=1,4
   
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
     
   
end do


			
CASE(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  E1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
 
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  
  count1=0
 do ii=1,5
   
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii)
	vvwg(count1)=(vvwpox(ii))
      
    
end do
	
			
CASE(6)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  F=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  F1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  
  count1=0
 do ii=1,6
    
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
   
end do
	
	
	
CASE(7,8,9)

  
  
  A1=0.3302393550012598	         ;a=0.0000000000000000
	B1=0.1806481606948574	;b=-0.8360311073266358
	C1=0.1806481606948574	;C=0.8360311073266358
	D1=0.0812743883615744	;D=-0.9681602395076261
	E1=0.0812743883615744	;E=0.9681602395076261
	F1=0.3123470770400029	;F=-0.3242534234038089
	G1=0.3123470770400029	;G=0.3242534234038089
	H1=0.2606106964029354	;H=-0.6133714327005904
	K1=0.2606106964029354	;K=0.6133714327005904
  
  
  
  
  
  
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e; vvnpox(6)=f ; vvnpox(7)=G ; vvnpox(8)=H ; vvnpox(9)=K
  
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1 ;vvwpox(7)=G1 ;vvwpox(8)=H1 ;vvwpox(9)=K1
  
  count1=0
 do ii=1,9
    
      
	count1=count1+1
	VVR1(count1)=vvnpox(ii)
	vvwg(count1)=vvwpox(ii)
      
   
end do	
	
	
			
END SELECT
		QPOINTS2D(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.5d0
! 		  WEQUA2D(:)=vvwg(:)
		do kk=1,qp_LINE
			WEQUA2D(kk)=vvwg(kk)
			R=VVR1(kk); 
				
			
			
			QPOINTS2D(:,kk)=((VEXT(1,1:2)+VEXT(2,1:2))/2.0D0)+(R*(VEXT(2,1:2)-VEXT(1,1:2))/2.0D0)
			
			
! 		
		END DO
		




END SUBROUTINE QUADRATURELINE



SUBROUTINE QUADRATURETETRA(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for a tetrahedral
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
INTEGER::Kk




WEQUA3D=0.0d0
QPOINTS=0.0d0
! R1=VEXT2(N,:)-VEXT1(N,:)
! R2=VEXT3(N,:)-VEXT1(N,:)
! R3=VEXT4(N,:)-VEXT1(N,:)

select case(IGQRULES)
case(1)

		vvwg(1) = 1.0d0
	    VVR1(1)=0.25d0	;VVR2(1)=0.25d0	;VVR3(1)=0.25d0; VVR4(1)=0.25d0
case(2)
	
  	

VVR1(1)=0.5854101966249680 ;VVR2(1)=0.1381966011250110 ;VVR3(1)=0.1381966011250110 ;VVR4(1)=0.1381966011250110 ;vvwg(1)=0.2500000000000000
VVR1(2)=0.1381966011250110 ;VVR2(2)=0.5854101966249680 ;VVR3(2)=0.1381966011250110 ;VVR4(2)=0.1381966011250110 ;vvwg(2)=0.2500000000000000
VVR1(3)=0.1381966011250110 ;VVR2(3)=0.1381966011250110 ;VVR3(3)=0.5854101966249680 ;VVR4(3)=0.1381966011250110 ;vvwg(3)=0.2500000000000000
VVR1(4)=0.1381966011250110 ;VVR2(4)=0.1381966011250110 ;VVR3(4)=0.1381966011250110 ;VVR4(4)=0.5854101966249680 ;vvwg(4)=0.2500000000000000

  	
case(3)

  
VVR1( 1)= 0.3108859192633006E+00;VVR2( 1)= 0.3108859192633006E+00;VVR3( 1)= 0.6734224221009816E-01;VVR4( 1)= 0.3108859192633006E+00;VVWG( 1)= 0.1126879257180159E+00
VVR1( 2)= 0.3108859192633006E+00;VVR2( 2)= 0.6734224221009816E-01;VVR3( 2)= 0.3108859192633006E+00;VVR4( 2)= 0.3108859192633005E+00;VVWG( 2)= 0.1126879257180159E+00
VVR1( 3)= 0.6734224221009816E-01;VVR2( 3)= 0.3108859192633006E+00;VVR3( 3)= 0.3108859192633006E+00;VVR4( 3)= 0.3108859192633005E+00;VVWG( 3)= 0.1126879257180159E+00
VVR1( 4)= 0.3108859192633006E+00;VVR2( 4)= 0.3108859192633006E+00;VVR3( 4)= 0.3108859192633006E+00;VVR4( 4)= 0.6734224221009810E-01;VVWG( 4)= 0.1126879257180159E+00
VVR1( 5)= 0.9273525031089125E-01;VVR2( 5)= 0.9273525031089125E-01;VVR3( 5)= 0.7217942490673264E+00;VVR4( 5)= 0.9273525031089114E-01;VVWG( 5)= 0.7349304311636196E-01
VVR1( 6)= 0.9273525031089125E-01;VVR2( 6)= 0.7217942490673264E+00;VVR3( 6)= 0.9273525031089125E-01;VVR4( 6)= 0.9273525031089114E-01;VVWG( 6)= 0.7349304311636196E-01
VVR1( 7)= 0.7217942490673264E+00;VVR2( 7)= 0.9273525031089125E-01;VVR3( 7)= 0.9273525031089125E-01;VVR4( 7)= 0.9273525031089114E-01;VVWG( 7)= 0.7349304311636196E-01
VVR1( 8)= 0.9273525031089125E-01;VVR2( 8)= 0.9273525031089125E-01;VVR3( 8)= 0.9273525031089125E-01;VVR4( 8)= 0.7217942490673263E+00;VVWG( 8)= 0.7349304311636196E-01
VVR1( 9)= 0.4550370412564964E-01;VVR2( 9)= 0.4544962958743504E+00;VVR3( 9)= 0.4544962958743504E+00;VVR4( 9)= 0.4550370412564964E-01;VVWG( 9)= 0.4254602077708147E-01
VVR1(10)= 0.4544962958743504E+00;VVR2(10)= 0.4550370412564964E-01;VVR3(10)= 0.4544962958743504E+00;VVR4(10)= 0.4550370412564964E-01;VVWG(10)= 0.4254602077708147E-01
VVR1(11)= 0.4550370412564964E-01;VVR2(11)= 0.4550370412564964E-01;VVR3(11)= 0.4544962958743504E+00;VVR4(11)= 0.4544962958743504E+00;VVWG(11)= 0.4254602077708147E-01
VVR1(12)= 0.4550370412564964E-01;VVR2(12)= 0.4544962958743504E+00;VVR3(12)= 0.4550370412564964E-01;VVR4(12)= 0.4544962958743504E+00;VVWG(12)= 0.4254602077708147E-01
VVR1(13)= 0.4544962958743504E+00;VVR2(13)= 0.4550370412564964E-01;VVR3(13)= 0.4550370412564964E-01;VVR4(13)= 0.4544962958743504E+00;VVWG(13)= 0.4254602077708147E-01
VVR1(14)= 0.4544962958743504E+00;VVR2(14)= 0.4544962958743504E+00;VVR3(14)= 0.4550370412564964E-01;VVR4(14)= 0.4550370412564964E-01;VVWG(14)= 0.4254602077708147E-01



! 
! 
! 
! 
! 
! VVR1(1)=0.7784952948213300 ;VVR2(1)=0.0738349017262234 ;VVR3(1)=0.0738349017262234 ;VVR4(1)=0.0738349017262234 ;vvwg(1)=0.0476331348432089
! VVR1(2)=0.0738349017262234 ;VVR2(2)=0.7784952948213300 ;VVR3(2)=0.0738349017262234 ;VVR4(2)=0.0738349017262234 ;vvwg(2)=0.0476331348432089
! VVR1(3)=0.0738349017262234 ;VVR2(3)=0.0738349017262234 ;VVR3(3)=0.7784952948213300 ;VVR4(3)=0.0738349017262234 ;vvwg(3)=0.0476331348432089
! VVR1(4)=0.0738349017262234 ;VVR2(4)=0.0738349017262234 ;VVR3(4)=0.0738349017262234 ;VVR4(4)=0.7784952948213300 ;vvwg(4)=0.0476331348432089
! VVR1(5)=0.4062443438840510 ;VVR2(5)=0.4062443438840510 ;VVR3(5)=0.0937556561159491 ;VVR4(5)=0.0937556561159491 ;vvwg(5)=0.1349112434378610
! VVR1(6)=0.4062443438840510 ;VVR2(6)=0.0937556561159491 ;VVR3(6)=0.4062443438840510 ;VVR4(6)=0.0937556561159491 ;vvwg(6)=0.1349112434378610
! VVR1(7)=0.4062443438840510 ;VVR2(7)=0.0937556561159491 ;VVR3(7)=0.0937556561159491 ;VVR4(7)=0.4062443438840510 ;vvwg(7)=0.1349112434378610
! VVR1(8)=0.0937556561159491 ;VVR2(8)=0.4062443438840510 ;VVR3(8)=0.4062443438840510 ;VVR4(8)=0.0937556561159491 ;vvwg(8)=0.1349112434378610
! VVR1(9)=0.0937556561159491 ;VVR2(9)=0.4062443438840510 ;VVR3(9)=0.0937556561159491 ;VVR4(9)=0.4062443438840510 ;vvwg(9)=0.1349112434378610
! VVR1(10)=0.0937556561159491 ;VVR2(10)=0.0937556561159491 ;VVR3(10)=0.4062443438840510 ;VVR4(10)=0.4062443438840510 ;vvwg(10)=0.1349112434378610


  












			
case(4)
VVR1( 1)= 0.4067395853461137E-01;VVR2( 1)= 0.4067395853461137E-01;VVR3( 1)= 0.8779781243961660E+00;VVR4( 1)= 0.4067395853461120E-01;VVWG( 1)= 0.1007721105532064E-01
VVR1( 2)= 0.4067395853461137E-01;VVR2( 2)= 0.8779781243961660E+00;VVR3( 2)= 0.4067395853461137E-01;VVR4( 2)= 0.4067395853461125E-01;VVWG( 2)= 0.1007721105532064E-01
VVR1( 3)= 0.8779781243961660E+00;VVR2( 3)= 0.4067395853461137E-01;VVR3( 3)= 0.4067395853461137E-01;VVR4( 3)= 0.4067395853461131E-01;VVWG( 3)= 0.1007721105532064E-01
VVR1( 4)= 0.4067395853461137E-01;VVR2( 4)= 0.4067395853461137E-01;VVR3( 4)= 0.4067395853461137E-01;VVR4( 4)= 0.8779781243961657E+00;VVWG( 4)= 0.1007721105532064E-01
VVR1( 5)= 0.3223378901422755E+00;VVR2( 5)= 0.3223378901422755E+00;VVR3( 5)= 0.3298632957317349E-01;VVR4( 5)= 0.3223378901422754E+00;VVWG( 5)= 0.5535718154365472E-01
VVR1( 6)= 0.3223378901422755E+00;VVR2( 6)= 0.3298632957317349E-01;VVR3( 6)= 0.3223378901422755E+00;VVR4( 6)= 0.3223378901422754E+00;VVWG( 6)= 0.5535718154365472E-01
VVR1( 7)= 0.3298632957317349E-01;VVR2( 7)= 0.3223378901422755E+00;VVR3( 7)= 0.3223378901422755E+00;VVR4( 7)= 0.3223378901422754E+00;VVWG( 7)= 0.5535718154365472E-01
VVR1( 8)= 0.3223378901422755E+00;VVR2( 8)= 0.3223378901422755E+00;VVR3( 8)= 0.3223378901422755E+00;VVR4( 8)= 0.3298632957317338E-01;VVWG( 8)= 0.5535718154365472E-01
VVR1( 9)= 0.2146028712591520E+00;VVR2( 9)= 0.2146028712591520E+00;VVR3( 9)= 0.3561913862225439E+00;VVR4( 9)= 0.2146028712591521E+00;VVWG( 9)= 0.3992275025816749E-01
VVR1(10)= 0.2146028712591520E+00;VVR2(10)= 0.3561913862225439E+00;VVR3(10)= 0.2146028712591520E+00;VVR4(10)= 0.2146028712591521E+00;VVWG(10)= 0.3992275025816749E-01
VVR1(11)= 0.3561913862225439E+00;VVR2(11)= 0.2146028712591520E+00;VVR3(11)= 0.2146028712591520E+00;VVR4(11)= 0.2146028712591521E+00;VVWG(11)= 0.3992275025816749E-01
VVR1(12)= 0.2146028712591520E+00;VVR2(12)= 0.2146028712591520E+00;VVR3(12)= 0.2146028712591520E+00;VVR4(12)= 0.3561913862225440E+00;VVWG(12)= 0.3992275025816749E-01
VVR1(13)= 0.6030056647916491E+00;VVR2(13)= 0.6366100187501750E-01;VVR3(13)= 0.2696723314583158E+00;VVR4(13)= 0.6366100187501761E-01;VVWG(13)= 0.4821428571428571E-01
VVR1(14)= 0.6030056647916491E+00;VVR2(14)= 0.6366100187501750E-01;VVR3(14)= 0.6366100187501750E-01;VVR4(14)= 0.2696723314583159E+00;VVWG(14)= 0.4821428571428571E-01
VVR1(15)= 0.6366100187501750E-01;VVR2(15)= 0.6366100187501750E-01;VVR3(15)= 0.6030056647916491E+00;VVR4(15)= 0.2696723314583160E+00;VVWG(15)= 0.4821428571428571E-01
VVR1(16)= 0.2696723314583158E+00;VVR2(16)= 0.6030056647916491E+00;VVR3(16)= 0.6366100187501750E-01;VVR4(16)= 0.6366100187501761E-01;VVWG(16)= 0.4821428571428571E-01
VVR1(17)= 0.6366100187501750E-01;VVR2(17)= 0.2696723314583158E+00;VVR3(17)= 0.6030056647916491E+00;VVR4(17)= 0.6366100187501766E-01;VVWG(17)= 0.4821428571428571E-01
VVR1(18)= 0.6366100187501750E-01;VVR2(18)= 0.6030056647916491E+00;VVR3(18)= 0.6366100187501750E-01;VVR4(18)= 0.2696723314583160E+00;VVWG(18)= 0.4821428571428571E-01
VVR1(19)= 0.2696723314583158E+00;VVR2(19)= 0.6366100187501750E-01;VVR3(19)= 0.6030056647916491E+00;VVR4(19)= 0.6366100187501766E-01;VVWG(19)= 0.4821428571428571E-01
VVR1(20)= 0.6366100187501750E-01;VVR2(20)= 0.2696723314583158E+00;VVR3(20)= 0.6366100187501750E-01;VVR4(20)= 0.6030056647916493E+00;VVWG(20)= 0.4821428571428571E-01
VVR1(21)= 0.6366100187501750E-01;VVR2(21)= 0.6366100187501750E-01;VVR3(21)= 0.2696723314583158E+00;VVR4(21)= 0.6030056647916493E+00;VVWG(21)= 0.4821428571428571E-01
VVR1(22)= 0.6366100187501750E-01;VVR2(22)= 0.6030056647916491E+00;VVR3(22)= 0.2696723314583158E+00;VVR4(22)= 0.6366100187501766E-01;VVWG(22)= 0.4821428571428571E-01
VVR1(23)= 0.2696723314583158E+00;VVR2(23)= 0.6366100187501750E-01;VVR3(23)= 0.6366100187501750E-01;VVR4(23)= 0.6030056647916493E+00;VVWG(23)= 0.4821428571428571E-01
VVR1(24)= 0.6030056647916491E+00;VVR2(24)= 0.2696723314583158E+00;VVR3(24)= 0.6366100187501750E-01;VVR4(24)= 0.6366100187501761E-01;VVWG(24)= 0.4821428571428571E-01




! VVR1(1)=0.9029422158182680 ;VVR2(1)=0.0323525947272439 ;VVR3(1)=0.0323525947272439 ;VVR4(1)=0.0323525947272439 ;vvwg(1)=0.0070670747944695
! VVR1(2)=0.0323525947272439 ;VVR2(2)=0.9029422158182680 ;VVR3(2)=0.0323525947272439 ;VVR4(2)=0.0323525947272439 ;vvwg(2)=0.0070670747944695
! VVR1(3)=0.0323525947272439 ;VVR2(3)=0.0323525947272439 ;VVR3(3)=0.9029422158182680 ;VVR4(3)=0.0323525947272439 ;vvwg(3)=0.0070670747944695
! VVR1(4)=0.0323525947272439 ;VVR2(4)=0.0323525947272439 ;VVR3(4)=0.0323525947272439 ;VVR4(4)=0.9029422158182680 ;vvwg(4)=0.0070670747944695
! VVR1(5)=0.2626825838877790 ;VVR2(5)=0.6165965330619370 ;VVR3(5)=0.0603604415251421 ;VVR4(5)=0.0603604415251421 ;vvwg(5)=0.0469986689718877
! VVR1(6)=0.6165965330619370 ;VVR2(6)=0.2626825838877790 ;VVR3(6)=0.0603604415251421 ;VVR4(6)=0.0603604415251421 ;vvwg(6)=0.0469986689718877
! VVR1(7)=0.2626825838877790 ;VVR2(7)=0.0603604415251421 ;VVR3(7)=0.6165965330619370 ;VVR4(7)=0.0603604415251421 ;vvwg(7)=0.0469986689718877
! VVR1(8)=0.6165965330619370 ;VVR2(8)=0.0603604415251421 ;VVR3(8)=0.2626825838877790 ;VVR4(8)=0.0603604415251421 ;vvwg(8)=0.0469986689718877
! VVR1(9)=0.2626825838877790 ;VVR2(9)=0.0603604415251421 ;VVR3(9)=0.0603604415251421 ;VVR4(9)=0.6165965330619370 ;vvwg(9)=0.0469986689718877
! VVR1(10)=0.6165965330619370 ;VVR2(10)=0.0603604415251421 ;VVR3(10)=0.0603604415251421 ;VVR4(10)=0.2626825838877790 ;vvwg(10)=0.0469986689718877
! VVR1(11)=0.0603604415251421 ;VVR2(11)=0.2626825838877790 ;VVR3(11)=0.6165965330619370 ;VVR4(11)=0.0603604415251421 ;vvwg(11)=0.0469986689718877
! VVR1(12)=0.0603604415251421 ;VVR2(12)=0.6165965330619370 ;VVR3(12)=0.2626825838877790 ;VVR4(12)=0.0603604415251421 ;vvwg(12)=0.0469986689718877
! VVR1(13)=0.0603604415251421 ;VVR2(13)=0.2626825838877790 ;VVR3(13)=0.0603604415251421 ;VVR4(13)=0.6165965330619370 ;vvwg(13)=0.0469986689718877
! VVR1(14)=0.0603604415251421 ;VVR2(14)=0.6165965330619370 ;VVR3(14)=0.0603604415251421 ;VVR4(14)=0.2626825838877790 ;vvwg(14)=0.0469986689718877
! VVR1(15)=0.0603604415251421 ;VVR2(15)=0.0603604415251421 ;VVR3(15)=0.2626825838877790 ;VVR4(15)=0.6165965330619370 ;vvwg(15)=0.0469986689718877
! VVR1(16)=0.0603604415251421 ;VVR2(16)=0.0603604415251421 ;VVR3(16)=0.6165965330619370 ;VVR4(16)=0.2626825838877790 ;vvwg(16)=0.0469986689718877
! VVR1(17)=0.3097693042728620 ;VVR2(17)=0.3097693042728620 ;VVR3(17)=0.3097693042728620 ;VVR4(17)=0.0706920871814129 ;vvwg(17)=0.1019369182898680
! VVR1(18)=0.3097693042728620 ;VVR2(18)=0.3097693042728620 ;VVR3(18)=0.0706920871814129 ;VVR4(18)=0.3097693042728620 ;vvwg(18)=0.1019369182898680
! VVR1(19)=0.3097693042728620 ;VVR2(19)=0.0706920871814129 ;VVR3(19)=0.3097693042728620 ;VVR4(19)=0.3097693042728620 ;vvwg(19)=0.1019369182898680
! VVR1(20)=0.0706920871814129 ;VVR2(20)=0.3097693042728620 ;VVR3(20)=0.3097693042728620 ;VVR4(20)=0.3097693042728620 ;vvwg(20)=0.1019369182898680



  











			
case(5)


VVR1( 1)= 0.2500000000000000E+00;VVR2( 1)= 0.2500000000000000E+00;VVR3( 1)= 0.2500000000000000E+00;VVR4( 1)= 0.2500000000000000E+00;VVWG( 1)= 0.9548528946413085E-01
VVR1( 2)= 0.3157011497782028E+00;VVR2( 2)= 0.3157011497782028E+00;VVR3( 2)= 0.5289655066539162E-01;VVR4( 2)= 0.3157011497782028E+00;VVWG( 2)= 0.4232958120996703E-01
VVR1( 3)= 0.3157011497782028E+00;VVR2( 3)= 0.5289655066539162E-01;VVR3( 3)= 0.3157011497782028E+00;VVR4( 3)= 0.3157011497782029E+00;VVWG( 3)= 0.4232958120996703E-01
VVR1( 4)= 0.5289655066539162E-01;VVR2( 4)= 0.3157011497782028E+00;VVR3( 4)= 0.3157011497782028E+00;VVR4( 4)= 0.3157011497782029E+00;VVWG( 4)= 0.4232958120996703E-01
VVR1( 5)= 0.3157011497782028E+00;VVR2( 5)= 0.3157011497782028E+00;VVR3( 5)= 0.3157011497782028E+00;VVR4( 5)= 0.5289655066539167E-01;VVWG( 5)= 0.4232958120996703E-01
VVR1( 6)= 0.5048982259839635E-01;VVR2( 6)= 0.4495101774016036E+00;VVR3( 6)= 0.4495101774016036E+00;VVR4( 6)= 0.5048982259839652E-01;VVWG( 6)= 0.3189692783285758E-01
VVR1( 7)= 0.4495101774016036E+00;VVR2( 7)= 0.5048982259839635E-01;VVR3( 7)= 0.4495101774016036E+00;VVR4( 7)= 0.5048982259839641E-01;VVWG( 7)= 0.3189692783285758E-01
VVR1( 8)= 0.5048982259839635E-01;VVR2( 8)= 0.5048982259839635E-01;VVR3( 8)= 0.4495101774016036E+00;VVR4( 8)= 0.4495101774016038E+00;VVWG( 8)= 0.3189692783285758E-01
VVR1( 9)= 0.5048982259839635E-01;VVR2( 9)= 0.4495101774016036E+00;VVR3( 9)= 0.5048982259839635E-01;VVR4( 9)= 0.4495101774016038E+00;VVWG( 9)= 0.3189692783285758E-01
VVR1(10)= 0.4495101774016036E+00;VVR2(10)= 0.5048982259839635E-01;VVR3(10)= 0.5048982259839635E-01;VVR4(10)= 0.4495101774016036E+00;VVWG(10)= 0.3189692783285758E-01
VVR1(11)= 0.4495101774016036E+00;VVR2(11)= 0.4495101774016036E+00;VVR3(11)= 0.5048982259839635E-01;VVR4(11)= 0.5048982259839646E-01;VVWG(11)= 0.3189692783285758E-01
VVR1(12)= 0.5751716375870000E+00;VVR2(12)= 0.1888338310260010E+00;VVR3(12)= 0.4716070036099790E-01;VVR4(12)= 0.1888338310260011E+00;VVWG(12)= 0.3720713072833462E-01
VVR1(13)= 0.5751716375870000E+00;VVR2(13)= 0.1888338310260010E+00;VVR3(13)= 0.1888338310260010E+00;VVR4(13)= 0.4716070036099795E-01;VVWG(13)= 0.3720713072833462E-01
VVR1(14)= 0.1888338310260010E+00;VVR2(14)= 0.1888338310260010E+00;VVR3(14)= 0.5751716375870000E+00;VVR4(14)= 0.4716070036099795E-01;VVWG(14)= 0.3720713072833462E-01
VVR1(15)= 0.4716070036099790E-01;VVR2(15)= 0.5751716375870000E+00;VVR3(15)= 0.1888338310260010E+00;VVR4(15)= 0.1888338310260010E+00;VVWG(15)= 0.3720713072833462E-01
VVR1(16)= 0.1888338310260010E+00;VVR2(16)= 0.4716070036099790E-01;VVR3(16)= 0.5751716375870000E+00;VVR4(16)= 0.1888338310260012E+00;VVWG(16)= 0.3720713072833462E-01
VVR1(17)= 0.1888338310260010E+00;VVR2(17)= 0.5751716375870000E+00;VVR3(17)= 0.1888338310260010E+00;VVR4(17)= 0.4716070036099795E-01;VVWG(17)= 0.3720713072833462E-01
VVR1(18)= 0.4716070036099790E-01;VVR2(18)= 0.1888338310260010E+00;VVR3(18)= 0.5751716375870000E+00;VVR4(18)= 0.1888338310260010E+00;VVWG(18)= 0.3720713072833462E-01
VVR1(19)= 0.1888338310260010E+00;VVR2(19)= 0.4716070036099790E-01;VVR3(19)= 0.1888338310260010E+00;VVR4(19)= 0.5751716375870001E+00;VVWG(19)= 0.3720713072833462E-01
VVR1(20)= 0.1888338310260010E+00;VVR2(20)= 0.1888338310260010E+00;VVR3(20)= 0.4716070036099790E-01;VVR4(20)= 0.5751716375870000E+00;VVWG(20)= 0.3720713072833462E-01
VVR1(21)= 0.1888338310260010E+00;VVR2(21)= 0.5751716375870000E+00;VVR3(21)= 0.4716070036099790E-01;VVR4(21)= 0.1888338310260011E+00;VVWG(21)= 0.3720713072833462E-01
VVR1(22)= 0.4716070036099790E-01;VVR2(22)= 0.1888338310260010E+00;VVR3(22)= 0.1888338310260010E+00;VVR4(22)= 0.5751716375870000E+00;VVWG(22)= 0.3720713072833462E-01
VVR1(23)= 0.5751716375870000E+00;VVR2(23)= 0.4716070036099790E-01;VVR3(23)= 0.1888338310260010E+00;VVR4(23)= 0.1888338310260011E+00;VVWG(23)= 0.3720713072833462E-01
VVR1(24)= 0.8108302410985486E+00;VVR2(24)= 0.2126547254148325E-01;VVR3(24)= 0.1466388138184849E+00;VVR4(24)= 0.2126547254148320E-01;VVWG(24)= 0.8110770829903342E-02
VVR1(25)= 0.8108302410985486E+00;VVR2(25)= 0.2126547254148325E-01;VVR3(25)= 0.2126547254148325E-01;VVR4(25)= 0.1466388138184849E+00;VVWG(25)= 0.8110770829903342E-02
VVR1(26)= 0.2126547254148325E-01;VVR2(26)= 0.2126547254148325E-01;VVR3(26)= 0.8108302410985486E+00;VVR4(26)= 0.1466388138184849E+00;VVWG(26)= 0.8110770829903342E-02
VVR1(27)= 0.1466388138184849E+00;VVR2(27)= 0.8108302410985486E+00;VVR3(27)= 0.2126547254148325E-01;VVR4(27)= 0.2126547254148325E-01;VVWG(27)= 0.8110770829903342E-02
VVR1(28)= 0.2126547254148325E-01;VVR2(28)= 0.1466388138184849E+00;VVR3(28)= 0.8108302410985486E+00;VVR4(28)= 0.2126547254148314E-01;VVWG(28)= 0.8110770829903342E-02
VVR1(29)= 0.2126547254148325E-01;VVR2(29)= 0.8108302410985486E+00;VVR3(29)= 0.2126547254148325E-01;VVR4(29)= 0.1466388138184849E+00;VVWG(29)= 0.8110770829903342E-02
VVR1(30)= 0.1466388138184849E+00;VVR2(30)= 0.2126547254148325E-01;VVR3(30)= 0.8108302410985486E+00;VVR4(30)= 0.2126547254148325E-01;VVWG(30)= 0.8110770829903342E-02
VVR1(31)= 0.2126547254148325E-01;VVR2(31)= 0.1466388138184849E+00;VVR3(31)= 0.2126547254148325E-01;VVR4(31)= 0.8108302410985485E+00;VVWG(31)= 0.8110770829903342E-02
VVR1(32)= 0.2126547254148325E-01;VVR2(32)= 0.2126547254148325E-01;VVR3(32)= 0.1466388138184849E+00;VVR4(32)= 0.8108302410985486E+00;VVWG(32)= 0.8110770829903342E-02
VVR1(33)= 0.2126547254148325E-01;VVR2(33)= 0.8108302410985486E+00;VVR3(33)= 0.1466388138184849E+00;VVR4(33)= 0.2126547254148320E-01;VVWG(33)= 0.8110770829903342E-02
VVR1(34)= 0.1466388138184849E+00;VVR2(34)= 0.2126547254148325E-01;VVR3(34)= 0.2126547254148325E-01;VVR4(34)= 0.8108302410985486E+00;VVWG(34)= 0.8110770829903342E-02
VVR1(35)= 0.8108302410985486E+00;VVR2(35)= 0.1466388138184849E+00;VVR3(35)= 0.2126547254148325E-01;VVR4(35)= 0.2126547254148320E-01;VVWG(35)= 0.8110770829903342E-02



! VVR1(1)=0.9197896733368800 ;VVR2(1)=0.0267367755543735 ;VVR3(1)=0.0267367755543735 ;VVR4(1)=0.0267367755543735 ;vvwg(1)=0.0021900463965388
! VVR1(2)=0.0267367755543735 ;VVR2(2)=0.9197896733368800 ;VVR3(2)=0.0267367755543735 ;VVR4(2)=0.0267367755543735 ;vvwg(2)=0.0021900463965388
! VVR1(3)=0.0267367755543735 ;VVR2(3)=0.0267367755543735 ;VVR3(3)=0.9197896733368800 ;VVR4(3)=0.0267367755543735 ;vvwg(3)=0.0021900463965388
! VVR1(4)=0.0267367755543735 ;VVR2(4)=0.0267367755543735 ;VVR3(4)=0.0267367755543735 ;VVR4(4)=0.9197896733368800 ;vvwg(4)=0.0021900463965388
! VVR1(5)=0.1740356302468940 ;VVR2(5)=0.7477598884818090 ;VVR3(5)=0.0391022406356488 ;VVR4(5)=0.0391022406356488 ;vvwg(5)=0.0143395670177665
! VVR1(6)=0.7477598884818090 ;VVR2(6)=0.1740356302468940 ;VVR3(6)=0.0391022406356488 ;VVR4(6)=0.0391022406356488 ;vvwg(6)=0.0143395670177665
! VVR1(7)=0.1740356302468940 ;VVR2(7)=0.0391022406356488 ;VVR3(7)=0.7477598884818090 ;VVR4(7)=0.0391022406356488 ;vvwg(7)=0.0143395670177665
! VVR1(8)=0.7477598884818090 ;VVR2(8)=0.0391022406356488 ;VVR3(8)=0.1740356302468940 ;VVR4(8)=0.0391022406356488 ;vvwg(8)=0.0143395670177665
! VVR1(9)=0.1740356302468940 ;VVR2(9)=0.0391022406356488 ;VVR3(9)=0.0391022406356488 ;VVR4(9)=0.7477598884818090 ;vvwg(9)=0.0143395670177665
! VVR1(10)=0.7477598884818090 ;VVR2(10)=0.0391022406356488 ;VVR3(10)=0.0391022406356488 ;VVR4(10)=0.1740356302468940 ;vvwg(10)=0.0143395670177665
! VVR1(11)=0.0391022406356488 ;VVR2(11)=0.1740356302468940 ;VVR3(11)=0.7477598884818090 ;VVR4(11)=0.0391022406356488 ;vvwg(11)=0.0143395670177665
! VVR1(12)=0.0391022406356488 ;VVR2(12)=0.7477598884818090 ;VVR3(12)=0.1740356302468940 ;VVR4(12)=0.0391022406356488 ;vvwg(12)=0.0143395670177665
! VVR1(13)=0.0391022406356488 ;VVR2(13)=0.1740356302468940 ;VVR3(13)=0.0391022406356488 ;VVR4(13)=0.7477598884818090 ;vvwg(13)=0.0143395670177665
! VVR1(14)=0.0391022406356488 ;VVR2(14)=0.7477598884818090 ;VVR3(14)=0.0391022406356488 ;VVR4(14)=0.1740356302468940 ;vvwg(14)=0.0143395670177665
! VVR1(15)=0.0391022406356488 ;VVR2(15)=0.0391022406356488 ;VVR3(15)=0.1740356302468940 ;VVR4(15)=0.7477598884818090 ;vvwg(15)=0.0143395670177665
! VVR1(16)=0.0391022406356488 ;VVR2(16)=0.0391022406356488 ;VVR3(16)=0.7477598884818090 ;VVR4(16)=0.1740356302468940 ;vvwg(16)=0.0143395670177665
! VVR1(17)=0.4547545999844830 ;VVR2(17)=0.4547545999844830 ;VVR3(17)=0.0452454000155172 ;VVR4(17)=0.0452454000155172 ;vvwg(17)=0.0250305395686746
! VVR1(18)=0.4547545999844830 ;VVR2(18)=0.0452454000155172 ;VVR3(18)=0.4547545999844830 ;VVR4(18)=0.0452454000155172 ;vvwg(18)=0.0250305395686746
! VVR1(19)=0.4547545999844830 ;VVR2(19)=0.0452454000155172 ;VVR3(19)=0.0452454000155172 ;VVR4(19)=0.4547545999844830 ;vvwg(19)=0.0250305395686746
! VVR1(20)=0.0452454000155172 ;VVR2(20)=0.4547545999844830 ;VVR3(20)=0.4547545999844830 ;VVR4(20)=0.0452454000155172 ;vvwg(20)=0.0250305395686746
! VVR1(21)=0.0452454000155172 ;VVR2(21)=0.4547545999844830 ;VVR3(21)=0.0452454000155172 ;VVR4(21)=0.4547545999844830 ;vvwg(21)=0.0250305395686746
! VVR1(22)=0.0452454000155172 ;VVR2(22)=0.0452454000155172 ;VVR3(22)=0.4547545999844830 ;VVR4(22)=0.4547545999844830 ;vvwg(22)=0.0250305395686746
! VVR1(23)=0.5031186450145980 ;VVR2(23)=0.2232010379623150 ;VVR3(23)=0.2232010379623150 ;VVR4(23)=0.0504792790607720 ;vvwg(23)=0.0479839333057554
! VVR1(24)=0.2232010379623150 ;VVR2(24)=0.5031186450145980 ;VVR3(24)=0.2232010379623150 ;VVR4(24)=0.0504792790607720 ;vvwg(24)=0.0479839333057554
! VVR1(25)=0.2232010379623150 ;VVR2(25)=0.2232010379623150 ;VVR3(25)=0.5031186450145980 ;VVR4(25)=0.0504792790607720 ;vvwg(25)=0.0479839333057554
! VVR1(26)=0.5031186450145980 ;VVR2(26)=0.2232010379623150 ;VVR3(26)=0.0504792790607720 ;VVR4(26)=0.2232010379623150 ;vvwg(26)=0.0479839333057554
! VVR1(27)=0.2232010379623150 ;VVR2(27)=0.5031186450145980 ;VVR3(27)=0.0504792790607720 ;VVR4(27)=0.2232010379623150 ;vvwg(27)=0.0479839333057554
! VVR1(28)=0.2232010379623150 ;VVR2(28)=0.2232010379623150 ;VVR3(28)=0.0504792790607720 ;VVR4(28)=0.5031186450145980 ;vvwg(28)=0.0479839333057554
! VVR1(29)=0.5031186450145980 ;VVR2(29)=0.0504792790607720 ;VVR3(29)=0.2232010379623150 ;VVR4(29)=0.2232010379623150 ;vvwg(29)=0.0479839333057554
! VVR1(30)=0.2232010379623150 ;VVR2(30)=0.0504792790607720 ;VVR3(30)=0.5031186450145980 ;VVR4(30)=0.2232010379623150 ;vvwg(30)=0.0479839333057554
! VVR1(31)=0.2232010379623150 ;VVR2(31)=0.0504792790607720 ;VVR3(31)=0.2232010379623150 ;VVR4(31)=0.5031186450145980 ;vvwg(31)=0.0479839333057554
! VVR1(32)=0.0504792790607720 ;VVR2(32)=0.5031186450145980 ;VVR3(32)=0.2232010379623150 ;VVR4(32)=0.2232010379623150 ;vvwg(32)=0.0479839333057554
! VVR1(33)=0.0504792790607720 ;VVR2(33)=0.2232010379623150 ;VVR3(33)=0.5031186450145980 ;VVR4(33)=0.2232010379623150 ;vvwg(33)=0.0479839333057554
! VVR1(34)=0.0504792790607720 ;VVR2(34)=0.2232010379623150 ;VVR3(34)=0.2232010379623150 ;VVR4(34)=0.5031186450145980 ;vvwg(34)=0.0479839333057554
! VVR1(35)=0.2500000000000000 ;VVR2(35)=0.2500000000000000 ;VVR3(35)=0.2500000000000000 ;VVR4(35)=0.2500000000000000 ;vvwg(35)=0.0931745731195340
	
	
CASE(6,7,8,9)
VVR1( 1)= 0.9551438045408220E+00;VVR2( 1)= 0.1495206515305920E-01;VVR3( 1)= 0.1495206515305920E-01;VVR4( 1)= 0.1495206515305920E-01;VVWG( 1)= 0.1037311233614000E-02
VVR1( 2)= 0.1495206515305920E-01;VVR2( 2)= 0.9551438045408220E+00;VVR3( 2)= 0.1495206515305920E-01;VVR4( 2)= 0.1495206515305920E-01;VVWG( 2)= 0.1037311233614000E-02
VVR1( 3)= 0.1495206515305920E-01;VVR2( 3)= 0.1495206515305920E-01;VVR3( 3)= 0.9551438045408220E+00;VVR4( 3)= 0.1495206515305920E-01;VVWG( 3)= 0.1037311233614000E-02
VVR1( 4)= 0.1495206515305920E-01;VVR2( 4)= 0.1495206515305920E-01;VVR3( 4)= 0.1495206515305920E-01;VVR4( 4)= 0.9551438045408220E+00;VVWG( 4)= 0.1037311233614000E-02
VVR1( 5)= 0.7799760084415400E+00;VVR2( 5)= 0.1518319491659370E+00;VVR3( 5)= 0.3409602119626150E-01;VVR4( 5)= 0.3409602119626150E-01;VVWG( 5)= 0.9601664539948001E-02
VVR1( 6)= 0.1518319491659370E+00;VVR2( 6)= 0.7799760084415400E+00;VVR3( 6)= 0.3409602119626150E-01;VVR4( 6)= 0.3409602119626150E-01;VVWG( 6)= 0.9601664539948001E-02
VVR1( 7)= 0.7799760084415400E+00;VVR2( 7)= 0.3409602119626150E-01;VVR3( 7)= 0.1518319491659370E+00;VVR4( 7)= 0.3409602119626150E-01;VVWG( 7)= 0.9601664539948001E-02
VVR1( 8)= 0.1518319491659370E+00;VVR2( 8)= 0.3409602119626150E-01;VVR3( 8)= 0.7799760084415400E+00;VVR4( 8)= 0.3409602119626150E-01;VVWG( 8)= 0.9601664539948001E-02
VVR1( 9)= 0.7799760084415400E+00;VVR2( 9)= 0.3409602119626150E-01;VVR3( 9)= 0.3409602119626150E-01;VVR4( 9)= 0.1518319491659370E+00;VVWG( 9)= 0.9601664539948001E-02
VVR1(10)= 0.1518319491659370E+00;VVR2(10)= 0.3409602119626150E-01;VVR3(10)= 0.3409602119626150E-01;VVR4(10)= 0.7799760084415400E+00;VVWG(10)= 0.9601664539948001E-02
VVR1(11)= 0.3409602119626150E-01;VVR2(11)= 0.7799760084415400E+00;VVR3(11)= 0.1518319491659370E+00;VVR4(11)= 0.3409602119626150E-01;VVWG(11)= 0.9601664539948001E-02
VVR1(12)= 0.3409602119626150E-01;VVR2(12)= 0.1518319491659370E+00;VVR3(12)= 0.7799760084415400E+00;VVR4(12)= 0.3409602119626150E-01;VVWG(12)= 0.9601664539948001E-02
VVR1(13)= 0.3409602119626150E-01;VVR2(13)= 0.7799760084415400E+00;VVR3(13)= 0.3409602119626150E-01;VVR4(13)= 0.1518319491659370E+00;VVWG(13)= 0.9601664539948001E-02
VVR1(14)= 0.3409602119626150E-01;VVR2(14)= 0.1518319491659370E+00;VVR3(14)= 0.3409602119626150E-01;VVR4(14)= 0.7799760084415400E+00;VVWG(14)= 0.9601664539948001E-02
VVR1(15)= 0.3409602119626150E-01;VVR2(15)= 0.3409602119626150E-01;VVR3(15)= 0.7799760084415400E+00;VVR4(15)= 0.1518319491659370E+00;VVWG(15)= 0.9601664539948001E-02
VVR1(16)= 0.3409602119626150E-01;VVR2(16)= 0.3409602119626150E-01;VVR3(16)= 0.1518319491659370E+00;VVR4(16)= 0.7799760084415400E+00;VVWG(16)= 0.9601664539948001E-02
VVR1(17)= 0.3549340560639790E+00;VVR2(17)= 0.5526556431060170E+00;VVR3(17)= 0.4620515041500170E-01;VVR4(17)= 0.4620515041500170E-01;VVWG(17)= 0.1644939767982320E-01
VVR1(18)= 0.5526556431060170E+00;VVR2(18)= 0.3549340560639790E+00;VVR3(18)= 0.4620515041500170E-01;VVR4(18)= 0.4620515041500170E-01;VVWG(18)= 0.1644939767982320E-01
VVR1(19)= 0.3549340560639790E+00;VVR2(19)= 0.4620515041500170E-01;VVR3(19)= 0.5526556431060170E+00;VVR4(19)= 0.4620515041500170E-01;VVWG(19)= 0.1644939767982320E-01
VVR1(20)= 0.5526556431060170E+00;VVR2(20)= 0.4620515041500170E-01;VVR3(20)= 0.3549340560639790E+00;VVR4(20)= 0.4620515041500170E-01;VVWG(20)= 0.1644939767982320E-01
VVR1(21)= 0.3549340560639790E+00;VVR2(21)= 0.4620515041500170E-01;VVR3(21)= 0.4620515041500170E-01;VVR4(21)= 0.5526556431060170E+00;VVWG(21)= 0.1644939767982320E-01
VVR1(22)= 0.5526556431060170E+00;VVR2(22)= 0.4620515041500170E-01;VVR3(22)= 0.4620515041500170E-01;VVR4(22)= 0.3549340560639790E+00;VVWG(22)= 0.1644939767982320E-01
VVR1(23)= 0.4620515041500170E-01;VVR2(23)= 0.3549340560639790E+00;VVR3(23)= 0.5526556431060170E+00;VVR4(23)= 0.4620515041500170E-01;VVWG(23)= 0.1644939767982320E-01
VVR1(24)= 0.4620515041500170E-01;VVR2(24)= 0.5526556431060170E+00;VVR3(24)= 0.3549340560639790E+00;VVR4(24)= 0.4620515041500170E-01;VVWG(24)= 0.1644939767982320E-01
VVR1(25)= 0.4620515041500170E-01;VVR2(25)= 0.3549340560639790E+00;VVR3(25)= 0.4620515041500170E-01;VVR4(25)= 0.5526556431060170E+00;VVWG(25)= 0.1644939767982320E-01
VVR1(26)= 0.4620515041500170E-01;VVR2(26)= 0.5526556431060170E+00;VVR3(26)= 0.4620515041500170E-01;VVR4(26)= 0.3549340560639790E+00;VVWG(26)= 0.1644939767982320E-01
VVR1(27)= 0.4620515041500170E-01;VVR2(27)= 0.4620515041500170E-01;VVR3(27)= 0.3549340560639790E+00;VVR4(27)= 0.5526556431060170E+00;VVWG(27)= 0.1644939767982320E-01
VVR1(28)= 0.4620515041500170E-01;VVR2(28)= 0.4620515041500170E-01;VVR3(28)= 0.5526556431060170E+00;VVR4(28)= 0.3549340560639790E+00;VVWG(28)= 0.1644939767982320E-01
VVR1(29)= 0.5381043228880020E+00;VVR2(29)= 0.2281904610687610E+00;VVR3(29)= 0.2281904610687610E+00;VVR4(29)= 0.5514754974477500E-02;VVWG(29)= 0.1537477665133100E-01
VVR1(30)= 0.2281904610687610E+00;VVR2(30)= 0.5381043228880020E+00;VVR3(30)= 0.2281904610687610E+00;VVR4(30)= 0.5514754974477500E-02;VVWG(30)= 0.1537477665133100E-01
VVR1(31)= 0.2281904610687610E+00;VVR2(31)= 0.2281904610687610E+00;VVR3(31)= 0.5381043228880020E+00;VVR4(31)= 0.5514754974477500E-02;VVWG(31)= 0.1537477665133100E-01
VVR1(32)= 0.5381043228880020E+00;VVR2(32)= 0.2281904610687610E+00;VVR3(32)= 0.5514754974477500E-02;VVR4(32)= 0.2281904610687610E+00;VVWG(32)= 0.1537477665133100E-01
VVR1(33)= 0.2281904610687610E+00;VVR2(33)= 0.5381043228880020E+00;VVR3(33)= 0.5514754974477500E-02;VVR4(33)= 0.2281904610687610E+00;VVWG(33)= 0.1537477665133100E-01
VVR1(34)= 0.2281904610687610E+00;VVR2(34)= 0.2281904610687610E+00;VVR3(34)= 0.5514754974477500E-02;VVR4(34)= 0.5381043228880020E+00;VVWG(34)= 0.1537477665133100E-01
VVR1(35)= 0.5381043228880020E+00;VVR2(35)= 0.5514754974477500E-02;VVR3(35)= 0.2281904610687610E+00;VVR4(35)= 0.2281904610687610E+00;VVWG(35)= 0.1537477665133100E-01
VVR1(36)= 0.2281904610687610E+00;VVR2(36)= 0.5514754974477500E-02;VVR3(36)= 0.5381043228880020E+00;VVR4(36)= 0.2281904610687610E+00;VVWG(36)= 0.1537477665133100E-01
VVR1(37)= 0.2281904610687610E+00;VVR2(37)= 0.5514754974477500E-02;VVR3(37)= 0.2281904610687610E+00;VVR4(37)= 0.5381043228880020E+00;VVWG(37)= 0.1537477665133100E-01
VVR1(38)= 0.5514754974477500E-02;VVR2(38)= 0.5381043228880020E+00;VVR3(38)= 0.2281904610687610E+00;VVR4(38)= 0.2281904610687610E+00;VVWG(38)= 0.1537477665133100E-01
VVR1(39)= 0.5514754974477500E-02;VVR2(39)= 0.2281904610687610E+00;VVR3(39)= 0.5381043228880020E+00;VVR4(39)= 0.2281904610687610E+00;VVWG(39)= 0.1537477665133100E-01
VVR1(40)= 0.5514754974477500E-02;VVR2(40)= 0.2281904610687610E+00;VVR3(40)= 0.2281904610687610E+00;VVR4(40)= 0.5381043228880020E+00;VVWG(40)= 0.1537477665133100E-01
VVR1(41)= 0.1961837595745600E+00;VVR2(41)= 0.3523052600879940E+00;VVR3(41)= 0.3523052600879940E+00;VVR4(41)= 0.9920572024945300E-01;VVWG(41)= 0.2935201183752300E-01
VVR1(42)= 0.3523052600879940E+00;VVR2(42)= 0.1961837595745600E+00;VVR3(42)= 0.3523052600879940E+00;VVR4(42)= 0.9920572024945300E-01;VVWG(42)= 0.2935201183752300E-01
VVR1(43)= 0.3523052600879940E+00;VVR2(43)= 0.3523052600879940E+00;VVR3(43)= 0.1961837595745600E+00;VVR4(43)= 0.9920572024945300E-01;VVWG(43)= 0.2935201183752300E-01
VVR1(44)= 0.1961837595745600E+00;VVR2(44)= 0.3523052600879940E+00;VVR3(44)= 0.9920572024945300E-01;VVR4(44)= 0.3523052600879940E+00;VVWG(44)= 0.2935201183752300E-01
VVR1(45)= 0.3523052600879940E+00;VVR2(45)= 0.1961837595745600E+00;VVR3(45)= 0.9920572024945300E-01;VVR4(45)= 0.3523052600879940E+00;VVWG(45)= 0.2935201183752300E-01
VVR1(46)= 0.3523052600879940E+00;VVR2(46)= 0.3523052600879940E+00;VVR3(46)= 0.9920572024945300E-01;VVR4(46)= 0.1961837595745600E+00;VVWG(46)= 0.2935201183752300E-01
VVR1(47)= 0.1961837595745600E+00;VVR2(47)= 0.9920572024945300E-01;VVR3(47)= 0.3523052600879940E+00;VVR4(47)= 0.3523052600879940E+00;VVWG(47)= 0.2935201183752300E-01
VVR1(48)= 0.3523052600879940E+00;VVR2(48)= 0.9920572024945300E-01;VVR3(48)= 0.1961837595745600E+00;VVR4(48)= 0.3523052600879940E+00;VVWG(48)= 0.2935201183752300E-01
VVR1(49)= 0.3523052600879940E+00;VVR2(49)= 0.9920572024945300E-01;VVR3(49)= 0.3523052600879940E+00;VVR4(49)= 0.1961837595745600E+00;VVWG(49)= 0.2935201183752300E-01
VVR1(50)= 0.9920572024945300E-01;VVR2(50)= 0.1961837595745600E+00;VVR3(50)= 0.3523052600879940E+00;VVR4(50)= 0.3523052600879940E+00;VVWG(50)= 0.2935201183752300E-01
VVR1(51)= 0.9920572024945300E-01;VVR2(51)= 0.3523052600879940E+00;VVR3(51)= 0.1961837595745600E+00;VVR4(51)= 0.3523052600879940E+00;VVWG(51)= 0.2935201183752300E-01
VVR1(52)= 0.9920572024945300E-01;VVR2(52)= 0.3523052600879940E+00;VVR3(52)= 0.3523052600879940E+00;VVR4(52)= 0.1961837595745600E+00;VVWG(52)= 0.2935201183752300E-01
VVR1(53)= 0.5965649956210169E+00;VVR2(53)= 0.1344783347929940E+00;VVR3(53)= 0.1344783347929940E+00;VVR4(53)= 0.1344783347929940E+00;VVWG(53)= 0.3662913664051080E-01
VVR1(54)= 0.1344783347929940E+00;VVR2(54)= 0.5965649956210169E+00;VVR3(54)= 0.1344783347929940E+00;VVR4(54)= 0.1344783347929940E+00;VVWG(54)= 0.3662913664051080E-01
VVR1(55)= 0.1344783347929940E+00;VVR2(55)= 0.1344783347929940E+00;VVR3(55)= 0.5965649956210169E+00;VVR4(55)= 0.1344783347929940E+00;VVWG(55)= 0.3662913664051080E-01
VVR1(56)= 0.1344783347929940E+00;VVR2(56)= 0.1344783347929940E+00;VVR3(56)= 0.1344783347929940E+00;VVR4(56)= 0.5965649956210169E+00;VVWG(56)= 0.3662913664051080E-01
 	
	
	
	
	
	
			
END select



		do kk=1,qp_tetra
			WEQUA3D(kk)=vvwg(kk)
			QPOINTS(:,kk)=(VVR1(kk)*VEXT(1,:))+(VVR2(kk)*VEXT(2,:))+(VVR3(kk)*VEXT(3,:))+(VVR4(kk)*VEXT(4,:))
			
		
		END DO


END SUBROUTINE QUADRATURETETRA

SUBROUTINE QUADRATUREPRISM(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for a prism
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
REAL::R,S,T,a,b,c,d,e,f
REAL::a1,b1,c1,d1,e1,f1,sumwe
INTEGER::Kk,J,ii,ij,ik,count1,alls

alls=igqrules*QP_TRIANGLE


 
 WEQUA3D=0.0d0
  QPOINTS=0.0d0
sumwe=0.0d0

SELECT CASE(IGQRULES)
 
  case(1)
		VVR1(1)=0.666666666666667/2.0d0 ;VVR2(1)=0.666666666666667/2.0d0 ;VVR3(1)=0.0d0
		
		  vvwg(1)=2.0d0

 
  	


 

 
 case(2)
		VVR1(1)=0.666666666666667 ;VVR2(1)=0.166666666666667 ;VVR3(1)=0.5773502691896257;vvwg(1)=0.33333333333333333333*1.0d0
		VVR1(2)=0.166666666666667 ;VVR2(2)=0.666666666666667 ;VVR3(2)=0.5773502691896257;vvwg(2)=0.33333333333333333333*1.0d0
		VVR1(3)=0.166666666666667 ;VVR2(3)=0.166666666666667 ;VVR3(3)=0.5773502691896257;vvwg(3)=0.33333333333333333333*1.0d0
		VVR1(4)=0.666666666666667 ;VVR2(4)=0.166666666666667 ;VVR3(4)=-0.5773502691896257;vvwg(4)=0.33333333333333333333*1.0d0
		VVR1(5)=0.166666666666667 ;VVR2(5)=0.666666666666667 ;VVR3(5)=-0.5773502691896257;vvwg(5)=0.33333333333333333333*1.0d0
		VVR1(6)=0.166666666666667 ;VVR2(6)=0.166666666666667 ;VVR3(6)=-0.5773502691896257;vvwg(6)=0.33333333333333333333*1.0d0
 

 
  	


  	
 CASE(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
 

		VVR1(1)=0.816847572980440 ;VVR2(1)=0.091576213509780 ;VVR3(1)=0.0 ;vvwg(1)=0.109951743655333*0.8888888888888888
		VVR1(2)=0.091576213509780 ;VVR2(2)=0.816847572980440 ;VVR3(2)=0.0 ;vvwg(2)=0.109951743655333*0.8888888888888888
		VVR1(3)=0.091576213509780 ;VVR2(3)=0.091576213509780 ;VVR3(3)=0.0 ;vvwg(3)=0.109951743655333*0.8888888888888888
		VVR1(4)=0.445948490915964 ;VVR2(4)=0.445948490915964 ;VVR3(4)=0.0 ;vvwg(4)=0.223381589678000*0.8888888888888888
		VVR1(5)=0.445948490915964 ;VVR2(5)=0.108103018168071 ;VVR3(5)=0.0 ;vvwg(5)=0.223381589678000*0.8888888888888888
		VVR1(6)=0.108103018168071 ;VVR2(6)=0.445948490915964 ;VVR3(6)=0.0 ;vvwg(6)=0.223381589678000*0.8888888888888888
		VVR1(7)=0.816847572980440 ;VVR2(7)=0.091576213509780 ;VVR3(7)=-0.7745966692414834 ;vvwg(7)=0.109951743655333*0.5555555555555556
		VVR1(8)=0.091576213509780 ;VVR2(8)=0.816847572980440 ;VVR3(8)=-0.7745966692414834 ;vvwg(8)=0.109951743655333*0.5555555555555556
		VVR1(9)=0.091576213509780 ;VVR2(9)=0.091576213509780 ;VVR3(9)=-0.7745966692414834 ;vvwg(9)=0.109951743655333*0.5555555555555556
		VVR1(10)=0.445948490915964 ;VVR2(10)=0.445948490915964 ;VVR3(10)=-0.7745966692414834 ;vvwg(10)=0.223381589678000*0.5555555555555556
		VVR1(11)=0.445948490915964 ;VVR2(11)=0.108103018168071 ;VVR3(11)=-0.7745966692414834 ;vvwg(11)=0.223381589678000*0.5555555555555556
		VVR1(12)=0.108103018168071 ;VVR2(12)=0.445948490915964 ;VVR3(12)=-0.7745966692414834 ;vvwg(12)=0.223381589678000*0.5555555555555556
		VVR1(13)=0.816847572980440 ;VVR2(13)=0.091576213509780 ;VVR3(13)=0.7745966692414834 ;vvwg(13)=0.109951743655333*0.5555555555555556
		VVR1(14)=0.091576213509780 ;VVR2(14)=0.816847572980440 ;VVR3(14)=0.7745966692414834 ;vvwg(14)=0.109951743655333*0.5555555555555556
		VVR1(15)=0.091576213509780 ;VVR2(15)=0.091576213509780 ;VVR3(15)=0.7745966692414834 ;vvwg(15)=0.109951743655333*0.5555555555555556
		VVR1(16)=0.445948490915964 ;VVR2(16)=0.445948490915964 ;VVR3(16)=0.7745966692414834 ;vvwg(16)=0.223381589678000*0.5555555555555556
		VVR1(17)=0.445948490915964 ;VVR2(17)=0.108103018168071 ;VVR3(17)=0.7745966692414834 ;vvwg(17)=0.223381589678000*0.5555555555555556
		VVR1(18)=0.108103018168071 ;VVR2(18)=0.445948490915964 ;VVR3(18)=0.7745966692414834 ;vvwg(18)=0.223381589678000*0.5555555555555556











			
 CASE(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
 
		VVR1(1)=0.816847572980440 ;VVR2(1)=0.091576213509780 ;VVR3(1)=-0.3399810435848563 ;vvwg(1)=0.109951743655333*0.6521451548625461
		VVR1(2)=0.091576213509780 ;VVR2(2)=0.816847572980440 ;VVR3(2)=-0.3399810435848563 ;vvwg(2)=0.109951743655333*0.6521451548625461
		VVR1(3)=0.091576213509780 ;VVR2(3)=0.091576213509780 ;VVR3(3)=-0.3399810435848563 ;vvwg(3)=0.109951743655333*0.6521451548625461
		VVR1(4)=0.445948490915964 ;VVR2(4)=0.445948490915964 ;VVR3(4)=-0.3399810435848563 ;vvwg(4)=0.223381589678000*0.6521451548625461
		VVR1(5)=0.445948490915964 ;VVR2(5)=0.108103018168071 ;VVR3(5)=-0.3399810435848563 ;vvwg(5)=0.223381589678000*0.6521451548625461
		VVR1(6)=0.108103018168071 ;VVR2(6)=0.445948490915964 ;VVR3(6)=-0.3399810435848563 ;vvwg(6)=0.223381589678000*0.6521451548625461
		VVR1(7)=0.816847572980440 ;VVR2(7)=0.091576213509780 ;VVR3(7)=0.3399810435848563 ;vvwg(7)=0.109951743655333*0.6521451548625461
		VVR1(8)=0.091576213509780 ;VVR2(8)=0.816847572980440 ;VVR3(8)=0.3399810435848563 ;vvwg(8)=0.109951743655333*0.6521451548625461
		VVR1(9)=0.091576213509780 ;VVR2(9)=0.091576213509780 ;VVR3(9)=0.3399810435848563 ;vvwg(9)=0.109951743655333*0.6521451548625461
		VVR1(10)=0.445948490915964 ;VVR2(10)=0.445948490915964 ;VVR3(10)=0.3399810435848563 ;vvwg(10)=0.223381589678000*0.6521451548625461
		VVR1(11)=0.445948490915964 ;VVR2(11)=0.108103018168071 ;VVR3(11)=0.3399810435848563 ;vvwg(11)=0.223381589678000*0.6521451548625461
		VVR1(12)=0.108103018168071 ;VVR2(12)=0.445948490915964 ;VVR3(12)=0.3399810435848563 ;vvwg(12)=0.223381589678000*0.6521451548625461
		VVR1(13)=0.816847572980440 ;VVR2(13)=0.091576213509780 ;VVR3(13)=-0.8611363115940526 ;vvwg(13)=0.109951743655333*0.3478548451374538
		VVR1(14)=0.091576213509780 ;VVR2(14)=0.816847572980440 ;VVR3(14)=-0.8611363115940526 ;vvwg(14)=0.109951743655333*0.3478548451374538
		VVR1(15)=0.091576213509780 ;VVR2(15)=0.091576213509780 ;VVR3(15)=-0.8611363115940526 ;vvwg(15)=0.109951743655333*0.3478548451374538
		VVR1(16)=0.445948490915964 ;VVR2(16)=0.445948490915964 ;VVR3(16)=-0.8611363115940526 ;vvwg(16)=0.223381589678000*0.3478548451374538
		VVR1(17)=0.445948490915964 ;VVR2(17)=0.108103018168071 ;VVR3(17)=-0.8611363115940526 ;vvwg(17)=0.223381589678000*0.3478548451374538
		VVR1(18)=0.108103018168071 ;VVR2(18)=0.445948490915964 ;VVR3(18)=-0.8611363115940526 ;vvwg(18)=0.223381589678000*0.3478548451374538
		VVR1(19)=0.816847572980440 ;VVR2(19)=0.091576213509780 ;VVR3(19)=0.8611363115940526 ;vvwg(19)=0.109951743655333*0.3478548451374538
		VVR1(20)=0.091576213509780 ;VVR2(20)=0.816847572980440 ;VVR3(20)=0.8611363115940526 ;vvwg(20)=0.109951743655333*0.3478548451374538
		VVR1(21)=0.091576213509780 ;VVR2(21)=0.091576213509780 ;VVR3(21)=0.8611363115940526 ;vvwg(21)=0.109951743655333*0.3478548451374538
		VVR1(22)=0.445948490915964 ;VVR2(22)=0.445948490915964 ;VVR3(22)=0.8611363115940526 ;vvwg(22)=0.223381589678000*0.3478548451374538
		VVR1(23)=0.445948490915964 ;VVR2(23)=0.108103018168071 ;VVR3(23)=0.8611363115940526 ;vvwg(23)=0.223381589678000*0.3478548451374538
		VVR1(24)=0.108103018168071 ;VVR2(24)=0.445948490915964 ;VVR3(24)=0.8611363115940526 ;vvwg(24)=0.223381589678000*0.3478548451374538


			
CASE(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  E1=0.2369268850561891
				VVR1(1)=0.888871894660413 ;VVR2(1)=0.055564052669793 ;VVR3(1)=0.0 ;vvwg(1)=0.041955512996649*0.5688888888888889
				VVR1(2)=0.055564052669793 ;VVR2(2)=0.888871894660413 ;VVR3(2)=0.0 ;vvwg(2)=0.041955512996649*0.5688888888888889
				VVR1(3)=0.055564052669793 ;VVR2(3)=0.055564052669793 ;VVR3(3)=0.0 ;vvwg(3)=0.041955512996649*0.5688888888888889
				VVR1(4)=0.295533711735893 ;VVR2(4)=0.634210747745723 ;VVR3(4)=0.0 ;vvwg(4)=0.112098412070887*0.5688888888888889
				VVR1(5)=0.295533711735893 ;VVR2(5)=0.070255540518384 ;VVR3(5)=0.0 ;vvwg(5)=0.112098412070887*0.5688888888888889
				VVR1(6)=0.070255540518384 ;VVR2(6)=0.295533711735893 ;VVR3(6)=0.0 ;vvwg(6)=0.112098412070887*0.5688888888888889
				VVR1(7)=0.634210747745723 ;VVR2(7)=0.295533711735893 ;VVR3(7)=0.0 ;vvwg(7)=0.112098412070887*0.5688888888888889
				VVR1(8)=0.634210747745723 ;VVR2(8)=0.070255540518384 ;VVR3(8)=0.0 ;vvwg(8)=0.112098412070887*0.5688888888888889
				VVR1(9)=0.070255540518384 ;VVR2(9)=0.634210747745723 ;VVR3(9)=0.0 ;vvwg(9)=0.112098412070887*0.5688888888888889
				VVR1(10)=0.333333333333333 ;VVR2(10)=0.333333333333333 ;VVR3(10)=0.0 ;vvwg(10)=0.201542988584730*0.5688888888888889
				VVR1(11)=0.888871894660413 ;VVR2(11)=0.055564052669793 ;VVR3(11)=-0.5384693101056831 ;vvwg(11)=0.041955512996649*0.4786286704993665
				VVR1(12)=0.055564052669793 ;VVR2(12)=0.888871894660413 ;VVR3(12)=-0.5384693101056831 ;vvwg(12)=0.041955512996649*0.4786286704993665
				VVR1(13)=0.055564052669793 ;VVR2(13)=0.055564052669793 ;VVR3(13)=-0.5384693101056831 ;vvwg(13)=0.041955512996649*0.4786286704993665
				VVR1(14)=0.295533711735893 ;VVR2(14)=0.634210747745723 ;VVR3(14)=-0.5384693101056831 ;vvwg(14)=0.112098412070887*0.4786286704993665
				VVR1(15)=0.295533711735893 ;VVR2(15)=0.070255540518384 ;VVR3(15)=-0.5384693101056831 ;vvwg(15)=0.112098412070887*0.4786286704993665
				VVR1(16)=0.070255540518384 ;VVR2(16)=0.295533711735893 ;VVR3(16)=-0.5384693101056831 ;vvwg(16)=0.112098412070887*0.4786286704993665
				VVR1(17)=0.634210747745723 ;VVR2(17)=0.295533711735893 ;VVR3(17)=-0.53846931010568314 ;vvwg(17)=0.112098412070887*0.4786286704993665
				VVR1(18)=0.634210747745723 ;VVR2(18)=0.070255540518384 ;VVR3(18)=-0.5384693101056831 ;vvwg(18)=0.112098412070887*0.4786286704993665
				VVR1(19)=0.070255540518384 ;VVR2(19)=0.634210747745723 ;VVR3(19)=-0.5384693101056831 ;vvwg(19)=0.112098412070887*0.4786286704993665
				VVR1(20)=0.333333333333333 ;VVR2(20)=0.333333333333333 ;VVR3(20)=-0.5384693101056831 ;vvwg(20)=0.201542988584730*0.4786286704993665
				VVR1(21)=0.888871894660413 ;VVR2(21)=0.055564052669793 ;VVR3(21)=0.5384693101056831 ;vvwg(21)=0.041955512996649*0.4786286704993665
				VVR1(22)=0.055564052669793 ;VVR2(22)=0.888871894660413 ;VVR3(22)=0.5384693101056831 ;vvwg(22)=0.041955512996649*0.4786286704993665
				VVR1(23)=0.055564052669793 ;VVR2(23)=0.055564052669793 ;VVR3(23)=0.5384693101056831 ;vvwg(23)=0.041955512996649*0.4786286704993665
				VVR1(24)=0.295533711735893 ;VVR2(24)=0.634210747745723 ;VVR3(24)=0.5384693101056831;vvwg(24)=0.112098412070887*0.4786286704993665
				VVR1(25)=0.295533711735893 ;VVR2(25)=0.070255540518384 ;VVR3(25)=0.5384693101056831 ;vvwg(25)=0.112098412070887*0.4786286704993665
				VVR1(26)=0.070255540518384 ;VVR2(26)=0.295533711735893 ;VVR3(26)=0.5384693101056831 ;vvwg(26)=0.112098412070887*0.4786286704993665
				VVR1(27)=0.634210747745723 ;VVR2(27)=0.295533711735893 ;VVR3(27)=0.5384693101056831 ;vvwg(27)=0.112098412070887*0.4786286704993665
				VVR1(28)=0.634210747745723 ;VVR2(28)=0.070255540518384 ;VVR3(28)=0.5384693101056831 ;vvwg(28)=0.112098412070887*0.4786286704993665
				VVR1(29)=0.070255540518384 ;VVR2(29)=0.634210747745723 ;VVR3(29)=0.5384693101056831 ;vvwg(29)=0.112098412070887*0.4786286704993665
				VVR1(30)=0.333333333333333 ;VVR2(30)=0.333333333333333 ;VVR3(30)=0.5384693101056831 ;vvwg(30)=0.201542988584730*0.4786286704993665
				VVR1(31)=0.888871894660413 ;VVR2(31)=0.055564052669793 ;VVR3(31)=-0.9061798459386640 ;vvwg(31)=0.041955512996649*0.2369268850561891
				VVR1(32)=0.055564052669793 ;VVR2(32)=0.888871894660413 ;VVR3(32)=-0.9061798459386640 ;vvwg(32)=0.041955512996649*0.2369268850561891
				VVR1(33)=0.055564052669793 ;VVR2(33)=0.055564052669793 ;VVR3(33)=-0.9061798459386640 ;vvwg(33)=0.041955512996649*0.2369268850561891
				VVR1(34)=0.295533711735893 ;VVR2(34)=0.634210747745723 ;VVR3(34)=-0.9061798459386640 ;vvwg(34)=0.112098412070887*0.2369268850561891
				VVR1(35)=0.295533711735893 ;VVR2(35)=0.070255540518384 ;VVR3(35)=-0.9061798459386640 ;vvwg(35)=0.112098412070887*0.2369268850561891
				VVR1(36)=0.070255540518384 ;VVR2(36)=0.295533711735893 ;VVR3(36)=-0.9061798459386640 ;vvwg(36)=0.112098412070887*0.2369268850561891
				VVR1(37)=0.634210747745723 ;VVR2(37)=0.295533711735893 ;VVR3(37)=-0.9061798459386640 ;vvwg(37)=0.112098412070887*0.2369268850561891
				VVR1(38)=0.634210747745723 ;VVR2(38)=0.070255540518384 ;VVR3(38)=-0.9061798459386640 ;vvwg(38)=0.112098412070887*0.2369268850561891
				VVR1(39)=0.070255540518384 ;VVR2(39)=0.634210747745723 ;VVR3(39)=-0.9061798459386640 ;vvwg(39)=0.112098412070887*0.2369268850561891
				VVR1(40)=0.333333333333333 ;VVR2(40)=0.333333333333333 ;VVR3(40)=-0.9061798459386640 ;vvwg(40)=0.201542988584730*0.2369268850561891
				VVR1(41)=0.888871894660413 ;VVR2(41)=0.055564052669793 ;VVR3(41)=0.9061798459386640 ;vvwg(41)=0.041955512996649*0.2369268850561891
				VVR1(42)=0.055564052669793 ;VVR2(42)=0.888871894660413 ;VVR3(42)=0.9061798459386640;vvwg(42)=0.041955512996649*0.2369268850561891
				VVR1(43)=0.055564052669793 ;VVR2(43)=0.055564052669793 ;VVR3(43)=0.9061798459386640 ;vvwg(43)=0.041955512996649*0.2369268850561891
				VVR1(44)=0.295533711735893 ;VVR2(44)=0.634210747745723 ;VVR3(44)=0.9061798459386640;vvwg(44)=0.112098412070887*0.2369268850561891
				VVR1(45)=0.295533711735893 ;VVR2(45)=0.070255540518384 ;VVR3(45)=0.9061798459386640 ;vvwg(45)=0.112098412070887*0.2369268850561891
				VVR1(46)=0.070255540518384 ;VVR2(46)=0.295533711735893 ;VVR3(46)=0.9061798459386640 ;vvwg(46)=0.112098412070887*0.2369268850561891
				VVR1(47)=0.634210747745723 ;VVR2(47)=0.295533711735893 ;VVR3(47)=0.9061798459386640 ;vvwg(47)=0.112098412070887*0.2369268850561891
				VVR1(48)=0.634210747745723 ;VVR2(48)=0.070255540518384 ;VVR3(48)=0.9061798459386640 ;vvwg(48)=0.112098412070887*0.2369268850561891
				VVR1(49)=0.070255540518384 ;VVR2(49)=0.634210747745723 ;VVR3(49)=0.9061798459386640 ;vvwg(49)=0.112098412070887*0.2369268850561891
				VVR1(50)=0.333333333333333 ;VVR2(50)=0.333333333333333 ;VVR3(50)=0.9061798459386640 ;vvwg(50)=0.201542988584730*0.2369268850561891
	
			
CASE(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  F=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  F1=0.1713244923791704
 
			      VVR1(1)=0.888871894660413 ;VVR2(1)=0.055564052669793 ;VVR3(1)=0.6612093864662645 ;vvwg(1)=0.041955512996649*0.3607615730481386
				VVR1(2)=0.055564052669793 ;VVR2(2)=0.888871894660413 ;VVR3(2)=0.6612093864662645 ;vvwg(2)=0.041955512996649*0.3607615730481386
				VVR1(3)=0.055564052669793 ;VVR2(3)=0.055564052669793 ;VVR3(3)=0.6612093864662645 ;vvwg(3)=0.041955512996649*0.3607615730481386
				VVR1(4)=0.295533711735893 ;VVR2(4)=0.634210747745723 ;VVR3(4)=0.6612093864662645 ;vvwg(4)=0.112098412070887*0.3607615730481386
				VVR1(5)=0.295533711735893 ;VVR2(5)=0.070255540518384 ;VVR3(5)=0.6612093864662645 ;vvwg(5)=0.112098412070887*0.3607615730481386
				VVR1(6)=0.070255540518384 ;VVR2(6)=0.295533711735893 ;VVR3(6)=0.6612093864662645 ;vvwg(6)=0.112098412070887*0.3607615730481386
				VVR1(7)=0.634210747745723 ;VVR2(7)=0.295533711735893 ;VVR3(7)=0.6612093864662645 ;vvwg(7)=0.112098412070887*0.3607615730481386
				VVR1(8)=0.634210747745723 ;VVR2(8)=0.070255540518384 ;VVR3(8)=0.6612093864662645 ;vvwg(8)=0.112098412070887*0.3607615730481386
				VVR1(9)=0.070255540518384 ;VVR2(9)=0.634210747745723 ;VVR3(9)=0.6612093864662645 ;vvwg(9)=0.112098412070887*0.3607615730481386
				VVR1(10)=0.333333333333333 ;VVR2(10)=0.333333333333333 ;VVR3(10)=0.6612093864662645 ;vvwg(10)=0.201542988584730*0.3607615730481386
				VVR1(11)=0.888871894660413 ;VVR2(11)=0.055564052669793 ;VVR3(11)=-0.6612093864662645 ;vvwg(11)=0.041955512996649*0.3607615730481386
				VVR1(12)=0.055564052669793 ;VVR2(12)=0.888871894660413 ;VVR3(12)=-0.6612093864662645 ;vvwg(12)=0.041955512996649*0.3607615730481386
				VVR1(13)=0.055564052669793 ;VVR2(13)=0.055564052669793 ;VVR3(13)=-0.6612093864662645 ;vvwg(13)=0.041955512996649*0.3607615730481386
				VVR1(14)=0.295533711735893 ;VVR2(14)=0.634210747745723 ;VVR3(14)=-0.6612093864662645 ;vvwg(14)=0.112098412070887*0.3607615730481386
				VVR1(15)=0.295533711735893 ;VVR2(15)=0.070255540518384 ;VVR3(15)=-0.6612093864662645 ;vvwg(15)=0.112098412070887*0.3607615730481386
				VVR1(16)=0.070255540518384 ;VVR2(16)=0.295533711735893 ;VVR3(16)=-0.6612093864662645;vvwg(16)=0.112098412070887*0.3607615730481386
				VVR1(17)=0.634210747745723 ;VVR2(17)=0.295533711735893 ;VVR3(17)=-0.6612093864662645 ;vvwg(17)=0.112098412070887*0.3607615730481386
				VVR1(18)=0.634210747745723 ;VVR2(18)=0.070255540518384 ;VVR3(18)=-0.6612093864662645;vvwg(18)=0.112098412070887*0.3607615730481386
				VVR1(19)=0.070255540518384 ;VVR2(19)=0.634210747745723 ;VVR3(19)=-0.6612093864662645 ;vvwg(19)=0.112098412070887*0.3607615730481386
				VVR1(20)=0.333333333333333 ;VVR2(20)=0.333333333333333 ;VVR3(20)=-0.6612093864662645 ;vvwg(20)=0.201542988584730*0.3607615730481386
				VVR1(21)=0.888871894660413 ;VVR2(21)=0.055564052669793 ;VVR3(21)=-0.2386191860831969 ;vvwg(21)=0.041955512996649*0.4679139345726910
				VVR1(22)=0.055564052669793 ;VVR2(22)=0.888871894660413 ;VVR3(22)=-0.2386191860831969 ;vvwg(22)=0.041955512996649*0.4679139345726910
				VVR1(23)=0.055564052669793 ;VVR2(23)=0.055564052669793 ;VVR3(23)=-0.2386191860831969;vvwg(23)=0.041955512996649*0.4679139345726910
				VVR1(24)=0.295533711735893 ;VVR2(24)=0.634210747745723 ;VVR3(24)=-0.2386191860831969;vvwg(24)=0.112098412070887*0.4679139345726910
				VVR1(25)=0.295533711735893 ;VVR2(25)=0.070255540518384 ;VVR3(25)=-0.2386191860831969 ;vvwg(25)=0.112098412070887*0.4679139345726910
				VVR1(26)=0.070255540518384 ;VVR2(26)=0.295533711735893 ;VVR3(26)=-0.2386191860831969 ;vvwg(26)=0.112098412070887*0.4679139345726910
				VVR1(27)=0.634210747745723 ;VVR2(27)=0.295533711735893 ;VVR3(27)=-0.2386191860831969 ;vvwg(27)=0.112098412070887*0.4679139345726910
				VVR1(28)=0.634210747745723 ;VVR2(28)=0.070255540518384 ;VVR3(28)=-0.2386191860831969 ;vvwg(28)=0.112098412070887*0.4679139345726910
				VVR1(29)=0.070255540518384 ;VVR2(29)=0.634210747745723 ;VVR3(29)=-0.2386191860831969 ;vvwg(29)=0.112098412070887*0.4679139345726910
				VVR1(30)=0.333333333333333 ;VVR2(30)=0.333333333333333 ;VVR3(30)=-0.2386191860831969 ;vvwg(30)=0.201542988584730*0.4679139345726910
				VVR1(31)=0.888871894660413 ;VVR2(31)=0.055564052669793 ;VVR3(31)=0.2386191860831969 ;vvwg(31)=0.041955512996649*0.4679139345726910
				VVR1(32)=0.055564052669793 ;VVR2(32)=0.888871894660413 ;VVR3(32)=0.2386191860831969 ;vvwg(32)=0.041955512996649*0.4679139345726910
				VVR1(33)=0.055564052669793 ;VVR2(33)=0.055564052669793 ;VVR3(33)=0.2386191860831969 ;vvwg(33)=0.041955512996649*0.4679139345726910
				VVR1(34)=0.295533711735893 ;VVR2(34)=0.634210747745723 ;VVR3(34)=0.2386191860831969 ;vvwg(34)=0.112098412070887*0.4679139345726910
				VVR1(35)=0.295533711735893 ;VVR2(35)=0.070255540518384 ;VVR3(35)=0.2386191860831969 ;vvwg(35)=0.112098412070887*0.4679139345726910
				VVR1(36)=0.070255540518384 ;VVR2(36)=0.295533711735893 ;VVR3(36)=0.2386191860831969 ;vvwg(36)=0.112098412070887*0.4679139345726910
				VVR1(37)=0.634210747745723 ;VVR2(37)=0.295533711735893 ;VVR3(37)=0.2386191860831969 ;vvwg(37)=0.112098412070887*0.4679139345726910
				VVR1(38)=0.634210747745723 ;VVR2(38)=0.070255540518384 ;VVR3(38)=0.2386191860831969 ;vvwg(38)=0.112098412070887*0.4679139345726910
				VVR1(39)=0.070255540518384 ;VVR2(39)=0.634210747745723 ;VVR3(39)=0.2386191860831969 ;vvwg(39)=0.112098412070887*0.4679139345726910
				VVR1(40)=0.333333333333333 ;VVR2(40)=0.333333333333333 ;VVR3(40)=0.2386191860831969 ;vvwg(40)=0.201542988584730*0.4679139345726910
				VVR1(41)=0.888871894660413 ;VVR2(41)=0.055564052669793 ;VVR3(41)=-0.9324695142031521 ;vvwg(41)=0.041955512996649*0.1713244923791704
				VVR1(42)=0.055564052669793 ;VVR2(42)=0.888871894660413 ;VVR3(42)=-0.9324695142031521;vvwg(42)=0.041955512996649*0.1713244923791704
				VVR1(43)=0.055564052669793 ;VVR2(43)=0.055564052669793 ;VVR3(43)=-0.9324695142031521 ;vvwg(43)=0.041955512996649*0.1713244923791704
				VVR1(44)=0.295533711735893 ;VVR2(44)=0.634210747745723 ;VVR3(44)=-0.9324695142031521;vvwg(44)=0.112098412070887*0.1713244923791704
				VVR1(45)=0.295533711735893 ;VVR2(45)=0.070255540518384 ;VVR3(45)=-0.9324695142031521 ;vvwg(45)=0.112098412070887*0.1713244923791704
				VVR1(46)=0.070255540518384 ;VVR2(46)=0.295533711735893 ;VVR3(46)=-0.9324695142031521 ;vvwg(46)=0.112098412070887*0.1713244923791704
				VVR1(47)=0.634210747745723 ;VVR2(47)=0.295533711735893 ;VVR3(47)=-0.9324695142031521 ;vvwg(47)=0.112098412070887*0.1713244923791704
				VVR1(48)=0.634210747745723 ;VVR2(48)=0.070255540518384 ;VVR3(48)=-0.9324695142031521 ;vvwg(48)=0.112098412070887*0.1713244923791704
				VVR1(49)=0.070255540518384 ;VVR2(49)=0.634210747745723 ;VVR3(49)=-0.9324695142031521 ;vvwg(49)=0.112098412070887*0.1713244923791704
				VVR1(50)=0.333333333333333 ;VVR2(50)=0.333333333333333 ;VVR3(50)=-0.9324695142031521 ;vvwg(50)=0.201542988584730*0.1713244923791704
				VVR1(51)=0.888871894660413 ;VVR2(51)=0.055564052669793 ;VVR3(51)=0.9324695142031521 ;vvwg(51)=0.041955512996649*0.1713244923791704
				VVR1(52)=0.055564052669793 ;VVR2(52)=0.888871894660413 ;VVR3(52)=0.9324695142031521;vvwg(52)=0.041955512996649*0.1713244923791704
				VVR1(53)=0.055564052669793 ;VVR2(53)=0.055564052669793 ;VVR3(53)=0.9324695142031521 ;vvwg(53)=0.041955512996649*0.1713244923791704
				VVR1(54)=0.295533711735893 ;VVR2(54)=0.634210747745723 ;VVR3(54)=0.9324695142031521;vvwg(54)=0.112098412070887*0.1713244923791704
				VVR1(55)=0.295533711735893 ;VVR2(55)=0.070255540518384 ;VVR3(55)=0.9324695142031521 ;vvwg(55)=0.112098412070887*0.1713244923791704
				VVR1(56)=0.070255540518384 ;VVR2(56)=0.295533711735893 ;VVR3(56)=0.9324695142031521 ;vvwg(56)=0.112098412070887*0.1713244923791704
				VVR1(57)=0.634210747745723 ;VVR2(57)=0.295533711735893 ;VVR3(57)=0.9324695142031521 ;vvwg(57)=0.112098412070887*0.1713244923791704
				VVR1(58)=0.634210747745723 ;VVR2(58)=0.070255540518384 ;VVR3(58)=0.9324695142031521 ;vvwg(58)=0.112098412070887*0.1713244923791704
				VVR1(59)=0.070255540518384 ;VVR2(59)=0.634210747745723 ;VVR3(59)=0.9324695142031521 ;vvwg(59)=0.112098412070887*0.1713244923791704
				VVR1(60)=0.333333333333333 ;VVR2(60)=0.333333333333333 ;VVR3(60)=0.9324695142031521 ;vvwg(60)=0.201542988584730*0.1713244923791704




	
			
END SELECT
		QPOINTS(:,:)=0.0d0
! 		
! 		  WEQUA3D(:)=vvwg(:)*0.5d0
		do kk=1,qp_prism
			WEQUA3D(kk)=vvwg(kk)*0.5d0
			
			R=VVR1(kk); S=VVR2(kk); T=VVR3(kk)
			VVnxi(1)=(0.5d0)*r*(1.0d0-t)
			VVnxi(2)=(0.5d0)*(s)*(1.0d0-t)
			VVnxi(3)=(0.5d0)*(1.0-R-s)*(1.0d0-t)
			VVnxi(4)=(0.5d0)*r*(1.0d0+t)
			VVnxi(5)=(0.5d0)*(s)*(1.0d0+t)
			VVnxi(6)=(0.5d0)*(1.0-R-s)*(1.0d0+t)
			
			DO J=1,6
			QPOINTS(:,kk)=QPOINTS(:,kk)+(VVNXI(j)*VEXT(j,:))
			END DO
		    
! 			

		END DO
! 		



END SUBROUTINE QUADRATUREPRISM


SUBROUTINE QUADRATUREPYRA(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for a pyramid
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
REAL::R,S,T,a,b,c,d,e,f,g
REAL::a1,b1,c1,d1,e1,f1,sumwe
INTEGER::Kk,J,ii,ij,ik,count1,alls

alls=QP_PYRA


 
 WEQUA3D=0.0d0
  QPOINTS=0.0d0
sumwe=0.0d0

SELECT CASE(IGQRULES)
 
  case(1)
		VVR1(1)=0.0d0 ;VVR2(1)=0.0d0 ;VVR3(1)=-0.5d0
		
		  vvwg(1)=8.0d0

 
  	


 

 
 case(2)
		VVR1(1)=-0.584237394672177188;VVR2(1)=-0.58423739467217718 ;VVR3(1)=-0.6666666666666666;vvwg(1)=0.81
		VVR1(2)=0.58423739467217718 ;VVR2(2)=-0.58423739467217718 ;VVR3(2)=-0.6666666666666666;vvwg(2)=0.81
		VVR1(3)=0.58423739467217718 ;VVR2(3)=0.58423739467217718 ;VVR3(3)=-0.6666666666666666;vvwg(3)=0.81
		VVR1(4)=-0.58423739467217718 ;VVR2(4)=0.58423739467217718 ;VVR3(4)=-0.6666666666666666;vvwg(4)=0.81
		VVR1(5)=0.0d0 ;VVR2(5)=0.0d0 ;VVR3(5)=0.4d0;vvwg(5)=3.76d0
		
 

 
  	


  	
 CASE(3)
  a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006d0*0.515003019323671498
  b1=1.104848006d0*0.2571837452420646589
  c1=1.104848006d0*2.474004977113405936
  d1=1.104848006d0*0.419515737191525950
 
 


		VVR1(1)=-a ;VVR2(1)=-a ;VVR3(1)=d ;vvwg(1)=a1
		VVR1(2)=a ;VVR2(2)=-a ;VVR3(2)=d ;vvwg(2)=a1
		VVR1(3)=a ;VVR2(3)=a ;VVR3(3)=d ;vvwg(3)=a1
		VVR1(4)=-a ;VVR2(4)=a ;VVR3(4)=d ;vvwg(4)=a1
		VVR1(5)=-b ;VVR2(5)=0.0 ;VVR3(5)=e ;vvwg(5)=b1
		VVR1(6)=b ;VVR2(6)=0.0 ;VVR3(6)=e ;vvwg(6)=b1
		VVR1(7)=0.0 ;VVR2(7)=-b ;VVR3(7)=e ;vvwg(7)=b1
		VVR1(8)=0.0 ;VVR2(8)=b ;VVR3(8)=e ;vvwg(8)=b1
		VVR1(9)=0.0 ;VVR2(9)=0.0 ;VVR3(9)=f ;vvwg(9)=c1
		VVR1(10)=-c ;VVR2(10)=-c ;VVR3(10)=g ;vvwg(10)=d1
		VVR1(11)=c ;VVR2(11)=-c ;VVR3(11)=g ;vvwg(11)=d1
		VVR1(12)=c ;VVR2(12)=c ;VVR3(12)=g ;vvwg(12)=d1
		VVR1(13)=-c ;VVR2(13)=c ;VVR3(13)=g ;vvwg(13)=d1










			
 CASE(4)
  a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006*0.515003019323671498
  b1=1.104848006*0.2571837452420646589
  c1=1.104848006*2.474004977113405936
  d1=1.104848006*0.419515737191525950
 
 


		VVR1(1)=-a ;VVR2(1)=-a ;VVR3(1)=d ;vvwg(1)=a1
		VVR1(2)=a ;VVR2(2)=-a ;VVR3(2)=d ;vvwg(2)=a1
		VVR1(3)=a ;VVR2(3)=a ;VVR3(3)=d ;vvwg(3)=a1
		VVR1(4)=-a ;VVR2(4)=a ;VVR3(4)=d ;vvwg(4)=a1
		VVR1(5)=-b ;VVR2(5)=0.0 ;VVR3(5)=e ;vvwg(5)=b1
		VVR1(6)=b ;VVR2(6)=0.0 ;VVR3(6)=e ;vvwg(6)=b1
		VVR1(7)=0.0 ;VVR2(7)=-b ;VVR3(7)=e ;vvwg(7)=b1
		VVR1(8)=0.0 ;VVR2(8)=b ;VVR3(8)=e ;vvwg(8)=b1
		VVR1(9)=0.0 ;VVR2(9)=0.0 ;VVR3(9)=f ;vvwg(9)=c1
		VVR1(10)=-c ;VVR2(10)=-c ;VVR3(10)=g ;vvwg(10)=d1
		VVR1(11)=c ;VVR2(11)=-c ;VVR3(11)=g ;vvwg(11)=d1
		VVR1(12)=c ;VVR2(12)=c ;VVR3(12)=g ;vvwg(12)=d1
		VVR1(13)=-c ;VVR2(13)=c ;VVR3(13)=g ;vvwg(13)=d1


			
CASE(5)

 a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006*0.515003019323671498
  b1=1.104848006*0.2571837452420646589
  c1=1.104848006*2.474004977113405936
  d1=1.104848006*0.419515737191525950
 
 


		VVR1(1)=-a ;VVR2(1)=-a ;VVR3(1)=d ;vvwg(1)=a1
		VVR1(2)=a ;VVR2(2)=-a ;VVR3(2)=d ;vvwg(2)=a1
		VVR1(3)=a ;VVR2(3)=a ;VVR3(3)=d ;vvwg(3)=a1
		VVR1(4)=-a ;VVR2(4)=a ;VVR3(4)=d ;vvwg(4)=a1
		VVR1(5)=-b ;VVR2(5)=0.0 ;VVR3(5)=e ;vvwg(5)=b1
		VVR1(6)=b ;VVR2(6)=0.0 ;VVR3(6)=e ;vvwg(6)=b1
		VVR1(7)=0.0 ;VVR2(7)=-b ;VVR3(7)=e ;vvwg(7)=b1
		VVR1(8)=0.0 ;VVR2(8)=b ;VVR3(8)=e ;vvwg(8)=b1
		VVR1(9)=0.0 ;VVR2(9)=0.0 ;VVR3(9)=f ;vvwg(9)=c1
		VVR1(10)=-c ;VVR2(10)=-c ;VVR3(10)=g ;vvwg(10)=d1
		VVR1(11)=c ;VVR2(11)=-c ;VVR3(11)=g ;vvwg(11)=d1
		VVR1(12)=c ;VVR2(12)=c ;VVR3(12)=g ;vvwg(12)=d1
		VVR1(13)=-c ;VVR2(13)=c ;VVR3(13)=g ;vvwg(13)=d1
	
			
CASE(6,7,8,9)

 a=0.673931986207731726
  b=0.610639618865075532
  c=0.580939660561084423
  d=-0.1428571428571428571
  e=-0.321428571428571429
  f=0.524394036075370072
  g=-0.830065359477124183
  a1=1.104848006*0.515003019323671498
  b1=1.104848006*0.2571837452420646589
  c1=1.104848006*2.474004977113405936
  d1=1.104848006*0.419515737191525950
 
 


		VVR1(1)=-a ;VVR2(1)=-a ;VVR3(1)=d ;vvwg(1)=a1
		VVR1(2)=a ;VVR2(2)=-a ;VVR3(2)=d ;vvwg(2)=a1
		VVR1(3)=a ;VVR2(3)=a ;VVR3(3)=d ;vvwg(3)=a1
		VVR1(4)=-a ;VVR2(4)=a ;VVR3(4)=d ;vvwg(4)=a1
		VVR1(5)=-b ;VVR2(5)=0.0 ;VVR3(5)=e ;vvwg(5)=b1
		VVR1(6)=b ;VVR2(6)=0.0 ;VVR3(6)=e ;vvwg(6)=b1
		VVR1(7)=0.0 ;VVR2(7)=-b ;VVR3(7)=e ;vvwg(7)=b1
		VVR1(8)=0.0 ;VVR2(8)=b ;VVR3(8)=e ;vvwg(8)=b1
		VVR1(9)=0.0 ;VVR2(9)=0.0 ;VVR3(9)=f ;vvwg(9)=c1
		VVR1(10)=-c ;VVR2(10)=-c ;VVR3(10)=g ;vvwg(10)=d1
		VVR1(11)=c ;VVR2(11)=-c ;VVR3(11)=g ;vvwg(11)=d1
		VVR1(12)=c ;VVR2(12)=c ;VVR3(12)=g ;vvwg(12)=d1
		VVR1(13)=-c ;VVR2(13)=c ;VVR3(13)=g ;vvwg(13)=d1



	
			
END SELECT
		QPOINTS(:,:)=0.0d0

! 		  WEQUA3D(:)=vvwg(:)*0.1250000000000
		do kk=1,qp_pyra
			WEQUA3D(kk)=vvwg(kk)*0.1250000000000
			
			R=VVR1(kk); S=VVR2(kk); T=VVR3(kk)
			VVnxi(1)=(0.1250000000000)*(1.0-R)*(1.0d0-s)*(1.0d0-t)
			VVnxi(2)=(0.1250000000000)*(1.0+R)*(1.0d0-s)*(1.0d0-t)
			VVnxi(3)=(0.1250000000000)*(1.0+R)*(1.0d0+s)*(1.0d0-t)
			VVnxi(4)=(0.1250000000000)*(1.0-R)*(1.0d0+s)*(1.0d0-t)
			VVnxi(5)=0.5d0*(1.0d0+t)
			
			
			DO J=1,5
			QPOINTS(:,kk)=QPOINTS(:,kk)+(VVNXI(j)*VEXT(j,:))
			END DO
		    
			

		END DO




END SUBROUTINE QUADRATUREPYRA


SUBROUTINE QUADRATUREHEXA(N,IGQRULES)
 !> @brief
!> This subroutine computes the quadrature points and weights for a hexahedral
IMPLICIT NONE
INTEGER,INTENT(IN)::IGQRULES,N
REAL::R,S,T,a,b,c,d,e,f
REAL::a1,b1,c1,d1,e1,f1
INTEGER::Kk,J,ii,ij,ik,count1,alls

 WEQUA3D=0.0d0
  QPOINTS=0.0d0

SELECT CASE(IGQRULES)
 

 case(1)

		vvwg(1) = 8.0d0
	    VVR1(1)=0.0d0	;VVR2(1)=0.0d0	;VVR3(1)=0.0d0


 case(2)

 a=-0.5773502691896257
  b=0.5773502691896257
  a1=1.0d0
  b1=1.0d0
  
  vvnpox(1)=a	;vvnpox(2)=b	
  vvnpoy(1)= a	;vvnpoy(2)=b	
  vvnpoz(1)= a	;vvnpoz(2)=b

  vvwpox(1)=a1	;vvwpox(2)=b1	
  vvwpoy(1)= a1	;vvwpoy(2)=b1	
  vvwpoz(1)= a1	;vvwpoz(2)=b1
  
  count1=0
 do ii=1,2
    do ij=1,2
      do ik=1,2
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) ;VVR3(count1)=vvnpoz(ik)   
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do
  	


  	
 CASE(3)
  a=0.0d0
  b=-0.7745966692414834
  c=0.7745966692414834
  a1=0.8888888888888888
  b1=0.5555555555555556
  c1=0.5555555555555556
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1
  count1=0
 do ii=1,3
    do ij=1,3
      do ik=1,3
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) ;VVR3(count1)=vvnpoz(ik)   
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do
 

			
 CASE(4)
  a=-0.3399810435848563
  b=0.3399810435848563
  c=-0.8611363115940526
  d=0.8611363115940526
   a1=0.6521451548625461
  b1=0.6521451548625461
  c1=0.3478548451374538
  d1=0.3478548451374538
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1
  count1=0
 do ii=1,4
    do ij=1,4
      do ik=1,4
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) ;VVR3(count1)=vvnpoz(ik)  
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do


			
CASE(5)

  a=0.0d0
  b=-0.5384693101056831
  c=0.5384693101056831
  d=-0.9061798459386640
  e=0.9061798459386640
    a1=0.5688888888888889
  b1=0.4786286704993665
  c1=0.4786286704993665
  d1=0.2369268850561891
  E1=0.2369268850561891
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1
  count1=0
 do ii=1,5
    do ij=1,5
      do ik=1,5
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) ;VVR3(count1)=vvnpoz(ik) 
	vvwg(count1)=(vvwpox(ii)*vvwpox(ij)*vvwpox(ik))
      end do
    end do
end do
	
			
CASE(6,7,8,9)

  a=0.6612093864662645
  b=-0.6612093864662645
  c=-0.2386191860831969
  d=0.2386191860831969
  e=-0.9324695142031521
  F=0.9324695142031521
  a1=0.3607615730481386
  b1=0.3607615730481386
  c1=0.4679139345726910
  d1=0.4679139345726910
  e1=0.1713244923791704
  F1=0.1713244923791704
  vvnpox(1)=a	;vvnpox(2)=b	;vvnpox(3)=c ;vvnpox(4)=d ;vvnpox(5)=e;vvnpox(6)=f
  vvnpoy(1)= a	;vvnpoy(2)=b	;vvnpoy(3)=c ;vvnpoy(4)=d ;vvnpoy(5)=e;vvnpoy(6)=f
  vvnpoz(1)= a	;vvnpoz(2)=b	;vvnpoz(3)=c ;vvnpoz(4)=d ;vvnpoz(5)=e;vvnpoz(6)=f
  vvwpox(1)=a1	;vvwpox(2)=b1	;vvwpox(3)=c1 ;vvwpox(4)=d1 ;vvwpox(5)=e1 ;vvwpox(6)=f1
  vvwpoy(1)= a1	;vvwpoy(2)=b1	;vvwpoy(3)=c1 ;vvwpoy(4)=d1 ;vvwpoy(5)=e1 ;vvwpoy(6)=f1
  vvwpoz(1)= a1	;vvwpoz(2)=b1	;vvwpoz(3)=c1 ;vvwpoz(4)=d1 ;vvwpoz(5)=e1 ;vvwpoz(6)=f1
  count1=0
 do ii=1,6
    do ij=1,6
      do ik=1,6
	count1=count1+1
	VVR1(count1)=vvnpox(ii);VVR2(count1)=vvnpoy(ij) ;VVR3(count1)=vvnpoz(ik) 
	vvwg(count1)=vvwpox(ii)*vvwpox(ij)*vvwpox(ik)
      end do
    end do
end do
	
			
END SELECT
		QPOINTS(:,:)=0.0d0
		
		  vvwg(:)=vvwg(:)*0.125d0
! 		  WEQUA3D(:)=vvwg(:)
		do kk=1,qp_hexa
			WEQUA3D(kk)=vvwg(kk)
			R=VVR1(kk); S=VVR2(kk); T=VVR3(kk)
			VVnxi(1)=(0.125d0)*(1.0d0-R)*(1.0d0-s)*(1.0d0-t)
			VVnxi(2)=(0.125d0)*(1.0d0+R)*(1.0d0-s)*(1.0d0-t)
			VVnxi(3)=(0.125d0)*(1.0d0+R)*(1.0d0+s)*(1.0d0-t)
			VVnxi(4)=(0.125d0)*(1.0d0-R)*(1.0d0+s)*(1.0d0-t)
			VVnxi(5)=(0.125d0)*(1.0d0-R)*(1.0d0-s)*(1.0d0+t)
			VVnxi(6)=(0.125d0)*(1.0d0+R)*(1.0d0-s)*(1.0d0+t)
			VVnxi(7)=(0.125d0)*(1.0d0+R)*(1.0d0+s)*(1.0d0+t)
			VVnxi(8)=(0.125d0)*(1.0d0-R)*(1.0d0+s)*(1.0d0+t)
			DO J=1,8
			QPOINTS(:,kk)=QPOINTS(:,kk)+(VVNXI(j)*VEXT(j,:))
			END DO
			

		END DO




 

END SUBROUTINE QUADRATUREHEXA

SUBROUTINE ROTATEF(N,TRI,ROTVECT,VECTCO,ANGLE1,ANGLE2)
 !> @brief
!> This subroutine rotates the vector of fluxes in the directions normal to the face in 3D 
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::TRI
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVECT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VECTCO
REAL,INTENT(IN)::ANGLE1,ANGLE2
REAL::sia1,coa1,coa2,sia2


!BUILD MATRIX OF ROTATION!
 coa1=COS(ANGLE1)
 sia1=sin(angle1)
 coa2=cos(angle2)
 sia2=sin(angle2)
 
tri=zero

TRI(1,1)=1.0d0

TRI(2,2)=coa1*sia2!COS(ANGLE1)*SIN(ANGLE2)
TRI(2,3)=sia1*sia2!SIN(ANGLE1)*SIN(ANGLE2)
TRI(2,4)=coa2!COS(ANGLE2)

TRI(3,2)=coa1*coa2!COS(ANGLE1)*COS(ANGLE2)
TRI(3,3)=sia1*coa2!SIN(ANGLE1)*COS(ANGLE2)
TRI(3,4)=-sia2!-SIN(ANGLE2)

TRI(4,2)=-sia1!-SIN(ANGLE1)
TRI(4,3)=coa1!COS(ANGLE1)
TRI(5,5)=1.0d0







ROTVECT(1:5)=MATMUL(TRI(1:5,1:5),VECTCO(1:5))
IF (MULTISPECIES.EQ.1)THEN
ROTVECT(6:nof_Variables)=VECTCO(6:nof_Variables)
END IF


END SUBROUTINE ROTATEF


SUBROUTINE ROTATEB(N,INVTRI,ROTVECT,VECTCO,ANGLE1,ANGLE2)
 !> @brief
!> This subroutine rotates back the vector of fluxes from the directions normal to the face to cartesian coordinates
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INVTRI
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVECT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VECTCO
REAL,INTENT(IN)::ANGLE1,ANGLE2
REAL::sia1,coa1,coa2,sia2

!BUILD MATRIX OF ROTATION!
 coa1=COS(ANGLE1)
 sia1=sin(angle1)
 coa2=cos(angle2)
 sia2=sin(angle2)




invtri=zero


!BUILD MATRIX OF ROTATION!
INVTRI(1,1)=1.0d0

INVTRI(2,2)=coa1*sia2!COS(ANGLE1)*SIN(ANGLE2)
INVTRI(2,3)=coa1*coa2!COS(ANGLE1)*COS(ANGLE2)
INVTRI(2,4)=-sia1!-SIN(ANGLE1)

INVTRI(3,2)=sia1*sia2!SIN(ANGLE1)*SIN(ANGLE2)
INVTRI(3,3)=sia1*coa2!SIN(ANGLE1)*COS(ANGLE2)
INVTRI(3,4)=coa1!COS(ANGLE1)

INVTRI(4,2)=coa2!COS(ANGLE2)
INVTRI(4,3)=-sia2!-SIN(ANGLE2)
INVTRI(5,5)=1.0d0

ROTVECT(1:5)=MATMUL(INVTRI(1:5,1:5),VECTCO(1:5))
IF (MULTISPECIES.EQ.1)THEN
ROTVECT(6:nof_Variables)=VECTCO(6:nof_Variables)
END IF


END SUBROUTINE ROTATEB


SUBROUTINE ROTATEF2d(N,TRI,ROTVECT,VECTCO,ANGLE1,ANGLE2)
 !> @brief
!> This subroutine rotates the vector of fluxes in the directions normal to the edge in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::TRI
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVECT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VECTCO
REAL,INTENT(IN)::ANGLE1,ANGLE2

ROTVECT(1)=VECTCO(1)
ROTVECT(2)=(ANGLE1*VECTCO(2))+(ANGLE2*VECTCO(3))
ROTVECT(3)=-(ANGLE2*VECTCO(2))+(ANGLE1*VECTCO(3))
ROTVECT(4)=VECTCO(4)

IF (MULTISPECIES.EQ.1)THEN
ROTVECT(5:nof_Variables)=VECTCO(5:nof_Variables)
END IF

END SUBROUTINE ROTATEF2d


SUBROUTINE ROTATEB2d(N,INVTRI,ROTVECT,VECTCO,ANGLE1,ANGLE2)
 !> @brief
!> This subroutine rotates back the vector of fluxes from the directions normal to the edge to cartesian coordinates
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::INVTRI
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ROTVECT
REAL,ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::VECTCO
REAL,INTENT(IN)::ANGLE1,ANGLE2




!BUILD MATRIX OF ROTATION!

ROTVECT(1)=VECTCO(1)
ROTVECT(2)=(ANGLE1*VECTCO(2))-(ANGLE2*VECTCO(3))
ROTVECT(3)=(ANGLE2*VECTCO(2))+(ANGLE1*VECTCO(3))
ROTVECT(4)=VECTCO(4)

IF (MULTISPECIES.EQ.1)THEN
ROTVECT(5:nof_Variables)=VECTCO(5:nof_Variables)
END IF

END SUBROUTINE ROTATEB2d




SUBROUTINE PROBEPOS(N,PROBEI)
 !> @brief
!> This subroutine establishes the cells where the probe positions belong to
IMPLICIT NONE
INTEGER,ALLOCATABLE,DIMENSION(:,:),INTENT(INOUT)::PROBEI
INTEGER,INTENT(IN)::N
INTEGER::I,J,K,L,KMAXE,INV
REAL::dist
REAL::DUMIN,DUMOUT,DELTA
ALLOCATE (PROBEI(N:N,NPROBES))
PROBEI(N:N,:)=0

KMAXE=XMPIELRANK(N)

	if (dimensiona.eq.3)then
	DO INV=1,NPROBES
	DELTA=TOLBIG
	DO I=1,KMAXE
		iconsi=i
		
		
		
		call COMPUTE_CENTRE3d(N,Iconsi)
		vext(1,1:3)=cords(1:3)
		vext(2,1:3)=PROBEC(inv,1:3)
		dist=distance3(n)
		
		
		IF (dist.LT.DELTA) THEN
			DELTA=dist
			L=I
		END IF
	END DO
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	DUMOUT=DELTA
		DUMIN=0.0d0
		CALL MPI_ALLREDUCE(DUMOUT,DUMIN,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	IF (abs(DUMIN-DELTA).le.1.0E-15) THEN
	PROBEI(N,INV)=L
	END IF
	END DO

	else
      	DO INV=1,NPROBES
	DELTA=TOLBIG
	DO I=1,KMAXE
		iconsi=i
		
		
		
		call COMPUTE_CENTRE2d(N,Iconsi)
		vext(1,1:2)=cords(1:2)
		vext(2,1:2)=PROBEC(inv,1:2)
		dist=distance2(n)
		
		
		IF (dist.Le.DELTA) THEN
			DELTA=dist
			L=I
		END IF
	END DO
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	
	DUMOUT=DELTA
		DUMIN=0.0d0
		CALL MPI_ALLREDUCE(DUMOUT,DUMIN,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)
		CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
	IF (abs(DUMIN-DELTA).le.1.0E-15) THEN
	PROBEI(N,INV)=L
	
	END IF
	END DO
	end if

	
END SUBROUTINE PROBEPOS


SUBROUTINE  ANGLEX(ANGLEFACEX)
 !> @brief
!> This subroutine computes the angles according to the quadrant sign
IMPLICIT NONE
REAL,INTENT(INOUT)::ANGLEFACEX

IF ((A_rot.NE.zero).AND.(b_rot.NE.zero))THEN
	IF ((A_rot.GT.zero).AND.(b_rot.GT.zero))THEN
		ANGLEFACEX=ATAN(b_rot/A_rot)
	END IF
	IF ((A_rot.LT.zero).AND.(b_rot.GT.zero))THEN
		ANGLEFACEX=PI+(ATAN(b_rot/A_rot))
	END IF
	IF ((A_rot.LT.zero).AND.(b_rot.LT.zero))THEN
		ANGLEFACEX=PI+(ATAN(b_rot/A_rot))
	END IF
	IF ((A_rot.GT.zero).AND.(b_rot.LT.zero))THEN
		ANGLEFACEX=(2.0d0*PI)+(ATAN(b_rot/A_rot))
	END IF
END IF
IF ((A_rot.EQ.zero).AND.(b_rot.NE.zero))THEN
	IF (B_rot.GT.zero) THEN
		ANGLEFACEX=(PI/2.0d0)
	END IF
	IF (B_rot.LT.zero) THEN
		ANGLEFACEX=3.0d0*(PI/2.0d0)
	END IF
END IF
IF ((A_rot.NE.zero).AND.(b_rot.EQ.zero))THEN
	IF (A_rot.GT.zero) THEN
		ANGLEFACEX=zero
	END IF
	IF (A_rot.LT.zero) THEN
		ANGLEFACEX=PI
	END IF
END IF
END SUBROUTINE ANGLEX
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !---------------------------------------------------------------------------------------------!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!FUNCTION TO CALCULATE THE ANGLES BASED ON THE COORDINATES !!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OF THE NORMAL PLANE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE ANGLEY(ANGLEFACEY)
 !> @brief
!> This subroutine computes the angles according to the quadrant sign
IMPLICIT NONE
REAL,INTENT(INOUT)::ANGLEFACEY

IF (C_rot.EQ.zero)THEN
	ANGLEFACEY=ACOS(zero)
END IF
IF (C_rot.NE.zero)THEN
	IF (C_rot.GT.zero)THEN
		ANGLEFACEY=ACOS(C_rot/ROOT_rot)
	END IF
	IF (C_rot.LT.zero)THEN
		ANGLEFACEY=ACOS(C_rot/ROOT_rot)
	END IF
END IF
END SUBROUTINE ANGLEY


SUBROUTINE ANGLE2D(ANGLEFACEX,ANGLEFACEY)
IMPLICIT NONE
REAL,INTENT(INOUT)::ANGLEFACEX,ANGLEFACEY
real::length

length=distance2(N)

ANGLEFACEX=(vext(2,2)-vext(1,2))/length
ANGLEFACEy=-(vext(2,1)-vext(1,1))/length



END SUBROUTINE ANGLE2D










END MODULE TRANSFORM
