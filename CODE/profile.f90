MODULE PROFILE
USE DECLARATION
USE LIBRARY
USE BASIS

IMPLICIT NONE


 CONTAINS
 
 
 
 REAL FUNCTION LINEAR_INIT3D(n)
 !> @brief
!> This function initialises the solution for linear advection in 3D,
!> various customisable profiles can be generated and assigned to each initcond code
IMPLICIT NONE
INTEGER,INTENT(IN)::N
!COORDINATES=POX(1),POY(1),POZ(1)

IF (INITCOND.EQ.0)THEN
IF(((POX(1).GE.0.25D0).AND.(POX(1).LE.0.75D0)).AND.((POZ(1).GE.0.25D0).AND.(POZ(1).LE.0.75D0)))THEN
	LINEAR_INIT3d=10.0D0
ELSE
	LINEAR_INIT3d=1.0D0
END IF
END IF
IF (INITCOND.EQ.2)THEN
LINEAR_INIT3d=(SIN((2.0D0*PI)*(POX(1))))*&
(SIN((2.0D0*PI)*(POY(1))))*(SIN((2.0D0*PI)*(POZ(1))))

END IF

END FUNCTION LINEAR_INIT3D


 REAL FUNCTION LINEAR_INIT2D(n)
  !> @brief
!> This function initialises the solution for linear advection in 2D,
!> various customisable profiles can be generated and assigned to each initcond code
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(90)::polyfun
REAL::AADX,AADY,SUMF,rd
integer::ixg
!COORDINATES=POX(1),POY(1),POZ(1)
SUMF=zero
IF (INITCOND.EQ.1)THEN
! IF(((POX(1).GE.0.25D0).AND.(POX(1).LE.0.75D0)).AND.((POy(1).GE.0.25D0).AND.(POy(1).LE.0.75D0)))THEN
! 	LINEAR_INIT2d=1.0D0
! ELSE
	LINEAR_INIT2d=10.0D0

END IF
IF (INITCOND.EQ.0)THEN
! LINEAR_INIT2d=(SIN((2.0D0*PI)*(POX(1))))*&
!  (SIN((2.0D0*PI)*(POY(1))))
AADX=POX(1)
AADY=POy(1)
 compwrt=0

   polyfun(1:ielem(n,iconsidered)%idegfree)=basis_rec2d(N,AADX,AADY,ielem(n,iconsidered)%iorder,Iconsidered,ielem(n,iconsidered)%idegfree)
 compwrt=0
   do ixg=1,ielem(n,iconsidered)%idegfree
   SUMF=SUMF+polyfun(ixg)
   end do
linear_init2d=SUMF

END IF

IF (INITCOND.EQ.3)THEN

LINEAR_INIT2d=0.0d0
if (sqrt(((pox(1)-0.25d0)**2)+((poy(1)-0.5d0)**2)).le.0.15)then
rd=(1.0d0/0.15d0)*sqrt(((pox(1)-0.25d0)**2)+((poy(1)-0.5d0)**2))

LINEAR_INIT2d=0.25d0*(1.0d0+cos(pi*min(rd,1.0d0)))
end if

if (sqrt(((pox(1)-0.5d0)**2)+((poy(1)-0.25d0)**2)).le.0.15)then

rd=(1.0d0/0.15d0)*sqrt(((pox(1)-0.5d0)**2)+((poy(1)-0.25d0)**2))
LINEAR_INIT2d=1.0d0-rd
end if

    if (sqrt(((pox(1)-0.5d0)**2)+((poy(1)-0.75d0)**2)).le.0.15)then

    rd=(1.0d0/0.15d0)*sqrt(((pox(1)-0.5d0)**2)+((poy(1)-0.75d0)**2))
	  if ((abs(pox(1)-0.5).GE.0.025d0).or.(poy(1).gt.0.85))then

	  LINEAR_INIT2d=1.0d0
	  else

	  LINEAR_INIT2d=0.0d0

	  end if
    end if

end if

IF (INITCOND.EQ.2)THEN
    LINEAR_INIT2d= (SIN((2.0D0*PI)*(POX(1)))) * (SIN((2.0D0*PI)*(POY(1))))
end if

IF (INITCOND.EQ.100)THEN
    LINEAR_INIT2d = 1.0d0
end if

END FUNCTION LINEAR_INIT2D
 
 
 SUBROUTINE INITIALISE_EULER3D(N)
 IMPLICIT NONE
  !> @brief
!> This function initialises the solution for EULER and NAVIER-STOKES equations in 3D,
!> various customisable profiles can be generated and assigned to each initcond code
INTEGER,INTENT(IN)::N
!COORDINATES=POX,POY,POZ
!SOLUTION=VECCOS
!COMPONENTS FROM DAT FILE GAMMA,UVEL,WVEL,VVEL,PRES,RRES
!INITCOND= PROFILE CHOICE FROM DATA FILE
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX,VHX,AMP,DVEL
integer::u_cond1,u_cond2,u_cond3,u_cond4
VECCOS(:)=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ INITIAL VALUES FROM DAT FILE
R1=RRES
P1=PRES
S1=SQRT((GAMMA*P1)/(R1))
U1=UVEL
V1=VVEL
W1=WVEL

!KINETIC ENERGY FIRST!
SKIN1=(oo2)*((U1**2)+(V1**2)+(W1**2))
!INTERNAL ENERGY 

IE1=((P1)/((GAMMA-1.0D0)*R1))

!TOTAL ENERGY
E1=R1*(SKIN1+IE1)

!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=R1*W1
VECCOS(5)=E1
IF (TURBULENCE.EQ.1)THEN

  IF (TURBULENCEMODEL.EQ.1)THEN

  VECCOS(6)=VISC*TURBINIT
  END IF
  IF (TURBULENCEMODEL.EQ.2)THEN
 
   if (zero_turb_init .eq. 0) then
  VECCOS(6)=(1.5D0*I_turb_inlet*(ufreestream**2))*R1
   VECCOS(7)=R1*veccos(6)/(10.0e-5*visc)	
  end if
  
  if (zero_turb_init .eq. 1) then
  VECCOS(6)=(1.5D0*I_turb_inlet*(ufreestream**2))*R1
   VECCOS(7)=R1*veccos(6)/(10.0e-5*visc)	
  end if
    
  END IF

  
END IF
IF (PASSIVESCALAR.GT.0)THEN

  VECCOS(5+TURBULENCEEQUATIONS+1:5+TURBULENCEEQUATIONS+PASSIVESCALAR)=ZERO

END IF




IF (INITCOND.EQ.95)THEN	!TAYLOR GREEN INITIAL PROFILE
R1=1.0D0

W1=0.0D0
P1=100.0D0+((R1/16.0D0)*((COS(2.0D0*POZ(1)))+2.0d0)*((COS(2.0D0*POX(1)))+(COS(2.0D0*POY(1)))))
!p1=100.0d0+((r1/16.0d0)*(COS(2.0D0*POZ(1))+(2.0d0*COS(2.0D0*POX(1)))+(COS(2.0D0*POY(1)))-2.0D0 ))
!p1=100.0d0+((r1/16.0d0)*(COS(2.0D0*POZ(1))+(2.0d0*COS(2.0D0*POX(1)))+(COS(2.0D0*POY(1)))-2.0D0 ))
!UU=(SQRT((GAMMA*P1)/(R1)))*0.28
!P1=1.0
u1=sin(POX(1))*COS(POY(1))*COS(POZ(1))
v1=-COS(POX(1))*SIN(POY(1))*COS(POZ(1))
!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2)+(W1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=R1*W1
VECCOS(5)=E1
END IF

IF (INITCOND.EQ.101)THEN	!shock density interaction
if (pox(1).lt.-4.0d0)then
r1=3.8571d0
u1=2.6294d0
v1=zero
w1=zero
p1=10.333d0
else
r1=(1.0d0+0.2d0*sin(5.0d0*pox(1)))
u1=zero
v1=zero
w1=zero
p1=1
end if
SKIN1=(OO2)*((U1**2)+(V1**2)+(W1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=R1*W1
VECCOS(5)=E1





end if



IF (INITCOND.EQ.405)THEN
!TEST CASE 4.5 OF CORALIC & COLONIUS

IF (POX(1).GT.0.1D0)THEN
MP_R(1)=0.166315789d0
MP_R(2)=1.658d0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=-114.49D0
V1= 0.0D0
w1=0.0D0
P1=159060.0d0


! SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))



R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION



if (sqrt(((pox(1)-0.05d0)**2)+((poy(1)-0.05d0)**2)+((poz(1)-0.05d0)**2)).LE.0.025d0)then
MP_R(1)=0.166315789d0
MP_R(2)=1.204D0
MP_A(1)=0.95d0
MP_A(2)=0.05D0
U1=0.0D0
V1=0.0D0
w1=0.0d0
P1=101325

! SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW
else



MP_R(1)=0.166315789d0
MP_R(2)=1.204D0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=0.0D0
V1=0.0D0
w1=0.0d0
P1=101325


! SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if






END IF












VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=R1*w1
VECCOS(5)=E1
VECCOS(6)=MP_R(1)*MP_A(1)
VECCOS(7)=MP_R(2)*MP_A(2)
VECCOS(8)=MP_A(1)




END IF



IF (INITCOND.EQ.408)THEN
!TEST CASE 4.5 OF CORALIC & COLONIUS

IF (POX(1).GT.0.10D0)THEN
MP_R(1)=6.03
MP_R(2)=1.658d0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=-114.49D0
V1= 0.0D0
w1=0.0D0
P1=159060.0d0


! SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION



if (sqrt(((pox(1)-0.079d0)**2)+((poy(1)-0.035d0)**2)+((poz(1)-0.035d0)**2)).LE.(0.0325d0/2.0d0))then
MP_R(1)=6.03
MP_R(2)=1.204D0
MP_A(1)=1.0d0
MP_A(2)=0.0D0
U1=0.0D0
V1=0.0D0
w1=0.0d0
P1=101325

! SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW
else



MP_R(1)=6.03
MP_R(2)=1.204D0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=0.0D0
V1=0.0D0
w1=0.0d0
P1=101325


! SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if






END IF












VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=R1*w1
VECCOS(5)=E1
VECCOS(6)=MP_R(1)*MP_A(1)
VECCOS(7)=MP_R(2)*MP_A(2)
VECCOS(8)=MP_A(1)




END IF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE INITIALISE_EULER3D



 SUBROUTINE INITIALISE_EULER2D(N)
 IMPLICIT NONE
   !> @brief
!> This function initialises the solution for EULER and NAVIER-STOKES equations in 2D,
!> various customisable profiles can be generated and assigned to each initcond code
INTEGER,INTENT(IN)::N
!COORDINATES=POX,POY
!SOLUTION=VECCOS
!COMPONENTS FROM DAT FILE GAMMA,UVEL,WVEL,VVEL,PRES,RRES
!INITCOND= PROFILE CHOICE FROM DATA FILE
real::acp,mscp,mvcp,vmcp,bcp,rcp,tcp,vfr,theta1
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX,VHX,AMP,DVEL,rgg,tt1
real::pr_Radius,pr_beta,pr_machnumberfree,pr_pressurefree,pr_temperaturefree,pr_gammafree,pr_Rgasfree,pr_xcenter,pr_ylength,pr_xlength,pr_ycenter,pr_densityfree,pr_cpconstant,pr_radiusvar,pr_velocityfree,pr_TemperatureVar
integer::u_cond1,u_cond2,u_cond3,u_cond4
VECCOS(:)=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! READ INITIAL VALUES FROM DAT FILE
R1=RRES
P1=PRES
S1=SQRT((GAMMA*P1)/(R1))
U1=UVEL
V1=VVEL

!KINETIC ENERGY FIRST!
SKIN1=(oo2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 

IE1=((P1)/((GAMMA-1.0D0)*R1))

!TOTAL ENERGY
E1=R1*(SKIN1+IE1)

!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
IF (TURBULENCE.EQ.1)THEN

  IF (TURBULENCEMODEL.EQ.1)THEN

  VECCOS(5)=VISC*TURBINIT/R1
  END IF
  IF (TURBULENCEMODEL.EQ.2)THEN
 
   if (zero_turb_init .eq. 0) then
  VECCOS(5)=(1.5D0*(I_turb_inlet*ufreestream)**2)*R1
  VECCOS(6)=ufreestream/L_turb_inlet
 VECCOS(6)=(C_MU_INLET**(-0.25D0))*SQRT(veccos(5))&
			/L_TURB_INLET*r1	
  end if
  
  if (zero_turb_init .eq. 1) then
  VECCOS(5)=ZERO
  VECCOS(6)=ufreestream/L_turb_inlet
  end if
    
  END IF

  
END IF
IF (PASSIVESCALAR.GT.0)THEN

  VECCOS(4+TURBULENCEEQUATIONS+1:nof_Variables+TURBULENCEEQUATIONS+PASSIVESCALAR)=ZERO

END IF




IF (INITCOND.EQ.95)THEN	!TAYLOR GREEN INITIAL PROFILE
IF ((POY(1).GE.0.25D0).AND.(POY(1).LE.0.75D0))then

r1=2.0d0
u1=-0.5d0
v1=0.01d0*sin(2.0d0*PI*(POX(1)-0.5))
P1=2.5


end if

IF ((POY(1).lt.0.25D0).or.(POY(1).gt.0.75D0))then
r1=1.0d0
u1=0.5d0
v1=0.01d0*sin(2.0d0*PI*(POX(1)-0.5))
P1=2.5
END IF



!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF


IF (INITCOND.EQ.75)THEN	!TAYLOR GREEN INITIAL PROFILE

if (((pox(1).ge.(2.0d0)).and.(pox(1).le.(4.0d0))).and.((poy(1).ge.(2.0d0)).and.(poy(1).le.(4.0d0))))then
R1=1.0D0
v1=0.0
u1=0.0
p1=10.0d0
else
R1=0.2D0
v1=0.0
u1=0.0
p1=1.0d0


end if


!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF




IF (INITCOND.EQ.25)THEN	!TAYLOR GREEN INITIAL PROFILE

if (sqrt(((pox(1)-1.0d0)**2)+((poy(1)-1.0d0)**2)).LE.0.4d0)then

R1=1.0D0
v1=0.0
u1=0.0
p1=1.0d0


else
R1=0.125D0
v1=0.0
u1=0.0
p1=0.1d0



end if


!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF

IF (INITCOND.EQ.31)THEN	!TAYLOR GREEN INITIAL PROFILE

if (sqrt(((pox(1)-1.0d0)**2)+((poy(1)-1.0d0)**2)).LE.0.4d0)then

R1=1.0D0
v1=0.0
u1=0.0
p1=1.0d0


else
R1=0.125D0
v1=0.0
u1=0.0
p1=0.1d0



end if


!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF



IF (INITCOND.EQ.65)THEN	!TAYLOR GREEN INITIAL PROFILE



R1=1.0D0
u1=1.0d0
v1=1.0d0
p1=1.0d0
!rgg=((pox(1)**2)+(poy(1)**2))
rgg=(((pox(1)-5.0d0)**2)+((poy(1)-5.0d0)**2))
u1=u1+(((5.0d0/(2.0d0*pi)))*(exp(((1.0d0-rgg)/(2.0d0))))*(5.0d0-poy(1)))
v1=v1+(((5.0d0/(2.0d0*pi)))*(exp(((1.0d0-rgg)/(2.0d0))))*(pox(1)-5.0d0))
tt1=1.0d0-(((gamma-1)*(25.0d0/(8.0d0*pi**2*gamma)))*exp(1.0d0-rgg))
r1=tt1**(1.0d0/(gamma-1.0d0))
P1=R1*tt1


    
!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF



IF (INITCOND.EQ.100)THEN	!TAYLOR GREEN INITIAL PROFILE
acp=0.075d0
bcp=0.175d0
mscp=1.5d0
mvcp=0.7d0
vmcp=sqrt(gamma)*mvcp

rcp=sqrt((pox(1)-0.25d0)**2+(poy(1)-0.5d0)**2)

if (rcp.le.acp)then
vfr=vmcp*rcp/acp
tcp=((rcp-0.175)*(gamma-1.0d0)*(vfr**2)/(rcp*gamma))+1.0d0
p1=tcp**(gamma/(gamma-1.0d0))
r1=tcp**(1.0d0/(gamma-1.0d0))
theta1=atan((poy(1)-0.5d0)/(pox(1)-0.25d0))
v1=vfr*COS(theta1)
u1=1.5d0*sqrt(gamma)-vfr*SIN(theta1)



ELSE

 if ((rcp.ge.acp).and.(rcp.le.bcp))then
vfr=vmcp*(acp/(acp**2-bcp**2))*(rcp-((bcp**2)/rcp))

vfr=vmcp*rcp/acp
tcp=((rcp-0.175)*(gamma-1.0d0)*(vfr**2)/(rcp*gamma))+1.0d0
p1=tcp**(gamma/(gamma-1.0d0))
r1=tcp**(1.0d0/(gamma-1.0d0))
theta1=atan((poy(1)-0.5d0)/(pox(1)-0.25d0))
v1=vfr*COS(theta1)
u1=1.5d0*sqrt(gamma)-vfr*SIN(theta1)








else
vfr=0.0d0
r1=1.0d0
u1=1.5d0*sqrt(gamma)
v1=0.0
p1=1.0d0
end if
END IF

if (pox(1).gt.0.5)then
r1=(9.0d0*gamma+9.0d0)/(9.0d0*gamma-1.0d0)
u1=sqrt(gamma)*((9.0d0*gamma-1.0d0)/(6.0d0*(gamma+1.0d0)))
p1=(7.0d0*gamma+2.0d0)/(2.0d0*gamma+2.0d0)
v1=0.0d0
end if


    
!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
END IF


IF (INITCOND.EQ.101)THEN	!shock density interaction
if (pox(1).lt.-4.0d0)then
r1=3.8571d0
u1=2.6294d0
v1=zero

p1=10.333d0
else
r1=(1.0d0+0.2d0*sin(5.0d0*pox(1)))
u1=zero
v1=zero

p1=1
end if
    
!KINETIC ENERGY FIRST!
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
end if


IF (INITCOND.EQ.102)THEN	!shock density interaction
if (pox(1).lt.((1.0d0/6.0d0)+(poy(1)/(sqrt(3.0d0)))))then
r1=8.0d0
u1=8.25*cos(pi/6.0d0)
v1=-8.25*sin(pi/6.0d0)
p1=116.5
else
r1=1.4d0
u1=zero
v1=zero
p1=1.0d0
end if
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1





end if

IF (INITCOND.EQ.104)THEN	!dmr_domain1
if (pox(1).lt.(1.0d0/6.0d0))then
r1=8.0d0
u1=8.25d0
v1=0.0d0
p1=116.5
else
r1=1.4d0
u1=zero
v1=zero
p1=1.0d0
end if
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1





end if





IF (INITCOND.EQ.30)THEN	!shock density interaction
if (pox(1).le.1.0d0)then
if (poy(1).le.1.0d0)then
r1=77.0d0/558.0d0
u1=4.0d0/sqrt(11.0d0)
v1=4.0d0/sqrt(11.0d0)
p1=9.0d0/310.0d0
end if
if (poy(1).gt.1.0d0)then
r1=33.0d0/62.0d0
u1=4.0d0/sqrt(11.0d0)
v1=0.0d0
p1=0.3d0
end if
end if
if (pox(1).gt.1.0d0)then
if (poy(1).le.1.0d0)then
r1=33.0d0/62.0d0
u1=0.0d0
v1=4.0d0/sqrt(11.0d0)
p1=0.3d0
end if
if (poy(1).gt.1.0d0)then
r1=1.5d0
u1=0.0d0
v1=0.0d0
p1=1.5d0
end if
end if
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1





end if

IF (INITCOND.EQ.3)THEN	!shock density interaction
if (pox(1).ge.zero)then
r1=1.1175d0
u1=0.0d0
v1=0.0d0
p1=95000.0d0
else
r1=1.7522
u1=166.34345
v1=0.0d0
p1=180219.75d0
end if
SKIN1=(OO2)*((U1**2)+(V1**2))
!INTERNAL ENERGY 
IE1=((P1)/((GAMMA-1.0D0)*R1))
!TOTAL ENERGY
E1=(P1/(GAMMA-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1





end if





IF (INITCOND.EQ.123)THEN	!fast inviscid vortex evolution



pr_Radius = 0.005
pr_beta = 0.2
pr_machnumberfree = 0.05
pr_pressurefree = 100000.0
pr_temperaturefree = 300.0
pr_gammafree = 1.4
pr_Rgasfree = 287.15
pr_xcenter = 0.05 ; pr_ycenter = 0.05
pr_xlength = 0.1 ;  pr_ylength = 0.1

pr_densityfree = pr_pressurefree/ (pr_Rgasfree*pr_temperaturefree)
 pr_cpconstant = pr_Rgasfree*(pr_gammafree/(pr_gammafree-1.0))

pr_radiusvar = sqrt( (POX(1)-pr_xcenter)**2 + (POY(1)-pr_ycenter)**2 )/ pr_Radius

pr_velocityfree = pr_machnumberfree * sqrt(pr_gammafree*pr_Rgasfree*pr_temperaturefree)


u1 = pr_velocityfree * (1.0 - (pr_beta* ((POY(1)-pr_ycenter)/pr_Radius )*exp((-pr_radiusvar**2)/2)  ) )
v1 = pr_velocityfree * ((pr_beta* ((POY(1)-pr_ycenter)/pr_Radius )*exp((-pr_radiusvar**2)/2)  ) )
pr_TemperatureVar = pr_temperaturefree - ( ((pr_velocityfree**2 *pr_beta**2.0)/(2.0 * pr_cpconstant))*exp(-pr_radiusvar**2) )

r1 = pr_densityfree * (pr_TemperatureVar/pr_temperaturefree)**(1.0/(pr_gammafree-1.0))

 p1 = r1*pr_Rgasfree*pr_TemperatureVar
! 
! W1=0.0D0
! P1=100.0D0+((R1/16.0D0)*((2.0D0*(COS(2.0D0*POX(1))))+(COS(2.0D0*POY(1)))-2.0D0))
! !UU=(SQRT((GAMMA*P1)/(R1)))*0.28
! !P1=1.0
! u1=sin(POX(1))*COS(POY(1))
! v1=-COS(POX(1))*SIN(POY(1))
! !KINETIC ENERGY FIRST!
 SKIN1=(OO2)*((U1**2)+(V1**2))
! !INTERNAL ENERGY 
 IE1=((P1)/((pr_gammafree-1.0D0)*R1))
! !TOTAL ENERGY
 E1=(P1/(pr_gammafree-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1




END IF



IF (INITCOND.EQ.124)THEN	!fast inviscid vortex evolution
pr_Radius = 0.005
pr_beta = 0.2
pr_machnumberfree = 0.2
pr_pressurefree = 100000.0
pr_temperaturefree = 300.0
pr_gammafree = 1.4
pr_Rgasfree = 287.15
pr_xcenter = 0.05 ; pr_ycenter = 0.05
pr_xlength = 0.1 ;  pr_ylength = 0.1

pr_densityfree = pr_pressurefree/ (pr_Rgasfree*pr_temperaturefree)
 pr_cpconstant = pr_Rgasfree*(pr_gammafree/(pr_gammafree-1.0))

pr_radiusvar = sqrt( (POX(1)-pr_xcenter)**2 + (POY(1)-pr_ycenter)**2 )/ pr_Radius

pr_velocityfree = pr_machnumberfree * sqrt(pr_gammafree*pr_Rgasfree*pr_temperaturefree)


u1 = pr_velocityfree * (1.0 - (pr_beta* ((POY(1)-pr_ycenter)/pr_Radius )*exp((-pr_radiusvar**2)/2)  ) )
v1 = pr_velocityfree * ((pr_beta* ((POY(1)-pr_ycenter)/pr_Radius )*exp((-pr_radiusvar**2)/2)  ) )
pr_TemperatureVar = pr_temperaturefree - ( ((pr_velocityfree**2 *pr_beta**2.0)/(2.0 * pr_cpconstant))*exp(-pr_radiusvar**2) )

r1 = pr_densityfree * (pr_TemperatureVar/pr_temperaturefree)**(1.0/(pr_gammafree-1.0))

 p1 = r1*pr_Rgasfree*pr_TemperatureVar
! 
! W1=0.0D0
! P1=100.0D0+((R1/16.0D0)*((2.0D0*(COS(2.0D0*POX(1))))+(COS(2.0D0*POY(1)))-2.0D0))
! !UU=(SQRT((GAMMA*P1)/(R1)))*0.28
! !P1=1.0
! u1=sin(POX(1))*COS(POY(1))
! v1=-COS(POX(1))*SIN(POY(1))
! !KINETIC ENERGY FIRST!
 SKIN1=(OO2)*((U1**2)+(V1**2))
! !INTERNAL ENERGY 
 IE1=((P1)/((pr_gammafree-1.0D0)*R1))
! !TOTAL ENERGY
 E1=(P1/(pr_gammafree-1))+(R1*SKIN1)
!VECTOR OF CONSERVED VARIABLES NOW
VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1




END IF






IF (INITCOND.EQ.401)THEN

if (sqrt(((pox(1)-1.0d0)**2)+((poy(1)-1.0d0)**2)).LE.0.4d0)then
MP_R(1)=1.0D0
MP_R(2)=0.125
MP_A(1)=1.0D0
MP_A(2)=0.0D0
U1=ZERO
V1=ZERO
P1=1.0D0

else
MP_R(1)=1.0D0
MP_R(2)=0.125
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=ZERO
V1=ZERO
P1=0.1D0

end if




R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW




VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)




END IF


IF (INITCOND.EQ.402)THEN

if ((pox(1).ge.0.25).and.(pox(1).lt.0.75))then
MP_R(1)=10.0D0
MP_R(2)=1.0D0
MP_A(1)=1.0D0
MP_A(2)=0.0D0
U1=0.5D0
V1=0.0D0
P1=1.0D0/1.4D0

else
MP_R(1)=10.0D0
MP_R(2)=1.0D0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=0.5D0
V1=0.0D0
P1=1.0D0/1.4D0

end if




R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW




VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)




END IF



IF (INITCOND.EQ.403)THEN
!TEST CASE 4.1 OF WANG, DEITERDING, PAN and REN


MP_R(1)=7.0D0
MP_R(2)=1.0D0
MP_A(1)=0.5D0+0.25d0*(SIN(PI*((pox(1)-1.0d0)+0.5d0)))
MP_A(2)=1.0D0-MP_A(1)
U1=1.0D0
V1=0.0D0
P1=1.0D0/1.4D0






R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW




VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)




END IF




IF (INITCOND.EQ.405)THEN
!TEST CASE 4.5 OF CORALIC & COLONIUS

IF (POX(1).GT.0.1D0)THEN
MP_R(1)=0.166315789d0
MP_R(2)=1.658d0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=-114.49D0
V1= 0.0D0
P1=159060.0d0


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION



if (sqrt(((pox(1)-0.05d0)**2)+((poy(1)-0.05d0)**2)).LE.0.025d0)then
MP_R(1)=0.166315789d0
MP_R(2)=1.204D0
MP_A(1)=0.95d0
MP_A(2)=0.05D0
U1=0.0D0
V1=0.0D0
P1=101325

R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW
else



MP_R(1)=0.166315789d0
MP_R(2)=1.204D0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=0.0D0
V1=0.0D0
P1=101325


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if






END IF












VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)




END IF







IF (INITCOND.EQ.406)THEN
!TEST CASE 4.5 OF CORALIC & COLONIUS

IF (POX(1).GT.0.1D0)THEN
MP_R(1)=0.166315789d0
MP_R(2)=1.658d0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=-114.49D0
V1= 0.0D0
P1=159060.0d0


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION
u_cond1=0;u_cond2=0;u_cond3=0; u_cond4=0

 if (((pox(1).ge.0.01d0).and.(pox(1).le.0.02)).or.((pox(1).ge.0.03d0).and.(pox(1).le.0.04)))then
    u_cond1=1
 end if
 if ((poy(1).ge.0.05d0).and.(poy(1).le.0.075))then
    u_cond2=1
 end if
 
 
 if ((poy(1).le.0.05d0).and.(poy(1).gt.0.02))then
    u_cond3=1
 end if
 if (u_cond3.eq.1)then
 if ((sqrt(((pox(1)-0.025d0)**2)+((poy(1)-0.05d0)**2)).LE.0.015d0).and.(sqrt(((pox(1)-0.025d0)**2)+((poy(1)-0.05d0)**2)).gE.0.005d0)) then
    u_cond4=1
 end if
 end if
 
 
 if (((u_cond1.eq.1).and.(u_cond2.eq.1)).or.(u_cond4.eq.1))then
 
MP_R(1)=0.166315789d0
MP_R(2)=1.204D0
MP_A(1)=0.95d0
MP_A(2)=0.05D0
U1=0.0D0
V1=0.0D0
P1=101325

R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW
else



MP_R(1)=0.166315789d0
MP_R(2)=1.204D0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=0.0D0
V1=0.0D0
P1=101325


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if






END IF












VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)




END IF





IF (INITCOND.EQ.408)THEN
!TEST CASE 4.5 OF CORALIC & COLONIUS

IF (POX(1).GT.0.10D0)THEN
MP_R(1)=6.03
MP_R(2)=1.658d0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=-114.49D0
V1= 0.0D0
w1=0.0D0
P1=159060.0d0


! SKIN1=(OO2)*((U1**2)+(V1**2)+(w1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION



if (sqrt(((pox(1)-0.079d0)**2)+((poy(1)-0.035d0)**2)).LE.(0.0325d0/2.0d0))then
MP_R(1)=6.03
MP_R(2)=1.204D0
MP_A(1)=1.0d0
MP_A(2)=0.0D0
U1=0.0D0
V1=0.0D0
w1=0.0d0
P1=101325.0d0

! SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW
else



MP_R(1)=6.03
MP_R(2)=1.204D0
MP_A(1)=0.0D0
MP_A(2)=1.0D0
U1=0.0D0
V1=0.0D0
P1=101325.0d0


! SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if



END IF



VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)


END IF










IF (INITCOND.EQ.410)THEN
!TEST CASE 6.2 of The Euler Equations for Multiphase Compressible Flow in Conservation Form Simulation of Shock–Bubble Interactions, Robin K. S. Hankin

!GAMMA_IN(1) = 7.15 ! Water
!GAMMA_IN(2) = 1.4  ! Air
!MP_PINF(1) = 3.43e8 !Water from Coralic and Colonius or 2.218e8(abgrall203)
!MP_PINF(2) = 0 ! Air 

IF (POX(1).LE.-0.007D0)THEN
MP_R(1)=1225.6d0 ! Water density
MP_R(2)=1.2d0 ! Air density
MP_A(1)=1.0D0 ! Water volume fraction (everything is water here)
MP_A(2)=0.0D0 ! Air volume fraction
U1=542.76D0	  ! m/s
V1= 0.0D0
P1=1.6e9      ! Pa


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION



if (sqrt(((pox(1))**2)+((poy(1))**2)).LE.0.003d0)then
MP_R(1)=1225.6d0 ! Water density
MP_R(2)=1.2d0 ! Air density
MP_A(1)=0.0D0 ! Water volume fraction (everything is water here)
MP_A(2)=1.0D0 ! Air volume fraction
U1= 0.0D0	  ! m/s
V1= 0.0D0
P1= 101325    ! Pa

R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW
else



MP_R(1)=1000.0d0 ! Water density
MP_R(2)=1.2d0 ! Air density
MP_A(1)=1.0D0 ! Water volume fraction (everything is water here)
MP_A(2)=0.0D0 ! Air volume fraction
U1= 0.0D0     ! m/s
V1= 0.0D0
P1=101325     ! Pa


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if

END IF

VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)

END IF




IF (INITCOND.EQ.411)THEN
!EXAMPLE VI Paper5.pdf

!GAMMA_IN(1) = 4.4 ! Water
!GAMMA_IN(2) = 1.4  ! Air
!MP_PINF(1) = 6e8 !Water from Coralic and Colonius or 2.218e8(abgrall203)
!MP_PINF(2) = 0 ! Air 

IF (POX(1).LE.0.0066D0)THEN
MP_R(2)=1323.65d0 	! Water density
MP_R(1)=1d0 		! Air density
MP_A(2)=1.0D0 		! Water volume fraction (everything is water here)
MP_A(1)=0.0D0 		! Air volume fraction
U1=681.058D0	  	! m/s
V1= 0.0D0
P1=1.9e9      		! Pa


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

!FIRST WITHIN BUBBLE REGION



if (sqrt(((pox(1)-0.012)**2)+((poy(1)-0.012)**2)).LE.0.003d0)then
MP_R(2)=1000.00d0 	! Water density
MP_R(1)=1d0 		! Air density
MP_A(2)=0.0D0 		! Water volume fraction (everything is water here)
MP_A(1)=1.0D0 		! Air volume fraction
U1= 0.0D0	  		! m/s
V1= 0.0D0
P1= 100000    			! Pa

R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW


else

MP_R(2)=1000.0d0 	! Water density
MP_R(1)=1d0 		! Air density
MP_A(2)=1.0D0 		! Water volume fraction (everything is water here)
MP_A(1)=0.0D0 		! Air volume fraction
U1= 0.0D0	  		! m/s
V1= 0.0D0
P1= 100000    			! Pa


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1

!VECTOR OF CONSERVED VARIABLES NOW
end if

END IF

VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)

END IF


IF (INITCOND.EQ.412)THEN
!EXAMPLE 4.3 Paper2.pdf

!GAMMA_IN(1) = 1.4 ! air
!GAMMA_IN(2) = 5.5  ! water
!MP_PINF(1) = 0 !Water from Coralic and Colonius or 2.218e8(abgrall203)
!MP_PINF(2) = 1.505 ! Air 

IF (POX(1).LE.0.00)THEN
MP_R(1)=1.241 	! air density
MP_R(2)=0.991 		! water density
MP_A(1)=1.0D0 		
MP_A(2)=0.0D0 		
U1=0.0	  	! m/s
V1= 0.0D0
P1=2.753     		! Pa


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

ELSE

MP_R(1)=1.241 	! air density
MP_R(2)=0.991 		! water density
MP_A(1)=0.0D0 		
MP_A(2)=1.0D0 		
U1=0.0	  	! m/s
V1= 0.0D0
P1=3.059*10e-4     		! Pa


R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW

END IF

VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)

END IF

IF (INITCOND.EQ.420)THEN

MP_R(1)=7.0 	! air density
MP_R(2)=1.0 		! water density
MP_A(1)=0.5d0+0.25*sin(pi*(((pox(1)-0.5d0)*2.0d0))) 		
MP_A(2)=1.0D0-MP_A(1)
U1=1.0D0	  	! m/s
V1= 0.0D0
P1=1.0d0/1.4d0    		! Pa




R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW



VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)

END IF


IF (INITCOND.EQ.421)THEN


IF (POX(1).LT.0.0D0)THEN

MP_R(1)=1.241 	! air density
MP_R(2)=0.991 		! water density
MP_A(1)=1.0D0		
MP_A(2)=0.0D0
U1=0.0D0	  	! m/s
V1= 0.0D0
P1=2.753D0    		! Pa
ELSE
MP_R(1)=1.241 	! air density
MP_R(2)=0.991 		! water density
MP_A(1)=0.0D0		
MP_A(2)=1.0D0
U1=0.0D0	  	! m/s
V1= 0.0D0
P1=3.059*(10.0e-4)    		! Pa


END IF


SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
! MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)*MP_R(1)))
! MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)*MP_R(2)))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))

IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW



VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)

END IF

IF (INITCOND.EQ.430)THEN


IF (POy(1).LT.0.386D0)THEN

MP_R(1)=3.483 	! air density
MP_R(2)=867 		! water density
MP_A(1)=0.0D0		
MP_A(2)=1.0D0
U1=0.0D0	  	! m/s
V1= 2
P1=300000    		! Pa
ELSE
MP_R(1)=3.483 	! air density
MP_R(2)=867		! water density
MP_A(1)=1.0D0		
MP_A(2)=0.0D0
U1=0.0D0	  	! m/s
V1= 0.0D0
P1=300000    		! Pa


END IF


SKIN1=(OO2)*((U1**2)+(V1**2))
R1=(MP_R(1)*MP_A(1))+(MP_R(2)*MP_A(2))
! MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)*MP_R(1)))
! MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)*MP_R(2)))
MP_IE(1)=((P1+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P1+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))

IE1=(MP_IE(1)*MP_A(1))+(MP_IE(2)*MP_A(2))
SKIN1=(OO2)*((U1**2)+(V1**2))
E1=(R1*SKIN1)+IE1
!VECTOR OF CONSERVED VARIABLES NOW



VECCOS(1)=R1
VECCOS(2)=R1*U1
VECCOS(3)=R1*V1
VECCOS(4)=E1
VECCOS(5)=MP_R(1)*MP_A(1)
VECCOS(6)=MP_R(2)*MP_A(2)
VECCOS(7)=MP_A(1)

END IF









!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE INITIALISE_EULER2D





END MODULE PROFILE
