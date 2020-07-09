MODULE PROFILE
USE DECLARATION
USE LIBRARY
USE BASIS

IMPLICIT NONE


 CONTAINS
 
 
 
 REAL FUNCTION LINEAR_INIT3D(n)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
!COORDINATES=POX(1),POY(1),POZ(1)

IF (INITCOND.EQ.0)THEN
IF(((POX(1).GE.0.25D0).AND.(POX(1).LE.0.75D0)).AND.((POZ(1).GE.0.25D0).AND.(POZ(1).LE.0.75D0)))THEN
	LINEAR_INIT3d=1.0D0
ELSE
	LINEAR_INIT3d=0.0D0
END IF
END IF
IF (INITCOND.EQ.2)THEN
LINEAR_INIT3d=(SIN((2.0D0*PI)*(POX(1))))*&
(SIN((2.0D0*PI)*(POY(1))))*(SIN((2.0D0*PI)*(POZ(1))))

END IF

END FUNCTION LINEAR_INIT3D


 REAL FUNCTION LINEAR_INIT2D(n)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(90)::polyfun
REAL::AADX,AADY,SUMF,rd
integer::ixg
!COORDINATES=POX(1),POY(1),POZ(1)
SUMF=zero
IF (INITCOND.EQ.1)THEN
IF(((POX(1).GE.0.25D0).AND.(POX(1).LE.0.75D0)).AND.((POy(1).GE.0.25D0).AND.(POy(1).LE.0.75D0)))THEN
	LINEAR_INIT2d=1.0D0
ELSE
	LINEAR_INIT2d=0.0D0
END IF
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
LINEAR_INIT2d=(SIN((2.0D0*PI)*(POX(1))))*&
 (SIN((2.0D0*PI)*(POY(1))))
end if


END FUNCTION LINEAR_INIT2D
 
 
 SUBROUTINE INITIALISE_EULER3D(N)
INTEGER,INTENT(IN)::N
!COORDINATES=POX,POY,POZ
!SOLUTION=VECCOS
!COMPONENTS FROM DAT FILE GAMMA,UVEL,WVEL,VVEL,PRES,RRES
!INITCOND= PROFILE CHOICE FROM DATA FILE
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX,VHX,AMP,DVEL
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
    IF (SRF.EQ.0) THEN
        VECCOS(6)=(1.5D0*I_turb_inlet*(ufreestream**2))*R1
        VECCOS(7)=R1*veccos(6)/(10.0e-5*visc)	
   ELSE
        VECCOS(6)=(1.5D0*I_turb_inlet*(KINIT_SRF**2))*R1
        VECCOS(7)=R1*veccos(6)/(10.0e-5*visc)
   END IF	
  end if
  
  if (zero_turb_init .eq. 1) then
    IF (SRF.EQ.0) THEN
        VECCOS(6)=(1.5D0*I_turb_inlet*(ufreestream**2))*R1
        VECCOS(7)=R1*veccos(6)/(10.0e-5*visc)	
   ELSE
        VECCOS(6)=(1.5D0*I_turb_inlet*(KINIT_SRF**2))*R1
        VECCOS(7)=R1*veccos(6)/(10.0e-5*visc)
   END IF		
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



!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE INITIALISE_EULER3D



 SUBROUTINE INITIALISE_EULER2D(N)
INTEGER,INTENT(IN)::N
!COORDINATES=POX,POY
!SOLUTION=VECCOS
!COMPONENTS FROM DAT FILE GAMMA,UVEL,WVEL,VVEL,PRES,RRES
!INITCOND= PROFILE CHOICE FROM DATA FILE
real::acp,mscp,mvcp,vmcp,bcp,rcp,tcp,vfr,theta1
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX,VHX,AMP,DVEL,rgg,tt1
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

  VECCOS(4+TURBULENCEEQUATIONS+1:4+TURBULENCEEQUATIONS+PASSIVESCALAR)=ZERO

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
if (pox(1).le.zero)then
if (poy(1).le.zero)then
r1=0.138
u1=1.206
v1=1.206
p1=0.029
end if
if (poy(1).gt.zero)then
r1=0.5323
u1=1.206
v1=0.0
p1=0.3
end if
end if
if (pox(1).gt.zero)then
if (poy(1).le.zero)then
r1=0.5323
u1=0.0
v1=1.206
p1=0.3
end if
if (poy(1).gt.zero)then
r1=1.5
u1=0.0
v1=0.0
p1=1.5
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE INITIALISE_EULER2D





END MODULE PROFILE
